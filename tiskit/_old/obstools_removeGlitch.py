#!env python3
""" Model and remove glitches from BBOBS data"""

#################################################
import obspy
from obspy.core import UTCDateTime
from obspy.core.stream import Stream
import diracComb_v2 as dc
import matplotlib.pyplot as plt
import numpy as np
#from scipy import signal
import math as M
import sys
import glob

#============================================================================
## CONVENIENCE CLASS FOR HOLDING GLITCH PARMETERS
#============================================================================
class GlitchParm:
  """ 
  Holds information about one periodic glitch signal
  
  starttime corresponds to onset of glitch, data will be sliced so that
    this arrives 1/3 of the way into each slice
  
  Attributes:
    period: best guess period for glitch (seconds)
    dp:      how much to change the period by when testing for better values
    clip_lo: lowest value that glitch attains (will clip lower values when
                calculating glitch)
    clip_hi: highest value that glitch attains (ditto)
    starttime: onset time of earliest glitch.  Can be a list, which allows
               a jump in onset time (for example, after a seismometer recenter
               for the 1h glitches).
               It is better to be a few seconds too early than too late, as
               the program will make glitch slices starting between here and
               1/3 of a glitch period earlier
  """
  def __init__(self,period,dp,clip_lo,clip_hi,starttime=None) :
    self.period=float(period)
    self.clip_lo=float(clip_lo)
    self.clip_hi=float(clip_hi)
    self.dp=float(dp)
    if isinstance(starttime,UTCDateTime):
      self.starttime=starttime
    elif starttime is None:
      self.starttime=None
    elif isinstance(starttime,list):
        self.starttime=starttime
            for s in self.starttime:
                s=UTCDateTime(s)    
    else:
      self.starttime=UTCDateTime(starttime)

#============================================================================
# The main routine
#============================================================================
def removeGlitch(stream,glitch_parms,eqfile='earthquakes_201210-201310.csv',
                    interactive=False, plots=False,
                    glitch=None,diracComb=None) :
    """
    calculate and remove periodic glitches from one channel of seismological
    data
    
    Usage: removeGlitch(trace,glitch_parms,eqfile,interactive,glitch=[])
    
    Inputs:
        trace = seismological channel
        glitch_parms = list of GlitchParm objects
        eqfile = csv file containing global earthquakes
                (column 1 = time (ISO8601), column4=magnitude)
        interactive (boolean): whether to let the user view and modify glitch 
                                interval and clip levels
        plots (boolean):  Whether to plot stuff
        #### THESE NEXT TWO ARE NOT YET IMPLEMENTED, MEANT FOR REMOVE
        #### GLITCHES BASED ON PRIOR INFORMATION
        glitches (list of numpy.ndarrays): glitches to use (will calculate the
                    best multiplication factor to fit the given trace) [None]
        diracComb (list of obspy traces): traces to convolve with glitches 
                to generate synthetic glitch time series [None]
                                
    Outputs:
        new_stream : (obspy.trace) glitch-corrected trace
        glitches   : (list of numpy.ndarray) list of average Glitches (one per
                     element in glitch_parms)
        diracCombs : (list of obspy traces): forms to convolve with
                      glitch to get synthetic glitch time series
        nGlitcheses: (list of ints): number of data glitches used to calculate
                                    average glitch(es)
    """
        
    #############################################
    ###                CONSTANTS              ###
    #############################################     
    match=True  # Match and cancel each pulse separately. In principal this is a
            # variable, but it works so well you should always do it!

    #############################################
    starttime=stream[0].stats.starttime
    endtime=stream[0].stats.endtime

    trace=stream.select(component='Z')[0]
    sps=trace.stats.sampling_rate
    station=trace.stats.station
    print('Station='+station)

    #######################################################################
    # WEILANDT METHOD TO REMOVE GLITCHES
    #######################################################################

    #======================
    # Bandpass to 2-30 mHz: demean; LP 30 mHz, order 6 (dt=1.6s); HP  2 mHz, order 2 (dt=1.6s)
    trace.detrend('demean')
    trace.detrend('linear')
    trace=trace.filter('lowpass',freq=0.03,corners=6)
    trace=trace.filter('highpass',freq=0.002,corners=2)
    
    #======================
    # Make a template to cut out EQ times, clip extreme values...)
    template=makeEQTemplate(trace,eqfile,plotit=plots)
  
    #===========================================================
    # Loop through different glitch periods (1 hour, 3 hours, ...) and remove the glitches
    #===========================================================
    iGlitchType=0
    glitches=list()
    dirac_combs=list()
    nGlitcheses=list()
    for gp in glitch_parms :
        iGlitchType=iGlitchType+1
        print('='*75)
        print(' '*15 + 'Fitting {:g}-second glitch'.format(gp.period))
        print('='*75)
        
        shifts=[]
        if isinstance(gp.starttime,list):
            # Search for glitch start corresponding to data start and 
            # make a list of glitch start shifts within data
           ref_time=gp.starttime
            if len(gp.starttime)>1:
                for test_time in gp.starttime[1:]:
                    if test_time < starttime:
                        ref_time=test_time
                    elif test_time < endtime:
                        # there is a break in glitch times within this data
                        shifts.append(test_time)
        elif gp.starttime:
            ref_time=gp.starttime
        else:
            ref_time=trace.stats.starttime
        
        if shifts:
            starttimes=shifts.copy().insert(0,starttime)
            endtimes=shifts.copy().append(endtime)
            reftimes=shifts.copy().insert(0,ref_time)
            print('Should have some way to cut out last period (because will be full of leveling)')
            for starttime,endtime,reftime in zip(starttimes,endtimes,reftimes):
                tr_slice=trace.cut(starttime=starttime,endtime=endtime)
                te_slice=template.cut(starttime=starttime,endtime=endtime)
                out, glitches, dirac_combs, nGlitcheses = \
                    _calcglitch(tr_slice,te_slice,ref_time,\
                    gp,interactive,plots)
            ## NEED TO STITCH TOGETHER DIFFERENT PARTS
            FILLER
        else:
            out, glitches, dirac_combs, nGlitcheses = \
                    _calcglitch(trace,template,ref_time, \
                    gp,interactive,plots)

    
    return (out, glitches, dirac_combs, nGlitcheses)
   
#=============================================================================
def _calcglitch(inp,template,ref_time, gp,interactive,plots) :
    # Calculate glitch for a given trace and (modifiable) glitch parameters
    #
    # Inputs:
    #   inp (obspy.trace):      input data trace (obspy.trace)
    #   template (obspy.trace): template of good data (same size as inp, 1s for good, 0s for bad
    #   ref_time (obspy.UTCDateTime): start time of one glitch (or a few seconds before)
    #   gp (removeGlitch.GlitchParm): glitch parameters
    #   interactive (Boolean):   ask user to verify/modify glitch parameters
    #   plots (Boolean):         plot results
     
    # Choose "slice" start time so that glitches will start no more than 1/3 of the way in:
    dp=gp.dp
    testper=gp.period
    clips=(gp.clip_lo,gp.clip_hi)
    first_glitch_offset=(ref_time-inp.stats.starttime)%testper
    if first_glitch_offset <= testper/3:
        startoffset=0
    else:
        print('='*72)
        print('\
According to the first glitch time you specified, the first glitch starts\n\
{:g} seconds ({:.0f}%) into slices based on the {:g} second glitch period.\n\
Shifting forward by {:g} seconds so that it starts 1/3 of the way in'.format(\
            first_glitch_offset,  100.*first_glitch_offset/testper,
            testper, first_glitch_offset- testper/3.) )
        print('='*72)

        startoffset=first_glitch_offset-testper/3
        newstarttime=inp.stats.starttime + startoffset
        #print(newstarttime,newstarttime-trace.stats.starttime)
  
    #======================
    ## INTERACTIVLY VERIFY GLITCH PARAMETERS
    #======================
    if interactive:
        # Set/verify clip levels
        clips=_ask_clip_levels(inp,template,startoffset,testper,clips)

        # SET/VERIFY GLITCH PERIOD
        cliptrace=inp.copy()
        cliptrace.data.clip(clips[0],clips[1], out=cliptrace.data)
        testper=_ask_test_period(cliptrace,template,startoffset,ref_time,testper)

    #======================
    # Calculate and Remove glitch 
    out,glit_1,glit_t,d_comb,nGlitch = dc.comb_clean(inp,testper,
                                                    dp,match,plots,
                                                    template,startoffset,
                                                    clips)
    out.stats.channel='CC{:d}'.format(iGlitchType)
    glit_t.stats.channel='GG{:d}'.format(iGlitchType)
    if plots:
        Stream([inp,out,glit_t]).plot()
    inp=out
    
    glitches.append(glit_1)
    dirac_combs.append(d_comb)
    nGlitcheses.append(nGlitch)
        
    return (out, glitches, dirac_combs, nGlitcheses)
#=============================================================================
def makeEQTemplate(trace,eqfile,mag_limit=5.85,dpm=1.5,verbose=True,plotit=True) :
    """
    Make a template to cut out times after large EQs
    
    Inputs:
        trace = seismological channel (obspy.trace)
        eqfile = csv file containing global earthquakes
                (column 1 = time (ISO8601), column4=magnitude)
        mag_limit: EQ Magnitude above which to cut out times [5.85]
        dpm: time (in days) to cut for each magnitude above mag_limit [1.5]
        verbose: print out list of EQs and the time that will be cut [True]
        plotit: plot traces with EQ cut out [True]
                                
    Outputs:
        template : (obspy.trace) ones for data to keep, zeros for data to remove
    """
    template=trace.copy()
    template.data=np.ones(trace.data.shape)
    template.stats.channel='TTT'

    f=open(eqfile,'r')
    for l in f:
        w=l.split(',')
        if w[0] == 'time' :
            continue
        starttime=UTCDateTime(w[0])
        mag=float(w[4])
        if starttime >= trace.stats.endtime :
            continue
        cutdays=(mag-mag_limit)*dpm
        # This is done before the magnitude limit check simply to allow
        #    printout of EQs within the acceptable time frame
        if cutdays<0:
            cutdays=0.01
        endtime=starttime+cutdays*86400
        if endtime <= trace.stats.starttime :
            continue
        print('M{:.1f} EQ of {} ...'.format(mag,w[0]),end='')
        if mag <= mag_limit :
            print('magnitude too small to worry about')
            continue
        print('blanking {:.1f} days'.format(cutdays))
        startaddr=int((starttime-trace.stats.starttime)*sps)
        if startaddr<0:
            startaddr=0
        if endtime > trace.stats.endtime :
            endtime = trace.stats.endtime
        endaddr = int((endtime-trace.stats.starttime)*sps)
        #template.data[startaddr:endaddr]=0
        template.data[startaddr:endaddr]=0
        
    # PLOT THE RESULT OF REMOVING BIG EQs
    if plotit:
        temp=template.copy()
        temp.data=temp.data*trace.data
        a=obspy.core.stream.Stream([trace,temp]).plot(color='blue')
        
    return template
    
#=============================================================================
#=============================================================================
def _ask_clip_levels(inp,template,st_off,testper,clip) :
  """
  Show clip levels and ask to update them until acceptable
  
  inp (obspy.trace) : seismological trace
  template (obspy.trace): signal of ones and zeros, zeros where trace is to be hidden
  st_off: offset from start of trace to start slicing
  testper: length of each slice (seconds)
  clip: 2-tuple of default clip values (lo, hi)
  """
  st=inp.stats.starttime
  sps=inp.stats.sampling_rate
  station=inp.stats.station
  stack_trace=inp.copy()
  stack_trace.data=stack_trace.data*template.data
  if st_off>0:
    stack_trace=stack_trace.slice(starttime=st + st_off)
  stack=dc.stack_data(stack_trace.data,testper*sps)
  nrows,ncols=stack.shape
  time=np.arange(nrows)/sps
  slicenums=np.arange(ncols)
  title='{} sliced at {:g}s, stacked'.format(station,testper)
  # Show clip levels and verify that they are ok
  while True :
    plt.plot(time,stack,linewidth=0.1)
    plt.plot([time[0],time[-1]],[clip[0],clip[0]],'k--',
                        label='clip_lo = {:g}'.format(clip[0]))
    plt.plot([time[0],time[-1]],[clip[1],clip[1]],'k-.',
                        label='clip_hi = {:g}'.format(clip[1]))
    plt.xlabel('Time (seconds)')
    clip_range=clip[1]-clip[0]
    plt.ylim((clip[0]-clip_range*.3,clip[1]+clip_range*.3))
    plt.legend()
    plt.title(title)
    plt.show()
    ## Ask for a new clip_lo,clip_hi tuple, continue if current value accepted
    newval=dc.input_float_tuple('Enter clip levels containing all glitches (RETURN to accept current value)',clip)
    if newval==clip:
      break
    else:
      clip=newval
  return clip


#======================================================================
def _ask_test_period(inp,template,st_off,ref_time,testper) :
  """
  Show glitch alignment and ask to update glitch period until acceptable
  
  Also allows to verify that the glitch_parm.starttime is ok, but not to
  modify it
  
  inp (obspy.trace) : seismological trace
  template (obspy.trace): signal of ones and zeros, zeros where trace is to be hidden
  st_off: offset from start of trace to start slicing
  ref_time (UTCDateTime) : time at start of (or just before) one of the glitches
  testper: length of each slice (seconds)
  """

  st=inp.stats.starttime
  sps=inp.stats.sampling_rate
  station=inp.stats.station
  stack_trace=inp.copy()
  stack_trace.data=stack_trace.data*template.data
  if st_off>0:
    stack_trace=stack_trace.slice(starttime=st + st_off)
  while True :
    stack=dc.stack_data(stack_trace.data,testper*sps)
    nrows,ncols=stack.shape
    time=np.arange(nrows)/sps
    slicenums=np.arange(ncols)
    if ref_time:
      ref_offset=(ref_time-stack_trace.stats.starttime)%testper
    else:
      ref_offset=0
    title='{} sliced at {:g}s, clipped, stacked'.format(station,testper)

    # Plot as stacked traces
#     plt.figure(1)
#     plt.plot(time,stack)
#     plt.plot([ref_offset,ref_offset],[stack.min().min(),stack.max().max()],
#             'k--')
#     plt.xlabel('Time (seconds)')
#     plt.title(title)

    # plot as pcolor 
    plt.figure(1)
    plt.pcolormesh(slicenums,time,stack)
    plt.plot([0,ncols],[ref_offset,ref_offset],'k--')
    plt.xlabel('Period slice #')
    plt.ylabel('Time (seconds)')
    plt.title(title)
    plt.grid(True,zorder=20)
    plt.show()
  
    # Ask for new test period, continue if current value accepted
    newval=dc.input_float(
            'Enter new test period (RETURN to accept current value)',
            testper)
    if newval==testper:
      break
    else:
      testper=newval
  return testper

#########################################################################
# if __name__ == "__main__":
# 	sys.exit(main())
