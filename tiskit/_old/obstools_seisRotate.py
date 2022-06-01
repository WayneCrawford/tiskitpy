#!env python3
""" 
Module to clean tilt noise from OBS vertical channel through
non-deforming rotation
"""

import sys
import math as M
import numpy as np
import scipy as sp
import matplotlib as mp
import obspy
from obspy import UTCDateTime
import matplotlib.pyplot as plt
from obspy.signal import rotate
from copy import deepcopy

DEBUG = True

def debug(s):
    if DEBUG:
        print(s)

#============================================================================
def cleanTiltNoise(stream,excludes,horiz_too=False,plotit=False) :
    """
        Minimizes tilt on vertical channel by rotating to minimize noise
    
        Usage: cleanTiltNoise(stream,exclude)
    
        Inputs:
            stream (obspy.stream) input data, must have a *Z, *[1|N] and *[2|E] channel
            excludes: list of dictionaries containing time periods to avoid in "start"
                     and "end" keys
            horiz_too (boolean): rotate horizontals also (True if you think
                        the channels are truly orthogonal) [False]
            plotit (boolean): Plot comparision of original and rotated vertical [False]
                
        Outputs:
            stream_rot (obspy.stream): rotated stream
            angles (tuple), angle,azimuth by which Z (or all three) were rotated
    """
    ##################################################################################
    # CONSTANTS
    ##################################################################################    
    quickTest=False     # Only run one day's data (will not save result)

  
    ################################################################################
    # Filter data to tilt noise band for best Angle calc)
    ################################################################################
    filtstream=stream.copy()
    for interval in excludes:
        filtstream.cutout( UTCDateTime(interval['start']),
                            UTCDateTime(interval['end'])   )
    filtstream.detrend('demean')
    filtstream.detrend('linear')
    filtstream.filter('lowpass',freq=0.01,corners=4,zerophase=True)
    filtstream.filter('highpass',freq=0.001,corners=4,zerophase=True)

    print('Creating seisRotate object')
    srData=SeisRotate(filtstream)   
    #srData=seisRotate.SeisRotate(filtstream)   

    print('Calculating rotation angle')
    (bestAngle,bestAzimuth)=srData.calc_zrotate_opt(verbose=True)
    print('Best angle,azimuth is ({:.2f},{:.2f}).  Rotating axes'.format(
                                                bestAngle,bestAzimuth))

    ################################################################################
    # Filter data (to make it easier to verify that rotation works on low freq noise)
    ################################################################################
    viewstream=stream.copy()
    print('Filtering data for viewing')
    viewstream.filter('lowpass',freq=0.005,corners=4,zerophase=True)
    viewstream.filter('highpass',freq=0.001,corners=4,zerophase=True)

    srData=SeisRotate(stream)   
    viewData=SeisRotate(viewstream)   
    #srData=seisRotate.SeisRotate(stream)   
    #viewData=seisRotate.SeisRotate(viewstream)   

    srData.zrotate(bestAngle,bestAzimuth,horiz_too)
    viewData.zrotate(bestAngle,bestAzimuth)

    print('Converting back to Stream')
    strm_rot=srData.stream()
    view_rot=viewData.stream()

    ##############################################################################
    # PLOT RESULTS (Z channels)
    ##############################################################################
    trace_view_Z=viewstream.select(component='Z')[0]
    trace_view_rot_Z=view_rot.select(component='Z')[0]
    trace_view_rot_Z.stats.channel=trace_view_rot_Z.stats.channel[0]+ \
                                   'R'+trace_view_rot_Z.stats.channel[2]
    compare_stream=obspy.core.stream.Stream([trace_view_Z,trace_view_rot_Z])
    # PLOT TWO DAYS
    #compare_stream.plot(equal_scale=True,starttime=trace_view_Z.stats.starttime,
    #                     endtime=trace_view_Z.stats.starttime+2*86400,method='full')
    if plotit:
        # PLOT ALL
        compare_stream.plot(equal_scale=True,method='full')

    return strm_rot,(bestAngle,bestAzimuth)


#============================================================================
class SeisRotate:
  """ 
  Class to clean tilt noise from OBS vertical channel through
  non-deforming rotation
  """
  def __init__(self,stream):
    """
    Create a seisRotate object from 3-component obsPy.stream.Stream data
    
    Channel names must end in Z, N and E or Z, 1, and 2 
    Z is up, 1 and 2 are horizontal orthogonal with 2 90° clockwise
    of "1"  (in other words, "1" corresponds to "N" and "2" to "E" , 
    except that they are not necessarily aligned with geographic
    cardinals)
    """
    # Read in data
    try:
      self.Z=stream.select(component='Z')[0].copy()
    except:
      print('No Z channel found')
      sys.exit(2)
    try :
      self.N=stream.select(component="N")[0].copy()
    except:
      try :
        self.N=stream.select(component="1")[0].copy()
      except:
        print('''No Z or "1" channel found''')
        sys.exit(2)
    try :
      self.E=stream.select(component="E")[0].copy()
    except:
      try :
        self.E=stream.select(component="2")[0].copy()
      except:
        print('''No E or "2" channel found''')
        sys.exit(2)
        
    # Verify all channels have same length
    if len(self.Z.data) != len(self.N.data):
      print('Z and N vectors have different lengths: {:d} & {:d}'.format(
          len(self.Z.data),len(self.N.data)))
      sys.exit(2)
    if len(self.Z.data) != len(self.E.data):
      print('Z and E vectors have different lengths')
      sys.exit(2)
      
    # Verify all channels have same sampling_rate
    self.fs=self.Z.stats.sampling_rate
    if (self.fs != self.E.stats.sampling_rate) or \
                  (self.fs != self.N.stats.sampling_rate) :
      print('Channels do not have the same sampling rate')
      sys.exit(2)
    
  def __len__(self):
    return len(self.Z.data)
    
  def plot(self):
    self.stream().plot()
    
  def copy(self):
        """
        Returns a copy of the seisRotate object.

        :return: Copy of seisRotate object.

        """
        #print(type(self.Z))
        return deepcopy(self)
        
  def stream(self):
    """ Return obspy stream() """    
    return obspy.core.stream.Stream(traces=[self.Z,self.N,self.E])
  
  def zrotate(self,angle,azimuth,horiz_also=True):
    """
    Rotate a seismometer's vertical component
    
    :param angle: angle (towards vertical) by which to rotate (degrees)
    :param azimuth: azimuth (CW from N) at which to rotate (degrees)
    :param horiz_also: rotate horizontal components as well, to retain
                       orthogonality
    :type  horiz_also: boolean
    
    I would set horiz_also to:
      False for the Guralp CMG-3T, which independently levels the three
               components
      True for the Trillium T240, which centers without rotation
               (according to Dr Wieland)
    """
    
    dip_Z=-90+angle
    #if not horiz_also:
    #  dip_N=0
    #  dip_E=0
    #else:
    # Calculated rotations for N and E using following logic
    # if azimuth=0 (aligned North), then Ndip=vertang, Edip=0
    # if azimuth=90 (aligned East), then Ndip=0, Edip=vertang
    # if azimuth=180 (aligned South), then Ndip=-vertang, Edip=0
    # if azimuth=270 (aligned West), then Ndip=0, Edip=-vertang
    #
    # Fits Ndip=vertang*cos(azimuth), Edip=vertang*sin(azimuth)
    dip_N=angle*M.cos(np.deg2rad(azimuth))
    dip_E=angle*M.sin(np.deg2rad(azimuth))
    
    [self.Z.data,N,E]=obspy.signal.rotate.rotate2zne( \
                              self.Z.data,azimuth, dip_Z, \
                              self.N.data, 0,      dip_N, \
                              self.E.data, 90,     dip_E)
    if horiz_also:
      self.N.data=N
      self.E.data=E
      
  def calc_zrotate_opt(self,lowcut=0.001,hicut=0.005,verbose=False):
    """
    Calculate the Z channel rotation angle that minimizes tilt noise
    
    :param lowfreq: low passband frequency in which to lood for energy
                   reduction
    :param hifreq:  high passpand frequency "" "" ""
    :param verbose: output information about the angles tested
    
    The default values (0.001, 0.005) correspond to the band where
    removing tilt noise generally has the greatest effect on the Z
    signal (at the low frequency end of the noise notch)
    
    You can input the returned angles to rotateSeis.zrotate() to get
    the properly rotated data
    """
    # Filter data, cut off edge effects and detrend
    
    filt=self.copy() 
    filt.Z.filter('bandpass',freqmin=lowcut,freqmax=hicut,corners=5)
    filt.N.filter('bandpass',freqmin=lowcut,freqmax=hicut,corners=5)
    filt.E.filter('bandpass',freqmin=lowcut,freqmax=hicut,corners=5)

    # Quick estimate of angles using Z/E and Z/N ratios. DOESN'T HELP: SKIP
    #(startAngle,startAzimuth)=filt._estimateAngles(verbose)
    startAngle,startAzimuth=(0,0)
    #if verbose:
    #  print("start angle, azimuth = {:.2f},{:.2f}".format(
    #                     startAngle,startAzimuth)
    
    # Search around quick estimated angles for the best angles
    (angle,azimuth)=filt._searchBestAngles(startAngle,startAzimuth,
                                           verbose)
    
    return (angle,azimuth)
    
  def _estimateAngles (self,verbose=False):
    """
    Estimate how far and in which direction Z is tilted from vertical
    by comparing with N and E
    
    azimuth is 0 towards N, 90 towards E
    dip is w.r.t. vertical
    
    Assumes data has already been filtered into the relevant band

    Doesn't seem to work AT ALL! I think there is too much other stuff
    (glitches, etc) for this to work.  Also doesn't account for any time
    lags: might work better if ratios calculated from cross-corellation
    
    :param verbose: display extra information
    """
    
    if verbose:
      print("Estimating preliminary angle based on signal ratios")
    #s=obspy.core.stream.Stream([self.Z,self.N])
    #s.plot()
    ZoverN= np.divide(self.Z.data,self.N.data)
    ZoverE= np.divide(self.Z.data,self.E.data)
    # Throw out not-a-numbers
    ZoverN=ZoverN[~np.isnan(ZoverN)]
    ZoverE=ZoverE[~np.isnan(ZoverE)]
    #if verbose:
    #  print("    ",np.median(ZoverN),np.median(ZoverE))
    
    ZNratio=M.copysign(M.sqrt(np.median(ZoverN**2)),np.median(ZoverN))
    ZEratio=M.copysign(M.sqrt(np.median(ZoverE**2)),np.median(ZoverE))
    ZNratio=np.median(ZoverN)
    ZEratio=np.median(ZoverE)
    plt.plot(ZoverN,'x')
    
    if verbose:
      print('    ZNratio = {:<10g}, ZEratio = {:<10g}'.format(ZNratio,ZEratio))
    
    ZNangle=(180/M.pi)*M.asin(ZNratio)
    ZEangle=(180/M.pi)*M.asin(ZEratio)
    azimuth=(180/M.pi)*M.atan(ZNratio/ZEratio)
    azimuth=90-(180/M.pi)*M.atan2(ZNratio,ZEratio)
    if verbose:
      print('    ZNangle = {:<10g}, ZEangle = {:<10g}'.format(ZNangle,ZEangle))
    # The angle from vertical below is just a guess, should developed  
    # mathematically (or at least derived empirically from the searched
    # results)
    angle=M.sqrt(ZNangle**2 + ZEangle**2)
    if verbose:
      print('    estAngle= {:<10g}, estAzim = {:<10g}'.format(angle,azimuth))

    # Sanity check: do the dip equations from zrotate() return ZN and ZE angles??
    dip_N=angle*M.cos(np.deg2rad(azimuth))
    dip_E=angle*M.sin(np.deg2rad(azimuth))
    print('dip_N= {:<10g}, dip_E={:<10g}'.format(dip_N,dip_E))
    
    return angle,azimuth
  
  def _searchBestAngles (self,startAngle=0,startAzimuth=0,
                         verbose=False):
    """
    Find best Z rotation angles
    
    :param startAngle: starting guess for the angle (degrees)
    :param startAzimuth: starting guess for the azimuth (degrees)
    :param verbose: display extra information
    :rtype: (angle, azimuth)
  
    Search for minimum Z energy as function of angle
    """
    if verbose:
      print("Calculating best angle based on variance minimization",end=' ')
      print("(uses sp.optimize.fmin)")
    #s=obspy.core.stream.Stream([self.Z,self.N])
    #s.plot()
    start_var=self._rotZ_variance([startAngle,startAzimuth])
    
    xopt,fopt,iter,funcalls,warnflag=sp.optimize.fmin( \
                                func=self._rotZ_variance,
                                x0=[startAngle,startAzimuth],
                                disp=verbose, 
                                full_output=True,retall=False)
    bestAngles=xopt
    if verbose:
      print(xopt,fopt,iter,funcalls)
      print("    best Angle,Azimuth = {:.2f},{:.2f}".format( \
                                        bestAngles[0],bestAngles[1]))
      print("    variance reduced from {:.2e} to {:.2e} ({:.1f}% lower)".format( \
                                        start_var,fopt,100*(1-fopt/start_var)))
    return (bestAngles[0],bestAngles[1])
    
    
  def _rotZ_variance(self,angles) :
    """
    Calculate the variance for a given rotation
    
    :param angles: 2 element list: angle and azimuth
    :type angles: list

    Assumes data are already filtered into relevant band
    """
    A=self.copy()
    #A.Z.stats.channel='XXX'
    A.zrotate(angles[0],angles[1])
    var=np.sum(A.Z.data**2)
    #print(angles[0],angles[1],var,A.Z.data[0:5],A.Z.data[0:5]**2,sum(A.Z.data[0:5]**2))
    return var
    
    
    