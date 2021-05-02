#!env python3
""" Resample, rotate to reduce tilt, and remove glitch"""

#################################################
import obspy
from obspy import UTCDateTime
from obspy.core import stream as Stream
import time
# The following are my own routines
import downSample as ds
import seisRotate as sr
import removeGlitch as rg
#import sys
    
##############################################################################
# Base parameters
##############################################################################
# Station problems: 
#   I07D has just a couple weeks of data (2016-03-08T13 to 2016-03-22T02); 
#   I20D has only BDH data
station='I07D'
interactiveGlitch=True  # Useful for tuning glitch parameters
starttime='2016-04-01T00:00:00.0Z'
endtime='2016-04-10T00:00:00.0Z'
starttime='2016-03-08T00:13:00.0Z'
endtime='2016-03-22T02:00:00.0Z'
decimates=[5,5,5]

baseDir="/Volumes/SISMOSMOOdata/Data_PiLAB/{}/".format(station)
filename="PILAB_XS_{}_00_*_nodriftcorr.mseed".format(station)
#filename="PILAB_XS_I14D_00_BDH_2016-070-130000_nodriftcorr.mseed"

##############################################################################
# Glitch parameters
##############################################################################
dp_1hr=0.05
dp_3hr=0.2
dp_24hr=1.
# 'glitch' is used by removeGlitch, provides glitch parameters
# 'exclude' is used by removeGlitch and cleanTiltNoise, provides time periods
#           with bad/noisy data that we don't want to use in the glitch calculation
#           or the tilt removal
station_parms=\
    {
    'I04D': 
        { 'glitch':
            [
            rg.GlitchParm(3620.27, dp_1hr,  -300, 200,'2016-04-01T00:28:30') ,
            rg.GlitchParm(9519.14, dp_3hr,  -110,  80,'2016-04-01T00:00:00') 
            ],
          'exclude':
            [
            ]
        },
    'I07D':
        { 'glitch':
            [
            rg.GlitchParm(3620.26, dp_1hr,  -700, 800,['2016-03-08T13:22:20','2016-03-10T17:00:00']) ,
            rg.GlitchParm(9519.1,  dp_3hr,   -60,  60,'2016-04-01T00:00:00') 
            ],
          'exclude':
            [
                {'start':'2012-12-01T00:00:00','end':'2012-12-12T14:35:00'},
                {'start':'2012-12-26T10:00:00','end':'2012-12-26T18:00:00'} 
            ]
        },
    'I12D':
        { 'glitch':
            [
            rg.GlitchParm(3619.73, dp_1hr,  -300, 300,'2016-04-01T00:00:00') ,
            ],
          'exclude':
            [
                {'start':'2012-12-01T00:00:00','end':'2012-12-12T14:35:00'},
                {'start':'2012-12-26T10:00:00','end':'2012-12-26T18:00:00'} 
            ]
        },
    'I14D':
        { 'glitch':
            [
            rg.GlitchParm(3619.73, dp_1hr,  -300, 300,'2016-04-01T00:00:00') ,
            ],
          'exclude':
            [
                {'start':'2012-12-01T00:00:00','end':'2012-12-12T14:35:00'},
                {'start':'2012-12-26T10:00:00','end':'2012-12-26T18:00:00'} 
            ]
        }
    }

##############################################################################
# Read data
##############################################################################
starttime=UTCDateTime(starttime)
endtime=UTCDateTime(endtime)
print('Reading {} from {} to {}... '.format(\
            filename,starttime.isoformat(),endtime.isoformat()),
            end='',flush=True)
tic=time.time()
# Use routine in downSample that allows read of miniSEED files > 2.1 GB
stream=ds.readMSEED(baseDir+filename,starttime=starttime,endtime=endtime)
print('Took {:.1f} seconds'.format(time.time()-tic))
print(stream)
print('Original data has {} samples'.format([tr.data.size for tr in stream]))
#stream.plot()
origdtype=stream[0].data.dtype

##############################################################################
# Decimate data
##############################################################################

print("Decimating data")
stream_decim = ds.decimate(stream,decimates)
stream_decim.verify()

##############################################################################
# Rotate vertical to minimize noise
##############################################################################

print("Rotating vertical")
stream_decim_rot,angles = sr.cleanTiltNoise(stream_decim,station_parms[station]['exclude'])
stream_decim_rot.verify()

##############################################################################
# Calculate and remove glitch(es)
##############################################################################

print("Removing glitches")
trace=stream.select(channel='*HZ')[0]
traceZ_deglitch,glitches,dirac_combs,nGlitcheses=rg.removeGlitch(\
                        stream_decim_rot,
                        station_parms[station]['glitch'],
                        eqfile='earthquakes_201210-201310.csv',
                        interactive=interactiveGlitch)
# can do the same for the horizontals using values from vertical

stream_deglitch=stream_decim_rot.copy()
stream_deglitch.select(channel='*HZ')[0]=traceZ_deglitch
# Change instrument code to "X" : Derived or Generate Channel
old_ch_name=stream_deglitch.select(channel='*HZ')[0].stats.channel
new_ch_name=old_ch_name[0]+'X'+old_ch_name[2]
stream_deglitch.select(channel='*HZ')[0].stats.channel=new_ch_name


##################################################################################
# Save result
##################################################################################
net=stream_deglitch[0].stats.network
# Convert traces back to original data type
print('Converting data back to {} format... '.format(origdtype),end='',flush=True)
tic=time.time()
for trace in stream_deglitch:
  trace.data=trace.data.astype(origdtype)
print('Took {:.1f} seconds'.format(time.time()-tic))

# Save data to single miniSEED file
outfname='./{}_{}_{}-{}_{}.mseed'.format(net,sta,new_ch_name,startdate.strftime('%Y%m%d'))
print('Saving downsampled, rotated, deglitched data to {}... '.format(outfname))
stream_deglitch.write(outfname,'MSEED')
# Save correction to another?