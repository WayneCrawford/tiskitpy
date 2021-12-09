rptransient
=========================

Remove periodic transient(s) from seismological data

Based on Matlab code by E Wielandt, used in Deen et al., 2017

Installation
=========================

Clone or download this repository, then from within the main repository directory, run:

```pip install .```

You can also install in editable mode (for developers), with:

```pip install -e .```

Overview
======================

Example run:
----------------------

Prepare decimated and clean-rotated data:
```python
from obspy import UTCDateTime
from obspy.clients.filesystem.sds import Client
from rptransient import decimate, rotate_clean

net = 'XX'
sta = 'I04D'
archive_dir = '/Volumes/Wayne_Data/_NoBackup/SDS_PILAB'
archive_dir = '/Users/crawford/_Work/Parc_OBS/7_Missions/2016.PiLAB/7_DataExtraction/2021.lc2ms/test/SDS'
eqfile = 'eqs_20160227-20170118.csv'

# READ DATA, DECIMATE AND ROTATE
starttime = UTCDateTime('2016-12-01T00:00:00')
endtime = starttime + 10*86400
decimates = [5, 5, 5]
sds = Client(archive_dir)

stream = sds.get_waveforms(net, sta, "*", "*", starttime, endtime)
# stream.merge()
print(endtime)
print(stream.__str__(extended=True))

dec_data = decimate(stream, decimates)
rot_data, ang, azi = rotate_clean(dec_data, remove_eq=True, plot=True)
print(f'{sta} Z rotated by {ang:.2f} degrees at azimuth {azi:.1f}')
time_str = f'{starttime.strftime("%Y%m%d")}_{endtime.strftime("%Y%m%d")}'
rot_fname = f'{sta}_dec_{time_str}.mseed'
rot_data.write(rot_fname, 'MSEED')
```

Calculate and remove transients
```python
from obspy.core import read
from rrtrans import Transients, PeriodicTransient as PT

station = 'I04D'
eqfile = 'earthquakes_201210-201310.csv'
transients = [PT("1h", 3620.27, 0.05, [-300, 200], ['2016-04-01T00:28:30']),
             [PT("3h", 9519.14, 0.2,  []-110,  80], ['2016-04-01T00:00:00'])]

rt = Transients(transients)
stream = read(f'{station}_dec.mseed','MSEED')
zdata = stream.select(channel='*Z')[0]
eq_template = EQTemplate(zdata, eqfile)
# Interactive tuning of transient parameters
rt.calc_timing(zdata, eq_template)   
trans = rt.calc_transient(zdata, eq_template)
cleaned = rt.clean_transient(zdata, trans)
cleaned.write(format='MSEED')
```

Functions:
----------------------
```python
decimate(stream, decimates, verbose=False):
    """
    Decimate obspy stream in some intelligent way

    Inputs:
        stream (obspy.stream)
        decimates: list of decimation factors (use 2,3 or 5)

    Outputs:
        stream_decim (obspy.stream) decimated data
    """
```

```python
def rotate_clean(stream, excludes=[], horiz_too=False, plot=False,
                 quickTest=False, remove_eq=True, verbose=True):
    """
    Rotates vertical channel to minimize noise
    Arguments:
        stream (Stream): input data, must have a *Z, *[1|N] and *[2|E] channel
        remove_eq (str, True or False): filename of catalog to use to remove
            earthquakes, will download catalog from USGS if True, not remove
            EQs if False
        excludes: list of dictionaries containing time periods to avoid,
                     using "start" and "end" keys
        horiz_too: rotate horizontals also (use True if you think
                      the channels are truly orthogonal)
        plot: Plot comparision of original and rotated vertical
        quickTest: Only run one day's data and do not save results
    Returns:
        (tuple): 3-tuple containing:
            (Stream): rotated stream
            (float): angle by which Z (or Z-X-Y) was rotated
            (float): azimuth by which Z (or Z-X-Y) was rotated
    """
```

Classes:
----------------------

```python
class Transients():
    """
    Periodic Transients class
    """
    def __init__(trans_parms=[]]):
        """
        :param trans_parms: list of transient parameters
        """
    def calc_times(stream, catalog=None):
        """
        Calculate transient times

        :param stream: data (trace?)
        :param catalog: earthquake catalog (to blank noisy times) 
        """
        
    def calc_trans(stream, times=[], catalog=None):
        """
        Return transient model

        :param stream: data (trace?)
        :param times: list of transient times
        :param catalog: earthquake catalog (to blank noisy times) 
        """
        
    def clean_trans(stream, trans_model, times):
        """
        Remove transient from data

        :param stream: data (trace?)
        :param trans_model: transient_model
        :param times: list of transient times
        :param catalog: earthquake catalog (to blank noisy times) 
        """
```

```python
class TransParms():
    """
    Information about one type of periodic transient
    """
    def __init__(self, name: str, period: float, dp: float, clips: tuple,
                 starttimes=[]):
        """
        :param name: name of this periodic transient (e.g., 'hourly')
        :param period: seconds between transients (seconds)
        :param dp: how many seconds to change the period by when testing for
                   better values
        :param clips: clip values outside of this range.  Should correspond to
                      max range of glitch
        :type clips: 2-tuple (lowest, highest)
        :param starttimes: onset times of earliest glitch(es).  Multiple values
                           allow for a change due to releveling, etc. It's
                           better to be a few seconds too early than too late.
                           The program will make transient slices starting
                           between a starttime and 1/3 of a transient period
                           earlier
        :type starttimes:  list of ~class `obspy.core.UTCDateTime`
        """

    def verify(self, trace, template, ref_time):
        """
        Interactively verify transient parameters
    
        :param trace:      input data trace
        :type trace: ~class `obspy.stream.Trace`
        :param template: template of good data (same size as inp, 1s for good,
                         0s for bad
        :type template: ~class `obspy.stream.Trace`
        :param ref_time: start time of first glitch (or a few seconds before)
        :type ref_time: ~class `obspy.core.UTCDateTime
        """
```

Other Useful Functions:
----------------------

```
decimate(): Decimates data in the most non-modifying fashion
```

```
rotate_clean(): Rotate the vertical axis to minimize current noise
```

```
plot_EQ_template(trace, catalog, mag_limit, dpm):
    """
    plot data before and after applying EQ masking template with provided
    maglimit and dpm
    """
```
