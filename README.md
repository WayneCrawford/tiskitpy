remove_regular_transients
=========================

Remove regularly-spaced transients from OBS data

Based on Matlab code by E Wielandt, used in Deen et al., 2017

Overview
======================

Example run:
----------------------

Prepare decimated and clean-rotated data:
```python
from obspy.clients.filesystem.sds import Client
from rrtrans import decimate, rotate_clean

net = 'YV'
sta = 'I04D'
noisys = []

# READ DATA, DECIMATE AND ROTATE
starttime = UTCDateTime('2012-12-01T00:00:00')
endtime = starttime + 86400
decimates = [5, 5, 5]
sds = Client('/SDS')

stream = sds.get_waveforms(net, sta, "", "??Z", starttime, endtime)

dec_data = decimate(stream, decimates)
rot_data = rotate_clean(dec_data, noisys)
rot_fname = f'{station}_dec.mseed'
rot_data.write(rot_fname, 'MSEED', encoding='STEIM1')
```

Calculate and remove transients
```python
from obspy.core import read
from rrtrans import Transients, PeriodicTransient as PT

station = 'I04D'
eqfile = 'earthquakes_201210-201310.csv'
transients = [PT("1h", 3620.27, 0.05, [-300, 200], ['2016-04-01T00:28:30']),
             [PT("3h", 9519.14, 0.2,  []-110,  80], ['2016-04-01T00:00:00'])]

stream = read(f'{station}_dec.mseed','MSEED')
rt = Transient(transients)
trans_times = rt.calc_times(rot_data.z, eq_file, delay)   # Interactive tuning of transient paramters
trans = rt.calc_trans(rot_data, trans_times, eq_file, delay)
cleaned = rt.clean_trans(rot_data, trans, trans_times)
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
rotate_clean(stream, excludes, horiz_too=False, plotit=False,
             quickTest=False):
    """
    Rotates vertical channel to minimize noise

    :param stream: input data, must have a *Z, *[1|N] and *[2|E] channel
    :type stream: ~class obspy.core.stream.Stream
    :param excludes: list of dictionaries containing time periods to avoid,
                     using "start" and "end" keys
    :param horiz_too: rotate horizontals also (use True if you think
                      the channels are truly orthogonal)
    :param plotit: Plot comparision of original and rotated vertical
    :param quickTest: Only run one day's data and do not save results

    :returns: stream_rot (obspy.stream): rotated stream
              (angle, azimuth) by which Z (or Z-X-Y) were rotated
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