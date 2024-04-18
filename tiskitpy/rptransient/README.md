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

### Prepare decimated and clean-rotated data:
```python
from obspy import UTCDateTime
from obspy.clients.filesystem.sds import Client
from rptransient import decimate, rotate_clean

net = 'XX'
sta = 'I04D'
archive_dir = '/Volumes/Wayne_Data/_NoBackup/SDS_PILAB_uncorrected'

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

### Calculate and remove transients
```python
from obspy.core import read
from rptransient import get_eq_spans, Transients, PeriodicTransient as PT

sta = 'I04D'
datafile = f'Data/{sta}_20161202_20161222_dec_rot.mseed'

transients = {'I04D': [PT("1h", 3620.26, 0.05, [-250, 140], '2016-12-02T00:46:00')],
              'I12D': [PT("1h", 3619.76, 0.05, [-320, 250], '2016-12-02T00:38:00')],
              'I14D': [PT("1h", 3619.73, 0.05, [-250, 180], '2016-12-02T00:50:00')],
              'I28D': [PT("1h", 3619.71, 0.05, [-450, 250], '2016-12-02T00:13:00')],
              'I34D': [PT("1h", 3619.67, 0.05, [-350, 400], '2016-12-02T00:11:00')]}
plot = True

rt = Transients(transients[sta])
stream = read(datafile,'MSEED')
zdata = stream.select(channel='*Z')[0]
eq_spans = get_eq_spans(zdata.stats.starttime, zdata.stats.endtime)
rt.calc_timing(zdata, eq_spans)   
rt.calc_transients(zdata, eq_spans, plot=plot)
cleaned = rt.remove_transients(zdata, plot=plot, match=False, prep_filter=True)
cleaned.write(datafile.split('.')[0] + '_cleaned.mseed', format='MSEED')
```

- `calc_timing()` plots two screens with parameters that you will have to
  validate or change on the command line.
  - `Select clip levels` presents the stacked glitches (using the glitch
    period you provided in `transients`) and the clip levels you provided in
    `transients`.  These should cover the min and max glitch levels, but not
    more.  If they do, hit <RETURN> on the command line.  If they don't, enter
    the clip levels, separated by a comma, on the command line.  The plot will
    be redisplayed with the new levels and you can validate or readjust, and
    so on (if you changed the clip levels, you should enter the new levels
    in `transients`, so you don't have to re-change next time)
  - `Select transient period` presents the stacked glitches and a dotted line
    that represents the glitch start time.  If the glitches are not aligned with
    the dotted line, then enter a new transient period on the command line, if
    they are, hit return.  Once the glitches are aligned, they should start just
    after (above) the dotted line.  If they aren't, you will need to change
    the start time in `transients` and rerun

- `calc_transients()` calculates the transients.  If `plot=True`, it will
  plot a few screens that you will need to close once you've looked at them.
  
- `remove_transients(plot=True)` removes the transients from the data.
  If `plot=True`, it plots a final screen comparing the signal
  before and after the correction, plus the correction itself. 

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

```python
def get_eq_spans(starttime, endtime, minmag=5.85,
                 days_per_magnitude=1.5, quiet=False):
    """
    Args:
        starttime (UTCDateTime): earliest data that will be presented
        endtime (UTCDateTime): latest data that will be presented
        minmag (float): EQ Magnitude above which to cut out times
        days_per_magnitude (float): days to cut per magnitude above
            min_magnitude
    Returns:
        eq_spans (~class `.TimeSpans`): time spans covering EQ signal
    """
```
Classes:
----------------------

```python
class Transients():
    """
    Periodic Transients class
    """
    def __init__(self, transients=[]):
        """
        Args:
            transients (list of PeriodicTransient): transients
        """

    def calc_timing(self, trace, eq_remover, prep_filter=True):
        """
        Calculate transient time parameters

        Args:
            trace (:class ~obspy.core.Trace): data
            eq_remover (:class ~EQRemover): periods to remove because of EQs
            prep_filter (bool): apply Wiedland prep filter before processing?
        """
        
    def calc_transients(self, trace, bad_times, plot=False,
                        prep_filter=True):
        """
        Calculate transients

        Args:
            trace (~class `obspy.core.trace.Trace`): data
            bad_times (~class `.TimeSpans`): time spans to zero out
            plot (bool): make information plots
            prep_filter (bool): apply Wieland prep filter before processing?
       """
        
    def remove_transients(self, trace, match=True, plot=False,
                          prep_filter=False):
        """
        Remove transient from data

        Args:
            trace (~class `obspy.core.trace.Trace`): data
            match (bool): individually match transients to data?
            plot (bool): make information plots
            prep_filter (bool): apply Wieland prep filter before processing?
        """
```

```python
class PeriodicTransient():
    """
    Information about one periodic transient
    """
    def __init__(self, name: str, period: float, dp: float, clips: tuple,
                 transient_starttime: UTCDateTime):
        """
        Args:
            name (str): name of this periodic transient (e.g., 'hourly')
            period (float): seconds between each transient
            dp (float): how many seconds to change the period by when testing for
                   better values
            clips (tuple): clip values outside of this range (low, high).
                Should contain the max range of the transient
            transient_starttime (~class `obspy.core.UTCDateTime`): onset
                time of earliest transient.
        
        The program will make transient slices starting between
        transient_starttime and 1/3 of self.period earlier
        
        transient_starttime must not be too late, a few seconds early is ok
        """
```
