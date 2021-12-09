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
from rptransient import EQRemover, Transients, PeriodicTransient as PT

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
eq_remover = EQRemover(zdata.stats.starttime, zdata.stats.endtime)
rt.calc_timing(zdata, eq_remover)   
rt.calc_transients(zdata, eq_remover, plot=plot)
cleaned = rt.remove_transients(zdata, plot=plot, match=False, prep_filter=True)
cleaned.write(datafile.split('.')[0] + '_cleaned.mseed', format='MSEED')
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

```python
class EQRemover:
    """
    A class for zeroing out data after big earthquakes

    Arguments:
        start_time (UTCDateTime): earliest data that will be presented
        end_time (UTCDateTime): latest data that will be presented
        min_magnitude (float): EQ Magnitude above which to cut out times
        days_per_magnitude (float): days to cut per magnitude above
            min_magnitude
    """

    def eq_filename(self, starttime, endtime, min_magnitude, format='quakeml'):
        """
        Create earthquake file filename
        
        Args:
            starttime (UTCDateTime): start time
            endtime (UTCDateTime): end time
            min_magmitude (float): minimum magnitude saved
            format (str): file format (only 'quakeml')
        Returns:
            filename (string):
        """

    def add_timespans(self, time_spans):
        """
        Add time spans to existing EQRemover object
        
        Args:
            time_spans (~class `.TimeSpans`): time spans to add
        """

    def has_zeros(self, starttime, endtime):
        """
        Does the trace's time span intersect an EQ zero?

        Arguments:
            starttime (UTCDateTime): start time
            endtime (UTCDateTime): end time

        Returns:
            (bool):
        """

    def zero(self, inp, plot=False):
        """
        Zero out data in trace  or stream when it is too close to an eq

        Returns an error if the EQ catalog does not bound the Trace time
        Returns a warning if the EQ catalog does not start far enough before
        the trace to cover an M9 EQ.

        Arguments:
            inp (obspy.core.Trace): seismological data
            plot: plot traces with EQs cut out

        Returns:
            trace with EQ-affected sections set to zero
        """
```
