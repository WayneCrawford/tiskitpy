"""
Routines for time series data processing

Use the obspy seismological Trace, Stream (data) and Inventory (metadata)
classes, but should work for non-seismology datasets as well if you stuff
them into those classes

Classes
=========================

- ``CleanRotator`` : Rotate data to minimize noise on vertical channel
- ``DataCleaner`` : Transfer_Function-based data cleaning
- ``Decimator`` : Decimate time series and update metadata with the
    decimator's response
- ``PeriodicTransient`` : Calculate and remove periodic transient (VERY manual!)
- ``SpectralDensity`` : Calculate and manipulate spectral density functions.
- ``TimeSpans`` : Specify time spans to be removed, kept, zeroed, etc.
- ``ResponseFunctions`` : Frequency response functions for a given input channel.
- ``CleanedStream`` : obspy `Stream` subclass that handles ``cleaned_sequence``
    information
              

Functions
=========================

- ``fir2caus`` : Transform zero-phase data to minimum phase (only works
    for LCHEAPO loggers, need to update to calculate/work for any
    zero-phase filter)
- ``read_MSEED`` : Read MSEED data, even if the file is too big (> 2 GB)
    for obspy's read() function
- ``Peterson_noise_model`` : Return the Peterson High and Low Noise Models
- ``stream_synchronize`` : Return a synchronized stream (all traces have
    same starttime and endtime).  Raises ValueError if not all streams have
    the same sample_rate
- ``stream_unmask`` : unmasks data in a stream, interpolating to fill any gaps

Command-line programs
=========================

Use the `-h` option for help

- ``tiskitpy_decimate_SDS`` : Decimate data stored in a SeisComp Data Structure
  database.
- ``tiskitpy_get_SDS_inventory``: Return the inventory corresponding to a
  SeisComp Data Structure database, using the FDSN Station webservice
"""
from .clean_rotator import CleanRotator
from .cleaned_stream import CleanedStream
from .data_cleaner import DataCleaner, RFList
from .decimate import Decimator
from .rptransient import PeriodicTransient
from .spectral_density import SpectralDensity, Peterson_noise_model
from .time_spans import TimeSpans
from .response_functions import ResponseFunctions
# Functions
from .read_mseed import read_MSEED
from .fir_corr import fir2caus
from .utils import stream_synchronize, stream_unmask
# These are only here for tests, there is probably a better way to access/hide them
from .logger import init_logger
from .utils import CleanSequence

# from .utils import remove_cleaner_string, CleanerString, CleanSequence
