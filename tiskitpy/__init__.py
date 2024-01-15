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

Command-line programs
=========================

- ``tiskitpy_decimate_SDS`` : Decimate data stored in a SeisComp Data Structure
    database.
    Inserts the data into the same database and creates a new StationXML file
    (based on an existing StationXML file for the input database)
"""
from .clean_rotator import CleanRotator
from .cleaned_stream import CleanedStream
from .logger import init_logger
from .data_cleaner import DataCleaner, RFList
from .decimate import Decimator
from .rptransient import PeriodicTransient
from .spectral_density import SpectralDensity, Peterson_noise_model
from .time_spans import TimeSpans
from .response_functions import ResponseFunctions
from .read_mseed import read_MSEED
from .fir_corr import fir2caus
from .utils import CleanSequence
# from .utils import remove_cleaner_string, CleanerString, CleanSequence
