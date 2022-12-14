"""
Routines for time series data processing

Use the obspy seismological Trace, Stream (data) and Inventory (metadata)
classes, but should work for non-seismology datasets as well if you stuff
them into those classes

Classes
=============

- ``CleanRotator`` : rotate data to minimize noise on vertical channel
- ``DataCleaner`` : Transfer_Function-based data cleaning
- ``Decimator`` : Decimate time series and update metadata with the decimator's response
- ``PeriodicTransient`` : calculate and remove periodic transient (VERY manual!)
- ``SpectralDensity`` : Calculate and manipulate spectral density functions.
- ``TimeSpans`` : Specify time spans to be removed, kept, zeroed, etc.
- ``TransferFunctions`` : Transfer functions for a given input channel.
               

Functions
=============

- ``fir2caus`` : transform zero-phase data to minimum phase (only works for LCHEAPO loggers, need to update to calculate/work for any zero-phase filter)
- ``read_MSEED`` : read in MSEED data, including if the file is too big (> 2 GB) for obspy's read() function
- ``Peterson_noise_model`` : return the Peterson High and Low Noise Models

"""
from .clean_rotator import CleanRotator
from .data_cleaner import DataCleaner, DCTF, DCTFs, CleanerString
from .decimate import Decimator
from .rptransient import PeriodicTransient
from .spectral_density import SpectralDensity, Peterson_noise_model
from .time_spans import TimeSpans
from .transfer_functions import TransferFunctions
from .read_mseed import read_MSEED
from .fir_corr import fir2caus
# from .seismo_tools import plot_response, plot_sensitivity
## `seismo_tools` submodule: seismology-specific functions
# - plot_response`: plot instrument response (command line?)
# - plot_sensitivity: plot instrument sensitivity (command line?)
