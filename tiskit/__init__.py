"""
Routines for time series data processing
Uses the obspy seismological Trace, Stream (data) and Inventory (metadata)
classes, but should work for non-seismology datasets as well

## Classes
- `CleanRotator`: rotate data to minimize noise on vertical channel
- `DataCleaner`: Transfer_Function-based data cleaning
- `Decimator`: Decimate time series and update metadata with the decimator's
               response
- `SpectralDensity`: Calculate and manipulate spectral density functions.
- `TimeSpans`: Specify time spans to be removed, kept, zeroed, etc.
- `TransferFunctions`: Transfer functions for a given input channel.
               

## Functions
- `FIR_corr`: transform zero-phase data to minimum phase (only works for
              LCHEAPO loggers, need to update to calculate/work for any
              zero-phase filter)
- `readMSEED`: read in MSEED data, including if the file is too big (> 2 GB)
               for obspy's read() function
- `rptransient`: calculate and remove periodic transient (VERY manual!)
- `PetersonNoiseModel``: return the Peterson High and Low Noise Models

"""
from .data_cleaner import DataCleaner
from .decimate import Decimator
from .read_mseed import read_MSEED
from .clean_rotator import CleanRotator
from .rptransient import PeriodicTransient
# from .seismo_tools import plot_response, plot_sensitivity
from .spectral_density import SpectralDensity, Peterson_noise_model
from .time_spans import TimeSpans
from .transfer_functions import TransferFunctions
## `seismo_tools` submodule: seismology-specific functions
# - plot_response`: plot instrument response (command line?)
# - plot_sensitivity: plot instrument sensitivity (command line?)
