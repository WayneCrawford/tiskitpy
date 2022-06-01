# TiSKit

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
- `rptransient`: calculate and remove periodic transient (VERY manual!).  
 	Based on Matlab code by E Wielandt, used in Deen et al., 2017

- `PetersonNoiseModel``: return the Peterson High and Low Noise Models

## `seismo_tools` submodule: seismology-specific functions
- plot_response`: plot instrument response (command line?)
- plot_sensitivity: plot instrument sensitivity (command line?)

## Installation

Clone or download this repository, then from within the main repository directory, run:

```pip install .```

You can also install in editable mode (for developers), with:

```pip install -e .```

## Creation

At first, combined different time series codes that I had in different projects:

- `crawtools`
    - `decimate` 
    - `spectral`

- `rptransient`: All except:
    - `decimate` (used one from `crawtools`)
    - `read_mseed` (changed to independent function

- `wayne_obstools`
    - `fir_corr` => `timeseries`?
    - `Peterson_model`, `plot_PPSD`, `plot_response`, `plot_sensitivity`

