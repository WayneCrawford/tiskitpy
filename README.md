# TiSKit

Routines for time series data processing

Uses the obspy seismological Trace, Stream (data) and Inventory (metadata)
classes, but should work for non-seismology datasets as well


[Documentation](https://tiskit.readthedocs.io/en/latest/index.html)


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

- `PetersonNoiseModel`: return the Peterson High and Low Noise Models


## `seismo_tools` submodule: seismology-specific functions

- `plot_response`: plot instrument response (command line?)
- `plot_sensitivity`: plot instrument sensitivity (command line?)


## Installation

First, install `obspy` using the instructions on their webpage.
Then, in the pip/conda environment that contains obspy...

### From this repository

Clone or download this repository, then from within the main repository directory, run:

`pip install .`

You can also install in editable mode (for developers), with:

`pip install -e .`

### Using `pip`

Type `pip install tiskit-py`

Note that I had to call this module `tiskit-py` on PyPI, rather than `tiskit`,
so you may have to change the `import` calls in your code from `import tiskit`
to `import tiskit-py as tiskit`
