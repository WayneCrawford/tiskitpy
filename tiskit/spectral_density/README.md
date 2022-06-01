# Spectral programs

## Active
- `spectral_density.py`: Defines `SpectraDensity` class, which calculates (or
  imports) cross- and auto-spectra.  Also calculates coherences  as derived
  properties, and coherence significance levels.
    - Currently has routines to calculate and apply cleaning transfer functions,
      though these might more logically sit in `transfer_function.py`
- `transfer_function.py`: Defines `TransferFunction` class, which calculates
  transfer functions from `SpectralDensity` objects
- `xf_cleaner.py`: Defines `XFCleaner` class, which calculates and applies
  cleaning transfer function lists
- `Peterson_noise_model.py`: defines function outputting Peterson high and low
  noise models.
  
## Grandfathered (?)
- `PSD.py`: calculate power spectral density for a single channel. Older than
  `SpectralDensity` and not linked to it.  It corrects for the instrument
  response during calculation, if an inventory was provided to the `calc()`
  method.
- `coherence.py`: Defines `Coherence` class, which is pretty much a wrapper for
  matplotlib.cohere or scipy.signal.coherence().  I don't know if there is
  any use not that we have `SpectralDensity`
- `utils.py`: contains some functions used by `PSD.py`?

## Remaining questions
- Is removing the transfer function times the spectrum really accurate?  It 
  seems that it would be more honest to remove the transfer function times
  each FFT.  In TiSKit, I started with the raw data and recalculated the
  spectra for every newly added transfer function (so applying to each FFT)