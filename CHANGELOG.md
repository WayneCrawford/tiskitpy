## v0.1:

combined time series codes that I had in different projects:

- `crawtools`
    - `decimate` 
    - `spectral`

- `rptransient`: All except:
    - `decimate` (used one from `crawtools`)
    - `read_mseed` (changed to independent function

- `wayne_obstools`
    - `fir_corr` => `timeseries`?
    - `Peterson_model`, `plot_PPSD`, `plot_response`, `plot_sensitivity`

## 0.2:

Made sure that all of the equations correspond to equations in Bendat&Piersol
(2010).  Verified that we now get as good of results as the old Matlab code.
Set up testing for the different sub-modules.  Probably changed some
method call parameters.  Created a readthedocs page.

## 0.3:

### MAJOR

- `SpectralDensity`:
    - Added required argument `window_s` to creator
- Renamed `TransferFunction` properties:
    - `values` to `frf` (frequency response function)
    - `uncert` to `uncertainty`
    - `channels` to `channel_names`
- `TimeSpans`:
    - Now created using `spans` by default
- `DataCleaner`
    - always renames cleaned channels, putting information in the
      location code slot (second from last element in string separated by '.'s)
    - Removed unused `fast_calc` parameter from `clean_sdf()`

### BUGFIXES

- Changed `SpectralDensity.from_stream()` z_threshold to apply to log10(spectra)
  rather than (spectra) (otherwise, "standard" zthreshold=3 rejects more than
  half of typical data)
- `TransferFunction()`
    - now provides a "frequency response function" (`frf`)
      and a `corrector` (for data cleaning): the two are the same only if
      all noise is on the output channel (see Bendat and Piersol chapter 6)
    - fixed (and simplified) calculation of uncertainty
- `DataCleaner` uses `TransferFunction`'s `corrector` property.
  
### MINOR

- `SpectralDensity`:
    - Added arguments `ts_starttime` and `ts_endtime` to creator
    - Added ``used()`, `unused()` methods and `window_seconds` property
    - `autospect()` and `.crossspect()` now accept wildcards in channel names
    - Added `channel_name()` method
    - Added `**fig_kw` to `plot_*()` methods
- `TimeSpans`:
    - Added `spans` property
    - Added `invert()` method
    -  `plot()` can include a trace or stream
- Updated documentation
- `DataCleaner`:
    - Created `CleanerString` class to help with channel names

## 0.3.1
    Fix `TransferFunction` errorbar plotting bug
## 0.3.2
    Change `MANIFEST.in` to recursively include subfiles/directories of decimate/
## 0.3.3
    - Decimator.decimate() now returns same data.dtype by default
    - TransferFunction now accepts wildcards for out_chan names
    - Added 'outfile' to TransferFunction.plot()
    - Added TransferFunction.put_response()
    
## 0.4

- Renamed TransferFunction to ResponseFunction, changed internal and method
  names to correspond:

    - .response -> instrument_response
    - .frf -> .value
    - .frf_wrt_counts -> .value_wrt_counts
    - SpectralDensity.channel_response -> SpectralDensity.channel_instrument_response
    - SpectralDensity.put_channel_response -> SpectralDensity.put_channel_instrument_response
    - SpectralDensity.__init__(response=) -> SpectralDensity.__init__(instrument_response)
- Internally renamed `tiskit` to `tiskitpy` 


## 0.4.1
    Rewrote tracking of cleaning steps, a lot of internal work, including
    new classes, but should be invisible when using the command-line codes
    and mostly invisible when using the API.

## 0.5

- Added channel identification by ``tiskitpy_id``, which includes cleaning
  information.
  The ``tiskitpy_id`` for uncleaned data is the ``seed_id``.
- Added class CleanedStream and revised the guts of several classes, including
  renaming ``SpectralDensity.channel_names`` to ``SpectralDensity.ids``.
- CleanRotator class now has a property `variance_reduction` which gives the
  variance reduction obtained during __init__()
- Added `SpectralDensity.plots()` and `.plots_coherences()`, to compare
  multiple `SpectralDensity` objects
  
## 0.5.1

- Fixed bug creating matrix of subplots for autospect