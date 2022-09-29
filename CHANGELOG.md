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
- `TransferFunction()` now provides a "frequency response function" (`frf`)
  and a `corrector` (for data cleaning): the two are the same only if
  all noise is on the output channel (see Bendat and Piersol chapter 6)
- `DataCleaner` uses `TransferFunction`'s `corrector` property.
  
### MINOR

- `SpectralDensity`:
    - Added arguments `ts_starttime` and `ts_endtime` to creator
    - Added ``used()`, `unused()` methods and `window_seconds` property
    - `autospect()` and `.crossspect()` now accept wildcards in channel names
    - Added `channel_name()` method
- `TimeSpans`:
    - Added `spans` property
    - Added `invert()` method
    -  `plot()` can include a trace or stream
- Updated documentation
- `DataCleaner`:
    - Created `CleanerString` class to help with channel names
