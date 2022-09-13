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

- Added required argument `window_s` to SpectralDensity() creator
- Renamed `TransferFunction.values` to `TransferFunction.frf` (frequency
  cresponse function)
- Renamed  `TransferFunction.channels` to `TransferFunction.channel_names`
- DataCleaner always renames cleaned channels, putting information in the
  location code slot (second from last element in string separated by '.'s)

### BUGFIXES

- Changed SpectralDensity.from_stream() z_threshold to apply to log10(spectra)
  rather than (spectra) (otherwise, "standard" zthreshold=3 rejects more than
  half of typical data)
- TransferFunction() now provides a "value" (transfer function) and "corrector"
  (for data cleaning) property, because the two are only the same in the case
  where all noise is on the output channel (see Bendat and Piersol chapter 6)
- DataCleaner uses TransferFunction's "corrector" property.
  
### MINOR

- Added arguments `ts_starttime` and `ts_endtime` to SpectralDensity() creator
- Added `TimeSpans.invert()` method
- Allowed `TimeSpans.plot()` to include a trace or stream
- Added `SpectralDensity` `used()`, `unused()` methods and `window_seconds`
  property
- Updated documentation
- `SpectralDensity.autospect()` and `.crossspect()` now accept wildcards
- New `SpectralDensity.channel_name()` method
- Created CleanerString class to help with DataCleaner cleaned channel names
