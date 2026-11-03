## 2.4.0
- Minor
    - ``FIRConverter`` now handles consecutive data in separate files
        - ``apply()`` takes preceding and following data (``prev_stream``,
          ``next_stream``) and interpolates small gaps or accepts edge
          effects (``gap_handling``)
        - Neighbours offset from the sample grid by up to ``max_misalignment``
          samples (default 0.5, e.g. clock-corrected day files) are treated
          as consecutive
        - ``pad_lengths()`` returns the number of neighbouring samples needed
        - ``output_stage_values()`` and ``convert_response_stage()`` give the
          FIR stage for the converted channels' metadata (minimum phase
          coefficients, no symmetry, zero delay and correction)
        - ``apply_SDS()`` fixed and tested: uses neighbouring day files
          (including across years), writes output files with the new
          location code and skips non-waveform files (e.g. ASCII logs)
        - ``sps`` is now optional in ``from_zeros_file()`` (it is only
          informational), which fixes ``tiskitpy FIR_convert_SDS --zeros_file``
        - ``update_inventory()`` and ``convert_channel()`` add the converted
          channels to an inventory
    - New ``tiskitpy FIR_convert_SDS`` subcommand, with the same interface as
      ``tiskitpy decimate_SDS``: takes both an SDS directory and a StationXML
      file, converts only the channels whose response contains the converter
      FIR, and writes an SDS and StationXML file with the new channels
    - ``Decimator.decimate()`` takes neighbouring data (``prev_stream``,
      ``next_stream``) and a decimation ``phase``, so that data decimated in
      separate pieces are identical to data decimated in one piece
    - ``decimate_SDS`` output is continuous across day files (including
      across years)
- Maintenance
    - Fixed a bug in the comment text in compliance output CSV files
    - ``decimate_SDS``: fixed ``tiskitpy decimate_SDS`` subcommand (crashed),
      re-runs into an existing output directory, and decimated channels
      that only change location code being rejected from the inventory
    - Shared SDS and inventory tools in ``util_functions/sds_tools.py`` and
      ``util_functions/inventory_tools.py``

### 2.3.1

No change, just to avoid pip install conflict

## 2.3.0

- Minor
    - Added ``Pressure_Noise_Model()`` function
    - Added ``FIRConverter`` class: converts zero-phase filtered output to its
      equivalent minimum phase form
    - Improvements to ``PeriodicTransient()``, including:
    	    - Automatic filtering of data for the modeling and matching part
    	    - Default calculation of EQ times (using ``TimeSpans.from_eqs``)
    	    - Using filtered data to improve alignment of diracs
    	- `TimeSpans.from_eqs()` now saves QML files in `.time_spans/` directory
- Maintenance
    - Added documentation on ``PeriodicTransient`` class
    - ``Decimator`` verifies that decimates contains only valid integers
    - Added (not uniformly) some plotting parameters:
	    - ``title`` to ``SpectralDensity.plot_coherences()`` and
	      ``ResponseFunctions.plot()``
	    - ``grid`` to `ResponseFunctions.plot()` and ``SpectralDensity.plot()``
    - Normalized order and names of the first three arguments to the
     ``SpectralDensity`` ``.plots()`` and ``.plots_coherences()`` methods 
    - Corrected legend line colors in ``plots_coherences()``

## 2.2.0

- Minor
    - ``ResponseFunctions`` class:
        - Added ``coherence`` property
    - ``Compliance`` class:
        - Added ``coherence`` property 
        - Added  ``coherence`` column to CSV file written by ``write()`` method
- Maintenance
    - ``TimeSpans`` class:
        - Added ``split()``, ``__and__()`` and ``__or__()`` methods
        - Grandfathered ``used_times()`` method, now named ``used_spans()``
    - ``SpectralDensity`` class:
        - Modified ``from_stream()`` to handle streams with data gaps
          (masked arrays or separate Traces with the same seed_id)
    - ``CleanRotator`` and ``DataCleaner`` classes:
        - ``apply()`` methods now return, by default, the same data dtype as
          the input

**post1 and post2**
Debugging of script paths, and debug documentation

### 2.1.2

- Fix error in SpectralDensity class when specifying time_spans (namespace conflict)

### 2.1.1

- Allow ``compliance.write_out()`` to accept ``None`` for p_response and z_response,
  in which case the compliance is assumed to already have been in counts/counts
- SeafloorSynthetic.streams() now puts `-self.water_depth` into the output inventory's
  station and channel elevations
  
## 2.1
- General
	- Closes figures that are not shown
	- Added ``plot_compliance_stack()`` function
- ``CleanRotator`` and ``SeisRotate`` classes:
	- Added optional ``H_over_Z`` parameter to account for inequal Horizontal and Z sensitivites
- ``SpectralDensity`` class:
	- Changed default windowtype from ``prol1pi`` to ``prol4pi``
- ``ResponseFunctions`` class:
	- removed ``to_norm_compliance()`` method (replaced by ``Compliance.from_response_functions()``
	- Added response function units to plots
- ``SeafloorSynthetic`` class:
	- Change default tilt noise parameters to more typical/reasonable values
- ``Compliance`` class:
	- Fixed Path bug
	- Added ``from_seafloor_synthetic`` class method
	- Added ``max_freq()`` static method and use it to limit frequencies by default
	- Added '_compliance' to ``write()`` method output filename
	- Added ``write_counts()`` method
	- Removed ``frf_to_compliance()`` (replaced by ``Compliance.from_from_response_functions()``)

# 2.0

Renamed ``ComplianceNoise`` to ``SeafloorSynthetic``

--------

## 1.1
- Fixed bug in ``ResponseFunctions`` plotting
- Added ``Compliance`` class
- ``ComplianceNoise`` now gives expected outputs
- ``SpectralDensity``: ``plots()`` and ``plot_coherences()`` methods now accept objects without identical seed_ids
- ``EarthModel1D``: Added ``plot()`` and ``calc_ncomp()`` methods

# 1.0

#### Added
-  ``ComplianceNoise`` class and documentation

#### Fixed
    - ``SpectralDensity`` plot overwriting axis labels
    - ``SpectralDensity`` ``remove_eqs=True`` now works

#### Renamed
    - ``DataCleaner.clean_stream()`` to ``DataCleaner.apply()``
    - ``DataCleaner.clean_sdf()`` to ``DataCleaner.apply_to_sdf()``

#### Other major changes
    - TimeSpans.from_eqs() now takes ``time_bounds`` input instead of ``starttime`` and ``endtime``
    - ``SpectralDensity``, ``DataCleaner`` and ``CleanRotate`` classes take ``remove_eqs``
      and ``avoid_spans`` as input and provide an ``avoided_spans`` property.

--------

## 0.1:

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

#### MAJOR

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

#### BUGFIXES

- Changed `SpectralDensity.from_stream()` z_threshold to apply to log10(spectra)
  rather than (spectra) (otherwise, "standard" zthreshold=3 rejects more than
  half of typical data)
- `TransferFunction()`
    - now provides a "frequency response function" (`frf`)
      and a `corrector` (for data cleaning): the two are the same only if
      all noise is on the output channel (see Bendat and Piersol chapter 6)
    - fixed (and simplified) calculation of uncertainty
- `DataCleaner` uses `TransferFunction`'s `corrector` property.
  
#### MINOR

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

### 0.3.1

Fix `TransferFunction` errorbar plotting bug
    
### 0.3.2
Change `MANIFEST.in` to recursively include subfiles/directories of `decimate/`

### 0.3.3

- `Decimator.decimate()` now returns same data.dtype by default
- `TransferFunction` now accepts wildcards for out_chan names
- Added 'outfile' parameter to `TransferFunction.plot()`
- Added `TransferFunction.put_response()`
    
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
  
## 0.5.2

- Fixed bug creating matrix of subplots for autospect
- Fixed bug overwriting axis labels for most plots
- Added ``remove_eqs`` option to SpectralDensity.from_stream()
- Plots of ResponseFunctions only show selected frequencies

## 0.5.3

- Fixed overwriting axis labels for plots_coherences()

