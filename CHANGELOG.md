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
