# TiSKit

Remove periodic transient(s) from seismological data

Based on Matlab code by E Wielandt, used in Deen et al., 2017

## Installation

Clone or download this repository, then from within the main repository directory, run:

```pip install .```

You can also install in editable mode (for developers), with:

```pip install -e .```

## Overview
- `timeseries`: functions that work on time series (probably all in obspy Stream or Trace format)
	- `transient`: identification and removal of transients
	- `decimate`: decimate data and provide updated metadata
- `spectral`: creation and manipulation of spectral data

## Creation

At first, just combines different time series codes that
I had in different directories:

- `crawtools`
    - `decimate` => `timeseries` (or `decimate`?)
    - `spectral` => `spectral`

- `rptransient`:
    - `decimate` `time_spans` and `eq_spans` =>  `timeseries`
    - `dirac_comb`, `transients` and `periodic_transient` => `transient` (or `timeseries`/`transient`?)
    - `rotate_clean` => `??`
    - `utils` => wherever appropriate
    - `read_mseed` => `timeseries`?

- `wayne_obstools`
    - `bbobs_clean` => `transient` and `timeseries` (how much overlap with `rptransient` functions?)
    - `fir_corr` => `timeseries`?
    - `Peterson_model`, `plot_PPSD`, `plot_response`, `plot_sensitivity` => `spectral` (how much overlap with crawtools/spectral?)
    - `instrument_tests` => **`lcheapo`**/`lctest`?

