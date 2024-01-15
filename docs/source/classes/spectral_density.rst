.. _SpectralDensity:

SpectralDensity
=======================

Calculates auto- and cross-spectral densities for a data stream.
Also outputs coherences as well.  Plots any of the above.   

Detailed information is in :ref:`tiskitpy.SpectralDensity`

The main methods are:

Constructor
---------------------

- ``SpectralDensity.from_stream(...)``: Make a :ref:`SpectralDensity` object from
  an ``obspy`` data stream

Properties
---------------------

- ``ids``: a list of the channel ids
- ``seed_ids``: a list of the channel seed_ids
- ``freqs``: the frequencies of the spectral density functions
- ``n_windows``: the number of data windows used to calculate the spectra
- ``window_type``: The type of tapering window used when calculating the
  spectral densities
- ``window_seconds``: Length of each window, in seconds
- ``starttimes``: get a list containing the starttimes for each data window
- ``used_times``: time spans used to calculate spectra
- ``unused_times``: time spans rejected or otherwise unused to calculate spectra


Methods
---------------------

Get Methods
^^^^^^^^^^^^^^^^^^

- ``autospect(id)``: the channel's auto-spectral density function
- ``coherence(in_id, out_id)``: the coherence between the given
  channels
- ``crossspect(in_id, out_id)``: the cross-spectral density function
  between the given channels
- ``channel_id(id)``: get the specified channel name, expanding
  wildcards and verifying that the result is unique
- ``seed_id(id)``: get the specified channel name, expanding
  wildcards and verifying that the result is unique
- ``channel_instrument_response(id)``: the channel's instrument response
- ``channel_units(id)``: the channel's input (physical) units
- ``units(in_id, out_id)``: get the units of the corresponding
  cross- or auto-spectra
- ``coh_signif(probability)``: get the coherence significance level

Other Methods
^^^^^^^^^^^^^^^^^^

- ``plot_autospectra(...)``: plot autospectra
- ``plot(...)``: shortcut for ``plot_autospectra()``
- ``plot_cross_spectra(...)``: plot cross- (and auto-) spectra
- ``plot_coherences(...)``: plot coherences
- ``plot_one_autospectra(channel, ...)``: plot autospectra for one channel
- ``plot_one_spectra(in_id, out_id, ...)``: plot cross-spectra
  for the given channels
- ``plot_one_coherence(in_id, out_id, ...)``: plot coherence
  for the given channels

Set Methods
^^^^^^^^^^^^^^^^^^

You probably won't ever use these (should I put a `_` before?)

- ``put_crossspect(in_id, out_id, spect)``: put a cross-spectral
  density in the given slot
- ``put_autospect(id, spect`)`: same as
  ``put_crossspect(id, channel, spect)``
- ``put_channel_instrument_response(id, response)``: put a channel
  response in the given slot
- ``replace_channel_id(id, replacement_id)``: change a channel id

Example
---------------------

:ref:`tiskitpy.SpectralDensity_example`
