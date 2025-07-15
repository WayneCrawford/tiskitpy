.. _SpectralDensity:

SpectralDensity
=======================

Calculates auto- and cross-spectral densities for a data stream.
Also outputs coherences as well.  Plots any of the above.   

Constructor
---------------------

- :meth:`SpectralDensity.from_stream`: Make a :ref:`SpectralDensity` object from
  a an obspy data :class:`Stream <obspy.core.stream.Stream>`

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

- :meth:`autospect <tiskitpy.SpectralDensity.autospect>`: Auto-spectral density function for a channel
- :meth:`coherence <tiskitpy.SpectralDensity.coherence>`: Coherence between two channels
- :meth:`crossspect <tiskitpy.SpectralDensity.crossspect>`: Cross-spectral density function between two channels
- :meth:`channel_id <tiskitpy.SpectralDensity.channel_id>`: Channel id, expanding
  wildcards and verifying that the result is unique
- :meth:`seed_id <tiskitpy.SpectralDensity.seed_id>`: The specified channel name, expanding
  wildcards and verifying that the result is unique
- :meth:`channel_instrument_response <tiskitpy.SpectralDensity.channel_instrument_response>`: Channel's instrument response
- :meth:`channel_units <tiskitpy.SpectralDensity.channel_units>`: A channel's input (physical) units
- :meth:`units <tiskitpy.SpectralDensity.units>`: Units of a cross- or auto-spectra
- :meth:`coh_signif <tiskitpy.SpectralDensity.coh_signif>`: The coherence significance level

Other Methods
^^^^^^^^^^^^^^^^^^

- :meth:`plot_autospectra <tiskitpy.SpectralDensity.plot_autospectra>`: plot autospectra
- :meth:`plot <tiskitpy.SpectralDensity.plot>`: shortcut for ``plot_autospectra()``
- :meth:`plot_cross_spectra <tiskitpy.SpectralDensity.plot_cross_spectra>`: plot cross- (and auto-) spectra
- :meth:`plot_coherences <tiskitpy.SpectralDensity.plot_coherences>`: plot coherences
- :meth:`plot_one_autospectra <tiskitpy.SpectralDensity.plot_one_autospectra>`: plot autospectra for one channel
- :meth:`plot_one_spectra <tiskitpy.SpectralDensity.plot_one_spectra>`: plot cross-spectra
  for the given channels
- :meth:`plot_one_coherence <tiskitpy.SpectralDensity.plot_one_coherence>`: plot coherence
  for the given channels
- :meth:`plots <tiskitpy.SpectralDensity.plots>`: overlay plot spectra specified in the list

Set Methods
^^^^^^^^^^^^^^^^^^

You probably won't ever use these (should I put a `_` before?)

- :meth:`put_crossspect <tiskitpy.SpectralDensity.put_crossspect>`: put a cross-spectral
  density in the given slot
- :meth:`put_autospect <tiskitpy.SpectralDensity.put_autospect>`: same as
  :meth:`put_crossspect <tiskitpy.SpectralDensity.put_crossspect>`
- :meth:`put_channel_instrument_response <tiskitpy.SpectralDensity.put_channel_instrument_response>`: put a channel
  response in the given slot
- :meth:`replace_channel_id <tiskitpy.SpectralDensity.replace_channel_id>`: change a channel id

Example
---------------------

:ref:`tiskitpy.SpectralDensity_example`
