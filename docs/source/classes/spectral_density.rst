SpectralDensity class
=======================

Calculates auto- and cross-spectral densities for a data stream.  Outputs coherences as well.  Plots any of the above.   The main methods are:

Constructor
---------------------

:`SpectralDensity.from_stream(...)`: Make a `SpectralDensity` object from a data stream

Properties
---------------------

:`channels`: a list of the channels
:`freqs`: get an array of the frequencies
:`n_windows`: the number of data windows used to calculate the spectra
:`starttimes`: get a list containing the starttimes for each data window
:`window_type`: The type of tapering window used when calculating the spectral densities


Get Methods
---------------------

:`autospect(channel)`: get the auto-spectral density function for the given channel
:`coherence(in_channel. out_channel)`: get the coherence between the given channels
:`crossspect(in_channel. out_channel)`: get the cross-spectral density between the given channels
:`channel_response(channel)`: get the instrument response for the given channel
:`channel_units(channel)`: get the input (physical) units for the given channel
:`units(in_channel, out_channel)`: get the units of the corresponding cross- or auto-spectra
:`coh_signal(probability)`: return the *probability*% coherence significance level

Set Methods
---------------------

You probably won't ever use these (should I put a `_` before?)

:`put_autospect(channel, spect)`: equivalent to `put_crossspect(channel, channel, spect)`
:`put_crossspect(in_channel, out_channel, spect)`: put a cross-spectral density in the given slot
:`put_channel_response(channel, response)`: put a channel response in the given slot
:`replace_channel_name(channel, replacement)`: change a channel name

Other
---------------------

:`plot_autospectra(...)`: plot autospectra
:`plot(...)`: shortcut for `plot_autospectra()`
:`plot_cross_spectra(...)`: plot cross- (and auto-) spectra
:`plot_coherences(...)`: plot coherences
:`plot_one_autospectra(channel, ...)`: plot autospectra for one channel
:`plot_one_spectra(in_channel, out_channel, ...)`: plot cross-spectra for the given channels
:`plot_one_coherence(in_channel, out_channel, ...)`: plot coherence for the given channels

Example
---------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import SpectralDensity
  
  stream = read('mydata', 'MSEED')
  inv = read_inventory('myinv', 'STATIONXML')
  sd = SpectralDensity.from_stream(stream, inv=inv)
  sd.plot()
  print(sd)