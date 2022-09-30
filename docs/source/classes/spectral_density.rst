.. _SpectralDensity:

SpectralDensity
=======================

Calculates auto- and cross-spectral densities for a data stream.
Also outputs coherences as well.  Plots any of the above.   

Detailed information is in :ref:`tiskit.SpectralDensity`

The main methods are:

Constructor
---------------------

:`SpectralDensity.from_stream(...)`: Make a `SpectralDensity` object from
    an `obspy` data stream

Properties
---------------------

:`channel_names`: a list of the channel names
:`freqs`: the frequencies of the spectral density functions
:`n_windows`: the number of data windows used to calculate the spectra
:`window_type`: The type of tapering window used when calculating the
    spectral densities
:`window_seconds`: Length of each window, in seconds
:`starttimes`: get a list containing the starttimes for each data window
:`used_times`: time spans used to calculate spectra
:`unused_times`: time spans rejected or otherwise unused to calculate spectra


Methods
---------------------

Get Methods
^^^^^^^^^^^^^^^^^^

:`autospect(channel)`: the channel's auto-spectral density function
:`coherence(in_channel. out_channel)`: the coherence between the given
    channels
:`crossspect(in_channel. out_channel)`: the cross-spectral density function
    between the given channels
:`channel_name(channel_name)`: get the specified channel name, expanding
    wildcards and verifying that the result is unique
:`channel_response(channel)`: the channel's instrument response
:`channel_units(channel)`: the channel's input (physical) units
:`units(in_channel, out_channel)`: get the units of the corresponding
    cross- or auto-spectra
:`coh_signif(probability)`: get the coherence significance level

Other Methods
^^^^^^^^^^^^^^^^^^

:`plot_autospectra(...)`: plot autospectra
:`plot(...)`: shortcut for `plot_autospectra()`
:`plot_cross_spectra(...)`: plot cross- (and auto-) spectra
:`plot_coherences(...)`: plot coherences
:`plot_one_autospectra(channel, ...)`: plot autospectra for one channel
:`plot_one_spectra(in_channel, out_channel, ...)`: plot cross-spectra for the given channels
:`plot_one_coherence(in_channel, out_channel, ...)`: plot coherence for the given channels

Set Methods
^^^^^^^^^^^^^^^^^^

You probably won't ever use these (should I put a `_` before?)

:`put_crossspect(in_channel, out_channel, spect)`: put a cross-spectral density in the given slot
:`put_autospect(channel, spect)`: same as `put_crossspect(channel, channel, spect)`
:`put_channel_response(channel, response)`: put a channel response in the given slot
:`replace_channel_name(channel, replacement)`: change a channel name

Example
---------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import SpectralDensity
  
  stream = read('mydata.mseed', 'MSEED')
  inv = read_inventory('myinv', 'STATIONXML')
  sd = SpectralDensity.from_stream(stream, inv=inv)
  sd.plot()
  print(sd)
