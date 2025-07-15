.. _ResponseFunctions:

ResponseFunctions
=======================

Calculate frequency response functions for a given input channel and
a range of output channels.

Constructor
---------------------

- :class:`ResponseFunctions <tiskitpy.ResponseFunctions>` 

Properties
---------------------

- ``freqs``: Frequency response function frequencies
- ``input_channel_id``: Frequency response function input channel
- ``output_channel_ids``: List of the output channels
- ``input_units``: Frequency response function input channel units
- ``noise_channels``: Noise channels for each xf
- ``input_clean_sequence``: Clean sequence text
- ``n_windows``: Number of time series data windows used

Methods
---------------------

Get
^^^^^^^^^^^^^^^^^^^^^

- :meth:`coh_signif <tiskitpy.ResponseFunctions.coh_signif>`: Coherence significance level with the given probability
- :meth:`corrector <tiskitpy.ResponseFunctions.corrector>`: Input channel's correction factor with
  respect to the given channel
- :meth:`corrector_wrt_counts <tiskitpy.ResponseFunctions.corrector_wrt_counts>`: As above, but with respect to
  raw data counts
- :meth:`value <tiskitpy.ResponseFunctions.value>`: Frequency response function
- :meth:`value_wrt_counts <tiskitpy.ResponseFunctions.value_wrt_counts>`: Frequency response function with
  respect to raw data counts
- :meth:`noise_channel <tiskitpy.ResponseFunctions.noise_channel>`: Return the channel ("input", "output" or
  "equal") assumed to have incoherent noise
- :meth:`output_units <tiskitpy.ResponseFunctions.output_units>`: Output channel units
- :meth:`instrument_response <tiskitpy.ResponseFunctions.instrument_response>`: Frequency response function's
  instrument response (use to switch from counts to physical units)
- :meth:`uncertainty <tiskitpy.ResponseFunctions.uncertainty>`: Uncertainty of the given ``value`` s
- :meth:`uncertainty_wrt_counts <tiskitpy.ResponseFunctions.uncertainty_wrt_counts>`:
  Uncertainty with respect to raw data counts

Other
^^^^^^^^^^^^^^^^^^^^^

- :meth:`to_norm_compliance <tiskitpy.ResponseFunctions.to_norm_compliance>`: Convert m/s^2 / Pa transfer functions to
  normalized compliance
- :meth:`plot <tiskitpy.ResponseFunctions.plot>`: Plot the frequency response functions
- :meth:`plot_one <tiskitpy.ResponseFunctions.plot_one>``: Plot one frequency response function

Example
---------------------

:ref:`tiskitpy.ResponseFunctions_example`
