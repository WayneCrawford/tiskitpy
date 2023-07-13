.. _ResponseFunctions:

ResponseFunctions
=======================

Calculate frequency response functions for a given input channel and
a range of output channels.

Detailed information is in :ref:`tiskitpy.ResponseFunctions`

The main methods are:

Constructor
---------------------

:`ResponseFunction(SpectralDensity, in_chan, ...)`: 

Properties
---------------------

:`freqs`: Frequency response function frequencies
:`n_windows`: Number of time series data windows used
:`input_channel`: Frequency response function input channel
:`input_units`: Frequency response function input channel units
:`output_channels`: List of the output channels
:`noise_channels`: Noise channels for each xf

Methods
---------------------

Get
^^^^^^^^^^^^^^^^^^^^^

:`coh_signif(prob)`: Coherence significance level with the given probability
:`corrector(output_channel)`: input channel's correction factor with respect
    to the given channel
:`corrector_wrt_counts(output_channel)`: as above, but with respect to raw
    data counts
:`value (output_channel)`: frequency response function
:`value_wrt_counts(output_channel)`: frequency response function with respect
    to raw data counts
:`noise_channel(output_channel)`: Return the channel ("input", "output" or "equal")
    assumed to have incoherent noise
:`output_units(output_channel)`: output channel units
:`instrument_response(output_channel)`: frequency response function's instrument respose 
    (output_channel_instrument_response / input_channel_instrument_response)
:`uncertainty(output_channel)`: uncertainty of the given `value` s
:`uncertainty_wrt_counts(output_channel)`: uncertainty with respect
    to raw data counts

Other
^^^^^^^^^^^^^^^^^^^^^

:`to_norm_compliance(water_depth)`: convert m/s^2 / Pa transfer functions to
    normalized compliance
:`plot`: plot the frequency response functions
:`plot_one(in_chan, out_chan, ...)`: plot one frequency response function

Example
---------------------


see :ref:`tiskitpy.ResponseFunctions_example`
