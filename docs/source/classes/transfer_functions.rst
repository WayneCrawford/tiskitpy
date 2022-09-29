TransferFunctions class
=======================

Calculate transfer functions (actually, frequency response functions, according
to [BP2010]_) for a given input channel and a range of output channels.

Detailed information is in :ref:`tiskit.TransferFunctions`

The main methods are:

Constructor
---------------------

:`TransferFunction(SpectralDensity, in_chan, ...)`: 

Properties
---------------------

:`freqs`: Transfer function frequencies
:`n_windows`: Number of time series data windows used
:`input_channel`: Transfer function input channel
:`input_units`: Transfer function input channel units
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
:`frf (output_channel)`: frequency response function
:`frf_wrt_counts(output_channel)`: frequency response function with respect
    to raw data counts
:`noise_channel(output_channel)`: Return the channel ("input", "output" or "equal")
    assumed to have incoherent noise
:`output_units(output_channel)`: output channel units
:`response(output_channel)`: transfer function's instrument respose 
    (output_channel_response / input_channel_response)
:`uncertainty(output_channel)`: uncertainty of the given `value` s
:`uncertainty_wrt_counts(output_channel)`: uncertainty with respect
    to raw data counts

Other
^^^^^^^^^^^^^^^^^^^^^

:`plot`: plot the transfer functions
:`plot_one(in_chan, out_chan, ...)`: plot one transfer function

Example
---------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import TransferFunction
  
