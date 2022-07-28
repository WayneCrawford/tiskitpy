TransferFunctions class
=======================

Calculate transfer functions (actually, frequency response functions, according
to [BP2010]_) for a given input channel and a range of output channels.

Detailed information is in :ref: `tiskit.TransferFunctions`

The main methods are:

Constructor
---------------------

:`TransferFunction(SpectralDensity, in_chan, ...)`: Calculate transfer functions
    based on the Spectral Density Matrix and for a given input (source)
    channel

Properties
---------------------

:`freqs`: Transfer function frequencies
:`n_windows`: Number of time series data windows used
:`input_channel`: Transfer function input channel
:`input_units`: Transfer function input channel units
:`output_channels`: List of the output channel
:`noise_channels`: List of which channel is assumed to have incoherent noise for each xf

Get Methods
---------------------

:`coh_signif(probability=0.95)`: Coherence significance level with the given
    probability
:`noise_channel(output_channel)`: Which channel ("input", "output" or "equal")
    is assumed to have incoherent noise for the given xf
:`output_units(output_channel)`: output channel units
:`values(output_channel)`: frequency response function
:`uncert(output_channel)`: uncertainty of the given `value` s
:`values_wrt_counts(output_channel)`: frequency response function with respect
    to raw data counts
:`uncert_wrt_counts(output_channel)`: uncertainty with respect
    to raw data counts
:`response(output_channel)`: transfer function's instrument respose 
    (output_channel_response / input_channel_response)



Other
---------------------

:`plot`: plot the transfer functions
:`plot_one(in_chan, out_chan, ...)`: plot one transfer function

Example
---------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import TransferFunction
  
