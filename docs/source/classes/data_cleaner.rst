DataCleaner class
=======================

Remove coherent noise seen on another channel

The main methods are:

Constructor
---------------------

:`DataCleaner(stream, remove_list,...)`: Calculate the DataCleaner object using
    the data stram and a list of channels to remove (in order)

Cleaning Methods
---------------------

:`clean_sdf(sdf)`: Clean an existing spectral density function.  Always an
    approximation, because a SpectralDensity object does not contain enough
    information for a full analysis
:`clean_stream(stream, ...)`: Clean a data stream
:`clean_stream_to_sdf(stream, ...)`: Calculate SpectralDensity function directly
    from the input stream, applying the DataCleaner values to each FFT


Other Methods
---------------------

:`plot()`: plot the transfer functions in the DataCleaner

Example
---------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import DataCleaner
  
