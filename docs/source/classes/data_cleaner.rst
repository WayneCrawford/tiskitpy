DataCleaner
=======================

Remove coherent noise seen on another channel

Detailed information is in :ref:`tiskitpy.DataCleaner`

The main methods are:

Constructor
---------------------

- ``DataCleaner(stream, remove_list,...)``: Calculate the DataCleaner object from
  a data stram and a list of channels to remove

Methods
---------------------

Cleaning
^^^^^^^^^^^^

- ``clean_sdf(sdf)``: Clean an existing spectral density function (approximation).
- ``clean_stream(stream, ...)``: Clean a data stream
- ``clean_stream_to_sdf(stream, ...)``: Calculate SpectralDensity function directly
  from the input stream

Other
^^^^^^^^^^^^

- ``plot()``: plot the transfer functions in the DataCleaner

Example
---------------------

:ref:`tiskitpy.DataCleaner_example`
