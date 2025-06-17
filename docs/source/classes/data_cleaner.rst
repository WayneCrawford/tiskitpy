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

- ``apply_to_sdf(sdf)``: Clean an existing spectral density function (approximation).
- ``apply(stream, ...)``: Clean a data stream
- ``apply_to_streams_sdf(stream, ...)``: Calculate SpectralDensity function directly
  from the input stream

Other
^^^^^^^^^^^^

- ``plot()``: plot the transfer functions in the DataCleaner

Example
---------------------

:ref:`tiskitpy.DataCleaner_example`
