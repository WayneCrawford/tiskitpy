DataCleaner
=======================

Remove coherent noise seen on another channel

Constructor
---------------------

- :class:`DataCleaner <tiskitpy.DataCleaner>`: Calculate the DataCleaner object from
  a data stream and a list of channels to remove

Methods
---------------------

Cleaning
^^^^^^^^^^^^

- :meth:`apply_to_sdf <tiskitpy.DataCleaner.apply_to_sdf>`: Clean an existing spectral density function (approximation).
- :meth:`apply <tiskitpy.DataCleaner.apply>`: Clean a data stream
- :meth:`apply_to_streams_sdf <tiskitpy.DataCleaner.apply_to_streams_sdf>`: Calculate SpectralDensity function directly
  from the input stream

Other
^^^^^^^^^^^^

- :meth:`plot <tiskitpy.DataCleaner.plot>`: plot the transfer functions in the DataCleaner

Example
---------------------

:ref:`tiskitpy.DataCleaner_example`
