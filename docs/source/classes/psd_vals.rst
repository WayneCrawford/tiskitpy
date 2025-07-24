.. _PSDVals:

PSDVals
=======================

Synthetic Power Spectral Densities and their waveforms

A helper class for :class:`tiskitpy.SeafloorSynthetic`


Constructor
---------------------

- :class:`tiskitpy.compliance.PSDVals`:
- :meth:`tiskitpy.compliance.PSDVals.from_loglog_slope`: Make an object with a constant slope in loglog space

Properties
---------------------

- ``freqs`` (list): monotonically increasing frequencies
- ``values`` (list): PSD values corresponding to each frequency
- ``value_units`` (str): The units of the PSD values

Dependent Properties
---------------------
- accel_as_vel: returns a PSDVals object with values divided by 2*pi*f

Methods
---------------------

- :meth:`copy <tiskitpy.compliance.PSDVals.copy>`: Deep copy of the object
- :meth:`resample <tiskitpy.compliance.PSDVals.resample>`: Resample the data using interpolation, modifying the object
- :meth:`resample_values <tiskitpy.compliance.PSDVals.resample_values>`: Resample the data using interpolation, returning just the resampled values
- :meth:`plot <tiskitpy.compliance.PSDVals.plot>`: Plot the PSD
- :meth:`as_trace <tiskitpy.compliance.PSDVals.as_trace>`: Return a :class:`obspy.Trace` with the given spectral shape

Static Methods
---------------------
- :meth:`sloped_freqs_and_values <tiskitpy.compliance.PSDVals.sloped_freqs_and_values>`: Create the freqs_and_values input to the constructor, for a loglog slope
- 

Example
---------------------

None for now
