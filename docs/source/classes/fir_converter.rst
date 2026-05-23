.. _FIRConverter:

FIRConverter
=======================

Convert zero-phase FIR to minimum phase.

Constructors
---------------------

- :class:`FIRConverter <tiskitpy.FIRConverter>`
- :meth:`FIRConverter.from_zeros <tiskitpy.FIRConverter.from_zeros>`:
  Calculate correction from a list of FIR zeros
- :meth:`FIRConverter.from_zeros <tiskitpy.FIRConverter.from_zeros>`:
  Calculate correction from a file containing a list of FIR zeros
- :meth:`FIRConverter.from_corr_file <tiskitpy.FIRConverter.from_zeros>`:
  Read in correction from a JSON file
- :meth:`FIRConverter.from_builtin <tiskitpy.FIRConverter.from_zeros>`:
  Set correction from one of the built-in cases

Properties
---------------------
- ``firzeros``: Linear-phase FIR zeros
- ``equivalent_zeros``: Equivalent minimum-phase FIR zeros
- ``sps``: samples per second output of the decimation phase associated with
  this FIR (0 means that it is used for more than one output sampling rate)
- ``timetag``: Number of samples that the linear-phase FIR must be shifted
  to make it zero-phase (point of symmetry: usually one-half of the FIR length)
- ``converter_type``: ???
- ``decimation_factor``: decimation associated with this FIR
- ``a``: AR coefficients
- ``b``: MA coefficients
- ``z``:` roots of ``firzeros``

Methods
---------------------

- :meth:`apply <tiskitpy.FIRConverter.apply>`: Apply correction to a data stream
- :meth:`apply_SDS <tiskitpy.FIRConverter.apply_SDS>`: Apply correction to an SDS repository
- :meth:`apply <tiskitpy.FIRConverter.plot_zplane>`: Plot z-plane representation of the FIR
- :meth:`apply <tiskitpy.FIRConverter.plot_impulse_parts>`: Plot the impulse response components of the FIR
- :meth:`apply <tiskitpy.FIRConverter.write_json>`: Write the correction to a JSON file
 
Example
---------------------

:ref:`tiskitpy.FIRConverter_example`
