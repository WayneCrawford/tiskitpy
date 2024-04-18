.. _Decimator:

Decimator
=======================

Decimate data and put decimation information in instrument responses.

Detailed information is in :ref:`tiskitpy.Decimator`

The main methods are:

Constructor
---------------------

- ``Decimate(decimates)``: Set up a series of decimations corresponding
  to the list ``decimates``.

Properties
---------------------
- ``decimates``: the list of decimation factors
- ``decimation_factor``: total decimation (product of ``decimates``)
- ``verbose``: True if object is chatty.


Methods
---------------------

- ``decimate(Stream or Trace)``: Decimate the data
- ``get_band_code(in_band_code, sample_rate)``: return the band code for a given
  sample rate.
- ``update_inventory(inv, ...)``: Return inventory with decimated channels added
- ``update_inventory_from_nslc(inv ...)``: Return inventory with only the
  specified network, station, channel, location(s) updated
 
 Command-line programs
---------------------

Use these programs' `-h` option for help

- ``tiskitpy_decimate_SDS``: Decimates the data in an SDS directory and updates
    the corresponding inventory file.

Example
---------------------

:ref:`tiskitpy.Decimator_example`
