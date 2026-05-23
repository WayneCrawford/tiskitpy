.. _Decimator:

Decimator
=======================

Decimate data and put decimation information in instrument responses.

Constructor
---------------------

- :class:`Decimator <tiskitpy.Decimator>`:  Set up a series of decimations corresponding
  to the list ``decimates``.

Properties
---------------------
- ``decimates``: the list of decimation factors
- ``decimation_factor``: total decimation (product of ``decimates``)
- ``verbose``: True if object is chatty.

Methods
---------------------

- :meth:`decimate <tiskitpy.Decimator.decimate>`: Decimate the data
- :meth:`update_inventory <tiskitpy.Decimator.update_inventory>`: Return inventory with decimated channels added
- :meth:`update_inventory_nslc <tiskitpy.Decimator.update_inventory_nslc>`: Return inventory with only
  specified networks, stations, channels and/or location updated
- :meth:`get_band_code <tiskitpy.Decimator.get_band_code>`: Return the band code for a given sampling rate
 
Command-line programs
---------------------

Use these programs' `-h` option for help

- ``tiskitpy_decimate_SDS``: Decimates the data in an SDS directory and updates
    the corresponding inventory file.

Example
---------------------

:ref:`tiskitpy.Decimator_example`
