Decimator class
=======================

Decimate data and put decimation information in instrument responses.

Detailed information is in :ref: `tiskit.Decimator`

The main methods are:

Constructor
---------------------

:`Decimate(decimates)`: Set up a series of decimations according to the list
    `decimates`, which can contain integers between 2 and 7.

Properties
---------------------
:`decimates`: returns the decimates list
:`decimation_factor`: the total decimation of the object (the product of
    the `decimates` values`)
:`verbose`: True if object is in verbose mode.


Methods
---------------------

:`decimate(Stream or Trace)`: Decimate the data
:`get_band_code(in_band_code, sample_rate)`: return the band code for a given
    sample rate.  `in_band_code` simply indicates if the instrument is broadband
    ("B" or other broadband band codes) or short period ("S" or other short-period
    band codes).
:`update_inventory(inv, ...)`: Return an inventory with decimated channels added
:`update_inventory_from_nslc(inv ...)`: Return an inventory with only the specified
 network, station, channel, location(s) updated
 

Example
---------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import Decimator
  
