CleanRotator class
=======================

Calculates angle and azimuth to rotate the vertical channel to minimize
vertical noise [#f1]_. 

Detailed information is in :ref:`tiskit.CleanRotator`

The main methods are:

Constructor
---------------------

:`CleanRotator(...)`:

Properties
---------------------

`angle`: angle by which the axes should be rotated
`azimuth`: azimuth by which the axes should be rotated

Get Methods
---------------------

`tfs()`: Return the Z-1 and Z-2 transfer functions equivalent to the
    rotation. *NOT FULLY TESTED*

Other
---------------------

:`apply(Stream, horiz_too=bool)`: Applies to vertical (and optionally
    horizontal, but see Note above) channels

Example
---------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import CleanRotator
  
.. [#f1]  THIS IS OK FOR THE VERTICAL CHANNEL, BUT WILL
   MIS-ALIGN THE HORIZONTAL CHANNELS, SHOULD CHANGE TO 3-value rotation

