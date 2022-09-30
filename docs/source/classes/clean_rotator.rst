.. _CleanRotator:

CleanRotator
=======================

Calculates the angle and azimuth to rotate the vertical channel to minimize
noise [#f1]_.

Detailed information is in :ref:`tiskit.CleanRotator`

The main methods are:

Constructor
---------------------

:`CleanRotator(stream,...)`: Calculate the CleanRotator object from
    a data stream

Properties
---------------------

:`angle`: angle (degrees) by which Z should be rotated
:`azimuth`: azimuth (degrees) by which Z should be rotated

Methods
---------------------

:`apply(stream, horiz_too=False)`: apply the rotation to the given stream
:`tfs()`: return the transfer functions equivalent to the rotation. *NOT
    FULLY TESTED/VALIDATED*

Example
---------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import CleanRotator

.. [#f1]  This only rotates
about 2 axes, so rotation of X and Y may lose their initial orientation (*need
to change to 3-value rotation that preserves X and Y azimuths*)]