.. _CleanRotator:

CleanRotator
=======================

Calculates the angle and azimuth to rotate the vertical channel to minimize
noise.

**BUG: The code only rotates about two axes, so it is not appropriate for
rotating the horizontal axes, which could lose their orientation**

Detailed information is in :ref:`tiskitpy.CleanRotator`

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


see :ref:`tiskitpy.CleanRotator_example`
