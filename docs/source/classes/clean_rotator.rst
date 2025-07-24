.. _CleanRotator:

CleanRotator
=======================

Calculates the angle and azimuth to rotate the vertical channel to minimize
noise.

**BUG: The code only rotates about two axes, so it is not appropriate for
rotating the horizontal axes, which could lose their orientation**

Constructor
---------------------

- :class:`CleanRotator <tiskitpy.CleanRotator>`: Calculate the CleanRotator object from
  a data stream

Properties
---------------------

- ``angle`` (float): angle (degrees) by which Z should be rotated
- ``azimuth`` (float): azimuth (degrees) by which Z should be rotated
- ``variance_reduction`` (float): amount by which variance was reduced during calculation
    (between 0. and 1.)
- ``avoided_spans`` (:class:`tiskitpy.TimeSpans`): time spans avoided during calculation

Methods
---------------------

- :meth:`apply <tiskitpy.CleanRotator.apply>`: Apply the rotation to the given stream
- :meth:`tfs <tiskitpy.CleanRotator.tfs>`: Return the transfer functions equivalent to the rotation. *NOT
  VALIDATED*

Example
---------------------

:ref:`tiskitpy.CleanRotator_example`
