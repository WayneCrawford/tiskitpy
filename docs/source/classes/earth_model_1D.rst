.. _EarthModel1D:

EarthModel1D
=======================

1D Earth model, used for compliance calculations


Constructor
---------------------

- :class:`tiskitpy.compliance.EarthModel1D`:

Properties
---------------------

- ``thicks`` (:class:`numpy.ndarray`): layer thicknesses (m)
- ``rhos`` (:class:`numpy.ndarray`): layer densities (kg/m^3)
- ``vps`` (:class:`numpy.ndarray`): layer compressional velocities (m/s)
- ``vss`` (:class:`numpy.ndarray`): layer shear velocities (m/s)

Methods
---------------------

- :meth:`plot <tiskitpy.compliance.EarthModel1D.plot>`: Plot the model
- :meth:`calc_ncompl <tiskitpy.compliance.EarthModel1D.calc_ncompl>`: Calculate normalized compliance

Example
---------------------

None for now
