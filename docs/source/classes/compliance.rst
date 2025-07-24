.. _Compliance:

Compliance
=======================

Data class for seafloor compliance, plus useful static functions


Constructor
---------------------

- :class:`Compliance <tiskitpy.Compliance>`
- :meth:`Compliance.from_response_functions <tiskitpy.Compliance.from_response_functions`:
  The most common way to create a Compliance object
- :meth:`Compliance.from_file <tiskitpy.Compliance.from_file`: read compliance from a file

Properties
---------------------

- ``freqs`` (:class:`numpy.ndarray`): Frequencies (Hz)
- ``values`` (:class:`numpy.ndarray`): Normalized compliance values (1/Pa)
- ``uncertainties`` (:class:`numpy.ndarray`): Normalized compliance uncertainties (1/Pa)
- ``water_depth`` (float): water depth in meters
- ``noise_channel`` (str or None): If a str, the compliance comes from data
  and this is the channel on which noise was assumed to dominate.  If None,
  the compliance comes from a calculation.
- ``gravity_corrected`` (bool): Has data-estimated compliance been corrected for
  gravitational attraction terms?

Dependent properties
^^^^^^^^^^^^^^^^^^^^^^


Methods
---------------------

- :meth:`correct_gravity_terms <tiskitpy.Compliance.correct_gravity_terms>`: Correct gravity terms
- :meth:`write <tiskitpy.Compliance.write>`: Write compliance to a text file
- :meth:`plot <tiskitpy.Compliance.plot>`: Plot the compliance

Static Methods
---------------------
- :meth:`plot_compliance_stack <tiskitpy.Compliance.plot_compliance_stack>`:
  Plot stacked PSDs, coherence and frequency response function used for
  compliance
- :meth:`calc_norm_compliance <tiskitpy.Compliance.calc_norm_compliance>`:
  Return normaliezed compliance of a :class:`tiskitpy.compliance.EarthModel1D` object
- :meth:`gravd <tiskitpy.Compliance.gravd>`: Return linear ocean surface
  gravity wave wavenumbers
- :meth:`raydep <tiskitpy.Compliance.raydep>`: Return surface motion and
  stress of P-SV waves for a layered 1D earth model and wave slownesses
- :meth:`calc_compliance <tiskitpy.Compliance.calc_compliance>`: Return compliance
  (m/Pa) of a :class:`tiskitpy.compliance.EarthModel1D` object

Example
---------------------

:ref:`tiskitpy.Compliance_example`
