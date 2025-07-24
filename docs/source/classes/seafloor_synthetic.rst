.. _ComplianceNoise:

SeafloorSynthetic
=======================

Generates synthetic seismological signals based on seafloor low frequency noise
and infragravity wave signal and noise levels.

Uses :class:`PSDVals <tiskitpy.PSDVals` class objects, which specify a
power spectral density using dB values at different frequenceis

Constructor
---------------------

- :class:`SeafloorSynthetic <tiskitpy.SeafloorSynthetic>`:

Properties
---------------------

- ``water_depth`` (float): water_depth in meters
- ``Z_offset_angles`` (tuple): (angle, azimuth) that Z is offset from vertical (degrees)
- ``IG_m_seasurface`` (:class:`tiskitpy.PSDVals`): IG sea surface wave height (dB ref m)
- ``noise_pressure`` (:class:`tiskitpy.PSDVals`): Pressure sensor self-noise (dB ref Pa)
- ``noise_seismo`` (:class:`tiskitpy.PSDVals`): Seismometer self noise (dB ref m/s^2)
- ``noise_tilt_max`` (:class:`tiskitpy.PSDVals`): Maximum tilt noise (dB ref nrad)
- ``noise_tilt_direction_limits`` (tuple): Minimum and maximum tilt directions (degrees)
- ``noise_tilt_variance`` (float): dBs below noise_tilt_max for minimum tilt noise
- ``earth_model`` (:class: `tiskitpy.EarthModel1D`): 1D Earth model used to calculate compliance
- ``IG_freqstep`` (float): maximum frequency step allowed between seafloor IG wave values

Dependent properties
^^^^^^^^^^^^^^^^^^^^^^

- ``IG_Pa_seafloor``: Infragravity wave seafloor pressure PSD
- ``PSDs``: Dictionary of all PSDs
- ``Z_angle_factor_DBs``: Rotation of horizontal tilt noise onto Z channel
- ``compliance_accel``: PSD of compliance * IG pressure

Methods
---------------------

- :meth:`make_tilt_ts <tiskitpy.SeafloorSynthetic.make_tilt_ts>`: Make a simple tilt time series
- :meth:`norm_compliance <tiskitpy.SeafloorSynthetic.norm_compliance>`: Return normalized compliance of objects ``EarthModel``
- :meth:`plot <tiskitpy.SeafloorSynthetic.plot>`: Plot spectral representation of the noise sources
- :meth:`save_compliance <tiskitpy.SeafloorSynthetic.save_compliance>`: save ``self.earth_model``'s compliance to a CSV file
- :meth:`streams <tiskitpy.SeafloorSynthetic.streams>`: Return streams generated from the noise model

Example
---------------------

:ref:`tiskitpy.ComplianceNoise_example`
