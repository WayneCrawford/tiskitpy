.. _PeriodicTransients

=============================
Removing Periodic Transients
=============================

Some data have periodic transients (instrument relevels, disk drive writing...)
that can interfere with signal processing.
The ``PeriodicTransient`` class lets you calculate and remove these transients,
if they are truly periodic and always have the same shape.
This is a labor-intensive, manual module, so if you use it
on a dataset, please share your results!
Properties and methods are:

Constructor
---------------------

:class:`PeriodicTransient <tiskitpy.PeriodicTransient>`


Properties
---------------------

- Input by the constructor 
    - ``name`` (str): name of this periodic transient (e.g., 'hourly')
    - ``period`` (float): seconds between each transient
    - ``dp`` (float): how many seconds to change the period by when testing for better values
    - ``clips`` (tuple): clip values outside of this range (low, high).
      Should include the max range of the transient
    - ``transient_starttime`` (`UTCDateTime``): onset time of earliest transient.
    - ``freq_HP`` (float or bool): highpass corner frequency used for training
      and matching.  If True, sets to 1/period.  If False, do not cut off
      low frequencies
    - ``freq_LP`` (float): lowpass corner frequency used for training and
      matching (0.05 is a good value to remove microseisms)
- Calculated by ``PeriodicTransient.calc_transient()``
    - ``transient_model`` (): Model of the periodic transient.
    - ``dirac_comb`` (): offsets of transients from the first one, in seconds
    - ``n_transients_used`` (): number of transients used to make the model
    - ``tm`` (): transient starting 1 sample earlier
    - ``tp``(): transient starting 1 sample later


Methods
---------------------

- :meth:`calc_timing <tiskitpy.PeriodTransient.calc_timing>`: Interactively
  calculate transient parameters
- :meth:`calc_transient <tiskitpy.PeriodTransient.calc_transient>`: Calculate
  the shape of the transient
- :meth:`remove_transient <tiskitpy.PeriodTransient.remove_transient>`: Remove
  the transient from a data trace

Example
---------------------

:ref:`tiskitpy.PeriodicTransient_example`