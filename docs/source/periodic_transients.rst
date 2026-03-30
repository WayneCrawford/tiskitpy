
=============================
Removing Periodic Transients
=============================

Some data has periodic transients (instrument relevels, disk drive writing...)
that can interfere with signal processing.  The `PeriodicTransients` class
allows one to minimize these transients, if they are truly periodic and always
have the same shape.  This is a labor-intensive, manual module, so if you use it
on a dataset, please share your results!  Properties and methods are:

Constructor
---------------------

- :class:`PeriodicTransient`:


Properties
---------------------

- Input by the constructor 
    - name (str): name of this periodic transient (e.g., 'hourly')
    - period (float): seconds between each transient
    - dp (float): how many seconds to change the period by when testing for better values
    - clips (tuple): clip values outside of this range (low, high). Should include the max range of the transient
    - transient_starttime (`UTCDateTime``): onset time of earliest transient.
- Calculated by `PeriodicTransient.calc_transient()`
    - transient_model (): Model of the periodic transient.  Created by `PeriodicTransient.calc_transient()`
    - dirac_comb (): offsets of transients from the first one, in seconds
    - n_transients_used (): number of transients used to make the model
    - tm (): transient starting 1 sample earlier
    - tp(): transient starting 1 sample later


Methods
---------------------

- :meth:`calc_timing <tiskitpy.PeriodTransient.calc_timing>`: Calculate the transient period
- :meth:`calc_transient <tiskitpy.PeriodTransient.calc_timing>`: Calculate the shape of the transient
- :meth:`calc_timing <tiskitpy.PeriodTransient.remove_transient>`: Remove the transient from a data trace

Example
---------------------

No example yet