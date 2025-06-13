**************************************************************
Avoiding earthquakes and other noisy data
**************************************************************

.. py:module:: tiskitpy

All of the class constructors that use time series data can remove unwanted
time spans, including the noisy time after large earthquakes.

The constructors are :py:meth:`SpectralDensity.from_stream`, 
:py:meth:`DataCleaner` and :py:meth:`CleanRotator`.
They share the following input parameters:

- ``avoid_spans`` (:py:class:`TimeSpans`)
- ``remove_eqs`` (``str`` or ``bool``)

and provide the following properties:

- ``avoided_spans```: time spans NOT used in processing
- ``used_spans```: inverse of ``avoided_spans``

Default behavior
=================

By default ``remove_eqs=True`` for each constructors, so they
will run :py:meth:`TimeSpans.from_eqs` to generate time spans to avoid based
on the USGS earthquake catalog (online or saved in a local file by a previous call).
The ``z_threshold`` parameter of
:py:meth:`SpectralDensity.from_stream` sets a :py:meth:`scipy.stats.zscore`
threshold for eliminating 
outliers in the calculated spectra: this may add to the time spans that are
"avoided".

The object property ``avoid_spans`` indicates all of the time spans avoided:
earthquake spans, spans added using ``avoid_spans`` and any time spans rejected
using the z_threshold (:py:class:`SpectralDensity` only)


Adding/modifying time spans to avoid
====================================

If the user sets ``remove_eqs`` to a string, this must correspond to the name
of a local file in QuakeML format.

The user can specify ``avoid_spans`` for other time spans to avoid.

The user can also specify another earthquake catalog to avoid (using the
default input paramters for :py:meth:`TimeSpans.from_eqs`) by setting
``remove_eqs=<filename>``, where ``<filename>`` is the name of a QuakeML file
in the current working directory.


If the user wants to calculate time spans to avoid based on earthquakes, but
using different parameters than the default values, they should run
:py:meth:`TimeSpans.from_eqs` themselves and enter the output to
the ``avoid_spans`` parameter.

If the user wants to calculate time spans to avoid themselves, they need
to convert these spans to :py:class:`TimeSpans` format and enter using
``avoid_spans``

Using the same time spans for different objects
=============================================================================

If you want to avoid exactly the same time spans for each class, you
should retrieve the ``avoided_spans`` parameter from one of them (preferably
:py:meth:`SpectralDensity.from_stream` since it uses the z-score to add time
spans to avoid) and use it as input to ``avoid_spans`` on the others.
