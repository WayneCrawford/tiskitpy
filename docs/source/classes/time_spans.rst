.. _TimeSpans:

TimeSpans
=======================

A class to cut out "undesired" data.

Detailed information is in :ref:`tiskitpy.TimeSpans`

The main methods are:

Constructor
----------------------

:``TimeSpans(spans)``: create a TimeSpans object using the
    provided lists of [start_time, end_time]s.
:``TimeSpans.from_eqs(start_time, end_time, ...)``: create a TimeSpans object from an earthquake catalog.  If the catalog is not provided, will download from USGS.

Properties
----------------------

:``start_times``: returns list of start_times
:``end_times``: returns list of end_times
:``spans``: returns list of spans

Methods
----------------------

Modify a Trace or Stream
^^^^^^^^^^^^^^^^^^^^^^^^^

:``cutout(Stream or Trace)``: cuts out data from the Stream or Trace (using
    ``obspy`` ``.cutout()`` method
:``interp(Stream or Trace)``: linearly interpolate values within the time spans from their value at the span start to their value at the span end
:``zero(Stream or Trace)``: set values within the time spans to zero

Other
^^^^^^^^^^^^^^^^^^^^^^^^^

:``append(new_TimeSpans)``: appends two :ref:`TimeSpans`` objects
:``invert(start_time, end_time)``: return inverted time spans
:``has_zeros(starttime, endtime)``: does the given time range intersect any of
    the TimeSpans?
:``plot(...)``: plot the total time range, stream or trace, higlighting
    the time spans.

Example
----------------------


see :ref:`tiskitpy.TimeSpans_example`
