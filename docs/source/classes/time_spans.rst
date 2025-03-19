.. _TimeSpans:

TimeSpans
=======================

A class to cut out "undesired" data.

Detailed information is in :ref:`tiskitpy.TimeSpans`

The main methods are:

Constructor
----------------------

- ``TimeSpans(spans)``: create a TimeSpans object using the
  provided lists of [start_time, end_time]s.
- ``TimeSpans.from_eqs(start_time, end_time, ...)``: create an object
  of TimeSpans to avoid, based on an earthquake catalog.
  If the catalog is not provided, ir will be downloaded from USGS.

Properties
----------------------

- ``start_times``: returns list of span start_times
- ``end_times``: returns list of span end_times
- ``spans``: returns list of spans

Methods
----------------------

Modify a Trace or Stream
^^^^^^^^^^^^^^^^^^^^^^^^^

- ``cutout(Stream or Trace)``: cuts out data in the time spans
- ``interp(Stream or Trace)``: linearly interpolate values within the time
  spans from their value at the span start to their value at the span end
- ``zero(Stream or Trace)``: set values within the time spans to zero

Other
^^^^^^^^^^^^^^^^^^^^^^^^^

- ``combine(new_TimeSpans)``: combines two :ref:`TimeSpans`` objects
- ``invert(start_time, end_time)``: return inverted time spans
- ``has_zeros(starttime, endtime)``: Returns a bool indicating if the given time
  range intersects any of the TimeSpans.
- ``plot(...)``: Make a stream or trace plot with highlighted time spans.

Example
----------------------

:ref:`tiskitpy.TimeSpans_example`
