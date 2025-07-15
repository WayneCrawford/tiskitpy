.. _TimeSpans:

TimeSpans
=======================

A class to cut out "undesired" data.


Constructor
----------------------

- :class:`TimeSpans <tiskitpy.TimeSpans>`: create a TimeSpans object using
  lists of [start_time, end_time]s.
- :meth:`TimeSpans.from_eqs <tiskitpy.TimeSpans.from_eqs>`: create an object
  of TimeSpans to avoid, based on an earthquake catalog.

Properties
----------------------

- ``start_times``: returns list of span start_times
- ``end_times``: returns list of span end_times
- ``spans``: returns list of spans

Methods
----------------------

Modify a Trace or Stream
^^^^^^^^^^^^^^^^^^^^^^^^^
TimeSpans
- :meth:`cutout <tiskitpy.ResponseFunctions.cutout>`: cuts out data in the time spans
- :meth:`interp <tiskitpy.TimeSpans.interp>`: linearly interpolate values within the time
  spans from their value at the span start to their value at the span end
- :meth:`zero <tiskitpy.TimeSpans.zero>`: set values within the time spans to zero

Other
^^^^^^^^^^^^^^^^^^^^^^^^^

- :meth:`combine <tiskitpy.TimeSpans.combine>`: combines two :ref:`TimeSpans`` objects
- :meth:`invert <tiskitpy.TimeSpans.invert>`: return inverted time spans
- :meth:`has_zeros <tiskitpy.TimeSpans.has_zeros>`: Returns a bool indicating if the given time
  range intersects any of the TimeSpans.
- :meth:`plot <tiskitpy.TimeSpans.plot>`: Make a stream or trace plot with highlighted time spans.

Example
----------------------

:ref:`tiskitpy.TimeSpans_example`
