TimeSpans class
=======================

A class to cut out "undesired" data.

Detailed information is in :ref:`tiskit.TimeSpans`

The main methods are:

Constructor
----------------------

:`TimeSpans(spans)`: create a TimeSpans object using the
    provided lists of [start_time, end_time]s.
:`TimeSpans.from_eqs(start_time, end_time, ...)`: create a TimeSpans
    object from an earthquake catalog.  If the catalog is not provided, will
    download from USGS.

Properties
----------------------

:`start_times`: returns list of start_times
:`end_times`: returns list of end_times
:`spans`: returns list of spans

Methods
----------------------

Modify a Trace or Stream
^^^^^^^^^^^^^^^^^^^^^^^^^

:`cutout(Stream or Trace)`: cuts out data from the Stream or Trace (using
    `obspy` `.cutout()` method
:`interp(Stream or Trace)`: linearly interpolate values within the time spans
    from their value at the span start to their value at the span end
:`zero(Stream or Trace)`: set values within the time spans to zero

Other
^^^^^^^^^^^^^^^^^^^^^^^^^

:`append(new_TimeSpans)`: appends two `.TimeSpan` objects
:`invert(start_time, end_time)`: return inverted time spans
:`has_zeros(starttime, endtime)`: does the given time range intersect any of
    the TimeSpans?
:`plot(...)`: plot the total time range, stream or trace, higlighting
    the time spans.

Example
----------------------

.. code-block:: python

  >>> from tiskit import TimeSpans
  
  >>> ts = TimeSpans([['2022-01-01T01', '2022-01-05'], ['2022-02-01', '2022-02-20']])

  >>> print(ts)
  TimeSpans: start             |            end
  =============================+===============================
   2022-01-01T01:00:00.000000Z | 2022-01-05T00:00:00.000000Z
   2022-02-01T10:00:00.000000Z | 2022-02-20T00:00:00.000000Z

  >>> print(ts.invert('2021-01-01', '2022-03-01'))
  TimeSpans: start             |            end
  =============================+===============================
   2021-01-01T00:00:00.000000Z | 2022-01-01T01:00:00.000000Z
   2022-01-05T00:00:00.000000Z | 2022-02-01T10:00:00.000000Z
   2022-02-20T00:00:00.000000Z | 2022-03-01T00:00:00.000000Z

  >>> print(ts.invert('2022-01-28', '2022-02-10'))
  TimeSpans: start              |            end
  ==============================+===============================
   2022-01-28T00:00:00.000000Z  | 2022-02-01T10:00:00.000000Z 

  >>> print(ts.invert(None, None))
  TimeSpans: start             |            end
  =============================+===============================
   2022-01-05T00:00:00.000000Z | 2022-02-01T10:00:00.000000Z

