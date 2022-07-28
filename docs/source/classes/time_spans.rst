TimeSpans class
=======================

A class to cut out "undesired" data.  Mostly used by the other classes, you
may never use this!

Detailed information is in :ref: `tiskit.TimeSpans`

The main methods are:

Constructor
----------------------

:`TimeSpans(start_times, end_times)`: create a TimeSpans object using the
    provided lists of start_times and end_times.
:`TimeSpans.from_eqs(start_time, end_time, ...)`: create a TimeSpans
    object from an earthquake catalog.  If the catalog is not provided, will
    download from USGS.

Properties
----------------------

`start_times`: returns list of start_times
`end_times`: returns list of end_times

Set Methods
----------------------

:`append(new_TimeSpans)`: appends two `.TimeSpan` objects


Trace/Stream modifying
----------------------

:`cutout(Stream or Trace)`: cuts out data from the Stream or Trace (using
    `obspy` `.cutout()` method
:`interp(Stream or Trace)`: linearly interpolate values within the time spans
    from their value at the span start to their value at the span end
:`zero(Stream or Trace)`: set values within the time spans to zero

Other
----------------------

:`has_zeros(starttime, endtime)`: does the given time range intersect any of
    the TimeSpans?
:`plot()`: plot a representation of the total time span, with individual
    time-spans highlighted in yellow

Example
----------------------

.. code-block:: python

  from obspy.core.stream import read
  from obspy.core.inventory import read_inventory
  from tiskit import CleanRotator
  
