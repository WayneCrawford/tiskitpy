.. _CleanedStream:

CleanedStream
=======================

Subclass of :class:`obspy.Stream` that embeds processing steps

Detailed information is in :ref:`tiskitpy.CleanedStream`

The main methods are:

Constructor
---------------------

- ``CleanedStream(stream)``: Works exactly like ``obspy.Stream``

Properties
---------------------

Same as ``obspy.Stream``

Methods
---------------------

Modified from ``obspy.Stream``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``__str__()``: outputs the :ref:`tiskitpy_id` instead of the seed_id
- ``plot()``: uses the :ref:`tiskitpy_id` instead of the seed_id
- ``select()``: selects on the :ref:`tiskitpy_id` if the seed_id doesn't work

New
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``tag()``: tags the stream with the given seed_id or transformation code

Example
---------------------

None for now
