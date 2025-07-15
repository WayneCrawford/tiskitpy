.. _CleanedStream:

CleanedStream
=======================

Subclass of :class:`obspy.Stream` that embeds processing steps

Constructor
---------------------

- :class:`CleanedStream <tiskitpy.CleanedStream>``: Works exactly like ``obspy.Stream``

Properties
---------------------

Same as ``obspy.Stream``

Methods
---------------------

Modified from ``obspy.Stream``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- :meth:`__str__ <tiskitpy.CleanedStream.__str__>`: outputs the :ref:`tiskitpy_id` instead of the seed_id
- :meth:`plot <tiskitpy.CleanedStream.plot>`: uses the :ref:`tiskitpy_id` instead of the seed_id
- :meth:`select <tiskitpy.CleanedStream.select>`: selects on the :ref:`tiskitpy_id` if the seed_id doesn't work

New
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- :meth:`tag <tiskitpy.CleanedStream.tag>`: tags the stream with the given seed_id or transformation code

Example
---------------------

None for now
