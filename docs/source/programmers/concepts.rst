*******************************
Concepts
**************
*****************

.. py:module:: tiskitpy

Clean sequences
=========================

Clean sequences are lists of the operations that have been performed on a data stream.
They are mostly handled by the :py:class:`tiskitpy.CleanSequence` class

``tiskitpy_id`` is a textual representation of the clean sequence list, with each
element separated by "-" and most of the processing codes shortened to their
first or last unambiguous letters.
The ``tiskitpy_id``` can be stuffed into a Trace's "location" code for
plotting purposes, but it shouldn't be stored there permanently because it can
hinder operations that are based on seed ids, such as finding an instrument
response in an inventory or comparing two channels.
The :py:method:`CleanSequence.tag` and :py:method:`tiskitpy.CleanSequence.untag` methods
are part of a suite of :py:class:`CleanSequence` methods to generate and remove
``tiskitpy_id`` from strings, seed_ids, etc.