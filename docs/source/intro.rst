
=============================
tiskitpy
=============================

A **TI**me **S**eries data processing tool**KIT**, designed to clean
data series and calculate the coherencies and frequence response functions between
them.  Most of the algorithms are based on [_BP2010].

Classes, functions and command_line programs are listed in the
:ref:`Functional overview`

The cleaning algorithms can create multiple instances of the 'same' data stream.
In order to keep track of these instances, tiskitpy stores clean_sequence
information, as described in :ref:`clean_sequences`

.. toctree::
   :maxdepth: 1

   intro/functional_overview
   intro/clean_sequences

.. [BP2010] Bendat J. S. and A. G. Piersol (1986), Random Data:
    Analysis and Measurement Procedures, 566 pp.