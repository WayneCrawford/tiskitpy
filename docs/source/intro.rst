
=============================
tiskitpy
=============================

TIme Series data processing toolKIT, for cleaning data series and calculating
the coherencies and frequence response functions between them.
Most of the algorithms are based on `BP2010`_.

Classes, functions and command_line programs are listed in the
:ref:`Module overview`

The cleaning algorithms can create multiple instances of the 'same' data stream.
In order to keep track of these instances, tiskitpy stores clean_sequence
information, as described in :ref:`clean_sequences`

.. toctree::
   :maxdepth: 1

   intro/module_overview
   intro/clean_sequences

.. [BP2010] Bendat J. S. and A. G. Piersol (1986), Random Data:
    Analysis and Measurement Procedures, 566 pp.