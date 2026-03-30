.. obsinfo-test documentation master file, created by
   sphinx-quickstart on Mon Jul 19 11:50:58 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TiSKitPy documentation
========================================

tiskitpy is the TIme Series data processing toolKIT, 
used for cleaning data series and calculating
the coherencies and frequence response functions between them.
Most of the algorithms are based on `BP2010`_.

Classes, functions and command_line programs are listed in the
:ref:`Overview`

The cleaning algorithms can create multiple instances of the 'same' data stream.
In order to keep track of these instances, tiskitpy stores clean_sequence
information, as described in :ref:`clean_sequences`

.. toctree::
  :maxdepth: 2
  :caption: Table of Contents:
  :glob:
  
  overview
  install
  classes
  functions
  periodic_transients
  avoid_eqs
  examples
  clean_sequences
  programmers


:ref:`genindex`

.. [BP2010] Bendat J. S. and A. G. Piersol (1986), Random Data:
    Analysis and Measurement Procedures, 566 pp.
