.. _tiskit.Decimate_SDS_example:

==============================
tiskit_decimate_SDS example
==============================

tiskit_decimate_SDS is a console script to decimate data and create associated
metadata, if the data is in an SDS directory and the metadata are in a
StationXML file.
The new data are added to the SDS directory, whereas the new and old metadata
are saved in a new StationXML file.

The program acts on data with a given sampling rate, so that the outputs will
all have the same sampling rate

The `-h` option provides basic documentation

.. code-block:: bash

    (obspy) friska:instrumentation 554> tiskit_decimate_SDS -h
    usage: tiskit_decimate_SDS [-h] [--of OUTPUT_FILE] [-q] SDS_root inv_file input_sample_rate {2,3,4,5,6,7} [{2,3,4,5,6,7} ...]

    Insert decimated channels and create new StationXML file

    positional arguments:
      SDS_root           SDS root directory
      inv_file           StationXML file
      input_sample_rate  Process only channels having this sample rate
      {2,3,4,5,6,7}      Sequence of decimation factors to use)

    optional arguments:
      -h, --help         show this help message and exit
      --of OUTPUT_FILE   Output StationXML filename (default = infile.replace('.xml', '_decim.xml')
      -q, --quiet        Suppress information messages


Here's an example with the data in a directory named "/Volumes/Data/SDS" and the
metadat in a file named ALPARRAY.INSU-IPGP.station.xml

.. code-block:: bash

    (obspy) friska:metadata 568> tiskit_decimate_SDS /Volumes/Data/SDS/ ALPARRAY-OBS.INSU-IPGP.station.xml 62.5 5 5 4
    INFO:root:output sampling rate will be 0.625 sps, band code will be L
    INFO:root:Working on year 2017
    INFO:root:    Working on net Z3
    INFO:root:        Working on station A401A
    INFO:root:            Working on channel BDH.D
    INFO:root:              Creating output channel dir "LDH.D"
    INFO:root:              22 files to process
    INFO:root:            Working on channel BH1.D
    INFO:root:              Creating output channel dir "LH1.D"
    INFO:root:              22 files to process
    INFO:root:            Working on channel BH2.D
    INFO:root:              Creating output channel dir "LH2.D"
    INFO:root:              22 files to process
    INFO:root:            Working on channel BHZ.D
    INFO:root:              Creating output channel dir "LHZ.D"
    INFO:root:              22 files to process
    INFO:root:        Working on station A413A
    INFO:root:            Working on channel BDH.D
    INFO:root:              Creating output channel dir "LDH.D"
    INFO:root:              192 files to process
    INFO:root:            Working on channel BH1.D
    INFO:root:              Creating output channel dir "LH1.D"
    INFO:root:              192 files to process
    INFO:root:            Working on channel BH2.D
    INFO:root:              Creating output channel dir "LH2.D"
    INFO:root:              192 files to process
    INFO:root:            Working on channel BHZ.D
    INFO:root:              Creating output channel dir "LHZ.D"
    INFO:root:              192 files to process
    INFO:root:        Working on station A416A
    INFO:root:          Working on channel BDH.D
    INFO:root:            Creating output channel dir "LDH.D"
    INFO:root:            189 files to process
    ...
    
The SDS directory's aborescence before:

.. code-block:: bash

    (base) friska:2017-18.AlpArray 508> tree -L 4 SDS
    SDS
    ├── 2017
    │   └── Z3
    │       ├── A401A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   ├── BHZ.D
    │       ├── A413A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   ├── BHZ.D
    │       ├── A416A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   ├── BHZ.D
    │       ├── A419A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   ├── BHZ.D
    │       ├── A422A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   └── BHZ.D
    │       ├── A425A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   └── BHZ.D
    │       └── A429A
    │           ├── BDH.D
    │           ├── BH1.D
    │           ├── BH2.D
    │           └── BHZ.D
    └── 2018
        └── Z3
            ├── A401A
            │   ├── BDH.D
            │   ├── BH1.D
            │   ├── BH2.D
            │   └── BHZ.D
            ├── A413A
            │   ├── BDH.D
            │   ├── BH1.D
            │   ├── BH2.D
            │   └── BHZ.D
            ├── A416A
            │   ├── BDH.D
            │   ├── BH1.D
            │   ├── BH2.D
            │   └── BHZ.D
            ├── A419A
            │   ├── BDH.D
            │   ├── BH1.D
            │   ├── BH2.D
            │   └── BHZ.D
            ├── A422A
            │   ├── BDH.D
            │   ├── BH1.D
            │   ├── BH2.D
            │   └── BHZ.D
            ├── A425A
            │   ├── BDH.D
            │   ├── BH1.D
            │   ├── BH2.D
            │   └── BHZ.D
            └── A429A
                ├── BDH.D
                ├── BH1.D
                ├── BH2.D
                └── BHZ.D

And after:
    
.. code-block:: bash

    (base) friska:2017-18.AlpArray 508> tree -L 4 SDS
    SDS
    ├── 2017
    │   └── Z3
    │       ├── A401A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   ├── BHZ.D
    │       │   ├── LDH.D
    │       │   ├── LH1.D
    │       │   ├── LH2.D
    │       │   └── LHZ.D
    │       ├── A413A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   ├── BHZ.D
    │       │   ├── LDH.D
    │       │   ├── LH1.D
    │       │   ├── LH2.D
    │       │   └── LHZ.D
    │       ├── A416A
    │       │   ├── BDH.D
    │       │   ├── BH1.D
    │       │   ├── BH2.D
    │       │   ├── BHZ.D
    │       │   ├── LDH.D
    │       │   ├── LH1.D
    │       │   ├── LH2.D
    │       │   └── LHZ.D
    ...
    
A new StationXML file is created, named ALPARRAY.INSU-IPGP_decim.xml