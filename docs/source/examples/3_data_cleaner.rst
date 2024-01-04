.. _tiskitpy.DataCleaner_example:

==============================
DataCleaner example code
==============================

from ``tiskitpy/_examples/3_DataCleaner.py``

.. code-block:: python

    from obspy.core.stream import read
    from obspy.core.inventory import read_inventory
    from tiskitpy import DataCleaner, SpectralDensity

    # Read data and inventory
    stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
    inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')

    # Calculate and apply Datacleaner from/to stream, subtracting *1, then *2, then *H
    dc = DataCleaner(stream, ['*1','*2','*H'])

.. code-block:: none

    [INFO] z_threshold=3 rejected 3 of 84 windows (4%)
    [INFO] z_threshold=3 rejected 3 of 84 windows (4%)
    [INFO] z_threshold=3 rejected 3 of 84 windows (4%)
    [INFO] z_threshold=3 rejected 3 of 84 windows (4%)

.. code-block:: python

    stream_cleaned = dc.clean_stream(stream)


.. code-block:: none

    [INFO] Correcting traces in the frequency domain

.. code-block:: python

    # Construct and plot a stream with original and cleaned Z channels
    z_compare = stream.select(channel='*Z') + stream_cleaned.select(channel='*Z')
    z_compare.plot(outfile='3_DataCleaner_timeseries.png')
    print(z_compare)

.. code-block:: none

    2 Trace(s) in Stream:
    XS.S11D..LHZ | 2016-12-10T23:59:59.992583Z - 2016-12-11T23:59:59.992583Z | 1.0 Hz, 86401 samples
    XS.S11D..LHZ | 2016-12-10T23:59:59.992583Z - 2016-12-11T23:59:59.992583Z | 1.0 Hz, 86401 samples

.. code-block:: python

    # Tag trace seed_ids with cleaning information so that it shows up in obspy plots
    z_compare_tagged = CS.seedid_tag(z_compare)
    print(z_compare_tagged)

.. code-block:: none

    2 Trace(s) in Stream:
    XS.S11D..LHZ       | 2016-12-10T23:59:59.992583Z - 2016-12-11T23:59:59.992583Z | 1.0 Hz, 86401 samples
    XS.S11D.-1-2-H.LHZ | 2016-12-10T23:59:59.992583Z - 2016-12-11T23:59:59.992583Z | 1.0 Hz, 86401 samples

.. code-block:: python

    z_compare_tagged.plot()

.. image:: images/3_DataCleaner_tagged_timeseries.png
   :width: 564
   

.. code-block:: python

    # compare spectral densities
    # No need to tag seed_ids, tiskitpy plot() methods do it automatically
    sd_compare = SpectralDensity.from_stream(z_compare, inv=inv)
    sd_compare.plot(overlay=True, outfile='3_DataCleaner_sd_overlay.png')

.. image:: images/3_DataCleaner_sd_overlay.png
   :width: 564
   
   
