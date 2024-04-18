.. _tiskitpy.TimeSpans_example:

==============================
TimeSpans example code
==============================

.. code-block:: python

    from obspy.core.stream import read
    from obspy.core.inventory import read_inventory
    from tiskitpy import SpectralDensity, TimeSpans

    # Read the data and metadata, select the Z channel
    stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
    inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
    Z = stream.select(channel='*Z')

    # Select a time span to avoid (has an earthquake)
    ts = TimeSpans([['2016-12-11T16:50', '2016-12-11T17:10']])
    ts.plot(Z)

.. image:: images/6_TimeSpans_tsplot.png
   :width: 564

.. code-block:: python

    # Zero the data directly in the stream and plot the result
    Zs = ts.zero(Z)  # Returns a CleanedStream object
    Zs += Z
    Zs.plot()

.. image:: images/6_TimeSpans_zeroed.png
   :width: 564

.. code-block:: python

    # Compare SpectralDensity objects calculated with and without this time span
    # (note: cancelled z_threshold selection to avoid automatically eliminating
    # the EQ!)
    kwargs={'inv': inv, 'z_threshold': None}
    sd =   SpectralDensity.from_stream(Z, **kwargs)
    sd_z = SpectralDensity.from_stream(Z, avoid_spans=ts, **kwargs)
    SpectralDensity.plots([sd, sd_z])

.. image:: images/6_TimeSpans_spect_both.png
   :width: 564

