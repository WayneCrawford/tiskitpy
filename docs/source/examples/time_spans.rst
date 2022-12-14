.. _tiskit.TimeSpans_example:

==============================
TimeSpans example code
==============================

.. code-block:: python

    from obspy.core.stream import read
    from obspy.core.inventory import read_inventory
    from tiskit import SpectralDensity, TimeSpans

    stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
    inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
    Z = stream.select(channel='*Z')
    ts = TimeSpans([['2016-12-11T16:50', '2016-12-11T17:10']])
    ts.plot(Z)

.. image:: images/TimeSpans_plot.png
   :width: 564

.. code-block:: python

    Zs = Z.copy()
    Zs += ts.zero(Z)
    Zs.plot()

.. image:: images/TimeSpans_zeroplot.png
   :width: 564

.. code-block:: python

    # Use large z_threshold so that SpectralDensity doesn't automatically
    # reject high-spectral data
    sd = SpectralDensity.from_stream(Z, inv=inv, z_threshold=1000)
    sd_zeroed = SpectralDensity.from_stream(Z, inv=inv, z_threshold=1000,
                                        avoid_spans=ts)
    sd.plot()

.. image:: images/TimeSpans_prespect.png
   :width: 564

.. code-block:: python

    sd_zeroed.plot()

.. image:: images/TimeSpans_postspect.png
   :width: 564

