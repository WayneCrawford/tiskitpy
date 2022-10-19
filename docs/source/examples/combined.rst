==============================
Example code
==============================

Here as an example of what you can do with TiSKit.  It uses several tiskit classes to clean seafloor seismological data.

.. code-block:: python

    from pathlib import Path
    import fnmatch

    import numpy as np
    from obspy.clients.fdsn import Client
    from obspy.core import UTCDateTime
    from obspy.core.stream import read
    from obspy.core.inventory import read_inventory
    from matplotlib import pyplot as plt

    from tiskit import Decimator, CleanRotator, SpectralDensity, DataCleaner

    client_address = 'IRIS'
    net = 'Z5'
    sta = 'BB650'  # Trillium T240. 100 sps
    loc = '*'
    cha = 'HH*,HDH'
    starttime = UTCDateTime('2014-01-01T') 
    endtime = UTCDateTime('2014-01-02T')

    # READ WAVEFORM DATA
    fname=f'data/{net}.{sta}.{starttime.strftime("%Y%m%d")}'
    client = Client(client_address)
    if Path(fname).is_file(): 
        print(f'reading data from file "{fname}"')
        stream = read(fname, "MSEED")
    else:
        print('reading data from FDSN server (should take a long time...)')
        stream = client.get_waveforms(network=net, station=sta, channel=cha,
                                      location=loc, starttime=starttime,
                                      endtime=endtime)
        Path(fname).parent.mkdir(exist_ok=True)
        stream.write(fname, "MSEED")  # COMMENT OUT ONCE YOU'VE READ FROM FDSN!!!!!

    # READ INVEMTORY (INSTRUMENT RESPONSE, etc)
    fname=f'data/{net}.{sta}.{starttime.strftime("%Y%m%d")}.station.xml'    
    if Path(fname).is_file(): 
        print(f'reading inventory from file "{fname}"')
        inv = read_inventory(fname, "STATIONXML")
    else:
        print('reading inventory from FDSN server (could take a long time...)')
        inv = client.get_stations(starttime=starttime, endtime=endtime,
                              network=net, station=sta, channel='*', location='*',
                              level='response')
        Path(fname).parent.mkdir(exist_ok=True)
        inv.write(fname, "STATIONXML")  # COMMENT OUT ONCE YOU'VE READ FROM FDSN!!!!!

    # DECIMATE DATA DOWN TO 1 Hz
    decim = Decimator([5, 5, 4])
    stream_decim = decim.decimate(stream)
    inv_decim = decim.update_inventory(inv, stream)

    sd_orig = SpectralDensity.from_stream(stream_decim, inv=inv_decim)

    # USE SIMPLE ROTATION TO REDUCE VERTICAL CHANNEL NOISE
    rotator = CleanRotator(stream_decim)
    rot_stream = rotator.apply(stream_decim)
    sd_rot = SpectralDensity.from_stream(rot_stream, inv=inv_decim)

    # USE TRANSFER FUNCTION BASED DATA CLEANER TO REDUCE VERTICAL CHANNEL NOISE
    dc = DataCleaner(rot_stream, ['*1', '*2', '*H'])
    # first clean the stream, then calculate the spectral density
    rot_stream_dc = dc.clean_stream(rot_stream)
    sd_rot_dc = SpectralDensity.from_stream(rot_stream_dc, inv=inv_decim)
    # directly calculate the spectral density, with the datacleaner as input
    sd_rot_sddc = dc.clean_stream_to_sdf(rot_stream, inv=inv_decim)

    # PLOT THE RESULTS
    fig, ax = plt.subplots()
    for sd, label in zip((sd_orig, sd_rot, sd_rot_dc, sd_rot_sddc),
                          ('original', 'rotated', 'rot + clean', 'rot+clean(sd)')
                        ):
        print(sd.channels)
        z_id = fnmatch.filter(sd.channels, '*.LHZ*')[0]
        ax.semilogx(sd.freqs, 10*np.log10(sd.autospect(z_id)), label=label)
    ax.legend()
    plt.show()
    print(stream)
    stream.plot()


The output plot should look like this:

.. image:: images/example_cleaning.png
   :width: 564