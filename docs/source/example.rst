# Example code

Here as an example of what you can do with TiSKit.  It uses several tiskit classes to clean seafloor seismological data.

.. code-block:: python

    from obspy.clients.fdsn import Client
    from obspy.core import UTCDateTime
    from obspy.core.stream import read
    from obspy.core.inventory import read_inventory
    from matplotlib import pyplot as plt
    
    client_address = 'IRIS'
    net = 'Z5'
    sta = 'BB870'
    loc = '*'
    cha = 'BH*,BDH'
    start_time = UTCDateTime('2014-08-01T')
    end_time = UTCDateTime('2014-08-02T')
    
    # READ WAVEFORM DATA
    client = Client(client_address)
    # This takes the longest time because you are reading XX bytes over the internet
    stream = client.get_waveforms(network=net, station=sta, channel=cha, location=loc,
                                  starttime=start_time, endtime=end_time)
    # Your next runs will be faster if you save and read locally:
    fname=f'{net}.{sta}.{starttime.strftime("%Y%m%d")}'    
    stream.write(fname, "MSEED")  # COMMENT OUT ONCE YOU'VE READ FROM FDSN!!!!!
    stream = read(fname, "MSEED")
    
    # READ INVEMTORY (INSTRUMENT RESPONSE, etc)
    inv = client.get_stations(starttime=start_time, endtime=end_time,
                              network=net, station=sta, channel=cha, location=loc)
    # Your next runs will be faster if you save and read locally
    fname=f'{net}.{sta}.{starttime.strftime("%Y%m%d")}.station.xml'    
    inv.write(fname, "STATIONXML")  # COMMENT OUT ONCE YOU'VE READ FROM FDSN!!!!!
    inv = read_inventory(fname, "STATIONXML")
    
    # DECIMATE DATA DOWN TO 1 Hz
    decim = Decimator([5, 5, 2])
    stream_decim = decim.decimate(stream)
    inv_decim = decim.update_inventory(inv, stream)
    
    sd_orig = SpectralDensity.from_stream(stream_decim, inv=inv_decim)
    
    # USE SIMPLE ROTATION TO REDUCE VERTICAL CHANNEL NOISE
    rot_stream = CleanRotator(stream_decim) # Downloads and stores a list of EQs
    rot_sd = SpectralDensity.from_stream(rot_stream, inv=inv_decim)
    
    # USE TRANSFER FUNCTION BASED DATA CLEANER TO REDUCE VERTICAL CHANNEL NOISE
    dc = DataCleaner(rot_stream, ['*1', '*2', '*H'])
    # first clean the stream, then calculate the spectral density
    rot_stream_dc = dc.clean_stream(rot_stream)
    sd_rot_dc = SpectralDensity.from_stream(rot_stream_dc, inv=inv_decim)
    # directly calculate the spectral density, with the datacleaner as input
    sd_rot_sddc = dc.clean_stream_to_sdf(rot_stream)
    
    # PLOT THE RESULTS
    fig, ax = plt.subplots()
    for sd, label in zip((sd_orig, sd_rot, sd_rot_dc, sd_rot_sddc),
                          ('original', 'rotated', 'rot + clean', 'rot+clean(sd)')
                        ):
        z_id = find('*Z', sd.channels)
        ax.semilogx(sd.freqs, sd.autospect(z_id), label=label)
    ax.le
    print(stream)
    stream.plot()
