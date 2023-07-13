"""
Calculate spectra and coherences for a given station/time period
"""
from pathlib import Path

from obspy.clients.fdsn.client import Client
from obspy.core import UTCDateTime
from obspy.core.stream import read
from obspy.core.inventory import read_inventory

#### EXAMPLE EXPERIMENTS ####
# - Pacific Array (2018XE) 
#     - SIO BBOBS stations CC01-12, EC01-05, EE01-04, WC01-05, WW01-04
#     - Recording period: 2018-04-29 to 2019-05-06, though some ended early
#     - Quiet periods: 2018/05/10-17, 5/19-6/20, 6/22-7/05,
#                       9/19-9/27, 2019/02/03-02/11


def read_data_inventory(network='XE', station='CC06', location='*',
                        channel='B*',
                        start_time = '2018-05-25T00:00:00',
                        end_time = '2018-05-26T00:00:00',
                        base_url='IRIS'):
    """
    Reads in data and an inventory for the spectral example codes
    
    Args:
        network (str): network to request
        station (str): station to request
        location (str): location to request
        channel (str): channels to request ('*' is a wildcard)
        start_time (str): start time for data
        end_time (str): end time of data to retrieve (YYYY-MM-DDTHH:MM:SS)
        base_url (str): FDSN data center or URL to interrogate
    """
    # Read the data
    t1, t2 = UTCDateTime(start_time), UTCDateTime(end_time)
    fbase = '{}.{}.{}.{}.{}-{}'.format(network.replace("*","-"),
                                       station.replace("*","-"),
                                       location.replace("*","-"),
                                       channel.replace("*","-"),
                                       t1.format_fissures().split(".")[0],
                                       t2.format_fissures().split(".")[0])
    if not Path(fbase+'.mseed').is_file() or not Path(fbase+'.xml').is_file():
        client = Client(base_url)
    if not Path(fbase+'.mseed').is_file():
        print('Reading waveforms from FDSN server...', end='', flush=True)
        stream = client.get_waveforms(network, station, location, channel, t1, t2)
        stream.merge(fill_value='latest')
        stream.write(fbase +'.mseed', format='MSEED')
    else:
        print(f'Reading waveforms from file {fbase +".mseed"} ...', end='', flush=True)
        stream = read(fbase+'.mseed', format='MSEED')
    print('Done reading')

    # Plot waveforms
    # stream.plot()

    # Read the inventory
    if not Path(fbase+'.xml').is_file():
        print('Reading inventory from FDSN server...', end='', flush=True)
        inv = client.get_stations(network=network, station=station, starttime=t1,
                                  endtime=t2, level='response')
        inv.write(fbase+'.xml', format='STATIONXML')
    else:
        print('Reading inventory from file...', end='', flush=True)
        inv = read_inventory(fbase+'.xml', format='STATIONXML')
    print('Done reading')
    
    return stream, inv
