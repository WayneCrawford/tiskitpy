from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime

from tiskit import Decimator

# Read data and inventory from IRIS server (takes a while)client_address = 'IRIS'
net = 'Z5'
sta = 'BB650'  # Trillium T240. 100 sps
loc = '*'
cha = 'HH*,HDH'
starttime = UTCDateTime('2014-01-01T00:52')
endtime = UTCDateTime('2014-01-01T00:55')
client = Client(client_address)
print('reading data from FDSN server...')
stream = client.get_waveforms(network=net, station=sta, channel=cha,
                              location=loc, starttime=starttime,
                              endtime=endtime)
print('reading inventory from FDSN server...')
inv = client.get_stations(starttime=starttime, endtime=endtime,
                          network=net, station=sta, channel='*', location='*',
                          level='response')

# DECIMATE DATA 100x
decim = Decimator([5, 5, 4])
stream_decim = decim.decimate(stream)
inv_decim = decim.update_inventory(inv, stream)
print(inv)
print(inv_decim)

compare_z = stream.select(channel='*Z') + stream_decim.select(channel='*Z')
print(compare_z)
compare_z.plot()
