from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime

from tiskitpy import Decimator

# Specify FDSN server dataset
client_address = 'IRIS'
net = 'Z5'
sta = 'BB650'  # Trillium T240. 100 sps
loc = '*'
cha = 'HH*,HDH'
starttime = UTCDateTime('2014-01-01T00:52')
endtime = UTCDateTime('2014-01-01T00:55')

# Read data and inventory
client = Client(client_address)
print('reading data from FDSN server...')
stream = client.get_waveforms(network=net, station=sta, channel=cha,
                              location=loc, starttime=starttime,
                              endtime=endtime)
print('reading inventory from FDSN server...')
inv = client.get_stations(starttime=starttime, endtime=endtime,
                          network=net, station=sta, channel='*', location='*',
                          level='response')
print(inv)

# DECIMATE DATA by 100 and update inventory
decim = Decimator([5, 5, 4])
stream_decim = decim.decimate(stream)

# Update the inventory to include the decimated channels
inv_decim = decim.update_inventory(inv, stream)

# Inventory now has 4 more channels (LDH, LHZ, LH1 and LH2)
print(inv_decim)

# Combine the original and decimated Z channels into one stream
compare_z = stream.select(channel='*Z') + stream_decim.select(channel='*Z')
print(compare_z)

# Plot the original and decimated Z channels
compare_z.plot()  # show=True, outfile='2_Decimator_time_series.png')
