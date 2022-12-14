from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskit import DataCleaner, SpectralDensity

stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
dc = DataCleaner(stream, ['*1','*2','*H'])
stream_cleaned = dc.clean_stream(stream)
z_compare = stream.select(channel='*Z') + stream_cleaned.select(channel='*Z')
z_compare.plot()
print(z_compare)
sd_compare = SpectralDensity.from_stream(z_compare, inv=inv)
sd_compare.plot(overlay=True)
