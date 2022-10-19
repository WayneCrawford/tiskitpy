from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskit import SpectralDensity

stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
sd = SpectralDensity.from_stream(stream, inv=inv)
sd.plot()
sd.plot(overlay=True)
print(sd)