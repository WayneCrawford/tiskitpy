from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskit import CleanRotator, SpectralDensity

stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
cr = CleanRotator(stream)
stream_rotated = cr.apply(stream)

# Change channel names so as not to conflict with existing
for tr in stream_rotated:
    tr.stats.channel = tr.stats.channel[0] + 'X' + tr.stats.channel[2]
inv_x = inv.copy()
for ch in inv_x[0][0]:
    ch.code = ch.code[0] + 'X' + ch.code[2]
    inv[0][0].channels.append(ch)
    
z_compare = stream.select(channel='*Z') + stream_rotated.select(channel='*Z')
z_compare.plot()
print(z_compare)
sd_compare = SpectralDensity.from_stream(z_compare, inv=inv)
sd_compare.plot(overlay=True)
