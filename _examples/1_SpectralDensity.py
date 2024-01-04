from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskitpy import SpectralDensity

# read data and inventory
stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')

# Calculate Spectral Density
sd = SpectralDensity.from_stream(stream, inv=inv)
print(sd)

# plot results
base = "1_SpectralDensity"
sd.plot(outfile=f'{base}_plot.png')
sd.plot(overlay=True, outfile=f'{base}_plot_overlay.png')
sd.plot_coherences(display="full", outfile=f'{base}_coher_full.png')
sd.plot_coherences(display="overlay", outfile=f'{base}_coher_overlay.png')
sd.plot_coherences(display="sparse", outfile=f'{base}_coher_sparse.png')
sd.plot_coherences(display="minimal", outfile=f'{base}_coher_minimal.png')
