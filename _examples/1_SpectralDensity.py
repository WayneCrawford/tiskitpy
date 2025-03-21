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
file_base = "1_SpectralDensity"
sd.plot()  #outfile=f'{file_base}_plot.png')

# plot results, overlaid
sd.plot(overlay=True)  # , outfile=f'{file_base}_plot_overlay.png')

# plot coherences
sd.plot_coherences(display="full")  #  , outfile=f'{file_base}_coher_full.png')

# plot coherences, overlaid
sd.plot_coherences(display="overlay")  # , outfile=f'{base}_coher_overlay.png')

# plot coherences, sparse
sd.plot_coherences(display="sparse")  # , outfile=f'{file_base}_coher_sparse.png')

# plot coherences, taking up minimal space
sd.plot_coherences(display="minimal")  # , outfile=f'{file_base}_coher_minimal.png')
