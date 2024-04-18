from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskitpy import CleanRotator, SpectralDensity

# Read the data and metadata, and calculate a CleanRotator
stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
cr = CleanRotator(stream)

# Apply the CleanRotator to the data
stream_rotated = cr.apply(stream)
    
# Combine the rotated and unrotated data, for comparison
z_compare = stream_rotated.select(channel='*Z') + stream.select(channel='*Z')

# Print and plot the combined data
print(z_compare)
z_compare.plot()  # outfile='5_CleanRotator_z_compare.png')

# Create and plot a SpectralDensity object for the rotated and unrotated data
sd_compare = SpectralDensity.from_stream(z_compare, inv=inv)
sd_compare.plot(overlay=True)  # , outfile='5_CleanRotator_spect_compare.png')
