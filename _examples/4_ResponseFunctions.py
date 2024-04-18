from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskitpy import SpectralDensity, ResponseFunctions

# Read the data and metadata, and calculate a SpectralDensity
stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
sd = SpectralDensity.from_stream(stream, inv=inv)

# Calculate the response functions with respect to the 'LDH' channel
rfs = ResponseFunctions(sd, '*H')

# Print the ResponseFunctions
print(rfs)

# Plot the ResponseFunctions
rfs.plot()  # outfile='4_ResponseFunctions.png')
