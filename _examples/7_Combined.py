import fnmatch

import numpy as np
from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from matplotlib import pyplot as plt

from tiskitpy import CleanRotator, SpectralDensity, DataCleaner

stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
sd_orig = SpectralDensity.from_stream(stream, inv=inv)

# USE SIMPLE ROTATION TO REDUCE VERTICAL CHANNEL NOISE
rotator = CleanRotator(stream)
rot_stream = rotator.apply(stream)
sd_rot = SpectralDensity.from_stream(rot_stream, inv=inv)

# USE TRANSFER FUNCTION BASED DATA CLEANER TO REDUCE VERTICAL CHANNEL NOISE
dc = DataCleaner(stream, ['*1', '*2', '*H'])
dc_rot = DataCleaner(rot_stream, ['*1', '*2', '*H'])

# Clean the stream, then calculate the spectral density
stream_dc = dc.clean_stream(stream)
rot_stream_dc = dc_rot.clean_stream(rot_stream)
# For now, must remove '-ROT' from location code
rot_stream_dc = CleanSequence.remove_from_loc(rot_stream_dc)
sd_dc = SpectralDensity.from_stream(stream_dc, inv=inv)
sd_rot_dc = SpectralDensity.from_stream(rot_stream_dc, inv=inv)

# Directly calculate the spectral density, with the DataCleaner as input
sd_rot_sddc = dc_rot.clean_stream_to_sdf(rot_stream, inv=inv)

# PLOT THE RESULTS
fig, ax = plt.subplots()
for sd, label in zip(
        (sd_orig, sd_rot, sd_dc, sd_rot_dc, sd_rot_sddc),
        ('original', 'rotated', 'cleaned', 'rot + clean', 'rot+clean(sd)')):
    # print(sd.channel_names)
    z_id = fnmatch.filter(sd.channel_names, '*.LHZ*')[0]
    ax.semilogx(sd.freqs, 10*np.log10(sd.autospect(z_id)), label=label)
ax.legend()
plt.show()
plot.outfile(f'7_Combined.png')
