from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskitpy import DataCleaner, SpectralDensity, CleanSequence as CS

# Read data and inventory
stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')

# Calculate and apply Datacleaner from/to stream, subtracting *1, then *2, then *H
print('='*80)
dc = DataCleaner(stream, ['*1','*2','*H'])
print('='*80)
stream_cleaned = dc.clean_stream(stream)
print('='*80)

# Construct and plot a stream with original and cleaned Z channels
z_compare = stream.select(channel='*Z') + stream_cleaned.select(channel='*Z')
print('='*80)
print(z_compare)

# Tag the seed_id with modifs so that it shows up in obspy plots
z_compare_tagged = CS.seedid_tag(z_compare)
print('='*80)
print(z_compare_tagged)
z_compare_tagged.plot(outfile='3_DataCleaner_tagged_timeseries.png')

# compare spectral densities
print('='*80)
sd_compare = SpectralDensity.from_stream(z_compare, inv=inv)
sd_compare.plot(overlay=True, outfile='3_DataCleaner_sd_overlay.png')
