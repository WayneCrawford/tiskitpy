from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskitpy import DataCleaner, SpectralDensity, CleanSequence

# Read data and inventory
stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')

# Calculate a Datacleaner that will subtract *1, then *2, then *H
print('='*80)
dc = DataCleaner(stream, ['*1','*2','*H'])

# Clean the data, then construct a stream with original and cleaned channels
print('='*80)
stream_cleaned = dc.clean_stream(stream)
z_compare = stream.select(channel='*Z') + stream_cleaned.select(channel='*Z')

# If you print and plot the stream "normally", can not see which is which
print('='*80)
print(z_compare)
z_compare.plot(outfile='3_DataCleaner_tagged_timeseries.png')

# The CleanSequence stream_plot() and stream_print() methods wrap the
# Stream functions after adding clean sequence info to the seed_id
CleanSequence.stream_print(z_compare)
CleanSequence.stream_plot(z_compare, outfile='3_DataCleaner_tagged_timeseries_csplot.png')

# compare spectral densities
# (tiskitpy plot() automatically include CleanSequence information)
print('='*80)
sd_compare = SpectralDensity.from_stream(z_compare, inv=inv)
sd_compare.plot(overlay=True, outfile='3_DataCleaner_sd_overlay.png')
