from tiskitpy import FIRConverter
from obspy.core import read, UTCDateTime

# Create a converter from a zero-phase FIR file
converter = FIRConverter.from_zeros_file('data/cs5322_fir3.json', 0)

# Plot technical figures about the input FIR
print('Plotting Z-plane representation of the input filter')
converter.plot_zplane()
print('Plotting decomposition of the input filter')
converter.plot_impulse_parts()

# Apply the converter to the data
datafile = 'data/LSVEI_20150827T0612.mseed'
plotstart = UTCDateTime('2015-08-27T06:12:37.5')
plot_seconds = 0.6

# Read and convert data
st = read(datafile, 'MSEED')
st_min, _ = converter.apply(st, new_loc_code='01')

# Plot
print(f'Plotting Lucky Strike event at {str(plotstart)}.  Locations are:\n'
      '    "00": the original data\n'
      '    "01": converted to minimum phase')
st_compare = st.select(component='Z') + st_min.select(component='Z')
st_compare.resample(500)  # Oversample, to smooth the waveform
st_compare.plot(starttime=plotstart,
                endtime=plotstart + plot_seconds)
