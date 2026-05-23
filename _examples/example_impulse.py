import numpy as np
from matplotlib import pyplot as plt

from tiskitpy import FIRConverter
from obspy.core import read, UTCDateTime, Trace, Stream

converter = FIRConverter.from_conv_file('../conv_coeffs/lc2000_fir3_250.conv.json')

n_data=140
i_impulse=n_data/2
sr = 100
impulse_offset = 70

impulse_data = np.zeros((n_data,), dtype=np.float64)
impulse_data[impulse_offset] = 1
impulse_data = np.convolve(impulse_data, converter.firzeros)
impulse_data = impulse_data[::converter.decimation_factor]    # decimate
tr = Trace(impulse_data, header={
    'sampling_rate': sr,
    'network': 'XX',
    'location': '00',
    'station': 'SSSSS',
    'channel': 'SHZ',
    'starttime': UTCDateTime('2008-01-01T00:00')-((impulse_offset + converter.timetag)/converter.decimation_factor)/sr
    })
tr.data = impulse_data
st = Stream([tr])
st_min, _ = converter.apply(st)
st_compare = st.select(component='Z') + st_min.select(component='Z')
st_compare.plot()
