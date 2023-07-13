"""
Calculate spectra and coherences for a given station/time period
"""
from matplotlib import pyplot as plt

from crawtools.spectral import SpectralDensity, TransferFunctions

from read_data_inventory import read_data_inventory


stream, inv = read_data_inventory()

# Calculate spectra
spect = SpectralDensity.from_stream(stream, inv=inv)
channels = spect.channels
# print(spect.channels)

# Calculate transfer functions w.r.t. pressure
tf = TransferFunctions(spect, '*BDH')
# Plot transfer functions
tf.plot()
# print(tf.xfs)

# Compare transfer functions calculated assuming noise on the input versus output channel
# Calculate transfer functions, assuming noise on the input channel
inp = '*BDH'
outp = '*BHZ'
tf_ni = TransferFunctions(spect, inp, noise_chan="input")
fig, ax = plt.subplots(1,1)
ax.loglog(tf.freqs, abs(tf.values(outp, zero_as_none=True)), color='red',
          label='assume noise on output channel')
ax.loglog(tf_ni.freqs, abs(tf_ni.values(outp, zero_as_none=True)), color='blue',
          label='assume noise on input channel')
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel(f'{outp}/{inp}')
ax.set_title('Transfer Functions')
ax.legend()
plt.show()

# help(TransferFunction)

# Print out data and coordinates
chan = '*BHZ'
print(f'{tf._ds=}')
print('OUTPUT CHANNEL INDEPENDENT')
print(f"{tf.freqs=}")
print(f"{tf.input_channel=}")
print(f"{tf.output_channels=}")
print(f"{tf.input_units=}")
print(f"{tf.n_windows=}")
print(f"{tf.coh_signif(0.95)=:.2f}")
print(f"")
print('OUTPUT CHANNEL DEPENDENT')
print(f"{chan=}")
print(f"{tf.values(chan)=}")
print(f"{tf.uncert(chan)=}")
print(f"{tf.output_units(chan)=}")
print(f"{tf.noise_channel(chan)=}")




# Compare original and cleaned vertical channels