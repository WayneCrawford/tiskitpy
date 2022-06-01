"""
Calculate spectra and coherences for a given station/time period
"""
from crawtools.spectral import SpectralDensity, PSD
from matplotlib import pyplot as plt
import numpy as np

from read_data_inventory import read_data_inventory


windowtype='prol1pi'
stream, inv = read_data_inventory()

# Calculate spectra
spect = SpectralDensity.from_stream(stream, inv=inv, windowtype=windowtype)

# Show different types of plots
spect.plot()
spect.plot(overlay=True)
spect.plot_cross_spectra()
spect.plot_cross_spectra(show_coherence=True)
spect.plot_coherences()
spect.plot_coherences(overlay=True)
