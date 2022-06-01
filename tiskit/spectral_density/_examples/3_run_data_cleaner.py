"""
Calculate spectra and coherences for a given station/time period
"""
import sys

from matplotlib import pyplot as plt
import numpy as np
import fnmatch

from crawtools.spectral import SpectralDensity, DataCleaner
from read_data_inventory import read_data_inventory

datasets = {'PacificArray_CC06': dict(network='XE', station='CC06', location='*',
                                      channel='B*',
                                      start_time = '2018-05-25T00:00:00',
                                      end_time = '2018-05-26T00:00:00',
                                      base_url='IRIS'),
            'RHUM-RUM_RR28': dict(network='YV', station='RR28', location='*',
                                  channel='B*',
                                      start_time = '2012-12-01T00:00:00',
                                      end_time = '2012-12-02T00:00:00',
                                      base_url='RESIF'),
            'RHUM-RUM_RR29': dict(network='YV', station='RR29', location='*',
                                  channel='B*',
                                      start_time = '2012-12-01T00:00:00',
                                      end_time = '2012-12-02T00:00:00',
                                      base_url='RESIF'),
            'AlpArray_A422A': dict(network='Z3', station='A422A', location='*',
                                  channel='B*',
                                      start_time = '2017-07-01T00:00:00',
                                      end_time = '2017-07-02T00:00:00',
                                      base_url='RESIF')
          }
for dataset in datasets.keys():
    print('\nReading data and metadata')
    stream, inv = read_data_inventory(**datasets[dataset])

    counter=0
    for clean_seq in (['*BDH'], ['*BH1']): # , ['*BDH', '*BH1', '*BH2']):
        # Clean vertical channel
        print('\nCalculating DataCleaner')
        dc = DataCleaner(stream, clean_seq, show_plots=False,
                         max_freq=0.1)

        orig_sdf = SpectralDensity.from_stream(stream)

        print('\nRunning DataCleaner.clean_stream(time)')
        clean_stream_time = dc.clean_stream(stream, in_time_domain=True)
        clean_stream_time_sdf = SpectralDensity.from_stream(clean_stream_time)

        print('\nRunning DataCleaner.clean_stream(freq)')
        clean_stream_freq = dc.clean_stream(stream, in_time_domain=False)
        clean_stream_freq_sdf = SpectralDensity.from_stream(clean_stream_freq)

        print('\nRunning DataCleaner.clean_stream_to_sdf')
        clean_sdf = dc.clean_stream_to_sdf(stream)

        sdf = SpectralDensity.from_stream(stream)
        print('\nRunning DataCleaner.clean_sdf')
        clean_sdf2 = dc.clean_sdf(sdf)
    
        fig, ax = plt.subplots(1)
        for dat, label in zip((orig_sdf, clean_stream_time_sdf, clean_stream_freq_sdf, clean_sdf, clean_sdf2),
                              ('original', 'clean_stream(time)', 'clean_stream(freq)', 'clean_stream_to_sdf', 'clean_sdf')):
            chan = fnmatch.filter(dat.channels, '*BHZ*')[0]
            print(f'{chan=}')
            ax.semilogx(dat.freqs, 10*np.log10(dat.autospect(chan)), label=label)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Raw PSD (counts^2/Hz)')
        ax.set_title(f'{dataset}: Comparison of data_cleaners, clean_seq={clean_seq}')
        ax.legend()
        plt.savefig(f'{dataset}_{counter}.png')
        counter += 1