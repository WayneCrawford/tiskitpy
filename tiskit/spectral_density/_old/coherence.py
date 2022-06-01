"""
Functions to calculate spectra, coherences and transfer functions
"""
import scipy.signal as ssig
from matplotlib import mlab
import matplotlib.pyplot as plt
import numpy as np
import math as m
import warnings

from .utils import choose_window_scipy, choose_window_mlab, seed_code

# Set variables
spect_library = 'scipy'  # 'mlab' or 'scipy': mlab gives weird coherences!


class Coherence:
    def __init__(self, data, chan_nums):
        self.data = data
        self.chan_nums = chan_nums


class Coherences:
    def __init__(self, freqs, cohs, num_windows, signif_level,
                 starttime, endtime, stats):
        self.freqs = np.array(freqs)
        self.cohs = cohs
        for coh in cohs:
            assert self.freqs.shape == self.cohs.shape
        self.num_window = num_windows
        self.signif_level = signif_level
        self.starttime = starttime
        self.endtime = endtime

    @classmethod
    def calc(cls, st, window_length=1000, window_type='prol1pi'):
        """
        Calculate coherences between channels of a data stream

        TO REMOVE PART COHERENT WITH ANOTHER CHANNEL, WILL HAVE TO DO
        WELCH MYSELF (MUST APPLY REMOVAL AT THE LEVEL OF EACH FFT)

        :type st: :class:`~obspy.core.stream.Stream`
        :param tr: Stream to be processed
        :type window_length: `numeric`
        :param window_length: minimum FFT window length in seconds
        :param window_type: 'fft_taper', 'prol1pi' or 'prol4pi'
        :returns: list of dictionaries containing freqs, data, units, name
        """
        cohers = []
        for i in range(len(st)-1):
            for j in range(i+1, len(st)):
                # print(i,j)
                tr_i = st[i]
                tr_j = st[j]
                tr_i.data = tr_i.data.astype(np.float64)
                tr_j.data = tr_j.data.astype(np.float64)

                # Verify that channels are compatible (same sampling rate,
                # length and starttime)
                tmp = 'Channels {:d} and {:d}'.format(i, j)
                if tr_i.stats.sampling_rate != tr_j.stats.sampling_rate:
                    warnings.warn(tmp + 'have different samp rates')
                    cohers.append([])
                    continue
                sampling_rate = tr_i.stats.sampling_rate
                if len(tr_i.data) != len(tr_j.data):
                    warnings.warn(tmp + 'have different lengths')
                    cohers.append([])
                    continue
                data_samples = len(tr_i.data)
                if abs(tr_i.stats.starttime - tr_j.stats.starttime) >\
                        1 / sampling_rate:
                    warnings.warn('tmp +  ' + 'are offset by > one sample')
                    cohers.append([])
                    continue
                starttime = tr_i.stats.starttime

                # Calculate Coherence
                nfft = 2**(m.ceil(m.log2(window_length * sampling_rate)))
                nlap = int(0.75 * nfft)
                if spect_library == 'mlab':   # MLAB GIVES STRANGE ANSWER
                    window = choose_window_mlab(window_type)
                    Cxy, _freq = mlab.cohere(tr_i.data, tr_j.data, nfft,
                                             sampling_rate,
                                             detrend=mlab.detrend_linear,
                                             window=window,
                                             noverlap=nlap, sides='onesided',
                                             scale_by_freq=False)
                # SCIPY GIVES SIMILAR ANSWER TO MATLAB
                elif spect_library == 'scipy':
                    window = choose_window_scipy(window_type, nfft)
                    _freq, Cxy = ssig.coherence(tr_i.data, tr_j.data,
                                                sampling_rate,
                                                window=window,
                                                nperseg=nfft,
                                                detrend="linear",
                                                noverlap=nlap)
                else:
                    warnings.warn('Unknown spectra library: "{}"'.format(
                        spect_library))
                    return False

                nW = 1 + np.floor((data_samples - nfft) / nlap)
                csl = np.sqrt(2. / nW)  # coherency 95% significance level

                # Make dictionary, leaving out first freq/Cxy entry (offset)
                cohers.append(Coherence(Cxy[1:], (i, j)))
        return cls(_freq[1:], cohers, nW, csl, starttime,
                   starttime + data_samples / sampling_rate,
                   [s.stats for s in st])

    def plot(self, outfile=None):
        """
        plot coherences calculated using calc_cohers (should become a class)
        """
        nCohers = len(self.cohers)
        nRows = m.ceil(m.sqrt(nCohers))
        nCols = m.ceil(nCohers / nRows)
        plt.figure(1)
        i = 0
        for coher in self.cohers:
            i += 1
            plt.subplot(nRows, nCols, i)
            plt.semilogx(self.freqs, np.absolute(coher.data))
            plt.title(seed_code(self.stats[coher.chan_nums[0]]) + '-' +
                      seed_code(self.stats[coher.chan_nums[1]]))
            plt.ylim([0, 1])
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
        return


if __name__ == "__main__":
    import doctest
    doctest.testmod()
