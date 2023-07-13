"""
Functions to calculate spectra, coherences and transfer functions
"""
import math as m
import warnings

from obspy.core.trace import Trace
from obspy.core.inventory import Inventory
import scipy.signal as ssig
from matplotlib import mlab
import matplotlib.pyplot as plt
import numpy as np

from .Peterson_noise_model import PetersonNoiseModel
from .utils import choose_window_scipy, choose_window_mlab, seed_code

# Set variables
spect_library = 'scipy'  # 'mlab' or 'scipy': mlab gives weird coherences!
use_PSDs = True


class PSD:
    """
    Power Spectral Density class
    """
    def __init__(self, freqs=None, data=None, units=None, seed_id=None):
        """
        Arguments:
            freqs (list or nparray): 1-D array of frequencies
            data (nparray): 1-D array of PSD values
            units (str): data units
            seed_id (str): SEED ID of the data
        """
        self.freqs = np.array(freqs)
        self.data = np.array(data)
        if data is not None:
            assert freqs.shape == data.shape
        self.units = units
        self.seed_id = seed_id

    def __repr__(self):
        """
        String describing object
        """
        s = f'PSD(freqs, data, units={self.units}, seed_id={self.seed_id}) '
        s += '<data.shape={}, freqs={:g}-{:g}>'.format(
            self.data.shape, self.freqs[0], self.freqs[-1])
        return s

    @classmethod
    def calc(cls, tr, window_length:float=1000., window_type='prol1pi', inv=None):
        """
        Calculate PSD of a data trace
        Based on obspy PPSD function

        REMOVING PART COHERENT WITH ANOTHER CHANNEL WILL REQUIRE DOING
        WELCH MYSELF (MUST APPLY REMOVAL AT THE LEVEL OF EACH FFT)

        Arguments:
            tr (obspy.core.stream.Trace): waveform, should have response attached
            window_length (float): minimum FFT window length in seconds
            window_type (str): 'fft_taper', 'prol1pi' or 'prol4pi'
            inv (obspy.core.inventory.Inventory): station inventory containing
                instrument response.  If 'None', response must be attached to
                the trace.
        """
        if not isinstance(tr, Trace):
            raise TypeError('trace is not an obspy Trace')
        if inv is not None and not isinstance(inv, Inventory):
            raise TypeError('inv is not an obspy Inventory')

        sampling_rate = tr.stats.sampling_rate
        data_len = tr.stats.endtime - tr.stats.starttime
        if window_length > data_len/2:
            window_length = int(data_len/3)
            print('data is only {:g}s long, reducing window_length to {:g}s'
                  .format(data_len, window_length))
        nfft = 2**(m.ceil(m.log2(window_length * sampling_rate)))
        nlap = int(0.75 * nfft)

        if spect_library == 'mlab':
            window = choose_window_mlab(window_type)
            spec, _freq = mlab.psd(tr.data, nfft, sampling_rate,
                                   detrend=mlab.detrend_linear,
                                   window=window, noverlap=nlap,
                                   sides='onesided', scale_by_freq=True)
        elif spect_library == 'scipy':
            window = choose_window_scipy(window_type, nfft)
            _freq, spec = ssig.welch(tr.data, sampling_rate, nperseg=nfft,
                                     window=window, detrend="linear",
                                     noverlap=nlap)
        else:
            warnings.warn('Unknown spectra library: "{}"'.format(
                spect_library))
            return False

        # leave out first entry (zero-freq)
        spec = spec[1:]
        _freq = _freq[1:]

        # Remove the response using the same conventions
        # Since the power is squared, square the sensitivity
        resp = None
        if inv is not None:
            try: 
                resp = inv.get_response(tr.id, tr.stats.starttime)
            except Exception as e:
                print(f'No response for {tr.id} found in inventory')
                print(e)
        if resp is None:
            try:
                resp = tr.stats.response
            except Exception:
                print('No response found in trace')

        if resp is not None:
            evalresp, freqs = resp.get_evalresp_response(
                t_samp=1 / sampling_rate, nfft=nfft, output="VEL")
            evalresp = evalresp[1:]
            freqs = freqs[1:]
            if not np.all(freqs == _freq):
                warnings.warn('psd freqs not same as evalresp freqs')
            # Get the amplitude response (squared)
            respamp = np.absolute(evalresp * np.conjugate(evalresp))
            # Make omega with the same conventions as spec
            w = 2.0 * m.pi * _freq
            # Remove response
            if resp.response_stages[0].input_units.upper()[:2] == "PA":
                print('Channel {} has input_units "{}": treating as hydrophone'.
                      format(tr.stats.channel,
                             resp.response_stages[0].input_units))
                spec = spec / respamp
                PSD_units = "dB ref 1 Pa^2/Hz"
            else:
                spec = (w**2) * spec / respamp
                PSD_units = "dB ref 1 (m/s^2)^2/Hz"
        else:
            PSD_units = "dB ref 1 count^2/Hz"
        return cls(_freq, spec, PSD_units, seed_code(tr.stats))

    def plot(self, ax=None, show=True, outfile=None, show_Peterson=True):
        """
        Plot a PSD
        """
        if ax:
            pass
            # plt.gca = ax
        else:
            plt.figure()
            ax = plt.gca()
        ax.semilogx(self.freqs, 10 * np.log10(self.data))
        if show_Peterson and not self.seed_id[-1] == 'H':
            lownoise, highnoise = PetersonNoiseModel(self.freqs, True)
            ax.semilogx(self.freqs, lownoise, '--')
            ax.semilogx(self.freqs, highnoise, '--')
        ax.set_ylabel(self.units)
        # ax.title()
        ax.set_title(self.seed_id)
        if outfile:
            plt.savefig(outfile)
        elif show:
            plt.show()


class PSDs:
    """
    List of PSD
    """
    def __init__(self, PSDs=None):
        assert isinstance(PSDs, list)
        for p in PSDs:
            assert isinstance(p, PSD)
        self.PSDs = PSDs

    def __str__(self):
        s = ''
        for PSD in self.PSDs:
            s += PSD.__str__() + '\n'

    def __len__(self):
        return len(self.PSDs)

    @classmethod
    def calc(cls, st, window_length=1000):
        """
        Calculate PSDs of a data stream

        :type st: :class:`~obspy.core.stream.Stream`
        :param st: Stream to be processed
        :type window_length: `numeric`
        :param window_length: minimum FFT window length in seconds
        :returns: list of dictionaries containing freqs, data, units, name
        """
        PSDs = []
        for tr in st:
            PSDs.append(PSD.calc(tr, window_length))
        return cls(PSDs=PSDs)

    def plot(self, outfile=None):
        """
        plot PSDs
        """
        nPSDs = len(self.PSDs)
        nRows = m.floor(m.sqrt(nPSDs))
        nCols = m.ceil(nPSDs/nRows)
        plt.figure(1)
        i = 0
        for p in self.PSDs:
            i += 1
            ax = plt.subplot(nRows, nCols, i)
            p.plot(ax=ax, show=False)
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
        return


if __name__ == "__main__":
    import doctest
    doctest.testmod()
