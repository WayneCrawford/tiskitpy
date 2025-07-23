from copy import deepcopy

from scipy.fft import irfft
import scipy.signal.windows as sp_windows
import numpy as np
from matplotlib import pyplot as plt
from obspy.core import Trace

# from ..spectral_density import SpectralDensity
from obspy.core.inventory.response import Response

class PSDVals():
    """
    Holds PSD frequencies and values
    
    Attributes:
        freqs (list): frequencies
        values (list): values in dB
        value_units (str): units of the values
    """
    def __init__(self, freqs_and_vals, value_units="unknown"):
        """
        Args:
            freqs_and_vals (tuple)): frequencies, PSD values, and is_dB (bool), entered as:
                ([[freq1, value1],
                  [freq2, value1],
                  ...
                  [freqN, valueN]],
                 is_dB)
                frequencies must be monotonically increasing
                is_DB: True: Input values are in dB ref 1 {value_units}^2/Hz
                       False: Input values are in {value_units}
            value_units (str): Units of values before converting to dB
        """
        freq_val_list = freqs_and_vals[0]
        is_dB = freqs_and_vals[1]
        if isinstance(freq_val_list, (list, tuple)):
            freq_val_list = np.array(freq_val_list)
        else:
            assert isinstance(freq_val_list, np.ndarray)
        assert np.all(np.diff(freq_val_list[:, 0]) > 0), 'freqs are not monotonically increasing'
        assert isinstance(value_units, str)
        self.freqs = freq_val_list[:, 0]
        self.values = freq_val_list[:, 1]
        if is_dB is not True:
            self.values = 20 * np.log10(self.values)
        if value_units.isalpha():
            # Purely alphanetical value_units don't need parentheses
            self.value_units = f'dB ref 1 {value_units}^2/Hz'
        elif value_units[0] == '(' and value_units[-1] == ')':
            # There are already parentheses around value_units, no need for more
            self.value_units = f'dB ref 1 {value_units}^2/Hz'
        else:
            # Add parenthesis around value_units
            self.value_units = f'dB ref 1 ({value_units})^2/Hz'

    @classmethod    
    def from_loglog_slope(cls, val_1Hz, slope, logf_low, logf_high, logf_step,
                          value_units="unknown"):
        """Create PSDVals object with a loglog slope
        
        Args:
            val_1Hz (float): PSD level (dBs)  at 1 Hz
            slope (float): PSD slope in log(dBs)/log(freq)
            logf_low (float): log10 of minimum frequency
            logf_high (float): log10 of maximum frequency
            logf_step (float): log10 frequency step
            value_units (str): Units of values before converting to dB
        """
        x = PSDVals.sloped_freqs_and_values(val_1Hz, slope, logf_low, logf_high,
                                            logf_step)
        return cls((x, True), value_units)

    def __str__(self):
        s = f"<PSDVals>:\n"
        s += f"        units={self.value_units}\n"
        s += "             freq  |   value  \n"
        s += "        ---------- | ---------\n"
        for f, v in zip(self.freqs, self.values):
            s += f"        {f:<10.4g} | {v:9.4g}\n"
        return s

    def __add__(self, other):
        """Add a scalar to each PSD value"""
        output = self.copy()
        output.values = [x + other for x in self.values]
        return output
        
    def __sub__(self, other):
        """Subtract a scalar from each PSD value"""
        output = self.copy()
        output.values = [x - other for x in self.values]
        return output
        
    def copy(self):
        """Return a deep copy of self"""
        return deepcopy(self)
        
    def resample(self, new_f, reference='loglog'):
        """
        Resample data using interpolation
        
        Args:
            new_f (list or :class:`np.array`): frequencies at which to evaluate
            reference (str): What x and y dimensions to use as reference:
                'loglog': log frequencies, log values (default)
                'semilogx': log frequencies, linear values
                'semilogy': linear frequencies, log values
                'linear': linear frequencies and values
        """
        self.values = self.resample_values(new_f, reference)
        self.freqs = new_f

    def plot(self):
        f, a = plt.subplots()
        a.semilogx(self.freqs, self.values)
        a.set_ylabel('Amplitude (dB)')
        a.set_title(self.value_units)
        plt.show()

    @staticmethod
    def _random_phases(n_values):
        rng = np.random.default_rng()
        return 360. * rng.random(n_values)

    def resample_values(self, freqs, reference='loglog'):
        """
        Resample values at the given frequencies

        Args:
            freqs (list or :class:`np.array`): frequencies at which to evaluate
            reference (str): What x and y dimensions to use as reference:
                'loglog': log frequencies, log values (default)
                'semilogx': log frequencies, linear values
                'semilogy': linear frequencies, log values
                'linear': linear frequencies and values
        """
        if reference=='linear':
            return np.log10(np.interp(freqs, self.freqs, np.pow(self.values, 10)))
        elif reference=='semilogx':
            return np.log10(np.interp(np.log10(freqs), np.log10(self.freqs), np.pow(self.values,10.)))
        elif reference=='semilogy':
            return np.interp(freqs, self.freqs, self.values)
        elif reference=='loglog':
            return np.interp(np.log10(freqs), np.log10(self.freqs), self.values)
        else:
            raise ValueError(f'{reference=} not in ("loglog", "semilogx", "semilogy", "linear")')

    @property
    def accel_as_vel(self):
        """
        Convert a PSD that was ref:acceleration to ref:velocity
        """
        if not self.value_units == 'dB ref 1 (m/s^2)^2/Hz':
            raise ValueError(f"{self.value_units=} are not '{ref}'")
        psd_list = [[f, v-20*np.log10(2*np.pi*f)]
                    for f, v in zip(self.freqs, self.values)]
        return PSDVals((psd_list, True), 'm/s')
        
    def _as_fft(self, freqs, left='taper', right='taper', phases=None,
               plotit=False):
        """
        Return an fft "equivalent" to the given Power Spectral Density
        Args:
            freqs (list, np.array or None): Resample at the given freqs
            left (None, float or 'taper'): how to handle values below the
                lowest self.freq:
                    - None: use np.interp() default (value at self.freqs[0])
                    - float: set to the given value
                    - 'taper': taper using Kaiser function
            right (None, float, or 'taper'): how to handle values above the
                highest self.freq
            phases (np.array or None): force phases to be the given values
                (must be same length as frequencies)
        """
        # VALIDATE INPUT PARAMETERS
        if not freqs[0] == 0:
            raise ValueError("Cannot create an fft without f[0] == 0")
        fdiffs = np.diff(freqs)
        if not np.all(np.abs(fdiffs - fdiffs[0]) < fdiffs[0] / 1e6):
            raise ValueError("freqs are not evenly spaced")

        # CREATE FFT FROM PSD
        # Using log(freqs) avoids bumps for widely spaced self.freqs
        fft = np.power(10., np.interp(np.log(freqs), np.log(self.freqs), self.values) / 20)
        np.nan_to_num(fft, copy=False)
        # Handle frequencies below/above min/maximum PSD frequency
        if left == 'taper':
            fft = self._add_left_taper(fft, freqs, self.freqs[0])
        elif left is not None:
            fft[freqs < self.freqs[0]] = left
        if right == 'taper':
            fft = self._add_right_taper(fft, freqs, self.freqs[-1])
        elif right is not None:
            fft[freqs > self.freqs[-1]] = right
        fft[0] = 0.  # DC = 0.

        if plotit is True:
            fig, ax = plt.subplots()
            ax.loglog(self.freqs, np.power(10., self.values/20), '+', freqs, fft)
            # ax.semilogx(self.freqs, self.values, '+', freqs, 20*np.log10(fft))
            plt.suptitle('PSDVals._as_fft()')
            plt.show()
        fft[0] = 0.  # Make sure the zero-frequency value is zero
        # Scale for sample rate and window length
        sampling_rate = 1 / (2 * freqs[-1])
        # Bendata&Piersol 1986 eqs 11.100 & 11.102
        mul_factor = np.sqrt(len(freqs) * sampling_rate / 2)
        fft *= mul_factor
        
        # ADD PHASES
        if phases is None:
            phases = np.radians(self._random_phases(len(fft)))
        return fft * np.exp(1j * phases)

    def as_trace(self, ref_trace, network=None, station=None, location=None,
                 channel=None, plotit=False, **kwargs):
        """
        Return a Trace with the given spectral shape
        
        Args:
            ref_trace (:class:`obspy.core.stream.Trace` or dict): trace whose
                parameters will be used, or dict with keys 'sampling_rate' (float),
                'starttime' (:class:`obspy.UTCDateTime`), 'endtime'
                (:class:`obspy.UTCDateTime`) and possibly 'response'
                (:class:`obspy.core.inventory.Response`).
            network (str): network code (default: value in ref Trace, or 'XX')
            station (str): station code (default: value in ref Trace, or 'STA')
            location (str): location code (default: value in ref Trace, or '00')
            channel (str): channel code (default: value in ref Trace, or 'CCC')
            plotit (bool or str): plot the components of the transformation
                (psd, fft, fresp).  If a string, use as the plot's title.
            **kwargs (dict): keyword arguments to pass to self._as_fft()
        
        Returns:
            tuple:
                trace (class:`obspy.stream.Trace`)
                phases (:class:`numpy.ndarray`): FFT phases (radians)
                    used to create this trace
                
        """
        if isinstance(ref_trace, dict):
            trace = self._make_trace(ref_trace)
        else:
            assert isinstance(ref_trace, Trace), 'ref_trace is not a Trace'
            trace = ref_trace.copy()
        sr = trace.stats.sampling_rate
        trace_pts = trace.stats.npts
        if network is not None:
            trace.stats.network = network
        if station is not None:
            trace.stats.station = station
        if location is not None:
            trace.stats.location = location
        if channel is not None:
            trace.stats.channel = channel
        npts = 2**int(np.ceil(np.log2(trace_pts)))
        f = np.linspace(0, trace.stats.sampling_rate / 2, npts)
        if 'response' in trace.stats.__dict__:
            fresp = trace.stats.response.get_evalresp_response_for_frequencies(f, output='DEF')
        else:
            fresp = np.ones(f.shape)
        if isinstance(plotit, str):
            title_text = plotit
            plotit = True
        else:
            title_text = None
        fft = self._as_fft(f, **kwargs)
        if plotit is True:
            fig, ax = plt.subplots()
            ax.loglog(f, np.abs(fft), 'b-', label='fft')
            ax.loglog(f, np.abs(fresp), 'g--', label='fresp')
            ax.loglog(f, np.abs(self._as_fft(f, **kwargs)*fresp), 'c--', label='fft*fresp')
            ax.loglog(self.freqs, np.power(10., self.values/20), 'r+', label='psd')
            ax.loglog(f, np.abs(fft/np.sqrt(len(f) * sr / 2)), 'b:', label='fft/sqrt(nf*sr/2)')
            if title_text is not None:
                ax.set_title(title_text)
            ax.legend()
            plt.show()
        trace.data = irfft(fft*fresp)[:trace_pts]
        return trace, np.angle(fft)

    @staticmethod    
    def sloped_freqs_and_values(val_1Hz, slope, logf_low, logf_high, logf_step):
        """Create freqs_and_values input for a loglog slope
        
        Args:
            val_1Hz (float): PSD level (dBs)  at 1 Hz
            slope (float): PSD slope in log(dBs)/log(freq)
            logf_low (float): log10 of minimum frequency
            logf_high (float): log10 of maximum frequency
            logf_step (float): log10 frequency step
        """
        x = [[np.power(10., logf), val_1Hz + logf * slope]
            for logf in np.arange(logf_low, logf_high, logf_step)]
        return (x, True)

    @staticmethod    
    def _make_trace(ref_trace, dtype=np.float64):
        """
        dtype (np.unit8 to np.float128) does not seem to affect output
        I'm guessing it's replaced downstream
        """
        for x in ("sampling_rate", "starttime", "endtime"):
            assert x in ref_trace, f'key "{x}" missing from ref_trace dict'
        sr = ref_trace["sampling_rate"]
        st = ref_trace["starttime"]
        trace_pts = 1 + int((ref_trace["endtime"] - st)/sr)
        stats = {"sampling_rate": sr, "starttime": st, "network": 'XX',
                 "station": 'STA', "location": '00', "channel": 'CCC'}
        if "response" in ref_trace:
            assert isinstance(ref_trace['response'], Response), 'ref_trace["response"] is not a Response object'
            stats["response"] = ref_trace["response"]
        trace = Trace(data=np.zeros(trace_pts, dtype=dtype), header=stats)
        return trace

    def _add_left_taper(self, fft, freqs, freq_lim, max_taper_len=100):
        n_zeros = len(fft[freqs < freq_lim])
        if n_zeros == 0:
            return fft
        fft[:n_zeros] = 0.
        if n_zeros > max_taper_len:
            taper_len = max_taper_len
        else:
            taper_len = n_zeros
        i_left = n_zeros - taper_len
        fft[i_left: n_zeros] = fft[n_zeros + 1] * self._taper_left(taper_len)
        return fft

    def _add_right_taper(self, fft, freqs, freq_lim, max_taper_len=100):
        n_zeros = len(fft[freqs > freq_lim])
        if n_zeros == 0:
            return fft
        fft[-n_zeros:] = 0.
        if n_zeros > max_taper_len:
            taper_len = max_taper_len
        else:
            taper_len = n_zeros
        i_right = -(n_zeros - taper_len + 1)
        # print(f'{taper_len=}, {i_right=}, {-n_zeros=}')
        fft[-(n_zeros + 1): i_right] = fft[-(n_zeros + 1)] * self._taper_right(taper_len)
        return fft

    @staticmethod
    def _taper_right(npts):
        return np.kaiser(2 * npts, 14)[npts:]

    @staticmethod
    def _taper_left(npts):
        return np.kaiser(2 * npts, 14)[:npts]

    @classmethod
    def from_amp_var_period(cls, amplitude, variance, period,
                            value_units="unknown"):
        """
        Generate a class object from an amplitude, a variance and a central period
        """
        raise ValueError('from_amp_var_period() not yet written')


if __name__ == "__main__":
    psd = PSDVals([0.001, 0.003, 0.006, 0.01, 0.02, 0.05, 0.1],
                  [-130, -160, -170, -175, -175, -180, -180],
                  'dB ref (m/s^2)^2/Hz')
    print(psd)
    psd.plot()
