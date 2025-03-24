import random

import numpy as np
from matplotlib import pyplot as plt


class PSDVals():
    """
    Holds PSD frequencies and values
    """
    def __init__(self, freqs, values, value_units="unknown"):
        """
        Args:
            freqs (list or :class:`np.array`): Frequencies (Hz)
            values (list or :class:`np.array`): PSD values in dB
            value_unites (str): Units of values
        """
        if isinstance(freqs, (list, tuple)):
            freqs = np.array(freqs)
        else:
            assert isinstance(freqs, np.ndarray)
        assert np.all(np.diff(freqs) > 0), 'freqs are not monotonically increasing'
        if isinstance(values, (list, tuple)):
            values = np.array(values)
        else:
            assert isinstance(values, np.ndarray)
        assert len(freqs) == len(values)
        assert isinstance(value_units, str)
        self.freqs = freqs
        self.values = values
        self.value_units = value_units
        
    def __str__(self):
        s = f"PSDVals (units={self.value_units}):\n"
        s += 20*" " + "     freq  |   value  \n"
        s += 20*" " + "---------- | ---------\n"
        for f, v in zip(self.freqs, self.values):
            s += 20*" " + f"{f:9.4g}  | {v:9.4g}\n"
        return s
    
    def resample(self, new_f):
        self.values = np.interp(new_f, self.freqs, self.values)
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
        return 360.*rng.random(n_values)

    def resample_values(self, freqs):
        """
        Resample values at the given frequencies
        """
        return np.interp(freqs, self.freqs, self.values)

    def as_fft(self, freqs, left='taper', right='taper', phases=None,
               plotit=False):
        """
        Return an fft "equivalent" to the given Power Spectral Density
        Args:
            freqs (list, np.array or None): Resample at the given freqs
            left (None, float or 'taper'): how to handle values below the
                lowest self.freq:
                    - None: use value at lowest input frequency
                    - float: set to the given value
                    - 'taper': taper using Kaiser function
            right (None, float, or 'taper'): how to handle values above the
                highest self.freq
            phases (np.array or None): force phases to be the given values
                (must be same length as frequencies)
        """
        if not freqs[0] == 0:
            raise ValueError(f"Cannot create an fft without f[0] == 0")
        fdiffs = np.diff(freqs)
        if not np.all(np.abs(fdiffs-fdiffs[0]) < fdiffs[0]/1e6):
            raise ValueError(f"freqs are not evenly spaced")

        fft = np.power(10., np.interp(freqs, self.freqs, self.values)/20)

        fft = np.nan_to_num(fft)        
        if left == 'taper':
            fft = self._add_left_taper(fft, freqs, self.freqs[0])
        elif left is not None:
            fft[freqs<self.freqs[0]] = left
        if right == 'taper':
            fft = self._add_right_taper(fft, freqs, self.freqs[-1])
        elif right is not None:
            fft[freqs>self.freqs[-1]] = right
            
        if plotit is True:
            f,a = plt.subplots()
            a.semilogx(freqs, fft, 'b')
            a.axvline(self.freqs[0], color='g', ls='--')
            a.axvline(self.freqs[-1], color='g', ls='--')
            a.semilogx(freqs, fft, 'r')
            plt.show()
        fft[0] = 0. # Make sure the zero-frequency value is zero
        # Scale for sample rate and window length
        sampling_rate = 1/(2*freqs[-1])
        # Bendata&Piersol 1986 eqs 11.100 & 11.102
        mul_factor = np.sqrt(len(freqs) * sampling_rate)  / 2    
        # mul_factor = (len(freqs)) / 2 
        # mul_factor = (sampling_rate) / 2 
        # div_factor = 2 / (sampling_rate) 
        # div_factor = 1    
        fft *= mul_factor
        if phases is None:
            phases = np.radians(self._random_phases(len(fft)))
        return fft * np.exp(1j * phases)


    def _add_left_taper(self, fft, freqs, freq_lim, max_taper_len=100):
        n_zeros = len(fft[freqs<freq_lim])
        if n_zeros == 0:
            return(fft)
        fft[:n_zeros] = 0.
        if n_zeros > max_taper_len:
            taper_len = max_taper_len
        else:
            taper_len = n_zeros
        i_left = n_zeros-taper_len
        fft[i_left: n_zeros] = fft[n_zeros+1]*self._taper_left(taper_len)
        return fft

    def _add_right_taper(self, fft, freqs, freq_lim, max_taper_len=100):
        n_zeros = len(fft[freqs > freq_lim])
        if n_zeros == 0:
            return(fft)
        fft[-n_zeros:] = 0.
        if n_zeros > max_taper_len:
            taper_len = max_taper_len
        else:
            taper_len = n_zeros
        i_right = -(n_zeros-taper_len+1)
        # print(f'{taper_len=}, {i_right=}, {-n_zeros=}')
        fft[-(n_zeros+1): i_right] = fft[-(n_zeros+1)]*self._taper_right(taper_len)
        return fft
    
    @staticmethod
    def _taper_right(npts):
        return np.kaiser(2*npts, 14)[npts:]

    @staticmethod
    def _taper_left(npts):
        return np.kaiser(2*npts, 14)[:npts]

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

