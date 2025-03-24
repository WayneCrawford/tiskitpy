from copy import deepcopy

from scipy.fft import irfft
import numpy as np
from matplotlib import pyplot as plt
from obspy.core.stream import Trace, Stream
from obspy.core import UTCDateTime

from ..spectral_density import SpectralDensity
from .compliance_functions import gravd, calc_norm_compliance
from .earth_model import EarthModel1D
from .tide_coefficients import TideCoefficients
from .psd_vals import PSDVals


class ComplianceNoise(object):
    """
    Generate synthetic seismological data based on environmental and noise factors
    """
    def __init__(self, water_depth, Z_offset_angles=(1, 10),
                 IG_m_seasurface=None, noise_pressure=None, noise_seismo=None,
                 noise_tilt_max=None, noise_tilt_variance=60, earth_model=None):
        """
        Return seismo and DPG time series corresponding to compliance plus noise

        Args:
            water_depth (numeric): water depth in meters
            Z_offset_angles (list): Seismometer's Z offset [angle, azimuth]
                from vertical, in degrees: (angle is the most important)
            IG_m_seasurface (:class:`PSDVals`, None): representation of
                infragravity wave PSD levels. Values are wave heights in
                dB ref m^2/Hz.
                If None, uses ComplianceNoise().default_IG_m_seasurface
            noise_pressure (:class:`PSDVals`, None): representation of DPG noise levels.
                Values are in dB ref Pa^2/Hz
                If None, uses ComplianceNoise().default_noise_pressure
            noise_seismo (:class:`PSDVals`, None): representation of seismometer noise levels.
                Values are in dB ref (m/s^2)^2/Hz
                If None, uses ComplianceNoise().default_noise_seismo
            noise_tilt_max (:class: `PSDVals`, None): maximum tilt noise levels (dB ref
                1 rad^2/Hz)
                If None, uses ComplianceNoise().default_tilt_max
            noise_tilt_variance (float): variance in dB of tilt noise levels
            earth_model (:class:`EarthModel1D`, None): 1D Earth model (last row
                treated as a half-space)
                If None, uses ComplianceNoise().default_earth_model
        """
        self.water_depth = water_depth
        self.Z_offset_angles = Z_offset_angles
        self.IG_m_seasurface = IG_m_seasurface
        self.noise_pressure = noise_pressure
        self.noise_seismo = noise_seismo
        self.noise_tilt_max = noise_tilt_max
        self.noise_tilt_variance = noise_tilt_variance
        self.earth_model = earth_model

        # Fill empty parameters with default values
        if self.IG_m_seasurface is None:
            self.IG_m_seasurface = self.default_IG_m_seasurface
        if self.noise_pressure is None:
            self.noise_pressure = self.default_noise_pressure
        if self.noise_seismo is None:
            self.noise_seismo = self.default_noise_seismo
        if self.noise_tilt_max is None:
            self.noise_tilt_max = self.default_tilt_max
        if self.earth_model is None:
            self.earth_model = self.default_earth_model

        # Validate variables
        assert isinstance(self.water_depth, (int, float))
        assert isinstance(self.IG_m_seasurface, PSDVals)
        assert isinstance(self.noise_pressure, PSDVals)
        assert isinstance(self.noise_seismo, PSDVals)
        assert isinstance(self.noise_tilt_max, PSDVals)
        assert isinstance(self.noise_tilt_variance, (int, float))
        assert isinstance(self.earth_model, EarthModel1D)

        self.IG_Pa_seafloor = self._calc_IG_Pa_seafloor()

    def __str__(self):
        s = '<ComplianceNoise>:'
        s += f'      water_depth={self.water_depth}'
        s += f'      Z_offset_angles={self.Z_offset_angles}\n'
        s += f'      IG_m_seasurface={self.IG_m_seasurface}\n'
        s += f'      noise_pressure={self.noise_pressure}\n'
        s += f'      noise_seismo={self.noise_seismo}\n'
        s += f'      noise_tilt_max={self.noise_tilt_max}\n'
        s += f'      tilt_min = noise_tilt_max - {self.noise_tilt_variance} dB\n'
        s += f'      earth_model={self.earth_model}'
        return s

    def _calc_IG_Pa_seafloor(self, min_f_spacing=0.003):
        psd = deepcopy(self.IG_m_seasurface)
        if np.any(np.diff(psd.freqs) > min_f_spacing):
            psd.resample(np.arange(psd.freqs[0],
                         psd.freqs[-1] + min_f_spacing * .999,
                         min_f_spacing))
        k = gravd(2 * np.pi * psd.freqs, self.water_depth)
        psd.values += 100 - self._cosh_dBs(k * self.water_depth)
        psd.value_units = 'dB ref 1 Pa^2/Hz'
        return psd

    def _cosh_dBs(self, x, max_input=700):
        # protect against values that are too big
        x[x > max_input] = max_input
        x[x < -max_input] = -max_input
        return to_DBs(np.cosh(x))

    @property
    def velPSD_compliance(self):
        om, k, ncompl = self._calc_ncompl()
        psd = deepcopy(self.IG_Pa_seafloor)
        psd.values += to_DBs(om * ncompl / k)
        return psd

    @property
    def default_IG_m_seasurface(self):
        """Default IG wave levels"""
        return PSDVals([0.001, 1], to_DBs([.002, .002]), "dB ref 1 m^2/Hz")

    @property
    def default_noise_pressure(self):
        """Default DPG noise levels"""
        return PSDVals([0.001, 0.003, 0.006, 0.01, 0.02, 0.05, 0.1],
                       [60, 30, -0, -10, -10, -10, -10],
                       'dB ref Pa^2/Hz')

    @property
    def default_noise_seismo(self):
        """Default seismometer noise levels"""
        return PSDVals([0.001, 0.003, 0.006, 0.01, 0.02, 0.05, 0.1],
                       [-130, -160, -170, -175, -175, -180, -180],
                       'dB ref (m/s^2)^2/Hz')

    @property
    def default_tilt_max(self):
        """Default noise_tilt_max values"""
        f = np.arange(0.001, 0.1, 0.001)
        hi = to_DBs(np.power(10., -6.5)*np.power(f, -1.5))
        return PSDVals(f.tolist(), hi.tolist(), 'dB ref 1 rad^2/Hz')

    @property
    def default_earth_model(self):
        """Default earth model"""
        return EarthModel1D([[1000, 3000, 3000, 1600],
                             [1000, 3000, 4000, 2300],
                             [1000, 3000, 5000, 2800],
                             [3000, 3000, 7500, 4300],
                             [3000, 3000, 8200, 4700]])

    def _default_tilt_ts(self, ref_trace, noise_tilt_variance):
        """
        returns default earth tides normalized to [0,1]
        and angle oscillating from 90 to 110

        Args:
            ref_trace (:class:`obspy.Trace`): a Trace used for sampling rate
                and time limits
            noise_tilt_variance (float): maximum dB variance.  Used to change tilt
                range from [0, 1] to [10^-[dB/20], 1]
        """
        assert isinstance(ref_trace, Trace)
        amp_trace, angle_trace = self.make_tilt_ts(ref_trace, TideCoefficients(), (90, 110))
        min_val = 10**(-noise_tilt_variance / 20)
        amp_trace.data[amp_trace.data < min_val] = min_val
        return amp_trace, angle_trace

    @property
    def accelPSD_compliance(self):
        om, k, ncompl = self._calc_ncompl()
        psd = deepcopy(self.IG_Pa_seafloor)
        psd.values = psd.values + to_DBs(om * om * ncompl / k)
        return psd

    def _calc_ncompl(self):
        f = self.IG_Pa_seafloor.freqs
        om = 2 * np.pi * f
        k = gravd(om, self.water_depth)
        ncompl = calc_norm_compliance(self.water_depth, f, self.earth_model)
        return om, k, ncompl
        
    def save_compliance(self, max_freq=None,
                        filename="model_compliance_Pa-1.csv", out_dir=None):
        """
        Saves self.earth_model's compliance to a BRUIT-FM CSV file
        
        Args:
            max_freq (float): only save up to this frequency (Hz)
            out_dir(str): output directory
            filename (str): output filename
        """
        oms, ks, ncompls = self._calc_ncompl()
        if out_dir is not None:
            filename = str(Path(out_dir) / filename)
        freqs = oms/(2*np.pi)
        if max_freq is not None:
            ncompls = ncompls[freqs <= max_freq]
            freqs = freqs[freqs <= max_freq]
        with open(filename, "w") as fid:
            fid.write('frequencies;compliance;uncertainty;phase\n')
            for freq, ncompl in zip(freqs, ncompls):
                fid.write('{:.5g};{:.5g};{:.5g};{:.5g}\n'
                          .format(freq, np.abs(ncompl), 0.000,
                                  np.angle(ncompl, deg=True)))
        

    def plot(self, fmin=0.001, fmax=0.1, fstep=0.001, outfile=None):
        """
        Plot the spectral representation of the noise sources in dB
        """
        f = np.arange(fmin, fmax + fstep / 2, fstep)
        # Plot
        fig, axs = plt.subplots(2, 1, sharex='col')
        # Plot the pressure signal
        axs[0].semilogx(f, self.IG_Pa_seafloor.resample_values(f), 'r', label='s_IG')
        axs[0].semilogx(f, self.noise_pressure.resample_values(f), 'b', label='n_press')
        axs[0].set_ylabel('dB ref 1 Pa^2/Hz')
        axs[0].legend()
        axs[0].set_ylim(-20, 60)
        # Plot the accel
        axs[1].semilogx(f, self.accelPSD_compliance.resample_values(f),
                        'r', label='s_IG')
        axs[1].semilogx(f, self.noise_seismo.resample_values(f),
                        'b', label='n_seismo')
        Z_angle_factor = to_DBs(np.sin(np.radians(self.Z_offset_angles[0])))
        axs[1].semilogx(f, self.noise_tilt_max.resample_values(f) + Z_angle_factor,
                        'g--', label='n_tilt_Z_max')
        axs[1].semilogx(f, self.noise_tilt_max.resample_values(f) + Z_angle_factor - self.noise_tilt_variance,
                        'g-.', label='n_tilt_Z_min')
        axs[1].set_ylabel('dB ref 1 (m/s^2)^2/Hz')
        axs[1].set_xlabel('Frequency (Hz)')
        axs[1].set_ylim(-200, -100)
        axs[1].legend()
        if outfile is not None:
            plt.savefig(outfile)
        plt.show()

    def streams(self, ref_trace, tilt_ts=None, s_sensitivity=3774870000,
                p_sensitivity=495, network='XX', station='SSSSS', plotit=False,
                forceInt32=False):
        """
        Return streams generated from to the noise and signal levels

        Args:
            ref_trace (:class: `obspy.Trace`): trace with time base to use
            tilt_ts (:class: `obspy.stream.Stream`): a tilt time series between the minimum
                and maximum values.  Must have the same sample rate and time span as
                the values defined by start_date, end_date and sps.
                The first trace is the tilt level modifier, which
                should vary between 0 and 1 (otherwise it will be shifted/compressed to
                do so).  Second trace is the direction of the tilt, in degrees clockwise
                from the seismometer's N channel.  Values must be between 0 and 360.
            s_sensitivity (float): Desired seismometer sensitivity (counts/m/s)
            p_sensitivity (float): Desired pressure gauge sensitivity (counts/Pa)

        Returns:
            streams (list):
                data (:class:`obspy.Stream): synthetic seafloor BB 4C data
                sources (:class:`obspy.Stream`): individual noises and signals
        """        
        trace_pts = ref_trace.stats.npts
        npts = 2**int(np.ceil(np.log2(trace_pts)))
        f = np.linspace(0, ref_trace.stats.sampling_rate / 2, npts)

        tr = []
        ref_trace = ref_trace.copy()    # Avoid overwriting original
        ref_trace.stats.station = station
        ref_trace.stats.network = network

        # IG wave pressure signal
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "IGP"
        p_fft = self.IG_Pa_seafloor.as_fft(f, plotit=False)
        tr[-1].data = irfft(p_fft)[:trace_pts]

        # DPG noise model
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "NOP"
        fft = self.noise_pressure.as_fft(f)
        tr[-1].data = irfft(fft)[:trace_pts]

        # Vertical compliance signal 
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "IGZ"
        z_fft = self.velPSD_compliance.as_fft(f, phases=np.angle(p_fft))
        tr[-1].data = irfft(z_fft)[:trace_pts]

        # Seismometer noise model
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "NOS"
        fft = self.noise_seismo.as_fft(f)
        tr[-1].data = irfft(fft)[:trace_pts]

        # Tilt noise model
        fft = self.noise_tilt_max.as_fft(f)
        noise_max = irfft(fft)[:trace_pts]
        dyntilt_amp, dyntilt_angle = self._default_tilt_ts(ref_trace, self.noise_tilt_variance)
        noise = noise_max * dyntilt_amp.data
        N_noise = np.cos(np.radians(dyntilt_amp)) * noise
        E_noise = np.sin(np.radians(dyntilt_amp)) * noise
        angfact = np.sin(np.radians(self.Z_offset_angles[0]))
        azefact_1 = np.sin(np.radians(self.Z_offset_angles[1]))
        azefact_2 = np.cos(np.radians(self.Z_offset_angles[1]))
        Z_noise = angfact * (azefact_1 * N_noise + azefact_2 * E_noise)
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "NO1"
        tr[-1].data = N_noise
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "NO2"
        tr[-1].data = E_noise
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "NOZ"
        tr[-1].data = Z_noise

        sources = Stream(traces=tr)
        if plotit is True:
            sources.plot(equal_scales=False)

        # Add signal and noise time series to make synthetic BBOBS channels
        tr = []
        # LDG: Differential pressure gauge
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "LDG"
        tr[-1].data = ((sources.select(channel='IGP')[0].data
                        + sources.select(channel='NOP')[0].data)
                       * p_sensitivity)
        # LH1: N-equivalent horizontal seismometer channel
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "LH1"
        tr[-1].data = ((sources.select(channel='NOS')[0].data
                        + sources.select(channel='NO1')[0].data)
                       * s_sensitivity)
        # LH2: E-equivalent horizontal seismometer channel
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "LH2"
        tr[-1].data = ((sources.select(channel='NOS')[0].data
                        + sources.select(channel='NO2')[0].data)
                       * s_sensitivity)
        # LHA: Vertical seismometer channel
        tr.append(ref_trace.copy())
        tr[-1].stats.channel = "LHZ"
        tr[-1].data = ((sources.select(channel='NOS')[0].data
                        + sources.select(channel='IGZ')[0].data
                        + sources.select(channel='NOZ')[0].data)
                       * s_sensitivity)

        data = Stream(traces=tr)
        if forceInt32 is True:
            for tr in data:
                tr.data  = np.require(tr.data, dtype=np.int32)
        if plotit is True:
            data.plot(equal_scales=False)

        return data, sources

    @staticmethod
    def make_tilt_ts(ref_trace, coefficients, angle_limits=(70, 100)):
        """
        make a simple tilt time series summing signals of given periods,
        amplitudes and starting phases

        Args:
            ref_trace (:class: `obspy.core.Trace`): A trace covering the
                desired period and with the desired sample rate
            coefficients (TideCoefficients): the tidal coefficients
            angle_limits (list): current direction angle limits ([min, max]
                in degrees).  Will be modulated by tide levels
        """
        tide_trace = coefficients.make_trace(ref_trace)
        # normalize between 0 and 1
        tide_trace.data -= np.min(tide_trace.data)
        tide_trace.data /= np.max(tide_trace.data)
        angles_trace = tide_trace.copy()
        angle_range = abs(angle_limits[1] - angle_limits[0])
        angles_trace.data = (angles_trace.data * angle_range) + min(angle_limits)
        return tide_trace, angles_trace


def to_DBs(inp):
    """
    Converts values to dBs ref 1

    Args:
        inp (float, list, or np.ndarray): values to convert
    """
    if isinstance(inp, (list, tuple)):
        inp = np.array(inp)
    return 20 * np.log10(inp)


def from_DBs(inp):
    """
    Converts values from dBs ref 1

    Args:
        inp (float, list or np.ndarray): values to convert
    """
    if isinstance(inp, (list, tuple)):
        inp = np.array(inp)
    return np.power(10., inp / 20)


if __name__ == "__main__":
    # Show an example
    wdepth = 2000
    noise_model = ComplianceNoise(wdepth)
    noise_model.plot(outfile='noise_model.png')
    noise_model.save_compliance(max_freq=0.07)
    resp_trace = Trace(np.zeros(86400*5),
                       header={'sample_rate': 1,
                               'starttime': UTCDateTime('2024-01-01T00:00:00')})
    data, sources = noise_model.streams(resp_trace)
    data.plot(equal_scale=False)
    sources.select(channel='ZIG').plot(method="full")
    data.write('synth_data.mseed', 'MSEED')
    sources.write('synth_sources.mseed', 'MSEED')
    sd_sources = SpectralDensity.from_stream(sources)
    sd_sources.plot()
    sd_data = SpectralDensity.from_stream(data)
    # sd_data.plot()
    sd_data.plot_coherences()
