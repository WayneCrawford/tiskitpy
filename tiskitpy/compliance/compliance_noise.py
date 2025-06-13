from scipy.fft import irfft
import numpy as np
from matplotlib import pyplot as plt
from obspy.core.stream import Trace, Stream
from obspy.core.inventory import Inventory, Network, Station, Channel, Response
from obspy.core import UTCDateTime
from pathlib import Path

from ..spectral_density import SpectralDensity
from .compliance_functions import gravd, calc_norm_compliance
from .earth_model import EarthModel1D
from .tide_coefficients import TideCoefficients
from .psd_vals import PSDVals

default_IG_m_seasurface = ([[0.001, .002], [1, .002]], False)
default_noise_pressure = ([[0.001, 60], [0.003, 30], [0.006, 0], [0.01, -10],
                           [0.02, -10], [0.05, -10], [0.1, -10], [1, -10]
                          ],
                          True)
default_noise_seismo = ([[0.001, -130], [0.003, -160], [0.006, -170], [0.01, -175],
                         [0.02, -175],  [0.05, -180], [ 0.1, -180],   [1, -170]
                        ],
                        True)
default_tilt_max = ([[f, np.power(10., -6.5) * np.power(f, -1.5)]
                     for f in np.power(10, np.arange(-3, 0.1, .25))
                    ],
                    False)
default_earth_model = [[1000, 3000, 3000, 1600],
                       [1000, 3000, 4000, 2300],
                       [1000, 3000, 5000, 2800],
                       [3000, 3000, 7500, 4300],
                       [3000, 3000, 8200, 4700]]


class ComplianceNoise(object):
    """
    Generate synthetic seismological data based on environmental and noise factors
    """
    def __init__(self, water_depth=2400,
                 Z_offset_angles=(2, 15),
                 IG_m_seasurface=default_IG_m_seasurface,
                 noise_pressure=default_noise_pressure,
                 noise_seismo=default_noise_seismo,
                 noise_tilt_max=default_tilt_max,
                 noise_tilt_variance=60,
                 noise_tilt_direction_limits=(100, 130),
                 earth_model=default_earth_model,
                 IG_freqstep=0.001):
        """
        Return seismo and DPG time series corresponding to compliance plus noise

        Args:
            water_depth (numeric): water depth in meters
            Z_offset_angles (list): Seismometer's Z offset [angle, azimuth]
                from vertical, in degrees: (angle is the most important)
            IG_m_seasurface (list, None): Infragravity wave PSD levels in format:
                ([[freq1, value1],
                 [freq2, value2],
                 ...
                 [freqN, valueN]],
                 is_dB)
                Values are wave heights in m (if is_dB is False) or in dB ref
                1 m^2/Hz.
                If None, uses ComplianceNoise().default_IG_m_seasurface
            noise_pressure (list, None): representation of DPG noise levels.
                Same format as for IG_m_seasurface, values are in Pa or dB
                equivalent
                If None, uses ComplianceNoise().default_noise_pressure
            noise_seismo (:class:`PSDVals`, None): representation of seismometer
                noise levels. Same format as for IG_m_seasurface, values are
                in m/s^2 or dB equivalent
                If None, uses ComplianceNoise().default_noise_seismo
            noise_tilt_max (:class: `PSDVals`, None): maximum tilt noise levels
                Same format as for IG_m_seasurface, values are in m/s^2 or dB
                equivalent.
                If None, uses ComplianceNoise().default_tilt_max
            noise_tilt_direction_limits (tuple): minimum and maximum tilt
                directions (degrees).
            noise_tilt_variance (float): variance in dB of tilt noise levels
            earth_model (list, None): 1D Earth model in the format:
                [[thick1, rho1, vp1, vs1]
                 [thick2, rho2, vp2, vs2]
                 ...
                 [thickN, rhoN, vpN, vsN]]
                wher units are meters, kg/m^2, m/s and m/s, and the last row
                is treated as a half-space)
                If None, uses ComplianceNoise().default_earth_model
            IG_freqstep (float): maximum frequency step for IG wave and compliance
                PSDs (must be small to capture shallow/deep water cutoff)
        """
        # Validate variables
        assert isinstance(water_depth, (int, float))
        assert isinstance(noise_tilt_variance, (int, float))

        self.water_depth = water_depth
        self.Z_offset_angles = Z_offset_angles
        self.IG_m_seasurface = PSDVals(IG_m_seasurface, "m")
        self.noise_pressure = PSDVals(noise_pressure, 'Pa')
        self.noise_seismo = PSDVals(noise_seismo, '(m/s^2)')
        self.noise_tilt_max = PSDVals(noise_tilt_max, '(m/s^2)')
        self.noise_tilt_direction_limits = noise_tilt_direction_limits
        self.noise_tilt_variance = noise_tilt_variance
        self.earth_model = EarthModel1D(earth_model)
        self.IG_freqstep = IG_freqstep

    @property
    def IG_Pa_seafloor(self):
        psd = self.IG_m_seasurface.copy()
        if np.any(np.diff(psd.freqs) > self.IG_freqstep):
            psd.resample(np.arange(psd.freqs[0],
                         psd.freqs[-1] + self.IG_freqstep * .999,
                         self.IG_freqstep))
        k = gravd(2 * np.pi * psd.freqs, self.water_depth)
        psd.values += 100 - self._cosh_dBs(k * self.water_depth)
        psd.value_units = 'dB ref 1 Pa^2/Hz'
        return psd

    @property
    def compliance_accel(self):
        """PSD of compliance * IG pressure, in (m/s^2)^2/Hz"""
        om, k, ncompl = self._calc_ncompl()
        ref = 'dB ref 1 Pa^2/Hz'
        assert self.IG_Pa_seafloor.value_units == ref, f"{self.IG_Pa_seafloor.value_units=} should be '{ref}'"
        psd = self.IG_Pa_seafloor.copy()
        psd.value_units = 'dB ref 1 (m/s^2)^2/Hz'
        psd.values = psd.values + to_DBs(om * om * abs(ncompl) / k)
        return psd

    @property
    def Z_angle_factor_DBs(self):
        """Rotation of horizontal tilt noise onto Z channel"""
        return to_DBs(np.sin(np.radians(self.Z_offset_angles[0])))

    @property
    def PSDs(self):
        """
        Return dictionary of all PSD values
        """
        return {'IGP': self.IG_Pa_seafloor,
                'NOP': self.noise_pressure,
                "IGZ": self.compliance_accel,
                "NOS": self.noise_seismo,
                "NOT_max": self.noise_tilt_max,
                "NOT_min": self.noise_tilt_max - self.noise_tilt_variance,
                "NOT_zmax": self.noise_tilt_max + self.Z_angle_factor_DBs,
                "NOT_zmin": self.noise_tilt_max + self.Z_angle_factor_DBs
                              - self.noise_tilt_variance}

    def __str__(self):
        s = '<ComplianceNoise>:\n'
        s += f'    water_depth={self.water_depth}\n'
        s += f'    Z_offset_angles={self.Z_offset_angles}\n'
        s += f'    IG_m_seasurface={self.IG_m_seasurface}\n'
        s += f'    noise_pressure={self.noise_pressure}\n'
        s += f'    noise_seismo={self.noise_seismo}\n'
        s += f'    noise_tilt_max={self.noise_tilt_max}\n'
        s += f'    tilt_min = noise_tilt_max - {self.noise_tilt_variance} dB\n'
        s += f'    earth_model={self.earth_model}'
        return s

    def _cosh_dBs(self, x, max_input=700):
        # protect against values that are too big
        x[x > max_input] = max_input
        x[x < -max_input] = -max_input
        return to_DBs(np.cosh(x))

    def norm_compliance(self, f=None):
        """
        Return normalized compliance of the object's EarthModel
        
        Args:
            f (list, np.array, None): frequencies at which to calculate.
                If None, then calculate at the frequencies of self.IG_Pa_seafloor
        """
        _, _, ncompl = self._calc_ncompl(f)
        return ncompl

    def _calc_ncompl(self, f=None):
        """
        Return normalized compliance of the object's EarthModel
        
        Args:
            f (list, np.array, None): frequencies at which to calculate.
                If None, then calculate at the frequencies of self.IG_Pa_seafloor
        Returns:
            (tuple): omega (np.array): angular frequencies
                     k (np.array): wavenumbers
                     ncompl (np.array): the normalized compliance
        """
        if f is None:
            f = self.IG_Pa_seafloor.freqs
        ncompl = calc_norm_compliance(self.water_depth, f, self.earth_model)
        # print(ncompl[~np.isnan(ncompl)])
        om = 2 * np.pi * f
        k = gravd(om, self.water_depth)
        # print(f'ComplianceNoise._calc_ncompl(): {self.water_depth=}, {om[:5]=}, {k[:5]=}')
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
        freqs = oms / (2 * np.pi)
        if max_freq is not None:
            ncompls = ncompls[freqs <= max_freq]
            freqs = freqs[freqs <= max_freq]
        with open(filename, "w") as fid:
            fid.write('frequencies;compliance;uncertainty;phase\n')
            for freq, ncompl in zip(freqs, ncompls):
                fid.write('{:.5g};{:.5g};{:.5g};{:.5g}\n'
                          .format(freq, np.abs(ncompl), 0.000,
                                  np.angle(ncompl, deg=True)))

    def plot(self, fmin=0.001, fmax=0.1, fstep=0.001, outfile=None, show=True):
        """
        Plot the spectral representation of the noise sources in dB
        """
        f = np.arange(fmin, fmax + fstep / 2, fstep)
        # Plot
        fig, axs = plt.subplots(2, 1, sharex='col')
        # Plot the pressure signal
        axs[0].semilogx(f, self.IG_Pa_seafloor.resample_values(f), 'r', label='IG_P')
        axs[0].semilogx(f, self.noise_pressure.resample_values(f), 'b', label='NO_P')
        axs[0].set_ylabel('dB ref 1 Pa^2/Hz')
        axs[0].legend()
        axs[0].set_ylim(-20, 60)
        # Plot the accel
        axs[1].semilogx(f, self.compliance_accel.resample_values(f),
                        'r', label='IG_S')
        axs[1].semilogx(f, self.noise_seismo.resample_values(f),
                        'b', label='NO_S')
        axs[1].semilogx(f,
                        self.noise_tilt_max.resample_values(f)+ self.Z_angle_factor_DBs,
                        'g--', label='NO_Z(max)')
        axs[1].semilogx(f,
                        self.noise_tilt_max.resample_values(f)
                        + self.Z_angle_factor_DBs - self.noise_tilt_variance,
                        'g-.', label='NO_Z(min)')
        axs[1].set_ylabel('dB ref 1 (m/s^2)^2/Hz')
        axs[1].set_xlabel('Frequency (Hz)')
        axs[1].set_ylim(-200, -100)
        axs[1].legend()
        if outfile is not None:
            plt.savefig(outfile)
        if show is True:
            plt.show()

    def streams(self, ref_trace, s_sensitivity=3774870000,
                p_sensitivity=495, network='XX', station='SSSSS', plotit=False,
                forceInt32=False):
        """
        Return streams generated from to the noise and signal levels
        Simply multiplies physical values by a sensitivity value, would be
        better to convolve with instrument response.

        Args:
            ref_trace (:class: `obspy.Trace`): trace with time base to use
            s_sensitivity (float): Desired seismometer sensitivity (counts/m/s)
            p_sensitivity (float): Desired pressure gauge sensitivity (counts/Pa)
            network (str): Network code (1-2 characters)
            station (str): Station code (1-5 characters)
            forceInt32 (bool): force output data to have dtype=np.int32

        Returns:
            streams (list):
                data (:class:`obspy.Stream): synthetic seafloor BB 4C data
                sources (:class:`obspy.Stream`): individual noises and signals
                inv (:class:`obspy.core.Inventory`): channel metadata
        """
        noise_to_vel = False    # Should be True, logically
        if noise_to_vel is True:
            print('noise PSDs converted from accel to vel before ifft')
        else:
            print('noise PSDs NOT converted from accel to vel')
        # SETUP
        # Validate inputs
        if not ref_trace.stats.channel[0] == 'L':
            raise ValueError("ref_trace channel code ({}) doesn't start with 'L'"
                .format(ref_trace.stats.channel))
        # Set up variables
        trace_pts = ref_trace.stats.npts
        npts = 2**int(np.ceil(np.log2(trace_pts)))
        sr = ref_trace.stats.sampling_rate
        location = ref_trace.stats.location
        channel = ref_trace.stats.channel
        f = np.linspace(0, ref_trace.stats.sampling_rate / 2, npts)
        velocity_response = Response.from_paz([], [], 1, input_units='m/s', output_units='count')
        accel_response = Response.from_paz([], [], 1, input_units='m/s**2', output_units='m/s**2')
        tr = []

        # Prepare a base_trace for the outputs
        base_trace = ref_trace.copy()    # Don't overwrite original
        base_trace.stats.station = station
        base_trace.stats.network = network
        if noise_to_vel is True:
            base_trace.stats.response = velocity_response
        else:
            base_trace.stats.response = accel_response

        # CREATE NOISE + IG/COMPLIANCE TRACES BY SOURCE
        # IG wave pressure signal
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "IGP"
        tr[-1].stats.response = None
        IG_fft = self.IG_Pa_seafloor.as_fft(f)
        tr[-1].data = irfft(IG_fft)[:trace_pts]

        # DPG noise model
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "NOP"
        tr[-1].stats.response = None
        tr[-1].data = irfft(self.noise_pressure.as_fft(f))[:trace_pts]

        # Vertical compliance signal
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "IGZ"
        tr[-1].stats.response = velocity_response
        # Phase_motion = Phase_pressure + 180° , Phase_velocity = Phase_motion + 90°
        # Ignores noise_to_vel because we validated accel_as_vel previously
        tr[-1].data = irfft(self.compliance_accel.accel_as_vel.as_fft(
            f, phases=np.angle(IG_fft)-np.pi/2))[:trace_pts]

        # Seismometer noise model
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "NOS"
        if noise_to_vel is True:
            tr[-1].data = irfft(self.noise_seismo.accel_as_vel.as_fft(f))[:trace_pts]
        else:
            tr[-1].data = irfft(self.noise_seismo.as_fft(f))[:trace_pts]
       
        # Tilt noise model
        if noise_to_vel is True:
            noise_max = irfft(self.noise_tilt_max.accel_as_vel.as_fft(f))[:trace_pts]
        else:
            noise_max = irfft(self.noise_tilt_max.as_fft(f))[:trace_pts]
        dyntilt_amp, dyntilt_angle = self.make_tilt_ts(base_trace, noise_max)
        angfact = np.sin(np.radians(self.Z_offset_angles[0]))
        azefact_1 = np.sin(np.radians(self.Z_offset_angles[1]))
        azefact_2 = np.cos(np.radians(self.Z_offset_angles[1]))
        N_noise = dyntilt_amp.data * np.cos(np.radians(dyntilt_angle.data))
        E_noise = dyntilt_amp.data * np.sin(np.radians(dyntilt_angle.data))
        Z_noise = angfact * (azefact_1 * N_noise + azefact_2 * E_noise)
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "NO1"
        tr[-1].data = N_noise
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "NO2"
        tr[-1].data = E_noise
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "NOZ"
        tr[-1].data = Z_noise

        sources = Stream(traces=tr)
        if plotit is True:
            sources.plot(equal_scales=False)

        # CREATE SYNTHETIC OBS CHANNELS WITH NOISE + IG/COMPLIANCE
        # Add signal and noise time series to make synthetic BBOBS channels
        tr = []
        base_trace.stats.response = None
        # LDG: Differential pressure gauge
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "LDG"
        tr[-1].data = ((sources.select(channel='IGP')[0].data
                        + sources.select(channel='NOP')[0].data)
                       * p_sensitivity)
        # LH1: N-equivalent horizontal seismometer channel
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "LH1"
        tr[-1].data = ((sources.select(channel='NOS')[0].data
                        + sources.select(channel='NO1')[0].data)
                       * s_sensitivity)
        # LH2: E-equivalent horizontal seismometer channel
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "LH2"
        tr[-1].data = ((sources.select(channel='NOS')[0].data
                        + sources.select(channel='NO2')[0].data)
                       * s_sensitivity)
        # LHZ: Vertical seismometer channel
        tr.append(base_trace.copy())
        tr[-1].stats.channel = "LHZ"
        tr[-1].data = ((sources.select(channel='NOS')[0].data
                        + sources.select(channel='IGZ')[0].data
                        + sources.select(channel='NOZ')[0].data)
                       * s_sensitivity)

        data = Stream(traces=tr)
        if forceInt32 is True:
            for tr in data:
                tr.data = np.require(tr.data, dtype=np.int32)
        if plotit is True:
            data.plot(equal_scales=False)

        # MAKE INVENTORY
        pressresp = Response.from_paz([], [], p_sensitivity, 1.0, 'm/s', 'count')
        seisresp = Response.from_paz([], [], s_sensitivity, input_units='m/s', output_units='count')
        # doesn't know 'Pa'
        pressresp.instrument_sensitivity.input_units='Pa' 
        pressresp.response_stages[0].input_units='Pa' 
        channels=[Channel('LHZ', location, 0, 0, 0, 0, response=seisresp, dip = -90),
                  Channel('LH1', location, 0, 0, 0, 0, response=seisresp),             
                  Channel('LH2', location, 0, 0, 0, 0, response=seisresp),             
                  Channel('LDG', location, 0, 0, 0, 0, response=pressresp)]
        stations = [Station(station, 0, 0, 0, channels=channels)]
        networks = [Network(network, stations=stations)]
        inv = Inventory(networks=networks)

        return data, sources, inv

    def make_tilt_ts(self, ref_trace, noise_max, coefficients=TideCoefficients()):
        """
        make a simple tilt time series summing signals of given periods,
        amplitudes and starting phases

        Args:
            ref_trace (:class: `obspy.core.Trace`): A trace covering the
                desired period and with the desired sample rate
            noise_max (:class:`obspy.core.Trace`): maximum tilt noise time series
            coefficients (TideCoefficients): the tidal coefficients
        """
        assert isinstance(ref_trace, Trace)
        tide_trace = coefficients.make_trace(ref_trace)
        # normalize between (-self.noise_tilt_variance dB) and 1
        tide_trace.data -= np.min(tide_trace.data)
        tide_trace.data /= np.max(tide_trace.data)
        min_val = 10**(-self.noise_tilt_variance / 20)
        tide_trace.data[tide_trace.data < min_val] = min_val

        amp_trace = tide_trace.copy()
        amp_trace.data *= noise_max

        angles_trace = tide_trace.copy()
        angle_range = abs(self.noise_tilt_direction_limits[1] - self.noise_tilt_direction_limits[0])
        angles_trace.data = (angles_trace.data * angle_range) + min(self.noise_tilt_direction_limits)

        return amp_trace, angles_trace


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
    resp_trace = Trace(np.zeros(86400 * 5),
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
