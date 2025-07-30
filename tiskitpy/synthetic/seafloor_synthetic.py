"""
Class to generate seafloor seismological noise models

Lacks:
- Microseisms
- 2nd order tilt effects on the vertical channel
- non-stationary IG wave spectra
- Gravitational attraction effects
"""
import numpy as np
from matplotlib import pyplot as plt
from obspy.core.stream import Trace, Stream
from obspy.core.inventory import Inventory, Network, Station, Channel, Response
from obspy.core import UTCDateTime
# from pathlib import Path

from ..spectral_density import SpectralDensity
from ..compliance import Compliance, EarthModel1D
from .tide_coefficients import TideCoefficients
from .psd_vals import PSDVals
from .functions import to_DBs  # , from_DBs

default_water_depth = 2400
default_Z_offset_angles = (2, 15)  # angle from vertical, azimuth from N
default_IG_m_seasurface = ([[0.001, .02], [1, .02]], False)
default_noise_pressure = ([[0.001, 60], [0.003, 30], [0.006, 0], [0.01, -10],
                           [0.02, -10], [0.05, -10], [0.1, -10], [1, -10]],
                          True)
default_noise_seismo = ([[0.001, -130], [0.003, -160], [0.006, -170],
                         [0.01, -175], [0.02, -175],  [0.05, -180],
                         [0.1, -180],   [1, -170]],
                        True)
default_tilt_max = PSDVals.sloped_freqs_and_values(-145, -30, -3, 0.1, .25)
default_tilt_variance = 50  # dB
default_tilt_direction_limits = (100, 130)  # degrees from "N"
default_earth_model = [[1000, 3000, 3000, 1600],
                       [1000, 3000, 4000, 2300],
                       [1000, 3000, 5000, 2800],
                       [3000, 3000, 7500, 4300],
                       [3000, 3000, 8200, 4700]]


class SeafloorSynthetic(object):
    """
    Generate synthetic seismological data based on environmental and noise factors

    Attributes:
        water_depth (float): water depth in meters
        Z_offset_angles (list): Seismometer's Z offset [angle, azimuth]
            from vertical, in degrees
        IG_m_seasurface (:class:`tiskitpy.PSDVals`): Infragravity wave PSD levels (ref m)
        noise_pressure (:class:`tiskitpy.PSDVals`): Pressure sensor noise levels (ref Pa)
        noise_seismo (:class:`tiskitpy.PSDVals`): Seismometer noise levels (ref m/s^2)
        noise_tilt_max (:class:`tiskitpy.PSDVals`): Maximum tilt noise (ref m/s^2))
        noise_tilt_direction_limits (tuple): minimum and maximum tilt
            directions (degrees).
        noise_tilt_variance (float): variance in dB of tilt noise levels
        earth_model (:class:`tiskitpy.EarthModel1D`): 1D Earth model
        IG_freqstep (float): maximum frequency step for IG wave and compliance
            PSDs
    """
    def __init__(self, water_depth=default_water_depth,
                 Z_offset_angles=default_Z_offset_angles,
                 IG_m_seasurface=default_IG_m_seasurface,
                 noise_pressure=default_noise_pressure,
                 noise_seismo=default_noise_seismo,
                 noise_tilt_max=default_tilt_max,
                 noise_tilt_variance=default_tilt_variance,
                 noise_tilt_direction_limits=default_tilt_direction_limits,
                 earth_model=default_earth_model,
                 IG_freqstep=0.001):
        """
        Return seismo and DPG time series corresponding to compliance plus noise

        Args:
            water_depth (numeric): water depth in meters
            Z_offset_angles (list): Seismometer's Z offset [angle, azimuth]
                from vertical, in degrees: (angle is the most important)
            IG_m_seasurface (tuple): Infragravity wave PSD levels in format:
                ([[freq1, value1],
                 [freq2, value2],
                 ...
                 [freqN, valueN]],
                 is_dB)
                Where is_dB is bool.  Values are wave heights in m (if is_dB
                is False) or in dB ref 1 (m^2)/Hz.
            noise_pressure (tuple): representation of DPG noise levels.
                Same format as for IG_m_seasurface, values are in Pa or dB
                equivalent
            noise_seismo (tuple): representation of seismometer
                noise levels. Same format as for IG_m_seasurface, values are
                in m/s^2 or dB equivalent
            noise_tilt_max (tuple): maximum tilt noise levels
                Same format as for IG_m_seasurface, values are in m/s^2 or dB
                equivalent.
            noise_tilt_direction_limits (tuple): minimum and maximum tilt
                directions (degrees).
            noise_tilt_variance (float): variance in dB of tilt noise levels
            earth_model (list, None): 1D Earth model in the format:
                [[thick1, rho1, vp1, vs1]
                 [thick2, rho2, vp2, vs2]
                 ...
                 [thickN, rhoN, vpN, vsN]]
                where units are meters, kg/m^2, m/s and m/s, and the last row
                is treated as a half-space)
            IG_freqstep (float): maximum frequency step for IG waves and
                compliance PSDs (must be small enough to capture shallow/deep
                water cutoff)
        """
        # Validate variables
        assert isinstance(water_depth, (int, float))
        assert isinstance(noise_tilt_variance, (int, float))

        self.water_depth = water_depth
        self.Z_offset_angles = Z_offset_angles
        self.IG_m_seasurface = PSDVals(IG_m_seasurface, "m")
        self.noise_pressure = PSDVals(noise_pressure, 'Pa')
        self.noise_seismo = PSDVals(noise_seismo, 'm/s^2')
        self.noise_tilt_max = PSDVals(noise_tilt_max, 'm/s^2')
        self.noise_tilt_direction_limits = noise_tilt_direction_limits
        self.noise_tilt_variance = noise_tilt_variance
        self.earth_model = EarthModel1D(earth_model)
        self.IG_freqstep = IG_freqstep

    @property
    def IG_Pa_seafloor(self):
        """
        Infragravity wave seafloor pressure PSD

        Based on self.IG_m_seasurface and self.water_depth)
        """
        # Seawater density is 1020-1029 at sea surface and up to 1050 at
        # the deep seafloor
        seawater_density = 1030
        g = 9.81  # 9.78 at equator, 9.83 at poles
        psd = self.IG_m_seasurface.copy()
        # Resample if frequency spacing is larger than the IG freqstep
        if np.any(np.diff(psd.freqs) > self.IG_freqstep):
            psd.resample(np.arange(psd.freqs[0],
                         psd.freqs[-1] + self.IG_freqstep * .999,
                         self.IG_freqstep))
        k = Compliance.gravd(2 * np.pi * psd.freqs, self.water_depth)
        psd.values += 20*np.log10(seawater_density*g)   # meters to Pascals
        psd.values -= self._cosh_dBs(k * self.water_depth)  # depth decay
        psd.value_units = 'dB ref 1 Pa^2/Hz'
        return psd

    @property
    def compliance_accel(self):
        """PSD of compliance * IG pressure, in (m/s^2)^2/Hz"""
        om, k, ncompl = self._calc_ncompl()
        ref = 'dB ref 1 Pa^2/Hz'
        if not self.IG_Pa_seafloor.value_units == ref:
            raise ValueError("IG_Pa_seafloor.value_units={} are not '{}'"
                             .format(self.IG_Pa_seafloor.value_units, ref))
        psd = self.IG_Pa_seafloor.copy()
        psd.value_units = 'dB ref 1 (m/s^2)^2/Hz'
        psd.values = psd.values + to_DBs(om * om * abs(ncompl) / k)
        return psd

    @property
    def Z_angle_factor_DBs(self):
        """Rotation of horizontal tilt noise onto Z channel"""
        return to_DBs(np.sin(np.radians(self.Z_offset_angles[0])))

    @property
    def stream_source_codes(self):
        """ return dict of streams and each trace's source codes """
        return {'LH1': ('NOS', 'NT1'),
                'LH2': ('NOS', 'NT2'),
                'LHZ': ('NOS', 'NTZ', 'IGZ'),
                'LDG': ('NOP', 'IGP')}

    @property
    def source_codes(self):
        """ return list of source codes """
        return ['IGP', 'NOP',  "IGZ", "NOS",
                "NTH_max", "NTH_min", "NTZ_max", "NTZ_min"]

    @property
    def trace_source_codes(self):
        """ return list of trace source codes """
        return ['IGP', 'NOP',  "IGZ", "NOS", "NT1", "NT2", "NTZ"]

    @property
    def PSDs(self):
        """
        Dictionary of all seafloor PSDs
        """
        return {k: self.source_by_code(k) for k in self.source_codes}

    def source_by_code(self, ch_code):
        """
        Return PSD by code

        Args:
            code (str): a 3-letter code that points to the given source.
                        If it's got more than three letters, its a combination
                        source
        Returns:
            :class:`PSDVals`: the source PSD
        """
        match ch_code:
            case 'IGP':
                return self.IG_Pa_seafloor
            case 'NOP':
                return self.noise_pressure
            case 'IGZ':
                return self.compliance_accel
            case 'NOS':
                return self.noise_seismo
            case 'NTH_max':
                return self.noise_tilt_max
            case 'NTH_min':
                return self.noise_tilt_max - self.noise_tilt_variance
            case 'NTZ_max':
                return self.noise_tilt_max + self.Z_angle_factor_DBs
            case 'NTZ_min':
                return (self.noise_tilt_max + self.Z_angle_factor_DBs
                        - self.noise_tilt_variance)
            case _:
                raise ValueError(f'"{ch_code}" is not a valide source channel code')
        return

    def __str__(self):
        s = '<SeafloorSynthetic>:\n'
        s += f'    water_depth={self.water_depth}\n'
        s += f'    Z_offset_angles={self.Z_offset_angles}\n'
        s += f'    IG_m_seasurface={self.IG_m_seasurface}\n'
        s += f'    noise_pressure={self.noise_pressure}\n'
        s += f'    noise_seismo={self.noise_seismo}\n'
        s += f'    noise_tilt_max={self.noise_tilt_max}\n'
        s += f'    noise_tilt_min = noise_tilt_max - {self.noise_tilt_variance} dB\n'
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
                If None, then calculate at the frequencies of
                self.IG_Pa_seafloor
        """
        _, _, ncompl = self._calc_ncompl(f)
        return ncompl

    def _calc_ncompl(self, f=None):
        """
        Return normalized compliance of the object's EarthModel

        Args:
            f (list, np.array, None): frequencies at which to calculate.
                If None, then calculate at the frequencies of
                self.IG_Pa_seafloor
        Returns:
            (tuple): omega (np.array): angular frequencies
                     k (np.array): wavenumbers
                     ncompl (np.array): the normalized compliance
        """
        if f is None:
            f = self.IG_Pa_seafloor.freqs
        ncompl = Compliance.calc_norm_compliance(self.water_depth, f,
                                                 self.earth_model)
        om = 2 * np.pi * f
        k = Compliance.gravd(om, self.water_depth)
        return om, k, ncompl

    def save_compliance(self, max_freq=True, base_name="model", out_dir=None):
        """
        GRANDFATHERED: use Compliance.write() instead
        Saves self.earth_model's compliance to a file

        Args:
            max_freq (float or bool): if float, only save up to given freq (Hz).
                If true, cut off at water depth based predicted compliance cutoff
                If false, use all of  IG_Pa_seafloor.freqs
            out_dir(str or :class:`Path`): output directory
            filename (str): output filename
        """
        ncompl = Compliance.from_seafloor_synthetic(self, max_freq)
        ncompl.write(base_name, out_dir=out_dir)

    def plot(self, fmin=0.001, fmax=0.1, fstep=0.001, outfile=None, show=True):
        """
        Plot the spectral representation of the noise sources in dB
        """
        f = np.arange(fmin, fmax + fstep / 2, fstep)
        # Plot
        fig, axs = plt.subplots(2, 1, sharex='col')
        # Plot the pressure signal
        for ch_code, color in zip(('IGP', 'NOP'), ('r', 'b')):
            axs[0].semilogx(f, self.source_by_code(ch_code).resample_values(f),
                            color, label=ch_code)
        axs[0].set_ylabel('dB ref 1 Pa^2/Hz')
        axs[0].legend()
        axs[0].set_ylim(-20, 60)
        axs[1].set_title('Pressure')
        # Plot the accel
        for ch_code, color in zip(
                ('IGZ', 'NOS', 'NTZ_max', 'NTZ_min', 'NTH_max', 'NTH_min'),
                ('r', 'b', 'g', 'g--', 'm', 'm--')):
            axs[1].semilogx(f, self.source_by_code(ch_code).resample_values(f),
                            color, label=ch_code)
        axs[1].set_ylabel('dB ref 1 (m/s^2)^2/Hz')
        axs[1].set_xlabel('Frequency (Hz)')
        axs[1].set_ylim(-200, -100)
        axs[1].legend()
        axs[1].set_title('Seismometer')
        plt.suptitle('SeafloorSynthetic components')
        if outfile is not None:
            plt.savefig(outfile)
        if show is True:
            plt.show()
        return axs

    def source_trace(self, code, trace_base, accel_to_vel=False, phases=None):
        """
        Return the :class:`obspy.stream.Trace` matching the given source code

        code (str): Valid source code
        trace_base (:class:`obspy.stream.Trace`): Trace to use as reference
            for dates, length, sampling rate, station, network and response
        accel_to_vel (bool): source PSDVals is an acceleration and should
            be converted to velocity.
        phases (np.array): list of phases to force fft to have (to correlate
            with another channel)
        """
        if accel_to_vel is False:
            return self.source_by_code(code).as_trace(trace_base, channel=code,
                                                      phases=phases)
        else:
            return self.source_by_code(code).accel_as_vel.as_trace(
                trace_base, channel=code, phases=phases)

    def streams(self, ref_trace, s_response=1., p_response=1.,
                network='XX', station='SSSSS', plotit=False, forceInt32=False):
        """
        Return streams generated from to the noise and signal levels
        Simply multiplies physical values by a sensitivity value, would be
        better to convolve with instrument response.

        Args:
            ref_trace (:class: `obspy.Trace`): trace with time base to use.
                band_code must be "L" and sampling rate near 1 sps.
            s_response (:class:`obspy.core.response.Response` or float):
                Seismometer response (counts/m/s).  If float, assumes flat
                response with this gain.
            p_response (class:`obspy.core.response.Response` or float):
                Pressure gauge response (counts/Pa).  If float, assumes flat
                response with this gain.
            network (str): Network code (1-2 characters)
            station (str): Station code (1-5 characters)
            forceInt32 (bool): force output data to have dtype=np.int32

        Returns:
            A tuple of (data, sources, inv), where
                data (:class:`obspy.Stream): synthetic seafloor BB 4C data
                sources (:class:`obspy.Stream`): individual noises and signals
                inv (:class:`obspy.core.Inventory`): channel metadata
        """
        # SETUP
        sr = ref_trace.stats.sampling_rate
        # Validate inputs
        if not ref_trace.stats.channel[0] == 'L':
            raise ValueError("ref_trace channel code ({}) doesn't start with L"
                             .format(ref_trace.stats.channel))
        if sr > 2 or sr < 0.5:
            raise ValueError(f'ref_trace {sr=} is not between 0.5 and 2 sps')
        # Set up variables
        trace_pts = ref_trace.stats.npts
        npts = 2**int(np.ceil(np.log2(trace_pts)))
        location = ref_trace.stats.location
        # channel = ref_trace.stats.channel
        _ = np.linspace(0, ref_trace.stats.sampling_rate / 2, npts)
        if not isinstance(p_response, Response):
            p_response = Response.from_paz([], [], p_response, 1.0, 'm/s',
                                           'count')
            # obspy doesn't understand Pa units, stuff them in afterwards
            p_response.instrument_sensitivity.input_units = 'Pa'
            p_response.response_stages[0].input_units = 'Pa'
        if not isinstance(s_response, Response):
            s_response = Response.from_paz([], [], s_response, 1.0, 'm/s',
                                           'count')

        # Prepare base seismo and pressure traces
        s_trace_base = ref_trace.copy()    # Don't overwrite original
        s_trace_base.stats.station = station
        s_trace_base.stats.network = network
        s_trace_base.stats.response = s_response
        p_trace_base = s_trace_base.copy()
        p_trace_base.stats.response = p_response

        # CREATE NOISE + IG/COMPLIANCE TRACES BY SOURCE
        sources = Stream([])

        # FOR THE PRESSURE CHANNEL
        # IG wave pressure signal
        IG_trace, IG_phases = self.source_trace("IGP", p_trace_base)
        sources += IG_trace
        sources += self.source_trace("NOP", p_trace_base)[0]

        # FOR THE SEISMOMETER CHANNELS
        # Vertical compliance signal
        # Phase_velocity = Phase_pressure + 270°
        sources += self.source_trace("IGZ", s_trace_base, True,
                                     phases=IG_phases-np.pi/2)[0]
        sources += self.source_trace("NOS", s_trace_base, True)[0]
        # Tilt noise model
        noise_max, _ = self.source_trace('NTH_max', s_trace_base, True)
        dyntilt_amp, dyntilt_angle = self.make_tilt_ts(noise_max)
        angfact = np.sin(np.radians(self.Z_offset_angles[0]))
        azefact_1 = np.sin(np.radians(self.Z_offset_angles[1]))
        azefact_2 = np.cos(np.radians(self.Z_offset_angles[1]))
        N_noise = dyntilt_amp.data * np.cos(np.radians(dyntilt_angle.data))
        E_noise = dyntilt_amp.data * np.sin(np.radians(dyntilt_angle.data))
        Z_noise = angfact * (azefact_1 * N_noise + azefact_2 * E_noise)
        sources += noise_max.copy()
        sources[-1].stats.channel = "NT1"
        sources[-1].data = N_noise
        sources += noise_max.copy()
        sources[-1].stats.channel = "NT2"
        sources[-1].data = E_noise
        sources += noise_max.copy()
        sources[-1].stats.channel = "NTZ"
        sources[-1].data = Z_noise

        if plotit is True:
            sources.plot(equal_scales=False)

        # CREATE SYNTHETIC BBOBS CHANNELS
        data = Stream([])
        for k, v in self.stream_source_codes.items():
            data += self._summed_channel(k, sources, v)
        if forceInt32 is True:
            for tr in data:
                tr.data = np.require(tr.data, dtype=np.int32)
        if plotit is True:
            data.plot(equal_scales=False)

        # Create Inventory
        channels = []
        for k, v in self.stream_source_codes.items():
            if k[1] == 'D':
                resp = p_response
                dip = 90.
            else:
                resp = s_response
                dip = 0.
                if k[2] == 'Z':
                    dip = -90.
            # Add BBOBS channels
            channels.append(Channel(k, location, 0, 0, 0, 0, response=resp,
                                    dip=dip))
            # Add source channels
            for x in v:
                channels.append(Channel(x, location, 0, 0, 0, 0, response=resp,
                                        dip=dip))
        stations = [Station(station, 0, 0, 0, channels=channels)]
        networks = [Network(network, stations=stations)]
        inv = Inventory(networks=networks)

        return data, sources, inv

    @staticmethod
    def _summed_channel(channel, source, source_chs):
        tr = source.select(channel=source_chs[0])[0].copy()
        tr.stats.channel = channel
        if len(source_chs) > 1:
            for c in source_chs[1:]:
                new_source = source.select(channel=c)[0]
                for key in ('station', 'network', 'location', 'response'):
                    assert new_source.stats[key] == tr.stats[key]
                tr.data += new_source.data
        return tr

    def make_tilt_ts(self, noise_max, coefficients=TideCoefficients(),
                     plotit=False):
        """
        make a simple tilt time series summing signals of given periods,
        amplitudes and starting phases

        Args:
            noise_max (:class:`obspy.core.Trace`): maximum tilt noise time series
            coefficients (TideCoefficients): the tidal coefficients
        """
        tide_trace = coefficients.make_trace(noise_max)
        tide_trace.stats.channel = 'TID'
        # normalize between (-self.noise_tilt_variance dB) and 1
        in_max = np.max(tide_trace.data)
        in_min = np.min(tide_trace.data)
        out_max = 1.
        out_min = 10**(-self.noise_tilt_variance / 20)
        tide_trace.data = (tide_trace.data - in_min)*(out_max-out_min)/(in_max-in_min) + out_min

        if plotit:
            tide_trace.plot()

        amp_trace = tide_trace.copy()
        amp_trace.stats.channel = 'AMP'
        amp_trace.data *= noise_max

        angles_trace = tide_trace.copy()
        angles_trace.stats.channel = 'ANG'
        angle_range = abs(self.noise_tilt_direction_limits[1]
                          - self.noise_tilt_direction_limits[0])
        angles_trace.data = ((angles_trace.data * angle_range)
                             + min(self.noise_tilt_direction_limits))

        if plotit:
            Stream([tide_trace, amp_trace, angles_trace]).plot(equal_scale=False)

        return amp_trace, angles_trace


if __name__ == "__main__":
    # Show an example
    wdepth = 2000
    noise_model = SeafloorSynthetic(wdepth)
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
    sd_data.plot()
    sd_data.plot_coherences()
