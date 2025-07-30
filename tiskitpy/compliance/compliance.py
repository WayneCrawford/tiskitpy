"""
Python compliance class

Authors:  W. Crawford, A. Doran
"""
from pathlib import Path

# import math
import numpy as np
from matplotlib import pyplot as plt

# from ..response_functions import ResponseFunctions


class Compliance(object):
    """
    Seafloor compliance class

    Attributes:
        freqs (:class:`numpy.ndarray`): Frequencies (Hz)
        values (:class:`numpy.ndarray`): Normalized compliance values (1/Pa)
        uncertainties (:class:`numpy.ndarray`): Normalized compliance
            uncertainties (1/Pa)
        water_depth (float): water depth in meters
        noise_channel (str or None): If a str, compliance comes from data and
            this is the channel on which noise was assumed to dominate
        gravity_corrected (bool): Has data-estimated compliance been corrected
            for gravitational attraction terms?
    """
    def __init__(self, freqs, values, uncertainties, water_depth,
                 noise_channel, gravity_corrected=False):
        """
        Seafloor compliance data class

        Args:
            freqs (:class:`numpy.ndarray`): Frequencies (Hz)
            values (:class:`numpy.ndarray`): Normalized compliance values
                (1/Pa)
            uncertainties (:class:`numpy.ndarray`): Normalized compliance
                uncertainties (1/Pa)
            water_depth (float): water depth in meters
            noise_channel (str): (str or None): If a str, compliance comes
                from data and this is the channel on which noise was assumed
                to dominate
            gravity_corrected (bool): Has data-estimated compliance been
                corrected for gravitational attraction terms?
        """
        self.freqs = freqs
        self.values = values
        self.uncertainties = uncertainties
        self.water_depth = water_depth
        self.noise_channel = noise_channel
        self.gravity_corrected = gravity_corrected

    # def frf_to_compliance(xf, wdepth, z_units='M/S'):
    @classmethod
    def from_response_functions(cls, rfs, wdepth, max_freq=None, z_str='*Z'):
        """
        Extracts compliance from ResponseFunctions object

        Args:
            rfs (:class:`tiskitpy.ResponseFunctions`): z/p transfer
                function(s). The input_channel should be the pressure channel
            wdepth (float): water depth (m)
            max_freq (float): maximum frequency to save.  If None, use
                sqrt(g/(2*pi*wdepth))
            z_str(str): channel_id to use for z channel (may include '*'
                wildcard)
        """
        # Validate input fields
        if not rfs.input_units.upper() == 'PA':
            raise ValueError(f'{rfs.input_units.upper()=}, not "PA"')
        try:
            _ = rfs.value(z_str)
        except Exception:
            raise ValueError(f'output channel "{z_str}" not in '
                             f'{rfs.output_channel_ids=}')
        if max_freq is None:
            max_freq = Compliance.max_freq(wdepth)  # about one wavelength

        f = rfs.freqs[rfs.freqs <= max_freq]
        zp = rfs.value(z_str)[rfs.freqs <= max_freq]
        zp_uncert = rfs.uncertainty(z_str)[rfs.freqs <= max_freq]
        z_units = rfs.output_units(z_str)
        return cls(f, Compliance._zp_to_ncompl(f, zp, wdepth, z_units),
                   Compliance._zp_to_ncompl(f, zp_uncert, wdepth, z_units),
                   wdepth, rfs.noise_channel)

    @classmethod
    def from_file(cls, filename):
        "NOT YET IMPLEMENTED"
        return

    @classmethod
    def from_earth_model_1D(cls, water_depth, freqs, earth_model,
                            limit_freqs=True):
        """
        Create object with the compliance of a 1D earth model

        Args:
            water_depth (float): water depth (m)
            freqs (list or :class:`numpy.ndarray`): frequencies
            earth_model (:class:`tiskitpy.compliance.EarthModel1D`): 1D earth
                model
            limit_freqs (bool): limit frequencies to below
                Compliance.max_freq()?
        """
        if isinstance(freqs, list):
            freqs = np.array(freqs)
        if limit_freqs is True:
            freqs = freqs[freqs < Compliance.max_freq(water_depth)]

        ncompl = Compliance.calc_norm_compliance(water_depth, freqs,
                                                 earth_model)
        uncert = np.zeros(ncompl.shape)
        return cls(freqs, ncompl, uncert, water_depth, None, True)

    @classmethod
    def from_seafloor_synthetic(cls, obj, max_freq=True):
        """
        Create object with the compliance of a SeafloorSynthetic object

        Uses objs earth_model, IG_Pa_seafloor.freqs and water_depth attributes

        Args:
            obj (:class:`tiskitpy.SeafloorSynthetic`): the object
            max_freq (bool or float): limit frequencies to:
                True: below Compliance.max_freq()
                float: below the value
                False: no limit applied
        """
        freqs = obj.IG_Pa_seafloor.freqs
        if not isinstance(max_freq, bool):
            freqs = freqs[freqs < max_freq]
            max_freq = False

        return cls.from_earth_model_1D(obj.water_depth, freqs, obj.earth_model,
                                       max_freq)

    @staticmethod
    def max_freq(water_depth):
        """
        Return estimated maximum compliance frequency for the given water depth

        About one wavelength, from Janiszewski et al. (2019?)
        """
        return np.sqrt(9.8/(2*np.pi*water_depth))

    def __str__(self):
        s = f"{self.__class__.__name__} object:\n"
        s +=  "  {} frequencies, from {} to {} Hz\n".format(
            len(self.freqs), np.min(self.freqs), np.max(self.freqs))
        s += f"  water_depth='{self.water_depth}'\n"
        s += f"  noise_channel={self.noise_channel}\n"
        s += f"  gravity_corrected={self.gravity_corrected}"
        return s

    def correct_gravity_terms(self):
        "NOT YET IMPLEMENTED"
        return

    def write(self, base_name, units='1/Pa', out_dir=None):
        """
        Save compliance to a CSV file

        Args:
            base_name (str): base filename.  "_{units}.csv" will be appended
            units (str): units in which to save compliance.  One of
                '1/Pa', 'm/Pa', 'm/s/Pa', 'm/s^2/Pa'
                ('1/Pa' is normalized compliance)
            out_dir (str or :class:`Path`): output directory (None: save to
                working directory)
        """
        if units == '1/Pa':
            filename = f'{base_name}_compliance_Pa-1.csv'
        elif units == 'm/Pa':
            filename = f'{base_name}_compliance_m.Pa-1.csv'
        elif units == 'm/s/Pa':
            filename = f'{base_name}_compliance_m.s-1.Pa-1.csv'
        elif units == 'm/s^2/Pa':
            filename = f'{base_name}_compliance_m.s-2.Pa-1.csv'
        else:
            raise ValueError(f"{base_name=} is not in ('1/Pa', 'm/Pa', "
                             "'m/s/Pa', 'm/s^2/Pa')")
        v, u = self._convert_compliance(units)
        if out_dir is not None:
            filename = str(Path(out_dir) / filename)
        with open(filename, "w") as fid:
            fid.write(f'# units={units}\n')
            fid.write(f'# water_depth={self.water_depth}\n')
            fid.write(f'# noise_channel={self.noise_channel}\n')
            fid.write(f'# gravity_corrected={self.gravity_corrected}\n')
            fid.write('frequencies;compliance;uncertainty;phase\n')
            for freq, ncompl, uncert in zip(self.freqs, self.values,
                                            self.uncertainties):
                fid.write('{:.5g};{:.5g};{:.5g};{:.5g}\n'
                          .format(freq, np.abs(ncompl), np.abs(uncert),
                                  np.angle(ncompl, deg=True)))

    def write_counts(self, base_name, z_response, p_response, out_dir=None):
        """
        Save compliance IN COUNTS to a CSV file

        Only useful for creating compliance values for people who don't
        know how to use inventories

        Args:
            base_name (str): base filename.  "_{units}.csv" will be appended
            z_reponse (:class:`obspy.core.inventory.Response`): response of
                the z channel
            p_reponse (:class:`obspy.core.inventory.Response`): response of
                the p channel
            out_dir (str or :class:`Path`): output directory (None: save to
                working directory)
        """
        filename = f'{base_name}_compliance_COUNTS.csv'

        assert p_response.instrument_sensitivity.input_units.lower() == 'pa'
        if z_response.instrument_sensitivity.input_units.lower() == 'm/s':
            v, u = self._convert_compliance('m/s/Pa')
        elif z_response.instrument_sensitivity.input_units.lower() == 'm/s^2':
            v, u = self._convert_compliance('m/s^2/Pa')
        elif z_response.instrument_sensitivity.input_units.lower() == 'm':
            v, u = self._convert_compliance('m/Pa')
        z_resp_f = z_response.get_evalresp_response_for_frequencies(self.freqs)
        p_resp_f = p_response.get_evalresp_response_for_frequencies(self.freqs)
        # responses are in counts/physical units, so multiply by z_resp/p_resp
        v = v * z_resp_f/p_resp_f
        u = u * z_resp_f/p_resp_f
        if out_dir is not None:
            filename = str(Path(out_dir) / filename)
        with open(filename, "w") as fid:
            fid.write('# units=COUNTS/COUNTS\n')
            fid.write(f'# water_depth={self.water_depth}\n')
            fid.write(f'# noise_channel={self.noise_channel}\n')
            fid.write(f'# gravity_corrected={self.gravity_corrected}\n')
            fid.write('frequencies;compliance;uncertainty;phase\n')
            for freq, ncompl, uncert in zip(self.freqs, v, u):
                fid.write('{:.5g};{:.5g};{:.5g};{:.5g}\n'.format(
                          freq, np.abs(ncompl), np.abs(uncert),
                          np.angle(ncompl, deg=True)))

    def _convert_compliance(self,  units):
        """
        Convert object's compliance and uncertainty values to the given units

        Args:
            units (str): units to convert to.  One of '1/Pa', 'm/Pa', 'm/s/Pa',
                or 'm/s^2/Pa' ('1/Pa' changes nothing)

        Returns:
            tuple:
                :class:`numpy.ndarray`: converted compliances
                :class:`numpy.ndarray`: converted uncertainties
        """
        if units == '1/Pa':
            return self.values, self.uncertainties
        else:
            k = Compliance.gravd(2*np.pi*self.freqs, self.water_depth)
            if units == 'm/Pa':
                multiplier = k**(-1)
            elif units == 'm/s/Pa':
                multiplier = 2 * np.pi * self.freqs / k
            elif units == 'm/s^2/Pa':
                multiplier = (2 * np.pi * self.freqs)**2 / k
            else:
                raise ValueError(f"{units=} is not in ('1/Pa', 'm/Pa', "
                                 "'m/s/Pa', 'm/s^2/Pa')")
            return self.values * multiplier, self.uncertainties * multiplier

    def plot(self, errorbars=True, show=True, outfile=None):
        """
        Plot the compliance

        Args:
            errorbars (bool): plot error bars
            show (bool): show on the screen
            outfile (str): save figure to this filename
        Returns:
            axa, axp: axis pair amplitude, phase
        """
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(3, 1, hspace=0)
        ax_a = fig.add_subplot(gs[:2, :])
        ax_p = fig.add_subplot(gs[2, :])
        # fig, axs = plt.subplots(2, 1, sharex=True)
        ncompl = self.values.copy()
        nuncert = self.uncertainties.copy()
        ibad = (ncompl == 0).nonzero()

        # Plot amplitude
        fig.suptitle("Compliance")
        ncompl[ncompl == 0] = np.nan
        nuncert[nuncert == 0] = np.nan
        if errorbars is True:
            ax_a.errorbar(self.freqs, np.abs(ncompl), np.abs(nuncert),
                          fmt='b_', ecolor='k', markersize=3, label='Estimated')
            if np.any(ncompl is not np.nan):
                ax_a.set_yscale('log')
            ax_a.set_xscale('log')
        else:
            ax_a.loglog(self.freqs, np.abs(ncompl + nuncert), color="blue",
                        linewidth=0.5)
            ax_a.loglog(self.freqs, np.abs(ncompl - nuncert), color="blue",
                        linewidth=0.5)
            ax_a.loglog(self.freqs, np.abs(ncompl), color="black", label='Estimated')
        # ax_a.set_xlim(self.freqs[1], self.freqs[-1])
        ax_a.set_ylabel("Norm Compliance (1/Pa)")
        ax_a.tick_params('x', which='both', direction="in")

        # Plot phase
        phases = np.angle(ncompl, deg=True)
        phases[ibad] = np.nan
        phase_lim = 220
        igood = np.invert(np.isnan(phases))
        wrap_phases = phases.copy()
        wrap_phases[igood] = np.unwrap(phases[igood], period=360)
        if (np.all(wrap_phases[igood] > -phase_lim)
                and np.all(wrap_phases[igood] < phase_lim)):
            ax_p.semilogx(self.freqs, wrap_phases)
        else:
            ax_p.semilogx(self.freqs, phases)
        ax_p.set_ylim(-phase_lim, phase_lim)
        # ax_p.set_xlim(self.freqs[1], self.freqs[-1])
        ax_p.set_yticks((-180, 0, 180))
        ax_p.set_ylabel("Phase")
        ax_p.set_xlabel("Frequency (Hz)")

        # Show and/or save plot
        # fig.tight_layout()
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()

        return ax_a, ax_p

    @staticmethod
    def _zp_to_ncompl(freqs, zp, wdepth, z_units):
        """
        Calculate compliance from the z/p ratio, freqs and water depth

        normalized compliance is defined as k*Z/P, with k in 1/m, Z in m and
        P in Pa. Its units are 1/Pa

        Args:
            freqs (:class:`numpy.nparray`): frequencies (1/s)
            zp (:class:`numpy.nparray`): vertical motion / pressure.  Pressure
                units are Pa, z_units are specified by the paramter z_units
            wdepth (float): water depth (m)
            z_units (str): z units, one of 'M', 'M/S' or 'M/S^2'
        """
        omega = 2 * np.pi * freqs
        k = Compliance.gravd(omega, wdepth)
        if z_units.upper() == 'M':
            omega_term = np.ones(omega.shape)
        elif z_units.upper() == 'M/S':
            omega_term = omega**(-1)
        elif z_units.upper() == 'M/S^2':
            omega_term = omega**(-2)
        else:
            raise ValueError(f'Z_units ({z_units}) is not in ("M", "M/S", "M/S^2")')
        return zp * k * omega_term

    @staticmethod
    def gravd(W, h):
        """
        Return linear ocean surface gravity wave wavenumbers

        Args:
            W (:class:`numpy.ndarray`): angular frequencies (rad/s)
            h (float): water depth (m)

        Returns:
            K (:class:`numpy.ndarray`): wavenumbers (rad/m)

        >>> f = np.array([0.0001, 0.001, 0.01, 0.1, 1, 10])
        >>> K = gravd(2 * np.pi * f, 2000)
        >>> wlen = 2 * np.pi * np.power(K, -1)
        >>> np.set_printoptions(precision=1)
        >>> print(f'{K=} rad/m')
        K=array([4.5e-06, 4.5e-05, 5.2e-04, 4.0e-02, 4.0e+00, 4.0e+02]) rad/m
        >>> print(f'{wlen=} m')
        wlen=array([1.4e+06, 1.4e+05, 1.2e+04, 1.6e+02, 1.6e+00, 1.6e-02]) m
        >>> print(f'c={f*wlen} m/s')
        c=[140.  139.8 121.1  15.6   1.6   0.2] m/s
        """
        # W must be array
        if not isinstance(W, np.ndarray):
            W = np.array([W])
        if np.any(W < 0):
            raise ValueError('there are omegas <= 0')
        G = 9.81  # 9.78 et equator, 9.83 at poles
        # N = len(W)
        W2 = W*W
        kDEEP = W2/G
        kSHAL = W/(np.sqrt(G*h))
        erDEEP = np.ones(np.shape(W)) - G*kDEEP*_dtanh(kDEEP*h)/W2
        one = np.ones(np.shape(W))
        d = np.copy(one)
        done = np.zeros(np.shape(W))
        done[W == 0] = 1   # if W==0, k is also zero
        nd = np.where(done == 0)

        k1 = np.copy(kDEEP)
        k2 = np.copy(kSHAL)
        e1 = np.copy(erDEEP)
        ktemp = np.copy(done)
        e2 = np.copy(done)
        e2[W == 0] = 0

        while True:
            e2[nd] = one[nd] - G*k2[nd] * _dtanh(k2[nd]*h)/W2[nd]
            d = e2*e2
            done = d < 1e-20
            if done.all():
                K = k2
                break
            nd = np.where(done == 0)
            ktemp[nd] = k1[nd]-e1[nd]*(k2[nd]-k1[nd])/(e2[nd]-e1[nd])
            k1[nd] = k2[nd]
            k2[nd] = ktemp[nd]
            e1[nd] = e2[nd]
        return K

    @staticmethod
    def raydep(P, om, d, ro, vp2, vs2):
        """
        Propagator matrix solutionn for P-SV waves, minor vector method

        Args:
            P (float): surface wave slowness (s/m)
            om (float): surface wave angular frequency (radians/sec)
            d (:class:`numpy.ndarray`): layer thicknesses (meters?)
            rho (:class:`numpy.ndarray`): layer densities (kg/m^3)
            vp2 (:class:`numpy.ndarray`): layer P velocities squared (m/s)^2
            vs2 (:class:`numpy.ndarray`): layer shear velocities squared
                (m/s)^2

        Returns:
            (list): Parameters, each value is at layer top
                v (:class:`numpy.ndarray`): vertical velocity (m/s?)
                u (:class:`numpy.ndarray`): horizontal velocity (m/s?)
                zz (:class:`numpy.ndarray`): vertical stress (Pa?)
                zx (:class:`numpy.ndarray`): horizontal stress (Pa?)

        Notes:
            d, rho, vp2 and vs2 have one value for each layer (top to bottom),
                must be same length
            (Normalized compliance = -k*v/(omega*sigzz) )

        >>> P = 1/140    # Corresponds to 2000m depth, low freqs
        >>> om = 2 * np.pi * 0.005
        >>> d = np.array([1000, 1000, 1000, 3000, 3000])
        >>> rho = np.array([3000, 3000, 3000, 3000, 3000])
        >>> vp2 = np.array([3000**2, 4000**2, 5000**2, 7500**2, 8200**2])
        >>> vs2 = np.array([1600**2, 2300**2, 2800**2, 4300**2, 4700**2])
        >>> np.set_printoptions(precision=1)
        >>> raydep(P, om, d, rho, vp2, vs2)
        (array([1. , 0.7, 0.5, 0.4, 0.3]), array([ 1.8e-01,  6.0e-02,  1.2e-02,  5.2e-05, -5.8e-02]), array([-2.8e+08, -2.8e+08, -2.7e+08, -2.6e+08, -1.9e+08]), array([-0.0e+00,  3.1e+07,  5.4e+07,  7.6e+07,  1.0e+08]))
        """
        mu = ro * vs2
        n = len(d)
        ist = n-1
        # ysav = 0
        psq = P*P
        r2 = 2 * mu[ist] * P
        # R and S are the "Wavenumbers" of compress and shear waves in botlayer
        # RoW and SoW are divided by ang freq
        RoW = np.sqrt(psq - 1/vp2[ist])
        SoW = np.sqrt(psq - 1/vs2[ist])
        ym = np.zeros((ist+1, 5))
        i = ist
        y = np.zeros((5, ))     # Minor vector matrix
        # Stress-displacement vector: (vert vel, hor vel, vert stress, hor stress)
        x = np.zeros((i+1, 4))

        y[2] =  RoW
        y[3] = -SoW
        y[0] = (RoW*SoW - psq) / ro[i]
        y[1] = r2*y[0] + P
        y[4] = ro[i] - r2*(P + y[1])
        ym[i, :] = y
        # *****PROPAGATE UP LAYERS*********
        while i > 0:
            i = i-1
            ha = psq - 1/vp2[i]
            ca, sa = _argdtray(om*d[i], ha)
            hb = psq - 1/vs2[i]
            cb, sb = _argdtray(om*d[i], hb)
            hbs = hb*sb
            has = ha*sa
            r1 = 1 / ro[i]
            r2 = 2 * mu[i] * P
            b1 = r2*y[0] - y[1]
            g3 = (y[4] + r2*(y[1]-b1)) * r1
            g1 = b1 + P*g3
            g2 = ro[i]*y[0] - P*(g1+b1)
            e1 = cb*g2 - hbs*y[2]
            e2 = -sb*g2 + cb*y[2]
            e3 = cb*y[3] + hbs*g3
            e4 = sb*y[3] + cb*g3
            y[2] = ca*e2 - has*e4
            y[3] = sa*e1 + ca*e3
            g3 = ca*e4 - sa*e2
            b1 = g1 - P*g3
            y[0] = (ca*e1 + has*e3 + P*(g1+b1))*r1
            y[1] = r2*y[0] - b1
            y[4] = ro[i]*g3 - r2*(y[1] - b1)
            ym[i, :] = y

        # de = y[4]/np.sqrt(y[0]*y[0] + y[1]*y[1])
        ynorm = 1/y[2]
        y[0: 4] = np.array([0, -ynorm,  0,  0])
        # *****PROPAGATE BACK DOWN LAYERS*********
        while i <= ist:
            x[i, 0] = -ym[i, 1]*y[0] - ym[i, 2]*y[1] + ym[i, 0]*y[3]
            x[i, 1] = -ym[i, 3]*y[0] + ym[i, 1]*y[1] - ym[i, 0]*y[2]
            x[i, 2] = -ym[i, 4]*y[1] - ym[i, 1]*y[2] - ym[i, 3]*y[3]
            x[i, 3] =  ym[i, 4]*y[0] - ym[i, 2]*y[2] + ym[i, 1]*y[3]
            ls = i
            if i >= 1:
                sum = abs(x[i, 0] + i*x[i, 1])
                # pbsq = 1 / vs2[i]
                if sum < 1e-4:
                    break

            ha = psq - 1/vp2[i]
            ca, sa = _argdtray(om*d[i], ha)
            hb = psq-1/vs2[i]
            cb, sb = _argdtray(om*d[i], hb)
            hbs = hb*sb
            has = ha*sa
            r2 = 2*P*mu[i]
            e2 = r2*y[1] - y[2]
            e3 = ro[i]*y[1] - P*e2
            e4 = r2*y[0] - y[3]
            e1 = ro[i]*y[0] - P*e4
            e6 = ca*e2 - sa*e1
            e8 = cb*e4 - sb*e3
            y[0] = (ca*e1 - has*e2+P*e8) / ro[i]
            y[1] = (cb*e3 - hbs*e4+P*e6) / ro[i]
            y[2] = r2*y[1] - e6
            y[3] = r2*y[0] - e8
            i = i+1
        #
        # if x(1,3) == 0
        #   error('vertical surface stress = 0 in DETRAY');
        # end
        ist = ls

        return x[:, 0], x[:, 1], x[:, 2], x[:, 3]


#     @staticmethod
#     def plot_compliance_stack(psd, zstr, pstr, water_depth,
#                               seawater_density=1030, show=True,
#                               outfile=None):
#         """
#         Plot from top to bottom: Z PSD, P PSD, coherence, Z/P
#
#         Args:
#             psd (SpectralDensity): PSDs including Z and P
#             zstr (str): channel id sub/string matching the Z channel (see
#                 :meth:`SpectralDensity.channel_id() documentation)
#             pstr (str): channel id sub/string matching the P channel
#             water_depth (float): water depth in meters
#             seawater_density (float): average water density overhead (kg/m^3)
#             show (bool): show the result on the screen
#             outfile (str): save the plot to the named file
#         """
#         # Validate inputs
#         try:
#             _ = psd.channel_id(zstr)
#         except Exception:
#             raise ValueError(f'{zstr=} invalid/unique channel id for {psd=}')
#         try:
#             _ = psd.channel_id(pstr)
#         except Exception:
#             raise ValueError(f'{pstr=} invalid/unique channel id for {psd=}')
#         assert water_depth > 0, f'{water_depth=} is not greater than 0'
#         for id, units in zip((zstr, pstr), ('m/s^2', 'Pa')):
#             assert psd.channel_units(id) == units, f'{id} units are {psd.channel_units(id)},  not "{units}"'
#
#         fig, axs = plt.subplots(4, 1, sharex=True)
#         fig.subplots_adjust(hspace=0)
#         # Plot Z PSD
#         axs[0].semilogx(psd.freqs, 20*np.log10(np.abs(psd.autospect(zstr))))
#         axs[0].set_ylabel(r'Z (dB ref 1 $m/s^2/\sqrt{Hz}$)', fontsize='small')
#         axs[0].set_title('Compliance Stack')
#         # Plot Pressure PSD
#         axs[1].semilogx(psd.freqs, 10*np.log10(np.abs(psd.autospect(pstr))))
#         axs[1].set_ylabel(r'P (dB ref 1 $Pa/\sqrt{Hz}$)', fontsize='small')
#         # Plot Pressure-Z coherence
#         psd.plot_one_coherence(pstr, zstr, fig=fig, ax_a=axs[2],
#                                show_phase=False, ylabel='')
#         axs[2].set_ylabel('Coherence', fontsize='small')
#         axs[2].set_ylim([0.001, 1])
#         # Plot Z/P ratios
#         for noise_channel, label, color in zip(('output', 'input'),
#                                                ('z_noise', 'p_noise'),
#                                                ('r', 'b')):
#             frf = ResponseFunctions(psd, pstr, [zstr],
#                                     noise_channel=noise_channel)
#             igood = np.abs(frf.uncertainty(zstr)) < np.abs(frf.value(zstr))
#             # igood = np.abs(frf.uncertainty(zstr)) == np.abs(frf.uncertainty(zstr))
#             # axs[3].plot(frf.freqs[igood], np.abs(frf.value(zstr))[igood],
#                           c=color, label=label)
#             axs[3].errorbar(frf.freqs[igood],
#                             y=np.abs(frf.value(zstr))[igood],
#                             yerr=np.abs(frf.uncertainty(zstr))[igood],
#                             fmt='.', ms=1, c=color, label=label)
#         # Overlay theoretical Z/P relation for LF Rayleigh waves (seafloor
#         # moving water column): m/s^2/Pa = 1/rho*H
#         axs[3].axhline(1/(seawater_density*water_depth), c='k', ls='--')
#         axs[3].text(frf.freqs[0], 1/(seawater_density*water_depth), r'1/$\rho H$',
#                     verticalalignment='bottom')
#         axs[3].set_ylabel(f'Z/P ({frf.output_units(zstr)}/{frf.input_units})',
#                           fontsize='small')
#         axs[3].set_xlabel('Frequency (Hz)')
#         axs[3].set_yscale('log')
#         # Put predicted compliance max frequency vertical line on each plot
#         IG_fmax = np.sqrt(9.8/(2*np.pi*water_depth)) # about one wavelength
#         for ax in axs:
#             ax.axvline(IG_fmax, c='k', ls='--')
#         axs[0].text(IG_fmax, np.max(20*np.log10(np.abs(psd.autospect(zstr)))),
#                    r'$\sqrt{\frac{g}{2 \pi H}}$', rotation='vertical',
#                    horizontalalignment='right', verticalalignment='top')
#
#         # Show or save plot
#         if show is False and outfile is None:
#             raise ValueError('Plot neither shown nor saved!')
#         if show is True:
#             plt.show()
#         if outfile is not None:
#             plot.savefig(outfile)

    @staticmethod
    def calc_compliance(wdepth, freq, model):
        """
        Return compliance of a 1D model

        Args:
            wdepth (float): water depth (m)
            freq (:class:`numpy.nparray`): frequencies (1/s)
            model (:class:`EarthModel1D`): 1D earth model
        """
        if np.any(freq <= 0):
            raise ValueError('At least one freq <= 0: cannot calculate compliance')
        vpsq = model.vps * model.vps
        vssq = model.vss * model.vss
        omega = 2 * np.pi * freq
        k = Compliance.gravd(omega, wdepth)
        ps = k / omega

        compl = np.zeros((len(ps)))
        for i in np.arange((len(ps))):
            v, _, sigzz, _ = Compliance.raydep(ps[i], omega[i], model.thicks,
                                               model.rhos, vpsq, vssq)
            # If raydep returned complex values, would need to divide by a
            # further 1j to go from (m/s)/Pa to m/Pa.  Returned value should
            # be negative because seafloor is lowest (DOWN) under maxixum
            # pressure, for quasi-static
            compl[i] = v[0] / (omega[i] * sigzz[0])
        return compl

    @staticmethod
    def calc_norm_compliance(wdepth, freq, model):
        """
        Return normalized compliance of a 1D model

        norm compliance == k(omega) * compliance

        Args:
            wdepth (float): water depth (m)
            freq (:class:`numpy.nparray`): frequencies (1/s)
            model (:class:`EarthModel1D`): 1D earth model

        >>> wdepth = 2000
        >>> freqs = np.array([0.001, 0.003, 0.005, 0.01, 0.03])
        >>> model = EarthModel1D([[1000, 3000, 3000, 1600],
        ...                      [1000, 3000, 4000, 2300],
        ...                      [1000, 3000, 5000, 2800],
        ...                      [3000, 3000, 7500, 4300],
        ...                      [3000, 3000, 8200, 4700]])
        >>> np.set_printoptions(precision=1)
        >>> calc_norm_compliance(wdepth, freqs, model)
        array([-1.4e-11, -2.0e-11, -2.6e-11, -4.2e-11, -9.0e-11])
        """
        k = Compliance.gravd(2 * np.pi * freq, wdepth)
        return k * Compliance.calc_compliance(wdepth, freq, model)


def _dtanh(x):
    """
    Stable hyperbolic tangent

    Args:
        x (:class:`numpy.ndarray`)
    """
    a = np.exp(x*(x <= 50))
    one = np.ones(np.shape(x))

    y = (abs(x) > 50) * (abs(x)/x) + (abs(x) <= 50)*((a-one/a) / (a+one/a))
    return y


def _argdtray(wd, h):
    hh = np.sqrt(abs(h))    # magnitude of wavenumber/freq
    th = wd * hh            # # of waves (or e-foldings) in layer (radians)
    if th >= 1.5e-14:
        if h <= 0:          # propagating wave
            c =  np.cos(th)
            s = -np.sin(th) / hh
        else:               # evenescent wave
            d = np.exp(th)
            c =  0.5*(d + 1/d)
            s = -0.5*(d - 1/d)/hh
    else:
        c = 1
        s = -wd
    return c, s


# Now that we have Compliance.from_response_functions, do we need this?
# def frf_to_compliance(xf, wdepth, z_units='M/S'):
#     """
#     Changes the response for each out_channel from z_units/Pa
#     to 1/Pa (normalized compliance)
#
#     Args:
#         xf (:class:`tiskitpy.ResponseFunctions`): z/p transfer function(s)
#         wdepth (float): water depth (m)
#         z_units (str): z units, one of 'M', 'M/S' or 'M/S^2'
#     """
#     compl = deepcopy(xf)
#     for oc in compl.output_channels:
#         if not compl.output_units(oc).upper() == z_units:
#             raise ValueError('output_units({}) ({}) != "{}"'.format(
#                 oc, compl.output_units(oc), z_units))
#         if not compl.input_units.upper() == 'PA':
#             raise ValueError(f'input_units ({compl.input_units}) != "PA"')
#         orig_resp = compl.response(oc)
#         new_resp = Compliance._zp_to_ncompl(compl.freqs, orig_resp,
#                                          wdepth, z_units)
#         # print(f'{new_resp/orig_resp=}')
#         # print(f'BEFORE {compl.response(oc)=}')
#         compl.put_response(new_resp, oc)
#         # print(f'AFTER {compl.response(oc)=}')
#         # compl._ds["response"].loc[dict(input=compl.input_channel,
#         #                             output=oc)] = new_resp
#     return compl


if __name__ == "__main__":
    import doctest
    doctest.testmod()
