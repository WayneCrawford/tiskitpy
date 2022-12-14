# Copyright 2021 Wayne Crawford
import pickle
import fnmatch  # Allows Unix filename pattern matching

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

from .spectral_density.utils import coherence_significance_level
from .spectral_density import SpectralDensity

np.seterr(all="ignore")
# np.set_printoptions(threshold=sys.maxsize)


class TransferFunctions(object):
    """
    Class of Transfer functions for a given input channel.

    From Bendat & Piersol, chapter 6.  The frequency response function is
    the relation between coherent parts of the signal: if the measured
    input x(t) = u(t) + m(t) and the measured output y(t) = v(t) + n(t),
    where u(t) and v(t) are coherent and m(t) and n(t) are not, then the
    frequency response function H(f) is such that v(t) = H(f)*u(t).  As
    to spectra, G_vv(f) = abs(H(f))^2 * G_uu(f)

    Args:
        sdm (:class:`.SpectralDensity`): Spectral density matrix objet
        in_chan (str): input channel.  Can use Unix wildcards (*, ?) but
            will return error if more than one string matches
        out_chans (list of str): output channels  (None => all but
            in_chan))
        noise_chan (str): 'input', 'output', 'equal', 'unknown'
        n_to_reject (int): number of neighboring frequencies for which the
            coherence must be above the 95% significance level in order
            to calculate transfer function (other values are set to 0,
            n_to_reject=0 means use all frequencies)
        min_freq (float or None): Return zero for frequencies below
            this value
        max_freq (float or None): Return zero for frequencies above
            this value
        quiet (bool): don't warn if creating a test object
        show_progress (bool): plot transfer functions and coherences
    """
    def __init__(self, sdm, in_chan, out_chans=None, noise_chan="output",
                 n_to_reject=3, min_freq=None, max_freq=None,
                 quiet=False, show_progress=False):
        """
        Attributes:
            _ds (:class: XArray.dataset): container for transfer functions and
                attributes
        """
        if sdm is None:
            if not quiet:
                print("Creating test TransferFunction object (no data)")
            if out_chans is None:
                raise ValueError("cannot have empty out_chans for test object")
            f = np.array([1])
            dims = ("input", "output", "f")
            shape = (1, len(out_chans), len(f))
            self._ds = xr.Dataset(
                data_vars=dict(value=(dims, np.zeros(shape, dtype="complex"))),
                coords=dict(input=[in_chan], output=out_chans, f=f),
            )

        else:
            if not isinstance(sdm, SpectralDensity):
                raise TypeError("sdm is not a SpectralDensity object")
            in_chan, out_chans = self._check_chans(sdm, in_chan, out_chans)
            in_units = sdm.channel_units(in_chan)
            out_units = [sdm.channel_units(x) for x in out_chans]
            f = sdm.freqs
            # Set properties
            dims = ("input", "output", "f")
            shape = (1, len(out_chans), len(f))
            # tf units are out_chan_units/in_chan_units
            # response units are in_chan_units/out_chan_units
            # tf * response gives tf w.r.t. data counts
            self._ds = xr.Dataset(
                data_vars={
                    "value": (dims, np.zeros(shape, dtype="complex")),
                    "uncert_mult": (dims, np.zeros(shape)),
                    "corr_mult": (dims, np.zeros(shape)),
                    "response": (dims, np.ones(shape, dtype="complex")),
                },
                coords= {
                    "input": [in_chan],
                    "output": out_chans,
                    "in_units": ("input", [in_units]),
                    "out_units": ("output", out_units),
                    "noise_chan": ("output", [noise_chan for x in out_chans]),
                    "f": f,
                },
                attrs=dict(n_winds=sdm.n_windows),
            )
            for out_chan in out_chans:
                xf, err_mult, corr_mult = self._calcxf(
                    sdm, in_chan, out_chan, noise_chan,
                    n_to_reject, min_freq, max_freq)
                self._ds["value"].loc[dict(input=in_chan,
                                           output=out_chan)] = xf
                self._ds["uncert_mult"].loc[dict(input=in_chan,
                                                 output=out_chan)] = err_mult
                self._ds["corr_mult"].loc[dict(input=in_chan,
                                               output=out_chan)] = corr_mult
                self._ds["response"].loc[
                    dict(input=in_chan, output=out_chan)
                ] = sdm.channel_response(out_chan) / sdm.channel_response(
                    in_chan
                )
            if show_progress:
                self._plot_progress(sdm)

    def __repr__(self):
        s = ("TransferFunctions(<SpectralDensity object>, "
             f"'{self.input_channel}', "
             f"{self.output_channels}, "
             f"{self.noise_channels[0]})")
        return s

    def __str__(self):
        s = "TransferFunctions object:\n"
        s += f"\tinput_channel='{self.input_channel}'\n"
        s += f"\toutput_channels={self.output_channels}\n"
        s += f"\tnoise_channels={self.noise_channels}\n"
        s += f"\tn_windows={self.n_windows}"
        return s

    @property
    def freqs(self):
        """Transfer function frequencies"""
        return self._ds.coords["f"].values

    @property
    def input_channel(self):
        """Transfer function input channel"""
        return str(self._ds.coords["input"].values[0])

    @property
    def output_channels(self):
        """Transfer function output channels"""
        return list(self._ds.coords["output"].values)

    @property
    def input_units(self):
        """Transfer function input channel units"""
        return str(self._ds.coords["in_units"].values[0])

    @property
    def n_windows(self):
        """Number of time series data windows used"""
        return self._ds.attrs["n_winds"]

    @property
    def noise_channels(self):
        """Names of the noise channel for each xf"""
        return list(self._ds.coords["noise_chan"].values)

    def coh_signif(self, prob=0.95):
        """
        Return coherence significance level

        Args:
            prob (float): significance level (between 0 and 1)
        """
        return coherence_significance_level(self.n_windows, prob)

    def output_units(self, output_channel):
        """Transfer function output channel units"""
        oc = self._match_out_chan(output_channel)
        return str(self._ds.sel(output=oc).coords["out_units"].values)

    def noise_channel(self, output_channel):
        """Transfer function noise channel string"""
        oc = self._match_out_chan(output_channel)
        return str(self._ds.sel(output=oc).coords["noise_chan"].values)

    def frf(self, output_channel, zero_as_none=False):
        """
        Return frequency response function for the given output channel

        Args:
            output_channel (str): output channel name
            zero_as_none (bool): return non-calculated values as Nones instead
                of zeros
        """
        oc = self._match_out_chan(output_channel)
        xf = np.squeeze(self._ds["value"].sel(output=oc).values)
        if zero_as_none:
            xf[xf == 0] = None
        return xf

    def frf_wrt_counts(self, output_channel, zero_as_none=False):
        """
        Return frequency response function with respect to raw data counts

        Args:
            output_channel (str): output channel name
            zero_as_none (bool): return non-calculated values as Nones instead
                of zeros
        """
        oc = output_channel
        return self.frf(oc, zero_as_none) * self.response(oc)

    def corrector(self, output_channel, zero_as_none=False):
        """
        Return coherent channel correction factor for the given output channel

        Args:
            output_channel (str): output channel name
            zero_as_none (bool): return non-calculated values as Nones instead
                of zeros
        """
        oc = self._match_out_chan(output_channel)
        corr_mult = np.squeeze(self._ds["corr_mult"].sel(output=oc).values)
        return self.frf(output_channel, zero_as_none) * corr_mult

    def corrector_wrt_counts(self, output_channel, zero_as_none=False):
        """
        Return transfer function with respect to raw data counts

        Args:
            output_channel (str): output channel name
            zero_as_none (bool): return non-calculated values as Nones instead
                of zeros
        """
        oc = output_channel
        return self.corrector(oc, zero_as_none) * self.response(oc)

    def response(self, output_channel, zero_as_none=False):
        """
        Return tf response (output_channel_response/input_channel_response)

        Divide count-based response by this to get unit-based response
        Multiply unit-based response by this to get count-based response
        """
        oc = self._match_out_chan(output_channel)
        return np.squeeze(self._ds["response"].sel(output=oc).values)

    def uncert_mult(self, output_channel):
        """Return transfer function uncertainty as a fraction of the transfer function"""
        oc = self._match_out_chan(output_channel)
        return np.squeeze(self._ds["uncert_mult"].sel(output=oc).values)

    def uncertainty(self, output_channel):
        """Return transfer function uncertainty for the given output channel"""
        return self.frf(output_channel) * self.uncert_mult(output_channel)

    def uncertainty_wrt_counts(self, output_channel):
        """
        Return transfer function uncertainty with respect to counts
        """
        return self.frf_wrt_counts(output_channel) * self.uncert_mult(output_channel)

    @staticmethod
    def _check_chans(sdm, in_chan, out_chans):
        """
        Verify that in_chan and out_chan are in the SpectralDensity object
        """
        # Validate in_chan
        if not isinstance(in_chan, str):
            raise TypeError("Error: in_chan is not a str")
        in_chans = fnmatch.filter(sdm.channel_names, in_chan)
        if len(in_chans) == 0:
            raise ValueError(f'No matches for "{in_chan}" in {sdm.channel_names}')
        elif len(in_chans) > 1:
            raise ValueError(
                f'Multiple channel matches for "{in_chan}": {in_chans}'
            )
        in_chan = in_chans[0]

        # Validate out_chan
        if out_chans is None:
            # Select all channels except in_chan
            out_chans = [x for x in sdm.channel_names if not x == in_chan]
        if not isinstance(out_chans, list):
            raise TypeError("Error: out_chans is not a list")
        return in_chan, out_chans

    def _match_out_chan(self, value):
        """
        Returns output_channel matching string (may have *,? wildcards)

        Error if there is not exactly one response
        """
        # Validate in_chan
        if not isinstance(value, str):
            raise TypeError("Error: value is not a str")
        out_chans = fnmatch.filter(self.output_channels, value)
        if len(out_chans) == 0:
            raise ValueError(f'No output channel matches "{value}"')
        elif len(out_chans) > 1:
            raise ValueError('Multiple output channels match "{}": {}'
                             .format(value, out_chans))
        return out_chans[0]

    def _calcxf(self, spect_density, input, output, noise_chan="output",
                n_to_reject=1, min_freq=None, max_freq=None):
        """
        Calculate frequency response function between two channels
        
        Returns 0 for values where coherence is beneath signif level
        
        Uses equations from Bendat&Piersol, 1986 (BP86)

        Args:
            spect_density(:class: ~SpectralDensity): cross-spectral density
                matrix
            input (str): input channel name
            output (str): output channel name
            noise_channel (str): which channel has the noise
            n_to_reject (int): only use values for which more than this
                many consecutive coherences are above the 95% significance
                level (0 = use all)
            min_freq (float): set to zero for frequencies below this value
            max_freq (float): set to zero for frequencies above this value

        Returns:
            (tuple):
                H (numpy.array): frequency response function
                H_err_mult (numpy.array): uncertainty multiplier
                corr_mult (numpy.array): value to multiply H by when correcting
                    spectra
        """
        coh = spect_density.coherence(input, output)
        f = spect_density.freqs
        # H = Gxy / Gxx  # B&P Equation 6.69
        # H = self._zero_bad(H, coh, n_to_reject, f, min_freq, max_freq)
        # if noise_chan == "output":
        #     xf = H * coh
        #     xferr = np.abs(xf) * errbase
        # elif noise_chan == "input":
        #     xf = H / coh
        #     xferr = np.abs(xf) * errbase
        # elif noise_chan == "equal":
        #     xf = H
        #     xferr = np.abs(xf) * errbase
        # elif noise_chan == "unknown":
        #     xf = H
        #     # VERY ad-hoc error guesstimate
        #     maxerr = np.abs(coh ** (-1)) + errbase
        #     minerr = np.abs(coh) - errbase
        #     xferr = np.abs(xf * (maxerr - minerr) / 2)
        # else:
        #     raise ValueError(f'unknown noise channel: "{noise_chan}"')
        # return xf, xferr
        
        # Calculate Frequency Response Function (Incorrectly called
        # transfer function)
        if noise_chan == "output" or noise_chan == "equal":
            Gxx = spect_density.autospect(input)
            Gxy = spect_density.crossspect(input, output)
            if noise_chan == "output":
                H = Gxy / Gxx  # BP86 eqn 6.37
                corr_mult = np.ones(H.shape)  # No change
            elif noise_chan == "equal":
                # Derived from BP86 eqns 6.48, 6.49, 6.51 and 6.52
                H = (Gxy / Gxx) / np.sqrt(coh)
                corr_mult = np.sqrt(coh)
        elif noise_chan == "input":
            Gyy = spect_density.autospect(output)
            Gyx = spect_density.crossspect(output, input)
            H = Gyy / Gyx # BP86 eqn 6.42
            corr_mult = coh  # derived from BP86 eqns 6.44 and 6.46
        # elif noise_chan == "unknown":
        #     xf = H
        #     # VERY ad-hoc error guesstimate
        #     maxerr = np.abs(coh ** (-1)) + errbase
        #     minerr = np.abs(coh) - errbase
        #     xferr = np.abs(xf * (maxerr - minerr) / 2)
        else:
            raise ValueError(f'unknown noise channel: "{noise_chan}"')
        H = self._zero_bad(H, coh, n_to_reject, f, min_freq, max_freq)

        # Calculate uncertainty
        # Crawford et al. 1991 eqn 4,  from BP2010 eqn 9.90
        H_err_mult = np.sqrt((np.ones(coh.shape) - coh) / (2*coh*self.n_windows))
        return H, H_err_mult, corr_mult

    def plot(self, errorbars=True, show=True):
        """
        Plot transfer functions

        Args:
            errorbars (bool): plot error bars
            show (bool): show on the screen
        Returns:
            (numpy.ndarray): array of axis pairs (amplitude, phase)
        """
        inp = self.input_channel
        outputs = self.output_channels
        rows = 1
        cols = len(outputs)
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        fig, axs = plt.subplots(rows, cols, sharex=True)
        in_suffix = self._find_str_suffix(inp, outputs)
        for out_chan, j in zip(outputs, range(len(outputs))):
            axa, axp = self.plot_one(inp, out_chan, fig, (rows, cols),
                                     (0, j), show_ylabel=j == 0,
                                     errorbars=errorbars,
                                     title=f"{out_chan}/{in_suffix}",
                                     show_xlabel=True)
        ax_array[0, j] = (axa, axp)
        if show:
            plt.show()
        return ax_array

    def _plot_progress(self, spect_dens, show=True):
        """
        Plot transfer functions and the coherences that made them

        Args:
            spect_dens (:class:`tiskit:SpectralDensity`): spectral density
                matrix
            show (bool): show plot
        Returns:
            (numpy.ndarray): array of axis pairs (amplitude, phase)
        """
        inp = self.input_channel
        outputs = self.output_channels
        shape = (2, len(outputs))
        ax_array = np.ndarray(shape, dtype=tuple)
        fig, axs = plt.subplots(shape[0], shape[1], sharex=True)
        in_suffix = self._find_str_suffix(inp, outputs)
        for out_chan, j in zip(outputs, range(len(outputs))):
            axa, axp = self.plot_one(
                inp,
                out_chan,
                fig,
                shape,
                (0, j),
                show_ylabel=j == 0,
                show_xlabel=True,
                title=f"{out_chan}/{in_suffix}",
            )
            ax_array[0, j] = (axa, axp)
            axa, axp = spect_dens.plot_one_coherence(
                inp,
                out_chan,
                fig,
                shape,
                (1, j),
                show_ylabel=j == 0,
                show_xlabel=True,
                show_phase=False,
            )
            ax_array[1, j] = axa
        if show:
            plt.show()
        return ax_array

    def _find_str_suffix(self, inp, outps):
        """
        Find longest non-common suffix of inp, outps

        Args:
            inp (str): string to reduce
            outps: (list of str): list of base strings to compare
        Returns:
            result (str): longest non-commmon suffix
        """
        ii_ind = len(inp) - 1
        for oc, j in zip(outps, range(len(outps))):
            for ii in range(0, len(oc)):
                if not inp[ii] == oc[ii]:
                    break
            if ii < ii_ind:
                ii_ind = ii
            return inp[ii:]

    def plot_one(self, in_chan, out_chan,
                 fig=None, fig_grid=(1, 1), plot_spot=(0, 0), errorbars=True,
                 label=None, title=None, show_xlabel=True, show_ylabel=True):
        """
        Plot one transfer function

        Args:
            in_chan (str): input channel
            out_chan (str): output channel
            fig (:class: ~matplotlib.figure.Figure): figure to plot on, if
                None this method will plot on the current figure or create
                a new figure.
            fig_grid (tuple): this plot sits in a grid of this many
                (rows, columns)
            plot_spot (tuple): put this plot at this (row, column) of
                the figure grid
            errorbars (bool): plot as vertical error bars
            label (str): string to put in legend
            title (str): string to put in title
            show_xlabel (bool): put an xlabel on this subplot
            show_ylabel (bool): put a y label on this subplot

        Returns:
            tuple:
                transfer function amplitude plot
                transfer function phase plot
        """
        xf = self.frf(out_chan).copy()
        xferr = self.uncertainty(out_chan).copy()
        f = self.freqs
        if fig is None:
            fig = plt.gcf()
        # Plot amplitude
        fig.suptitle("Transfer Functions")
        ax_a = plt.subplot2grid(
            (3 * fig_grid[0], 1 * fig_grid[1]),
            (3 * plot_spot[0] + 0, plot_spot[1] + 0),
            rowspan=2,
        )
        xf[xf == 0] = None
        xferr[xf == 0] = None
        if errorbars is True:
            ax_a.errorbar(f, np.abs(xf), np.abs(xferr), fmt='none')
            if np.any(xf is not None):
                ax_a.set_yscale('log')
            ax_a.set_xscale('log')
        else:
            ax_a.loglog(f, np.abs(xf + xferr), color="blue", linewidth=0.5)
            ax_a.loglog(f, np.abs(xf - xferr), color="blue", linewidth=0.5)
            ax_a.loglog(f, np.abs(xf), color="black", label=label)
        # ax_a.loglog(f, np.abs(xf), label=f"'{out_chan}' / '{in_chan}'")
        ax_a.set_xlim(f[1], f[-1])
        if title is not None:
            ax_a.set_title(title)
        if label is not None:
            ax_a.legend()

        if show_ylabel:
            ax_a.set_ylabel("TF")
        ax_a.set_xticklabels([])
        # Plot phase
        ax_p = plt.subplot2grid(
            (3 * fig_grid[0], 1 * fig_grid[1]),
            (3 * plot_spot[0] + 2, plot_spot[1] + 0),
            sharex=ax_a,
        )
        ax_p.semilogx(f, np.degrees(np.angle(xf)))
        ax_p.set_ylim(-180, 180)
        ax_p.set_xlim(f[1], f[-1])
        ax_p.set_yticks((-180, 0, 180))
        if show_ylabel:
            ax_p.set_ylabel("Phase")
        else:
            ax_p.set_yticklabels([])
        if show_xlabel:
            ax_p.set_xlabel("Frequency (Hz)")
        return ax_a, ax_p

    def _zero_bad(self, H, coh, n_to_reject, f, min_freq, max_freq):
        """
        Set non-significant elements to zero

        Args:
            H (np.ndarray): one-D array
            coh (np.ndarray): absolute coherence (1D)
            n_to_reject: how many consecutive "coherent" values are needed to
                accept the value at thise indice?
            f (np.array): frequencies
            min_freq: set values to zero for frequencies below this
            max_freq: set values to zero for frequencies above this
        """
        assert isinstance(H, np.ndarray)
        assert isinstance(coh, np.ndarray)
        if min_freq is not None and max_freq is not None:
            if min_freq >= max_freq:
                raise ValueError("min_freq >= max_freq")
        goods = np.full(coh.shape, True)
        if min_freq is not None:
            goods[f < min_freq] = False
        if max_freq is not None:
            goods[f > max_freq] = False
        if n_to_reject > 0:
            goods = np.logical_and(goods, coh > self.coh_signif(0.95))
            # goods = coh > self.coh_signif(0.95)

            # for n == 1, should do nothing, for n == 2, shift once, etc
            if n_to_reject > 1:
                goods_orig = goods.copy()
                # calc n values for both sides: 2=>0, 3,4=>1, 5,6->2, etc
                both_sides = int(np.floor((n_to_reject - 1) / 2))
                # Shift to both sides
                for n in range(1, both_sides + 1):
                    goods = np.logical_and(
                        goods, np.roll(goods_orig, n)
                    )
                    goods = np.logical_and(
                        goods, np.roll(goods_orig, -n)
                    )
                if n_to_reject % 2 == 0:  # IF EVEN, roll in one from above
                    goods = np.logical_and(
                        goods,
                        goods_orig.roll(
                            f=-both_sides - 1, fill_value=goods_orig[-1]
                        ),
                    )
            H[~goods] = 0
        return H
