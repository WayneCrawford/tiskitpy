"""
Spectral Density Functions
"""
import logging
import fnmatch

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt

# from obstools.atacr import DayNoise, StaNoise  # removed for readthedocs
from .utils import _prol1pi, _prol4pi, coherence_significance_level
from obspy.core.stream import Stream
from obspy.core import UTCDateTime
from scipy import signal, stats

from .Peterson_noise_model import Peterson_noise_model
from ..time_spans import TimeSpans

np.seterr(all="ignore")

logging.basicConfig(level=logging.INFO)


class SpectralDensity:
    """
    Class for spectral density functions.

    The standard constructor is rarely used, generate objects using
    `SpectralDensity.from_stream()`

    No public attributes, access data through provided methods
    """
    def __init__(self, chan_names, freqs, chan_units, n_windows, window_type,
                 window_s, ts_starttime=None, ts_endtime=None, starttimes=None,
                 data=None, responses=None):
        """
        Args:
            chan_names (list of str): channel names
            freqs (np.ndarray): frequencies
            chan_units (list of str): channel physical units (e.g m/s^2, Pa)
            n_windows (int): of windows used to calculate spectra
            window_type (str): type of window used
            window_s (float): length of data windows in seconds
            ts_starttime (:class:`obspy.core.UTCDateTime`): start of time
                series used for this object
            ts_endtime (:class:`obspy.core.UTCDateTime`): end of time series
                used for this object
            starttimes (list of UTCDateTime): starttime for each window
            data (:class:`np.ndarray`):
                one-sided spectral density functions.
                shape = (len(chan_names), len(chan_names), len(freqs)
                units = chan_units(i)*chan_units(j)/Hz
            responses (:class:`np.ndarray`):
                instrument response for each channel.
                shape=(n_spects,n_freqs)
                units=(counts/chan_units)
        """
        # Validate Dimensions
        n_ch = len(chan_names)
        n_f = len(freqs)
        shape = (n_ch, n_ch, n_f)
        dims = ("input", "output", "f")
        responses_shape = (n_ch, n_f)
        responses_dims = ("input", "f")
        assert freqs.size == (n_f)  # Make sure it's one dimensional
        assert len(chan_units) == n_ch
        if starttimes is not None:
            for x in starttimes:
                assert isinstance(x, UTCDateTime)
        if data is not None:
            assert data.size() == shape
            assert data.dtype == "complex"
        else:
            data = np.zeros(shape, dtype="complex")
        if responses is not None:
            assert responses.shape == responses_shape
            assert responses.dtype == "complex"
        else:
            responses = np.ones((n_ch, n_f), dtype="complex")

        self._ds = xr.Dataset(
            data_vars={
                "spectra": (dims, np.zeros(shape, dtype="complex")),
                "response": (responses_dims, responses),
            },
            # dims=dims,
            coords={
                "input": chan_names,
                "output": chan_names,
                "f": freqs,
                "in_units": ("input", chan_units),
                "out_units": ("output", chan_units),
            },
            attrs={
                "n_windows": n_windows,
                "window_type": window_type,
                "ts_starttime": ts_starttime,
                "ts_endtime": ts_endtime,
                "starttimes": starttimes,
                "window_s": window_s,
                "long_name": "spectral density function",
                "units": "input units * output units / Hz",
                "description": "One-sided spectral density functions",
            },
        )

    def __str__(self):
        s = "SpectralDensity object:\n"
        s += f"\tchannel_names={self.channel_names}\n"
        s += "\tchannel_units={}\n".format([self.channel_units(ch)
                                            for ch in self.channel_names])
        f = self.freqs
        s += f"\t{len(f)} frequencies, from {f[0]:.3g} to {f[-1]:.3g}Hz\n"
        s += f"\tn_windows={self.n_windows}\n"
        s += f"\twindow_type={self.window_type}"
        return s

    def __eq__(self, other):
        return self._ds == other._ds

    @property
    def channel_names(self):
        """
        Channel names

        Returns:
            (list of str):
        """
        assert list(self._ds.coords["input"].values) == list(
            self._ds.coords["output"].values
        )
        return list(self._ds.coords["input"].values)

    @property
    def freqs(self):
        """
        Frequencies of the spectral density functions

        Returns:
            (:class:`numpy.ndarray`):
        """
        return self._ds.coords["f"].values

    @property
    def window_type(self):
        """
        The type of window used to calculate the spectral densities

        Returns:
            (str):
        """
        return self._ds.window_type

    @property
    def window_seconds(self):
        """
        length of each window, in seconds

        Returns:
            (float):
        """
        return self._ds.window_s

    @property
    def starttimes(self):
        """
        Start times for each data window used to calculate spectra

        Returns:
            (list of :class:`obspy.UTCDateTimes`):
        """
        return self._ds.starttimes

    @property
    def used_times(self):
        """
        TimeSpans object containing time_spans used during processing

        Returns:
            (:class:`obspy.TimeSpans`):
        """
        spans = [[x, x + self.window_seconds] for x in self.starttimes]

        return TimeSpans([[x, x + self.window_seconds]
                           for x in self.starttimes])

    @property
    def unused_times(self):
        """
        TimeSpans object containing time_spans rejected during processing

        Returns:
            (:class:`obspy.TimeSpans`):
        """
        if self._ds.ts_starttime is not None:
            ts_start = self._ds.ts_starttime
        else:
            ts_start = self.starttimes[0]
            logging.info('no starttime information, using first used window')
        if self._ds.ts_endtime is not None:
            ts_end = self._ds.ts_endtime
        else:
            ts_start = self.starttimes[-1] + self.window_seconds
            logging.info('no endtime information, using end of last used window')

        return self.used_times.invert(ts_start, ts_end)

    @property
    def n_windows(self):
        """
        The number of data windows used to calculate spectra

        Returns:
            (int):
        """
        return self._ds.n_windows

    @classmethod
    def from_stream(cls, stream, window_s=1000, windowtype="prol1pi",
                    inv=None, data_cleaner=None, starttimes=None,
                    time_spans=None, avoid_spans=None, z_threshold=3):
        """
        Calculate spectral density functions from the provided stream

        Should add a window selection algorithm, for now just steps by
        the window length

        Args:
            stream (:class:`obspy.core.stream.Stream`): data
            window_s (float): desired window length in seconds
            windowtype (str): window type, must be a valid
            inv (:class:`obspy.core.inventory.Inventory`): inventory containing
                instrument responses.  If none is found for the given channel,
                will look in the channel's stats.response object
            data_cleaner (:class:`DataCleaner`): Data cleaner to
                apply to channels as ffts are calculated
            starttimes (list of :class:`obspy.UTCDateTime`): Use provided
                start window times (ignores z_threshold). Incompatible with
                `time_spans`
            time_spans (:class:`TimeSpans`): only use windows in the provided
                time spans.  Incompatible with `starttimes`
            avoid_spans (:class:`TimeSpans`): Do NOT use data from the provided
                time spans.  Incompatible with `starttimes` and "time_spans"
            subtract_tf_suffix (str): suffix to add to channel names if tf
                is subtracted
            z_threshold (float or None): reject windows with z-score greater than this
                value.  None: no rejection
        """
        if not isinstance(stream, Stream):
            raise ValueError(f"stream is a {type(stream)}, not obspy Stream")
        if avoid_spans is not None:
            if time_spans is not None:
                raise RuntimeError("Provided both time_spans and avoid_spans")
            time_spans = avoid_spans.invert(stream[0].stats.starttime,
            stream[0].stats.endtime)
        if starttimes is not None and time_spans is not None:
            if avoid_spans is not None:
                raise RuntimeError("Provided both starttimes and avoid_spans")
            else:
                raise RuntimeError("Provided both starttimes and time_spans")

        stream = stream.copy()  # avoid modifying original stream
        stream = _align_traces(stream)

        # Select windows
        sr = stream[0].stats.sampling_rate
        if (window_s * sr > stream[0].stats.npts):
            raise ValueError(
                'Requested window size > data length ({:g} > {:g} s)'
                .format(window_s, stream[0].stats.npts/sr))
        ws = int(_npow2(window_s * sr))
        if (ws > stream[0].stats.npts):
            logging.warning(f'window pts > data pts ({ws:d} > '
                            f'{stream[0].stats.npts:d} pts), reducing ...')
            while ws > stream[0].stats.npts:
                ws /= 2
            ws = int(ws)
            logging.warning(f'New window size={ws:d} pts')
        multfac = 2 / (ws * sr)  # Bendata&Piersol 1986 eqs 11.100 & 11.102
        # window_starts = WindowSelect(stream, ws, windowtype)

        # Calculate FFTs
        ft, evalresps, units = {}, {}, []
        ids = [tr.id for tr in stream]
        if not len(ids) == len(set(ids)):
            raise ValueError("stream has duplicate IDs")
        for id in ids:  # Calculate Fourier transforms
            tr = stream.select(id=id)[0]
            ft[id], f, sts = SpectralDensity._windowed_rfft(
                tr, ws, ws, windowtype, starttimes, time_spans)
            # Transform fft to physical units
            ft[id], resp, evalresp, ft_units = _correct_response(
                ft[id], f, id, tr.stats, inv)
            units.append(ft_units)
            evalresps[id] = evalresp
            if resp is not None:
                tr.stats.response = resp
        # Remove outliers if starttimes not forced and z_threshold specified
        if starttimes is None and z_threshold is not None:
            n_winds_orig = len(sts)
            ft, sts = cls._remove_outliers(ft, sts, z_threshold)
            n_winds = len(sts)
            if n_winds < n_winds_orig:
                rejected = n_winds_orig - n_winds
                logging.info(f'{z_threshold=} rejected {rejected:d} of '
                             f'{n_winds_orig:d} windows '
                             f'({100.*rejected/n_winds_orig:.0f}%)')

        if data_cleaner is not None:  # clean data using correlated noise
            dctfs = data_cleaner.DCTFs
            old_ids = ids
            ft = dctfs.ft_subtract_tfs(ft, evalresps)
            ids = dctfs.update_channel_names(old_ids)
            evalresps = dctfs.update_channel_keys(evalresps)
        # Create DataArray
        obj = cls(ids, f, units, n_winds, windowtype, ws,
                  ts_starttime=min([x.stats.starttime for x in stream]),
                  ts_endtime=max([x.stats.endtime for x in stream]),
                  starttimes=sts)
        # Fill DataArray with Cross-Spectral Density Functions
        for inp in ids:
            if evalresps[inp] is not None:
                obj.put_channel_response(inp, evalresps[inp])
            for outp in ids:
                # (p 547, Bendat and Piersol, 2010)
                obj.put_crossspect(inp, outp, 2 * multfac * np.mean(
                    np.conj(ft[inp])*ft[outp], axis=0))
        return obj

    @staticmethod
    def _remove_outliers(ft, sts, z_threshold=3, recursive=True):
        """Remove  windows with z-score above z_threshold

        Args:
            ft (dict): Fourier transforms.  Each value is an np.array where
                axis 0 = windows and axis 1 = frequencies m windows x n
                frequencies
            sts (list of :class:`UTCDateTime``): starttimes for each Fourier
                transform
            z_threshold (float): threshold.
            recursive (bool): repeat the test until there are no outliers
        """
        ft = ft.copy()
        sts = sts.copy()
        n_reject = 1  # Just starting the motor...
        while n_reject > 0:
            # norm should have as many columns as windows, as many rows as keys
            norm = np.array([np.mean(np.log10(np.abs(ft[id]*np.conj(ft[id]))),
                                     axis=1) for id in ft.keys()])
            # calculate a score for each channel (row) and window(column)
            z_scores = stats.zscore(norm, axis=1)
            # flatten down to one score for each window
            z_scores = np.max(np.abs(z_scores), axis=0)
            n_reject = np.count_nonzero(np.abs(z_scores) > z_threshold)
            if n_reject > 0:
                logging.debug('{:d} of {:d} had z_score > {:g}: rejected'
                              .format(n_reject, len(z_scores), z_threshold))
                keepers = z_scores <= z_threshold
                sts = [x for x, keep in zip(sts, keepers.tolist())
                       if keep is True]
                for id in ft.keys():
                    ft[id] = ft[id][keepers, :]
            if recursive is False:
                break
        return ft, sts

    def autospect(self, channel):
        """
        Get auto-spectral_density function for the channel

        Args:
            channel (str): channel name
        Returns:
            (:class:`numpy.ndarray`): auto-spectral density function
        """
        ch_name = self.channel_name(channel, "in_channel")
        return np.abs(self._ds["spectra"].sel(input=ch_name,
                                              output=ch_name).values.flatten())

    def crossspect(self, in_channel, out_channel):
        """
        Get cross-spectral density function for the channels

        Args:
            in_channel (str): input channel name
            out_channel (str): output channel name
        Returns:
            (:class:`numpy.ndarray`): cross-spectral density function
        """
        ichn = self.channel_name(in_channel, "in_channel")
        ochn = self.channel_name(out_channel, "out_channel")
        return (self._ds["spectra"].sel(input=ichn,
                                        output=ochn).values.flatten())

    def channel_name(self, channel, ch_identifier='chname'):
        """
        Return channel name, verifying that it exists and is unique

        Can expand wildcards, if they match only one channel
        
        Args:
            channel (str): channel name to search for
            ch_identifier (str): description of the kind of channel (input,
                output...), useful for error messages
        """
        if not isinstance(channel, str):
            raise TypeError(f"{ch_identifier} is a {type(channel)}, not a str")
        names = fnmatch.filter(self.channel_names, channel)
        if len(names) == 0:
            raise ValueError(f'{ch_identifier} "{channel}" not in '
                             f'{self.channel_names=}')
        elif len(names) > 1:
            raise ValueError(f'{ch_identifier} "{channel}" matches more than '
                             f'one {self.channel_names=}')
        else:
            return names[0]

    def put_autospect(self, channel, auto_spect):
        """
        Equivalent to put_cross_spect(channel, channel, auto_spect)

        Args:
            channel (str): auto-spectra channel
            auto_spect (:class:`numpy.ndarray`): the auto-spectral density
        """
        if not auto_spect.shape == self.freqs.shape:
            raise ValueError(
                "auto_spect has different shape than freqs ("
                f"{auto_spect.shape} vs {self.freqs.shape})"
            )
        if not auto_spect.dtype == "complex":
            try:
                auto_spect = auto_spect.astype(dtype=complex)
            except Exception:
                raise ValueError(
                    'auto_spect could not be converted to dtype=complex')
        if channel not in self.channel_names:
            raise ValueError(
                'channel "{}" is not in channel_names {}'.format(
                    channel, self.channel_names
                )
            )
        self._ds["spectra"].loc[
            dict(input=channel, output=channel)
        ] = auto_spect

    def replace_channel_name(self, channel, replacement):
        """
        Args:
            channel (str): original channel name
            replacement (str): replacement channel name
        """
        channel_names = self.channel_names
        channel_names[channel_names.index(channel)] = replacement
        self._ds["input"] = channel_names
        self._ds["output"] = channel_names

    def put_crossspect(self, in_channel, out_channel, cross_spect):
        """
        Put data into one of the cross-spectra.  Also puts the complex
        conjugate in the symmetric index

        Args:
            in_channel (str): cross-spectra input channel
            out_channel (str): cross-spectra output channel
            cross_spect (:class:`numpy.ndarray`): a cross-spectral density
        """
        assert cross_spect.shape == self.freqs.shape
        assert cross_spect.dtype == "complex"
        assert in_channel in self.channel_names
        assert out_channel in self.channel_names
        self._ds["spectra"].loc[
            dict(input=in_channel, output=out_channel)
        ] = cross_spect
        if not in_channel == out_channel:
            self._ds["spectra"].loc[
                dict(input=out_channel, output=in_channel)
            ] = np.conj(cross_spect)

    def channel_response(self, channel):
        """
        Get channel's instrument response

        Args:
            channel (str): channel name
        Returns:
            (:class:`numpy.ndarray`):
        """
        return self._ds["response"].sel(input=channel)

    def put_channel_response(self, channel, response):
        """
        Put a channel's instrument response into the object

        Verifies that the response has the same shape as the object's
        `frequency` property and that it is of type=`complex`

        Args:
            channel (str): the channel name
            response (:class:`numpy.ndarray`): the response
        """
        assert response.shape == self.freqs.shape
        assert response.dtype == "complex"
        assert channel in self.channel_names
        self._ds["response"].loc[dict(input=channel)] = response

    def channel_units(self, channel):
        """
        Get channel's input (physical) units
        
        Args:
            channel (str): the channel name
        Returns:
            (str): Channel's input units
        """
        return str(
            self._ds["spectra"].sel(input=channel).coords["in_units"].values
        )

    def units(self, in_channel, out_channel):
        """
        The units of the given cross-  or auto-spectra

        Args:
            in_channel (str): input channel
            out_channel (str): output channel
        Returns:
            (str): the units
        """
        in_units = self.channel_units(in_channel)
        out_units = self.channel_units(out_channel)
        if in_units == out_units:
            return f"({in_units})^2/Hz"
        return f"({in_units})*({out_units})/Hz"

    #     @classmethod
    #     def from_ATACR(cls, objnoise, horizontal_format='separate'):
    #         """
    #         Initiate class from ATACR DayNoise or StaNoise class
    #
    #         Args:
    #             objnoise (:class:`.DayNoise` or :class:`.StaNoise`):
    #                 noise spectra and frequencies
    #             horizontal_format (str): which type of horizontal channels
    #                 to use:
    #                     'aligned': one horizontal channel, aligned with
    #                         highest noise
    #                     'separate': two orthogonal horizontal channels
    #         """
    #         if (not objnoise and not isinstance(objnoise, DayNoise) and
    #                 not isinstance(objnoise, StaNoise)):
    #             raise TypeError("Error: A TFNoise object must be initialized"
    #                             " with only one of type DayNoise or "
    #                             "StaNoise object")
    #
    #         if not objnoise.av:
    #             raise(Exception("Error: Noise object has not been processed"
    #                             " (QC and averaging) - aborting"))
    #
    #         if horizontal_format == 'separate':
    #             chans = ['1', '2', 'Z', 'P']
    #             units = ['m/s^2', 'm/s^2', 'm/s^2', 'Pa']
    #         elif horizontal_format == 'aligned':
    #             chans = ['L', 'Z', 'P']
    #             units = ['m/s^2', 'm/s^2', 'Pa']
    #         else:
    #             raise ValueError('horizontal_format not "separate" or '
    #                              '"aligned"')
    #         # shape = (len(chans), len(chans), len(objnoise.f))
    #         if hasattr(objnoise, 'nwins'):
    #             n_winds = np.sum(objnoise.nwins)
    #         else:
    #             n_winds = np.sum(objnoise.goodwins)
    #
    #         obj = cls(chans, objnoise.f, units, n_winds, 'hanning')
    #         if horizontal_format == 'separate':
    #             obj.put_crossspect('1', '1', objnoise.power.c11)
    #             obj.put_crossspect('1', '2', objnoise.cross.c12)
    #             obj.put_crossspect('1', 'Z', objnoise.cross.c1Z)
    #             obj.put_crossspect('1', 'P', objnoise.cross.c1P)
    #             obj.put_crossspect('2', '1', np.conj(objnoise.cross.c12))
    #             obj.put_crossspect('2', '2', objnoise.power.c22)
    #             obj.put_crossspect('2', 'Z', objnoise.cross.c2Z)
    #             obj.put_crossspect('2', 'P', objnoise.cross.c2P)
    #             obj.put_crossspect('Z', '1', np.conj(objnoise.cross.c1Z))
    #             obj.put_crossspect('Z', '2', np.conj(objnoise.cross.c2Z))
    #             obj.put_crossspect('Z', 'Z', objnoise.power.cZZ)
    #             obj.put_crossspect('Z', 'P', objnoise.cross.cZP)
    #             obj.put_crossspect('P', '1', np.conj(objnoise.cross.c1P))
    #             obj.put_crossspect('P', '2', np.conj(objnoise.cross.c2P))
    #             obj.put_crossspect('P', 'Z', np.conj(objnoise.cross.cZP))
    #             obj.put_crossspect('P', 'P', objnoise.power.cPP)
    #         elif horizontal_format == 'aligned':
    #             obj.put_crossspect('L', 'L', objnoise.rotation.cHH)
    #             obj.put_crossspect('L', 'Z', objnoise.rotation.cHZ)
    #             obj.put_crossspect('L', 'P', objnoise.rotation.cHP)
    #             obj.put_crossspect('Z', 'L', np.conj(objnoise.rotation.cHZ))
    #             obj.put_crossspect('Z', 'Z', objnoise.power.cZZ)
    #             obj.put_crossspect('Z', 'P', objnoise.cross.cZP)
    #             obj.put_crossspect('P', 'L', np.conj(objnoise.rotation.cHP))
    #             obj.put_crossspect('P', 'Z', np.conj(objnoise.cross.cZP))
    #             obj.put_crossspect('P', 'P', objnoise.power.cPP)
    #         return obj

    def coherence(self, in_chan, out_chan):
        """
        Get coherence between the channels

        Args:
            in_chan (str): input channel.  Must match one of the
                coordinates in _ds
            out_chan (str): output channel.  Must match one of the
                coordinates in _ds
        Returns:
            (:class:`numpy.ndarray`): Coherence absolute value

        Coherence is a real-valued quantity, for the cross-spectral phase,
        use the cross-spectral density function.
        From Bendat & Piersol (1986), Appendix B, gamma_xy^2 (f)
        """
        in_chan = self.channel_name(in_chan, 'in_chan')
        out_chan = self.channel_name(out_chan, 'out_chan')
        if in_chan not in self._ds.input:
            raise ValueError(f'{in_chan=} not in spectral density matrix {self._ds.input}')
        if out_chan not in self._ds.output:
            raise ValueError(f'{out_chan=} not in spectral density matrix {self._ds.output}')
        coherence = (np.abs(self.crossspect(in_chan, out_chan))**2
                     / (self.autospect(in_chan) * self.autospect(out_chan)))
        return np.abs(coherence)  # shouldn't be necessary

    def coh_signif(self, prob=0.95):
        """
        Get the coherence significance level

        Args:
            prob (float): probability of an incoherent frequency
                passing this value (between 0 and 1)
        Returns:
            (float):
        """
        return coherence_significance_level(self.n_windows, prob)

    def plot(self, **kwargs):
        """Shortcut for `plot_autospectra()`"""
        self.plot_autospectra(**kwargs)

    def plot_autospectra(
        self,
        x=None,
        overlay=False,
        plot_peterson=True,
        show=True,
        outfile=None,
        title=None,
        **fig_kw
    ):
        """
        Plot autospectra

        Args:
            x (list of str): limit to the listed channels
            overlay (bool): put all spect on one axis
            plot_peterson(bool): plot Peterson Noise model if any channel has
                units of (m/s^2)^2/Hz
            show (bool): show on desktop
            outfile (str): save figure to this filename
            title (str): custom plot title
            **fig_kw (**dict): all additional keyword arguments (such as `figsize`
                and `dpi`) are passed to the `pyplot.figure` call
        Returns:
            (:class:`numpy.ndarray`): array of axis pairs (amplitude, phase)
        """
        x = self._get_validate_channel_names(x)
        if not overlay:
            rows, cols = _squarish_grid(len(x))
        else:
            rows, cols = 1, 1
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
        if title is None:
            title = "Auto-spectra"
        fig.suptitle(title)
        if not overlay:
            for key, i in zip(x, range(len(x))):
                i_row = int(i / cols)
                i_col = i - cols * i_row
                axa, axp = self.plot_one_spectra(
                    key,
                    key,
                    fig,
                    (rows, cols),
                    (i_row, i_col),
                    show_ylabel=i_col == 0,
                    show_xlabel=i_row == rows - 1,
                    show_phase=False,
                    plot_peterson=plot_peterson,
                )
                ax_array[i_row, i_col] = (axa, axp)
        else:
            axa, axp = None, None
            for key, i in zip(x, range(len(x))):
                axa, axp = self.plot_one_spectra(
                    key,
                    key,
                    fig,
                    (1, 1),
                    (0, 0),
                    show_ylabel=i == len(x) - 1,
                    show_xlabel=i == len(x) - 1,
                    ax_a=axa,
                    ax_p=axp,
                    show_phase=False,
                    plot_peterson=plot_peterson,
                )
            ax_array[0, 0] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_cross_spectra(
        self,
        x=None,
        show=True,
        show_coherence=False,
        outfile=None,
        plot_peterson=False,
        **fig_kw
    ):
        """
        Plot cross (and auto) spectra

        Args:
            x (list of str): limit to the listed channels
            show (bool): show on desktop
            plot_peterson(bool): plot Peterson Noise model if any channel has
                units of (m/s^2)^2/Hz
            show_coherence (bool): show coherence as well
            fig_kw (**dict): all additional keyword arguments (such as `figsize`
                and `dpi`) are passed to the `pyplot.figure` call
        Returns:
            :class:`numpy.ndarray`: array of axis pairs (amplitude, phase)
        """
        x = self._get_validate_channel_names(x)
        n_subkeys = len(x)
        rows, cols = n_subkeys, n_subkeys
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
        fig.suptitle("Cross-spectra (dB ref UNITS/Hz)")
        for in_chan, i in zip(x, range(len(x))):
            for out_chan, j in zip(x, range(len(x))):
                title = out_chan if i == 0 else None
                axa, axp = self.plot_one_spectra(
                    in_chan,
                    out_chan,
                    fig,
                    (rows, cols),
                    (i, j),
                    show_ylabel=j == 0,
                    show_xlabel=i == rows - 1,
                    ylabel=in_chan,
                    label="units",
                    title=title,
                    show_coherence=show_coherence,
                )
                ax_array[i, j] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_one_autospectra(self, key, **kwargs):
        """
        Plot one autospectral density

        Arguments are the same as for `plot_one_spectra()`, except
        there is no `subkey` argument and `show_phase` is ignored
        """
        kwargs["show_phase"] = False
        self.plot_one_spectra(key, key, **kwargs)

    def plot_one_spectra(
        self,
        key,
        subkey,
        fig=None,
        fig_grid=(1, 1),
        plot_spot=(0, 0),
        show_xlabel=True,
        show_ylabel=None,
        ax_a=None,
        ax_p=None,
        ylabel=None,
        label=None,
        title=None,
        show_coherence=False,
        show_phase=True,
        plot_peterson=True,
        outfile=None,
    ):
        """
        Plot one spectral density

        Args:
            key (str): input (driving) channel
            subkey (str): output (response) channel
            fig (:class:`matplotlib.figure.Figure`): figure to plot on, if
                None this method will plot on the current figure or create
                a new figure.
            fig_grid (tuple): this plot sits in a grid of this many
                              (rows, columns)
            subplot_spot (tuple): put this plot at this (row,column) of
                                  the figure grid
            show_xlabel (bool): put an xlabel on this subplot
            show_ylabel (bool): put a ylabel on this subplot
            ylabel (str): label to put on y axis (if show_label).  If not
                speficied, will use 'dB ref UNITS/Hz'
            label (str): 'units': print units without channel name in legend
            ax_a (Axis): use an existing axis for the amplitude plot
            ax_p (Axis): use this existing axis for the phase plot
            title (str): title to put on this subplot
            show_coherence (bool): draw coherence on the same plot
            show_phase (bool): show phase as well as amplitude
            plot_peterson(bool): plot Peterson Noise model if channel has
                units of (m/s^2)^2/Hz
            outfile (str): save figure to this filename

        Returns:
            (tuple): tuple containing
                - :class:`matplotlib.axes.axis`: amplitude plot axis
                - :class:`matplotlib.axes.axis`: phase plot axis
        """
        # da = self._ds["spectra"].sel(input=key, output=subkey)
        # in_units = da.coords['in_units'].values
        # out_units = da.coords['out_units'].values
        psd = self.crossspect(key, subkey)
        in_units = self.channel_units(key)
        out_units = self.channel_units(subkey)
        if in_units == out_units:
            PSD_units = f"({in_units})^2"
        else:
            PSD_units = f"{in_units}*{out_units}"
        f = self.freqs
        if fig is None:
            fig = plt.gcf()
        # Plot amplitude
        if ax_a is None:
            if show_phase:
                ax_a = plt.subplot2grid(
                    (3 * fig_grid[0], 1 * fig_grid[1]),
                    (3 * plot_spot[0] + 0, plot_spot[1] + 0),
                    rowspan=2,
                )
            else:
                ax_a = plt.subplot2grid(
                    (fig_grid[0], fig_grid[1]), (plot_spot[0], plot_spot[1])
                )
        if show_coherence:
            ax2 = ax_a.twinx()
            ax2.semilogx(
                f,
                np.abs(self.coherence(key, subkey)),
                color="red",
                linewidth=0.5,
                alpha=0.8,
            )
            ax2.axhline(
                self.coh_signif(0.95),
                color="red",
                linewidth=0.5,
                alpha=0.8,
                ls="--",
            )
            ax2.set_ylim(0, 1)
            if plot_spot[1] == fig_grid[1] - 1:  # Rightmost column
                ax2.set_ylabel("Coher", color="red")
            else:
                ax2.set_yticklabels([])
        psd[psd == 0] = None

        # Plot amplitude
        if label is None:
            label = f"{subkey} ({PSD_units})"
        elif label == "units":
            label = f"{PSD_units}"
        ax_a.semilogx(f, 10 * np.log10(np.abs(psd)), label=label)
        ax_a.set_xlim(f[1], f[-1])
        if plot_peterson is True and PSD_units.lower() == "(m/s^2)^2":
            lownoise, highnoise = Peterson_noise_model(f, True)
            ax_a.semilogx(f, lownoise, "k--")
            ax_a.semilogx(f, highnoise, "k--")

        if label is not None:
            legend_1 = ax_a.legend()
            if show_coherence:
                legend_1.remove()
                ax2.add_artist(legend_1)
        if show_ylabel:
            if ylabel is None:
                ylabel = "dB ref UNITS/Hz"
            ax_a.set_ylabel(ylabel)
        if title:
            ax_a.set_title(title)

        # Plot phase
        if show_phase:
            if ax_p is None:
                ax_p = plt.subplot2grid(
                    (3 * fig_grid[0], 1 * fig_grid[1]),
                    (3 * plot_spot[0] + 2, plot_spot[1] + 0),
                )
            ax_p.semilogx(f, np.degrees(np.angle(psd)))
            ax_p.set_ylim(-180, 180)
            ax_p.set_xlim(f[1], f[-1])
            ax_p.set_yticks((-180, 0, 180))
            if show_ylabel:
                # ax_p.set_ylabel('Phase')
                pass
            else:
                ax_p.set_yticklabels([])
            bottom_axis = ax_p
        else:
            ax_p = None
            bottom_axis = ax_a
        if show_xlabel:
            bottom_axis.set_xlabel("Frequency (Hz)")
        else:
            bottom_axis.set_xticklabels([])
        if outfile:
            plt.savefig(outfile)
        return ax_a, ax_p

    @staticmethod
    def _seedid_chan(inp):
        return inp.split('.')[-1]

    @staticmethod
    def _seedid_full(inp):
        return inp

    @staticmethod
    def _seedid_loc_chan(inp):
        return '.'.join(inp.split('.')[-2:])

    def _seedid_strfun(self, inp_type):
        """Choose the function that will be used to modify a seedid string"""
        allowed_inps = ('full', 'chan', 'loc-chan')
        if inp_type not in allowed_inps:
            raise ValueError(f'{inp_type} not in {allowed_inps}')
        if inp_type == 'chan':
            return self._seedid_chan
        elif inp_type == 'loc-chan':
            return self._seedid_loc_chan
        else:
            return self._seedid_full

    def plot_coherences(self, x=None, y=None, overlay=False, show=True,
                        outfile=None, labels="full", sort_by="full",
                        **fig_kw):
        """
        Plot coherences

        Args:
            x (list of str): limit to the listed input channels
            y (list of str): limit to the listed output channels
            overlay (bool): put all coherences on one plot
            show (bool): show on desktop
            outfile (str): save to the named file
            labels (str): labels to put on x and y axes ('full', 'chan' or
                'loc-chan')
            sort_by (str): how to sort x and y axes ('full', 'chan' or
                'loc-chan')
            fig_kw (**dict): all additional keyword arguments (such as `figsize`
                and `dpi`) are passed to the `pyplot.figure` call

        Returns:
            (:class:`numpy.ndarray`): array of axis pairs (amplitude, phase)
        """
        strfun = self._seedid_strfun(sort_by)
        x = sorted(self._get_validate_channel_names(x), key=strfun)
        y = sorted(self._get_validate_channel_names(y), key=strfun)
        if not overlay:
            rows, cols = len(x), len(y)
            ax_array = np.ndarray((rows, cols), dtype=tuple)
            fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
            fig.suptitle("Coherences")
            strfun = self._seedid_strfun(labels)
            for in_chan, i in zip(x, range(len(x))):
                for out_chan, j in zip(y, range(len(y))):
                    in_chan_label = strfun(in_chan)
                    out_chan_label = strfun(out_chan)
                    title = out_chan_label if i == 0 else None
                    axa, axp = self.plot_one_coherence(
                        in_chan,
                        out_chan,
                        fig,
                        (rows, cols),
                        (i, j),
                        show_ylabel=j == 0,
                        show_xlabel=i == rows - 1,
                        ylabel=in_chan_label,
                        title=title,
                    )
                    ax_array[i, j] = (axa, axp)
        else:
            ax_array = np.ndarray((1, 1), dtype=tuple)
            fig, axs = plt.subplots(1, 1, sharex=True, **fig_kw)
            fig.suptitle("Coherences")
            labels = []
            axa, axp = None, None
            for in_chan, i in zip(x, range(len(x))):
                for out_chan, j in zip(y, range(len(y))):
                    if out_chan == in_chan:
                        break
                    if f"{out_chan}-{in_chan}" in labels:
                        break
                    label = f"{in_chan}-{out_chan}"
                    axa, axp = self.plot_one_coherence(
                        in_chan,
                        out_chan,
                        fig,
                        (1, 1),
                        (0, 0),
                        ylabel="Coherence",
                        label=label,
                        ax_a=axa,
                        ax_p=axp,
                    )
                    labels.append(label)
            ax_array[0, 0] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_one_coherence(
        self,
        in_chan,
        out_chan,
        fig=None,
        fig_grid=(1, 1),
        plot_spot=(0, 0),
        show_xlabel=True,
        show_ylabel=True,
        ax_a=None,
        ax_p=None,
        ylabel=None,
        label=None,
        title=None,
        show_phase=True,
        **kwargs
    ):
        """
        Plot one coherence

        Args:
            in_chan (str): input (driving) channel
            out_chan (str): output (response) channel
            fig (:class:`matplotlib.figure.Figure`): figure to plot on, if
                None this method will plot on the current figure or create
                a new figure.
            fig_grid (tuple): this plot sits in a grid of this many
                              (rows, columns)
            plot_spot (tuple): put this plot at this (row,column) of
                                  the figure grid
            show_xlabel (bool): put an xlabel on this subplot
            show_ylabel (bool): put a ylabel on this subplot
            ylabel (str): label to put on y axis (if show_label).  If not
                speficied, will use 'dB ref UNITS/Hz'
            label (str): text to put in legend
            ax_a (Axis): use an existing axis for the amplitude plot
            ax_p (Axis): use this existing axis for the phase plot
            title (str): title to put on this subplot
            show_phase (bool): show phase as well as amplitude
            kwargs (**dict): values to pass on to plotting routines

        Returns:
            (tuple): tuple containing:
                - (:class:`matplotlib.axes.axis`): amplitude plot axis
                - (:class:`matplotlib.axes.axis`): phase plot axis
        """
        in_chan = self.channel_name(in_chan,'in_chan')
        out_chan = self.channel_name(out_chan, 'out_chan')
        ds = self._ds["spectra"].sel(input=in_chan, output=out_chan)
        f = self._ds.coords["f"].values
        if fig is None:
            fig = plt.gcf()
        # Plot amplitude
        if ax_a is None:
            if show_phase:
                ax_a = plt.subplot2grid(
                    (3 * fig_grid[0], 1 * fig_grid[1]),
                    (3 * plot_spot[0] + 0, plot_spot[1] + 0),
                    rowspan=2,
                )
            else:
                ax_a = plt.subplot2grid(
                    (fig_grid[0], fig_grid[1]), (plot_spot[0], plot_spot[1])
                )
        ax_a.semilogx(
            f, np.abs(self.coherence(in_chan, out_chan)), label=label, **kwargs
        )
        ax_a.axhline(
            self.coh_signif(0.95),
            color="red",
            linewidth=0.5,
            alpha=0.8,
            ls="--",
        )
        ax_a.set_ylim(0, 1)
        # ax_a.set_yticklabels([])
        # Plot amplitude
        if label is not None:
            ax_a.legend()
        if show_ylabel:
            if ylabel is None:
                ylabel = "Coherence"
            ax_a.set_ylabel(ylabel)
        if title:
            ax_a.set_title(title)

        # Plot phase
        if show_phase:
            if ax_p is None:
                ax_p = plt.subplot2grid(
                    (3 * fig_grid[0], 1 * fig_grid[1]),
                    (3 * plot_spot[0] + 2, plot_spot[1] + 0),
                )
            ax_p.semilogx(f, np.angle(ds, deg=True), **kwargs)
            ax_p.set_ylim(-180, 180)
            ax_p.set_xlim(f[1], f[-1])
            ax_p.set_yticks((-180, 0, 180))
            if show_ylabel:
                # ax_p.set_ylabel('Phase')
                pass
            else:
                ax_p.set_yticklabels([])
            bottom_axis = ax_p
        else:
            ax_p = None
            bottom_axis = ax_a
        if show_xlabel:
            bottom_axis.set_xlabel("Frequency (Hz)")
        else:
            bottom_axis.set_xticklabels([])
        return ax_a, ax_p

    def _get_validate_channel_names(self, x):
        """
        If x is None, return list of all channel names
        If x is a list, validate all of the names
        """
        if x is None:
            return list(self._ds.coords["input"].values)
        for key in x:
            if key not in list(self._ds.coords["input"].values):
                ValueError('key "{key}" not in channel list')
        return x

    @staticmethod
    def _remove_subtracted_loc(id):
        """
        Remove loc code characters including and after first "-"

        Allows the use of "-?" in the loc code to specify removed coherent
        noise

        Args:
            id (str): seed ID code

        Example:
            >>> SD._remove_subtracted_loc('hello')
            'hello'
            >>> SD._remove_subtracted_loc('NN.SSSS.LL.CCC')
            'NN.SSSS.LL.CCC'
            >>> SD._remove_subtracted_loc('NN.SSSS.-LL.CCC')
            'NN.SSSS..CCC'
            >>> SD._remove_subtracted_loc('NN.SSSS.L-LL.CCC')
            'NN.SSSS.L.CCC'
        """
        comps = id.split(".")
        if not len(comps) == 4:
            return id
        if "-" not in comps[2]:
            return id
        comps[2] = comps[2].partition("-")[0]
        return ".".join(comps)

    @staticmethod
    def _remove_subtracted_chan(id):
        """
        Remove chan code characters including and after first "-"
        """
        comps = id.split(".")
        if not len(comps) == 4:
            return id
        if "-" not in comps[3]:
            return id
        comps[3] = comps[3].partition("-")[0]
        return ".".join(comps)

    @staticmethod
    def _windowed_rfft(trace, ws, ss=None, win_taper="hanning",
                       starttimes=None, time_spans=None):
        """
        Calculates windowed Fourier transform of real data

        Args:
            trace (:class:`obspy.core.Trace`): Input trace data
            ws (int): Window size, in number of samples
            ss (int): Step size, or number of samples until next window
            win_taper (str): taper to apply to data  ['hanning', 'prol4pi',
                'prol1pi', 'bartlett', 'blackman', 'hamming']
            starttimes (list of :class:`obspy.UTCDateTime`): Use provided
                start window times (ignores z_threshold). Incompatible with
                `time_spans`
            time_spans (:class:`TimeSpans`): only use windows in the provided
                time spans.  Incompatible with `starttimes`

        Returns:
            (tuple): tuple containing:
                ft (:class:`numpy.ndarray`): Fourier transforms of trace, first
                    index corresponds to each window, second index to each freq
                f (:class:`numpy.ndarray`): Frequency axis in Hz
                t (list of UTCDateTime): Times of window starts
        """
        # Extract data windows
        a, starttimes = SpectralDensity._make_windows(trace, ws, ss, win_taper,
                                           starttimes, time_spans)
        sr = trace.stats.sampling_rate
        # Fourier transform
        n2 = _npow2(ws)
        ft = np.fft.rfft(a, n=n2)
        f = np.fft.rfftfreq(ws, 1.0 / sr)
        # f = np.linspace(0., 1., int(n2/2) + 1) * trace.stats.sampling_rate/2.
        # Don't return zero frequency
        return ft[:, 1:], f[1:], starttimes

    @staticmethod
    def _make_windows(trace, ws, ss, win_taper, starttimes, time_spans):
        """
        Returns:
            (tuple):
                windows (numpy.array): array of tapered windows (n_wind x ws)
                starttimes (list of :class:`obspy.UTCDateTime`): starttimes
                    for each window
        """
        st = trace.stats.starttime
        et = trace.stats.endtime
        sr = trace.stats.sampling_rate
        npts = trace.stats.npts

        # Calculate offsets
        if starttimes is not None:
            if time_spans is not None:
                raise RuntimeError("cannot provide time_spans AND starttimes")
            # Verify that times don't overrun data
            if np.any(starttimes < st):
                raise ValueError('provided starttimes before data start time')
            if np.any([x + ws/sr for x in starttimes] > et):
                raise ValueError('provided starttimes+ws after data end time')
            offsets = [int((x-st)*sr) for x in starttimes]
        elif time_spans is not None:
            offsets = []
            for s, e in zip(time_spans.start_times, time_spans.end_times):
                spanoffsets = SpectralDensity._sliding_window(int((e-s)*sr),
                                                              ws, ss)
                reloffsets = [int(x + (s-st)*sr) for x in spanoffsets]
                offsets.extend([x for x in reloffsets if x >= 0 and x+ws<=npts])
        else:
            offsets = SpectralDensity._sliding_window(trace.stats.npts, ws, ss)

        # make window taper
        if win_taper in ["hanning", "hamming", "blackman", "bartlett"]:
            taper = eval(f"np.{win_taper}(ws)")
        elif win_taper == "prol1pi":
            taper = _prol1pi(ws)
        elif win_taper == "prol4pi":
            taper = _prol4pi(ws)
        else:
            raise ValueError(f'Unknown taper type "{win_taper}"')

        if len(offsets)==0:
            logging.warning('No offsets returned')
            return None, None
        # Make tapered windows
        a = np.ndarray((len(offsets), ws), dtype=trace.data.dtype)
        for offset, i in zip(offsets, range(len(offsets))):
            if offset < 0:
                raise ValueError(f'{(offset+ws)=} > {trace.stats.npts=}')
            if offset+ws > trace.stats.npts:
                raise ValueError(f'{(offset+ws)=} > {trace.stats.npts=}')
            a[i] = signal.detrend(trace.data[offset:offset+ws]) * taper
        if starttimes is None:
            starttimes = [st + x/sr for x in offsets]

        return a, starttimes

    @staticmethod
    def _sliding_window(npts, ws, ss=None):
        """
        Return offsets for overlapping sub-windows

        Args:
            npts (int): number of data samples
            ws (int): Window size in samples
            ss (int): Step size in samples. If not provided, ss=ws

        Returns:
            offsets (list): list of sample offsets for window starts
        """
        if ss is None:
            ss = ws
        ws = int(ws)
        if (ws > npts):
            return []
        # Calculate the number of windows to return, ignoring leftover samples,
        nd = 1 + (npts - ws) // ss
        offsets = []
        if nd == 0:
            offsets.append(0)
        for i in range(nd):
            # "slide" the window along the samples
            offsets.append(i*ss)
        return offsets

def _align_traces(stream):
    """Trim stream so that all traces are aligned and same length"""
    # Determine last starttime and first endtime, and verify that
    # all traces have the same sample rate
    first_start = last_start = stream[0].stats.starttime
    first_end = last_end = stream[0].stats.endtime
    sampling_rate = stream[0].stats.sampling_rate
    for tr in stream[1:]:
        if tr.stats.starttime > last_start:
            last_start = tr.stats.starttime
        elif tr.stats.starttime < first_start:
            first_start = tr.stats.starttime
        if tr.stats.endtime < first_end:
            first_end = tr.stats.endtime
        elif tr.stats.endtime > last_end:
            last_end = tr.stats.endtime
        if not tr.stats.sampling_rate == sampling_rate:
            raise ValueError("not all traces have same sample rate")

    # Make sure that the array is not masked
    if np.any([np.ma.count_masked(x.data) for x in stream]):
        logging.warning('Unmasking masked data (usually a gap or overlap)')
        stream = stream.split().merge(fill_value='interpolate')

    if last_start >= first_end:
        raise ValueError("There are non-overlapping traces")
    if last_start - first_start > 1 / sampling_rate:
        logging.debug(f"Cutting up to {last_start-first_start} seconds "
                      "from trace starts")
    if last_end - first_end > 1 / sampling_rate:
        logging.debug(f"Cutting up to {last_end-first_end} seconds "
                      "from trace ends")
    stream.trim(last_start, first_end)
    min_len = min([tr.stats.npts for tr in stream])
    max_len = max([tr.stats.npts for tr in stream])
    if not max_len == min_len:
        for tr in stream:
            tr.data = tr.data[:min_len]
    return stream


def _npow2(x):
    """ Returns power of 2 >= x"""
    return 1 if x == 0 else 2 ** int(x - 1).bit_length()


def _squarish_grid(n_elems):
    """
    Returns as close to a square grid as fits the data without waste

    Args:
        n_elems (int): Number of elements in the grid

    Returns:
        (tuple): tuple containing:
            rows (int): number of rows in grid
            cols (int): number of columns in grid
    """
    n_elems = int(n_elems)
    if n_elems == 1:
        return 1, 1
    elif n_elems == 2:
        return 1, 2
    elif n_elems <= 4:
        return 2, 2
    elif n_elems <= 6:
        return 2, 3
    else:
        cols = np.ceil(np.sqrt(n_elems))
        rows = np.ceil(n_elems / cols)
    return rows, cols


def _subtract_tfs(fts, subtract_tfs):
    """
    Args:
        fts (dict): dictionary containing Fourier transforms for each channel.
            Each Fourier transform is N*ws, where ws is the window size and N
            is the mumber of windows
        subtract_tfs (list of DCTF): transfer functions to
            subtract from channels as ffts are calculated (NOT SURE IF
            THIS IS MORE USEFUL THAN SUBTRACTING FROM THE FINAL
            SPECTRALDENSITY (LIKE ATACR), BUT NEED TO TEST)
    Returns:
        fts (dict): dictionary containg corrected Fourier transforms for each
            channel.
    """
    for stf in subtract_tfs:
        tfs = stf.tfs
        in_chan = tfs.input_channel
        fts_ic = in_chan.split("-")[0]  # take off any '-*' tail
        if not tfs.freqs.shape == fts[fts_ic].shape:
            ValueError("transfer function and ft have different shapes "
                       f"({tfs.freqs.shape} vs {fts[fts_ic].shape})")
        for out_chan in tfs.output_channels:
            fts_oc = out_chan.split("-")[0]
            fts[fts_oc] -= fts[fts_ic] * tfs.corrector(out_chan)
    return fts


def _correct_response(ft, f, id, stats, inv=None):
    """
    Convert Fourier transform into channel in_units

    Args:
        ft (:class:`numpy.ndarray`): fourier transforms for one channel
        f (:class:`numpy.ndarray`): frequencies
        id (str): channel id
        stats (:class:`obspy.core.stream.stats`): trace statistics
        inv (:class:`obspy.core.inventory.Inventory`): station inventory
    Returns:
        ft: corrected fourier transforms
        response: the channel response
        evalresp: the channel response evaluated at frequencies f
        units: the channel response units
    """
    resp, evalresp, units = None, None, "Counts"
    if inv is not None:
        try:
            resp = inv.get_response(id, stats.starttime)
        except Exception:
            new_id = SpectralDensity._remove_subtracted_loc(id)
            try:
                resp = inv.get_response(new_id, stats.starttime)
            except Exception:
                new_id = SpectralDensity._remove_subtracted_chan(id)
                try:
                    resp = inv.get_response(new_id, stats.starttime)
                except Exception:
                    raise ValueError(f'No match found for "{new_id}" in inv')
    if resp is None and "response" in stats:
        resp = stats.response
    if resp is not None:
        if "pa" in resp.instrument_sensitivity.input_units.lower():
            evalresp = resp.get_evalresp_response_for_frequencies(f, "VEL")
            units = "Pa"
        else:
            evalresp = resp.get_evalresp_response_for_frequencies(f, "ACC")
            units = "m/s^2"
        ft /= evalresp
    return ft, resp, evalresp, units


if __name__ == "__main__":
    import doctest

    doctest.testmod()
