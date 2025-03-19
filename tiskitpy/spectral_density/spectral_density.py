"""
Spectral Density Class

*The ``from_ATACR()`` method is disabled because readthedocs had
a problem importing the `obstools` package using `pip`.  (the method was never
validated anyway)*
"""
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
from ..logger import init_logger, change_level
from ..cleaned_stream import CleanedStream
from ..utils import match_one_str, CleanSequence as CS

logger = init_logger()
np.seterr(all="ignore")


class SpectralDensity:
    """
    Class for spectral density functions.

    Warning:
        The standard constructor is rarely used, generate objects using
        :meth:`SpectralDensity.from_stream()`

    No public attributes, access data through provided methods
    """
    def __init__(self, freqs, seed_ids, chan_units, n_windows,
                 window_type, window_s, chan_clean_sequences=None,
                 ts_starttime=None, ts_endtime=None, starttimes=None,
                 data=None, instrument_responses=None):
        """
        Warning:
            This constructor is rarely used, use instead
            :meth:`SpectralDensity.from_stream()`
        Args:
            freqs (np.ndarray): frequencies
            seed_ids (list of str): seed_ids for each channel
            chan_units (list of str): channel physical units (e.g m/s^2, Pa)
            n_windows (int): of windows used to calculate spectra
            window_type (str): type of window used
            window_s (float): length of data windows in seconds
            chan_clean_sequences (list of list of str): clean_sequences
                applied to each channel.
            ts_starttime (:class:`obspy.core.UTCDateTime`): start of time
                series used for this object
            ts_endtime (:class:`obspy.core.UTCDateTime`): end of time series
                used for this object
            starttimes (list of UTCDateTime): starttime for each window
            data (:class:`np.ndarray`):
                one-sided spectral density functions.
                shape = (len(ids), len(ids), len(freqs)
                units = chan_units(i)*chan_units(j)/Hz
                The diagonal is used for "autospectra", the others for
                coherency and frequency response functions
            instrument_responses (:class:`np.ndarray`):
                instrument response for each channel.
                shape=(n_spects,n_freqs)
                units=(counts/chan_units)
        """
        n_f, n_ch = len(freqs), len(seed_ids)
        if data is None:
            data = np.zeros((n_ch, n_ch, n_f), dtype="complex")
        if instrument_responses is None:
            instrument_responses = np.ones((n_ch, n_f), dtype="complex")
        assert isinstance(chan_clean_sequences, list)
        for x in chan_clean_sequences:
            assert isinstance(x, list)
            for y in x:
                assert isinstance(y, str)
        _validate_dimensions(freqs, seed_ids, chan_units, starttimes,
                             chan_clean_sequences, data, instrument_responses)
        ids = [CS.tiskitpy_id(s,cs) for (s, cs) in zip(seed_ids, chan_clean_sequences)]
        self._ds = xr.Dataset(
            data_vars={
                "spectra": (("input", "output", "f"), data),
                "instrument_response": (("input", "f"), instrument_responses)
            },
            coords={
                "input": ids,
                "output": ids,
                "f": freqs,
                "in_units": ("input", chan_units),
                "out_units": ("output", chan_units),
                "seed_ids": ("input", seed_ids)
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
        self._clean_sequences = {k: v for k, v in zip(ids,
                                                      chan_clean_sequences)}

    def __str__(self):
        s = "SpectralDensity object:\n"
        s += f"\tids={self.ids}\n"
        s += "\tchannel_units={}\n".format([self.channel_units(ch)
                                            for ch in self.ids])
        f = self.freqs
        s += f"\t{len(f)} frequencies, from {f[0]:.3g} to {f[-1]:.3g}Hz\n"
        s += f"\tn_windows={self.n_windows}\n"
        s += f"\twindow_type={self.window_type}"
        return s

    def __eq__(self, other):
        return self._ds == other._ds

    @property
    def ids(self):
        """
        List of Channel ids: seed_ids + any cleaning information
        ** Previously called _channel_names **
        """
        outp = list(self._ds.coords["input"].values)
        assert outp == list(self._ds.coords["output"].values)
        return outp

    @property
    def seed_ids(self):
        """list of seed_ids"""
        return list(self._ds.coords["seed_ids"].values)

    @property
    def clean_sequences(self):
        """
        Returns:
            (list of list):
        """
        return [self._clean_sequences[key] for key in self.ids]

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
            logger.info('no starttime information, using first used window')
        if self._ds.ts_endtime is not None:
            ts_end = self._ds.ts_endtime
        else:
            ts_start = self.starttimes[-1] + self.window_seconds
            logger.info('no endtime information, using end of last window')

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
                    time_spans=None, avoid_spans=None, remove_eqs=False,
                    z_threshold=3, quiet=False):
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
            subtract_rf_suffix (str): suffix to add to channel names if rf
                is subtracted
            remove_eqs (bool or str): if str, filename of QuakeML file
                containing earthquakes to remove using default parameters of
                TimeSpans.remove_eqs().  If True, download and use earthquakes
                from USGS website.  If False, do not remove earthquakes. See
                TimeSpans.remove_eqs() for details.
            z_threshold (float or None): reject windows with z-score greater
                than this value.  None: no rejection
            quiet (bool): only output errors and beyond to console
        """
        if quiet is True:
            change_level(logger, 'console','ERROR')
        if not isinstance(stream, Stream):
            raise ValueError(f"stream is a {type(stream)}, not obspy Stream")
        stream = stream.copy()  # avoid modifying original stream

        if time_spans is not None:
            stream = CleanedStream(stream).tag('SPANS')
        if avoid_spans is not None:
            if time_spans is not None:
                raise RuntimeError("Provided both time_spans and avoid_spans")
            stream = CleanedStream(stream).tag('AVOID')
            time_spans = avoid_spans.invert(stream[0].stats.starttime,
                                            stream[0].stats.endtime)
        if remove_eqs is not False:
            avoid_eqs = TimeSpans.remove_eqs(remove_eqs)
            ts_avoided = avoid_eqs.invert(stream[0].stats.starttime,
                                          stream[0].stats.endtime)
            if time_spans is None:
                time_spans = ts_avoided
            else:
                time_spans.combine(time_spans, ts_avoided)
        if starttimes is not None and time_spans is not None:
            raise RuntimeError("Provided both starttimes and time spans")

        stream = _align_traces(stream)

        # Select windows
        sr = stream[0].stats.sampling_rate
        if (window_s * sr > stream[0].stats.npts):
            raise ValueError(
                'Requested window size > data length ({:g} > {:g} s)'
                .format(window_s, stream[0].stats.npts/sr))
        ws = int(_npow2(window_s * sr))
        if (ws > stream[0].stats.npts):
            logger.warning(f'window pts > data pts ({ws:d} > '
                           f'{stream[0].stats.npts:d} pts), reducing ...')
            while ws > stream[0].stats.npts:
                ws /= 2
            ws = int(ws)
            logger.warning(f'New window size={ws:d} pts')
        multfac = 2 / (ws * sr)  # Bendata&Piersol 1986 eqs 11.100 & 11.102
        # window_starts = WindowSelect(stream, ws, windowtype)

        # Calculate FFTs
        ft, evalresps, units = {}, {}, []
        tagged_stream = CS.seedid_tag(stream)
        ids = [tr.id for tr in tagged_stream]
        clean_seq_dict = {tr.id: tr.stats.get('clean_sequence',[]) for tr in tagged_stream}
        seed_ids = [tr.id for tr in stream]
        if not len(ids) == len(set(ids)):
            raise ValueError(f"stream has duplicate IDs: {ids}")
        for id in ids:  # Calculate Fourier transforms
            tr_st = tagged_stream.select(id=id)
            if len(tr_st) == 0:
                raise ValueError(f'{id=} not found in tagged stream = {tagged_stream.__str__()}')
            elif not len(tr_st) == 1:
                raise ValueError(f'{len(tr_st)} {id=}s found in tagged stream = {tagged_stream.__str__()}')
            tr = tr_st[0]
            ft[id], f, sts = SpectralDensity._windowed_rfft(
                tr, ws, ws, windowtype, starttimes, time_spans)
            # Transform fft to physical units
            ft[id], resp, evalresp, ft_units = _correct_instrument_response(
                ft[id], f, id, tr.stats, inv)
            units.append(ft_units)
            evalresps[id] = evalresp
            if resp is not None:
                tr.stats.response = resp

        n_winds = len(sts)
        # Remove outliers
        if starttimes is None and z_threshold is not None:
            n_winds_orig = len(sts)
            ft, sts = cls._remove_outliers(ft, sts, z_threshold)
            n_winds_new = len(sts)
            if not n_winds_new == n_winds:
                rejected = n_winds - n_winds_new
                logger.info(f'{z_threshold=}, rejected '
                            f'{100.*rejected/n_winds:.0f}% '
                            f'of windows ({rejected:d}/{n_winds:d})')
                n_winds = n_winds_new

        # Clean data
        new_clean_seq_dict = {}
        if data_cleaner is not None:  # clean data using correlated noise
            rf_list = data_cleaner.RFList  # the ids here are beyond those in ft
            ft, new_clean_seq_dict = rf_list.ft_subtract_rfs(ft, evalresps)
            clean_seq_dict = {k: clean_seq_dict.get(k,[]) + new_clean_seq_dict.get(k,[]) for k in ids}

        clean_seq_list = [clean_seq_dict.get(x, []) for x in ids]
        # Create object
        obj = cls(f,
                  seed_ids,
                  units,
                  n_winds,
                  windowtype,
                  ws,
                  clean_seq_list,
                  ts_starttime=min([x.stats.starttime for x in stream]),
                  ts_endtime=max([x.stats.endtime for x in stream]),
                  starttimes=sts)
        # Fill object with Cross-Spectral Density Functions
        for inp in ids:
            in_id = CS.tiskitpy_id(CS.seed_id(inp), clean_seq_dict.get(inp,[]))
            if evalresps[inp] is not None:
                obj.put_channel_instrument_response(in_id, evalresps[inp])
            for outp in ids:
                out_id = CS.tiskitpy_id(CS.seed_id(outp), clean_seq_dict.get(outp,[]))
                # (p 547, Bendat and Piersol, 2010)
                obj.put_crossspect(in_id, out_id, 2 * multfac * np.mean(
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
                logger.debug('{:d} of {:d} had z_score > {:g}: rejected'
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
            channel (str): channel id
        Returns:
            (:class:`numpy.ndarray`): auto-spectral density function
        """
        ch_id = self.channel_id(channel, "in_id")
        return np.abs(self._ds["spectra"].sel(input=ch_id,
                                              output=ch_id).values.flatten())

    def crossspect(self, in_id, out_id):
        """
        Get cross-spectral density function for the channels

        Args:
            in_id (str): input channel id
            out_id (str): output channel id
        Returns:
            (:class:`numpy.ndarray`): cross-spectral density function
        """
        if in_id == out_id:
            return self.autospect(in_id)
        ichn = self.channel_id(in_id, "in_id")
        ochn = self.channel_id(out_id, "out_id")
        return (self._ds["spectra"].sel(input=ichn,
                                        output=ochn).values.flatten())

    def channel_id(self, test_id, ch_identifier='id'):
        """
        Return channel id, verifying that it exists and is unique

        Can expand wildcards, if they match only one channel

        Args:
            test_id (str): channel id to search for
            ch_identifier (str): description of the kind of channel (input,
                output...), useful for error messages
        """
        if not isinstance(test_id, str):
            raise TypeError(f"{ch_identifier} is a {type(test_id)}, not a str")
        name = match_one_str(test_id, self.ids,
                             "test_id", "self.ids")
        return name

    def seed_id(self, id):
        return self._ds["seed_ids"].sel(input=id)

    def put_autospect(self, channel, auto_spect):
        """
        Equivalent to put_cross_spect(channel, channel, auto_spect)

        Args:
            channel (str): auto-spectra channel id
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
        if channel not in self.ids:
            raise ValueError(
                'channel "{}" is not in ids {}'.format(
                    channel, self.ids
                )
            )
        self._ds["spectra"].loc[
            dict(input=channel, output=channel)
        ] = auto_spect

    def replace_channel_id(self, channel, replacement):
        """
        Args:
            channel (str): original channel id
            replacement (str): replacement channel id
        """
        ids = self.ids
        ids[ids.index(channel)] = replacement
        self._ds["input"] = ids
        self._ds["output"] = ids
        self._clean_sequences[replacement] = self._clean_sequences.pop(channel)

    def put_crossspect(self, in_id, out_id, cross_spect):
        """
        Put data into one of the cross-spectra.  Also puts the complex
        conjugate in the symmetric index

        Args:
            in_id (str): cross-spectra input channel
            out_id (str): cross-spectra output channel
            cross_spect (:class:`numpy.ndarray`): a cross-spectral density
        """
        assert cross_spect.shape == self.freqs.shape
        assert cross_spect.dtype == "complex"
        assert in_id in self.ids
        assert out_id in self.ids
        self._ds["spectra"].loc[
            dict(input=in_id, output=out_id)
        ] = cross_spect
        if not in_id == out_id:
            self._ds["spectra"].loc[
                dict(input=out_id, output=in_id)
            ] = np.conj(cross_spect)

    def channel_instrument_response(self, channel):
        """
        Get channel's instrument response

        Args:
            channel (str): channel name
        Returns:
            (:class:`numpy.ndarray`):
        """
        return self._ds["instrument_response"].sel(input=channel)

    def put_channel_instrument_response(self, channel, instrument_response):
        """
        Put a channel's instrument response into the object

        Verifies that the instrument_response has the same shape as the
        object's `frequency` property and that it is of type=`complex`

        Args:
            channel (str): the channel name
            instrument_response (:class:`numpy.ndarray`): the instrument
                response
        """
        assert instrument_response.shape == self.freqs.shape
        assert instrument_response.dtype == "complex"
        assert channel in self.ids
        self._ds["instrument_response"].loc[dict(input=channel)] = instrument_response

    def channel_units(self, channel):
        """
        Get channel's input (physical) units

        Args:
            channel (str): the channel name
        Returns:
            (str): Channel's input units
        """
        if channel not in self._ds["spectra"].coords["input"]:
            raise ValueError("channel {} not found in spectra.input {}"
                             .format(
                                channel,
                                self._ds["spectra"].coords["input"].values))
        return str(
            self._ds["spectra"].sel(input=channel).coords["in_units"].values
        )

    def clean_sequence(self, channel):
        """
        Clean sequence applied to this channel

        Args:
            channel (str): a channel name
        Returns:
            (list): List of cleaners applied, in order
        """
        if channel not in self._clean_sequences.keys():
            raise ValueError("channel {} not found in clean_sequence keys {}"
                             .format(channel, self._clean_sequences.keys()))
        return self._clean_sequences[channel]

    def put_clean_sequence(self, channel, clean_sequence):
        """
        Put a channel's clean_sequence into the object

        Args:
            channel (str): the channel name
            clean_sequence (list of str): clean channels/methods, in order
        """
        if channel in self._clean_sequences:
            self._clean_sequences[channel] = clean_sequence
        else:
            logger.warning('tried to put a clean_sequence in non-existent '
                           f'channel "{channel}"')

    def units(self, in_id, out_id):
        """
        The units of the given cross-  or auto-spectra

        Args:
            in_id (str): input channel
            out_id (str): output channel
        Returns:
            (str): the units
        """
        in_units = self.channel_units(in_id)
        out_units = self.channel_units(out_id)
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
        in_chan = self.channel_id(in_chan, 'in_id')
        out_chan = self.channel_id(out_chan, 'out_id')
        if in_chan not in self._ds.input:
            raise ValueError('in_chan={} not in spectral density matrix {}'
                             .format(in_chan, self._ds.input))
        if out_chan not in self._ds.output:
            raise ValueError('out_chan={} not in spectral density matrix {}'
                             .format(out_chan, self._ds.output))
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

    @staticmethod
    def plots(sds,
              channel=None,
              line_kws=None,
              labels=None,
              x=None,
              overlay=True,
              plot_peterson=True,
              show=True,
              outfile=None,
              title=None,
              **fig_kw):
        """
        Plot overlaid autospectra of multiple SpectralDensity objects

        Args:
            sds (list): SpectralDensity functions to plot
            channel (str): Limit to the given channel
            line_kws(list of dict): Line keywords for each SpectralDensity function
            labels(list of dict): labels for each sd
        Other Properties:
            **kwargs: any arguments used in plot_autospectra, except
                overlay (always true)
        Returns:
            figinfo (list):
                fig
                axa (): amplitude axis
        """
        # Validate inputs
        if overlay is not True:
            logger.warning('You requested overlay=False, ignored!')
        line_kws, labels = _validate_plots_args(sds, line_kws, labels)

        rows, cols = 1, 1
        # fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
        fig = plt.figure(**fig_kw)
        if title is None:
            title = "Auto-spectra, multiple SpectralDensities"
        fig.suptitle(title)
        axa, axp = None, None
        first_time = True
        for sd, l_kw, label in zip(sds, line_kws, labels):
            assert isinstance(l_kw, dict)
            new_x = sd._get_validate_ids(x)
            for key, i in zip(new_x, range(len(new_x))):
                if channel is not None:
                    if not key.split('.')[-1] == channel:
                        continue
                axa, axp = sd.plot_one_spectra(
                    key,
                    key,
                    fig,
                    (1, 1),
                    (0, 0),
                    show_ylabel=True,   # Would be better to do only LAST time
                    show_xlabel=True,   # Would be better to do only LAST time
                    ax_a=axa,
                    ax_p=axp,
                    show_phase=False,
                    plot_peterson=plot_peterson,
                    annotate=False,
                    label=label,
                    **fig_kw,
                    **l_kw
                )
                first_time = False
        if channel is not None and first_time is True:
            logger.error(f'No channel matching {channel} found, nothing plotted')
        plt.legend(fontsize='x-small')
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return fig, axa

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
                units of :math:`(m/s^2)^2/Hz`
            show (bool): show on desktop
            outfile (str): save figure to this filename
            title (str): custom plot title
            fig_kw (dict): all additional keyword arguments (such as `figsize`
                and `dpi`) are passed to the `pyplot.figure` call
        Returns:
            (:class:`numpy.ndarray`): array of axis pairs (amplitude, phase)
        """
        x = self._get_validate_ids(x)
        if not overlay:
            rows, cols = _squarish_grid(len(x))
        else:
            rows, cols = 1, 1
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        # fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
        fig = plt.figure(**fig_kw)
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
                    annotate=False
                )
            ax_array[0, 0] = (axa, axp)
            plt.legend(fontsize='small')
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
                units of :math:`(m/s^2)^2/Hz`
            show_coherence (bool): show coherence as well
            fig_kw (dict): all additional keyword arguments (such as `figsize`
                and `dpi`) are passed to the `pyplot.figure` call
        Returns:
            :class:`numpy.ndarray`: array of axis pairs (amplitude, phase)
        """
        x = self._get_validate_ids(x)
        n_subkeys = len(x)
        rows, cols = n_subkeys, n_subkeys
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        # fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
        fig = plt.figure(**fig_kw)
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
        annotate=True,
        **plot_kws
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
            subplot_spot (tuple): put this plot at this (row, column) of
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
                units of :math:`(m/s^2)^2/Hz`
            outfile (str): save figure to this filename
            **plot_kws (dict): keywords to pass on to plot command

        Returns:
            (tuple): tuple containing
                - :class:`matplotlib.axes.axis`: amplitude plot axis
                - :class:`matplotlib.axes.axis`: phase plot axis
        """
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
        ax_a.semilogx(f, 10 * np.log10(np.abs(psd)), label=label, **plot_kws)
        ax_a.set_xlim(f[1], f[-1])
        if plot_peterson is True and PSD_units.lower() == "(m/s^2)^2":
            lownoise, highnoise = Peterson_noise_model(f, True)
            ax_a.semilogx(f, lownoise, "k--")
            ax_a.semilogx(f, highnoise, "k--")

        if label is not None and annotate is True:
            ax_a.annotate(label, (0.5, 0.98),  xycoords="axes fraction",
                          ha='center', va='top', fontsize='small',
                          bbox=dict(boxstyle='square', fc='w', ec='k',
                                    alpha=0.5))
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
            ax_p.semilogx(f, np.degrees(np.angle(psd), **plot_kws))
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

    @staticmethod
    def plots_coherences(sds,
              sds_names=None,
              line_kws=None,
              labels=None,
              x=None,
              y=None,
              display='sparse',
              overlay=False,
              show=True,
              label_by="chan",
              sort_by="chan",
              outfile=None,
              title=None,
              **fig_kw):
        """
        Plot overlaid coherences of multiple SpectralDensity objects

        Args:
            sds (list): SpectralDensity functions to plot.  Each must have
                same seed_ids
            sds_names (list): Names to give to each sd in the plot legend
            line_kws(list of dict): Line keywords for each SpectralDensity function
            labels(list of dict): Labels for each SpectralDensity function
        Other Properties:
            **kwargs: any arguments used in plot_coherences, except
                overlay (always true)
        Returns:
            figinfo (list):
                fig
                axa (): amplitude axis
        """
        # Validate inputs
        if overlay is True:
            logger.error("Can't overlay multiple coherences")
            overlay=False
        assert display=='sparse'
        line_kws, labels = _validate_plots_args(sds, line_kws, labels)
        if sds_names is None:
            sds_names = [max(x.clean_sequences, key=len) for x in sds]
        else:
            assert len(sds_names) == len(sds)

        x_inp, y_inp = x, y
        strfun = sds[0]._seedid_strfun(sort_by)
        x = sorted(sds[0]._get_validate_ids(x_inp), key=strfun)
        y = sorted(sds[0]._get_validate_ids(y_inp), key=strfun)
        # Copied from plot_coherences "sparse" option, there's got
        # to be a way not to repeat code
        rows, cols = len(x), len(y)
        reduce_display = False
        if x[0] == y[0]:
            reduce_display = True   # Can get rid of one row and column
            rows -= 1
            cols -= 1
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        ax_array.fill((None, None))
        fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
        net_sta = '.'.join(sds[0].seed_ids[0].split('.')[:2])
        fig.suptitle(f"{net_sta} Coherences, multiple SpectralDensitys")
        first_time = True
        for sd, l_kw in zip(sds, line_kws):
            sort_strfun = sd._seedid_strfun(sort_by)
            x = sorted(sd._get_validate_ids(x_inp), key=sort_strfun)
            y = sorted(sd._get_validate_ids(y_inp), key=sort_strfun)
            label_strfun = sd._seedid_strfun(label_by)
            assert isinstance(l_kw, dict)
            # new_x = sd._get_validate_ids(x)
            plotted = []
            for in_chan, i in zip(x, range(len(x))):
                new_row = True
                for out_chan, jbase in zip(y, range(len(y))):
                    j = jbase
                    if reduce_display==True:
                        j -= 1
                    if in_chan == out_chan or (out_chan, in_chan) in plotted:
                        if i < rows and j >= 0:
                            axs[i, j].axis('off')
                        continue
                    plotted.append((in_chan, out_chan))
                    in_chan_label = label_strfun(in_chan)
                    out_chan_label = label_strfun(out_chan)
                    title = out_chan_label if i == 0 else None
                    axa, axp = sd.plot_one_coherence(
                        in_chan,
                        out_chan,
                        fig,
                        (rows, cols),
                        (i, j),
                        ax_a=ax_array[i,j][0],
                        ax_p=ax_array[i,j][1],
                        show_ylabel=new_row & first_time,
                        show_xlabel=(i == j)  & first_time,
                        ylabel=in_chan_label,
                        title=title,
                        **l_kw
                    )
                    new_row = False
                    ax_array[i, j] = (axa, axp)
        i, j = rows-1, 0
        axs[i,j].axis('on')
        for sds_name, l_kw, label in zip(sds_names, line_kws, labels):
            if label is None:
                label = sds_name
            axs[i,j].plot([1,1],[1,1], label=label, **l_kw)
        axs[i,j].legend(fontsize='x-small')
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return fig, axa

    def plot_coherences(self, x=None, y=None, display='full', show=True,
                        outfile=None, label_by="full", sort_by="full",
                        overlay=False, **fig_kw):
        """
        Plot coherences

        Args:
            x (list of str or None): limit to the listed input channels
            y (list of str or None): limit to the listed output channels
            display (str): how to arrange plots:
                - "full": a row for every channel, a column for every channel,
                  every cell filled
                - "sparse": Only plot the upper diagonal
                - "minimal": Plot upper diagonal elements in the least
                  number of cells possible
                - "overlay": One plot with all upper diagonal elemetns overlain
            overlay (bool): [GRANDFATHERED]: same as display="overlay"
            show (bool): show on desktop
            outfile (str): save to the named file
            label_by (str): labels to put on x and y axes ('full', 'chan' or
                'loc-chan')
            sort_by (str): how to sort x and y axes ('full', 'chan' or
                'loc-chan')
            fig_kw (dict): all additional keyword arguments (such as `figsize`
                and `dpi`) are passed to the `pyplot.figure` call

        Returns:
            (:class:`numpy.ndarray`): array of axis pairs (amplitude, phase)
        """
        if display not in ('full', 'sparse', 'minimal', 'overlay'):
            raise ValueError(f'Unknown display value: "{display}"')
        strfun = self._seedid_strfun(sort_by)
        x = sorted(self._get_validate_ids(x), key=strfun)
        y = sorted(self._get_validate_ids(y), key=strfun)
        if overlay is True:
            logger.warning('parameter `overlay` is grandfathered, '
                           'use display="overlay"')
            if display == 'full':
                display = 'overlay'
        if display == 'full':
            rows, cols = len(x), len(y)
            ax_array = np.ndarray((rows, cols), dtype=tuple)
            # fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
            fig = plt.figure(**fig_kw)
            fig.suptitle("Coherences")
            strfun = self._seedid_strfun(label_by)
            for in_chan, i in zip(x, range(rows)):
                for out_chan, j in zip(y, range(cols)):
                    in_chan_label = strfun(in_chan)
                    out_chan_label = strfun(out_chan)
                    title = out_chan_label if i == 0 else None
                    axa, axp = self.plot_one_coherence(
                        in_chan, out_chan,
                        fig, (rows, cols), (i, j),
                        show_ylabel=j == 0,
                        show_xlabel=i == rows - 1,
                        ylabel=in_chan_label,
                        title=title,
                    )
                    ax_array[i, j] = (axa, axp)
        if display == 'sparse':
            rows, cols = len(x), len(y)
            reduce_display = False
            if x[0] == y[0]:
                reduce_display = True   # Can get rid of one row and column
                rows -= 1
                cols -= 1
            ax_array = np.ndarray((rows, cols), dtype=tuple)
            fig, axs = plt.subplots(rows, cols, sharex=True, **fig_kw)
            fig.suptitle("Coherences")
            strfun = self._seedid_strfun(label_by)
            plotted = []
            for in_chan, i in zip(x, range(len(x))):
                new_row = True
                for out_chan, jbase in zip(y, range(len(y))):
                    j = jbase
                    if reduce_display==True:
                        j -= 1
                    if in_chan == out_chan or (out_chan, in_chan) in plotted:
                        if i < rows and j >= 0:
                            axs[i, j].axis('off')
                        continue
                    plotted.append((in_chan, out_chan))
                    in_chan_label = strfun(in_chan)
                    out_chan_label = strfun(out_chan)
                    title = out_chan_label if i == 0 else None
                    axa, axp = self.plot_one_coherence(
                        in_chan, out_chan,
                        fig, (rows, cols), (i, j),
                        show_ylabel=new_row,
                        # show_xlabel=i == rows - 1,
                        show_xlabel=i == j,
                        ylabel=in_chan_label,
                        title=title,
                    )
                    new_row = False
                    ax_array[i, j] = (axa, axp)
        if display == 'minimal':
            combis = []
            for i in x:
                for o in y:
                    if i == o or (o, i) in combis:
                        continue
                    combis.append((i, o))
            rows = int(np.sqrt(len(combis)))
            if rows == 0:
                rows = 1
            cols = int(np.ceil(len(combis)/rows))
            ax_array = np.ndarray((rows, cols), dtype=tuple)
            # fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True,
            #                         **fig_kw)
            fig = plt.figure(**fig_kw)
            fig.suptitle("Coherencies")
            strfun = self._seedid_strfun(label_by)
            i, j = 0, 0
            for combi in combis:
                in_chan = combi[0]
                out_chan = combi[1]
                in_chan_label = strfun(in_chan)
                out_chan_label = strfun(out_chan)
                title = out_chan_label if i == 0 else None
                # Get unique part of in_chan_label
                for ctr in range(len(in_chan_label)):
                    if ctr > len(out_chan_label):
                        raise ValueError("strings match all of out_chan_label")
                    if in_chan_label[ctr] == out_chan_label[ctr]:
                        continue
                    in_chan_sublabel = in_chan_label[ctr:]
                axa, axp = self.plot_one_coherence(
                    in_chan, out_chan,
                    fig, (rows, cols), (i, j),
                    show_ylabel=j == 0,
                    show_xlabel=i == rows - 1,
                    ylabel='Coherence',
                    title=None
                )
                label = f'{out_chan_label}/{in_chan_sublabel}'
                axa.annotate(label, (0.5, 0.98),  xycoords="axes fraction",
                             ha='center', va='top', fontsize='small',
                             bbox=dict(boxstyle='square', fc='w', ec='k',
                                       alpha=0.5))
                ax_array[i, j] = (axa, axp)
                j += 1
                if j >= cols:
                    j = 0
                    i += 1
        elif display == 'overlay':
            ax_array = np.ndarray((1, 1), dtype=tuple)
            # fig, axs = plt.subplots(1, 1, sharex=True, **fig_kw)
            fig = plt.figure(**fig_kw)
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
                        in_chan, out_chan,
                        fig, (1, 1), (0, 0),
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
        self, in_chan, out_chan,
        fig=None, fig_grid=(1, 1), plot_spot=(0, 0),
        show_xlabel=True, show_ylabel=True,
        ax_a=None, ax_p=None,
        ylabel=None, label=None, title=None,
        show_phase=True,
        outfile=None, show=False, **kwargs):
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
            outfile (str): plot to the named file
            show (bool): show on the screen (False by default: parent function shows)
            kwargs (dict): values to pass on to plotting routines

        Returns:
            (tuple): tuple containing:
                - (:class:`matplotlib.axes.axis`): amplitude plot axis
                - (:class:`matplotlib.axes.axis`): phase plot axis
        """
        in_chan = self.channel_id(in_chan, 'in_id')
        out_chan = self.channel_id(out_chan, 'out_id')
        ds = self._ds["spectra"].sel(input=in_chan, output=out_chan)
        f = self._ds.coords["f"].values
        if fig is None:
            fig = plt.gcf()
        plot_shape = (3*fig_grid[0], fig_grid[1])
        # Plot amplitude
        if ax_a is None:
            if show_phase:
                ax_a = plt.subplot2grid(
                    plot_shape,
                    (3 * plot_spot[0], plot_spot[1]),
                    rowspan=2
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
            ax_a.legend(fontsize='small')
        if not show_ylabel:
            ax_a.set_yticklabels([])
        if title:
            ax_a.set_title(title, fontsize='small')

        # Plot phase
        if show_phase:
            if ax_p is None:
                ax_p = plt.subplot2grid(plot_shape,
                                        (3 * plot_spot[0] + 2, plot_spot[1]))
            ax_p.semilogx(f, np.angle(ds, deg=True), **kwargs)
            ax_p.set_ylim(-180, 180)
            ax_p.set_xlim(f[1], f[-1])
            ax_p.set_yticks((-180, 0, 180))
            if show_ylabel:
                # ax_p.set_ylabel('Phase')
                pass
            else:
                ax_p.set_yticklabels([])
            ax_a.set_xticklabels([])
            bottom_axis = ax_p
        else:
            ax_p = None
            bottom_axis = ax_a

        if show_xlabel:
            bottom_axis.set_xlabel("Frequency (Hz)")
        else:
            bottom_axis.set_xticklabels([])

        # Put ylabel across amplitude and phase plots
        if show_ylabel:
            rows, cols = fig_grid
            row, col = plot_spot
            fig.add_subplot(rows, cols, row * cols + col + 1, frameon=False)
            if ylabel is None:
                ylabel = "Coherence"
            plt.tick_params(labelcolor='none', which='both', top=False,
                            bottom=False, left=False, right=False)
            plt.ylabel(ylabel, fontsize='small')

        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_a, ax_p

    def _get_validate_ids(self, x):
        """
        If x is None, return list of all channel ids
        If x is a list, validate each id

        Args:
            x (list or None): channel ids to validate
        """
        if x is None:
            x = list(self._ds.coords["input"].values)
            x = self.ids
        else:
            for key in x:
                # if key not in list(self._ds.coords["input"].values):
                if key not in self.ids:
                    ValueError('key "{key}" not in id list')
        return x

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
                reloffs = [int(x + (s-st)*sr) for x in spanoffsets]
                offsets.extend([x for x in reloffs
                                if x >= 0 and x + ws <= npts])
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

        if len(offsets) == 0:
            logger.warning('No offsets returned')
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


def _validate_dimensions(freqs, seed_ids, chan_units, starttimes,
                         chans_cleaned, data, instrument_responses):
    """Validate dimensions of __init__() variables"""
    n_f, n_ch = len(freqs), len(seed_ids)
    assert freqs.size == (n_f)  # Make sure it's one dimensional
    assert freqs.dtype == "float"
    assert len(chan_units) == n_ch
    assert len(chans_cleaned) == n_ch
    for x in chan_units:
        assert isinstance(x, str)
    for x in seed_ids:
        assert isinstance(x, str)
    for x in chans_cleaned:
        if not isinstance(x, list):
            raise TypeError(f'chans_cleaned element is a {type(x)}, '
                            'not a list')
        for y in x:
            if not isinstance(y, str):
                raise TypeError('chans_cleaned subelement is a '
                                f'{type(y)}, not a str')
    if starttimes is not None:
        for x in starttimes:
            assert isinstance(x, UTCDateTime)
    assert data.shape == (n_ch, n_ch, n_f)
    assert data.dtype == "complex"
    assert instrument_responses.shape == (n_ch, n_f)
    assert instrument_responses.dtype == "complex"


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
        logger.warning('Unmasking masked data (usually a gap or overlap)')
        stream = stream.split().merge(fill_value='interpolate')

    if last_start >= first_end:
        raise ValueError("There are non-overlapping traces")
    if last_start - first_start > 1 / sampling_rate:
        logger.debug("Cutting up to {} seconds from trace starts"
                     .format(last_start-first_start))
    if last_end - first_end > 1 / sampling_rate:
        logger.debug("Cutting up to {} seconds from trace ends"
                     .format(last_end-first_end))
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
        cols = int(np.ceil(np.sqrt(n_elems)))
        rows = int(np.ceil(n_elems / cols))
    return rows, cols


def _subtract_rfs(fts, subtract_rfs):
    """
    Args:
        fts (dict): dictionary containing Fourier transforms for each channel.
            Each Fourier transform is N*ws, where ws is the window size and N
            is the mumber of windows
        subtract_rfs (list of :class:`.ResponseFunctions``): frequency response
            functions to subtract from channels as ffts are calculated
            (NOT SURE IF THIS IS MORE USEFUL THAN SUBTRACTING FROM THE FINAL
            SPECTRALDENSITY (LIKE ATACR), BUT NEED TO TEST)
    Returns:
        fts (dict): dictionary containg corrected Fourier transforms for each
            channel.
    """
    for srf in subtract_rfs:
        rfs = srf.rfs
        in_chan = rfs.input_channel
        fts_ic = in_chan.split("-")[0]  # take off any '-*' tail
        if not rfs.freqs.shape == fts[fts_ic].shape:
            ValueError("frequency response function and ft have different "
                       f"shapes ({rfs.freqs.shape} vs {fts[fts_ic].shape})")
        for out_chan in rfs.output_channels:
            fts_oc = out_chan.split("-")[0]
            fts[fts_oc] -= fts[fts_ic] * rfs.corrector(out_chan)
    return fts


def _correct_instrument_response(ft, f, id, stats, inv=None):
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
        instrument_response: the channel instrument_response
        evalresp: the channel instrument_response evaluated at frequencies f
        units: the channel instrument_response units
    """
    resp, evalresp, units = None, None, "Counts"
    if inv is not None:
        try:
            resp = inv.get_response(id, stats.starttime)
        except Exception:
            # remove subtraction codes from location code
            new_id = _remove_subtracted_loc(id)
            try:
                resp = inv.get_response(new_id, stats.starttime)
            except Exception:
                new_id = _remove_subtracted_chan(id)
                try:
                    resp = inv.get_response(new_id, stats.starttime)
                except Exception:
                    raise ValueError(f'No match found for "{new_id}" in inv')
    if resp is None and "response" in stats:
        resp = stats.response
    if resp is not None:
        if resp.response_stages is None:
            raise ValueError(f'channel {id} has no response_stages')
        elif len(resp.response_stages) == 0:
            raise ValueError(f'channel {id} has zero response_stages')
        if "pa" in resp.instrument_sensitivity.input_units.lower():
            evalresp = resp.get_evalresp_response_for_frequencies(f, "VEL")
            units = "Pa"
        else:
            evalresp = resp.get_evalresp_response_for_frequencies(f, "ACC")
            units = "m/s^2"
        ft /= evalresp
    return ft, resp, evalresp, units


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


def _validate_plots_args(sds, line_kws, labels):
    """ validate arguments passed to plots() and to plots_coherences()"""
    if not isinstance(sds, (list, tuple)):
        raise ValueError('sds is not a list or tuple')
    for i, sd in zip(range(len(sds)), sds):
        if not isinstance(sd, SpectralDensity):
            raise ValueError(f'sds[{i}] is not a SpectralDensity object')
        if i == 0:
            seed_ids = sorted(sds[0].seed_ids)
        else:
            if not (seed_ids == sorted(sd.seed_ids)):
                raise ValueError(f"sds[{i}].seed_ids={sd.seed_ids}"
                                 f" does not match {sds[0].seed_ids=}")
    if line_kws is not None:
        if not isinstance(line_kws, (list, tuple)):
            raise ValueError('line_kws is not a list or tuple')
        if not len(line_kws) == len(sds):
            raise ValueError(f'{len(line_kws)=} != {len(sds)=}')
    else:
        line_kws = [{} for x in sds]
    if labels is not None:
        if not isinstance(labels, (list, tuple)):
            raise ValueError('labels is not a list or tuple')
        if not len(labels) == len(sds):
            raise ValueError(f'{len(labels)=} != {len(sds)=}')
    else:
        labels = [None for x in sds]
    return line_kws, labels


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


if __name__ == "__main__":
    import doctest

    doctest.testmod()
