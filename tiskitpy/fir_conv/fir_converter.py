import copy
import json
import math
from typing import Tuple
from pathlib import Path

# from pprint import pprint

import numpy as np
from numpy.typing import ArrayLike
# from numpy.fft import fft
import matplotlib.pyplot as plt
import scipy.signal as sig
from obspy import UTCDateTime
from obspy.clients.filesystem.sds import Client as SDS_Client
from obspy.core.inventory import FIRResponseStage
from obspy.core.inventory.util import FloatWithUncertaintiesAndUnit

from ..util_functions.sds_tools import (sds_channel_files, sds_file_path,
                                        sds_neighbour_streams,
                                        get_neighbour_pad)
from ..util_functions.inventory_tools import find_channel, insert_channel

# from ..logger import init_logger

# logger = init_logger()


class FIRConverter:
    """
    Calculate and perform FIR conversion from zero phase to minimum phase

    Parameters:
        firzeros (ArrayLike) Input FIR filter zeros (polynomial coefficients).
        sps (float or None): Sampling rate (Hz) at which the FIR runs
            (before decimation).  Informational only (FDIG_EFF in
            Scherbaum's .prt files): not used in the conversion.
        timetag (float): Time-tag correction in samples.
        input_type (str): 'linear' or 'minphase'
        output_type (str): 'minphase' or 'linear'
        decimation_factor (int): decimation associated with this FIR
        z (np.ndarray): z-plane poles of firzeros
        a (np.ndarray):
        b (np.ndarray):

    Implements FIR filter conversion described in chapter 8 of Scherbaum, 1996,
    "Of poles and zeros"

    BASED ON Scherbaum C code:
    VERSION: 2.0
    DATE: 1996-10-18 (Frank Scherbaum)
    DESCRIPTION: This program implements the FIR filter conversion described in
    chapter 8 of Scherbaum, F., Of poles and zeros: fundamentals of digital
    seismology, Kluwer Academic Publishers, 1996. The key equation
    is (8.15) on page 120. The ARMA coefficients are assumed to be
    stored in a '*.prt' file produced by Frank Scherbaum's analyse.m
    Mathematica program. This is the  <filter coefficient file
    for conversion filter> above.

    Filter x[] which is the reversed input sequence x2[]
    using the difference equation:

            mx                mx
            --                -----
    y'[i] =  > a[k]*y'[i-k]   + > b[l] x[i-k]
            __                __
            k=1               l=0

    This corresponds to equ. (8.15) in 'Scherbaum, F: Of poles and zeros,
    Fundamentals of Digital Seismology, Kluwer Academic Publ., 1996'
    mx = number of AR coefficients
    b[l] = MA coefficients for l = 0, mx
    a[k] = AR coefficients for k = 1, mx

    x[i] = reversed input sequence for i = 0 .....
    y'[] = output sequence

    Reverse the output sequence y'[] in time again to obtain the
    corrected sequence y[n]!
    """

    def __init__(self, firzeros: ArrayLike, sps: float, timetag: float,
                 decimation_factor: int, input_type: str, 
                 a: np.ndarray, b: np.ndarray) -> None:
        """
        Constructor method:
        
        Args:
            firzeros: 
            sps (float or None): FIR sampling rate (informational only)
            timetag (float)
            decimation_factor (int)
            input_type (str)
        """
        self.firzeros: np.ndarray = np.asarray(firzeros, dtype=float)
        self.sps = None if sps is None else float(sps)
        self.timetag: float = float(timetag)
        self.input_type: str = input_type
        self.decimation_factor: int = decimation_factor
        self.a: np.ndarray = a
        self.b: np.ndarray = b
        self.z = np.roots(self.firzeros)
        self.equivalent_zeros = self._compute_minimum_phase_filter()

    # ------------------------------------------------------------
    # Dynamic properties
    # ------------------------------------------------------------
    @property
    def z_max(self):
        """Maximum phase z poles"""
        return self.z[np.abs(self.z) < 0.999]

    @property
    def z_min(self):
        """Minimum phase z poles"""
        return self.z[np.abs(self.z) > 1.001]

    @property
    def z_UC(self):
        """Unit Circle z poles"""
        return self.z[np.logical_and(np.abs(self.z) >= 0.999,
                                     np.abs(self.z) <= 1.001)]

    @property
    def output_type(self):
        """Equivalent FIR output type"""
        if self.input_type == 'linear':
            return 'minphase'
        elif self.input_type == 'minphase':
            return 'linear'
        else:
            raise ValueError(f'{self.input_type=} not in ("linear", "minphase")')

    # ------------------------------------------------------------
    # Alternative constructors
    # ------------------------------------------------------------
    @classmethod
    def from_zeros(cls, firzeros: ArrayLike, sps: float, timetag: float,
                   decimation_factor: int, input_firtype='linear'):
        """
        Calculate converter from provided FIR zeros
        """
        firzeros = np.asarray(firzeros, dtype=float)

        z = np.roots(firzeros)
        if input_firtype == 'linear':
            # fir_mp = cls._compute_minimum_phase_filter(z)
            a, b = cls._compute_AR_MA(z)
        elif input_firtype == 'minphase':
            # fir_mp = cls._compute_minimum_phase_filter(z)
            a, b = cls._compute_AR_MA(z)
        else:
            raise ValueError(f'Unknown {input_firtype=}')
        return cls(firzeros, sps, timetag, decimation_factor, input_firtype,
                   a, b)

    @classmethod
    def from_zeros_file(cls, filename: str, sps: float = None):
        """
        Read FIR coefficients from a JSON file, and calculate the conversion

        The JSON file must have the elements 'fir', 'timetag',
        'decimation_factor', and 'type'

        Args:
            filename (str): JSON file name
            sps (float or None): FIR sampling rate (informational only)
        """
        with open(filename, "r") as f:
            book = json.load(f)
        firzeros = np.array(book["fir"], dtype="double")
        if 'symmetry' in book:
            if book['symmetry'] == 'NONE':
                firzeros = firzeros
            elif book['symmetry'] == 'EVEN':
                firzeros = np.concatenate((firzeros, firzeros[-1::-1]))
            elif book['symmetry'] == 'ODD':
                firzeros =  np.concatenate((firzeros, firzeros[-2::-1]))
            else:
                raise ValueError(f'{filename}: Unknown {book["symmetry"]=}')
        timetag = float(book["timetag"])
        decimation_factor = int(book["decimation_factor"])
        input_firtype = book["type"]
        return cls.from_zeros(firzeros, sps, timetag, decimation_factor,
                              input_firtype)

    @classmethod
    def from_conv_file(cls, filename: str):
        """Read in values from a JSON conversion coefficients file"""
        with open(filename, "r") as f:
            book = json.load(f)
        firzeros = np.array(book["firdata_eff"], dtype="double")
        sps = book.get("sps")
        timetag = float(book["timetag"])
        input_type = book["input_type"]
        decimation_factor = int(book["decimation_factor"])
        a = np.array(book["ar"], dtype="double")
        b = np.array(book["ma"], dtype="double")
        return cls(firzeros, sps, timetag, decimation_factor, input_type,
                   a, b)

    @classmethod
    def from_builtin(cls, name=None):
        """
        Read in values from builtin conversion coefficients

        Args:
            name (str or None): if str, uses the specified builtin
                                if None, prints a list of builtins
        """
        corr_dir = Path(__file__).parent.resolve() / 'conv_coeffs'
        if name is None:
            children = list(corr_dir.glob('*.json'))
            if len(children) == 0:
                print('No built-in conversion coefficients')
            else:
                print('Built-in conversion coefficients:')
                for child in children:
                    print('\t' + child.stem)
            return
        return cls.from_conv_file(corr_dir / (name + '.json'))

    # ------------------------------------------------------------
    # Computational methods
    # ------------------------------------------------------------
    def _compute_minimum_phase_filter(self) -> np.ndarray:
        """
        Construct the equivalent minimum-phase FIR filter by
        reflecting maximum-phase roots inside the unit circle.

        Returns:
            (np.ndarray):  Minimum-phase FIR coefficients, normalized to
                sum to 1.
        """
        z_mp = self.z.copy()
        i_max = np.abs(self.z) < 0.999  # maximum phase zs
        # Invert max phase zs, making them min phase zs
        z_mp[i_max] = np.power(z_mp[i_max], -1)

        fir_mp = np.poly(z_mp)[::-1]
        fir_mp /= np.sum(fir_mp) # Normalize to sum = 1
        return fir_mp

    @staticmethod
    def _compute_AR_MA(z) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute AR and MA coefficients from the maximum-phase roots.

        Returns
        -------
        (np.ndarray, np.ndarray)
            Tuple of AR coefficients `a` and MA coefficients `b`.
        """
        F_max = z[np.abs(z) > 1.001]  # minimum phase zs
        fmax = np.poly(F_max)  # calculate coefficients of minimum phase part
        mx = len(F_max)
        # print(f'{F_max=}, {fmax=}, {mx=}')

        # a is -fmax[mx-1: -1: -1] / fmax[mx]
        a = -fmax[(mx - np.arange(1, mx + 1))] / fmax[mx]
        # b is fmax[0: mx+1] / fmax[mx]
        b = fmax[np.arange(0, mx + 1)] / fmax[mx]
        return a, b

    # ------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------
    def plot_zplane(self) -> None:
        """Plot the z-plane roots."""
        ZPlanePlotter.plot(self.z, self.z_max, self.z_min, self.z_UC,
                           "Input filter")

    def plot_output_zplane(self) -> None:
        """Plot the z-plane roots of the output filter"""
        # ATTENTION RECURSION!!!
        if self.input_type == 'linear':
            timetag = 0
        else:
            timetag = len(self.z)/2
        outp = FIRConverter.from_zeros(
            self.equivalent_zeros, self.sps, timetag,
            self.decimation_factor, self.output_type)
        ZPlanePlotter.plot(outp.z, outp.z_max, outp.z_min, outp.z_UC,
                           "Output filter")

    def plot_impulse_parts(self) -> None:
        """Plot impulse responses for each root category."""
        ImpulsePlotter.plot(self.firzeros, self.z, self.z_max, self.z_min,
                            self.z_UC, self.equivalent_zeros)

    # ------------------------------------------------------------
    # Data conversion
    # ------------------------------------------------------------
    def apply(self, stream, prev_stream=None, FIRdecim=None,
              new_loc_code='01', detrend=None, next_stream=None,
              gap_handling='interpolate', max_interp_gap=None,
              rel_tol=2.**-32, max_misalignment=0.5):
        """
        Convert a zero-phase stream to minimum-phase
        (if the object has the correct zero-phase information)

        The prev_stream, next_stream, gap_handling, max_interp_gap and
        rel_tol parameters are only used if you want to apply to consecutive
        data contained in different streams (or in different traces of the
        same stream).
        Each output sample depends on input samples from n_before samples
        before it to (in theory) infinitely after it, as the conversion
        filter is an anti-causal IIR filter (see `pad_lengths()`).  To get
        the same result as if the data had been converted in one piece,
        provide the data preceding (`prev_stream`) and following
        (`next_stream`) the input stream: the required number of samples
        will be taken from each and the output will be cut to the same
        samples as the input.  Other traces in `stream` with the same id are
        also used as neighbours.

        A neighbour whose samples are offset from the trace's sample grid by
        no more than `max_misalignment` samples (for example, day files
        whose start times were shifted by a clock correction) is shifted
        onto the trace's grid, so that its samples are used unchanged.
        If a neighbour is not consecutive (gap or misaligned samples), it is
        either linearly interpolated across the gap
        (`gap_handling='interpolate'` and the gap is no longer than
        `max_interp_gap` samples) or ignored (`gap_handling='edge'`).
        Without a valid neighbour, the start is padded with the first data
        value and the end is not padded, giving edge effects.

        Args:
            stream (:class:`obspy.core.stream.Stream`): input waveforms
            prev_stream (:class:`obspy.core.stream.Stream`): waveforms
                preceding the input (for example the buffer returned by a
                previous call)
            FIRdecim (int): decimation factor associated with FIR to correct
                If None, use self value
            new_loc_code (str or None): Use provided loc code for the output
                traces
            detrend (bool or None): detrend data before conversion.  The
                trend is calculated over the padded data, so results only
                match across neighbouring streams if the data are not
                detrended.  If None, only detrend traces that have no valid
                neighbour on either side
            next_stream (:class:`obspy.core.stream.Stream`): waveforms
                following the input
            gap_handling (str): how to handle non-consecutive neighbours:
                'interpolate' or 'edge'
            max_interp_gap (int): longest gap (in samples) to interpolate.
                If None, the FIR length divided by FIRdecim
            rel_tol (float): tolerance used to calculate n_after (see
                `pad_lengths()`)
            max_misalignment (float): largest offset (in samples) of a
                neighbour from the trace's sample grid for which it is
                shifted onto the grid.  Larger offsets are resampled by
                linear interpolation.
        Returns:
            :class:`obspy.core.stream.Stream`: output waveforms
            :class:`obspy.core.stream.Stream`: buffer to use as prev_stream
                for the next stream
        """
        if detrend not in (True, False, None):
            raise ValueError(f'{detrend=} not in (True, False, None)')
        if gap_handling not in ('interpolate', 'edge'):
            raise ValueError(f'{gap_handling=} not in ("interpolate", "edge")')
        s = stream.copy()
        if FIRdecim is None:
            FIRdecim = self.decimation_factor
        buffer_size = self._buffer_size(FIRdecim)
        n_before, n_after = self.pad_lengths(FIRdecim, rel_tol)
        if max_interp_gap is None:
            max_interp_gap = int(math.ceil(len(self.firzeros) / FIRdecim))
        for tr, tr_orig in zip(s, stream):
            original_dtype = tr.data.dtype
            sr_out = tr.stats.sampling_rate

            # Gather neighbours and add pads
            neighbours = [x for x in stream
                          if x is not tr_orig and x.id == tr.id]
            for st in (prev_stream, next_stream):
                if st is not None:
                    neighbours.extend(st.select(id=tr.id))
            for nb in neighbours:
                if not nb.stats.sampling_rate == sr_out:
                    raise ValueError(f'{tr.id}: trace and neighbour have '
                                     'different sampling rates')
            pre = get_neighbour_pad(tr, neighbours, n_before, 'before',
                                gap_handling, max_interp_gap,
                                max_misalignment)
            post = get_neighbour_pad(tr, neighbours, n_after, 'after',
                                 gap_handling, max_interp_gap,
                                 max_misalignment)
            if detrend is None:
                # Only detrend isolated traces, to keep continuity
                tr_detrend = pre is None and post is None
            else:
                tr_detrend = detrend
            if pre is None:
                pre = np.full(buffer_size, tr.data[0])
            if post is None:
                post = np.array([], dtype=original_dtype)
            data = np.concatenate((pre, tr.data, post))

            # Convert, cut to the input samples and force back to
            # original dtype
            data_corr = self._convert_array(data, len(pre), FIRdecim,
                                            tr_detrend, original_dtype)
            tr.data = np.array(data_corr[:tr.stats.npts],
                               dtype=original_dtype)

            # Adjust starttime if timetag/FIRdecim is not an integer
            offset_samps = self.timetag/FIRdecim
            if offset_samps < buffer_size:
                tr.stats.starttime -= (buffer_size - offset_samps) / sr_out
            elif offset_samps > buffer_size:
                raise ValueError(f'f{offset_samps=} > {buffer_size=}')

            if new_loc_code is not None:
                tr.stats.location = new_loc_code

        # Create a new buffer from the end of the original stream traces
        newbuffer = stream.copy()
        for tr in newbuffer:
            tr.trim(starttime=tr.stats.endtime
                              - (n_before-1)/tr.stats.sampling_rate,
                    nearest_sample=True)
        return s, newbuffer

    def pad_lengths(self, FIRdecim=None, rel_tol=2.**-32):
        """
        Number of neighbouring samples needed to convert data exactly

        The output sample at index j depends on input samples j - n_before
        to infinity, because the conversion filter is an anti-causal IIR
        filter whose impulse response decays as |p|^n, where p is its
        largest pole (the reflected FIR zero closest to the unit circle).
        n_before is therefore exact, whereas n_after is the smallest number
        of samples for which the neglected impulse response sums to no more
        than rel_tol, so that the error is at most rel_tol * max(|data|).
        The default (2^-32) gives an error < 0.5 for any 32-bit integer data.

        Both values are calculated from the impulse response of the full
        conversion (interpolation, filtering and decimation), so they
        adapt to any FIR length and decimation factor.

        Args:
            FIRdecim (int): decimation factor associated with FIR to correct
                If None, use self value
            rel_tol (float): maximum sum of the neglected impulse response
        Returns:
            (tuple):
                n_before (int): number of samples needed before the data
                n_after (int): number of samples needed after the data
        """
        if FIRdecim is None:
            FIRdecim = self.decimation_factor
        buffer_size = self._buffer_size(FIRdecim)
        p_max = np.max(np.abs(np.roots(np.insert(-self.a, 0, 1.))))
        if p_max >= 1:
            raise ValueError(f'conversion filter is unstable ({p_max=})')
        # Probe far enough that the unprobed response is negligible
        m = int(math.ceil(np.log(rel_tol * 1e-6) / np.log(p_max) / FIRdecim)
                ) + buffer_size + 100
        impulse = np.zeros(buffer_size + 2*m + 1)
        impulse[buffer_size + m] = 1.
        h = np.abs(self._convert_array(impulse, buffer_size, FIRdecim,
                                       False, impulse.dtype))
        n_before = max(int(np.nonzero(h)[0].max()) - m, buffer_size)
        # h[m-k] is the weight of the input sample k samples after output
        h_cum = np.cumsum(h)
        i_ok = np.searchsorted(h_cum, rel_tol, side='right') - 1
        if i_ok < 0:
            raise ValueError('impulse response probe is too short')
        n_after = max(m - 1 - i_ok, 0)
        return n_before, n_after

    def _buffer_size(self, FIRdecim):
        """Number of samples that conversion shifts the data"""
        return int(math.ceil(self.timetag/FIRdecim))

    def _convert_array(self, data, n_pre, FIRdecim, detrend, original_dtype):
        """
        Convert a padded 1D array

        Args:
            data (np.ndarray): data, including n_pre samples of padding
                at the start
            n_pre (int): number of padding samples at the start (must be at
                least self._buffer_size(FIRdecim))
            FIRdecim (int): decimation factor associated with FIR to correct
            detrend (bool): detrend data before conversion
            original_dtype (np.dtype): dtype of the input data
        Returns:
            (np.ndarray): converted data, starting at data[n_pre]
        """
        buffer_size = self._buffer_size(FIRdecim)
        if n_pre < buffer_size:
            raise ValueError(f'{n_pre=} < {buffer_size=}')

        # DETREND AND INTERPOLATE
        if detrend is True:
            data = sig.detrend(data, type='linear')
        if original_dtype == np.float32 and data.dtype != np.float32:
            data = np.require(data, dtype=np.float32)
        data_interp = Interpolator(FIRdecim).apply(data, False)

        # CHANGE FIR FROM ZERO PHASE TO MIN PHASE
        data_interp_corr, _ = self._firconvert(data_interp, False)

        # RESAMPLE TO ORIGINAL RATE and remove shift
        data_corr = data_interp_corr[::FIRdecim]
        return data_corr[n_pre - buffer_size: len(data) - buffer_size]

    # ------------------------------------------------------------
    # Data conversion
    # ------------------------------------------------------------
    def apply_SDS(self, SDS_input, SDS_output=None, startdate=None,
                  enddate=None, FIRdecim=None, loc_code='01', verbose=False,
                  max_misalignment=0.5, channels=None):
        """
        Convert an entire SDS directory

        Each file is converted using the previous and next day files of the
        same channel, so that the output is continuous across files (see
        `apply()`).  Neighbouring files are used even if they are outside
        of the requested date range.

        Args:
            SDS_input (str or Path): input SDS filepath
            SDS_output (str, Path or None): output SDS filepath.  If None,
                write to SDS_input
            startdate (str, UTCDateTime or None): limit to files starting
                at or after this date.  If string, must be ISO8601-standard
            enddate (str, UTCDateTime or None): limit to files ending at or
                before this date.  If string, must be ISO8601-standard
            FIRdecim (int): decimation factor associated with FIR to correct
                If None, use self value
            loc_code (str or None): location code to give to output traces.
                If None, keep the input location code
            verbose (bool): print the name of each converted file
            max_misalignment (float): largest offset (in samples) of a
                neighbouring file from a file's sample grid for which it is
                treated as on the grid (see `apply()`)
            channels (list or None): only convert these channels, given as
                (net, sta, loc, chan) tuples.  If None, convert all channels
        Returns:
            (list): (net, sta, loc, chan) of the converted channels
        Raises:
            ValueError: if SDS_output exists and loc_code is None or already
                used in SDS_output
        """
        if loc_code is not None and len(loc_code) > 2:
            raise ValueError(f'{loc_code=} has more than two letters')

        SDS_input = Path(SDS_input)
        SDS_output = SDS_input if SDS_output is None else Path(SDS_output)
        if startdate is not None:
            startdate = UTCDateTime(startdate)
        if enddate is not None:
            enddate = UTCDateTime(enddate)
        if SDS_output.exists():
            self._check_output_valid(SDS_output, loc_code)

        converted = []
        sds_channels = sds_channel_files(SDS_input)
        for key in sorted(sds_channels):
            if channels is not None and key[:4] not in channels:
                continue
            for m, prev_st, st, next_st in sds_neighbour_streams(
                    sds_channels[key], startdate, enddate, verbose):
                st_corr, _ = self.apply(st, prev_st, FIRdecim=FIRdecim,
                                        new_loc_code=loc_code,
                                        next_stream=next_st,
                                        max_misalignment=max_misalignment)
                # Write to output file, with the new location code
                out_fpath = sds_file_path(SDS_output, m, loc=loc_code)
                out_fpath.parent.mkdir(parents=True, exist_ok=True)
                st_corr.write(str(out_fpath), format='MSEED')  # Retains encoding
                if key[:4] not in converted:
                    converted.append(key[:4])
        return converted

    def update_inventory(self, inv, channels, loc_code='01',
                         inv_force_overwrite=False):
        """
        Add the converted channels to an inventory

        Each new channel is a copy of its source channel, with the
        converter's FIR stage replaced by its minimum phase equivalent
        (see `convert_response_stage()`) and the location code changed.

        Args:
            inv (:class:`obspy.core.inventory.Inventory`): inventory
            channels (list): source channels, as (net, sta, loc, chan) tuples
                (for example, the output of `apply_SDS()`)
            loc_code (str): location code of the converted channels
            inv_force_overwrite (bool): If the new channel's seed_id already
                exists, overwrite it
        Returns:
            (:class:`obspy.core.inventory.Inventory`): new inventory
        Raises:
            ValueError: if a channel is not in the inventory or has no FIR
                stage matching the converter
        """
        inv = inv.copy()  # Don't overwrite input inventory
        for net, sta, loc, chan in channels:
            new_cha = self.convert_channel(
                find_channel(inv, net, sta, loc, chan), loc_code)
            inv, _ = insert_channel(inv, net, sta, new_cha, loc, chan,
                                    force_overwrite=inv_force_overwrite)
        return inv

    def convert_channel(self, channel, loc_code='01'):
        """
        Return a copy of an inventory channel for the converted data

        Args:
            channel (:class:`obspy.core.inventory.Channel`): source channel
            loc_code (str): location code of the converted channel
        Returns:
            (:class:`obspy.core.inventory.Channel`)
        Raises:
            ValueError: if no FIR stage of the channel matches the converter
        """
        new_cha = copy.deepcopy(channel)
        stages = new_cha.response.response_stages
        for i, stage in enumerate(stages):
            if not isinstance(stage, FIRResponseStage):
                continue
            try:
                stages[i] = self.convert_response_stage(stage)
                break
            except ValueError:
                continue
        else:
            raise ValueError(f'No FIR stage of {channel.code} matches the '
                             'converter FIR')
        new_cha.location_code = loc_code
        return new_cha

    def matches_channel(self, channel):
        """True if one of the channel's FIR stages is the converter FIR"""
        try:
            self.convert_channel(channel)
        except ValueError:
            return False
        return True

    @staticmethod
    def _check_output_valid(SDS_output, loc_code):
        """
        Verify that a location code is specified and not in output_dir
        """
        # ERROR if no loc_code specified
        if loc_code is None:
            raise ValueError("No location code specified when writing to "
                             "an existing SDS directory")
        # ERROR if loc_code already exists in output_dir
        client = SDS_Client(str(SDS_output))
        bad_ids = []
        for nslc in client.get_all_nslc():
            if nslc[2] == loc_code:
                bad_ids.append('.'.join(nslc))
        if not len(bad_ids) == 0:
            raise ValueError(f"{loc_code=} already exists in {SDS_output=}: {bad_ids}")
        return

    def _firconvert(self, intrace, returnOrigDType=True, plot=False):
        """
        Implements FIR filter conversion
        Args:
            intrace (:class:`numpy.ndarray`): 1d float array
            returnOrigDtype (bool): return output trace as same type as input
                trace (True) or as float double (False)

        Returns:
            (tuple):
                outtrace (:class:`numpy.ndarray`): causal version of trace
                timetag (int): number of samples that outtrace delta response
                    is BEFORE intrace delta response
        """
        if plot is True:
            fig, ax = plt.subplots()
            ax.plot(self.a, color='r', label='a')
            ax.plot(self.b, color='b', label='b')
            plt.legend()
            plt.show()
        original_dtype = intrace.dtype
        # Convert to double for manipulations
        x = np.array(intrace, dtype="double")

        x = x[::-1]  # flip in time
        # because lfilter subtracts "a"s, multiply them by -1
        ainv = -self.a.copy()
        # Insert 1. at beginning of array.  Needed to match results of
        # Scherbaum, even though his code inserts 0.
        ainv = np.insert(ainv, 0, 1.)

        y = sig.lfilter(self.b, ainv, x)
        y = y[::-1]  # flip back in time

        if returnOrigDType:
            y = np.array(y, dtype=original_dtype)

        return y, self.timetag

    # ------------------------------------------------------------
    # Response metadata
    # ------------------------------------------------------------
    def output_stage_values(self, FIRdecim=None):
        """
        FIR stage values describing the converted (minimum phase) data

        `apply()` removes the time-tag shift, so the converted data are the
        signal filtered by the causal minimum-phase FIR, starting at each
        sample's time.  The FIR stage of the output channels should
        therefore have the minimum-phase coefficients, no symmetry and no
        delay or correction (whereas the input channels' linear-phase stage
        has delay = correction = timetag / FIR sampling rate).

        Args:
            FIRdecim (int): decimation factor associated with the FIR.
                If None, use self value
        Returns:
            (dict): 'coefficients' (list, normalized to sum to 1),
                'symmetry' ('NONE'), 'delay' (0.), 'correction' (0.),
                'decimation_factor' and 'input_sample_rate' (self.sps,
                may be None)
        """
        if not self.input_type == 'linear':
            raise ValueError('Output stage values are only defined for '
                             'linear phase input FIRs')
        if FIRdecim is None:
            FIRdecim = self.decimation_factor
        return {'coefficients': self.equivalent_zeros.tolist(),
                'symmetry': 'NONE',
                'delay': 0.,
                'correction': 0.,
                'decimation_factor': FIRdecim,
                'input_sample_rate': self.sps}

    def convert_response_stage(self, stage):
        """
        Return a copy of an input channel's FIR stage for the converted data

        Replaces the coefficients by their minimum phase equivalent (scaled
        to the same sum as the input coefficients, so that the stage gain is
        unchanged), the symmetry by 'NONE' and the decimation delay and
        correction by 0 (see `output_stage_values()`).

        Args:
            stage (:class:`obspy.core.inventory.response.FIRResponseStage`):
                linear phase FIR stage of the input channel
        Returns:
            (:class:`obspy.core.inventory.response.FIRResponseStage`)
        Raises:
            ValueError: if the stage's coefficients are not this converter's
                input FIR
        """
        values = self.output_stage_values()
        coeffs = np.array([float(x) for x in stage.coefficients])
        if stage.symmetry == 'EVEN':
            coeffs = np.concatenate((coeffs, coeffs[::-1]))
        elif stage.symmetry == 'ODD':
            coeffs = np.concatenate((coeffs, coeffs[-2::-1]))
        firzeros = self.firzeros / self.firzeros.sum()
        if not (len(coeffs) == len(firzeros)
                and np.allclose(coeffs / coeffs.sum(), firzeros,
                                rtol=0, atol=1e-4 * np.abs(firzeros).max())):
            raise ValueError(f'{stage.name=} coefficients do not match the '
                             'converter FIR')
        new_stage = copy.deepcopy(stage)
        new_stage.coefficients = (np.array(values['coefficients'])
                                  * coeffs.sum()).tolist()
        new_stage.symmetry = values['symmetry']
        # obspy needs these as FloatWithUncertaintiesAndUnit to write them
        new_stage.decimation_delay = FloatWithUncertaintiesAndUnit(
            values['delay'], unit='SECONDS')
        new_stage.decimation_correction = FloatWithUncertaintiesAndUnit(
            values['correction'], unit='SECONDS')
        if new_stage.name is not None:
            new_stage.name += ' (minimum phase)'
        return new_stage

    # ------------------------------------------------------------
    # File writing
    # ------------------------------------------------------------
    def write_prt(self, filename: str) -> None:
        """
        Write the conversion filter to a .prt file.

        Parameters
        ----------
        filename : str
            Output file path.
        """
        PRTWriter.write(
            filename=filename,
            firzeros=self.firzeros,
            sps=self.sps,
            timetag=self.timetag,
            a=self.a,
            b=self.b
        )

    def write_json(self, filename):
        """
        Make a file of coefficients for converting from zero- to minimum phase

        Args:
            filename: the name of the output file
        """
        JSONWriter.write(filename, self.firzeros, self.sps, self.timetag,
                         self.input_type, self.decimation_factor,
                         self.a, self.b)


class ZPlanePlotter:
    """Utility for plotting z-plane roots."""

    @staticmethod
    def plot(z: np.ndarray, z_max, z_min, z_UC, title_head) -> None:
        """
        Plot complex roots on the z-plane.

        Args:
            z (np.ndarray):  Complex roots (all of them).
            z_max (np.ndarray):  Maximum phase roots
            z_min (np.ndarray):  Minimium phase roots
            z_UC (np.ndarray):  Unit circle roots.
            title_head (str): Text to prepend title
        """
        plt.plot(z.real, z.imag, "ko")
        plt.plot(z_max.real, z_max.imag, "bo", label="maximum phase")
        plt.plot(z_min.real, z_min.imag, "ro", label="minimum phase")
        plt.plot(z_UC.real, z_UC.imag, "go", label="unit circle")
        # Plot unit circle
        theta = np.linspace(0, 2*np.pi, 200)
        plt.plot(np.cos(theta), np.sin(theta), ":")
        plt.legend()
        plt.axis("equal")
        plt.xlabel("Real")
        plt.ylabel("Imag")
        plt.title(f"{title_head} Z-plane roots")
        plt.show()


class ImpulsePlotter:
    """Utility for plotting impulse responses of root subsets."""

    @staticmethod
    def plot(forig: ArrayLike, z_total: np.ndarray,
             z_max: np.ndarray, z_min: np.ndarray, z_UC: np.ndarray,
             f_equiv: np.ndarray, figsize=(6,8)) -> None:
        """
        Plot impulse responses for original, total, max-phase,
        min-phase, and unit-circle components.

        Args:
            forig (ArrayLike):  Original FIR coefficients.
            z_total (np.ndarray): roots of FIR coefficients
            z_max (np.ndarray): roots inside unit circle
            z_min (np.ndarray): roots outside unit circle
            z_UC (np.ndarray): roots on unit circle
            f_equiv (np.ndarray): Equivalent FIR coefficients
            figsize (tuple): figure (x, y) size in inches
        """
        ftotal = np.poly(z_total)
        fmax = np.poly(z_max)
        fmin = np.poly(z_min)
        fUC = np.poly(z_UC)

        titles = [
            "original FIR response",
            "maximum phase part",
            "minimum phase part",
            "unit circle part",
            "calculated FIR response",
            "equivalent minimum phase response"
        ]
        data = [forig, fmax, fmin, fUC, ftotal, f_equiv]

        plt.figure(figsize=figsize)
        i, j = 1, 1
        for (title, arr) in zip(titles, data):
            if i == 2:
                plt.subplot(5, 2, 2+j)
                if j == 1:
                    j = 2
                else:
                    i += 1
            else:
                plt.subplot(5, 1, i)
                i += 1
            plt.plot(arr, "-o")
            plt.title(title)

        plt.tight_layout()
        plt.show()


class PRTWriter:
    """Utility for writing .prt conversion files."""

    @staticmethod
    def write(
        filename: str,
        firzeros: np.ndarray,
        sps: float,
        timetag: float,
        a: np.ndarray,
        b: np.ndarray
    ) -> None:
        """
        Write FIR conversion data to a .prt file.

        Parameters
        ----------
        filename : str
            Output file path.
        firzeros : np.ndarray
            Original FIR zeros.
        sps : float or None
            FIR sampling rate (informational only, written as 0 if None).
        timetag : float
            Time-tag correction.
        a : np.ndarray
            AR coefficients.
        b : np.ndarray
            MA coefficients.
        """
        with open(filename, "w") as fid:
            fid.write("#METHOD: POLYNOMIAL ROOTING\n")
            fid.write("#FIRNAME file name effective FIR filter coefficients\n")
            fid.write("test.001\n")
            fid.write("#NO_EFF:\n")
            fid.write(f"{len(firzeros)}\n")
            fid.write("#FDIG_EFF:\n")
            fid.write(f"{0 if sps is None else sps}\n")
            fid.write("#FIRDATA_EFF:\n")
            fid.write(" ".join(f"{x:19.12g}" for x in firzeros) + "\n\n")

            fid.write("#CORRECTION FILTER:\n")
            fid.write("#CORR_AR:\n")
            fid.write(f"{len(a)}\n")
            fid.write("#CORR_AR_DATA:\n")
            fid.write(" ".join(f"{x:19.12g}" for x in a) + "\n\n")

            fid.write("#CORR_MA:\n")
            fid.write(f"{len(b)}\n")
            fid.write("#CORR_MA_DATA:\n")
            fid.write(" ".join(f"{x:19.12g}" for x in b) + "\n\n")

            fid.write("#TIMETAG:\n")
            fid.write(f"{timetag}\n")


class JSONWriter:
    """Utility for writing JSON conversion files."""

    @staticmethod
    def write(filename: str, firzeros: np.ndarray, sps: float,
              timetag: float, input_type: str, decimation_factor: int,
              a: np.ndarray, b: np.ndarray) -> None:
        """
        Write FIR conversion data to a JSON file.

        Args:
            filename (str): Output file path.
            firzeros (np.ndarray): Original FIR zeros.
            sps (float or None): FIR sampling rate (informational only)
            timetag (float): Time-tag correction.
            input_type (str): input FIR type
            decimation_factor (int): the decimation factor associated with
                the FIR
            a (np.ndarray): AR coefficients.
            b (np.ndarray):  MA coefficients.
        """
        mydict = {'method': 'POLYNOMIAL ROOTING',
                  'firname': "test.001",
                  'fdig': "0",
                  'input_type': input_type,
                  'decimation_factor': decimation_factor,
                  'firdata_eff': firzeros.tolist(),
                  'ar': a.tolist(),
                  'ma': b.tolist(),
                  'sps': sps,
                  'timetag': timetag
                  }
        with open(filename, "w") as fid:
            json.dump(mydict, fid, indent=4)


class Interpolator():
    def __init__(self, ipol_fac):
        """
        Args:
            ipol_fac (int): factor to interpolate by.  Must be a multiple
                of 2, 3 and/or 5
        """
        # Validate value
        conv_fac = ipol_fac
        for mult in (2, 3, 5):
            while conv_fac % mult == 0:
                conv_fac /= mult
        if conv_fac != 1:
            raise ValueError(f"Factor {ipol_fac:d} cannot be separated "
                             "into factors 2, 3, and 5\n")

        self.ipol_fac = ipol_fac

    def apply(self, trace, returnOrigDType=True):
        """
        Interpolate data by an integer multiple of 2, 3, and/or 5

        Args:
            trace (:class:`numpy.ndarray`): data to interpolate
            returnOrigDType (bool): return output trace as same type as input
                trace (True) or as float double (False)
        """

        if self.ipol_fac == 1:
            print("Nothing to interpolate, exit...\n")
            return trace

        conv_fac = int(self.ipol_fac)
        # number of divisions by 2
        pow2 = 2
        no_it2 = 0
        while conv_fac / pow2 == int(conv_fac / pow2):
            pow2 *= 2
            no_it2 += 1
        if no_it2 > 0:
            conv_fac = 2 * conv_fac / pow2

        # number of divisions by 3 */
        pow3 = 3
        no_it3 = 0
        while conv_fac / pow3 == int(conv_fac / pow3):
            pow3 *= 3
            no_it3 += 1
        if no_it3 > 0:
            conv_fac = 3 * conv_fac / pow3

        #  number of divisions by 5 */
        pow5 = 5
        no_it5 = 0
        while conv_fac / pow5 == int(conv_fac / pow5):
            pow5 *= 5
            no_it5 += 1
        if no_it5 > 0:
            conv_fac = 5 * conv_fac / pow5

        if conv_fac != 1:
            raise ValueError(f"Factor {self.ipol_fac:d} cannot be separated "
                             "into factors 2, 3, and 5\n")

        original_dtype = trace.dtype
        # interpolations by factors of  2
        for k in range(0, no_it2):  # (k=0;k<no_it2;k++) :
            trace = self._ipol2(trace)
        # interpolations by factors of 3
        for k in range(0, no_it3):  # (k=0;k<no_it3;k++) :
            trace = self._ipol3(trace)
        # interpolations by factors of  5
        for k in range(0, no_it5):  # (k=0;k<no_it5;k++)
            trace = self._ipol5(trace)

        if returnOrigDType:
            # Force back to original data type:
            return np.array(trace, dtype=original_dtype)
        else:
            return np.array(trace)

    # E. Wielandt's filter coefficients for interpolation ratios 2, 3 and 5 */

    @staticmethod
    def _ipol2(xin, debug=False):
        "interpolation by factor 2"

        g2_1 = np.array(
            [
                -0.0002616,
                0.0009302,
                -0.0023258,
                0.0049142,
                -0.0092537,
                0.0161355,
                -0.0266150,
                0.0424554,
                -0.0670612,
                0.1091686,
                -0.2008408,
                0.6327541,
                0.6327542,
                -0.2008408,
                0.1091686,
                -0.0670612,
                0.0424554,
                -0.0266150,
                0.0161355,
                -0.0092537,
                0.0049142,
                -0.0023258,
                0.0009302,
                -0.0002616,
            ],
            dtype="double",
        )

        # xx = interpolated sample at center between index i and i+1
        xx = np.convolve(xin, g2_1[::-1], mode="full")
        if debug:
            print(f"len(g2_1)={len(g2_1)}")
        xx = xx[int(len(g2_1) / 2): (len(xin) + int(len(g2_1) / 2))]
        # print("{:d} {:d}".format(len(xin),len(xx)))
        # WEAVE data together
        xout = np.column_stack((xin, xx))  # 2D
        xout = xout.flatten()  # To 1D
        xout = xout[0:-1]  # get rid of extrapolations
        # Should verify that xout[0]==xin[0] and that xout[-1]==xin[-1]
        return xout

    @staticmethod
    def _ipol3(xin):
        """interpolation by factor 3"""
        ifac = 3  # interpolation factor
        g3_1 = np.array(
            [
                -0.0002324,
                0.0008256,
                -0.0020654,
                0.0043684,
                -0.0082410,
                0.0144087,
                -0.0238643,
                0.0383080,
                -0.0611438,
                0.1015046,
                -0.1959616,
                0.8226883,
                0.4111037,
                -0.1564917,
                0.0885536,
                -0.0553522,
                0.0353777,
                -0.0223077,
                0.0135759,
                -0.0078056,
                0.0041522,
                -0.0019670,
                0.0007871,
                -0.0002211,
            ],
            dtype="double",
        )

        xx1 = np.convolve(xin, g3_1[::-1], mode="full")
        xx2 = np.convolve(xin, g3_1, mode="full")
        xx1 = xx1[len(g3_1) / 2: (len(xin) + len(g3_1) / 2)]
        xx2 = xx2[len(g3_1) / 2: (len(xin) + len(g3_1) / 2)]
        # WEAVE xx1 and xx2 into xin
        xout = np.column_stack((xin, xx1, xx2))
        xout = xout.flatten()
        xout = xout[0: (-ifac + 1)]  # get rid of extrapolations
        # Should verify that xout[0]==xin[0] and that xout[-1]==xin[-1]
        return xout

    @staticmethod
    def _ipol5(xin):
        """interpolation by factor 5"""
        ifac = 5  # interpolation factor
        g5_1 = np.array(
            [
                -0.0001611,
                0.0005720,
                -0.0014316,
                0.0030307,
                -0.0057263,
                0.0100352,
                -0.0166788,
                0.0269186,
                -0.0433567,
                0.0732492,
                -0.1480766,
                0.9320452,
                0.2327664,
                -0.0984035,
                0.0572469,
                -0.0362360,
                0.0233232,
                -0.0147709,
                0.0090151,
                -0.0051932,
                0.0027660,
                -0.0013112,
                0.0005248,
                -0.0001473,
            ],
            dtype="double",
        )
        g5_2 = np.array(
            [
                -0.0002526,
                0.0008977,
                -0.0022452,
                0.0047467,
                -0.0089479,
                0.0156272,
                -0.0258392,
                0.0413715,
                -0.0657512,
                0.1082703,
                -0.2048060,
                0.7525200,
                0.5015040,
                -0.1790148,
                0.0997642,
                -0.0619420,
                0.0394431,
                -0.0248145,
                0.0150789,
                -0.0086611,
                0.0046043,
                -0.0021804,
                0.0008723,
                -0.0002452,
            ],
            dtype="double",
        )

        xx1 = np.convolve(xin, g5_1[::-1], mode="full")
        xx2 = np.convolve(xin, g5_2[::-1], mode="full")
        xx3 = np.convolve(xin, g5_2, mode="full")
        xx4 = np.convolve(xin, g5_1, mode="full")
        xx1 = xx1[len(g5_1) / 2: (len(xin) + len(g5_1) / 2)]
        xx2 = xx2[len(g5_2) / 2: (len(xin) + len(g5_2) / 2)]
        xx3 = xx3[len(g5_2) / 2: (len(xin) + len(g5_2) / 2)]
        xx4 = xx4[len(g5_1) / 2: (len(xin) + len(g5_1) / 2)]

        # WEAVE data together
        xout = np.column_stack((xin, xx1, xx2, xx3, xx4))
        xout = xout.flatten()
        xout = xout[0: (-ifac + 1)]  # get rid of extrapolations
        # Should verify that xout[0]==xin[0] and that xout[-1]==xin[-1]
        return xout
