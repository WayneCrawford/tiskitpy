import json
import math
from typing import Tuple
from pathlib import Path
import re

# from pprint import pprint

import numpy as np
from numpy.typing import ArrayLike
# from numpy.fft import fft
import matplotlib.pyplot as plt
import scipy.signal as sig
from obspy.core.stream import read as stream_read
from obspy import UTCDateTime
from obspy.clients.filesystem.sds import Client as SDS_Client

# from ..logger import init_logger

# logger = init_logger()


class FIRConverter:
    """
    Calculate and perform FIR conversion from zero phase to minimum phase

    Parameters:
        firzeros (ArrayLike) Input FIR filter zeros (polynomial coefficients).
        sps (float): Sampling rate in Hz.
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
            sps (float):
            timetag (float)
            decimation_factor (int)
            input_type (str)
        """
        self.firzeros: np.ndarray = np.asarray(firzeros, dtype=float)
        self.sps: float = float(sps)
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
    def from_zeros_file(cls, filename: str, sps: float):
        """
        Read FIR coefficients from a JSON file, and calculate the conversion

        The JSON file must have the elements 'fir', 'timetag',
        'decimation_factor', and 'type'
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
                return ValueError(f'{filename}: Unknown {book("symmetry")=}')
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
        sps = float(book["sps"])
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
    def apply(self, stream, buffer=None, FIRdecim=None, new_loc_code='01'):
        """
        Convert a zero-phase stream to minimum-phase
        (if the object has the correct zero-phase information)

        To keep the same starttime as the input stream, adds
        self.timetag/2 samples to start and cuts off the same number of
        samples at the end.  This is done by adding a self.timetag/2 - sample
        buffer to the data before converting.  The buffer may be provided in
        the function call (it should be data just preceding the current data:
        this will be verified by comparing the buffer's endtime to the stream's
        starttime): otherwise a buffer will be created with the same values as
        the first sample in each stream trace

        Args:
            stream (:class:`obspy.core.stream.Stream`): input waveforms
            buffer (:class:`obspy.core.stream.Stream`): waveform preceding the
                input (must have same trace_ids as stream)
            FIRdecim (int): decimation factor associated with FIR to correct
                If None, use self value
            new_loc_code (str or None): Use provided loc code for the output
                traces
        Returns:
            :class:`obspy.core.stream.Stream`: output waveforms
            :class:`obspy.core.stream.Stream`: buffer for the next stream
        """
        ##################################################

        s = stream.copy()
        newbuffer = stream.copy()
        if FIRdecim is None:
            FIRdecim = self.decimation_factor
        buffer_size = int(math.ceil(self.timetag/FIRdecim))
        for tr in s:
            original_dtype = tr.data.dtype
            sr_out = tr.stats.sampling_rate

            # Extract data (1D numpy array) and add buffer
            data = self._add_buffer(tr, buffer, buffer_size)

            # DETREND AND INTERPOLATE
            data = sig.detrend(data, type='linear')
            if original_dtype == np.float32 and data.dtype != np.float32:
                data = np.require(data, dtype=np.float32)
            data_interp = Interpolator(FIRdecim).apply(data, False)

            # CHANGE FIR FROM ZERO PHASE TO MIN PHASE
            data_interp_corr, offset_npts = self._firconvert(
                data_interp, False)

            # RESAMPLE TO ORIGINAL RATE and remove extras from end
            data_corr = data_interp_corr[::FIRdecim]

            # Force data back to original dtype and put back in trace
            tr.data = np.array(data_corr[:-buffer_size], dtype=original_dtype)

            # Adjust starttime if offset_npts/FIRdecim is not an integer
            offset_samps = offset_npts/FIRdecim
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
                              - (buffer_size+1)/tr.stats.sampling_rate,
                    nearest_sample=True)
        return s, newbuffer

    def _add_buffer(self, trace, buffer, buffer_size):
        """
        Add a buffer to data, so that it won't be shifted after converting
        to minimum phase
        """
        data = trace.data
        if buffer is not None:
            buffer = buffer.copy().select(id=trace.id)
            if len(buffer) == 0:
                buffer = None
            elif len(buffer) > 1:
                raise ValueError(f'more than one buffer trace matches {trace.id=}')
            else:
                buffer_trace = buffer[0]
                sr = trace.stats.sampling_rate
                if not sr == buffer_trace.stats.sampling_rate:
                    raise ValueError('trace and buffer have different sampling rates')
                if not math.abs(trace.stats.starttime
                                - buffer_trace.stats.endtime
                                - 1/sr) < 0.5/sr:
                    print('trace and buffer are not continuous, '
                          'using generic buffer')
                    buffer = None
                elif len(buffer_trace.data) < buffer_size:
                    print(f'{len(buffer_trace.data)=} < {buffer_size=}, '
                          'using generic buffer')
                    buffer = None
                else:
                    buffer_data = buffer_trace.data[:buffer_size]

        if buffer is None:
            buffer_data = np.full(buffer_size, data[0])

        data = np.concatenate((buffer_data, data))
        return data

    # ------------------------------------------------------------
    # Data conversion
    # ------------------------------------------------------------
    def apply_SDS(self, SDS_input, SDS_output, startdate=None, enddate=None,
                  FIRdecim=None, loc_code='01'):
        """
        Convert an entire SDS directory

        UNTESTED!!!!!

        Args:
            SDS_input (str or Path): input SDS filepath
            SDS_output (str or Path): output SDS filepath
            startdate (str, UTCDateTime or None): limit to files starting after
                this date.  If string, must be ISO8601-standard
            enddate (str, UTCDateTime or None): limit to files ending before
                this date.  If string, must be ISO8601-standard
            FIRdecim (int): decimation factor associated with FIR to correct
                If None, use self value
            loc_code (char): location_code to give to output traces
        Raises:
            ValueError: if SDS_output directory exists already
        """
        if len(loc_code) > 2:
            raise ValueError(f'{loc_code=} has more than two letters')

        SDS_input = Path(SDS_input)
        SDS_output = Path(SDS_output)
        if startdate is not None:
            startdate = UTCDateTime(startdate)
        if enddate is not None:
            enddate = UTCDateTime(enddate)
        if SDS_output.exists():
            self._check_output_dir_loc_validity(SDS_output, loc_code)
        year_paths = [x for x in SDS_input.iterdir() if x.is_dir()]
        for year_path in year_paths:
            net_paths = [x for x in year_path.iterdir() if x.is_dir()]
            for net_path in net_paths:
                sta_paths = [x for x in net_path.iterdir() if x.is_dir()]
                for sta_path in sta_paths:
                    chan_paths = [x for x in sta_path.iterdir() if x.is_dir()]
                    for chan_path in chan_paths:
                        buffer = None
                        # Must sort files, to pass buffer for next day
                        file_paths = sorted(chan_path.iterdir(),
                                            key=lambda x: x.name)
                        for file_path in file_paths:
                            buffer = self._read_write_SDS_data(file_path,
                                                               startdate,
                                                               enddate,
                                                               SDS_input,
                                                               SDS_output,
                                                               buffer,
                                                               loc_code)

    def _check_output_valid(self, SDS_output, loc_code):
        """
        Verify that a location code is specified and not in output_dir
        """
        # ERROR if no loc_code specified
        if loc_code == None:
            raise ValueError("No location code specified when writing to "
                             "an existing SDS directory")
        # ERROR if loc_code already exists in output_dir
        client = SDS_Client(SDS_output)
        bad_ids = []
        for nslc in client.get_all_nslc():
            if nslc[2] == loc_code:
                bad_ids.extend('.'.join(nslc))
        if not len(bad_ids) == 0:
            raise ValueError(f"{loc_code=} already exists in {SDS_output=}: {bad_ids}")
        return

    def _read_write_SDS_data(self, file_path, startdate, enddate,
                             SDS_input, SDS_output, buffer, loc_code='01'):
        # VERIFY THAT IT'S AN SDS FILE AND IN THE DATA RANGE
        if not re.match('\d{4}-\d{2}-\d{2}-\d{4}-\d{2}M\..........', file_path):
            return None
        filedate = UTCDateTime(file_path[:18])
        # TEST AGAINST DATA BOUNDS
        if startdate is not None:
            if filedate < startdate:
                return None
        if enddate is not None:
            if filedate > enddate:
                return None
        print(file_path)

        # Read and convert
        st = stream_read(file_path)
        st_corr, buffer = self.apply(st, new_loc_code=loc_code)
        # Create output file and write to it
        out_fname = str(file_path).replace(SDS_input, SDS_output)
        out_fpath = Path(out_fname)
        out_fpath.parent.mkdir(parents=True, exist_ok=True)
        st_corr.write(out_fname, format='MSEED')  # Retains encoding
        return buffer

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
        sps : float
            Sampling rate.
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
            fid.write(f"{sps}\n")
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
            sps (float): Sampling rate.
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
