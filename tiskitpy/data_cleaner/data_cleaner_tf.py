"""
Clean data using frequency response functions between channels

Test using spectra (ATACR-style) and time-series (TisKit style):  is there a
difference?  If not, will be a lot faster to remove noise after calculating
spectra!
"""
from copy import deepcopy

import numpy as np
from scipy import signal
from obspy.core.stream import Stream
from matplotlib import pyplot as plt

from ..response_functions import ResponseFunctions
from ..spectral_density import SpectralDensity
from .rf_list import RFList
from ..cleaned_stream import CleanedStream
from ..utils import CleanSequence as CS, stream_synchronize, stream_unmask
from tiskitpy.logger import init_logger

logger = init_logger()


class DataCleaner:
    """
    Calculate and apply ResponseFunction-based data cleaner

    Works on raw data, without instrument corrections

    By default, calculates corrected intermediate and final spectra directly
    from the stream data, applying the frequency response functions to each FFT.

    """

    def __init__(self, stream, remove_list, noise_channel="output",
                 n_to_reject=3, min_freq=None, max_freq=None,
                 show_progress=False, fast_calc=False, **kwargs):
        """
        Args:
            stream (:class:`obspy.core.stream.Stream`): time series data used
                to calculate cleaner
            remove_list (list): list of channels to remove, in order.  Can use
                "*" and "?" as wildcards
            noise_channel (str): which channel the noise is on.  Choices are:
                "input", "output", "equal", "unknown", "model"
            n_to_reject (int): Number of neighboring frequencies for which
                the coherence must be above the 95% signifcance level in order
                to use for cleaning (0 means use all frequencies)
            min_freq (float): Do not process frequencies below this value
            max_freq (float): Do not process frequencies above this value
            show_progress (bool): Show progress plots
            fast_calc (bool): Calculate corrected spectra directly from
                previous spectra (ATACR-style).
            kwargs (:class:`SpectralDensity.from_stream()` properties):
                 window_s, windowtype, z_threshold, avoid_spans, ...
        Attributes:
            RFList (list of :class:`ResponseFunctions`): list of data cleaner
                frequency response functions
            starttimes (list of :class:`obspy.UTCDateTime`): start times for
                each spectra used
        """
        if not isinstance(stream, Stream):
            raise ValueError("stream is a {type(stream)}, not an obspy Stream")
        stream = stream.copy()   # Do not change values in original stream
        # Remove any instrument response from stream
        for tr in stream:
            tr.stats.response = None
        sdfs = [SpectralDensity.from_stream(stream, **kwargs)]
        # Make sure we're using the same starttimes throughout the process
        self.starttimes = sdfs[0].starttimes
        self.RFList = RFList()
        # Calculate removal frequency response functions
        # out_chan_ids = sdfs[0].ids
        clean_sequence = []
        # in_ids = [sdfs[0].channel_id(x) for x in remove_list]
        # for in_id in in_ids:
        removed_codes=[]
        for remove_code in remove_list:
            in_chan_id = sdfs[-1].channel_id(remove_code)
            removed_ids = [sdfs[-1].channel_id(x) for x in removed_codes]
            out_chan_ids = [CS.tiskitpy_id(x, y)
                       for x, y in zip(sdfs[-1].seed_ids, sdfs[-1].clean_sequences)
                       if not x == CS.seed_id(in_chan_id) and not x in [CS.seed_id(y) for y in removed_ids]]
            rf = ResponseFunctions(
                sdfs[-1],
                in_chan_id,
                out_chan_ids,
                noise_channel,
                n_to_reject,
                min_freq,
                max_freq,
                show_progress=show_progress,
            )
            removed_codes.append(remove_code)
            # Apply data cleaner and update channel names with removed channels
            clean_sequence.append(in_chan_id)
            self.RFList.append(deepcopy(rf))
            if fast_calc:
                new_sdf = self._removeRF_SDF(sdfs[-1], rf, show_progress)
            else:
                new_sdf = SpectralDensity.from_stream(
                    stream, data_cleaner=self, **kwargs)
            sdfs.append(new_sdf)
        if show_progress:
            logger.info("Plotting sdfs after data cleaner applied")
            self._plot_sdfs(sdfs)

    def __str__(self):
        s = "DataCleaner object:\n"
        s += f" {len(self.starttimes)} starttimes\n"
        s += f" response_functions = {self.RFList.__str__()}"
        return s

    def clean_sdf(self, sdf):
        """
        Correct spectral density functions

        Can only apply fast_calc method, because SpectralDensity object
        does not contain enough information to do otherwise

        Args:
            sdf (:class:`.SpectralDensity`): data to apply to
        Return:
            sdf (:class:`.SpectralDensity`): corrected spectral densities
        """
        if not isinstance(sdf, SpectralDensity):
            raise TypeError('sdf is not a SpectralDensity object')
        for rf in self.RFList:
            sdf = self._removeRF_SDF(sdf, rf)
        return sdf

    def clean_stream_to_sdf(self, stream, fast_calc=False, **kwargs):
        """
        Calculate corrected spectral density functions from an input stream

        Applys the DataCleaner values to each FFT

        Args:
            stream (:class:`obspy.core.stream.Stream`): data to apply to
            fast_calc (bool): Calculate corrected spectra directly from
                previous spectra.
            kwargs (dict): keyword arguments for `SpectralDensity.from_stream()`
        Return:
            sdf (:class:`.SpectralDensity`): corrected spectral density
                functions
        """
        assert isinstance(stream, Stream)
        if fast_calc:
            sdf = SpectralDensity.from_stream(stream, **kwargs)
            sdf = self.clean_sdf(sdf)
        else:
            sdf = SpectralDensity.from_stream(stream, data_cleaner=self,
                                              **kwargs)
        return sdf

    def clean_stream(self, stream, in_time_domain=False):
        """
        Calculate corrected data stream

        Args:
            stream (Stream): list of channels to remove, in order
            in_time_domain(bool): do work in time domain (default is freq
                domain, time_domain is much slower)
        Returns:
            stream (:class:`obspy.core.stream.Stream`): corrected data

        Assumes frequency response function is for counts, not resp_corrected data.
        """
        assert isinstance(stream, Stream)
        
        stream = stream_unmask(stream)
        # Make sure that the array is not masked
        # if np.any([np.ma.count_masked(tr.data) for tr in stream]):
        #     logger.warning('Unmasking masked data (usually a gap or overlap)')
        #     stream = stream.split().merge(fill_value='interpolate')

        out_stream = CleanedStream(stream)

        if in_time_domain is True:
            logger.warning("Correcting traces in the time domain: VERY SLOW")
        else:
            logger.info("Correcting traces in the frequency domain")

        for rfs in self.RFList:
            in_id = rfs.input_channel_id
            for out_id in rfs.output_channel_ids:
                # print(f'{out_stream=}, {in_id=}, {out_id=}')
                in_trace = out_stream.select(id=in_id)[0]
                out_trace = out_stream.select(id=out_id)[0]
                out_stream.remove(out_trace)
                out_trace = self._correct_trace(
                    in_trace, out_trace,
                    rfs.freqs, rfs.corrector_wrt_counts(out_id),
                    in_time_domain
                )
                out_trace = CS.tag(out_trace, in_trace.id)
                out_stream += out_trace
        return out_stream

    def plot(self):
        """
        Plot data cleaning frequency response functions
        """
        self.RFList.plot()

    def _correct_trace(self, in_trace, out_trace, f, rf, in_time_domain=False,
                       max_reject_sync=0.01):
        """
        Correct a trace using an input trace and a frequency response function

        Note that the frequency response function should be between the trace's units

        Args:
            in_trace (:class: `obspy.core.trace.Trace`): input trace
            out_trace (:class: `obspy.core.trace.Trace`): original output trace
            f (:class:`numpy.ndarray`): frequencies
            rf (:class:`numpy.ndarray`): frequency response function between the input
                and output traces (counts/count)
            in_time_domain (bool): do correction in time domain
            max_reject_sync (float): max_reject for stream_synchronize()

        Returns:
            out_trace_corrected (:class:`obspy.core.trace.Trace`): corrected
                output trace
        """
        stream = stream_synchronize(Stream([in_trace, out_trace]), max_reject_sync)
        in_trace, out_trace = stream[0], stream[1]
        # self._validate_streams_synchronized(in_trace, out_trace)
        # in_trace = in_trace.copy()
        # out_trace = out_trace.copy()
        in_trace.detrend("linear")
        out_trace.detrend("linear")
        if in_time_domain:
            return self._correct_trace_time(in_trace, out_trace, f, rf)
        else:
            return self._correct_trace_freq(in_trace, out_trace, f, rf)

    @staticmethod
    def _correct_trace_time(in_trace, out_trace, f, rf):
        """
        Correct trace in the time domain

        Args:
            in_trace (:class:`obspy.core.trace.Trace`): input (noise) trace
            out_trace (:class:`obspy.core.trace.Trace`): trace to correct
            f (:class:`numpy.ndarray`): frequencies
            rf (:class:`numpy.ndarray`): out_trace/in_trace frequency response function
        """
        if not 2 * f[-1] == in_trace.stats.sampling_rate:
            raise ValueError(
                "different rf ({}) & trace ({}) sample rates".format(
                    2 * f[-1], in_trace.stats.sampling_rate
                )
            )
        # Time domain transform of frequency response function
        # Why don't we need to use the complex conjugate of rf?
        irf = np.fft.ifftshift(np.fft.irfft(rf))
        out_trace_corr = out_trace.copy()
        corr_data = np.convolve(in_trace.data.copy(), irf, mode="same")
        out_trace_corr.data -= signal.detrend(corr_data)
        return out_trace_corr

    @staticmethod
    def _correct_trace_freq(in_trace, out_trace, f, rf):
        """
        Correct trace in the frequency domain

        Args:
            in_trace (:class:`obspy.core.trace.Trace`): input (noise) trace
            out_trace (:class:`obspy.core.trace.Trace`): trace to correct
            f (:class:`numpy.ndarray`): frequencies
            rf (:class:`numpy.ndarray`): out_trace/in_trace frequency response function
        """
        npts = in_trace.stats.npts
        trace_sr = in_trace.stats.sampling_rate
        fft_len = 2 ** int(np.ceil(npts)).bit_length()
        buff = np.zeros(fft_len - npts)
        # Transform all to frequency domain
        in_rfft = np.fft.rfft(np.concatenate((in_trace.data, buff)))
        out_rfft = np.fft.rfft(np.concatenate((out_trace.data, buff)))
        f_rfft = np.fft.rfftfreq(fft_len, 1.0 / trace_sr)
        rf_interp = np.interp(f_rfft, f, rf)
        # Subract coherent part and convert to time domain
        # Why isn't the complex conjugate of the rf used?
        new_out = np.fft.irfft(out_rfft - in_rfft * rf_interp)
        # new_out = np.fft.irfft(out_rfft - in_rfft * np.conj(rf_interp))
        # Stuff into new trace
        out_trace_corr = out_trace.copy()
        out_trace_corr.data = new_out[:npts]
        return out_trace_corr

    @staticmethod
    def _validate_streams_synchronized(in_trace, out_trace):
        its, ots = in_trace.stats, out_trace.stats
        errhdr = "in_trace, out_trace have different "
        if its.npts != ots.npts:
            raise ValueError(errhdr + f"lengths ({its.npts}, {ots.npts})")
        if its.sampling_rate != ots.sampling_rate:
            raise ValueError(
                errhdr
                + "samp rates ({}, {})".format(
                    its.sampling_rate, ots.sampling_rate
                )
            )
        if abs(its.starttime - ots.starttime) > 1 / its.sampling_rate:
            raise ValueError(
                errhdr
                + "start times ({}, {})".format(its.starttime, ots.starttime)
            )

    def _removeRF_SDF(self, sdf, rf, show_progress=False):
        """
        Directly clean a SpectralDensity object using frequency response functions

        This is the ATACR method.  It works on the overall spectra, not
        individual windows.

        Used for the "fast_calc" method and for data already in
        SpectralDensity format

        Args:
            sdf (:class:`.SpectralDensity`): one-sided spectral density
                functions
            rf (:class:`.ResponseFunctions`): Data Cleaner Frequency Response Functions
                w.r.t. one input channel
            show_progress (bool): show progress plots for each step

        Returns:
            sdf (:class:`.SpectralDensity`): cleaned spectral density functions
        """
        ic = rf.input_channel_id
        out_ids = rf.output_channel_ids

        sdf = deepcopy(sdf)  # Avoids overwriting input object
        in_auto = sdf.autospect(ic)
        for oc in out_ids:
            rf_oc = rf.corrector(oc)
            out_auto = sdf.autospect(oc)
            crossspect = sdf.crossspect(ic, oc)
            if show_progress:
                self._plot_removeRF_SDF(
                    rf, sdf, in_auto, out_auto, rf_oc, ic, oc
                )
            # Bendat & Piersol 1986 eqs 6.35 and 6.41
            sdf.put_autospect(oc, out_auto - in_auto * np.abs(rf_oc) ** 2)
            # Bendat & Piersol eq6.36 with
            sdf.put_crossspect(ic, oc, crossspect - in_auto * rf_oc)
            # Add to clean_sequence:
            # sdf.put_clean_sequence(oc, CS.tag(sdf.clean_sequence(oc), ic))
            new_clean_sequence = sdf.clean_sequence(oc) + [ic]
            sdf.put_clean_sequence(oc, new_clean_sequence)
            sdf.replace_channel_id(oc, CS.tiskitpy_id(CS.seed_id(oc),
                                                      new_clean_sequence))
        return sdf

    @staticmethod
    def _plot_removeRF_SDF(rf, sdf, in_auto, out_auto, rf_oc, ic, oc):
        """Plot elements of frequency response function removal using SDFs"""
        f = rf.freqs
        fig, ax = plt.subplots(1, 1)
        ax2 = ax.twinx()
        ax2.semilogx(f, np.abs(sdf.coherence(ic, oc)),
                     "k--", alpha=0.5, label="coherence")
        ax2.axhline(np.abs(sdf.coh_signif(0.95)), color="k", ls=":",
                    alpha=0.5, label="coherence significance level")
        ax.loglog(f, np.abs(out_auto), label=oc)
        ax.loglog(f, np.abs(in_auto), alpha=0.5, label=ic)
        # ax.loglog(f, np.abs(rf_oc), label='rf^2')
        ax.loglog(f, np.abs(in_auto * np.abs(rf_oc) ** 2), alpha=0.5,
                  label=ic + "*rf^2")
        ax.loglog(f, np.abs(out_auto - in_auto * np.abs(rf_oc) ** 2),
                  alpha=0.5, label=f"{oc} - ({ic}*rf^2)")
        ax.set_title(f"in={ic}, out={oc}")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Autospectral Density (counts^2/Hz)")
        ax.legend()
        ax2.set_ylabel("Coherence")
        ax2.set_ylim(0, 1)
        plt.show()

    def _plot_sdfs(self, sdfs):
        """
        Plot comparison of SpectralDensities for each channel

        sdfs (list of :class:`SpectralDensity`): spectral density functions
            for each step of input channel removal
        dont_plot (list of str): list of channel not to plot
        """
        if len(sdfs) == len(self.RFList) + 1:
            ValueError("len(sdfs) isn't len(self.RFList)+1")
        inputs = [""] + [t.rfs.input_channel for t in self.RFList]
        n_channels = len(sdfs[0].ids)
        rows = int(np.floor(np.sqrt(n_channels)))
        cols = int(np.ceil(n_channels / rows))
        fig, axs = plt.subplots(rows, cols)
        irow, icol = 0, 0
        for chan_id in sdfs[0].ids:
            ax = axs[irow, icol]
            plot_sdfs = sdfs
            if id in inputs:
                plot_sdfs = sdfs[:inputs.index(chan_id)]
            for sdf in plot_sdfs:
                ax.loglog(sdf.freqs,
                          np.abs(sdf.autospect(chan_id)),
                          alpha=0.75,
                          label=chan_id + sdf.cleaners(chan_id))
            if irow == rows - 1:
                ax.set_xlabel("Frequency (Hz)")
            ax.legend()
            icol += 1
            if icol >= cols:
                irow += 1
                icol = 0
        plt.show()


if __name__ == "__main__":
    import doctest

    doctest.testmod()
