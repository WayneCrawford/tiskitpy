"""
Clean data using transfer functions between channels

Test using spectra (ATACR-style) and time-series (TisKit style):  is there a
difference?  If not, will be a lot faster to remove noise after calculating
spectra!
"""
import fnmatch
from copy import deepcopy
import logging

import numpy as np
from scipy import signal
from obspy.core.stream import Stream
from matplotlib import pyplot as plt

from ..transfer_functions import TransferFunctions
from ..spectral_density import SpectralDensity
from .dctfs import DCTF, DCTFs
from .cleaner_string import CleanerString as CS


class DataCleaner:
    """
    Calculate and applying Transfer_Function-based data cleaner

    Works on raw data, without instrument corrections

    By default, calculates corrected intermediate and final spectra directly
    from the stream data, applying the transfer functions to each FFT.

    Note:
        Channel names are updated with each removal of noise correlated
        with another channel, in the "location" slot of the seed_id
        (NET.STA.LOC.CHAN) by adding '-x' to LOC, where 'x' is the shortest
        unique string from the end of the input channel's name.  If the channel
        name has no '.', one is added at the beginning to allow the cleaning
        codes to be added before

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
    """

    def __init__(self, stream, remove_list, noise_channel="output",
                 n_to_reject=3, min_freq=None, max_freq=None,
                 show_progress=False, fast_calc=False, **kwargs):
        """
        Attributes:
            DCTFs (list of :class:`DCTF`): list of data cleaner transfer
                functions
            starttimes (list of :class:`obspy.UTCDateTime`): start times for
                each spectra used
        """
        if not isinstance(stream, Stream):
            raise ValueError("stream is a {type(stream)}, not an obspy Stream")
        stream = stream.copy()   # Do not change values in original stream
        # Remove any response from stream
        for tr in stream:
            tr.stats.response = None
        sdfs = [SpectralDensity.from_stream(stream, **kwargs)]
        # Make sure we're using the same starttimes throughout the process
        self.starttimes = sdfs[0].starttimes
        self.DCTFs = DCTFs()
        # Calculate removal transfer functions
        out_chans = sdfs[0].channel_names
        # in_list = [_list_match_pattern(x, out_chans) for x in remove_list]
        in_list = [sdfs[0].channel_name(x) for x in remove_list]
        remove_seq, remove_new = "", ""
        for in_chan in in_list:
            # ic = in_chan + remove_seq
            ic = CS.insert(remove_seq, in_chan)
            out_chans = [x for x in out_chans if not x == ic]
            tf = TransferFunctions(
                sdfs[-1],
                ic,
                out_chans,
                noise_channel,
                n_to_reject,
                min_freq,
                max_freq,
                show_progress=show_progress,
            )
            # Apply data cleaner and update channel names with removed channels
            remove_new = CS.make(in_chan, remove_seq)
            remove_seq += remove_new
            new_dctf = DCTF(remove_channel=ic, remove_sequence=remove_seq,
                            tfs=tf)
            self.DCTFs.append(new_dctf)
            if fast_calc:
                sdfs.append(
                    self._removeDCTF_SDF(sdfs[-1], new_dctf, show_progress)
                )
            else:
                sdfs.append(SpectralDensity.from_stream(
                    stream, data_cleaner=self, **kwargs))
            out_chans = [CS.insert(remove_new, x) for x in out_chans]
            # out_chans = [x + remove_new for x in out_chans]
        if show_progress:
            logging.info("Plotting sdfs after data cleaner applied")
            self._plot_sdfs(sdfs)

    def __str__(self):
        s = "DataCleaner object:\n"
        s += "   Input channel      | Output channels\n"
        s += "   ================== | ===============\n"
        for tf in self.DCTFs:
            s += "   {:18s} | {}\n".format(
                tf.tfs.input_channel,
                [CS.insert(tf.remove_sequence, CS.strip(x))
                 for x in tf.tfs.output_channels]
            )
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
        for dctf in self.DCTFs:
            sdf = self._removeDCTF_SDF(sdf, dctf)
        return sdf

    def clean_stream_to_sdf(self, stream, fast_calc=False, **kwargs):
        """
        Calculate corrected spectral density functions from an input stream

        Applys the DataCleaner values to each FFT
        Args:
            stream (:class:`obspy.core.stream.Stream`): data to apply to
            fast_calc (bool): Calculate corrected spectra directly from
                previous spectra.
            **kwargs (dict): keyword arguments for
                SpectralDensity.from_stream()
        Return:
            sdf (:class:`.SpectralDensity`): corrected spectral density
                functions
        """
        assert isinstance(stream, Stream)
        if fast_calc:
            sdf = SpectralDensity.from_stream(stream, **kwargs)
            sdf = self.clean_sdf(sdf)
        else:
            # sdf = SpectralDensity.from_stream(stream, dctfs=self.DCTFs,
            #                                   **kwargs)
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

        Assumes transfer function is for counts, not resp_corrected data.
        """
        assert isinstance(stream, Stream)

        # Make sure that the array is not masked
        if np.any([np.ma.count_masked(tr.data) for tr in stream]):
            logging.warning('Unmasking masked data (usually a gap or overlap)')
            stream = stream.split().merge(fill_value='interpolate')
        out_stream = stream.copy()
        seed_ids = [tr.id for tr in out_stream]
        remove_seqs = {k: "" for k in seed_ids}

        if in_time_domain is True:
            logging.warning("Correcting traces in the time domain: VERY SLOW")
        else:
            logging.info("Correcting traces in the frequency domain")

        for dctf in self.DCTFs:
            tfs = dctf.tfs
            in_chan = tfs.input_channel
            for out_chan in tfs.output_channels:
                in_trace = out_stream.select(id=in_chan)[0]
                out_trace = out_stream.select(id=out_chan)[0]
                out_stream.remove(out_trace)
                out_trace = self._correct_trace(
                    in_trace,
                    out_trace,
                    tfs.freqs,
                    tfs.corrector_wrt_counts(out_chan),
                    in_time_domain,
                )
                out_trace.id = CS.insert(dctf.remove_sequence,
                                         CS.strip(out_trace.id))
                out_stream += out_trace
        return out_stream

    def plot(self):
        """
        Plot data cleaning transfer functions
        """
        self.DCTFs.plot()

    def _correct_trace(self, in_trace, out_trace, f, tf, in_time_domain=False):
        """
        Correct a trace using an input trace and a transfer function

        Note that the transfer function should be between the trace's units

        Args:
            in_trace (:class: `obspy.core.trace.Trace`): input trace
            out_trace (:class: `obspy.core.trace.Trace`): original output trace
            f (:class:`numpy.ndarray`): frequencies
            tf (:class:`numpy.ndarray`): transfer function between the input
                and output traces (counts/count)
            in_time_domain (bool): do correction in time domain

        Returns:
            out_trace_corrected (:class:`obspy.core.trace.Trace`): corrected
                output trace
        """
        self._validate_streams_synchronized(in_trace, out_trace)
        in_trace = in_trace.copy()
        out_trace = out_trace.copy()
        in_trace.detrend("linear")
        out_trace.detrend("linear")
        if in_time_domain:
            return self._correct_trace_time(in_trace, out_trace, f, tf)
        else:
            return self._correct_trace_freq(in_trace, out_trace, f, tf)

    @staticmethod
    def _correct_trace_time(in_trace, out_trace, f, tf):
        """
        Correct trace in the time domain

        Args:
            in_trace (:class:`obspy.core.trace.Trace`): input (noise) trace
            out_trace (:class:`obspy.core.trace.Trace`): trace to correct
            f (:class:`numpy.ndarray`): frequencies
            tf (:class:`numpy.ndarray`): out_trace/in_trace transfer function
        """
        if not 2 * f[-1] == in_trace.stats.sampling_rate:
            raise ValueError(
                "different tf ({}) & trace ({}) sample rates".format(
                    2 * f[-1], in_trace.stats.sampling_rate
                )
            )
        # Time domain transform of transfer function
        # Why don't we need to use the complex conjugate of tf?
        itf = np.fft.ifftshift(np.fft.irfft(tf))
        out_trace_corr = out_trace.copy()
        corr_data = np.convolve(in_trace.data.copy(), itf, mode="same")
        out_trace_corr.data -= signal.detrend(corr_data)
        return out_trace_corr

    @staticmethod
    def _correct_trace_freq(in_trace, out_trace, f, tf):
        """
        Correct trace in the frequency domain

        Args:
            in_trace (:class:`obspy.core.trace.Trace`): input (noise) trace
            out_trace (:class:`obspy.core.trace.Trace`): trace to correct
            f (:class:`numpy.ndarray`): frequencies
            tf (:class:`numpy.ndarray`): out_trace/in_trace transfer function
        """
        npts = in_trace.stats.npts
        trace_sr = in_trace.stats.sampling_rate
        fft_len = 2 ** int(np.ceil(npts)).bit_length()
        buff = np.zeros(fft_len - npts)
        # Transform all to frequency domain
        in_rfft = np.fft.rfft(np.concatenate((in_trace.data, buff)))
        out_rfft = np.fft.rfft(np.concatenate((out_trace.data, buff)))
        f_rfft = np.fft.rfftfreq(fft_len, 1.0 / trace_sr)
        tf_interp = np.interp(f_rfft, f, tf)
        # Subract coherent part and convert to time domain
        # Why isn't the complex conjugate of the tf used?
        new_out = np.fft.irfft(out_rfft - in_rfft * tf_interp)
        # new_out = np.fft.irfft(out_rfft - in_rfft * np.conj(tf_interp))
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

    def _removeDCTF_SDF(self, sdf, dctf, show_progress=False):
        """
        Directly clean a SpectralDensity object using transfer functions

        This is the ATACR method.  It works on the overall spectra, not
        individual windows.

        Used for the "fast_calc" method and for data already in
        SpectralDensity format

        Args:
            sdf (:class:`.SpectralDensity`): one-sided spectral density
                functions
            dctf (:class:`.DCTF`): Data Cleaner Transfer Functions w.r.t. one
                input channel
            show_progress (bool): show progress plots for each step

        Returns:
            sdf (:class:`.SpectralDensity`): cleaned spectral density functions
        """
        tf = dctf.tfs
        ic = tf.input_channel
        out_chans = tf.output_channels

        sdf = deepcopy(sdf)  # Avoids overwriting input object
        in_auto = sdf.autospect(ic)
        for oc in out_chans:
            # tf_oc = tf.values_counts(oc)
            tf_oc = tf.corrector(oc)
            # out_auto = sdf.autospect(strip_remove_one(oc))
            # crossspect = sdf.crossspect(ic, strip_remove_one(oc))
            out_auto = sdf.autospect(oc)
            crossspect = sdf.crossspect(ic, oc)
            if show_progress:
                self._plot_removeDCTF_SDF(
                    tf, sdf, in_auto, out_auto, tf_oc, ic, oc
                )
            # Bendat & Piersol 1986 eqs 6.35 and 6.41
            sdf.put_autospect(oc, out_auto - in_auto * np.abs(tf_oc) ** 2)
            # Bendat & Piersol eq6.36 with
            sdf.put_crossspect(ic, oc, crossspect - in_auto * tf_oc)
            sdf.replace_channel_name(oc, CS.insert(dctf.remove_sequence,
                                                   CS.strip(oc)))
        return sdf

    @staticmethod
    def _plot_removeDCTF_SDF(tf, sdf, in_auto, out_auto, tf_oc, ic, oc):
        """Plot elements of transfer function removal using SDFs"""
        f = tf.freqs
        fig, ax = plt.subplots(1, 1)
        ax2 = ax.twinx()
        ax2.semilogx(f, np.abs(sdf.coherence(ic, oc)),
                     "k--", alpha=0.5, label="coherence")
        ax2.axhline(np.abs(sdf.coh_signif(0.95)), color="k", ls=":",
                    alpha=0.5, label="coherence significance level")
        ax.loglog(f, np.abs(out_auto), label=oc)
        ax.loglog(f, np.abs(in_auto), alpha=0.5, label=ic)
        # ax.loglog(f, np.abs(tf_oc), label='tf^2')
        ax.loglog(f, np.abs(in_auto * np.abs(tf_oc) ** 2), alpha=0.5,
                  label=ic + "*tf^2")
        ax.loglog(f, np.abs(out_auto - in_auto * np.abs(tf_oc) ** 2),
                  alpha=0.5, label=f"{oc} - ({ic}*tf^2)")
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
        if len(sdfs) == len(self.DCTFs) + 1:
            ValueError("len(sdfs) isn't len(self.DCTFs)+1")
        remove_seqs = [""] + [t.remove_sequence for t in self.DCTFs]
        inputs = [""] + [t.remove_channel for t in self.DCTFs]
        n_channels = len(sdfs[0].channel_names)
        rows = int(np.floor(np.sqrt(n_channels)))
        cols = int(np.ceil(n_channels / rows))
        fig, axs = plt.subplots(rows, cols)
        irow, icol = 0, 0
        for chan in sdfs[0].channel_names:
            ax = axs[irow, icol]
            # make one line for each input
            inputs_base = [CS.strip(i) for i in inputs]
            ch_base = CS.strip(chan)
            if chan in inputs_base:
                i = inputs_base.index(chan)
                for sdf, rs in zip(sdfs[:i], remove_seqs[:i]):
                    ax.loglog(
                        sdf.freqs,
                        np.abs(sdf.autospect(ch_base + rs)),
                        alpha=0.75,
                        label=ch_base + rs,
                    )
            else:
                for sdf, rs in zip(sdfs, remove_seqs):
                    ax.loglog(
                        sdf.freqs,
                        np.abs(sdf.autospect(ch_base + rs)),
                        alpha=0.75,
                        label=ch_base + rs,
                    )
            if irow == rows - 1:
                ax.set_xlabel("Frequency (Hz)")
            ax.legend()
            icol += 1
            if icol >= cols:
                irow += 1
                icol = 0
        plt.show()


# def _list_match_pattern(pattern, chan_list):
#     """
#     Return list element matching pattern.  Pattern can include file wildcards
# 
#     Returns error if there are no matches, or more than one match
# 
#     Args:
#         pattern (str): pattern to match, can include "*" and "?" wildcards
#         chan_list (list of str): list of channels
#     """
#     selected = [x for x in chan_list if fnmatch.fnmatch(x, pattern)]
#     if len(selected) == 0:
#         raise ValueError(f'input pattern "{pattern}" not matched')
#     elif len(selected) > 1:
#         raise ValueError(
#             'more than one match for input pattern "{}": {}'.format(
#                 pattern, selected
#             )
#         )
#     return selected[0]


if __name__ == "__main__":
    import doctest

    doctest.testmod()
