#!env python3
""" Model and remove transients from BBOBS data"""
from scipy.signal import butter, sosfiltfilt
from obspy.core import UTCDateTime
from obspy.core import Stream  # , Trace
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from ..time_spans import TimeSpans

from .dirac_comb import comb_calc, comb_remove
from .utils import stack_data, input_float, input_float_tuple
from ..logger import init_logger

logger = init_logger()
def_mag_limit = 5.85
def_days_per_magnitude = 1.5


class PeriodicTransient:
    """
    Class to determine parameters for and remove a periodic transient

    The program will make transient slices starting between
    transient_starttime and 1/3 of self.period earlier

    To include the entire transient in each trace slice,
    transient_starttime must not be too late. A few seconds early is OK.

    Parameters:
        name (str): name of this periodic transient (e.g., 'hourly')
        period (float): seconds between each transient
        dp (float): how many seconds to change the period by when testing
                   for better values.  Should be equal to or slightly greater
                   than the uncertainty in your estimate of `period`
        clips (tuple): clip values outside of this range (low, high).
                Should contain the max range of the transient
        transient_starttime (str or ~class `obspy.core.UTCDateTime`): onset
                time of earliest transient.
        freq_HP (float or bool): highpass corner frequency used for training
            and matching.  If True, sets to 1/period.  If False, do not cut off
            low frequencies
        freq_LP (float): lowpass corner frequency used for training and
            matching (0.05 is a good value to remove microseisms)
        transient_model: Model of the transient
        dirac_comb: the dirac comb to convolve with the transient
        n_transients_used: number of transients used to calculate the model
        self.tm: transient model shifted one left
        self.tp: transient model shifted one right
    """

    def __init__(self, name, period, dp, clips, transient_starttime,
                 freq_HP=True, freq_LP=0.05):
        """
        Constructor

        Args:
            name (str): name of this periodic transient (e.g., 'hourly')
            period (float): seconds between each transient
            dp (float): how many seconds to change the period by when testing
                       for better values.  Should be equal to or slightly greater
                       than the uncertainty in your estimate of `period`
            clips (tuple): clip values outside of this range (low, high).
                    Should contain the max range of the transient
            transient_starttime (str or ~class `obspy.core.UTCDateTime`): onset
                    time of earliest transient.
            freq_HP (float or bool): highpass corner frequency used for training
                and matching.  If True, sets to 1/period.  If False, do not cut off
                low frequencies
            freq_LP (float): lowpass corner frequency used for training and
                matching (0.05 is a good value to remove microseisms)
        """
        self.name = name
        self.period = float(period)
        self.dp = float(dp)
        self.clips = (float(clips[0]), float(clips[1]))
        if isinstance(transient_starttime, str):
            transient_starttime = UTCDateTime(transient_starttime)
        assert isinstance(transient_starttime, UTCDateTime)
        self.transient_starttime = transient_starttime
        self.freq_HP = freq_HP
        self.freq_LP = freq_LP
        # Values to be calculated
        self.transient_model = None
        self.dirac_comb = None
        self.n_transients_used = 0
        # These should probably just be dependent parameters based on
        # self.transient_model
        self.tm = None
        self.tp = None

    def __str__(self):
        s = f'"{self.name}": {self.period:.2f}s+-{self.dp}, clips={self.clips}'
        s += f", training and matching freq bounds=({self.freq_HP}, {self.freq_LP})"
        s += f", transient_starttime={self.transient_starttime}"
        return s

    def calc_timing(self, trace, eq_spans=True):
        """
        Interactively calculate transient timing parameters.

        The user chooses the transient period, clip levels and starttime
        best fitting the data.
        
        Args:
            trace (class obspy.core.Trace): data
            eq_spans (:class:`TimeSpans` or bool): 
                if TimeSpans: times spans in which to zero data
                if True: use TimeSpans.from_eqs()
        """
        if eq_spans is True:
            eq_spans = TimeSpans.from_eqs(trace)
        elif eq_spans is False:
            eq_spans = TimeSpans(None)

        slice_starttime = self._calc_slice_starttime(trace)

        # Set/verify clip levels
        cliptrace = self._filtered_trace(trace)
        self._ask_clips(cliptrace, eq_spans, slice_starttime)
        print(f'clips level set to {self.clips}')

        # Set/verify transient period
        cliptrace.data.clip(self.clips[0], self.clips[1], out=cliptrace.data)
        self._ask_period_or_time(cliptrace, eq_spans, slice_starttime, True)
        print(f'period set to {self.period}')
        
        # Set/verify transient_starttime
        self._ask_period_or_time(cliptrace, eq_spans, slice_starttime, False)
        print(f'transient starttime set to {self.transient_starttime}')

    def calc_transient(self, trace, eq_spans=True, match=True, plots=False):
        """
        Calculate transient for a given trace and transient parameters

        Args:
            trace (:class:`obspy.core.stream.Trace`): input data trace
            eq_spans (:class:`TimeSpans` or bool): 
                if TimeSpans: times spans in which to zero data
                if True: use TimeSpans.from_eqs()
            plots (bool): plot results
        """
        if eq_spans is True:
            eq_spans = TimeSpans.from_eqs(trace)
        elif eq_spans is False:
            eq_spans = TimeSpans(None)

        slice_starttime = self._calc_slice_starttime(trace)
        if self.transient_starttime < trace.stats.starttime:
            print("\tshifting transient startime to first within data")
            self.transient_starttime = (
                trace.stats.starttime + self._transient_offset(trace)
            )
        filttrace = self._filtered_trace(trace)
        transient, dc, nG, tm, tp, cbuff = comb_calc(
            filttrace, self, plots, eq_spans, slice_starttime
        )
        transient.stats.channel = f"TR{trace.stats.channel[-1]}"
        if plots:
            transient.plot()
        self.transient_model = transient
        self.tm, self.tp = tm, tp
        self.dirac_comb = dc
        self.n_transients_used = nG
        # print(f"calc_transient(): {cbuff=}")
        self.comb_buffer = cbuff

    def _filtered_trace(self, tr, zerophase=True):
        filt_trace = tr.copy()
        if self.freq_HP is False:
            filt_trace.filter("lowpass", freq=self.freq_LP, zerophase=zerophase)
        elif self.freq_HP is True:
            filt_trace.filter("bandpass", freqmin=1./self.period, freqmax=self.freq_LP, zerophase=zerophase)
        else:
            filt_trace.filter("bandpass", freqmin=self.freq_HP, freqmax=self.freq_LP, zerophase=zerophase)
        return filt_trace

    def remove_transient(self, trace, match=True, plots=False):
        """
        Remove transient from trace

        Args:
            trace(:class:`obspy.Trace`): input data
            match (bool): match each transient individually.  This is useful
                if the transient dominates the trace, not if the transient is
                subtle.
        Returns:
            (:class:`obspy.Trace`): output data
        """
        assert self.transient_model is not None

        slice_starttime = self._calc_slice_starttime(trace)
        out, synth = comb_remove(trace, self._filtered_trace(trace),
                                 self, match, slice_starttime)
        if plots:
            # Change channel names for the plot
            out.stats.channel = "CLN"
            synth.stats.channel = "SYN"
            Stream([trace, out, synth]).plot(method="full")
            # Revert output channel name
            out.stats.channel=trace.stats.channel
        return out

    def _ask_clips(self, trace, eq_spans, slice_starttime):
        """
        Show clip levels and ask to update them until acceptable

        Args:
            trace (~class `obspy.core.stream.Trace``): seismological trace
            slice_starttime (UTCDateTime): first slice starttime
            eq_spans (:class:`TimeSpans` or bool): times spans in which to
                zero data
        """
        stt = trace.stats.starttime
        sps = trace.stats.sampling_rate
        seed_id = trace.id

        stack_trace = trace.copy()
        stack_trace = eq_spans.zero(stack_trace)
        if slice_starttime > stt:
            stack_trace = stack_trace.slice(starttime=slice_starttime)
        stack = stack_data(stack_trace.data, self.period * sps)
        nrows, ncols = stack.shape
        time = np.arange(nrows) / sps
        # slicenums = np.arange(ncols)
        title = f"{seed_id}, sliced at {self.period:g}s, stacked, clips={self.clips}"

        # Show clip levels and verify that they are ok
        fig, ax = plt.subplots(1, 1, num="Select clip levels")
        ax.plot(time, stack, linewidth=0.1)
        c1, c2 = self.clips
        llo = ax.axhline(c1, c="b", ls="--", label="clip_lo")
        lhi = ax.axhline(c2, c="r", ls="--", label="clip_hi")
        ax.set_xlabel("Time (seconds)")
        ax.set_title(title)
        ax.set_ylim((c1 - (c2 - c1) * 0.3, c2 + (c2 - c1) * 0.3))
        ax.legend(loc="upper left")
        plt.ion()
        plt.show()
        while True:
            newval = input_float_tuple(
                "Enter clip levels containing all transients", self.clips
            )
            if newval == self.clips:
                break
            else:
                self.clips = newval
                c1, c2 = self.clips
                llo.set_ydata([c1, c1])
                lhi.set_ydata([c2, c2])
                ax.set_ylim((c1 - (c2 - c1) * 0.3, c2 + (c2 - c1) * 0.3))
                plt.draw()
        plt.close(fig)
        plt.ioff()

    def _ask_period_or_time(self, trace, eq_spans, slice_starttime, ask_period=True):
        """
        Show transient alignment and ask to update period until acceptable

        Also allows to verify that self.transient_starttime is ok, but not to
        modify it

        Args:
            trace (~class obspy.core.stream.Trace): seismological trace
            eq_spans (:class:`TimeSpans`): times spans in which to zero data
            slice_starttime (UTCDateTime): first slice starttime
            ask_period (bool): True: ask for period, False: ask for time
        """
        stt = trace.stats.starttime
        sps = trace.stats.sampling_rate
        seed_id = trace.id

        stack_trace = trace.copy()
        stack_trace = eq_spans.zero(stack_trace)
        if slice_starttime > stt:
            stack_trace = stack_trace.slice(starttime=slice_starttime)
        fig, ax = plt.subplots(1, 1, num="Select transient period")
        ax.set_xlabel("Slice starttime")
        ax.set_ylabel("Time (seconds)")
        locator = mdates.AutoDateLocator()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
        ax.grid(True, zorder=20)
        hline = ax.axhline(0, c="k", ls="--")
        fig.autofmt_xdate()
        plt.ion()
        plt.show()
        
        if ask_period is True:
            print('\nIf the slope is positive, INCREASE the period')
        else:
            print('\nIf the dotted line is beneath the transient, enter a POSTIVE offset')
        while True:
            stack = stack_data(stack_trace.data, self.period * sps)
            nrows, ncols = stack.shape
            timey = np.arange(nrows) / sps
            x_offset = np.arange(ncols) * self.period
            timex = (
                stack_trace.stats.starttime.matplotlib_date + x_offset / 86400
            )
            ref_offs = self._transient_offset(stack_trace)
            if ask_period is True:
                title = f"{seed_id}, period={self.period:g}s"
            else:
                title = f'{seed_id}, transient_starttime={self.transient_starttime.strftime("%Y%m%dT%H%M%S")}'
            # plot as pcolor
            ax.set_title(title, size="medium")
            # ax.pcolormesh(slicenums, time, stack, shading='auto')
            ax.pcolormesh(timex, timey, stack, shading="auto")
            hline.set_ydata([ref_offs, ref_offs])
            plt.draw()

            if ask_period is True:
                # Ask for new test period, continue if current value accepted
                newval = input_float("Enter new test period", self.period)
                if newval == self.period:
                    break
                else:
                    self.period = newval
            else:
                # Ask for transient_starttime offset
                print(f'transient_starttime = {self.transient_starttime}')
                offset = input_float("Enter offset seconds", 0.)
                if offset == 0:
                    break
                else:
                    self.transient_starttime += offset
        plt.close(fig)
        plt.ioff()

    def _calc_slice_starttime(self, trace, verbose=True):
        """
        Choose first "slice" starttime so that transients start no more
        than 1/3 of the way in

        :param trace: input data
        """
        slice_starttime = trace.stats.starttime
        transient_offset = self._transient_offset(trace)
        max_offset = self.period / 3
        if transient_offset > max_offset:
            if verbose:
                print(f"\told {slice_starttime=}")
            slice_starttime += transient_offset - max_offset
            if verbose:
                print(f"\t{self.transient_starttime=}")
                print(f"\t{self.period=}")
                print(
                    "\tFirst transient starts {:.0f}% of period into data".format(
                        100.0 * transient_offset / self.period
                    )
                )
                print(f"\tReducing to 33% by shifting forward {transient_offset - max_offset:g}s")
                print(f"\tnew {slice_starttime=}")
        return slice_starttime

    def _transient_offset(self, trace):
        """
        Return offset of first transient in the trace
        """
        return (self.transient_starttime - trace.stats.starttime) % self.period


#########################################################################
# if __name__ == "__main__":
# 	sys.exit(main())
