#!env python3
""" Model and remove transients from BBOBS data"""

# import math as M
# import sys
# import glob
# import obspy

from obspy.core import UTCDateTime
from obspy.core import Stream, Trace
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
# from scipy import signal

from .dirac_comb import comb_calc, comb_remove
from .utils import stack_data, input_float, input_float_tuple

def_mag_limit = 5.85
def_days_per_magnitude = 1.5


class PeriodicTransient():
    """
    Information about one periodic transient
    """
    def __init__(self, name: str, period: float, dp: float, clips: tuple,
                 transient_starttime: UTCDateTime):
        """
        :param name: name of this periodic transient (e.g., 'hourly')
        :param period: seconds between transients (seconds)
        :param dp: how many seconds to change the period by when testing for
                   better values
        :param clips: clip values outside of this range.  Should correspond to
                      max range of transient
        :type clips: 2-tuple (lowest, highest)
        :param transient_starttime: onset time of earliest transient.
        :type transient_starttime:  ~class `obspy.core.UTCDateTime`
        
        The program will make transient slices starting between
        transient_starttime and 1/3 of self.period earlier
        
        transient_starttime must not be too late, a few seconds early is ok
        """
        self.name = name
        self.period = float(period)
        self.dp = float(dp)
        self.clips = [None, None]
        self.clips[0] = float(clips[0])
        self.clips[1] = float(clips[1])
        if isinstance(transient_starttime, str): 
                transient_starttime = UTCDateTime(transient_starttime)
        assert isinstance(transient_starttime, UTCDateTime)
        self.transient_starttime = transient_starttime
        # Values to be calculated
        self.transient_model = None
        self.dirac_comb = None
        self.n_transients_used = 0
        # These should probably just be dependent parameters based on
        # self.transient_model
        self.tm = None
        self.tp = None

    def __str__(self):
        s = f'"{self.name}": {self.period:.2f}s+-{self.dp}, clips={self.clips}, '
        s += f'transient_starttime={self.transient_starttime}'
        return s

    def calc_timing(self, trace, eq_template):
        """
        Calculate transient time parameters

        :param trace: data
        :param eq_template: EQTemplate() object
        """
        self.verify(trace, eq_template)

    def calc_transient(self, trace, eq_template, match=True, plots=False):
        """
        Calculate transient for a given trace and transient parameters

        :param trace:      input data trace
        :type trace: ~class `obspy.stream.Trace`
        :param eq_template: template of good data (same size as trace,
                            1 for good, 0 for bad
        :type eq_template: ~class `EQTemplate` or `obspy.core.Trace`
        :param match: match and cancel each pulse separately
        :type match: bool
        :param plots: plot results
        :type plots: bool
        :returns:  cleaned data
        """
        slice_starttime = self._calc_slice_starttime(trace)
        if self.transient_starttime < trace.stats.starttime:
            print("\tshifting transient startime to first within data")
            self.transient_starttime = trace.stats.starttime + self._transient_offset(trace)
        transient, dc, nG, tm, tp, cbuff = comb_calc(trace, self, plots,
                                                     eq_template,
                                                     slice_starttime)
        transient.stats.channel = f'TR{trace.stats.channel[-1]}'
        if plots:
            transient.plot()
        self.transient_model = transient
        self.tm, self.tp = tm, tp
        self.dirac_comb = dc
        self.n_transients_used = nG
        # print(f"calc_transient(): {cbuff=}")
        self.comb_buffer = cbuff

    def remove_transient(self, trace, match=True, plots=False):
        """
        Remove transient from trace

        :param trace: data
        :param match: match each transient individually
        """
        assert self.transient_model is not None

        slice_starttime = self._calc_slice_starttime(trace)
        out, synth = comb_remove(trace, self, match, slice_starttime)
        if plots:
            out.stats.channel = 'CLN'
            synth.stats.channel = 'SYN'
            Stream([trace, out, synth]).plot()
        return out

    def verify(self, trace, template):
        """
        Interactively verify transient parameters

        :param trace:      input data trace
        :type trace: ~class `obspy.stream.Trace`
        :param template: template of good data (same size as trace,
                         1 for good, 0 for bad
        :type template: ~class `obspy.stream.Trace`
        """
        slice_starttime = self._calc_slice_starttime(trace)

        # Set/verify clip levels
        self._ask_clips(trace, template, slice_starttime)

        # SET/VERIFY TRANSIENT PERIOD
        cliptrace = trace.copy()
        cliptrace.data.clip(self.clips[0], self.clips[1], out=cliptrace.data)
        self._ask_period(cliptrace, template, slice_starttime)

    def _ask_clips(self, trace, template, slice_starttime):
        """
        Show clip levels and ask to update them until acceptable

        :param trace: seismological trace
        :type trace: obspy.core.stream.Trace
        :param template: signal of ones and zeros, with zeros where trace is to
                         be hidden
        :type template: obspy.core.stream.Trace
        :param slice_starttime: first slice starttime
        :param testper: length of each slice (seconds)
        :param clip: 2-tuple of default clip values (lo, hi)
        """
        stt = trace.stats.starttime
        sps = trace.stats.sampling_rate
        sta = trace.stats.station
        
        stack_trace = trace.copy()
        stack_trace.data = stack_trace.data * template.data
        if slice_starttime > stt:
            stack_trace = stack_trace.slice(starttime=slice_starttime)
        stack = stack_data(stack_trace.data, self.period*sps)
        nrows, ncols = stack.shape
        time = np.arange(nrows) / sps
        # slicenums = np.arange(ncols)
        title = '{} sliced at {:g}s, stacked'.format(sta, self.period)
        
        # Show clip levels and verify that they are ok
        fig, ax = plt.subplots(1,1, num="Select clip levels")
        ax.plot(time, stack, linewidth=0.1)
        c1, c2 = self.clips
        llo = ax.axhline(c1, c='b', ls='--', label=f'clip_lo')
        lhi = ax.axhline(c2, c='r', ls='--', label=f'clip_hi')
        ax.set_xlabel('Time (seconds)')
        ax.set_title(title)
        ax.set_ylim((c1 - (c2-c1)*.3, c2 + (c2-c1)*.3))
        ax.legend()
        plt.ion()
        plt.show()
        while True:
            newval = input_float_tuple(
                'Enter clip levels containing all transients', self.clips)
            if newval == self.clips:
                break
            else:
                self.clips = newval
                c1, c2 = self.clips
                llo.set_ydata([c1, c1])
                lhi.set_ydata([c2, c2])
                ax.set_ylim((c1 - (c2-c1)*.3, c2 + (c2-c1)*.3))
                plt.draw()
        plt.close(fig)
        plt.ioff()

    def _ask_period(self, trace, template, slice_starttime):
        """
        Show transient alignment and ask to update period until acceptable

        Also allows to verify that self.transient_starttime is ok, but not to
        modify it

        :param trace: seismological trace
        :type trace: obspy.core.stream.Trace
        :param template: signal of 1s and 0s, 0s where trace is to be hidden
        :type template: obspy.core.stream.Trace
        :param slice_starttime: offset from start of trace to start slicing
        """
        stt = trace.stats.starttime
        sps = trace.stats.sampling_rate
        sta = trace.stats.station
        
        stack_trace = trace.copy()
        stack_trace.data = stack_trace.data * template.data
        if slice_starttime > stt:
            stack_trace = stack_trace.slice(starttime=slice_starttime)
        fig, ax = plt.subplots(1,1, num="Select transient period")
        ax.set_xlabel('Slice starttime')
        ax.set_ylabel('Time (seconds)')
        locator = mdates.AutoDateLocator()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
        ax.grid(True, zorder=20)
        hline = ax.axhline(0, c='k', ls='--')
        fig.autofmt_xdate()
        plt.ion()
        plt.show()
        while True:
            stack = stack_data(stack_trace.data, self.period*sps)
            nrows, ncols = stack.shape
            timey = np.arange(nrows) / sps
            x_offset  = np.arange(ncols)*self.period
            timex = stack_trace.stats.starttime.matplotlib_date + x_offset/86400
            slicenums = np.arange(ncols)
            ref_offs = self._transient_offset(stack_trace)
            title = '{}, transient, data start={},{}, period={:g}s'.format(
                sta, self.transient_starttime.strftime('%Y%m%dT%H%M%S'),
                stt.strftime('%Y%m%dT%H%M%S'), self.period)
            # plot as pcolor
            ax.set_title(title, size='medium')
            # ax.pcolormesh(slicenums, time, stack, shading='auto')
            ax.pcolormesh(timex, timey, stack, shading='auto')
            hline.set_ydata([ref_offs, ref_offs])
            plt.draw()

            # Ask for new test period, continue if current value accepted
            newval = input_float('Enter new test period', self.period)
            if newval == self.period:
                break
            else:
                self.period = newval
        plt.close(fig)
        plt.ioff()

    def _calc_slice_starttime(self, trace):
        """
        Choose first "slice" starttime so that transients start no more
        than 1/3 of the way in
        
        :param trace: input data
        """
        slice_starttime = trace.stats.starttime
        transient_offset = self._transient_offset(trace)
        if transient_offset > self.period / 3:
            shift = transient_offset - self.period/3.
            slice_starttime += shift
            print(f'\t{slice_starttime=}')
            print(f'\t{self.transient_starttime=}')
            print(f'\t{self.period=}')
            print('\tFirst transient starts {:.0f}% of period into data'
                  .format(100. * transient_offset / self.period))
            print(f'\tReducing to 33% by shifting forward {shift:g}s')
        return slice_starttime

    def _transient_offset(self, trace):
        """
        Return offset of first transient in the trace
        """
        return (self.transient_starttime - trace.stats.starttime) % self.period

#########################################################################
# if __name__ == "__main__":
# 	sys.exit(main())
