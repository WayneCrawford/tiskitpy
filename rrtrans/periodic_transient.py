#!env python3
""" Model and remove glitches from BBOBS data"""

#################################################
# import obspy
from obspy.core import UTCDateTime
from obspy.core.stream import Stream
import diracComb_v2 as dc
import matplotlib.pyplot as plt
import numpy as np
# from scipy import signal
# import math as M
# import sys
# import glob

from .utils import stack_data, input_float, input_float_tuple

def_mag_limit = 5.85
def_days_per_magnitude = 1.5

class Periodic_Transient():
    """
    Information about one periodic transient
    """
    def __init__(self, name: str, period: float, dp: float, clips: tuple,
                 starttimes: list=[]):
        """
        :param name: name of this periodic transient (e.g., 'hourly')
        :param period: seconds between transients (seconds)
        :param dp: how many seconds to change the period by when testing for
                   better values
        :param clips: clip values outside of this range.  Should correspond to
                      max range of glitch
        :type clips: 2-tuple (lowest, highest)
        :param starttimes: onset times of earliest glitch(es).  Multiple values
                           allow for a change due to releveling, etc. It's
                           better to be a few seconds too early than too late.
                           The program will make transient slices starting
                           between a starttime and 1/3 of a transient period
                           earlier
        :type starttimes:  list of ~class `obspy.core.UTCDateTime`
        """
        self.name = name
        self.period = float(period)
        self.dp = float(dp)
        self.clips = [None, None]
        self.clips[0] = float(clips[0])
        self.clips[1] = float(clips[1])
        assert isinstance(starttimes, list)
        for t in starttimes:
            assert t is None or isinstance(t, UTCDateTime)
        self.starttimes = starttimes

    def __str__(self):
        s= '"{self.name}": {self.period}s+-{self.dp}, clips={self.clips}, '
        s += 'starttimes={self.starttimes}'
        return s

    def verify(self, trace, template, ref_time):
        """
        Interactively verify transient parameters
    
        :param trace:      input data trace
        :type trace: ~class `obspy.stream.Trace`
        :param template: template of good data (same size as trace,
                         1 for good, 0 for bad
        :type template: ~class `obspy.stream.Trace`
        :param ref_time: start time of first glitch (or a few seconds before)
        :type ref_time: ~class `obspy.core.UTCDateTime
        """
        startoffset = self._calc_start_offset(ref_time, trace.stats.starttime)

        # Set/verify clip levels
        self._ask_clips(trace, template, startoffset)

        # SET/VERIFY GLITCH PERIOD
        cliptrace = trace.copy()
        cliptrace.data.clip(self.clips[0], self.clips[1], out=cliptrace.data)
        self._ask_period(cliptrace, template, startoffset, ref_time)

    def _ask_clips(self, trace, template, st_off):
        """
        Show clip levels and ask to update them until acceptable

        :param trace: seismological trace
        :type trace: obspy.core.stream.Trace
        :param template: signal of ones and zeros, with zeros where trace is to
                         be hidden
        :type template: obspy.core.stream.Trace
        :param st_off: offset from start of trace to start slicing
        :param testper: length of each slice (seconds)
        :param clip: 2-tuple of default clip values (lo, hi)
        """
        st = trace.stats.starttime
        sps = trace.stats.sampling_rate
        station = trace.stats.station
        stack_trace = trace.copy()
        stack_trace.data = stack_trace.data * template.data
        if st_off > 0:
            stack_trace = stack_trace.slice(starttime=st + st_off)
        stack = stack_data(stack_trace.data, self.period*sps)
        nrows, ncols = stack.shape
        time = np.arange(nrows) / sps
        # slicenums = np.arange(ncols)
        title = '{} sliced at {:g}s, stacked'.format(station, self.period)
        # Show clip levels and verify that they are ok
        while True:
            plt.plot(time, stack, linewidth=0.1)
            plt.plot([time[0], time[-1]], [self.clips[0], self.clips[0]],
                     'k--', label='clip_lo = {:g}'.format(self.clips[0]))
            plt.plot([time[0], time[-1]], [self.clips[1], self.clips[1]],
                     'k-.', label='clip_hi = {:g}'.format(self.clips[1]))
            plt.xlabel('Time (seconds)')
            clip_range = self.clips[1] - self.clips[0]
            plt.ylim((self.clips[0] - clip_range*.3,
                      self.clips[1] + clip_range*.3))
            plt.legend()
            plt.title(title)
            plt.show()
            # Ask for a new clip_lo,clip_hi tuple, continue if current value
            # accepted
            newval = input_float_tuple('Enter clip levels containing all '
                                       'glitches (RETURN to accept current '
                                       'value)', self.clips)
            if newval == self.clips:
                break
            else:
                self.clips = newval

    def _ask_test_period(self, trace, template, st_off, ref_time):
        """
        Show glitch alignment and ask to update glitch period until acceptable

        Also allows to verify that the glitch_parm.starttime is ok, but not to
        modify it

        :param trace: seismological trace
        :type trace: obspy.core.stream.Trace
        :param template: signal of 1s and 0s, 0s where trace is to be hidden
        :type template: obspy.core.stream.Trace
        :param st_off: offset from start of trace to start slicing
        :param ref_time: time at start of (or just before) one of the glitches
        :type ref_time: ~class `obspy.core.UTCDateTime`
        """
        st = trace.stats.starttime
        sps = trace.stats.sampling_rate
        station = trace.stats.station
        stack_trace = trace.copy()
        stack_trace.data = stack_trace.data * template.data
        if st_off > 0:
            stack_trace = stack_trace.slice(starttime=st + st_off)
        while True:
            stack = stack_data(stack_trace.data, self.period*sps)
            nrows, ncols = stack.shape
            time = np.arange(nrows) / sps
            slicenums = np.arange(ncols)
            if ref_time:
                ref_offs = (ref_time-stack_trace.stats.starttime) % self.period
            else:
                ref_offs = 0
            title = '{} sliced at {:g}s, clipped, stacked'.format(station,
                                                                  self.period)

            # plot as pcolor
            plt.figure(1)
            plt.pcolormesh(slicenums, time, stack)
            plt.plot([0, ncols], [ref_offs, ref_offs], 'k--')
            plt.xlabel('Period slice #')
            plt.ylabel('Time (seconds)')
            plt.title(title)
            plt.grid(True, zorder=20)
            plt.show()

            # Ask for new test period, continue if current value accepted
            newval = input_float('Enter new test period (RETURN to accept '
                                 'current value)', self.period)
            if newval == self.period:
                break
            else:
                self.period = newval

    def _calc_start_offset(self, ref_time, starttime):
        """
        Choose "slice" start time so that glitches will start no more than 1/3
        of the way in
        """
        first_glitch_offset = (ref_time - starttime) % self.period
        if first_glitch_offset <= self.period / 3:
            offset = 0
        else:
            offset = first_glitch_offset - self.period/3.
            print('='*72)
            print('According to the first glitch time you specified, the '
                  'first glitch starts\n {:g} seconds ({:.0f}%) into slices '
                  'based on the {:g} second glitch period.\n Shifting forward by {:g} '
                  'seconds so that it starts 1/3 of the way in'
                  .format(first_glitch_offset,
                          100. * first_glitch_offset / self.period,
                          self.period, offset))
            print('='*72)
        return offset

    def calcglitch(self, trace, template, ref_time, iTrans, match, plots):
        """
        Calculate transient for a given trace and transient parameters

        :param trace:      input data trace
        :type trace: ~class `obspy.stream.Trace`
        :param template: template of good data (same size as trace, 1 for good,
                         0 for bad
        :type template: ~class `obspy.stream.Trace`
        :param ref_time: start time of first glitch (or a few seconds before)
        :type ref_time: ~class `obspy.core.UTCDateTime
        :param match: match and cancel each pulse separately
        :type match: bool
        :param iTrans: counter of which transient we are removing
        :param plots: plot results
        :type plots: bool
        """
        startoffset = self._calc_start_offset(ref_time, trace.stats.starttime)

        # Calculate and Remove glitch
        out, glit_1, glit_t, d_comb, nGlitch = dc.comb_clean(trace, self,
                                                             match, plots,
                                                             template,
                                                             startoffset)
        out.stats.channel = f'CC{iTrans:d}'
        glit_t.stats.channel = f'GG{iTrans:d}'
        if plots:
            Stream([trace, out, glit_t]).plot()
        trace = out
        return (out, glit_1, d_comb, nGlitch)

    def _get_transient_shifts(self, trace):
        """
        Search for glitch start corresponding to data start and
        make a list of glitch start shifts within data
        """
        shifts = []
        ref_time = self.starttime
        if len(self.starttime) > 1:
            for test_time in self.starttime[1:]:
                if test_time < trace.stats.starttime:
                    ref_time = test_time
                elif test_time < trace.stats.endtime:
                    # there is a break in glitch times within this data
                    shifts.append(test_time)
            return ref_time, shifts
        elif self.starttime:
            return self.starttime, shifts
        else:
            return trace.stats.starttime, shifts


#########################################################################
# if __name__ == "__main__":
# 	sys.exit(main())
