#!env python3
"""Class of time spans to remove, keep, zero, etc. in Trace or Stream data"""

import numpy as np

# import obspy
from obspy.core import UTCDateTime
from obspy.core.stream import Stream, Trace
from dataclasses import dataclass


@dataclass
class TimeSpans:
    """
    A class specifying time spans, to be removed, kept, zeroed, etc.
    """
    start_times: list
    end_times: list

    def __post_init__(self):
        if not self.validate():
            self.organize()

    def __len__(self):
        return len(self.start_times)

    def __eq__(self, other):
        if not len(self.start_times) == len(other.start_times):
            return False
        for sta, eta, stb, etb in zip(self.start_times, self.end_times,
                                      other.start_times, other.end_times):
            if not sta == stb:
                return False
            if not eta == etb:
                return False
        return True

    def validate(self):
        """
        Returns false if start_times are not ordered or there are overlaps
        """
        if not len(self.start_times) == len(self.end_times):
            raise ValueError("starttimes[] and end_times[] are not the "
                             "same length")
        for x in self.start_times:
            if not isinstance(x, UTCDateTime):
                return False
        for x in self.end_times:
            if not isinstance(x, UTCDateTime):
                return False
        if not sorted(self.start_times) == self.start_times:
            return False
        for startafter, endbefore in zip(self.start_times[1:],
                                         self.end_times[:-2]):
            if startafter < endbefore:
                return False

    def organize(self):
        """
        Order starttimes and endtimes by increasing starttime, and consolidate
        overlapping time spans
        """
        # sort by time
        self.end_times = [x for _, x in sorted(zip(self.start_times,
                                                   self.end_times))]
        self.start_times = sorted(self.start_times)

        # remove any overlaps
        start_times = [self.start_times[0]]
        end_times = [self.end_times[0]]
        for st, et in zip(self.start_times[1:], self.end_times[1:]):
            if st > end_times[-1]:
                start_times.append(st)
                end_times.append(et)
            # otherwise extend the previous entry
            else:
                end_times[-1] = et
        self.start_times = start_times
        self.end_times = end_times

    def zero(self, inp, plot=False):
        """
        Zero out data in the time spans

        Arguments:
            inp (Trace or Stream): seismological data
            plot: plot traces with spans cut out

        Returns:
            trace with spans set to zero
        """
        if isinstance(inp, Trace):
            stream = Stream([inp])
        elif isinstance(inp, Stream):
            stream = inp.copy()
        else:
            raise(ValueError, 'inp is not an obspy Trace or Stream')
        for tr in stream:
            tr.stats.channel = 'XX' + tr.stats.channel[2]

            for st, et in zip(self.start_times, self.end_times):
                start_addr, end_addr = self._get_addrs(st, et, tr)
                if start_addr == None:
                    continue
                tr.data[start_addr:end_addr+1] = 0.
        if plot:
            (stream + inp).plot(color='blue', equal_scale=False)
        if isinstance(inp, Trace):
            return stream[0]
        return stream

    def interp(self, inp, plot=False):
        """
        Interpolate data from the start to end values in each time span

        Arguments:
            inp (Trace or Stream): seismological data
            plot: plot traces with spans cut out

        Returns:
            trace with spans set to zero
        """
        if isinstance(inp, Trace):
            stream = Stream([inp])
        elif isinstance(inp, Stream):
            stream = inp.copy()
        else:
            raise(ValueError, 'inp is not an obspy Trace or Stream')
        for tr in stream:
            tr.stats.channel = 'XX' + tr.stats.channel[2]

            for st, et in zip(self.start_times, self.end_times):
                start_addr, end_addr = self._get_addrs(st, et, tr)
                if start_addr == None:
                    continue
                tr.data[start_addr:end_addr+1] = np.linspace(
                    tr.data[start_addr], tr.data[end_addr],
                    end_addr - start_addr + 1)
        if plot:
            (stream + inp).plot(color='blue', equal_scale=False)
        if isinstance(inp, Trace):
            return stream[0]
        return stream

    def _get_addrs(self, start_time, end_time, trace):
        if end_time < trace.stats.starttime or start_time > trace.stats.endtime:
            return None, None
        start_addr = max(np.floor((start_time - trace.stats.starttime)
                         * trace.stats.sampling_rate), 0)
        end_addr = min(np.ceil((end_time - trace.stats.starttime)
                       * trace.stats.sampling_rate), trace.stats.npts)
        # print(st, et, start_sample, end_sample)
        return int(start_addr), int(end_addr)
        

    def cutout(self, inp, plot=False):
        """
        Cut out data in the time spans (using Trace/Stream.cutout)

        Arguments:
            inp (Trace or Stream): seismological data

        Returns:
            new Trace or Stream
        """
        if not isinstance(inp, Trace) and not isinstance(inp, Stream):
            raise(ValueError, 'inp is not an obspy Trace or Stream')
        outp = inp.copy()
        for tr in outp:
            tr.stats.channel = 'XX' + tr.stats.channel[2]

            for st, et in zip(self.start_times, self.end_times):
                # Skip events that don't cover the trace time range
                if et < tr.stats.starttime:
                    continue
                if st > tr.stats.endtime:
                    continue
                start_sample = max(st - tr.stats.starttime, 0)
                end_sample = min(et - tr.stats.starttime, tr.stats.npts)
                tr.data[start_sample:end_sample] = 0.
        if plot:
            (outp + inp).plot(color='blue', equal_scale=False)
        if isinstance(inp, Trace):
            return outp[0]
        return outp