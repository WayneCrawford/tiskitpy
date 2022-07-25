#!env python3
"""Class of time spans to remove, keep, zero, etc. in Trace or Stream data"""
from pathlib import Path
from dataclasses import dataclass

from obspy.clients.fdsn import Client
from obspy.core.event import Catalog, read_events
from obspy.core import UTCDateTime
from obspy.core.stream import Stream, Trace
import numpy as np
from matplotlib import pyplot as plt


@dataclass
class TimeSpans:
    """
    A class specifying time spans, to be removed, kept, zeroed, etc.
    """

    start_times: list
    end_times: list

    def __post_init__(self):
        self.organize()

    def __len__(self):
        return len(self.start_times)

    def __eq__(self, other):
        if not len(self.start_times) == len(other.start_times):
            return False
        for sta, eta, stb, etb in zip(
            self.start_times,
            self.end_times,
            other.start_times,
            other.end_times,
        ):
            if not sta == stb:
                return False
            if not eta == etb:
                return False
        return True

    def __str__(self):
        s = "TimeSpans: start            |            end\n"
        s += "===========================+===============================\n"
        for st, et in zip(self.start_times, self.end_times):
            s += f" {st} | {et}\n"
        return s

    @classmethod
    def from_eqs(
        cls,
        starttime,
        endtime,
        minmag=5.85,
        days_per_magnitude=1.5,
        eq_file=None,
        save_eq_file=True,
        quiet=False,
    ):
        """
        Generate timespans to avoid because of earthquakes

        Will read earthquakes from the USGS online catalog the first time,
        saving the information to a file that can be subsequently used

        Args:
            starttime (UTCDateTime or str): earliest data that will be
                presented.  If a str, must by ISO8601 compatible
            endtime (UTCDateTime or str): latest data that will be presented.
                If a str, must by ISO8601 compatible
            minmag (float): EQ Magnitude above which to cut out times
            days_per_magnitude (float): days to cut per magnitude above
                min_magnitude
            eq_file (str): the eq filename (otherwise, generates it)
            save_eq_file (bool): save the catalog file for future use
        Returns:
            eq_spans (~class `.TimeSpans`): time spans covering EQ signal
        """
        if isinstance(starttime, str):
            try:
                starttime = UTCDateTime(starttime)
            except Exception:
                raise ValueError(f"UTCDateTime() could not read {starttime=}")
        if isinstance(endtime, str):
            try:
                endtime = UTCDateTime(endtime)
            except Exception:
                raise ValueError(f"UTCDateTime() could not read {endtime=}")
        if eq_file is None:
            eq_file = _eq_filename(starttime, endtime, minmag)
        if Path(eq_file).is_file():
            cat = read_events(eq_file, format="quakeml")
        else:
            if not quiet:
                print("Reading EQs from USGS online catalog...",
                      end="", flush=True)
            cat = Client("USGS").get_events(
                starttime=starttime
                - _calc_eq_cut(9, minmag, days_per_magnitude),
                endtime=endtime,
                minmagnitude=minmag,
                orderby="time-asc",
            )
            if not quiet:
                print("Done", flush=True)
                print(f'writing catalog to "{eq_file}"')
            if save_eq_file:
                cat.write(eq_file, format="quakeml")

        new_cat = Catalog(
            events=[x for x in cat if x.preferred_magnitude().mag >= minmag]
        )
        start_times = [x.preferred_origin().time for x in new_cat]
        end_times = [
            x.preferred_origin().time
            + _calc_eq_cut(
                x.preferred_magnitude().mag, minmag, days_per_magnitude
            )
            for x in new_cat
        ]
        return cls(start_times, end_times)

    def validate(self):
        """
        Returns false if there is a problem

        Problems:
            - start_times are not ordered
            - there are overlaps
            - len(start_times) ~= len(end_times)
            - not every object in start_times and end_times is a UTCDateTime
        """
        if not len(self.start_times) == len(self.end_times):
            raise ValueError(
                "starttimes[] and end_times[] are not the " "same length"
            )
        for x in self.start_times:
            if not isinstance(x, UTCDateTime):
                raise ValueError("There is a non-UTCDateTime starttime value")
        for x in self.end_times:
            if not isinstance(x, UTCDateTime):
                raise ValueError("There is a non-UTCDateTime endtime value")
        if not sorted(self.start_times) == self.start_times:
            return False
        for startafter, endbefore in zip(
            self.start_times[1:], self.end_times[:-2]
        ):
            if startafter < endbefore:
                return False

    def organize(self):
        """
        Order starttimes and endtimes by increasing starttime, and consolidate
        overlapping time spans
        """
        if self.validate() is True:
            return

        # sort by time
        self.end_times = [
            x for _, x in sorted(zip(self.start_times, self.end_times))
        ]
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

    def append(self, new_time_spans):
        """
        Appends TimeSpan object to self

        Args:
            new_time_spans (~class `TimeSpans`): time spans to append
        """
        new_time_spans.validate()
        self.validate()
        self.start_times.extend(new_time_spans.start_times)
        self.end_times.extend(new_time_spans.end_times)
        self.organize()

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
            raise (ValueError, "inp is not an obspy Trace or Stream")
        for tr in stream:
            tr.stats.channel = "XX" + tr.stats.channel[2]

            for st, et in zip(self.start_times, self.end_times):
                start_addr, end_addr = self._get_addrs(st, et, tr.stats)
                if start_addr is None:
                    continue
                tr.data[start_addr : end_addr + 1] = 0.0
        if plot:
            (stream + inp).plot(color="blue", equal_scale=False)
        if isinstance(inp, Trace):
            return stream[0]
        return stream

    def has_zeros(self, starttime, endtime):
        """
        Does a trace's time span intersect any of the TimeSpans?

        Arguments:
            starttime (UTCDateTime): start time
            endtime (UTCDateTime): end time

        Returns:
            (bool):
        """
        for st, et in zip(self.start_times, self.end_times):
            if st < starttime and et > starttime:
                return True
            elif st >= starttime and st < endtime:
                return True
        return False

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
            raise (ValueError, "inp is not an obspy Trace or Stream")
        for tr in stream:
            tr.stats.channel = "XX" + tr.stats.channel[2]

            for st, et in zip(self.start_times, self.end_times):
                start_addr, end_addr = self._get_addrs(st, et, tr.stats)
                if start_addr is None:
                    continue
                tr.data[start_addr : end_addr + 1] = np.linspace(
                    tr.data[start_addr],
                    tr.data[end_addr],
                    end_addr - start_addr + 1,
                )
        if plot:
            (stream + inp).plot(color="blue", equal_scale=False)
        if isinstance(inp, Trace):
            return stream[0]
        return stream

    def _get_addrs(self, starttime, endtime, stats):
        if not isinstance(starttime, UTCDateTime):
            raise TypeError(f"starttime is {type(starttime)}, not UTCDateTime")
        if not isinstance(endtime, UTCDateTime):
            raise TypeError(f"endtime is {type(endtime)}, not UTCDateTime")
        if endtime < stats.starttime or starttime > stats.endtime:
            return None, None
        start_addr = max(
            np.floor((starttime - stats.starttime) * stats.sampling_rate), 0
        )
        end_addr = min(
            np.ceil((endtime - stats.starttime) * stats.sampling_rate),
            stats.npts,
        )
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
            raise (ValueError, "inp is not an obspy Trace or Stream")
        outp = inp.copy()
        for tr in outp:
            tr.stats.channel = "XX" + tr.stats.channel[2]

            for st, et in zip(self.start_times, self.end_times):
                # Skip events that don't cover the trace time range
                if et < tr.stats.starttime:
                    continue
                if st > tr.stats.endtime:
                    continue
                start_sample = max(st - tr.stats.starttime, 0)
                end_sample = min(et - tr.stats.starttime, tr.stats.npts)
                tr.data[start_sample:end_sample] = 0.0
        if plot:
            (outp + inp).plot(color="blue", equal_scale=False)
        if isinstance(inp, Trace):
            return outp[0]
        return outp

    def plot(self, trace=None, ax=None, show=None):
        """
        Plot representation of selected time periods as yellow highlights

        Args:
            trace (:class:`obspy.core.trace.Trace`): trace to plot along with
                selected time periods
            show (bool): show the plot on the screen
            ax (:class:`matplotlib.pyplot.Axis`): axis to plot on (plot with
                any existing plot).  If specifed, will not "show" by default
        """
        if ax is None:
            ax_existed = False
            _, ax = plt.subplots(1, 1)
            if show is None:
                show = True
        else:
            ax_existed = True
            if show is None:
                show = False
        for st, et in zip(self.start_times, self.end_times):
            ax.axvspan(st.datetime, et.datetime, facecolor="red", alpha=0.25)
        if trace is not None:
            assert isinstance(trace, Trace), "trace is not an obspy Trace"
            ax.plot(trace.times(), trace.data)
        if show:
            plt.show()


def _calc_eq_cut(mag, minmag, days_per_magnitude):
    if mag < minmag:
        return None
    return (mag - minmag) * 86400 * days_per_magnitude


def _eq_filename(starttime, endtime, minmag):
    """
    Create earthquake file filename

    Args:
        starttime (UTCDateTime): start time
        endtime (UTCDateTime): end time
        min_magmitude (float): minimum magnitude saved
    Returns:
        filename (string):
    """
    tfmt = "%Y%m%dT%H%M%S"
    return "{}-{}_MM{:g}_eqcat.qml".format(
        starttime.strftime(tfmt), endtime.strftime(tfmt), minmag
    )
