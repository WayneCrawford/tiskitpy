#!env python3
"""Class of time spans to remove, keep, zero, etc. in Trace or Stream data"""
from pathlib import Path

from obspy.clients.fdsn import Client
from obspy.core.event import Catalog, read_events
from obspy.core import UTCDateTime
from obspy.core.stream import Stream, Trace
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.dates import date2num

from .utils import CleanSequence
from .cleaned_stream import CleanedStream
from .logger import init_logger

logger = init_logger()
ZERO_TAG = 'ZEROS'

class TimeSpans:
    """
    A class specifying time spans, to be removed, kept, zeroed, etc.
    
    Times are obspy UTCDateTimes
    
    You can initialize using:
    
        A list of (start_time, end_time) pairs (start_time and end_time can
        be any input to UTCDateTime, including UTCDateTime)
    - or -
        A list of UTCDateTime start_times and a list of UTCDateTime end_times
        (the lists must be the same length)
    """
    def __init__(self, spans: list=None, start_times: list=None, end_times: list=None):
        if spans is not None:
            if start_times is not None or end_times is not None:
                raise ValueError('You provided spans AND start_times '
                                 'and/or end_times')
            try:
                self._start_times = [UTCDateTime(x[0]) for x in spans]
                self._end_times = [UTCDateTime(x[1]) for x in spans]
            except Exception as err:
                raise TypeError('Could not convert at least one of the span '
                                f'input values to UTCDateTime')
        elif start_times is not None and end_times is not None:
            if not len(start_times) == len(end_times):
                raise ValueError(
                    f'{len(start_times)=} != {len(end_times)=}')
            try:
                self._start_times = [UTCDateTime(x) for x in start_times]
            except Exception as err:
                raise TypeError('Could not convert at least one of the '
                                f'start_times to UTCDateTime(): {err}')
            try:
                self._end_times = [UTCDateTime(x) for x in end_times]
            except Exception as err:
                raise TypeError('Could not convert at least one of the '
                                f'end_times to UTCDateTime(): {err}')
        else:
            raise ValueError('You must provide spans or start_times '
                                 'and end_times')
        self._organize()

    @property
    def start_times(self):
        """list of span start times"""
        return self._start_times.copy()

    @property
    def end_times(self):
        """list of span start times"""
        return self._end_times.copy()

    @property
    def spans(self):
        """list of spans [start_time, end_time]"""
        return [[x, y] for x, y in zip(self._start_times, self._end_times)]

    @classmethod
    def from_eqs(cls, starttime, endtime, minmag=5.85, days_per_magnitude=1.5,
                 eq_file=None, save_eq_file=True):
        """
        Generate timespans to avoid because of earthquakes

        Will read earthquakes from the USGS online catalog the first time,
        saving the information to a file that can be subsequently used

        Args:
            starttime (:class:`UTCDateTime` or str): earliest data that will be
                presented.  If a str, must by ISO8601 compatible.  Forced
                to the beginning of the day
            endtime (:class:`UTCDateTime` or str): latest data that will be presented.
                If a str, must by ISO8601 compatible.  Forced to the
                end of the day
            minmag (float): EQ Magnitude above which to cut out times
            days_per_magnitude (float): days to cut per magnitude above
                min_magnitude
            eq_file (str): the eq filename (otherwise, generates it)
            save_eq_file (bool): save the catalog file for future use
        Returns:
            eq_spans (:class:`TimeSpans`): time spans covering EQ signal
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
        starttime = starttime.replace(hour=0, minute=0, second=0, microsecond=0)
        endtime = endtime.replace(hour=23, minute=59, second=59, microsecond=999999)
        if eq_file is None:
            eq_file = _eq_filename(starttime, endtime, minmag)
        if Path(eq_file).is_file():
            cat = read_events(eq_file, format="quakeml")
        else:
            logger.info(f"Didn't find local EQ file '{eq_file}', reading from USGS online catalog...")
            try:
                cat = Client("USGS").get_events(
                    starttime=starttime
                    - _calc_eq_cut(9, minmag, days_per_magnitude),
                    endtime=endtime,
                    minmagnitude=minmag,
                    orderby="time-asc",
                )
                logger.info("Done")
                logger.info(f'writing catalog to "{eq_file}"')
                if save_eq_file:
                    cat.write(eq_file, format="quakeml")
            except Exception as e:  # except FDSNNoServiceException as e:
                logger.warning(e)
                logger.warning('!!!Continuing without removing EQs!!!')
                return cls([])
            

        new_cat = Catalog(
            events=[x for x in cat if x.preferred_magnitude().mag >= minmag]
        )
        spans = [[x.preferred_origin().time,
                  x.preferred_origin().time
                  + _calc_eq_cut(x.preferred_magnitude().mag,
                                 minmag, days_per_magnitude)]
                  for x in new_cat if _calc_eq_cut(x.preferred_magnitude().mag,
                                                   minmag, days_per_magnitude)
                                   > 0]
        return cls(spans)

    # INFORMATION METHODS
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

    def __repr__(self):
        return f"TimeSpans({len(self._start_times):d} spans)"

    def __str__(self):
        s =  "TimeSpans: start              |            end\n"
        s += "==============================+===============================\n"
        for st, et in zip(self.start_times, self.end_times):
            s += f" {str(st):28} | {str(et):28}\n"
        return s

    def __add__(self, other):
        """
        Appends TimeSpan object to self

        Args:
            other (~class `TimeSpans`): time spans to append
        """
        if not isinstance(other, TimeSpans):
            raise TypeError(f"Tried to add a {type(other)} to a TimeSpans")
        return TimeSpans(start_times=self.start_times + other.start_times,
                         end_times=self.end_times + other.end_times)

    def __radd__(self, other):  # Allows us to use sum()
        if other == 0:
            return self
        else:
            return self.__add__(other)

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

    def _validate(self):
        """
        Raises error if there is a problem, returns false if not organized

        Raises:
            (ValueError): if
                - len(start_times) ~= len(end_times)
                - not every object in start_times and end_times is a UTCDateTime
                - start_times[i] >= end_times[i]
        Returns:
            (bool): False if start times are not ordered, or there are overlaps
                 between time spans.  True otherwise
        """
        # Look for errors
        if not len(self.start_times) == len(self.end_times):
            raise ValueError(
                f"{len(self.start_times)=} != {len(self.end_times)=}")
        for x in self.start_times:
            if not isinstance(x, UTCDateTime):
                raise ValueError("There is a non-UTCDateTime starttime value")
        for x in self.end_times:
            if not isinstance(x, UTCDateTime):
                raise ValueError("There is a non-UTCDateTime endtime value")
        for (s, e) in self.spans:
            if s >= e:
                raise ValueError(f"start_time >= end_time ({s}, {e})")
        
        # Look for unorganized data
        # start_times not sorted
        if not sorted(self.start_times) == self.start_times:
            return False
        # overlapping windows
        for startafter, endbefore in zip(self.start_times[1:],
                                         self.end_times[:-2]):
            if startafter <= endbefore:
                return False
        return True

    # METHODS OPERATING ON THE TIMESPANS OBJECT
    def invert(self, ts_starttime, ts_endtime):
        """
        Return inverted time spans

        If input TimeSpan had spans from A->A'', B->B'and C->C'', output
        TimeSpan will have spans from ts_starttime->A, A'->B, B'->C and
        C'->ts_endtime  if ts_starttime < A and ts_endtime > B.
        If ts_start_time -> ts_endtime, is outside of the time spans, will
        simply use them for the time span

        Args:
            ts_startime (:class:`obspy.core.UTCDateTime`): timeseries start
                time.  If None, inverted TimeSpans will start at the first
                end_time
            ts_startime (:class:`obspy.core.UTCDateTime`): timeseries end time.
                If None, inverted TimeSpans will end at the last start_time
        Returns:
            (:class:`TimeSpans`)
        """
        self._organize()
        
        if isinstance(ts_starttime, str):
            try:
                ts_starttime = UTCDateTime(ts_starttime)
            except Exception:
                raise RunTimeError(f'Could not convert  {ts_starttime=} to UTCDateTime')
        elif not isinstance(ts_starttime, (UTCDateTime, type(None))):
            raise TypeError('ts_starttime is not a UTCDateTime, str or None')
        if isinstance(ts_endtime, str):
            try:
                ts_endtime = UTCDateTime(ts_endtime)
            except Exception:
                raise RunTimeError(f'Could not convert {ts_endtime=} to UTCDateTime')
        if not isinstance(ts_endtime, (UTCDateTime, type(None))):
            raise TypeError('ts_endtime is not a UTCDateTime, str or None')
        
        if ts_starttime is None:
            ts_starttime = self.end_times[0]
        if ts_endtime is None:
            ts_endtime = self.start_times[-1]
        if ts_starttime >= ts_endtime:
            raise ValueError('ts_starttime > ts_endtime')
        new_spans = [x for x in self.spans if x[1] > ts_starttime and x[0] < ts_endtime]
        if len(new_spans) == 0:
            if ts_starttime > self.end_times[-1]:
                logger.info(f'{ts_starttime=} after end of TimeSpans')
            elif ts_endtime < self.start_times[0]:
                logger.info(f'{ts_endtime=} before start of TimeSpans')
            return TimeSpans([[ts_starttime, ts_endtime]])
        if ts_starttime < new_spans[0][0]:
            new_spans = [[ts_starttime, ts_starttime]] + new_spans
        if ts_endtime > new_spans[-1][1]:
            new_spans.append([ts_endtime, ts_endtime])
        new_st = [x[1] for x in new_spans[:-1]]
        new_et = [x[0] for x in new_spans[1:]]

        return TimeSpans(start_times=new_st, end_times=new_et)

    def append(self, new_time_spans):
        """
        Combines two TimeSpan objects (REPLACED BY COMBINE)

        Args:
            new_time_spans (:class:`TimeSpans`): object to combine with self
        """
        logger.warning('TimeSpans.append() has been renamed '
                       'TimeSpans.combine(), change your code before append()'
                       ' is deleted!')
        self.combine(new_time_spans)

    def combine(self, new_time_spans):
        """
        Combines two TimeSpan objects

        Args:
            new_time_spans (:class:`TimeSpans`): object to combine with self
        """
        new_time_spans._validate()
        self._validate()
        self._start_times.extend(new_time_spans.start_times)
        self._end_times.extend(new_time_spans.end_times)
        self._organize()

    def _organize(self):
        """
        Order starttimes and endtimes by increasing starttime, and consolidate
        overlapping time spans
        """
        if self._validate() is True:
            return

        # sort by time
        self._end_times = [
            x for _, x in sorted(zip(self.start_times, self.end_times))
        ]
        self._start_times = sorted(self.start_times)

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
        self._start_times = start_times
        self._end_times = end_times

    # METHODS OPERATING ON OUTSIDE OBJECTS
    def zero(self, inp, plot=False):
        """
        Zero out data in the time spans

        Arguments:
            inp (Trace or Stream): seismological data
            plot: plot traces with spans cut out

        Returns:
            Trace or Stream with spans set to zero
        """
        if isinstance(inp, Trace):
            tr = inp.copy()  # Do not destroy original
            for st, et in zip(self.start_times, self.end_times):
                start_addr, end_addr = self._get_addrs(st, et, tr.stats)
                if start_addr is not None:
                    tr.data[start_addr: end_addr + 1] = 0.0
                    tr = CleanSequence.tag(tr, ZERO_TAG)
                    # tr.stats.channel = "XX" + tr.stats.channel[2]  # Mark channel code
            if plot:
                Stream([trace,inp]).plot(color="blue", equal_scale=False)
            return tr
        elif isinstance(inp, Stream):
            stream = CleanedStream(inp)
            for tr in stream:
                new_tr = self.zero(tr)
                stream.remove(tr)
                stream.append(new_tr)
            # stream = stream.tag(ZERO_TAG)
            if plot:
                (stream + inp).plot(color="blue", equal_scale=False)
            return stream
        else:
            raise (ValueError, "inp is not an obspy Trace or Stream")

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
            trace with interpolated data in each time span
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
                tr.data[start_addr: end_addr + 1] = np.linspace(
                    tr.data[start_addr],
                    tr.data[end_addr],
                    end_addr - start_addr + 1,
                )
        if plot:
            (stream + inp).plot(color="blue", equal_scale=False)
        if isinstance(inp, Trace):
            return stream[0]
        return stream

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

    def plot(self, stream=None, color="red", alpha=0.25, title=None, **kwargs):
        """
        Make a stream or trace plot with highlighted time spans

        Args:
            stream (:class:`obspy.core.trace.Stream` or :class:`obspy.core.trace.Trace`):
                obspy stream/trace to plot
            color (str): highlight color
            alpha (float): highlight transparency alpha (1=opaque, 0 =
                invisible)
            title (str): figure title
            kwargs (dict): arguments to stream/trace plot program
        Returns:
            (tuple):
                fig (matplotlib Figure):
                ax (list): list of axes
        """
        if alpha < 0:
            raise ValueError("alpha < 0")
        elif alpha > 1:
            raise ValueError("alpha > 1")
        # Delay plotting or saving until after we decorate the graph
        outfile = kwargs.pop('outfile', None)
        show = kwargs.pop('show', True)
        block = kwargs.pop('block', True)
        kwargs['handle'] = True

        if stream is not None:
            # Plot the stream or trace
            fig = stream.plot(**kwargs)
        else:
            fig, ax = plt.subplots()
            ax.xaxis.axis_date()

        # Add colors for time spans
        for ax in fig.get_axes():
            for st, et in zip(self.start_times, self.end_times):
                xmin = date2num(st)
                xmax = date2num(et)
                ax.axvspan(xmin, xmax, facecolor=color, alpha=0.25)
        
        if title is not None:
            fig.suptitle(title)

        # Plot to screen or save to file
        if outfile:
            fig.savefig(outfile)
        elif show is True:
            try:
                plt.show(block=block)
            except Exception:
                plt.show()
        return fig, fig.get_axes()

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
        filename (string): includes startday to endday information
    """
    tfmt = "%Y%m%d"
    return "{}-{}_MM{:g}_eqcat.qml".format(
        starttime.strftime(tfmt), endtime.strftime(tfmt), minmag
    )
