#!env python3
"""Create a template of time ranges to skip after earthquakes"""

#################################################
from dataclasses import dataclass
from pathlib import Path

# import obspy
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.core.event import Catalog, read_events
from obspy.core.stream import Stream, Trace
# import numpy as np
# import matplotlib.pyplot as plt

from .time_spans import TimeSpans


@dataclass
class EQRemover:
    """
    A class for zeroing out data after big earthquakes

    Arguments:
        start_time (UTCDateTime): earliest data that will be presented
        end_time (UTCDateTime): latest data that will be presented
        min_magnitude (float): EQ Magnitude above which to cut out times
        days_per_magnitude (float): days to cut per magnitude above
            min_magnitude
    """
    start_time: UTCDateTime
    end_time: UTCDateTime
    min_magnitude: float = 5.85
    days_per_magnitude: float = 1.5
    quiet: bool = False

    def __post_init__(self):
        self.eqfile = self.eq_filename(self.start_time, self.end_time,
                                       self.min_magnitude)
        if Path(self.eqfile).is_file():
            cat = read_events(self.eqfile, format='quakeml')
        else:
            if not self.quiet:
                print('Reading EQs from USGS online catalog...', end='',
                      flush=True)
            cat = Client("USGS").get_events(
                starttime=self.start_time - self._calc_cut(9),
                endtime=self.end_time,
                minmagnitude=self.min_magnitude,
                orderby='time-asc')
            if not self.quiet:
                print('Done', flush=True)
            print(f'writing catalog to "{self.eqfile}"')
            cat.write(self.eqfile, format='quakeml')

        new_cat = Catalog(events=[x for x in cat
                                  if x.preferred_magnitude().mag
                                  >= self.min_magnitude])
        start_times = [x.preferred_origin().time for x in new_cat]
        end_times = [x.preferred_origin().time
                     + self._calc_cut(x.preferred_magnitude().mag)
                     for x in new_cat]
        self.spans = TimeSpans(start_times, end_times)

    def _calc_cut(self, mag):
        if mag < self.min_magnitude:
            return None
        return (mag - self.min_magnitude) * 86400 * self.days_per_magnitude

    def eq_filename(self, starttime, endtime, min_magnitude, format='quakeml'):
        """
        Create earthquake file filename
        
        Args:
            starttime (UTCDateTime): start time
            endtime (UTCDateTime): end time
            min_magmitude (float): minimum magnitude saved
            format (str): file format (only 'quakeml')
        Returns:
            filename (string):
        """
        if format == 'quakeml':
            return '{}-{}_MM{:g}_eqcat.qml'.format(
                starttime.strftime('%Y%m%dT%H%M%S'),
                endtime.strftime('%Y%m%dT%H%M%S'),
                min_magnitude)
        raise ValueError('format is not "quakeml"')

    def add_timespans(self, time_spans):
        """
        Add time spans to existing EQRemover object
        
        Args:
            time_spans (~class `.TimeSpans`): time spans to add
        """
        self.spans.append(time_spans)

    def has_zeros(self, starttime, endtime):
        """
        Does the trace's time span intersect an EQ zero?

        Arguments:
            starttime (UTCDateTime): start time
            endtime (UTCDateTime): end time

        Returns:
            (bool):
        """
        return self.spans.has_zeros(starttime, endtime)

    def zero(self, inp, plot=False):
        """
        Zero out data in trace  or stream when it is in a flagged time span

        Arguments:
            inp (obspy.core.Trace): seismological data
            plot: plot traces with time spans set to zero

        Returns:
            trace with time spans set to zero
        """
        if not isinstance(inp, Trace) and not isinstance(inp, Stream):
            raise(ValueError, 'inp is not a Trace or Stream object')
        out = self.spans.zero(inp, plot=plot)
        return out
