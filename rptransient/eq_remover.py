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
import numpy as np
import matplotlib.pyplot as plt

from .time_spans import TimeSpans


@dataclass
class EQRemover:
    """
    A class for zeroing out data after big earthquakes

    Arguments:
        start_time: earliest data that will be presented
        end_time: latest data that will be presented
        min_magnitude (float): EQ Magnitude above which to cut out times
        days_per_magnitude (float): days to cut per magnitude above
            min_magnitude
        eqfile: name of  csv file containing global earthquakes
            (column 1 = time (ISO8601), column4=magnitude)
            If `None`, will download a catalog starting early enough to cover
            a M9 EQ from the USGS FDSN server and save this catalog to a
            QuakeML file, or it will read this previously created file if all
            other arguments are the same.
    """
    start_time: UTCDateTime
    end_time: UTCDateTime
    min_magnitude: float = 5.85
    days_per_magnitude: float = 1.5
    quiet: bool = False
    eqfile: str = ''

    def __post_init__(self):
        if not self.eqfile:
            self.eqfile = '{}-{}_M{:g}_DM{:g}_eqcat.qml'.format(
                self.start_time.strftime('%Y%m%dT%H%M%S'),
                self.end_time.strftime('%Y%m%dT%H%M%S'),
                self.min_magnitude, self.days_per_magnitude)
            if Path(self.eqfile).is_file():
                cat = read_events(self.eqfile, format='quakeml')
            else:
                if not self.quiet:
                    print('Reading EQs from USGS online catalog...', end = '',
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
        else:
            cat = read_events(self.eqfile, format='csv')

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

    def zero(self, inp, plot=False):
        """
        Zero out data in trace  or stream when it is too close to an eq

        Returns an error if the EQ catalog does not bound the Trace time
        Returns a warning if the EQ catalog does not start far enough before
        the trace to cover an M9 EQ.

        Arguments:
            inp (obspy.core.Trace): seismological data
            plot: plot traces with EQs cut out

        Returns:
            trace with EQ-affected sections set to zero
        """
        if not isinstance(inp, Trace) and not isinstance(inp, Stream):
            raise(ValueError, 'inp is not a Trace or Stream object')
        out = self.spans.zero(inp, plot=plot)
        return out