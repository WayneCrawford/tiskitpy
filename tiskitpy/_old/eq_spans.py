#!env python3
"""Create a template of time spans to skip after earthquakes"""

#################################################
from pathlib import Path

# import obspy
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog, read_events
# import numpy as np
# import matplotlib.pyplot as plt

from .time_spans import TimeSpans


def get_eq_spans(starttime, endtime, minmag=5.85,
                 days_per_magnitude=1.5, quiet=False):
    """
    Return timespans to avoid because of earthquakes
    
    Will read earthquakes from the USGS online catalog the first time, saving
    the information to a file that can be subsequently used

    Args:
        starttime (UTCDateTime): earliest data that will be presented
        endtime (UTCDateTime): latest data that will be presented
        minmag (float): EQ Magnitude above which to cut out times
        days_per_magnitude (float): days to cut per magnitude above
            min_magnitude
    Returns:
        eq_spans (~class `.TimeSpans`): time spans covering EQ signal
    """
    eqfile = eq_filename(starttime, endtime, minmag)
    if Path(eqfile).is_file():
        cat = read_events(eqfile, format='quakeml')
    else:
        if not quiet:
            print('Reading EQs from USGS online catalog...', end='',
                  flush=True)
        cat = Client("USGS").get_events(
            starttime=starttime - _calc_cut(9, minmag, days_per_magnitude),
            endtime=endtime,
            minmagnitude=minmag,
            orderby='time-asc')
        if not quiet:
            print('Done', flush=True)
        print(f'writing catalog to "{eqfile}"')
        cat.write(eqfile, format='quakeml')

    new_cat = Catalog(events=[x for x in cat
                              if x.preferred_magnitude().mag >= minmag])
    start_times = [x.preferred_origin().time for x in new_cat]
    end_times = [x.preferred_origin().time
                 + _calc_cut(x.preferred_magnitude().mag, minmag,
                             days_per_magnitude)
                 for x in new_cat]
    return TimeSpans(start_times, end_times)


def _eq_filename(starttime, endtime, minmag, format='quakeml'):
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
            minmag)
    raise ValueError('format is not "quakeml"')


def _calc_cut(mag, minmag, days_per_magnitude):
    if mag < minmag:
        return None
    return (mag - minmag) * 86400 * days_per_magnitude
