"""
Tools for processing SeisComp Data Structure (SDS) archives file by file,
using neighbouring data to keep the results continuous across files
"""
import re
from pathlib import Path

import numpy as np
from obspy import UTCDateTime
from obspy.core.stream import read as stream_read

# SDS file name: NET.STA.LOC.CHAN.TYPE.YEAR.DAY
SDS_FNAME = re.compile(r'^(?P<net>[^.]*)\.(?P<sta>[^.]*)\.(?P<loc>[^.]*)'
                       r'\.(?P<chan>[^.]*)\.(?P<type>[^.])'
                       r'\.(?P<year>\d{4})\.(?P<day>\d{3})$')


def sds_channel_files(SDS_root):
    """
    Group the files of an SDS archive by channel (across years)

    Args:
        SDS_root (str or Path): SDS root directory
    Returns:
        (dict): keys are (net, sta, loc, chan, type) tuples, values are
            lists of (filedate, file_path, match) tuples, sorted by date,
            where match is the `SDS_FNAME` match of the file name
    """
    channels = {}
    for file_path in Path(SDS_root).glob('*/*/*/*/*'):
        m = SDS_FNAME.match(file_path.name)
        if m is None or not file_path.is_file():
            continue
        filedate = UTCDateTime(year=int(m['year']), julday=int(m['day']))
        key = (m['net'], m['sta'], m['loc'], m['chan'], m['type'])
        channels.setdefault(key, []).append((filedate, file_path, m))
    for key in channels:
        channels[key].sort(key=lambda x: x[0])
    return channels


def sds_file_path(SDS_root, match, loc=None, chan=None):
    """
    Path of an SDS file, optionally with a new location and/or channel code

    Args:
        SDS_root (str or Path): SDS root directory
        match: `SDS_FNAME` match of the original file name
        loc (str or None): new location code (None: keep original)
        chan (str or None): new channel code (None: keep original)
    """
    m = match
    loc = m['loc'] if loc is None else loc
    chan = m['chan'] if chan is None else chan
    return (Path(SDS_root) / m['year'] / m['net'] / m['sta']
            / f"{chan}.{m['type']}"
            / '.'.join((m['net'], m['sta'], loc, chan, m['type'], m['year'],
                        m['day'])))


def sds_neighbour_streams(files, startdate=None, enddate=None, verbose=False):
    """
    Read the files of one SDS channel, with their neighbouring files

    Files without waveform data (e.g. ASCII log records) are skipped.
    Neighbouring files are returned even if they are outside of the date
    range.  Each file is read only once.

    Args:
        files (list): one value of the `sds_channel_files()` dict
        startdate (UTCDateTime or None): skip files starting before this
        enddate (UTCDateTime or None): skip files ending after this
        verbose (bool): print the name of each file
    Yields:
        (tuple): match, previous stream (or None), stream, next stream (or
            None)
    """
    streams = {}
    for i, (filedate, file_path, m) in enumerate(files):
        if startdate is not None and filedate < startdate:
            continue
        if enddate is not None and filedate + 86400 > enddate:
            continue
        if verbose:
            print(f'{file_path=}')
        window = [x[1] for x in files[max(i - 1, 0): i + 2]]
        for path in window:
            if path not in streams:
                streams[path] = stream_read(str(path))
        streams = {p: streams[p] for p in window}
        if not has_waveforms(streams[file_path]):
            print(f'{file_path.name}: no waveform data, skipped')
            continue
        prev_st = streams[files[i - 1][1]] if i > 0 else None
        next_st = streams[files[i + 1][1]] if i < len(files) - 1 else None
        yield m, prev_st, streams[file_path], next_st


def has_waveforms(stream):
    """True if any trace in the stream has samples at a non-zero rate"""
    return any(tr.stats.sampling_rate > 0 and tr.stats.npts > 0
               for tr in stream)


def neighbour_gap(trace, neighbour, side, max_misalignment=0.5):
    """
    Number of missing samples between a trace and a neighbouring trace

    Args:
        trace (:class:`obspy.core.trace.Trace`): reference trace
        neighbour (:class:`obspy.core.trace.Trace`): trace before or after
        side (str): 'before' or 'after': which side the neighbour is on
        max_misalignment (float): largest offset (in samples) from the
            trace's sample grid for which the neighbour is considered on
            the grid
    Returns:
        (tuple):
            aligned (bool): neighbour is on the trace's sample grid
            positions (np.ndarray): neighbour sample positions (in samples)
                relative to the trace's first sample (rounded if aligned)
            gap (float): missing samples between the trace and the neighbour
                (0: consecutive, < 0: overlap)
    """
    sr = trace.stats.sampling_rate
    offset = (neighbour.stats.starttime - trace.stats.starttime) * sr
    aligned = abs(offset - round(offset)) <= max_misalignment
    if aligned:
        offset = round(offset)
    pos = offset + np.arange(neighbour.stats.npts)
    if side == 'before':
        gap = -1 - pos[-1]
    else:
        gap = pos[0] - trace.stats.npts
    return aligned, pos, gap


def is_consecutive(prev_trace, trace, max_misalignment=0.5):
    """
    True if trace follows prev_trace with no gap or overlap

    Sample offsets of up to max_misalignment samples from the trace's
    sample grid are accepted (see `neighbour_gap()`)
    """
    if not (prev_trace.id == trace.id and prev_trace.stats.sampling_rate
            == trace.stats.sampling_rate):
        return False
    aligned, _, gap = neighbour_gap(trace, prev_trace, 'before',
                                    max_misalignment)
    return aligned and gap == 0


def get_neighbour_pad(trace, neighbours, n, side, gap_handling='interpolate',
                      max_gap=0, max_misalignment=0.5):
    """
    Return n samples preceding or following a trace, from its neighbours

    Args:
        trace (:class:`obspy.core.trace.Trace`): trace to pad
        neighbours (list of :class:`obspy.core.trace.Trace`): traces with
            the same id and sampling rate
        n (int): number of samples
        side (str): 'before' or 'after'
        gap_handling (str): 'interpolate' or 'edge'
        max_gap (int): longest gap (in samples) to interpolate across
        max_misalignment (float): largest offset (in samples) from the
            trace's sample grid for which the neighbour is shifted onto
            the grid
    Returns:
        (np.ndarray or None): pad, or None if there is no valid neighbour.
            If the neighbour is shorter than n, its last value is
            repeated.
    """
    npts = trace.stats.npts
    if side == 'before':
        cands = [x for x in neighbours
                 if x.stats.starttime < trace.stats.starttime]
        if len(cands) == 0:
            return None
        nb = max(cands, key=lambda x: x.stats.endtime)
        want = np.arange(-n, 0)
    else:
        cands = [x for x in neighbours
                 if x.stats.endtime > trace.stats.endtime]
        if len(cands) == 0:
            return None
        nb = min(cands, key=lambda x: x.stats.starttime)
        want = np.arange(npts, npts + n)

    # Neighbour sample positions, relative to the trace's first sample
    aligned, pos, _ = neighbour_gap(trace, nb, side, max_misalignment)
    nb_data = nb.data.astype(float)

    # Keep neighbour samples outside of the trace, add the trace's edge
    # sample (used to interpolate across a gap)
    if side == 'before':
        keep = pos < -0.5
        if not np.any(keep):
            return None
        xp = np.append(pos[keep], 0)
        fp = np.append(nb_data[keep], trace.data[0])
        gap = -1 - pos[keep][-1]
    else:
        keep = pos > npts - 0.5
        if not np.any(keep):
            return None
        xp = np.insert(pos[keep], 0, npts - 1)
        fp = np.insert(nb_data[keep], 0, trace.data[-1])
        gap = pos[keep][0] - npts

    if not (aligned and gap == 0):
        if gap_handling == 'edge' or gap > max_gap:
            print(f'{trace.id}: {side} neighbour not consecutive '
                  f'({gap=:g} samples), accepting edge effects')
            return None
    return np.interp(want, xp, fp)
