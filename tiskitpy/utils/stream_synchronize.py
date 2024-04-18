"""
Copyright 2023 Wayne Crawford
"""
from ..logger import init_logger

logger = init_logger()


def stream_synchronize(in_stream, max_reject=0.5):
    """
    Ensures that all traces in a stream are same sample rate, starttime
    and number of samples.  If stream is already synchronized, returns the
    original stream object

    Args:
        in_stream (:class:`obspy.core.Stream`): stream
        max_reject (float): maximum portion to accept rejecting (0-1)
    Returns:
        stream (:class:`obspy.core.Stream`): stream
    """
    if _is_synced(in_stream) is True:
        return in_stream

    if max_reject < 0 or max_reject > 1.:
        raise ValueError('max_reject must be between 0 and 1')

    # Try merging the input stream
    out_stream = in_stream.copy()
    out_stream.merge(0, fill_value='interpolate')

    if _is_synced(out_stream) is True:
        return out_stream

    # Synchronize start times
    srate = _get_sampling_rate(out_stream)
    starts = [x.stats.starttime for x in out_stream]
    if srate is False:
        raise ValueError("stream has more than one sampling rate")
    last_start = max(starts)
    for tr in out_stream:
        if last_start - tr.stats.starttime > 1 / srate:
            offset = int((last_start - tr.stats.starttime)*srate)
            tr.data = tr.data[offset:]
            tr.stats.starttime += offset/srate
    # Equalize data lengths
    min_nsamps = min([len(x) for x in out_stream])
    for tr in out_stream:
        if len(tr) > min_nsamps:
            tr.data = tr.data[:min_nsamps]
    # Verify that did not throw out too much
    for tr_in, tr_out in zip(in_stream, out_stream):
        if len(tr_out)/len(tr_in) < max_reject:
            logger.warning('Cancel synchronization, would need to cut more than {max_reject*100}%')
            return None

    return out_stream

def _is_synced(stream):
    """
    Returns True is stream is synchronized, False otherwise
    """
    # If only 1 (or zero) channels, it's synced!
    if len(stream) < 2:
        return True

    # All channels must have the same sampling rate
    srate = _get_sampling_rate(stream)
    if srate is False:
        return False

    starts = [x.stats.starttime for x in stream]
    nsamps = [len(x) for x in stream]
    for start, nsamp in zip(starts[1:], nsamps[1:]):
        if abs(start-starts[0]) > 1 / srate:
            return False
        if not nsamp == nsamps[0]:
            return False
    return True

def _get_sampling_rate(stream):
    """
    return the stream's sampling rate, or False if there is more than one
    """
    # All channels must have the same sampling rate
    srate = stream[0].stats.sampling_rate
    for sr in [x.stats.sampling_rate for x in stream[1:]]:
        if not sr == srate:
            return False
    return srate
