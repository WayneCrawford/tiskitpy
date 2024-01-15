"""
Remove information about the cleaning sequence from the location code of
a Trace, Stream, or SEED id string
"""
from obspy.core import Trace, Stream

from .cleaner_string import CleanerString
from ..logger import init_logger

logger = init_logger()



def remove_cleaner_string(input):
    """
    Remove a cleaning sequence (everything starting with "-") from a Trace,
    Stream, or SEED id's location code
    
    Args:
        input (str, :class:`Trace`, or :class:`Stream`): input string, Trace,
            or Stream
    Returns:
        output (str, :class:`Trace`, or :class:`Stream`)
    """
    if isinstance(input, str):
        return CleanerString.strip(input)
    elif isinstance(input, Trace):
        new_trace = input.copy()
        new_trace.id = CleanerString.strip(new_trace.get_id())
        return new_trace
    elif isinstance(input, Stream):
        new_stream = input.copy()
        for tr in new_stream:
            new_tr = remove_cleaner_string(tr)
            new_stream.remove(tr)
            new_stream += new_tr
        return new_stream
    else:
        raise TypeError('input is a "{type(input)}", must be a str, Trace or Stream')