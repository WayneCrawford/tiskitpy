"""
Copyright 2026 Wayne Crawford
"""
import numpy as np
from ..logger import init_logger

logger = init_logger()


def get_stream_dtype(stream):
    """
    Get the dtype of the traces in a stream
    
    Args:
        stream (:class:`obspy.core.stream.Stream`): the input stream
        dtype (:class:`numpy.dtype`): the dtype
    Returns:
        :class:`numpy.dtype`: the dtype
    Raises:
        ValueError if not all traces have the same dtype
    """
    dtypes = [tr.data.dtype for tr in stream]
    if not len(list(set(dtypes))) == 1:
        logger.warning(f"Multiple dtypes in the same stream: {set(dtypes)}. "
                       f"Returning the first one: {dtypes[0]}")
    return dtypes[0]


def set_stream_dtype(stream, dtype):
    """
    Set all traces of a stream to the given dtype
    
    Args:
        stream (:class:`obspy.core.stream.Stream`): the input stream
        dtype (:class:`numpy.dtype`): the dtype
    Returns:
        :class:`obspy.core.stream.Stream`: the dtype-modified stream
    """
    for tr in stream:
        tr.data = tr.data.astype(dtype)
    return stream


if __name__ == '__main__':
    print('not a command line code')
    sys.exit(1)
