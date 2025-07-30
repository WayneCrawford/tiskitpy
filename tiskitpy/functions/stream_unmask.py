"""
Copyright 2023 Wayne Crawford
"""
import numpy as np
from ..logger import init_logger

logger = init_logger()

def stream_unmask(stream):
    """
    Check if a stream is masked and, if so, unmask it
    Interpolates data in gaps
        """
    if np.any([np.ma.count_masked(tr.data) for tr in stream]):
        logger.warning('Unmasking masked data (usually a gap or overlap)')
        return stream.split().merge(fill_value='interpolate')
    return stream


if __name__ == '__main__':
    print('not a command line code')
    sys.exit(1)
