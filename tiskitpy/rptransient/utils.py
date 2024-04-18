""""
Utilities used for transient detection/removal
"""
import sys
import math as M

import numpy as np

from ..logger import init_logger

logger = init_logger()



def stack_data(data, slice_len):
    """
    Return data sliced and stacked into columns

    If the offset is non-integer, finds the best compromise for each slice
    (may lose a sample here or there)

    Args:
        data (:class:`numpy.ndarray`): 1-D array of data
        slice_len (float): the number of samples per slice (may be non-integer)
    Returns:
        (:class:`numpy.ndarray`): a 2-D array with one slice in each column
    """
    n = data.size
    nCols = int(M.floor(n/slice_len))
    nRows = int(M.floor(slice_len))
    stack = np.zeros((nRows, nCols))
    for i in range(nCols):
        off = int(np.round(i*slice_len))
        n1 = off
        n2 = n1 + nRows
        stack[:, i] = data[n1:n2]
    return stack


def prep_filter(trace, filtparms=[0.002, 2, 0.03, 6]):
    """
    Prefilter parameters used by Wielandt (Bandpass and demean)

    Args:
        trace (:class:`obspy.core.stream.Trace`): input trace
        filtparms (tuple): hipass freq, npoles, lowpass freq, npoles
    Returns:
        (:class:`obspy.core.stream.Trace`): new trace, filtered
    """
    lpf, lpc = filtparms[2:]
    hpf, hpc = filtparms[:2]

    f = trace.copy()
    f.detrend('demean')
    f.detrend('linear')
    f = f.filter('lowpass', freq=lpf, corners=lpc)
    f = f.filter('highpass', freq=hpf, corners=hpc)
    return f


def input_float(text, default):
    """
    Return a floating point number
    """
    try:
        _ = float(default)
    except Exception:
        logger.critical('The default value is not numeric!')
        sys.exit(2)
    while True:
        inp = input(f'{text}  (RETURN uses current value) [{default:g}]: ')
        if len(inp) == 0:
            return default
        try:
            f = float(inp)
            return f
        except Exception:
            print("Oops! That was no valid number. Try again...")


def _is_float_tuple(inp, nvals=None):
    if not isinstance(inp, (tuple, list)):
        print(f'input type = {type(inp)}')
        return False
    else:
        if nvals is not None:
            assert len(inp) == nvals
        try:
            _ = (float(v) for v in inp)
        except Exception:
            print(inp)
            return False
    return True


def input_float_tuple(text, default):
    """
    Return a tuple of floating point numbers
    """
    assert _is_float_tuple(default, 2), \
        'The values in the default tuple are not floating point!'
    nElements = len(default)
    while True:
        inp = input(f'{text} (RETURN for current value) {default}: ')
        if len(inp) == 0:
            return default
        try:
            newval = eval(inp)
        except Exception:
            print("Oops! That was invalid input. Try again...")
            continue
        if not _is_float_tuple(newval):
            print("You did not enter comma-separated numbers!  Try again...")
        elif len(newval) != nElements:
            print("You did not enter {:d} comma-separated #s! Try again..."
                  .format(nElements))
        else:
            break
    return newval
