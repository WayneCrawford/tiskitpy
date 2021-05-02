""""
Dirac Comb routines to detect and remove glitches

Versions:
  0: straight implementation of Wielandt's Matlab routines
  1: Add features to handle earthquakes and noisy data:
    - Remove "teeth" from the Dirac comb wherever the "template" time series
      covering the same samples as the data contains a zero within the slice
      corresponding to that tooth
    - Also remove teeth for slices that have significantly higher variance than
      the others
    - Use this damaged comb and a clipped version of the data (using the clips
      variable) to generate the average glitch.
    - Apply this average glitch and the full comb to clean the data
  2:
    - In individual matching, allow further shifting of comb tooth:
        - In time
        - In height? (shouldn't allow too much liberty or may create errors)
    - Recalculate glitch based on this comb?
    - Return a sparse dirac comb for the time series, which could be combined
      with others to create a dirac comb for the whole dataset
    - Return the number of teeth used to calculate the master glitch
"""
import sys
import math as M
# import time

# import obspy.core
import numpy as np
# import scipy as sp
# import matplotlib.pyplot as plt
# from scipy.signal import convolve
# from scipy.signal import correlate, deconvolve
# from scipy.fftpack import ifft

# DEBUG = False


def stack_data(data, slice_len):
    """
    Return data sliced and stacked into columns

    If the offset is non-integer, finds the best compromise for each slice
    (may lose a sample here or there)

    :param data: 1-D array of data
    :type data: numpy.ndarray
    :param offset: the number of samples per slice (may be non-integer)
    :returns: a 2-D array with one slice in each column
    :rtype: numpy.ndarray
    """
    n = data.size
    nCols = M.floor(n/slice_len)
    nRows = M.floor(slice_len)
    stack = np.zeros((nRows, nCols))
    for i in range(nCols):
        off = np.round(i*slice_len)
        n1 = off
        n2 = n1 + nRows
        if n2 <= n:
            stack[:, i] = data[n1:n2]
    return stack


def prep_filter(trace):
    """
    Prefilter parameters used by Wielandt (Bandpass and demean)

    :param trace: input trace
    :returns: new trace, filtered
    """
    # LP 30 mHz, order 6 (dt=1.6s)
    lpf, lpc = 0.03, 6
    # HP 2 mHz, order 2 (dt=1.6s)
    hpf, hpc = 0.002, 2

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
        print('The default value is not numeric!')
        sys.exit(2)
    while True:
        inp = input(f'{text} [{default:g}]: ')
        if len(inp) == 0:
            return default
        try:
            f = float(inp)
            return f
        except Exception:
            print("Oops! That was no valid number. Try again...")


def _is_float_tuple(inp):
    if not isinstance(inp, tuple):
        return False
    else:
        try:
            _ = (float(v) for v in inp)
        except Exception:
            return False
    return True


def input_float_tuple(text, default):
    """
    Return a tuple of floating point numbers
    """
    assert _is_float_tuple(default), \
        'The values in the default tuple are not floating point!'
    nElements = len(default)
    while True:
        inp = input(f'{text} [{default}]: ')
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
