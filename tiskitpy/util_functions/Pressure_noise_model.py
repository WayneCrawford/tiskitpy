#!/usr/bin/env python3
"""
Pressure Low and High noise models

From Brown et al [2014] (IDC_2010_LH and _HH).
Should add Webb LF bounds
"""
import sys

import numpy as np

from ..logger import init_logger

logger = init_logger()


#     Freq      dBs ref uPa2/Hz
LH = [[0.001,   162.5],
      [0.005,   145],
      [0.01,    135],
      [0.02,    115],
      [0.03,    105],
      [0.05,    100],
      [0.1,     100],
      [0.15,    110.5],
      [0.2,     130],
      [0.3,     121],
      [0.4,     117.5],
      [0.5,     112],
      [1.0,     95],
      [1.5,     91],
      [2.0,     90],
      [3.0,     81],
      [4.0,     76],
      [5.0,     72],
      [6.0,     71],
      [8.0,     72.5],
      [10.0,    71.5],
      [20.0,    71.0],
      [30.0,    73.0],
      [40.0,    74.0],
      [50.0,    73.0],
      [60.0,    72.0],
      [70.0,    70.5],
      [80.0,    68.5],
      [100.0,   63]]

HH = [[0.001,   172],
      [0.002,   166],
      [0.003,   166],
      [0.004,   168],
      [0.006,   168.5],
      [0.01,    167],
      [0.013,   160],
      [0.02,    145],
      [0.03,    119],
      [0.05,    117],
      [0.1,     125],
      [0.11,    139],
      [0.13,    142.5],
      [0.18,    156],
      [0.2,     155],
      [0.3,     149.5],
      [0.4,     143.5],
      [0.5,     137.5],
      [1.0,     119],
      [2.0,     106],
      [3.0,     102],
      [4.0,     100.5],
      [5.0,     100],
      [10.0,    95.5],
      [13.0,    93],
      [15.0,    99],
      [20.0,    96],
      [30.0,    91],
      [40.0,    90],
      [50.0,    89.5],
      [60.0,    89],
      [70.0,    87.5],
      [80.0,    87],
      [90.0,    84],
      [100.0,   81]]


def Pressure_noise_model(periods, as_freqs=False):
    """
    Return Brown et al [2014] low and high marine pressure noise models

    Args:
        periods (list): periods to use (should be increasing).
        as_freqs (bool): interpret "periods" as frequencies instead
    Returns:
        tuple: (lownoise, highnoise) values in dB ref to 1 (m/s^2)^2/Hz
    """
    if as_freqs:
        freqs = periods
        lownoise = _fit_points(freqs, LH)
        highnoise = _fit_points(freqs, HH)
    else:
        freqs = np.power(periods[::-1], -1)
        lownoise = _fit_points(freqs, LH)
        highnoise = _fit_points(freqs, HH)
        lownoise = lownoise[::-1]
        highnoise = highnoise[::-1]
    return lownoise, highnoise


def _fit_points(freqs, model):
    """
    Fit points to a noise model

    :param freqs: freqs in increasing order
    :type freqs: list
    :param model: list of [freq, value]
    :type model: list of lists
    """
    x = np.log10(freqs)
    xp = np.log10([x[0] for x in model])
    yp = [x[1]-120 for x in model] # Convert to ref Pa^2/Hz
    assert np.all(np.diff(x) > 0), 'x is not increasing'
    assert np.all(np.diff(xp) > 0), 'xp is not increasing'
    return np.interp(x, xp, yp, left=np.nan, right=np.nan)


if __name__ == '__main__':
    print('not a command line code')
    sys.exit(1)
