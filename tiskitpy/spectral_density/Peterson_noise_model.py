#!/usr/bin/env python3
"""
Peterson Low and High noise models
"""
import sys

import numpy as np

from ..logger import init_logger

logger = init_logger()


#        Period    Level       Slope
LPAB = [[0.10,    -162.36,     5.64],
        [0.17,    -166.7,      0.00],
        [0.40,    -170.0,     -8.30],
        [0.80,    -166.4,     28.90],
        [1.24,    -168.60,    52.48],
        [2.40,    -159.98,    29.81],
        [4.30,    -141.10,     0.00],
        [5.00,     -71.36,   -99.77],
        [6.00,     -97.26,   -66.49],
        [10.00,   -132.18,   -31.57],
        [12.00,   -205.27,    36.16],
        [15.60,    -37.65,  -104.33],
        [21.90,   -114.37,   -47.10],
        [31.60,   -160.58,   -16.28],
        [45.00,   -187.50,     0.00],
        [70.00,   -216.47,    15.70],
        [101.00,  -185.00,     0.00],
        [154.00,  -168.34,    -7.61],
        [328.00,  -217.43,    11.90],
        [600.00,  -258.28,    26.60],
        [10000.0, -346.88,    48.75],
        [100000,  -346.88,    48.75]]

HPAB = [[0.10,    -108.73,   -17.23],
        [0.22,    -150.34,   -80.50],
        [0.32,    -122.31,   -23.87],
        [0.80,    -116.85,    32.51],
        [3.80,    -108.48,    18.08],
        [4.60,     -74.66,   -32.95],
        [6.30,       0.66,  -127.18],
        [7.90,     -93.37,   -22.42],
        [15.40,     73.54,  -162.98],
        [20.00,   -151.52,    10.01],
        [354.80,  -206.66,    31.63],
        [100000,  -206.66,    31.63]]


def Peterson_noise_model(periods, as_freqs=False):
    """
    Return Peterson low and high seismological noise models

    returns the acceleration noise models in dB ref to 1 (m/s^2)^2/Hz

    :param periods: periods to use (should be increasing).
    :type freqs: list
    :param as_freqs: interpret "periods" as frequencies instead
    :type as_freqs: bool, opt
    """
    if not as_freqs:
        lownoise = _fit_points(periods, LPAB)
        highnoise = _fit_points(periods, HPAB)
    else:
        periods = np.power(periods[::-1], -1)
        lownoise = _fit_points(periods, LPAB)
        highnoise = _fit_points(periods, HPAB)
        lownoise = lownoise[::-1]
        highnoise = highnoise[::-1]
    return lownoise, highnoise


def _fit_points(periods, model):
    """
    Fit points to a noise model

    :param periods: periods in increasing order
    :type periods: list
    :param model: list of [period, value, slope]
    :type model: list of lists
    """
    x = np.log10(periods)
    xp = np.log10([x[0] for x in model])
    yp = [x[1] + x[2]*np.log10(x[0]) for x in model]
    assert np.all(np.diff(x) > 0), 'x is not increasing'
    assert np.all(np.diff(xp) > 0), 'xp is not increasing'
    return np.interp(x, xp, yp, left=np.nan, right=np.nan)


if __name__ == '__main__':
    print('not a command line code')
    sys.exit(1)
