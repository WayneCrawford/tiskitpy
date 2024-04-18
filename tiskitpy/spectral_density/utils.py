"""
Functions to calculate spectra, coherences and transfer functions
"""
from obspy.signal.invsim import cosine_taper
import numpy as np

from ..logger import init_logger

logger = init_logger()

# Set variables
# spect_library = 'scipy'  # 'mlab' or 'scipy': mlab gives weird coherences!


def choose_window_scipy(wind_type, n):
    if wind_type == 'fft_taper':
        return _fft_taper(n)
    elif wind_type == 'prol1pi':
        return _prol1pi(n)
    elif wind_type == 'prol4pi':
        return _prol4pi(n)
    else:
        raise NameError('wind_type must be fft_taper, prol1pi or prol4pi')


def choose_window_mlab(wind_type):
    if wind_type == 'fft_taper':
        return _fft_taper_times
    elif wind_type == 'prol1pi':
        return _prol1pi_times
    elif wind_type == 'prol4pi':
        return _prol4pi_times
    else:
        raise NameError('wind_type must be fft_taper, prol1pi or prol4pi')


def seed_code(stats):
    """
    Returns SEED code from obspy stats container

    >>> from obspy.core.trace import Stats
    >>> stats = Stats()
    >>> stats.network = 'YV'
    >>> stats.station = 'TEST'
    >>> stats.channel = 'BHZ'
    >>> _seed_code(stats)
    'YV.TEST..BHZ'
    """
    try:
        return '{}.{}.{}.{}'.format(stats.network, stats.station,
                                    stats.location, stats.channel)
    except AttributeError:
        return ''


def _fft_taper(n):
    """
    Cosine taper, 10 percent at each end (like done by [McNamara2004]_).
    """
    return cosine_taper(n, 0.2)


def _fft_taper_times(data):
    return data * _fft_taper(len(data))


def _prol1pi(n):
    """
    1-pi prolate spheroidal window
    """
    fOnE = np.ones(n)
    x = (2*np.arange(1, n+1) - fOnE) / n
    u = (2*fOnE - x) * x
    w = np.sqrt(2.) * ((((((((((
        5.3476939016920851e-11 * u + 2.2654256220146656e-9*fOnE) * u +
        7.8075102004229667e-8*fOnE) * u + 2.1373409644281953e-6*fOnE) * u +
        4.5094847544714943e-5*fOnE) * u + 7.0498957221483167e-4*fOnE) * u +
        7.7412693304064753e-3*fOnE) * u + 5.5280627452077586e-2*fOnE) * u +
        2.2753754228751827e-1*fOnE) * u + 4.3433904277546202e-1*fOnE) * u +
        2.2902051859068017e-1*fOnE)
    return w


def _prol1pi_times(data):
    return data * _prol1pi(len(data))


def _prol4pi(n):
    """
    4-pi prolate spheroidal window.
    """
    fOnE = np.ones(n)
    x = (2*np.arange(1, n+1) - fOnE) / n
    u = (2*fOnE - x) * x
    w = np.sqrt(2./.508125548147497) * (((((((((((((((((((((
        2.6197747176990866e-11 * u + 2.9812025862125737e-10*fOnE) * u +
        3.0793023552299688e-9*fOnE) * u + 2.8727486379692354e-8*fOnE) * u +
        2.4073904863499725e-7*fOnE) * u + 1.8011359410323110e-6*fOnE) * u +
        1.1948784162527709e-5*fOnE) * u + 6.9746276641509466e-5*fOnE) * u +
        3.5507361197109845e-4*fOnE) * u + 1.5607376779150113e-3*fOnE) * u +
        5.8542015072142441e-3*fOnE) * u + 1.8482388295519675e-2*fOnE) * u +
        4.8315671140720506e-2*fOnE) * u + 1.0252816895203814e-1*fOnE) * u +
        1.7233583271499150e-1*fOnE) * u + 2.2242525852102708e-1*fOnE) * u +
        2.1163435697968192e-1*fOnE) * u + 1.4041394473085307e-1*fOnE) * u +
        5.9923940532892353e-2*fOnE) * u + 1.4476509897632850e-2*fOnE) * u +
        1.5672417352380246e-3*fOnE) * u + 4.2904633140034110e-5*fOnE)
    return w


def _prol4pi_times(data):
    return data * _prol4pi(len(data))


def coherence_significance_level(n_windows, prob=0.95):
    """
    Definition: L_1(alpha, q) = sqrt(1-alpha**(1/q))

    where alpha = 1-prob and 2(q+1) = nwinds (degree of freedom)

    For nwinds >> 1, L1 ~ sqrt(1-alpha**(2/nwinds))
    For a 95% signif level this comes out to
        sqrt(1-.05**(2/nwinds)) for nwinds >> 1.
    I previously used sqrt(2/nwinds) for the 95% signif level (alpha=0.05),
    but L1 is much closer to sqrt(6/nwinds).

    Args:
        n_windows (int): number of windows
        prob (float): significance level (between 0 and 1)
    """
    assert prob >= 0 and prob <= 1
    alpha = 1 - prob
    q = n_windows/2 - 1
    return np.sqrt(1 - alpha ** (1. / q))


if __name__ == "__main__":
    import doctest
    doctest.testmod()
