"""
Calculate spectra and coherences for a given station/time period
"""
from tiskit import SpectralDensity  # , PetersonNoiseModel
from matplotlib import pyplot as plt
import numpy as np
import scipy.signal as ssig
from obspy.signal.invsim import cosine_taper

from read_data_inventory import read_data_inventory


def main():
    stream, inv = read_data_inventory()
    for windowtype in ['prol1pi', 'prol4pi']:
        BHZstream = stream.select(channel='BHZ')
        spect = SpectralDensity.from_stream(BHZstream, inv=inv,
                                            windowtype=windowtype)
        f_sp, s_sp = scipy_spect(BHZstream[0], windowtype)
        #  lownoise, highnoise = PetersonNoiseModel(spect.freqs, as_freqs=True)
        plt.semilogx(spect.freqs,
                     10*np.log10(np.abs(spect.autospect('XE.CC06..BHZ'))),
                     label='SpectralDensity')
        plt.semilogx(f_sp, 10*np.log10(np.abs(s_sp)),
                     label='scipy.signal.welch)')
        # plt.plot(spect.freqs, highnoise, 'g--')
        # plt.plot(spect.freqs, lownoise, 'r--')
        plt.ylabel('Autospectral density')
        plt.ylabel('Frequencies (Hz)')
        plt.title(f'{windowtype} window')
        plt.legend()
        plt.show()


def scipy_spect(tr, windowtype, window_length=1000.):
    sampling_rate = tr.stats.sampling_rate
    nfft = int(2**(np.ceil(np.log2(window_length * sampling_rate))))
    nlap = int(0.75 * nfft)

    window = choose_window_scipy(windowtype, nfft)
    freq, spec = ssig.welch(tr.data, sampling_rate, nperseg=nfft,
                            window=window, detrend="linear",
                            noverlap=nlap)
    # Requires response in trace (works because done after SpectralDensity)
    resp = tr.stats.response
    evalresp, freqs = resp.get_evalresp_response(
        t_samp=1 / sampling_rate, nfft=nfft, output="ACC")
    assert np.all(freq == freqs)
    # Get the amplitude response (squared)
    spec /= np.absolute(evalresp * np.conjugate(evalresp))
    # leave out first entry (zero-freq)
    return freq[1:], spec[1:]


def choose_window_scipy(wind_type, n):
    if wind_type == 'fft_taper':
        return _fft_taper(n)
    elif wind_type == 'prol1pi':
        return _prol1pi(n)
    elif wind_type == 'prol4pi':
        return _prol4pi(n)
    else:
        raise NameError('wind_type must be fft_taper, prol1pi or prol4pi')


def _fft_taper(n):
    """
    Cosine taper, 10 percent at each end (like done by [McNamara2004]_).
    """
    return cosine_taper(n, 0.2)


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


if __name__ == "__main__":
    main()
