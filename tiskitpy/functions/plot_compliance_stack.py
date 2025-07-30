"""
Python compliance class

Authors:  W. Crawford, A. Doran
"""
import numpy as np
from matplotlib import pyplot as plt

from ..response_functions import ResponseFunctions


def plot_compliance_stack(psd, zstr, pstr, water_depth, seawater_density=1030,
                          show=True, outfile=None):
    """
    Plot from top to bottom: Z PSD, P PSD, coherence, Z/P

    Args:
        psd (SpectralDensity): PSDs including Z and P
        zstr (str): channel id sub/string matching the Z channel (see
            :meth:`SpectralDensity.channel_id() documentation)
        pstr (str): channel id sub/string matching the P channel
        water_depth (float): water depth in meters
        seawater_density (float): average water density overhead (kg/m^3)
        show (bool): show the result on the screen
        outfile (str): save the plot to the named file
    """
    # Validate inputs
    try:
        _ = psd.channel_id(zstr)
    except Exception:
        raise ValueError(f'{zstr=} not a valid/unique channel id for {psd=}')
    try:
        _ = psd.channel_id(pstr)
    except Exception:
        raise ValueError(f'{pstr=} not a valid/unique channel id for {psd=}')
    assert water_depth > 0, f'{water_depth=} is not greater than 0'
    for id, units in zip((zstr, pstr), ('m/s^2', 'Pa')):
        assert psd.channel_units(id) == units, f'{id} units are {psd.channel_units(id)},  not "{units}"'

    fig, axs = plt.subplots(4, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    # Plot Z PSD
    axs[0].semilogx(psd.freqs, 20*np.log10(np.abs(psd.autospect(zstr))))
    axs[0].set_ylabel(r'Z (dB ref 1 $m/s^2/\sqrt{Hz}$)', fontsize='small')
    axs[0].set_title('Compliance Stack')
    # Plot Pressure PSD
    axs[1].semilogx(psd.freqs, 10*np.log10(np.abs(psd.autospect(pstr))))
    axs[1].set_ylabel(r'P (dB ref 1 $Pa/\sqrt{Hz}$)', fontsize='small')
    # Plot Pressure-Z coherence
    psd.plot_one_coherence(pstr, zstr, fig=fig, ax_a=axs[2], show_phase=False,
                           ylabel='')
    axs[2].set_ylabel('Coherence', fontsize='small')
    axs[2].set_ylim([0.001, 1])
    # Plot Z/P ratios
    for noise_channel, label, color in zip(('output', 'input'),
                                           ('z_noise', 'p_noise'),
                                           ('r', 'b')):
        frf = ResponseFunctions(psd, pstr, [zstr], noise_channel=noise_channel)
        igood = np.abs(frf.uncertainty(zstr)) < np.abs(frf.value(zstr))
        axs[3].errorbar(frf.freqs[igood],
                        y=np.abs(frf.value(zstr))[igood],
                        yerr=np.abs(frf.uncertainty(zstr))[igood],
                        fmt='.', ms=1, c=color, label=label)
    # Overlay theoretical Z/P relation for LF Rayleigh waves (seafloor moving
    # water column): m/s^2/Pa = 1/rho*H
    axs[3].axhline(1/(seawater_density*water_depth), c='k', ls='--')
    axs[3].text(frf.freqs[0], 1/(seawater_density*water_depth), r'1/$\rho H$',
                verticalalignment='bottom')
    axs[3].set_ylabel(f'Z/P ({frf.output_units(zstr)}/{frf.input_units})',
                      fontsize='small')
    axs[3].set_xlabel('Frequency (Hz)')
    axs[3].set_yscale('log')
    # Put predicted compliance max frequency vertical line on each plot
    IG_fmax = np.sqrt(9.8/(2*np.pi*water_depth))  # about one wavelength
    for ax in axs:
        ax.axvline(IG_fmax, c='k', ls='--')
    axs[0].text(IG_fmax, np.max(20*np.log10(np.abs(psd.autospect(zstr)))),
                r'$\sqrt{\frac{g}{2 \pi H}}$', rotation='vertical',
                horizontalalignment='right', verticalalignment='top')

    # Show or save plot
    if show is False and outfile is None:
        raise ValueError('Plot neither shown nor saved!')
    if show is True:
        plt.show()
    if outfile is not None:
        plt.savefig(outfile)


if __name__ == '__main__':
    print('not a command line code')
    sys.exit(1)
