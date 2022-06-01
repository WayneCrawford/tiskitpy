"""
Calculate spectral transfer functions
"""
import matplotlib.pyplot as plt
import numpy as np
import warnings

# Set variables
# spect_library = 'scipy'  # 'mlab' or 'scipy': mlab gives weird coherences!


class TransferFunction:
    def __init__(self, freqs, data, uncerts, drive_stats, resp_stats,
                 gooddata=None, noisechan='response'):
        """
        :parm freqs: 1-D array of frequencies
        :parm data: 1-D array of transfer function
        :parm uncerts: 1-D array of uncertainties
        :parm units: data units (str)
        :parm noisechan: 'response', 'driving', 'equal' or 'unknown'
        :type drive_stats, resp_stats: :class:`~obspy.core.trace.Stats`
        """
        self.freqs = freqs
        self.data = data
        self.uncerts = uncerts
        self.drive_stats = drive_stats
        self.resp_stats = resp_stats
        self.gooddata = gooddata
        self.noisechan = noisechan

    @classmethod
    def calc(cls, spects, cohers, drivechan='H', respchan='Z',
             noisechan='response', verbose=True):
        """
        Calculate the transfer function between 2 channels

        XF = calcXF(spect,cohers,drivechan,respchan,noisechan)
        input:
            spects    (list)  contains PSDs calculated using calc_PSDs()
            cohers    (list)  contains coherences
            drivechan (str)   the driving channel name (or last character(s))
            respchan  (str)   the response channel name (or last character(s))
            noisechan (str)	  is which channel to assume contains the noise
                'response'	[default] Assume all noise on the response channel
                'driving'	Assume all noise on the driving channel
                'equal'		Assume the same signal/noise on the driving  and
                            response channels
                'unknown'	Make no assumption about noise
        output:
            dictionary containing the transfer function
        """

        if verbose:
            print(f'Calculating XF between "{drivechan}" (driving) and'
                  f'"{respchan}" (response) channels, assume noise is on'
                  f'{noisechan}')

        # FIND PSD and coherence channels matching drivechan and respchan
        coh = None
        for c in cohers:
            chan_i = cohers.stats[c.ch_nums[0]].channel
            chan_j = cohers.stats[c.ch_nums[1]].channel
            if (chan_i.endswith(drivechan) and chan_j.endswith(respchan)) or \
               (chan_i.endswith(respchan) and chan_j.endswith(drivechan)):
                coh = c
                break
        if not coh:
            warnings.warn('Did no t find a coherence with channels {}, {}'.
                          format(drivechan, respchan))
            return False
        drive_spect = None
        resp_spect = None
        for s in spects:
            if s['chan'].endswith(drivechan):
                drive_spect = s
                if verbose:
                    print('drive_spect channel is "{}"'.format(s['chan']))
            elif s['chan'].endswith(respchan):
                resp_spect = s
                if verbose:
                    print('resp_spect channel is "{}"'.format(s['chan']))
        if not drive_spect:
            warnings.warn('Did not find a spectra with channel {}'.
                          format(drivechan))
            return False
        if not resp_spect:
            warnings.warn('Did not find a spectra with channel {}'.
                          format(respchan))
            return False

        # CALCULATE RESP/DRIVE
        RespOverDrive = np.sqrt(np.divide(resp_spect['data'],
                                          drive_spect['data']))

        # CALCULATE TRANSFER FUNCTION
        # Equations from Bendat&Piersol "Random Data" 1986, pp 176-181 (xfs)
        # & pp 317, Table 9.6
        cohmagsq = np.multiply(np.absolute(coh['data']),
                               np.absolute(coh['data']))
        errbase = np.divide(np.sqrt(np.ones(cohmagsq.shape) - cohmagsq),
                            2 * coh['num_windows'] * cohmagsq)
        if noisechan == 'response':
            xf = np.multiply(RespOverDrive, coh['data'])
            xferr = np.multiply(np.abs(xf), errbase)
        elif noisechan == 'driving':
            xf = np.divide(RespOverDrive, coh['data'])
            xferr = np.multiply(np.abs(xf), errbase)
        elif noisechan == 'equal':
            xf = RespOverDrive
            xferr = np.abs(np.multiply(xf, errbase))
        elif noisechan == 'unknown':
            xf = RespOverDrive
            # Ad-hoc error guesstimate
            maxerr = np.abs(np.power(coh['data'], -1.)) + errbase
            minerr = np.abs(coh['data']) - errbase
            xferr = np.abs(np.multiply(xf, (maxerr-minerr) / 2))
        else:
            warnings.warn('Invalid: noisechan = ' + noisechan)
            return False
        return cls(drive_spect.freqs, xf, xferr,
                   drive_spect.stats, resp_spect.stats,
                   np.absolute(coh['data']) > coh['signif_level'],
                   noisechan)
        # XF=dict(drive_info=dict(channel=drive_spect['chan'],
        #                         units=drive_spect['units']),
        #         resp_info= dict(channel=resp_spect['chan'],
        #                         units=resp_spect['units']),
        #         freq=     drive_spect['freq'],
        #         data=      xf,
        #         error=     xferr,
        #         noisechan=noisechan,
        #         gooddata=  np.absolute(coh['data']) > coh['signif_level'])

    def clean(self, sdf):
        """
        Clean a SpectralDensity object using the given transfer function
        """
    
    def plot(self, outfile=None, debug=False):
        """
        plot transfer function
        """
        plt.figure(1)
        plt.clf()
        if debug:
            print(self.freqs[0])
            print(self.data[0])
            print(self.uncert[0])
        # plt.errorbar(np.log10(XF['freq']), np.absolute(XF['data']),
        #              yerr=np.absolute(XF['error']))
        # plt.ylim(0 , np.max(np.absolute(XF['data'])))
        plt.loglog(self.freqs, np.absolute(self.data), marker='o',
                   linestyle='')
        plt.loglog(self.freqs, np.absolute(self.data) + self.uncert)
        plt.loglog(self.freqs, np.absolute(self.data) - self.uncert)
        plt.ylim(np.min(np.absolute(self.data)),
                 np.max(np.absolute(self.data)))
        plt.title(f'Transfer function, noise channel: ({self.noisechan})')
        plt.ylabel(self.resp_stats.channel + '/' + self.drive_stats.channel)
        plt.xlabel('log(freq)')
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
        return


if __name__ == "__main__":
    import doctest
    doctest.testmod()
