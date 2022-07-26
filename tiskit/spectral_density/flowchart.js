op2=>operation: '\nSpectral Density Functions\n'
op4=>operation: import logging
op6=>operation: import xarray as xr
op8=>operation: import numpy as np
op10=>operation: from matplotlib import pyplot as plt
op12=>operation: from .utils import _prol1pi, _prol4pi, coherence_significance_level
op14=>operation: from obspy.core.stream import Stream
op16=>operation: from obspy.core import UTCDateTime
op18=>operation: from scipy import signal, stats
op20=>operation: from .Peterson_noise_model import Peterson_noise_model
sub22=>subroutine: np.seterr(all='ignore')
op24=>operation: class SpectralDensity():
    '\n    Class for spectral density functions.\n\n    The standard constructor is rarely used, generate objects using\n    `SpectralDensity.from_stream()`\n\n    No public attributes, access data through provided methods\n    '

    def __init__(self, chan_names, freqs, chan_units, n_windows, window_type, starttimes=None, data=None, responses=None):
        '\n        Args:\n            chan_names (list of str): channel names\n            freqs (np.ndarray): frequencies\n            chan_units (list of str): channel physical units (e.g m/s^2, Pa)\n            n_windows (int): of windows used to calculate spectra\n            windpw_type (str): type of window used\n            starttimes (list of UTCDateTime): starttime for each window\n            data (:class:`np.ndarray`):\n                one-sided spectral density functions.\n                shape = (len(chan_names), len(chan_names), len(freqs)\n                units = chan_units(i)*chan_units(j)/Hz\n            responses (:class:`np.ndarray`):\n                instrument response for each channel.\n                shape=(n_spects,n_freqs)\n                units=(counts/chan_units)\n        '
        n_ch = len(chan_names)
        n_f = len(freqs)
        shape = (n_ch, n_ch, n_f)
        dims = ('input', 'output', 'f')
        responses_shape = (n_ch, n_f)
        responses_dims = ('input', 'f')
        assert (freqs.size == n_f)
        assert (len(chan_units) == n_ch)
        if (starttimes is not None):
            for x in starttimes:
                assert isinstance(x, UTCDateTime)
        if (data is not None):
            assert (data.size() == shape)
            assert (data.dtype == 'complex')
        else:
            data = np.zeros(shape, dtype='complex')
        if (responses is not None):
            assert (responses.shape == responses_shape)
            assert (responses.dtype == 'complex')
        else:
            responses = np.ones((n_ch, n_f), dtype='complex')
        self._ds = xr.Dataset(data_vars={'spectra': (dims, np.zeros(shape, dtype='complex')), 'response': (responses_dims, responses)}, coords={'input': chan_names, 'output': chan_names, 'f': freqs, 'in_units': ('input', chan_units), 'out_units': ('output', chan_units)}, attrs={'n_windows': n_windows, 'window_type': window_type, 'starttimes': starttimes, 'long_name': 'spectral density function', 'units': 'input units * output units / Hz', 'description': 'One-sided spectral density functions'})

    def __str__(self):
        s = 'SpectralDensity object:\n'
        s += f'''	channels={self.channels}
'''
        s += '\tchannel_units={}\n'.format([self.channel_units(ch) for ch in self.channels])
        f = self.freqs
        s += f'''	{len(f)} frequencies, from {f[0]:.3g} to {f[(- 1)]:.3g}Hz
'''
        s += f'''	n_windows={self.n_windows}
'''
        s += f'	window_type={self.window_type}'
        return s

    def __eq__(self, other):
        return (self._ds == other._ds)

    @property
    def channels(self):
        '\n        Channel names\n\n        Returns:\n            (list of str):\n        '
        assert (list(self._ds.coords['input'].values) == list(self._ds.coords['output'].values))
        return list(self._ds.coords['input'].values)

    @property
    def freqs(self):
        '\n        Frequencies of the spectral density functions\n\n        Returns:\n            (:class:`numpy.ndarray`):\n        '
        return self._ds.coords['f'].values

    @property
    def window_type(self):
        '\n        The type of window used to calculate the spectral densities\n\n        Returns:\n            (str):\n        '
        return self._ds.window_type

    @property
    def starttimes(self):
        '\n        Start times for each data window used to calculate spectra\n\n        Returns:\n            (list of :class:`obspy.UTCDateTimes`):\n        '
        return self._ds.starttimes

    @property
    def n_windows(self):
        '\n        The number of data windows used to calculate spectra\n\n        Returns:\n            (int):\n        '
        return self._ds.n_windows

    @classmethod
    def from_stream(cls, stream, window_s=1000, windowtype='prol1pi', inv=None, data_cleaner=None, z_threshold=3):
        "\n        Calculate spectral density functions from the provided stream\n\n        Should add a window selection algorithm, for now just steps by\n        the window length\n\n        Args:\n            stream (:class:`obspy.core.stream.Stream`): data\n            window_s (float): desired window length in seconds\n            windowtype (str): window type, must be a valid\n            inv (:class:`obspy.core.inventory.Inventory`): inventory containing\n                instrument responses.  If none is found for the given channel,\n                will look in the channel's stats.response object\n            data_cleaner (:class:`DataCleaner`): Data cleaner to\n                apply to channels as ffts are calculated\n            subtract_tf_suffix (str): suffix to add to channel names if tf\n                is subtracted\n            z_threshold (float): reject windows with z-score greater than this\n                                 value\n        "
        if (not isinstance(stream, Stream)):
            raise ValueError(f'stream is a {type(stream)}, not obspy Stream')
        stream = _align_traces(stream)
        sr = stream[0].stats.sampling_rate
        if ((window_s * sr) > stream[0].stats.npts):
            raise ValueError('Requested window size > data length ({:g} > {:g} s)'.format(window_s, (stream[0].stats.npts / sr)))
        ws = int(_npow2((window_s * sr)))
        if (ws > stream[0].stats.npts):
            logging.warning(f'window pts > data pts ({ws:d} > {stream[0].stats.npts:d} pts), reducing ...')
            while (ws > stream[0].stats.npts):
                ws /= 2
            ws = int(ws)
            logging.warning(f'New window size={ws:d} pts')
        multfac = (2 / (ws * sr))
        (ft, evalresps, units) = ({}, {}, [])
        ids = [tr.id for tr in stream]
        if (not (len(ids) == len(set(ids)))):
            raise ValueError('stream has duplicate IDs')
        for id in ids:
            tr = stream.select(id=id)[0]
            (ft[id], f, sts) = _calculate_windowed_rfft(tr, ws, ws, windowtype)
            (ft[id], resp, evalresp, ft_units) = _correct_response(ft[id], f, id, tr.stats, inv)
            units.append(ft_units)
            evalresps[id] = evalresp
            if (resp is not None):
                tr.stats.response = resp
        (ft, sts) = cls._remove_outliers(ft, sts, z_threshold)
        n_winds = len(sts)
        if (data_cleaner is not None):
            dctfs = data_cleaner.DCTFs
            old_ids = ids
            ft = dctfs.ft_subtract_tfs(ft)
            ids = dctfs.update_channel_names(old_ids)
            evalresps = dctfs.update_channel_keys(evalresps)
        obj = cls(ids, f, units, n_winds, windowtype, starttimes=sts)
        for inp in ids:
            if (evalresps[inp] is not None):
                obj.put_channel_response(inp, evalresps[inp])
            for outp in ids:
                obj.put_crossspect(inp, outp, ((2 * np.mean((np.conj(ft[inp]) * ft[outp]), axis=0)) * multfac))
        return obj

    @staticmethod
    def _remove_outliers(ft, sts, z_threshold=3, recursive=True):
        'Remove  windows with z-score above z_threshold\n        \n        Args:\n            ft (dict): each value is an np.array where axis 0 = windows and\n                axis 1 = frequencies m windows x n frequencies\n            recursive (bool): repeat the test until there are no outliers\n        '
        ft = ft.copy()
        sts = sts.copy()
        n_reject = 1
        while (n_reject > 0):
            norm = np.array([np.mean(np.abs((ft[id] * np.conj(ft[id]))), axis=1) for id in ft.keys()])
            z_scores = stats.zscore(norm, axis=1)
            z_scores = np.max(np.abs(z_scores), axis=0)
            n_reject = np.count_nonzero((np.abs(z_scores) > z_threshold))
            if (n_reject > 0):
                logging.info('{:d} of {:d} had z_score > {:g}: rejected'.format(n_reject, len(z_scores), z_threshold))
                keepers = (z_scores <= z_threshold)
                sts = [x for (x, keep) in zip(sts, keepers.tolist()) if (keep is True)]
                for id in ft.keys():
                    ft[id] = ft[id][(keepers, :)]
            if (recursive is False):
                break
        return (ft, sts)

    def autospect(self, channel):
        '\n        Auto-spectral_density function for the given channel\n\n        Args:\n            channel (str): channel name\n        Returns:\n            (:class:`numpy.ndarray`): auto-spectral density function\n        '
        self._verify_channel(channel, 'in_channel')
        return np.abs(self._ds['spectra'].sel(input=channel, output=channel).values.flatten())

    def crossspect(self, in_channel, out_channel):
        '\n        Cross-spectral density function for the given channels\n\n        Args:\n            in_channel (str): input channel name\n            out_channel (str): output channel name\n        Returns:\n            (:class:`numpy.ndarray`): cross-spectral density function\n        '
        self._verify_channel(in_channel, 'in_channel')
        self._verify_channel(out_channel, 'out_channel')
        return self._ds['spectra'].sel(input=in_channel, output=out_channel).values.flatten()

    def _verify_channel(self, channel, ch_identifier):
        if (not isinstance(channel, str)):
            raise TypeError(f'{ch_identifier} is a {type(channel)}, not a str')
        if (channel not in self.channels):
            raise ValueError(f'{ch_identifier} "{channel}" not in channels {self.channels}')

    def put_autospect(self, channel, auto_spect):
        '\n        Equivalent to put_cross_spect(channel, channel, auto_spect)\n\n        Args:\n            channel (str): auto-spectra channel\n            auto_spect (:class:`numpy.ndarray`): the auto-spectral density\n        '
        if (not (auto_spect.shape == self.freqs.shape)):
            raise ValueError(f'auto_spect has different shape than freqs ({auto_spect.shape} vs {self.freqs.shape})')
        if (not (auto_spect.dtype == 'complex')):
            try:
                auto_spect = auto_spect.astype(dtype=complex)
            except Exception:
                raise ValueError('auto_spect could not be converted to dtype=complex')
        if (channel not in self.channels):
            raise ValueError('channel "{}" is not in channels {}'.format(channel, self.channels))
        self._ds['spectra'].loc[dict(input=channel, output=channel)] = auto_spect

    def replace_channel_name(self, channel, replacement):
        '\n        Args:\n            channel (str): original channel name\n            replacement (str): replacement channel name\n        '
        channel_names = self.channels
        channel_names[channel_names.index(channel)] = replacement
        self._ds['input'] = channel_names
        self._ds['output'] = channel_names

    def put_crossspect(self, in_channel, out_channel, cross_spect):
        '\n        Put data into one of the cross-spectra.  Also puts the complex\n        conjugate in the symmetric index\n\n        Args:\n            in_channel (str): cross-spectra input channel\n            out_channel (str): cross-spectra output channel\n            cross_spect (:class:`numpy.ndarray`): a cross-spectral density\n        '
        assert (cross_spect.shape == self.freqs.shape)
        assert (cross_spect.dtype == 'complex')
        assert (in_channel in self.channels)
        assert (out_channel in self.channels)
        self._ds['spectra'].loc[dict(input=in_channel, output=out_channel)] = cross_spect
        if (not (in_channel == out_channel)):
            self._ds['spectra'].loc[dict(input=out_channel, output=in_channel)] = np.conj(cross_spect)

    def channel_response(self, channel):
        "\n        A channel's instrument response\n\n        Args:\n            channel (str): channel name\n        Returns:\n            (:class:`numpy.ndarray`):\n        "
        return self._ds['response'].sel(input=channel)

    def put_channel_response(self, channel, response):
        "\n        Put a channel's instrument response into the object\n\n        Verifies that the response has the same shape as the object's\n        `frequency` property and that it is of type=`complex`\n\n        Args:\n            channel (str): the channel name\n            response (:class:`numpy.ndarray`): the response\n        "
        assert (response.shape == self.freqs.shape)
        assert (response.dtype == 'complex')
        assert (channel in self.channels)
        self._ds['response'].loc[dict(input=channel)] = response

    def channel_units(self, channel):
        '\n        Args:\n            channel (str): the channel name\n        Returns:\n            (str): Input (physical) units of the given channel\n        '
        return str(self._ds['spectra'].sel(input=channel).coords['in_units'].values)

    def units(self, in_channel, out_channel):
        '\n        The units of the given cross-  or auto-spectra\n\n        Args:\n            in_channel (str): input channel\n            out_channel (str): output channel\n        Returns:\n            (str): the units\n        '
        in_units = self.channel_units(in_channel)
        out_units = self.channel_units(out_channel)
        if (in_units == out_units):
            return f'({in_units})^2/Hz'
        return f'({in_units})*({out_units})/Hz'

    def coherence(self, in_chan, out_chan):
        '\n        The coherence for the given input and output channels\n\n        Args:\n            in_chan (str): input channel.  Must match one of the\n                coordinates in _ds\n            out_chan (str): output channel.  Must match one of the\n                coordinates in _ds\n        Returns:\n            (:class:`numpy.ndarray`): Coherence absolute value\n\n        Coherence is a real-valued quantity, for the cross-spectral phase,\n        use the cross-spectral density function.\n        From Bendat & Piersol (1986), Appendix B, gamma_xy^2 (f)\n        '
        if (in_chan not in self._ds.input):
            raise ValueError('"in_chan" not in spectral density matrix')
        if (out_chan not in self._ds.output):
            raise ValueError('"out_chan" not in spectral density matrix')
        coherence = ((np.abs(self.crossspect(in_chan, out_chan)) ** 2) / (self.autospect(in_chan) * self.autospect(out_chan)))
        return np.abs(coherence)

    def coh_signif(self, prob=0.95):
        '\n        The coherence significance level\n\n        Args:\n            prob (float): significance level (between 0 and 1)\n        Returns:\n            (float):\n        '
        return coherence_significance_level(self.n_windows, prob)

    def plot(self, **kwargs):
        'Shortcut for `plot_autospectra()`'
        self.plot_autospectra(**kwargs)

    def plot_autospectra(self, x=None, overlay=False, plot_peterson=True, show=True, outfile=None, title=None):
        '\n        Plot autospectra\n\n        Args:\n            x (list of str): limit to the listed channels\n            overlay (bool): put all spect on one axis\n            plot_peterson(bool): plot Peterson Noise model if any channel has\n                units of (m/s^2)^2/Hz\n            show (bool): show on desktop\n            outfile (str): save figure to this filename\n            title (str): custom plot title\n        Returns:\n            (:class:`numpy.ndarray`): array of axis pairs (amplitude, phase)\n        '
        x = self._get_validate_channel_names(x)
        if (not overlay):
            (rows, cols) = _squarish_grid(len(x))
        else:
            (rows, cols) = (1, 1)
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        (fig, axs) = plt.subplots(rows, cols, sharex=True)
        if (title is None):
            title = 'Auto-spectra'
        fig.suptitle(title)
        if (not overlay):
            for (key, i) in zip(x, range(len(x))):
                i_row = int((i / cols))
                i_col = (i - (cols * i_row))
                (axa, axp) = self.plot_one_spectra(key, key, fig, (rows, cols), (i_row, i_col), show_ylabel=(i_col == 0), show_xlabel=(i_row == (rows - 1)), show_phase=False, plot_peterson=plot_peterson)
                ax_array[(i_row, i_col)] = (axa, axp)
        else:
            (axa, axp) = (None, None)
            for (key, i) in zip(x, range(len(x))):
                (axa, axp) = self.plot_one_spectra(key, key, fig, (1, 1), (0, 0), show_ylabel=(i == (len(x) - 1)), show_xlabel=(i == (len(x) - 1)), ax_a=axa, ax_p=axp, show_phase=False, plot_peterson=plot_peterson)
            ax_array[(0, 0)] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_cross_spectra(self, x=None, show=True, show_coherence=False, outfile=None, plot_peterson=False):
        '\n        Plot cross (and auto) spectra\n\n        Args:\n            x (list of str): limit to the listed channels\n            show (bool): show on desktop\n            plot_peterson(bool): plot Peterson Noise model if any channel has\n                units of (m/s^2)^2/Hz\n            show_coherence (bool): show coherence as well\n        Returns:\n            :class:`numpy.ndarray`: array of axis pairs (amplitude, phase)\n        '
        x = self._get_validate_channel_names(x)
        n_subkeys = len(x)
        (rows, cols) = (n_subkeys, n_subkeys)
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        (fig, axs) = plt.subplots(rows, cols, sharex=True)
        fig.suptitle('Cross-spectra (dB ref UNITS/Hz)')
        for (in_chan, i) in zip(x, range(len(x))):
            for (out_chan, j) in zip(x, range(len(x))):
                title = (out_chan if (i == 0) else None)
                (axa, axp) = self.plot_one_spectra(in_chan, out_chan, fig, (rows, cols), (i, j), show_ylabel=(j == 0), show_xlabel=(i == (rows - 1)), ylabel=in_chan, label='units', title=title, show_coherence=show_coherence)
                ax_array[(i, j)] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_one_autospectra(self, key, **kwargs):
        '\n        Plot one autospectral density\n\n        Arguments are the same as for `plot_one_spectra()`, except\n        there is no `subkey` argument and `show_phase` is ignored\n        '
        kwargs['show_phase'] = False
        self.plot_one_spectra(key, key, **kwargs)

    def plot_one_spectra(self, key, subkey, fig=None, fig_grid=(1, 1), plot_spot=(0, 0), show_xlabel=True, show_ylabel=None, ax_a=None, ax_p=None, ylabel=None, label=None, title=None, show_coherence=False, show_phase=True, plot_peterson=True, outfile=None):
        "\n        Plot one spectral density\n\n        Args:\n            key (str): input (driving) channel\n            subkey (str): output (response) channel\n            fig (:class:`matplotlib.figure.Figure`): figure to plot on, if\n                None this method will plot on the current figure or create\n                a new figure.\n            fig_grid (tuple): this plot sits in a grid of this many\n                              (rows, columns)\n            subplot_spot (tuple): put this plot at this (row,column) of\n                                  the figure grid\n            show_xlabel (bool): put an xlabel on this subplot\n            show_ylabel (bool): put a ylabel on this subplot\n            ylabel (str): label to put on y axis (if show_label).  If not\n                speficied, will use 'dB ref UNITS/Hz'\n            label (str): 'units': print units without channel name in legend\n            ax_a (Axis): use an existing axis for the amplitude plot\n            ax_p (Axis): use this existing axis for the phase plot\n            title (str): title to put on this subplot\n            show_coherence (bool): draw coherence on the same plot\n            show_phase (bool): show phase as well as amplitude\n            plot_peterson(bool): plot Peterson Noise model if channel has\n                units of (m/s^2)^2/Hz\n            outfile (str): save figure to this filename\n\n        Returns:\n            (tuple): tuple containing\n                - :class:`matplotlib.axes.axis`: amplitude plot axis\n                - :class:`matplotlib.axes.axis`: phase plot axis\n        "
        psd = self.crossspect(key, subkey)
        in_units = self.channel_units(key)
        out_units = self.channel_units(subkey)
        if (in_units == out_units):
            PSD_units = f'({in_units})^2'
        else:
            PSD_units = f'{in_units}*{out_units}'
        f = self.freqs
        if (fig is None):
            fig = plt.gcf()
        if (ax_a is None):
            if show_phase:
                ax_a = plt.subplot2grid(((3 * fig_grid[0]), (1 * fig_grid[1])), (((3 * plot_spot[0]) + 0), (plot_spot[1] + 0)), rowspan=2)
            else:
                ax_a = plt.subplot2grid((fig_grid[0], fig_grid[1]), (plot_spot[0], plot_spot[1]))
        if show_coherence:
            ax2 = ax_a.twinx()
            ax2.semilogx(f, np.abs(self.coherence(key, subkey)), color='red', linewidth=0.5, alpha=0.8)
            ax2.axhline(self.coh_signif(0.95), color='red', linewidth=0.5, alpha=0.8, ls='--')
            ax2.set_ylim(0, 1)
            if (plot_spot[1] == (fig_grid[1] - 1)):
                ax2.set_ylabel('Coher', color='red')
            else:
                ax2.set_yticklabels([])
        psd[(psd == 0)] = None
        if (label is None):
            label = f'{subkey} ({PSD_units})'
        elif (label == 'units'):
            label = f'{PSD_units}'
        ax_a.semilogx(f, (10 * np.log10(np.abs(psd))), label=label)
        ax_a.set_xlim(f[1], f[(- 1)])
        if ((plot_peterson is True) and (PSD_units.lower() == '(m/s^2)^2')):
            (lownoise, highnoise) = Peterson_noise_model(f, True)
            ax_a.semilogx(f, lownoise, 'k--')
            ax_a.semilogx(f, highnoise, 'k--')
        if (label is not None):
            legend_1 = ax_a.legend()
            if show_coherence:
                legend_1.remove()
                ax2.add_artist(legend_1)
        if show_ylabel:
            if (ylabel is None):
                ylabel = 'dB ref UNITS/Hz'
            ax_a.set_ylabel(ylabel)
        if title:
            ax_a.set_title(title)
        if show_phase:
            if (ax_p is None):
                ax_p = plt.subplot2grid(((3 * fig_grid[0]), (1 * fig_grid[1])), (((3 * plot_spot[0]) + 2), (plot_spot[1] + 0)))
            ax_p.semilogx(f, np.degrees(np.angle(psd)))
            ax_p.set_ylim((- 180), 180)
            ax_p.set_xlim(f[1], f[(- 1)])
            ax_p.set_yticks(((- 180), 0, 180))
            if show_ylabel:
                pass
            else:
                ax_p.set_yticklabels([])
            bottom_axis = ax_p
        else:
            ax_p = None
            bottom_axis = ax_a
        if show_xlabel:
            bottom_axis.set_xlabel('Frequency (Hz)')
        else:
            bottom_axis.set_xticklabels([])
        if outfile:
            plt.savefig(outfile)
        return (ax_a, ax_p)

    def plot_coherences(self, x=None, y=None, overlay=False, show=True, outfile=None):
        '\n        Plot coherences\n\n        Args:\n            x (list of str): limit to the listed input channels\n            y (list of str): limit to the listed output channels\n            overlay (bool): put all coherences on one plot\n            show (bool): show on desktop\n            outfile (str): save to the named file\n\n        Returns:\n            (:class:`numpy.ndarray`): array of axis pairs (amplitude, phase)\n        '
        x = self._get_validate_channel_names(x)
        y = self._get_validate_channel_names(y)
        if (not overlay):
            (rows, cols) = (len(x), len(y))
            ax_array = np.ndarray((rows, cols), dtype=tuple)
            (fig, axs) = plt.subplots(rows, cols, sharex=True)
            fig.suptitle('Coherences')
            for (in_chan, i) in zip(x, range(len(x))):
                for (out_chan, j) in zip(y, range(len(y))):
                    title = (out_chan if (i == 0) else None)
                    (axa, axp) = self.plot_one_coherence(in_chan, out_chan, fig, (rows, cols), (i, j), show_ylabel=(j == 0), show_xlabel=(i == (rows - 1)), ylabel=in_chan, title=title)
                    ax_array[(i, j)] = (axa, axp)
        else:
            ax_array = np.ndarray((1, 1), dtype=tuple)
            (fig, axs) = plt.subplots(1, 1, sharex=True)
            fig.suptitle('Coherences')
            labels = []
            (axa, axp) = (None, None)
            for (in_chan, i) in zip(x, range(len(x))):
                for (out_chan, j) in zip(y, range(len(y))):
                    if (out_chan == in_chan):
                        break
                    if (f'{out_chan}-{in_chan}' in labels):
                        break
                    label = f'{in_chan}-{out_chan}'
                    (axa, axp) = self.plot_one_coherence(in_chan, out_chan, fig, (1, 1), (0, 0), ylabel='Coherence', label=label, ax_a=axa, ax_p=axp)
                    labels.append(label)
            ax_array[(0, 0)] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_one_coherence(self, in_chan, out_chan, fig=None, fig_grid=(1, 1), plot_spot=(0, 0), show_xlabel=True, show_ylabel=True, ax_a=None, ax_p=None, ylabel=None, label=None, title=None, show_phase=True):
        "\n        Plot one coherence\n\n        Args:\n            in_chan (str): input (driving) channel\n            out_chan (str): output (response) channel\n            fig (:class:`matplotlib.figure.Figure`): figure to plot on, if\n                None this method will plot on the current figure or create\n                a new figure.\n            fig_grid (tuple): this plot sits in a grid of this many\n                              (rows, columns)\n            plot_spot (tuple): put this plot at this (row,column) of\n                                  the figure grid\n            show_xlabel (bool): put an xlabel on this subplot\n            show_ylabel (bool): put a ylabel on this subplot\n            ylabel (str): label to put on y axis (if show_label).  If not\n                speficied, will use 'dB ref UNITS/Hz'\n            label (str): text to put in legend\n            ax_a (Axis): use an existing axis for the amplitude plot\n            ax_p (Axis): use this existing axis for the phase plot\n            title (str): title to put on this subplot\n            show_phase (bool): show phase as well as amplitude\n\n        Returns:\n            (tuple): tuple containing:\n                - (:class:`matplotlib.axes.axis`): amplitude plot axis\n                - (:class:`matplotlib.axes.axis`): phase plot axis\n        "
        ds = self._ds['spectra'].sel(input=in_chan, output=out_chan)
        f = self._ds.coords['f'].values
        if (fig is None):
            fig = plt.gcf()
        if (ax_a is None):
            if show_phase:
                ax_a = plt.subplot2grid(((3 * fig_grid[0]), (1 * fig_grid[1])), (((3 * plot_spot[0]) + 0), (plot_spot[1] + 0)), rowspan=2)
            else:
                ax_a = plt.subplot2grid((fig_grid[0], fig_grid[1]), (plot_spot[0], plot_spot[1]))
        ax_a.semilogx(f, np.abs(self.coherence(in_chan, out_chan)), label=label)
        ax_a.axhline(self.coh_signif(0.95), color='red', linewidth=0.5, alpha=0.8, ls='--')
        ax_a.set_ylim(0, 1)
        if (label is not None):
            ax_a.legend()
        if show_ylabel:
            if (ylabel is None):
                ylabel = 'Coherence'
            ax_a.set_ylabel(ylabel)
        if title:
            ax_a.set_title(title)
        if show_phase:
            if (ax_p is None):
                ax_p = plt.subplot2grid(((3 * fig_grid[0]), (1 * fig_grid[1])), (((3 * plot_spot[0]) + 2), (plot_spot[1] + 0)))
            ax_p.semilogx(f, np.degrees(np.angle(ds)))
            ax_p.set_ylim((- 180), 180)
            ax_p.set_xlim(f[1], f[(- 1)])
            ax_p.set_yticks(((- 180), 0, 180))
            if show_ylabel:
                pass
            else:
                ax_p.set_yticklabels([])
            bottom_axis = ax_p
        else:
            ax_p = None
            bottom_axis = ax_a
        if show_xlabel:
            bottom_axis.set_xlabel('Frequency (Hz)')
        else:
            bottom_axis.set_xticklabels([])
        return (ax_a, ax_p)

    def _get_validate_channel_names(self, x):
        '\n        If x is None, return list of all channel names\n        If x is a list, validate all of the names\n        '
        if (x is None):
            return list(self._ds.coords['input'].values)
        for key in x:
            if (key not in list(self._ds.coords['input'].values)):
                ValueError('key "{key}" not in channel list')
        return x

    @staticmethod
    def _remove_subtracted_loc(id):
        '\n        Remove loc code characters including and after first "-"\n\n        Allows the use of "-?" in the loc code to specify removed coherent\n        noise\n\n        Args:\n            id (str): seed ID code\n\n        Example:\n            >>> SD._remove_subtracted_loc(\'hello\')\n            \'hello\'\n            >>> SD._remove_subtracted_loc(\'NN.SSSS.LL.CCC\')\n            \'NN.SSSS.LL.CCC\'\n            >>> SD._remove_subtracted_loc(\'NN.SSSS.-LL.CCC\')\n            \'NN.SSSS..CCC\'\n            >>> SD._remove_subtracted_loc(\'NN.SSSS.L-LL.CCC\')\n            \'NN.SSSS.L.CCC\'\n        '
        comps = id.split('.')
        if (not (len(comps) == 4)):
            return id
        if ('-' not in comps[2]):
            return id
        comps[2] = comps[2].partition('-')[0]
        return '.'.join(comps)

    @staticmethod
    def _sliding_window(a, ws, ss=None, win_taper='hanning'):
        "\n        Split a data array into overlapping, tapered sub-windows\n\n        Args:\n            a (:class:`numpy.ndarray`): 1D array of data to split\n            ws (int): Window size in samples\n            ss (int): Step size in samples. If not provided, window and step\n                size are equal.\n            win_taper (str): taper to apply to data ['hanning', 'prol4pi',\n                'prol1pi', 'bartlett', 'blackman', 'hamming']\n\n        Returns:\n            (tuple): tuple conaining\n                out (:class:`numpy.ndarray`): 1D array of windowed data\n                nd (int): Number of windows\n                offsets (list): list of sample offsets for window starts\n        "
        if (ss is None):
            ss = ws
        ws = int(ws)
        if (ws > len(a)):
            raise ValueError('window size > data length ({} > {})'.format(ws, len(a)))
        nd = (1 + ((len(a) - ws) // ss))
        out = np.ndarray((nd, ws), dtype=a.dtype)
        if (win_taper in ['hanning', 'hamming', 'blackman', 'bartlett']):
            taper = eval(f'np.{win_taper}(ws)')
        elif (win_taper == 'prol1pi'):
            taper = _prol1pi(ws)
        elif (win_taper == 'prol4pi'):
            taper = _prol4pi(ws)
        else:
            raise ValueError(f'Unknown taper type "{win_taper}"')
        offsets = []
        if (nd == 0):
            out = (signal.detrend(a) * taper)
            offsets.append(0)
        for i in range(nd):
            start = (i * ss)
            stop = (start + ws)
            out[i] = (signal.detrend(a[start:stop]) * taper)
            offsets.append(start)
        return (out, nd, offsets)
st27=>start: start _align_traces
io29=>inputoutput: input: stream
op32=>operation: 'Trim stream so that all traces are aligned and same length'
op34=>operation: first_start = last_start = stream[0].stats.starttime
op36=>operation: first_end = last_end = stream[0].stats.endtime
op38=>operation: sampling_rate = stream[0].stats.sampling_rate
cond41=>condition: for tr in stream[1:]
cond100=>condition: if (tr.stats.starttime > last_start)
op104=>operation: last_start = tr.stats.starttime
cond121=>condition: if (tr.stats.endtime < first_end)
op125=>operation: first_end = tr.stats.endtime
cond142=>operation: raise ValueError('not all traces have same sample rate') if  (not (tr.stats.sampling_rate == sampling_rate))
cond130=>operation: last_end = tr.stats.endtime if  (tr.stats.endtime > last_end)
cond109=>operation: first_start = tr.stats.starttime if  (tr.stats.starttime < first_start)
cond155=>condition: if np.any([np.ma.count_masked(x.data) for x in stream])
sub159=>subroutine: logging.warning('Unmasking masked data (usually a gap or overlap)')
op161=>operation: stream = stream.split().merge(fill_value='interpolate')
cond167=>operation: raise ValueError('There are non-overlapping traces') if  (last_start >= first_end)
cond178=>operation: logging.info('Cutting up to {}s from trace starts'.format((last_start - first_start))) if  ((last_start - first_start) > (1 / sampling_rate))
cond189=>operation: logging.info('Cutting up to {last_start-first_start}s from trace ends') if  ((last_end - first_end) > (1 / sampling_rate))
sub199=>subroutine: stream.trim(last_start, first_end)
op201=>operation: min_len = min([tr.stats.npts for tr in stream])
op203=>operation: max_len = max([tr.stats.npts for tr in stream])
cond206=>condition: if (not (max_len == min_len))
cond211=>operation: tr.data = tr.data[:min_len] while  tr in stream
io229=>inputoutput: output:  stream
e227=>end: end function return

op2->op4
op4->op6
op6->op8
op8->op10
op10->op12
op12->op14
op14->op16
op16->op18
op18->op20
op20->sub22
sub22->op24
op24->st27
st27->io29
io29->op32
op32->op34
op34->op36
op36->op38
op38->cond41
cond41(yes)->cond100
cond100(yes)->op104
op104->cond121
cond121(yes)->op125
op125->cond142
cond142->cond41
cond121(no)->cond130
cond130->cond142
cond100(no)->cond109
cond109->cond121
cond41(no)->cond155
cond155(yes)->sub159
sub159->op161
op161->cond167
cond167->cond178
cond178->cond189
cond189->sub199
sub199->op201
op201->op203
op203->cond206
cond206(yes)->cond211
cond211->io229
io229->e227
cond206(no)->io229
cond155(no)->cond167

