"""
Spectral Density Functions
"""
import pickle
# from dataclasses import dataclass

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
# from obstools.atacr import DayNoise, StaNoise  # removed for readthedocs
from .utils import _prol1pi, _prol4pi, coherence_significance_level
from obspy.core.stream import Stream
from obspy.core import UTCDateTime
from scipy import signal

from .Peterson_noise_model import Peterson_noise_model

np.seterr(all='ignore')
# np.set_printoptions(threshold=sys.maxsize)


class SpectralDensity:
    """
    Class for spectral density functions.
    
    The standard constructor is rarely used, generate objects using
    `SpectralDensity.from_stream()`
    
    No public attributes, access data through provided methods
    """
    def __init__(self, chan_names, freqs, chan_units, n_windows, window_type,
		 starttimes=None, data=None, responses=None):
		"""
		    Args:
            chan_names (list of str): channel names
            freqs (np.ndarray): frequencies
            chan_units (list of str): channel physical units (e.g m/s^2, Pa)
            n_windows (int): of windows used to calculate spectra
            windpw_type (str): type of window used
            starttimes (list of UTCDateTime): starttime for each window
            data (:class:`np.ndarray`):
                one-sided spectral density functions.
                shape = (len(chan_names), len(chan_names), len(freqs)
                units = chan_units(i)*chan_units(j)/Hz
            responses (:class:`np.ndarray`):
                instrument response for each channel.
                shape=(n_spects,n_freqs)
                units=(counts/chan_units) 
		"""
        # Validate Dimensions
        n_ch = len(chan_names)
        n_f = len(freqs)
        shape = (n_ch, n_ch, n_f)
        dims = ('input', 'output', 'f')
        responses_shape = (n_ch, n_f)
        responses_dims = ('input', 'f')
        assert freqs.size == (n_f)   # Make sure it's one dimensional
        assert len(chan_units) == n_ch
        if starttimes is not None:
            for x in starttimes:
                assert isinstance(x, UTCDateTime)
        if data is not None:
            assert data.size() == shape
            assert data.dtype == 'complex'
        else:
            data = np.zeros(shape, dtype='complex')
        if responses is not None:
            assert responses.shape == responses_shape
            assert response.dtype == 'complex'
        else:
            responses = np.ones((n_ch, n_f), dtype='complex')

        self._ds = xr.Dataset(
            data_vars={"spectra": (dims,  np.zeros(shape, dtype='complex')),
                       "response": (responses_dims, responses)},
            # dims=dims,
            coords={"input": chan_names,
                    "output": chan_names,
                    "f": freqs,
                    "in_units": ("input", chan_units),
                    "out_units": ("output", chan_units)},
            attrs={"n_windows": n_windows,
                   "window_type": window_type,
                   "starttimes": starttimes,
                   "long_name": "spectral density function",
                   "units": "input units * output units / Hz",
                   "description": "One-sided spectral density functions"}
            )

    def __str__(self):
        s = 'SpectralDensity object:\n'
        s += f'\tchannels={self.channels}\n'
        s += f'\tchannel_units={[self.channel_units(ch) for ch in self.channels]}\n'
        f = self.freqs
        s += f'\t{len(f)} frequencies, from {f[0]:.3g} to {f[-1]:.3g}Hz\n'
        s += f'\tn_windows={self.n_windows}\n'
        s += f'\twindow_type={self.window_type}'
        return s

    @property
    def channels(self):
        """
	Channel names
	
        Returns:
            (list of str):
        """
        assert (list(self._ds.coords['input'].values)
                == list(self._ds.coords['output'].values))
        return list(self._ds.coords['input'].values)

    @property
    def freqs(self):
        """
        Frequencies of the spectral density functions

        Returns:
            (:class:`numpy.ndarray`): 
        """
        return self._ds.coords['f'].values

    def autospect(self, channel):
        """
        Auto-spectral_density function for the given channel

        Args:
            channel (str): channel name
        Returns:
            (:class:`numpy.ndarray`): auto-spectral density function
        """
        self._verify_channel(channel, "in_channel")
        return np.abs(self._ds["spectra"].sel(
            input=channel, output=channel).values.flatten())

    def crossspect(self, in_channel, out_channel):
        """
        Cross-spectral density function for the given channels

        Args:
            in_channel (str): input channel name
            out_channel (str): output channel name
        Returns:
            (:class:`numpy.ndarray`): cross-spectral density function
        """
        self._verify_channel(in_channel, "in_channel")
        self._verify_channel(out_channel, "out_channel")
        return self._ds["spectra"].sel(
            input=in_channel, output=out_channel).values.flatten()

    def _verify_channel(self, channel, ch_name):
        if not isinstance(channel, str):
            raise TypeError(f'{ch_name} is a {type(channel)}, not a str')
        if channel not in self.channels:
            raise ValueError(f'{ch_name} "{channel}" not in channels {self.channels}')

    def put_autospect(self, channel, auto_spect):
        """
        Equivalent to put_cross_spect(channel, channel, auto_spect)

        Args:
            channel (str): auto-spectra channel
            auto_spect (:class:`numpy.ndarray`): the auto-spectral density
        """
        if not auto_spect.shape == self.freqs.shape:
            raise ValueError('auto_spect has different shape than freqs ('
                             f'{auto_spect.shape} vs {self.freqs.shape})')
        assert auto_spect.dtype == 'complex'
        if channel not in self.channels:
            raise ValueError('channel "{channel}" is not in channels {self.channels}')
        self._ds["spectra"].loc[dict(input=channel, output=channel)] = \
            auto_spect

    def replace_channel_name(self, channel, replacement):
        """
        Args:
            channel (str): original channel name
            replacement (str): replacement channel name
        Example:
            >>> obj = SpectralDensity(('a','b','c'), np.array([0,1,2]),
                                      ('un', 'un', 'un'), 10, 'no')
            >>> obj.replace_channel_name('a', 'hello')
            >>> print(obj)
                SpectralDensity object:
	                channels=['a', 'b', 'c']
	                channel_units=['un', 'un', 'un']
	                3 frequencies, from 0 to 2Hz
	                n_windows=10
	                window_type=no

        """
        channel_names = self.channels
        channel_names[channel_names.index(channel)] = replacement
        self._ds['input']=channel_names
        self._ds['output']=channel_names

    def put_crossspect(self, in_channel, out_channel, cross_spect):
        """
        Put data into one of the cross-spectra.  Also puts the complex
        conjugate in the symmetric index

        Args:
            in_channel (str): cross-spectra input channel
            out_channel (str): cross-spectra output channel
            cross_spect (:class:`numpy.ndarray`): a cross-spectral density
        """
        assert cross_spect.shape == self.freqs.shape
        assert cross_spect.dtype == 'complex'
        assert in_channel in self.channels
        assert out_channel in self.channels
        self._ds["spectra"].loc[dict(input=in_channel, output=out_channel)] = \
            cross_spect
        if not in_channel == out_channel:
            self._ds["spectra"].loc[dict(input=out_channel, output=in_channel)] = \
                np.conj(cross_spect)


    def channel_response(self, channel):
        """
        A channel's instrument response
        
        Args:
            channel (str): channel name
        Returns:
            (:class:`numpy.ndarray`):
        """
        return self._ds["response"].sel(input=channel)

    def put_channel_response(self, channel, response):
        """
        Put a channel's instrument response into the object

        Verifies that the response has the same shape as the object's
        `frequency` property and that it is of type=`complex`
        
        Args:
            channel (str): the channel name
            response (:class:`numpy.ndarray`): the response
        """
        assert response.shape == self.freqs.shape
        assert response.dtype == 'complex'
        assert channel in self.channels
        self._ds["response"].loc[dict(input=channel)] = response

    def channel_units(self, channel):
        """
        Args:
            channel (str): the channel name
        Returns:
            (str): The input (physical) units of the given input or output channel
        """
        return str(self._ds["spectra"].sel(input=channel).coords['in_units'].values)

    def units(self, in_channel, out_channel):
        """
        The units of the given cross-  or auto-spectra

        Args:
            in_channel (str): input channel
            out_channel (str): output channel
        Returns:
            (str): the units
        """
        in_units = self.channel_units(in_channel)
        out_units = self.channel_units(out_channel)
        if in_units == out_units:
            return('({in_units})^2/Hz')
        return '({in_units})*({out_units})/Hz'

    @property
    def window_type(self):
        """
        The type of window used to calculate the spectral densities
        
        Returns:
            (str):
        """
        return(self._ds.window_type)

    @property
    def starttimes(self):
        """
        Start times for each data window used to calculate spectra
        
        Returns:
            (list of :class:`obspy.UTCDateTimes`): 
        """
        return(self._ds.starttimes)

    @property
    def n_windows(self):
        """
	    The number of data windows used to calculate spectra
	
        Returns:
            (int):
        """
        return(self._ds.n_windows)

#     @classmethod
#     def from_ATACR(cls, objnoise, horizontal_format='separate'):
#         """
#         Initiate class from ATACR DayNoise or StaNoise class
# 
#         Args:
#             objnoise (:class:`.DayNoise` or :class:`.StaNoise`):
#                 noise spectra and frequencies
#             horizontal_format (str): which type of horizontal channels to use:
#                 'aligned': one horizontal channel, aligned with highest noise
#                 'separate': two orthogonal horizontal channels
#         """
#         if (not objnoise and not isinstance(objnoise, DayNoise) and
#                 not isinstance(objnoise, StaNoise)):
#             raise TypeError("Error: A TFNoise object must be initialized with"
#                             " only one of type DayNoise or StaNoise object")
# 
#         if not objnoise.av:
#             raise(Exception("Error: Noise object has not been processed (QC "
#                             "and averaging) - aborting"))
# 
#         if horizontal_format == 'separate':
#             chans = ['1', '2', 'Z', 'P']
#             units = ['m/s^2', 'm/s^2', 'm/s^2', 'Pa']
#         elif horizontal_format == 'aligned':
#             chans = ['L', 'Z', 'P']
#             units = ['m/s^2', 'm/s^2', 'Pa']
#         else:
#             raise ValueError('horizontal_format not "separate" or "aligned"')
#         # shape = (len(chans), len(chans), len(objnoise.f))
#         if hasattr(objnoise, 'nwins'):
#             n_winds = np.sum(objnoise.nwins)
#         else:
#             n_winds = np.sum(objnoise.goodwins)
# 
#         obj = cls(chans, objnoise.f, units, n_winds, 'hanning')
#         if horizontal_format == 'separate':
#             obj.put_crossspect('1', '1', objnoise.power.c11)
#             obj.put_crossspect('1', '2', objnoise.cross.c12)
#             obj.put_crossspect('1', 'Z', objnoise.cross.c1Z)
#             obj.put_crossspect('1', 'P', objnoise.cross.c1P)
#             obj.put_crossspect('2', '1', np.conj(objnoise.cross.c12))
#             obj.put_crossspect('2', '2', objnoise.power.c22)
#             obj.put_crossspect('2', 'Z', objnoise.cross.c2Z)
#             obj.put_crossspect('2', 'P', objnoise.cross.c2P)
#             obj.put_crossspect('Z', '1', np.conj(objnoise.cross.c1Z))
#             obj.put_crossspect('Z', '2', np.conj(objnoise.cross.c2Z))
#             obj.put_crossspect('Z', 'Z', objnoise.power.cZZ)
#             obj.put_crossspect('Z', 'P', objnoise.cross.cZP)
#             obj.put_crossspect('P', '1', np.conj(objnoise.cross.c1P))
#             obj.put_crossspect('P', '2', np.conj(objnoise.cross.c2P))
#             obj.put_crossspect('P', 'Z', np.conj(objnoise.cross.cZP))
#             obj.put_crossspect('P', 'P', objnoise.power.cPP)
#         elif horizontal_format == 'aligned':
#             obj.put_crossspect('L', 'L', objnoise.rotation.cHH)
#             obj.put_crossspect('L', 'Z', objnoise.rotation.cHZ)
#             obj.put_crossspect('L', 'P', objnoise.rotation.cHP)
#             obj.put_crossspect('Z', 'L', np.conj(objnoise.rotation.cHZ))
#             obj.put_crossspect('Z', 'Z', objnoise.power.cZZ)
#             obj.put_crossspect('Z', 'P', objnoise.cross.cZP)
#             obj.put_crossspect('P', 'L', np.conj(objnoise.rotation.cHP))
#             obj.put_crossspect('P', 'Z', np.conj(objnoise.cross.cZP))
#             obj.put_crossspect('P', 'P', objnoise.power.cPP)
#         return obj

    @classmethod
    def from_stream(cls, stream, window_s=1000, windowtype='prol1pi',
                    inv=None, data_cleaner=None):
        """
        Calculate spectral density functions from the provided stream

        Should add a window selection algorithm, for now just steps by
        the window length

        Args:
            stream (:class:`obspy.core.stream.Stream`): data
            window_s (float): desired window length in seconds
            windowtype (str): window type, must be a valid
            inv (:class:`obspy.core.inventory.Inventory`): inventory containing
                instrument responses.  If none is found for the given channel,
                will look in the channel's stats.response object
            data_cleaner (:class:`DataCleaner`): Data cleaner to
                apply to channels as ffts are calculated
            subtract_tf_suffix (str): suffix to add to channel names if tf
                is subtracted
        """
        if not isinstance(stream, Stream):
            raise ValueError(f'stream is a {type(stream)}, not obspy Stream')
        stream = _align_traces(stream)

        # Select windows
        sr = stream[0].stats.sampling_rate
        ws = _npow2(window_s*sr)
        multfac = 2/(ws*sr)  # Bendata&Piersol 1986 eqs 11.100 & 11.102
        # window_starts = WindowSelect(stream, ws, windowtype)

        # Calculate FFTs
        ft, evalresps, units = {}, {}, []
        ids = [tr.id for tr in stream]
        if not len(ids) == len(set(ids)):
            raise ValueError('stream has duplicate IDs')
        for id in ids:  # Calculate Fourier transforms
            tr = stream.select(id=id)[0]
            ft[id], f = _calculate_windowed_rfft(tr, ws, ws, windowtype)
            n_winds = ft[id].shape[0]
	    	# Transform fft to physical units
            ft[id], resp, evalresp, ft_units = _correct_response(
				ft[id], f, id, tr.stats, inv)
            units.append(ft_units)
            evalresps[id] = evalresp
            if resp is not None:
                tr.stats.response = resp
        if data_cleaner is not None:  # clean data using correlated noise
            dctfs = data_cleaner.DCTFs
            old_ids = ids
            ft = dctfs.ft_subtract_tfs(ft)
            ids = dctfs.update_channel_names(old_ids)
            evalresps = dctfs.update_channel_keys(evalresps)
        # Create DataArray
        obj = cls(ids, f, units, n_winds, windowtype,
                  starttimes=[stream[0].stats.starttime])
        # Fill DataArray
        for inp in ids:
            if evalresps[inp] is not None:
                obj.put_channel_response(inp, evalresps[inp])
            for outp in ids:
                # This also puts the autospectra
                obj.put_crossspect(
					inp, outp, np.mean(ft[inp]*np.conj(ft[outp]), axis=0)*multfac)
        return obj

    def coherence(self, in_chan, out_chan):
        """
        The coherence for the given input and output channels

        Args:
            in_chan (str): input channel.  Must match one of the
                coordinates in _ds
            out_chan (str): output channel.  Must match one of the
                coordinates in _ds
        Returns:
            (:class:`numpy.ndarray`): Coherence absolute value

        Coherence is a real-valued quantity, for the cross-spectral phase,
        use the cross-spectral density function.
        From Bendat & Piersol (1986), eq 6.27, gamma^2_xy(f)
        """
        if in_chan not in self._ds.input:
            raise ValueError('"in_chan" not in spectral density matrix')
        if out_chan not in self._ds.output:
            raise ValueError('"out_chan" not in spectral density matrix')
        return np.abs(self.crossspect(in_chan, out_chan)**2
                      / (self.autospect(in_chan) * self.autospect(out_chan)))

    def coh_signif(self, prob=0.95):
        """
        The coherence significance level
        
        Args:
            prob (float): significance level (between 0 and 1)
        Returns:
            (float):
        """
        return coherence_significance_level(self.n_windows, prob)

    def plot(self, **kwargs):
        """Shortcut for `plot_autospectra()`"""
        self.plot_autospectra(**kwargs)

    def plot_autospectra(self, x=None, overlay=False, plot_peterson=True,
                         show=True, outfile=None, title=None):
        """
        Plot autospectra

        Args:
            x (list of str): limit to the listed channels
            overlay (bool): put all spect on one axis
            plot_peterson(bool): plot Peterson Noise model if any channel has
                units of (m/s^2)^2/Hz
            show (bool): show on desktop
            outfile (str): save figure to this filename
            title (str): custom plot title
        Returns:
            (:class:`numpy.ndarray`): array of axis pairs (amplitude, phase)
        """
        x = self._get_validate_channel_names(x)
        if not overlay:
            rows, cols = _squarish_grid(len(x))
        else:
            rows, cols = 1, 1
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        fig, axs = plt.subplots(rows, cols, sharex=True)
        if title is None:
            title = 'Auto-spectra'
        fig.suptitle(title)
        if not overlay:
            for key, i in zip(x, range(len(x))):
                i_row = int(i/cols)
                i_col = i - cols*i_row
                axa, axp = self.plot_one_spectra(key, key, fig, (rows, cols),
                                                 (i_row, i_col),
                                                 show_ylabel=i_col == 0,
                                                 show_xlabel=i_row == rows-1,
                                                 show_phase=False,
                                                 plot_peterson=plot_peterson)
                ax_array[i_row, i_col] = (axa, axp)
        else:
            axa, axp = None, None
            for key, i in zip(x, range(len(x))):
                axa, axp = self.plot_one_spectra(key, key, fig, (1, 1),
                                                 (0, 0),
                                                 show_ylabel=i == len(x)-1,
                                                 show_xlabel=i == len(x)-1,
                                                 ax_a=axa, ax_p=axp,
                                                 show_phase=False,
                                                 plot_peterson=plot_peterson)
            ax_array[0, 0] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_cross_spectra(self, x=None, show=True, show_coherence=False,
                           outfile=None, plot_peterson=False):
        """
        Plot cross (and auto) spectra

        Args:
            x (list of str): limit to the listed channels
            show (bool): show on desktop
            plot_peterson(bool): plot Peterson Noise model if any channel has
                units of (m/s^2)^2/Hz
            show_coherence (bool): show coherence as well
        Returns:
            :class:`numpy.ndarray`: array of axis pairs (amplitude, phase)
        """
        x = self._get_validate_channel_names(x)
        n_subkeys = len(x)
        rows, cols = n_subkeys, n_subkeys
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        fig, axs = plt.subplots(rows, cols, sharex=True)
        fig.suptitle('Cross-spectra (dB ref UNITS/Hz)')
        for in_chan, i in zip(x, range(len(x))):
            for out_chan, j in zip(x, range(len(x))):
                title = out_chan if i == 0 else None
                axa, axp = self.plot_one_spectra(in_chan, out_chan, fig,
                                                 (rows, cols), (i, j),
                                                 show_ylabel=j == 0,
                                                 show_xlabel=i == rows-1,
                                                 ylabel=in_chan,
                                                 label='units',
                                                 title=title,
                                                 show_coherence=show_coherence)
                ax_array[i, j] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_one_autospectra(self, key, **kwargs):
        """
        Plot one autospectral density
        
        Arguments are the same as for `plot_one_spectra()`, except
        there is no `subkey` argument and `show_phase` is ignored
        """
        kwargs['show_phase'] = False
        self.plot_one_spectra(key, key, **kwargs)
        
    def plot_one_spectra(self, key, subkey, fig=None, fig_grid=(1, 1),
                         plot_spot=(0, 0), show_xlabel=True, show_ylabel=None,
                         ax_a=None, ax_p=None, ylabel=None, label=None,
                         title=None, show_coherence=False,
                         show_phase=True, plot_peterson=True, outfile=None):
        """
        Plot one spectral density

        Args:
            key (str): input (driving) channel
            subkey (str): output (response) channel
            fig (:class:`matplotlib.figure.Figure`): figure to plot on, if
                None this method will plot on the current figure or create
                a new figure.
            fig_grid (tuple): this plot sits in a grid of this many
                              (rows, columns)
            subplot_spot (tuple): put this plot at this (row,column) of
                                  the figure grid
            show_xlabel (bool): put an xlabel on this subplot
            show_ylabel (bool): put a ylabel on this subplot
            ylabel (str): label to put on y axis (if show_label).  If not
                speficied, will use 'dB ref UNITS/Hz'
            label (str): 'units': print units without channel name in legend
            ax_a (Axis): use an existing axis for the amplitude plot
            ax_p (Axis): use this existing axis for the phase plot
            title (str): title to put on this subplot
            show_coherence (bool): draw coherence on the same plot
            show_phase (bool): show phase as well as amplitude
            plot_peterson(bool): plot Peterson Noise model if channel has
                units of (m/s^2)^2/Hz
            outfile (str): save figure to this filename

        Returns:
            (tuple): tuple containing
                - :class:`matplotlib.axes.axis`: amplitude plot axis
                - :class:`matplotlib.axes.axis`: phase plot axis
        """
        # da = self._ds["spectra"].sel(input=key, output=subkey)
        # in_units = da.coords['in_units'].values
        # out_units = da.coords['out_units'].values
        psd = self.crossspect(key, subkey)
        in_units = self.channel_units(key)
        out_units = self.channel_units(subkey)
        if in_units == out_units:
            PSD_units = f'({in_units})^2'
        else:
            PSD_units = f'{in_units}*{out_units}'
        f = self.freqs
        if fig is None:
            fig = plt.gcf()
        # Plot amplitude
        if ax_a is None:
            if show_phase:
                ax_a = plt.subplot2grid((3*fig_grid[0], 1*fig_grid[1]),
                                        (3*plot_spot[0]+0, plot_spot[1]+0),
                                        rowspan=2)
            else:
                ax_a = plt.subplot2grid((fig_grid[0], fig_grid[1]),
                                        (plot_spot[0], plot_spot[1]))
        if show_coherence:
            ax2 = ax_a.twinx()
            ax2.semilogx(f, np.abs(self.coherence(key, subkey)),
                         color='red', linewidth=0.5, alpha=0.8)
            ax2.axhline(self.coh_signif(0.95), color='red',
                        linewidth=0.5, alpha=0.8, ls='--')
            ax2.set_ylim(0, 1)
            if plot_spot[1] == fig_grid[1]-1:  # Rightmost column
                ax2.set_ylabel('Coher', color='red')
            else:
                ax2.set_yticklabels([])
        psd[psd == 0] = None

        # Plot amplitude
        if label is None:
            label = f'{subkey} ({PSD_units})'
        elif label == 'units':
            label = f'{PSD_units}'
        ax_a.semilogx(f, 10*np.log10(np.abs(psd)), label=label)
        ax_a.set_xlim(f[1], f[-1])
        if plot_peterson is True and PSD_units.lower() == '(m/s^2)^2':
            lownoise, highnoise = Peterson_noise_model(f, True)
            ax_a.semilogx(f, lownoise, 'k--')
            ax_a.semilogx(f, highnoise, 'k--')

        if label is not None:
            legend_1 = ax_a.legend()
            if show_coherence:
                legend_1.remove()
                ax2.add_artist(legend_1)
        if show_ylabel:
            if ylabel is None:
                ylabel = 'dB ref UNITS/Hz'
            ax_a.set_ylabel(ylabel)
        if title:
            ax_a.set_title(title)

        # Plot phase
        if show_phase:
            if ax_p is None:
                ax_p = plt.subplot2grid((3*fig_grid[0], 1*fig_grid[1]),
                                        (3*plot_spot[0]+2, plot_spot[1]+0))
            ax_p.semilogx(f, np.degrees(np.angle(psd)))
            ax_p.set_ylim(-180, 180)
            ax_p.set_xlim(f[1], f[-1])
            ax_p.set_yticks((-180, 0, 180))
            if show_ylabel:
                # ax_p.set_ylabel('Phase')
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
        return ax_a, ax_p

    def plot_coherences(self, x=None, y=None, overlay=False, show=True,
                        outfile=None):
        """
	    Plot coherences

        Args:
            x (list of str): limit to the listed input channels
            y (list of str): limit to the listed output channels
            overlay (bool): put all coherences on one plot
            show (bool): show on desktop
            outfile (str): save to the named file

        Returns:
            (:class:`numpy.ndarray`): array of axis pairs (amplitude, phase)
        """
        x = self._get_validate_channel_names(x)
        y = self._get_validate_channel_names(y)
        if not overlay:
            rows, cols = len(x), len(y)
            ax_array = np.ndarray((rows, cols), dtype=tuple)
            fig, axs = plt.subplots(rows, cols, sharex=True)
            fig.suptitle('Coherences')
            for in_chan, i in zip(x, range(len(x))):
                for out_chan, j in zip(y, range(len(y))):
                    title = out_chan if i == 0 else None
                    axa, axp = self.plot_one_coherence(in_chan, out_chan, fig,
                                                       (rows, cols), (i, j),
                                                       show_ylabel=j == 0,
                                                       show_xlabel=i == rows-1,
                                                       ylabel=in_chan,
                                                       title=title)
                    ax_array[i, j] = (axa, axp)
        else:
            ax_array = np.ndarray((1, 1), dtype=tuple)
            fig, axs = plt.subplots(1, 1, sharex=True)
            fig.suptitle('Coherences')
            labels = []
            axa, axp = None, None
            for in_chan, i in zip(x, range(len(x))):
                for out_chan, j in zip(y, range(len(y))):
                    if out_chan == in_chan:
                        break
                    if f'{out_chan}-{in_chan}' in labels:
                        break
                    label = f'{in_chan}-{out_chan}'
                    axa, axp = self.plot_one_coherence(in_chan, out_chan, fig,
                                                       (1, 1), (0, 0),
                                                       ylabel='Coherence',
                                                       label=label,
                                                       ax_a=axa, ax_p=axp)
                    labels.append(label)
            ax_array[0, 0] = (axa, axp)
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        return ax_array

    def plot_one_coherence(self, in_chan, out_chan, fig=None, fig_grid=(1, 1),
                           plot_spot=(0, 0), show_xlabel=True,
                           show_ylabel=True, ax_a=None, ax_p=None,
                           ylabel=None, label=None,
                           title=None, show_phase=True):
        """
        Plot one coherence

        Args:
            in_chan (str): input (driving) channel
            out_chan (str): output (response) channel
            fig (:class:`matplotlib.figure.Figure`): figure to plot on, if
                None this method will plot on the current figure or create
                a new figure.
            fig_grid (tuple): this plot sits in a grid of this many
                              (rows, columns)
            plot_spot (tuple): put this plot at this (row,column) of
                                  the figure grid
            show_xlabel (bool): put an xlabel on this subplot
            show_ylabel (bool): put a ylabel on this subplot
            ylabel (str): label to put on y axis (if show_label).  If not
                speficied, will use 'dB ref UNITS/Hz'
            label (str): text to put in legend
            ax_a (Axis): use an existing axis for the amplitude plot
            ax_p (Axis): use this existing axis for the phase plot
            title (str): title to put on this subplot
            show_phase (bool): show phase as well as amplitude

        Returns:
            (tuple): tuple containing:
                - (:class:`matplotlib.axes.axis`): amplitude plot axis
                - (:class:`matplotlib.axes.axis`): phase plot axis
        """
        ds = self._ds['spectra'].sel(input=in_chan, output=out_chan)
        f = self._ds.coords['f'].values
        if fig is None:
            fig = plt.gcf()
        # Plot amplitude
        if ax_a is None:
            if show_phase:
                ax_a = plt.subplot2grid((3*fig_grid[0], 1*fig_grid[1]),
                                        (3*plot_spot[0]+0, plot_spot[1]+0),
                                        rowspan=2)
            else:
                ax_a = plt.subplot2grid((fig_grid[0], fig_grid[1]),
                                        (plot_spot[0], plot_spot[1]))
        ax_a.semilogx(f, np.abs(self.coherence(in_chan, out_chan)),
                      label=label)
        ax_a.axhline(self.coh_signif(0.95), color='red', linewidth=0.5,
                     alpha=0.8, ls='--')
        ax_a.set_ylim(0, 1)
        # ax_a.set_yticklabels([])
        # Plot amplitude
        if label is not None:
            ax_a.legend()
        if show_ylabel:
            if ylabel is None:
                ylabel = 'Coherence'
            ax_a.set_ylabel(ylabel)
        if title:
            ax_a.set_title(title)

        # Plot phase
        if show_phase:
            if ax_p is None:
                ax_p = plt.subplot2grid((3*fig_grid[0], 1*fig_grid[1]),
                                        (3*plot_spot[0]+2, plot_spot[1]+0))
            ax_p.semilogx(f, np.degrees(np.angle(ds)))
            ax_p.set_ylim(-180, 180)
            ax_p.set_xlim(f[1], f[-1])
            ax_p.set_yticks((-180, 0, 180))
            if show_ylabel:
                # ax_p.set_ylabel('Phase')
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
        return ax_a, ax_p

    def _get_validate_channel_names(self, x):
        """
        If x is None, return list of all channel names
        If x is a list, validate all of the names
        """
        if x is None:
            return list(self._ds.coords['input'].values)
        for key in x:
            if key not in list(self._ds.coords['input'].values):
                ValueError('key "{key}" not in channel list')
        return x

    @staticmethod
    def _remove_subtracted_loc(id):
        """Remove loc code characters including and after first "-"
    
        Allows the use of "-?" in the loc code to specify removed coherent noise

        Args:
            id (str): seed ID code
    
        Example:
            >>> SD._remove_subtracted_loc('hello')
            'hello'
            >>> SD._remove_subtracted_loc('NN.SSSS.LL.CCC')
            'NN.SSSS.LL.CCC'
            >>> SD._remove_subtracted_loc('NN.SSSS.-LL.CCC')
            'NN.SSSS..CCC'
            >>> SD._remove_subtracted_loc('NN.SSSS.L-LL.CCC')
            'NN.SSSS.L.CCC'
        """
        comps = id.split('.')
        if not len(comps) == 4:
            return id
        if not '-' in comps[2]:
            return id
        comps[2]=comps[2].partition('-')[0]
        return '.'.join(comps)
    

def _align_traces(stream):
    """Trim stream so that all traces are aligned and same length"""
    # Determine last starttime and first endtime, and verify that
    # all traces have the same sample rate
    first_start = last_start = stream[0].stats.starttime
    first_end = last_end = stream[0].stats.endtime
    sampling_rate = stream[0].stats.sampling_rate
    for tr in stream[1:]:
        if tr.stats.starttime > last_start:
            last_start = tr.stats.starttime
        elif tr.stats.starttime < first_start:
            first_start = tr.stats.starttime
        if tr.stats.endtime < first_end:
            first_end = tr.stats.endtime
        elif tr.stats.endtime > last_end:
            last_end = tr.stats.endtime
        if not tr.stats.sampling_rate == sampling_rate:
            raise ValueError("not all traces have same sample rate")

    if last_start >= first_end:
        raise ValueError("There are non-overlapping traces")
    if last_start - first_start > 1/sampling_rate:
        print("Cutting up to {last_start-first_start}s from trace starts")
    if last_end - first_end > 1/sampling_rate:
        print("Cutting up to {last_start-first_start}s from trace ends")
    stream.trim(last_start, first_end)
    min_len = min([tr.stats.npts for tr in stream])
    max_len = max([tr.stats.npts for tr in stream])
    if not max_len == min_len:
        for tr in stream:
            tr.data = tr.data[:min_len]
    return stream


# COPIED FROM ATACR, but other tapers added/allowed and assumes data are real
def _calculate_windowed_rfft(trace, ws, ss=None, win_taper='hanning'):
    """
    Calculates windowed Fourier transform of real data

    Args:
        trace (:class:`obspy.core.Trace`): Input trace data
        ws (int): Window size, in number of samples
        ss (int): Step size, or number of samples until next window
        win_taper (str): taper to apply to data ['hanning', 'prol4pi',
            'prol1pi']

    Returns:
        (tuple): tuple containing:
            ft (:class:`numpy.ndarray`): Fourier transform of trace
            f (:class:`numpy.ndarray`): Frequency axis in Hz
    """
    # Extract sliding windows
    a, nd = _sliding_window(trace.data, ws, ss, win_taper)
    # Fourier transform
    n2 = _npow2(ws)
    ft = np.fft.rfft(a, n=n2)
    f = np.fft.rfftfreq(ws, 1./trace.stats.sampling_rate)
    # f = np.linspace(0., 1., int(n2/2) + 1) * trace.stats.sampling_rate/2.
    # Don't return zero frequency
    return ft[:, 1:], f[1:]


def _sliding_window(a, ws, ss=None, win_taper='hanning'):
    """
    Split a data array into overlapping, tapered sub-windows

    Args:
        a (:class:`numpy.ndarray`): 1D array of data to split
        ws (int): Window size in samples
        ss (int): Step size in samples. If not provided, window and step size
            are equal.
        win_taper (str): taper to apply to data ['hanning', 'prol4pi',
            'prol1pi', 'bartlett', 'blackman', 'hamming']

    Returns:
        (tuple): tuple conaining
            out (:class:`numpy.ndarray`): 1D array of windowed data
            nd (int): Number of windows
    """
    if ss is None:
        # no step size was provided. Return non-overlapping windows
        ss = ws
    ws = int(ws)
    # Calculate the number of windows to return, ignoring leftover samples, and
    # allocate memory to contain the samples
    valid = len(a) - ss
    nd = (valid) // ss
    out = np.ndarray((nd, ws), dtype=a.dtype)
    if win_taper in ['hanning', 'hamming', 'blackman', 'bartlett']:
        taper = eval(f'np.{win_taper}(ws)')
    elif win_taper == 'prol1pi':
        taper = _prol1pi(ws)
    elif win_taper == 'prol4pi':
        taper = _prol4pi(ws)
    else:
        raise ValueError(f'Unknown taper type "{win_taper}"')
    if nd == 0:
        out = signal.detrend(a) * taper
    for i in range(nd):
        # "slide" the window along the samples
        start = i * ss
        stop = start + ws
        out[i] = signal.detrend(a[start: stop]) * taper
    return out, nd


def _npow2(x):
    return 1 if x == 0 else 2**int(x-1).bit_length()


def _squarish_grid(n_elems):
    """
    Returns as close to a square grid as fits the data without waste

    Args:
        n_elems (int): Number of elements in the grid

    Returns:
        (tuple): tuple containing:
            rows (int): number of rows in grid
            cols (int): number of columns in grid
    """
    n_elems = int(n_elems)
    if n_elems == 1:
        return 1, 1
    elif n_elems == 2:
        return 1, 2
    elif n_elems <= 4:
        return 2, 2
    elif n_elems <= 6:
        return 2, 3
    else:
        cols = np.ceil(np.sqrt(n_elems))
        rows = np.ceil(n_elems/cols)
    return rows, cols

def _subtract_tfs(fts, subtract_tfs):
    """
    Args:
        fts (dict): dictionary containing Fourier transforms for each channel.
            Each Fourier transform is N*ws, where ws is the window size and N
            is the mumber of windows
        subtract_tfs (list of DCTF): transfer functions to
            subtract from channels as ffts are calculated (NOT SURE IF
            THIS IS MORE USEFUL THAN SUBTRACTING FROM THE FINAL
            SPECTRALDENSITY (LIKE ATACR), BUT NEED TO TEST)
    Returns:
        fts (dict): dictionary containg corrected Fourier transforms for each
            channel.
    """
    for stf in subtract_tfs:
        tfs = stf.tfs
        in_chan = tfs.input_channel
        fts_ic = in_chan.split('-')[0] # take off any '-*' tail
        if not tfs.freqs.shape == fts[fts_ic].shape:
            ValueError('transfer function and ft have different shapes ({} vs {})'
                       .format(tfs.freqs.shape, fts[fts_ic].shape))
        for out_chan in tfs.output_channels:
            fts_oc = out_chan.split('-')[0]
            fts[fts_oc] -= fts[fts_ic] * tfs.values(out_chan)
    return fts

def _correct_response(ft, f, id, stats, inv=None):
    """
    Convert fourier transform into channel in_units
    
    Args:
        ft (:class:`numpy.ndarray`): fourier transforms for one channel
        f (:class:`numpy.ndarray`): frequencies
        id (str): channel id
        tr (:class:`obspy.core.stream.stats`): trace statistics
        inv (:class:`obspy.core.inventory.Inventory`): station inventory
    Returns:
        ft: corrected fourier transforms
        response: the channel response
        evalresp: the channel response evaluated at frequencies f
        units: the channel response units
    """
    resp, evalresp, units = None, None, 'Counts'
    if inv is not None:
        try:
            resp = inv.get_response(id, stats.starttime)
        except Exception:
            new_id = SpectralDensity._remove_subtracted_loc(id)
            try:
                resp = inv.get_response(new_id, stats.starttime)
            except Exception:
                raise ValueError(f'No match found for "{new_id}" in inv')
    if resp is None and 'response' in stats:
        resp = stats.response
    if resp is not None:
        if "pa" in resp.instrument_sensitivity.input_units.lower():
            evalresp = resp.get_evalresp_response_for_frequencies(f, 'VEL')
            units = 'Pa'
        else:
            evalresp = resp.get_evalresp_response_for_frequencies(f, 'ACC')
            units = 'm/s^2'
        ft /= evalresp
    return ft, resp, evalresp, units

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
