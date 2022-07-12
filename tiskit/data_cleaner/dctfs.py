from dataclasses import dataclass
from collections import UserList

import numpy as np
from matplotlib import pyplot as plt

from ..transfer_functions import TransferFunctions


class DCTFs(UserList):
    """
    list of DCTFs
    """

    def ft_subtract_tfs(self, fts):
        """
        Args:
            fts (dict): dictionary containing Fourier transforms for each
                channel. Each Fourier transform is N*ws, where ws is the
                window size and N is the mumber of windows
        Returns:
            fts (dict): dictionary containg corrected Fourier transforms for
                each channel.
        """
        for dctf in self:
            tfs = dctf.tfs
            # mapping = dict()
            in_chan = tfs.input_channel
            if not tfs.freqs.shape == fts[in_chan].shape:
                ValueError(
                    "transfer function and ft have different shapes "
                    f"({tfs.freqs.shape} vs {fts[in_chan].shape})"
                )
            for out_chan in tfs.output_channels:
                new_out_chan = (
                    strip_remove_str(out_chan) + dctf.remove_sequence
                )
                # EQ 8, Crawford & Webb 2000
                fts[new_out_chan] = fts[out_chan] - fts[in_chan] * np.conj(
                    tfs.values(out_chan)
                )
        return fts

    def update_channel_names(self, channel_names):
        """
        Change an array's channel names to final channel removal codes

        Args:
            channel_names (list): list of valid channel names.

        Returns:
            new_channel_names (dict): dict whose keys are updated channel names
        """
        mapping = self._channel_map()
        new_channel_names = channel_names.copy()
        for k, v in mapping.items():
            if k not in new_channel_names:
                raise ValueError(f'"{k}" not in channel_names')
            new_channel_names[channel_names.index(k)] = v
        return new_channel_names

    def update_channel_keys(self, in_dict):
        """
        Adds final channel removal codes to a dict's keys

        Args:
            in_dict (dict): dictionary whose keys are valid channel names

        Returns:
            out_dict (dict): dict whose keys are updated channel names
        """
        mapping = self._channel_map()
        out_dict = in_dict.copy()
        for k, v in mapping.items():
            if k not in out_dict.keys():
                raise ValueError(f'"{k}" not in in_dict keys')
            out_dict[v] = out_dict[k]
        return out_dict

    def plot(self, outfile=None, title=None):
        """
        Plot transfer functions

        Args:
            outfile (str): output file name.  If None, will plot to screen
            title (str): plot title.  If None, will use self[0].remove_channel
        """
        ncols = max([len(dctf.tfs.output_channels) for dctf in self])
        nrows = len(self)
        # fig, axs = plt.subplots(nrows, ncols)
        fig, axs = plt.subplots()
        row, col = 0, 0
        for dctf in self:
            tfs = dctf.tfs
            for oc in tfs.output_channels:
                ic = tfs.input_channel
                tfs.plot_one(
                    ic,
                    oc,
                    fig=fig,
                    fig_grid=(nrows, ncols),
                    plot_spot=(row, col),
                    label=oc.split(".")[-1] + "/" + ic.split(".")[-1],
                    show_ylabel=col == 0,
                    show_xlabel=row == nrows - 1,
                )
                col += 1
            row += 1
            col = 0
        seed_idl = self[0].remove_channel.split(".")
        if title is not None:
            fig.suptitle(title)
        elif len(seed_idl) > 1:
            title = ".".join(seed_idl[:2])
        else:
            title = self[0].remove_channel
        fig.suptitle(f"{title} DataCleaner Transfer Functions")
        if outfile:
            fig.savefig(outfile)
        else:
            plt.show()

    def _channel_map(self):
        """
        Make a mapping dict between original and output channel names

        Args:
            channel_names (list): list of valid channel names

        Returns:
            mapping (dict): keys= old channel names, values = new channel names
        """
        mapping = dict()
        for dctf in self:
            tfs = dctf.tfs
            # in_chan = tfs.input_channel
            for out_chan in tfs.output_channels:
                orig_out = strip_remove_str(out_chan)
                mapping[orig_out] = orig_out + dctf.remove_sequence
        return mapping


@dataclass
class DCTF:
    """
    Data Cleaner Transfer Function class

    Organizes and names TransferFunctions for DataCleaner
    """

    remove_channel: str
    remove_sequence: str
    tfs: TransferFunctions


def strip_remove_str(in_str):
    """
    Strips transfer function removal string from a string

    Just takes out everything after first '-'

    Also handles lists of strings

    Example:
        >>> text = 'BHZ-1-2-H'
        >>> strip_remove_str(text)
        'BHZ'
        >>> text = 'BHZ'
        >>> strip_remove_str(text)
        'BHZ'
    """
    if isinstance(in_str, list):
        return [strip_remove_str(x) for x in in_str]
    return in_str.split("-")[0]


def strip_remove_one(in_str):
    """
    Strips last transfer function removal string from a string

    Also handles lists of strings

    Example:
        >>> text = 'BHZ-1-2-H'
        >>> strip_remove_one(text)
        'BHZ-1-2'
        >>> text = 'BHZ'
        >>> strip_remove_one(text)
        'BHZ'
    """
    if isinstance(in_str, list):
        return [strip_remove_one(x) for x in in_str]
    return in_str.rsplit("-", 1)[0]


def remove_str(in_chan, prev_remove_seq=""):
    """
    Generate a string describing the TransferFunction removal chain

    Args:
        in_chan: name of channel to remove
        prev_remove_seq (str): previous removal sequence string, used to
            avoid repeating an already-given channel code
    Returns:
        remove_new (str): string to add to end of removal sequence

    Example:
        >>> remove_str('BHZ', '-X-Y')
        '-Z'
        >>> remove_str('BHX', '-X-Y')
        '-HX'
        >>> remove_str('X', '-X-Y')
        Traceback (most recent call last):
        ...
        ValueError: new_name('X') already in prev_remove_seq('-X-Y')
    """
    prev_components = prev_remove_seq.split("-")
    offset = 1
    while (new_name := in_chan[-offset:]) in prev_components:
        offset += 1
        if offset > len(in_chan):
            raise ValueError(
                f"new_name('{new_name}') already in "
                f"prev_remove_seq('{prev_remove_seq}')"
            )
    remove_elem = "-" + new_name
    return remove_elem
