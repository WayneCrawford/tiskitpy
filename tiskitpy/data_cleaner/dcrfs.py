from dataclasses import dataclass
from collections import UserList

import numpy as np
from matplotlib import pyplot as plt

from ..response_functions import ResponseFunctions
from .cleaner_string import CleanerString as CS


@dataclass
class DCRF:
    """
    Data Cleaner Response Function class

    Organizes and names ResponseFunctions for DataCleaner
    
    Attributes:
        remove_channel: input channel name
        remove_sequence: 
    """

    remove_channel: str
    remove_sequence: str
    rfs: ResponseFunctions


class DCRFs(UserList):
    """
    list of DCRFs

    Most of the code involves channel renaming based on the removed channels
    """

    def ft_subtract_rfs(self, fts, evalresps=None):
        """
        Remove effect of channels on each other.  Modifies channel names
        to reflect the response function sequence.
        Args:
            fts (dict): Fourier transforms for each channel.
                        Each Fourier transform is N*ws, where ws is the
                        window size and N is the number of windows
            evalresps (dict): instrument responses for each channel (if
                fts are response-corrected)
        Returns:
            (tuple):
                fts (dict): corrected Fourier transforms for
                    each channel.
                evalresps (dict): dict with modified keys
        """
        for dcrf in self:
            rfs = dcrf.rfs
            # mapping = dict()
            id_in = rfs.input_channel
            if not rfs.freqs.shape == fts[id_in].shape:
                ValueError("frequency response function and ft have different shapes "
                           f"({rfs.freqs.shape} vs {fts[id_in].shape})")
            for id_out in rfs.output_channels:
                if not fts[id_out].shape == fts[id_in].shape:
                    ValueError(f"ft[{id_in=}] and ft[{id_out=}] "
                               "have different shapes "
                               f"({fts[id_in].shape} vs {fts[id_out].shape})")
                new_id_out = CS.insert(dcrf.remove_sequence, CS.strip(id_out))
                # new_id_out = (CS.strip(id_out)
                #               + dcrf.remove_sequence)
                # EQ 8, Crawford & Webb 2000
                multiplier = np.ones(fts[id_in].shape)
                if evalresps is not None:
                    oid = CS.strip(id_out)
                    iid = CS.strip(id_in)
                    if evalresps[iid] is not None:
                        if evalresps[oid] is not None:
                            multiplier = evalresps[iid] / evalresps[oid]
                fts[new_id_out] = (fts[id_out] - fts[id_in]
                                   * multiplier * rfs.corrector(id_out))
        return fts

    def update_channel_names(self, channel_names):
        """
        Change an array's channel names to include final channel removal codes

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
        Plot frequency response functions

        Args:
            outfile (str): output file name.  If None, will plot to screen
            title (str): plot title.  If None, will use self[0].remove_channel
        """
        ncols = max([len(dcrf.rfs.output_channels) for dcrf in self])
        nrows = len(self)
        # fig, axs = plt.subplots(nrows, ncols)
        fig, axs = plt.subplots()
        row, col = 0, 0
        for dcrf in self:
            rfs = dcrf.rfs
            for oc in rfs.output_channels:
                ic = rfs.input_channel
                rfs.plot_one(
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
        fig.suptitle(f"{title} DataCleaner Frequency Response Functions")
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
        for dcrf in self:
            rfs = dcrf.rfs
            # in_chan = rfs.input_channel
            for out_chan in rfs.output_channels:
                orig_out = CS.strip(out_chan)
                # mapping[orig_out] = orig_out + dcrf.remove_sequence
                mapping[orig_out] = CS.insert(dcrf.remove_sequence, orig_out)
        return mapping
