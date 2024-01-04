from dataclasses import dataclass
from collections import UserList

import numpy as np
from matplotlib import pyplot as plt

from ..response_functions import ResponseFunctions
# from ..utils import CleanSequence as CS
from ..logger import init_logger

logger = init_logger()


class RFList(UserList):
    """
    list of response_functions, used by Data Cleaner

    Most of the code is about channel renaming based on the removed channels
    """

    def ft_subtract_rfs(self, fts, evalresps=None):
        """
        Remove effect of channels on each other.
        
        Args:
            fts (dict): Fourier transforms for each channel.
                        Each Fourier transform is N*ws, where ws is the
                        window size and N is the number of windows
                        The key for each item is the channel name
            evalresps (dict): instrument responses for each channel (if
                fts are response-corrected)
        Returns:
            (tuple):
                fts (dict): corrected Fourier transforms for each channel.
                clean_sequence_dict (dict of list): clean_sequence list for each channel
        """
        clean_sequence_dict = {}
        for rfs in self:
            id_in = rfs.input_channel
            if not rfs.freqs.shape == fts[id_in].shape:
                ValueError("frequency response function and ft have different shapes "
                           f"({rfs.freqs.shape} vs {fts[id_in].shape})")
            for id_out in rfs.output_channels:
                # Verify that the fourier transforms have the same shape
                if not fts[id_out].shape == fts[id_in].shape:
                    ValueError(f"ft[{id_in=}] and ft[{id_out=}] "
                               "have different shapes ({} vs {})".response(
                               fts[id_in].shape, fts[id_out].shape))
                # Clean coherent noise (EQ 8, Crawford & Webb 2000)
                multiplier = np.ones(fts[id_in].shape)
                if evalresps is not None:
                    if evalresps[id_in] is not None and evalresps[id_out] is not None:
                            multiplier = evalresps[id_in] / evalresps[id_out]
                fts[id_out] -= fts[id_in] * multiplier * rfs.corrector(id_out)
                # Add input channel to clean_sequence for output channel
                if id_out not in clean_sequence_dict:
                    clean_sequence_dict[id_out] = [id_in]
                else:
                    clean_sequence_dict[id_out].append(id_in)
                # print(f'clean_sequence_dict[{id_out}]={clean_sequence_dict[id_out]}')
        # print(f'{clean_sequence_dict=}')
        return fts, clean_sequence_dict

    def __str__(self):
        s = "RFList object:\n"
        s += "   Step | Channel to remove  | Channels to remove from\n"
        s += "   ==== | ================== | ========================================\n"
        for rfs, i in zip(self, range(len(self))):
            s += " {:>6d} | {:18s} | {}\n".format(
                i+1, rfs.input_channel,
                "['" +"', '".join(rfs.output_channels) + "']")
        return s

    def plot(self, outfile=None, title=None):
        """
        Plot frequency response functions

        Args:
            outfile (str): output file name.  If None, will plot to screen
            title (str): plot title.  If None, will use self[0].rfs.input_channel
        """
        ncols = max([len(rfs.output_channels) for rfs in self])
        nrows = len(self)
        # fig, axs = plt.subplots(nrows, ncols)
        fig, axs = plt.subplots()
        row, col = 0, 0
        for rfs in self:
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
        seed_idl = self[0].rfs.input_channel.split(".")
        if title is not None:
            fig.suptitle(title)
        elif len(seed_idl) > 1:
            title = ".".join(seed_idl[:2])
        else:
            title = self[0].rfs.input_channel
        fig.suptitle(f"{title} DataCleaner Frequency Response Functions")
        if outfile:
            fig.savefig(outfile)
        else:
            plt.show()
