#!env python3
""" FIR filter class"""

import inspect
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve, freqz, lfilter
from obspy.core.inventory.response import FIRResponseStage

from ..logger import init_logger

logger = init_logger()


@dataclass
class FIRFilter:
    """Create and implement FIR filters"""

    coeffs: list
    offset: int = 0
    decim: int = 1
    name: str = ""

    @classmethod
    def from_SAC(cls, decim, normalize=True):
        """
        Make SAC FIR filter

        Args:
            decim (int): decimation factor (must be one of 2, 3, 4, 5, 6, 7)
            renormalize (bool): renormalize filter to have unit sum (avoids
                                obspy warning, but is it right?)
        """
        if decim not in (2, 3, 4, 5, 6, 7):
            raise ValueError(f"No SAC decimation filter for {decim=}")
        SAC_path = (
            Path(inspect.getfile(inspect.currentframe())).resolve().parent
            / "sac_fir"
        )
        with open(str(SAC_path / f"dec{decim}")) as f:
            f.readline()
            A = f.readline().split()
            divisor, n_values = float(A[0]), int(A[1])
            coeffs = [
                float(x) / divisor
                for line in f.readlines()
                for x in line.split()
            ]
        if not len(coeffs) == n_values:
            raise RuntimeError(
                "# FIR values is different from that stated: "
                f"{len(coeffs)} != {n_values}"
            )
        # SAC only saves right half of symmetrical odd filter
        offset = len(coeffs) - 1
        coeffs = coeffs[-1:0:-1] + coeffs
        if normalize:
            if abs(np.sum(coeffs) - 1) > 0.001:
                divisor = np.sum(coeffs)
                coeffs = [x / divisor for x in coeffs]
        return cls(coeffs, offset, decim, f"SAC_decim{decim}")

    def __sub__(self, y):
        if not len(self.coeffs) == len(y.coeffs):
            raise ValueError("Cannot subtract filters with different lengths")
        return FIRFilter([x - y for x, y in zip(self.coeffs, y.coeffs)])

    def __eq__(self, y):
        if not len(self.coeffs) == len(y.coeffs):
            return False
        if not self.coeffs == y.coeffs:
            return False
        if not self.offset == y.offset:
            return False
        if not self.decim == y.decim:
            return False
        if not self.name == y.name:
            return False
        return True
        return FIRFilter([x - y for x, y in zip(self.coeffs, y.coeffs)])

    def convolve(self, data):
        """
        convolve data with FIR filter and return result
        """
        con_data = convolve(data, self.coeffs, mode="full")
        cut_data = con_data[
            self.offset: -(len(self.coeffs) - self.offset - 1)
        ]
        assert len(data) == len(cut_data)
        return cut_data

    def _apply_symmetry(self, symmetry_code):
        """
        Apply symmetry_code to make a full filter with "Symmetry=NONE"
        """
        if symmetry_code == "NONE":
            return
        elif symmetry_code == "ODD":
            self.coeffs = self.coeffs + self.coeffs[-2::-1]
        elif symmetry_code == "EVEN":
            self.coeffs = self.coeffs + self.coeffs[-1::-1]
        else:
            raise NameError(f"FIR symmetry = {symmetry_code} not implemented")
        return

    def to_obspy(self, sampling_rate, stage_sequence_number=0,
                 prior_output_units=False):
        """
        Return an obspy FIRResponseStage
        """
        units = 'count'     # StationXML recommendation
        if prior_output_units is not False:
            if prior_output_units[:5].lower() == 'count':
                units = prior_output_units
        return FIRResponseStage(
            stage_sequence_number=stage_sequence_number,
            stage_gain=1,
            stage_gain_frequency=0,
            input_units=units,
            input_units_description="digital counts",
            output_units=units,
            output_units_description="digital counts",
            symmetry="NONE",
            name=self.name,
            coefficients=self.coeffs,
            decimation_input_sample_rate=sampling_rate,
            decimation_factor=self.decim,
            decimation_offset=self.offset,
            decimation_delay=self.offset / sampling_rate,
            decimation_correction=self.offset / sampling_rate,
        )

    def plot(self, show=False, savefig=None):
        """Plot FIR filter and important parameters"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 6))

        plot_freqz(ax1, self.coeffs, 1, 1 / self.decim)
        plot_phasez(ax3, self.coeffs, 1, 1 / self.decim, self.offset)
        plot_impz(ax2, self.coeffs, 1, delay=self.offset)
        plot_stepz(ax4, self.coeffs, 1, delay=self.offset)

        fig.subplots_adjust(
            hspace=0.4, wspace=0.2, left=0.07, right=0.99, bottom=0.1, top=0.95
        )
        if show:
            plt.show()
        if savefig is not None:
            plt.savefig(savefig)

        return fig


# This script is based on code in a really nice reference:
# http://pyageng.mpastell.com/book/dsp.html
COLOR = "black"


# WCC added "delay" to plot_impz, plot_stepz and plot_responsez
#     fixed freqz y scale from -100 to 5
def plot_freqz(ax, b, a=1, desired_nyquist=None):
    w, h = freqz(b, a)
    h_dB = 20 * np.log10(np.abs(h))

    ax.plot(w / np.pi, h_dB, color=COLOR)

    if desired_nyquist:
        ax.vlines(desired_nyquist, -200, 200, lw=0.5, linestyle="--")
    ax.hlines(
        [0, -48.16, -96.33, -144.5, -192.7], 0.0, 1.0, lw=0.5, linestyle="--"
    )

    ax.set_ylim([-200, 10])

    ax.set_ylabel("Magnitude (db)")
    ax.set_title(r"Frequency response")


def plot_phasez(ax, b, a=1, desired_nyquist=None, delay=0):
    w, h = freqz(b, a)
    h_phase = np.unwrap(np.arctan2(h.imag, h.real))
    h_phase = (h_phase + (w * delay) + np.pi / 2) % np.pi - np.pi / 2  # WCC

    ax.plot(w / np.pi, h_phase, color=COLOR)
    ax.set_ylim([min(min(h_phase), -np.pi), max(max(h_phase), np.pi)])

    if desired_nyquist:
        ylim = ax.get_ylim()
        ax.vlines(desired_nyquist, -200, 200, lw=0.5, linestyle="--")
        ax.set_ylim(*ylim)
    ax.hlines(0, 1.0, 0.0, lw=0.5, linestyle="--")

    ax.set_ylabel("Phase (radians)")
    ax.set_xlabel(r"Normalized Frequency (x$\pi$rad/sample)")


def plot_impz(ax, b, a=1, delay=0.0):
    # FIR filter
    if isinstance(a, int):
        lb = len(b)
    # IIR filter.
    else:
        lb = 100

    impulse = np.repeat(0.0, lb)
    impulse[0] = 1.0
    x = np.linspace(0, lb, lb)

    response = lfilter(b, a, impulse)
    # ax.stem(x, response, linefmt='k-', basefmt='k-', markerfmt='ko')
    ax.plot(x, response, color=COLOR)
    ax.text(
        0.96,
        0.8,
        f"sum={np.sum(response):.2f}\ndelay = {delay:g} samples",
        horizontalalignment="right",
        transform=ax.transAxes,
    )

    ylim = ax.get_ylim()
    # ax.vlines(l / 2.0 , -2, 2, lw=0.5, linestyle="--")
    ax.vlines(delay + 0.5, -2, 2, lw=0.5, linestyle="--")
    ax.set_ylim(*ylim)
    # ax.set_ylim([-.1,1.1])

    ax.set_ylabel("Amplitude")
    ax.set_xlabel("n (samples)")
    ax.set_title("Impulse response")


def plot_stepz(ax, b, a=1, delay=0):
    # FIR filter
    if isinstance(a, int):
        lb = len(b)
    # IIR filter.
    else:
        lb = 100

    impulse = np.repeat(0.0, lb)
    impulse[0] = 1.0
    x = np.linspace(0, lb, lb)

    response = lfilter(b, a, impulse)
    step = np.cumsum(response)

    ax.plot(x, step, color=COLOR)
    # ylim = ax.get_ylim()
    ax.vlines(delay + 0.5, -2, 2, lw=0.5, linestyle="--")
    ax.set_ylim([min(min(step), -0.15), max(max(step), 1.15)])

    ax.set_ylabel("Amplitude")
    ax.set_xlabel(r"n (samples)")
    ax.set_title(r"Step response")
