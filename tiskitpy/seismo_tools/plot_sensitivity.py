#!/usr/bin/env python3
"""
Plot instrument sensitivity


1: Overall instrument response
2: Instrument sensitivity
"""
import sys
import argparse

from obspy import read_inventory  # , UTCDateTime
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import numpy as np
import matplotlib.pyplot as plt

from ..logger import init_logger

logger = init_logger()



def main():
    arguments = _get_arguments()

    inv = read_inventory(arguments.sta_file)
    net, sta, loc = _get_default(arguments, inv)
    SEED_id = f"{net}.{sta}.{loc}.{arguments.component}"
    try:
        channel = inv.select(
            network=net, station=sta, location=loc, channel=arguments.component
        )[0][0][0]
    except Exception:
        logger.error(f"{SEED_id} not found in {arguments.sta_file}")
    SEED_id = f"{net}.{sta}.{loc}.{channel.code}"
    if not arguments.min_freq:
        logger.info("Calculating minimum frequency from lowest frequency pole:")
        paz = channel.response.get_paz()
        min_pole = np.min(np.abs(paz.poles)) / (2 * np.pi)
        logger.info(f"    LF pole found at {min_pole:.2g} Hz ({1/min_pole:.2g} s)")
        arguments.min_freq = (min_pole) / 10

    outf_base = f"{SEED_id}"
    if arguments.gain_mult:
        outf_base = outf_base + f"_x{arguments.gain_mult:g}"

    plt.figure(1)
    _plot_sensitivity(channel, SEED_id, arguments)
    plt.savefig(outf_base + "_sensitivity.png")

    return 0


def _get_arguments():
    """
    Get command line arguments
    """
    # Set up
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("sta_file", help="StationXML file")
    p.add_argument("component", help="component to plot")
    p.add_argument("max_counts", type=int, help="A/D maximum value (counts)")
    p.add_argument("dyn_range", type=float, help="A/D dynamic range (dB)")
    p.add_argument(
        "--min_counts",
        default=False,
        action="store_true",
        dest="dyn_is_min_counts",
        help="Interpret 4th argument as AD noise level (counts)",
    )
    p.add_argument(
        "-d",
        "--clip_disp",
        type=float,
        default=None,
        help="clip level as displacement (m)",
    )
    p.add_argument(
        "-v",
        "--clip_vel",
        type=float,
        default=None,
        help="clip level as velocity (m/s)",
    )
    p.add_argument(
        "-f",
        "--min_freq",
        type=float,
        default=None,
        help="approx lower frequency limit",
    )
    p.add_argument(
        "-g",
        "--gain_mult",
        type=float,
        default=None,
        help="multiply instrument gain by this value",
    )
    p.add_argument(
        "-o",
        "--output_type",
        default="VEL",
        choices=["VEL", "ACC", "DISP"],
        help="output type",
    )
    p.add_argument("-n", "--network", default="*", help="network code")
    p.add_argument("-s", "--station", default="*", help="station code")
    p.add_argument("-l", "--location", default="*", help="location code")

    # Parse
    results = p.parse_args()
    return results


def _get_default(arguments, inv):
    """Gets first net,sta,loc values from inventory, if not already set"""
    net = arguments.network
    sta = arguments.station
    loc = arguments.location
    if net == "*":
        net = inv[0].code
    if sta == "*":
        sta = inv.select(network=net)[0][0].code
    if loc == "*":
        loc = inv.select(network=net, station=sta)[0][0][0].location_code
    return net, sta, loc


def _get_bradley():
    """Return "HF" vertical noise

    (100m depth, Bradley et al. 1997 JGR)
    """
    f = [0.1, 0.15, 0.3, 1, 6, 10, 15]
    vals = [-122, -122, -103, -118, -146, -135, -132]
    return np.reciprocal(f), vals


def _get_wolin_high():
    """High HF vertical noise model from Wolin and McNamara, 2020 BSSA"""
    f = [10, 15.38, 25, 33.33, 100]
    vals = [-91.5, -88.0, -88.0, -90, -90]
    return np.reciprocal(f), vals


def _get_wolin_lo():
    """Low HF vertical noise model from Wolin and McNamara, 2020 BSSA"""
    f = [1.25, 100]
    vals = [-169.38, -139.43]
    return np.reciprocal(f), vals


def _plot_sensitivity(channel, seed_id, args):
    """Plot instrument's minimum and maximum sensitivity"""

    if args.dyn_is_min_counts:
        min_counts = args.dyn_range
        args.dyn_range = 20 * np.log10(args.max_counts / min_counts)
        logger.info(f"minimum A/D counts provided:  {min_counts:d}")
        logger.info(f"A/D Dynamic range calculated: {args.dyn_range:.0f}dB")
    else:
        min_counts = args.max_counts / np.power(10, args.dyn_range / 20)
        logger.info(f"A/D dynamic range provided:  {args.dyn_range:.0f}dB")
        logger.info(f"min_counts calculated:       {min_counts:.1f}")

    # Plot reference values
    t, vals = get_nlnm()
    plt.semilogx(t, vals, "--", label="NLNM")
    t, vals = get_nhnm()
    plt.semilogx(t, vals, "--", label="NHNM")
    t, vals = _get_bradley()
    plt.semilogx(t, vals, "-.", label="Seafloor noise")
    t, vals = _get_wolin_lo()
    plt.semilogx(t, vals, "r-.", label="Wolin Low HF noise")
    t, vals = _get_wolin_lo()
    plt.semilogx(t, vals, "b-.", label="Wolin High HF noise")

    # Calculate and plot sensitivity
    sr = channel.sample_rate
    t_samp = 1.0 / sr
    nfft = 2 ** (np.ceil(np.log2(sr / args.min_freq)))
    resp = channel.response
    r, f = resp.get_evalresp_response(t_samp, nfft, output="ACC")
    if args.gain_mult:
        r *= args.gain_mult
    plt.semilogx(
        f[1:] ** -1,
        20 * np.log10(min_counts * np.reciprocal(abs(r[1:]))),
        label="Min sensitivity",
    )
    plt.semilogx(
        f[1:] ** -1,
        20 * np.log10(args.max_counts * np.reciprocal(abs(r[1:]))),
        label="Max sensitivity",
    )

    # Plot clip levels, if specified
    if args.clip_vel:
        lev = args.clip_vel * (2 * np.pi * f)
        plt.semilogx(
            f[1:] ** -1, 20 * np.log10(lev[1:]), ":", label="sensor clip level"
        )
    if args.clip_disp:
        lev = args.clip_disp * np.power(2 * np.pi * f, 2)
        plt.semilogx(
            f[1:] ** -1, 20 * np.log10(lev[1:]), ":", label="sensor clip level"
        )

    # Finish up plot
    plt.xlabel("Period (s)")
    plt.ylabel("PSD (dB ref 1 (m/s^2)^2/Hz)")
    # plt.xlim([f[1],f[-1]])
    plt.xlim([0.01, 200])
    plt.ylim([-200, 50])
    title = seed_id
    if channel.equipment:
        if channel.equipment.description:
            title = f"{seed_id} ({channel.equipment.description})"
    if args.gain_mult:
        title = title + f", gain multiplied by {args.gain_mult:g}"
    plt.title(title)
    plt.legend(loc="upper right", fontsize="small", labelspacing=0.1)
    # plt.show()


if __name__ == "__main__":
    sys.exit(main())
