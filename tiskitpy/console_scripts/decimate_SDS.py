"""
Script to decimate SDS data, stuff new channels into the SDS structure
and return the modified inventory
"""
import argparse
from pathlib import Path

from obspy import UTCDateTime
from obspy.core.stream import Stream, read
from obspy.core.inventory import read_inventory

from ..decimate import Decimator
from ..util_functions.sds_tools import (sds_channel_files, sds_file_path,
                                        sds_neighbour_streams, is_consecutive)
from ..util_functions.inventory_tools import write_inventory_unique
from ..logger import init_logger, change_console_level

logger = init_logger()


def add_arguments(p):
    """Add decimate_SDS arguments to an ArgumentParser"""
    p.add_argument("SDS_root", help='SDS root directory')
    p.add_argument("inv_file", help="StationXML file")
    p.add_argument("input_sample_rate", type=float,
                   help="Process channels having this sample rate")
    p.add_argument("decim_factors", type=int, nargs="+",
                   choices=[2, 3, 4, 5, 6, 7],
                   help="Sequence of decimation factors to use")
    p.add_argument("--out_file", dest="output_file", default=None,
                   help="Output StationXML filename "
                        "(default = infile.replace('.xml', '_decim.xml')")
    p.add_argument("--out_dir", dest="output_dir", default=None,
                   help="Output data to a separate SDS directory")
    p.add_argument("--station_only", action="store_true", default=False,
                   help="Only create a new StationXML, not new data")
    p.add_argument("--inv_dont_overwrite", action="store_true", default=False,
                   help="Don't overwrite an existing inventory channel")
    p.add_argument("-q", "--quiet", action="store_true",
                   default=False, help="Suppress information messages")


def decimate_SDS(args, inv):
    """
    Decimate SDS channels and add the decimated channels to the inventory

    Each file is decimated using the end of the previous file and the start
    of the next one (if they are consecutive), and the decimation phase is
    carried over from one file to the next, so that the output is the same
    as decimating the whole channel in one piece.  If a trace does not
    follow on from the previous one, its first output sample is the first
    input sample closest to a multiple of the output sampling interval
    (counted from 1970-01-01).

    Only works for broadband channels

    Args:
        args (:class: argparse): Command line arguments
        inv (:class:`obspy.core.inventory.Inventory`): station inventory
    Returns:
        (:class:`obspy.core.inventory.Inventory`): inventory with the
            decimated channels added
    """
    logging_default = 'WARNING' if args.quiet is True else 'INFO'
    change_console_level(logger, logging_default)
    SDS_root = Path(args.SDS_root)
    decimator = Decimator(args.decim_factors)
    D = decimator.decimation_factor
    output_sample_rate = args.input_sample_rate / D
    in_band_code = Decimator.get_band_code('B', args.input_sample_rate)
    out_band_code = Decimator.get_band_code('B', output_sample_rate)
    logger.info('output sampling rate will be {:g} sps, band code will be {}'
                .format(output_sample_rate, out_band_code))
    if in_band_code == out_band_code:
        raise ValueError(f'identical input & output band codes: {in_band_code}')

    out_SDS_root = SDS_root
    if args.output_dir is not None:
        out_SDS_root = Path(args.output_dir)

    sds_channels = sds_channel_files(SDS_root)
    for key in sorted(sds_channels):
        net, sta, loc, cha, typ = key
        seed_id = f'{net}.{sta}.{loc}.{cha}'
        if not cha[0] == in_band_code:
            logger.info(f'{seed_id}: band_code != {in_band_code}, skipping')
            continue
        elif not _channel_in_inv(inv, net, sta, loc, cha):
            logger.warning(f'{seed_id}: not in inventory, skipping')
            continue
        elif not _ch_sample_rate(inv, net, sta, loc, cha) == args.input_sample_rate:
            logger.info(f'{seed_id}: inv sampling rate != '
                        f'{args.input_sample_rate}, skipping')
            continue
        logger.info(f'{seed_id}: processing {len(sds_channels[key])} files')
        och = out_band_code + cha[1:]
        if args.station_only is True:
            stream = read(str(sds_channels[key][0][1]), 'MSEED',
                          headonly=True)
        else:
            stream = None
            last_tr, last_phase = None, None
            for m, prev_st, st, next_st in sds_neighbour_streams(
                    sds_channels[key]):
                d_stream = Stream()
                for tr in sorted(st, key=lambda x: x.stats.starttime):
                    if tr.stats.sampling_rate != args.input_sample_rate:
                        logger.warning(
                            f"{tr.id} {tr.stats.starttime}: sampling rate != "
                            f"{args.input_sample_rate}, skipping...")
                        continue
                    if (last_tr is not None
                            and is_consecutive(last_tr, tr)):
                        phase = (last_phase - last_tr.stats.npts) % D
                    else:
                        phase = _grid_phase(tr, D)
                    change_console_level(logger, 'WARNING')
                    d_stream += decimator.decimate(
                        Stream([tr]), prev_stream=_others(prev_st, st, tr),
                        next_stream=next_st, phase=phase)
                    change_console_level(logger, logging_default)
                    last_tr, last_phase = tr, phase
                if len(d_stream) == 0:
                    continue
                if stream is None:
                    stream = st
                out_fpath = sds_file_path(out_SDS_root, m, chan=och)
                out_fpath.parent.mkdir(parents=True, exist_ok=True)
                d_stream.write(str(out_fpath), 'MSEED')
            if stream is None:
                logger.warning(f'{seed_id}: no data decimated')
                continue
        inv = decimator.update_inventory(
            inv, stream.select(id=seed_id)[:1],
            inv_force_overwrite=args.inv_dont_overwrite is False)
    return inv


def _grid_phase(trace, D):
    """
    Index of the first sample closest to a multiple of D samples from 1970
    """
    n0 = round((trace.stats.starttime - UTCDateTime(0))
               * trace.stats.sampling_rate)
    return int(-n0 % D)


def _others(prev_st, st, tr):
    """prev_st plus the other traces of st (neighbours of tr)"""
    others = Stream([x for x in st if x is not tr])
    return others if prev_st is None else prev_st + others


def _channel_in_inv(inv, net, sta, loc, cha):
    return len(inv.select(network=net, station=sta, location=loc,
                          channel=cha)) > 0


def _ch_sample_rate(inv, net, sta, loc, cha):
    chinfo = inv.select(network=net, station=sta, location=loc,
                        channel=cha)[0][0][0]
    return chinfo.sample_rate


def decimate_SDS_StationXML(args):
    """
    Applies decimate_SDS when the inventory is in a StationXML file

    Args:
        args (:class: argparse): Command line arguments

    """
    inv = read_inventory(args.inv_file, 'STATIONXML')
    inv = decimate_SDS(args, inv)
    of = args.output_file
    if of is None:
        of = args.inv_file.replace('.xml', '_decim.xml')
    of = write_inventory_unique(inv, of)
    logger.info(f'Wrote {of}')


def main(args=None):
    """
    Decimate SDS data and write the updated StationXML file

    Args:
        args (:class: argparse): Command line arguments.  If None, read
            them from the command line
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description="Insert decimated channels and create a new "
                        "StationXML file")
        add_arguments(parser)
        args = parser.parse_args()
    decimate_SDS_StationXML(args)
