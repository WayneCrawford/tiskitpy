"""
Script to convert SDS data from zero-phase to minimum phase, stuff the new
channels into the SDS structure and return the modified inventory
"""
import argparse

from obspy.core.inventory import read_inventory

from tiskitpy import FIRConverter
from ..util_functions.sds_tools import sds_channel_files
from ..util_functions.inventory_tools import (find_channel,
                                              write_inventory_unique)
from ..logger import init_logger, change_console_level

logger = init_logger()


def add_arguments(p):
    """Add FIR_convert_SDS arguments to an ArgumentParser"""
    p.add_argument("SDS_root", help='SDS root directory')
    p.add_argument("inv_file", help="StationXML file")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--zeros_file',
                       help='specify a zero-phase filter zeros file')
    group.add_argument('--conv_file',
                       help='specify a zero-phase to minimum-phase conversion file')
    group.add_argument('--builtin', help='specify a built-in conversion')
    p.add_argument("--loc_code", default='01',
                   help="Location code for output data (default: %(default)s)")
    p.add_argument("--out_file", dest="output_file", default=None,
                   help="Output StationXML filename "
                        "(default = infile.replace('.xml', '_FIRconv.xml')")
    p.add_argument("--out_dir", dest="output_dir", default=None,
                   help="Output data to a separate SDS directory")
    p.add_argument("--station_only", action="store_true", default=False,
                   help="Only create a new StationXML, not new data")
    p.add_argument("--inv_dont_overwrite", action="store_true", default=False,
                   help="Don't overwrite an existing inventory channel")
    p.add_argument("-q", "--quiet", action="store_true",
                   default=False, help="Suppress information messages")


def FIR_convert_SDS(args, inv):
    """
    Convert the SDS channels whose response contains the converter's FIR

    Args:
        args (:class: argparse): Command line arguments
        inv (:class:`obspy.core.inventory.Inventory`): station inventory
    Returns:
        (:class:`obspy.core.inventory.Inventory`): inventory with the
            converted channels added
    """
    change_console_level(logger, 'WARNING' if args.quiet else 'INFO')
    if args.zeros_file is not None:
        converter = FIRConverter.from_zeros_file(args.zeros_file)
    elif args.conv_file is not None:
        converter = FIRConverter.from_conv_file(args.conv_file)
    else:
        converter = FIRConverter.from_builtin(args.builtin)

    # Select channels that are in the inventory and contain the FIR
    channels = []
    for net, sta, loc, cha, _ in sorted(sds_channel_files(args.SDS_root)):
        seed_id = f'{net}.{sta}.{loc}.{cha}'
        try:
            inv_cha = find_channel(inv, net, sta, loc, cha)
        except ValueError as e:
            logger.info(f'{seed_id}: {e}, skipping')
            continue
        if not converter.matches_channel(inv_cha):
            logger.info(f'{seed_id}: no FIR stage matches the converter, '
                        'skipping')
            continue
        logger.info(f'{seed_id}: processing')
        channels.append((net, sta, loc, cha))
    if len(channels) == 0:
        logger.warning('No channels to convert')
        return inv

    if args.station_only is False:
        converter.apply_SDS(args.SDS_root, args.output_dir,
                            loc_code=args.loc_code, channels=channels,
                            verbose=not args.quiet)
    return converter.update_inventory(
        inv, channels, loc_code=args.loc_code,
        inv_force_overwrite=args.inv_dont_overwrite is False)


def main(args=None):
    """
    Convert SDS data and write the updated StationXML file

    Args:
        args (:class: argparse): Command line arguments.  If None, read
            them from the command line
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description='Convert waveforms in an SDS directory from '
                        'zero_phase to minimum phase response, and add the '
                        'new channels to a StationXML file')
        add_arguments(parser)
        args = parser.parse_args()
    inv = read_inventory(args.inv_file, 'STATIONXML')
    inv = FIR_convert_SDS(args, inv)
    of = args.output_file
    if of is None:
        of = args.inv_file.replace('.xml', '_FIRconv.xml')
    of = write_inventory_unique(inv, of)
    logger.info(f'Wrote {of}')
