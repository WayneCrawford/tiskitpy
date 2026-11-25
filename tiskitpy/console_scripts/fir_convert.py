"""
Script to convert SDS data from zero-phase to minimum phase, stuff the new
channels into the SDS structure and return the modified inventory
"""
import argparse
import shutil
from pathlib import Path

from obspy.core.inventory import read_inventory

from tiskitpy import FIRConverter
from ..util_functions.sds_tools import (sds_channel_files, sds_file_path,
                                        sds_inventory_comparison,
                                        log_sds_inventory_comparison)
from ..util_functions.inventory_tools import (find_channel,
                                              write_inventory_unique)
from ..logger import init_logger, change_console_level

logger = init_logger()


def add_arguments(p):
    """Add FIR_convert_SDS arguments to an ArgumentParser"""
    p.add_argument("SDS_root", help='SDS root directory')
    p.add_argument("inv_file", help="StationXML file")
    conv_group = p.add_argument_group('FIR converter (required, choose one)')
    # List it just after the positional arguments, before the options
    p._action_groups.remove(conv_group)
    p._action_groups.insert(1, conv_group)
    conv = conv_group.add_mutually_exclusive_group(required=True)
    conv.add_argument('--zeros_file', help='zero-phase filter zeros file')
    conv.add_argument('--conv_file',
                      help='zero-phase to minimum-phase conversion file')
    builtins = FIRConverter.builtin_names()
    conv.add_argument('--builtin', choices=builtins, metavar='NAME',
                      help='built-in conversion: ' + ', '.join(builtins))
    p.add_argument("--loc_code", default='01',
                   help="Location code for output data (default: %(default)s)")
    p.add_argument("--out_file", dest="output_file", default=None,
                   help="Output StationXML filename "
                        "(default = infile.replace('.xml', '_FIRconv.xml')")
    p.add_argument("--out_dir", dest="output_dir", default=None,
                   help="Output data to a separate SDS directory, which "
                        "must not already exist (see --out_dir_exists_ok)")
    p.add_argument("--no_copy", dest="copy", action="store_false",
                   default=True,
                   help="Don't copy the channels that are not converted "
                        "to --out_dir")
    p.add_argument("--out_dir_exists_ok", action="store_true", default=False,
                   help="Allow --out_dir to be an existing directory.  "
                        "Channels (NSLCs) are added to it, but channels "
                        "that already exist in it are not overwritten")
    p.add_argument("--station_only", action="store_true", default=False,
                   help="Only create a new StationXML, not new data")
    p.add_argument("--inv_dont_overwrite", action="store_true", default=False,
                   help="Don't overwrite an existing inventory channel")
    p.add_argument("-q", "--quiet", action="store_true",
                   default=False, help="Suppress information messages")


def FIR_convert_SDS(args, inv):
    """
    Convert the SDS channels whose response contains the converter's FIR

    Channels (NSLCs) that already exist in the output SDS directory are
    never overwritten: their conversion or copy is aborted.  A summary of
    the action for each SDS channel is logged at the end.

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

    writing_data = args.station_only is False
    writing_out_dir = args.output_dir is not None and writing_data
    if (writing_out_dir and Path(args.output_dir).exists()
            and not args.out_dir_exists_ok):
        raise ValueError(f'--out_dir "{args.output_dir}" already exists '
                         '(use --out_dir_exists_ok to add to it)')
    out_root = Path(args.SDS_root if args.output_dir is None
                    else args.output_dir)
    copying = (args.copy and writing_out_dir
               and out_root.resolve() != Path(args.SDS_root).resolve())

    sds_channels = sds_channel_files(args.SDS_root)
    if not writing_data:
        out_existing = set()
    elif out_root.resolve() == Path(args.SDS_root).resolve():
        out_existing = {key[:4] for key in sds_channels}
    else:
        out_existing = {key[:4] for key in sds_channel_files(out_root)}

    # Sort the SDS channels into those to convert and those to copy/skip
    comparison = sds_inventory_comparison(sds_channels, inv)
    reasons = {x: 'not in StationXML' for x in comparison['sds_only']}
    to_convert = []
    for net, sta, loc, cha in comparison['both']:
        try:
            inv_cha = find_channel(inv, net, sta, loc, cha)
        except ValueError as e:
            reasons[(net, sta, loc, cha)] = str(e)
            continue
        if not converter.matches_channel(inv_cha):
            reasons[(net, sta, loc, cha)] = ('no FIR stage matches the '
                                             'converter')
            continue
        to_convert.append((net, sta, loc, cha))

    # Decide each channel's action, never overwriting an existing NSLC
    actions = {}   # nslc: (action, detail)
    for nslc, reason in reasons.items():
        if not copying:
            actions[nslc] = ('skipped', reason)
        elif nslc in out_existing:
            actions[nslc] = ('aborted', f'{reason}, but already exists in '
                                        f'{out_root}')
        else:
            actions[nslc] = ('copied', reason)
    taken = out_existing | {x for x, a in actions.items()
                            if a[0] == 'copied'}
    channels = []
    for nslc in to_convert:
        target = converter.output_nslc(nslc, args.loc_code)
        out_id = '.'.join(target)
        if target in taken:
            actions[nslc] = ('aborted', f'{out_id} already exists in '
                                        f'{out_root}')
        else:
            actions[nslc] = ('converted', f'to {out_id}')
            taken.add(target)
            channels.append(nslc)
    if len(channels) == 0:
        for nslc, (action, detail) in actions.items():
            if action == 'copied':
                actions[nslc] = ('skipped', detail + ', no channel converted')

    log_sds_inventory_comparison(
        comparison, logger, 'copying' if copying and channels else 'skipping')
    for nslc in sorted(actions):
        seed_id = '.'.join(nslc)
        action, detail = actions[nslc]
        if nslc in comparison['sds_only'] and action != 'aborted':
            continue   # Already in the warning above
        if action == 'converted':
            logger.info(f'{seed_id}: processing')
        elif action == 'copied':
            logger.info(f'{seed_id}: {detail}, copying')
        elif action == 'skipped':
            logger.info(f'{seed_id}: {detail}, skipping')
        else:
            logger.warning(f'{seed_id}: {detail}, aborting')

    if len(channels) == 0:
        logger.warning('No channels to convert' + (', nothing copied'
                                                   if copying else ''))
        _log_summary(actions)
        return inv

    if writing_data:
        converter.apply_SDS(args.SDS_root, args.output_dir,
                            loc_code=args.loc_code, channels=channels,
                            sds_channels=sds_channels,
                            verbose=not args.quiet)
    if copying:
        _copy_channels(sds_channels,
                       [x for x, a in actions.items() if a[0] == 'copied'],
                       out_root)
    _log_summary(actions)
    return converter.update_inventory(
        inv, channels, loc_code=args.loc_code,
        inv_force_overwrite=args.inv_dont_overwrite is False)


def _log_summary(actions):
    """
    Log the action taken for each SDS channel

    Args:
        actions (dict): keys are (net, sta, loc, cha), values are
            (action, detail) tuples
    """
    counts = {a: sum(1 for x in actions.values() if x[0] == a)
              for a in ('converted', 'copied', 'skipped', 'aborted')}
    lines = [f'  {".".join(nslc)}: {actions[nslc][0]} ({actions[nslc][1]})'
             for nslc in sorted(actions)]
    logger.info('Summary: ' + ', '.join(f'{n} {a}' for a, n in counts.items())
                + '\n' + '\n'.join(lines))


def _copy_channels(sds_channels, nslcs, out_root):
    """
    Copy the SDS files of the given channels to another SDS directory

    Args:
        sds_channels (dict): output of `sds_channel_files()`
        nslcs (list): (net, sta, loc, cha) tuples of the channels to copy
        out_root (str or Path): output SDS root directory
    """
    nslcs = set(nslcs)
    for key, files in sds_channels.items():
        if key[:4] not in nslcs:
            continue
        for _, file_path, m in files:
            out_path = sds_file_path(out_root, m)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(file_path, out_path)


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
