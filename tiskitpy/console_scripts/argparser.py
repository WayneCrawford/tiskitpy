"""
 Main functions for obsinfo-makeStationXML

 Creates obsinfo objects starting with a network object in a hierarchy which
 strongly follows the hierarchy of StationXML files.
 Then converts the objects to a StationXML file using obspy.
"""
import sys
import warnings
from pathlib import Path  # , PurePath
from argparse import ArgumentParser

from .decimate_SDS import main as decimate_main
from .get_SDS_inventory import run as inventory_main
from .fir_convert import main as fir_convert_main

warnings.simplefilter("once")
warnings.filterwarnings("ignore", category=DeprecationWarning)

basedir = Path(__file__).parent.parent


def _print_version(args):
    file = Path(__file__).parent.parent.joinpath("version.py")
    version = {}
    with open(file) as fp:
        exec(fp.read(), version)
    version = version['__version__']

    print(f"{version=}")


def main():
    """
    Entry point all obsinfo sub-commands
    """
    # create the top-level parser
    parser = ArgumentParser(prog="tiskitpy")

    subparsers = parser.add_subparsers(title='subcommands')

    # "version" subcommand
    parser_version = subparsers.add_parser(
        'version', help='Print tiskitpy version',
        description='Print tiskitpy version')
    parser_version.set_defaults(func=_print_version)

    # "decimate_SDS" subcommand
    p = subparsers.add_parser(
        'decimate_SDS',
        description='Insert decimated channels and create a new StationXML file',
        help='decimate SDS data and update StationXML')
    # flags
    p.add_argument("SDS_root", help='SDS root directory')
    p.add_argument("inv_file", help="StationXML file")
    p.add_argument("input_sample_rate", type=float,
                        help="Process channels having this sample rate")
    p.add_argument("decim_factors", type=int, nargs="+",
                        choices=[2,3,4,5,6,7],
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
    p.set_defaults(func=decimate_main)

    # "SDS_inventory" subcommand
    p = subparsers.add_parser(
        'SDS_inventory',
        description='Create an inventory file corresponding to the stations in an SDS directory',
        help='Create an inventory file for an SDS repository')
    # flags
    p.add_argument("SDS_root", help='SDS root directory')
    p.add_argument("-s", "--servers", dest="servers", nargs='+',
                   default=['IRIS', 'EIDA', 'RESIF'],
                   help="shortcut names for FDSN web service providers "
                        "(see https://docs.obspy.org/packages/obspy.clients.fdsn.html)"
                        " (default: %(default)s)")
    p.add_argument('--dryrun', action="store_true", help="Just print list of station-channels found")
    p.add_argument("--of", dest="output_file", default='SDS_directory.xml',
                   help="Set the StationXML filename (default = SDS_directory.xml")
    p.add_argument("-q", "--quiet", action="store_true",
                        default=False, help="Suppress information messages")
    p.set_defaults(func=inventory_main)


    # "FIR_convert" subcommand
    p = subparsers.add_parser(
        'FIR_convert',
        description='Convert waveforms in an SDS directory from zero_phase to minimum phase reponse',
        epilog='You must enter one of --zeros_file, --conv_file, or --builtin. '
               'If the output directory exists, you must specify an unused --loc_code.',
        help='Convert waveforms to minimum phase equivalent')
    # flags
    p.add_argument("SDS_input", help='SDS input directory')
    p.add_argument("SDS_output", help='SDS output directory')
    group = p.add_mutually_exclusive_group()
    group.add_argument('--zeros_file', help='specify a zero-phase filter zeros file')
    group.add_argument('--conv_file', help='specify a zero-phase to minimum-phase conversion file')
    group.add_argument('--builtin', help='specify a built-in conversion')
    p.add_argument("--loc_code", help="Set the location code for output data")
    p.add_argument("-q", "--quiet", action="store_true",
                        default=False, help="Suppress information messages")
    p.set_defaults(func=fir_convert_main)

    args = parser.parse_args()
    if not "func" in args:
        parser.print_help()
    else:
        args.func(args)  # run the appropriate function
    sys.exit(0)


if __name__ == '__main__':
    main()
