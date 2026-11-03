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
from .decimate_SDS import add_arguments as add_decimate_arguments
from .get_SDS_inventory import run as inventory_main
from .fir_convert import main as fir_convert_main
from .fir_convert import add_arguments as add_fir_convert_arguments

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
    add_decimate_arguments(p)
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


    # "FIR_convert_SDS" subcommand
    p = subparsers.add_parser(
        'FIR_convert_SDS',
        description='Convert waveforms in an SDS directory from zero_phase '
                    'to minimum phase response, and add the new channels to '
                    'a StationXML file',
        epilog='You must enter one of --zeros_file, --conv_file, or --builtin. '
               'Only channels in the StationXML file whose response contains '
               'the corresponding FIR are converted.',
        help='Convert SDS waveforms to minimum phase equivalent and update '
             'StationXML')
    add_fir_convert_arguments(p)
    p.set_defaults(func=fir_convert_main)

    args = parser.parse_args()
    if not "func" in args:
        parser.print_help()
    else:
        args.func(args)  # run the appropriate function
    sys.exit(0)


if __name__ == '__main__':
    main()
