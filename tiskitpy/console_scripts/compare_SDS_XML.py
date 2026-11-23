"""
List the channels (NSLCs) found in an SDS directory and/or a StationXML file
"""
import argparse

from obspy.core.inventory import read_inventory

from ..util_functions.sds_tools import sds_inventory_comparison


def add_arguments(p):
    """Add compare_SDS_XML arguments to an ArgumentParser"""
    p.add_argument("SDS_root", help='SDS root directory')
    p.add_argument("inv_file", help="StationXML file")


def compare_SDS_XML(SDS_root, inv_file):
    """
    Print the NSLCs in both, only in the StationXML file or only in the SDS

    Args:
        SDS_root (str or Path): SDS root directory
        inv_file (str): StationXML file
    Returns:
        (dict): output of `sds_inventory_comparison()`
    """
    comparison = sds_inventory_comparison(
        SDS_root, read_inventory(inv_file, 'STATIONXML'))
    for key, title in (('both', 'In both'),
                       ('inv_only', 'Only in StationXML file'),
                       ('sds_only', 'Only in SDS directory')):
        print(f'{title} ({len(comparison[key])}):')
        for nslc in comparison[key]:
            print('    ' + '.'.join(nslc))
    return comparison


def main(args=None):
    """
    Compare the NSLCs in an SDS directory and a StationXML file

    Args:
        args (:class: argparse): Command line arguments.  If None, read
            them from the command line
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description='List the NSLCs in both, or only one of, an SDS '
                        'directory and a StationXML file')
        add_arguments(parser)
        args = parser.parse_args()
    compare_SDS_XML(args.SDS_root, args.inv_file)
