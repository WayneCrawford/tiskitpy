"""
Create an inventory file corresponding to the stations in an SDS directory
"""
import argparse
from pathlib import Path
from dataclasses import dataclass
import sys

from obspy.core.inventory import Inventory
from obspy.clients.fdsn import Client

from ..logger import init_logger, change_console_level

logger = init_logger()


def run(args):
    """
    The main function

    Args:
        args (:class: argparse): Command line arguments
    """
    if args.quiet is True:
        logging_default = 'WARNING'
    else:
        logging_default = 'DEBUG'
    # print(f'{logging_default=}')
    change_console_level(logger, logging_default)
    SDS_root = Path(args.SDS_root)
    inv = Inventory()
    station_channels = _get_station_channels(SDS_root)
    if args.dryrun is True:
        print('Dry run.  Station_channels found are:')
        for s in station_channels:
            print(f'    {s}')
        sys.exit()
    clients = [Client(server) for server in args.servers]
    if args.quiet is False:
        for s, x in zip(args.servers, clients):
            print(f'\nserver "{s}" client:')
            print('\n'.join(x.__str__().splitlines()[:2]))
        print()
    
    for station in station_channels:
        print(station)
        inv_temp = None
        for client in clients:
            try:
                inv_temp = client.get_stations(network=station.network,
                                               station=station.station,
                                               channel=station.channels,
                                               level='response')
            except Exception as e:
                continue  # try next client
            break
        if inv_temp is None:
            logger.warning(f'Station_channels not found on any client')
        else:
            net_exists = False
            for inv_network in inv:
                if inv_temp[0].code == inv_network.code:
                    sta = inv_temp.networks[0].stations[0]
                    # print(sta)
                    inv_network.stations.append(sta)
                    net_exists == True
                    break
            if net_exists is False:
                inv.extend(inv_temp)

    inv.write(args.output_file, format="STATIONXML")


def _get_station_channels(SDS_root):
    station_channels = []
    for year_dir in [x for x in SDS_root.iterdir() if x.is_dir()]:
        for net_dir in [x for x in year_dir.iterdir() if x.is_dir()]:
            for sta_dir in [x for x in net_dir.iterdir() if x.is_dir()]:
                cha_names = [y.name.split('.')[0] for y in
                             [x for x in sta_dir.iterdir() if x.is_dir()]]
                info = StationInfo(net_dir.name, sta_dir.name,
                                   ','.join(cha_names))
                if info not in station_channels:
                    station_channels.append(info)
    return station_channels


@dataclass
class StationInfo:
    """Class for storing station information"""
    network: str
    station: str
    channels: str


def main():
    parser = argparse.ArgumentParser(
        description="Create an inventory file corresponding to the stations in an SDS directory"
    )
    parser.add_argument("SDS_root", help='SDS root directory')
    parser.add_argument("-s", "--servers", dest="servers", nargs='+',
                        default=['IRIS', 'EIDA', 'RESIF'],
                        help="shortcut names for FDSN web service providers "
                             "(see https://docs.obspy.org/packages/obspy.clients.fdsn.html)"
                             " (default: %(default)s)")
    parser.add_argument('--dryrun', action="store_true", help="Just print list of station-channels found")
    parser.add_argument("--of", dest="output_file", default='SDS_directory.xml',
                        help="Output StationXML filename (default = SDS_directory.xml")
    parser.add_argument("-q", "--quiet", action="store_true",
                        default=False, help="Suppress information messages")
    args = parser.parse_args()
    run(args)
