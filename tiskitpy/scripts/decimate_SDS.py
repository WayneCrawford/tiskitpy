"""
Script to decimate SDS data, stuff new channels into the SDS structure
and return the modified inventory
"""
import argparse
from pathlib import Path

from obspy.core.stream import read
from obspy.core.inventory import read_inventory

from ..decimate import Decimator
from ..logger import init_logger, change_console_level

logger = init_logger()


def decimate_SDS(args, inv):
    """
    The main function

    Should ensure continuity across day and year boundaries.
    Currently just checks if the requested decimation is an integral divisor
    of the current day's samples and returns an error if not
    Only works for broadband channels

    Args:
        args (:class: argparse): Command line arguments
        inv (:class:`obspy.core.inventory.Inventory`): station inventory
    """
    if args.quiet is True:
        logging_default ='WARNING'
    else:
        logging_default = 'DEBUG'
    change_console_level(logger, logging_default)
    print(f'{logging_default=}')
    SDS_root = Path(args.SDS_root)
    decimator = Decimator(args.decim_factors)
    output_sample_rate = args.input_sample_rate/decimator.decimation_factor
    in_band_code = Decimator.get_band_code('B', args.input_sample_rate)
    out_band_code = Decimator.get_band_code('B', output_sample_rate)
    logger.info('output sampling rate will be {:g} sps, band code will be {}'
                 .format(output_sample_rate, out_band_code))
    if in_band_code == out_band_code:
        raise ValueError(f'identical input & output band codes: {in_band_code}')

    out_SDS_root = SDS_root
    if args.output_dir is not None:
        out_SDS_root = Path(args.output_dir)
        
    for year_dir in [x for x in SDS_root.iterdir() if x.is_dir()]:
        year_name = year_dir.name
        logger.info(f' {year_name=}')
        for net_dir in [x for x in year_dir.iterdir() if x.is_dir()]:
            net_name = net_dir.name
            logger.info(f'\t{net_name=}')
            for sta_dir in [x for x in net_dir.iterdir() if x.is_dir()]:
                sta_name=sta_dir.name
                logger.info(f'\t\t{sta_name=}')
                for cha_dir in [x for x in sta_dir.iterdir() if x.is_dir()]:
                    cha_name = cha_dir.name
                    if not cha_name[0] == in_band_code:
                        logger.info(f'\t\t\t{cha_name=}: band_code != {in_band_code}, skipping')
                        continue
                    elif not _channel_in_inv(inv, net_name, sta_name, cha_name.split('.')[0]):
                        logger.warning(f'\t\t\t{cha_name=}: not in inventory, skipping')
                        continue
                    elif not _ch_sample_rate(inv, net_name, sta_name, cha_name.split('.')[0]) == args.input_sample_rate:
                        logger.info(f'\t\t\t{cha_name=}: inv sampling rate != {args.input_sample_rate}, skipping')
                        continue
                    else:
                        logger.info(f'\t\t\t{cha_name=}: processing')
                    if args.station_only is True:
                        stream = read(str(list(cha_dir.glob(f'*.{cha_name}.*'))[0]), 'MSEED')
                    else:
                        out_cha_dir = out_SDS_root / year_name / net_name / sta_name / (out_band_code + cha_name[1:])
                        out_cha_dir.mkdir(parents=True)
                        logger.info('\t\t\t\tCreating output channel dir "{}"'
                                    .format(out_cha_dir.name))
                        files = list(cha_dir.glob(f'*.{cha_name}.*'))
                        logger.info(f'\t\t\t\t{len(files):d} files to process')
                        for f in files:
                            stream = read(str(f), 'MSEED')
                            if stream[0].stats.npts % decimator.decimation_factor == 0:
                                ValueError(
                                    f"day's stream length ({stream[0].stats.npts})"
                                    " is not divisible by decimator "
                                    f"({decimator.decimation_factor})")
                            if stream[0].stats.sampling_rate != args.input_sample_rate:
                                logger.warning(
                                    f"{str(f)} first block's sampling rate != "
                                    f"{args.input_sample_rate}, skipping...")
                            net, sta, loc, ich, typ, yr, dy = str(f.name).split('.')
                            # logging.getLogger().setLevel(logging.WARN)
                            change_console_level(logger, 'WARNING')
                            d_stream = decimator.decimate(stream)
                            change_console_level(logger, logging_default)
                            # logger.setLevel(logging_default)
                            och = out_band_code + ich[1:]
                            outfname = f'{net}.{sta}.{loc}.{och}.{typ}.{yr}.{dy}'
                            d_stream.write(out_cha_dir / outfname, 'MSEED')
                    inv = decimator.update_inventory(
                        inv, stream,
                        inv_force_overwrite=args.inv_dont_overwrite is False)
    return inv

def _channel_in_inv(inv, net, sta, cha):
    return len(inv.select(network=net, station=sta, channel=cha))> 0

def _ch_sample_rate(inv, net, sta, cha):
    chinfo = inv.select(network=net, station=sta, channel=cha)[0][0][0]
    return chinfo.sample_rate


def decimate_SDS_StationXML(args):
    """
    Applies decimate_SDS when the inventory is in a StationXML file

    Args:
        args (:class: argparse): Command line arguments

    """
    inv = read_inventory(args.inv_file, 'STATIONXML')
    inv = decimate_SDS(args, inv)
    if args.output_file is None:
        of = args.inv_file.replace('.xml', '_decim.xml')
    else:
        of = args.output_file
    startfile = of
    counter = 1
    while Path(of).exists():
        of = startfile.replace('.xml', str(counter)+'.xml')
        counter +=1
    inv.write(of, format="STATIONXML")


def main():
    parser = argparse.ArgumentParser(
        description="Insert decimated channels and create a new StationXML file"
    )
    parser.add_argument("SDS_root", help='SDS root directory')
    parser.add_argument("inv_file", help="StationXML file")
    parser.add_argument("input_sample_rate", type=float,
                        help="Process channels having this sample rate")
    parser.add_argument("decim_factors", type=int, nargs="+",
                        choices=[2,3,4,5,6,7],
                        help="Sequence of decimation factors to use")
    parser.add_argument("--out_file", dest="output_file", default=None,
                        help="Output StationXML filename "
                             "(default = infile.replace('.xml', '_decim.xml')")
    parser.add_argument("--out_dir", dest="output_dir", default=None,
                        help="Output data to a separate SDS directory")
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument("--out_dir", dest="output_dir", default=None,
    #                    help="Output data to a separate SDS directory. Triggers "
    #                          "--out_only_decim")
    # group.add_argument("--out_only_decim", action="store_true",  default=False,
    #                    help="Only output decimated channels to the new StationXML file")
    parser.add_argument("--station_only", action="store_true", default=False,
                        help="Only create a new StationXML, not new data")
    parser.add_argument("--inv_dont_overwrite", action="store_true", default=False,
                        help="Don't overwrite an existing inventory channel")
    parser.add_argument("-q", "--quiet", action="store_true",
                        default=False, help="Suppress information messages")
    args = parser.parse_args()
    # if args.out_dir is True:
    #      args.out_only_decim = True
    decimate_SDS_StationXML(args)
