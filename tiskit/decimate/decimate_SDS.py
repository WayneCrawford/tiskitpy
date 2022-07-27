"""
Script to decimate SDS data, stuff new channels into the SDS structure
and return the modified inventory
"""
import logging
import argparse
from pathlib import Path

from obspy.core.stream import read
from obspy.core.inventory import read_inventory

from tiskit import Decimator


def decimate_SDS(SDS_root, inv, input_sample_rate, decim_list):
    """
    The main function

    Should ensure continuity across day and year boundaries.
    Currently just checks if the requested decimation is an integral divisor
    of the current day's samples and returns an error if not
    Only works for broad-band channels

    Args:
        SDS_root (str or Path): SDS root directory
        inv (:class:`obspy.core.inventory.Inventory`): station inventory
        input_sample_rate (float): sample rate of data to process
        decim_list (list): list of decimation factors (integers between
            2 and 7)
    """
    SDS_root = Path(SDS_root)
    decimator = Decimator(decim_list)
    output_sample_rate = input_sample_rate/decimator.decimation_factor
    in_band_code = Decimator.get_band_code('B', input_sample_rate)
    out_band_code = Decimator.get_band_code('B', output_sample_rate)
    logging.info('output sampling rate will be {:g} sps, band code will be {}'
                 .format(output_sample_rate, out_band_code))
    if in_band_code == out_band_code:
        raise ValueError('identical input & output band codes: {in_band_code}')

    for year_dir in [x for x in SDS_root.iterdir() if x.is_dir()]:
        logging.info('Working on year {str(year_dir.name)}')
        for net_dir in [x for x in year_dir.iterdir() if x.is_dir()]:
            logging.info('\tWorking on net {str(net_dir.name)}')
            for sta_dir in [x for x in net_dir.iterdir() if x.is_dir()]:
                logging.info('\t\tWorking on station {str(sta_dir.name)}')
                for cha_dir in [x for x in sta_dir.iterdir() if x.is_dir()]:
                    logging.info('\t\t\tWorking on channel {}'
                                 .format(str(sta_dir.name)))
                    cha_name = cha_dir.name
                    if not cha_name[0] == in_band_code:
                        continue
                    out_cha_dir = (sta_dir / (out_band_code + cha_name[1:]))
                    out_cha_dir.mkdir()
                    logging.info('\t\t\t\tCreating output channel dir "{}"'
                                 .format(out_cha_dir.name))
                    files = cha_dir.glob(f'*.{cha_name}.*')
                    logging.info('\t\t\t\t{len(files):d} files to process')
                    for f in files:
                        stream = read(str(f), 'MSEED')
                        if stream[0].stats.npts % decimator.decimation_factor == 0:
                            ValueError(
                                f"day's stream length ({stream[0].stats.npts})"
                                " is not divisible by decimator "
                                f"({decimator.decimation_factor})")
                        if stream[0].stats.sampling_rate != input_sample_rate:
                            logging.warning(
                                f"{str(f)} first block's sampling rate != "
                                "{input_sample_rate}, skipping...")
                        net, sta, loc, ich, typ, yr, dy = str(f.name).split('.')
                        d_stream = decimator.decimate(stream)
                        och = out_band_code + ich[1:]
                        outfname = f'{net}.{sta}.{loc}.{och}.{typ}.{yr}.{dy}'
                        d_stream.write(out_cha_dir / outfname, 'MSEED')
                    inv = decimator.update_inventory(inv, stream)
    return inv


def decimate_SDS_StationXML(SDS_root, inv_file, input_sample_rate, decim_list,
                            output_file=None):
    """
    Applies decimate_SDS when the inventory is in a StationXML file

    Args:
        SDS_root (str or Path): SDS root directory
        inv_file (str or Path): Path to StationXML file
        input_sample_rate (float): sample rate of data to process
        decim_list (list): list of decimation factors (integers between
            2 and 7)
        output_file (str): output StationXML filename (None => 
            infile.replace('.xml', '_decim.xml'))

    The output StationXML file will have the suffix  `_decim.xml`
    """
    inv = read_inventory(inv_file, 'STATIONXML')
    inv = decimate_SDS(SDS_root, inv, input_sample_rate, decim_list)
    if output_file is None:
        output_file = inv_file.replace('.xml', '_decim.xml')
    inv.write(output_file, format="STATIONXML")


def main():
    parser = argparse.ArgumentParser(
        description="Insert decimated channels and create new StationXML file"
    )
    parser.add_argument("SDS_root", help='SDS root directory')
    parser.add_argument("inv_file",
                        help="StationXML file")
    parser.add_argument("input_sample_rate", type=float,
                        help="Process only channels having this sample rate")
    parser.add_argument("decim_factor", type=int, nargs="+",
                        choices=[2,3,4,5,6,7],
                        help="Sequence of decimation factors to use)")
    parser.add_argument("--of", dest="output_file", default=None,
                        help="Output StationXML filename "
                             "(default = infile.replace('.xml', '_decim.xml')")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress information messages")
    args = parser.parse_args()
    if not args.quiet:
        logging.disable(logging.DEBUG)
    decimate_SDS_StationXML(args.SDS_root, args.inv_file,
                            args.input_sample_rate, args.decim_factor,
                            args.output_file)
