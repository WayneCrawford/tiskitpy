"""
Script to decimate SDS data, stuff new channels into the SDS structure
and return the modified inventory
"""
from tiskitpy import FIRConverter
from ..logger import init_logger, change_console_level

logger = init_logger()


def main(args):
    """
    The main function

    Just calls the appropriate FIRConverter constructor, then the
    apply_SDS method

    Args:
        args (:class: argparse): Command line arguments
    """
    if args.quiet is True:
        logging_default ='WARNING'
    else:
        logging_default = 'DEBUG'

    if args.zeros_file is not None:
        converter = FIRConverter.from_zeros_file(args.zeros_file)
    elif args.conv_file is not None:
        converter = FIRConverter.from_conv_file(args.conv_file)
    elif args.builtin is not None:
        converter = FIRConverter.from_builtin(args.builtin)
    else:
        logger.error('You must provide one of --zeros_file, --conv_file, or --builtin')
        raise ValueError('Missing FIR filter information')
    
    converter.apply_SDS(args.SDS_input, args.SDS_output, loc_code=args.loc_code)