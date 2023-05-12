import argparse

from . import __version__

def cli():
    """
    Generate command-line argument parser.

    :return:
        The parser.
    """

    usage = """
=================   VeloxChem   =================

%(prog)s input_file [output_file]

or

python -m veloxchem input_file [output_file]
    """
    parser = argparse.ArgumentParser(prog="vlx", usage=usage)
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=__version__,
    )
    parser.add_argument(
        'input_file',
        type=str,
        help='Input file',
    )
    parser.add_argument(
        'output_file',
        type=str,
        nargs="?",
        default=".",
        help='Output file (default: STDOUT)',
    )

    return parser
