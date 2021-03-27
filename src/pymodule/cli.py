from mpi4py import MPI
import argparse
import sys
import os

from . import __version__
from .veloxchemlib import mpi_master


def cli():
    """
    Parses argument strings.

    :return:
        The populated namespace.
    """

    cli = argparse.ArgumentParser(description='Front-end CLI for VeloxChem')
    cli.add_argument('-v', '--version', action='version', version=__version__)
    cli.add_argument(
        'input_output_files',
        type=str,
        nargs=argparse.REMAINDER,
        help='Input/Output files',
    )

    return cli.parse_args()


def print_help():
    """
    Prints help text.
    """

    info_txt = [
        '',
        '=================   VeloxChem   =================',
        '',
        'Usage:',
        '    python3 -m veloxchem input_file [output_file]',
        '',
    ]
    if MPI.COMM_WORLD.Get_rank() == mpi_master():
        print(os.linesep.join(info_txt), file=sys.stdout)
    sys.exit(0)
