import argparse

from . import __version__


def cli():
    cli = argparse.ArgumentParser(description="Front-end CLI for VeloxChem")
    cli.add_argument("-v", "--version", action="version", version=__version__)
    cli.add_argument(
        "input_output_files",
        type=str,
        nargs=argparse.REMAINDER,
        help="Input/Output files",
    )

    return cli.parse_args()
