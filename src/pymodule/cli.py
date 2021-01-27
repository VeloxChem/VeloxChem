import argparse

from . import __version__


def cli():
    cli = argparse.ArgumentParser(description="Front-end CLI for VeloxChem")
    cli.add_argument("-v", "--version", action="version", version=__version__)
    cli.add_argument("input_file", type=str, help="Input file")
    cli.add_argument(
        "output_file",
        type=str,
        nargs=argparse.REMAINDER,
        default=None,
        help="Output file",
    )

    return cli.parse_args()
