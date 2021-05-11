#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
    parser.add_argument('-v', '--version', action='version', version=__version__)
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
