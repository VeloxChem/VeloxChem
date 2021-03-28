#
#                           VELOXCHEM 1.0-RC
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
