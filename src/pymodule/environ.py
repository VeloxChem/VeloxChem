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

from multiprocessing import cpu_count
from os import environ
from pathlib import Path
from sys import stdout


def get_basis_path():
    """Returns location of basis files within module."""
    return Path(__file__).parent / "basis"


def set_vlxbasispath():
    if 'VLXBASISPATH' not in environ:
        environ['VLXBASISPATH'] = str(get_basis_path())


def set_omp_num_threads(ncores=None):
    if ncores is not None:
        environ['OMP_NUM_THREADS'] = ncores
    else:
        if 'OMP_NUM_THREADS' not in environ:
            ncores = cpu_count()
            environ['OMP_NUM_THREADS'] = str(ncores)
            print('* Warning * Environment variable OMP_NUM_THREADS not set.',
                  file=stdout)
            print('* Warning * Setting OMP_NUM_THREADS to {:d}.'.format(ncores),
                  file=stdout)
            stdout.flush()
