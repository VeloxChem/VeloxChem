#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

from os import environ, cpu_count
from pathlib import Path

from .mklconf import configure_mkl_rt


def get_basis_path():
    """
    Returns location of basis files within module.

    :return:
        The location of basis files within module.
    """

    return Path(__file__).parent / "basis"


def set_vlxbasispath():
    """
    Sets location of basis files.
    """

    if 'VLXBASISPATH' not in environ:
        environ['VLXBASISPATH'] = str(get_basis_path())


def get_data_path():
    """
    Returns location of data files within module.

    :return:
        The location of data files within module.
    """

    return Path(__file__).parent / "database"


def set_vlxdatapath():
    """
    Sets location of data files.
    """

    if 'VLXDATAPATH' not in environ:
        environ['VLXDATAPATH'] = str(get_data_path())


def set_omp_num_threads(ncores=None):
    """
    Sets number of OpenMP threads.

    :param ncores:
        The number of cores available for OpenMP threads.
    """

    if ncores is not None:
        environ['OMP_NUM_THREADS'] = str(ncores)
    else:
        if 'OMP_NUM_THREADS' not in environ:
            ncores = cpu_count()
            if ncores is None:
                ncores = 1
            environ['OMP_NUM_THREADS'] = str(ncores)
