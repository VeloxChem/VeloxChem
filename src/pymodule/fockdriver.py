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

from .veloxchemlib import mpi_master
from .veloxchemlib import FockDriver


def _FockDriver_mpi_compute(self, comm, screener, density, label,
                            exchange_factor, omega, ithreshold):
    """
    Computes Fock matrix/matrices for given density matrix/matrices.

    :param comm:
        The MPI communicator.
    :param screener:
        The sreener with ERIs data.
    :param density:
        The density matrix/matrices to compute Fock matrix.
    :param label:
        The standard label/labels of Fock matrix type.
    :param exchange_factor:
        The exchange scaling factor.
    :param omega:
        The range separation factor.
    :param ithreshold:
        The threshold of integrals.

    :return:
        Fock matrix/matrices.
    """

    # compute local Fock matrix
    loc_mat = self.compute(screener, comm.Get_rank(), comm.Get_size(), density,
                           label, exchange_factor, omega, ithreshold)

    # reduce Fock matrix
    return loc_mat.reduce(comm, mpi_master())


FockDriver.mpi_compute = _FockDriver_mpi_compute
