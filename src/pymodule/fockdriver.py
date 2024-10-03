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

from mpi4py import MPI
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import _FockDriver
from .outputstream import OutputStream


class FockDriver:
    """
    Implements Fock driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes Fock driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self._fock_drv = _FockDriver()

    def update_block_size_factor(self, nao):
        """
        Updates block size factor for Fock driver.

        :param nao:
            The number of AOs.
        """

        if self.rank == mpi_master():
            if nao < 900:
                block_size_factor = 16
            elif nao < 1800:
                block_size_factor = 8
            elif nao < 3600:
                block_size_factor = 4
            else:
                block_size_factor = 2
        else:
            block_size_factor = None

        block_size_factor = self.comm.bcast(block_size_factor,
                                            root=mpi_master())

        self._fock_drv.set_block_size_factor(block_size_factor)

    def compute(self, screener, *args):

        return self._fock_drv._compute_local_fock(screener, self.rank,
                                                  self.nodes, *args)

    def _compute_fock_omp(self, *args):

        return self._fock_drv._compute_fock_omp(*args)
