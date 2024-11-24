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
from os import environ
import math
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import _FockDriver
from .veloxchemlib import T4CScreener
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


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

    def compute(self, screener, *args):

        return self._fock_drv._compute_local_fock(screener, self.rank,
                                                  self.nodes, *args)

    def compute_on_subcomm(self, subcomm, screener, *args):

        assert_msg_critical(
            self._check_subcomm(subcomm),
            'FockDriver.compute_on_subcomm: subcomm must be a ' +
            'sub-communicator of self.comm')

        return self._fock_drv._compute_local_fock(screener, subcomm.Get_rank(),
                                                  subcomm.Get_size(), *args)

    def _check_subcomm(self, subcomm):
        """
        Checks that subcomm is a sub-communicator of self.comm.
        """

        translated_ranks = subcomm.allgather(self.comm.Get_rank())

        return all((rank in range(self.nodes)) for rank in translated_ranks)

    def compute_eri(self, molecule, basis, eri_thresh=1.0e-12):
        """
        Computes ERI as a 4D array.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param eri_thresh:
            The ERI screening threshold.

        :return:
            ERI as a 4D array.
        """

        t4c_drv = T4CScreener()
        t4c_drv.partition(basis, molecule, 'eri')

        naos = basis.get_dimensions_of_basis()
        thresh_int = int(-math.log10(eri_thresh))

        return self._fock_drv.compute_eri(t4c_drv, naos, thresh_int)

    def _compute_fock_omp(self, *args):

        return self._fock_drv._compute_fock_omp(*args)

    def _set_block_size_factor(self, factor):

        assert_msg_critical(
            factor in [1, 2, 4, 8, 16, 32, 64, 128],
            'FockDriver._set_block_size_factor: Invalid factor')

        total_cores = self.nodes * int(environ['OMP_NUM_THREADS'])

        if total_cores < 1024:
            extra_factor = 1
        elif total_cores < 4096:
            extra_factor = 2
        else:
            extra_factor = 4

        self._fock_drv._set_block_size_factor(factor * extra_factor)
