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

from mpi4py import MPI
import sys

from .veloxchemlib import mpi_master
from .gradientdriver import GradientDriver
from .outputstream import OutputStream


class OpenMMGradientDriver(GradientDriver):
    """
    Implements OpenMM gradient driver.

    :param openmm_drv:
        The OpenMM driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - openmm_drv: The OpenMM driver.
    """

    def __init__(self, openmm_drv, comm=None, ostream=None):
        """
        Initializes OpenMM gradient driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        self.flag = 'OpenMM Gradient Driver'
        self.openmm_drv = openmm_drv

    def compute(self, molecule):
        """
        Performs calculation of OpenMM analytical gradient.

        :param molecule:
            The molecule.
        """

        self.print_header()

        self.gradient = self.openmm_drv.get_gradient()
        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())

        self.print_geometry(molecule)
        self.print_gradient(molecule, molecule.get_labels())

        self.ostream.print_blank()
        self.ostream.flush()
