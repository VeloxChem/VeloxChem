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

from .veloxchemlib import mpi_master
from .gradientdriver import GradientDriver


class XTBGradientDriver(GradientDriver):
    """
    Implements XTB gradient driver.

    :param xtb_drv:
        The XTB driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - xtb_drv: The XTB driver.
    """

    def __init__(self, xtb_drv, comm=None, ostream=None):
        """
        Initializes XTB gradient driver.
        """

        super().__init__(xtb_drv, comm, ostream)

        self.flag = 'XTB Gradient Driver'
        self.xtb_drv = xtb_drv

    def compute(self, molecule):
        """
        Performs calculation of XTB analytical gradient.

        :param molecule:
            The molecule.
        """

        self.print_header()

        self.gradient = self.xtb_drv.get_gradient()
        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())

        self.print_geometry(molecule)
        self.print_gradient(molecule) 

        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule):
        """
        Performs calculation of XTB energy.

        :param molecule:
            The molecule.
        """

        self.xtb_drv.compute(molecule, self.ostream)
        return self.xtb_drv.get_energy()
