#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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


class OpenMMGradientDriver(GradientDriver):
    """
    Implements OpenMM gradient driver.

    :param openmm_drv:
        The OpenMM driver.

    Instance variables
        - flag: The driver flag.
    """

    def __init__(self, openmm_drv):
        """
        Initializes OpenMM gradient driver.
        """

        super().__init__(openmm_drv.comm, openmm_drv.ostream)

        self.openmm_driver = openmm_drv
        self.flag = 'OpenMM Gradient Driver'

    def compute(self, molecule):
        """
        Performs calculation of OpenMM analytical gradient.

        :param molecule:
            The molecule.
        """

        self.print_header()

        self.openmm_driver.compute(molecule)

        self.gradient = self.openmm_driver.get_gradient()
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

        self.openmm_driver.compute(molecule)

        return self.openmm_driver.get_energy()
