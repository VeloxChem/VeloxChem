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


class MMGradientDriver(GradientDriver):
    """
    Implements MM gradient driver.

    :param mm_drv:
        The MM driver.

    Instance variables
        - flag: The driver flag.
    """

    def __init__(self, mm_drv):
        """
        Initializes MM gradient driver.
        """

        super().__init__(mm_drv.comm, mm_drv.ostream)

        self.mm_driver = mm_drv
        self.flag = 'MM Gradient Driver'

    def compute(self, molecule):
        """
        Performs calculation of MM analytical gradient.

        :param molecule:
            The molecule.
        """

        self.print_header()

        self.mm_driver.compute(molecule)

        self.gradient = self.mm_driver.get_gradient()
        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())

        self.print_geometry(molecule)
        self.print_gradient(molecule)

        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule):
        """
        Performs calculation of MM energy.

        :param molecule:
            The molecule.
        """

        self.mm_driver.compute(molecule)

        return self.mm_driver.get_energy()
