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

from mpi4py import MPI
from copy import deepcopy
import time as tm

from .veloxchemlib import XCFunctional, MolecularGrid
from .outputstream import OutputStream
from .gradientdriver import GradientDriver


class ScfGradientDriver(GradientDriver):
    """
    Implements SCF gradient driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - flag: The driver flag.
        - delta_h: The displacement for finite difference.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'SCF Gradient Driver'

        self.numerical = True
        self.delta_h = 0.001

    def compute(self, molecule, ao_basis, scf_drv):
        """
        Performs calculation of gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        """

        start_time = tm.time()
        self.print_header()

        scf_drv.ostream.mute()

        # Currently, only numerical gradients activated
        self.compute_numerical(molecule, ao_basis, scf_drv)

        scf_drv.ostream.unmute()

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule, ao_basis, scf_drv):
        """
        Computes the energy at current geometry.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.

        :return:
            The energy.
        """

        scf_drv.ostream.mute()

        scf_drv.restart = False
        scf_drv.compute(molecule, ao_basis)

        scf_drv.ostream.unmute()

        return scf_drv.get_scf_energy()

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_grad_drv = ScfGradientDriver(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                pass
            elif isinstance(val, XCFunctional):
                new_grad_drv.key = XCFunctional(val)
            elif isinstance(val, MolecularGrid):
                new_grad_drv.key = MolecularGrid(val)
            else:
                new_grad_drv.key = deepcopy(val)

        return new_grad_drv
