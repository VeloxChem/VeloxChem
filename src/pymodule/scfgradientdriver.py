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
import time as tm
import sys

from .gradientdriver import GradientDriver
from .outputstream import OutputStream


class ScfGradientDriver(GradientDriver):
    """
    Implements SCF gradient driver.

    :param scf_drv:
        The SCF driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - delta_h: The displacement for finite difference.
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        super().__init__(comm, ostream)

        self.flag = 'SCF Gradient Driver'
        self.scf_drv = scf_drv
        self.delta_h = 0.001

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Performs calculation of gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        start_time = tm.time()
        self.print_header()

        # atom labels
        labels = molecule.get_labels()

        # Currently, only numerical gradients activated
        self.compute_numerical(molecule, ao_basis, min_basis)

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule, labels)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def pre_numerical_gradient(self):
        """
        Silences the energy driver and save the current ostream state.

        :return:
            The ostream state.
        """

        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False
        return scf_ostream_state

    def compute_energy(self, molecule, ao_basis, min_basis=None):
        """
        Computes the energy at current geometry.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.

        :return:
            The energy.
        """

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        return self.scf_drv.get_scf_energy()

    def post_numerical_gradient(self,
                                ostream_state,
                                molecule,
                                ao_basis,
                                min_basis=None):
        """
        Restores the energy driver to initial state.

        :param ostream_state:
            The state of the energy driver's ostream.
        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        self.scf_drv.ostream.state = ostream_state
