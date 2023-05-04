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

from copy import deepcopy
import time as tm
import numpy as np

from .veloxchemlib import mpi_master
from .errorhandler import assert_msg_critical

from .gradientdriver import GradientDriver


class TddftGradientDriver(GradientDriver):
    """
    Implements the TDDFT gradient driver.

    :param scf_drv:
        The SCF driver.
    :param rsp_drv:
        The linear response driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - scf_drv: The SCF driver.
        - rsp_drv: The linear response driver.
        - delta_h: The displacement for finite diference.
        - state_deriv_index: The index of the excited state of interest.
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes the TDDFT gradient driver.
        """

        super().__init__(scf_drv, comm, ostream)

        self.flag = 'RPA Gradient Driver'
        self.tamm_dancoff = False

        self.scf_drv = scf_drv

        self.delta_h = 0.001
        self.state_deriv_index = 1

        self.numerical = True

    def update_settings(self, grad_dict, rsp_dict, method_dict=None):
        """
        Updates settings in gradient driver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param rsp_dict:
            The input dictionary of response settings  group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        if method_dict is None:
            method_dict = {}

        # basic settings from parent class
        super().update_settings(grad_dict, method_dict)

        if 'tamm_dancoff' in rsp_dict:
            key = rsp_dict['tamm_dancoff'].lower()
            self.tamm_dancoff = True if key in ['yes', 'y'] else False

        if self.tamm_dancoff:
            self.flag = 'TDA Gradient Driver'
        else:
            self.flag = 'RPA Gradient Driver'

        # Excited state of interest
        # NOTE: the indexing starts at 1.
        if 'state_deriv_index' in grad_dict:
            self.state_deriv_index = int(grad_dict['state_deriv_index'])

    def compute(self, molecule, basis, rsp_drv):
        """
        Performs calculation of analytical or numerical gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver.
        """

        # sanity check
        # TODO automatically determine nstates for response driver
        if self.rank == mpi_master():
            error_message = 'TddftGradientDriver: some of the '
            error_message += 'selected states have not been calculated.'
            assert_msg_critical(self.state_deriv_index <= rsp_drv.nstates,
                                error_message)

        self.print_header(self.state_deriv_index)

        start_time = tm.time()

        self.scf_drv.ostream.mute()

        # Currently, only numerical gradients are available
        self.compute_numerical(molecule, basis, rsp_drv)

        self.scf_drv.ostream.unmute()

        # print gradient
        if self.rank == mpi_master():
            self.print_geometry(molecule)
            self.print_gradient(molecule)

            valstr = '*** Time spent in gradient calculation: '
            valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_energy(self, molecule, basis, rsp_drv):
        """
        Computes the energy at the current position.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param rsp_drv:
            The linear response driver.
        :param min_basis:
            The reduced basis set.
        """

        scf_results = self.scf_drv.compute(molecule, basis)

        rsp_drv._is_converged = False
        rsp_results = rsp_drv.compute(molecule, basis, scf_results)

        if self.rank == mpi_master():
            scf_ene = scf_results['scf_energy']
            exc_ene = rsp_results['eigenvalues'][self.state_deriv_index - 1]
            total_ene = scf_ene + exc_ene
        else:
            total_ene = None
        total_ene = self.comm.bcast(total_ene, root=mpi_master())

        return total_ene
