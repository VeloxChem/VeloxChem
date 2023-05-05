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

import time as tm

from .veloxchemlib import mpi_master
from .gradientdriver import GradientDriver
from .errorhandler import assert_msg_critical


class TddftGradientDriver(GradientDriver):
    """
    Implements the TDDFT gradient driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - flag: The driver flag.
        - delta_h: The displacement for finite diference.
        - state_deriv_index: The index of the excited state of interest.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the TDDFT gradient driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'RPA Gradient Driver'
        self.tamm_dancoff = False
        self.state_deriv_index = 1

        self.numerical = True
        self.delta_h = 0.001

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

    def compute(self, molecule, basis, scf_drv, rsp_drv):
        """
        Performs calculation of analytical or numerical gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        :param rsp_drv:
            The RPA or TDA driver.
        """

        if self.tamm_dancoff:
            self.flag = 'TDA Gradient Driver'
        else:
            self.flag = 'RPA Gradient Driver'

        if self.rank == mpi_master():
            error_message = 'TddftGradientDriver: The state of interest is '
            error_message += 'beyond the number of solved states.'
            assert_msg_critical(self.state_deriv_index <= rsp_drv.nstates,
                                error_message)

        self.print_header(self.state_deriv_index)

        start_time = tm.time()

        scf_drv.ostream.mute()

        # Currently, only numerical gradients are available
        self.compute_numerical(molecule, basis, scf_drv, rsp_drv)

        scf_drv.ostream.unmute()

        if self.rank == mpi_master():
            self.print_geometry(molecule)
            self.print_gradient(molecule)

            valstr = '*** Time spent in gradient calculation: '
            valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_energy(self, molecule, basis, scf_drv, rsp_drv):
        """
        Computes the energy at the current position.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param scf_drv:
            The SCF driver.
        :param rsp_drv:
            The RPA or TDA driver.
        """

        scf_drv.restart = False
        scf_results = scf_drv.compute(molecule, basis)
        assert_msg_critical(scf_drv.is_converged,
                            'TddftGradientDriver: SCF did not converge')

        rsp_drv.restart = False
        rsp_results = rsp_drv.compute(molecule, basis, scf_results)
        assert_msg_critical(rsp_drv.is_converged,
                            'TddftGradientDriver: response did not converge')

        if self.rank == mpi_master():
            scf_ene = scf_results['scf_energy']
            exc_ene = rsp_results['eigenvalues'][self.state_deriv_index - 1]
            total_ene = scf_ene + exc_ene
        else:
            total_ene = None
        total_ene = self.comm.bcast(total_ene, root=mpi_master())

        return total_ene
