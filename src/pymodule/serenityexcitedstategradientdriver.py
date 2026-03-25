#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np

from .veloxchemlib import mpi_master
from .errorhandler import assert_msg_critical
from .gradientdriver import GradientDriver
from .serenityscfdriver import SerenityScfDriver
from .serenitylrrspeigensolver import SerenityLinearResponseSolver

try:
    from qcserenity import serenipy as spy
except ImportError:
    pass


class SerenityExcitedStateGradientDriver(GradientDriver):
    """
    Implements Serenity excited-state gradient driver.

    :param serenity_scf_drv:
        The Serenity SCF driver.
    :param serenity_rsp_drv:
        Optional Serenity LR response solver.

    Instance variables
        - state_deriv_index: Excited-state index of interest (1-based).
        - exc_method: Excited-state method (`tda` or `tddft`).
    """

    def __init__(self, serenity_scf_drv, serenity_rsp_drv=None):
        """
        Initializes Serenity excited-state gradient driver.
        """

        errmsg = 'SerenityExcitedStateGradientDriver: invalid Serenity SCF '
        errmsg += 'driver.'
        assert_msg_critical(isinstance(serenity_scf_drv, SerenityScfDriver),
                            errmsg)

        super().__init__(serenity_scf_drv.comm, serenity_scf_drv.ostream)

        self.serenity_driver = serenity_scf_drv
        self.rsp_driver = (serenity_rsp_drv if serenity_rsp_drv is not None else
                           SerenityLinearResponseSolver(serenity_scf_drv))

        self.state_deriv_index = 1
        self.exc_method = 'tda'

        self.excited_state_energy = None
        self.total_energy = None

        self.flag = 'Serenity Excited-State Gradient Driver'

    def set_state_deriv_index(self, state_deriv_index):
        """
        Sets the excited-state index of interest (1-based).
        """

        state = int(state_deriv_index)
        assert_msg_critical(
            state > 0,
            'SerenityExcitedStateGradientDriver: state index must be > 0')
        self.state_deriv_index = state

    def set_exc_method(self, exc_method):
        """
        Sets excited-state method (`tda` or `tddft`).
        """

        method = str(exc_method).strip().lower()
        if method in ('tda', 'cis'):
            self.exc_method = 'tda'
        elif method in ('tddft', 'rpa'):
            self.exc_method = 'tddft'
        else:
            errmsg = 'SerenityExcitedStateGradientDriver: invalid exc_method '
            errmsg += f'"{exc_method}"'
            assert_msg_critical(False, errmsg)

        self.rsp_driver.set_exc_method(self.exc_method)

    def update_settings(self,
                        grad_dict,
                        rsp_dict=None,
                        method_dict=None):
        """
        Updates settings in excited-state gradient driver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param rsp_dict:
            The input dictionary of response settings group.
        :param method_dict:
            The input dictionary of method settings group.
        """

        if grad_dict is None:
            grad_dict = {}

        if rsp_dict is None:
            rsp_dict = {}

        if method_dict is None:
            method_dict = {}

        # parse common gradient settings: numerical, do_four_point, delta_h
        super().update_settings(grad_dict, method_dict={})

        if 'state_deriv_index' in grad_dict:
            state_val = grad_dict['state_deriv_index']
            if isinstance(state_val, (list, tuple, np.ndarray)):
                assert_msg_critical(
                    len(state_val) > 0,
                    'SerenityExcitedStateGradientDriver: empty state list')
                state_val = state_val[0]
            self.set_state_deriv_index(state_val)

        if 'exc_method' in grad_dict:
            self.set_exc_method(grad_dict['exc_method'])

        if 'method' in rsp_dict:
            self.set_exc_method(rsp_dict['method'])

        if 'tamm_dancoff' in rsp_dict:
            self.set_exc_method('tda' if bool(rsp_dict['tamm_dancoff']) else
                                'tddft')

        self.rsp_driver.update_settings(rsp_dict, method_dict)

    def compute(self, molecule):
        """
        Performs calculation of Serenity excited-state gradient.

        :param molecule:
            The molecule.
        """

        
        self.ostream.mute()

        if self.numerical:
            if self.rank == mpi_master():
                self.compute_numerical(molecule)
            else:
                self.gradient = None
        else:
            if self.rank == mpi_master():
                self.gradient = self._compute_analytical_master(molecule)
            else:
                self.gradient = None

        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())

        self.print_geometry(molecule)
        self.print_gradient(molecule, [self.state_deriv_index])

        self.serenity_driver._invalidate_cache()
        self.rsp_driver._invalidate_rsp_cache()

        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule):
        """
        Computes excited-state total energy at current geometry.

        :param molecule:
            The molecule.

        :return:
            The excited-state total energy.
        """

        if self.rank != mpi_master():
            return None

        rsp_results = self.rsp_driver.compute(molecule, broadcast=False)

        exc_ene = self._extract_excitation_energy(rsp_results)
        self.excited_state_energy = float(exc_ene)
        self.total_energy = float(self.serenity_driver.get_energy() +
                                  self.excited_state_energy)

        return self.total_energy

    def _compute_analytical_master(self, molecule):

        rsp_results = self.rsp_driver.compute(molecule, broadcast=False)

        self.excited_state_energy = self._extract_excitation_energy(rsp_results)

        # Ensure SCF/system for current geometry is available and synchronized.
        self.serenity_driver._compute_energy_master(molecule)

        mode = self.serenity_driver._current_scf_mode
        if mode == 'restricted':
            grad_task = spy.GradientTask_R(self.serenity_driver._system)
        else:
            grad_task = spy.GradientTask_U(self.serenity_driver._system)

        self._configure_excited_gradient_task(grad_task)

        # with self.serenity_driver._serenity_output_context():
        grad_task.run()

        gradient = np.array(self.serenity_driver._system.getGeometry().getGradients(),
                            dtype=float)
        self.total_energy = float(self.serenity_driver.get_energy() +
                                  self.excited_state_energy)

        return gradient

    def _configure_excited_gradient_task(self, grad_task):
        if hasattr(grad_task, 'generalSettings'):
            grad_task.generalSettings.printLevel = (
                spy.GLOBAL_PRINT_LEVELS.MINIMUM)

        grad_task.settings.gradType = 'analytical'
        grad_task.settings.excMethod = self.exc_method
        grad_task.settings.excGradList = [int(self.state_deriv_index)]

        nstates_req = max(int(self.state_deriv_index),
                          int(self.rsp_driver.nstates))

        grad_task.settings.lrscfSettings.method = self.exc_method
        grad_task.settings.lrscfSettings.nEigen = int(nstates_req)

        if self.rsp_driver.conv_thresh is not None:
            grad_task.settings.lrscfSettings.conv = float(
                self.rsp_driver.conv_thresh)

        if self.rsp_driver.max_cycles is not None:
            grad_task.settings.lrscfSettings.maxCycles = int(
                self.rsp_driver.max_cycles)

        if self.rsp_driver.max_subspace_dimension is not None:
            grad_task.settings.lrscfSettings.maxSubspaceDimension = int(
                self.rsp_driver.max_subspace_dimension)

        if self.rsp_driver.densfit_j is not None:
            grad_task.settings.lrscfSettings.densFitJ = self.rsp_driver.densfit_j

        if self.rsp_driver.grid_accuracy is not None:
            grad_task.settings.lrscfSettings.grid.accuracy = int(
                self.rsp_driver.grid_accuracy)

        if self.rsp_driver.small_grid_accuracy is not None:
            grad_task.settings.lrscfSettings.grid.smallGridAccuracy = int(
                self.rsp_driver.small_grid_accuracy)
        if self.rsp_driver.spinflip:
            grad_task.settings.excMethod = 'sftda'

    def _extract_excitation_energy(self, rsp_results):
        eig = rsp_results['eigenvalues']

        state = int(self.state_deriv_index)
        nst = len(eig)

        errmsg = 'SerenityExcitedStateGradientDriver: requested state index '
        errmsg += f'{state} but only {nst} state(s) are available.'
        assert_msg_critical(state <= nst, errmsg)

        return float(eig[state - 1])
