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
from .openqpscfdriver import OpenQPScfDriver
from .openqpexcitedstatesdriver import OpenQPExcitedStatesDriver


class OpenQPExcitedStateGradientDriver(GradientDriver):
    """
    Implements the OpenQP excited-state (MRSF-/SF-/TD-DFT) gradient driver.

    :param openqp_scf_drv:
        The OpenQP SCF driver providing the OpenQP engine.
    :param openqp_exc_drv:
        Optional OpenQP excited-states driver (for shared settings).

    Instance variables
        - state_deriv_index: OpenQP state/root index of interest.  This matches
          the OpenQP ``properties.grad`` convention: state 1 is the ground
          state (S0), state 2 is S1, state 3 is S2, and so on.
        - tddft_type: OpenQP TD type (`mrsf`, `sf`, `tda`, or `rpa`).
        - nstates: Number of excited-state roots to compute.
    """

    def __init__(self, openqp_scf_drv, openqp_exc_drv=None):
        """
        Initializes the OpenQP excited-state gradient driver.
        """

        errmsg = 'OpenQPExcitedStateGradientDriver: invalid OpenQP SCF driver.'
        assert_msg_critical(isinstance(openqp_scf_drv, OpenQPScfDriver), errmsg)

        super().__init__(openqp_scf_drv.comm, openqp_scf_drv.ostream)

        self.openqp_driver = openqp_scf_drv
        self.exc_driver = (openqp_exc_drv if openqp_exc_drv is not None else
                           OpenQPExcitedStatesDriver(openqp_scf_drv))

        self.state_deriv_index = 3
        self.tddft_type = 'mrsf'
        self.nstates = 10

        self.excited_state_energy = None
        self.total_energy = None

        self.flag = 'OpenQP Excited-State Gradient Driver'

    def set_state_deriv_index(self, state_deriv_index):
        """
        Sets the OpenQP state/root index of interest.

        The index follows the OpenQP ``properties.grad`` convention: state 1 is
        the ground state (S0), 2 is S1, 3 is S2, etc.
        """

        state = int(state_deriv_index)
        assert_msg_critical(
            state > 0,
            'OpenQPExcitedStateGradientDriver: state index must be > 0')
        self.state_deriv_index = state

    def set_tddft_type(self, tddft_type):
        """
        Sets the OpenQP TD type (`mrsf`, `sf`, `tda`, or `rpa`).
        """

        ttype = str(tddft_type).strip().lower()
        assert_msg_critical(
            ttype in ('mrsf', 'sf', 'tda', 'rpa'),
            f'OpenQPExcitedStateGradientDriver: invalid tddft_type '
            f'"{tddft_type}"')
        self.tddft_type = ttype
        self.exc_driver.set_tddft_type(ttype)

    def set_nstates(self, nstates):
        """
        Sets the number of excited-state roots to compute.
        """

        self.nstates = int(nstates)
        self.exc_driver.set_nstates(nstates)

    def update_settings(self, grad_dict, rsp_dict=None, method_dict=None):
        """
        Updates settings in the excited-state gradient driver.

        :param grad_dict:
            The input dictionary of gradient settings.
        :param rsp_dict:
            The input dictionary of response settings.
        :param method_dict:
            The input dictionary of method settings.
        """

        if grad_dict is None:
            grad_dict = {}
        if rsp_dict is None:
            rsp_dict = {}
        if method_dict is None:
            method_dict = {}

        # parse common gradient settings (numerical, do_four_point, delta_h)
        super().update_settings(grad_dict, method_dict={})

        if 'state_deriv_index' in grad_dict:
            state_val = grad_dict['state_deriv_index']
            if isinstance(state_val, (list, tuple, np.ndarray)):
                assert_msg_critical(
                    len(state_val) > 0,
                    'OpenQPExcitedStateGradientDriver: empty state list')
                state_val = state_val[0]
            self.set_state_deriv_index(state_val)

        if 'tddft_type' in grad_dict:
            self.set_tddft_type(grad_dict['tddft_type'])
        if 'tddft_type' in rsp_dict:
            self.set_tddft_type(rsp_dict['tddft_type'])
        if 'method' in rsp_dict:
            self.set_tddft_type(rsp_dict['method'])
        if 'nstates' in rsp_dict:
            self.set_nstates(rsp_dict['nstates'])
        if 'nstate' in rsp_dict:
            self.set_nstates(rsp_dict['nstate'])

    def compute(self, molecule, *args):
        """
        Performs the calculation of the OpenQP excited-state gradient.

        :param molecule:
            The molecule.
        """

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

        self.openqp_driver._invalidate_cache()

        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule, *args):
        """
        Computes the excited-state total energy at the current geometry.

        :param molecule:
            The molecule.

        :return:
            The excited-state total energy.
        """

        if self.rank != mpi_master():
            return None

        mol = self._run_state(molecule, need_gradient=False)
        self.total_energy = float(mol.energies[self.state_deriv_index])
        self.excited_state_energy = float(
            mol.energies[self.state_deriv_index] - mol.energies[1])

        return self.total_energy

    def _compute_analytical_master(self, molecule):
        mol = self._run_state(molecule, need_gradient=True)

        state = self.state_deriv_index
        self.total_energy = float(mol.energies[state])
        self.excited_state_energy = float(mol.energies[state] - mol.energies[1])

        gradient = np.array(mol.grads[state], dtype=float).reshape(
            (molecule.number_of_atoms(), 3))

        return gradient

    def _run_state(self, molecule, need_gradient):
        drv = self.openqp_driver

        scf_type = 'rohf' if self.tddft_type in ('mrsf', 'sf') else \
            drv._effective_scf_type(molecule)
        multiplicity = 3 if self.tddft_type in ('mrsf', 'sf') else \
            drv._effective_multiplicity(molecule)
        functional = drv.functional if drv.method == 'dft' else ''

        nstate = max(int(self.nstates), int(self.state_deriv_index))

        return drv.run_oqp(
            molecule,
            method='tdhf',
            functional=functional,
            scf_type=scf_type,
            multiplicity=multiplicity,
            tdhf_type=self.tddft_type,
            nstate=nstate,
            grad_state=self.state_deriv_index,
            need_gradient=need_gradient,
        )
