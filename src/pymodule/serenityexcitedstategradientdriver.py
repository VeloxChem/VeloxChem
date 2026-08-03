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

from copy import deepcopy

import numpy as np

from .veloxchemlib import mpi_master, hartree_in_ev

from .errorhandler import assert_msg_critical
from .gradientdriver import GradientDriver
from .serenityscfdriver import SerenityScfDriver
from .serenitylrrspeigensolver import SerenityLinearResponseSolver
from .transitiondensitytracker import (StateTrackingError,
                                       StateTrackingResult,
                                       StateTrackingStatus,
                                       TransitionDensityTracker)

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
        self.enforce_same_multiplicity = True
        self.target_multiplicity = None
        self.s2_tolerance = 0.5

        self.excited_state_energy = None
        self.selected_excitation_energy = None
        self.total_energy = None
        self.computes_energy_with_gradient = True

        self.reference_s2 = None
        self.delta_s2 = None
        self.state_s2 = None
        self.state_multiplicities = None
        self.ideal_state_s2 = None
        self.s2_deviation = None
        self.multiplicity_valid_mask = None
        self.selected_s2 = None
        self.selected_multiplicity = None
        self.selected_s2_deviation = None

        self.flag = 'Serenity Excited-State Gradient Driver'

        self._grad_task = None
        self.state_tracker = None
        self.tracking_info = None
        self.last_lrscf_controller = None
        self._last_tracking_system = None
        self._last_tracking_mode = None
        self._last_tracking_molecule = None
        self._tracking_applied_in_compute = False
        # OptimizationEngine enables this while bounded LOW_OVERLAP recovery
        # is still available.  Non-strict fallbacks are applied only after
        # those retries have been exhausted.
        self._tracking_recovery_active = False
        self._current_gradient_task_roots = []
        self.tracking_history = []
        self._tracking_evaluation_counter = 0

        self._input_keywords['gradient'].update({
            'enforce_same_multiplicity':
                ('bool', 'restrict spin-flip root tracking by multiplicity'),
            'target_multiplicity':
                ('int', 'target spin multiplicity; omit to infer initially'),
            's2_tolerance':
                ('float', 'maximum deviation from ideal spin squared'),
        })


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

        assert_msg_critical(
            self.target_multiplicity is None or self.target_multiplicity > 0,
            'SerenityExcitedStateGradientDriver: target_multiplicity must be '
            'positive.')
        assert_msg_critical(
            self.s2_tolerance >= 0.0,
            'SerenityExcitedStateGradientDriver: s2_tolerance must be '
            'non-negative.')

    def compute(self, molecule):
        """
        Performs calculation of Serenity excited-state gradient.

        :param molecule:
            The molecule.
        """

        self._tracking_applied_in_compute = False
        self._current_gradient_task_roots = []

        tracking_error = None
        try:
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
        except StateTrackingError as error:
            tracking_error = (str(error), error.result)

        tracking_error = self.comm.bcast(
            tracking_error, root=mpi_master())
        if tracking_error is not None:
            message, result = tracking_error
            raise StateTrackingError(message, result)

        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())
        if self.rank == mpi_master():
            state_payload = {
                'state_deriv_index': int(self.state_deriv_index),
                'target_multiplicity': self.target_multiplicity,
                'excited_state_energy': self.excited_state_energy,
                'selected_excitation_energy':
                    self.selected_excitation_energy,
                'total_energy': self.total_energy,
                'reference_s2': self.reference_s2,
                'delta_s2': self.delta_s2,
                'state_s2': self.state_s2,
                'state_multiplicities': self.state_multiplicities,
                'ideal_state_s2': self.ideal_state_s2,
                's2_deviation': self.s2_deviation,
                'multiplicity_valid_mask':
                    self.multiplicity_valid_mask,
                'selected_s2': self.selected_s2,
                'selected_multiplicity': self.selected_multiplicity,
                'selected_s2_deviation': self.selected_s2_deviation,
                'tracking_info': self.tracking_info,
                'tracking_applied': self._tracking_applied_in_compute,
            }
        else:
            state_payload = None

        state_payload = self.comm.bcast(
            state_payload, root=mpi_master())
        self.state_deriv_index = int(state_payload['state_deriv_index'])
        self.target_multiplicity = state_payload['target_multiplicity']
        self.excited_state_energy = state_payload['excited_state_energy']
        self.selected_excitation_energy = (
            state_payload['selected_excitation_energy'])
        self.total_energy = state_payload['total_energy']
        self.reference_s2 = state_payload['reference_s2']
        self.delta_s2 = state_payload['delta_s2']
        self.state_s2 = state_payload['state_s2']
        self.state_multiplicities = (
            state_payload['state_multiplicities'])
        self.ideal_state_s2 = state_payload['ideal_state_s2']
        self.s2_deviation = state_payload['s2_deviation']
        self.multiplicity_valid_mask = (
            state_payload['multiplicity_valid_mask'])
        self.selected_s2 = state_payload['selected_s2']
        self.selected_multiplicity = state_payload['selected_multiplicity']
        self.selected_s2_deviation = (
            state_payload['selected_s2_deviation'])
        self.tracking_info = state_payload['tracking_info']
        self._tracking_applied_in_compute = bool(
            state_payload['tracking_applied'])
        if (self._tracking_applied_in_compute and
                self.tracking_info is not None):
            self.tracking_history.append(deepcopy(self.tracking_info))

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
        eigenvalues = np.asarray(
            rsp_results['eigenvalues'], dtype=float).reshape(-1)

        if self.rsp_driver.spinflip:
            self._set_spin_metadata(rsp_results)
            if self.enforce_same_multiplicity:
                self.state_deriv_index = (
                    self._select_energy_root_by_multiplicity(
                        self.state_deriv_index))
            self._update_selected_spin_diagnostics()

        self.excited_state_energy = eigenvalues
        self.selected_excitation_energy = self._extract_excitation_energy(
            rsp_results)
        self.total_energy = float(
            self.serenity_driver.get_energy() +
            self.selected_excitation_energy)

        return self.total_energy

    def _compute_analytical_master(self, molecule):
        # Ensure SCF/system for current geometry is available and synchronized.
        self.serenity_driver._compute_energy_master(molecule)

        mode = self.serenity_driver._current_scf_mode
        requested_state = int(self.state_deriv_index)
        grad_task = self._run_excited_gradient_task(mode)

        if self.rsp_driver.spinflip:
            self._update_spin_metadata(
                grad_task.getLRSCFController())

            if self.enforce_same_multiplicity:
                selected_state = self._select_same_multiplicity_state(
                    molecule, grad_task, mode, requested_state)
                if selected_state != requested_state:
                    self.state_deriv_index = selected_state
                    grad_task = self._run_excited_gradient_task(mode)
                    self._update_spin_metadata(
                        grad_task.getLRSCFController())
            else:
                selected_state = requested_state
        else:
            selected_state = requested_state

        self.state_deriv_index = int(selected_state)
        if self.rsp_driver.spinflip:
            self._update_selected_spin_diagnostics()
            self.ostream.print_info(
                'Spin-flip state selection: '
                f'root {self.state_deriv_index}, '
                f'{self._format_selected_spin_label()}, '
                f'<S^2> = {self.selected_s2:.6f}, '
                f'deviation = {self.selected_s2_deviation:.6f}.')
            self.ostream.flush()

        controller = grad_task.getLRSCFController()
        excitation_energies = np.asarray(
            controller.getExcitationEnergies('isolated'),
            dtype=float).reshape(-1) / hartree_in_ev()

        assert_msg_critical(
            self.state_deriv_index <= excitation_energies.size,
            'SerenityExcitedStateGradientDriver: selected state index '
            f'{self.state_deriv_index} but only {excitation_energies.size} '
            'state(s) are available.')

        self.excited_state_energy = excitation_energies
        self.selected_excitation_energy = float(
            excitation_energies[self.state_deriv_index - 1])
        self.total_energy = float(
            self.serenity_driver.get_energy() +
            self.selected_excitation_energy)
        gradient = np.array(
            self.serenity_driver._system.getGeometry().getGradients(),
            dtype=float)

        self._grad_task = grad_task
        self.last_lrscf_controller = controller
        self._last_tracking_system = self.serenity_driver._system
        self._last_tracking_mode = mode
        self._last_tracking_molecule = molecule

        if (self.rsp_driver.spinflip and
                self.enforce_same_multiplicity and
                self.state_tracker is not None):
            if (self.tracking_info is not None and
                    self.tracking_info.get(
                        'reference_update_recommended', False) and
                    not self.tracking_info.get('initialized', False)):
                # The selected-root gradient and its controller snapshot are
                # staged together. geomeTRIC's accepted-step hook commits both;
                # a rejected trial discards the proposal.
                self.state_tracker.propose_reference(
                    self._last_tracking_system,
                    self.last_lrscf_controller,
                    self._last_tracking_mode,
                    self.state_deriv_index,
                )
                self.tracking_info['reference_staged'] = True
            self._tracking_applied_in_compute = True

        if self.tracking_info is not None:
            self.tracking_info.update({
                'step': int(len(self.tracking_history)),
                'evaluation_id': int(self._tracking_evaluation_counter),
                'gradient_root': int(self.state_deriv_index),
                'gradient_task_roots':
                    list(self._current_gradient_task_roots),
                'gradient_recomputed':
                    len(self._current_gradient_task_roots) > 1,
                'selected_excitation_energy':
                    float(self.selected_excitation_energy),
                'total_energy': float(self.total_energy),
                'gradient_rms': float(np.sqrt(np.mean(gradient**2))),
                'state_excitation_energies':
                    np.asarray(self.excited_state_energy, dtype=float).copy(),
                'state_multiplicities':
                    np.asarray(self.state_multiplicities, dtype=int).copy(),
                'state_s2': np.asarray(self.state_s2, dtype=float).copy(),
                's2_deviation':
                    np.asarray(self.s2_deviation, dtype=float).copy(),
            })
            self._tracking_evaluation_counter += 1
            self._print_tracking_diagnostics(self.tracking_info)

        return gradient

    def _run_excited_gradient_task(self, mode):
        """
        Runs one Serenity gradient task for the current raw response root.
        """

        self._current_gradient_task_roots.append(
            int(self.state_deriv_index))

        if mode == 'restricted':
            grad_task = spy.GradientTask_R(self.serenity_driver._system)
        else:
            grad_task = spy.GradientTask_U(self.serenity_driver._system)

        self._configure_excited_gradient_task(grad_task)

        with self.serenity_driver._serenity_output_context():
            grad_task.run()

        return grad_task

    def _update_spin_metadata(self, controller):
        metadata = self.rsp_driver.get_spinflip_metadata(controller)
        self._set_spin_metadata(metadata)

    def _set_spin_metadata(self, metadata):
        self.reference_s2 = metadata['reference_s2']
        self.delta_s2 = np.asarray(metadata['delta_s2'], dtype=float)
        self.state_s2 = np.asarray(metadata['state_s2'], dtype=float)
        self.state_multiplicities = np.asarray(
            metadata['state_multiplicities'], dtype=int)
        self.ideal_state_s2 = np.asarray(
            metadata['ideal_state_s2'], dtype=float)
        self.s2_deviation = np.asarray(
            metadata['s2_deviation'], dtype=float)

    def _update_selected_spin_diagnostics(self):
        selected_index = int(self.state_deriv_index) - 1
        assert_msg_critical(
            0 <= selected_index < self.state_s2.size,
            'SerenityExcitedStateGradientDriver: selected spin-flip root is '
            'outside the spin metadata range.')

        self.selected_s2 = float(self.state_s2[selected_index])
        self.selected_multiplicity = int(
            self.state_multiplicities[selected_index])
        self.selected_s2_deviation = float(
            self.s2_deviation[selected_index])

    def _format_selected_spin_label(self):
        """
        Formats the spin assignment of the selected spin-flip root.

        Nearest-multiplicity classification always returns some multiplicity,
        so an SF root with ``<S^2>`` near 1 -- equidistant from the singlet
        (0) and triplet (2) values -- is reported as a clean singlet or
        triplet depending only on which side of the tie it falls.  Those roots
        are spin-incomplete rather than spin-pure, and printing a bare
        "multiplicity 3" for ``<S^2> = 1.02`` claims a triplet that is not
        there.  Qualify the label whenever the root sits further from its
        nearest physical ``<S^2>`` than the configured tolerance allows.
        """

        if self.selected_s2_deviation > self.s2_tolerance:
            return ('spin-incomplete (nearest multiplicity '
                    f'{self.selected_multiplicity})')

        return f'multiplicity {self.selected_multiplicity}'

    def _select_energy_root_by_multiplicity(self, requested_state):
        """
        Validates a raw root and records spin diagnostics for energy-only use.

        Serenity SF-TDA roots are not guaranteed spin-pure, so the approximate
        ``<S^2>`` classification is diagnostic here rather than a hard mask.
        """

        requested_state = int(requested_state)
        nstates = int(self.state_s2.size)
        assert_msg_critical(
            1 <= requested_state <= nstates,
            'SerenityExcitedStateGradientDriver: requested state index is '
            'outside the spin-flip spectrum.')

        if self.target_multiplicity is None:
            self.target_multiplicity = int(
                self.state_multiplicities[requested_state - 1])

        self.multiplicity_valid_mask = (
            (self.state_multiplicities == int(self.target_multiplicity)) &
            (self.s2_deviation <= float(self.s2_tolerance))
        )
        return requested_state

    def _select_same_multiplicity_state(self, molecule, grad_task, mode,
                                        requested_state):
        """
        Tracks all valid roots, then applies the configured failure policy.

        Approximate ``<S^2>`` never removes a finite response root before the
        overlap comparison.  A multiplicity mismatch is reported separately as
        ``SPIN_CONFLICT``.
        """

        del molecule

        nstates = int(np.asarray(self.state_s2).size)
        assert_msg_critical(
            1 <= requested_state <= nstates,
            'SerenityExcitedStateGradientDriver: requested state index is '
            'outside the spin-flip spectrum.')

        if self.target_multiplicity is None:
            self.target_multiplicity = int(
                self.state_multiplicities[requested_state - 1])

        same_inferred_multiplicity = (
            (np.asarray(self.state_multiplicities, dtype=int) ==
             int(self.target_multiplicity))
        )
        within_s2_tolerance = (
            np.asarray(self.s2_deviation, dtype=float) <=
            float(self.s2_tolerance)
        )
        self.multiplicity_valid_mask = (
            same_inferred_multiplicity & within_s2_tolerance)

        if self.state_tracker is None:
            self.state_tracker = TransitionDensityTracker(
                self.rsp_driver,
                target_state=requested_state,
            )

        controller = grad_task.getLRSCFController()
        # Spin compatibility must combine the inferred multiplicity with the
        # <S^2> tolerance.  Nearest-multiplicity classification alone labels a
        # spin-incomplete SF root (<S^2> ~ 1, equidistant from 0 and 2) as a
        # singlet, so following it would be reported as CONFIDENT instead of
        # SPIN_CONFLICT.
        candidate_metadata = self._tracking_candidate_metadata(
            controller, self.multiplicity_valid_mask)
        result = self.state_tracker.track(
            self.serenity_driver._system,
            controller,
            mode,
            active_reference_state=requested_state,
            # Only invalid response vectors are hard-excluded by the tracker.
            allowed_states=None,
            candidate_metadata=candidate_metadata,
        )

        if (result['initialized'] and
                not bool(self.multiplicity_valid_mask[requested_state - 1])):
            diagnostics = result.to_dict()
            diagnostics.pop('status', None)
            diagnostics['assignment_confident'] = False
            diagnostics['warnings'] = list(diagnostics['warnings']) + [
                'initial raw root conflicts with target multiplicity'
            ]
            result = StateTrackingResult(
                StateTrackingStatus.SPIN_CONFLICT, diagnostics)

        selected_state, result = self._apply_tracking_failure_policy(result)
        selected_index = selected_state - 1
        info = result.to_dict()
        info.update({
            'target_multiplicity': int(self.target_multiplicity),
            'selected_state': int(selected_state),
            'selected_multiplicity':
                int(self.state_multiplicities[selected_index]),
            'selected_s2': float(self.state_s2[selected_index]),
            'selected_s2_deviation':
                float(self.s2_deviation[selected_index]),
            'spin_diagnostic_roots': [
                int(state) for state in
                np.flatnonzero(same_inferred_multiplicity) + 1
            ],
            's2_tolerance_roots': [
                int(state) for state in
                np.flatnonzero(self.multiplicity_valid_mask) + 1
            ],
            'reference_updated': False,
            'reference_staged': False,
        })
        self.tracking_info = info

        return int(selected_state)

    def _tracking_candidate_metadata(self, controller, spin_compatible):
        """Builds per-root diagnostics, including optional solver metadata."""

        nstates = int(np.asarray(self.state_s2).size)
        residuals = self._optional_controller_vector(
            controller,
            ('getResidualNorms', 'getResiduals', 'getEigenpairResiduals'),
            nstates,
        )
        converged = self._optional_controller_vector(
            controller,
            ('getConvergedRoots', 'getRootConvergence'),
            nstates,
        )
        return {
            'solver_residual': residuals,
            'solver_converged': converged,
            's2': np.asarray(self.state_s2, dtype=float),
            'inferred_multiplicity': np.asarray(
                self.state_multiplicities, dtype=int),
            's2_deviation': np.asarray(self.s2_deviation, dtype=float),
            'spin_compatible': np.asarray(spin_compatible, dtype=bool),
        }

    @staticmethod
    def _optional_controller_vector(controller, method_names, nstates):
        for method_name in method_names:
            method = getattr(controller, method_name, None)
            if method is None:
                continue
            try:
                values = np.asarray(method(), dtype=float).reshape(-1)
            except Exception:
                continue
            if values.size == nstates:
                return values
        return None

    def _apply_tracking_failure_policy(self, result):
        status = result.status
        selected_state = result.get('new_state')
        if status is StateTrackingStatus.CONFIDENT:
            return int(selected_state), result

        if status in (StateTrackingStatus.NO_ELIGIBLE_ROOT,
                      StateTrackingStatus.INVALID_RESPONSE):
            raise self._tracking_error(result)

        policy = self.state_tracker.failure_policy
        if (policy == 'strict' or
                (status is StateTrackingStatus.LOW_OVERLAP and
                 getattr(self, '_tracking_recovery_active', False))):
            raise self._tracking_error(result)

        diagnostics = result.to_dict()
        diagnostics.pop('status', None)
        diagnostics['provisional'] = True
        diagnostics['reference_update_recommended'] = False
        diagnostics['policy_action'] = policy

        if policy == 'best_effort':
            if status is StateTrackingStatus.AMBIGUOUS:
                raise self._tracking_error(result)
            diagnostics['warnings'] = list(diagnostics['warnings']) + [
                'BEST_EFFORT: using the unique argmax provisionally; the last '
                'valid reference is preserved'
            ]
            provisional = StateTrackingResult(status, diagnostics)
            self._print_provisional_warning(provisional)
            return int(selected_state), provisional

        # Adiabatic fallback: energy-order the inferred target-spin roots.  No
        # rotated diabatic state is constructed; a real Serenity raw-root
        # gradient is always requested.
        candidates = [
            row for row in diagnostics['candidate_table']
            if row['eligible'] and row['spin_compatible'] is True and
            row['excitation_energy_ev'] is not None
        ]
        if not candidates:
            raise self._tracking_error(result)
        adiabatic = min(
            candidates, key=lambda row: row['excitation_energy_ev'])
        selected_state = int(adiabatic['raw_root'])
        diagnostics.update({
            'new_state': selected_state,
            'selected_by_overlap': int(result['new_state']),
            'reference_update_recommended': True,
        })
        diagnostics['warnings'] = list(diagnostics['warnings']) + [
            'ADIABATIC: following the lowest-energy inferred target-spin '
            f'raw root {selected_state}'
        ]
        provisional = StateTrackingResult(status, diagnostics)
        self._print_provisional_warning(provisional)
        return selected_state, provisional

    def _tracking_error(self, result):
        message = (
            'Serenity state tracking failed with status '
            f'{result.status.value}: reference root {result.get("old_state")}, '
            f'candidate root {result.get("new_state")}, normalized overlap '
            f'{result.get("max_overlap")}, second/best ratio '
            f'{result.get("overlap_ratio")}.')
        self.tracking_info = result.to_dict()
        self._print_tracking_diagnostics(self.tracking_info)
        return StateTrackingError(message, result)

    def _print_provisional_warning(self, result):
        self.ostream.print_info(
            'WARNING: Serenity state assignment is provisional '
            f'({result.status.value}, policy='
            f'{self.state_tracker.failure_policy}).')
        self.ostream.flush()

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
        grad_task.settings.lrscfSettings.restart = True
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
            if hasattr(grad_task.settings.lrscfSettings, 'scfstab'):
                grad_task.settings.lrscfSettings.scfstab = 'spinflip'
            elif hasattr(grad_task.settings, 'scfstab'):
                grad_task.settings.scfstab = 'spinflip'

    def _extract_excitation_energy(self, rsp_results):
        eig = rsp_results['eigenvalues']

        state = int(self.state_deriv_index)
        nst = len(eig)

        errmsg = 'SerenityExcitedStateGradientDriver: requested state index '
        errmsg += f'{state} but only {nst} state(s) are available.'
        assert_msg_critical(state <= nst, errmsg)

        return float(eig[state - 1])
    
    def set_state_tracker(self, tracker):
        assert_msg_critical(
            tracker is None or isinstance(tracker, TransitionDensityTracker),
            "SerenityExcitedStateGradientDriver: invalid transition-density tracker.")
        self.state_tracker = tracker
        self.tracking_info = None
        self.tracking_history = []

    def _print_tracking_diagnostics(self, info):
        """Prints one Serenity transition-density tracking decision."""

        self.ostream.print_header('Serenity Spin-Flip State Tracking')
        self.ostream.print_info(
            f"Status                     : {info.get('status', 'CONFIDENT')}")
        self.ostream.print_info(
            f"Previous raw root          : {info['old_state']}")
        self.ostream.print_info(
            f"Selected candidate root    : {info.get('new_state')}")
        self.ostream.print_info(
            'Hard-eligible roots         : ' +
            ' '.join(str(state)
                     for state in info.get('allowed_states', [])))
        if info.get('spin_diagnostic_roots') is not None:
            self.ostream.print_info(
                'Target-spin diagnostics : ' +
                ' '.join(str(state)
                         for state in info['spin_diagnostic_roots']))
        if info.get('s2_tolerance_roots') is not None:
            # The roots that are actually spin-compatible: target multiplicity
            # and within the <S^2> tolerance.  Roots listed above but missing
            # here are spin-incomplete SF roots.
            self.ostream.print_info(
                'Spin-compatible roots   : ' +
                ' '.join(str(state)
                         for state in info['s2_tolerance_roots']))
        if info.get('initialized', False):
            self.ostream.print_info(
                'Tracking reference initialized at the current geometry.')
        else:
            if info.get('max_overlap') is not None:
                self.ostream.print_info(
                    f"Target overlap             : "
                    f"{info['max_overlap']:.8f}")
            if info.get('second_overlap') is not None:
                self.ostream.print_info(
                    f"Second eligible overlap    : "
                    f"{info['second_overlap']:.8f}")
            if info.get('overlap_ratio') is not None:
                self.ostream.print_info(
                    f"Ambiguity ratio            : "
                    f"{info['overlap_ratio']:.8f}")
            self.ostream.print_info(
                f"Assignment confident       : "
                f"{info.get('assignment_confident', False)}")
            if info.get('global_state') is not None:
                self.ostream.print_info(
                    f"Hungarian target root      : {info['global_state']}")
        self.ostream.print_info(
            f"Raw root changed           : {info.get('swapped', False)}")
        if info.get('gradient_task_roots') is not None:
            self.ostream.print_info(
                'Gradient task root sequence: ' +
                ' '.join(str(state)
                         for state in info['gradient_task_roots']))

        if info.get('candidate_table'):
            self.ostream.print_info(
                'Candidates: raw  energy/eV residual conv       <S^2> mult '
                'raw-overlap norm-overlap eligible exclusion')
            for row in info['candidate_table']:
                def value(key, fmt):
                    item = row.get(key)
                    return 'n/a' if item is None else format(item, fmt)

                self.ostream.print_info(
                    f"  {row['raw_root']:3d} "
                    f"{value('excitation_energy_ev', '.8f'):>12s} "
                    f"{value('solver_residual', '.3e'):>8s} "
                    f"{str(row.get('solver_converged')):>5s} "
                    f"{value('s2', '.6f'):>10s} "
                    f"{str(row.get('inferred_multiplicity')):>4s} "
                    f"{value('raw_overlap', '.8f'):>11s} "
                    f"{value('normalized_overlap', '.8f'):>12s} "
                    f"{str(row['eligible']):>8s} "
                    f"{row.get('exclusion_reason') or '-'}")
        self.ostream.print_blank()
        self.ostream.flush()

    def commit_tracking_reference(self):
        """Commits the staged reference after an accepted nuclear step."""

        if self.state_tracker is None:
            return False
        if self.rank == mpi_master():
            committed = self.state_tracker.commit()
            reference_state = int(self.state_tracker.reference_state)
            if committed and self.tracking_info is not None:
                self.tracking_info['reference_updated'] = True
                self.tracking_info['reference_staged'] = False
                if self.tracking_history:
                    self.tracking_history[-1]['reference_updated'] = True
                    self.tracking_history[-1]['reference_staged'] = False
        else:
            committed = None
            reference_state = None
        committed, reference_state = self.comm.bcast(
            (committed, reference_state), root=mpi_master())
        if committed:
            self.state_deriv_index = int(reference_state)
        return committed

    def rollback_tracking_reference(self):
        """Rolls back a rejected trial and its selected raw-root identity."""

        if self.state_tracker is None:
            return False
        if self.rank == mpi_master():
            rolled_back = self.state_tracker.rollback()
            reference_state = int(self.state_tracker.reference_state)
        else:
            rolled_back = None
            reference_state = None
        rolled_back, reference_state = self.comm.bcast(
            (rolled_back, reference_state), root=mpi_master())
        self.state_deriv_index = int(reference_state)
        self._tracking_applied_in_compute = False
        return rolled_back

    def begin_tracking_retry(self):
        """Restores the committed root and starts midpoint subdivision."""

        if self.state_tracker is None:
            return False
        if self.rank == mpi_master():
            self.state_tracker.begin_retry_chain()
            reference_state = int(self.state_tracker.reference_state)
        else:
            reference_state = None
        reference_state = self.comm.bcast(
            reference_state, root=mpi_master())
        self.state_deriv_index = int(reference_state)
        return True

    def promote_tracking_reference_for_retry(self):
        """Promotes a successful midpoint only within the retry transaction."""

        if self.state_tracker is None:
            return False
        if self.rank == mpi_master():
            promoted = self.state_tracker.promote_pending_for_retry()
            reference_state = int(self.state_tracker.reference_state)
        else:
            promoted = None
            reference_state = None
        promoted, reference_state = self.comm.bcast(
            (promoted, reference_state), root=mpi_master())
        if promoted:
            self.state_deriv_index = int(reference_state)
        return promoted

    def track_state(self, molecule, recompute_on_switch=True,
                update_reference=True):
        if self.state_tracker is None:
            self.tracking_info = None
            return None

        # Spin-flip calculations with multiplicity enforcement are tracked as
        # part of compute(), so energy, gradient, root and spin diagnostics are
        # selected atomically. OpenMM may still call track_state() afterwards;
        # return the already completed decision instead of tracking twice.
        if (self.rsp_driver.spinflip and
                self.enforce_same_multiplicity and
                self._tracking_applied_in_compute):
            if update_reference:
                self.commit_tracking_reference()
            return self.tracking_info

        assert_msg_critical(
            self.last_lrscf_controller is not None,
            "track_state: call compute() before track_state()."
        )
        old_state = int(self.state_deriv_index)
        result = self.state_tracker.track(
            self._last_tracking_system,
            self.last_lrscf_controller,
            self._last_tracking_mode,
            active_reference_state=old_state,
        )
        info = result.to_dict()
        new_state = int(info["new_state"])
        info["gradient_recomputed"] = False
        info["reference_updated"] = False
        if new_state != old_state:
            self.set_state_deriv_index(new_state)
            if recompute_on_switch:
                self.gradient = self._compute_analytical_master(molecule)
                info["gradient_recomputed"] = True
        if update_reference:
            self.state_tracker.accept_reference(
                self._last_tracking_system,
                self.last_lrscf_controller,
                self._last_tracking_mode,
                int(self.state_deriv_index),
            )
            info["reference_updated"] = True
        self.tracking_info = info

        return self.tracking_info
