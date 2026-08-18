#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2026 VeloxChem developers
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

"""Serenity transition-density state tracking.

Serenity returns a normalized transition-density overlap matrix with current
roots in rows and reference roots in columns.  This module owns assignment,
confidence classification, diagnostics, and transactional reference lifetime;
the electronic-structure backend remains responsible for constructing the
cross-geometry AO metric.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum

import numpy as np

from .errorhandler import assert_msg_critical
from .serenitylrrspeigensolver import SerenityLinearResponseSolver
from .surfacehoppingtracker import solve_assignment
from .veloxchemlib import hartree_in_ev

try:
    from qcserenity import serenipy as spy
except ImportError:
    pass


class StateTrackingStatus(str, Enum):
    """Mutually exclusive outcome of one state-assignment attempt."""

    CONFIDENT = 'CONFIDENT'
    LOW_OVERLAP = 'LOW_OVERLAP'
    AMBIGUOUS = 'AMBIGUOUS'
    SPIN_CONFLICT = 'SPIN_CONFLICT'
    GROUND_STATE_COLLISION = 'GROUND_STATE_COLLISION'
    NO_ELIGIBLE_ROOT = 'NO_ELIGIBLE_ROOT'
    INVALID_RESPONSE = 'INVALID_RESPONSE'


@dataclass(frozen=True)
class StateTrackingResult(Mapping):
    """Structured tracking result with mapping compatibility for older code."""

    status: StateTrackingStatus
    diagnostics: dict

    def to_dict(self):
        result = dict(self.diagnostics)
        result['status'] = self.status.value
        return result

    def __getitem__(self, key):
        if key == 'status':
            return self.status.value
        return self.diagnostics[key]

    def __iter__(self):
        yield 'status'
        yield from self.diagnostics

    def __len__(self):
        return len(self.diagnostics) + 1

    @property
    def is_confident(self):
        return self.status is StateTrackingStatus.CONFIDENT


class StateTrackingError(RuntimeError):
    """Typed failure carrying the complete state-tracking diagnostics."""

    def __init__(self, message, result):
        super().__init__(message)
        self.result = result
        self.status = result.status
        self.diagnostics = result.to_dict()
        self.recoverable = result.status is StateTrackingStatus.LOW_OVERLAP


class TransitionDensityTracker:
    """Tracks a Serenity LR root by normalized transition-density overlap.

    ``failure_policy`` controls what the gradient driver does after a
    non-confident assignment.  ``strict`` is the default, ``best_effort`` uses
    a unique argmax provisionally, and ``adiabatic`` selects the lowest-energy
    target-spin root.  Low-overlap subdivision is bounded by
    ``max_subdivisions`` and is performed by :class:`OptimizationEngine`.
    """

    def __init__(self,
                 lr_resp_driver,
                 target_state=1,
                 min_overlap=0.5,
                 overlap_ratio_bounds=(0.3, 0.6),
                 degeneracy_threshold=None,
                 lrscf_type='isolated',
                 failure_policy='strict',
                 max_subdivisions=2,
                 max_root_expansions=0,
                 root_window_increment=5,
                 guard_ground_state=True,
                 ground_state_root=None):

        assert_msg_critical(
            isinstance(lr_resp_driver, SerenityLinearResponseSolver),
            'TransitionDensityTracker: invalid LR response driver.')

        self.lrresp_driver = lr_resp_driver
        self.target_state = int(target_state)
        self.min_overlap = float(min_overlap)
        self.overlap_ratio_bounds = tuple(overlap_ratio_bounds)
        self.degeneracy_threshold = degeneracy_threshold
        self.lrscf_type = lrscf_type
        self.failure_policy = str(failure_policy).strip().lower()
        self.max_subdivisions = int(max_subdivisions)
        self.max_root_expansions = int(max_root_expansions)
        self.root_window_increment = int(root_window_increment)
        self.guard_ground_state = bool(guard_ground_state)
        self.ground_state_root = (None if ground_state_root is None else
                                  int(ground_state_root))

        assert_msg_critical(
            self.target_state > 0,
            'TransitionDensityTracker: target_state must be positive.')
        assert_msg_critical(
            0.0 <= self.min_overlap <= 1.0,
            'TransitionDensityTracker: min_overlap must be in [0, 1].')
        assert_msg_critical(
            len(self.overlap_ratio_bounds) == 2 and
            0.0 <= float(self.overlap_ratio_bounds[0]) <=
            float(self.overlap_ratio_bounds[1]) <= 1.0,
            'TransitionDensityTracker: invalid overlap_ratio_bounds.')
        assert_msg_critical(
            self.failure_policy in ('strict', 'best_effort', 'adiabatic'),
            'TransitionDensityTracker: failure_policy must be strict, '
            'best_effort, or adiabatic.')
        assert_msg_critical(
            self.max_subdivisions >= 0 and self.max_root_expansions >= 0 and
            self.root_window_increment > 0,
            'TransitionDensityTracker: invalid retry settings.')

        self._ref_exc = None
        self._ref_geom = None
        self._ref_coeff = None
        self._ref_nstates = None
        self._ref_energies_ev = None
        self.reference_state = self.target_state
        self.last_info = None

        self._pending = None
        self._retry_base = None
        self._reference_revision = 0
        self.rejected_commits = 0

    def has_reference(self):
        return self._ref_exc is not None

    @property
    def has_pending(self):
        return self._pending is not None

    @property
    def reference_revision(self):
        return self._reference_revision

    def set_target_state(self, state):
        self.target_state = int(state)

    def _capture_reference(self, system, lrscf_controller, mode,
                           reference_state):
        del mode  # The reference objects themselves encode restricted/unrestricted.
        excitation_vectors = lrscf_controller.getExcitationVectors(
            self.lrscf_type)
        return {
            'exc': excitation_vectors,
            'nstates': self._number_of_states(excitation_vectors),
            'geom': system.getGeometry().copy(),
            'coeff': lrscf_controller.getCoefficients(),
            'energies_ev': np.asarray(
                lrscf_controller.getExcitationEnergies(self.lrscf_type),
                dtype=float).reshape(-1).copy(),
            'reference_state': int(reference_state),
        }

    def _install_reference(self, snapshot):
        self._ref_exc = snapshot['exc']
        self._ref_nstates = int(snapshot['nstates'])
        self._ref_geom = snapshot['geom']
        self._ref_coeff = snapshot['coeff']
        self._ref_energies_ev = snapshot['energies_ev']
        self.reference_state = int(snapshot['reference_state'])

    def _current_reference_snapshot(self):
        if not self.has_reference():
            return None
        return {
            'exc': self._ref_exc,
            'nstates': self._ref_nstates,
            'geom': self._ref_geom,
            'coeff': self._ref_coeff,
            'energies_ev': self._ref_energies_ev,
            'reference_state': self.reference_state,
        }

    def initialize_reference(self, system, lrscf_controller, mode,
                             reference_state=None):
        state = (self.target_state if reference_state is None else
                 int(reference_state))
        self._install_reference(
            self._capture_reference(system, lrscf_controller, mode, state))
        self._reference_revision += 1

    def propose_reference(self, system, lrscf_controller, mode,
                          reference_state=None):
        state = (self.reference_state if reference_state is None else
                 int(reference_state))
        self._pending = self._capture_reference(
            system, lrscf_controller, mode, state)
        return True

    def commit(self):
        if self._pending is None:
            return False
        self._install_reference(self._pending)
        self._pending = None
        self._retry_base = None
        self._reference_revision += 1
        return True

    def rollback(self):
        changed = self._pending is not None or self._retry_base is not None
        self._pending = None
        if self._retry_base is not None:
            self._install_reference(self._retry_base)
            self._retry_base = None
        return changed

    def begin_retry_chain(self):
        """Starts an ephemeral midpoint-reference chain."""

        self.rollback()
        self._retry_base = self._current_reference_snapshot()

    def promote_pending_for_retry(self):
        """Uses a successful midpoint as a working, not committed, reference."""

        if self._pending is None:
            return False
        if self._retry_base is None:
            self._retry_base = self._current_reference_snapshot()
        self._install_reference(self._pending)
        self._pending = None
        return True

    def accept_reference(self,
                         system,
                         lrscf_controller,
                         mode,
                         reference_state=None):
        """Immediately accepts a reference for non-transactional callers."""

        self.propose_reference(system, lrscf_controller, mode, reference_state)
        return self.commit()

    def track(self,
              system,
              lrscf_controller,
              mode,
              active_reference_state=None,
              allowed_states=None,
              candidate_metadata=None):
        """Assigns the committed reference character to one current root."""

        current_exc = lrscf_controller.getExcitationVectors(self.lrscf_type)
        current_nstates = self._number_of_states(current_exc)
        energies_ev = np.asarray(
            lrscf_controller.getExcitationEnergies(self.lrscf_type),
            dtype=float).reshape(-1)
        response_valid, vector_norms, response_reasons = \
            self._response_validity(current_exc, energies_ev)
        response_valid, response_reasons = self._apply_solver_validity(
            response_valid, response_reasons, candidate_metadata)

        if not self.has_reference():
            ref_state = (self.target_state if active_reference_state is None
                         else int(active_reference_state))
            assert_msg_critical(
                1 <= ref_state <= current_nstates,
                'TransitionDensityTracker: initial reference state is outside '
                'the current excitation-vector range.')
            allowed = self._allowed_mask(allowed_states, current_nstates)
            eligible = allowed & response_valid
            status = (StateTrackingStatus.CONFIDENT
                      if eligible[ref_state - 1] else
                      StateTrackingStatus.INVALID_RESPONSE)
            candidates = self._candidate_table(
                energies_ev, np.full(current_nstates, np.nan),
                np.full(current_nstates, np.nan), eligible, response_valid,
                response_reasons, vector_norms, candidate_metadata)
            diagnostics = {
                'initialized': True,
                'swapped': False,
                'old_state': ref_state,
                'new_state': ref_state,
                'reference_state': ref_state,
                'max_overlap': 1.0,
                'second_overlap': 0.0,
                'overlap_ratio': 0.0,
                'overlap_matrix': None,
                'assignment_confident': status is StateTrackingStatus.CONFIDENT,
                'reference_update_recommended': False,
                'degenerate': False,
                'allowed_states': [
                    int(i) for i in np.flatnonzero(eligible) + 1
                ],
                'candidate_table': candidates,
                'greedy_state': ref_state,
                'global_state': ref_state,
                'global_assignment': list(range(1, current_nstates + 1)),
                'global_assignment_disagrees': False,
                'subspace': None,
                'warnings': [],
                'provisional': False,
            }
            result = StateTrackingResult(status, diagnostics)
            self.last_info = result
            if result.is_confident:
                # The optimizer's initial point is the accepted baseline; no
                # accepted-step callback precedes its first trial geometry.
                self.initialize_reference(
                    system, lrscf_controller, mode, ref_state)
            return result

        # Once initialized, the committed reference owns the column identity.
        # A caller's most recently evaluated raw root may belong to a rejected
        # trial and must never redirect the reference column.
        ref_state = int(self.reference_state)
        assert_msg_critical(
            1 <= ref_state <= int(self._ref_nstates),
            'TransitionDensityTracker: committed reference state is outside '
            'the stored excitation-vector range.')

        tracker_cls = (spy.TransitionDensityTracking_R
                       if mode == 'restricted' else
                       spy.TransitionDensityTracking_U)
        degeneracy_threshold = self._degeneracy_threshold()
        try:
            tden = tracker_cls(
                system,
                lrscf_controller,
                self._ref_exc,
                self._ref_geom,
                self._ref_coeff,
                spy.LRSCF_TYPE.ISOLATED,
                ref_state,
                float(degeneracy_threshold),
            )
            tden.track()
            overlap = np.asarray(tden.getOvlpMatrix())
        except Exception as exc:
            result = self._invalid_response_result(
                current_nstates, energies_ev, response_valid,
                response_reasons, vector_norms, candidate_metadata,
                f'Serenity transition-density tracking failed: {exc}')
            self.last_info = result
            return result

        result = self.evaluate_overlap_matrix(
            overlap,
            reference_state=ref_state,
            current_energies_ev=energies_ev,
            allowed_states=allowed_states,
            response_valid_mask=response_valid,
            response_reasons=response_reasons,
            vector_norms=vector_norms,
            candidate_metadata=candidate_metadata,
            degeneracy_threshold=degeneracy_threshold,
        )
        self.last_info = result
        return result

    def evaluate_overlap_matrix(self,
                                overlap_matrix,
                                reference_state=None,
                                current_energies_ev=None,
                                allowed_states=None,
                                response_valid_mask=None,
                                response_reasons=None,
                                vector_norms=None,
                                candidate_metadata=None,
                                degeneracy_threshold=None):
        """Classifies a supplied matrix; useful for deterministic replays."""

        overlap = np.asarray(overlap_matrix)
        if overlap.ndim != 2:
            return self._minimal_invalid_result(
                'Serenity returned a non-matrix transition-density overlap.')

        n_current, n_reference = overlap.shape
        ref_state = (self.reference_state if reference_state is None else
                     int(reference_state))
        if ref_state < 1 or ref_state > n_reference or n_current == 0:
            return self._minimal_invalid_result(
                'The reference state is outside the overlap matrix.')

        if np.iscomplexobj(overlap) and np.any(
                np.abs(np.imag(overlap)) > 1.0e-12):
            return self._minimal_invalid_result(
                'Serenity returned complex transition-density overlaps.')
        overlap = np.asarray(np.real(overlap), dtype=float)

        if current_energies_ev is None:
            energies_ev = np.full(n_current, np.nan)
        else:
            energies_ev = np.asarray(
                current_energies_ev, dtype=float).reshape(-1)
        if energies_ev.size != n_current:
            return self._minimal_invalid_result(
                'The energy and overlap root windows have different sizes.')

        allowed = self._allowed_mask(allowed_states, n_current)
        response_valid = (np.ones(n_current, dtype=bool)
                          if response_valid_mask is None else
                          np.asarray(response_valid_mask,
                                     dtype=bool).reshape(-1))
        if response_valid.size != n_current:
            return self._minimal_invalid_result(
                'The response-validity mask has the wrong size.')
        reasons = ([''] * n_current if response_reasons is None else
                   list(response_reasons))
        norms = (np.full(n_current, np.nan) if vector_norms is None else
                 np.asarray(vector_norms, dtype=float).reshape(-1))
        response_valid, reasons = self._apply_solver_validity(
            response_valid, reasons, candidate_metadata)

        raw_column = overlap[:, ref_state - 1]
        finite_overlap = np.isfinite(raw_column)
        eligible = allowed & response_valid & finite_overlap
        normalized_column = np.abs(raw_column)

        # Keep an excited-state run off S0.  The guard is skipped when the
        # ground state is what the tracker was asked to follow.
        ground_index = self._ground_state_index(energies_ev, n_current)
        guarding_ground_state = (
            self.guard_ground_state and
            ground_index is not None and
            int(self.target_state) != ground_index + 1 and
            bool(eligible[ground_index]))
        ground_state_collision = False
        if guarding_ground_state:
            unguarded_scores = np.where(eligible, normalized_column, -np.inf)
            ground_state_collision = (
                int(np.argmax(unguarded_scores)) == ground_index)
            eligible = eligible.copy()
            eligible[ground_index] = False
        candidates = self._candidate_table(
            energies_ev, raw_column, normalized_column, eligible,
            response_valid & finite_overlap, reasons, norms,
            candidate_metadata)

        if not np.any(eligible):
            status = (StateTrackingStatus.INVALID_RESPONSE
                      if np.any(allowed) else
                      StateTrackingStatus.NO_ELIGIBLE_ROOT)
            diagnostics = self._base_diagnostics(
                ref_state, None, overlap, candidates, eligible)
            diagnostics['warnings'] = [
                'No valid current response root is available for assignment.'
            ]
            return StateTrackingResult(status, diagnostics)

        scores = np.where(eligible, normalized_column, -np.inf)
        new_index = int(np.argmax(scores))
        new_state = new_index + 1
        best = float(normalized_column[new_index])
        eligible_indices = np.flatnonzero(eligible)
        alternatives = eligible_indices[eligible_indices != new_index]
        if alternatives.size:
            second_index = int(
                alternatives[np.argmax(normalized_column[alternatives])])
            second = float(normalized_column[second_index])
            second_state = second_index + 1
        else:
            second = 0.0
            second_state = None
        ratio = float(second / best) if best > 0.0 else np.inf

        threshold = (self._degeneracy_threshold()
                     if degeneracy_threshold is None else
                     float(degeneracy_threshold))
        degenerate = self._energies_are_degenerate(
            energies_ev, new_state, second_state, threshold)

        _, maximum_ratio = self.overlap_ratio_bounds
        if ground_state_collision:
            # Reported ahead of the other diagnostics: once the tracked
            # character sits on S0 the remaining overlap numbers describe a
            # fallback assignment, not the state that was asked for.
            status = StateTrackingStatus.GROUND_STATE_COLLISION
        elif ratio > float(maximum_ratio) or degenerate:
            status = StateTrackingStatus.AMBIGUOUS
        elif not np.isfinite(best) or best < self.min_overlap:
            status = StateTrackingStatus.LOW_OVERLAP
        elif not self._candidate_spin_compatible(
                candidate_metadata, new_index):
            status = StateTrackingStatus.SPIN_CONFLICT
        else:
            status = StateTrackingStatus.CONFIDENT

        global_assignment = None
        global_state = None
        if n_current == n_reference and np.all(np.isfinite(overlap)):
            global_assignment = solve_assignment(np.abs(overlap).T) + 1
            global_state = int(global_assignment[ref_state - 1])

        subspace = self._subspace_diagnostics(
            overlap, energies_ev, ref_state, new_state, threshold)
        warnings = self._status_warnings(
            status, best, ratio, new_state,
            None if ground_index is None else ground_index + 1)
        diagnostics = {
            'initialized': False,
            'ground_state_root': (None if ground_index is None else
                                  ground_index + 1),
            'ground_state_guarded': bool(guarding_ground_state),
            'ground_state_collision': bool(ground_state_collision),
            'swapped': new_state != ref_state,
            'old_state': ref_state,
            'new_state': new_state,
            'reference_state': ref_state,
            'max_overlap': best,
            'second_overlap': second,
            'second_state': second_state,
            'overlap_ratio': ratio,
            'overlap_matrix': overlap,
            'assignment_confident': status is StateTrackingStatus.CONFIDENT,
            'reference_update_recommended':
                status is StateTrackingStatus.CONFIDENT,
            'degenerate': degenerate,
            'allowed_states': [
                int(i) for i in np.flatnonzero(eligible) + 1
            ],
            'candidate_table': candidates,
            'greedy_state': new_state,
            'global_state': global_state,
            'global_assignment': (None if global_assignment is None else
                                  global_assignment.astype(int).tolist()),
            'global_assignment_disagrees':
                global_state is not None and global_state != new_state,
            'subspace': subspace,
            'warnings': warnings,
            'provisional': False,
        }
        return StateTrackingResult(status, diagnostics)

    def _invalid_response_result(self, nstates, energies_ev, response_valid,
                                 response_reasons, vector_norms,
                                 candidate_metadata, message):
        candidates = self._candidate_table(
            energies_ev, np.full(nstates, np.nan),
            np.full(nstates, np.nan), np.zeros(nstates, dtype=bool),
            response_valid, response_reasons, vector_norms,
            candidate_metadata)
        diagnostics = self._base_diagnostics(
            self.reference_state, None, None, candidates,
            np.zeros(nstates, dtype=bool))
        diagnostics['warnings'] = [message]
        return StateTrackingResult(
            StateTrackingStatus.INVALID_RESPONSE, diagnostics)

    def _minimal_invalid_result(self, message):
        diagnostics = self._base_diagnostics(
            self.reference_state, None, None, [], np.zeros(0, dtype=bool))
        diagnostics['warnings'] = [message]
        return StateTrackingResult(
            StateTrackingStatus.INVALID_RESPONSE, diagnostics)

    @staticmethod
    def _base_diagnostics(ref_state, new_state, overlap, candidates,
                          eligible):
        return {
            'initialized': False,
            'swapped': new_state is not None and new_state != ref_state,
            'old_state': int(ref_state),
            'new_state': new_state,
            'reference_state': int(ref_state),
            'max_overlap': None,
            'second_overlap': None,
            'second_state': None,
            'overlap_ratio': None,
            'overlap_matrix': overlap,
            'assignment_confident': False,
            'reference_update_recommended': False,
            'degenerate': False,
            'allowed_states': [
                int(i) for i in np.flatnonzero(eligible) + 1
            ],
            'candidate_table': candidates,
            'greedy_state': new_state,
            'global_state': None,
            'global_assignment': None,
            'global_assignment_disagrees': False,
            'subspace': None,
            'warnings': [],
            'provisional': False,
        }

    def _candidate_table(self, energies_ev, raw_overlap, normalized_overlap,
                         eligible, response_valid, response_reasons,
                         vector_norms, metadata):
        nstates = len(energies_ev)
        table = []
        for index in range(nstates):
            response_ok = bool(response_valid[index])
            hard_allowed = bool(eligible[index])
            if not response_ok:
                reason = (response_reasons[index]
                          if index < len(response_reasons) else
                          'invalid response')
            elif not hard_allowed:
                reason = 'excluded by exact eligibility mask'
            else:
                reason = None
            row = {
                'raw_root': index + 1,
                'excitation_energy_ev': self._finite_or_none(
                    energies_ev[index]),
                'excitation_energy_hartree': self._finite_or_none(
                    energies_ev[index] / hartree_in_ev()),
                'solver_residual': self._metadata_value(
                    metadata, 'solver_residual', index),
                'solver_converged': self._metadata_value(
                    metadata, 'solver_converged', index),
                'response_vector_norm': self._finite_or_none(
                    vector_norms[index] if index < len(vector_norms) else
                    np.nan),
                's2': self._metadata_value(metadata, 's2', index),
                'inferred_multiplicity': self._metadata_value(
                    metadata, 'inferred_multiplicity', index),
                's2_deviation': self._metadata_value(
                    metadata, 's2_deviation', index),
                'spin_compatible': self._metadata_value(
                    metadata, 'spin_compatible', index),
                # The native API exposes one normalized matrix only.  Preserve
                # its signed/raw entry and the phase-invariant magnitude.
                'raw_overlap': self._finite_or_none(raw_overlap[index]),
                'normalized_overlap': self._finite_or_none(
                    normalized_overlap[index]),
                'eligible': hard_allowed,
                'exclusion_reason': reason,
            }
            table.append(row)
        return table

    @staticmethod
    def _metadata_value(metadata, key, index):
        if metadata is None or key not in metadata:
            return None
        values = metadata[key]
        if values is None:
            return None
        flattened = np.asarray(values).reshape(-1)
        if index >= flattened.size:
            return None
        value = flattened[index]
        if isinstance(value, (np.bool_, bool)):
            return bool(value)
        if isinstance(value, (np.integer, int)):
            return int(value)
        if isinstance(value, (np.floating, float)):
            return TransitionDensityTracker._finite_or_none(value)
        return value

    @staticmethod
    def _finite_or_none(value):
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        return number if np.isfinite(number) else None

    @staticmethod
    def _candidate_spin_compatible(metadata, index):
        if metadata is None or metadata.get('spin_compatible') is None:
            return True
        return bool(np.asarray(metadata['spin_compatible']).reshape(-1)[index])

    @staticmethod
    def _allowed_mask(allowed_states, nstates):
        if allowed_states is None:
            return np.ones(nstates, dtype=bool)
        mask = np.zeros(nstates, dtype=bool)
        try:
            indices = [int(state) for state in allowed_states]
        except (TypeError, ValueError):
            indices = []
        assert_msg_critical(
            all(1 <= state <= nstates for state in indices),
            'TransitionDensityTracker: allowed state index is outside the '
            'current root range.')
        if indices:
            mask[np.asarray(indices, dtype=int) - 1] = True
        return mask

    @staticmethod
    def _response_validity(excitation_vectors, energies_ev):
        nstates = energies_ev.size
        valid = np.isfinite(energies_ev)
        squared_norms = np.zeros(nstates, dtype=float)
        reasons = ['' for _ in range(nstates)]
        for component in excitation_vectors:
            array = np.asarray(component)
            if array.ndim != 2 or array.shape[1] != nstates:
                return (np.zeros(nstates, dtype=bool),
                        np.full(nstates, np.nan),
                        ['inconsistent response-vector shape'] * nstates)
            complex_columns = (np.any(np.abs(np.imag(array)) > 1.0e-12,
                                      axis=0)
                               if np.iscomplexobj(array) else
                               np.zeros(nstates, dtype=bool))
            finite_columns = np.all(np.isfinite(array), axis=0)
            valid &= finite_columns & ~complex_columns
            squared_norms += np.sum(np.abs(array)**2, axis=0)
            for index in np.flatnonzero(complex_columns):
                reasons[index] = 'complex response vector'
            for index in np.flatnonzero(~finite_columns):
                reasons[index] = 'nonfinite response vector'
        norms = np.sqrt(squared_norms)
        zero_norm = ~np.isfinite(norms) | (norms <= 1.0e-14)
        valid &= ~zero_norm
        for index in np.flatnonzero(zero_norm):
            reasons[index] = 'zero or invalid response-vector norm'
        for index in np.flatnonzero(~np.isfinite(energies_ev)):
            reasons[index] = 'nonfinite excitation energy'
        return valid, norms, reasons

    def _apply_solver_validity(self, valid, reasons, metadata):
        """Hard-excludes explicitly non-converged or invalid eigenpairs."""

        valid = np.asarray(valid, dtype=bool).copy()
        reasons = list(reasons)
        nstates = valid.size
        if metadata is None:
            return valid, reasons

        converged = metadata.get('solver_converged')
        if converged is not None:
            values = np.asarray(converged).reshape(-1)
            if values.size == nstates:
                failed = ~values.astype(bool)
                valid &= ~failed
                for index in np.flatnonzero(failed):
                    reasons[index] = 'response root did not converge'

        residuals = metadata.get('solver_residual')
        if residuals is not None:
            values = np.asarray(residuals, dtype=float).reshape(-1)
            if values.size == nstates:
                invalid = ~np.isfinite(values)
                threshold = float(self.lrresp_driver.conv_thresh or 1.0e-6)
                unconverged = values > threshold
                valid &= ~(invalid | unconverged)
                for index in np.flatnonzero(invalid):
                    reasons[index] = 'nonfinite response residual'
                for index in np.flatnonzero(unconverged & ~invalid):
                    reasons[index] = (
                        f'response residual exceeds {threshold:.3e}')
        return valid, reasons

    def _degeneracy_threshold(self):
        if self.degeneracy_threshold is not None:
            return float(self.degeneracy_threshold)
        return 10.0 * float(self.lrresp_driver.conv_thresh or 1.0e-6)

    @staticmethod
    def _energies_are_degenerate(energies_ev, state_a, state_b, threshold):
        if state_b is None or state_a == state_b:
            return False
        energy_a = energies_ev[state_a - 1] / hartree_in_ev()
        energy_b = energies_ev[state_b - 1] / hartree_in_ev()
        return (np.isfinite(energy_a) and np.isfinite(energy_b) and
                abs(energy_a - energy_b) < float(threshold))

    def _subspace_diagnostics(self, overlap, current_energies_ev, ref_state,
                              new_state, threshold):
        if self._ref_energies_ev is None or threshold <= 0.0:
            return None
        ref_ha = np.asarray(self._ref_energies_ev) / hartree_in_ev()
        current_ha = np.asarray(current_energies_ev) / hartree_in_ev()
        if ref_state > ref_ha.size or new_state > current_ha.size:
            return None
        ref_indices = np.flatnonzero(
            np.abs(ref_ha - ref_ha[ref_state - 1]) <= threshold)
        current_indices = np.flatnonzero(
            np.abs(current_ha - current_ha[new_state - 1]) <= threshold)
        if ref_indices.size < 2 and current_indices.size < 2:
            return None
        # Row/column phase changes are unitary transformations and therefore
        # leave these singular values invariant. Elementwise magnitudes would
        # destroy that property for a genuinely rotated subspace.
        block = overlap[np.ix_(current_indices, ref_indices)]
        singular_values = np.linalg.svd(block, compute_uv=False)
        clipped = np.clip(singular_values, 0.0, 1.0)
        return {
            'reference_roots': (ref_indices + 1).astype(int).tolist(),
            'current_roots': (current_indices + 1).astype(int).tolist(),
            'singular_values': singular_values.astype(float).tolist(),
            'principal_angles_radians':
                np.arccos(clipped).astype(float).tolist(),
            'gradient_is_adiabatic_root_only': True,
        }

    def _ground_state_index(self, energies_ev, n_current):
        """Zero-based index of the root that is the electronic ground state.

        A spin-flip response spectrum contains the closed-shell ground state
        as an ordinary root -- normally the lowest one -- so an excited-state
        run can be assigned onto S0 and then simply minimize the ground state.
        The ground state is identified by energy rather than by a fixed number
        because the raw root order changes along a path; ``ground_state_root``
        pins it explicitly when a spectrum needs that. Ordinary TDA/TDDFT
        spectra contain excited roots only, so their lowest excitation must not
        be mistaken for S0.

        :param energies_ev:
            Excitation energies of the current roots.
        :param n_current:
            Number of current roots.

        :return:
            The zero-based ground-state index, or ``None`` when it cannot be
            identified.
        """

        if self.ground_state_root is not None:
            index = self.ground_state_root - 1
            return index if 0 <= index < n_current else None

        if not bool(getattr(self.lrresp_driver, 'spinflip', False)):
            return None

        energies = np.asarray(energies_ev, dtype=float).reshape(-1)
        if energies.size != n_current or not np.any(np.isfinite(energies)):
            # Without energies the SF convention (lowest root is S0) is the
            # only information available.
            return 0 if n_current else None

        return int(np.nanargmin(energies))

    def _status_warnings(self, status, best, ratio, new_state,
                         ground_state_root=None):
        warnings = []
        if status is StateTrackingStatus.GROUND_STATE_COLLISION:
            warnings.append(
                f'the tracked excited state has collapsed onto ground-state '
                f'root {ground_state_root}; S0 was excluded from the '
                f'assignment and root {new_state} was used instead')
        elif status is StateTrackingStatus.LOW_OVERLAP:
            warnings.append(
                f'unique root {new_state} has normalized overlap {best:.6f}, '
                f'below minimum {self.min_overlap:.6f}')
        elif status is StateTrackingStatus.AMBIGUOUS:
            warnings.append(
                f'assignment is ambiguous (second/best={ratio:.6f})')
        elif status is StateTrackingStatus.SPIN_CONFLICT:
            warnings.append(
                f'best physical continuation root {new_state} conflicts with '
                'the inferred target multiplicity')
        return warnings

    @staticmethod
    def _number_of_states(excitation_vectors):
        assert_msg_critical(
            excitation_vectors is not None and len(excitation_vectors) > 0,
            'TransitionDensityTracker: empty excitation-vector container.')
        first_component = np.asarray(excitation_vectors[0])
        assert_msg_critical(
            first_component.ndim == 2,
            'TransitionDensityTracker: expected a two-dimensional '
            'excitation-vector component.')
        return int(first_component.shape[1])
