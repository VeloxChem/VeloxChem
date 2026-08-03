#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2026 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.

"""Persistent-state assignment for OpenQP MRSF calculations.

The assignment consumes a state-overlap matrix with previous *raw* roots in
its rows and current raw roots in its columns.  After the first root exchange,
however, a persistent state identity is no longer the same thing as a previous
raw row.  This module owns that bookkeeping and deliberately contains no
OpenQP imports, which keeps the assignment logic independently testable.

Two conventions matter and are both implemented here rather than at the call
site, so that they can be tested without a working OpenQP installation:

* :func:`native_overlap_to_previous_current` converts OpenQP's tagarray view
  of ``OQP::td_states_overlap`` into the row/column order documented above.
* :class:`OpenQPMRSFStateTracker` is transactional.  A finished electronic
  evaluation only ever *proposes* a new reference; the reference advances when
  the consumer (for a relaxed optimization, geomeTRIC's accepted-step hook)
  commits it.  A rejected nuclear trial rolls the proposal back, so the next
  geometry is still compared against the last accepted one.
"""

from copy import deepcopy
from dataclasses import dataclass, field
import json

import numpy as np

from .surfacehoppingtracker import solve_assignment


_REQUIRED_SNAPSHOT_KEYS = (
    'OQP::VEC_MO_A',
    'OQP::VEC_MO_B',
    'OQP::E_MO_A',
    'OQP::E_MO_B',
    'OQP::td_bvec_mo',
    'OQP::td_energies',
)

# An overlap of a normalized state with another normalized state cannot exceed
# one.  A small allowance absorbs the TLF approximation and finite Davidson
# convergence; anything beyond it means the native quantity is not what the
# assignment assumes and must surface rather than be silently clipped.
_OVERLAP_MAGNITUDE_TOLERANCE = 1.0e-3


def native_overlap_to_previous_current(state_overlap):
    """Reinterprets an OpenQP tagarray overlap as (previous, current).

    OpenQP builds ``OQP::td_states_overlap`` in ``get_states_overlap.F90`` as
    ``s_st(oi, ni)`` with ``oi`` indexing the *previous* geometry's states and
    ``ni`` the *current* ones.  Fortran stores that column-major, but the
    tagarray hands Python a NumPy array that reports the Fortran shape while
    being read in C order.  The delivered array is therefore the index
    transpose of the Fortran matrix.

    The same wrapper behaviour is directly observable on the rectangular
    ``OQP::td_bvec_mo`` (declared ``bvec_mo(xvec_dim, nstate)``): its
    per-state amplitude vectors are unit-norm only after a C-order
    ``reshape(nstate, xvec_dim)``, never as ``raw[:, state]``.  OpenQP's own
    ``NACME.align_x`` relies on exactly that reshape.

    For a square ``(nstates, nstates)`` matrix the correction reduces to a
    transpose.

    :param state_overlap:
        The array as delivered by ``mol.data['OQP::td_states_overlap']``.

    :return:
        A copy with previous raw roots in rows and current raw roots in
        columns.
    """

    overlap = np.asarray(state_overlap, dtype=float)

    if overlap.ndim != 2:
        raise ValueError(
            'native_overlap_to_previous_current: expected a two-dimensional '
            f'state overlap, got {overlap.ndim} dimension(s).')
    if overlap.shape[0] != overlap.shape[1]:
        raise ValueError(
            'native_overlap_to_previous_current: expected a square state '
            f'overlap, got shape {overlap.shape}.')

    # reshape(d2, d1).T is the general Fortran-order recovery; for the square
    # state-overlap matrix it is the plain transpose, written explicitly here
    # so the intent survives future edits.
    return np.ascontiguousarray(overlap.T)


@dataclass
class OpenQPStateTrackingResult:
    """Result of one MRSF persistent-state assignment.

    Root numbers and the ``tracked_to_raw`` permutation are one-based, matching
    OpenQP's ``properties.grad`` convention. ``is_reliable`` applies to the
    configured target identity; ``global_confidence`` separately exposes the
    weakest assignment among all calculated roots.
    """

    step: int
    target_tracked_root: int
    previous_raw_root: int
    selected_raw_root: int
    tracked_to_raw: np.ndarray
    assignment_scores: np.ndarray
    target_score: float
    second_score: float
    ambiguity_ratio: float
    global_confidence: float
    is_reliable: bool
    reordered: bool
    root_changed: bool
    state_energies: np.ndarray
    similarity_matrix: np.ndarray
    signed_overlap_matrix: np.ndarray = None
    warnings: list = field(default_factory=list)
    reference_initialized: bool = False
    seeded_from_reference: bool = False
    ground_state_raw_root: int = None
    ground_state_collision: bool = False

    def to_dict(self):
        """Returns a JSON-serializable diagnostics dictionary."""

        return {
            'step': int(self.step),
            'target_tracked_root': int(self.target_tracked_root),
            'previous_raw_root': int(self.previous_raw_root),
            'selected_raw_root': int(self.selected_raw_root),
            'tracked_to_raw': self.tracked_to_raw.astype(int).tolist(),
            'assignment_scores':
                self.assignment_scores.astype(float).tolist(),
            'target_score': float(self.target_score),
            'second_score': float(self.second_score),
            'ambiguity_ratio': float(self.ambiguity_ratio),
            'global_confidence': float(self.global_confidence),
            'is_reliable': bool(self.is_reliable),
            'reordered': bool(self.reordered),
            'root_changed': bool(self.root_changed),
            'state_energies': self.state_energies.astype(float).tolist(),
            'similarity_matrix':
                self.similarity_matrix.astype(float).tolist(),
            'signed_overlap_matrix':
                (None if self.signed_overlap_matrix is None else
                 self.signed_overlap_matrix.astype(float).tolist()),
            'warnings': list(self.warnings),
            'reference_initialized': bool(self.reference_initialized),
            'seeded_from_reference': bool(self.seeded_from_reference),
            'ground_state_raw_root':
                (None if self.ground_state_raw_root is None else
                 int(self.ground_state_raw_root)),
            'ground_state_collision': bool(self.ground_state_collision),
        }


class OpenQPMRSFStateTracker:
    """Tracks persistent MRSF states using OpenQP native state overlaps.

    :param nstates:
        Number of physical MRSF states in every overlap calculation.
    :param target_state:
        One-based persistent state identity to follow.  The identity is
        initialized from the raw energy ordering at the first geometry.
    :param minimum_overlap:
        Minimum accepted absolute native overlap for the target assignment.
    :param maximum_ambiguity_ratio:
        Maximum accepted ratio of the second-largest target-row overlap to the
        assigned overlap.  Values near one indicate ambiguous character.
    :param failure_mode:
        ``"error"`` stops before evaluating a gradient when tracking is
        unreliable. ``"warn"`` continues with the globally optimal assignment
        and records a prominent warning.
    """

    method_name = 'openqp_mrsf_native_overlap'
    state_version = 1

    def __init__(self,
                 nstates,
                 target_state,
                 minimum_overlap=0.5,
                 maximum_ambiguity_ratio=0.8,
                 failure_mode='error',
                 guard_ground_state=True):

        self.nstates = int(nstates)
        self.target_state = int(target_state)
        self.minimum_overlap = float(minimum_overlap)
        self.maximum_ambiguity_ratio = float(maximum_ambiguity_ratio)
        self.failure_mode = str(failure_mode).strip().lower()
        self.guard_ground_state = bool(guard_ground_state)

        if self.nstates < 1:
            raise ValueError(
                'OpenQPMRSFStateTracker: nstates must be positive.')
        if not 1 <= self.target_state <= self.nstates:
            raise ValueError(
                'OpenQPMRSFStateTracker: target_state must be between 1 and '
                'nstates.')
        if not 0.0 <= self.minimum_overlap <= 1.0:
            raise ValueError(
                'OpenQPMRSFStateTracker: minimum_overlap must be in [0, 1].')
        if not 0.0 <= self.maximum_ambiguity_ratio <= 1.0:
            raise ValueError(
                'OpenQPMRSFStateTracker: maximum_ambiguity_ratio must be in '
                '[0, 1].')
        if self.failure_mode not in ('error', 'warn'):
            raise ValueError(
                'OpenQPMRSFStateTracker: failure_mode must be "error" or '
                '"warn".')

        self._previous_xyz = None
        self._previous_data = None
        self._tracked_to_raw = np.arange(self.nstates, dtype=int)
        self._step = -1
        self.history = []
        self._pending = None
        self._reference_revision = 0
        self.rejected_commits = 0

    @property
    def has_reference(self):
        """Whether a complete committed previous-geometry payload exists."""

        return (self._previous_xyz is not None and
                self._previous_data is not None)

    @property
    def has_pending(self):
        """Whether an evaluation is staged but not committed."""

        return self._pending is not None

    @property
    def reference_revision(self):
        """Counter that increments on every committed reference.

        Consumers that cache results by geometry must include this value in
        the cache key: the same coordinates evaluated against a different
        committed reference are a different calculation.
        """

        return self._reference_revision

    def committed_guess_payload(self):
        """Returns the committed OpenQP snapshot for SCF initialization.

        Only the committed reference is exposed.  A pending or rolled-back
        proposal must never seed an SCF, otherwise a rejected nuclear trial
        would still steer the orbitals of the next geometry.
        """

        if self._previous_data is None:
            return None

        return deepcopy(self._previous_data)

    @property
    def tracked_to_raw(self):
        """One-based persistent-to-raw root permutation."""

        return self._tracked_to_raw.copy() + 1

    @property
    def selected_raw_root(self):
        """One-based raw root currently carrying the target identity."""

        return int(self._tracked_to_raw[self.target_state - 1] + 1)

    def reset(self):
        """Discards the reference payload, any proposal, and the history."""

        self._previous_xyz = None
        self._previous_data = None
        self._tracked_to_raw = np.arange(self.nstates, dtype=int)
        self._step = -1
        self.history = []
        self._pending = None
        self._reference_revision = 0
        self.rejected_commits = 0

    def get_reference_payload(self):
        """Returns a defensive copy for ``OpenQP Molecule.back_door``."""

        if not self.has_reference:
            return None
        return deepcopy(self._previous_xyz), deepcopy(self._previous_data)

    def initialize_result(self, state_energies):
        """Defines persistent identities from the first raw root ordering."""

        energies = self._validate_energies(state_energies)
        identity = np.arange(self.nstates, dtype=int) + 1
        return OpenQPStateTrackingResult(
            step=0,
            target_tracked_root=self.target_state,
            previous_raw_root=self.target_state,
            selected_raw_root=self.target_state,
            tracked_to_raw=identity,
            assignment_scores=np.ones(self.nstates),
            target_score=1.0,
            second_score=0.0,
            ambiguity_ratio=0.0,
            global_confidence=1.0,
            is_reliable=True,
            reordered=False,
            root_changed=False,
            state_energies=energies,
            similarity_matrix=np.eye(self.nstates),
            signed_overlap_matrix=None,
            warnings=['OpenQP MRSF tracking reference initialized.'],
            reference_initialized=True)

    def evaluate(self, state_energies, state_overlap):
        """Assigns current raw roots to the previous persistent identities.

        ``state_overlap`` must have previous raw roots in rows and current raw
        roots in columns, as returned by ``OQP::td_states_overlap``.
        """

        if not self.has_reference:
            raise RuntimeError(
                'OpenQPMRSFStateTracker.evaluate: reference is not initialized.')

        energies = self._validate_energies(state_energies)
        signed = self._validate_overlap(state_overlap)
        # The metric itself is unchanged: absolute native overlap.  Magnitudes
        # above one are rejected in _validate_overlap rather than clipped.
        raw_similarity = np.abs(signed)

        # OpenQP rows are in the previous raw order.  Reorder those rows into
        # persistent-state order before solving the global assignment.
        similarity = raw_similarity[self._tracked_to_raw, :]
        assignment = solve_assignment(similarity)
        scores = similarity[np.arange(self.nstates), assignment]

        target_idx = self.target_state - 1
        previous_raw = int(self._tracked_to_raw[target_idx])
        selected_raw = int(assignment[target_idx])
        target_score = float(scores[target_idx])

        alternatives = np.delete(similarity[target_idx], selected_raw)
        second_score = (float(np.max(alternatives))
                        if alternatives.size else 0.0)
        ambiguity_ratio = (second_score / target_score
                           if target_score > np.finfo(float).eps else np.inf)
        global_confidence = float(np.min(scores))

        reliable_overlap = target_score >= self.minimum_overlap
        reliable_ratio = (
            ambiguity_ratio <= self.maximum_ambiguity_ratio)
        is_reliable = bool(reliable_overlap and reliable_ratio)

        warnings = []
        if not reliable_overlap:
            warnings.append(
                f'target overlap {target_score:.6f} is below the configured '
                f'minimum {self.minimum_overlap:.6f}')
        if not reliable_ratio:
            warnings.append(
                f'target ambiguity ratio {ambiguity_ratio:.6f} exceeds the '
                f'configured maximum {self.maximum_ambiguity_ratio:.6f}')

        one_based_assignment = assignment + 1
        reordered = bool(
            np.any(one_based_assignment !=
                   np.arange(self.nstates, dtype=int) + 1))
        root_changed = selected_raw != previous_raw
        if root_changed:
            warnings.append(
                f'persistent root {self.target_state} moved from raw root '
                f'{previous_raw + 1} to {selected_raw + 1}')

        # An MRSF spectrum carries S0 as an ordinary root, so a tracked
        # excited identity can be assigned onto the ground state at a conical
        # intersection and then quietly minimize S0 instead.  The ground state
        # is the lowest-energy raw root at this geometry; the raw order is not
        # fixed along a path, so it is located by energy.
        ground_raw = int(np.argmin(energies))
        ground_state_collision = bool(
            self.guard_ground_state and
            self.target_state != 1 and
            selected_raw == ground_raw)
        if ground_state_collision:
            is_reliable = False
            warnings.append(
                f'persistent root {self.target_state} was assigned to raw '
                f'root {ground_raw + 1}, which is the ground state at this '
                'geometry; the tracked excited state has collapsed onto S0')

        return OpenQPStateTrackingResult(
            step=self._step + 1,
            target_tracked_root=self.target_state,
            previous_raw_root=previous_raw + 1,
            selected_raw_root=selected_raw + 1,
            tracked_to_raw=one_based_assignment,
            assignment_scores=scores,
            target_score=target_score,
            second_score=second_score,
            ambiguity_ratio=float(ambiguity_ratio),
            global_confidence=global_confidence,
            is_reliable=is_reliable,
            reordered=reordered,
            root_changed=root_changed,
            state_energies=energies,
            similarity_matrix=similarity,
            signed_overlap_matrix=signed,
            warnings=warnings,
            ground_state_raw_root=ground_raw + 1,
            ground_state_collision=ground_state_collision)

    def propose(self, result, xyz, data):
        """Stages a finished evaluation as a *candidate* reference.

        The committed reference is untouched, so a nuclear trial that the
        optimizer later rejects cannot become the comparison point for the
        next geometry.  Call :meth:`commit` or :meth:`rollback` afterwards.

        :param result:
            The :class:`OpenQPStateTrackingResult` for this geometry.
        :param xyz:
            Geometry that produced ``data``, in bohr.
        :param data:
            OpenQP electronic snapshot taken before the gradient step.
        """

        if not isinstance(result, OpenQPStateTrackingResult):
            raise TypeError(
                'OpenQPMRSFStateTracker.propose: invalid tracking result.')
        if result.tracked_to_raw.size != self.nstates:
            raise ValueError(
                'OpenQPMRSFStateTracker.propose: invalid root permutation.')

        snapshot = deepcopy(data)
        self._validate_snapshot(snapshot)

        coordinates = np.asarray(xyz, dtype=float)
        if coordinates.size == 0 or not np.all(np.isfinite(coordinates)):
            raise ValueError(
                'OpenQPMRSFStateTracker.propose: invalid geometry snapshot.')

        self._pending = {
            'result': result,
            'xyz': coordinates.copy(),
            'data': snapshot,
        }

        return result

    def commit(self):
        """Promotes a staged proposal to the committed reference.

        Committing is refused for an unreliable assignment: an ambiguous or
        low-overlap decision must not become the identity that every later
        geometry is measured against.

        :return:
            ``True`` when the committed reference advanced.
        """

        if self._pending is None:
            # Nothing was evaluated since the last commit (for example a
            # cached optimizer step).  Committing must stay idempotent.
            return False

        result = self._pending['result']
        if not bool(result.is_reliable):
            self._pending = None
            self.rejected_commits += 1
            return False

        self._previous_xyz = self._pending['xyz']
        self._previous_data = self._pending['data']
        self._tracked_to_raw = (
            np.asarray(result.tracked_to_raw, dtype=int) - 1)
        self._step = int(result.step)
        self.history.append(result)
        self._reference_revision += 1
        self._pending = None

        return True

    def rollback(self):
        """Discards a staged proposal, keeping the committed reference.

        :return:
            ``True`` when a proposal was actually discarded.
        """

        had_pending = self._pending is not None
        self._pending = None

        return had_pending

    def accept(self, result, xyz, data):
        """Proposes and immediately commits, for non-transactional callers.

        Molecular dynamics has no rejected-step concept: every evaluated
        geometry is on the trajectory, so propose/commit collapse into one
        call.  Returns ``True`` when the reference advanced.
        """

        self.propose(result, xyz, data)

        return self.commit()

    def state_dict(self):
        """Returns complete restart state as plain Python containers."""

        return {
            'version': self.state_version,
            'method': self.method_name,
            'nstates': self.nstates,
            'target_state': self.target_state,
            'minimum_overlap': self.minimum_overlap,
            'maximum_ambiguity_ratio': self.maximum_ambiguity_ratio,
            'failure_mode': self.failure_mode,
            'step': self._step,
            # Only committed state is serialized; a pending proposal belongs
            # to an optimizer step that has not been accepted yet.
            'reference_revision': self._reference_revision,
            'tracked_to_raw': self.tracked_to_raw.tolist(),
            'previous_xyz':
                (None if self._previous_xyz is None else
                 self._previous_xyz.tolist()),
            'previous_data': _json_compatible(self._previous_data),
            'history': [
                deepcopy(result) if isinstance(result, dict) else
                result.to_dict()
                for result in self.history
            ],
        }

    def save_state(self, filename):
        """Writes a complete JSON restart payload."""

        with open(filename, 'w', encoding='utf-8') as state_file:
            json.dump(self.state_dict(), state_file, indent=2)

    def load_state(self, filename):
        """Restores a restart payload written by :meth:`save_state`."""

        with open(filename, encoding='utf-8') as state_file:
            state = json.load(state_file)
        self.load_state_dict(state)

    def load_state_dict(self, state):
        """Restores and validates a plain restart dictionary."""

        if int(state.get('version', -1)) != self.state_version:
            raise ValueError(
                'OpenQPMRSFStateTracker: unsupported restart-state version.')
        if state.get('method') != self.method_name:
            raise ValueError(
                'OpenQPMRSFStateTracker: restart method does not match.')
        if int(state.get('nstates', -1)) != self.nstates:
            raise ValueError(
                'OpenQPMRSFStateTracker: restart nstates does not match.')
        if int(state.get('target_state', -1)) != self.target_state:
            raise ValueError(
                'OpenQPMRSFStateTracker: restart target_state does not match.')

        permutation = np.asarray(state.get('tracked_to_raw'), dtype=int)
        expected = np.arange(self.nstates, dtype=int) + 1
        if (permutation.shape != (self.nstates,) or
                not np.array_equal(np.sort(permutation), expected)):
            raise ValueError(
                'OpenQPMRSFStateTracker: restart permutation is invalid.')

        xyz = state.get('previous_xyz')
        data = state.get('previous_data')
        if (xyz is None) != (data is None):
            raise ValueError(
                'OpenQPMRSFStateTracker: incomplete restart reference.')
        if data is not None:
            self._validate_snapshot(data)
            xyz = np.asarray(xyz, dtype=float)
            if xyz.size == 0 or not np.all(np.isfinite(xyz)):
                raise ValueError(
                    'OpenQPMRSFStateTracker: invalid restart geometry.')

        self._tracked_to_raw = permutation - 1
        self._step = int(state.get('step', -1))
        self._previous_xyz = None if xyz is None else xyz.copy()
        self._previous_data = deepcopy(data)
        # Payloads written before the transactional API carry no revision.
        self._reference_revision = int(state.get('reference_revision', 0))
        self._pending = None
        self.rejected_commits = 0
        # Historic diagnostics are intentionally retained as dictionaries.
        # New in-process entries remain strongly typed result objects.
        self.history = list(state.get('history', []))

    def _validate_energies(self, state_energies):
        energies = np.asarray(state_energies, dtype=float).reshape(-1)
        if energies.shape != (self.nstates,):
            raise ValueError(
                'OpenQPMRSFStateTracker: expected '
                f'{self.nstates} state energies, got {energies.size}.')
        if not np.all(np.isfinite(energies)):
            raise ValueError(
                'OpenQPMRSFStateTracker: nonfinite state energies.')
        return energies.copy()

    def _validate_overlap(self, state_overlap):
        overlap = np.asarray(state_overlap, dtype=float)
        if overlap.size != self.nstates * self.nstates:
            raise ValueError(
                'OpenQPMRSFStateTracker: native overlap has size '
                f'{overlap.size}, expected {self.nstates * self.nstates}.')
        overlap = overlap.reshape((self.nstates, self.nstates))
        if not np.all(np.isfinite(overlap)):
            raise ValueError(
                'OpenQPMRSFStateTracker: nonfinite native state overlap.')

        # Never clip: a magnitude above one means the overlap is not a
        # normalized state overlap, and hiding it would silently corrupt the
        # assignment instead of reporting a broken native quantity.
        largest = float(np.max(np.abs(overlap)))
        if largest > 1.0 + _OVERLAP_MAGNITUDE_TOLERANCE:
            raise ValueError(
                'OpenQPMRSFStateTracker: native state overlap has magnitude '
                f'{largest:.6f}, which exceeds one by more than the accepted '
                f'tolerance {_OVERLAP_MAGNITUDE_TOLERANCE:.1e}. The overlap '
                'is not a normalized state overlap.')

        return overlap.copy()

    def _validate_snapshot(self, data):
        if not isinstance(data, dict):
            raise TypeError(
                'OpenQPMRSFStateTracker: OpenQP snapshot must be a dictionary.')
        missing = [key for key in _REQUIRED_SNAPSHOT_KEYS if key not in data]
        if missing:
            raise ValueError(
                'OpenQPMRSFStateTracker: OpenQP snapshot is missing ' +
                ', '.join(missing) + '.')

        td_energies = np.asarray(data['OQP::td_energies']).reshape(-1)
        if td_energies.size < self.nstates:
            raise ValueError(
                'OpenQPMRSFStateTracker: snapshot contains fewer TD states '
                'than requested.')


def _json_compatible(value):
    """Recursively converts NumPy containers for a JSON restart payload."""

    if value is None:
        return None
    if isinstance(value, dict):
        return {str(key): _json_compatible(val)
                for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_compatible(val) for val in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value
