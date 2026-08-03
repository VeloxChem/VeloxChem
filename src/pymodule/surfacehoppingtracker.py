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
"""
Electronic-state tracking for surface-hopping dynamics.

The root index returned by a TDA/TDDFT eigensolver is an energy ordering, not
a persistent physical state identity: two roots that cross between consecutive
geometries swap indices without anything physical happening.  Feeding raw root
indices into a Landau-Zener gap history therefore produces spurious events.

This module assigns a persistent identity to each root by maximizing a
state-similarity matrix

.. math::

    O_{ij} = \\left| \\langle \\Psi_i(t) | \\Psi_j(t+\\Delta t) \\rangle
             \\right|^2

or a documented approximation to it, subject to a one-to-one constraint solved
with the Hungarian algorithm.

The module is deliberately free of VeloxChem driver imports so that trackers
can be unit tested against synthetic descriptors.
"""

from dataclasses import dataclass, field
import numpy as np

try:
    from scipy.optimize import linear_sum_assignment
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False


@dataclass
class StateDescriptors:
    """
    Quantities describing the character of each tracked electronic state at
    one geometry.

    Only ``excitation_energies`` is mandatory.  Trackers use whichever of the
    optional fields are populated by the electronic-structure provider.

    :param excitation_energies:
        Array of shape (n_states,), in Hartree.  Index 0 corresponds to the
        ground state and is conventionally 0.0.
    :param oscillator_strengths:
        Array of shape (n_states,), dimensionless.
    :param transition_dipoles:
        Array of shape (n_states, 3), in atomic units.
    :param excitation_vectors:
        Array of shape (n_states, n_amplitudes) with response vectors in a
        common basis.  Only meaningful when the molecular-orbital bases at the
        two geometries have been aligned.
    :param overlap_matrix:
        Optional precomputed |<Psi_i|Psi_j>|^2 matrix supplied directly by a
        rigorous provider, of shape (n_states, n_states).
    """

    excitation_energies: np.ndarray
    oscillator_strengths: np.ndarray = None
    transition_dipoles: np.ndarray = None
    excitation_vectors: np.ndarray = None
    overlap_matrix: np.ndarray = None

    @property
    def n_states(self):
        """
        :return:
            The number of tracked states described.
        """

        return int(np.asarray(self.excitation_energies).size)


@dataclass
class TrackingResult:
    """
    Outcome of one state-tracking step.

    :param permutation:
        Array of shape (n_states,) mapping tracked (persistent) index to the
        raw root index at the current geometry.
    :param scores:
        Array of shape (n_states,) with the similarity score assigned to each
        tracked state.
    :param confidence:
        Overall tracking confidence in [0, 1]; the minimum assignment score.
    :param is_reliable:
        Whether the confidence meets the configured threshold.
    :param swapped:
        Whether the assignment differs from the identity permutation.
    :param method:
        Name of the tracker that produced the result.
    :param warnings:
        List of diagnostic warning strings.
    :param similarity_matrix:
        Full matrix used for the one-to-one assignment.  Rows are persistent
        tracked states at the previous geometry and columns are raw roots at
        the current geometry.
    """

    permutation: np.ndarray
    scores: np.ndarray
    confidence: float
    is_reliable: bool
    swapped: bool
    method: str
    warnings: list = field(default_factory=list)
    similarity_matrix: np.ndarray = None


def solve_assignment(similarity):
    """
    Solves the maximum-weight one-to-one assignment problem

    .. math:: \\pi = \\arg\\max_\\pi \\sum_i O_{i, \\pi(i)}

    using the Hungarian algorithm when SciPy is available and a greedy
    fallback otherwise.

    :param similarity:
        Square similarity matrix of shape (n, n); rows are previous states,
        columns are current roots.

    :return:
        Array of shape (n,) with the assigned column index for each row.
    """

    similarity = np.asarray(similarity, dtype=float)

    if similarity.ndim != 2 or similarity.shape[0] != similarity.shape[1]:
        raise ValueError('solve_assignment: similarity must be square.')

    n = similarity.shape[0]

    if not np.all(np.isfinite(similarity)):
        raise ValueError('solve_assignment: nonfinite similarity matrix.')

    if _SCIPY_AVAILABLE:
        rows, cols = linear_sum_assignment(-similarity)
        assignment = np.empty(n, dtype=int)
        assignment[rows] = cols
        return assignment

    # Greedy fallback: repeatedly take the largest remaining entry.  Exact for
    # well-separated states, which is the regime where tracking is easy
    # anyway; ambiguous cases are caught by the confidence threshold.
    assignment = -np.ones(n, dtype=int)
    remaining_rows = set(range(n))
    remaining_cols = set(range(n))

    while remaining_rows:
        best = None
        for i in remaining_rows:
            for j in remaining_cols:
                if best is None or similarity[i, j] > similarity[best]:
                    best = (i, j)
        assignment[best[0]] = best[1]
        remaining_rows.discard(best[0])
        remaining_cols.discard(best[1])

    return assignment


class ElectronicStateTracker:
    """
    Base class for electronic-state trackers.

    A tracker maps the raw root ordering at each new geometry onto persistent
    tracked identities established at the previous geometry.

    :param minimum_overlap:
        Minimum accepted per-state similarity for the assignment to be
        considered reliable.
    """

    method_name = 'base'

    def __init__(self, minimum_overlap=0.5):

        self.minimum_overlap = float(minimum_overlap)
        self._reference = None

    @property
    def has_reference(self):
        """
        :return:
            Whether a reference geometry has been recorded.
        """

        return self._reference is not None

    def reset(self):
        """
        Discards the stored reference descriptors.
        """

        self._reference = None

    def similarity_matrix(self, previous, current):
        """
        Builds the state-similarity matrix between two sets of descriptors.

        :param previous:
            :class:`StateDescriptors` at the previous geometry.
        :param current:
            :class:`StateDescriptors` at the current geometry.

        :return:
            Array of shape (n_states, n_states) with entries in [0, 1].
        """

        raise NotImplementedError

    def track(self, descriptors):
        """
        Assigns persistent identities to the roots at the current geometry.

        On the first call the current ordering defines the reference and the
        identity permutation is returned.

        :param descriptors:
            :class:`StateDescriptors` at the current geometry.

        :return:
            A :class:`TrackingResult`.
        """

        n_states = descriptors.n_states

        if not self.has_reference:
            self._reference = descriptors
            return TrackingResult(permutation=np.arange(n_states),
                                  scores=np.ones(n_states),
                                  confidence=1.0,
                                  is_reliable=True,
                                  swapped=False,
                                  method=self.method_name,
                                  warnings=['tracking reference initialized'],
                                  similarity_matrix=np.eye(n_states))

        if self._reference.n_states != n_states:
            raise ValueError(
                'ElectronicStateTracker: number of tracked states changed '
                f'from {self._reference.n_states} to {n_states}.')

        similarity = self.similarity_matrix(self._reference, descriptors)
        permutation = solve_assignment(similarity)
        scores = similarity[np.arange(n_states), permutation]

        confidence = float(np.min(scores)) if n_states else 1.0
        is_reliable = confidence >= self.minimum_overlap
        swapped = bool(np.any(permutation != np.arange(n_states)))

        warnings = []
        if not is_reliable:
            warnings.append(
                f'state tracking confidence {confidence:.3f} below threshold '
                f'{self.minimum_overlap:.3f}')
        if swapped:
            warnings.append(f'root reordering detected: {permutation.tolist()}')

        # The reference is stored in *tracked* order, so that identities
        # propagate across an arbitrary number of reorderings.
        self._reference = _permute_descriptors(descriptors, permutation)

        return TrackingResult(permutation=permutation,
                              scores=scores,
                              confidence=confidence,
                              is_reliable=is_reliable,
                              swapped=swapped,
                              method=self.method_name,
                              warnings=warnings,
                              similarity_matrix=similarity)


class OverlapStateTracker(ElectronicStateTracker):
    """
    Rigorous tracker driven by a wavefunction-overlap matrix supplied by the
    electronic-structure provider.

    The provider must populate :attr:`StateDescriptors.overlap_matrix` with
    ``|<Psi_i(t)|Psi_j(t+dt)>|^2`` evaluated with properly aligned molecular
    orbital bases, or populate :attr:`StateDescriptors.excitation_vectors`
    with response vectors already expressed in a common basis.
    """

    method_name = 'transition_density'

    def similarity_matrix(self, previous, current):
        """
        See :func:`ElectronicStateTracker.similarity_matrix`.
        """

        if current.overlap_matrix is not None:
            overlap = np.abs(np.asarray(current.overlap_matrix, dtype=float))
            return np.clip(overlap, 0.0, 1.0)

        if (previous.excitation_vectors is not None and
                current.excitation_vectors is not None):
            prev_vec = _normalize_rows(
                np.asarray(previous.excitation_vectors, dtype=float))
            curr_vec = _normalize_rows(
                np.asarray(current.excitation_vectors, dtype=float))
            return np.clip(np.abs(prev_vec @ curr_vec.T)**2, 0.0, 1.0)

        raise ValueError(
            'OverlapStateTracker: neither an overlap matrix nor aligned '
            'excitation vectors were provided. Use the descriptor tracker '
            'instead, or supply rigorous overlaps.')


class DescriptorStateTracker(ElectronicStateTracker):
    """
    Fallback tracker based on continuity of state descriptors.

    .. warning::

        This tracker is **not** a wavefunction overlap.  It combines
        continuity of excitation energy, oscillator strength and transition
        dipole direction into a heuristic similarity score.  It is reliable
        for well-separated states with distinct character and degrades
        precisely where nonadiabatic dynamics is most interesting, namely near
        conical intersections.  The reported confidence must be monitored, and
        the tracker should be replaced by :class:`OverlapStateTracker` as soon
        as rigorous overlaps are available for the chosen electronic-structure
        backend.

    The similarity between previous state ``i`` and current root ``j`` is

    .. math::

        O_{ij} = w_E s^E_{ij} + w_f s^f_{ij} + w_\\mu s^\\mu_{ij}

    with Lorentzian energy and oscillator-strength similarities and a squared
    transition-dipole direction cosine.  Weights of unavailable descriptors
    are redistributed over the available ones.

    :param minimum_overlap:
        Minimum accepted similarity for a reliable assignment.
    :param energy_scale:
        Lorentzian width for the excitation-energy similarity, in Hartree.
    :param oscillator_scale:
        Lorentzian width for the oscillator-strength similarity.
    """

    method_name = 'descriptor'

    def __init__(self,
                 minimum_overlap=0.5,
                 energy_scale=0.02,
                 oscillator_scale=0.1):

        super().__init__(minimum_overlap=minimum_overlap)

        self.energy_scale = float(energy_scale)
        self.oscillator_scale = float(oscillator_scale)

        if self.energy_scale <= 0.0 or self.oscillator_scale <= 0.0:
            raise ValueError(
                'DescriptorStateTracker: similarity scales must be positive.')

    def similarity_matrix(self, previous, current):
        """
        See :func:`ElectronicStateTracker.similarity_matrix`.
        """

        e_prev = np.asarray(previous.excitation_energies, dtype=float)
        e_curr = np.asarray(current.excitation_energies, dtype=float)
        n = e_prev.size

        ones = np.ones((n, n))

        # Each component carries a per-element validity mask.  Where a
        # descriptor is undefined for a given pair - most importantly the null
        # transition dipole of the ground state - its weight is redistributed
        # over the remaining components instead of scoring the pair as
        # dissimilar.
        components = []

        # Excitation-energy continuity (always available).
        d_energy = np.abs(e_prev[:, None] - e_curr[None, :])
        components.append(
            (2.0, 1.0 / (1.0 + (d_energy / self.energy_scale)**2), ones))

        # Oscillator-strength continuity.
        if (previous.oscillator_strengths is not None and
                current.oscillator_strengths is not None):
            f_prev = np.asarray(previous.oscillator_strengths, dtype=float)
            f_curr = np.asarray(current.oscillator_strengths, dtype=float)
            d_osc = np.abs(f_prev[:, None] - f_curr[None, :])
            components.append(
                (1.0, 1.0 / (1.0 + (d_osc / self.oscillator_scale)**2), ones))

        # Transition-dipole direction continuity.
        if (previous.transition_dipoles is not None and
                current.transition_dipoles is not None):
            mu_prev = np.asarray(previous.transition_dipoles, dtype=float)
            mu_curr = np.asarray(current.transition_dipoles, dtype=float)
            norm_prev = np.linalg.norm(mu_prev, axis=1)
            norm_curr = np.linalg.norm(mu_curr, axis=1)
            cosine = _normalize_rows(mu_prev) @ _normalize_rows(mu_curr).T
            mask = np.outer(norm_prev > 0.0, norm_curr > 0.0).astype(float)
            components.append((1.5, np.clip(cosine**2, 0.0, 1.0), mask))

        numerator = np.zeros((n, n))
        denominator = np.zeros((n, n))

        for weight, component, mask in components:
            numerator += weight * mask * component
            denominator += weight * mask

        # denominator is never zero: the energy component is always valid.
        similarity = numerator / denominator

        return np.clip(similarity, 0.0, 1.0)


class IdentityStateTracker(ElectronicStateTracker):
    """
    Non-tracker that keeps the raw energy ordering.

    Provided only for reference calculations and for reproducing the legacy
    fixed-root behaviour.  It cannot detect root reordering and will therefore
    generate spurious Landau-Zener events at state crossings.
    """

    method_name = 'none'

    def similarity_matrix(self, previous, current):
        """
        See :func:`ElectronicStateTracker.similarity_matrix`.
        """

        return np.eye(previous.n_states)


TRACKERS = {
    'transition_density': OverlapStateTracker,
    'descriptor': DescriptorStateTracker,
    'none': IdentityStateTracker,
}


def create_state_tracker(method, minimum_overlap=0.5, **kwargs):
    """
    Instantiates a state tracker by name.

    :param method:
        One of the keys of :data:`TRACKERS`.
    :param minimum_overlap:
        Minimum accepted tracking confidence.

    :return:
        An :class:`ElectronicStateTracker` instance.
    """

    key = str(method).lower()

    if key not in TRACKERS:
        valid = ', '.join(sorted(TRACKERS))
        raise ValueError(
            f'create_state_tracker: unknown tracker {method}; valid: {valid}')

    return TRACKERS[key](minimum_overlap=minimum_overlap, **kwargs)


def _normalize_rows(matrix):
    """
    Normalizes each row of a matrix to unit length, leaving null rows
    untouched.

    :param matrix:
        Array of shape (n, m).

    :return:
        The row-normalized array.
    """

    matrix = np.asarray(matrix, dtype=float)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms = np.where(norms > 0.0, norms, 1.0)

    return matrix / norms


def _permute_descriptors(descriptors, permutation):
    """
    Reorders descriptors into tracked order.

    :param descriptors:
        :class:`StateDescriptors` in raw root order.
    :param permutation:
        Array mapping tracked index to raw root index.

    :return:
        A new :class:`StateDescriptors` in tracked order.
    """

    permutation = np.asarray(permutation, dtype=int)

    def take(array):
        if array is None:
            return None
        return np.asarray(array)[permutation]

    return StateDescriptors(
        excitation_energies=take(descriptors.excitation_energies),
        oscillator_strengths=take(descriptors.oscillator_strengths),
        transition_dipoles=take(descriptors.transition_dipoles),
        excitation_vectors=take(descriptors.excitation_vectors),
        overlap_matrix=None)
