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
Backend adapters that give surface hopping *target-state* semantics.

Conventional linear-response TDA/TDDFT answers a question that neither
MRSF-TDDFT nor SF-TDDFT answers: "what is the ground-state energy, and what
are the excitations above it?".  Both spin-flip families instead produce a
*working reference* of a different spin manifold together with a set of
physical target states that may lie on either side of it.  Forcing them
through the conventional contract is what produced the two defects this
module repairs:

* the OpenQP MRSF driver exposes a **scalar** excitation energy of the
  selected root, so a multi-state ladder could not be built from it at all;
* the Serenity SF ladder was built by prepending the high-spin working
  reference, exposing a state of the wrong spin manifold as dynamics state 0.

The repair is a single normalized contract, :class:`ElectronicSnapshot`,
which carries **only target-state total energies** to the hopping controller.
The working reference survives as provenance and diagnostics and is never a
dynamics surface.

Index spaces
------------

``raw target index r``
    ``0 .. N-1``.  The backend's own ordering of the *physical target states*
    at one geometry.  The working reference is not part of this space.

``backend derivative selector``
    What the backend must be told to differentiate raw target ``r``.  Both
    installed backends use ``r + 1``, but the arithmetic is performed once,
    inside the adapter, and is carried explicitly on the snapshot in
    :attr:`ElectronicSnapshot.derivative_selectors`.

``tracked dynamics state t``
    Persistent physical identity owned by the surface-hopping controller and
    its state tracker.  Adapters never see it.

Measured backend contracts
--------------------------

OpenQP 1.2.0, MRSF-TDDFT, ``nstate = N``::

    len(mol.energies) == N + 1
    mol.energies[0]       high-spin ROHF working reference (diagnostics only)
    mol.energies[1:N+1]   the N total MRSF target-state energies
    properties.grad = [r + 1]     gradient of raw target r

Serenity 1.7, SF-TDA/SF-TDDFT, ``nEigen = N``::

    E_reference               high-spin SCF energy, system.getEnergy()
    omega[r]                  SF response energy of raw target r, r = 0..N-1
    E_target[r]               E_reference + omega[r]
    excGradList = [r + 1]     gradient of raw target r

The Serenity gradient was verified by central finite differences to be the
**total** target-state gradient ``d(E_reference + omega)/dR``; no reference
gradient is added by the adapter.  See
``tests/test_surface_hopping_backend_gates.py``.
"""

from dataclasses import dataclass, field
import hashlib
import json
import os
import shutil
import uuid

import numpy as np

from .errorhandler import assert_msg_critical
from .molecule import Molecule


#: Revision of the normalized snapshot contract.  Any change of the meaning
#: of a snapshot field must increment this, because it participates in every
#: cache key and in the checkpoint payload.
SNAPSHOT_CONTRACT_REVISION = '1'

#: Number of decimals used to canonicalize a geometry into a fingerprint.
#: 1e-10 bohr is far below any meaningful nuclear displacement and far above
#: the noise of the OpenMM position round trip.
GEOMETRY_FINGERPRINT_DECIMALS = 10


class BackendCapabilityError(RuntimeError):
    """
    Raised when a backend cannot satisfy a production requirement.

    Distinct from a numerical failure: it means the installed backend, its
    version, or the requested method does not provide something surface
    hopping needs (multi-state target energies, a state-specific gradient, or
    a cross-geometry electronic descriptor).
    """


def geometry_fingerprint(geometry):
    """
    Canonicalizes a geometry into a stable content fingerprint.

    Atom order is part of the fingerprint by construction: the array is
    hashed in row-major order and never sorted.

    :param geometry:
        Array of shape (n_atoms, 3) in bohr.

    :return:
        A hexadecimal digest string.
    """

    coordinates = np.asarray(geometry, dtype=float)

    assert_msg_critical(
        coordinates.ndim == 2 and coordinates.shape[1] == 3,
        'geometry_fingerprint: expected an (n_atoms, 3) geometry.')
    assert_msg_critical(
        np.all(np.isfinite(coordinates)),
        'geometry_fingerprint: nonfinite geometry.')

    rounded = np.round(coordinates, GEOMETRY_FINGERPRINT_DECIMALS)
    # +0.0 collapses the two signed zeros so that -0.0 and 0.0 agree.
    rounded = rounded + 0.0

    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(rounded, dtype='<f8').tobytes())

    return digest.hexdigest()


def settings_fingerprint(payload):
    """
    Canonicalizes every result-affecting setting into a stable fingerprint.

    :param payload:
        A JSON-serializable mapping.  Keys are sorted, so the caller need not
        care about insertion order.

    :return:
        A hexadecimal digest string.
    """

    text = json.dumps(payload, sort_keys=True, separators=(',', ':'),
                      default=str)

    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def _frozen(array, dtype=float):
    """
    Returns a read-only copy of an array, so a snapshot cannot be mutated by
    whoever received it.
    """

    if array is None:
        return None

    frozen = np.array(array, dtype=dtype, copy=True)
    frozen.setflags(write=False)

    return frozen


def assignment_similarity(raw_overlap):
    """
    The single documented assignment metric shared by every backend and by
    both the geometry-optimization and the surface-hopping entry points.

    The metric is the **absolute** state overlap ``|<Psi_i|Psi_j>|``, not its
    square.  Taking the modulus is what makes the assignment invariant under
    the arbitrary sign (phase) of a response vector: a root whose eigenvector
    is returned with the opposite sign at the next geometry is still the same
    physical state.

    Magnitudes above one are refused rather than clipped: a normalized state
    overlap cannot exceed one, so a larger value means the supplied quantity
    is not the assumed one and silently clipping it would corrupt the
    assignment.

    :param raw_overlap:
        Square array with previous raw targets in rows and current raw
        targets in columns.

    :return:
        The similarity matrix used for the one-to-one assignment.
    """

    overlap = np.asarray(raw_overlap, dtype=float)

    if overlap.ndim != 2 or overlap.shape[0] != overlap.shape[1]:
        raise ValueError(
            'assignment_similarity: expected a square (previous, current) '
            f'overlap, got shape {overlap.shape}.')
    if not np.all(np.isfinite(overlap)):
        raise ValueError('assignment_similarity: nonfinite state overlap.')

    similarity = np.abs(overlap)
    largest = float(np.max(similarity)) if similarity.size else 0.0

    if largest > 1.0 + 1.0e-3:
        raise ValueError(
            'assignment_similarity: state overlap has magnitude '
            f'{largest:.6f}, which exceeds one by more than the accepted '
            'tolerance 1.0e-03; the supplied quantity is not a normalized '
            'state overlap.')

    return similarity


@dataclass(frozen=True)
class ElectronicSnapshot:
    """
    Immutable, backend-neutral electronic result at one nuclear geometry.

    Only :attr:`target_energies` reaches the hopping controller.  The working
    reference is retained in :attr:`reference_energy` for provenance and is
    never a dynamics surface.

    :param backend:
        Backend identifier, e.g. ``'openqp'`` or ``'serenity'``.
    :param method:
        Electronic method, e.g. ``'mrsf-tddft'`` or ``'sf-tda'``.
    :param backend_version:
        Installed backend version string.
    :param adapter_revision:
        Revision of the adapter that produced this snapshot.
    :param calculation_id:
        Unique identifier of the backend calculation that produced it.  Used
        to key gradients and native payloads; never reused.
    :param geometry:
        Read-only array of shape (n_atoms, 3) in bohr.
    :param atom_labels:
        Tuple of element labels in the same order as ``geometry``.
    :param geometry_fingerprint:
        Content fingerprint of ``geometry``.
    :param settings_fingerprint:
        Content fingerprint of every result-affecting setting.
    :param charge:
        Molecular charge.
    :param reference_multiplicity:
        Spin multiplicity of the working reference (3 for both installed
        spin-flip backends).
    :param target_manifold:
        Spin manifold requested for the target states.
    :param functional:
        Exchange-correlation functional label.
    :param basis:
        Basis-set label.
    :param target_energies:
        Read-only array of shape (n_states,) with **total** target-state
        energies in Hartree, in raw target order.
    :param derivative_selectors:
        Tuple of length n_states; entry ``r`` is what the backend must be
        told to differentiate raw target ``r``.
    :param reference_energy:
        Working-reference total energy in Hartree.  Provenance only.
    :param response_energies:
        Read-only array of shape (n_states,) with the backend response
        energies relative to the working reference.  Provenance only.
    :param overlap_to_previous:
        Read-only (n_prev, n_states) array with previous raw targets in rows
        and current raw targets in columns, or ``None`` at the first
        geometry.
    :param previous_calculation_id:
        Identifier of the snapshot ``overlap_to_previous`` was measured
        against.  Reversing the pair is a different quantity.
    :param overlap_source:
        Name of the conditioned descriptor selected for assignment.
    :param state_overlap_method:
        Backend algorithm that produced the native cross-geometry overlap.
    :param native_overlap_matrix:
        Signed native OpenQP MRSF overlap retained for provenance.
    :param response_overlap_matrix:
        Signed normalized MRSF response-vector overlap, when available.
    :param native_spectral_norm:
        Largest singular value of the native overlap.
    :param selected_spectral_norm:
        Largest singular value of :attr:`overlap_to_previous`.
    :param overlap_is_conditioned:
        Whether the selected overlap satisfies the shared OpenQP conditioning
        criterion.
    :param spin_metadata:
        Backend spin diagnostics per raw target (``<S^2>``, inferred
        multiplicity), or ``None``.
    :param scf_converged:
        Whether the reference SCF converged.
    :param response_converged:
        Whether the response calculation converged.
    :param reference_stability:
        Structured diagnostics for a requested OpenQP reference-stability
        check on this snapshot, or ``None``.
    :param reference_continuity:
        Compact occupied-subspace singular-value diagnostics comparing this
        high-spin reference with the previous accepted reference.
    :param provider_generation:
        Generation counter of the owning provider; bumped on cache clear so
        that a stale entry can never be mistaken for a current one.
    """

    backend: str
    method: str
    backend_version: str
    adapter_revision: str
    calculation_id: str

    geometry: np.ndarray
    atom_labels: tuple

    geometry_fingerprint: str
    settings_fingerprint: str

    charge: float
    reference_multiplicity: int
    target_manifold: str
    functional: str
    basis: str

    target_energies: np.ndarray
    derivative_selectors: tuple

    reference_energy: float
    response_energies: np.ndarray

    overlap_to_previous: np.ndarray = None
    previous_calculation_id: str = None
    overlap_source: str = None
    state_overlap_method: str = None
    native_overlap_matrix: np.ndarray = None
    response_overlap_matrix: np.ndarray = None
    native_spectral_norm: float = None
    selected_spectral_norm: float = None
    overlap_is_conditioned: bool = None

    spin_metadata: dict = None
    scf_converged: bool = True
    response_converged: bool = True
    reference_stability: dict = None
    reference_continuity: dict = None

    energy_unit: str = 'hartree'
    coordinate_unit: str = 'bohr'
    gradient_unit: str = 'hartree/bohr'

    contract_revision: str = SNAPSHOT_CONTRACT_REVISION
    provider_generation: int = 0
    warnings: tuple = field(default_factory=tuple)

    @property
    def n_states(self):
        """
        :return:
            Number of physical target states.
        """

        return int(np.asarray(self.target_energies).size)

    @property
    def raw_target_indices(self):
        """
        :return:
            Tuple ``(0, ..., n_states - 1)``.
        """

        return tuple(range(self.n_states))

    def selector_for(self, raw_target):
        """
        Maps a raw target index to its backend derivative selector.

        This is the only place where a caller may obtain a selector; no
        surface-hopping module performs selector arithmetic of its own.

        :param raw_target:
            Raw target index in ``0 .. n_states - 1``.

        :return:
            The backend derivative selector.
        """

        index = int(raw_target)

        assert_msg_critical(
            0 <= index < self.n_states,
            f'ElectronicSnapshot: raw target index {index} is outside the '
            f'{self.n_states} computed target states.')

        return int(self.derivative_selectors[index])

    def is_valid(self):
        """
        :return:
            True when the snapshot is complete, converged and finite.
        """

        if not (self.scf_converged and self.response_converged):
            return False
        if not np.all(np.isfinite(np.asarray(self.target_energies,
                                             dtype=float))):
            return False
        if not np.isfinite(float(self.reference_energy)):
            return False

        return True

    def matches_geometry(self, geometry):
        """
        :param geometry:
            Array of shape (n_atoms, 3) in bohr.

        :return:
            True when the geometry is bit-identical after canonicalization.
        """

        return geometry_fingerprint(geometry) == self.geometry_fingerprint

    def to_restart_dict(self):
        """
        Serializes the portable part of the snapshot for a checkpoint.

        Native backend handles (an OpenQP molecule, a Serenity LRSCF
        controller) are deliberately excluded: they are not serializable.  A
        restart recomputes them from ``geometry`` and ``settings_fingerprint``
        before the first tracked step, so tracking never restarts from an
        empty history.

        :return:
            A JSON-compatible dictionary.
        """

        return {
            'contract_revision': str(self.contract_revision),
            'backend': str(self.backend),
            'method': str(self.method),
            'backend_version': str(self.backend_version),
            'adapter_revision': str(self.adapter_revision),
            'calculation_id': str(self.calculation_id),
            'geometry': np.asarray(self.geometry, dtype=float).tolist(),
            'atom_labels': list(self.atom_labels),
            'geometry_fingerprint': str(self.geometry_fingerprint),
            'settings_fingerprint': str(self.settings_fingerprint),
            'charge': float(self.charge),
            'reference_multiplicity': int(self.reference_multiplicity),
            'target_manifold': str(self.target_manifold),
            'functional': str(self.functional),
            'basis': str(self.basis),
            'target_energies':
                np.asarray(self.target_energies, dtype=float).tolist(),
            'derivative_selectors': [int(s) for s in self.derivative_selectors],
            'reference_energy': float(self.reference_energy),
            'response_energies':
                np.asarray(self.response_energies, dtype=float).tolist(),
            'energy_unit': str(self.energy_unit),
            'coordinate_unit': str(self.coordinate_unit),
            'gradient_unit': str(self.gradient_unit),
            'reference_stability': (None if self.reference_stability is None
                                    else dict(self.reference_stability)),
            'reference_continuity': (
                None if self.reference_continuity is None else
                dict(self.reference_continuity)),
            'state_overlap_method': self.state_overlap_method,
        }

    @classmethod
    def from_restart_dict(cls, payload):
        """
        Rebuilds the portable part of a snapshot written by
        :func:`to_restart_dict`.

        The result carries no native handles and no overlap, so it can only be
        used to *locate* the accepted geometry that must be recomputed.
        """

        assert_msg_critical(
            str(payload.get('contract_revision')) ==
            SNAPSHOT_CONTRACT_REVISION,
            'ElectronicSnapshot: the checkpoint was written with an '
            'incompatible snapshot contract revision.')

        return cls(
            backend=str(payload['backend']),
            method=str(payload['method']),
            backend_version=str(payload['backend_version']),
            adapter_revision=str(payload['adapter_revision']),
            calculation_id=str(payload['calculation_id']),
            geometry=_frozen(payload['geometry']),
            atom_labels=tuple(payload['atom_labels']),
            geometry_fingerprint=str(payload['geometry_fingerprint']),
            settings_fingerprint=str(payload['settings_fingerprint']),
            charge=float(payload['charge']),
            reference_multiplicity=int(payload['reference_multiplicity']),
            target_manifold=str(payload['target_manifold']),
            functional=str(payload['functional']),
            basis=str(payload['basis']),
            target_energies=_frozen(payload['target_energies']),
            derivative_selectors=tuple(
                int(s) for s in payload['derivative_selectors']),
            reference_energy=float(payload['reference_energy']),
            response_energies=_frozen(payload['response_energies']),
            state_overlap_method=payload.get('state_overlap_method'),
            reference_stability=(
                None if payload.get('reference_stability') is None else
                dict(payload['reference_stability'])),
            reference_continuity=(
                None if payload.get('reference_continuity') is None else
                dict(payload['reference_continuity'])),
            energy_unit=str(payload['energy_unit']),
            coordinate_unit=str(payload['coordinate_unit']),
            gradient_unit=str(payload['gradient_unit']))


def build_snapshot(*,
                   backend,
                   method,
                   backend_version,
                   adapter_revision,
                   geometry,
                   atom_labels,
                   settings_digest,
                   charge,
                   reference_multiplicity,
                   target_manifold,
                   functional,
                   basis,
                   reference_energy,
                   target_energies,
                   derivative_selectors,
                   response_energies=None,
                   overlap_to_previous=None,
                   previous_calculation_id=None,
                   overlap_source=None,
                   state_overlap_method=None,
                   native_overlap_matrix=None,
                   response_overlap_matrix=None,
                   native_spectral_norm=None,
                   selected_spectral_norm=None,
                   overlap_is_conditioned=None,
                   spin_metadata=None,
                   scf_converged=True,
                   response_converged=True,
                   reference_stability=None,
                   reference_continuity=None,
                   provider_generation=0,
                   warnings=()):
    """
    Validates and freezes one electronic snapshot.

    Every production invariant that can be checked from the data alone is
    checked here, once, rather than at each consumer:

    * the working reference is not among the target energies by construction,
      because it is a separate field;
    * target energies are finite and complete;
    * derivative selectors are unique and cover every raw target;
    * an overlap, when present, is square, finite, and of the right size.

    :return:
        A frozen :class:`ElectronicSnapshot`.
    """

    coordinates = np.asarray(geometry, dtype=float)
    energies = np.asarray(target_energies, dtype=float).reshape(-1)
    n_states = int(energies.size)

    assert_msg_critical(
        n_states >= 1,
        'build_snapshot: at least one physical target state is required.')
    assert_msg_critical(
        np.all(np.isfinite(energies)),
        f'build_snapshot: the {backend} {method} calculation returned a '
        'nonfinite target-state energy.')
    assert_msg_critical(
        np.isfinite(float(reference_energy)),
        f'build_snapshot: the {backend} {method} calculation returned a '
        'nonfinite working-reference energy.')

    selectors = tuple(int(s) for s in derivative_selectors)
    assert_msg_critical(
        len(selectors) == n_states,
        'build_snapshot: one derivative selector is required per target '
        f'state; got {len(selectors)} for {n_states} states.')
    assert_msg_critical(
        len(set(selectors)) == n_states,
        'build_snapshot: derivative selectors must be distinct, otherwise '
        'two target states would alias onto the same backend gradient.')

    if response_energies is None:
        response_energies = energies - float(reference_energy)
    response = np.asarray(response_energies, dtype=float).reshape(-1)
    assert_msg_critical(
        response.size == n_states,
        'build_snapshot: response-energy count does not match the target '
        'states.')

    if overlap_to_previous is not None:
        overlap = np.asarray(overlap_to_previous, dtype=float)
        assert_msg_critical(
            overlap.ndim == 2 and overlap.shape[1] == n_states,
            'build_snapshot: the cross-geometry overlap must have one column '
            f'per current target state; got shape {overlap.shape} for '
            f'{n_states} states.')
        assert_msg_critical(
            np.all(np.isfinite(overlap)),
            'build_snapshot: the cross-geometry overlap contains nonfinite '
            'entries.')
        assert_msg_critical(
            previous_calculation_id is not None,
            'build_snapshot: an overlap must name the snapshot it was '
            'measured against.')

        for name, matrix in (
                ('native overlap', native_overlap_matrix),
                ('response-vector overlap', response_overlap_matrix)):
            if matrix is None:
                continue
            candidate = np.asarray(matrix, dtype=float)
            assert_msg_critical(
                candidate.shape == overlap.shape and
                np.all(np.isfinite(candidate)),
                f'build_snapshot: the {name} provenance must be finite and '
                'have the selected overlap shape.')

        for name, value in (
                ('native spectral norm', native_spectral_norm),
                ('selected spectral norm', selected_spectral_norm)):
            if value is not None:
                assert_msg_critical(
                    np.isfinite(float(value)) and float(value) >= 0.0,
                    f'build_snapshot: the {name} is invalid.')

    return ElectronicSnapshot(
        backend=str(backend),
        method=str(method),
        backend_version=str(backend_version),
        adapter_revision=str(adapter_revision),
        calculation_id=uuid.uuid4().hex,
        geometry=_frozen(coordinates),
        atom_labels=tuple(str(label) for label in atom_labels),
        geometry_fingerprint=geometry_fingerprint(coordinates),
        settings_fingerprint=str(settings_digest),
        charge=float(charge),
        reference_multiplicity=int(reference_multiplicity),
        target_manifold=str(target_manifold),
        functional=str(functional),
        basis=str(basis),
        target_energies=_frozen(energies),
        derivative_selectors=selectors,
        reference_energy=float(reference_energy),
        response_energies=_frozen(response),
        overlap_to_previous=_frozen(overlap_to_previous),
        previous_calculation_id=(None if previous_calculation_id is None else
                                 str(previous_calculation_id)),
        overlap_source=(None if overlap_source is None else
                        str(overlap_source)),
        state_overlap_method=(None if state_overlap_method is None else
                              str(state_overlap_method)),
        native_overlap_matrix=_frozen(native_overlap_matrix),
        response_overlap_matrix=_frozen(response_overlap_matrix),
        native_spectral_norm=(None if native_spectral_norm is None else
                              float(native_spectral_norm)),
        selected_spectral_norm=(None if selected_spectral_norm is None else
                                float(selected_spectral_norm)),
        overlap_is_conditioned=(None if overlap_is_conditioned is None else
                                bool(overlap_is_conditioned)),
        spin_metadata=(None if spin_metadata is None else dict(spin_metadata)),
        scf_converged=bool(scf_converged),
        response_converged=bool(response_converged),
        reference_stability=(None if reference_stability is None else
                             dict(reference_stability)),
        reference_continuity=(None if reference_continuity is None else
                              dict(reference_continuity)),
        provider_generation=int(provider_generation),
        warnings=tuple(warnings))


class ElectronicBackendAdapter:
    """
    Base class for backend adapters.

    An adapter owns exactly one responsibility: turn a nuclear geometry into a
    validated :class:`ElectronicSnapshot`, and turn a snapshot plus a raw
    target index into that target's gradient.  It owns no persistent state
    identity, no cache policy, and no hopping logic.

    :param molecule_template:
        A :class:`Molecule` supplying atom labels, charge and the working
        reference multiplicity.  Its coordinates are replaced at every
        geometry.
    :param number_of_states:
        Number of physical target states the dynamics tracks.
    """

    backend = 'abstract'
    method = 'abstract'
    adapter_revision = '1'
    target_manifold = 'unspecified'

    def __init__(self, molecule_template, number_of_states):

        assert_msg_critical(
            isinstance(molecule_template, Molecule),
            f'{type(self).__name__}: molecule_template must be a Molecule.')
        assert_msg_critical(
            int(number_of_states) >= 2,
            f'{type(self).__name__}: number_of_states must be >= 2 for '
            'surface hopping.')

        self.molecule_template = molecule_template
        self.number_of_states = int(number_of_states)
        self.atom_labels = tuple(molecule_template.get_labels())
        self.charge = float(molecule_template.get_charge())
        self.reference_multiplicity = int(
            molecule_template.get_multiplicity())

        self.n_scf_calls = 0
        self.n_response_calls = 0
        self.n_gradient_calls = 0
        self.n_overlap_calls = 0

        self._settings_digest = None

    # -- required interface ------------------------------------------------

    def compute_snapshot(self, geometry, previous_snapshot=None,
                         gradient_hint=None):
        """
        Runs the backend at one geometry.

        :param geometry:
            Array of shape (n_atoms, 3) in bohr.
        :param previous_snapshot:
            The previously *accepted* snapshot, or ``None`` at the first
            geometry.  When given, the returned snapshot must carry
            ``overlap_to_previous``.
        :param gradient_hint:
            Raw target index whose gradient the caller will most likely want.
            Backends that produce energies and one gradient in a single task
            use it to avoid a second task; others may ignore it.

        :return:
            A tuple ``(snapshot, gradient_or_None)``.
        """

        raise NotImplementedError

    def compute_gradient(self, snapshot, raw_target):
        """
        Evaluates the gradient of one raw target at the snapshot's geometry.

        :param snapshot:
            The :class:`ElectronicSnapshot` the gradient must belong to.
        :param raw_target:
            Raw target index in ``0 .. n_states - 1``.

        :return:
            Array of shape (n_atoms, 3) in Hartree/bohr.
        """

        raise NotImplementedError

    def validate_startup(self):
        """
        Checks that the installed backend can satisfy every production
        requirement, before a trajectory is started.

        :return:
            A dictionary of validated capabilities, for the run log.
        """

        raise NotImplementedError

    def release(self, calculation_id):
        """
        Drops the native payload of a snapshot that is no longer reachable.
        """

    # -- shared helpers ----------------------------------------------------

    def build_molecule(self, geometry):
        """
        Builds a :class:`Molecule` at a geometry from the template.

        Atom order, charge and the working-reference multiplicity always come
        from the template, so they cannot drift between geometries.
        """

        coordinates = np.asarray(geometry, dtype=float)

        assert_msg_critical(
            coordinates.shape == (len(self.atom_labels), 3),
            f'{type(self).__name__}: geometry shape {coordinates.shape} does '
            f'not match the {len(self.atom_labels)}-atom template.')

        molecule = Molecule(list(self.atom_labels), coordinates, units='au')
        molecule.set_charge(self.charge)
        molecule.set_multiplicity(self.reference_multiplicity)

        return molecule

    def settings_digest(self):
        """
        :return:
            The fingerprint of every result-affecting setting of this adapter.
        """

        if self._settings_digest is None:
            self._settings_digest = settings_fingerprint(
                self.describe_settings())

        return self._settings_digest

    def describe_settings(self):
        """
        :return:
            A JSON-serializable mapping of every setting that can change a
            result.  Subclasses extend it; the base contributes what every
            backend shares.
        """

        return {
            'contract_revision': SNAPSHOT_CONTRACT_REVISION,
            'backend': self.backend,
            'method': self.method,
            'adapter_revision': self.adapter_revision,
            'number_of_states': self.number_of_states,
            'atom_labels': list(self.atom_labels),
            'charge': self.charge,
            'reference_multiplicity': self.reference_multiplicity,
            'target_manifold': self.target_manifold,
            'coordinate_unit': 'bohr',
        }

    def _validate_snapshot_pair(self, snapshot, raw_target):
        """
        Refuses a gradient request that does not belong to this adapter's
        current configuration or to a valid raw target.
        """

        assert_msg_critical(
            isinstance(snapshot, ElectronicSnapshot),
            f'{type(self).__name__}: a gradient requires an '
            'ElectronicSnapshot.')
        assert_msg_critical(
            snapshot.backend == self.backend and
            snapshot.method == self.method,
            f'{type(self).__name__}: refusing to differentiate a '
            f'{snapshot.backend}/{snapshot.method} snapshot.')
        assert_msg_critical(
            snapshot.settings_fingerprint == self.settings_digest(),
            f'{type(self).__name__}: the snapshot was produced with '
            'different settings than the adapter is configured with now.')

        return snapshot.selector_for(raw_target)


class OpenQPMRSFAdapter(ElectronicBackendAdapter):
    """
    OpenQP MRSF-TDDFT adapter.

    One ``compute_snapshot`` performs a triplet ROHF reference, the MRSF
    response, and - once a previous accepted snapshot exists - OpenQP's own
    ``OQP::td_states_overlap`` plus an independent normalized response-vector
    overlap.  The optimizer's shared conditioning selector chooses the matrix
    used for assignment. Gradients of individual target states are then taken
    from the *same* OpenQP molecule, so a gradient can never belong to a
    different response calculation than its energy.

    :param molecule_template:
        Molecule supplying labels, charge and the triplet working reference.
    :param scf_driver:
        A configured :class:`OpenQPScfDriver` (basis, functional, scratch).
    :param number_of_states:
        Number of physical MRSF target states.
    :param exc_multiplicity:
        Multiplicity of the requested MRSF target manifold; 1 for singlets.
    :param native_cache_size:
        Number of OpenQP molecules kept alive for gradient evaluation.
    """

    backend = 'openqp'
    method = 'mrsf-tddft'
    adapter_revision = '3'
    target_manifold = 'singlet'

    #: OpenQP MRSF requires a high-spin triplet ROHF working reference.
    REQUIRED_REFERENCE_MULTIPLICITY = 3

    def __init__(self,
                 molecule_template,
                 scf_driver,
                 number_of_states,
                 exc_multiplicity=1,
                 native_cache_size=4):

        super().__init__(molecule_template, number_of_states)

        self.scf_driver = scf_driver
        self.exc_multiplicity = int(exc_multiplicity)
        self.native_cache_size = max(2, int(native_cache_size))
        self.target_manifold = ('singlet' if self.exc_multiplicity == 1 else
                                f'multiplicity-{self.exc_multiplicity}')

        assert_msg_critical(
            self.reference_multiplicity ==
            self.REQUIRED_REFERENCE_MULTIPLICITY,
            'OpenQPMRSFAdapter: MRSF-TDDFT requires a high-spin triplet '
            'working reference; set the template multiplicity to 3.  The '
            'triplet is a working reference only and is never a dynamics '
            'surface.')

        self._native = {}
        self._stability_policy_active = False
        self._reference_continuity_threshold = None

        # The native snapshot cache keeps whole OpenQP molecules alive, and
        # write_native_provenance() re-reads their scratch log and result
        # bundle.  The driver's scratch retention window must therefore cover
        # at least the cached snapshots plus the one being computed, or a
        # cached molecule would point at files the driver already retired.
        retention = getattr(self.scf_driver, 'set_scratch_retention', None)
        if callable(retention):
            retention(max(
                int(getattr(self.scf_driver, 'scratch_retention', 0)),
                self.native_cache_size + 1))

    def configure_reference_stability_policy(self, enabled):
        """Lets the trajectory policy selectively override driver stability."""

        self._stability_policy_active = bool(enabled)

    def configure_reference_continuity_policy(self, minimum_singular_value):
        """Sets an optional explicit occupied-subspace continuity cutoff."""

        self._reference_continuity_threshold = (
            None if minimum_singular_value is None else
            float(minimum_singular_value))

    def describe_settings(self):
        """
        See :func:`ElectronicBackendAdapter.describe_settings`.
        """

        driver = self.scf_driver
        payload = super().describe_settings()
        payload.update({
            'basis': str(driver.basis),
            'functional': str(getattr(driver, 'functional', '')),
            'scf_method': str(getattr(driver, 'method', '')),
            'scf_type': 'rohf',
            'max_cycles': int(getattr(driver, 'max_cycles', 0) or 0),
            'd4': bool(getattr(driver, 'd4', False)),
            'exc_multiplicity': self.exc_multiplicity,
            'tddft_type': 'mrsf',
            'tracking_overlap_algorithm': str(
                driver.tracking_overlap_algorithm),
        })

        return payload

    @property
    def basis_label(self):
        """
        :return:
            The configured basis-set label.
        """

        return str(self.scf_driver.basis)

    @property
    def functional_label(self):
        """
        :return:
            The configured functional label, empty for Hartree-Fock.
        """

        driver = self.scf_driver

        return str(driver.functional) if driver.method == 'dft' else ''

    def backend_version(self):
        """
        :return:
            The installed OpenQP version, or ``'unknown'``.
        """

        try:
            from importlib.metadata import version

            return str(version('openqp'))
        except Exception:
            return 'unknown'

    def validate_startup(self):
        """
        See :func:`ElectronicBackendAdapter.validate_startup`.
        """

        from .openqpscfdriver import OpenQPScfDriver

        if not OpenQPScfDriver.is_available():
            raise BackendCapabilityError(
                'OpenQPMRSFAdapter: the oqp module is not importable; OpenQP '
                'MRSF surface hopping is unavailable.')

        if not hasattr(self.scf_driver, 'run_oqp_mrsf_states'):
            raise BackendCapabilityError(
                'OpenQPMRSFAdapter: the installed OpenQPScfDriver does not '
                'provide run_oqp_mrsf_states(); multi-state MRSF energies '
                'cannot be obtained.')

        return {
            'backend': self.backend,
            'backend_version': self.backend_version(),
            'method': self.method,
            'basis': self.basis_label,
            'functional': self.functional_label,
            'reference_multiplicity': self.reference_multiplicity,
            'target_manifold': self.target_manifold,
            'number_of_states': self.number_of_states,
            'multi_state_total_energies': True,
            'state_specific_gradients': True,
            'cross_geometry_descriptor':
                'conditioned_openqp_mrsf_'
                f'{self.scf_driver.tracking_overlap_algorithm}_overlap_'
                'with_response_fallback',
            'scf_reference_stability': True,
            'scf_reference_continuity':
                'docc_somo_subspace_singular_values',
        }

    def compute_snapshot(self, geometry, previous_snapshot=None,
                         gradient_hint=None, stability_request=None,
                         reference_guess_payload=None):
        """
        See :func:`ElectronicBackendAdapter.compute_snapshot`.
        """

        molecule = self.build_molecule(geometry)
        previous_payload = None

        if previous_snapshot is not None:
            previous_payload = self._native_payload(previous_snapshot)

        driver_request = stability_request
        if self._stability_policy_active and driver_request is None:
            driver_request = {'stability_requested': False}

        driver_kwargs = {
            'functional': self.functional_label,
            'nstate': self.number_of_states,
            'previous_payload': previous_payload,
            'exc_multiplicity': self.exc_multiplicity,
        }
        if driver_request is not None:
            driver_kwargs['stability_request'] = driver_request
        if reference_guess_payload is not None:
            driver_kwargs['reference_guess_payload'] = reference_guess_payload
        if self._reference_continuity_threshold is not None:
            driver_kwargs['reference_continuity_threshold'] = (
                self._reference_continuity_threshold)

        try:
            native = self.scf_driver.run_oqp_mrsf_states(
                molecule, **driver_kwargs)
        finally:
            self.n_scf_calls += 1

        self.n_response_calls += 1

        energies = np.asarray(native['energies'], dtype=float).reshape(-1)

        # The measured OpenQP 1.2.0 contract.  A different length means the
        # installed OpenQP does not follow it, and guessing which slots are
        # target states would silently corrupt every potential.
        if energies.size != self.number_of_states + 1:
            raise BackendCapabilityError(
                'OpenQPMRSFAdapter: OpenQP returned '
                f'{energies.size} energies for nstate='
                f'{self.number_of_states}; the supported contract is exactly '
                f'{self.number_of_states + 1} entries, one working reference '
                'followed by one total energy per MRSF target state.')

        overlap = native.get('overlap', None)
        overlap_selection = None
        if previous_snapshot is not None:
            if overlap is None:
                raise BackendCapabilityError(
                    'OpenQPMRSFAdapter: OpenQP did not return a native state '
                    'overlap against the previous accepted geometry; '
                    'production state tracking cannot proceed.')

            # Apply the exact selector used by the optimization tracker even
            # when a test/downgraded driver only returns the candidate matrices.
            # This also verifies the provenance returned by the production
            # driver before a gradient can be evaluated.
            from .openqpstatetracker import select_mrsf_assignment_overlap

            native_overlap = native.get('native_overlap', None)
            if native_overlap is None:
                native_overlap = overlap
            overlap_selection = select_mrsf_assignment_overlap(
                native_overlap,
                response_overlap=native.get('response_overlap', None),
                nstates=self.number_of_states)
            overlap = overlap_selection.selected_overlap
            if not overlap_selection.is_conditioned:
                details = '; '.join(overlap_selection.warnings)
                raise BackendCapabilityError(
                    'OpenQPMRSFAdapter: no conditioned cross-geometry MRSF '
                    f'descriptor is available; refusing state assignment: '
                    f'{details}')
            self.n_overlap_calls += 1

        snapshot = build_snapshot(
            backend=self.backend,
            method=self.method,
            backend_version=self.backend_version(),
            adapter_revision=self.adapter_revision,
            geometry=geometry,
            atom_labels=self.atom_labels,
            settings_digest=self.settings_digest(),
            charge=self.charge,
            reference_multiplicity=self.reference_multiplicity,
            target_manifold=self.target_manifold,
            functional=self.functional_label,
            basis=self.basis_label,
            reference_energy=float(energies[0]),
            target_energies=energies[1:self.number_of_states + 1],
            derivative_selectors=tuple(
                range(1, self.number_of_states + 1)),
            response_energies=native.get('response_energies', None),
            overlap_to_previous=overlap,
            previous_calculation_id=(
                None if previous_snapshot is None else
                previous_snapshot.calculation_id),
            overlap_source=(None if overlap_selection is None else
                            overlap_selection.overlap_source),
            state_overlap_method=native.get('state_overlap_method', None),
            native_overlap_matrix=(None if overlap_selection is None else
                                   overlap_selection.native_overlap),
            response_overlap_matrix=(None if overlap_selection is None else
                                     overlap_selection.response_overlap),
            native_spectral_norm=(None if overlap_selection is None else
                                  overlap_selection.native_spectral_norm),
            selected_spectral_norm=(None if overlap_selection is None else
                                    overlap_selection.selected_spectral_norm),
            overlap_is_conditioned=(None if overlap_selection is None else
                                    overlap_selection.is_conditioned),
            scf_converged=bool(native.get('scf_converged', True)),
            response_converged=bool(native.get('response_converged', True)),
            reference_stability=native.get('reference_stability', None),
            reference_continuity=native.get('reference_continuity', None),
            warnings=(overlap_selection.warnings if
                      overlap_selection is not None else ()))

        self._store_native(snapshot, native)

        gradient = None
        if gradient_hint is not None:
            gradient = self.compute_gradient(snapshot, gradient_hint)

        return snapshot, gradient

    def compute_gradient(self, snapshot, raw_target):
        """
        See :func:`ElectronicBackendAdapter.compute_gradient`.
        """

        selector = self._validate_snapshot_pair(snapshot, raw_target)
        payload = self._native_payload(snapshot, for_gradient=True)

        gradient = self.scf_driver.compute_oqp_state_gradient(
            payload['molecule'], selector, len(self.atom_labels))
        self.n_gradient_calls += 1

        return np.asarray(gradient, dtype=float)

    def state_overlap(self, previous_snapshot, current_snapshot):
        """
        Returns the stored conditioned overlap of an ordered snapshot pair.

        OpenQP builds the state overlap *inside* the response calculation, by
        injecting the previous MO/response payload between the reference SCF
        and the current MRSF step.  It therefore cannot be recomputed from two
        finished snapshots, and the ordered pair is fixed at computation time.

        :return:
            Array with previous raw targets in rows, current in columns.
        """

        assert_msg_critical(
            current_snapshot.previous_calculation_id ==
            previous_snapshot.calculation_id,
            'OpenQPMRSFAdapter: the requested overlap pair does not match the '
            'pair the conditioned overlap was computed for; reversing or '
            're-pairing snapshots is a different quantity.')

        return np.asarray(current_snapshot.overlap_to_previous, dtype=float)

    def release(self, calculation_id):
        """
        See :func:`ElectronicBackendAdapter.release`.
        """

        self._native.pop(str(calculation_id), None)

    def reference_restart_payload(self, snapshot):
        """Returns the accepted SCF orbitals/densities for a checkpoint."""

        from .openqpscfdriver import scf_guess_payload

        payload = self._native_payload(snapshot)
        guess = scf_guess_payload(payload['data'])
        if 'OQP::VEC_MO_A' not in guess:
            return None
        return guess

    def write_native_provenance(self, snapshot, directory, step):
        """Archives OpenQP's own result bundle for one accepted snapshot.

        Returns compact SCF fields used by the provider for
        ``scf_reference.json``.  Test doubles and older drivers without the
        native serializer return ``None`` so generic surface-hopping tests do
        not create filesystem output.
        """

        saver = getattr(self.scf_driver, '_save_openqp_results', None)
        if not callable(saver):
            return None

        payload = self._native_payload(snapshot)
        molecule = payload.get('molecule', None)
        if molecule is None:
            return None

        saved = saver(molecule, force=True)
        if saved is None:
            raise OSError(
                'OpenQPMRSFAdapter: OpenQP did not create its requested '
                'provenance bundle.')

        directory = os.path.abspath(str(directory))
        result_name = f'openqp_step_{int(step):08d}.json'
        copies = {
            'openqp_results_file': (
                saved['openqp_results_file'], result_name),
            'openqp_settings_file': (
                saved['openqp_settings_file'], 'openqp_settings.json'),
            'openqp_log_file': (saved['openqp_log_file'], 'openqp.log'),
        }
        archived_files = {}
        for label, (source, name) in copies.items():
            destination = os.path.join(directory, name)
            shutil.copy2(source, destination)
            archived_files[label] = name

        data = payload.get('data', {})
        n_alpha = int(data.get('nelec_A', 0))
        n_beta = int(data.get('nelec_B', 0))

        def values(key):
            try:
                return np.asarray(data[key], dtype=float).reshape(-1).tolist()
            except (KeyError, TypeError, ValueError):
                return None

        n_orbitals = 0
        alpha_energies = values('OQP::E_MO_A')
        beta_energies = values('OQP::E_MO_B')
        if alpha_energies is not None:
            n_orbitals = len(alpha_energies)
        elif beta_energies is not None:
            n_orbitals = len(beta_energies)

        continuity = snapshot.reference_continuity or {}
        return {
            'archived_openqp_files': archived_files,
            'reference_identifier': continuity.get(
                'current_reference_identifier'),
            'reference_energy': float(snapshot.reference_energy),
            'scf_converged': bool(snapshot.scf_converged),
            'scf_iteration_count': None,
            'n_alpha_electrons': n_alpha,
            'n_beta_electrons': n_beta,
            'alpha_occupations': (
                [1.0] * n_alpha + [0.0] * max(0, n_orbitals - n_alpha)),
            'beta_occupations': (
                [1.0] * n_beta + [0.0] * max(0, n_orbitals - n_beta)),
            'alpha_orbital_energies': alpha_energies,
            'beta_orbital_energies': beta_energies,
            'reference_multiplicity': int(snapshot.reference_multiplicity),
            'stability_result': (
                None if snapshot.reference_stability is None else
                dict(snapshot.reference_stability)),
        }

    def _store_native(self, snapshot, native):
        """
        Retains the OpenQP molecule and the aligned response payload of a
        snapshot, bounded so that a long trajectory cannot accumulate one
        OpenQP molecule per nuclear step.
        """

        self._native[snapshot.calculation_id] = native

        while len(self._native) > self.native_cache_size:
            self._native.pop(next(iter(self._native)))

    def _native_payload(self, snapshot, for_gradient=False):
        """
        Looks up the native payload of a snapshot.

        A miss is an error rather than a silent recomputation: returning data
        from another calculation is precisely the stale-result failure the
        provider contract forbids.
        """

        payload = self._native.get(snapshot.calculation_id, None)

        if payload is None:
            what = ('a state-specific gradient' if for_gradient else
                    'a native state overlap')
            raise BackendCapabilityError(
                f'OpenQPMRSFAdapter: the native OpenQP payload of calculation '
                f'{snapshot.calculation_id} is no longer retained, so {what} '
                'cannot be produced without reusing another calculation. '
                'Increase native_cache_size or recompute the snapshot.')

        return payload


class SerenitySFAdapter(ElectronicBackendAdapter):
    """
    Serenity SF-TDA / SF-TDDFT adapter.

    The high-spin determinant is the *working reference* that generates the
    spin-flip targets.  It is never a dynamics surface, and the physical S0
    generally lies **below** it, so a negative SF response energy is
    physically correct and is accepted.

    Target-state total energies are ``E_reference + omega[r]``.  The Serenity
    gradient task was verified by central finite differences to return the
    total target gradient, so no reference gradient is added.

    :param molecule_template:
        Molecule supplying labels, charge and the high-spin reference
        multiplicity.
    :param gradient_driver:
        A configured :class:`SerenityExcitedStateGradientDriver` whose
        response driver has ``spinflip`` enabled.
    :param number_of_states:
        Number of physical SF target states.
    :param degeneracy_threshold:
        Threshold handed to Serenity's transition-density tracking.
    """

    backend = 'serenity'
    method = 'sf-tda'
    adapter_revision = '1'
    target_manifold = 'spin-flip'

    #: Serenity SF requires a high-spin working reference.
    MINIMUM_REFERENCE_MULTIPLICITY = 3

    def __init__(self,
                 molecule_template,
                 gradient_driver,
                 number_of_states,
                 degeneracy_threshold=0.0,
                 native_cache_size=4,
                 gradient_identity_threshold=0.9,
                 gradient_ambiguity_ratio=0.8,
                 gradient_energy_tolerance=1.0e-6):

        super().__init__(molecule_template, number_of_states)

        self.gradient_driver = gradient_driver
        self.response_driver = gradient_driver.rsp_driver
        self.scf_driver = gradient_driver.serenity_driver
        self.degeneracy_threshold = float(degeneracy_threshold)
        self.native_cache_size = max(2, int(native_cache_size))
        self.gradient_identity_threshold = float(
            gradient_identity_threshold)
        self.gradient_ambiguity_ratio = float(gradient_ambiguity_ratio)
        self.gradient_energy_tolerance = float(gradient_energy_tolerance)

        assert_msg_critical(
            0.0 <= self.gradient_identity_threshold <= 1.0,
            'SerenitySFAdapter: gradient_identity_threshold must be in '
            '[0, 1].')
        assert_msg_critical(
            0.0 <= self.gradient_ambiguity_ratio <= 1.0,
            'SerenitySFAdapter: gradient_ambiguity_ratio must be in [0, 1].')
        assert_msg_critical(
            self.gradient_energy_tolerance > 0.0,
            'SerenitySFAdapter: gradient_energy_tolerance must be positive.')

        assert_msg_critical(
            self.reference_multiplicity >=
            self.MINIMUM_REFERENCE_MULTIPLICITY,
            'SerenitySFAdapter: SF-TDDFT requires a high-spin working '
            'reference; set the template multiplicity to at least 3.  The '
            'high-spin determinant is a working reference only and is never '
            'a dynamics surface.')
        assert_msg_critical(
            bool(getattr(self.response_driver, 'spinflip', False)),
            'SerenitySFAdapter: the Serenity response driver must have '
            'spinflip enabled; without it the roots are ordinary excitations '
            'above the high-spin reference and the ladder is meaningless.')
        assert_msg_critical(
            int(self.response_driver.nstates) >= self.number_of_states,
            'SerenitySFAdapter: the Serenity response driver computes '
            f'{self.response_driver.nstates} roots but {self.number_of_states} '
            'SF target states are tracked.')

        self.method = ('sf-tda' if gradient_driver.exc_method == 'tda' else
                       'sf-tddft')

        # Backend spin selection must never silently renumber the ladder: the
        # adapter owns the raw-target -> selector map, and a driver-side
        # multiplicity search would break it.
        self.gradient_driver.enforce_same_multiplicity = False

        self._native = {}
        self._active_geometry_fingerprint = None

    def describe_settings(self):
        """
        See :func:`ElectronicBackendAdapter.describe_settings`.
        """

        scf = self.scf_driver
        rsp = self.response_driver
        payload = super().describe_settings()
        payload.update({
            'basis': str(scf.basis),
            'functional': str(scf.dft_functional),
            'scf_method': str(scf.method),
            'scf_mode': str(scf.scf_mode),
            'max_cycles': int(scf.max_cycles),
            'densfit_j': str(scf.densfit_j),
            'grid_accuracy': int(scf.grid_accuracy),
            'small_grid_accuracy': int(scf.small_grid_accuracy),
            'exc_method': str(self.gradient_driver.exc_method),
            'spinflip': True,
            'response_nstates': int(rsp.nstates),
            'response_conv_thresh': rsp.conv_thresh,
            'response_max_cycles': rsp.max_cycles,
            'gradient_identity_threshold': self.gradient_identity_threshold,
            'gradient_ambiguity_ratio': self.gradient_ambiguity_ratio,
            'gradient_energy_tolerance': self.gradient_energy_tolerance,
        })

        return payload

    @property
    def basis_label(self):
        """
        :return:
            The configured basis-set label.
        """

        return str(self.scf_driver.basis)

    @property
    def functional_label(self):
        """
        :return:
            The configured functional label.
        """

        return str(self.scf_driver.dft_functional)

    def backend_version(self):
        """
        :return:
            The installed Serenity version, or ``'unknown'``.
        """

        try:
            from importlib.metadata import version

            return str(version('qcserenity'))
        except Exception:
            return 'unknown'

    def validate_startup(self):
        """
        See :func:`ElectronicBackendAdapter.validate_startup`.
        """

        from .serenityscfdriver import SerenityScfDriver

        if not SerenityScfDriver.is_available():
            raise BackendCapabilityError(
                'SerenitySFAdapter: qcserenity.serenipy is not importable; '
                'Serenity SF surface hopping is unavailable.')

        try:
            from qcserenity import serenipy as spy
        except ImportError as exc:
            raise BackendCapabilityError(
                f'SerenitySFAdapter: cannot import serenipy: {exc}')

        for name in ('TransitionDensityTracking_R',
                     'TransitionDensityTracking_U'):
            if not hasattr(spy, name):
                raise BackendCapabilityError(
                    f'SerenitySFAdapter: serenipy lacks {name}; no validated '
                    'cross-geometry electronic descriptor is available and '
                    'production surface hopping is refused.')

        return {
            'backend': self.backend,
            'backend_version': self.backend_version(),
            'method': self.method,
            'basis': self.basis_label,
            'functional': self.functional_label,
            'reference_multiplicity': self.reference_multiplicity,
            'target_manifold': self.target_manifold,
            'number_of_states': self.number_of_states,
            'multi_state_total_energies': True,
            'state_specific_gradients': True,
            'cross_geometry_descriptor': 'serenity_transition_density_overlap',
        }

    def compute_snapshot(self, geometry, previous_snapshot=None,
                         gradient_hint=None):
        """
        See :func:`ElectronicBackendAdapter.compute_snapshot`.
        """

        from .veloxchemlib import hartree_in_ev

        molecule = self.build_molecule(geometry)
        raw_hint = 0 if gradient_hint is None else int(gradient_hint)
        selector = raw_hint + 1

        assert_msg_critical(
            0 <= raw_hint < self.number_of_states,
            f'SerenitySFAdapter: gradient hint {raw_hint} is outside the '
            f'{self.number_of_states} tracked SF target states.')

        task = self._run_gradient_task(molecule, selector)
        controller = task.getLRSCFController()

        omegas = np.asarray(controller.getExcitationEnergies('isolated'),
                            dtype=float).reshape(-1) / hartree_in_ev()

        if omegas.size < self.number_of_states:
            raise BackendCapabilityError(
                f'SerenitySFAdapter: Serenity converged {omegas.size} SF '
                f'roots but {self.number_of_states} target states are '
                'required.')

        reference_energy = self.scf_driver.get_energy()

        assert_msg_critical(
            reference_energy is not None,
            'SerenitySFAdapter: Serenity did not expose the high-spin '
            'working-reference energy; the SF target totals cannot be built.')
        reference_energy = float(reference_energy)

        omegas = omegas[:self.number_of_states]
        # A negative SF response energy is physical: the lowest SF target
        # normally lies BELOW the high-spin working reference.  It is kept.
        target_energies = reference_energy + omegas

        spin_metadata = self._spin_metadata(controller)
        overlap = None

        if previous_snapshot is not None:
            overlap = self._transition_density_overlap(previous_snapshot,
                                                       controller)
            self.n_overlap_calls += 1

        self.n_scf_calls += 1
        self.n_response_calls += 1
        self.n_gradient_calls += 1

        snapshot = build_snapshot(
            backend=self.backend,
            method=self.method,
            backend_version=self.backend_version(),
            adapter_revision=self.adapter_revision,
            geometry=geometry,
            atom_labels=self.atom_labels,
            settings_digest=self.settings_digest(),
            charge=self.charge,
            reference_multiplicity=self.reference_multiplicity,
            target_manifold=self.target_manifold,
            functional=self.functional_label,
            basis=self.basis_label,
            reference_energy=reference_energy,
            target_energies=target_energies,
            derivative_selectors=tuple(
                range(1, self.number_of_states + 1)),
            response_energies=omegas,
            overlap_to_previous=overlap,
            previous_calculation_id=(
                None if previous_snapshot is None else
                previous_snapshot.calculation_id),
            spin_metadata=spin_metadata)

        gradient = np.asarray(
            self.scf_driver._system.getGeometry().getGradients(), dtype=float)

        self._store_native(snapshot, {
            'controller': controller,
            'system': self.scf_driver._system,
            'mode': self.scf_driver._current_scf_mode,
            'geometry_object': self.scf_driver._system.getGeometry().copy(),
            'coefficients': controller.getCoefficients(),
            'excitation_vectors': controller.getExcitationVectors('isolated'),
            'molecule': molecule,
            'gradients': {raw_hint: gradient.copy()},
        })

        return snapshot, gradient.copy()

    def compute_gradient(self, snapshot, raw_target):
        """
        See :func:`ElectronicBackendAdapter.compute_gradient`.
        """

        selector = self._validate_snapshot_pair(snapshot, raw_target)
        payload = self._native_payload(snapshot)
        index = int(raw_target)

        cached = payload['gradients'].get(index, None)
        if cached is not None:
            return np.array(cached, dtype=float, copy=True)

        # Serenity's gradient task rebuilds its own LRSCF controller.  The
        # snapshot's descriptor payload must keep belonging to the response
        # calculation that produced the energies, so the new controller is
        # deliberately discarded rather than replacing it.
        task = self._run_gradient_task(payload['molecule'], selector)
        self.n_gradient_calls += 1
        mapped_selector = self._validate_recomputed_gradient_identity(
            snapshot, index, task.getLRSCFController())

        if mapped_selector != selector:
            # The response roots permuted between the energy snapshot and the
            # independent gradient task. The first gradient therefore belongs
            # to the wrong raw root and is discarded. Recompute once with the
            # descriptor-mapped selector, then prove that calculation retained
            # the same assignment before accepting its gradient.
            task = self._run_gradient_task(payload['molecule'],
                                           mapped_selector)
            self.n_gradient_calls += 1
            confirmed = self._validate_recomputed_gradient_identity(
                snapshot, index, task.getLRSCFController())
            if confirmed != mapped_selector:
                raise BackendCapabilityError(
                    'SerenitySFAdapter: the SF response roots changed again '
                    'while correcting a state-specific gradient; refusing an '
                    'unvalidated target force.')

        gradient = np.asarray(
            self.scf_driver._system.getGeometry().getGradients(), dtype=float)

        payload['gradients'][index] = gradient.copy()

        return gradient.copy()

    def state_overlap(self, previous_snapshot, current_snapshot):
        """
        Returns the transition-density overlap of an ordered snapshot pair.

        Unlike OpenQP, Serenity can measure this after the fact from the two
        stored response payloads, so the ordered pair can also be evaluated on
        demand.  Rows are previous raw targets, columns current raw targets.
        """

        current_payload = self._native_payload(current_snapshot)

        return self._transition_density_overlap(previous_snapshot,
                                                current_payload['controller'])

    def release(self, calculation_id):
        """
        See :func:`ElectronicBackendAdapter.release`.
        """

        self._native.pop(str(calculation_id), None)

    # -- internals ---------------------------------------------------------

    def _run_gradient_task(self, molecule, selector):
        """
        Runs one Serenity excited-state gradient task for a raw selector.

        The driver's own ``compute`` is deliberately bypassed so that a second
        state gradient at the *same* geometry does not pay for a fresh SCF.
        Its cache reset must therefore be reproduced here, and it must be
        reproduced on every geometry change: Serenity caches
        basis-function-on-grid controllers on the system object, and
        ``setCoordinates`` alone leaves them tied to the previous geometry.
        Reusing the system across a nuclear step was measured to shift the
        whole SF ladder by about 0.03 Hartree while leaving the response
        energies unchanged - a silent, uniform error in every potential.

        :param molecule:
            The molecule at the geometry to evaluate.
        :param selector:
            Serenity ``excGradList`` selector; raw target ``r`` uses ``r + 1``.

        :return:
            The finished Serenity gradient task.
        """

        driver = self.gradient_driver
        driver.set_state_deriv_index(int(selector))

        fingerprint = geometry_fingerprint(
            molecule.get_coordinates_in_bohr())

        if self._active_geometry_fingerprint != fingerprint:
            self.scf_driver._invalidate_cache()
            self.response_driver._invalidate_rsp_cache()
            self._active_geometry_fingerprint = fingerprint

        self.scf_driver._compute_energy_master(molecule)
        mode = self.scf_driver._current_scf_mode

        return driver._run_excited_gradient_task(mode)

    def _spin_metadata(self, controller):
        """
        Collects Serenity's spin diagnostics for every raw target.

        Serenity's SF roots interleave spin manifolds, so ``<S^2>`` and the
        inferred multiplicity are recorded per raw target and travel with the
        snapshot.  They are diagnostics: the raw target numbering is *not*
        renumbered by spin, because the derivative selector must stay
        ``r + 1``.
        """

        try:
            metadata = self.response_driver.get_spinflip_metadata(controller)
        except Exception:
            return None

        n = self.number_of_states

        def take(name):
            values = metadata.get(name, None)
            if values is None:
                return None
            return np.asarray(values).reshape(-1)[:n].tolist()

        return {
            'reference_s2': float(metadata.get('reference_s2', np.nan)),
            'state_s2': take('state_s2'),
            'state_multiplicities': take('state_multiplicities'),
            's2_deviation': take('s2_deviation'),
            's2_tolerance': float(self.gradient_driver.s2_tolerance),
        }

    def _validate_recomputed_gradient_identity(
            self, snapshot, raw_target, controller):
        """Matches an independent Serenity gradient task to its snapshot.

        A Serenity gradient task owns a fresh LRSCF controller. Root ordering
        is therefore not assumed to reproduce the energy snapshot. Native
        same-geometry transition-density overlaps provide the primary state
        identity; mapped total target energies and spin sectors are independent
        consistency checks.

        :return:
            The one-based selector in the recomputed controller that represents
            ``raw_target`` of ``snapshot``.
        """

        from .surfacehoppingtracker import (assignment_ambiguity,
                                            solve_assignment)
        from .veloxchemlib import hartree_in_ev

        overlap = self._transition_density_overlap(snapshot, controller)
        similarity = assignment_similarity(overlap)
        assignment = solve_assignment(similarity)
        scores = similarity[np.arange(self.number_of_states), assignment]
        ambiguity = assignment_ambiguity(similarity, assignment)

        if (float(np.min(scores)) < self.gradient_identity_threshold or
                ambiguity > self.gradient_ambiguity_ratio):
            raise BackendCapabilityError(
                'SerenitySFAdapter: the independently recomputed gradient '
                'roots cannot be assigned reliably to the energy snapshot '
                f'(minimum overlap {float(np.min(scores)):.6f}, ambiguity '
                f'ratio {ambiguity:.6f}).')

        current_omegas = np.asarray(
            controller.getExcitationEnergies('isolated'),
            dtype=float).reshape(-1) / hartree_in_ev()
        if current_omegas.size < self.number_of_states:
            raise BackendCapabilityError(
                'SerenitySFAdapter: the gradient task returned too few SF '
                'roots to validate the selected state.')

        reference_energy = self.scf_driver.get_energy()
        if reference_energy is None or not np.isfinite(float(reference_energy)):
            raise BackendCapabilityError(
                'SerenitySFAdapter: the gradient task returned no finite '
                'working-reference energy for state validation.')

        current_totals = (float(reference_energy) +
                          current_omegas[:self.number_of_states])
        mapped_totals = current_totals[assignment]
        deviations = np.abs(
            mapped_totals - np.asarray(snapshot.target_energies, dtype=float))
        largest = float(np.max(deviations))
        if largest > self.gradient_energy_tolerance:
            raise BackendCapabilityError(
                'SerenitySFAdapter: the recomputed gradient roots disagree '
                'with the energy snapshot after descriptor assignment by '
                f'{largest:.3e} Hartree (limit '
                f'{self.gradient_energy_tolerance:.3e}).')

        previous_spin = snapshot.spin_metadata
        current_spin = self._spin_metadata(controller)
        if previous_spin is not None and current_spin is not None:
            previous_mult = np.asarray(
                previous_spin.get('state_multiplicities', []),
                dtype=int).reshape(-1)
            current_mult = np.asarray(
                current_spin.get('state_multiplicities', []),
                dtype=int).reshape(-1)
            if (previous_mult.size != self.number_of_states or
                    current_mult.size != self.number_of_states or
                    not np.array_equal(current_mult[assignment],
                                       previous_mult)):
                raise BackendCapabilityError(
                    'SerenitySFAdapter: the recomputed gradient root maps to '
                    'an incompatible spin sector.')

        return int(assignment[int(raw_target)]) + 1

    def _transition_density_overlap(self, previous_snapshot, controller):
        """
        Evaluates Serenity's native transition-density overlap.

        Serenity returns ``|<Psi_i(current)|Psi_j(reference)>|`` with current
        roots in rows and reference roots in columns, already phase-invariant
        by construction.  The transpose restores the module-wide
        ``(previous, current)`` orientation, and the block is truncated to the
        tracked target states.
        """

        from qcserenity import serenipy as spy

        payload = self._native_payload(previous_snapshot)
        mode = self.scf_driver._current_scf_mode

        if mode != payload['mode']:
            raise BackendCapabilityError(
                'SerenitySFAdapter: the Serenity SCF mode changed from '
                f'{payload["mode"]} to {mode} between geometries; the '
                'transition-density overlap is not defined across that '
                'change.')

        tracker_cls = (spy.TransitionDensityTracking_R if mode == 'restricted'
                       else spy.TransitionDensityTracking_U)

        try:
            tracking = tracker_cls(
                self.scf_driver._system,
                controller,
                payload['excitation_vectors'],
                payload['geometry_object'],
                payload['coefficients'],
                spy.LRSCF_TYPE.ISOLATED,
                1,
                float(self.degeneracy_threshold))
            tracking.track()
            overlap = np.asarray(tracking.getOvlpMatrix())
        except Exception as exc:
            raise BackendCapabilityError(
                'SerenitySFAdapter: Serenity transition-density tracking '
                f'failed between the accepted and current geometry: {exc}')

        if np.iscomplexobj(overlap):
            if np.any(np.abs(np.imag(overlap)) > 1.0e-12):
                raise BackendCapabilityError(
                    'SerenitySFAdapter: Serenity returned complex '
                    'transition-density overlaps.')
            overlap = np.real(overlap)

        overlap = np.asarray(overlap, dtype=float)

        if overlap.ndim != 2:
            raise BackendCapabilityError(
                'SerenitySFAdapter: Serenity returned a non-matrix '
                'transition-density overlap.')

        n = self.number_of_states

        if overlap.shape[0] < n or overlap.shape[1] < n:
            raise BackendCapabilityError(
                'SerenitySFAdapter: the transition-density overlap has shape '
                f'{overlap.shape} but {n} tracked target states are required.')

        # (current, reference) -> (previous, current), truncated to the
        # tracked block.
        return np.ascontiguousarray(overlap[:n, :n].T)

    def _store_native(self, snapshot, native):
        """
        Retains the Serenity response payload of a snapshot, bounded.
        """

        self._native[snapshot.calculation_id] = native

        while len(self._native) > self.native_cache_size:
            self._native.pop(next(iter(self._native)))

    def _native_payload(self, snapshot):
        """
        Looks up the native payload of a snapshot; a miss is an error.
        """

        payload = self._native.get(snapshot.calculation_id, None)

        if payload is None:
            raise BackendCapabilityError(
                'SerenitySFAdapter: the Serenity response payload of '
                f'calculation {snapshot.calculation_id} is no longer '
                'retained. Producing a gradient or descriptor from another '
                'calculation is refused; increase native_cache_size or '
                'recompute the snapshot.')

        return payload


#: Backend identifiers accepted by the production driver.
BACKEND_ADAPTERS = {
    'openqp_mrsf': OpenQPMRSFAdapter,
    'serenity_sf': SerenitySFAdapter,
}


def create_backend_adapter(backend,
                           molecule_template,
                           number_of_states,
                           *,
                           basis,
                           functional,
                           exc_method='tda',
                           response_states=None,
                           exc_multiplicity=1,
                           scf_max_cycles=None,
                           response_max_cycles=150,
                           scratch_dir=None,
                           verbose=False):
    """
    Builds a configured backend adapter by name.

    This is an optional convenience factory.  Input scripts that need the
    usual VeloxChem driver-level control should construct the drivers
    themselves and instantiate :class:`OpenQPMRSFAdapter` or
    :class:`SerenitySFAdapter` directly.  The adapters and
    ``BackendStateProvider`` never require this factory.

    This is the single place a production driver selects its electronic
    backend.  There is no default and no fallback: an unknown name is an
    error, and a name that is installed but unusable raises
    :class:`BackendCapabilityError` rather than degrading to native VeloxChem
    TDA.

    :param backend:
        ``'openqp_mrsf'`` or ``'serenity_sf'``.
    :param molecule_template:
        Molecule supplying labels, charge and the high-spin working-reference
        multiplicity, which must be 3.
    :param number_of_states:
        Number of physical target states the dynamics tracks.
    :param basis:
        Basis-set label, in the backend's own spelling.
    :param functional:
        Exchange-correlation functional label, in the backend's own spelling.
    :param exc_method:
        Serenity response method, ``'tda'`` or ``'tddft'``.  Ignored by
        OpenQP, whose MRSF response is fixed.
    :param response_states:
        Number of roots the response solver computes; defaults to
        ``number_of_states``.  A larger window can stabilize the solver
        without changing which states the dynamics tracks.
    :param exc_multiplicity:
        Multiplicity of the requested MRSF target manifold (OpenQP only).
    :param scf_max_cycles:
        Maximum SCF iterations, or ``None`` for the driver default.
    :param response_max_cycles:
        Maximum response iterations.
    :param scratch_dir:
        Scratch directory.  Each backend calculation still receives its own
        unique output identity inside it, so concurrent or restarted runs
        cannot parse another calculation's files.
    :param verbose:
        Echo the native backend log.

    :return:
        A configured :class:`ElectronicBackendAdapter`.
    """

    key = str(backend).strip().lower()

    if key not in BACKEND_ADAPTERS:
        raise BackendCapabilityError(
            f'create_backend_adapter: unknown electronic backend {backend!r}; '
            'valid: ' + ', '.join(sorted(BACKEND_ADAPTERS)) + '.')

    roots = int(response_states or number_of_states)

    assert_msg_critical(
        roots >= int(number_of_states),
        'create_backend_adapter: response_states must be at least '
        'number_of_states.')

    if key == 'openqp_mrsf':
        from .openqpscfdriver import OpenQPScfDriver

        scf_driver = OpenQPScfDriver()
        scf_driver.set_basis(basis)
        scf_driver.set_functional(functional)
        scf_driver.openqp_verbose = bool(verbose)
        scf_driver.set_print_openqp_log(bool(verbose))
        scf_driver.set_save_openqp_json(False)
        if scf_max_cycles is not None:
            scf_driver.set_max_cycles(int(scf_max_cycles))
        if scratch_dir is not None:
            scf_driver.set_scratch_dir(scratch_dir)
        if not verbose:
            scf_driver.ostream.mute()

        # OpenQP's MRSF response returns exactly ``nstate`` target states, so
        # a wider window would change the snapshot contract rather than only
        # stabilizing the solver.
        assert_msg_critical(
            roots == int(number_of_states),
            'create_backend_adapter: OpenQP MRSF computes exactly '
            'number_of_states target states; response_states cannot differ.')

        adapter = OpenQPMRSFAdapter(molecule_template, scf_driver,
                                    number_of_states,
                                    exc_multiplicity=exc_multiplicity)
    else:
        from .serenityscfdriver import SerenityScfDriver
        from .serenityexcitedstategradientdriver import (
            SerenityExcitedStateGradientDriver)
        from .serenitylrrspeigensolver import SerenityLinearResponseSolver

        scf_driver = SerenityScfDriver()
        scf_driver.set_basis(basis)
        scf_driver.set_functional(functional)
        scf_driver.serenity_verbose = bool(verbose)
        if scf_max_cycles is not None:
            scf_driver.max_cycles = int(scf_max_cycles)
        if scratch_dir is not None:
            scf_driver.scratch_dir = scratch_dir
        if not verbose:
            scf_driver.ostream.mute()

        response_driver = SerenityLinearResponseSolver(scf_driver)
        response_driver.update_settings({
            'method': str(exc_method),
            'nstates': roots,
            'spinflip': True,
            'max_cycles': int(response_max_cycles),
        })

        gradient_driver = SerenityExcitedStateGradientDriver(scf_driver,
                                                             response_driver)
        gradient_driver.set_state_deriv_index(1)
        if not verbose:
            gradient_driver.ostream.mute()

        adapter = SerenitySFAdapter(molecule_template, gradient_driver,
                                    number_of_states)

    # Fail before a trajectory starts rather than at its first hop.
    adapter.validate_startup()

    return adapter


def create_backend_provider(backend, molecule_template, number_of_states,
                            cache_size=8, **kwargs):
    """
    Builds a production :class:`BackendStateProvider` by backend name.

    This convenience path constructs drivers internally.  For explicit
    dependency injection, construct the native driver stack in the input
    script, pass it to the matching adapter, and call
    ``BackendStateProvider(adapter)``.  Both setup paths produce the same
    runtime provider and do not change the surface-hopping controller.

    See :func:`create_backend_adapter` for the keyword arguments.

    :return:
        The tuple ``(provider, capabilities)``, where ``capabilities`` is the
        validated startup report suitable for the run log.
    """

    from .surfacehoppingprovider import BackendStateProvider

    adapter = create_backend_adapter(backend, molecule_template,
                                     number_of_states, **kwargs)
    provider = BackendStateProvider(adapter, cache_size=cache_size)

    return provider, provider.validate_startup()


__all__ = [
    'SNAPSHOT_CONTRACT_REVISION',
    'BackendCapabilityError',
    'ElectronicSnapshot',
    'ElectronicBackendAdapter',
    'OpenQPMRSFAdapter',
    'SerenitySFAdapter',
    'BACKEND_ADAPTERS',
    'create_backend_adapter',
    'create_backend_provider',
    'assignment_similarity',
    'build_snapshot',
    'geometry_fingerprint',
    'settings_fingerprint',
]
