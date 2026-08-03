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
Electronic-state providers for surface-hopping dynamics.

A provider owns the electronic-structure drivers and answers a single
question: given a nuclear geometry and an active state, what are the tracked
adiabatic state energies, the gradient of the active state, and the
descriptors needed for state tracking?

State-index convention
----------------------

Three index conventions coexist and are converted in exactly one place, in
:func:`StateIndexConverter`:

``tracked dynamics-state index``
    Persistent physical identity used by the controller, gap detector and hop
    history.  0 is initially the ground state and 1 the initially selected
    first excited state, but an excited tracked state may map to a different
    raw root after reordering.

``current raw response-root index``
    The state ordering returned at one geometry: 0 is the ground state and
    raw index ``i >= 1`` is response root ``i``.  The tracker permutation maps
    tracked indices to these current raw indices.

``driver index`` (``state_deriv_index``)
    Used by the VeloxChem excited-state gradient drivers.  It is one-based and
    counts excited roots only, so current raw state ``i >= 1`` corresponds to
    derivative index ``i``.  The ground state has no derivative index and is
    handled by the ground-state gradient driver instead.

For linear-response TDA/TDDFT the state energies are

.. math:: E_0 = E_\\mathrm{SCF}, \\qquad E_i = E_\\mathrm{SCF} + \\omega_i
"""

from dataclasses import dataclass
import numpy as np

from .errorhandler import assert_msg_critical
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .surfacehoppingtracker import StateDescriptors


class StateIndexConverter:
    """
    Single point of conversion among tracked dynamics-state indices, current
    raw response-root indices and excited-state-gradient driver indices.

    No ``+1`` or ``-1`` state-index arithmetic may appear anywhere else in the
    surface-hopping layer.
    """

    @staticmethod
    def validate_permutation(permutation):
        """
        Validates a tracked-to-raw state permutation.

        :param permutation:
            One-dimensional array whose entry ``i`` is the current raw root
            representing persistent tracked state ``i``.  Raw index 0 is the
            ground state; raw indices 1..N are response roots 1..N.

        :return:
            The validated integer permutation.
        """

        permutation = np.asarray(permutation, dtype=int)

        assert_msg_critical(
            permutation.ndim == 1,
            'StateIndexConverter: state permutation must be one-dimensional.')
        assert_msg_critical(
            np.array_equal(np.sort(permutation), np.arange(permutation.size)),
            'StateIndexConverter: state permutation must contain every raw '
            'state index exactly once.')

        return permutation

    @staticmethod
    def tracked_to_raw(tracked_state, permutation):
        """
        Maps a persistent tracked state to its current raw response root.
        """

        permutation = StateIndexConverter.validate_permutation(permutation)
        state = int(tracked_state)

        assert_msg_critical(
            0 <= state < permutation.size,
            'StateIndexConverter: tracked state index is out of range.')

        return int(permutation[state])

    @staticmethod
    def raw_to_tracked(raw_state, permutation):
        """
        Maps a current raw response root to its persistent tracked state.
        """

        permutation = StateIndexConverter.validate_permutation(permutation)
        state = int(raw_state)

        assert_msg_critical(
            0 <= state < permutation.size,
            'StateIndexConverter: raw state index is out of range.')

        inverse = np.empty_like(permutation)
        inverse[permutation] = np.arange(permutation.size)

        return int(inverse[state])

    @staticmethod
    def is_ground_state(dynamics_state):
        """
        :param dynamics_state:
            Current raw state index (0 is the ground state).

        :return:
            True when the state is the electronic ground state.
        """

        return int(dynamics_state) == 0

    @staticmethod
    def to_driver_index(dynamics_state):
        """
        :param dynamics_state:
            Current raw state index; must be >= 1.

        :return:
            The one-based excited-state driver index.
        """

        state = int(dynamics_state)

        assert_msg_critical(
            state >= 1,
            'StateIndexConverter: the ground state has no excited-state '
            'driver index.')

        return state

    @staticmethod
    def to_dynamics_index(driver_index):
        """
        :param driver_index:
            One-based excited-state driver index.

        :return:
            The corresponding dynamics state index.
        """

        index = int(driver_index)

        assert_msg_critical(
            index >= 1, 'StateIndexConverter: driver state index must be >= 1.')

        return index

    @staticmethod
    def excitation_index(dynamics_state):
        """
        :param dynamics_state:
            Current raw state index; must be >= 1.

        :return:
            The zero-based index into the array of excitation energies.
        """

        return StateIndexConverter.to_driver_index(dynamics_state) - 1

    @staticmethod
    def raw_to_driver_index(raw_state):
        """
        Maps a current raw state to the VeloxChem derivative-state index.

        :return:
            ``None`` for the ground state, otherwise the one-based response
            root index used by ``state_deriv_index``.
        """

        state = int(raw_state)

        if StateIndexConverter.is_ground_state(state):
            return None

        return StateIndexConverter.to_driver_index(state)

    @staticmethod
    def tracked_to_driver_index(tracked_state, permutation):
        """
        Maps a persistent tracked state directly to the VeloxChem
        derivative-state index at the current geometry.
        """

        raw_state = StateIndexConverter.tracked_to_raw(tracked_state,
                                                       permutation)

        return StateIndexConverter.raw_to_driver_index(raw_state)


@dataclass
class ElectronicStateResult:
    """
    All electronic-structure information belonging to one nuclear geometry.

    :param geometry:
        Array of shape (n_atoms, 3) with the QM geometry in bohr.
    :param ground_energy:
        The SCF ground-state energy, in Hartree.
    :param excitation_energies:
        Array of shape (n_states,) with excitation energies in Hartree; the
        first entry is 0.0 for the ground state.
    :param state_energies:
        Array of shape (n_states,) with total state energies in Hartree, in
        raw root order.
    :param active_state:
        Current raw state index whose gradient was evaluated.  Retained for
        compatibility; :attr:`active_raw_state` carries the explicit name.
    :param active_gradient:
        Array of shape (n_atoms, 3) with the active-state gradient, in Hartree
        per bohr.
    :param state_descriptors:
        :class:`StateDescriptors` for state tracking, in raw root order.
    :param state_permutation:
        Array mapping tracked index to raw root index; the identity until the
        tracker has run.
    :param tracking_scores:
        Per-state tracking similarity scores.
    :param scf_converged:
        Whether the SCF calculation converged.
    :param response_converged:
        Whether the response calculation converged.
    :param gradient_converged:
        Whether the gradient evaluation succeeded.
    :param active_raw_state:
        Current raw state index whose gradient was evaluated.
    :param active_tracked_state:
        Persistent tracked identity represented by ``active_raw_state``.
    :param derivative_state_index:
        One-based VeloxChem excited-state derivative index, or ``None`` for
        the ground-state gradient driver.
    """

    geometry: np.ndarray
    ground_energy: float
    excitation_energies: np.ndarray
    state_energies: np.ndarray
    active_state: int
    active_gradient: np.ndarray
    state_descriptors: object = None
    state_permutation: np.ndarray = None
    tracking_scores: np.ndarray = None
    scf_converged: bool = True
    response_converged: bool = True
    gradient_converged: bool = True
    active_raw_state: int = None
    active_tracked_state: int = None
    derivative_state_index: int = None

    def __post_init__(self):
        """
        Initializes explicit index metadata for provider results.
        """

        if self.active_raw_state is None:
            self.active_raw_state = int(self.active_state)

        if (self.derivative_state_index is None and
                not StateIndexConverter.is_ground_state(self.active_raw_state)):
            self.derivative_state_index = (
                StateIndexConverter.raw_to_driver_index(self.active_raw_state))

    def tracked_state_energies(self):
        """
        :return:
            The state energies reordered into tracked identity order.
        """

        if self.state_permutation is None:
            return np.asarray(self.state_energies, dtype=float)

        return np.asarray(self.state_energies,
                          dtype=float)[np.asarray(self.state_permutation,
                                                  dtype=int)]

    def raw_to_tracked_permutation(self):
        """
        :return:
            Array mapping each current raw state index to its tracked identity.
        """

        permutation = (np.arange(np.asarray(self.state_energies).size)
                       if self.state_permutation is None else
                       StateIndexConverter.validate_permutation(
                           self.state_permutation))
        inverse = np.empty_like(permutation)
        inverse[permutation] = np.arange(permutation.size)

        return inverse

    def tracked_state_descriptors(self):
        """
        :return:
            The state descriptors reordered into persistent tracked order.
        """

        if self.state_descriptors is None:
            return None

        permutation = (np.arange(self.state_descriptors.n_states)
                       if self.state_permutation is None else
                       StateIndexConverter.validate_permutation(
                           self.state_permutation))
        descriptors = self.state_descriptors

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

    def is_valid(self):
        """
        :return:
            True when all convergence flags are set and no energy or gradient
            component is nonfinite.
        """

        if not (self.scf_converged and self.response_converged and
                self.gradient_converged):
            return False

        if not np.all(np.isfinite(np.asarray(self.state_energies,
                                             dtype=float))):
            return False

        if not np.all(np.isfinite(np.asarray(self.active_gradient,
                                             dtype=float))):
            return False

        return True


class ElectronicStateProvider:
    """
    Base class for electronic-state providers.

    Subclasses implement :func:`_evaluate` and optionally
    :func:`compute_state_gradient`.  The base class owns geometry-keyed
    caching so that repeated queries at the same geometry - which happen
    during checkpoint-and-replay - never trigger a second SCF or response
    calculation.

    :param number_of_states:
        Number of tracked adiabatic states including the ground state.
    :param geometry_tolerance:
        Two geometries closer than this, in bohr, are treated as identical for
        caching purposes.
    :param cache_size:
        Maximum number of retained results.  The controller only ever queries
        the current frame, the trial frame and the restored central frame, so a
        small bound suffices; an unbounded cache would retain one full set of
        SCF and response results per nuclear step for the whole trajectory.
    """

    #: Smallest cache that still covers the checkpoint-and-replay chronology:
    #: the central frame, the discarded trial endpoint and the replayed
    #: endpoint, plus one entry of slack.
    MINIMUM_CACHE_SIZE = 4

    def __init__(self,
                 number_of_states,
                 geometry_tolerance=1.0e-10,
                 cache_size=8):

        assert_msg_critical(
            int(number_of_states) >= 2,
            'ElectronicStateProvider: number_of_states must be >= 2.')

        assert_msg_critical(
            int(cache_size) >= self.MINIMUM_CACHE_SIZE,
            'ElectronicStateProvider: cache_size must be at least '
            f'{self.MINIMUM_CACHE_SIZE} to cover checkpoint-and-replay.')

        self.number_of_states = int(number_of_states)
        self.geometry_tolerance = float(geometry_tolerance)
        self.cache_size = int(cache_size)

        self._cache = {}
        self.n_scf_calls = 0
        self.n_response_calls = 0
        self.n_gradient_calls = 0

    @property
    def number_of_roots(self):
        """
        :return:
            Number of excited-state roots to request from the response driver.
        """

        return self.number_of_states - 1

    def clear_cache(self):
        """
        Drops all cached results.
        """

        self._cache.clear()
        electronic_cache = getattr(self, '_electronic_data_cache', None)
        if electronic_cache is not None:
            electronic_cache.clear()

    def compute(self, geometry, active_state):
        """
        Returns the electronic-state result at a geometry, reusing the cache
        when the same geometry and active state were evaluated before.

        :param geometry:
            Array of shape (n_atoms, 3) with the QM geometry in bohr.
        :param active_state:
            Current raw state index whose gradient is requested.  The
            controller converts from its persistent tracked state before
            calling the provider and corrects the request after the current
            state permutation is known.

        :return:
            An :class:`ElectronicStateResult`.
        """

        geometry = np.asarray(geometry, dtype=float)
        key = self._cache_key(geometry, active_state)

        if key in self._cache:
            return self._cache[key]

        result = self._evaluate(geometry, int(active_state))

        assert_msg_critical(
            isinstance(result, ElectronicStateResult),
            'ElectronicStateProvider: _evaluate must return an '
            'ElectronicStateResult.')

        self._cache[key] = result
        self._evict_oldest(self._cache)

        return result

    def _evict_oldest(self, cache):
        """
        Drops the oldest entries of a cache until it fits :attr:`cache_size`.

        Python dictionaries preserve insertion order, so the first keys are the
        least recently inserted.

        :param cache:
            The dictionary to bound in place.
        """

        while len(cache) > self.cache_size:
            cache.pop(next(iter(cache)))

    def compute_state_gradient(self, geometry, state):
        """
        Evaluates the gradient of an arbitrary state on demand.

        Only called when ``gap_gradient`` momentum rescaling is enabled and a
        hop has already been proposed, so that ordinary steps never pay for an
        extra gradient.

        :param geometry:
            Array of shape (n_atoms, 3) with the QM geometry in bohr.
        :param state:
            Current raw state index at this geometry.

        :return:
            Array of shape (n_atoms, 3) with the gradient in Hartree/bohr.
        """

        raise NotImplementedError(
            'ElectronicStateProvider: this provider does not support '
            'on-demand gradients for arbitrary states; use global_qm '
            'momentum rescaling instead.')

    def _evaluate(self, geometry, active_state):
        """
        Performs the actual electronic-structure calculation.

        :param geometry:
            Array of shape (n_atoms, 3) with the QM geometry in bohr.
        :param active_state:
            Dynamics index of the active state.

        :return:
            An :class:`ElectronicStateResult`.
        """

        raise NotImplementedError

    def _cache_key(self, geometry, active_state):
        """
        Builds a cache key from a rounded geometry and the active state.

        :param geometry:
            Array of shape (n_atoms, 3) in bohr.
        :param active_state:
            Dynamics state index.

        :return:
            A hashable key.
        """

        decimals = max(0, int(-np.log10(self.geometry_tolerance)))
        rounded = np.round(np.asarray(geometry, dtype=float), decimals)

        return (int(active_state), rounded.tobytes())


class VeloxChemStateProvider(ElectronicStateProvider):
    """
    Provider for the native VeloxChem TDA/TDDFT stack.

    The target workflow per nuclear step is exactly one SCF calculation, one
    multiroot response calculation, and one gradient of the active state.
    Converged orbitals from the previous geometry are reused as the SCF
    starting guess whenever the driver supports it.

    :param molecule_template:
        A :class:`Molecule` supplying the atom labels, charge and
        multiplicity; its coordinates are replaced at each step.
    :param basis_label:
        Basis-set label, e.g. ``'def2-svp'``.
    :param scf_driver:
        An SCF driver instance, e.g. :class:`ScfRestrictedDriver`.
    :param rsp_driver:
        A response driver instance, e.g. :class:`TdaEigenSolver` or
        :class:`LinearResponseEigenSolver`.
    :param scf_gradient_driver:
        Ground-state gradient driver, e.g. :class:`ScfGradientDriver`.
    :param excited_gradient_driver:
        Excited-state gradient driver, e.g. :class:`TddftGradientDriver`.
    :param number_of_states:
        Number of tracked states including the ground state.
    """

    def __init__(self, molecule_template, basis_label, scf_driver, rsp_driver,
                 scf_gradient_driver, excited_gradient_driver,
                 number_of_states):

        super().__init__(number_of_states)

        self.molecule_template = molecule_template
        self.basis_label = basis_label
        self.scf_driver = scf_driver
        self.rsp_driver = rsp_driver
        self.scf_gradient_driver = scf_gradient_driver
        self.excited_gradient_driver = excited_gradient_driver

        # Request enough roots to cover every tracked excited state.
        if hasattr(self.rsp_driver, 'nstates'):
            self.rsp_driver.nstates = max(int(self.rsp_driver.nstates or 0),
                                          self.number_of_roots)

        self._last_scf_results = None
        self._last_rsp_results = None
        self._last_basis = None
        self._last_molecule = None
        self._electronic_data_cache = {}

    @staticmethod
    def is_tamm_dancoff(rsp_results):
        """
        Decides whether a response result comes from a Tamm-Dancoff (TDA)
        eigensolver or from a full RPA/TDDFT eigensolver.

        This is inferred from the result dictionary itself rather than from
        the driver class: the TDA solver returns dense ``eigenvectors`` while
        the RPA solver returns ``eigenvectors_distributed``.  Inferring it
        here cannot drift out of sync with the solvers.

        :param rsp_results:
            The response driver result dictionary.

        :return:
            True for TDA, False for RPA.
        """

        return 'eigenvectors_distributed' not in rsp_results

    def _configure_gradient_driver(self, rsp_results):
        """
        Propagates the TDA/RPA choice to the excited-state gradient driver.

        The orbital-response step reads the excitation vectors from different
        keys of the response-result dictionary depending on this flag
        (``eigenvectors`` for TDA, ``eigenvectors_distributed`` for RPA), so a
        mismatch fails with an opaque ``KeyError``.  Keeping the two drivers
        consistent is done here once rather than being left to the user.

        :param rsp_results:
            The response driver result dictionary.
        """

        driver = self.excited_gradient_driver

        if driver is None or not hasattr(driver, 'update_settings'):
            return

        tda_flag = 'yes' if self.is_tamm_dancoff(rsp_results) else 'no'

        driver.update_settings({}, {'tamm_dancoff': tda_flag})

    def build_molecule(self, geometry):
        """
        Builds a :class:`Molecule` at the given geometry from the template.

        :param geometry:
            Array of shape (n_atoms, 3) in bohr.

        :return:
            A new :class:`Molecule` in bohr.
        """

        labels = self.molecule_template.get_labels()
        geometry = np.asarray(geometry, dtype=float)

        assert_msg_critical(
            geometry.shape == (len(labels), 3),
            'VeloxChemStateProvider: geometry shape does not match the '
            'molecule template.')

        molecule = Molecule(labels, geometry, units='au')
        molecule.set_charge(self.molecule_template.get_charge())
        molecule.set_multiplicity(self.molecule_template.get_multiplicity())

        return molecule

    def _evaluate(self, geometry, active_state):
        """
        See :func:`ElectronicStateProvider._evaluate`.
        """

        molecule = self.build_molecule(geometry)
        basis = MolecularBasis.read(molecule, self.basis_label)

        # --- SCF -----------------------------------------------------------
        scf_results = self.scf_driver.compute(molecule, basis)
        self.n_scf_calls += 1

        scf_converged = bool(getattr(self.scf_driver, 'is_converged', True))
        ground_energy = float(self.scf_driver.get_scf_energy())

        # --- response ------------------------------------------------------
        rsp_results = self.rsp_driver.compute(molecule, basis, scf_results)
        self.n_response_calls += 1

        response_converged = bool(getattr(self.rsp_driver, 'is_converged',
                                          True))
        response_converged = response_converged and bool(rsp_results)

        eigenvalues = np.asarray(rsp_results.get('eigenvalues', []),
                                 dtype=float)

        assert_msg_critical(
            eigenvalues.size >= self.number_of_roots,
            'VeloxChemStateProvider: the response driver returned '
            f'{eigenvalues.size} roots but {self.number_of_roots} are '
            'required.')

        # yapf: disable
        excitation_energies = np.concatenate(
            [[0.0], eigenvalues[:self.number_of_roots]])
        # yapf: enable
        state_energies = ground_energy + excitation_energies

        # --- descriptors ---------------------------------------------------
        descriptors = self._build_descriptors(rsp_results, excitation_energies)

        # SCF/response data are keyed by geometry independently of the
        # requested gradient state.  Checkpoint-and-replay may request a
        # target gradient at the restored central frame after a trial
        # endpoint has already been evaluated.
        geometry_key = self._geometry_key(geometry)
        self._electronic_data_cache[geometry_key] = (molecule, basis,
                                                     scf_results, rsp_results)
        self._evict_oldest(self._electronic_data_cache)

        # --- active-state gradient -----------------------------------------
        gradient = self._compute_gradient(molecule, basis, scf_results,
                                          rsp_results, active_state)

        self._last_scf_results = scf_results
        self._last_rsp_results = rsp_results
        self._last_basis = basis
        self._last_molecule = molecule

        return ElectronicStateResult(
            geometry=np.asarray(geometry, dtype=float),
            ground_energy=ground_energy,
            excitation_energies=excitation_energies,
            state_energies=state_energies,
            active_state=int(active_state),
            active_gradient=gradient,
            state_descriptors=descriptors,
            state_permutation=np.arange(self.number_of_states),
            tracking_scores=np.ones(self.number_of_states),
            scf_converged=scf_converged,
            response_converged=response_converged,
            gradient_converged=bool(np.all(np.isfinite(gradient))))

    def _build_descriptors(self, rsp_results, excitation_energies):
        """
        Assembles :class:`StateDescriptors` from a response-result dictionary.

        The ground state is prepended with zero oscillator strength and a null
        transition dipole.

        :param rsp_results:
            The response driver result dictionary.
        :param excitation_energies:
            Array of shape (n_states,) with excitation energies in Hartree.

        :return:
            A :class:`StateDescriptors`.
        """

        n_roots = self.number_of_roots

        oscillator = rsp_results.get('oscillator_strengths', None)
        if oscillator is not None:
            oscillator = np.concatenate([[0.0],
                                         np.asarray(oscillator,
                                                    dtype=float)[:n_roots]])

        dipoles = rsp_results.get('electric_transition_dipoles', None)
        if dipoles is not None:
            dipoles = np.vstack(
                [np.zeros((1, 3)),
                 np.asarray(dipoles, dtype=float)[:n_roots]])

        return StateDescriptors(excitation_energies=excitation_energies,
                                oscillator_strengths=oscillator,
                                transition_dipoles=dipoles)

    def _compute_gradient(self, molecule, basis, scf_results, rsp_results,
                          state):
        """
        Evaluates the gradient of one state.

        Dispatches to the ground-state gradient driver for current raw state 0
        and to the excited-state gradient driver otherwise.

        :param molecule:
            The :class:`Molecule` at the current geometry.
        :param basis:
            The :class:`MolecularBasis`.
        :param scf_results:
            The SCF result dictionary.
        :param rsp_results:
            The response result dictionary.
        :param state:
            Current raw state index.

        :return:
            Array of shape (n_atoms, 3) with the gradient in Hartree/bohr.
        """

        self.n_gradient_calls += 1

        if StateIndexConverter.is_ground_state(state):
            self.scf_gradient_driver.compute(molecule, basis, scf_results)
            return np.asarray(self.scf_gradient_driver.get_gradient(),
                              dtype=float)

        driver_index = StateIndexConverter.to_driver_index(state)

        # The gradient driver must agree with the response driver on TDA/RPA.
        self._configure_gradient_driver(rsp_results)

        assert_msg_critical(
            driver_index <= self.number_of_roots,
            'VeloxChemStateProvider: requested gradient for state '
            f'{state} but only {self.number_of_roots} roots are tracked.')

        self.excited_gradient_driver.state_deriv_index = (driver_index,)
        self.excited_gradient_driver.compute(molecule, basis, self.scf_driver,
                                             self.rsp_driver, rsp_results)

        gradient = np.asarray(self.excited_gradient_driver.get_gradient(),
                              dtype=float)

        # TddftGradientDriver returns one gradient per requested state, with
        # shape (n_states, n_atoms, 3); exactly one state is requested here.
        if gradient.ndim == 3:
            assert_msg_critical(
                gradient.shape[0] == 1,
                'VeloxChemStateProvider: expected a single excited-state '
                f'gradient but received {gradient.shape[0]}.')
            gradient = gradient[0]

        return gradient

    def compute_state_gradient(self, geometry, state):
        """
        See :func:`ElectronicStateProvider.compute_state_gradient`.

        Reuses the cached SCF and response results at the same geometry, so
        only the gradient itself is recomputed.
        """

        geometry = np.asarray(geometry, dtype=float)

        geometry_key = self._geometry_key(geometry)

        assert_msg_critical(
            geometry_key in self._electronic_data_cache,
            'VeloxChemStateProvider: no cached SCF/response data are '
            'available for the requested gradient geometry.')

        molecule, basis, scf_results, rsp_results = (
            self._electronic_data_cache[geometry_key])

        return self._compute_gradient(molecule, basis, scf_results, rsp_results,
                                      state)

    def _geometry_key(self, geometry):
        """
        Builds the state-independent key used for cached SCF/response data.
        """

        decimals = max(0, int(-np.log10(self.geometry_tolerance)))
        rounded = np.round(np.asarray(geometry, dtype=float), decimals)

        return rounded.tobytes()


class ExternalDriverStateProvider(ElectronicStateProvider):
    """
    Provider for the external excited-state backends that follow the
    VeloxChem ``state_deriv_index`` / ``total_energy`` / ``get_gradient``
    driver contract, namely the Serenity and OpenQP excited-state gradient
    drivers.

    These drivers evaluate the reference, the excitations and the gradient of
    a single selected root in one ``compute`` call, so the provider sets the
    driver's state index, calls ``compute`` once, and reads the excitation
    spectrum back out of the driver.

    :param molecule_template:
        A :class:`Molecule` supplying atom labels, charge and multiplicity.
    :param excited_gradient_driver:
        A Serenity or OpenQP excited-state gradient driver.
    :param number_of_states:
        Number of tracked states including the ground state.
    :param ground_gradient_driver:
        Optional ground-state gradient driver, required only when the ground
        state is reachable by hopping.
    :param scf_driver:
        Optional SCF driver used for the ground-state energy.
    """

    def __init__(self,
                 molecule_template,
                 excited_gradient_driver,
                 number_of_states,
                 ground_gradient_driver=None,
                 scf_driver=None):

        super().__init__(number_of_states)

        self.molecule_template = molecule_template
        self.excited_gradient_driver = excited_gradient_driver
        self.ground_gradient_driver = ground_gradient_driver
        self.scf_driver = scf_driver

    def build_molecule(self, geometry):
        """
        Builds a :class:`Molecule` at the given geometry from the template.

        :param geometry:
            Array of shape (n_atoms, 3) in bohr.

        :return:
            A new :class:`Molecule`.
        """

        labels = self.molecule_template.get_labels()
        molecule = Molecule(labels,
                            np.asarray(geometry, dtype=float),
                            units='au')
        molecule.set_charge(self.molecule_template.get_charge())
        molecule.set_multiplicity(self.molecule_template.get_multiplicity())

        return molecule

    def _evaluate(self, geometry, active_state):
        """
        See :func:`ElectronicStateProvider._evaluate`.
        """

        molecule = self.build_molecule(geometry)
        driver = self.excited_gradient_driver

        assert_msg_critical(
            not StateIndexConverter.is_ground_state(active_state) or
            self.ground_gradient_driver is not None,
            'ExternalDriverStateProvider: a ground-state gradient driver is '
            'required to run on the ground state.')

        if StateIndexConverter.is_ground_state(active_state):
            self.ground_gradient_driver.compute(molecule)
            gradient = np.asarray(self.ground_gradient_driver.get_gradient(),
                                  dtype=float)
            self.n_gradient_calls += 1
            # The ground-state gradient driver leaves the excitation spectrum
            # cached on the excited-state driver at whatever geometry it last
            # saw.  Reading it here would freeze the gaps of an S0 trajectory
            # at a stale geometry, so refresh it explicitly.
            self._refresh_excitation_spectrum(molecule)
        else:
            driver.set_state_deriv_index(
                StateIndexConverter.to_driver_index(active_state))
            driver.compute(molecule)
            gradient = np.asarray(driver.get_gradient(), dtype=float)
            self.n_scf_calls += 1
            self.n_response_calls += 1
            self.n_gradient_calls += 1

        excitation_energies = self._read_excitation_energies(driver)
        ground_energy = self._read_ground_energy(driver, active_state,
                                                 excitation_energies)

        state_energies = ground_energy + excitation_energies

        descriptors = StateDescriptors(excitation_energies=excitation_energies,
                                       oscillator_strengths=self._read_optional(
                                           driver, 'oscillator_strengths'),
                                       transition_dipoles=self._read_optional(
                                           driver,
                                           'electric_transition_dipoles'))

        return ElectronicStateResult(
            geometry=np.asarray(geometry, dtype=float),
            ground_energy=ground_energy,
            excitation_energies=excitation_energies,
            state_energies=state_energies,
            active_state=int(active_state),
            active_gradient=gradient,
            state_descriptors=descriptors,
            state_permutation=np.arange(self.number_of_states),
            tracking_scores=np.ones(self.number_of_states),
            scf_converged=True,
            response_converged=True,
            gradient_converged=bool(np.all(np.isfinite(gradient))))

    def _refresh_excitation_spectrum(self, molecule):
        """
        Recomputes the excitation spectrum at the current geometry without
        asking for an excited-state gradient.

        Needed on ground-state steps, where the excited-state driver is not
        invoked and its cached spectrum belongs to an earlier geometry.

        :param molecule:
            The :class:`Molecule` at the current geometry.
        """

        driver = self.excited_gradient_driver
        compute_energy = getattr(driver, 'compute_energy', None)

        assert_msg_critical(
            callable(compute_energy),
            'ExternalDriverStateProvider: running on the ground state requires '
            'the excited-state driver to expose compute_energy(molecule) so '
            'that the excitation spectrum can be evaluated at the current '
            'geometry; otherwise the tracked gaps would be stale.')

        compute_energy(molecule)
        self.n_scf_calls += 1
        self.n_response_calls += 1

    def compute_state_gradient(self, geometry, state):
        """
        See :func:`ElectronicStateProvider.compute_state_gradient`.

        The external single-root drivers evaluate the reference, the
        excitations and one selected root's gradient in a single ``compute``
        call, so an on-demand gradient for another root repeats that whole
        calculation.  This is only requested after the tracker has moved the
        active state onto a different raw root, or after a hop has been
        accepted, so it costs at most one extra evaluation per hopping event.

        :param geometry:
            Array of shape (n_atoms, 3) with the QM geometry in bohr.
        :param state:
            Current raw response-root index at this geometry; 0 is the ground
            state.

        :return:
            Array of shape (n_atoms, 3) with the gradient in Hartree/bohr.
        """

        molecule = self.build_molecule(geometry)

        if StateIndexConverter.is_ground_state(state):
            assert_msg_critical(
                self.ground_gradient_driver is not None,
                'ExternalDriverStateProvider: a ground-state gradient driver '
                'is required to evaluate the ground-state gradient.')
            self.ground_gradient_driver.compute(molecule)
            self.n_gradient_calls += 1

            return np.asarray(self.ground_gradient_driver.get_gradient(),
                              dtype=float)

        driver = self.excited_gradient_driver
        driver_index = StateIndexConverter.to_driver_index(state)

        assert_msg_critical(
            driver_index <= self.number_of_roots,
            'ExternalDriverStateProvider: requested the gradient of raw state '
            f'{state} but only {self.number_of_roots} roots are tracked.')

        driver.set_state_deriv_index(driver_index)
        driver.compute(molecule)
        self.n_scf_calls += 1
        self.n_response_calls += 1
        self.n_gradient_calls += 1

        return np.asarray(driver.get_gradient(), dtype=float)

    def _read_excitation_energies(self, driver):
        """
        Extracts the excitation spectrum from an external driver.

        :param driver:
            The excited-state gradient driver.

        :return:
            Array of shape (n_states,) with the ground state prepended.
        """

        energies = getattr(driver, 'excited_state_energy', None)

        assert_msg_critical(
            energies is not None,
            'ExternalDriverStateProvider: the excited-state driver did not '
            'expose an excitation spectrum.')

        energies = np.asarray(energies, dtype=float).ravel()

        assert_msg_critical(
            energies.size >= self.number_of_roots,
            'ExternalDriverStateProvider: the driver returned '
            f'{energies.size} roots but {self.number_of_roots} are required.')

        return np.concatenate([[0.0], energies[:self.number_of_roots]])

    def _read_ground_energy(self, driver, active_state, excitation_energies):
        """
        Recovers the reference (ground-state) energy.

        Prefers an explicit SCF driver, and otherwise back-solves it from the
        driver's total energy of the active state.

        :param driver:
            The excited-state gradient driver.
        :param active_state:
            Dynamics index of the active state.
        :param excitation_energies:
            Array of excitation energies in Hartree.

        :return:
            The ground-state energy in Hartree.
        """

        if self.scf_driver is not None:
            energy = self.scf_driver.get_energy()
            if energy is not None:
                return float(energy)

        total_energy = getattr(driver, 'total_energy', None)

        assert_msg_critical(
            total_energy is not None,
            'ExternalDriverStateProvider: neither an SCF energy nor a total '
            'energy is available.')

        return float(total_energy) - float(
            excitation_energies[int(active_state)])

    def _read_optional(self, driver, name):
        """
        Reads an optional descriptor array from a driver, prepending the
        ground-state entry.

        :param driver:
            The excited-state gradient driver.
        :param name:
            Attribute name to look up.

        :return:
            The descriptor array, or None when unavailable.
        """

        values = getattr(driver, name, None)

        if values is None:
            return None

        values = np.asarray(values, dtype=float)

        if values.ndim == 1:
            if values.size < self.number_of_roots:
                return None
            return np.concatenate([[0.0], values[:self.number_of_roots]])

        if values.ndim == 2 and values.shape[1] == 3:
            if values.shape[0] < self.number_of_roots:
                return None
            return np.vstack([np.zeros((1, 3)), values[:self.number_of_roots]])

        return None
