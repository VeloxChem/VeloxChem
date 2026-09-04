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
Adiabatic Landau-Zener surface-hopping dynamics controller.

The controller composes - rather than inherits from - :class:`OpenMMDynamics`,
so that all existing OpenMM workflows remain untouched.  It owns the
chronology of a surface-hopping trajectory:

1. save the complete OpenMM state at frame ``n``;
2. propagate a trial step from ``n`` to ``n+1`` on the current surface;
3. evaluate the energy-ordered adiabatic states at ``n+1`` and diagnose their
   character continuity;
4. detect whether an adiabatic energy gap had a local minimum at frame ``n``;
5. evaluate the Landau-Zener probability and draw a random number;
6. keep the trial frame when no hop is selected;
7. otherwise restore frame ``n``, rescale momenta, switch the active state,
   install the new force, and replay the step.

Because a three-point minimum at frame ``n`` can only be recognized once frame
``n+1`` exists, the hop is applied at the central frame ``n`` and the trial
endpoint is discarded.  The fitted event time ``t_n + tau_c`` is recorded in
the diagnostics so that a later implementation can split the nuclear timestep
into pre-hop and post-hop segments without changing this controller's
interface.
"""

from dataclasses import dataclass
from pathlib import Path
import copy
import json
import pickle
import sys
import numpy as np

from .veloxchemlib import hartree_in_kjpermol, bohr_in_angstrom
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .landauzener import LandauZenerEventDetector
from .surfacehoppingprovider import (ElectronicStateProvider,
                                     StateIndexConverter)
from .surfacehoppingrescaler import SurfaceHopMomentumRescaler
from .surfacehoppingsettings import (SurfaceHoppingSettings, AU_IN_FS,
                                     NM_IN_BOHR, NM_PER_PS_IN_AU_VELOCITY,
                                     AU_VELOCITY_IN_NM_PER_PS,
                                     AMU_IN_ELECTRON_MASS)
from .surfacehoppingtracker import (create_state_tracker,
                                    is_production_tracker)

try:
    import openmm as mm
    import openmm.unit as unit
except ImportError:
    pass


class SurfaceHoppingError(RuntimeError):
    """
    Raised when a surface-hopping trajectory must be aborted.

    The controller writes a restartable checkpoint before raising whenever
    checkpointing is enabled.
    """


@dataclass
class TrajectoryFrame:
    """
    Complete restorable state of the trajectory at one frame.

    :param step:
        Integer frame index.
    :param time_fs:
        Simulation time in femtoseconds.
    :param positions_nm:
        Array of shape (n_particles, 3) with all OpenMM positions in nm.
    :param velocities_nm_ps:
        Array of shape (n_particles, 3) with all OpenMM velocities in nm/ps.
    :param box_vectors:
        Periodic box vectors in nm, or None for a non-periodic system.
    :param active_state:
        Dynamics index of the active state.
    :param electronic_result:
        The :class:`ElectronicStateResult` at this frame.
    :param rng_state:
        State dictionary of the surface-hopping random number generator.
    :param tracker_reference:
        Tracked-order descriptor reference at this frame.
    :param detector_history:
        Gap histories and duplicate-event bookkeeping at this frame.
    :param state_permutation:
        Mapping from persistent character labels to raw roots at this frame.
    :param current_frame_in_detector:
        Whether the detector snapshot already contains this frame.
    :param provider_reference:
        The provider's accepted electronic reference at this frame.  A
        backend that measures a cross-geometry descriptor needs the same
        accepted/trial chronology as the tracker, so the reference is
        captured and restored together with it.
    """

    step: int
    time_fs: float
    positions_nm: np.ndarray
    velocities_nm_ps: np.ndarray
    box_vectors: object
    active_state: int
    electronic_result: object
    rng_state: object
    tracker_reference: object
    detector_history: object
    state_permutation: np.ndarray
    current_frame_in_detector: bool
    provider_reference: object = None


class LandauZenerSurfaceHoppingDynamics:
    """
    High-level controller for adiabatic Landau-Zener surface-hopping dynamics.

    :param openmm_dynamics:
        A configured :class:`OpenMMDynamics` instance whose system, topology
        and QM force have already been built, e.g. through
        :func:`OpenMMDynamics.create_qmmm_system_from_files`.
    :param provider:
        An :class:`ElectronicStateProvider` supplying state energies,
        gradients and tracking descriptors.
    :param settings:
        A :class:`SurfaceHoppingSettings`, or a dictionary of surface-hopping
        keywords.
    :param ostream:
        Output stream for human-readable logging.
    """

    _ACTIVE_STATE_INDEX_SPACE = 'adiabatic_raw'

    def __init__(self, openmm_dynamics, provider, settings=None, ostream=None):

        if isinstance(settings, SurfaceHoppingSettings):
            self.settings = settings
            self.settings.validate()
        else:
            self.settings = SurfaceHoppingSettings.from_dict(settings)

        assert_msg_critical(
            isinstance(provider, ElectronicStateProvider),
            'LandauZenerSurfaceHoppingDynamics: provider must be an '
            'ElectronicStateProvider.')

        assert_msg_critical(
            provider.number_of_states == self.settings.number_of_states,
            'LandauZenerSurfaceHoppingDynamics: provider and settings '
            'disagree on number_of_states.')

        self.md = openmm_dynamics
        self.provider = provider

        self.ostream = ostream if ostream is not None else OutputStream(sys.stdout)

        self._scf_stability_mode = str(
            self.settings.scf_stability_mode).strip().lower()
        if self._scf_stability_mode != 'off':
            configure_stability = getattr(
                self.provider, 'configure_reference_stability_policy', None)
            assert_msg_critical(
                callable(configure_stability),
                'LandauZenerSurfaceHoppingDynamics: the selected provider '
                'does not implement the requested SCF reference-stability '
                'policy.')
            configure_stability(True)

        continuity_threshold = (
            self.settings.scf_reference_continuity_min_singular_value)
        configure_continuity = getattr(
            self.provider, 'configure_reference_continuity_policy', None)
        if callable(configure_continuity):
            configure_continuity(continuity_threshold)
        elif continuity_threshold is not None:
            assert_msg_critical(
                False,
                'LandauZenerSurfaceHoppingDynamics: the selected provider '
                'does not implement the requested SCF-reference continuity '
                'threshold.')

        self._validate_tracking_policy()

        tracker_options = {}
        if str(self.settings.state_tracking_method).lower() == \
                'backend_overlap':
            tracker_options['maximum_ambiguity_ratio'] = float(
                self.settings.maximum_tracking_ambiguity)

        self.tracker = create_state_tracker(
            self.settings.state_tracking_method,
            minimum_overlap=self.settings.minimum_tracking_overlap,
            **tracker_options)

        # The controller saves exactly one frame per nuclear interval and
        # applies a hop at that frame.  The central frame of a gap window of
        # length L lies L // 2 samples behind the newest sample, so the saved
        # frame coincides with the detector's central frame only for L = 3.
        # A longer window would require rolling back several frames, which the
        # single-frame checkpoint-and-replay chronology cannot express; refuse
        # loudly instead of hopping at the wrong geometry.
        assert_msg_critical(
            int(self.settings.gap_fit_points) == 3,
            'LandauZenerSurfaceHoppingDynamics: gap_fit_points must be 3. '
            'The checkpoint-and-replay controller rolls back a single nuclear '
            'frame, which is the central frame of a three-point gap window '
            'only. A five-point window peaks two frames back and would apply '
            'the hop at the wrong geometry; use LandauZenerEventDetector with '
            'the five-point fitter for offline gap analysis instead.')

        self.detector = LandauZenerEventDetector(
            timestep=self.settings.timestep_in_au,
            fitter=('five_point'
                    if self.settings.gap_fit_points == 5 else 'three_point'),
            gap_screening_threshold=self.settings.gap_screening_threshold,
            minimum_curvature=self.settings.minimum_curvature,
            exponent_clip=self.settings.maximum_lz_exponent,
            minimum_fit_quality=self.settings.minimum_fit_quality,
            multistate_warning_gap=self.settings.multistate_warning_gap)

        self.rescaler = SurfaceHopMomentumRescaler(
            mode=self.settings.momentum_rescaling,
            remove_translation=self.settings.remove_translation,
            remove_rotation=self.settings.remove_rotation,
            energy_tolerance=self.settings.hop_energy_tolerance)

        self.rng = np.random.default_rng(self.settings.random_seed)

        self.active_state = int(self.settings.initial_active_state)
        self.current_step = 0
        self.diagnostics = []
        self._state_permutation = np.arange(self.settings.number_of_states)
        self._tracked_spin_multiplicities = None
        self._current_frame_in_detector = False

        self._diagnostics_handle = None
        self._diagnostics_path = None
        self._diagnostics_started = False
        self._restart_diagnostics_path = None
        self._restart_diagnostics_offset = None
        self._n_hops_accepted = 0
        self._n_hops_frustrated = 0
        self._n_hops_rejected = 0
        self.reference_stability_diagnostics = []
        self.reference_continuity_diagnostics = []
        self._scf_stability_check_count = 0
        self._scf_stability_initial_completed = (
            self._scf_stability_mode == 'off')
        self._scf_stability_last_accepted_step = None
        self._recorded_scf_stability_checks = set()
        self._recorded_reference_continuity_checks = set()
        self._electronic_provenance_calculation_index = 0
        # While one nuclear interval is provisional, keep the last fully
        # accepted frame available for fail-closed rollback.  In particular,
        # an overlap/tracking failure occurs only after OpenMM has already
        # propagated to the rejected endpoint; checkpointing that endpoint
        # together with the previous electronic reference would create a
        # nuclear/electronic mismatch on restart.
        self._abort_frame = None

    def _validate_tracking_policy(self):
        """
        Enforces the production state-tracking policy before anything runs.

        A production spin-flip backend exposes neither an oscillator strength
        nor a transition dipole, so the descriptor tracker collapses onto
        continuity of the excitation energy and the identity tracker is the
        raw energy ordering by construction.  Using either as the sole
        state-identity evidence is the silent energy-order fallback that this
        controller must never perform, so it is refused unless the run is
        explicitly marked as unsafe diagnostics.
        """

        method = str(self.settings.state_tracking_method).lower()
        backend = str(self.settings.electronic_backend).lower()
        production_backend = (
            backend in SurfaceHoppingSettings.PRODUCTION_BACKENDS)
        production_provider = bool(self.provider.is_production_ready())

        if production_backend:
            assert_msg_critical(
                production_provider,
                'LandauZenerSurfaceHoppingDynamics: electronic_backend '
                f'{backend} was requested but the supplied provider does not '
                'implement a production backend. Refusing to run native '
                'VeloxChem TDA under a spin-flip backend label.')

        if not (production_backend or production_provider):
            return

        assert_msg_critical(
            is_production_tracker(method) or
            bool(self.settings.allow_unsafe_state_tracking),
            'LandauZenerSurfaceHoppingDynamics: production surface hopping '
            f'refuses state_tracking_method={method!r}. Neither spin-flip '
            'backend exposes an oscillator strength or a transition dipole, '
            'so this tracker would assign state identity from excitation-'
            "energy order alone. Use 'backend_overlap'.")

        if not is_production_tracker(method):
            self.ostream.print_warning(
                'UNSAFE STATE TRACKING: allow_unsafe_state_tracking is set, '
                f'so persistent identities come from {method!r} tracking, '
                'which for this backend is excitation-energy continuity '
                'alone. This run is diagnostics only and must not be used '
                'for production data.')
            self.ostream.flush()

    def _validate_serenity_spin_policy(self, result, tracking):
        """Validates spin sectors without renumbering Serenity raw roots.

        Serenity's SF response contains singlet, triplet and sometimes
        spin-incomplete roots in one energy-ordered ladder. The ladder and
        its derivative selectors remain untouched, but without spin-orbit
        coupling a production trajectory may only connect adiabatic roots of
        the same multiplicity. The first accepted frame also fixes the spin
        sector of every persistent character label; later frames must carry
        the same sectors through their descriptor-based permutation.

        :return:
            Proposed character-ordered multiplicities, or ``None`` for a
            backend whose target manifold is already spin-pure.
        """

        snapshot = getattr(result, 'snapshot', None)
        if snapshot is None or str(snapshot.backend).lower() != 'serenity':
            return None

        metadata = snapshot.spin_metadata
        if not isinstance(metadata, dict):
            raise SurfaceHoppingError(
                'Serenity SF returned no per-root spin metadata; production '
                'dynamics cannot exclude spin-incompatible state identity or '
                'hopping without spin-orbit coupling.')

        n_states = int(self.settings.number_of_states)
        multiplicities = np.asarray(
            metadata.get('state_multiplicities', []), dtype=int).reshape(-1)
        deviations = np.asarray(metadata.get('s2_deviation', []),
                                dtype=float).reshape(-1)
        tolerance = float(metadata.get('s2_tolerance', 0.5))

        if (multiplicities.size != n_states or
                deviations.size != n_states or
                not np.all(np.isfinite(deviations))):
            raise SurfaceHoppingError(
                'Serenity SF returned incomplete or nonfinite per-root spin '
                'metadata; production dynamics is refused.')

        invalid = np.flatnonzero(deviations > tolerance)
        if invalid.size:
            roots = (invalid + 1).tolist()
            raise SurfaceHoppingError(
                'Serenity SF target root(s) '
                f'{roots} are spin-incomplete: their <S^2> deviations exceed '
                f'the configured tolerance {tolerance:.6f}. Increase the '
                'response window or choose a compatible target ladder; these '
                'roots cannot enter production surface hopping.')

        permutation = np.asarray(tracking.permutation, dtype=int)
        tracked = multiplicities[permutation]

        incompatible = set()
        for source in range(n_states):
            for target in self.settings.candidate_states(source):
                if multiplicities[source] != multiplicities[target]:
                    incompatible.add(tuple(sorted((source, target))))

        if incompatible:
            pairs = sorted(incompatible)
            raise SurfaceHoppingError(
                'Serenity SF allowed_transitions contains spin-changing '
                f'pair(s) {pairs} in adiabatic-root space, but no spin-orbit '
                'coupling is implemented. Use an explicit list containing '
                'only same-multiplicity adiabatic-root pairs.')

        if (self._tracked_spin_multiplicities is not None and
                not np.array_equal(tracked,
                                   self._tracked_spin_multiplicities)):
            raise SurfaceHoppingError(
                'Serenity SF state tracking assigned a persistent character '
                'label to an incompatible spin sector: accepted '
                'multiplicities '
                f'{self._tracked_spin_multiplicities.tolist()}, current '
                f'{tracked.tolist()}. The electronic step is refused rather '
                'than falling back to energy order.')

        return np.array(tracked, dtype=int, copy=True)

    # ------------------------------------------------------------------
    # OpenMM state access
    # ------------------------------------------------------------------

    @property
    def context(self):
        """
        :return:
            The OpenMM context of the underlying simulation.
        """

        assert_msg_critical(
            getattr(self.md, 'simulation', None) is not None,
            'LandauZenerSurfaceHoppingDynamics: the OpenMM simulation has '
            'not been created.')

        return self.md.simulation.context

    @property
    def qm_atoms(self):
        """
        :return:
            List of OpenMM particle indices belonging to the QM region.
        """

        return list(self.md.qm_atoms)

    @property
    def rescaling_atoms(self):
        """
        :return:
            List of OpenMM particle indices whose velocities are rescaled on a
            hop.
        """

        spec = self.settings.rescaling_atoms

        if isinstance(spec, str):
            assert_msg_critical(
                spec.lower() == 'qm',
                'LandauZenerSurfaceHoppingDynamics: rescaling_atoms must be '
                '"qm" or an explicit list of indices.')
            return self.qm_atoms

        return [int(idx) for idx in spec]

    def get_masses_in_au(self, indices):
        """
        Reads particle masses from the OpenMM system.

        :param indices:
            Iterable of particle indices.

        :return:
            Array of shape (n,) with masses in electron masses.
        """

        masses = [
            self.md.system.getParticleMass(int(idx)).value_in_unit(unit.dalton)
            for idx in indices
        ]

        return np.asarray(masses, dtype=float) * AMU_IN_ELECTRON_MASS

    def get_qm_geometry_in_bohr(self):
        """
        :return:
            Array of shape (n_qm, 3) with the QM-region geometry in bohr.
        """

        state = self.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True).value_in_unit(
            unit.nanometer)

        return np.asarray(positions[self.qm_atoms], dtype=float) * NM_IN_BOHR

    def save_frame(self, step, electronic_result):
        """
        Captures the complete OpenMM and surface-hopping state at a frame.

        :param step:
            Integer frame index.
        :param electronic_result:
            The :class:`ElectronicStateResult` at this frame.

        :return:
            A :class:`TrajectoryFrame`.
        """

        state = self.context.getState(getPositions=True, getVelocities=True)

        try:
            box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(
                unit.nanometer)
            box = np.asarray(box, dtype=float)
        except Exception:
            box = None
        # yapf: disable
        return TrajectoryFrame(
            step=int(step),
            time_fs=float(state.getTime().value_in_unit(unit.femtosecond)),
            positions_nm=np.asarray(
                state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
                dtype=float),
            velocities_nm_ps=np.asarray(
                state.getVelocities(asNumpy=True).value_in_unit(
                    unit.nanometer / unit.picosecond),
                dtype=float),
            box_vectors=box,
            active_state=int(self.active_state),
            electronic_result=electronic_result,
            rng_state=copy.deepcopy(self.rng.bit_generator.state),
            tracker_reference=copy.deepcopy(self.tracker._reference),
            detector_history=copy.deepcopy(
                self.detector.get_history_snapshot()),
            state_permutation=np.array(self._state_permutation, copy=True),
            current_frame_in_detector=bool(self._current_frame_in_detector),
            # Electronic snapshots are immutable, so the accepted reference is
            # captured by reference rather than copied.
            provider_reference=self.provider.get_reference_state())
        # yapf: enable

    def restore_frame(self, frame):
        """
        Restores a previously saved frame into the OpenMM context.

        :param frame:
            A :class:`TrajectoryFrame`.
        """

        try:
            if frame.box_vectors is not None:
                self.context.setPeriodicBoxVectors(*[
                    mm.Vec3(*row) * unit.nanometer for row in frame.box_vectors
                ])

            self.context.setPositions(frame.positions_nm * unit.nanometer)
            self.context.setVelocities(frame.velocities_nm_ps * unit.nanometer /
                                       unit.picosecond)
            self.context.setTime(frame.time_fs * unit.femtosecond)
        except Exception as exc:
            raise SurfaceHoppingError(
                f'failed to restore OpenMM state at step {frame.step}: {exc}')

        self.active_state = int(frame.active_state)
        self.rng.bit_generator.state = copy.deepcopy(frame.rng_state)
        self.tracker._reference = copy.deepcopy(frame.tracker_reference)
        self.detector.set_history_snapshot(copy.deepcopy(
            frame.detector_history))
        self._state_permutation = np.array(frame.state_permutation, copy=True)
        self._current_frame_in_detector = bool(frame.current_frame_in_detector)
        # The discarded trial endpoint must not remain the electronic
        # reference the replayed endpoint is measured against.
        self.provider.set_reference_state(frame.provider_reference)
        self.install_force(frame.electronic_result.active_gradient)

    def install_force(self, gradient):
        """
        Writes the active-state force into the QM ``CustomExternalForce``.

        Uses the same sign and unit convention as
        :func:`OpenMMDynamics.update_forces`: the force is ``-grad(E_a)``
        converted from Hartree/bohr to kJ/mol/nm.

        :param gradient:
            Array of shape (n_qm, 3) with the gradient in Hartree/bohr.
        """

        gradient = np.asarray(gradient, dtype=float)

        assert_msg_critical(
            np.all(np.isfinite(gradient)),
            'LandauZenerSurfaceHoppingDynamics: nonfinite gradient cannot be '
            'installed as an OpenMM force.')

        conversion = hartree_in_kjpermol() * 10.0 / bohr_in_angstrom()
        force_array = -gradient * conversion

        custom_force = self.md.system.getForce(self.md.qm_force_index)

        for i, atom_idx in enumerate(self.qm_atoms):
            custom_force.setParticleParameters(i, atom_idx, force_array[i])

        custom_force.updateParametersInContext(self.context)

    # ------------------------------------------------------------------
    # electronic evaluation and tracking
    # ------------------------------------------------------------------

    def _reference_stability_reason(self, step):
        """Returns the policy reason for checking this accepted-step index."""

        mode = getattr(self, '_scf_stability_mode', 'off')
        if mode == 'off':
            return None
        if not getattr(self, '_scf_stability_initial_completed', False):
            return 'initial'
        if (mode == 'periodic' and int(step) > 0 and
                int(step) % int(self.settings.scf_stability_interval) == 0 and
                getattr(self, '_scf_stability_last_accepted_step', None) !=
                int(step)):
            return 'periodic'
        return None

    def _request_reference_stability(self, step):
        """Stages one due stability check on the OpenQP backend provider."""

        reason = self._reference_stability_reason(step)
        if reason is None:
            return None

        request_stability = getattr(
            self.provider, 'request_reference_stability', None)
        assert_msg_critical(
            callable(request_stability),
            'LandauZenerSurfaceHoppingDynamics: the provider cannot accept an '
            'SCF reference-stability request.')

        try:
            state = self.context.getState()
            time_fs = float(
                state.getTime().value_in_unit(unit.femtosecond))
        except Exception:
            time_fs = float(step) * float(self.settings.electronic_timestep)

        self._scf_stability_check_count += 1
        request = {
            'stability_mode': self._scf_stability_mode,
            'stability_requested': True,
            'reason_for_check': reason,
            'step': int(step),
            'time': time_fs,
            'check_index': int(self._scf_stability_check_count),
            'diagnostics_dir':
                str(self.settings.scf_stability_diagnostics_dir),
        }
        request_stability(request)
        return request

    def _provider_reference_stability_record(self):
        """Returns the latest backend record after an interrupted evaluation."""

        getter = getattr(
            self.provider, 'last_reference_stability_record', None)
        if callable(getter):
            return getter()
        return None

    def _provider_reference_continuity_record(self):
        """Returns the latest backend continuity record after a failure."""

        getter = getattr(
            self.provider, 'last_reference_continuity_record', None)
        if callable(getter):
            return getter()
        return None

    def _write_reference_stability_diagnostics(self, record, error=None):
        """Persists one requested-check event in the trajectory JSONL."""

        if record is None:
            return
        payload = copy.deepcopy(record)
        check_index = int(payload.get('stability_check_index', 0))
        if check_index in self._recorded_scf_stability_checks:
            return

        payload.update({
            'record_type': 'scf_reference_stability',
            'trajectory_id': self.settings.trajectory_id,
            'reference_compatible': bool(
                payload.get('reason_for_check') == 'initial' or
                not payload.get('reference_changed', False)),
            'eligible_for_lz_history': bool(
                payload.get('stability_performed', False) and
                payload.get('stability_converged', False) and
                (payload.get('reason_for_check') == 'initial' or
                 not payload.get('reference_changed', False))),
            'failure': None if error is None else str(error),
        })
        self._recorded_scf_stability_checks.add(check_index)
        self.reference_stability_diagnostics.append(payload)
        self._write_diagnostics_record(payload)

    def _record_reference_continuity(self, record, error=None):
        """Retains one compact continuity decision without changing tracking."""

        if record is None:
            return
        payload = copy.deepcopy(record)
        key = (payload.get('previous_reference_identifier'),
               payload.get('current_reference_identifier'))
        if key in self._recorded_reference_continuity_checks:
            return
        self._recorded_reference_continuity_checks.add(key)
        self.reference_continuity_diagnostics.append(payload)

        # Successful records are part of the ordinary per-step diagnostics and
        # accepted-step archive.  A rejected calculation has neither, so write
        # a dedicated JSONL event before aborting.
        if error is not None:
            payload.update({
                'record_type': 'scf_reference_continuity_rejection',
                'trajectory_id': self.settings.trajectory_id,
                'trajectory_step': int(self.current_step),
                'eligible_for_lz_history': False,
                'failure': str(error),
            })
            self._write_diagnostics_record(payload)

    def evaluate_electronic_state(self, step):
        """
        Evaluates the electronic states at the current geometry, applies state
        tracking, and returns the tracked result.

        :param step:
            Integer frame index, used for diagnostics.

        :return:
            The tuple (result, tracking_result).
        """

        geometry = self.get_qm_geometry_in_bohr()
        accepted_tracker_reference = copy.deepcopy(self.tracker._reference)
        get_provider_reference = getattr(
            self.provider, 'get_reference_state', None)
        accepted_provider_reference = (
            get_provider_reference() if callable(get_provider_reference) else
            None)
        accepted_permutation = np.array(self._state_permutation, copy=True)
        # ``active_state`` is an adiabatic raw-root index.  Character tracking
        # must never redirect the gradient: doing so switches potential-energy
        # surfaces without a Landau-Zener hop.
        raw_state_hint = int(self.active_state)
        stability_request = self._request_reference_stability(step)

        def rollback_electronic_transaction():
            """Restores both owners of accepted electronic-state history."""

            self.provider.rollback_reference()
            set_provider_reference = getattr(
                self.provider, 'set_reference_state', None)
            if callable(set_provider_reference):
                set_provider_reference(accepted_provider_reference)
            self.tracker._reference = copy.deepcopy(
                accepted_tracker_reference)
            self._state_permutation = np.array(accepted_permutation,
                                               copy=True)

        try:
            result = self.provider.compute(geometry, raw_state_hint)
        except Exception as exc:
            # A failed electronic step commits nothing: no accepted reference,
            # no tracking history, no force.
            rollback_electronic_transaction()
            stability_record = getattr(exc, 'stability_record', None)
            if stability_record is None and stability_request is not None:
                stability_record = self._provider_reference_stability_record()
            self._write_reference_stability_diagnostics(
                stability_record, error=exc)
            continuity_record = getattr(exc, 'continuity_record', None)
            if continuity_record is None:
                continuity_record = (
                    self._provider_reference_continuity_record())
            self._record_reference_continuity(
                continuity_record, error=exc)
            raise SurfaceHoppingError(
                f'electronic-structure calculation failed at step {step}: '
                f'{exc}')

        stability_record = None
        snapshot = getattr(result, 'snapshot', None)
        if snapshot is not None:
            stability_record = getattr(snapshot, 'reference_stability', None)
        self._write_reference_stability_diagnostics(stability_record)
        continuity_record = (None if snapshot is None else
                             getattr(snapshot, 'reference_continuity', None))

        tracking = None
        try:
            if not result.scf_converged:
                raise SurfaceHoppingError(
                    f'SCF did not converge at step {step}.')
            if not result.response_converged:
                raise SurfaceHoppingError(
                    'the response calculation did not converge at step '
                    f'{step}.')
            if not result.gradient_converged:
                raise SurfaceHoppingError(
                    f'the active-state gradient failed at step {step}.')
            if not result.is_valid():
                raise SurfaceHoppingError(
                    f'nonfinite electronic energies or gradient at step '
                    f'{step}.')

            if (continuity_record is not None and
                    continuity_record.get(
                        'reference_continuity_checked', False) and
                    not continuity_record.get('reference_continuous', False)):
                raise SurfaceHoppingError(
                    'the high-spin SCF reference is discontinuous with the '
                    f'previous accepted reference at step {step} '
                    f"({continuity_record.get('reference_continuity_reason', 'unknown')}).")

            tracking = self.tracker.track(result.state_descriptors)

            if not tracking.is_reliable:
                detail = '; '.join(tracking.warnings) or (
                    f'confidence {tracking.confidence:.3f} is below the '
                    f'threshold {self.settings.minimum_tracking_overlap:.3f}')
                raise SurfaceHoppingError(
                    f'state tracking became ambiguous at step {step} '
                    f'({tracking.method}/{tracking.descriptor}): {detail}.')

            proposed_spin_multiplicities = (
                self._validate_serenity_spin_policy(result, tracking))
        except SurfaceHoppingError as exc:
            self._archive_rejected_electronic_trial(
                step, result, tracking, exc)
            rollback_electronic_transaction()
            if (continuity_record is not None and
                    not continuity_record.get('reference_continuous', True)):
                self._record_reference_continuity(
                    continuity_record, error=exc)
            raise

        self._record_reference_continuity(continuity_record)

        result.state_permutation = tracking.permutation
        result.tracking_scores = tracking.scores
        self._state_permutation = np.asarray(tracking.permutation, dtype=int)

        active_raw_state = int(self.active_state)

        if int(result.active_raw_state) != active_raw_state:
            try:
                result.active_gradient = np.asarray(
                    self.provider.compute_state_gradient(
                        geometry, active_raw_state),
                    dtype=float)
            except Exception as exc:
                rollback_electronic_transaction()
                raise SurfaceHoppingError(
                    'failed to obtain the active adiabatic-state gradient at '
                    f'step {step}: raw root {active_raw_state}: '
                    f'{exc}')

        result.active_state = active_raw_state
        result.active_raw_state = active_raw_state
        result.active_tracked_state = StateIndexConverter.raw_to_tracked(
            active_raw_state, self._state_permutation)
        result.derivative_state_index = result.selector_for(active_raw_state)
        result.gradient_converged = bool(
            np.all(np.isfinite(result.active_gradient)))

        if not result.gradient_converged:
            rollback_electronic_transaction()
            raise SurfaceHoppingError(
                f'the active adiabatic-state gradient failed at step {step}.')

        # Commit only after the active adiabatic-root gradient is finite and
        # belongs to the same validated snapshot. Until this point the
        # provider snapshot and the tracker's permutation are a single trial
        # transaction; a failure in either half restores both accepted
        # histories.
        try:
            self.provider.commit_reference()
        except Exception as exc:
            rollback_electronic_transaction()
            raise SurfaceHoppingError(
                f'failed to commit the electronic state at step {step}: '
                f'{exc}')

        archive = getattr(
            self.provider, 'archive_accepted_provenance', None)
        provenance_path = None
        self._electronic_provenance_calculation_index = int(getattr(
            self, '_electronic_provenance_calculation_index', 0)) + 1
        if callable(archive):
            try:
                provenance_path = archive(
                    step=int(step),
                    result=result,
                    tracking=tracking,
                    directory=self.settings.electronic_provenance_dir,
                    calculation_index=
                        self._electronic_provenance_calculation_index,
                    wavefunction_retention=
                        self.settings
                        .electronic_provenance_wavefunction_steps)
            except Exception as exc:
                rollback_electronic_transaction()
                raise SurfaceHoppingError(
                    'failed to archive the accepted electronic provenance at '
                    f'step {step}: {exc}')
        result.electronic_provenance_path = provenance_path

        if proposed_spin_multiplicities is not None:
            self._tracked_spin_multiplicities = proposed_spin_multiplicities

        if (stability_record is not None and
                stability_record.get('reason_for_check') == 'initial'):
            self._scf_stability_initial_completed = True

        return result, tracking

    def _archive_rejected_electronic_trial(
            self, step, result, tracking, error):
        """Archives a failed trial before its provider transaction is rolled back."""

        archive = getattr(
            self.provider, 'archive_rejected_provenance', None)
        if not callable(archive) or getattr(result, 'snapshot', None) is None:
            return None

        history = self.detector.get_history_snapshot()
        pairs = set(history.get('gap_history', {}))
        pairs.update(history.get('step_history', {}))
        serializable_history = {
            'pairs': [
                {
                    'states': [int(value) for value in pair],
                    'gaps': [float(value) for value in
                             history.get('gap_history', {}).get(pair, ())],
                    'steps': [int(value) for value in
                              history.get('step_history', {}).get(pair, ())],
                }
                for pair in sorted(pairs)
            ],
            'last_event_steps': [
                {
                    'states': [int(value) for value in pair],
                    'step': int(value),
                }
                for pair, value in sorted(
                    history.get('last_event_step', {}).items())
            ],
        }

        self._electronic_provenance_calculation_index = int(getattr(
            self, '_electronic_provenance_calculation_index', 0)) + 1
        try:
            path = archive(
                step=int(step),
                result=result,
                tracking=tracking,
                directory=self.settings.electronic_provenance_dir,
                calculation_index=(
                    self._electronic_provenance_calculation_index),
                lz_history_before_rejection=serializable_history,
                wavefunction_retention=(
                    self.settings.electronic_provenance_wavefunction_steps))
        except Exception as archive_error:
            self._write_diagnostics_record({
                'record_type': 'rejected_electronic_provenance_error',
                'trajectory_id': self.settings.trajectory_id,
                'trajectory_step': int(step),
                'accepted': False,
                'eligible_for_lz': False,
                'failure': str(error),
                'archive_failure': str(archive_error),
            })
            return None

        self._write_diagnostics_record({
            'record_type': 'rejected_electronic_trial',
            'trajectory_id': self.settings.trajectory_id,
            'trajectory_step': int(step),
            'accepted': False,
            'eligible_for_lz': False,
            'active_raw_state': int(result.active_raw_state),
            'active_tracked_state': int(result.active_tracked_state),
            'failure': str(error),
            'electronic_provenance_path': path,
            'lz_history_before_rejection': serializable_history,
        })
        return path

    # ------------------------------------------------------------------
    # hopping
    # ------------------------------------------------------------------

    def select_event(self, tracking):
        """
        Runs the detector and returns the highest-priority valid event.

        The deterministic selection policy is: discard invalid events, sort by
        fitted event time, break ties by smallest minimum gap, and evaluate
        only the first event.  At most one hop is processed per nuclear
        interval.

        :param tracking:
            The :class:`TrackingResult` of the current frame.

        :return:
            The tuple (event, n_candidates); event is None when no valid
            candidate exists.
        """

        events = self.detector.detect(
            tracking_quality=tracking.confidence,
            minimum_tracking_quality=self.settings.minimum_tracking_overlap)

        if not events:
            return None, 0

        return events[0], len(events)

    def attempt_hop(self, event, frame):
        """
        Attempts a hop at the central frame.

        :param event:
            The :class:`LandauZenerEvent` being processed.
        :param frame:
            The :class:`TrajectoryFrame` at the central frame ``n``.

        :return:
            A dictionary describing the outcome, with keys ``accepted``,
            ``frustrated``, ``kinetic_before``, ``kinetic_after``,
            ``energy_residual`` and ``message``.
        """

        result = frame.electronic_result
        adiabatic_energies = result.adiabatic_state_energies()

        source = int(event.source_state)
        target = int(event.target_state)
        source_raw = source
        target_raw = target
        delta_energy = float(adiabatic_energies[target_raw] -
                             adiabatic_energies[source_raw])

        indices = self.rescaling_atoms
        masses = self.get_masses_in_au(indices)

        positions = frame.positions_nm[indices] * NM_IN_BOHR
        velocities = (frame.velocities_nm_ps[indices] *
                      NM_PER_PS_IN_AU_VELOCITY)

        direction = None
        if self.rescaler.requires_target_gradient():
            direction = self._gap_gradient_direction(frame, source, target,
                                                     indices)

        rescaling = self.rescaler.rescale(masses, positions, velocities,
                                          delta_energy, direction)

        outcome = {
            'accepted': bool(rescaling.accepted),
            'frustrated': bool(rescaling.frustrated),
            'kinetic_before': float(rescaling.kinetic_before),
            'kinetic_after': float(rescaling.kinetic_after),
            'energy_residual': float(rescaling.energy_residual),
            'scaling_factor': rescaling.scaling_factor,
            'gamma': rescaling.gamma,
            'delta_energy': delta_energy,
            'source_raw_state': source_raw,
            'target_raw_state': target_raw,
            'target_derivative_state_index': result.selector_for(target_raw),
            'message': rescaling.message,
        }

        if not rescaling.accepted:
            return outcome

        # Write the rescaled velocities back into the restored frame.
        new_velocities = np.array(frame.velocities_nm_ps, copy=True)
        new_velocities[indices] = (rescaling.velocities *
                                   AU_VELOCITY_IN_NM_PER_PS)

        try:
            self.context.setVelocities(new_velocities * unit.nanometer /
                                       unit.picosecond)
        except Exception as exc:
            raise SurfaceHoppingError(
                f'failed to apply rescaled velocities at step {frame.step}: '
                f'{exc}')

        self.active_state = target

        return outcome

    def _gap_gradient_direction(self, frame, source, target, indices):
        """
        Builds the gap-gradient direction ``grad(E_b) - grad(E_a)`` restricted
        to the rescaling atom set.

        Only evaluated when ``gap_gradient`` rescaling is enabled and a hop
        has already been proposed.

        :param frame:
            The :class:`TrajectoryFrame` at the central frame.
        :param source:
            Adiabatic raw index of the source state.
        :param target:
            Adiabatic raw index of the target state.
        :param indices:
            Particle indices of the rescaling atom set.

        :return:
            Array of shape (n, 3) with the direction in Hartree/bohr.
        """

        geometry = frame.electronic_result.geometry
        source_raw = int(source)
        target_raw = int(target)

        assert_msg_critical(
            int(frame.electronic_result.active_raw_state) == source_raw,
            'LandauZenerSurfaceHoppingDynamics: the saved source gradient '
            'does not belong to the adiabatic source state.')

        # yapf: disable
        try:
            grad_source = np.asarray(frame.electronic_result.active_gradient,
                                     dtype=float)
            grad_target = np.asarray(
                self.provider.compute_state_gradient(geometry, target_raw),
                dtype=float)
        except NotImplementedError as exc:
            raise SurfaceHoppingError(
                f'gap_gradient rescaling is not supported by this provider: '
                f'{exc}')
        # yapf: enable

        direction_qm = grad_target - grad_source

        qm_atoms = self.qm_atoms
        lookup = {atom: i for i, atom in enumerate(qm_atoms)}

        assert_msg_critical(
            all(idx in lookup for idx in indices),
            'LandauZenerSurfaceHoppingDynamics: gap_gradient rescaling '
            'requires the rescaling atoms to be a subset of the QM atoms.')

        return np.asarray([direction_qm[lookup[idx]] for idx in indices],
                          dtype=float)

    # ------------------------------------------------------------------
    # main loop
    # ------------------------------------------------------------------

    def run(self,
            nsteps,
            checkpoint_file=None,
            diagnostics_file=None,
            traj_file=None,
            save_freq=1,
            print_freq=1):
        """
        Runs a surface-hopping trajectory and writes only retained endpoints
        to the optional PDB trajectory.

        :param nsteps:
            Number of nuclear steps to propagate.
        :param checkpoint_file:
            Overrides :attr:`SurfaceHoppingSettings.checkpoint_file`.
        :param diagnostics_file:
            Overrides :attr:`SurfaceHoppingSettings.diagnostics_file`.
        :param traj_file:
            Path to the PDB trajectory file. If provided, accepted endpoints
            are written at ``save_freq`` frequency.
        :param save_freq:
            Frequency (in steps) for writing PDB frames. Default: 1 (every step).
        :param print_freq:
            Frequency (in steps) for formatted console output. Default: 1 (every step).

        :return:
            The list of per-step diagnostics dictionaries.
        """
        checkpoint_file = checkpoint_file or self.settings.checkpoint_file
        diagnostics_file = (diagnostics_file or self.settings.diagnostics_file)

        # Keep trajectory output under the controller's transaction.  An
        # OpenMM PDBReporter attached to ``simulation.reporters`` fires during
        # every provisional trial/replay step and therefore archives frames
        # that the electronic transaction later discards.  Write only after
        # _advance_one_step() has accepted the endpoint.
        trajectory_handle = None
        trajectory_app = None
        if traj_file:
            try:
                import openmm.app as app
                trajectory_handle = open(traj_file, 'w', encoding='utf-8')
                trajectory_app = app
                app.PDBFile.writeHeader(self.md.topology, trajectory_handle)
                self.ostream.print_info(
                    f"PDB trajectory will be written to: {traj_file}")
            except ImportError:
                self.ostream.print_warning(
                    'PDBReporter requires OpenMM.app module. '
                    'Trajectory writing disabled.')
            except Exception as exc:
                self.ostream.print_warning(
                    f'Failed to initialize PDB trajectory output: {exc}. '
                    'Trajectory writing disabled.')
                if trajectory_handle is not None:
                    trajectory_handle.close()
                    trajectory_handle = None
                    trajectory_app = None

        self._check_ensemble()
        self._check_nuclear_timestep()
        self._open_diagnostics(diagnostics_file)

        try:
            result, tracking = self.evaluate_electronic_state(self.current_step)

            if not self._current_frame_in_detector:
                self.detector.push(
                    self.current_step, self.active_state,
                    result.adiabatic_state_energies(),
                    self.settings.candidate_states(self.active_state))
                self._current_frame_in_detector = True

            # Initial snapshot output
            if print_freq > 0:
                initial_record = self._new_record(-1, None, result, tracking,
                                                  None, 0)
                self.print_snapshot(initial_record)

            for _ in range(int(nsteps)):
                result, tracking = self._advance_one_step(
                    result, tracking, checkpoint_file)

                record = self.diagnostics[-1]

                if (trajectory_handle is not None and
                        self.current_step % int(save_freq) == 0):
                    try:
                        endpoint = self.context.getState(getPositions=True)
                        trajectory_app.PDBFile.writeModel(
                            self.md.topology,
                            endpoint.getPositions(),
                            trajectory_handle,
                            modelIndex=self.current_step)
                        trajectory_handle.flush()
                    except Exception as exc:
                        self.ostream.print_warning(
                            f'Failed to write PDB trajectory endpoint: {exc}. '
                            'Trajectory writing disabled.')
                        trajectory_handle.close()
                        trajectory_handle = None
                        trajectory_app = None

                # Periodic formatted output
                if print_freq > 0 and self.current_step % print_freq == 0:
                    self.print_snapshot(record)

        except SurfaceHoppingError:
            if self._abort_frame is not None:
                abort_frame = self._abort_frame
                self._abort_frame = None
                self.restore_frame(abort_frame)
                self.current_step = int(abort_frame.step)
            if self.settings.checkpoint_interval:
                self._write_checkpoint(checkpoint_file, best_effort=True)
            raise
        finally:
            self._close_diagnostics()
            if trajectory_handle is not None:
                try:
                    trajectory_app.PDBFile.writeFooter(
                        self.md.topology, trajectory_handle)
                finally:
                    trajectory_handle.close()

        return self.diagnostics

    def print_snapshot(self, record):
        """
        Prints formatted snapshot information to output stream.

        :param record:
            The diagnostics record dictionary from _new_record().
        """
        step = record.get('step', self.current_step)
        time_fs = record.get('time', 0.0)
        time_ps = time_fs / 1000.0

        active_before = record.get('active_state', self.active_state)
        active_after = record.get('active_state_after', active_before)

        lines = [
            f"\n{'='*70}",
            f"  SNAPSHOT {step:6d}  |  Time: {time_ps:8.2f} ps  |  State: {active_before} -> {active_after}",
            f"{'='*70}"
        ]

        # Energy section
        if 'kinetic_energy' in record:
            lines.append(
                f"  Kinetic Energy:     {record['kinetic_energy']:12.8f} Eh")
        if 'active_state_energy' in record:
            lines.append(
                f"  Active State Energy: {record['active_state_energy']:12.8f} Eh"
            )
        if 'total_energy' in record:
            lines.append(
                f"  Total Energy:       {record['total_energy']:12.8f} Eh")
        if 'temperature_kelvin' in record:
            lines.append(
                f"  Temperature:        {record['temperature_kelvin']:8.2f} K")

        # State energies
        if 'tracked_state_energies' in record:
            energies = record['tracked_state_energies']
            lines.append("\n  State Energies (Eh):")
            for i, energy in enumerate(energies):
                marker = " <-- ACTIVE" if i == active_after else ""
                lines.append(f"    State {i}: {energy:12.8f}{marker}")

        # Energy gaps
        if 'energy_gaps' in record and len(record['energy_gaps']) > 0:
            lines.append("\n  Energy Gaps (Eh):")
            for i, gap in enumerate(record['energy_gaps']):
                lines.append(f"    Gap {i}-{i+1}: {gap:12.8f}")
            if 'minimum_energy_gap' in record and record['minimum_energy_gap']:
                lines.append(
                    f"    Min Gap:    {record['minimum_energy_gap']:12.8f}")

        # Hopping information
        if record.get('hop_proposed', False):
            lines.append(f"\n  {'-' * 66}")
            lines.append("  HOP PROPOSED:")
            lines.append(
                f"    {record.get('candidate_source_state', '?')} -> {record.get('candidate_target_state', '?')}"
            )
            lines.append(
                f"    Min Gap:       {record.get('minimum_gap', 0.0):12.8f} Eh")
            lines.append(
                f"    Gap Curvature: {record.get('gap_curvature', 0.0):12.8f} Eh/a.u.²"
            )
            lines.append(
                f"    LZ Exponent:   {record.get('landau_zener_exponent', 0.0):12.4f}"
            )
            lines.append(
                f"    LZ Probability: {record.get('landau_zener_probability', 0.0):12.6f}"
            )
            lines.append(
                f"    Random Num:    {record.get('random_number', 0.0):12.6f}")

            if record.get('hop_accepted'):
                lines.append("    -> HOP ACCEPTED")
                if 'hop_energy_residual' in record:
                    lines.append(
                        f"    Energy Residual: {record['hop_energy_residual']:12.8f} Eh"
                    )
            elif record.get('hop_frustrated'):
                lines.append(
                    f"    -> HOP FRUSTRATED: {record.get('hop_message', 'N/A')}"
                )
            else:
                lines.append("    -> HOP REJECTED")

        # Convergence flags
        lines.append("\n  Convergence:")
        lines.append(f"    SCF: {record.get('scf_converged', True)} | "
                     f"Response: {record.get('response_converged', True)} | "
                     f"Gradient: {record.get('gradient_converged', True)}")

        lines.append(f"{'='*70}\n")
        output = "\n".join(lines)
        self.ostream.print_info(output)
        self.ostream.flush()

    def write_trajectory_pdb(self,
                             traj_file,
                             start_step=None,
                             end_step=None,
                             step_stride=1,
                             include_all_states=False):
        """
        Writes a PDB trajectory file from the stored frames.

        :param traj_file:
            Path to the output PDB trajectory file.
        :param start_step:
            First step to include (default: 0).
        :param end_step:
            Last step to include (default: final step).
        :param step_stride:
            Stride for selecting steps (default: 1, every step).
        :param include_all_states:
            If True, writes a separate PDB file for each state with the geometry
            at which that state was active (default: False).
        """
        try:
            import openmm.app as app
        except ImportError:
            self.ostream.print_warning(
                'PDB trajectory writing requires OpenMM.app module')
            return

        if start_step is None:
            start_step = 0
        if end_step is None:
            end_step = self.current_step

        # For now, we'll write the current trajectory from the simulation
        # In a full implementation, you would need to store all positions
        # or be able to regenerate them from checkpoints

        try:
            with open(traj_file, 'w') as f:
                for step in range(start_step, end_step + 1, step_stride):
                    # Get current state
                    state = self.context.getState(getPositions=True)
                    positions = state.getPositions()

                    # Write frame
                    app.PDBFile.writeFile(
                        self.md.topology,
                        positions,
                        f,
                        append=True if step > start_step else False)

            self.ostream.print_info(f"PDB trajectory written to {traj_file}")

        except Exception as exc:
            self.ostream.print_warning(f'Failed to write PDB trajectory: {exc}')

    def _advance_one_step(self, result, tracking, checkpoint_file):
        """
        Propagates one nuclear step with checkpoint-and-replay hop handling.

        :param result:
            The :class:`ElectronicStateResult` at the current frame.
        :param tracking:
            The :class:`TrackingResult` at the current frame.
        :param checkpoint_file:
            Path of the checkpoint file.

        :return:
            The tuple (result, tracking) at the new current frame.
        """

        step = self.current_step

        # 1. save the complete state at frame n
        frame = self.save_frame(step, result)
        self._abort_frame = frame

        # 2. trial propagation from n to n+1 on the current surface
        self.install_force(result.active_gradient)
        self.md.simulation.step(1)
        self.current_step = step + 1

        # 3. electronic states at n+1
        trial_result, trial_tracking = self.evaluate_electronic_state(
            self.current_step)

        # 4. gap history and event detection; the central frame is n
        self.detector.push(self.current_step, self.active_state,
                           trial_result.adiabatic_state_energies(),
                           self.settings.candidate_states(self.active_state))
        self._current_frame_in_detector = True

        event, n_candidates = self.select_event(trial_tracking)

        record = self._new_record(step, frame, trial_result, trial_tracking,
                                  event, n_candidates)
        if event is None:
            self.install_force(trial_result.active_gradient)
            self._record_force_state(record, trial_result, 'after_replay')
            self._finish_step(record, checkpoint_file)
            self._abort_frame = None
            return trial_result, trial_tracking

        # The hop is applied at the saved frame, so the detector's central
        # frame must be that same frame.  A mismatch means the gap window and
        # the rollback depth disagree, which would rescale momenta and switch
        # surfaces at the wrong geometry.
        assert_msg_critical(
            int(event.central_step) == int(frame.step),
            'LandauZenerSurfaceHoppingDynamics: the Landau-Zener event has '
            f'central frame {event.central_step} but the restorable frame is '
            f'{frame.step}; the gap window and the rollback depth disagree.')

        self.detector.mark_processed(event)

        # 5. stochastic decision
        xi = float(self.rng.random())
        decision_rng_state = copy.deepcopy(self.rng.bit_generator.state)
        record['random_number'] = xi
        record['hop_proposed'] = bool(xi < event.probability)

        if not record['hop_proposed']:
            # 6. no hop selected: retain the trial frame unchanged
            self._n_hops_rejected += 1
            record['hop_outcome'] = 'rejected'
            self.install_force(trial_result.active_gradient)
            self._record_force_state(record, trial_result, 'after_replay')
            self._finish_step(record, checkpoint_file)
            self._abort_frame = None
            return trial_result, trial_tracking

        # 7. restore frame n and attempt the hop there
        self.restore_frame(frame)
        self.rng.bit_generator.state = decision_rng_state
        self.current_step = step

        outcome = self.attempt_hop(event, frame)

        record.update({
            'hop_accepted': bool(outcome['accepted']),
            'hop_frustrated': bool(outcome['frustrated']),
            'kinetic_energy_before': outcome['kinetic_before'],
            'kinetic_energy_after': outcome['kinetic_after'],
            'hop_energy_residual': outcome['energy_residual'],
            'momentum_rescaling_scale': outcome['scaling_factor'],
            'momentum_rescaling_gamma': outcome['gamma'],
            'hop_message': outcome['message'],
            'target_raw_root_at_event': outcome['target_raw_state'],
            'target_derivative_state_index_at_event':
                outcome['target_derivative_state_index'],
        })

        if not outcome['accepted']:
            # Frustrated or rejected rescaling: the active state and the
            # velocities are unchanged, so simply replay the original step.
            if outcome['frustrated']:
                self._n_hops_frustrated += 1
                record['hop_outcome'] = 'frustrated'
            else:
                self._n_hops_rejected += 1
                record['hop_outcome'] = 'rescaling_rejected'
            self.install_force(result.active_gradient)
            self.md.simulation.step(1)
            self.current_step = step + 1

            replay_result, replay_tracking = self.evaluate_electronic_state(
                self.current_step)
            self.detector.push(
                self.current_step, self.active_state,
                replay_result.adiabatic_state_energies(),
                self.settings.candidate_states(self.active_state))
            self._current_frame_in_detector = True
            self.detector.mark_processed(event)
            self.install_force(replay_result.active_gradient)
            self._replace_trial_electronic_record(record, replay_result,
                                                  replay_tracking)
            self._record_force_state(record, replay_result, 'after_replay')
            self._finish_step(record, checkpoint_file)

            self._abort_frame = None
            return replay_result, replay_tracking

        # Accepted hop: install the target-state force at frame n and replay.
        self._n_hops_accepted += 1
        record['hop_outcome'] = 'accepted'

        record['electronic_energy_before'] = float(
            result.adiabatic_state_energies()[event.source_state])
        record['electronic_energy_after'] = float(
            result.adiabatic_state_energies()[event.target_state])

        target_raw = outcome['target_raw_state']
        # yapf: disable
        try:
            target_gradient = np.asarray(
                self.provider.compute_state_gradient(
                    frame.electronic_result.geometry, target_raw),
                dtype=float)
        except Exception as exc:
            raise SurfaceHoppingError(
                'failed to evaluate the accepted-hop target gradient at '
                f'step {step}: adiabatic raw root {target_raw}: {exc}')
        # yapf: enable

        hopped_result = copy.deepcopy(frame.electronic_result)
        hopped_result.active_gradient = target_gradient
        hopped_result.active_state = target_raw
        hopped_result.active_raw_state = target_raw
        hopped_result.active_tracked_state = (
            StateIndexConverter.raw_to_tracked(
                target_raw, hopped_result.state_permutation))
        hopped_result.derivative_state_index = hopped_result.selector_for(
            target_raw)
        hopped_result.gradient_converged = bool(
            np.all(np.isfinite(target_gradient)))

        if not hopped_result.gradient_converged:
            raise SurfaceHoppingError(
                f'nonfinite accepted-hop target gradient at step {step}.')

        self.install_force(hopped_result.active_gradient)
        self._record_force_state(record, hopped_result, 'before_replay')

        # Histories involving the old active state are invalid on the target
        # surface.  The restored tracker reference remains at frame n.
        self.detector.reset_pairs_involving(event.source_state)
        self.md.simulation.step(1)
        self.current_step = step + 1

        replay_result, replay_tracking = self.evaluate_electronic_state(
            self.current_step)

        self.detector.push(self.current_step, self.active_state,
                           replay_result.adiabatic_state_energies(),
                           self.settings.candidate_states(self.active_state))
        self._current_frame_in_detector = True

        self.install_force(replay_result.active_gradient)
        self._replace_trial_electronic_record(record, replay_result,
                                              replay_tracking)
        self._record_force_state(record, replay_result, 'after_replay')
        self._finish_step(record, checkpoint_file)

        self._abort_frame = None
        return replay_result, replay_tracking

    @staticmethod
    def _electronic_provenance(result, tracking):
        """
        Backend and descriptor provenance of one electronic frame.

        Every field needed to prove after the fact that the potentials were
        target-state totals, that the working reference was excluded, and
        that the persistent identity followed an electronic descriptor rather
        than energy order.
        """

        snapshot = getattr(result, 'snapshot', None)

        payload = {
            'tracking_descriptor': getattr(tracking, 'descriptor',
                                           'unspecified'),
            'tracking_ambiguity_ratio': float(
                getattr(tracking, 'ambiguity_ratio', 0.0)),
            'tracking_warnings': list(getattr(tracking, 'warnings', [])),
            'state_selectors': (None if result.state_selectors is None else
                                [int(s) for s in result.state_selectors]),
            'working_reference_energy': float(result.ground_energy),
            'state_tracking_reliable': bool(
                getattr(tracking, 'is_reliable', False)),
            'electronic_provenance_path': getattr(
                result, 'electronic_provenance_path', None),
        }

        if snapshot is None:
            payload.update({
                'electronic_backend': 'veloxchem',
                'electronic_method': 'linear_response',
                'backend_version': None,
                'calculation_id': None,
                'geometry_fingerprint': None,
                'settings_fingerprint': None,
                'energy_unit': 'hartree',
                'gradient_unit': 'hartree/bohr',
                'spin_metadata': None,
                'overlap_source': None,
                'selected_overlap_matrix': None,
                'native_overlap_matrix': None,
                'response_overlap_matrix': None,
                'native_spectral_norm': None,
                'selected_spectral_norm': None,
                'overlap_is_conditioned': None,
                'scf_reference_stability': None,
                'scf_reference_continuity': None,
                'reference_continuity_checked': False,
                'reference_continuous': None,
                'reference_continuity_reason': 'not_available',
                'reference_overlap_method': None,
                'docc_overlap_sigma_min': None,
                'docc_overlap_sigma_values': [],
                'somo_overlap_sigma_min': None,
                'somo_overlap_sigma_values': [],
                'electronic_warnings': [],
            })
            return payload

        def matrix_or_none(matrix):
            if matrix is None:
                return None
            return np.asarray(matrix, dtype=float).tolist()

        continuity = ({} if snapshot.reference_continuity is None else
                      dict(snapshot.reference_continuity))
        payload.update({
            'electronic_backend': snapshot.backend,
            'electronic_method': snapshot.method,
            'backend_version': snapshot.backend_version,
            'calculation_id': snapshot.calculation_id,
            'previous_calculation_id': snapshot.previous_calculation_id,
            'geometry_fingerprint': snapshot.geometry_fingerprint,
            'settings_fingerprint': snapshot.settings_fingerprint,
            'energy_unit': snapshot.energy_unit,
            'gradient_unit': snapshot.gradient_unit,
            'spin_metadata': snapshot.spin_metadata,
            'response_energies':
                np.asarray(snapshot.response_energies, dtype=float).tolist(),
            'overlap_source': snapshot.overlap_source,
            'selected_overlap_matrix':
                matrix_or_none(snapshot.overlap_to_previous),
            'native_overlap_matrix':
                matrix_or_none(snapshot.native_overlap_matrix),
            'response_overlap_matrix':
                matrix_or_none(snapshot.response_overlap_matrix),
            'native_spectral_norm': snapshot.native_spectral_norm,
            'selected_spectral_norm': snapshot.selected_spectral_norm,
            'overlap_is_conditioned': snapshot.overlap_is_conditioned,
            'scf_reference_stability': (
                None if snapshot.reference_stability is None else
                dict(snapshot.reference_stability)),
            'scf_reference_continuity': (
                None if snapshot.reference_continuity is None else
                dict(snapshot.reference_continuity)),
            'reference_continuity_checked': bool(
                continuity.get('reference_continuity_checked', False)),
            'reference_continuous': continuity.get(
                'reference_continuous', None),
            'reference_continuity_reason': continuity.get(
                'reference_continuity_reason', 'not_available'),
            'reference_overlap_method': continuity.get(
                'reference_overlap_method', None),
            'docc_overlap_sigma_min': continuity.get(
                'docc_overlap_sigma_min', None),
            'docc_overlap_sigma_values': continuity.get(
                'docc_overlap_sigma_values', []),
            'somo_overlap_sigma_min': continuity.get(
                'somo_overlap_sigma_min', None),
            'somo_overlap_sigma_values': continuity.get(
                'somo_overlap_sigma_values', []),
            'electronic_warnings': list(snapshot.warnings),
        })

        return payload

    @staticmethod
    def _record_force_state(record, result, suffix):
        """
        Records which adiabatic root owns the gradient installed in OpenMM.

        ``suffix`` is ``before_replay`` for the central-frame force and
        ``after_replay`` for the retained endpoint force.
        """

        record[f'force_tracked_state_{suffix}'] = (
            None if result.active_tracked_state is None else int(
                result.active_tracked_state))
        record[f'force_raw_state_{suffix}'] = int(result.active_raw_state)
        record[f'force_derivative_state_index_{suffix}'] = (
            None if result.derivative_state_index is None else int(
                result.derivative_state_index))

    def _replace_trial_electronic_record(self, record, result, tracking):
        """
        Replaces fields belonging to a discarded trial endpoint with fields
        from the endpoint retained after replay.
        """

        tracked_energies = result.tracked_state_energies()
        raw_to_tracked = result.raw_to_tracked_permutation()

        record.update({
            'raw_state_energies': np.asarray(result.state_energies,
                                             dtype=float).tolist(),
            'tracked_state_energies': tracked_energies.tolist(),
            'state_permutation': np.asarray(tracking.permutation,
                                            dtype=int).tolist(),
            'tracked_to_raw': np.asarray(tracking.permutation,
                                         dtype=int).tolist(),
            'raw_to_tracked': np.asarray(raw_to_tracked, dtype=int).tolist(),
            'tracking_scores': np.asarray(tracking.scores,
                                          dtype=float).tolist(),
            'tracking_similarity_matrix': np.asarray(tracking.similarity_matrix,
                                                     dtype=float).tolist(),
            'tracking_confidence': float(tracking.confidence),
            'tracking_method': tracking.method,
            'active_tracked_state': int(result.active_tracked_state),
            'active_raw_root': int(result.active_raw_state),
            'derivative_state_index':
                (None if result.derivative_state_index is None else int(
                    result.derivative_state_index)),
            'scf_converged': bool(result.scf_converged),
            'response_converged': bool(result.response_converged),
            'gradient_converged': bool(result.gradient_converged),
            'n_scf_calls': int(self.provider.n_scf_calls),
            'n_response_calls': int(self.provider.n_response_calls),
            'n_gradient_calls': int(self.provider.n_gradient_calls),
        })
        record.update(self._electronic_provenance(result, tracking))
        self._update_endpoint_observables(record, result)

    def _kinetic_energy_au(self, all_velocities_nm_ps):
        """Returns the rescaling atom set's kinetic energy in Hartree."""

        indices = self.rescaling_atoms
        masses_au = self.get_masses_in_au(indices)
        velocities_au = (np.asarray(all_velocities_nm_ps, dtype=float)[indices]
                         * NM_PER_PS_IN_AU_VELOCITY)

        return float(0.5 * np.sum(
            masses_au * np.sum(velocities_au**2, axis=1)))

    def _temperature_kelvin(self, kinetic_energy_au):
        """Converts the reported atom set's kinetic energy to temperature."""

        n_atoms = len(self.rescaling_atoms)
        n_dof = 3 * n_atoms - 6
        if n_atoms <= 2:
            n_dof = 3 * n_atoms - 5
        if n_dof <= 0:
            return 0.0

        ke_joule = (float(kinetic_energy_au) * hartree_in_kjpermol() *
                    1000.0 / 6.02214076e23)
        k_B = 1.380649e-23

        return float(2.0 * ke_joule / (n_dof * k_B))

    def _update_endpoint_observables(self, record, result):
        """Synchronizes all reported energy terms at the retained endpoint.

        Electronic energies in a step record belong to the propagated or
        replayed endpoint.  Kinetic energy must be read from that same OpenMM
        state; combining it with the saved central-frame velocity produces a
        quantity that is not the total energy at either time, especially
        across an accepted hop.
        """

        state = self.context.getState(getVelocities=True)
        velocity_unit = unit.nanometer / unit.picosecond
        velocities_nm_ps = np.asarray(
            state.getVelocities(asNumpy=True).value_in_unit(velocity_unit),
            dtype=float)
        kinetic_energy = self._kinetic_energy_au(velocities_nm_ps)
        adiabatic_energies = result.adiabatic_state_energies()
        active_energy = float(adiabatic_energies[self.active_state])

        record.update({
            'endpoint_step': int(self.current_step),
            'endpoint_time': float(
                state.getTime().value_in_unit(unit.femtosecond)),
            'endpoint_active_state': int(self.active_state),
            'energy_observation_frame': 'retained_endpoint',
            'kinetic_energy': kinetic_energy,
            'temperature_kelvin': self._temperature_kelvin(kinetic_energy),
            'active_state_energy': active_energy,
            'total_energy': kinetic_energy + active_energy,
            'energy_gaps': [
                abs(float(adiabatic_energies[i + 1] -
                          adiabatic_energies[i]))
                for i in range(len(adiabatic_energies) - 1)
            ],
        })
        record['minimum_energy_gap'] = (
            float(min(record['energy_gaps']))
            if record['energy_gaps'] else None)

    def _new_record(self, step, frame, trial_result, trial_tracking, event,
                    n_candidates):
        """
        Builds the machine-readable diagnostics record for one step.

        :param step:
            The central frame index.
        :param frame:
            The saved :class:`TrajectoryFrame` at the central frame.
        :param trial_result:
            The :class:`ElectronicStateResult` at the trial frame.
        :param trial_tracking:
            The :class:`TrackingResult` at the trial frame.
        :param event:
            The selected :class:`LandauZenerEvent`, or None.
        :param n_candidates:
            Number of valid candidate events at this step.

        :return:
            A JSON-serializable dictionary.
        """
        # Without a saved frame this is the pre-trajectory snapshot: the active
        # state is the current one and the energies are those just evaluated.
        time_fs = 0.0
        active_state = int(self.active_state)
        if frame is not None:

            time_fs = frame.time_fs
            active_state = frame.active_state

        # yapf: disable
        record = {
            'trajectory_id': self.settings.trajectory_id,
            'step': int(step),
            'time': float(time_fs),
            'central_step': (None if frame is None else int(frame.step)),
            'central_time': (None if frame is None else float(frame.time_fs)),
            'active_state': int(active_state),
            'raw_state_energies':
                np.asarray(trial_result.state_energies).tolist(),
            'tracked_state_energies':
                trial_result.tracked_state_energies().tolist(),
            'state_permutation':
                np.asarray(trial_tracking.permutation).tolist(),
            'tracked_to_raw': np.asarray(trial_tracking.permutation).tolist(),
            'raw_to_tracked':
                trial_result.raw_to_tracked_permutation().tolist(),
            'tracking_scores': np.asarray(trial_tracking.scores).tolist(),
            'tracking_similarity_matrix': np.asarray(
                trial_tracking.similarity_matrix).tolist(),
            'tracking_confidence': float(trial_tracking.confidence),
            'tracking_method': trial_tracking.method,
            'active_tracked_state': int(trial_result.active_tracked_state),
            'active_raw_root': int(trial_result.active_raw_state),
            'derivative_state_index':
                (None if trial_result.derivative_state_index is None else int(
                    trial_result.derivative_state_index)),
            'candidate_source_state': None,
            'candidate_target_state': None,
            'event_tracked_to_raw': None,
            'event_raw_to_tracked': None,
            'source_raw_root_at_event': None,
            'source_derivative_state_index_at_event': None,
            'target_raw_root_at_event': None,
            'target_derivative_state_index_at_event': None,
            'minimum_gap': None,
            'gap_curvature': None,
            'fitted_event_time': None,
            'fit_quality': None,
            'landau_zener_exponent': None,
            'landau_zener_probability': None,
            'random_number': None,
            'hop_proposed': False,
            'hop_accepted': False,
            'hop_frustrated': False,
            'hop_outcome': 'none',
            'n_candidate_events': int(n_candidates),
            'multistate_warning': bool(self.detector.last_multistate_warning),
            'momentum_rescaling_mode': self.settings.momentum_rescaling,
            'momentum_rescaling_scale': None,
            'momentum_rescaling_gamma': None,
            'kinetic_energy_before': None,
            'kinetic_energy_after': None,
            'electronic_energy_before': None,
            'electronic_energy_after': None,
            'hop_energy_residual': None,
            'hop_message': '',
            'force_tracked_state_before_replay': None,
            'force_raw_state_before_replay': None,
            'force_derivative_state_index_before_replay': None,
            'force_tracked_state_after_replay': None,
            'force_raw_state_after_replay': None,
            'force_derivative_state_index_after_replay': None,
            'scf_converged': bool(trial_result.scf_converged),
            'response_converged': bool(trial_result.response_converged),
            'gradient_converged': bool(trial_result.gradient_converged),
            'n_scf_calls': int(self.provider.n_scf_calls),
            'n_response_calls': int(self.provider.n_response_calls),
            'n_gradient_calls': int(self.provider.n_gradient_calls),
        }
        # yapf: enable

        # Preserve the central kinetic energy as explicit event provenance.
        # The legacy/public kinetic_energy field is populated below from the
        # retained endpoint so it is synchronized with the electronic terms.
        if frame is not None:
            record['central_kinetic_energy'] = self._kinetic_energy_au(
                frame.velocities_nm_ps)
        else:
            record['central_kinetic_energy'] = None

        self._update_endpoint_observables(record, trial_result)

        if event is not None:
            event_permutation = np.asarray(
                frame.electronic_result.state_permutation, dtype=int)
            event_inverse = (
                frame.electronic_result.raw_to_tracked_permutation())
            source_raw = int(event.source_state)
            target_raw = int(event.target_state)
            record.update({
                'candidate_source_state': int(event.source_state),
                'candidate_target_state': int(event.target_state),
                'event_tracked_to_raw': event_permutation.tolist(),
                'event_raw_to_tracked': np.asarray(event_inverse,
                                                   dtype=int).tolist(),
                'source_raw_root_at_event': int(source_raw),
                'source_derivative_state_index_at_event':
                    frame.electronic_result.selector_for(source_raw),
                'target_raw_root_at_event': int(target_raw),
                'target_derivative_state_index_at_event':
                    frame.electronic_result.selector_for(target_raw),
                'minimum_gap': float(event.minimum_gap),
                'gap_curvature': float(event.gap_curvature),
                'fitted_event_time': float(event.event_time * AU_IN_FS),
                'fit_quality': float(event.fit_quality),
                'landau_zener_exponent': float(event.exponent),
                'landau_zener_probability': float(event.probability),
            })

        if frame is not None:
            self._record_force_state(record, frame.electronic_result,
                                     'before_replay')

        record.update(self._electronic_provenance(trial_result,
                                                  trial_tracking))

        return record

    def _finish_step(self, record, checkpoint_file):
        """
        Appends a diagnostics record and writes a checkpoint when due.

        :param record:
            The diagnostics dictionary of the completed step.
        :param checkpoint_file:
            Path of the checkpoint file.
        """

        record['active_state_after'] = int(self.active_state)

        self.diagnostics.append(record)
        self._write_diagnostics_record(record)

        stability_record = record.get('scf_reference_stability', None)
        if (stability_record is not None and
                stability_record.get('reason_for_check') == 'periodic'):
            self._scf_stability_last_accepted_step = int(self.current_step)

        interval = int(self.settings.checkpoint_interval)

        if interval and self.current_step % interval == 0:
            self._write_checkpoint(checkpoint_file, best_effort=True)

    def _format_step_output(self, record):
        """
        Formats a comprehensive human-readable output for a dynamics step.

        :param record:
            The diagnostics record dictionary from _new_record().

        :return:
            Formatted string with step information.
        """
        step = record.get('step', self.current_step)
        time_fs = record.get('time', 0.0)
        time_ps = time_fs / 1000.0

        active_before = record.get('active_state', self.active_state)
        active_after = record.get('active_state_after', active_before)

        lines = [
            f"\n{'='*70}",
            f"SURFACE HOPPING STEP {step:6d}  Time: {time_ps:8.2f} ps",
            f"{'='*70}", f"Active State: {active_before} -> {active_after}"
        ]

        # Energy information
        if 'active_state_energy' in record:
            active_energy = record['active_state_energy']
            lines.append(f"Active State Energy: {active_energy:.8f} Eh")

        if 'kinetic_energy' in record:
            ke = record['kinetic_energy']
            lines.append(f"Kinetic Energy:      {ke:.8f} Eh")

        if 'total_energy' in record:
            total = record['total_energy']
            lines.append(f"Total Energy:        {total:.8f} Eh")

        if 'temperature_kelvin' in record:
            temp = record['temperature_kelvin']
            lines.append(f"Temperature:          {temp:.2f} K")

        # State energies and gaps
        if 'tracked_state_energies' in record:
            energies = record['tracked_state_energies']
            lines.append("\nState Energies (Eh):")
            for i, energy in enumerate(energies):
                marker = " *" if i == active_after else ""
                lines.append(f"  State {i:2d}: {energy:12.8f}{marker}")

        if 'energy_gaps' in record:
            gaps = record['energy_gaps']
            min_gap = record.get('minimum_energy_gap')
            lines.append("\nEnergy Gaps (Eh):")
            for i, gap in enumerate(gaps):
                lines.append(f"  Gap {i}->{i+1}: {gap:.8f}")
            if min_gap is not None:
                lines.append(f"  Minimum Gap:  {min_gap:.8f}")

        # Hopping information
        if record.get('hop_proposed', False):
            lines.append("\nHOP PROPOSED:")
            lines.append(
                f"  Source State:      {record.get('candidate_source_state', 'N/A')}"
            )
            lines.append(
                f"  Target State:      {record.get('candidate_target_state', 'N/A')}"
            )
            lines.append(
                f"  Minimum Gap:       {record.get('minimum_gap', 'N/A'):.8f} Eh"
            )
            lines.append(
                f"  Gap Curvature:     {record.get('gap_curvature', 'N/A'):.8f} Eh/a.u.²"
            )
            lines.append(
                f"  Event Time:        {record.get('fitted_event_time', 'N/A'):.2f} fs"
            )
            lines.append(
                f"  LZ Exponent:       {record.get('landau_zener_exponent', 'N/A'):.4f}"
            )
            lines.append(
                f"  LZ Probability:    {record.get('landau_zener_probability', 'N/A'):.6f}"
            )
            lines.append(
                f"  Random Number:     {record.get('random_number', 'N/A'):.6f}"
            )

            if record.get('hop_accepted', False):
                lines.append("  -> HOP ACCEPTED")
                if 'hop_energy_residual' in record:
                    lines.append(
                        f"  Energy Residual:   {record['hop_energy_residual']:.8f} Eh"
                    )
                if 'momentum_rescaling_scale' in record:
                    lines.append(
                        f"  Momentum Scale:    {record['momentum_rescaling_scale']:.6f}"
                    )
            elif record.get('hop_frustrated', False):
                lines.append(
                    f"  -> HOP FRUSTRATED: {record.get('hop_message', 'N/A')}")
            else:
                lines.append("  -> HOP REJECTED")

        # SCF convergence
        lines.append(
            f"\nSCF Converged:      {record.get('scf_converged', True)}")
        lines.append(
            f"Response Converged: {record.get('response_converged', True)}")
        lines.append(
            f"Gradient Converged: {record.get('gradient_converged', True)}")

        # Performance counters
        lines.append("\nCounters:")
        lines.append(f"  SCF Calls:         {record.get('n_scf_calls', 0)}")
        lines.append(
            f"  Response Calls:    {record.get('n_response_calls', 0)}")
        lines.append(
            f"  Gradient Calls:    {record.get('n_gradient_calls', 0)}")

        lines.append(f"{'='*70}\n")

        return "\n".join(lines)

    def print_step_summary(self, record):
        """
        Prints a formatted summary of the current dynamics step.

        :param record:
            The diagnostics record dictionary.
        """
        output = self._format_step_output(record)
        self.ostream.print_info(output)
        self.ostream.flush()

    def _write_pdb_frame(self, filename, step=None, state=None):
        """
        Writes the current system geometry to a PDB file.

        :param filename:
            Path to the PDB file. If contains '{step}' or '{state}', these will be
            replaced with the current step and active state.
        :param step:
            Current step number (defaults to self.current_step).
        :param state:
            Current active state (defaults to self.active_state).
        """
        try:
            import openmm.app as app
        except ImportError:
            self.ostream.print_warning('PDB writing requires OpenMM.app module')
            return

        if step is None:
            step = self.current_step
        if state is None:
            state = self.active_state

        # Format filename with step and state if placeholders exist
        formatted_filename = filename
        if '{step}' in filename:
            formatted_filename = filename.format(step=step)
        if '{state}' in filename:
            formatted_filename = formatted_filename.format(state=state)

        state = self.context.getState(getPositions=True)
        positions = state.getPositions()

        try:
            app.PDBFile.writeFile(self.md.topology, positions,
                                  open(formatted_filename, 'w'))
        except Exception as exc:
            self.ostream.print_warning(
                f'Failed to write PDB frame to {formatted_filename}: {exc}')

    def _check_ensemble(self):
        """
        Warns about thermostats and barostats, which have no scientifically
        defensible chronology relative to the hopping algorithm implemented
        here.
        """

        ensemble = str(getattr(self.md, 'ensemble', 'NVE')).upper()

        if ensemble != 'NVE':
            self.ostream.print_warning(
                f'Surface hopping is running in the {ensemble} ensemble. '
                'Momentum rescaling and thermostatting compete for the same '
                'kinetic energy, and no defensible chronology between them '
                'is implemented. NVE is strongly recommended for production '
                'trajectories.')
            self.ostream.flush()

    def get_openmm_timestep_in_fs(self):
        """
        Reads the nuclear timestep actually used by the OpenMM integrator.

        :return:
            The integrator step size in femtoseconds, or None when the
            underlying dynamics object does not expose one.
        """

        integrator = getattr(getattr(self.md, 'simulation', None), 'integrator',
                             None)

        for quantity in (getattr(integrator, 'getStepSize', lambda: None)(),
                         getattr(self.md, 'timestep', None)):
            if quantity is None:
                continue
            try:
                return float(quantity.value_in_unit(unit.femtosecond))
            except Exception:
                continue

        return None

    def _check_nuclear_timestep(self):
        """
        Verifies that the declared electronic timestep is the timestep the
        nuclei are actually propagated with.

        The gap curvature is a finite difference over consecutive nuclear
        frames, ``Zddot = (Z_{n-1} - 2 Z_n + Z_{n+1}) / h**2``, and ``h`` is
        taken from :attr:`SurfaceHoppingSettings.electronic_timestep`.  If that
        setting disagrees with the OpenMM integrator step size the curvature is
        wrong by the square of the ratio and every Landau-Zener probability is
        silently wrong, so this is refused rather than warned about.
        """

        openmm_timestep = self.get_openmm_timestep_in_fs()

        if openmm_timestep is None:
            self.ostream.print_warning(
                'Surface hopping could not read the OpenMM integrator '
                'timestep; the declared electronic_timestep of '
                f'{float(self.settings.electronic_timestep):.6g} fs is assumed '
                'to be the nuclear propagation timestep.')
            self.ostream.flush()
            return

        declared = float(self.settings.electronic_timestep)

        assert_msg_critical(
            abs(openmm_timestep - declared) <= 1.0e-8 * max(1.0, abs(declared)),
            'LandauZenerSurfaceHoppingDynamics: electronic_timestep is '
            f'{declared:.6g} fs but the OpenMM integrator propagates '
            f'{openmm_timestep:.6g} fs per step. The gap curvature is a finite '
            'difference over consecutive nuclear frames, so the two must be '
            'equal; otherwise every Landau-Zener probability is wrong by the '
            'square of the ratio.')

    # ------------------------------------------------------------------
    # diagnostics and restart
    # ------------------------------------------------------------------

    def _open_diagnostics(self, filename):
        """
        Opens the machine-readable diagnostics file for retained records.

        A checkpoint records the flushed byte offset of its diagnostics file.
        When that exact file is reused for a restart, any records written after
        the checkpoint are discarded before replay appends their replacements.
        This prevents stale trial history and duplicate stochastic decisions
        from surviving a checkpoint rollback.

        A trajectory that does not continue a checkpoint starts its diagnostics
        file from scratch.  Appending to records of an unrelated earlier run
        would produce a file with repeated step indices and repeated random
        numbers that no longer describes one trajectory.  Successive
        :func:`run` calls on the same controller do continue one trajectory and
        keep appending.

        :param filename:
            Path of the JSON-lines diagnostics file.
        """

        if not filename:
            return

        path = Path(filename)
        restart_path = self._restart_diagnostics_path
        restart_offset = self._restart_diagnostics_offset
        continues_trajectory = (self._diagnostics_started or
                                restart_path is not None)

        if restart_path is not None and restart_offset is not None:
            same_file = path.resolve() == Path(restart_path).resolve()

            if same_file:
                offset = int(restart_offset)

                if offset < 0:
                    raise SurfaceHoppingError(
                        'checkpoint contains a negative diagnostics offset.')

                if not path.exists():
                    if offset != 0:
                        raise SurfaceHoppingError(
                            'the checkpoint diagnostics file is missing; '
                            'use a different diagnostics_file for the '
                            'continuation.')
                else:
                    current_size = path.stat().st_size

                    if current_size < offset:
                        raise SurfaceHoppingError(
                            'the diagnostics file is shorter than the offset '
                            'stored in the checkpoint; refusing to append.')

                    if current_size > offset:
                        with open(path, 'r+b') as handle:
                            handle.truncate(offset)

        # Consume the pending rollback metadata even when the caller selected
        # a different diagnostics path intentionally.
        self._restart_diagnostics_path = None
        self._restart_diagnostics_offset = None

        self._diagnostics_path = str(path.resolve())
        self._diagnostics_handle = open(path,
                                        'a' if continues_trajectory else 'w',
                                        encoding='utf-8')
        self._diagnostics_started = True

    def _write_diagnostics_record(self, record):
        """
        Writes one JSON-lines diagnostics record.

        :param record:
            The diagnostics dictionary.
        """

        if self._diagnostics_handle is None:
            return

        self._diagnostics_handle.write(json.dumps(record) + '\n')
        self._diagnostics_handle.flush()

    def _close_diagnostics(self):
        """
        Closes the diagnostics file.
        """

        if self._diagnostics_handle is not None:
            self._diagnostics_handle.close()
            self._diagnostics_handle = None

    def _write_checkpoint(self, filename, best_effort=False):
        """
        Writes a restart checkpoint.

        :param filename:
            Path of the checkpoint file.
        :param best_effort:
            When True, checkpoint failures are reported as warnings instead of
            aborting the trajectory.
        """

        if not filename:
            return

        # yapf: disable
        try:
            state = self.context.getState(getPositions=True, getVelocities=True)

            diagnostics_path = self._diagnostics_path
            diagnostics_offset = None

            if diagnostics_path is not None:
                if self._diagnostics_handle is not None:
                    self._diagnostics_handle.flush()

                path = Path(diagnostics_path)
                if path.is_file():
                    diagnostics_offset = int(path.stat().st_size)

            try:
                box = np.asarray(
                    state.getPeriodicBoxVectors(
                        asNumpy=True).value_in_unit(unit.nanometer),
                    dtype=float)
            except Exception:
                box = None

            payload = {
                'step': int(self.current_step),
                'time_fs':
                    float(state.getTime().value_in_unit(unit.femtosecond)),
                'positions_nm':
                    np.asarray(
                        state.getPositions(
                            asNumpy=True).value_in_unit(unit.nanometer),
                        dtype=float),
                'velocities_nm_ps': np.asarray(
                    state.getVelocities(asNumpy=True).value_in_unit(
                        unit.nanometer / unit.picosecond),
                    dtype=float),
                'box_vectors': box,
                'active_state': int(self.active_state),
                'active_state_index_space': self._ACTIVE_STATE_INDEX_SPACE,
                'rng_state': self.rng.bit_generator.state,
                'detector_history': self.detector.get_history_snapshot(),
                'tracker_reference': self.tracker._reference,
                # Native backend handles are not serializable, so only the
                # portable description of the accepted electronic reference is
                # written.  The provider recomputes it at restart, which keeps
                # the descriptor chain unbroken instead of resuming with an
                # empty tracking history.
                'provider_reference':
                    self.provider.reference_restart_payload()
                    if hasattr(self.provider, 'reference_restart_payload')
                    else None,
                'state_permutation': np.asarray(self._state_permutation,
                                                dtype=int),
                'tracked_spin_multiplicities':
                    None if self._tracked_spin_multiplicities is None else
                    np.asarray(self._tracked_spin_multiplicities, dtype=int),
                'current_frame_in_detector': bool(
                    self._current_frame_in_detector),
                'diagnostics_path': diagnostics_path,
                'diagnostics_offset': diagnostics_offset,
                'settings_fingerprint': self.settings.fingerprint(),
                'settings': self.settings.as_dict(),
                'n_hops_accepted': self._n_hops_accepted,
                'n_hops_frustrated': self._n_hops_frustrated,
                'n_hops_rejected': self._n_hops_rejected,
                'scf_stability_check_count':
                    int(self._scf_stability_check_count),
                'scf_stability_initial_completed':
                    bool(self._scf_stability_initial_completed),
                'scf_stability_last_accepted_step':
                    self._scf_stability_last_accepted_step,
                'electronic_provenance_calculation_index':
                    int(self._electronic_provenance_calculation_index),
            }

            with open(filename, 'wb') as fh:
                pickle.dump(payload, fh)

        except Exception as exc:
            if not best_effort:
                raise
            self.ostream.print_warning(
                f'Failed to write surface-hopping checkpoint: {exc}')
            self.ostream.flush()
        # yapf: enable

    def restart_from_checkpoint(self, filename):
        """
        Restores a trajectory from a checkpoint.

        The checkpoint is refused when the surface-hopping settings that
        affect the trajectory have changed, since the restarted trajectory
        would then not continue the original one.

        :param filename:
            Path of the checkpoint file.
        """

        path = Path(filename)

        assert_msg_critical(
            path.is_file(),
            f'LandauZenerSurfaceHoppingDynamics: checkpoint {filename} does '
            'not exist.')

        with open(path, 'rb') as fh:
            payload = pickle.load(fh)

        settings_compatible = (
            payload.get('settings_fingerprint') == self.settings.fingerprint())

        # Checkpoints written before all trajectory-changing settings were
        # included in the fingerprint carry a legacy digest.  They also carry
        # the full settings dictionary, which can be re-fingerprinted using
        # the corrected key set without weakening compatibility checks.
        if not settings_compatible and isinstance(payload.get('settings'),
                                                  dict):
            try:
                checkpoint_settings = SurfaceHoppingSettings.from_dict(
                    payload['settings'])
                settings_compatible = (checkpoint_settings.fingerprint() ==
                                       self.settings.fingerprint())
            except Exception:
                settings_compatible = False

        assert_msg_critical(
            settings_compatible,
            'LandauZenerSurfaceHoppingDynamics: refusing to restart from a '
            'checkpoint written with incompatible surface-hopping settings.')

        checkpoint_index_space = payload.get('active_state_index_space')
        checkpoint_permutation = np.asarray(
            payload.get('state_permutation',
                        np.arange(self.settings.number_of_states)),
            dtype=int)
        if checkpoint_index_space is None:
            # Legacy checkpoints stored active_state in persistent-character
            # space. The two conventions are equivalent only while the
            # tracked-to-raw mapping is the identity.
            assert_msg_critical(
                np.array_equal(
                    checkpoint_permutation,
                    np.arange(self.settings.number_of_states)),
                'LandauZenerSurfaceHoppingDynamics: refusing a legacy '
                'checkpoint with a non-identity state permutation because '
                'its active_state is in the obsolete character-label index '
                'space. Restart from an earlier identity-mapped checkpoint.')
        else:
            assert_msg_critical(
                checkpoint_index_space == self._ACTIVE_STATE_INDEX_SPACE,
                'LandauZenerSurfaceHoppingDynamics: checkpoint active-state '
                f'index space {checkpoint_index_space!r} is incompatible '
                f'with {self._ACTIVE_STATE_INDEX_SPACE!r}.')

        if payload.get('box_vectors') is not None:
            self.context.setPeriodicBoxVectors(*[
                mm.Vec3(*row) * unit.nanometer for row in payload['box_vectors']
            ])

        self.context.setPositions(payload['positions_nm'] * unit.nanometer)
        self.context.setVelocities(payload['velocities_nm_ps'] *
                                   unit.nanometer / unit.picosecond)
        self.context.setTime(payload['time_fs'] * unit.femtosecond)

        self.current_step = int(payload['step'])
        self.active_state = int(payload['active_state'])
        self.rng.bit_generator.state = payload['rng_state']
        self.detector.set_history_snapshot(payload['detector_history'])
        self.tracker._reference = payload['tracker_reference']

        if hasattr(self.provider, 'restart_from'):
            self.provider.restart_from(payload.get('provider_reference'))

        # yapf: disable
        self._state_permutation = checkpoint_permutation
        # yapf: enable
        stored_spin = payload.get('tracked_spin_multiplicities')
        self._tracked_spin_multiplicities = (
            None if stored_spin is None else
            np.asarray(stored_spin, dtype=int))
        self._current_frame_in_detector = bool(
            payload.get('current_frame_in_detector', True))
        self._restart_diagnostics_path = payload.get('diagnostics_path')
        self._restart_diagnostics_offset = payload.get('diagnostics_offset')

        self._n_hops_accepted = int(payload.get('n_hops_accepted', 0))
        self._n_hops_frustrated = int(payload.get('n_hops_frustrated', 0))
        self._n_hops_rejected = int(payload.get('n_hops_rejected', 0))
        self._scf_stability_check_count = int(
            payload.get('scf_stability_check_count', 0))
        self._scf_stability_initial_completed = bool(payload.get(
            'scf_stability_initial_completed',
            self._scf_stability_mode == 'off' or self.current_step > 0))
        stored_stability_step = payload.get(
            'scf_stability_last_accepted_step', None)
        self._scf_stability_last_accepted_step = (
            None if stored_stability_step is None else
            int(stored_stability_step))
        self._electronic_provenance_calculation_index = int(payload.get(
            'electronic_provenance_calculation_index', 0))

    def get_summary(self):
        """
        :return:
            A dictionary summarizing the hopping statistics of the
            trajectory.
        """

        return {
            'trajectory_id': self.settings.trajectory_id,
            'steps': int(self.current_step),
            'final_active_state': int(self.active_state),
            'hops_accepted': self._n_hops_accepted,
            'hops_frustrated': self._n_hops_frustrated,
            'hops_rejected': self._n_hops_rejected,
            'scf_calls': int(self.provider.n_scf_calls),
            'response_calls': int(self.provider.n_response_calls),
            'gradient_calls': int(self.provider.n_gradient_calls),
            'scf_stability_checks': int(self._scf_stability_check_count),
        }
