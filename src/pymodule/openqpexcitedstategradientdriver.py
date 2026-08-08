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
from .openqpstatetracker import OpenQPMRSFStateTracker


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
        self.target_state_index = 3
        self.tddft_type = 'mrsf'
        self.nstates = 10

        self.excited_state_energy = None
        self.total_energy = None

        # ``excited_state_energy`` is a *scalar* excitation energy of the
        # selected root relative to the lowest MRSF target.  It can never be
        # used to build a multi-state ladder.  The two attributes below carry
        # what a multi-state consumer actually needs and are refreshed by
        # every evaluation: the working-reference energy (diagnostics only)
        # and the total target-state energies in raw target order.
        self.reference_energy = None
        self.target_state_energies = None

        # State tracking is opt-in so that existing fixed-root calculations
        # retain their exact control flow.  A tracked calculation provides its
        # energy and gradient atomically.
        self.state_tracking = False
        self.computes_energy_with_gradient = False
        self.state_tracker = None
        self.last_tracking_result = None
        self.tracking_history = []
        self._tracked_geometry_signature = None
        self._input_keywords['gradient'].update({
            'state_tracking':
                ('bool', 'enable native OpenQP MRSF state tracking'),
            'minimum_tracking_overlap':
                ('float', 'minimum accepted target-state overlap'),
            'maximum_ambiguity_ratio':
                ('float', 'maximum accepted second/best overlap ratio'),
            'tracking_failure':
                ('str_lower', 'tracking failure policy: error or warn'),
        })

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
        self.target_state_index = state
        if self.state_tracker is not None:
            self.state_tracker = self._new_state_tracker(
                target_state=state)
            self.last_tracking_result = None
            self.tracking_history = []
            self._tracked_geometry_signature = None

    def set_tddft_type(self, tddft_type):
        """
        Sets the OpenQP TD type (`mrsf`, `sf`, `tda`, or `rpa`).
        """

        ttype = str(tddft_type).strip().lower()
        assert_msg_critical(
            ttype in ('mrsf', 'sf', 'tda', 'rpa'),
            f'OpenQPExcitedStateGradientDriver: invalid tddft_type '
            f'"{tddft_type}"')
        assert_msg_critical(
            not self.state_tracking or ttype == 'mrsf',
            'OpenQPExcitedStateGradientDriver: native state tracking is '
            'available only for MRSF.')
        self.tddft_type = ttype
        self.exc_driver.set_tddft_type(ttype)

    def set_nstates(self, nstates):
        """
        Sets the number of excited-state roots to compute.
        """

        new_nstates = int(nstates)
        assert_msg_critical(
            new_nstates > 0,
            'OpenQPExcitedStateGradientDriver: nstates must be positive.')
        assert_msg_critical(
            not self.state_tracking or
            new_nstates >= int(self.target_state_index),
            'OpenQPExcitedStateGradientDriver: nstates must include the '
            'persistent target state.')
        changed = new_nstates != self.nstates
        self.nstates = new_nstates
        self.exc_driver.set_nstates(nstates)
        if changed and self.state_tracker is not None:
            self.state_tracker = self._new_state_tracker()
            self.last_tracking_result = None
            self.tracking_history = []
            self._tracked_geometry_signature = None

    def enable_state_tracking(self,
                              minimum_overlap=0.5,
                              maximum_ambiguity_ratio=0.8,
                              failure_mode='error',
                              guard_ground_state=True,
                              max_stale_reference_steps=10):
        """Enables native OpenQP MRSF persistent-state tracking.

        The currently configured ``state_deriv_index`` establishes the
        persistent identity at the first geometry.  Subsequent
        ``state_deriv_index`` values report the raw OpenQP root carrying that
        identity.

        ``guard_ground_state`` keeps an excited identity off the lowest raw
        root; see :class:`OpenQPMRSFStateTracker`.  Leave it enabled for every
        excited-state optimization -- without it a torsion scan that reaches
        an S1/S0 conical intersection silently continues on S0.
        """

        assert_msg_critical(
            self.tddft_type == 'mrsf',
            'OpenQPExcitedStateGradientDriver: native state tracking is '
            'available only for MRSF.')
        assert_msg_critical(
            int(self.nstates) >= int(self.target_state_index),
            'OpenQPExcitedStateGradientDriver: nstates must include the '
            'persistent target state.')

        self.state_tracker = OpenQPMRSFStateTracker(
            nstates=self.nstates,
            target_state=self.target_state_index,
            minimum_overlap=minimum_overlap,
            maximum_ambiguity_ratio=maximum_ambiguity_ratio,
            failure_mode=failure_mode,
            guard_ground_state=guard_ground_state,
            max_stale_reference_steps=max_stale_reference_steps)
        self.state_tracking = True
        self.computes_energy_with_gradient = True
        self.last_tracking_result = None
        self.tracking_history = []
        self._tracked_geometry_signature = None

    def disable_state_tracking(self):
        """Disables tracking and restores the configured persistent root."""

        self.state_tracking = False
        self.computes_energy_with_gradient = False
        self.state_tracker = None
        self.last_tracking_result = None
        self.tracking_history = []
        self._tracked_geometry_signature = None
        self.state_deriv_index = self.target_state_index

    def set_state_tracker(self, tracker):
        """Installs a configured :class:`OpenQPMRSFStateTracker`."""

        assert_msg_critical(
            isinstance(tracker, OpenQPMRSFStateTracker),
            'OpenQPExcitedStateGradientDriver: invalid OpenQP state tracker.')
        assert_msg_critical(
            tracker.nstates == self.nstates,
            'OpenQPExcitedStateGradientDriver: tracker nstates mismatch.')
        assert_msg_critical(
            tracker.target_state == self.target_state_index,
            'OpenQPExcitedStateGradientDriver: tracker target-state mismatch.')
        assert_msg_critical(
            self.tddft_type == 'mrsf',
            'OpenQPExcitedStateGradientDriver: native state tracking is '
            'available only for MRSF.')

        self.state_tracker = tracker
        self.state_tracking = True
        self.computes_energy_with_gradient = True
        self.last_tracking_result = None
        self.tracking_history = list(tracker.history)
        self._tracked_geometry_signature = None
        if tracker.has_reference:
            self.state_deriv_index = tracker.selected_raw_root

    def save_tracking_state(self, filename):
        """Writes the complete OpenQP tracking restart state to JSON."""

        assert_msg_critical(
            self.state_tracker is not None,
            'OpenQPExcitedStateGradientDriver: state tracking is not enabled.')
        if self.rank == mpi_master():
            self.state_tracker.save_state(filename)
        self.comm.barrier()

    def load_tracking_state(self, filename):
        """Loads a JSON tracking restart into the current driver."""

        if self.rank == mpi_master():
            master_tracker = self._new_state_tracker()
            master_tracker.load_state(filename)
            state = master_tracker.state_dict()
        else:
            state = None

        state = self.comm.bcast(state, root=mpi_master())
        tracker = self._new_state_tracker()
        tracker.load_state_dict(state)
        self.set_state_tracker(tracker)

    def _new_state_tracker(self, target_state=None):
        """Rebuilds a tracker while retaining configured confidence options."""

        old = self.state_tracker
        return OpenQPMRSFStateTracker(
            nstates=self.nstates,
            target_state=(self.target_state_index if target_state is None else
                          int(target_state)),
            minimum_overlap=(0.5 if old is None else old.minimum_overlap),
            maximum_ambiguity_ratio=(
                0.8 if old is None else old.maximum_ambiguity_ratio),
            failure_mode=('error' if old is None else old.failure_mode),
            guard_ground_state=(
                True if old is None else old.guard_ground_state),
            max_stale_reference_steps=(
                10 if old is None else old.max_stale_reference_steps))

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

        tracking_value = grad_dict.get('state_tracking',
                                       grad_dict.get('track_state', None))
        if tracking_value is not None:
            enabled = (tracking_value if isinstance(tracking_value, bool) else
                       str(tracking_value).strip().lower() in
                       ('1', 'true', 'yes', 'on'))
            if enabled:
                self.enable_state_tracking(
                    minimum_overlap=float(
                        grad_dict.get('minimum_tracking_overlap', 0.5)),
                    maximum_ambiguity_ratio=float(
                        grad_dict.get('maximum_ambiguity_ratio', 0.8)),
                    failure_mode=grad_dict.get('tracking_failure', 'error'),
                    guard_ground_state=bool(
                        grad_dict.get('guard_ground_state', True)),
                    max_stale_reference_steps=grad_dict.get(
                        'max_stale_reference_steps', 10))
            else:
                self.disable_state_tracking()

    def compute(self, molecule, *args):
        """
        Performs the calculation of the OpenQP excited-state gradient.

        :param molecule:
            The molecule.
        """

        assert_msg_critical(
            not (self.state_tracking and self.numerical),
            'OpenQPExcitedStateGradientDriver: numerical gradients cannot '
            'update a sequential state-tracking reference safely.')

        computed_new = False
        if self.numerical:
            if self.rank == mpi_master():
                self.compute_numerical(molecule)
            else:
                self.gradient = None
        else:
            # The reuse decision must be identical on every rank, otherwise
            # the ranks disagree about entering the collective section below.
            # Only the master owns the tracker, so it decides and broadcasts.
            if self.rank == mpi_master():
                reuse = self._can_reuse_tracked_result(molecule)
            else:
                reuse = None
            reuse = self.comm.bcast(reuse, root=mpi_master())

            if self.rank == mpi_master():
                if reuse:
                    # compute_energy() already requested this exact tracked
                    # evaluation against the same committed reference.
                    pass
                else:
                    self.gradient = self._compute_analytical_master(molecule)
                    computed_new = True
                    if self.state_tracking:
                        self._tracked_geometry_signature = (
                            self._tracked_cache_key(molecule))
            else:
                self.gradient = None

        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())

        if self.state_tracking:
            payload = ((self.total_energy, self.excited_state_energy,
                        self.state_deriv_index, self.last_tracking_result,
                        computed_new)
                       if self.rank == mpi_master() else None)
            payload = self.comm.bcast(payload, root=mpi_master())
            (self.total_energy, self.excited_state_energy,
             self.state_deriv_index, self.last_tracking_result,
             computed_new) = payload
            if computed_new:
                self.tracking_history.append(self.last_tracking_result)

        self.print_geometry(molecule)
        self.print_gradient(molecule, [self.state_deriv_index])
        if self.state_tracking:
            self._print_tracking_diagnostics(self.last_tracking_result)

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

        if self.state_tracking:
            if self.rank == mpi_master():
                reuse = self._can_reuse_tracked_result(molecule)
            else:
                reuse = None
            reuse = self.comm.bcast(reuse, root=mpi_master())
            if not reuse:
                self.compute(molecule, *args)
            return (self.total_energy
                    if self.rank == mpi_master() else None)

        if self.rank != mpi_master():
            return None

        mol = self._run_state(molecule, need_gradient=False)
        self._record_state_energies(mol)
        self.total_energy = float(mol.energies[self.state_deriv_index])
        self.excited_state_energy = float(
            mol.energies[self.state_deriv_index] - mol.energies[1])

        return self.total_energy

    def _record_state_energies(self, mol):
        """Stores the working reference and the total target-state energies.

        OpenQP returns ``nstate + 1`` entries: slot 0 is the high-spin ROHF
        working reference and slots ``1 .. nstate`` are the total MRSF
        target-state energies.  The working reference is diagnostics only and
        is never a dynamics surface.
        """

        energies = np.asarray(mol.energies, dtype=float).reshape(-1)
        n_targets = max(0, energies.size - 1)

        self.reference_energy = float(energies[0]) if energies.size else None
        self.target_state_energies = energies[1:1 + n_targets].copy()

    def _tracked_cache_key(self, molecule):
        """Cache key for a tracked evaluation: geometry *and* reference.

        Coordinates alone are not a sufficient key.  After a rejected
        optimizer trial the committed reference is rolled back, so the same
        coordinates evaluated against a different committed reference are a
        different calculation and must not return the earlier proposal.
        """

        geometry_signature = self.openqp_driver._get_geometry_signature(
            molecule)
        revision = (self.state_tracker.reference_revision
                    if self.state_tracker is not None else -1)

        return (geometry_signature, int(revision))

    def _can_reuse_tracked_result(self, molecule):
        """Whether the cached tracked energy and gradient are still valid."""

        if not self.state_tracking:
            return False
        if self.gradient is None or self.total_energy is None:
            return False
        if self._tracked_geometry_signature is None:
            return False
        # A staged-but-uncommitted proposal belongs to the cached result; once
        # it is committed or rolled back the key changes with the revision.
        return self._tracked_cache_key(molecule) == \
            self._tracked_geometry_signature

    def commit_tracking_reference(self):
        """Commits the staged tracking proposal (accepted nuclear step).

        Called from :class:`OptimizationEngine` through geomeTRIC's
        accepted-step hook.  Returns ``True`` when the reference advanced.

        An unreliable proposal is refused by the tracker, so a low-overlap or
        ambiguous assignment never becomes the identity that later geometries
        are measured against.
        """

        if not self.state_tracking or self.state_tracker is None:
            return False

        if self.rank == mpi_master():
            committed = self.state_tracker.commit()
            if not committed and self.state_tracker.rejected_commits:
                self.ostream.print_info(
                    'OpenQP MRSF tracking: refused to commit an unreliable '
                    'assignment; the reference stays at the last accepted '
                    'geometry.')
                self.ostream.flush()
        else:
            committed = None

        return self.comm.bcast(committed, root=mpi_master())

    def rollback_tracking_reference(self):
        """Discards the staged proposal (rejected nuclear trial).

        Called from :class:`OptimizationEngine` through geomeTRIC's
        rejected-step hook.  Returns ``True`` when a proposal was discarded.
        """

        if not self.state_tracking or self.state_tracker is None:
            return False

        if self.rank == mpi_master():
            rolled_back = self.state_tracker.rollback()
            if rolled_back:
                # The cached proposal was tied to the discarded reference.
                self._tracked_geometry_signature = None
                self.ostream.print_info(
                    'OpenQP MRSF tracking: nuclear trial rejected; the '
                    'tracking reference was rolled back.')
                self.ostream.flush()
        else:
            rolled_back = None

        return self.comm.bcast(rolled_back, root=mpi_master())

    def _compute_analytical_master(self, molecule):
        if self.state_tracking:
            mol, tracking_result = self._run_tracked_state(molecule)
            state = int(tracking_result.selected_raw_root)
            self.state_deriv_index = state
            # Set by OptimizationEngine on a constrained scan.  The tracker
            # holds this same object in its history, so stamping it here is
            # what later lets the history be aligned with the grid points.
            tracking_result.scan_point = getattr(self, 'current_scan_point',
                                                 None)
            self.last_tracking_result = tracking_result
        else:
            mol = self._run_state(molecule, need_gradient=True)
            state = self.state_deriv_index
        self._record_state_energies(mol)
        self.total_energy = float(mol.energies[state])
        self.excited_state_energy = float(mol.energies[state] - mol.energies[1])

        gradient = np.array(mol.grads[state], dtype=float).reshape(
            (molecule.number_of_atoms(), 3))

        return gradient

    def _run_tracked_state(self, molecule):
        assert_msg_critical(
            self.tddft_type == 'mrsf',
            'OpenQPExcitedStateGradientDriver: native tracking requires MRSF.')
        assert_msg_critical(
            self.state_tracker is not None,
            'OpenQPExcitedStateGradientDriver: missing state tracker.')

        drv = self.openqp_driver
        functional = drv.functional if drv.method == 'dft' else ''
        return drv.run_oqp_tracked_mrsf(
            molecule,
            functional=functional,
            multiplicity=3,
            nstate=self.nstates,
            tracker=self.state_tracker,
            exc_multiplicity=getattr(self.exc_driver,
                                     'tddft_multiplicity', 1))

    def _print_tracking_diagnostics(self, result):
        """Prints the complete root decision at the current geometry."""

        if self.rank != mpi_master() or result is None:
            return

        hartree_to_ev = 27.211386245988
        energies = np.asarray(result.state_energies, dtype=float)
        reference_energy = float(energies[0])

        self.ostream.print_header('OpenQP MRSF State Tracking')
        self.ostream.print_info(
            f'Persistent target identity : {result.target_tracked_root}')
        self.ostream.print_info(
            f'Previous raw root          : {result.previous_raw_root}')
        self.ostream.print_info(
            f'Selected gradient root     : {result.selected_raw_root}')
        self.ostream.print_info(
            'Tracked -> raw permutation : ' +
            ' '.join(str(root) for root in result.tracked_to_raw))
        self.ostream.print_info(
            f'Target |overlap|           : {result.target_score:.8f}')
        self.ostream.print_info(
            f'Second-best |overlap|      : {result.second_score:.8f}')
        self.ostream.print_info(
            f'Ambiguity ratio            : {result.ambiguity_ratio:.8f}')
        self.ostream.print_info(
            f'Minimum assignment score   : {result.global_confidence:.8f}')
        self.ostream.print_info(
            f'Target assignment reliable : {result.is_reliable}')
        self.ostream.print_info(
            f'Raw root changed           : {result.root_changed}')

        if result.ground_state_raw_root is not None:
            guard_text = ('on (S0 excluded from the assignment)'
                          if result.ground_state_guarded else 'OFF')
            self.ostream.print_info(
                f'Ground-state raw root      : '
                f'{result.ground_state_raw_root}  [guard {guard_text}]')
        if result.ground_state_gap_ev is not None:
            self.ostream.print_info(
                f'Selected root above S0     : '
                f'{result.ground_state_gap_ev:.6f} eV')
        if result.ground_state_collision:
            self.ostream.print_info(
                'S1/S0 seam reached         : the tracked character sits on '
                f'raw root {result.unguarded_raw_root} (the ground state)')

        # Where this geometry's starting orbitals came from, and how far the
        # comparison point has fallen behind.  Both are silent by default and
        # both change what the numbers above mean.
        self.ostream.print_info(
            'SCF guess for this geometry: ' +
            ('previous accepted geometry (committed tracking reference)'
             if result.seeded_from_reference else
             'OpenQP default (no reusable reference)'))
        guess_info = getattr(self.openqp_driver, 'last_scf_guess_info', None)
        if guess_info is not None and guess_info.get('seeded'):
            # The reported number is the violation the Loewdin step removed,
            # i.e. how far the previous geometry's orbitals were from
            # orthonormal in *this* geometry's AO metric -- a measure of how
            # far the geometry moved, not a residual error.
            self.ostream.print_info(
                f"  {guess_info['basis_functions']} orbitals re-orthonormalized"
                ' (Loewdin) in this AO metric; largest violation removed: '
                f"{guess_info['orthonormality_error']:.3e}")
        if result.stale_reference_steps:
            self.ostream.print_info(
                'Reference age              : '
                f'{result.stale_reference_steps} accepted step(s) behind '
                '(commits refused); overlaps decay while it stays behind')
        self.ostream.print_blank()
        self.ostream.print_info(
            'Raw  Tracked       Energy / Eh       Relative to S0 / eV')
        for raw_idx, energy in enumerate(energies, start=1):
            tracked_idx = int(
                np.where(result.tracked_to_raw == raw_idx)[0][0] + 1)
            relative_ev = (float(energy) - reference_energy) * hartree_to_ev
            marker = ' <- gradient' if raw_idx == \
                result.selected_raw_root else ''
            self.ostream.print_info(
                f'{raw_idx:3d}  {tracked_idx:7d}  {energy:18.10f}'
                f'  {relative_ev:20.8f}{marker}')

        if self.print_level >= 3 and not result.reference_initialized:
            self.ostream.print_blank()
            self.ostream.print_info(
                'Absolute overlap matrix '
                '(rows=persistent, columns=current raw):')
            for row in np.asarray(result.similarity_matrix):
                self.ostream.print_info(
                    '  ' + ' '.join(f'{value:10.6f}' for value in row))

        for warning in result.warnings:
            self.ostream.print_info('Tracking note: ' + warning)
        self.ostream.print_blank()
        self.ostream.flush()

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
