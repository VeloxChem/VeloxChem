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
from contextlib import nullcontext
import os

import hashlib

from .veloxchemlib import mpi_master, hartree_in_ev
from .errorhandler import assert_msg_critical
from .serenityscfdriver import (SerenityScfDriver, SerenityCalculationError,
                                parse_serenity_lr_output)

from .resultsio import write_rsp_results_to_hdf5
from .resultsio import write_rsp_full_solution_to_hdf5
from .resultsio import write_scf_results_to_hdf5

try:
    from qcserenity import serenipy as spy
    import qcserenity as qc
except ImportError:
    pass


class SerenityLinearResponseSolver:
    """
    Implements Serenity linear-response solver.

    :param serenity_scf_drv:
        The Serenity SCF driver.

    Instance variables
        - exc_method: Excited-state method (`tda` or `tddft`).
        - nstates: Number of requested states.
        - densfit_j: Density fitting setup for Coulomb terms.
        - grid_accuracy: Main LR integration grid accuracy.
        - small_grid_accuracy: Pre-optimization LR integration grid accuracy.
        - lr_restart_policy: ``'same_geometry'`` (default) or ``'never'``.
        - lr_convergence_retries: Fresh re-solves with doubled Davidson
          cycles before a non-converged spectrum is rejected.

    LRSCF restart versus SCF warm start versus state tracking
        Serenity's LRSCF restart loads the stored excitation vectors of the
        System as raw coefficients in the *current* MO basis and then follows
        those roots.  Between two geometries the MOs change phase, order and
        shape, while the stored vectors do not (X_ai must become s_a s_i X_ai
        under phi -> s phi), so a cross-geometry restart can steer the
        Davidson solver onto the wrong roots and lose states (T1 vanished at
        bend 130 / torsion 10).  A restart is therefore only requested when
        the stored solution was obtained at the same geometry, on the same
        System, from the same SCF solution and with the same response
        settings.  Reusing the previous SCF orbitals (SCF warm start) is a
        separate matter decided by the SCF driver, and comparing old and new
        excited states (transition-density tracking) never feeds a previous
        geometry's vectors into the eigensolver.
    """

    def __init__(self, serenity_scf_drv):
        """
        Initializes Serenity linear-response solver.
        """

        errmsg = 'SerenityLinearResponseSolver: invalid Serenity SCF driver.'
        assert_msg_critical(isinstance(serenity_scf_drv, SerenityScfDriver),
                            errmsg)

        self.scf_driver = serenity_scf_drv
        self.comm = serenity_scf_drv.comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.exc_method = 'tda'
        self.spinflip = False
        self.nstates = 3
        self.nroots = 3
        self.conv_thresh = None
        self.max_cycles = None
        self.max_subspace_dimension = None

        # 'none' avoids missing RI auxiliary basis failures in minimal setups.
        self.densfit_j = serenity_scf_drv.densfit_j
        self.grid_accuracy = serenity_scf_drv.grid_accuracy
        self.small_grid_accuracy = serenity_scf_drv.small_grid_accuracy

        self.lr_restart_policy = 'same_geometry'
        self.lr_convergence_retries = 1

        self._lr_task = None
        self._rsp_results = None
        self._rsp_results_key = None
        self.last_lr_provenance = None

        # Ledger of the LR solution currently stored in the System's restart
        # file.  Updated only after a converged LRSCF solve (by this solver or
        # by the LRSCF step of an excited-state GradientTask) and cleared on
        # any failure, so a restart can never pick up a foreign solution.
        # _last_system_signature is (System name, SCF revision).
        self._last_rsp_geom_signature = None
        self._last_rsp_settings_signature = None
        self._last_system_signature = None

    @staticmethod
    def is_available():
        """
        Returns if Serenity python driver is available.
        """

        return SerenityScfDriver.is_available()

    def set_exc_method(self, exc_method):
        """
        Sets excited-state method (`tda` or `tddft`).
        """

        if self.rank != mpi_master():
            return

        method = str(exc_method).strip().lower()
        if method in ('tda', 'cis'):
            self.exc_method = 'tda'
        elif method in ('tddft', 'rpa'):
            self.exc_method = 'tddft'
        else:
            errmsg = 'SerenityLinearResponseSolver: invalid exc_method '
            errmsg += f'"{exc_method}"'
            assert_msg_critical(False, errmsg)

        self._invalidate_rsp_cache()

    def set_nstates(self, nstates):
        """
        Sets number of requested excited states.
        """

        if self.rank != mpi_master():
            return

        nst = int(nstates)
        assert_msg_critical(nst > 0,
                            'SerenityLinearResponseSolver: nstates must be > 0')
        self.nstates = nst
        self.nroots = nst
        self._invalidate_rsp_cache()

    # Backward-compatible alias.
    def set_nroots(self, nroots):
        self.set_nstates(nroots)

    def update_settings(self, rsp_dict=None, method_dict=None):
        """
        Updates settings in linear-response solver.

        :param rsp_dict:
            The input dictionary of response settings group.
        :param method_dict:
            The input dictionary of method settings group.
        """

        if self.rank != mpi_master():
            return

        if rsp_dict is None:
            rsp_dict = {}

        if method_dict is None:
            method_dict = {}

        if 'tamm_dancoff' in rsp_dict:
            self.set_exc_method('tda' if bool(rsp_dict['tamm_dancoff']) else
                                'tddft')

        if 'method' in rsp_dict:
            self.set_exc_method(rsp_dict['method'])

        if 'spinflip' in rsp_dict:
            self.spinflip = rsp_dict.get('spinflip', False)
            self._invalidate_rsp_cache()

        for key in ('nstates', 'nroots', 'n_eigen', 'neigen', 'nEigen'):
            if key in rsp_dict:
                self.set_nstates(rsp_dict[key])
                break

        if 'conv' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv'])
            self._invalidate_rsp_cache()

        if 'max_cycles' in rsp_dict:
            self.max_cycles = int(rsp_dict['max_cycles'])
            self._invalidate_rsp_cache()

        if 'max_subspace_dimension' in rsp_dict:
            self.max_subspace_dimension = int(rsp_dict['max_subspace_dimension'])
            self._invalidate_rsp_cache()

        if 'densfit_j' in rsp_dict:
            self.densfit_j = str(rsp_dict['densfit_j'])
            self._invalidate_rsp_cache()

        if 'grid_accuracy' in rsp_dict:
            self.grid_accuracy = int(rsp_dict['grid_accuracy'])
            self._invalidate_rsp_cache()

        if 'small_grid_accuracy' in rsp_dict:
            self.small_grid_accuracy = int(rsp_dict['small_grid_accuracy'])
            self._invalidate_rsp_cache()

        for key in ('lr_restart_policy', 'restart_policy'):
            if key in rsp_dict:
                self.set_lr_restart_policy(rsp_dict[key])
                break

        if 'lr_convergence_retries' in rsp_dict:
            self.lr_convergence_retries = int(rsp_dict['lr_convergence_retries'])

        if 'grid_level' in method_dict:
            lvl = int(method_dict['grid_level'])
            self.grid_accuracy = lvl
            self.small_grid_accuracy = lvl
            self._invalidate_rsp_cache()

    def compute(self, molecule, broadcast=True):
        """
        Performs Serenity LR-SCF calculation.

        :param molecule:
            The molecule.
        :param broadcast:
            Broadcast response results from master to all ranks?

        :return:
            The response results dictionary (all ranks if broadcast is True).
        """

        errmsg = 'SerenityLinearResponseSolver: qcserenity is not available. '
        errmsg += 'Please install/build Serenity python bindings.'
        assert_msg_critical(self.is_available(), errmsg)

        failure = None
        if self.rank == mpi_master():
            try:
                rsp_results = self._compute_master(molecule)
            except SerenityCalculationError as error:
                # Master-only callers (broadcast=False) handle it themselves.
                if self.nodes == 1 or not broadcast:
                    raise
                rsp_results = None
                failure = (str(error), error.stage, error.details)
        else:
            rsp_results = None

        if broadcast and self.nodes > 1:
            failure = self.comm.bcast(failure, root=mpi_master())
            if failure is not None:
                raise SerenityCalculationError(*failure)

        if broadcast:
            rsp_results = self.comm.bcast(rsp_results, root=mpi_master())
            self._rsp_results = self._copy_rsp_results(rsp_results)
            return self._copy_rsp_results(rsp_results)

        if self.rank == mpi_master():
            self._rsp_results = self._copy_rsp_results(rsp_results)
            return self._copy_rsp_results(rsp_results)
        return None

    def get_results(self):
        """
        Gets the latest LR-SCF results.
        """

        if self._rsp_results is None:
            return None

        return self._copy_rsp_results(self._rsp_results)

    def get_excitation_energies(self):
        """
        Gets excitation energies in Hartree.
        """

        if self._rsp_results is None:
            return None

        return self._rsp_results['eigenvalues'].copy()

    def get_excitation_energy(self, state_deriv_index):
        """
        Gets one excitation energy by 1-based state index.
        """

        errmsg = 'SerenityLinearResponseSolver: response results are missing.'
        assert_msg_critical(self._rsp_results is not None, errmsg)

        state = int(state_deriv_index)
        assert_msg_critical(state > 0,
                            'SerenityLinearResponseSolver: state index must be > 0')

        nstates = len(self._rsp_results['eigenvalues'])
        errmsg = 'SerenityLinearResponseSolver: requested state index '
        errmsg += f'{state} but only {nstates} state(s) are available.'
        assert_msg_critical(state <= nstates, errmsg)

        return float(self._rsp_results['eigenvalues'][state - 1])

    # Backward-compatible helper.
    def compute_lrresp_master(self, molecule):
        return self._compute_master(molecule)

    def set_lr_restart_policy(self, policy):
        """
        Sets when a stored LRSCF solution may seed the Davidson solver.

        :param policy:
            ``'same_geometry'``: only when the stored solution belongs to the
            same geometry, System, SCF solution and response settings, e.g.
            the LRSCF step of a gradient task right after the spectrum was
            solved at that geometry.  ``'never'``: always solve from scratch.
        """

        label = str(policy).strip().lower()
        assert_msg_critical(
            label in ('same_geometry', 'never'),
            'SerenityLinearResponseSolver: lr_restart_policy must be '
            '"same_geometry" or "never"; cross-geometry restarts are not '
            'supported because Serenity does not transform stored vectors '
            'to the MO basis of the new geometry.')
        self.lr_restart_policy = label

    def _invalidate_rsp_cache(self):
        self._lr_task = None
        self._rsp_results = None
        self._rsp_results_key = None
        self.invalidate_lr_restart_ledger()

    def invalidate_lr_restart_ledger(self):
        """Forgets the stored LR solution; the next solve starts fresh."""

        self._last_rsp_geom_signature = None
        self._last_rsp_settings_signature = None
        self._last_system_signature = None

    def _current_system_signature(self):
        """(System name, SCF revision): identifies the MO basis in use."""

        return (getattr(self.scf_driver, '_system_name', None),
                int(getattr(self.scf_driver, '_scf_revision', 0)))

    def lr_restart_file(self):
        """Path of the System's isolated LRSCF restart file, or ``None``."""

        name = getattr(self.scf_driver, '_system_name', None)
        root = getattr(self.scf_driver, 'scratch_dir', None)
        if not name or not root:
            return None
        suffix = ('res.h5' if self.scf_driver._current_scf_mode == 'restricted'
                  else 'unres.h5')
        return os.path.join(root, name, f'{name}_lrscf.iso.{suffix}')

    def decide_lr_restart(self, settings_signature=None):
        """
        Decides whether the next LRSCF solve may restart from the stored one.

        :param settings_signature:
            Response-settings signature of the solve about to run; defaults
            to this solver's settings.

        :return:
            Dictionary with ``restart`` (bool), ``reason`` and ``policy``.
        """

        if settings_signature is None:
            settings_signature = self._get_rsp_signature()
        current_system = self._current_system_signature()
        geometry = self.scf_driver._active_geom_signature

        restart = False
        if self.lr_restart_policy == 'never':
            reason = 'restart policy is "never"'
        elif (self._last_rsp_geom_signature is None or
              self._last_system_signature is None):
            reason = 'no converged LR solution is stored for this System'
        elif self._last_system_signature[0] != current_system[0]:
            reason = ('the stored LR solution belongs to a different '
                      'Serenity System')
        elif self._last_rsp_geom_signature != geometry:
            reason = 'geometry changed since the stored LR solution'
        elif self._last_system_signature[1] != current_system[1]:
            reason = ('the SCF was re-run since the stored LR solution '
                      '(the MO phases may differ)')
        elif self._last_rsp_settings_signature != settings_signature:
            reason = 'response settings differ from the stored LR solution'
        else:
            path = self.lr_restart_file()
            if path is None or not os.path.isfile(path):
                reason = 'the LRSCF restart file is missing'
            else:
                restart = True
                reason = ('same geometry, System, SCF solution and response '
                          'settings as the stored LR solution')

        return {
            'restart': restart,
            'reason': reason,
            'policy': self.lr_restart_policy,
        }

    def record_lr_solution(self, settings_signature):
        """
        Registers the converged LR solution now held in the restart file.

        :param settings_signature:
            Response-settings signature the solution was obtained with.
        """

        self._last_rsp_geom_signature = self.scf_driver._active_geom_signature
        self._last_rsp_settings_signature = settings_signature
        self._last_system_signature = self._current_system_signature()

    def get_lr_controller(self):
        """Returns the LRSCF controller of the latest converged solve."""

        assert_msg_critical(
            self._lr_task is not None,
            'SerenityLinearResponseSolver: no LRSCF task is available.')
        return self._lr_task.getLRSCFControllers()[0]

    @staticmethod
    def signature_digest(signature):
        """Short stable digest of a settings signature for the records."""

        return hashlib.sha1(repr(signature).encode('utf-8')).hexdigest()[:16]

    def _compute_master(self, molecule):
        # Ensure SCF/system are up-to-date for this geometry.
        self.scf_driver._compute_energy_master(molecule)

        geom_signature = self.scf_driver._active_geom_signature
        system_signature = self._current_system_signature()
        rsp_signature = self._get_rsp_signature()
        results_key = (geom_signature, system_signature, rsp_signature)

        if (self._lr_task is not None and self._rsp_results is not None and
                self._rsp_results_key == results_key):
            results = self._copy_rsp_results(self._rsp_results)
            provenance = dict(results.get('lr_provenance') or {})
            provenance['cache_hit'] = True
            results['lr_provenance'] = provenance
            return results

        # Cross-geometry LR vectors are never used as a Davidson guess; see
        # the class docstring and decide_lr_restart().
        decision = self.decide_lr_restart(rsp_signature)
        attempts = []
        max_cycles = None
        for attempt in range(1 + max(0, int(self.lr_convergence_retries))):
            # A retry after non-convergence always starts from scratch.
            restart = bool(decision['restart']) and attempt == 0
            run = self._run_lr_task(restart, max_cycles)
            attempts.append(run)
            # converged is None when Serenity printed no convergence status;
            # more Davidson cycles cannot fix that.
            if run['converged'] is not False:
                break
            max_cycles = 2 * int(run['max_cycles'])

        if attempts[-1]['converged'] is not True:
            # A non-converged solve may have overwritten the restart file.
            self.invalidate_lr_restart_ledger()
            self._rsp_results = None
            self._rsp_results_key = None
            if attempts[-1]['converged'] is None:
                reason = ('Serenity printed no convergence status for the '
                          'LRSCF solve, so convergence cannot be verified')
            else:
                reason = ('Serenity LRSCF did not converge ("Convergence '
                          'criterion not reached") in '
                          f'{len(attempts)} attempt(s)')
            raise SerenityCalculationError(
                f'{reason}; the spectrum is not usable.', stage='response',
                details={'attempts': attempts})

        transitions = self._get_serenity_transitions()
        eigvecs = self._get_serenity_excitation_vectors()
        self._validate_lr_solution(transitions, eigvecs)

        rsp_results = self._build_rsp_results(transitions)
        rsp_results["exc_method"] = self.exc_method
        rsp_results["eigenvectors"] = eigvecs

        controller = self._lr_task.getLRSCFControllers()[0]
        if self.spinflip:
            self._add_spinflip_metadata(rsp_results, controller)

        self._add_native_rsp_metadata(rsp_results, eigvecs)

        final_run = attempts[-1]
        rsp_results['lr_provenance'] = {
            'geometry_signature': geom_signature,
            'settings_signature': self.signature_digest(rsp_signature),
            'system_name': system_signature[0],
            'scf_revision': system_signature[1],
            'restart_policy': self.lr_restart_policy,
            'restart_requested': bool(final_run['restart_requested']),
            'restart_used': final_run['restart_used'],
            'restart_reason': decision['reason'],
            'converged': True,
            'davidson_iterations': final_run['davidson_iterations'],
            'attempts': attempts,
            'nstates': int(self.nstates),
            'exc_method': self.exc_method,
            'cache_hit': False,
            'scf': self.scf_driver.get_scf_provenance(),
        }
        self.last_lr_provenance = dict(rsp_results['lr_provenance'])

        # Only a converged, validated solution enters the ledger.
        self.record_lr_solution(rsp_signature)
        self._rsp_results = rsp_results
        self._rsp_results_key = results_key

        self._write_final_hdf5(molecule, self._rsp_results)
        self._write_response_vectors(self.scf_driver.get_final_h5py_file(),
                                     eigvecs)

        return self._copy_rsp_results(self._rsp_results)

    def _run_lr_task(self, restart, max_cycles=None):
        """
        Runs one LRSCF solve on the current System and parses its output.

        :param restart:
            Explicit Serenity restart flag for this solve.
        :param max_cycles:
            Optional Davidson cycle limit overriding ``max_cycles``.

        :return:
            Dictionary describing the solve (restart requested/used,
            convergence, Davidson iterations, cycle limit, warnings).
        """

        mode = self.scf_driver._current_scf_mode
        with self.scf_driver._serenity_output_context():
            if mode == 'restricted':
                self._lr_task = spy.LRSCFTask_R(self.scf_driver._system)
            else:
                self._lr_task = spy.LRSCFTask_U(self.scf_driver._system)

        self._configure_lr_task(restart=restart, max_cycles=max_cycles)

        capture = self.scf_driver.capture_serenity_output('response')
        try:
            with capture:
                self._lr_task.run()
        except Exception as error:
            self.invalidate_lr_restart_ledger()
            raise SerenityCalculationError(
                f'Serenity LRSCF failed: {error}', stage='response',
                details={'serenity_output_tail': capture.text[-4000:]}
            ) from error

        parsed = parse_serenity_lr_output(capture.text)
        return {
            'restart_requested': bool(restart),
            # None if Serenity printed no restart message.
            'restart_used': parsed['restart_loaded'] if restart else False,
            'converged': parsed['converged'],
            'davidson_iterations': parsed['davidson_iterations'],
            'max_cycles': int(self._lr_task.settings.maxCycles),
            'warnings': parsed['warnings'],
        }

    def _validate_lr_solution(self, transitions, eigvecs):
        """Rejects spectra with missing, nonfinite or empty roots."""

        nroots = int(transitions.shape[0])
        energies = transitions[:, 0] if nroots else np.zeros(0)
        vectors = np.asarray(eigvecs, dtype=float)
        problems = []
        if nroots < int(self.nstates):
            problems.append(
                f'{nroots} root(s) returned but {self.nstates} requested')
        if not np.all(np.isfinite(energies)):
            problems.append('nonfinite excitation energies')
        if vectors.ndim != 2 or vectors.shape[1] != nroots:
            problems.append('inconsistent excitation-vector shape')
        else:
            if not np.all(np.isfinite(vectors)):
                problems.append('nonfinite excitation vectors')
            elif np.any(np.linalg.norm(vectors, axis=0) <= 1.0e-12):
                problems.append('zero-norm excitation vector')
        if problems:
            self.invalidate_lr_restart_ledger()
            raise SerenityCalculationError(
                'Serenity LRSCF returned an invalid spectrum: ' +
                '; '.join(problems), stage='response',
                details={'problems': problems})

    def _configure_lr_task(self, restart=False, max_cycles=None):
        """
        Configures the LRSCF task.

        :param restart:
            Explicit Serenity restart flag.  Must only be True when
            ``decide_lr_restart`` allows it; never set it for a new geometry.
        :param max_cycles:
            Optional Davidson cycle limit overriding ``max_cycles``.
        """

        # Print level NORMAL: at MINIMUM Serenity silences "Iterative solver
        # converged in N iterations" and the restart messages, which are
        # the only evidence of convergence and restart use.  The output is
        # captured and only echoed in verbose mode.
        if hasattr(self._lr_task, 'generalSettings'):
            self._lr_task.generalSettings.printLevel = (
                spy.GLOBAL_PRINT_LEVELS.NORMAL)

        self._lr_task.settings.method = self.exc_method
        self._lr_task.settings.nEigen = int(self.nstates)
        self._lr_task.settings.restart = bool(restart)
        if self.spinflip:
            self._lr_task.settings.scfstab = 'spinflip'
        if self.conv_thresh is not None:
            self._lr_task.settings.conv = float(self.conv_thresh)

        if max_cycles is None:
            max_cycles = self.max_cycles
        if max_cycles is not None:
            # Serenity treats maxCycles == 1 as an FDEc step and reports
            # neither convergence nor failure.
            assert_msg_critical(
                int(max_cycles) >= 2,
                'SerenityLinearResponseSolver: max_cycles must be at least '
                '2; Serenity reports no convergence status for 1 cycle')
            self._lr_task.settings.maxCycles = int(max_cycles)

        if self.max_subspace_dimension is not None:
            self._lr_task.settings.maxSubspaceDimension = int(
                self.max_subspace_dimension)

        if self.densfit_j is not None:
            self._lr_task.settings.densFitJ = self.densfit_j

        if self.grid_accuracy is not None:
            self._lr_task.settings.grid.accuracy = int(self.grid_accuracy)

        if self.small_grid_accuracy is not None:
            self._lr_task.settings.grid.smallGridAccuracy = int(
                self.small_grid_accuracy)

    def _get_rsp_signature(self, nstates=None, exc_method=None):
        """
        Returns the response-settings signature.

        :param nstates:
            Optional root count overriding ``nstates`` (gradient tasks).
        :param exc_method:
            Optional method overriding ``exc_method`` (SF gradients are
            always SF-TDA).
        """

        return (
            self.exc_method if exc_method is None else str(exc_method),
            bool(self.spinflip),
            int(self.nstates if nstates is None else nstates),
            self.conv_thresh,
            self.max_cycles,
            self.max_subspace_dimension,
            self.densfit_j,
            self.grid_accuracy,
            self.small_grid_accuracy,
        )
    
    def _get_serenity_transitions(self):
        transitions = np.array(self._lr_task.getTransitions(), dtype=float)

        if transitions.size == 0:
            return np.zeros((0, 6), dtype=float)

        if transitions.ndim == 1:
            transitions = transitions.reshape(1, -1)

        return transitions
    
    # def _get_serenity_excitation_vectors(self):
    #     controller = self._lr_task.getLRSCFControllers()[0]

    #     try:
    #         eigvecs_ser = np.array(
    #             controller.getExcitationVectors(spy.LRSCF_TYPE.ISOLATED),
    #             dtype=float,
    #         )
    #     except Exception:
    #         eigvecs_ser = np.array(controller.getExcitationVectors(), dtype=float)

    #     if eigvecs_ser.ndim == 3:
    #         if self.exc_method in ("tda", "cis"):
    #             return eigvecs_ser[0, :, :]
    #         return np.vstack((eigvecs_ser[0, :, :], eigvecs_ser[1, :, :]))

    #     return eigvecs_ser
    
    def _get_serenity_excitation_vectors(self):
        controller = self._lr_task.getLRSCFControllers()[0]

        try:
            eigvecs_ser = np.array(
                controller.getExcitationVectors(spy.LRSCF_TYPE.ISOLATED),
                dtype=float,
            )
        except Exception:
            eigvecs_ser = np.array(controller.getExcitationVectors(), dtype=float)

        if eigvecs_ser.ndim == 3:
            tda_like = self.spinflip or self.exc_method in ("tda", "cis")
            if tda_like:
                return eigvecs_ser[0, :, :]
            return np.vstack((eigvecs_ser[0, :, :], eigvecs_ser[1, :, :]))

        return eigvecs_ser

    def get_spinflip_metadata(self, controller=None):
        """
        Gets spin diagnostics for the current spin-flip roots.

        Serenity exposes the change in spin squared relative to the
        unrestricted reference.  This method adds the reference value and
        classifies every root by the nearest physically compatible
        multiplicity.

        :param controller:
            Optional LRSCF controller. Uses the controller of the latest
            response task when omitted.

        :return:
            Dictionary containing ``delta_s2``, ``reference_s2``,
            ``state_s2``, ``state_multiplicities`` and ``s2_deviation``.
        """

        assert_msg_critical(
            self.spinflip,
            'SerenityLinearResponseSolver: spin diagnostics requested for a '
            'non-spin-flip calculation.')

        if controller is None:
            assert_msg_critical(
                self._lr_task is not None,
                'SerenityLinearResponseSolver: no LRSCF task is available for '
                'spin diagnostics.')
            controller = self._lr_task.getLRSCFControllers()[0]

        assert_msg_critical(
            hasattr(controller, 'getSpinSquared'),
            'SerenityLinearResponseSolver: the Serenity Python bindings do '
            'not expose LRSCFController.getSpinSquared().')

        delta_s2 = np.asarray(
            controller.getSpinSquared('isolated'), dtype=float).reshape(-1)
        excitation_energies = np.asarray(
            controller.getExcitationEnergies('isolated'),
            dtype=float).reshape(-1)
        assert_msg_critical(
            delta_s2.size == excitation_energies.size,
            'SerenityLinearResponseSolver: the number of Delta <S^2> values '
            'does not match the number of spin-flip roots.')

        reference_s2 = self._compute_reference_s2()
        state_s2 = reference_s2 + delta_s2
        assert_msg_critical(
            np.all(np.isfinite(state_s2)),
            'SerenityLinearResponseSolver: nonfinite spin-squared value in '
            'the spin-flip spectrum.')

        scf_results = self.scf_driver._scf_results
        electron_count = None
        if scf_results is not None:
            occ_alpha = np.asarray(
                scf_results.get('occ_alpha', []), dtype=float)
            occ_beta = np.asarray(
                scf_results.get('occ_beta', []), dtype=float)
            if occ_alpha.size and occ_beta.size:
                electron_count = int(
                    np.rint(np.sum(occ_alpha) + np.sum(occ_beta)))

        multiplicities, ideal_s2 = self._classify_multiplicities(
            state_s2, electron_count)
        s2_deviation = np.abs(state_s2 - ideal_s2)

        return {
            'delta_s2': delta_s2,
            'reference_s2': float(reference_s2),
            'state_s2': state_s2,
            'state_multiplicities': multiplicities,
            'ideal_state_s2': ideal_s2,
            's2_deviation': s2_deviation,
        }

    def _add_spinflip_metadata(self, rsp_results, controller):
        metadata = self.get_spinflip_metadata(controller)
        rsp_results.update(metadata)

        energies_ev = np.asarray(
            controller.getExcitationEnergies('isolated'),
            dtype=float).reshape(-1)
        delta_s2 = metadata['delta_s2']
        state_s2 = metadata['state_s2']
        multiplicities = metadata['state_multiplicities']

        self.scf_driver.ostream.print_blank()
        self.scf_driver.ostream.print_header(
            'Serenity Spin-Flip State Analysis')
        self.scf_driver.ostream.print_header(34 * '=')
        self.scf_driver.ostream.print_info(
            f"Reference <S^2>: {metadata['reference_s2']:.6f}")
        self.scf_driver.ostream.print_header(
            ' State   Exc. energy / eV   Delta <S^2>      <S^2>   Mult.')

        for index, (energy, delta, total, multiplicity) in enumerate(
                zip(energies_ev, delta_s2, state_s2, multiplicities), start=1):
            self.scf_driver.ostream.print_header(
                f'{index:6d} {energy:18.6f} {delta:13.6f} '
                f'{total:12.6f} {multiplicity:7d}')

        self.scf_driver.ostream.print_blank()
        self.scf_driver.ostream.flush()

    def _compute_reference_s2(self):
        """
        Computes orbital-based spin squared for the unrestricted reference.
        """

        scf_results = self.scf_driver._scf_results
        assert_msg_critical(
            scf_results is not None,
            'SerenityLinearResponseSolver: SCF results are unavailable for '
            'the reference <S^2> calculation.')

        required = ('C_alpha', 'C_beta', 'S', 'occ_alpha', 'occ_beta')
        assert_msg_critical(
            all(key in scf_results for key in required),
            'SerenityLinearResponseSolver: incomplete SCF results for the '
            'reference <S^2> calculation.')

        coeff_alpha = np.asarray(scf_results['C_alpha'], dtype=float)
        coeff_beta = np.asarray(scf_results['C_beta'], dtype=float)
        overlap = np.asarray(scf_results['S'], dtype=float)
        occ_alpha = np.asarray(scf_results['occ_alpha'], dtype=float)
        occ_beta = np.asarray(scf_results['occ_beta'], dtype=float)

        alpha_indices = np.flatnonzero(occ_alpha > 0.5)
        beta_indices = np.flatnonzero(occ_beta > 0.5)
        n_alpha = int(alpha_indices.size)
        n_beta = int(beta_indices.size)
        spin_z = 0.5 * (n_alpha - n_beta)

        if n_alpha == 0 or n_beta == 0:
            overlap_term = 0.0
        else:
            occupied_overlap = (
                coeff_alpha[:, alpha_indices].T @ overlap @
                coeff_beta[:, beta_indices])
            overlap_term = float(np.sum(np.abs(occupied_overlap)**2))

        # This symmetric UHF expression is equivalent to
        # Sz(Sz + 1) + N_beta - overlap_term when N_alpha >= N_beta.
        return float(spin_z**2 + 0.5 * (n_alpha + n_beta) - overlap_term)

    @staticmethod
    def _classify_multiplicities(state_s2, electron_count=None):
        """
        Classifies spin values by the nearest parity-compatible multiplicity.
        """

        state_s2 = np.asarray(state_s2, dtype=float).reshape(-1)
        first_multiplicity = (
            1 if electron_count is None or electron_count % 2 == 0 else 2)
        candidates = np.arange(first_multiplicity,
                               first_multiplicity + 10,
                               2,
                               dtype=int)
        spins = 0.5 * (candidates - 1)
        candidate_s2 = spins * (spins + 1.0)

        distances = np.abs(state_s2[:, None] - candidate_s2[None, :])
        nearest = np.argmin(distances, axis=1)

        return candidates[nearest], candidate_s2[nearest]
    
    def _build_spinflip_effective_scf_results(self, eigvecs):
        scf = self.scf_driver._scf_results.copy()

        C_a = np.asarray(scf["C_alpha"])
        C_b = np.asarray(scf["C_beta"])
        E_a = np.asarray(scf["E_alpha"])
        E_b = np.asarray(scf["E_beta"])
        occ_a = np.asarray(scf["occ_alpha"])
        occ_b = np.asarray(scf["occ_beta"])

        nocc_a = int(np.count_nonzero(occ_a > 0.0))
        nocc_b = int(np.count_nonzero(occ_b > 0.0))
        nvir_a = C_a.shape[1] - nocc_a
        nvir_b = C_b.shape[1] - nocc_b

        vector_dim = eigvecs.shape[0]

        if vector_dim == nocc_a * nvir_b:
            # alpha -> beta spin flip
            C_eff = np.hstack((C_a[:, :nocc_a], C_b[:, nocc_b:]))
            E_eff = np.hstack((E_a[:nocc_a], E_b[nocc_b:]))
            nocc_eff = nocc_a
            nvir_eff = nvir_b
            direction = "alpha_to_beta"

        elif vector_dim == nocc_b * nvir_a:
            # beta -> alpha spin flip
            C_eff = np.hstack((C_b[:, :nocc_b], C_a[:, nocc_a:]))
            E_eff = np.hstack((E_b[:nocc_b], E_a[nocc_a:]))
            nocc_eff = nocc_b
            nvir_eff = nvir_a
            direction = "beta_to_alpha"

        else:
            raise RuntimeError(
                "Serenity spin-flip vector size does not match either "
                "alpha->beta or beta->alpha transition space."
            )

        occ_eff = np.zeros(nocc_eff + nvir_eff)
        occ_eff[:nocc_eff] = 1.0

        scf["scf_type"] = "restricted"
        scf["C_alpha"] = C_eff
        scf["C_beta"] = C_eff.copy()
        scf["E_alpha"] = E_eff
        scf["E_beta"] = E_eff.copy()
        scf["occ_alpha"] = occ_eff
        scf["occ_beta"] = occ_eff.copy()

        return scf, {
            "num_core": 0,
            "num_valence": nocc_eff,
            "num_virtual": nvir_eff,
            "scf_type": "restricted",
            "response_vector_layout": "tda",
            "serenity_spinflip_direction": direction,
        }
    
    def _add_native_rsp_metadata(self, rsp_results, eigvecs):
        scf_results = self.scf_driver._scf_results
        if scf_results is None:
            return

        occ_alpha = np.array(scf_results["occ_alpha"], dtype=float)
        nocc = int(np.count_nonzero(occ_alpha > 0.0))
        nmo = int(len(occ_alpha))
        nvir = nmo - nocc

        rsp_results["num_core"] = 0
        rsp_results["num_valence"] = nocc
        rsp_results["num_virtual"] = nvir

        nstates = int(rsp_results["number_of_states"])
        details = []

        for state in range(min(nstates, eigvecs.shape[1])):
            if self.scf_driver._current_scf_mode == "restricted" and not self.spinflip:
                details.append(
                    self._get_excitation_details(eigvecs[:, state], nocc, nvir)
                )
            else:
                details.append([])

        rsp_results["excitation_details"] = details
        
    def _get_excitation_details(self, eigvec, nocc, nvir, coef_thresh=0.2):
        n_ov = nocc * nvir
        if eigvec.size not in (n_ov, 2 * n_ov):
            return []

        excitations = []
        de_excitations = []

        for i in range(nocc):
            homo = "HOMO" if i == nocc - 1 else f"HOMO-{nocc - 1 - i}"
            for a in range(nvir):
                lumo = "LUMO" if a == 0 else f"LUMO+{a}"
                ia = i * nvir + a

                c = eigvec[ia]
                if abs(c) > coef_thresh:
                    excitations.append(
                        (abs(c), f"{homo:<8s} -> {lumo:<8s} {c:10.4f}")
                    )

                if eigvec.size == 2 * n_ov:
                    y = eigvec[n_ov + ia]
                    if abs(y) > coef_thresh:
                        de_excitations.append(
                            (abs(y), f"{homo:<8s} <- {lumo:<8s} {y:10.4f}")
                        )

        return [x[1] for x in sorted(excitations, reverse=True)] + [
            x[1] for x in sorted(de_excitations, reverse=True)
        ]

    @staticmethod
    def _build_rsp_results(transitions):
        if transitions.size == 0:
            transitions = np.zeros((0, 6), dtype=float)
        elif transitions.ndim == 1:
            transitions = transitions.reshape(1, -1)

        nroots = transitions.shape[0]
        ncols = transitions.shape[1]

        eigenvalues = transitions[:, 0].copy() if ncols >= 1 else np.zeros(
            nroots)
        osc_len = transitions[:, 1].copy() if ncols >= 2 else np.zeros(nroots)
        osc_vel = transitions[:, 2].copy() if ncols >= 3 else np.zeros(nroots)

        if ncols >= 6:
            rot = transitions[:, 3:6].copy()
        else:
            rot = np.zeros((nroots, 3), dtype=float)

        
        
        return {
            'eigenvalues': eigenvalues,
            'eigenvalues_ev': eigenvalues * hartree_in_ev(),
            'oscillator_strengths': osc_len,
            'oscillator_strengths_velocity': osc_vel,
            'rotatory_strengths': rot,
            'transitions': transitions,
            'number_of_states': int(nroots),
        }

    def _serenity_output_context(self):
        if self.scf_driver.serenity_verbose and not self.scf_driver.ostream.is_muted:
            return nullcontext()
        return qc.redirectOutputToFile(os.devnull)
    
    @staticmethod
    def _copy_rsp_results(rsp_results):
        if rsp_results is None:
            return None

        copied = {}
        for key, val in rsp_results.items():
            if isinstance(val, np.ndarray):
                copied[key] = val.copy()
            else:
                copied[key] = val
        return copied

    def _write_response_vectors(self, fname, eigvecs):
        if fname is None:
            return

        eigvecs = np.array(eigvecs, dtype=float)

        if eigvecs.ndim == 1:
            eigvecs = eigvecs.reshape(-1, 1)

        for state in range(eigvecs.shape[1]):
            write_rsp_full_solution_to_hdf5(fname, eigvecs[:, state], state + 1, self.nstates)
        
    def _write_final_hdf5(self, molecule, rsp_results):
        final_h5_fname = self.scf_driver.get_final_h5py_file()
        if final_h5_fname is None:
            return

        # # Guarantees root metadata and /scf exist before /rsp is appended.
        # self.scf_driver.write_final_hdf5(final_h5_fname, molecule)

        # h5_results = self._build_hdf5_rsp_results(rsp_results)
        # write_lr_rsp_results_to_hdf5(final_h5_fname, h5_results)
        
        if self.spinflip:
            effective_scf_results, sf_rsp_metadata = (
                self._build_spinflip_effective_scf_results(rsp_results["eigenvectors"])
            )

            # write pseudo-restricted /scf for the visualizer
            self.scf_driver.write_final_hdf5(final_h5_fname, molecule)
            write_scf_results_to_hdf5(final_h5_fname, effective_scf_results)

            h5_results = self._build_hdf5_rsp_results(rsp_results)
            h5_results.update(sf_rsp_metadata)
        else:
            self.scf_driver.write_final_hdf5(final_h5_fname, molecule)
            h5_results = self._build_hdf5_rsp_results(rsp_results)

        write_rsp_results_to_hdf5(final_h5_fname, h5_results)

    def _build_hdf5_rsp_results(self, rsp_results):
        transitions = np.array(rsp_results.get('transitions', []), dtype=float)
        eigenvalues = np.array(rsp_results['eigenvalues'], dtype=float)
        oscillator_strengths = np.array(
            rsp_results['oscillator_strengths'], dtype=float)

        nstates = int(rsp_results.get('number_of_states', len(eigenvalues)))

        scf_results = self.scf_driver._scf_results
        if scf_results is not None and 'occ_alpha' in scf_results:
            occ_alpha = np.array(scf_results['occ_alpha'], dtype=float)
            nocc = int(np.count_nonzero(occ_alpha > 0.0))
            nmo = int(len(occ_alpha))
        else:
            nocc = 0
            nmo = 0

        h5_results = {
            'eigenvalues': eigenvalues,
            'oscillator_strengths': oscillator_strengths,
            'number_of_states': nstates,
            'num_core': 0,
            'num_valence': nocc,
            'num_virtual': max(nmo - nocc, 0),
            'serenity_transitions': transitions,
            'exc_method': self.exc_method,
        }
        
        transition_dipole_keys = (
            "electric_transition_dipoles",
            "velocity_transition_dipoles",
            "magnetic_transition_dipoles",
        )
        missing_transition_dipoles = []
        for key in transition_dipole_keys:
            if key in rsp_results:
                h5_results[key] = np.array(rsp_results[key], dtype=float)
            else:
                h5_results[key] = np.zeros((nstates, 3), dtype=float)
                missing_transition_dipoles.append(key)

        h5_results["serenity_transition_dipoles_are_placeholders"] = bool(
            missing_transition_dipoles)

        if "excitation_details" in rsp_results:
            h5_results["excitation_details"] = rsp_results["excitation_details"]

        if 'eigenvalues_ev' in rsp_results:
            h5_results['eigenvalues_ev'] = np.array(
                rsp_results['eigenvalues_ev'], dtype=float)

        if 'oscillator_strengths_velocity' in rsp_results:
            h5_results['oscillator_strengths_velocity'] = np.array(
                rsp_results['oscillator_strengths_velocity'], dtype=float)

        for key in ('delta_s2', 'state_s2', 'state_multiplicities',
                    'ideal_state_s2', 's2_deviation'):
            if key in rsp_results:
                h5_results[key] = np.asarray(rsp_results[key])

        if 'reference_s2' in rsp_results:
            h5_results['reference_s2'] = float(rsp_results['reference_s2'])

        # Native VeloxChem expects 1D rotatory_strengths. Serenity currently
        # exposes extra transition columns; verify their physical meaning before
        # using them for production ECD analysis.
        if transitions.ndim == 2 and transitions.shape[1] >= 6:
            h5_results["serenity_rotatory_strengths_length"] = transitions[:, 3]
            h5_results["serenity_rotatory_strengths_velocity"] = transitions[:, 4]
            h5_results["serenity_rotatory_strengths_modified"] = transitions[:, 5]

            # Native VeloxChem uses one rotatory_strengths vector.
            h5_results["rotatory_strengths"] = transitions[:, 4]

        return h5_results
