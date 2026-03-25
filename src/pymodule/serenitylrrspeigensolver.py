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

from .veloxchemlib import mpi_master, hartree_in_ev
from .errorhandler import assert_msg_critical
from .serenityscfdriver import SerenityScfDriver

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
        self.ostream = serenity_scf_drv.ostream

        self.exc_method = 'tda'
        self.spinflip = False
        self.nstates = 3
        self.nroots = 3
        self.conv_thresh = None
        self.max_cycles = None
        self.max_subspace_dimension = None

        # 'none' avoids missing RI auxiliary basis failures in minimal setups.
        self.densfit_j = 'none'
        self.grid_accuracy = 7
        self.small_grid_accuracy = 7

        self._lr_task = None
        self._rsp_results = None
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
        self.ostream.mute()
    
        if self.rank == mpi_master():
            rsp_results = self._compute_master(molecule)
        else:
            rsp_results = None

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

    def _invalidate_rsp_cache(self):
        self._lr_task = None
        self._rsp_results = None
        self._last_rsp_geom_signature = None
        self._last_rsp_settings_signature = None
        self._last_system_signature = None

    def _compute_master(self, molecule):
        # Ensure SCF/system are up-to-date for this geometry.
        self.scf_driver._compute_energy_master(molecule)

        geom_signature = self.scf_driver._active_geom_signature
        system_signature = self.scf_driver._system_signature
        rsp_signature = self._get_rsp_signature()

        recompute_lr = (self._lr_task is None or
                        self._last_rsp_geom_signature != geom_signature or
                        self._last_rsp_settings_signature != rsp_signature or
                        self._last_system_signature != system_signature)

        if recompute_lr:
            mode = self.scf_driver._current_scf_mode
            with self.scf_driver._serenity_output_context():
                if mode == 'restricted':
                    self._lr_task = spy.LRSCFTask_R(self.scf_driver._system)
                else:
                    self._lr_task = spy.LRSCFTask_U(self.scf_driver._system)

            self._configure_lr_task()

            
            with self.scf_driver._serenity_output_context():
                self._lr_task.run()
            
            transitions = np.array(self._lr_task.getTransitions(), dtype=float)
            self._rsp_results = self._build_rsp_results(transitions)
            self._rsp_results['exc_method'] = self.exc_method

            self._last_rsp_geom_signature = geom_signature
            self._last_rsp_settings_signature = rsp_signature
            self._last_system_signature = system_signature

        return self._copy_rsp_results(self._rsp_results)

    def _configure_lr_task(self):
        if hasattr(self._lr_task, 'generalSettings'):
            self._lr_task.generalSettings.printLevel = (
                spy.GLOBAL_PRINT_LEVELS.MINIMUM)

        self._lr_task.settings.method = self.exc_method
        self._lr_task.settings.nEigen = int(self.nstates)
        
        if self.spinflip:
            self._lr_task.settings.scfstab = 'spinflip'
        if self.conv_thresh is not None:
            self._lr_task.settings.conv = float(self.conv_thresh)

        if self.max_cycles is not None:
            self._lr_task.settings.maxCycles = int(self.max_cycles)

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

    def _get_rsp_signature(self):
        return (
            self.exc_method,
            int(self.nstates),
            self.conv_thresh,
            self.max_cycles,
            self.max_subspace_dimension,
            self.densfit_j,
            self.grid_accuracy,
            self.small_grid_accuracy,
        )

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
        if self.scf_driver.serenity_verbose and not self.ostream.is_muted:
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
