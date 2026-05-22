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
from .serenitylrrspeigensolver import SerenityLinearResponseSolver 

try:
    from qcserenity import serenipy as spy
    import qcserenity as qc
except ImportError:
    pass


class TransitionDensityTracker:
    def __init__(self, lr_resp_driver, target_state=1,
                 min_overlap=0.5, overlap_ratio_bounds=(0.3, 0.6),
                 degeneracy_threshold=None, lrscf_type="isolated"):
        assert_msg_critical(isinstance(lr_resp_driver, SerenityLinearResponseSolver),
                            "TransitionDensityTracker: invalid LR response driver.")
        self.lrresp_driver = lr_resp_driver
        self.target_state = int(target_state)
        self.min_overlap = float(min_overlap)
        self.overlap_ratio_bounds = tuple(overlap_ratio_bounds)
        self.degeneracy_threshold = degeneracy_threshold
        self.lrscf_type = lrscf_type
        self._ref_exc = None
        self._ref_geom = None
        self._ref_coeff = None
        self.last_info = None

    def has_reference(self):
        return self._ref_exc is not None

    def set_target_state(self, state):
        self.target_state = int(state)

    def initialize_reference(self, system, lrscf_controller, mode):
        self._ref_exc = lrscf_controller.getExcitationVectors(self.lrscf_type)
        self._ref_geom = system.getGeometry().copy()
        self._ref_coeff = lrscf_controller.getCoefficients()

    def track(self, system, lrscf_controller, mode, active_reference_state=None):
        """
        Assigns the stored reference state character to one current root.

        Parameters
        ----------
        system
            Current Serenity SystemController.
        lrscf_controller
            Current LRSCFController from the just-finished GradientTask/LRSCFTask.
        mode
            "restricted" or "unrestricted".
        active_reference_state
            1-based root index in the stored reference geometry. If None, use
            self.reference_state.

        Returns
        -------
        dict
            Tracking decision and diagnostics.
        """

        if active_reference_state is None:
            active_reference_state = getattr(self, "reference_state", None)

        if active_reference_state is None:
            active_reference_state = self.target_state

        ref_state = int(active_reference_state)
        assert_msg_critical(ref_state > 0,
                            "TransitionDensityTracker: reference state must be > 0.")

        # First call: there is no previous transition density yet.
        # Store the current state character as the reference.
        if not self.has_reference():
            self.reference_state = ref_state
            self.initialize_reference(system, lrscf_controller, mode)

            self.last_info = {
                "initialized": True,
                "swapped": False,
                "old_state": ref_state,
                "new_state": ref_state,
                "reference_state": ref_state,
                "max_overlap": None,
                "second_overlap": None,
                "overlap_ratio": None,
                "overlap_matrix": None,
                "assignment_confident": True,
                "reference_update_recommended": False,
                "degenerate": False,
            }
            return dict(self.last_info)

        tracker_cls = (spy.TransitionDensityTracking_R
                    if mode == "restricted"
                    else spy.TransitionDensityTracking_U)

        degeneracy_threshold = self.degeneracy_threshold
        if degeneracy_threshold is None:
            degeneracy_threshold = 10.0 * float(self.lrresp_driver.conv_thresh
                                                or 1.0e-6)

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

        # We call Serenity's tracker to build the transition-density overlap matrix.
        # We do not rely on getNewState(), because the Python-side assignment below
        # follows a stored reference character by column.
        swapped: bool = tden.track()

        # ovlp = np.array(tden.getOvlpMatrix(), dtype=float)

        # n_current, n_reference = ovlp.shape
        # assert_msg_critical(
        #     ref_state <= n_reference,
        #     "TransitionDensityTracker: reference state is outside overlap matrix."
        # )

        # # Serenity matrix convention:
        # # rows    = current roots
        # # columns = reference roots
        # #
        # # To follow reference character ref_state, search the corresponding column.
        # column = ovlp[:, ref_state - 1]

        # new_state = int(np.argmax(column)) + 1
        # max_overlap = float(column[new_state - 1])
        # second_overlap = self._second_largest(column)

        # if max_overlap > 0.0:
        #     ratio = float(second_overlap / max_overlap)
        # else:
        #     ratio = np.inf

        # swapped = (new_state != ref_state)

        # degenerate = self._states_are_degenerate(
        #     lrscf_controller,
        #     ref_state,
        #     new_state,
        #     degeneracy_threshold,
        # )

        # if degenerate:
        #     # Same as Serenity's own conservative behavior: do not relabel roots
        #     # inside the degeneracy window.
        #     new_state = ref_state
        #     swapped = False

        # lo, hi = self.overlap_ratio_bounds

        # assignment_confident = (
        #     max_overlap >= self.min_overlap and
        #     ratio <= hi and
        #     not degenerate
        # )

        # # Serenity's geometry-optimization logic updates the reference only in a
        # # controlled overlap-ratio window, and not directly after a switch.
        # reference_update_recommended = (
        #     not swapped and
        #     not degenerate and
        #     max_overlap >= self.min_overlap and
        #     lo <= ratio <= hi
        # )

        self.last_info = {
            "initialized": False,
            "swapped": swapped,
            "old_state": ref_state,
            "new_state": int(tden.getNewState()),
            "reference_state": ref_state,
        }

        return dict(self.last_info)

    def accept_reference(self, system, lrscf_controller, mode,
                     reference_state=None):
        if reference_state is not None:
            self.reference_state = int(reference_state)

        self._ref_exc = lrscf_controller.getExcitationVectors(self.lrscf_type)
        self._ref_geom = system.getGeometry().copy()
        self._ref_coeff = lrscf_controller.getCoefficients()
    
    def initialize_reference(self, system, lrscf_controller, mode):
        self._ref_exc = lrscf_controller.getExcitationVectors(self.lrscf_type)
        self._ref_geom = system.getGeometry().copy()
        self._ref_coeff = lrscf_controller.getCoefficients()

    def _second_largest(self, values):
        vals = np.asarray(values, dtype=float)

        if vals.size < 2:
            return 0.0

        imax = int(np.argmax(vals))
        remaining = np.delete(vals, imax)

        return float(np.max(remaining))


    def _states_are_degenerate(self, lrscf_controller, state_a, state_b, threshold):
        if state_a == state_b:
            return False

        # Python LRSCFController returns eV, so convert to Hartree.
        exc_ha = (np.array(lrscf_controller.getExcitationEnergies("isolated"),
                        dtype=float) / hartree_in_ev())

        return abs(exc_ha[state_a - 1] - exc_ha[state_b - 1]) < float(threshold)

    # def _get_current_coefficients(self, system, mode):
    #     if mode == "restricted":
    #         #Here get the coefficients from LRSCF controller
    #         return system.getElectronicStructure_R().coeff()
    #     es = system.getElectronicStructure_U()
    #     return np.stack((es.alphaCoeff(), es.betaCoeff()))
        