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
"""Geometry-fidelity RMSD helpers: Kabsch shape, angle multiset, assignment."""
import math

import numpy as np
import pytest

from veloxchem.autocdam import Evaluator, FitTarget, CalculationResult
from veloxchem.autocdam.geometry import angle_rmsd_deg, assign_to_reference, kabsch_rmsd, pairwise_angles_deg
from veloxchem.autocdam.models import AQUEOUS_PROPERTIES, GEOMETRY_RMSD_PROPERTIES
from veloxchem.errorhandler import VeloxChemError

OCTAHEDRON = np.array([
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, -1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
])


def _rot_z(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


# ---------------------------------------------------------------------------
# role: a factor is either a fitting target (optimized) or analysis-only.
# ---------------------------------------------------------------------------


class TestAutoCdamGeometryRmsd:

    def test_pairwise_angles_of_octahedron(self):
        angles = np.sort(pairwise_angles_deg(OCTAHEDRON))
        # 12 edges at 90 deg, 3 diagonals at 180 deg
        assert len(angles) == 15
        assert np.allclose(angles[:12], 90.0, atol=1e-6)
        assert np.allclose(angles[12:], 180.0, atol=1e-6)

    def test_angle_rmsd_zero_for_ideal(self):
        ideal = np.sort(pairwise_angles_deg(OCTAHEDRON))
        assert angle_rmsd_deg(OCTAHEDRON, ideal) == pytest.approx(0.0, abs=1e-9)

    def test_kabsch_rmsd_invariant_to_rotation(self):
        rotated = OCTAHEDRON @ _rot_z(0.7).T
        # pure rotation -> identical shape -> ~0 RMSD after superposition
        assert kabsch_rmsd(rotated, OCTAHEDRON) == pytest.approx(0.0, abs=1e-9)

    def test_kabsch_rmsd_detects_distortion(self):
        distorted = OCTAHEDRON.copy()
        distorted[0] *= 1.5  # stretch one vertex outward
        rmsd = kabsch_rmsd(distorted, OCTAHEDRON)
        assert rmsd > 0.1

    def test_assignment_recovers_permutation(self):
        perm = np.array([2, 0, 5, 1, 4, 3])
        shuffled = OCTAHEDRON[perm]
        got = assign_to_reference(shuffled, OCTAHEDRON)
        # each shuffled[i] should map back to its original reference row perm[i]
        assert list(got) == list(perm)

    def test_geometry_targets_registered(self):
        for name in GEOMETRY_RMSD_PROPERTIES:
            assert AQUEOUS_PROPERTIES[name].job == "md"

    def test_analysis_role_target_excluded_from_objective(self):
        calc = CalculationResult(values={"IOD": 2.0, "SHELL_SHAPE_RMSD": 0.5})
        fit_only = [FitTarget("IOD", 1.9, 0.05, role="objective")]
        plus_analysis = fit_only + [
            FitTarget("SHELL_SHAPE_RMSD", 0.0, 0.1, role="analysis")
        ]
        obj_fit = Evaluator().evaluate(calc, fit_only).objective
        obj_both = Evaluator().evaluate(calc, plus_analysis).objective
        # the analysis factor is measured but must not change the objective
        assert obj_both == pytest.approx(obj_fit)
        res = Evaluator().evaluate(calc, plus_analysis)
        assert "SHELL_SHAPE_RMSD" in res.metrics  # still reported
        assert res.target_results["SHELL_SHAPE_RMSD"]["role"] == "analysis"

    def test_objective_role_target_included_in_objective(self):
        calc = CalculationResult(values={"IOD": 2.0, "SHELL_SHAPE_RMSD": 0.5})
        base = [FitTarget("IOD", 1.9, 0.05, role="objective")]
        plus_fit = base + [
            FitTarget("SHELL_SHAPE_RMSD", 0.0, 0.1, role="objective")
        ]
        assert Evaluator().evaluate(calc, plus_fit).objective > \
            Evaluator().evaluate(calc, base).objective

    def test_default_role_is_objective_and_bad_role_rejected(self):
        assert FitTarget("IOD", 1.9, 0.05).is_objective is True
        with pytest.raises(VeloxChemError):
            FitTarget("IOD", 1.9, 0.05, role="diagnostic")
