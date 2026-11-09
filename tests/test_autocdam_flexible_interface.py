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
"""Interface contracts for flexible multisystem fitting (weighted multisystem fits).

Tests marked ``xfail(strict=True)`` describe behavior implemented in a later
development step; strictness makes them fail loudly once they start passing,
so the mark is removed in the step that implements them.
"""
import json
from pathlib import Path

import pytest

from veloxchem.autocdam import AutoCDAM, CDAMCandidate, CalculationResult, Evaluator, FitTarget
from veloxchem.errorhandler import VeloxChemError

STRUCTURES = Path(__file__).parent / "data"


def step2(function):
    """Marks a test of behavior implemented in Step 2."""
    return function


def step3(function):
    """Marks a test of behavior implemented in Step 3."""
    return function


class TestAutoCdamFlexibleInterface:

    def test_node_structure_and_charges_are_consistent(self):
        lines = (STRUCTURES /
                 "autocdam_al_node_fragment.xyz").read_text().splitlines()
        charges = json.loads(
            (STRUCTURES / "autocdam_al_node_fragment_charges.json").read_text())
        assert int(lines[0]) == len(charges["charges_e"]) == 19
        assert sum(charges["charges_e"]) == pytest.approx(
            charges["molecular_charge"])

    @step2
    def test_parameters_accept_bounds_fixed_values_and_aliases(self):
        parameters = {
            "sigma": (1.4, 2.6),
            "epsilon": (0.02, 0.5),
            "metal_dummy_distance": 0.9
        }
        fit = AutoCDAM()
        fit.setup(CDAMCandidate("Al", 3),
                  {"IOD": {
                      "value": 1.9,
                      "tolerance": 0.03
                  }},
                  parameters=parameters)
        specs = {spec.name: spec for spec in fit.parameter_specs}
        assert specs["metal_lj_sigma"].bounds == (
            1.4, 2.6) and specs["metal_lj_sigma"].optimizable
        assert specs["metal_lj_epsilon"].bounds == (0.02, 0.5)
        assert not specs["metal_dummy_distance"].optimizable
        assert specs["metal_dummy_distance"].initial == 0.9

    @step2
    def test_unknown_parameter_is_rejected_before_calculation(self):
        with pytest.raises(VeloxChemError, match="unsupported"):
            AutoCDAM().setup(CDAMCandidate("Al", 3),
                             {"IOD": {
                                 "value": 1.9,
                                 "tolerance": 0.03
                             }},
                             parameters={"neighbor_charge": (0.0, 1.0)})

    @step2
    def test_target_roles_and_normalized_weighted_loss(self):
        targets = [
            FitTarget("A", 1.0, tolerance=0.5, weight=1.0, role="weighted"),
            FitTarget("B", 2.0, tolerance=1.0, weight=3.0, role="weighted"),
            FitTarget("M", 0.0, tolerance=1.0, role="monitor")
        ]
        result = Evaluator().evaluate(
            CalculationResult(values={
                "A": 2.0,
                "B": 2.0,
                "M": 50.0
            }), targets)
        # ((2-1)/0.5)^2 * 1 + 0 * 3, divided by the weight sum 4; the monitor is excluded.
        assert result.objective == pytest.approx(1.0)
        assert result.feasible

    @step2
    def test_required_target_cannot_be_compensated_by_weighted_score(self):
        targets = [
            FitTarget("IOD", 1.90, tolerance=0.03, role="required"),
            FitTarget("CN", 6.0, tolerance=0.2, weight=100.0)
        ]
        result = Evaluator().evaluate(
            CalculationResult(values={
                "IOD": 1.99,
                "CN": 6.0
            }), targets)
        assert not result.feasible
        assert set(result.constraint_violations) == {"IOD"}

    @step3
    def test_framework_system_accepts_veloxchem_structure_and_target_label(
            self):
        import veloxchem as vlx
        from veloxchem.autocdam import FrameworkSystem

        molecule = vlx.Molecule.read_xyz_file(
            str(STRUCTURES / "autocdam_al_node_fragment.xyz"))
        molecule.set_charge(-1)
        charges = json.loads(
            (STRUCTURES /
             "autocdam_al_node_fragment_charges.json").read_text())["charges_e"]
        system = FrameworkSystem(structure=molecule,
                                 charges=charges,
                                 target_metal="Al")
        baseline = system.prepare_baseline()
        assert baseline is system.prepare_baseline()  # prepared once and reused
        assert system.target_index == 0
