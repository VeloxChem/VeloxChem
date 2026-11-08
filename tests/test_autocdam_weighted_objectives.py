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
"""Step 2: target roles, normalized weighted loss, constraints, and parameter selection."""
import math

import pytest

from veloxchem.autocdam import AutoCDAM, CalculationResult, CDAMCandidate, Evaluator, FitTarget
from veloxchem.autocdam.evaluator import within_tolerance


class SurfaceCalculator:
    """IOD rises with sigma and CN with dummy charge; records every request."""

    def __init__(self):
        self.calls = []

    def calculate(self, batch, parameters, targets):
        self.calls.append((dict(parameters.values), [t.name for t in targets]))
        sigma = parameters["metal_lj_sigma"]
        values = {
            "IOD": 0.9 + 0.5 * sigma,
            "CN": 4.0 + 4.0 * parameters["dummy_charge"],
            "SHELL": 0.01
        }
        return CalculationResult(
            values={t.name: values[t.name] for t in targets})


def fit_with(targets, parameters=None, **kwargs):
    calculator = SurfaceCalculator()
    parameters = parameters or {
        "sigma": (1.4, 2.6),
        "dummy_charge": (0.2, 0.7),
        "metal_dummy_distance": 0.9
    }
    search = {"method": "staged", "max_iterations": 2, "samples_per_stage": 6}
    fit = AutoCDAM(candidate=CDAMCandidate("Al", 3),
                   parameters=parameters,
                   fit_targets=targets,
                   calculator=calculator,
                   search=search,
                   **kwargs)
    return fit, calculator


class TestAutoCdamWeightedObjectives:

    def test_scale_normalizes_error_but_tolerance_decides_convergence(self):
        target = FitTarget("D", 1.0, tolerance=0.1, scale=0.5)
        evaluation = Evaluator().evaluate(CalculationResult(values={"D": 1.08}),
                                          [target])
        assert evaluation.errors["D"] == pytest.approx(0.16)
        assert evaluation.objective == pytest.approx(0.16**2)
        assert within_tolerance(evaluation, [target])

    def test_loss_is_unit_independent_when_tolerance_uses_the_same_unit(self):
        angstrom = [FitTarget("D", 1.90, tolerance=0.03, unit="angstrom")]
        nanometer = [FitTarget("D", 0.190, tolerance=0.003, unit="nm")]
        in_angstrom = Evaluator().evaluate(
            CalculationResult(values={"D": 1.95}), angstrom)
        in_nm = Evaluator().evaluate(CalculationResult(values={"D": 0.195}),
                                     nanometer)
        assert in_angstrom.objective == pytest.approx(in_nm.objective)

    def test_zero_weight_required_target_constrains_without_entering_loss(self):
        targets = [
            FitTarget("IOD", 1.9, tolerance=0.03, weight=0.0, role="required"),
            FitTarget("CN", 6.0, tolerance=0.2)
        ]
        evaluation = Evaluator().evaluate(
            CalculationResult(values={
                "IOD": 2.5,
                "CN": 6.1
            }), targets)
        assert evaluation.objective == pytest.approx(0.25)
        assert not evaluation.feasible and "IOD" in evaluation.constraint_violations

    def test_legacy_roles_still_work(self):
        targets = [
            FitTarget("A", 0.0, 1.0, role="objective"),
            FitTarget("B", 0.0, 1.0, role="analysis")
        ]
        assert targets[0].is_objective and targets[0].kind == "weighted"
        assert targets[1].is_monitor and not targets[1].is_objective

    def test_failed_calculation_is_infeasible_with_finite_penalty(self):
        evaluation = Evaluator().evaluate(
            CalculationResult(failed=True, failure_reason="MD exploded"),
            [FitTarget("IOD", 1.9, 0.03, role="required")])
        assert evaluation.failed and not evaluation.feasible
        assert math.isfinite(
            evaluation.objective) and "MD exploded" in evaluation.failure_reason

    def test_missing_monitor_value_is_not_a_failure(self):
        targets = [
            FitTarget("IOD", 1.9, 0.03),
            FitTarget("SHELL", 0.0, 0.1, role="monitor")
        ]
        evaluation = Evaluator().evaluate(
            CalculationResult(values={"IOD": 1.9}), targets)
        assert not evaluation.failed and "SHELL" not in evaluation.metrics

    @pytest.mark.parametrize("bad, message", [
        ({
            "role": "preferred"
        }, "role"),
        ({
            "weight": -1.0
        }, "weight"),
        ({
            "tolerance": 0.0
        }, "tolerance"),
        ({
            "scale": -2.0
        }, "scale"),
        ({
            "atomz": [0, 1]
        }, "unsupported keys"),
    ])
    def test_invalid_target_inputs_are_rejected(self, bad, message):
        with pytest.raises(ValueError, match=message):
            FitTarget.from_mapping("T", {"value": 1.0, "tolerance": 0.1, **bad})

    def test_target_selection_keys_are_kept_in_metadata(self):
        target = FitTarget.from_mapping(
            "metal_oxygen", {
                "system": "mof",
                "property": "distance",
                "atoms": ["AL1", "O3"],
                "value": 1.93,
                "tolerance": 0.05
            })
        assert target.system == "mof" and target.property_name == "distance"
        assert target.metadata["atoms"] == ["AL1", "O3"]
        assert FitTarget("IOD", 1.9, 0.03).property_name == "IOD"

    @pytest.mark.parametrize("parameters, message", [
        ({
            "sigma": (2.6, 1.4)
        }, "lower < upper"),
        ({
            "sigma": (1.0, 2.0, 3.0)
        }, "lower, upper"),
        ({
            "sigma": 2.0
        }, "at least one parameter"),
        ({
            "sigma": (1.4, 2.6),
            "metal_lj_sigma": (1.5, 2.5)
        }, "more than once"),
    ])
    def test_invalid_parameter_inputs_are_rejected(self, parameters, message):
        with pytest.raises(ValueError, match=message):
            AutoCDAM(candidate=CDAMCandidate("Al", 3),
                     parameters=parameters,
                     fit_targets={"IOD": {
                         "value": 1.9,
                         "tolerance": 0.03
                     }})

    def test_parameters_and_fit_parameters_are_exclusive(self):
        with pytest.raises(ValueError, match="either"):
            AutoCDAM(candidate=CDAMCandidate("Al", 3),
                     parameters={"sigma": (1.4, 2.6)},
                     fit_parameters={"metal_lj_sigma": {
                         "bounds": [1.0, 3.0]
                     }},
                     fit_targets={"IOD": {
                         "value": 1.9,
                         "tolerance": 0.03
                     }})

    def test_bounded_start_uses_template_value_or_midpoint(self):
        targets = {"IOD": {"value": 1.9, "tolerance": 0.03}}
        parameters = {"sigma": (1.4, 2.6), "dummy_charge": (0.6, 0.7)}
        fit, _ = fit_with(targets, parameters=parameters)
        specs = {spec.name: spec for spec in fit.parameter_specs}
        assert specs[
            "metal_lj_sigma"].initial == 2.0  # template default is inside
        assert specs["dummy_charge"].initial == pytest.approx(
            0.65)  # template 0.5 is outside

    def test_fixed_parameters_never_change_and_monitors_are_not_searched(self):
        fit, calculator = fit_with({
            "IOD": {
                "value": 1.9,
                "tolerance": 0.03,
                "role": "required"
            },
            "CN": {
                "value": 6.0,
                "tolerance": 0.2,
                "weight": 2.0
            },
            "SHELL": {
                "value": 0.0,
                "tolerance": 0.1,
                "role": "monitor"
            },
        })
        result = fit.fit()
        search_calls = calculator.calls[:-1]
        assert all(
            call[0]["metal_dummy_distance"] == 0.9 for call in calculator.calls)
        assert all(call[1] == ["IOD", "CN"] for call in search_calls)
        # one final monitor calculation
        assert calculator.calls[-1][1] == ["SHELL"]
        assert result.parameters["metal_lj_sigma"] == calculator.calls[-1][0][
            "metal_lj_sigma"]
        assert result.diagnostics["monitor"]["SHELL"]["value"] == 0.01
        assert result.diagnostics["monitor"]["SHELL"][
            "within_tolerance"] is True

    def test_monitor_with_evaluate_always_is_calculated_for_every_candidate(
            self):
        fit, calculator = fit_with({
            "IOD": {
                "value": 1.9,
                "tolerance": 0.03
            },
            "SHELL": {
                "value": 0.0,
                "tolerance": 0.1,
                "role": "monitor",
                "evaluate": "always"
            },
        })
        fit.fit()
        assert all(call[1] == ["IOD", "SHELL"] for call in calculator.calls)

    def test_only_monitors_cannot_be_fitted(self):
        with pytest.raises(ValueError, match="required or weighted"):
            fit_with(
                {"SHELL": {
                    "value": 0.0,
                    "tolerance": 0.1,
                    "role": "monitor"
                }})

    def test_unreachable_required_target_is_never_reported_converged(self):
        fit, _ = fit_with({
            "IOD": {
                "value": 5.0,
                "tolerance": 0.01,
                "role": "required"
            },  # beyond the bounds
            "CN": {
                "value": 6.0,
                "tolerance": 0.5
            },
        })
        result = fit.fit()
        assert not result.converged and result.status != "success"
        assert result.diagnostics["feasible"] is False
        assert "IOD" in result.diagnostics["constraint_violations"]
