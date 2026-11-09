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
"""Behavioral tests: real candidate regeneration, refitting and observed stopping."""
import math

import pytest

from veloxchem.autocdam import (AutoCDAM, CalculationResult, CDAMCandidate,
                                CoordinateGenerator, Evaluator, FitTarget,
                                ParameterSet, ParameterSpec)
from veloxchem.autocdam.models import RESULT_GROUP, read_h5_group
from veloxchem.autocdam.search import (CHECKPOINT_GROUP, CalculationInput,
                                       FittingDriver)


def candidate():
    return CDAMCandidate("Al",
                         3,
                         water_model="tip3p",
                         metal_lj_sigma=2.0,
                         metal_lj_epsilon=0.2)


class SyntheticCalculator:
    """Known property surface, independent of the requested reference."""

    def __init__(self):
        self.batches = []

    def calculate(self, batch, parameters, targets):
        self.batches.append(batch)
        sigma = batch.system.nonbonded_parameters["metal_lj_sigma_angstrom"]
        assert sigma == parameters["metal_lj_sigma"]
        return CalculationResult(values={"IOD": sigma, "CN": 2 * sigma})


def make_driver(*,
                max_iterations=3,
                checkpoint_file=None,
                tolerance=1e-12,
                **options):
    generator = CoordinateGenerator()
    calculator = SyntheticCalculator()
    specs = [ParameterSpec("metal_lj_sigma", 1.0, (0.5, 4.0))]
    targets = [
        FitTarget("IOD", 2.345, tolerance),
        FitTarget("CN", 4.69, 2 * tolerance)
    ]

    def inputs(proposals, ids):
        return [
            CalculationInput(b) for b in generator.generate_parameter_batch(
                candidate(), proposals, sample_ids=ids)
        ]

    return FittingDriver(specs,
                         targets,
                         calculator,
                         Evaluator(),
                         inputs,
                         max_iterations=max_iterations,
                         samples_per_stage=8,
                         seed=42,
                         checkpoint_file=checkpoint_file,
                         **options), calculator


class TestAutoCdamStagedFitting:

    def test_stages_calculate_new_candidates_refit_and_shrink(self):
        driver, calculator = make_driver()
        result = driver.fit()
        assert result.iterations == 3
        assert result.diagnostics["stop_reason"] == "max_iterations"
        assert not result.converged
        assert result.status == "warning"
        assert result.diagnostics["n_calculations"] == 24
        assert len(driver.surrogate.models) == 2
        assert result.history[1]["predicted_candidate"] is not None
        assert result.history[-1]["best_objective"] <= result.history[0][
            "best_objective"]
        for before, after in zip(result.history, result.history[1:]):
            old_lo, old_hi = before["region_bounds"]["metal_lj_sigma"]
            new_lo, new_hi = after["region_bounds"]["metal_lj_sigma"]
            assert new_hi - new_lo == pytest.approx((old_hi - old_lo) * 0.5)
        points = [b.parameters["metal_lj_sigma"] for b in calculator.batches]
        assert len(set(points)) == len(points)
        assert result.diagnostics["best_properties"][
            "IOD"] == result.parameters["metal_lj_sigma"]

    def test_threshold_uses_observations_and_completes_whole_stage(self):
        driver, calculator = make_driver(max_iterations=10, tolerance=10)
        result = driver.fit()
        assert result.converged
        assert result.iterations == 1
        assert len(calculator.batches) == 8
        assert result.diagnostics["stop_reason"] == "threshold"

    def test_resume_keeps_samples_and_replays_same_next_stage(self, tmp_path):
        path = tmp_path / "fit"
        first, _ = make_driver(max_iterations=1, checkpoint_file=path)
        first.fit()
        resumed, calculator = make_driver(max_iterations=3,
                                          checkpoint_file=path)
        result = resumed.fit(resume=True)
        uninterrupted, _ = make_driver(max_iterations=3)
        expected = uninterrupted.fit()
        assert result.parameters.values == expected.parameters.values
        assert result.history == expected.history
        assert len(calculator.batches) == 16
        incompatible, _ = make_driver(max_iterations=4,
                                      checkpoint_file=path,
                                      tolerance=2)
        with pytest.raises(ValueError, match="settings differ"):
            incompatible.fit(resume=True)

    def test_resume_without_a_checkpoint_is_an_error(self, tmp_path):
        missing, _ = make_driver(checkpoint_file=tmp_path / "absent")
        with pytest.raises(ValueError,
                           match="restart requires an existing checkpoint"):
            missing.fit(resume=True)

    def test_checkpoint_groups_share_one_run_file(self, tmp_path):
        path = tmp_path / "run"
        group = "autocdam/checkpoint_a"
        first, _ = make_driver(max_iterations=1,
                               checkpoint_file=path,
                               checkpoint_group=group)
        result = first.fit()
        assert result.checkpoint == str(tmp_path / "run.h5")
        assert result.checkpoint_group == group
        assert read_h5_group(path,
                             group)["history"] == result.to_dict()["history"]
        assert read_h5_group(path, CHECKPOINT_GROUP) is None
        assert read_h5_group(path, RESULT_GROUP) is None  # driver-level save
        other, _ = make_driver(max_iterations=1, checkpoint_file=path)
        with pytest.raises(ValueError,
                           match="restart requires an existing checkpoint"):
            other.fit(resume=True)
        resumed, calculator = make_driver(max_iterations=1,
                                          checkpoint_file=path,
                                          checkpoint_group=group)
        assert resumed.fit(resume=True).iterations == 1
        assert calculator.batches == []

    def test_failed_and_unreviewed_results_never_win_or_train(self):

        class FailedCalculator:

            def calculate(self, batch, parameters, targets):
                if parameters["x"] < 0.3:
                    raise RuntimeError("unstable candidate")
                return CalculationResult(
                    values={"IOD": 0.0},
                    metadata={"IOD": {
                        "quality_status": "unreviewed"
                    }})

        driver = FittingDriver(
            [ParameterSpec("x", 0.0, (0.0, 1.0))], [FitTarget("IOD", 0, 1)],
            FailedCalculator(),
            Evaluator(),
            lambda points, ids: [
                CalculationInput(candidate().to_generated_batch())
                for p in points
            ],
            max_iterations=2,
            samples_per_stage=6)
        result = driver.fit()
        assert not result.converged
        assert result.status == "failed"
        assert math.isinf(result.objective)
        assert result.diagnostics["quality_counts"]["unreviewed"] > 0
        assert result.diagnostics["quality_counts"]["rejected"] > 0
        assert not driver.surrogate.models

    def test_candidate_generation_error_does_not_discard_neighbors(self):
        calculator = SyntheticCalculator()
        targets = {"IOD": {"value": 2.0, "tolerance": 100}}
        parameters = {"metal_lj_sigma": {"initial": 1.0, "bounds": [-2, 4]}}
        cdam = AutoCDAM()
        cdam.calculator = calculator
        cdam.update_settings({
            "search_method": "staged",
            "max_iterations": 1,
            "samples_per_stage": 12
        })
        result = cdam.fit(candidate(), targets, fit_parameters=parameters)
        assert result.converged
        assert result.diagnostics["quality_counts"]["rejected"] > 0
        assert result.diagnostics["quality_counts"]["usable"] > 0

    def test_array_properties_have_forward_models(self):

        class ArrayCalculator:

            def calculate(self, batch, parameters, targets):
                x = parameters["x"]
                return CalculationResult(values={"energy": [x, 2 * x]})

        driver = FittingDriver(
            [ParameterSpec("x", 0.1, (0, 1))],
            [FitTarget("energy", [0.345, 0.69], 1e-10)],
            ArrayCalculator(),
            Evaluator(),
            lambda points, ids: [
                CalculationInput(candidate().to_generated_batch())
                for p in points
            ],
            max_iterations=2,
            samples_per_stage=5)
        result = driver.fit()
        assert result.iterations == 2
        predictions = driver.surrogate.predict([ParameterSet({"x": 0.5})])
        prediction = predictions["energy"][0]
        assert len(prediction.mean) == 2
        assert len(prediction.standard_deviation) == 2

    def test_unknown_live_parameter_is_rejected(self):
        with pytest.raises(ValueError, match="unsupported CDAM parameters"):
            CoordinateGenerator().generate_parameter_batch(
                candidate(), [{
                    "epsilon_typo": 1
                }])

    def test_tolerance_solution_is_selected_even_if_weighted_loss_is_larger(
            self):
        from veloxchem.autocdam.search import ParameterPropertySample
        driver, _ = make_driver()
        driver.targets = [FitTarget("IOD", 0, 1), FitTarget("CN", 0, 1)]
        driver.table.add(
            ParameterPropertySample(
                "low_loss", ParameterSet({"metal_lj_sigma": 1}),
                CalculationResult(values={
                    "IOD": 1.1,
                    "CN": 0
                }), "usable"))
        driver.table.add(
            ParameterPropertySample(
                "all_tolerances", ParameterSet({"metal_lj_sigma": 2}),
                CalculationResult(values={
                    "IOD": 0.9,
                    "CN": 0.9
                }), "usable"))
        result = driver.fit()
        assert result.converged
        assert result.diagnostics["best_sample_id"] == "all_tolerances"
        assert result.iterations == 0

    def test_warning_does_not_block_fitting_and_is_retained(self):
        driver, calculator = make_driver(tolerance=10)
        original_calculate = calculator.calculate

        def calculate_with_warning(*args):
            result = original_calculate(*args)
            result.warnings.append(
                "Short sampling: inspect equilibration separately.")
            return result

        calculator.calculate = calculate_with_warning
        result = driver.fit()
        assert result.converged
        assert len(driver.surrogate.models) == 2
        assert result.diagnostics["quality_counts"]["usable"] == 8
        assert len(result.warnings) == 8
        assert all(
            sample.calculation.warnings for sample in driver.table.samples)

    @pytest.mark.parametrize(
        "value", [float("nan"), float("inf"), [1.0, float("nan")]])
    def test_nonfinite_observations_are_rejected_even_with_warnings(
            self, value):
        from veloxchem.autocdam.search import default_quality_check

        result = CalculationResult(values={"IOD": value},
                                   warnings=["Diagnostic warning"])
        assert default_quality_check(result) == "rejected"

    def test_surrogate_learns_surface_and_predicts_an_unsampled_point(self):
        from veloxchem.autocdam.search import ParameterPropertySample, ParameterPropertyTable
        from veloxchem.autocdam.search import SklearnGaussianProcessSurrogate

        table = ParameterPropertyTable()
        for index in range(9):
            x = index / 8
            table.add(
                ParameterPropertySample(
                    str(index),
                    ParameterSet({"x": x}),
                    CalculationResult(values={"IOD": 2 + math.sin(2 * x)}),
                    "usable",
                ))
        surrogate = SklearnGaussianProcessSurrogate()
        surrogate.fit(table, ["IOD"])
        prediction = surrogate.predict([ParameterSet({"x": 0.37})])["IOD"][0]
        assert prediction.mean == pytest.approx(2 + math.sin(0.74), abs=0.02)
        assert math.isfinite(prediction.standard_deviation)
        assert prediction.standard_deviation >= 0
