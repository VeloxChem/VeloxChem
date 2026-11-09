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
"""Step 6: Gaussian-process Bayesian search on synthetic surfaces with known optima."""
import math

import numpy as np
import pytest

from veloxchem.autocdam import (AutoCDAM, CalculationResult, CDAMCandidate,
                                samples_table)
from veloxchem.autocdam.models import read_h5_group
from veloxchem.errorhandler import VeloxChemError

PARAMETERS = {
    "sigma": (1.4, 2.6),
    "dummy_charge": (0.2, 0.7),
    "metal_dummy_distance": 0.9
}


class Surface:
    """IOD from sigma, CN from dummy charge; optimum at sigma 2.31, q 0.29."""

    def __init__(self, noise=0.0, fail_above_sigma=None, seed=0):
        self.noise, self.fail_above_sigma = noise, fail_above_sigma
        self.rng = np.random.default_rng(seed)
        self.points = []

    def calculate(self, batch, parameters, targets, **kwargs):
        sigma, q = parameters["metal_lj_sigma"], parameters.values.get(
            "dummy_charge", 0.5)
        self.points.append((sigma, q))
        if self.fail_above_sigma is not None and sigma > self.fail_above_sigma:
            return CalculationResult(failed=True,
                                     failure_reason="synthetic instability")
        values = {
            "IOD": 0.9 + 0.5 * sigma + self.rng.normal(0, self.noise),
            "CN": 4.0 + 4.0 * q + self.rng.normal(0, 5 * self.noise),
            "D": 3.0 - 0.6 * sigma
        }
        return CalculationResult(
            values={t.name: values[t.name] for t in targets})


TARGETS = {
    "IOD": {
        "value": 2.055,
        "tolerance": 0.01,
        "role": "required"
    },
    "CN": {
        "value": 5.16,
        "tolerance": 0.04
    }
}


def bayesian_driver(search, calculator, **settings):
    driver = AutoCDAM()
    driver.calculator = calculator
    driver.update_settings({"search_method": "bayesian", **search, **settings})
    return driver


def fit(search,
        targets=TARGETS,
        calculator=None,
        parameters=PARAMETERS,
        **settings):
    calculator = calculator or Surface()
    driver = bayesian_driver(search, calculator, **settings)
    result = driver.fit(CDAMCandidate("Al", 3),
                        targets,
                        parameters=parameters)
    return driver, result, calculator


def _coupled_surface(sigma, q):
    """Nonlinear, coupled synthetic IOD and CN with a known optimum."""
    iod = 1.30 + 0.32 * sigma - 0.45 * (q - 0.45)**2 + 0.05 * math.sin(
        4.0 * sigma)
    cn = 4.2 + 5.0 * q + 0.8 * (sigma - 2.0)**2 - 0.6 * q * (sigma - 1.8)
    return {"IOD": iod, "CN": cn}


class _SyntheticCalculator:
    """Counts calls and returns the coupled surface instead of running MD."""

    def __init__(self):
        self.calls = 0

    def calculate(self, batch, parameters, targets, **kwargs):
        self.calls += 1
        values = _coupled_surface(parameters["metal_lj_sigma"],
                                  parameters["dummy_charge"])
        return CalculationResult(
            values={t.name: values[t.name] for t in targets})


def run(strategy, seed, budget):
    """One fit of the coupled surface (optimum sigma 1.62, q 0.61) with ``budget`` calls."""
    reference = _coupled_surface(1.62, 0.61)
    targets = {
        "IOD": {
            "value": reference["IOD"],
            "tolerance": 0.01,
            "role": "required"
        },
        "CN": {
            "value": reference["CN"],
            "tolerance": 0.04
        }
    }
    if strategy == "bayesian":
        search = {
            "search_method": "bayesian",
            "mode": "constrained",
            "max_evaluations": budget,
            "seed": seed
        }
    else:
        search = {
            "max_iterations": budget // 5,
            "samples_per_stage": 5,
            "seed": seed
        }
    parameters = {"sigma": (1.4, 2.6), "dummy_charge": (0.2, 0.7)}
    calculator = _SyntheticCalculator()
    driver = AutoCDAM()
    driver.calculator = calculator
    driver.update_settings(search)
    result = driver.fit(CDAMCandidate("Al", 3),
                        targets,
                        parameters=parameters)
    return {"solved": bool(result.converged), "evaluations": calculator.calls}


class TestAutoCdamBayesianSearch:

    def test_weighted_search_reaches_the_known_optimum(self):
        _, result, calculator = fit({
            "mode": "weighted",
            "max_evaluations": 25,
            "seed": 3
        })
        assert result.converged and result.diagnostics[
            "stop_reason"] == "threshold"
        assert result.parameters["metal_lj_sigma"] == pytest.approx(2.31,
                                                                    abs=0.03)
        assert result.parameters["dummy_charge"] == pytest.approx(0.29,
                                                                  abs=0.02)
        assert len(calculator.points) <= 25
        assert all(
            math.isclose(p, 0.9)
            for p in [result.parameters["metal_dummy_distance"]])
        assert set(result.diagnostics["surrogate_at_best"]) == {"IOD", "CN"}

    def test_same_seed_reproduces_and_another_seed_differs(self):
        first = fit({"max_evaluations": 10, "seed": 5})[2].points
        again = fit({"max_evaluations": 10, "seed": 5})[2].points
        other = fit({"max_evaluations": 10, "seed": 6})[2].points
        assert first == again and first != other

    def test_constrained_mode_finds_the_feasible_region(self):
        targets = {
            "IOD": {
                "value": 2.055,
                "tolerance": 0.01,
                "role": "required"
            },
            "CN": {
                "value": 5.16,
                "tolerance": 0.01,
                "weight": 1.0
            }
        }
        _, result, _ = fit(
            {
                "mode": "constrained",
                "max_evaluations": 20,
                "seed": 1
            }, targets)
        assert result.diagnostics["feasible"]
        assert abs(result.diagnostics["best_properties"]["IOD"] - 2.055) <= 0.01

    def test_unreachable_constraint_is_never_reported_as_converged(self):
        targets = {
            "IOD": {
                "value": 9.0,
                "tolerance": 0.01,
                "role": "required"
            },
            "CN": {
                "value": 5.16,
                "tolerance": 0.04
            }
        }
        _, result, calculator = fit(
            {
                "mode": "constrained",
                "max_evaluations": 9,
                "seed": 2
            }, targets)
        assert not result.converged and not result.diagnostics["feasible"]
        assert result.diagnostics["stop_reason"] in ("max_evaluations",
                                                     "acquisition_exhausted")
        assert len(calculator.points) <= 9

    def test_pareto_mode_reports_a_nondominated_tradeoff(self):
        targets = {
            "IOD": {
                "value": 1.80,
                "tolerance": 0.01
            },
            "D": {
                "value": 2.0,
                "tolerance": 0.01
            }
        }
        search = {"mode": "pareto", "max_evaluations": 14, "seed": 4}
        _, result, _ = fit(search, targets, parameters={"sigma": (1.4, 2.6)})
        front = result.diagnostics["pareto_front"]
        assert len(front) >= 2
        errors = [
            np.array(list(point["squared_errors"].values())) for point in front
        ]
        for i, a in enumerate(errors):
            assert not any(
                np.all(b <= a) and np.any(b < a)
                for j, b in enumerate(errors)
                if j != i)

    def test_batches_are_distinct_and_never_repeat_points(self):
        driver, result, calculator = fit({
            "batch_size": 3,
            "max_evaluations": 14,
            "seed": 7,
            "initial_samples": 5
        })
        phases = [
            entry for entry in driver.fitting_driver.history
            if entry["phase"] == "acquisition"
        ]
        assert phases and all(entry["n_calculations"] <= 3 for entry in phases)
        assert any(entry["n_calculations"] == 3
                   for entry in phases) or result.converged
        assert len(set(calculator.points)) == len(calculator.points)

    def test_resume_continues_exactly_like_an_uninterrupted_run(self, tmp_path):
        search = {"max_evaluations": 9, "seed": 8, "initial_samples": 5}
        targets = {
            "IOD": {
                "value": 2.055,
                "tolerance": 1e-6
            },
            "CN": {
                "value": 5.16,
                "tolerance": 1e-6
            }
        }
        full = fit({**search, "max_evaluations": 12}, targets)[2].points
        filename = str(tmp_path / "bo")
        first = fit(search, targets, filename=filename,
                    cache=False)[2].points
        assert first == full[:9]
        resumed_calc = Surface()
        longer_search = {**search, "max_evaluations": 12}
        driver = bayesian_driver(longer_search,
                                 resumed_calc,
                                 filename=filename,
                                 cache=False,
                                 restart=True)
        with pytest.raises(
                ValueError,
                match="differ"):  # changed budget is a different study
            driver.fit(CDAMCandidate("Al", 3), targets, parameters=PARAMETERS)
        driver = bayesian_driver(search,
                                 resumed_calc,
                                 filename=filename,
                                 cache=False,
                                 restart=True)
        result = driver.fit(CDAMCandidate("Al", 3),
                            targets,
                            parameters=PARAMETERS)
        assert resumed_calc.points == [] and result.diagnostics[
            "n_evaluations"] == 9

    def test_failed_calculations_stay_visible_and_do_not_train(self):
        search = {"max_evaluations": 15, "seed": 9}
        _, result, calculator = fit(search,
                                    calculator=Surface(fail_above_sigma=2.45))
        assert result.diagnostics["quality_counts"]["rejected"] >= 1
        assert result.parameters["metal_lj_sigma"] <= 2.45
        assert result.converged

    def test_noise_aware_gp_learns_noise_and_still_finds_the_region(self):
        driver, result, _ = fit(
            {
                "max_evaluations": 25,
                "seed": 10,
                "noise": True
            }, {
                "IOD": {
                    "value": 2.055,
                    "tolerance": 0.03
                },
                "CN": {
                    "value": 5.16,
                    "tolerance": 0.12
                }
            },
            calculator=Surface(noise=0.005, seed=1))
        from veloxchem.autocdam.search import _learned_noise
        model = driver.fitting_driver.surrogate.models["IOD"][0]
        assert _learned_noise(model) > 0
        assert result.converged

    def test_negligible_acquisition_stops_the_search_honestly(self):
        targets = {
            "IOD": {
                "value": 9.0,
                "tolerance": 0.01
            },
            "CN": {
                "value": 5.16,
                "tolerance": 0.04
            }
        }
        _, result, calculator = fit(
            {
                "max_evaluations": 30,
                "seed": 11,
                "initial_samples": 5,
                "acquisition_tolerance": 1e9,
                "patience": 1
            }, targets)
        assert result.diagnostics["stop_reason"] == "acquisition_exhausted"
        assert not result.converged and len(calculator.points) == 6

    def test_discrete_geometries_are_searched_separately(self, tmp_path):

        class GeometrySurface(Surface):

            def calculate(self, batch, parameters, targets, **kwargs):
                result = super().calculate(batch, parameters, targets)
                if batch.metadata["geometry"] == "tetrahedral":
                    result.values[
                        "IOD"] += 0.4  # tetrahedral cannot reach the target
                return result

        driver, result, _ = fit(
            {
                "max_evaluations": 12,
                "seed": 12,
                "geometries": ["tetrahedral", "octahedral"]
            },
            calculator=GeometrySurface(),
            filename=str(tmp_path / "run"))
        assert result.diagnostics["geometry"] == "octahedral"
        assert set(result.diagnostics["geometry_results"]) == {
            "tetrahedral", "octahedral"
        }
        assert driver.candidate.geometry == "octahedral"
        # One checkpoint group per geometry in the same run file.
        assert result.checkpoint_group == "autocdam/checkpoint_octahedral"
        for geometry in ("tetrahedral", "octahedral"):
            assert read_h5_group(tmp_path / "run",
                                 f"autocdam/checkpoint_{geometry}")["samples"]
        assert read_h5_group(tmp_path / "run", "autocdam/checkpoint") is None
        rows = samples_table(tmp_path / "run.h5")
        assert rows == samples_table(result)
        assert len(rows) == len(driver.fitting_driver.table.samples)

    @pytest.mark.parametrize("search, message", [
        ({
            "mode": "greedy"
        }, "mode"),
        ({
            "initial_samples": 40,
            "max_evaluations": 10
        }, "exceed"),
        ({
            "initial_samples": 1
        }, "at least 2"),
        ({
            "geometries": ["cubic"]
        }, "geometry"),
    ])
    def test_invalid_search_settings_are_rejected(self, search, message):
        with pytest.raises(VeloxChemError, match=message):
            fit(search)

    def test_constrained_mode_needs_a_required_target(self):
        with pytest.raises(VeloxChemError, match="required"):
            fit({"mode": "constrained"},
                {"IOD": {
                    "value": 2.0,
                    "tolerance": 0.01
                }})

    def test_unknown_search_settings_are_rejected(self):
        with pytest.raises(VeloxChemError, match="unknown settings"):
            AutoCDAM().update_settings({"acquisition": "ucb"})

    def test_unknown_search_method_is_rejected(self):
        driver = AutoCDAM()
        driver.update_settings({"search_method": "matrix"})
        with pytest.raises(VeloxChemError, match="search_method"):
            driver.fit(CDAMCandidate("Al", 3), TARGETS, parameters=PARAMETERS)

    def test_vector_targets_are_left_to_the_staged_search(self):
        with pytest.raises(VeloxChemError, match="scalar"):
            fit({"max_evaluations": 6},
                {"IOD": {
                    "value": [2.0, 2.1],
                    "tolerance": 0.01
                }})

    def test_bayesian_search_beats_staged_search_on_equal_budget(self):
        """Coupled synthetic surface; deterministic seeds, equal 20-evaluation budget."""
        bayes = [run("bayesian", seed, 20) for seed in range(3)]
        staged = [run("staged", seed, 20) for seed in range(3)]
        assert sum(r["solved"] for r in bayes) >= sum(
            r["solved"] for r in staged)
        assert sum(r["evaluations"] for r in bayes) < sum(
            r["evaluations"] for r in staged)

    @pytest.mark.parametrize("metric", ["absolute", "squared"])
    def test_acquisition_scores_match_the_evaluator(self, metric):
        from veloxchem.autocdam import Evaluator, FitTarget

        targets = [
            FitTarget("IOD", 1.9, 0.05, metric=metric, role="required"),
            FitTarget("CN", 6.0, 0.4, metric=metric, scale=0.2)
        ]
        draws = np.array([[1.93, 6.1], [1.70, 6.9], [2.0, 5.5]])
        terms, within = Evaluator.score_draws(draws, targets)
        for row, term, ok in zip(draws, terms, within):
            evaluation = Evaluator().evaluate(
                CalculationResult(values={
                    "IOD": row[0],
                    "CN": row[1]
                }), targets)
            assert term == pytest.approx(
                [evaluation.errors["IOD"]**2, evaluation.errors["CN"]**2])
            assert list(ok) == [
                evaluation.target_results[t.name]["within_tolerance"]
                for t in targets
            ]
