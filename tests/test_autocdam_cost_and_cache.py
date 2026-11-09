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
"""Step 7: property caching, cost tracking, screening, and failure visibility."""
import h5py
import numpy as np
import pytest

from veloxchem.autocdam import (AutoCDAM, CalculationResult, CDAMCandidate,
                                FitResult, FitTarget, ParameterSet,
                                ParameterSpec)
from veloxchem.autocdam.cache import CACHE_GROUP, CachedCalculator, load_cache
from veloxchem.autocdam.models import RESULT_GROUP, read_h5_group
from veloxchem.autocdam.bayesian import BayesianFittingDriver
from veloxchem.autocdam.search import ParameterPropertySample
from veloxchem.errorhandler import VeloxChemError

PARAMETERS = {"sigma": (1.4, 2.6), "dummy_charge": (0.2, 0.7)}


class Counting:
    """IOD/CN surface; records each call's targets and options."""

    def __init__(self, fail=False, interrupt_at=None):
        self.calls, self.fail, self.interrupt_at = [], fail, interrupt_at

    def calculate(self, batch, parameters, targets, **kwargs):
        self.calls.append(([t.name for t in targets], dict(kwargs)))
        if self.interrupt_at is not None and len(
                self.calls) >= self.interrupt_at:
            raise KeyboardInterrupt
        if self.fail:
            return CalculationResult(failed=True, failure_reason="diverged")
        sigma = batch.system.nonbonded_parameters["metal_lj_sigma_angstrom"]
        q = batch.system.charges[1]
        offset = 0.05 if kwargs.get("fidelity") == "low" else 0.0
        values = {
            "IOD": 0.9 + 0.5 * sigma + offset,
            "CN": 4.0 + 4.0 * q - 3 * offset,
            "HFE": -4500.0
        }
        metadata = {
            "HFE": {
                "quality_status": "unreviewed"
            }
        } if any(t.name == "HFE" for t in targets) else {}
        return CalculationResult(
            values={t.name: values[t.name] for t in targets}, metadata=metadata)


def batch(sigma=2.0, q=0.5):
    return CDAMCandidate("Al", 3, metal_lj_sigma=sigma,
                         dummy_charge=q).to_generated_batch()


IOD = FitTarget("IOD", 1.9, 0.03)
CN = FitTarget("CN", 6.0, 0.2)


def run(search, calculator, targets=None, **settings):
    fit = AutoCDAM()
    fit.calculator = calculator
    fit.update_settings({"search_method": "bayesian", **search, **settings})
    if targets is None:
        targets = {
            "IOD": {
                "value": 2.055,
                "tolerance": 1e-4
            },
            "CN": {
                "value": 5.16,
                "tolerance": 1e-4
            }
        }
    return fit, fit.fit(CDAMCandidate("Al", 3),
                        targets,
                        parameters=PARAMETERS)


class TestAutoCdamCostAndCache:

    def test_identical_candidates_reuse_properties_even_when_objectives_change(
            self):
        inner = Counting()
        cached = CachedCalculator(inner)
        first = cached.calculate(batch(), ParameterSet(), [IOD, CN])
        changed = [
            FitTarget("IOD", 2.5, 0.5, weight=9.0, role="required"),
            FitTarget("CN", 5.0, 1.0, role="monitor")
        ]
        again = cached.calculate(batch(), ParameterSet(), changed)
        assert len(inner.calls) == 1
        assert again.values == first.values and again.metadata["IOD"]["cached"]
        assert cached.stats["cache_hits"] == 2

    @pytest.mark.parametrize("change", ["parameter", "option", "unit"])
    def test_a_different_candidate_option_or_unit_misses_the_cache(
            self, change):
        inner = Counting()
        cached = CachedCalculator(inner)
        cached.calculate(batch(),
                         ParameterSet(), [IOD],
                         md_options={"steps": 10})
        request = {
            "batch": batch(),
            "targets": [IOD],
            "options": {
                "md_options": {
                    "steps": 10
                }
            }
        }
        if change == "parameter":
            request["batch"] = batch(sigma=2.1)
        elif change == "option":
            request["options"] = {"md_options": {"steps": 20}}
        else:
            request["targets"] = [FitTarget("IOD", 1.9, 0.03, unit="nm")]
        cached.calculate(request["batch"], ParameterSet(), request["targets"],
                         **request["options"])
        assert len(inner.calls) == 2

    def test_cache_key_is_readable(self):
        key = CachedCalculator.key(batch(), IOD, {"work_dir": "/a"})
        assert key.startswith("IOD [")
        assert '"metal_lj_sigma_angstrom": 2.0' in key and "work_dir" not in key

    def test_work_dir_does_not_change_the_identity(self):
        inner = Counting()
        cached = CachedCalculator(inner)
        cached.calculate(batch(), ParameterSet(), [IOD], work_dir="/a")
        cached.calculate(batch(), ParameterSet(), [IOD], work_dir="/b")
        assert len(inner.calls) == 1

    def test_only_missing_properties_are_calculated(self):
        inner = Counting()
        cached = CachedCalculator(inner)
        cached.calculate(batch(), ParameterSet(), [IOD])
        result = cached.calculate(batch(), ParameterSet(), [IOD, CN])
        assert inner.calls[-1][0] == ["CN"] and set(
            result.values) == {"IOD", "CN"}

    def test_failures_are_never_cached(self):
        inner = Counting(fail=True)
        cached = CachedCalculator(inner)
        for _ in range(2):
            assert cached.calculate(batch(), ParameterSet(), [IOD]).failed
        assert len(
            inner.calls) == 2 and cached.stats["failed_calculations"] == 2

    def test_persistent_cache_lives_in_the_run_file_and_survives_a_new_process(
            self, tmp_path):
        path = tmp_path / "run"
        first = CachedCalculator(Counting(), path)
        first.calculate(batch(), ParameterSet(), [IOD, CN])
        assert first.path == tmp_path / "run.h5"
        with h5py.File(first.path, "r") as h5f:
            assert CACHE_GROUP in h5f
        assert len(load_cache(path)) == 2
        inner = Counting()
        second = CachedCalculator(inner, path)
        result = second.calculate(batch(), ParameterSet(), [IOD, CN])
        assert inner.calls == [] and result.values["IOD"] == pytest.approx(1.9)
        second.calculate(batch(sigma=2.1), ParameterSet(), [IOD])
        assert len(load_cache(path)) == 3  # old entries + the new one

    def test_failed_calculations_are_not_written_to_the_run_file(self, tmp_path):
        path = tmp_path / "run"
        cached = CachedCalculator(Counting(fail=True), path)
        assert cached.calculate(batch(), ParameterSet(), [IOD]).failed
        assert load_cache(path) == {}
        assert CachedCalculator(Counting(), path).entries == {}

    def test_missing_run_file_gives_an_empty_cache(self, tmp_path):
        assert load_cache(tmp_path / "absent") == {}

    def test_new_weights_reuse_earlier_calculations_through_a_persistent_cache(
            self, tmp_path):
        filename = str(tmp_path / "properties")
        search = {"max_evaluations": 8, "initial_samples": 6, "seed": 3}
        first = Counting()
        run(search, first, cache=True, filename=filename)
        assert len(load_cache(filename)) == 2 * len(first.calls)
        reweighted = Counting()
        _, result = run(search,
                        reweighted,
                        cache=True,
                        filename=filename,
                        targets={
                            "IOD": {
                                "value": 2.055,
                                "tolerance": 1e-4,
                                "weight": 5.0
                            },
                            "CN": {
                                "value": 5.16,
                                "tolerance": 1e-4,
                                "weight": 0.2
                            }
                        })
        assert len(reweighted.calls) <= 2  # the six design points are reused
        assert result.diagnostics["cost"]["cache_hits"] >= 12

    def test_cost_report_counts_calculations_and_wall_time(self):
        calculator = Counting()
        _, result = run({
            "max_evaluations": 7,
            "initial_samples": 5,
            "seed": 4
        }, calculator)
        cost = result.diagnostics["cost"]
        assert cost["calculations"] == len(
            calculator.calls) == cost["n_candidates"] == 7
        assert cost["wall_time_s"] >= 0 and "default" in cost["per_system"]
        assert cost["per_system"]["default"]["calculations"] == 7

    def test_interrupted_search_resumes_from_its_checkpoint_and_cache(
            self, tmp_path):
        search = {"max_evaluations": 10, "initial_samples": 5, "seed": 5}
        _, reference = run(search, Counting())
        filename = str(tmp_path / "run")
        with pytest.raises(KeyboardInterrupt):
            run(search,
                Counting(interrupt_at=8),
                filename=filename,
                cache=True)
        assert len(load_cache(filename)) == 2 * 7  # the finished calculations
        assert read_h5_group(filename, RESULT_GROUP) is None
        resumed = Counting()
        targets = {
            "IOD": {
                "value": 2.055,
                "tolerance": 1e-4
            },
            "CN": {
                "value": 5.16,
                "tolerance": 1e-4
            }
        }
        fit = AutoCDAM()
        fit.calculator = resumed
        fit.update_settings({
            "search_method": "bayesian",
            **search,
            "filename": filename,
            "cache": True,
            "restart": True,
        })
        result = fit.fit(CDAMCandidate("Al", 3), targets, parameters=PARAMETERS)
        assert result.diagnostics["n_evaluations"] == 10
        assert result.parameters.to_dict() == reference.parameters.to_dict()
        assert len(resumed.calls) <= 3  # only the unfinished steps
        assert FitResult.load(filename).parameters.to_dict(
        ) == result.parameters.to_dict()

    def test_failures_are_listed_with_reasons_in_history_and_result(self):

        class Flaky(Counting):

            def calculate(self, batch, parameters, targets, **kwargs):
                if parameters["metal_lj_sigma"] > 2.3:
                    return CalculationResult(failed=True,
                                             failure_reason="box collapsed")
                return super().calculate(batch, parameters, targets, **kwargs)

        fit, result = run(
            {
                "max_evaluations": 9,
                "initial_samples": 6,
                "seed": 6
            }, Flaky())
        reasons = {f["reason"] for f in result.diagnostics["failures"]}
        assert reasons == {"box collapsed"}
        logged = [
            f for entry in fit.fitting_driver.history for f in entry["failures"]
        ]
        assert len(logged) == len(result.diagnostics["failures"])

    def test_screening_ranks_candidates_but_never_trains_the_surrogate(self):
        calculator = Counting()
        search = {
            "max_evaluations": 9,
            "initial_samples": 5,
            "seed": 7,
            "batch_size": 2,
            "screening": {
                "calculator_options": {
                    "fidelity": "low"
                },
                "oversample": 3
            }
        }
        fit, result = run(search, calculator)
        full = [c for c in calculator.calls if c[1].get("fidelity") != "low"]
        low = [c for c in calculator.calls if c[1].get("fidelity") == "low"]
        driver = fit.fitting_driver
        assert len(full) == len(driver.table.samples) == 9
        acquisition_steps = [
            e for e in driver.history if e["phase"] == "acquisition"
        ]
        assert len(low) == sum(
            len(e["screening"])
            for e in acquisition_steps) == 6 * len(acquisition_steps)
        screening = result.diagnostics["screening"]
        assert screening["used_in_surrogate"] is False
        assert screening["agreement"]["IOD"]["n_pairs"] == 4
        assert screening["agreement"]["IOD"][
            "mean_difference"] == pytest.approx(-0.05)

    @pytest.mark.parametrize("screening, message", [
        ({
            "oversample": 3
        }, "cheaper"),
        ({
            "calculator_options": {
                "x": 1
            },
            "oversample": 1
        }, "oversample"),
        ({
            "calculator_options": {
                "x": 1
            },
            "fidelity": "low"
        }, "unsupported"),
    ])
    def test_invalid_screening_settings(self, screening, message):
        with pytest.raises(VeloxChemError, match=message):
            run({"max_evaluations": 6, "screening": screening}, Counting())

    def test_cost_model_predicts_expensive_regions(self):
        specs = [ParameterSpec("metal_lj_sigma", 2.0, (1.4, 2.6))]
        search = {"max_evaluations": 10, "cost_aware": True}
        driver = BayesianFittingDriver(specs, [IOD],
                                       Counting(),
                                       None,
                                       lambda *a: [],
                                       search=search)
        for i, sigma in enumerate(np.linspace(1.4, 2.6, 6)):
            driver.table.add(
                ParameterPropertySample(f"s{i}",
                                        ParameterSet({"metal_lj_sigma": sigma}),
                                        CalculationResult(values={"IOD": 1.9}),
                                        "usable",
                                        runtime_seconds=float(
                                            np.exp(3 * (sigma - 1.4)))))
        cheap, expensive = driver._cost_factor(np.array([[0.05], [0.95]]))
        assert expensive > 5 * cheap

    def test_unreviewed_hfe_never_trains_or_wins(self):
        calculator = Counting()
        search = {"max_evaluations": 6, "initial_samples": 5, "seed": 8}
        targets = {"HFE": {"value": -4500.0, "tolerance": 50.0}}
        _, result = run(search, calculator, targets=targets)
        assert result.diagnostics["quality_counts"]["unreviewed"] == 6
        assert not result.converged and result.status == "failed"

    def test_cached_default_calculator_still_records_the_md_protocol(self):
        from veloxchem.autocdam.calculator import Calculator
        from veloxchem.autocdam.md import md_settings
        fit = AutoCDAM()
        fit.calculator = Calculator()
        fit.setup(CDAMCandidate("Al", 3),
                  {"IOD": {
                      "value": 2.055,
                      "tolerance": 0.01
                  }},
                  parameters=PARAMETERS)
        assert isinstance(fit.calculator, CachedCalculator)
        assert fit.calculator_options["md_options"] == md_settings(None)
