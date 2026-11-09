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
"""Step 8: fit results, diagnostics, export on request, and plotting."""
import json
from pathlib import Path

import pytest
import veloxchem as vlx

from veloxchem.autocdam import AutoCDAM, CalculationResult, CDAMCandidate, FitResult, FrameworkSystem, plot_fit

STRUCTURES = Path(__file__).parent / "data"
CHARGES = json.loads(
    (STRUCTURES /
     "autocdam_al_node_fragment_charges.json").read_text())["charges_e"]


class Surface:

    def calculate(self, batch, parameters, targets, **kwargs):
        sigma = batch.system.nonbonded_parameters["metal_lj_sigma_angstrom"]
        q = batch.system.charges[1]
        if sigma > 2.5:
            return CalculationResult(failed=True, failure_reason="unstable")
        values = {
            "IOD": 0.9 + 0.5 * sigma,
            "CN": 4.0 + 4.0 * q,
            "SHELL": 0.02,
            "D": 3.0 - 0.6 * sigma
        }
        return CalculationResult(
            values={t.name: values[t.name] for t in targets})


def fit(targets=None, search=None, parameters=None, **kwargs):
    parameters = parameters or {
        "sigma": (1.4, 2.6),
        "dummy_charge": (0.2, 0.7),
        "metal_dummy_distance": 0.9
    }
    targets = targets or {
        "IOD": {
            "value": 2.055,
            "tolerance": 0.01,
            "role": "required"
        },
        "CN": {
            "value": 5.16,
            "tolerance": 0.04,
            "weight": 2.0
        },
        "SHELL": {
            "value": 0.0,
            "tolerance": 0.1,
            "role": "monitor"
        }
    }
    search = search or {"max_evaluations": 14, "seed": 1}
    driver = AutoCDAM(candidate=CDAMCandidate("Al", 3),
                      parameters=parameters,
                      fit_targets=targets,
                      calculator=Surface(),
                      search=search,
                      **kwargs)
    return driver, driver.fit()


def result_parameters():
    from veloxchem.autocdam import ParameterSet
    return ParameterSet({"metal_lj_sigma": 2.0})


class TestAutoCdamResultsDiagnostics:

    def test_result_reports_definitions_targets_constraints_and_status(self):
        _, result = fit()
        d = result.diagnostics
        definitions = {p["name"]: p for p in d["parameter_definitions"]}
        assert definitions["metal_lj_sigma"]["adjustable"] and definitions[
            "metal_lj_sigma"]["bounds"] == [1.4, 2.6]
        assert definitions["metal_dummy_distance"] == {
            "name": "metal_dummy_distance",
            "adjustable": False,
            "bounds": None,
            "fixed_value": 0.9,
            "unit": None
        }
        roles = {o["name"]: o["role"] for o in d["objectives"]}
        assert roles == {
            "IOD": "required",
            "CN": "weighted",
            "SHELL": "monitor"
        }
        assert d["target_metal"] == "Al" and d["feasible"]
        assert set(d["target_results"]) == {"IOD", "CN"}
        assert "SHELL" in d["monitor"]
        assert d["optimization_converged"] is True and d[
            "scientifically_validated"] is False
        assert "not" in d["physical_validation"]
        rows = d["evaluations"]
        assert len(rows) == d["n_evaluations"] and {
            "parameters", "values", "feasible"
        } <= set(rows[0])
        assert d["surrogate_at_best"]["IOD"]["std"] >= 0
        assert d["cost"]["calculations"] >= len(rows)

    def test_failures_are_visible_in_the_evaluation_record(self):
        _, result = fit(search={
            "max_evaluations": 12,
            "seed": 2,
            "initial_samples": 10
        })
        failed = [r for r in result.diagnostics["evaluations"] if r["failed"]]
        assert failed and all(
            r["failure_reason"] == "unstable" and r["objective"] is None
            for r in failed)

    def test_solution_at_a_bound_is_flagged(self):
        targets = {"IOD": {"value": 1.60, "tolerance": 0.001}}
        search = {"max_evaluations": 12, "seed": 3}
        _, result = fit(targets=targets,
                        parameters={"sigma": (1.4, 2.4)},
                        search=search)
        assert result.parameters["metal_lj_sigma"] == pytest.approx(1.4,
                                                                    abs=0.03)
        assert any(
            "lower bound" in w for w in result.diagnostics["bound_warnings"])
        assert any("lower bound" in w for w in result.warnings)

    def test_summary_distinguishes_convergence_from_validation(self):
        _, result = fit()
        text = result.summary()
        assert "optimization converged: True" in text and "Physical validation: not performed" in text
        assert "metal_dummy_distance = 0.9 (fixed)" in text and "[required]" in text
        assert "SHELL [monitor]" in text

    def test_nothing_is_written_unless_requested(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        fit()
        assert list(tmp_path.iterdir()) == []

    def test_result_json_round_trip_and_plot_from_file(self, tmp_path):
        _, result = fit()
        path = result.save(tmp_path / "result.json")
        data = json.loads(path.read_text())
        assert data["diagnostics"]["evaluations"] and data["converged"]
        figure = plot_fit(path, output_path=tmp_path / "from_file.png")
        assert len(figure.axes) == 4 and (tmp_path / "from_file.png").exists()

    def test_diagnostic_plot_has_convergence_residual_bound_and_history_panels(
            self, tmp_path):
        _, result = fit()
        output = tmp_path / "diagnostics.png"
        figure = result.plot(output_path=output)
        assert output.exists()
        assert [axis.get_title(loc="left") for axis in figure.axes] == [
            "Optimization convergence",
            "Target residuals at best (orange: required)",
            "Parameter bounds (orange: at a bound)", "Recorded objective"
        ]

    def test_pareto_tradeoff_plot(self):
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
        search = {"mode": "pareto", "max_evaluations": 12, "seed": 4}
        _, result = fit(targets=targets,
                        parameters={"sigma": (1.4, 2.4)},
                        search=search)
        assert result.diagnostics["pareto_front"]
        assert result.plot().axes[3].get_title(loc="left") == "Trade-off"

    def test_legacy_results_without_evaluations_still_plot(self):
        legacy = FitResult(parameters=result_parameters(),
                           history=[{
                               "iteration": 1,
                               "objective": 3.0,
                               "errors": {
                                   "IOD": 2.0
                               }
                           }, {
                               "iteration": 2,
                               "objective": 1.0,
                               "errors": {
                                   "IOD": 0.8
                               }
                           }])
        assert plot_fit(legacy) is not None

    def test_mof_fit_reports_preserved_framework_and_exports_on_request(
            self, tmp_path):
        mol = vlx.Molecule.read_xyz_file(
            str(STRUCTURES / "autocdam_al_node_fragment.xyz"))
        mol.set_charge(-1)
        mof = FrameworkSystem(structure=mol,
                              charges=CHARGES,
                              target_metal="Al",
                              representation="cdam",
                              name="mof")
        parameters = {"sigma": (1.5, 2.5), "dummy_charge": (0.3, 0.7)}
        targets = {
            "al_o": {
                "property": "distance",
                "atoms": ["Al", "O3"],
                "value": 1.93,
                "tolerance": 0.05
            }
        }
        search = {"max_evaluations": 6, "initial_samples": 5, "seed": 5}
        driver = AutoCDAM(candidate=CDAMCandidate("Al", 3),
                          systems={"mof": mof},
                          parameters=parameters,
                          fit_targets=targets,
                          search=search)
        result = driver.fit()
        report = result.diagnostics["systems"]["mof"]
        assert report["non_target_preserved"] and "VeloxChem" in report[
            "baseline_source"]
        assert "neighbouring parameters preserved: True" in result.summary()
        written = driver.export_fitted_model(tmp_path / "model")
        names = {Path(p).name for p in written}
        assert {
            "fitted_cdam.xml", "fitted_cdam.pdb", "mof_system.xml",
            "mof_positions.xyz", "fitted_model.json"
        } <= names
        lines = (tmp_path / "model" /
                 "mof_positions.xyz").read_text().splitlines()
        assert int(lines[0]) == 19 + 6  # framework atoms plus six dummies
        provenance = json.loads(
            (tmp_path / "model" / "fitted_model.json").read_text())
        assert provenance["scientifically_validated"] is False
        assert provenance["systems"]["mof"][
            "non_target_parameters"] == "baseline values, unchanged"

    def test_export_requires_a_fit(self):
        driver = AutoCDAM(
            candidate=CDAMCandidate("Al", 3),
            parameters={"sigma": (1.4, 2.6)},
            calculator=Surface(),
            fit_targets={"IOD": {
                "value": 2.0,
                "tolerance": 0.01
            }})
        with pytest.raises(ValueError, match="fit"):
            driver.export_fitted_model("unused")
