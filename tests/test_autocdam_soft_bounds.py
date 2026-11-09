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
"""UFF-referenced soft bounds and their widening during the Bayesian search."""
import json

import pytest

from veloxchem.autocdam import AutoCDAM, CalculationResult, CDAMCandidate, ParameterSpec
from veloxchem.autocdam.models import HARD_LIMITS, reference_bounds, uff_lj


class Surface:
    """IOD rises with sigma; CN with dummy charge."""

    def calculate(self, batch, parameters, targets, **kwargs):
        sigma = batch.system.nonbonded_parameters["metal_lj_sigma_angstrom"]
        q = batch.system.charges[1]
        values = {"IOD": 0.9 + 0.5 * sigma, "CN": 4.0 + 4.0 * q}
        return CalculationResult(
            values={t.name: values[t.name] for t in targets})


def fit(parameters, iod, search=None, **settings):
    targets = {"IOD": {"value": iod, "tolerance": 0.002}}
    driver = AutoCDAM()
    driver.calculator = Surface()
    driver.update_settings({
        "search_method": "bayesian",
        "max_evaluations": 10,
        "seed": 1,
        **(search or {}),
        **settings
    })
    return driver, driver.fit(CDAMCandidate("Al", 3),
                              targets,
                              parameters=parameters)


class TestAutoCdamSoftBounds:

    def test_uff_reference_scales_bounds_by_element(self):
        sigma, epsilon = uff_lj("Al")
        assert sigma == pytest.approx(
            4.008, abs=1e-3) and epsilon == pytest.approx(0.505, abs=1e-3)
        bounds = reference_bounds("Al")
        assert bounds["metal_lj_sigma"]["bounds"] == pytest.approx(
            (0.3 * sigma, 1.1 * sigma))
        assert bounds["metal_lj_epsilon"]["hard_bounds"] == HARD_LIMITS[
            "metal_lj_epsilon"]
        assert reference_bounds("Zn")["metal_lj_sigma"]["bounds"][1] < bounds[
            "metal_lj_sigma"]["bounds"][1]

    def test_uff_keyword_gives_soft_bounds_and_explicit_bounds_stay_hard(self):
        parameters = {"sigma": "uff", "dummy_charge": (0.2, 0.7)}
        driver = AutoCDAM()
        driver.calculator = Surface()
        driver.setup(CDAMCandidate("Al", 3), {"IOD": {
            "value": 2.0
        }},
                     parameters=parameters)
        specs = {s.name: s for s in driver.parameter_specs}
        assert specs["metal_lj_sigma"].expandable and specs[
            "metal_lj_sigma"].hard_bounds
        assert "UFF Al" in specs["metal_lj_sigma"].metadata["bounds_source"]
        assert not specs["dummy_charge"].expandable

    def test_default_parameters_use_uff_soft_bounds(self):
        driver = AutoCDAM()
        driver.calculator = Surface()
        driver.setup(CDAMCandidate("Al", 3), {"IOD": {"value": 2.0}})
        names = {s.name for s in driver.parameter_specs}
        assert names == {"metal_lj_sigma", "metal_lj_epsilon"}
        assert all(s.expandable for s in driver.parameter_specs)

    @pytest.mark.parametrize("value, message",
                             [("gaff", "unknown bound keyword"),
                              ("uff", "no reference bounds")])
    def test_bad_reference_requests(self, value, message):
        parameters = {"sigma": (1.4, 2.6), "dummy_charge": value}
        with pytest.raises(ValueError, match=message):
            AutoCDAM().setup(CDAMCandidate("Al", 3), {"IOD": {
                "value": 2.0
            }},
                             parameters=parameters)

    def test_hard_bounds_must_contain_bounds(self):
        with pytest.raises(ValueError, match="inside hard_bounds"):
            ParameterSpec("metal_lj_sigma",
                          2.0, (1.0, 3.0),
                          hard_bounds=(1.5, 3.0))

    def test_solution_beyond_a_soft_bound_is_reached_by_widening(self):
        # IOD 1.40 needs sigma = 1.0, below the first soft bound [1.4, 2.6].
        soft = ParameterSpec("metal_lj_sigma",
                             2.0, (1.4, 2.6),
                             expandable=True,
                             hard_bounds=(0.5, 8.0))
        driver, result = fit({"sigma": soft}, iod=1.40)
        d = result.diagnostics
        assert d["bound_expansions"] and d["bound_expansions"][0][
            "side"] == "lower"
        assert d["final_bounds"]["metal_lj_sigma"][0] < 1.0
        assert result.parameters["metal_lj_sigma"] == pytest.approx(1.0,
                                                                    abs=0.01)
        assert result.converged and not d["bound_warnings"]
        definition = d["parameter_definitions"][0]
        assert definition["initial_bounds"] == [1.4, 2.6]
        assert definition["expandable"]
        n_events = len({c["step"] for c in d["bound_expansions"]})
        assert d["n_evaluations"] <= 10 + 3 * n_events
        assert "bound widened" in result.summary()

    def test_hard_bound_is_never_widened_and_is_flagged(self):
        _, result = fit({"sigma": (1.4, 2.6)}, iod=1.40)
        d = result.diagnostics
        assert d["bound_expansions"] == [] and d["final_bounds"][
            "metal_lj_sigma"] == [1.4, 2.6]
        assert any("lower bound" in w for w in d["bound_warnings"])

    def test_widening_stops_at_the_hard_limit_and_expansion_count(self):
        soft = ParameterSpec("metal_lj_sigma",
                             2.0, (1.4, 2.6),
                             expandable=True,
                             hard_bounds=(1.2, 8.0))
        _, result = fit({"sigma": soft},
                        iod=1.40,
                        search={"max_bound_expansions": 2})
        d = result.diagnostics
        assert d["final_bounds"]["metal_lj_sigma"][0] == pytest.approx(1.2)
        assert len({c["step"] for c in d["bound_expansions"]}) <= 2
        assert any(
            "soft bound already widened" in w for w in d["bound_warnings"])

    def test_design_points_near_a_bound_do_not_widen_it(self):
        soft = ParameterSpec("metal_lj_sigma",
                             2.0, (1.4, 2.6),
                             expandable=True,
                             hard_bounds=(0.5, 8.0))
        _, result = fit({"sigma": soft},
                        iod=1.40,
                        search={"initial_samples": 4})
        d = result.diagnostics
        assert d["bound_expansions"]
        assert all(c["step"] > 1 for c in d["bound_expansions"])

    def test_widened_search_resumes_with_its_bounds(self, tmp_path):
        soft = ParameterSpec("metal_lj_sigma",
                             2.0, (1.4, 2.6),
                             expandable=True,
                             hard_bounds=(0.5, 8.0))
        checkpoint = tmp_path / "bo.json"
        _, first = fit({"sigma": soft},
                       iod=1.40,
                       checkpoint_file=str(checkpoint))
        saved = json.loads(checkpoint.read_text())
        assert saved["history"][-1]["bounds"]["metal_lj_sigma"][0] < 1.4
        _, resumed = fit({"sigma": soft},
                         iod=1.40,
                         checkpoint_file=str(checkpoint),
                         restart=True)
        assert resumed.parameters.to_dict() == first.parameters.to_dict()
        assert resumed.diagnostics["final_bounds"] == first.diagnostics[
            "final_bounds"]
