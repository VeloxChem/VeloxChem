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

from veloxchem.autocdam import AutoCDAM, CDAMCandidate


class TestAutoCdamDriverSmoke:

    def test_autocdam_physical_fit_runs_end_to_end(self, tmp_path):
        import math
        import random
        from veloxchem import OutputStream

        random.seed(42)
        ostream = OutputStream(str(tmp_path / "fit.log"))
        try:
            cdam = AutoCDAM(ostream=ostream)
            cdam.update_settings({
                "search_method": "staged",
                "max_iterations": 1,
                "samples_per_stage": 2,
                "pilot": True,
                "work_dir": str(tmp_path),
                "checkpoint_file": str(tmp_path / "checkpoint.json"),
            })
            result = cdam.fit(CDAMCandidate("Zn", 2), {
                "IOD": {
                    "value": 2.09,
                    "tolerance": 0.03
                },
                "CN": {
                    "value": 6.0,
                    "tolerance": 0.25
                }
            })
            assert result.iterations == 1
            assert result.diagnostics["n_calculations"] == 2
            assert math.isfinite(result.objective)
            assert set(result.diagnostics["best_properties"]) == {"IOD", "CN"}
            best = next(
                sample for sample in cdam.fitting_driver.table.samples
                if sample.sample_id == result.diagnostics["best_sample_id"])
            assert best.parameters.values == result.parameters.values
            assert best.calculation.values == result.diagnostics[
                "best_properties"]
            assert best.calculation.raw_data["md_analysis"]["n_frames"] >= 2
            assert best.calculation.metadata["IOD"][
                "method"] == "OpenMM MD trajectory analysis"
            assert best.calculation.metadata["IOD"]["md_options"][
                "timestep_fs"] == 0.25
            assert best.calculation.metadata["IOD"]["pilot"] is True
            assert (tmp_path / "checkpoint.json").exists()
        finally:
            ostream.close()

    def test_custom_evaluator_receives_unmodified_calculator_result(
            self) -> None:
        from veloxchem.autocdam import CalculationResult, EvaluationResult, ParameterSet

        calculation = CalculationResult(
            values={"IOD": 2.1},
            metadata={"IOD": {
                "uncertainty": 0.05
            }},
            raw_data={"trajectory": "caller-owned"},
        )

        class SuppliedCalculator:

            def calculate(self, batch, parameters, targets):
                assert isinstance(parameters, ParameterSet)
                return calculation

        class SuppliedEvaluator:

            def evaluate(self, result, targets):
                assert result is calculation
                assert targets[0].name == "IOD"
                return EvaluationResult(objective=7.0,
                                        metrics={"custom_loss": 7.0})

        cdam = AutoCDAM()
        cdam.calculator = SuppliedCalculator()
        cdam.evaluator = SuppliedEvaluator()
        result = cdam.fit(
            CDAMCandidate("Zn", 2, metal_lj_sigma=2.0, metal_lj_epsilon=0.2),
            {"IOD": {
                "reference": 2.0,
                "tolerance": 0.1
            }},
            fit_parameters={
                "dummy_charge": {
                    "initial": 0.5,
                    "bounds": [0.0, 1.0]
                }
            })
        assert result.objective == 7.0
        assert result.metrics == {"custom_loss": 7.0}
