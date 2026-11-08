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
"""Check physical calculation planning and explicit failure reporting."""
import pytest

from veloxchem.autocdam import Calculator, CDAMCandidate, Evaluator, FitTarget, ParameterSet


class TestAutoCdamCalculatorInterface:

    def test_calculation_plan_reuses_md_job_for_iod_and_cn(self):
        plan = Calculator().plan([
            FitTarget("IOD", reference=2.0, tolerance=0.1),
            FitTarget("CN", reference=6.0, tolerance=0.5),
        ])
        assert plan.jobs["md"] == ["IOD", "CN"]
        assert "OpenMM MD" in plan.methods["IOD"]

    def test_calculator_uses_injected_veloxchem_output_stream(self, tmp_path):
        from veloxchem import OutputStream

        path = tmp_path / "calculation.log"
        ostream = OutputStream(str(path))
        try:
            Calculator(ostream=ostream).plan(
                [FitTarget("IOD", reference=2.0, tolerance=0.1)])
            text = path.read_text()
            assert "AutoCDAM calculation plan" in text
            assert "OpenMM MD trajectory analysis" in text
        finally:
            ostream.close()

    @pytest.mark.parametrize("name", ["UNKNOWN", "iod"])
    def test_unsupported_target_fails_before_system_construction(
            self, name, monkeypatch):

        def unexpected_platform_probe(self):
            raise AssertionError(
                "Unsupported targets must fail before platform initialization")

        monkeypatch.setattr(Calculator, "_initialize_platform",
                            unexpected_platform_probe)
        with pytest.raises(ValueError, match="unsupported target"):
            Calculator().calculate(None, ParameterSet(),
                                   [FitTarget(name, 0, 1)])

    def test_missing_fep_api_fails_before_system_construction(
            self, monkeypatch):
        import veloxchem.autocdam.md as module

        monkeypatch.setattr(module, "SolvationFepDriver",
                            type("UnavailableFepDriver", (), {}))

        def unexpected_platform_probe(self):
            raise AssertionError(
                "FEP capability should be checked before initialization")

        monkeypatch.setattr(Calculator, "_initialize_platform",
                            unexpected_platform_probe)
        candidate = CDAMCandidate("Zn",
                                  2,
                                  metal_lj_sigma=2.0,
                                  metal_lj_epsilon=0.2)
        batch = candidate.to_generated_batch()
        with pytest.raises(RuntimeError, match="lacks skip_vacuum"):
            Calculator().calculate(batch, ParameterSet(),
                                   [FitTarget("HFE", -10.0, 1.0)])

    def test_failed_calculation_reports_failure_to_evaluator(self):
        target = FitTarget("IOD", reference=2.0, tolerance=0.1)
        result = Calculator()._generate_failed_result([target], ParameterSet())
        evaluation = Evaluator(failure_penalty=123.0).evaluate(result, [target])
        assert result.failed
        assert evaluation.failed
        assert evaluation.objective == 123.0
