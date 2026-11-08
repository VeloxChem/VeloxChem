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

import math

from veloxchem.autocdam import CalculationResult, Evaluator, FitTarget


class TestAutoCdamEvaluatorInterface:

    def test_evaluator_handles_scalar_and_array_targets(self) -> None:
        evaluator = Evaluator()
        result = evaluator.evaluate(
            CalculationResult(values={
                "IOD": 2.1,
                "QM_ENERGY": [0.0, 1.4, 3.0]
            }),
            [
                FitTarget("IOD", reference=2.0, tolerance=0.1, weight=2.0),
                FitTarget("QM_ENERGY", reference=[0.0, 1.0, 3.0],
                          tolerance=0.2),
            ],
        )
        assert math.isclose(result.metrics["IOD"], 0.1)
        assert result.metrics["QM_ENERGY"] > 0
        assert result.objective > 0
        assert result.failed is False

    def test_failed_calculation_gets_finite_penalty_and_reason(self) -> None:
        result = Evaluator(failure_penalty=123.0).evaluate(
            CalculationResult(failed=True,
                              failure_reason="unstable simulation"),
            [FitTarget("HFE", reference=-10.0, tolerance=1.0)],
        )
        assert result.objective == 123.0
        assert math.isfinite(result.objective)
        assert result.failed is True
        assert result.failure_reason == "unstable simulation"
        assert "unstable simulation" in result.warnings[0]

    def test_nonfinite_prediction_is_not_a_successful_objective(self) -> None:
        for value in (float("nan"), float("inf"), -float("inf")):
            result = Evaluator(failure_penalty=123.0).evaluate(
                CalculationResult(values={"IOD": value}),
                [FitTarget("IOD", reference=2.0, tolerance=0.1)],
            )
            assert result.failed
            assert result.objective == 123.0
