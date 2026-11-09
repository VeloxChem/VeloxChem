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
"""Score calculated properties against scalar or vector reference targets.

Errors are divided by each target's normalization scale (its tolerance by
default) and combined into the normalized weighted loss
``sum_i w_i e_i^2 / sum_i w_i`` over required and weighted targets. Required
targets are also checked as hard constraints, independently of the loss, so a
good weighted score never hides a violated constraint. Monitored targets are
reported only. Failed or invalid calculations receive a finite penalty and
retain their failure reason. Raw calculated values are never modified.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np

from ..outputstream import OutputStream
from .models import CalculationResult, EvaluationResult, FitTarget, _is_sequence, default_ostream


class Evaluator:
    """
    Converts calculator values into a normalized scalar objective.

    The ``CalculationResult -> EvaluationResult`` handoff compares scalar and
    array targets, applies the target tolerance and weight, and assigns a large
    finite penalty when a calculation fails or cannot be compared.  It does
    not know whether values came from MD, FEP, QM, or another calculator.

    ``ostream`` uses the same output contract as the surrounding VeloxChem
    package, so warnings have the familiar VeloxChem formatting and are still
    retained in ``EvaluationResult.warnings`` for programmatic inspection.

    :param failure_penalty:
        The objective of a failed proposal (positive and finite).
    :param ostream:
        The output stream.
    """

    def __init__(
        self,
        failure_penalty: float = 1.0e6,
        ostream: OutputStream | None = None,
    ) -> None:
        """
        Initializes objective normalization and warning output.

        ``failure_penalty`` must be finite so an optimizer can continue after
        an unstable proposal while still recording the failure in its history.
        The same VeloxChem output stream used by the driver receives the human-
        readable warning.
        """

        if not math.isfinite(failure_penalty) or failure_penalty <= 0:
            raise ValueError("failure_penalty must be a positive finite number")
        self.failure_penalty = float(failure_penalty)
        self.ostream = default_ostream(ostream)

    def evaluate(
        self,
        calculation: CalculationResult,
        targets: Iterable[FitTarget],
    ) -> EvaluationResult:
        """
        Compares scalar and array targets using one shared objective formula.

        :param calculation:
            The calculation result.
        :param targets:
            The targets.

        :return:
            The evaluation result.
        """

        targets = list(targets)
        if calculation.failed:
            reason = calculation.failure_reason or "calculator reported failure"
            return self._failure(reason, calculation.warnings)

        weighted_sum = 0.0
        weight_total = 0.0
        metrics: dict[str, float] = {}
        errors: dict[str, float] = {}
        target_results: dict[str, dict[str, object]] = {}
        violations: dict[str, float] = {}

        for target in targets:
            if target.name not in calculation.values:
                if target.is_monitor:
                    # Monitors are usually calculated for the final candidate only.
                    continue
                reason = f"missing calculated value for target {target.name}"
                return self._failure(reason, calculation.warnings)
            prediction = calculation.values[target.name]
            try:
                raw_error, metric_name = self._error(prediction,
                                                     target.reference,
                                                     target.metric)
                if not math.isfinite(raw_error):
                    raise ValueError("calculated error is not finite")
            except (TypeError, ValueError) as exc:
                reason = f"target {target.name} could not be compared: {exc}"
                return self._failure(reason, calculation.warnings)
            normalized = raw_error / target.normalization_scale
            within_tolerance = raw_error <= target.tolerance
            if target.is_objective:
                weighted_sum += target.weight * normalized**2
                weight_total += target.weight
            if target.is_required and not within_tolerance:
                violations[target.name] = raw_error / target.tolerance
            metrics[target.name] = raw_error
            errors[target.name] = normalized
            target_results[target.name] = {
                "prediction": prediction,
                "reference": target.reference,
                "metric": metric_name,
                "raw_error": raw_error,
                "normalized_error": normalized,
                "within_tolerance": within_tolerance,
                "role": target.role,
                "kind": target.kind,
                "weight": target.weight,
            }

        objective = weighted_sum / weight_total if weight_total > 0 else 0.0
        return EvaluationResult(
            objective=float(objective),
            metrics=metrics,
            errors=errors,
            target_results=target_results,
            warnings=list(calculation.warnings),
            feasible=not violations,
            constraint_violations=violations,
        )

    def _failure(self, reason: str, warnings: list[str]) -> EvaluationResult:
        """Report one consistent penalty while preserving calculator diagnostics."""
        warning = f"WARNING: Evaluation failed: {reason} A finite penalty was assigned."
        self.ostream.print_warning(warning)
        self.ostream.flush()
        return EvaluationResult(
            objective=self.failure_penalty,
            failed=True,
            failure_reason=reason,
            warnings=[*warnings, warning],
            feasible=False,
        )

    def evaluate_predictions(self, values: dict,
                             targets: Iterable[FitTarget]) -> EvaluationResult:
        """
        Scores surrogate means for choosing where to calculate next.

        This score cannot stop the fit. Custom evaluators that require raw
        trajectories can provide their own prediction-scoring method.

        :param values:
            The predicted values by target name.
        :param targets:
            The targets.

        :return:
            The evaluation result.
        """

        return self.evaluate(CalculationResult(values=values), targets)

    @staticmethod
    def score_draws(
            draws: np.ndarray,
            targets: Sequence[FitTarget]) -> tuple[np.ndarray, np.ndarray]:
        """
        Vectorized ``evaluate`` for scalar targets, for many predicted draws at once.

        ``draws`` has shape (..., len(targets)), one column per target. Uses the
        same scalar metrics as ``_error``.

        :param draws:
            The predicted draws.
        :param targets:
            The scalar targets.

        :return:
            The squared normalized errors (the loss terms) and whether each raw
            error is within tolerance.
        """

        reference = np.array([float(t.reference) for t in targets])
        scale = np.array([t.normalization_scale for t in targets])
        tolerance = np.array([t.tolerance for t in targets])
        squared = np.array([t.metric.lower() == "squared" for t in targets])
        difference = draws - reference
        raw = np.where(squared, difference**2, np.abs(difference))
        return (raw / scale)**2, raw <= tolerance

    @staticmethod
    def _error(prediction, reference, metric: str) -> tuple[float, str]:
        """Raw error and the metric name actually used.

        Arrays use ``mae``, ``sse`` or (default) ``rmse``; scalars use
        ``squared`` or (default) ``absolute``. An unknown metric falls back to
        the default. Raises ValueError for arrays of unequal or zero length.
        """
        if _is_sequence(prediction) or _is_sequence(reference):
            predicted = [float(value) for value in prediction]
            expected = [float(value) for value in reference]
            if len(predicted) != len(expected) or not predicted:
                raise ValueError(
                    "array predictions and references must have equal non-zero length"
                )
            differences = [
                actual - target for actual, target in zip(predicted, expected)
            ]
            if metric.lower() == "mae":
                return sum(abs(value)
                           for value in differences) / len(differences), "mae"
            if metric.lower() == "sse":
                return sum(value**2 for value in differences), "sse"
            return math.sqrt(
                sum(value**2 for value in differences) /
                len(differences)), "rmse"
        difference = float(prediction) - float(reference)
        if metric.lower() == "squared":
            return difference**2, "squared"
        return abs(difference), "absolute"


def within_tolerance(evaluation: EvaluationResult,
                     targets: Iterable[FitTarget]) -> bool:
    """True when every required and weighted target is inside its tolerance.

    Monitored targets never gate convergence. Raw errors are compared with the
    tolerance directly, so a custom normalization ``scale`` cannot change it.
    """
    if evaluation is None or evaluation.failed:
        return False
    return all(target.name in evaluation.metrics and
               evaluation.metrics[target.name] <= target.tolerance
               for target in targets
               if target.is_objective)
