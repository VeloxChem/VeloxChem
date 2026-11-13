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
"""Parameter searches: observations, surrogates, shared bookkeeping, staged search.

A search proposes parameter points, calculates them into one observation
table, and judges them with the ``Evaluator``. Only calculated results can
stop a search; predictions only guide it. Sections, in order:

* sampling: ``ParameterPropertyTable`` records every attempt, failed ones
  included; ``default_quality_check`` decides which observations may train a
  surrogate or win a fit. Raw values are stored without target weights so a
  table can be rescored with different targets;
* surrogate: one Gaussian process per raw property;
* ``SearchDriver``: the table, measured best, stopping criterion, checkpoints
  and the final ``FitResult``;
* ``FittingDriver``, the staged (coarse-to-fine) search: each stage calculates
  ``samples_per_stage`` points (the surrogate's predicted best, Latin-hypercube
  samples of the current region and, after stage 1, some global samples), then
  the region shrinks by ``shrink_fraction`` around the best measured point.

``bayesian.py`` adds the Bayesian search on the same base class.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict, replace
from typing import Any, Iterable, Literal, Sequence, Callable
import json
import math
import time
from pathlib import Path

import importlib.util

import numpy as np
from mpi4py import MPI

from ..errorhandler import assert_msg_critical
from ..veloxchemlib import mpi_master
from .evaluator import within_tolerance
from .models import (CalculationResult, EvaluationResult, FitTarget,
                     GeneratedBatch, ParameterSet, ParameterSpec, all_finite,
                     FitResult, format_value, h5_filename, print_settings,
                     read_h5_group, write_h5_groups)

# ==========================================================================
# sampling
# ==========================================================================

QualityStatus = Literal["unreviewed", "usable", "rejected"]


@dataclass
class ParameterPropertySample:
    """One parameter point and its calculated properties, quality, and protocol."""

    sample_id: str
    parameters: ParameterSet
    calculation: CalculationResult
    quality_status: QualityStatus = "unreviewed"
    protocol_id: str | None = None
    runtime_seconds: float | None = None
    notes: list[str] = field(default_factory=list)


@dataclass
class ParameterPropertyTable:
    """
    Accumulated observations, including failed attempts, independent of scoring.

    :param samples:
        The samples.
    :param study_metadata:
        The metadata that defines the study.
    """

    samples: list[ParameterPropertySample] = field(default_factory=list)
    study_metadata: dict[str, Any] = field(default_factory=dict)

    def add(self, sample: ParameterPropertySample) -> None:
        """
        Appends a sample.

        Raises ``ValueError`` if its ``sample_id`` already exists.

        :param sample:
            The sample.
        """

        if any(row.sample_id == sample.sample_id for row in self.samples):
            raise ValueError(f"duplicate sample_id: {sample.sample_id}")
        self.samples.append(sample)

    def rescore(self, evaluator,
                targets: Iterable[FitTarget]) -> dict[str, EvaluationResult]:
        """
        Evaluates every usable sample against the targets.

        :param evaluator:
            The evaluator.
        :param targets:
            The targets.

        :return:
            The evaluations by sample id.
        """

        targets = list(targets)
        return {
            row.sample_id: evaluator.evaluate(row.calculation, targets)
            for row in self.samples
            if row.quality_status == "usable"
        }


def latin_hypercube(specs: Sequence[ParameterSpec], n_samples: int,
                    seed: int) -> list[ParameterSet]:
    """Seeded Latin-hypercube points over the bounds; fixed parameters keep ``initial``."""
    rng = np.random.default_rng(seed)
    points = [
        dict((s.name, s.initial) for s in specs) for _ in range(n_samples)
    ]
    for spec in specs:
        if not spec.optimizable:
            continue
        lo, hi = spec.bounds
        fractions = (rng.permutation(n_samples) +
                     rng.random(n_samples)) / n_samples
        for point, fraction in zip(points, fractions):
            point[spec.name] = float(lo + fraction * (hi - lo))
    return [ParameterSet(point) for point in points]


@dataclass
class CalculationInput:
    """Generated candidate and calculator options for one physical calculation."""

    batch: GeneratedBatch
    calculator_kwargs: dict[str, Any] = field(default_factory=dict)
    protocol_id: str | None = None


def default_quality_check(calculation: CalculationResult) -> QualityStatus:
    """Accept finite results and honor explicit calculator quality decisions.

    Warnings remain attached to the observation but do not exclude it from
    fitting. A calculator can mark an estimate ``unreviewed`` or ``rejected``;
    a study-specific check can apply additional sampling requirements.
    """
    values = calculation.values
    if calculation.failed or not values or not all(
            map(all_finite, values.values())):
        return "rejected"
    states = [
        v.get("quality_status")
        for v in calculation.metadata.values()
        if isinstance(v, dict)
    ]
    if "rejected" in states:
        return "rejected"
    if "unreviewed" in states:
        return "unreviewed"
    return "usable"


# ==========================================================================
# surrogate
# ==========================================================================


@dataclass(frozen=True)
class PropertyPrediction:
    """Predicted property mean and model uncertainty at one parameter point."""

    mean: Any
    standard_deviation: Any


class SklearnGaussianProcessSurrogate:
    """
    Learns one Gaussian process per property from usable observations.

    Inputs and outputs are standardized, and kernel parameters are learned from
    the data. Vector properties keep the same component order in every sample.
    Predicted uncertainty describes the model, rather than MD/FEP convergence.

    :param noise:
        Whether a learned white-noise term is added for noisy MD properties.
    """

    def __init__(self, noise: bool = False):
        """
        Initializes the surrogate; scikit-learn is required.

        With ``noise=True`` predicted standard deviations describe the latent
        property, not a new noisy observation, which is what expected improvement
        needs.
        """

        assert_msg_critical(
            importlib.util.find_spec("sklearn") is not None,
            "SklearnGaussianProcessSurrogate: scikit-learn is required")
        self.noise = noise
        self.models: dict[str, Any] = {}
        self.parameter_names: list[str] = []

    def fit(self, table: ParameterPropertyTable,
            property_names: Sequence[str]) -> None:
        """
        Fits one GP per property from usable, finite observations.

        Properties with fewer than two such observations are skipped (and absent
        from ``predict``). Raises ``ValueError`` for properties that are neither
        scalars nor fixed-length vectors.

        :param table:
            The observation table.
        :param property_names:
            The properties to fit.
        """

        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import Matern, WhiteKernel
        from sklearn.preprocessing import StandardScaler

        usable = [
            sample for sample in table.samples if
            sample.quality_status == "usable" and not sample.calculation.failed
        ]
        self.models = {}
        self.parameter_names = list(usable[0].parameters) if usable else []
        for name in property_names:
            rows = [
                sample for sample in usable
                if name in sample.calculation.values and
                all_finite(sample.calculation.values[name])
            ]
            if len(rows) < 2:
                continue
            parameter_values = np.array(
                [[sample.parameters[key]
                  for key in self.parameter_names]
                 for sample in rows])
            property_values = np.asarray(
                [sample.calculation.values[name] for sample in rows],
                dtype=float)
            is_vector = property_values.ndim == 2
            if property_values.ndim not in (1, 2):
                raise ValueError(
                    f"{name}: expected a scalar or fixed-length vector property"
                )
            input_scaler = StandardScaler()
            output_scaler = StandardScaler()
            scaled_inputs = input_scaler.fit_transform(parameter_values)
            scaled_outputs = output_scaler.fit_transform(
                property_values if is_vector else property_values[:, None])
            # Standardized outputs already have unit variance. Learn the spatial
            # length scale; constant observations contain no signal to tune it.
            constant_output = np.all(np.ptp(scaled_outputs, axis=0) == 0)
            kernel = Matern(length_scale=1.0, nu=2.5)
            if self.noise and not constant_output:
                kernel = kernel + WhiteKernel(noise_level=1e-2,
                                              noise_level_bounds=(1e-8, 1.0))
            model = GaussianProcessRegressor(
                kernel=kernel,
                alpha=1e-6,  # Stabilize nearly identical parameter points.
                optimizer=None if constant_output else "fmin_l_bfgs_b",
                n_restarts_optimizer=2 if self.noise else 0,
                random_state=0,
            )
            model.fit(scaled_inputs,
                      scaled_outputs if is_vector else scaled_outputs[:, 0])
            self.models[name] = (model, input_scaler, output_scaler, is_vector)

    def predict(
        self, parameters: Sequence[ParameterSet]
    ) -> dict[str, list[PropertyPrediction]]:
        """
        Predicts each fitted property in its original units.

        Standard deviations exclude the learned white noise, so they describe the
        latent property rather than a new noisy observation.

        :param parameters:
            The parameter points.

        :return:
            The predictions (mean and standard deviation) per point, by property.
        """

        if not parameters:
            return {name: [] for name in self.models}
        parameter_values = np.array(
            [[point[name]
              for name in self.parameter_names]
             for point in parameters])
        predictions = {}
        for name, (model, input_scaler, output_scaler,
                   is_vector) in self.models.items():
            means, deviations = model.predict(
                input_scaler.transform(parameter_values), return_std=True)
            means = output_scaler.inverse_transform(
                means if is_vector else means[:, None])
            deviations = np.asarray(deviations)
            noise = _learned_noise(model)
            if noise:
                deviations = np.sqrt(np.maximum(deviations**2 - noise, 0.0))
            if deviations.ndim == 1:
                deviations = deviations[:, None]
            deviations *= output_scaler.scale_
            predictions[name] = [
                PropertyPrediction(
                    mean=mean.tolist() if is_vector else float(mean[0]),
                    standard_deviation=deviation.tolist()
                    if is_vector else float(deviation[0]),
                ) for mean, deviation in zip(means, deviations)
            ]
        return predictions


def _learned_noise(model) -> float:
    """White-noise variance (standardized units) learned by a fitted GP, else 0."""
    from sklearn.gaussian_process.kernels import WhiteKernel

    noise = getattr(model.kernel_, "k2",
                    None)  # the kernel is Matern [+ WhiteKernel]
    return float(noise.noise_level) if isinstance(noise, WhiteKernel) else 0.0


# ==========================================================================
# search base
# ==========================================================================

CHECKPOINT_VERSION = 1
CHECKPOINT_GROUP = "autocdam/checkpoint"


def _json_value(value):
    """Serialize scientific summaries without converting numbers to strings."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(
        f"checkpoint contains unsupported value {type(value).__name__}")


class SearchDriver:
    """Shared state and bookkeeping of a parameter search (not used directly).

    ``make_inputs(proposals, sample_ids)`` yields one generated input (or an
    exception) per proposal. ``quality_check`` decides which observations count.
    ``objective_threshold`` replaces "all targets within tolerance" as the
    success criterion; required targets must hold in either case.
    """

    def __init__(self,
                 specs: Sequence[ParameterSpec],
                 targets: Sequence[FitTarget],
                 calculator,
                 evaluator,
                 make_inputs: Callable,
                 *,
                 seed: int = 42,
                 objective_threshold: float | None = None,
                 quality_check: Callable[[CalculationResult],
                                         QualityStatus] = default_quality_check,
                 noise: bool = False,
                 study_metadata: dict[str, Any] | None = None,
                 checkpoint_file: str | Path | None = None,
                 checkpoint_group: str = CHECKPOINT_GROUP,
                 comm=None,
                 ostream=None):
        """Validate the parameters and targets and create an empty observation table.

        It is a critical error if parameters or targets are missing or
        duplicated, an optimizable parameter has no bounds, or
        ``objective_threshold`` is negative or not finite. ``noise`` lets the property GPs learn observation noise.
        """
        self.specs = list(specs)
        self.targets = list(targets)
        assert_msg_critical(
            bool(self.specs) and bool(self.targets),
            "fitting requires parameters and targets")
        assert_msg_critical(
            len({s.name for s in self.specs}) == len(self.specs) and
            len({t.name for t in self.targets}) == len(self.targets),
            "parameter and target names must be unique")
        assert_msg_critical(
            objective_threshold is None or
            (math.isfinite(objective_threshold) and objective_threshold >= 0),
            "objective_threshold must be finite and nonnegative")
        for spec in self.specs:
            assert_msg_critical(
                not spec.optimizable or spec.bounds is not None,
                f"bounded sampling requires bounds for {spec.name}")
        self.calculator = calculator
        self.evaluator = evaluator
        self.make_inputs = make_inputs
        self.seed = seed
        self.objective_threshold = objective_threshold
        self.quality_check = quality_check
        self.surrogate = SklearnGaussianProcessSurrogate(noise=noise)
        self.table = ParameterPropertyTable(
            study_metadata=dict(study_metadata or {}))
        self.checkpoint_file = (h5_filename(checkpoint_file)
                                if checkpoint_file else None)
        self.checkpoint_group = checkpoint_group
        # only the MPI master rank writes the checkpoint
        self.comm = comm if comm is not None else MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.ostream = ostream
        self.history: list[dict[str, Any]] = []

    @property
    def property_names(self) -> list[str]:
        """Names of the fitted targets, which are also the surrogate's property names."""
        return [t.name for t in self.targets]

    def _info(self, message: str) -> None:
        """Print an information line to ``ostream``, if one was given."""
        if self.ostream is not None:
            self.ostream.print_info(message)
            self.ostream.flush()

    # ------------------------------------------------------------- evaluation

    def _observed_best(
            self
    ) -> tuple[ParameterPropertySample | None, EvaluationResult | None]:
        """Prefer a measured tolerance solution, then minimize the observed loss."""
        scores = self.table.rescore(self.evaluator, self.targets)
        rows = [
            r for r in self.table.samples
            if r.sample_id in scores and not scores[r.sample_id].failed and
            math.isfinite(scores[r.sample_id].objective)
        ]
        if not rows:
            return None, None
        # A tolerance solution wins; otherwise required constraints outrank the
        # weighted loss, so a low score cannot compensate for a violated target.
        satisfying = [
            r for r in rows if self._meets_threshold(scores[r.sample_id])
        ]
        row = min(
            satisfying or rows,
            key=lambda r:
            (not scores[r.sample_id].feasible, scores[r.sample_id].objective))
        return row, scores[row.sample_id]

    def _meets_threshold(self, evaluation: EvaluationResult | None) -> bool:
        """Required targets always gate success; then threshold or tolerances."""
        if evaluation is None or evaluation.failed or not evaluation.feasible:
            return False
        if self.objective_threshold is not None:
            return evaluation.objective <= self.objective_threshold
        return within_tolerance(evaluation, self.targets)

    # ------------------------------------------------------------ calculation

    def _calculate(self, inputs: CalculationInput | Exception,
                   parameters: ParameterSet, **overrides) -> CalculationResult:
        """Calculate one generated input; any error becomes a failed result.

        KeyboardInterrupt and SystemExit are not intercepted.
        """
        try:
            if isinstance(inputs, Exception):
                raise inputs
            return self.calculator.calculate(
                inputs.batch, parameters, self.targets, **{
                    **inputs.calculator_kwargs,
                    **overrides
                })
        except Exception as exc:
            return CalculationResult(
                failed=True, failure_reason=f"{type(exc).__name__}: {exc}")

    def _collect_stage(self, stage: int,
                       proposals: Sequence[ParameterSet]) -> list[str]:
        """Generate each proposal lazily, calculate it, and keep individual failures.

        Raises ``ValueError`` if ``quality_check`` returns an unknown status.
        """
        offset = len(self.table.samples)
        sample_ids = [
            f"stage-{stage:03d}-{offset + index:06d}"
            for index in range(len(proposals))
        ]
        inputs = iter(self.make_inputs(proposals, sample_ids))
        for parameters, sample_id in zip(proposals, sample_ids):
            started = time.monotonic()
            item = next(inputs)
            calculation = self._calculate(item, parameters)
            quality = self.quality_check(calculation)
            if quality not in ("usable", "unreviewed", "rejected"):
                raise ValueError(
                    "quality_check must return usable, unreviewed, or rejected")
            self.table.add(
                ParameterPropertySample(
                    sample_id=sample_id,
                    parameters=parameters,
                    calculation=calculation,
                    quality_status="rejected"
                    if calculation.failed else quality,
                    protocol_id=getattr(item, "protocol_id", None),
                    runtime_seconds=time.monotonic() - started))
        return sample_ids

    def _history_entry(self, best, evaluation) -> dict[str, Any]:
        """History fields every search records: the measured best so far."""
        objective = evaluation.objective if evaluation is not None else None
        return {
            "parameters": best.parameters.to_dict() if best is not None else {},
            "objective": objective,
            "best_objective": objective,
            "metrics": evaluation.metrics if evaluation is not None else {},
            "errors": evaluation.errors if evaluation is not None else {},
        }

    # ------------------------------------------------------------ checkpoints

    def _signature(self) -> dict[str, Any]:
        """Settings that must agree when continuing a saved study."""
        return {
            "specs": [asdict(s) for s in self.specs],
            "targets": [asdict(t) for t in self.targets],
            "study_metadata": self.table.study_metadata,
            "seed": self.seed,
            "objective_threshold": self.objective_threshold,
        }

    def _checkpoint_state(self) -> dict[str, Any]:
        """Search-specific checkpoint entries besides signature, history and samples."""
        return {}

    def _checkpoint_label(self) -> str | None:
        """The checkpoint location for printing: ``file.h5:group``, or None."""
        if self.checkpoint_file is None:
            return None
        return f"{self.checkpoint_file}:{self.checkpoint_group}"

    def _save(self) -> None:
        """Replace the checkpoint only after writing the complete next version."""
        if self.checkpoint_file is None or self.rank != mpi_master():
            return
        data = {
            "version": CHECKPOINT_VERSION,
            "signature": self._signature(),
            **self._checkpoint_state(), "history": self.history,
            "samples": [asdict(r) for r in self.table.samples]
        }
        # Plain values (as JSON would store them), then one HDF5 group.
        data = json.loads(json.dumps(data, default=_json_value))
        write_h5_groups(self.checkpoint_file, {self.checkpoint_group: data})

    def _load_checkpoint(self, differ_message: str) -> dict[str, Any]:
        """Check a saved study matches this one, restore its samples, return the data."""
        data = (read_h5_group(self.checkpoint_file, self.checkpoint_group)
                if self.checkpoint_file is not None else None)
        if data is None:
            raise ValueError("restart requires an existing checkpoint in the "
                             "run file")
        signature = json.loads(
            json.dumps(self._signature(), default=_json_value))
        if (data.get("version") != CHECKPOINT_VERSION or
                data["signature"] != signature):
            raise ValueError(differ_message)
        self.history = data["history"]
        self.table.samples = []
        for item in data["samples"]:
            item["parameters"] = ParameterSet(**item["parameters"])
            item["calculation"] = CalculationResult(**item["calculation"])
            self.table.add(ParameterPropertySample(**item))
        return data

    # ----------------------------------------------------------------- result

    def _build_result(self, best, evaluation, converged,
                      stop_reason) -> FitResult:
        """Measured properties of the best sample, with quality and protocol diagnostics."""
        samples = self.table.samples
        metadata = [
            item for sample in samples
            for item in sample.calculation.metadata.values()
            if isinstance(item, dict)
        ]
        parameters = (best.parameters.to_dict() if best is not None else {
            spec.name: spec.initial for spec in self.specs
        })
        status = "success" if converged else (
            "warning" if best is not None else "failed")
        found = evaluation is not None
        best_properties = best.calculation.values if best is not None else {}
        violations = dict(evaluation.constraint_violations) if found else {}
        return FitResult(
            parameters=ParameterSet(parameters),
            metrics=evaluation.metrics if evaluation is not None else {},
            objective=evaluation.objective
            if evaluation is not None else math.inf,
            status=status,
            converged=converged,
            iterations=len(self.history),
            history=list(self.history),
            warnings=[
                warning for sample in samples
                for warning in sample.calculation.warnings
            ],
            diagnostics={
                "stop_reason": stop_reason,
                "n_calculations": len(samples),
                "best_sample_id": best.sample_id if best is not None else None,
                "best_properties": best_properties,
                "feasible": bool(found and evaluation.feasible and
                                 not evaluation.failed),
                "constraint_violations": violations,
                "target_results": evaluation.target_results if found else {},
                "quality_counts": {
                    quality: sum(s.quality_status == quality for s in samples)
                    for quality in ("usable", "unreviewed", "rejected")
                },
                "fast_mode": any(item.get("fast_mode") for item in metadata),
                "scientifically_validated": False,
            },
            checkpoint=(str(self.checkpoint_file)
                        if self.checkpoint_file is not None else None),
            checkpoint_group=(self.checkpoint_group
                              if self.checkpoint_file is not None else None),
        )


# ==========================================================================
# staged
# ==========================================================================

CANDIDATE_POOL_SIZE = 256  # points scored by the surrogate per stage
EXPLORATION_FRACTION = 0.2  # share of each later stage sampled from the original bounds


class FittingDriver(SearchDriver):
    """
    Generates, calculates, refits, and refines stages until tolerance or budget.

    :param specs:
        The parameter specifications.
    :param targets:
        The fitted targets.
    :param calculator:
        The property calculator.
    :param evaluator:
        The evaluator.
    :param make_inputs:
        The function ``make_inputs(proposals, sample_ids)`` that yields one
        generated input (or an exception) per proposal.
    :param max_iterations:
        The number of stages.
    :param samples_per_stage:
        The number of candidates per stage.
    :param shrink_fraction:
        The width kept per refinement, in (0, 1).
    :param options:
        The shared search options (seed, objective_threshold, quality_check,
        noise, study_metadata, checkpoint settings, comm, ostream).
    """

    def __init__(self,
                 specs,
                 targets,
                 calculator,
                 evaluator,
                 make_inputs,
                 *,
                 max_iterations: int = 3,
                 samples_per_stage: int = 10,
                 shrink_fraction: float = 0.5,
                 **options):
        """
        Sets the stage budget; the region starts at the full bounds.

        ``shrink_fraction`` is the width kept per refinement. Out-of-range
        sizes or fractions are a critical error.
        """

        super().__init__(specs, targets, calculator, evaluator, make_inputs,
                         **options)
        assert_msg_critical(
            max_iterations >= 1 and samples_per_stage >= 2,
            "max_iterations >= 1 and samples_per_stage >= 2 are required")
        assert_msg_critical(0 < shrink_fraction < 1,
                            "shrink_fraction must lie in (0, 1)")
        self.max_iterations = max_iterations
        self.samples_per_stage = samples_per_stage
        self.shrink_fraction = shrink_fraction
        self.region = {s.name: s.bounds for s in self.specs if s.optimizable}

    # ------------------------------------------------------------------ fit

    def fit(self, *, resume: bool = False) -> FitResult:
        """
        Calculates complete stages until a measured solution or the budget limit.

        :param resume:
            Whether to resume from the checkpoint.

        :return:
            The fit result.
        """

        if resume:
            self._resume()
        self._print_setup(resume)
        best, evaluation = self._observed_best()
        converged = self._meets_threshold(evaluation)
        stop_reason = "threshold" if converged else "max_iterations"
        if self.table.samples:
            self.surrogate.fit(self.table, self.property_names)

        for stage in range(len(self.history) + 1, self.max_iterations + 1):
            if converged:
                break
            self._info(
                f"Fitting stage {stage}/{self.max_iterations}: generate, calculate, refit"
            )
            used_bounds = dict(self.region)
            proposals, predicted = self._propose_stage(stage)
            if not proposals:
                stop_reason = "search_exhausted"
                break
            sample_ids = self._collect_stage(stage, proposals)
            best, evaluation = self._observed_best()
            self.surrogate.fit(self.table, self.property_names)
            converged = self._meets_threshold(evaluation)
            if best is not None and not converged:
                self._refine(best.parameters)
            self._record_stage(stage, used_bounds, sample_ids, predicted, best,
                               evaluation)
            if converged:
                stop_reason = "threshold"
                break
        return self._build_result(best, evaluation, converged, stop_reason)

    # ------------------------------------------------------------ proposals

    def _region_specs(self) -> list[ParameterSpec]:
        """The specs restricted to the current region (fixed values untouched)."""
        specs = []
        for spec in self.specs:
            if spec.optimizable:
                lower, upper = self.region[spec.name]
                spec = replace(spec,
                               initial=min(max(spec.initial, lower), upper),
                               bounds=(lower, upper))
            specs.append(spec)
        return specs

    def _predicted_best(self,
                        pool: Sequence[ParameterSet]) -> ParameterSet | None:
        """Pool point with the best predicted score (evaluators opt in via evaluate_predictions)."""
        score_predictions = getattr(self.evaluator, "evaluate_predictions",
                                    None)
        if score_predictions is None:
            return None
        predictions = self.surrogate.predict(pool)
        if any(t.name not in predictions for t in self.targets):
            return None
        best, value = None, (True, math.inf)
        for i, parameters in enumerate(pool):
            result = score_predictions(
                {t.name: predictions[t.name][i].mean for t in self.targets},
                self.targets)
            key = (not result.feasible, result.objective)
            if not result.failed and math.isfinite(
                    result.objective) and key < value:
                best, value = parameters, key
        return best

    def _propose_stage(
            self, stage: int) -> tuple[list[ParameterSet], ParameterSet | None]:
        """Predicted point + local samples + global exploration, without repeats."""
        region_specs = self._region_specs()
        stage_seed = self.seed + stage * 3
        predicted = None
        if self.table.samples:
            predicted = self._predicted_best(
                latin_hypercube(region_specs, CANDIDATE_POOL_SIZE, stage_seed))

        proposals = latin_hypercube(region_specs, self.samples_per_stage,
                                    stage_seed + 1)
        proposals[0] = predicted if predicted is not None else ParameterSet(
            {spec.name: spec.initial for spec in region_specs})
        if stage > 1:
            count = int(self.samples_per_stage * EXPLORATION_FRACTION)
            if count:
                proposals[-count:] = latin_hypercube(self.specs, count,
                                                     stage_seed + 2)

        seen = {
            self._parameter_key(sample.parameters)
            for sample in self.table.samples
        }
        fresh = []
        for proposal in proposals:
            key = self._parameter_key(proposal)
            if key not in seen:
                seen.add(key)
                fresh.append(proposal)
        return fresh, predicted

    def _parameter_key(self, parameters: ParameterSet) -> tuple[float, ...]:
        """Hashable tuple of all parameter values, used to skip repeated points."""
        return tuple(parameters[spec.name] for spec in self.specs)

    def _refine(self, parameters: ParameterSet) -> None:
        """Center a smaller region on the measured best, inside the global bounds."""
        for spec in self.specs:
            if not spec.optimizable:
                continue
            lo, hi = self.region[spec.name]
            global_lo, global_hi = spec.bounds
            width = max((hi - lo) * self.shrink_fraction,
                        (global_hi - global_lo) * 1e-8)
            lower = max(
                global_lo,
                min(parameters[spec.name] - width / 2, global_hi - width))
            self.region[spec.name] = (lower, lower + width)

    # --------------------------------------------------------- bookkeeping

    def _record_stage(self, stage, bounds, sample_ids, predicted, best,
                      evaluation) -> None:
        """Append the stage's region, samples, best and predicted point; checkpoint."""
        predicted = predicted.to_dict() if predicted is not None else None
        self.history.append({
            "iteration": stage,
            "region_bounds": bounds,
            "next_region_bounds": dict(self.region),
            "sample_ids": sample_ids,
            "n_calculations": len(sample_ids),
            **self._history_entry(best, evaluation),
            "predicted_candidate": predicted,
        })
        self._save()

    def _signature(self) -> dict[str, Any]:
        """Settings that must match to resume a saved staged study."""
        return {
            **super()._signature(), "samples_per_stage": self.samples_per_stage,
            "shrink_fraction": self.shrink_fraction,
            "exploration_fraction": EXPLORATION_FRACTION,
            "candidate_pool_size": CANDIDATE_POOL_SIZE
        }

    def _checkpoint_state(self) -> dict[str, Any]:
        """Save the current refinement region with the checkpoint."""
        return {"region": self.region}

    def _resume(self) -> None:
        """Restore samples, history and region from the checkpoint.

        Raises ``ValueError`` if the settings differ or more stages were already
        completed than ``max_iterations`` allows.
        """
        data = self._load_checkpoint(
            "checkpoint parameter/target/protocol settings differ from this study"
        )
        self.region = {k: tuple(v) for k, v in data["region"].items()}
        if len(self.history) > self.max_iterations:
            raise ValueError(
                "max_iterations is smaller than the checkpoint's completed stage count"
            )
        for entry in self.history:
            for name in ("region_bounds", "next_region_bounds"):
                entry[name] = {
                    key: tuple(value) for key, value in entry[name].items()
                }

    def _print_setup(self, resume: bool) -> None:
        """Print the parameters, targets and stage settings to ``ostream``."""
        rows = []
        for spec in self.specs:
            unit = f" {spec.unit}" if spec.unit else ""
            if spec.optimizable:
                lower, upper = self.region[spec.name]
                text = (
                    f"[{format_value(lower)}, {format_value(upper)}]{unit}, "
                    f"initial {format_value(spec.initial)}")
            else:
                text = f"fixed at {format_value(spec.initial)}{unit}"
            rows.append((f"Parameter {spec.name}", text))
        rows += [(
            f"Target {t.name}",
            f"{format_value(t.reference)}{' ' + t.unit if t.unit else ''}"
            f" +/- {format_value(t.tolerance)}, weight {format_value(t.weight)}, {t.metric}"
        ) for t in self.targets]
        quality_check = getattr(self.quality_check, "__name__",
                                type(self.quality_check).__name__)
        print_settings(self.ostream, "CDAM Staged Fitting Setup", [
            *rows,
            ("Maximum Stages", self.max_iterations),
            ("Samples per Stage", self.samples_per_stage),
            ("Region Shrink Fraction", self.shrink_fraction),
            ("Global Exploration Fraction", EXPLORATION_FRACTION),
            ("Surrogate Candidate Pool", CANDIDATE_POOL_SIZE),
            ("Surrogate Model", type(self.surrogate).__name__),
            ("Stopping Criterion",
             f"objective <= {format_value(self.objective_threshold)}"
             if self.objective_threshold is not None else
             "all targets within tolerance"),
            ("Quality Check", quality_check),
            ("Random Seed", self.seed),
            ("Checkpoint", self._checkpoint_label()),
            ("Resumed Stages",
             f"{len(self.history)} stages, {len(self.table.samples)} samples"
             if resume else None),
        ])
