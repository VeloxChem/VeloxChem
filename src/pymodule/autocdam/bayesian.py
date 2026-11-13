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
"""Gaussian-process Bayesian search over the bounded CDAM parameters.

A quasi-random design is followed by acquisition steps. Each property has its
own GP, so the loss is computed from predicted properties and target weights
can change without recalculating. Acquisitions use Monte Carlo draws of the
property predictions with fixed random numbers per step (deterministic for a
seed):

* ``weighted``: noisy expected improvement of the normalized loss over the
  best predicted loss at an evaluated point;
* ``constrained``: expected improvement times the probability that every
  required target is within tolerance (that probability alone until a
  feasible point has been calculated);
* ``pareto``: ParEGO-style random augmented-Chebyshev scalarizations of the
  per-target squared errors; the non-dominated observations are reported.

Batches avoid duplicates and use local penalization around points already
chosen. The search stops at the tolerance criterion, the evaluation budget,
or when the acquisition stays below ``acquisition_tolerance`` for
``patience`` steps.

Soft bounds (``ParameterSpec.expandable``, e.g. UFF-scaled defaults) widen
when the best candidate of an acquisition step lies within 2% of one: that
side moves out by ``bound_expansion_factor`` times the width, never past
``hard_bounds``, at most ``max_bound_expansions`` times, and each widening
adds ``bound_expansion_evaluations`` to the budget. Hard bounds never move;
a solution there is reported as a bound warning.

``cost_aware`` divides positive acquisition values by the predicted wall time
(a GP on log runtime), i.e. expected improvement per second. ``screening``
``{"calculator_options": {...}, "framework_protocol": {...}, "oversample": 4}``
calculates ``oversample * batch_size`` candidates with a cheap protocol and
sends only the best ``batch_size`` to the full calculation. Screening values
never train the surrogate; their agreement with full results is reported.
"""
from __future__ import annotations

import warnings
from dataclasses import asdict, replace
from typing import Any

import numpy as np

from ..errorhandler import assert_msg_critical
from .evaluator import Evaluator
from .models import (FitResult, ParameterSet, _is_sequence, format_value,
                     print_settings, widen_soft_bounds)
from .search import SearchDriver

SEARCH_DEFAULTS = {
    "mode": "weighted",  # weighted | constrained | pareto
    "max_evaluations": 30,  # total calculated candidates, including the design
    "initial_samples": None,  # default: max(5, 2 * n_parameters + 1)
    "initial_design": "sobol",  # sobol | lhs
    "batch_size": 1,  # candidates proposed per acquisition step
    "seed": 42,
    "pool_size": 512,  # quasi-random acquisition candidates per step
    "n_mc": 256,  # Monte Carlo draws of the property predictions
    "noise": True,  # learn observation noise in each GP
    "local_refinement": True,  # polish the best pool point with L-BFGS-B
    "duplicate_tolerance": 1e-3,  # minimum distance (unit cube) to evaluated points
    "acquisition_tolerance": 1e-6,  # stop after `patience` steps below this value
    "patience": 3,
    "cost_aware": False,  # expected improvement per predicted second
    "cost_exponent": 1.0,
    "screening": None,  # cheap pre-ranking; see the module docstring
    "max_bound_expansions": 3,  # widenings of soft (expandable) bounds per search
    "bound_expansion_factor": 1.0,  # a widened side moves by this fraction of the width
    "bound_expansion_evaluations": 3,  # extra budget granted per widening
}
SCREENING_KEYS = {"calculator_options", "framework_protocol", "oversample"}
LOCAL_POOL_SPREAD = 0.08  # unit-cube spread of extra pool points around the best
PENALTY_RADIUS = 0.10  # local penalization radius between batch members
PARETO_RHO = 0.05  # augmentation weight of the Chebyshev scalarization


def resolve_search_settings(search, specs, targets) -> dict[str, Any]:
    """Defaults plus ``search``, validated against the parameters and targets."""
    unknown = set(search or {}) - set(SEARCH_DEFAULTS)
    assert_msg_critical(not unknown,
                        f"unsupported search settings: {sorted(unknown)}")
    settings = {**SEARCH_DEFAULTS, **(search or {})}
    if settings["initial_samples"] is None:
        settings["initial_samples"] = max(
            5, 2 * sum(spec.optimizable for spec in specs) + 1)
    assert_msg_critical(
        settings["mode"] in ("weighted", "constrained", "pareto"),
        "search mode must be 'weighted', 'constrained', or 'pareto'")
    assert_msg_critical(settings["initial_design"] in ("sobol", "lhs"),
                        "initial_design must be 'sobol' or 'lhs'")
    for key in ("max_evaluations", "initial_samples", "batch_size", "pool_size",
                "n_mc", "patience"):
        assert_msg_critical(
            int(settings[key]) >= 1,
            f"search setting {key} must be a positive integer")
    assert_msg_critical(
        int(settings["max_bound_expansions"]) >= 0 and
        int(settings["bound_expansion_evaluations"]) >= 0,
        "bound expansion counts must be nonnegative")
    assert_msg_critical(
        float(settings["bound_expansion_factor"]) > 0,
        "bound_expansion_factor must be positive")
    assert_msg_critical(
        settings["initial_samples"] >= 2,
        "initial_samples must be at least 2 to fit a surrogate")
    assert_msg_critical(
        settings["initial_samples"] <= settings["max_evaluations"],
        "initial_samples cannot exceed max_evaluations")
    assert_msg_critical(
        settings["mode"] != "constrained" or
        any(t.is_required for t in targets),
        "constrained mode needs at least one target with role='required'")
    assert_msg_critical(
        settings["mode"] != "pareto" or
        sum(t.is_objective for t in targets) >= 2,
        "pareto mode needs at least two required or weighted targets")
    screening = settings["screening"]
    if screening is not None:
        unknown = set(screening) - SCREENING_KEYS
        assert_msg_critical(
            not unknown, f"unsupported screening settings: {sorted(unknown)}")
        assert_msg_critical(
            int(screening.get("oversample", 4)) >= 2,
            "screening oversample must be at least 2")
        assert_msg_critical(
            bool(screening.get("calculator_options") or
                 screening.get("framework_protocol")),
            "screening needs cheaper calculator_options or framework_protocol")
    return settings


def _quasi_random(d: int, n: int, seed: int):
    """Return ``n`` scrambled Sobol points in the ``d``-dimensional unit cube, shape (n, d)."""
    from scipy.stats import qmc
    with warnings.catch_warnings(
    ):  # Sobol balance warning for non-powers of two
        warnings.simplefilter("ignore")
        return qmc.Sobol(d=d, scramble=True, seed=seed).random(n)


class BayesianFittingDriver(SearchDriver):
    """
    Bayesian search; see the module docstring for modes, bounds and screening.

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
    :param search:
        The search settings (see ``SEARCH_DEFAULTS``).
    :param options:
        The shared search options (objective_threshold, quality_check,
        study_metadata, checkpoint settings, comm, ostream).
    """

    def __init__(self,
                 specs,
                 targets,
                 calculator,
                 evaluator,
                 make_inputs,
                 *,
                 search: dict[str, Any] | None = None,
                 **options):
        """
        Resolves and validates ``search`` and sets up the driver.

        One GP per property learns observation noise when ``search["noise"]`` is
        true. Unknown or inconsistent search settings are a critical error.
        """

        settings = resolve_search_settings(search, specs, targets)
        super().__init__(specs,
                         targets,
                         calculator,
                         evaluator,
                         make_inputs,
                         seed=int(settings["seed"]),
                         noise=bool(settings["noise"]),
                         **options)
        self.search = settings
        self.initial_specs = list(self.specs)
        self.bound_expansions: list[dict[str, Any]] = []
        self.free = [spec for spec in self.specs if spec.optimizable]
        self.objective_targets = [t for t in self.targets if t.is_objective]

    # ------------------------------------------------------------------ fit

    def fit(self, *, resume: bool = False) -> FitResult:
        """
        Runs the design, then acquisition steps, until tolerance, budget, or
        stagnation.

        :param resume:
            Whether to resume from the checkpoint.

        :return:
            The fit result.
        """

        if resume:
            self._resume()
        self._print_setup()
        for target in self.objective_targets:
            assert_msg_critical(
                not _is_sequence(target.reference),
                f"Bayesian search supports scalar targets; {target.name} is a "
                "vector (use the staged search)")
        best, evaluation = self._observed_best()
        converged = self._meets_threshold(evaluation)
        stop_reason = "threshold" if converged else "max_evaluations"
        budget = self._budget()
        design = self._design()
        patience, tolerance = int(
            self.search["patience"]), self.search["acquisition_tolerance"]
        quiet = sum(1 for entry in self.history[-patience:]
                    if entry.get("max_acquisition") is not None and
                    entry["max_acquisition"] < tolerance)

        while not converged and len(self.table.samples) < budget:
            step = len(self.history) + 1
            n_done = len(self.table.samples)
            if n_done < len(design):
                phase, acquisition, screened = "design", [], None
                proposals = [self._from_unit(u) for u in design[n_done:]]
            else:
                phase = "acquisition"
                proposals, acquisition, screened = self._acquisition_step(step)
                if not proposals:
                    stop_reason = "search_exhausted"
                    break
            proposals = proposals[:budget - n_done]
            self._info(
                f"Bayesian step {step}: {phase}, {len(proposals)} candidate(s), "
                f"{n_done}/{budget} evaluated")

            sample_ids = self._collect_stage(step, proposals)
            best, evaluation = self._observed_best()
            converged = self._meets_threshold(evaluation)
            max_acquisition = max(acquisition) if acquisition else None
            # Widen only from the search's own proposals: a space-filling design
            # point that lands near a bound says nothing about where the optimum is.
            expansions = ([] if converged or phase == "design" else
                          self._expand_bounds(best, step))
            if expansions:
                budget = self._budget()
                quiet = 0
            self._record_step(step, phase, sample_ids, acquisition,
                              max_acquisition, best, evaluation, expansions,
                              screened)
            if converged:
                stop_reason = "threshold"
                break
            if max_acquisition is not None and not expansions:
                quiet = quiet + 1 if max_acquisition < tolerance else 0
                if quiet >= patience:
                    stop_reason = "acquisition_exhausted"
                    break
        if self.table.samples:
            self.surrogate.fit(self.table, self.property_names)
        return self._build_result(best, evaluation, converged, stop_reason)

    def _acquisition_step(self, step):
        """Refit the surrogate and propose (and optionally screen) the next batch."""
        self.surrogate.fit(self.table, self.property_names)
        screening = self.search["screening"]
        count = (int(screening.get("oversample", 4)) *
                 int(self.search["batch_size"])
                 if screening is not None else None)
        proposals, acquisition = self._propose_batch(step, count)
        screened = None
        if screening is not None and proposals:
            proposals, acquisition, screened = self._screen(
                step, proposals, acquisition)
        return proposals, acquisition, screened

    def _record_step(self, step, phase, sample_ids, acquisition,
                     max_acquisition, best, evaluation, expansions,
                     screened) -> None:
        """Append one step (samples, acquisition, best, bounds, failures) and checkpoint."""
        failed = [
            r for r in self.table.samples
            if r.sample_id in sample_ids and r.calculation.failed
        ]
        failures = [{
            "sample_id": r.sample_id,
            "reason": r.calculation.failure_reason
        } for r in failed]
        self.history.append({
            "iteration": step,
            "phase": phase,
            "sample_ids": sample_ids,
            "n_calculations": len(sample_ids),
            "n_evaluated": len(self.table.samples),
            "acquisition": acquisition,
            "max_acquisition": max_acquisition,
            **self._history_entry(best, evaluation),
            "feasible": bool(evaluation is not None and evaluation.feasible),
            "failures": failures,
            "bounds": self.current_bounds(),
            "bound_expansions": expansions,
            "screening":
                (self._link_screening(screened, sample_ids)
                 if phase == "acquisition" and screened is not None else None),
        })
        self._save()

    # ----------------------------------------------------------- soft bounds

    def _set_bounds(self, bounds: dict[str, tuple[float, float]]) -> None:
        """Replace the bounds of the named parameters and refresh the list of free specs."""
        self.specs = [
            replace(s, bounds=tuple(bounds[s.name])) if s.name in bounds else s
            for s in self.specs
        ]
        self.free = [spec for spec in self.specs if spec.optimizable]

    def current_bounds(self) -> dict[str, tuple[float, float]]:
        """
        Gets the current (possibly widened) bounds of the optimizable parameters.

        :return:
            The bounds by parameter name.
        """

        return {s.name: tuple(s.bounds) for s in self.free}

    def _n_widenings(self) -> int:
        """Number of acquisition steps that widened at least one soft bound."""
        return len({e["step"] for e in self.bound_expansions})

    def _progress_total(self) -> int:
        """The evaluation budget, including widenings so far."""
        return self._budget()

    def _budget(self) -> int:
        """``max_evaluations`` plus the extra evaluations granted per widening."""
        return (int(self.search["max_evaluations"]) + self._n_widenings() *
                int(self.search["bound_expansion_evaluations"]))

    def _expand_bounds(self, best, step: int) -> list[dict[str, Any]]:
        """Widen the soft bounds the best candidate has reached; return the changes."""
        if best is None or self._n_widenings() >= int(
                self.search["max_bound_expansions"]):
            return []
        _, changes = widen_soft_bounds(self.free, best.parameters,
                                       float(self.search["bound_expansion_factor"]),
                                       step)
        bounds = {c["parameter"]: tuple(c["new"]) for c in changes}
        if changes:
            self._set_bounds(bounds)
            self.bound_expansions.extend(changes)
            for c in changes:
                self._info(
                    f"Best {c['parameter']} = {c['value']:.4g} reached its {c['side']} "
                    f"bound; widened [{c['old'][0]:.4g}, {c['old'][1]:.4g}] -> "
                    f"[{c['new'][0]:.4g}, {c['new'][1]:.4g}]")
        return changes

    # ------------------------------------------------- unit-cube coordinates

    def _to_unit(self, parameters: ParameterSet) -> np.ndarray:
        """Map the free parameters to unit-cube coordinates within the current bounds."""
        return np.array([
            (parameters[s.name] - s.bounds[0]) / (s.bounds[1] - s.bounds[0])
            for s in self.free
        ])

    def _from_unit(self, unit: np.ndarray) -> ParameterSet:
        """Map unit-cube coordinates (clipped to [0, 1]) to a full parameter set."""
        values = {spec.name: spec.initial for spec in self.specs}
        for spec, fraction in zip(self.free, np.clip(unit, 0.0, 1.0)):
            lower, upper = spec.bounds
            values[spec.name] = float(lower + fraction * (upper - lower))
        return ParameterSet(values)

    def _evaluated_units(self) -> np.ndarray:
        """Unit-cube coordinates of every calculated sample, failed ones included."""
        return np.array(
            [self._to_unit(sample.parameters) for sample in self.table.samples])

    def _design(self) -> np.ndarray:
        """Initial space-filling design in the unit cube (Sobol or Latin hypercube), seeded."""
        n = int(self.search["initial_samples"])
        if self.search["initial_design"] == "lhs":
            from scipy.stats import qmc
            return qmc.LatinHypercube(d=len(self.free),
                                      seed=self.seed).random(n)
        return _quasi_random(len(self.free), n, self.seed)

    # ----------------------------------------------------------- acquisition

    def _predict(self, units: np.ndarray):
        """Property means and latent standard deviations, shape (n, n_targets)."""
        points = [self._from_unit(unit) for unit in units]
        predictions = self.surrogate.predict(points)
        if any(t.name not in predictions for t in self.objective_targets):
            return None
        columns = [predictions[t.name] for t in self.objective_targets]
        rows = range(len(points))
        mean = np.array([[column[i].mean for column in columns] for i in rows],
                        dtype=float)
        sd = np.array(
            [[c[i].standard_deviation for c in columns] for i in rows],
            dtype=float)
        return mean, np.maximum(sd, 1e-12)

    def _draw_scores(self, units, normals, weights_lambda=None):
        """Loss (or scalarized loss) and feasibility for every MC draw."""
        predicted = self._predict(units)
        if predicted is None:
            return None
        mean, sd = predicted
        targets = self.objective_targets
        draws = mean[:, None, :] + sd[:, None, :] * normals[None, :, :]
        # The evaluator's own error metric, so acquisitions optimize what is reported.
        errors, within = Evaluator.score_draws(draws, targets)
        weights = np.array([t.weight for t in targets])
        if weights_lambda is not None:
            loss = np.max(errors * weights_lambda,
                          axis=2) + PARETO_RHO * (errors @ weights_lambda)
        elif weights.sum() > 0:
            loss = errors @ weights / weights.sum()
        else:
            loss = np.zeros(errors.shape[:2])
        required = np.array([t.is_required for t in targets])
        feasible = np.all(within | ~required, axis=2)
        return loss, feasible

    def _incumbent(self, normals, weights_lambda=None) -> float | None:
        """Best predicted (noise-free) score at an evaluated usable point."""
        scores = self.table.rescore(self.evaluator, self.targets)
        rows = [
            r for r in self.table.samples
            if r.sample_id in scores and not scores[r.sample_id].failed
        ]
        if self.search["mode"] == "constrained":
            rows = [r for r in rows if scores[r.sample_id].feasible]
        if not rows:
            return None
        units = np.array([self._to_unit(r.parameters) for r in rows])
        result = self._draw_scores(units, np.zeros((1, normals.shape[1])),
                                   weights_lambda)
        if result is None:
            return None
        return float(result[0][:, 0].min())

    def _acquisition(self,
                     units,
                     normals,
                     incumbent,
                     weights_lambda=None) -> np.ndarray:
        """Acquisition value of each unit-cube point for the current mode.

        Expected improvement over ``incumbent`` averaged over the MC draws
        (times feasibility in constrained mode). Without an incumbent it falls back
        to the feasibility probability or the negative mean loss, and without a
        usable surrogate to the distance from evaluated points.
        """
        result = self._draw_scores(units, normals, weights_lambda)
        if result is None:
            return self._space_filling(units)
        loss, feasible = result
        if self.search["mode"] == "constrained":
            if incumbent is None:
                return feasible.mean(axis=1)
            return (np.maximum(incumbent - loss, 0.0) * feasible).mean(axis=1)
        if incumbent is None:
            return -loss.mean(axis=1)
        return np.maximum(incumbent - loss, 0.0).mean(axis=1)

    def _space_filling(self, units) -> np.ndarray:
        """Distance of each point to its nearest evaluated point (ones if none)."""
        evaluated = self._evaluated_units()
        if not len(evaluated):
            return np.ones(len(units))
        return np.min(np.linalg.norm(units[:, None] - evaluated[None], axis=2),
                      axis=1)

    def _cost_factor(self, pool: np.ndarray) -> np.ndarray | None:
        """Predicted wall time of pool candidates from earlier runtimes, or None."""
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import Matern, WhiteKernel

        rows = [
            s for s in self.table.samples
            if s.runtime_seconds and s.runtime_seconds > 0
        ]
        if len(rows) < 3:
            return None
        x = np.array([self._to_unit(s.parameters) for s in rows])
        y = np.log([s.runtime_seconds for s in rows])
        if np.ptp(y) == 0:
            return None
        model = GaussianProcessRegressor(Matern(length_scale=0.3, nu=2.5) +
                                         WhiteKernel(1e-2),
                                         normalize_y=True,
                                         random_state=0).fit(x, y)
        return np.exp(model.predict(pool))

    def _propose_batch(self, step: int, count: int | None = None):
        """Up to ``count`` (default batch_size) new points and their acquisition values."""
        from scipy.optimize import minimize

        rng = np.random.default_rng(self.seed + 7919 * step)
        n_pool = int(self.search["pool_size"])
        normals = rng.standard_normal(
            (int(self.search["n_mc"]), len(self.objective_targets)))
        pool = _quasi_random(len(self.free), n_pool, self.seed + step)
        best, _ = self._observed_best()
        if best is not None:  # local candidates around the measured best
            local = self._to_unit(
                best.parameters) + LOCAL_POOL_SPREAD * rng.standard_normal(
                    (n_pool // 4, len(self.free)))
            pool = np.vstack([pool, np.clip(local, 0.0, 1.0)])
        evaluated = self._evaluated_units()
        duplicate = self.search["duplicate_tolerance"]
        cost = self._cost_factor(pool) if self.search["cost_aware"] else None
        chosen, values, lambdas = [], [], None

        for _ in range(int(count or self.search["batch_size"])):
            if self.search["mode"] == "pareto":
                lambdas = rng.dirichlet(np.ones(len(self.objective_targets)))
            incumbent = self._incumbent(normals, lambdas)
            acquisition = self._acquisition(pool, normals, incumbent, lambdas)
            if cost is not None:  # improvement per predicted second
                acquisition = np.where(
                    acquisition > 0,
                    acquisition / cost**self.search["cost_exponent"],
                    acquisition)
            spread = np.ptp(acquisition) if np.ptp(acquisition) > 0 else 1.0
            for point in chosen:  # local penalization always lowers nearby values
                closeness = np.exp(
                    -0.5 *
                    (np.linalg.norm(pool - point, axis=1) / PENALTY_RADIUS)**2)
                acquisition = acquisition - closeness * (np.abs(acquisition) +
                                                         spread)
            for point in [*evaluated, *chosen]:  # never recalculate a point
                acquisition = np.where(
                    np.linalg.norm(pool - point, axis=1) < duplicate, -np.inf,
                    acquisition)
            if not np.isfinite(acquisition).any():
                break
            index = int(np.argmax(acquisition))
            point, value = pool[index], float(acquisition[index])
            if (self.search["local_refinement"] and cost is None and
                    self._predict(point[None]) is not None):

                def negative(x):
                    """Negated acquisition at ``x`` (clipped to the unit cube), for minimization."""
                    return -float(
                        self._acquisition(
                            np.clip(x, 0, 1)[None], normals, incumbent,
                            lambdas)[0])

                refined = minimize(negative,
                                   point,
                                   method="L-BFGS-B",
                                   bounds=[(0.0, 1.0)] * len(self.free),
                                   options={"maxiter": 50})
                far_enough = all(
                    np.linalg.norm(refined.x - other) >= duplicate
                    for other in [*evaluated, *chosen])
                if refined.success and -refined.fun > value and far_enough:
                    point, value = np.clip(refined.x, 0.0,
                                           1.0), float(-refined.fun)
            chosen.append(point)
            values.append(value)
        return [self._from_unit(point) for point in chosen], values

    # -------------------------------------------------------------- screening

    def _screen(self, step, proposals, acquisition):
        """Calculate candidates cheaply; keep the best ``batch_size`` for the full run."""
        screening = self.search["screening"]
        ids = [f"screen-{step:03d}-{i:03d}" for i in range(len(proposals))]
        records = []
        options = dict(screening.get("calculator_options", {}))
        if screening.get("framework_protocol"):
            options["framework_protocol"] = screening["framework_protocol"]
        for parameters, value, inputs in zip(proposals, acquisition,
                                             self.make_inputs(proposals, ids)):
            calculation = self._calculate(inputs, parameters, **options)
            evaluation = self.evaluator.evaluate(calculation, self.targets)
            records.append({
                "parameters": parameters.to_dict(),
                "acquisition": value,
                "values": dict(calculation.values),
                "failed": calculation.failed,
                "failure_reason": calculation.failure_reason,
                "objective": evaluation.objective,
                "feasible": bool(evaluation.feasible and not evaluation.failed)
            })

        def rank(i):
            record = records[i]
            return (record["failed"], not record["feasible"],
                    record["objective"])

        order = sorted(range(len(records)), key=rank)
        keep = sorted(order[:int(self.search["batch_size"])])
        for i, record in enumerate(records):
            record["selected"] = i in keep
        selected = [proposals[i] for i in keep]
        return selected, [acquisition[i] for i in keep], records

    @staticmethod
    def _link_screening(records, sample_ids):
        """Attach full-calculation sample ids, in order, to the selected records."""
        for record, sample_id in zip([r for r in records if r["selected"]],
                                     sample_ids):
            record["sample_id"] = sample_id
        return records

    def screening_agreement(self) -> dict[str, Any] | None:
        """
        Gets the correlation between screening and full values for doubly
        calculated points.

        :return:
            The agreement per target, or None without screening.
        """

        if self.search["screening"] is None:
            return None
        full = {
            r.sample_id: r
            for r in self.table.samples
            if not r.calculation.failed
        }
        pairs = {t.name: [] for t in self.targets}
        screened = 0
        for entry in self.history:
            for record in entry.get("screening") or []:
                screened += 1
                row = full.get(record.get("sample_id"))
                if row is None or record["failed"]:
                    continue
                for name in pairs:
                    if name in record[
                            "values"] and name in row.calculation.values:
                        pairs[name].append((record["values"][name],
                                            row.calculation.values[name]))
        agreement = {}
        for name, values in pairs.items():
            entry = {"n_pairs": len(values)}
            if len(values) >= 3:
                low, high = np.array(values, dtype=float).T
                if np.ptp(low) > 0 and np.ptp(high) > 0:
                    from scipy.stats import spearmanr
                    entry["pearson"] = float(np.corrcoef(low, high)[0, 1])
                    entry["spearman"] = float(spearmanr(low, high).statistic)
                entry["mean_difference"] = float(np.mean(high - low))
            agreement[name] = entry
        return {
            "n_screened": screened,
            "agreement": agreement,
            "used_in_surrogate": False,
            "note": "screening values only rank candidates; they never train the surrogate"
        }

    # ----------------------------------------------------------------- result

    def pareto_front(self) -> list[dict[str, Any]]:
        """
        Finds the non-dominated usable observations in per-target normalized
        squared error.

        Feasible observations are preferred when any exist.

        :return:
            The front, one dictionary per observation.
        """

        scores = self.table.rescore(self.evaluator, self.targets)
        rows = [(r, scores[r.sample_id])
                for r in self.table.samples
                if r.sample_id in scores and not scores[r.sample_id].failed]
        rows = [(r, e) for r, e in rows if e.feasible] or rows
        names = [t.name for t in self.objective_targets]
        vectors = [
            np.array([e.errors[name]**2 for name in names]) for _, e in rows
        ]
        front = []
        for i, (row, evaluation) in enumerate(rows):
            dominated = any(
                np.all(other <= vectors[i]) and np.any(other < vectors[i])
                for j, other in enumerate(vectors)
                if j != i)
            if not dominated:
                front.append({
                    "sample_id": row.sample_id,
                    "parameters": row.parameters.to_dict(),
                    "squared_errors": dict(zip(names, vectors[i].tolist())),
                    "feasible": evaluation.feasible
                })
        return front

    def _build_result(self, best, evaluation, converged,
                      stop_reason) -> FitResult:
        """Base result plus bounds, surrogate-at-best, Pareto, screening diagnostics."""
        result = super()._build_result(best, evaluation, converged, stop_reason)
        diagnostics = result.diagnostics
        diagnostics["strategy"] = "bayesian"
        diagnostics["search"] = dict(self.search)
        diagnostics["n_evaluations"] = len(self.table.samples)
        diagnostics["final_bounds"] = {
            k: list(v) for k, v in self.current_bounds().items()
        }
        diagnostics["bound_expansions"] = list(self.bound_expansions)
        if best is not None:
            predicted = self._predict(self._to_unit(best.parameters)[None])
            if predicted is not None:
                diagnostics["surrogate_at_best"] = {
                    t.name: {
                        "mean": float(predicted[0][0, i]),
                        "std": float(predicted[1][0, i])
                    } for i, t in enumerate(self.objective_targets)
                }
        if self.search["mode"] == "pareto":
            diagnostics["pareto_front"] = self.pareto_front()
        if self.search["screening"] is not None:
            diagnostics["screening"] = self.screening_agreement()
        diagnostics["failures"] = [{
            "sample_id": r.sample_id,
            "reason": r.calculation.failure_reason
        } for r in self.table.samples if r.calculation.failed]
        return result

    # ------------------------------------------------------------ checkpoints

    def _signature(self) -> dict[str, Any]:
        """Settings that must match to resume, using the bounds before any widening."""
        return {
            **super()._signature(),
            # bounds before any widening
            "specs": [asdict(s) for s in self.initial_specs],
            "search": dict(self.search),
        }

    def _resume(self) -> None:
        """Restore samples, history, bound widenings and latest bounds from the checkpoint."""
        self._load_checkpoint(
            "checkpoint parameter/target/search settings differ from this study"
        )
        self.bound_expansions = [
            c for entry in self.history
            for c in entry.get("bound_expansions") or []
        ]
        if self.history and self.history[-1].get("bounds"):
            self._set_bounds(self.history[-1]["bounds"])

    def _print_setup(self) -> None:
        """Print parameters, targets and every search setting to ``ostream``."""
        rows = []
        for spec in self.specs:
            text = (
                f"[{format_value(spec.bounds[0])}, {format_value(spec.bounds[1])}]"
                if spec.optimizable else
                f"fixed at {format_value(spec.initial)}")
            rows.append((f"Parameter {spec.name}", text))
        rows += [(
            f"Target {t.name}", f"{format_value(t.reference)} +/- "
            f"{format_value(t.tolerance)}, {t.kind}, weight {format_value(t.weight)}"
        ) for t in self.targets]
        rows += [(f"Search {key}", value) for key, value in self.search.items()]
        rows.append(("Checkpoint", self._checkpoint_label()))
        print_settings(self.ostream, "CDAM Bayesian Fitting Setup", rows)
