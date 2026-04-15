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

"""
Trust-radius grouping and diagnostics utilities.

This module is intentionally self-contained so it can be integrated in phases.

Integration intent
------------------
This file provides two components intended for
``IMDatabasePointCollecter.determine_trust_radius_gradient``:

1) Improved grouping algorithm:
   - Adaptive geometric cutoff from kNN distance statistics.
   - Mutual-kNN graph construction to reduce single-link chaining.
   - Connected-components basin extraction (same conceptual output shape as the
     legacy grouping, but with additional diagnostics).
   - Automatic fallback to threshold-only grouping if the stricter graph
     disconnects all mixed (dp + reference) basins.

2) Analysis tool for trust-radius optimization quality:
   - Basin-level convergence and boundary diagnostics.
   - Per-datapoint shrinkage, support, and sensitivity metrics.
   - Ranking of problematic datapoints where large radius contraction is paired
     with high local objective sensitivity or weak support.
   - Export helpers for HDF5/CSV.

Minimal integration outline
---------------------------
Inside ``determine_trust_radius_gradient`` (collector):

1. Build ``D_db, D_ref, D_dpref`` exactly as today.
2. Replace ``group_by_connected_components(...)`` with
   ``improved_group_by_connected_components(...)`` from this module.
3. Keep your current basin loop and optimizer calls.
4. Collect optimizer outputs in a list (one item per basin).
5. After optimization, call ``TrustRadiusOptimizationAnalyzer.analyze(...)``.
6. Optionally persist with ``write_trust_radius_report_h5`` and/or
   ``write_trust_radius_report_csv``.

Notes
-----
- The analyzer can run with or without an objective callback.
  If no callback is supplied, sensitivity terms are left as NaN and only
  geometry/support/shrink diagnostics are used.
- All data structures are plain dict/list/ndarray compatible for easy handoff
  to existing code paths.

Worked Example (How To Read The Grouping Output)
------------------------------------------------
The output from ``improved_group_by_connected_components(...)`` has three top
keys: ``groups``, ``diagnostics``, ``labels``.

Given a concrete example:

- ``diagnostics['n_datapoints'] = 9``
- ``diagnostics['n_references'] = 191``
- ``diagnostics['n_total_nodes'] = 200``

Interpretation:

1) Node indexing in the combined graph:
   - Global nodes ``0..8`` are datapoints (9 total).
   - Global nodes ``9..199`` are references (191 total).
   - This comes from block assembly in ``improved_group_by_connected_components``:
     ``D_global = [[D_db, D_dpref], [D_dpref.T, D_ref]]``.

2) Why cutoff is exactly ``1.0`` even though ``quantile_raw`` can be smaller:
   - ``cutoff_info['quantile_raw']`` is the raw adaptive quantile from kNN
     distances.
   - The final cutoff is clipped to ``[min_cutoff_distance, max_cutoff_distance]``.
   - If ``quantile_raw = 0.6749`` and ``min_cutoff_distance = 1.0``, then
     final ``cutoff = 1.0``.

3) Why ``groups`` can include only a subset of datapoints:
   - ``labels`` is the connected-component id for every global node.
   - Example ``labels`` prefix ``[0,1,2,3,4,5,6,7,0,...]`` means:
     datapoint 0 is in component 0, datapoints 1..7 are in their own isolated
     components, datapoint 8 is also in component 0.
   - All references in this example are in component 0.
   - A group is kept only if a component has BOTH:
     at least one datapoint AND at least one reference.
   - Therefore only datapoints ``[0, 8]`` appear in ``groups[0]``.

4) Why count values in ``group_stats`` look like ``2``, ``36290``, ``382``:
   - For one basin with ``n_dp = 2`` and ``n_ref = 191``:
   - ``dp_dp.count`` upper bound (off-diagonal directed entries) is
     ``n_dp*(n_dp-1) = 2*1 = 2``.
   - ``ref_ref.count`` upper bound is ``n_ref*(n_ref-1) = 191*190 = 36290``.
   - ``dp_ref.count`` upper bound is ``n_dp*n_ref = 2*191 = 382``.
   - When counts hit those maxima, it means all those pair distances are finite
     and strictly positive after filtering.

5) How ``adjacency_density`` is formed:
   - The adjacency matrix is binary (1=edge, 0=no edge), symmetric, with
     diagonal removed.
   - Density is computed as:
     ``sum(adjacency) / (n_total_nodes * (n_total_nodes - 1))``.
   - With 200 nodes, denominator is ``200*199 = 39800``.

Use ``explain_grouping_result_for_non_experts(...)`` in this module to generate
a narrative explanation directly from any run output.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import csv
import math

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

try:
    import h5py
except ImportError:
    h5py = None


@dataclass
class GroupingConfig:
    """Configuration for improved grouping."""

    knn_k: int = 6
    adaptive_quantile: float = 0.80
    min_cutoff_distance: float = 1.0
    max_cutoff_distance: float = 6.0
    fixed_cutoff_distance: float | None = None
    include_threshold_edges: bool = True
    require_mutual_knn: bool = True
    fallback_to_threshold_grouping: bool = True
    threshold_fallback_cutoff: float = 3.0


@dataclass
class SensitivityConfig:
    """Configuration for per-datapoint sensitivity analysis."""

    rel_alpha_step: float = 0.10
    abs_alpha_step_min: float = 1.0e-3
    alpha_lower_bound: float = 1.0e-6
    alpha_upper_bound: float = 1.5
    shrink_ratio_problematic: float = 0.50
    support_distance_cutoff: float | None = None
    min_support_for_stability: int = 2
    sensitivity_percentile: float = 0.75


def _finite_knn_values(
    distance_matrix: np.ndarray,
    k: int,
) -> np.ndarray:
    """Collect the k-th nearest finite neighbor distance per row."""

    if distance_matrix.size == 0:
        return np.array([], dtype=np.float64)

    n = distance_matrix.shape[0]
    kth = []
    for i in range(n):
        row = distance_matrix[i]
        finite = row[np.isfinite(row)]
        finite = finite[finite > 0.0]
        if finite.size == 0:
            continue
        finite.sort()
        idx = min(k - 1, finite.size - 1)
        kth.append(float(finite[idx]))
    if len(kth) == 0:
        return np.array([], dtype=np.float64)
    return np.asarray(kth, dtype=np.float64)


def _compute_adaptive_cutoff(
    d_global: np.ndarray,
    cfg: GroupingConfig,
) -> tuple[float, dict[str, Any]]:
    """Compute a robust adaptive cutoff from kNN statistics."""

    if cfg.fixed_cutoff_distance is not None:
        cutoff = float(cfg.fixed_cutoff_distance)
        return cutoff, {
            "mode": "fixed",
            "cutoff": cutoff,
            "knn_k": int(cfg.knn_k),
            "adaptive_quantile": float(cfg.adaptive_quantile),
            "kth_values_size": 0,
        }

    kth = _finite_knn_values(d_global, max(1, int(cfg.knn_k)))
    if kth.size == 0:
        cutoff = float(cfg.threshold_fallback_cutoff)
        return cutoff, {
            "mode": "fallback_empty_knn",
            "cutoff": cutoff,
            "knn_k": int(cfg.knn_k),
            "adaptive_quantile": float(cfg.adaptive_quantile),
            "kth_values_size": 0,
        }

    q = float(np.quantile(kth, np.clip(cfg.adaptive_quantile, 0.05, 0.99)))
    cutoff = float(np.clip(q, cfg.min_cutoff_distance, cfg.max_cutoff_distance))
    return cutoff, {
        "mode": "adaptive",
        "cutoff": cutoff,
        "quantile_raw": q,
        "knn_k": int(cfg.knn_k),
        "adaptive_quantile": float(cfg.adaptive_quantile),
        "kth_values_size": int(kth.size),
        "kth_values_median": float(np.median(kth)),
    }


def _build_neighbor_lists(
    d_global: np.ndarray,
    cutoff: float,
    k: int,
) -> list[np.ndarray]:
    """Build top-k neighbor indices per node under a cutoff."""

    n = d_global.shape[0]
    neighbors: list[np.ndarray] = []
    for i in range(n):
        row = d_global[i]
        valid = np.where(np.isfinite(row) & (row > 0.0) & (row <= cutoff))[0]
        if valid.size == 0:
            neighbors.append(np.array([], dtype=np.int64))
            continue
        order = np.argsort(row[valid])
        chosen = valid[order[: max(1, int(k))]]
        neighbors.append(chosen.astype(np.int64))
    return neighbors


def _adjacency_from_neighbors(
    d_global: np.ndarray,
    neighbors: list[np.ndarray],
    cutoff: float,
    include_threshold_edges: bool,
    require_mutual_knn: bool,
) -> np.ndarray:
    """Create symmetric adjacency from neighbor lists."""

    n = d_global.shape[0]
    adjacency = np.zeros((n, n), dtype=np.int8)
    neighbor_sets = [set(arr.tolist()) for arr in neighbors]

    for i, nbrs in enumerate(neighbors):
        for j in nbrs:
            if require_mutual_knn and (i not in neighbor_sets[j]):
                continue
            adjacency[i, j] = 1
            adjacency[j, i] = 1

    if include_threshold_edges:
        threshold_edges = (d_global <= cutoff).astype(np.int8)
        np.fill_diagonal(threshold_edges, 0)
        adjacency = np.maximum(adjacency, threshold_edges)

    np.fill_diagonal(adjacency, 0)
    return adjacency


def _extract_groups_from_labels(
    labels: np.ndarray,
    n_dp: int,
    n_ref: int,
) -> list[dict[str, Any]]:
    """Extract mixed-component groups with dp/ref index lists."""

    groups: list[dict[str, Any]] = []
    n_components = int(labels.max() + 1) if labels.size > 0 else 0
    ref_global_indices = np.arange(n_dp, n_dp + n_ref)

    for cluster_id in range(n_components):
        refs_in_cluster = [
            int(global_idx - n_dp)
            for global_idx in ref_global_indices
            if labels[global_idx] == cluster_id
        ]
        if len(refs_in_cluster) <= 10:
            continue

        dps_in_cluster = [int(i) for i in range(n_dp) if labels[i] == cluster_id]
        if len(dps_in_cluster) <= 1:
            continue

        groups.append(
            {
                "basin_id": int(cluster_id),
                "datapoint_indices": dps_in_cluster,
                "reference_indices": refs_in_cluster,
            }
        )
    return groups


def _group_diagnostics(
    groups: list[dict[str, Any]],
    d_db: np.ndarray,
    d_ref: np.ndarray,
    d_dpref: np.ndarray,
) -> list[dict[str, Any]]:
    """Compute compact per-basin distance diagnostics."""

    diag = []
    for group in groups:
        dp_idx = np.asarray(group["datapoint_indices"], dtype=np.int64)
        ref_idx = np.asarray(group["reference_indices"], dtype=np.int64)

        d_dp_dp = d_db[np.ix_(dp_idx, dp_idx)] if dp_idx.size else np.array([])
        d_rf_rf = d_ref[np.ix_(ref_idx, ref_idx)] if ref_idx.size else np.array([])
        d_dp_rf = d_dpref[np.ix_(dp_idx, ref_idx)] if (dp_idx.size and ref_idx.size) else np.array([])

        def finite_stats(mat: np.ndarray) -> dict[str, float | int]:
            if mat.size == 0:
                return {"count": 0, "min": math.nan, "median": math.nan, "max": math.nan}
            finite = mat[np.isfinite(mat)]
            finite = finite[finite > 0.0]
            if finite.size == 0:
                return {"count": 0, "min": math.nan, "median": math.nan, "max": math.nan}
            return {
                "count": int(finite.size),
                "min": float(np.min(finite)),
                "median": float(np.median(finite)),
                "max": float(np.max(finite)),
            }

        diag.append(
            {
                "basin_id": int(group["basin_id"]),
                "n_datapoints": int(dp_idx.size),
                "n_references": int(ref_idx.size),
                "dp_dp": finite_stats(d_dp_dp),
                "ref_ref": finite_stats(d_rf_rf),
                "dp_ref": finite_stats(d_dp_rf),
            }
        )
    return diag


def improved_group_by_connected_components(
    d_db: np.ndarray,
    d_ref: np.ndarray,
    d_dpref: np.ndarray,
    config: GroupingConfig | None = None,
) -> dict[str, Any]:
    """
    Improved grouping for trust-radius optimization.

    Returns a dict with:
    - ``groups``: compatible with legacy grouping output.
    - ``diagnostics``: cutoff mode, adjacency density, and per-basin stats.
    - ``labels``: connected-component label per global node.
    """

    cfg = config or GroupingConfig()

    n_dp = int(d_db.shape[0])
    n_ref = int(d_ref.shape[0])
    n_tot = n_dp + n_ref

    d_global = np.block([[d_db, d_dpref], [d_dpref.T, d_ref]])

    cutoff, cutoff_info = _compute_adaptive_cutoff(d_global, cfg)
    neighbors = _build_neighbor_lists(d_global, cutoff, max(1, int(cfg.knn_k)))
    adjacency = _adjacency_from_neighbors(
        d_global,
        neighbors,
        cutoff,
        include_threshold_edges=bool(cfg.include_threshold_edges),
        require_mutual_knn=bool(cfg.require_mutual_knn),
    )
    graph = csr_matrix(adjacency)
    _, labels = connected_components(csgraph=graph, directed=False, return_labels=True)
    groups = _extract_groups_from_labels(labels, n_dp=n_dp, n_ref=n_ref)

    used_fallback = False
    if len(groups) == 0 and cfg.fallback_to_threshold_grouping:
        used_fallback = True
        adjacency_fallback = (d_global <= float(cfg.threshold_fallback_cutoff)).astype(np.int8)
        np.fill_diagonal(adjacency_fallback, 0)
        graph_fb = csr_matrix(adjacency_fallback)
        _, labels = connected_components(csgraph=graph_fb, directed=False, return_labels=True)
        groups = _extract_groups_from_labels(labels, n_dp=n_dp, n_ref=n_ref)
        adjacency = adjacency_fallback
        cutoff = float(cfg.threshold_fallback_cutoff)

    adjacency_density = float(np.sum(adjacency) / max(1, n_tot * (n_tot - 1)))
    diagnostics = {
        "n_datapoints": n_dp,
        "n_references": n_ref,
        "n_total_nodes": n_tot,
        "cutoff": float(cutoff),
        "cutoff_info": cutoff_info,
        "used_fallback": bool(used_fallback),
        "n_groups": int(len(groups)),
        "adjacency_density": adjacency_density,
        "group_stats": _group_diagnostics(groups, d_db, d_ref, d_dpref),
    }
    return {
        "groups": groups,
        "diagnostics": diagnostics,
        "labels": labels,
    }


def _optimizer_result_to_dict(result: Any) -> dict[str, Any]:
    """Convert scipy/minimizer-like result object to serializable dict."""

    if result is None:
        return {
            "success": False,
            "status": -999,
            "message": "no_result",
            "fun": math.nan,
            "nit": -1,
            "nfev": -1,
            "njev": -1,
        }

    data = {
        "success": bool(getattr(result, "success", False)),
        "status": int(getattr(result, "status", -1)),
        "message": str(getattr(result, "message", "")),
        "fun": float(getattr(result, "fun", math.nan)),
        "nit": int(getattr(result, "nit", -1)),
        "nfev": int(getattr(result, "nfev", -1)),
        "njev": int(getattr(result, "njev", -1)),
    }
    return data


class TrustRadiusOptimizationAnalyzer:
    """Analyzer for trust-radius optimization stability and datapoint sensitivity."""

    def __init__(self, config: SensitivityConfig | None = None):
        self.config = config or SensitivityConfig()

    def _support_counts(
        self,
        groups: list[dict[str, Any]],
        d_dpref: np.ndarray,
        n_dp: int,
        support_cutoff: float,
    ) -> np.ndarray:
        support = np.zeros(n_dp, dtype=np.int64)
        for group in groups:
            dp_idx = np.asarray(group["datapoint_indices"], dtype=np.int64)
            ref_idx = np.asarray(group["reference_indices"], dtype=np.int64)
            if dp_idx.size == 0 or ref_idx.size == 0:
                continue
            local = d_dpref[np.ix_(dp_idx, ref_idx)]
            count = np.sum(local <= support_cutoff, axis=1)
            support[dp_idx] = np.maximum(support[dp_idx], count.astype(np.int64))
        return support

    def _per_datapoint_distance_stats(
        self,
        groups: list[dict[str, Any]],
        d_dpref: np.ndarray,
        n_dp: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        min_dist = np.full(n_dp, np.nan, dtype=np.float64)
        mean_dist = np.full(n_dp, np.nan, dtype=np.float64)
        for group in groups:
            dp_idx = np.asarray(group["datapoint_indices"], dtype=np.int64)
            ref_idx = np.asarray(group["reference_indices"], dtype=np.int64)
            if dp_idx.size == 0 or ref_idx.size == 0:
                continue
            local = d_dpref[np.ix_(dp_idx, ref_idx)]
            min_dist[dp_idx] = np.min(local, axis=1)
            mean_dist[dp_idx] = np.mean(local, axis=1)
        return min_dist, mean_dist

    def _sensitivity_scan(
        self,
        initial_alpha: np.ndarray,
        final_alpha: np.ndarray,
        objective_fn: Callable[[np.ndarray], float] | None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Estimate local sensitivity by finite differences around optimized alphas.

        Returns:
            - gradient proxy dL/dalpha_i
            - max perturbation penalty max(L(a+), L(a-) ) - L(a*)
        """

        n = final_alpha.size
        grad = np.full(n, np.nan, dtype=np.float64)
        penalty = np.full(n, np.nan, dtype=np.float64)

        if objective_fn is None or n == 0:
            return grad, penalty

        base = final_alpha.copy()
        base_loss = float(objective_fn(base))

        for i in range(n):
            step = max(
                self.config.abs_alpha_step_min,
                self.config.rel_alpha_step * max(abs(base[i]), abs(initial_alpha[i]), 1.0e-8),
            )
            a_plus = base.copy()
            a_minus = base.copy()
            a_plus[i] = float(np.clip(a_plus[i] + step, self.config.alpha_lower_bound, self.config.alpha_upper_bound))
            a_minus[i] = float(np.clip(a_minus[i] - step, self.config.alpha_lower_bound, self.config.alpha_upper_bound))

            f_plus = float(objective_fn(a_plus))
            f_minus = float(objective_fn(a_minus))
            denom = (a_plus[i] - a_minus[i])
            if abs(denom) > 1.0e-16:
                grad[i] = (f_plus - f_minus) / denom
            penalty[i] = max(f_plus, f_minus) - base_loss

        return grad, penalty

    def analyze(
        self,
        *,
        initial_alphas: np.ndarray,
        final_alphas: np.ndarray,
        groups: list[dict[str, Any]],
        d_dpref: np.ndarray,
        grouping_diagnostics: dict[str, Any] | None = None,
        basin_optimizer_results: list[dict[str, Any]] | None = None,
        datapoint_labels: list[str] | None = None,
        objective_fn: Callable[[np.ndarray], float] | None = None,
    ) -> dict[str, Any]:
        """
        Produce basin-level and datapoint-level diagnostics.

        Parameters
        ----------
        objective_fn:
            Optional callable that accepts full alpha vector and returns scalar
            objective. If omitted, sensitivity fields are NaN.
        """

        a0 = np.asarray(initial_alphas, dtype=np.float64).reshape(-1)
        a1 = np.asarray(final_alphas, dtype=np.float64).reshape(-1)
        if a0.shape != a1.shape:
            raise ValueError("initial_alphas and final_alphas must have identical shape.")

        n_dp = int(a1.size)
        if datapoint_labels is None:
            datapoint_labels = [f"dp_{i}" for i in range(n_dp)]
        if len(datapoint_labels) != n_dp:
            raise ValueError("datapoint_labels length mismatch.")

        cutoff = self.config.support_distance_cutoff
        if cutoff is None:
            cutoff = 3.0
            if grouping_diagnostics is not None:
                cutoff = float(grouping_diagnostics.get("cutoff", cutoff))

        support = self._support_counts(groups, d_dpref, n_dp, cutoff)
        min_dist, mean_dist = self._per_datapoint_distance_stats(groups, d_dpref, n_dp)

        sens_grad, sens_penalty = self._sensitivity_scan(a0, a1, objective_fn)

        rel_change = (a1 - a0) / np.maximum(np.abs(a0), 1.0e-12)
        shrink_ratio = a1 / np.maximum(a0, 1.0e-12)
        at_lower = np.isclose(a1, self.config.alpha_lower_bound, rtol=0.0, atol=1.0e-6)
        at_upper = np.isclose(a1, self.config.alpha_upper_bound, rtol=0.0, atol=1.0e-6)

        finite_pen = sens_penalty[np.isfinite(sens_penalty)]
        if finite_pen.size > 0:
            sensitivity_threshold = float(np.quantile(finite_pen, self.config.sensitivity_percentile))
        else:
            sensitivity_threshold = math.nan

        per_dp = []
        for i in range(n_dp):
            heavy_shrink = bool(shrink_ratio[i] <= self.config.shrink_ratio_problematic)
            low_support = bool(support[i] < self.config.min_support_for_stability)
            high_sensitivity = bool(
                np.isfinite(sens_penalty[i]) and np.isfinite(sensitivity_threshold) and sens_penalty[i] >= sensitivity_threshold
            )
            boundary_issue = bool(at_lower[i] or at_upper[i])
            problematic = bool(heavy_shrink and (high_sensitivity or low_support or boundary_issue))

            issue_score = 0.0
            issue_score += max(0.0, (self.config.shrink_ratio_problematic - shrink_ratio[i]))
            issue_score += 0.5 if low_support else 0.0
            issue_score += 0.5 if boundary_issue else 0.0
            if np.isfinite(sens_penalty[i]) and np.isfinite(sensitivity_threshold) and sensitivity_threshold > 0.0:
                issue_score += max(0.0, sens_penalty[i] / sensitivity_threshold - 1.0)

            per_dp.append(
                {
                    "index": int(i),
                    "label": str(datapoint_labels[i]),
                    "alpha_initial": float(a0[i]),
                    "alpha_final": float(a1[i]),
                    "alpha_delta": float(a1[i] - a0[i]),
                    "alpha_rel_change": float(rel_change[i]),
                    "alpha_shrink_ratio": float(shrink_ratio[i]),
                    "support_count": int(support[i]),
                    "min_dp_ref_distance": float(min_dist[i]) if np.isfinite(min_dist[i]) else math.nan,
                    "mean_dp_ref_distance": float(mean_dist[i]) if np.isfinite(mean_dist[i]) else math.nan,
                    "sensitivity_gradient_proxy": float(sens_grad[i]) if np.isfinite(sens_grad[i]) else math.nan,
                    "sensitivity_penalty": float(sens_penalty[i]) if np.isfinite(sens_penalty[i]) else math.nan,
                    "hits_lower_bound": bool(at_lower[i]),
                    "hits_upper_bound": bool(at_upper[i]),
                    "heavy_shrink": heavy_shrink,
                    "low_support": low_support,
                    "high_sensitivity": high_sensitivity,
                    "problematic": problematic,
                    "issue_score": float(issue_score),
                }
            )

        per_dp.sort(key=lambda row: row["issue_score"], reverse=True)
        problematic = [row for row in per_dp if row["problematic"]]

        basin_results = []
        if basin_optimizer_results is not None:
            for item in basin_optimizer_results:
                basin_id = int(item.get("basin_id", -1))
                result_obj = item.get("result")
                result_dict = _optimizer_result_to_dict(result_obj)
                if result_obj is not None:
                    history = getattr(result_obj, "history", None)
                    if history is not None:
                        result_dict["history"] = history
                        result_dict["history_len"] = int(len(history))
                    comp = getattr(result_obj, "component_summary", None)
                    if comp is not None:
                        result_dict["component_summary"] = comp
                basin_results.append(
                    {
                        "basin_id": basin_id,
                        "n_datapoints": int(item.get("n_datapoints", 0)),
                        "n_references": int(item.get("n_references", 0)),
                        "result": result_dict,
                    }
                )

        summary = {
            "n_datapoints": n_dp,
            "n_problematic": int(len(problematic)),
            "problematic_fraction": float(len(problematic) / max(1, n_dp)),
            "mean_alpha_initial": float(np.mean(a0)) if n_dp > 0 else math.nan,
            "mean_alpha_final": float(np.mean(a1)) if n_dp > 0 else math.nan,
            "median_alpha_shrink_ratio": float(np.median(shrink_ratio)) if n_dp > 0 else math.nan,
            "sensitivity_penalty_threshold": sensitivity_threshold,
        }

        return {
            "summary": summary,
            "grouping_diagnostics": grouping_diagnostics or {},
            "basin_optimizer_results": basin_results,
            "datapoints_ranked": per_dp,
            "problematic_datapoints": problematic,
        }


def write_trust_radius_report_h5(
    report: dict[str, Any],
    h5_path: str | Path,
    *,
    state: int,
    run_id: int | None = None,
) -> None:
    """Persist analysis report in HDF5 under /trust_radius_analysis/state_<state>/run_<id>."""

    if h5py is None:
        raise ImportError(
            "h5py is required for write_trust_radius_report_h5 but is not available."
        )

    path = Path(h5_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(path, "a") as h5:
        root = h5.require_group("trust_radius_analysis")
        state_group = root.require_group(f"state_{int(state)}")

        if run_id is None:
            existing = sorted(
                [
                    int(name.split("_")[-1])
                    for name in state_group.keys()
                    if name.startswith("run_") and name.split("_")[-1].isdigit()
                ]
            )
            run_id = (existing[-1] + 1) if existing else 0

        run_group = state_group.require_group(f"run_{int(run_id)}")

        # Overwrite helper for scalar attributes.
        for key, value in report.get("summary", {}).items():
            run_group.attrs[str(key)] = value

        # Grouping diagnostics as attributes where possible.
        gdiag = run_group.require_group("grouping_diagnostics")
        for key, value in report.get("grouping_diagnostics", {}).items():
            if isinstance(value, (dict, list, tuple)):
                continue
            gdiag.attrs[str(key)] = value

        # Datapoint table as column datasets.
        dp_rows = report.get("datapoints_ranked", [])
        table = run_group.require_group("datapoints_ranked")

        if len(dp_rows) > 0:
            keys = sorted(dp_rows[0].keys())
            n = len(dp_rows)
            for key in keys:
                if key in table:
                    del table[key]

                col = [row.get(key) for row in dp_rows]
                if isinstance(col[0], str):
                    dt = h5py.string_dtype(encoding="utf-8")
                    table.create_dataset(key, data=np.asarray(col, dtype=object), dtype=dt)
                elif isinstance(col[0], bool):
                    table.create_dataset(key, data=np.asarray(col, dtype=np.int8))
                elif isinstance(col[0], int):
                    table.create_dataset(key, data=np.asarray(col, dtype=np.int64))
                else:
                    table.create_dataset(key, data=np.asarray(col, dtype=np.float64))
            table.attrs["n_rows"] = int(n)
        else:
            table.attrs["n_rows"] = 0


def write_trust_radius_report_csv(
    report: dict[str, Any],
    csv_path: str | Path,
) -> None:
    """Write ranked datapoint analysis table to CSV."""

    rows = report.get("datapoints_ranked", [])
    path = Path(csv_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if len(rows) == 0:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["index", "label", "problematic", "issue_score"])
        return

    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_report_for_console(report: dict[str, Any], top_k: int = 10) -> str:
    """Render a compact human-readable summary string."""

    summary = report.get("summary", {})
    rows = report.get("datapoints_ranked", [])
    lines = []
    lines.append(
        "Trust-radius analysis: "
        f"n_dp={summary.get('n_datapoints', 0)}, "
        f"n_problematic={summary.get('n_problematic', 0)}, "
        f"problematic_fraction={summary.get('problematic_fraction', math.nan):.3f}"
    )
    lines.append("Top problematic datapoints:")
    for row in rows[: max(0, int(top_k))]:
        lines.append(
            "  "
            f"{row.get('label')} "
            f"(idx={row.get('index')}, "
            f"R0={row.get('alpha_initial'):.4g}, R*={row.get('alpha_final'):.4g}, "
            f"ratio={row.get('alpha_shrink_ratio'):.3f}, "
            f"support={row.get('support_count')}, "
            f"sens_pen={row.get('sensitivity_penalty')}, "
            f"problematic={row.get('problematic')})"
        )
    return "\n".join(lines)


def _preview_int_list(values: list[int], max_items: int = 12) -> str:
    """Compact preview for index lists."""

    if len(values) <= max_items:
        return str(values)
    head = values[: max_items // 2]
    tail = values[-max(1, max_items // 2):]
    return f"{head} ... {tail} (len={len(values)})"