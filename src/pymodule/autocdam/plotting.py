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
"""Standalone plotting for CDAM fits and real-atom metal-complex LJ fits.

``plot_fit`` (also ``FitResult.plot``) accepts a FitResult, a result
dictionary, or a saved result JSON path and picks the figure: search
diagnostics, CDAM stage history, or the metal-complex QM/MM scan. Plotting
never runs QM, MM, MD, or fitting. Matplotlib is imported only when plotting.
It returns a Figure and optionally saves it; it does not call show().
"""
from __future__ import annotations

from collections.abc import Mapping
import importlib.util
import json
from pathlib import Path
import textwrap

import numpy as np

from ..errorhandler import assert_msg_critical

BLUE = "#236A9A"
ORANGE = "#C96B28"
GREY = "#62666B"
COLORS = [BLUE, ORANGE, "#788244", "#9A5781", "#49586B"]


def _plt():
    """Import and return ``matplotlib.pyplot``."""
    assert_msg_critical(
        importlib.util.find_spec("matplotlib") is not None,
        "AutoCDAM plotting: matplotlib is required")
    import matplotlib.pyplot as plt

    return plt


def _result(result) -> dict:
    """Normalize a FitResult, dict or run file to a result dictionary."""
    if isinstance(result, (str, Path)):
        from .models import RESULT_GROUP, read_h5_group
        result = read_h5_group(result, RESULT_GROUP)
        assert_msg_critical(result is not None,
                            "AutoCDAM plotting: no result in the run file")
    if isinstance(result, Mapping):
        # Batch example summaries wrap the shared FitResult under 'result'.
        return result.get("result", result)
    return {
        name: getattr(result, name)
        for name in ("history", "diagnostics", "parameters", "converged",
                     "warnings")
    }


def _dataset(dataset):
    """Normalize a QMScanDataset, dict or JSON path to ``(metadata, record dicts)``."""
    if isinstance(dataset, (str, Path)):
        dataset = json.loads(Path(dataset).read_text())
    if isinstance(dataset, Mapping):
        return dataset.get("metadata", {}), dataset["records"]
    if dataset is not None:
        from dataclasses import asdict
        return dataset.metadata, [asdict(record) for record in dataset.records]
    return {}, []


def _style(axis):
    """Apply the shared axis style: white background, full frame, inward ticks, no grid."""
    axis.set_facecolor("white")
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("#27272A")
    axis.tick_params(labelsize=10,
                     direction="in",
                     top=True,
                     right=True,
                     width=0.8)
    axis.grid(False, which="both")


def _finish(figure,
            data,
            title,
            output_path,
            dpi,
            *,
            note="",
            details=False,
            summary=""):
    """Add the title (and with ``details`` a status line and footer), lay out, and optionally save.

    The footer holds ``note`` and the first recorded warning. Returns the figure.
    """
    figure.set_facecolor("white")
    if details:
        figure.set_size_inches(max(11, figure.get_figwidth()),
                               max(7.3, figure.get_figheight()))
    figure.suptitle(title, fontsize=13, fontweight="normal", y=0.98)
    if not details:
        figure.tight_layout(rect=(0.015, 0.015, 0.99, 0.94),
                            h_pad=1.4,
                            w_pad=1.8)
    else:
        _status_and_footer(figure, data, note, summary)
    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path,
                       dpi=dpi,
                       facecolor="white",
                       **({} if details else {
                           "bbox_inches": "tight"
                       }))
    return figure


def _status_and_footer(figure, data, note, summary):
    """Status line (convergence, stop reason, ``summary``) and footer (``note``, first warning)."""
    diag = data.get("diagnostics", {})
    state = data.get("converged")
    status = "Saved stage history" if state is None else f"Converged: {'yes' if state else 'no'}"
    if diag.get("stop_reason"):
        status += "  |  Stop: " + str(diag["stop_reason"])
    figure.text(0.08,
                1 - 0.52 / figure.get_figheight(),
                status.replace("_", " "),
                fontsize=10,
                color=GREY,
                va="top")
    if summary:
        figure.text(0.08,
                    1 - 0.86 / figure.get_figheight(),
                    summary,
                    fontsize=10,
                    color=GREY,
                    va="top")
    warnings = data.get("warnings", [])
    footer = note
    if warnings:
        footer += ("\n" if footer else "") + "Warning: " + str(warnings[0])
    if footer:
        footer = "\n".join(
            textwrap.fill(line, 115) for line in footer.splitlines())
        figure.text(0.08, 0.025, footer, fontsize=9, color=GREY, va="bottom")
    figure.tight_layout(rect=(0.02, 0.12, 0.99, 0.88), h_pad=2.0, w_pad=2.5)


def _draw_history(axis, data, xlabel):
    """Draw the recorded objective per history entry on ``axis``."""
    history = data.get("history", [])
    x = [entry.get("iteration", i + 1) for i, entry in enumerate(history)]
    y = [
        np.nan if entry.get("objective") is None else float(entry["objective"])
        for entry in history
    ]
    axis.plot(x, y, marker="o", markersize=4, linewidth=1.5, color=BLUE)
    if not history:
        axis.text(0.5,
                  0.5,
                  "No recorded objective history",
                  transform=axis.transAxes,
                  ha="center")
    if x and len(x) <= 12:
        axis.set_xticks(x)
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Weighted normalized objective")
    axis.set_title("Recorded objective", fontsize=12, loc="left")
    _style(axis)


def _plot_cdam_fit(data,
                   *,
                   output_path=None,
                   title=None,
                   dpi=200,
                   details=False):
    """Plot CDAM stage objectives and each target's tolerance-normalized error.

    Target errors and objectives are those actually recorded in stage history;
    raw properties with different physical units are not put on a shared scale.
    """
    history = data.get("history", [])
    names = sorted(
        {name for entry in history for name in entry.get("errors", {})})
    figure, axes = _plt().subplots(1,
                                   2 if names else 1,
                                   figsize=(11, 4.9) if names else (9, 4.7),
                                   squeeze=False)
    _draw_history(axes[0, 0], data, "Fitting stage")
    if names:
        axis = axes[0, 1]
        stages = [
            entry.get("iteration", i + 1) for i, entry in enumerate(history)
        ]
        for i, name in enumerate(names):
            values = [
                entry.get("errors", {}).get(name, np.nan) for entry in history
            ]
            axis.plot(stages,
                      values,
                      marker=["o", "s", "^", "D", "v"][i % 5],
                      linestyle="-" if i < 5 else "--",
                      color=COLORS[i % 5],
                      label=name,
                      linewidth=1.5,
                      markersize=4)
        axis.axhline(1,
                     color=GREY,
                     linestyle="--",
                     linewidth=1,
                     label="Tolerance = 1")
        axis.set_ylim(bottom=0)
        if len(stages) <= 12:
            axis.set_xticks(stages)
        axis.set_xlabel("Fitting stage")
        axis.set_ylabel("Error / target tolerance")
        axis.set_title("Recorded target errors", fontsize=12, loc="left")
        axis.legend(frameon=False, fontsize=9)
        _style(axis)
    note = "Calculated stage results; no surrogate predictions are shown. Numerical convergence is separate from scientific validation."
    return _finish(figure,
                   data,
                   title or "CDAM fitting",
                   output_path,
                   dpi,
                   note=note,
                   details=details)


def _evaluation_rows(data):
    """Per-candidate evaluation rows from the diagnostics; raises ValueError if absent."""
    rows = data.get("diagnostics", {}).get("evaluations")
    if not rows:
        raise ValueError(
            "this result has no per-candidate evaluations; refit with the current "
            "AutoCDAM to record them")
    return rows


def _draw_convergence(axis, data):
    """Draw each candidate's loss (feasible, infeasible, failed) and the running best."""
    rows = _evaluation_rows(data)
    x = np.arange(1, len(rows) + 1)
    objectives = np.array(
        [np.nan if r["objective"] is None else r["objective"] for r in rows])
    feasible = np.array([r["feasible"] for r in rows])
    failed = np.array(
        [r["failed"] or r["quality_status"] != "usable" for r in rows])
    ok = ~failed
    best, running = np.inf, []
    for objective, usable in zip(objectives, ok):
        if usable and np.isfinite(objective):
            best = min(best, objective)
        running.append(best if np.isfinite(best) else np.nan)
    axis.plot(x[ok & feasible],
              objectives[ok & feasible],
              "o",
              color=BLUE,
              markersize=5,
              label="Feasible")
    axis.plot(x[ok & ~feasible],
              objectives[ok & ~feasible],
              "o",
              markerfacecolor="white",
              color=ORANGE,
              markersize=5,
              label="Violates a required target")
    finite = objectives[ok & np.isfinite(objectives)]
    top = float(np.nanmax(finite)) if finite.size else 1.0
    if failed.any():
        axis.plot(x[failed],
                  np.full(failed.sum(), top),
                  "x",
                  color=GREY,
                  markersize=6,
                  label="Failed / not usable")
    axis.step(x,
              running,
              where="post",
              color="#27272A",
              linewidth=1.2,
              label="Best so far")
    if finite.size and np.nanmin(finite) > 0 and top / max(
            np.nanmin(finite), 1e-12) > 50:
        axis.set_yscale("log")
    axis.set_xlabel("Calculated candidate")
    axis.set_ylabel("Weighted normalized loss")
    axis.set_title("Optimization convergence", fontsize=12, loc="left")
    axis.legend(frameon=False, fontsize=8)
    _style(axis)


def _draw_residuals(axis, data):
    """Draw signed (calculated - target) / tolerance bars for the best candidate's targets."""
    diag = data.get("diagnostics", {})
    tolerances = {o["name"]: o["tolerance"] for o in diag.get("objectives", [])}
    bars = []  # (label, signed residual / tolerance, style)
    for name, row in diag.get("target_results", {}).items():
        try:
            value = (float(row["prediction"]) -
                     float(row["reference"])) / tolerances.get(name, 1.0)
        except (TypeError,
                ValueError):  # vector targets: unsigned normalized error
            value = row["normalized_error"]
        bars.append(
            (name, value,
             "required" if row.get("kind") == "required" else "weighted"))
    for name, row in diag.get("monitor", {}).items():
        if row.get("value") is not None:
            value = (float(row["value"]) -
                     float(row["reference"])) / tolerances.get(name, 1.0)
            bars.append((f"{name} (monitor)", value, "monitor"))
    if not bars:
        axis.text(0.5,
                  0.5,
                  "No evaluated best candidate",
                  transform=axis.transAxes,
                  ha="center")
        _style(axis)
        return
    styles = {
        "required": {
            "color": ORANGE
        },
        "weighted": {
            "color": BLUE
        },
        "monitor": {
            "color": "white",
            "edgecolor": GREY,
            "hatch": "//"
        }
    }
    axis.axvspan(-1, 1, color="#E5E7EB", zorder=0)
    for i, (_, value, style) in enumerate(bars):
        axis.barh(i, value, height=0.6, **styles[style])
    axis.axvline(0, color="#27272A", linewidth=0.8)
    axis.set_yticks(range(len(bars)))
    axis.set_yticklabels([label for label, _, _ in bars], fontsize=9)
    axis.set_xlabel(
        "(calculated - target) / tolerance; grey band = within tolerance")
    axis.set_title("Target residuals at best (orange: required)",
                   fontsize=12,
                   loc="left")
    _style(axis)


def _draw_bounds(axis, data):
    """Draw sampled and best parameter values normalized to their bounds."""
    diag = data.get("diagnostics", {})
    rows = _evaluation_rows(data)
    specs = [
        p for p in diag.get("parameter_definitions", []) if p["adjustable"]
    ]
    best = data.get("parameters", {})
    best = best.to_dict() if hasattr(best, "to_dict") else dict(best)
    for i, spec in enumerate(specs):
        lower, upper = spec["bounds"]
        sampled = [(r["parameters"][spec["name"]] - lower) / (upper - lower)
                   for r in rows]
        axis.plot(sampled,
                  np.full(len(sampled), i),
                  "|",
                  color=GREY,
                  markersize=10,
                  alpha=0.6)
        value = (best[spec["name"]] - lower) / (upper - lower)
        near = min(value, 1 - value) <= 0.02
        axis.plot([value], [i],
                  "D",
                  color=ORANGE if near else BLUE,
                  markersize=8)
        axis.text(1.02,
                  i,
                  f"{best[spec['name']]:.4g}",
                  va="center",
                  fontsize=8,
                  transform=axis.get_yaxis_transform())
    axis.set_yticks(range(len(specs)))
    axis.set_yticklabels([spec["name"] for spec in specs], fontsize=9)
    axis.set_xlim(-0.03, 1.03)
    axis.set_xticks([0, 0.5, 1])
    axis.set_xticklabels(["lower", "mid", "upper"])
    axis.set_title("Parameter bounds (orange: at a bound)",
                   fontsize=12,
                   loc="left")
    _style(axis)


def _draw_tradeoff(axis, data):
    """Draw squared errors of the first two Pareto targets for all candidates and the front."""
    diag = data.get("diagnostics", {})
    front = diag.get("pareto_front") or []
    names = list(front[0]["squared_errors"])[:2] if front else []
    rows = [r for r in _evaluation_rows(data) if not r["failed"]]
    if len(names) < 2:
        axis.text(0.5,
                  0.5,
                  "No trade-off recorded",
                  transform=axis.transAxes,
                  ha="center")
        _style(axis)
        return
    points = np.array(
        [[r["errors"].get(n, np.nan)**2 for n in names] for r in rows])
    axis.plot(points[:, 0],
              points[:, 1],
              "o",
              color=GREY,
              markersize=4,
              alpha=0.6,
              label="Calculated")
    front_points = np.array(
        sorted([[f["squared_errors"][n] for n in names] for f in front]))
    axis.plot(front_points[:, 0],
              front_points[:, 1],
              "o-",
              color=ORANGE,
              markersize=6,
              label="Non-dominated")
    axis.set_xlabel(f"{names[0]} squared normalized error")
    axis.set_ylabel(f"{names[1]} squared normalized error")
    axis.set_title("Trade-off", fontsize=12, loc="left")
    axis.legend(frameon=False, fontsize=8)
    _style(axis)


def _plot_search_diagnostics(data,
                             *,
                             output_path=None,
                             title=None,
                             dpi=200,
                             details=False):
    """Convergence, residuals, parameter bounds, and (if recorded) the trade-off."""
    tradeoff = bool(data.get("diagnostics", {}).get("pareto_front"))
    plt = _plt()
    figure, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    _draw_convergence(axes[0, 0], data)
    _draw_residuals(axes[0, 1], data)
    _draw_bounds(axes[1, 0], data)
    if tradeoff:
        _draw_tradeoff(axes[1, 1], data)
    else:
        _draw_history(axes[1, 1], data, "Search step")
    note = (
        "Calculated results only. Optimization convergence is not physical validation of "
        "the force field.")
    return _finish(figure,
                   data,
                   title or "CDAM fitting diagnostics",
                   output_path,
                   dpi,
                   note=note,
                   details=details)


def _plot_metal_complex_fit(data,
                            *,
                            dataset=None,
                            zoom_range=None,
                            output_path=None,
                            title=None,
                            dpi=200,
                            show_residuals=False,
                            show_rejected=False,
                            details=False):
    """Plot QM/MM whole-complex interaction energies and MM-minus-QM residuals.

    Scan distances and directions come from the result's diagnostics; curves
    are separated by approach direction and lines only connect sampled
    distances. ``zoom_range=(lo, hi)`` adds an expanded energy panel in
    angstrom. The default is a clean framed fit: QM points and the fitted MM
    curve. Enable ``show_residuals``, ``show_rejected`` (rejected points of the
    QM ``dataset``) or ``details`` for diagnostics.
    """
    diag = data.get("diagnostics", {})
    ids = diag["record_ids"]
    qm = np.asarray(diag["qm_interaction_kjmol"], dtype=float)
    mm = np.asarray(diag["mm_interaction_kjmol"], dtype=float)
    metadata, records = _dataset(dataset)
    distances = np.array(diag["scan_distances_angstrom"], dtype=float)
    directions = diag["scan_directions"]
    if len(ids) != len(distances) or len(ids) != len(qm) or len(ids) != len(
            mm) or len(ids) != len(directions):
        raise ValueError(
            "scan coordinates, QM energies, and MM energies must have equal lengths"
        )
    if not len(ids) or not np.isfinite(np.r_[distances, qm, mm]).all():
        raise ValueError("plotting requires finite, nonempty scan data")
    if zoom_range is not None:
        lo, hi = zoom_range
        if not np.isfinite([lo, hi]).all() or lo >= hi:
            raise ValueError(
                "zoom_range must be finite increasing distances in angstrom")

    width = 10 if zoom_range is not None else 6.5
    figure = _plt().figure(figsize=(width, 6.4 if show_residuals else 4.5))
    grid = figure.add_gridspec(
        2 if show_residuals else 1,
        2,
        height_ratios=[1.5, 1] if show_residuals else [1])
    overview = figure.add_subplot(grid[
        0, 0] if zoom_range is not None else grid[0, :])
    expanded = figure.add_subplot(grid[0,
                                       1]) if zoom_range is not None else None
    residual = figure.add_subplot(grid[1, :]) if show_residuals else None
    groups = {}
    for i, direction in enumerate(directions):
        groups.setdefault(tuple(np.round(direction, 8)), []).append(i)
    for group_index, indices in enumerate(groups.values()):
        indices = np.array(indices)[np.argsort(distances[indices])]
        suffix = f", direction {group_index+1}" if len(groups) > 1 else ""
        marker = ["o", "s", "^", "D", "v"][group_index % 5]
        for axis in [overview] + ([expanded] if expanded is not None else []):
            axis.plot(distances[indices],
                      qm[indices],
                      color=BLUE,
                      marker=marker,
                      markersize=5,
                      markerfacecolor="white",
                      linestyle="none",
                      label="QM" + suffix)
            axis.plot(distances[indices],
                      mm[indices],
                      color=ORANGE,
                      linestyle=["-", "--", "-.", ":"][group_index % 4],
                      linewidth=1.6,
                      label="Fitted MM" + suffix)
        if residual is not None:
            residual.plot(distances[indices], (mm - qm)[indices],
                          marker=marker,
                          color=ORANGE,
                          linewidth=1.3,
                          markersize=4,
                          label="MM − QM" + suffix)
    for axis in [overview] + ([expanded] if expanded is not None else []):
        if details:
            axis.axhline(0, color=GREY, linewidth=0.6)
        axis.set_xlabel("Metal–probe distance (Å)")
        axis.set_ylabel("Interaction energy (kJ/mol)")
        axis.legend(frameon=False, fontsize=9)
        _style(axis)
    if zoom_range is not None:
        overview.set_title("Full scan", fontsize=11)
    if expanded is not None:
        mask = (distances >= lo) & (distances <= hi)
        if not mask.any():
            _plt().close(figure)
            raise ValueError("zoom_range contains no sampled distances")
        energy = np.r_[qm[mask], mm[mask], 0]
        padding = max(float(np.ptp(energy)) * 0.12, 0.01)
        expanded.set_xlim(lo, hi)
        expanded.set_ylim(np.min(energy) - padding, np.max(energy) + padding)
        expanded.set_title("Expanded energy scale", fontsize=11)
    rejected = [row for row in records if not row.get("accepted", True)]
    if rejected and show_rejected:
        overview.plot([row["distance_angstrom"] for row in rejected],
                      [0.035] * len(rejected),
                      linestyle="none",
                      marker="|",
                      markersize=8,
                      color=GREY,
                      transform=overview.get_xaxis_transform(),
                      label="Rejected point (no energy)")
        overview.legend(frameon=False, fontsize=9)
    if residual is not None:
        tolerance = diag.get("energy_tolerance_kjmol")
        if tolerance is not None and np.isfinite(tolerance) and tolerance > 0:
            residual.axhline(tolerance,
                             color=GREY,
                             linestyle=":",
                             linewidth=0.8,
                             label="Per-point tolerance")
            residual.axhline(-tolerance,
                             color=GREY,
                             linestyle=":",
                             linewidth=0.8)
        residual.axhline(0, color=GREY, linewidth=0.6)
        residual.set_xlabel("Metal–probe distance (Å)")
        residual.set_ylabel("MM − QM (kJ/mol)")
        residual.legend(frameon=False, fontsize=9)
        _style(residual)
    rmse = float(np.sqrt(np.mean((mm - qm)**2)))
    params = data.get("parameters", {})
    if hasattr(params, "to_dict"):
        params = params.to_dict()
    if "values" in params:
        params = params["values"]
    annotation = f"RMSE = {rmse:.3g} kJ/mol   |   {len(ids)} fitted points"
    if "metal_lj_sigma" in params and "metal_lj_epsilon" in params:
        annotation += f"   |   σ = {params['metal_lj_sigma']:.4g} Å   |   ε = {params['metal_lj_epsilon']:.4g} kJ/mol"
    probe = diag.get("probe", metadata.get("study", {}).get("probe"))
    default_title = "QM–MM interaction fit" + (f" · {probe} probe"
                                               if probe else "")
    note = "Whole-complex interaction energies; lines connect sampled points. Numerical agreement alone does not establish scientific validity."
    return _finish(figure,
                   data,
                   title or default_title,
                   output_path,
                   dpi,
                   note=note,
                   details=details,
                   summary=annotation)


def plot_fit(result, *, dataset=None, **kwargs):
    """Plot a FitResult, result dict or saved JSON; returns the Matplotlib Figure.

    Keyword options: ``output_path``, ``title``, ``dpi``, ``details``; for a
    metal-complex fit also ``dataset``, ``zoom_range``, ``show_residuals`` and
    ``show_rejected``.
    """
    data = _result(result)
    if "qm_interaction_kjmol" in data.get("diagnostics", {}):
        return _plot_metal_complex_fit(data, dataset=dataset, **kwargs)
    if data.get("diagnostics", {}).get("evaluations"):
        return _plot_search_diagnostics(data, **kwargs)
    return _plot_cdam_fit(data, **kwargs)
