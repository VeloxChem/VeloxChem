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
"""Post-fit diagnostics, flat result tables, and fitted-model export.

Diagnostics functions return plain dictionaries or lists for
``FitResult.diagnostics``, so results stay plain, storable values.

``samples_table`` accepts a ``FitResult``, a search driver, a
``ParameterPropertyTable``, or a run file (``<filename>.h5``, with the search
checkpoint and the result) and returns one row per candidate::

    sample_id, step, parameter.<name>..., property.<name>...,
    uncertainty.<name>... (e.g. HFE), spread.<name>... (e.g. IOD, CN),
    objective, feasible (when scores are available),
    quality_status, calculation_failed, failure_reason, protocol_id, runtime_seconds

``save_table`` writes the rows as CSV; ``format_table`` renders them as aligned
text. ``export_cdam_model`` writes a fitted CDAM model: candidate PDB/XML,
framework systems and provenance JSON.
"""

from __future__ import annotations

from typing import Any, Sequence, Iterable
import csv
import json
import math
from dataclasses import asdict
from pathlib import Path

try:
    import openmm as mm
except ImportError:
    pass

from .models import (RESULT_GROUP, FitResult, FitTarget, ParameterSpec,
                     h5_filename, read_h5_group)

# ==========================================================================
# reporting
# ==========================================================================

BOUND_FRACTION = 0.02  # a fitted value within this fraction of the width is "at the bound"
PHYSICAL_VALIDATION_NOTE = (
    "not performed: meeting the targets shows the optimizer converged, not that the "
    "force field is transferable or physically validated")


def parameter_definitions(specs: Sequence[ParameterSpec],
                          final_bounds: dict[str, Any]) -> list[dict[str, Any]]:
    """Each parameter with its searched bounds (soft bounds: initial, final, hard)."""
    definitions = []
    for spec in specs:
        definition = {
            "name": spec.name,
            "adjustable": spec.optimizable,
            "bounds": (list(final_bounds.get(spec.name, spec.bounds))
                       if spec.bounds else None),
            "initial" if spec.optimizable else "fixed_value": spec.initial,
            "unit": spec.unit
        }
        if spec.expandable:
            definition.update({
                "initial_bounds": list(spec.bounds),
                "expandable": True,
                "hard_bounds": list(spec.hard_bounds or ()),
                "bounds_source": spec.metadata.get("bounds_source")
            })
        definitions.append(definition)
    return definitions


def objective_definitions(targets: Sequence[FitTarget]) -> list[dict[str, Any]]:
    """Each target with its system, property, reference, tolerance, scale, weight and role."""
    return [{
        "name": t.name,
        "system": t.system,
        "property": t.property_name,
        "reference": t.reference,
        "tolerance": t.tolerance,
        "scale": t.normalization_scale,
        "weight": t.weight,
        "role": t.kind,
        "unit": t.unit
    } for t in targets]


def bound_warnings(specs: Sequence[ParameterSpec], parameters: dict[str, float],
                   final_bounds: dict[str, Any]) -> list[str]:
    """Warnings for fitted values that ended at a (final) bound."""
    warnings = []
    for spec in specs:
        if not spec.optimizable or spec.name not in parameters:
            continue
        lower, upper = final_bounds.get(spec.name, spec.bounds)
        margin = BOUND_FRACTION * (upper - lower)
        value = parameters[spec.name]
        if value - lower <= margin or upper - value <= margin:
            side = "lower" if value - lower <= margin else "upper"
            advice = (
                "soft bound already widened to its limit or the maximum number of "
                "times; check identifiability" if spec.expandable else
                "widen the bounds or check identifiability")
            warnings.append(
                f"WARNING: fitted {spec.name} = {value:.4g} lies at its {side} bound "
                f"[{lower:.4g}, {upper:.4g}]; {advice}")
    return warnings


def evaluation_records(samples, evaluator, targets) -> list[dict[str, Any]]:
    """Every calculated candidate with its values, score, and quality."""
    records = []
    for sample in samples:
        score = evaluator.evaluate(sample.calculation, targets)
        records.append({
            "sample_id": sample.sample_id,
            "parameters": sample.parameters.to_dict(),
            "values": dict(sample.calculation.values),
            "objective": None if score.failed else score.objective,
            "errors": dict(score.errors),
            "feasible": bool(score.feasible and not score.failed),
            "quality_status": sample.quality_status,
            "failed": sample.calculation.failed,
            "failure_reason": sample.calculation.failure_reason,
            "runtime_seconds": sample.runtime_seconds,
        })
    return records


def cost_report(samples, calculator) -> dict[str, Any]:
    """Candidate wall time, plus calculation counts and cache reuse if cached."""
    cost = {
        "n_candidates": len(samples),
        "candidate_wall_time_s": float(
            sum(s.runtime_seconds or 0.0 for s in samples))
    }
    report = getattr(calculator, "cost_report", None)
    if callable(report):
        cost.update(report())
    return cost


def system_report(systems, samples, water_model: str) -> dict[str, Any]:
    """Each system's identity; for frameworks, the baseline and that its surroundings held."""
    report = {}
    for name, system in systems.systems.items():
        if system is None:
            report[name] = {"type": "aqueous", "water_model": water_model}
            continue
        baseline = system.prepare_baseline()
        checked = [
            meta.get("non_target_preserved")
            for sample in samples
            for meta in sample.calculation.metadata.values()
            if isinstance(meta, dict) and meta.get("system") == name
        ]
        report[name] = {
            "type": "framework",
            "baseline_source": baseline.source,
            "baseline_fingerprint": baseline.fingerprint,
            "target_metal": system.labels[system.target_index],
            "representation": system.representation,
            "periodic": baseline.periodic,
            "removed_metal_bonded_terms": dict(baseline.removed_metal_terms),
            "non_target_preserved": bool(checked) and all(checked),
            "n_checked_calculations": len(checked),
        }
    return report


def monitor_report(calculation, evaluation, targets) -> dict[str, Any]:
    """Calculated value, reference and tolerance status of each monitored target."""
    return {
        t.name: {
            "value": calculation.values.get(t.name),
            "reference": t.reference,
            "within_tolerance": evaluation.target_results.get(
                t.name, {}).get("within_tolerance")
        } for t in targets
    }


def add_validation_status(result: FitResult) -> None:
    """Record in place that the optimizer converged (or not) but nothing was physically validated."""
    result.diagnostics["optimization_converged"] = bool(result.converged)
    result.diagnostics["scientifically_validated"] = False
    result.diagnostics["physical_validation"] = PHYSICAL_VALIDATION_NOTE


# ==========================================================================
# tables
# ==========================================================================


def samples_table(source) -> list[dict[str, Any]]:
    """
    Creates rows for every calculated candidate of a source.

    See the module docstring for the columns.

    :param source:
        A FitResult, a search driver, a ParameterPropertyTable, or a run file.

    :return:
        The list of rows.
    """

    samples, scores = _load(source)
    return [_row(s, scores.get(s["sample_id"])) for s in samples]


def _row(sample: dict[str, Any], score: dict | None) -> dict[str, Any]:
    """One flat row; per-property metadata adds uncertainty and spread columns."""
    row = {"sample_id": sample["sample_id"], "step": _step(sample["sample_id"])}
    row.update({f"parameter.{k}": v for k, v in sample["parameters"].items()})
    row.update({f"property.{k}": v for k, v in sample["values"].items()})
    for name, meta in (sample.get("metadata") or {}).items():
        if isinstance(meta, dict):
            for key in ("uncertainty", "spread"):
                if meta.get(key) is not None:
                    row[f"{key}.{name}"] = meta[key]
    if score is not None:
        row["objective"] = score.get("objective")
        row["feasible"] = score.get("feasible")
    row.update({
        "quality_status": sample.get("quality_status"),
        "calculation_failed": sample.get("failed", False),
        "failure_reason": sample.get("failure_reason"),
        "protocol_id": sample.get("protocol_id"),
        "runtime_seconds": sample.get("runtime_seconds")
    })
    return row


def save_table(rows: Iterable[dict[str, Any]], path) -> Path:
    """
    Writes rows as CSV (columns in first-seen order).

    :param rows:
        The rows.
    :param path:
        The CSV file.

    :return:
        The path of the CSV file.
    """

    rows = list(rows)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def format_table(rows: Iterable[dict[str, Any]],
                 columns=None,
                 digits: int = 4) -> str:
    """
    Formats rows as an aligned plain-text table.

    Floats are shown with ``digits`` significant figures, missing cells are blank.

    :param rows:
        The rows.
    :param columns:
        The selected columns, in order (default: all).
    :param digits:
        The number of significant figures.

    :return:
        The table text.
    """

    rows = list(rows)
    if columns is None:
        columns = list(dict.fromkeys(key for row in rows for key in row))
    cells = [[_cell(row.get(c), digits) for c in columns] for row in rows]
    widths = [
        max([len(c)] + [len(r[i])
                        for r in cells])
        for i, c in enumerate(columns)
    ]
    lines = [
        "  ".join(c.ljust(w) for c, w in zip(columns, widths)),
        "  ".join("-" * w for w in widths)
    ]
    lines += ["  ".join(v.rjust(w) for v, w in zip(r, widths)) for r in cells]
    return "\n".join(lines)


def _cell(value, digits: int) -> str:
    """Compact text for one table cell."""
    if value is None:
        return ""
    if isinstance(value, float):
        return "nan" if math.isnan(value) else f"{value:.{digits}g}"
    return str(value)


# ------------------------------------------------------------------ loading


def _step(sample_id) -> int | None:
    """Step number from a ``<prefix>-<step>-<index>`` sample id, else None."""
    parts = str(sample_id).split("-")
    return int(parts[1]) if len(parts) == 3 and parts[1].isdigit() else None


def _flatten(s: dict[str, Any]) -> dict[str, Any]:
    """Plain sample dictionary from a saved (or ``asdict``) ``ParameterPropertySample``."""
    calculation = s["calculation"]
    return {
        "sample_id": s["sample_id"],
        "parameters": s["parameters"]["values"],
        "values": calculation["values"],
        "metadata": calculation["metadata"],
        "quality_status": s["quality_status"],
        "failed": calculation["failed"],
        "failure_reason": calculation["failure_reason"],
        "protocol_id": s["protocol_id"],
        "runtime_seconds": s["runtime_seconds"]
    }


def _from_objects(samples) -> list[dict[str, Any]]:
    """Plain sample dictionaries from in-memory sample objects."""
    return [_flatten(asdict(s)) for s in samples]


def _from_checkpoint(filename, group=None) -> list[dict[str, Any]]:
    """Plain sample dictionaries read from the checkpoint group of a run file."""
    from .search import CHECKPOINT_GROUP
    data = read_h5_group(filename, group or CHECKPOINT_GROUP)
    return [_flatten(s) for s in (data or {}).get("samples", [])]


def _from_evaluations(
        evaluations) -> tuple[list[dict[str, Any]], dict[str, dict]]:
    """Samples and scores from ``FitResult.diagnostics["evaluations"]`` records (no metadata)."""
    samples = [{
        "sample_id": e["sample_id"],
        "parameters": e["parameters"],
        "values": e["values"],
        "quality_status": e.get("quality_status"),
        "failed": e.get("failed", False),
        "failure_reason": e.get("failure_reason"),
        "runtime_seconds": e.get("runtime_seconds")
    } for e in evaluations]
    return samples, _scores(evaluations)


def _scores(evaluations) -> dict[str, dict]:
    """``{sample_id: {"objective", "feasible"}}`` from evaluation records."""
    return {
        e["sample_id"]: {
            "objective": e.get("objective"),
            "feasible": e.get("feasible")
        } for e in evaluations or []
    }


def _load(source) -> tuple[list[dict[str, Any]], dict[str, dict]]:
    """Samples and scores from any accepted source.

    ``source`` is a FitResult, a search table or driver, or a run file. A
    checkpoint is preferred over saved evaluations because it also carries
    per-property metadata. Raises FileNotFoundError for a missing run file and
    ValueError for a run file without samples or result.
    """
    from .models import FitResult
    from .search import ParameterPropertyTable

    if isinstance(source, ParameterPropertyTable):
        return _from_objects(source.samples), {}
    if hasattr(source, "table") and isinstance(source.table,
                                               ParameterPropertyTable):
        return _from_objects(source.table.samples), {}
    if isinstance(source, FitResult):
        evaluations = source.diagnostics.get("evaluations", [])
        if source.checkpoint and h5_filename(source.checkpoint).exists():
            samples = _from_checkpoint(source.checkpoint,
                                       source.checkpoint_group)
            if samples:
                return samples, _scores(evaluations)
        return _from_evaluations(evaluations)

    # a run file: the checkpoint carries per-property metadata, so prefer it
    path = h5_filename(source)
    if not path.exists():
        raise FileNotFoundError(f"no AutoCDAM run file {path}")
    result = read_h5_group(path, RESULT_GROUP)
    evaluations = (result or {}).get("diagnostics", {}).get("evaluations", [])
    samples = _from_checkpoint(path, (result or {}).get("checkpoint_group"))
    if samples:
        return samples, _scores(evaluations)
    if evaluations:
        return _from_evaluations(evaluations)
    raise ValueError(f"no AutoCDAM samples or result in {path}")


# ==========================================================================
# export
# ==========================================================================


def export_cdam_model(directory,
                      candidate,
                      parameters: dict[str, float],
                      systems=None) -> list[Path]:
    """Write ``candidate`` with ``parameters`` into ``directory``; return the written paths.

    Per framework system: the serialized OpenMM system and its coordinates
    (dummies included, angstrom). ``fitted_model.json`` records the provenance.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    point = candidate.with_parameters(parameters)
    written = list(point.write_files(directory, stem="fitted_cdam"))
    provenance = {
        "parameters": parameters,
        "candidate": point.metadata(),
        "scientifically_validated": False,
        "systems": {}
    }
    for name, system in (systems.systems.items() if systems else []):
        if system is None:
            continue
        omm, positions, site = system.build(point)
        xml_path = directory / f"{name}_system.xml"
        xml_path.write_text(mm.XmlSerializer.serialize(omm))
        labels = list(
            system.labels) + [f"D{i + 1}" for i in range(len(site["dummies"]))]
        xyz_path = directory / f"{name}_positions.xyz"
        xyz_path.write_text(
            f"{len(labels)}\n{name} fitted model, angstrom\n" +
            "\n".join(f"{label} {x * 10:.6f} {y * 10:.6f} {z * 10:.6f}"
                      for label, (x, y, z) in zip(labels, positions)) + "\n")
        baseline = system.prepare_baseline()
        provenance["systems"][name] = {
            "baseline_source": baseline.source,
            "baseline_fingerprint": baseline.fingerprint,
            "target_metal": system.labels[system.target_index],
            "representation": system.representation,
            "site": dict(site),
            "non_target_parameters": "baseline values, unchanged"
        }
        written += [xml_path, xyz_path]
    provenance_path = directory / "fitted_model.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, default=str))
    written.append(provenance_path)
    return written
