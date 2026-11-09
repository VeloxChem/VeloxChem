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
"""Per-property cache in front of any calculator, with cost accounting.

With a run file, the cache is the ``autocdam/cache`` group of
``<filename>.h5``: it is loaded when the file exists and rewritten after each
new calculation, so a later run with the same file reuses the properties.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from .models import (CalculationResult, FitTarget, GeneratedBatch, all_finite,
                     h5_filename, read_h5_group, write_h5_groups)

CACHE_GROUP = "autocdam/cache"


def _canonical(value):
    """JSON-stable representation with floats rounded to 12 significant digits."""
    if isinstance(value, dict):
        return {
            str(k): _canonical(v)
            for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_canonical(v) for v in value]
    if isinstance(value, np.ndarray):
        return _canonical(value.tolist())
    if isinstance(value, (float, np.floating)):
        return float(f"{float(value):.12g}")
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, Path):
        return str(value)
    return value


class CachedCalculator:
    """
    Reuses property results for candidates that were already calculated.

    Each target is cached under its name and unit, the calculator options,
    and the candidate model (metal, geometry, water model, charges, every CDAM
    parameter); see ``key``. System settings (framework baseline, protocol) are
    not part of the key, so a persistent cache should only be shared by runs of
    the same study. Only
    missing properties are recalculated. Failed or nonfinite results are never
    cached, so they stay visible as failures. ``path`` makes the cache
    persistent (the run file, see the module docstring). Wall time and
    calculation counts are tracked per system.

    :param calculator:
        The wrapped calculator.
    :param filename:
        The run file of the persistent cache.
    :param write:
        Whether the run file is written (False on MPI ranks other than the
        master).
    """

    def __init__(self,
                 calculator,
                 filename: str | Path | None = None,
                 write: bool = True):
        """
        Wraps the calculator and loads the cache of the run file, if any.

        ``write=False`` (MPI ranks other than the master) never writes the file.
        """

        self.inner = calculator
        self.write = write
        self.entries: dict[str, dict[str, Any]] = {}
        self.path = h5_filename(filename) if filename is not None else None
        if self.path is not None:
            self.entries = load_cache(self.path)
        self.stats = {
            "requests": 0,
            "target_requests": 0,
            "cache_hits": 0,
            "calculations": 0,
            "failed_calculations": 0,
            "wall_time_s": 0.0,
            "per_system": {}
        }

    def __getattr__(self, name):
        """Delegate unknown attributes to the wrapped calculator."""
        return getattr(self.inner, name)

    @staticmethod
    def key(batch: GeneratedBatch, target: FitTarget,
            options: dict[str, Any]) -> str:
        """
        Creates the readable cache key: target name and unit, calculator options,
        candidate model.

        Values are stored in the target's unit, and options such as a cheaper
        screening protocol change the values, so both are part of the key.
        References, tolerances, weights, roles and ``work_dir`` are left out, so
        re-weighting reuses results.

        :param batch:
            The generated batch.
        :param target:
            The target.
        :param options:
            The calculator options.

        :return:
            The cache key.
        """

        model = {k: v for k, v in batch.metadata.items() if k != "sample_id"}
        model.update(batch.system.nonbonded_parameters)
        model.update(batch.system.bonded_parameters)
        options = {k: v for k, v in options.items() if k != "work_dir"}
        identity = _canonical({"options": options, "model": model})
        text = json.dumps(identity, sort_keys=True, default=str)
        return f"{target.name} [{target.unit or ''}] {text}"

    def calculate(self, batch, parameters, targets,
                  **kwargs) -> CalculationResult:
        """
        Returns cached values where available and calculates only the missing
        targets.

        All missing targets go to the inner calculator in one call. Finite fresh values
        are stored and the cache file is rewritten; a failed inner result is returned
        unchanged. Reused entries carry ``cached: True`` in their metadata.

        :param batch:
            The generated batch.
        :param parameters:
            The parameter set.
        :param targets:
            The targets.
        :param kwargs:
            The calculator options.

        :return:
            The calculation result.
        """

        targets = list(targets)
        self.stats["requests"] += 1
        self.stats["target_requests"] += len(targets)
        keys = {t.name: self.key(batch, t, kwargs) for t in targets}
        missing = [t for t in targets if keys[t.name] not in self.entries]
        missing_names = {t.name for t in missing}
        result = CalculationResult()
        if missing:
            started = time.monotonic()
            fresh = self.inner.calculate(batch, parameters, missing, **kwargs)
            elapsed = time.monotonic() - started
            self.stats["calculations"] += 1
            self.stats["wall_time_s"] += elapsed
            self._record_systems(fresh, missing, elapsed)
            if fresh.failed:
                self.stats["failed_calculations"] += 1
                return fresh
            for target in missing:
                value = fresh.values.get(target.name)
                if all_finite(value):
                    self.entries[keys[target.name]] = {
                        "value": _canonical(value),
                        "metadata": _canonical(
                            fresh.metadata.get(target.name, {})),
                    }
            self._save()
            if len(missing) == len(targets):
                return fresh  # nothing reused: hand back the calculator's own result
            result.raw_data, result.warnings = dict(fresh.raw_data), list(
                fresh.warnings)
            result.values.update(fresh.values)
            result.metadata.update(fresh.metadata)
        for target in targets:
            if target.name in missing_names:
                continue
            entry = self.entries[keys[target.name]]
            self.stats["cache_hits"] += 1
            result.values[target.name] = entry["value"]
            result.metadata[target.name] = {**entry["metadata"], "cached": True}
        return result

    def _record_systems(self, fresh, targets, elapsed) -> None:
        """Attribute wall time to systems (router timings, else the whole call)."""
        per_system = self.stats["per_system"]
        timings = {
            name: data.get("runtime_seconds")
            for name, data in fresh.raw_data.items()
            if isinstance(data, dict) and
            data.get("runtime_seconds") is not None
        }
        if not timings:
            system_of = getattr(self.inner, "system_of", None)
            if system_of:
                names = {system_of(t) for t in targets}
            else:
                names = {"default"}
            timings = {name: elapsed / len(names) for name in names}
        for name, seconds in timings.items():
            entry = per_system.setdefault(name, {
                "calculations": 0,
                "wall_time_s": 0.0
            })
            entry["calculations"] += 1
            entry["wall_time_s"] += float(seconds)

    def _save(self) -> None:
        """Rewrite the cache group of the run file; no-op without a run file."""
        if self.path is None or not self.write:
            return
        write_h5_groups(
            self.path, {
                CACHE_GROUP: {
                    "keys": list(self.entries),
                    "values": [entry["value"] for entry in self.entries.values()],
                    "metadata":
                        [entry["metadata"] for entry in self.entries.values()],
                }
            })

    def cost_report(self) -> dict[str, Any]:
        """
        Gets the request/hit/calculation counters and wall times, plus cache size
        and path.

        :return:
            A copy of the cost counters.
        """

        report = json.loads(json.dumps(self.stats))
        report["cached_properties"] = len(self.entries)
        report["cache_path"] = str(self.path) if self.path else None
        return report


def load_cache(filename) -> dict[str, dict[str, Any]]:
    """
    Gets the cached properties of a run file.

    :param filename:
        The run file (``.h5`` added when missing).

    :return:
        The entries ``{key: {"value", "metadata"}}``; empty without a cache.
    """

    data = read_h5_group(filename, CACHE_GROUP)
    if not data:
        return {}
    return {
        key: {
            "value": value,
            "metadata": metadata
        } for key, value, metadata in zip(data["keys"], data["values"],
                                          data["metadata"])
    }
