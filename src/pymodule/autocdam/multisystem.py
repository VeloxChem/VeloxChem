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
"""Route each target to its system (aqueous CDAM or a framework) and merge the results."""

from __future__ import annotations

import time
from dataclasses import replace
from typing import Iterable

from ..outputstream import OutputStream
from .calculator import Calculator
from .framework import FrameworkSystem
from .candidate import CDAMCandidate
from .models import (CalculationResult, FitTarget, GeneratedBatch, ParameterSet,
                     default_ostream, AQUEOUS_PROPERTIES, convert_units)


class MultiSystemCalculator:
    """Calculate one parameter set in several systems and merge the results.

    Each target names its system (``FitTarget.system``) and property
    (``FitTarget.property_name``). Aqueous targets go to the CDAM
    ``Calculator`` with the shared candidate batch; framework targets go to
    their ``FrameworkSystem``, which receives the same CDAM parameter point.
    Values are keyed by target name, so several targets may share a property
    in different systems. A failure in any system fails the whole proposal
    and its reason names the system.
    """

    def __init__(self,
                 systems,
                 aqueous: Calculator | None = None,
                 *,
                 default_system: str | None = None,
                 ostream: OutputStream | None = None):
        """Store the systems (``None`` marks the aqueous one) and the aqueous ``Calculator``."""
        self.systems = dict(systems)
        self.aqueous = aqueous
        self.default_system = default_system
        self.ostream = default_ostream(ostream)
        self.aqueous_names = {
            name for name, system in self.systems.items() if system is None
        }
        self.requires_openmm_systems = False

    def system_of(self, target: FitTarget) -> str:
        """Name of the system a target belongs to; raises ``ValueError`` if missing or unknown."""
        name = target.system or self.default_system
        if name is None:
            raise ValueError(f"target {target.name!r} must name its system "
                             f"({sorted(self.systems)})")
        if name not in self.systems:
            raise ValueError(
                f"target {target.name!r} uses unknown system {name!r}")
        return name

    def check_targets(self, targets) -> None:
        """Validate routing and properties before any calculation starts."""
        targets = list(targets)
        for name, system in self.systems.items():
            group = [t for t in targets if self.system_of(t) == name]
            if system is None:
                for target in group:
                    if target.property_name not in AQUEOUS_PROPERTIES:
                        raise ValueError(
                            f"target {target.name!r}: aqueous systems support "
                            f"{sorted(AQUEOUS_PROPERTIES)}")
            elif group:
                system.check_targets(group)
        self.requires_openmm_systems = any(
            self.systems[self.system_of(t)] is None and
            t.property_name != "QM_ENERGY" for t in targets)

    def calculate(self, batch: GeneratedBatch, parameters: ParameterSet,
                  targets: Iterable[FitTarget], **kwargs) -> CalculationResult:
        """Calculate every system that has targets and merge the results under target names.

        Aqueous targets are deduplicated by property so one aqueous calculation serves
        several targets; it runs in the canonical unit and each target's value is then
        converted to that target's own unit. ``framework_protocol`` is forwarded only to framework systems.
        Per-system raw data records ``runtime_seconds``. Any system failure marks the
        merged result failed, with reasons prefixed by the system name.
        """
        targets = list(targets)
        merged = CalculationResult()
        failures = []
        framework_protocol = kwargs.pop("framework_protocol", None)
        for name, system in self.systems.items():
            group = [t for t in targets if self.system_of(t) == name]
            if not group:
                continue
            started = time.monotonic()
            if system is None:
                # One aqueous calculation per property, in its canonical unit; each
                # target's own unit is applied below, so targets may differ in unit.
                inner = {}
                for target in group:
                    inner.setdefault(
                        target.property_name,
                        replace(target, name=target.property_name, unit=None))
                result = self.aqueous.calculate(batch, parameters,
                                                list(inner.values()), **kwargs)
            else:
                candidate = CDAMCandidate.from_generated_batch(batch)
                result = system.calculate(candidate,
                                          group,
                                          protocol=framework_protocol)
            for target in group:
                key = target.property_name if system is None else target.name
                if key not in result.values:
                    if isinstance(result.metadata.get(key), dict):
                        merged.metadata[target.name] = {
                            **result.metadata[key], "system": name
                        }
                    continue
                values = {key: result.values[key]}
                metadata = {key: dict(result.metadata.get(key) or {})}
                if system is None:
                    convert_units(values, metadata, [replace(target, name=key)])
                merged.values[target.name] = values[key]
                if isinstance(result.metadata.get(key), dict):
                    merged.metadata[target.name] = {
                        **metadata[key], "system": name
                    }
            merged.raw_data[name] = {
                **result.raw_data, "runtime_seconds": time.monotonic() - started
            }
            merged.warnings.extend(
                f"[{name}] {warning}" for warning in result.warnings)
            if result.failed:
                failures.append(f"[{name}] {result.failure_reason}")
        if failures:
            merged.failed = True
            merged.failure_reason = "; ".join(failures)
        return merged


def build_systems(systems, targets, calculator,
                  ostream) -> MultiSystemCalculator:
    """Create the named systems and a router that sends each target to its system.

    ``{"type": "aqueous"}`` is the solvated CDAM ion (calculated by
    ``calculator``); ``"framework"`` (aliases ``cluster``/``molecular``) takes
    ``FrameworkSystem`` keywords. Targets without ``system`` use the only
    system, or the aqueous one.
    """
    built = {}
    for name, spec in systems.items():
        if isinstance(spec, FrameworkSystem):
            built[name] = spec
            continue
        spec = dict(spec)
        kind = spec.pop("type", "framework")
        if kind == "aqueous":
            if spec:
                raise ValueError(
                    f"aqueous system {name!r} takes no options; use "
                    "calculator_options and the candidate's water_model")
            built[name] = None
        elif kind in ("framework", "cluster", "molecular"):
            built[name] = FrameworkSystem(name=name, ostream=ostream, **spec)
        else:
            raise ValueError(
                f"system {name!r}: type must be 'aqueous' or 'framework'")
    aqueous = [name for name, system in built.items() if system is None]
    if len(aqueous) > 1:
        raise ValueError("define at most one aqueous system")
    default = next(
        iter(built)) if len(built) == 1 else (aqueous[0] if aqueous else None)
    if aqueous and not isinstance(calculator, Calculator):
        raise ValueError(
            "an injected calculator cannot be combined with an aqueous system")
    router = MultiSystemCalculator(built,
                                   calculator if aqueous else None,
                                   default_system=default,
                                   ostream=ostream)
    router.check_targets(targets)
    return router
