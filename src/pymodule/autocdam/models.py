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
"""Shared data contracts, units, reference bounds, aqueous properties and printing.

Python objects in this module are the source of truth for the workflow. XML,
PDB, and checkpoint files are output artifacts and must not become the primary
internal representation. Sections, in order:

* printing: VeloxChem-style output streams and centered ``Label : value``
  settings blocks (streams default to ``sys.stdout`` at call time);
* units: CDAM (angstrom, kcal/mol) <-> OpenMM (nm, kJ/mol) conversions;
* reference: UFF Lennard-Jones values (Rappé et al., J. Am. Chem. Soc. 1992,
  114, 10024; table from VeloxChem) that scale wide, soft first bounds. The
  Bayesian search widens a soft side the best candidate reaches, up to the hard
  physical limits; explicit ``(lower, upper)`` bounds stay hard;
* models: parameters, targets, batches and results;
* properties: which job computes each aqueous property, in which canonical
  unit. Adding one means an entry in ``AQUEOUS_PROPERTIES`` plus its
  measurement (``geometry.py`` for MD properties); ``convert_units`` converts
  to the unit a target asks for.
"""

from __future__ import annotations

import math
import sys
import textwrap
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping, Sequence

from mpi4py import MPI

from ..errorhandler import assert_msg_critical
from ..outputstream import OutputStream
from ..veloxchemlib import hartree_in_kcalpermol, hartree_in_kjpermol, mpi_master

# ==========================================================================
# printing
# ==========================================================================

SETTINGS_WIDTH = 84
LABEL_WIDTH = 32


def default_ostream(ostream: OutputStream | None = None,
                    comm=None) -> OutputStream:
    """Return ``ostream``, or the VeloxChem default stream for ``comm``.

    The default writes to the current ``sys.stdout`` on the MPI master rank and
    nowhere on the other ranks.
    """
    if ostream is not None:
        return ostream
    if comm is None:
        comm = MPI.COMM_WORLD
    return OutputStream(
        sys.stdout) if comm.Get_rank() == mpi_master() else OutputStream(None)


def format_value(value) -> str:
    """Format numbers, sequences, and mappings compactly for a settings line."""
    if value is None:
        return "None"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, float):
        if not math.isfinite(value):
            return str(value)
        return f"{value:.6g}" if value == 0.0 or 1e-4 <= abs(
            value) < 1e6 else f"{value:.4e}"
    if hasattr(value, "item") and getattr(value, "ndim", None) == 0:
        return format_value(value.item())
    if isinstance(value, Mapping):
        return ", ".join(
            f"{key}={format_value(item)}" for key, item in value.items())
    if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        return ", ".join(format_value(item) for item in value)
    return str(value)


def print_settings(ostream: OutputStream | None, title: str, rows) -> None:
    """Print a titled ``Label : value`` block; ``rows`` holds (label, value) pairs.

    Rows whose value is ``None`` are skipped so optional settings stay quiet.
    """
    if ostream is None:
        return
    ostream.print_blank()
    ostream.print_header(title)
    ostream.print_header("=" * (len(title) + 2))
    ostream.print_blank()
    for label, value in rows:
        if value is None:
            continue
        # Wrap long values (paths, lists) under the value column so every
        # line keeps the block's centered alignment.
        indent = " " * (LABEL_WIDTH + 2)
        lines = textwrap.wrap(f"{label:<{LABEL_WIDTH}s}: {format_value(value)}",
                              SETTINGS_WIDTH,
                              subsequent_indent=indent,
                              break_on_hyphens=False) or [label]
        for line in lines:
            ostream.print_header(line.ljust(SETTINGS_WIDTH))
    ostream.print_blank()
    ostream.flush()


# ==========================================================================
# files
# ==========================================================================


RESULT_GROUP = "autocdam/result"


def h5_filename(filename) -> Path:
    """
    Gets the HDF5 file of a run (``.h5`` added when missing).

    :param filename:
        The file name, with or without the ``.h5`` suffix.

    :return:
        The path of the HDF5 file.
    """

    path = Path(filename)
    return path if path.suffix == '.h5' else path.with_name(path.name + '.h5')


def _check_h5_keys(value, where='value') -> None:
    """
    Checks that dictionary keys can be HDF5 names (no '/' and not '.').

    :param value:
        The value to be written.
    :param where:
        The location used in the error message.
    """

    if isinstance(value, Mapping):
        for key, item in value.items():
            if isinstance(key, str) and ('/' in key or key in ('', '.')):
                raise ValueError(
                    f"key {key!r} in {where} cannot be stored in HDF5")
            _check_h5_keys(item, f"{where}[{key!r}]")
    elif isinstance(value, (list, tuple)):
        for item in value:
            _check_h5_keys(item, where)


def write_h5_groups(filename, groups: Mapping[str, Any]) -> Path:
    """
    Writes (or replaces) named groups of a run file through a temporary copy.

    Values are written with VeloxChem's results serializer and keep their
    dictionary order; other groups of the file are left unchanged. The file
    is replaced only after the complete new version is written.

    :param filename:
        The run file (``.h5`` added when missing).
    :param groups:
        The values by group path, such as ``{'autocdam/result': {...}}``.

    :return:
        The path of the HDF5 file.
    """

    import shutil

    import h5py

    from ..resultsio import _write_value_to_hdf5

    for name, value in groups.items():
        _check_h5_keys(value, name)
    path = h5_filename(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    if path.exists():
        shutil.copyfile(path, temporary)
    config = h5py.get_config()
    track_order = config.track_order
    config.track_order = True
    try:
        with h5py.File(temporary, 'a') as handle:
            for name, value in groups.items():
                if name in handle:
                    del handle[name]
                parent, _, key = name.rpartition('/')
                group = handle.require_group(parent) if parent else handle
                _write_value_to_hdf5(group, key, value, 'AutoCDAM value')
    finally:
        config.track_order = track_order
    temporary.replace(path)
    return path


def read_h5_group(filename, name: str, default: Any = None) -> Any:
    """
    Reads one named group of a run file.

    :param filename:
        The run file (``.h5`` added when missing).
    :param name:
        The group path, such as ``'autocdam/result'``.
    :param default:
        The value returned when the file or group does not exist.

    :return:
        The stored value, or ``default``.
    """

    import h5py

    from ..resultsio import _read_value_from_hdf5

    path = h5_filename(filename)
    if not path.exists():
        return default
    with h5py.File(path, 'r') as handle:
        if name not in handle:
            return default
        return _read_value_from_hdf5(handle[name], 'AutoCDAM value')


# ==========================================================================
# units
# ==========================================================================

# Exact by definition of the thermochemical calorie; the ratio of VeloxChem's
# hartree conversions differs in the last digits, which perturbs MD trajectories.
KCAL_TO_KJ = 4.184
ANGSTROM_TO_NM = 0.1


def site_parameters_openmm(*, sigma, epsilon, c4, distance, bond_k,
                           angle_k) -> dict[str, float]:
    """CDAM site parameters in OpenMM units.

    Input: sigma and distance in angstrom, epsilon in kcal/mol, C4 in
    kcal/mol*angstrom^4, bond k in kcal/mol/angstrom^2, angle k in kcal/mol/rad^2.
    Output: nm, kJ/mol, kJ/mol*nm^4, kJ/mol/nm^2 and kJ/mol/rad^2.
    """
    return {
        "sigma": float(sigma) * ANGSTROM_TO_NM,
        "epsilon": float(epsilon) * KCAL_TO_KJ,
        "c4": float(c4) * KCAL_TO_KJ * ANGSTROM_TO_NM**4,
        "distance": float(distance) * ANGSTROM_TO_NM,
        "bond_k": bond_k * KCAL_TO_KJ / ANGSTROM_TO_NM**2,
        "angle_k": angle_k * KCAL_TO_KJ,
    }


# ==========================================================================
# reference
# ==========================================================================

# Bounds as multiples of the UFF value (sigma in angstrom, epsilon in kcal/mol).
UFF_BOUND_FACTORS = {
    "metal_lj_sigma": (0.3, 1.1),
    "metal_lj_epsilon": (0.1, 10.0),
}
# Limits that no widening may cross.
HARD_LIMITS = {
    "metal_lj_sigma": (0.5, 8.0),
    "metal_lj_epsilon": (1e-4, 100.0),
}
REFERENCE_KEYWORDS = ("uff", "auto")


def uff_lj(element: str) -> tuple[float, float]:
    """UFF Lennard-Jones (sigma / angstrom, epsilon / kcal/mol) of one element."""
    from ..uffparameters import get_uff_parameters

    table = get_uff_parameters()
    assert_msg_critical(element in table,
                        f"no UFF parameters for element {element!r}")
    values = table[element]
    return 10.0 * values["sigma"], values["epsilon"] / KCAL_TO_KJ


def reference_bounds(element: str) -> dict[str, dict]:
    """Soft first bounds per parameter name, derived from UFF for ``element``."""
    references = dict(
        zip(("metal_lj_sigma", "metal_lj_epsilon"), uff_lj(element)))
    result = {}
    for name, value in references.items():
        low, high = UFF_BOUND_FACTORS[name]
        hard_low, hard_high = HARD_LIMITS[name]
        result[name] = {
            "bounds": (max(hard_low, low * value), min(hard_high,
                                                       high * value)),
            "hard_bounds": HARD_LIMITS[name],
            "reference": value,
            "source": f"UFF {element} x [{low:g}, {high:g}]",
        }
    return result


# ==========================================================================
# models
# ==========================================================================

CoordinateFrame = list[list[float]]
FitStatus = Literal["success", "warning", "failed"]


@dataclass
class SystemData:
    """
    Topology and basic force-field data used by a generated batch.

    The containers are deliberately permissive because topology may originate
    from an external ForceFieldGenerator.  A generator may create a safely
    copied, augmented instance when dummy or probe atoms are required.

    :param atoms:
        The atoms.
    :param residues:
        The residues.
    :param bonds:
        The bonds.
    :param atom_types:
        The atom types (same length as ``atoms`` when given).
    :param charges:
        The charges (same length as ``atoms`` when given).
    :param bonded_parameters:
        The bonded force-field parameters.
    :param nonbonded_parameters:
        The nonbonded force-field parameters.
    :param metadata:
        The metadata.
    """

    atoms: list[Any] = field(default_factory=list)
    residues: list[Any] = field(default_factory=list)
    bonds: list[Any] = field(default_factory=list)
    atom_types: list[Any] = field(default_factory=list)
    charges: list[float] = field(default_factory=list)
    bonded_parameters: dict[str, Any] = field(default_factory=dict)
    nonbonded_parameters: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Apply only low-risk consistency checks at the shared boundary."""

        if self.atoms and self.atom_types and len(self.atom_types) != len(
                self.atoms):
            raise ValueError("atom_types must have the same length as atoms")
        if self.atoms and self.charges and len(self.charges) != len(self.atoms):
            raise ValueError("charges must have the same length as atoms")


@dataclass
class OpenMMSystems:
    """
    Candidate-owned OpenMM systems and starting positions (in nm).

    Topologies can be shared read-only. Systems and positions must be copied
    before updating a proposal; live Context objects are never shared.

    :param vacuum_system:
        The vacuum OpenMM system.
    :param vacuum_topology:
        The vacuum topology.
    :param vacuum_positions:
        The vacuum positions.
    :param solvated_system:
        The solvated OpenMM system.
    :param solvated_topology:
        The solvated topology.
    :param solvated_positions:
        The solvated positions.
    :param model_metadata:
        The model metadata.
    """

    vacuum_system: Any
    vacuum_topology: Any
    vacuum_positions: Any
    solvated_system: Any
    solvated_topology: Any
    solvated_positions: Any
    model_metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class GeneratedBatch:
    """
    One shared topology plus coordinate frames for multiple configurations.

    Every frame must use identical atom ordering and the same topology.  A
    frame is represented as ``n_atoms`` rows of ``[x, y, z]`` coordinates.

    :param system:
        The shared topology.
    :param configurations:
        The coordinate frames.
    :param metadata:
        The metadata.
    :param parameters:
        The parameter set the batch was generated for.
    :param openmm_systems:
        The candidate-owned OpenMM systems.
    """

    system: SystemData
    configurations: list[CoordinateFrame]
    metadata: dict[str, Any] = field(default_factory=dict)
    parameters: "ParameterSet | None" = None
    openmm_systems: OpenMMSystems | None = None

    def __post_init__(self) -> None:
        """Check that every frame is non-empty, finite, 3-D and of one shape; raise ValueError."""
        if not self.configurations:
            raise ValueError(
                "GeneratedBatch requires at least one configuration")

        expected_atoms = len(self.system.atoms)
        first_shape: tuple[int, int] | None = None
        for frame_index, frame in enumerate(self.configurations):
            if expected_atoms and len(frame) != expected_atoms:
                raise ValueError(
                    f"configuration {frame_index} does not have the same frame shape: "
                    f"{len(frame)} atoms; expected {expected_atoms}")
            for atom_index, coordinates in enumerate(frame):
                if len(coordinates) != 3:
                    raise ValueError(
                        f"configuration {frame_index}, atom {atom_index} "
                        "must have three coordinates")
                if not all(
                        math.isfinite(float(value)) for value in coordinates):
                    raise ValueError("configuration coordinates must be finite")
            shape = (len(frame), 3)
            if first_shape is None:
                first_shape = shape
            elif shape != first_shape:
                raise ValueError(
                    "all configurations must have the same frame shape")

    @property
    def n_configurations(self) -> int:
        """
        Gets the number of coordinate frames in the batch.

        :return:
            The number of frames.
        """

        return len(self.configurations)


@dataclass(frozen=True)
class ParameterSpec:
    """
    Definition of one search parameter and its optional optimization bound.

    ``initial`` is the starting value and ``bounds`` the ``(lower, upper)``
    search interval, both in the parameter's own ``unit`` (CDAM parameters use
    angstrom and kcal/mol). A spec with ``optimizable=False`` is held fixed at
    ``initial``. ``expandable`` bounds are soft and may be widened up to
    ``hard_bounds``.

    :param name:
        The parameter name.
    :param initial:
        The starting value.
    :param bounds:
        The ``(lower, upper)`` search interval.
    :param unit:
        The unit.
    :param optimizable:
        Whether the parameter is adjusted.
    :param metadata:
        The metadata.
    :param expandable:
        Whether the bounds are soft.
    :param hard_bounds:
        The limits soft bounds are never widened past.
    """

    name: str
    initial: float = 0.0
    bounds: tuple[float, float] | None = None
    unit: str | None = None
    optimizable: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)
    # Soft bounds may be widened by the Bayesian search, never past hard_bounds.
    expandable: bool = False
    hard_bounds: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        """Validate name, finite initial value and bounds nesting (critical error)."""
        assert_msg_critical(bool(self.name), "parameter name cannot be empty")
        assert_msg_critical(math.isfinite(self.initial),
                            "initial parameter value must be finite")
        if self.bounds is not None:
            lower, upper = self.bounds
            assert_msg_critical(
                math.isfinite(lower) and math.isfinite(upper),
                "parameter bounds must be finite")
            assert_msg_critical(lower < upper,
                                "parameter bounds must satisfy lower < upper")
            assert_msg_critical(lower <= self.initial <= upper,
                                "initial parameter value must be inside bounds")
        if self.hard_bounds is not None:
            assert_msg_critical(self.bounds is not None,
                                "hard_bounds need bounds")
            assert_msg_critical(
                self.hard_bounds[0] <= self.bounds[0] < self.bounds[1] <=
                self.hard_bounds[1], "bounds must lie inside hard_bounds")

    @classmethod
    def from_mapping(cls, name: str, values: Mapping[str,
                                                     Any]) -> "ParameterSpec":
        """
        Creates a spec from the compact user-facing configuration mapping.

        :param name:
            The parameter name.
        :param values:
            The mapping with ``ParameterSpec`` keys (``value`` is an alias of
            ``initial``).

        :return:
            The parameter specification.
        """

        raw_bounds = values.get("bounds")
        bounds = tuple(
            float(value)
            for value in raw_bounds) if raw_bounds is not None else None
        initial = float(
            values.get("initial",
                       values.get("value", bounds[0] if bounds else 0.0)))
        return cls(
            name=name,
            initial=initial,
            bounds=bounds,  # type: ignore[arg-type]
            unit=values.get("unit"),
            optimizable=bool(values.get("optimizable", True)),
            metadata=dict(values.get("metadata", {})),
            expandable=bool(values.get("expandable", False)),
            hard_bounds=(tuple(float(v) for v in values["hard_bounds"])
                         if values.get("hard_bounds") is not None else None),
        )


BOUND_MARGIN = 0.02  # "at a bound" = within this fraction of the width


def widen_soft_bounds(specs: Sequence[ParameterSpec],
                      values: Mapping[str, float], factor: float,
                      step: int) -> tuple[list[ParameterSpec], list[dict]]:
    """
    Widens the soft bounds that the given parameter values have reached.

    A value within ``BOUND_MARGIN`` of the width from a side of an expandable
    spec moves that side out by ``factor`` times the width, never past
    ``hard_bounds``. Fixed and non-expandable specs are returned unchanged.

    :param specs:
        The parameter specifications.
    :param values:
        The parameter values (the current best point) by name.
    :param factor:
        The widening as a fraction of the current width.
    :param step:
        The search step recorded with each change.

    :return:
        The specs with widened bounds, and one record per widened side.
    """

    widened, changes = [], []
    for spec in specs:
        if not (spec.optimizable and spec.expandable):
            widened.append(spec)
            continue
        lower, upper = spec.bounds
        hard_lower, hard_upper = spec.hard_bounds or (-math.inf, math.inf)
        width, value = upper - lower, values[spec.name]
        margin = BOUND_MARGIN * width
        if value - lower <= margin and lower > hard_lower:
            side, new = "lower", (max(hard_lower,
                                      lower - factor * width), upper)
        elif upper - value <= margin and upper < hard_upper:
            side, new = "upper", (lower,
                                  min(hard_upper, upper + factor * width))
        else:
            widened.append(spec)
            continue
        widened.append(replace(spec, bounds=new))
        changes.append({
            "step": step,
            "parameter": spec.name,
            "side": side,
            "value": value,
            "old": [lower, upper],
            "new": list(new)
        })
    return widened, changes


@dataclass
class ParameterSet:
    """
    Current parameter values passed to one calculation.

    :param values:
        The parameter values by name.
    """

    values: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Coerce keys to str and values to float."""
        self.values = {
            str(name): float(value) for name, value in self.values.items()
        }

    def __getitem__(self, name: str) -> float:
        """Value of the named parameter; raises KeyError if absent."""
        return self.values[name]

    def __iter__(self):
        """Iterate over parameter names."""
        return iter(self.values)

    def __len__(self) -> int:
        """Number of parameters."""
        return len(self.values)

    def items(self):
        """
        Gets the (name, value) pairs, as for a dict.

        :return:
            The (name, value) pairs.
        """

        return self.values.items()

    def to_dict(self) -> dict[str, float]:
        """
        Gets a plain copy of the values.

        :return:
            The ``{name: value}`` dictionary.
        """

        return dict(self.values)


TARGET_ROLES = ("required", "weighted", "monitor")
LEGACY_ROLES = {"objective": "weighted", "analysis": "monitor"}
# Target keys that select where and how a property is calculated. They are kept
# in ``metadata`` so FitTarget stays a small, calculator-independent contract.
TARGET_SELECTION_KEYS = ("system", "property", "atoms", "selection",
                         "reference_structure", "method", "evaluate",
                         "cutoff_angstrom", "shell_size")


@dataclass
class FitTarget:
    """
    Reference data and objective settings for one named target.

    ``reference`` may be a scalar or a sequence.  Names such as HFE, IOD, CN,
    and QM_ENERGY are target labels, not architectural special cases.

    ``role`` is ``required`` (a hard tolerance constraint that also enters the
    loss), ``weighted`` (loss only) or ``monitor`` (reported, never optimized).
    The legacy names ``objective`` and ``analysis`` mean weighted and monitor.
    ``scale`` normalizes the error; it defaults to the tolerance.

    :param name:
        The target name.
    :param reference:
        The reference value (scalar or sequence).
    :param tolerance:
        The tolerance (positive).
    :param weight:
        The weight in the loss (nonnegative).
    :param metric:
        The error metric.
    :param unit:
        The unit of the reference.
    :param role:
        The role (required, weighted or monitor).
    :param metadata:
        The selection keys (system, property, atoms, ...) and other metadata.
    :param scale:
        The normalization scale of the error.
    """

    name: str
    reference: Any
    tolerance: float
    weight: float = 1.0
    metric: str = "rmse"
    unit: str | None = None
    role: str = "objective"
    metadata: dict[str, Any] = field(default_factory=dict)
    scale: float | None = None

    def __post_init__(self) -> None:
        """Validate name, positive tolerance/scale, non-negative weight and role (critical error)."""
        assert_msg_critical(bool(self.name), "target name cannot be empty")
        assert_msg_critical(
            math.isfinite(self.tolerance) and self.tolerance > 0,
            "target tolerance must be positive")
        assert_msg_critical(
            math.isfinite(self.weight) and self.weight >= 0,
            "target weight cannot be negative")
        assert_msg_critical(
            self.role in TARGET_ROLES or self.role in LEGACY_ROLES,
            f"target role must be one of {TARGET_ROLES}, got {self.role!r}")
        assert_msg_critical(
            self.scale is None or
            (math.isfinite(self.scale) and self.scale > 0),
            "target scale must be positive")

    @property
    def kind(self) -> str:
        """
        Gets the role with legacy names mapped to required/weighted/monitor.

        :return:
            The kind of the target.
        """

        return LEGACY_ROLES.get(self.role, self.role)

    @property
    def is_objective(self) -> bool:
        """
        Checks whether this target contributes to the weighted loss.

        :return:
            True unless the target is a monitor.
        """

        return self.kind != "monitor"

    @property
    def is_required(self) -> bool:
        """
        Checks whether this target is a hard tolerance constraint.

        :return:
            True for a required target.
        """

        return self.kind == "required"

    @property
    def is_monitor(self) -> bool:
        """
        Checks whether the target is only reported, never optimized.

        :return:
            True for a monitor target.
        """

        return self.kind == "monitor"

    @property
    def normalization_scale(self) -> float:
        """
        Gets the divisor for the error.

        :return:
            ``scale`` if set, otherwise ``tolerance``.
        """

        return self.scale if self.scale is not None else self.tolerance

    @property
    def system(self) -> str | None:
        """
        Gets the name of the system the target is calculated on.

        :return:
            The system name, or None if no system was selected.
        """

        return self.metadata.get("system")

    @property
    def property_name(self) -> str:
        """
        Gets the calculated property.

        :return:
            The property name; defaults to the target name (IOD, CN, ...).
        """

        return self.metadata.get("property", self.name)

    @classmethod
    def from_mapping(cls, name: str, values: Mapping[str, Any]) -> "FitTarget":
        """
        Creates a target from the intended user-facing dictionary syntax.

        :param name:
            The target name.
        :param values:
            The target settings: ``reference`` (or ``value``), ``tolerance``,
            ``weight``, ``metric``, ``unit``, ``role``, ``metadata``, ``scale`` and
            the selection keys.

        :return:
            The fit target.
        """

        if "reference" in values:
            reference = values["reference"]
        elif "value" in values:
            reference = values["value"]
        else:
            assert_msg_critical(
                False, f"target {name!r} needs a 'reference' or 'value'")
        known = {
            "reference", "value", "tolerance", "weight", "metric", "unit",
            "role", "metadata", "scale", *TARGET_SELECTION_KEYS
        }
        unknown = set(values) - known
        assert_msg_critical(
            not unknown,
            f"target {name!r} has unsupported keys {sorted(unknown)}")
        default_metric = "rmse" if _is_sequence(reference) else "absolute"
        metadata = dict(values.get("metadata", {}))
        metadata.update(
            {key: values[key] for key in TARGET_SELECTION_KEYS if key in values})
        return cls(
            name=name,
            reference=reference,
            tolerance=float(values.get("tolerance", 1.0)),
            weight=float(values.get("weight", 1.0)),
            metric=str(values.get("metric", default_metric)),
            unit=values.get("unit"),
            role=str(values.get("role", "objective")),
            metadata=metadata,
            scale=float(values["scale"])
            if values.get("scale") is not None else None,
        )


@dataclass
class CalculationResult:
    """
    Values produced by a calculator, independent of how they were obtained.

    ``values`` maps target names to calculated scalars or arrays; ``raw_data``
    holds calculator-specific intermediates. ``failed`` with
    ``failure_reason`` marks a calculation the evaluator must penalize.

    :param values:
        The calculated values by target name.
    :param metadata:
        The per-target metadata.
    :param raw_data:
        The calculator-specific intermediates.
    :param warnings:
        The warnings.
    :param failed:
        Whether the calculation failed.
    :param failure_reason:
        Why the calculation failed.
    """

    values: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    raw_data: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    failed: bool = False
    failure_reason: str | None = None


@dataclass
class EvaluationResult:
    """
    Normalized objective and target-level diagnostics for one proposal.

    ``objective`` is the normalized weighted loss (or the failure penalty).
    ``metrics`` holds each target's raw error in the target's unit, ``errors``
    the error divided by its normalization scale, and ``target_results`` a
    per-target row with prediction, reference and tolerance status.
    ``constraint_violations`` maps each required target outside tolerance to
    raw error / tolerance; ``feasible`` is False when any exist.

    :param objective:
        The normalized weighted loss.
    :param metrics:
        The raw errors.
    :param errors:
        The normalized errors.
    :param target_results:
        The per-target rows.
    :param failed:
        Whether the evaluation failed.
    :param failure_reason:
        Why the evaluation failed.
    :param warnings:
        The warnings.
    :param feasible:
        Whether all required targets are within tolerance.
    :param constraint_violations:
        The violated required targets.
    """

    objective: float
    metrics: dict[str, float] = field(default_factory=dict)
    errors: dict[str, float] = field(default_factory=dict)
    target_results: dict[str, dict[str, Any]] = field(default_factory=dict)
    failed: bool = False
    failure_reason: str | None = None
    warnings: list[str] = field(default_factory=list)
    feasible: bool = True
    constraint_violations: dict[str, float] = field(default_factory=dict)


@dataclass
class FitResult:
    """
    Final optimizer output, including history and diagnostic information.

    ``parameters`` and ``metrics`` describe the selected best candidate;
    ``history`` holds one record per evaluated step and ``diagnostics`` the
    target results, parameter definitions, cost and stop reason read by
    ``summary``. ``checkpoint`` is the run file and ``checkpoint_group`` the
    group of the search checkpoint in it, if one was written.

    :param parameters:
        The best parameters.
    :param metrics:
        The metrics of the best candidate.
    :param objective:
        The objective of the best candidate.
    :param status:
        The status (success, warning or failed).
    :param converged:
        Whether the optimization converged.
    :param iterations:
        The number of iterations.
    :param history:
        The per-step records.
    :param warnings:
        The warnings.
    :param diagnostics:
        The diagnostics.
    :param checkpoint:
        The run file.
    :param checkpoint_group:
        The group of the search checkpoint.
    """

    parameters: ParameterSet
    metrics: dict[str, float] = field(default_factory=dict)
    objective: float = math.inf
    status: FitStatus = "failed"
    converged: bool = False
    iterations: int = 0
    history: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    checkpoint: str | None = None
    checkpoint_group: str | None = None

    def plot(self, **kwargs):
        """
        Delegates CDAM or metal-complex plotting to the standalone module.

        :param kwargs:
            The keyword options of ``plot_fit``.

        :return:
            The Matplotlib figure.
        """

        from .plotting import plot_fit
        return plot_fit(self, **kwargs)

    def to_dict(self) -> dict[str, Any]:
        """
        Converts the complete result to a JSON-ready dictionary.

        :return:
            The dictionary.
        """

        import json
        objective = self.objective if math.isfinite(self.objective) else None
        checkpoint = str(self.checkpoint) if self.checkpoint else None
        group = self.checkpoint_group if checkpoint else None
        data = {
            "parameters": self.parameters.to_dict(),
            "metrics": self.metrics,
            "objective": objective,
            "status": self.status,
            "converged": self.converged,
            "iterations": self.iterations,
            "history": self.history,
            "warnings": self.warnings,
            "diagnostics": self.diagnostics,
            "checkpoint": checkpoint,
            "checkpoint_group": group,
        }
        return json.loads(json.dumps(data, default=_json_default))

    def samples_table(self) -> list[dict[str, Any]]:
        """
        Creates one flat row per calculated candidate.

        See ``veloxchem.autocdam.reporting``.

        :return:
            The list of rows.
        """

        from .reporting import samples_table
        return samples_table(self)

    def save(self, filename) -> Path:
        """
        Writes the result to the ``autocdam/result`` group of a run file.

        :param filename:
            The run file (``.h5`` added when missing).

        :return:
            The path of the HDF5 file.
        """

        return write_h5_groups(filename, {RESULT_GROUP: self.to_dict()})

    @classmethod
    def load(cls, filename) -> "FitResult":
        """
        Reads a result saved with ``save``.

        :param filename:
            The run file (``.h5`` added when missing).

        :return:
            The fit result.
        """

        data = read_h5_group(filename, RESULT_GROUP)
        if data is None:
            raise FileNotFoundError(
                f"no AutoCDAM result in {h5_filename(filename)}")
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "FitResult":
        """
        Rebuilds a result from ``to_dict`` output.

        :param data:
            The dictionary of a saved result.

        :return:
            The fit result.
        """

        data = dict(data)
        objective = data.get("objective")
        return cls(parameters=ParameterSet(dict(data["parameters"])),
                   metrics=dict(data.get("metrics") or {}),
                   objective=math.inf if objective is None else objective,
                   status=data.get("status", "failed"),
                   converged=bool(data.get("converged", False)),
                   iterations=int(data.get("iterations", 0)),
                   history=list(data.get("history") or []),
                   warnings=list(data.get("warnings") or []),
                   diagnostics=dict(data.get("diagnostics") or {}),
                   checkpoint=data.get("checkpoint"),
                   checkpoint_group=data.get("checkpoint_group"))

    def __repr__(self) -> str:
        """
        Creates the formatted result shown in notebooks and by ``print``.

        :return:
            The formatted result.
        """

        d = self.diagnostics
        title = "AutoCDAM Fit Result"
        lines = [title, "=" * len(title)]
        cost = d.get("cost", {})
        rows = [
            ("Status", f"{self.status} (converged: {self.converged}, "
             f"stop: {d.get('stop_reason')})"),
            ("Objective", f"{self.objective:.6g} (feasible: "
             f"{d.get('feasible')})"),
            ("Candidates", cost.get("n_candidates",
                                    d.get("n_calculations"))),
            ("Physical validation", d.get("physical_validation",
                                          "not performed")),
        ]
        lines += [f"{label:<20s}: {value}" for label, value in rows
                  if value is not None]
        if d.get("study_id"):
            lines.insert(2, f"{'Study ID':<20s}: {d['study_id']}")
        units = {
            row["name"]: row.get("unit") or ""
            for row in d.get("parameter_definitions", [])
        }
        lines += ["", f"{'Parameter':<22s}{'Value':>14s}  Unit"]
        lines += [f"{name:<22s}{value:>14.6g}  {units.get(name, '')}"
                  for name, value in self.parameters.items()]
        units = {
            row["name"]: row.get("unit") or ""
            for row in d.get("objectives", [])
        }
        results = d.get("target_results", {})
        if results or d.get("monitor"):
            lines += [
                "", f"{'Target':<16s}{'Calculated':>14s}{'Reference':>14s}"
                f"{'Error/Tol':>11s}  {'Met':<9s}Unit"
            ]
        for name, row in results.items():
            met = 'yes' if row.get('within_tolerance') else 'no'
            lines.append(f"{name:<16s}{_short(row['prediction']):>14s}"
                         f"{_short(row['reference']):>14s}"
                         f"{row['normalized_error']:>11.3g}  {met:<9s}"
                         f"{units.get(name, '')}")
        for name, row in d.get("monitor", {}).items():
            lines.append(f"{name:<16s}{_short(row['value']):>14s}"
                         f"{_short(row['reference']):>14s}{'':>11s}  "
                         f"{'monitor':<9s}{units.get(name, '')}")
        warnings = list(dict.fromkeys(self.warnings))
        if warnings:
            lines += [""] + [f"Warning: {warning}" for warning in warnings]
        return "\n".join(lines)

    def summary(self) -> str:
        """
        Creates a plain-text report.

        The report lists parameters, targets, constraints, cost, and status.

        :return:
            The report text.
        """

        d = self.diagnostics
        lines = [
            f"Status: {self.status} | optimization converged: {self.converged} | "
            f"stop: {d.get('stop_reason')}", "Physical validation: " +
            str(d.get("physical_validation", "not performed")), "Parameters:"
        ]
        definitions = {p["name"]: p for p in d.get("parameter_definitions", [])}
        for name, value in self.parameters.items():
            spec = definitions.get(name, {})
            kind = (f"adjustable {spec['bounds']}"
                    if spec.get("adjustable") else "fixed")
            if spec.get("expandable"):
                source = spec.get("bounds_source")
                kind += f", soft from {spec['initial_bounds']}" + (
                    f" ({source})" if source else "")
            lines.append(f"  {name} = {value:.6g} ({kind})")
        for change in d.get("bound_expansions", []):
            lines.append(
                f"  bound widened at step {change['step']}: {change['parameter']} "
                f"{change['side']} {change['old']} -> {change['new']}")
        lines.append(
            f"Objective: {self.objective:.6g} | feasible: {d.get('feasible')}")
        for name, row in d.get("target_results", {}).items():
            lines.append(
                f"  {name} [{row.get('kind', row.get('role'))}] = "
                f"{row['prediction']!s} vs {row['reference']!s}, "
                f"error/tolerance-scale {row['normalized_error']:.3g}, "
                f"within tolerance: {row.get('within_tolerance')}")
        for name, row in d.get("monitor", {}).items():
            lines.append(
                f"  {name} [monitor] = {row['value']} vs {row['reference']}")
        for name, system in d.get("systems", {}).items():
            if system.get("type") == "framework":
                lines.append(
                    f"System {name}: {system['baseline_source']}; target "
                    f"{system['target_metal']} ({system['representation']}); neighbouring "
                    f"parameters preserved: {system['non_target_preserved']}")
        cost = d.get("cost", {})
        if cost:
            lines.append(
                f"Cost: {cost.get('calculations', cost.get('n_candidates'))} "
                f"calculations, {cost.get('wall_time_s', cost.get('candidate_wall_time_s', 0)):.1f} s, "
                f"{cost.get('cache_hits', 0)} cached properties reused")
        lines.extend(
            f"Warning: {warning}" for warning in dict.fromkeys(self.warnings))
        return "\n".join(lines)


def _json_default(value):
    """``json.dumps`` fallback: NumPy to Python, non-finite floats to None, else ``str``."""
    import numpy as np
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return str(value)


def _is_sequence(value: Any) -> bool:
    """True for a non-string sequence or an array with at least one dimension."""
    return (isinstance(value, Sequence) and
            not isinstance(value,
                           (str, bytes))) or getattr(value, "ndim", 0) > 0


def all_finite(value: Any) -> bool:
    """True if ``value`` (a number or nested sequence) is numeric and finite everywhere."""
    import numpy as np
    try:
        return bool(np.isfinite(np.asarray(value, dtype=float)).all())
    except (TypeError, ValueError):
        return False


def write_atomic(path: Path, text: str) -> None:
    """Write ``text`` through a temporary file, so ``path`` is never left half written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def coerce_targets(
    targets: Mapping[str, FitTarget | Mapping[str, Any]] | Iterable[FitTarget],
) -> list[FitTarget]:
    """Normalize compact mappings and already-constructed targets."""

    if isinstance(targets, Mapping):
        return [
            t if isinstance(t, FitTarget) else FitTarget.from_mapping(name, t)
            for name, t in targets.items()
        ]
    return list(targets)


def parse_parameters(
    parameters: Mapping[str, Any],
    *,
    current: Mapping[str, float],
    allowed: Iterable[str],
    aliases: Mapping[str, str] | None = None,
    references: Mapping[str, Mapping[str, Any]] | None = None,
) -> list[ParameterSpec]:
    """Turn ``{"name": (lo, hi) | value | "uff" | mapping}`` into parameter specs.

    A two-value sequence is an adjustable hard bound, a scalar is a fixed value,
    and a mapping uses the ``ParameterSpec`` keys. ``"uff"`` (or ``"auto"``)
    takes soft, expandable bounds from ``references`` (see
    ``veloxchem.autocdam.models``). Bounded parameters start from the template's
    ``current`` value when it lies inside the bounds, otherwise from the
    midpoint. Unknown names are rejected so no proposal can be ignored.
    """

    aliases = dict(aliases or {})
    references = dict(references or {})
    allowed = set(allowed)

    def start(name, lower, upper):
        """The template's current value if inside the bounds, else the midpoint."""
        value = current.get(name)
        return float(
            value) if value is not None and lower <= value <= upper else 0.5 * (
                lower + upper)

    specs: list[ParameterSpec] = []
    for raw_name, value in parameters.items():
        name = aliases.get(raw_name, raw_name)
        assert_msg_critical(
            name in allowed,
            f"unsupported parameter {raw_name!r}; adjustable parameters are "
            f"{sorted(allowed)} (aliases {sorted(aliases)})")
        assert_msg_critical(
            all(spec.name != name for spec in specs),
            f"parameter {name!r} is given more than once")
        if isinstance(value, str):
            assert_msg_critical(
                value.lower() in REFERENCE_KEYWORDS,
                f"unknown bound keyword {value!r} for {raw_name!r}; use "
                f"(lower, upper), a fixed value, or one of {REFERENCE_KEYWORDS}")
            assert_msg_critical(
                name in references,
                f"no reference bounds for {raw_name!r}; give (lower, upper) "
                f"(references exist for {sorted(references)})")
            reference = references[name]
            lower, upper = reference["bounds"]
            specs.append(
                ParameterSpec(name,
                              initial=start(name, lower, upper),
                              bounds=(lower, upper),
                              expandable=True,
                              hard_bounds=reference["hard_bounds"],
                              metadata={"bounds_source": reference["source"]}))
        elif isinstance(value, ParameterSpec):
            specs.append(value)
        elif isinstance(value, Mapping):
            specs.append(ParameterSpec.from_mapping(name, value))
        elif _is_sequence(value):
            assert_msg_critical(
                len(value) == 2,
                f"bounds for {raw_name!r} must be (lower, upper)")
            lower, upper = float(value[0]), float(value[1])
            specs.append(
                ParameterSpec(name,
                              initial=start(name, lower, upper),
                              bounds=(lower, upper)))
        else:
            specs.append(
                ParameterSpec(name, initial=float(value), optimizable=False))
    assert_msg_critical(
        any(spec.optimizable for spec in specs),
        "at least one parameter needs (lower, upper) bounds to be fitted")
    return specs


def coerce_parameter_specs(
    parameters: Mapping[str, ParameterSpec | Mapping[str, Any]] |
    Iterable[ParameterSpec],
) -> list[ParameterSpec]:
    """Normalize compact parameter configuration into shared specs."""

    if isinstance(parameters, Mapping):
        return [
            p if isinstance(p, ParameterSpec) else ParameterSpec.from_mapping(
                name, p) for name, p in parameters.items()
        ]
    return list(parameters)


# ==========================================================================
# properties
# ==========================================================================


@dataclass(frozen=True)
class AqueousProperty:
    """One computable property: the job that produces it, its method label and canonical unit."""

    name: str
    job: str  # "md", "fep" or "qm": properties sharing a job share one simulation
    method: str  # description recorded in the result metadata
    unit: str  # canonical unit of the calculated value


def _register(*properties: AqueousProperty) -> dict[str, AqueousProperty]:
    """Index properties by name."""
    return {p.name: p for p in properties}


AQUEOUS_PROPERTIES = _register(
    AqueousProperty("HFE", "fep", "OpenMM FEP", "kJ/mol"),
    AqueousProperty("IOD", "md", "OpenMM MD trajectory analysis", "angstrom"),
    AqueousProperty("CN", "md", "OpenMM MD trajectory analysis",
                    "dimensionless"),
    # First hydration shell vs the ideal polyhedron.
    AqueousProperty("SHELL_DIST_RMSD", "md",
                    "first-shell metal-O distance RMSD (shell uniformity)",
                    "angstrom"),
    AqueousProperty("SHELL_ANGLE_RMSD", "md",
                    "first-shell O-M-O angle RMSD vs ideal geometry", "degree"),
    AqueousProperty("SHELL_SHAPE_RMSD", "md",
                    "first-shell Kabsch shape RMSD vs ideal geometry",
                    "angstrom"),
    # CDAM dummy star vs its equilibrium geometry.
    AqueousProperty("DUMMY_DIST_RMSD", "md",
                    "CDAM dummy metal-D bond RMSD vs equilibrium", "angstrom"),
    AqueousProperty("DUMMY_ANGLE_RMSD", "md",
                    "CDAM dummy D-M-D angle RMSD vs ideal", "degree"),
    AqueousProperty("DUMMY_SHAPE_RMSD", "md",
                    "CDAM dummy Kabsch shape RMSD vs ideal", "angstrom"),
    AqueousProperty("QM_ENERGY", "qm", "VeloxChem single-point SCF energy",
                    "hartree"),
)

GEOMETRY_RMSD_PROPERTIES = tuple(
    name for name in AQUEOUS_PROPERTIES if name.endswith("_RMSD"))

# Accepted spellings (lower case) and the factor from the canonical unit.
UNIT_FACTORS = {
    "kJ/mol": {
        "kj/mol": 1.0,
        "kcal/mol": 1.0 / KCAL_TO_KJ
    },
    "angstrom": {
        "angstrom": 1.0,
        "a": 1.0,
        "å": 1.0,
        "nm": 0.1
    },
    "dimensionless": {
        "dimensionless": 1.0,
        "1": 1.0
    },
    "degree": {
        "degree": 1.0,
        "deg": 1.0,
        "degrees": 1.0,
        "radian": math.pi / 180.0,
        "rad": math.pi / 180.0
    },
    "hartree": {
        "hartree": 1.0,
        "eh": 1.0,
        "au": 1.0,
        "a.u.": 1.0,
        "kj/mol": hartree_in_kjpermol(),
        "kcal/mol": hartree_in_kcalpermol()
    },
}


def aqueous_property(name: str) -> AqueousProperty:
    """Registered property for a target name; raises ``ValueError`` if unsupported."""
    if name not in AQUEOUS_PROPERTIES:
        raise ValueError(
            f"unsupported target {name!r}; Calculator supports HFE, IOD, CN, "
            "QM_ENERGY, and the geometry RMSD targets "
            "SHELL_DIST_RMSD/SHELL_ANGLE_RMSD/SHELL_SHAPE_RMSD and "
            "DUMMY_DIST_RMSD/DUMMY_ANGLE_RMSD/DUMMY_SHAPE_RMSD")
    return AQUEOUS_PROPERTIES[name]


def convert_units(values: dict, metadata: dict,
                  targets: Iterable[FitTarget]) -> None:
    """Convert calculated values (and their spread/uncertainty) to the target units in place."""
    for target in targets:
        if target.name not in AQUEOUS_PROPERTIES or target.name not in values:
            continue
        canonical = AQUEOUS_PROPERTIES[target.name].unit
        requested = target.unit or canonical
        factors = UNIT_FACTORS[canonical]
        key = requested.strip().lower()
        if key not in factors:
            raise ValueError(f"unsupported {target.name} unit: {requested}")
        factor = factors[key]
        values[target.name] *= factor
        metadata[target.name]["unit"] = requested
        metadata[target.name]["calculation_unit"] = canonical
        for name in ("spread", "uncertainty"):
            if name in metadata[target.name]:
                metadata[target.name][name] *= factor


def _short(value) -> str:
    """A value for the formatted result (vectors by their length)."""
    try:
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return f"[{len(value)} values]"
