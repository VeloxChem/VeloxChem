#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#

"""Schema-4 reusable local-factor domain objects and pure validation."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import math
from typing import Any

import numpy as np

from .imcorefamilies import (
    _FrozenRecord,
    _as_frozen_record,
    _dataclass_payload,
    _float_tuple,
    _nonnegative,
    _string_tuple,
    _unique_by_id,
    _validated_id,
    canonical_json,
    stable_msi_id,
)


FACTOR_BINDING_MODES = frozenset({"native", "shared", "disabled"})


def _int_tuple(value):
    return tuple(int(item) for item in (value or ()))


def _int_pairs(value, field_name):
    result = tuple(tuple(int(item) for item in pair) for pair in (value or ()))
    if any(len(pair) != 2 for pair in result):
        raise ValueError(f"{field_name} must contain pairs.")
    return result


def _string_pairs(value, field_name):
    result = tuple(tuple(str(item) for item in pair) for pair in (value or ()))
    if any(len(pair) != 2 for pair in result):
        raise ValueError(f"{field_name} must contain pairs.")
    return result


@dataclass(frozen=True)
class ProvenanceFingerprint:
    """Compatibility-critical molecular, coordinate, and QM provenance."""

    molecular_graph: str
    charge: int
    multiplicity: int
    electronic_state: str
    qm_method: str
    qm_backend: str
    basis: str
    dispersion: str
    coordinate_policy: str
    symmetry_policy: str
    units: str
    derivative_order: int
    extra: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        for field_name in (
            "molecular_graph",
            "electronic_state",
            "qm_method",
            "qm_backend",
            "basis",
            "coordinate_policy",
            "symmetry_policy",
            "units",
        ):
            value = str(getattr(self, field_name)).strip()
            if not value:
                raise ValueError(f"{field_name} must not be empty.")
            object.__setattr__(self, field_name, value)
        object.__setattr__(self, "dispersion", str(self.dispersion).strip())
        if int(self.multiplicity) < 1:
            raise ValueError("multiplicity must be positive.")
        if int(self.derivative_order) not in (0, 1, 2):
            raise ValueError("derivative_order must be 0, 1, or 2.")
        object.__setattr__(self, "charge", int(self.charge))
        object.__setattr__(self, "multiplicity", int(self.multiplicity))
        object.__setattr__(self, "derivative_order", int(self.derivative_order))
        object.__setattr__(
            self, "extra", _as_frozen_record(self.extra, "extra")
        )

    @property
    def digest(self):
        return hashlib.sha256(
            canonical_json(self.to_payload()).encode("utf-8")
        ).hexdigest()

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


@dataclass(frozen=True)
class FactorClassSpec:
    """Identity of one reusable local physical motion."""

    factor_class_id: str
    kind: str
    canonical_graph_fingerprint: str
    owned_atom_roles: tuple[tuple[str, int], ...]
    rotor_axes: tuple[tuple[int, int], ...]
    symmetry_orders: tuple[int, ...]
    canonical_coordinate_signature: tuple[str, ...]
    default_support_policy: str = "local_response_v1"
    default_relaxation_policy: str = "anchored_local_v1"
    provenance: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        object.__setattr__(
            self,
            "factor_class_id",
            _validated_id(self.factor_class_id, "factor_class_id"),
        )
        for field_name in (
            "kind",
            "canonical_graph_fingerprint",
            "default_support_policy",
            "default_relaxation_policy",
        ):
            value = str(getattr(self, field_name)).strip()
            if not value:
                raise ValueError(f"{field_name} must not be empty.")
            object.__setattr__(self, field_name, value)
        roles = tuple(
            (str(role), int(atom)) for role, atom in self.owned_atom_roles
        )
        if len({role for role, _ in roles}) != len(roles):
            raise ValueError("owned_atom_roles must have unique role names.")
        if len({atom for _, atom in roles}) != len(roles):
            raise ValueError("owned_atom_roles must map bijectively to atoms.")
        object.__setattr__(self, "owned_atom_roles", roles)
        axes = _int_pairs(self.rotor_axes, "rotor_axes")
        if any(left == right or left < 0 or right < 0 for left, right in axes):
            raise ValueError("rotor_axes must contain distinct non-negative atoms.")
        object.__setattr__(self, "rotor_axes", axes)
        orders = _int_tuple(self.symmetry_orders)
        if len(orders) != len(axes) or any(order < 1 for order in orders):
            raise ValueError(
                "symmetry_orders must contain one positive order per rotor axis."
            )
        object.__setattr__(self, "symmetry_orders", orders)
        object.__setattr__(
            self,
            "canonical_coordinate_signature",
            _string_tuple(self.canonical_coordinate_signature),
        )
        object.__setattr__(
            self, "provenance", _as_frozen_record(self.provenance, "provenance")
        )

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


@dataclass(frozen=True)
class FactorDefinitionSpec:
    """A stored native residual library and its canonical source anchor."""

    factor_definition_id: str
    factor_class_id: str
    source_core_family_id: str
    source_core_point_id: str
    state_datapoint_labels: tuple[str, ...]
    anchor_state_ids: tuple[str, ...]
    residual_derivative_order: int
    projector_policy: str
    support_policy: str
    symmetry_chart_policy: str
    coordinate_qm_fingerprint: str
    source_environment_descriptor: tuple[float, ...] = ()
    validated_source_domain_radius: float = 0.0
    provenance: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        for field_name in (
            "factor_definition_id",
            "factor_class_id",
            "source_core_family_id",
            "source_core_point_id",
        ):
            object.__setattr__(
                self,
                field_name,
                _validated_id(getattr(self, field_name), field_name),
            )
        labels = _string_tuple(self.state_datapoint_labels)
        anchors = _string_tuple(self.anchor_state_ids)
        if not labels:
            raise ValueError("state_datapoint_labels must not be empty.")
        if not anchors:
            raise ValueError("anchor_state_ids must not be empty.")
        object.__setattr__(self, "state_datapoint_labels", labels)
        object.__setattr__(self, "anchor_state_ids", anchors)
        order = int(self.residual_derivative_order)
        if order not in (1, 2):
            raise ValueError("residual_derivative_order must be 1 or 2.")
        object.__setattr__(self, "residual_derivative_order", order)
        for field_name in (
            "projector_policy",
            "support_policy",
            "symmetry_chart_policy",
            "coordinate_qm_fingerprint",
        ):
            value = str(getattr(self, field_name)).strip()
            if not value:
                raise ValueError(f"{field_name} must not be empty.")
            object.__setattr__(self, field_name, value)
        object.__setattr__(
            self,
            "source_environment_descriptor",
            _float_tuple(
                self.source_environment_descriptor,
                "source_environment_descriptor",
            ),
        )
        object.__setattr__(
            self,
            "validated_source_domain_radius",
            _nonnegative(
                self.validated_source_domain_radius,
                "validated_source_domain_radius",
            ),
        )
        object.__setattr__(
            self, "provenance", _as_frozen_record(self.provenance, "provenance")
        )

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


@dataclass(frozen=True)
class FactorBindingSpec:
    """How a target core family consumes a native or shared definition."""

    factor_binding_id: str
    target_core_family_id: str
    factor_class_id: str
    factor_definition_id: str
    mode: str
    atom_transport_map: tuple[tuple[int, int], ...] = ()
    row_transport_map: tuple[tuple[int, int], ...] = ()
    target_anchor_mapping: tuple[tuple[str, str], ...] = ()
    validation_record_ids: tuple[str, ...] = ()
    applicability_domain: Mapping[str, Any] = _FrozenRecord()
    precedence: int = 0
    provenance: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        for field_name in (
            "factor_binding_id",
            "target_core_family_id",
            "factor_class_id",
        ):
            object.__setattr__(
                self,
                field_name,
                _validated_id(getattr(self, field_name), field_name),
            )
        mode = str(self.mode).strip().lower()
        if mode not in FACTOR_BINDING_MODES:
            raise ValueError(
                f"mode must be one of {sorted(FACTOR_BINDING_MODES)}; correction "
                "surfaces are not implemented in schema-4 v1."
            )
        object.__setattr__(self, "mode", mode)
        definition_id = _validated_id(
            self.factor_definition_id,
            "factor_definition_id",
            allow_empty=(mode == "disabled"),
        )
        if mode != "disabled" and not definition_id:
            raise ValueError(f"A {mode} binding requires a factor definition.")
        object.__setattr__(self, "factor_definition_id", definition_id)
        object.__setattr__(
            self,
            "atom_transport_map",
            _int_pairs(self.atom_transport_map, "atom_transport_map"),
        )
        object.__setattr__(
            self,
            "row_transport_map",
            _int_pairs(self.row_transport_map, "row_transport_map"),
        )
        object.__setattr__(
            self,
            "target_anchor_mapping",
            _string_pairs(self.target_anchor_mapping, "target_anchor_mapping"),
        )
        object.__setattr__(
            self,
            "validation_record_ids",
            _string_tuple(self.validation_record_ids),
        )
        object.__setattr__(
            self,
            "applicability_domain",
            _as_frozen_record(self.applicability_domain, "applicability_domain"),
        )
        object.__setattr__(self, "precedence", int(self.precedence))
        object.__setattr__(
            self, "provenance", _as_frozen_record(self.provenance, "provenance")
        )

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


@dataclass(frozen=True)
class FactorValidationRecord:
    """Measured target-family evidence for a proposed shared definition."""

    validation_record_id: str
    factor_definition_id: str
    target_core_family_id: str
    probe_labels: tuple[str, ...]
    residual_energy_error: float
    projected_gradient_error: float
    omitted_gradient_error: float
    projected_hessian_error: float
    cross_hessian_error: float
    environment_distance: float
    threshold_policy: Mapping[str, Any]
    passed: bool
    qm_provenance_fingerprint: str
    created_at: str
    validator_version: str
    evidence: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        for field_name in (
            "validation_record_id",
            "factor_definition_id",
            "target_core_family_id",
        ):
            object.__setattr__(
                self,
                field_name,
                _validated_id(getattr(self, field_name), field_name),
            )
        labels = _string_tuple(self.probe_labels)
        if not labels:
            raise ValueError("probe_labels must not be empty.")
        object.__setattr__(self, "probe_labels", labels)
        for field_name in (
            "residual_energy_error",
            "projected_gradient_error",
            "omitted_gradient_error",
            "projected_hessian_error",
            "cross_hessian_error",
            "environment_distance",
        ):
            object.__setattr__(
                self,
                field_name,
                _nonnegative(getattr(self, field_name), field_name),
            )
        object.__setattr__(
            self,
            "threshold_policy",
            _as_frozen_record(self.threshold_policy, "threshold_policy"),
        )
        for field_name in (
            "qm_provenance_fingerprint",
            "created_at",
            "validator_version",
        ):
            value = str(getattr(self, field_name)).strip()
            if not value:
                raise ValueError(f"{field_name} must not be empty.")
            object.__setattr__(self, field_name, value)
        object.__setattr__(
            self, "evidence", _as_frozen_record(self.evidence, "evidence")
        )

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


def factor_fingerprint(factor_class, provenance):
    """Fingerprint the physical factor identity and compatibility provenance."""

    payload = factor_class.to_payload()
    payload.pop("factor_class_id", None)
    payload["coordinate_qm_provenance"] = provenance.to_payload()
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def provenance_compatibility(left, right):
    """Return ``(compatible, reasons)`` for strict schema-4 reuse fields."""

    fields = (
        "molecular_graph",
        "charge",
        "multiplicity",
        "electronic_state",
        "qm_method",
        "qm_backend",
        "basis",
        "dispersion",
        "coordinate_policy",
        "symmetry_policy",
        "units",
        "derivative_order",
    )
    reasons = tuple(
        f"{name}: {getattr(left, name)!r} != {getattr(right, name)!r}"
        for name in fields
        if getattr(left, name) != getattr(right, name)
    )
    return not reasons, reasons


def validate_transport_map(
        binding, factor_class=None, n_target_atoms=None, n_target_rows=None):
    """Validate bijective canonical-to-target atom and row transport."""

    if binding.mode == "disabled":
        return

    atom_sources = tuple(source for source, _ in binding.atom_transport_map)
    atom_targets = tuple(target for _, target in binding.atom_transport_map)
    row_sources = tuple(source for source, _ in binding.row_transport_map)
    row_targets = tuple(target for _, target in binding.row_transport_map)

    for name, values in (
        ("atom sources", atom_sources),
        ("atom targets", atom_targets),
        ("row sources", row_sources),
        ("row targets", row_targets),
    ):
        if len(set(values)) != len(values):
            raise ValueError(f"Transport {name} must be unique.")
        if any(value < 0 for value in values):
            raise ValueError(f"Transport {name} must be non-negative.")

    if factor_class is not None:
        required_atoms = {atom for _, atom in factor_class.owned_atom_roles}
        if required_atoms and set(atom_sources) != required_atoms:
            raise ValueError(
                "atom_transport_map sources do not match factor owned atoms."
            )

    if binding.mode == "shared":
        if not binding.atom_transport_map:
            raise ValueError("A shared binding requires an atom transport map.")
        if not binding.target_anchor_mapping:
            raise ValueError("A shared binding requires a target anchor mapping.")

    if n_target_atoms is not None and any(
            atom >= int(n_target_atoms) for atom in atom_targets):
        raise ValueError("atom_transport_map contains an out-of-range target atom.")
    if n_target_rows is not None and any(
            row >= int(n_target_rows) for row in row_targets):
        raise ValueError("row_transport_map contains an out-of-range target row.")


def resolve_factor_binding(bindings, target_core_family_id, factor_class_id):
    """Resolve one deterministic static family-level binding by precedence."""

    candidates = [
        binding
        for binding in bindings
        if binding.target_core_family_id == str(target_core_family_id)
        and binding.factor_class_id == str(factor_class_id)
    ]
    if not candidates:
        return None
    maximum = max(binding.precedence for binding in candidates)
    winners = [
        binding for binding in candidates if binding.precedence == maximum
    ]
    identities = {
        (binding.mode, binding.factor_definition_id) for binding in winners
    }
    if len(identities) != 1:
        raise ValueError(
            "Conflicting factor bindings have equal precedence for target "
            f"family {target_core_family_id!r}, class {factor_class_id!r}."
        )
    return min(winners, key=lambda binding: binding.factor_binding_id)


def validate_factor_domain(
        factor_classes,
        factor_definitions,
        factor_bindings,
        validation_records,
        core_family_ids=(),
        core_point_ids=()):
    """Validate entity references and evidence requirements for schema-4 v1."""

    class_by_id = _unique_by_id(
        factor_classes, "factor_class_id", "factor class"
    )
    definition_by_id = _unique_by_id(
        factor_definitions, "factor_definition_id", "factor definition"
    )
    binding_by_id = _unique_by_id(
        factor_bindings, "factor_binding_id", "factor binding"
    )
    record_by_id = _unique_by_id(
        validation_records, "validation_record_id", "validation record"
    )
    del binding_by_id
    family_ids = {str(item) for item in core_family_ids}
    point_ids = {str(item) for item in core_point_ids}

    for definition in factor_definitions:
        if definition.factor_class_id not in class_by_id:
            raise ValueError(
                f"Factor definition {definition.factor_definition_id!r} "
                f"references missing class {definition.factor_class_id!r}."
            )
        if family_ids and definition.source_core_family_id not in family_ids:
            raise ValueError(
                f"Factor definition {definition.factor_definition_id!r} "
                "references a missing source family."
            )
        if point_ids and definition.source_core_point_id not in point_ids:
            raise ValueError(
                f"Factor definition {definition.factor_definition_id!r} "
                "references a missing source point."
            )

    for record in validation_records:
        if record.factor_definition_id not in definition_by_id:
            raise ValueError(
                f"Validation record {record.validation_record_id!r} references "
                "a missing factor definition."
            )
        if family_ids and record.target_core_family_id not in family_ids:
            raise ValueError(
                f"Validation record {record.validation_record_id!r} references "
                "a missing target family."
            )

    seen_targets = set()
    for binding in factor_bindings:
        if binding.factor_class_id not in class_by_id:
            raise ValueError(
                f"Binding {binding.factor_binding_id!r} references a missing class."
            )
        if family_ids and binding.target_core_family_id not in family_ids:
            raise ValueError(
                f"Binding {binding.factor_binding_id!r} references a missing family."
            )
        definition = None
        if binding.mode != "disabled":
            definition = definition_by_id.get(binding.factor_definition_id)
            if definition is None:
                raise ValueError(
                    f"Binding {binding.factor_binding_id!r} references a missing "
                    "definition."
                )
            if definition.factor_class_id != binding.factor_class_id:
                raise ValueError(
                    f"Binding {binding.factor_binding_id!r} and definition use "
                    "different factor classes."
                )
        if binding.mode == "native" and (
                definition.source_core_family_id != binding.target_core_family_id):
            raise ValueError(
                f"Native binding {binding.factor_binding_id!r} must use a "
                "definition sourced from its target family."
            )
        if binding.mode == "shared":
            if definition.source_core_family_id == binding.target_core_family_id:
                raise ValueError(
                    f"Shared binding {binding.factor_binding_id!r} must use a "
                    "definition from a different family."
                )
            records = []
            for record_id in binding.validation_record_ids:
                record = record_by_id.get(record_id)
                if record is None:
                    raise ValueError(
                        f"Shared binding {binding.factor_binding_id!r} references "
                        f"missing validation record {record_id!r}."
                    )
                records.append(record)
            if not any(
                record.passed
                and record.factor_definition_id == binding.factor_definition_id
                and record.target_core_family_id == binding.target_core_family_id
                for record in records
            ):
                raise ValueError(
                    f"Shared binding {binding.factor_binding_id!r} has no passing "
                    "compatible validation evidence."
                )
        validate_transport_map(binding, class_by_id[binding.factor_class_id])

        target = (binding.target_core_family_id, binding.factor_class_id)
        if target not in seen_targets:
            resolve_factor_binding(
                factor_bindings, binding.target_core_family_id, binding.factor_class_id
            )
            seen_targets.add(target)


@dataclass(frozen=True)
class CanonicalResidualState:
    """Anchor-relative local residual jet used by native and shared bindings."""

    state_id: str
    relative_coordinates: tuple[float, ...]
    residual_energy: float
    residual_gradient: tuple[float, ...]
    residual_hessian: tuple[tuple[float, ...], ...]
    source_datapoint_label: str
    residual_cardinal: bool = True

    def __post_init__(self):
        object.__setattr__(self, "state_id", _validated_id(self.state_id, "state_id"))
        object.__setattr__(
            self,
            "relative_coordinates",
            _float_tuple(self.relative_coordinates, "relative_coordinates"),
        )
        energy = float(self.residual_energy)
        if not math.isfinite(energy):
            raise ValueError("residual_energy must be finite.")
        object.__setattr__(self, "residual_energy", energy)
        gradient = _float_tuple(self.residual_gradient, "residual_gradient")
        hessian = tuple(
            _float_tuple(row, "residual_hessian") for row in self.residual_hessian
        )
        if len(gradient) != len(self.relative_coordinates):
            raise ValueError("Residual gradient and coordinate dimensions differ.")
        if len(hessian) != len(gradient) or any(
                len(row) != len(gradient) for row in hessian):
            raise ValueError("Residual Hessian must be square in local coordinates.")
        object.__setattr__(self, "residual_gradient", gradient)
        object.__setattr__(self, "residual_hessian", hessian)
        if not str(self.source_datapoint_label).strip():
            raise ValueError("source_datapoint_label must not be empty.")

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


@dataclass(frozen=True)
class CanonicalResidualLibrary:
    """A source-anchor-relative residual library in canonical local rows."""

    factor_definition_id: str
    source_rows: tuple[int, ...]
    periodic_source_rows: tuple[int, ...]
    anchor_state_id: str
    states: tuple[CanonicalResidualState, ...]
    coordinate_qm_fingerprint: str

    def __post_init__(self):
        object.__setattr__(
            self,
            "factor_definition_id",
            _validated_id(self.factor_definition_id, "factor_definition_id"),
        )
        rows = _int_tuple(self.source_rows)
        if not rows or len(set(rows)) != len(rows) or any(row < 0 for row in rows):
            raise ValueError("source_rows must be unique non-negative rows.")
        object.__setattr__(self, "source_rows", rows)
        periodic = _int_tuple(self.periodic_source_rows)
        if not set(periodic).issubset(set(rows)):
            raise ValueError("periodic_source_rows must be a subset of source_rows.")
        object.__setattr__(self, "periodic_source_rows", periodic)
        object.__setattr__(
            self, "anchor_state_id", _validated_id(
                self.anchor_state_id, "anchor_state_id"
            )
        )
        states = tuple(self.states)
        state_ids = {state.state_id for state in states}
        if len(state_ids) != len(states):
            raise ValueError("Canonical residual state IDs must be unique.")
        if self.anchor_state_id not in state_ids:
            raise ValueError("anchor_state_id is missing from residual states.")
        if any(
            len(state.relative_coordinates) != len(rows) for state in states
        ):
            raise ValueError("Residual state coordinate dimensions are inconsistent.")
        object.__setattr__(self, "states", states)
        fingerprint = str(self.coordinate_qm_fingerprint).strip()
        if not fingerprint:
            raise ValueError("coordinate_qm_fingerprint must not be empty.")
        object.__setattr__(self, "coordinate_qm_fingerprint", fingerprint)
        validate_anchor_zero_residual(self.anchor_state)

    @property
    def anchor_state(self):
        return next(
            state for state in self.states if state.state_id == self.anchor_state_id
        )

    def to_payload(self):
        return {
            "factor_definition_id": self.factor_definition_id,
            "source_rows": list(self.source_rows),
            "periodic_source_rows": list(self.periodic_source_rows),
            "anchor_state_id": self.anchor_state_id,
            "states": [state.to_payload() for state in self.states],
            "coordinate_qm_fingerprint": self.coordinate_qm_fingerprint,
        }

    @classmethod
    def from_payload(cls, payload):
        payload = dict(payload)
        payload["states"] = tuple(
            CanonicalResidualState.from_payload(item)
            for item in payload.get("states", ())
        )
        return cls(**payload)


def build_canonical_residual_library(
        definition, datapoints_by_label, source_rows, periodic_source_rows=()):
    """Subtract one stored source anchor E/G/H jet from every factor state."""

    source_rows = _int_tuple(source_rows)
    if not source_rows:
        raise ValueError("A canonical residual library requires source rows.")
    state_datapoints = {}
    for label in definition.state_datapoint_labels:
        datapoint = datapoints_by_label.get(label)
        if datapoint is None:
            raise ValueError(
                f"Factor definition {definition.factor_definition_id!r} is "
                f"missing datapoint {label!r}."
            )
        state_datapoints[label] = datapoint

    state_id_to_label = definition.provenance.get("state_id_to_label", {})
    anchor_label = None
    for anchor_id in definition.anchor_state_ids:
        if anchor_id in state_datapoints:
            anchor_label = anchor_id
            break
        candidate = state_id_to_label.get(str(anchor_id))
        if candidate in state_datapoints:
            anchor_label = str(candidate)
            break
        try:
            candidate = definition.state_datapoint_labels[int(anchor_id)]
        except (ValueError, IndexError):
            continue
        if candidate in state_datapoints:
            anchor_label = candidate
            break
    if anchor_label is None:
        raise ValueError(
            f"Factor definition {definition.factor_definition_id!r} has no "
            "resolvable source anchor."
        )

    anchor = state_datapoints[anchor_label]
    anchor_coordinates = np.asarray(
        anchor.internal_coordinates_values, dtype=np.float64
    )
    anchor_gradient = np.asarray(anchor.internal_gradient, dtype=np.float64)
    anchor_hessian = np.asarray(anchor.internal_hessian, dtype=np.float64)
    row_array = np.asarray(source_rows, dtype=np.int64)
    periodic = set(int(row) for row in periodic_source_rows)
    states = []

    for label in definition.state_datapoint_labels:
        datapoint = state_datapoints[label]
        coordinates = np.asarray(
            datapoint.internal_coordinates_values, dtype=np.float64
        )
        relative = coordinates[row_array] - anchor_coordinates[row_array]
        for local_index, source_row in enumerate(source_rows):
            if source_row in periodic:
                relative[local_index] = (
                    relative[local_index] + np.pi
                ) % (2.0 * np.pi) - np.pi

        gradient = (
            np.asarray(datapoint.internal_gradient, dtype=np.float64)
            - anchor_gradient
        )[row_array]
        if int(definition.residual_derivative_order) >= 2:
            hessian = (
                np.asarray(datapoint.internal_hessian, dtype=np.float64)
                - anchor_hessian
            )[np.ix_(row_array, row_array)]
        else:
            hessian = np.zeros((len(source_rows), len(source_rows)))
        energy = float(datapoint.energy) - float(anchor.energy)
        state_id = stable_msi_id(
            "residual_state",
            {
                "factor_definition_id": definition.factor_definition_id,
                "datapoint_label": label,
            },
        )
        states.append(
            CanonicalResidualState(
                state_id=state_id,
                relative_coordinates=tuple(float(value) for value in relative),
                residual_energy=energy,
                residual_gradient=tuple(float(value) for value in gradient),
                residual_hessian=tuple(
                    tuple(float(value) for value in row) for row in hessian
                ),
                source_datapoint_label=label,
            )
        )

    anchor_state = next(
        state for state in states if state.source_datapoint_label == anchor_label
    )
    return CanonicalResidualLibrary(
        factor_definition_id=definition.factor_definition_id,
        source_rows=source_rows,
        periodic_source_rows=tuple(sorted(periodic)),
        anchor_state_id=anchor_state.state_id,
        states=tuple(states),
        coordinate_qm_fingerprint=definition.coordinate_qm_fingerprint,
    )


def select_high_information_probe(library, environment_distances=None):
    """Select the non-anchor state with the largest curvature/response score."""

    environment_distances = environment_distances or {}
    scored = []
    for state in library.states:
        if state.state_id == library.anchor_state_id:
            continue
        gradient = np.asarray(state.residual_gradient, dtype=np.float64)
        hessian = np.asarray(state.residual_hessian, dtype=np.float64)
        score = (
            abs(float(state.residual_energy))
            + float(np.linalg.norm(gradient))
            + float(np.linalg.norm(hessian, ord="fro"))
            + float(environment_distances.get(state.source_datapoint_label, 0.0))
        )
        scored.append((score, state.source_datapoint_label, state))
    if not scored:
        raise ValueError("A validation probe requires a non-anchor residual state.")
    return max(scored, key=lambda item: (item[0], item[1]))[2]


def calculate_transfer_validation_metrics(
        donor_state,
        target_state,
        omitted_gradient_error=0.0,
        cross_hessian_error=0.0,
        environment_distance=0.0):
    """Calculate all schema-4 transfer metrics for one target-family probe."""

    donor_gradient = np.asarray(donor_state.residual_gradient, dtype=np.float64)
    target_gradient = np.asarray(target_state.residual_gradient, dtype=np.float64)
    donor_hessian = np.asarray(donor_state.residual_hessian, dtype=np.float64)
    target_hessian = np.asarray(target_state.residual_hessian, dtype=np.float64)
    if donor_gradient.shape != target_gradient.shape:
        raise ValueError("Donor and target projected gradients have different shapes.")
    if donor_hessian.shape != target_hessian.shape:
        raise ValueError("Donor and target projected Hessians have different shapes.")
    return {
        "residual_energy_error": abs(
            float(donor_state.residual_energy)
            - float(target_state.residual_energy)
        ),
        "projected_gradient_error": float(
            np.sqrt(np.mean((donor_gradient - target_gradient)**2))
        ),
        "omitted_gradient_error": _nonnegative(
            omitted_gradient_error, "omitted_gradient_error"
        ),
        "projected_hessian_error": float(
            np.sqrt(np.mean((donor_hessian - target_hessian)**2))
        ),
        "cross_hessian_error": _nonnegative(
            cross_hessian_error, "cross_hessian_error"
        ),
        "environment_distance": _nonnegative(
            environment_distance, "environment_distance"
        ),
    }


def build_transfer_validation_record(
        validation_record_id,
        factor_definition_id,
        target_core_family_id,
        probe_labels,
        metrics,
        thresholds,
        qm_provenance_fingerprint,
        created_at,
        validator_version="transfer.v1"):
    """Persist measured metrics and pass only when every threshold is met."""

    required = (
        "residual_energy_error",
        "projected_gradient_error",
        "omitted_gradient_error",
        "projected_hessian_error",
        "cross_hessian_error",
        "environment_distance",
    )
    missing = [name for name in required if name not in metrics]
    missing_thresholds = [name for name in required if name not in thresholds]
    if missing or missing_thresholds:
        raise ValueError(
            "Transfer validation requires every energy, gradient, Hessian, and "
            "environment metric and threshold."
        )
    passed = all(
        float(metrics[name]) <= float(thresholds[name]) for name in required
    )
    return FactorValidationRecord(
        validation_record_id=validation_record_id,
        factor_definition_id=factor_definition_id,
        target_core_family_id=target_core_family_id,
        probe_labels=tuple(probe_labels),
        residual_energy_error=metrics["residual_energy_error"],
        projected_gradient_error=metrics["projected_gradient_error"],
        omitted_gradient_error=metrics["omitted_gradient_error"],
        projected_hessian_error=metrics["projected_hessian_error"],
        cross_hessian_error=metrics["cross_hessian_error"],
        environment_distance=metrics["environment_distance"],
        threshold_policy=dict(thresholds),
        passed=passed,
        qm_provenance_fingerprint=qm_provenance_fingerprint,
        created_at=created_at,
        validator_version=validator_version,
        evidence={"all_required_metrics_present": True},
    )


def choose_shared_or_native_binding(
        provisional_shared_binding, validation_record, native_fallback_binding):
    """Accept passing sharing evidence or require a complete native fallback."""

    shared = provisional_shared_binding
    if shared.mode != "shared":
        raise ValueError("The provisional binding must use shared mode.")
    compatible_record = (
        validation_record.factor_definition_id == shared.factor_definition_id
        and validation_record.target_core_family_id
        == shared.target_core_family_id
    )
    if compatible_record and validation_record.passed:
        return shared
    native = native_fallback_binding
    if native is None or native.mode != "native":
        raise ValueError(
            "Failed sharing validation requires a complete native target-family "
            "factor library."
        )
    if (
        native.target_core_family_id != shared.target_core_family_id
        or native.factor_class_id != shared.factor_class_id
    ):
        raise ValueError("Native fallback targets the wrong family or factor class.")
    return native


def validate_anchor_zero_residual(
        residual_state,
        energy_tolerance=1.0e-10,
        gradient_tolerance=1.0e-12,
        hessian_tolerance=1.0e-12):
    """Reject a transported residual anchor whose E/G/H jet is not zero."""

    maximum_gradient = max(
        (abs(value) for value in residual_state.residual_gradient), default=0.0
    )
    maximum_hessian = max(
        (
            abs(value)
            for row in residual_state.residual_hessian
            for value in row
        ),
        default=0.0,
    )
    failures = []
    if abs(residual_state.residual_energy) > float(energy_tolerance):
        failures.append("energy")
    if maximum_gradient > float(gradient_tolerance):
        failures.append("gradient")
    if maximum_hessian > float(hessian_tolerance):
        failures.append("hessian")
    if failures:
        raise ValueError(
            "Residual anchor has a non-zero " + "/".join(failures) + " jet."
        )
