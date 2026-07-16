#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#

"""Schema-4 core-family domain objects and promotion records."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
import hashlib
import json
import math
import re
from typing import Any

import numpy as np


_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
_LOCALITY_CLASSES = frozenset(
    {"local", "boundary_response", "core_candidate", "ambiguous"}
)
_PROMOTION_DECISIONS = frozenset(
    {
        "no_action",
        "add_core_point_to_existing_family",
        "create_core_family",
        "add_factor_validation_probe",
        "build_native_factor_library",
    }
)


DEFAULT_LOCALITY_THRESHOLD_POLICY = {
    "policy_id": "schema4_locality_v1",
    "local_omitted_displacement": 0.15,
    "local_omitted_gradient_fraction": 0.10,
    "local_cross_hessian_fraction": 0.10,
    "local_environment_distance": 0.20,
    "core_omitted_displacement": 0.50,
    "core_omitted_gradient_fraction": 0.30,
    "core_cross_hessian_fraction": 0.25,
    "core_environment_distance": 0.50,
    "significant_response_fraction": 0.05,
    "epsilon": 1.0e-14,
}


@dataclass(frozen=True)
class _FrozenRecord(Mapping):
    """A small recursively immutable, deterministically ordered JSON mapping."""

    items_tuple: tuple[tuple[str, Any], ...] = ()

    def __iter__(self):
        return (key for key, _ in self.items_tuple)

    def __len__(self):
        return len(self.items_tuple)

    def __getitem__(self, key):
        for item_key, value in self.items_tuple:
            if item_key == key:
                return value
        raise KeyError(key)


def _freeze_json(value):
    if isinstance(value, _FrozenRecord):
        return value
    if isinstance(value, Mapping):
        return _FrozenRecord(
            tuple(
                (str(key), _freeze_json(item))
                for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
            )
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_json(item) for item in value)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Schema-4 JSON payloads require finite numbers.")
        return float(value)
    if hasattr(value, "item"):
        return _freeze_json(value.item())
    raise TypeError(
        "Schema-4 JSON payloads support mappings, sequences, strings, booleans, "
        f"finite numbers, and null; got {type(value).__name__}."
    )


def _thaw_json(value):
    if isinstance(value, _FrozenRecord):
        return {key: _thaw_json(item) for key, item in value.items_tuple}
    if isinstance(value, tuple):
        return [_thaw_json(item) for item in value]
    return value


def _as_frozen_record(value, field_name):
    if value is None:
        return _FrozenRecord()
    frozen = _freeze_json(value)
    if not isinstance(frozen, _FrozenRecord):
        raise TypeError(f"{field_name} must be a mapping.")
    return frozen


def _json_safe(value):
    if isinstance(value, _FrozenRecord):
        return _thaw_json(value)
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    return value


def _dataclass_payload(instance):
    return {
        field.name: _json_safe(getattr(instance, field.name))
        for field in fields(instance)
    }


def canonical_json(payload):
    """Return a stable compact JSON encoding used for schema-4 fingerprints."""

    return json.dumps(
        _json_safe(_freeze_json(payload)),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )


def stable_msi_id(namespace, payload, digest_length=20):
    """Construct a readable deterministic ID from a canonical JSON payload."""

    namespace = _validated_id(namespace, "namespace")
    digest = hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()
    return f"{namespace}:{digest[:int(digest_length)]}"


def _validated_id(value, field_name, allow_empty=False):
    value = str(value).strip()
    if not value and allow_empty:
        return value
    if not value or _ID_PATTERN.fullmatch(value) is None:
        raise ValueError(
            f"{field_name} must be a non-empty stable ID containing only "
            "letters, digits, '.', '_', ':', or '-'."
        )
    return value


def _string_tuple(value):
    return tuple(str(item) for item in (value or ()))


def _float_tuple(value, field_name):
    result = tuple(float(item) for item in (value or ()))
    if not all(math.isfinite(item) for item in result):
        raise ValueError(f"{field_name} must contain finite values.")
    return result


def _nonnegative(value, field_name):
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{field_name} must be finite and non-negative.")
    return value


@dataclass(frozen=True)
class EnvironmentDescriptorSpec:
    """Versioned definition of a differentiable environment descriptor."""

    descriptor_spec_id: str
    version: str
    feature_definitions: tuple[Any, ...]
    feature_scales: tuple[float, ...]
    distance_coefficient: float = 1.0
    invariance_policy: str = "rigid_motion_invariant"
    descriptor_kind: str = "contact_torsion_v1"
    provenance: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        object.__setattr__(
            self,
            "descriptor_spec_id",
            _validated_id(self.descriptor_spec_id, "descriptor_spec_id"),
        )
        object.__setattr__(self, "version", str(self.version).strip())
        if not self.version:
            raise ValueError("version must not be empty.")
        features = tuple(_freeze_json(item) for item in self.feature_definitions)
        if not all(isinstance(item, _FrozenRecord) for item in features):
            raise TypeError("Each feature definition must be a mapping.")
        object.__setattr__(self, "feature_definitions", features)
        scales = _float_tuple(self.feature_scales, "feature_scales")
        if len(scales) != len(features):
            raise ValueError(
                "feature_scales must have one entry per feature definition."
            )
        if any(scale <= 0.0 for scale in scales):
            raise ValueError("feature_scales must be strictly positive.")
        object.__setattr__(self, "feature_scales", scales)
        object.__setattr__(
            self,
            "distance_coefficient",
            _nonnegative(self.distance_coefficient, "distance_coefficient"),
        )
        if not str(self.invariance_policy).strip():
            raise ValueError("invariance_policy must not be empty.")
        if not str(self.descriptor_kind).strip():
            raise ValueError("descriptor_kind must not be empty.")
        object.__setattr__(
            self, "provenance", _as_frozen_record(self.provenance, "provenance")
        )

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


@dataclass(frozen=True)
class CoreFamilySpec:
    """A physical environment regime containing one or more outer points."""

    core_family_id: str
    root: str
    descriptor_spec_id: str
    prototype_descriptor: tuple[float, ...]
    descriptor_domain_radius: float
    environment_group_ids: tuple[str, ...] = ()
    factor_binding_ids: tuple[str, ...] = ()
    construction_policy_id: str = "schema4_default"
    provenance: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        for field_name in (
            "core_family_id",
            "root",
            "descriptor_spec_id",
            "construction_policy_id",
        ):
            object.__setattr__(
                self,
                field_name,
                _validated_id(getattr(self, field_name), field_name),
            )
        object.__setattr__(
            self,
            "prototype_descriptor",
            _float_tuple(self.prototype_descriptor, "prototype_descriptor"),
        )
        object.__setattr__(
            self,
            "descriptor_domain_radius",
            _nonnegative(
                self.descriptor_domain_radius, "descriptor_domain_radius"
            ),
        )
        object.__setattr__(
            self, "environment_group_ids", _string_tuple(self.environment_group_ids)
        )
        object.__setattr__(
            self, "factor_binding_ids", _string_tuple(self.factor_binding_ids)
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
class CorePointSpec:
    """An individual outer full-PES expansion point in a core family."""

    core_point_id: str
    core_family_id: str
    datapoint_label: str
    reference_descriptor: tuple[float, ...]
    confidence_radius: float
    full_pes_cardinal: bool = True
    factor_binding_ids: tuple[str, ...] = ()
    second_order_cardinal: bool = False
    provenance: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        for field_name in ("core_point_id", "core_family_id"):
            object.__setattr__(
                self,
                field_name,
                _validated_id(getattr(self, field_name), field_name),
            )
        label = str(self.datapoint_label).strip()
        if not label:
            raise ValueError("datapoint_label must not be empty.")
        object.__setattr__(self, "datapoint_label", label)
        object.__setattr__(
            self,
            "reference_descriptor",
            _float_tuple(self.reference_descriptor, "reference_descriptor"),
        )
        object.__setattr__(
            self,
            "confidence_radius",
            _nonnegative(self.confidence_radius, "confidence_radius"),
        )
        object.__setattr__(
            self, "factor_binding_ids", _string_tuple(self.factor_binding_ids)
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
class LocalityAuditRecord:
    """Persisted evidence used to classify a provisional local state."""

    locality_audit_id: str
    datapoint_label: str
    factor_class_id: str
    omitted_displacement: float
    omitted_gradient_fraction: float
    cross_hessian_fraction: float
    environment_distance: float
    contact_regime_changed: bool
    classification: str
    threshold_policy_id: str
    metrics: Mapping[str, Any] = _FrozenRecord()

    def __post_init__(self):
        for field_name in (
            "locality_audit_id",
            "factor_class_id",
            "threshold_policy_id",
        ):
            object.__setattr__(
                self,
                field_name,
                _validated_id(getattr(self, field_name), field_name),
            )
        if not str(self.datapoint_label).strip():
            raise ValueError("datapoint_label must not be empty.")
        for field_name in (
            "omitted_displacement",
            "omitted_gradient_fraction",
            "cross_hessian_fraction",
            "environment_distance",
        ):
            object.__setattr__(
                self,
                field_name,
                _nonnegative(getattr(self, field_name), field_name),
            )
        classification = str(self.classification).strip().lower()
        if classification not in _LOCALITY_CLASSES:
            raise ValueError(
                f"classification must be one of {sorted(_LOCALITY_CLASSES)}."
            )
        object.__setattr__(self, "classification", classification)
        object.__setattr__(
            self, "metrics", _as_frozen_record(self.metrics, "metrics")
        )

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


@dataclass(frozen=True)
class PromotionRecord:
    """A construction-time decision for an audited candidate state."""

    promotion_record_id: str
    datapoint_label: str
    locality_audit_id: str
    decision: str
    source_core_family_id: str = ""
    target_core_family_id: str = ""
    reason_codes: tuple[str, ...] = ()
    descriptor_distance: float = 0.0
    policy_id: str = "schema4_default"

    def __post_init__(self):
        for field_name in (
            "promotion_record_id",
            "locality_audit_id",
            "policy_id",
        ):
            object.__setattr__(
                self,
                field_name,
                _validated_id(getattr(self, field_name), field_name),
            )
        for field_name in ("source_core_family_id", "target_core_family_id"):
            object.__setattr__(
                self,
                field_name,
                _validated_id(
                    getattr(self, field_name), field_name, allow_empty=True
                ),
            )
        if not str(self.datapoint_label).strip():
            raise ValueError("datapoint_label must not be empty.")
        decision = str(self.decision).strip().lower()
        if decision not in _PROMOTION_DECISIONS:
            raise ValueError(
                f"decision must be one of {sorted(_PROMOTION_DECISIONS)}."
            )
        object.__setattr__(self, "decision", decision)
        object.__setattr__(self, "reason_codes", _string_tuple(self.reason_codes))
        object.__setattr__(
            self,
            "descriptor_distance",
            _nonnegative(self.descriptor_distance, "descriptor_distance"),
        )

    def to_payload(self):
        return _dataclass_payload(self)

    @classmethod
    def from_payload(cls, payload):
        return cls(**dict(payload))


def validate_core_domain(descriptor_specs, core_families, core_points):
    """Validate schema-4 core-domain IDs, references, and descriptor shapes."""

    descriptor_by_id = _unique_by_id(
        descriptor_specs, "descriptor_spec_id", "descriptor spec"
    )
    family_by_id = _unique_by_id(
        core_families, "core_family_id", "core family"
    )
    _unique_by_id(core_points, "core_point_id", "core point")

    for family in core_families:
        descriptor = descriptor_by_id.get(family.descriptor_spec_id)
        if descriptor is None:
            raise ValueError(
                f"Core family {family.core_family_id!r} references missing "
                f"descriptor {family.descriptor_spec_id!r}."
            )
        if len(family.prototype_descriptor) != len(descriptor.feature_scales):
            raise ValueError(
                f"Core family {family.core_family_id!r} prototype has the wrong "
                "descriptor dimension."
            )

    for point in core_points:
        family = family_by_id.get(point.core_family_id)
        if family is None:
            raise ValueError(
                f"Core point {point.core_point_id!r} references missing family "
                f"{point.core_family_id!r}."
            )
        if len(point.reference_descriptor) != len(family.prototype_descriptor):
            raise ValueError(
                f"Core point {point.core_point_id!r} has the wrong descriptor "
                "dimension."
            )


def _unique_by_id(items, attribute, entity_name):
    result = {}
    for item in items:
        item_id = str(getattr(item, attribute))
        if item_id in result:
            raise ValueError(f"Duplicate {entity_name} ID {item_id!r}.")
        result[item_id] = item
    return result


BOHR_IN_ANGSTROM = 0.529177210903


def evaluate_environment_descriptor(spec, cartesian_coordinates):
    """Evaluate ``contact_torsion_v1`` and its analytical Cartesian Jacobian."""

    coordinates = np.asarray(cartesian_coordinates, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("cartesian_coordinates must have shape (natoms, 3).")
    if spec.descriptor_kind == "legacy_none":
        return np.zeros(0, dtype=np.float64), np.zeros(
            (0, *coordinates.shape), dtype=np.float64
        )
    if spec.descriptor_kind != "contact_torsion_v1":
        raise ValueError(
            f"No analytical Jacobian implementation for descriptor kind "
            f"{spec.descriptor_kind!r}."
        )

    values = []
    jacobians = []
    for frozen_feature in spec.feature_definitions:
        feature = _thaw_json(frozen_feature)
        kind = str(feature.get("kind", ""))
        if kind == "smooth_contact":
            value, jacobian = _smooth_contact_feature(
                coordinates,
                int(feature["atom_i"]),
                int(feature["atom_j"]),
                float(feature["r0_bohr"]),
                int(feature.get("exponent", 6)),
            )
        elif kind == "smooth_contact_pool":
            pairs = tuple(tuple(int(atom) for atom in pair) for pair in feature["pairs"])
            reduction = str(feature.get("reduction", "sum"))
            value = 0.0
            jacobian = np.zeros_like(coordinates)
            for atom_i, atom_j in pairs:
                pair_value, pair_jacobian = _smooth_contact_feature(
                    coordinates,
                    atom_i,
                    atom_j,
                    float(feature["r0_bohr"]),
                    int(feature.get("exponent", 6)),
                )
                value += pair_value
                jacobian += pair_jacobian
            if reduction == "mean" and pairs:
                value /= len(pairs)
                jacobian /= len(pairs)
            elif reduction != "sum":
                raise ValueError(
                    "smooth_contact_pool reduction must be 'sum' or 'mean'."
                )
        elif kind in ("distance_mean", "distance_rms"):
            pairs = tuple(tuple(int(atom) for atom in pair) for pair in feature["pairs"])
            value, jacobian = _distance_pool_feature(
                coordinates, pairs, rms=(kind == "distance_rms")
            )
        elif kind in ("torsion_sin", "torsion_cos"):
            torsion, torsion_jacobian = _torsion_value_jacobian(
                coordinates,
                tuple(int(atom) for atom in feature["atoms"]),
            )
            if kind == "torsion_sin":
                value = float(np.sin(torsion))
                jacobian = np.cos(torsion) * torsion_jacobian
            else:
                value = float(np.cos(torsion))
                jacobian = -np.sin(torsion) * torsion_jacobian
        else:
            raise ValueError(
                f"Unsupported contact_torsion_v1 feature kind {kind!r}."
            )
        values.append(float(value))
        jacobians.append(np.asarray(jacobian, dtype=np.float64))

    return np.asarray(values), np.asarray(jacobians)


def descriptor_squared_distance(spec, value, reference, jacobian=None):
    """Return the scaled descriptor distance and optional Cartesian derivative."""

    value = np.asarray(value, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    scales = np.asarray(spec.feature_scales, dtype=np.float64)
    if value.shape != reference.shape or value.shape != scales.shape:
        raise ValueError("Descriptor value/reference/scale dimensions differ.")
    delta = value - reference
    coefficient = float(spec.distance_coefficient)
    distance_squared = coefficient * float(np.dot(scales * delta, scales * delta))
    if jacobian is None:
        return distance_squared
    jacobian = np.asarray(jacobian, dtype=np.float64)
    if jacobian.shape[0] != value.size:
        raise ValueError("Descriptor Jacobian has the wrong feature dimension.")
    factors = 2.0 * coefficient * scales * scales * delta
    gradient = np.tensordot(factors, jacobian, axes=1)
    return distance_squared, gradient


def build_contact_torsion_v1_spec(
        descriptor_spec_id,
        labels,
        connectivity,
        local_group_model,
        z_matrix,
        distance_coefficient=1.0,
        contact_r0_angstrom=4.0,
        contact_exponent=6):
    """Build descriptor features from detected groups and graph topology."""

    labels = tuple(str(label) for label in labels)
    connectivity = np.asarray(connectivity)
    if connectivity.shape != (len(labels), len(labels)):
        raise ValueError("connectivity must be square with one row per atom.")
    heavy_atoms = {
        atom for atom, label in enumerate(labels)
        if label.strip().upper() != "H"
    }
    ring_atoms = _cycle_atoms(connectivity) & heavy_atoms
    moving_kinds = {
        "amide_side_chain",
        "anchored_linker",
        "alkyl_chain",
        "side_chain",
    }
    selected_clusters = [
        cluster for cluster in local_group_model.clusters.values()
        if str(cluster.family_label) in moving_kinds
    ]
    moving_atoms = {
        int(atom)
        for cluster in selected_clusters
        for atom in (
            tuple(getattr(cluster, "owned_atoms", ()))
            + tuple(getattr(cluster, "active_atoms", ()))
        )
        if int(atom) in heavy_atoms
    }
    core_atoms = set(int(atom) for atom in local_group_model.core_atoms)
    ring_core_atoms = (ring_atoms & core_atoms) or ring_atoms
    moving_atoms.difference_update(ring_core_atoms)
    pairs = tuple(
        (moving_atom, ring_atom)
        for moving_atom in sorted(moving_atoms)
        for ring_atom in sorted(ring_core_atoms)
        if _shortest_graph_distance(
            connectivity, moving_atom, ring_atom, maximum=3
        ) > 3
    )
    if not pairs:
        raise ValueError(
            "contact_torsion_v1 detection found no non-bonded side-chain/ring "
            "heavy-atom pairs."
        )

    feature_definitions = [
        {
            "kind": "smooth_contact_pool",
            "pairs": pairs,
            "r0_bohr": float(contact_r0_angstrom) / BOHR_IN_ANGSTROM,
            "exponent": int(contact_exponent),
            "reduction": "sum",
        },
        {"kind": "distance_mean", "pairs": pairs},
    ]
    selected_rotor_ids = {
        str(rotor_id)
        for cluster in selected_clusters
        for rotor_id in cluster.rotor_ids
    }
    torsions = []
    for rotor_id in sorted(selected_rotor_ids):
        rotor = local_group_model.rotors.get(rotor_id)
        phase_coordinate = tuple(
            int(atom) for atom in getattr(rotor, "phase_coordinate", ())
        )
        if len(phase_coordinate) == 4:
            torsions.append(phase_coordinate)
    if not torsions:
        flat_z_matrix = []
        for section in ("bonds", "angles", "dihedrals", "impropers"):
            flat_z_matrix.extend(z_matrix.get(section, ()))
        selected_axes = {
            tuple(sorted(tuple(int(atom) for atom in rotor.axis)))
            for rotor_id, rotor in local_group_model.rotors.items()
            if str(rotor_id) in selected_rotor_ids
        }
        torsions.extend(
            tuple(int(atom) for atom in coordinate)
            for coordinate in flat_z_matrix
            if len(coordinate) == 4
            and tuple(sorted((int(coordinate[1]), int(coordinate[2]))))
            in selected_axes
        )
    torsions = tuple(dict.fromkeys(torsions))
    for torsion in torsions:
        feature_definitions.extend((
            {"kind": "torsion_sin", "atoms": torsion},
            {"kind": "torsion_cos", "atoms": torsion},
        ))

    feature_scales = tuple(
        1.0 / max(len(pairs), 1)
        if feature["kind"] == "smooth_contact_pool"
        else (1.0 / (4.0 / BOHR_IN_ANGSTROM)
              if feature["kind"] == "distance_mean" else 1.0)
        for feature in feature_definitions
    )
    return EnvironmentDescriptorSpec(
        descriptor_spec_id=descriptor_spec_id,
        version="contact_torsion_v1",
        feature_definitions=tuple(feature_definitions),
        feature_scales=feature_scales,
        distance_coefficient=distance_coefficient,
        invariance_policy="pair_distance_and_periodic_torsion",
        descriptor_kind="contact_torsion_v1",
        provenance={
            "moving_heavy_atoms": sorted(moving_atoms),
            "ring_core_heavy_atoms": sorted(ring_core_atoms),
            "detected_cluster_ids": sorted(
                str(cluster.cluster_id) for cluster in selected_clusters
            ),
            "contact_r0_angstrom": float(contact_r0_angstrom),
            "contact_exponent": int(contact_exponent),
        },
    )


def build_contact_torsion_v1_spec_from_factor_topology(
        descriptor_spec_id,
        labels,
        connectivity,
        family_bank,
        z_matrix,
        distance_coefficient=1.0,
        contact_r0_angstrom=4.0,
        contact_exponent=6):
    """Build the graph descriptor when adapting an already stored factor bank.

    The migration path derives moving atoms and linker torsions from persisted
    detector topology.  It intentionally does not contain molecule-specific
    atom numbers.
    """

    model = local_group_model_from_factor_topology(family_bank, z_matrix)
    return build_contact_torsion_v1_spec(
        descriptor_spec_id=descriptor_spec_id,
        labels=labels,
        connectivity=connectivity,
        local_group_model=model,
        z_matrix=z_matrix,
        distance_coefficient=distance_coefficient,
        contact_r0_angstrom=contact_r0_angstrom,
        contact_exponent=contact_exponent,
    )


def local_group_model_from_factor_topology(family_bank, z_matrix):
    """Rehydrate the local-group topology required by migration/probe tools."""

    from .imlocalgroups import LocalCluster, LocalGroupModel, LocalRotor

    factor_info = family_bank.get("factor_info", {})
    terms = family_bank.get("terms", {})
    flat_z_matrix = []
    for section in ("bonds", "angles", "dihedrals", "impropers"):
        flat_z_matrix.extend(z_matrix.get(section, ()))
    rotors = {}
    clusters = {}
    moving_kinds = ("anchored_linker", "alkyl_chain", "amide_side_chain")
    for term_id, term in terms.items():
        policy = str(term.get("relaxation_policy_id", ""))
        identity = f"{str(term_id)} {policy}".lower()
        family_label = next(
            (kind for kind in moving_kinds if kind in identity),
            str(term_id).split("_", 1)[0],
        )
        rotor_ids = tuple(str(value) for value in term.get("rotor_ids", ()))
        grouped_rows = tuple(term.get("grouped_signature_rows", ()))
        for rotor_index, rotor_id in enumerate(rotor_ids):
            signature_rows = tuple(
                int(row) for row in (
                    grouped_rows[rotor_index]
                    if rotor_index < len(grouped_rows) else ()
                )
            )
            torsion_rows = tuple(
                row for row in signature_rows
                if 0 <= row < len(flat_z_matrix)
                and len(flat_z_matrix[row]) == 4
            )
            phase_coordinate = (
                tuple(int(atom) for atom in flat_z_matrix[torsion_rows[0]])
                if torsion_rows else None
            )
            axis = (
                (phase_coordinate[1], phase_coordinate[2])
                if phase_coordinate is not None else (0, 1)
            )
            rotors[rotor_id] = LocalRotor(
                rotor_id=rotor_id,
                kind=family_label,
                axis=axis,
                symmetry_order=1,
                owned_atoms=tuple(int(atom) for atom in term.get(
                    "active_atoms", ())),
                unit_atom_sets=(),
                torsion_rows=torsion_rows,
                phase_coordinate=phase_coordinate,
                signature_rows=signature_rows,
            )
        clusters[str(term_id)] = LocalCluster(
            cluster_id=str(term_id),
            family_label=family_label,
            rotor_ids=rotor_ids,
            owned_atoms=tuple(int(atom) for atom in term.get(
                "active_atoms", ())),
            active_atoms=tuple(int(atom) for atom in term.get(
                "active_atoms", ())),
            active_rows=tuple(int(row) for row in term.get(
                "active_rows", ())),
        )
    return LocalGroupModel(
        version=int(factor_info.get("schema_version", 3)),
        enabled=bool(clusters),
        rotors=rotors,
        clusters=clusters,
        core_atoms=tuple(int(atom) for atom in factor_info.get(
            "core_atoms", ())),
        core_rows=tuple(int(row) for row in factor_info.get("core_rows", ())),
    )


def select_descriptor_medoids(descriptors, count, feature_scales=None):
    """Select deterministic centroid seed plus farthest-point coverage medoids."""

    descriptors = np.asarray(descriptors, dtype=np.float64)
    if descriptors.ndim != 2 or descriptors.shape[0] == 0:
        raise ValueError("descriptors must be a non-empty two-dimensional array.")
    count = min(max(int(count), 1), descriptors.shape[0])
    scales = np.ones(descriptors.shape[1], dtype=np.float64)
    if feature_scales is not None:
        scales = np.asarray(feature_scales, dtype=np.float64)
        if scales.shape != (descriptors.shape[1],):
            raise ValueError("feature_scales has the wrong descriptor dimension.")
    scaled = descriptors * scales
    centroid = np.mean(scaled, axis=0)
    first = int(np.argmin(np.sum((scaled - centroid)**2, axis=1)))
    selected = [first]
    nearest_squared = np.sum((scaled - scaled[first])**2, axis=1)
    while len(selected) < count:
        candidate = int(np.argmax(nearest_squared))
        selected.append(candidate)
        candidate_squared = np.sum(
            (scaled - scaled[candidate])**2, axis=1
        )
        nearest_squared = np.minimum(nearest_squared, candidate_squared)
    return tuple(selected)


def audit_local_factor_state(
        datapoint_label,
        factor_class_id,
        state_internal_coordinates,
        anchor_internal_coordinates,
        state_internal_gradient,
        anchor_internal_gradient,
        state_internal_hessian,
        projector_rows,
        environment_distance=0.0,
        contact_regime_changed=False,
        boundary_rows=(),
        coordinate_scales=None,
        periodic_rows=(),
        mapping_compatible=True,
        threshold_policy=None):
    """Measure and classify state response outside a proposed factor support.

    The three primary metrics implement the schema-4 :math:`eta_q`,
    :math:`eta_g`, and :math:`eta_H` definitions.  Classification is deliberately
    driven by a persisted, versioned policy instead of implicit constants.
    """

    policy = dict(DEFAULT_LOCALITY_THRESHOLD_POLICY)
    if threshold_policy is not None:
        policy.update(dict(threshold_policy))
    policy_id = _validated_id(policy.pop("policy_id"), "policy_id")
    epsilon = float(policy["epsilon"])

    state_q = np.asarray(state_internal_coordinates, dtype=np.float64)
    anchor_q = np.asarray(anchor_internal_coordinates, dtype=np.float64)
    state_g = np.asarray(state_internal_gradient, dtype=np.float64).reshape(-1)
    anchor_g = np.asarray(anchor_internal_gradient, dtype=np.float64).reshape(-1)
    state_h = np.asarray(state_internal_hessian, dtype=np.float64)
    if state_q.shape != anchor_q.shape or state_q.ndim != 1:
        raise ValueError("State and anchor internal coordinates must be 1-D peers.")
    if state_g.shape != state_q.shape or anchor_g.shape != state_q.shape:
        raise ValueError("Internal gradients must match the coordinate dimension.")
    if state_h.shape != (state_q.size, state_q.size):
        raise ValueError("The internal Hessian has the wrong shape.")

    support = np.zeros(state_q.size, dtype=bool)
    rows = np.asarray(tuple(int(row) for row in projector_rows), dtype=np.int64)
    if rows.size and (np.min(rows) < 0 or np.max(rows) >= state_q.size):
        raise ValueError("projector_rows contains an out-of-range coordinate.")
    support[rows] = True
    omitted = ~support

    delta_q = state_q - anchor_q
    periodic = np.asarray(tuple(int(row) for row in periodic_rows), dtype=np.int64)
    if periodic.size:
        if np.min(periodic) < 0 or np.max(periodic) >= state_q.size:
            raise ValueError("periodic_rows contains an out-of-range coordinate.")
        delta_q[periodic] = (
            (delta_q[periodic] + np.pi) % (2.0 * np.pi) - np.pi
        )
    scales = np.ones(state_q.size, dtype=np.float64)
    if coordinate_scales is not None:
        scales = np.asarray(coordinate_scales, dtype=np.float64)
        if scales.shape != state_q.shape or np.any(scales <= 0.0):
            raise ValueError(
                "coordinate_scales must be positive and match the coordinates."
            )
    omitted_displacement = float(np.linalg.norm(scales[omitted] * delta_q[omitted]))

    delta_g = state_g - anchor_g
    gradient_norm = float(np.linalg.norm(delta_g))
    omitted_gradient_norm = float(np.linalg.norm(delta_g[omitted]))
    omitted_gradient_fraction = omitted_gradient_norm / (gradient_norm + epsilon)

    cross_block = state_h[np.ix_(support, omitted)]
    hessian_norm = float(np.linalg.norm(state_h, ord="fro"))
    cross_hessian_fraction = float(
        np.linalg.norm(cross_block, ord="fro") / (hessian_norm + epsilon)
    )

    significant_rows = set()
    if omitted_gradient_norm > epsilon:
        significant_cutoff = (
            float(policy["significant_response_fraction"])
            * omitted_gradient_norm
        )
        significant_rows.update(
            int(row) for row in np.flatnonzero(
                omitted & (np.abs(delta_g) >= significant_cutoff)
            )
        )
    if omitted_displacement > epsilon:
        displacement_norm = float(np.linalg.norm(scales * delta_q))
        significant_cutoff = (
            float(policy["significant_response_fraction"])
            * displacement_norm
        )
        significant_rows.update(
            int(row) for row in np.flatnonzero(
                omitted & (np.abs(scales * delta_q) >= significant_cutoff)
            )
        )
    boundary_set = {int(row) for row in boundary_rows}
    boundary_localized = bool(significant_rows) and significant_rows <= boundary_set

    metrics = {
        "mapping_compatible": bool(mapping_compatible),
        "boundary_localized": bool(boundary_localized),
        "significant_omitted_rows": sorted(significant_rows),
        "projector_rows": sorted(int(row) for row in rows),
        "boundary_rows": sorted(boundary_set),
        "gradient_residual_norm": gradient_norm,
        "omitted_gradient_norm": omitted_gradient_norm,
        "hessian_norm": hessian_norm,
        "thresholds": policy,
    }
    classification = classify_locality_metrics(
        omitted_displacement=omitted_displacement,
        omitted_gradient_fraction=omitted_gradient_fraction,
        cross_hessian_fraction=cross_hessian_fraction,
        environment_distance=environment_distance,
        contact_regime_changed=contact_regime_changed,
        boundary_localized=boundary_localized,
        mapping_compatible=mapping_compatible,
        threshold_policy=policy,
    )
    audit_id = stable_msi_id(
        "locality_audit",
        {
            "datapoint_label": str(datapoint_label),
            "factor_class_id": str(factor_class_id),
            "threshold_policy_id": policy_id,
            "metrics": {
                "eta_q": omitted_displacement,
                "eta_g": omitted_gradient_fraction,
                "eta_H": cross_hessian_fraction,
                "environment_distance": float(environment_distance),
                "contact_regime_changed": bool(contact_regime_changed),
                "mapping_compatible": bool(mapping_compatible),
            },
        },
    )
    return LocalityAuditRecord(
        locality_audit_id=audit_id,
        datapoint_label=str(datapoint_label),
        factor_class_id=str(factor_class_id),
        omitted_displacement=omitted_displacement,
        omitted_gradient_fraction=omitted_gradient_fraction,
        cross_hessian_fraction=cross_hessian_fraction,
        environment_distance=float(environment_distance),
        contact_regime_changed=bool(contact_regime_changed),
        classification=classification,
        threshold_policy_id=policy_id,
        metrics=metrics,
    )


def classify_locality_metrics(
        omitted_displacement,
        omitted_gradient_fraction,
        cross_hessian_fraction,
        environment_distance,
        contact_regime_changed=False,
        boundary_localized=False,
        mapping_compatible=True,
        threshold_policy=None):
    """Classify persisted locality metrics with a versioned threshold policy."""

    policy = dict(DEFAULT_LOCALITY_THRESHOLD_POLICY)
    if threshold_policy is not None:
        policy.update(dict(threshold_policy))
    if not mapping_compatible:
        return "ambiguous"
    values = (
        float(omitted_displacement),
        float(omitted_gradient_fraction),
        float(cross_hessian_fraction),
        float(environment_distance),
    )
    if not all(math.isfinite(value) and value >= 0.0 for value in values):
        return "ambiguous"
    if contact_regime_changed or any((
        values[0] > float(policy["core_omitted_displacement"]),
        values[1] > float(policy["core_omitted_gradient_fraction"]),
        values[2] > float(policy["core_cross_hessian_fraction"]),
        values[3] > float(policy["core_environment_distance"]),
    )):
        return "core_candidate"
    exceeds_local = any((
        values[0] > float(policy["local_omitted_displacement"]),
        values[1] > float(policy["local_omitted_gradient_fraction"]),
        values[2] > float(policy["local_cross_hessian_fraction"]),
        values[3] > float(policy["local_environment_distance"]),
    ))
    if not exceeds_local:
        return "local"
    if boundary_localized:
        return "boundary_response"
    return "ambiguous"


def assign_core_candidate(
        audit,
        candidate_descriptor,
        core_families,
        descriptor_spec,
        source_core_family_id="",
        policy_id="schema4_promotion_v1"):
    """Decide whether a promoted candidate extends or creates a family."""

    if audit.classification != "core_candidate":
        decision = "no_action"
        target_family_id = ""
        reason_codes = (f"locality:{audit.classification}",)
        nearest_distance = 0.0
    else:
        candidate = np.asarray(candidate_descriptor, dtype=np.float64)
        compatible = [
            family for family in core_families
            if family.descriptor_spec_id == descriptor_spec.descriptor_spec_id
        ]
        distances = []
        for family in compatible:
            distance_squared = descriptor_squared_distance(
                descriptor_spec, candidate, family.prototype_descriptor
            )
            distances.append((math.sqrt(max(distance_squared, 0.0)), family))
        distances.sort(key=lambda item: (item[0], item[1].core_family_id))
        nearest_distance, nearest_family = (
            distances[0] if distances else (math.inf, None)
        )
        inside_family = (
            nearest_family is not None
            and nearest_distance <= nearest_family.descriptor_domain_radius
            and not audit.contact_regime_changed
        )
        if inside_family:
            decision = "add_core_point_to_existing_family"
            target_family_id = nearest_family.core_family_id
            reason_codes = ("inside_descriptor_domain", "outer_coverage_gap")
        else:
            decision = "create_core_family"
            target_family_id = ""
            reason_codes = (
                ("contact_regime_changed", "outside_descriptor_domain")
                if audit.contact_regime_changed else
                ("outside_descriptor_domain",)
            )
        if not math.isfinite(nearest_distance):
            nearest_distance = 0.0

    promotion_id = stable_msi_id(
        "promotion",
        {
            "datapoint_label": audit.datapoint_label,
            "locality_audit_id": audit.locality_audit_id,
            "decision": decision,
            "target_core_family_id": target_family_id,
            "policy_id": policy_id,
        },
    )
    return PromotionRecord(
        promotion_record_id=promotion_id,
        datapoint_label=audit.datapoint_label,
        locality_audit_id=audit.locality_audit_id,
        decision=decision,
        source_core_family_id=source_core_family_id,
        target_core_family_id=target_family_id,
        reason_codes=reason_codes,
        descriptor_distance=nearest_distance,
        policy_id=policy_id,
    )


def _shortest_graph_distance(connectivity, start, target, maximum=None):
    """Return graph-bond distance, stopping once ``maximum`` is exceeded."""

    if int(start) == int(target):
        return 0
    connectivity = np.asarray(connectivity)
    frontier = {int(start)}
    visited = set(frontier)
    distance = 0
    while frontier:
        distance += 1
        if maximum is not None and distance > int(maximum):
            return distance
        next_frontier = set()
        for atom in frontier:
            next_frontier.update(
                int(index) for index in np.flatnonzero(connectivity[atom])
                if int(index) not in visited
            )
        if int(target) in next_frontier:
            return distance
        visited.update(next_frontier)
        frontier = next_frontier
    return math.inf


def _smooth_contact_feature(coordinates, atom_i, atom_j, r0, exponent):
    if r0 <= 0.0 or exponent < 2 or exponent % 2:
        raise ValueError("Smooth contacts require r0 > 0 and an even exponent >= 2.")
    displacement = coordinates[atom_i] - coordinates[atom_j]
    distance = float(np.linalg.norm(displacement))
    jacobian = np.zeros_like(coordinates)
    if distance <= 1.0e-14:
        return 1.0, jacobian
    ratio = distance / float(r0)
    power = ratio**int(exponent)
    value = 1.0 / (1.0 + power)
    derivative = (
        -int(exponent)
        * ratio**(int(exponent) - 1)
        / (float(r0) * (1.0 + power)**2)
    )
    direction = displacement / distance
    jacobian[atom_i] = derivative * direction
    jacobian[atom_j] = -jacobian[atom_i]
    return float(value), jacobian


def _distance_pool_feature(coordinates, pairs, rms=False):
    if not pairs:
        raise ValueError("A distance pool requires at least one atom pair.")
    distances = []
    distance_jacobians = []
    for atom_i, atom_j in pairs:
        displacement = coordinates[atom_i] - coordinates[atom_j]
        distance = float(np.linalg.norm(displacement))
        if distance <= 1.0e-14:
            raise ValueError("Distance descriptor is undefined for coincident atoms.")
        jacobian = np.zeros_like(coordinates)
        jacobian[atom_i] = displacement / distance
        jacobian[atom_j] = -jacobian[atom_i]
        distances.append(distance)
        distance_jacobians.append(jacobian)
    distances = np.asarray(distances)
    jacobians = np.asarray(distance_jacobians)
    if not rms:
        return float(np.mean(distances)), np.mean(jacobians, axis=0)
    value = float(np.sqrt(np.mean(distances * distances)))
    weights = distances / (len(distances) * value)
    return value, np.tensordot(weights, jacobians, axes=1)


def _torsion_value_jacobian(coordinates, atoms):
    """Evaluate a torsion with a small analytical forward-mode vector algebra."""

    if len(atoms) != 4 or len(set(atoms)) != 4:
        raise ValueError("A torsion feature requires four distinct atoms.")
    ncart = coordinates.size
    dual_points = []
    for atom in atoms:
        value = coordinates[atom].copy()
        jacobian = np.zeros((3, ncart), dtype=np.float64)
        for component in range(3):
            jacobian[component, 3 * atom + component] = 1.0
        dual_points.append((value, jacobian))

    b0 = _dual_sub(dual_points[0], dual_points[1])
    b1 = _dual_sub(dual_points[2], dual_points[1])
    b2 = _dual_sub(dual_points[3], dual_points[2])
    axis = _dual_normalize(b1)
    v = _dual_sub(b0, _dual_scalar_vector(_dual_dot(b0, axis), axis))
    w = _dual_sub(b2, _dual_scalar_vector(_dual_dot(b2, axis), axis))
    x = _dual_dot(v, w)
    y = _dual_dot(_dual_cross(axis, v), w)
    denominator = x[0] * x[0] + y[0] * y[0]
    if denominator <= 1.0e-24:
        raise ValueError("Torsion descriptor is singular for collinear atoms.")
    value = float(np.arctan2(y[0], x[0]))
    derivative = (x[0] * y[1] - y[0] * x[1]) / denominator
    return value, derivative.reshape(coordinates.shape)


def _dual_sub(left, right):
    return left[0] - right[0], left[1] - right[1]


def _dual_dot(left, right):
    return (
        float(np.dot(left[0], right[0])),
        np.matmul(left[0], right[1]) + np.matmul(right[0], left[1]),
    )


def _dual_cross(left, right):
    value = np.cross(left[0], right[0])
    jacobian = np.column_stack([
        np.cross(left[1][:, index], right[0])
        + np.cross(left[0], right[1][:, index])
        for index in range(left[1].shape[1])
    ])
    return value, jacobian


def _dual_normalize(vector):
    norm = float(np.linalg.norm(vector[0]))
    if norm <= 1.0e-12:
        raise ValueError("Cannot normalize a zero-length descriptor vector.")
    value = vector[0] / norm
    projector = np.eye(3) / norm - np.outer(vector[0], vector[0]) / norm**3
    return value, np.matmul(projector, vector[1])


def _dual_scalar_vector(scalar, vector):
    value = scalar[0] * vector[0]
    jacobian = (
        vector[0][:, None] * scalar[1][None, :]
        + scalar[0] * vector[1]
    )
    return value, jacobian


def _cycle_atoms(connectivity):
    adjacency = {
        atom: {
            int(neighbor)
            for neighbor in np.flatnonzero(connectivity[atom])
        }
        for atom in range(connectivity.shape[0])
    }
    cycle_atoms = set()
    for atom_i, neighbors in adjacency.items():
        for atom_j in neighbors:
            if atom_j <= atom_i:
                continue
            stack = [atom_i]
            visited = {atom_i}
            parent = {atom_i: None}
            found = False
            path_endpoint = None
            while stack and not found:
                current = stack.pop()
                for neighbor in adjacency[current]:
                    if {current, neighbor} == {atom_i, atom_j}:
                        continue
                    if neighbor == atom_j:
                        found = True
                        path_endpoint = current
                        break
                    if neighbor not in visited:
                        visited.add(neighbor)
                        parent[neighbor] = current
                        stack.append(neighbor)
            if found:
                current = path_endpoint
                while current is not None:
                    cycle_atoms.add(current)
                    current = parent[current]
                cycle_atoms.add(atom_j)
    return cycle_atoms
