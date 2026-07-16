#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#

"""Validated schema-4 MSI registry I/O and legacy in-memory adaptation."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import json
import uuid

import h5py
import numpy as np

from .imcorefamilies import (
    CoreFamilySpec,
    CorePointSpec,
    EnvironmentDescriptorSpec,
    _FrozenRecord,
    _as_frozen_record,
    _dataclass_payload,
    stable_msi_id,
    validate_core_domain,
)
from .imlocalfactorlibrary import (
    FactorBindingSpec,
    FactorClassSpec,
    FactorDefinitionSpec,
    FactorValidationRecord,
    validate_factor_domain,
)


MSI_SCHEMA_VERSION = 4
MSI_SEMANTIC_VERSION = "multicore_shared_residual_v1"
MSI_REGISTRY_GROUP = "msi_model_registry"

_JSON_DATASETS = (
    "model_info_json",
    "descriptor_specs_json",
    "core_families_json",
    "core_points_json",
    "factor_classes_json",
    "factor_definitions_json",
    "factor_bindings_json",
    "validation_records_json",
    "datapoint_index_json",
)

_MODEL_INFO_FIELDS = (
    "schema_version",
    "semantic_version",
    "root",
    "coordinate_conventions",
    "unit_conventions",
    "global_qm_provenance_fingerprint",
    "outer_distance_policy",
    "cardinality_contract",
    "created_by",
)


@dataclass(frozen=True)
class MSIModelRegistry:
    """Immutable, JSON-safe schema-4 model registry."""

    model_info: Mapping
    descriptor_specs: tuple[EnvironmentDescriptorSpec, ...]
    core_families: tuple[CoreFamilySpec, ...]
    core_points: tuple[CorePointSpec, ...]
    factor_classes: tuple[FactorClassSpec, ...] = ()
    factor_definitions: tuple[FactorDefinitionSpec, ...] = ()
    factor_bindings: tuple[FactorBindingSpec, ...] = ()
    validation_records: tuple[FactorValidationRecord, ...] = ()
    datapoint_index: Mapping = _FrozenRecord()

    def __post_init__(self):
        object.__setattr__(
            self,
            "model_info",
            _as_frozen_record(self.model_info, "model_info"),
        )
        for field_name in (
            "descriptor_specs",
            "core_families",
            "core_points",
            "factor_classes",
            "factor_definitions",
            "factor_bindings",
            "validation_records",
        ):
            object.__setattr__(self, field_name, tuple(getattr(self, field_name)))
        object.__setattr__(
            self,
            "datapoint_index",
            _as_frozen_record(self.datapoint_index, "datapoint_index"),
        )

    @property
    def root(self):
        return str(self.model_info["root"])

    def to_payloads(self):
        return {
            "model_info_json": dict(_dataclass_payload(self)["model_info"]),
            "descriptor_specs_json": [
                item.to_payload() for item in self.descriptor_specs
            ],
            "core_families_json": [
                item.to_payload() for item in self.core_families
            ],
            "core_points_json": [item.to_payload() for item in self.core_points],
            "factor_classes_json": [
                item.to_payload() for item in self.factor_classes
            ],
            "factor_definitions_json": [
                item.to_payload() for item in self.factor_definitions
            ],
            "factor_bindings_json": [
                item.to_payload() for item in self.factor_bindings
            ],
            "validation_records_json": [
                item.to_payload() for item in self.validation_records
            ],
            "datapoint_index_json": dict(
                _dataclass_payload(self)["datapoint_index"]
            ),
        }

    @classmethod
    def from_payloads(cls, payloads):
        return cls(
            model_info=payloads["model_info_json"],
            descriptor_specs=tuple(
                EnvironmentDescriptorSpec.from_payload(item)
                for item in payloads["descriptor_specs_json"]
            ),
            core_families=tuple(
                CoreFamilySpec.from_payload(item)
                for item in payloads["core_families_json"]
            ),
            core_points=tuple(
                CorePointSpec.from_payload(item)
                for item in payloads["core_points_json"]
            ),
            factor_classes=tuple(
                FactorClassSpec.from_payload(item)
                for item in payloads["factor_classes_json"]
            ),
            factor_definitions=tuple(
                FactorDefinitionSpec.from_payload(item)
                for item in payloads["factor_definitions_json"]
            ),
            factor_bindings=tuple(
                FactorBindingSpec.from_payload(item)
                for item in payloads["factor_bindings_json"]
            ),
            validation_records=tuple(
                FactorValidationRecord.from_payload(item)
                for item in payloads["validation_records_json"]
            ),
            datapoint_index=payloads["datapoint_index_json"],
        )


def validate_msi_registry(registry):
    """Validate a complete schema-4 model before write and after load."""

    missing_info = [
        field for field in _MODEL_INFO_FIELDS if field not in registry.model_info
    ]
    if missing_info:
        raise ValueError(
            "Schema-4 model_info is missing required fields: "
            + ", ".join(missing_info)
        )
    if int(registry.model_info["schema_version"]) != MSI_SCHEMA_VERSION:
        raise ValueError("MSI model registry schema_version must be 4.")
    if str(registry.model_info["semantic_version"]) != MSI_SEMANTIC_VERSION:
        raise ValueError(
            "Unsupported schema-4 MSI semantic_version "
            f"{registry.model_info['semantic_version']!r}."
        )
    if not registry.root:
        raise ValueError("Schema-4 model root must not be empty.")

    validate_core_domain(
        registry.descriptor_specs,
        registry.core_families,
        registry.core_points,
    )
    family_ids = tuple(item.core_family_id for item in registry.core_families)
    point_ids = tuple(item.core_point_id for item in registry.core_points)
    validate_factor_domain(
        registry.factor_classes,
        registry.factor_definitions,
        registry.factor_bindings,
        registry.validation_records,
        family_ids,
        point_ids,
    )

    family_by_id = {
        family.core_family_id: family for family in registry.core_families
    }
    binding_by_id = {
        binding.factor_binding_id: binding
        for binding in registry.factor_bindings
    }
    for family in registry.core_families:
        for binding_id in family.factor_binding_ids:
            binding = binding_by_id.get(binding_id)
            if binding is None:
                raise ValueError(
                    f"Core family {family.core_family_id!r} references missing "
                    f"binding {binding_id!r}."
                )
            if binding.target_core_family_id != family.core_family_id:
                raise ValueError(
                    f"Core family {family.core_family_id!r} references binding "
                    f"{binding_id!r} owned by another family."
                )

    for point in registry.core_points:
        if point.core_family_id not in family_by_id:
            continue
        for binding_id in point.factor_binding_ids:
            binding = binding_by_id.get(binding_id)
            if binding is None:
                raise ValueError(
                    f"Core point {point.core_point_id!r} references missing "
                    f"binding {binding_id!r}."
                )
            if binding.target_core_family_id != point.core_family_id:
                raise ValueError(
                    f"Core point {point.core_point_id!r} references a binding "
                    "for another family."
                )

    index = registry.datapoint_index
    for point in registry.core_points:
        record = index.get(point.datapoint_label)
        if record is None:
            raise ValueError(
                f"Core point {point.core_point_id!r} references missing datapoint "
                f"label {point.datapoint_label!r}."
            )
        roles = _datapoint_model_roles(record)
        if not roles.intersection(("outer_core", "outer_global")):
            raise ValueError(
                f"Core datapoint {point.datapoint_label!r} is not indexed as an "
                "outer point."
            )

    for definition in registry.factor_definitions:
        for label in definition.state_datapoint_labels:
            record = index.get(label)
            if record is None:
                raise ValueError(
                    f"Factor definition {definition.factor_definition_id!r} "
                    f"references missing state label {label!r}."
                )
            if "residual_state" not in _datapoint_model_roles(record):
                raise ValueError(
                    f"Factor state {label!r} is not indexed as residual_state."
                )


def _datapoint_model_roles(record):
    """Return all schema-4 roles, including promotion-time dual ownership."""

    roles = record.get("model_roles", ())
    if isinstance(roles, str):
        roles = (roles,)
    result = {str(role) for role in roles if str(role)}
    role = str(record.get("model_role", ""))
    if role:
        result.add(role)
    return frozenset(result)


def has_msi_registry(filename, root):
    """Return whether a complete active schema-4 group exists for ``root``."""

    path = f"{MSI_REGISTRY_GROUP}/root_{root}"
    with h5py.File(filename, "r") as h5f:
        return path in h5f and all(
            name in h5f[path] for name in _JSON_DATASETS
        )


def list_msi_registry_roots(filename):
    """List complete active roots while ignoring staging and backup groups."""

    with h5py.File(filename, "r") as h5f:
        if MSI_REGISTRY_GROUP not in h5f:
            return ()
        base = h5f[MSI_REGISTRY_GROUP]
        roots = []
        for name, group in base.items():
            if not name.startswith("root_"):
                continue
            if all(dataset in group for dataset in _JSON_DATASETS):
                roots.append(name.replace("root_", "", 1))
        return tuple(sorted(roots, key=str))


def write_msi_registry_atomic(filename, registry):
    """Stage, re-read, validate, and replace one active schema-4 root group."""

    validate_msi_registry(registry)
    payloads = registry.to_payloads()
    root_name = f"root_{registry.root}"
    token = uuid.uuid4().hex
    staging_name = f"__staging_{root_name}_{token}"
    backup_name = f"__backup_{root_name}_{token}"

    with h5py.File(filename, "a") as h5f:
        base = h5f.require_group(MSI_REGISTRY_GROUP)
        staging = base.create_group(staging_name)
        try:
            for dataset_name in _JSON_DATASETS:
                _write_json_dataset(
                    staging, dataset_name, payloads[dataset_name]
                )
            h5f.flush()
            staged_registry = _registry_from_group(staging)
            validate_msi_registry(staged_registry)
            if staged_registry != registry:
                raise RuntimeError(
                    "Schema-4 staged registry changed during JSON round-trip."
                )

            had_active = root_name in base
            if had_active:
                base.move(root_name, backup_name)
            try:
                base.move(staging_name, root_name)
            except Exception:
                if had_active and backup_name in base and root_name not in base:
                    base.move(backup_name, root_name)
                raise
            if backup_name in base:
                del base[backup_name]
            h5f.flush()
        except Exception:
            if staging_name in base:
                del base[staging_name]
            if backup_name in base and root_name not in base:
                base.move(backup_name, root_name)
            raise


def read_msi_registry(filename, root):
    """Read and validate only the active root group, ignoring staging groups."""

    path = f"{MSI_REGISTRY_GROUP}/root_{root}"
    with h5py.File(filename, "r") as h5f:
        if path not in h5f:
            raise KeyError(f"No active schema-4 MSI registry at /{path}.")
        registry = _registry_from_group(h5f[path])
    validate_msi_registry(registry)
    if registry.root != str(root):
        raise ValueError(
            f"Schema-4 registry root {registry.root!r} does not match requested "
            f"root {str(root)!r}."
        )
    return registry


def _write_json_dataset(group, name, payload):
    data = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )
    string_type = h5py.string_dtype(encoding="utf-8")
    group.create_dataset(
        name, data=np.array(data, dtype=object), dtype=string_type
    )


def _read_json_dataset(group, name):
    if name not in group:
        raise ValueError(
            f"Schema-4 registry group {group.name!r} is missing {name!r}."
        )
    raw = group[name][()]
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    return json.loads(str(raw))


def _registry_from_group(group):
    return MSIModelRegistry.from_payloads(
        {name: _read_json_dataset(group, name) for name in _JSON_DATASETS}
    )


def adapt_legacy_banks_to_schema4(legacy_banks, root):
    """Create native-only schema-4 entities without inferring factor sharing."""

    root = str(root)
    descriptor = EnvironmentDescriptorSpec(
        descriptor_spec_id=f"legacy.descriptor.root_{root}",
        version="legacy_adapter_v1",
        feature_definitions=(),
        feature_scales=(),
        distance_coefficient=0.0,
        invariance_policy="legacy_outer_distance_only",
        descriptor_kind="legacy_none",
        provenance={"migration_status": "unvalidated_native_adapter"},
    )

    families = []
    points = []
    classes = []
    definitions = []
    bindings = []
    datapoint_index = {}

    for family_label, family_bank in sorted(
            (legacy_banks or {}).items(), key=lambda item: str(item[0])):
        family_label = str(family_label)
        family_id = stable_msi_id(
            "family", {"root": root, "legacy_family_label": family_label}
        )
        core = family_bank.get("core")
        core_label = _datapoint_label(
            core,
            family_bank.get("point_index", {}).get("core_label"),
            fallback=f"{family_label}_core",
        )
        point_id = stable_msi_id(
            "core_point", {"family": family_id, "label": core_label}
        )
        datapoint_index[core_label] = {
            "model_role": "outer_core",
            "full_pes_cardinal": True,
            "residual_cardinal": False,
            "legacy_family_label": family_label,
            "core_family_id": family_id,
            "core_point_id": point_id,
        }

        family_binding_ids = []
        for term_id, term_bank in sorted(
                family_bank.get("terms", {}).items(),
                key=lambda item: str(item[0])):
            term_id = str(term_id)
            class_id = stable_msi_id(
                "factor_class",
                {
                    "root": root,
                    "legacy_family_label": family_label,
                    "legacy_term_id": term_id,
                },
            )
            active_atoms = tuple(
                int(atom) for atom in term_bank.get("active_atoms", ())
            )
            owned_roles = tuple(
                (f"legacy_atom_{index}", atom)
                for index, atom in enumerate(active_atoms)
            )
            signature_rows = tuple(
                int(row)
                for rows in term_bank.get("grouped_signature_rows", ())
                for row in rows
            )
            active_rows = tuple(
                int(row) for row in term_bank.get("active_rows", ())
            )
            projector_rows = tuple(
                int(row) for row in term_bank.get("projector_rows", ())
            )
            transport_rows = tuple(
                sorted(set(active_rows) | set(signature_rows) | set(projector_rows))
            )
            classes.append(
                FactorClassSpec(
                    factor_class_id=class_id,
                    kind=str(term_bank.get("role", "legacy_factor")),
                    canonical_graph_fingerprint=(
                        f"legacy:{root}:{family_label}:{term_id}"
                    ),
                    owned_atom_roles=owned_roles,
                    rotor_axes=(),
                    symmetry_orders=(),
                    canonical_coordinate_signature=tuple(
                        f"legacy_row:{row}" for row in signature_rows
                    ),
                    default_support_policy="legacy_projector_unchanged",
                    default_relaxation_policy=str(
                        term_bank.get("relaxation_policy_id", "legacy")
                    ),
                    provenance={
                        "legacy_family_label": family_label,
                        "legacy_term_id": term_id,
                    },
                )
            )

            state_items = sorted(
                term_bank.get("expected_states", {}).items(),
                key=lambda item: str(item[0]),
            )
            state_labels = tuple(
                _datapoint_label(
                    datapoint,
                    fallback=f"{family_label}_{term_id}_state_{state_id}",
                )
                for state_id, datapoint in state_items
            )
            if not state_labels:
                raise ValueError(
                    f"Legacy term {family_label!r}/{term_id!r} has no states."
                )
            for (state_id, _), label in zip(state_items, state_labels):
                datapoint_index[label] = {
                    "model_role": "residual_state",
                    "full_pes_cardinal": False,
                    "residual_cardinal": True,
                    "legacy_family_label": family_label,
                    "legacy_term_id": term_id,
                    "legacy_state_id": str(state_id),
                    "factor_class_id": class_id,
                }

            definition_id = stable_msi_id(
                "factor_definition",
                {"family": family_id, "term": term_id, "states": state_labels},
            )
            factor_info = family_bank.get("factor_info", {})
            fingerprint = hashlib.sha256(
                json.dumps(
                    factor_info,
                    sort_keys=True,
                    separators=(",", ":"),
                    default=str,
                ).encode("utf-8")
            ).hexdigest()
            anchor_ids = tuple(
                str(state_id)
                for state_id in term_bank.get("anchor_state_ids", ())
            )
            if not anchor_ids:
                anchor_ids = (str(state_items[0][0]),)
            definitions.append(
                FactorDefinitionSpec(
                    factor_definition_id=definition_id,
                    factor_class_id=class_id,
                    source_core_family_id=family_id,
                    source_core_point_id=point_id,
                    state_datapoint_labels=state_labels,
                    anchor_state_ids=anchor_ids,
                    residual_derivative_order=2,
                    projector_policy=str(
                        term_bank.get("projector_policy_id", "legacy_unchanged")
                    ),
                    support_policy="legacy_unchanged",
                    symmetry_chart_policy="legacy_mapping_masks_unchanged",
                    coordinate_qm_fingerprint=fingerprint,
                    provenance={
                        "migration_status": "unvalidated_native_adapter",
                        "legacy_schema_version": int(
                            factor_info.get("schema_version", 0)
                        ),
                        "legacy_family_label": family_label,
                        "legacy_term_id": term_id,
                        "state_id_to_label": {
                            str(state_id): label
                            for (state_id, _), label in zip(
                                state_items, state_labels
                            )
                        },
                    },
                )
            )
            binding_id = stable_msi_id(
                "factor_binding",
                {
                    "target_family": family_id,
                    "definition": definition_id,
                    "mode": "native",
                },
            )
            family_binding_ids.append(binding_id)
            bindings.append(
                FactorBindingSpec(
                    factor_binding_id=binding_id,
                    target_core_family_id=family_id,
                    factor_class_id=class_id,
                    factor_definition_id=definition_id,
                    mode="native",
                    atom_transport_map=tuple(
                        (atom, atom) for atom in active_atoms
                    ),
                    row_transport_map=tuple(
                        (row, row) for row in transport_rows
                    ),
                    target_anchor_mapping=tuple(
                        (anchor_id, anchor_id) for anchor_id in anchor_ids
                    ),
                    applicability_domain={"policy": "legacy_static_family"},
                    provenance={
                        "migration_status": "unvalidated_native_adapter"
                    },
                )
            )

        families.append(
            CoreFamilySpec(
                core_family_id=family_id,
                root=root,
                descriptor_spec_id=descriptor.descriptor_spec_id,
                prototype_descriptor=(),
                descriptor_domain_radius=0.0,
                factor_binding_ids=tuple(family_binding_ids),
                construction_policy_id="legacy_native_adapter",
                provenance={
                    "legacy_family_label": family_label,
                    "migration_status": "unvalidated_native_adapter",
                },
            )
        )
        points.append(
            CorePointSpec(
                core_point_id=point_id,
                core_family_id=family_id,
                datapoint_label=core_label,
                reference_descriptor=(),
                confidence_radius=float(
                    getattr(core, "confidence_radius", 0.0) or 0.0
                ),
                full_pes_cardinal=True,
                second_order_cardinal=True,
                provenance={"legacy_family_label": family_label},
            )
        )

    registry = MSIModelRegistry(
        model_info={
            "schema_version": MSI_SCHEMA_VERSION,
            "semantic_version": MSI_SEMANTIC_VERSION,
            "root": root,
            "coordinate_conventions": "legacy_unchanged",
            "unit_conventions": "legacy_atomic_units",
            "global_qm_provenance_fingerprint": "legacy_unavailable",
            "outer_distance_policy": "legacy_unchanged",
            "cardinality_contract": {
                "outer_core": "full_pes_cardinal",
                "factor_state": "residual_cardinal",
            },
            "created_by": "schema2_3_native_adapter_v1",
            "migration_status": "unvalidated_native_adapter",
        },
        descriptor_specs=(descriptor,),
        core_families=tuple(families),
        core_points=tuple(points),
        factor_classes=tuple(classes),
        factor_definitions=tuple(definitions),
        factor_bindings=tuple(bindings),
        validation_records=(),
        datapoint_index=datapoint_index,
    )
    validate_msi_registry(registry)
    return registry


def adapt_legacy_file_to_schema4(
        filename, root, z_matrix, interpolation_settings):
    """Load schema 2/3 through its unchanged reader and adapt only in memory."""

    from .imlocalfactorregistry import load_signed_factor_banks_for_root

    banks = load_signed_factor_banks_for_root(
        filename, root, z_matrix, interpolation_settings
    )
    return adapt_legacy_banks_to_schema4(banks, root)


def _datapoint_label(datapoint, candidate=None, fallback=None):
    for value in (
        getattr(datapoint, "point_label", None),
        candidate,
        fallback,
    ):
        if value is not None and str(value).strip():
            return str(value)
    raise ValueError("A legacy datapoint has no stable label.")
