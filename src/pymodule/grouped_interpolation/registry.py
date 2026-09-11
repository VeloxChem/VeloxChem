"""Versioned, restartable HDF5 registry for grouped interpolation models."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json

import h5py
import numpy as np

from .coordinates import canonical_json
from .models import GroupedCalculationResult


GROUPED_SCHEMA_VERSION = "1.1.0"


def _read_text(dataset) -> str:
    value = dataset[()]
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _write_text(group, name: str, value: str) -> None:
    temporary = f"__new_{name}"
    if temporary in group:
        del group[temporary]
    group.create_dataset(temporary, data=value, dtype=h5py.string_dtype("utf-8"))
    if name in group:
        del group[name]
    group.move(temporary, name)


def _write_json(group, name: str, value) -> None:
    _write_text(group, name, canonical_json(value))


def _read_json(group, name: str):
    return json.loads(_read_text(group[name]))


class GroupedModelRegistry:
    """Own grouped metadata without relying on point-name suffix semantics."""

    def __init__(self, filename, model_id):
        self.filename = str(filename)
        self.model_id = str(model_id)

    @property
    def path(self):
        return f"grouped_models/{self.model_id}"

    def initialize(self, *, plan, coordinate_definition, rotors, policy, metadata):
        with h5py.File(self.filename, "a") as h5file:
            models = h5file.require_group("grouped_models")
            if self.model_id in models:
                group = models[self.model_id]
                existing = _read_json(group, "coordinate_definition_json")
                if existing["fingerprint"] != coordinate_definition.fingerprint:
                    raise ValueError(
                        "Cannot resume grouped model with a different coordinate chart."
                    )
                return _read_text(group["status"])

            group = models.create_group(self.model_id)
            model_metadata = dict(metadata)
            model_metadata.update(
                {
                    "schema_version": GROUPED_SCHEMA_VERSION,
                    "model_id": self.model_id,
                    "family_id": plan.family_id,
                    "coordinate_fingerprint": coordinate_definition.fingerprint,
                    "construction_mode": "grouped_symmetry_aware",
                    "construction_policy": "signed_rigid_residual",
                    "created_at_utc": datetime.now(timezone.utc).isoformat(),
                }
            )
            _write_json(group, "metadata_json", model_metadata)
            _write_json(group, "coordinate_definition_json", asdict(coordinate_definition))
            _write_json(group, "motion_registry_json", [asdict(rotor) for rotor in rotors])
            _write_json(group, "symmetry_registry_json", [])
            _write_json(group, "coupling_registry_json", [])
            _write_json(group, "interaction_registry_json", [])
            _write_json(group, "subset_registry_json", [asdict(subset) for subset in plan.subsets])
            manifest = [
                {"request": asdict(request), "status": "pending"}
                for request in plan.requests
            ]
            _write_json(group, "request_manifest_json", manifest)
            _write_json(group, "point_index_json", {})
            _write_json(group, "validation_json", {"status": "pending"})
            _write_json(group, "policy_json", asdict(policy))
            _write_json(
                group,
                "factor_registry_json",
                {
                    subset_id: list(request_ids)
                    for subset_id, request_ids in plan.factor_state_requests
                },
            )
            _write_text(group, "status", "staging")
            _write_text(group, "content_hash", "")
            group.create_group("points")
            h5file.flush()
        return "staging"

    def status(self) -> str:
        with h5py.File(self.filename, "r") as h5file:
            return _read_text(h5file[self.path]["status"])

    def manifest(self):
        with h5py.File(self.filename, "r") as h5file:
            return _read_json(h5file[self.path], "request_manifest_json")

    def point_index(self):
        with h5py.File(self.filename, "r") as h5file:
            return _read_json(h5file[self.path], "point_index_json")

    def request_is_complete(self, request_id: str) -> bool:
        return any(
            item["request"]["request_id"] == request_id and item["status"] == "complete"
            for item in self.manifest()
        )

    def write_result(self, request, result: GroupedCalculationResult) -> None:
        if result.coordinate_fingerprint != request.coordinate_fingerprint:
            raise ValueError("Refusing to store a result from a different coordinate chart.")
        if result.request_id != request.request_id:
            raise ValueError("Grouped calculation result/request IDs do not agree.")
        with h5py.File(self.filename, "a") as h5file:
            group = h5file[self.path]
            points = group["points"]
            if request.request_id in points:
                del points[request.request_id]
            point_group = points.create_group(request.request_id)
            point_group.attrs["coordinate_fingerprint"] = result.coordinate_fingerprint
            point_group.attrs["is_virtual_image"] = bool(request.is_virtual_image)
            point_group.attrs["converged"] = bool(result.converged)
            point_group.create_dataset("energy", data=np.float64(result.energy))
            for name, value in (
                ("cartesian_coordinates", result.cartesian_coordinates),
                ("cartesian_gradient", result.cartesian_gradient),
                ("internal_coordinates", result.internal_coordinates),
                ("internal_gradient", result.internal_gradient),
                ("measured_phase_signature", result.measured_phase_signature),
                ("constraint_errors", result.constraint_errors),
            ):
                point_group.create_dataset(name, data=np.asarray(value, dtype=float), compression="gzip" if np.asarray(value).size else None)
            if result.cartesian_hessian is not None:
                point_group.create_dataset(
                    "cartesian_hessian",
                    data=np.asarray(result.cartesian_hessian, dtype=float),
                    compression="gzip",
                )
            if result.internal_hessian is not None:
                point_group.create_dataset(
                    "internal_hessian",
                    data=np.asarray(result.internal_hessian, dtype=float),
                    compression="gzip",
                )
            _write_text(point_group, "source_point_label", result.source_point_label)
            _write_text(point_group, "provenance", result.provenance)

            manifest = _read_json(group, "request_manifest_json")
            found = False
            for item in manifest:
                if item["request"]["request_id"] == request.request_id:
                    item["status"] = "complete"
                    found = True
                    break
            if not found:
                manifest.append({"request": asdict(request), "status": "complete"})
            point_index = _read_json(group, "point_index_json")
            point_index[request.request_id] = f"{self.path}/points/{request.request_id}"
            _write_json(group, "request_manifest_json", manifest)
            _write_json(group, "point_index_json", point_index)
            h5file.flush()

    def read_result(self, request_id: str) -> dict:
        with h5py.File(self.filename, "r") as h5file:
            point = h5file[self.path]["points"][request_id]
            result = {
                "request_id": request_id,
                "energy": float(point["energy"][()]),
                "coordinate_fingerprint": str(point.attrs["coordinate_fingerprint"]),
                "is_virtual_image": bool(point.attrs["is_virtual_image"]),
                "converged": bool(point.attrs["converged"]),
                "source_point_label": _read_text(point["source_point_label"]),
                "provenance": _read_text(point["provenance"]),
            }
            for name in (
                "cartesian_coordinates",
                "cartesian_gradient",
                "cartesian_hessian",
                "internal_coordinates",
                "internal_gradient",
                "internal_hessian",
                "measured_phase_signature",
                "constraint_errors",
            ):
                result[name] = None if name not in point else np.asarray(point[name], dtype=float)
            return result

    def write_symmetry_registry(self, operations) -> None:
        with h5py.File(self.filename, "a") as h5file:
            _write_json(
                h5file[self.path],
                "symmetry_registry_json",
                [asdict(operation) for operation in operations],
            )

    def write_coupling_registry(self, coupling_records) -> None:
        with h5py.File(self.filename, "a") as h5file:
            _write_json(
                h5file[self.path],
                "coupling_registry_json",
                coupling_records,
            )

    def write_interaction_registry(self, interaction_records) -> None:
        with h5py.File(self.filename, "a") as h5file:
            _write_json(
                h5file[self.path],
                "interaction_registry_json",
                interaction_records,
            )

    def write_factor_registry(self, factor_registry) -> None:
        with h5py.File(self.filename, "a") as h5file:
            _write_json(h5file[self.path], "factor_registry_json", factor_registry)

    def write_validation(self, validation) -> None:
        with h5py.File(self.filename, "a") as h5file:
            _write_json(h5file[self.path], "validation_json", validation)

    def update_metadata(self, updates) -> None:
        with h5py.File(self.filename, "a") as h5file:
            group = h5file[self.path]
            metadata = _read_json(group, "metadata_json")
            metadata.update(dict(updates))
            _write_json(group, "metadata_json", metadata)

    def _calculate_content_hash(self, group) -> str:
        digest = hashlib.sha256()
        for name in (
            "metadata_json",
            "coordinate_definition_json",
            "motion_registry_json",
            "symmetry_registry_json",
            "coupling_registry_json",
            "interaction_registry_json",
            "subset_registry_json",
            "request_manifest_json",
            "point_index_json",
            "validation_json",
            "policy_json",
            "factor_registry_json",
        ):
            if name not in group:
                continue
            digest.update(name.encode("utf-8"))
            if name == "validation_json":
                validation = _read_json(group, name)
                validation.pop("content_hash", None)
                digest.update(canonical_json(validation).encode("utf-8"))
            elif name == "metadata_json":
                metadata = _read_json(group, name)
                metadata.pop("created_at_utc", None)
                digest.update(canonical_json(metadata).encode("utf-8"))
            else:
                digest.update(_read_text(group[name]).encode("utf-8"))
        for request_id in sorted(group["points"].keys()):
            point = group["points"][request_id]
            digest.update(request_id.encode("utf-8"))
            for name in sorted(point.keys()):
                value = point[name][()]
                if isinstance(value, bytes):
                    digest.update(value)
                else:
                    digest.update(np.asarray(value).tobytes())
        return digest.hexdigest()

    def publish(self) -> str:
        with h5py.File(self.filename, "a") as h5file:
            group = h5file[self.path]
            manifest = _read_json(group, "request_manifest_json")
            incomplete = [
                item["request"]["request_id"]
                for item in manifest
                if item["status"] != "complete"
            ]
            if incomplete:
                raise RuntimeError(f"Cannot publish with incomplete requests: {incomplete}")
            validation = _read_json(group, "validation_json")
            if validation.get("status") != "passed":
                raise RuntimeError("Cannot publish a grouped model that did not pass validation.")
            content_hash = self._calculate_content_hash(group)
            validation["content_hash"] = content_hash
            _write_json(group, "validation_json", validation)
            _write_text(group, "content_hash", content_hash)
            # Publication state is deliberately the final write.
            _write_text(group, "status", "published")
            h5file.flush()
            return content_hash

    def load_bundle(self, *, require_published=True):
        with h5py.File(self.filename, "r") as h5file:
            group = h5file[self.path]
            status = _read_text(group["status"])
            if require_published and status != "published":
                raise RuntimeError(
                    f"Grouped model {self.model_id} is {status}, not published."
                )
            point_ids = sorted(group["points"].keys())
            json_names = (
                    "metadata_json",
                    "coordinate_definition_json",
                    "motion_registry_json",
                    "symmetry_registry_json",
                    "coupling_registry_json",
                    "interaction_registry_json",
                    "subset_registry_json",
                    "request_manifest_json",
                    "point_index_json",
                    "validation_json",
                    "policy_json",
                    "factor_registry_json",
                )
            bundle = {
                name: _read_json(group, name)
                for name in json_names
                if name in group
            }
            bundle.setdefault("coupling_registry_json", [])
            bundle.setdefault("interaction_registry_json", [])
            bundle["status"] = status
            bundle["content_hash"] = _read_text(group["content_hash"])
            bundle["points"] = {}
            for request_id in point_ids:
                point = group["points"][request_id]
                bundle["points"][request_id] = {
                    "request_id": request_id,
                    "energy": float(point["energy"][()]),
                    "internal_coordinates": np.asarray(point["internal_coordinates"], dtype=float),
                    "internal_gradient": np.asarray(point["internal_gradient"], dtype=float),
                    "internal_hessian": None
                    if "internal_hessian" not in point
                    else np.asarray(point["internal_hessian"], dtype=float),
                    "cartesian_coordinates": np.asarray(point["cartesian_coordinates"], dtype=float),
                    "cartesian_gradient": np.asarray(point["cartesian_gradient"], dtype=float),
                    "is_virtual_image": bool(point.attrs["is_virtual_image"]),
                    "coordinate_fingerprint": str(point.attrs["coordinate_fingerprint"]),
                }
            return bundle


def find_published_model(filename, model_id=None):
    with h5py.File(filename, "r") as h5file:
        if "grouped_models" not in h5file:
            raise KeyError(f"No grouped_models registry exists in {filename}.")
        models = h5file["grouped_models"]
        candidates = []
        for candidate in models:
            status = _read_text(models[candidate]["status"])
            if status == "published" and (model_id is None or candidate == model_id):
                candidates.append(candidate)
        if not candidates:
            requested = "any model" if model_id is None else model_id
            raise KeyError(f"No published grouped model matching {requested} exists in {filename}.")
        if model_id is None and len(candidates) > 1:
            raise ValueError(
                "Multiple published grouped models exist; set grouped_model_id explicitly."
            )
        return candidates[0]


def read_grouped_coordinate_z_matrix(filename, model_id=None):
    """Read the canonical chart from one grouped model, including staging models."""

    with h5py.File(filename, "r") as h5file:
        if "grouped_models" not in h5file:
            raise KeyError(f"No grouped_models registry exists in {filename}.")
        models = h5file["grouped_models"]
        candidates = [
            candidate for candidate in models
            if model_id is None or candidate == model_id
        ]
        if len(candidates) != 1:
            raise ValueError(
                "A unique grouped model is required to recover its coordinate chart."
            )
        definition = _read_json(
            models[candidates[0]], "coordinate_definition_json"
        )
        z_matrix = {"bonds": [], "angles": [], "dihedrals": [], "impropers": []}
        section_by_type = {
            "bond": "bonds",
            "angle": "angles",
            "proper_torsion": "dihedrals",
            "improper_torsion": "impropers",
        }
        for row, row_type in zip(definition["rows"], definition["row_types"]):
            z_matrix[section_by_type[row_type]].append(tuple(int(atom) for atom in row))
        return candidates[0], z_matrix
