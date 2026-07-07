#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json

import h5py
import numpy as np

from .interpolationdatapoint import InterpolationDatapoint


def _restore_id(value):
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return str(value)


def _json_key(value):
    return str(value)


def _read_scalar_string(ds):
    val = ds[()]
    if isinstance(val, bytes):
        return val.decode("utf-8")
    return str(val)


def _write_string_dataset(group, name, value):
    if name in group:
        del group[name]
    dt = h5py.string_dtype(encoding="utf-8")
    group.create_dataset(name, data=np.array(value, dtype=object), dtype=dt)


def _as_float_list(value):
    if value is None:
        return []
    return [float(x) for x in np.asarray(value, dtype=np.float64).reshape(-1)]


def _as_int_tuple(value):
    if value is None:
        return ()
    return tuple(int(x) for x in np.asarray(value, dtype=np.int64).reshape(-1))


def _rotor_key(rotor_ids):
    return tuple(sorted((str(rotor_id) for rotor_id in rotor_ids), key=str))


def _phase_for_record(record, n_rotors):
    phase = np.asarray(
        record.get("phase_signature", ()),
        dtype=np.float64,
    ).reshape(-1)
    out = np.zeros(int(n_rotors), dtype=np.float64)
    ncopy = min(out.size, phase.size)
    if ncopy:
        out[:ncopy] = phase[:ncopy]
    return out


def _rotor_rows(local_group_model, rotor_id):
    rotor = local_group_model.rotors.get(str(rotor_id))
    if rotor is None:
        return ()
    signature_rows = getattr(rotor, "signature_rows", ())
    if signature_rows:
        return tuple(int(row) for row in signature_rows)
    return tuple(int(row) for row in rotor.torsion_rows)


def _rotor_row_types(local_group_model, rotor_id):
    rotor = local_group_model.rotors.get(str(rotor_id))
    if rotor is None:
        return ()

    rows = _rotor_rows(local_group_model, rotor_id)
    row_types = tuple(
        str(row_type).strip().lower()
        for row_type in getattr(rotor, "signature_row_types", ())
    )
    if len(row_types) == len(rows):
        return row_types

    return tuple("torsion" for _ in rows)


def _rotor_row_scales(local_group_model, rotor_id):
    rotor = local_group_model.rotors.get(str(rotor_id))
    if rotor is None:
        return ()

    rows = _rotor_rows(local_group_model, rotor_id)
    row_scales = tuple(
        float(scale)
        for scale in getattr(rotor, "signature_row_scales", ())
    )
    if len(row_scales) == len(rows):
        return row_scales

    return tuple(
        0.35 if row_type == "angle" else 1.0
        for row_type in _rotor_row_types(local_group_model, rotor_id)
    )


def _signature_rows(local_group_model, rotor_ids):
    return tuple(
        tuple(int(row) for row in _rotor_rows(local_group_model, rotor_id))
        for rotor_id in rotor_ids
    )


def _signature_row_types(local_group_model, rotor_ids):
    return tuple(
        tuple(str(row_type) for row_type in _rotor_row_types(
            local_group_model, rotor_id))
        for rotor_id in rotor_ids
    )


def _signature_row_scales(local_group_model, rotor_ids):
    return tuple(
        tuple(float(scale) for scale in _rotor_row_scales(
            local_group_model, rotor_id))
        for rotor_id in rotor_ids
    )


def _active_rows_for_term(local_group_model, rotor_ids, cluster=None):
    rows = []
    if cluster is not None:
        rows.extend(int(row) for row in getattr(cluster, "active_rows", ()))
    if not rows:
        for rotor_id in rotor_ids:
            rows.extend(_rotor_rows(local_group_model, rotor_id))
    return tuple(sorted(set(int(row) for row in rows)))


def _active_atoms_for_term(local_group_model, rotor_ids, cluster=None):
    atoms = []
    if cluster is not None:
        atoms.extend(int(atom) for atom in getattr(cluster, "active_atoms", ()))
    if not atoms:
        for rotor_id in rotor_ids:
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is not None:
                atoms.extend(int(atom) for atom in rotor.owned_atoms)
    return tuple(sorted(set(int(atom) for atom in atoms)))


def _overlap_active_rows(local_group_model, rotor_ids):
    rotor_set = set(str(rotor_id) for rotor_id in rotor_ids)
    containing_rows = [
        set(int(row) for row in cluster.active_rows)
        for cluster in local_group_model.clusters.values()
        if rotor_set.issubset(set(str(rotor_id) for rotor_id in cluster.rotor_ids))
    ]

    rows = set.intersection(*containing_rows) if containing_rows else set()
    for rotor_id in rotor_ids:
        rows.update(_rotor_rows(local_group_model, rotor_id))
    return tuple(sorted(rows))


def _validate_local_group_row_partition(local_group_model):
    core_rows = set(int(row) for row in local_group_model.core_rows)
    active_rows = set()

    for cluster_id, cluster in local_group_model.clusters.items():
        cluster_rows = set(int(row) for row in cluster.active_rows)
        active_rows.update(cluster_rows)

        for rotor_id in cluster.rotor_ids:
            rotor_rows = set(_rotor_rows(local_group_model, rotor_id))
            missing = rotor_rows - cluster_rows
            if missing:
                raise RuntimeError(
                    f"Local cluster {cluster_id} does not own torsion rows "
                    f"{sorted(missing)} for rotor {rotor_id}."
                )

    overlap = core_rows & active_rows
    if overlap:
        raise RuntimeError(
            "Local-group core and factor rows overlap: "
            f"{sorted(overlap)}."
        )

    assigned_rows = core_rows | active_rows
    if assigned_rows:
        expected_rows = set(range(max(assigned_rows) + 1))
        missing = expected_rows - assigned_rows
        if missing:
            raise RuntimeError(
                "Local-group row partition is incomplete; unassigned rows: "
                f"{sorted(missing)}."
            )


def _records_for_cluster(point_labels_by_cluster, cluster_id):
    records = point_labels_by_cluster.get(str(cluster_id), {})
    if isinstance(records, dict):
        return [records[key] for key in sorted(records, key=lambda x: int(x))]
    return sorted(records, key=lambda rec: int(rec.get("state_id", 0)))


def _state_payload(state_id, rotor_ids, phase_signature, label, is_anchor=False):
    return {
        "state_id": int(state_id),
        "rotor_ids": [str(rotor_id) for rotor_id in rotor_ids],
        "phase_signature": _as_float_list(phase_signature),
        "is_anchor": bool(is_anchor),
        "label": str(label),
    }


def _find_source_factor_for_overlap(factors, overlap_key):
    overlap_set = set(str(rotor_id) for rotor_id in overlap_key)
    for cluster_id, cluster in factors:
        cluster_rotors = tuple(str(rotor_id) for rotor_id in cluster.rotor_ids)
        if overlap_set.issubset(set(cluster_rotors)):
            return str(cluster_id), cluster
    return None, None


def build_signed_factor_registry_payload(
        local_group_model,
        family_label,
        core_label,
        point_labels_by_cluster):
    """
    Build an inclusion-exclusion registry for possibly overlapping clusters.

    Each sampled coupled cluster is a positive factor. Intersections of rotor
    sets become signed overlap factors. Empty intersections become a signed
    coefficient on the core Taylor model. This preserves a single copy of the
    shared core and subtracts shared local rotors such as a tert-butyl parent
    rotor that appears in two direct factors.
    """

    _validate_local_group_row_partition(local_group_model)

    factors = sorted(
        ((str(cluster_id), cluster)
         for cluster_id, cluster in local_group_model.clusters.items()),
        key=lambda item: str(item[0]),
    )

    coefficient_by_key = Counter()
    core_coefficient = 0.0

    for order in range(1, len(factors) + 1):
        sign = 1.0 if order % 2 == 1 else -1.0
        for combo in combinations(factors, order):
            intersection = set(str(rotor_id) for rotor_id in combo[0][1].rotor_ids)
            for _, cluster in combo[1:]:
                intersection.intersection_update(
                    str(rotor_id) for rotor_id in cluster.rotor_ids
                )

            if intersection:
                coefficient_by_key[_rotor_key(intersection)] += sign
            else:
                core_coefficient += sign

    explicit_key_to_cluster = {
        _rotor_key(cluster.rotor_ids): (str(cluster_id), cluster)
        for cluster_id, cluster in factors
    }

    primitive_signature_rows = {}
    for _, cluster in factors:
        for rotor_id in cluster.rotor_ids:
            rotor_id = str(rotor_id)
            primitive_signature_rows.setdefault(
                rotor_id,
                tuple(int(row) for row in _rotor_rows(local_group_model, rotor_id)),
            )

    terms = {}
    term_state_labels = {}
    term_angle_library = {}
    overlap_counter = 0

    for rotor_key, coefficient in sorted(
            coefficient_by_key.items(), key=lambda item: item[0]):
        coefficient = float(coefficient)
        if abs(coefficient) <= 1.0e-12:
            continue

        explicit = explicit_key_to_cluster.get(rotor_key)
        if explicit is not None:
            cluster_id, cluster = explicit
            term_id = cluster_id
            term_rotor_ids = tuple(str(rotor_id) for rotor_id in cluster.rotor_ids)
            source_records = _records_for_cluster(point_labels_by_cluster, cluster_id)

            state_payloads = []
            state_labels = {}
            for record in source_records:
                state_id = int(record.get("state_id", len(state_payloads)))
                label = str(record["label"])
                phase = _phase_for_record(record, len(term_rotor_ids))
                state_payloads.append(
                    _state_payload(
                        state_id,
                        term_rotor_ids,
                        phase,
                        label,
                        is_anchor=bool(record.get("is_anchor", False)),
                    )
                )
                state_labels[int(state_id)] = label

            if not state_payloads:
                raise RuntimeError(
                    "No local-group datapoints were registered for sampled "
                    f"cluster {cluster_id}."
                )

            active_rows = _active_rows_for_term(
                local_group_model, term_rotor_ids, cluster=cluster)
            active_atoms = _active_atoms_for_term(
                local_group_model, term_rotor_ids, cluster=cluster)
            role = "factor"

        else:
            source_cluster_id, source_cluster = _find_source_factor_for_overlap(
                factors, rotor_key)
            if source_cluster_id is None:
                raise RuntimeError(
                    "Could not find a sampled source cluster for local-factor "
                    f"overlap {rotor_key}."
                )

            source_rotors = tuple(str(rotor_id) for rotor_id in source_cluster.rotor_ids)
            rotor_key_set = set(rotor_key)
            term_rotor_ids = tuple(
                rotor_id for rotor_id in source_rotors if rotor_id in rotor_key_set
            )
            keep_indices = [source_rotors.index(rotor_id) for rotor_id in term_rotor_ids]
            drop_indices = [
                idx for idx, rotor_id in enumerate(source_rotors)
                if rotor_id not in rotor_key_set
            ]

            state_payloads = []
            state_labels = {}
            seen_phases = set()

            for record in _records_for_cluster(
                    point_labels_by_cluster, source_cluster_id):
                source_phase = _phase_for_record(record, len(source_rotors))
                if any(abs(float(source_phase[idx])) > 1.0e-10
                       for idx in drop_indices):
                    continue

                phase = tuple(float(source_phase[idx]) for idx in keep_indices)
                phase_key = tuple(round(value, 12) for value in phase)
                if phase_key in seen_phases:
                    continue
                seen_phases.add(phase_key)

                state_id = len(state_payloads)
                label = str(record["label"])
                state_payloads.append(
                    _state_payload(
                        state_id,
                        term_rotor_ids,
                        phase,
                        label,
                        is_anchor=all(abs(value) <= 1.0e-12 for value in phase),
                    )
                )
                state_labels[state_id] = label

            if not state_payloads:
                raise RuntimeError(
                    "No compatible source states were found for local-factor "
                    f"overlap {rotor_key}."
                )

            term_id = f"overlap_{overlap_counter}"
            overlap_counter += 1
            active_rows = _overlap_active_rows(
                local_group_model, term_rotor_ids)
            active_atoms = _active_atoms_for_term(local_group_model, term_rotor_ids)
            role = "overlap"

        term_key = _json_key(term_id)
        terms[term_key] = {
            "term_id": str(term_id),
            "role": role,
            "coefficient": coefficient,
            "rotor_ids": [str(rotor_id) for rotor_id in term_rotor_ids],
            "active_rows": [int(row) for row in active_rows],
            "active_atoms": [int(atom) for atom in active_atoms],
            "grouped_signature_rows": [
                [int(row) for row in rows]
                for rows in _signature_rows(local_group_model, term_rotor_ids)
            ],
            "grouped_signature_row_types": [
                [str(row_type) for row_type in row_types]
                for row_types in _signature_row_types(
                    local_group_model, term_rotor_ids)
            ],
            "grouped_signature_row_scales": [
                [float(scale) for scale in row_scales]
                for row_scales in _signature_row_scales(
                    local_group_model, term_rotor_ids)
            ],
        }
        term_state_labels[term_key] = {
            str(int(state_id)): str(label)
            for state_id, label in state_labels.items()
        }
        term_angle_library[term_key] = state_payloads

    factor_info = {
        "schema_version": 2,
        "topology": "signed_local_factors",
        "family_label": str(family_label),
        "core_label": str(core_label),
        "core_rows": [int(row) for row in local_group_model.core_rows],
        "core_atoms": [int(atom) for atom in local_group_model.core_atoms],
        "core_coefficient": float(core_coefficient),
        "primitive_signature_rows": {
            _json_key(rotor_id): [int(row) for row in rows]
            for rotor_id, rows in primitive_signature_rows.items()
        },
        "primitive_signature_row_types": {
            _json_key(rotor_id): [
                str(row_type)
                for row_type in _rotor_row_types(local_group_model, rotor_id)
            ]
            for rotor_id in primitive_signature_rows
        },
        "primitive_signature_row_scales": {
            _json_key(rotor_id): [
                float(scale)
                for scale in _rotor_row_scales(local_group_model, rotor_id)
            ]
            for rotor_id in primitive_signature_rows
        },
        "terms": terms,
    }
    factor_index = {
        "family_label": str(family_label),
        "core_label": str(core_label),
        "term_state_labels": term_state_labels,
    }

    return factor_info, term_angle_library, factor_index


def write_signed_factor_registry_for_family(
        imff_file,
        root,
        family_label,
        local_group_model,
        core_label,
        point_labels_by_cluster):
    factor_info, factor_angle_library, factor_index = (
        build_signed_factor_registry_payload(
            local_group_model,
            family_label,
            core_label,
            point_labels_by_cluster,
        )
    )

    with h5py.File(imff_file, "a") as h5f:
        prefix = f"local_factor_registry/root_{root}/family_{family_label}"
        group = h5f.require_group(prefix)
        _write_string_dataset(
            group,
            "factor_info_json",
            json.dumps(factor_info, sort_keys=True),
        )
        _write_string_dataset(
            group,
            "factor_angle_library_json",
            json.dumps(factor_angle_library, sort_keys=True),
        )
        _write_string_dataset(
            group,
            "point_index_json",
            json.dumps(factor_index, sort_keys=True),
        )


def _iter_signed_factor_family_refs(imff_file, root=None):
    with h5py.File(imff_file, "r") as h5f:
        base = "local_factor_registry"
        if base not in h5f:
            return []

        if root is None:
            root_keys = sorted(h5f[base].keys())
        else:
            root_key = f"root_{root}"
            root_keys = [root_key] if root_key in h5f[base] else []

        refs = []
        for root_key in root_keys:
            root_group = h5f[base][root_key]
            for family_key in sorted(root_group.keys()):
                if family_key.startswith("family_"):
                    family = family_key.replace("family_", "", 1)
                else:
                    family = family_key
                refs.append((root_key.replace("root_", "", 1), family))
        return refs


def read_signed_factor_registry_for_family(imff_file, root, family_label):
    with h5py.File(imff_file, "r") as h5f:
        prefix = f"local_factor_registry/root_{root}/family_{family_label}"
        info_json = _read_scalar_string(h5f[prefix + "/factor_info_json"])
        library_json = _read_scalar_string(
            h5f[prefix + "/factor_angle_library_json"])
        index_json = _read_scalar_string(h5f[prefix + "/point_index_json"])

    return json.loads(info_json), json.loads(library_json), json.loads(index_json)


def _local_metadata_labels_from_hdf5(imff_file):
    labels = []
    seen = set()
    suffixes = ("_rinv_dihedral", "_eq_dihedral", "_r_dihedral")
    with h5py.File(imff_file, "r") as h5f:
        for key in h5f.keys():
            if not key.endswith("_bank_role"):
                continue

            label = key[:-len("_bank_role")]
            for suffix in suffixes:
                if label.endswith(suffix):
                    label = label[:-len(suffix)]
                    break

            if label not in seen:
                seen.add(label)
                labels.append(label)
    return labels


def _phase_for_datapoint(datapoint, n_rotors):
    phase = getattr(datapoint, "phase_signature", None)
    if phase is None:
        return np.zeros(int(n_rotors), dtype=np.float64)

    phase = np.asarray(phase, dtype=np.float64).reshape(-1)
    out = np.zeros(int(n_rotors), dtype=np.float64)
    ncopy = min(out.size, phase.size)
    if ncopy:
        out[:ncopy] = phase[:ncopy]
    return out


def _flat_z_matrix(z_matrix):
    if isinstance(z_matrix, dict):
        return InterpolationDatapoint.flatten_z_matrix(z_matrix)
    return list(z_matrix or [])


def _find_dihedral_row(flat_z_matrix, dihedral):
    dihedral = tuple(int(x) for x in np.asarray(dihedral).reshape(-1))
    reversed_dihedral = tuple(reversed(dihedral))
    for row_idx, coord in enumerate(flat_z_matrix):
        coord = tuple(int(x) for x in coord)
        if coord == dihedral or coord == reversed_dihedral:
            return int(row_idx)
    return None


def _axis_torsion_rows(flat_z_matrix, row_idx):
    coord = tuple(int(x) for x in flat_z_matrix[int(row_idx)])
    if len(coord) != 4:
        return (int(row_idx),)

    axis = set(coord[1:3])
    rows = []
    for idx, candidate in enumerate(flat_z_matrix):
        candidate = tuple(int(x) for x in candidate)
        if len(candidate) == 4 and set(candidate[1:3]) == axis:
            rows.append(int(idx))
    return tuple(rows) if rows else (int(row_idx),)


def _dihedrals_to_rotate(datapoint):
    dihedrals = getattr(datapoint, "dihedrals_to_rotate", None)
    if dihedrals is None:
        return []
    arr = np.asarray(dihedrals, dtype=np.int64)
    if arr.size == 0:
        return []
    return [tuple(int(x) for x in row) for row in arr.reshape(-1, 4)]


def _infer_rotor_signature_rows(cluster_infos, flat_z_matrix):
    rows_by_rotor = {}

    for cluster_info in cluster_infos.values():
        rotor_ids = tuple(str(rotor_id) for rotor_id in cluster_info["rotor_ids"])
        for rotor_id in rotor_ids:
            rows_by_rotor.setdefault(rotor_id, set())

        for datapoint in cluster_info["states"].values():
            phase = _phase_for_datapoint(datapoint, len(rotor_ids))
            nonzero_indices = [
                idx for idx, value in enumerate(phase)
                if abs(float(value)) > 1.0e-10
            ]
            dihedrals = _dihedrals_to_rotate(datapoint)
            for rotor_index, dihedral in zip(nonzero_indices, dihedrals):
                row_idx = _find_dihedral_row(flat_z_matrix, dihedral)
                if row_idx is None:
                    continue
                rows_by_rotor[rotor_ids[rotor_index]].update(
                    _axis_torsion_rows(flat_z_matrix, row_idx)
                )

    for cluster_info in cluster_infos.values():
        active_rows = tuple(int(row) for row in cluster_info.get("active_rows", ()))
        for rotor_id in cluster_info["rotor_ids"]:
            rotor_id = str(rotor_id)
            if not rows_by_rotor.get(rotor_id):
                rows_by_rotor.setdefault(rotor_id, set()).update(active_rows)

    return {
        rotor_id: tuple(sorted(int(row) for row in rows))
        for rotor_id, rows in rows_by_rotor.items()
    }


def _build_inferred_factor_bank(family_label, core_dp, cluster_infos, flat_z_matrix):
    factors = sorted(cluster_infos.items(), key=lambda item: str(item[0]))
    coefficient_by_key = Counter()
    core_coefficient = 0.0

    for order in range(1, len(factors) + 1):
        sign = 1.0 if order % 2 == 1 else -1.0
        for combo in combinations(factors, order):
            intersection = set(str(rotor_id) for rotor_id in combo[0][1]["rotor_ids"])
            for _, cluster_info in combo[1:]:
                intersection.intersection_update(
                    str(rotor_id) for rotor_id in cluster_info["rotor_ids"]
                )
            if intersection:
                coefficient_by_key[_rotor_key(intersection)] += sign
            else:
                core_coefficient += sign

    explicit_key_to_cluster = {
        _rotor_key(info["rotor_ids"]): (str(cluster_id), info)
        for cluster_id, info in factors
    }
    rows_by_rotor = _infer_rotor_signature_rows(cluster_infos, flat_z_matrix)

    terms = {}
    overlap_counter = 0
    for rotor_key, coefficient in sorted(
            coefficient_by_key.items(), key=lambda item: item[0]):
        coefficient = float(coefficient)
        if abs(coefficient) <= 1.0e-12:
            continue

        explicit = explicit_key_to_cluster.get(rotor_key)
        if explicit is not None:
            term_id, source_info = explicit
            term_rotor_ids = tuple(str(rotor_id) for rotor_id in source_info["rotor_ids"])
            expected_states = dict(source_info["states"])
            active_rows = set(
                int(row) for row in source_info.get("active_rows", ()))
            for rotor_id in term_rotor_ids:
                active_rows.update(rows_by_rotor.get(str(rotor_id), ()))
            active_rows = tuple(sorted(active_rows))
            active_atoms = tuple(int(atom) for atom in source_info.get("active_atoms", ()))
            role = "factor"
        else:
            source_cluster_id, source_info = None, None
            overlap_set = set(rotor_key)
            for candidate_id, candidate_info in factors:
                candidate_rotors = tuple(str(rotor_id) for rotor_id in candidate_info["rotor_ids"])
                if overlap_set.issubset(set(candidate_rotors)):
                    source_cluster_id = str(candidate_id)
                    source_info = candidate_info
                    break
            if source_info is None:
                continue

            source_rotors = tuple(str(rotor_id) for rotor_id in source_info["rotor_ids"])
            term_rotor_ids = tuple(
                rotor_id for rotor_id in source_rotors if rotor_id in overlap_set
            )
            keep_indices = [source_rotors.index(rotor_id) for rotor_id in term_rotor_ids]
            drop_indices = [
                idx for idx, rotor_id in enumerate(source_rotors)
                if rotor_id not in overlap_set
            ]

            expected_states = {}
            seen_phases = set()
            for datapoint in source_info["states"].values():
                source_phase = _phase_for_datapoint(datapoint, len(source_rotors))
                if any(abs(float(source_phase[idx])) > 1.0e-10
                       for idx in drop_indices):
                    continue
                phase = tuple(float(source_phase[idx]) for idx in keep_indices)
                phase_key = tuple(round(value, 12) for value in phase)
                if phase_key in seen_phases:
                    continue
                seen_phases.add(phase_key)
                expected_states[len(expected_states)] = datapoint

            if not expected_states:
                continue

            term_id = f"overlap_{overlap_counter}"
            overlap_counter += 1
            active_rows = tuple(
                sorted({
                    int(row)
                    for rotor_id in term_rotor_ids
                    for row in rows_by_rotor.get(str(rotor_id), ())
                })
            )
            active_atoms = ()
            role = "overlap"

        terms[str(term_id)] = {
            "term_id": str(term_id),
            "role": role,
            "coefficient": coefficient,
            "rotor_ids": term_rotor_ids,
            "active_rows": active_rows,
            "active_atoms": active_atoms,
            "grouped_signature_rows": tuple(
                tuple(int(row) for row in rows_by_rotor.get(str(rotor_id), ()))
                for rotor_id in term_rotor_ids
            ),
            "grouped_signature_row_types": tuple(
                tuple("torsion" for _ in rows_by_rotor.get(str(rotor_id), ()))
                for rotor_id in term_rotor_ids
            ),
            "grouped_signature_row_scales": tuple(
                tuple(1.0 for _ in rows_by_rotor.get(str(rotor_id), ()))
                for rotor_id in term_rotor_ids
            ),
            "expected_states": expected_states,
            "angle_library": [],
        }

    return {
        "topology": "signed_local_factors",
        "root": None,
        "factor_info": {"schema_version": 0, "source": "inferred_metadata"},
        "factor_angle_library": {},
        "point_index": {"family_label": str(family_label)},
        "core": core_dp,
        "core_coefficient": float(core_coefficient),
        "core_rows": _as_int_tuple(getattr(core_dp, "active_rows", None)),
        "primitive_signature_rows": rows_by_rotor,
        "terms": terms,
    }


def infer_signed_factor_banks_from_datapoint_metadata(imff_file, z_matrix, im_settings):
    labels = _local_metadata_labels_from_hdf5(imff_file)
    if not labels:
        return {}

    flat_z_matrix = _flat_z_matrix(z_matrix)
    families = {}

    for label in labels:
        datapoint = InterpolationDatapoint(z_matrix)
        datapoint.update_settings(im_settings)
        datapoint.read_hdf5(imff_file, label)

        family_label = getattr(datapoint, "family_label", None)
        if family_label is None:
            continue

        family = families.setdefault(
            str(family_label), {"core": None, "clusters": {}})
        bank_role = getattr(datapoint, "bank_role", "global")
        if bank_role == "core":
            family["core"] = datapoint
        elif bank_role == "cluster":
            cluster_id = str(getattr(datapoint, "cluster_id", ""))
            if not cluster_id:
                continue
            cluster_info = family["clusters"].setdefault(
                cluster_id,
                {
                    "rotor_ids": tuple(str(rotor_id) for rotor_id in (
                        getattr(datapoint, "cluster_rotor_ids", ()) or ())),
                    "active_rows": _as_int_tuple(getattr(datapoint, "active_rows", None)),
                    "active_atoms": _as_int_tuple(getattr(datapoint, "active_atoms", None)),
                    "states": {},
                },
            )
            state_id = int(getattr(datapoint, "cluster_state_id", 0))
            cluster_info["states"][state_id] = datapoint

    banks = {}
    for family_label, family in families.items():
        core_dp = family.get("core")
        clusters = family.get("clusters", {})
        if core_dp is None or not clusters:
            continue
        banks[family_label] = _build_inferred_factor_bank(
            family_label,
            core_dp,
            clusters,
            flat_z_matrix,
        )

    return banks


def load_signed_factor_banks_for_root(imff_file, root, z_matrix, im_settings):
    banks = {}
    if imff_file is None:
        return banks

    registry_refs = _iter_signed_factor_family_refs(imff_file, root)
    if not registry_refs:
        return infer_signed_factor_banks_from_datapoint_metadata(
            imff_file, z_matrix, im_settings)

    for registry_root, family in registry_refs:
        factor_info, factor_angle_library, factor_index = (
            read_signed_factor_registry_for_family(
                imff_file, registry_root, family)
        )

        core_dp = InterpolationDatapoint(z_matrix)
        core_dp.update_settings(im_settings)
        core_dp.read_hdf5(imff_file, factor_index["core_label"])

        terms = {}
        for term_key, term_payload in factor_info.get("terms", {}).items():
            state_labels = factor_index.get("term_state_labels", {}).get(
                term_key, {})
            if not state_labels:
                state_labels = {
                    str(state["state_id"]): state["label"]
                    for state in factor_angle_library.get(term_key, [])
                }

            expected_states = {}
            for state_id, label in state_labels.items():
                dp = InterpolationDatapoint(z_matrix)
                dp.update_settings(im_settings)
                dp.read_hdf5(imff_file, label)
                expected_states[int(state_id)] = dp

            grouped_signature_rows = tuple(
                tuple(int(row) for row in rows)
                for rows in term_payload.get("grouped_signature_rows", ())
            )
            grouped_signature_row_types = tuple(
                tuple(str(row_type).strip().lower() for row_type in row_types)
                for row_types in term_payload.get(
                    "grouped_signature_row_types", ())
            )
            if len(grouped_signature_row_types) != len(grouped_signature_rows):
                grouped_signature_row_types = tuple(
                    tuple("torsion" for _ in rows)
                    for rows in grouped_signature_rows
                )

            grouped_signature_row_scales = tuple(
                tuple(float(scale) for scale in row_scales)
                for row_scales in term_payload.get(
                    "grouped_signature_row_scales", ())
            )
            if len(grouped_signature_row_scales) != len(grouped_signature_rows):
                grouped_signature_row_scales = tuple(
                    tuple(
                        0.35 if row_type == "angle" else 1.0
                        for row_type in row_types
                    )
                    for row_types in grouped_signature_row_types
                )
            active_rows = set(
                int(row) for row in term_payload.get("active_rows", ())
            )
            for rows in grouped_signature_rows:
                active_rows.update(rows)

            terms[str(term_payload.get("term_id", term_key))] = {
                "term_id": str(term_payload.get("term_id", term_key)),
                "role": str(term_payload.get("role", "factor")),
                "coefficient": float(term_payload.get("coefficient", 1.0)),
                "rotor_ids": tuple(
                    _restore_id(rotor_id)
                    for rotor_id in term_payload.get("rotor_ids", ())
                ),
                "active_rows": tuple(sorted(active_rows)),
                "active_atoms": tuple(
                    int(atom) for atom in term_payload.get("active_atoms", ())
                ),
                "grouped_signature_rows": grouped_signature_rows,
                "grouped_signature_row_types": grouped_signature_row_types,
                "grouped_signature_row_scales": grouped_signature_row_scales,
                "expected_states": expected_states,
                "angle_library": factor_angle_library.get(term_key, []),
            }

        banks[str(family)] = {
            "topology": "signed_local_factors",
            "root": str(registry_root),
            "factor_info": factor_info,
            "factor_angle_library": factor_angle_library,
            "point_index": factor_index,
            "core": core_dp,
            "core_coefficient": float(factor_info.get("core_coefficient", 0.0)),
            "core_rows": tuple(int(row) for row in factor_info.get("core_rows", ())),
            "primitive_signature_rows": {
                _restore_id(rotor_id): tuple(int(row) for row in rows)
                for rotor_id, rows in factor_info.get(
                    "primitive_signature_rows", {}).items()
            },
            "primitive_signature_row_types": {
                _restore_id(rotor_id): tuple(
                    str(row_type).strip().lower() for row_type in row_types
                )
                for rotor_id, row_types in factor_info.get(
                    "primitive_signature_row_types", {}).items()
            },
            "primitive_signature_row_scales": {
                _restore_id(rotor_id): tuple(float(scale) for scale in row_scales)
                for rotor_id, row_scales in factor_info.get(
                    "primitive_signature_row_scales", {}).items()
            },
            "terms": terms,
        }

    return banks


def iter_signed_factor_datapoints(family_bank):
    core = family_bank.get("core")
    if core is not None:
        yield core
    for term_bank in family_bank.get("terms", {}).values():
        for dp in term_bank.get("expected_states", {}).values():
            if dp is not None:
                yield dp


def _coerce_mapping_masks(datapoint, masks=None):
    if masks is None:
        masks = getattr(datapoint, "mapping_masks", None)

    n_ic = len(datapoint.internal_coordinates_values)
    if masks is None or len(masks) == 0:
        return np.arange(n_ic, dtype=np.int64).reshape(1, -1)

    masks_arr = np.asarray(masks, dtype=np.int64)
    if masks_arr.ndim == 1:
        masks_arr = masks_arr.reshape(1, -1)
    return masks_arr


def _normalize_symmetry_candidate_weight(
        driver,
        d2,
        grad_d2,
        confidence_radius=1.0,
        power=2,
        eps=1.0e-12):
    rho2 = max(float(confidence_radius)**2, eps)
    x = float(d2) / rho2 + eps
    w = 1.0 / (2.0 * (x**power))
    dw_dd2 = -(power / (2.0 * rho2)) * (x**(-(power + 1)))
    grad_w = dw_dd2 * np.asarray(grad_d2, dtype=np.float64)

    natm = driver.impes_coordinate.cartesian_coordinates.shape[0]
    return float(w), grad_w.reshape(natm, 3)


def _normalize_mask_weight_pool(driver, raw, raw_grad, eps=1.0e-14):
    raw = np.asarray(raw, dtype=np.float64)
    raw_grad = np.asarray(raw_grad, dtype=np.float64)

    natm = driver.impes_coordinate.cartesian_coordinates.shape[0]
    if raw.size == 0:
        return (
            np.zeros((0,), dtype=np.float64),
            np.zeros((0, natm, 3), dtype=np.float64),
        )

    total = float(raw.sum())
    if total < eps:
        idx = int(np.argmax(raw))
        weights = np.zeros_like(raw)
        weights[idx] = 1.0
        grad_weights = np.zeros_like(raw_grad)
        return weights, grad_weights

    sum_grad = raw_grad.sum(axis=0)
    weights = raw / total
    grad_weights = (raw_grad * total - raw[:, None, None] * sum_grad) / (total**2)
    return weights, grad_weights


def _candidate_current_chart(
        driver,
        org_int_coords,
        b_matrix,
        datapoint=None,
        mask=None):
    org_values = np.asarray(org_int_coords, dtype=np.float64)
    b_values = np.asarray(b_matrix, dtype=np.float64)

    if (
        datapoint is None
        or mask is None
        or not driver.impes_coordinate.use_eq_bond_length
        or driver._get_eq_symmetry_mode() == 'symmetrized'
    ):
        return org_values, b_values

    eq_reference = getattr(datapoint, 'eq_bond_lengths', None)
    if eq_reference is None:
        raise RuntimeError(
            "Local-factor datapoint is missing eq_bond_lengths: "
            f"{getattr(datapoint, 'point_label', None)}")

    base_cache = getattr(driver, '_eq_candidate_base_cache', None)
    if base_cache is None:
        base_cache = driver.impes_coordinate.prepare_eq_candidate_base_cache()
        driver._eq_candidate_base_cache = base_cache

    return (
        driver.impes_coordinate
        .build_masked_current_chart_from_reference_eq_fast(
            reference_eq_bond_lengths=eq_reference,
            mask=mask,
            base_cache=base_cache,
        )
    )


def _append_local_rmsd_diagnostics(driver, delta_eff, bounds):
    if getattr(driver, "bond_rmsd", None) is not None:
        bend = int(bounds["bond_end"])
        value = 0.0 if bend == 0 else np.sqrt(np.mean(delta_eff[:bend]**2))
        driver.bond_rmsd.append(float(value))

    if getattr(driver, "angle_rmsd", None) is not None:
        bend = int(bounds["bond_end"])
        aend = int(bounds["angle_end"])
        value = 0.0 if aend <= bend else np.sqrt(np.mean(delta_eff[bend:aend]**2))
        driver.angle_rmsd.append(float(value))

    if getattr(driver, "dihedral_rmsd", None) is not None:
        dstart = int(bounds["dihedral_start"])
        dend = int(bounds["dihedral_end"])
        value = 0.0 if dend <= dstart else np.sqrt(
            np.mean(delta_eff[dstart:dend]**2))
        driver.dihedral_rmsd.append(float(value))


def _restricted_taylor_eval(
        driver,
        *,
        energy,
        grad,
        hess,
        ref_coords,
        org_int_coords,
        b_matrix,
        projector):
    bounds = driver._get_internal_coordinate_partitions()
    dihedral_start = int(bounds["dihedral_start"])
    dihedral_end = int(bounds["dihedral_end"])

    delta_raw = (
        np.asarray(org_int_coords, dtype=np.float64)
        - np.asarray(ref_coords, dtype=np.float64)
    )
    delta_eff = delta_raw.copy()
    chain = np.ones_like(delta_raw)

    d_prop = driver._principal_torsion_delta(
        delta_raw[dihedral_start:dihedral_end])
    delta_eff[dihedral_start:dihedral_end] = 2.0 * np.sin(0.5 * d_prop)
    chain[dihedral_start:dihedral_end] = np.cos(0.5 * d_prop)

    d_imp = driver._principal_torsion_delta(delta_raw[dihedral_end:])
    imp_slice = slice(dihedral_end, len(delta_raw))
    delta_eff[imp_slice] = 2.0 * np.tan(0.5 * d_imp)
    chain[imp_slice] = 1.0 / np.maximum(
        np.cos(0.5 * d_imp)**2,
        1.0e-12,
    )

    _append_local_rmsd_diagnostics(driver, delta_eff, bounds)

    projector = np.asarray(projector, dtype=np.float64)
    delta_projected = projector * delta_eff

    energy = float(energy)
    grad = np.asarray(grad, dtype=np.float64)
    hess = np.asarray(hess, dtype=np.float64)

    value = float(
        energy
        + np.matmul(delta_projected.T, grad)
        + 0.5 * np.linalg.multi_dot([delta_projected.T, hess, delta_projected])
    )
    grad_internal_eff = projector * (grad + np.matmul(hess, delta_projected))
    grad_internal = chain * grad_internal_eff

    b_matrix = np.asarray(b_matrix, dtype=np.float64)
    natm = b_matrix.shape[1] // 3
    grad_cart = np.matmul(b_matrix.T, grad_internal).reshape(natm, 3)
    return value, grad_cart


def _restricted_taylor_eval_for_mask(
        driver,
        *,
        datapoint,
        mask,
        projector,
        org_int_coords,
        precomputed_current=None):
    mask = np.asarray(mask, dtype=np.int64)
    n_ic = len(datapoint.internal_coordinates_values)
    mask0 = np.arange(n_ic, dtype=np.int64)

    grad_eff = np.asarray(datapoint.internal_gradient, dtype=np.float64).copy()
    hess_eff = np.asarray(datapoint.internal_hessian, dtype=np.float64).copy()
    ref_masked = np.asarray(
        datapoint.internal_coordinates_values, dtype=np.float64)[mask]

    grad_eff[mask0] = grad_eff[mask]
    hess_eff[np.ix_(mask0, mask0)] = hess_eff[np.ix_(mask, mask)]

    if precomputed_current is None:
        org_eval, b_eval = _candidate_current_chart(
            driver,
            org_int_coords,
            driver.impes_coordinate.b_matrix,
            datapoint=datapoint,
            mask=mask,
        )
    else:
        org_eval, b_eval = precomputed_current

    return _restricted_taylor_eval(
        driver,
        energy=float(datapoint.energy),
        grad=grad_eff,
        hess=hess_eff,
        ref_coords=ref_masked,
        org_int_coords=org_eval,
        b_matrix=b_eval,
        projector=np.asarray(projector, dtype=np.float64),
    )


def _rows_are_nonempty(rows):
    if rows is None:
        return False
    return np.asarray(rows).size > 0


def _factor_core_rows(family_bank, n_ic):
    occupied = set()
    for term_bank in family_bank.get("terms", {}).values():
        occupied.update(int(row) for row in term_bank.get("active_rows", ()))

    core_dp = family_bank.get("core")
    stored_rows = family_bank.get("core_rows", ())
    if not stored_rows and core_dp is not None and _rows_are_nonempty(
            getattr(core_dp, "active_rows", None)):
        stored_rows = tuple(
            int(row) for row in np.asarray(core_dp.active_rows).reshape(-1))

    schema_version = int(
        family_bank.get("factor_info", {}).get("schema_version", 0))
    if schema_version >= 2:
        rows = set(int(row) for row in stored_rows)
        overlap = rows & occupied
        missing = set(range(int(n_ic))) - rows - occupied
        if overlap or missing:
            raise RuntimeError(
                "Invalid schema-2 local-factor row partition: "
                f"core/factor overlap={sorted(overlap)}, "
                f"unassigned={sorted(missing)}."
            )
        return tuple(sorted(rows))

    # Schema 0/1 compatibility: old registries did not classify boundary rows.
    # Signature rows have already been restored to their factors while genuinely
    # unassigned rows retain the historical core fallback.
    return tuple(sorted(set(range(int(n_ic))) - occupied))


def _factor_core_projector(family_bank, n_ic):
    projector = np.zeros(n_ic, dtype=np.float64)
    rows = _factor_core_rows(family_bank, n_ic)
    if rows:
        projector[list(rows)] = 1.0
    return projector


def _factor_term_projector(family_bank, term_bank, n_ic):
    rows = set(_factor_core_rows(family_bank, n_ic))
    rows.update(int(row) for row in term_bank.get("active_rows", ()))
    projector = np.zeros(n_ic, dtype=np.float64)
    if rows:
        projector[list(sorted(rows))] = 1.0
    return projector


def _build_factor_signature(term_bank, int_coords):
    int_coords = np.asarray(int_coords, dtype=np.float64)
    return [
        int_coords[list(rows)].copy()
        for rows in term_bank.get("grouped_signature_rows", ())
    ]


def _factor_state_weight_mode(driver):
    mode = getattr(driver, "local_factor_state_weight_mode", "signature")
    mode = "signature" if mode is None else str(mode).strip().lower()
    mode = mode.replace("-", "_")

    aliases = {
        "": "signature",
        "signature": "signature",
        "torsion": "signature",
        "torsion_signature": "signature",
        "rotor_signature": "signature",
        "active": "active_internal",
        "active_rows": "active_internal",
        "active_internal": "active_internal",
        "active_atoms": "active_atoms_internal",
        "active_atom_rows": "active_atoms_internal",
        "active_atoms_internal": "active_atoms_internal",
        "local_internal": "active_atoms_internal",
        "signature_active": "signature_active_internal",
        "signature_active_internal": "signature_active_internal",
        "hybrid_active": "signature_active_internal",
        "hybrid_active_internal": "signature_active_internal",
        "signature_active_atoms": "signature_active_atoms_internal",
        "signature_active_atom_rows": "signature_active_atoms_internal",
        "signature_active_atoms_internal": "signature_active_atoms_internal",
        "hybrid_active_atoms": "signature_active_atoms_internal",
        "hybrid_active_atoms_internal": "signature_active_atoms_internal",
        "all": "all_internal",
        "all_rows": "all_internal",
        "all_internal": "all_internal",
        "cartesian": "cartesian_active",
        "cartesian_active": "cartesian_active",
        "active_cartesian": "cartesian_active",
    }
    if mode not in aliases:
        raise ValueError(
            "Unsupported local_factor_state_weight_mode="
            f"'{getattr(driver, 'local_factor_state_weight_mode', mode)}'."
        )
    return aliases[mode]


def _factor_signature_row_set(term_bank):
    rows = set()
    for group_rows in term_bank.get("grouped_signature_rows", ()):
        rows.update(int(row) for row in group_rows)
    return rows


def _factor_internal_metric_rows(driver, term_bank, mode, n_ic):
    signature_rows = _factor_signature_row_set(term_bank)

    if mode == "active_internal":
        rows = set(int(row) for row in term_bank.get("active_rows", ()))
        rows.update(signature_rows)

    elif mode == "active_atoms_internal":
        z_matrix = getattr(driver.impes_coordinate, "z_matrix", None)
        active_atoms = set(int(atom) for atom in term_bank.get("active_atoms", ()))
        if z_matrix is None or not active_atoms:
            rows = set(int(row) for row in term_bank.get("active_rows", ()))
        else:
            rows = {
                idx for idx, coord in enumerate(z_matrix[:int(n_ic)])
                if set(int(atom) for atom in coord) & active_atoms
            }
        rows.update(signature_rows)

    elif mode == "all_internal":
        rows = set(range(int(n_ic)))

    else:
        rows = set()

    return tuple(sorted(row for row in rows if 0 <= int(row) < int(n_ic)))


def _factor_bond_angle_rows(driver, rows):
    selected = []
    for row in rows:
        section = _factor_row_section(driver, row)
        if section in ("bonds", "angles"):
            selected.append(int(row))
    return tuple(selected)


def _factor_row_section(driver, row_idx):
    z_matrix = getattr(driver.impes_coordinate, "z_matrix", None)
    if z_matrix is None:
        return None

    coord = tuple(int(atom) for atom in z_matrix[int(row_idx)])
    if len(coord) == 2:
        return "bonds"
    if len(coord) == 3:
        return "angles"

    bounds = driver._get_internal_coordinate_partitions()
    if int(row_idx) >= int(bounds["dihedral_end"]):
        return "impropers"
    return "dihedrals"


def _factor_internal_rows_distance_and_gradient(
        driver,
        rows,
        org_eval,
        ref_eval,
        b_matrix):
    rows = tuple(int(row) for row in rows)
    b_matrix = np.asarray(b_matrix, dtype=np.float64)
    ncart = b_matrix.shape[1]
    if not rows:
        return 0.0, np.zeros(ncart, dtype=np.float64)

    q_cur = np.asarray(org_eval, dtype=np.float64)
    q_ref = np.asarray(ref_eval, dtype=np.float64)

    d2_total = 0.0
    grad_total = np.zeros(ncart, dtype=np.float64)

    for row in rows:
        section = _factor_row_section(driver, row)
        if section is None:
            continue

        sigma = driver._imp_coordinate_sigma(section, row, q_ref)
        sigma2 = max(float(sigma) * float(sigma), 1.0e-24)
        delta_raw = float(q_cur[row] - q_ref[row])

        if section == "dihedrals":
            delta = driver._principal_torsion_delta(delta_raw)
            metric = 2.0 * (1.0 - np.cos(delta)) / sigma2
            grad_coeff = 2.0 * np.sin(delta) / sigma2

        elif section == "impropers":
            delta = driver._principal_torsion_delta(delta_raw)
            eta = driver._improper_dihedral_displacement(delta)
            chain = driver._improper_dihedral_chain(delta)
            metric = (eta * eta) / sigma2
            grad_coeff = 2.0 * eta * chain / sigma2

        else:
            metric = (delta_raw * delta_raw) / sigma2
            grad_coeff = 2.0 * delta_raw / sigma2

        d2_total += float(metric)
        grad_total += grad_coeff * b_matrix[row, :]

    count = max(len(rows), 1)
    return float(d2_total / count), grad_total / count


def _factor_cartesian_active_distance_and_gradient(driver, term_bank, datapoint):
    active_atoms = tuple(int(atom) for atom in term_bank.get("active_atoms", ()))
    natm = driver.impes_coordinate.cartesian_coordinates.shape[0]
    if not active_atoms:
        active_atoms = tuple(range(natm))

    active_atoms_arr = np.asarray(active_atoms, dtype=np.int64)
    current = np.asarray(
        driver.impes_coordinate.cartesian_coordinates,
        dtype=np.float64,
    )[active_atoms_arr]
    reference = np.asarray(
        datapoint.cartesian_coordinates,
        dtype=np.float64,
    )[active_atoms_arr]

    current_centered = current - np.mean(current, axis=0)
    reference_centered = reference - np.mean(reference, axis=0)
    delta = current_centered - reference_centered

    d2 = float(np.dot(delta.reshape(-1), delta.reshape(-1)))
    grad_sub = 2.0 * (delta - np.mean(delta, axis=0))

    inv_sqrt = getattr(driver.impes_coordinate, "inv_sqrt_masses", None)
    if inv_sqrt is not None:
        inv_sqrt = np.asarray(inv_sqrt, dtype=np.float64).reshape(-1)
        active_dofs = np.array(
            [3 * atom + comp for atom in active_atoms_arr for comp in range(3)],
            dtype=np.int64,
        )
        grad_sub = grad_sub.reshape(-1) * inv_sqrt[active_dofs]
        grad_sub = grad_sub.reshape(active_atoms_arr.size, 3)

    grad = np.zeros((natm, 3), dtype=np.float64)
    grad[active_atoms_arr] = grad_sub

    return d2, grad.reshape(-1)


def _factor_signature_distance_and_gradient(
        term_bank,
        current_signature,
        reference_signature,
        b_matrix):
    rows_by_rotor = term_bank.get("grouped_signature_rows", ())
    row_types_by_rotor = term_bank.get("grouped_signature_row_types", ())
    row_scales_by_rotor = term_bank.get("grouped_signature_row_scales", ())
    b_matrix = np.asarray(b_matrix, dtype=np.float64)

    ncart = b_matrix.shape[1]
    d2_total = 0.0
    grad_total = np.zeros(ncart, dtype=np.float64)
    n_groups = len(rows_by_rotor)

    for rotor_index, rows in enumerate(rows_by_rotor):
        if len(rows) == 0:
            continue
        rows_arr = np.asarray(rows, dtype=np.int64)
        delta = (
            np.asarray(current_signature[rotor_index], dtype=np.float64)
            - np.asarray(reference_signature[rotor_index], dtype=np.float64)
        )

        row_types = ()
        if rotor_index < len(row_types_by_rotor):
            row_types = tuple(
                str(row_type).strip().lower()
                for row_type in row_types_by_rotor[rotor_index]
            )
        if len(row_types) != len(rows):
            row_types = tuple("torsion" for _ in rows)

        row_scales = ()
        if rotor_index < len(row_scales_by_rotor):
            row_scales = tuple(
                max(float(scale), 1.0e-12)
                for scale in row_scales_by_rotor[rotor_index]
            )
        if len(row_scales) != len(rows):
            row_scales = tuple(
                0.35 if row_type == "angle" else 1.0
                for row_type in row_types
            )

        row_metric = np.zeros(len(rows), dtype=np.float64)
        row_grad_coeff = np.zeros(len(rows), dtype=np.float64)
        for idx, (row_type, scale) in enumerate(zip(row_types, row_scales)):
            if row_type in ("torsion", "dihedral", "periodic"):
                row_metric[idx] = 2.0 * (1.0 - np.cos(delta[idx])) / (scale**2)
                row_grad_coeff[idx] = 2.0 * np.sin(delta[idx]) / (scale**2)
            else:
                row_metric[idx] = (delta[idx] / scale)**2
                row_grad_coeff[idx] = 2.0 * delta[idx] / (scale**2)

        d2_total += float(np.mean(row_metric))
        grad_total += np.mean(
            row_grad_coeff[:, None] * b_matrix[rows_arr, :],
            axis=0,
        )

    if n_groups > 0:
        d2_total /= n_groups
        grad_total /= n_groups

    return d2_total, grad_total


def _factor_mask_signature_distance_and_gradient(
        driver,
        *,
        term_bank,
        datapoint,
        mask,
        org_int_coords,
        b_matrix=None,
        precomputed_current=None,
        mode=None):
    mask = np.asarray(mask, dtype=np.int64)
    mode = "signature" if mode is None else str(mode)

    if mode == "cartesian_active":
        return _factor_cartesian_active_distance_and_gradient(
            driver, term_bank, datapoint)

    if precomputed_current is None:
        org_eval, b_eval = _candidate_current_chart(
            driver,
            org_int_coords,
            driver.impes_coordinate.b_matrix if b_matrix is None else b_matrix,
            datapoint=datapoint,
            mask=mask,
        )
    else:
        org_eval, b_eval = precomputed_current

    cur_sig = _build_factor_signature(term_bank, org_eval)
    ref_masked = np.asarray(
        datapoint.internal_coordinates_values, dtype=np.float64)[mask]

    if mode == "signature":
        ref_sig = _build_factor_signature(term_bank, ref_masked)
        return _factor_signature_distance_and_gradient(
            term_bank,
            cur_sig,
            ref_sig,
            b_matrix=b_eval,
        )

    if mode in ("signature_active_internal", "signature_active_atoms_internal"):
        ref_sig = _build_factor_signature(term_bank, ref_masked)
        sig_d2, sig_grad = _factor_signature_distance_and_gradient(
            term_bank,
            cur_sig,
            ref_sig,
            b_matrix=b_eval,
        )
        internal_mode = {
            "signature_active_internal": "active_internal",
            "signature_active_atoms_internal": "active_atoms_internal",
        }[mode]
        rows = _factor_bond_angle_rows(
            driver,
            _factor_internal_metric_rows(
                driver,
                term_bank,
                internal_mode,
                len(ref_masked),
            ),
        )
        row_d2, row_grad = _factor_internal_rows_distance_and_gradient(
            driver,
            rows,
            org_eval,
            ref_masked,
            b_eval,
        )
        return sig_d2 + row_d2, sig_grad + row_grad

    rows = _factor_internal_metric_rows(
        driver,
        term_bank,
        mode,
        len(ref_masked),
    )
    return _factor_internal_rows_distance_and_gradient(
        driver,
        rows,
        org_eval,
        ref_masked,
        b_eval,
    )


def _factor_state_signature_metric(driver, term_bank, datapoint, org_int_coords):
    masks = _coerce_mapping_masks(datapoint)
    mode = _factor_state_weight_mode(driver)

    best_d2 = None
    best_grad = None
    for mask in masks:
        d2, grad_d2 = _factor_mask_signature_distance_and_gradient(
            driver,
            term_bank=term_bank,
            datapoint=datapoint,
            mask=mask,
            org_int_coords=org_int_coords,
            mode=mode,
        )
        if best_d2 is None or d2 < best_d2:
            best_d2, best_grad = d2, grad_d2

    return best_d2, best_grad


def _compute_factor_weights(driver, term_bank, org_int_coords):
    populated = [
        (state_id, dp)
        for state_id, dp in sorted(term_bank["expected_states"].items())
        if dp is not None
    ]

    natm = driver.impes_coordinate.cartesian_coordinates.shape[0]
    if not populated:
        return (
            np.zeros((0,), dtype=np.float64),
            np.zeros((0, natm, 3), dtype=np.float64),
            [],
        )

    metric_results = [
        _factor_state_signature_metric(driver, term_bank, dp, org_int_coords)
        for _, dp in populated
    ]
    d2_values = np.asarray([item[0] for item in metric_results], dtype=np.float64)
    exact_idx = np.where(d2_values <= 1.0e-10)[0]

    if exact_idx.size > 0:
        weights = np.zeros(len(populated), dtype=np.float64)
        weights[exact_idx] = 1.0 / exact_idx.size
        grad_weights = np.zeros((len(populated), natm, 3), dtype=np.float64)
        return weights, grad_weights, populated

    raw = []
    raw_grad = []
    p = float(driver.exponent_p if driver.exponent_p is not None else 2.0)
    q = float(driver.exponent_q if driver.exponent_q is not None else p / 2.0)

    for (_, dp), (d2_i, grad_d2_i) in zip(populated, metric_results):
        radius = getattr(dp, "confidence_radius", 1.0)
        if radius is None:
            radius = 1.0
        radius = float(np.asarray(radius, dtype=np.float64).reshape(-1)[0])
        radius = max(radius, 1.0e-8)

        rho2 = radius * radius
        u = max(float(d2_i) / rho2, 1.0e-12)
        denom = u**p + u**q
        w_i = 1.0 / denom
        ddenom_du = p * u**(p - 1.0) + q * u**(q - 1.0)
        dw_dd2 = -(ddenom_du / rho2) / (denom * denom)

        raw.append(w_i)
        raw_grad.append(
            (dw_dd2 * np.asarray(grad_d2_i, dtype=np.float64)).reshape(natm, 3)
        )

    weights, grad_weights = _normalize_mask_weight_pool(driver, raw, raw_grad)
    return weights, grad_weights, populated


def _evaluate_masked_factor_state(
        driver,
        term_bank,
        datapoint,
        projector,
        org_int_coords,
        exact_tol=1.0e-10):
    masked_energy = []
    masked_gradient = []
    raw = []
    raw_grad = []
    exact_flags = []

    natm = driver.impes_coordinate.cartesian_coordinates.shape[0]

    for mask in _coerce_mapping_masks(datapoint):
        mask_arr = np.asarray(mask, dtype=np.int64)
        current = _candidate_current_chart(
            driver,
            org_int_coords,
            driver.impes_coordinate.b_matrix,
            datapoint=datapoint,
            mask=mask_arr,
        )
        energy_m, gradient_m = _restricted_taylor_eval_for_mask(
            driver,
            datapoint=datapoint,
            mask=mask_arr,
            projector=projector,
            org_int_coords=org_int_coords,
            precomputed_current=current,
        )
        d2_m, grad_d2_m = _factor_mask_signature_distance_and_gradient(
            driver,
            term_bank=term_bank,
            datapoint=datapoint,
            mask=mask_arr,
            org_int_coords=org_int_coords,
            b_matrix=current[1],
            precomputed_current=current,
        )

        masked_energy.append(energy_m)
        masked_gradient.append(gradient_m)

        if d2_m <= exact_tol:
            exact_flags.append(True)
            raw.append(0.0)
            raw_grad.append(np.zeros((natm, 3), dtype=np.float64))
        else:
            exact_flags.append(False)
            w_m, grad_w_m = _normalize_symmetry_candidate_weight(
                driver, d2_m, grad_d2_m)
            raw.append(w_m)
            raw_grad.append(grad_w_m)

    masked_energy = np.asarray(masked_energy, dtype=np.float64)
    masked_gradient = np.asarray(masked_gradient, dtype=np.float64)
    exact_idx = np.where(np.asarray(exact_flags, dtype=bool))[0]
    if exact_idx.size > 0:
        weights = np.zeros(masked_energy.size, dtype=np.float64)
        weights[exact_idx] = 1.0 / exact_idx.size
        grad_weights = np.zeros((masked_energy.size, natm, 3), dtype=np.float64)
    else:
        weights, grad_weights = _normalize_mask_weight_pool(driver, raw, raw_grad)

    energy = float(np.dot(weights, masked_energy))
    gradient = (
        np.tensordot(weights, masked_gradient, axes=1)
        + np.tensordot(masked_energy - energy, grad_weights, axes=1)
    )
    return energy, gradient


def _all_factor_mask_signature_distance_and_gradient(
        driver,
        *,
        family_bank,
        datapoint,
        mask,
        org_int_coords,
        b_matrix=None,
        precomputed_current=None):
    mask = np.asarray(mask, dtype=np.int64)

    if precomputed_current is None:
        org_eval, b_eval = _candidate_current_chart(
            driver,
            org_int_coords,
            driver.impes_coordinate.b_matrix if b_matrix is None else b_matrix,
            datapoint=datapoint,
            mask=mask,
        )
    else:
        org_eval, b_eval = precomputed_current

    ref_masked = np.asarray(
        datapoint.internal_coordinates_values, dtype=np.float64)[mask]
    b_eval = np.asarray(b_eval, dtype=np.float64)

    d2_total = 0.0
    grad_total = np.zeros(b_eval.shape[1], dtype=np.float64)
    rows_by_rotor = family_bank.get("primitive_signature_rows", {})

    for rows in rows_by_rotor.values():
        if len(rows) == 0:
            continue
        rows_arr = np.asarray(rows, dtype=np.int64)
        delta = np.asarray(org_eval, dtype=np.float64)[rows_arr] - ref_masked[rows_arr]
        d2_total += float(np.mean(2.0 * (1.0 - np.cos(delta))))
        grad_total += np.mean(
            2.0 * np.sin(delta)[:, None] * b_eval[rows_arr, :],
            axis=0,
        )

    n_groups = len(rows_by_rotor)
    if n_groups > 0:
        d2_total /= n_groups
        grad_total /= n_groups

    return d2_total, grad_total


def _evaluate_masked_factor_core_state(
        driver,
        family_bank,
        core_dp,
        core_projector,
        org_int_coords):
    masked_energy = []
    masked_gradient = []
    raw = []
    raw_grad = []
    exact_flags = []
    natm = driver.impes_coordinate.cartesian_coordinates.shape[0]

    for mask in _coerce_mapping_masks(core_dp):
        mask_arr = np.asarray(mask, dtype=np.int64)
        current = _candidate_current_chart(
            driver,
            org_int_coords,
            driver.impes_coordinate.b_matrix,
            datapoint=core_dp,
            mask=mask_arr,
        )
        energy_m, gradient_m = _restricted_taylor_eval_for_mask(
            driver,
            datapoint=core_dp,
            mask=mask_arr,
            projector=core_projector,
            org_int_coords=org_int_coords,
            precomputed_current=current,
        )
        d2_m, grad_d2_m = _all_factor_mask_signature_distance_and_gradient(
            driver,
            family_bank=family_bank,
            datapoint=core_dp,
            mask=mask_arr,
            org_int_coords=org_int_coords,
            b_matrix=current[1],
            precomputed_current=current,
        )

        masked_energy.append(energy_m)
        masked_gradient.append(gradient_m)
        if d2_m <= 1.0e-10:
            exact_flags.append(True)
            raw.append(0.0)
            raw_grad.append(np.zeros((natm, 3), dtype=np.float64))
        else:
            exact_flags.append(False)
            w_m, grad_w_m = _normalize_symmetry_candidate_weight(
                driver, d2_m, grad_d2_m)
            raw.append(w_m)
            raw_grad.append(grad_w_m)

    masked_energy = np.asarray(masked_energy, dtype=np.float64)
    masked_gradient = np.asarray(masked_gradient, dtype=np.float64)
    exact_idx = np.where(np.asarray(exact_flags, dtype=bool))[0]
    if exact_idx.size > 0:
        weights = np.zeros(masked_energy.size, dtype=np.float64)
        weights[exact_idx] = 1.0 / exact_idx.size
        grad_weights = np.zeros((masked_energy.size, natm, 3), dtype=np.float64)
    else:
        weights, grad_weights = _normalize_mask_weight_pool(driver, raw, raw_grad)

    energy = float(np.dot(weights, masked_energy))
    gradient = (
        np.tensordot(weights, masked_gradient, axes=1)
        + np.tensordot(masked_energy - energy, grad_weights, axes=1)
    )
    return energy, gradient


def evaluate_signed_local_factor_model(driver, family_bank, org_int_coords):
    """
    Evaluate a signed local-factor model for one local-group family.

    The returned energy and gradient already include the core, all positive
    sampled factors, and inclusion-exclusion overlap terms.
    """

    core_dp = family_bank["core"]
    org_coords = np.asarray(org_int_coords, dtype=np.float64)
    n_ic = len(org_coords)

    core_projector = _factor_core_projector(family_bank, n_ic)
    core_coeff = float(family_bank.get("core_coefficient", 0.0))

    total_energy = 0.0
    total_gradient = np.zeros_like(driver.impes_coordinate.cartesian_coordinates)

    if abs(core_coeff) > 1.0e-12 or not family_bank.get("terms"):
        core_energy, core_gradient = _evaluate_masked_factor_core_state(
            driver,
            family_bank,
            core_dp,
            core_projector,
            org_coords,
        )
        total_energy += core_coeff * core_energy
        total_gradient += core_coeff * core_gradient

    evaluated_terms = 0
    for _, term_bank in sorted(
            family_bank.get("terms", {}).items(), key=lambda item: str(item[0])):
        projector = _factor_term_projector(family_bank, term_bank, n_ic)
        weights, grad_weights, states = _compute_factor_weights(
            driver, term_bank, org_coords)
        if not states:
            continue

        state_energy = []
        state_gradient = []
        for _, state_dp in states:
            energy_s, gradient_s = _evaluate_masked_factor_state(
                driver,
                term_bank,
                state_dp,
                projector,
                org_coords,
            )
            state_energy.append(energy_s)
            state_gradient.append(gradient_s)

        state_energy = np.asarray(state_energy, dtype=np.float64)
        state_gradient = np.asarray(state_gradient, dtype=np.float64)
        term_energy = float(np.dot(weights, state_energy))
        term_gradient = (
            np.tensordot(weights, state_gradient, axes=1)
            + np.tensordot(state_energy - term_energy, grad_weights, axes=1)
        )

        coefficient = float(term_bank.get("coefficient", 1.0))
        total_energy += coefficient * term_energy
        total_gradient += coefficient * term_gradient
        evaluated_terms += 1

    if evaluated_terms == 0 and abs(core_coeff) <= 1.0e-12:
        return _evaluate_masked_factor_core_state(
            driver,
            family_bank,
            core_dp,
            core_projector,
            org_coords,
        )

    return float(total_energy), total_gradient
