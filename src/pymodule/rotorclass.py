from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np


@dataclass(frozen=True)
class RotorDefinition:
    rotor_id: int
    center: tuple[int, int]
    torsion_rows: tuple[int, ...]
    torsion_coords: tuple[tuple[int, int, int, int], ...]
    symmetry_order: int
    atom_group: tuple[int, ...]


@dataclass(frozen=True)
class RotorClusterDefinition:
    cluster_id: int
    rotor_ids: tuple[int, ...]
    cluster_type: Literal["independent", "coupled"]
    torsion_rows: tuple[int, ...]


@dataclass
class RotorClusterInformation:
    dihedral_start: int
    dihedral_end: int
    rotors: dict[int, RotorDefinition] = field(default_factory=dict)
    clusters: dict[int, RotorClusterDefinition] = field(default_factory=dict)
    rotor_to_cluster: dict[int, int] = field(default_factory=dict)
    coupling_matrix: np.ndarray | None = None


@dataclass(frozen=True)
class RotorClusterStateDefinition:
    cluster_id: int
    state_id: int
    cluster_type: Literal["independent", "coupled"]
    rotor_ids: tuple[int, ...]
    angle_assignment: dict[tuple[int, int, int, int], float]
    dihedrals_to_rotate: tuple[tuple[int, int, int, int], ...] | None = None
    phase_signature: np.ndarray | None = None
    is_anchor: bool = False
    label_suffix: str = ""


@dataclass
class RotorClusterAngleLibrary:
    state_banks: dict[int, tuple[RotorClusterStateDefinition, ...]] = field(
        default_factory=dict
    )

    def __getitem__(self, cluster_id: int) -> tuple[RotorClusterStateDefinition, ...]:
        return self.state_banks[cluster_id]


@dataclass(frozen=True)
class RotorClusterRuntime:
    cluster_id: int
    cluster_type: Literal["independent", "coupled"]
    rotor_ids: tuple[int, ...]
    cluster_rows: tuple[int, ...]
    projector_rows: tuple[int, ...]
    signature_rows: tuple[int, ...]
    grouped_signature_rows: tuple[tuple[int, ...], ...]
    projector_mask: np.ndarray
    cluster_only_mask: np.ndarray
    signature_mask: np.ndarray


@dataclass
class RotorClusterRuntimeSet:
    n_ic: int
    core_rows: tuple[int, ...]
    core_mask: np.ndarray
    clusters: dict[int, RotorClusterRuntime] = field(default_factory=dict)


def rotor_coupling_score(H, rotor_a_rows, rotor_b_rows, eps=1.0e-12):
    H_ab = H[np.ix_(rotor_a_rows, rotor_b_rows)]
    H_aa = H[np.ix_(rotor_a_rows, rotor_a_rows)]
    H_bb = H[np.ix_(rotor_b_rows, rotor_b_rows)]
    num = np.linalg.norm(H_ab, ord="fro")
    den = np.sqrt(np.linalg.norm(H_aa, ord="fro") * np.linalg.norm(H_bb, ord="fro") + eps)
    return float(num / den)


def build_rotor_clusters(rotors, coupling_matrix, threshold):
    unvisited = set(rotors.keys())
    clusters = []
    while unvisited:
        seed = unvisited.pop()
        component = {seed}
        stack = [seed]
        while stack:
            r = stack.pop()
            neighbors = {s for s in unvisited if coupling_matrix[r, s] >= threshold}
            component.update(neighbors)
            stack.extend(neighbors)
            unvisited.difference_update(neighbors)
        clusters.append(tuple(sorted(component)))
    return clusters


def build_runtime_cluster_set(cluster_info: RotorClusterInformation, n_ic: int):
    all_cluster_rows = []
    for cluster in cluster_info.clusters.values():
        all_cluster_rows.extend(cluster.torsion_rows)

    core_rows = tuple(sorted(set(range(n_ic)) - set(all_cluster_rows)))
    core_mask = np.zeros(n_ic, dtype=np.float64)
    core_mask[list(core_rows)] = 1.0

    runtime = RotorClusterRuntimeSet(
        n_ic=n_ic,
        core_rows=core_rows,
        core_mask=core_mask,
        clusters={},
    )

    for cluster_id, cluster in cluster_info.clusters.items():
        grouped_signature_rows = []
        cluster_rows = []
        for rotor_id in sorted(cluster.rotor_ids):
            rows_r = tuple(cluster_info.rotors[rotor_id].torsion_rows)
            grouped_signature_rows.append(rows_r)
            cluster_rows.extend(rows_r)

        cluster_rows = tuple(cluster_rows)
        projector_rows = tuple(sorted(set(core_rows).union(cluster_rows)))
        signature_rows = cluster_rows

        projector_mask = np.zeros(n_ic, dtype=np.float64)
        projector_mask[list(projector_rows)] = 1.0

        cluster_only_mask = np.zeros(n_ic, dtype=np.float64)
        cluster_only_mask[list(cluster_rows)] = 1.0

        signature_mask = np.zeros(n_ic, dtype=np.float64)
        signature_mask[list(signature_rows)] = 1.0

        runtime.clusters[cluster_id] = RotorClusterRuntime(
            cluster_id=cluster_id,
            cluster_type=cluster.cluster_type,
            rotor_ids=tuple(sorted(cluster.rotor_ids)),
            cluster_rows=cluster_rows,
            projector_rows=projector_rows,
            signature_rows=signature_rows,
            grouped_signature_rows=tuple(grouped_signature_rows),
            projector_mask=projector_mask,
            cluster_only_mask=cluster_only_mask,
            signature_mask=signature_mask,
        )
    return runtime