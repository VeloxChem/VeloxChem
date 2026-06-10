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


from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

"""
This class is designed to be a helper function for the rotor-framework
within the interpolation framework for symmetry considerations of methyl
groups. The framework is intended to be furhter extended to more symmetry
equivalent groups which is feasible.
"""

@dataclass(frozen=True)
class RotorDefinition:
    """
    Object class used for handling rotor specific variables
    for correct integration within the interpolation driver
    """
    # unique identifier of the rotor
    rotor_id: int 
    # rotatable bond
    center: tuple[int, int]
    # all indices of participating dohedrals within the z-matrix
    torsion_rows: tuple[int, ...]
    # acutal dihedrals
    torsion_coords: tuple[tuple[int, int, int, int], ...]
    # what kind of periodicity --> currently always 3 for methyl rotors
    symmetry_order: int
    # atom indices defining the symmetry-equivalent group attached to the rotor
    atom_group: tuple[int, ...]


@dataclass(frozen=True)
class RotorClusterDefinition:
    """
    Possible categorization of rotor combinations into clusters
    independent --> Cluster with a single rotor
    coupled --> Cluster containing all single rotor combinations due to coupling
    """
    # unique identifier of this rotor cluster
    cluster_id: int
    # all rotor IDs grouped into this cluster
    rotor_ids: tuple[int, ...]
    # cluster category: one rotor (independent) or multiple coupled rotors
    cluster_type: Literal["independent", "coupled"]
    # internal-coordinate rows covered by the rotors in this cluster
    torsion_rows: tuple[int, ...]


@dataclass
class RotorClusterInformation:
    """
    Container holding all rotor definitions, cluster assignments, and
    coupling metadata for one interpolation family/state.
    """
    # start index of the dihedral block in internal-coordinate vectors
    dihedral_start: int
    # end index (exclusive) of the dihedral block in internal-coordinate vectors
    dihedral_end: int
    # lookup table: rotor ID -> rotor definition
    rotors: dict[int, RotorDefinition] = field(default_factory=dict)
    # lookup table: cluster ID -> cluster definition
    clusters: dict[int, RotorClusterDefinition] = field(default_factory=dict)
    # lookup table: rotor ID -> owning cluster ID
    rotor_to_cluster: dict[int, int] = field(default_factory=dict)
    # pairwise rotor-coupling strengths used to build cluster components
    coupling_matrix: np.ndarray | None = None


@dataclass(frozen=True)
class RotorClusterStateDefinition:
    """
    Defines one discrete state in a cluster angle library, including the
    assigned torsion values and metadata used for labeling/evaluation.
    """
    # identifier of the cluster this state belongs to
    cluster_id: int
    # state index inside the cluster bank (state 0 is the anchor state)
    state_id: int
    # cluster category copied for fast branching in downstream logic
    cluster_type: Literal["independent", "coupled"]
    # rotor IDs participating in this state
    rotor_ids: tuple[int, ...]
    # target torsion values used to generate/evaluate this cluster state
    angle_assignment: dict[tuple[int, int, int, int], float]
    # dihedrals actively rotated when constructing this state from the anchor
    dihedrals_to_rotate: tuple[tuple[int, int, int, int], ...] | None = None
    # optional grouped phase descriptor for signature-based matching
    phase_signature: np.ndarray | None = None
    # marks the unrotated reference state for the cluster
    is_anchor: bool = False
    # optional suffix used when building persistent point labels
    label_suffix: str = ""


@dataclass
class RotorClusterAngleLibrary:
    """
    Registry of pre-enumerated cluster states keyed by cluster ID.
    """
    # mapping from cluster ID to the ordered state bank for that cluster
    state_banks: dict[int, tuple[RotorClusterStateDefinition, ...]] = field(
        default_factory=dict
    )

    def __getitem__(self, cluster_id: int) -> tuple[RotorClusterStateDefinition, ...]:
        return self.state_banks[cluster_id]


@dataclass(frozen=True)
class RotorClusterRuntime:
    """
    Precomputed runtime masks and row-partitions for one rotor cluster,
    used during projected Taylor interpolation and signature matching.
    """
    # cluster identifier for this runtime view
    cluster_id: int
    # cluster category: independent or coupled
    cluster_type: Literal["independent", "coupled"]
    # rotor IDs represented by this runtime cluster
    rotor_ids: tuple[int, ...]
    # all internal-coordinate rows owned by this cluster
    cluster_rows: tuple[int, ...]
    # rows retained in restricted Taylor expansion (core + cluster rows)
    projector_rows: tuple[int, ...]
    # rows used to build periodic cluster-signature distances
    signature_rows: tuple[int, ...]
    # signature rows grouped rotor-by-rotor for per-rotor phase metrics
    grouped_signature_rows: tuple[tuple[int, ...], ...]
    # dense projector mask over all internal coordinates
    projector_mask: np.ndarray
    # dense mask containing only rows that belong to this cluster
    cluster_only_mask: np.ndarray
    # dense mask containing only rows used for signature matching
    signature_mask: np.ndarray


@dataclass
class RotorClusterRuntimeSet:
    """
    Runtime bundle containing the global core projector and per-cluster
    runtime projectors/signature masks for one internal-coordinate size.
    """
    # total number of internal coordinates for which masks are defined
    n_ic: int
    # internal-coordinate rows that are not part of any rotor cluster
    core_rows: tuple[int, ...]
    # dense mask selecting the non-cluster (core) rows
    core_mask: np.ndarray
    # mapping from cluster ID to runtime masks/signature metadata
    clusters: dict[int, RotorClusterRuntime] = field(default_factory=dict)


def rotor_coupling_score(H, rotor_a_rows, rotor_b_rows, eps=1.0e-12):
    """
    Computes a normalized Hessian coupling score between two rotor blocks.
    """
    # Hessian cross-block between rotor A rows and rotor B rows.
    H_ab = H[np.ix_(rotor_a_rows, rotor_b_rows)]
    # Hessian self-block for rotor A.
    H_aa = H[np.ix_(rotor_a_rows, rotor_a_rows)]
    # Hessian self-block for rotor B.
    H_bb = H[np.ix_(rotor_b_rows, rotor_b_rows)]
    # Frobenius norm of inter-rotor coupling.
    num = np.linalg.norm(H_ab, ord="fro")
    # Normalization from intra-rotor block norms (regularized by eps).
    den = np.sqrt(np.linalg.norm(H_aa, ord="fro") * np.linalg.norm(H_bb, ord="fro") + eps)
    return float(num / den)


def build_rotor_clusters(rotors, coupling_matrix, threshold):
    """
    Builds connected rotor components from a coupling matrix threshold.
    """
    # Rotor IDs that are not yet assigned to any cluster.
    unvisited = set(rotors.keys())
    # List of discovered rotor clusters (connected components).
    clusters = []
    while unvisited:
        # Seed rotor used to start exploring a new connected component.
        seed = unvisited.pop()
        # Current connected component being grown.
        component = {seed}
        # DFS/LIFO stack of frontier rotor IDs.
        stack = [seed]
        while stack:
            # Rotor ID currently being expanded.
            r = stack.pop()
            # Unvisited rotors coupled strongly enough to be linked with r.
            neighbors = {s for s in unvisited if coupling_matrix[r, s] >= threshold}
            # Add discovered neighboring rotors to current component.
            component.update(neighbors)
            # Continue exploring from newly found neighbors.
            stack.extend(neighbors)
            # Remove assigned neighbors so they are not revisited.
            unvisited.difference_update(neighbors)
        # Store cluster as a sorted tuple for stable downstream indexing.
        clusters.append(tuple(sorted(component)))
    return clusters


def build_runtime_cluster_set(cluster_info: RotorClusterInformation, n_ic: int):
    """
    Precomputes core/cluster masks and row groupings used at runtime.
    """
    # Flattened list of all internal-coordinate rows participating in any cluster.
    all_cluster_rows = []
    for cluster in cluster_info.clusters.values():
        all_cluster_rows.extend(cluster.torsion_rows)

    # Internal-coordinate rows not owned by any rotor cluster.
    core_rows = tuple(sorted(set(range(n_ic)) - set(all_cluster_rows)))
    # Dense selector mask for core rows in projected Taylor evaluation.
    core_mask = np.zeros(n_ic, dtype=np.float64)
    core_mask[list(core_rows)] = 1.0

    # Runtime container holding global core projection and per-cluster runtime metadata.
    runtime = RotorClusterRuntimeSet(
        n_ic=n_ic,
        core_rows=core_rows,
        core_mask=core_mask,
        clusters={},
    )

    for cluster_id, cluster in cluster_info.clusters.items():
        # Per-rotor grouped row sets used by cluster phase-signature distances.
        grouped_signature_rows = []
        # Flattened row list for all rotors in this cluster.
        cluster_rows = []
        for rotor_id in sorted(cluster.rotor_ids):
            # Rows associated with this single rotor.
            rows_r = tuple(cluster_info.rotors[rotor_id].torsion_rows)
            grouped_signature_rows.append(rows_r)
            cluster_rows.extend(rows_r)

        # Stable tuple view of all rows belonging to the cluster.
        cluster_rows = tuple(cluster_rows)
        # Rows retained during restricted evaluation: core rows + cluster rows.
        projector_rows = tuple(sorted(set(core_rows).union(cluster_rows)))
        # Rows used when building periodic signatures for this cluster.
        signature_rows = cluster_rows

        # Dense projector mask selecting core + cluster rows.
        projector_mask = np.zeros(n_ic, dtype=np.float64)
        projector_mask[list(projector_rows)] = 1.0

        # Dense mask selecting only rows that belong to this cluster.
        cluster_only_mask = np.zeros(n_ic, dtype=np.float64)
        cluster_only_mask[list(cluster_rows)] = 1.0

        # Dense mask selecting rows that contribute to signature comparison.
        signature_mask = np.zeros(n_ic, dtype=np.float64)
        signature_mask[list(signature_rows)] = 1.0

        # Register runtime metadata for this cluster ID.
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
