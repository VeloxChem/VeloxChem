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

@dataclass(frozen=True)
class LocalRotor:
    rotor_id: str
    kind: str                 # methyl, tertbutyl_parent, phenyl, generic
    axis: tuple[int, int]
    symmetry_order: int       # 1 for non-symmetric
    owned_atoms: tuple[int, ...]
    unit_atom_sets: tuple[tuple[int, ...], ...]
    torsion_rows: tuple[int, ...]
    phase_coordinate: tuple[int, int, int, int] | None = None
    phase_values: tuple[float, ...] | None = None
    signature_rows: tuple[int, ...] = ()
    signature_row_types: tuple[str, ...] = ()
    signature_row_scales: tuple[float, ...] = ()

@dataclass(frozen=True)
class LocalCluster:
    cluster_id: str
    family_label: str
    rotor_ids: tuple[str, ...]
    owned_atoms: tuple[int, ...]
    active_atoms: tuple[int, ...]
    active_rows: tuple[int, ...]
    anchor_atoms: tuple[int, ...] = ()
    relaxation_atoms: tuple[int, ...] = ()
    role: str = "factor"
    canonical_subset_key: tuple[str, ...] = ()
    parent_cluster_ids: tuple[str, ...] = ()
    response_rows: tuple[int, ...] = ()
    projector_rows: tuple[int, ...] = ()
    projector_policy_id: str = "legacy"
    relaxation_policy_id: str = "default"
    anchor_state_ids: tuple[int, ...] = ()

@dataclass
class LocalGroupModel:
    version: int = 1
    enabled: bool = False
    rotors: dict[str, LocalRotor] = field(default_factory=dict)
    clusters: dict[str, LocalCluster] = field(default_factory=dict)
    core_atoms: tuple[int, ...] = ()
    core_rows: tuple[int, ...] = ()
