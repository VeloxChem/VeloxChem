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

import sys
from pathlib import Path

import numpy as np
import networkx as nx
from mpi4py import MPI
import h5py
import re

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:
    pass

from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from ...molecule import Molecule

from ..io.basic import nn, nl, pname, is_list_A_in_B, lname, arr_dimension
from ..utils.geometry import (unit_cell_to_cartesian_matrix,
                              fractional_to_cartesian, cartesian_to_fractional,
                              locate_min_idx, reorthogonalize_matrix,
                              find_optimal_pairings, find_edge_pairings,
                              Carte_points_generator)
from .other import fetch_X_atoms_ind_array, find_pair_x_edge_fc, order_edge_array
from .superimpose import superimpose_rotation_only


class SupercellBuilder:
    """Build supercell graph from primitive cell graph sG and optionally add virtual edges.

    For ditopic linkers, expands sG by supercell and checks virtual edges. For multitopic,
    bundles multiedges (CV nodes), then expands and updates bundle in supercell.
    """

    def __init__(self, comm=None, ostream=None):
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)

        #need to be set before use
        self.sG = None
        self.linker_connectivity = None
        self.supercell = [1, 1, 1]
        self.cell_info = None  #a list of 6 elements, a,b,c,alpha,beta,gamma
        #virtual edge settings for bridge type nodes
        self.add_virtual_edge = False
        self.vir_edge_range = 0.5  #in fractional coordinate should be less than 0.5
        self.vir_edge_max_neighbor = 2

        #will be generated during the process
        self.multiedge_bundlings = None
        self.superG = None

        self._debug = False

    def _is_ditopic_linker(self):
        """Return True if linker is ditopic (connectivity 2)."""
        return self.linker_connectivity == 2

    def build_supercellGraph(self):
        """Build superG from sG: expand by supercell, (for multitopic) bundle multiedges, optionally add virtual edges. Returns superG."""
        self.ostream.print_info("building supercell graph...")
        self.ostream.flush()
        self.superG_cell_info = [
            self.cell_info[0] * self.supercell[0],
            self.cell_info[1] * self.supercell[1],
            self.cell_info[2] * self.supercell[2],
            self.cell_info[3],
            self.cell_info[4],
            self.cell_info[5],
        ]
        self.ostream.print_info(
            f"supercell cell info: {self.superG_cell_info}")
        self.ostream.flush()
        sG = self.sG.copy()
        if not self._is_ditopic_linker():
            self.multiedge_bundlings = self._bundle_multiedge(sG)

        superG = self._update_supercell_node_fpoints_loose(sG, self.supercell)
        superG = self._update_supercell_edge_fpoints(sG, superG,
                                                     self.supercell)
        if self._debug:
            self.ostream.print_info(
                f"nodes in supercell graph superG: {len(superG.nodes())}")
            self.ostream.print_info(
                f"edges in supercell graph superG: {len(superG.edges())}")
            self.ostream.flush()

        if self._is_ditopic_linker():
            self.superG = self._check_virtual_edge(superG.copy())
            return self.superG
        # for multitopic linker MOF
        self.prim_multiedge_bundlings = self.multiedge_bundlings
        self.super_multiedge_bundlings = self._make_super_multiedge_bundlings(
            self.prim_multiedge_bundlings, self.supercell)
        superG = self._update_supercell_bundle(superG,
                                               self.super_multiedge_bundlings)
        new_superG = self._check_multiedge_bundlings_insuperG(
            self.super_multiedge_bundlings, superG)
        self.superG = self._check_virtual_edge(new_superG.copy())
        self.ostream.print_info("supercell graph built")
        self.ostream.flush()
        return self.superG

    def _check_virtual_edge(self, supG):
        if self._debug:
            self.ostream.print_info(
                "checking virtual edge for bridge nodes...")
            self.ostream.print_info(
                f"current edges in superG: {len(supG.edges())}")
            self.ostream.flush()

        if self.add_virtual_edge:
            superG = self._add_virtual_edge(
                self.sc_unit_cell,
                supG.copy(),
                self.vir_edge_range,
                self.vir_edge_max_neighbor,
            )
            if self._debug:
                self.ostream.print_info("added virtual edge for bridge nodes")
                self.ostream.print_info(
                    f"current edges in superG after adding virtual edge:{len(superG.edges())}"
                )
                self.ostream.flush()
        else:
            superG = supG.copy()
        return superG

    def _bundle_multiedge(self, sG):
        multiedge_bundlings = []
        for n in sG.nodes:
            if "CV" in n and sG.nodes[n]["type"] == "V":
                edges = []
                for con_n in list(sG.neighbors(n)):
                    edges.append(sG.edges[n, con_n]["coords"])
                multiedge_bundlings.append((n, list(sG.neighbors(n))))
        return multiedge_bundlings

    # check if node is at the boundary of the supercell
    # if at boundary, then add the diff_e to the node name and add the diff_e to the f_points
    def _update_supercell_node_fpoints_loose(self, sG, supercell):
        # boundary_node_res = []
        # incell_node_res = []
        superG = nx.Graph()
        for n in sG.nodes():
            if sG.nodes[n]["type"] != "V":  # get rid of SV, sv will be pnode+i
                superG.add_node(
                    n,
                    f_points=sG.nodes[n]["f_points"],
                    fcoords=sG.nodes[n]["fcoords"],
                    type="SV",
                    note=sG.nodes[n]["note"],
                )

                continue

            # add the node to superG, if lname(n) is np.array([0,0,0])
            superG.add_node(
                pname(n) + "_" + str(lname(n)),
                f_points=sG.nodes[n]["f_points"],
                fcoords=sG.nodes[n]["fcoords"],
                type="V",
                note=sG.nodes[n]["note"],
            )

            diffs = (np.mod(sG.nodes[n]["fcoords"], 1) -
                     sG.nodes[n]["fcoords"] + np.asarray(supercell))
            diffs = diffs.astype(int)
            diff_ele = Carte_points_generator(diffs)
            diff_ele = diff_ele[1:]  # remove the first element [0,0,0]

            # if len(diff_ele) > supercell_Carte.shape[0]:
            #    boundary_node_res.append(n)
            # else:
            #    incell_node_res.append(n)

            for diff_e in diff_ele:
                diff_e = np.asarray(diff_e)
                if (pname(n) + "_" + str(lname(n) + diff_e)) in superG.nodes():
                    self.ostream.print_info(
                        f"node already in superG {pname(n) + '_' + str(lname(n) + diff_e)}"
                    )
                    self.ostream.flush()
                    continue
                superG.add_node(
                    (pname(n) + "_" + str(lname(n) + diff_e)),
                    f_points=np.hstack((
                        sG.nodes[n]["f_points"][:, 0:2],
                        sG.nodes[n]["f_points"][:, 2:5].astype(float) + diff_e,
                    )),  # NOTE:modified because of extra column of atom type
                    fcoords=sG.nodes[n]["fcoords"] + diff_e,
                    type="SV",
                    note=sG.nodes[n]["note"],
                )

        return superG

    """
    sG_node_note_set:{'CV', 'V'}
    sG_node_type_set:{'DV', 'V'}
    sG_edge_type_set:{'DE', 'E'}
    """

    # check if edge is at the boundary of the supercell
    def _update_supercell_edge_fpoints(self, sG, superG, supercell):
        # boundary edge is DE
        # incell edge is E
        supercell_Carte = Carte_points_generator(supercell)
        for e in sG.edges():
            for i in supercell_Carte:
                s_edge = (
                    pname(e[0]) + "_" + str(i + lname(e[0])),
                    pname(e[1]) + "_" + str(i + lname(e[1])),
                )

                # check if node e[0]+'_'+str(diff_e) and e[1]+'_'+str(diff_e) in superG
                if (s_edge[0] in superG.nodes()) and (s_edge[1]
                                                      in superG.nodes()):
                    superG.add_edge(
                        s_edge[0],
                        s_edge[1],
                        f_points=np.hstack((
                            sG.edges[e]["f_points"][:, 0:2],
                            sG.edges[e]["f_points"][:, 2:5].astype(float) + i,
                        )),
                        fcoords=sG.edges[e]["fcoords"] + i,
                        type=sG.edges[e]["type"],
                    )

                elif (s_edge[0] in superG.nodes()) or (s_edge[1]
                                                       in superG.nodes()):
                    if s_edge[0] in superG.nodes():
                        superG.add_node(
                            s_edge[1],
                            f_points=np.hstack((
                                sG.nodes[e[1]]["f_points"][:, 0:2],
                                sG.nodes[e[1]]["f_points"][:,
                                                           2:5].astype(float) +
                                i,
                            )),
                            fcoords=sG.nodes[e[1]]["fcoords"] + i,
                            type="DSV",
                            note=sG.nodes[e[1]]["note"],
                        )

                    else:
                        superG.add_node(
                            s_edge[0],
                            f_points=np.hstack((
                                sG.nodes[e[0]]["f_points"][:, 0:2],
                                sG.nodes[e[0]]["f_points"][:,
                                                           2:5].astype(float) +
                                i,
                            )),
                            fcoords=sG.nodes[e[0]]["fcoords"] + i,
                            type="DSV",
                            note=sG.nodes[e[0]]["note"],
                        )

                    superG.add_edge(
                        s_edge[0],
                        s_edge[1],
                        f_points=np.hstack(
                            (
                                sG.edges[e]["f_points"][:, 0:2],
                                sG.edges[e]["f_points"][:, 2:5].astype(float) +
                                i,
                            )
                        ),  # NOTE:modified because of extra column of atom type
                        fcoords=sG.edges[e]["fcoords"] + i,
                        type="DSE",
                    )

                else:
                    if self._debug:
                        self.ostream.print_info(
                            f"edge not in superG: {s_edge[0]}, {s_edge[1]}")
                        self.ostream.flush()
        return superG

    ########## the below is to process the multiedge bundling in superG###########
    # need to combine with the multiedge_bundling.py
    # replace bundle dvnode with vnode+diff_e
    def _replace_bundle_dvnode_with_vnode(self, dv_v_pairs,
                                          multiedge_bundlings):
        for dv, v in dv_v_pairs:
            for bund in multiedge_bundlings:
                if dv in bund[1]:
                    if pname(bund[1][bund[1].index(dv)]) == pname(dv):
                        bund[1][bund[1].index(dv)] = v
                if dv in bund[0]:
                    if pname(bund[0][bund[0].index(dv)]) == pname(dv):
                        bund[0][bund[0].index(dv)] = v
        # update v if no list then add [0,0,0]
        # convert tuple to list
        updated_bundlings = []
        for bund in multiedge_bundlings:
            ec_node = pname(bund[0]) + "_" + str(lname(bund[0]))
            con_nodes = [pname(i) + "_" + str(lname(i)) for i in bund[1]]
            updated_bundlings.append((ec_node, con_nodes))
        return updated_bundlings

    # loop bundle and check if any element in the bundle is in the superG, if not, add the element to the superG
    def _make_super_multiedge_bundlings(self, prim_multiedge_bundlings,
                                        supercell):
        super_multiedge_bundlings = {}
        for i in Carte_points_generator(supercell):
            for bund in prim_multiedge_bundlings:
                ec_node = pname(bund[0]) + "_" + str(i + lname(bund[0]))
                con_nodes = [
                    pname(n) + "_" + str(i + lname(n)) for n in bund[1]
                ]
                super_multiedge_bundlings[ec_node] = con_nodes
        return super_multiedge_bundlings

    def _update_supercell_bundle(self, superG, super_multiedge_bundlings):
        if self._debug:
            self.ostream.print_info("updating supercell bundle...")
            self.ostream.print_info(
                f"current nodes number in superG: {len(superG.nodes())}")
            self.ostream.flush()
        for ec_node in super_multiedge_bundlings.keys():
            con_nodes = super_multiedge_bundlings[ec_node]
            # order the con_nodes by th x-x pair of the ecnode X atoms
            prim_ecname = pname(ec_node) + "_" + str(np.array([0.0, 0.0, 0.0]))
            if ec_node not in superG.nodes():
                trans = lname(ec_node)
                superG.add_node(
                    ec_node,
                    f_points=np.hstack(
                        (superG.nodes[prim_ecname]["f_points"][:, 0:2],
                         superG.nodes[prim_ecname]["f_points"][:, 2:5].astype(
                             float) + trans)),
                    fcoords=superG.nodes[prim_ecname]["fcoords"] + trans,
                    type="SV",
                    note=superG.nodes[prim_ecname]["note"],
                )
            for j in range(len(con_nodes)):
                cn = con_nodes[j]
                prim_cnname = super_multiedge_bundlings[prim_ecname][
                    j]  # find prim_ecname in super_multiedge_bundlings and then get the corresponding prim_cnname
                trans = lname(cn) - lname(prim_cnname)
                if cn not in superG.nodes():
                    superG.add_node(
                        cn,
                        f_points=np.hstack(
                            (superG.nodes[prim_cnname]["f_points"][:, 0:2],
                             superG.nodes[prim_cnname]["f_points"]
                             [:, 2:5].astype(float) + trans)),
                        fcoords=superG.nodes[prim_cnname]["fcoords"] + trans,
                        type="SV",
                        note=superG.nodes[prim_cnname]["note"],
                    )
                superG.add_edge(
                    ec_node,
                    cn,
                    f_points=np.hstack(
                        (superG.edges[prim_ecname,
                                      prim_cnname]["f_points"][:, 0:2],
                         superG.edges[prim_ecname, prim_cnname]["f_points"]
                         [:, 2:5].astype(float) + trans)),
                    fcoords=superG.edges[prim_ecname, prim_cnname]["fcoords"] +
                    trans,
                    type="DSE",
                )

        return superG

    def _check_multiedge_bundlings_insuperG(self, super_multiedge_bundlings,
                                            superG):
        super_multiedge_bundlings_edges = []
        for ec_node in super_multiedge_bundlings:
            # check is all CV node in superG are in the super_multiedge_bundlings_edges first element
            cvnodes = [
                n for n in superG.nodes() if superG.nodes[n]["note"] == "CV"
            ]
            # use set to check if all cvnodes are in the super_multiedge_bundlings_edges
            if set(cvnodes) == set(
                [i[0] for i in super_multiedge_bundlings_edges]):
                return superG
            else:
                self.ostream.print_info(
                    "not all CV nodes in super_multiedge_bundlings_edges")
                self.ostream.flush()
                diff_element = set(cvnodes).difference(
                    set(list(super_multiedge_bundlings)))
                if self._debug:
                    self.ostream.print_info(
                        f"removing diff_element: {diff_element}")
                    self.ostream.flush()
                # remove the diff_element from the superG
                for n in diff_element:
                    superG.remove_node(n)
                    # remove all edges linked to the node
                    edges = [e for e in superG.edges(n)]
                    for e in edges:
                        superG.remove_edge(e[0], e[1])

                return superG

    def _add_virtual_edge(self,
                          unit_cell,
                          superG,
                          bridge_node_distance,
                          max_neighbor=2):
        # add pillar nodes virtual edges
        nodes_list = [
            n for n in superG.nodes() if superG.nodes[n]["note"] == "V"
        ]
        n_n_distance_matrix = np.zeros((len(nodes_list), len(nodes_list)))

        for i in range(len(nodes_list)):
            for j in range(len(nodes_list)):
                n_n_distance_matrix[i, j] = np.linalg.norm(
                    np.dot(
                        unit_cell,
                        superG.nodes[nodes_list[i]]["fcoords"] -
                        superG.nodes[nodes_list[j]]["fcoords"],
                    ))
            n_n_distance_matrix[i, i] = 1000

        # find the shortest path between all nodes
        for i in range(len(nodes_list)):
            neighbor_count = 0
            while neighbor_count < max_neighbor:

                def add_v_e(i, n_n_distance_matrix, superG,
                            bridge_node_distance, count):
                    n_n_min_distance = np.min(n_n_distance_matrix[i:i + 1, :])
                    if n_n_min_distance < bridge_node_distance:
                        _, n_j = locate_min_idx(n_n_distance_matrix[i:i +
                                                                    1, :])
                        superG.add_edge(nodes_list[i],
                                        nodes_list[n_j],
                                        type="virtual")
                        if self._debug:
                            self.ostream.print_info(
                                f"add virtual edge between {nodes_list[i]}, {nodes_list[n_j]}"
                            )
                            self.ostream.flush()
                        n_n_distance_matrix[i, n_j] = 1000
                        return True, count + 1, n_n_distance_matrix, superG
                    else:
                        return False, count, n_n_distance_matrix, superG

                added, neighbor_count, n_n_distance_matrix, superG = add_v_e(
                    i, n_n_distance_matrix, superG, bridge_node_distance,
                    neighbor_count)
                if not added:
                    break

        return superG

    def _add_virtual_edge_for_bridge_node(self, superG):
        """
        after setting the virtual edge search, add the virtual edge to the target supercell superG MOF
        """
        if self.add_virtual_edge:
            add_superG = self._add_virtual_edge(
                self.sc_unit_cell,
                superG,
                self.vir_edge_range,
                self.vir_edge_max_neighbor,
            )
            if self._debug():
                self.ostream.print_info("add virtual edge")
                self.ostream.flush()
            return add_superG
        else:
            return superG


class EdgeGraphBuilder:
    """Build edge graph (eG) from supercell graph superG: V/EDGE nodes, XOO attachment, and optional cleaving.

    Converts superG to eG (ditopic or multitopic), adds XOO from nodes to edges,
    cleaves to custom_fbox/supercell range, and sets matched_vnode_xind and unsaturated lists.
    """

    def __init__(self, comm=None, ostream=None):
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)

        #need to be set before use
        self.superG = None
        self.linker_connectivity = None
        self.sc_unit_cell = None  #3x3 matrix of the supercell unit cell
        self.supercell = None  #same as the supercell in supercellbuilder
        self.node_connectivity = None
        #specific buffer range for cleaving the supercell
        self.custom_fbox = None  #in fractional coordinate

        #will be generated during the process
        self.eG = None
        self.cleaved_eG = None
        self.cleaved_edges = None
        self.cleaved_nodes = None
        self.matched_vnode_xind = None
        self.xoo_dict = None
        self.unsaturated_nodes = []
        self.unsaturated_linkers = []

        self._debug = False

    def build_edgeG_from_superG(self):
        """Build eG from superG, add XOO to edges, cleave to range, set unsaturated nodes/linkers."""
        if self.linker_connectivity == 2:
            eG, eG_index_name_dict = self._superG_to_eG_ditopic(self.superG)
        else:
            eG, eG_index_name_dict = self._superG_to_eG_multitopic(
                self.superG, self.sc_unit_cell)

        # only keep the main fragment of the target MOF cell, remove the other fragments, to avoid the disconnected fragments
        self.ostream.print_info("main fragment of the eG kept")
        self.ostream.flush()
        self.eG = [eG.subgraph(c).copy() for c in nx.connected_components(eG)
                   ][0]  #keep the main fragment of the eG
        if self.linker_connectivity == 2:
            self.eG, unsaturated_linker, self.matched_vnode_xind, self.xoo_dict = self._addxoo2edge_ditopic(
                self.eG, self.sc_unit_cell)

        else:
            self.eG, unsaturated_linker, self.matched_vnode_xind, self.xoo_dict = self._addxoo2edge_multitopic(
                self.eG, self.sc_unit_cell)

        #cleave the range of supercell buffer
        self._make_supercell_range_cleaved_eG(self.eG, self.supercell,
                                              self.custom_fbox)
        self.eG_index_name_dict = eG_index_name_dict
        #will generate the cleaved_eG, cleaved_edges, cleaved_nodes
        self.unsaturated_linkers = unsaturated_linker

        self.unsaturated_nodes = self._find_unsaturated_node(
            self.eG, self.node_connectivity)

    def _find_unsaturated_node(self, eG, node_connectivity):
        # find unsaturated node V in eG
        unsaturated_node = []
        for n in eG.nodes():
            if pname(n) != "EDGE":
                real_neighbor = []
                for cn in eG.neighbors(n):
                    if eG.edges[(n, cn)]["type"] == "real":
                        real_neighbor.append(cn)
                if len(real_neighbor) < node_connectivity:
                    unsaturated_node.append(n)
        return unsaturated_node

    def _superG_to_eG_multitopic(self, superG, sc_unit_cell):
        """
        Convert a multitopic supercell graph (superG) into an edge graph (eG).

        - V nodes become NODE entries (keeps index on superG for later reference).
        - CV nodes become EDGE entries: incident edge f_points are merged and
          X->x pairings are resolved via make_paired_Xto_x.
        - Virtual edges in superG are preserved in the resulting eG.
        """
        eG = nx.Graph()
        edge_count = 1  #2n
        node_count = 0  #2n+1
        eG_index_name_dict = {}

        # First pass: add NODE entries for V nodes and EDGE nodes for CV centers
        for n in superG.nodes():
            node_data = superG.nodes[n]
            note = node_data.get("note")

            if note == "V":
                eG.add_node(
                    n,
                    f_points=node_data.get("f_points"),
                    fcoords=node_data.get("fcoords"),
                    type=node_data.get("type", "V"),
                    name="NODE",
                    note="V",
                    index=2 * node_count + 1,  # odd indices for NODEs
                )
                # keep index on superG for cross-referencing downstream
                superG.nodes[n]["index"] = 2 * node_count + 1
                eG_index_name_dict[2 * node_count + 1] = n  # map index to name
                node_count += 1

            elif note == "CV":
                neighbors = list(superG.neighbors(n))
                if not neighbors:
                    continue

                # Collect edge f_points from all incident edges
                fp_list = []
                for nb in neighbors:
                    fp = superG.edges[n, nb].get("f_points")
                    if fp is None:
                        continue
                    fp_list.append(np.asarray(fp))

                merged_edges = np.vstack(fp_list) if fp_list else np.empty(
                    (0, 0))

                # Pair X atoms between the CV center and neighbor edges, converting matched X -> x
                ec_merged_edges = self._make_paired_Xto_x(
                    node_data.get("f_points"), merged_edges, len(neighbors),
                    sc_unit_cell)

                edge_name = "EDGE_" + str(
                    2 * edge_count)  # EDGE_2, EDGE_4, ...
                eG.add_node(
                    edge_name,
                    f_points=ec_merged_edges,
                    fcoords=node_data.get("fcoords"),
                    type="Edge",
                    name="EDGE",
                    note="E",
                    index=2 * edge_count,  # even indices for EDGEs
                )
                eG_index_name_dict[2 *
                                   edge_count] = edge_name  # map index to name

                # connect EDGE to each neighbor with a 'real' half-edge
                for nb in neighbors:
                    eG.add_edge(edge_name,
                                nb,
                                index="E_" + str(2 * edge_count),
                                type="real")
                edge_count += 1

        # Preserve virtual edges present in superG
        for u, v in superG.edges():
            if superG.edges[u, v].get("type") == "virtual":
                # avoid duplicating edges already created
                if not eG.has_edge(u, v):
                    eG.add_edge(u, v, type="virtual")

        return eG, eG_index_name_dict

    def _superG_to_eG_ditopic(self, superG):
        """
        Convert a ditopic supercell graph (superG) into an edge graph (eG).

        - Each V node becomes a NODE in eG (keeps index on superG for later reference).
        - Each unique V-V connection (non-virtual) is represented by:
            * a new EDGE node (named "EDGE_<neg_index>") containing the edge f_points/fcoords
            * a 'real' edge between the two V nodes (keeps the connectivity)
            * two 'half' edges connecting the EDGE node to each V node
        - Virtual edges in superG are preserved (added as edges of type "virtual").
        """
        eG = nx.Graph()

        edge_count = 1  # 2n will be decremented for EDGE naming (EDGE_2, EDGE_4, ...)
        node_count = 0  # 2n+1 will be incremented for NODE naming (NODE_1, NODE_3, ...)
        eG_index_name_dict = {}

        # Use a set of frozensets to record handled undirected node pairs (fast membership test)
        handled_pairs = set()

        # Iterate over nodes once, adding V nodes and creating EDGE nodes for each unique V-V connection
        for n, ndata in superG.nodes(data=True):
            if ndata.get("note") != "V":
                continue

            # Add V node as NODE in eG and store index back to superG for cross-reference

            eG.add_node(
                n,
                f_points=ndata.get("f_points"),
                fcoords=ndata.get("fcoords"),
                type=ndata.get("type", "V"),
                note=ndata.get("note"),
                name="NODE",
                index=2 * node_count + 1,
            )
            superG.nodes[n]["index"] = 2 * node_count + 1
            eG_index_name_dict[2 * node_count + 1] = n  # map index to name

            node_count += 1

            # Process neighbors of this V node
            for ne in superG.neighbors(n):
                pair_key = frozenset((n, ne))

                # Skip if already created for the undirected pair
                if pair_key in handled_pairs:
                    continue

                edge_data = superG.edges[n, ne]
                if edge_data.get("type") == "virtual":
                    # Virtual edges are handled later (preserve original behavior)
                    continue

                # Mark pair as handled
                handled_pairs.add(pair_key)

                # Create EDGE node representing this V-V connection

                edge_name = f"EDGE_{2*edge_count}"
                eG.add_node(
                    edge_name,
                    f_points=edge_data.get("f_points"),
                    fcoords=edge_data.get("fcoords"),
                    type="Edge",
                    name="EDGE",
                    note="E",
                    index=2 * edge_count,  # even indices for EDGEs
                )
                eG_index_name_dict[2 *
                                   edge_count] = edge_name  # map index to name

                # Preserve the connectivity:
                # - a 'real' edge between the two original V nodes
                # - two 'half' edges connecting the EDGE node to each V node
                eG.add_edge(n, ne, index=f"E_{2*edge_count}", type="real")
                eG.add_edge(edge_name,
                            ne,
                            index=f"E_{2*edge_count}",
                            type="half")
                eG.add_edge(edge_name,
                            n,
                            index=f"E_{2*edge_count}",
                            type="half")
                edge_count += 1

        # Preserve virtual edges from superG (add without duplicating)
        for u, v, edata in superG.edges(data=True):
            if edata.get("type") == "virtual":
                if not eG.has_edge(u, v):
                    eG.add_edge(u, v, type="virtual")

        return eG, eG_index_name_dict

    def _make_paired_Xto_x(self, ec_arr, merged_arr, neighbor_number,
                           sc_unit_cell):
        """
        ec_arr: the array of the linker center node
        merged_arr: the array of the merged edges
        neighbor_number: the number of the neighbor nodes of the linker center node
        """
        ec_indices, ec_fpoints = fetch_X_atoms_ind_array(ec_arr, 0, "X")
        if len(ec_indices) < neighbor_number:
            # duplicate the cv_xatoms
            ec_fpoints = np.vstack([ec_fpoints] * neighbor_number)
        # extract only the X atoms in the neighbor edges but not the cv nodes: [len(ec_arr) :]
        nei_indices, nei_fcpoints = fetch_X_atoms_ind_array(merged_arr, 0, "X")

        #skip atom type column
        row_ind, col_ind = find_pair_x_edge_fc(
            ec_fpoints[:, 2:5].astype(float),
            nei_fcpoints[:, 2:5].astype(float), sc_unit_cell)

        # according col_ind order  to reorder the connected edge points
        # switch the X to x for nei_fcpoints
        for i in col_ind:
            if nn(merged_arr[nei_indices[i], 0]) == "X":
                merged_arr[nei_indices[i],
                           0] = "x" + nl(merged_arr[nei_indices[i], 0])
        for k in row_ind:
            if ec_indices[k] < len(ec_arr):
                if nn(ec_arr[ec_indices[k], 0]) == "X":
                    ec_arr[ec_indices[k],
                           0] = "x" + nl(ec_arr[ec_indices[k], 0])

        ordered_edges_points_follow_ecXatoms = order_edge_array(
            row_ind, col_ind, merged_arr)
        # remove the duplicated cv_xatoms
        ec_merged_arr = np.vstack(
            (ec_arr, ordered_edges_points_follow_ecXatoms))
        return ec_merged_arr

    def _find_nearest_neighbor(self, i, n_n_distance_matrix):
        n_n_min_distance = np.min(n_n_distance_matrix[i:i + 1, :])
        _, n_j = locate_min_idx(n_n_distance_matrix[i:i + 1, :])
        # set the column to 1000 to avoid the same atom being selected again
        n_n_distance_matrix[:, n_j] = 1000
        return n_j, n_n_min_distance, n_n_distance_matrix

    def _find_surrounding_points(self, ind, n_n_distance_matrix, max_number):
        stop = 0  # if while loop is too long, stop it
        nearest_neighbor = {}
        nearest_neighbor[ind] = []
        while len(nearest_neighbor[ind]) < max_number:
            stop += 1
            if stop > 100:
                break
            n_j, _, n_n_distance_matrix = self._find_nearest_neighbor(
                ind, n_n_distance_matrix)
            nearest_neighbor[ind].append(n_j)
        return nearest_neighbor

    # Function to find 'XOO' pairs for a specific node
    def _xoo_pair_ind_node(self, single_node_fc, sc_unit_cell):
        # if the node x is not surrounded by two o atoms,
        # then modify the fetch_X_atoms_ind_array(single_node, 0, 'O') find_surrounding_points(k, xs_os_dist_matrix, 2)
        # this function is to find the XOO pairs in a specific node(by node_id),
        # this xoo pair is the indice of x and nearest two o atoms in the same node
        # return the indice of x and nearest two o atoms in the same node, which can be convert to a dict with x_index as key and o_indices as value
        # the distance is in cartesian coordinates
        # single_node_fc: coordinates of any node in the main fragment
        # sc_unit_cell: supercell unit cell matrix
        single_node = np.hstack((
            single_node_fc[:, 0:1],
            np.dot(sc_unit_cell, single_node_fc[:, 2:5].astype(float).T).T,
        ))  # NOTE: modified to skip atom type
        xind, xs_coords = fetch_X_atoms_ind_array(single_node, 0, "X")
        oind, os_coords = fetch_X_atoms_ind_array(single_node, 0, "O")
        xs_os_dist_matrix = np.zeros((len(xs_coords), len(os_coords)))
        for i in range(len(xs_coords)):
            for j in range(len(os_coords)):
                xs_os_dist_matrix[
                    i, j] = np.linalg.norm(xs_coords[i, 1:4].astype(float) -
                                           os_coords[j, 1:4].astype(float))
        xoo_ind_list = []
        for k in range(len(xind)):
            nearest_dict = self._find_surrounding_points(
                k, xs_os_dist_matrix, 2)
            for key in nearest_dict.keys():
                xoo_ind_list.append(
                    [xind[key],
                     sorted([oind[m] for m in nearest_dict[key]])])
        return xoo_ind_list

    def _get_xoo_dict_of_node(self, eG, sc_unit_cell):
        # quick check the order of xoo in every node are same, select n0 and n1, if xoo_ind_node0 == xoo_ind_node1, then xoo_dict is the same
        # return xoo dict of every node, key is x index, value is o index
        n0 = [i for i in eG.nodes() if pname(i) != "EDGE"][0]
        n1 = [i for i in eG.nodes() if pname(i) != "EDGE"][1]
        xoo_ind_node0 = self._xoo_pair_ind_node(
            eG.nodes[n0]["f_points"],
            sc_unit_cell)  # pick node one and get xoo_ind pair
        xoo_ind_node1 = self._xoo_pair_ind_node(
            eG.nodes[n1]["f_points"],
            sc_unit_cell)  # pick node two and get xoo_ind pair
        if xoo_ind_node0 == xoo_ind_node1:
            xoo_dict = {}
            for xoo in xoo_ind_node0:
                xoo_dict[xoo[0]] = xoo[1]
        else:
            self.ostream.print_warning(
                "the order of xoo in every node are not same, please check the input"
            )
            self.ostream.print_info(f"xoo_ind_node0: {xoo_ind_node0}")
            self.ostream.print_info(f"xoo_ind_node1: {xoo_ind_node1}")
            self.ostream.flush()

        return xoo_dict

    def _addxoo2edge_multitopic(self, eG, sc_unit_cell):
        xoo_dict = self._get_xoo_dict_of_node(eG, sc_unit_cell)
        matched_vnode_X = []
        unsaturated_linker = []
        # for every X atom in the EDGE node, search for the paired(nearest) X atom in the connected V node
        # and then use the xoo_dict of the connected V node to extract the xoos of the connected V node
        # and then add the xoos to the EDGE node
        # all xoo_node for the V node is the same
        EDGE_nodes = [n for n in eG.nodes() if pname(n) == "EDGE"]
        for n in EDGE_nodes:
            eG.nodes[n]["xoo_f_points"] = np.zeros((0, 5))
            Xs_edge_indices, Xs_edge_fpoints = fetch_X_atoms_ind_array(
                eG.nodes[n]["f_points"], 0, "X")
            Xs_edge_ccpoints = np.hstack((
                Xs_edge_fpoints[:, 0:2],
                np.dot(sc_unit_cell, Xs_edge_fpoints[:,
                                                     2:5].astype(float).T).T,
            ))
            V_nodes = [i for i in eG.neighbors(n) if pname(i) != "EDGE"]
            if len(V_nodes) == 0:
                # unsaturated_linker.append(n)
                # print(
                #    "no V node connected to this edge node, this linker is a isolated linker, will be ignored",
                #    n,
                # ) # debug
                if self._debug:
                    self.ostream.print_warning(
                        "no V node connected to this edge node, this linker is a isolated linker, will be ignored",
                        n,
                    )
                    self.ostream.flush()
                continue
            all_Xs_vnodes_ind = []
            all_Xs_vnodes_ccpoints = np.zeros((0, 5))
            for v in V_nodes:
                # find the connected V node and its X atoms
                Xs_vnode_indices, Xs_vnode_fpoints = fetch_X_atoms_ind_array(
                    eG.nodes[v]["f_points"], 0, "X")
                Xs_vnode_ccpoints = np.hstack((
                    Xs_vnode_fpoints[:, 0:2],
                    np.dot(sc_unit_cell,
                           Xs_vnode_fpoints[:, 2:5].astype(float).T).T,
                ))
                for ind in Xs_vnode_indices:
                    all_Xs_vnodes_ind.append([v, ind, n])
                all_Xs_vnodes_ccpoints = np.vstack(
                    (all_Xs_vnodes_ccpoints, Xs_vnode_ccpoints))
            edgeX_vnodeX_dist_matrix = np.zeros(
                (len(Xs_edge_ccpoints), len(all_Xs_vnodes_ccpoints)))
            for i in range(len(Xs_edge_ccpoints)):
                for j in range(len(all_Xs_vnodes_ccpoints)):
                    edgeX_vnodeX_dist_matrix[i, j] = np.linalg.norm(
                        Xs_edge_ccpoints[i, 2:5].astype(float) -
                        all_Xs_vnodes_ccpoints[j, 2:5].astype(float))

            for k in range(len(Xs_edge_fpoints)):
                n_j, min_dist, edgeX_vnodeX_dist_matrix = self._find_nearest_neighbor(
                    k, edgeX_vnodeX_dist_matrix)

                if min_dist > 4:
                    unsaturated_linker.append(n)
                    if self._debug:
                        self.ostream.print_warning(
                            f"{min_dist} no xoo for edge node, this linker is a dangling unsaturated linker{n}"
                        )
                        self.ostream.flush()
                    continue
                # add the xoo to the edge node

                nearest_vnode = all_Xs_vnodes_ind[n_j][0]
                nearest_X_ind_in_vnode = all_Xs_vnodes_ind[n_j][1]
                matched_vnode_X.append(all_Xs_vnodes_ind[n_j])
                corresponding_o_indices = xoo_dict[nearest_X_ind_in_vnode]
                xoo_ind_in_vnode = [[nearest_X_ind_in_vnode] +
                                    corresponding_o_indices]
                xoo_fpoints_in_vnode = [
                    eG.nodes[nearest_vnode]["f_points"][i]
                    for i in xoo_ind_in_vnode
                ]
                xoo_fpoints_in_vnode = np.vstack(xoo_fpoints_in_vnode)
                eG.nodes[n]["xoo_f_points"] = np.vstack(
                    (eG.nodes[n]["xoo_f_points"], xoo_fpoints_in_vnode))
                if self._debug:
                    self.ostream.print_info(f"add xoo to edge node {n}")
                    self.ostream.flush()
        return eG, unsaturated_linker, matched_vnode_X, xoo_dict

    def _addxoo2edge_ditopic(self, eG, sc_unit_cell):
        """
        Add XOO groups from V nodes to ditopic EDGE nodes in eG.

        For each EDGE node:
        - find X atoms on the edge and on connected V nodes;
        - compute cartesian coordinates for matching (using sc_unit_cell);
        - for each X on the edge, find the nearest unmatched X on the connected V nodes
          (use _find_nearest_neighbor which also marks matched columns to prevent reuse);
        - if the nearest distance is reasonable (< 4 Å) copy the XOO (X + its two O indices
          from xoo_dict) from the matched V node into the EDGE node's "xoo_f_points";
        - collect unsaturated linkers (no nearby XOO at all sites) and matched_vnode_X .

        Returns:
            eG, unsaturated_linker, matched_vnode_X, xoo_dict
        """
        xoo_dict = self._get_xoo_dict_of_node(eG, sc_unit_cell)
        matched_vnode_X = []
        unsaturated_linker = []

        EDGE_nodes = [n for n in eG.nodes() if pname(n) == "EDGE"]

        for edge in EDGE_nodes:
            # Initialize/clear xoo_f_points for this edge
            eG.nodes[edge]["xoo_f_points"] = np.zeros((0, 5))

            # Fetch X atoms on the edge
            Xs_edge_indices, Xs_edge_fpoints = fetch_X_atoms_ind_array(
                eG.nodes[edge]["f_points"], 0, "X")
            if len(Xs_edge_indices) == 0:
                # nothing to match for this edge
                continue
            elif len(Xs_edge_indices) == 1:
                #only one dot in edge
                #duplicate this dot
                Xs_edge_fpoints = np.vstack([Xs_edge_fpoints] * 2)
                Xs_edge_indices = np.hstack([Xs_edge_indices] * 2)

            # Cartesian coordinates for edge Xs
            Xs_edge_ccpoints = np.hstack((
                Xs_edge_fpoints[:, 0:2],
                np.dot(sc_unit_cell, Xs_edge_fpoints[:,
                                                     2:5].astype(float).T).T,
            ))

            # Collect X atoms from all connected V nodes
            V_nodes = [v for v in eG.neighbors(edge) if pname(v) != "EDGE"]
            if len(V_nodes) == 0:
                if self._debug:
                    self.ostream.print_warning(
                        "no V node connected to this edge node, this linker is a isolated linker, will be ignored",
                        edge,
                    )
                    self.ostream.flush()
                continue

            all_Xs_vnodes_ind = []  # list of [vnode, x_index_in_vnode, edge]
            all_Xs_vnodes_ccpoints = np.zeros((0, 5))
            for v in V_nodes:
                Xs_vnode_indices, Xs_vnode_fpoints = fetch_X_atoms_ind_array(
                    eG.nodes[v]["f_points"], 0, "X")
                if len(Xs_vnode_indices) == 0:
                    continue
                Xs_vnode_ccpoints = np.hstack((
                    Xs_vnode_fpoints[:, 0:2],
                    np.dot(sc_unit_cell,
                           Xs_vnode_fpoints[:, 2:5].astype(float).T).T,
                ))
                for ind in Xs_vnode_indices:
                    all_Xs_vnodes_ind.append([v, ind, edge])
                all_Xs_vnodes_ccpoints = np.vstack(
                    (all_Xs_vnodes_ccpoints, Xs_vnode_ccpoints))

            # If no Xs on connected vnodes, cannot match
            if all_Xs_vnodes_ccpoints.shape[0] == 0:
                if self._debug:
                    self.ostream.print_warning(
                        "no X atoms found on connected V nodes for edge", edge)
                    self.ostream.flush()
                continue

            # Build distance matrix between edge Xs and vnode Xs (cartesian)
            # safe broadcasting-based distance computation
            edge_coords = Xs_edge_ccpoints[:, 2:5].astype(float)
            vnode_coords = all_Xs_vnodes_ccpoints[:, 2:5].astype(float)
            diff = edge_coords[:, None, :] - vnode_coords[None, :, :]
            edgeX_vnodeX_dist_matrix = np.linalg.norm(diff, axis=2)

            # For efficient appending collect xoo fpoints in a list, then stack once
            collected_xoo = []

            # For each edge X find nearest vnode X (while marking matched vnode X columns
            # to avoid reusing the same vnode X multiple times)
            for k in range(len(Xs_edge_fpoints)):
                n_j, min_dist, edgeX_vnodeX_dist_matrix = self._find_nearest_neighbor(
                    k, edgeX_vnodeX_dist_matrix)

                if min_dist > 4.0:
                    # dangling/unsaturated linker
                    unsaturated_linker.append(edge)
                    if self._debug:
                        self.ostream.print_warning(
                            f"{min_dist} no xoo for edge node, this linker is a dangling unsaturated linker{edge}"
                        )
                        self.ostream.flush()
                    continue

                # Add the matched vnode X and its XOO group
                nearest_vnode = all_Xs_vnodes_ind[n_j][0]
                nearest_X_ind_in_vnode = all_Xs_vnodes_ind[n_j][1]
                matched_vnode_X.append(all_Xs_vnodes_ind[n_j])

                # Find corresponding O indices from xoo_dict
                corresponding_o_indices = xoo_dict.get(nearest_X_ind_in_vnode,
                                                       [])
                xoo_inds = [nearest_X_ind_in_vnode] + corresponding_o_indices

                # Extract xoo f_points from the vnode and append
                vnode_fpoints = eG.nodes[nearest_vnode]["f_points"]
                try:
                    xoo_fpoints_in_vnode = vnode_fpoints[xoo_inds]
                except Exception:
                    # Fallback: build by list comprehension if direct indexing fails
                    xoo_fpoints_in_vnode = np.vstack(
                        [vnode_fpoints[idx] for idx in xoo_inds])

                collected_xoo.append(xoo_fpoints_in_vnode)

                if self._debug:
                    self.ostream.print_info(f"add xoo to edge node {edge}")
                    self.ostream.flush()

            # After matching all edge Xs, stack collected xoo points (if any) into node attribute
            if collected_xoo:
                eG.nodes[edge]["xoo_f_points"] = np.vstack(collected_xoo)

        return eG, unsaturated_linker, matched_vnode_X, xoo_dict

    def _make_supercell_range_cleaved_eG(self, eG, supercell, custom_fbox):
        new_eG = eG.copy()
        removed_edges = []
        removed_nodes = []

        if custom_fbox is None:
            #self.cleaved_eG = new_eG.copy()
            #self.cleaved_edges = removed_edges
            #self.cleaved_nodes = removed_nodes
            #return
            custom_fbox = [[
                0,
                supercell[0],
            ], [0, supercell[1]], [0, supercell[2]]]
        custom_fbox = np.array(custom_fbox)  # [[xlo,xhi],[ylo,yhi],[zlo,zhi]]

        def check_supercell_box_range(fcoords, supercell, custom_fbox):
            # check if the fcoords is in the supercell box range with buffer
            for i in range(3):
                if fcoords[i] < custom_fbox[i,
                                            0] or fcoords[i] > custom_fbox[i,
                                                                           1]:
                    return False
            return True

        for n in eG.nodes():
            if pname(n) != "EDGE":
                if check_supercell_box_range(eG.nodes[n]["fcoords"], supercell,
                                             custom_fbox):
                    pass
                else:
                    new_eG.remove_node(n)
                    removed_nodes.append(n)
            elif pname(n) == "EDGE":
                # ditopic linker have two points in the fcoords
                if (arr_dimension(eG.nodes[n]["fcoords"]) == 2):
                    edge_coords = np.mean(eG.nodes[n]["fcoords"], axis=0)
                # multitopic linker have one point in the fcoords from EC
                elif (arr_dimension(eG.nodes[n]["fcoords"]) == 1):
                    edge_coords = eG.nodes[n]["fcoords"]

                if check_supercell_box_range(edge_coords, supercell,
                                             custom_fbox):
                    pass
                else:
                    new_eG.remove_node(n)
                    removed_edges.append(n)

        self.matched_vnode_xind = self._update_matched_nodes_xind(
            removed_nodes, removed_edges, self.matched_vnode_xind)

        self.cleaved_eG = new_eG.copy()
        self.cleaved_edges = removed_edges
        self.cleaved_nodes = removed_nodes

    def _update_matched_nodes_xind(self, removed_nodes_list,
                                   removed_edges_list, matched_vnode_xind):
        """
        Update the matched_vnode_xind list after nodes/edges have been removed by cleaving.

        Each entry in matched_vnode_xind is expected to be a tuple/list of the form:
            (node, xind, edge)
        - node: the vnode identifier (string) that was matched
        - xind: the index (or indices) of the X atom(s) in that node
        - edge: the edge identifier (string) that was matched

        This function removes entries where:
        - the matched edge has been removed (edge in removed_edges_list) but the vnode itself remains,
          because that match is no longer valid; or
        - the vnode itself has been removed (node in removed_nodes_list).

        Parameters:
            removed_nodes_list (list): node identifiers removed during cleaving
            removed_edges_list (list): edge identifiers removed during cleaving
            matched_vnode_xind (list): list of (node, xind, edge) matches to be filtered

        Returns:
            list: filtered matched_vnode_xind with invalid entries removed
        """
        # indices of entries to delete from matched_vnode_xind
        to_remove_row = []

        # Iterate through the matched list and mark entries that should be removed.
        # We use indices because we will build a new list by excluding these indices.
        for i in range(len(matched_vnode_xind)):
            node, xind, edge = matched_vnode_xind[i]

            # If the edge was removed but the node is still present, the match is invalid.
            if edge in removed_edges_list and node not in removed_nodes_list:
                if self._debug:
                    self.ostream.print_info(
                        f"remove edge {edge}, {matched_vnode_xind[i]} from matched_vnode_xind"
                    )
                to_remove_row.append(i)

            # If the node itself was removed, remove the match regardless of the edge.
            elif node in removed_nodes_list:
                if self._debug:
                    self.ostream.print_info(
                        f"remove node {node}, {matched_vnode_xind[i]} from matched_vnode_xind"
                    )
                to_remove_row.append(i)

        # Build a new list excluding the marked indices. This preserves original ordering
        # for entries that remain valid.
        update_matched_vnode_xind = [
            entry for j, entry in enumerate(matched_vnode_xind)
            if j not in to_remove_row
        ]

        return update_matched_vnode_xind


def remove_node_by_index(eG, remove_node_list, remove_edge_list):
    """Remove from eG all nodes whose index is in remove_node_list and all EDGEs whose index is in remove_edge_list (by -index). Modifies eG in place."""
    for n in eG.nodes():
        if pname(n) != "EDGE":
            if eG.nodes[n]["index"] in remove_node_list:
                eG.remove_node(n)
        if pname(n) == "EDGE":
            if -1 * eG.nodes[n]["index"] in remove_edge_list:
                eG.remove_node(n)
    return eG
