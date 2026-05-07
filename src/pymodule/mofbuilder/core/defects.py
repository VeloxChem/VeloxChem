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
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import networkx as nx
from mpi4py import MPI
import h5py
import re

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
from .superimpose import superimpose_rotation_only, superimpose


class TerminationDefectGenerator:
    """Generate defects by finding unsaturated nodes, removing nodes/linkers, and adding terminations.

    Can remove specified nodes/linkers (by index), optionally clean unsaturated linkers,
    add termination groups to unsaturated nodes (using termination_data), and replace
    nodes or linkers by superimposing new geometry onto the existing graph.

    Attributes:
        comm: MPI communicator.
        rank: MPI rank of this process.
        nodes: MPI size (number of processes).
        ostream: Output stream for logging.
        cleaved_eG: Edge graph (eG) after cleaving; set before use.
        node_connectivity: Expected number of real neighbors per node.
        linker_connectivity: Expected number of connections per linker (edge).
        eG_index_name_dict: Mapping from eG index to node/edge name.
        sc_unit_cell: 3x3 supercell unit cell matrix.
        sc_unit_cell_inv: Inverse of sc_unit_cell.
        res_idx2rm: List of eG indices to remove (nodes/edges).
        xoo_dict: Dict mapping X index to list of O indices (XOO groups) in node.
        matched_vnode_xind: List of (node, xind, edge) for matched V-node X to edges.
        updated_matched_vnode_xind: Same after removal/termination.
        unsaturated_nodes: List of unsaturated node names.
        updated_unsaturated_nodes: Same after removal.
        unsaturated_linkers: List of unsaturated linker (EDGE) names.
        termination_data: Full termination atom array (XOO).
        termination_X_data: Termination X atoms.
        termination_Y_data: Termination Y atoms (e.g. O-O center).
        use_termination: If True, add terminations to unsaturated nodes.
        clean_unsaturated_linkers: If True, remove unsaturated linkers from graph.
        update_node_termination: If True, update matched nodes after removal.
        nodes_idx2rp: Indices of nodes to replace.
        linkers_idx2rp: Indices of linkers to replace.
        new_node_data: New node atom data for replacement (set by caller).
        new_node_X_data: New node X-atom data for replacement.
        new_linker_data: New linker atom data for replacement.
        new_linker_X_data: New linker X-atom data for replacement.
        defectG: eG after removing nodes/linkers (no terminations).
        termG: eG after adding terminations.
        finalG: eG after removing XOO from nodes.
        _debug: If True, print extra debug messages.
    """

    def __init__(
        self,
        comm: Optional[Any] = None,
        ostream: Optional[Any] = None,
    ) -> None:
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)

        #need to be set before use
        self.cleaved_eG = None
        self.node_connectivity = None
        self.linker_connectivity = None
        self.eG_index_name_dict = None
        self.sc_unit_cell = None
        self.sc_unit_cell_inv = None
        #for remove
        self.res_idx2rm = []
        self.xoo_dict = None

        self.matched_vnode_xind = None
        self.updated_matched_vnode_xind = None  #will be set after remove
        self.unsaturated_nodes = []  # list of unsaturated node name
        self.updated_unsaturated_nodes = []  #will be set after remove
        self.unsaturated_linkers = []  # list of unsaturated linker name
        self.termination_data = []  # list of termination XOO vectors
        self.termination_X_data = []  # list of termination X atoms
        self.termination_Y_data = []  # list of termination Y atoms

        self.use_termination = True

        self.clean_unsaturated_linkers = False
        self.update_node_termination = True

        #for replace
        self.nodes_idx2rp = []
        self.linkers_idx2rp = []
        self.new_node_data = None
        self.new_node_X_data = None
        self.new_linker_data = None
        self.new_linker_X_data = None

        #will be set after use
        #self.saved_eG_matched_vnode_xind = None  #saved matched_vnode_xind before remove items
        self.defectG = None  #eG after removing nodes or linkers
        self.termG = None  #eG after adding terminations
        self.finalG = None  #eG after removing xoo from node

        #debug
        self._debug = False

    def remove_items_or_terminate(
        self,
        res_idx2rm: Optional[List[int]] = None,
        cleaved_eG: Optional[nx.Graph] = None,
    ) -> nx.Graph:
        """Remove nodes/linkers by index and optionally add terminations to unsaturated nodes.

        If res_idx2rm is empty, only finds unsaturated linkers and optionally adds
        terminations to unsaturated nodes. Otherwise removes the given nodes/linkers,
        optionally removes unsaturated linkers, then adds terminations to unsaturated
        nodes using termination_data.

        Args:
            res_idx2rm: List of eG indices (nodes/edges) to remove.
            cleaved_eG: Edge graph (eG) to modify; a copy is used.

        Returns:
            Defect graph after removal and termination (termG or defectG).
        """
        if res_idx2rm is None:
            res_idx2rm = []
        #just terminate nodes
        if not res_idx2rm:
            defectG = cleaved_eG.copy()
            new_unsaturated_linkers = self._find_unsaturated_linkers(
                defectG, self.linker_connectivity)
            if (not self.use_termination) and (
                    not self.clean_unsaturated_linkers):
                self.unsaturated_linkers = new_unsaturated_linkers
                return defectG
            elif (not self.use_termination) and (
                    self.clean_unsaturated_linkers):
                for linker_name in new_unsaturated_linkers:
                    if linker_name in defectG.nodes():
                        defectG.remove_node(linker_name)
                        self.ostream.print_info(
                            f"unsaturated linker {linker_name} removed")
                self.updated_unsaturated_nodes = self.unsaturated_nodes
                self.updated_matched_vnode_xind = self.matched_vnode_xind
                return defectG
            elif (self.use_termination) and (
                    not self.clean_unsaturated_linkers):
                new_unsaturated_nodes = self._find_unsaturated_nodes(
                    defectG, self.node_connectivity)
                termG, _ = self._add_terminations_to_unsaturated_nodes(
                    defectG, new_unsaturated_nodes, self.matched_vnode_xind)
                self.updated_unsaturated_nodes = new_unsaturated_nodes
                self.updated_matched_vnode_xind = self.matched_vnode_xind
                self.unsaturated_linkers = new_unsaturated_linkers
                return termG
            else:  #self.use_termination and self.clean_unsaturated_linkers
                res_idx2rm = []
                for k, v in self.eG_index_name_dict.items():
                    if v in new_unsaturated_linkers:
                        res_idx2rm.append(k)

        self.ostream.print_info(f"trying removing items: {res_idx2rm}")
        nodes_names2rm = self._extract_node_name_from_eG_dict(
            res_idx2rm, self.eG_index_name_dict)
        self.ostream.print_info(f"removing nodes/linkers: {nodes_names2rm}")
        self.ostream.flush()

        defectG = cleaved_eG.copy()

        new_unsaturated_nodes = self._find_unsaturated_nodes(
            defectG, self.node_connectivity)
        for node_name in nodes_names2rm:
            if node_name in defectG.nodes():
                defectG.remove_node(node_name)
                self.ostream.print_info(f"node {node_name} removed")
                self.ostream.flush()

        #remove all unsaturated linkers
        if self.clean_unsaturated_linkers:
            new_unsaturated_linkers = self._find_unsaturated_linkers(
                defectG, self.linker_connectivity)
            for linker_name in new_unsaturated_linkers:
                if linker_name in defectG.nodes():
                    defectG.remove_node(linker_name)
                    self.ostream.print_info(
                        f"unsaturated linker {linker_name} removed")
                    self.ostream.flush()
            self.unsaturated_linkers = new_unsaturated_linkers

        if self._debug:
            self.ostream.print_info(
                f"new unsaturated nodes: {new_unsaturated_nodes}")
            self.ostream.print_info(
                f"new unsaturated linkers: {new_unsaturated_linkers}")
            self.ostream.flush()

        #old_unsaturated_nodes = self.unsaturated_nodes
        #old_unsaturated_linkers = self.unsaturated_linkers

        if not self.use_termination:
            if self._debug:
                self.ostream.print_info("no termination, return the defectG")
                self.ostream.flush()
            return defectG

        if self.update_node_termination:
            #update unsaturated nodes
            self.ostream.print_info("update unsaturated nodes")
            self.ostream.flush()

            updated_matched_vnode_xind = self._update_matched_nodes_xind(
                nodes_names2rm, self.matched_vnode_xind)
            #self.saved_eG_matched_vnode_xind = self.matched_vnode_xind
            #add termination to the new unsaturated node
            termG, _ = self._add_terminations_to_unsaturated_nodes(
                defectG, new_unsaturated_nodes, updated_matched_vnode_xind)
            self.updated_unsaturated_nodes = new_unsaturated_nodes
            self.updated_matched_vnode_xind = updated_matched_vnode_xind
            return termG

        else:
            #self.matched_vnode_xind = self.saved_eG_matched_vnode_xind
            #add termination to the old unsaturated node
            termG, _ = self._add_terminations_to_unsaturated_nodes(
                defectG, self.unsaturated_nodes, self.matched_vnode_xind)
            self.updated_unsaturated_nodes = self.unsaturated_nodes
            self.updated_matched_vnode_xind = self.matched_vnode_xind
            return termG

    def replace_items(
        self, res_idx2rp: List[int], G: nx.Graph
    ) -> nx.Graph:
        """Replace specified nodes or linkers in G with new_node_data/new_linker_data by superimposition.

        Args:
            res_idx2rp: List of eG indices to replace (nodes or EDGE names).
            G: Edge graph to modify (copy is used).

        Returns:
            Graph with replaced node/linker coordinates; unchanged if data not set.
        """
        nodes_name2rp = self._extract_node_name_from_eG_dict(
            res_idx2rp, self.eG_index_name_dict)
        #split the node or edge name
        new_nodes_name = [n for n in nodes_name2rp if pname(n) != "EDGE"]
        replace_edges_name = [n for n in nodes_name2rp if pname(n) == "EDGE"]
        if (not new_nodes_name) and (not replace_edges_name):
            if self._debug:
                self.ostream.print_info(
                    "no nodes or linkers to replace, return the original G")
                self.ostream.flush()
            return G
        rpG = G.copy()
        if new_nodes_name:
            self.ostream.print_info(f"replace nodes: {new_nodes_name}")
            #check if the new_node_data is set
            if self.new_node_data is None:
                self.ostream.print_warning(
                    "new_node_data is not set, skip replacing nodes")
                #skip replacing nodes
                rpG = G.copy()
            else:
                rpG = self._replace_items_in_G(new_nodes_name, G,
                                               self.new_node_data,
                                               self.new_node_X_data,
                                               self.sc_unit_cell_inv)
        if replace_edges_name:
            self.ostream.print_info(f"replace linkers: {replace_edges_name}")
            #check if the new_linker_data is set
            if self.new_linker_data is None:
                self.ostream.print_warning(
                    "new_linker_data is not set, skip replacing linkers")
                #skip replacing linkers
                rpG = rpG.copy()
            else:
                rpG = self._replace_items_in_G(replace_edges_name, rpG,
                                               self.new_linker_data,
                                               self.new_linker_X_data,
                                               self.sc_unit_cell_inv)
        return rpG

    def _find_unsaturated_nodes(
        self, eG: nx.Graph, node_connectivity: int
    ) -> List[str]:
        """Return list of node names in eG that have fewer real neighbors than node_connectivity."""
        unsaturated_nodes = []
        for n in eG.nodes():
            if pname(n) != "EDGE":
                real_neighbor = []
                for cn in eG.neighbors(n):
                    if eG.edges[(n, cn)]["type"] == "real":
                        real_neighbor.append(cn)
                if len(real_neighbor) < node_connectivity:
                    unsaturated_nodes.append(n)
        return unsaturated_nodes

    def _find_unsaturated_linkers(
        self, eG: nx.Graph, linker_topics: int
    ) -> List[str]:
        """Return list of EDGE node names in eG that have fewer neighbors than linker_topics."""
        unsaturated_linkers = []
        for n in eG.nodes():
            if pname(n) == "EDGE" and len(list(
                    eG.neighbors(n))) < linker_topics:
                unsaturated_linkers.append(n)
        return unsaturated_linkers

    def _extract_node_name_from_eG_dict(
        self, idx_list: List[int], eG_dict: Dict[int, str]
    ) -> List[str]:
        """Map eG indices to node/edge names using eG_index_name_dict."""
        return [eG_dict[i] for i in idx_list if i in eG_dict]

    def _update_matched_nodes_xind(
        self,
        nodes_name_list: List[str],
        old_matched_vnode_xind: List[Tuple[str, int, str]],
    ) -> List[Tuple[str, int, str]]:
        """Drop (node, xind, edge) entries where node or edge was removed from the graph."""
        update_matched_vnode_xind = []
        for i in range(len(old_matched_vnode_xind)):
            # format: old_matched_vnode_xind[i]: old_matched_vnode_xind[i]
            if (old_matched_vnode_xind[i][0]
                    in nodes_name_list) or (old_matched_vnode_xind[i][-1]
                                            in nodes_name_list):
                if self._debug:
                    self.ostream.print_info(
                        f"matched_vnode_xind {old_matched_vnode_xind[i]} removed"
                    )
                continue
            update_matched_vnode_xind.append(old_matched_vnode_xind[i])

        return update_matched_vnode_xind

    def _add_terminations_to_unsaturated_nodes(self,
                                               eG,
                                               unsaturated_nodes: Optional[Sequence[Any]] = None,
                                               matched_vnode_xind: Optional[Sequence[Any]] = None):
        """
        use the node terminations to add terminations to the unsaturated nodes
        """
        unsaturated_nodes = [] if unsaturated_nodes is None else unsaturated_nodes
        matched_vnode_xind = [] if matched_vnode_xind is None else matched_vnode_xind

        #generate term_xoovecs
        term_xoovecs = np.vstack((self.termination_X_data,
                                  self.termination_Y_data))[:,
                                                            5:8].astype(float)
        term_coords = self.termination_data[:, 5:8].astype(float)
        term_atoms = self.termination_data[:, 0:2]

        termG = eG.copy()
        unsaturated_nodes = [
            n for n in unsaturated_nodes if n in termG.nodes()
        ]

        if not unsaturated_nodes:
            self.ostream.print_warning(
                "no unsaturated nodes found in the graph, skip adding terminations"
            )
            return termG, {}

        (
            unsaturated_vnode_xind_dict,
            unsaturated_vnode_xoo_dict,  #will be updated after adding terminations
            matched_vnode_xind_dict,
        ) = self._make_unsaturated_vnode_xoo_dict(unsaturated_nodes,
                                                  self.xoo_dict,
                                                  matched_vnode_xind, termG,
                                                  self.sc_unit_cell)
        # term_file: path to the termination file
        # ex_node_cxo_cc: exposed node coordinates

        node_oovecs_record = []
        for n in termG.nodes():
            termG.nodes[n]["term_c_points"] = {}
        for exvnode_xind_key in unsaturated_vnode_xoo_dict.keys():
            exvnode_x_ccoords = unsaturated_vnode_xoo_dict[exvnode_xind_key][
                "x_cpoints"]
            exvnode_oo_ccoords = unsaturated_vnode_xoo_dict[exvnode_xind_key][
                "oo_cpoints"]
            node_xoo_ccoords = np.vstack(
                [exvnode_x_ccoords, exvnode_oo_ccoords])

            # make the beginning point of the termination to the center of the oo atoms
            node_oo_center_cvec = np.mean(
                exvnode_oo_ccoords[:, 2:5].astype(float), axis=0)
            node_xoo_cvecs = (node_xoo_ccoords[:, 2:5].astype(float) -
                              node_oo_center_cvec)
            node_xoo_cvecs = node_xoo_cvecs.astype("float")
            # use record to record the rotation matrix for get rid of the repeat calculation
            indices = [
                index for index, value in enumerate(node_oovecs_record)
                if is_list_A_in_B(node_xoo_cvecs, value[0])
            ]
            if len(indices) == 1:
                rot = node_oovecs_record[indices[0]][1]
            else:
                _, rot, _ = superimpose(term_xoovecs, node_xoo_cvecs)
                node_oovecs_record.append((node_xoo_cvecs, rot))

            adjusted_term_vecs = np.dot(term_coords, rot) + node_oo_center_cvec
            adjusted_term = np.hstack((
                term_atoms,
                adjusted_term_vecs,
            ))
            # add the adjusted term to the terms, add index, add the node name
            unsaturated_vnode_xoo_dict[exvnode_xind_key][
                "node_term_c_points"] = (adjusted_term)
            termG.nodes[exvnode_xind_key[0]]["term_c_points"][
                exvnode_xind_key[1]] = (adjusted_term)
            if self._debug:
                self.ostream.print_info(
                    f"node {exvnode_xind_key[0]} xind {exvnode_xind_key[1]} added termination"
                )
                self.ostream.flush()

        return termG, unsaturated_vnode_xoo_dict

    def _make_unsaturated_vnode_xoo_dict(self, unsaturated_nodes, xoo_dict,
                                         matched_vnode_xind, eG, sc_unit_cell):
        """Build data structures for unsaturated nodes:
              - unsaturated_vnode_xind_dict: for each unsaturated node -> list of exposed X indices
              - unsaturated_vnode_xoo_dict: for each (node, xind) -> metadata including fractional and cartesian points
              - matched_vnode_xind_dict: mapping node -> list of matched x indices
            """
        # Normalize matched_vnode_xind into a dict: node -> [xind, ...]
        matched_vnode_xind = matched_vnode_xind or []
        matched_vnode_xind_dict = {}
        for entry in matched_vnode_xind:
            try:
                n, xind, _ = entry  # entry expected as [node, xind, edge] (edge ignored here)
            except Exception:
                # malformed entry: skip
                continue
            matched_vnode_xind_dict.setdefault(n, []).append(xind)

        # Prepare exposed x-index list for each unsaturated vnode
        xoo_keys = list(xoo_dict.keys())
        unsaturated_vnode_xind_dict = {}
        for vnode in unsaturated_nodes:
            if vnode not in eG.nodes():
                # node not present in graph (may have been removed) -> treat as no exposed indices
                unsaturated_vnode_xind_dict[vnode] = []
                continue

            if vnode in matched_vnode_xind_dict:
                # exposed indices are those X keys not already matched for this node
                matched_list = matched_vnode_xind_dict[vnode]
                unsaturated_vnode_xind_dict[vnode] = [
                    i for i in xoo_keys if i not in matched_list
                ]
            else:
                # no matched indices -> all X keys are exposed
                # copy to avoid accidental mutation of xoo_keys
                unsaturated_vnode_xind_dict[vnode] = xoo_keys.copy()
                self.ostream.print_warning(
                    f"Node is isolated {vnode} whose index is {eG.nodes[vnode]['index']}"
                )

        # Build the detailed X/OO dictionary for each exposed (node, xind)
        unsaturated_vnode_xoo_dict = {}
        for vnode, exposed_x_indices in unsaturated_vnode_xind_dict.items():
            # skip nodes without exposed indices
            if not exposed_x_indices:
                continue

            fpoints = eG.nodes[vnode].get("f_points")
            if fpoints is None:
                # missing data: skip this node
                continue

            try:
                fpoints_len = len(fpoints)
            except Exception:
                fpoints_len = 0

            for xind in exposed_x_indices:
                # validate x index
                if not (0 <= xind < fpoints_len):
                    # invalid index -> skip
                    continue

                # fractional f_points for this X
                x_fpoints = fpoints[xind]
                # build Cartesian coords for X
                x_cpoints = np.hstack((
                    x_fpoints[0:2],
                    fractional_to_cartesian(x_fpoints[2:5], sc_unit_cell),
                ))

                # get associated OO indices (may be missing -> empty list)
                oo_ind_in_vnode = xoo_dict.get(xind, [])
                if oo_ind_in_vnode:
                    # gather fractional arrays for those OO atoms and compute Cartesian coords in one call
                    try:
                        oo_fpoints_in_vnode = np.vstack(
                            [fpoints[i] for i in oo_ind_in_vnode])
                    except Exception:
                        # if any index invalid, skip this xind
                        continue
                    oo_cpoints = np.hstack((
                        oo_fpoints_in_vnode[:, 0:2],
                        fractional_to_cartesian(oo_fpoints_in_vnode[:, 2:5],
                                                sc_unit_cell),
                    ))
                else:
                    # no OO partners found
                    oo_fpoints_in_vnode = np.empty((0, 0))
                    oo_cpoints = np.empty((0, 3))

                unsaturated_vnode_xoo_dict[(vnode, xind)] = {
                    "xind": xind,
                    "oo_ind": oo_ind_in_vnode,
                    "x_fpoints": x_fpoints,
                    "x_cpoints": x_cpoints,
                    "oo_fpoints": oo_fpoints_in_vnode,
                    "oo_cpoints": oo_cpoints,
                }

        return (
            unsaturated_vnode_xind_dict,
            unsaturated_vnode_xoo_dict,
            matched_vnode_xind_dict,
        )

    def _remove_xoo_from_node(self, G: nx.Graph) -> nx.Graph:
        """Remove XOO atom rows from each node's f_points and set noxoo_f_points on each node in eG."""

        xoo_dict = self.xoo_dict
        eG = G.copy()

        all_xoo_indices = []
        for x_ind, oo_ind in xoo_dict.items():
            all_xoo_indices.append(x_ind)
            all_xoo_indices.extend(oo_ind)

        for n in eG.nodes():
            if pname(n) != "EDGE":
                all_f_points = eG.nodes[n]["f_points"]
                noxoo_f_points = np.delete(all_f_points,
                                           all_xoo_indices,
                                           axis=0)
                eG.nodes[n]["noxoo_f_points"] = noxoo_f_points

        return eG

    def _replace_items_in_G(
        self,
        edge_n_list: List[str],
        G: nx.Graph,
        new_n_data: np.ndarray,
        new_n_X_data: np.ndarray,
        sc_unit_cell_inv: np.ndarray,
        newname: Optional[str] = None,
    ) -> nx.Graph:
        """Replace coordinates of nodes in edge_n_list with new_n_data, aligned by X atoms via superimpose."""
        new_n_atoms = new_n_data[:, 0:2]
        new_n_ccoords = new_n_data[:, 5:8].astype(float)
        new_n_x_ccoords = new_n_X_data[:, 5:8].astype(float)

        for n in edge_n_list:
            if n not in G.nodes():
                self.ostream.print_warning(
                    f"this linker is not in MOF, will be skipped {n}")
                continue
            n_f_points = G.nodes[n]["f_points"]
            x_indices = [
                i for i in range(len(n_f_points))
                if nn(n_f_points[i][0]) == "X"
            ]
            #if x atoms in new linker is not equal to the x atoms in the processed edge print warning and skip
            if len(x_indices) != len(new_n_x_ccoords):
                self.ostream.print_warning(
                    f"the number of X atoms in the new fragment does not match the origin connectivity of {n}, skipping"
                )
                continue
            old_n_x_points = n_f_points[x_indices]
            old_n_com = np.mean(old_n_x_points[:, 2:5].astype(float), axis=0)
            old_n_x_fcoords = old_n_x_points[:, 2:5].astype(float) - old_n_com

            new_n_x_fcoords = cartesian_to_fractional(new_n_x_ccoords,
                                                      sc_unit_cell_inv)
            new_n_fcoords = cartesian_to_fractional(new_n_ccoords,
                                                    sc_unit_cell_inv)

            _, rot, trans = superimpose(new_n_x_fcoords, old_n_x_fcoords)
            replaced_linker_fcoords = np.dot(new_n_fcoords, rot) + old_n_com
            replaced_linker_f_points = np.hstack(
                (new_n_atoms, replaced_linker_fcoords))

            G.nodes[n]["f_points"] = replaced_linker_f_points
            G.nodes[n][
                "name"] = newname if newname else "R" + G.nodes[n]["name"]
            self.ostream.print_info(f"linker {n} replaced")

        return G
