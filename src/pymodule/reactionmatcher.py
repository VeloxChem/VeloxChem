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

from mpi4py import MPI
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.algorithms.isomorphism import categorical_node_match
import networkx as nx
import numpy as np
import time
import sys

from .outputstream import OutputStream
from .veloxchemlib import mpi_master


class ReactionMatcher:

    def __init__(self, comm=None, ostream=None):
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # output stream
        self.ostream = ostream

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self._iso_count = 0
        self._mono_count = 0
        self._start_time = 0
        self.max_time = 600
        self._breaking_depth = 2  # how many breaking edges to try
        self._check_monomorphic = False

    def get_mapping(self, reactant_ff, rea_elems, product_ff, pro_elems,
                    breaking_bonds):

        self.rea_graph, self.pro_graph = self._create_reaction_graphs(
            reactant_ff,
            rea_elems,
            product_ff,
            pro_elems,
            breaking_bonds,
        )

        map, breaking_edges, forming_edges = self._find_mapping(
            self.rea_graph.copy(), self.pro_graph.copy(), breaking_bonds)
        if map is None:

            self._check_monomorphic = True
            self.ostream.print_info(
                "No mapping found by checking subgraph isomorphism, trying to find mapping with subgraph monomorphism."
            )
            self.ostream.print_info(f"Total subgraph isomorphism checks: {self._iso_count}")
            self.ostream.flush()
            self._iso_count = 0
            map, breaking_edges, forming_edges = self._find_mapping(
                self.rea_graph.copy(), self.pro_graph.copy(), breaking_bonds)
        if map is None:
            self.ostream.print_warning(
                "No subgraph solution with two broken bonds found, aborting. Try suggesting more breaking bonds."
            )
        self.ostream.flush()
        self.ostream.print_info(f"Total subgraph isomorphism checks: {self._iso_count}, total subgraph monomorphism checks: {self._mono_count}")
        return map, breaking_edges, forming_edges

    def _create_reaction_graphs(self, reactant_ff, rea_elems, product_ff,
                                pro_elems, breaking_bonds):
        rea_graph = nx.Graph()
        reactant_bonds = list(reactant_ff.bonds.keys())
        # Remove the bonds that are being broken, so that these segments get treated as seperate reactants

        rea_graph.add_nodes_from(reactant_ff.atoms.keys())
        rea_graph.add_edges_from(reactant_bonds)

        for i, elem in enumerate(rea_elems):
            rea_graph.nodes[i]['elem'] = elem

        pro_graph = nx.Graph()
        pro_graph.add_nodes_from(product_ff.atoms.keys())
        pro_graph.add_edges_from(list(product_ff.bonds.keys()))
        for i, elem in enumerate(pro_elems):
            pro_graph.nodes[i]['elem'] = elem

        return rea_graph, pro_graph

    def _find_mapping(self, A, B, forced_breaking_edges=set()):
        """
        Find a mapping between the connected components of A and B, while avoiding reconnecting broken edges.
        """
        # make copies of A and B to avoid modifying the original graphs
        for edge in forced_breaking_edges:
            A.remove_edge(*edge)
        # loop progressively over more and more broken bonds till all connected components of A are subgraphs of B
        # can be done based on elements
        swapped = False

        if len(A.edges()) > len(B.edges()):
            self.ostream.print_info(
                "A (reactant) has more edges than B (product), swapping A and B."
            )
            self.ostream.flush()
            A, B = B, A
            swapped = True

        self._start_time = time.time()
        breaking_edges = None
        for depth in range(1, self._breaking_depth):
            self.ostream.print_info(
                f"Finding breaking edges at depth {depth}"
            )
            self.ostream.flush()
            breaking_edges = self._find_breaking_edges(A, B, depth)
            if breaking_edges is not None:
                break
            
        if breaking_edges is None:
            return None, None, None

        insert = "monomorphisms" if self._check_monomorphic else "isomorphisms"
        self.ostream.print_info(
            f"Found subgraph {insert} with breaking bonds: {list(breaking_edges)} and forced breaking bonds: {list(forced_breaking_edges)}, continuing to find forming bonds"
        )
        self.ostream.flush()
        forming_edges = self._find_forming_edges(
            A, B, breaking_edges | forced_breaking_edges)
        if forming_edges is None:
            return None, None, None

        self.ostream.print_info(
            f"Found forming bonds: {forming_edges}. Finding mapping from isomorphism."
        )
        self.ostream.flush()
        # print(forming_edges)
        GM = GraphMatcher(A, B, categorical_node_match('elem', ''))
        map = next(GM.isomorphisms_iter())

        self._check_time(
            f"finding mapping.")
        if swapped:
            # if we swapped A and B, we need to swap the mapping back
            self.ostream.print_info(
                "Inverting mapping because A and B were swapped.")
            map = {v: k for k, v in map.items()}
        return map, breaking_edges, forming_edges


    def _find_breaking_edges(self,A,B, depth):
        if not self._check_time():
            return None

        if self._connected_components_are_subgraphs(A, B):
            return set()
            
        if depth > 0:
            for edge in A.edges():
                A.remove_edge(*edge)
                # self.ostream.print_info(
                #     f"Trying to find breaking edge {edge} at depth {depth}"
                # )
                self.ostream.flush()
                edges = self._find_breaking_edges(A, B, depth-1)
                if edges is not None:
                    edges.add(edge)
                    self.ostream.print_info(
                        f"Found breaking edges: {edges} at depth {depth}."
                    )
                    self.ostream.flush()
                    return edges
                A.add_edge(*edge)
        return None


    def _find_forming_edges(self, A, B, breaking_edges):
        """
        Find forming edges in A that are not in B, while avoiding reconnecting broken edges.
        """
        bond_count = 0
        forming_edges = []
        nodes = list(A.nodes())
        active_nodes = set()
        for edge in breaking_edges:
            active_nodes.add(edge[0])
            active_nodes.add(edge[1])
        # Move every element of active_nodes in nodes to the front of the list
        nodes = [n for n in nodes if n in active_nodes] + [n for n in nodes if n not in active_nodes]

        elem_combinations_B = {}
        for edge in B.edges():
            comb = self._get_elem_comb(edge[0], edge[1], B)
            if comb not in elem_combinations_B.keys():
                elem_combinations_B[comb] = 1
            else:
                elem_combinations_B[comb] += 1

        elem_combinations_A = {}
        for edge in A.edges():
            comb = self._get_elem_comb(edge[0], edge[1], A)
            if comb not in elem_combinations_A.keys():
                elem_combinations_A[comb] = 1
            else:
                elem_combinations_A[comb] += 1

        #todo only look through the list of node combinations that form bonds of element pairs that are missing
        while len(A.edges()) < len(B.edges()):
            bond_count += 1
            edge = None
            for i,node_i in enumerate(nodes):
                for j,node_j in enumerate(A.nodes()):
                    if j<= i:
                        continue
                    edge_in_A = (node_i, node_j) in A.edges() or (node_j, node_i) in A.edges()
                    edge_in_broken = (node_i, node_j) in breaking_edges or (
                        node_j, node_i) in breaking_edges
                    
                    # HH_bond = A.nodes[node_i]['elem'] == 1.0 and A.nodes[node_j]['elem'] == 1.0
                    if edge_in_A or edge_in_broken:
                        continue
                    
                    comb = self._get_elem_comb(node_i, node_j, A)
                    if elem_combinations_B.get(comb, 0) <= elem_combinations_A.get(comb, 0):
                        continue

                    A.add_edge(node_i, node_j)
                    self.ostream.print_info(
                        f"Trying to find bond {bond_count}: ({node_i}, {node_j})"
                    )
                    self.ostream.flush()
                    if not self._connected_components_are_subgraphs(A, B):
                        A.remove_edge(node_i, node_j)
                        continue
                    if not self._check_time():
                        return None

                    edge = (node_i, node_j)
                    forming_edges.append(edge)
                    self._check_time(f"finding bond {bond_count}: {edge}")
                    self.ostream.flush()
                    break
                if edge is not None:
                    break
            if edge is None:
                self.ostream.print_info(f"Bond not found, aborting")
                self.ostream.flush()
                return None
        return forming_edges
    
    @staticmethod
    def _get_elem_comb(node1,node2,A):
        """
        Get the element combination of two nodes in a graph.
        """
        elem1 = A.nodes[node1]['elem']
        elem2 = A.nodes[node2]['elem']
        return tuple(sorted((elem1, elem2)))

    
    def _connected_components_are_subgraphs(self, A, B):
        """
        Check if all connected components of A are subgraphs of B. If true, return a list of GraphMatchers for each connected component and their corresponding subgraphs. 
        Otherwise return None.
        """
        cc_A = list(nx.connected_components(A))
        cc_A.sort(key=len) # Sorting by length speeds up the algorithm for larger systems
        cc_B = list(nx.connected_components(B))
        cc_B.sort(key=len)
        max_nodes = len(cc_B[-1])
        # self.ostream.print_info(f"number of connected components in A: {len(connected_components)}")
        for g in cc_A:
            if len(g) > max_nodes:
                return False
            A_sub = A.subgraph(g)
            GM = GraphMatcher(B, A_sub, categorical_node_match('elem', ''))
            self._iso_count += 1
            is_iso = GM.subgraph_is_isomorphic()
            if not is_iso:
                if self._check_monomorphic:
                    self._mono_count += 1
                    is_mono = GM.subgraph_is_monomorphic()
                    if not is_mono:

                        return False
                else:
                    return False
        
        # If every connected component is a subgraph, we need to verify it is a unique subgraph. 
        # This is done later, because we need to start with the largest component this time.
        cc_A.reverse()
        B = nx.Graph(B)
        for g in cc_A:
            A_sub = A.subgraph(g)
            GM = GraphMatcher(B, A_sub, categorical_node_match('elem', ''))
            self._iso_count += 1
            if not GM.subgraph_is_isomorphic():
                return False
            map = next(GM.subgraph_isomorphisms_iter())
            
            B.remove_nodes_from(map.keys())
        return True

    def _check_time(self, msg=None):
        if msg is not None:
            self.ostream.print_info(
                f"Total time = {time.time() - self._start_time:.3f} s after {msg}")
            self.ostream.flush()
            return True

        if time.time() - self._start_time > self.max_time:
            self.ostream.print_warning(
                f"Spent more then {self.max_time:.3f} s finding mapping, aborting"
            )
            return False
        return True