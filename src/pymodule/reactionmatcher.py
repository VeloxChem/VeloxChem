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
import math

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
        self._verbose = False
        self.max_time = 600

        self.max_break_attempts_guess = 1e5
        self._breaking_depth = 1  # how many breaking edges to try
        self._check_monomorphic = False

    def get_mapping(
        self,
        reactant_ff,
        rea_elems,
        product_ff,
        pro_elems,
        breaking_bonds,
        forming_bonds,
    ):

        self.rea_graph, self.pro_graph = self._create_reaction_graphs(
            reactant_ff,
            rea_elems,
            product_ff,
            pro_elems,
        )

        N = len(self.rea_graph.edges)
        while math.factorial(N) / math.factorial(
                N - self._breaking_depth) < self.max_break_attempts_guess:
            self._breaking_depth += 1
            if self._breaking_depth == N:
                break

        self.ostream.print_info(
            f"Set max breaking depth to {self._breaking_depth} from max breaking attempts {self.max_break_attempts_guess}"
        )
        self.ostream.flush()

        map, breaking_edges, forming_edges = self._find_mapping(
            self.rea_graph.copy(),
            self.pro_graph.copy(),
            breaking_bonds,
            forming_bonds,
        )

        self.ostream.flush()
        self.ostream.print_info(
            f"Total subgraph isomorphism checks: {self._iso_count}, total subgraph monomorphism checks: {self._mono_count}"
        )
        return map, breaking_edges, forming_edges

    def _create_reaction_graphs(self, reactant_ff, rea_elems, product_ff,
                                pro_elems):
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

    def _find_mapping(self,
                      A,
                      B,
                      forced_breaking_edges=set(),
                      forced_forming_edges=set()):
        """
        Find a mapping between the connected components of A and B, while avoiding reconnecting broken edges.
        """
        for edge in forced_breaking_edges:
            A.remove_edge(*edge)

        for edge in forced_forming_edges:
            A.add_edge(*edge)

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
        breaking_edges = self._find_breaking_edges(
            A,
            B,
            forced_forming_edges,
        )
        if breaking_edges is None:
            return None, None, None

        insert = "monomorphisms" if self._check_monomorphic else "isomorphisms"
        self.ostream.print_info(
            f"Found subgraph {insert} with breaking bonds: {self._print_bond_list(breaking_edges)} and forced breaking bonds: {self._print_bond_list(forced_breaking_edges)}, continuing to find forming bonds"
        )
        self.ostream.flush()

        forming_edges = self._find_forming_edges(
            A, B, breaking_edges | forced_breaking_edges)
        if forming_edges is None:
            if not self._check_monomorphic:
                self.ostream.print_info(
                    f"Could not find forming edges from subgraph isomorphism, attempting subgraph monomorphism"
                )
                self._check_monomorphic = True
                forming_edges = self._find_forming_edges(
                    A, B, breaking_edges | forced_breaking_edges)
        if forming_edges is None:
            return None, None, None

        self.ostream.print_info(
            f"Found forming bonds: {self._print_bond_list(forming_edges)} with forced forming bonds: {self._print_bond_list(forced_forming_edges)}. Finding mapping from isomorphism."
        )
        self.ostream.flush()
        # print(forming_edges)
        GM = GraphMatcher(A, B, categorical_node_match('elem', ''))
        map = next(GM.isomorphisms_iter())

        self._check_time(f"finding mapping.")
        if swapped:
            # if we swapped A and B, we need to swap the mapping back
            self.ostream.print_info(
                "Inverting mapping because A and B were swapped.")
            map = {v: k for k, v in map.items()}
        return map, breaking_edges, forming_edges

    def _find_breaking_edges(self, A, B, forced_forming_edges):
        breaking_edges = None
        B_composition = self._get_bond_element_composition(B)
        for depth in range(0, self._breaking_depth):
            self._check_monomorphic = False
            self.ostream.print_info(
                f"Finding isomorphic breaking edges at depth {depth}")
            self.ostream.flush()
            breaking_edges = self._recurse_breaking_edges(
                A,
                B,
                depth,
                B_composition,
                forced_forming_edges,
            )
            if breaking_edges is not None:
                break
            self._check_monomorphic = True
            self.ostream.print_info(
                f"Finding monomorphic breaking edges at depth {depth}")
            self.ostream.flush()
            breaking_edges = self._recurse_breaking_edges(
                A,
                B,
                depth,
                B_composition,
                forced_forming_edges,
            )
            if breaking_edges is not None:
                break
            else:
                self.ostream.print_info(
                    f"No breaking edges found at depth {depth}, increasing depth."
                )
                self.ostream.flush()
        return breaking_edges

    def _recurse_breaking_edges(self, A, B, depth, B_composition,
                                forced_forming_edges):
        """
        Search the minimal amount of edges that need to be removed from A so that all connected components of A are subgraphs of B. 
        Searches recursively up to given depth. Returns None if no edges are found at a given depth, or when the time limit is exceeded.
        In looping, it prioritizes breaking bonds involving combinations of elements that are more abundant in A then in B. 
        After that, it prioritizes breaking bonds involving Hydrogen atoms.
        """

        if not self._check_time():
            return None

        if self._connected_components_are_subgraphs(A, B):
            return set()

        if depth > 0:

            # B will have more bonds then A, if there is a type of bond in A of which there are more in A then in B, then those bonds for sure need to be removed
            A_composition = self._get_bond_element_composition(A)
            required_element_combinations = set()
            for key, val in A_composition.items():
                if B_composition.get(key, 0) < val:
                    required_element_combinations.add(key)

            # Sorted edges has all edges involving H first, prioritizing a solution that breaks bonds with Hydrogen atoms
            sorted_edges = self._sort_edges(A)
            for edge in sorted_edges:
                if len(required_element_combinations) > 0:
                    comb = self._get_elem_comb(edge[0], edge[1], A)
                    if comb not in required_element_combinations:
                        continue
                if edge in forced_forming_edges or (
                        edge[1], edge[0]) in forced_forming_edges:
                    continue

                A.remove_edge(*edge)
                if self._verbose:
                    self.ostream.print_info(
                        f"Trying to find breaking edge {self._print_bond(edge)} at depth {depth}"
                    )
                self.ostream.flush()
                edges = self._recurse_breaking_edges(
                    A,
                    B,
                    depth - 1,
                    B_composition,
                    forced_forming_edges,
                )
                if edges is not None:
                    edges.add(edge)
                    self.ostream.print_info(
                        f"Found breaking edges: {self._print_bond_list(edges)} at depth {depth}."
                    )
                    self.ostream.flush()
                    return edges
                A.add_edge(*edge)
        return None

    def _find_forming_edges(self, A, B, breaking_edges):
        """
        Keep adding edges to A until it is isomorphic to B, while avoiding reconnecting broken edges.
        """
        bond_count = 0
        forming_edges = []
        nodes = list(A.nodes())
        active_nodes = set()
        for edge in breaking_edges:
            active_nodes.add(edge[0])
            active_nodes.add(edge[1])

        # Move every element of active_nodes in nodes to the front of the list
        nodes = [n for n in nodes if n in active_nodes
                 ] + [n for n in nodes if n not in active_nodes]

        elem_combinations_B = self._get_bond_element_composition(B)
        elem_combinations_A = self._get_bond_element_composition(A)

        while len(A.edges()) < len(B.edges()):
            bond_count += 1
            edge = None
            for i, node_i in enumerate(nodes):
                for j, node_j in enumerate(nodes):
                    if j <= i:
                        continue
                    # If an edge is in A, it cannot be formed
                    edge_in_A = (node_i,
                                 node_j) in A.edges() or (node_j,
                                                          node_i) in A.edges()
                    # If an edge was broken, it should not be re-formed
                    edge_in_broken = (node_i, node_j) in breaking_edges or (
                        node_j, node_i) in breaking_edges

                    # HH_bond = A.nodes[node_i]['elem'] == 1.0 and A.nodes[node_j]['elem'] == 1.0
                    if edge_in_A or edge_in_broken:
                        continue

                    comb = self._get_elem_comb(node_i, node_j, A)
                    if elem_combinations_B.get(comb,
                                               0) <= elem_combinations_A.get(
                                                   comb, 0):
                        continue

                    A.add_edge(node_i, node_j)
                    if self._verbose:
                        self.ostream.print_info(
                            f"Trying to find bond {bond_count}: {self._print_bond((node_i, node_j))}"
                        )
                        self.ostream.flush()
                    if not self._connected_components_are_subgraphs(A, B):
                        A.remove_edge(node_i, node_j)
                        continue
                    if not self._check_time():
                        return None

                    edge = (node_i, node_j)
                    forming_edges.append(edge)
                    self._check_time(
                        f"finding bond {bond_count}: {self._print_bond(edge)}")
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
    def _get_bond_element_composition(A):
        composition = {}
        for edge in A.edges():
            comb = ReactionMatcher._get_elem_comb(edge[0], edge[1], A)
            if comb not in composition.keys():
                composition[comb] = 1
            else:
                composition[comb] += 1
        return composition

    @staticmethod
    def _sort_edges(A):
        H_edges = []
        other_edges = []
        # TODO: add more sorting, for example based on metals
        for edge in A.edges():
            comb = ReactionMatcher._get_elem_comb(edge[0], edge[1], A)
            if 1.0 in comb:
                H_edges.append(edge)
            else:
                other_edges.append(edge)
        sorted_edges = H_edges + other_edges
        return sorted_edges

    @staticmethod
    def _get_elem_comb(node1, node2, A):
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
        cc_A.sort(
            key=len
        )  # Sorting by length speeds up the algorithm for larger systems
        cc_B = list(nx.connected_components(B))
        cc_B.sort(key=len)
        max_nodes = len(cc_B[-1])
        for g in cc_A:
            if len(g) > max_nodes:
                return False
            A_sub = A.subgraph(g)
            GM = GraphMatcher(B, A_sub, categorical_node_match('elem', ''))
            if self._check_monomorphic:
                self._mono_count += 1
                is_mono = GM.subgraph_is_monomorphic()
                if not is_mono:
                    return False
            else:
                self._iso_count += 1
                is_iso = GM.subgraph_is_isomorphic()
                if not is_iso:
                    return False

        # If every connected component is a subgraph, we need to verify it is a unique subgraph.
        # This is done after the first routine, because we need to start with the largest component this time.
        cc_A.reverse()
        B = nx.Graph(B)
        
        # Loop over every connected subgraph in A in decreasing size
        # Look for all maps from A_sub to B and find the mapping that leaves the least amount of connected components in B when removing the mapped nodes
        for g in cc_A:
            A_sub = A.subgraph(g)
            GM = GraphMatcher(B, A_sub, categorical_node_match('elem', ''))
            if self._check_monomorphic:
                self._mono_count += 1
                if not GM.subgraph_is_monomorphic():
                    return False
            else:
                self._iso_count += 1
                if not GM.subgraph_is_isomorphic():
                    return False

            # find mapping that leaves least amount of connected components
            if self._check_monomorphic:
                iterator = GM.subgraph_monomorphisms_iter()
                self._mono_count += 1
            else:
                iterator = GM.subgraph_isomorphisms_iter()
                self._iso_count += 1
                
            best_map = next(iterator)
            
            B_copy = nx.Graph(B)
            B_copy.remove_nodes_from(best_map.keys())
            best_cc_count = len(list(nx.connected_components(B_copy)))
            for map in iterator:
                if self._check_monomorphic:
                    self._mono_count += 1
                else:
                    self._iso_count += 1
                    
                B_copy = nx.Graph(B)
                B_copy.remove_nodes_from(map.keys())
                cc_count = len(list(nx.connected_components(B_copy)))
                if cc_count < best_cc_count:
                    best_cc_count = cc_count
                    best_map = map
                    
                # The least amount of connected components is 1, or 0 if there is nothing left
                # TODO are there ways to predict that the minimum needs to be higher?
                if best_cc_count <= 1:
                    break


            B.remove_nodes_from(best_map.keys())
        return True

    def _check_time(self, msg=None):
        if msg is not None:
            self.ostream.print_info(
                f"Total time = {time.time() - self._start_time:.3f} s after {msg}"
            )
            self.ostream.flush()
            return True

        if time.time() - self._start_time > self.max_time:
            self.ostream.print_warning(
                f"Spent more then {self.max_time:.3f} s finding mapping, aborting"
            )
            return False
        return True

    @staticmethod
    def _print_bond_list(bonds):
        bonds = list(bonds)
        return f"[{', '.join([ReactionMatcher._print_bond(bond) for bond in bonds])}]"

    @staticmethod
    def _print_bond(bond):
        return f"({bond[0]+1}, {bond[1]+1})"
