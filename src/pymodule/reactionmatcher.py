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
            A, B = B, A
            swapped = True

        self._start_time = time.time()
        breaking_edges = self._find_breaking_edges(A, B, forced_breaking_edges)
        if breaking_edges is None:
            return None, None, None

        insert = "monomorphisms" if self._check_monomorphic else "isomorphisms"
        self.ostream.print_info(
            f"Found subgraph {insert} with breaking bonds: {list(breaking_edges)} and forced breaking bonds: {list(forced_breaking_edges)}, continuing to find forming bonds"
        )
        self.ostream.flush()
        forming_edges = self._find_forming_edges(
            A, B, breaking_edges, forced_breaking_edges)
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

    def _find_breaking_edges(self,A,B, forced_breaking_edges):
        breaking_edges = set()
        if self._connected_components_are_subgraphs(A, B):
            return breaking_edges

        for edge in A.edges():
            A.remove_edge(*edge)
            if self._connected_components_are_subgraphs(A, B):
                breaking_edges.add(edge)
                break
            else:
                A.add_edge(*edge)
            if not self._check_time():
                return None

        self._check_time("finding one breaking bond")
        if len(breaking_edges) > 0:
            return breaking_edges
        
        self.ostream.print_info(
            "No subgraph solution with one broken bond found, trying to find two breaking bonds"
        )
        self.ostream.flush()
        for edge_a in A.edges():
            if edge_a not in A.edges():
                continue
            A.remove_edge(*edge_a)
            for edge_b in A.edges():
                A.remove_edge(*edge_b)
                if self._connected_components_are_subgraphs(A, B):
                    breaking_edges.add(edge_a)
                    breaking_edges.add(edge_b)
                    return breaking_edges
                else:
                    A.add_edge(*edge_b)
                if not self._check_time():
                    return None
            A.add_edge(*edge_a)
        self.ostream.flush()
        return None
        

    def _find_forming_edges(self, A, B, breaking_edges, forced_breaking_edges):
        """
        Find forming edges in A that are not in B, while avoiding reconnecting broken edges.
        """
        bond_count = 0
        forming_edges = []
        while len(A.edges()) < len(B.edges()):
            bond_count += 1
            edge = None
            for i in A.nodes():
                for j in A.nodes():
                    if j <= i:
                        continue
                    edge_in_A = (i, j) in A.edges() or (j, i) in A.edges()
                    edge_in_broken = (i, j) in forced_breaking_edges or (
                        j, i) in breaking_edges
                    if edge_in_A or edge_in_broken:
                        continue
                    A.add_edge(i, j)
                    if not self._connected_components_are_subgraphs(A, B):
                        A.remove_edge(i, j)
                        continue
                    if not self._check_time():
                        return None

                    edge = (i, j)
                    forming_edges.append(edge)
                    self._check_time(f"finding bond {bond_count}: {edge}")
                    self.ostream.flush()
                    break
                if edge is not None:
                    break
            if edge is None:
                self.ostream.print_info(f"Bond {i} not found, aborting")
                self.ostream.flush()
                return None
        return forming_edges

    
    def _connected_components_are_subgraphs(self, A, B):
        """
        Check if all connected components of A are subgraphs of B. If true, return a list of GraphMatchers for each connected component and their corresponding subgraphs. 
        Otherwise return None.
        """
        connected_components = list(nx.connected_components(A))
        for g in connected_components:
            A_sub = A.subgraph(g)
            GM = GraphMatcher(B, A_sub, categorical_node_match('elem', ''))
            self._iso_count += 1
            if not GM.subgraph_is_isomorphic():
                if self._check_monomorphic:
                    self._mono_count += 1
                    if not GM.subgraph_is_monomorphic():
                        return False
                else:
                    return False
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