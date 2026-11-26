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
from itertools import permutations, product
import networkx as nx
import numpy as np
import time
import sys
import copy
import math

from .outputstream import OutputStream
from .veloxchemlib import mpi_master
from .errorhandler import assert_msg_critical
from collections import Counter


class _TimeLimitExceeded(Exception):
    """
    Raised once the deadline of ReactionMatcher.get_mapping has passed.
    """


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
        self._deadline = 0
        self._time_limit_active = False
        self.verbose = False
        self.print_starting_index = 1
        self.max_time = 600
        self.force_hydrogen_inclusion = True
        self.brute_force_lim = 1000
        self.rdFMCS_assist_assign = False
        self.rdRascalMCES_assist_assign = False
        self._reduce_hydrogen = False
        self._assist_min_depth = 3

        self.max_break_attempts_guess = 1e7
        self._breaking_depth = -1  # how many breaking edges to try
        self._check_monomorphic = True
        self._assisting_map = {}

    def get_mapping(
        self,
        reactant,
        product,
        breaking_bonds=None,
        forming_bonds=None,
    ):
        """
        Find the mapping from the reactant atoms to the product atoms, with the breaking and forming bonds.
        Returns None for all results if no mapping exists. Exceeding max_time, which bounds the entire search
        including the retry with hydrogens, is a critical error.
        """
        self._start_time = time.monotonic()
        self._deadline = self._start_time + self.max_time
        self._time_limit_active = True
        finished_in_time = True
        try:
            map, breaking_edges, forming_edges = self._get_mapping(
                reactant, product, breaking_bonds, forming_bonds)
        except _TimeLimitExceeded:
            finished_in_time = False
        finally:
            self._time_limit_active = False

        assert_msg_critical(
            finished_in_time,
            f"ReactionMatcher: time limit of {self.max_time:.3f} s exceeded, "
            + "the search for a mapping was aborted. This does not mean "
            + "that no mapping exists. Increase max_time to search longer.")
        return map, breaking_edges, forming_edges

    def _get_mapping(self, reactant, product, breaking_bonds, forming_bonds):
        map, breaking_edges, forming_edges = self._match_reaction_graphs(
            reactant,
            product,
            breaking_bonds,
            forming_bonds,
            self.force_hydrogen_inclusion,
        )
        if map is None and self._reduce_hydrogen:
            self.ostream.print_info(
                "No mapping found, retrying with hydrogens included.")
            self.ostream.flush()
            map, breaking_edges, forming_edges = self._match_reaction_graphs(
                reactant,
                product,
                breaking_bonds,
                forming_bonds,
                True,
            )

        self.ostream.print_info(
            f"Total subgraph isomorphism checks: {self._iso_count}, total subgraph monomorphism checks: {self._mono_count}"
        )
        self.ostream.flush()
        return map, breaking_edges, forming_edges

    def _match_reaction_graphs(self, reactant, product, breaking_bonds,
                               forming_bonds, force_hydrogen_inclusion):
        self.rea_graph, self.pro_graph = self._create_reaction_graphs(
            reactant, product, breaking_bonds, forming_bonds,
            force_hydrogen_inclusion)

        self.rea_graph, self.pro_graph, self.assisting_map = self._assign_assist_ids(
            reactant, product, self.rea_graph, self.pro_graph)

        if len(self.assisting_map) == len(self.rea_graph.nodes):
            map = dict(self.assisting_map)
            breaking_edges, forming_edges = set(), set()
        else:
            if self._breaking_depth == -1:
                self._breaking_depth = self._decide_breaking_depth()

            map, breaking_edges, forming_edges = self._find_mapping(
                self.rea_graph.copy(),
                self.pro_graph.copy(),
                breaking_bonds,
                forming_bonds,
            )

        if map is not None and self._reduce_hydrogen:
            self._restore_hydrogen(map, self.rea_graph, self.pro_graph)
        return map, breaking_edges, forming_edges

    def _decide_breaking_depth(self):
        N = len(self.rea_graph.edges)
        _breaking_depth = 1

        if N == 1:
            return 1

        # Equivalent to N! / (N - depth)!
        def estimate_attempts(N, depth):
            lg = math.lgamma(N + 1) - math.lgamma(N - depth + 1)
            return math.exp(lg)

        while estimate_attempts(
                N, _breaking_depth) < self.max_break_attempts_guess:
            _breaking_depth += 1
            if _breaking_depth == N:
                break

        self.ostream.print_info(
            f"Set max breaking depth to {_breaking_depth} from max breaking attempts {self.max_break_attempts_guess}"
        )
        self.ostream.flush()
        return _breaking_depth

    def _create_reaction_graphs(self, reactant, product, breaking_bonds,
                                forming_bonds, force_hydrogen_inclusion):

        rea_graph = self._prepare_graph(reactant, breaking_bonds, forming_bonds)

        pro_graph = self._prepare_graph(product)

        self._reduce_hydrogen = False
        if not force_hydrogen_inclusion:
            self._reduce_hydrogen = self._decide_hydrogen_reduction(
                rea_graph, pro_graph)
            if self._reduce_hydrogen:
                self.ostream.print_info("Reducing hydrogens into parent atoms")
                self.ostream.flush()
                rea_hydrogens = [
                    n for n in rea_graph.nodes
                    if rea_graph.nodes[n]['elem'] == 1.0
                ]
                pro_hydrogens = [
                    n for n in pro_graph.nodes
                    if pro_graph.nodes[n]['elem'] == 1.0
                ]
                rea_graph.remove_nodes_from(rea_hydrogens)
                pro_graph.remove_nodes_from(pro_hydrogens)

        return rea_graph, pro_graph

    def _prepare_graph(self, molecule, breaking_bonds=None, forming_bonds=None):
        graph = nx.Graph()
        mat = molecule.get_connectivity_matrix()
        bonds = set()
        for i in range(len(mat)):
            for j in range(i, len(mat)):
                if mat[i, j]:
                    bonds.add((i, j))
        # Remove the bonds that are being broken, so that these segments get treated as seperate reactants

        graph.add_nodes_from(list(range(len(mat))))
        graph.add_edges_from(bonds)

        if breaking_bonds is not None:
            for edge in breaking_bonds:
                graph.remove_edge(*edge)

        if forming_bonds is not None:
            for edge in forming_bonds:
                graph.add_edge(*edge)

        elements = [float(id) for id in molecule.get_element_ids()]
        for i, elem in enumerate(elements):
            graph.nodes[i]['elem'] = elem

        for i in graph.nodes:
            connected_indices, connected_elements = self._get_connected_atoms(
                graph, i)
            H_count = sum(1 for elem in connected_elements if elem == 1.0)
            graph.nodes[i]['H_bond_count'] = H_count
            if H_count > 0:
                graph.nodes[i]['H_indices'] = [
                    idx
                    for idx, elem in zip(connected_indices, connected_elements)
                    if elem == 1.0
                ]

        return graph

    def _get_connected_atoms(self, graph, node):
        indices = list(graph.neighbors(node))
        elements = [graph.nodes[n]['elem'] for n in indices]
        return indices, elements

    def _get_range_subgraph(self, graph, node, cutoff):
        nodes = [
            node for node, dist in nx.single_source_shortest_path_length(
                graph, node, cutoff=cutoff).items()
        ]
        return graph.subgraph(nodes)

    def _decide_hydrogen_reduction(self, A, B):
        # Small enough graphs don't need hydrogen reduction
        if len(A.nodes) < 30:
            return False

        # Returns true if the types of bonds that are involving hydrogens are balanced between A and B
        A_bond_composition = self._get_bond_element_composition(A)
        B_bond_composition = self._get_bond_element_composition(B)
        hydrogen_combs = {
            comb
            for comb in set(A_bond_composition) | set(B_bond_composition)
            if 1.0 in comb
        }
        for comb in hydrogen_combs:
            if A_bond_composition.get(comb,
                                      0) != B_bond_composition.get(comb, 0):
                return False
        return True

    def _assign_assist_ids(self, reactant, product, rea_graph, pro_graph):
        # Tries to figure out the 'obvious' atom mappings before starting the full graph matching
        # This speeds up the matching significantly for larger molecules
        # Subgraphs of decreasing size in the reactant and product and product are itterated over
        # If they are isomorphic, it is assumed that these structures are the same in both the reactant and the product

        self.ostream.print_info(
            "Searching for matching subgraphs to assign assist ids")
        self.ostream.flush()
        if self.rdFMCS_assist_assign or self.rdRascalMCES_assist_assign:
            return self._get_rdkit_assist_ids(reactant, product, rea_graph,
                                              pro_graph)
        else:
            return self._calculate_assist_ids(rea_graph, pro_graph)

    def _get_rdkit_assist_ids(self, reactant, product, rea_graph, pro_graph):
        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolfiles
            from rdkit.Chem import rdRascalMCES
            from rdkit.Chem import rdFMCS

            from rdkit.Chem import rdDetermineBonds

            rd_mol1 = Chem.RWMol(
                rdmolfiles.MolFromXYZBlock(reactant.get_xyz_string()))
            rd_mol2 = Chem.RWMol(
                rdmolfiles.MolFromXYZBlock(product.get_xyz_string()))

            rdDetermineBonds.DetermineConnectivity(rd_mol1)
            rdDetermineBonds.DetermineBondOrders(rd_mol1)
            rdDetermineBonds.DetermineConnectivity(rd_mol2)
            rdDetermineBonds.DetermineBondOrders(rd_mol2)
            if self.rdRascalMCES_assist_assign:

                res = rdRascalMCES.FindMCES(rd_mol1, rd_mol2)
                if len(res) > 0:
                    assist_map = {i: j for i, j in res[0].atomMatches()}
                    for rea_node, pro_node in assist_map.items():
                        rea_graph.nodes[rea_node]['assist_id'] = rea_node
                        pro_graph.nodes[pro_node]['assist_id'] = rea_node
                    self.ostream.print_info(str(res[0].atomMatches()))
                    self.ostream.print_info(
                        f"Assigned {len(assist_map)} assist ids: {self._print_mapping(assist_map)}"
                    )
                    self.ostream.flush()
                    self.assisting_map = assist_map
                    return rea_graph, pro_graph, assist_map
                return rea_graph, pro_graph, {}

            elif self.rdFMCS_assist_assign:

                mcs = rdFMCS.FindMCS([rd_mol1, rd_mol2])
                if mcs.numAtoms > 0:
                    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                    match1 = rd_mol1.GetSubstructMatch(mcs_mol)
                    match2 = rd_mol2.GetSubstructMatch(mcs_mol)
                    # target_atm1 = []
                    assist_map = {}
                    for atom1, atom2 in zip(match1, match2):
                        assist_map[atom1] = atom2
                        rea_graph.nodes[atom1]['assist_id'] = atom1
                        pro_graph.nodes[atom2]['assist_id'] = atom1

                    self.ostream.print_info(
                        f"Assigned {len(assist_map)} assist ids: {self._print_mapping(assist_map)}"
                    )
                    self.ostream.flush()
                    self.assisting_map = assist_map
                    return rea_graph, pro_graph, assist_map
                return rea_graph, pro_graph, {}
        except ImportError:
            self.ostream.print_warning(
                "rdkit not available, skipping rdkit based assist id assignment. Turn off rdFMCS_assist_assign and rdRascalMCES_assist_assign to use VeloxChem assist id assignment."
            )
            self.ostream.flush()
            return rea_graph, pro_graph, {}

    def _calculate_assist_ids(self, rea_graph, pro_graph):
        # Tries to figure out the 'obvious' atom mappings before starting the full graph matching
        # This speeds up the matching significantly for larger molecules
        # Neighbourhoods of decreasing range are compared through Weisfeiler-Lehman hashes. A heavy atom whose
        # neighbourhood occurs exactly once among the unassigned atoms of both the reactant and the product is a
        # centre. Inside the matching neighbourhoods of a centre, the atoms that can be told apart are mapped as well.
        # Atoms that only look equivalent locally are never mapped, because picking one of several local
        # isomorphisms depends on the search order of the graph matcher and can be chemically wrong.
        # The hashes include the assist ids, so every assigned atom can make the atoms around it unique.

        min_depth_positive = self._assist_min_depth >= 1
        assert_msg_critical(
            min_depth_positive,
            f"ReactionMatcher: the minimum assist depth must be at least 1, got {self._assist_min_depth}"
        )

        assisting_map = {}

        depth = 1
        for cc in nx.connected_components(rea_graph):
            dia = nx.diameter(rea_graph.subgraph(cc))
            r = math.ceil(dia / 2)
            if r > depth:
                depth = r
        depth += 1

        # Enough iterations for hashes that see an entire reactant molecule
        full_depth = 2 * depth

        while depth >= self._assist_min_depth:
            self._check_time()
            centres = self._find_unique_centres(rea_graph, pro_graph, depth)
            new_count = self._assign_centres(rea_graph, pro_graph, depth,
                                             centres, assisting_map)

            # Once no atom is unique anymore, break the symmetry of one set of equivalent atoms by assigning one of
            # them. This makes the other atoms of the set unique for the next round. Only one set is assigned per
            # round, because assigning one set can break the symmetry of another.
            if new_count == 0:
                for centre in self._find_symmetric_centres(
                        rea_graph, pro_graph, depth, full_depth):
                    new_count = self._assign_centres(rea_graph, pro_graph,
                                                     depth, [centre],
                                                     assisting_map)
                    if new_count > 0:
                        break

            # Lower the range only once no new atoms can be assigned at this one
            if new_count == 0:
                depth -= 1
            else:
                self.ostream.print_info(
                    f"Assigned {new_count} heavy atom assist ids at depth {depth}."
                )
                self.ostream.flush()

        sorted_assisting_map = {
            k: assisting_map[k]
            for k in sorted(assisting_map)
        }
        self.ostream.print_info(
            f"Assigned {len(assisting_map)} assist ids: {self._print_mapping(sorted_assisting_map)}"
        )
        self._assisting_map = sorted_assisting_map
        self.ostream.flush()
        return rea_graph, pro_graph, sorted_assisting_map

    def _find_unique_centres(self, rea_graph, pro_graph, depth):
        """
        Find the unassigned heavy atoms whose neighbourhood within depth occurs exactly once in both graphs.
        Returns a sorted list of (reactant node, product node) pairs.
        """
        rea_hashes = self._get_neighbourhood_hashes(rea_graph, depth)
        pro_hashes = self._get_neighbourhood_hashes(pro_graph, depth)
        rea_unique = self._get_unique_hashes(rea_graph, rea_hashes)
        pro_unique = self._get_unique_hashes(pro_graph, pro_hashes)
        pairs = [(rea_node, pro_unique[h]) for h, rea_node in rea_unique.items()
                 if h in pro_unique]
        pairs.sort()
        return pairs

    def _map_neighbourhoods(self, rea_graph, rea_id, pro_graph, pro_id,
                            depth):
        """
        Map the unassigned heavy atoms within depth - 1 of rea_id and pro_id that can be told apart inside the
        neighbourhoods within depth. The mapping is only returned if the neighbourhoods are isomorphic with all of
        its atoms pinned, otherwise an empty mapping is returned.
        """
        rea_subgraph = self._get_range_subgraph(rea_graph, rea_id,
                                                depth).copy()
        pro_subgraph = self._get_range_subgraph(pro_graph, pro_id,
                                                depth).copy()
        if len(rea_subgraph.nodes) != len(pro_subgraph.nodes):
            return {}

        # Pin the centres to each other through an assist id that no other atom has
        rea_subgraph.nodes[rea_id]['assist_id'] = -2
        pro_subgraph.nodes[pro_id]['assist_id'] = -2

        # 2 * depth iterations let the hashes see the entire neighbourhood
        rea_hashes = self._get_neighbourhood_hashes(rea_subgraph, 2 * depth)
        pro_hashes = self._get_neighbourhood_hashes(pro_subgraph, 2 * depth)
        rea_inner = self._get_range_subgraph(rea_graph, rea_id, depth - 1).nodes
        pro_inner = self._get_range_subgraph(pro_graph, pro_id, depth - 1).nodes

        map = {rea_id: pro_id}
        pro_unique = self._get_unique_hashes(pro_subgraph, pro_hashes)
        for h, rea_node in self._get_unique_hashes(rea_subgraph,
                                                   rea_hashes).items():
            pro_node = pro_unique.get(h, None)
            if rea_node in rea_inner and pro_node in pro_inner:
                map[rea_node] = pro_node

        for rea_node, pro_node in map.items():
            rea_subgraph.nodes[rea_node]['assist_id'] = rea_node
            pro_subgraph.nodes[pro_node]['assist_id'] = rea_node
        GM = self.get_graph_matcher(rea_subgraph, pro_subgraph)
        if not GM.is_isomorphic():
            return {}
        return map

    def _assign_centres(self, rea_graph, pro_graph, depth, centres,
                        assisting_map):
        """
        Assign assist ids to the centres and to the atoms that can be told apart in their neighbourhoods.
        Returns the number of assigned heavy atoms.
        """
        count = 0
        for rea_id, pro_id in centres:

            # An earlier centre can have assigned this one already
            if rea_graph.nodes[rea_id].get('assist_id') is not None:
                continue
            if pro_graph.nodes[pro_id].get('assist_id') is not None:
                continue

            map = self._map_neighbourhoods(rea_graph, rea_id, pro_graph,
                                           pro_id, depth)
            if len(map) == 0:
                continue

            if self.verbose:
                self.ostream.print_info(
                    f"Found matching neighbourhood between reactant atom {rea_id + self.print_starting_index} and product atom {pro_id + self.print_starting_index} at depth {depth} with mapping: {self._print_mapping(map)}."
                )
                self.ostream.flush()

            for rea_node, pro_node in map.items():
                rea_graph.nodes[rea_node]['assist_id'] = rea_node
                pro_graph.nodes[pro_node]['assist_id'] = rea_node
                assisting_map.update({rea_node: pro_node})

                # Reduced hydrogens are not in the graphs, _restore_hydrogen maps them at the end
                if self._reduce_hydrogen:
                    continue
                rea_H_indices = rea_graph.nodes[rea_node].get('H_indices')
                pro_H_indices = pro_graph.nodes[pro_node].get('H_indices')
                if rea_H_indices is not None and pro_H_indices is not None:
                    if len(rea_H_indices) == len(pro_H_indices):
                        for i, j in zip(rea_H_indices, pro_H_indices):
                            rea_graph.nodes[i]['assist_id'] = i
                            pro_graph.nodes[j]['assist_id'] = i
                            assisting_map.update({i: j})
            count += len(map)
        return count

    def _find_symmetric_centres(self, rea_graph, pro_graph, depth,
                                full_depth):
        """
        Find the sets of unassigned heavy atoms that are equivalent within the entire reactant and the entire product,
        and that are the only atoms with their neighbourhood within depth in both graphs.
        Yields one (reactant node, product node) pair of each set, ordered by the lowest reactant node. Assigning one
        atom of such a set does not exclude any mapping, because the equivalent atoms can be swapped by a symmetry of
        the graph.
        """
        rea_groups = self._group_unassigned_by_hash(
            rea_graph, self._get_neighbourhood_hashes(rea_graph, depth))
        pro_groups = self._group_unassigned_by_hash(
            pro_graph, self._get_neighbourhood_hashes(pro_graph, depth))

        rea_full = None
        pro_full = None
        for h in sorted(rea_groups, key=lambda h: min(rea_groups[h])):
            rea_nodes = rea_groups[h]
            pro_nodes = pro_groups.get(h, [])
            if len(rea_nodes) < 2 or len(rea_nodes) != len(pro_nodes):
                continue

            # Hashes over the entire graph are a cheap first test of a symmetry beyond the neighbourhood
            if rea_full is None:
                rea_full = self._get_neighbourhood_hashes(rea_graph, full_depth)
                pro_full = self._get_neighbourhood_hashes(pro_graph, full_depth)
            if len({rea_full[n] for n in rea_nodes}) != 1:
                continue
            if len({pro_full[n] for n in pro_nodes}) != 1:
                continue

            if not self._are_symmetric(rea_graph, rea_nodes):
                continue
            if not self._are_symmetric(pro_graph, pro_nodes):
                continue
            yield rea_nodes[0], pro_nodes[0]

    def _are_symmetric(self, graph, nodes):
        """
        Check that a symmetry of graph, keeping all assigned atoms in place, maps the first node on each other node.
        Equal Weisfeiler-Lehman hashes do not guarantee this in general.
        """
        first = graph.copy()
        first.nodes[nodes[0]]['assist_id'] = -2
        for n in nodes[1:]:
            other = graph.copy()
            other.nodes[n]['assist_id'] = -2
            GM = self.get_graph_matcher(first, other)
            if not GM.is_isomorphic():
                return False
        return True

    def _get_neighbourhood_hashes(self, graph, depth):
        """
        Get the Weisfeiler-Lehman hash of the neighbourhood within depth of every node in graph.
        The hashes take the same node attributes into account as the graph matcher.
        """
        # The brackets keep the concatenated labels of the first iteration unambiguous
        labelled = nx.Graph()
        for n, data in graph.nodes(data=True):
            label = f"[{data['elem']}|{data.get('assist_id', -1)}"
            if self._reduce_hydrogen:
                label += f"|{data.get('H_bond_count', 0)}"
            labelled.add_node(n, label=label + "]")
        labelled.add_edges_from(graph.edges)

        hashes = nx.weisfeiler_lehman_subgraph_hashes(labelled,
                                                      node_attr='label',
                                                      iterations=depth)
        return {n: h[-1] for n, h in hashes.items()}

    @staticmethod
    def _group_unassigned_by_hash(graph, hashes):
        """
        Group the unassigned heavy atoms of graph by their hash.
        """
        groups = {}
        for n in graph.nodes:
            if graph.nodes[n]['elem'] == 1.0:
                continue
            if graph.nodes[n].get('assist_id', None) is not None:
                continue
            groups.setdefault(hashes[n], []).append(n)
        return groups

    @staticmethod
    def _get_unique_hashes(graph, hashes):
        """
        Get the hashes that belong to exactly one unassigned heavy atom of graph, mapped to that atom.
        """
        groups = ReactionMatcher._group_unassigned_by_hash(graph, hashes)
        return {h: nodes[0] for h, nodes in groups.items() if len(nodes) == 1}

    def _find_mapping(self, A, B, forced_breaking_edges, forced_forming_edges):
        """
        Find a mapping between the connected components of A and B, while avoiding reconnecting broken edges.
        """

        if forced_breaking_edges is None:
            forced_breaking_edges = set()
        if forced_forming_edges is None:
            forced_forming_edges = set()

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

            self.assisting_map = {v: k for k, v in self.assisting_map.items()}

        # self.assisting_map = {v: k for k, v in self.assisting_map.items()}
        A_assist = self._get_assisting_graph(A)
        B_assist = self._get_assisting_graph(B)
        # B_assist = nx.relabel_nodes(B_assist, self.assisting_map)

        unassigned_nodes_A = A.nodes - A_assist.nodes
        unassigned_elems_count = list(
            Counter([A.nodes[i]['elem'] for i in unassigned_nodes_A]).values())
        total_expected_perms = 1
        for elem_count in unassigned_elems_count:
            total_expected_perms *= math.factorial(elem_count)

        if total_expected_perms <= self.brute_force_lim:
            unassigned_nodes_B = B.nodes - B_assist.nodes

            # Only nodes of the same element are permuted, because other permutations are not chemically valid
            elems = sorted({A.nodes[n]['elem'] for n in unassigned_nodes_A})
            ordered_nodes_A = [
                n for elem in elems for n in sorted(unassigned_nodes_A)
                if A.nodes[n]['elem'] == elem
            ]
            elem_permutations = [
                permutations(
                    sorted(n for n in unassigned_nodes_B
                           if B.nodes[n]['elem'] == elem)) for elem in elems
            ]
            permuted_unassigned_nodes_B = [
                sum(perm, ()) for perm in product(*elem_permutations)
            ]
            self.ostream.print_info(
                f"Trying all {len(permuted_unassigned_nodes_B)} permutations of unassigned nodes of the same element. Expecting to consider {total_expected_perms} permutations based on element counts."
            )
            self.ostream.flush()
            return self._find_brute_force_mapping(A, B, forced_breaking_edges,
                                                  forced_forming_edges, swapped,
                                                  ordered_nodes_A,
                                                  permuted_unassigned_nodes_B)

        breaking_edges = self._find_breaking_edges(
            A,
            A_assist,
            B,
            B_assist,
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
            A, A_assist, B, B_assist, breaking_edges | forced_breaking_edges)
        if forming_edges is None:
            if not self._check_monomorphic:
                self.ostream.print_info(
                    "Could not find forming edges from subgraph isomorphism, attempting subgraph monomorphism"
                )
                self.ostream.flush()
                self._check_monomorphic = True
                forming_edges = self._find_forming_edges(
                    A, A_assist, B, B_assist,
                    breaking_edges | forced_breaking_edges)
        if forming_edges is None:
            return None, None, None

        self.ostream.print_info(
            f"Found forming bonds: {self._print_bond_list(forming_edges)} with forced forming bonds: {self._print_bond_list(forced_forming_edges)}. Finding mapping from isomorphism."
        )
        self.ostream.flush()
        # print(forming_edges)
        GM = self.get_graph_matcher(A, B)
        map = next(GM.isomorphisms_iter())

        self._print_elapsed("finding mapping.")
        if swapped:
            # if we swapped A and B, we need to swap the mapping back
            self.ostream.print_info(
                "Inverting mapping because A and B were swapped.")
            map = {v: k for k, v in map.items()}
        return map, breaking_edges, forming_edges

    @staticmethod
    def _get_assisting_graph(A, relabel=False):
        A_assisting = nx.Graph(A)
        nodes = list(A_assisting.nodes().keys())
        for n in nodes:
            assist_id = A_assisting.nodes[n].get('assist_id', None)
            if assist_id is None:
                A_assisting.remove_node(n)

        return A_assisting

    def _find_brute_force_mapping(self, A, B, forced_breaking_edges,
                                  forced_forming_edges, swapped,
                                  unassigned_nodes_A,
                                  permuted_unassigned_nodes_B):
        maps = []
        breaking_edges_list = []
        forming_edges_list = []
        N_changing_bonds = []
        N_changing_H_bonds = []

        # Loop through all permutations
        for perm in permuted_unassigned_nodes_B:
            self._check_time()
            temp_map = {k: v for k, v in zip(unassigned_nodes_A, perm)}
            skip = False

            # Check if the map is chemically valid
            for k, v in temp_map.items():
                if A.nodes[k]['elem'] != B.nodes[v]['elem']:
                    skip = True
                    break
            if skip:
                continue

            total_temp_map = copy.copy(temp_map)
            total_temp_map.update(self.assisting_map)

            total_temp_map_inverted = {v: k for k, v in total_temp_map.items()}

            B_copy = nx.Graph(B)
            B_copy = nx.relabel_nodes(B_copy, total_temp_map_inverted)

            breaking_edges = A.edges - B_copy.edges
            forming_edges = B_copy.edges - A.edges
            # If the forced breaking or forced forming edges are not respected, skip this mapping
            if (forced_breaking_edges.issubset(forming_edges)
                    and len(forced_breaking_edges)
                    > 0) or (forced_forming_edges.issubset(breaking_edges)
                             and len(forced_forming_edges) > 0):
                continue

            changing_bonds = breaking_edges | forming_edges
            n_changing_H_bonds = sum([
                1 for u, v in changing_bonds
                if A.nodes[u]['elem'] == 1.0 or A.nodes[v]['elem'] == 1.0
            ])
            if self.verbose:
                self.ostream.print_info(
                    f"Permutation {self._print_mapping(temp_map)} has a total of {len(changing_bonds)} changing bonds, forming: {self._print_bond_list(forming_edges)}, breaking: {self._print_bond_list(breaking_edges)} of which {n_changing_H_bonds} involve hydrogens."
                )
                self.ostream.flush()

            maps.append(copy.copy(total_temp_map))
            breaking_edges_list.append(breaking_edges)
            forming_edges_list.append(forming_edges)
            N_changing_bonds.append(len(changing_bonds))
            N_changing_H_bonds.append(n_changing_H_bonds)

        # The best mapping is the one with the least amount of changing bonds
        # If there is a tie, the tie-breaker is the one with the most changing H-bonds
        # If there is still a tie, pick the first one
        least_changing_bonds_indices = np.where(
            N_changing_bonds == np.array(N_changing_bonds).min())[0]
        if len(least_changing_bonds_indices) > 1:
            self.ostream.print_info(
                "Multiple mappings with same least amount of changing bonds found, picking one the most changing H-bonds."
            )
            least_changing_bonds_H_bond_counts = [
                N_changing_H_bonds[i] for i in least_changing_bonds_indices
            ]
            best_H_bond_count = np.array(
                least_changing_bonds_H_bond_counts).max()
            least_changing_bonds_index = np.where(
                least_changing_bonds_H_bond_counts == best_H_bond_count)[0]

            best_index = least_changing_bonds_indices[
                least_changing_bonds_index[0]]
            if len(least_changing_bonds_index):
                self.ostream.print_info(
                    f"{len(least_changing_bonds_index)} equally good mappings found, picking first one."
                )
                for i in least_changing_bonds_index:
                    temp_i = least_changing_bonds_indices[i]
                    self.ostream.print_info(
                        f"Mapping {maps[temp_i]} with breaking bonds: {breaking_edges_list[temp_i]} and forming bonds: {forming_edges_list[temp_i]}"
                    )
                self.ostream.flush()
        else:
            best_index = least_changing_bonds_indices[0]
        self.ostream.print_info(
            f"Picking mapping {maps[best_index]} with breaking bonds: {breaking_edges_list[best_index]} and forming bonds: {forming_edges_list[best_index]}"
        )
        self.ostream.flush()
        to_return_map = maps[best_index]
        if swapped:
            # if we swapped A and B, we need to swap the mapping back
            self.ostream.print_info(
                "Inverting mapping because A and B were swapped.")
            to_return_map = {v: k for k, v in to_return_map.items()}
        return to_return_map, breaking_edges_list[
            best_index], forming_edges_list[best_index]

    def _find_breaking_edges(self, A, A_assist, B, B_assist,
                             forced_forming_edges):
        breaking_edges = None
        B_composition = self._get_bond_element_composition(B)
        for depth in range(0, self._breaking_depth + 1):
            self._check_monomorphic = False
            self.ostream.print_info(
                f"Finding isomorphic breaking edges at depth {depth}")
            self.ostream.flush()
            breaking_edges = self._recurse_breaking_edges(
                A,
                A_assist,
                B,
                B_assist,
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
                A_assist,
                B,
                B_assist,
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

    def _recurse_breaking_edges(self, A, A_assist, B, B_assist, depth,
                                B_composition, forced_forming_edges):
        """
        Search the minimal amount of edges that need to be removed from A so that all connected components of A are subgraphs of B.
        Searches recursively up to given depth. Returns None if no edges are found at a given depth.
        In looping, it prioritizes breaking bonds involving combinations of elements that are more abundant in A then in B.
        After that, it prioritizes breaking bonds involving Hydrogen atoms.
        """

        self._check_time()

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

                if A_assist.has_edge(*edge) and B_assist.has_edge(*edge):
                    continue

                A.remove_edge(*edge)
                if self.verbose:
                    self.ostream.print_info(
                        f"Trying to find breaking edge {self._print_bond(edge)} at depth {depth}"
                    )
                    self.ostream.flush()
                edges = self._recurse_breaking_edges(
                    A,
                    A_assist,
                    B,
                    B_assist,
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

    def _find_forming_edges(self, A, A_assist, B, B_assist, breaking_edges):
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
                    edge_in_A = A.has_edge(node_i, node_j)
                    # If an edge was broken, it should not be re-formed
                    edge_in_broken = (node_i, node_j) in breaking_edges or (
                        node_j, node_i) in breaking_edges

                    # HH_bond = A.nodes[node_i]['elem'] == 1.0 and A.nodes[node_j]['elem'] == 1.0
                    if edge_in_A or edge_in_broken:
                        continue

                    comb = self._get_elem_comb(node_i, node_j, A)
                    elem_comb_A = elem_combinations_A.get(comb, 0)
                    elem_comb_B = elem_combinations_B.get(comb, 0)
                    if elem_comb_B <= elem_comb_A:
                        continue

                    if node_i in A_assist.nodes and node_j in A_assist.nodes:
                        if node_i in B_assist.nodes and node_j in B_assist.nodes:
                            if not A_assist.has_edge(
                                    node_i, node_j) and not B_assist.has_edge(
                                        node_i, node_j):
                                continue

                    A.add_edge(node_i, node_j)
                    if self.verbose:
                        self.ostream.print_info(
                            f"Trying to find forming bond {bond_count}: {self._print_bond((node_i, node_j))}"
                        )
                        self.ostream.flush()
                    # self.ostream.print_info(f"Checking bond {self._print_bond((node_i, node_j))}")
                    if not self._connected_components_are_subgraphs(A, B):
                        A.remove_edge(node_i, node_j)
                        continue
                    self._check_time()

                    edge = (node_i, node_j)
                    forming_edges.append(edge)
                    self._print_elapsed(
                        f"finding bond {bond_count}: {self._print_bond(edge)}")
                    self.ostream.flush()
                    break
                if edge is not None:
                    break
            if edge is None:
                self.ostream.print_info("Bond not found, aborting")
                self.ostream.flush()
                return None
        return forming_edges

    @staticmethod
    def _get_bond_element_composition(A):
        gen = (ReactionMatcher._get_elem_comb(u, v, A) for u, v in A.edges())
        composition = Counter(gen)
        return dict(composition)

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
            GM = self.get_graph_matcher(B, A_sub)
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
            GM = self.get_graph_matcher(B, A_sub)
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

    def get_graph_matcher(self, A, B):

        def node_match(a, b):
            # The graph matcher calls this for every candidate pair, which makes it the place to stop a long search
            self._check_time()

            # the elem key is always expected, the assist_id key and H_bond_count key are optional
            elem_match = a['elem'] == b['elem']
            h_count_match = a.get('H_bond_count', 0) == b.get('H_bond_count', 0)
            assist_match = a.get('assist_id', -1) == b.get('assist_id', -1)
            if self._reduce_hydrogen:
                return elem_match and assist_match and h_count_match
            else:
                return elem_match and assist_match

        GM = GraphMatcher(A, B, node_match=node_match)
        return GM

    def _restore_hydrogen(self, map, rea_graph, pro_graph):
        for rea_node, rea_info in rea_graph.nodes.items():
            if rea_info.get('H_bond_count', 0) > 0:
                rea_H = rea_info['H_indices']
                pro_node = map[rea_node]
                pro_info = pro_graph.nodes[pro_node]
                pro_H = pro_info['H_indices']

                for h1, h2 in zip(rea_H, pro_H):
                    map.update({h1: h2})
                self.ostream.flush()

    def _check_time(self):
        """
        Raise _TimeLimitExceeded once the deadline of get_mapping has passed.
        """
        if self._time_limit_active and time.monotonic() > self._deadline:
            raise _TimeLimitExceeded()

    def _print_elapsed(self, msg):
        self.ostream.print_info(
            f"Total time = {time.monotonic() - self._start_time:.3f} s after {msg}"
        )
        self.ostream.flush()

    def _print_mapping(self, mapping):
        items = list(mapping.items())
        items.sort()
        return "{" + ", ".join([
            f"{k+self.print_starting_index}: {v+self.print_starting_index}"
            for k, v in items
        ]) + "}"

    def _print_bond_list(self, bonds):
        bonds = list(bonds)
        return f"[{', '.join([self._print_bond(bond) for bond in bonds])}]"

    def _print_bond(self, bond):
        return f"({bond[0]+self.print_starting_index}, {bond[1]+self.print_starting_index})"
