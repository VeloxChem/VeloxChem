from mpi4py import MPI
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.algorithms.isomorphism import categorical_node_match
import networkx as nx
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

        self.max_time = 600

    def match_reaction_graphs(self, A: nx.Graph, B: nx.Graph):
        # figure out how many bonds need to change -> difficult to impossible with multiple graphs right?
        largest_in_A = len(list(nx.connected_components(A))[0])
        largest_in_B = len(list(nx.connected_components(B))[0])
        swapped = False
        if largest_in_A > largest_in_B:
            A, B = B, A
            swapped = True
        
        mapping, forming_bonds, breaking_bonds = ReactionMatcher._match_subgraph(A, B)
        start_time = time.time()

        if not mapping == {}:
            spent_time = time.time() - start_time
            self.ostream.print_info(f"Found mapping without removing bonds.")
        else:
            # self.ostream.print_info("No mapping found without removing bonds, trying to remove one bond")
            Am1_mappings = []
            for edge in A.edges:
                Am1 = nx.Graph(A)
                Am1.remove_edge(*edge)
                bond_elems = (A.nodes[edge[0]]["elem"], A.nodes[edge[1]]["elem"])

                mapping, forming_bonds, breaking_bonds = ReactionMatcher._match_subgraph(Am1, B)
                if mapping != {}:
                    Am1_mappings.append({"mapping": mapping, "forming_bonds": forming_bonds, "breaking_bonds": breaking_bonds})

            if len(Am1_mappings) > 0:
                self.ostream.print_info(f"Found {len(Am1_mappings)} mappings with removing 1 bond")
                mapping = self._sort_broken_mappings(Am1_mappings)
            else:
                Am2_mappings = []
                remaining_edges = set(A.edges)
                for edge1 in A.edges:
                    if time.time() - start_time > self.max_time:
                        self.ostream.print_info(f"Spent more then {self.max_time} seconds finding mapping, aborting")
                        mapping = {}
                        break
                    edge_set = set()
                    edge_set.add(edge1)
                    remaining_edges = remaining_edges - edge_set
                    Am1 = nx.Graph(A)
                    Am1.remove_edge(*edge1)
                    for edge2 in remaining_edges:
                        Am2 = nx.Graph(Am1)
                        Am2.remove_edge(*edge2)
                        bond_elems1 = {A.nodes[edge1[0]]["elem"], A.nodes[edge1[1]]["elem"]}
                        bond_elems2 = {A.nodes[edge2[0]]["elem"], A.nodes[edge2[1]]["elem"]}
                        mapping, forming_bonds, breaking_bonds = ReactionMatcher._match_subgraph(Am2, B)
                        if mapping != {}:
                            Am2_mappings.append({"mapping": mapping, "forming_bonds": forming_bonds, "breaking_bonds": breaking_bonds})
                if len(Am2_mappings) > 0:
                    mapping = self._sort_broken_mappings(Am2_mappings)
                    self.ostream.print_info(f"Found {len(Am2_mappings)} mappings with removing 2 bonds")
        if mapping == {}:
            self.ostream.print_info("No mapping found with removing two bonds, removing 3 bonds is not implemented. Try suggesting some broken bonds.")
        if swapped:
            mapping = {v: k for k, v in mapping.items()}

        spent_time = time.time() - start_time
        self.ostream.print_info(f"Spent {spent_time:.3f} seconds finding mapping")
        self.ostream.flush()
        return mapping

    def _sort_broken_mappings(self, mappings):
        sorted_mappings = sorted(mappings, key=lambda x: x["forming_bonds"]+x["breaking_bonds"])
        total_mappings = len(sorted_mappings)
        least_changing_bonds = mappings[0]["forming_bonds"]+mappings[0]["breaking_bonds"]
        best_mappnigs = [mapping for mapping in mappings if mapping["forming_bonds"]+mapping["breaking_bonds"] == least_changing_bonds]
        
        self.ostream.print_info(f"Found {total_mappings} mappings with removing 2 bonds, of which {len(best_mappnigs)} have {least_changing_bonds}+2 changing bonds which is the minimal amount")
        return mappings[0]["mapping"]
        
    @staticmethod
    def _match_subgraph(A, B):
        A = nx.Graph(A)
        B = nx.Graph(B)
        A_copy = nx.Graph(A)
        B_copy = nx.Graph(B)
        total_mapping = {}
        while A.number_of_nodes() > 0:
            # if find largest subgraph
            mapping, res, size = ReactionMatcher._find_largest_subgraph(A, B)
            if mapping is None:
                return {}, -1, -1
            # remove mapping from both graphs
            
            for key, value in mapping.items():
                A.remove_node(key)
                B.remove_node(value)
            
            total_mapping.update(mapping)
    
        forming_bonds, breaking_bonds = ReactionMatcher._count_changing_bonds(A_copy, B_copy, total_mapping)
        return total_mapping, forming_bonds, breaking_bonds

    @staticmethod
    def _count_changing_bonds(A, B, mapping):
        # apply mapping to B
        B_mapped = nx.Graph(B)
        mapping = {v: k for k, v in mapping.items()}
        B_mapped = nx.relabel_nodes(B_mapped, mapping)
        
        A_bonds = set(tuple(sorted(bond)) for bond in A.edges)
        B_bonds = set(tuple(sorted(bond)) for bond in B_mapped.edges)
        breaking_bonds = A_bonds - B_bonds
        forming_bonds = B_bonds - A_bonds
        return len(forming_bonds), len(breaking_bonds)
        
    @staticmethod
    def _find_largest_subgraph(a: nx.Graph, b: nx.Graph):
        swapped = False
        # Sort the molecules by size
        A = ReactionMatcher._sort_graph_by_size(ReactionMatcher._split_graphs(a))
        B = ReactionMatcher._sort_graph_by_size(ReactionMatcher._split_graphs(b))

        # If the largest graph is in B, swap A and B, the reactant is now being treated as the product or vice versa
        if len(A[0].nodes) < len(B[0].nodes):
            A, B = B, A

            if swapped is False:
                swapped = True

        # Find the next largest subgraph isomorphism of the elements of B in A[0]
        # move swapping into here
        mapping, res = ReactionMatcher._find_next_subgraph(A[0], B)

        # Save the obtained mapping
        # If the reactant is being treated as the product, the mapping needs to be inverted before saving
        if mapping is not None:
            if swapped:
                mapping = {v: k for k, v in mapping.items()}

            # self.ostream.print_info(f"Found mapping through largest subgraph: {mapping}")
            size = len(mapping)
        else:
            size = None

        return mapping, res, size

    # Split a graph into a list of connected graphs
    @staticmethod
    def _split_graphs(graph: nx.Graph | list[nx.Graph]) -> list[nx.Graph]:
        graphs = []
        if type(graph) is nx.Graph:
            graph = [graph]
        for g in graph:
            graphs.extend(list(g.subgraph(c) for c in nx.connected_components(g)))
        return graphs

    
    @staticmethod
    def _sort_graph_by_size(graphs: list[nx.Graph]) -> list[nx.Graph]:
        return sorted(graphs, key=lambda graph: len(graph.nodes), reverse=True)

    @staticmethod
    def _find_next_subgraph(a: nx.Graph, B: list[nx.Graph]):
        # Loops through B and finds the largest subgraph isomorphism in A that leaves the least scattered atoms
        # Returns the mapping and the index of the subgraph in B
        B = sorted(B, key=lambda graph: len(graph.nodes), reverse=True)

        for b in B:
            GM = GraphMatcher(a, b, node_match=categorical_node_match("elem", 0))

            # Get all subgraph isomorphisms
            sub_graph_mappings = [sub for sub in GM.subgraph_isomorphisms_iter()]

            # If there are any subgraph isomorphisms, find the one with the smallest number of connected components in the remaining graph
            if len(sub_graph_mappings) > 0:
                best_mapping = None
                best_res = -1

                for mapping in sub_graph_mappings:

                    temp_a = nx.Graph(a)  # Unfreeze the graph to allow for node removal
                    temp_a.remove_nodes_from(mapping.keys())
                    res = nx.number_connected_components(temp_a)
                    if best_res == -1 or res < best_res:
                        best_res = res
                        best_mapping = mapping

                return best_mapping, best_res
        return None, None
