import networkx as nx
import itertools
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.algorithms.isomorphism import categorical_node_match
import numpy as np

from collections import Counter


class ReactionMatcher:

    @staticmethod
    def match_reaction_graphs(A: nx.Graph, B: nx.Graph):
        # figure out how many bonds need to change -> difficult to impossible with multiple graphs right?
        
        breaking_bonds = 0
        #make all graphs with one bond broken in A, Am1, and one bond broken in B, Bm1
        # Am1 = []
        # for edge in A.edges:
        #     temp = nx.Graph(A)
        #     temp.remove_edge(*edge)
        #     Am1.append(temp)
        # Bm1 = []
        # for edge in B.edges:
        #     temp = nx.Graph(B)
        #     temp.remove_edge(*edge)
        #     Bm1.append(temp)
        #check if from any in Am1 there's a subgraph isomorphism to B, or from Bm1 to A, 
        #pick the one with the minimal residue = unconnected leftover parts
        #then from Am1 to Bm1, again pick the one with 
        total_mapping = {}
        while A.number_of_nodes() > 0:
            # if find largest subgraph
            mapping = ReactionMatcher._find_largest_subgraph(A, B)

            # elif find unique groups
            if mapping is None:
                mapping = ReactionMatcher._map_unique_groups(A, B)
            # else find largest matching chain
            if mapping is None:
                mapping = ReactionMatcher._find_largest_matching_chain(A, B)
            if mapping is None:
                raise Exception("No mapping found. You broke the algorithm. Please send an e-mail with your input structures to bvh@kth.se")
            # remove mapping from both graphs
            
            for key, value in mapping.items():
                A.remove_node(key)
                B.remove_node(value)
            
            total_mapping.update(mapping)
        
        return total_mapping


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
        mapping = ReactionMatcher._find_next_subgraph(A[0], B)

        # Save the obtained mapping
        # If the reactant is being treated as the product, the mapping needs to be inverted before saving
        if mapping is not None:
            if swapped:
                mapping = {v: k for k, v in mapping.items()}

            print(f"Found mapping through largest subgraph: {mapping}")

        return mapping

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

                return best_mapping
        return None

    @staticmethod
    def _map_unique_groups(A: nx.Graph, B: nx.Graph):
        filters = [None, [6] ,[6,1]]
        mapping = {}
        
        for filter in filters:
            # Make a list of all non-hydrogen node id's their own element and the connected elements for both A and B
            
            A_unique, A_all = ReactionMatcher._find_unique_groups(A, filter)
            B_unique, B_all = ReactionMatcher._find_unique_groups(B, filter)
            a_contradictions = []
            b_contradictions = []
            group_ids = []
            pass
            for ida, groupa in A_unique.items():
                for idb, groupb in B_unique.items():
                    if groupa['element'] == groupb['element'] and groupa['connected_elements'] == groupb['connected_elements']:
                        
                        a_ids = groupa["neighbour_ids"]+[ida]
                        b_ids = groupb["neighbour_ids"]+[idb]
                        for a, b in zip(a_ids,b_ids):
                            if (a not in mapping.keys() and a not in a_contradictions and b not in mapping.values() and b not in b_contradictions):
                                mapping.update({a: b})
                            else:
                                if a in mapping.keys():
                                    mapping.pop(a)
                                if b in mapping.values():
                                    mapping = {k: v for k, v in mapping.items() if v != b}
                                if a not in a_contradictions:
                                    a_contradictions.append(a)
                                if b not in b_contradictions:
                                    b_contradictions.append(b)

            if len(mapping) > 0:
                print(f"Found mapping  {mapping} through unique groups with reactant ids {group_ids} and with filter {filter}")
                return mapping
        return None

    @staticmethod
    def _find_unique_groups(A: nx.Graph, filter):
        groups = {}
        for id in A.nodes:
            elem = A.nodes[id]["elem"]
            if elem != 1:

                neighbour_ids = [n for n in A.neighbors(id)]
                neighbour_elements = [A.nodes[n]["elem"] for n in neighbour_ids]

                filtered_ids = []
                filtered_elements = []
                if filter is None:
                    filtered_ids = neighbour_ids
                    filtered_elements = neighbour_elements
                else:
                    for neigh_id, neigh_elem in zip(neighbour_ids, neighbour_elements):
                        if neigh_elem not in filter:
                            filtered_ids.append(neigh_id)
                            filtered_elements.append(neigh_elem)

                if len(filtered_elements) > 0:
                    filtered_elements, filtered_ids = zip(*sorted(zip(filtered_elements, filtered_ids)))

                groups.update({
                    id: {
                        "element": elem,
                        "connected_elements": list(filtered_elements),
                        "neighbour_ids": list(filtered_ids)
                    }
                })

        unique = {}
        #Loop htrough all 
        for i, (id, groupi) in enumerate(groups.items()):
            found_match = False
            for j, groupj in enumerate(groups.values()):
                if not found_match and i != j:
                    if groupi['element'] == groupj['element'] and groupi['connected_elements'] == groupj['connected_elements']:
                        found_match = True
            if not found_match:
                unique.update({id: groupi})
        
        return unique, groups
