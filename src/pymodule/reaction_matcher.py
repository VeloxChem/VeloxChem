import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.algorithms.isomorphism import categorical_node_match


class ReactionMatcher:

    @staticmethod
    def _match_reaction_graphs(rea_graph: nx.Graph, pro_graph: nx.Graph):
        # while rea_graph.number_of_nodes() > 0:
        #     # if find largest subgraph
        #     mapping = ReactionMatcher._find_largest_subgraph(rea_graph, pro_graph)
        #     # elif find unique groups
        #     if mapping is None:
        #         mapping = ReactionMatcher._find_unique_group(rea_graph, pro_graph)
        #     # else find largest matching chain
        #     if mapping is None:
        #         mapping = ReactionMatcher._find_largest_matching_chain(rea_graph, pro_graph)
        #     if mapping is None:
        #         raise Exception("No mapping found. You broke the algorithm. Please send an e-mail with your input structures to bvh@kth.se")
        #     # remove mapping from both graphs
        #     rea_graph, pro_graph = ReactionMatcher._remove_mapping_from_graphs(rea_graph, pro_graph, mapping)
        #     # add mapping to total mapping

        # Seperate all molecules into seperate graphs
        product_graphs = ReactionMatcher._split_graphs(pro_graph)
        reactant_graphs = ReactionMatcher._split_graphs(rea_graph)
        print(f"{len(reactant_graphs)} reactant molecule(s) and {len(product_graphs)} product molecule(s)")

        A = reactant_graphs
        B = product_graphs

        swapped = False
        total_mapping = {}

        while len(A) > 0 and len(B) > 0:
            # Sort the molecules by size
            A = ReactionMatcher._sort_graph_by_size(ReactionMatcher._split_graphs(A))
            B = ReactionMatcher._sort_graph_by_size(ReactionMatcher._split_graphs(B))

            # If the largest graph is in B, swap A and B, the reactant is now being treated as the product or vice versa
            if len(A[0].nodes) < len(B[0].nodes):
                A, B = B, A

                if swapped is False:
                    swapped = True
                else:
                    swapped = False

            # Find the next largest subgraph isomorphism of the elements of B in A[0]
            # move swapping into here
            new_mapping, mapping_index = ReactionMatcher._find_next_subgraph(A[0], B)

            # Save the obtained mapping
            # If the reactant is being treated as the product, the mapping needs to be inverted before saving
            if swapped:
                total_mapping.update({v: k for k, v in new_mapping.items()})
            else:
                total_mapping.update(new_mapping)

            # Remove the matched subgraph from both A and B

            # B can be removed as an element from the array
            B.pop(mapping_index)
            # Remove the corresponding nodes from A
            A[0] = nx.Graph(A[0])  # Unfreeze the graph to allow for node removal
            A[0].remove_nodes_from(new_mapping.values())
        return total_mapping


    @staticmethod
    def _find_largest_subgraph(A: nx.Graph, B: nx.Graph):
        # Take the largest connected graph in either A or B
        # Go through the other graph for connected graphs by size, and find the largest subgraph isomorphism
        # Return the mapping 
        return {}

    @staticmethod
    def _find_unique_group(A: nx.Graph, B: nx.Graph):
        # Make a list of all non-hydrogen atoms and their connected atoms in both A and B
        # find if there are any groups in A and B that are unique in at least either A or B
        # do the same but excluding carbons -> should give OH groups
        # do the same excluding carbonds and hydrogens -> should give O and N groups if unique
        # for the first category with hits, return all groups as mappings
        return {}

    @staticmethod
    def _find_largest_matching_chain(A: nx.Graph, B: nx.Graph):
        # todo
        return {}

    @staticmethod
    def _remove_mapping_from_graphs(A: nx.Graph, B: nx.Graph, mapping: dict[int, int]):
        return A, B


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

    # Loops through B and finds the largest subgraph isomorphism in A that leaves the least scattered atoms
    # Returns the mapping and the index of the subgraph in B
    @staticmethod
    def _find_next_subgraph(a: nx.Graph, B: list[nx.Graph]):
        B = sorted(B, key=lambda graph: len(graph.nodes), reverse=True)

        for graph_index, b in enumerate(B):
            GM = GraphMatcher(a, b, node_match=categorical_node_match("mass", 0))

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


                inverted_mapping = {v: k for k, v in best_mapping.items()}  # type: ignore
                return inverted_mapping, graph_index
        raise Exception("No subgraph isomorphism found.") #todo remove this error