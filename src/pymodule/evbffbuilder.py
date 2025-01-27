import os

import shutil

from .molecule import Molecule
from .molecularbasis import MolecularBasis

from .scfunrestdriver import ScfUnrestrictedDriver
from .respchargesdriver import RespChargesDriver
from .xtbdriver import XtbDriver
from .xtbhessiandriver import XtbHessianDriver
from .optimizationdriver import OptimizationDriver
from .forcefieldgenerator import ForceFieldGenerator

import networkx as nx
import numpy as np
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.algorithms.isomorphism import categorical_node_match


class EvbForceFieldBuilder():

    def __init__(self):
        self.optimise: bool = False
        self.reparameterise: bool = True

        self.input_folder: str = "input_files"

        self.reactant: ForceFieldGenerator
        self.products: list[ForceFieldGenerator]
        self.gaff_path = None
        pass

    def build_forcefields(
        self,
        reactant_input: dict,
        product_input: list[dict],
        reactant_charge: int = 0,
        product_charge: list[int] = [0],
        reactant_multiplicity: int = 1,
        product_multiplicity: list[int] = [1],
        ordered_input: bool = False,
        breaking_bonds: list[tuple[int, int]] | None = None,
    ):

        self.reactant = self.get_forcefield(
            reactant_input,
            reactant_charge,
            reactant_multiplicity,
            self.reparameterise,
            self.optimise,
        )
        self.reactant.ostream.flush()

        products: list[ForceFieldGenerator] = []

        for i, input in enumerate(product_input):
            products.append(
                self.get_forcefield(
                    input,
                    product_charge[i],
                    product_multiplicity[i],  # type:ignore
                    self.reparameterise,
                    self.optimise,
                ))

        #todo all this stuff should get removed, at least moved to a higher class, don't want to be managing files here
        # Get all files and directories in the current directory
        # items = os.listdir(".")

        # for item in items:
        #     # If the item starts with 'vlx_'
        #     if item.startswith("vlx_"):
        #         # Construct full item path
        #         item_path = os.path.join(".", item)
        #         # If it's a file, remove it
        #         if os.path.isfile(item_path):
        #             os.remove(item_path)
        #         # If it's a directory, remove it including all its content
        #         elif os.path.isdir(item_path):
        #             shutil.rmtree(item_path)

        # Never merge the reactant forcefield generators, for these we actually need positions
        self.product= self._create_combined_forcefield(products)
        
        if not ordered_input:
            self.product = self._match_reactant_and_product(self.reactant, self.product, breaking_bonds)
        self.product.ostream.flush()
        self._summarise_reaction(self.reactant, self.product)
        
        return self.reactant, self.product

    def get_forcefield(
        self,
        input: dict,
        charge: int,
        multiplicity: int,
        reparameterise: bool,
        optimise: bool,
    ) -> ForceFieldGenerator:

        molecule = input["molecule"]
        molecule.set_multiplicity(multiplicity)
        molecule.set_charge(charge)

        # # If charges exist, load them into the forcefield object before creating the topology

        # The topology creation will calculate charges if they're not already set

        if input["forcefield"] is not None:
            forcefield = input["forcefield"]
            forcefield.molecule = molecule
        else:
            if input["optimise"] and optimise:
                print("Optimising the geometry with xtb.")
                scf_drv = XtbDriver()
                opt_drv = OptimizationDriver(scf_drv)
                opt_drv.hessian = "last"
                opt_results = opt_drv.compute(molecule)
                molecule = Molecule.from_xyz_string(opt_results["final_geometry"])

            forcefield = ForceFieldGenerator()
            if self.gaff_path is not None:
                forcefield.force_field_data = self.gaff_path

            #Load or calculate the charges
            
            if input["charges"] is not None:
                forcefield.partial_charges = input["charges"]
                print("Creating topology")
                forcefield.create_topology(molecule)
            else:
                if max(molecule.get_masses()) > 84:
                    basis = MolecularBasis.read(molecule, "STO-6G", ostream=None)
                    print(
                        f"Heavy ({max(molecule.get_masses())}) atom found. Using STO-6G basis (only comes in for RESP calculation)."
                    )
                else:
                    basis = MolecularBasis.read(molecule, "6-31G*", ostream=None)
                scf_drv = ScfUnrestrictedDriver()
                scf_results = scf_drv.compute(molecule, basis)
                if not scf_drv.is_converged:
                    scf_drv.conv_thresh = 1.0e-4
                    scf_drv.max_iter = 200
                    scf_results = scf_drv.compute(molecule, basis)
                assert scf_drv.is_converged, f"SCF calculation for RESP charges on compound {filename} did not converge, aborting"

                resp_drv = RespChargesDriver()
                print("Calculating RESP charges")
                forcefield.partial_charges = resp_drv.compute(molecule, basis,scf_results,'resp')
                print("Creating topology")
                forcefield.create_topology(molecule, basis, scf_result=scf_results)

            # The atomtypeidentifier returns water with no lj on the hydrogens, this leads to unstable simulations
            for atom in forcefield.atoms.values():
                if atom['type'] == 'ow':
                    sigma = 1.8200 * 2**(-1 / 6) * 2 / 10
                    epsilon = 0.0930 * 4.184
                    atom['type'] = 'oh'
                    atom['sigma'] = sigma
                    atom['epsilon'] = epsilon
                    atom['comment'] = "Reaction-water oxygen"
                elif atom['type'] == 'hw':
                    
                    sigma = 0.3019 * 2**(-1 / 6) * 2 / 10
                    epsilon = 0.0047 * 4.184
                    atom['type'] = 'ho'
                    atom['sigma'] = sigma
                    atom['epsilon'] = epsilon
                    atom['comment'] = "Reaction-water hydrogen"         

            #Reparameterise the forcefield if necessary and requested
            if reparameterise:
                print("Reparameterising force field.")

                if input["hessian"] is not None:
                    hessian = input["hessian"]
                else:
                    xtb_drv = XtbDriver()
                    xtb_hessian_drv = XtbHessianDriver(xtb_drv)
                    xtb_hessian_drv.compute(molecule)
                    hessian = np.copy(xtb_hessian_drv.hessian)  # type: ignore
                forcefield.reparametrize(hessian=hessian)
        return forcefield


    #Match the indices of the reactant and product forcefield generators
    @staticmethod
    def _match_reactant_and_product(
        reactant: ForceFieldGenerator,
        product: ForceFieldGenerator,
        breaking_bonds: list[tuple[int, int]] | None = None,
    ) -> ForceFieldGenerator:
        # Turn the reactand and product into graphs
        rea_graph = nx.Graph()
        reactant_bonds = list(reactant.bonds.keys())
        # Remove the bonds that are being broken, so that these segments get treated as seperate reactants
        if breaking_bonds is not None:
            reactant_bonds = [bond for bond in reactant_bonds if bond not in breaking_bonds]

        rea_graph.add_nodes_from(reactant.atoms.keys())
        rea_graph.add_edges_from(reactant.bonds.keys())

        for atom in reactant.atoms.items():
            rea_graph.nodes[atom[0]]["mass"] = atom[1]["mass"]

        pro_graph = nx.Graph()
        pro_graph.add_nodes_from(product.atoms.keys())
        pro_graph.add_edges_from(list(product.bonds.keys()))
        for atom in product.atoms.items():
            pro_graph.nodes[atom[0]]["mass"] = atom[1]["mass"]

        # Seperate all molecules into seperate graphs
        product_graphs = EvbForceFieldBuilder._split_graphs(pro_graph)
        reactant_graphs = EvbForceFieldBuilder._split_graphs(rea_graph)
        print(f"{len(reactant_graphs)} reactant molecule(s) and {len(product_graphs)} product molecule(s)")

        A = reactant_graphs
        B = product_graphs

        swapped = False
        total_mapping = {}

        while len(A) > 0 and len(B) > 0:
            # Sort the molecules by size
            A = EvbForceFieldBuilder._sort_graph_by_size(EvbForceFieldBuilder._split_graphs(A))
            B = EvbForceFieldBuilder._sort_graph_by_size(EvbForceFieldBuilder._split_graphs(B))

            # If the largest graph is in B, swap A and B, the reactant is now being treated as the product or vice versa
            if len(A[0].nodes) < len(B[0].nodes):
                A, B = B, A

                if swapped is False:
                    swapped = True
                else:
                    swapped = False

            # Find the next largest subgraph isomorphism of the elements of B in A[0]
            new_mapping, mapping_index = EvbForceFieldBuilder._find_next_subgraph(A[0], B)

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
        print(f"Mapping: {total_mapping}")

        EvbForceFieldBuilder._apply_mapping_to_forcefield(product, total_mapping)
        return product

    # Merge a list of forcefield generators into a single forcefield generator while taking care of the atom indices
    @staticmethod
    def _create_combined_forcefield(forcefields: list[ForceFieldGenerator]) -> ForceFieldGenerator:
        forcefield = ForceFieldGenerator()
        forcefield.atoms = {}
        forcefield.bonds = {}
        forcefield.angles = {}
        forcefield.dihedrals = {}
        forcefield.impropers = {}
        atom_count = 0
        
        for product_ffgen in forcefields:
            # Shift all atom keys by the current atom count so that every atom has a unique ID
            # todo make this more robust, check if shifting is necessary
            shift = atom_count
            mapping = {atom_key: atom_key + shift for atom_key in product_ffgen.atoms}
            EvbForceFieldBuilder._apply_mapping_to_forcefield(product_ffgen, mapping)
            atom_count += len(product_ffgen.atoms)
            forcefield.atoms.update(product_ffgen.atoms)
            forcefield.bonds.update(product_ffgen.bonds)
            forcefield.angles.update(product_ffgen.angles)
            forcefield.dihedrals.update(product_ffgen.dihedrals)
            forcefield.impropers.update(product_ffgen.impropers)
            
        return forcefield

    # Split a graph into a list of connected graphs
    @staticmethod
    def _split_graphs(graph: nx.Graph | list[nx.Graph]) -> list[nx.Graph]:
        graphs = []
        if type(graph) is nx.Graph:
            graph = [graph]
        for g in graph:
            graphs.extend(list(g.subgraph(c) for c in nx.connected_components(g)))
        return graphs

    # Remap indices in the forcefield to the new indices
    @staticmethod
    def _apply_mapping_to_forcefield(forcefield: ForceFieldGenerator, mapping: dict[int, int]) -> ForceFieldGenerator:
        new_product_atoms = {}
        for atom_key in forcefield.atoms:
            key = mapping[atom_key]
            val = forcefield.atoms[atom_key]
            new_product_atoms.update({key: val})

        # Sort the atoms by index
        forcefield.atoms = dict(sorted(new_product_atoms.items()))

        forcefield.bonds = EvbForceFieldBuilder._apply_mapping_to_parameters(forcefield.bonds, mapping)
        forcefield.angles = EvbForceFieldBuilder._apply_mapping_to_parameters(forcefield.angles, mapping)
        forcefield.dihedrals = EvbForceFieldBuilder._apply_mapping_to_parameters(forcefield.dihedrals, mapping)
        forcefield.impropers = EvbForceFieldBuilder._apply_mapping_to_parameters(forcefield.impropers, mapping)
        return forcefield

    #Remap the indices in a specific set of parameters
    @staticmethod
    def _apply_mapping_to_parameters(old_parameters: dict[tuple, dict], mapping: dict[int, int]) -> dict[tuple, dict]:
        new_parameters = {}
        for old_key in old_parameters:
            new_key = tuple([mapping[atom_key] for atom_key in old_key])
            val = old_parameters[old_key]

            # Make sure that the new key is still properly ordered§
            if new_key[-1] < new_key[0]:
                new_key = new_key[::-1]
            new_parameters.update({new_key: val})
        return new_parameters

    @staticmethod
    def _sort_graph_by_size(graphs: list[nx.Graph]) -> list[nx.Graph]:
        return sorted(graphs, key=lambda graph: len(graph.nodes), reverse=True)

    # Loops through B and finds the largest subgraph isomorphism in A that leaves the least scattered atoms
    # Returns the mapping and the index of the subgraph in B
    @staticmethod
    def _find_next_subgraph(a: nx.Graph, B: list[nx.Graph]):
        # find largest graph B in other that is subgraph isomorphic
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

        raise Exception("No subgraph isomorphism found")

    def _summarise_reaction(self, reactant, product):
        """
        Summarises the reaction by printing the bonds that are being broken and formed.

        Returns:
            None
        """
        reactant_bonds = set(reactant.bonds)
        product_bonds = set(product.bonds)
        formed_bonds = product_bonds - reactant_bonds
        broken_bonds = reactant_bonds - product_bonds
        print(f"{len(broken_bonds)} breaking bonds:")
        if len(broken_bonds) > 0:
            print("ReaType, ProType, ID - ReaType, ProType, ID")
        for bond_key in broken_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0]
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1]
            print(f"{reactant_type0:<9}{product_type0:<9}{id0:<2} - {reactant_type1:<9}{product_type1:<9}{id1:<2}")

        print(f"{len(formed_bonds)} forming bonds:")
        if len(formed_bonds) > 0:
            print("ReaType, ProType, ID\t - ReaType, ProType, ID")
        for bond_key in formed_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0]
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1]
            print(f"{reactant_type0:<9}{product_type0:<9}{id0:<2} - {reactant_type1:<9}{product_type1:<9}{id1:<2}")
