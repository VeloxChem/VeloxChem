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
from .reaction_matcher import ReactionMatcher

import numpy as np
import networkx as nx
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
        assert len(reactant.atoms) == len(product.atoms), "The number of atoms in the reactant and product do not match"
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

        total_mapping = ReactionMatcher._match_reaction_graphs(rea_graph, pro_graph)
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

        return formed_bonds, broken_bonds
        
