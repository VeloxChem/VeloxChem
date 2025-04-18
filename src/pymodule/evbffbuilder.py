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
import numpy as np
import networkx as nx
import sys

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestdriver import ScfRestrictedDriver
from .respchargesdriver import RespChargesDriver
from .xtbdriver import XtbDriver
from .xtbhessiandriver import XtbHessianDriver
from .optimizationdriver import OptimizationDriver
from .mmforcefieldgenerator import MMForceFieldGenerator
from .reactionmatcher import ReactionMatcher
from .outputstream import OutputStream


class EvbForceFieldBuilder():

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

        self.optimize: bool = False
        self.reparameterize: bool = True

        self.input_folder: str = "input_files"

        self.reactant: MMForceFieldGenerator = None
        self.products: list[MMForceFieldGenerator] = None

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
            self.reparameterize,
            self.optimize,
        )
        self.reactant.ostream.flush()

        products: list[MMForceFieldGenerator] = []

        

        for i, input in enumerate(product_input):
            products.append(
                self.get_forcefield(
                    input,
                    product_charge[i],
                    product_multiplicity[i],  # type:ignore
                    self.reparameterize,
                    self.optimize
                ))

        rea_elems = self.reactant.molecule.get_element_ids()
        pro_elems = [
            element_id
            for pro_ff in products
            for element_id in pro_ff.molecule.get_element_ids()
        ]

        # Never merge the reactant forcefield generators, for these we actually need positions
        self.ostream.print_info("Creating combined product force field")
        self.ostream.flush()
        self.product = self._create_combined_forcefield(products)
        
        if not ordered_input:
            self.ostream.print_info("Matching reactant and product force fields")
            self.ostream.flush()
            self.product = self._match_reactant_and_product(self.reactant, rea_elems,self.product, pro_elems, breaking_bonds)
        self.product.ostream.flush()
        self._summarise_reaction(self.reactant, self.product)
        
        return self.reactant, self.product

    def get_forcefield(
        self,
        input: dict,
        charge: int,
        multiplicity: int,
        reparameterize: bool,
        optimize: bool,
    ) -> MMForceFieldGenerator:

        molecule = input["molecule"]
        molecule.set_multiplicity(multiplicity)
        molecule.set_charge(charge)

        # # If charges exist, load them into the forcefield object before creating the topology

        # The topology creation will calculate charges if they're not already set

        if input["forcefield"] is not None:
            forcefield = input["forcefield"]
            forcefield.molecule = molecule
        else:
            if input["optimize"] and optimize:
                self.ostream.print_info("Optimising the geometry with xtb.")
                scf_drv = XtbDriver()
                opt_drv = OptimizationDriver(scf_drv)
                opt_drv.hessian = "last"
                opt_results = opt_drv.compute(molecule)
                molecule = Molecule.from_xyz_string(opt_results["final_geometry"])

            forcefield = MMForceFieldGenerator()

            #Load or calculate the charges
            
            if input["charges"] is not None:
                assert len(input["charges"]) == molecule.number_of_atoms(), "The number of provided charges does not match the number of atoms in the molecule"
                charge_sum = sum(input["charges"])
                if charge_sum - round(charge_sum) > 0.001:
                    self.ostream.print_warning(f"Sum of charges is {charge_sum} is not close to an integer. Confirm that the input is correct.")  
                forcefield.partial_charges = input["charges"]
                self.ostream.print_info("Creating topology")
                self.ostream.flush()
                forcefield.create_topology(molecule)
            else:
                if max(molecule.get_masses()) > 84:
                    basis = MolecularBasis.read(molecule, "STO-6G", ostream=None)
                    self.ostream.print_info(
                        f"Heavy ({max(molecule.get_masses())}) atom found. Using STO-6G basis (only comes in for RESP calculation)."
                    )
                else:
                    basis = MolecularBasis.read(molecule, "6-31G*", ostream=None)
                if multiplicity == 1:
                    scf_drv = ScfRestrictedDriver()
                else:
                    scf_drv = ScfUnrestrictedDriver()
                self.ostream.flush()
                scf_results = scf_drv.compute(molecule, basis)
                if not scf_drv.is_converged:
                    scf_drv.conv_thresh = 1.0e-4
                    scf_drv.max_iter = 200
                    scf_results = scf_drv.compute(molecule, basis)
                assert scf_drv.is_converged, f"SCF calculation for RESP charges did not converge, aborting"
                
                resp_drv = RespChargesDriver()
                self.ostream.print_info("Calculating RESP charges")
                self.ostream.flush()
                forcefield.partial_charges = resp_drv.compute(molecule, basis,scf_results,'resp')
                self.ostream.flush()
                self.ostream.print_info("Creating topology")
                forcefield.create_topology(molecule, basis, scf_results=scf_results)

            # The atomtypeidentifier returns water with no Lennard-Jones on the hydrogens, which leads to unstable simulations
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

            #Reparameterize the forcefield if necessary and requested
            if reparameterize:
                self.ostream.print_info("Reparameterising force field.")
                
                if input["hessian"] is not None:
                    hessian = input["hessian"]
                else:
                    xtb_drv = XtbDriver()
                    xtb_hessian_drv = XtbHessianDriver(xtb_drv)
                    self.ostream.flush()
                    xtb_hessian_drv.compute(molecule)
                    hessian = np.copy(xtb_hessian_drv.hessian)  # type: ignore
                self.ostream.flush()
                forcefield.reparameterize(hessian=hessian)
        return forcefield

    #Match the indices of the reactant and product forcefield generators
    def _match_reactant_and_product(
        self,
        reactant_ff: MMForceFieldGenerator,
        rea_elems: list,
        product_ff: MMForceFieldGenerator,
        pro_elems: list,
        breaking_bonds: list[tuple[int, int]] | None = None,
    ) -> MMForceFieldGenerator:
        assert len(reactant_ff.atoms) == len(product_ff.atoms), "The number of atoms in the reactant and product do not match"
        # Turn the reactand and product into graphs
        rea_graph = nx.Graph()
        reactant_bonds = list(reactant_ff.bonds.keys())
        # Remove the bonds that are being broken, so that these segments get treated as seperate reactants
        if breaking_bonds is not None:
            reactant_bonds = [bond for bond in reactant_bonds if bond not in breaking_bonds]

        rea_graph.add_nodes_from(reactant_ff.atoms.keys())
        rea_graph.add_edges_from(reactant_bonds)

        for i, elem in enumerate(rea_elems):
            rea_graph.nodes[i]['elem'] = elem

        pro_graph = nx.Graph()
        pro_graph.add_nodes_from(product_ff.atoms.keys())
        pro_graph.add_edges_from(list(product_ff.bonds.keys()))
        for i, elem in enumerate(pro_elems):
            pro_graph.nodes[i]['elem'] = elem
            
        rm = ReactionMatcher()
        total_mapping = rm.match_reaction_graphs(rea_graph, pro_graph)
        total_mapping = {v: k for k, v in total_mapping.items()}
        self.ostream.print_info(f"Mapping: {total_mapping}")
        EvbForceFieldBuilder._apply_mapping_to_forcefield(product_ff, total_mapping)
        return product_ff

        # Merge a list of forcefield generators into a single forcefield generator while taking care of the atom indices
    @staticmethod
    def _create_combined_forcefield(forcefields: list[MMForceFieldGenerator]) -> MMForceFieldGenerator:
        forcefield = MMForceFieldGenerator()
        forcefield.atoms = {}
        forcefield.bonds = {}
        forcefield.angles = {}
        forcefield.dihedrals = {}
        forcefield.impropers = {}
        atom_count = 0
        
        for product_ffgen in forcefields:
            # Shift all atom keys by the current atom count so that every atom has a unique ID
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
    def _apply_mapping_to_forcefield(forcefield: MMForceFieldGenerator, mapping: dict[int, int]) -> MMForceFieldGenerator:
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
        self.ostream.print_info(f"{len(broken_bonds)} breaking bonds:")
        if len(broken_bonds) > 0:
            self.ostream.print_info("ReaType, ProType, ID - ReaType, ProType, ID")
        for bond_key in broken_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0]
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1]
            self.ostream.print_info(f"{reactant_type0:<9}{product_type0:<9}{id0:<2} - {reactant_type1:<9}{product_type1:<9}{id1:<2}")

        self.ostream.print_info(f"{len(formed_bonds)} forming bonds:")
        if len(formed_bonds) > 0:
            self.ostream.print_info("ReaType, ProType, ID\t - ReaType, ProType, ID")
        for bond_key in formed_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0]
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1]
            self.ostream.print_info(f"{reactant_type0:<9}{product_type0:<9}{id0:<2} - {reactant_type1:<9}{product_type1:<9}{id1:<2}")

        return formed_bonds, broken_bonds
        
