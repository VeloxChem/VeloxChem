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
from .veloxchemlib import Point

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


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
        self.product: MMForceFieldGenerator = None

        self.optimize_ff: bool = True

    def build_forcefields(
        self,
        reactant_input: list[dict],
        product_input: list[dict],
        ordered_input: bool = False,
        breaking_bonds: list[tuple[int, int]] | None = None,
    ):

        reactants: list[MMForceFieldGenerator] = []
        for input in reactant_input:
            reactants.append(
                self.get_forcefield(input, self.reparameterize, self.optimize))

        self.ostream.print_info("Creating combined reactant force field")
        self.ostream.flush()
        reamol = molecule = self._combine_molecule(
            [rea['molecule'] for rea in reactant_input],
            reactants,
            self.optimize_ff,
        )
        self.reactant = self._combine_forcefield(reactants)
        self.reactant.molecule = reamol

        products: list[MMForceFieldGenerator] = []

        for input in product_input:
            products.append(
                self.get_forcefield(input, self.reparameterize, self.optimize))

        self.ostream.print_info("Creating combined product force field")
        self.ostream.flush()
        promol = self._combine_molecule(
            [pro['molecule'] for pro in product_input],
            products,
            self.optimize_ff,
        )
        self.product = self._combine_forcefield(products)
        self.product.molecule = promol

        if not ordered_input:
            self.ostream.print_info(
                "Matching reactant and product force fields")
            self.ostream.flush()
            self.product = self._match_reactant_and_product(
                self.reactant, reamol.get_element_ids(), self.product, promol.get_element_ids(),
                breaking_bonds)

        formed_bonds, broken_bonds = self._summarise_reaction(
            self.reactant, self.product)

        return self.reactant, self.product, formed_bonds, broken_bonds

    @staticmethod
    def _combine_molecule(molecules,
                          forcefields,
                          optimise_mm=True,
                          name='MOL',):
        
        combined_molecule = molecules[0]
        # pos = []
        for mol in molecules[1:]:

            max_x = max(combined_molecule.get_coordinates_in_angstrom()[:, 0])
            min_x = min(combined_molecule.get_coordinates_in_angstrom()[:, 0])
            shift = max_x - min_x + 2

            for elem, coord in zip(mol.get_element_ids(),
                                   mol.get_coordinates_in_angstrom()):
                coord[0] += shift
                # pos.append(coord)
                combined_molecule.add_atom(int(elem), Point(coord), 'angstrom')

        # merged_forcefield.molecule = combined_molecule

        if optimise_mm:

            # ffs = []
            ffnames = []
            for i, ff in enumerate(forcefields):
                ffname = f'{name}_{i}'
                ffnames.append(f'{ffname}.xml')
                ff.write_openmm_files(ffname,ffname)
                pdb = mmapp.PDBFile(f'{ffname}.pdb')
                if i == 0:
                    modeller = mmapp.Modeller(pdb.topology,pdb.positions)
                else:
                    modeller.add(pdb.topology, pdb.positions)

            top = modeller.getTopology()
            pos = modeller.getPositions()

            ff = mmapp.ForceField(*ffnames)

            mmsys = ff.createSystem(
                top,
                nonbondedMethod=mmapp.CutoffNonPeriodic,
                nonbondedCutoff=1.0 * mmunit.nanometers,
            )
            integrator = mm.VerletIntegrator(0.001)
            sim = mmapp.Simulation(top, mmsys, integrator)
            sim.context.setPositions(pos)
            sim.minimizeEnergy()

            state = sim.context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True).value_in_unit(mmunit.angstrom)

            elements = combined_molecule.get_element_ids()
            new_molecule = Molecule()
            for i in range(len(elements)):
                point = Point(pos[i])
                new_molecule.add_atom(int(elements[i]), point, 'angstrom')

            combined_molecule = new_molecule
        return combined_molecule

    def get_forcefield(
        self,
        input: dict,
        reparameterize: bool,
        optimize: bool,
    ) -> MMForceFieldGenerator:

        molecule = input["molecule"]

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
                molecule = Molecule.from_xyz_string(
                    opt_results["final_geometry"])

            forcefield = MMForceFieldGenerator()
            forcefield.eq_param = False
            #Load or calculate the charges

            if input["charges"] is not None:
                assert len(input["charges"]) == molecule.number_of_atoms(
                ), "The number of provided charges does not match the number of atoms in the molecule"
                charge_sum = sum(input["charges"])
                if charge_sum - round(charge_sum) > 0.001:
                    self.ostream.print_warning(
                        f"Sum of charges is {charge_sum} is not close to an integer. Confirm that the input is correct."
                    )
                forcefield.partial_charges = input["charges"]
                self.ostream.print_info("Creating topology")
                self.ostream.flush()
                forcefield.create_topology(molecule)
            else:
                if max(molecule.get_masses()) > 84:
                    basis = MolecularBasis.read(molecule,
                                                "STO-6G",
                                                ostream=None)
                    self.ostream.print_info(
                        f"Heavy ({max(molecule.get_masses())}) atom found. Using STO-6G basis (only comes in for RESP calculation)."
                    )
                else:
                    basis = MolecularBasis.read(molecule,
                                                "6-31G*",
                                                ostream=None)
                if molecule.get_multiplicity() == 1:
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
                forcefield.partial_charges = resp_drv.compute(
                    molecule, basis, scf_results, 'resp')
                self.ostream.flush()
                self.ostream.print_info("Creating topology")
                forcefield.create_topology(molecule,
                                           basis,
                                           scf_results=scf_results)

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
            unknowns_params = 'Guessed' in [
                par['comment'] for par in list(forcefield.bonds.values()) +
                list(forcefield.angles.values())
            ]
            if reparameterize and unknowns_params:
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
        assert len(reactant_ff.atoms) == len(
            product_ff.atoms
        ), "The number of atoms in the reactant and product do not match"
        # Turn the reactand and product into graphs
        rea_graph = nx.Graph()
        reactant_bonds = list(reactant_ff.bonds.keys())
        # Remove the bonds that are being broken, so that these segments get treated as seperate reactants
        if breaking_bonds is not None:
            reactant_bonds = [
                bond for bond in reactant_bonds if bond not in breaking_bonds
            ]

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
        self.ostream.flush()
        product_ff = EvbForceFieldBuilder._apply_mapping_to_forcefield(
            product_ff, total_mapping)

        product_ff.molecule = EvbForceFieldBuilder._apply_mapping_to_molecule(
            product_ff.molecule, total_mapping)
        return product_ff

        # Merge a list of forcefield generators into a single forcefield generator while taking care of the atom indices
    @staticmethod
    def _combine_forcefield(
            forcefields: list[MMForceFieldGenerator]) -> MMForceFieldGenerator:
        forcefield = MMForceFieldGenerator()
        forcefield.atoms = {}
        forcefield.bonds = {}
        forcefield.angles = {}
        forcefield.dihedrals = {}
        forcefield.impropers = {}
        atom_count = 0

        for i, ff in enumerate(forcefields):
            # Shift all atom keys by the current atom count so that every atom has a unique ID
            shift = atom_count
            mapping = {
                atom_key: atom_key + shift
                for atom_key in ff.atoms
            }
            EvbForceFieldBuilder._apply_mapping_to_forcefield(
                ff, mapping)
            atom_count += len(ff.atoms)
            # for atom in ff.atoms.values():
            #     atom['name'] = f"{atom['name']}_{i}"
            forcefield.atoms.update(ff.atoms)
            forcefield.bonds.update(ff.bonds)
            forcefield.angles.update(ff.angles)
            forcefield.dihedrals.update(ff.dihedrals)
            forcefield.impropers.update(ff.impropers)

        return forcefield

    # Remap indices in the forcefield to the new indices
    @staticmethod
    def _apply_mapping_to_forcefield(
            forcefield: MMForceFieldGenerator,
            mapping: dict[int, int]) -> MMForceFieldGenerator:
        new_product_atoms = {}
        for atom_key in forcefield.atoms:
            key = mapping[atom_key]
            val = forcefield.atoms[atom_key]
            new_product_atoms.update({key: val})

        # Sort the atoms by index
        forcefield.atoms = dict(sorted(new_product_atoms.items()))

        forcefield.bonds = EvbForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.bonds, mapping)
        forcefield.angles = EvbForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.angles, mapping)
        forcefield.dihedrals = EvbForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.dihedrals, mapping)
        forcefield.impropers = EvbForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.impropers, mapping)
        return forcefield

    @staticmethod
    def _apply_mapping_to_molecule(molecule, mapping):
        new_molecule = Molecule()
        positions = molecule.get_coordinates_in_angstrom()
        element_ids = molecule.get_element_ids()
        sorted_ids = dict(sorted(mapping.items(),
                                 key=lambda item: item[1])).keys()
        for id in sorted_ids:
            new_molecule.add_atom(int(element_ids[id]), Point(positions[id]),
                                  'angstrom')
        return new_molecule
        # int(elem), Point(coord), 'angstrom'

    #Remap the indices in a specific set of parameters
    @staticmethod
    def _apply_mapping_to_parameters(
            old_parameters: dict[tuple, dict],
            mapping: dict[int, int]) -> dict[tuple, dict]:
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
            self.ostream.print_info(
                "ReaType, ProType, ID - ReaType, ProType, ID")
        for bond_key in broken_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0]
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1]
            self.ostream.print_info(
                f"{reactant_type0:<9}{product_type0:<9}{id0:<2} - {reactant_type1:<9}{product_type1:<9}{id1:<2}"
            )

        self.ostream.print_info(f"{len(formed_bonds)} forming bonds:")
        if len(formed_bonds) > 0:
            self.ostream.print_info(
                "ReaType, ProType, ID\t - ReaType, ProType, ID")
        for bond_key in formed_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0]
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1]
            self.ostream.print_info(
                f"{reactant_type0:<9}{product_type0:<9}{id0:<2} - {reactant_type1:<9}{product_type1:<9}{id1:<2}"
            )

        return formed_bonds, broken_bonds
