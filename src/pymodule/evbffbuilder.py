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
import os
import copy
from .veloxchemlib import mpi_master
from .sanitychecks import molecule_sanity_check
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
from .waterparameters import get_water_parameters

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

        self.mute_scf: bool = True

        self.optimize_ff: bool = True
        self.water_model: str

    def build_forcefields(
            self,
            reactant_input: list[dict],
            product_input: list[dict],
            reactant_total_multiplicity: int,
            product_total_multiplicity: int,
            ordered_input: bool = False,
            breaking_bonds: set[tuple[int, int]] = set(),
    ):

        reactants: list[MMForceFieldGenerator] = []
        for input in reactant_input:
            reactants.append(
                self.get_forcefield(input, self.reparameterize, self.optimize))

        self.ostream.print_info("Creating combined reactant force field")
        self.ostream.flush()
        reamol = self._combine_molecule(
            [rea['molecule'] for rea in reactant_input],
            reactant_total_multiplicity)
        self.ostream.print_info(
            f"Combined reactant with total charge {reamol.get_charge()} and multiplicity {reamol.get_multiplicity()}"
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
            product_total_multiplicity)
        self.ostream.print_info(
            f"Combined product with total charge {promol.get_charge()} and multiplicity {promol.get_multiplicity()}"
        )
        self.product = self._combine_forcefield(products)
        self.product.molecule = promol

        if not ordered_input:
            self.ostream.print_info(
                "Matching reactant and product force fields")
            self.ostream.flush()
            self.product = self._match_reactant_and_product(
                self.reactant, reamol.get_element_ids(), self.product,
                promol.get_element_ids(), breaking_bonds)

        formed_bonds, broken_bonds = self._summarise_reaction(
            self.reactant, self.product)

        if self.optimize_ff:
            self.reactant.molecule = self._optimize_molecule(
                self.reactant.molecule.get_element_ids(),
                self.reactant,
                formed_bonds,
                note='reactant',
            )
            self.product.molecule = self._optimize_molecule(
                self.product.molecule.get_element_ids(),
                self.product,
                broken_bonds,
                note='product',
            )

        return self.reactant, self.product, formed_bonds, broken_bonds, reactants, products

    @staticmethod
    def _combine_molecule(molecules, total_multiplicity):

        combined_molecule = Molecule()
        # pos = []
        charge = 0
        Sm1 = 0
        for mol in molecules:
            charge += mol.get_charge()
            Sm1 += mol.get_multiplicity() - 1
            if combined_molecule.number_of_atoms() > 0:
                max_x = max(combined_molecule.get_coordinates_in_angstrom()[:,
                                                                            0])
                min_x = min(combined_molecule.get_coordinates_in_angstrom()[:,
                                                                            0])
                shift = max_x - min_x + 2
            else:
                shift = 0

            for elem, coord in zip(mol.get_element_ids(),
                                   mol.get_coordinates_in_angstrom()):
                coord[0] += shift
                # pos.append(coord)
                combined_molecule.add_atom(int(elem), Point(coord), 'angstrom')
        combined_molecule.set_charge(charge)
        if total_multiplicity > -1:
            combined_molecule.set_multiplicity(total_multiplicity)
        else:
            combined_molecule.set_multiplicity(Sm1 + 1)

        molecule_sanity_check(combined_molecule)
        return combined_molecule

    def _optimize_molecule(self,
                           elemental_ids,
                           forcefield,
                           changing_bonds,
                           name='MOL',
                           note=None):
        # for i, ff in enumerate(forcefields):
        changing_atoms = []
        for bond in changing_bonds:
            s1 = forcefield.atoms[bond[0]]['sigma']
            s2 = forcefield.atoms[bond[1]]['sigma']
            if bond[0] not in changing_atoms:
                changing_atoms.append(bond[0])
            if bond[1] not in changing_atoms:
                changing_atoms.append(bond[1])
            s = (s1 + s2) / 2
            rmin = s * 2**(1 / 6)
            # s*=0.8
            forcefield.bonds.update(
                {bond: {
                    'equilibrium': rmin,
                    'force_constant': 2500000
                }})
            self.ostream.print_info(
                f"Adding bond {bond} with r {rmin:.3f} for equilibration {note}"
            )

        self.ostream.flush()


        for bond in changing_bonds:
            forcefield.bonds.pop(bond)

        #merging of systems through openmm files is shaky, as it depends on the atom naming working. See atom renaming in combine_forcefield
        #todo find a better way to do this
        forcefield.write_openmm_files(name, name)

        pdb = mmapp.PDBFile(f'{name}.pdb')
        ff = mmapp.ForceField(f'{name}.xml')

        modeller = mmapp.Modeller(pdb.topology, pdb.positions)

        top = modeller.getTopology()
        pos = modeller.getPositions()

        mmsys = ff.createSystem(
            top,
            nonbondedMethod=mmapp.CutoffNonPeriodic,
            nonbondedCutoff=1.0 * mmunit.nanometers,
        )

        nbforce = [
            force for force in mmsys.getForces()
            if isinstance(force, mm.NonbondedForce)
        ][0]
        existing_exceptions = {}
        for i in range(nbforce.getNumExceptions()):
            params = nbforce.getExceptionParameters(i)
            if params[0] in changing_atoms or params[1] in changing_atoms:
                sorted_tuple = (params[0], params[1])
                existing_exceptions.update({sorted_tuple: i})
        added_exceptions = {}
        for i in changing_atoms:
            for j in changing_atoms:
                if i > j:
                    if (i, j) in existing_exceptions.keys():
                        nbforce.setExceptionParameters(
                            existing_exceptions[(i, j)], i, j, 0, 1, 0)
                    elif (j, i) in existing_exceptions.keys():
                        nbforce.setExceptionParameters(
                            existing_exceptions[(j, i)], i, j, 0, 1, 0)
                    else:

                        index = nbforce.addException(i, j, 0, 1, 0)
                        added_exceptions.update({(i, j): index})
        with open(f'{name}_sys.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(mmsys))
        os.unlink(f'{name}.xml')
        os.unlink(f'{name}.pdb')
        os.unlink(f'{name}_sys.xml')
        integrator = mm.VerletIntegrator(0.001)
        sim = mmapp.Simulation(top, mmsys, integrator)
        sim.context.setPositions(pos)
        sim.minimizeEnergy()
        sim.context.setVelocitiesToTemperature(300 * mmunit.kelvin)
        sim.step(100)
        sim.minimizeEnergy()

        state = sim.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(mmunit.angstrom)

        new_molecule = Molecule()
        for i, elem in enumerate(elemental_ids):
            point = Point(pos[i])
            new_molecule.add_atom(int(elem), point, 'angstrom')
        new_molecule.set_charge(forcefield.molecule.get_charge())
        new_molecule.set_multiplicity(forcefield.molecule.get_multiplicity())
        return new_molecule

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
                scf_drv = XtbDriver(ostream=self.ostream)
                opt_drv = OptimizationDriver(scf_drv)
                opt_drv.hessian = "last"
                if self.mute_scf:
                    self.ostream.print_info("Optimising the geometry with xtb.")
                    self.ostream.mute()
                opt_results = opt_drv.compute(molecule)
                self.ostream.unmute()
                molecule = Molecule.from_xyz_string(
                    opt_results["final_geometry"])

            forcefield = MMForceFieldGenerator(ostream=self.ostream)

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
                forcefield.create_topology(molecule,
                                           water_model=self.water_model)
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
                    scf_drv = ScfRestrictedDriver(ostream=self.ostream)
                else:
                    scf_drv = ScfUnrestrictedDriver(ostream=self.ostream)

                if self.mute_scf:
                    self.ostream.print_info("Calculating SCF for RESP charges")
                    self.ostream.mute()
                scf_results = scf_drv.compute(molecule, basis)
                if not scf_drv.is_converged:
                    scf_drv.conv_thresh = 1.0e-4
                    scf_drv.max_iter = 200
                    scf_results = scf_drv.compute(molecule, basis)
                self.ostream.unmute()
                assert scf_drv.is_converged, f"SCF calculation for RESP charges did not converge, aborting"
                resp_drv = RespChargesDriver(ostream=self.ostream)
                self.ostream.flush()
                if self.mute_scf:
                    self.ostream.print_info("Calculating RESP charges")
                    self.ostream.mute()
                forcefield.partial_charges = resp_drv.compute(
                    molecule, basis, scf_results, 'resp')

                self.ostream.unmute()
                self.ostream.flush()
                self.ostream.print_info("Creating topology")
                forcefield.create_topology(
                    molecule,
                    basis,
                    scf_results=scf_results,
                    water_model=self.water_model,
                )

            # The atomtypeidentifier returns water with no Lennard-Jones on the hydrogens, which leads to unstable simulations
            atom_types = [atom['type'] for atom in forcefield.atoms.values()]
            if 'ow' in atom_types and 'hw' in atom_types and len(
                    atom_types) == 3:
                water_model = get_water_parameters()[self.water_model]
                for atom_id, atom in forcefield.atoms.items():
                    forcefield.atoms[atom_id] = copy.copy(
                        water_model[atom['type']])
                    forcefield.atoms[atom_id]['name'] = forcefield.atoms[
                        atom_id]['name'][0] + str(atom_id)
                for bond_id in forcefield.bonds.keys():
                    forcefield.bonds[bond_id] = water_model['bonds']
                for ang_id in forcefield.angles.keys():
                    forcefield.angles[ang_id] = water_model['angles']

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
        breaking_bonds: set[tuple[int, int]],
    ) -> MMForceFieldGenerator:
        assert len(reactant_ff.atoms) == len(
            product_ff.atoms
        ), "The number of atoms in the reactant and product do not match"
        # Turn the reactand and product into graphs

        rm = ReactionMatcher(ostream=self.ostream)
        total_mapping, breaking_bonds, forming_bonds = rm.get_mapping(
            reactant_ff, rea_elems, product_ff, pro_elems, breaking_bonds)
        if total_mapping is None:
            raise ValueError(
                "Could not find a mapping between the reactant and product force fields."
            )
        total_mapping = {v: k for k, v in total_mapping.items()}
        self.ostream.print_info(f"Mapping: {total_mapping}")
        self.ostream.flush()
        product_ff = EvbForceFieldBuilder._apply_mapping_to_forcefield(
            product_ff,
            total_mapping,
        )

        product_ff.molecule = EvbForceFieldBuilder._apply_mapping_to_molecule(
            product_ff.molecule,
            total_mapping,
        )
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
        forcefield.unique_atom_types = []
        forcefield.pairs = {}

        for i, ff in enumerate(forcefields):
            # Shift all atom keys by the current atom count so that every atom has a unique ID
            shift = atom_count
            mapping = {atom_key: atom_key + shift for atom_key in ff.atoms}
            EvbForceFieldBuilder._apply_mapping_to_forcefield(ff, mapping)
            atom_count += len(ff.atoms)
            for atom in ff.atoms.values():
                atom['name'] = f"{atom['name']}{i+1}"
            forcefield.atoms.update(ff.atoms)
            forcefield.bonds.update(ff.bonds)
            forcefield.angles.update(ff.angles)
            forcefield.dihedrals.update(ff.dihedrals)
            forcefield.impropers.update(ff.impropers)
            if hasattr(ff, 'unique_atom_types'):
                forcefield.unique_atom_types += ff.unique_atom_types
            if hasattr(ff, 'pairs'):
                forcefield.pairs.update(ff.pairs)

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
