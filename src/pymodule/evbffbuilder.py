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
import json
import sys
import os
import copy
import math
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
from .errorhandler import assert_msg_critical

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

        self.reactant: MMForceFieldGenerator = None
        self.product: MMForceFieldGenerator = None

        self.optimize_mol: bool = False
        self.reparameterize: bool = True
        self.optimize_ff: bool = True
        self.mm_opt_constrain_bonds: bool = False
        self.water_model: str = 'spce'
        self.reactant_total_multiplicity: int = -1
        self.product_total_multiplicity: int = -1
        self.breaking_bonds: set[tuple[int, int]] | tuple = set()
        self.reactant_partial_charges: list[float] | list[
            list[float]] | None = None
        self.product_partial_charges: list[float] | list[
            list[float]] | None = None
        self.reactant_hessians: np.ndarray | list[np.ndarray
                                                  | None] | None = None
        self.product_hessians: np.ndarray | list[np.ndarray
                                                 | None] | None = None
        # todo what to do with this option?
        self.mute_scf: bool = True

        self.keywords = {
            "optimize_mol": bool,
            "reparameterize": bool,
            "optimize_ff": bool,
            "mm_opt_constrain_bonds": bool,
            "water_model": str,
            "reactant_total_multiplicity": int,
            "product_total_multiplicity": int,
            "breaking_bonds": set | tuple,
            "reactant_partial_charges": list | None,
            "product_partial_charges": list | None,
            "reactant_hessians": np.ndarray | list | None,
            "product_hessians": np.ndarray | list | None,
            "mute_scf": bool
        }

    def read_keywords(self, **kwargs):
        for key, value in kwargs.items():
            if key in self.keywords.keys():
                if isinstance(value, self.keywords[key]):
                    setattr(self, key, value)
                else:
                    raise ValueError(
                        f"Type for given keyword {key} is {type(value)} but should be {self.keywords[key]}"
                    )
            else:
                raise ValueError(
                    f"Unknown keyword {key} in EvbForceFieldBuilder")

    def build_forcefields(
        self,
        reactant: Molecule | list[Molecule],
        product: Molecule | list[Molecule],
    ):

        self.reactant, reactants, reactant_total_charge = self._create_combined_forcefield(
            reactant,
            self.reactant_partial_charges,
            self.reactant_hessians,
            "REA",
        )

        self.product, products, product_total_charge = self._create_combined_forcefield(
            product,
            self.product_partial_charges,
            self.product_hessians,
            "PRO",
        )

        assert reactant_total_charge == product_total_charge, f"Total charge of reactants {reactant_total_charge} and products {product_total_charge} must match"

        self.ostream.print_info("Matching reactant and product force fields")
        self.ostream.flush()

        self.product, product_mapping = self._match_reactant_and_product(
            self.reactant,
            self.reactant.molecule.get_element_ids(),
            self.product,
            self.product.molecule.get_element_ids(),
            self.breaking_bonds,
        )

        formed_bonds, broken_bonds = self._summarise_reaction(
            self.reactant, self.product)

        if self.optimize_ff:
            if len(reactants) > 1:
                self.reactant.molecule = self._optimize_molecule(
                    self.reactant.molecule.get_element_ids(),
                    self.reactant,
                    formed_bonds,
                    note='reactant',
                )
            if len(products) > 1:
                self.product.molecule = self._optimize_molecule(
                    self.product.molecule.get_element_ids(),
                    self.product,
                    broken_bonds,
                    note='product',
                )

        return self.reactant, self.product, formed_bonds, broken_bonds, reactants, products, product_mapping

    def _create_combined_forcefield(
        self,
        molecules: list[Molecule],
        partial_charges: list[list[float]]
        | list[float] | None,
        hessians: np.ndarray | list[np.ndarray | None] | None,
        name: str,
    ):
        if isinstance(molecules, Molecule):
            molecules = [molecules]

        if partial_charges is None:
            partial_charges = [None] * len(molecules)  # type: ignore
        elif isinstance(partial_charges[0], float) or isinstance(
                partial_charges[0], int):
            partial_charges = [partial_charges]  # type: ignore
        assert isinstance(partial_charges, list)

        if isinstance(hessians, np.ndarray) or hessians is None:
            hessians = [hessians] * len(molecules)

        assert_msg_critical(
            len(molecules) == len(partial_charges),
            "Amount of input molecules and lists of partial charges must match")
        assert_msg_critical(
            len(molecules) == len(hessians),
            "Amount of input molecules and lists of hessians must match")

        single_ffs: list[MMForceFieldGenerator] = []
        for i, (molecule, partial_charge,
                hessian) in enumerate(zip(molecules, partial_charges,
                                          hessians)):
            molecule_sanity_check(molecule)
            if partial_charge is not None:
                # Casting to float is necessary for json serialization
                assert isinstance(partial_charge, list)
                partial_charge = [float(x) for x in partial_charge]
                cond = len(partial_charge) == molecule.number_of_atoms()
                msg = f"Number of partial charges {len(partial_charge)} must match the number of atoms {molecule.number_of_atoms()} in the molecule."
                assert_msg_critical(cond, msg)

            self.ostream.print_blank()
            self.ostream.print_header(f"Building force field for {name}_{i+1}")
            self.ostream.print_blank()
            self.ostream.flush()
            ff = self._create_single_forcefield(
                molecule,
                partial_charge,
                hessian,
            )
            single_ffs.append(ff)

        total_charge = sum([mol.get_charge() for mol in molecules])

        self.ostream.print_info("Creating combined reactant force field")
        self.ostream.flush()
        combined_mol = self._combine_molecule(
            molecules,
            self.reactant_total_multiplicity,
        )
        self.ostream.print_info(
            f"Combined reactant with total charge {combined_mol.get_charge()} and multiplicity {combined_mol.get_multiplicity()}"
        )
        self.ostream.flush()
        combined_ff = self._combine_forcefield(single_ffs)
        combined_ff.molecule = combined_mol

        return combined_ff, single_ffs, total_charge

    # todo add hessian option
    # Transforms input difctionary with keys 'molecule', 'charges', 'forcefield', 'optimize' into a forcefield generator
    # If the forcefield is not provided, it will be created from the molecule and charges
    # Calculates resp charges with some fallbacks, does optimization if specified, reparameterizes the forcefield if necessary
    def _create_single_forcefield(
        self,
        molecule,
        partial_charges: list[float] | None,
        hessian: np.ndarray | None = None,
    ) -> MMForceFieldGenerator:

        # # If charges exist, load them into the forcefield object before creating the topology

        # The topology creation will calculate charges if they're not already set
        if self.optimize_mol:
            scf_drv = XtbDriver()
            opt_drv = OptimizationDriver(scf_drv)
            opt_drv.hessian = "last"
            if self.mute_scf:
                self.ostream.print_info("Optimising the geometry with xtb.")
                self.ostream.flush()
                scf_drv.ostream.mute()
            opt_results = opt_drv.compute(molecule)

            molecule = Molecule.from_xyz_string(opt_results["final_geometry"])
            molecule.set_charge(molecule.get_charge())
            molecule.set_multiplicity(molecule.get_multiplicity())

        forcefield = MMForceFieldGenerator(ostream=self.ostream)

        forcefield.eq_param = False
        #Load or calculate the charges

        if partial_charges is not None:
            assert len(partial_charges) == molecule.number_of_atoms(
            ), "The number of provided charges does not match the number of atoms in the molecule"
            charge_sum = sum(partial_charges)
            if charge_sum - round(charge_sum) > 0.001:
                self.ostream.print_warning(
                    f"Sum of charges is {charge_sum} is not close to an integer. Confirm that the input is correct."
                )
            forcefield.partial_charges = partial_charges
            self.ostream.print_info("Creating topology")
            self.ostream.flush()
            forcefield.create_topology(molecule, water_model=self.water_model)
        else:
            if max(molecule.get_masses()) > 84:
                basis = MolecularBasis.read(molecule, "STO-6G", ostream=None)
                self.ostream.print_info(
                    f"Heavy ({max(molecule.get_masses())}) atom found. Using STO-6G basis (only comes in for RESP calculation)."
                )
            else:
                basis = MolecularBasis.read(molecule, "6-31G*", ostream=None)
            if molecule.get_multiplicity() == 1:
                scf_drv = ScfRestrictedDriver()
            else:
                scf_drv = ScfUnrestrictedDriver()

            self.ostream.print_info("Calculating SCF for RESP charges")
            self.ostream.flush()
            if self.mute_scf:
                scf_drv.ostream.mute()
            scf_results = scf_drv.compute(molecule, basis)
            if not scf_drv.is_converged:
                self.ostream.print_warning(
                    "SCF did not converge, increasing convergence threshold to 1.0e-4 and maximum itterations to 200."
                )
                self.ostream.flush()
                scf_drv.conv_thresh = 1.0e-4
                scf_drv.max_iter = 200
                scf_results = scf_drv.compute(molecule, basis)
            # self.ostream.unmute()
            assert scf_drv.is_converged, f"SCF calculation for RESP charges did not converge, aborting"
            resp_drv = RespChargesDriver()
            self.ostream.flush()
            if self.mute_scf:
                resp_drv.ostream.mute()
                self.ostream.print_info("Calculating RESP charges")
                self.ostream.flush()
            forcefield.partial_charges = resp_drv.compute(
                molecule, basis, scf_results, 'resp')
            self.ostream.print_info(
                f"RESP charges: {forcefield.partial_charges}")
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
        if 'ow' in atom_types and 'hw' in atom_types and len(atom_types) == 3:
            water_model = get_water_parameters()[self.water_model]
            for atom_id, atom in forcefield.atoms.items():
                forcefield.atoms[atom_id] = copy.copy(water_model[atom['type']])
                forcefield.atoms[atom_id]['name'] = forcefield.atoms[atom_id][
                    'name'][0] + str(atom_id)
            for bond_id in forcefield.bonds.keys():
                forcefield.bonds[bond_id] = water_model['bonds']
            for ang_id in forcefield.angles.keys():
                forcefield.angles[ang_id] = water_model['angles']

        #Reparameterize the forcefield if necessary and requested
        unknowns_params = 'Guessed' in [
            par['comment'] for par in list(forcefield.bonds.values()) +
            list(forcefield.angles.values())
        ]
        if self.reparameterize and unknowns_params:
            self.ostream.print_info("Reparameterising force field.")
            self.ostream.flush()
            if hessian is None:
                self.ostream.print_info(
                    "Calculating Hessian with XTB to reparameterise the force field."
                )
                self.ostream.flush()
                xtb_drv = XtbDriver()
                xtb_hessian_drv = XtbHessianDriver(xtb_drv)
                xtb_hessian_drv.ostream.mute()
                self.ostream.flush()
                xtb_hessian_drv.compute(molecule)
                hessian = np.copy(xtb_hessian_drv.hessian)  # type: ignore
            else:
                cond = np.shape(hessian) == (molecule.number_of_atoms() * 3,
                                             molecule.number_of_atoms() * 3)
                msg = f"Hessian shape {np.shape(hessian)} should be square with width 3 times the number"
                msg += f" of atoms 3*{molecule.number_of_atoms()}={3*molecule.number_of_atoms()} in the molecule."
                assert_msg_critical(cond, msg)
            self.ostream.flush()
            forcefield.reparameterize(hessian)
        return forcefield

    # Guesses how to combine molecular structures without overlapping them
    def _combine_molecule(self, molecules, total_multiplicity):

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
            self.ostream.print_info(
                f"Setting multiplicity of the combined molecule to {Sm1 + 1} based on the multiplicities of the provided molecules."
            )
            self.ostream.flush()

            combined_molecule.set_multiplicity(Sm1 + 1)

        molecule_sanity_check(combined_molecule)
        return combined_molecule

    def _add_reaction_bonds(self, forcefield, mmsys, changing_bonds, note):
        self.ostream.print_info(
            "Guessing intra-molecular constraints based on breaking bonds. This option can be turned off with mm_opt_constrain_bonds."
        )
        morse_expr = "D*(1-exp(-a*(r-re)))^2 + b*(r-rb)/(1-exp(-k*(r-rb)));"
        morse_force = mm.CustomBondForce(morse_expr)
        morse_force.setName("Reaction morse bond")
        morse_force.addPerBondParameter("D")
        morse_force.addPerBondParameter("a")
        morse_force.addPerBondParameter("re")
        morse_force.addPerBondParameter("b")
        morse_force.addPerBondParameter("rb")
        morse_force.addPerBondParameter("k")

        changing_atoms = []
        for bond in changing_bonds:
            s1 = forcefield.atoms[bond[0]]['sigma']
            s2 = forcefield.atoms[bond[1]]['sigma']
            e1 = forcefield.atoms[bond[0]]['epsilon']
            e2 = forcefield.atoms[bond[1]]['epsilon']
            if bond[0] not in changing_atoms:
                changing_atoms.append(bond[0])
            if bond[1] not in changing_atoms:
                changing_atoms.append(bond[1])
            s = (s1 + s2) / 2
            e = math.sqrt(e1 * e2)

            rmin = s * (2**(1 / 6))  # rmin = sigma * 2^(1/6)
            fc = 250000 * e  # force constant in kJ/mol/angstrom^2
            D = 500
            a = math.sqrt(fc / (2 * D))
            b = 100
            rb = 3 * rmin
            k = 100
            morse_force.addBond(bond[0], bond[1], [D, a, rmin, b, rb, k])
            self.ostream.print_info(
                f"Replacing nonbonded interaction for atoms {bond} with Morse bond  with r {rmin:.3f} and fc {fc:.3f} for MM equilibration of the {note}"
            )
            self.ostream.flush()
        mmsys.addForce(morse_force)

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
        for bond in changing_bonds:
            self.ostream.flush()
            if bond in existing_exceptions.keys():
                nbforce.setExceptionParameters(existing_exceptions[bond],
                                               bond[0], bond[1], 0, 1, 0)
            elif (bond[1], bond[0]) in existing_exceptions.keys():
                nbforce.setExceptionParameters(
                    existing_exceptions[(bond[1], bond[0])], bond[0], bond[1],
                    0, 1, 0)
            else:
                index = nbforce.addException(bond[0], bond[1], 0, 1, 0)
                added_exceptions.update({bond: index})
        return mmsys

    #Match the indices of the reactant and product forcefield generators
    def _match_reactant_and_product(
        self,
        reactant_ff: MMForceFieldGenerator,
        rea_elems: list,
        product_ff: MMForceFieldGenerator,
        pro_elems: list,
        breaking_bonds: set[tuple[int, int]],
    ):
        assert len(reactant_ff.atoms) == len(
            product_ff.atoms
        ), "The number of atoms in the reactant and product do not match"
        # Turn the reactand and product into graphs

        rm = ReactionMatcher(ostream=self.ostream)
        total_mapping, breaking_bonds, forming_bonds = rm.get_mapping(
            reactant_ff,
            rea_elems,
            product_ff,
            pro_elems,
            breaking_bonds,
        )  # type: ignore
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
        return product_ff, total_mapping

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

        for l, ff in enumerate(forcefields):
            # Shift all atom keys by the current atom count so that every atom has a unique ID
            shift = atom_count
            mapping = {atom_key: atom_key + shift for atom_key in ff.atoms}
            EvbForceFieldBuilder._apply_mapping_to_forcefield(ff, mapping)
            atom_count += len(ff.atoms)
            for atom in ff.atoms.values():
                atom['name'] = f"{atom['name']}{l+1}"
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
                "ReaType, ProType, ID - ReaType, ProType, ID")
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

    # Does an FF optimization of the molecule.
    def _optimize_molecule(
        self,
        elemental_ids,
        forcefield,
        changing_bonds,
        name='MOL',
        note=None,
    ):

        #merging of systems through openmm files is shaky, as it depends on the atom naming working. See atom renaming in combine_forcefield
        #todo find a better way to do this
        forcefield.write_openmm_files(name, name)

        # for bond in changing_bonds:
        #     forcefield.bonds.pop(bond)

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
        mmsys_bak = copy.deepcopy(mmsys)
        if self.mm_opt_constrain_bonds:
            mmsys = self._add_reaction_bonds(forcefield, mmsys, changing_bonds,
                                             note)

        with open(f'{name}_sys.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(mmsys))
        os.unlink(f'{name}.xml')
        os.unlink(f'{name}.pdb')
        os.unlink(f'{name}_sys.xml')
        integrator = mm.VerletIntegrator(0.001)
        sim = mmapp.Simulation(top, mmsys, integrator)
        sim.context.setPositions(pos)
        self.ostream.print_info(f"Minimizing {note} molecule.")
        sim.minimizeEnergy()
        sim.context.setVelocitiesToTemperature(600 * mmunit.kelvin)
        self.ostream.print_info(
            f"Running 1000 NVE steps with initial T at 600 K for {note} molecule."
        )
        self.ostream.flush()
        sim.step(1000)
        self.ostream.print_info(f"Minimizing {note} molecule again.")
        self.ostream.flush()
        sim.minimizeEnergy()
        state = sim.context.getState(getPositions=True)
        if self.mm_opt_constrain_bonds:
            integrator = mm.VerletIntegrator(0.001)
            sim = mmapp.Simulation(top, mmsys_bak, integrator)
            sim.context.setPositions(state.getPositions())
            self.ostream.print_info(
                f"Minimizing {note} molecule without bond constraints.")
            self.ostream.flush()
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
