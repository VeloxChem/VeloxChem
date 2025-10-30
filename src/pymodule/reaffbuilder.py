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
from .scfhessiandriver import ScfHessianDriver
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
from .openmmdynamics import OpenMMDynamics

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


class ReactionForceFieldBuilder():

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

        self.optimize_mol: bool = False
        self.reparameterize_bonds: bool = False
        self.reparameterize_angles: bool = False
        self.optimize_ff: bool = True
        self.optimize_steps: int = 1000
        self.optimize_conformer_snapshots: int = 10
        self.optimize_temp: int = 600
        self.optimize_dist_restraint_offset = 0.5  # Angstrom
        self.optimize_dist_restraint_k = 250000.0  # kJ mol^-1 nm^-2
        self.mm_opt_constrain_bonds: bool = True
        self.water_model: str = 'cspce'
        self.product_mapping: dict[int, int] | None = None  # one-indexed
        self.mute_scf: bool = True
        self.skip_reaction_matching: bool = False
        #Todo get a better functional and basis set from here https://pubs.acs.org/doi/10.1021/acs.jctc.3c00558
        self.hessian_xc_fun: str = 'B3LYP'
        #Todo get better basis set once we have f-functionals
        # Can (should?) be scaled up to def2-TZVPPD, and if only we had our ECP's by now
        self.hessian_basis = 'def2-SV_P_'

    def build_forcefields(
        self,
        reactant: Molecule | list[Molecule],
        product: Molecule | list[Molecule],
        reactant_partial_charges: list[float] | list[list[float]] | None = None,
        product_partial_charges: list[float] | list[list[float]] | None = None,
        reactant_hessians: np.ndarray | list[np.ndarray | None] | None = None,
        product_hessians: np.ndarray | list[np.ndarray | None] | None = None,
        reactant_total_multiplicity: int = -1,
        product_total_multiplicity: int = -1,
        forced_breaking_bonds: set[tuple[int, int]] | tuple = (),
        forced_forming_bonds: set[tuple[int, int]] | tuple = (),
        product_mapping: dict[int, int] | None = None,
    ):

        reactant_ff, reactant_ffs, reactant_total_charge = self._create_combined_forcefield(
            reactant,
            reactant_partial_charges,
            reactant_hessians,
            "REA",
            reactant_total_multiplicity,
        )

        product_ff, product_ffs, product_total_charge = self._create_combined_forcefield(
            product,
            product_partial_charges,
            product_hessians,
            "PRO",
            product_total_multiplicity,
        )

        assert reactant_total_charge == product_total_charge, f"Total charge of reactants {reactant_total_charge} and products {product_total_charge} must match"

        if not self.skip_reaction_matching and product_mapping is None:
            breaking_bonds_insert = "no forced breaking bonds"
            if len(forced_breaking_bonds) > 0:
                breaking_bonds_insert = f"forced breaking bonds: {forced_breaking_bonds}"

            forming_bonds_insert = "no forced forming bonds"
            if len(forced_forming_bonds) > 0:
                forming_bonds_insert = f"forced forming bonds: {forced_forming_bonds}"

            msg = "Matching reactant and product force fields with "
            msg += breaking_bonds_insert + " and " + forming_bonds_insert + "."
            self.ostream.print_info(msg)
            self.ostream.flush()

            # adjust for 1-indexed input of breaking bonds
            forced_breaking_bonds = {(bond[0] - 1, bond[1] - 1)
                                     for bond in forced_breaking_bonds}
            forced_forming_bonds = {(bond[0] - 1, bond[1] - 1)
                                    for bond in forced_forming_bonds}

            product_ff, product_mapping = self._match_reactant_and_product(
                reactant_ff,
                reactant_ff.molecule.get_element_ids(),
                product_ff,
                product_ff.molecule.get_element_ids(),
                forced_breaking_bonds,
                forced_forming_bonds,
            )
        elif product_mapping is not None:
            self.ostream.print_info(
                f"Skipping reaction matching because the mapping {product_mapping} is already provided"
            )
            product_mapping = {k - 1: v - 1 for k, v in product_mapping.items()}
        else:
            self.ostream.print_info("Skipping reaction matching")
        self.ostream.flush()

        if product_mapping is not None:
            product_ff = ReactionForceFieldBuilder._apply_mapping_to_forcefield(
                product_ff,
                product_mapping,
            )

            product_ff.molecule = ReactionForceFieldBuilder._apply_mapping_to_molecule(
                product_ff.molecule,
                product_mapping,
            )

        forming_bonds, breaking_bonds = self._summarise_reaction(
            reactant_ff, product_ff, self.ostream)

        for bond in breaking_bonds:
            reactant_ff.bonds[bond]['comment'] += ', broken in reaction'
        for bond in forming_bonds:
            product_ff.bonds[bond]['comment'] += ', formed in reaction'

        self.ostream.flush()
        if self.optimize_ff:
            # TODO this optimisation can likely be taken care of by the openmmdynamics class
            reactant_ff.molecule = self._optimize_molecule(
                reactant_ff.molecule.get_element_ids(),
                reactant_ff,
                forming_bonds,
                note='reactant',
            )
            product_ff.molecule = self._optimize_molecule(
                product_ff.molecule.get_element_ids(),
                product_ff,
                breaking_bonds,
                note='product',
            )
            # if len(reactant_ffs) > 1:
            # if len(product_ffs) > 1:

        return reactant_ff, product_ff, forming_bonds, breaking_bonds, reactant_ffs, product_ffs, product_mapping

    def _create_combined_forcefield(
        self,
        molecules: list[Molecule],
        partial_charges: list[list[float]]
        | list[float] | None,
        hessians: np.ndarray | list[np.ndarray | None] | None,
        name: str,
        total_multiplicity: int,
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
            total_multiplicity,
            name,
        )
        self.ostream.print_info(
            f"Combined reactant with total charge {combined_mol.get_charge()} and multiplicity {combined_mol.get_multiplicity()}"
        )
        self.ostream.flush()
        combined_ff = self._combine_forcefield(single_ffs)
        combined_ff.molecule = combined_mol

        return combined_ff, single_ffs, total_charge

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

        #Reparameterize the forcefield if necessary and requested
        unknown_pairs = set()
        unknown_params = set()
        if self.reparameterize_bonds:
            for key, bond in forcefield.bonds.items():
                if bond['comment'] == 'Guessed':
                    sorted_tuple = tuple(sorted(key))
                    unknown_pairs.add(sorted_tuple)
                    unknown_params.add(key)
        if self.reparameterize_angles:
            for key, angle in forcefield.angles.items():
                if angle['comment'] == 'Guessed':
                    sorted_tuple = tuple(sorted((key[0], key[1])))
                    unknown_pairs.add(sorted_tuple)
                    sorted_tuple = tuple(sorted((key[1], key[2])))
                    unknown_pairs.add(sorted_tuple)
                    unknown_params.add(key)

        if len(unknown_pairs) > 0:
            self.ostream.print_info("Reparameterising force field.")
            self.ostream.flush()
            if hessian is None:
                self.ostream.print_info(
                    f"Calculating hessian submatrices for atom pairs {unknown_pairs} to reparameterise the force field."
                )
                if molecule.get_multiplicity() == 1:
                    scf_drv = ScfRestrictedDriver()
                else:
                    scf_drv = ScfUnrestrictedDriver()
                self.ostream.flush()
                if self.mute_scf:
                    scf_drv.ostream.mute()
                basis = MolecularBasis.read(molecule, self.hessian_basis)
                scf_drv.xcfun = self.hessian_xc_fun
                scf_drv.dispersion = True
                scf_drv.compute(molecule, basis)
                if scf_drv.is_converged is False:
                    self.ostream.print_warning(
                        "SCF did not converge, increasing convergence threshold to 1.0e-4 and maximum itterations to 200."
                    )
                    self.ostream.flush()
                    scf_drv.conv_thresh = 1.0e-4
                    scf_drv.max_iter = 200
                    scf_drv.compute(molecule, basis)
                assert scf_drv.is_converged, f"SCF calculation for Hessian did not converge, aborting"

                hess_drv = ScfHessianDriver(scf_drv)
                if self.mute_scf:
                    hess_drv.ostream.mute()
                hess_drv.atom_pairs = list(unknown_pairs)
                hess_drv.compute(molecule, basis)
                hessian = np.copy(hess_drv.hessian)  # type: ignore
            else:
                cond = np.shape(hessian) == (molecule.number_of_atoms() * 3,
                                             molecule.number_of_atoms() * 3)
                msg = f"Hessian shape {np.shape(hessian)} should be square with width 3 times the number"
                msg += f" of atoms 3*{molecule.number_of_atoms()}={3*molecule.number_of_atoms()} in the molecule."
                assert_msg_critical(cond, msg)
            self.ostream.flush()
            forcefield.reparameterize(hessian,
                                      reparameterize_keys=unknown_params)
        return forcefield

    # Guesses how to combine molecular structures without overlapping them
    def _combine_molecule(self, molecules, total_multiplicity, name='MOL'):

        combined_molecule = Molecule()
        # pos = []
        charge = 0
        Sm1 = 0
        for i, mol in enumerate(molecules):
            charge += mol.get_charge()
            Sm1 += mol.get_multiplicity() - 1
            if combined_molecule.number_of_atoms() > 0:
                coords = combined_molecule.get_coordinates_in_angstrom()
                max_x = max(coords[:, 0])
                min_x = min(coords[:, 0])
                shift = max_x - min_x + 2
                self.ostream.print_info(
                    f"max_x: {max_x}, min_x: {min_x}, shifting next molecule by {shift} angstrom"
                )
            else:
                shift = 0
            coords = mol.get_coordinates_in_angstrom()
            min_x = min(coords[:, 0])
            coords[:, 0] -= min_x
            for elem, coord in zip(mol.get_element_ids(), coords):
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

    #Match the indices of the reactant and product forcefield generators
    def _match_reactant_and_product(
        self,
        reactant_ff: MMForceFieldGenerator,
        rea_elems: list,
        product_ff: MMForceFieldGenerator,
        pro_elems: list,
        breaking_bonds: set[tuple[int, int]],
        forming_bonds: set[tuple[int, int]],
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
            forming_bonds,
        )  # type: ignore
        if total_mapping is None:
            raise ValueError(
                "Could not find a mapping between the reactant and product force fields."
            )
        total_mapping = {v: k for k, v in total_mapping.items()}
        print_mapping = {k + 1: v + 1 for k, v in total_mapping.items()}
        self.ostream.print_info(f"Mapping: {print_mapping}")
        self.ostream.flush()
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
        forcefield.atom_info_dict = {}
        for l, ff in enumerate(forcefields):
            # Shift all atom keys by the current atom count so that every atom has a unique ID
            shift = atom_count
            mapping = {atom_key: atom_key + shift for atom_key in ff.atoms}
            ReactionForceFieldBuilder._apply_mapping_to_forcefield(ff, mapping)
            atom_count += len(ff.atoms)
            for atom in ff.atoms.values():
                atom['name'] = f"{atom['name']}{l+1}"
            forcefield.atom_info_dict.update(ff.atom_info_dict)
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
        new_ff_atoms = {}
        for atom_key in forcefield.atoms:
            key = mapping[atom_key]
            val = forcefield.atoms[atom_key]
            new_ff_atoms.update({key: val})

        new_atom_info = {}
        for atom_info_key in forcefield.atom_info_dict.keys():
            key = mapping[atom_info_key - 1] + 1
            val = forcefield.atom_info_dict[atom_info_key]
            val['ConnectedAtomsNumbers'] = [
                mapping[key - 1] + 1 for key in val['ConnectedAtomsNumbers']
            ]
            val['AtomNumber'] = mapping[val['AtomNumber'] - 1] + 1
            new_atom_info.update({key: val})

        # Sort the atoms by index
        forcefield.atoms = dict(sorted(new_ff_atoms.items()))
        forcefield.atom_info_dict = dict(sorted(new_atom_info.items()))

        forcefield.bonds = ReactionForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.bonds, mapping)
        forcefield.angles = ReactionForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.angles, mapping)
        forcefield.dihedrals = ReactionForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.dihedrals, mapping)
        forcefield.impropers = ReactionForceFieldBuilder._apply_mapping_to_parameters(
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

    def _summarise_reaction(self, reactant, product, ostream):
        """
        Summarises the reaction by printing the bonds that are being broken and formed.

        Returns:
            None
        """
        reactant_bonds = set(reactant.bonds)
        product_bonds = set(product.bonds)
        formed_bonds = product_bonds - reactant_bonds
        broken_bonds = reactant_bonds - product_bonds
        ostream.print_header("Reaction summary")
        ostream.print_header(f"{len(broken_bonds)} breaking bonds:")

        if len(broken_bonds) > 0:
            ostream.print_header(f"ReaType  ProType  ID - ReaType  ProType  ID")
        for bond_key in broken_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0] + 1
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1] + 1
            ostream.print_header(
                f"{reactant_type0:^9}{product_type0:^9}{id0:^2} - {reactant_type1:^9}{product_type1:^9}{id1:^2}"
            )
        ostream.print_blank()
        ostream.print_header(f"{len(formed_bonds)} forming bonds:")
        if len(formed_bonds) > 0:
            ostream.print_header("ReaType  ProType  ID - ReaType  ProType  ID")
        for bond_key in formed_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0] + 1
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1] + 1
            ostream.print_header(
                f"{reactant_type0:^9}{product_type0:^9}{id0:^2} - {reactant_type1:^9}{product_type1:^9}{id1:^2}"
            )
        ostream.print_blank()
        ostream.flush()
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

        opm_dyn = OpenMMDynamics()
        opm_dyn.ostream.mute()
        opm_dyn.openmm_platform = "CPU"

        opm_dyn.pdb = pdb
        opm_dyn.system = mmsys

        self.ostream.print_info(
            f"Running conformational sampling with {self.optimize_steps*self.optimize_conformer_snapshots} steps and {self.optimize_conformer_snapshots} snapshots at {self.optimize_temp} K for {note} molecule."
        )
        self.ostream.flush()
        try:
            conformers_dict = opm_dyn.conformational_sampling(
                ensemble='NVT',
                nsteps=self.optimize_steps * self.optimize_conformer_snapshots,
                snapshots=self.optimize_conformer_snapshots,
                temperature=self.optimize_temp,
            )

            min_arg = np.argmin(conformers_dict['energies'])
            new_molecule = conformers_dict['molecules'][min_arg]
            self.ostream.print_info(
                f"Found {len(conformers_dict['molecules'])} conformers with energies {conformers_dict['energies']} during optimization of the {note} molecule."
            )
            self.ostream.flush()
        except Exception as e:
            self.ostream.print_warning(
                f"OpenMM optimization of the {note} molecule failed with error: {e}. Reverting to original geometry."
            )
            self.ostream.flush()
            new_molecule = forcefield.molecule
            
        new_molecule.set_charge(forcefield.molecule.get_charge())
        new_molecule.set_multiplicity(forcefield.molecule.get_multiplicity())
        os.unlink(f'{name}.xml')
        os.unlink(f'{name}.pdb')
        os.unlink(f'{name}_sys.xml')
        return new_molecule

    def _add_reaction_bonds(self, forcefield, mmsys, changing_bonds, note):
        self.ostream.print_info(
            "Guessing intra-molecular constraints based on breaking bonds. This option can be turned off with mm_opt_constrain_bonds."
        )

        dist_restraint_expr = "0.5 * k * (r - r0)^2*step(r-r0)"
        dist_restraint_force = mm.CustomBondForce(dist_restraint_expr)
        dist_restraint_force.addPerBondParameter("r0")
        dist_restraint_force.addGlobalParameter("k",
                                                self.optimize_dist_restraint_k)

        for bond in changing_bonds:
            s1 = forcefield.atoms[bond[0]]['sigma']
            s2 = forcefield.atoms[bond[1]]['sigma']
            s = (s1 + s2) / 2

            r0 = s * (
                2**(1 / 6)
            ) + self.optimize_dist_restraint_offset * 0.1  # rmin = sigma * 2^(1/6)
            dist_restraint_force.addBond(bond[0], bond[1], [r0])
            self.ostream.print_info(
                f"Adding distance restraint for atoms {tuple(b+1 for b in bond)} with r0 {r0:.3f} and k {self.optimize_dist_restraint_k} for MM equilibration of the {note}"
            )
            self.ostream.flush()

        mmsys.addForce(dist_restraint_force)
        return mmsys
