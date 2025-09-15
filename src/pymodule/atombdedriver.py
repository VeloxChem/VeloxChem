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
import time
import numpy as np
import sys

from .veloxchemlib import mpi_master,hartree_in_kcalpermol, hartree_in_kjpermol
from .outputstream import OutputStream
from .molecule import Molecule
from .atomtypeidentifier import AtomTypeIdentifier
from .errorhandler import assert_msg_critical
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .optimizationdriver import OptimizationDriver


class AtomBdeDriver:
    """
    This class implements a driver for calculating bond dissociation energies (BDEs) of hydrogen or other target atoms in molecules.

    Overview:
    ---------
    - Calculates the BDE for each unique hydrogen (or target) atom in a molecule.
    - Target atoms must be singly connected (i.e., bonded to only one other atom).
    - Use `self.compute()` to calculate BDEs for each unique hydrogen atom.
    - Use `self.show()` to visualize the results in kJ/mol (default: only hydrogens bonded to sp3 carbons).
    - Requires two sets of basis sets and exchange-correlation functionals: one for geometry optimization, one for single-point energy calculation.

    Workflow:
    ---------
    1. Optimize the whole molecule geometry using `basis_set1` and `functional1`.
    2. Perform a single-point energy calculation with `basis_set2` and `functional2`.
    3. Analyze the molecule to:
        a. Determine atom types, equivalent atom groups, and connectivity.
        b. Filter unique hydrogen/target atoms based on equivalence and bonding.
    4. For each unique hydrogen/target atom:
        a. Remove it from the optimized molecule to generate the corresponding radical.
        b. Optimize the radical geometry and perform a single-point energy calculation.
    5. Compute the single-point energy of the isolated hydrogen/target atom.
    6. Calculate the BDE:
        BDE = (radical/leftover molecule energy) + (hydrogen/target atom energy) - (whole molecule energy)

    Defaults:
    ---------
    - Atom(default hydrogen) radical: doublet, neutral, bonded to sp3 carbon.
    - Optimization: def2-svp/blyp/RI/grid_level2, SCF convergence 1e-3.
    - Single-point: def2-tzvp/b3lyp, SCF convergence 1e-3.
    - Level shifting of 0.5 applied to radical SCF drivers.
    - If SCF for a radical does not converge, BDE is set to 0.0 and displayed as "fail".

    Usage:
    ------
    - Use `_generate_radical_molecules()` to generate radical candidates from a molecule.
    - Use `compute(molecules)` to calculate BDEs for a single molecule or a list of molecules.
    - Use `show()` to visualize results.
    - Input to `compute()` should be a `Molecule` object (or a list thereof) with correct charge and multiplicity set.
    - Set the target atom type, multiplicity, and charge via:
        `self.target_atom = "H"`
        `self.target_atom_multiplicity = 2` #hydrogen radical
        `self.target_atom_charge = 0` #neutral hydrogen radical
    - Set radical/leftover molecule multiplicity and charge via:
        `self.mol_rad_multiplicity = 2` #radical molecule when removing H radical
        `self.mol_rad_charge = 0` #neutral radical molecule when removing H radical
    - Mute SCF and optimization output by setting `self.mute_output = True`.
    - Set energy units via `self.energy_unit = "kcal"`, `"kJ"`, or `"hartree"` or  `"au"`.
    - Customize basis sets and functionals via:
        `self.basis_sets = ["def2-svp", "def2-tzvp"]`
        `self.xcfunctionals = ["blyp", "b3lyp"]`
    """
    def __init__(self, comm=None, ostream=None):
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self._comm = comm
        self._rank = comm.Get_rank()
        self._size = comm.Get_size()
        self.ostream = ostream

        #two basis sets needed for optimization and final single point energy calculation
        self.basis_sets = ["def2-svp", "def2-tzvp"]  
        self.xcfunctionals = ["blyp", "b3lyp"]
        self.radical_level_shifting = 0.5
        #whole molecule optimization workflow
        self.mol_scf_drv = ScfRestrictedDriver()
        self.mol_scf_drv.conv_thresh = 1e-3
        self.mol_scf_drv.ri_coulomb = True
        self.mol_scf_drv.grid_level = 2
        self.mol_opt_drv = OptimizationDriver(self.mol_scf_drv)
        self.mol_opt_drv.conv_energy = 1e-04
        self.mol_opt_drv.conv_drms = 1e-02
        self.mol_opt_drv.conv_dmax = 2e-02
        self.mol_opt_drv.conv_grms = 4e-03
        self.mol_opt_drv.conv_gmax = 8e-03
        self.mol_final_scf_drv = ScfRestrictedDriver()
        self.mol_final_scf_drv.conv_thresh = 1e-3
        #radical optimization workflow
        self.radical_scf_drv = ScfUnrestrictedDriver()
        self.radical_scf_drv.conv_thresh = 1e-3
        self.radical_scf_drv.ri_coulomb = True
        self.radical_scf_drv.grid_level = 2
        self.radical_opt_drv = OptimizationDriver(self.radical_scf_drv)
        self.radical_opt_drv.conv_energy = 1e-04
        self.radical_opt_drv.conv_drms = 1e-02
        self.radical_opt_drv.conv_dmax = 2e-02
        self.radical_opt_drv.conv_grms = 4e-03
        self.radical_opt_drv.conv_gmax = 8e-03
        self.radical_final_scf_drv = ScfUnrestrictedDriver()
        self.radical_final_scf_drv.conv_thresh = 1e-3
        #Target atom radical scf driver
        self.target_atom_final_scf_drv = ScfUnrestrictedDriver()
        self.target_atom_final_scf_drv.conv_thresh = 1e-3
        #analyze all atoms in the molecule, setting for _atom_analyzer
        self.analyze_allatoms = False
        self.target_atom = "H"
        self.target_atom_multiplicity = 2
        self.target_atom_charge = 0
        self.mol_rad_multiplicity = 2
        self.mol_rad_charge = 0
        self.show_mol = False
        
        #save files
        self.save_files = True

        #attributes for atom analysis
        self.atom_idx = None
        self.labels = None
        self.connectivity_matrix = None
        self.atoms_types = None
        self.atom_info_dict = {}
        self.target_atom_info_dict = {}
        self.use_equiv = True
        self.only_sp3_carbon_connections = True
        self.only_hartree_fock = False
        self.mute_output = False
        self.energy_unit = "kJ"  #kcal, kJ, hartree

    def _update_drv_input_attributes(self, old_drv, new_drv):
        #if attribute is in input keywords, then update it but skip _scf_type of old_drv
        update_keywords = []
        if hasattr(new_drv, '_input_keywords'):
            for key in new_drv._input_keywords.keys():# scf driver
                update_keywords.extend(new_drv._input_keywords[key].keys())
        if hasattr(new_drv, 'input_keywords'): # optimization driver
            for key in new_drv.input_keywords.keys():
                update_keywords.extend(new_drv.input_keywords[key].keys())

        for attr in vars(old_drv):
            # skip attributes of _scf_type
            if attr == "_scf_type":
                continue
            if ((attr in update_keywords) & hasattr(old_drv, attr)):
                setattr(new_drv, attr, getattr(old_drv, attr))
        return new_drv


    def _check_proper_drv(self,whole_molecule):
        #reset scf drv for 'molecule', 'radical', 'target_atom' based on multiplicity
        #default target atom (hydrogen) is doublet with ScfUnrestrictedDriver
        if self.target_atom_multiplicity == 1: #singlet
            self.target_atom_final_scf_drv = self._update_drv_input_attributes(self.target_atom_final_scf_drv, ScfRestrictedDriver())
        elif self.target_atom_multiplicity > 1: #multiplet
            self.target_atom_final_scf_drv = self._update_drv_input_attributes(self.target_atom_final_scf_drv, ScfUnrestrictedDriver())

        #default radical is doublet with ScfUnrestrictedDriver
        if self.mol_rad_multiplicity == 1: #singlet
            self.radical_scf_drv = self._update_drv_input_attributes(self.radical_scf_drv, ScfRestrictedDriver())
            self.radical_final_scf_drv = self._update_drv_input_attributes(self.radical_final_scf_drv, ScfRestrictedDriver())
            self.radical_opt_drv = self._update_drv_input_attributes(self.radical_opt_drv, OptimizationDriver(self.radical_scf_drv))
        elif self.mol_rad_multiplicity > 1: #multiplet
            self.radical_scf_drv = self._update_drv_input_attributes(self.radical_scf_drv, ScfUnrestrictedDriver())
            self.radical_final_scf_drv = self._update_drv_input_attributes(self.radical_final_scf_drv, ScfUnrestrictedDriver())
            self.radical_opt_drv = self._update_drv_input_attributes(self.radical_opt_drv, OptimizationDriver(self.radical_scf_drv))

        #default molecule is closed shell
        if whole_molecule.get_multiplicity() > 1: #multiplet
            self.mol_scf_drv = self._update_drv_input_attributes(self.mol_scf_drv, ScfUnrestrictedDriver())
            self.mol_final_scf_drv = self._update_drv_input_attributes(self.mol_final_scf_drv, ScfUnrestrictedDriver())
            self.mol_opt_drv = self._update_drv_input_attributes(self.mol_opt_drv, OptimizationDriver(self.mol_scf_drv))
        elif whole_molecule.get_multiplicity() == 1: #singlet
            self.mol_scf_drv = self._update_drv_input_attributes(self.mol_scf_drv, ScfRestrictedDriver())
            self.mol_final_scf_drv = self._update_drv_input_attributes(self.mol_final_scf_drv, ScfRestrictedDriver())
            self.mol_opt_drv = self._update_drv_input_attributes(self.mol_opt_drv, OptimizationDriver(self.mol_scf_drv))
    
    def _check_scf_mute(self):
        if self.mute_output:
            self.mol_scf_drv.ostream.mute()
            self.mol_final_scf_drv.ostream.mute()
            self.radical_scf_drv.ostream.mute()
            self.radical_final_scf_drv.ostream.mute()
            self.target_atom_final_scf_drv.ostream.mute()
            self.mol_opt_drv.ostream.mute()
            self.radical_opt_drv.ostream.mute()
    

    def _check_radical_level_shifting(self):
        if self.radical_level_shifting is not None:
            if self.mol_rad_multiplicity <2:
                self.ostream.print_warning(
                    "level shifting is only applied to open shell calculations, "
                    "will skip setting level shifting for radical scf drivers"
                    )
                self.ostream.flush()
                return
            if hasattr(self.radical_scf_drv, 'level_shifting'):
                self.radical_scf_drv.level_shifting = self.radical_level_shifting
                self.radical_final_scf_drv.level_shifting = self.radical_level_shifting

    def _check_functionals(self):
        """
        check if exchange-correlation functionals as list are set, suppose to be [functional1,functional2]
        if not, then check if the user set the only_hartree_fock flag or extract the scf_driver functionals
        """
        #if only hartree fock, then set functionals to None
        if self.only_hartree_fock:
            self.xcfunctionals = [] #reset xcfunctionals to an empty list
            self.ostream.print_info(
                "Using HF for all calculations")
            self.ostream.flush()
            return

        #if not setting HF only, then check if functionals are provided in the scf drivers
        if len(self.xcfunctionals) == 0 and not self.only_hartree_fock:
            self.ostream.print_warning(
                "No functionals provided, checking scf drv functionals")
            self.ostream.flush()
            if self.mol_scf_drv.xcfun is not None:
                self.xcfunctionals.append(self.mol_scf_drv.xcfun)
                if self.mol_final_scf_drv.xcfun is not None:
                    self.xcfunctionals.append(self.mol_final_scf_drv.xcfun)
                else:
                    self.xcfunctionals.append(self.mol_scf_drv.xcfun)
            elif self.radical_scf_drv.xcfun is not None:
                self.xcfunctionals.append(self.radical_scf_drv.xcfun)
                if self.radical_final_scf_drv.xcfun is not None:
                    self.xcfunctionals.append(
                        self.radical_final_scf_drv.xcfun)
                else:
                    self.xcfunctionals.append(self.radical_scf_drv.xcfun)
            else:
                assert_msg_critical(
                    "No exchange-correlation functionals provided and no functional found in scf drivers, "
                    "please provide functionals or set only_hartree_fock=True"
                )
        elif len(self.xcfunctionals) == 1:
            self.xcfunctionals.append(self.xcfunctionals[0])

        self.ostream.print_info(
            f"Using provided exchange-correlation functionals{self.xcfunctionals}"
        )
        self.ostream.flush()
        #set functionals to drivers
        self.mol_scf_drv.xcfun = self.xcfunctionals[0]
        self.mol_final_scf_drv.xcfun = self.xcfunctionals[1]
        self.radical_scf_drv.xcfun = self.xcfunctionals[0]
        self.radical_final_scf_drv.xcfun = self.xcfunctionals[1]
        self.target_atom_final_scf_drv.xcfun = self.xcfunctionals[1]

    def _get_equiv(self, mol):
        """
        Get the equivalent atom types for a molecule.

        :param mol: 
            The molecule to analyze.

        :return: 
            tuple (a list of atom types, a list of equivalent atom types).
        """
        idtf = AtomTypeIdentifier()
        idtf.ostream.mute()
        atom_type = idtf.generate_gaff_atomtypes(mol)
        idtf.identify_equivalences()
        equivalent_charges = idtf.equivalent_charges
        self.atom_type = atom_type
        self.equivalent_charges = equivalent_charges
        return atom_type, equivalent_charges

    def _analyze_equiv(self, molecule):
        """
        analyze and process the equivalence of atom types in a molecule.

        :param molecule: The molecule object

        :return: 
            A tuple (a list of atom types, a list of equivalent atom types).
        """
        atom_type, equiv = self._get_equiv(molecule)
        equiv_atoms_groups = []
        # Split the string by comma
        substrings = equiv.split(",")
        # Split each substring by "="
        for substr in substrings:
            unit = substr.split("=")
            equiv_atoms_groups.append(unit)
        # map str to int
        equiv_atoms_groups = [list(map(int, x)) for x in equiv_atoms_groups]
        for i in range(len(equiv_atoms_groups)
                      ):  # this is to make the atom index start from 0
            equiv_atoms_groups[i] = [j - 1 for j in equiv_atoms_groups[i]]
        self.atom_type = atom_type
        self.equiv_atoms_groups = equiv_atoms_groups
        return atom_type, equiv_atoms_groups

    def _search_in_equiv_atoms_groups(self, equiv_atoms_groups, atom_idx):
        """
        Search for the equivalent atom group containing the specified atom index.

        :param equiv_atoms_groups: The list of equivalent atom groups.
        :param atom_idx: The atom index to search for.

        :return: The equivalent atom group containing the specified atom index, or an empty list if not found.
        """
        for group in equiv_atoms_groups:
            if atom_idx in group:
                return group
        return [atom_idx]

    def _atoms_analyzer(self, molecule):
        """
        Analyze and process the atom types in a molecule.
        This function is aiming to create a dictionary of atom information.
        including atom types, coordinates, equivalent atom types, and connectivity information.

        keys are "atom_type","equiv_group","H_connected_atom","coord"
        the equiv_group is used to identify equivalent atoms in the molecule and can be used for skipping redundant calculations.
        the H_connected_atom is used to identify if the atom is connected to a sp3 carbon atom.
        
        :param molecule: The molecule object

        :return: A dictionary containing atom information.
        if setting self.analyze_allatoms= True, then return the all atoms info dictionary 
        if setting self.analyze_allatoms= False, then return the target atom info dictionary

        """
        if self.show_mol:
            molecule.show(atom_indices=True)
        atom_type, equiv_atoms_groups = self._analyze_equiv(molecule)
        con_matrix = molecule.get_connectivity_matrix()
        labels = molecule.get_labels()
        coords = molecule.get_coordinates_in_angstrom()
        atom_info_dict = {}
        for i in range(len(labels)):
            atom_info_dict[labels[i] + "_" + str(i)] = {
                "atom_type": atom_type[i],
                "equiv_group": self._search_in_equiv_atoms_groups(
                    equiv_atoms_groups, i),
                "H_connected_atom": self._add_H_connected_atom_info(
                    labels, atom_type, i, con_matrix),
                "coord": coords[i]
            }
        if self.analyze_allatoms:
            self.atom_info_dict = atom_info_dict
            return atom_info_dict
        else:
            target_atom_info_dict = {}
            for key in atom_info_dict.keys():
                if self.target_atom == key.split("_")[0]:
                    target_atom_info_dict[key] = atom_info_dict[key]
            self.target_atom_info = target_atom_info_dict
            return target_atom_info_dict

    def _add_H_connected_atom_info(self, labels, atoms_types, atom_idx,
                                   connectivity_matrix):
        """
        this function is to use connectivity matrix to find the H connected atom information

        :param labels: The labels of the atoms in the molecule
        :param atoms_types: The types of the atoms in the molecule
        :param atom_idx: The index of the atom to analyze
        :param connectivity_matrix: The connectivity matrix of the molecule

        :return: A string containing the H connected atom information,this string is 
        in the format of "atomlabel_atomindex_atomtype".
        """
        # if the atom is not H or target_atom_type, then return ''
        if labels[atom_idx] != self.target_atom:
            return ""
        else:
            # search it in connectivity matrix to get the connected atom index
            con_info = connectivity_matrix[atom_idx]
            # if the atom is H, then it should have only one connected atom, which value should be 1
            connected_atom_idx = np.where(con_info == 1)[0]
            if len(connected_atom_idx) != 1:
                assert_msg_critical(
                    f"{self.target_atom} atom should have only one connected atom")
            connected_atom = labels[connected_atom_idx[0]]
            connected_atom_type = atoms_types[connected_atom_idx[0]]
            h_connected_atom_info = (str(connected_atom) + "_" +
                                     str(connected_atom_idx[0]) + "_" +
                                     str(connected_atom_type))
            return h_connected_atom_info

    def _fetch_unique_H(self,
                        target_atoms_dict,
                        use_equiv=False,
                        only_sp3_carbon_connections=False):
        """
        this function is used as filter for unique (hydrogen) atoms or only hydrogen connected to sp3 carbon

        :param target_atoms_dict: The dictionary containing target atom information.
        :param use_equiv: Whether to use equivalent hydrogen atoms.
        :param only_sp3_carbon_hydrogen: Whether to only include hydrogen connected to sp3 carbon.

        :return: A tuple containing the filtered hydrogen keys and their indices.

        """
        if not use_equiv:
            for key in target_atoms_dict.keys():
                target_atoms_dict[key]["equiv_group"] = []

        target_atoms_record = []
        unique_target_atom_keys = []
        for key in target_atoms_dict.keys():
            if key.split("_")[0] == self.target_atom and int(
                    key.split("_")[1]) not in target_atoms_record:
                target_atoms_record.extend(target_atoms_dict[key]["equiv_group"])
                unique_target_atom_keys.append(key)
        unique_target_atom_indices = [
            int(x.split("_")[1]) for x in unique_target_atom_keys
        ]
        if not only_sp3_carbon_connections:
            return unique_target_atom_keys, unique_target_atom_indices
        else:
            sp3_carbon_unique_target_atoms_indices = []
            sp3_carbon_unique_target_atoms_keys = []
            for key in unique_target_atom_keys:
                if target_atoms_dict[key]["H_connected_atom"] != "":
                    connected_atom = target_atoms_dict[key][
                        "H_connected_atom"].split("_")[0]
                    connected_atom_type = target_atoms_dict[key][
                        "H_connected_atom"].split("_")[2]
                    if connected_atom == "C" and connected_atom_type == "c3":
                        sp3_carbon_unique_target_atoms_indices.append(
                            int(key.split("_")[1]))
                        sp3_carbon_unique_target_atoms_keys.append(key)
            return sp3_carbon_unique_target_atoms_keys, sp3_carbon_unique_target_atoms_indices

    def _update_equiv_target_atoms_dissociation_energy(
            self, unique_target_atoms_dissociation_energies, unique_target_atoms_keys,
            target_atoms_dict):
        """
        this function is used to add the dissociation energy of equivalent target atoms in the target_atoms_dict

        :param unique_target_atoms_dissociation_energies: The bond dissociation energies of the unique target atoms.
        :param unique_target_atoms_keys: The keys of the unique target atoms.
        :param target_atoms_dict: The dictionary containing target atom information.

        :return: tuple. The updated target_atoms_dict, a list of tuple, 
        like [(bond dissociation energy, coordinates),(bond dissociation energy, coordinates)].
        """
        hartree2kcal = hartree_in_kcalpermol()
        hartree2kj = hartree_in_kjpermol()
        target_atoms_bdes_kj_coords = []
        unique_target_atoms_bdes_kj_coords = []
        for i in range(len(unique_target_atoms_keys)):
            key = unique_target_atoms_keys[i]
            energy_hartree = unique_target_atoms_dissociation_energies[i]
            if "kj" in self.energy_unit.lower():
                unique_target_atoms_bdes_kj_coords.append(
                    (energy_hartree * hartree2kj, target_atoms_dict[key]["coord"]))
            elif "kcal" in self.energy_unit.lower():
                unique_target_atoms_bdes_kj_coords.append(
                    (energy_hartree * hartree2kcal, target_atoms_dict[key]["coord"]))
            elif ("hartree" in self.energy_unit.lower() or "au" in self.energy_unit.lower()):
                unique_target_atoms_bdes_kj_coords.append(
                    (energy_hartree, target_atoms_dict[key]["coord"]))
            equiv_group = target_atoms_dict[key]["equiv_group"]
            for j in equiv_group:
                target_atoms_dict[
                    str(self.target_atom)+"_" + str(j)]["dissociation_energy_hartree"] = energy_hartree
                target_atoms_dict[str(self.target_atom)+"_" + str(j)]["dissociation_energy_kcal"] = (
                    energy_hartree * hartree2kcal)
                target_atoms_dict[str(self.target_atom)+"_" + str(j)]["dissociation_energy_kj"] = (
                    energy_hartree * hartree2kj)
                target_atoms_bdes_kj_coords.append(
                    (energy_hartree * hartree2kj,
                     target_atoms_dict[str(self.target_atom)+"_" + str(j)]["coord"]))
        return (target_atoms_dict, target_atoms_bdes_kj_coords,
                unique_target_atoms_bdes_kj_coords)

    def _print_target_atoms_bond_dissociation_energy(self,
                                                 target_atoms_dict,
                                                 unit="kcal"):
        """
        this function is used to print the bond dissociation energy of target atoms, with targeting energy units
        """
        self.ostream.print_info("-" * 50)
        self.ostream.print_info(f"bond dissociation energy of {self.target_atom} atoms:")
        self.ostream.flush()
        for key in target_atoms_dict.keys():
            if "dissociation_energy_hartree" in target_atoms_dict[key].keys():
                #skip the SCF not converged ones, whose dissociation_energy_hartree is set to 0.0
                if target_atoms_dict[key]["dissociation_energy_hartree"] - 0.0 < 1.0e-6:
                    continue
                if "kcal" in unit.lower():
                    self.ostream.print_info(
                        f"{key} {round(target_atoms_dict[key]['dissociation_energy_kcal'], 1)} kcal/mol"
                    )
                    self.ostream.flush()
                elif "kj" in unit.lower():
                    self.ostream.print_info(
                        f"{key} {round(target_atoms_dict[key]['dissociation_energy_kj'], 1)} kJ/mol"
                    )
                    self.ostream.flush()
                elif ("hartree" in unit.lower() or "au" in unit.lower()):
                    self.ostream.print_info(
                        f"{key} {target_atoms_dict[key]['dissociation_energy_hartree']} hartree"
                    )
                    self.ostream.flush()

    def _compute_whole_mol_scf_energy(self, molecule, mol_idx):
        """
        this function is used to compute the whole molecule SCF energy with optimization and single-point energy calculation
        close shell with ScfRestrictedDriver

        :param molecule: The molecule object.

        :return: The optimized molecule object.

        """
        self.ostream.print_info(
            "Optimizing geometry of the molecule before removing target atoms...")
        self.ostream.flush()
        if self.save_files:
            self.mol_scf_drv.filename = f'bde_mol_{mol_idx+1}'
            self.mol_opt_drv.filename = f'bde_mol_{mol_idx+1}_opt'
            self.mol_final_scf_drv.filename = f'bde_mol_{mol_idx+1}_final'
        else:
            self.mol_scf_drv.filename = None
            self.mol_opt_drv.filename = None
            self.mol_final_scf_drv.filename = None

        basis_set1 = MolecularBasis.read(molecule, self.basis_sets[0])
        basis_set2 = MolecularBasis.read(molecule, self.basis_sets[1])
        scf_results = self.mol_scf_drv.compute(molecule, basis_set1)
        opt_results = self.mol_opt_drv.compute(molecule, basis_set1, scf_results)

        opt_molecule = Molecule.read_xyz_string(opt_results["final_geometry"])
        opt_molecule.set_charge(molecule.get_charge())
        opt_molecule.set_multiplicity(molecule.get_multiplicity())

        # final energy
        final_single_point_scf_result = self.mol_final_scf_drv.compute(opt_molecule, basis_set2)
        if not self.mol_final_scf_drv.is_converged:
            self.new_method()
            self.ostream.flush()
        self.whole_mol_single_point_scf_energy = self.mol_final_scf_drv.get_scf_energy()
        return opt_molecule

    def new_method(self):
        self.ostream.print_warning(
                "SCF of final single point calculation did not converge for the whole molecule, "
                "the result is not reliable, please modify the settings if needed"
            )

    def _compute_target_atoms_radical_scf_energy(self, mol_idx):
        """
        this function is used to compute the target atoms radical SCF energy
        open shell with ScfUnrestrictedDriver
        """

        target_atom = Molecule.read_str(f"{self.target_atom} 0.0 0.0 0.0")
        target_atom.set_multiplicity(self.target_atom_multiplicity)
        target_atom.set_charge(self.target_atom_charge)
        basis_set2 = MolecularBasis.read(target_atom, self.basis_sets[1])
        if self.save_files:
            self.target_atom_final_scf_drv.filename = f'bde_mol_{mol_idx+1}_target_{self.target_atom}_final'
        else:
            self.target_atom_final_scf_drv.filename = None
        scf_resultsH = self.target_atom_final_scf_drv.compute(target_atom, basis_set2)
        if not self.target_atom_final_scf_drv.is_converged:
            self.ostream.print_warning(
                f"SCF of final single point calculation did not converge for {self.target_atom} radical or ion, "
                f"the result is not reliable, please modify the settings if needed"
            )
            self.ostream.flush()
        self.target_atom_single_point_scf_energy = self.target_atom_final_scf_drv.get_scf_energy()

    def _remove_atom_by_idx(self, mol, atom_indices_to_remove, carbon_indices):
        """
        this function is used to remove atoms from a molecule by using their indices 
        and find the exposed carbon for guess

        :param mol: The molecule object.
        :param atom_indices_to_remove: The list of atom indices to remove.
        :param carbon_indices: The list of carbon atom indices.

        :returns: A list of radical molecule objects and a list of exposed carbon indices.
        """

        mol_string = mol.get_xyz_string()
        number_of_atoms = mol.number_of_atoms()
        mol_stringlist = mol_string.split("\n")
        # Identify the lines that start with atom and save the positions
        allmolecules = []
        radical_carbon_indices = []
        for idx, carbon_idx in zip(atom_indices_to_remove, carbon_indices):
            new_mol = mol_stringlist.copy()
            # remove the index+2 line, the first line is the number of atoms, the second line is the comment
            new_mol.pop(idx + 2)
            # Update the number of atoms
            new_mol[0] = str(number_of_atoms - 1)
            new_mol = "\n".join(new_mol)
            allmolecules.append(Molecule.read_xyz_string(new_mol))
            if idx > carbon_idx:
                radical_carbon_indices.append(carbon_idx)
            else:
                radical_carbon_indices.append(carbon_idx - 1)
        return allmolecules, radical_carbon_indices

    def _compute_mol_rad_scf_energy(self, mol, radical_carbon_idx, run_idx, mol_idx):
        """
        this function is used to compute the SCF energy for the given radical molecule 
        optimization and single point energy calculation.
        open shell with ScfUnrestrictedDriver

        :param mol: The molecule object.
        :param radical_carbon_idx: The index of the radical carbon atom for guess.
        :param run_idx: The index of the current radical molecule among all radical molecules. Can be removed because it is only used to save file

        :return: The SCF energy of the radical molecule. or None if the SCF did not converge

        """
        self.ostream.print_info("-" * 50)
        if self.show_mol:
            mol.show(atom_indices=True)
            
        step_start = time.time()
        mol.set_multiplicity(self.mol_rad_multiplicity)
        mol.set_charge(self.mol_rad_charge)
        basis_set1 = MolecularBasis.read(mol, self.basis_sets[0])
        basis_set2 = MolecularBasis.read(mol, self.basis_sets[1])
        if self.save_files:
            self.radical_scf_drv.filename = f'bde_mol_{mol_idx+1}_{run_idx+1}'
            self.radical_opt_drv.filename = f'bde_mol_{mol_idx+1}_{run_idx+1}_opt'
            self.radical_final_scf_drv.filename = f'bde_mol_{mol_idx+1}_{run_idx+1}_final'
        else:
            self.radical_scf_drv.filename = None
            self.radical_opt_drv.filename = None
            self.radical_final_scf_drv.filename = None

        #set guess for radical optimization
        if self.mol_rad_multiplicity != 1:
            self.radical_scf_drv.guess_unpaired_electrons = f'{radical_carbon_idx+1}({self.mol_rad_multiplicity-1}.0)'
        
        try:
            scf_resultsmol = self.radical_scf_drv.compute(mol, basis_set1)
            opt_results_rad = self.radical_opt_drv.compute(mol, basis_set1, scf_resultsmol)
        except Exception as e:
            if "ScfGradientDriver: SCF did not converge" in str(e):
                self.ostream.print_warning(
                    f"SCF did not converge for radical molecule {run_idx+1},will skip this radical, please modify the settings if needed"
                )
                self.ostream.flush()
                return None
        mol = Molecule.read_xyz_string(opt_results_rad["final_geometry"])
        mol.set_multiplicity(self.mol_rad_multiplicity)

        if self.mol_rad_multiplicity != 1:
            self.radical_final_scf_drv.guess_unpaired_electrons = f'{radical_carbon_idx+1}({self.mol_rad_multiplicity-1}.0)'
        mol.set_charge(self.mol_rad_charge)
        scf_results_rad_big = self.radical_final_scf_drv.compute(mol, basis_set2)
 
        if not self.radical_final_scf_drv.is_converged:
            self.ostream.print_warning(
                f"SCF of final single point calculation did not converge for radical molecule {run_idx+1}, "
                f"will skip this radical or molecule, please modify the settings if needed"
            )
            self.ostream.flush()
            return None

        radical_single_point_scf_energy = self.radical_final_scf_drv.get_scf_energy()
        step_end = time.time()
        self.ostream.print_info(
            f"SCF energy for radical molecule {run_idx+1} is {radical_single_point_scf_energy} hartree"
        )
        self.ostream.print_info(
            f"Time cost : {round(step_end - step_start, 2)} seconds"
        )
        self.ostream.flush()
        return radical_single_point_scf_energy

    def _show_bde_on_atom(self,
                          molecule,
                          width=400,
                          height=300,
                          atom_indices=False,
                          atom_labels=False,
                          one_indexed=True,
                          bdes_coords=None):
        """
        Creates a 3D view with py3dmol.

        :param width:
            The width.
        :param height:
            The height.
        :param atom_indices:
            The flag for showing atom indices (1-based).
        :param atom_labels:
            The flag for showing atom labels.
        :param one_indexed:
            The flag for using one-based indexing. Supposed to be True.
        :param bdes_coords:
            The list of tuples of (bde, coordinates) of target atoms.
        """

        try:
            import py3Dmol
            viewer = py3Dmol.view(width=width, height=height)
            viewer.addModel(molecule.get_xyz_string())
            viewer.setViewStyle({"style": "outline", "width": 0.05})
            viewer.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
            if atom_indices or atom_labels:
                coords = molecule.get_coordinates_in_angstrom()
                labels = molecule.get_labels()
                for i in range(coords.shape[0]):
                    text = ''
                    if atom_labels:
                        text += f'{labels[i]}'
                    if atom_indices:
                        if one_indexed:
                            text += f'{i + 1}'
                        else:
                            text += f'{i}'
                    viewer.addLabel(
                        text, {
                            'position': {
                                'x': coords[i, 0],
                                'y': coords[i, 1],
                                'z': coords[i, 2],
                            },
                            'alignment': 'center',
                            'fontColor': 0x000000,
                            'backgroundColor': 0xffffff,
                            'backgroundOpacity': 0.0,
                        })
            #add bde based on coords and unique_bde_hartree
            for i in range(len(bdes_coords)):
                #bde coords is a list of tuple [(bde, (x, y, z)),(bde, (x, y, z))]
                #default in kJ/mol but changed by self.energy_unit
                #skip the scf not converged ones, whose bde is set to 0.0
                if bdes_coords[i][0] - 0.0 < 1.0e-6:
                    bde_kj = 'Fail'
                else:
                    bde_kj = round(bdes_coords[i][0],
                                0)  #only show the integer part
                    bde_kj = int(bde_kj)
                viewer.addLabel(
                    f'{bde_kj}', {
                        'position': {
                            'x': bdes_coords[i][1][0],
                            'y': bdes_coords[i][1][1],
                            'z': bdes_coords[i][1][2],
                        },
                        'alignment': 'top',
                        'fontColor': 'red',
                        'backgroundColor': 0xffffff,
                        'backgroundOpacity': 0.5,
                    })
            viewer.zoomTo()
            viewer.show()

        except ImportError:
            raise ImportError('Unable to import py3Dmol')

    def _generate_radical_molecules(self, whole_molecule: object):
        """
        only applied to single molecule
        Generate radical molecules from the whole molecule by removing hydrogen atoms or target atoms.

        :param whole_molecule: The whole molecule object to generate radicals from.
        :return: A tuple, containing the radical molecules, their corresponding carbon indices, and guess messages.

        """
        self._check_proper_drv(whole_molecule)
        self._check_scf_mute()
        self._check_functionals()
        self._check_radical_level_shifting()
        target_atoms_dict = self._atoms_analyzer(whole_molecule)
        unique_target_atoms_keys, unique_target_atoms_indices = self._fetch_unique_H(
            target_atoms_dict,
            use_equiv=self.use_equiv,
            only_sp3_carbon_connections=self.only_sp3_carbon_connections)
        #self.ostream.print_info(f'unique_target_atoms_indices {unique_target_atoms_indices}')
        carbon_indices = []
        # Generate radical molecules by removing target atoms
        # check connectivity matrix to find the connected carbon index for each target atom
        # if connected to more than one atom, then throw error
        conn_mat = whole_molecule.get_connectivity_matrix()
        for x in unique_target_atoms_indices:
            assert list(conn_mat[x]).count(1) == 1
            carbon_indices.append(list(conn_mat[x]).index(1))
        #self.ostream.print_info(f'connected carbon_indices {carbon_indices}')

        opt_whole_molecule = self._compute_whole_mol_scf_energy(whole_molecule)
        radical_molecules, radical_carbon_indices = self._remove_atom_by_idx(
            opt_whole_molecule, unique_target_atoms_indices, carbon_indices)
        guess_msg = []
        for i in range(len(radical_molecules)):
            guess = f"guess_unpaired_electrons = {radical_carbon_indices[i]+1}(1.0)"
            guess_msg.append(guess)
        return radical_molecules, radical_carbon_indices, guess_msg

    def _compute_single_molecule(self, whole_molecule: object, mol_idx: int):
        """
        used to compute hydrogen/target atom bond dissociation energies (BDEs) for a single molecule.
        """
        self._check_proper_drv(whole_molecule)
        self._check_scf_mute()
        self._check_functionals()
        self._check_radical_level_shifting()

        target_atoms_dict = self._atoms_analyzer(whole_molecule)
        unique_target_atoms_keys, unique_target_atoms_indices = self._fetch_unique_H(
            target_atoms_dict,
            use_equiv=self.use_equiv,
            only_sp3_carbon_connections=self.only_sp3_carbon_connections)

        #self.ostream.print_info(f'unique_target_atoms_indices {unique_target_atoms_indices}')
        carbon_indices = []
        conn_mat = whole_molecule.get_connectivity_matrix()
        for x in unique_target_atoms_indices:
            assert list(conn_mat[x]).count(1) == 1
            carbon_indices.append(list(conn_mat[x]).index(1))
        #self.ostream.print_info(f'connected carbon_indices {carbon_indices}')

        opt_whole_molecule = self._compute_whole_mol_scf_energy(whole_molecule,mol_idx)
        self._compute_target_atoms_radical_scf_energy(mol_idx)
        self.opt_whole_molecule = opt_whole_molecule
        molecules_rads, radical_carbon_indices = self._remove_atom_by_idx(
            opt_whole_molecule, unique_target_atoms_indices, carbon_indices)

        #calculate BDEs for each radical molecule
        unique_BDEs_hartree = []
        count = 1
        for run_idx, (mol_rad, radical_carbon_idx) in enumerate(
                zip(molecules_rads, radical_carbon_indices)):
            self.ostream.print_info("-" * 50)
            self.ostream.print_info(
                f"Computing energy of structure : {count} of {len(molecules_rads)}"
            )
            self.ostream.flush()
            mol_rad_scf_energy = None #reset to None for each radical molecule
            mol_rad_scf_energy = self._compute_mol_rad_scf_energy(
                mol_rad, radical_carbon_idx, run_idx, mol_idx)
            if mol_rad_scf_energy is None:
                bde_hartree = 0.0
                self.ostream.print_warning(
                    f"SCF did not converge for radical molecule {count}, setting BDE to 0.0, please modify the settings if needed"
                )
                self.ostream.flush()
            else:
                bde_hartree = mol_rad_scf_energy - self.whole_mol_single_point_scf_energy + self.target_atom_single_point_scf_energy
            unique_BDEs_hartree.append(bde_hartree)
            count += 1
        self.unique_BDEs_hartree = unique_BDEs_hartree
        self.unique_target_atoms_keys = unique_target_atoms_keys
        # loop the unique_target_atoms_indices to remove the H atoms from the molecule and
        # calculate the dissociation energy but save the energy for all equivalent H atoms
        # print the dissociation energy for each H atom
        self.target_atoms_dict, self.bdes_coords, self.unique_target_atoms_bdes_coords = self._update_equiv_target_atoms_dissociation_energy(
            unique_BDEs_hartree, unique_target_atoms_keys, target_atoms_dict)
        self._print_target_atoms_bond_dissociation_energy(self.target_atoms_dict,
                                                      unit=self.energy_unit)

    def compute(self, mol_list: list):
        """
        for computing bond dissociation energies (BDEs) for a list of molecules.
        can accept one molecule object as well.
        will save self.mols_bdes_list
        """
        self.mols_bdes_list = []
        if not isinstance(mol_list, list):
            mol_list = [mol_list]
        for mol_idx, mol in enumerate(mol_list):
            self._compute_single_molecule(mol,mol_idx)
            #use a dictionary to store all bde and target atom information
            mol_bde_dict = {
                'target_atoms_dict': self.target_atoms_dict,
                'bdes_coords': self.bdes_coords,
                'unique_target_atoms_keys': self.unique_target_atoms_keys,
                'unique_BDEs_hartree': self.unique_BDEs_hartree,
                'unique_target_atoms_bdes_coords': self.unique_target_atoms_bdes_coords,
                'opt_whole_molecule': self.opt_whole_molecule,
            }
            self.mols_bdes_list.append(mol_bde_dict)

    def show(self,
             atom_indices=False,
             atom_labels=False,
             width=400,
             height=300,
             unique_target_atoms=True):
        """
        this function is to visualize the bond dissociation energies (BDEs) on the target atoms in the molecule.
        allow atom indices and labels to be shown.
        """
        for item in self.mols_bdes_list:
            mol = item['opt_whole_molecule']

            if unique_target_atoms:
                bdes_coords = item['unique_target_atoms_bdes_coords']
            else:
                bdes_coords = item['bdes_coords']

            self._show_bde_on_atom(mol,
                                   width=width,
                                   height=height,
                                   atom_indices=atom_indices,
                                   atom_labels=atom_labels,
                                   one_indexed=True,
                                   bdes_coords=bdes_coords)
