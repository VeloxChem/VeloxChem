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

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .molecule import Molecule
from .atomtypeidentifier import AtomTypeIdentifier
from .errorhandler import assert_msg_critical
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .optimizationdriver import OptimizationDriver


class HydrogenBdeDriver:

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

        #set basis
        self.basis_sets = [
        ]  #two basis sets needed for optimization and final single point energy calculation
        self.xcfunctionals = []
        #whole molecule optimization workflow
        self.mol_scf_drv = ScfRestrictedDriver()
        self.mol_opt_drv = OptimizationDriver(self.mol_scf_drv)
        self.mol_final_scf_drv = ScfRestrictedDriver()
        #radical optimization workflow
        self.radical_scf_drv = ScfUnrestrictedDriver()
        self.radical_opt_drv = OptimizationDriver(self.radical_scf_drv)
        self.radical_final_scf_drv = ScfUnrestrictedDriver()
        #hydrogen radical scf driver
        self.hydrogen_final_scf_drv = ScfUnrestrictedDriver()
        #analyze all atoms in the molecule, setting for _atom_analyzer
        self.analyze_allatoms = False
        self.target_atom = "H"
        self.show_mol = False
        #attributes for atom analysis
        self.atom_idx = None
        self.labels = None
        self.connectivity_matrix = None
        self.atoms_types = None
        self.atom_info_dict = {}
        self.target_atom_info_dict = {}
        self.use_equiv = True
        self.only_sp3_carbon_hydrogen = True
        self.only_hartree_fock = False
        self.mute_output = False
        self.energy_unit = "kj"  #kcal, kj, au

    def check_scf_mute(self):
        if self.mute_output:
            self.mol_scf_drv.ostream.mute()
            self.mol_final_scf_drv.ostream.mute()
            self.radical_scf_drv.ostream.mute()
            self.radical_final_scf_drv.ostream.mute()
            self.hydrogen_final_scf_drv.ostream.mute()
            self.mol_opt_drv.ostream.mute()
            self.radical_opt_drv.ostream.mute()

    def check_functionals(self):
        """
        check if exchange-correlation functionals as list are set, suppose to be [functional1,functional2]
        if not, then check if the user set the only_hartree_fock flag or extract the scf_driver functionals
        """
        #if no functionals then all HF
        if len(self.xcfunctionals) == 0:
            if self.only_hartree_fock:
                #warning
                self.ostream.print_info("No functionals provided, using HF for all calculations")
                self.ostream.flush()
                return
            else:
                self.ostream.print_warning("No functionals provided, checking scf drv functionals")
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
                        "No exchange-correlation functionals provided and no functional found in scf drivers, please provide functionals or set only_hartree_fock=True"
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
        self.hydrogen_final_scf_drv.xcfun = self.xcfunctionals[1]

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
        return []

    def _atoms_analyzer(self, molecule):
        """
        Analyze and process the atom types in a molecule.
        This function is aiming to create a dictionary of atom information.
        including atom types, coordinates, equivalent atom types, and connectivity information.

        keys are "atom_type","equiv_group","H_connected_atom","coord"
        the equiv_group is used to identify equivalent atoms in the molecule and can be used for skip redundant calculations.
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

        :return: A string containing the H connected atom information,this string is in the format of "atomlabel_atomindex_atomtype".
        """
        # if the atom is not H, then return ''
        if labels[atom_idx] != "H":
            return ""
        else:
            # search it in connectivity matrix to get the connected atom index
            con_info = connectivity_matrix[atom_idx]
            # if the atom is H, then it should have only one connected atom, which value should be 1
            connected_atom_idx = np.where(con_info == 1)[0]
            if len(connected_atom_idx) != 1:
                assert_msg_critical(
                    "H atom should have only one connected atom")
            connected_atom = labels[connected_atom_idx[0]]
            connected_atom_type = atoms_types[connected_atom_idx[0]]
            h_connected_atom_info = (str(connected_atom) + "_" +
                                     str(connected_atom_idx[0]) + "_" +
                                     str(connected_atom_type))
            return h_connected_atom_info

    def _fetch_unique_H(self,
                        hydrogen_atoms_dict,
                        use_equiv=False,
                        only_sp3_carbon_hydrogen=False):
        """
        this function is used as filter for unique hydrogen atoms or only hydrogen connected to sp3 carbon

        :param hydrogen_atoms_dict: The dictionary containing hydrogen atom information.
        :param use_equiv: Whether to use equivalent hydrogen atoms.
        :param only_sp3_carbon_hydrogen: Whether to only include hydrogen connected to sp3 carbon.

        :return: A tuple containing the filtered hydrogen keys and their indices.

        """
        if not use_equiv:
            for key in hydrogen_atoms_dict.keys():
                hydrogen_atoms_dict[key]["equiv_group"] = []

        hydrogen_record = []
        unique_hydrogen_keys = []
        for key in hydrogen_atoms_dict.keys():
            if key.split("_")[0] == "H" and int(
                    key.split("_")[1]) not in hydrogen_record:
                hydrogen_record.extend(hydrogen_atoms_dict[key]["equiv_group"])
                unique_hydrogen_keys.append(key)
        unique_hydrogen_indices = [
            int(x.split("_")[1]) for x in unique_hydrogen_keys
        ]
        if not only_sp3_carbon_hydrogen:
            return unique_hydrogen_keys, unique_hydrogen_indices
        else:
            sp3_carbon_unique_hydrogen_indices = []
            sp3_carbon_unique_hydrogen_keys = []
            for key in unique_hydrogen_keys:
                if hydrogen_atoms_dict[key]["H_connected_atom"] != "":
                    connected_atom = hydrogen_atoms_dict[key][
                        "H_connected_atom"].split("_")[0]
                    connected_atom_type = hydrogen_atoms_dict[key][
                        "H_connected_atom"].split("_")[2]
                    if connected_atom == "C" and connected_atom_type == "c3":
                        sp3_carbon_unique_hydrogen_indices.append(
                            int(key.split("_")[1]))
                        sp3_carbon_unique_hydrogen_keys.append(key)
            return sp3_carbon_unique_hydrogen_keys, sp3_carbon_unique_hydrogen_indices

    def _update_equiv_hydrogen_dissociation_energy(
            self, unique_hydrogen_dissociation_energies, unique_hydrogen_keys,
            hydrogen_atoms_dict):
        """
        this function is used to add the dissociation energy of equivalent hydrogen atoms in the hydrogen_atoms_dict

        :param unique_hydrogen_dissociation_energies: The bond dissociation energies of the unique hydrogen atoms.
        :param unique_hydrogen_keys: The keys of the unique hydrogen atoms.
        :param hydrogen_atoms_dict: The dictionary containing hydrogen atom information.

        :return: tuple. The updated hydrogen_atoms_dict, a list of tuple, like [(bond dissociation energy, coordinates),(bond dissociation energy, coordinates)].
        """
        au2kcal = 627.509
        au2kj = 2625.5
        hydrogen_bdes_kj_coords = []
        for i in range(len(unique_hydrogen_keys)):
            key = unique_hydrogen_keys[i]
            energy_au = unique_hydrogen_dissociation_energies[i]
            equiv_group = hydrogen_atoms_dict[key]["equiv_group"]
            for j in equiv_group:
                hydrogen_atoms_dict[
                    "H_" + str(j)]["dissociation_energy_au"] = energy_au
                hydrogen_atoms_dict["H_" +
                                    str(j)]["dissociation_energy_kcal"] = (
                                        energy_au * au2kcal)
                hydrogen_atoms_dict["H_" + str(j)]["dissociation_energy_kj"] = (
                    energy_au * au2kj)
                hydrogen_bdes_kj_coords.append(
                    (energy_au * au2kj,
                     hydrogen_atoms_dict["H_" + str(j)]["coord"]))
        return (hydrogen_atoms_dict, hydrogen_bdes_kj_coords)

    def _print_hydrogen_bond_dissociation_energy(self,
                                                 hydrogen_atoms_dict,
                                                 unit="kcal"):
        """
        this function is used to print the bond dissociation energy of hydrogen atoms, with targeting energy units
        """
        self.ostream.print_info("-" * 50)
        self.ostream.print_info("bond dissociation energy of hydrogen atoms:")
        self.ostream.flush()
        for key in hydrogen_atoms_dict.keys():
            if "dissociation_energy_au" in hydrogen_atoms_dict[key].keys():
                if unit == "kcal":
                    self.ostream.print_info(
                        f"{key} {round(hydrogen_atoms_dict[key]['dissociation_energy_kcal'], 1)} kcal/mol"
                    )
                    self.ostream.flush()
                elif unit == "kj":
                    self.ostream.print_info(
                        f"{key} {round(hydrogen_atoms_dict[key]['dissociation_energy_kj'], 1)} kj/mol"
                    )
                    self.ostream.flush()
                elif unit == "au":
                    self.ostream.print_info(
                        f"{key} {hydrogen_atoms_dict[key]['dissociation_energy_au']} au"
                    )
                    self.ostream.flush()

    def _compute_whole_mol_scf_energy(self, molecule):
        """
        this function is used to compute the whole molecule SCF energy with optimization and single-point energy calculation
        close shell with ScfRestrictedDriver

        :param molecule: The molecule object.

        :return: The optimized molecule object.

        """
        self.ostream.print_info(
            "Optimizing geometry of the molecule before removing hydrogens")
        self.ostream.flush()
        basis_set1 = MolecularBasis.read(molecule, self.basis_sets[0])
        basis_set2 = MolecularBasis.read(molecule, self.basis_sets[1])
        scf_results = self.mol_scf_drv.compute(molecule, basis_set1)
        opt_results = self.mol_opt_drv.compute(molecule, basis_set1,
                                               scf_results)

        opt_molecule = Molecule.read_xyz_string(opt_results["final_geometry"])

        # final energy
        final_single_point_scf_result = self.mol_final_scf_drv.compute(
            opt_molecule, basis_set2)
        self.whole_mol_single_point_scf_energy = self.mol_final_scf_drv.get_scf_energy(
        )
        return opt_molecule

    def _compute_hydrogen_radical_scf_energy(self):
        """
        this function is used to compute the hydrogen radical SCF energy
        open shell with ScfUnrestrictedDriver
        """

        hydrogen = Molecule.read_str(""" H 0.0 0.0 0.0 """)
        hydrogen.set_multiplicity(2)
        basis_set2 = MolecularBasis.read(hydrogen, self.basis_sets[1])
        scf_resultsH = self.hydrogen_final_scf_drv.compute(hydrogen, basis_set2)
        self.hydrogen_single_point_scf_energy = self.hydrogen_final_scf_drv.get_scf_energy(
        )

    def _remove_atom_by_idx(self, mol, atom_indices_to_remove, carbon_indices):
        """
        this function is used to remove atoms from a molecule by using their indices and find the exposed carbon for guess

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

    def _compute_mol_rad_scf_energy(self, mol, radical_carbon_idx, run_idx):
        """
        this function is used to compute the SCF energy for the given radical molecule 
        optimization and single point energy calculation.
        open shell with ScfUnrestrictedDriver

        :param mol: The molecule object.
        :param radical_carbon_idx: The index of the radical carbon atom for guess.
        :param run_idx: The index of the current radical molecule among all radical molecules. can be removed because it is only used to save file

        :return: The SCF energy of the radical molecule.
        """
        step_start = time.time()
        mol.set_multiplicity(2)
        basis_set1 = MolecularBasis.read(mol, self.basis_sets[0])
        basis_set2 = MolecularBasis.read(mol, self.basis_sets[1])
        self.radical_scf_drv.filename = f'bde_{run_idx+1}'
        self.radical_scf_drv.guess_unpaired_electrons = f'{radical_carbon_idx+1}(1.0)'
        scf_resultsmol = self.radical_scf_drv.compute(mol, basis_set1)
        opt_results_rad = self.radical_opt_drv.compute(mol, basis_set1,
                                                       scf_resultsmol)
        mol = Molecule.read_xyz_string(opt_results_rad["final_geometry"])
        mol.set_multiplicity(2)
        self.radical_final_scf_drv.guess_unpaired_electrons = f'{radical_carbon_idx+1}(1.0)'
        scf_results_rad_big = self.radical_final_scf_drv.compute(
            mol, basis_set2)

        step_end = time.time()
        self.ostream.print_info("-" * 50)
        self.ostream.print_info(
            f"time cost : {round(step_end - step_start, 2)} seconds")
        self.ostream.flush()
        radical_single_point_scf_energy = self.radical_final_scf_drv.get_scf_energy(
        )
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
            The list of tuples of (bde, coordinates) of hydrogen atoms.
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
            #add bde based on coords and unique_bde_au
            for i in range(len(bdes_coords)):
                #bde coords is a list of tuple [(bde, (x, y, z)),(bde, (x, y, z))]
                bde_kj = round(bdes_coords[i][0], 1)
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

    def _compute_single_molecule(self, whole_molecule):
        """
        used to compute hydrogen bond dissociation energies (BDEs) for a single molecule.
        """
        self.check_scf_mute()
        self.check_functionals()
        hydrogen_atoms_dict = self._atoms_analyzer(whole_molecule)
        unique_hydrogen_keys, unique_hydrogen_indices = self._fetch_unique_H(
            hydrogen_atoms_dict,
            use_equiv=self.use_equiv,
            only_sp3_carbon_hydrogen=self.only_sp3_carbon_hydrogen)

        #self.ostream.print_info(f'unique_hydrogen_indices {unique_hydrogen_indices}')
        carbon_indices = []
        conn_mat = whole_molecule.get_connectivity_matrix()
        for x in unique_hydrogen_indices:
            assert list(conn_mat[x]).count(1) == 1
            carbon_indices.append(list(conn_mat[x]).index(1))
        #self.ostream.print_info(f'connected carbon_indices {carbon_indices}')
        opt_whole_molecule = self._compute_whole_mol_scf_energy(whole_molecule)
        self._compute_hydrogen_radical_scf_energy()
        self.opt_whole_molecule = opt_whole_molecule
        molecules_rads, radical_carbon_indices = self._remove_atom_by_idx(
            opt_whole_molecule, unique_hydrogen_indices, carbon_indices)

        unique_BDEs_au = []
        count = 1
        for run_idx, (mol_rad, radical_carbon_idx) in enumerate(
                zip(molecules_rads, radical_carbon_indices)):
            self.ostream.print_info(
                f"Computing energy of structure : {count} of {len(molecules_rads)}"
            )
            self.ostream.flush()
            mol_rad_scf_energy = self._compute_mol_rad_scf_energy(
                mol_rad, radical_carbon_idx, run_idx)
            bde_au = mol_rad_scf_energy - self.whole_mol_single_point_scf_energy + self.hydrogen_single_point_scf_energy
            unique_BDEs_au.append(bde_au)
            count += 1
        self.unique_BDEs_au = unique_BDEs_au
        self.unique_hydrogen_keys = unique_hydrogen_keys
        # loop the unique_hydrogen_indices to remove the H atoms from the molecule and calulate the dissciation energy but save the energy for all equivalent H atoms
        # print the dissociation energy for each H atom
        self.hydrogen_atoms_dict, self.bdes_coords = self._update_equiv_hydrogen_dissociation_energy(
            unique_BDEs_au, unique_hydrogen_keys, hydrogen_atoms_dict)
        self._print_hydrogen_bond_dissociation_energy(self.hydrogen_atoms_dict,
                                                      unit=self.energy_unit)

    def compute(self, mol_list: list):
        """
        for computing hydrogen bond dissociation energies (BDEs) for a list of molecules.
        can accept one molecule object as well.
        will save self.mols_bdes_list
        """
        self.mols_bdes_list = []
        if not isinstance(mol_list, list):
            mol_list = [mol_list]
        for mol in mol_list:
            self._compute_single_molecule(mol)
            #use a dictionary to store all bde and hydrogen atom information
            mol_bde_dict = {
                'hydrogen_atoms_dict': self.hydrogen_atoms_dict,
                'bdes_coords': self.bdes_coords,
                'unique_hydrogen_keys': self.unique_hydrogen_keys,
                'unique_BDEs_au': self.unique_BDEs_au,
                'opt_whole_molecule': self.opt_whole_molecule,
            }
            self.mols_bdes_list.append(mol_bde_dict)

    def show(self,
             atom_indices=False,
             atom_labels=False,
             width=400,
             height=300):
        """
        this function is to visualize the hydrogen bond dissociation energies (BDEs) on the hydrogen atoms in the molecule.
        allow atom indices and labels to be shown.
        """
        for item in self.mols_bdes_list:
            mol = item['opt_whole_molecule']
            bdes_coords = item['bdes_coords']
            self._show_bde_on_atom(mol,
                                   width=width,
                                   height=height,
                                   atom_indices=atom_indices,
                                   atom_labels=atom_labels,
                                   one_indexed=True,
                                   bdes_coords=bdes_coords)
