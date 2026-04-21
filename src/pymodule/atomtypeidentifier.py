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
from collections import defaultdict
import numpy as np
import networkx as nx
import sys
import re

from .veloxchemlib import mpi_master
from .veloxchemlib import chemical_element_identifier
from .outputstream import OutputStream
from .mathutils import safe_arccos


class AtomTypeIdentifier:
    """
    A class to identify atom types in a molecule using GAFF (General Amber
    Force Field) atom types based on a VeloxChem molecule object.

    The class processes a molecule object containing the atomic coordinates of
    a molecule to determine the types of atoms according to the GAFF. It
    involves several steps including reading the file, determining covalent
    radii, creating a connectivity matrix, identifying cyclic structures, and
    assigning atom types.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - atomic_symbols: A list of atomic symbols for each atom in the
          molecule.
        - coordinates: A 2D numpy array of atomic coordinates for each atom in
          the molecule.
        - covalent_radii: A list of covalent radii for each atom in the
          molecule.
        - connectivity_matrix: A 2D numpy array indicating which atoms are
          bonded.
        - distance_matrix: A 2D numpy array containing the distances between
          atoms.
        - cyclic_atoms: A set of atom indices that are part of a cyclic
          structure.
        - aromaticity: A list of aromaticity classifications for each cyclic
          structure.
        - atom_cycle_info: A dictionary containing cycle information for each
          atom.
        - atom_info_dict: A dictionary containing detailed information for each
          atom.
        - atom_types_dict: A dictionary containing the GAFF atom types for each
          atom.

    OPLS-AA type annotations
    ------------------------
    The 'opls' key in each atom type dictionary entry stores the corresponding
    OPLS-AA atom type number for reference.  Type numbers follow the canonical
    OPLS-AA 2001 numbering (Jorgensen et al., JACS 1996) extended with the
    commonly used GROMACS oplsaa.ff additions.  Where no standard OPLS-AA
    equivalent exists (e.g. hypervalent sulfur, most phosphorus types, strained
    ring nitrogen, C=S) the value is set to None so that downstream code can
    detect the gap programmatically rather than silently using a wrong type.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the AtomTypeIdentifier instance.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # GAFF version
        self.gaff_version = None

    def plot_connectivity_map_3D(self):
        """
        Plots a 3D connectivity map of the molecule using Matplotlib.
        """

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('Unable to import Matplotlib')

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot bonds based on the connectivity matrix
        for i in range(len(self.atomic_symbols)):
            for j in range(i + 1, len(self.atomic_symbols)):
                if self.connectivity_matrix[i][j] == 1:
                    x = [self.coordinates[i, 0], self.coordinates[j, 0]]
                    y = [self.coordinates[i, 1], self.coordinates[j, 1]]
                    z = [self.coordinates[i, 2], self.coordinates[j, 2]]
                    ax.plot(x, y, z, 'k-', linewidth=2)

        # Add atom labels
        for i, symbol in enumerate(self.atomic_symbols):
            ax.text(self.coordinates[i, 0],
                    self.coordinates[i, 1],
                    self.coordinates[i, 2],
                    f'{symbol} ({i + 1})',
                    fontsize=12,
                    ha='center',
                    va='center')

        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Connectivity Map in 3D')
        plt.show()

    def is_sp2_carbon(self, atom_idx):
        """
        Determines if a C atom, identified by its index, is sp2 hybridized.

        :param atom_idx:
            Index of the atom in the molecule.

        :return:
            True if the atom is sp2 hybridized, False otherwise.
        """

        return (self.atomic_symbols[atom_idx] == 'C'
                and len(list(self.graph.neighbors(atom_idx))) == 3)

    def is_sp2_nitrogen(self, atom_idx):
        """
        Determines if a N atom, identified by its index, is sp2 hybridized.

        :param atom_idx:
            Index of the atom in the molecule.

        :return:
            True if the atom is sp2 hybridized, False otherwise.
        """

        neighbor_identifiers = [
            chemical_element_identifier(self.atomic_symbols[x].upper())
            for x in list(self.graph.neighbors(atom_idx))
        ]

        nonmetal_neighbor_count = 0
        for elem_id in neighbor_identifiers:
            if not self.element_id_is_metal(elem_id):
                nonmetal_neighbor_count += 1

        return (self.atomic_symbols[atom_idx] == 'N'
                and nonmetal_neighbor_count == 2)

    def element_is_metal(self, element):
        return self.element_id_is_metal(chemical_element_identifier(element))

    def element_id_is_metal(self, elem_id):
        return ((elem_id in [3, 4, 11, 12, 13])
                or (19 <= elem_id and elem_id <= 32)
                or (37 <= elem_id and elem_id <= 51)
                or (55 <= elem_id and elem_id <= 84)
                or (87 <= elem_id and elem_id <= 108))

    def detect_closed_cyclic_structures(self):
        """
        Detects closed cyclic structures in a molecule and determines their
        aromaticity.
        """

        self.graph = nx.Graph()
        for i in range(len(self.atomic_symbols)):
            self.graph.add_node(i)
            for j in range(i + 1, len(self.atomic_symbols)):
                if self.connectivity_matrix[i][j] == 1:
                    self.graph.add_edge(i, j)

        all_cycles = list(nx.simple_cycles(self.graph, length_bound=10))

        # Sort by cycle length
        all_cycles = sorted(all_cycles, key=len)

        # Filter only the cycles of size 3 to 10
        filtered_cycles = [
            cycle for cycle in all_cycles if 3 <= len(cycle) <= 10
        ]

        # Remove super-cycles (cycles that contain smaller cycles)
        cycles_to_remove = set()
        for cycle in filtered_cycles:
            for other_cycle in filtered_cycles:
                if (len(cycle) < len(other_cycle)
                        and set(cycle).issubset(set(other_cycle))):
                    cycles_to_remove.add(tuple(other_cycle))

        # Build reduced_cycles by excluding the ones in cycles_to_remove
        self.reduced_cycles = [
            cycle for cycle in filtered_cycles
            if tuple(cycle) not in cycles_to_remove
        ]

        self.cyclic_atoms = set()
        self.aromaticity = []
        self.atom_cycle_info = {}

        # Assignation of aromaticity to all the reduced cycles

        for cycle in self.reduced_cycles:
            self.cyclic_atoms.update(cycle)

            cycle_elements = [self.atomic_symbols[i] for i in cycle]

            cc_bond_distances = []
            for ind, a_i in enumerate(cycle):
                a_j = cycle[(ind + 1) % len(cycle)]
                if (self.atomic_symbols[a_i] == 'C'
                        and self.atomic_symbols[a_j] == 'C'):
                    cc_bond_distances.append(self.distance_matrix[a_i][a_j])

            max_distance, min_distance = 0.0, 0.0
            if cc_bond_distances:
                max_distance = max(cc_bond_distances)
                min_distance = min(cc_bond_distances)

            only_carbon_nitrogen_in_cycle = all(
                [elem in ['C', 'N'] for elem in cycle_elements])

            all_carbons_sp2 = all([
                self.is_sp2_carbon(atom_idx) for atom_idx in cycle
                if self.atomic_symbols[atom_idx] == 'C'
            ])

            all_nitrogens_sp2 = all([
                self.is_sp2_nitrogen(atom_idx) for atom_idx in cycle
                if self.atomic_symbols[atom_idx] == 'N'
            ])

            all_nitrogens_likely_sp2 = True
            for atom_idx in cycle:
                if self.atomic_symbols[atom_idx] == 'N':
                    n_neighbors = list(self.graph.neighbors(atom_idx))
                    if len(n_neighbors) == 2:
                        continue
                    elif len(n_neighbors) == 3:
                        for n_neighbor_idx in n_neighbors:
                            if n_neighbor_idx in cycle:
                                if (self.atomic_symbols[n_neighbor_idx] == 'C'
                                        and self.is_sp2_carbon(n_neighbor_idx)):
                                    continue
                                elif (self.atomic_symbols[n_neighbor_idx] == 'N'
                                      and self.is_sp2_nitrogen(n_neighbor_idx)):
                                    continue
                                else:
                                    all_nitrogens_likely_sp2 = False
                                    break
                    else:
                        all_nitrogens_likely_sp2 = False
                        break

            # TODO: add aromaticity detection for rings with 7-10 atoms

            if len(cycle) == 6 and only_carbon_nitrogen_in_cycle:

                if all_carbons_sp2 and all_nitrogens_sp2:
                    if max_distance - min_distance <= 0.08:
                        aro = 'pure_aromatic'
                    else:
                        aro = 'non_pure_aromatic'
                elif all_carbons_sp2 and (max_distance - min_distance <= 0.08):
                    aro = 'non_pure_aromatic'
                elif all_carbons_sp2 and all_nitrogens_likely_sp2:
                    aro = 'non_pure_aromatic'
                else:
                    aro = 'non_aromatic'

            elif len(cycle) == 5 and all_carbons_sp2:

                if 'S' in cycle_elements or 'O' in cycle_elements:
                    has_ss_so_oo_bond = False
                    max_s_or_o_neighbors = 0
                    for atom_idx, elem in zip(cycle, cycle_elements):
                        if elem not in ['S', 'O']:
                            continue
                        s_or_o_neighbors = [
                            self.atomic_symbols[i]
                            for i in list(self.graph.neighbors(atom_idx))
                        ]
                        if 'S' in s_or_o_neighbors or 'O' in s_or_o_neighbors:
                            has_ss_so_oo_bond = True
                        if len(s_or_o_neighbors) > max_s_or_o_neighbors:
                            max_s_or_o_neighbors = len(s_or_o_neighbors)
                    if (not has_ss_so_oo_bond) and max_s_or_o_neighbors == 2:
                        if 'N' in cycle_elements:
                            aro = ('non_pure_aromatic' if
                                   all_nitrogens_likely_sp2 else 'non_aromatic')
                        else:
                            aro = 'non_pure_aromatic'
                    else:
                        aro = 'non_aromatic'
                else:
                    if 'N' in cycle_elements:
                        aro = ('non_pure_aromatic'
                               if all_nitrogens_likely_sp2 else 'non_aromatic')
                    else:
                        aro = 'non_aromatic'

            elif len(cycle) == 4:

                if set(cycle_elements) == {'C'} and all_carbons_sp2:
                    aro = 'non_pure_aromatic'
                else:
                    aro = 'non_aromatic'

            else:

                aro = 'non_aromatic'

            self.aromaticity.append(aro)

            for atom in cycle:
                if atom not in self.atom_cycle_info:
                    self.atom_cycle_info[atom] = {
                        'sizes': [],
                        'aromaticities': []
                    }
                self.atom_cycle_info[atom]['sizes'].append(len(cycle))
                self.atom_cycle_info[atom]['aromaticities'].append(aro)

        # Additional logic for reassignment of aromaticity in special cases
        # where 3 atoms are shared with aromatic rings.
        for index, cycle in enumerate(self.reduced_cycles):
            # Check if all carbons in the cycle are sp2
            all_carbons_sp2 = all(
                self.is_sp2_carbon(atom_idx) for atom_idx in cycle
                if self.atomic_symbols[atom_idx] == "C")
            if (len(cycle) == 5 and all_carbons_sp2
                    and self.aromaticity[index] == 'non_aromatic'):
                count_pure_aromatic_atoms = sum([
                    1 for atom in cycle if 'pure_aromatic' in
                    self.atom_cycle_info[atom]['aromaticities']
                ])
                if count_pure_aromatic_atoms >= 3:
                    self.aromaticity[index] = 'non_pure_aromatic'
                    for atom in cycle:
                        self.atom_cycle_info[atom]['aromaticities'] = [
                            'non_pure_aromatic' if a == 'non_aromatic' else a
                            for a in self.atom_cycle_info[atom]['aromaticities']
                        ]

    def create_atom_info_dict(self):
        """
        Creates a dictionary containing detailed information for each atom in
        the molecule.
        """

        atom_info_dict = {}

        for i, symbol in enumerate(self.atomic_symbols):
            num_connected_atoms = np.sum(self.connectivity_matrix[i])
            connected_atoms_numbers = [
                j + 1 for j in range(len(self.atomic_symbols))
                if self.connectivity_matrix[i][j] == 1
            ]
            connected_atoms_symbols = [
                self.atomic_symbols[j] for j in range(len(self.atomic_symbols))
                if self.connectivity_matrix[i][j] == 1
            ]
            connected_atoms_distances = [
                self.distance_matrix[i][j]
                for j in range(len(self.atomic_symbols))
                if self.connectivity_matrix[i][j] == 1
            ]

            info = {
                "AtomicSymbol": symbol,
                "AtomNumber": i + 1,
                "NumConnectedAtoms": num_connected_atoms,
                "ConnectedAtoms": connected_atoms_symbols,
                "ConnectedAtomsNumbers": connected_atoms_numbers,
                "ConnectedAtomsDistances": connected_atoms_distances
            }

            if i in self.cyclic_atoms:
                info["CyclicStructure"] = "cycle"
                info["CycleSize"] = self.atom_cycle_info[i]['sizes']
                info["Aromaticity"] = self.atom_cycle_info[i]['aromaticities']
                info["CycleNumber"] = [
                    c_ind for c_ind, c in enumerate(self.reduced_cycles)
                    if i in c
                ]
            else:
                info["CyclicStructure"] = "none"

            atom_info_dict[i + 1] = info
        return atom_info_dict

    def decide_atom_type(self, atom_info_dict):
        """
        Analyzes the molecular structure information to assign atom types to
        each atom in the molecule.
        """

        atom_types_dict = {}

        using_gaff_220 = False
        if self.gaff_version is not None:
            gaff_version_major = self.gaff_version.split('.')[0]
            gaff_version_minor = self.gaff_version.split('.')[1]
            # here we only compare the first digit of gaff_version_minor...
            if gaff_version_minor:
                gaff_version_minor = gaff_version_minor[0]
            if gaff_version_major.isdigit() and gaff_version_minor.isdigit():
                using_gaff_220 = (int(gaff_version_major) >= 2
                                  and int(gaff_version_minor) >= 2)

        for atom_number, info in atom_info_dict.items():
            # Chemical environment information
            connected_symbols = set(info['ConnectedAtoms'])
            name = f"{info['AtomicSymbol']}{info['AtomNumber']}"
            # Carbon type decision
            atom_type = None

            if info['AtomicSymbol'] == 'C':
                atom_type = self.assign_carbon_type(atom_info_dict,
                                                    using_gaff_220, info,
                                                    connected_symbols)

                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    hydrogen_type = None
                    connected_atom_info, hydrogen_type = self.assign_hydrogen_carbon_type(
                        atom_info_dict, atom_type, connected_atom_number, info)

                    if hydrogen_type is not None:
                        H_name = f"H{connected_atom_info['AtomNumber']}"
                        # TODO: rewrite to avoid over-condensed one-liner
                        atom_types_dict.setdefault(H_name,
                                                   {}).update(hydrogen_type)

            elif info['AtomicSymbol'] == 'O':
                atom_type = self.assign_oxygen_type(atom_info_dict, info,
                                                    connected_symbols)

                # Hydrogen type assignment based on oxygen type
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    hydrogen_type = None
                    connected_atom_info, hydrogen_type = self.assign_hydrogen_oxygen_type(
                        atom_info_dict, connected_atom_number, atom_type)

                    if hydrogen_type is not None:
                        H_name = f"H{connected_atom_info['AtomNumber']}"
                        atom_types_dict.setdefault(H_name,
                                                   {}).update(hydrogen_type)

            elif info['AtomicSymbol'] == 'N':

                atom_type = self.assign_nitrogen_type(atom_info_dict, info,
                                                      connected_symbols)

                # Hydrogen type assignment based on nitrogen type
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    hydrogen_type = None
                    connected_atom_info, hydrogen_type = self.assign_hydrogen_nitrogen_type(
                        atom_info_dict, connected_atom_number)

                    if hydrogen_type is not None:
                        H_name = f"H{connected_atom_info['AtomNumber']}"
                        atom_types_dict.setdefault(H_name,
                                                   {}).update(hydrogen_type)

            elif info['AtomicSymbol'] == 'P':

                # Disclaimer: Only non-cyclic structures and non-conjugated cases are considered for now.
                # TODO: Add cyclic and conjugated cases.
                if info.get('CyclicStructure') == 'none':

                    atom_type = self.assign_phosphorus_type(
                        info, connected_symbols)

                    # Hydrogens based on the phosphorus type
                    for connected_atom_number in info['ConnectedAtomsNumbers']:
                        hydrogen_type = None
                        connected_atom_info, hydrogen_type = self.assign_hydrogen_phosphorus_type(
                            atom_info_dict, connected_atom_number, atom_type)

                        if hydrogen_type is not None:
                            H_name = f"H{connected_atom_info['AtomNumber']}"
                            atom_types_dict.setdefault(H_name,
                                                       {}).update(hydrogen_type)

            elif info['AtomicSymbol'] == 'S':

                atom_type = self.assign_sulfur_type(atom_info_dict, info,
                                                    connected_symbols)

                # Hydrogen assignment in the case of thiols
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    hydrogen_type = None
                    connected_atom_info, hydrogen_type = self.assign_hydrogen_sulfur_type(
                        atom_info_dict, connected_atom_number, atom_type)

                    if hydrogen_type is not None:
                        H_name = f"H{connected_atom_info['AtomNumber']}"
                        atom_types_dict.setdefault(H_name,
                                                   {}).update(hydrogen_type)

            # Decision for halogens
            # Halogens (F, Cl, Br, I): OPLS-AA has specific types but they depend
            # on the bonding context (e.g. opls_135 for CH2Cl2, opls_145 for
            # chlorobenzene).  The generic gaff type is used here; a full OPLS
            # halogen mapping requires context that is not yet implemented.
            elif info['AtomicSymbol'] in ['F', 'Cl', 'Br', 'I']:
                atom_type = {'opls': None, 'gaff': info['AtomicSymbol'].lower()}

            uff_type = {'uff': info['AtomicSymbol']}
            if atom_type is not None:  # An atom type was found in the decision tree
                atom_type.update(uff_type)
                atom_types_dict[name] = atom_type
            else:
                atom_types_dict.setdefault(name, {}).update(uff_type)

        # Sort the atom_types_dict into the order matching the input molecule
        # atomtypeidentifier.identify_equivalences and mmforcefieldgenerator.populate_atoms rely on this ordering
        atom_types_dict = self.sort_atom_types(atom_types_dict)
        return atom_types_dict

    def sort_atom_types(self, atom_types_dict):
        """Sorts the entries of the atom_types_dict based on the number occurring in the key
        """
        new_dict = {}
        keys = list(atom_types_dict.keys())
        sorted_keys = sorted(keys, key=self.get_atom_number)
        for key in sorted_keys:
            new_dict[key] = atom_types_dict[key]

        return new_dict

    @staticmethod
    def get_atom_number(atom_type_str):
        """
        Extracts the numeric part from an atom type string.

        This static method uses regular expression to find the first sequence
        of digits in the atom type string, which typically represents the atom
        number.

        :param atom_type_str:
            The atom type string containing a numeric suffix.

        :return:
            The numeric part extracted from the atom type string. Returns 0 if
            no number is found.
        """

        at_match = re.search(r'\d+', atom_type_str)

        return (int(at_match.group()) if at_match else 0)

    def assign_hydrogen_sulfur_type(self, atom_info_dict, connected_atom_number,
                                    sulfur_type):

        hydrogen_type = None

        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):
            if sulfur_type == {'opls': 'opls_200', 'gaff': 'sh'}:
                # opls_204: H on thiol sulfur (Jorgensen OPLS-AA)
                hydrogen_type = {'opls': 'opls_204', 'gaff': 'hs'}

        return connected_atom_info, hydrogen_type

    def assign_sulfur_type(self, atom_info_dict, info, connected_symbols):

        sulfur_type = None

        # S with one connected atom
        if info['NumConnectedAtoms'] == 1:
            # Terminal S (e.g. thioketone C=S): no standard OPLS-AA equivalent
            sulfur_type = {'opls': None, 'gaff': 's'}

        # S with two connected atoms
        elif info['NumConnectedAtoms'] == 2:

            if 'H' in connected_symbols:
                # opls_200: thiol sulfur –SH (Jorgensen OPLS-AA)
                sulfur_type = {'opls': 'opls_200', 'gaff': 'sh'}

            elif all([
                    atom_info_dict[num]['AtomicSymbol'] in ['C', 'N', 'S', 'O']
                    for num in info['ConnectedAtomsNumbers']
            ]):
                # Thio-ether or thio-ester: opls_202 (Jorgensen OPLS-AA)
                sulfur_type = {'opls': 'opls_202', 'gaff': 'ss'}

            else:
                # S2 (disulfide-like): opls_203 (Jorgensen OPLS-AA)
                sulfur_type = {'opls': 'opls_203', 'gaff': 's2'}

        # S with three connected atoms
        elif info['NumConnectedAtoms'] == 3:

            has_sp2_carbon = any([
                atom_info_dict[num]['AtomicSymbol'] == 'C'
                and atom_info_dict[num]['NumConnectedAtoms'] == 3
                for num in info['ConnectedAtomsNumbers']
            ])
            has_sp2_nitrogen = any([
                atom_info_dict[num]['AtomicSymbol'] == 'N'
                and atom_info_dict[num]['NumConnectedAtoms'] == 2
                for num in info['ConnectedAtomsNumbers']
            ])

            if has_sp2_carbon or has_sp2_nitrogen:
                # sx: S bonded to sp2 C/N; no standard OPLS-AA equivalent
                sulfur_type = {'opls': None, 'gaff': 'sx'}

            else:
                # s4: hypervalent S with 3 substituents; no standard OPLS-AA equivalent
                sulfur_type = {'opls': None, 'gaff': 's4'}

        # S with four connected atoms
        elif info['NumConnectedAtoms'] == 4:

            if any(atom_info_dict[num]['AtomicSymbol'] == 'C'
                   and atom_info_dict[num]['NumConnectedAtoms'] == 3
                   for num in info['ConnectedAtomsNumbers']):
                # sy: S bonded to sp2 C with 4 total bonds; no standard OPLS-AA equivalent
                sulfur_type = {'opls': None, 'gaff': 'sy'}

            else:
                # s6: hexavalent S (e.g. SF6, DMSO2); no standard OPLS-AA equivalent
                sulfur_type = {'opls': None, 'gaff': 's6'}

            # TODO: Sp3 S connected with hydrogen

        return sulfur_type

    def assign_hydrogen_phosphorus_type(self, atom_info_dict,
                                        connected_atom_number, phosphorus_type):

        if phosphorus_type is None:
            self.ostream.print_warning(
                f"Phosphorus type is None for atom {connected_atom_number}, "
                "assigning hp to connected hydrogen.")

        hydrogen_type = None

        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):
            # opls_440: H on phosphorus (GROMACS oplsaa.ff extension for P–H)
            hydrogen_type = {'opls': 'opls_440', 'gaff': 'hp'}

        return connected_atom_info, hydrogen_type

    def assign_phosphorus_type(self, info, connected_symbols):

        phosphorus_type = None

        # Hypervalent phosphorus, 4 substituents (e.g. phosphate ester)
        if info['NumConnectedAtoms'] == 4:
            connected_to_metal = any([
                self.element_is_metal(elem.upper())
                for elem in connected_symbols
            ])
            if 'O' in connected_symbols or connected_to_metal:
                # p5: no standard OPLS-AA equivalent
                phosphorus_type = {'opls': None, 'gaff': 'p5'}

        # sp3 phosphorus, 3 substituents
        elif info['NumConnectedAtoms'] == 3:
            if 'O' in connected_symbols:
                oxygen_count = info['ConnectedAtoms'].count('O')
            else:
                oxygen_count = 0

            if oxygen_count == 0:
                # p3: regular sp3 P (e.g. phosphine PR3); no standard OPLS-AA equivalent
                phosphorus_type = {'opls': None, 'gaff': 'p3'}
            else:
                # p4: hypervalent P with 3 substituents inc. O; no standard OPLS-AA equivalent
                phosphorus_type = {'opls': None, 'gaff': 'p4'}

        # sp2 phosphorus (C=P etc.)
        elif info['NumConnectedAtoms'] == 2:
            # p2: no standard OPLS-AA equivalent
            phosphorus_type = {'opls': None, 'gaff': 'p2'}

        return phosphorus_type

    def assign_hydrogen_nitrogen_type(self, atom_info_dict,
                                      connected_atom_number):

        hydrogen_type = None

        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):
            # opls_240: H on N (amine/amide N–H; Jorgensen OPLS-AA)
            hydrogen_type = {'opls': 'opls_240', 'gaff': 'hn'}

        return connected_atom_info, hydrogen_type

    def assign_hydrogen_oxygen_type(self, atom_info_dict, connected_atom_number,
                                    oxygen_type):

        hydrogen_type = None

        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):

            if oxygen_type == {'opls': 'opls_111', 'gaff': 'ow'}:
                # opls_112: H in TIP3P/SPC water
                hydrogen_type = {'opls': 'opls_112', 'gaff': 'hw'}

            elif oxygen_type == {'opls': 'opls_154', 'gaff': 'oh'}:
                # opls_155: H on alcohol or phenol O–H
                hydrogen_type = {'opls': 'opls_155', 'gaff': 'ho'}

        return connected_atom_info, hydrogen_type

    def assign_nitrogen_type(self, atom_info_dict, info, connected_symbols):

        nitrogen_type = None

        connected_atoms = info['ConnectedAtoms']
        connected_atoms_numbers = info['ConnectedAtomsNumbers']

        num_hydrogens = sum([1 for symbol in connected_atoms if symbol == 'H'])

        # Non-cyclic

        if info.get('CyclicStructure') == 'none':
            if info['NumConnectedAtoms'] == 4:
                if num_hydrogens == 4:
                    # n+: NH4+; opls_287 (Jorgensen OPLS-AA)
                    nitrogen_type = {'opls': 'opls_287', 'gaff': 'n+'}

                elif num_hydrogens == 3:
                    # nz: RNH3+ (protonated primary amine); opls_288
                    nitrogen_type = {'opls': 'opls_288', 'gaff': 'nz'}

                elif num_hydrogens == 2:
                    # ny: R2NH2+ (protonated secondary amine); opls_289
                    nitrogen_type = {'opls': 'opls_289', 'gaff': 'ny'}

                elif num_hydrogens == 1:
                    # nx: R3NH+ (protonated tertiary amine); opls_290
                    nitrogen_type = {'opls': 'opls_290', 'gaff': 'nx'}

                else:
                    # n4: NR4+ (quaternary N, no H); opls_291
                    nitrogen_type = {'opls': 'opls_291', 'gaff': 'n4'}

            elif info['NumConnectedAtoms'] == 3:
                # This case is highly dependent on the environment of
                # the nitrogen.
                # Create flags to check for specific cases.
                found_nitro = False
                found_amide = False
                found_aromatic = False
                found_sp2_carbon = False

                if (sorted(info['ConnectedAtoms']) in [
                        sorted(['C', 'O', 'O']),
                        sorted(['N', 'O', 'O']),
                        sorted(['O', 'O', 'O']),
                ]):
                    found_nitro = True

                elif 'C' in connected_symbols:
                    for idx, atom in enumerate(connected_atoms_numbers):
                        if connected_atoms[idx] != 'C':
                            continue

                            # Check for aromaticity
                        if ('cycle' in atom_info_dict[atom]['CyclicStructure']
                                and 'pure_aromatic'
                                in self.atom_info_dict[atom]['Aromaticity']):
                            found_aromatic = True

                        elif atom_info_dict[atom]['NumConnectedAtoms'] == 3:
                            # Check for amide
                            for connected_to_carbon in atom_info_dict[atom][
                                    'ConnectedAtomsNumbers']:
                                atom_symbol = atom_info_dict[
                                    connected_to_carbon]['AtomicSymbol']
                                atom_connectivity = atom_info_dict[
                                    connected_to_carbon]['NumConnectedAtoms']

                                if ((atom_symbol == 'O'
                                     and atom_connectivity == 1)
                                        or (atom_symbol == 'S'
                                            and atom_connectivity == 1)):
                                    found_amide = True
                                    break

                            if not found_amide:
                                found_sp2_carbon = True

                # Assign nitrogen types based on the flags using the hierarchy:
                # 1. Nitro  2. Amide  3. Aromatic / sp2 carbon

                if found_nitro:
                    # no: nitro N (–NO2); opls_760 (Jorgensen OPLS-AA)
                    nitrogen_type = {'opls': 'opls_760', 'gaff': 'no'}

                elif found_amide:
                    if num_hydrogens == 1:
                        # ns: secondary amide N (1H, –CONH–); opls_241
                        nitrogen_type = {'opls': 'opls_241', 'gaff': 'ns'}

                    elif num_hydrogens == 2:
                        # nt: primary amide N (2H, –CONH2); opls_237
                        nitrogen_type = {'opls': 'opls_237', 'gaff': 'nt'}

                    else:
                        # n: tertiary amide N (0H, –CON<); opls_238
                        nitrogen_type = {'opls': 'opls_238', 'gaff': 'n'}

                elif found_aromatic or found_sp2_carbon:
                    if num_hydrogens == 1:
                        # nu: amine N adjacent to sp2 C, 1H (aniline-like); opls_901
                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nu'}

                    elif num_hydrogens == 2:
                        # nv: amine N adjacent to sp2 C, 2H; opls_901
                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nv'}

                    else:
                        # nh: amine N adjacent to sp2 C, 0H; opls_901
                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nh'}

                else:
                    if num_hydrogens == 1:
                        # n7: sp3 N with 1H; opls_300
                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n7'}

                    elif num_hydrogens == 2:
                        # n8: sp3 N with 2H; opls_300
                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n8'}

                    else:
                        # n3: sp3 N with 0H; opls_300
                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n3'}

            elif info['NumConnectedAtoms'] == 2:
                has_sp1_nitrogen = any([
                    atom_info_dict[atom]['AtomicSymbol'] == 'N'
                    and atom_info_dict[atom]['NumConnectedAtoms'] == 1
                    for atom in connected_atoms_numbers
                ])

                sp2_carbon_count = sum(
                    1 for num in connected_atoms_numbers
                    if atom_info_dict[num]['AtomicSymbol'] == 'C'
                    and atom_info_dict[num]['NumConnectedAtoms'] == 3)
                sp1_carbon_count = sum(
                    1 for num in connected_atoms_numbers
                    if atom_info_dict[num]['AtomicSymbol'] == 'C'
                    and atom_info_dict[num]['NumConnectedAtoms'] == 2)
                n2_count = sum(
                    1 for num in connected_atoms_numbers
                    if atom_info_dict[num]['AtomicSymbol'] == 'N'
                    and atom_info_dict[num]['NumConnectedAtoms'] == 2)

                # Check if the Nitrogen is connected to another
                # Nitrogen with sp1 hybridization
                if has_sp1_nitrogen:
                    # n1: sp1 N (terminal, nitrile-end-like); opls_753
                    nitrogen_type = {'opls': 'opls_753', 'gaff': 'n1'}

                elif sp2_carbon_count + sp1_carbon_count + n2_count == 2:
                    # ne: sp2 N in conjugated non-ring system (imine C=N–);
                    # opls_296 (closest OPLS-AA equivalent)
                    nitrogen_type = {'opls': 'opls_296', 'gaff': 'ne'}

                else:
                    # Check bond angle to distinguish sp1 from sp2
                    vec_ji = (self.coordinates[connected_atoms_numbers[0] - 1] -
                              self.coordinates[info['AtomNumber'] - 1])
                    vec_jk = (self.coordinates[connected_atoms_numbers[1] - 1] -
                              self.coordinates[info['AtomNumber'] - 1])
                    theta_ijk = safe_arccos(
                        np.dot(vec_ji, vec_jk) /
                        (np.linalg.norm(vec_ji) * np.linalg.norm(vec_jk)))

                    if theta_ijk * 180.0 / np.pi > 170.0:
                        # Near-linear geometry → sp1 N; opls_753
                        nitrogen_type = {'opls': 'opls_753', 'gaff': 'n1'}
                    else:
                        # n2: sp2 N, non-aromatic, 2-connected; opls_531
                        nitrogen_type = {'opls': 'opls_531', 'gaff': 'n2'}

            elif info['NumConnectedAtoms'] == 1:
                # n1: terminal sp1 N (e.g. isocyanide); opls_753
                nitrogen_type = {'opls': 'opls_753', 'gaff': 'n1'}

        # Cyclic

        elif info.get('CyclicStructure') == 'cycle':
            if (info['NumConnectedAtoms'] == 2
                    and 'pure_aromatic' in info.get('Aromaticity')):
                # nb: sp2 N in pure aromatic ring (pyridine-type); opls_520
                nitrogen_type = {'opls': 'opls_520', 'gaff': 'nb'}

            elif (info['NumConnectedAtoms'] == 2
                  and 'non_pure_aromatic' in info.get('Aromaticity')):
                # nc: sp2 N in non-pure aromatic ring; opls_520 (best match)
                nitrogen_type = {'opls': 'opls_520', 'gaff': 'nc'}

            elif (info['NumConnectedAtoms'] == 3
                  and 'pure_aromatic' in info.get('Aromaticity')):
                # nb: pyridine N as ligand in organometallic; opls_520
                nitrogen_type = {'opls': 'opls_520', 'gaff': 'nb'}

            elif (info['NumConnectedAtoms'] == 3
                  and 'non_pure_aromatic' in info.get('Aromaticity')):
                # na: pyrrole-type N with H in non-pure aromatic ring;
                # opls_534 (GROMACS oplsaa.ff pyrrole N)
                nitrogen_type = {'opls': 'opls_534', 'gaff': 'na'}

                if 'C' in connected_symbols:
                    # Check for amides / sulfamides via connected carbons
                    found_CO = False

                    for idx, atom in enumerate(connected_atoms_numbers):
                        if connected_atoms[idx] != 'C':
                            continue

                        for connected_to_carbon in atom_info_dict[atom][
                                'ConnectedAtomsNumbers']:
                            atom_symbol = atom_info_dict[connected_to_carbon][
                                'AtomicSymbol']
                            atom_connectivity = atom_info_dict[
                                connected_to_carbon]['NumConnectedAtoms']

                            if ((atom_symbol == 'O' and atom_connectivity == 1)
                                    or (atom_symbol == 'S'
                                        and atom_connectivity == 1)):
                                if num_hydrogens == 1:
                                    # ns: cyclic secondary amide N; opls_241
                                    nitrogen_type = {
                                        'opls': 'opls_241',
                                        'gaff': 'ns'
                                    }

                                elif num_hydrogens == 2:
                                    # nt: cyclic primary amide N; opls_237
                                    nitrogen_type = {
                                        'opls': 'opls_237',
                                        'gaff': 'nt'
                                    }

                                else:
                                    # n: cyclic tertiary amide N; opls_238
                                    nitrogen_type = {
                                        'opls': 'opls_238',
                                        'gaff': 'n'
                                    }

                                found_CO = True
                                break

                        if found_CO:
                            break

            # Nitrogens in non-aromatic cycles
            elif (info['NumConnectedAtoms'] == 3
                  and 'non_aromatic' in info.get('Aromaticity')):
                if 'C' in connected_symbols:
                    # Check for amides / sulfamides via connected carbons
                    found_CO = False

                    for idx, atom in enumerate(connected_atoms_numbers):
                        if connected_atoms[idx] != 'C':
                            continue

                        for connected_to_carbon in atom_info_dict[atom][
                                'ConnectedAtomsNumbers']:
                            atom_symbol = atom_info_dict[connected_to_carbon][
                                'AtomicSymbol']
                            atom_connectivity = atom_info_dict[
                                connected_to_carbon]['NumConnectedAtoms']

                            if ((atom_symbol == 'O' and atom_connectivity == 1)
                                    or (atom_symbol == 'S'
                                        and atom_connectivity == 1)):
                                if num_hydrogens == 1:
                                    # ns: secondary amide N; opls_241
                                    nitrogen_type = {
                                        'opls': 'opls_241',
                                        'gaff': 'ns'
                                    }

                                elif num_hydrogens == 2:
                                    # nt: primary amide N; opls_237
                                    nitrogen_type = {
                                        'opls': 'opls_237',
                                        'gaff': 'nt'
                                    }

                                else:
                                    # n: tertiary amide N; opls_238
                                    nitrogen_type = {
                                        'opls': 'opls_238',
                                        'gaff': 'n'
                                    }

                                found_CO = True
                                break

                        if found_CO:
                            break

                    if not found_CO:
                        # Check if the nitrogen is adjacent to a sp2 C or N
                        connected_to_sp2_carbon = any([
                            atom_info_dict[atom]['AtomicSymbol'] == 'C'
                            and atom_info_dict[atom]['NumConnectedAtoms'] == 3
                            for atom in connected_atoms_numbers
                        ])
                        connected_to_sp2_nitrogen = any([
                            atom_info_dict[atom]['AtomicSymbol'] == 'N'
                            and atom_info_dict[atom]['NumConnectedAtoms'] == 2
                            for atom in connected_atoms_numbers
                        ])

                        if connected_to_sp2_carbon or connected_to_sp2_nitrogen:
                            if num_hydrogens == 1:
                                # nu: amine N next to sp2 C, 1H; opls_901
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nu'
                                }

                            elif num_hydrogens == 2:
                                # nv: amine N next to sp2 C, 2H; opls_901
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nv'
                                }

                            else:
                                # nh: amine N next to sp2 C, 0H; opls_901
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nh'
                                }

                        else:
                            # n3/n7/n8 and small-ring special cases
                            if num_hydrogens == 1:
                                if 3 in info.get('CycleSize'):
                                    # n5: N in aziridine (3-ring); no OPLS-AA equivalent
                                    nitrogen_type = {'opls': None, 'gaff': 'n5'}

                                elif 4 in info.get('CycleSize'):
                                    # n6: N in azetidine (4-ring); no OPLS-AA equivalent
                                    nitrogen_type = {'opls': None, 'gaff': 'n6'}

                                else:
                                    # n7: sp3 N, 1H; opls_300
                                    nitrogen_type = {
                                        'opls': 'opls_300',
                                        'gaff': 'n7'
                                    }

                            elif num_hydrogens == 2:
                                # n8: sp3 N, 2H; opls_300
                                nitrogen_type = {
                                    'opls': 'opls_300',
                                    'gaff': 'n8'
                                }

                            else:
                                # n3: sp3 N, 0H; opls_300
                                nitrogen_type = {
                                    'opls': 'opls_300',
                                    'gaff': 'n3'
                                }

                # No C neighbour — check for special RG3/RG4 cases
                else:
                    if num_hydrogens == 1:
                        if 3 in info.get('CycleSize'):
                            # n5: aziridine N; no OPLS-AA equivalent
                            nitrogen_type = {'opls': None, 'gaff': 'n5'}

                        elif 4 in info.get('CycleSize'):
                            # n6: azetidine N; no OPLS-AA equivalent
                            nitrogen_type = {'opls': None, 'gaff': 'n6'}

                        else:
                            # n7: sp3 N, 1H; opls_300
                            nitrogen_type = {'opls': 'opls_300', 'gaff': 'n7'}

                    elif num_hydrogens == 0:
                        if 3 in info.get('CycleSize'):
                            # np: bridging N in 3-ring, no H; no OPLS-AA equivalent
                            nitrogen_type = {'opls': None, 'gaff': 'np'}

                        elif 4 in info.get('CycleSize'):
                            # nq: bridging N in 4-ring, no H; no OPLS-AA equivalent
                            nitrogen_type = {'opls': None, 'gaff': 'nq'}

                        else:
                            # n3: sp3 N, 0H; opls_300
                            nitrogen_type = {'opls': 'opls_300', 'gaff': 'n3'}

            elif (info['NumConnectedAtoms'] == 2
                  and 'non_aromatic' in info.get('Aromaticity')):
                sp2_carbon_count = sum(
                    1 for num in connected_atoms_numbers
                    if atom_info_dict[num]['AtomicSymbol'] == 'C'
                    and atom_info_dict[num]['NumConnectedAtoms'] == 3)

                sp1_carbon_count = sum(
                    1 for num in connected_atoms_numbers
                    if atom_info_dict[num]['AtomicSymbol'] == 'C'
                    and atom_info_dict[num]['NumConnectedAtoms'] == 2)

                n2_count = sum(
                    1 for num in connected_atoms_numbers
                    if atom_info_dict[num]['AtomicSymbol'] == 'N'
                    and atom_info_dict[num]['NumConnectedAtoms'] == 2)

                if sp2_carbon_count + sp1_carbon_count + n2_count == 2:
                    # ne: conjugated sp2 N in non-aromatic ring; opls_296
                    nitrogen_type = {'opls': 'opls_296', 'gaff': 'ne'}

                else:
                    # n2: sp2 N in non-aromatic ring, 2-connected; opls_531
                    nitrogen_type = {'opls': 'opls_531', 'gaff': 'n2'}

        return nitrogen_type

    def assign_carbon_type(self, atom_info_dict, using_gaff_220, info,
                           connected_symbols):

        carbon_type = None

        # Pure aromatic cycles

        if (info.get('CyclicStructure') == 'cycle'
                and 'pure_aromatic' in info.get('Aromaticity')):
            if info['NumConnectedAtoms'] == 3:
                # Check for identifying biphenyls:
                # connected_carbons_in_diff_cycle_and_pure_aromatic
                connected_carbon_atom = None

                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    connected_atom_info = atom_info_dict[connected_atom_number]

                    if connected_atom_info.get('AtomicSymbol') != 'C':
                        continue

                    if not connected_atom_info.get('CycleNumber'):
                        continue

                    if 'pure_aromatic' not in connected_atom_info.get(
                            'Aromaticity'):
                        continue

                    common_cycle_numbers = (
                        set(connected_atom_info.get('CycleNumber'))
                        & set(info.get('CycleNumber')))

                    # Note: Here we count both pure_aromatic and
                    # non_pure_aromatic rings
                    common_aromatic_cycles = [
                        cycle_num for cycle_num in common_cycle_numbers
                        if self.aromaticity[cycle_num] in
                        ['pure_aromatic', 'non_pure_aromatic']
                    ]

                    if not common_aromatic_cycles:
                        connected_carbon_atom = connected_atom_number
                        break

                if connected_symbols == {'C', 'H'}:
                    # ca: aromatic C–H (benzene-like); opls_145
                    carbon_type = {'opls': 'opls_145', 'gaff': 'ca'}

                elif connected_symbols == {'C', 'N', 'H'}:
                    # ca in heteroaromatic ring (e.g. pyridine C adjacent to N);
                    # opls_521 (GROMACS oplsaa.ff pyridine C)
                    carbon_type = {'opls': 'opls_521', 'gaff': 'ca'}

                elif connected_carbon_atom is not None:
                    # cp: aromatic C at biphenyl-type inter-ring junction;
                    # opls_521 (same environment as heteroaromatic C)
                    carbon_type = {'opls': 'opls_521', 'gaff': 'cp'}

                else:
                    # ca: generic aromatic C (no H, no special environment);
                    # opls_145
                    carbon_type = {'opls': 'opls_145', 'gaff': 'ca'}

        # Non-pure aromatic cycles

        elif (info.get('CyclicStructure') == 'cycle'
              and 'non_pure_aromatic' in info.get('Aromaticity')):
            has_terminal_oxygen = False
            has_terminal_sulfur = False
            for connected_atom_number in info['ConnectedAtomsNumbers']:
                if ((atom_info_dict[connected_atom_number]['AtomicSymbol']
                     == 'O') and
                    (atom_info_dict[connected_atom_number]['CyclicStructure']
                     == 'none') and
                    (atom_info_dict[connected_atom_number]['NumConnectedAtoms']
                     == 1)):
                    has_terminal_oxygen = True
                elif (
                    (atom_info_dict[connected_atom_number]['AtomicSymbol']
                     == 'S') and
                    (atom_info_dict[connected_atom_number]['CyclicStructure']
                     == 'none') and
                    (atom_info_dict[connected_atom_number]['NumConnectedAtoms']
                     == 1)):
                    has_terminal_sulfur = True

            if has_terminal_oxygen:
                # c: carbonyl C in non-pure aromatic ring (e.g. lactam C=O);
                # opls_235
                carbon_type = {'opls': 'opls_235', 'gaff': 'c'}

            elif has_terminal_sulfur:
                # cs: C=S in non-pure aromatic ring; no OPLS-AA equivalent
                carbon_type = {'opls': None, 'gaff': 'cs'}

            else:
                # cc: sp2 C in non-pure aromatic ring;
                # opls_142 (generic sp2 alkene C is the closest OPLS-AA type)
                carbon_type = {'opls': 'opls_142', 'gaff': 'cc'}

        # Non-aromatic cycles

        elif (info.get('CyclicStructure') == 'cycle'
              and 'non_aromatic' in info.get('Aromaticity')):
            if info['NumConnectedAtoms'] == 4:
                if 3 in info['CycleSize']:
                    # cx: sp3 C in cyclopropane; opls_352 (Jorgensen OPLS-AA)
                    carbon_type = {'opls': 'opls_352', 'gaff': 'cx'}

                elif 4 in info['CycleSize']:
                    # cy: sp3 C in cyclobutane; opls_356 (Jorgensen OPLS-AA)
                    carbon_type = {'opls': 'opls_356', 'gaff': 'cy'}

                elif 5 in info['CycleSize'] and using_gaff_220:
                    # c5: sp3 C in 5-membered ring; OPLS-AA does not distinguish
                    # 5-ring sp3 C from generic sp3 C → opls_136
                    carbon_type = {'opls': 'opls_136', 'gaff': 'c5'}

                elif 6 in info['CycleSize'] and using_gaff_220:
                    # c6: sp3 C in 6-membered ring; opls_136 (same reasoning)
                    carbon_type = {'opls': 'opls_136', 'gaff': 'c6'}

                else:
                    # c3: generic sp3 C; opls_135
                    carbon_type = {'opls': 'opls_135', 'gaff': 'c3'}

            elif info['NumConnectedAtoms'] == 3:
                has_terminal_oxygen = False
                has_terminal_sulfur = False
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    if ((atom_info_dict[connected_atom_number]['AtomicSymbol']
                         == 'O') and (atom_info_dict[connected_atom_number]
                                      ['NumConnectedAtoms'] == 1)):
                        has_terminal_oxygen = True
                    elif ((atom_info_dict[connected_atom_number]['AtomicSymbol']
                           == 'S') and (atom_info_dict[connected_atom_number]
                                        ['NumConnectedAtoms'] == 1)):
                        has_terminal_sulfur = True

                if has_terminal_oxygen:
                    # c: carbonyl C in non-aromatic ring; opls_235
                    carbon_type = {'opls': 'opls_235', 'gaff': 'c'}

                elif has_terminal_sulfur:
                    # cs: C=S in non-aromatic ring; no OPLS-AA equivalent
                    carbon_type = {'opls': None, 'gaff': 'cs'}

                elif 3 in info['CycleSize']:
                    # cu: sp2 C in cyclopropene-type ring; opls_350
                    carbon_type = {'opls': 'opls_350', 'gaff': 'cu'}

                elif 4 in info['CycleSize']:
                    # cv: sp2 C in cyclobutene-type ring; opls_354
                    carbon_type = {'opls': 'opls_354', 'gaff': 'cv'}

                elif 'C' in connected_symbols or 'N' in connected_symbols:
                    sp2_carbon_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'C'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 3)
                    sp1_carbon_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'C'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 2)
                    n2_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'N'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 2)

                    if (sp2_carbon_count + sp1_carbon_count == 2
                            or sp2_carbon_count + sp1_carbon_count == 3
                            or sp2_carbon_count + sp1_carbon_count + n2_count
                            == 2):
                        # ce: conjugated sp2 C in non-aromatic ring; opls_142
                        carbon_type = {'opls': 'opls_142', 'gaff': 'ce'}

                    else:
                        # c2: generic sp2 C in non-aromatic ring; opls_141
                        carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}

                else:
                    # c2: generic sp2 C in ring; opls_141
                    carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}

            elif info['NumConnectedAtoms'] == 2:
                # c1: sp C in ring; opls_235 (best available OPLS-AA match)
                carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}

        # Non-cyclic

        elif info.get('CyclicStructure') == 'none':
            if info['NumConnectedAtoms'] == 4:
                # c3: generic sp3 C; opls_135
                carbon_type = {'opls': 'opls_135', 'gaff': 'c3'}

            elif info['NumConnectedAtoms'] == 3:
                has_terminal_oxygen = False
                has_terminal_sulfur = False
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    if ((atom_info_dict[connected_atom_number]['AtomicSymbol']
                         == 'O') and (atom_info_dict[connected_atom_number]
                                      ['NumConnectedAtoms'] == 1)):
                        has_terminal_oxygen = True
                    elif ((atom_info_dict[connected_atom_number]['AtomicSymbol']
                           == 'S') and (atom_info_dict[connected_atom_number]
                                        ['NumConnectedAtoms'] == 1)):
                        has_terminal_sulfur = True

                if has_terminal_oxygen:
                    # c: non-cyclic carbonyl C; opls_235
                    carbon_type = {'opls': 'opls_235', 'gaff': 'c'}

                elif has_terminal_sulfur:
                    # cs: non-cyclic C=S; no OPLS-AA equivalent
                    carbon_type = {'opls': None, 'gaff': 'cs'}

                elif 'C' in connected_symbols or 'N' in connected_symbols:
                    sp2_carbon_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'C'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 3)
                    sp1_carbon_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'C'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 2)
                    n2_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'N'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 2)

                    if (sp2_carbon_count + sp1_carbon_count == 2
                            or sp2_carbon_count + sp1_carbon_count == 3
                            or sp2_carbon_count + sp1_carbon_count + n2_count
                            == 2):
                        # ce: conjugated sp2 C (non-ring); opls_142
                        carbon_type = {'opls': 'opls_142', 'gaff': 'ce'}

                    else:
                        # c2: generic sp2 C (non-ring); opls_141
                        carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}

                else:
                    # c2: generic sp2 C (non-ring); opls_141
                    carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}

            elif info['NumConnectedAtoms'] == 2:
                if 'O' in connected_symbols:
                    # c1: C in carbonyl group or acid anhydride; opls_235
                    carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}

                else:
                    sp2_carbon_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'C'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 3)
                    sp1_carbon_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'C'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 2)
                    n_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'N'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 1)
                    n2_count = sum(
                        1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] == 'N'
                        and atom_info_dict[num]['NumConnectedAtoms'] == 2)

                    if sp2_carbon_count == 2:
                        # c1: allenic / cumulene central C; opls_235
                        carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}

                    elif (sp2_carbon_count + sp1_carbon_count + n_count == 2 or
                          sp2_carbon_count + sp1_carbon_count + n2_count == 2):
                        # Distinguish sp vs sp2 by bond length
                        vec_ji = (
                            self.coordinates[info['ConnectedAtomsNumbers'][0] -
                                             1] -
                            self.coordinates[info['AtomNumber'] - 1])
                        vec_jk = (
                            self.coordinates[info['ConnectedAtomsNumbers'][1] -
                                             1] -
                            self.coordinates[info['AtomNumber'] - 1])
                        bond_length_ji = np.linalg.norm(vec_ji)
                        bond_length_jk = np.linalg.norm(vec_jk)

                        if (bond_length_ji <= 1.3475
                                and bond_length_jk <= 1.3475):
                            # c1: sp C (alkyne/nitrile end); opls_235
                            carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}
                        else:
                            # cg: sp C (internal alkyne attached to sp2 C);
                            # opls_157 (GROMACS oplsaa.ff internal alkyne C)
                            carbon_type = {'opls': 'opls_157', 'gaff': 'cg'}

                    else:
                        # c1: sp C (default 2-connected case); opls_235
                        carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}

        return carbon_type

    def assign_hydrogen_carbon_type(self, atom_info_dict, carbon_type,
                                    connected_atom_number, info):

        hydrogen_type = None

        # Assignment for the Hydrogens linked to the carbons
        # ewd = Electron withdrawing atoms

        ewd_atoms = ['N', 'Br', 'Cl', 'I', 'F', 'S', 'O']

        ewd_count = sum(1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] in ewd_atoms)
        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):
            # sp1 carbon
            if carbon_type == {'opls': 'opls_235', 'gaff': 'c1'}:
                # ha: H on sp1 C (terminal alkyne ≡C–H); opls_146
                hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}

            # sp2 carbon — covers all sp2 C types including ring variants
            elif carbon_type == {
                    'opls': 'opls_350',
                    'gaff': 'cu'
            } or carbon_type == {
                    'opls': 'opls_354',
                    'gaff': 'cv'
            } or carbon_type == {
                    'opls': 'opls_141',
                    'gaff': 'c2'
            } or carbon_type == {
                    'opls': 'opls_145',
                    'gaff': 'ca'
            } or carbon_type == {
                    'opls': 'opls_142',
                    'gaff': 'ce'
            } or carbon_type == {
                    'opls': 'opls_142',
                    'gaff': 'cc'
            } or carbon_type == {
                    'opls': 'opls_142',
                    'gaff': 'cf'
            } or carbon_type == {
                    'opls': 'opls_142',
                    'gaff': 'cd'
            } or carbon_type == {
                    'opls': None,
                    'gaff': 'cs'
            } or carbon_type == {
                    'opls': 'opls_235',
                    'gaff': 'c'
            }:
                if ewd_count == 1:
                    # h4: H on sp2 C with 1 EWD neighbour; opls_146
                    hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}

                elif ewd_count == 2:
                    # h5: H on sp2 C with 2 EWD neighbours;
                    # OPLS-AA has no distinct type for this case → opls_146
                    hydrogen_type = {'opls': 'opls_146', 'gaff': 'h5'}

                else:
                    # ha: H on sp2 C (no EWD); opls_146
                    hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}

            # Hydrogens on C in heteroaromatic ring (pyridine C–N–C pattern)
            elif carbon_type == {'opls': 'opls_521', 'gaff': 'ca'}:
                # h4: H on heteroaromatic C adjacent to N; opls_146
                hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}

            # sp3 carbon — covers generic and ring-specific sp3 C types
            elif carbon_type == {
                    'opls': 'opls_135',
                    'gaff': 'c3'
            } or carbon_type == {
                    'opls': 'opls_136',
                    'gaff': 'c5'
            } or carbon_type == {
                    'opls': 'opls_136',
                    'gaff': 'c6'
            } or carbon_type == {
                    'opls': 'opls_352',
                    'gaff': 'cx'
            } or carbon_type == {
                    'opls': 'opls_356',
                    'gaff': 'cy'
            }:
                if ewd_count == 1:
                    # h1: H on sp3 C with 1 EWD neighbour; opls_140
                    hydrogen_type = {'opls': 'opls_140', 'gaff': 'h1'}

                elif ewd_count == 2:
                    # h2: H on sp3 C with 2 EWD neighbours; opls_139
                    hydrogen_type = {'opls': 'opls_139', 'gaff': 'h2'}

                elif ewd_count == 3:
                    # h3: H on sp3 C with 3 EWD neighbours; opls_138
                    hydrogen_type = {'opls': 'opls_138', 'gaff': 'h3'}

                else:
                    # hc: H on sp3 C (no EWD); opls_140
                    hydrogen_type = {'opls': 'opls_140', 'gaff': 'hc'}

        return connected_atom_info, hydrogen_type

    def assign_oxygen_type(self, atom_info_dict, info, connected_symbols):

        oxygen_type = None

        # Non-cyclic

        if info.get('CyclicStructure') == 'none':

            if (info['NumConnectedAtoms'] == 2 and connected_symbols == {'H'}):
                # ow: water O (TIP3P/SPC); opls_111
                oxygen_type = {'opls': 'opls_111', 'gaff': 'ow'}

            elif (info['NumConnectedAtoms'] == 2 and 'H' in connected_symbols):
                # oh: alcohol or phenol O (R–OH); opls_154
                oxygen_type = {'opls': 'opls_154', 'gaff': 'oh'}

            elif info['NumConnectedAtoms'] == 2:
                # os: ether / ester O (sp3, 2-connected, no H); opls_180
                oxygen_type = {'opls': 'opls_180', 'gaff': 'os'}

            elif info['NumConnectedAtoms'] == 1:

                # If necessary we could check if the carbon connected
                # to the oxygen is connected to another oxygen. It is
                # useful to identify carboxylic acids and esters.

                # carbons = [
                #     atom for atom in info['ConnectedAtomsNumbers']
                #     if atom_info_dict[atom]['AtomicSymbol'] == 'C'
                # ]
                # if any('O' in atom_info_dict[carbon]['ConnectedAtoms']
                #        for carbon in carbons):
                #     ...

                if connected_symbols == {'H'}:
                    # oxygen in OH
                    oxygen_type = {'opls': 'opls_154', 'gaff': 'oh'}
                else:
                    oxygen_type = {'opls': 'opls_XXX', 'gaff': 'o'}

        # Cyclic

        elif info.get('CyclicStructure') == 'cycle':

            if (info['NumConnectedAtoms'] == 2 and connected_symbols == {'H'}):
                # oh: cyclic O–H (e.g. sugar hydroxyl); opls_154
                oxygen_type = {'opls': 'opls_154', 'gaff': 'oh'}

            elif info['NumConnectedAtoms'] == 2:
                # os: cyclic ether O (furan, THF, etc.); opls_180
                oxygen_type = {'opls': 'opls_180', 'gaff': 'os'}

            elif info['NumConnectedAtoms'] == 1:
                # o: cyclic carbonyl O (lactone, lactam); opls_236
                oxygen_type = {'opls': 'opls_236', 'gaff': 'o'}

        return oxygen_type

    def extract_gaff_atom_types(self, atom_types_dict):
        """
        Extracts GAFF atom types from the atom types dictionary.
        """

        # Initialize the list of gaff atom types
        gaff_atom_types = []

        # Append the gaff atom types to the list
        for atom_type in atom_types_dict.keys():
            if isinstance(atom_types_dict[atom_type], dict):
                gaff_type = atom_types_dict[atom_type].get('gaff', None)
                uff_type = atom_types_dict[atom_type].get('uff', None)
                if gaff_type:
                    gaff_atom_types.append(gaff_type)
                else:
                    gaff_atom_types.append(f'{uff_type}_unknown')
        return gaff_atom_types

    def check_alternating_atom_types(self):
        """
        Checks alternating atom types in GAFF, including cc, cd, ce, cf, cg,
        ch, nc, nd, ne and nf.
        """

        atom_types = list(self.gaff_atom_types)

        # Look for bonds formed between cc, ce, cg, nc, ne

        assigned_bonds = []
        counted_atom_ids = []

        while True:

            counted_atom_ids_prev = list(counted_atom_ids)

            # attempt to update counted_atom_ids and assigned_bonds by
            # following the chain of bonds (note that the ordering in
            # counted_atom_ids is important)

            for i in counted_atom_ids:
                for j in range(len(atom_types)):
                    if (j not in counted_atom_ids
                            and atom_types[j] in ['cc', 'ce', 'cg', 'nc', 'ne']
                            and self.connectivity_matrix[i][j] == 1):
                        assigned_bonds.append((i, j))
                        counted_atom_ids.append(j)

            # in case counted_atom_ids remains unchanged, look for a new chain
            # of bonds

            if counted_atom_ids_prev == counted_atom_ids:

                all_found = True
                new_index = None

                for i in range(len(atom_types)):
                    if (atom_types[i] in ['cc', 'ce', 'cg', 'nc', 'ne']
                            and i not in counted_atom_ids):
                        all_found = False
                        new_index = i
                        break

                if all_found:
                    break
                else:
                    counted_atom_ids.append(new_index)

        # Check distances for alternation

        conjugated_atom_type_pairs = {
            'cc': ('cc', 'cd'),
            'cd': ('cc', 'cd'),
            'ce': ('ce', 'cf'),
            'cf': ('ce', 'cf'),
            'cg': ('cg', 'ch'),
            'ch': ('cg', 'ch'),
            'nc': ('nc', 'nd'),
            'nd': ('nc', 'nd'),
            'ne': ('ne', 'nf'),
            'nf': ('ne', 'nf'),
        }

        for i, j in assigned_bonds:

            is_c_n_bond = ((atom_types[i][0] == 'c' and atom_types[j][0] == 'n')
                           or (atom_types[i][0] == 'n'
                               and atom_types[j][0] == 'c'))

            is_n_n_bond = (atom_types[i][0] == 'n' and atom_types[j][0] == 'n')

            has_sp1_carbon = (atom_types[i] in ['cg', 'ch']
                              or atom_types[j] in ['cg', 'ch'])

            if is_c_n_bond or is_n_n_bond or has_sp1_carbon:
                single_bond_thresh = 1.3475
            else:
                single_bond_thresh = 1.3985

            if self.distance_matrix[i][j] <= single_bond_thresh:
                if atom_types[i] in ['cc', 'ce', 'cg', 'nc', 'ne']:
                    atom_types[j] = conjugated_atom_type_pairs[atom_types[j]][1]
                else:
                    atom_types[j] = conjugated_atom_type_pairs[atom_types[j]][0]
            else:
                if atom_types[i] in ['cc', 'ce', 'cg', 'nc', 'ne']:
                    atom_types[j] = conjugated_atom_type_pairs[atom_types[j]][0]
                else:
                    atom_types[j] = conjugated_atom_type_pairs[atom_types[j]][1]

        # Check and correct e.g. ce-ce-ce-ce chain
        # Need to find the central pair of ce-ce-ce-ce chain
        # Then correct central pair if both atoms have two ce-ce/cf bonds

        while True:
            found_correction = False
            for i, j in assigned_bonds:
                if atom_types[i] not in ['cc', 'ce', 'cg', 'nc', 'ne']:
                    continue
                if atom_types[j] not in ['cc', 'ce', 'cg', 'nc', 'ne']:
                    continue
                i_count = {'ee_and_ef': 1, 'ee': 1}
                j_count = {'ee_and_ef': 1, 'ee': 1}
                for k, l in assigned_bonds:
                    if i == k and j != l:
                        i_count['ee_and_ef'] += 1
                        if atom_types[l] in ['cc', 'ce', 'cg', 'nc', 'ne']:
                            i_count['ee'] += 1
                    if i == l and j != k:
                        i_count['ee_and_ef'] += 1
                        if atom_types[k] in ['cc', 'ce', 'cg', 'nc', 'ne']:
                            i_count['ee'] += 1
                    if j == k and i != l:
                        j_count['ee_and_ef'] += 1
                        if atom_types[l] in ['cc', 'ce', 'cg', 'nc', 'ne']:
                            j_count['ee'] += 1
                    if j == l and i != k:
                        j_count['ee_and_ef'] += 1
                        if atom_types[k] in ['cc', 'ce', 'cg', 'nc', 'ne']:
                            j_count['ee'] += 1
                if (i_count['ee_and_ef'] == 2 and i_count['ee'] > 1
                        and j_count['ee_and_ef'] == 2 and j_count['ee'] > 1):
                    atom_types[i] = conjugated_atom_type_pairs[atom_types[i]][1]
                    atom_types[j] = conjugated_atom_type_pairs[atom_types[j]][1]
                    found_correction = True
                    break
            if not found_correction:
                break

        # Look for bonds formed between cp

        for i, at_i in enumerate(atom_types):
            for j, at_j in enumerate(atom_types):
                if (j > i and self.connectivity_matrix[i][j] == 1
                        and at_i in ['cp', 'cq'] and at_j in ['cp', 'cq']):
                    if self.get_common_cycles(i, j, 'pure_aromatic'):
                        atom_types[j] = 'cq' if atom_types[i] == 'cp' else 'cp'
                    else:
                        atom_types[j] = atom_types[i]

        self.gaff_atom_types = atom_types

        for atom_key, gaff_type in zip(self.atom_types_dict,
                                       self.gaff_atom_types):
            self.atom_types_dict[atom_key]['gaff'] = gaff_type

    def get_common_cycles(self, i, j, cycle_type='any'):
        """
        Gets number of common cycles shared by a pair of atoms.

        :param i:
            The index of the first atom.
        :param j:
            The index of the second atom.
        :param cycle_type:
            The type of cycles.

        :return:
            The number of common cycles.
        """

        common_cycle_numbers = (
            set(self.atom_info_dict[i + 1].get('CycleNumber'))
            & set(self.atom_info_dict[j + 1].get('CycleNumber')))

        if cycle_type in ['pure_aromatic', 'non_pure_aromatic', 'non_aromatic']:
            return [
                cycle_num for cycle_num in common_cycle_numbers
                if self.aromaticity[cycle_num] == cycle_type
            ]
        else:
            return list(common_cycle_numbers)

    def generate_gaff_atomtypes(self, molecule, connectivity_matrix=None):
        """
        Generates GAFF (General Amber Force Field) atom types for a given molecule.

        :param molecule:
            A VeloxChem molecule object.

        :return:
            A list of GAFF atom types for each atom in the molecule.
        """

        # Workflow of the method
        self.coordinates = molecule.get_coordinates_in_angstrom()
        self.atomic_symbols = molecule.get_labels()

        if connectivity_matrix is None:
            self.connectivity_matrix = molecule.get_connectivity_matrix()
        else:
            # For the add bond feature
            self.connectivity_matrix = connectivity_matrix.copy()

        self.distance_matrix = molecule.get_distance_matrix_in_angstrom()

        self.detect_closed_cyclic_structures()
        self.atom_info_dict = self.create_atom_info_dict()
        self.atom_types_dict = self.decide_atom_type(self.atom_info_dict)

        self.gaff_atom_types = self.extract_gaff_atom_types(
            self.atom_types_dict)
        self.check_alternating_atom_types()

        # Printing output
        self.ostream.print_info("VeloxChem Atom Type Identification")
        self.ostream.print_info("-" * 40)  # Dashed line

        # Detected number of atoms
        num_atoms = molecule.number_of_atoms()
        self.ostream.print_info(f"Detected number of atoms: {num_atoms}")

        # Print table header
        self.ostream.print_info("{:<30} {:<20}".format(
            "Symbol (id)", "GAFF atom type assigned"))

        # Print atom symbol, atom number, and GAFF atom type for each atom
        for i, (symbol, gaff_type) in enumerate(zip(self.atomic_symbols,
                                                    self.gaff_atom_types),
                                                start=1):
            self.ostream.print_info("{:<30} {:<20}".format(
                f"{symbol} ({i})", gaff_type))

        # Print cycle information (aromaticity) for each cycle
        cycle_sizes = [len(cycle) for cycle in self.reduced_cycles]

        if cycle_sizes:
            cycle_suffix = 's' if len(cycle_sizes) > 1 else ''
            self.ostream.print_blank()
            self.ostream.print_info(
                f"Detected {len(cycle_sizes)} cycle{cycle_suffix}:")

        for size, aromaticity in zip(cycle_sizes, self.aromaticity):
            if aromaticity == "non_aromatic":
                self.ostream.print_info(
                    f"Cycle size {size}: Non-aromatic Cycle")
            elif aromaticity == "non_pure_aromatic":
                self.ostream.print_info(
                    f"Cycle size {size}: Non-pure Aromatic Cycle")
            elif aromaticity == "pure_aromatic":
                self.ostream.print_info(
                    f"Cycle size {size}: Pure Aromatic Cycle")

        self.ostream.flush()

        return self.gaff_atom_types

    def identify_equivalences(self, depth=10):
        """
        Identifies equivalent atoms in the molecule.
        The depth parameter specifies how many bonds are considered for the equivalence.
        The default value is 10.

        :param depth:
            Depth of the equivalence search
        """

        def gather_neighbors(atom_index, current_depth=0, path=()):
            """
            Gather the paths to neighbors of an atom up to a certain depth.
            """

            if current_depth == depth:
                return []

            if current_depth == 0 and not path:
                path = (atom_index, )

            neighbors = []

            for i, connected in enumerate(connectivity_matrix[atom_index]):
                if connected and i not in path:
                    new_path = path + (i, )
                    neighbors.append(new_path)
                    if current_depth < depth - 1:
                        neighbors.extend(
                            gather_neighbors(i, current_depth + 1, new_path))

            return neighbors

        # Main logic for identifying equivalences

        conjugated_atomtype_mapping = {
            'cc': 'cd',
            'cd': 'cc',
            'ce': 'cf',
            'cf': 'ce',
            'cg': 'ch',
            'ch': 'cg',
            'nc': 'nd',
            'nd': 'nc',
            'ne': 'nf',
            'nf': 'ne',
            'cp': 'cq',
            'cq': 'cp',
        }

        conjugated_cc_cd_mapping = {
            'cc': 'cd',
            'cd': 'cc',
        }

        conjugated_nc_nd_mapping = {
            'nc': 'nd',
            'nd': 'nc',
        }

        self.equivalent_atoms = []
        atom_types_for_equil = []
        for name, type in self.atom_types_dict.items():
            if 'gaff' in type:
                self.equivalent_atoms.append(f'{type["gaff"]}_00')
                atom_types_for_equil.append(type['gaff'])
            else:
                self.equivalent_atoms.append(f'{name}_00')
                atom_types_for_equil.append(name)

        connectivity_matrix = self.connectivity_matrix

        for atom_type in list(set(atom_types_for_equil)):

            # skip cd/cf/ch/nd/nf/cq since they will be counted by
            # cc/ce/cg/nc/ne/cp
            if atom_type == 'cd' and 'cc' in atom_types_for_equil:
                continue
            if atom_type == 'cf' and 'ce' in atom_types_for_equil:
                continue
            if atom_type == 'ch' and 'cg' in atom_types_for_equil:
                continue
            if atom_type == 'nd' and 'nc' in atom_types_for_equil:
                continue
            if atom_type == 'nf' and 'ne' in atom_types_for_equil:
                continue
            if atom_type == 'cq' and 'cp' in atom_types_for_equil:
                continue

            swapped_atom_type = conjugated_atomtype_mapping.get(
                atom_type, atom_type)

            atom_type_count = atom_types_for_equil.count(atom_type)
            if swapped_atom_type != atom_type:
                atom_type_count += atom_types_for_equil.count(swapped_atom_type)

            if atom_type_count == 1:
                continue

            atom_paths = defaultdict(set)

            for idx, at in enumerate(atom_types_for_equil):
                if at in [atom_type, swapped_atom_type]:
                    paths = gather_neighbors(idx)
                    path_types = [
                        tuple(atom_types_for_equil[step] for step in path)
                        for path in paths
                    ]
                    path_types = tuple(sorted(path_types))

                    # this is full swapping using conjugated_atomtype_mapping
                    swapped_path_types = [
                        tuple(
                            conjugated_atomtype_mapping.get(at, at)
                            for at in path) for path in path_types
                    ]
                    swapped_path_types = tuple(sorted(swapped_path_types))

                    # this is just cc-cd swapping using conjugated_cc_cd_mapping
                    swapped_cc_cd_path_types = [
                        tuple(
                            conjugated_cc_cd_mapping.get(at, at) for at in path)
                        for path in path_types
                    ]
                    swapped_cc_cd_path_types = tuple(
                        sorted(swapped_cc_cd_path_types))

                    # this is just nc-nd swapping using conjugated_nc_nd_mapping
                    swapped_nc_nd_path_types = [
                        tuple(
                            conjugated_nc_nd_mapping.get(at, at) for at in path)
                        for path in path_types
                    ]
                    swapped_nc_nd_path_types = tuple(
                        sorted(swapped_nc_nd_path_types))

                    # we should include more single/double/... pair swapping
                    # but for now doing only cc-cd and nc-nd seems to suffice

                    if swapped_path_types in atom_paths:
                        atom_paths[swapped_path_types].add(idx)
                    elif swapped_cc_cd_path_types in atom_paths:
                        atom_paths[swapped_cc_cd_path_types].add(idx)
                    elif swapped_nc_nd_path_types in atom_paths:
                        atom_paths[swapped_nc_nd_path_types].add(idx)
                    else:
                        atom_paths[path_types].add(idx)

            for eq_id, (path_set, indices) in enumerate(atom_paths.items()):
                for idx in indices:
                    self.equivalent_atoms[idx] = f'{atom_type}_{eq_id:02d}'

        equal_charges_str_list = []

        for eq_at in list(set(self.equivalent_atoms)):
            if self.equivalent_atoms.count(eq_at) > 1:
                eq_str_list = [
                    str(i + 1) for i, a in enumerate(self.equivalent_atoms)
                    if a == eq_at
                ]
                eq_str = ' = '.join(eq_str_list)
                equal_charges_str_list.append(eq_str)

        self.equivalent_charges = ', '.join(equal_charges_str_list)
