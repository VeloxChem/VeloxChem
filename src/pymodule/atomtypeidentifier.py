#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
import re
import networkx as nx
import math

from .molecule import Molecule
from .veloxchemlib import bohr_in_angstrom
from .veloxchemlib import mathconst_pi


class AtomTypeIdentifier:
    """
    A class to identify atom types in a molecule using GAFF (General Amber Force Field) atom types based on 
    a VeloxChem molecule object.

    The class processes a molecule object containing the atomic coordinates of a molecule to determine the types
    of atoms according to the GAFF. It involves several steps including reading the file, determining covalent 
    radii, creating a connectivity matrix, identifying cyclic structures, and assigning atom types.

    Instance variables
        - molecule: A VeloxChem molecule object.
        - atomic_symbols: A list of atomic symbols for each atom in the molecule.
        - coordinates: A 2D numpy array of atomic coordinates for each atom in the molecule.
        - covalent_radii: A list of covalent radii for each atom in the molecule.
        - connectivity_matrix: A 2D numpy array indicating which atoms are bonded.
        - distance_matrix: A 2D numpy array containing the distances between atoms.
        - cyclic_atoms: A set of atom indices that are part of a cyclic structure.
        - cycle_sizes: A list of the sizes of each cyclic structure.
        - aromaticity: A list of aromaticity classifications for each cyclic structure.
        - cycles: A list of cyclic structures.
        - cycle_ids: A list of cyclic structure IDs.
        - atom_cycle_info: A dictionary containing cycle information for each atom.
        - atom_info_dict: A dictionary containing detailed information for each atom.
        - atom_types_dict: A dictionary containing the GAFF atom types for each atom.

    """

    def __init__(self):
        """
        Initializes the AtomTypeIdentifier instance.

        Args:
            self
        """

    def create_connectivity_matrix(self, factor=1.3):
        """
        Creates a connectivity matrix for the molecule based on the atomic coordinates
        and covalent radii, determining which atoms are bonded.

        This method iterates through pairs of atoms and calculates the distance between
        them. If this distance is less than or equal to the sum of their covalent radii
        scaled by a factor, it is considered a bond, and the connectivity matrix is updated
        accordingly. The method also constructs a corresponding distance matrix with the
        actual distances between connected atoms.

        :param factor: 
            A scaling factor for the covalent radii to account for the bond threshold.
            Default value is 1.3.
        :return:
            tuple: A tuple containing two 2D numpy arrays:
                   - The first array is the connectivity matrix with 1s indicating bonded atom pairs.
                   - The second array is the distance matrix with actual distances between atoms.
        """
        num_atoms = len(self.coordinates)
        self.connectivity_matrix = np.zeros((num_atoms, num_atoms), dtype=int)
        self.distance_matrix = np.zeros((num_atoms, num_atoms))

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                distance = self.measure_length(self.coordinates[i],
                                               self.coordinates[j])
                covalent_distance = self.covalent_radii[
                    i] + self.covalent_radii[j]
                adjusted_threshold = covalent_distance * float(factor)

                if distance <= adjusted_threshold:
                    self.connectivity_matrix[i][j] = 1
                    self.connectivity_matrix[j][i] = 1
                    self.distance_matrix[i][j] = distance
                    self.distance_matrix[j][i] = distance

        return self.connectivity_matrix, self.distance_matrix

    def plot_connectivity_map_3D(self):
        """
        Plots a 3D connectivity map of the molecule using Matplotlib.

        """
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
                    symbol,
                    fontsize=12,
                    ha='center',
                    va='center')

        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Connectivity Map in 3D')
        plt.show()

        pass

    def is_sp2_carbon(self, atom_idx):
        """
        Determines if a C atom, identified by its index, is sp2 hybridized.

        :param atom_idx:
            Index of the atom in the molecule.

        :return:
            True if the atom is sp2 hybridized, False otherwise.
        """
        return self.atomic_symbols[atom_idx] == "C" and len(
            list(self.graph.neighbors(atom_idx))) == 3

    def is_sp2_nitrogen(self, atom_idx):
        """
        Determines if a N atom, identified by its index, is sp2 hybridized.

        :param atom_idx:
            Index of the atom in the molecule.

        :return:
            True if the atom is sp2 hybridized, False otherwise.
        """
        return self.atomic_symbols[atom_idx] == "N" and len(
            list(self.graph.neighbors(atom_idx))) == 2

    def detect_closed_cyclic_structures(self):
        """
        Detects closed cyclic structures in a molecule and determines their aromaticity.

        This method analyzes the graph of atoms and their connectivity to identify cycles,
        determine the size of each cycle, and classify them based on aromaticity criteria.

        """
        self.graph = nx.Graph()
        for i in range(len(self.atomic_symbols)):
            self.graph.add_node(i)
            for j in range(i + 1, len(self.atomic_symbols)):
                if self.connectivity_matrix[i][j] == 1:
                    self.graph.add_edge(i, j)

        all_cycles = list(nx.simple_cycles(self.graph))

        # Sort by cycle length
        all_cycles = sorted(all_cycles, key=len)

        # Filter only the cycles of size 3 to 6
        filtered_cycles = [
            cycle for cycle in all_cycles if 3 <= len(cycle) <= 6
        ]

        # Remove super-cycles (cycles that contain smaller cycles)
        self.reduced_cycles = filtered_cycles[:]
        cycles_to_remove = set()
        for i, cycle in enumerate(filtered_cycles):
            for j, larger_cycle in enumerate(filtered_cycles):
                if len(cycle) < len(larger_cycle) and set(cycle).issubset(
                        set(larger_cycle)):
                    cycles_to_remove.add(tuple(larger_cycle))

        # Convert cycles_to_remove to list of lists and then to a set of tuples for faster lookup
        cycles_to_remove = {tuple(c) for c in cycles_to_remove}

        # Rebuild reduced_cycles excluding the ones in cycles_to_remove
        self.reduced_cycles = [
            cycle for cycle in self.reduced_cycles
            if tuple(cycle) not in cycles_to_remove
        ]

        self.cyclic_atoms = set()
        self.cycle_sizes = []
        self.aromaticity = []
        self.cycle_ids = []
        self.atom_cycle_info = {}

        # Assignation of aromaticity to all the reduced cycles

        for cycle_id, cycle in enumerate(self.reduced_cycles):
            self.cyclic_atoms.update(cycle)

            size = len(cycle)
            self.cycle_sizes.append(size)

            cycle_elements = [self.atomic_symbols[i] for i in cycle]

            cc_bond_distances = [
                self.distance_matrix[cycle[i]][cycle[(i + 1) % size]]
                for i in range(size)
                if self.atomic_symbols[cycle[i]] == "C" and
                self.atomic_symbols[cycle[(i + 1) % size]] == "C"
            ]

            max_distance = max(cc_bond_distances) if cc_bond_distances else 0
            min_distance = min(cc_bond_distances) if cc_bond_distances else 0

            if size == 6 and all(elem in ["C", "N"] for elem in cycle_elements):
                if (all(
                        self.is_sp2_carbon(atom_idx)
                        for atom_idx in cycle
                        if self.atomic_symbols[atom_idx] == "C") and all(
                            self.is_sp2_nitrogen(atom_idx)
                            for atom_idx in cycle
                            if self.atomic_symbols[atom_idx] == "N")):
                    aro = "pure_aromatic" if max_distance - min_distance <= 0.08 else "non_pure_aromatic"
                elif (all(
                        self.is_sp2_carbon(atom_idx)
                        for atom_idx in cycle
                        if self.atomic_symbols[atom_idx] == "C") and
                      max_distance - min_distance <= 0.08):
                    aro = "non_pure_aromatic"
                else:
                    aro = "non_aromatic"
            elif (size == 5 and all(
                    self.is_sp2_carbon(atom_idx)
                    for atom_idx in cycle
                    if self.atomic_symbols[atom_idx] == "C")):

                if "S" in cycle_elements:
                    # Check if the S is connected to 2 atoms in the cycle, if so, it is non_pure_aromatic
                    if len(
                            list(
                                self.graph.neighbors(
                                    cycle[cycle_elements.index("S")]))) == 2:
                        aro = "non_pure_aromatic"
                    else:
                        aro = "non_aromatic"
                else:
                    if 'N' in cycle_elements or 'O' in cycle_elements:
                        aro = "non_pure_aromatic"
                    else:
                        aro = "non_aromatic"

            elif size == 4:
                if all(
                        self.is_sp2_carbon(atom_idx)
                        for atom_idx in cycle
                        if self.atomic_symbols[atom_idx] == "C"):
                    aro = "non_pure_aromatic"
                else:
                    aro = "non_aromatic"
            else:
                aro = "non_aromatic"

            self.aromaticity.append(aro)

            for atom in cycle:
                if atom not in self.atom_cycle_info:
                    self.atom_cycle_info[atom] = {
                        'sizes': [],
                        'aromaticities': []
                    }
                self.atom_cycle_info[atom]['sizes'].append(size)
                self.atom_cycle_info[atom]['aromaticities'].append(aro)

            self.cycle_ids.append(cycle_id)

        # Additional logic for reassignment of aromaticity in special cases where 3 atoms are shared with aromatic rings.
        for index, cycle in enumerate(self.reduced_cycles):
            # Check if all carbons in the cycle are sp2
            all_carbons_sp2 = all(
                self.is_sp2_carbon(atom_idx)
                for atom_idx in cycle
                if self.atomic_symbols[atom_idx] == "C")
            if self.cycle_sizes[index] == 5 and self.aromaticity[
                    index] == 'non_aromatic' and all_carbons_sp2:
                count_pure_aromatic_atoms = sum(
                    1 for atom in cycle if 'pure_aromatic' in
                    self.atom_cycle_info[atom]['aromaticities'])
                if count_pure_aromatic_atoms >= 3:
                    self.aromaticity[index] = 'non_pure_aromatic'
                    for atom in cycle:
                        self.atom_cycle_info[atom]['aromaticities'] = [
                            'non_pure_aromatic' if a == 'non_aromatic' else a
                            for a in self.atom_cycle_info[atom]['aromaticities']
                        ]

    def create_atom_info_dict(self):
        """
        Creates a dictionary containing detailed information for each atom in the molecule.

        This method compiles the atomic symbol, atom number, number of connected atoms, symbols of connected atoms,
        atom numbers of connected atoms, distances to connected atoms, and cycle information into a structured dictionary.

        :return:
            The dictionary where each key is an atom number and each value is another dictionary of atom information.
        """
        self.atom_info_dict = {}
        for i, symbol in enumerate(self.atomic_symbols):
            num_connected_atoms = np.sum(self.connectivity_matrix[i])
            connected_atoms_numbers = [
                j + 1
                for j in range(len(self.atomic_symbols))
                if self.connectivity_matrix[i][j] == 1
            ]
            connected_atoms_symbols = [
                self.atomic_symbols[j]
                for j in range(len(self.atomic_symbols))
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
                    self.cycle_ids[self.reduced_cycles.index(c)]
                    for c in self.reduced_cycles
                    if i in c
                ]
            else:
                info["CyclicStructure"] = "none"

            self.atom_info_dict[i + 1] = info

    def decide_atom_type(self):
        """
        Analyzes the molecular structure information to assign atom types to each atom in the molecule.

        This method traverses through the atom information dictionary created from the molecular structure data
        and applies a series of rules to determine the appropriate atom type for each atom. The rules consider
        factors such as the atom's chemical environment, its connectivity to other atoms, and whether it is part of
        a cyclic structure.

        :return:
            A dictionary where each key corresponds to an atom identifier (e.g., "C1" for the first carbon atom),
            and each value is another dictionary containing the 'opls' and 'gaff' force field identifiers for the atom.
        """

        self.bad_hydrogen = False

        self.atom_types_dict = {}

        for atom_number, info in self.atom_info_dict.items():
            
            if info['AtomicSymbol'] == 'C':

                # If this carbon was previously assigned, skip the rest of the loop
                key = f"C{info['AtomNumber']}"
                if key in self.atom_types_dict:
                    continue

                # Chemical environment information
                connected_symbols = set(info['ConnectedAtoms'])
                connected_distances = info['ConnectedAtomsDistances']

                # Cyclic and pure aromatic
                if info.get('CyclicStructure'
                           ) == 'cycle' and 'pure_aromatic' in info.get(
                               'Aromaticity'):

                    if info['NumConnectedAtoms'] == 3:
                        # Check for identifying biphenyls
                        # TODO: Clean up this code / Check 0091.xyz. Use normal for loop
                        connected_carbons_in_diff_cycle_and_pure_aromatic = []

                        for connected_atom_number in info[
                                'ConnectedAtomsNumbers']:
                            connected_atom_info = self.atom_info_dict[
                                connected_atom_number]

                            if connected_atom_info.get('AtomicSymbol') == 'C' and connected_atom_info.get('CycleNumber') and \
                                not set(connected_atom_info.get('CycleNumber')) & set(info.get('CycleNumber')) and \
                                'pure_aromatic' in connected_atom_info.get('Aromaticity'):

                                connected_carbons_in_diff_cycle_and_pure_aromatic.append(
                                    connected_atom_number)

                        # If the list is not empty, set the connected_carbon_atom to the first element. Else, set it to None.
                        connected_carbon_atom = connected_carbons_in_diff_cycle_and_pure_aromatic[
                            0] if connected_carbons_in_diff_cycle_and_pure_aromatic else None

                        if connected_symbols == {'C', 'H'}:
                            carbon_type = {'opls': 'opls_145', 'gaff': 'ca'}
                        elif connected_symbols == {'C', 'N', 'H'}:
                            carbon_type = {'opls': 'opls_521', 'gaff': 'ca'}

                        elif connected_carbon_atom:
                            carbon_type = {'opls': 'opls_521', 'gaff': 'cp'}

                            index_in_distances = info[
                                'ConnectedAtomsNumbers'].index(
                                    connected_carbon_atom)
                            d = info['ConnectedAtomsDistances'][
                                index_in_distances]

                            if d >= 1.4685:
                                biphenyl_carbon = {
                                    'opls': 'opls_521',
                                    'gaff': 'cp'
                                }
                            elif d <= 1.4685:
                                biphenyl_carbon = {
                                    'opls': 'opls_CQ',
                                    'gaff': 'cq'
                                }

                            # Store the atomtype in the dictionary for the biphenyl carbon
                            self.atom_types_dict[
                                f"C{connected_carbon_atom}"] = biphenyl_carbon
                        else:
                            carbon_type = {'opls': 'opls_145', 'gaff': 'ca'}

                    elif info[
                            'NumConnectedAtoms'] == 4 and connected_symbols == {
                                'C', 'H'
                            }:
                        carbon_type = {'opls': 'opls_135', 'gaff': 'c3'}
                    else:
                        carbon_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'cx{info["AtomNumber"]}'
                        }

                # Default assignation for carbons in non-pure aromatic cyclic structures
                elif info.get('CyclicStructure'
                             ) == 'cycle' and 'non_pure_aromatic' in info.get(
                                 'Aromaticity'):

                    if 'O' in connected_symbols:
                        # Directly loop through the connected atom numbers
                        for connected_atom_number in info[
                                'ConnectedAtomsNumbers']:
                            if self.atom_info_dict[connected_atom_number][
                                    'AtomicSymbol'] == 'O' and self.atom_info_dict[
                                        connected_atom_number][
                                            'CyclicStructure'] == 'none':
                                if self.atom_info_dict[connected_atom_number][
                                        'NumConnectedAtoms'] == 1:
                                    carbon_type = {
                                        'opls': 'opls_235',
                                        'gaff': 'c'
                                    }  # Carbonyl Carbon
                                    break  # Exit the loop once the carbonyl carbon is identified
                            else:
                                carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}

                    elif 'C' in connected_symbols and info[
                            'NumConnectedAtoms'] == 3 and any(
                                size in info['CycleSize']
                                for size in {4, 5, 6}):

                        carbon_type = {
                            'opls': 'opls_508',
                            'gaff': 'cc'
                        }  # Default assignment

                        # TODO: This part of the code will be deprecated in the future.
                        # The method check_alternating_atoms() will be used instead.

                        # Immediately assign C atom types based on distances
                        connected_atoms_numbers = info['ConnectedAtomsNumbers']

                        # Current atom's cycle number (will be 'none' if not in a cycle)
                        current_cycle_number = info.get('CycleNumber')

                        for neighboring_atom_number, distance in zip(
                                connected_atoms_numbers, connected_distances):
                            # Check if the neighbour is also a C
                            if self.atom_info_dict[neighboring_atom_number][
                                    'AtomicSymbol'] != 'C':
                                continue
                            key_neigh = f"C{neighboring_atom_number}"
                            # Check if the C atom has been previously assigned
                            if key_neigh in self.atom_types_dict:
                                existing_type = self.atom_types_dict[key_neigh]

                                # For neighboring atom type {'opls': 'opls_XXX', 'gaff': 'cc'}
                                if existing_type == {
                                        'opls': 'opls_508',
                                        'gaff': 'cc'
                                }:
                                    if distance <= 1.4:
                                        carbon_type = {
                                            'opls': 'opls_XXX',
                                            'gaff': 'cd'
                                        }
                                    elif distance > 1.4:
                                        carbon_type = {
                                            'opls': 'opls_508',
                                            'gaff': 'cc'
                                        }

                                # For neighboring atom type {'opls': 'opls_XXX', 'gaff': 'cd'}
                                elif existing_type == {
                                        'opls': 'opls_XXX',
                                        'gaff': 'cd'
                                }:
                                    if distance <= 1.4:
                                        carbon_type = {
                                            'opls': 'opls_508',
                                            'gaff': 'cc'
                                        }
                                    elif distance > 1.4:
                                        carbon_type = {
                                            'opls': 'opls_XXX',
                                            'gaff': 'cd'
                                        }

                                else:
                                    continue
                            else:

                                neighboring_cycle_number = self.atom_info_dict[
                                    neighboring_atom_number].get(
                                        'CycleNumber', 'none')

                                # Check if neighboring atom is not in a cycle or belongs to another cycle ID
                                if neighboring_cycle_number == 'none' or neighboring_cycle_number != current_cycle_number:
                                    continue
                                else:
                                    if distance <= 1.4:  # Double bonded
                                        conn_carbon_type = {
                                            'opls': 'opls_XXX',
                                            'gaff': 'cd'
                                        }

                                    elif distance >= 1.4:  # Single bonded
                                        conn_carbon_type = {
                                            'opls': 'opls_508',
                                            'gaff': 'cc'
                                        }

                                    # Store the atomtype in the dictionary
                                    self.atom_types_dict[
                                        f"C{neighboring_atom_number}"] = conn_carbon_type

                                    # If there is an H connected to this carbon, assign it now
                                    for connected_atom_number in self.atom_info_dict[
                                            neighboring_atom_number][
                                                'ConnectedAtomsNumbers']:
                                        connected_atom_info = self.atom_info_dict[
                                            connected_atom_number]
                                        if connected_atom_info[
                                                'AtomicSymbol'] == 'H' and connected_atom_info[
                                                    'NumConnectedAtoms'] == 1:
                                            ewd_atoms = [
                                                'N', 'Br', 'Cl', 'I', 'F', 'S',
                                                'O'
                                            ]
                                            ewd_count = sum(
                                                1
                                                for num in self.atom_info_dict[
                                                    neighboring_atom_number]
                                                ['ConnectedAtomsNumbers']
                                                if self.atom_info_dict[num]
                                                ['AtomicSymbol'] in ewd_atoms)
                                            if ewd_count == 1:
                                                hydrogen_type = {
                                                    'opls': 'opls_146',
                                                    'gaff': 'h4'
                                                }  # 1 EWD atom
                                            elif ewd_count == 2:
                                                hydrogen_type = {
                                                    'opls': 'opls_h5',
                                                    'gaff': 'h5'
                                                }  # 2 EWD atom
                                            else:
                                                hydrogen_type = {
                                                    'opls': 'opls_146',
                                                    'gaff': 'ha'
                                                }  # Default Aliphatic sp2 Hydrogen for c2

                                            # Store it in the dictionary
                                            self.atom_types_dict[
                                                f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

                                    # Next iteration
                                    continue
                    else:
                        carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}

                # Aliphatic cycles
                elif info.get('CyclicStructure'
                             ) == 'cycle' and 'non_aromatic' in info.get(
                                 'Aromaticity'):

                    if info['NumConnectedAtoms'] == 4 and 3 in info['CycleSize']:
                        carbon_type = {'opls': 'opls_CX', 'gaff': 'cx'}
                    elif info['NumConnectedAtoms'] == 4 and 4 in info[
                            'CycleSize']:
                        carbon_type = {'opls': 'opls_CY', 'gaff': 'cy'}
                    elif info['NumConnectedAtoms'] == 4 and 5 in info[
                            'CycleSize']:
                        carbon_type = {'opls': 'opls_c5', 'gaff': 'c5'}
                    elif info['NumConnectedAtoms'] == 4 and 6 in info[
                            'CycleSize']:
                        carbon_type = {'opls': 'opls_c6', 'gaff': 'c6'}
                    elif info['NumConnectedAtoms'] == 3 and 3 in info[
                            'CycleSize']:
                        carbon_type = {'opls': 'opls_CU', 'gaff': 'cu'}
                    elif info['NumConnectedAtoms'] == 3 and 4 in info[
                            'CycleSize']:
                        carbon_type = {'opls': 'opls_CV', 'gaff': 'cv'}

                    # Cases for general Non-Aromatic cycles bigger than 4
                    elif info['NumConnectedAtoms'] == 3:
                        if 'O' in connected_symbols:
                            # Directly loop through the connected atom numbers
                            for connected_atom_number in info[
                                    'ConnectedAtomsNumbers']:
                                if self.atom_info_dict[connected_atom_number][
                                        'AtomicSymbol'] == 'O':
                                    if self.atom_info_dict[
                                            connected_atom_number][
                                                'NumConnectedAtoms'] == 1:
                                        carbon_type = {
                                            'opls': 'opls_235',
                                            'gaff': 'c'
                                        }  # Carbonyl Carbon
                                        break  # Exit the loop once the carbonyl carbon is identified
                                    else:
                                        carbon_type = {
                                            'opls': 'opls_141',
                                            'gaff': 'c2'
                                        }  # Alcohol

                        elif 'C' in connected_symbols:
                            # Count the number of sp2-hybridized carbons connected to the current carbon
                            sp2_carbon_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'C' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 3)
                            sp1_carbon_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'C' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 2)
                            if sp2_carbon_count + sp1_carbon_count == 2 or sp2_carbon_count + sp1_carbon_count == 3:  # If the current carbon is connected to 2 sp2 carbons
                                carbon_type = {
                                    'opls': 'opls_XXX',
                                    'gaff': 'ce'
                                }  # Inner Sp2 carbons in conjugated systems
                            else:
                                carbon_type = {
                                    'opls': 'opls_141',
                                    'gaff': 'c2'
                                }  # Aliphatic sp2 Carbon
                        else:
                            carbon_type = {
                                'opls': 'opls_141',
                                'gaff': 'c2'
                            }  #Generic sp2 C

                    else:
                        carbon_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'cx{info["AtomNumber"]}'
                        }

                elif info.get('CyclicStructure') == 'none':
                    if info['NumConnectedAtoms'] == 4:
                        if 'C' and 'H' in connected_symbols:
                            carbon_type = {
                                'opls': 'opls_135',
                                'gaff': 'c3'
                            }  # Aliphatic sp3 Carbon
                        elif 'O' in connected_symbols and connected_symbols != {
                                'C', 'O', 'O', 'O'
                        }:
                            carbon_type = {
                                'opls': 'opls_135',
                                'gaff': 'c3'
                            }  # Carbon bonded to Oxygen (but not in carboxylate)
                        else:
                            carbon_type = {
                                'opls': 'opls_135',
                                'gaff': 'c3'
                            }  # General case, anything connected to it

                    elif info['NumConnectedAtoms'] == 3:

                        if 'O' in connected_symbols:
                            # Directly loop through the connected atom numbers
                            for connected_atom_number in info[
                                    'ConnectedAtomsNumbers']:
                                if self.atom_info_dict[connected_atom_number][
                                        'AtomicSymbol'] == 'O':
                                    if self.atom_info_dict[
                                            connected_atom_number][
                                                'NumConnectedAtoms'] == 1:
                                        carbon_type = {
                                            'opls': 'opls_235',
                                            'gaff': 'c'
                                        }  # Carbonyl Carbon
                                        break  # Exit the loop once the carbonyl carbon is identified
                                    else:
                                        carbon_type = {
                                            'opls': 'opls_141',
                                            'gaff': 'c2'
                                        }

                        elif 'S' in connected_symbols:  # Double bond with a Sulfur
                            for connected_atom_number in info[
                                    'ConnectedAtomsNumbers']:
                                if self.atom_info_dict[connected_atom_number][
                                        'AtomicSymbol'] == 'S':
                                    if self.atom_info_dict[
                                            connected_atom_number][
                                                'NumConnectedAtoms'] == 1:
                                        carbon_type = {
                                            'opls': 'opls_cs',
                                            'gaff': 'cs'
                                        }  # Carbon double bonded to Sulfur
                                        break  # Exit the loop as soon as we find a sulfur atom with only one connected atom
                                    else:
                                        carbon_type = {
                                            'opls': 'opls_141',
                                            'gaff': 'c2'
                                        }

                        elif 'C' in connected_symbols:
                            # Count the number of sp2-hybridized carbons connected to the current carbon
                            sp2_carbon_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'C' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 3)
                            sp1_carbon_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'C' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 2)
                            n_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'N' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 1)
                            n2_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'N' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 2)
                            if sp2_carbon_count + sp1_carbon_count == 2 or sp2_carbon_count + sp1_carbon_count == 3:  
                                carbon_type = {
                                    'opls': 'opls_XXX',
                                    'gaff': 'ce'
                                }  # Inner Sp2 carbons in conjugated systems
                            elif sp2_carbon_count + sp1_carbon_count + n2_count == 2:  
                                carbon_type = {
                                    'opls': 'opls_XXX',
                                    'gaff': 'ce'
                                }  # Inner Sp2 carbons in conjugated systems
                            else:
                                carbon_type = {
                                    'opls': 'opls_141',
                                    'gaff': 'c2'
                                }  # Aliphatic sp2 Carbon

                        else:
                            carbon_type = {
                                'opls': 'opls_141',
                                'gaff': 'c2'
                            }  #Generic sp2 C

                    elif info['NumConnectedAtoms'] == 2:
                        if 'O' in connected_symbols:
                            carbon_type = {
                                'opls': 'opls_235',
                                'gaff': 'c1'
                            }  # Carbon in carbonyl group or acid anhydride
                        else:
                            # Count the number of sp2-hybridized carbons connected to the current carbon
                            sp2_carbon_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'C' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 3)
                            sp1_carbon_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'C' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 2)
                            n_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'N' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 1)
                            n2_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] ==
                                'N' and self.atom_info_dict[num]
                                ['NumConnectedAtoms'] == 2)
                            if sp2_carbon_count == 2:
                                carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}
                            # If the current carbon is connected to 2 sp2 carbons
                            elif sp2_carbon_count + sp1_carbon_count + n_count == 2:
                                carbon_type = {'opls': 'opls_XXX', 'gaff': 'cg'}
                            elif sp2_carbon_count + sp1_carbon_count + n2_count == 2:
                                carbon_type = {'opls': 'opls_XXX', 'gaff': 'ch'}
                            else:
                                carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}
                    elif info['NumConnectedAtoms'] == 1:
                        carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}
                    else:
                        carbon_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'cx{info["AtomNumber"]}'
                        }
                else:
                    carbon_type = {
                        'opls': f'opls_x{info["AtomNumber"]}',
                        'gaff': f'cx{info["AtomNumber"]}'
                    }

                # Assignment for the Hydrogens linked to the carbons
                # ewd = Electron withdrawing atoms
                ewd_atoms = ['N', 'Br', 'Cl', 'I', 'F', 'S', 'O']

                # TODO: Simplify the rules that can be in the same types
                # Check also remove long comments on the right side
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    connected_atom_info = self.atom_info_dict[
                        connected_atom_number]
                    if connected_atom_info[
                            'AtomicSymbol'] == 'H' and connected_atom_info[
                                'NumConnectedAtoms'] == 1:
                        # Sp1 Carbon
                        if carbon_type == {'opls': 'opls_235', 'gaff': 'c1'}:
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}
                        # Carbonyl
                        elif carbon_type == {'opls': 'opls_235', 'gaff': 'c'}:
                            # Count the number of connected 'EWD' atoms for this carbon
                            ewd_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] in
                                ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {
                                    'opls': 'opls_146',
                                    'gaff': 'h4'
                                }  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {
                                    'opls': 'opls_xxx',
                                    'gaff': 'h5'
                                }  # 2 EWD atoms
                        # Sp3 Carbon types
                        elif carbon_type == {
                                'opls': 'opls_135',
                                'gaff': 'c3'
                        } or carbon_type == {
                                'opls': 'opls_c5',
                                'gaff': 'c5'
                        } or carbon_type == {
                                'opls': 'opls_c6',
                                'gaff': 'c6'
                        }:
                            # Count the number of connected 'EWD' atoms for this carbon
                            ewd_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] in
                                ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {
                                    'opls': 'opls_140',
                                    'gaff': 'h1'
                                }  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {
                                    'opls': 'opls_xxx',
                                    'gaff': 'h2'
                                }  # 2 EWD atoms
                            elif ewd_count == 3:
                                hydrogen_type = {
                                    'opls': 'opls_xxx',
                                    'gaff': 'h3'
                                }  # 3 EWD atoms
                            else:
                                hydrogen_type = {
                                    'opls': 'opls_140',
                                    'gaff': 'hc'
                                }  # Aliphatic sp3 Hydrogen as default for c3

                        # Sp3 triangular and squared cycles
                        elif carbon_type == {
                                'opls': 'opls_CX',
                                'gaff': 'cx'
                        } or carbon_type == {
                                'opls': 'opls_CY',
                                'gaff': 'cy'
                        }:
                            # Count the number of connected 'EWD' atoms for this carbon
                            ewd_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] in
                                ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {
                                    'opls': 'opls_140',
                                    'gaff': 'h1'
                                }  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {
                                    'opls': 'opls_xxx',
                                    'gaff': 'h2'
                                }  # 2 EWD atoms
                            elif ewd_count == 3:
                                hydrogen_type = {
                                    'opls': 'opls_xxx',
                                    'gaff': 'h3'
                                }  # 3 EWD atoms
                            else:
                                hydrogen_type = {
                                    'opls': 'opls_140',
                                    'gaff': 'hc'
                                }  # Aliphatic sp3 Hydrogen as default for c3

                        # Sp2 triangular and squared cycles
                        elif carbon_type == {
                                'opls': 'opls_CU',
                                'gaff': 'cu'
                        } or carbon_type == {
                                'opls': 'opls_CV',
                                'gaff': 'cv'
                        }:
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}

                        # Sp2 carbons in aromatic or alkenes or non-aromatic cycles
                        elif carbon_type == {
                                'opls': 'opls_141',
                                'gaff': 'c2'
                        } or carbon_type == {
                                'opls': 'opls_145',
                                'gaff': 'ca'
                        }:
                            # Count the number of connected 'EWD' atoms for this carbon
                            ewd_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] in
                                ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {
                                    'opls': 'opls_146',
                                    'gaff': 'h4'
                                }  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {
                                    'opls': 'opls_xxx',
                                    'gaff': 'h5'
                                }  # 2 EWD atoms
                            else:
                                hydrogen_type = {
                                    'opls': 'opls_146',
                                    'gaff': 'ha'
                                }  # Default sp2 Hydrogen for c2 and ca types

                        elif carbon_type == {'opls': 'opls_XXX', 'gaff': 'ce'}:
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}
                        elif carbon_type == {'opls': 'opls_521', 'gaff': 'ca'}:
                            hydrogen_type = {
                                'opls': 'opls_146',
                                'gaff': 'h4'
                            }  # Hydrogens connected to C in heterocycle with N as in pyridine C-N-C

                        elif carbon_type == {
                                'opls': 'opls_XXX',
                                'gaff': 'ce'
                        } or carbon_type == {
                                'opls': 'opls_508',
                                'gaff': 'cc'
                        } or carbon_type == {
                                'opls': 'opls_XXX',
                                'gaff': 'cf'
                        } or carbon_type == {
                                'opls': 'opls_XXX',
                                'gaff': 'cd'
                        }:  # Hydrogens connected to a non-pure aromatic cycle or non_aromatic cycle
                            ewd_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] in
                                ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {
                                    'opls': 'opls_146',
                                    'gaff': 'h4'
                                }  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {
                                    'opls': 'opls_h5',
                                    'gaff': 'h5'
                                }  # 2 EWD atom
                            else:
                                hydrogen_type = {
                                    'opls': 'opls_146',
                                    'gaff': 'ha'
                                }  # Default Aliphatic sp2 Hydrogen for c2

                        # Hydrogens connected to a carbon double bonded to a sulfur
                        elif carbon_type == {'opls': 'opls_cs', 'gaff': 'cs'}:
                            ewd_count = sum(
                                1 for num in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[num]['AtomicSymbol'] in
                                ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {
                                    'opls': 'opls_146',
                                    'gaff': 'h4'
                                }
                            elif ewd_count == 2:
                                hydrogen_type = {
                                    'opls': 'opls_h5',
                                    'gaff': 'h5'
                                }
                            else:
                                hydrogen_type = {
                                    'opls': 'opls_146',
                                    'gaff': 'ha'
                                }

                        else:
                            hydrogen_type = {
                                'opls': f'opls_x{info["AtomNumber"]}',
                                'gaff': f'hx{info["AtomNumber"]}'
                            }

                        self.atom_types_dict[
                            f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

                self.atom_types_dict[f"C{info['AtomNumber']}"] = carbon_type

            # Oxygen type decision
            elif info['AtomicSymbol'] == 'O':
                if info.get('CyclicStructure') == 'none':
                    connected_symbols = set(info['ConnectedAtoms'])

                    if info['NumConnectedAtoms'] == 2 and connected_symbols == {
                            'H', 'H'
                    }:
                        oxygen_type = {
                            'opls': 'opls_154',
                            'gaff': 'ow'
                        }  # Water

                    elif info[
                            'NumConnectedAtoms'] == 2 and 'H' in connected_symbols:
                        oxygen_type = {'opls': 'opls_154', 'gaff': 'oh'}

                    elif info[
                            'NumConnectedAtoms'] == 2 and 'C' in connected_symbols:
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'O'
                               for atom in info['ConnectedAtomsNumbers']):
                            oxygen_type = {
                                'opls': 'opls_160',
                                'gaff': 'os'
                            }  # Ether group or secondary amide
                        else:
                            oxygen_type = {
                                'opls': 'opls_156',
                                'gaff': 'os'
                            }  # Carbonyl group
                    elif info['NumConnectedAtoms'] == 2:
                        oxygen_type = {
                            'opls': 'opls_156',
                            'gaff': 'os'
                        }  # General case

                    elif info['NumConnectedAtoms'] == 1:
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'C'
                               for atom in info['ConnectedAtomsNumbers']):
                            # This checks if the carbon connected to the oxygen is connected to another oxygen.
                            # It is useful to identify carboxylic acids and esters.
                            carbons = [
                                atom for atom in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[atom]['AtomicSymbol'] ==
                                'C'
                            ]
                            if any('O' in self.atom_info_dict[carbon]
                                   ['ConnectedAtoms'] for carbon in carbons):
                                oxygen_type = {
                                    'opls': 'opls_157',
                                    'gaff': 'o'
                                }  # Carboxylic acid or ester
                            else:
                                oxygen_type = {
                                    'opls': 'opls_158',
                                    'gaff': 'o'
                                }  # Aldehyde
                        elif any(
                                self.atom_info_dict[atom]['AtomicSymbol'] == 'O'
                                for atom in info['ConnectedAtomsNumbers']):
                            oxygen_type = {
                                'opls': 'opls_159',
                                'gaff': 'oo'
                            }  # Peroxide
                        else:
                            oxygen_type = {'opls': 'opls_157', 'gaff': 'o'}

                    else:
                        oxygen_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'ox{info["AtomNumber"]}'
                        }

                    self.atom_types_dict[f"O{info['AtomNumber']}"] = oxygen_type

                    # Hydrogen type assignment based on oxygen type
                    for connected_atom_number in info['ConnectedAtomsNumbers']:
                        connected_atom_info = self.atom_info_dict[
                            connected_atom_number]
                        if connected_atom_info[
                                'AtomicSymbol'] == 'H' and connected_atom_info[
                                    'NumConnectedAtoms'] == 1:
                            if oxygen_type == {
                                    'opls': 'opls_154',
                                    'gaff': 'oh'
                            }:
                                hydrogen_type = {
                                    'opls': 'opls_240',
                                    'gaff': 'ho'
                                }  # Alcohol Hydrogen
                            else:
                                hydrogen_type = {
                                    'opls': f'opls_x{info["AtomNumber"]}',
                                    'gaff': f'hox{info["AtomNumber"]}'
                                }

                            self.atom_types_dict[
                                f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

                if info.get('CyclicStructure') == 'cycle':
                    connected_symbols = set(info['ConnectedAtoms'])

                    if info['NumConnectedAtoms'] == 2 and connected_symbols == {
                            'H'
                    }:
                        oxygen_type = {
                            'opls': 'opls_154',
                            'gaff': 'oh'
                        }  # Alcohol group

                    elif info['NumConnectedAtoms'] == 2 and (
                            'C' in connected_symbols or 'S' in connected_symbols
                            or 'N' in connected_symbols):  # Ether and thioether
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'O'
                               for atom in info['ConnectedAtomsNumbers']):
                            oxygen_type = {
                                'opls': 'opls_160',
                                'gaff': 'os'
                            }  # Ether group or secondary amide
                        else:
                            oxygen_type = {
                                'opls': 'opls_156',
                                'gaff': 'os'
                            }  # Carbonyl group

                    elif info['NumConnectedAtoms'] == 1:
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'C'
                               for atom in info['ConnectedAtomsNumbers']):
                            # This checks if the carbon connected to the oxygen is connected to another oxygen.
                            # It is useful to identify carboxylic acids and esters.
                            carbons = [
                                atom for atom in info['ConnectedAtomsNumbers']
                                if self.atom_info_dict[atom]['AtomicSymbol'] ==
                                'C'
                            ]
                            if any('O' in self.atom_info_dict[carbon]
                                   ['ConnectedAtoms'] for carbon in carbons):
                                oxygen_type = {
                                    'opls': 'opls_157',
                                    'gaff': 'o'
                                }  # Carboxylic acid or ester
                            else:
                                oxygen_type = {
                                    'opls': 'opls_158',
                                    'gaff': 'o'
                                }  # Aldehyde
                        elif any(
                                self.atom_info_dict[atom]['AtomicSymbol'] == 'O'
                                for atom in info['ConnectedAtomsNumbers']):
                            oxygen_type = {
                                'opls': 'opls_159',
                                'gaff': 'oo'
                            }  # Peroxide
                        else:
                            oxygen_type = {'opls': 'opls_157', 'gaff': 'o'}

                    else:
                        oxygen_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'ox{info["AtomNumber"]}'
                        }

                    self.atom_types_dict[f"O{info['AtomNumber']}"] = oxygen_type

                    # Hydrogen type assignment based on oxygen type
                    for connected_atom_number in info['ConnectedAtomsNumbers']:
                        connected_atom_info = self.atom_info_dict[
                            connected_atom_number]
                        if connected_atom_info[
                                'AtomicSymbol'] == 'H' and connected_atom_info[
                                    'NumConnectedAtoms'] == 1:
                            if oxygen_type == {
                                    'opls': 'opls_154',
                                    'gaff': 'oh'
                            }:
                                hydrogen_type = {
                                    'opls': 'opls_240',
                                    'gaff': 'ho'
                                }  # Alcohol Hydrogen
                            else:
                                hydrogen_type = {
                                    'opls': f'opls_x{info["AtomNumber"]}',
                                    'gaff': f'hox{info["AtomNumber"]}'
                                }

                            self.atom_types_dict[
                                f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

            # Nitrogen type decision
            elif info['AtomicSymbol'] == 'N':
                connected_symbols = set(info['ConnectedAtoms'])
                connected_atoms = info['ConnectedAtoms']
                connected_atoms_numbers = info['ConnectedAtomsNumbers']

                if info.get('CyclicStructure') == 'none':

                    if info['NumConnectedAtoms'] == 4:

                        num_hydrogens = sum(
                            [1 for symbol in connected_atoms if symbol == 'H'])

                        if num_hydrogens == 4:
                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n+'}
                        elif num_hydrogens == 3:
                            nitrogen_type = {
                                'opls': 'opls_XXX',
                                'gaff': 'nz'
                            }  # Sp3 N with three hydrogen atoms
                        elif num_hydrogens == 2:
                            nitrogen_type = {
                                'opls': 'opls_XXX',
                                'gaff': 'ny'
                            }  # Sp3 N with two hydrogen atoms
                        elif num_hydrogens == 1:
                            nitrogen_type = {
                                'opls': 'opls_XXX',
                                'gaff': 'nx'
                            }  # Sp3 N with one hydrogen atom
                        else:
                            nitrogen_type = {
                                'opls': 'opls_XXX',
                                'gaff': 'n4'
                            }  # Sp3 N with four connected atoms, but no hydrogens

                    elif info['NumConnectedAtoms'] == 3:

                        # This case is highly dependent on the environment of the nitrogen
                        # Create flags to check for specific cases

                        # List of flags
                        found_nitro = False
                        found_aromatic = False
                        found_amide = False
                        found_idine = False
                        found_sp2_carbon = False
                        found_sp1_carbon = False
                        found_sp3_carbon = False

                        if sorted(info['ConnectedAtoms']) == sorted([
                                'C', 'O', 'O'
                        ]) or sorted(info['ConnectedAtoms']) == sorted(
                            ['N', 'O', 'O']):
                            found_nitro = True

                        elif 'C' in connected_symbols:
                            for atom in connected_atoms_numbers:

                                # Check for aromaticity
                                if ('cycle' in self.atom_info_dict[atom]
                                    ['CyclicStructure'] and
                                        'pure_aromatic' in self.
                                        atom_info_dict[atom]['Aromaticity']):
                                    found_aromatic = True

                                elif self.atom_info_dict[atom][
                                        'NumConnectedAtoms'] == 4:  # Connected to an sp3 carbon
                                    found_sp3_carbon = True

                                elif self.atom_info_dict[atom][
                                        'NumConnectedAtoms'] == 3:  # Connected to an sp2 carbon
                                    # Check the atoms connected to the Carbon
                                    for connected_to_carbon in self.atom_info_dict[
                                            atom]['ConnectedAtomsNumbers']:
                                        atom_symbol = self.atom_info_dict[
                                            connected_to_carbon]['AtomicSymbol']
                                        atom_connectivity = self.atom_info_dict[
                                            connected_to_carbon][
                                                'NumConnectedAtoms']

                                        # Amides and Sulfamides
                                        if (atom_symbol == 'O' and
                                                atom_connectivity == 1) or (
                                                    atom_symbol == 'S' and
                                                    atom_connectivity == 1):
                                            found_amide = True

                                        # Idines
                                        elif atom_symbol == 'N' and atom_connectivity == 2:
                                            found_idine = True

                                        # Connected to C sp1
                                        elif atom_symbol == 'C' and atom_connectivity == 2:
                                            found_sp1_carbon = True

                                        # Connected to C sp2
                                        elif atom_symbol == 'C' and atom_connectivity == 3:
                                            found_sp2_carbon = True
                                else:
                                    # n3 in GAFF but special cases n7 and n8 added in GAFF2
                                    num_hydrogens = sum([
                                        1 for symbol in connected_atoms
                                        if symbol == 'H'
                                    ])

                                    if num_hydrogens == 1:
                                        nitrogen_type = {
                                            'opls': 'opls_300',
                                            'gaff': 'n7'
                                        }  # Like n3, but with 1 attached hydrogen atom
                                    elif num_hydrogens == 2:
                                        nitrogen_type = {
                                            'opls': 'opls_300',
                                            'gaff': 'n8'
                                        }  # Like n3, but with 2 attached hydrogen atoms
                                    else:
                                        nitrogen_type = {
                                            'opls': 'opls_300',
                                            'gaff': 'n3'
                                        }  # Special case of H2-N-C Sp1

                        # General Sp3 N with three connected atoms
                        else:
                            nitrogen_type = {'opls': 'opls_300', 'gaff': 'n3'}

                        # Append all the true flags to a dictionary
                        flags = {
                            'found_nitro': found_nitro,
                            'found_aromatic': found_aromatic,
                            'found_amide': found_amide,
                            'found_idine': found_idine,
                            'found_sp2_carbon': found_sp2_carbon,
                            'found_sp1_carbon': found_sp1_carbon,
                            'found_sp3_carbon': found_sp3_carbon
                        }

                        # Now assign nitrogen types based on the flags using the following hierarchy:
                        # 1. Nitro
                        # 2. Amide
                        # 3. Idine
                        # 4. Aromatic
                        # 5. Sp2 carbon
                        # 6. Sp1 carbon
                        # 7. Sp3 carbon

                        if flags['found_nitro']:
                            nitrogen_type = {
                                'opls': 'opls_XXX',
                                'gaff': 'no'
                            }  # Nitro N
                        elif flags['found_amide']:
                            num_hydrogens = sum([
                                1 for symbol in connected_atoms if symbol == 'H'
                            ])
                            if num_hydrogens == 1:
                                nitrogen_type = {
                                    'opls': 'opls_XXX',
                                    'gaff': 'ns'
                                }
                            elif num_hydrogens == 2:
                                nitrogen_type = {
                                    'opls': 'opls_XXX',
                                    'gaff': 'nt'
                                }
                            else:
                                nitrogen_type = {
                                    'opls': 'opls_XXX',
                                    'gaff': 'n'
                                }
                        elif flags['found_idine']:
                            num_hydrogens = sum([
                                1 for symbol in connected_atoms if symbol == 'H'
                            ])
                            if num_hydrogens == 1:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nu'
                                }
                            elif num_hydrogens == 2:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nv'
                                }
                            else:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nh'
                                }
                        elif flags['found_aromatic']:
                            num_hydrogens = sum([
                                1 for symbol in connected_atoms if symbol == 'H'
                            ])
                            if num_hydrogens == 1:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nu'
                                }
                            elif num_hydrogens == 2:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nv'
                                }
                            else:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nh'
                                }
                        elif flags['found_sp2_carbon']:
                            num_hydrogens = sum([
                                1 for symbol in connected_atoms if symbol == 'H'
                            ])
                            if num_hydrogens == 1:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nu'
                                }
                            elif num_hydrogens == 2:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nv'
                                }
                            else:
                                nitrogen_type = {
                                    'opls': 'opls_901',
                                    'gaff': 'nh'
                                }
                        elif flags['found_sp1_carbon']:
                            num_hydrogens = sum([
                                1 for symbol in connected_atoms if symbol == 'H'
                            ])
                            if num_hydrogens == 1:
                                nitrogen_type = {
                                    'opls': 'opls_300',
                                    'gaff': 'n7'
                                }
                            elif num_hydrogens == 2:
                                nitrogen_type = {
                                    'opls': 'opls_300',
                                    'gaff': 'n8'
                                }
                            else:
                                nitrogen_type = {
                                    'opls': 'opls_300',
                                    'gaff': 'n3'
                                }
                        elif flags['found_sp3_carbon']:
                            num_hydrogens = sum([
                                1 for symbol in connected_atoms if symbol == 'H'
                            ])
                            if num_hydrogens == 1:
                                nitrogen_type = {
                                    'opls': 'opls_300',
                                    'gaff': 'n7'
                                }
                            elif num_hydrogens == 2:
                                nitrogen_type = {
                                    'opls': 'opls_300',
                                    'gaff': 'n8'
                                }
                            else:
                                nitrogen_type = {
                                    'opls': 'opls_300',
                                    'gaff': 'n3'
                                }

                    elif info['NumConnectedAtoms'] == 2:
                        # Check if the Nitrogen is connected to another Nitrogen with sp1 hybridization
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'N'
                               and self.atom_info_dict[atom]
                               ['NumConnectedAtoms'] == 1
                               for atom in connected_atoms_numbers):
                            nitrogen_type = {
                                'opls': 'opls_753',
                                'gaff': 'n1'
                            }  # N triple bond
                        else:
                            nitrogen_type = {
                                'opls': 'opls_903',
                                'gaff': 'n2'
                            }  # General case

                    elif info['NumConnectedAtoms'] == 1:
                        nitrogen_type = {
                            'opls': 'opls_753',
                            'gaff': 'n1'
                        }  # N triple bond

                    else:
                        nitrogen_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'nx{info["AtomNumber"]}'
                        }

                # Cyclic structures containing N
                elif info.get('CyclicStructure') == 'cycle':

                    if info['NumConnectedAtoms'] == 2 and 'pure_aromatic' in info.get(
                            'Aromaticity'):
                        nitrogen_type = {
                            'opls': 'opls_520',
                            'gaff': 'nb'
                        }  # Sp2 N in pure aromatic systems

                    elif info[
                            'NumConnectedAtoms'] == 2 and 'non_pure_aromatic' in info.get(
                                'Aromaticity'):
                        nitrogen_type = {
                            'opls': 'opls_520',
                            'gaff': 'nc'
                        }  # Sp2 N in non-pure aromatic systems

                    elif info[
                            'NumConnectedAtoms'] == 3 and 'pure_aromatic' in info.get(
                                'Aromaticity'):
                        nitrogen_type = {
                            'opls': 'opls_183',
                            'gaff': 'nb'
                        }  # Pyridine as a ligand in an organometallic complex

                    elif info[
                            'NumConnectedAtoms'] == 3 and 'non_pure_aromatic' in info.get(
                                'Aromaticity'):
                        #nitrogen_type = {'opls': 'opls_na', 'gaff': 'na'}
                        if 'C' in connected_symbols:
                            # Check for amides and sulfamides by checking the connected atoms to the carbon
                            found_CO = False  # Flag variable

                            for atom in connected_atoms_numbers:
                                for connected_to_carbon in self.atom_info_dict[
                                        atom]['ConnectedAtomsNumbers']:
                                    atom_symbol = self.atom_info_dict[
                                        connected_to_carbon]['AtomicSymbol']
                                    atom_connectivity = self.atom_info_dict[
                                        connected_to_carbon][
                                            'NumConnectedAtoms']
                                    if (atom_symbol == 'O' and atom_connectivity
                                            == 1) or (atom_symbol == 'S' and
                                                      atom_connectivity == 1):
                                        num_hydrogens = sum([
                                            1 for symbol in connected_atoms
                                            if symbol == 'H'
                                        ])
                                        if num_hydrogens == 1:
                                            nitrogen_type = {
                                                'opls': 'opls_XXX',
                                                'gaff': 'ns'
                                            }
                                        elif num_hydrogens == 2:
                                            nitrogen_type = {
                                                'opls': 'opls_XXX',
                                                'gaff': 'nt'
                                            }
                                        else:
                                            nitrogen_type = {
                                                'opls': 'opls_XXX',
                                                'gaff': 'n'
                                            }
                                        found_CO = True  # Set the flag to True
                                        break  # This will break the inner loop
                                if found_CO:  # Check the flag after each iteration of the outer loop
                                    break  # If the flag is True, break the outer loop
                            else:
                                nitrogen_type = {
                                    'opls': 'opls_na',
                                    'gaff': 'na'
                                }  # General n3 case
                        else:
                            nitrogen_type = {
                                'opls': f'opls_x{info["AtomNumber"]}',
                                'gaff': f'nx{info["AtomNumber"]}'
                            }

                    # Nitrogens in Non-aromatic cycles
                    elif info[
                            'NumConnectedAtoms'] == 3 and 'non_aromatic' in info.get(
                                'Aromaticity'):
                        if 'C' in connected_symbols:
                            # Check for amides and sulfamides by checking the connected atoms to the carbon
                            found_CO = False  # Flag variable

                            for atom in connected_atoms_numbers:
                                for connected_to_carbon in self.atom_info_dict[
                                        atom]['ConnectedAtomsNumbers']:
                                    atom_symbol = self.atom_info_dict[
                                        connected_to_carbon]['AtomicSymbol']
                                    atom_connectivity = self.atom_info_dict[
                                        connected_to_carbon][
                                            'NumConnectedAtoms']
                                    if (atom_symbol == 'O' and atom_connectivity
                                            == 1) or (atom_symbol == 'S' and
                                                      atom_connectivity == 1):
                                        num_hydrogens = sum([
                                            1 for symbol in connected_atoms
                                            if symbol == 'H'
                                        ])
                                        if num_hydrogens == 1:
                                            nitrogen_type = {
                                                'opls': 'opls_XXX',
                                                'gaff': 'ns'
                                            }
                                        elif num_hydrogens == 2:
                                            nitrogen_type = {
                                                'opls': 'opls_XXX',
                                                'gaff': 'nt'
                                            }
                                        else:
                                            nitrogen_type = {
                                                'opls': 'opls_XXX',
                                                'gaff': 'n'
                                            }
                                        found_CO = True  # Set the flag to True
                                        break  # This will break the inner loop
                                if found_CO:  # Check the flag after each iteration of the outer loop
                                    break  # If the flag is True, break the outer loop

                            if not found_CO:
                                # Check if the nitrogen is connected to a carbon sp2
                                if any(self.atom_info_dict[atom]['AtomicSymbol']
                                       == 'C' and self.atom_info_dict[atom]
                                       ['NumConnectedAtoms'] == 3
                                       for atom in connected_atoms_numbers):
                                    num_hydrogens = sum([
                                        1 for symbol in connected_atoms
                                        if symbol == 'H'
                                    ])
                                    if num_hydrogens == 1:
                                        nitrogen_type = {
                                            'opls': 'opls_901',
                                            'gaff': 'nu'
                                        }
                                    elif num_hydrogens == 2:
                                        nitrogen_type = {
                                            'opls': 'opls_901',
                                            'gaff': 'nv'
                                        }
                                    else:
                                        nitrogen_type = {
                                            'opls': 'opls_901',
                                            'gaff': 'nh'
                                        }
                                else:
                                    # n3 in GAFF but special cases n7 and n8 added in GAFF2
                                    num_hydrogens = sum([
                                        1 for symbol in connected_atoms
                                        if symbol == 'H'
                                    ])

                                    if num_hydrogens == 1:
                                        nitrogen_type = {
                                            'opls': 'opls_300',
                                            'gaff': 'n7'
                                        }  # Like n3, but with 1 attached hydrogen atom
                                    elif num_hydrogens == 2:
                                        nitrogen_type = {
                                            'opls': 'opls_300',
                                            'gaff': 'n8'
                                        }  # Like n3, but with 2 attached hydrogen atoms
                                    else:
                                        nitrogen_type = {
                                            'opls': 'opls_300',
                                            'gaff': 'n3'
                                        }  # Special case of H2-N-C Sp1

                        # Check for special RG3 and RG4 cases
                        else:
                            num_hydrogens = sum([
                                1 for symbol in connected_atoms if symbol == 'H'
                            ])
                            if num_hydrogens == 1:
                                if 3 in info.get('CycleSize'):
                                    nitrogen_type = {
                                        'opls': 'opls_n5',
                                        'gaff': 'n5'
                                    }
                                elif 4 in info.get('CycleSize'):
                                    nitrogen_type = {
                                        'opls': 'opls_n6',
                                        'gaff': 'n6'
                                    }
                                else:
                                    nitrogen_type = {
                                        'opls': 'opls_300',
                                        'gaff': 'n7'
                                    }
                            elif num_hydrogens == 0:
                                if 3 in info.get('CycleSize'):
                                    nitrogen_type = {
                                        'opls': 'opls_np',
                                        'gaff': 'np'
                                    }
                                elif 4 in info.get('CycleSize'):
                                    nitrogen_type = {
                                        'opls': 'opls_nq',
                                        'gaff': 'nq'
                                    }
                                else:
                                    nitrogen_type = {
                                        'opls': 'opls_300',
                                        'gaff': 'n3'
                                    }
                            else:
                                nitrogen_type = {
                                    'opls': f'opls_x{info["AtomNumber"]}',
                                    'gaff': f'nx{info["AtomNumber"]}'
                                }

                    elif info[
                            'NumConnectedAtoms'] == 2 and 'non_aromatic' in info.get(
                                'Aromaticity'):
                        nitrogen_type = {'opls': 'opls_903', 'gaff': 'n2'}

                    else:
                        # Add other conditions for cyclic nitrogen atoms if needed or add a default case
                        nitrogen_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'nx{info["AtomNumber"]}'
                        }

                else:
                    nitrogen_type = {
                        'opls': f'opls_x{info["AtomNumber"]}',
                        'gaff': f'nx{info["AtomNumber"]}'
                    }

                self.atom_types_dict[f"N{info['AtomNumber']}"] = nitrogen_type

                # Hydrogen type assignment based on nitrogen type
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    connected_atom_info = self.atom_info_dict[
                        connected_atom_number]
                    # Assign hydrogen type
                    if connected_atom_info[
                            'AtomicSymbol'] == 'H' and connected_atom_info[
                                'NumConnectedAtoms'] == 1:
                        hydrogen_type = {'opls': 'opls_240', 'gaff': 'hn'}

                        self.atom_types_dict[
                            f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

            # Phosphorus type decision
            elif info['AtomicSymbol'] == 'P':
                if info.get('CyclicStructure') == 'none':
                    connected_symbols = set(info['ConnectedAtoms'])

                    # Phosphate groups or phosphoric acid
                    if info['NumConnectedAtoms'] == 4 and 'O' in connected_symbols:
                        oxygen_count = list(connected_symbols).count('O')
                        if oxygen_count == 4:
                            phosphorus_type = {
                                'opls': 'opls_900P',
                                'gaff': 'p5'
                            }  # Simplified, it could be tetrahedral phosphate
                        else:
                            phosphorus_type = {
                                'opls': f'opls_x{info["AtomNumber"]}',
                                'gaff': f'px{info["AtomNumber"]}'
                            }

                    # Phosphine
                    elif info[
                            'NumConnectedAtoms'] == 3 and 'H' in connected_symbols:
                        hydrogen_count = list(connected_symbols).count('H')
                        if hydrogen_count == 3:
                            phosphorus_type = {
                                'opls': 'opls_901P',
                                'gaff': 'ph3'
                            }
                        else:
                            phosphorus_type = {
                                'opls': f'opls_x{info["AtomNumber"]}',
                                'gaff': f'px{info["AtomNumber"]}'
                            }

                    # Phosphine oxides
                    elif info[
                            'NumConnectedAtoms'] == 4 and 'O' in connected_symbols:
                        hydrogen_count = list(connected_symbols).count('H')
                        if hydrogen_count == 3:
                            phosphorus_type = {
                                'opls': 'opls_902P',
                                'gaff': 'po'
                            }
                        else:
                            phosphorus_type = {
                                'opls': f'opls_x{info["AtomNumber"]}',
                                'gaff': f'px{info["AtomNumber"]}'
                            }

                    # Phosphonates and Phosphites
                    elif info[
                            'NumConnectedAtoms'] == 3 and 'O' in connected_symbols:
                        oxygen_count = list(connected_symbols).count('O')
                        if oxygen_count == 2:
                            phosphorus_type = {
                                'opls': 'opls_903P',
                                'gaff': 'p3'
                            }  # Again simplified, could distinguish between phosphonates and phosphites
                        else:
                            phosphorus_type = {
                                'opls': f'opls_x{info["AtomNumber"]}',
                                'gaff': f'px{info["AtomNumber"]}'
                            }

                    else:
                        phosphorus_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'px{info["AtomNumber"]}'
                        }

                    self.atom_types_dict[
                        f"P{info['AtomNumber']}"] = phosphorus_type

                    # Note: Hydrogens in phosphine groups are less commonly parameterized in force fields and hence not added here.

            # Sulfur type decision
            elif info['AtomicSymbol'] == 'S':
                connected_symbols = set(info['ConnectedAtoms'])

                # S with one connected atom
                if info['NumConnectedAtoms'] == 1:
                    sulfur_type = {'opls': 'opls_920S', 'gaff': 's'}

                # S with two connected atom, involved at least one double bond
                elif info['NumConnectedAtoms'] == 2:
                    if 'H' in connected_symbols:
                        sulfur_type = {'opls': 'opls_924S', 'gaff': 'sh'}
                    elif all(
                            self.atom_info_dict[num]['AtomicSymbol'] in
                        ['C', 'N', 'S', 'O']
                            for num in info['ConnectedAtomsNumbers']
                    ):  # Both connected atoms are carbons or one carbon and one nitrogen
                        sulfur_type = {
                            'opls': 'opls_SS',
                            'gaff': 'ss'
                        }  # Thio-ether or Thio-ester
                    else:
                        sulfur_type = {'opls': 'opls_921S', 'gaff': 's2'}

                # S with three connected atoms
                elif info['NumConnectedAtoms'] == 3:
                    if any(self.atom_info_dict[num]['AtomicSymbol'] == 'C' and
                           self.atom_info_dict[num]['NumConnectedAtoms'] == 3
                           for num in info['ConnectedAtomsNumbers']):
                        sulfur_type = {'opls': 'opls_922X', 'gaff': 'sx'}
                    else:
                        sulfur_type = {'opls': 'opls_922S', 'gaff': 's4'}

                # S with four connected atoms
                elif info['NumConnectedAtoms'] == 4:
                    if any(self.atom_info_dict[num]['AtomicSymbol'] == 'C' and
                           self.atom_info_dict[num]['NumConnectedAtoms'] == 3
                           for num in info['ConnectedAtomsNumbers']):
                        sulfur_type = {'opls': 'opls_922X', 'gaff': 'sy'}
                    else:
                        sulfur_type = {'opls': 'opls_923S', 'gaff': 's6'}

                # Sp3 S connected with hydrogen
                elif info['NumConnectedAtoms'] == 4 and {
                        'H', 'H', 'H'
                } <= connected_symbols:  # use <= to check if a set is a subset
                    sulfur_type = {'opls': 'opls_924S', 'gaff': 'sh'}

                else:
                    sulfur_type = {
                        'opls': f'opls_x{info["AtomNumber"]}',
                        'gaff': f'sx{info["AtomNumber"]}'
                    }

                self.atom_types_dict[f"S{info['AtomNumber']}"] = sulfur_type

                # Hydrogen assignment in the case of thiols
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    connected_atom_info = self.atom_info_dict[
                        connected_atom_number]
                    if (connected_atom_info['AtomicSymbol'] == 'H' and
                            connected_atom_info['NumConnectedAtoms'] == 1):
                        if sulfur_type == {'opls': 'opls_924S', 'gaff': 'sh'}:
                            hydrogen_type = {'opls': 'opls_926H', 'gaff': 'hs'}
                        else:
                            hydrogen_type = {
                                'opls': f'opls_x{info["AtomNumber"]}',
                                'gaff': f'hsx{info["AtomNumber"]}'
                            }

                        self.atom_types_dict[
                            f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

            # Decision for halogens

            elif info['AtomicSymbol'] == 'Br':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'br'}
                self.atom_types_dict[
                    f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type
                
                if 'H' in info['ConnectedAtoms']:
                    hydrogen_type = {'opls': 'opls_h_x', 'gaff': 'h_x'}
                    self.atom_types_dict[
                        f"H{info['ConnectedAtomsNumbers'][info['ConnectedAtoms'].index('H')]}"
                    ] = hydrogen_type
                    self.bad_hydrogen = True

            elif info['AtomicSymbol'] == 'Cl':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'cl'}
                self.atom_types_dict[
                    f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type
                if 'H' in info['ConnectedAtoms']:
                    hydrogen_type = {'opls': 'opls_h_x', 'gaff': 'h_x'}
                    self.atom_types_dict[
                        f"H{info['ConnectedAtomsNumbers'][info['ConnectedAtoms'].index('H')]}"
                    ] = hydrogen_type
                    self.bad_hydrogen = True


            elif info['AtomicSymbol'] == 'F':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'f'}
                self.atom_types_dict[
                    f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type
                if 'H' in info['ConnectedAtoms']:
                    hydrogen_type = {'opls': 'opls_h_x', 'gaff': 'h_x'}
                    self.atom_types_dict[
                        f"H{info['ConnectedAtomsNumbers'][info['ConnectedAtoms'].index('H')]}"
                    ] = hydrogen_type
                    self.bad_hydrogen = True
        

            elif info['AtomicSymbol'] == 'I':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'i'}
                self.atom_types_dict[
                    f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type
                if 'H' in info['ConnectedAtoms']:
                    hydrogen_type = {'opls': 'opls_h_x', 'gaff': 'h_x'}
                    self.atom_types_dict[
                        f"H{info['ConnectedAtomsNumbers'][info['ConnectedAtoms'].index('H')]}"
                    ] = hydrogen_type
                    self.bad_hydrogen = True

            # Decision for Transition Metals

            elif info['AtomicSymbol'] not in [
                    'C', 'H', 'O', 'N', 'P', 'Br', 'Cl', 'F', 'I'
            ]:

                # Assign atom types based on AtomicSymbol and AtomNumber
                atom_type = {
                    'opls': f'opls_{info["AtomicSymbol"]}{info["AtomNumber"]}',
                    'gaff': f'{info["AtomicSymbol"]}{info["AtomNumber"]}'
                }
                self.atom_types_dict[
                    f"{info['AtomicSymbol']}{info['AtomNumber']}"] = atom_type
                
                if 'H' in info['ConnectedAtoms']:
                    hydrogen_type = {'opls': 'opls_h_x', 'gaff': 'h_x'}
                    self.atom_types_dict[
                        f"H{info['ConnectedAtomsNumbers'][info['ConnectedAtoms'].index('H')]}"
                    ] = hydrogen_type
                    self.bad_hydrogen = True

            else:
                # Else for the cases falling off the decision tree
                # The Hydrogen are assigned outside the main branches of the decision tree
                # Therefore, they need to be out of the else case.

                if info['AtomicSymbol'] != 'H':
                    print(
                        'Warning:',
                        f"{info['AtomicSymbol']}{info['AtomNumber']}",
                        'Has not been found in the decision tree, check it carefully'
                    )

        return self.atom_types_dict

    def extract_gaff_atom_types(self, atom_type):
        """
        Extracts GAFF atom types from the atom types dictionary.

        """

        # Initialize the list of gaff atom types
        self.gaff_atom_types = []

        # Sort atom types based on the number after the atomic symbol
        sorted_atom_types = sorted(self.atom_types_dict.keys(),
                                   key=self.get_atom_number)

        # Append the gaff atom types to the list
        for atom_type in sorted_atom_types:
            if isinstance(self.atom_types_dict[atom_type], dict):
                gaff_type = self.atom_types_dict[atom_type].get('gaff', None)
                if gaff_type:
                    self.gaff_atom_types.append(gaff_type)


    def check_alternating_atom_types(self):
        """
        This method checks the alternating atom types in GAFF2.
        Those are cc-cd and ce-cf, cg-ch and nc-nd.
        The decide_atom_types method does assign a default atom type
        as cc, ce, cg and nc, but it does not assign the other atom type.

        :return:
            The dictionary with the atom types updated with the
            alternating atom types
        """

        # Create a dictionary with the assigned bond types
        assigned_bonds = {}
        atom_types = self.gaff_atom_types
        distances = self.distance_matrix

        # Look for cc-cc, ce-ce, cg-cg and nc-nc bonds and append the distance
        # to the dictionary
        for i in range(len(atom_types)):
            for j in range(i + 1, len(atom_types)):
                if atom_types[i] == atom_types[j]:
                    if atom_types[i] in ['cc', 'ce', 'cg', 'nc']:
                        assigned_bonds[(i, j)] = distances[i][j]

        # Check distances for alternation
        for bond in assigned_bonds:
            if atom_types[bond[0]] == 'cc':
                if distances[bond[0]][bond[1]] < 1.4:
                    atom_types[bond[1]] = 'cd'
                else:
                    atom_types[bond[1]] = 'cc'
            elif atom_types[bond[0]] == 'ce':

                if distances[bond[0]][bond[1]] < 1.4:
                    atom_types[bond[1]] = 'cf'
                else:
                    atom_types[bond[1]] = 'ce'
            elif atom_types[bond[0]] == 'cg':
                if distances[bond[0]][bond[1]] < 1.4:
                    atom_types[bond[1]] = 'ch'
                else:
                    atom_types[bond[1]] = 'cg'
            elif atom_types[bond[0]] == 'nc':
                if distances[bond[0]][bond[1]] < 1.4:
                    atom_types[bond[1]] = 'nd'
                else:
                    atom_types[bond[1]] = 'nc'

        # Return an array with the new atom types
        self.gaff_atom_types = atom_types

        return self.gaff_atom_types

    # TODO: Fetch the force field file from the URL: 
    # https://github.com/openmm/openmmforcefields/blob/main/openmmforcefields/ffxml/amber/gaff/dat/gaff-2.11.dat

    def check_for_bad_assignations(self, gaff_force_field):
        '''
        Method that checks if there are any atom types that have not been
        assigned properly. To do that, the method will check if any of the 
        bonds assigned do not exist in the GAFF force field. 

        :param gaff_force_field:
            The dat file containing the force field information

        '''

        # Create a dictionary with the assigned bond types
        assigned_bonds = {}
        atom_types = self.gaff_atom_types
        distances = self.distance_matrix

        for i in range(len(atom_types)):
            for j in range(
                    i + 1,
                    len(atom_types)):  # Avoid checking the same bond twice
                # Append the bond and the atom ids to the dictionary
                if distances[i][j] != 0:
                    atomtype1 = atom_types[i]
                    atom_id1 = i + 1
                    atomtype2 = atom_types[j]
                    atom_id2 = j + 1
                    bond = (atomtype1, atomtype2)
                    assigned_bonds[bond] = (atom_id1, atom_id2)

        with open(gaff_force_field, 'r') as f:
            lines = f.readlines()

        start_bonds_section = False
        # Save a list with the bonds in the force field
        bonds = {}

        for line in lines:
            if 'hn  ho  hs  n   na  nc  nd  ne  nf  n2  n3  n4' in line:
                start_bonds_section = True
                continue
            elif line.strip() == '':
                start_bonds_section = False
            elif start_bonds_section:
                # Bonds sections always are xy-zw or x -zw, so, the atom
                # types are always in the same position
                atomtype1 = line[0:2].strip()
                element1 = atomtype1[0]
                atomtype2 = line[3:5].strip()
                element2 = atomtype2[0]
                bond = (atomtype1, atomtype2)
                bonds[bond] = (element1, element2)

        # Check if the assigned bonds are in the force field
        # TODO: Propose alternatives
        bad_bonds = {
        }  # Dictionary with the bad bonds format {(atomtype1, atomtype2): (atom_id1, atom_id2)}
        for bond in assigned_bonds:
            # Check if the bond or the reverse bond is in the force field
            if bond not in bonds and bond[::-1] not in bonds:
                print(
                    f'Warning: the bond {bond[0]}-{bond[1]} is not in the force field'
                )
                print(
                    f'Which means that {bond[0]} or {bond[1]} may not have been assigned properly'
                )
                print('Check the atom type', bond[0], 'assigned to id',
                      assigned_bonds[bond][0])
                print('Check the atom type', bond[1], 'assigned to id',
                      assigned_bonds[bond][1])
                bad_bonds[bond] = assigned_bonds[bond]
            else:
                continue

        if bad_bonds == {}:  # If there are no bad bonds, print a message
            print('All atom types have been assigned properly')

        # Print the definition of the atom types that have been assigned wrong
        first_atoms = []
        second_atoms = []

        for key, value in bad_bonds.items():
            first_atom = key[0]
            first_atoms.append(first_atom)
            second_atom = key[1]
            second_atoms.append(second_atom)
            found_first_atom = found_second_atom = False
            for line in lines[1:100]:
                # Look for first atom type
                if not found_first_atom and first_atom in line:
                    parts = line.split()
                    comment = ' '.join(parts[3:])
                    print(f'{first_atom} is {comment}')
                    found_first_atom = True
                elif not found_second_atom and second_atom in line:
                    parts = line.split()
                    comment = ' '.join(parts[3:])
                    print(f'{second_atom} is {comment}')
                    found_second_atom = True
                if found_first_atom and found_second_atom:
                    break
    
    # Main wrapper method of the class
    def generate_gaff_atomtypes(self, molecule):
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
        self.num_atoms = len(self.atomic_symbols)
        self.covalent_radii = molecule.covalent_radii_to_numpy(
        ) * bohr_in_angstrom()
        self.create_connectivity_matrix()
        self.detect_closed_cyclic_structures()
        self.create_atom_info_dict()
        self.atom_types_dict = self.decide_atom_type()
        self.extract_gaff_atom_types(self.atom_types_dict)
        self.gaff_atom_types = self.check_alternating_atom_types()

        # Printing output
        print("VeloxChem Atom Type Identification\n")
        print("-" * 40)  # Dashed line

        # Detected number of atoms
        num_atoms = len(self.atomic_symbols)
        print(f"Detected number of atoms: {num_atoms}\n")

        # Print table header
        print("{:<30} {:<20}".format("Symbol (id)", "GAFF atom type assigned"))

        # Print atom symbol, atom number, and GAFF atom type for each atom
        for i, (symbol, gaff_type) in enumerate(zip(self.atomic_symbols,
                                                    self.gaff_atom_types),
                                                start=1):
            print("{:<30} {:<20}".format(f"{symbol} ({i})", gaff_type))

        if len(self.cycle_ids) == 1:
            print(f"\nDetected {len(self.cycle_ids)} cycle of size:")
        else:
            print(f"\nDetected {len(self.cycle_ids)} cycles of sizes:")

        # Print cycle information (aromaticity) for each cycle
        for size, aromaticity in zip(self.cycle_sizes, self.aromaticity):
            if aromaticity == "non_aromatic":
                print(f"Cycle size {size}: Non-aromatic Cycle")
            elif aromaticity == "non_pure_aromatic":
                print(f"Cycle size {size}: Non-pure Aromatic Cycle")
            elif aromaticity == "pure_aromatic":
                print(f"Cycle size {size}: Pure Aromatic Cycle")
        
        if self.bad_hydrogen:
            print('\nWarning: Hydrogen type not defined in GAFF')

        return self.gaff_atom_types

    def compute_structural_features(self):
        '''
        Computes the structural features of a molecule based on its connectivity and distance matrices. 

        :return:
            A dictionary containing the lists of bonds, angles, dihedrals, impropers, 
            and pairs with their respective structural information.
        '''

        # Initialize the dictionary structure for structural features

        structural_features = {
            'bonds': [],
            'angles': [],
            'dihedrals': [],
            'impropers': [],
            'pairs': []
        }

        # Compute bonds based on the connectivity matrix and the distance matrix
        for y in range(len(self.connectivity_matrix)):  #Loop over every row
            for x in range(y + 1, len(self.connectivity_matrix)
                          ):  #Loop over every upper triangular element
                if self.connectivity_matrix[y][x] == 1:
                    bond = {
                        'atoms': (y + 1, x + 1),
                        'types':
                            (self.gaff_atom_types[y], self.gaff_atom_types[x]),
                        'distance': self.distance_matrix[y][x]
                    }
                    structural_features['bonds'].append(bond)

        # Compute angles based on the connectivity matrix and the distance matrix

        def add_angle(i, j, k):
            angle = {
                'atoms': (i + 1, j + 1, k + 1),
                'types': (self.gaff_atom_types[i], self.gaff_atom_types[j],
                          self.gaff_atom_types[k]),
                'angle': self.measure_angle(self.coordinates[i],
                                            self.coordinates[j],
                                            self.coordinates[k])
            }
            structural_features['angles'].append(angle)

        def add_dihedral(i, j, k, l):
            dihedral = {
                'atoms': (i + 1, j + 1, k + 1, l + 1),
                'types': (self.gaff_atom_types[i], self.gaff_atom_types[j],
                          self.gaff_atom_types[k], self.gaff_atom_types[l]),
                'dihedral_angle': self.measure_dihedral(self.coordinates[i],
                                                        self.coordinates[j],
                                                        self.coordinates[k],
                                                        self.coordinates[l])
            }
            structural_features['dihedrals'].append(dihedral)

        def add_improper(i, j, k, l):

            improper_dihedral = {
                'atoms': (i + 1, j + 1, k + 1, l + 1),
                'types': (self.gaff_atom_types[i], self.gaff_atom_types[j],
                          self.gaff_atom_types[k], self.gaff_atom_types[l]),
                'improper_angle': self.measure_dihedral(self.coordinates[i],
                                                        self.coordinates[j],
                                                        self.coordinates[k],
                                                        self.coordinates[l])
            }
            structural_features['impropers'].append(improper_dihedral)

        def check_improper(center, neighbours):
            #todo add better functionality for recognising sp2 atom-types
            sp2_atom_types = [
                'c ', 'cs', 'c2', 'ca', 'cp', 'cq', 'cc', 'cd', 'ce', 'cf',
                'cu', 'cv', 'cz', 'n ', 'n2', 'na', 'nb', 'nc', 'nd', 'ne',
                'nf', 'pb', 'pc', 'pd', 'pe', 'pf'
            ]
            if self.gaff_atom_types[center] in sp2_atom_types:
                # check for cycle
                cycle = False
                if self.gaff_atom_types[neighbours[
                        0]] in sp2_atom_types and self.gaff_atom_types[
                            neighbours[1]] in sp2_atom_types:
                    add_improper(center, neighbours[0], neighbours[1],
                                 neighbours[2])
                    cycle = True
                if self.gaff_atom_types[neighbours[
                        1]] in sp2_atom_types and self.gaff_atom_types[
                            neighbours[2]] in sp2_atom_types:
                    add_improper(center, neighbours[1], neighbours[2],
                                 neighbours[0])
                    cycle = True
                if self.gaff_atom_types[neighbours[
                        2]] in sp2_atom_types and self.gaff_atom_types[
                            neighbours[0]] in sp2_atom_types:
                    add_improper(center, neighbours[2], neighbours[0],
                                 neighbours[1])
                    cycle = True
                if cycle:
                    return

                #check for double bond
                if self.gaff_atom_types[neighbours[0]] in sp2_atom_types:
                    add_improper(center, neighbours[1], neighbours[2],
                                 neighbours[0])
                if self.gaff_atom_types[neighbours[1]] in sp2_atom_types:
                    add_improper(center, neighbours[2], neighbours[0],
                                 neighbours[1])
                if self.gaff_atom_types[neighbours[2]] in sp2_atom_types:
                    add_improper(center, neighbours[0], neighbours[1],
                                 neighbours[2])
                return

        for y in range(len(self.connectivity_matrix)):  #Loop over every row
            x_ids = []
            #Loop over the row, and collect the appearances of 1s
            for x in range(len(self.connectivity_matrix)):
                if self.connectivity_matrix[y][
                        x] == 1:  # If i and j are connected
                    x_ids.append(x)

            #if there are 2 bonds, that is an angle
            if len(x_ids) >= 2:
                add_angle(x_ids[0], y, x_ids[1])

            #if there are 3 bonds, it is a junction, which is made up of 3 angles and possibly an improper dihedral
            if len(x_ids) >= 3:
                add_angle(x_ids[0], y, x_ids[2])
                add_angle(x_ids[1], y, x_ids[2])
                check_improper(y, x_ids)

            #if there are 4 bonds, it is a junction, which is made up of 3 angles and no impropers
            if len(x_ids) == 4:
                add_angle(x_ids[0], y, x_ids[3])
                add_angle(x_ids[1], y, x_ids[3])
                add_angle(x_ids[2], y, x_ids[3])

            #Add dihedrals
            for x_id in x_ids:  #for all bonds in this row
                if x_id > y:  #if it is in the top triangular part
                    for y2 in range(
                            len(self.connectivity_matrix)
                    ):  #Loop over the whole column corresponding to x_id
                        if y2 != y and self.connectivity_matrix[
                                x_id,
                                y2] == 1:  #Find all the bonds (and thus skip the current row)
                            #Loop over all other x_ids
                            for x_id2 in x_ids:
                                if x_id2 != x_id:
                                    add_dihedral(y2, x_id, y, x_id2)

        return structural_features

    def generate_force_field_dict(self, ff_file_path):
        '''
        Generates a force field dictionary for a molecule based on its structural features 
        and a given GAFF force field file.

        :param ff_file_path:
            The file path to the GAFF force field data file.

        :return:
            The dictionary containing force field parameters for atom types, 
            bonds, angles, dihedrals, impropers, and non-bonded pairs.

        '''

        # Required structural features
        self.structural_features = self.compute_structural_features()

        # Initialize the dictionary to hold force field data
        force_field_data = {
            'atomtypes': {},
            'bonds': {},
            'angles': {},
            'dihedrals': {},
            'impropers': {},
            'pairs': {}
        }

        # Conversion factors
        angstrom_to_nm = 0.1
        kcalmol_to_kjmol = 4.184

        # Read the force field file
        with open(ff_file_path, 'r') as ff_file:
            ff_lines = ff_file.readlines()

        # Create a dictionary from the force field file for atom types
        atomtype_data = {
            line.split()[0]: line
            for line in ff_lines
            if re.match(r'^\S+\s+\d+\.\d+\s+\d+\.\d+\s+.*$', line)
        }

        # Create a dictionary from the force field file for sigma and epsilon values
        sigma_epsilon_data = {
            line.split()[0]: line
            for line in ff_lines
            if re.match(r'^\s*\S+\s+\d+\.\d+\s+\d+\.\d+', line)
        }

        # Parse atomtypes section and initialize force field data with unique IDs
        for index, gaff_atom_type in enumerate(
                self.gaff_atom_types,
                start=1):  # Start index at 1 since IDs start at 1
            force_field_data['atomtypes'][index] = {
                'element': self.atomic_symbols[index - 1],
                'type': gaff_atom_type,
                'mass': 0,
                'sigma': 0,
                'epsilon': 0,
                'charge': 0,
                'comment': 'undefined'
            }

            comment = ""
            # If atom type data is found in the force field, update mass and info
            if gaff_atom_type in atomtype_data:
                atom_type_match = re.match(
                    r'^(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(.+)$',
                    atomtype_data[gaff_atom_type])
                if atom_type_match:
                    mass = float(atom_type_match.group(2))
                    info = atom_type_match.group(4).strip()
                    force_field_data['atomtypes'][index].update({
                        'mass': mass,
                        'comment': 'GAFF2 ' + info
                    })
                    comment += 'GAFF2 ' + info

            # If sigma and epsilon data is found in the force field, update those values
            if gaff_atom_type in sigma_epsilon_data:
                sigma_epsilon_match = re.match(
                    r'^\s*(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)',
                    sigma_epsilon_data[gaff_atom_type])
                if sigma_epsilon_match:
                    sigma = float(sigma_epsilon_match.group(2)) * 2**(
                        -1 / 6) * 2 * angstrom_to_nm
                    epsilon = float(
                        sigma_epsilon_match.group(3)) * kcalmol_to_kjmol
                    force_field_data['atomtypes'][index].update({
                        'sigma': sigma,
                        'epsilon': epsilon
                    })
            force_field_data['atomtypes'][index].update({'comment': comment})

        # Parse bonds section
        for bond in self.structural_features['bonds']:
            bond_type_pattern = r'\s?-'.join(sorted(bond['types']))
            bond_regex = bond_type_pattern + r'\s+(\d+\.\d+)\s+(\d+\.\d+)'
            match_found = False  # Flag to check if we find a match in the force field file
            bond_key = tuple(sorted(bond['atoms']))

            for line in ff_lines:
                bond_match = re.search(bond_regex, line)
                if bond_match:
                    fc = float(bond_match.group(
                        1)) * 2 * kcalmol_to_kjmol / angstrom_to_nm**2  #unit
                    eq_distance = float(bond_match.group(2)) * angstrom_to_nm
                    # The bond key will include atom types and IDs for uniqueness
                    force_field_data['bonds'][bond_key] = {
                        # 'ids': bond['atoms'],  # Include the atom IDs
                        'types': bond['types'],
                        'fc': fc,
                        'eq': eq_distance,
                        'comment': f"GAFF2, {bond['types']}".replace("'", "")
                    }
                    match_found = True
                    break  # Exit the loop after finding the match

            if not match_found:
                # Default values if no match is found, using the computed distance from structural_features
                computed_distance = bond[
                    'distance']  # Assuming 'distance' is already in the correct units
                # The bond key will include atom types and IDs for uniqueness
                force_field_data['bonds'][bond_key] = {
                    # 'ids': bond['atoms'],  # Include the atom IDs
                    'types': bond['types'],
                    'fc': 25000.000,  # Default force constant converted to kJ/mol
                    'eq': computed_distance * angstrom_to_nm,
                    'comment': f"unknown, {bond['types']}".replace("'", "")
                }

        # Parse angles section
        for angle in self.structural_features['angles']:
            angle_types = angle['types']
            if angle_types[2] < angle_types[0]:
                angle_types = angle_types[::-1]
            angle_type_pattern = r'\s?-'.join(angle_types)
            angle_regex = angle_type_pattern + r'\s+(\d+\.\d+)\s+(\d+\.\d+)'
            match_found = False  # Flag to check if we find a match in the force field file
            angle_key = tuple(angle['atoms'])
            #Order the key so that the first index is always lower then the last one
            if angle_key[0] > angle_key[2]:
                angle_key = angle_key[::-1]

            for line in ff_lines:
                angle_match = re.search(angle_regex, line)
                if angle_match:
                    fc = float(
                        angle_match.group(1)) * 2 * kcalmol_to_kjmol  #unit
                    eq_angle_radians = float(angle_match.group(2))
                    eq_angle_degrees = eq_angle_radians  # Assuming angle is in radians, convert to degrees if necessary
                    # The angle key will include atom types and IDs for uniqueness
                    force_field_data['angles'][angle_key] = {
                        # 'ids': angle['atoms'],  # Include the atom IDs
                        'types': angle['types'],
                        'fc': fc,
                        'eq': eq_angle_degrees,
                        'comment': f"GAFF2, {angle['types']}".replace("'", "")
                    }
                    match_found = True
                    break  # Exit the loop after finding the match

            if not match_found:
                # Default values if no match is found, using the computed angle from structural_features
                computed_angle = angle[
                    'angle']  # Assuming 'angle' is already computed and stored in degrees
                # The angle key will include atom types and IDs for uniqueness
                force_field_data['angles'][angle_key] = {
                    # 'ids': angle['atoms'],  # Include the atom IDs
                    'types': angle['types'],
                    'fc': 500.00,  # Default force constant converted to kJ/mol #unit
                    'eq': computed_angle,
                    'comment': f"unknown, {angle['types']}".replace("'", "")
                }

        #todo get this to use the same code as in forcefieldgenerator, or other way around
        for dihedral in self.structural_features['dihedrals']:
            dihedral_key = tuple(dihedral['atoms'])
            #Order the key so that the first index is always lower then the last one
            if dihedral_key[0] > dihedral_key[3]:
                dihedral_key = dihedral_key[::-1]

            i, j, k, l = dihedral['atoms']
            at_1, at_2, at_3, at_4 = dihedral['types']

            patterns = [
                re.compile(r'\A' + f'{at_1}-{at_2}-{at_3}-{at_4} '),
                re.compile(r'\A' + f'{at_4}-{at_3}-{at_2}-{at_1} '),
            ]
            # Find the dihedral
            dihedral_found = False

            dihedral_ff_lines = []
            for line in ff_lines:
                matches = [re.search(p, line) for p in patterns]
                if any(matches):
                    dihedral_ff = line[11:60].strip().split()
                    if len(dihedral_ff) == 4:
                        dihedral_ff_lines.append(line)
                        dihedral_found = True

            if not dihedral_found:
                patterns = [
                    re.compile(r'\A' + f'X -{at_2}-{at_3}-X  '),
                    re.compile(r'\A' + f'X -{at_3}-{at_2}-X  '),
                ]

                dihedral_ff_lines = []
                for line in ff_lines:
                    matches = [re.search(p, line) for p in patterns]
                    if any(matches):
                        dihedral_ff = line[11:60].strip().split()
                        if len(dihedral_ff) == 4:
                            dihedral_ff_lines.append(line)
                            dihedral_found = True

            errmsg = 'ForceFieldGenerator: proper dihedral'
            errmsg += f' {at_1}-{at_2}-{at_3}-{at_4} is not available.'
            # assert_msg_critical(dihedral_found, errmsg)

            if dihedral_found:
                for line in dihedral_ff_lines:
                    dihedral_ff = line[11:60].strip().split()

                    multiplicity = int(dihedral_ff[0])
                    fc = float(dihedral_ff[1]) * 4.184 / multiplicity
                    eq = float(dihedral_ff[2])
                    # Note: negative periodicity implies multitermed dihedral
                    # See https://ambermd.org/FileFormats.php
                    try:
                        periodicity = int(dihedral_ff[3])
                    except ValueError:
                        periodicity = int(float(dihedral_ff[3]))

                    force_field_data['dihedrals'][dihedral_key] = {
                        # 'ids': dihedral['atoms'],  # Include the atom IDs
                        'types': dihedral['types'],
                        'multiplicity': multiplicity,
                        'fc': fc,
                        'eq': eq,
                        'periodicity': periodicity,
                        'comment': f"GAFF2, {dihedral['types']}".replace(
                            "'", ""),
                        'type': 9
                    }

                    if periodicity > 0:
                        break

            else:
                computed_angle = dihedral['dihedral_angle']
                # The dihedral key will include atom types and IDs for uniqueness
                force_field_data['dihedrals'][dihedral_key] = {
                    # 'ids': dihedral['atoms'],  # Include the atom IDs
                    'types': dihedral['types'],
                    'multiplicity': -1,  # Default periodicity
                    'fc': 0,  # Default force constant
                    'eq': computed_angle,
                    'periodicity': 1,  # Default phase
                    'comment': f"unknown, {dihedral['types']}".replace("'", ""),
                    'type': 9
                }

        # improper dihedrals
        for improper in self.structural_features['impropers']:
            at_1 = dihedral['types'][0]
            at_2 = dihedral['types'][1]
            at_3 = dihedral['types'][2]
            at_4 = dihedral['types'][3]

            patterns = [
                re.compile(r'\A' + f'{at_4}-{at_1}-{at_2}-{at_3} '),
                re.compile(r'\A' + f'{at_4}-{at_3}-{at_2}-{at_1} '),
                re.compile(r'\A' + f'{at_1}-{at_3}-{at_2}-{at_4} '),
                re.compile(r'\A' + f'{at_1}-{at_4}-{at_2}-{at_3} '),
                re.compile(r'\A' + f'{at_3}-{at_1}-{at_2}-{at_4} '),
                re.compile(r'\A' + f'{at_3}-{at_4}-{at_2}-{at_1} '),
            ]

            dihedral_found = False

            for line in ff_lines:
                matches = [re.search(p, line) for p in patterns]
                if any(matches):
                    dihedral_ff = line[11:60].strip().split()
                    if len(dihedral_ff) == 3:
                        dihedral_found = True
                        break

            if not dihedral_found:
                patterns = [
                    re.compile(r'\A' + f'X -{at_1}-{at_2}-{at_3} '),
                    re.compile(r'\A' + f'X -{at_3}-{at_2}-{at_1} '),
                    re.compile(r'\A' + f'X -{at_3}-{at_2}-{at_4} '),
                    re.compile(r'\A' + f'X -{at_4}-{at_2}-{at_3} '),
                    re.compile(r'\A' + f'X -{at_1}-{at_2}-{at_4} '),
                    re.compile(r'\A' + f'X -{at_4}-{at_2}-{at_1} '),
                ]

                for line in ff_lines:
                    matches = [re.search(p, line) for p in patterns]
                    if any(matches):
                        dihedral_ff = line[11:60].strip().split()
                        if len(dihedral_ff) == 3:
                            dihedral_found = True
                            break

            if not dihedral_found:
                patterns = [
                    re.compile(r'\A' + f'X -X -{at_2}-{at_3} '),
                    re.compile(r'\A' + f'X -X -{at_2}-{at_1} '),
                    re.compile(r'\A' + f'X -X -{at_2}-{at_4} '),
                ]

                for line in ff_lines:
                    matches = [re.search(p, line) for p in patterns]
                    if any(matches):
                        dihedral_ff = line[11:60].strip().split()
                        if len(dihedral_ff) == 3:
                            dihedral_found = True
                            break

            fc = float(dihedral_ff[0]) * 4.184
            eq = float(dihedral_ff[1])
            periodicity = abs(int(float(dihedral_ff[2])))

            improper_key = tuple(
                improper['atoms'])  # Use atom IDs instead of types for the key
            if improper_key[0] > improper_key[3]:
                improper_key = improper_key[::-1]

            if dihedral_found:
                force_field_data['impropers'][improper_key] = {
                    # 'ids': improper['ids'],  # Include the atom IDs
                    'types': improper['types'],
                    'periodicity': periodicity,
                    'fc': fc,
                    'eq': eq,
                    'comment': f"GAFF2, {improper['types']}".replace("'", ""),
                    'type': 4
                }
            else:
                computed_angle = improper['improper_angle']
                force_field_data['impropers'][improper_key] = {
                    # 'ids': improper['ids'],  # Include the atom IDs
                    'types': improper['types'],
                    'periodicity': 2,
                    'fc': 0,
                    'eq': computed_angle,
                    'comment': f"unknown, {improper['types']}".replace("'", ""),
                    'type': 4
                }

        # Add all dihedral pairs
        for dihedral in self.structural_features['dihedrals']:
            # Extract first and last atoms from the dihedral by IDs
            dihedral_pair_ids = (dihedral['atoms'][0], dihedral['atoms'][-1])
            # Initialize the pair in the force field data
            force_field_data['pairs'][dihedral_pair_ids] = {
                'comment': 'computed from dihedrals'
            }

        # Exclude pairs that are already present in bonds
        for par in self.structural_features['bonds']:
            bond_pair_ids = (par['atoms'][0], par['atoms'][1])
            # Remove the bond pair from the force field data if it exists
            force_field_data['pairs'].pop(bond_pair_ids, None)

        # Exclude terminal atoms present in angles
        for angle in self.structural_features['angles']:
            # Remove the first and last atom pair from pairs if it exists
            angle_pair_ids = (angle['atoms'][0], angle['atoms'][-1])
            force_field_data['pairs'].pop(angle_pair_ids, None)

        return force_field_data

    @staticmethod
    def get_atom_number(atom_type_str):
        """
        Extracts the numeric part from an atom type string.

        This static method uses regular expression to find the first sequence of digits
        in the atom type string, which typically represents the atom number.

        :param atom_type_str: 
            The atom type string containing a numeric suffix.
            
        :return:
            The numeric part extracted from the atom type string. Returns 0 if no number is found.
        """
        match = re.search(r'\d+', atom_type_str)

        return int(match.group()) if match else 0

    @staticmethod
    def measure_length(v1, v2):
        "Calculates the distance between v1 and v2"

        return np.linalg.norm(np.array(v1) - np.array(v2))

    @staticmethod
    def measure_angle(v1, v2, v3):
        "Calculates the angle v1-v2-v3 in degrees"
        rad_to_deg = 180 / mathconst_pi()
        v12 = np.array(v1) - np.array(v2)
        v32 = np.array(v3) - np.array(v2)

        return math.acos(
            np.dot(v12, v32) /
            (np.linalg.norm(v12) * np.linalg.norm(v32))) * rad_to_deg

    @staticmethod
    def measure_dihedral(v1, v2, v3, v4):
        """
        Calculates the dihedral v1-v2-v3-v4 in degrees
        https://en.wikipedia.org/wiki/Dihedral_angle#Mathematical_background
        """
        rad_to_deg = 180 / mathconst_pi()
        #Get normal vectors
        u1 = np.array(v2) - np.array(v1)
        u2 = np.array(v3) - np.array(v2)
        u3 = np.array(v4) - np.array(v3)
        n1 = np.cross(u1, u2)
        n2 = np.cross(u2, u3)

        x = np.dot(u2, np.cross(n1, n2))
        y = np.linalg.norm(u2,) * np.dot(n1, n2)

        return math.atan2(x, y) * rad_to_deg
