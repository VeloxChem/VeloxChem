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

from .molecule import Molecule
from .veloxchemlib import bohr_in_angstrom


# TODO add a small test in python_tests for AtomIdentifier

class AtomTypeIdentifier:

    """
    A class to identify atom types in a molecule using GAFF (General Amber Force Field) atom types based on 
    a VeloxChem molecule object.

    The class processes a molecule object containing the atomic coordinates of a molecule to determine the types
    of atoms according to the GAFF. It involves several steps including reading the file, determining covalent 
    radii, creating a connectivity matrix, identifying cyclic structures, and assigning atom types.

    Attributes:
        molecule (vlx.Molecule object): The molecule object from VeloxChem.
        num_atoms (int): The number of atoms in the molecule.
        atomic_symbols (list): The list of atomic symbols for the atoms in the molecule.
        coordinates (numpy.ndarray): The coordinates of the atoms in the molecule.
        covalent_radii (numpy.ndarray): The covalent radii of the atoms.
        connectivity_matrix (numpy.ndarray): The connectivity matrix representing bonded atoms in the molecule.
        distances (numpy.ndarray): The distances between bonded atoms in the molecule.
        cyclic_atoms (set): A set of atom indices that are part of cyclic structures.
        cycle_sizes (list): The list of sizes of the detected cycles.
        aromaticity (list): The list indicating the aromaticity of each cycle.
        cycles (list): The list of cycles detected in the molecule.
        cycle_ids (list): The list of unique identifiers for each cycle.
        atom_cycle_info (dict): A dictionary containing cycle information for each atom.
        self.atom_info_dict (dict): A dictionary holding detailed atom information used for atom type decision.
        self.atom_types_dict (dict): A dictionary with assigned atom types for each atom.
        gaff_atom_types (list): A list of GAFF atom types extracted from `self.atom_types_dict`.

    Methods:
        create_connectivity_matrix: Creates a matrix that represents bonded atoms in the molecule.
        detect_closed_cyclic_structures: Identifies cyclic structures within the molecule.
        create_self.atom_info_dict: Constructs a dictionary with information about each atom.
        decide_atom_type: Assigns atom types based on the gathered information.
        extract_gaff_atom_types: Extracts the GAFF atom types from the assigned atom types.
        get_atom_number: Utility method to extract the atom number from a string.
        to_array: Generates an array with the GAFF atom types.

    The class should be initialized with the path to an XYZ file, after which it will automatically process 
    the file and set the corresponding attributes.

    Example:
        >>> atom_identifier = AtomTypeIdentifier('molecule.xyz')
        >>> print(atom_identifier.gaff_atom_types)
    """

    def __init__(self):
        """
        Initializes the AtomTypeIdentifier instance.

        Args:
            self
        """

    def create_connectivity_matrix(self, coordinates, covalent_radii, factor=1.3):
        """
        Creates a connectivity matrix for the molecule based on the atomic coordinates
        and covalent radii, determining which atoms are bonded.

        This method iterates through pairs of atoms and calculates the distance between
        them. If this distance is less than or equal to the sum of their covalent radii
        scaled by a factor, it is considered a bond, and the connectivity matrix is updated
        accordingly. The method also constructs a corresponding distance matrix with the
        actual distances between connected atoms.

        Parameters:
            factor (float, optional): A scaling factor for the covalent radii to account for
                                      the bond threshold. Defaults to 1.3.

        Returns:
            tuple: A tuple containing two 2D numpy arrays:
                   - The first array is the connectivity matrix with 1s indicating bonded atom pairs.
                   - The second array is the distance matrix with actual distances between atoms.

        Side Effects:
            - Updates the 'connectivity_matrix' and 'distance_matrix' attributes of the class.

        Example:
            >>> atom_identifier = AtomTypeIdentifier('molecule.xyz')
            >>> atom_identifier.create_connectivity_matrix()
            >>> print(atom_identifier.connectivity_matrix)
            >>> print(atom_identifier.distance_matrix)
        """
        num_atoms = len(self.coordinates)
        self.connectivity_matrix = np.zeros((num_atoms, num_atoms), dtype=int)
        self.distance_matrix = np.zeros((num_atoms, num_atoms))

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                distance = self.calculate_3d_distance(self.coordinates[i], self.coordinates[j])
                covalent_distance = self.covalent_radii[i] + self.covalent_radii[j]
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

        This method uses the atomic symbols, connectivity matrix, and coordinates
        to plot a 3D representation of the molecule, showing bonds between atoms.

        Side Effects:
            - Generates a 3D plot showing the connected atoms in the molecule.

        Example:
            >>> atom_identifier = AtomTypeIdentifier('molecule.xyz')
            >>> atom_identifier.plot_connectivity_map_3D()
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
            ax.text(self.coordinates[i, 0], self.coordinates[i, 1], self.coordinates[i, 2], symbol, fontsize=12, ha='center', va='center')

        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Connectivity Map in 3D')
        plt.show()

        pass


    def is_sp2_carbon(self, atom_idx):
        """
        Determines if a given atom, identified by its index, is sp2 hybridized.

        Args:
            atom_idx (int): Index of the atom in the molecule.

        Returns:
            bool: True if the atom is sp2 hybridized, False otherwise.
        """
        return self.atomic_symbols[atom_idx] == "C" and len(list(self.graph.neighbors(atom_idx))) == 3

    def detect_closed_cyclic_structures(self,atomic_symbols, connectivity_matrix, distance_matrix):
        """
        Detects closed cyclic structures in a molecule and determines their aromaticity.

        This method analyzes the graph of atoms and their connectivity to identify cycles,
        determine the size of each cycle, and classify them based on aromaticity criteria.

        Returns:
            tuple: Contains sets and lists of cyclic atoms, cycle sizes, aromaticity,
            unique cycles, cycle ids, and atom cycle information.
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
        filtered_cycles = [cycle for cycle in all_cycles if 3 <= len(cycle) <= 6]

        # Remove super-cycles (cycles that contain smaller cycles) 
        reduced_cycles = filtered_cycles[:]  
        cycles_to_remove = set()
        for i, cycle in enumerate(filtered_cycles):
            for j, larger_cycle in enumerate(filtered_cycles):
                if len(cycle) < len(larger_cycle) and set(cycle).issubset(set(larger_cycle)):
                    cycles_to_remove.add(tuple(larger_cycle))

        # Convert cycles_to_remove to list of lists and then to a set of tuples for faster lookup
        cycles_to_remove = {tuple(c) for c in cycles_to_remove}

        # Rebuild reduced_cycles excluding the ones in cycles_to_remove
        reduced_cycles = [cycle for cycle in reduced_cycles if tuple(cycle) not in cycles_to_remove]


        self.cyclic_atoms = set()
        self.cycle_sizes = []
        self.aromaticity = []
        self.cycle_ids = []
        self.atom_cycle_info = {}

        # Assignation of aromaticity to all the reduced cycles
        for cycle_id, cycle in enumerate(reduced_cycles):
            self.cyclic_atoms.update(cycle)

            size = len(cycle)
            self.cycle_sizes.append(size)

            cycle_elements = [self.atomic_symbols[i] for i in cycle]
            cc_bond_distances = [self.distances[cycle[i]][cycle[(i + 1) % size]] for i in range(size) if self.atomic_symbols[cycle[i]] == "C" and self.atomic_symbols[cycle[(i + 1) % size]] == "C"]

            max_distance = max(cc_bond_distances) if cc_bond_distances else 0
            min_distance = min(cc_bond_distances) if cc_bond_distances else 0

            # Aromaticity determination logic
            if size == 6 and all(elem in ["C", "N"] for elem in cycle_elements):
                if all(self.is_sp2_carbon(atom_idx) for atom_idx in cycle if self.atomic_symbols[atom_idx] == "C"):
                    aro = "pure_aromatic" if max_distance - min_distance <= 0.08 else "non_pure_aromatic"
                else:
                    aro = "non_aromatic"
            elif size == 5 and all(self.is_sp2_carbon(atom_idx) for atom_idx in cycle if self.atomic_symbols[atom_idx] == "C"):
                if "S" in cycle_elements:
                    aro = "non_pure_aromatic" if len(list(self.graph.neighbors(cycle_elements.index("S")))) == 2 else "non_aromatic"
                else:
                    aro = "non_pure_aromatic" if "N" in cycle_elements or "O" in cycle_elements else "non_aromatic"
            elif size == 4:
                aro = "non_pure_aromatic" if all(self.is_sp2_carbon(atom_idx) for atom_idx in cycle if self.atomic_symbols[atom_idx] == "C") else "non_aromatic"
            else:
                aro = "non_aromatic"

            self.aromaticity.append(aro)

            for atom in cycle:
                if atom not in self.atom_cycle_info:
                    self.atom_cycle_info[atom] = {'sizes': [], 'aromaticities': []}
                self.atom_cycle_info[atom]['sizes'].append(size)
                self.atom_cycle_info[atom]['aromaticities'].append(aro)

            self.cycle_ids.append(cycle_id)


        # Additional logic for reassignment of aromaticity in special cases where 3 atoms are shared with aromatic rings.
        for index, cycle in enumerate(reduced_cycles):
            # Check if all carbons in the cycle are sp2
            all_carbons_sp2 = all(self.is_sp2_carbon(atom_idx) for atom_idx in cycle if self.atomic_symbols[atom_idx] == "C")
            if self.cycle_sizes[index] == 5 and self.aromaticity[index] == 'non_aromatic' and all_carbons_sp2:
                count_pure_aromatic_atoms = sum(1 for atom in cycle if 'pure_aromatic' in self.atom_cycle_info[atom]['aromaticities'])
                if count_pure_aromatic_atoms >= 3:
                    self.aromaticity[index] = 'non_pure_aromatic'
                    for atom in cycle:
                        self.atom_cycle_info[atom]['aromaticities'] = ['non_pure_aromatic' if a == 'non_aromatic' else a for a in self.atom_cycle_info[atom]['aromaticities']]

        return self.cyclic_atoms, self.cycle_sizes, self.aromaticity, reduced_cycles, self.cycle_ids, self.atom_cycle_info

    def create_atom_info_dict(self, atomic_symbols, connectivity_matrix, distances, cyclic_atoms, cycle_sizes, aromaticity, cycles, cycle_ids, atom_cycle_info):
        """
        Creates a dictionary containing detailed information for each atom in the molecule.

        This method compiles the atomic symbol, atom number, number of connected atoms, symbols of connected atoms,
        atom numbers of connected atoms, distances to connected atoms, and cycle information into a structured dictionary.

        Returns:
            dict: A dictionary where each key is an atom number and each value is another dictionary of atom information.
        """
        self.atom_info = {}
        for i, symbol in enumerate(self.atomic_symbols):
            num_connected_atoms = np.sum(self.connectivity_matrix[i])
            connected_atoms_numbers = [j + 1 for j in range(len(self.atomic_symbols)) if self.connectivity_matrix[i][j] == 1]
            connected_atoms_symbols = [self.atomic_symbols[j] for j in range(len(self.atomic_symbols)) if self.connectivity_matrix[i][j] == 1]
            connected_atoms_distances = [self.distances[i][j] for j in range(len(self.atomic_symbols)) if self.connectivity_matrix[i][j] == 1]

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
                info["CycleNumber"] = [self.cycle_ids[self.cycles.index(c)] for c in self.cycles if i in c]
            else:
                info["CyclicStructure"] = "none"

            self.atom_info[i + 1] = info

        return self.atom_info

    def decide_atom_type(self, atom_info):
        """
        Analyzes the molecular structure information to assign atom types to each atom in the molecule.

        This method traverses through the atom information dictionary created from the molecular structure data
        and applies a series of rules to determine the appropriate atom type for each atom. The rules consider
        factors such as the atom's chemical environment, its connectivity to other atoms, and whether it is part of
        a cyclic structure.

        The atom types are determined for both the OPLS and GAFF force fields, which are used in various molecular
        dynamics simulations and computational chemistry analyses. The assignment process is crucial for the accurate
        representation of the molecule in these simulations.

        Returns:
            dict: A dictionary where each key corresponds to an atom identifier (e.g., "C1" for the first carbon atom),
                  and each value is another dictionary containing the 'opls' and 'gaff' force field identifiers for the atom.
        """

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
                if info.get('CyclicStructure') == 'cycle' and 'pure_aromatic' in info.get('Aromaticity'):
                    
                    if info['NumConnectedAtoms'] == 3: 
                        # Check for identifying biphenyls
                        connected_carbons_in_diff_cycle_and_pure_aromatic = [
                            connected_atom_number for connected_atom_number in info['ConnectedAtomsNumbers']
                            if self.atom_info_dict[connected_atom_number].get('AtomicSymbol') == 'C' and 
                            self.atom_info_dict[connected_atom_number].get('CycleNumber') and  # Ensure it's not empty
                            not set(self.atom_info_dict[connected_atom_number].get('CycleNumber')) & set(info.get('CycleNumber')) and 
                            'pure_aromatic' in self.atom_info_dict[connected_atom_number].get('Aromaticity')
                        ]

                        # If the list is not empty, set the connected_carbon_atom to the first element. Else, set it to None.
                        connected_carbon_atom = connected_carbons_in_diff_cycle_and_pure_aromatic[0] if connected_carbons_in_diff_cycle_and_pure_aromatic else None

                        if connected_symbols == {'C', 'H'}: 
                            carbon_type = {'opls': 'opls_145', 'gaff': 'ca'}
                        elif connected_symbols == {'C', 'N', 'H'}:
                            carbon_type = {'opls': 'opls_521', 'gaff': 'ca'}
                                    
                        elif connected_carbon_atom:
                            carbon_type = {'opls': 'opls_521', 'gaff': 'cp'}

                            index_in_distances = info['ConnectedAtomsNumbers'].index(connected_carbon_atom)
                            d = info['ConnectedAtomsDistances'][index_in_distances]
                            
                            if d >= 1.4685:
                                biphenyl_carbon = {'opls': 'opls_521', 'gaff': 'cp'}
                            elif d <= 1.4685:
                                biphenyl_carbon = {'opls': 'opls_CQ', 'gaff': 'cq'}
                            
                            # Store the atomtype in the dictionary for the biphenyl carbon
                            self.atom_types_dict[f"C{connected_carbon_atom}"] = biphenyl_carbon
                        else:
                            carbon_type = {'opls': 'opls_145', 'gaff': 'ca'}

                    elif info['NumConnectedAtoms'] == 4 and connected_symbols == {'C', 'H'}:
                        carbon_type = {'opls': 'opls_135', 'gaff': 'c3'}
                    else:
                        carbon_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'cx{info["AtomNumber"]}'}
                        
                # Default assignation for carbons in non-pure aromatic cyclic structures
                elif info.get('CyclicStructure') == 'cycle' and 'non_pure_aromatic' in info['Aromaticity']:
                #elif info.get('CyclicStructure') == 'cycle' and info['Aromaticity'] == ['non_pure_aromatic']:
            
                    if 'O' in connected_symbols:
                    # Directly loop through the connected atom numbers
                        for connected_atom_number in info['ConnectedAtomsNumbers']:
                            if self.atom_info_dict[connected_atom_number]['AtomicSymbol'] == 'O' and self.atom_info_dict[connected_atom_number]['CyclicStructure'] == 'none':
                                if self.atom_info_dict[connected_atom_number]['NumConnectedAtoms'] == 1:
                                    carbon_type = {'opls': 'opls_235', 'gaff': 'c'}  # Carbonyl Carbon
                                    break  # Exit the loop once the carbonyl carbon is identified
                            else:
                                carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}
                                    
                    elif 'C' in connected_symbols and info['NumConnectedAtoms'] == 3 and any(size in info['CycleSize'] for size in {4, 5, 6}):  

                        carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}  # Default assignment

                        # Immediately assign C atom types based on distances
                        connected_atoms_numbers = info['ConnectedAtomsNumbers']

                        # Current atom's cycle number (will be 'none' if not in a cycle)
                        current_cycle_number = info.get('CycleNumber')

                        for neighboring_atom_number, distance in zip(connected_atoms_numbers, connected_distances):
                            # Check if the neighbour is also a C
                            if self.atom_info_dict[neighboring_atom_number]['AtomicSymbol'] != 'C':
                                continue
                            key_neigh = f"C{neighboring_atom_number}"
                            # Check if the C atom has been previously assigned
                            if key_neigh in self.atom_types_dict:
                                existing_type = self.atom_types_dict[key_neigh]
                                #carbon_type = None
                                # For neighboring atom type {'opls': 'opls_XXX', 'gaff': 'cc'}
                                if existing_type == {'opls': 'opls_508', 'gaff': 'cc'}:
                                    if distance <= 1.4:
                                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'cd'}
                                    elif distance > 1.4:
                                        carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}
                                
                                # For neighboring atom type {'opls': 'opls_XXX', 'gaff': 'cf'}
                                elif existing_type == {'opls': 'opls_XXX', 'gaff': 'cd'}:
                                    if distance <= 1.4:
                                        carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}
                                    elif distance > 1.4:
                                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'cd'}
                                
                                else:
                                    continue
                            else:

                                neighboring_cycle_number = self.atom_info_dict[neighboring_atom_number].get('CycleNumber', 'none')

                                # Check if neighboring atom is not in a cycle or belongs to another cycle ID
                                if neighboring_cycle_number == 'none' or neighboring_cycle_number != current_cycle_number:
                                    continue
                                else:
                                    if distance <= 1.4: # Double bonded
                                        conn_carbon_type = {'opls': 'opls_XXX', 'gaff': 'cd'}
                                        
                                    elif distance >= 1.4: # Single bonded
                                        conn_carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}
                                    
                                    # Store the atomtype in the dictionary 
                                    self.atom_types_dict[f"C{neighboring_atom_number}"] = conn_carbon_type
                                    
                                    # If there is an H connected to this carbon, assign it now
                                    for connected_atom_number in self.atom_info_dict[neighboring_atom_number]['ConnectedAtomsNumbers']:
                                        connected_atom_info = self.atom_info_dict[connected_atom_number]
                                        if connected_atom_info['AtomicSymbol'] == 'H' and connected_atom_info['NumConnectedAtoms'] == 1:
                                            ewd_atoms = ['N', 'Br', 'Cl', 'I', 'F', 'S', 'O']
                                            ewd_count = sum(1 for num in self.atom_info_dict[neighboring_atom_number]['ConnectedAtomsNumbers'] 
                                                            if self.atom_info_dict[num]['AtomicSymbol'] in ewd_atoms)
                                            if ewd_count == 1:
                                                hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}  # 1 EWD atom
                                            elif ewd_count == 2:
                                                hydrogen_type = {'opls': 'opls_h5', 'gaff': 'h5'}  # 2 EWD atom
                                            else:
                                                hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}  # Default Aliphatic sp2 Hydrogen for c2
                                            
                                            # Store it in the dictionary
                                            self.atom_types_dict[f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type
                                    
                                    # Next iteration
                                    continue
                    else:
                        carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}

                # Chain structures
                elif info.get('CyclicStructure') == 'cycle' and 'non_aromatic' in info['Aromaticity']:

                    if info['NumConnectedAtoms'] == 4 and 3 in info['CycleSize']:
                        carbon_type = {'opls': 'opls_CX', 'gaff': 'cx'}
                    elif info['NumConnectedAtoms'] == 4 and 4 in info['CycleSize']:
                        carbon_type = {'opls': 'opls_CY', 'gaff': 'cy'}
                    elif info['NumConnectedAtoms'] == 4 and 5 in info['CycleSize']:
                        carbon_type = {'opls': 'opls_c5', 'gaff': 'c5'}
                    elif info['NumConnectedAtoms'] == 4 and 6 in info['CycleSize']:
                        carbon_type = {'opls': 'opls_c6', 'gaff': 'c6'}
                    elif info['NumConnectedAtoms'] == 3 and 3 in info['CycleSize']:
                        carbon_type = {'opls': 'opls_CU', 'gaff': 'cu'}
                    elif info['NumConnectedAtoms'] == 3 and 4 in info['CycleSize']:
                        carbon_type = {'opls': 'opls_CV', 'gaff': 'cv'}
                    elif info['NumConnectedAtoms'] == 3 and 5 in info['CycleSize']:

                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'ce'} # Inner Sp2 carbons in conjugated systems - Default

                        # Immediately assign C atom types based on distances
                        connected_atoms_numbers = info['ConnectedAtomsNumbers']

                        # Current atom's cycle number (will be 'none' if not in a cycle)
                        current_cycle_number = info.get('CycleNumber')

                        for neighboring_atom_number, distance in zip(connected_atoms_numbers, connected_distances):
                            # Check if the neighbour is also a C
                            if self.atom_info_dict[neighboring_atom_number]['AtomicSymbol'] != 'C':
                                continue
                            key_neigh = f"C{neighboring_atom_number}"
                            
                            # Check if the C atom has been previously assigned
                            if key_neigh in self.atom_types_dict:
                                existing_type = self.atom_types_dict[key_neigh]
                                #carbon_type = None
                                # For neighboring atom type {'opls': 'opls_XXX', 'gaff': 'ce'}
                                if existing_type == {'opls': 'opls_XXX', 'gaff': 'ce'}:
                                    if distance <= 1.4:
                                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'cf'}
                                    elif distance > 1.4:
                                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'ce'}
                                
                                # For neighboring atom type {'opls': 'opls_XXX', 'gaff': 'cf'}
                                elif existing_type == {'opls': 'opls_XXX', 'gaff': 'cf'}:
                                    if distance <= 1.4:
                                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'ce'}
                                    elif distance > 1.4:
                                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'cf'}
                                
                                else:
                                    continue
                            else:
                                neighboring_cycle_number = self.atom_info_dict[neighboring_atom_number].get('CycleNumber', 'none')

                                # Check if neighboring atom is not in a cycle or belongs to another cycle ID
                                if neighboring_cycle_number == 'none' or neighboring_cycle_number != current_cycle_number:
                                    continue
                                else:
                                    if distance <= 1.4: # Double bonded
                                        conn_carbon_type = {'opls': 'opls_XXX', 'gaff': 'cf'}
                                        
                                    elif distance > 1.4: # Single bonded
                                        conn_carbon_type = {'opls': 'opls_XXX', 'gaff': 'ce'}
                                    
                                    # Store the atomtype in the dictionary 
                                    self.atom_types_dict[f"C{neighboring_atom_number}"] = conn_carbon_type

                                    # If there is an H connected to this carbon, assign it now
                                    for connected_atom_number in self.atom_info_dict[neighboring_atom_number]['ConnectedAtomsNumbers']:
                                        connected_atom_info = self.atom_info_dict[connected_atom_number]
                                        if connected_atom_info['AtomicSymbol'] == 'H' and connected_atom_info['NumConnectedAtoms'] == 1:
                                            ewd_atoms = ['N', 'Br', 'Cl', 'I', 'F', 'S', 'O']
                                            ewd_count = sum(1 for num in self.atom_info_dict[neighboring_atom_number]['ConnectedAtomsNumbers'] 
                                                            if self.atom_info_dict[num]['AtomicSymbol'] in ewd_atoms)
                                            if ewd_count == 1:
                                                hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}  # 1 EWD atom
                                            elif ewd_count == 2:
                                                hydrogen_type = {'opls': 'opls_h5', 'gaff': 'h5'}  # 2 EWD atom
                                            else:
                                                hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}  # Default Aliphatic sp2 Hydrogen for c2
                                            
                                            # Store it in the dictionary
                                            self.atom_types_dict[f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type
                                    
                                    # Next iteration
                                    continue

                    # Cases for general Non-Aromatic cycles bigger than 5
                    elif info['NumConnectedAtoms'] == 3:
                        if 'O' in connected_symbols:
                            # Directly loop through the connected atom numbers
                            for connected_atom_number in info['ConnectedAtomsNumbers']:
                                if self.atom_info_dict[connected_atom_number]['AtomicSymbol'] == 'O':
                                    if self.atom_info_dict[connected_atom_number]['NumConnectedAtoms'] == 1:
                                        carbon_type = {'opls': 'opls_235', 'gaff': 'c'}  # Carbonyl Carbon
                                        break  # Exit the loop once the carbonyl carbon is identified
                                    else:
                                        carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}  # Alcohol

                        elif 'C' in connected_symbols:
                            # Count the number of sp2-hybridized carbons connected to the current carbon
                            sp2_carbon_count = sum(1 for num in info['ConnectedAtomsNumbers'] 
                                                if self.atom_info_dict[num]['AtomicSymbol'] == 'C' and self.atom_info_dict[num]['NumConnectedAtoms'] == 3)
                            sp1_carbon_count = sum(1 for num in info['ConnectedAtomsNumbers'] 
                                                if self.atom_info_dict[num]['AtomicSymbol'] == 'C' and self.atom_info_dict[num]['NumConnectedAtoms'] == 2)
                            if sp2_carbon_count + sp1_carbon_count == 2 or sp2_carbon_count + sp1_carbon_count == 3: # If the current carbon is connected to 2 sp2 carbons
                                carbon_type = {'opls': 'opls_XXX', 'gaff': 'ce'} # Inner Sp2 carbons in conjugated systems
                            else:
                                carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}  # Aliphatic sp2 Carbon
                        else:
                            carbon_type = {'opls': 'opls_141', 'gaff': 'c2'} #Generic sp2 C
                    
                    else:
                        carbon_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'cx{info["AtomNumber"]}'}

                elif info.get('CyclicStructure') == 'none':
                    if info['NumConnectedAtoms'] == 4:
                        if 'C' and 'H' in connected_symbols:
                            carbon_type = {'opls': 'opls_135', 'gaff': 'c3'}  # Aliphatic sp3 Carbon
                        elif 'O' in connected_symbols and connected_symbols != {'C', 'O', 'O', 'O'}:
                            carbon_type = {'opls': 'opls_135', 'gaff': 'c3'}  # Carbon bonded to Oxygen (but not in carboxylate)
                        else:
                            carbon_type = {'opls': 'opls_135', 'gaff': 'c3'} # General case, anything connected to it
                            
                    elif info['NumConnectedAtoms'] == 3:
                        if 'O' in connected_symbols:
                            # Directly loop through the connected atom numbers
                            for connected_atom_number in info['ConnectedAtomsNumbers']:
                                if self.atom_info_dict[connected_atom_number]['AtomicSymbol'] == 'O':
                                    if self.atom_info_dict[connected_atom_number]['NumConnectedAtoms'] == 1:
                                        carbon_type = {'opls': 'opls_235', 'gaff': 'c'}  # Carbonyl Carbon
                                        break  # Exit the loop once the carbonyl carbon is identified
                                    else:
                                        carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}  # Alcohol

                        elif 'C' in connected_symbols:
                            # Count the number of sp2-hybridized carbons connected to the current carbon
                            sp2_carbon_count = sum(1 for num in info['ConnectedAtomsNumbers'] 
                                                if self.atom_info_dict[num]['AtomicSymbol'] == 'C' and self.atom_info_dict[num]['NumConnectedAtoms'] == 3)
                            sp1_carbon_count = sum(1 for num in info['ConnectedAtomsNumbers'] 
                                                if self.atom_info_dict[num]['AtomicSymbol'] == 'C' and self.atom_info_dict[num]['NumConnectedAtoms'] == 2)
                            if sp2_carbon_count + sp1_carbon_count == 2 or sp2_carbon_count + sp1_carbon_count == 3: # If the current carbon is connected to 2 sp2 carbons
                                carbon_type = {'opls': 'opls_XXX', 'gaff': 'ce'} # Inner Sp2 carbons in conjugated systems
                            else:
                                carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}  # Aliphatic sp2 Carbon
                        else:
                            carbon_type = {'opls': 'opls_141', 'gaff': 'c2'} #Generic sp2 C
                            
                    elif info['NumConnectedAtoms'] == 2:
                        if 'O' in connected_symbols:
                            carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}  # Carbon in carbonyl group or acid anhydride
                        else:
                            # Count the number of sp2-hybridized carbons connected to the current carbon
                            sp2_carbon_count = sum(1 for num in info['ConnectedAtomsNumbers'] 
                                                if self.atom_info_dict[num]['AtomicSymbol'] == 'C' and self.atom_info_dict[num]['NumConnectedAtoms'] == 3)
                            sp1_carbon_count = sum(1 for num in info['ConnectedAtomsNumbers'] 
                                                if self.atom_info_dict[num]['AtomicSymbol'] == 'C' and self.atom_info_dict[num]['NumConnectedAtoms'] == 2)
                            n_count = sum(1 for num in info['ConnectedAtomsNumbers'] 
                                                if self.atom_info_dict[num]['AtomicSymbol'] == 'N' and self.atom_info_dict[num]['NumConnectedAtoms'] == 1)
                            n2_count = sum(1 for num in info['ConnectedAtomsNumbers'] 
                                                if self.atom_info_dict[num]['AtomicSymbol'] == 'N' and self.atom_info_dict[num]['NumConnectedAtoms'] == 2)
                            if sp2_carbon_count + sp1_carbon_count + n_count == 2: # If the current carbon is connected to 2 sp2 carbons
                                carbon_type = {'opls': 'opls_XXX', 'gaff': 'cg'} 
                            elif sp2_carbon_count + sp1_carbon_count + n2_count == 2:
                                carbon_type = {'opls': 'opls_XXX', 'gaff': 'ch'}
                            else:
                                carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}
                    elif info['NumConnectedAtoms'] == 1:
                        carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}
                    else:
                        carbon_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'cx{info["AtomNumber"]}'}
                else:
                    carbon_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'cx{info["AtomNumber"]}'}
                    
                # Assignment for the Hydrogens linked to the carbons
                # ewd = Electron withdrawing atoms
                ewd_atoms = ['N', 'Br', 'Cl', 'I', 'F', 'S', 'O']

                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    connected_atom_info = self.atom_info_dict[connected_atom_number]
                    if connected_atom_info['AtomicSymbol'] == 'H' and connected_atom_info['NumConnectedAtoms'] == 1:
                        # Sp1 Carbon
                        if carbon_type == {'opls': 'opls_235', 'gaff': 'c1'}:
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}
                        # Carbonyl
                        elif carbon_type == {'opls': 'opls_235', 'gaff': 'c'}:
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'} # Hydrogen in a Carbonyl
                        # Sp3 Carbon types
                        elif carbon_type == {'opls': 'opls_135', 'gaff': 'c3'} or carbon_type == {'opls': 'opls_c5', 'gaff': 'c5'} or carbon_type == {'opls': 'opls_c6', 'gaff': 'c6'}:
                            # Count the number of connected 'EWD' atoms for this carbon
                            ewd_count = sum(1 for num in info['ConnectedAtomsNumbers'] if self.atom_info_dict[num]['AtomicSymbol'] in ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {'opls': 'opls_140', 'gaff': 'h1'}  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {'opls': 'opls_xxx', 'gaff': 'h2'}  # 2 EWD atoms
                            elif ewd_count == 3:
                                hydrogen_type = {'opls': 'opls_xxx', 'gaff': 'h3'}  # 3 EWD atoms
                            else:
                                hydrogen_type = {'opls': 'opls_140', 'gaff': 'hc'}  # Aliphatic sp3 Hydrogen as default for c3

                        # Sp3 triangular and squared cycles
                        elif carbon_type == {'opls': 'opls_CX', 'gaff': 'cx'} or carbon_type == {'opls': 'opls_CY', 'gaff': 'cy'}:
                            # Count the number of connected 'EWD' atoms for this carbon
                            ewd_count = sum(1 for num in info['ConnectedAtomsNumbers'] if self.atom_info_dict[num]['AtomicSymbol'] in ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {'opls': 'opls_140', 'gaff': 'h1'}  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {'opls': 'opls_xxx', 'gaff': 'h2'}  # 2 EWD atoms
                            elif ewd_count == 3:
                                hydrogen_type = {'opls': 'opls_xxx', 'gaff': 'h3'}  # 3 EWD atoms
                            else:
                                hydrogen_type = {'opls': 'opls_140', 'gaff': 'hc'}  # Aliphatic sp3 Hydrogen as default for c3
                                
                        # Sp2 triangular and squared cycles
                        elif carbon_type == {'opls': 'opls_CU', 'gaff': 'cu'} or carbon_type == {'opls': 'opls_CV', 'gaff': 'cv'}:
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}

                        # Sp2 carbons in aromatic or alkenes or non-aromatic cycles
                        elif carbon_type == {'opls': 'opls_141', 'gaff': 'c2'} or carbon_type == {'opls': 'opls_145', 'gaff': 'ca'}: 
                            # Count the number of connected 'EWD' atoms for this carbon
                            ewd_count = sum(1 for num in info['ConnectedAtomsNumbers'] if self.atom_info_dict[num]['AtomicSymbol'] in ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {'opls': 'opls_xxx', 'gaff': 'h5'}  # 2 EWD atoms
                            else:
                                hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}  # Default sp2 Hydrogen for c2 and ca types

                        elif carbon_type == {'opls': 'opls_XXX', 'gaff': 'ce'}:
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}
                        elif carbon_type == {'opls': 'opls_521', 'gaff': 'ca'}: 
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}  # Hydrogens connected to C in heterocycle with N as in pyridine C-N-C

                        elif carbon_type == {'opls': 'opls_XXX', 'gaff': 'ce'} or carbon_type == {'opls': 'opls_508', 'gaff': 'cc'} or carbon_type == {'opls': 'opls_XXX', 'gaff': 'cf'} or carbon_type == {'opls': 'opls_XXX', 'gaff': 'cd'}: # Hydrogens connected to a non-pure aromatic cycle or non_aromatic cycle
                            ewd_count = sum(1 for num in info['ConnectedAtomsNumbers'] if self.atom_info_dict[num]['AtomicSymbol'] in ewd_atoms)
                            if ewd_count == 1:
                                hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}  # 1 EWD atom
                            elif ewd_count == 2:
                                hydrogen_type = {'opls': 'opls_h5', 'gaff': 'h5'}  # 2 EWD atom
                            else:
                                hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}  # Default Aliphatic sp2 Hydrogen for c2
                        else:
                            hydrogen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'hx{info["AtomNumber"]}'}

                        self.atom_types_dict[f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

                self.atom_types_dict[f"C{info['AtomNumber']}"] = carbon_type
            
            # Oxygen type decision
            elif info['AtomicSymbol'] == 'O':
                if info.get('CyclicStructure') == 'none':
                    connected_symbols = set(info['ConnectedAtoms'])
                    
                    if info['NumConnectedAtoms'] == 2 and connected_symbols == {'H','H'}:
                        oxygen_type = {'opls': 'opls_154', 'gaff': 'ow'}  # Water
                    
                    elif info['NumConnectedAtoms'] == 2 and 'H' in connected_symbols:
                        oxygen_type = {'opls': 'opls_154', 'gaff': 'oh'}
                    
                    elif info['NumConnectedAtoms'] == 2 and 'C' in connected_symbols:
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'O' for atom in info['ConnectedAtomsNumbers']):
                            oxygen_type = {'opls': 'opls_160', 'gaff': 'os'}  # Ether group or secondary amide
                        else:
                            oxygen_type = {'opls': 'opls_156', 'gaff': 'os'}  # Carbonyl group
                    elif info['NumConnectedAtoms'] == 2:
                        oxygen_type = {'opls': 'opls_156', 'gaff': 'os'} # General case

                    elif info['NumConnectedAtoms'] == 1:
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'C' for atom in info['ConnectedAtomsNumbers']):
                            # This checks if the carbon connected to the oxygen is connected to another oxygen. 
                            # It is useful to identify carboxylic acids and esters.
                            carbons = [atom for atom in info['ConnectedAtomsNumbers'] if self.atom_info_dict[atom]['AtomicSymbol'] == 'C']
                            if any('O' in self.atom_info_dict[carbon]['ConnectedAtoms'] for carbon in carbons):
                                oxygen_type = {'opls': 'opls_157', 'gaff': 'o'}  # Carboxylic acid or ester
                            else:
                                oxygen_type = {'opls': 'opls_158', 'gaff': 'o'}  # Aldehyde
                        elif any(self.atom_info_dict[atom]['AtomicSymbol'] == 'O' for atom in info['ConnectedAtomsNumbers']):
                            oxygen_type = {'opls': 'opls_159', 'gaff': 'oo'}  # Peroxide
                        else:
                            oxygen_type = {'opls': 'opls_157', 'gaff': 'o'}
                    
                    else:
                        oxygen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'ox{info["AtomNumber"]}'}
                    
                    self.atom_types_dict[f"O{info['AtomNumber']}"] = oxygen_type
                    
                    # Hydrogen type assignment based on oxygen type
                    for connected_atom_number in info['ConnectedAtomsNumbers']:
                        connected_atom_info = self.atom_info_dict[connected_atom_number]
                        if connected_atom_info['AtomicSymbol'] == 'H' and connected_atom_info['NumConnectedAtoms'] == 1:
                            if oxygen_type == {'opls': 'opls_154', 'gaff': 'oh'}:
                                hydrogen_type = {'opls': 'opls_240', 'gaff': 'ho'}  # Alcohol Hydrogen
                            else:
                                hydrogen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'hox{info["AtomNumber"]}'}
                            
                            self.atom_types_dict[f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

                if info.get('CyclicStructure') == 'cycle':
                    connected_symbols = set(info['ConnectedAtoms'])
                    
                    if info['NumConnectedAtoms'] == 2 and connected_symbols == {'H'}:
                        oxygen_type = {'opls': 'opls_154', 'gaff': 'oh'}  # Alcohol group
                    
                    elif info['NumConnectedAtoms'] == 2 and 'C' in connected_symbols:
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'O' for atom in info['ConnectedAtomsNumbers']):
                            oxygen_type = {'opls': 'opls_160', 'gaff': 'os'}  # Ether group or secondary amide
                        else:
                            oxygen_type = {'opls': 'opls_156', 'gaff': 'os'}  # Carbonyl group
                    
                    elif info['NumConnectedAtoms'] == 1:
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'C' for atom in info['ConnectedAtomsNumbers']):
                            # This checks if the carbon connected to the oxygen is connected to another oxygen. 
                            # It is useful to identify carboxylic acids and esters.
                            carbons = [atom for atom in info['ConnectedAtomsNumbers'] if self.atom_info_dict[atom]['AtomicSymbol'] == 'C']
                            if any('O' in self.atom_info_dict[carbon]['ConnectedAtoms'] for carbon in carbons):
                                oxygen_type = {'opls': 'opls_157', 'gaff': 'o'}  # Carboxylic acid or ester
                            else:
                                oxygen_type = {'opls': 'opls_158', 'gaff': 'o'}  # Aldehyde
                        elif any(self.atom_info_dict[atom]['AtomicSymbol'] == 'O' for atom in info['ConnectedAtomsNumbers']):
                            oxygen_type = {'opls': 'opls_159', 'gaff': 'oo'}  # Peroxide
                        else:
                            oxygen_type = {'opls': 'opls_157', 'gaff': 'o'}
                    
                    else:
                        oxygen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'ox{info["AtomNumber"]}'}
                    
                    self.atom_types_dict[f"O{info['AtomNumber']}"] = oxygen_type
                    
                    # Hydrogen type assignment based on oxygen type
                    for connected_atom_number in info['ConnectedAtomsNumbers']:
                        connected_atom_info = self.atom_info_dict[connected_atom_number]
                        if connected_atom_info['AtomicSymbol'] == 'H' and connected_atom_info['NumConnectedAtoms'] == 1:
                            if oxygen_type == {'opls': 'opls_154', 'gaff': 'oh'}:
                                hydrogen_type = {'opls': 'opls_240', 'gaff': 'ho'}  # Alcohol Hydrogen
                            else:
                                hydrogen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'hox{info["AtomNumber"]}'}
                            
                            self.atom_types_dict[f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

            # Nitrogen type decision
            elif info['AtomicSymbol'] == 'N':
                connected_symbols = set(info['ConnectedAtoms'])
                connected_atoms = info['ConnectedAtoms']
                connected_atoms_numbers = info['ConnectedAtomsNumbers']

                if info.get('CyclicStructure') == 'none':

                    if info['NumConnectedAtoms'] == 4:

                        num_hydrogens = sum([1 for symbol in connected_atoms if symbol == 'H'])

                        if num_hydrogens == 4:
                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n+'}
                        elif num_hydrogens == 3:
                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'nz'}  # Sp3 N with three hydrogen atoms
                        elif num_hydrogens == 2:
                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'ny'}  # Sp3 N with two hydrogen atoms
                        elif num_hydrogens == 1:
                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'nx'}  # Sp3 N with one hydrogen atom
                        else:
                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n4'}  # Sp3 N with four connected atoms, but no hydrogens

                    elif info['NumConnectedAtoms'] == 3:

                        # Check for Nitro N
                        if connected_symbols == {'C', 'O', 'O'}:
                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'no'}  # Nitro N

                        elif 'C' in connected_symbols:
                            for atom in connected_atoms_numbers:
                                    
                                # Check for aromaticity
                                if 'cycle' in self.atom_info_dict[atom]['CyclicStructure'] and 'pure_aromatic' in self.atom_info_dict[atom]['Aromaticity']:
                                    num_hydrogens = sum([1 for symbol in connected_atoms if symbol == 'H'])
                                    if num_hydrogens == 1:
                                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nu'}  # Like nh, but with 1 attached hydrogen atom
                                    elif num_hydrogens == 2:
                                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nv'}  # Like nh, but with 2 attached hydrogen atoms
                                    else:
                                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nh'}  # Special case of idine
                                    break

                                elif self.atom_info_dict[atom]['NumConnectedAtoms'] == 4: # Connected to an sp3 carbon
                                    num_hydrogens = sum([1 for symbol in connected_atoms if symbol == 'H'])

                                    if num_hydrogens == 1:
                                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n7'}  # Like n3, but with 1 attached hydrogen atom
                                    elif num_hydrogens == 2:
                                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n8'}  # Like n3, but with 2 attached hydrogen atoms
                                    else:
                                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n3'}  # Special case of H2-N-C Sp1
                                    break

                                # Check the atoms connected to the Carbon
                                for connected_to_carbon in self.atom_info_dict[atom]['ConnectedAtomsNumbers']:
                                    atom_symbol = self.atom_info_dict[connected_to_carbon]['AtomicSymbol']
                                    atom_connectivity = self.atom_info_dict[connected_to_carbon]['NumConnectedAtoms']

                                    # Amides and Sulfamides
                                    if (atom_symbol == 'O' and atom_connectivity == 1) or (atom_symbol == 'S' and atom_connectivity == 1):
                                        num_hydrogens = sum([1 for symbol in connected_atoms if symbol == 'H'])
                                        if num_hydrogens == 1:
                                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'ns'}  # Like n, but with 1 attached hydrogen atom
                                        elif num_hydrogens == 2:
                                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'nt'}  # Like n, but with 2 attached hydrogen atoms
                                        else:
                                            nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n'}  # Amide and sulfamide
                                        break
                                    # Idines
                                    elif atom_symbol == 'N' and atom_connectivity == 2:
                                        num_hydrogens = sum([1 for symbol in connected_atoms if symbol == 'H'])
                                        if num_hydrogens == None:
                                            continue
                                        if num_hydrogens == 1:
                                            nitrogen_type = {'opls': 'opls_901', 'gaff': 'nu'}  # Like n3, but with 1 attached hydrogen atom
                                        elif num_hydrogens == 2:
                                            nitrogen_type = {'opls': 'opls_901', 'gaff': 'nv'}  # Like n3, but with 2 attached hydrogen atoms
                                        else:
                                            nitrogen_type = {'opls': 'opls_901', 'gaff': 'nh'}  # Special case of idine
                                        break
                                    # Connected to C sp1
                                    elif atom_symbol == 'C' and atom_connectivity == 2:
                                        # n3 in GAFF but special cases n7 and n8 added in GAFF2
                                        num_hydrogens = sum([1 for symbol in connected_atoms if symbol == 'H'])
                                        if num_hydrogens == None:
                                            continue
                
                                        if num_hydrogens == 1:
                                            nitrogen_type = {'opls': 'opls_300', 'gaff': 'n7'}  # Like n3, but with 1 attached hydrogen atom
                                        elif num_hydrogens == 2:
                                            nitrogen_type = {'opls': 'opls_300', 'gaff': 'n8'}  # Like n3, but with 2 attached hydrogen atoms
                                        else:
                                            nitrogen_type = {'opls': 'opls_300', 'gaff': 'n3'}  # Special case of H2-N-C Sp1
                                        break

                                    elif atom_symbol == 'C' and atom_connectivity == 3:
                                        nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'nh'}  # Special case of H2-N-C Sp2
                                        continue
                                else:
                                    # n3 in GAFF but special cases n7 and n8 added in GAFF2
                                    num_hydrogens = sum([1 for symbol in connected_atoms if symbol == 'H'])

                                    if num_hydrogens == 1:
                                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n7'}  # Like n3, but with 1 attached hydrogen atom
                                    elif num_hydrogens == 2:
                                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n8'}  # Like n3, but with 2 attached hydrogen atoms
                                    else:
                                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n3'}  # Special case of H2-N-C Sp1
                                    break
                                break

                        # General Sp3 N with three connected atoms
                        else:
                            nitrogen_type = {'opls': 'opls_300', 'gaff': 'n3'}

                    elif info['NumConnectedAtoms'] == 2:
                        # Check if the Nitrogen is connected to another Nitrogen with sp1 hybridization
                        if any(self.atom_info_dict[atom]['AtomicSymbol'] == 'N' and self.atom_info_dict[atom]['NumConnectedAtoms'] == 1 for atom in connected_atoms_numbers):
                            nitrogen_type = {'opls': 'opls_753', 'gaff': 'n1'}  # N triple bond
                        else:
                            nitrogen_type = {'opls': 'opls_903', 'gaff': 'n2'}  # General case


                    elif info['NumConnectedAtoms'] == 1:
                        nitrogen_type = {'opls': 'opls_753', 'gaff': 'n1'}  # N triple bond

                    else:
                        nitrogen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'nx{info["AtomNumber"]}'}
            
                # Cyclic structures containing N
                elif info.get('CyclicStructure') == 'cycle':
                    
                    if info['NumConnectedAtoms'] == 2 and 'pure_aromatic' in info.get('Aromaticity'):
                        nitrogen_type = {'opls': 'opls_520', 'gaff': 'nb'}  # Sp2 N in pure aromatic systems 

                    elif info['NumConnectedAtoms'] == 2 and 'non_pure_aromatic' in info.get('Aromaticity'):
                        nitrogen_type = {'opls': 'opls_520', 'gaff': 'nc'}  # Sp2 N in non-pure aromatic systems

                    elif info['NumConnectedAtoms'] == 3 and 'pure_aromatic' in info.get('Aromaticity'):
                        nitrogen_type = {'opls': 'opls_183', 'gaff': 'nb'}  # Pyridine as a ligand in an organometallic complex

                    elif info['NumConnectedAtoms'] == 3 and 'non_pure_aromatic' in info.get('Aromaticity'):
                        nitrogen_type = {'opls': 'opls_na', 'gaff': 'na'}
                
                    # Nitrogens in Non aromatic cycles
                    elif info['NumConnectedAtoms'] == 3 and 'non_aromatic' in info.get('Aromaticity'):
                        if 3 in info.get('CycleSize'):
                            nitrogen_type = {'opls': 'opls_np', 'gaff': 'np'}
                        if 4 in info.get('CycleSize'):
                            nitrogen_type = {'opls': 'opls_nq', 'gaff': 'nq'}
                        else:
                            nitrogen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'nx{info["AtomNumber"]}'}
                    else:
                        # Add other conditions for cyclic nitrogen atoms if needed or add a default case
                        nitrogen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'nx{info["AtomNumber"]}'}

                else:
                    nitrogen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'nx{info["AtomNumber"]}'}   
                
                self.atom_types_dict[f"N{info['AtomNumber']}"] = nitrogen_type

                # Hydrogen type assignment based on nitrogen type
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    connected_atom_info = self.atom_info_dict[connected_atom_number]
                    # Assign hydrogen type
                    if connected_atom_info['AtomicSymbol'] == 'H' and connected_atom_info['NumConnectedAtoms'] == 1:
                        hydrogen_type = {'opls': 'opls_240', 'gaff': 'hn'}

                        self.atom_types_dict[f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

            # Phosphorus type decision
            elif info['AtomicSymbol'] == 'P':
                if info.get('CyclicStructure') == 'none':
                    connected_symbols = set(info['ConnectedAtoms'])
                    
                    # Phosphate groups or phosphoric acid
                    if info['NumConnectedAtoms'] == 4 and 'O' in connected_symbols:
                        oxygen_count = list(connected_symbols).count('O')
                        if oxygen_count == 4:
                            phosphorus_type = {'opls': 'opls_900P', 'gaff': 'p5'}  # Simplified, it could be tetrahedral phosphate
                        else:
                            phosphorus_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'px{info["AtomNumber"]}'}
                    
                    # Phosphine
                    elif info['NumConnectedAtoms'] == 3 and 'H' in connected_symbols:
                        hydrogen_count = list(connected_symbols).count('H')
                        if hydrogen_count == 3:
                            phosphorus_type = {'opls': 'opls_901P', 'gaff': 'ph3'}
                        else:
                            phosphorus_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'px{info["AtomNumber"]}'}
                    
                    # Phosphine oxides
                    elif info['NumConnectedAtoms'] == 4 and 'O' in connected_symbols:
                        hydrogen_count = list(connected_symbols).count('H')
                        if hydrogen_count == 3:
                            phosphorus_type = {'opls': 'opls_902P', 'gaff': 'po'}
                        else:
                            phosphorus_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'px{info["AtomNumber"]}'}
                    
                    # Phosphonates and Phosphites
                    elif info['NumConnectedAtoms'] == 3 and 'O' in connected_symbols:
                        oxygen_count = list(connected_symbols).count('O')
                        if oxygen_count == 2:
                            phosphorus_type = {'opls': 'opls_903P', 'gaff': 'p3'}  # Again simplified, could distinguish between phosphonates and phosphites
                        else:
                            phosphorus_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'px{info["AtomNumber"]}'}
                    
                    else:
                        phosphorus_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'px{info["AtomNumber"]}'}
                    
                    self.atom_types_dict[f"P{info['AtomNumber']}"] = phosphorus_type
                    
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
                    elif all(self.atom_info_dict[num]['AtomicSymbol'] in ['C', 'N', 'S'] for num in info['ConnectedAtomsNumbers']):  # Both connected atoms are carbons or one carbon and one nitrogen
                        sulfur_type = {'opls': 'opls_SS', 'gaff': 'ss'}  # Thio-ether or Thio-ester
                    else:
                        sulfur_type = {'opls': 'opls_921S', 'gaff': 's2'}

                # S with three connected atoms 
                elif info['NumConnectedAtoms'] == 3:
                    if any(self.atom_info_dict[num]['AtomicSymbol'] == 'C' and self.atom_info_dict[num]['NumConnectedAtoms'] == 3 for num in info['ConnectedAtomsNumbers']):
                        sulfur_type = {'opls': 'opls_922X', 'gaff': 'sx'}
                    else:
                        sulfur_type = {'opls': 'opls_922S', 'gaff': 's4'}

                # S with four connected atoms 
                elif info['NumConnectedAtoms'] == 4:
                    if any(self.atom_info_dict[num]['AtomicSymbol'] == 'C' and self.atom_info_dict[num]['NumConnectedAtoms'] == 3 for num in info['ConnectedAtomsNumbers']):
                        sulfur_type = {'opls': 'opls_922X', 'gaff': 'sy'}
                    else:
                        sulfur_type = {'opls': 'opls_923S', 'gaff': 's6'}

                # Sp3 S connected with hydrogen 
                elif info['NumConnectedAtoms'] == 4 and {'H', 'H', 'H'} <= connected_symbols:  # use <= to check if a set is a subset
                    sulfur_type = {'opls': 'opls_924S', 'gaff': 'sh'}

                else:
                    sulfur_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'sx{info["AtomNumber"]}'}

                self.atom_types_dict[f"S{info['AtomNumber']}"] = sulfur_type

                    
                # Hydrogen assignment in the case of thiols
                for connected_atom_number in info['ConnectedAtomsNumbers']:
                    connected_atom_info = self.atom_info_dict[connected_atom_number]
                    if connected_atom_info['AtomicSymbol'] == 'H' and connected_atom_info['NumConnectedAtoms'] == 1:
                        if sulfur_type == {'opls': 'opls_924S', 'gaff': 'sh'}:
                            hydrogen_type = {'opls': 'opls_926H', 'gaff': 'hs'}
                        else:
                            hydrogen_type = {'opls': f'opls_x{info["AtomNumber"]}', 'gaff': f'hsx{info["AtomNumber"]}'}

                        self.atom_types_dict[f"H{connected_atom_info['AtomNumber']}"] = hydrogen_type

            # Decision for halogens

            elif info['AtomicSymbol'] == 'Br':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'br'}
                self.atom_types_dict[f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type

            elif info['AtomicSymbol'] == 'Cl':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'cl'}
                self.atom_types_dict[f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type

            elif info['AtomicSymbol'] == 'F':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'f'}
                self.atom_types_dict[f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type

            elif info['AtomicSymbol'] == 'I':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'i'}
                self.atom_types_dict[f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type


            # Decision for Transition Metals

            elif info['AtomicSymbol'] not in ['C', 'H', 'O', 'N', 'P', 'Br', 'Cl', 'F', 'I']: 
                
                # Assign atom types based on AtomicSymbol and AtomNumber
                atom_type = {'opls': f'opls_{info["AtomicSymbol"]}{info["AtomNumber"]}', 
                            'gaff': f'{info["AtomicSymbol"]}{info["AtomNumber"]}'}
                self.atom_types_dict[f"{info['AtomicSymbol']}{info['AtomNumber']}"] = atom_type

            else:
                # Else for the cases falling off the decision tree 
                # The Hydrogen are assigned outside the main branches of the decision tree
                # Therefore, they need to be out of the else case.

                if info['AtomicSymbol'] != 'H':
                    print('Warning:',f"{info['AtomicSymbol']}{info['AtomNumber']}", 'Has not been found in the decision tree, check it carefully')


        return self.atom_types_dict

    def extract_gaff_atom_types(self, atom_type):
        """
        Extracts GAFF atom types from the atom types dictionary.

        This method sorts the atom types based on their numerical suffix and then extracts
        the GAFF atom types from the dictionary that contains atom type information for
        both OPLS and GAFF force fields.

        Returns:
            list: A list of GAFF atom types in the order determined by the sorted atom keys.
        """

        self.gaff_atom_types = []

        # Sort atom types based on the number after the atomic symbol
        sorted_atom_types = sorted(self.atom_types_dict.keys(), key=self.get_atom_number)

        for atom_type in sorted_atom_types:
            if isinstance(self.atom_types_dict[atom_type], dict):
                gaff_type = self.atom_types_dict[atom_type].get('gaff', None)
                if gaff_type:
                    self.gaff_atom_types.append(gaff_type)

        return self.gaff_atom_types


    def generate_gaff_atomtypes(self, molecule):
        """
        Generates GAFF (General Amber Force Field) atom types for a given molecule.

        This method takes a VeloxChem molecule object as input and performs several
        steps to determine GAFF atom types for each atom in the molecule.

        Steps:
        1. Extracts molecular information such as coordinates, atomic symbols, and
        the number of atoms from the input molecule.
        2. Calculates covalent radii for each atom and stores them.
        3. Creates a connectivity matrix and computes pairwise distances between atoms.
        4. Detects closed cyclic structures within the molecule and gathers information
        about cyclic atoms, cycle sizes, aromaticity, and cycle IDs.
        5. Generates an atom information dictionary using the gathered data.
        6. Determines atom types for each atom based on the atom information dictionary.
        7. Extracts GAFF atom types from the determined atom types.

        Args:
            molecule: A VeloxChem molecule object.

        Returns:
            list: A list of GAFF atom types for each atom in the molecule.

        Note:
            The method populates various attributes of the class instance with
            intermediate data during its execution, such as 'coordinates', 'atomic_symbols',
            'covalent_radii', 'connectivity_matrix', 'distances', 'cyclic_atoms',
            'cycle_sizes', 'aromaticity', 'cycles', 'cycle_ids', 'atom_cycle_info',
            'atom_info_dict', and 'atom_types_dict'.

        Example:
            To generate GAFF atom types for a VeloxChem molecule 'my_molecule', you can
            call this method as follows:
            gaff_atom_types = my_instance.generate_gaff_atomtypes(my_molecule)
        """

        # Workflow of the method
        self.coordinates = molecule.get_coordinates_in_angstrom()
        self.atomic_symbols = molecule.get_labels()
        self.num_atoms = len(self.atomic_symbols)
        self.covalent_radii = molecule.covalent_radii_to_numpy() * bohr_in_angstrom()
        self.connectivity_matrix, self.distances = self.create_connectivity_matrix(self.coordinates, self.covalent_radii)
        self.cyclic_atoms, self.cycle_sizes, self.aromaticity, self.cycles, self.cycle_ids, self.atom_cycle_info = self.detect_closed_cyclic_structures(self.atomic_symbols, self.connectivity_matrix, self.distances)
        self.atom_info_dict = self.create_atom_info_dict(self.atomic_symbols, self.connectivity_matrix, self.distances, self.cyclic_atoms, self.cycle_sizes, self.aromaticity, self.cycles, self.cycle_ids, self.atom_cycle_info)
        self.atom_types_dict = self.decide_atom_type(self.atom_info_dict)
        self.gaff_atom_types = self.extract_gaff_atom_types(self.atom_types_dict)

        # Printing output
        print("VeloxChem Atom Type Identification\n")
        print("-" * 40)  # Dashed line

        # Detected number of atoms
        num_atoms = len(self.atomic_symbols)
        print(f"Detected number of atoms: {num_atoms}\n")

        # Print table header
        print("{:<30} {:<20}".format("Symbol (id)", "GAFF atom type assigned"))

        # Print atom symbol, atom number, and GAFF atom type for each atom
        for i, (symbol, gaff_type) in enumerate(zip(self.atomic_symbols, self.gaff_atom_types), start=1):
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

        return self.gaff_atom_types
        
    def compute_structural_features(self):
        '''
        Computes the structural features of a molecule based on its connectivity and distance matrices. The features include bonds, angles, dihedrals, and improper dihedrals. Each feature is represented by the involved atoms' IDs, the corresponding atom types, and the geometric parameters like distances and angles.

        This method populates a dictionary with the following keys and values:
        - 'bonds': A list of dictionaries, each representing a bond. Each dictionary contains:
        - 'atoms': A tuple of atom IDs that form the bond.
        - 'types': A tuple of GAFF atom types for the atoms in the bond.
        - 'distance': The distance between the atoms, derived from the distance matrix.
        
        - 'angles': A list of dictionaries, each representing an angle formed by three connected atoms. Each dictionary contains:
        - 'atoms': A tuple of atom IDs that form the angle.
        - 'types': A tuple of GAFF atom types for the atoms in the angle.
        - 'angle': The calculated angle between the atoms, in degrees.

        - 'dihedrals': A list of dictionaries, each representing a dihedral angle formed by four sequentially bonded atoms. Each dictionary contains:
        - 'atoms': A tuple of atom IDs that form the dihedral.
        - 'types': A tuple of GAFF atom types for the atoms in the dihedral.
        - 'dihedral': The calculated dihedral angle between the atoms, in degrees.
        
        - 'impropers': A list of dictionaries, each representing an improper dihedral angle where one atom is connected to three others, forming a pyramid structure. Each dictionary contains:
        - 'atoms': A tuple of atom IDs, with the first being the central atom.
        - 'types': A tuple of GAFF atom types for the atoms in the improper dihedral.
        - 'improper_angle': The calculated improper dihedral angle, in degrees.

        - 'pairs': Currently initialized as an empty list, to be populated with non-bonded atom pairs that need to be considered during force field generation.

        Returns:
        - A dictionary containing the lists of bonds, angles, dihedrals, impropers, and pairs with their respective structural information.

        Note:
        - This method assumes that `self.connectivity_matrix` and `self.distance_matrix` are already computed and available as attributes of the class instance.
        - Atom IDs are assumed to be 1-indexed, corresponding to their order in the `self.gaff_atom_types` list.
        - The method `AtomTypeIdentifier.calculate_angle` is used to calculate the angles, and `AtomTypeIdentifier.calculate_dihedral` is used for dihedral angles. These methods should be defined elsewhere in the AtomTypeIdentifier class.
        - The `self.coordinates` attribute is assumed to hold the coordinates of each atom, used for geometric calculations.
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
        for i in range(len(self.connectivity_matrix)):
            for j in range(i + 1, len(self.connectivity_matrix)):
                if self.connectivity_matrix[i][j] == 1:
                    bond = {
                        'atoms': (i + 1, j + 1),
                        'types': (self.gaff_atom_types[i], self.gaff_atom_types[j]),
                        'distance': self.distance_matrix[i][j]
                    }
                    structural_features['bonds'].append(bond)

        # Compute angles based on the connectivity matrix and the distance matrix
        for i in range(len(self.connectivity_matrix)):
            for j in range(len(self.connectivity_matrix)):
                if self.connectivity_matrix[i][j] == 1:  # If i and j are connected
                    for k in range(j + 1, len(self.connectivity_matrix)):
                        if self.connectivity_matrix[j][k] == 1 and self.connectivity_matrix[i][k] == 0:  # If j and k are connected, but i and k are not (forming an angle)
                            angle = {
                                'atoms': (i + 1, j + 1, k + 1),
                                'types': (self.gaff_atom_types[i], self.gaff_atom_types[j], self.gaff_atom_types[k]),
                                'angle': AtomTypeIdentifier.calculate_angle(
                                    self.coordinates[i], self.coordinates[j], self.coordinates[k]
                                )
                            }
                            structural_features['angles'].append(angle)

        # Compute dihedrals based on the connectivity matrix and the distance matrix
        for i in range(len(self.connectivity_matrix)):
            for j in range(len(self.connectivity_matrix)):
                if self.connectivity_matrix[i][j] == 1:  # If i and j are connected
                    for k in range(len(self.connectivity_matrix)):
                        if self.connectivity_matrix[j][k] == 1 and i != k:  # If j and k are connected, and i and k are different
                            for l in range(len(self.connectivity_matrix)):
                                # Ensure k and l are connected, l is not j, and i and l are not directly connected (to form a chain)
                                if self.connectivity_matrix[k][l] == 1 and l != j and self.connectivity_matrix[i][l] == 0:
                                    dihedral = {
                                        'atoms': (i + 1, j + 1, k + 1, l + 1),
                                        'types': (
                                            self.gaff_atom_types[i], self.gaff_atom_types[j],
                                            self.gaff_atom_types[k], self.gaff_atom_types[l]
                                        ),
                                        'dihedral': AtomTypeIdentifier.calculate_dihedral(
                                            self.coordinates[i], self.coordinates[j],
                                            self.coordinates[k], self.coordinates[l]
                                        )
                                    }
                                    structural_features['dihedrals'].append(dihedral)

        # Compute improper dihedrals based on the connectivity matrix
        for i in range(len(self.connectivity_matrix)):
            # Find the central atom bonded to three other atoms
            bonded_atoms = [j for j in range(len(self.connectivity_matrix)) if self.connectivity_matrix[i][j] == 1]
            if len(bonded_atoms) == 3:
                # Indices of the three bonded atoms
                j, k, l = bonded_atoms
                improper_dihedral = {
                    'atoms': (i+1, j+1, k+1, l+1),
                    'types': (self.gaff_atom_types[i], self.gaff_atom_types[j], self.gaff_atom_types[k], self.gaff_atom_types[l]),
                    'improper_angle': AtomTypeIdentifier.calculate_improper_dihedral(
                        self.coordinates[i], self.coordinates[j], self.coordinates[k], self.coordinates[l]
                    )
                }
                structural_features['impropers'].append(improper_dihedral)

        return structural_features
    
    def generate_force_field(self, ff_file_path):
        '''
        Generates a force field dictionary for a molecule based on its structural features and a given GAFF force field file.

        Arguments:
        - ff_file_path (str): The file path to the GAFF force field data file.

        Returns:
        - dict: A dictionary containing force field parameters for atom types, bonds, angles, dihedrals, impropers, and non-bonded pairs.

        Overview:
        This method constructs a force field data structure that includes parameters for each atom type, bond, angle, dihedral, and improper torsion in the molecule. It also calculates non-bonded pairs that are not explicitly covered by bonds or angles. The method relies on pre-computed structural features such as connectivity and distance matrices provided by `self.compute_structural_features()`.

        Details:
        - The force field data for each atom type includes its mass, sigma, epsilon values, and any additional info as defined in the GAFF force field file.
        - Bonds are defined by the atoms they connect and include force constants and equilibrium distances.
        - Angles are determined by three connected atoms and include force constants and equilibrium angles.
        - Dihedrals are defined by four sequentially bonded atoms and include force constants, equilibrium angles, periodicity, and phase.
        - Improper torsions are defined by a central atom connected to three other atoms, including force constants and equilibrium angles.
        - Non-bonded pairs are computed based on dihedrals but exclude atoms already connected by bonds or as terminal atoms in angles.

        The output dictionary keys are as follows:
        - 'atomtypes': Each entry contains the ID, type, mass, sigma, epsilon, and additional info for each atom type.
        - 'bonds': A list of bond dictionaries, each with atom IDs, types, force constants, and equilibrium distances.
        - 'angles': A list of angle dictionaries, each with atom IDs, types, force constants, and equilibrium angles.
        - 'dihedrals': A list of dihedral dictionaries, each with atom IDs, types, force constants, equilibrium angles, periodicity, and phase.
        - 'impropers': A list of improper dictionaries, each with atom IDs, types, force constants, and equilibrium angles.
        - 'pairs': A dictionary of non-bonded atom pairs with comments describing their computed origin.

        Exceptions:
        - FileNotFoundError: If the provided force field file path does not exist or is inaccessible.
        - ValueError or IndexError: If there are issues parsing the force field file or if the structural features do not align with the expected format.

        Usage:
        The method is designed to be used within an instance of a molecule class that contains methods and properties for structural feature computation and has GAFF atom types assigned to its atoms.

        Example:
        ```python
        molecule = MoleculeClass(...)  # Instance of a molecule class with GAFF atom types and structural features.
        force_field_data = molecule.generate_force_field('./gaff-2.20.dat')
        '''
        
        # Required structural features
        self.structural_features = self.compute_structural_features()

        # Initialize the dictionary to hold force field data
        force_field_data = {
            'atomtypes': {},
            'bonds': {},
            'angles': {},
            'dihedrals':{},
            'impropers':{},
            'pairs':{}
        }

        # Conversion factors
        angstrom_to_nm = 0.1
        kcalmol_to_kjmol = 4.184

        # Read the force field file
        with open(ff_file_path, 'r') as ff_file:
            ff_lines = ff_file.readlines()

        # Create a dictionary from the force field file for atom types
        atomtype_data = {line.split()[0]: line for line in ff_lines if re.match(r'^\S+\s+\d+\.\d+\s+\d+\.\d+\s+.*$', line)}

        # Create a dictionary from the force field file for sigma and epsilon values
        sigma_epsilon_data = {line.split()[0]: line for line in ff_lines if re.match(r'^\s*\S+\s+\d+\.\d+\s+\d+\.\d+', line)}

        # Parse atomtypes section and initialize force field data with unique IDs
        for index, gaff_atom_type in enumerate(self.gaff_atom_types, start=1):  # Start index at 1 since IDs start at 1
            force_field_data['atomtypes'][index] = {
                'type': gaff_atom_type,
                'id': index,
                'mass': 0,
                'sigma': 0,
                'epsilon': 0,
                'info': 'undefined'
            }

            # If atom type data is found in the force field, update mass and info
            if gaff_atom_type in atomtype_data:
                atom_type_match = re.match(r'^(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(.+)$', atomtype_data[gaff_atom_type])
                if atom_type_match:
                    mass = float(atom_type_match.group(2))
                    info = atom_type_match.group(4).strip()
                    force_field_data['atomtypes'][index].update({
                        'mass': mass,
                        'info': info
                    })

            # If sigma and epsilon data is found in the force field, update those values
            if gaff_atom_type in sigma_epsilon_data:
                sigma_epsilon_match = re.match(r'^\s*(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)', sigma_epsilon_data[gaff_atom_type])
                if sigma_epsilon_match:
                    sigma = float(sigma_epsilon_match.group(2)) * angstrom_to_nm
                    epsilon = float(sigma_epsilon_match.group(3)) * kcalmol_to_kjmol
                    force_field_data['atomtypes'][index].update({
                        'sigma': sigma,
                        'epsilon': epsilon
                    })

            # Parse bonds section
            for bond in self.structural_features['bonds']:
                bond_type_pattern = '-'.join(sorted(bond['types']))
                bond_regex = re.escape(bond_type_pattern) + r'\s+(\d+\.\d+)\s+(\d+\.\d+)'
                match_found = False  # Flag to check if we find a match in the force field file

                for line in ff_lines:
                    bond_match = re.search(bond_regex, line)
                    if bond_match:
                        force_constant = float(bond_match.group(1)) * kcalmol_to_kjmol
                        eq_distance = float(bond_match.group(2)) * angstrom_to_nm
                        # The bond key will include atom types and IDs for uniqueness
                        bond_key = tuple(sorted((bond['atoms'][0], bond['atoms'][1])))
                        force_field_data['bonds'][bond_key] = {
                            'ids': bond['atoms'],  # Include the atom IDs
                            'types': bond['types'],
                            'force_constant': force_constant,
                            'eq_distance': eq_distance,
                            'comment': 'GAFF2'
                        }
                        match_found = True
                        break  # Exit the loop after finding the match

                if not match_found:
                    # Default values if no match is found, using the computed distance from structural_features
                    computed_distance = bond['distance']  # Assuming 'distance' is already in the correct units
                    # The bond key will include atom types and IDs for uniqueness
                    bond_key = tuple(sorted((bond['atoms'][0], bond['atoms'][1])))
                    force_field_data['bonds'][bond_key] = {
                        'ids': bond['atoms'],  # Include the atom IDs
                        'types': bond['types'],
                        'force_constant': 250000.000 * kcalmol_to_kjmol,  # Default force constant converted to kJ/mol
                        'eq_distance': computed_distance,
                        'comment': 'unknown'
                    }


        # Parse angles section
        for angle in self.structural_features['angles']:
            angle_type_pattern = '-'.join(sorted(angle['types']))
            angle_regex = re.escape(angle_type_pattern) + r'\s+(\d+\.\d+)\s+(\d+\.\d+)'
            match_found = False  # Flag to check if we find a match in the force field file

            for line in ff_lines:
                angle_match = re.search(angle_regex, line)
                if angle_match:
                    force_constant = float(angle_match.group(1)) * kcalmol_to_kjmol
                    eq_angle_radians = float(angle_match.group(2))
                    eq_angle_degrees = eq_angle_radians  # Assuming angle is in radians, convert to degrees if necessary
                    # The angle key will include atom types and IDs for uniqueness
                    angle_key = tuple(sorted((angle['atoms'][0], angle['atoms'][1], angle['atoms'][2])))
                    force_field_data['angles'][angle_key] = {
                        'ids': angle['atoms'],  # Include the atom IDs
                        'types': angle['types'],
                        'force_constant': force_constant,
                        'eq_angle': eq_angle_degrees,
                        'comment': 'GAFF2'
                    }
                    match_found = True
                    break  # Exit the loop after finding the match

            if not match_found:
                # Default values if no match is found, using the computed angle from structural_features
                computed_angle = angle['angle']  # Assuming 'angle' is already computed and stored in degrees
                # The angle key will include atom types and IDs for uniqueness
                angle_key = tuple(sorted((angle['atoms'][0], angle['atoms'][1], angle['atoms'][2])))
                force_field_data['angles'][angle_key] = {
                    'ids': angle['atoms'],  # Include the atom IDs
                    'types': angle['types'],
                    'force_constant': 2000.000 * kcalmol_to_kjmol,  # Default force constant converted to kJ/mol
                    'eq_angle': computed_angle,
                    'comment': 'unknown'
                }

        
        # Parse dihedrals section
        for dihedral in self.structural_features['dihedrals']:
            # Use 'X' as a wildcard for the terminal atoms in the dihedral pattern
            dihedral_type_pattern = 'X-' + '-'.join(dihedral['types'][1:3]) + '-X'
            dihedral_regex = re.escape(dihedral_type_pattern) + r'\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+)'
            match_found = False  # Flag to check if we find a match in the force field file

            for line in ff_lines:
                dihedral_match = re.search(dihedral_regex, line)
                if dihedral_match:
                    periodicity = int(dihedral_match.group(1))
                    force_constant = float(dihedral_match.group(2)) * kcalmol_to_kjmol
                    eq_angle = float(dihedral_match.group(3))  # Assuming this angle does not require conversion
                    phase = float(dihedral_match.group(4))
                    # The dihedral key will include atom types and IDs for uniqueness
                    dihedral_key = tuple(sorted((dihedral['atoms'][0], dihedral['atoms'][1], dihedral['atoms'][2], dihedral['atoms'][3])))
                    force_field_data['dihedrals'][dihedral_key] = {
                        'ids': dihedral['atoms'],  # Include the atom IDs
                        'types': dihedral['types'],
                        'periodicity': periodicity,
                        'force_constant': force_constant,
                        'eq_angle': eq_angle,
                        'phase': phase,
                        'comment': 'GAFF2'
                    }
                    match_found = True
                    break  # Exit the loop after finding the match

            if not match_found:
                # Default values if no match is found, using the computed angle from structural_features
                computed_angle = dihedral.get('angle', 0)  # Default to 0 if not provided
                # The dihedral key will include atom types and IDs for uniqueness
                dihedral_key = tuple(sorted((dihedral['atoms'][0], dihedral['atoms'][1], dihedral['atoms'][2], dihedral['atoms'][3])))
                force_field_data['dihedrals'][dihedral_key] = {
                    'ids': dihedral['atoms'],  # Include the atom IDs
                    'types': dihedral['types'],
                    'periodicity': 2,  # Default periodicity
                    'force_constant': 0,  # Default force constant
                    'eq_angle': computed_angle,
                    'phase': 0,  # Default phase
                    'comment': 'unknown'
                }


        # Parse impropers section
        for improper in self.structural_features['impropers']:
            # Assuming 'types' and 'ids' are lists of the four atom types and IDs involved in the improper
            improper_type_pattern = '-'.join(improper['types'])
            improper_regex = re.escape(improper_type_pattern) + r'\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
            match_found = False  # Flag to check if we find a match in the force field file

            for line in ff_lines:
                improper_match = re.search(improper_regex, line)
                if improper_match:
                    periodicity = int(improper_match.group(1))
                    force_constant = float(improper_match.group(2)) * kcalmol_to_kjmol
                    eq_angle = float(improper_match.group(3))  # Assuming this angle does not require conversion
                    # The improper key will include atom types and IDs for uniqueness
                    improper_key = tuple(improper['ids'])  # Use atom IDs instead of types for the key
                    force_field_data['impropers'][improper_key] = {
                        'ids': improper['ids'],  # Include the atom IDs
                        'types': improper['types'],
                        'periodicity': periodicity,
                        'force_constant': force_constant,
                        'eq_angle': eq_angle,
                        'comment': 'GAFF2'
                    }
                    match_found = True
                    break  # Exit the loop after finding the match

            if not match_found:
                # Default values if no match is found, using the computed angle from structural_features
                computed_angle = improper.get('angle', 0)  # Default to 0 if not provided
                # The improper key will include atom types and IDs for uniqueness
                improper_key = tuple(improper['ids'])  # Use atom IDs instead of types for the key
                force_field_data['impropers'][improper_key] = {
                    'ids': improper['ids'],  # Include the atom IDs
                    'types': improper['types'],
                    'periodicity': 2,  # Default periodicity
                    'force_constant': 0,  # Default force constant
                    'eq_angle': computed_angle,
                    'comment': 'unknown'
                }

        # Add all dihedral pairs
        for dihedral in self.structural_features['dihedrals']:
            # Extract first and last atoms from the dihedral by IDs
            dihedral_pair_ids = (dihedral['atoms'][0], dihedral['atoms'][-1])
            # Initialize the pair in the force field data
            force_field_data['pairs'][dihedral_pair_ids] = {'comment': 'computed from dihedrals'}

        # Exclude pairs that are already present in bonds
        for bond in self.structural_features['bonds']:
            bond_pair_ids = (bond['atoms'][0], bond['atoms'][1])
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

        Args:
            atom_type_str (str): The atom type string containing a numeric suffix.

        Returns:
            int: The numeric part extracted from the atom type string. Returns 0 if no number is found.
        """
        match = re.search(r'\d+', atom_type_str)
        return int(match.group()) if match else 0
    
    @staticmethod
    def calculate_3d_distance(coord1, coord2):
        return np.linalg.norm(np.array(coord1) - np.array(coord2))

    @staticmethod
    def calculate_angle(coord1, coord2, coord3):
        """
        Calculate the angle formed by three points given their coordinates.

        Args:
            coord1 (np.array): The coordinates of the first point.
            coord2 (np.array): The coordinates of the second point (vertex of the angle).
            coord3 (np.array): The coordinates of the third point.

        Returns:
            float: The angle in degrees.
        """
        # Vectors from coord2 to coord1 and coord3
        vector1 = coord1 - coord2
        vector2 = coord3 - coord2

        # Compute the dot product and the magnitudes of the vectors
        dot_product = np.dot(vector1, vector2)
        mag_vector1 = np.linalg.norm(vector1)
        mag_vector2 = np.linalg.norm(vector2)

        # Ensure the cosine value falls within the valid range of -1 to 1
        cosine_angle = np.clip(dot_product / (mag_vector1 * mag_vector2), -1.0, 1.0)

        # Calculate the angle using the law of cosines
        angle_rad = np.arccos(cosine_angle)

        # Convert the angle from radians to degrees
        angle_deg = np.degrees(angle_rad)

        return angle_deg

    @staticmethod
    def calculate_dihedral(point1, point2, point3, point4):
        """
        Calculate the dihedral angle between four points.

        Args:
            point1 (array_like): The xyz coordinates of the first atom.
            point2 (array_like): The xyz coordinates of the second atom.
            point3 (array_like): The xyz coordinates of the third atom.
            point4 (array_like): The xyz coordinates of the fourth atom.

        Returns:
            float: The dihedral angle in degrees.
        """
        # Convert points to numpy arrays if they aren't already
        p0 = np.asarray(point1)
        p1 = np.asarray(point2)
        p2 = np.asarray(point3)
        p3 = np.asarray(point4)

        # Vectors between points
        b0 = p0 - p1
        b1 = p1 - p2
        b2 = p2 - p3

        # Normal vectors to planes formed by the vectors
        n1 = np.cross(b0, b1)
        n2 = np.cross(b1, b2)
        
        # Normalize the normal vectors
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)

        # Unit vector along b1
        m1 = np.cross(n1, b1/np.linalg.norm(b1))

        # Calculate the dihedral angle
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)

        # Compute the angle using arctan2 to get the correct quadrant
        # TODO double check if we have made sure that we don't get e.g. NaN
        angle_rad = np.arctan2(y, x)
        angle_deg = np.degrees(angle_rad)

        return angle_deg

    @staticmethod
    def calculate_improper_dihedral(point1, point2, point3, point4):
        """
        Calculate the improper dihedral angle for a given set of four points.

        Args:
            point1 (array_like): The xyz coordinates of the central atom.
            point2 (array_like): The xyz coordinates of the first outer atom.
            point3 (array_like): The xyz coordinates of the second outer atom.
            point4 (array_like): The xyz coordinates of the third outer atom.

        Returns:
            float: The improper dihedral angle in degrees.
        """
        # Convert points to numpy arrays if they aren't already
        p0 = np.asarray(point1)
        p1 = np.asarray(point2)
        p2 = np.asarray(point3)
        p3 = np.asarray(point4)

        # Vectors from central atom to outer atoms
        v1 = p1 - p0
        v2 = p2 - p0
        v3 = p3 - p0

        # Normal vector of the plane formed by the first three atoms
        n1 = np.cross(v1, v2)

        # Normal vector of the plane formed by the last three atoms
        n2 = np.cross(v2, v3)

        # Normalize the normal vectors
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)

        # Calculate the angle between the normals
        cos_theta = np.dot(n1, n2)
        sin_theta = np.linalg.norm(np.cross(n1, n2))

        # Compute the angle using arctan2 to get the correct quadrant
        # TODO double check if we have made sure that we don't get e.g. NaN
        angle_rad = np.arctan2(sin_theta, cos_theta)
        angle_deg = np.degrees(angle_rad)

        # The improper dihedral angle is the complement of the angle between normals
        improper_angle_deg = 180.0 - angle_deg

        return improper_angle_deg

