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

from mpi4py import MPI
import numpy as np
import networkx as nx
import sys
import re

from .veloxchemlib import mpi_master, bohr_in_angstrom
from .outputstream import OutputStream


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
        - molecule: A VeloxChem molecule object.
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
        - cycle_sizes: A list of the sizes of each cyclic structure.
        - aromaticity: A list of aromaticity classifications for each cyclic
          structure.
        - cycles: A list of cyclic structures.
        - cycle_ids: A list of cyclic structure IDs.
        - atom_cycle_info: A dictionary containing cycle information for each
          atom.
        - atom_info_dict: A dictionary containing detailed information for each
          atom.
        - atom_types_dict: A dictionary containing the GAFF atom types for each
          atom.

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

    def create_connectivity_matrix(self, factor=1.3):
        """
        Creates a connectivity matrix for the molecule based on the atomic
        coordinates and covalent radii, determining which atoms are bonded.

        This method iterates through pairs of atoms and calculates the distance
        between them. If this distance is less than or equal to the sum of
        their covalent radii scaled by a factor, it is considered a bond, and
        the connectivity matrix is updated accordingly. The method also
        constructs a corresponding distance matrix with the actual distances
        between connected atoms.

        :param factor:
            A scaling factor for the covalent radii to account for the bond
            threshold.  Default value is 1.3.

        :return:
            tuple: A tuple containing two 2D numpy arrays:
                   - The first array is the connectivity matrix with 1s
                     indicating bonded atom pairs.
                   - The second array is the distance matrix with actual
                     distances between atoms.
        """

        num_atoms = len(self.coordinates)
        self.connectivity_matrix = np.zeros((num_atoms, num_atoms), dtype='int32')
        self.distance_matrix = np.zeros((num_atoms, num_atoms))

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                distance = self.measure_length(self.coordinates[i],
                                               self.coordinates[j])
                adjusted_threshold = (self.covalent_radii[i] +
                                      self.covalent_radii[j]) * float(factor)

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

    def is_sp2_carbon(self, atom_idx):
        """
        Determines if a C atom, identified by its index, is sp2 hybridized.

        :param atom_idx:
            Index of the atom in the molecule.

        :return:
            True if the atom is sp2 hybridized, False otherwise.
        """

        return (self.atomic_symbols[atom_idx] == "C" and
                len(list(self.graph.neighbors(atom_idx))) == 3)

    def is_sp2_nitrogen(self, atom_idx):
        """
        Determines if a N atom, identified by its index, is sp2 hybridized.

        :param atom_idx:
            Index of the atom in the molecule.

        :return:
            True if the atom is sp2 hybridized, False otherwise.
        """

        return (self.atomic_symbols[atom_idx] == "N" and
                len(list(self.graph.neighbors(atom_idx))) == 2)

    def detect_closed_cyclic_structures(self):
        """
        Detects closed cyclic structures in a molecule and determines their
        aromaticity.

        This method analyzes the graph of atoms and their connectivity to
        identify cycles, determine the size of each cycle, and classify them
        based on aromaticity criteria.
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

        # Convert cycles_to_remove to list of lists and then to a set of tuples
        # for faster lookup
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
                    if max_distance - min_distance <= 0.08:
                        aro = "pure_aromatic"
                    else:
                        aro = "non_pure_aromatic"
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
                    # Check if the S is connected to 2 atoms in the cycle, if
                    # so, it is non_pure_aromatic
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

        # Additional logic for reassignment of aromaticity in special cases
        # where 3 atoms are shared with aromatic rings.
        for index, cycle in enumerate(self.reduced_cycles):
            # Check if all carbons in the cycle are sp2
            all_carbons_sp2 = all(
                self.is_sp2_carbon(atom_idx)
                for atom_idx in cycle
                if self.atomic_symbols[atom_idx] == "C")
            if (self.cycle_sizes[index] == 5 and
                    self.aromaticity[index] == 'non_aromatic' and
                    all_carbons_sp2):
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

        This method compiles the atomic symbol, atom number, number of
        connected atoms, symbols of connected atoms, atom numbers of connected
        atoms, distances to connected atoms, and cycle information into a
        structured dictionary.

        :return:
            The dictionary where each key is an atom number and each value is
            another dictionary of atom information.
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
        Analyzes the molecular structure information to assign atom types to
        each atom in the molecule.

        This method traverses through the atom information dictionary created
        from the molecular structure data and applies a series of rules to
        determine the appropriate atom type for each atom. The rules consider
        factors such as the atom's chemical environment, its connectivity to
        other atoms, and whether it is part of a cyclic structure.

        :return:
            A dictionary where each key corresponds to an atom identifier
            (e.g., "C1" for the first carbon atom), and each value is another
            dictionary containing the 'opls' and 'gaff' force field identifiers
            for the atom.
        """

        self.bad_hydrogen = False

        self.atom_types_dict = {}

        for atom_number, info in self.atom_info_dict.items():

            if info['AtomicSymbol'] == 'C':

                # If this carbon was previously assigned, skip the rest of the
                # loop
                key = f"C{info['AtomNumber']}"
                if key in self.atom_types_dict:
                    continue

                # Chemical environment information
                connected_symbols = set(info['ConnectedAtoms'])
                connected_distances = info['ConnectedAtomsDistances']

                # Cyclic and pure aromatic
                if (info.get('CyclicStructure') == 'cycle' and
                        'pure_aromatic' in info.get('Aromaticity')):

                    if info['NumConnectedAtoms'] == 3:
                        # Check for identifying biphenyls
                        connected_carbons_in_diff_cycle_and_pure_aromatic = []

                        for connected_atom_number in info[
                                'ConnectedAtomsNumbers']:
                            connected_atom_info = self.atom_info_dict[
                                connected_atom_number]

                            if (connected_atom_info.get('AtomicSymbol') == 'C'
                                    and
                                    connected_atom_info.get('CycleNumber') and
                                    not set(
                                        connected_atom_info.get('CycleNumber'))
                                    and set(info.get('CycleNumber')) and
                                    'pure_aromatic'
                                    in connected_atom_info.get('Aromaticity')):

                                connected_carbons_in_diff_cycle_and_pure_aromatic.append(
                                    connected_atom_number)

                        # If the list is not empty, set the
                        # connected_carbon_atom to the first element. Else, set
                        # it to None.
                        connected_carbon_atom = (
                            connected_carbons_in_diff_cycle_and_pure_aromatic[0]
                            if connected_carbons_in_diff_cycle_and_pure_aromatic
                            else None)

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

                            # Store the atomtype in the dictionary for the
                            # biphenyl carbon
                            self.atom_types_dict[
                                f"C{connected_carbon_atom}"] = biphenyl_carbon
                        else:
                            carbon_type = {'opls': 'opls_145', 'gaff': 'ca'}

                    elif (info['NumConnectedAtoms'] == 4 and
                          connected_symbols == {'C', 'H'}):
                        carbon_type = {'opls': 'opls_135', 'gaff': 'c3'}

                    else:
                        carbon_type = {
                            'opls': f'opls_x{info["AtomNumber"]}',
                            'gaff': f'cx{info["AtomNumber"]}'
                        }

                # Default assignation for carbons in non-pure aromatic cyclic
                # structures
                elif (info.get('CyclicStructure') == 'cycle' and
                      'non_pure_aromatic' in info.get('Aromaticity')):

                    if 'O' in connected_symbols:
                        # Directly loop through the connected atom numbers
                        for connected_atom_number in info[
                                'ConnectedAtomsNumbers']:
                            if ((self.atom_info_dict[connected_atom_number]
                                 ['AtomicSymbol'] == 'O') and
                                (self.atom_info_dict[connected_atom_number]
                                 ['CyclicStructure'] == 'none')):
                                if self.atom_info_dict[connected_atom_number][
                                        'NumConnectedAtoms'] == 1:
                                    # Carbonyl Carbon
                                    carbon_type = {
                                        'opls': 'opls_235',
                                        'gaff': 'c'
                                    }
                                    # Exit the loop once the carbonyl carbon is
                                    # identified
                                    break
                            else:
                                carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}

                    elif ('C' in connected_symbols and
                          info['NumConnectedAtoms'] == 3 and
                          any(size in info['CycleSize'] for size in {4, 5, 6})):

                        # Default assignment
                        carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}

                        # TODO: This part of the code will be deprecated in the
                        # future. The method check_alternating_atoms() will be
                        # used instead.

                        # Immediately assign C atom types based on distances
                        connected_atoms_numbers = info['ConnectedAtomsNumbers']

                        # Current atom's cycle number (will be 'none' if not in
                        # a cycle)
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

                                # Check if neighboring atom is not in a cycle
                                # or belongs to another cycle ID
                                if (neighboring_cycle_number == 'none' or
                                        neighboring_cycle_number
                                        != current_cycle_number):
                                    continue
                                else:
                                    if distance <= 1.4:
                                        # Double bonded
                                        conn_carbon_type = {
                                            'opls': 'opls_XXX',
                                            'gaff': 'cd'
                                        }
                                    else:
                                        # Single bonded
                                        conn_carbon_type = {
                                            'opls': 'opls_508',
                                            'gaff': 'cc'
                                        }

                                    # Store the atomtype in the dictionary
                                    atom_num = neighboring_atom_number
                                    self.atom_types_dict[
                                        f"C{atom_num}"] = conn_carbon_type

                                    # If there is an H connected to this
                                    # carbon, assign it now
                                    for connected_atom_number in self.atom_info_dict[
                                            neighboring_atom_number][
                                                'ConnectedAtomsNumbers']:
                                        connected_atom_info = self.atom_info_dict[
                                            connected_atom_number]
                                        if (connected_atom_info['AtomicSymbol']
                                                == 'H' and connected_atom_info[
                                                    'NumConnectedAtoms'] == 1):
                                            ewd_atoms = [
                                                'N', 'Br', 'Cl', 'I', 'F', 'S',
                                                'O'
                                            ]
                                            ewd_count = sum([
                                                1
                                                for num in self.atom_info_dict[
                                                    neighboring_atom_number]
                                                ['ConnectedAtomsNumbers']
                                                if self.atom_info_dict[num]
                                                ['AtomicSymbol'] in ewd_atoms
                                            ])
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
                                                # Default Aliphatic sp2
                                                # Hydrogen for c2
                                                hydrogen_type = {
                                                    'opls': 'opls_146',
                                                    'gaff': 'ha'
                                                }

                                            # Store it in the dictionary
                                            atom_num = connected_atom_info[
                                                'AtomNumber']
                                            self.atom_types_dict[
                                                f'H{atom_num}'] = hydrogen_type

                                    # Next iteration
                                    continue
                    else:
                        carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}

                # Aliphatic cycles
                elif (info.get('CyclicStructure') == 'cycle' and
                      'non_aromatic' in info.get('Aromaticity')):

                    if (info['NumConnectedAtoms'] == 4 and
                            3 in info['CycleSize']):
                        carbon_type = {'opls': 'opls_CX', 'gaff': 'cx'}
                    elif (info['NumConnectedAtoms'] == 4 and
                          4 in info['CycleSize']):
                        carbon_type = {'opls': 'opls_CY', 'gaff': 'cy'}
                    elif (info['NumConnectedAtoms'] == 4 and
                          5 in info['CycleSize']):
                        carbon_type = {'opls': 'opls_c5', 'gaff': 'c5'}
                    elif (info['NumConnectedAtoms'] == 4 and
                          6 in info['CycleSize']):
                        carbon_type = {'opls': 'opls_c6', 'gaff': 'c6'}
                    elif (info['NumConnectedAtoms'] == 3 and
                          3 in info['CycleSize']):
                        carbon_type = {'opls': 'opls_CU', 'gaff': 'cu'}
                    elif (info['NumConnectedAtoms'] == 3 and
                          4 in info['CycleSize']):
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
                                        # Carbonyl Carbon
                                        carbon_type = {
                                            'opls': 'opls_235',
                                            'gaff': 'c'
                                        }
                                        # Exit the loop once the carbonyl
                                        # carbon is identified
                                        break
                                    else:
                                        carbon_type = {
                                            'opls': 'opls_141',
                                            'gaff': 'c2'
                                        }  # Alcohol

                        elif 'C' in connected_symbols:
                            # Count the number of sp2-hybridized carbons
                            # connected to the current carbon
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
                            # If the current carbon is connected to 2 sp2 carbons
                            if (sp2_carbon_count + sp1_carbon_count == 2 or
                                    sp2_carbon_count + sp1_carbon_count == 3):
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
                            }  # Generic sp2 C

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
                                        # Carbonyl Carbon
                                        carbon_type = {
                                            'opls': 'opls_235',
                                            'gaff': 'c'
                                        }
                                        # Exit the loop once the carbonyl
                                        # carbon is identified
                                        break
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
                                        # Carbon double bonded to Sulfur
                                        carbon_type = {
                                            'opls': 'opls_cs',
                                            'gaff': 'cs'
                                        }
                                        # Exit the loop as soon as we find a
                                        # sulfur atom with only one connected
                                        # atom
                                        break
                                    else:
                                        carbon_type = {
                                            'opls': 'opls_141',
                                            'gaff': 'c2'
                                        }

                        elif 'C' in connected_symbols:
                            # Count the number of sp2-hybridized carbons
                            # connected to the current carbon
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
                            if (sp2_carbon_count + sp1_carbon_count == 2 or
                                    sp2_carbon_count + sp1_carbon_count == 3):
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
                            }  # Generic sp2 C

                    elif info['NumConnectedAtoms'] == 2:
                        if 'O' in connected_symbols:
                            # Carbon in carbonyl group or acid anhydride
                            carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}
                        else:
                            # Count the number of sp2-hybridized carbons
                            # connected to the current carbon
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

                # TODO: Simplify the rules that can be in the same types
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
                            # Hydrogens connected to C in heterocycle with N as
                            # in pyridine C-N-C
                            hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}

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
                        }:
                            # Hydrogens connected to a non-pure aromatic cycle
                            # or non_aromatic cycle
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

                        # Hydrogens connected to a carbon double bonded to a
                        # sulfur
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
                            # This checks if the carbon connected to the oxygen
                            # is connected to another oxygen.
                            # It is useful to identify carboxylic acids and
                            # esters.
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
                            # This checks if the carbon connected to the oxygen
                            # is connected to another oxygen.
                            # It is useful to identify carboxylic acids and
                            # esters.
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

                        # This case is highly dependent on the environment of
                        # the nitrogen
                        # Create flags to check for specific cases

                        # List of flags
                        found_nitro = False
                        found_aromatic = False
                        found_amide = False
                        found_idine = False
                        found_sp2_carbon = False
                        found_sp1_carbon = False
                        found_sp3_carbon = False

                        if (sorted(info['ConnectedAtoms']) == sorted(
                            ['C', 'O', 'O']) or sorted(info['ConnectedAtoms'])
                                == sorted(['N', 'O', 'O'])):
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
                                        'NumConnectedAtoms'] == 4:
                                    # Connected to an sp3 carbon
                                    found_sp3_carbon = True

                                elif self.atom_info_dict[atom][
                                        'NumConnectedAtoms'] == 3:
                                    # Connected to an sp2 carbon
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
                                        elif (atom_symbol == 'N' and
                                              atom_connectivity == 2):
                                            found_idine = True

                                        # Connected to C sp1
                                        elif (atom_symbol == 'C' and
                                              atom_connectivity == 2):
                                            found_sp1_carbon = True

                                        # Connected to C sp2
                                        elif (atom_symbol == 'C' and
                                              atom_connectivity == 3):
                                            found_sp2_carbon = True
                                else:
                                    # n3 in GAFF but special cases n7 and n8
                                    # added in GAFF2
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

                        # Now assign nitrogen types based on the flags using
                        # the following hierarchy:
                        # 1. Nitro
                        # 2. Amide
                        # 3. Idine
                        # 4. Aromatic
                        # 5. Sp2 carbon
                        # 6. Sp1 carbon
                        # 7. Sp3 carbon

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
                        # Check if the Nitrogen is connected to another
                        # Nitrogen with sp1 hybridization
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

                    if (info['NumConnectedAtoms'] == 2 and
                            'pure_aromatic' in info.get('Aromaticity')):
                        nitrogen_type = {
                            'opls': 'opls_520',
                            'gaff': 'nb'
                        }  # Sp2 N in pure aromatic systems

                    elif (info['NumConnectedAtoms'] == 2 and
                          'non_pure_aromatic' in info.get('Aromaticity')):
                        nitrogen_type = {
                            'opls': 'opls_520',
                            'gaff': 'nc'
                        }  # Sp2 N in non-pure aromatic systems

                    elif (info['NumConnectedAtoms'] == 3 and
                          'pure_aromatic' in info.get('Aromaticity')):
                        nitrogen_type = {
                            'opls': 'opls_183',
                            'gaff': 'nb'
                        }  # Pyridine as a ligand in an organometallic complex

                    elif (info['NumConnectedAtoms'] == 3 and
                          'non_pure_aromatic' in info.get('Aromaticity')):
                        if 'C' in connected_symbols:
                            # Check for amides and sulfamides by checking the
                            # connected atoms to the carbon
                            found_CO = False

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
                                        found_CO = True
                                        break
                                if found_CO:
                                    break
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
                    elif (info['NumConnectedAtoms'] == 3 and
                          'non_aromatic' in info.get('Aromaticity')):
                        if 'C' in connected_symbols:
                            # Check for amides and sulfamides by checking the
                            # connected atoms to the carbon
                            found_CO = False

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
                                        found_CO = True
                                        break
                                if found_CO:
                                    break

                            if not found_CO:
                                # Check if the nitrogen is connected to a carbon sp2
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
                                    # n3 in GAFF but special cases n7 and n8
                                    # added in GAFF2
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

                        # Check for special RG3 and RG4 cases
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
                        # Add other conditions for cyclic nitrogen atoms if
                        # needed or add a default case
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
                            # Again simplified, could distinguish between
                            # phosphonates and phosphites
                            phosphorus_type = {
                                'opls': 'opls_903P',
                                'gaff': 'p3'
                            }
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

                    # Note: Hydrogens in phosphine groups are less commonly
                    # parameterized in force fields and hence not added here.

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
                    elif all([
                            self.atom_info_dict[num]['AtomicSymbol']
                            in ['C', 'N', 'S', 'O']
                            for num in info['ConnectedAtomsNumbers']
                    ]):
                        # Both connected atoms are carbons or one carbon and
                        # one nitrogen
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
                    atom_num = info['ConnectedAtomsNumbers'][
                        info['ConnectedAtoms'].index('H')]
                    self.atom_types_dict[f'H{atom_num}'] = hydrogen_type
                    self.bad_hydrogen = True

            elif info['AtomicSymbol'] == 'Cl':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'cl'}
                self.atom_types_dict[
                    f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type
                if 'H' in info['ConnectedAtoms']:
                    hydrogen_type = {'opls': 'opls_h_x', 'gaff': 'h_x'}
                    atom_num = info['ConnectedAtomsNumbers'][
                        info['ConnectedAtoms'].index('H')]
                    self.atom_types_dict[f'H{atom_num}'] = hydrogen_type
                    self.bad_hydrogen = True

            elif info['AtomicSymbol'] == 'F':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'f'}
                self.atom_types_dict[
                    f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type
                if 'H' in info['ConnectedAtoms']:
                    hydrogen_type = {'opls': 'opls_h_x', 'gaff': 'h_x'}
                    atom_num = info['ConnectedAtomsNumbers'][
                        info['ConnectedAtoms'].index('H')]
                    self.atom_types_dict[f'H{atom_num}'] = hydrogen_type
                    self.bad_hydrogen = True

            elif info['AtomicSymbol'] == 'I':
                halogen_type = {'opls': 'opls_XXX', 'gaff': 'i'}
                self.atom_types_dict[
                    f"{info['AtomicSymbol']}{info['AtomNumber']}"] = halogen_type
                if 'H' in info['ConnectedAtoms']:
                    hydrogen_type = {'opls': 'opls_h_x', 'gaff': 'h_x'}
                    atom_num = info['ConnectedAtomsNumbers'][
                        info['ConnectedAtoms'].index('H')]
                    self.atom_types_dict[f'H{atom_num}'] = hydrogen_type
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
                    atom_num = info['ConnectedAtomsNumbers'][
                        info['ConnectedAtoms'].index('H')]
                    self.atom_types_dict[f'H{atom_num}'] = hydrogen_type
                    self.bad_hydrogen = True

            else:
                # Else for the cases falling off the decision tree
                # The Hydrogen are assigned outside the main branches of the
                # decision tree
                # Therefore, they need to be out of the else case.

                if info['AtomicSymbol'] != 'H':
                    self.ostream.print_warning(
                        f"{info['AtomicSymbol']}{info['AtomNumber']}" +
                        'Has not been found in the decision tree, check it carefully'
                    )

        return self.atom_types_dict

    def extract_gaff_atom_types(self, atom_type):
        """
        Extracts GAFF atom types from the atom types dictionary.

        """

        # Initialize the list of gaff atom types
        self.gaff_atom_types = []

        # Sort atom types based on the number after the atomic symbol
        sorted_atom_types = sorted(self.atom_types_dict.keys(),
                                   key=self.get_atom_number)

        # Append the gaff atom types to the list
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

        # Check distances for alternation
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

        # Return an array with the new atom types
        self.gaff_atom_types = atom_types

        return self.gaff_atom_types

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
        self.covalent_radii = (molecule.covalent_radii_to_numpy() *
                               bohr_in_angstrom())
        self.create_connectivity_matrix()
        self.detect_closed_cyclic_structures()
        self.create_atom_info_dict()
        self.atom_types_dict = self.decide_atom_type()
        self.extract_gaff_atom_types(self.atom_types_dict)
        self.gaff_atom_types = self.check_alternating_atom_types()

        # Printing output
        self.ostream.print_info("VeloxChem Atom Type Identification")
        self.ostream.print_info("-" * 40)  # Dashed line

        # Detected number of atoms
        num_atoms = len(self.atomic_symbols)
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

        if self.cycle_ids:
            cycle_suffix = 's' if len(self.cycle_ids) > 1 else ''
            self.ostream.print_blank()
            self.ostream.print_info(
                f"Detected {len(self.cycle_ids)} cycle{cycle_suffix}:")

        # Print cycle information (aromaticity) for each cycle
        for size, aromaticity in zip(self.cycle_sizes, self.aromaticity):
            if aromaticity == "non_aromatic":
                self.ostream.print_info(
                    f"Cycle size {size}: Non-aromatic Cycle")
            elif aromaticity == "non_pure_aromatic":
                self.ostream.print_info(
                    f"Cycle size {size}: Non-pure Aromatic Cycle")
            elif aromaticity == "pure_aromatic":
                self.ostream.print_info(
                    f"Cycle size {size}: Pure Aromatic Cycle")

        if self.bad_hydrogen:
            self.ostream.print_warning('Hydrogen type not defined in GAFF')

        return self.gaff_atom_types

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

    @staticmethod
    def measure_length(v1, v2):
        """
        Calculates the distance between v1 and v2
        """

        return np.linalg.norm(np.array(v1) - np.array(v2))
