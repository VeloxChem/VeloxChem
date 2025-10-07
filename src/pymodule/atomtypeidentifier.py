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
from .outputstream import OutputStream
from .errorhandler import safe_arccos


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

        return (self.atomic_symbols[atom_idx] == 'N'
                and len(list(self.graph.neighbors(atom_idx))) == 2)

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
                    connected_atom_info, hydrogen_type = self.assign_hydrogen_carbon_type(
                        atom_info_dict, atom_type, connected_atom_number, info)
                    if hydrogen_type is not None:
                        H_name = f"H{connected_atom_info['AtomNumber']}"
                        atom_types_dict.setdefault(H_name,
                                                   {}).update(hydrogen_type)

            elif info['AtomicSymbol'] == 'O':
                atom_type = self.assign_oxygen_type(atom_info_dict, info,
                                                    connected_symbols)

                # Hydrogen type assignment based on oxygen type
                for connected_atom_number in info['ConnectedAtomsNumbers']:
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
                    # If px then hx, else hp
                    for connected_atom_number in info['ConnectedAtomsNumbers']:
                        hydrogen_type = None
                        connected_atom_info, hydrogen_type = self.assign_hydrogen_phosphorus_type(
                            atom_info_dict, connected_atom_number)

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

            elif info['AtomicSymbol'] in ['F', 'Cl', 'Br', 'I']:
                atom_type = {
                    'opls': 'opls_XXX',
                    'gaff': info['AtomicSymbol'].lower()
                }

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
            if sulfur_type == {'opls': 'opls_924S', 'gaff': 'sh'}:
                hydrogen_type = {'opls': 'opls_926H', 'gaff': 'hs'}
        return connected_atom_info, hydrogen_type

    def assign_sulfur_type(self, atom_info_dict, info, connected_symbols):
        sulfur_type = None
        # S with one connected atom
        if info['NumConnectedAtoms'] == 1:

            sulfur_type = {'opls': 'opls_920S', 'gaff': 's'}

        # S with two connected atom, involved at least one double bond
        elif info['NumConnectedAtoms'] == 2:

            if 'H' in connected_symbols:
                sulfur_type = {'opls': 'opls_924S', 'gaff': 'sh'}

            elif all([
                    atom_info_dict[num]['AtomicSymbol'] in ['C', 'N', 'S', 'O']
                    for num in info['ConnectedAtomsNumbers']
            ]):
                # Both connected atoms are carbons or one carbon and
                # one nitrogen
                # Thio-ether or Thio-ester
                sulfur_type = {'opls': 'opls_SS', 'gaff': 'ss'}

            else:
                sulfur_type = {'opls': 'opls_921S', 'gaff': 's2'}

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
                sulfur_type = {'opls': 'opls_922X', 'gaff': 'sx'}

            else:
                sulfur_type = {'opls': 'opls_922S', 'gaff': 's4'}

        # S with four connected atoms
        elif info['NumConnectedAtoms'] == 4:

            if any(atom_info_dict[num]['AtomicSymbol'] == 'C'
                   and atom_info_dict[num]['NumConnectedAtoms'] == 3
                   for num in info['ConnectedAtomsNumbers']):
                sulfur_type = {'opls': 'opls_922X', 'gaff': 'sy'}

            else:
                sulfur_type = {'opls': 'opls_923S', 'gaff': 's6'}

            # TODO: Sp3 S connected with hydrogen

        return sulfur_type

    def assign_hydrogen_phosphorus_type(self, atom_info_dict,
                                        connected_atom_number):
        hydrogen_type = None
        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):
            hydrogen_type = {'opls': 'opls_XXX', 'gaff': 'hp'}
        return connected_atom_info, hydrogen_type

    def assign_phosphorus_type(self, info, connected_symbols):
        phosphorus_type = None
        # hypervalent phosphorus, 4 subst.
        if (info['NumConnectedAtoms'] == 4 and 'O' in connected_symbols):
            phosphorus_type = {'opls': 'opls_900P', 'gaff': 'p5'}

            # sp3 phosphorus, 3 subst.
        elif info['NumConnectedAtoms'] == 3:
            # Oxygen determines if the phosphorus is hypervalent or not
            if 'O' in connected_symbols:
                oxygen_count = info['ConnectedAtoms'].count('O')
            else:
                oxygen_count = 0

                # Regular sp3 P with three connected atoms, such as PH3
            if oxygen_count == 0:
                phosphorus_type = {'opls': 'opls_901P', 'gaff': 'p3'}

                #  hypervalent phosphorus, 3 subst.
            else:
                phosphorus_type = {'opls': 'opls_900P', 'gaff': 'p4'}

                # sp2 phosphorus (C=P, etc.)
        elif info['NumConnectedAtoms'] == 2:
            phosphorus_type = {'opls': 'opls_900P', 'gaff': 'p2'}
        return phosphorus_type

    def assign_hydrogen_nitrogen_type(self, atom_info_dict,
                                      connected_atom_number):
        hydrogen_type = None
        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):
            hydrogen_type = {'opls': 'opls_240', 'gaff': 'hn'}
        return connected_atom_info, hydrogen_type

    def assign_hydrogen_oxygen_type(self, atom_info_dict, connected_atom_number,
                                    oxygen_type):
        hydrogen_type = None
        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):
            if oxygen_type == {'opls': 'opls_111', 'gaff': 'ow'}:
                hydrogen_type = {'opls': 'opls_112', 'gaff': 'hw'}

            elif oxygen_type == {'opls': 'opls_154', 'gaff': 'oh'}:
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
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n+'}

                elif num_hydrogens == 3:
                    # Sp3 N with three hydrogen atoms
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'nz'}

                elif num_hydrogens == 2:
                    # Sp3 N with two hydrogen atoms
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'ny'}

                elif num_hydrogens == 1:
                    # Sp3 N with one hydrogen atom
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'nx'}

                else:
                    # Sp3 N with four connected atoms, but no hydrogens
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n4'}

            elif info['NumConnectedAtoms'] == 3:
                # This case is highly dependent on the environment of
                # the nitrogen
                # Create flags to check for specific cases
                # List of flags
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

                        # Now assign nitrogen types based on the flags using
                        # the following hierarchy:
                        # 1. Nitro
                        # 2. Amide
                        # 3. Aromatic / sp2 carbon

                if found_nitro:
                    # Nitro N
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'no'}

                elif found_amide:
                    if num_hydrogens == 1:
                        nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'ns'}

                    elif num_hydrogens == 2:
                        nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'nt'}

                    else:
                        nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n'}

                elif found_aromatic or found_sp2_carbon:
                    if num_hydrogens == 1:
                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nu'}

                    elif num_hydrogens == 2:
                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nv'}

                    else:
                        nitrogen_type = {'opls': 'opls_901', 'gaff': 'nh'}

                else:
                    if num_hydrogens == 1:
                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n7'}

                    elif num_hydrogens == 2:
                        nitrogen_type = {'opls': 'opls_300', 'gaff': 'n8'}

                    else:
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
                    nitrogen_type = {'opls': 'opls_753', 'gaff': 'n1'}

                elif sp2_carbon_count + sp1_carbon_count + n2_count == 2:
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'ne'}

                else:
                    # Check bond angle
                    vec_ji = (self.coordinates[connected_atoms_numbers[0] - 1] -
                              self.coordinates[info['AtomNumber'] - 1])
                    vec_jk = (self.coordinates[connected_atoms_numbers[1] - 1] -
                              self.coordinates[info['AtomNumber'] - 1])
                    theta_ijk = safe_arccos(
                        np.dot(vec_ji, vec_jk) /
                        (np.linalg.norm(vec_ji) * np.linalg.norm(vec_jk)))

                    if theta_ijk * 180.0 / np.pi > 170.0:
                        nitrogen_type = {'opls': 'opls_753', 'gaff': 'n1'}
                    else:
                        nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n2'}

            elif info['NumConnectedAtoms'] == 1:
                nitrogen_type = {'opls': 'opls_753', 'gaff': 'n1'}

                # Cyclic

        elif info.get('CyclicStructure') == 'cycle':
            if (info['NumConnectedAtoms'] == 2
                    and 'pure_aromatic' in info.get('Aromaticity')):
                # Sp2 N in pure aromatic systems
                nitrogen_type = {'opls': 'opls_520', 'gaff': 'nb'}

            elif (info['NumConnectedAtoms'] == 2
                  and 'non_pure_aromatic' in info.get('Aromaticity')):
                # Sp2 N in non-pure aromatic systems
                nitrogen_type = {'opls': 'opls_520', 'gaff': 'nc'}

            elif (info['NumConnectedAtoms'] == 3
                  and 'pure_aromatic' in info.get('Aromaticity')):
                # Pyridine as a ligand in an organometallic complex
                nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'nb'}

            elif (info['NumConnectedAtoms'] == 3
                  and 'non_pure_aromatic' in info.get('Aromaticity')):
                # Default assignment
                # General n3 case
                nitrogen_type = {'opls': 'opls_na', 'gaff': 'na'}

                if 'C' in connected_symbols:
                    # Check for amides and sulfamides by checking the
                    # connected atoms to the carbon
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

                    # Nitrogens in Non-aromatic cycles
            elif (info['NumConnectedAtoms'] == 3
                  and 'non_aromatic' in info.get('Aromaticity')):
                if 'C' in connected_symbols:
                    # Check for amides and sulfamides by checking the
                    # connected atoms to the carbon
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
                        # Check if the nitrogen is connected to a sp2 C or N
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
                            # n3 and special cases
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

                        # Check for special RG3 and RG4 cases
                else:
                    if num_hydrogens == 1:
                        if 3 in info.get('CycleSize'):
                            nitrogen_type = {'opls': 'opls_n5', 'gaff': 'n5'}

                        elif 4 in info.get('CycleSize'):
                            nitrogen_type = {'opls': 'opls_n6', 'gaff': 'n6'}

                        else:
                            nitrogen_type = {'opls': 'opls_300', 'gaff': 'n7'}

                    elif num_hydrogens == 0:
                        if 3 in info.get('CycleSize'):
                            nitrogen_type = {'opls': 'opls_np', 'gaff': 'np'}

                        elif 4 in info.get('CycleSize'):
                            nitrogen_type = {'opls': 'opls_nq', 'gaff': 'nq'}

                        else:
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
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'ne'}

                else:
                    nitrogen_type = {'opls': 'opls_XXX', 'gaff': 'n2'}
        return nitrogen_type

    def assign_carbon_type(self, atom_info_dict, using_gaff_220, info,
                           connected_symbols):
        carbon_type = None
        # Pure aromatic cycles

        if (info.get('CyclicStructure') == 'cycle'
                and 'pure_aromatic' in info.get('Aromaticity')):
            if info['NumConnectedAtoms'] == 3:
                # Check for identifying biphenyls
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
                    carbon_type = {'opls': 'opls_145', 'gaff': 'ca'}

                elif connected_symbols == {'C', 'N', 'H'}:
                    carbon_type = {'opls': 'opls_521', 'gaff': 'ca'}

                elif connected_carbon_atom is not None:
                    carbon_type = {'opls': 'opls_521', 'gaff': 'cp'}

                else:
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
                # Carbonyl Carbon
                carbon_type = {'opls': 'opls_235', 'gaff': 'c'}

            elif has_terminal_sulfur:
                # Carbon double bonded to Sulfur
                carbon_type = {'opls': 'opls_cs', 'gaff': 'cs'}

            else:
                carbon_type = {'opls': 'opls_508', 'gaff': 'cc'}

                # Non-aromatic cycles

        elif (info.get('CyclicStructure') == 'cycle'
              and 'non_aromatic' in info.get('Aromaticity')):
            if info['NumConnectedAtoms'] == 4:
                if 3 in info['CycleSize']:
                    carbon_type = {'opls': 'opls_CX', 'gaff': 'cx'}

                elif 4 in info['CycleSize']:
                    carbon_type = {'opls': 'opls_CY', 'gaff': 'cy'}

                elif 5 in info['CycleSize'] and using_gaff_220:
                    carbon_type = {'opls': 'opls_c5', 'gaff': 'c5'}

                elif 6 in info['CycleSize'] and using_gaff_220:
                    carbon_type = {'opls': 'opls_c6', 'gaff': 'c6'}

                else:
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
                    # Carbonyl Carbon
                    carbon_type = {'opls': 'opls_235', 'gaff': 'c'}

                elif has_terminal_sulfur:
                    # Carbon double bonded to Sulfur
                    carbon_type = {'opls': 'opls_cs', 'gaff': 'cs'}

                elif 3 in info['CycleSize']:
                    carbon_type = {'opls': 'opls_CU', 'gaff': 'cu'}

                elif 4 in info['CycleSize']:
                    carbon_type = {'opls': 'opls_CV', 'gaff': 'cv'}

                elif 'C' in connected_symbols or 'N' in connected_symbols:
                    # Count the number of sp2/sp1 hybridized C and N
                    # connected to the current carbon
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
                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'ce'}

                    else:
                        carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}

                else:
                    carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}

            elif info['NumConnectedAtoms'] == 2:
                carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}

                # Non-cyclic

        elif info.get('CyclicStructure') == 'none':
            if info['NumConnectedAtoms'] == 4:
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
                    # Carbonyl Carbon
                    carbon_type = {'opls': 'opls_235', 'gaff': 'c'}

                elif has_terminal_sulfur:
                    # Carbon double bonded to Sulfur
                    carbon_type = {'opls': 'opls_cs', 'gaff': 'cs'}

                elif 'C' in connected_symbols or 'N' in connected_symbols:
                    # Count the number of sp2/sp1 hybridized C and N
                    # connected to the current carbon
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
                        carbon_type = {'opls': 'opls_XXX', 'gaff': 'ce'}

                    else:
                        carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}

                else:
                    carbon_type = {'opls': 'opls_141', 'gaff': 'c2'}

            elif info['NumConnectedAtoms'] == 2:
                if 'O' in connected_symbols:
                    # Carbon in carbonyl group or acid anhydride
                    carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}

                else:
                    # Count the number of sp2/sp1 hybridized C and N
                    # connected to the current carbon
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
                        carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}

                    elif (sp2_carbon_count + sp1_carbon_count + n_count == 2 or
                          sp2_carbon_count + sp1_carbon_count + n2_count == 2):
                        # Check bond length
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
                            carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}
                        else:
                            carbon_type = {'opls': 'opls_XXX', 'gaff': 'cg'}

                    else:
                        carbon_type = {'opls': 'opls_235', 'gaff': 'c1'}

                # Assignment for the Hydrogens linked to the carbons
                # ewd = Electron withdrawing atoms

        return carbon_type

    def assign_hydrogen_carbon_type(self, atom_info_dict, carbon_type,
                                    connected_atom_number, info):
        hydrogen_type = None

        ewd_atoms = ['N', 'Br', 'Cl', 'I', 'F', 'S', 'O']

        ewd_count = sum(1 for num in info['ConnectedAtomsNumbers']
                        if atom_info_dict[num]['AtomicSymbol'] in ewd_atoms)
        connected_atom_info = atom_info_dict[connected_atom_number]

        if (connected_atom_info['AtomicSymbol'] == 'H'
                and connected_atom_info['NumConnectedAtoms'] == 1):
            # sp1 carbon
            if carbon_type == {'opls': 'opls_235', 'gaff': 'c1'}:
                hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}

                # sp2 carbon
            elif carbon_type == {
                    'opls': 'opls_CU',
                    'gaff': 'cu'
            } or carbon_type == {
                    'opls': 'opls_CV',
                    'gaff': 'cv'
            } or carbon_type == {
                    'opls': 'opls_141',
                    'gaff': 'c2'
            } or carbon_type == {
                    'opls': 'opls_145',
                    'gaff': 'ca'
            } or carbon_type == {
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
            } or carbon_type == {
                    'opls': 'opls_cs',
                    'gaff': 'cs'
            } or carbon_type == {
                    'opls': 'opls_235',
                    'gaff': 'c'
            }:
                if ewd_count == 1:
                    hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}

                elif ewd_count == 2:
                    hydrogen_type = {'opls': 'opls_xxx', 'gaff': 'h5'}

                else:
                    hydrogen_type = {'opls': 'opls_146', 'gaff': 'ha'}

                    # Hydrogens connected to C in heterocycle with N as in
                    # pyridine C-N-C
            elif carbon_type == {'opls': 'opls_521', 'gaff': 'ca'}:
                hydrogen_type = {'opls': 'opls_146', 'gaff': 'h4'}

                # sp3 carbon
            elif carbon_type == {
                    'opls': 'opls_135',
                    'gaff': 'c3'
            } or carbon_type == {
                    'opls': 'opls_c5',
                    'gaff': 'c5'
            } or carbon_type == {
                    'opls': 'opls_c6',
                    'gaff': 'c6'
            } or carbon_type == {
                    'opls': 'opls_CX',
                    'gaff': 'cx'
            } or carbon_type == {
                    'opls': 'opls_CY',
                    'gaff': 'cy'
            }:
                if ewd_count == 1:
                    hydrogen_type = {'opls': 'opls_140', 'gaff': 'h1'}

                elif ewd_count == 2:
                    hydrogen_type = {'opls': 'opls_xxx', 'gaff': 'h2'}

                elif ewd_count == 3:
                    hydrogen_type = {'opls': 'opls_xxx', 'gaff': 'h3'}

                else:
                    hydrogen_type = {'opls': 'opls_140', 'gaff': 'hc'}

        return connected_atom_info, hydrogen_type

    def assign_oxygen_type(self, atom_info_dict, info, connected_symbols):
        # Non-cyclic
        oxygen_type = None
        if info.get('CyclicStructure') == 'none':

            if (info['NumConnectedAtoms'] == 2 and connected_symbols == {'H'}):

                oxygen_type = {'opls': 'opls_111', 'gaff': 'ow'}

            elif (info['NumConnectedAtoms'] == 2 and 'H' in connected_symbols):

                oxygen_type = {'opls': 'opls_154', 'gaff': 'oh'}

            elif info['NumConnectedAtoms'] == 2:

                oxygen_type = {'opls': 'opls_XXX', 'gaff': 'os'}

            elif info['NumConnectedAtoms'] == 1:

                # If necessary we could check if the carbon connected
                # to the oxygen is connected to another oxygen. It is
                # useful to identify carboxylic acids and esters.

                oxygen_type = {'opls': 'opls_XXX', 'gaff': 'o'}

        # Cyclic

        elif info.get('CyclicStructure') == 'cycle':

            if (info['NumConnectedAtoms'] == 2 and connected_symbols == {'H'}):

                oxygen_type = {'opls': 'opls_154', 'gaff': 'oh'}

            elif info['NumConnectedAtoms'] == 2:

                oxygen_type = {'opls': 'opls_XXX', 'gaff': 'os'}

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

                oxygen_type = {'opls': 'opls_XXX', 'gaff': 'o'}

            # else:

            #     oxygen_type = {
            #         'opls': f'opls_x{info["AtomNumber"]}',
            #         'gaff': f'ox{info["AtomNumber"]}'
            #     }

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
        # Need to the central pair of ce-ce-ce-ce chain
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
                if (i_count['ee_and_ef'] == 2 and i_count['ee'] > 1 and
                        j_count['ee_and_ef'] == 2 and j_count['ee'] > 1):
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
            if atom_type in ['cd', 'cf', 'ch', 'nd', 'nf', 'cq']:
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

                    swapped_path_types = [
                        tuple(
                            conjugated_atomtype_mapping.get(at, at)
                            for at in path) for path in path_types
                    ]
                    swapped_path_types = tuple(sorted(swapped_path_types))

                    if swapped_path_types in atom_paths:
                        atom_paths[swapped_path_types].add(idx)
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
