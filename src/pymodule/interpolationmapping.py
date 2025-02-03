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
import itertools
from pathlib import Path
from sys import stdout
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom
import random
import re
from contextlib import redirect_stderr
from io import StringIO
from scipy.optimize import linear_sum_assignment
from .errorhandler import assert_msg_critical

from .molecularbasis import MolecularBasis
# from .imdatabasedriver import IMDatabaseDriver
#from .impesforcefieldgenerator_parallel import ImpesForceFieldGeneratorParallel
from .forcefieldgenerator import ForceFieldGenerator


# Drivers
from .scfrestdriver import ScfRestrictedDriver
from .molecularbasis import MolecularBasis
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .externalqmdriver import ExternalQMDriver
from .externalqmgradientdriver import ExternalQMGradientDriver
from .externalqmhessiandriver import ExternalQMHessianDriver
from .impesdriver import ImpesDriver
from .impescoordinates import ImpesCoordinates

from .atomtypeidentifier import AtomTypeIdentifier
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from .molecule import Molecule
from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom

with redirect_stderr(StringIO()) as fg_err:
    import geometric

class InterpolationMapping:


    def __init__(self, molecule, non_core_symmetry_groups=None):

        self.init_variable = None
        self.ff_data = None
        self.molecule = molecule
        
        if non_core_symmetry_groups is None:
            non_core_symmetry_groups = self.symmetric_groups(molecule)
        self.non_core_symmetry_groups = non_core_symmetry_groups
    
    def calculate_translation_coordinates_analysis(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates

    def cost_matrix_analysis(self, current_coordinates, datapoint_coordinate):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                ImpesCoordinates object
        """
        
        non_core_symmetry_groups = self.non_core_symmetry_groups
        if non_core_symmetry_groups is None:
            non_core_symmetry_groups = self.symmetric_groups(self.molecule)

        # First, translate the cartesian coordinates to zero
        target_coordinates = self.calculate_translation_coordinates_analysis(datapoint_coordinate)
        reference_coordinates = self.calculate_translation_coordinates_analysis(current_coordinates)
        # remove the non_core atoms (currently just H)
        target_coordinates_core = np.delete(target_coordinates, [], axis=0)
        reference_coordinates_core = np.delete(reference_coordinates, [], axis=0)
        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)     
        rotation_matrix_core = geometric.rotate.get_rot(target_coordinates_core,
                                                reference_coordinates_core)
        

        # Rotate the data point
        rotated_coordinates_core = np.dot(rotation_matrix_core, target_coordinates.T).T
        
        if len(non_core_symmetry_groups) == 0:

            # Calculate the Cartesian distance
            ref_structure_check = reference_coordinates.copy()
            distance_core = (np.linalg.norm(rotated_coordinates_core - ref_structure_check))
        
        else:

        
            reference_group_idx = [idx for idx, _ in enumerate(self.molecule.get_labels()) if idx not in non_core_symmetry_groups]
            reference_group = np.delete(np.array(reference_coordinates), reference_group_idx, axis=0)
            current_group = np.delete(np.array(rotated_coordinates_core), reference_group_idx, axis=0)
        
            
            row_ind, col_ind = self.assign_atoms(current_group, reference_group)
            #define the overall cost of the matrix

            swapped_indices = {}
            swapped_elements_map = {}
            swapped_indices_list = []
            reference_indices_list = []
            indices_that_need_to_be_swapped = []
            element = 0
            current_element_list = []
            for i, j in zip(row_ind, col_ind):
                swapped_indices[non_core_symmetry_groups[i]] = non_core_symmetry_groups[j]
                swapped_indices_list.append(non_core_symmetry_groups[i])
                current_element_list.append(non_core_symmetry_groups[j])
                reference_indices_list.append(non_core_symmetry_groups[j])
                indices_that_need_to_be_swapped.append(non_core_symmetry_groups[j])
                if len(current_element_list) == 3:
                    swapped_elements_map[element] = tuple(current_element_list)
                    current_element_list = []
                    element += 1

            y_assigned = reference_coordinates[reference_indices_list]
            ref_structure_check = reference_coordinates.copy()
            ref_structure_check[swapped_indices_list] = y_assigned

            # Calculate the Cartesian distance
            distance_core = (np.linalg.norm(rotated_coordinates_core - ref_structure_check))
        
        #print('Distance norm', distance_core)
        # Calculate the gradient of the interpolation weights
        # (required for energy gradient interpolation)
        distance_vector_core = (ref_structure_check - rotated_coordinates_core)
        distance_vector_core_norm = np.zeros(reference_coordinates_core.shape[0])

        for i in range(len(distance_vector_core_norm)):
            distance_vector_core_norm[i] += np.linalg.norm(distance_vector_core[i])

        return ref_structure_check, distance_core, distance_vector_core
    
    def assign_atoms(self, x, y):
    
        cost_matrix = np.linalg.norm(x[:, np.newaxis, :] - y[np.newaxis, :, :], axis=2)
        
        # Apply the Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        
        return row_ind, col_ind
    
    
    def symmetric_groups(self, molecule=None, ff_gen=None):
        

        if self.ff_data is None and ff_gen is None:
            assert_msg_critical(self.ff_data is None, 'FF_data not defined: Can not determine symmetric groups.')
        
        if molecule is not None:
            self.molecule = molecule
        

        if ff_gen is None:
            ff_gen = ForceFieldGenerator()
            ff_gen.force_field_data = self.ff_data
            ff_gen.partial_charges = self.molecule.get_partial_charges(self.molecule.get_charge())
            ff_gen.create_topology(self.molecule)
        
        connectivity_matrix = self.molecule.get_connectivity_matrix()
        rotatable_bonds = ff_gen.rotatable_bonds
        dihedral_side_atom_types = {}            
        atom_types = ff_gen.atom_types
        
        

        for rot_bond in rotatable_bonds:
        
            atom_1_connections = [(index, atom_types[index]) for index, connected in enumerate(connectivity_matrix[rot_bond[0] - 1]) if connected == 1 and index != rot_bond[1] - 1]
            atom_2_connections = [(index, atom_types[index])  for index, connected in enumerate(connectivity_matrix[rot_bond[1] - 1]) if connected == 1 and index != rot_bond[0] - 1]

            if rot_bond[0] not in dihedral_side_atom_types:

                dihedral_side_atom_types[rot_bond[0]-1] = atom_1_connections
                    
            if rot_bond[1] not in dihedral_side_atom_types:

                dihedral_side_atom_types[rot_bond[1]-1] = atom_2_connections
            
        index_dict = {}

        for keys, diehdral_group in dihedral_side_atom_types.items():

            index_dict[keys] = {}

            for tuple in diehdral_group:

                if tuple[1] not in index_dict[keys]:
                    index_dict[keys][tuple[1]] = []
                
                index_dict[keys][tuple[1]].append(tuple[0])
        
        non_core_symmetric_groups = {}
        non_core_symmetric_groups_combinded = []
        for key, groups in index_dict.items():

            for _, atom_group in groups.items():
                
                if len(atom_group) <= 2:
                    continue
                
                if key not in non_core_symmetric_groups:
                     non_core_symmetric_groups[key] = atom_group
                                    
                non_core_symmetric_groups_combinded += atom_group
        
        return non_core_symmetric_groups_combinded