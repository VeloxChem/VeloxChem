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
from math import ceil
import concurrent.futures
import pandas as pd

from .molecularbasis import MolecularBasis
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
from .interpolationmapping import InterpolationMapping

from .atomtypeidentifier import AtomTypeIdentifier
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from .molecule import Molecule
from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class FindBestCombination:

    def __init__(self, molecule, ff_datafile):

        self.init_variable = 1.0
        self.orginal_molecule = molecule
        self.ff_datafile = ff_datafile

        forcefield_generator = ForceFieldGenerator()
        forcefield_generator.force_field_data = ff_datafile
        forcefield_generator.partial_charges = molecule.get_partial_charges(molecule.get_charge())
        forcefield_generator.create_topology(molecule)
        
        self.ff_gen = forcefield_generator
        self.atom_types = forcefield_generator.atom_types
        print('define the forcefield datafile (gaff-2.11.dat) in self.ff_datafile !!!')

    # parse category is used to combine different atom types to form internal bond, angle, dihedral types
    
    def parse_category(self, category, current_atom_types):
        pattern = '|'.join(sorted(current_atom_types, key=len, reverse=True)) 
        elements = re.findall(pattern, category)
        return elements

    # get_internal_coordinate_atom_type_mapping is used to determine the internal bond, angle, dihedral types
    def get_internal_coordinate_atom_type_mapping(self, molecule, atom_types):

        z_matrix = molecule.get_z_matrix()

        grouping = {}
        for internal_coord in z_matrix:

            current_atom_types = [atom_types[index] for index in internal_coord]
            
            combined_atom_type = ''.join([atom_types[index] for index in internal_coord])
            elements = self.parse_category(combined_atom_type, current_atom_types)
            reordered_category = ''.join(sorted(elements))
            if reordered_category in grouping:
                grouping[reordered_category].append(internal_coord)
                    
            
            else:
                grouping[reordered_category] = [internal_coord]
            
        return grouping
    
    # define_internal_coordinates_types() is there to determines the distance of the same internal coordinate feature types (bond, angle, dihedral) formed from the same atom types of the molecule
    # for that the internal coordinate indices the corresponding atom types and the respected values determined by VeloxChem are extracted from the functions within the notebook and the forcefield generator Class previously explained!

    def define_internal_coordinates_types(self, reference_molecules, distorted_molecule_list, atom_types=None):

        atom_types = self.atom_types
        z_mat = self.orginal_molecule.get_z_matrix()
        datapoint_distances = []
        
        im_mapping_obj = InterpolationMapping(self.orginal_molecule)
        cartesian_distances_dict = {i:[] for i in range(len(distorted_molecule_list))}

        internal_coord_grouping = self.get_internal_coordinate_atom_type_mapping(self.orginal_molecule, atom_types)
                
        element_to_category = {item: category for category, items in internal_coord_grouping.items() for item in items}
        
        for j, dist_mol in enumerate(distorted_molecule_list):

            internal_distance_grouping = {}
            internal_coordinate_reference_groups = {}

            for i, ref_mol in enumerate(reference_molecules):

                impes_coordinate = ImpesCoordinates(z_mat)
                
                impes_coordinate.use_inverse_bond_length = False
                impes_coordinate.cartesian_coordinates = ref_mol.get_coordinates_in_bohr()
                impes_coordinate.define_internal_coordinates() 
                impes_coordinate.compute_internal_coordinates_values()

                reference_internal_coordinate_values = impes_coordinate.internal_coordinates_values
                internal_difference_vectors = []
    

                for index, element in enumerate(z_mat):

                    current_category = element_to_category.get(element)
                    
                    if current_category not in internal_distance_grouping:
                        
                        internal_distance_grouping[current_category] = [[] for i in range(len(reference_molecules))]
                
                mapped_structure, distance_norm, _ = im_mapping_obj.cost_matrix_analysis(dist_mol.get_coordinates_in_bohr(), ref_mol.get_coordinates_in_bohr())
                cartesian_distances_dict[j].append(distance_norm)
                impes_coordinate.cartesian_coordinates = mapped_structure
                impes_coordinate.define_internal_coordinates()
                impes_coordinate.compute_internal_coordinates_values()

                internal_difference_vectors.append(reference_internal_coordinate_values - impes_coordinate.internal_coordinates_values)
                for index, element in enumerate(z_mat):
                    current_category = element_to_category.get(element)
                    internal_distance_grouping[current_category][i].append(internal_difference_vectors[-1][index])  
                    internal_coordinate_reference_groups[current_category] = reference_internal_coordinate_values[index]
                
                # for index, element in enumerate(z_mat):
                #     current_category = element_to_category.get(element)
                    # internal_distance_grouping[current_category][i] = np.linalg.norm(internal_distance_grouping[current_category][i])

            datapoint_distances.append(internal_distance_grouping)
        return internal_coord_grouping, datapoint_distances, cartesian_distances_dict


    # The Greddy algorithm is defined for the problem at hand in the greedy_select_reference_points_combination() and works as follows:
    # 1. a while loop is established which is running as long as the size of the set with the best combination of labels is smaller then the size of the general label set
    # 2. A for loop increases the current best combination of labels with 1 more label in each iteration -> calculating the energy and checking if this combination has to be established as the new best combination of labels
    # 3. After the loop the solution is compared to the original interpolation energy and checked if an impovment has occured -> if nothing has improved and the energy is lower then the org value the search is finished and the while loop breaks
    # 4. if the performance has improved the while loop starts from the begining and the new best combination paired with new labels unless the existing best combination is below a certain threshold

    # This algorithm makes sure that a minimum is found even though it might not be the gloabal minimum rather a local one.


    def greedy_select_reference_points_combination(
        self,
        molecule,
        interpolation_settings,
        labels,
        distances,
        datafile,
        reference_value,
        interpolated_energy,
        non_core_symmetry_groups,
        must_include=None,
        tolerance=3e-1,
        alpha=1.0,
        beta=0.1
    ):
        progress_log = []
        
        # Copy the lists to avoid mutating the originals
        available_labels = labels[:]
        max_points = len(labels)

        # Initialize best_combination and best_distances
        best_combination = []
        best_distances = []
        best_weights = None
        best_difference = float('inf')
        start = False

        # If there are labels that must be included, add them right away
        if must_include is None:
            must_include = []
        
        for forced_label in must_include:
            if forced_label in available_labels:
                best_combination.append(forced_label)
                # Find and remove the corresponding distance too
                idx = available_labels.index(forced_label)
                forced_label_distance = distances[forced_label]
                best_distances.append(forced_label_distance)
                
                del available_labels[idx]
            else:
                # Optionally handle the case where forced_label wasn't found
                pass

        # Now your best_combination already includes the must_include labels.
        # Proceed with the normal greedy logic for the *other* points:
        
        print('Here is best combination and max_points', best_combination, len(best_combination), max_points)

        if len(best_combination) == max_points:
            start = True
            progress_log.append({
                "step": len(best_combination),
                "difference": abs(interpolated_energy - reference_value),
                "distance": best_distances,
                'weights': [1.0 for i in range(len(best_combination))],
                "combination": best_combination.copy(),
            })
        
        while len(best_combination) < max_points and start == False:
            improvement_found = False
            best_new_label = None
            best_new_diff = best_difference
            best_new_score = float('inf')

            for i, label in enumerate(available_labels):
                candidate_combo = best_combination + [label]

                # Calculate the interpolated energy for the candidate combo
                im_energies_kj_mol, important_weights, _ = self.calc_interpolation_energies_combination(
                    [molecule],
                    interpolation_settings,
                    candidate_combo,
                    datafile,
                    non_core_symmetry_groups
                )
                error = abs(im_energies_kj_mol[0] - reference_value)
                distance = distances[label]

                # Let the distance influence the current score
                score = alpha * error - beta * distance

                if score < best_new_score + 0.01:
                    best_new_score = score
                    best_new_diff = error
                    best_new_label = label
                    best_new_distance = distance
                    best_weights = important_weights
                    improvement_found = True

            # If we found no improvement, but current difference is better than baseline,
            # we can break early
            if not improvement_found and best_new_diff < abs(reference_value - interpolated_energy):
                break

            # Append the chosen label to the combination
            best_combination.append(best_new_label)
            best_difference = best_new_diff
            best_distances.append(best_new_distance)
            best_weights = best_weights[:len(best_distances)]
            # Remove from the pool
            index_to_remove = available_labels.index(best_new_label)
            del available_labels[index_to_remove]

            # Log the progress
            progress_log.append({
                "step": len(best_combination),
                "added_label": best_new_label,
                "difference": best_difference,
                "distance": best_distances,
                'weights': best_weights,
                "combination": best_combination.copy(),
                "score": best_new_score
            })

            if best_difference < tolerance:
                break

        return progress_log



    # find_the_best_combination() sets up the necessary variables for the Greddy algorithm and extracts the best solutin from the progress_log of the search.

    def find_best_combination(self, molecule, reference_value, interpolated_energy, interpolation_settings, datafile, distance_dict, non_core_symmetry_groups, labels_sorted=None, must_include=None):
        
        z_mat = molecule.get_z_matrix()
        
        impes_driver = ImpesDriver(z_mat)

        impes_driver.checkpoint_file_name = datafile
        if labels_sorted is None:
            labels = impes_driver.read_labels() 
            labels_sorted = sorted(labels, key=lambda x: int(x.split('_')[1]))

        progress_log = self.greedy_select_reference_points_combination(molecule, interpolation_settings, labels_sorted[:], distance_dict, datafile, reference_value, interpolated_energy, non_core_symmetry_groups, must_include)

        best_entry = min(progress_log, key=lambda x: x["difference"])

        # Extract the best_difference and combination
        best_difference = best_entry["difference"]
        best_combination = best_entry["combination"]
        best_distance = best_entry['distance']
        best_weights = best_entry['weights']

        return best_combination, best_difference, best_distance, best_weights

    # The main function of this cell is the data_frame_creator() which is used to generate the data frame in pandas from the train-, calibaration-, testset are defined. For that different distance features are used:
    # 1. The distances between the current structure X adn the reference structures X_i
    # 2. Additionally the distance of identical bond-, angle- and dihedral types between the current structure X and the reference structures X_i (which have been determined from the atom types itself by define_internal_coordinates_types())
    # -> It is possible that there are multiple bonds, angle, dihedrals of the same type... For those the mean is used as a value for that internal coordinate
    # 3. The last columns are the labels of the reference poionts that the model has to predict, which can take binary values (0 -> not to include or 1 -> to include)

    # Due to the large number of randomly selected structures and the time consuming Greedy algorithm the process for identifying the best combination of reference points has been parallelised using concurrent.features.
    # For that the process_batch() and process_single_structure(). The whole structure set is first divided into batches for which "workers" are assigned that perform calculations within a batch parallel.

    # The rsulting dataframe is used for the next steps

    def process_single_structure(self, index, key, structural_differences_dict,
                                strucutre_coordinats, qm_energies, interpolation_energies,
                                impes_dictionary, database_file, internal_distances_list,
                                labels_sorted, feature_labels, non_core_symmetry_groups):
        
        try:
            current_molecule = strucutre_coordinats[index]

            best_combination, best_difference, best_distances, best_weights = self.find_best_combination(
                current_molecule,
                qm_energies[index],
                interpolation_energies[index],
                impes_dictionary,
                database_file,
                structural_differences_dict[key],
                non_core_symmetry_groups
            )
            
            print('INDEX', index, best_combination, best_difference, abs(qm_energies[index] - interpolation_energies[index]))

            internal_distances = internal_distances_list[index]
            current_distances = []
            for key_internal, entries in internal_distances.items():
                for entry in entries:
                    current_distances.append(entry)

            combination_binary = [1 if label in best_combination else 0 for label in labels_sorted]

            row_data = []
            # row_data.extend(structural_differences_dict[key])
            row_data.extend(current_distances)
            row_data.extend(combination_binary)

            
            return index, row_data
    
        except Exception as e:
            print(f"Error in process_single_structure: {e}")
            raise 

    def process_batch(self, batch, structural_differences_dict, strucutre_coordinats, qm_energies, 
                    interpolation_energies, impes_dictionary, database_file, 
                    internal_distances_list, labels_sorted, feature_labels, non_core_symmetry_groups):
        """
        Process a batch of structures.
        Returns a list of (index, row_data) tuples.
        """
        results = []
        for index, key in batch:
            try:
                result = self.process_single_structure(
                    index,
                    key,
                    structural_differences_dict,
                    strucutre_coordinats,
                    qm_energies,
                    interpolation_energies,
                    impes_dictionary,
                    database_file,
                    internal_distances_list,
                    labels_sorted,
                    feature_labels,
                    non_core_symmetry_groups
                )
                results.append(result)
                
            except Exception as e:
                print(f"Error at idx={index}, key={key}: {e}")
                raise  # or handle

        return results

    def data_frame_creator_parallel(self, strucutre_coordinats, structural_differences_dict, internal_distances_labels, internal_distances_list, number_of_dp, labels_sorted, interpolation_energies, qm_energies, impes_dictionary, database_file, non_core_symmetry_groups):

        
        # dist_feature_labels = [f'dist_ref_{index}' for index in range(len(datapoint_energies))] 
        # energy_feature_labels = [f'energy_ref_{index}' for index in range(len(datapoint_energies))] 
        # gradient_feature_labels = [f'gradient_ref_{index}' for index in range(len(datapoint_energies))] 
        # hessian_feature_labels_0 = [f'hessian_min_{i}' for i in range(len(datapoint_energies))] 
        # hessian_feature_labels_1 = [f'hessian_max_{i}' for i in range(len(datapoint_energies))] 
        # hessian_feature_labels_2 = [f'hessian_trace_{i}' for i in range(len(datapoint_energies))] 
        internal_distance_features = [f'{key}_dp_{i}' for i in range(number_of_dp) for key in internal_distances_labels.keys()]
        
        print(internal_distance_features)
        # org_interpolation_energy_feature = ['org_interpolation_energy']
        # best_combination_interpolation_energy_feature = ['best_comb_interpolation_energy']

        
        
        feature_labels = internal_distance_features + labels_sorted #+ energy_feature_labels + gradient_feature_labels + hessian_feature_labels_0 + hessian_feature_labels_1 + hessian_feature_labels_2 + sorted_labels


        max_workers = 2
        items_to_process = list(enumerate(structural_differences_dict.keys()))
        batch_size = ceil(len(items_to_process) / max_workers)  # Adjust size dynamically
        batches = [items_to_process[i:i + batch_size] for i in range(0, len(items_to_process), batch_size)]
        # batch_size = 8  # !!!!!!!!! This value needs to be carefulluy choosen depending on the sieze of the data and the number of CPUs included because otherwise one can get problems with communication and splitting the data!!!!!!!!!!!!!!!
        # batches = [items_to_process[i:i+batch_size] for i in range(0, len(items_to_process), batch_size)]

        
        print('batches', len(batches), len(batches[0]))
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all batches as tasks
            futures = [
                executor.submit(self.process_batch, batch, structural_differences_dict, strucutre_coordinats, qm_energies, interpolation_energies, impes_dictionary, database_file, internal_distances_list, labels_sorted, feature_labels, non_core_symmetry_groups) for batch in batches ]

            # Collect the results
            final_dataframe = {}
            for future in concurrent.futures.as_completed(futures):
                print('Iam in future concurent result')
                batch_results = future.result()
                for idx, row_data in batch_results:
                    final_dataframe[idx] = row_data
                print('After the funal loop')
        
        df = pd.DataFrame.from_dict(final_dataframe, orient='index', columns=feature_labels)
        dataframe = df.sort_index()
        dataframe.to_csv('dataset.csv', index=False)
        return final_dataframe


    def database_extracter(self, datafile, molecule):
        
        im_driver = ImpesDriver(molecule.get_z_matrix()) # -> implemented Class in VeloxChem that is capable to perform interpolation calculations for a given molecule and provided z_matrix and database
        im_driver.checkpoint_file_name = datafile
        labels = im_driver.read_labels()
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        z_mat = molecule.get_z_matrix()
        mol_labels = molecule.get_labels()

        impes_coordinate = ImpesCoordinates(z_mat) # -> implemented Class in VeloxChem that handles all transformations and database changes concerning the interpolation
        data_point_molecules = []


        for label in sorted_labels:
            impes_coordinate.read_hdf5(datafile, label) # -> read in function from the ImpesDriver object
            coordinates_in_angstrom = impes_coordinate.cartesian_coordinates * bohr_in_angstrom()
            current_molecule = Molecule(mol_labels, coordinates_in_angstrom, 'angstrom') # -> creates a VeloxChem Molecule object
            data_point_molecules.append(current_molecule)

        return data_point_molecules, sorted_labels

    
    def compute_best_combination(self, current_molecule, sorted_labels, impes_dict, reference_calc_dict, distance_dict, non_core_symmetry_groups, must_include=None):


        datafile = impes_dict['checkpoint_file_name']
        im_energy, important_weights, _ = self.calc_interpolation_energies_combination([current_molecule], impes_dict, sorted_labels, datafile, non_core_symmetry_groups)
        reference_energy, _ = self.qm_energy_calculator(current_molecule, reference_calc_dict['driver'], reference_calc_dict['basis_set'], reference_calc_dict['xc_func'])
        
        atom_types = self.ff_gen.atom_types
        im_mapping_obj = InterpolationMapping(self.orginal_molecule)
        non_core_symmetry_groups = im_mapping_obj.non_core_symmetry_groups
        datapoint_molecules, sorted_labels = self.database_extracter(datafile, self.orginal_molecule)
        internal_coordinates_grouping, internal_distance_dict_list, distances_to_references = self.define_internal_coordinates_types(datapoint_molecules, [current_molecule], atom_types)
        
        print('internal distances:', internal_coordinates_grouping, '   ', internal_distance_dict_list, '   ', distances_to_references)

        best_combination, best_difference, best_distances, best_weights = self.find_best_combination(current_molecule, reference_energy[0], im_energy[0], impes_dict, datafile, distance_dict, non_core_symmetry_groups, sorted_labels, must_include)

        return best_combination, best_difference, best_distances, best_weights, im_energy[0], reference_energy[0]


    def compute_dataframe(self, dynamics_molecules, impes_dict, reference_calc_dict, number_of_labels):
        
        forcefield_generator = ForceFieldGenerator()
        forcefield_generator.force_field_data = self.ff_datafile
        forcefield_generator.partial_charges = self.orginal_molecule.get_partial_charges(self.orginal_molecule.get_charge())
        forcefield_generator.create_topology(self.orginal_molecule)
        atom_types = forcefield_generator.atom_types
        datafile = impes_dict['checkpoint_file_name']
        im_mapping_obj = InterpolationMapping(self.orginal_molecule)
        non_core_symmetry_groups = im_mapping_obj.non_core_symmetry_groups
        datapoint_molecules, sorted_labels = self.database_extracter(datafile, self.orginal_molecule)
        internal_coordinates_grouping, internal_distance_dict_list, distances_to_references = self.define_internal_coordinates_types(datapoint_molecules, dynamics_molecules, atom_types)

        sorted_labels = sorted_labels[:number_of_labels]
        
        reference_energies = []
        im_energies = []
        for mol in dynamics_molecules:
            
            im_energy, _, _ = self.calc_interpolation_energies_combination([mol], impes_dict, sorted_labels, datafile, non_core_symmetry_groups)
            reference_energy, _ = self.qm_energy_calculator(mol, reference_calc_dict['driver'], reference_calc_dict['basis_set'], reference_calc_dict['xc_func'])
            reference_energies.append(reference_energy[0])
            im_energies.append(im_energy[0])

        dataframe = self.data_frame_creator_parallel(dynamics_molecules, distances_to_references, internal_coordinates_grouping, internal_distance_dict_list, len(datapoint_molecules), sorted_labels, im_energies, reference_energies, impes_dict, datafile, non_core_symmetry_groups)

        return dataframe

    def qm_energy_calculator(self, mol, QM_Driver='XTB', basis_set_label='sto-3g', xc_func='b3lyp', gradient=False):

        qm_energies_kcal_mol = []
        gradients = []
        if QM_Driver == 'XTB':
            
            xtb_drv = XtbDriver()
            xtb_grad_drv = XtbGradientDriver(xtb_drv)

                
            xtb_drv.compute(mol)
            qm_energies_kcal_mol.append(xtb_drv.get_energy() * 627.509474)
            if gradient is True:
                xtb_grad_drv.compute(mol)
                gradients.append(xtb_grad_drv.get_gradient())
        

        ### This section uses more elaborate and expensive methods and have not been used further due to time constained
        if QM_Driver == 'Veloxchem':

            vlx_driver = ScfRestrictedDriver() # -> QM - Driver implemented within Veloxchem to calculate the energy for a given vlx.Molecule based on the self-consistent-field formalism
            vlx_driver.xcfun = xc_func # -> exchange-correlation functionals stored within veloxchem 
            vlx_gradient_driver = ScfGradientDriver(vlx_driver) # -> QM - GradientDriver implemented within Veloxchem to calculate the gradient for a given vlx.Molecule based on the self-consistent-field formalism
            vlx_driver.ostream.mute()
            basis = MolecularBasis.read(mol, basis_set_label,) # -> neccessary basis set for making the SCF calculations possible
            print('Basis', basis) 
            scf_result = vlx_driver.compute(mol, basis)
            qm_energies_kcal_mol.append(vlx_driver.scf_energy * 627.509474)
                
            if gradient is True:

                vlx_gradient_driver.compute(scf_result)
                gradients.append(vlx_gradient_driver.gradient)

        
        if QM_Driver == 'Orca':
            
            qm_driver = ExternalQMDriver('ORCA', 0, 0, 1, 16, hostname=False) # ->  ExternalDriver implemented within Veloxchem to communicate with another QM-program ORCA to calculate the energy for a given vlx.Molecule based on the self-consistent-field formalism
            qm_driver.spin_flip = False
            qm_driver.basis_set_label = basis_set_label
            qm_driver.xc_func = xc_func
            qm_grad_driver = ExternalQMGradientDriver(qm_driver) # ->  ExternalGradientDriver implemented within Veloxchem to communicate with another QM-program ORCA to calculate the gradient for a given vlx.Molecule based on the self-consistent-field formalism
            
            energy = qm_driver.compute(mol)
            qm_energies_kcal_mol.append(energy[0] * 627.509474)
            if gradient is True:
                qm_grad_driver.compute(mol)
                gradient = qm_grad_driver.extract_gradients()
                gradients.append(gradient[0])
        
        return qm_energies_kcal_mol, gradients
    
    # calc_interpolation_energies_combination() is used to claculate the interpolated energy for a given combination of labels -> used for the combination given by the Greedy algorithm

    def calc_interpolation_energies_combination(self, molecule_list, interpolation_settings, labels, datafile, non_core_symmetry_group, calc_gradient=False):
        
        print('labels', labels)
        impes_driver = ImpesDriver(self.orginal_molecule.get_z_matrix())
        impes_driver.update_settings(interpolation_settings)
        impes_driver.checkpoint_file_name = datafile
        impes_driver.distance_thrsh = 100.0
        impes_driver.non_core_symmetry_group = non_core_symmetry_group

        im_energies_kj_mol = []
        gradients = []

        for mol in molecule_list:
            impes_driver.compute(mol, None, None, labels=labels)
            important_weights = impes_driver.individual_weights
            im_energies_kj_mol.append(impes_driver.impes_coordinate.energy * 627.509474) # transofrmation in kcal/mol

            if calc_gradient is True:

                gradients.append(impes_driver.impes_coordinate.gradient)
        
        return im_energies_kj_mol, important_weights, gradients