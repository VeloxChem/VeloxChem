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

from .molecularbasis import MolecularBasis
from .imdatabasedriver import IMDatabaseDriver
#from .impesforcefieldgenerator_parallel import ImpesForceFieldGeneratorParallel
from .forcefieldgenerator import ForceFieldGenerator
from .openmmdynamics import OpenMMDynamics

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

class DatabaseDriver:
    """
    Class to set up and control the construction of the Interpolation Dynamics database.
    
    Instances variables:

        - density_of_datapoints: variable that keeps track of the number of datapoints currently in the database.
        - qm_data_points: database point object with cartesian_coordinates, internal_coordinate_values, qm_energy,
          internal gradient and internal hessian.
        - qmlabels: label of the each database point.
        - molecules_along_rp: (currently) moleculear structures along a given internal coordinate.
        - impes_dict: dictionary map with settings for the IMPES driver object.
        - molecule: initital molecule (defined by the user!).
        - datafile: the database.
        - angle_threshold: defines a range in which a current structure stays within the system (only used if dihedral for rotation is given)
        - z_matrix: original list of defined internal coordinates (bonds, angles and dihedrals),
        - parallel: (not really implemented attempt to allow a split of the database calculation)
        - add_bias_force: steering the dynmacis in the database construction to a certain internal coordinate if defined  
    """
    
    def __init__(self):
        """
        Initializes the DATABASE driver.
        """
        
        #self.impes_forcefield_generator = TestNewImpesFFGen()

        self.density_of_datapoints = None
        self.qm_data_points = None
        self.qmlabels = None
        self.molecules_along_rp = None
        self.impes_dict = None
        self.reference_calc_dict = None
        self.molecule = None
        self.datafile = None
        self.dihedrals = None
        self.distance_threshold = None
        self.z_matrix = None
        self.paralell = False
        self.add_bias_force = False

        # individual run
        self.qm_energies = []
        self.total_energies = []
        self.molecules = []
        self.kinetic_energies = []
        self.point_added_molecules = []
        self.unique_molecules = []

    def set_up_the_system(self, molecule, impes_dict, dynamics_settings, target_dihedral=None, sampling_structures=1):

        """
        Assign the neccessary variables with respected values. 

        :param molecule: original molecule

        :param impes_dict: user defined dictionary for the interpolation cacluations (see IMPES driver)

        :param dynamics_settings: sets all important parameters for the dynamics:

            
        """


        self.qm_data_points = None
        self.molecule = molecule
        self.impes_dict = impes_dict
        if dynamics_settings:
            self.dynamics_settings = dynamics_settings
        if self.z_matrix is None:
            self.z_matrix = molecule.get_z_matrix()
        if 'checkpoint_file_name' in self.impes_dict:
            data_base_file = self.impes_dict['checkpoint_file_name']
            print('checkpoint file name', self.impes_dict, data_base_file)
            impes_driver = ImpesDriver(self.z_matrix)
            impes_driver.update_settings(self.impes_dict)
            self.qmlabels = impes_driver.read_labels()

            self.qm_data_points = []
            self.qm_energies = []
            for label in self.qmlabels:
                qm_data_point = ImpesCoordinates(self.z_matrix)
                qm_data_point.read_hdf5(data_base_file, label)
                self.qm_energies.append(qm_data_point.energy)
                self.qm_data_points.append(qm_data_point)        
        
        if target_dihedral:
            self.density_of_datapoints, self.molecules_along_rp, self.distance_threshold = self.database_density_check_with_molecule(molecule, self.qm_data_points, target_dihedral, nsampling=sampling_structures)
            self.dihedrals = target_dihedral
            self.impes_dict['dihedrals'] = target_dihedral
        elif self.molecules_along_rp is None:
            self.density_of_datapoints, self.molecules_along_rp, self.distance_threshold = self.database_density_check_with_molecule(molecule, self.qm_data_points, rotate=False, distance_threshold=np.inf)
        
        print('here is the density', self.density_of_datapoints)

    def run_data_base_construction(self):
        
        for counter, mol in enumerate(self.molecules_along_rp):

            #structure = structure * bohr_in_angstrom()
            #current_molecule = Molecule(self.molecule.get_labels(), structure, units="angstrom")
            forcefield_generator = ForceFieldGenerator()
            forcefield_generator.force_field_data = self.dynamics_settings['FF_datafile']
            scf_result = 'vlx_20241016_097c295c.scf.h5 '
            basis = MolecularBasis.read(mol,
                                        '6-31G*',
                                        ostream=None)

            self.dynamics_settings['trajectory_file'] = f'trajectory_{counter}.pdb'
            forcefield_generator.partial_charges = mol.get_partial_charges(mol.get_charge())
            forcefield_generator.create_topology(mol)

            # Identify the core that needs to be mapped and identify symmetric groups

            core_structure = self.define_core_system_for_mapping(self.molecule)
            non_core_symmetric_groups_map, non_core_symmetric_groups_combined = self.symmetric_groups(self.molecule, forcefield_generator)
            mapping_of_elements = self.map_swapping_elements(non_core_symmetric_groups_map)
            print('non_core_symmetric', non_core_symmetric_groups_map, non_core_symmetric_groups_combined, '\n\n', mapping_of_elements)
            self.impes_dict['symmetry_groups'] = non_core_symmetric_groups_combined
            different_bins = False
            if counter > 0:
                relaxed = True
                if different_bins:
                    self.impes_dict['basename'] = f'im_collcection_run{counter}'
                else:
                    if 'checkpoint_file_name' not in self.impes_dict:
                        self.impes_dict['checkpoint_file_name'] = f'{self.impes_dict['basename']}.h5'
                    
                    self.datafile = self.impes_dict['checkpoint_file_name']
                    impes_driver = ImpesDriver(self.z_matrix)
                    impes_driver.update_settings(self.impes_dict)
                    self.qmlabels = impes_driver.read_labels()

                    self.qm_data_points = []
                    self.qm_energies = []
                    for label in self.qmlabels:
                        qm_data_point = ImpesCoordinates(self.z_matrix)
                        qm_data_point.read_hdf5(self.datafile, label)
                        self.qm_energies.append(qm_data_point.energy)
                        self.qm_data_points.append(qm_data_point)
                    new_density_around_data_points, _, _ =  self.database_density_check_with_molecule(mol, self.qm_data_points, rotate=False, distance_threshold=self.distance_threshold)
                    print('here is the denisty after first run', new_density_around_data_points)
                    self.density_of_datapoints[counter] = new_density_around_data_points[-1]
                
            im_database_driver = IMDatabaseDriver()
            # if len(self.molecule.get_labels()) > 25 and self.paralell == True:
            #     im_database_driver = ImpesForceFieldGeneratorParallel()            IMDatabaseDriver
            # im_database_driver.core_structure = core_structure
            # im_database_driver.non_core_symmetry_groups = non_core_symmetric_groups
            im_database_driver.platform = 'CUDA'
            solvent_settings = 'gas'
            if 'solvent' in self.dynamics_settings:
                solvent_settings = self.dynamics_settings['solvent']
            im_database_driver.system_from_molecule(mol, self.z_matrix, forcefield_generator, solvent=solvent_settings, qm_atoms='all')
            if self.add_bias_force is True:
                im_database_driver.add_bias_force([4, 5, 6, 7], 40, -30)
            nroots = int(self.dynamics_settings['roots']) + 1
            desiered_point_density = int(self.dynamics_settings['desired_datapoint_density'])
            im_database_driver.update_settings(self.dynamics_settings, self.impes_dict)
            
            print('here is the seccond density', self.density_of_datapoints, self.density_of_datapoints[counter], desiered_point_density)
            if self.density_of_datapoints[counter] / nroots < desiered_point_density:   
                #density_around_data_points = im_database_driver.run_qmmm(self.density_of_datapoints[counter], desiered_point_density, self.distance_threshold, unadded_cycles, nroots, calc_NAC, checkpoint, ensemble, nsteps=nsteps, temperature=temperature, timestep=timestep, snapshots=snapshots, out_file=f'{counter}_run_{out}')            
                im_database_driver.density_around_data_point = self.density_of_datapoints[counter]
                im_database_driver.distance_threshold = self.distance_threshold
                im_database_driver.reference_calc_dict = self.reference_calc_dict
                im_database_driver.run_qmmm()
                self.density_of_datapoints[counter] = im_database_driver.density_around_data_point

                # individual impes run objects
                self.qm_energies.append(im_database_driver.qm_potentials)
                self.total_energies.append(im_database_driver.total_energies)
                self.kinetic_energies.append(im_database_driver.kinetic_energies)
                self.molecules.append(im_database_driver.molecules)
                self.point_added_molecules.append(im_database_driver.point_adding_molecule)
                self.unique_molecules.append(im_database_driver.allowed_molecules)


            counter += 1
        
        print('The construction of the database was sucessfull')    

    
    def database_density_check_with_molecule(self, molecule, qm_datapoints, specific_dihedral=None, rotate=True, nsampling=None, distance_threshold=1.0):
        """
        loads the molecule and specific features (dihedral) that is changed in order to sample different molecular conformations.

        :param filename: 
            Path to the PDB file. Example: 'system.pdb'.
        """
        sampled_molecules = []
        if rotate == True:
            rotation_values = np.linspace(0, 360, nsampling)
            print(rotation_values, specific_dihedral)
            for theta in rotation_values:
                molecule.set_dihedral_in_degrees(specific_dihedral, theta)
                new_molecule = Molecule(molecule.get_labels(), molecule.get_coordinates_in_bohr(), 'bohr')
                sampled_molecules.append(new_molecule)
        
        else:
            sampled_molecules.append(molecule)

        # check for the point densities

        point_densities = []
        print('qm data_points', qm_datapoints)
        if qm_datapoints:
            print('I am going into the qm_dataPoint')
            impes_driver = ImpesDriver(self.z_matrix)
            impes_driver.update_settings(self.impes_dict)
            for i in sampled_molecules:
                impes_driver.define_impes_coordinate(i.get_coordinates_in_bohr())
                density = 0
                for point in qm_datapoints:

                    distance, _, _ = impes_driver.cartesian_distance(point)
                    if distance * bohr_in_angstrom() <= distance_threshold:
                        density += 1
                point_densities.append(density)
        else:
            for i in sampled_molecules:
                density = 0
                point_densities.append(density)

        return point_densities, sampled_molecules, distance_threshold
    
    
    def define_core_system_for_mapping(self, molecule):

        
        connectivity = molecule.get_connectivity_matrix()
        
        mol_labels = molecule.get_labels()

        atoms_to_exclude = []
        for atom, mol_label in enumerate(mol_labels):

            connections = np.sum(connectivity[atom, :])

            if connections == 1 and mol_label == 'H':

                atoms_to_exclude.append(atom)
        
        core_structure = [idx for idx, _ in enumerate(mol_labels) if idx not in atoms_to_exclude]

        return core_structure
    
    def symmetric_groups(self, molecule, ff_gen):
        
        connectivity_matrix = molecule.get_connectivity_matrix()
        rotatable_bonds = ff_gen.rotatable_bonds
        dihedral_side_atom_types = {}
        atom_types = ff_gen.atom_types
        

        print('rotatable bonds', rotatable_bonds)

        for rot_bond in rotatable_bonds:
            
            print('current_rot_bond', rot_bond)
            atom_1_connections = [(index, atom_types[index]) for index, connected in enumerate(connectivity_matrix[rot_bond[0] - 1]) if connected == 1 and index != rot_bond[1] - 1]
            atom_2_connections = [(index, atom_types[index])  for index, connected in enumerate(connectivity_matrix[rot_bond[1] - 1]) if connected == 1 and index != rot_bond[0] - 1]
            print(atom_1_connections, '\n', atom_2_connections)

            if rot_bond[0] - 1  not in dihedral_side_atom_types:

                dihedral_side_atom_types[rot_bond[0]-1] = atom_1_connections
                    
            if rot_bond[1] - 1 not in dihedral_side_atom_types:

                dihedral_side_atom_types[rot_bond[1]-1] = atom_2_connections
            
        index_dict = {}
        print('Here', dihedral_side_atom_types.items())
        for keys, diehdral_group in dihedral_side_atom_types.items():

            index_dict[keys] = {}

            for tuple in diehdral_group:

                if tuple[1] not in index_dict[keys]:
                    index_dict[keys][tuple[1]] = []
                
                index_dict[keys][tuple[1]].append(tuple[0])

        print(index_dict)

        non_core_symmetric_groups = {}
        non_core_symmetric_groups_combinded = []
        for key, groups in index_dict.items():

            for atom_type, atom_group in groups.items():
                
                if len(atom_group) <= 2:
                    continue
                
                if key not in non_core_symmetric_groups:
                     non_core_symmetric_groups[key] = atom_group
                                    
                #non_core_symmetric_groups[key].append(atom_group)
                non_core_symmetric_groups_combinded += atom_group
        

        return non_core_symmetric_groups, non_core_symmetric_groups_combinded
    
    def map_swapping_elements(self, non_core_symmetry_groups_map):

    # find element of interest in the z_matrix

        mapping_between_elements = {}
        
        for key, entries in non_core_symmetry_groups_map.items():
            
            possible_indices_perm = list(itertools.permutations(entries))

            print(possible_indices_perm)
            
            mapping_between_elements.update({current_perm:{org_index: new for new, org_index in zip(entries, current_perm)} for current_perm in possible_indices_perm})
        

        return mapping_between_elements
    
    # database_extracter is a function that that extracts from a given interpolation database the information stored (molecule, energy, gradient, hessian)

    def database_extracter(self, datafile, molecule):
        
        im_driver = ImpesDriver(molecule.get_z_matrix()) # -> implemented Class in VeloxChem that is capable to perform interpolation calculations for a given molecule and provided z_matrix and database
        im_driver.checkpoint_file_name = datafile
        labels = im_driver.read_labels()
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        z_mat = molecule.get_z_matrix()
        mol_labels = molecule.get_labels()

        impes_coordinate = ImpesCoordinates(z_mat) # -> implemented Class in VeloxChem that handles all transformations and database changes concerning the interpolation
        data_point_molecules = []
        data_point_energy = []
        data_point_gradient = []
        data_point_hessian = {}

        data_point_hessian['min'] = []
        data_point_hessian['max'] = []
        data_point_hessian['trace'] = []


        for label in sorted_labels:
            impes_coordinate.read_hdf5(datafile, label) # -> read in function from the ImpesDriver object
            coordinates_in_angstrom = impes_coordinate.cartesian_coordinates * bohr_in_angstrom()
            current_molecule = Molecule(mol_labels, coordinates_in_angstrom, 'angstrom') # -> creates a VeloxChem Molecule object

            eigenvalues, _ = np.linalg.eigh(impes_coordinate.internal_hessian)

            min_eigenvalue = np.min(eigenvalues)
            max_eigenvalue = np.max(eigenvalues)

            trace_hess = np.trace(impes_coordinate.internal_hessian)
            
            data_point_molecules.append(current_molecule)
            data_point_energy.append(float(impes_coordinate.energy))
            data_point_gradient.append(float(np.linalg.norm(impes_coordinate.internal_gradient)))
            data_point_hessian['min'].append(min_eigenvalue)
            data_point_hessian['max'].append(max_eigenvalue)
            data_point_hessian['trace'].append(trace_hess)



        return data_point_molecules, data_point_energy, data_point_gradient, data_point_hessian
    
    def datapoint_distribution(self, datapoint_molecules_1, molecule, ff_gen=None, datapoint_molecules_2=None):
        
        distances_from_ref_to_ref_1 = {}
        distances_from_ref_to_ref_2 = {}
        
        _, non_core_symmetric_groups_combined = self.symmetric_groups(molecule, ff_gen)
        
        for i in range(len(datapoint_molecules_1)):
            distances_from_ref_to_ref_1[i] = [] 
            for j in range(0, len(datapoint_molecules_1)):
                if i == j:
                    continue     
            
                _, distance_norm_1 = self.cost_matrix_analysis(datapoint_molecules_1[i].get_coordinates_in_bohr(), datapoint_molecules_1[j].get_coordinates_in_bohr(), non_core_symmetric_groups_combined, molecule)

                distances_from_ref_to_ref_1[i].append(distance_norm_1)

        if datapoint_molecules_2 is not None:
            for i in range(len(datapoint_molecules_2)):
                distances_from_ref_to_ref_2[i] = [] 
                for j in range(0, len(datapoint_molecules_2)):
                    if i == j:
                        continue       

                    _, distance_norm_2 = self.cost_matrix_analysis(datapoint_molecules_2[i].get_coordinates_in_bohr(), datapoint_molecules_2[j].get_coordinates_in_bohr(), non_core_symmetric_groups_combined, molecule)

                    distances_from_ref_to_ref_2[i].append(distance_norm_2)


        return distances_from_ref_to_ref_1, distances_from_ref_to_ref_2

    def binning_molecules_into_distance(self, molecules, datapoint_molecules, ff_gen=None, start=0, input=30):
        
        # Define the bins as tuples
        bin_system = np.linspace(start, input, input * 4 + 1)

        # Generate bins dynamically
        bins = [(bin_system[i], bin_system[i + 1]) for i in range(0, len(bin_system) - 1)]

        # Add the last bin with a very large upper limit
        bins[-1] = (bins[-1][0], np.inf)
        bins[0] = (-np.inf, bins[0][1])

        # Initialize the dictionary with empty lists for each bin

        # Function to find the bin for a value
        def get_bin(value, bins):
            for bin_range in bins:
                if bin_range[0] <= value < bin_range[1]:  # Check if value falls in the range
                    return bin_range
            return None  # If the value doesn't fit any bin

        non_core_symmetric_groups_combined = [] 
        if ff_gen is not None:

            _, non_core_symmetric_groups_combined = self.symmetric_groups(molecules[0], ff_gen)
        
        dict_for_each_molecule = []
        for i in range(len(molecules)):
            bin_molecule_dict = {bin_range: [0, []] for bin_range in bins}
            for j in range(len(datapoint_molecules)):
                        
                _, distance_norm = self.cost_matrix_analysis(molecules[i].get_coordinates_in_bohr(), datapoint_molecules[j].get_coordinates_in_bohr(), non_core_symmetric_groups_combined, molecules[0])

                bin_range = get_bin(distance_norm, bins)

                if bin_range:
                    bin_molecule_dict[bin_range][0] += 1
                    bin_molecule_dict[bin_range][1].append(distance_norm)
            
            dict_for_each_molecule.append(bin_molecule_dict)

        return dict_for_each_molecule
    
    def simple_run_dynamics(self, molecule, forcefield_generator, dynamics_settings, scanning_method, interpolation_settings=None, reference_calc_dict=None):

        all_structures = None

        if scanning_method == 'MM':
            
            # define OpenMMDriver object and perform a dynamical sampling
            openmmdyn = OpenMMDynamics()
            openmmdyn.system_from_molecule(molecule, forcefield_generator)
            # openmmdyn.add_bias_force([21, 5, 2, 1], 2, -210)
            _, conformation_structures = openmmdyn.conformational_sampling(snapshots=dynamics_settings['snapshots'], nsteps=dynamics_settings['nsteps'], temperature=dynamics_settings['temperature'], minimize=False)
            all_structures = [Molecule.from_xyz_string(structure) for structure in conformation_structures]
            

        if scanning_method == 'IM':
            openmmdyn = OpenMMDynamics()
            openmmdyn.reference_calc_dict = reference_calc_dict
            openmmdyn.system_from_molecule(molecule, ff_gen=forcefield_generator, solvent=dynamics_settings['solvent'], qm_atoms='all')
            openmmdyn.platform = 'CUDA'
            openmmdyn.update_settings(dynamics_settings, interpolation_settings)
            openmmdyn.add_bias_force([21, 5, 2, 1], 2, -210)
            openmmdyn.run_qmmm()
            
            all_structures = openmmdyn.unique_molecules
        elif scanning_method == 'QM':
            openmmdyn = OpenMMDynamics()
            openmmdyn.system_from_molecule(molecule, forcefield_generator, solvent=dynamics_settings['solvent'], qm_atoms='all')
            openmmdyn.platform = 'CUDA'
            openmmdyn.update_settings(dynamics_settings, interpolation_settings)
            openmmdyn.run_qmmm()
            all_structures = openmmdyn.unique_molecules

        return all_structures


    def confirm_database_quality(self, 
                                 molecule, 
                                 interpolation_database, 
                                 interpolation_settings, 
                                 dynamics_settings, 
                                 reference_calculation_settings, 
                                 scanning_method='GF',
                                 given_structure_file=None,
                                 number_of_structures=10,
                                 threshold=2.5,
                                 extend=False,
                                 second_interpolation_database=None):

        # For all Methods a ForceField of the molecule is requiered
        forcefield_generator = ForceFieldGenerator()
        forcefield_generator.force_field_data = dynamics_settings['FF_datafile']
        forcefield_generator.partial_charges = molecule.get_partial_charges(molecule.get_charge())
        forcefield_generator.create_topology(molecule)
        
        database_quality = False 
        z_matrix = molecule.get_z_matrix()
        while database_quality is False:

            all_structures = self.simple_run_dynamics(molecule, forcefield_generator, dynamics_settings, scanning_method, interpolation_settings, reference_calculation_settings)

            # if scanning_method == 'MM':

                
            #     # define OpenMMDriver object and perform a dynamical sampling

            #     openmmdyn = OpenMMDynamics()
            #     openmmdyn.system_from_molecule(molecule, forcefield_generator)
            #     # openmmdyn.add_bias_force([21, 5, 2, 1], 2, -210)
            #     _, conformation_structures = openmmdyn.conformational_sampling(snapshots=dynamics_settings['snapshots'], nsteps=dynamics_settings['nsteps'], temperature=dynamics_settings['temperature'], minimize=False)

            #     all_structures = [Molecule.from_xyz_string(structure) for structure in conformation_structures]
                

    
            # if scanning_method == 'IM':

            #     openmmdyn = OpenMMDynamics()
            #     openmmdyn.system_from_molecule(molecule, ff_gen=forcefield_generator, solvent=dynamics_settings['solvent'], qm_atoms='all')
            #     openmmdyn.platform = 'CUDA'

            #     openmmdyn.update_settings(dynamics_settings, interpolation_settings)
            #     openmmdyn.add_bias_force([21, 5, 2, 1], 2, -210)
            #     openmmdyn.run_qmmm()
                
            #     all_structures = openmmdyn.unique_molecules

            # elif scanning_method == 'QM':

            #     openmmdyn = OpenMMDynamics()
            #     openmmdyn.system_from_molecule(molecule, forcefield_generator, solvent=dynamics_settings['solvent'], qm_atoms='all')
            #     openmmdyn.platform = 'CUDA'

            #     openmmdyn.update_settings(dynamics_settings, interpolation_settings)
            #     openmmdyn.run_qmmm()

            #     all_structures = openmmdyn.unique_molecules

            # elif scanning_method == 'GF':
                
            #     all_structures, _, _, _, _, _ =self.read_in_xyz_traj(given_structure_file)
                

            datapoint_molecules, _, data_points_gradient, data_points_hessian = self.database_extracter(interpolation_database, molecule)

            datapoint_molecules_2, data_points_gradient_2, data_points_hessian_2 = None, None, None
            if second_interpolation_database is not None:
                datapoint_molecules_2, _, data_points_gradient_2, data_points_hessian_2 = self.database_extracter(second_interpolation_database, molecule)

            distances_1, distances_2 = self.datapoint_distribution(datapoint_molecules, molecule, forcefield_generator, datapoint_molecules_2)
            
            distances_list_1_min = []
            for distances_from_point in distances_1.values():
                distances_within_range = [distance for distance in distances_from_point if abs(distance - min(distances_from_point)) < 0.2]
                if np.mean(distances_within_range) < 1e-5:
                    distances_list_1_min.append(0.05)
                else:
                    distances_list_1_min.append(np.mean(distances_within_range))
            
            distances_list_2_min = []
            for distances_from_point in distances_2.values():
                distances_within_range = [distance for distance in distances_from_point if abs(distance - min(distances_from_point)) < 0.2]
                if np.mean(distances_within_range) < 1e-5:
                    distances_list_2_min.append(1.0)
                else:
                    distances_list_2_min.append(np.mean(distances_within_range))
            
            print('here is the len', len(distances_list_1_min), len(distances_1), len(distances_list_2_min))
            non_core_symmetric_groups_map, non_core_symmetric_groups_combined = self.symmetric_groups(molecule, forcefield_generator)
            
            # mapping_of_elements = self.map_swapping_elements(non_core_symmetric_groups_map)

            rmsd = - np.inf
            random_structure_choices = None
            counter = 0
            if scanning_method == 'GF':
                random_structure_choices = all_structures
            while rmsd < 0.3 and scanning_method != 'GF' or counter > 20 and scanning_method != 'GF' :
                
                individual_distances = []
                random_structure_choices = random.choices(all_structures, k=number_of_structures)

                for datapoint_molecule in datapoint_molecules:

                    for random_struc in random_structure_choices:
                        
                        _, distance_norm = self.cost_matrix_analysis(random_struc.get_coordinates_in_bohr(), datapoint_molecule.get_coordinates_in_bohr(), non_core_symmetric_groups_combined, molecule)

                        individual_distances.append(distance_norm)
                
                rmsd = min(individual_distances)
                counter += 1

                if rmsd >= 0.3:
                    print(f'The overall RMSD is {rmsd} -> The current structures are well seperated from the database conformations! loop is discontinued')
                else:
                    print(f'The overall RMSD is {rmsd} -> The current structures are not all well seperated from the database conformations! loop is continued')        
            
            
            qm_energies = []
            im_energies = []
            im_energies_2 = []

            im_energies_no_mapping = []
            im_energies_2_no_mapping = []


            qm_driver = dynamics_settings['qm_driver']

            impes_driver = ImpesDriver(molecule.get_z_matrix())
            impes_driver.update_settings(interpolation_settings)
            impes_driver.checkpoint_file_name = interpolation_database
            impes_driver.distance_thrsh = 1.0
            impes_driver.data_points_gradient = distances_list_1_min
            impes_driver.data_points_hessian = data_points_hessian

            impes_driver.non_core_structure = non_core_symmetric_groups_combined
            impes_driver.non_core_symmetry_group = non_core_symmetric_groups_combined

            impes_driver_no_mapping = ImpesDriver(molecule.get_z_matrix())
            impes_driver_no_mapping.update_settings(interpolation_settings)
            impes_driver_no_mapping.checkpoint_file_name = interpolation_database
            impes_driver_no_mapping.distance_thrsh = 1.0

            impes_driver_no_mapping.non_core_structure = []
            impes_driver_no_mapping.non_core_symmetry_group = []

            labels = impes_driver.read_labels()
            sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

            for i, mol in enumerate(random_structure_choices):
                
                print('Angle of the dihedral', mol.get_dihedral_in_degrees([21, 5, 2, 1]))
                impes_driver.compute(mol, labels=sorted_labels)
                impes_driver_no_mapping.compute(mol, labels=sorted_labels)
                reference_energy, reference_gradient = self.qm_energy_calculator(mol, 
                                                                                QM_Driver=reference_calculation_settings['driver'], 
                                                                                basis_set_label=reference_calculation_settings['basis_set'],
                                                                                xc_func=reference_calculation_settings['xc_func'],
                                                                                gradient=reference_calculation_settings['gradient'])
                
                if second_interpolation_database is not None:
                    
                    impes_driver_2 = ImpesDriver(molecule.get_z_matrix())
                    impes_driver_2.update_settings(interpolation_settings)
                    impes_driver_2.checkpoint_file_name = second_interpolation_database
                    impes_driver_2.distance_thrsh = 1.0
                    # impes_driver_2.data_points_gradient = distances_list_2_min
                    # impes_driver_2.data_points_hessian = data_points_hessian_2

                    impes_driver_2.non_core_structure = non_core_symmetric_groups_combined
                    impes_driver_2.non_core_symmetry_group = non_core_symmetric_groups_combined
                    
                    impes_driver_2_no_mapping = ImpesDriver(molecule.get_z_matrix())
                    impes_driver_2_no_mapping.update_settings(interpolation_settings)
                    impes_driver_2_no_mapping.checkpoint_file_name = second_interpolation_database
                    impes_driver_2_no_mapping.distance_thrsh = 1.0

                    impes_driver_2_no_mapping.non_core_structure = []
                    impes_driver_2_no_mapping.non_core_symmetry_group = []

                    labels_2 = impes_driver_2.read_labels()
                    sorted_labels_2 = sorted(labels_2, key=lambda x: int(x.split('_')[1]))
                    impes_driver_2.compute(mol, labels=sorted_labels_2)
                    impes_driver_2_no_mapping.compute(mol, labels=sorted_labels_2)

                qm_energies.append(reference_energy[0])
                im_energies.append(impes_driver.impes_coordinate.energy * 627.509474)
                im_energies_2.append(impes_driver_2.impes_coordinate.energy * 627.509474)
                im_energies_no_mapping.append(impes_driver_no_mapping.impes_coordinate.energy * 627.509474)
                im_energies_2_no_mapping.append(impes_driver_2_no_mapping.impes_coordinate.energy * 627.509474)

                
                print('Energies', qm_energies[-1], im_energies[-1], im_energies_2[-1], im_energies_no_mapping[-1], im_energies_2_no_mapping[-1])
                
                print(f'\n\n ########## Step {i} ######### \n\n')
                if abs(qm_energies[-1] - im_energies[-1]) > threshold and extend is True:
                    
                    labels = []
                    labels.append("point_{0}".format((len(sorted_labels) + 1)))

                    external_qm_driver = ExternalQMDriver('ORCA', 0, 0, 1, 8, hostname=False)
                    external_grad_driver = ExternalQMGradientDriver(external_qm_driver)
                    external_hess_driver = ExternalQMHessianDriver(external_qm_driver)

                    self.add_point(mol, labels, interpolation_database, datapoint_molecules, interpolation_settings, z_matrix, external_qm_driver, external_grad_driver, external_hess_driver, non_core_symmetric_groups_combined)
                    break

            if len(qm_energies) == len(random_structure_choices):
                database_quality = True

        self.structures_to_xyz_file(random_structure_choices, 'random_xyz_traj.xyz', im_energies, qm_energies, im_energies_2, im_energies_no_mapping, im_energies_2_no_mapping)


        return random_structure_choices, qm_energies, im_energies, im_energies_2, im_energies_no_mapping, im_energies_2_no_mapping

        
    def calculate_translation_coordinates_analysis(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates

    def cost_matrix_analysis(self, current_coordinates, datapoint_coordinate, non_core_symmetry_groups, molecule):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                ImpesCoordinates object
        """
            
        # First, translate the cartesian coordinates to zero
        target_coordinates = self.calculate_translation_coordinates_analysis(datapoint_coordinate)
        reference_coordinates = self.calculate_translation_coordinates_analysis(current_coordinates)
        # remove the non_core atoms (currently just H)
        target_coordinates_core = np.delete(target_coordinates, non_core_symmetry_groups, axis=0)
        reference_coordinates_core = np.delete(reference_coordinates, non_core_symmetry_groups, axis=0)
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

        
            reference_group_idx = [idx for idx, _ in enumerate(molecule.get_labels()) if idx not in non_core_symmetry_groups]
            reference_group = np.delete(np.array(reference_coordinates), reference_group_idx, axis=0)
            current_group = np.delete(np.array(rotated_coordinates_core), reference_group_idx, axis=0)
        
            
            row_ind, col_ind, cost_matrix = self.assign_atoms(current_group, reference_group)
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

        return ref_structure_check, distance_core
    
    # This function is used to calculate the reference energy for different possible QM-Methods (from which xTB was used for the paper) where the veloxchem functionalities are used to caclualte the energy and the gradient
    # Orca and VeloxChem are more expensive methods which have been tried -- but due to time it was decided to use xTB as a first trial

    def qm_energy_calculator(self, mol, QM_Driver='XTB', basis_set_label='sto-3g', xc_func='b3lyp', gradient=False):

        qm_energies_kcal_mol = 0
        gradients = 0
        if QM_Driver == 'XTB':
            
            xtb_drv = XtbDriver()
            xtb_grad_drv = XtbGradientDriver(xtb_drv)

                
            xtb_drv.compute(mol)
            qm_energies_kcal_mol = xtb_drv.get_energy() * 627.509474
            if gradient is True:
                xtb_grad_drv.compute(mol)
                gradients = xtb_grad_drv.get_gradient()
        

        ### This section uses more elaborate and expensive methods and have not been used further due to time constained
        if QM_Driver == 'Veloxchem':

            vlx_driver = ScfRestrictedDriver() # -> QM - Driver implemented within Veloxchem to calculate the energy for a given vlx.Molecule based on the self-consistent-field formalism
            vlx_driver.xcfun = xc_func # -> exchange-correlation functionals stored within veloxchem 
            vlx_gradient_driver = ScfGradientDriver(vlx_driver) # -> QM - GradientDriver implemented within Veloxchem to calculate the gradient for a given vlx.Molecule based on the self-consistent-field formalism
            
            basis = basis = MolecularBasis.read(mol, basis_set_label,) # -> neccessary basis set for making the SCF calculations possible
                
            scf_result = vlx_driver.compute(mol, basis)
            qm_energies_kcal_mol = vlx_driver.scf_energy * 627.509474
                
            if gradient is True:

                vlx_gradient_driver.compute(scf_result)
                gradients = vlx_gradient_driver.gradient

        
        if QM_Driver == 'Orca':
            
            qm_driver = ExternalQMDriver('ORCA', 0, 0, 1, 2, hostname=False) # ->  ExternalDriver implemented within Veloxchem to communicate with another QM-program ORCA to calculate the energy for a given vlx.Molecule based on the self-consistent-field formalism
            qm_driver.spin_flip = False
            qm_driver.basis_set_label = basis_set_label
            qm_driver.xc_func = xc_func
            qm_grad_driver = ExternalQMGradientDriver(qm_driver) # ->  ExternalGradientDriver implemented within Veloxchem to communicate with another QM-program ORCA to calculate the gradient for a given vlx.Molecule based on the self-consistent-field formalism
            
            energy = qm_driver.compute(mol)
            qm_energies_kcal_mol = energy[0] * 627.509474
            if gradient is True:
                qm_grad_driver.compute(mol)
                gradient = qm_grad_driver.extract_gradients()
                gradients = gradient[0]
        
        return qm_energies_kcal_mol, gradients
    

    def add_point(self, molecule, labels, filename, datapoint_molecules, impes_dict, z_matrix, qm_driver, grad_driver, hess_driver, non_core_symmetry_group, basis=None, scf_results=None, mapping=True):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param label:
                the label for the new point to be added.
            :param basis:
                the basis set (if required)
        """
        if qm_driver is None:
            raise ValueError("No energy driver defined.")
        if grad_driver is None:
            raise ValueError("No gradient driver defined.")
        if hess_driver is None:
            raise ValueError("No Hessian driver defined.")

        energy = None
        gradient = None
        hessian = None
        
        symmetry_adjusted_coordinates = molecule.get_coordinates_in_bohr()
        if len(datapoint_molecules) > 0 and mapping is True:
            symmetry_adjusted_coordinates = self.database_symmetry_group_swapper(molecule.get_coordinates_in_bohr(), datapoint_molecules, non_core_symmetry_group=non_core_symmetry_group, molecule=molecule)

        symmetry_adjusted_molecule = Molecule(molecule.get_labels(), symmetry_adjusted_coordinates, units='bohr')

        energy, scf_results = self.compute_energy(qm_driver, symmetry_adjusted_molecule, basis)
        gradient, NACs = self.compute_gradient(grad_driver, symmetry_adjusted_molecule, basis, scf_results)
        hessian = self.compute_hessian(hess_driver, symmetry_adjusted_molecule, basis, scf_results)

        change = False
        # if self.adiabatic_basis:
        #     NACs_deriv = NACs
        #     ADTMatrix = NACs
        #     energy, gradient, hessian = self.adiabatic_transformation(energy, gradient, hessian, NACs, NACs_deriv, ADTMatrix)

        for root in range(1):
            impes_coordinate = ImpesCoordinates(z_matrix)
            impes_coordinate.update_settings(impes_dict)
            impes_coordinate.cartesian_coordinates = symmetry_adjusted_molecule.get_coordinates_in_bohr()

            impes_coordinate.energy = energy[root]
            impes_coordinate.gradient = gradient[root]

            impes_coordinate.hessian = hessian[root]
            impes_coordinate.transform_gradient_and_hessian()
            # if NACs != None and root < self.roots - 1:
            #     impes_coordinate.NAC = NACs[root]
            #     impes_coordinate.transform_cartesian_NACs_to_internal_NACs()
            #     self.NAC.append(NACs)
            # else:
            #     if self.calc_NAC == True:
            #         self.calc_NAC = False
            #         change = True

            impes_coordinate.write_hdf5(filename, labels[root], False)
           
            print(f"Database expansion with {', '.join(labels)}")
            for l, e in zip(labels, energy):
                print(f"{l}: Energy = {e:.4f} hatree")

        # if change == True:
        #     self.calc_NAC = True
        #     change = False


    def compute_energy(self, qm_driver, molecule, basis=None):
        """ Computes the QM energy using self.qm_driver.

            :param molecule:
                The molecule.
            :param basis:
                The basis set.

            :returns the QM energy.
        """
        # Dtermine the type of energy driver, to be able to
        # call it correctly.

        qm_energy = None
        scf_tensors = None

        # XTB
        if isinstance(qm_driver, XtbDriver):

            qm_driver.compute(molecule)
            qm_energy = qm_driver.get_energy()
            qm_energy = np.array([qm_energy])
            print('qm_energy', qm_energy, qm_energy[0])

        elif isinstance(qm_driver, ExternalQMDriver):
            qm_energy = qm_driver.compute(molecule)
            print('qm_energy', qm_energy)

        # restricted SCF
        elif isinstance(qm_driver, ScfRestrictedDriver):
            qm_driver.ostream.mute()
            scf_tensors = qm_driver.compute(molecule, basis)
            qm_energy = qm_driver.scf_energy
            qm_energy = np.array([qm_energy])
            print('qm_energy', qm_energy)
            qm_driver.ostream.unmute()

        # TDDFT / TDHF
        elif ( isinstance(self.qm_driver, TDAExciDriver) or
             isinstance(self.qm_driver, LinearResponseEigenSolver)):
            if self.grad_driver.scf_drv is None:
                raise ValueError("No SCF driver defined.")
            self.grad_driver = scf_drv.mute()
            scf_tensors = self.grad_driver.scf_drv.compute(molecule, basis)
            self.qm_driver.ostream.mute()
            self.qm_driver._is_converged = False
            self.rsp_results = self.qm_driver.compute(molecule, basis,
                                                      scf_tensors)
            qm_energy = ( self.grad_driver.scf_drv.scf_energy
                   + self.rsp.results['eigenvalues'][self.excited_state_index] )

            self.scf_drv.unmute()
            self.qm_driver.ostream.unmute()

        if qm_energy is None:
            error_txt = "Could not compute the QM energy. "
            error_txt += "Please define a QM driver."
            raise ValueError(error_txt)

        return qm_energy, scf_tensors

    def compute_derivatives(self, molecule, NACs=False):
        qm_gradient = None
        qm_hessian = None
        qm_nacs = None
        if isinstance(self.grad_driver, ExternalQMGradientDriver) and isinstance(self.hess_driver, ExternalQMHessianDriver):
            self.grad_driver.compute_gradient(molecule)
            jobs_finished = self.hess_driver.compute_analytical_hessian(molecule)
            if jobs_finished == 0:
                qm_gradient = self.grad_driver.extract_gradients()
                qm_hessian = self.hess_driver.extract_hessians()
            else:
                print('jobs have not been finshed or there was and error!')
                exit()

        return qm_gradient, qm_nacs, qm_hessian

    def compute_gradient(self, grad_driver, molecule, basis=None, scf_results=None):
        """ Computes the QM gradient using self.grad_driver.

            :param molecule:
                The molecule.
            :param basis:
                The basis set.

            :returns the QM gradient.
        """

        qm_gradient = None

        if isinstance(grad_driver, XtbGradientDriver):
            grad_driver.ostream.mute()
            grad_driver.compute(molecule)
            grad_driver.ostream.unmute()
            qm_gradient = grad_driver.gradient
            qm_gradient = np.array([qm_gradient])

        elif isinstance(grad_driver, ExternalQMGradientDriver):
            print('here is the gradient object', grad_driver)
            grad_driver.compute(molecule)
            qm_gradient = grad_driver.extract_gradients()


        elif isinstance(grad_driver, ScfGradientDriver):
            grad_driver.ostream.mute()
            grad_driver.compute(molecule, basis, scf_results)
            qm_gradient = grad_driver.gradient
            qm_gradient = np.array([qm_gradient])
            self.grad_driver.ostream.unmute()

        elif isinstance(self.grad_driver, TddftGradientDriver):
            self.grad_driver.ostream.mute()
            self.grad_driver.compute(molecule, basis, self.qm_driver,
                                         self.rsp_results)
            qm_gradient = self.grad_driver.gradient[self.excited_state_index]
            self.grad_driver.ostream.unmute()

        if qm_gradient is None:
            error_txt = "Could not compute the QM gradient. "
            error_txt += "Please define a QM gradient driver."
            raise ValueError(error_txt)
        
        print('Gradient', qm_gradient)
        return qm_gradient, None


    # TODO: mute outside to save time?
    def compute_hessian(self, hess_driver, molecule, basis=None, scf_results=None):
        """ Computes the QM Hessian using self.hess_driver.

            :param molecule:
                The molecule.
            :param basis:
                The basis set.

            :returns the QM Hessian matrix.
        """

        qm_hessian = None

        if isinstance(hess_driver, XtbHessianDriver):
            hess_driver.ostream.mute()
            hess_driver.compute(molecule)
            qm_hessian = hess_driver.hessian
            hess_driver.ostream.unmute()
            qm_hessian = np.array([qm_hessian])

        elif isinstance(hess_driver, ExternalQMHessianDriver):
            if hess_driver.analytical:
                hess_driver.compute(molecule)
                qm_hessian = hess_driver.extract_hessians()
            else:
                qm_hessian = hess_driver.compute_numerical_hessian(molecule)

        elif isinstance(hess_driver, ScfHessianDriver):
            hess_driver.ostream.mute()
            hess_driver.compute(molecule, basis, scf_results)
            qm_hessian = hess_driver.hessian
            qm_hessian = np.array([qm_hessian])
            hess_driver.ostream.unmute()

        elif isinstance(self.hesian_driver, TddftHessianDriver):
            self.hess_driver.ostream.mute()
            self.hess_driver.compute(molecule, basis)
            qm_hessian = self.hess_driver.hessian
            self.hess_driver.ostream.unmute()

        if qm_hessian is None:
            error_txt = "Could not compute the QM Hessian. "
            error_txt += "Please define a QM Hessian driver."
            raise ValueError(error_txt)


        return qm_hessian
    

    def database_symmetry_group_swapper(self, molecule_coordinates, ref_qm_data_points, non_core_symmetry_group, molecule):

        current_coordinates = molecule_coordinates.copy()

        
        
        new_molecular_coordinates = self.find_best_mapping(current_coordinates, ref_qm_data_points, non_core_symmetry_group, molecule)



        return new_molecular_coordinates
    
    def find_best_mapping(self, current_coordinate, datapoint_structures, non_core_symmetry_groups, molecule, data_point_idx=None):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                ImpesCoordinates object
        """

        #print('QM_Datapoint internal coordinates \n\n')
        
        best_cost = np.inf
        best_cost_matrix = None
        best_row_ind = None
        best_col_ind = None 
        for i, datapoint_mol in enumerate(datapoint_structures):
            
            datapoint_coord = datapoint_mol.get_coordinates_in_bohr()
            # First, translate the cartesian coordinates to zero
            target_coordinates = self.calculate_translation_coordinates_analysis(datapoint_coord)
            reference_coordinates = self.calculate_translation_coordinates_analysis(current_coordinate)

            # remove the non_core atoms (currently just H)
            target_coordinates_core = np.delete(target_coordinates, non_core_symmetry_groups, axis=0)
            reference_coordinates_core = np.delete(reference_coordinates, non_core_symmetry_groups, axis=0)

            # Then, determine the rotation matrix which
            # aligns data_point (target_coordinates)
            # to self.impes_coordinate (reference_coordinates)     

            rotation_matrix_core = geometric.rotate.get_rot(target_coordinates_core,
                                                    reference_coordinates_core)
            
            # Rotate the data point
            rotated_coordinates_core = np.dot(rotation_matrix_core, target_coordinates.T).T

            reference_group_idx = [idx for idx, _ in enumerate(molecule.get_labels()) if idx not in non_core_symmetry_groups]
            reference_group = np.delete(np.array(reference_coordinates), reference_group_idx, axis=0)
            current_group = np.delete(np.array(rotated_coordinates_core), reference_group_idx, axis=0)

        
            
            row_ind, col_ind, cost_matrix = self.assign_atoms(current_group, reference_group)

            #define the overall cost of the matrix
            total_cost = cost_matrix[row_ind, col_ind].sum()

            if total_cost < best_cost:
                
                best_cost = total_cost
                best_cost_matrix = cost_matrix
                best_row_ind = row_ind
                best_col_ind = col_ind

        swapped_indices = {}
        swapped_elements_map = {}
        swapped_indices_list = []
        reference_indices_list = []
        indices_that_need_to_be_swapped = []
        element = 0
        current_element_list = []
        for i, j in zip(best_row_ind, best_col_ind):
            swapped_indices[non_core_symmetry_groups[i]] = non_core_symmetry_groups[j]
            swapped_indices_list.append(non_core_symmetry_groups[i])
            current_element_list.append(non_core_symmetry_groups[j])
            reference_indices_list.append(non_core_symmetry_groups[j])
            indices_that_need_to_be_swapped.append(non_core_symmetry_groups[j])
            if len(current_element_list) == 3:
                swapped_elements_map[element] = tuple(current_element_list)
                current_element_list = []
                element += 1

        y_assigned = current_coordinate[reference_indices_list]
        ref_structure_check = current_coordinate.copy()
        ref_structure_check[swapped_indices_list] = y_assigned

        return ref_structure_check
    
    def assign_atoms(self, x, y):
    
        cost_matrix = np.linalg.norm(x[:, np.newaxis, :] - y[np.newaxis, :, :], axis=2)
        
        # Apply the Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        
        return row_ind, col_ind, cost_matrix
    
    def structures_to_xyz_file(self, molecules_for_xyz, structure_filename, im_energies=None, qm_energies=None, im_energies_2=None, im_energies_no_mapping=None, im_energies_2_no_mapping=None):

        with open(structure_filename, 'w') as file:
            pass

        for i, dyn_mol in enumerate(molecules_for_xyz):

            current_xyz_string = dyn_mol.get_xyz_string()

            xyz_lines = current_xyz_string.splitlines()

            if len(xyz_lines) >= 2 and im_energies_2 is not None and im_energies is not None:

                xyz_lines[1] += f'Energies  QM: {qm_energies[i]}  IM: {im_energies[i]}  IM_2: {im_energies_2[i]}  IM_1_NM: {im_energies_no_mapping[i]}  IM_2_NM: {im_energies_2_no_mapping[i]}   Diff_1: {abs(qm_energies[i] - im_energies[i])}'


            updated_xyz_string = "\n".join(xyz_lines)

            with open(structure_filename, 'a') as file:
                file.write(f"{updated_xyz_string}\n\n")
    
    def read_in_xyz_traj(self, xyz_traj):

        molecules_along_traj = []

        xyz_strings = []
        qm_energies = []
        im_energies_1 = []
        im_energies_2 = []

        im_energies_1_nm = []
        im_energies_2_nm = []

        with open(xyz_traj, 'r') as file:
            lines = file.readlines()

        current_xyz = []
        for line in lines: # read in lines from the file
            line = line.strip()
            if line.isdigit():  
                if current_xyz: 
                    xyz_strings.append('\n'.join(current_xyz))
                    current_xyz = []  
                current_xyz.append(line) 
            elif "Energies" in line: # Extract the information of the second line that holds energy information

                match = re.search(r"QM:\s*([\-\d\.]+)\s+IM:\s*([\-\d\.]+)\s+IM_2:\s*([\-\d\.]+)\s+IM_1_NM:\s*([\-\d\.]+)\s+IM_2_NM:\s*([\-\d\.]+)\s+Diff_1:\s*([\-\d\.]+)", line) # group the individual energies to append them in a list
                if match:
                    qm_energies.append(float(match.group(1)))
                    im_energies_1.append(float(match.group(2)))
                    im_energies_2.append(float(match.group(3)))
                    im_energies_1_nm.append(float(match.group(4)))
                    im_energies_2_nm.append(float(match.group(4)))

                current_xyz.append(line) # still need to be stored as string because vlx reads it in as xyz string
            else:
                current_xyz.append(line)
        


        if current_xyz:
            xyz_strings.append('\n'.join(current_xyz))

        for xyz_string in xyz_strings:
            molecules_along_traj.append(Molecule.from_xyz_string(xyz_string))

        molecules_along_traj[0].show()
        return molecules_along_traj, qm_energies, im_energies_1, im_energies_2, im_energies_1_nm, im_energies_2_nm
    
    def create_database_from_trajextory(self, structures_as_xyz, interpolation_settings, ff_data, qm_driver, grad_driver, hess_driver, database_name='interpolation_database.h5', mapping=True):
        

        molecule_traj, _, _, _, _, _ = self.read_in_xyz_traj(structures_as_xyz)

        forcefield_generator = ForceFieldGenerator()
        forcefield_generator.force_field_data = ff_data
        forcefield_generator.partial_charges = molecule_traj[0].get_partial_charges(molecule_traj[0].get_charge())
        forcefield_generator.create_topology(molecule_traj[0])
        _, non_core_symmetric_groups_combined = self.symmetric_groups(molecule_traj[0], forcefield_generator)

        for i, mol in enumerate(molecule_traj):
            
            label = [f'point_{i}']

            self.add_point(mol, 
                           label, 
                           database_name, 
                           molecule_traj[:i], 
                           interpolation_settings,
                           molecule_traj[0].get_z_matrix(), 
                           qm_driver, 
                           grad_driver, 
                           hess_driver, 
                           non_core_symmetric_groups_combined, mapping=mapping)
            
    
    def create_new_database_from_old_db(self, old_db, ref_mol, interpolation_settings, ff_data, qm_driver, grad_driver, hess_driver, database_name='interpolation_database.h5'):
        
        molecule_traj, _, _, _ = self.database_extracter(old_db, ref_mol)

        forcefield_generator = ForceFieldGenerator()
        forcefield_generator.force_field_data = ff_data
        forcefield_generator.partial_charges = molecule_traj[0].get_partial_charges(molecule_traj[0].get_charge())
        forcefield_generator.create_topology(molecule_traj[0])
        _, non_core_symmetric_groups_combined = self.symmetric_groups(molecule_traj[0], forcefield_generator)

        for i, mol in enumerate(molecule_traj):
            
            label = [f'point_{i}']

            self.add_point(mol, 
                           label, 
                           database_name, 
                           molecule_traj[:i], 
                           interpolation_settings,
                           molecule_traj[0].get_z_matrix(), 
                           qm_driver, 
                           grad_driver, 
                           hess_driver, 
                           non_core_symmetric_groups_combined)