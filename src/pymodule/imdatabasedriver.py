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
from contextlib import redirect_stderr
from io import StringIO

from .impesdriver import ImpesDriver
from .impescoordinates import ImpesCoordinates
from .imdatabasepointcollecter import IMDatabasePointCollecter
from .mmforcefieldgenerator import MMForceFieldGenerator
from .molecule import Molecule

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class ImDatabaseDriver:
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
        
        self.density_of_datapoints = None
        self.qm_data_points = None
        self.qmlabels = None
        self.molecules_along_rp = None
        self.impes_dict = None
        self.molecule = None
        self.datafile = None
        self.dihedrals = None
        self.allowed_devation = None
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

    
    def define_z_matrix(self, molecule):
        """
        Creates the z-matrix of redundant internal coordinates based on the
        topology from geomeTRIC.

        :return:
            a list of 2-tuples, 3-tuples, and 4-tuples corresponding to all bonds,
            bond agles, and respectively dihedral angles in the molecule.
        """

        g_molecule = geometric.molecule.Molecule()
        g_molecule.elem = molecule.get_labels()
        g_molecule.xyzs = [molecule.get_coordinates_in_bohr() * geometric.nifty.bohr2ang]

        g_molecule.build_topology()
        g_molecule.build_bonds()

        bonds = g_molecule.Data['bonds']
        angles = g_molecule.find_angles()
        dihedrals = g_molecule.find_dihedrals()

        z_matrix = []
        for bond in bonds:
            z_matrix.append(bond)
        for angle in angles:
            z_matrix.append(angle)
        for dihedral in dihedrals:
            z_matrix.append(dihedral)

        return z_matrix

    def set_up_the_system(self, molecule, impes_dict, dynamics_settings, target_dihedrals=None, sampling_structures=1):

        """
        Assign the neccessary variables with respected values. 

        :param molecule: original molecule

        :param impes_dict: user defined dictionary for the interpolation cacluations (see IMPES driver)

        :param dynamics_settings: sets all important parameters for the dynamics:

        :param target_dihedrals: is a list of dihedrals that should be scanned during the dynamics
            
        """


        self.qm_data_points = None
        self.molecule = molecule
        self.impes_dict = impes_dict
        if dynamics_settings:
            self.dynamics_settings = dynamics_settings
        if self.z_matrix is None:
            self.z_matrix = self.define_z_matrix(molecule)
        
        # Read the database to determine the  
        
        if 'checkpoint_file_name' in self.impes_dict:
            data_base_file = self.impes_dict['checkpoint_file_name']
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
        
        self.density_of_datapoints, self.molecules_along_rp, self.allowed_devation = self.database_density_check_with_molecule(molecule, self.qm_data_points, specific_dihedrals=target_dihedrals, nsampling=sampling_structures)
        self.dihedrals = target_dihedrals
        self.impes_dict['dihedrals'] = target_dihedrals
        
        print('here is the density', self.density_of_datapoints)

    def run_data_base_construction(self):
        
        for counter, entry in enumerate(self.molecules_along_rp.items()):
            key, molecules = entry
            for mol in molecules:
                
                if self.dihedrals is not None:
                    key = (key, int(round(mol.get_dihedral_in_degrees(key))))
                    print('key', key, self.density_of_datapoints)
                #structure = structure * bohr_in_angstrom()
                #current_molecule = Molecule(self.molecule.get_labels(), structure, units="angstrom")
                forcefield_generator = MMForceFieldGenerator()
                self.dynamics_settings['trajectory_file'] = f'trajectory_{counter}.pdb'
                forcefield_generator.partial_charges = mol.get_partial_charges(mol.get_charge())
                
                forcefield_generator.create_topology(mol)
                
                # different_bins = False
                # if counter > 0:

                #     if different_bins:
                #         self.impes_dict['basename'] = f'im_collcection_run{counter}'
                #     else:
                #         if 'checkpoint_file_name' not in self.impes_dict:
                #             self.impes_dict['checkpoint_file_name'] = f'{self.impes_dict['basename']}.h5'
                        
                #         self.datafile = self.impes_dict['checkpoint_file_name']
                #         impes_driver = ImpesDriver(self.z_matrix)
                #         impes_driver.update_settings(self.impes_dict)
                #         self.qmlabels = impes_driver.read_labels()

                #         self.qm_data_points = []
                #         self.qm_energies = []
                #         for label in self.qmlabels:
                #             qm_data_point = ImpesCoordinates(self.z_matrix)
                #             qm_data_point.read_hdf5(self.datafile, label)
                #             self.qm_energies.append(qm_data_point.energy)
                #             self.qm_data_points.append(qm_data_point)
                #         new_density_around_data_points, _ =  self.database_density_check_with_molecule(mol, self.qm_data_points)
                #         print('here is the denisty after first run', new_density_around_data_points)
                #         self.density_of_datapoints[counter] = new_density_around_data_points
                    
                im_database_driver = IMDatabasePointCollecter()
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
                    im_database_driver.add_bias_force([2, 3, 4, 8], 12, 180)
                nroots = int(self.dynamics_settings['roots']) + 1
                desiered_point_density = int(self.dynamics_settings['desired_datapoint_density'])

                if self.density_of_datapoints[key] / nroots < desiered_point_density: 
                    #density_around_data_points = im_database_driver.run_qmmm(self.density_of_datapoints[counter], desiered_point_density, self.distance_threshold, unadded_cycles, nroots, calc_NAC, checkpoint, ensemble, nsteps=nsteps, temperature=temperature, timestep=timestep, snapshots=snapshots, out_file=f'{counter}_run_{out}')            
                    im_database_driver.density_around_data_point = [self.density_of_datapoints[key], key[0]]
                    
                    print('Allowed deviation', self.allowed_devation)
                    
                    im_database_driver.allowed_molecule_deviation = [key[1] - self.allowed_devation, -key[1] + self.allowed_devation]
                    im_database_driver.update_settings(self.dynamics_settings, self.impes_dict)
                    im_database_driver.run_qmmm()
                    self.density_of_datapoints[key] = im_database_driver.density_around_data_point

                    # individual impes run objects
                    self.qm_energies.append(im_database_driver.qm_potentials)
                    self.total_energies.append(im_database_driver.total_energies)
                    self.kinetic_energies.append(im_database_driver.kinetic_energies)
                    self.molecules.append(im_database_driver.molecules)
                    self.point_added_molecules.append(im_database_driver.point_adding_molecule)
                    self.unique_molecules.append(im_database_driver.allowed_molecules)


                counter += 1
        
        print('The construction of the database was sucessfull')    

    
    def database_density_check_with_molecule(self, molecule, qm_datapoints, specific_dihedrals=None, nsampling=None):
        """
        loads the molecule and specific features (dihedral) that is changed in order to sample different molecular conformations.

        :param filename: 
            Path to the PDB file. Example: 'system.pdb'.
        """
        sampled_molecules = {}
        point_densities = {}
        rotation_values = np.linspace(0, 360, nsampling, endpoint=False)
        if specific_dihedrals is not None:
            for specific_dihedral in specific_dihedrals:
                sampled_molecules[tuple(specific_dihedral)] = []
                for theta in rotation_values:
                    molecule.set_dihedral_in_degrees(specific_dihedral, theta)
                    print('get_the molecule dihedral', molecule.get_dihedral_in_degrees(specific_dihedral))
                    new_molecule = Molecule(molecule.get_labels(), molecule.get_coordinates_in_bohr(), 'bohr')
                    sampled_molecules[tuple(specific_dihedral)].append(new_molecule)

                    key = (tuple(specific_dihedral), int(molecule.get_dihedral_in_degrees(specific_dihedral)))
                    point_densities[key] = 0
        
        else:
            sampled_molecules[None, 180] = [molecule]
            point_densities[None, 180] = 0
            
            if qm_datapoints is not None:
                point_densities[None, 180] = len(qm_datapoints)
            
    
        def dihedral_to_vector(angle):
            """
            Converts a dihedral angle in degrees to a 2D vector on the unit circle.
            """
            rad = np.radians(angle)
            return np.array([np.cos(rad), np.sin(rad)])

        def structure_to_vector(dihedrals):
            """
            Converts a list of dihedral angles to a concatenated vector.
            For N dihedrals, returns a vector of length 2N.
            """
            return np.concatenate([dihedral_to_vector(angle) for angle in dihedrals])

        def dihedral_distance_vectorized(dihedrals1, dihedrals2):
            """
            Computes the Euclidean distance between two sets of dihedrals by mapping each
            angle to a unit circle vector.
            
            Parameters:
                dihedrals1, dihedrals2: Lists or tuples of angles (in degrees).
            
            Returns:
                A scalar distance.
            """
            vec1 = structure_to_vector(dihedrals1)
            vec2 = structure_to_vector(dihedrals2)
            return np.linalg.norm(vec1 - vec2)

        # Example usage:

        if qm_datapoints and specific_dihedrals is not None:
            impes_driver = ImpesDriver(self.z_matrix)
            impes_driver.update_settings(self.impes_dict)
            for specific_dihedral in specific_dihedrals:
                for point in qm_datapoints:
                    min_distance = np.inf
                    key = None
                    for i, mol in enumerate(sampled_molecules[tuple(specific_dihedral)]):
                        datapoint_molecule = Molecule(self.molecule.get_labels(), point.cartesian_coordinates, 'bohr')
                        dihedrals_of_mol = [mol.get_dihedral_in_degrees(specific_dihedral)]
                        dihedrals_of_dp = [datapoint_molecule.get_dihedral_in_degrees(specific_dihedral)]
                        distance_vectorized = dihedral_distance_vectorized(dihedrals_of_mol, dihedrals_of_dp)

                        if abs(distance_vectorized) < min_distance:
                            min_distance = abs(distance_vectorized)
                            key = (tuple(specific_dihedral), int(mol.get_dihedral_in_degrees(specific_dihedral)))
                    
                    point_densities[key] += 1

        print(point_densities)
        
        allowed_devation = 360 / nsampling
        normalized_angle = ((allowed_devation + 180) % 360) - 180

        print('normalized angle', normalized_angle)

        return point_densities, sampled_molecules, normalized_angle
    

 
