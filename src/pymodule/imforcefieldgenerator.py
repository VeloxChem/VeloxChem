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
import os
import random
from contextlib import redirect_stderr
from io import StringIO

from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .molecularbasis import MolecularBasis
from .openmmdynamics import OpenMMDynamics
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .imdatabasepointcollecter import IMDatabasePointCollecter
from .mmforcefieldgenerator import MMForceFieldGenerator
from .molecule import Molecule
from .errorhandler import assert_msg_critical
from. veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class IMForceFieldGenerator:
    """
    Class to set up and control the construction of the Interpolation Dynamics (IM) database.
    
    This class handles multiple aspects of database creation for force field generation, including
    molecule sampling, quantum mechanical calculations, and interpolation settings. It manages 
    data such as molecular structures, energy calculations, and control parameters for dynamics.

    Instance Variables:

        - density_of_datapoints: Tracks the number of current interpolation data points in the database. Used to monitor 
                                 database size and data density.

        - qm_data_points: Stores interpolation data points.

        - qmlabels: Labels assigned to each interpolation data point in the database.

        - molecules_along_rp: Represents molecular structures along a predefined reaction path or internal coordinate 
                              pathway.

        - dihedrals: A list of dihedral angles to be rotated or scanned during simulations. Used to determine a pre-
                        defined path that should be sampled within the interpolation database construction.

        - sampling_structures: Specifies how many structures to generate around rotatable dihedral angles for database 
                               population.

        - molecule: The initial molecular structure. (This mus be provided by the user)

        - datafile: Represents the database file (interpolation forcefield) used to store molecular structures and data points. 
                    Typically initialized as `IMDatabase.h5`.

        - z_matrix: The original Z-matrix (internal coordinate definition) for the molecule, specifying bonds, angles, 
                    and dihedrals. It serves as the basis for internal coordinate transformations.

        - allowed_deviation: A threshold for how much a generated or sampled 
                            structure is allowed to deviate from an expected configuration. Ensures structural integrity 
                            during sampling or dynamics.

        - angle_threshold: Defines the range within which dihedral angles can vary during sampling and dynamics.
                           Making sure the smapling for 1 structure stays within a certain constained space.

        - interpolation_settings: A dictionary containing settings for the interpolation.

        - interpolation_type: Defines the type of interpolation used in the force field generation. The default value is 
                              'shepard', likely referring to Shepard interpolation, a type of distance-weighted interpolation.

        - qm_driver, qm_grad_driver, qm_hess_driver: Instances of drivers for quantum mechanical (QM) calculations, including 
                                                     single-point energy (qm_driver), gradient (qm_grad_driver), and Hessian 
                                                     (qm_hess_driver) calculations. These drivers are used to perform QM tasks 
                                                     during database construction. (Currently given from the user)

        - dynamics_settings: Holds various settings for molecular dynamics simulations, such as temperature, pressure, 
                             force constants, and timestep.

        - basis_set_label: Specifies the basis set for QM calculations, with a default value of 'def2-svp'.

        - xcfun: Specifies the exchange-correlation functional for QM calculations. Default is 'b3lyp'.

        - duration: Specifies the total steps that can pass without a point being added to the database to determine
                    early breaks for a current structure.

        - temperature: Temperature for molecular dynamics simulations, defaulting to 150.15 K.

        - pressure: Pressure value for molecular dynamics, defaulting to 1.0 atm.

        - force_constant: The force constant used in simulations or sampling dynamics, with a default value of 1.0.

        - ensemble: Specifies the statistical ensemble used in the dynamics. Default: 'NVE'.

        - timestep: The time increment (in femtoseconds) between simulation steps. Default: 0.5 fs.

        - nsteps: Number of steps to be performed in the dynamics simulation.

        - snapshots: Specifies how many snapshots of the trajectory to record. Defaults: nsteps.

        - trajectory_file: Name of the file to store the trajectory of the molecular simulation. Default: 'trajectory.pdb'.

        - desired_point_density: Defines the desired density of data points for 1 structure in the database. Default: 50.

        - converged_cycle: Defines the number of cycles required for a simulation or database sampling to be considered 
                           converged. Default: 4.

        - energy_threshold: Specifies an energy threshold to determine when a structure is necessary to be added into 
                            the interpolation database. Default: 1.5 kcal/mol.

        - start_collect: Specifies at which step in the simulation interpolation datapoint collection should begin. Default: 0.

        - solvent: Specifies the solvent environment for the dynamics. Default: 'gas'. (Should always be gas for the construction)

        - FF_datafile: Specifies the file path for a force field data file, defaulting to a specific location.

        - qm_energies: A list to store QM energy from individual simulations or calculations (kj/mol).

        - total_energies: A list to store total energy (kj/mol).

        - molecules: A list to store molecular structures sampled during simulations (kj/mol).

        - kinetic_energies: A list to store kinetic energy from simulations (kj/mol).

        - point_added_molecules: A list of molecules for which new data points were added to the database.

        - unique_molecules: A list of unique molecular structures identified during database construction.

        - dynamics_method: Determines the method to generate molecular structures for the database quality conformation.
        
        - nstruc_to_confirm_database_quality: Number of randomly selected strucutures for the database quality check.
    """
    
    def __init__(self):

        self.density_of_datapoints = None
        self.qm_data_points = None
        self.qmlabels = None
        self.molecules_along_rp = None
        self.dihedrals = None
        self.sampling_structures = 1
        

        self.molecule = None
        self.datafile = None
        self.dihedrals = None
        self.allowed_deviation = None
        self.z_matrix = None

        # should be necessary to initialize
        self.qm_driver = ScfRestrictedDriver()
        self.qm_grad_driver = ScfGradientDriver(self.qm_driver)
        self.qm_hess_driver = ScfHessianDriver(self.qm_driver)

        # variables for the interpolation
        self.interpolation_settings = None
        self.interpolation_type = 'shepard'
        self.exponent_p = '12'
        self.exponent_q = '2'
        self.confidence_radius = '0.5'
        self.imforcefieldfile = None
        file_exists = 'IMDatabase.h5' in os.listdir(os.getcwd())
        if file_exists:
            print('IMPORTANT: IM ForceFieldFile is initalized from the current directory as IMDatabase.h5')
            self.imforcefieldfile = 'IMDatabase.h5'
        # variables for the forcefield generation and database expansion
        self.dynamics_settings = None
        self.basis_set_label = 'def2-svp'
        self.xcfun = 'b3lyp'
        self.duration = 2000
        self.temperature = 150.15
        self.pressure = 1.0
        self.force_constant = 1.0
        self.ensemble = 'NVE'
        self.timestep = 0.5
        self.friction = 1.0
        self.nsteps = 1000
        self.snapshots = self.nsteps
        self.trajectory_file = 'trajectory.pdb'
        self.desired_point_density = 50
        self.converged_cycle = 5
        self.energy_threshold = 1.5
        self.start_collect = 0
        self.solvent = 'gas'

        self.FF_datafile = '/home/vlind06/phd_work/interpolation/gaff-2.11.dat'

        # individual run
        self.qm_energies = []
        self.total_energies = []
        self.molecules = None
        self.kinetic_energies = []
        self.point_added_molecules = []
        self.unique_molecules = []


        # confirm database quality
        self.dynamics_method = 'MM'
        self.nstruc_to_confirm_database_quality = 50

        self.interpolation_settings = { 'interpolation_type':self.interpolation_type, 
                            'exponent_p':self.exponent_p,
                            'exponent_q':self.exponent_q, 
                            'confidence_radius':self.confidence_radius,
                            'imforcefield_file':self.imforcefieldfile,
                            'use_inverse_bond_length':True
                          }

        self.dynamics_settings = {  'qm_driver': self.qm_driver,
                                    'grad_driver': self.qm_grad_driver,
                                    'hess_driver': self.qm_hess_driver,
                                    'basis_set':self.basis_set_label,
                                    'xc_fun':self.xcfun,
                                    'duration':self.duration, 'temperature':self.temperature, 'solvent':self.solvent,
                                    'pressure':self.force_constant, 'force_constant': self.force_constant, 'ensemble':self.ensemble,
                                    'timestep': self.timestep, 'nsteps': self.nsteps, 'friction':self.friction,
                                    'snapshots':self.snapshots, 'trajectory_file':self.trajectory_file,
                                    'desired_datapoint_density':self.desired_point_density, 'converged_cycle': self.converged_cycle, 'energy_threshold':self.energy_threshold,
                                    'NAC':False, 'load_system': None, 'collect_qm_points_from':self.start_collect,
                                    'FF_datafile':self.FF_datafile}


    def set_up_the_system(self, molecule, target_dihedrals=None, sampling_structures=1):

        """
        Assign the neccessary variables with respected values. 

        :param molecule: original molecule

        :param target_dihedrals: is a list of dihedrals that should be scanned during the dynamics

        :param sampling_structures: devides the searchspace around given rotatbale dihedrals
            
        """


        self.qm_data_points = None
        self.molecule = molecule
        # self.z_matrix = sorted(self.z_matrix, key=len)

        if self.z_matrix is None and self.imforcefieldfile is None:
            self.z_matrix = self.define_z_matrix(molecule)
        # Read the database to determine the  
        
        if self.imforcefieldfile is not None:

            impes_driver = InterpolationDriver(self.z_matrix)
            impes_driver.update_settings(self.interpolation_settings)
            self.qmlabels, self.z_matrix = impes_driver.read_labels()

            self.qm_data_points = []
            self.qm_energies = []
            for label in self.qmlabels:
                qm_data_point = InterpolationDatapoint(self.z_matrix)
                qm_data_point.read_hdf5(self.imforcefieldfile, label)
                self.qm_energies.append(qm_data_point.energy)
                self.qm_data_points.append(qm_data_point)        
        
        self.density_of_datapoints, self.molecules_along_rp, self.allowed_deviation = self.database_density_check_with_molecule(molecule, self.qm_data_points, specific_dihedrals=target_dihedrals, nsampling=sampling_structures)

        self.interpolation_settings['dihedrals'] = target_dihedrals

    def compute(self, molecule):

        """
        Construct the interpolation dynamics database by generating molecular structures, 
        performing QM calculations, and collecting data points.

        :param molecule: The input molecular structure for which the database is constructed.
                        This molecule serves as a reference for sampling and simulation tasks.

        The method sets up quantum mechanical drivers, initializes interpolation and dynamics settings,
        and iterates over molecular structures along a predefined reaction path. It runs simulations
        to expand/generate the interpolation forcefield with new data points.

        """
        
        # First set up the system for which the database needs to be constructed
        if self.qm_driver is ScfRestrictedDriver:
            assert_msg_critical(self.qm_driver, 'IMForceFieldGenerator: QM-Driver/ QM-Gradient-Driver / QM-Hessian-Driver were initialized as ScF Method!.')
        
        print('MOlecule labels', molecule.get_labels())
        self.interpolation_settings = { 'interpolation_type':self.interpolation_type, 
                            'exponent_p':self.exponent_p,
                            'exponent_q':self.exponent_q, 
                            'confidence_radius':self.confidence_radius,
                            'imforcefield_file':self.imforcefieldfile,
                            'use_inverse_bond_length':True
                          }

        self.dynamics_settings = {  'qm_driver': self.qm_driver,
                                    'grad_driver': self.qm_grad_driver,
                                    'hess_driver': self.qm_hess_driver,
                                    'basis_set':self.basis_set_label,
                                    'xc_fun':self.xcfun,
                                    'duration':self.duration, 'temperature':self.temperature, 'solvent':self.solvent,
                                    'pressure':self.force_constant, 'force_constant': self.force_constant, 'ensemble':self.ensemble,
                                    'timestep': self.timestep, 'nsteps': self.nsteps, 'friction':self.friction,
                                    'snapshots':self.snapshots, 'trajectory_file':self.trajectory_file,
                                    'desired_datapoint_density':self.desired_point_density, 'converged_cycle': self.converged_cycle, 'energy_threshold':self.energy_threshold,
                                    'NAC':False, 'load_system': None, 'collect_qm_points_from':self.start_collect,
                                    'FF_datafile':self.FF_datafile}
        

        self.set_up_the_system(molecule, self.dihedrals, self.sampling_structures)

        for counter, entry in enumerate(self.molecules_along_rp.items()):
            key_old, molecules = entry
            for mol in molecules:
                
                key = key_old
                if self.dihedrals is not None:
                    
                    key = (key_old, int(round(mol.get_dihedral_in_degrees(key_old))))
                

                forcefield_generator = MMForceFieldGenerator()
                forcefield_generator.force_field_data = self.dynamics_settings['FF_datafile']
                self.dynamics_settings['trajectory_file'] = f'trajectory_{counter}.pdb'
                forcefield_generator.partial_charges = mol.get_partial_charges(mol.get_charge())
                
                forcefield_generator.create_topology(mol)
                    
                im_database_driver = IMDatabasePointCollecter()
                im_database_driver.platform = 'CUDA'

                im_database_driver.system_from_molecule(mol, self.z_matrix, forcefield_generator, solvent=self.solvent, qm_atoms='all')  
                desiered_point_density = int(self.dynamics_settings['desired_datapoint_density'])

                if self.density_of_datapoints[key] < desiered_point_density: 
                    im_database_driver.density_around_data_point = [self.density_of_datapoints[key], key[0]]
                    
                    if key[0] is None:
                        im_database_driver.allowed_molecule_deviation = [key[1] - self.allowed_deviation, -key[1] + self.allowed_deviation]
                    else:
                        im_database_driver.allowed_molecule_deviation = [((key[1] - self.allowed_deviation + 180) % 360) - 180, ((key[1] + self.allowed_deviation + 180) % 360) - 180]
                    
                    im_database_driver.update_settings(self.dynamics_settings, self.interpolation_settings)
                    if self.interpolation_settings['imforcefield_file'] is None:
                        self.interpolation_settings['imforcefield_file'] = im_database_driver.qm_datafile

                    im_database_driver.run_qmmm()
                    self.density_of_datapoints[key] = im_database_driver.density_around_data_point

                    # individual impes run objects
                    self.qm_energies.append(im_database_driver.qm_potentials)
                    self.total_energies.append(im_database_driver.total_energies)
                    self.kinetic_energies.append(im_database_driver.kinetic_energies)
                    self.molecules = im_database_driver.molecules
                    self.point_added_molecules.append(im_database_driver.point_adding_molecule)
                    self.unique_molecules.append(im_database_driver.allowed_molecules)

                
                entries = list(self.molecules_along_rp.values())

                if self.imforcefieldfile is None:
                    self.imforcefieldfile = 'IMDatabase.h5'

                self.confirm_database_quality(entries[0][0], self.imforcefieldfile, given_molecular_strucutres=self.molecules)
    
                counter += 1
        
        print('The construction of the database was sucessfull')

    
    def database_density_check_with_molecule(self, molecule, qm_datapoints, specific_dihedrals=None, nsampling=None):
        """

        Sample molecular structures by rotating specific dihedrals if defined and determine the current density of 
        datapoints for given structure (or structures).

        :param molecule: The original molecule object on which rotations are performed.

        :param qm_datapoints: A list of interpolation data points used to track how many
                            structures already exist (if specific_dihedrals: for certain dihedral configurations).

        :param specific_dihedrals: A list of dihedral angle definitions (as tuples of atoms) that 
                                will be scanned. If not provided, no specific dihedrals are rotated.

        :param nsampling: The number of samples to generate by rotating each dihedral from 0 to 360 degrees.
                        The rotation values are evenly spaced based on this parameter.

        The method creates sampled molecular structures by setting each specified dihedral angle to different 
        rotation values. The rotated structures are stored in `sampled_molecules`. Additionally, it initializes
        or updates `point_densities`, which tracks how many data points exist for each dihedral configuration.
        If no specific dihedrals are provided, the method uses a default structure at a 180-degree dihedral.

        :returns:

        - sampled_molecules: A dictionary where each key is a specific dihedral (or `None` if no dihedral is given),
                            and the value is a list of sampled molecular structures.
        
        - point_densities: A dictionary where keys are tuples of (dihedral, angle) and values represent the 
                            number of existing quantum mechanical data points for that configuration.

        - normalized_angle: determines the normalized dihedral angle how much the the angle is allowed to change within
                            the dynamics
        """

        sampled_molecules = {}
        point_densities = {}
        rotation_values = np.linspace(0, 360, nsampling, endpoint=False)
        if specific_dihedrals is not None:
            for specific_dihedral in specific_dihedrals:
                sampled_molecules[tuple(specific_dihedral)] = []
                for theta in rotation_values:
                    molecule.set_dihedral_in_degrees(specific_dihedral, theta)
                    print('get_the molecule dihedral', int(round(molecule.get_dihedral_in_degrees(specific_dihedral))))
                    new_molecule = Molecule(molecule.get_labels(), molecule.get_coordinates_in_bohr(), 'bohr')
                    sampled_molecules[tuple(specific_dihedral)].append(new_molecule)

                    key = (tuple(specific_dihedral), int(round(molecule.get_dihedral_in_degrees(specific_dihedral))))
                    point_densities[key] = 0
        
        else:
            sampled_molecules[None, 180] = [molecule]
            point_densities[None, 180] = 0
            
            if qm_datapoints is not None:
                point_densities[None, 180] = len(qm_datapoints)
            
    
        def dihedral_to_vector(angle):
            """
            Converts a dihedral angle in degrees to a 2D vector on the unit circle.

            :param angle:
                Angle to be transformed into a sinus and cosinus basis.
            
            :returns:
                Angle in vector form.
            """
            rad = np.radians(angle)
            return np.array([np.cos(rad), np.sin(rad)])

        def structure_to_vector(dihedrals):
            """
            Converts a list of dihedral angles to a concatenated vector.
            For N dihedrals, returns a vector of length 2N.

            :param dihedrals:
                List of dihedrals which will be transformed into vector form.
            
            :return:
                A concetenate list of the vectors.

            """
            return np.concatenate([dihedral_to_vector(angle) for angle in dihedrals])

        def dihedral_distance_vectorized(dihedrals1, dihedrals2):
            """
            Computes the Euclidean distance between two sets of dihedrals by mapping each
            angle to a unit circle vector.
            
            :param dihedrals1:
                Lists or tuples of angles (in degrees).
            
            :param dihedrals2:
                Lists or tuples of angles (in degrees).
            
            :returns:
                Norm of the distance.
            """
            vec1 = structure_to_vector(dihedrals1)
            vec2 = structure_to_vector(dihedrals2)
            return np.linalg.norm(vec1 - vec2)

        if qm_datapoints and specific_dihedrals is not None:
            impes_driver = InterpolationDriver(self.z_matrix)
            impes_driver.update_settings(self.interpolation_settings)
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
                            key = (tuple(specific_dihedral), int(round(mol.get_dihedral_in_degrees(specific_dihedral))))
                    
                    point_densities[key] += 1


        allowed_deviation = 360 / nsampling
        normalized_angle = ((allowed_deviation + 180) % 360) - 180

        return point_densities, sampled_molecules, normalized_angle
    

    def calculate_translation_coordinates_analysis(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0).
        
           :param given_coordinates:
                Coordinate that is translated to the center.
        """
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates
    
    def calculate_distance_to_ref(self, current_coordinates, datapoint_coordinate):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param current_coordinates:
                current molecular coordinates.
           
           :param data_point:
                InterpolationDatapoint object.

           :returns:
              Norm of the distance between 2 structures.
        """

        # First, translate the cartesian coordinates to zero
        target_coordinates = self.calculate_translation_coordinates_analysis(datapoint_coordinate)
        reference_coordinates = self.calculate_translation_coordinates_analysis(current_coordinates)

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)     
        rotation_matrix_core = geometric.rotate.get_rot(target_coordinates,
                                                reference_coordinates)
        

        # Rotate the data point
        rotated_coordinates_core = np.dot(rotation_matrix_core, target_coordinates.T).T
        # Calculate the Cartesian distance
        ref_structure_check = reference_coordinates.copy()
        distance_core = (np.linalg.norm(rotated_coordinates_core - ref_structure_check))

        return distance_core
    
    def database_extracter(self, datafile, mol_labels):
        """Extracts molecular structures from a given database file.

        :param datafile:
            Database file containing interpolation data.
        
        :param mol_labels:
            List of molecular labels.

        :returns:
            A list of VeloxChem Molecule objects extracted from the database.
        
        """
        
        im_driver = InterpolationDriver() # -> implemented Class in VeloxChem that is capable to perform interpolation calculations for a given molecule and provided z_matrix and database
        im_driver.imforcefield_file = datafile
        labels, z_matrix = im_driver.read_labels()
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        impes_coordinate = InterpolationDatapoint(z_matrix) # -> implemented Class in VeloxChem that handles all transformations and database changes concerning the interpolation
        data_point_molecules = []

        for label in sorted_labels:
            impes_coordinate.read_hdf5(datafile, label) # -> read in function from the ImpesDriver object
            coordinates_in_angstrom = impes_coordinate.cartesian_coordinates * bohr_in_angstrom()
            current_molecule = Molecule(mol_labels, coordinates_in_angstrom, 'angstrom') # -> creates a VeloxChem Molecule object
            
            data_point_molecules.append(current_molecule)

        return data_point_molecules
    
    def simple_run_dynamics(self, molecule, forcefield_generator):
        """Performs a molecular dynamics simulation using the selected method.
        This function runs a molecular dynamics simulation using either:
        - OpenMM for molecular mechanics (MM) sampling.
        - Interpolation-based QM/MM sampling.

        :param molecule:
            A VeloxChem Molecule object representing the target system.

        :param forcefield_generator:
            A defined MMForceFieldGenerator object.

        :returns:
            A list of molecular structures obtained from the dynamics simulation.
        """

        all_structures = None

        if self.dynamics_method == 'MM':
            
            # define OpenMMDriver object and perform a dynamical sampling
            openmmdyn = OpenMMDynamics()

            openmmdyn.create_system_from_molecule(molecule, ff_gen=forcefield_generator, 
                                                  solvent=self.dynamics_settings['solvent'], 
                                                  qm_atoms='all')
    
            _, conformation_structures = openmmdyn.conformational_sampling(ensemble=self.ensemble,
                                                                           snapshots=self.snapshots, 
                                                                           nsteps=self.timestep, 
                                                                           temperature=self.temperature, 
                                                                           minimize=True)

            all_structures = conformation_structures
            

        if self.dynamics_method == 'IM':
            interpolation_driver = InterpolationDriver()
            interpolation_driver.update_settings(self.interpolation_settings)
            interpolation_driver.imforcefield_file = self.imforcefieldfile
            qmlabels, z_matrix = interpolation_driver.read_labels()
            sorted_labels = sorted(qmlabels, key=lambda x: int(x.split('_')[1]))

            interpolation_driver.impes_coordinate.z_matrix = z_matrix

            qm_data_points = []

            for label in sorted_labels:
                print('label', label)
                qm_data_point = InterpolationDatapoint(z_matrix)
                qm_data_point.read_hdf5(self.imforcefieldfile, label)
                qm_data_points.append(qm_data_point)   
            interpolation_driver.labels = sorted_labels
            interpolation_driver.qm_data_points = qm_data_points
            openmmdyn = OpenMMDynamics()
            openmmdyn.create_system_from_molecule(molecule, ff_gen=forcefield_generator, solvent=self.solvent, qm_atoms='all')
            openmmdyn.platform = 'CUDA'

            openmmdyn.run_qmmm(interpolation_driver, interpolation_driver, ensemble=self.ensemble, temperature=self.temperature,
                               pressure=self.pressure, friction=self.friction, timestep=self.timestep, nsteps=self.nsteps,
                               snapshots=self.snapshots, traj_file='Database_quality_conformation.pdb')
            
            all_structures = openmmdyn.dynamic_molecules

        if self.dynamics_method == 'IM_Driver':
            interpolation_driver = InterpolationDriver()
            interpolation_driver.update_settings(self.interpolation_settings)
            interpolation_driver.imforcefield_file = self.imforcefieldfile
            self.qmlabels, z_matrix = interpolation_driver.read_labels()

            interpolation_driver.impes_coordinate.z_matrix = z_matrix
            im_database_driver = IMDatabasePointCollecter()
            im_database_driver.platform = 'CUDA'

            im_database_driver.system_from_molecule(molecule, z_matrix, forcefield_generator, solvent=self.solvent, qm_atoms='all')
            self.interpolation_settings['imforcefield_file'] = self.imforcefieldfile
            im_database_driver.update_settings(self.dynamics_settings, self.interpolation_settings)
            im_database_driver.collect_qm_points = self.nsteps
            im_database_driver.run_qmmm()

        return all_structures


    def confirm_database_quality(self, molecule, im_database_file, given_molecular_strucutres=None):
        """Validates the quality of an interpolation database for a given molecule.

       This function assesses the quality of the provided interpolation database 
       comparing the interpolated energy with a QM-reference energy.

       :param molecule:
           A VeloxChem molecule object representing the reference molecular system.

       :param im_database_file:
           Interpolation database file.

       :param given_molecular_strucutres:
           An optional list of additional molecular structures that will be used for the validation.

       :returns:
           List of IM-energies, QM-energies.
        """

        # For all Methods a ForceField of the molecule is requiered
        forcefield_generator = MMForceFieldGenerator()
        forcefield_generator.force_field_data = self.FF_datafile
        forcefield_generator.partial_charges = molecule.get_partial_charges(molecule.get_charge())
        forcefield_generator.create_topology(molecule)
        
        self.imforcefieldfile = im_database_file
        if self.interpolation_settings is None:
            self.interpolation_settings = { 'interpolation_type':self.interpolation_type, 
                            'exponent_p':self.exponent_p,
                            'exponent_q':self.exponent_q, 
                            'confidence_radius':self.confidence_radius,
                            'imforcefield_file':self.imforcefieldfile,
                            'use_inverse_bond_length':True
                          }

        if self.dynamics_method == 'MM':
            rot_bonds = forcefield_generator.rotatable_bonds
            forcefield_generator.reparameterize_dihedrals(rot_bonds[0], scan_range=[180, 360], n_points=7, visualize=True)

        database_quality = False


        while database_quality is False:
            
            database_expanded = False
            all_structures = given_molecular_strucutres
            if given_molecular_strucutres is None:
                all_structures = self.simple_run_dynamics(molecule, forcefield_generator)

            datapoint_molecules = self.database_extracter(im_database_file, molecule.get_labels())
                
            rmsd = -np.inf
            random_structure_choices = None
            counter = 0
            # if given_molecular_strucutres is not None:
            #     random_structure_choices = given_molecular_strucutres

    
            
            while rmsd < 0.3 and counter <= 20:
                
                individual_distances = []
                random_structure_choices = random.choices(all_structures, k=self.nstruc_to_confirm_database_quality)
                for datapoint_molecule in datapoint_molecules:
                    for random_struc in random_structure_choices:
                        
                        distance_norm = self.calculate_distance_to_ref(random_struc.get_coordinates_in_bohr(), datapoint_molecule.get_coordinates_in_bohr())
                        individual_distances.append(distance_norm / np.sqrt(len(molecule.get_labels())) * bohr_in_angstrom())
                
                rmsd = min(individual_distances)

                counter += 1
                if rmsd >= 0.3:
                    print(f'The overall RMSD is {rmsd} -> The current structures are well seperated from the database conformations! loop is discontinued')
                else:
                    print(f'The overall RMSD is {rmsd} -> The current structures are not all well seperated from the database conformations! loop is continued')        
            
            
            qm_energies = []
            im_energies = []
            qm_driver = self.qm_driver
            basis = MolecularBasis.read(molecule, self.basis_set_label)
            impes_driver = InterpolationDriver()
            impes_driver.update_settings(self.interpolation_settings)
            impes_driver.imforcefield_file = im_database_file
            labels, z_matrix = impes_driver.read_labels()
            impes_driver.impes_coordinate.z_matrix = z_matrix
            sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))
            print('labels used for the interpolation', sorted_labels)
            for i, mol in enumerate(random_structure_choices):

                impes_driver.compute(mol, labels=sorted_labels)
                reference_energy = self.compute_energy(qm_driver, mol, basis)
                qm_energies.append(reference_energy[0])
                im_energies.append(impes_driver.impes_coordinate.energy)
                
                print('Energies', qm_energies[-1], im_energies[-1])
                
                print(f'\n\n ########## Step {i} ######### \n')
                print(f'delta_E:   {abs(qm_energies[-1] - im_energies[-1]) * hartree_in_kcalpermol()} kcal/mol \n')
                if abs(qm_energies[-1] - im_energies[-1]) * hartree_in_kcalpermol() > self.energy_threshold:
                    print(mol.get_xyz_string())
                    labels = []
                    labels.append("point_{0}".format((len(sorted_labels) + 1)))
                    self.add_point(mol, im_database_file)
                    database_expanded = True
                    print('The interpolation quality was to low! Structre as been added to the database')

            if not database_expanded:
                database_quality = True

        
        self.plot_final_energies(qm_energies, im_energies)
        
        self.structures_to_xyz_file(random_structure_choices, 'random_xyz_traj.xyz', im_energies, qm_energies)



        return qm_energies, im_energies
    
    def plot_final_energies(self, qm_energies, im_energies):
        """Plots the final potential energies of QM and IM methods.

        :param qm_energies:
            A list or NumPy array of QM potential energies.

        :param im_energies:
            A list or NumPy array of IM potential energies.

        """
        
        import matplotlib.pyplot as plt
        
        qm_energies_plot = np.array(qm_energies)
        im_energies_org = np.array(im_energies)


        sort_indices = np.argsort(qm_energies_plot) 
        qm_energies_sorted = qm_energies_plot[sort_indices]  
        im_energies_org_sorted = im_energies_org[sort_indices]

        ss_total = np.sum((qm_energies_sorted - np.mean(qm_energies_sorted)) ** 2)
        r2_org = 1 - np.sum((qm_energies_sorted - im_energies_org_sorted) ** 2) / ss_total

        plt.figure(figsize=(10, 6))

        plt.plot(qm_energies_sorted - qm_energies_sorted[0], qm_energies_sorted - qm_energies_sorted[0], label="Reference (y = x)", color="black", linestyle="-")
        plt.scatter(qm_energies_sorted - qm_energies_sorted[0], qm_energies_sorted - qm_energies_sorted[0], label="QM Energy", color="black", alpha=1.0)

        plt.scatter(qm_energies_sorted - qm_energies_sorted[0], im_energies_org_sorted - qm_energies_sorted[0], label=f"IM Energy (R$^2$ = {r2_org})", color="red", alpha=0.7)
        plt.xlabel("Potential Energies (kcal mol$^{-1}$)", fontsize=12)
        plt.ylabel("Potential Energies (kcal mol$^{-1}$)", fontsize=12)
        plt.legend()
        plt.grid(alpha=0.5)

        plt.savefig("correlation_energy_plot.svg")
        plt.show(block=True)

    def structures_to_xyz_file(self, molecules_for_xyz, structure_filename, im_energies=None, qm_energies=None):
        """Writes molecular structures to an XYZ file.

        :param molecules_for_xyz:
            A list of VeloxChem molecular objects.

        :param structure_filename:
            The name of the output file where XYZ structures will be stored.

        :param im_energies:
            An optional list of interpolation energies corresponding to each molecule.

        :param qm_energies:
            An optional list of quantum mechanical energies corresponding to each molecule.

        """

        with open(structure_filename, 'w') as file:
            pass

        for i, dyn_mol in enumerate(molecules_for_xyz):

            current_xyz_string = dyn_mol.get_xyz_string()

            xyz_lines = current_xyz_string.splitlines()

            if len(xyz_lines) >= 2 and im_energies is not None:

                xyz_lines[1] += f'Energies  QM: {qm_energies[i]}  IM: {im_energies[i]}  delta_E: {abs(qm_energies[i] - im_energies[i])}'


            updated_xyz_string = "\n".join(xyz_lines)

            with open(structure_filename, 'a') as file:
                file.write(f"{updated_xyz_string}\n\n")
        
    
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
    

    def add_point(self, molecule, imforcefielddatafile):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param imforcefielddatafile:
                Datafile containing the information of the IM forcefield.

        """
        if self.qm_driver is None:
            raise ValueError("No energy driver defined.")
        if self.qm_grad_driver is None:
            raise ValueError("No gradient driver defined.")
        if self.qm_hess_driver is None:
            raise ValueError("No Hessian driver defined.")

        energy = None
        gradient = None
        hessian = None
           
        # define impesdriver to determine if stucture should be added:
        interpolation_driver = InterpolationDriver()
        interpolation_driver.update_settings(self.interpolation_settings)
        interpolation_driver.imforcefield_file = imforcefielddatafile
        labels, z_matrix = interpolation_driver.read_labels()
        interpolation_driver.impes_coordinate.z_matrix = z_matrix
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        interpolation_driver.compute(molecule, labels=sorted_labels)

        basis = MolecularBasis.read(molecule, self.basis_set_label)

        energy, scf_results = self.compute_energy(self.qm_driver, molecule, basis)
        
        # if (abs(energy[0] - interpolation_driver.impes_coordinate.energy)) * hartree_in_kcalpermol() < 0.95:
        #     assert_msg_critical((abs(energy[0] - interpolation_driver.impes_coordinate.energy)) * hartree_in_kcalpermol(), 
        #                         'The potential energy for the given structure is already approximated within chemical accuracy (less then 1.0 kcal/mol)')
        
        print((abs(energy[0] - interpolation_driver.impes_coordinate.energy)) * hartree_in_kcalpermol()) 

        gradient = self.compute_gradient(self.qm_grad_driver, molecule, basis, scf_results)
        hessian = self.compute_hessian(self.qm_hess_driver, molecule, basis)

        natoms = molecule.number_of_atoms()
        elem = molecule.get_labels()
        coords = molecule.get_coordinates_in_bohr().reshape(natoms * 3)

        vib_frequencies, normal_modes_vec, gibbs_energy = (
            geometric.normal_modes.frequency_analysis(
                coords,
                hessian[0],
                elem,
                energy=energy[0],
                temperature=self.temperature,
                pressure=self.pressure,
                outfnm=f'vibrational_point_{energy[0]}',
                normalized=False))
    
        print('Vibrational frequencies', vib_frequencies, gibbs_energy)

        impes_coordinate = InterpolationDatapoint(z_matrix)
        impes_coordinate.update_settings(self.interpolation_settings)
        impes_coordinate.cartesian_coordinates = molecule.get_coordinates_in_bohr()
        impes_coordinate.energy = energy[0]
        impes_coordinate.gradient = gradient[0]
        impes_coordinate.hessian = hessian[0]
        impes_coordinate.transform_gradient_and_hessian()

        impes_coordinate.write_hdf5(imforcefielddatafile, f'point_{len(sorted_labels) + 1}')
        
        print(f"Database expansion with {', '.join(labels)}")
        for l, e in zip(labels, energy):
            print(f"{l}: Energy = {e:.4f} hatree")

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

        # restricted SCF
        elif isinstance(qm_driver, ScfRestrictedDriver):
            qm_driver.ostream.mute()
            scf_tensors = qm_driver.compute(molecule, basis)
            qm_energy = qm_driver.scf_energy
            qm_driver.ostream.unmute()

        if qm_energy is None:
            error_txt = "Could not compute the QM energy. "
            error_txt += "Please define a QM driver."
            raise ValueError(error_txt)

        return qm_energy, scf_tensors

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

        elif isinstance(grad_driver, ScfGradientDriver):
            grad_driver.ostream.mute()
            grad_driver.compute(molecule, basis, scf_results)
            qm_gradient = grad_driver.gradient
            qm_gradient = np.array([qm_gradient])
            grad_driver.ostream.unmute()

        if qm_gradient is None:
            error_txt = "Could not compute the QM gradient. "
            error_txt += "Please define a QM gradient driver."
            raise ValueError(error_txt)
        
        return qm_gradient


    # TODO: mute outside to save time?
    def compute_hessian(self, hess_driver, molecule, basis=None):
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

        elif isinstance(hess_driver, ScfHessianDriver):
            hess_driver.ostream.mute()
            hess_driver.compute(molecule, basis)
            qm_hessian = hess_driver.hessian
            qm_hessian = np.array([qm_hessian])
            hess_driver.ostream.unmute()

        if qm_hessian is None:
            error_txt = "Could not compute the QM Hessian. "
            error_txt += "Please define a QM Hessian driver."
            raise ValueError(error_txt)


        return qm_hessian

 


