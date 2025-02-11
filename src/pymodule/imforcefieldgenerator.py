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
from pathlib import Path
from sys import stdout
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom

from contextlib import redirect_stderr
from io import StringIO
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .imdatabasepointcollecter import IMDatabasePointCollecter
# from .impesdatabasebuilder import ImpesDatabaseBuilder
# from .imdatabasedriver import IMDatabaseDriver
#from .impesforcefieldgenerator_parallel import ImpesForceFieldGeneratorParallel
from .forcefieldgenerator import ForceFieldGenerator
from .atomtypeidentifier import AtomTypeIdentifier
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from .molecule import Molecule
from .errorhandler import assert_msg_critical

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

        - impes_dict: A dictionary containing settings for the interpolation.

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
        self.qm_driver = None
        self.qm_grad_driver = None
        self.qm_hess_driver = None

        # variables for the interpolation
        self.impes_dict = None
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
        self.temperature = '150.15'
        self.pressure = 1.0
        self.force_constant = 1.0
        self.ensemble = 'NVE'
        self.timestep = 0.5
        self.nsteps = 1000
        self.snapshots = self.nsteps
        self.trajectory_file = 'trajectory.pdb'
        self.desired_point_density = 50
        self.converged_cycle = 4
        self.energy_threshold = 1.5
        self.start_collect = 0
        self.solvent = 'gas'

        self.FF_datafile = '/home/vlind06/phd_work/interpolation/gaff-2.11.dat'

        # individual run
        self.qm_energies = []
        self.total_energies = []
        self.molecules = []
        self.kinetic_energies = []
        self.point_added_molecules = []
        self.unique_molecules = []


    def set_up_the_system(self, molecule, target_dihedrals=None, sampling_structures=1):

        """
        Assign the neccessary variables with respected values. 

        :param molecule: original molecule

        :param target_dihedrals: is a list of dihedrals that should be scanned during the dynamics

        :param sampling_structures: devides the searchspace around given rotatbale dihedrals
            
        """


        self.qm_data_points = None
        self.molecule = molecule

        if self.z_matrix is None:
            self.z_matrix = self.define_z_matrix(molecule)
        
        # Read the database to determine the  
        
        if self.imforcefieldfile is not None:

            impes_driver = InterpolationDriver(self.z_matrix)
            impes_driver.update_settings(self.impes_dict)
            self.qmlabels = impes_driver.read_labels()

            self.qm_data_points = []
            self.qm_energies = []
            for label in self.qmlabels:
                qm_data_point = InterpolationDatapoint(self.z_matrix)
                qm_data_point.read_hdf5(self.imforcefieldfile, label)
                self.qm_energies.append(qm_data_point.energy)
                self.qm_data_points.append(qm_data_point)        
        
        self.density_of_datapoints, self.molecules_along_rp, self.allowed_deviation = self.database_density_check_with_molecule(molecule, self.qm_data_points, specific_dihedrals=target_dihedrals, nsampling=sampling_structures)

        self.impes_dict['dihedrals'] = target_dihedrals

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
        if self.qm_driver is None:
            assert_msg_critical(self.qm_driver is not None, 'ImForceFieldGenerator: No QM-Driver/ QM-Gradient-Driver / QM-Hessian-Driver were initialized!ß.')
        
        print('MOlecule labels', molecule.get_labels())
        self.impes_dict = { 'interpolation_type':self.interpolation_type, 
                            'exponent_p':self.exponent_p,
                            'exponent_q':self.exponent_q, 
                            'confidence_radius':self.confidence_radius,
                            'imforcefield_file':self.imforcefieldfile,
                          }

        self.dynamics_settings = {  'qm_driver': self.qm_driver,
                                    'grad_driver': self.qm_grad_driver,
                                    'hess_driver': self.qm_hess_driver,
                                    'basis_set':self.basis_set_label,
                                    'duration':self.duration, 'temperature':self.temperature, 'solvent':self.solvent,
                                    'pressure':self.force_constant, 'force_constant': self.force_constant, 'ensemble':self.ensemble,
                                    'timestep': self.timestep, 'nsteps': self.nsteps,
                                    'snapshots':self.snapshots, 'trajectory_file':self.trajectory_file,
                                    'desired_datapoint_density':self.desired_point_density, 'converged_cycle': self.converged_cycle, 'energy_threshold':self.energy_threshold,
                                    'NAC':False, 'load_system': None, 'collect_qm_points_from':self.start_collect,
                                    'FF_datafile':self.FF_datafile}
        

        self.set_up_the_system(molecule, self.dihedrals, self.sampling_structures)


        for counter, entry in enumerate(self.molecules_along_rp.items()):
            key, molecules = entry
            for mol in molecules:
                
                if self.dihedrals is not None:
                    key = (key, int(round(mol.get_dihedral_in_degrees(key))))
                    print('key', key, self.density_of_datapoints)

                forcefield_generator = ForceFieldGenerator()
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
                    im_database_driver.allowed_molecule_deviation = [key[1] - self.allowed_deviation, -key[1] + self.allowed_deviation]
                    im_database_driver.update_settings(self.dynamics_settings, self.impes_dict)
                    if self.impes_dict['imforcefield_file'] is None:
                        self.impes_dict['imforcefield_file'] = im_database_driver.qm_datafile

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

        Returns:

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

        if qm_datapoints and specific_dihedrals is not None:
            impes_driver = InterpolationDriver(self.z_matrix)
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


        allowed_deviation = 360 / nsampling
        normalized_angle = ((allowed_deviation + 180) % 360) - 180

        return point_densities, sampled_molecules, normalized_angle
    
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

 
