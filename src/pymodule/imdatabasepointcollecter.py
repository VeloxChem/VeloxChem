#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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
from scipy.optimize import minimize
import scipy
import h5py
import itertools
import re
import os
from pathlib import Path
from sys import stdout
import sys
import random
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom
from contextlib import redirect_stderr
from io import StringIO

from .molecule import Molecule
from .veloxchemlib import mpi_master
from. veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .solvationbuilder import SolvationBuilder
from .optimizationdriver import OptimizationDriver

# Drivers
from .scfrestdriver import ScfRestrictedDriver
from .molecularbasis import MolecularBasis
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .externalscfdriver import ExternalScfDriver
from .externalgradientdriver import ExternalGradientDriver
from .externalhessiandriver import ExternalHessianDriver
from .externalexcitedstatedriver import ExternalExcitedStatesScfDriver
from .externalexcitedstategradientdriver import ExternalExcitedStatesGradientDriver
from .externalexcitedstatehessiandriver import ExternalExcitedStatesHessianDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .localbayesresidual import LocalBayesResidual

with redirect_stderr(StringIO()) as fg_err:
    import geometric

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass


class IMDatabasePointCollecter:
    """
    Implements the OpenMMDynamics.
    Performs classical MD and QM/MM simulations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - platform: The platform for OpenMM. Default is 'reference'
                    Options are 'CPU', 'CUDA' and 'OpenCL'.
        - ensemble: The thermodynamic ensemble used in the simulation.
                    Options are: 'NVE', 'NVT' and 'NPT'.
        - temperature: The temperature.
        - friction: The friction coefficient which couples the system to the heat bath.
        - timestep: The timestep for the integrator in fs.
        - nsteps: The number of steps.
        - parent_ff: The Force Field to be used in standard residues.
        - water_ff: The Force Field to be used for water.
        - box_size: The size of the box for periodic simulations.
        - padding: The padding for the box in nm.
        - cutoff: The cutoff for nonbonded interactions in nm.
        - integrator: The integrator object.
        - system: The OpenMM system object.
        - topology: The OpenMM topology object.
        - simulation: The OpenMM simulation object.
        - constraints: The constraints for the system.
        - pdb: The PDB file object.
        - modeller: The modeller object.
        - labels: The atom labels.
        - molecule: The VeloxChem molecule object.
        - unique_residues: The list of unique residues in the system.
        - unique_molecules: The list of unique molecules in the system.
        - qm_driver: The VeloxChem driver object. Options are XtbDriver and InterpolationDriver.
        - grad_driver: The VeloxChem gradient driver object.
        - qm_atoms: The list of atom indices for the QM region.
        - mm_subregion: The list of atom indices for the MM subregion (part of a molecule).
        - linking_atoms: The list of atom indices for the linking atoms.
        - qm_force_index: The index of the QM force in the system.
        - driver_flag: The flag to determine the driver to be used.
        - scaling_factor: The scaling factor for the QM stabilizer.
        - linking_atom_distance: The distance between the QM and MM regions in angstroms.
    """
    
    def __init__(self, comm=None, ostream=None):
        """
        Initializes the class with default simulation parameters.
        """
        np.set_printoptions(threshold=sys.maxsize)
        # MPI and output stream
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # Instance variables
        # Simulation parameters
        self.platform = None
        self.ensemble = None
        self.temperature = None
        self.friction = None
        self.timestep = None
        self.nsteps = None
        self.parent_ff = 'amber03.xml'
        self.water_ff = 'spce.xml'
        self.box_size = 2.0 
        self.padding = 1.0
        self.cutoff = 1.0
        self.scaling_factor = 0.0
        self.integrator = None

        # OpenMM objects
        self.system = None
        self.topology = None
        self.integrator = None
        self.simulation = None
        self.constraints = None
        self.pdb = None
        self.modeller = None
        self.labels = []

        # VeloxChem objects
        self.molecule = None
        self.ghost_atom = (False, None)
        self.energy_gabs = {}
        self.state_energies = {}
        self.unique_residues = []
        self.unique_molecules = []

        self.current_im_choice = None
        
        # QM Region parameters
        self.drivers = None
        self.qm_atoms = None
        self.mm_subregion = None
        self.linking_atoms = None
        self.qm_force_index = None
        self.driver_flag = None
        self.swapped = False
        self.state_swtiched = False
        self.velo_switch = False
        self.excitation_pulse = None

        self.roots_to_follow = []

        self.snapshots = None
        self.load_system = None
        self.pressure = None
        self.output_file = None
        self.adiabatic_basis = False
        self.density_around_data_point = None
        self.impes_drivers = None
        self.im_labels = None
        self.sorted_im_labels = []
        self.qm_energies = None
        self.qm_data_points = None
        self.point_adding_molecule = {}
        self.energy_threshold = None
        self.collect_qm_points = None
        self.previous_energy_list = []
        self.non_core_symmetry_groups = []
        self.interpolation_settings = None
        self.allowed_molecules = None
        self.reference_struc_energies_file = None

        self.starting_temperature = None

        self.current_state = None
        self.starting_state = 0
        self.distance_thrsh = 0.1
        self.current_im_choice = None
        self.current_gradient = 0
        self.current_energy = 0
        self.point_checker = 1
        self.allowed_molecule_deviation = None
        self.last_point_added = None
        self.im_labels = None
        self.add_a_point = False
        self.check_a_point = False
        self.cluster_run = None
        self.skipping_value = 0
        self.basis_set_label = None
        self.molecule = None
        self.expansion = True
        self.expansion_molecules = []
        self.dynamics_settings_interpolation_run = None
        self.sampled_molecules = None
        self.bayes_models = None


        # output_file variables that will be written into self.general_variable_output
        self.gradients = None
        self.velocities = None
        self.start_velocities = None
        self.coordinates = None
        self.coordinates_xyz = None
        self.state_specific_molecules = None
        self.all_gradients = []
        
        self.identfy_relevant_int_coordinates = True
        self.use_symmetry = True
        self.use_opt_confidence_radius = True
        self.qm_symmetry_dict = None
        self.add_bayes_model = True
        
        self.summary_output = 'summary_output.h5'
        self.coordinates_xyz = None

        self.velocities = []

        # Default value for the C-H linker distance
        self.linking_atom_distance = 1.0705 
            
    # Method to generate OpenMM system from VeloxChem objects
    def system_from_molecule(self,
                             molecule,
                             z_matrix, 
                             ff_gen, 
                             solvent='gas', 
                             qm_atoms=None, 
                             filename='residue',
                             trust_radius=False, 
                             residue_name='MOL'):
        """
        Generates an OpenMM system from a VeloxChem molecule and a forcefield generator.
        
        :param molecule:
            VeloxChem molecule object.
        :param ff_gen:
            VeloxChem forcefield generator object.
        :param solvent:
            Available options:'gas', 'spce', 'tip3p', 'ethanol', 'methanol', 'acetone', 
            'chloroform', 'hexane', 'toluene', 'dcm', 'benzene', 'dmso', 'thf', 
            'acetonitrile', 'other' or 'itself'.
        :param qm_atoms:
            Options: None, 'all', or list of atom indices for QM region.
        :param filename:
            Base name for the generated files. Default is 'residue'.
        :param residue_name:
            Name of the residue. Default is 'MOL'.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Store the molecule object and generate OpenMM compatible files
        self.molecule = molecule
        self.positions = molecule.get_coordinates_in_angstrom()
        self.labels = molecule.get_labels()
        self.z_matrix = z_matrix

        # Options for the QM region if it's required.
        # TODO: Take this if else tree to a separate method.
        if qm_atoms:
            filename = 'qm_region'
            if trust_radius:
                filename += 'trust_radius'
            # Create a QM region topology template
            if qm_atoms == 'all':
                msg = 'Full molecule as QM region'
                self.ostream.print_info(msg)
                self.ostream.flush()
                qm_atoms = list(range(len(self.labels)))
                self._create_QM_residue(ff_gen, 
                                        qm_atoms, 
                                        filename)

            elif isinstance (qm_atoms, list):
                if qm_atoms == list(range(len(self.labels))):
                    msg = 'Full molecule as QM region'
                    self.ostream.print_info(msg)
                    self.ostream.flush()
                    self._create_QM_residue(ff_gen,
                                            qm_atoms,
                                            filename)
                msg = 'QM/MM partition inside the molecule'
                self.ostream.print_info(msg)
                self.ostream.flush()
                self._create_QM_subregion(ff_gen,
                                          qm_atoms, 
                                          molecule, 
                                          filename)

            ff_gen.write_pdb(f'{filename}.pdb', 'QMR')

        elif qm_atoms is None:
            ff_gen.write_openmm_files(filename, residue_name)

        # TODO: Incorporate the SystemBuilder here to avoid using the Modeller.
        if solvent == 'gas':
            phase = 'gas'
            self.pdb = app.PDBFile(f'{filename}.pdb')     
            # Common forcefield loading, modified according to phase specifics
            forcefield_files = [f'{filename}.xml']

        if solvent != 'gas':
            # Solvate the molecule using the SystemBuilder
            phase = 'periodic'
            sys_builder = SolvationBuilder()
            sys_builder.solvate(solute=molecule, 
                                solvent=solvent,
                                padding=self.padding,
                                equilibrate=True
                                )
            #sys_builder.write_openmm_files(solute_ff=ff_gen)
            if solvent == 'spce':
                self.water_ff = 'spce.xml'
                forcefield_files = [f'{filename}.xml', self.parent_ff, self.water_ff]
            elif solvent == 'tip3p':
                self.water_ff = 'tip3p.xml'
                forcefield_files = [f'{filename}.xml', self.parent_ff, self.water_ff]
            elif solvent == 'itself':
                # Solute and solvent are the same
                forcefield_files = [f'{filename}.xml', self.parent_ff]
            else:
                # Any other solvent, including custom ones
                solvent_ff = 'solvent_1.xml'
                forcefield_files = [f'{filename}.xml', self.parent_ff, solvent_ff]

            # Load the PDB from the SystemBuilder
            self.pdb = app.PDBFile('equilibrated_system.pdb')

        # Create the ForceField object        
        forcefield = app.ForceField(*forcefield_files)

        # Create the System object
        if phase == 'gas':
            self.system = forcefield.createSystem(self.pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)

        else:
            self.system = forcefield.createSystem(self.pdb.topology, 
                                                  nonbondedMethod=app.PME, 
                                                  nonbondedCutoff=self.cutoff * unit.nanometer, 
                                                  constraints=app.HBonds)
        
        # Modify the system to include QM and QM/MM forces.
        if qm_atoms:
            self.set_qm_mm_system(phase, 
                                  ff_gen)
            # self.qm_stabilizer(ff_gen)
        
        # Write the system to a xml file (for debugging purposes)
        with open(f'{filename}_system.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(self.system))
            msg = f'System parameters written to {filename}_system.xml'
            self.ostream.print_info(msg)
            self.ostream.flush()

        # Write the system to a pdb file (for debugging purposes)
        if phase == 'gas':
            app.PDBFile.writeFile(self.pdb.topology, 
                                  self.positions, 
                                  open(f'{filename}_system.pdb', 'w'))
            msg = f'System coordinates written to {filename}_system.pdb'
            self.ostream.print_info(msg)
            self.ostream.flush()

        elif phase == 'periodic':
            app.PDBFile.writeFile(self.pdb.topology, 
                                  self.pdb.positions, 
                                  open(f'{filename}_system.pdb', 'w'))
            msg = f'System coordinates written to {filename}_system.pdb'
            self.ostream.print_info(msg)
            self.ostream.flush()
        
        self.phase = phase

    def add_bias_force(self, atoms, force_constant, target):
        """
        Method to add a biasing force to the system.

        :param atoms:
            List of atom indices to apply the force.
        :param force_constant:
            Force constant for the biasing force.
        :param target:
            Target equilibrium parameter (distance (nm), angle and torsion (deg))
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        if len(atoms) == 2:
            msg = f'Adding stretch force between atoms {atoms[0]} and {atoms[1]} with force constant {force_constant}.'
            self.ostream.print_info(msg)
            self.ostream.flush()
            force = mm.CustomBondForce('0.5*k*(r-r0)^2')
            force.addGlobalParameter('k', force_constant)
            force.addGlobalParameter('r0', target)
            force.addBond(atoms[0], atoms[1])
            self.system.addForce(force)
        elif len(atoms) == 3:
            target_rad = target * np.pi / 180
            msg = f'Adding bend force between atoms {atoms[0]}, {atoms[1]}, and {atoms[2]} with force constant {force_constant}.'
            self.ostream.print_info(msg)
            self.ostream.flush()
            force = mm.CustomAngleForce('0.5*k*(theta-theta0)^2')
            force.addGlobalParameter('k', force_constant)
            force.addGlobalParameter('theta0', target_rad)
            force.addAngle(atoms[0], atoms[1], atoms[2])
            self.system.addForce(force)
        elif len(atoms) == 4:
            target_rad = target * np.pi / 180
            msg = f'Adding torsion force between atoms {atoms[0]}, {atoms[1]}, {atoms[2]}, and {atoms[3]} with force constant {force_constant}.'        
            self.ostream.print_info(msg)
            self.ostream.flush()
            force = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
            force.addGlobalParameter('k', force_constant)
            force.addGlobalParameter('theta0', target_rad)
            force.addTorsion(atoms[0], atoms[1], atoms[2], atoms[3])
            self.system.addForce(force)
        else:
            raise ValueError('Invalid number of atoms for the biasing force.')
     
    def conformational_sampling(self, 
                                ensemble='NVT', 
                                temperature=700, 
                                timestep=2.0, 
                                nsteps=10000, 
                                snapshots=10,
                                minimize=True):

        """
        Runs a high-temperature MD simulation to sample conformations and minimize the energy of these conformations.

        :param ensemble:
            Type of ensemble. Options are 'NVE', 'NVT', 'NPT'. Default is 'NVT'.
        :param temperature:
            Temperature of the system in Kelvin. Default is 700 K.
        :param timestep:
            Timestep of the simulation in femtoseconds. Default is 2.0 fs.
        :param nsteps:
            Number of steps in the simulation. Default is 1000.
        :param snapshots:
            The number of snapshots to save. Default is 10.


        :return:
            Tuple containing the minimized potential energies and the XYZ format strings of the relaxed coordinates.

        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        if self.system is None:
            raise RuntimeError('System has not been created!')
        if self.molecule is None:
            raise RuntimeError('Molecule object does not exist!')

        self.ensemble = ensemble
        self.temperature = temperature * unit.kelvin

        self.timestep = timestep * unit.femtosecond
        self.nsteps = nsteps

        self.integrator = self._create_integrator()
        topology = self.modeller.topology if self.phase in ['water', 'periodic'] else self.pdb.topology
        self.positions = self.modeller.positions if self.phase in ['water', 'periodic'] else self.pdb.positions

        self.simulation = app.Simulation(topology, self.system, self.integrator)
        self.simulation.context.setPositions(self.positions)
        self.simulation.context.setVelocitiesToTemperature(self.temperature)

        save_freq = nsteps // snapshots if snapshots else nsteps
        energies = []
        opt_coordinates = []

        for step in range(nsteps):
            self.simulation.step(1)
            if step % save_freq == 0:
                state = self.simulation.context.getState(getPositions=True, getEnergy=True)
                energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                coordinates = state.getPositions()
                self.simulation.context.setPositions(coordinates) 

                if minimize:

                    print(f'Step: {step}, Potential energy: {energy}')
                    self.simulation.minimizeEnergy()
                    minimized_state = self.simulation.context.getState(getPositions=True, getEnergy=True)
                    minimized_energy = minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                    minimized_coordinates = minimized_state.getPositions()

                    energies.append(minimized_energy)
                    print(f'Minimized energy: {minimized_energy}')
                    xyz = f"{len(self.labels)}\n\n"
                    for label, coord in zip(self.labels, minimized_coordinates):
                        xyz += f"{label} {coord.x * 10} {coord.y * 10} {coord.z * 10}\n"  # Convert nm to Angstroms
                    print('Saved coordinates for step', step)
                    opt_coordinates.append(xyz)

                # Create a xyz with the coordinates
                print('Energy of the conformation:', energy)
                xyz = f"{len(self.labels)}\n\n"
                for label, coord in zip(self.labels, coordinates):
                    xyz += f"{label} {coord.x * 10} {coord.y * 10} {coord.z * 10}\n"
                opt_coordinates.append(xyz)
                energies.append(energy)
                print('Saved coordinates for step', step)

        print('Conformational sampling completed!')
        print(f'Number of conformations: {len(opt_coordinates)}')

        return energies, opt_coordinates

    def sort_points_with_association(self, points, associated_list):
        '''
        This function is sorting the points from a given database
        starting from the first entry to the last (point_0,..., point_n)

        :param points:
            list of database points
        :param associated_list:
            list of labels that are associated to the individual datapoints

        '''
        
        # Define a custom key function to extract the numerical part
        def extract_number(point):
            return int(point.split('_')[1])

        # Combine the points and associated list into tuples
        combined = list(zip(points, associated_list))

        # Sort the combined list using the custom key function
        sorted_combined = sorted(combined, key=lambda x: extract_number(x[0]))

        # Unzip the sorted combined list back into two lists
        sorted_points, sorted_associated_list = zip(*sorted_combined)

        return list(sorted_points), list(sorted_associated_list)
    
    
    def update_settings(self, dynamics_settings, interpolation_settings=None):
        """
        Updates settings in the ImpesDynamicsDriver.
        :param molecule:
            The Molecule object.
        :param impes_dict:
            The input dictionary of settings for IMPES.
        """
        
        self.drivers = dynamics_settings['drivers']

        self.interpolation_settings = interpolation_settings
        self.dynamics_settings_interpolation_run = dynamics_settings


        if 'print_step' in dynamics_settings:
            self.print_step = int(dynamics_settings['print_step'])

        if 'duration' in dynamics_settings:
            self.duration = float(dynamics_settings['duration'])

        if 'temperature' in dynamics_settings:
            self.temperature = dynamics_settings['temperature']
            self.starting_temperature = dynamics_settings['temperature']

        if 'pressure' in dynamics_settings:
            self.pressure = float(dynamics_settings['pressure'])

        if 'force_constant' in dynamics_settings:
            self.force_constant = float(dynamics_settings['force_constant'])

        if 'friction' in dynamics_settings:
            self.friction = float(dynamics_settings['friction'])

        if 'timestep' in dynamics_settings:
            self.timestep = float(dynamics_settings['timestep'])

        if 'nsteps' in dynamics_settings:
            self.nsteps = int(dynamics_settings['nsteps'])

        if 'snapshots' in dynamics_settings:
            self.snapshots = int(dynamics_settings['snapshots'])

        # Start the simulation form a given state
        if 'load_system' in dynamics_settings:
            self.load_system = dynamics_settings['load_system']

        if 'trajectory_file' in dynamics_settings:
            self.out_file = dynamics_settings['trajectory_file']
        else:
            self.out_file = 'trajectory.pdb'

        if 'reference_struc_energy_file' in dynamics_settings:
            self.reference_struc_energies_file = dynamics_settings['reference_struc_energy_file']

        # Determines the ensemble in order to set the correct simulation set_up
        if 'ensemble' in dynamics_settings:
            self.ensemble = dynamics_settings['ensemble']

        #################################### DATABASE construciton inputs #############################

        # The desired density around a given starting structure/datapoint
        if 'desired_datapoint_density' in dynamics_settings:
            self.desired_datpoint_density = int(dynamics_settings['desired_datapoint_density'])

        # The number of iteration cycles without database expansion after which the description around the reference point has converged
        if 'converged_cycle' in dynamics_settings:
            self.unadded_cycles = int(dynamics_settings['converged_cycle'])
        
        if 'basis_set_label' in dynamics_settings:
            basis_set_label = dynamics_settings['basis_set_label']
            self.basis_set_label = basis_set_label
    
        # Dertermines if non-adiabatic couplings should be calculated
        if 'NAC' in dynamics_settings:
            self.calc_NAC = dynamics_settings['NAC']
        else:
            self.calc_NAC = False

        if 'ADT' in dynamics_settings:
            self.adiabatic_basis = True

        if 'cluster_run' in dynamics_settings:
            self.cluster_run = dynamics_settings['cluster_run']
        # TODO: these need to be different dictionary entries
        if 'collect_qm_points_from' in dynamics_settings:
            self.collect_qm_points = dynamics_settings['collect_qm_points_from']

        # time step at which QM data point collection should stop
        if 'qmc_stop' in dynamics_settings:
            self.qmc_stop = float(dynamics_settings["qmc_stop"])

        # index of the excited state (in case of TDDFT QM data points)
        if "roots_to_follow" in dynamics_settings:
            self.nstates = len(list(dynamics_settings["roots_to_follow"]))
            self.roots_to_follow = list(dynamics_settings['roots_to_follow'])

        if 'excitation_pulse' in dynamics_settings:
            self.excitation_pulse = list(dynamics_settings['excitation_pulse'])
        
        # keywords used to reverse the dynamics
        # to a previous iteration.
        if 'reverse' in dynamics_settings:
            self.reverse = True
        else:
            self.reverse = False

        # how many iterations to reverse
        if 'reverse_by' in dynamics_settings:
            self.reverse_by = int(dynamics_settings['reverse_by'])
            self.reverse = True
        
        if 'energy_threshold' in dynamics_settings:
            self.energy_threshold = dynamics_settings['energy_threshold']

        # threshold to determine if a bond is breaking
        if 'bond_threshold' in dynamics_settings:
            self.bond_threshold = float(dynamics_settings['bond_threshold'])

    
    
    def run_qmmm(self):
        """
        Runs a QM/MM simulation using OpenMM, storing the trajectory and simulation data.

        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        if self.system is None:
            raise RuntimeError('System has not been created!')
        
        # temperature = self.temperature
        self.temperature_number = self.temperature
        self.temperature = self.temperature * unit.kelvin
        
        # friction = self.friction
        self.friction = self.friction / unit.picosecond
        
        timestep = self.timestep
        self.timestep = timestep * unit.femtoseconds

        self.qm_potentials = []
        self.qm_mm_interaction_energies = []
        self.mm_potentials = []
        self.total_potentials = []
        self.kinetic_energies = []
        self.temperatures = []
        self.total_energies = []
        self.coordinates_xyz = []

        save_freq = self.nsteps // self.snapshots if self.snapshots else self.nsteps

        # Create or update the integrator
        if self.integrator is None:
            new_integrator = self._create_integrator()
        else:
            new_integrator = self.integrator

        self.topology = self.pdb.topology

        self.positions = self.pdb.positions
        
        self.simulation = app.Simulation(self.topology, self.system, new_integrator)

        # Load the state if a restart file is provided
        if self.load_system is not None:
            self.simulation.loadState(self.load_system)
        
        else:
            self.simulation.context.setPositions(self.positions)

            # Set initial velocities if the ensemble is NVT or NPT
            if self.ensemble in ['NVT', 'NPT']:
                self.simulation.context.setVelocitiesToTemperature(self.temperature)
        
        self.start_velocities = self.simulation.context.getState(getVelocities=True).getVelocities()
        # Set up reporting
        self.simulation.reporters.clear()
        self.simulation.reporters.append(app.PDBReporter(self.out_file, save_freq))

        # Print header
        print('QM/MM Simulation Parameters')
        print('=' * 60)
        print('QM Driver:', self.driver_flag)
        print('Ensemble:', self.ensemble)
        print('Integration method:', new_integrator.__class__.__name__)
        if self.ensemble in ['NVT', 'NPT']:
            print('Temperature:', self.temperature, 'K')
            if self.ensemble == 'NPT':
                print('Pressure:', self.pressure, 'atm')
        print('Friction:', self.friction, '1/ps')
        print('Timestep:', self.timestep, 'fs')
        print('Total simulation time in ns:', self.nsteps * timestep / 1e6)
        print('=' * 60)

        self.last_point_added = 0
        self.cycle_iteration = self.unadded_cycles
        
        print('current datapoints around the given starting structure', self.desired_datpoint_density, self.density_around_data_point[0], self.density_around_data_point[1], '\n allowed derivation from the given structure', self.allowed_molecule_deviation, '\n ---------------------------------------')
        self.allowed_molecules = {root: {'molecules': [], 'qm_energies': [], 'im_energies':[]} for root in self.roots_to_follow}
        openmm_coordinate = self.simulation.context.getState(getPositions=True).getPositions()
        self.coordinates = [openmm_coordinate]
        self.coordinates_xyz = []
        self.gradients = []
        self.velocities = []
        self.velocities_np = []
        self.velocities_np.append(self.simulation.context.getState(getVelocities=True).getVelocities(True))

        self.impes_drivers = []
        
        state_spec_dict = {'pot_energies':[], 'gradients':[]}


        self.gloabal_sim_informations = {f'state_{root}':state_spec_dict for root in self.roots_to_follow}
        self.gloabal_sim_informations['coordinates_ang'] = []
        self.gloabal_sim_informations['state'] = []
        self.gloabal_sim_informations['temperatures'] = []

        self.state_specific_molecules = {root: [] for root in self.roots_to_follow}
        self.sorted_state_spec_im_labels = {root: [] for root in self.roots_to_follow}

        self.bayes_models = {root: [] for root in self.roots_to_follow}
        self.qm_data_point_dict = {root: [] for root in self.roots_to_follow}
        self.qm_symmetry_datapoint_dict = {root: {} for root in self.roots_to_follow}
        self.qm_energies_dict = {root: [] for root in self.roots_to_follow}
        self.impes_drivers = {root: None for root in self.roots_to_follow}
        self.sampled_molecules = {root: {'molecules': [], 'im_energies': [], 'qm_energies': [], 'distances': []} for root in self.roots_to_follow}
        print(self.reference_struc_energies_file)
        last_added = 0
        if self.reference_struc_energies_file is not None:
            self.sampled_molecules = self.extract_reference_structures(self.reference_struc_energies_file, self.roots_to_follow)
            self.allowed_molecules = self.extract_reference_structures(self.reference_struc_energies_file, self.roots_to_follow)
            last_added += len(self.sampled_molecules[self.roots_to_follow[0]]['molecules'])
        else:
             self.reference_struc_energies_file = 'ref_struc_energy.xyz'

        for root in self.roots_to_follow:
            # Dynamically create an attribute name
            attribute_name = f'impes_driver_{root}'
            # Initialize the object
            driver_object = InterpolationDriver(self.z_matrix)
            driver_object.update_settings(self.interpolation_settings[root])
            if root == 0:
                driver_object.symmetry_information = self.non_core_symmetry_groups['gs']
            else:
                driver_object.symmetry_information = self.non_core_symmetry_groups['es']
            driver_object.use_symmetry = self.use_symmetry

            im_labels, _ = driver_object.read_labels()
            print('beginning labels', im_labels, self.add_bayes_model)
            self.qm_data_points = []
            self.qm_energies = []
            old_label = None

            for label in im_labels:
                if '_symmetry' not in label:
                    qm_data_point = InterpolationDatapoint(self.z_matrix)
                    qm_data_point.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)

                    self.qm_energies_dict[root].append(qm_data_point.energy)
                    self.qm_data_point_dict[root].append(qm_data_point)

                    self.sorted_state_spec_im_labels[root].append(label)
                    
                    if self.add_bayes_model:
                        bayes_model = LocalBayesResidual(self.z_matrix)
                        bayes_model.compute_internal_coordinates_values(qm_data_point.cartesian_coordinates)
                    
                        if root == 0:
                            bayes_model.symmetry_information = self.non_core_symmetry_groups['gs']
                        else:
                            bayes_model.symmetry_information = self.non_core_symmetry_groups['es']
                        
                        bayes_model.init_bayesian()

                        self.bayes_models[root].append(bayes_model)
                    
                    old_label = qm_data_point.point_label
                    driver_object.qm_symmetry_data_points[old_label] = [qm_data_point]
                    self.qm_symmetry_datapoint_dict[root][old_label] = [qm_data_point]

                else:
                    symmetry_data_point = InterpolationDatapoint(self.z_matrix)
                    symmetry_data_point.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)
                    
                    driver_object.qm_symmetry_data_points[old_label].append(symmetry_data_point)
                    self.qm_symmetry_datapoint_dict[root][old_label].append(symmetry_data_point)
                    
                    # driver_object.qm_symmetry_data_points_1 = {old_label: [self.qm_data_point_dict[root][2], self.qm_data_point_dict[root][4]]}
                    # driver_object.qm_symmetry_data_points_2 = {old_label: [self.qm_data_point_dict[root][3], self.qm_data_point_dict[root][5]]}

            print(driver_object.qm_symmetry_data_points)       
            # Set the object as an attribute of the instance
            setattr(self, attribute_name, driver_object)
            # Append the object to the list
            self.impes_drivers[root] = driver_object

        self.current_state = self.starting_state

        start_time = time()
        self.step = 0 
        for step in range(self.nsteps):

            self.update_forces(self.simulation.context)
            
            # Potential energies
            # QM region
            qm = self.get_qm_potential_energy()
            self.qm_potentials.append(qm)

            # QM/MM interactions
            qm_mm = self.simulation.context.getState(getEnergy=True, groups={1,2}).getPotentialEnergy()
            self.qm_mm_interaction_energies.append(qm_mm.value_in_unit(unit.kilojoules_per_mole))

            # MM region 
            mm = self.simulation.context.getState(getEnergy=True, groups={3,4,5,6,7}).getPotentialEnergy()
            self.mm_potentials.append(mm.value_in_unit(unit.kilojoules_per_mole))

            # Total potential energy
            pot = qm * unit.kilojoules_per_mole + qm_mm + mm
            self.total_potentials.append(pot.value_in_unit(unit.kilojoules_per_mole))

            # Kinetic energy
            kinetic = self.simulation.context.getState(getEnergy=True).getKineticEnergy()
            self.kinetic_energies.append(kinetic.value_in_unit(unit.kilojoules_per_mole))

            # Temperature
            dof = self.system.getNumParticles() * 3 - self.system.getNumConstraints()
            if hasattr(self.integrator, 'computeSystemTemperature'):
                temp = self.integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
            else:
                temp = (2 * kinetic / (dof * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)
            self.temperatures.append(temp)
    
            # Total energy
            total = pot + kinetic
            self.total_energies.append(total.value_in_unit(unit.kilojoules_per_mole))

            self.gloabal_sim_informations['temperatures'].append(temp)
            self.gloabal_sim_informations['state'].append(self.current_state)
            self.gloabal_sim_informations['coordinates_ang'].append(self.coordinates_xyz[-1])

            # Information output
            if step % save_freq == 0:
                
                print(f"Step: {step} / {self.nsteps} Time: {round((step * timestep) / 1000, 2)} ps")
                print('Potential Energy QM region:', qm, 'kJ/mol')
                print('Potential Energy MM region:', mm)
                print('QM/MM Interaction Energy:', qm_mm)
                print('Total Potential Energy:', pot)
                print('Kinetic Energy:', kinetic)
                print('Temperature:', temp, 'K')
                print('Total Energy:', total, '±', np.std(self.total_energies), 'kJ/mol')
                for root in self.roots_to_follow:
                    print('Current Density', self.density_around_data_point[0][root], '-->', self.desired_datpoint_density, self.unadded_cycles)   
                print('Current State (PES):', self.current_state) 
                print('-' * 60)   
            # self.output_file_writer(self.summary_output)
            self.step += 1

            self.simulation.step(1)
            
            if step % 100 == 0 and step != 0:
                self.simulation.saveCheckpoint('checkpoint')
                
                # self.output_file_writer(self.summary_output)
                state_spec_dict = {'pot_energies':[], 'gradients':[]}
                self.gloabal_sim_informations = {f'state_{root}':state_spec_dict for root in self.roots_to_follow}
                self.gloabal_sim_informations['coordinates_ang'] = []
                self.gloabal_sim_informations['state'] = []
                self.gloabal_sim_informations['temperatures'] = []

                for root in self.roots_to_follow:

                    self.write_qm_energy_determined_points(self.sampled_molecules[root]['molecules'][last_added: ],
                                                           self.sampled_molecules[root]['qm_energies'][last_added:],
                                                           self.sampled_molecules[root]['im_energies'][last_added:],
                                                           [],
                                                           root)
                    last_added = len(self.sampled_molecules[root]['molecules'])
                #print('cooridnates', simulation.context.getState(getPositions=True).getPositions())

            # if step == self.nsteps and self.density_around_data_point[0] != self.desired_datpoint_density:
            #     step = 0
            
            for root in self.roots_to_follow:
                if self.density_around_data_point[0][root] >= self.desired_datpoint_density and self.expansion:
                    step = self.nsteps
                    break
            
            if not self.expansion:
                if len(self.expansion_molecules) >= 20:
                    break

            if self.unadded_cycles == 0:
                step = self.nsteps
                break

        end_time = time()
        elapsed_time = end_time - start_time
        elapsed_time_days = elapsed_time / (24 * 3600)
        performance = (self.nsteps * timestep / 1e6) / elapsed_time_days

        print('QM/MM simulation completed!')
        print(f'Number of steps: {self.nsteps}')
        print(f'Trajectory saved as {self.out_file}')
        print('=' * 60)
        print('Simulation Averages:')
        print('=' * 60)
        print('QM Potential Energy:', np.mean(self.qm_potentials), '±', np.std(self.qm_potentials), 'kJ/mol')
        print('QM/MM Interaction Energy:', np.mean(self.qm_mm_interaction_energies), '±', np.std(self.qm_mm_interaction_energies), 'kJ/mol')
        print('MM Potential Energy:', np.mean(self.mm_potentials), '±', np.std(self.mm_potentials), 'kJ/mol')
        print('Total Potential Energy:', np.mean(self.total_potentials), '±', np.std(self.total_potentials), 'kJ/mol')
        print('Kinetic Energy:', np.mean(self.kinetic_energies), '±', np.std(self.kinetic_energies), 'kJ/mol')
        print('Temperature:', np.mean(self.temperatures), '±', np.std(self.temperatures), 'K')
        print('Total Energy:', np.mean(self.total_energies), '±', np.std(self.total_energies), 'kJ/mol')
        print('=' * 60)
        print(f'Elapsed time: {int(elapsed_time // 60)} minutes, {int(elapsed_time % 60)} seconds')
        print(f'Performance: {performance:.2f} ns/day')
        print(f'Trajectory saved as {self.out_file}')


    # Post-simulation analysis methods
    def plot_energy(self,
                    components=['total'],
                    filename='energy_plot.png'):
        """
        Plots the energy components of the simulation.

        :param components:
            List of energy components to plot. Default is ['total'].
            It can include 'total', 'kinetic', 'potential', 'temperature', 'qm', 'mm', 'qm_mm'.
        :param filename:
            Name of the output file. Default is 'energy_plot.png'.
        """
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot also the trend of the energy components as a linear regression
        if 'total' in components:
            ax.plot(self.total_energies, label='Total Energy', marker='o')
            # Linear regression excluding the first 10% of the data
            x = np.arange(len(self.total_energies))
            m, b = np.polyfit(x[int(0.1 * len(x)):], self.total_energies[int(0.1 * len(x)):], 1)
            ax.plot(x, m * x + b, label=f'Total Energy Trend ({m:.4e} kJ/mol/step)', linestyle='--')

        if 'kinetic' in components:
            ax.plot(self.kinetic_energies, label='Kinetic Energy', marker='o')

        if 'potential' in components:
            ax.plot(self.total_potentials, label='Potential Energy', marker='o')
        if 'temperature' in components:
            ax.plot(self.temperatures, label='Temperature', marker='o')
            ax.set_ylabel('Temperature (K)')
        if 'qm' in components:
            ax.plot(self.qm_potentials, label='QM Energy', marker='o')
        if 'mm' in components:
            ax.plot(self.mm_potentials, label='MM Energy', marker='o')
        if 'qm_mm' in components:
            ax.plot(self.qm_mm_interaction_energies, label='QM/MM Interaction Energy', marker='o')

        ax.set_title('Energy Components of the Simulation')

        ax.set_xlabel('Step')
        ax.set_ylabel('Energy (kJ/mol)')
        ax.legend()
        plt.tight_layout()
        plt.savefig(filename)
        plt.show()
        
    def visualize_trajectory(self, 
                             trajectory_file='trajectory.pdb', 
                             interval=1):
        """
        Visualizes the trajectory of the simulation using py3Dmol.

        :param trajectory_file:
            Path to the PDB file containing the trajectory. Default is 'trajectory.pdb'.
        """
        try:
            import py3Dmol

        except ImportError:
            raise ImportError("py3Dmol is not installed. Please install it using `pip install py3Dmol`.")
        
        viewer = py3Dmol.view(width=800, height=600)

        viewer.addModelsAsFrames(open(trajectory_file, 'r').read(),'pdb', {'keepH': True})
        viewer.animate({'interval': interval, 'loop': "forward", 'reps': 10})
        viewer.setStyle({"stick":{},"sphere": {"scale":0.25}})
        viewer.zoomTo()

        viewer.show()

    # Private methods
    def _save_output(self, output_file):
        """
        Save the simulation output to an out file.
        Reports all the energy components for every snapshot.

        :param output_file:
            Path to the output file.
        """

        with open(output_file + '.out', 'w', encoding='utf-8') as out:
            # Write the header
            out.write('VeloxChem/OpenMM Simulation Output\n')
            out.write('=' * 60 + '\n')
            if self.qm_atoms is not None:
                out.write(f'Simulation type: {self.driver_flag}\n')
            else:
                out.write('Simulation type: Classical MD\n')
            out.write(f'Ensemble: {self.ensemble}\n')
            out.write(f'Temperature: {self.temperature}\n')
            out.write(f'Friction: {self.friction}\n')
            out.write(f'Timestep: {self.timestep} \n')
            out.write(f'Number of steps: {self.nsteps}\n')
            out.write('=' * 60 + '\n')

            # Write the energy components from the lists
            for i in range(len(self.total_energies)):
                out.write(f'Step: {i}\n')
                out.write(f'Potential Energy (kJ/mol): {self.total_potentials[i]:.4f}\n')
                out.write(f'Kinetic Energy (kJ/mol): {self.kinetic_energies[i]:.4f}\n')
                out.write(f'Temperature (K): {self.temperatures[i]:.4f}\n')
                if self.qm_atoms is not None:
                    out.write(f'QM Potential Energy (kJ/mol): {self.qm_potentials[i]:.4f}\n')
                    out.write(f'QM/MM Interaction Energy (kJ/mol): {self.qm_mm_interaction_energies[i]:.4f}\n')
                    out.write(f'MM Potential Energy (kJ/mol): {self.mm_potentials[i]:.4f}\n')
                out.write(f'Total Energy (kJ/mol): {self.total_energies[i]:.4f}\n')
                out.write('-' * 60 + '\n')
            # Write the averages
            out.write('Summary\n')
            out.write('=' * 60 + '\n')
            out.write(f'Average Potential Energy (kJ/mol): {np.mean(self.total_potentials):.4f} ± {np.std(self.total_potentials):.4f}\n')
            out.write(f'Average Kinetic Energy (kJ/mol): {np.mean(self.kinetic_energies):.4f} ± {np.std(self.kinetic_energies):.4f}\n')
            out.write(f'Average Temperature (K): {np.mean(self.temperatures):.4f} ± {np.std(self.temperatures):.4f}\n')
            if self.qm_atoms is not None:
                out.write(f'Average QM Potential Energy (kJ/mol): {np.mean(self.qm_potentials):.4f} ± {np.std(self.qm_potentials):.4f}\n')
                out.write(f'Average QM/MM Interaction Energy (kJ/mol): {np.mean(self.qm_mm_interaction_energies):.4f} ± {np.std(self.qm_mm_interaction_energies):.4f}\n')
                out.write(f'Average MM Potential Energy (kJ/mol): {np.mean(self.mm_potentials):.4f} ± {np.std(self.mm_potentials):.4f}\n')
            out.write(f'Average Total Energy (kJ/mol): {np.mean(self.total_energies):.4f} ± {np.std(self.total_energies):.4f}\n')
            out.write('=' * 60 + '\n')      

    def _format_PDB_file(self, filename):
        """
        Check the format of the PDB file and correct it if necessary.

        :param pdb_file:
            Path to the PDB file to be checked.
        :return:
            Path to the corrected PDB file.
        """

        pdb_path = Path(filename)
        if not pdb_path.is_file():
            raise FileNotFoundError(f"{filename} does not exist.")
        
        pdbstr = pdb_path.read_text()

        # Formatting issues flags
        label_guess_warning = False
        conect_warning = False

        for line in pdbstr.strip().splitlines():
            if line.startswith(('ATOM', 'HETATM')):
                atom_label = line[76:78].strip()
                if not atom_label:
                    label_guess_warning = True

            # Check if the CONECT records are present and skip them
            # If they are not present activate the warning flag
            if line.startswith('CONECT'):
                conect_warning = False
            else:
                conect_warning = True

        # Overwrite the PDB file with the guessed atom labels in columns 77-78
        # if the labels are missing
        if label_guess_warning:

            with open(filename, 'w') as f:
                for line in pdbstr.strip().splitlines():
                    if line.startswith(('ATOM', 'HETATM')):
                        atom_name = line[12:16].strip()
                        atom_label = atom_name[0]
                        new_line = line[:76] + f"{atom_label:>2}" + line[78:] + '\n'
                        f.write(new_line)
                    else:
                        f.write(line + '\n')
            
            print('Warning: Atom labels were guessed based on atom names (first character).')
            print(f'Please verify the atom labels in the {filename} PDB file.')

        if conect_warning:

            # Create a molecule from the PDB file
            molecule = Molecule.read_PDB_file(filename)
            connectivity_matrix = molecule.get_connectivity_matrix()
            # Determine all the bonds in the molecule
            with open(filename, 'a') as f:
                for i in range(connectivity_matrix.shape[0]):
                    for j in range(i + 1, connectivity_matrix.shape[1]):
                        if connectivity_matrix[i, j] == 1:
                            # Convert indices to 1-based index for PDB format and ensure proper column alignment
                            i_index = i + 1
                            j_index = j + 1
                            # Align to the right 
                            con_string = "{:6s}{:>5d}{:>5d}".format('CONECT', i_index, j_index)
                            f.write(con_string + '\n')
       
            print('The CONECT records were not found in the PDB file.')
            print('The connectivity matrix was used to determine the bonds.')

    def _create_integrator(self):
        """
        Creates an OpenMM integrator object based on the specified ensemble type.

        Returns:
            OpenMM Integrator: Configured integrator for the simulation.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Common parameters for Langevin integrators

        if self.ensemble in ['NVT', 'NPT']:
            integrator = mm.LangevinIntegrator(self.temperature, self.friction, self.timestep)
            integrator.setConstraintTolerance(1e-5)
            if self.ensemble == 'NPT':
                # Add a barostat for pressure control in NPT ensemble
                barostat = mm.MonteCarloBarostat(1 * unit.atmospheres, self.temperature, 25)
                self.system.addForce(barostat)
        elif self.ensemble == 'NVE':
            integrator = mm.VerletIntegrator(self.timestep)

        else:
            raise ValueError("Unsupported ensemble type. Please choose 'NVE', 'NVT', or 'NPT'.")

        return integrator
    
    # Methods to create a QM region/subregion in the system
    def _create_QM_residue(self, 
                           ff_gen,
                           qm_atoms, 
                           filename='qm_region', 
                           residue_name='QMR'):
        """
        This method creates an xml file for a QM region.
        The xml file only contains atomtypes, residue, and nonbonded parameters.

        :param ff_gen:
            VeloxChem forcefield generator object
        :param qm_atoms:
            List of atom indices in the QM region.
        :param filename:
            Name of the files to be generated. Default is 'qm_region'
        :param residue_name:
            Name of the residue. Default is 'QMR'
        """
        self.qm_atoms = qm_atoms
        atoms = ff_gen.atoms
        bonds = ff_gen.bonds

        # Create the root element of the XML file
        ForceField = ET.Element("ForceField")
        
        # AtomTypes section
        AtomTypes = ET.SubElement(ForceField, "AtomTypes")

        for i, atom in atoms.items():
            element = ''.join([i for i in atom['name'] if not i.isdigit()])  
            attributes = {
                # Name is the atom type_molname
                "name": atom['name'] + '_' + residue_name,
                "class": str(i + 1),
                "element": element,
                "mass": str(atom['mass']) 
            }
            ET.SubElement(AtomTypes, "Type", **attributes)

        # Residues section
        Residues = ET.SubElement(ForceField, "Residues")
        Residue = ET.SubElement(Residues, "Residue", name=residue_name)
        for atom_id, atom_data in atoms.items():
            ET.SubElement(Residue, "Atom", name=atom_data['name'], type=atom_data['name'] + '_' + residue_name, charge=str(atom_data['charge']))
        for bond_id, bond_data in bonds.items():
            ET.SubElement(Residue, "Bond", atomName1=atoms[bond_id[0]]['name'], atomName2=atoms[bond_id[1]]['name'])

        # NonbondedForce section
        NonbondedForce = ET.SubElement(ForceField, "NonbondedForce", coulomb14scale=str(ff_gen.fudgeQQ), lj14scale=str(ff_gen.fudgeLJ))
        for atom_id, atom_data in atoms.items():
            attributes = {
                "type": atom_data['name'] + '_' + residue_name,
                "charge": str(0.0),
                "sigma": str(0.0),
                "epsilon": str(0.0)
            }
            ET.SubElement(NonbondedForce, "Atom", **attributes)

        # Generate the tree and write to file
        tree = ET.ElementTree(ForceField)
        rough_string = ET.tostring(ForceField, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        indented_string = reparsed.toprettyxml(indent="    ")  

        with open(filename + '.xml', 'w') as output_file:
            output_file.write(indented_string)

        self.ostream.print_info(f'QM region parameters written to {filename}.xml')
        self.ostream.flush()
        
    def _create_QM_subregion(self, ff_gen, qm_atoms, molecule, filename='qm_region', residue_name='QMR'):
        """
        Creates the xml file for a molecule with a QM region.

        :param ff_gen:
            VeloxChem forcefield generator object
        :param qm_atoms:
            List of atom indices in the QM region.
        :param molecule:
            VeloxChem molecule object
        :param filename:
            Name of the files to be generated. Default is 'qm_subregion'
        :param residue_name:
            Name of the residue. Default is 'QMR'
        """

        # Atoms belonging to the QM region do not have bonded parameters
        # Atoms belonging to the MM region have bonded parameters

        # List of atom indices in the molecule
        atom_indices = list(range(molecule.number_of_atoms()))

        connectivity_matrix = molecule.get_connectivity_matrix()

        # Get the atoms connected to the QM region
        qm_connected = []
        for i in qm_atoms:
            qm_connected.extend(np.where(connectivity_matrix[i])[0])
        qm_connected = list(set(qm_connected) - set(qm_atoms))
        self.linking_atoms = qm_connected
        print('Linking atoms indices:', qm_connected)

        # Get the broken bonds between the QM and MM regions
        self.broken_bonds = []
        for i in qm_connected:
            for j in qm_atoms:
                if connectivity_matrix[i, j]:
                    self.broken_bonds.append((i, j))
        # Set with the atoms in broken bonds
        self.broken_bonds_atoms = list(set([atom for bond in self.broken_bonds for atom in bond]))

        # The linked atoms are part of the QM region since a customforce is used.
        mm_atoms = list(set(atom_indices) - set(qm_atoms) - set(qm_connected))
        # QM atoms are the specified QM region + the linking atoms
        qm_atoms.extend(qm_connected)
        # Order the atoms
        qm_atoms.sort()
        self.qm_atoms = qm_atoms
        print('QM subregion (self.qm_atoms):', self.qm_atoms)
        self.mm_subregion = mm_atoms
        print('MM subregion:(self.mm_subregion)', self.mm_subregion)

        atoms = ff_gen.atoms
        bonds = ff_gen.bonds
        angles = ff_gen.angles
        dihedrals = ff_gen.dihedrals
        impropers = ff_gen.impropers

        # Create the root element of the XML file
        ForceField = ET.Element("ForceField")

        # AtomTypes section
        AtomTypes = ET.SubElement(ForceField, "AtomTypes")

        for i, atom in atoms.items():
            element = ''.join([i for i in atom['name'] if not i.isdigit()])  
            attributes = {
                # Name is the atom type_molname
                "name": atom['name'] + '_' + residue_name,
                "class": str(i + 1),
                "element": element,
                "mass": str(atom['mass']) 
            }
            ET.SubElement(AtomTypes, "Type", **attributes)

        # Residues section
        Residues = ET.SubElement(ForceField, "Residues")
        Residue = ET.SubElement(Residues, "Residue", name=residue_name)
        for atom_id, atom_data in atoms.items():
            ET.SubElement(Residue, "Atom", name=atom_data['name'], type=atom_data['name'] + '_' + residue_name, charge=str(atom_data['charge']))
        for bond_id, bond_data in bonds.items():
            ET.SubElement(Residue, "Bond", atomName1=atoms[bond_id[0]]['name'], atomName2=atoms[bond_id[1]]['name'])

        # NonbondedForce section
        NonbondedForce = ET.SubElement(ForceField, "NonbondedForce", coulomb14scale=str(ff_gen.fudgeQQ), lj14scale=str(ff_gen.fudgeLJ))
        for atom_id, atom_data in atoms.items():
            if atom_id in qm_atoms:
                charge = 0.0
                sigma = 0.0
                epsilon = 0.0
            else:
                charge = atom_data['charge']
                sigma = atom_data['sigma']
                epsilon = atom_data['epsilon']
            attributes = {
                "type": atom_data['name'] + '_' + residue_name,
                "charge": str(charge),
                "sigma": str(sigma),
                "epsilon": str(epsilon)
            }
            ET.SubElement(NonbondedForce, "Atom", **attributes)

        bonded_atoms = mm_atoms + self.broken_bonds_atoms

        # BondForce section
        BondForce = ET.SubElement(ForceField, "HarmonicBondForce")
        for bond_id, bond_data in bonds.items():
            if (bond_id[0] in bonded_atoms and 
                bond_id[1] in bonded_atoms):
                attributes = {
                    "class1": str(bond_id[0] + 1),
                    "class2": str(bond_id[1] + 1),
                    "length": str(bond_data['equilibrium']),
                    "k": str(bond_data['force_constant'])
                }
                ET.SubElement(BondForce, "Bond", **attributes)

        # AngleForce section
        AngleForce = ET.SubElement(ForceField, "HarmonicAngleForce")
        for angle_id, angle_data in angles.items():
            if (angle_id[0] in bonded_atoms and 
                angle_id[1] in bonded_atoms and 
                angle_id[2] in bonded_atoms):
                attributes = {
                    "class1": str(angle_id[0] + 1),
                    "class2": str(angle_id[1] + 1),
                    "class3": str(angle_id[2] + 1),
                    "angle": str(angle_data['equilibrium'] * np.pi / 180),
                    "k": str(angle_data['force_constant'])
                }
                ET.SubElement(AngleForce, "Angle", **attributes)

        # DihedralForce section
        DihedralForce = ET.SubElement(ForceField, "PeriodicTorsionForce")
        for dihedral_id, dihedral_data in dihedrals.items():
            if (dihedral_id[0] in bonded_atoms and
                dihedral_id[1] in bonded_atoms and
                dihedral_id[2] in bonded_atoms and
                dihedral_id[3] in bonded_atoms):
                # Skip RB dihedrals
                if dihedral_data['type'] == 'RB':
                    continue
                attributes = {
                    "class1": str(dihedral_id[0] + 1),
                    "class2": str(dihedral_id[1] + 1),
                    "class3": str(dihedral_id[2] + 1),
                    "class4": str(dihedral_id[3] + 1),
                    "periodicity1": str(dihedral_data['periodicity']),
                    "phase1": str(dihedral_data['phase'] * np.pi / 180),
                    "k1": str(dihedral_data['barrier'])
                }
                ET.SubElement(DihedralForce, "Proper", **attributes)

        # RB Dihedrals section
        RBForce = ET.SubElement(ForceField, "RBTorsionForce")
        for dihedral_id, dihedral_data in dihedrals.items():
            if (dihedral_id[0] in bonded_atoms and
                dihedral_id[1] in bonded_atoms and
                dihedral_id[2] in bonded_atoms and
                dihedral_id[3] in bonded_atoms):
                # Skip Fourier dihedrals
                if dihedral_data['type'] == 'Fourier':
                    continue
                attributes = {
                    "class1": str(dihedral_id[0] + 1),
                    "class2": str(dihedral_id[1] + 1),
                    "class3": str(dihedral_id[2] + 1),
                    "class4": str(dihedral_id[3] + 1),
                    "c0": str(dihedral_data['RB_coefficients'][0]),
                    "c1": str(dihedral_data['RB_coefficients'][1]),
                    "c2": str(dihedral_data['RB_coefficients'][2]),
                    "c3": str(dihedral_data['RB_coefficients'][3]),
                    "c4": str(dihedral_data['RB_coefficients'][4]),
                    "c5": str(dihedral_data['RB_coefficients'][5])
                }
                ET.SubElement(RBForce, "Torsion", **attributes)

        # ImproperForce section
        ImproperForce = ET.SubElement(ForceField, "PeriodicTorsionForce")
        for improper_id, improper_data in impropers.items():
            if (improper_id[0] in bonded_atoms and
                improper_id[1] in bonded_atoms and
                improper_id[2] in bonded_atoms and
                improper_id[3] in bonded_atoms):
                attributes = {
                    "class1": str(improper_id[1] + 1),
                    "class2": str(improper_id[0] + 1),
                    "class3": str(improper_id[2] + 1),
                    "class4": str(improper_id[3] + 1),
                    "periodicity1": str(improper_data['periodicity']),
                    "phase1": str(improper_data['phase'] * np.pi / 180),
                    "k1": str(improper_data['barrier'])
                }
                ET.SubElement(ImproperForce, "Torsion", **attributes)

        # Generate the tree and write to file
        # tree = ET.ElementTree(ForceField)
        rough_string = ET.tostring(ForceField, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        indented_string = reparsed.toprettyxml(indent="    ")  

        with open(filename + '.xml', 'w') as output_file:
            output_file.write(indented_string)

        msg = f'QM subregion parameters written to {filename}.xml'
        self.ostream.print_info(msg)
        self.ostream.flush()

    # Method to set a qm_mm system:
    def set_qm_mm_system(self, phase, ff_gen):
        """
        Configures the system for QM/MM calculations.

        :param qm_atoms:
            List of atom indices to be included in the QM region.
        :param phase:
            Phase of the system ('gas', 'water', 'periodic').
        :param ff_gen:
            MMForceFieldGenerator object from VeloxChem.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Set the QM/MM Interaction Groups
        total_atoms = self.system.getNumParticles()
        
        # The MM subregion is counted as regular MM atoms
        qm_group = set(self.qm_atoms)
        print('QM Group:', qm_group)
        mm_group = set(range(total_atoms)) - qm_group
        print('MM Group:', mm_group)
        if not mm_group:
            print('No external MM atoms found in the system')

        # Add custom hessian forces
        force_expression = mm.CustomExternalForce("-fx*x-fy*y-fz*z")
        self.system.addForce(force_expression)
        force_expression.addPerParticleParameter("fx")
        force_expression.addPerParticleParameter("fy")
        force_expression.addPerParticleParameter("fz")

        for i in self.qm_atoms:
            force_expression.addParticle(i, [0, 0, 0])

        # QM Hessian Force Group
        force_expression.setForceGroup(0)

        # If a MM region is present define the interactions
        if mm_group:
            
            nonbonded_force = None
            for force in self.system.getForces():
                if isinstance(force, mm.NonbondedForce):
                    nonbonded_force = force
                    break
    
            if nonbonded_force is None:
                raise ValueError("NonbondedForce not found in the system")

            # CustomNonbondedForce for QM/MM interactions
            vdw = mm.CustomNonbondedForce("4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
            if phase in ['water', 'periodic']:
                vdw.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
                vdw.setCutoffDistance(self.cutoff)
            vdw.addPerParticleParameter("sigma")
            vdw.addPerParticleParameter("epsilon")


            if phase in ['water', 'periodic']:
                # OpenMM uses Reaction Field method for CustomNB with PBC.
                # DEBUG Only for testing purposes
                print('Using Reaction Field method for long-range electrostatics!')
                rfDielectric = nonbonded_force.getReactionFieldDielectric()
                krf = (1 / (self.cutoff**3)) * (rfDielectric - 1) / (2*rfDielectric + 1)
                crf = (1 /  self.cutoff) * (3*rfDielectric) / (2*rfDielectric + 1)
                coulomb_rf = f"(138.935456*charge1*charge2)*(1/r + {krf}*r*r - {crf});"
                coulomb = mm.CustomNonbondedForce(coulomb_rf)
                coulomb.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
                coulomb.setCutoffDistance(self.cutoff)
                # Disable long-range electrostatics for nonbonded force
                nonbonded_force.setUseDispersionCorrection(False)
            else:
                coulomb = mm.CustomNonbondedForce("138.935456*charge1*charge2/r;")
                coulomb.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
            
            coulomb.addPerParticleParameter("charge")
            
            # Apply the same exclusions to the custom forces as in the NonbondedForce
            # This is needed to avoid errors (all forces must have identical exclusions)
            for i in range(nonbonded_force.getNumExceptions()):
                p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
                vdw.addExclusion(p1, p2)
                coulomb.addExclusion(p1, p2)
            
            self.system.addForce(vdw)
            self.system.addForce(coulomb)

            # Add particles to the custom forces
            # QM region
            for i in qm_group:
                vdw.addParticle([ff_gen.atoms[i]['sigma']* unit.nanometer, 
                                 ff_gen.atoms[i]['epsilon']] * unit.kilojoules_per_mole)
                coulomb.addParticle([ff_gen.atoms[i]['charge']] * unit.elementary_charge)

            # MM region
            # Obtain the sigma, epsilon, and charge values from the system
            nonbonded_force = [f for f in self.system.getForces() if isinstance(f, mm.NonbondedForce)][0]

            for i in mm_group:
                # The charges, sigmas, and epsilons are taken from the system
                charge, sigma, epsilon = nonbonded_force.getParticleParameters(i)
                sigma = sigma * unit.nanometer
                epsilon = epsilon * unit.kilojoules_per_mole
                charge = charge * unit.elementary_charge
                vdw.addParticle([sigma, epsilon])
                coulomb.addParticle([charge])

            vdw.addInteractionGroup(qm_group, mm_group)
            coulomb.addInteractionGroup(qm_group, mm_group)

            # Set force groups
            vdw.setForceGroup(1)
            coulomb.setForceGroup(2)

            # Set a force group for the MM region
            for force in self.system.getForces():
                # Non bonded MM
                if isinstance(force, mm.NonbondedForce):
                    force.setForceGroup(3)
                # Bonded
                elif isinstance(force, mm.HarmonicBondForce):
                    force.setForceGroup(4)
                elif isinstance(force, mm.HarmonicAngleForce):
                    force.setForceGroup(5)
                elif isinstance(force, mm.PeriodicTorsionForce):
                    force.setForceGroup(6)
                elif isinstance(force, mm.RBTorsionForce):
                    force.setForceGroup(7)
                elif isinstance(force, mm.CMMotionRemover):
                    force.setForceGroup(8)
                
        # Determine the force index for the QM region
        # it is an instance of openmm.openmm.CustomExternalForce
        for i, force in enumerate(self.system.getForces()):
            if isinstance(force, mm.CustomExternalForce):
                self.qm_force_index = i
                break


    def cartesian_just_distance(self, coordinate_1, coordinate_2, non_core_atoms):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.
           :param data_point:
                InterpolationDatapoint object
        """
        target_coordinates_core = np.delete(coordinate_1, non_core_atoms, axis=0)
        reference_coordinates_core = np.delete(coordinate_2, non_core_atoms, axis=0)
        # First, translate the cartesian coordinates to zero
        target_coordinates = self.calculate_translation_coordinates(target_coordinates_core)
        reference_coordinates = (
            self.calculate_translation_coordinates(reference_coordinates_core))
        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)
        rotation_matrix = geometric.rotate.get_rot(target_coordinates,
                                                   reference_coordinates)
        # Rotate the data point
        rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
        # Calculate the Cartesian distance
        distance_vector = (reference_coordinates - rotated_coordinates)

        return np.linalg.norm(distance_vector)
    
    def qm_stabilizer(self, ff_gen_qm):
        
        """
        Implements a MM potential to stabilize the QM region.
        The forces are 1% of the regular MM forces.

        :param qm_atoms: 
            List of atom indices to be included in the QM region.
        :param ff_gen_qm: 
            MMForceFieldGenerator object from VeloxChem.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Harmonic bond contribution. Parameters are read from ff_gen_qm
        bonds = ff_gen_qm.bonds
        bond_force = mm.HarmonicBondForce()
        for bond, params in bonds.items():
            bond_force.addBond(*bond,
                            params['equilibrium'] * unit.nanometer,
                            params['force_constant'] * unit.kilojoule_per_mole / unit.nanometer**2 * self.scaling_factor)
        self.system.addForce(bond_force)

        # Harmonic angle contribution. Parameters are read from ff_gen_qm
        angles = ff_gen_qm.angles
        angle_force = mm.HarmonicAngleForce()
        for angle, params in angles.items():
            angle_force.addAngle(*angle,
                                params['equilibrium'] * np.pi / 180 * unit.radian,
                                params['force_constant'] * unit.kilojoule_per_mole / unit.radian**2 * self.scaling_factor)
        self.system.addForce(angle_force)

        # Periodic torsion contribution. Parameters are read from ff_gen_qm
        torsions = ff_gen_qm.dihedrals
        torsion_force = mm.PeriodicTorsionForce()
        rb_torsion_force = mm.RBTorsionForce()
        for torsion, params in torsions.items():
            if params['type'] == 'Fourier':
                torsion_force.addTorsion(*torsion,
                                         params['periodicity'],
                                         params['phase'] * np.pi / 180 * unit.radian,
                                         params['barrier'] * unit.kilojoule_per_mole * self.scaling_factor)
            elif params['type'] == 'RB':
                rb_torsion_force.addTorsion(*torsion,
                                            params['RB_coefficients'][0] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][1] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][2] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][3] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][4] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][5] * unit.kilojoule_per_mole * self.scaling_factor)
        self.system.addForce(torsion_force)
        self.system.addForce(rb_torsion_force)

        # Improper torsion contribution. Parameters are read from ff_gen_qm
        impropers = ff_gen_qm.impropers
        improper_force = mm.PeriodicTorsionForce()
        for improper, params in impropers.items():
            improper_force.addTorsion(*improper,
                                    params['periodicity'],
                                    params['phase'] * np.pi / 180 * unit.radian,
                                    params['barrier'] * unit.kilojoule_per_mole * self.scaling_factor)
        self.system.addForce(improper_force)

    def update_gradient_and_energy(self, new_positions):
        """
        Updates and returns the gradient and potential energy of the QM region.

        :param new_positions:
            The new positions of the atoms in the QM region.
        :return:
            The gradient and potential energy of the QM region.
        """

        # self.coordinates_xyz.append(new_positions * 10)
        positions_ang = new_positions * 10 

        new_molecule = None
        # Check if there is a QM/MM partition in the system
        if self.mm_subregion is not None:
            # Create a molecule with a link atom (H)
            # The linking atom is defined in self.linking_atoms
            # It need to be changed to a H atom at 1.0 angstrom
            # The QM region is defined in self.qm_atoms

            # Create a molecule with the new positions
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]

            qm_atom_labels = []
            for i in self.qm_atoms:
                # Change the linking atom to H
                if i in self.linking_atoms:
                    qm_atom_labels.append('H')
                else:
                    qm_atom_labels.append(atom_labels[i])

            
            # Change the positions of the linking atoms to 1.0 angstrom
            for atom1, atom2 in self.broken_bonds:

                current_distance = np.linalg.norm(positions_ang[atom1] - positions_ang[atom2])
                # Change the position of the linking atom to 1.0 angstrom
                # By construction, the first atom is the linking atom
                direction = (positions_ang[atom2] - positions_ang[atom1]) / current_distance
                positions_ang[atom1] = positions_ang[atom2] - direction * self.linking_atom_distance

            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")

        else:
            # Atom labels for the QM region
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            self.unique_molecules.append(new_molecule)
        
        for root in self.roots_to_follow:
            self.impes_drivers[root].qm_data_points = self.qm_data_point_dict[root]
            self.impes_drivers[root].compute(new_molecule)

        transitions = []
        if self.nstates > 1:
            for root_1 in range(0, self.nstates):
                for root_2 in range(root_1 + 1, self.nstates):

                    potential_kjmol = self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
                    potential_kjmol_2 = self.impes_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
                    transitions.append((self.roots_to_follow[root_1], self.roots_to_follow[root_2], potential_kjmol - potential_kjmol_2))
                    print(f'compare the energies between roots: {self.roots_to_follow[root_1]} -> {self.roots_to_follow[root_2]}', potential_kjmol_2 - potential_kjmol)

                    if self.calc_NAC == True and np.linalg.norm(self.velocities_np[-1]) > 0.0:
                        current_NAC = self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.NAC.flatten()
                        current_velocitites = self.velocities_np[-1].flatten() * 4.566180e-4

                        hopping_potential = np.exp(-abs((np.pi/(4)) * (( self.impes_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy - self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy ) / np.linalg.multi_dot([current_NAC, current_velocitites]))))
                        print('#######################', '\n\n', hopping_potential, potential_kjmol_2 - potential_kjmol, '\n\n', '#######################')

                    if abs(potential_kjmol_2 - potential_kjmol) < 20:
                        # Choose a random integer between 0 and 1
                        random_integer = random.randint(0, 1)
                        if random_integer == 1 and self.current_state == self.roots_to_follow[root_1]:
                            self.current_state = self.roots_to_follow[root_2]
                        elif random_integer == 1 and self.current_state == self.roots_to_follow[root_2]:
                            self.current_state = self.roots_to_follow[root_1]
                        break
        
        else:
            self.current_state = self.roots_to_follow[0]

        if len(self.roots_to_follow) > 1 and self.step == self.excitation_pulse[0]:
            self.current_state = self.excitation_pulse[1]
            

        potential_kjmol = self.impes_drivers[self.current_state].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
        self.current_gradient = self.impes_drivers[self.current_state].impes_coordinate.gradient

        self.current_energy = potential_kjmol
                        
            # Potential energy is in Hartree, convert to kJ/mol
            

        return self.current_gradient, potential_kjmol


    def update_gradient(self, new_positions):
        """
        Updates and returns the gradient of the QM region.

        :param new_positions:
            The new positions of the atoms in the QM region.
        :return:
            The gradient of the QM region.
        """
        gradient, _ = self.update_gradient_and_energy(new_positions)

        return gradient

    def update_forces(self, context):
        """
        Updates the forces in the system based on a new gradient.

        Args:
            context: The OpenMM context object.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        conversion_factor = (4.184 * hartree_in_kcalpermol() * 10.0 / bohr_in_angstrom()) * unit.kilojoule_per_mole / unit.nanometer
        new_positions = context.getState(getPositions=True).getPositions()
        
        # Update the forces of the QM region
        qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])

        self.velocities_np.append(context.getState(getVelocities=True).getVelocities(True))
        gradient = self.update_gradient(qm_positions)

        positions_ang = (qm_positions) * 10

        atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
        qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]

        new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")

        force = -np.array(gradient) * conversion_factor

        self.all_gradients.append(gradient)
        ############################################################
        #################### Correlation Check #####################
        ############################################################
        allowed = True 
        self.add_a_point = False

        if self.density_around_data_point[1] is not None and self.current_state == self.density_around_data_point[2]:
            current_dihedral = new_molecule.get_dihedral_in_degrees([self.density_around_data_point[1][0], self.density_around_data_point[1][1], self.density_around_data_point[1][2], self.density_around_data_point[1][3]])%360.0
            lower, upper = self.allowed_molecule_deviation[self.current_state][self.density_around_data_point[1]][self.density_around_data_point[3]]

            # Case 1: If boundaries do not wrap (e.g., [-60, 60])
            if lower < upper:
                allowed = lower <= current_dihedral <= upper
            else:
                # Case 2: If boundaries wrap around (e.g., [120, -120])
                allowed = current_dihedral >= lower or current_dihedral <= upper

        if allowed:

            openmm_coordinate = context.getState(getPositions=True).getPositions()
            self.coordinates.append(openmm_coordinate)
            self.velocities.append(context.getState(getVelocities=True).getVelocities())
            self.gradients.append(gradient)
            if self.skipping_value == 0 and self.step + 1 > self.collect_qm_points:

                if self.add_bayes_model and all(len(self.bayes_models[self.current_state][key].observation_queue) > 3 for key in self.impes_drivers[self.current_state].weights):

                    # current_basis = MolecularBasis.read(new_molecule, self.basis_set_label)

                    # qm_energy, scf_tensors = self.compute_energy(self.drivers['ground_state'][0], new_molecule, current_basis)

                    # energy_difference = (abs(qm_energy[0] - self.impes_drivers[self.current_state].impes_coordinate.energy))

                    current_int_coordinates = self.impes_drivers[0].impes_coordinate.internal_coordinates_values.copy()
                    r_correction = 0.0
                    variance = 0.0

                    for key, weight in self.impes_drivers[self.current_state].weights.items():
                        pred, var = self.bayes_models[self.current_state][key].bayes_predict(current_int_coordinates)
                        r_correction += weight * pred
                        variance += weight * var
                    # print('Len of obs', len(self.bayes_models[self.current_state][0].observation_queue))
                    print('Bayes Error', r_correction, 'variance', 2.0 * np.sqrt(variance))
                    if abs(r_correction) > self.energy_threshold * 0.5 or 2.0 * np.sqrt(variance) > 0.4:
                        self.add_a_point = True

                else:

                    print('rmsd bond', np.mean(np.array(self.impes_drivers[0].bond_rmsd)), '\n')
                    print('rmsd angle', np.mean(np.array(self.impes_drivers[0].angle_rmsd)), '\n')
                    print('rmsd dihedral', np.mean(np.array(self.impes_drivers[0].dihedral_rmsd)), '\n')

                    mean_bond_rmsd = np.mean(np.array(self.impes_drivers[0].bond_rmsd))
                    mean_angle_rmsd = np.mean(np.array(self.impes_drivers[0].angle_rmsd))
                    mean_dihedral_rmsd = np.mean(np.array(self.impes_drivers[0].dihedral_rmsd))

                    if mean_bond_rmsd > 0.1 or mean_angle_rmsd > 0.1 or mean_dihedral_rmsd > 1.0:

                        K = 5  # Number of closest previous matches to cache
                        threshold = self.distance_thrsh - 0.05
                        new_coords = new_molecule.get_coordinates_in_bohr()
                        n_atoms = len(new_molecule.get_labels())
                        sym_group = self.non_core_symmetry_groups['gs' if self.current_state == 0 else 'es'][4]
                        molecule_list = self.allowed_molecules[self.current_state]['molecules']

                        if not hasattr(self, "previous_candidate_indices"):
                            self.previous_candidate_indices = []

                        scanned = False
                        closest_indices = []
                        index_added = []

                        # --- 1. Try cached candidates first
                        for idx in self.previous_candidate_indices:
                            if idx >= len(molecule_list):
                                continue
                            checked_molecule = molecule_list[idx]
                            checked_coords = checked_molecule.get_coordinates_in_bohr()
                            checked_distance = self.cartesian_just_distance(checked_coords, new_coords, sym_group)
                            normed_dist = (np.linalg.norm(checked_distance) / np.sqrt(n_atoms)) * bohr_in_angstrom()

                            if idx not in index_added:
                                print("cached candidate idx", idx, "→", normed_dist)
                                closest_indices.append((idx, normed_dist))
                                index_added.append(idx)
                            if normed_dist <= threshold:
                                scanned = True
                                break

                        # --- 2. Fallback: full scan only if no match
                        if not scanned:
                            for idx, checked_molecule in enumerate(molecule_list):
                                checked_coords = checked_molecule.get_coordinates_in_bohr()
                                checked_distance = self.cartesian_just_distance(checked_coords, new_coords, sym_group)
                                normed_dist = (np.linalg.norm(checked_distance) / np.sqrt(n_atoms)) * bohr_in_angstrom()
                                closest_indices.append((idx, normed_dist))
                                if normed_dist <= threshold:
                                    scanned = True
                                    break

                        # --- 3. Update cached candidate indices for the next step
                        closest_indices.sort(key=lambda x: x[1])  # Sort by distance
                        self.previous_candidate_indices = [i for i, _ in closest_indices[:K]]
                        if not scanned:
                            for i, qm_data_point in enumerate(self.qm_data_point_dict[self.current_state], start=1):

                                length_vectors = (self.impes_drivers[self.current_state].impes_coordinate.cartesian_distance_vector(qm_data_point))

                                if (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels()))) * bohr_in_angstrom() <= self.distance_thrsh - 0.05:
                                    self.add_a_point = False
                                    break

                                if i == len(self.qm_data_point_dict[self.current_state]):
                                    self.add_a_point = True
                    # else:
                    #     self.add_a_point = True

            else:
                self.skipping_value -= 1
                if self.skipping_value < 0:
                    self.skipping_value = 0
            
            # if self.step % 100 == 0:
                # self.add_a_point = True
            self.point_checker += 1 
            if self.add_a_point == True:
                self.point_correlation_check(new_molecule)
            if self.point_checker == 0:            
                
                self.point_checker += 1

                context.setPositions(self.coordinates[0])
                self.coordinates = [self.coordinates[0]]
                # context.setVelocities(self.velocities[0])


            # Set initial velocities if the ensemble is NVT or NPT
                if self.ensemble in ['NVT', 'NPT']:

                    context.setVelocitiesToTemperature(self.temperature)
                else:
                    context.setVelocities(self.start_velocities)
                self.velocities = [context.getState(getVelocities=True).getVelocities()]
                
                new_positions = context.getState(getPositions=True).getPositions()
                qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])
                positions_ang = (qm_positions) * 10
                atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
                qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
                new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
                gradient_2 = self.update_gradient(qm_positions)
                force = -np.array(gradient_2) * conversion_factor
            
            if (self.point_checker + self.last_point_added) % self.duration == 0 and self.point_checker != 0:
                print('start when last point was added', self.point_checker + self.last_point_added)
                self.last_point_added = 0
                self.unadded_cycles -= 1

            if self.point_checker < 500 and self.cycle_iteration != self.unadded_cycles:
                self.unadded_cycles += 1
            
            self.coordinates_xyz.append(qm_positions * 10)

            ############################################################
            self.state_specific_molecules[self.current_state].append(new_molecule)
            for i, atom_idx in enumerate(self.qm_atoms):
          
                self.system.getForce(self.qm_force_index).setParticleParameters(i, atom_idx, force[i])
            self.system.getForce(self.qm_force_index).updateParametersInContext(context)


            for root in self.roots_to_follow:
                self.gloabal_sim_informations[f'state_{root}']['pot_energies'].append(self.impes_drivers[root].get_energy() * hartree_in_kcalpermol())
                self.gloabal_sim_informations[f'state_{root}']['gradients'].append(self.impes_drivers[root].get_gradient() * hartree_in_kcalpermol() * (1.0 / bohr_in_angstrom()))
                
        
        else:
            context.setPositions(self.coordinates[0])
            self.coordinates = [self.coordinates[0]]
            if self.ensemble in ['NVT', 'NPT']:

                context.setVelocitiesToTemperature(self.temperature)
            else:
                context.setVelocities(self.start_velocities)
            self.velocities = [context.getState(getVelocities=True).getVelocities()]
            new_positions = context.getState(getPositions=True).getPositions()
            qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])
            positions_ang = (qm_positions) * 10
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            gradient_2 = self.update_gradient(qm_positions)
            force = -np.array(gradient_2) * conversion_factor
            
            self.state_specific_molecules[self.current_state].append(new_molecule)
            self.coordinates_xyz.append(qm_positions * 10)
            for i, atom_idx in enumerate(self.qm_atoms):
                
                print(atom_idx, force[i])
                self.system.getForce(self.qm_force_index).setParticleParameters(i, atom_idx, force[i])
            self.system.getForce(self.qm_force_index).updateParametersInContext(context)
    
    
    ####################################################################
    ################ Functions to expand the database ##################
    ####################################################################

    def point_correlation_check(self, molecule):
        """ Takes the current point on the PES and checks with a QM-energy
            calculation is necessary based on the current difference to the
            interpolation. Based on the difference the step_size of the next
            step is determined.

            :param molecule:
                The molecule corresponding to the current conformation
                in the IM dynamics.
        """
        
        for state, key in enumerate(self.drivers.keys()):
            
            drivers = None
            
            if state == 0 and self.drivers[key] is not None and 0 in self.roots_to_follow and self.current_state == state:
                drivers = self.drivers[key]
   
            elif state == 1 and self.drivers['excited_state'] is not None and any(x > 0 for x in self.roots_to_follow) and self.current_state >= 1:
                drivers = self.drivers[key]
            else:
                continue

            qm_energy = 0
            print('############# Energy is QM claculated ############')
            current_basis = MolecularBasis.read(molecule, self.basis_set_label)
            qm_energy, scf_tensors = self.compute_energy(drivers[0], molecule, current_basis)
            energy_difference = (abs(qm_energy[0] - self.impes_drivers[self.current_state].impes_coordinate.energy))
            
            print('Energy difference', energy_difference * hartree_in_kcalpermol(), 'kcal/mol')
            if self.add_bayes_model:
                gradients = self.compute_gradient(drivers[1], molecule, current_basis, scf_tensors)
                gradient_difference = gradients[0] - self.impes_drivers[self.current_state].impes_coordinate.gradient

                masses = molecule.get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
                grad_mw = gradient_difference / np.sqrt(masses)[:, None]


                impes_coordinate = InterpolationDatapoint(self.z_matrix)
                impes_coordinate.cartesian_coordinates = molecule.get_coordinates_in_bohr()
                impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                impes_coordinate.gradient = grad_mw
                impes_coordinate.transform_gradient()

                current_int_coordinates = self.impes_drivers[self.current_state].impes_coordinate.internal_coordinates_values.copy()

                r_correction = sum(weight * self.bayes_models[self.current_state][key].bayes_predict(current_int_coordinates)[0]
                    for key, weight in self.impes_drivers[self.current_state].weights.items())
                variance = sum(weight * self.bayes_models[self.current_state][key].bayes_predict(current_int_coordinates)[1]
                    for key, weight in self.impes_drivers[self.current_state].weights.items())
                print('correction', r_correction, 'energy differences', energy_difference * hartree_in_kcalpermol(), 'difference', abs(r_correction - energy_difference * hartree_in_kcalpermol()), 'variance', np.sqrt(variance))

                if abs(r_correction - energy_difference * hartree_in_kcalpermol()) > 0.2 or r_correction == 0.0:
                    print('point_would be added')

                    for key, weight in self.impes_drivers[self.current_state].weights.items():

                        # current_int_coordinates[self.non_core_symmetry_groups[5]] = qm_data_point.internal_coordinates_values[self.non_core_symmetry_groups[5]]
                        self.bayes_models[self.current_state][key].bayes_update_sliding(current_int_coordinates, abs(qm_energy[0] - self.impes_drivers[self.current_state].impes_coordinate.energy) * hartree_in_kcalpermol(), impes_coordinate.internal_gradient * hartree_in_kcalpermol(), weight)
                        # mu_r, var_r = self.bayes_models[self.current_state].bayes_predict(current_int_coordinates)

                        # print("Original dE:", energy_difference * hartree_in_kcalpermol())
                        # mu, var = self.bayes_models.bayes_predict(current_int_coordinates)
                        # print("Predicted mu_r:", mu)
                        # print("Residual:", abs(mu - energy_difference * hartree_in_kcalpermol()))
                        # print("Relative error %:", 100 * abs(mu - energy_difference * hartree_in_kcalpermol()) / abs(energy_difference * hartree_in_kcalpermol()))

            # calcualte energy gradient
            self.previous_energy_list.append(energy_difference)

            if len(self.previous_energy_list) > 3:
                # Compute gradient (rate of change of energy difference)
                grad1 = self.previous_energy_list[-2] - self.previous_energy_list[-3]
                grad2 = self.previous_energy_list[-1] - self.previous_energy_list[-2]

                # Base skipping value calculation
                base_skip = min(round(abs(self.energy_threshold / (energy_difference * hartree_in_kcalpermol())**2)), 20) - 1

                # Adjust skipping value based on gradient
                if grad2 > grad1:  # Energy difference is increasing
                    self.skipping_value = max(1, base_skip - 1)  # Reduce skipping for more frequent checks
                else:  # Energy difference is decreasing
                    self.skipping_value = base_skip + 10  # Increase skipping to check less often

                print(f"Energy Difference: {(energy_difference * hartree_in_kcalpermol()):.6f}, Gradient: {grad2 - grad1:.6f}, Skipping Value: {self.skipping_value}")

            else:
                self.skipping_value = min(round(abs(self.energy_threshold / (energy_difference * hartree_in_kcalpermol())**2)), 20)

            print('len of molecules', len(self.sampled_molecules[self.current_state]['molecules']), self.use_opt_confidence_radius)
            if energy_difference * hartree_in_kcalpermol() > self.energy_threshold:
                
                
                
                if self.use_opt_confidence_radius[0] and len(self.sampled_molecules[self.current_state]['molecules']) > 20:
                    self.add_a_point = True
                    if not self.expansion:
                        length_vectors = (self.impes_drivers[-1].impes_coordinate.cartesian_distance_vector(self.qm_data_points[0]))
                        rmsd = (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels())) * bohr_in_angstrom())
                        self.expansion_molecules.append((molecule, energy_difference * hartree_in_kcalpermol(), rmsd, (np.linalg.norm(length_vectors))))
                        self.last_point_added = self.point_checker - 1
                        self.point_checker = 0
                
                elif energy_difference * hartree_in_kcalpermol() > self.energy_threshold and len(self.sampled_molecules[self.current_state]['molecules']) < 100 and self.use_opt_confidence_radius[0]:
  
                    self.sampled_molecules[self.current_state]['molecules'].append(molecule)
                    self.sampled_molecules[self.current_state]['im_energies'].append(self.impes_drivers[self.current_state].impes_coordinate.energy)
                    self.sampled_molecules[self.current_state]['qm_energies'].append(qm_energy[0])
                    self.last_point_added = self.point_checker - 1
                    self.point_checker = 0
                    self.add_a_point = False
                
                else:
                    self.add_a_point = True
                    if not self.expansion:
                        length_vectors = (self.impes_drivers[-1].impes_coordinate.cartesian_distance_vector(self.qm_data_points[0]))
                        rmsd = (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels())) * bohr_in_angstrom())
                        self.expansion_molecules.append((molecule, energy_difference * hartree_in_kcalpermol(), rmsd, (np.linalg.norm(length_vectors))))
                        self.last_point_added = self.point_checker - 1
                        self.point_checker = 0

            else:
                self.allowed_molecules[self.current_state]['molecules'].append(molecule)
                self.allowed_molecules[self.current_state]['im_energies'].append(self.impes_drivers[self.current_state].impes_coordinate.energy)
                self.allowed_molecules[self.current_state]['qm_energies'].append(qm_energy[0])
                self.sampled_molecules[self.current_state]['molecules'].append(molecule)
                self.sampled_molecules[self.current_state]['im_energies'].append(self.impes_drivers[self.current_state].impes_coordinate.energy)
                self.sampled_molecules[self.current_state]['qm_energies'].append(qm_energy[0])
                
                self.add_a_point = False
            if self.add_a_point and self.expansion:
                print('✨ A point is added! ✨', self.point_checker)
                print(molecule.get_xyz_string())

                ############# Implement constraint optimization ############

                if self.identfy_relevant_int_coordinates:
                    opt_qm_driver = ScfRestrictedDriver()
                    opt_qm_driver.xcfun = 'b3lyp'
                    _, scf_tensors = self.compute_energy(opt_qm_driver, molecule, current_basis)
                    opt_drv = OptimizationDriver(opt_qm_driver)
                    opt_drv.ostream.mute()
                    
                    # interpolation_driver = InterpolationDriver(self.z_matrix)
                    # interpolation_driver.update_settings(self.interpolation_settings[self.current_state])
                    # if self.current_state == 0:
                    #     interpolation_driver.symmetry_information = self.non_core_symmetry_groups['gs']
                    # else:
                    #     interpolation_driver.symmetry_information = self.non_core_symmetry_groups['es']

                    # interpolation_driver.distance_thrsh = 1000
                    # interpolation_driver.exponent_p = 2
                    # interpolation_driver.store_weights = True
                    # interpolation_driver.qm_
                    # interpolation_driver.compute(molecule)
                    self.impes_drivers[self.current_state].compute(molecule)
                    current_weights = self.impes_drivers[self.current_state].weights

                    weights = [value for _, value in current_weights.items()]
                    used_labels = [label_idx for label_idx, _ in current_weights.items()]
                    internal_coordinate_datapoints = []
                    for i, weight in enumerate(weights):
                        
                        if weight >= max(weights) - 0.2:
                            
                            internal_coordinate_datapoints.append(self.qm_data_point_dict[self.current_state][used_labels[i]])
                        else:
                            break

                    # qm_datapoints_weighted = [qm_datapoint for qm_datapoint in enumerate if ]
                    constraints = self.impes_drivers[self.current_state].determine_important_internal_coordinates(qm_energy, molecule, self.z_matrix, internal_coordinate_datapoints)
                    
                    
                    print('FINAL WIEGHTS', current_weights)
                    print('CONSTRAINTS', constraints)

                    opt_constraint_list = []
                    for constraint in constraints:
                        if len(constraint) == 2:
                            opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                            opt_constraint_list.append(opt_constraint)
                        
                        elif len(constraint) == 3:
                            opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                            opt_constraint_list.append(opt_constraint)
                    
                        else:
                            opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                            opt_constraint_list.append(opt_constraint)
                    opt_drv.constraints = opt_constraint_list
                    opt_results = opt_drv.compute(molecule, current_basis, scf_tensors)
                    optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                    opt_current_basis = MolecularBasis.read(optimized_molecule, current_basis.get_main_basis_label())
                    qm_energy, scf_tensors = self.compute_energy(drivers[0], optimized_molecule, opt_current_basis)
                    print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                    
                    label = f"point_{len(self.sorted_state_spec_im_labels[self.current_state]) + 1}"
                    self.add_point(optimized_molecule, label, current_basis, self.non_core_symmetry_groups)
                    self.last_point_added = self.point_checker - 1
                    self.point_checker = 0
                    self.point_adding_molecule[self.step] = (molecule, qm_energy, label)
                
                else:
                    label = f"point_{len(self.sorted_state_spec_im_labels[self.current_state]) + 1}"
                    self.add_point(molecule, label, current_basis, self.non_core_symmetry_groups)
                    self.last_point_added = self.point_checker - 1
                    self.point_checker = 0
                    self.point_adding_molecule[self.step] = (molecule, qm_energy, label)


    def add_point(self, molecule, label, basis, symmetry_information):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param label:
                the label for the new point to be added.
            :param energy:
                the energy of the previous QM calcualtion.
            :param basis:
                the basis set (if required).
            :scf_results:
                the scf_results of previous QM calculation (if required).
        """

        adjusted_molecule = {'gs': [], 'es': []}
        symmetry_mapping_groups = []
        symmetry_exclusion_groups = []
        new_label = label
        category_label = None
 
        if len(symmetry_information) != 0:

            for key, sym_inf in symmetry_information.items():
                if 'gs' == key and 0 in self.roots_to_follow and len(sym_inf[2]) != 0:
                    symmetry_mapping_groups = [item for item in range(len(molecule.get_labels()))]
                    symmetry_exclusion_groups = [item for element in sym_inf[1] for item in element]
                    sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(sym_inf[1], sym_inf[5])
                    
                    # Generate all combinations
                    keys = list(sym_dihedrals.keys())
                    values = list(sym_dihedrals.values())
                    combinations = list(itertools.product(*values))

                    # Convert to list of dictionaries
                    molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)
                    print('MOlecule configs', molecule_configs)
                    for i, molecule_config in enumerate(molecule_configs):
                        cur_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
                        dihedral_to_change = []
                        for dihedral, angle in molecule_config.items():
                            cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], angle, 'radian')
                            dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])

                        current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())
                        adjusted_molecule['gs'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change))
                elif 'es' == key and any(x > 0 for x in self.roots_to_follow) and len(sym_inf[2]) != 0:
                    symmetry_mapping_groups = [item for item in range(len(molecule.get_labels()))]
                    symmetry_exclusion_groups = [item for element in sym_inf[1] for item in element]
                    sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(sym_inf[1], sym_inf[5])
                    
                    # Generate all combinations
                    keys = list(sym_dihedrals.keys())
                    values = list(sym_dihedrals.values())
                    combinations = list(itertools.product(*values))

                    # Convert to list of dictionaries
                    molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)

                    for i, molecule_config in enumerate(molecule_configs):
                        cur_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
                        dihedral_to_change = []
                        for dihedral, angle in molecule_config.items():
                            cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], angle, 'radian')
                            dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])

                        current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())
                        adjusted_molecule['es'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change))
                else:
                    adjusted_molecule[key].append((molecule, basis, 1, None))
        
        else:
            adjusted_molecule['gs'].append((molecule, basis, 1, None))  
            adjusted_molecule['es'].append((molecule, basis, 1, None))
        
        for state, key in enumerate(self.drivers.keys()):
            
            drivers = None
            
            if state == 0 and self.drivers[key] is not None and 0 in self.roots_to_follow and self.current_state == 0:
                drivers = self.drivers[key]
   
            elif state == 1 and self.drivers['excited_state'] is not None and any(x > 0 for x in self.roots_to_follow) and self.current_state > 0:
                drivers = self.drivers[key]
            else:
                continue
            
            for mol_state, entries in adjusted_molecule.items():
                
                if state == 0 and mol_state != 'gs':
                    continue
                if state == 1 and mol_state != 'es':
                    continue
                
                new_label = label
                for label_counter, mol_basis in enumerate(entries):

                    energies, scf_results = self.compute_energy(drivers[0], mol_basis[0], mol_basis[1])

                    gradients = self.compute_gradient(drivers[1], mol_basis[0], mol_basis[1], scf_results)
                    hessians = self.compute_hessian(drivers[2], mol_basis[0], mol_basis[1])

                    masses = molecule.get_masses().copy()
                    masses_cart = np.repeat(masses, 3)
                    inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)

                    for number in range(len(energies)):
                        
                        if number != self.current_state:
                            continue

                            
                        grad = gradients[number].copy()
                        hess = hessians[number].copy()


                        if self.ghost_atom[0]:

                            # All DOFs: 0 to 17
                            n_atoms = len(molecule.get_labels()) + 1
                            full_cart_indices = np.arange(3 * n_atoms)

                            # Cartesian indices to remove: 3 * ghost_atom_index + [0, 1, 2]
                            ghost_cartesian_indices = np.arange(3 * self.ghost_atom[1], 3 * self.ghost_atom[1] + 3)

                            # Remaining indices
                            indices = np.setdiff1d(full_cart_indices, ghost_cartesian_indices)

                            grad = grad.reshape(-1)  # if it's (6, 3), flatten to (18,)
                            grad = grad[indices].reshape(-1, 3) 
                            hess = hess[np.ix_(indices, indices)]  # shape will now be (15, 15)

                           

                        grad_mw = grad / np.sqrt(masses)[:, None]
                        H_mw = hess * inv_sqrt_masses[:, np.newaxis] * inv_sqrt_masses[np.newaxis, :]

                        
                        impes_coordinate = InterpolationDatapoint(self.z_matrix)
                        impes_coordinate.cartesian_coordinates = mol_basis[0].get_coordinates_in_bohr()
                        impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                        impes_coordinate.energy = energies[number]
                        impes_coordinate.gradient = grad_mw
                        impes_coordinate.hessian = H_mw
                        impes_coordinate.transform_gradient_and_hessian()
                        # trust_radius = impes_coordinate.determine_trust_radius(molecule, self.interpolation_settings, self.dynamics_settings, label, specific_coordiante, self.allowed_deviation, symmetry_groups=self.symmetry_groups)
                        sampled_space = False
                        sampled_distances = []
                        impes_coordinate.confidence_radius = float(self.use_opt_confidence_radius[2])

                        if label_counter > 0:
                            new_label = f'{label}_symmetry_{label_counter}'

                        print('Org Molecules', mol_basis[0].get_xyz_string())

                        if mol_basis[2] > 1:
                                
                            rotation_combinations = None
                            
                            if mol_basis[2] == 3:
                                from itertools import product
                                rotations = [0.0, 2*np.pi/3, 4*np.pi/3]
                                dihedrals = mol_basis[3]  # your list of rotatable dihedrals
                                rotation_combinations = list(product(rotations, repeat=len(dihedrals)))
                            
                            if mol_basis[2] == 2:
                                from itertools import product
                                rotations = [0.0, np.pi]
                                dihedrals = mol_basis[3]  # your list of rotatable dihedrals
                                rotation_combinations = list(product(rotations, repeat=len(dihedrals)))
                            print('rot comb', rotation_combinations)
                            
                            org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                            masks = [org_mask]
                            for combo in rotation_combinations:
                                if all(r == 0 for r in combo):
                                    continue  # skip the base geometry if desired
                                
                        
                                rot_mol = Molecule(self.molecule.get_labels(), mol_basis[0].get_coordinates_in_bohr(), 'bohr')
                                for angle, dihedral in zip(combo, dihedrals):
                                    
                                    rot_mol.set_dihedral(dihedral, mol_basis[0].get_dihedral(dihedral, 'radian') - angle, 'radian')
                                print('rot_molecules', combo, rot_mol.get_xyz_string())
                                target_coordinates = self.calculate_translation_coordinates(rot_mol.get_coordinates_in_bohr())
                                reference_coordinates = self.calculate_translation_coordinates(mol_basis[0].get_coordinates_in_bohr())
                                
                                target_coordinates_core = target_coordinates.copy()
                                reference_coordinates_core = reference_coordinates.copy()
                                    
                                target_coordinates_core = np.delete(target_coordinates, symmetry_exclusion_groups, axis=0)
                                reference_coordinates_core = np.delete(reference_coordinates, symmetry_exclusion_groups, axis=0)
                                rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                                        reference_coordinates_core)
                                # Rotate the data point
                                rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
                                rotated_coordinates = np.ascontiguousarray(rotated_coordinates)
        
                                mapping_dict_12 = self.perform_symmetry_assignment(symmetry_mapping_groups, symmetry_exclusion_groups, 
                                                                                    mol_basis[0].get_coordinates_in_bohr()[symmetry_exclusion_groups], rotated_coordinates[symmetry_exclusion_groups])
                                z_matrix_dict = {tuple(sorted(element)): i 
                                    for i, element in enumerate(impes_coordinate.z_matrix)}
                                
                                mapping_dict = [mapping_dict_12]
                                
                                
                                reorded_int_coords = [impes_coordinate.internal_coordinates_values.copy()]
            
                                for ord in range(len(mapping_dict)):
                                    mask = []
                                    reorded_int_coord = np.zeros_like(impes_coordinate.internal_coordinates_values)
                                    for i, element in enumerate(impes_coordinate.z_matrix):
                                        # Otherwise, reorder the element
                                        reordered_element = [mapping_dict[ord].get(x, x) for x in element]
                                        key = tuple(sorted(reordered_element))
                                        z_mat_index = z_matrix_dict.get(key)
                                        mask.append(z_mat_index)
                                        reorded_int_coord[i] = (float(reorded_int_coords[0][z_mat_index]))
                                    reorded_int_coords.append(reorded_int_coord)
                                    masks.append(mask)

                            impes_coordinate.mapping_masks = masks
                        else:
                            org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                            masks = [org_mask]
                            impes_coordinate.mapping_masks = masks
                    
                    # while not sampled_space and self.use_opt_confidence_radius[0]:
                    #     # selected_molecules = random.sample(self.state_specific_molecules[state + i], 30 - len(self.allowed_molecules[state + i]))
                        
                    #     for mol in self.state_specific_molecules[state + i]:
                    #         _, distance, _ = self.calculate_distance_to_ref(mol.get_coordinates_in_bohr(), molecule.get_coordinates_in_bohr())
                    #         sampled_distances.append(distance)

                    #     sampled_distances_np = np.array(sampled_distances)
                    #     mu = np.mean(sampled_distances_np)
                    #     sigma = np.std(sampled_distances_np)

                    #     weights = np.exp(-((sampled_distances_np - mu) ** 2) / (2 * sigma ** 2))
                    #     # kT = 4.0  # You can tune this value
                    #     # weights = np.exp(-sampled_distances_np / kT)
                        
                    #     probabilities = weights / sum(weights)
                        
                    #     num_samples = 20
                    #     sampled_indices = np.random.choice(len(sampled_distances), size=num_samples, replace=False, p=probabilities)
                    #     print(sampled_indices)
                    #     chosen_structures = np.array(self.state_specific_molecules[state + i])[sampled_indices]
                        
                    #     for sel_mol in chosen_structures:
                    #         energies, scf_results = self.compute_energy(drivers[0], sel_mol, basis)
                    #         self.impes_drivers[state + i].qm_data_points = self.qm_data_point_dict[state + i]
                    #         self.impes_drivers[state + i].compute(sel_mol)
                    #         if abs(energies[i] - self.impes_drivers[state + i].get_energy()) * hartree_in_kcalpermol() < self.energy_threshold:
                    #             self.sampled_molecules[state + i]['molecules'].append(sel_mol)
                    #             self.sampled_molecules[state + i]['im_energies'].append(self.impes_drivers[state + i].get_energy())
                    #             self.sampled_molecules[state + i]['qm_energies'].append(energies[i])
                    #     sampled_space = True
                    
                    
                        if label_counter == 0:
                            category_label = label
                            if impes_coordinate.use_inverse_bond_length:
                                category_label += "_rinv"
                            else:
                                category_label += "_r"
                            if impes_coordinate.use_cosine_dihedral:
                                category_label += "_cosine"
                            else:
                                category_label += "_dihedral"
                            print('dict and key', self.impes_drivers[state + number].qm_symmetry_data_points, category_label)
                            self.qm_data_point_dict[state + number].append(impes_coordinate)
                            self.sorted_state_spec_im_labels[state + number].append(new_label)
                            self.qm_symmetry_datapoint_dict[state + number][category_label] = [impes_coordinate]
                            self.impes_drivers[state + number].qm_symmetry_data_points[category_label] = [impes_coordinate]
                            impes_coordinate.point_label = category_label
                            self.qm_energies_dict[state + number].append(energies[number])
                            
                            self.impes_drivers[state + number].qm_data_points = self.qm_data_point_dict[self.current_state]
                            self.impes_drivers[state + number].labels = self.sorted_state_spec_im_labels[self.current_state]
                            print('I am writing the file', label, self.interpolation_settings[state + number]['imforcefield_file'])
                        
                            self.density_around_data_point[0][state + number] += 1
                            if self.add_bayes_model:
                                bayes_model = LocalBayesResidual(self.z_matrix)
                                bayes_model.compute_internal_coordinates_values(impes_coordinate.cartesian_coordinates)
                                if state + number == 0:
                                    bayes_model.symmetry_information = self.non_core_symmetry_groups['gs']
                                else:
                                    bayes_model.symmetry_information = self.non_core_symmetry_groups['es']
                                bayes_model.init_bayesian()
                                self.bayes_models[state + number].append(bayes_model)
                        else:
                            
                            self.qm_symmetry_datapoint_dict[state + number][category_label].append(impes_coordinate)
                            self.impes_drivers[state + number].qm_symmetry_data_points[category_label].append(impes_coordinate)
                        
                        impes_coordinate.confidence_radius = self.use_opt_confidence_radius[2]
                        # impes_coordinate.write_hdf5(self.interpolation_settings[state + number]['imforcefield_file'], new_label)
                        # trust_radius = self.determine_trust_radius_single(self.sampled_molecules['molecules'], self.sampled_molecules['qm_energies'], self.sampled_molecules['im_energies'], self.qm_data_point_dict[state + i][:-1], impes_coordinate, self.interpolation_settings[state + i])
                        
                        print('data list', len(self.sampled_molecules[state + number]['molecules']))
                        # if len(self.qm_data_point_dict[state + i]) > 2:
                        #     remove_datapoints = False
                        #     cut_off = 0
                        #     for counter in range(3, len(self.qm_data_point_dict[state + i])):
            
                        #         if not remove_datapoints:
                        #             trust_radius_multi = self.determine_trust_radius(self.sampled_molecules['molecules'], self.sampled_molecules['qm_energies'], self.sampled_molecules['im_energies'], self.qm_data_point_dict[state + i][:counter], self.interpolation_settings[state + i])
                                    
                        #             for idx, trust_radius in enumerate(trust_radius_multi):
                        #                 self.qm_data_point_dict[state + i][idx].confidence_radius = trust_radius
                        #             interpolation_driver = InterpolationDriver(self.z_matrix)
                        #             interpolation_driver.update_settings(self.interpolation_settings[state + i])
                        #             interpolation_driver.symmetry_sub_groups = self.non_core_symmetry_groups
                        #             interpolation_driver.distance_thrsh = 1000
                        #             interpolation_driver.exponent_p = 2
                        #             interpolation_driver.print = False
                        #             interpolation_driver.qm_data_points = self.qm_data_point_dict[state + i][:len(trust_radius_multi) - 1]
                        #             curr_dp_molecule = Molecule(molecule.get_labels(), self.qm_data_point_dict[state + i][len(trust_radius_multi)].cartesian_coordinates, 'bohr')
                        #             interpolation_driver.compute(curr_dp_molecule)
                        #             new_im_energy = interpolation_driver.get_energy()
                        #             diff = abs(new_im_energy - self.qm_data_point_dict[state + i][len(trust_radius_multi)].energy) * hartree_in_kcalpermol()
                        #             print('DIFF for removeale', diff)
                        #             if diff < self.energy_threshold - self.energy_threshold / 3.0:
                        #                 remove_datapoints = True
                                
                        #         if remove_datapoints:
                        #             self.qm_data_point_dict[state + i][counter - 1].remove_point_from_hdf5(self.interpolation_settings[state + i]['imforcefield_file'], self.sorted_state_spec_im_labels[state + i][counter - 1], use_inverse_bond_length=True, use_cosine_dihedral=False)
                        #             cut_off = counter
                        #     if remove_datapoints:
                        #         self.qm_data_point_dict[state + i] = self.qm_data_point_dict[state + i][:cut_off - 1]
                        #         exit()  
                        impes_coordinate.write_hdf5(self.interpolation_settings[state + number]['imforcefield_file'], new_label)
                        if label_counter == 0:
                            if self.use_opt_confidence_radius[0]:
                                
                                trust_radius = None
                                if self.use_opt_confidence_radius[1] == 'single':
                                    sym_dict = self.non_core_symmetry_groups['gs']
                                    if state + number > 0:
                                        sym_dict = self.non_core_symmetry_groups['es']
                                        
                                    
                                    print(len(self.qm_data_point_dict[state + number][:-1]), len(self.qm_data_point_dict[state + number][:]))
                                    trust_radius = self.determine_trust_radius_single(self.sampled_molecules[state + number]['molecules'], 
                                                                                    self.sampled_molecules[state + number]['qm_energies'],
                                                                                    self.sampled_molecules[state + number]['im_energies'], 
                                                                                    self.qm_data_point_dict[state + number][:-1], 
                                                                                    self.qm_data_point_dict[state + number][-1], 
                                                                                    self.interpolation_settings[state + number],
                                                                                    self.qm_symmetry_datapoint_dict[state + number],
                                                                                    sym_dict)
                                elif self.use_opt_confidence_radius[1] == 'multi':
                                    
                                    sym_dict = self.non_core_symmetry_groups['gs']
                                    if state + number > 0:
                                        sym_dict = self.non_core_symmetry_groups['es']
                                    trust_radius = self.determine_trust_radius(self.sampled_molecules[state + number]['molecules'], 
                                                                            self.sampled_molecules[state + number]['qm_energies'], 
                                                                            self.sampled_molecules[state + number]['im_energies'], 
                                                                            self.qm_data_point_dict[state + number], 
                                                                            self.interpolation_settings[state + number],
                                                                            self.qm_symmetry_datapoint_dict[state + number],
                                                                            sym_dict)
                                
                                else:
                                    trust_radius = self.determine_beysian_trust_radius(self.sampled_molecules[state + number]['molecules'], 
                                                                                    self.sampled_molecules[state + number]['qm_energies'], 
                                                                                    self.qm_data_point_dict[state + number], 
                                                                                    self.interpolation_settings[state + number], 
                                                                                    self.qm_symmetry_datapoint_dict[state + number],
                                                                                    sym_dict)
                            
                                # exit()
                                for idx, trust_radius in enumerate(trust_radius):
                                    print(self.sorted_state_spec_im_labels[state + number][idx])
                                    self.qm_data_point_dict[state + number][idx].update_confidence_radius(self.interpolation_settings[state + number]['imforcefield_file'], self.sorted_state_spec_im_labels[state + number][idx], trust_radius)
                                    self.qm_data_point_dict[state + number][idx].confidence_radius = trust_radius
                            
                    

    def determine_beysian_trust_radius(self, molecules, qm_energies, current_datapoints, interpolation_setting, sym_datapoints, sym_dict):
    
        trust_radii = []
        for dp in current_datapoints:   
            sum_sq_error = 0.0
            combined_datapoints = [dp]
            


            for i, mol in enumerate(molecules):
                _, distance, _ = self.calculate_distance_to_ref(mol.get_coordinates_in_bohr(), dp.cartesian_coordinates)

                interpolation_driver = InterpolationDriver(self.z_matrix)
                interpolation_driver.update_settings(interpolation_setting)
                interpolation_driver.symmetry_information = sym_dict
                interpolation_driver.qm_symmetry_data_points = sym_datapoints
                interpolation_driver.distance_thrsh = 1000
                interpolation_driver.exponent_p = 2
                interpolation_driver.print = False
                interpolation_driver.qm_data_points = combined_datapoints
                
                interpolation_driver.compute(mol)
                new_im_energy = interpolation_driver.get_energy()
                diff = (new_im_energy - qm_energies[i]) * hartree_in_kcalpermol()
                sum_sq_error += (diff)**2 / (0.1**2 * distance**6)
                # interpolation_driver.determine_important_internal_coordinates(qm_e[i], mol, self.z_matrix, combined_datapoints)

            bey_trust_radius = (1/sum_sq_error)**(1/6)

            trust_radii.append(bey_trust_radius)
            
        return trust_radii
    
    def determine_beysian_grad_trust_radius(self, molecules, qm_energies, current_datapoint, interpolation_setting, sym_datapoints, sym_dict):
    
                
        sum_sq_error = 0.0
        combined_datapoints = [current_datapoint]
        for i, mol in enumerate(molecules):
            _, distance, _ = self.calculate_distance_to_ref(mol.get_coordinates_in_bohr(), current_datapoint.cartesian_coordinates)

            interpolation_driver = InterpolationDriver(self.z_matrix)
            interpolation_driver.update_settings(interpolation_setting)
            interpolation_driver.symmetry_information = sym_dict
            interpolation_driver.qm_symmetry_data_points = sym_datapoints
            interpolation_driver.distance_thrsh = 1000
            interpolation_driver.exponent_p = 2
            interpolation_driver.print = False
            interpolation_driver.qm_data_points = combined_datapoints
            
            interpolation_driver.compute(mol)
            new_im_energy = interpolation_driver.get_energy()
            diff = (new_im_energy - qm_energies[i]) * hartree_in_kcalpermol()
            sum_sq_error += (diff)**2 / (0.1**2 * distance**6)

        bey_trust_radius = (1/sum_sq_error)**(1/6)

        return bey_trust_radius

    def determine_trust_radius(self, molecules, qm_energies, im_energies, datapoints, interpolation_setting, sym_datapoints, sym_dict):
        
        def obj_energy_function(alphas, structure_list, qm_e, im_e, dps, impes_dict, sym_datapoints, sym_dict):
                
            sum_sq_error = 0.

            for i, dp in enumerate(dps):
                dp.confidence_radius = alphas[i]

            for i, mol in enumerate(structure_list):

                interpolation_driver = InterpolationDriver(self.z_matrix)
                interpolation_driver.update_settings(impes_dict)
                interpolation_driver.symmetry_information = sym_dict
                interpolation_driver.qm_symmetry_data_points = sym_datapoints
                interpolation_driver.distance_thrsh = 1000
                interpolation_driver.exponent_p = 2
                interpolation_driver.print = False
                interpolation_driver.qm_data_points = dps
                
                interpolation_driver.compute(mol)
                new_im_energy = interpolation_driver.get_energy()
                diff = (new_im_energy * hartree_in_kcalpermol() - qm_e[i] * hartree_in_kcalpermol())
                sum_sq_error += (diff)**2

            print('sum_of_sqaure', sum_sq_error, alphas)
            return sum_sq_error
        
        def obj_gradient_function(alphas, structure_list, qm_e, im_e, dps, impes_dict, sym_datapoints, sym_dict):
            
            n_points = len(dps)
            dF_dalphas = np.zeros(n_points, dtype=float)

            for i, dp in enumerate(dps):
                dp.confidence_radius = alphas[i]

            for i, mol in enumerate(structure_list):
                interpolation_driver = InterpolationDriver(self.z_matrix)
                interpolation_driver.update_settings(impes_dict)
                interpolation_driver.symmetry_information = sym_dict
                interpolation_driver.qm_symmetry_data_points = sym_datapoints
                interpolation_driver.distance_thrsh = 1000
                interpolation_driver.exponent_p = 2
                interpolation_driver.qm_data_points = dps
                
                interpolation_driver.compute(mol)
                
                sum_of_weights = interpolation_driver.sum_of_weights
                new_im_energy = interpolation_driver.get_energy()
                
                diff = (new_im_energy * hartree_in_kcalpermol() - qm_e[i] * hartree_in_kcalpermol())
                    
                S     = interpolation_driver.sum_of_weights      # scalar
                E_hat = interpolation_driver.get_energy()        # \hat{E}_m
    
                for j, dp in enumerate(dps):
                    # weight derivative for THIS alpha_j
                    dw = interpolation_driver.trust_radius_weight_gradient(dp)   # w'_{mj}

                    P_j = interpolation_driver.potentials[j]                     # P_j
                    # ---- one-liner from (★) ----
                    dE_dalpha = dw * (P_j * hartree_in_kcalpermol() - new_im_energy * hartree_in_kcalpermol()) / S                           # hartree
                    # -------------------------------------------

                    # objective derivative (★★) – ONE unit conversion only
                    dF_dalphas[j] += 2.0 * diff * dE_dalpha   
                    

            

            return dF_dalphas
        
        def optimize_trust_radius(alphas, geom_list, E_ref_list, E_im_list, dps, impes_dict, sym_datapoints, sym_dict):
            """
            Perform the gradient-based optimization to find R*
            that minimizes the sum of squared errors to reference QM energies.
            """

            # We pass the data as extra arguments to avoid repeated global references:
            args = (geom_list, E_ref_list, E_im_list, dps, impes_dict, sym_datapoints, sym_dict)

            # Minimizer expects a function returning just the objective,
            # plus a function returning the gradient if we pass `jac=True`.
            res = minimize(
                fun = obj_energy_function,
                x0 = alphas,
                args = args,
                method='L-BFGS-B',
                jac = obj_gradient_function,
                bounds=[(0.01, 1.5)],
                options={'disp': True}
            )

            # res.x is the optimal R, res.fun is F(R) at optimum
            return res
        

        inital_alphas = [dp.confidence_radius for dp in datapoints]

        print('INPUT Trust radius', inital_alphas)

        trust_radius = optimize_trust_radius(inital_alphas, molecules, qm_energies, im_energies, datapoints, interpolation_setting, sym_datapoints, sym_dict)
        
        return trust_radius['x']
    

    def determine_trust_radius_single(self, molecules, qm_energies, im_energies, old_datapoints, current_datapoint, interpolation_setting, sym_datapoints, sym_dict):
        
        def obj_energy_function(alpha, structure_list, qm_e, im_e, old_dps, curr_dp, impes_dict, sym_datapoints, sym_dict):
                
            
            # alphas = np.linspace(0.01, 10.0, 100)
            # min_diff = np.inf
            # global_minimum = []
            # for alpha in alphas:
            
            alpha = [alpha]
            sum_sq_error = 0.0
            new_datapoint = curr_dp
            new_datapoint.confidence_radius = float(alpha[0])
            combined_datapoints = [dp for dp in old_dps]
            combined_datapoints.append(new_datapoint)
            
            
            for i, mol in enumerate(structure_list):
                _, distance, _ = self.calculate_distance_to_ref(mol.get_coordinates_in_bohr(), new_datapoint.cartesian_coordinates)
                interpolation_driver = InterpolationDriver(self.z_matrix)
                interpolation_driver.update_settings(impes_dict)
                interpolation_driver.symmetry_information = sym_dict
                interpolation_driver.qm_symmetry_data_points = sym_datapoints
                interpolation_driver.distance_thrsh = 1000
                interpolation_driver.exponent_p = 2
                interpolation_driver.print = False
                interpolation_driver.use_symmetry = False
                interpolation_driver.qm_data_points = combined_datapoints
                
                interpolation_driver.compute(mol)
                new_im_energy = interpolation_driver.get_energy()
                diff = (new_im_energy - qm_e[i]) * hartree_in_kcalpermol()
                sum_sq_error += (diff)**2
                old_diff = abs(im_e[i] - qm_e[i]) * hartree_in_kcalpermol()

                # interpolation_driver.determine_important_internal_coordinates(qm_e[i], mol, self.z_matrix, combined_datapoints)

            #     if sum_sq_error < min_diff:
            #         min_diff = sum_sq_error
            #         global_minimum = [min_diff, alpha]
            
            # print('Global minimum', global_minimum)
            # exit()

            return sum_sq_error
        
        def obj_gradient_function(alpha, structure_list, qm_e, im_e, old_dps, curr_dp, impes_dict, sym_datapoints, sym_dict):
            
            dF_dalpha = 0.0
            new_datapoint = curr_dp
            new_datapoint.confidence_radius  = alpha[0]
            combined_datapoints = [dp for dp in old_dps]
            combined_datapoints.append(new_datapoint)

            for i, mol in enumerate(structure_list):
                interpolation_driver = InterpolationDriver(self.z_matrix)
                interpolation_driver.update_settings(impes_dict)
                interpolation_driver.symmetry_information = sym_dict
                interpolation_driver.qm_symmetry_data_points = sym_datapoints
                interpolation_driver.distance_thrsh = 1000
                interpolation_driver.exponent_p = 2
                interpolation_driver.use_symmetry = False
                interpolation_driver.qm_data_points = combined_datapoints
                
                interpolation_driver.compute(mol)
                
                sum_of_weights = interpolation_driver.sum_of_weights
                gradient_weight = interpolation_driver.trust_radius_weight_gradient(new_datapoint)
                new_im_energy = interpolation_driver.get_energy()

                diff = (new_im_energy * hartree_in_kcalpermol() - qm_e[i] * hartree_in_kcalpermol())                
                
                energy_gradient = (gradient_weight * (interpolation_driver.potentials[-1] * hartree_in_kcalpermol()
                            - new_im_energy * hartree_in_kcalpermol()) / sum_of_weights)
                


                energy_gradient = energy_gradient

                dF_dalpha += 2 * diff * energy_gradient 

            print('dF_alpha', dF_dalpha)
            return np.array([dF_dalpha])
        
        def optimize_trust_radius(alpha, geom_list, E_ref_list, E_im_list, old_dps, curr_dp, impes_dict, sym_datapoints, sym_dict):
            """
            Perform the gradient-based optimization to find R*
            that minimizes the sum of squared errors to reference QM energies.
            """

            # We pass the data as extra arguments to avoid repeated global references:
            args = (geom_list, E_ref_list, E_im_list, old_dps, curr_dp, impes_dict, sym_datapoints, sym_dict)

            # Minimizer expects a function returning just the objective,
            # plus a function returning the gradient if we pass `jac=True`.
            res = minimize(
                fun = obj_energy_function,
                x0 = alpha,
                args = args,
                method='L-BFGS-B',
                jac = obj_gradient_function,
                bounds=[(0.01, 10.0)],
                options={'disp': True}
            )

            # res.x is the optimal R, res.fun is F(R) at optimum
            return res
        
                

        print('INPUT Trust radius', current_datapoint.confidence_radius)
        current_dps_list = old_datapoints.copy()
        trust_radius = optimize_trust_radius(current_datapoint.confidence_radius, molecules, qm_energies, im_energies, current_dps_list, current_datapoint, interpolation_setting, sym_datapoints, sym_dict)
        print('final trust radius', trust_radius)
        return trust_radius['x']

    
    
    def perform_symmetry_assignment(self, atom_map, sym_group, reference_group, datapoint_group):
        """ Performs the atom mapping. """
        from scipy.optimize import linear_sum_assignment
        new_map = np.array(atom_map.copy())
        mapping_dict = {}
        # cost = self.get_dihedral_cost(atom_map, sym_group, non_group_atoms)
        cost = np.linalg.norm(datapoint_group[:, np.newaxis, :] - reference_group[np.newaxis, :, :], axis=2)
        row, col = linear_sum_assignment(cost)
        assigned = False
        if not np.equal(row, col).all():
            assigned = True
            
            # atom_maps = self.linear_assignment_solver(cost)

            reordred_arr = np.array(sym_group)[col]
            new_map[sym_group] = new_map[reordred_arr]

            mapping_dict = {org: new for org, new in zip(np.array(sym_group), reordred_arr)}
        
        return mapping_dict

    def adjust_symmetry_dihedrals(self, symmetry_groups, rot_bonds):
        
        def symmetry_group_dihedral(reference_set, dihedrals, rot_bonds):
            rot_bond_set = {frozenset(bond) for bond in rot_bonds}

            filtered_dihedrals = []
            for d in dihedrals:
                middle_bond = frozenset([d[1], d[2]])

                if middle_bond in rot_bond_set:
                    common_elements = [x for x in [d[0], d[3]] if x in reference_set]
                    if len(common_elements) == 1:
                        filtered_dihedrals.append(d)
            return filtered_dihedrals
        
        all_dihedrals = [element for element in self.z_matrix if len(element) == 4]

        symmetry_group_dihedral_dict = {} 
        angles_to_set = {}
        periodicities = {}
        dihedral_groups = {2: [], 3: []}

        for symmetry_group in symmetry_groups:

            symmetry_group_dihedral_list = symmetry_group_dihedral(symmetry_group, all_dihedrals, rot_bonds)

            symmetry_group_dihedral_dict[tuple(symmetry_group)] = symmetry_group_dihedral_list
            if len(symmetry_group) == 3:

                # angles_to_set[symmetry_group_dihedral_list[0]] = ([0.0, np.pi/3.0])
                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = ([0.0, np.pi/3.0])

                periodicities[tuple(symmetry_group_dihedral_list[0])] = 3
                dihedral_groups[3].extend([tuple(sorted(element)) for element in symmetry_group_dihedral_list])

            elif len(symmetry_group) == 2:

                # angles_to_set[symmetry_group_dihedral_list[0]] = ([0.0, np.pi/3.0])
                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = ([0.0, np.pi/2.0])

                periodicities[tuple(symmetry_group_dihedral_list[0])] = 2
                dihedral_groups[2].extend([tuple(sorted(element)) for element in symmetry_group_dihedral_list])

        return angles_to_set, periodicities, symmetry_group_dihedral_dict, dihedral_groups       
    
    
    
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
            qm_energy = np.array([qm_energy])
            qm_driver.ostream.unmute()
        
        elif isinstance(qm_driver, ExternalScfDriver):
            qm_energy = qm_driver.compute_energy(molecule, basis.get_main_basis_label())
            print('qm_energy', qm_energy)
        
        elif isinstance(qm_driver, ExternalExcitedStatesScfDriver):
            qm_energy = qm_driver.compute_energy(molecule, basis.get_main_basis_label())

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
        
        elif isinstance(grad_driver, ExternalGradientDriver):
            grad_driver.compute_gradient(molecule)
            qm_gradient = grad_driver.extract_gradients()
            qm_gradient = qm_gradient
        
        elif isinstance(grad_driver, ExternalExcitedStatesGradientDriver):
            grad_driver.compute_gradient(molecule)
            qm_gradient = grad_driver.extract_gradients()

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
        
        elif isinstance(hess_driver, ExternalHessianDriver):
            hess_driver.compute_analytical_hessian(molecule)
            qm_hessian = hess_driver.extract_hessians()

            
        elif isinstance(hess_driver, ExternalExcitedStatesHessianDriver):
            hess_driver.compute_analytical_hessian(molecule)
            qm_hessian = hess_driver.extract_hessians()

        if qm_hessian is None:
            error_txt = "Could not compute the QM Hessian. "
            error_txt += "Please define a QM Hessian driver."
            raise ValueError(error_txt)


        return qm_hessian



    def get_qm_potential_energy(self):
        """
        Returns the potential energy of the QM region.

        Args:
            context: The OpenMM context object.
        Returns:
            The potential energy of the QM region.
        """

        potential_energy = self.current_energy

        return potential_energy
    
    def output_file_writer(self, outputfile):
        """
        Writes the current simulation summary data (stored in self.gloabal_sim_informations)
        into the HDF5 file and then resets the dictionary.
        """
        
        valid_checkpoint = (outputfile and isinstance(outputfile, str))

        if valid_checkpoint:
            try:
                h5f = h5py.File(outputfile, 'a')
            except IOError:
                h5f = h5py.File(outputfile, 'w')
        

        # 1. Write temperatures
        ds_name = 'mol_labels'
        if ds_name in h5f:
            print('labels already exist')
        else:
            labels = self.molecule.get_labels()
            dt = h5py.string_dtype(encoding='utf-8')
            label_ds = h5f.create_dataset(
                ds_name,
                shape=(0,),
                maxshape=(None,),
                dtype=dt,
                chunks=True
            )
            old_size = label_ds.shape[0]
            new_size = old_size + labels.shape[0]
            label_ds.resize((new_size,))
            label_ds[old_size:new_size] = labels

        ds_name = "temperatures"
        if ds_name in h5f:
            temp_ds = h5f[ds_name]
        else:
            temp_ds = h5f.create_dataset(
                ds_name,
                shape=(0,),
                maxshape=(None,),
                dtype='float64',
                chunks=True
            )
        # Convert the list to a NumPy array and append
        new_temps = np.array(self.gloabal_sim_informations['temperatures'], dtype='float64')
        old_size = temp_ds.shape[0]
        new_size = old_size + new_temps.shape[0]
        temp_ds.resize((new_size,))
        temp_ds[old_size:new_size] = new_temps
        # 2. Write state labels
        ds_name = "state"
        # Adjust dtype as needed (here assuming integer states; if strings, you'll need a special dtype)
        if ds_name in h5f:
            state_ds = h5f[ds_name]
        else:
            state_ds = h5f.create_dataset(
                ds_name,
                shape=(0,),
                maxshape=(None,),
                dtype='int32',
                chunks=True
            )
        new_states = np.array(self.gloabal_sim_informations['state'], dtype='int32')
        old_size = state_ds.shape[0]
        new_size = old_size + new_states.shape[0]
        state_ds.resize((new_size,))
        state_ds[old_size:new_size] = new_states
        
        # 3. Write coordinates (assumed to be np.array of shape (n_atoms, 3))
        ds_name = "coordinates_ang"
        coord_list = self.gloabal_sim_informations['coordinates_ang']
        if len(coord_list) > 0:
            # Determine the number of atoms from the first snapshot.
            n_atoms = coord_list[0].shape[0]
            new_coords = np.stack(coord_list, axis=0)  # shape: (n_snapshots, n_atoms, 3)
            if ds_name in h5f:
                coord_ds = h5f[ds_name]
            else:
                coord_ds = h5f.create_dataset(
                    ds_name,
                    data=new_coords,
                    maxshape=(None, n_atoms, 3),
                    dtype=new_coords.dtype,
                    chunks=True
                )
                new_coords = None  # Already stored
            if new_coords is not None:
                old_size = coord_ds.shape[0]
                new_size = old_size + new_coords.shape[0]
                coord_ds.resize((new_size, n_atoms, 3))
                coord_ds[old_size:new_size, :, :] = new_coords
        # 4. Write state-specific data (for each state in self.roots_to_follow)
        for root in self.roots_to_follow:
            key = f'state_{root}'
            state_dict = self.gloabal_sim_informations[key]
            state_group = h5f.require_group(key)
            
            # 1. Store potential energies (assumed to be 1D arrays)
            for subkey in ['pot_energies']:
                if subkey in state_group:
                    ds = state_group[subkey]
                else:
                    ds = state_group.create_dataset(
                        subkey,
                        shape=(0,),
                        maxshape=(None,),
                        dtype='float64',
                        chunks=True
                    )
                new_data = np.array(state_dict[subkey], dtype='float64')
                old_size = ds.shape[0]
                new_size = old_size + new_data.shape[0]
                ds.resize((new_size,))
                ds[old_size:new_size] = new_data
            
            # 2. Store gradients (each entry is an array of shape (n_atoms, 3))
            subkey = 'gradients'
            new_gradients = state_dict[subkey]
            if new_gradients:  # Only proceed if there is new data
                # Determine number of atoms from the first gradient array.
                n_atoms = new_gradients[0].shape[0]
                # Stack new gradients along a new axis: shape becomes (n_snapshots, n_atoms, 3)
                grad_data = np.stack(new_gradients, axis=0)
                
                if subkey in state_group:
                    ds = state_group[subkey]
                else:
                    # Create the dataset with initial data
                    ds = state_group.create_dataset(
                        subkey,
                        data=grad_data,
                        maxshape=(None, n_atoms, 3),
                        dtype=grad_data.dtype,
                        chunks=True
                    )
                    grad_data = None  # Already stored
                # If the dataset already existed, append the new data.
                if grad_data is not None:
                    old_size = ds.shape[0]
                    new_size = old_size + grad_data.shape[0]
                    ds.resize((new_size, n_atoms, 3))
                    ds[old_size:new_size, :, :] = grad_data
    
    def read_simulation_summary(self, fname, roots_to_follow):
        """
        Reads a simulation summary from an HDF5 file.

        :param fname: Name of the HDF5 file.
        :param roots_to_follow: List of root states (e.g., [0, 1]) that were written.
        :return: A dictionary containing all read data.
        """
        data = {}

        with h5py.File(fname, 'r') as h5f:
            # 1. Read temperatures
            if 'temperatures' in h5f:
                data['temperatures'] = np.array(h5f['temperatures'])

            # 2. Read states
            if 'state' in h5f:
                data['state'] = np.array(h5f['state'])

            # 3. Read recoordinates
            if 'coordinates_ang' in h5f:
                data['coordinates_ang'] = np.array(h5f['coordinates_ang'])

            # 4. Read state-specific groups
            for root in roots_to_follow:
                state_key = f'state_{root}'
                if state_key in h5f:
                    state_group = h5f[state_key]
                    state_data = {}

                    if 'pot_energies' in state_group:
                        state_data['pot_energies'] = np.array(state_group['pot_energies'])

                    if 'gradients' in state_group:
                        state_data['gradients'] = np.array(state_group['gradients'])  # shape: (n_snapshots, n_atoms, 3)

                    data[state_key] = state_data

        return data
    
    def write_qm_energy_determined_points(self, molecules, qm_energies, im_energies, distances, state):           
        
        mode = 'a' if os.path.exists(self.reference_struc_energies_file) and os.path.getsize(self.reference_struc_energies_file) > 0 else 'w'
        print(f"Opening file in mode: {mode}")

        with open(self.reference_struc_energies_file, mode) as file:
            
            for i, dyn_mol in enumerate(molecules):
                print('i', i)
                current_xyz_string = dyn_mol.get_xyz_string()
                xyz_lines = current_xyz_string.splitlines()

                # Add energy info to comment line (second line)
                xyz_lines[1] += f' Energies  QM: {qm_energies[i]} IM: {im_energies[i]} Distance: {1.0}  State: {state}'
                updated_xyz_string = "\n".join(xyz_lines)

                file.write(f"{updated_xyz_string}\n\n")

    def extract_reference_structures(self, filename, roots):

        xyz_strings = []  
        qm_energies = []  
        im_energies = []
        distances = []
        states = []  
        
        reference_molecules_check = {root: {'molecules': [], 'qm_energies': [], 'im_energies': [], 'distances': []} for root in roots}
    
        with open(filename, 'r') as file:
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

                match = re.search(r"QM:\s*([\-\d\.]+)\s+IM:\s*([\-\d\.]+)\s+Distance:\s*([\-\d\.]+)\s+State:\s*([\-\d\.]+)", line) # group the individual energies to append them in a list
                if match:
                    qm_energies.append(float(match.group(1)))
                    im_energies.append(float(match.group(2)))
                    distances.append(float(match.group(3)))
                    states.append(float(match.group(4)))

                current_xyz.append(line) # still need to be stored as string because vlx reads it in as xyz string
            else:
                current_xyz.append(line)


        if current_xyz:
            xyz_strings.append('\n'.join(current_xyz))
        
        for i, state in enumerate(states):
            mol = Molecule.from_xyz_string(xyz_strings[i])

            reference_molecules_check[state]['molecules'].append(mol)
            reference_molecules_check[state]['qm_energies'].append(qm_energies[i])
            reference_molecules_check[state]['im_energies'].append(im_energies[i])
            reference_molecules_check[state]['distances'].append(distances[i])

        



        return reference_molecules_check

    def calculate_translation_coordinates(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates
    

    def calculate_distance_to_ref(self, current_coordinates, datapoint_coordinate):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                InterpolationDatapoint object
        """

        # First, translate the cartesian coordinates to zero
        target_coordinates = self.calculate_translation_coordinates(datapoint_coordinate)
        reference_coordinates = self.calculate_translation_coordinates(current_coordinates)

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

        distance_vector_core = (ref_structure_check - rotated_coordinates_core)
        distance_vector_core_norm = np.zeros(reference_coordinates.shape[0])

        for i in range(len(distance_vector_core_norm)):
            distance_vector_core_norm[i] += np.linalg.norm(distance_vector_core[i])

        return ref_structure_check, distance_core, distance_vector_core
        