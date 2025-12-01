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

from math import cos

from networkx import node_clique_number
from mpi4py import MPI
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import scipy
import h5py
import itertools
import re
import os, copy, math
import torch
import gpytorch
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
from. veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom, hartree_in_kjpermol
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .solvationbuilder import SolvationBuilder
from .optimizationdriver import OptimizationDriver

# Drivers
from .scfrestdriver import ScfRestrictedDriver
from .molecularbasis import MolecularBasis
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .tdaeigensolver import TdaEigenSolver
from .lreigensolver import LinearResponseEigenSolver
from .tddftgradientdriver import TddftGradientDriver
from .tddfthessiandriver import TddftHessianDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .localbayesresidual import LocalBayesResidual
from .alphaoptimizer import AlphaOptimizer
from .gprinterpolationdriver import GPRInterpolationDriver

with redirect_stderr(StringIO()) as fg_err:
    import geometric

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass


# --- put the custom reporter here ---
class GhostSafePDBReporter(app.PDBReporter):
    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True)
        n_atoms = simulation.topology.getNumAtoms()
        trimmed_positions = positions[:n_atoms]
        app.PDBFile.writeModel(simulation.topology, trimmed_positions, self._out, self._nextModel)
        self._nextModel += 1


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
        self.bias_force_reaction_idx = None
        self.bias_force_reaction_prop = None

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
        self.prev_dens_of_points = None

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
        self.gradient_rmsd_thrsh = None
        self.force_orient_thrsh = None
        self.collect_qm_points = None
        self.previous_energy_list = []
        self.non_core_symmetry_groups = []
        self.interpolation_settings = None
        self.allowed_molecules = None
        self.reference_struc_energies_file = None
        self.all_rot_bonds = None

        self.starting_temperature = None

        self.current_state = None
        self.starting_state = 0
        self.distance_thrsh = 0.1
        self.current_im_choice = None
        self.current_gradient = 0
        self.current_energy = 0
        self.point_checker = 1
        self.last_added = 0
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
        
        self.identfy_relevant_int_coordinates = (True, [])
        self.use_symmetry = True
        self.use_opt_confidence_radius = True
        self.qm_symmetry_dict = None
        self.add_gpr_model = True
        
        self.summary_output = 'summary_output.h5'
        self.coordinates_xyz = None

        self.velocities = []

        # Default value for the C-H linker distance
        self.linking_atom_distance = 1.0705 

        self.eq_bond_force_constants = None
            
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
                                equilibrate=True,
                                constraint_solute=True
                                )
            sys_builder.write_openmm_files(solute_ff=ff_gen)
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
                print('I am going in here', forcefield_files)

            # Load the PDB from the SystemBuilder
            self.pdb = app.PDBFile('system.pdb')

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


    def compute_tetrahedral_anchor(self, positions, Cbeta, Calpha, Hbeta1, Hbeta2, bond_length=0.109):
        def normalize(v):
            return v / np.linalg.norm(v)

        rC  = positions[Cbeta].value_in_unit(unit.nanometer)
        rCa = positions[Calpha].value_in_unit(unit.nanometer)
        rH1 = positions[Hbeta1].value_in_unit(unit.nanometer)
        rH2 = positions[Hbeta2].value_in_unit(unit.nanometer)

        v1 = normalize(rCa - rC)
        v2 = normalize(rH1 - rC)
        v3 = normalize(rH2 - rC)

        # Compute plane normal (sp2 plane)
        normal = normalize(np.cross(v2 - v1, v3 - v1))

        # Average substituent direction (points roughly into the sp2 plane)
        avg = normalize(v1 + v2 + v3)

        # The ghost should point mostly along -avg, but tilted out of the plane
        # Combine avg and normal with weights to get ~109.5° separation
        alpha = np.cos(np.deg2rad(109.5))   # ~ -0.333
        beta  = np.sqrt(1 - alpha**2)      # ~ 0.943
        vghost = normalize(alpha * avg + beta * normal)

        anchor = rC + bond_length * vghost
        return anchor


    def add_bias_force(self, atoms, force_constant, target=0.109, anchor=False):
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

        if len(atoms) == 2 and anchor:
            msg = f'Adding stretch force between atom {atoms[0]} and anchor {atoms[1]} with force constant {force_constant}.'
            self.ostream.print_info(msg)
            self.ostream.flush()
            for force in self.system.getForces():
                if isinstance(force, mm.NonbondedForce):
                    force.addParticle(0.0, 1.0, 0.0)  

            force = mm.CustomBondForce("0.5*k*(max(0, r-r0))^2")
            force.addGlobalParameter("k", force_constant)
            force.addGlobalParameter("r0", 0.109)  # nm
            force.addBond(atoms[0], atoms[1])
            self.system.addForce(force)
        
        elif len(atoms) == 2:
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

        if 'excitation_pulse' in dynamics_settings and dynamics_settings['excitation_pulse'] is not None:
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
        
        if 'grad_rmsd_thrsh' in dynamics_settings:
            self.gradient_rmsd_thrsh = dynamics_settings['grad_rmsd_thrsh']
        
        if 'force_orient_thrsh' in dynamics_settings:
            self.force_orient_thrsh = dynamics_settings['force_orient_thrsh']

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

        # Compute anchor from geometry
        if self.bias_force_reaction_prop is not None:

            if self.bias_force_reaction_prop[3]:

            
                ghost_index = self.system.addParticle(0.0)

                self.add_bias_force((self.bias_force_reaction_prop[0][0], ghost_index), self.bias_force_reaction_prop[1], anchor=self.bias_force_reaction_prop[3])

                
                Cbeta, Calpha, H1beta, H2beta = self.bias_force_reaction_idx[0], self.bias_force_reaction_idx[1], self.bias_force_reaction_idx[2],self.bias_force_reaction_idx[3]
                anchor = self.compute_tetrahedral_anchor(self.positions, Cbeta, Calpha, H1beta, H2beta)
    
            
                anchor_pos = np.array(anchor) * unit.nanometer

                self.positions = list(self.positions) + [anchor_pos]

            else:
                self.add_bias_force(self.bias_force_reaction_prop[0], self.bias_force_reaction_prop[1], self.bias_force_reaction_prop[2])   
        
        self.simulation = app.Simulation(self.topology, self.system, new_integrator)

        # Load the state if a restart file is provided
        if self.load_system is not None:
            self.simulation.loadState(self.load_system)
        
        else:
            self.simulation.context.setPositions(self.positions)

            # Set initial velocities if the ensemble is NVT or NPT
            if self.ensemble in ['NVT', 'NPT']:
                self.simulation.context.setVelocitiesToTemperature(self.temperature)
            # else:
            #     self.simulation.context.setVelocitiesToTemperature(150 * unit.kelvin)
        self.start_velocities = self.simulation.context.getState(getVelocities=True).getVelocities()
        # Set up reporting
        self.simulation.reporters.clear()

        if self.bias_force_reaction_prop is not None and self.bias_force_reaction_prop[3]:
            self.simulation.reporters.append(GhostSafePDBReporter(self.out_file, save_freq))
        else:
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
        self.allowed_molecules = {root: {'molecules': [], 'qm_energies': [], 'im_energies':[], 'qm_gradients':[], 'distances': []} for root in self.roots_to_follow}
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

        self.transition_information = []
        self.prev_dens_of_points = {root: 0 for root in self.roots_to_follow}
        self.bayes_models = {root: [] for root in self.roots_to_follow}
        self.qm_data_point_dict = {root: [] for root in self.roots_to_follow}
        self.qm_symmetry_datapoint_dict = {root: {} for root in self.roots_to_follow}
        self.qm_energies_dict = {root: [] for root in self.roots_to_follow}
        self.impes_drivers = {root: None for root in self.roots_to_follow}
        self.sampled_molecules = {root: {'molecules': [], 'im_energies': [], 'qm_energies': [], 'qm_gradients': [], 'distances': []} for root in self.roots_to_follow}
        
        self.allowed_molecules_to_check = {root: [] for root in self.roots_to_follow}

        if self.reference_struc_energies_file is not None:
            print('Extracting reference structures from', self.reference_struc_energies_file)
            self.sampled_molecules = self.extract_reference_structures(self.reference_struc_energies_file, self.roots_to_follow)
            self.allowed_molecules = self.extract_reference_structures(self.reference_struc_energies_file, self.roots_to_follow)
            self.last_added += len(self.sampled_molecules[self.roots_to_follow[0]]['molecules'])

        else:
             self.reference_struc_energies_file = 'ref_struc_energy.xyz'

        masses = self.molecule.get_masses().copy()
        masses_cart = np.repeat(masses, 3)
        inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)

        for root in self.roots_to_follow:
            # Dynamically create an attribute name
            attribute_name = f'impes_driver_{root}'
            # Initialize the object
            driver_object = InterpolationDriver(self.z_matrix)
            driver_object.update_settings(self.interpolation_settings[root])
            driver_object.eq_bond_force_constants = self.eq_bond_force_constants
            driver_object.impes_coordinate.eq_bond_force_constants = self.eq_bond_force_constants
            # print('Interpolation driver settings updated for root', root, self.eq_bond_force_constants)
            
            driver_object.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
            if root == 0:
                driver_object.symmetry_information = self.non_core_symmetry_groups['gs']
            elif root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
                driver_object.symmetry_information = self.non_core_symmetry_groups['gs']
            else:
                driver_object.symmetry_information = self.non_core_symmetry_groups['es']
            

            driver_object.use_symmetry = self.use_symmetry

            if len(self.sampled_molecules[root]['molecules']) != 0 and self.add_gpr_model:
                
                internal_coordaintes = []
                dE_gpr = []


                # Last molecule only
                x_list, y_list = [], []
                x_groups = None
                for idx, mol in enumerate(self.allowed_molecules[root]['molecules'][:]):
                    # 1) ICs + groups (groups should be identical for every mol)
                    X_single, groups, types = driver_object.compute_ic_and_groups(mol.get_coordinates_in_bohr())
                    if x_groups is None:
                        x_groups = groups
                    else:
                        # sanity: same layout across all mols
                        assert all(np.array_equal(x_groups[k], groups[k]) for k in x_groups.keys()), "IC groups changed!"
                    x_list.append(np.asarray(X_single, dtype=np.float64).reshape(1, -1))
                    # 2) energies -> per-atom ΔE (use mol, not molecule!)
                    per_atom = len(mol.get_labels())
                    y_list.append(
                        np.asarray(
                            (self.allowed_molecules[root]['qm_energies'][idx] - self.allowed_molecules[root]['im_energies'][idx])
                            * hartree_in_kcalpermol() / per_atom,
                            dtype=np.float64
                        ).reshape(-1)
                    )
                # 3) stack batch
                X_np = np.vstack(x_list)              # (N,D)
                y_np = np.concatenate(y_list, axis=0) # (N,)
                    

                gpr_driver = GPRInterpolationDriver(X_np, x_groups, y_np)
                gpr_driver.refit(steps="auto", verbose=False)
                base = gpr_driver.model.covar_module.base_kernel
                groups = base.kernels if hasattr(base, "kernels") else [base]
                ells = [float(k.lengthscale.detach().cpu().view(-1)[0]) for k in groups]
                print("group ℓ:", ells)

                # outputscale & noise
                print("outputscale:", float(gpr_driver.model.covar_module.outputscale.detach().cpu()))
                print("noise:", float(gpr_driver.likelihood.noise.detach().cpu()))
                driver_object.gpr_intdriver = gpr_driver 

            im_labels, _ = driver_object.read_labels()
            print('beginning labels', im_labels, self.add_gpr_model)
            self.qm_data_points = []
            self.qm_energies = []
            old_label = None

            for label in im_labels:
                if '_symmetry' not in label:
                    qm_data_point = InterpolationDatapoint(self.z_matrix)
                    qm_data_point.update_settings(self.interpolation_settings[root])
                    qm_data_point.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)
                    qm_data_point.inv_sqrt_masses = inv_sqrt_masses
                    self.qm_energies_dict[root].append(qm_data_point.energy)
                    self.qm_data_point_dict[root].append(qm_data_point)

                    self.sorted_state_spec_im_labels[root].append(label)
                    self.prev_dens_of_points[root] += 1
                    old_label = qm_data_point.point_label
                    driver_object.qm_symmetry_data_points[old_label] = [qm_data_point]
                    self.qm_symmetry_datapoint_dict[root][old_label] = [qm_data_point]

                else:
                    symmetry_data_point = InterpolationDatapoint(self.z_matrix)
                    symmetry_data_point.update_settings(self.interpolation_settings[root])
                    symmetry_data_point.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)
                    
                    driver_object.qm_symmetry_data_points[old_label].append(symmetry_data_point)
                    self.qm_symmetry_datapoint_dict[root][old_label].append(symmetry_data_point)

            print(driver_object.qm_symmetry_data_points)       
            # Set the object as an attribute of the instance
            setattr(self, attribute_name, driver_object)
            # Append the object to the list
            self.impes_drivers[root] = driver_object


        self.current_state = self.roots_to_follow[0]
 
        self.swap_back = True
        self.prev_state = self.current_state
        start_time = time()
        self.step = 0 
        self.prev_dE_gpr = 0.0
        self.last_gpr_addition = 0

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

            self.step += 1

            # self._update_induction_embedding(phase == 'periodic')
            self.simulation.step(1)
            self.prev_state = self.current_state
            
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

        if not hasattr(self, "force_IM_embed"):
            self.force_IM_embed = mm.CustomExternalForce("-fx*x - fy*y - fz*z")
            self.force_IM_embed.addPerParticleParameter("fx")
            self.force_IM_embed.addPerParticleParameter("fy")
            self.force_IM_embed.addPerParticleParameter("fz")
            for a in qm_group:
                self.force_IM_embed.addParticle(a, (0.0, 0.0, 0.0))
            self.force_IM_embed.setForceGroup(9)
            self.system.addForce(self.force_IM_embed)

        # 2) Constant-force container for MM back-reaction (group 10)
        if mm_group and not hasattr(self, "force_MM_embed"):
            self.force_MM_embed = mm.CustomExternalForce("-fx*x - fy*y - fz*z")
            self.force_MM_embed.addPerParticleParameter("fx")
            self.force_MM_embed.addPerParticleParameter("fy")
            self.force_MM_embed.addPerParticleParameter("fz")
            for b in mm_group:
                self.force_MM_embed.addParticle(b, (0.0, 0.0, 0.0))
            self.force_MM_embed.setForceGroup(10)
            self.system.addForce(self.force_MM_embed)

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
    
    def _update_induction_embedding(self, periodic: bool):
        # --- fetch positions (nm) ---
        state = self.simulation.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)

        # --- get charges for MM atoms (e) once and cache ---
        if not hasattr(self, "_mm_charge_vec"):
            nb = [f for f in self.system.getForces() if isinstance(f, mm.NonbondedForce)][0]
            q = np.zeros(self.system.getNumParticles())
            for i in range(nb.getNumParticles()):
                qi, sigma, eps = nb.getParticleParameters(i)
                q[i] = qi._value   # elementary charge
            self._mm_charge_vec = q

        # --- Coulomb constant in kJ/mol·nm·e^-2 ---
        Kc = 138.935456

        # --- Reaction Field constants if periodic & using your RF coupling ---
        if periodic:
            nb = [f for f in self.system.getForces() if isinstance(f, mm.NonbondedForce)][0]
            rfDielectric = nb.getReactionFieldDielectric()
            rc = float(self.cutoff)  # nm
            krf = (1.0/(rc**3)) * (rfDielectric - 1.0) / (2.0*rfDielectric + 1.0)
            # crf drops out of the field (it’s a constant in potential), see derivation below
        else:
            krf = 0.0

        # --- work arrays ---
        F_im = np.zeros((len(self.qm_atoms), 3))                # forces to add on IM atoms
        F_mm = np.zeros((len(self.mm_subregion), 3))            # forces to add on MM atoms
        U_ind = 0.0

        # --- user-supplied isotropic polarizabilities for IM sites (nm^3 units) ---
        # Suggestion: start with element-wise constants or a single scalar.
        # e.g., self.alpha_iso[a] already defined in nm^3 (1 Å^3 = 1e-3 nm^3)
        alpha = self.alpha_iso   # dict or array indexed by atom id

        # Helper to map mm_subregion index -> row id in F_mm
        mm_index_of = {b: k for k, b in enumerate(self.mm_subregion)}

        # --- loop over IM sites, accumulate field and forces ---
        for ia, a in enumerate(self.qm_atoms):
            ra = np.array(pos[a], dtype=float)
            Ea = np.zeros(3, dtype=float)

            # field from all MM charges
            for b in self.mm_subregion:
                qb = self._mm_charge_vec[b]
                if qb == 0.0:
                    continue
                R = ra - np.array(pos[b], dtype=float)
                r2 = float(np.dot(R, R)) + 1e-24
                r = np.sqrt(r2)
                R3 = r2 * r

                # Electric field (Reaction Field consistent):
                # E = Kc * q * ( R/r^3  -  2*krf * R )
                Ea += Kc * qb * (R / R3 - 2.0 * krf * R)

            # induction energy contribution at site a
            alpha_a = float(alpha[a])  # nm^3
            U_ind += -0.5 * alpha_a * float(np.dot(Ea, Ea))

            # Now build forces from d/dR of U = -1/2 alpha |E|^2
            for b in self.mm_subregion:
                qb = self._mm_charge_vec[b]
                if qb == 0.0:
                    continue
                R = ra - np.array(pos[b], dtype=float)
                r2 = float(np.dot(R, R)) + 1e-24
                r = np.sqrt(r2)
                RR = np.outer(R, R)

                # J(R)·E = ∂E/∂R · E  with Reaction Field:
                # ∂E/∂R = Kc*q * [ (r^2 I - 3 RR^T)/r^5  -  2*krf I ]
                JdotE = Kc * qb * (((r2*np.eye(3) - 3.0*RR) @ Ea) / (r2*r2*r) - 2.0*krf * Ea)

                # forces: F_a += +alpha * J·E ; F_b += -alpha * J·E
                Fa =  alpha_a * JdotE
                Fb = -alpha_a * JdotE

                F_im[ia] += Fa
                F_mm[mm_index_of[b]] += Fb

        # --- push to context (don’t forget to call updateParametersInContext) ---
        for k, a in enumerate(self.qm_atoms):
            self.force_IM_embed.setParticleParameters(k, a, F_im[k])
        for k, b in enumerate(self.mm_subregion):
            self.force_MM_embed.setParticleParameters(k, b, F_mm[k])

        self.force_IM_embed.updateParametersInContext(self.simulation.context)
        self.force_MM_embed.updateParametersInContext(self.simulation.context)

        # store energy for your own total-energy bookkeeping
        self.last_U_induction = U_ind     # kJ/mol
    
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
            new_molecule.set_charge(self.molecule.get_charge())
            new_molecule.set_multiplicity(self.molecule.get_multiplicity())

        else:
            # Atom labels for the QM region
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            new_molecule.set_charge(self.molecule.get_charge())
            new_molecule.set_multiplicity(self.molecule.get_multiplicity())
            self.unique_molecules.append(new_molecule)
        
        for root in self.roots_to_follow:
            self.impes_drivers[root].qm_data_points = self.qm_data_point_dict[root]
            self.impes_drivers[root].compute(new_molecule)


        if self.nstates > 1:
            for root_1 in range(0, self.nstates):
                for root_2 in range(root_1 + 1, self.nstates):

                    potential_kjmol = self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
                    potential_kjmol_2 = self.impes_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
                    self.transition_information.append((self.roots_to_follow[root_1], self.roots_to_follow[root_2], potential_kjmol - potential_kjmol_2))
                    print(f'compare the energies between roots: {self.roots_to_follow[root_1]} -> {self.roots_to_follow[root_2]}', potential_kjmol_2 - potential_kjmol)

                    if self.calc_NAC == True and np.linalg.norm(self.velocities_np[-1]) > 0.0:
                        current_NAC = self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.NAC.flatten()
                        current_velocitites = self.velocities_np[-1].flatten() * 4.566180e-4

                        hopping_potential = np.exp(-abs((np.pi/(4)) * (( self.impes_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy - self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy ) / np.linalg.multi_dot([current_NAC, current_velocitites]))))
                        print('#######################', '\n\n', hopping_potential, potential_kjmol_2 - potential_kjmol, '\n\n', '#######################')

                    if potential_kjmol_2 - potential_kjmol < 5 and self.swap_back:
                        # Choose a random integer between 0 and 1
                        random_integer = random.randint(0, 1)
                        self.swap_back = False
                        random_integer = 1.0
                        if random_integer == 1 and self.current_state == self.roots_to_follow[root_1]:
                            self.current_state = self.roots_to_follow[root_2]
                        elif random_integer == 1 and self.current_state == self.roots_to_follow[root_2]:
                            self.current_state = self.roots_to_follow[root_1]
                        break

        
        else:
            self.current_state = self.roots_to_follow[0]

        if self.excitation_pulse is not None and len(self.roots_to_follow) > 1 and self.point_checker == self.excitation_pulse[0] and self.current_state != self.excitation_pulse[1] and self.swap_back:
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
        new_molecule.set_charge(self.molecule.get_charge())
        new_molecule.set_multiplicity(self.molecule.get_multiplicity())
        
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
            
            if self.skipping_value == 0 and self.point_checker + 1 > self.collect_qm_points:
        
                mean_bond_rmsd = np.mean(np.array(self.impes_drivers[self.current_state].bond_rmsd))
                mean_angle_rmsd = np.mean(np.array(self.impes_drivers[self.current_state].angle_rmsd))
                mean_dihedral_rmsd = np.mean(np.array(self.impes_drivers[self.current_state].dihedral_rmsd))

                if mean_bond_rmsd > 0.1 or mean_angle_rmsd > 0.1 or mean_dihedral_rmsd > 1.0:

                    K = 5  # Number of closest previous matches to cache
                    threshold = self.distance_thrsh - 0.05
                    new_coords = new_molecule.get_coordinates_in_bohr()
                    n_atoms = len(new_molecule.get_labels())
                    sym_group = self.non_core_symmetry_groups['gs' if self.current_state == 0 or self.drivers['es'] is not None and self.drivers['es'][0].spin_flip and self.current_state == 1 else 'es'][4]
                    molecule_list = self.allowed_molecules_to_check[self.current_state]

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

                if self.impes_drivers[self.current_state].gpr_intdriver is not None and len(self.allowed_molecules[self.current_state]['molecules']) > 5:
                
                    X_single_list, groups, type = self.impes_drivers[self.current_state].compute_ic_and_groups(new_molecule.get_coordinates_in_bohr())
                    
                    X_new = np.asarray(
                                X_single_list,
                                dtype=np.float64
                            ).reshape(1, -1)
                    
                    mean, std = self.impes_drivers[self.current_state].gpr_intdriver.predict(X_new, return_std=True) 
                    mean_latent, std_latent = self.impes_drivers[self.current_state].gpr_intdriver.predict_latent(X_new, return_std=True) 

                    ### second check

                    trigger, information = self.impes_drivers[self.current_state].gpr_intdriver.should_qm_check(X_new, T=self.energy_threshold)
                    
                    error_p_variance = abs(mean) + std
                    
                    error_proxy = float(np.abs(mean).ravel()[0] + std.ravel()[0])
                    # print('Error proxy', error_proxy, trigger, information)
                    Xt, yt = self.impes_drivers[self.current_state].gpr_intdriver._transform(self.impes_drivers[self.current_state].gpr_intdriver.X_raw, self.impes_drivers[self.current_state].gpr_intdriver.y_raw)
                      

                    self.add_a_point = False
                    if trigger:
                        self.add_a_point = True 
                    self.prev_dE_gpr = error_p_variance


            else:
                self.skipping_value -= 1
                if self.skipping_value < 0:
                    self.skipping_value = 0
            
            # if self.step % 1 == 0:
            #     self.add_a_point = True
            
            if self.current_state != self.prev_state:
                self.add_a_point = True

            self.point_checker += 1 
            self.last_gpr_addition += 1

            print('current state', self.current_state)
            if self.add_a_point == True:
                self.point_correlation_check(new_molecule)
                self.last_gpr_addition = 0
            if self.point_checker == 0:            
                
                
                self.swap_back = True
                self.current_state = self.starting_state
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
                new_molecule.set_charge(self.molecule.get_charge())
                new_molecule.set_multiplicity(self.molecule.get_multiplicity())
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
            new_molecule.set_charge(self.molecule.get_charge())
            new_molecule.set_multiplicity(self.molecule.get_multiplicity())
            gradient_2 = self.update_gradient(qm_positions)
            force = -np.array(gradient_2) * conversion_factor
            
            self.state_specific_molecules[self.current_state].append(new_molecule)
            self.coordinates_xyz.append(qm_positions * 10)
            for i, atom_idx in enumerate(self.qm_atoms):
                self.system.getForce(self.qm_force_index).setParticleParameters(i, atom_idx, force[i])
            self.system.getForce(self.qm_force_index).updateParametersInContext(context)

        if self.current_state < self.prev_state and 1 == 2:
            
            qm = self.get_qm_potential_energy()
            self.qm_potentials.append(qm)

            # QM/MM interactions
            qm_mm = context.getState(getEnergy=True, groups={1,2}).getPotentialEnergy()

            # MM region 
            mm = context.getState(getEnergy=True, groups={3,4,5,6,7}).getPotentialEnergy()

            # Total potential energy
            pot = qm * unit.kilojoules_per_mole + qm_mm + mm

            # Kinetic energy
            kinetic = context.getState(getEnergy=True).getKineticEnergy()
    
            # Total energy
            total = pot + kinetic
            E_tot_diff = np.sqrt(abs(self.total_energies[-1] - total.value_in_unit(unit.kilojoules_per_mole)) / kinetic.value_in_unit(unit.kilojoules_per_mole))
    
            # Assume you have a Simulation object called 'simulation'
            state = context.getState(getVelocities=True)
            velocities = state.getVelocities(asNumpy=True)  # unit: nanometers/picosecond

            scaling_factor = E_tot_diff  # for example, to reduce velocities by 20%
            velocities *= scaling_factor

            context.setVelocities(velocities)
        
        if len(self.roots_to_follow) > 1 and self.swap_back == False:# and self.temperatures[-1] > 50:
            state = context.getState(getVelocities=True)
            velocities = state.getVelocities(asNumpy=True)  # unit: nanometers/picosecond

         
            velocities *= 0.01

            context.setVelocities(velocities)
            print('I am setting back the velocities')

        elif len(self.roots_to_follow) > 1 and self.swap_back == False and self.temperatures[-1] < self.temperatures[0]:
            self.swap_back = True
            self.point_checker = 0
        
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
        
        current_state_difference = {}
        state_specific_energies = {}
        state_specific_gradients = {}
        current_state_to_state_difference = {}
        self.add_a_point = False

        all_combinations = list(itertools.combinations(self.roots_to_follow, 2))

        for comb in all_combinations:
            current_state_to_state_difference[comb] = [None, None]
        for root in self.roots_to_follow:
            state_specific_energies[root] = None
            state_specific_gradients[root] = None
            current_state_difference[root] = [None, None, None]

        addition_of_state_specific_points = []

        for state, key in enumerate(self.drivers.keys()):
            
            drivers = None
            identification_state = state
            
            if state == 0 and self.drivers[key] is not None and 0 in self.roots_to_follow:
                drivers = self.drivers[key]



   
            elif state == 1 and self.drivers['es'] is not None and any(x > 0 for x in self.roots_to_follow):
                drivers = self.drivers[key]
                if 0 not in self.roots_to_follow:
                    identification_state = 0

            else:
                continue

            natms = len(molecule.get_labels())
            print('############# Energy is QM claculated ############')
            print('identifaction state', identification_state, 'current state', self.current_state)
            current_basis = None
            if self.current_state == 0:
                current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])
            else:
                current_basis = MolecularBasis.read(molecule, self.basis_set_label['es'])
            
   
            qm_energy, scf_results, rsp_results = self.compute_energy(drivers[0], molecule, current_basis)
            print(qm_energy)
            if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                qm_energy = [qm_energy[identification_state + state]]


            gradients = self.compute_gradient(drivers[1], molecule, current_basis, scf_results, rsp_results)

            for e_idx in range(len(qm_energy)):
                grad = gradients[e_idx]         # (3N,)

                print('weights', self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].weights)
                energy_difference = (abs(qm_energy[e_idx] - self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.energy))

     
                gradient_difference = (grad - self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.gradient) * hartree_in_kcalpermol() * bohr_in_angstrom()
                # print('predicted Grad', self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.gradient)
                rmsd_gradient    = np.sqrt((gradient_difference**2).mean())
                
                cos_theta = 1.0
                if rmsd_gradient > 0.5 and energy_difference > 0.01:
                    
                    gq = grad.ravel()
                    gi = self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.gradient.ravel()
                    
                    # print('gq', gq, 'gi', gi)
                    cos_theta = 0.0
                    significant_mask = (np.abs(gq) > 1e-4) | (np.abs(gi) > 1e-4)
                    if np.sum(significant_mask) == 0:
                        cos_theta = 1.0
                    # print('significant components', np.sum(significant_mask))
                    gq_clean = gq[significant_mask]
                    gi_clean = gi[significant_mask]

                    norm_gq_clean = np.linalg.norm(gq_clean)
                    norm_gi_clean = np.linalg.norm(gi_clean)
                    
                    if norm_gq_clean == 0 and norm_gi_clean == 0:
                        cos_theta = 1.0
                    
                    cos_theta = np.dot(gq_clean, gi_clean) / (norm_gq_clean * norm_gi_clean)
                # cos_theta_old = np.dot(gq, gi) / (np.linalg.norm(gq) * np.linalg.norm(gi))
                # print('gq', cos_theta, cos_theta_old)
                
                print('Energy difference', energy_difference, energy_difference * hartree_in_kcalpermol(), 'kcal/mol', 'energy differences rmsd', energy_difference / natms * hartree_in_kcalpermol())
                print('gradients alignment', cos_theta, 'rmsd gradient', rmsd_gradient, 'kcal/mol/angstrom')


                state_specific_energies[self.roots_to_follow[identification_state + e_idx]] = [qm_energy[e_idx], self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.energy]
                current_state_difference[self.roots_to_follow[identification_state + e_idx]][0] = energy_difference / natms * hartree_in_kcalpermol()
                current_state_difference[self.roots_to_follow[identification_state + e_idx]][1] = rmsd_gradient
                current_state_difference[self.roots_to_follow[identification_state + e_idx]][2] = cos_theta
                state_specific_gradients[self.roots_to_follow[identification_state + e_idx]] = [grad, self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.gradient]

        if current_state_difference[self.current_state][0] > self.energy_threshold or current_state_difference[self.current_state][1] > self.gradient_rmsd_thrsh or current_state_difference[self.current_state][2] < self.force_orient_thrsh:
            print('The point would be added due to the thresholds')
            addition_of_state_specific_points.append(self.current_state)

        for comb in all_combinations:
            current_state_to_state_difference[comb][0] = abs((state_specific_energies[comb[0]][0] - state_specific_energies[comb[1]][0]) * hartree_in_kjpermol())
            current_state_to_state_difference[comb][1] = abs((state_specific_energies[comb[0]][1] - state_specific_energies[comb[1]][1]) * hartree_in_kjpermol())
            print('current state to state', comb, current_state_to_state_difference[comb])
            if current_state_to_state_difference[comb][0] < 30.0 and abs(current_state_to_state_difference[comb][1] - current_state_to_state_difference[comb][0]) > 10.0:
                if comb[0] not in addition_of_state_specific_points:
                    addition_of_state_specific_points.append(comb[0])
                if comb[1] not in addition_of_state_specific_points:
                    addition_of_state_specific_points.append(comb[1])

        if len(addition_of_state_specific_points) > 0:

            for state in self.roots_to_follow:
                if self.use_opt_confidence_radius[0] and len(self.allowed_molecules[state]['molecules']) > 20 and 1 == 2:
                    self.add_a_point = True
                    if not self.expansion:
                        length_vectors = (self.impes_drivers[-1].impes_coordinate.cartesian_distance_vector(self.qm_data_points[0]))
                        rmsd = (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels())) * bohr_in_angstrom())
                        self.expansion_molecules.append((molecule, energy_difference / natms * hartree_in_kcalpermol(), rmsd, (np.linalg.norm(length_vectors))))
                        self.last_point_added = self.point_checker - 1
                        self.point_checker = 0
                
                else:
                    self.add_a_point = True
                    #TODO: Check if GPR implementation is correctly implemented and updates correctly considering the GPRinterpolationdriver class
                    if self.add_gpr_model:

                       for root in self.roots_to_follow:

                            if len(self.allowed_molecules[state]['molecules']) < 2:
                                continue

                            x_list, y_list = [], []
                            interpolation_driver = InterpolationDriver(self.z_matrix)
                            interpolation_driver.update_settings(self.interpolation_settings[root])
                            interpolation_driver.symmetry_information = self.impes_drivers[root].symmetry_information
                            interpolation_driver.qm_symmetry_data_points = self.impes_drivers[root].qm_symmetry_data_points
                            interpolation_driver.distance_thrsh = 1000
                            interpolation_driver.exponent_p = self.impes_drivers[root].exponent_p
                            interpolation_driver.print = False
                            interpolation_driver.qm_data_points = self.impes_drivers[root].qm_data_points
                            interpolation_driver.calc_optim_trust_radius = True

                            x_groups = None

                            for idx, mol in enumerate(self.allowed_molecules[root]['molecules']):
                                X_single, groups, types = self.impes_drivers[root].compute_ic_and_groups(mol.get_coordinates_in_bohr())
                                if x_groups is None:
                                    x_groups = groups
                                else:                           
                                    assert all(np.array_equal(x_groups[k], groups[k]) for k in x_groups.keys()), "IC groups changed!"

                                x_list.append(np.asarray(X_single, dtype=np.float64).reshape(1, -1))

                 
                                interpolation_driver.compute(mol)
                                im_energy = interpolation_driver.get_energy()
                                self.allowed_molecules[root]['im_energies'][idx] = im_energy

                                per_atom = len(mol.get_labels())
                                y_list.append(
                                    np.asarray(
                                        (self.allowed_molecules[root]['qm_energies'][idx] - im_energy)
                                        * hartree_in_kcalpermol() / per_atom,
                                        dtype=np.float64
                                    ).reshape(-1)
                                )

   
                            X_np = np.vstack(x_list)           
                            y_np = np.concatenate(y_list, axis=0) 
                            new_idx = X_np.shape[0] - 1

                            drv = self.impes_drivers[root].gpr_intdriver
                            
                            if drv is not None and len(self.allowed_molecules[root]['molecules']) > 2:
                                X_t = torch.from_numpy(X_np).to(drv.dtype).to(drv.device)
                                y_t = torch.from_numpy(y_np).to(drv.dtype).to(drv.device)
                                # Current training set (exclude the new point)
                                Xtr = X_t[:-1, :]
                                ytr = y_t[:-1]
                                Xq  = X_t[new_idx:new_idx+1, :]    # (1,D)
                                y_true = float(y_t[new_idx].cpu().numpy())
                                # Transform with the CURRENT drv scalers (not refit yet)
                                Xtr_std, ytr_std = drv._transform(Xtr, ytr)
                                Xq_std           = drv._transform(Xq)
                                # Latent prediction out-of-sample (using current drv trained on X_raw/y_raw)
                                mu_lat, sig_lat = drv.predict_latent(Xq, return_std=True)  # add this method if you haven’t
                                mu_lat = float(mu_lat.squeeze()); sig_lat = float(sig_lat.squeeze())
                                # Effective distance of query to current training set
                                d_norm = drv._compute_d_norm(Xq_std, Xtr_std, groups)  # ← normalized distance
                                d_thresh = drv.get_distance_threshold(1.0)
                                s_max       = drv.max_kernel_similarity(Xq_std, Xtr_std, topk=5)
                                # Scalar absolute residual for THIS query (QM label already available)
                                abs_res = abs(y_true - mu_lat)
                                
                                # Update the dynamic threshold with this single pair
                                drv.record_distance_event(d_norm, True, abs_res, T=self.energy_threshold)
                                drv.update_distance_threshold()
                                
                                drv.record_similarity_event(s_max, True, abs_res, T=self.energy_threshold)
                                drv.update_similarity_threshold()
                                # ---------- 2) NOW add the point and refit (warm start) ----------
                                drv.replace_data(X_t, y_t, clone=False)
                                drv.refit(steps="auto", verbose=False)
                                # (optional) log query-centric diagnostics after refit
                                print(f"[update] d_norm={d_norm:.3f}, |residual|={abs_res:.4g}, "
                                    f"mu_lat={mu_lat:.4g}, sigma_lat={sig_lat:.4g}, T={self.energy_threshold}")
                                
                            # 4) predict for the last structure (shape (1,D))
                            mean, std = self.impes_drivers[root].gpr_intdriver.predict(X_np[-1:,:], return_std=True)
                            print('### Predictions:', mean.ravel().tolist(), std.ravel().tolist())
                            # 5) (optional) inspect kernel hyperparameters to judge reasonableness
                                              
            if self.add_a_point:    
            
                if not self.expansion:
                    length_vectors = (self.impes_drivers[-1].impes_coordinate.cartesian_distance_vector(self.qm_data_points[0]))
                    rmsd = (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels())) * bohr_in_angstrom())
                    self.expansion_molecules.append((molecule, energy_difference * hartree_in_kcalpermol(), rmsd, (np.linalg.norm(length_vectors))))
                    self.last_point_added = self.point_checker - 1
                    self.point_checker = 0
        else:
            for root in self.roots_to_follow:
                self.allowed_molecules[root]['molecules'].append(molecule)
                self.allowed_molecules[root]['im_energies'].append(self.impes_drivers[root].impes_coordinate.energy)
                self.allowed_molecules[root]['qm_energies'].append(state_specific_energies[root][0])
                self.allowed_molecules[root]['qm_gradients'].append(state_specific_gradients[root][0])
                self.allowed_molecules_to_check[root].append(molecule)

                if self.add_gpr_model:
                    # Last molecule only
                    x_list, y_list = [], []
                    interpolation_driver = InterpolationDriver(self.z_matrix)
                    interpolation_driver.update_settings(self.interpolation_settings[root])
                    interpolation_driver.symmetry_information = self.impes_drivers[root].symmetry_information
                    interpolation_driver.qm_symmetry_data_points = self.impes_drivers[root].qm_symmetry_data_points
                    interpolation_driver.distance_thrsh = 1000
                    interpolation_driver.eq_bond_force_constants = self.eq_bond_force_constants
                    interpolation_driver.impes_coordinate.eq_bond_force_constants = self.eq_bond_force_constants
                    interpolation_driver.exponent_p = self.impes_drivers[root].exponent_p
                    interpolation_driver.print = False
                    interpolation_driver.qm_data_points = self.impes_drivers[root].qm_data_points
                    interpolation_driver.calc_optim_trust_radius = False
                    x_groups = None
                    for idx, mol in enumerate(self.allowed_molecules[root]['molecules']):
                        # 1) ICs + groups (groups should be identical for every mol)
                        X_single, groups, types = self.impes_drivers[root].compute_ic_and_groups(mol.get_coordinates_in_bohr())
                        if x_groups is None:
                            x_groups = groups
                        else:
                            # sanity: same layout across all mols
                            assert all(np.array_equal(x_groups[k], groups[k]) for k in x_groups.keys()), "IC groups changed!"
                        x_list.append(np.asarray(X_single, dtype=np.float64).reshape(1, -1))
                        # 2) energies -> per-atom ΔE (use mol, not molecule!)
                        interpolation_driver.compute(mol)
                        im_energy = interpolation_driver.get_energy()
                        self.allowed_molecules[root]['im_energies'][idx] = im_energy
                        per_atom = len(mol.get_labels())
                        y_list.append(
                            np.asarray(
                                (self.allowed_molecules[root]['qm_energies'][idx] - im_energy)
                                * hartree_in_kcalpermol() / per_atom,
                                dtype=np.float64
                            ).reshape(-1)
                        )
                    # 3) stack batch
                    X_np = np.vstack(x_list)              # (N,D)
                    y_np = np.concatenate(y_list, axis=0) # (N,)
                    new_idx = X_np.shape[0] - 1
                    drv = self.impes_drivers[root].gpr_intdriver
                    if drv is None:
                        # First time: create once, with grouped kernel
    
                        gpr_driver = GPRInterpolationDriver(X_np, x_groups, y_np)  # your class builds grouped kernel inside __init__
                        # Fit hypers at least once so lengthscales/noise move off defaults
                        gpr_driver.refit()
                        self.impes_drivers[root].gpr_intdriver = gpr_driver
                        
                    else:
                        X_t = torch.from_numpy(X_np).to(drv.dtype).to(drv.device)
                        y_t = torch.from_numpy(y_np).to(drv.dtype).to(drv.device)
                        # Current training set (exclude the new point)
                        Xtr = X_t[:-1, :]
                        ytr = y_t[:-1]
                        Xq  = X_t[new_idx:new_idx+1, :]    # (1,D)
                        y_true = float(y_t[new_idx].cpu().numpy())
                        # Transform with the CURRENT drv scalers (not refit yet)
                        Xtr_std, ytr_std = drv._transform(Xtr, ytr)
                        Xq_std           = drv._transform(Xq)
                        # Latent prediction out-of-sample (using current drv trained on X_raw/y_raw)
                        mu_lat, sig_lat = drv.predict_latent(Xq, return_std=True)  # add this method if you haven’t
                        mu_lat = float(mu_lat.squeeze()); sig_lat = float(sig_lat.squeeze())
                        # Effective distance of query to current training set
                        d_norm = drv._compute_d_norm(Xq_std, Xtr_std, groups)  # ← normalized distance
                        d_thresh = drv.get_distance_threshold(1.0)
                        s_max       = drv.max_kernel_similarity(Xq_std, Xtr_std, topk=5)
                        # Scalar absolute residual for THIS query (QM label already available)
                        abs_res = abs(y_true - mu_lat)
                        
                        # Update the dynamic threshold with this single pair
                        drv.record_distance_event(d_norm, False, abs_res, T=self.energy_threshold)
                        drv.update_distance_threshold()
                        
                        drv.record_similarity_event(s_max, False, abs_res, T=self.energy_threshold)
                        drv.update_similarity_threshold()
                        # ---------- 2) NOW add the point and refit (warm start) ----------
                        drv.replace_data(X_t, y_t, clone=False)
                        drv.refit(steps="auto", verbose=False)
                        # (optional) log query-centric diagnostics after refit
                        print(f"[update] d_norm={d_norm:.3f}, |residual|={abs_res:.4g}, "
                            f"mu_lat={mu_lat:.4g}, sigma_lat={sig_lat:.4g}, T={self.energy_threshold}")
                        
                    # 4) predict for the last structure (shape (1,D))
                    mean, std = self.impes_drivers[root].gpr_intdriver.predict(X_np[-1:,:], return_std=True)
                    print('### Predictions:', mean.ravel().tolist(), std.ravel().tolist())
                    # 5) (optional) inspect kernel hyperparameters to judge reasonableness

                self.write_qm_energy_determined_points(self.allowed_molecules[root]['molecules'][self.last_added: ],
                                                    self.allowed_molecules[root]['qm_energies'][self.last_added:],
                                                    self.allowed_molecules[root]['qm_gradients'][self.last_added:],
                                                    self.allowed_molecules[root]['im_energies'][self.last_added:],
                                                    [],
                                                    root)
                self.last_added = len(self.allowed_molecules[root]['molecules'])


            if self.use_opt_confidence_radius[0] and len(self.allowed_molecules[self.current_state]['molecules']) >= 2 and self.density_around_data_point[0][self.current_state] > 1 and self.density_around_data_point[0][self.current_state] % 5 == 0 and self.prev_dens_of_points[self.current_state] != self.density_around_data_point[0][self.current_state]:
                self.prev_dens_of_points[self.current_state] = self.density_around_data_point[0][self.current_state]        
                trust_radius = None
                sym_dict = self.non_core_symmetry_groups['gs']
                if self.current_state > 0 or root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
                    sym_dict = self.non_core_symmetry_groups['es']
                 
                for run_through in range(1):
                    indices = [i for i in range(len(self.allowed_molecules[self.current_state]['molecules']) - 1)]
                    chosen_structures = [self.allowed_molecules[self.current_state]['molecules'][idx] for idx in indices]
                    chosen_qm_energies = [self.allowed_molecules[self.current_state]['qm_energies'][idx] for idx in indices]
                    chosen_im_energies = [self.allowed_molecules[self.current_state]['im_energies'][idx] for idx in indices]
                    chosen_qm_gradients = [self.allowed_molecules[self.current_state]['qm_gradients'][idx] for idx in indices]
                    
                        
                    if self.use_opt_confidence_radius[1] == 'multi_grad':
                        
                        trust_radius = self.determine_trust_radius_gradient(chosen_structures, 
                                                                chosen_qm_energies,
                                                                chosen_qm_gradients,
                                                                chosen_im_energies, 
                                                                self.qm_data_point_dict[self.current_state], 
                                                                self.interpolation_settings[self.current_state],
                                                                self.qm_symmetry_datapoint_dict[self.current_state],
                                                                sym_dict,
                                                                exponent_p_q = (self.impes_drivers[self.current_state].exponent_p, self.impes_drivers[self.current_state].exponent_q))
                    
                    
                    elif self.use_opt_confidence_radius[1] == 'bayes':
                        trust_radius = self.determine_beysian_trust_radius(self.allowed_molecules[self.current_state]['molecules'], 
                                                                        self.allowed_molecules[self.current_state]['qm_energies'], 
                                                                        self.qm_data_point_dict[self.current_state], 
                                                                        self.interpolation_settings[self.current_state], 
                                                                        self.qm_symmetry_datapoint_dict[self.current_state],
                                                                        sym_dict)
                        

                    for idx, trust_radius in enumerate(trust_radius):
                        print(self.sorted_state_spec_im_labels[self.current_state][idx])
                        self.qm_data_point_dict[self.current_state][idx].update_confidence_radius(self.interpolation_settings[self.current_state]['imforcefield_file'], self.sorted_state_spec_im_labels[self.current_state][idx], trust_radius)
                        self.qm_data_point_dict[self.current_state][idx].confidence_radius = trust_radius

                
                self.add_a_point = False
        
        # calcualte energy gradient
        if self.add_a_point is False:
            self.previous_energy_list.append(current_state_difference[self.current_state][0])
            curr_state_diff = (current_state_difference[self.current_state][0])
            if curr_state_diff == 0:
                curr_state_diff = 1e-8
            if len(self.previous_energy_list) > 3:
                # Compute gradient (rate of change of energy difference)
                grad1 = self.previous_energy_list[-2] - self.previous_energy_list[-3]
                grad2 = self.previous_energy_list[-1] - self.previous_energy_list[-2]
                # Base skipping value calculation
                
                base_skip = min(round(abs(self.energy_threshold / (curr_state_diff)**2)), 20) - 1
                # Adjust skipping value based on gradient
                if grad2 > grad1:  # Energy difference is increasing
                    self.skipping_value = max(1, base_skip - 1)  # Reduce skipping for more frequent checks
                else:  # Energy difference is decreasing
                    self.skipping_value = base_skip + 10  # Increase skipping to check less often
                print(f"Energy Difference: {(current_state_difference[self.current_state][0]):.6f}, Gradient: {grad2 - grad1:.6f}, Skipping Value: {self.skipping_value}")
            else:
                self.skipping_value = min(round(abs(self.energy_threshold / (curr_state_diff)**2)), 20)
            print('len of molecules', len(self.allowed_molecules[self.current_state]['molecules']), self.use_opt_confidence_radius)
        
        self.skipping_value = 5
        if self.add_a_point and self.expansion:
            print('✨ A point is added! ✨', self.point_checker)
            print(molecule.get_xyz_string())
            ############# Implement constraint optimization ############
            
            state_specific_molecules = []
            
            
            if self.identfy_relevant_int_coordinates[0]:  
                for state_to_optim in addition_of_state_specific_points:                  
                    optimized_molecule = molecule
                    self.impes_drivers[state_to_optim].compute(molecule)
                    current_weights = self.impes_drivers[state_to_optim].weights
                    weights = [value for _, value in current_weights.items()]
                    used_labels = [label_idx for label_idx, _ in current_weights.items()]
                    # Sort labels and weights by descending weight
                    sorted_items = sorted(zip(used_labels, weights), key=lambda x: x[1], reverse=True)
                    total_weight = sum(weights)
                    cumulative_weight = 0.0
                    internal_coordinate_datapoints = []
                    for label, weight in sorted_items:
                        cumulative_weight += weight
                        internal_coordinate_datapoints.append(self.qm_data_point_dict[state_to_optim][label])
                        if cumulative_weight >= 0.8 * total_weight:
                            break
                    # qm_datapoints_weighted = [qm_datapoint for qm_datapoint in enumerate if ]
                    print('Items', sorted_items, len(internal_coordinate_datapoints), state_specific_energies[state_to_optim][0])
                    constraints = self.impes_drivers[state_to_optim].determine_important_internal_coordinates(state_specific_energies[state_to_optim][0], state_specific_gradients[state_to_optim][0], molecule, self.z_matrix, internal_coordinate_datapoints)
                    print('CONSTRAINTS', constraints)
                    if state_to_optim == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                        _, scf_tensors, rsp_results = self.compute_energy(self.drivers['gs'][0], molecule, current_basis)
     
                        opt_drv = OptimizationDriver(self.drivers['gs'][0])
                        opt_drv.ostream.mute()
                        
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
                        
                        for constraint in self.identfy_relevant_int_coordinates[1]:
                            if constraint in opt_constraint_list:
                                continue
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
                        optimized_molecule.set_charge(molecule.get_charge())
                        optimized_molecule.set_multiplicity(molecule.get_multiplicity())
                        current_basis = MolecularBasis.read(optimized_molecule, current_basis.get_main_basis_label())

                        # qm_energy, scf_tensors = self.compute_energy(drivers[0], optimized_molecule, opt_current_basis)
                        print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                    
                    elif state_to_optim >= 1 and isinstance(self.drivers['es'][0], LinearResponseEigenSolver) or state_to_optim >= 1 and isinstance(self.drivers['es'][0], TdaEigenSolver):
                            
                        self.drivers['es'][1].state_deriv_index = [opt_roots]
                        opt_drv = OptimizationDriver(self.drivers['es'][1])
                        current_basis = MolecularBasis.read(molecule, states_basis['es'])
                        _, _, rsp_results = self.compute_energy(self.drivers['es'][0], molecule, current_basis)
                        opt_drv.ostream.mute()
                        
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
                        
                        for constraint in self.identfy_relevant_int_coordinates[1]:
                            if constraint in opt_constraint_list:
                                continue
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
                        self.drivers['es'][3].ostream.mute()
                        self.drivers['es'][0].ostream.mute()
                        opt_results = opt_drv.compute(molecule, current_basis, self.drivers['es'][3], self.drivers['es'][0], rsp_results)
                        excitated_roots = [root for root in self.roots_to_follow if root != 0]
                        self.drivers['es'][1].state_deriv_index = excitated_roots
                        energy = opt_results['opt_energies'][-1]
                        optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                        print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                        molecule = optimized_molecule          
                    
                    
                    same = False
                    counter = 1
                    sum_of_values = 0
                    old_sorted_shorted = []
                    for item, vlaue in sorted_items:
                        sum_of_values += vlaue
                        old_sorted_shorted.append(item)
                        if sum_of_values > 0.8:       
                            break

                    old_list_change = []
                    while not same:
                        from typing import Iterable, Tuple, Set

                        def _norm_pair(pair: Tuple[int, int]) -> Tuple[int, int]:
                            """Order-insensitive bond key."""
                            a, b = pair
                            return (a, b) if a < b else (b, a)

                        def _to_rotset(rotatable_pairs: Iterable[Tuple[int, int]]) -> Set[Tuple[int, int]]:
                            """Normalize a list/set of rotatable (j,k) pairs to an order-insensitive set."""
                            return {_norm_pair(p) for p in rotatable_pairs}

                        def dihedral_center_pair(dihedral: Tuple[int, int, int, int]) -> Tuple[int, int]:
                            """Return the central bond (j,k) of a dihedral (i,j,k,l), order-insensitive."""
                            _, j, k, _ = dihedral
                            return _norm_pair((j, k))

                        def is_dihedral_rotatable(
                            dihedral: Tuple[int, int, int, int],
                            rotatable_pairs: Iterable[Tuple[int, int]]
                        ) -> bool:
                            """True if the dihedral's central bond is in the rotatable list."""
                            rotset = _to_rotset(rotatable_pairs)
                            return dihedral_center_pair(dihedral) in rotset

                        def filter_rotatable_dihedrals(
                            dihedrals: Iterable[Tuple[int, int, int, int]],
                            rotatable_pairs: Iterable[Tuple[int, int]]
                        ) -> list[Tuple[int, int, int, int]]:
                            """Keep only dihedrals whose central bond is rotatable."""
                            rotset = _to_rotset(rotatable_pairs)
                            return [d for d in dihedrals if dihedral_center_pair(d) in rotset]

                        interpolation_driver = InterpolationDriver(self.z_matrix)
                        interpolation_driver.update_settings(self.interpolation_settings[state_to_optim])
                        interpolation_driver.symmetry_information = self.impes_drivers[state_to_optim].symmetry_information
                        interpolation_driver.qm_symmetry_data_points = self.impes_drivers[state_to_optim].qm_symmetry_data_points
                        interpolation_driver.distance_thrsh = 1000
                        interpolation_driver.exponent_p = self.impes_drivers[state_to_optim].exponent_p
                        interpolation_driver.eq_bond_force_constants = self.eq_bond_force_constants
                        interpolation_driver.impes_coordinate.eq_bond_force_constants = self.eq_bond_force_constants
                        interpolation_driver.print = False
                        interpolation_driver.qm_data_points = self.impes_drivers[state_to_optim].qm_data_points
                        interpolation_driver.calc_optim_trust_radius = True

                        interpolation_driver.compute(optimized_molecule)

                        newly_weights = interpolation_driver.weights
                        new_weights = [value for _, value in newly_weights.items()]
                        new_used_labels = [label_idx for label_idx, _ in newly_weights.items()]
                        # Sort labels and weights by descending weight
                        new_sorted_items = sorted(zip(new_used_labels, new_weights), key=lambda x: x[1], reverse=True)
                    

                        new_order = [lbl for lbl, _ in new_sorted_items[:len(old_sorted_shorted)]]

                        same = (new_order == old_sorted_shorted)
                        # quick check: exact same ordering?
            
                        print(new_sorted_items, sorted_items, same)
                        int_diff = []
                        elem_list = []
                        if not same:
                            for element_idx, element in enumerate(self.z_matrix):

                                if len(element) == 4 and element not in constraints and is_dihedral_rotatable(element, self.all_rot_bonds):

                                    diff = np.sin(interpolation_driver.impes_coordinate.internal_coordinates_values[element_idx] - self.impes_drivers[state_to_optim].impes_coordinate.internal_coordinates_values[element_idx])
                                    if abs(diff) > 1e-1:
                                        elem_list.append(element)
                                        int_diff.append(diff)
                                
                            ordered_diif = sorted(zip(elem_list, int_diff), key=lambda x:x[1], reverse=True)
                            
                            print('Len ordered list', ordered_diif)
                            
                            if len(ordered_diif) > 0:
                                for spec_elem in range(counter):
                                    if spec_elem == len(ordered_diif):
                                        break
                                    
                                    dihedral_val_to_reset = molecule.get_dihedral_in_degrees((ordered_diif[spec_elem][0][0] + 1, ordered_diif[spec_elem][0][1] + 1, ordered_diif[spec_elem][0][2] + 1, ordered_diif[spec_elem][0][3] + 1))
                                    # print(dihedral_val_to_reset, ordered_diif[spec_elem][0])
                                    optimized_molecule.set_dihedral_in_degrees((ordered_diif[spec_elem][0][0] + 1, ordered_diif[spec_elem][0][1] + 1, ordered_diif[spec_elem][0][2] + 1, ordered_diif[spec_elem][0][3] + 1), dihedral_val_to_reset)
                                    
                                    print(optimized_molecule.get_dihedral_in_degrees((ordered_diif[spec_elem][0][0] + 1, ordered_diif[spec_elem][0][1] + 1, ordered_diif[spec_elem][0][2] + 1, ordered_diif[spec_elem][0][3] + 1)))      
                                    print('Orderd diff element', ordered_diif[spec_elem], old_list_change)
                                    if ordered_diif[spec_elem][0] in old_list_change:
                                        same = True
                            else:
                                same = True
                            
                            counter +=1
                            if len(old_list_change) == 0:
                                old_list_change = elem_list.copy()
                    state_specific_molecules.append((optimized_molecule, current_basis, [state_to_optim]))

                print('New optimized molecule \n', optimized_molecule.get_xyz_string())
                self.add_point(state_specific_molecules, self.non_core_symmetry_groups)
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0
                # self.point_adding_molecule[self.step] = (molecule, qm_energy, label)

            else:
                    
                current_basis = None
                if self.current_state == 0:
                    current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])
                else:
                    current_basis = MolecularBasis.read(molecule, self.basis_set_label['es'])
                state_specific_molecules.append((molecule, current_basis, addition_of_state_specific_points))
                    
                
                self.add_point(state_specific_molecules, self.non_core_symmetry_groups)
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0
                # self.point_adding_molecule[self.step] = (molecule, qm_energy, label)


    def add_point(self, state_specific_molecules, symmetry_information):
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
        category_label = None
 

        for entries in state_specific_molecules:
            symmetry_point = False
            if 0 in entries[2] and len(symmetry_information['gs']) != 0 and len(symmetry_information['gs'][2]) != 0:
                symmetry_mapping_groups = [item for item in range(len(entries[0].get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['gs'][1] for item in element]
                sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(symmetry_information['gs'][1], symmetry_information['gs'][5])
                
                # Generate all combinations
                keys = list(sym_dihedrals.keys())
                values = list(sym_dihedrals.values())
                combinations = list(itertools.product(*values))
                # Convert to list of dictionaries
                molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)
                print('MOlecule configs', molecule_configs)
                for i, molecule_config in enumerate(molecule_configs):
                    cur_molecule = Molecule.from_xyz_string(entries[0].get_xyz_string())
                    cur_molecule.set_charge(entries[0].get_charge())
                    cur_molecule.set_multiplicity(entries[0].get_multiplicity())
                    dihedral_to_change = []
                    for dihedral, angle in molecule_config.items():
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                    current_basis = MolecularBasis.read(cur_molecule, entries[1].get_main_basis_label())
                    if i > 0:
                        symmetry_point = True
                    adjusted_molecule['gs'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, symmetry_information['gs'], entries[2], symmetry_point))
            elif any(x > 0 for x in entries[2]) and len(symmetry_information['es']) != 0 and len(symmetry_information['es'][2]) != 0:
                symmetry_mapping_groups = [item for item in range(len(entries[0].get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['es'][1] for item in element]
                sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(symmetry_information['es'][1], symmetry_information['es'][5])
                
                # Generate all combinations
                keys = list(sym_dihedrals.keys())
                values = list(sym_dihedrals.values())
                combinations = list(itertools.product(*values))
                # Convert to list of dictionaries
                molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)
                for i, molecule_config in enumerate(molecule_configs):
                    cur_molecule = Molecule.from_xyz_string(entries[0].get_xyz_string())
                    cur_molecule.set_charge(entries[0].get_charge())
                    cur_molecule.set_multiplicity(entries[0].get_multiplicity())
                    dihedral_to_change = []
                    for dihedral, angle in molecule_config.items():
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                    if i > 0:
                        symmetry_point = True   
                    current_basis = MolecularBasis.read(cur_molecule, entries[1].get_main_basis_label())
                    adjusted_molecule['es'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, entries[2], symmetry_point))
            else:
                if 0 in entries[2]:
                    adjusted_molecule['gs'].append((entries[0], entries[1], 1, None, [0], symmetry_point)) 
                if any(x > 0 for x in entries[2]): 
                    states = [state for state in entries[2] if state > 0]
                    adjusted_molecule['es'].append((entries[0], entries[1], 1, None, states, symmetry_point)) 

        
        for key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue
            print(key, entries)
            drivers = self.drivers[key]
            
            org_roots = self.roots_to_follow[0]
            if any(root > 0 for root in self.roots_to_follow): 
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                    org_roots = drivers[1].state_deriv_index  

            label_counter = 0
            for mol_basis in entries:
                if any(root > 0 for root in self.roots_to_follow):
                    if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers, TdaEigenSolver):
                        root_to_follow_calc = mol_basis[4]           
                        drivers[1].state_deriv_index = root_to_follow_calc 
                
                energies, scf_results, rsp_results = self.compute_energy(drivers[0], mol_basis[0], mol_basis[1])
                
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    energies = energies[mol_basis[4]]

                print(energies)

                gradients = self.compute_gradient(drivers[1], mol_basis[0], mol_basis[1], scf_results, rsp_results)
                
                hessians = self.compute_hessian(drivers[2], mol_basis[0], mol_basis[1])
                masses = mol_basis[0].get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
                for number in range(len(energies)):
                    print(gradients[number], hessians[number])
                    label = f"point_{len(self.qm_data_point_dict[mol_basis[4][number]])}"
                    if not mol_basis[5]:
                        label = f"point_{len(self.qm_data_point_dict[mol_basis[4][number]]) +1}"
                    new_label = label

                    grad_vec = gradients[number].reshape(-1)         # (3N,)
                    hess_mat = hessians[number].reshape(grad_vec.size, grad_vec.size)
                    mw_grad_vec = inv_sqrt_masses * grad_vec
                    mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]
                    
                    impes_coordinate = InterpolationDatapoint(self.z_matrix)
                    impes_coordinate.update_settings(self.interpolation_settings[mol_basis[4][number]])
                    impes_coordinate.cartesian_coordinates = mol_basis[0].get_coordinates_in_bohr()
                    impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                    impes_coordinate.energy = energies[number]
                    impes_coordinate.gradient = mw_grad_vec.reshape(gradients[number].shape)
                    impes_coordinate.hessian = mw_hess_mat.reshape(hessians[number].shape)
                    impes_coordinate.transform_gradient_and_hessian()

                    
                    sampled_space = False
                    sampled_distances = []
                    impes_coordinate.confidence_radius = float(self.use_opt_confidence_radius[2])

                    if mol_basis[5]:
                            new_label = f'{label}_symmetry_{label_counter}'

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
                    print('Label counter', label_counter, len(self.qm_data_point_dict[mol_basis[4][number]]), new_label)
                    if not mol_basis[5]:
                        category_label = label
                        if impes_coordinate.use_inverse_bond_length:
                            category_label += "_rinv"
                        else:
                            category_label += "_r"
                        if impes_coordinate.use_cosine_dihedral:
                            category_label += "_cosine"
                        else:
                            category_label += "_dihedral"
                        print('dict and key', self.impes_drivers[mol_basis[4][number]].qm_symmetry_data_points, category_label)
                        self.qm_data_point_dict[mol_basis[4][number]].append(impes_coordinate)
                        self.sorted_state_spec_im_labels[mol_basis[4][number]].append(new_label)
                        self.qm_symmetry_datapoint_dict[mol_basis[4][number]][category_label] = [impes_coordinate]
                        self.impes_drivers[mol_basis[4][number]].qm_symmetry_data_points[category_label] = [impes_coordinate]
                        impes_coordinate.point_label = category_label
                        self.qm_energies_dict[mol_basis[4][number]].append(energies[number])
                        
                        self.impes_drivers[mol_basis[4][number]].qm_data_points = self.qm_data_point_dict[self.current_state]
                        self.impes_drivers[mol_basis[4][number]].labels = self.sorted_state_spec_im_labels[self.current_state]
                        print('I am writing the file', category_label, self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'])
                    
                        self.density_around_data_point[0][mol_basis[4][number]] += 1

                    else:
                        
                        self.qm_symmetry_datapoint_dict[mol_basis[4][number]][category_label].append(impes_coordinate)
                        self.impes_drivers[mol_basis[4][number]].qm_symmetry_data_points[category_label].append(impes_coordinate)

                    
                    impes_coordinate.confidence_radius = self.use_opt_confidence_radius[2]

                    print('data list', len(self.allowed_molecules[mol_basis[4][number]]['molecules']))
                    
                    impes_coordinate.write_hdf5(self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'], new_label)
                    if mol_basis[5]:
                        if self.use_opt_confidence_radius[0] and len(self.allowed_molecules[self.current_state]['molecules']) >= 20:
                            
                            trust_radius = None
                            sym_dict = self.non_core_symmetry_groups['gs']
                            if mol_basis[4][number] > 0 or root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
                                sym_dict = self.non_core_symmetry_groups['es']
                            
                            
                            for run_through in range(1):
                                indices = [i for i in range(len(self.allowed_molecules[mol_basis[4][number]]['molecules']) - 1)]
                                chosen_structures = [self.allowed_molecules[mol_basis[4][number]]['molecules'][idx] for idx in indices]
                                chosen_qm_energies = [self.allowed_molecules[mol_basis[4][number]]['qm_energies'][idx] for idx in indices]
                                chosen_im_energies = [self.allowed_molecules[mol_basis[4][number]]['im_energies'][idx] for idx in indices]
                                chosen_qm_gradients = [self.allowed_molecules[mol_basis[4][number]]['qm_gradients'][idx] for idx in indices]
                                
                                    
                                if self.use_opt_confidence_radius[1] == 'multi_grad':
                                    
                                    trust_radius = self.determine_trust_radius_gradient(chosen_structures, 
                                                                            chosen_qm_energies,
                                                                            chosen_qm_gradients,
                                                                            chosen_im_energies, 
                                                                            self.qm_data_point_dict[mol_basis[4][number]], 
                                                                            self.interpolation_settings[mol_basis[4][number]],
                                                                            self.qm_symmetry_datapoint_dict[mol_basis[4][number]],
                                                                            sym_dict,
                                                                            exponent_p_q = (self.impes_drivers[mol_basis[4][number]].exponent_p, self.impes_drivers[mol_basis[4][number]].exponent_q))
                                
                                
                                elif self.use_opt_confidence_radius[1] == 'bayes':
                                    trust_radius = self.determine_beysian_trust_radius(self.allowed_molecules[mol_basis[4][number]]['molecules'], 
                                                                                    self.allowed_molecules[mol_basis[4][number]]['qm_energies'], 
                                                                                    self.qm_data_point_dict[mol_basis[4][number]], 
                                                                                    self.interpolation_settings[mol_basis[4][number]], 
                                                                                    self.qm_symmetry_datapoint_dict[mol_basis[4][number]],
                                                                                    sym_dict)
                            
                                for idx, trust_radius in enumerate(trust_radius):
                                    print(self.sorted_state_spec_im_labels[mol_basis[4][number]][idx])
                                    self.qm_data_point_dict[mol_basis[4][number]][idx].update_confidence_radius(self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'], self.sorted_state_spec_im_labels[mol_basis[4][number]][idx], trust_radius)
                                    self.qm_data_point_dict[mol_basis[4][number]][idx].confidence_radius = trust_radius
                                self.simulation.saveCheckpoint('checkpoint')
                                
                                # self.output_file_writer(self.summary_output)
                                state_spec_dict = {'pot_energies':[], 'gradients':[]}
                                self.gloabal_sim_informations = {f'state_{root}':state_spec_dict for root in self.roots_to_follow}
                                self.gloabal_sim_informations['coordinates_ang'] = []
                                self.gloabal_sim_informations['state'] = []
                                self.gloabal_sim_informations['temperatures'] = []
                    # exit()
                    for root in mol_basis[4]:
                        self.write_qm_energy_determined_points(self.allowed_molecules[root]['molecules'][self.last_added: ],
                                                            self.allowed_molecules[root]['qm_energies'][self.last_added:],
                                                            self.allowed_molecules[root]['qm_gradients'][self.last_added:],
                                                            self.allowed_molecules[root]['im_energies'][self.last_added:],
                                                            [],
                                                            root)
                        self.last_added = len(self.allowed_molecules[root]['molecules'])
                label_counter += 1
            if any(root > 0 for root in self.roots_to_follow):
                                   
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                    
                    drivers[1].state_deriv_index = org_roots
    
    
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
                interpolation_driver.print = False
                interpolation_driver.qm_data_points = combined_datapoints
                
                interpolation_driver.compute(mol)
                new_im_energy = interpolation_driver.get_energy()
                diff = (new_im_energy - qm_energies[i]) * hartree_in_kcalpermol()
                sum_sq_error += (diff)**2 / (0.1**2 * distance**6)

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

    
    def determine_trust_radius_gradient(self, molecules, qm_energies, qm_gradients, im_energies, datapoints, interpolation_setting, sym_datapoints, sym_dict, exponent_p_q):
        
        def optimize_trust_radius(alphas, geom_list, E_ref_list, G_ref_list, E_im_list, dps, impes_dict, sym_datapoints, sym_dict, exponent_p_q):
            """
            Perform the gradient-based optimization to find R*
            that minimizes the sum of squared errors to reference QM energies.
            """

            bounds = [(0.01, 1.5)] * len(dps)
             
            opt = AlphaOptimizer(self.z_matrix, impes_dict, sym_dict, sym_datapoints, dps,
                 geom_list, E_ref_list, G_ref_list, exponent_p_q,
                 e_x=self.use_opt_confidence_radius[2],
                 beta=0.8, n_workers=os.cpu_count())  # pick sensible n_workers

            minimizer_kwargs = {"method": "L-BFGS-B", "jac": opt.jac, "bounds": bounds, "options": {"disp": True, "gtol": 1e-4, "ftol": 1e-9, "maxls": 10}}
            res = basinhopping(opt.fun, x0=alphas, minimizer_kwargs=minimizer_kwargs, niter=10)

            return res
        
        inital_alphas = [dp.confidence_radius for dp in datapoints]

        print('INPUT Trust radius', inital_alphas)

        additional_molecules = []
        additional_energies = []
        additional_gradients = []
        for dp in datapoints:
            additional_molecules.append(Molecule(molecules[0].get_labels(), dp.cartesian_coordinates, 'bohr'))
            additional_energies.append(dp.energy)
            additional_gradients.append(dp.gradient)
        molecules.extend(additional_molecules)
        qm_energies.extend(additional_energies)
        qm_gradients.extend(additional_gradients)
        trust_radius = optimize_trust_radius(inital_alphas, molecules[:], qm_energies[:], qm_gradients[:], im_energies[:], datapoints, interpolation_setting, sym_datapoints, sym_dict, exponent_p_q)

        print('FINAL Trust radius', trust_radius['x'])

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
        scf_results = None
        rsp_results = None

        # XTB
        if isinstance(qm_driver, XtbDriver):
            qm_driver.ostream.mute()
            qm_driver.compute(molecule)
            qm_energy = qm_driver.get_energy()
            qm_energy = np.array([qm_energy])

        # restricted SCF
        elif isinstance(qm_driver, ScfRestrictedDriver):
            qm_driver.ostream.mute()
            scf_results = qm_driver.compute(molecule, basis)
            qm_energy = qm_driver.scf_energy
            qm_energy = np.array([qm_energy])
            qm_driver.ostream.unmute()
            qm_driver.filename = None
            qm_driver.checkpoint_file = None
        
        elif isinstance(qm_driver, LinearResponseEigenSolver) or isinstance(qm_driver, TdaEigenSolver):
            self.drivers['es'][0].ostream.mute()
            scf_results = self.drivers['es'][3].compute(molecule, basis)
            self.drivers['es'][3].filename = None
            self.drivers['es'][3].checkpoint_file = None
            scf_energy = self.drivers['es'][3].scf_energy
            qm_driver.ostream.mute()
            rsp_results = qm_driver.compute(molecule, basis, scf_results)
            
            qm_energy = np.insert(rsp_results['eigenvalues'] + scf_energy, 0, scf_energy)

        if qm_energy is None:
            error_txt = "Could not compute the QM energy. "
            error_txt += "Please define a QM driver."
            raise ValueError(error_txt)

        return qm_energy, scf_results, rsp_results

    def compute_gradient(self, grad_driver, molecule, basis=None, scf_results=None, rsp_results=None):
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

        elif isinstance(grad_driver, TddftGradientDriver):
            grad_driver.ostream.mute()
            grad_driver._rsp_results = None
            grad_driver.compute(molecule, basis, self.drivers['es'][3], self.drivers['es'][0], rsp_results)
            qm_gradient = grad_driver.gradient

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

        qm_hessians = None

        if isinstance(hess_driver, XtbHessianDriver):
            hess_driver.ostream.mute()
            hess_driver.compute(molecule)
            qm_hessian = hess_driver.hessian
            hess_driver.ostream.unmute()
            qm_hessians = np.array([qm_hessian])

        elif isinstance(hess_driver, ScfHessianDriver):
            hess_driver.ostream.mute()
            hess_driver.compute(molecule, basis)
            qm_hessian = hess_driver.hessian
            qm_hessians = np.array([qm_hessian])
            hess_driver.ostream.unmute()
        
        elif isinstance(hess_driver, TddftHessianDriver):
            roots = self.drivers['es'][1].state_deriv_index

            qm_hessians = []
            for root in roots:
                print('In hessian clac', root, hess_driver)
                self.drivers['es'][1].state_deriv_index = root
                hess_driver.compute(molecule, basis)
                qm_hessian = hess_driver.hessian
                qm_hessians.append(qm_hessian)

        if qm_hessians is None:
            error_txt = "Could not compute the QM Hessian. "
            error_txt += "Please define a QM Hessian driver."
            raise ValueError(error_txt)


        return qm_hessians



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
    
    def write_qm_energy_determined_points(self, molecules, qm_energies, qm_gradients, im_energies, distances, state):           
        
        mode = 'a' if os.path.exists(self.reference_struc_energies_file) and os.path.getsize(self.reference_struc_energies_file) > 0 else 'w'
        print(f"Opening file in mode: {mode}")

        with open(self.reference_struc_energies_file, mode) as file:
            
            for i, dyn_mol in enumerate(molecules):
                current_xyz_string = dyn_mol.get_xyz_string()
                xyz_lines = current_xyz_string.splitlines()

                # Add energy info to comment line (second line)
                xyz_lines[1] += f' Energies  QM: {qm_energies[i]} QM_G: {qm_gradients[i]} IM: {im_energies[i]} Distance: {1.0}  State: {state}'
                updated_xyz_string = "\n".join(xyz_lines)

                file.write(f"{updated_xyz_string}\n\n")

    def extract_reference_structures(self, filename, roots):

        def parse_gradient(grad_str):
            # Remove outer [[ and ]]
            grad_str = grad_str.strip()
            grad_str = grad_str.lstrip("[").rstrip("]")  # remove only one layer
            # Split by rows: each line with numbers
            rows = grad_str.strip().split("\n")
            data = []
            for row in rows:
                # Remove any inner brackets and extra spaces
                row_clean = row.replace("[", "").replace("]", "").strip()
                if row_clean:  # skip empty lines
                    data.append([float(x) for x in row_clean.split()])
            return np.array(data)

        xyz_strings = []  
        qm_energies = []
        qm_gradients = []  
        im_energies = []
        distances = []
        states = []  
        
        reference_molecules_check = {root: {'molecules': [], 'qm_energies': [], 'qm_gradients': [], 'im_energies': [], 'distances': []} for root in roots}
    
        with open(filename, 'r') as file:
            lines = file.readlines()

        atom_line_pattern = re.compile(r"^[A-Z][a-z]?\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+")

        current_xyz = []
        skip_gradient = False

        for line in lines:
            line = line.strip()

            # new molecule start (atom count)
            if line.isdigit():  
                if current_xyz: 
                    xyz_strings.append('\n'.join(current_xyz))
                    current_xyz = []  
                current_xyz.append(line)
                current_xyz.append(' ')  # Add an empty line for the comment

            elif "Energies" in line:
                skip_gradient = True

            elif skip_gradient:
                # Check if this line is an atom line → resume adding
                if atom_line_pattern.match(line):
                    skip_gradient = False
                    current_xyz.append(line)
                # else: we are still in QM_G gradient → skip it
            else:
                current_xyz.append(line)

        content = "\n".join(lines)
        pattern = (
                    r"QM:\s*([-0-9\.]+)\s*"
                    r"QM_G:\s*(\[\[[\s\S]*?\]\](?=\s*IM:))\s*"  # stop right before IM:
                    r"IM:\s*([-0-9\.]+)\s*"
                    r"Distance:\s*([-0-9\.]+)\s*"
                    r"State:\s*([0-9]+)"
                )
        matches = list(re.finditer(pattern, content, re.DOTALL))
        for match in matches:
            qm_energy = float(match.group(1))
            grad = parse_gradient(match.group(2))
            im_energy = float(match.group(3))
            distance = float(match.group(4))
            state = int(match.group(5))


            qm_energies.append(qm_energy)
            im_energies.append(im_energy)
            distances.append(distance)
            states.append(state)
            qm_gradients.append(grad)

        if current_xyz:
            xyz_strings.append('\n'.join(current_xyz))

        for i, state in enumerate(states):
            
            mol = Molecule.from_xyz_string(xyz_strings[i])

            reference_molecules_check[state]['molecules'].append(mol)
            reference_molecules_check[state]['qm_energies'].append(qm_energies[i])
            reference_molecules_check[state]['qm_gradients'].append(qm_gradients[i])
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
        