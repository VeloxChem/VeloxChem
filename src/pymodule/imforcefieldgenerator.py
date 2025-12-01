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

import numpy as np
import math
import itertools
import os
import random
from contextlib import redirect_stderr
from io import StringIO
from copy import deepcopy
from itertools import combinations


from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .tdaeigensolver import TdaEigenSolver
from .lreigensolver import LinearResponseEigenSolver
from .tddftgradientdriver import TddftGradientDriver
from .tddfthessiandriver import TddftHessianDriver

from .molecularbasis import MolecularBasis
from .openmmdynamics import OpenMMDynamics
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .imdatabasepointcollecter import IMDatabasePointCollecter
# from .impesdatabasebuilder import ImpesDatabaseBuilder
# from .imdatabasedriver import IMDatabaseDriver
#from .impesforcefieldgenerator_parallel import ImpesForceFieldGeneratorParallel
from .mmforcefieldgenerator import MMForceFieldGenerator
from .conformergenerator import ConformerGenerator
from .optimizationdriver import OptimizationDriver
from .atommapper import AtomMapper
from .tsguesser import TransitionStateGuesser

import openmm as mm
import openmm.app as app
import openmm.unit as unit
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
                    Typically initialized as `im_database.h5`.

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

        - qm_energies: A list to store QM energy from individual simulations or calculations (kj/mol).

        - total_energies: A list to store total energy (kj/mol).

        - molecules: A list to store molecular structures sampled during simulations (kj/mol).

        - kinetic_energies: A list to store kinetic energy from simulations (kj/mol).

        - point_added_molecules: A list of molecules for which new data points were added to the database.

        - unique_molecules: A list of unique molecular structures identified during database construction.

        - dynamics_method: Determines the method to generate molecular structures for the database quality conformation.
        
        - nstruc_to_confirm_database_quality: Number of randomly selected strucutures for the database quality check.
    """
    
    def __init__(self, ground_state_driver=None, excited_state_driver=None, roots_to_follow=[0], rsp_method='tda'):

        self.density_of_datapoints = None
        self.qm_data_points = None
        self.qmlabels = None
        self.molecules_along_rp = None
        self.conformal_structures = None
        self.sampling_structures = 1
        
        self.excitation_pulse = None
        self.state_specific_molecules = None
        self.datafile = None
        self.dihedrals_dict = None
        self.allowed_deviation = None
        self.z_matrix = None
        self.int_coord_bond_information = None
        self.symmetry_information = None
        self.symmetry_dihedral_lists = {}
        self.all_rotatable_bonds = None
        self.gs_basis_set_label = 'def2-svp'
        self.es_basis_set_label = '6-31g*'

        self.drivers = {'gs': None, 'es':None}

        if isinstance(ground_state_driver, ScfRestrictedDriver):
        # should be necessary to initialize
       
            qm_grad_driver = ScfGradientDriver(ground_state_driver)
            qm_hess_driver = ScfHessianDriver(ground_state_driver)
            
            self.drivers['gs'] = (ground_state_driver, qm_grad_driver, qm_hess_driver)
        ##########################################################
        ################# External Settings ######################
        ##########################################################

        if isinstance(excited_state_driver, ScfRestrictedDriver):
            
            excited_state_response_driver = TdaEigenSolver()
            if rsp_method == 'tddft':
                excited_state_response_driver = LinearResponseEigenSolver()
            excited_state_response_driver.nstates = 10
            excited_state_grad_driver = TddftGradientDriver(excited_state_driver)
            excitated_roots = [root for root in roots_to_follow if root != 0]
            excited_state_grad_driver.state_deriv_index = excitated_roots
            excited_state_hess_driver = TddftHessianDriver(excited_state_driver, excited_state_response_driver, excited_state_grad_driver)

            self.drivers['es'] = (excited_state_response_driver, excited_state_grad_driver, excited_state_hess_driver, excited_state_driver)


        if isinstance(ground_state_driver, XtbDriver):
            qm_grad_driver = XtbGradientDriver(ground_state_driver)
            qm_hess_driver = XtbHessianDriver(ground_state_driver)
            self.drivers['gs'] = (ground_state_driver, qm_grad_driver, qm_hess_driver)
        
        self.states_interpolation_settings = {root: None for root in roots_to_follow}
        self.states_data_point_density = {root: None for root in roots_to_follow}
        self.roots_to_follow = roots_to_follow
        ##########################################################
        # variables for the interpolation
        self.interpolation_settings = None
        self.interpolation_type = 'shepard'
        self.weightfunction_type = 'cartesian'
        self.exponent_p = 2
        self.exponent_q = 2
        self.confidence_radius = 0.5
        self.imforcefieldfiles = None
        self.use_inverse_bond_length = True
        self.use_cosine_dihedral = False

        self.reaction_coordinates = None
        self.optim_constraints = None
        self.reaction_molecules_dict = None

        # variables for the forcefield generation and database expansion
        self.dynamics_settings = None
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
        self.reference_struc_energy_file = None   
        self.desired_point_density = 50
        self.converged_cycle = 5
        self.energy_threshold = 1.5
        self.gradient_rmsd_thrsh = 1.0
        self.force_orient_thrsh = 0.93
        self.distance_thrsh = 0.1
        self.start_collect = 0
        self.solvent = 'gas'
        self.add_bias_force = None
        self.bias_force_reaction_idx = None
        self.bias_force_reaction_prop = None

        self.atom_transfer_reaction_path = None

        # individual run
        self.qm_energies = []
        self.total_energies = []
        self.molecules = None
        self.kinetic_energies = []
        self.point_added_molecules = []
        self.unique_molecules = []

        # In here I want to store Number_of_dp, exponent_p, exponent_q
        self.im_results = {'n_datapoints': None, 'RMSD': None, '|D|':None} 


        # confirm database quality
        self.dynamics_method = 'MM'
        self.nstruc_to_confirm_database_quality = 50

        
        # set boolean for the optimization features of the metho
        self.use_minimized_structures = [True, []],
        self.add_conformal_structures = True
        self.add_structures_along_rcs = False
        self.use_symmetry = True
        self.identfy_relevant_int_coordinates = True
        self.add_gpr_model = True
        self.use_opt_confidence_radius = [False, 'single', 0.5, 0.3]
    
        self.eq_bond_force_constants = None


    def set_up_the_system(self, molecule, extract_z_matrix=None):

        """
        Assign the neccessary variables with respected values. 

        :param molecule: original molecule

        :param target_dihedrals: is a list of dihedrals that should be scanned during the dynamics

        :param sampling_structures: devides the searchspace around given rotatbale dihedrals
            
        """
        
        
        def determine_dimer_fragments(molecule):
            

            natms = len(molecule.get_labels())
            con_matrix = molecule.get_connectivity_matrix()

            dimers = []
            checked_indices = []

            for i in range(natms):
                cont_i = False
                for dimer in dimers:
                    if i in dimer:
                        cont_i = True
                        break
                if cont_i:
                    continue  
                if i not in checked_indices:
                    checked_indices.append(i)      
                expanded = False
                for counter, conn_j in enumerate(checked_indices):

                    for row_idx in range(natms):

                        connectivity = con_matrix[conn_j, row_idx]

                        if row_idx not in checked_indices and connectivity == 1:
                            
                            checked_indices.append(row_idx)
                            expanded = True
                    

                    
                if not expanded:
                
                    dimers.append(checked_indices)
                    checked_indices = []

            bond_indices = []

            # go through all unique pairs of dimers
            for d1, d2 in combinations(dimers, 2):
                # make all cross-pairs between them
                bond_indices.extend((i, j) for i in d1 for j in d2)


            return bond_indices


        def regroup_by_rotatable_connection(molecule, groups, rotatable_bonds, conn):
            new_groups = {'gs': [], 'es': [], 'non_rotatable': []}
            rot_groups = {'gs': [], 'es': []}
            labels = molecule.get_labels()

            def determine_state(a1, a2):
                neighbors_a1 = sum(conn[a1])
                neighbors_a2 = sum(conn[a2])
                element_a1 = labels[a1][0]
                element_a2 = labels[a2][0]

                # NH2 rule
                if (element_a1 == 'N' and neighbors_a1 == 2) or (element_a2 == 'N' and neighbors_a2 == 2):
                    return 'gs'

                # Oxygen + sp2 carbon → es
                if (element_a1 == 'O' and element_a2 == 'C' and neighbors_a2 == 3) or \
                (element_a2 == 'O' and element_a1 == 'C' and neighbors_a1 == 3):
                    return 'es'

                # sp2-sp2 → es
                if neighbors_a1 == 3 and neighbors_a2 == 3:
                    return 'es'

                # default → gs
                return 'gs'

            for group in groups:
                if len(group) > 3:
                    continue
                connected_subgroups = {}  # key: (state, atom in bond), value: atoms in group connected to it
                connected_rotatable_bonds = set()
                for atom in group:
                    for a1, a2 in rotatable_bonds:
                        state = determine_state(a1, a2)

                        if conn[atom, a1] or conn[atom, a2]:
                            connected_rotatable_bonds.add((a1, a2))

                        if conn[atom, a1]:
                            connected_subgroups.setdefault((state, a1), []).append(atom)
                        if conn[atom, a2]:
                            connected_subgroups.setdefault((state, a2), []).append(atom)
                
                if len(connected_rotatable_bonds) > 1:
                    new_groups['non_rotatable'].append(group)
                    continue

                if not connected_subgroups:
                    new_groups['non_rotatable'].append(group)
                else:
                    for (state, _), subgroup in connected_subgroups.items():
                        subgroup = list(set(subgroup))
                        if len(subgroup) > 1 and subgroup not in rot_groups[state]:
                            for atom in subgroup:
                                print(molecule.get_labels()[atom], sum(conn[atom]))
                            rot_groups[state].append(sorted(subgroup))
                            new_groups[state].append(sorted(subgroup))
  
            return new_groups, rot_groups


        if self.z_matrix is None and extract_z_matrix is None:
            dimer_coordinates = determine_dimer_fragments(molecule)
            if self.reaction_coordinates is not None:
                self.reaction_coordinates.extend(dimer_coordinates)

            elif len(dimer_coordinates) > 0:
                self.reaction_coordinates = dimer_coordinates

            self.z_matrix = self.define_z_matrix(molecule, self.reaction_coordinates)
        elif self.z_matrix is None:
            int_driver = InterpolationDriver()
            int_driver.update_settings({ 'interpolation_type':self.interpolation_type,
                                    'weightfunction_type':self.weightfunction_type,
                                    'exponent_p':self.exponent_p,
                                    'exponent_q':self.exponent_q, 
                                    'confidence_radius':self.confidence_radius,
                                    'imforcefield_file':self.imforcefieldfiles[self.roots_to_follow[0]],
                                    'use_inverse_bond_length':self.use_inverse_bond_length,
                                    'use_cosine_dihedral':self.use_cosine_dihedral
                                })

            _, self.z_matrix = int_driver.read_labels()

        angle_index = next((i for i, x in enumerate(self.z_matrix) if len(x) == 3), len(self.z_matrix))
        dihedral_index = next((i for i, x in enumerate(self.z_matrix) if len(x) == 4), len(self.z_matrix))

        self.qm_data_points = None
        self.molecule = molecule

        atom_mapper = AtomMapper(molecule, molecule)
        symmetry_groups = atom_mapper.determine_symmetry_group()

        if not self.use_symmetry:
            symmetry_groups = (symmetry_groups[0], [], symmetry_groups[2])
        ff_gen = MMForceFieldGenerator()
        ff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
        ff_gen.create_topology(molecule)

        rotatable_bonds = deepcopy(ff_gen.rotatable_bonds)
        es_rotatable_bonds = deepcopy(ff_gen.excited_states_rot_bond)

        rotatable_bonds.extend(es_rotatable_bonds)

        rotatable_bonds_zero_based = [(i - 1, j - 1) for (i, j) in rotatable_bonds]
        self.all_rotatable_bonds = rotatable_bonds_zero_based
        all_exclision = [element for rot_bond in rotatable_bonds_zero_based for element in rot_bond]

        symmetry_groups_ref = [groups for groups in symmetry_groups[1] if not any(item in all_exclision for item in groups)]

        # reduce the symmetry to only CH3 or CH2 symmetry groups for the time being
        regrouped, rot_groups = regroup_by_rotatable_connection(molecule, symmetry_groups_ref, rotatable_bonds_zero_based, molecule.get_connectivity_matrix())
        

        self.symmetry_information = {'gs': (), 'es': ()}
        for root in self.roots_to_follow:
            if root == 0 or root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
        
                non_core_atoms = [element for group in regrouped['gs'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                angles_to_set, _, _, self.symmetry_dihedral_lists = self.adjust_symmetry_dihedrals(molecule, rot_groups['gs'], rotatable_bonds_zero_based)
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.z_matrix[dihedral_index:], start=dihedral_index):

                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)

                self.symmetry_information['gs'] = (symmetry_groups[0], rot_groups['gs'], regrouped['gs'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [angle_index, dihedral_index])


            if root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip is False:

                non_core_atoms = [element for group in regrouped['gs'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                angles_to_set, _, _, self.symmetry_dihedral_lists = self.adjust_symmetry_dihedrals(molecule, rot_groups['gs'], rotatable_bonds_zero_based)
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.z_matrix[dihedral_index:], start=dihedral_index):

                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)

                self.symmetry_information['es'] = (symmetry_groups[0], rot_groups['es'], regrouped['es'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [angle_index, dihedral_index])

            if root >= 2 and len(self.symmetry_information['es']) == 0 and self.drivers['es']:
                non_core_atoms = [element for group in regrouped['es'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                print(rot_groups['es'], rotatable_bonds_zero_based)
                angles_to_set, _, _, self.symmetry_dihedral_lists = self.adjust_symmetry_dihedrals(molecule, rot_groups['es'], rotatable_bonds_zero_based)
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.z_matrix[dihedral_index:], start=dihedral_index):

                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)

                self.symmetry_information['es'] = (symmetry_groups[0], rot_groups['es'], regrouped['es'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [angle_index, dihedral_index])


        if self.reaction_molecules_dict is not None and extract_z_matrix is None:
            for entry in self.reaction_molecules_dict:
                if entry['root'] not in self.atom_transfer_reaction_path:

                    self.atom_transfer_reaction_path = {entry['root']: []}
                self.atom_transfer_reaction_path[entry['root']].append(self.determine_atom_transfer_reaction_path(entry['reactants'], entry['products']))


        if self.add_conformal_structures:

            rotatable_dihedrals_dict = {}
        
            conformer_generator = ConformerGenerator()
            conformal_structures = conformer_generator.generate(molecule)
            dihedral_canditates = conformer_generator.dihedral_candidates

            
            conf_molecule_xyz = conformal_structures['molecules'][0].get_xyz_string()

            for entry_idx, entries in enumerate(dihedral_canditates):
                
                dihedral = entries[0]
                periodicity = len(entries[1])
                conformal_structures['molecules'][entry_idx] = (conformal_structures['molecules'][entry_idx], [dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])

                for i in range(len(entries[1])):
                    current_molecule = Molecule.from_xyz_string(conf_molecule_xyz)
                    if i + 1 < len(entries[1]):
                        new_angle = entries[1][i] - (entries[1][i + 1]/periodicity)
                        current_molecule.set_dihedral_in_degrees([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], new_angle)
                        constriant = [dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1]
                        conformal_structures['molecules'].append((current_molecule, constriant))
                        print(new_angle)
                    else:
                        new_angle = entries[1][i] - (entries[1][i]/2)
                        current_molecule.set_dihedral_in_degrees([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], new_angle)
                        constriant = [dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1]
                        conformal_structures['molecules'].append((current_molecule, constriant))
                        print(new_angle)
            
            self.conformal_structures = {None : conformal_structures['molecules']}

    def compute(self, molecule, states_basis=None):

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
        states_basis = {'gs':self.gs_basis_set_label, 'es':self.es_basis_set_label}
        root_extract_z_matrix = None
        if self.imforcefieldfiles is not None:
            root_extract_z_matrix = {}
            for i, root in enumerate(self.roots_to_follow):
                if root not in self.imforcefieldfiles:
                    self.imforcefieldfiles[self.roots_to_follow[i]] = f'im_database_{root}.h5'
                else:
                    root_extract_z_matrix[self.roots_to_follow[i]] = True

        else:
            self.imforcefieldfiles = {}
            standard_files = [f'im_database_{root}.h5' for root in self.roots_to_follow]
            for i, file in enumerate(standard_files):
                
                if os.path.exists(file):
                    if root_extract_z_matrix is None:
                        root_extract_z_matrix = {}
                        root_extract_z_matrix[i] = True               
                self.imforcefieldfiles[self.roots_to_follow[i]] = file
        
        print(f'IMPORTANT: IM ForceFieldFile is initalized from the current directory as {file}')
        self.set_up_the_system(molecule, extract_z_matrix=root_extract_z_matrix)

        print('Set up the system is here')
        if len(self.roots_to_follow) > 1:
            print('The molecule that is used for the database construction is minimized always based on the ground state Energy -> to find the local minimum based on the '
            'potential energy surface of the ground state!')

            for root in self.roots_to_follow:
                    
                imforcefieldfile = self.imforcefieldfiles[root]
                self.states_interpolation_settings[root] = { 'interpolation_type':self.interpolation_type,
                                    'weightfunction_type':self.weightfunction_type,
                                    'exponent_p':self.exponent_p,
                                    'exponent_q':self.exponent_q, 
                                    'confidence_radius':self.confidence_radius,
                                    'imforcefield_file':imforcefieldfile,
                                    'use_inverse_bond_length':self.use_inverse_bond_length,
                                    'use_cosine_dihedral':self.use_cosine_dihedral
                                }
                
            self.dynamics_settings = {  'drivers':self.drivers,
                                        'basis_set_label': states_basis,
                                        'duration':self.duration, 'temperature':self.temperature, 'solvent':self.solvent,
                                        'pressure':self.force_constant, 'force_constant': self.force_constant, 'ensemble':self.ensemble,
                                        'timestep': self.timestep, 'nsteps': self.nsteps, 'friction':self.friction,
                                        'snapshots':self.snapshots, 'trajectory_file':self.trajectory_file, 'reference_struc_energy_file':self.reference_struc_energy_file,
                                        'desired_datapoint_density':self.desired_point_density, 'converged_cycle': self.converged_cycle, 
                                        'energy_threshold':self.energy_threshold, 'grad_rmsd_thrsh': self.gradient_rmsd_thrsh, 'force_orient_thrsh':self.force_orient_thrsh,
                                        'NAC':False, 'load_system': None, 'collect_qm_points_from':self.start_collect, 'roots_to_follow':self.roots_to_follow, 'excitation_pulse':self.excitation_pulse}
            
            files_to_add_conf = []
            molecules_to_add_info = []
            for root in self.roots_to_follow:

                if not os.path.exists(self.imforcefieldfiles[root]):
                    files_to_add_conf.append(root)
            
            if self.atom_transfer_reaction_path is not None and len(files_to_add_conf) != 0:

                for ts_root, reaction_mol in self.atom_transfer_reaction_path.items():
                    
                    if self.use_minimized_structures[0]:       
                        
                        if ts_root == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                        
                            opt_drv = OptimizationDriver(opt_qm_driver)
                            current_basis = MolecularBasis.read(reaction_mol, states_basis['gs'])
                            _, scf_results, _ = self.compute_energy(opt_qm_driver, reaction_mol, current_basis)
                            opt_drv.ostream.mute()
                            
                            opt_constraint_list = []
                            for constraint in self.use_minimized_structures[1]:
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

                            opt_results = opt_drv.compute(reaction_mol, current_basis, scf_results)
                            optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                            reaction_mol = optimized_molecule
                            print(optimized_molecule.get_xyz_string())

                        
                    if ts_root == 0:
                        current_basis = MolecularBasis.read(reaction_mol, states_basis['gs'])
                        molecules_to_add_info.append((reaction_mol, current_basis, [ts_root]))
                    else:
                        current_basis = MolecularBasis.read(reaction_mol, states_basis['es'])
                        molecules_to_add_info.append((reaction_mol, current_basis, [ts_root]))
                self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)
            
            if self.add_conformal_structures and len(files_to_add_conf) != 0:
                for counter, entry in enumerate(self.conformal_structures.items()):
                    key_old, molecules = entry
                    for i, mol in enumerate(molecules):
                    
                        key = key_old
                        for opt_root in self.roots_to_follow:
                            if self.use_minimized_structures[0]:
                                if opt_roots not in files_to_add_conf or opt_root not in self.use_minimized_structures[2]:
                                    if opt_root == 0:
                                        current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                                        molecules_to_add_info.append((molecule, current_basis, [opt_root]))
                                    else:
                                        current_basis = MolecularBasis.read(molecule, states_basis['es'])
                                        molecules_to_add_info.append((molecule, current_basis, [opt_root]))
                                    continue
                                root_to_add = opt_root
                                if  root_to_add == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                                
                                    opt_drv = OptimizationDriver(self.drivers['gs'][0])
                                    current_basis = MolecularBasis.read(mol, states_basis['gs'])
                                    _, scf_results, _ = self.compute_energy(self.drivers['gs'][0], mol, current_basis)
                                    opt_drv.ostream.mute()
                                    
                                    opt_constraint_list = []
                                    for constraint in self.use_minimized_structures[1]:
                                        if len(constraint) == 2:
                                            opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                                            opt_constraint_list.append(opt_constraint)
                                        
                                        elif len(constraint) == 3:
                                            opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                                            opt_constraint_list.append(opt_constraint)
                                    
                                        else:
                                            opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                                            opt_constraint_list.append(opt_constraint)
                                    opt_constraint_list.extend(constraint_conf)
                                    opt_drv.constraints = opt_constraint_list
                                    
                                    opt_results = opt_drv.compute(mol, current_basis, scf_results)
                                    energy = opt_results['opt_energies'][-1]
                                    optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])

                                    current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                                    molecules_to_add_info.append((optimized_molecule, current_basis, [root_to_add]))

                                    print(optimized_molecule.get_xyz_string())

                            else:
                                if opt_root == 0:
                                    current_basis = MolecularBasis.read(mol, states_basis['gs'])
                                    molecules_to_add_info.append((mol, current_basis, [opt_root]))
                                else:
                                    current_basis = MolecularBasis.read(mol, states_basis['es'])
                                    molecules_to_add_info.append((mol, current_basis, [opt_root]))
                        self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)

            elif not self.add_conformal_structures and len(files_to_add_conf) != 0:
                
                molecules_to_add_info = []
                if self.use_minimized_structures[0]:       
                    for opt_roots in self.roots_to_follow:
                        if opt_roots not in files_to_add_conf or opt_roots not in self.use_minimized_structures[2]:
                            if opt_roots == 0:
                                current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                                molecules_to_add_info.append((molecule, current_basis, [opt_roots]))
                            else:
                                current_basis = MolecularBasis.read(molecule, states_basis['es'])
                                molecules_to_add_info.append((molecule, current_basis, [opt_roots]))
                            continue
                        root_to_add = opt_roots
                        if opt_roots == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                        
                            opt_drv = OptimizationDriver(self.drivers['gs'][0])
                            current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                            _, scf_results, _ = self.compute_energy(self.drivers['gs'][0], molecule, current_basis)
                            opt_drv.ostream.mute()
                            
                            opt_constraint_list = []
                            for constraint in self.use_minimized_structures[1]:
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
                            
                            opt_results = opt_drv.compute(molecule, current_basis, scf_results)
                            energy = opt_results['opt_energies'][-1]
                            optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                            print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                            molecule = optimized_molecule

                        elif opt_roots >= 1 and isinstance(self.drivers['es'][0], LinearResponseEigenSolver) or opt_roots >= 1 and isinstance(self.drivers['es'][0], TdaEigenSolver):
                            
                            self.drivers['es'][1].state_deriv_index = [opt_roots]
                            opt_drv = OptimizationDriver(self.drivers['es'][1])
                            current_basis = MolecularBasis.read(molecule, states_basis['es'])
                            _, _, rsp_results = self.compute_energy(self.drivers['es'][0], molecule, current_basis)
                            opt_drv.ostream.mute()
                            
                            opt_constraint_list = []
                            for constraint in self.use_minimized_structures[1]:
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

                        if root_to_add == 0:
                            current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                            molecules_to_add_info.append((molecule, current_basis, [root_to_add]))
                        else:
                            current_basis = MolecularBasis.read(molecule, states_basis['es'])
                            molecules_to_add_info.append((molecule, current_basis, [root_to_add]))

                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)
                    
                else:
                    for opt_root in self.roots_to_follow:
                        if opt_root not in files_to_add_conf:
                            continue
                        if opt_root == 0:
                            current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                            molecules_to_add_info.append((molecule, current_basis, [opt_root]))
                        else:
                            current_basis = MolecularBasis.read(molecule, states_basis['es'])
                            molecules_to_add_info.append((molecule, current_basis, [opt_root]))
               
                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)
            
            else:
                for opt_root in self.roots_to_follow:
                    if opt_root == 0:
                        current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                        molecules_to_add_info.append((molecule, current_basis, [opt_root]))
                    else:
                        current_basis = MolecularBasis.read(molecule, states_basis['es'])
                        molecules_to_add_info.append((molecule, current_basis, [opt_root]))
            
            self.density_of_datapoints, self.molecules_along_rp, self.allowed_deviation = self.determine_molecules_along_dihedral_scan(molecules_to_add_info, self.roots_to_follow, specific_dihedrals=self.dihedrals_dict)    
            density_of_datapoints = self.determine_datapoint_density(self.density_of_datapoints, self.states_interpolation_settings)
            self.states_data_point_density = density_of_datapoints
            
            for counter, (state, dihedral_dict) in enumerate(self.molecules_along_rp.items()):
                for key, mol_info in dihedral_dict.items():
                    molecules, start = mol_info
                    print(f"State: {state}, Dihedral: {key}, Molecules: {molecules}, Start: {start}")
                    for i, mol in enumerate(molecules):
                        
                        if i < start:
                            continue
                        current_dihedral_angle = list(self.allowed_deviation[state][key].keys())[i]
                        forcefield_generator = MMForceFieldGenerator()
                        self.dynamics_settings['trajectory_file'] = f'trajectory_{state}_{key}_{i}.pdb'
                        forcefield_generator.partial_charges = mol.get_partial_charges(mol.get_charge())
                        forcefield_generator.create_topology(mol)
                            
                        im_database_driver = IMDatabasePointCollecter()
                        im_database_driver.distance_thrsh = self.distance_thrsh
                        im_database_driver.non_core_symmetry_groups = self.symmetry_information
                        im_database_driver.starting_state = state
                        im_database_driver.all_rot_bonds = self.all_rotatable_bonds
                    
                        im_database_driver.platform = 'CUDA'

                        im_database_driver.identfy_relevant_int_coordinates = (self.identfy_relevant_int_coordinates, self.use_minimized_structures[1])
                        im_database_driver.use_symmetry = self.use_symmetry
                        im_database_driver.ghost_atom = self.ghost_atom

                        im_database_driver.add_gpr_model = self.add_gpr_model

                        im_database_driver.use_opt_confidence_radius = self.use_opt_confidence_radius
        
                        im_database_driver.system_from_molecule(mol, self.z_matrix, forcefield_generator, solvent=self.solvent, qm_atoms='all')  

                        desired_density = False
                        current_structure_density = {}
                        
                        desiered_point_density = int(self.dynamics_settings['desired_datapoint_density'])
                        density_of_datapoints = self.determine_datapoint_density(self.density_of_datapoints, self.states_interpolation_settings)
                        self.density_of_datapoints = density_of_datapoints
                        print('density of points', self.states_data_point_density)
                        for root in density_of_datapoints.keys():

                            current_structure_density[root] = density_of_datapoints[root][key][current_dihedral_angle]
                            if density_of_datapoints[root][key][current_dihedral_angle] >= desiered_point_density:
                                print('database is already converged for state:', root)
                                break
                        print('density of points', self.states_data_point_density)

                        if desired_density is False:
                            im_database_driver.density_around_data_point = [current_structure_density, key, state, current_dihedral_angle]
                            print(im_database_driver.density_around_data_point)
                            if key is None:
                                im_database_driver.allowed_molecule_deviation = self.allowed_deviation
                            else:
                                im_database_driver.allowed_molecule_deviation = self.allowed_deviation

                            im_database_driver.update_settings(self.dynamics_settings, self.states_interpolation_settings)
                            im_database_driver.run_qmmm()
                            self.density_of_datapoints[key] = im_database_driver.density_around_data_point
                            # individual impes run objects
                            self.qm_energies.append(im_database_driver.qm_potentials)
                            self.total_energies.append(im_database_driver.total_energies)
                            self.kinetic_energies.append(im_database_driver.kinetic_energies)
                            self.state_specific_molecules = im_database_driver.state_specific_molecules
                            self.point_added_molecules.append(im_database_driver.point_adding_molecule)
                            self.unique_molecules.append(im_database_driver.allowed_molecules)
                    
                    entries = list(self.molecules_along_rp.values())
                    # self.confirm_database_quality(molecule, self.imforcefieldfiles, basis=basis, given_molecular_strucutres=self.state_specific_molecules)
                
            for root in self.roots_to_follow:
                density_of_datapoints = self.determine_datapoint_density(self.density_of_datapoints, self.states_interpolation_settings)
                self.states_data_point_density[root] = density_of_datapoints
            print('The construction of the database was sucessfull', self.states_data_point_density)
            self.im_results['n_datapoints'] = self.states_data_point_density
                    
        else:
            print('The single-state collection considers only 1 state which can get problematic if the states are not sufficiently seperated!')
            assert_msg_critical(len(self.roots_to_follow) == 1, 'ImForceFieldGenerator: The root to follow is not defined!')

            imforcefieldfile = self.imforcefieldfiles[self.roots_to_follow[0]]
            self.states_interpolation_settings[self.roots_to_follow[0]] = { 'interpolation_type':self.interpolation_type,
                                'weightfunction_type':self.weightfunction_type, 
                                'exponent_p':self.exponent_p,
                                'exponent_q':self.exponent_q, 
                                'confidence_radius':self.confidence_radius,
                                'imforcefield_file':imforcefieldfile,
                                'use_inverse_bond_length':self.use_inverse_bond_length,
                                'use_cosine_dihedral':self.use_cosine_dihedral
                            }

            self.dynamics_settings = {  'drivers':self.drivers,
                                        'basis_set_label': states_basis,
                                        'duration':self.duration, 'temperature':self.temperature, 'solvent':self.solvent,
                                        'pressure':self.force_constant, 'force_constant': self.force_constant, 'ensemble':self.ensemble,
                                        'timestep': self.timestep, 'nsteps': self.nsteps, 'friction':self.friction,
                                        'snapshots':self.snapshots, 'trajectory_file':self.trajectory_file, 'reference_struc_energy_file':self.reference_struc_energy_file,
                                        'desired_datapoint_density':self.desired_point_density, 'converged_cycle': self.converged_cycle, 
                                        'energy_threshold':self.energy_threshold, 'grad_rmsd_thrsh': self.gradient_rmsd_thrsh, 'force_orient_thrsh':self.force_orient_thrsh,
                                        'NAC':False, 'load_system': None, 'collect_qm_points_from':self.start_collect, 'roots_to_follow':self.roots_to_follow}
            
            files_to_add_conf = []
            molecules_to_add_info = []
            # if os.path.exists(self.imforcefieldfiles[self.roots_to_follow[0]]):
            #     if self.roots_to_follow[0] == 0:
            #         molecules_to_add_info.append((molecule, states_basis['gs'], self.roots_to_follow[0]))
            #     else:
            #         molecules_to_add_info.append((molecule, states_basis['gs'], self.roots_to_follow[0]))


                

            if not os.path.exists(self.imforcefieldfiles[self.roots_to_follow[0]]):
                    files_to_add_conf.append(self.roots_to_follow[0])

            if self.atom_transfer_reaction_path is not None and not os.path.exists(imforcefieldfile):

                for reaction_mol in self.atom_transfer_reaction_path:
                    
                    if self.use_minimized_structures[0]:       
                        optimized_molecule = None
                        if self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                        
                            opt_drv = OptimizationDriver(self.drivers['gs'][0])
                            current_basis = MolecularBasis.read(reaction_mol, states_basis['gs'])
                            _, scf_results, _ = self.compute_energy(self.drivers['gs'][0], reaction_mol, current_basis)
                            opt_drv.ostream.mute()
                            
                            opt_constraint_list = []
                            for constraint in self.use_minimized_structures[1]:
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

                            opt_results = opt_drv.compute(reaction_mol, current_basis, scf_results)
                            optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])

                            print(optimized_molecule.get_xyz_string())                           

                        current_basis = MolecularBasis.read(optimized_molecule, basis.get_main_basis_label())
                        self.add_point(optimized_molecule, current_basis, self.states_interpolation_settings, symmetry_information=self.symmetry_information, files_to_add=files_to_add_conf)



            if self.add_conformal_structures and not os.path.exists(imforcefieldfile):
                for counter, entry in enumerate(self.conformal_structures.items()):

                    key, molecules = entry
                    molecules_to_add_info = []
                    for i, mol_const in enumerate(molecules):
                        mol = mol_const[0]
                        constraint_conf = [f"freeze dihedral {mol_const[1][0]} {mol_const[1][1]} {mol_const[1][2]} {mol_const[1][3]}"]

                        if self.use_minimized_structures[0]:       
                            optimized_molecule = None
                            if self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                            

                                opt_drv = OptimizationDriver(self.drivers['gs'][0])
                                current_basis = MolecularBasis.read(mol, states_basis['gs'])
                                _, scf_results, _ = self.compute_energy(self.drivers['gs'][0], mol, current_basis)
                                opt_drv.ostream.mute()
                                
                                opt_constraint_list = []
                                for constraint in self.use_minimized_structures[1]:
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
                                
                                opt_results = opt_drv.compute(mol, current_basis, scf_results)
                                energy = opt_results['opt_energies'][-1]
                                optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])

                                current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                                molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow))

                                print(optimized_molecule.get_xyz_string())
                            
                            elif self.roots_to_follow[0] >= 1 and isinstance(self.drivers['es'][0], LinearResponseEigenSolver) or self.roots_to_follow[0] >= 1 and isinstance(self.drivers['es'][0], TdaEigenSolver):

                                    opt_drv = OptimizationDriver(self.drivers['es'][1])
                                    current_basis = MolecularBasis.read(mol, states_basis['es'])
                                    _, _, rsp_results  = self.compute_energy(self.drivers['es'][0], mol, current_basis)
                                    opt_drv.ostream.mute()
                                    
                                    opt_constraint_list = []
                                    for constraint in self.use_minimized_structures[1]:
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
                                    energy = opt_results['opt_energies'][-1]
                                    optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                                    print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                                    current_basis = MolecularBasis.read(optimized_molecule, states_basis['es'])
                                    molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow))
                                
                            self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)


                        else:
                            if self.roots_to_follow[0] == 0:
                                current_basis = MolecularBasis.read(mol, states_basis['gs'])
                                molecules_to_add_info.append((mol, current_basis, self.roots_to_follow))
                            else:
                                current_basis = MolecularBasis.read(mol, states_basis['es'])
                                molecules_to_add_info.append((mol, current_basis, self.roots_to_follow))
                            
                            self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)

            elif not self.add_conformal_structures and not os.path.exists(imforcefieldfile):
        
                molecules_to_add_info = []
                if self.use_minimized_structures[0]:       
                    optimized_molecule = None
                    if self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                    

                        opt_drv = OptimizationDriver(self.drivers['gs'][0])
                        current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                        _, scf_results, _ = self.compute_energy(self.drivers['gs'][0], molecule, current_basis)
                        opt_drv.ostream.mute()
                        
                        opt_constraint_list = []
                        for constraint in self.use_minimized_structures[1]:
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
                        
                        opt_results = opt_drv.compute(molecule, current_basis, scf_results)
                        energy = opt_results['opt_energies'][-1]
                        optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])

                        current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                        molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow))

                        print(optimized_molecule.get_xyz_string())
                    
                    elif self.roots_to_follow[0] >= 1 and isinstance(self.drivers['es'][0], LinearResponseEigenSolver) or self.roots_to_follow[0] >= 1 and isinstance(self.drivers['es'][0], TdaEigenSolver):

                            opt_drv = OptimizationDriver(self.drivers['es'][1])
                            current_basis = MolecularBasis.read(molecule, states_basis['es'])
                            _, _, rsp_results  = self.compute_energy(self.drivers['es'][0], molecule, current_basis)
                            opt_drv.ostream.mute()
                            
                            opt_constraint_list = []
                            for constraint in self.use_minimized_structures[1]:
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
                            energy = opt_results['opt_energies'][-1]
                            optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                            print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                            current_basis = MolecularBasis.read(optimized_molecule, states_basis['es'])
                            molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow))
                        
                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)


                else:
                    if self.roots_to_follow[0] == 0:
                        current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                        molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow))
                    else:
                        current_basis = MolecularBasis.read(molecule, states_basis['es'])
                        molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow))
                    
                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)

            else:
                
                if self.roots_to_follow[0] == 0:
                    current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                    molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow))
                else:
                    current_basis = MolecularBasis.read(molecule, states_basis['es'])
                    molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow))

            self.density_of_datapoints, self.molecules_along_rp, self.allowed_deviation = self.determine_molecules_along_dihedral_scan(molecules_to_add_info, self.roots_to_follow, specific_dihedrals=self.dihedrals_dict)
            density_of_datapoints = self.determine_datapoint_density(self.density_of_datapoints, self.states_interpolation_settings)
            self.states_data_point_density = density_of_datapoints
            
            if self.add_structures_along_rcs:
                
                for counter, (states, dihedral_dict) in enumerate(self.molecules_along_rp.items()):
                    for key, mol_info in dihedral_dict.items():
                        molecules, start = mol_info
                        molecules_to_add = []
                        print(f"State: {states}, Dihedral: {key}, Molecules: {molecules}, Start: {start}")
                        # Do your processing here
                        for i, mol in enumerate(molecules):
                            optimized_molecule = mol
                            if self.use_minimized_structures[0]:   
                                energy = None
                                self.use_minimized_structures[1].append(key)
                                if self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver):

                                    opt_drv = OptimizationDriver(self.drivers['gs'][0])
                                    current_basis = MolecularBasis.read(mol, 'def2-svp')
                                    _, scf_results, _ = self.compute_energy(self.drivers['gs'][0], mol, current_basis)
                                    opt_drv.ostream.mute()
                                    
                                    opt_constraint_list = []
                                    for constraint in self.use_minimized_structures[1]:
                                        if len(constraint) == 2:
                                            opt_constraint = f"freeze distance {constraint[0]} {constraint[1]}"
                                            opt_constraint_list.append(opt_constraint)
                                        
                                        elif len(constraint) == 3:
                                            opt_constraint = f"freeze angle {constraint[0]} {constraint[1]} {constraint[2]}"
                                            opt_constraint_list.append(opt_constraint)
                                    
                                        else:
                                            opt_constraint = f"freeze dihedral {constraint[0]} {constraint[1]} {constraint[2]} {constraint[3]}"
                                            opt_constraint_list.append(opt_constraint)
                                    opt_drv.constraints = opt_constraint_list
                                    
                                    opt_results = opt_drv.compute(mol, current_basis, scf_results)
                                    optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                                    energy = opt_results['opt_energies'][-1]
                                    print(optimized_molecule.get_xyz_string())
                            
                            if self.roots_to_follow[0] == 0 and energy is None:

                                current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                                energy, scf_results, _ = self.compute_energy(self.drivers['gs'][0], mol, current_basis)
                            elif self.roots_to_follow[0] > 0 and energy is None:
                                current_basis = MolecularBasis.read(molecule, states_basis['es'])
                                energy, scf_results, _ = self.compute_energy(self.drivers['es'][0], mol, current_basis)
                            impes_driver = InterpolationDriver(self.z_matrix)
                            impes_driver.update_settings(self.states_interpolation_settings[self.roots_to_follow[0]])
                            if self.roots_to_follow[0] == 0:
                                impes_driver.symmetry_information = self.symmetry_information['gs']
                            else:
                                impes_driver.symmetry_information = self.symmetry_information['es']
                            
                            old_label = None
                            im_labels, _ = impes_driver.read_labels()
                            impes_driver.qm_data_points = []
                            for label in im_labels:
                                if '_symmetry' not in label:
                                    qm_data_point = InterpolationDatapoint(self.z_matrix)
                                    qm_data_point.read_hdf5(self.states_interpolation_settings[self.roots_to_follow[0]]['imforcefield_file'], label)
                                    
                                    old_label = qm_data_point.point_label
                                    impes_driver.qm_symmetry_data_points[old_label] = [qm_data_point]
                                    impes_driver.qm_data_points.append(qm_data_point)

                                else:
                                    symmetry_data_point = InterpolationDatapoint(self.z_matrix)
                                    symmetry_data_point.read_hdf5(self.states_interpolation_settings[self.roots_to_follow[0]]['imforcefield_file'], label)
                                    
                                    impes_driver.qm_symmetry_data_points[old_label].append(symmetry_data_point)

                            impes_driver.compute(optimized_molecule)
                            energy_difference = (abs(energy - impes_driver.impes_coordinate.energy))

                            if energy_difference * hartree_in_kcalpermol() > self.energy_threshold:
                                print(f'Energy difference {energy_difference} is larger than the threshold {self.energy_threshold}, adding point to the database')

                                self.add_point(optimized_molecule, current_basis, self.states_interpolation_settings, symmetry_information=self.symmetry_information)
            
            for counter, (state, dihedral_dict) in enumerate(self.molecules_along_rp.items()):
                for key, mol_info in dihedral_dict.items():
                    molecules, start = mol_info
                    print(f"State: {state}, Dihedral: {key}, Molecules: {molecules}, Start: {start}")
                    # Do your processing here
                    for i, mol in enumerate(molecules):
                        
                        if i < start:
                            continue

                        current_dihedral_angle = list(self.allowed_deviation[state][key].keys())[i]

                        forcefield_generator = MMForceFieldGenerator()
                        self.dynamics_settings['trajectory_file'] = f'trajectory_{state}_{key}_{i}.pdb'
                        forcefield_generator.partial_charges = mol.get_partial_charges(mol.get_charge())
                        
                        forcefield_generator.create_topology(mol)
    
                        im_database_driver = IMDatabasePointCollecter()
                        im_database_driver.distance_thrsh = self.distance_thrsh
                        im_database_driver.non_core_symmetry_groups = self.symmetry_information
                        im_database_driver.platform = 'CUDA'
                        im_database_driver.all_rot_bonds = self.all_rotatable_bonds
                        
                        # set optimization features in the construction run

                        im_database_driver.identfy_relevant_int_coordinates = (self.identfy_relevant_int_coordinates, self.use_minimized_structures[1])
                        im_database_driver.use_symmetry = self.use_symmetry

                        im_database_driver.add_gpr_model = self.add_gpr_model
                        im_database_driver.ghost_atom = self.ghost_atom
                        im_database_driver.use_opt_confidence_radius = self.use_opt_confidence_radius
                        im_database_driver.eq_bond_force_constants = self.eq_bond_force_constants
                        
                        
                        im_database_driver.system_from_molecule(mol, self.z_matrix, forcefield_generator, solvent=self.solvent, qm_atoms='all')  

                        if self.bias_force_reaction_prop is not None:
                            im_database_driver.bias_force_reaction_idx = self.bias_force_reaction_idx
                            im_database_driver.bias_force_reaction_prop = self.bias_force_reaction_prop
                        desired_density = False
                        current_structure_density = {}
                        
                        desiered_point_density = int(self.dynamics_settings['desired_datapoint_density'])
                        density_of_datapoints = self.determine_datapoint_density(self.states_data_point_density, self.states_interpolation_settings)
                        self.density_of_datapoints = density_of_datapoints
                        print('density of points', self.states_data_point_density)
                        for root in density_of_datapoints.keys():

                            current_structure_density[root] = density_of_datapoints[root][key][current_dihedral_angle]
                            if density_of_datapoints[root][key][current_dihedral_angle] >= desiered_point_density:
                                desiered_point_density = True
                                print('database is already converged for state:', root)
                                continue
                        
        
                        if desired_density is False:
                            im_database_driver.density_around_data_point = [current_structure_density, key, state, current_dihedral_angle]
                            if key is None:
                                im_database_driver.allowed_molecule_deviation = self.allowed_deviation
                            else:
                                im_database_driver.allowed_molecule_deviation = self.allowed_deviation

                            im_database_driver.update_settings(self.dynamics_settings, self.states_interpolation_settings)
                            im_database_driver.run_qmmm()
                            self.density_of_datapoints[key] = im_database_driver.density_around_data_point
                            # individual impes run objects
                            self.qm_energies.append(im_database_driver.qm_potentials)
                            self.total_energies.append(im_database_driver.total_energies)
                            self.kinetic_energies.append(im_database_driver.kinetic_energies)
                            self.state_specific_molecules = im_database_driver.state_specific_molecules
                            self.point_added_molecules.append(im_database_driver.point_adding_molecule)
                            self.unique_molecules.append(im_database_driver.allowed_molecules)
                        
                        entries = list(self.molecules_along_rp.values())
                        # self.confirm_database_quality(molecule, self.states_interpolation_settings, basis=basis, given_molecular_strucutres=self.state_specific_molecules)

            
            desiered_point_density = int(self.dynamics_settings['desired_datapoint_density'])
            density_of_datapoints = self.determine_datapoint_density(self.states_data_point_density, self.states_interpolation_settings)
            self.states_data_point_density = density_of_datapoints
            
            print('The construction of the database was sucessfull', self.states_data_point_density)
            self.im_results['n_datapoints'] = self.states_data_point_density

        return self.im_results 

    
    def determine_atom_transfer_reaction_path(self, reactants, products, scf=True):

        reactants_partial_charges_list = []
        products_partial_charges_list = []
        
        for reactant in reactants:
            
            current_charge = reactant.get_partial_charges(reactant.get_charge())
            reactants_partial_charges_list.append(current_charge)
        
        for product in products:

            current_charge = product.get_partial_charges(product.get_charge())
            products_partial_charges_list.append(current_charge)

        ts_guesser = TransitionStateGuesser()
        TS_mol, results = ts_guesser.find_TS(reactants, products, reactant_partial_charges=reactants_partial_charges_list, product_partial_charges=products_partial_charges_list, scf=True)

        molecules_along_reaction_path = []
        for xyz_string in results['xyz_geometries']:
            molecules_along_reaction_path.append(Molecule.from_xyz_string(xyz_string))


        return molecules_along_reaction_path

    def determine_molecules_along_dihedral_scan(self, molecules_to_add_info, roots_to_follow, specific_dihedrals=None):
        """

        Sample molecular structures by rotating specific dihedrals if defined and determine the current density of 
        datapoints for given structure (or structures).

        :param molecule: The original molecule object on which rotations are performed.


        :param specific_dihedrals: A list of dihedral angle definitions (as tuples of atoms) that 
                                will be scanned. If not provided, no specific dihedrals are rotated.


        The method creates sampled molecular structures by setting each specified dihedral angle to different 
        rotation values. The rotated structures are stored in `sampled_molecules`. Additionally, it initializes
        or updates `point_densities`, which tracks how many data points exist for each dihedral configuration.
        If no specific dihedrals are provided, the method uses a default structure at a 180-degree dihedral.

        :returns:

        - sampled_molecules: A dictionary where each key is a specific dihedral (or `None` if no dihedral is given),
                            and the value is a list of sampled molecular structures.
        
        - point_densities: A dictionary where keys are tuples of (dihedral, angle) and values represent the 
                            number of existing quantum mechanical data points for that configuration.

        - allowed_deviations: The allowed angle deviation within the dynamics for the given conformer.
        """

        

        
        sampled_molecules = {root: {} for root in roots_to_follow}
        point_densities = {root: {} for root in roots_to_follow}
        allowed_deviation = {root: {} for root in roots_to_follow}
        
        molecules_info = {states[0]: molecule for molecule, _, states in molecules_to_add_info}

        if specific_dihedrals is not None:
            for entries in specific_dihedrals:
                specific_dihedral = entries[0]
                n_sampling = entries[1]
                state = entries[2]
                start = entries[3]
                molecule = molecules_info[state]
 
                rotation_values = np.linspace(0, 360, n_sampling, endpoint=False)

                sampled_molecules[state][specific_dihedral] = ([], start)
                normalized_angle = (360) / (2 * n_sampling)
                dihedral_in_deg = molecule.get_dihedral_in_degrees([specific_dihedral[0], specific_dihedral[1], specific_dihedral[2], specific_dihedral[3]])
                
                allowed_deviation[state][specific_dihedral] = {int(rotation_values[i] + dihedral_in_deg): (
                                                        ((rotation_values[i] + dihedral_in_deg) - normalized_angle)%360.0,
                                                        ((rotation_values[i] + dihedral_in_deg) + normalized_angle)%360.0
                                                        )
                                                        for i in range(len(rotation_values))}
                point_densities[state][specific_dihedral] = {int(rotation_values[i] + dihedral_in_deg): 0 for i in range(len(rotation_values))}
                for theta in rotation_values:
                    rotation_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
                    rotation_molecule.set_dihedral_in_degrees([specific_dihedral[0], specific_dihedral[1], specific_dihedral[2], specific_dihedral[3]], dihedral_in_deg + theta)
                    new_molecule = Molecule.from_xyz_string(rotation_molecule.get_xyz_string())
                    sampled_molecules[state][specific_dihedral][0].append(new_molecule)
        

        else:
            for root in roots_to_follow:
                molecule = molecules_info[root]
                sampled_molecules[root][None] = ([molecule], 0)
                point_densities[root][None] = {360: 0}
                
                allowed_deviation[root][None] = {360: (0.0, 360.0)}

        return point_densities, sampled_molecules, allowed_deviation
    
    def determine_conformal_structures(self, molecule, specific_dihedrals=None):
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

        for specific_dihedral, periodicity, n_sampling in specific_dihedrals:
            rotation_values = np.linspace(0, 360 / periodicity, n_sampling, endpoint=False)
            sampled_molecules[specific_dihedral] = []
            for theta in rotation_values:
                new_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
                new_molecule.set_dihedral_in_degrees([specific_dihedral[0] + 1, specific_dihedral[1] + 1, specific_dihedral[2] + 1, specific_dihedral[3] + 1], theta)
                sampled_molecules[specific_dihedral].append(new_molecule)

        return sampled_molecules

    def determine_datapoint_density(self, point_densities_dict, imforcefieldfile):
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
        
        reseted_point_densities_dict = {outer_key: {inner_key: {key: 0 for key in point_densities_dict[outer_key][inner_key].keys()} for inner_key in point_densities_dict[outer_key].keys()} for outer_key in point_densities_dict.keys()}

        for state in point_densities_dict.keys():
            qm_datapoints = []
            if imforcefieldfile[state]['imforcefield_file'] in os.listdir(os.getcwd()):
                impes_driver = InterpolationDriver(self.z_matrix)
                impes_driver.update_settings(imforcefieldfile[state])
                self.qmlabels, self.z_matrix = impes_driver.read_labels()
                for label in self.qmlabels:
                    if '_symmetry' not in label:
                        qm_data_point = InterpolationDatapoint(self.z_matrix)
                        qm_data_point.update_settings(imforcefieldfile[state])
                        qm_data_point.read_hdf5(imforcefieldfile[state]['imforcefield_file'], label)
                        qm_datapoints.append(qm_data_point)
            for specific_dihedral in point_densities_dict[state].keys():
                for point in qm_datapoints:
                    if specific_dihedral is None:
                        reseted_point_densities_dict[state][specific_dihedral][360] += 1
                    else:
                        min_distance = np.inf
                        key = None
                        for dihedral in point_densities_dict[state][specific_dihedral].keys():
                            datapoint_molecule = Molecule(self.molecule.get_labels(), point.cartesian_coordinates, 'bohr')
                            dihedrals_of_dp = [datapoint_molecule.get_dihedral_in_degrees([specific_dihedral[0], specific_dihedral[1], specific_dihedral[2], specific_dihedral[3]])]
                            distance_vectorized = dihedral_distance_vectorized([dihedral], dihedrals_of_dp)
                            if abs(distance_vectorized) < min_distance:
                                min_distance = abs(distance_vectorized)
                                key = dihedral
                    
                        reseted_point_densities_dict[state][specific_dihedral][key] += 1

        return reseted_point_densities_dict


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
    
    def database_extracter(self, datafile, mol_labels, im_settings):
        """Extracts molecular structures from a given database file.

        :param datafile:
            Database file containing interpolation data.
        
        :param mol_labels:
            List of molecular labels.

        :returns:
            A list of VeloxChem Molecule objects extracted from the database.
        
        """

        
        im_driver = InterpolationDriver() # -> implemented Class in VeloxChem that is capable to perform interpolation calculations for a given molecule and provided z_matrix and database
        im_driver.update_settings(im_settings)
        # im_driver.imforcefield_file = datafile
        labels, z_matrix = im_driver.read_labels()
    
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        # impes_coordinate = InterpolationDatapoint(z_matrix) # -> implemented Class in VeloxChem that handles all transformations and database changes concerning the interpolation
        data_point_molecules = []
        datapoints = []

        for label in sorted_labels:
            impes_coordinate = InterpolationDatapoint(z_matrix)
            impes_coordinate.update_settings(im_settings)
            impes_coordinate.read_hdf5(datafile, label) # -> read in function from the ImpesDriver object

            coordinates_in_angstrom = impes_coordinate.cartesian_coordinates * bohr_in_angstrom()

           
            current_molecule = Molecule(mol_labels, coordinates_in_angstrom, 'angstrom') # -> creates a VeloxChem Molecule object

            datapoints.append(impes_coordinate)
            data_point_molecules.append(current_molecule)

        return data_point_molecules, datapoints, z_matrix
    
    def confirm_database_quality(self, molecule, states_interpolation_settings, basis, given_molecular_strucutres=None, improve=True):
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
           List of QM-energies, IM-energies.
        """

        # For all Methods a ForceField of the molecule is requiered
        forcefield_generator = MMForceFieldGenerator()
        forcefield_generator.create_topology(molecule)
        
        # self.imforcefieldfile = im_database_files
        # if self.interpolation_settings is None:
        #     self.interpolation_settings = { 'interpolation_type':self.interpolation_type, 
        #                     'exponent_p':self.exponent_p,
        #                     'exponent_q':self.exponent_q, 
        #                     'confidence_radius':self.confidence_radius,
        #                     'imforcefield_file':self.imforcefieldfile,
        #                     'use_inverse_bond_length':True
        #                   }

        if self.dynamics_method == 'MM':
            rot_bonds = forcefield_generator.rotatable_bonds
            forcefield_generator.reparameterize_dihedrals(rot_bonds[0], scan_range=[180, 360], n_points=7, visualize=True)

        for root in self.roots_to_follow:
            database_quality = False
            drivers = None
            if root == 0:
                drivers = self.drivers['gs']
            else:
                drivers = self.drivers['es']
            all_structures = given_molecular_strucutres[root]
            current_datafile = states_interpolation_settings[root]['imforcefield_file']
            datapoint_molecules, _ = self.database_extracter(current_datafile, molecule.get_labels())
            
            if len(all_structures) == 0:
                continue

            while database_quality is False:
                
                database_expanded = False
                if given_molecular_strucutres is None:
                    print('provide reference structures!')
                    exit()
                    
                rmsd = -np.inf
                random_structure_choices = None
                counter = 0
                # if given_molecular_strucutres is not None:
                #     random_structure_choices = given_molecular_strucutres
                dist_ok = False
                while dist_ok == False and counter <= 20:
                    print('Dihedrals dict is given here', self.dihedrals_dict)
                    if self.dihedrals_dict is not None:
                        selected_molecules = []
                        for entries in self.dihedrals_dict:
                            specific_dihedral = entries[0]
                            n_sampling = entries[1]
                            state = entries[2]
                            if state != root:
                                continue
                            start = entries[3]
                            desired_angles = np.linspace(0, 360, 36)
                            angles_mols = {int(angle):[] for angle in desired_angles}

                            keys = list(angles_mols.keys())

                            for mol in all_structures:  
                                mol_angle = (mol.get_dihedral_in_degrees(specific_dihedral) + 360) % 360
                                
                                
                                for i in range(len(desired_angles) - 1):
                                    
                                    if keys[i] <= mol_angle < keys[i + 1]:
                                        angles_mols[keys[i]].append(mol)
                                        break
                        
                            # List to hold selected molecules
                            
                            total_molecules = sum(len(mols) for mols in angles_mols.values())

                            # Adaptive selection based on bin sizes
                
                            for angle_bin, molecules_in_bin in angles_mols.items():
                                num_mols_in_bin = len(molecules_in_bin)

                                if num_mols_in_bin == 0:
                                    continue  # Skip empty bins

                                # If 2 or fewer molecules, take all
                                elif num_mols_in_bin <= 2:
                                    selected_molecules.extend(molecules_in_bin)

                                else:
                                    # Calculate proportional number of molecules to select
                                    proportion = num_mols_in_bin / total_molecules
                                    num_to_select = max(1, math.ceil(proportion * self.nstruc_to_confirm_database_quality))

                                    # Randomly select the proportional number
                                    selected_mols = random.sample(molecules_in_bin, min(num_to_select, num_mols_in_bin))
                                    selected_molecules.extend(selected_mols)
                    else:
                        selected_molecules = random.sample(all_structures, min(self.nstruc_to_confirm_database_quality, len(all_structures)))
                          
                    individual_distances = []
                    random_structure_choices = selected_molecules
                    for datapoint_molecule in datapoint_molecules:
                        for random_struc in random_structure_choices:
                            
                            distance_norm = self.calculate_distance_to_ref(random_struc.get_coordinates_in_bohr(), datapoint_molecule.get_coordinates_in_bohr())
                            individual_distances.append(distance_norm / np.sqrt(len(molecule.get_labels())) * bohr_in_angstrom())
                    
                    rmsd = min(individual_distances)
                    counter += 1
                    if rmsd >= 0.1:
                        print(f'The overall RMSD is {rmsd} -> The current structures are well seperated from the database conformations! loop is discontinued')
                        dist_ok = True
                    else:
                        print(f'The overall RMSD is {rmsd} -> The current structures are not all well seperated from the database conformations! loop is continued')        
                
                # if self.roots_to_follow[0] == 0 and energy is None:
                #     energy, scf_results = self.compute_energy(self.drivers['gs'][0], optimized_molecule, current_basis)
                # elif self.roots_to_follow[0] > 0 and energy is None:
                #     energy, scf_results = self.compute_energy(self.drivers['es'][0], optimized_molecule, current_basis)
                impes_driver = InterpolationDriver(self.z_matrix)
                impes_driver.update_settings(states_interpolation_settings[root])
                if root == 0:
                    impes_driver.symmetry_information = self.symmetry_information['gs']
                else:
                    impes_driver.symmetry_information = self.symmetry_information['es']
                
                old_label = None
                im_labels, _ = impes_driver.read_labels()
                impes_driver.qm_data_points = []
                for label in im_labels:
                    if '_symmetry' not in label:
                        qm_data_point = InterpolationDatapoint(self.z_matrix)
                        qm_data_point.read_hdf5(current_datafile, label)
                        
                        old_label = qm_data_point.point_label
                        impes_driver.qm_symmetry_data_points[old_label] = [qm_data_point]
                        impes_driver.qm_data_points.append(qm_data_point)

                    else:
                        symmetry_data_point = InterpolationDatapoint(self.z_matrix)
                        symmetry_data_point.read_hdf5(current_datafile, label)
                        
                        impes_driver.qm_symmetry_data_points[old_label].append(symmetry_data_point)


                qm_energies = []
                im_energies = []

                for i, mol in enumerate(random_structure_choices):

                    current_basis = MolecularBasis.read(mol, basis.get_main_basis_label())
                    impes_driver.compute(mol)
                    reference_energy, scf_results, _ = self.compute_energy(drivers[0], mol, current_basis)
                    
                    current_element = 0
                    if root >= 1:
                        current_element = root - 1
                    
                    qm_energies.append(reference_energy[current_element])
                    im_energies.append(impes_driver.impes_coordinate.energy)
                    
                    print('Energies', qm_energies[-1], im_energies[-1])
                    
                    print(f'\n\n ########## Step {i} ######### \n')
                    print(f'delta_E:   {abs(qm_energies[-1] - im_energies[-1]) * hartree_in_kcalpermol()} kcal/mol \n')
                    if abs(qm_energies[-1] - im_energies[-1]) * hartree_in_kcalpermol() > self.energy_threshold and improve == True:

                        print(mol.get_xyz_string())

                        if self.use_minimized_structures[0]:
                            current_weights = impes_driver.weights

                            weights = [value for _, value in current_weights.items()]
                            used_labels = [label_idx for label_idx, _ in current_weights.items()]

                            # Sort labels and weights by descending weight
                            sorted_items = sorted(zip(used_labels, weights), key=lambda x: x[1], reverse=True)

                            total_weight = sum(weights)
                            cumulative_weight = 0.0
                            internal_coordinate_datapoints = []

                            for label, weight in sorted_items:
                                cumulative_weight += weight
                                internal_coordinate_datapoints.append(impes_driver.qm_data_points[label])
                                if cumulative_weight >= 0.8 * total_weight:
                                    break

                            # qm_datapoints_weighted = [qm_datapoint for qm_datapoint in enumerate if ]
                            constraints = impes_driver.determine_important_internal_coordinates(qm_energies[-1], mol, self.z_matrix, internal_coordinate_datapoints)

                            print('CONSTRAINTS', constraints)

                            
                            if isinstance(drivers[0], ScfRestrictedDriver):

                                _, scf_tensors = self.compute_energy(drivers[0], molecule, current_basis)
                                opt_drv = OptimizationDriver(drivers[0])
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
                                
                                for constraint in self.use_minimized_structures[1]:
                                    if constraint in opt_constraint_list:
                                        continue
                                    if len(constraint) == 2:
                                        opt_constraint = f"freeze distance {constraint[0]} {constraint[1]}"
                                        opt_constraint_list.append(opt_constraint)
                                    
                                    elif len(constraint) == 3:
                                        opt_constraint = f"freeze angle {constraint[0]} {constraint[1]} {constraint[2]}"
                                        opt_constraint_list.append(opt_constraint)
                                
                                    else:
                                        opt_constraint = f"freeze dihedral {constraint[0]} {constraint[1]} {constraint[2]} {constraint[3]}"
                                        opt_constraint_list.append(opt_constraint)
                                opt_drv.constraints = opt_constraint_list
                                opt_results = opt_drv.compute(mol, current_basis, scf_tensors)
                                optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                                optimized_molecule.set_charge(mol.get_charge())
                                optimized_molecule.set_multiplicity(mol.get_multiplicity())
                                
                                # qm_energy, scf_tensors = self.compute_energy(drivers[0], optimized_molecule, opt_current_basis)
                                print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', mol.get_xyz_string())
                            
                            elif isinstance(drivers[0], ExternalScfDriver):
                                
                                for constraint in self.use_minimized_structures[1]:
                                    if constraint in constraints:
                                        continue
                                    constraints.append(constraint)
                                optim_driver = ExternalOptimDriver(drivers[0])
                                optim_driver.constraints = constraints
                                opt_mol_string, energy = optim_driver.optimize(molecule)
                                optimized_molecule = Molecule.from_xyz_string(opt_mol_string)

                                print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                            current_basis = MolecularBasis.read(optimized_molecule, current_basis.get_main_basis_label())      
                            labels = []
                            labels.append("point_{0}".format((len(im_labels) + 1)))
                            self.add_point(optimized_molecule, current_basis, states_interpolation_settings, symmetry_information=self.symmetry_information)
                            database_expanded = True
                            print('The interpolation quality was to low! Structre as been added to the database')
                        
                        else:
                            labels = []
                            labels.append("point_{0}".format((len(im_labels) + 1)))
                            current_basis = MolecularBasis.read(mol, current_basis.get_main_basis_label())
                            self.add_point(mol, current_basis, states_interpolation_settings, symmetry_information=self.symmetry_information)
                            database_expanded = True

                if not database_expanded:
                    database_quality = True

            
            # self.plot_final_energies(qm_energies, im_energies)
            self.structures_to_xyz_file(all_structures, 'full_xyz_traj.xyz')
            self.structures_to_xyz_file(random_structure_choices, 'random_xyz_structures.xyz', im_energies, qm_energies)


    
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

            if len(xyz_lines) >= 2 and im_energies is not None or len(xyz_lines) >= 2 and qm_energies is not None:
                
                if im_energies is not None and qm_energies is None:
                    xyz_lines[1] += f'Energies  IM: {im_energies[i]}'
                elif im_energies is None and qm_energies is not None:
                    xyz_lines[1] += f'Energies  IM: {qm_energies[i]}'
                else:
                    xyz_lines[1] += f'Energies  QM: {qm_energies[i]}  IM: {im_energies[i]}  delta_E: {abs(qm_energies[i] - im_energies[i])}'


            updated_xyz_string = "\n".join(xyz_lines)

            with open(structure_filename, 'a') as file:
                file.write(f"{updated_xyz_string}\n\n")
        
    
    def define_z_matrix(self, molecule, add_coordinates=None):
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
        
        # primitives = geometric.internal.PrimitiveInternalCoordinates(g_molecule)
        
        # for Internal in primitives:
        #     print(Internal.value())
        # exit()


        z_matrix = []
        for bond in bonds:
            z_matrix.append(bond)        
        if add_coordinates is not None:
            for coord in add_coordinates:
                if len(coord) == 2:
                    z_matrix.append(coord)
        for angle in angles:
            z_matrix.append(angle)
        if add_coordinates is not None:
            for coord in add_coordinates:
                if len(coord) == 3:
                    z_matrix.append(coord)
        for dihedral in dihedrals:
            z_matrix.append(dihedral)
            if self.use_cosine_dihedral:
                z_matrix.append(dihedral)
        if add_coordinates is not None:
            for coord in add_coordinates:
                if len(coord) == 4:
                    z_matrix.append(coord)
                    if self.use_cosine_dihedral:
                        z_matrix.append(coord)

        
        return z_matrix
    

    def add_point(self, molecule_specific_information, interpolation_settings, symmetry_information={}):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param imforcefielddatafile:
                Datafile containing the information of the IM forcefield.

        """
        if len(self.drivers) == 0:
            raise ValueError("No energy driver defined.")
           
        # define impesdriver to determine if stucture should be added:
        
        ## For symmetry groups of periodicty of 3 it is crucial for the interpolation to set the dihedral to the position between 2 extreme points in order to account
        ## for the symmetry correclty using only one reference point!
        
        # create all molecule combinations

        adjusted_molecule = {'gs': [], 'es': []}
        symmetry_mapping_groups = []
        symmetry_exclusion_groups = []

        for entries in molecule_specific_information:
            molecule = entries[0]
            basis = entries[1]
            states = entries[2]
            symmetry_point = False
            if 0 in states and len(symmetry_information['gs']) != 0 and len(symmetry_information['gs'][2]) != 0 or 1 in states and len(symmetry_information['gs']) != 0 and self.drivers['es'][0].spin_flip and len(symmetry_information['gs'][2]) != 0:
                symmetry_mapping_groups = [item for item in range(len(molecule.get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['gs'][1] for item in element]
                sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(molecule, symmetry_information['gs'][1], symmetry_information['gs'][5])
                
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
                    constraints = []
                    for dihedral, angle in molecule_config.items():
                        opt_dihedral_angle = molecule.get_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], 'radian')
                        print(molecule.get_xyz_string())
                        print('Dihedral creation', dihedral, opt_dihedral_angle, angle)
                        
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], opt_dihedral_angle + angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                        constraint = f"freeze dihedral {dihedral[0] + 1} {dihedral[1] + 1} {dihedral[2] + 1} {dihedral[3] + 1}"
                        constraints.append(constraint)
                        if self.symmetry_information is not None:
                            print('opt_dihedral_angle', opt_dihedral_angle)
                            self.symmetry_information['gs'][8][dihedral].append(opt_dihedral_angle)
                    if i > 0:
                        symmetry_point = True
                    current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())
                    adjusted_molecule['gs'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, entries[2], symmetry_point))
            if 1 in states and len(symmetry_information['es']) != 0 and len(symmetry_information['es'][2]) != 0 and self.drivers['es'][0].spin_flip is False:
                symmetry_mapping_groups = [item for item in range(len(molecule.get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['es'][1] for item in element]
                sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(molecule, symmetry_information['es'][1], symmetry_information['es'][5])
                
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
                        print('Dihedral', dihedral)
                        opt_dihedral_angle = cur_molecule.get_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], 'radian')
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], opt_dihedral_angle + angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                        print('påt_dihedral', opt_dihedral_angle + angle)
                        if self.symmetry_information is not None:
                            self.symmetry_information['es'][8][dihedral].append(opt_dihedral_angle + angle)
                    if i > 0:
                        symmetry_point = True
                    current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())
                    adjusted_molecule['es'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, entries[2], symmetry_point))
            elif any(x > 1 for x in states) and len(symmetry_information['es']) != 0 and len(symmetry_information['es'][2]) != 0 and 1 in states and self.drivers['es'][0].spin_flip is True:
                symmetry_mapping_groups = [item for item in range(len(molecule.get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['es'][1] for item in element]
                sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(molecule, symmetry_information['es'][1], symmetry_information['es'][5])
                
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
                        print('Dihedral', dihedral)
                        opt_dihedral_angle = cur_molecule.get_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], 'radian')
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], opt_dihedral_angle + angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                        print('påt_dihedral', opt_dihedral_angle + angle)
                        if self.symmetry_information is not None:
                            self.symmetry_information['es'][8][dihedral].append(opt_dihedral_angle + angle)
                    if i > 0:
                        symmetry_point = True
                    current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())
                    adjusted_molecule['es'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, entries[2], symmetry_point))
            
            else:
                if 0 in entries[2]:
                    adjusted_molecule['gs'].append((entries[0], entries[1], 1, None, [0], symmetry_point)) 
                if any(x > 0 for x in entries[2]): 
                    states = [state for state in entries[2] if state > 0]
                    adjusted_molecule['es'].append((entries[0], entries[1], 1, None, states, symmetry_point)) 
        # else:
        #     adjusted_molecule['gs'].append((molecule, basis, 1, None))  
        #     adjusted_molecule['es'].append((molecule, basis, 1, None))  


        for key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue
            print(key, entries)
            drivers = self.drivers[key]


            org_roots = None

            if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                    org_roots = drivers[1].state_deriv_index  

            label_counter = 0
            for mol_basis in entries:
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    root_to_follow_calc = mol_basis[4]           
                    drivers[1].state_deriv_index = root_to_follow_calc 
                    # drivers[2].roots_to_follow = root_to_follow_calc
                
                
                energies, scf_results, rsp_results = self.compute_energy(drivers[0], mol_basis[0], mol_basis[1])
                print(energies)
                
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    energies = energies[mol_basis[4]]

                gradients = self.compute_gradient(drivers[1], mol_basis[0], mol_basis[1], scf_results, rsp_results)
                print(gradients)
                hessians = self.compute_hessian(drivers[2], mol_basis[0], mol_basis[1])
                print(hessians)
                masses = mol_basis[0].get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
                for number in range(len(energies)):
    
                    interpolation_driver = InterpolationDriver()
                    interpolation_driver.update_settings(interpolation_settings[mol_basis[4][number]])
                    interpolation_driver.imforcefield_file = interpolation_settings[mol_basis[4][number]]['imforcefield_file']
                    if self.z_matrix is None:
                        self.z_matrix = self.define_z_matrix(molecule, self.reaction_coordinates)
                    sorted_labels = []
                    if interpolation_settings[mol_basis[4][number]]['imforcefield_file'] in os.listdir(os.getcwd()):
                        org_labels, self.z_matrix = interpolation_driver.read_labels()
                        labels = [label for label in org_labels if '_symmetry' not in label]
                        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))
            
                    label = None
                    grad = gradients[number].copy()
                    hess = hessians[number].copy()

                    grad_vec = grad.reshape(-1)         # (3N,)
                    hess_mat = hess.reshape(grad_vec.size, grad_vec.size)
                    mw_grad_vec = inv_sqrt_masses * grad_vec
                    mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]
                    
    
                    if mol_basis[5] == False:
                        label = f'point_{len(sorted_labels) + 1}'
                        old_label = f'point_{len(sorted_labels) + 1}'
                        current_label = f'point_{len(sorted_labels) + 1}'
                    else:
                        label = f'{old_label}_symmetry_{label_counter}'


                    print('force dict', self.eq_bond_force_constants)

                    impes_coordinate = InterpolationDatapoint(self.z_matrix)
                    impes_coordinate.eq_bond_force_constants = self.eq_bond_force_constants
                    impes_coordinate.update_settings(interpolation_settings[mol_basis[4][number]])
                    impes_coordinate.cartesian_coordinates = mol_basis[0].get_coordinates_in_bohr()
                    impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                    impes_coordinate.energy = energies[number]
                    impes_coordinate.gradient =  mw_grad_vec.reshape(grad.shape)
                    impes_coordinate.hessian = mw_hess_mat.reshape(hess.shape)
                    impes_coordinate.transform_gradient_and_hessian()
                    
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
                                print('dihedral', angle, dihedral, (mol_basis[0].get_dihedral(dihedral, 'radian')))
                                rot_mol.set_dihedral(dihedral, mol_basis[0].get_dihedral(dihedral, 'radian') - angle, 'radian')
                                print('dihedral', np.cos(3.0 * mol_basis[0].get_dihedral(dihedral, 'radian')), np.cos(3.0 * (mol_basis[0].get_dihedral(dihedral, 'radian') - angle)))
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
                                inverse_mapping = {v: k for k, v in mapping_dict[ord].items()}
                                reorded_int_coord = np.zeros_like(impes_coordinate.internal_coordinates_values)
                                for i, element in enumerate(impes_coordinate.z_matrix):
                                    # Otherwise, reorder the element
                                    reordered_element = [inverse_mapping.get(x, x) for x in element]
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
                            
                    trust_radius = self.use_opt_confidence_radius[2]
                    impes_coordinate.confidence_radius = trust_radius
                    
                    impes_coordinate.write_hdf5(interpolation_settings[mol_basis[4][number]]['imforcefield_file'], label)
                    interpolation_driver.imforcefield_file = interpolation_settings[mol_basis[4][number]]['imforcefield_file']
                    
                    labels, z_matrix = interpolation_driver.read_labels()
                    
                    print(f"Database expansion with {', '.join(labels)}")
                
                label_counter += 1
                if any(root > 0 for root in self.roots_to_follow): 
                    if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                        drivers[1].state_deriv_index = org_roots
                    

    def calculate_translation_coordinates(self, cart_coord):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(cart_coord, axis=0)
        translated_coordinates = cart_coord - center

        return translated_coordinates

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

    def adjust_symmetry_dihedrals(self, molecule, symmetry_groups, rot_bonds):
        
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
        
        if self.z_matrix is None:
            self.z_matrix = self.define_z_matrix(molecule, self.reaction_coordinates)
        all_dihedrals = [element for element in self.z_matrix if len(element) == 4]

        symmetry_group_dihedral_dict = {} 
        angles_to_set = {}
        periodicities = {}
        dihedral_groups = {2: [], 3: []}

        for symmetry_group in symmetry_groups:
             
            symmetry_group_dihedral_list = symmetry_group_dihedral(symmetry_group, all_dihedrals, rot_bonds)
            symmetry_group_dihedral_dict[tuple(symmetry_group)] = symmetry_group_dihedral_list
            print('dihedral list', symmetry_group_dihedral_list)
            
            if len(symmetry_group) == 3:

                # angles_to_set[symmetry_group_dihedral_list[0]] = ([0.0, np.pi/3.0])
                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = ([0.0, np.pi/3.0])

                periodicities[tuple(symmetry_group_dihedral_list[0])] = 3
                dihedral_groups[3].extend([tuple(sorted(element, reverse=False)) for element in symmetry_group_dihedral_list])

            elif len(symmetry_group) == 2:

                # angles_to_set[symmetry_group_dihedral_list[0]] = ([0.0, np.pi/3.0])
                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = ([0.0, np.pi/2.0])

                periodicities[tuple(symmetry_group_dihedral_list[0])] = 2
                dihedral_groups[2].extend([tuple(sorted(element, reverse=False)) for element in symmetry_group_dihedral_list])

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
            qm_energy = np.array([qm_driver.scf_energy])
            qm_driver.ostream.unmute()
            qm_driver.filename = None
            qm_driver.checkpoint_file = None

            print('qm_energy in SCF driver', qm_energy)

        elif isinstance(qm_driver, LinearResponseEigenSolver) or isinstance(qm_driver, TdaEigenSolver):
            self.drivers['es'][0].ostream.mute()
            scf_results = self.drivers['es'][3].compute(molecule, basis)
            self.drivers['es'][3].filename = None
            self.drivers['es'][3].checkpoint_file = None
            scf_energy = self.drivers['es'][3].scf_energy
            qm_driver.ostream.unmute()
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
                self.drivers['es'][1].state_deriv_index = root
                hess_driver.compute(molecule, basis)
                qm_hessian = hess_driver.hessian
                qm_hessians.append(qm_hessian)

        if qm_hessians is None:
            error_txt = "Could not compute the QM Hessian. "
            error_txt += "Please define a QM Hessian driver."
            raise ValueError(error_txt)


        return qm_hessians