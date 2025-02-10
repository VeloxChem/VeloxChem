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

from mpi4py import MPI
import numpy as np
from scipy.optimize import linear_sum_assignment
from pathlib import Path
from sys import stdout
import sys
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom
import random
from scipy.spatial import Delaunay
from scipy.spatial import distance
import alphashape
import re

from contextlib import redirect_stderr
from io import StringIO
with redirect_stderr(StringIO()) as fg_err:
    import geometric

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from .molecule import Molecule
from .veloxchemlib import mpi_master
from. veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .forcefieldgenerator import ForceFieldGenerator
# from .systembuilder import SystemBuilder
from .profiler import Profiler
# from .interpolationmapping import InterpolationMapping
# from .findbestcombination import FindBestCombination


# Drivers
from .scfrestdriver import ScfRestrictedDriver
from .molecularbasis import MolecularBasis
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .vibrationalanalysis import VibrationalAnalysis
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
# from .externalqmdriver import ExternalQMDriver
# from .externalqmgradientdriver import ExternalQMGradientDriver
# from .externalqmhessiandriver import ExternalQMHessianDriver
from .impesdriver import ImpesDriver
from .impescoordinates import ImpesCoordinates



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
        - qm_driver: The VeloxChem driver object. Options are XtbDriver and ImpesDriver.
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
        self.platform = 'reference'
        self.ensemble = 'NVE'
        self.temperature = 298.15 
        self.friction = 1.0 
        self.timestep = 2.0
        self.nsteps = 1000
        self.parent_ff = 'amber03.xml'
        self.water_ff = 'spce.xml'
        self.box_size = 2.0 
        self.padding = 1.0
        self.cutoff = 1.0
        self.integrator = None
        self.scaling_time = 5

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
        self.energy_gabs = {}
        self.state_energies = {}
        self.unique_residues = []
        self.unique_molecules = []

        self.starting_state = 0
        self.current_im_choice = None
        self.excitation_step = 10000
        
        # QM Region parameters
        self.qm_driver = None
        self.grad_driver = None
        self.qm_atoms = None
        self.mm_subregion = None
        self.linking_atoms = None
        self.qm_force_index = None
        self.driver_flag = None
        self.swapped = False
        self.state_swtiched = False
        self.velo_switch = False

        self.snapshots = None
        self.calc_NAC = False
        self.roots = 1
        self.load_system = None
        self.pressure = None
        self.output_file = None
        self.adiabatic_basis = False
        self.density_around_data_point = None
        self.impes_drivers = None
        self.im_labels = None
        self.qm_energies = None
        self.qm_data_points = None
        self.point_adding_molecule = {}
        self.energy_threshold = 1.0
        self.collect_qm_points = 0
        self.core_structure = None
        self.non_core_structure = None
        self.non_core_symmetry_groups = None

        self.ff_datafile = None
        self.starting_temperature = None

        self.reference_calc_dict = None
        self.current_state = None
        self.current_im_choice = None
        self.current_gradient = 0
        self.current_energy = 0
        self.point_checker = 1
        self.allowed_molecule_deviation = None
        self.last_point_added = None
        self.im_labels = []
        self.add_a_point = False
        self.check_a_point = False
        self.cluster_run = None
        
        self.normal_modes_vec = {}
        self.distplacement_threshold = {}
        self.distance_threshold_dict = {}
        self.skipping_value = 0



        self.basis = None

        # output_file variables that will be written into self.general_variable_output
        self.gradients = None
        self.velocities = None
        self.coordinates = None
        self.coordinates_xyz = None
        self.molecules = []
        self.all_gradients = []
        
        
        self.summary_output = 'summary_output.txt'
        with open(self.summary_output, 'w') as file:
            file.write("########## Summaray Ouput of Structures and Energies ##########\n\n")
        self.coordinates_xyz = None
        
        self.compare_current_points = []
        self.compare_alligned_points = []
        self.compare_nearest_rot_points = []
        self.compare_nearest_points = []
        self.compare_nearest_2_points = []
        self.restoring_geometry = None
        self.velocities = []

        self.determine_convex_hull = False
        self.scaling_factor = 0.0
        # Default value for the C-H linker distance
        self.linking_atom_distance = 1.0705 
        
    # Loading methods
    # TODO: Integrate the guess with the read_PDB_file in Molecule.
    def load_system_PDB(self, filename):
        """
        Loads a system from a PDB file, extracting unique residues and their details.

        :param filename: 
            Path to the PDB file. Example: 'system.pdb'.
        """
        
        pdb_path = Path(filename)
        if not pdb_path.is_file():
            raise FileNotFoundError(f"{filename} does not exist.")
        
        pdbstr = pdb_path.read_text()
        residues = {}

        # Formatting issues flags
        label_guess_warning = False
        conect_warning = True

        self.unique_residues = []  # Initialize unique_residues earlier

        for line in pdbstr.strip().splitlines():
            # Skip the header lines
            if line.startswith(('REMARK', 'HEADER', 'TITLE', 'COMPND', 'SOURCE')):
                continue
            # Skip the CRYST1 record it will be extracted from the System object
            if line.startswith('CRYST1'):
                continue
            # Add the ATOM and HETATM records to the residues dictionary
            if line.startswith(('ATOM', 'HETATM')):
                atom_label = line[76:78].strip()
                if not atom_label:
                    label_guess_warning = True
                    atom_name = line[12:16].strip()
                    atom_label = atom_name[0]

                residue = line[17:20].strip()
                residue_number = int(line[22:26])
                coordinates = [float(line[i:i+8]) for i in range(30, 54, 8)]

                residue_identifier = (residue, residue_number)

                if residue_identifier not in residues:
                    residues[residue_identifier] = {
                        'labels': [],
                        'coordinates': []
                    }
                
                residues[residue_identifier]['labels'].append(atom_label)
                residues[residue_identifier]['coordinates'].append(coordinates)
            
            # If a line starts with TER, skip it
            if line.startswith('TER'):
                continue
            # If any line starts with CONECT, set the conect_warning flag to False
            if line.startswith('CONECT'):
                conect_warning = False

            # If a line starts with END, break the loop
            if line.startswith('END'):
                break
            

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
            msg = 'Atom labels were guessed based on atom names (first character).'
            msg += f'Please verify the atom labels in the {filename} PDB file.'
            self.ostream.print_warning(msg)
            self.ostream.flush()

        if conect_warning:
            msg = 'CONECT records not found in the PDB file.'
            msg += 'The connectivity matrix will be used to determine the bonds.'
            self.ostream.print_warning(msg)
            self.ostream.flush()
            # Create a molecule from the PDB file
            molecule = Molecule.read_pdb_file(filename)
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

        # Create VeloxChem Molecule objects for each unique residue

        molecules = []
        unq_residues = []

        for (residue, number), data in residues.items():
            coordinates_array = np.array(data['coordinates'])
            mol = Molecule(data['labels'], coordinates_array, "angstrom")
            molecules.append(mol)
            unq_residues.append((residue, number))


        # Initialize a set to track the first occurrence of each residue name
        seen_residue_names = set()

        # Lists to store the filtered residues and molecules
        self.unique_residues = []
        self.unique_molecules = []

        for index, (residue_name, number) in enumerate(unq_residues):
            if residue_name not in seen_residue_names:
                seen_residue_names.add(residue_name)
                self.unique_residues.append((residue_name, number))
                self.unique_molecules.append(molecules[index])

        # Print results
        info_msg = f"Unique Residues: {self.unique_residues}, saved as molecules."
        self.ostream.print_info(info_msg)
        self.ostream.flush()
    
    def load_system_from_files(self, system_xml, system_pdb):
        """
        Loads a preexisting system from an XML file and a PDB file.

        :param system_xml: 
            XML file containing the system parameters.
        :param system_pdb: 
            PDB file containing the system coordinates.
        """
        self.pdb = app.PDBFile(system_pdb)

        with open(system_xml, 'r') as file:
            self.system = mm.XmlSerializer.deserialize(file.read())

        self.phase = 'gas'
        self.qm_atoms = []
        self.qm_force_index = -1  # Initialize with an invalid index

        # Iterate over all forces in the system to find the CustomExternalForce and detect QM atoms
        for i, force in enumerate(self.system.getForces()):
            if isinstance(force, mm.CustomExternalForce):
                self.qm_force_index = i  # Set the index of the QM force
                # Assuming that the particle indices are directly stored in the force
                for j in range(force.getNumParticles()):
                    # Get the actual particle index (usually force index and particle index are the same)
                    index, parameters = force.getParticleParameters(j)
                    self.qm_atoms.append(index)
                msg = f"QM Atoms detected: {self.qm_atoms}"
                self.ostream.print_info(msg)
                self.ostream.flush()
                break

        if self.qm_force_index == -1:
            msg = "No CustomExternalForce found, no QM atoms detected."
            self.ostream.print_info(msg)
            self.ostream.flush()
        else:
            msg = f"CustomExternalForce found at index {self.qm_force_index}."
            self.ostream.print_info(msg)
            self.ostream.flush()
            
    # Method to generate OpenMM system from VeloxChem objects
    def system_from_molecule(self,
                             molecule,
                             z_matrix, 
                             ff_gen, 
                             solvent='gas', 
                             qm_atoms=None, 
                             filename='residue', 
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

        # Store the molecule object and generate OpenMM compatible files
        self.molecule = molecule
        self.positions = molecule.get_coordinates_in_angstrom()
        self.labels = molecule.get_labels()
        self.z_matrix = z_matrix

        # Options for the QM region if it's required.
        # TODO: Take this if else tree to a separate method.
        if qm_atoms:
            filename = 'qm_region'
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

            ff_gen.write_pdb('qm_region.pdb', 'QMR')

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
            sys_builder = SystemBuilder()
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

    # Methods to build a custom system from a PDB file and custom XML files
    def md_system_from_files(self, 
                            pdb_file, 
                            xml_file,
                            filename='custom'):
        """
        Builds a system from a PDB file containing multiple residues and custom XML files.
        
        :param system_pdb:
            PDB file containing the system or a list of PDB files.
        :param xml_file:
            XML file containing the forcefield parameters or a list of XML files.
        :param filename:
            Base filename for output and intermediate files. Default is 'custom_system'.
        """
        
        # Load the PDB file and format it correctly
        # This method already prints the unique residues info!
        self.load_system_PDB(pdb_file)

        # Load the generated QM region topology and system PDB
        self.pdb = app.PDBFile(pdb_file)
        # Check if the system is periodic from the PDB file
        if self.pdb.topology.getUnitCellDimensions() is not None:
            periodic = True
            msg = 'PBC detected from the PDB file.'
            self.ostream.print_info(msg)
            self.ostream.flush()

        # Combine XML files into a forcefield
        xml_files = [xml_file] if isinstance(xml_file, str) else xml_file
        msg = f'Added XML Files: {xml_files}'
        self.ostream.print_info(msg)
        self.ostream.flush()

        # Combine XML files into a forcefield
        forcefield = app.ForceField(*xml_files, self.parent_ff, self.water_ff)

        # Create the system with appropriate nonbonded settings
        if periodic:
            nonbondedMethod = app.PME  # Particle Mesh Ewald method for periodic systems
            nonbondedCutoff = self.cutoff * unit.nanometer  # Assuming self.cutoff is defined elsewhere appropriately
            constraints = app.HBonds
        else:
            nonbondedMethod = app.NoCutoff  # No cutoff for non-periodic systems
            nonbondedCutoff = None
            constraints = app.HBonds

        # Check if nonbondedCutoff is defined before using it in createSystem
        system_arguments = {
            'nonbondedMethod': nonbondedMethod,
            'constraints': constraints
        }   
        
        if nonbondedCutoff is not None:
            system_arguments['nonbondedCutoff'] = nonbondedCutoff

        self.system = forcefield.createSystem(self.pdb.topology, 
                                              **system_arguments)

        # Save system to XML and PDB for inspection and reuse
        with open(f'{filename}_system.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(self.system))
            msg = f'System parameters written to {filename}_system.xml'
            self.ostream.print_info(msg)
            self.ostream.flush()

        with open(f'{filename}_system.pdb', 'w') as pdb_file:
            app.PDBFile.writeFile(self.pdb.topology, 
                                  self.pdb.positions, 
                                  pdb_file)
            msg = f'System coordinates written to {filename}_system.pdb'
            self.ostream.print_info(msg)
            self.ostream.flush()

        self.phase = 'gas'

    def qmmm_system_from_files(self, 
                               pdb_file, 
                               xml_file, 
                               qm_residue, 
                               ff_gen_qm=None,
                               qm_atoms='all', 
                               filename='qmmm'):
        """
        Builds a QM/MM system from a PDB file containing multiple residues and custom XML files.
        
        :param pdb_file: 
            PDB file containing the system.
        :param xml_file: 
            XML file containing the forcefield parameters. Can be a list of XML files.
        :param qm_residue: 
            Tuple containing the name and number of the QM residue. (e.g. ('MOL', 1))
        :param ff_gen_qm: 
            Optional custom forcefield generator for the QM region. If None, a new one will be created.
        :param qm_atoms:
            Options: 'all' or a list of atom indices for QM region.
        :param filename: 
            Base filename for output and intermediate files.
        """

        # Load the PDB file and format it correctly
        self.load_system_PDB(pdb_file)

        # Extracting the residue name and number for clarity and correct usage
        qm_residue_name, qm_residue_number = qm_residue
        msg = f'Unique Residues: {self.unique_residues}'
        self.ostream.print_info(msg)
        self.ostream.flush()
        qm_residue_index = self.unique_residues.index((qm_residue_name, qm_residue_number))
        msg = f'QM Residue: {qm_residue_name}'
        self.ostream.print_info(msg)
        self.ostream.flush()

        qm_molecule = self.unique_molecules[qm_residue_index]
        self.z_matrix = qm_molecule.get_z_matrix()
        # Generate or use an existing forcefield generator for the QM region
        if ff_gen_qm is None:
            ff_gen_qm = ForceFieldGenerator()
            ff_gen_qm.force_field_data = 'gaff-2.11.dat'
            ff_gen_qm.create_topology(qm_molecule)

        # Determine qm_atoms based on the QM residue
        if qm_atoms == 'all':
            msg = 'Full molecule as QM region'
            self.ostream.print_info(msg)
            qm_atoms = list(range(qm_molecule.number_of_atoms()))
            self._create_QM_residue(ff_gen_qm,
                        qm_atoms, 
                        filename=filename, 
                        residue_name=qm_residue_name)
        elif isinstance(qm_atoms, list):
            if qm_atoms == list(range(qm_molecule.number_of_atoms())): 
                msg = 'Full molecule as QM region'
                self.ostream.print_info(msg)
                self._create_QM_residue(ff_gen_qm,
                        qm_atoms, 
                        filename=filename, 
                        residue_name=qm_residue_name)
            msg = 'QM/MM partition inside the molecule'
            self.ostream.print_info(msg)
            qm_atoms = qm_atoms
            self._create_QM_subregion(ff_gen_qm,
                                      qm_atoms, 
                                      qm_molecule, 
                                      filename,
                                      residue_name=qm_residue_name)
        else:
            raise ValueError('Invalid value for qm_atoms. Please use "all" or a list of atom indices.')

        # Load the generated QM region topology and system PDB
        self.pdb = app.PDBFile(pdb_file)
        # Check if the system is periodic from the PDB file
        if self.pdb.topology.getUnitCellDimensions() is not None:
            periodic = True
            msg = 'PBC detected from the PDB file.'
            self.ostream.print_info(msg)
            self.ostream.flush()

        # Get self.positions from the topology
        self.positions = self.pdb.positions

        # Combine XML files into a forcefield
        xml_files = [f'{filename}.xml'] + ([xml_file] if isinstance(xml_file, str) else xml_file)
        msg = f'Added XML Files: {xml_files}'
        self.ostream.print_info(msg)
        self.ostream.flush()

        # Combine XML files into a forcefield
        forcefield = app.ForceField(*xml_files)

        # Create the system with appropriate nonbonded settings
        if periodic:
            nonbondedMethod = app.PME  # Particle Mesh Ewald method for periodic systems
            nonbondedCutoff = self.cutoff * unit.nanometer  # Assuming self.cutoff is defined elsewhere appropriately
            constraints = app.HBonds
        else:
            nonbondedMethod = app.NoCutoff  # No cutoff for non-periodic systems
            nonbondedCutoff = None
            constraints = app.HBonds

        # Check if nonbondedCutoff is defined before using it in createSystem
        system_arguments = {
            'nonbondedMethod': nonbondedMethod,
            'constraints': constraints
        }   

        if nonbondedCutoff is not None:
            system_arguments['nonbondedCutoff'] = nonbondedCutoff

        self.system = forcefield.createSystem(self.pdb.topology, **system_arguments)

        # Set QM atoms based on the QM residue
        #self.qm_atoms = [atom.index for atom in self.pdb.topology.atoms() if (atom.residue.name == qm_residue_name and atom.residue.id == str(qm_residue_number))]
        self.qm_atoms = qm_atoms
        msg = f'QM Atom indices: {self.qm_atoms}'
        self.ostream.print_info(msg)
        self.ostream.flush()

        # Setting QM/MM system
        self.set_qm_mm_system('periodic' if periodic else 'gas', ff_gen_qm)

        # Adding the stabilizer force to the QM region
        self.qm_stabilizer(ff_gen_qm)

        # Save system to XML and PDB for inspection and reuse
        with open(f'{filename}_system.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(self.system))
            msg = f'System parameters written to {filename}_system.xml'
            self.ostream.print_info(msg)
            self.ostream.flush()

        with open(f'{filename}_system.pdb', 'w') as pdb_file:
            app.PDBFile.writeFile(self.pdb.topology, self.pdb.positions, pdb_file)
            msg = f'System coordinates written to {filename}_system.pdb'
            self.ostream.print_info(msg)
            self.ostream.flush()

        # Correct phase setting
        self.phase = 'gas'

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

    # Simulation methods
    def energy_minimization(self, max_iter=0, tol=10.0):
        """
        Minimizes the energy of the system using the specified parameters.

        Args:
            max_iter (int): Maximum number of iterations for the minimization. Default is 0 (no limit).
            tol (float): Tolerance for the energy minimization, in kJ/mol. Default is 10.0.

        Raises:
            RuntimeError: If the system has not been created prior to the call.

        Returns:
            float: The minimized potential energy of the system.
            str: XYZ format string of the relaxed coordinates.
        """
        if self.system is None:
            raise RuntimeError('System has not been created!')

        # Create an integrator and simulation object
        self.integrator = self._create_integrator()

        self.simulation = app.Simulation(self.modeller.topology if 
                                         self.phase in ['water', 'custom', 'periodic'] else
                                         self.pdb.topology,
                                        self.system, self.integrator)
        
        self.simulation.context.setPositions(self.modeller.positions if 
                                             self.phase in ['water', 'custom', 'periodic'] else
                                             self.pdb.positions)

        # Perform energy minimization
        self.simulation.minimizeEnergy(tolerance=tol * unit.kilojoules_per_mole / unit.nanometer, maxIterations=max_iter)
        
        # Retrieve and process the final state of the system
        state = self.simulation.context.getState(getEnergy=True, getPositions=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        coordinates = np.array(state.getPositions().value_in_unit(unit.nanometer)) * 10  # Convert nm to Angstroms

        # Construct XYZ format string for the coordinates
        xyz = f"{len(self.labels)}\n\n"
        xyz += "\n".join(f"{label} {x} {y} {z}" for label, (x, y, z) in zip(self.labels, coordinates))
        
        return energy, xyz
    
    def steepest_descent(self, max_iter=1000, learning_rate=0.0001, convergence_threshold=1.0e-2):
        """
        Performs a steepest descent minimization on the system.

        :param max_iter: 
            Maximum number of iterations. Default is 1000.
        :param learning_rate:
            Initial learning rate for the optimization. Default is 0.01.
        :param convergence_threshold:
            Convergence threshold for the optimization. Default is 1e-3.

        :return:
            Tuple containing the XYZ format string of the relaxed coordinates and the final energy.
        """
            
        if self.system is None:
            raise RuntimeError('System has not been created!')
        
        # Create an integrator and simulation object
        self.integrator = mm.VerletIntegrator(0.002 * unit.picoseconds)

        self.simulation = app.Simulation(self.modeller.topology if 
                                         self.phase in ['water', 'custom', 'periodic'] else
                                         self.pdb.topology,
                                        self.system, self.integrator)
        
        # Get initial positions
        positions = self.modeller.positions if self.phase in ['water', 'custom', 'periodic'] else self.pdb.positions
        positions = np.array(positions.value_in_unit(unit.nanometer))

        # Define the energy and gradient functions
        def compute_energy_and_gradient(positions):

            self.simulation.context.setPositions(positions)
            state = self.simulation.context.getState(
                getPositions=True,
                getEnergy=True, 
                getForces=True)
            
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            gradient = -1.0 * state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)

            return energy, gradient

        for iteration in range(max_iter):
            energy, gradients = compute_energy_and_gradient(positions)
            print("Gradient norm:", np.linalg.norm(gradients))

            # Optional: Normalize gradients
            norm = np.linalg.norm(gradients)
            if norm > 0:
                gradients /= norm

            # Update positions
            new_positions = positions - learning_rate * gradients

            # Adaptive learning rate adjustment
            if iteration > 0 and energy > previous_energy:
                learning_rate *= 0.5
                print("Reducing learning rate to:", learning_rate)

            previous_energy = energy

            # Check for convergence
            if np.linalg.norm(new_positions - positions) < convergence_threshold:
                print(f"Convergence reached after {iteration+1} iterations.")
                break

            positions = new_positions
            print(f"Iteration {iteration+1}: Energy = {energy}")

        # Once converged, return the final energy and positions
        # The positions shall be written in XYZ format in angstroms
        final_positions = positions * 10
        # Construct XYZ format string for the coordinates
        xyz = f"{len(self.labels)}\n\n"
        xyz += "\n".join(f"{label} {x} {y} {z}" for label, (x, y, z) in zip(self.labels, final_positions))

        return energy, xyz
    
    def steepest_gradient(self):
        """
        Minimize the energy using the steepest descent algorithm.
        Requires additional openmmtools package.
        """

        try:
            from openmmtools.integrators import GradientDescentMinimizationIntegrator
        except ImportError:
            raise ImportError('The openmmtools package is required for this method.')
        
        if self.system is None:
            raise RuntimeError('System has not been created!')

        self.integrator = GradientDescentMinimizationIntegrator()
        self.simulation = app.Simulation(self.modeller.topology if self.phase in ['water', 'custom', 'periodic'] else self.pdb.topology,
                                         self.system, self.integrator)
        self.simulation.context.setPositions(self.modeller.positions if self.phase in ['water', 'custom', 'periodic'] else self.pdb.positions)

        self.simulation.minimizeEnergy()

        state = self.simulation.context.getState(getEnergy=True, getPositions=True)

        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        coordinates = np.array(state.getPositions().value_in_unit(unit.nanometer)) * 10  # Convert nm to Angstroms

        xyz = f"{len(self.labels)}\n\n"
        xyz += "\n".join(f"{label} {x} {y} {z}" for label, (x, y, z) in zip(self.labels, coordinates))

        return energy, xyz
     
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

    def run_md(self, 
               restart_file=None,
               ensemble='NVE',
               temperature=298.15, 
               pressure=1.0, 
               friction=1.0,
               timestep=2.0, 
               nsteps=1000, 
               snapshots=100, 
               traj_file='trajectory.pdb',
               state_file='output.xml',
               save_last_frame=True,
               output_file='output'):
        """
        Runs an MD simulation using OpenMM, storing the trajectory and simulation data.

        :param restart_file:
            Last state of the simulation to restart from. Default is None.
        :param ensemble:
            Type of ensemble. Options are 'NVE', 'NVT', 'NPT'. Default is 'NVE'.
        :param temperature:
            Temperature of the system in Kelvin. Default is 298.15 K.
        :param pressure:
            Pressure of the system in atmospheres. Default is 1.0 atm.
        :param friction:
            Friction coefficient in 1/ps. Default is 1.0.
        :param timestep:
            Timestep of the simulation in femtoseconds. Default is 2.0 fs.
        :param nsteps:
            Number of steps in the simulation. Default is 1000.
        :param snapshots:
            Frequency of snapshots. Default is 100.
        :param traj_file:
            Output file name for the trajectory. Default is 'trajectory.pdb'.
        :param state_file:
            Output file name for the simulation state. Default is 'output.xml'.
        """

        if self.system is None:
            raise RuntimeError('System has not been created!')

        self.ensemble = ensemble
        self.temperature = temperature * unit.kelvin
        self.friction = friction / unit.picosecond
        self.timestep = timestep * unit.femtoseconds
        self.nsteps = nsteps

        self.total_potentials = []
        self.kinetic_energies = []
        self.temperatures = []
        self.total_energies = []


        # Create or update the integrator
        if self.integrator is None:
            new_integrator = self._create_integrator()
        else:
            new_integrator = self.integrator

        self.topology = self.pdb.topology

        self.positions = self.pdb.positions
        
        self.simulation = app.Simulation(self.topology, self.system, new_integrator)
        self.simulation.context.setPositions(self.positions)
        
        # Load the state if a restart file is provided
        if restart_file is not None:
            self.simulation.loadState(restart_file)
            msg = f'Restarting simulation from {restart_file}'
            self.ostream.print_info(msg)
            self.ostream.flush()
        
        else:
            self.simulation.context.setPositions(self.positions)

            # Set initial velocities if the ensemble is NVT or NPT
            if self.ensemble in ['NVT', 'NPT']:
                self.simulation.context.setVelocitiesToTemperature(self.temperature)

        # Minimize the energy before starting the simulation

        # Set up reporting
        save_freq = max(nsteps // snapshots, 1)
        self.simulation.reporters.clear()  
        self.simulation.reporters.append(app.PDBReporter(traj_file, save_freq))

        # Print header
        print('MD Simulation parameters:')
        print('=' * 50)
        print('Ensemble:', ensemble)
        if ensemble in ['NVT', 'NPT']:
            print('Temperature:', temperature, 'K')
            if ensemble == 'NPT':
                print('Pressure:', pressure, 'atm')
        print('Friction:', friction, '1/ps')
        print('Timestep:', timestep, 'fs')
        print('Total simulation time in ns:', nsteps * timestep / 1e6)
        print('=' * 50)


        start_time = time()
        for step in range(nsteps):
            self.simulation.step(1)
            if step % save_freq == 0:
                state = self.simulation.context.getState(getEnergy=True)
                print(f"Step: {step} / {nsteps} Time: {round((step * timestep)/1000, 2)} ps")

                # Potential energy
                pot = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                print('Potential Energy', pot, 'kJ/mol')
                self.total_potentials.append(pot)

                # Kinetic energy
                kin = state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
                print('Kinetic Energy:', kin, 'kJ/mol')
                self.kinetic_energies.append(kin)

                # Temperature
                R_kjmol = 0.00831446261815324
                particles = self.system.getNumParticles() * 1.5
                temp = kin / (particles * R_kjmol)
                self.temperatures.append(temp)
                print('Temperature:', temp, 'K')

                # Total energy
                total = pot + kin
                self.total_energies.append(total)
                print('Total Energy:', total, 'kJ/mol')
                print('-' * 60)


        end_time = time()
        elapsed_time = end_time - start_time
        elapsed_time_days = elapsed_time / (24 * 3600)
        performance = (nsteps * timestep / 1e6) / elapsed_time_days

        print('MD simulation completed!')
        print(f'Number of steps: {nsteps}')
        print(f'Trajectory saved as {traj_file}')
        print('=' * 60)
        print('Simulation Averages:')
        print('=' * 60)
        print('Total Potential Energy:', np.mean(self.total_potentials), '±', np.std(self.total_potentials), 'kJ/mol')
        print('Kinetic Energy:', np.mean(self.kinetic_energies), '±', np.std(self.kinetic_energies), 'kJ/mol')
        print('Temperature:', np.mean(self.temperatures), '±', np.std(self.temperatures), 'K')
        print('Total Energy:', np.mean(self.total_energies), '±', np.std(self.total_energies), 'kJ/mol')
        print('=' * 60)
        print(f'Elapsed time: {int(elapsed_time // 60)} minutes, {int(elapsed_time % 60)} seconds')
        print(f'Performance: {performance:.2f} ns/day')
        print(f'Trajectory saved as {traj_file}')

        # Save the state of the simulation
        self.simulation.saveState(state_file)
        print(f'Simulation state saved as {state_file}')

        if save_last_frame:
            # Save the last frame of the trajectory as last_frame.pdb
            with open('last_frame.pdb', 'w') as f:
                app.PDBFile.writeFile(self.simulation.topology, self.simulation.context.getState(getPositions=True).getPositions(), f)
            print('Last frame saved as last_frame.pdb')

        # Write the output to a file
        self._save_output(output_file)
        self.ostream.print_info(f'Simulation output saved as {output_file}')
        self.ostream.flush()

    def sort_points_with_association(self, points, associated_list):
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
    
    
    def update_settings(self, dynamics_settings, impes_dict=None):
        """
        Updates settings in the ImpesDynamicsDriver.
        :param molecule:
            The Molecule object.
        :param impes_dict:
            The input dictionary of settings for IMPES.
        """

        assert_msg_critical(
        dynamics_settings['qm_driver'] is not None or dynamics_settings['grad_driver'] is not None or dynamics_settings['hess_driver'] is not None,
        'dynamics_settings: Not all drivers (energy, gradient, hessian) are defined please add the missing driver/s in settings.')

        self.qm_driver = dynamics_settings['qm_driver']
        self.grad_driver = dynamics_settings['grad_driver']
        self.hess_driver = dynamics_settings['hess_driver']

        self.qm_driver.ostream.mute()
        self.impes_dict = impes_dict


        #if 'checkpoint_file_name' not in impes_dict:
        
        #    #positions_ang = self.molecule.get_coordinates()
        #    #atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
        #    #qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
        #    #new_molecule = Molecule(qm_atom_labels, positions_ang, units="au")
        #    self.roots = int(dynamics_settings['roots'])
        #    new_molecule = self.molecule
        #    if 'calc_NAC' in dynamics_settings:
        #        self.calc_NAC = dynamics_settings['calc_NAC']
        #    else:
        #        self.calc_NAC = False
        #    self.qm_data_points = []
        #    potential_kjmol = self.compute_energy(new_molecule)
        #    self.im_labels = []
        #    self.qm_energies = []
        #    basename = impes_dict['basename']
        #    self.qm_datafile = f'{basename}.h5'
        #    if self.roots > 1:
        #        self.add_point(new_molecule, ['point_0', 'point_1'], potential_kjmol, self.qm_datafile)
        #    else:
        #        self.add_point(new_molecule, ['point_0'], potential_kjmol, self.qm_datafile)
        #    impes_dict['checkpoint_file_name'] = self.qm_datafile


        if 'print_step' in dynamics_settings:
            self.print_step = int(dynamics_settings['print_step'])

        if 'duration' in dynamics_settings:
            self.duration = float(dynamics_settings['duration'])

        if 'temperature' in dynamics_settings:
            self.temperature = float(dynamics_settings['temperature'])
            self.starting_temperature = float(dynamics_settings['temperature'])

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

        if 'equilibrium_coordinate' in dynamics_settings:
            coords = self.parse_coordinate(
                            dynamics_settings['equilibrium_coordinate'])
            self.define_equilibrium_coordinate(coords)

        # Start the simulation form a given state
        if 'load_system' in dynamics_settings:
            self.load_system = dynamics_settings['load_system']

        if 'trajectory_file' in dynamics_settings:
            self.out_file = dynamics_settings['trajectory_file']
        else:
            self.out_file = 'trajectory.pdb'

        # Determines the ensemble in order to set the correct simulation set_up
        if 'ensemble' in dynamics_settings:
            self.ensemble = dynamics_settings['ensemble']
        
        if 'FF_datafile' in dynamics_settings:
            self.ff_datafile = dynamics_settings['FF_datafile']

        #################################### DATABASE construciton inputs #############################

        if 'checkpoint_file_name' in impes_dict:
            self.qm_datafile = impes_dict['checkpoint_file_name']

        # The desired density around a given starting structure/datapoint
        if 'desired_datapoint_density' in dynamics_settings:
            self.desired_datpoint_density = int(dynamics_settings['desired_datapoint_density'])
        else:
            self.desired_datpoint_density = 30

        # The number of iteration cycles without database expansion after which the description around the reference point has converged
        if 'converged_cycle' in dynamics_settings:
            self.unadded_cycles = int(dynamics_settings['converged_cycle'])
        else:
            self.unadded_cycles = 3

        # Dertermines the number of roots that should be calculated
        if 'roots' in dynamics_settings:
            self.roots = int(dynamics_settings['roots']) + 1
        
        if 'basis_set' in dynamics_settings:
            basis_label = dynamics_settings['basis_set']
            self.basis = MolecularBasis.read(self.molecule, basis_label)


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

        # time step at which QM data point collection should start
        if 'qmc_start' in dynamics_settings:
            self.qmc_start = float(dynamics_settings["qmc_start"])
            # To prevent collecting points before qmc_start,
            # set collect_qm_points to False
            self.collect_qm_points = False

        # time step at which QM data point collection should stop
        if 'qmc_stop' in dynamics_settings:
            self.qmc_stop = float(dynamics_settings["qmc_stop"])

        # index of the excited state (in case of TDDFT QM data points)
        if "excited_state_index" in dynamics_settings:
            self.excited_state_index = int(dynamics_settings["excited_state_index"])

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
        
        if 'symmetry_groups' in impes_dict:
            self.non_core_symmetry_groups = impes_dict['symmetry_groups']

        if 'checkpoint_file_name' not in impes_dict:

            #positions_ang = self.molecule.get_coordinates()
            #atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            #qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            #new_molecule = Molecule(qm_atom_labels, positions_ang, units="au")
            new_molecule = self.molecule
            self.qm_data_points = []
            print('roots', self.roots, self.basis)
            qm_energy, scf_tensors = self.compute_energy(new_molecule, basis=self.basis)
            
            self.im_labels = []
            self.qm_energies = []
            basename = impes_dict['basename']
            self.qm_datafile = f'{basename}.h5'
            label_list = [f'point_{i}' for i in range(self.roots)]
            self.add_point(new_molecule, label_list, qm_energy, self.qm_datafile, self.basis, scf_results=scf_tensors)
            impes_dict['checkpoint_file_name'] = self.qm_datafile

    
    
    def run_qmmm(self, state_file='output.xml'):
        """
        Runs a QM/MM simulation using OpenMM, storing the trajectory and simulation data.

        :param qm_driver:
            QM driver object from VeloxChem.
        :param grad_driver:
            Gradient driver object from VeloxChem.
        :param restart_file:
            Last state of the simulation to restart from. Default is None.
        :param ensemble:
            Type of ensemble. Options are 'NVE', 'NVT', 'NPT'. Default is 'NVE'.
        :param temperature:
            Temperature of the system in Kelvin. Default is 298.15 K.
        :param pressure:
            Pressure of the system in atmospheres. Default is 1.0 atm.
        :param friction:
            Friction coefficient in 1/ps. Default is 1.0.
        :param timestep:
            Timestep of the simulation in femtoseconds. Default is 0.5 fs.
        :param nsteps:
            Number of steps in the simulation. Default is 1000.
        :param snapshots:
            Frequency of snapshots. Default is 100.
        :param traj_file:
            Output file name for the trajectory. Default is 'trajectory.pdb'.
        :param state_file:
            Output file name for the simulation state. Default is 'output.xml'.
        """
        if self.system is None:
            raise RuntimeError('System has not been created!')
        
        self.temperature = self.temperature * unit.kelvin
        self.friction = self.friction / unit.picosecond
        timestep = self.timestep
        self.timestep = timestep * unit.femtoseconds

        # Driver flag 
        if isinstance(self.qm_driver, XtbDriver):
            self.driver_flag = 'XTb Driver'
        elif isinstance(self.qm_driver, ImpesDriver):
            self.driver_flag = 'Impes Driver'

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
            #else:
                #self.simulation.context.setVelocitiesToTemperature(10 * unit.kelvin)

        # There is no minimization step for QM/MM simulations
        # It causes instabilities in the MM region!
        
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

        # load datafile if there is one
        impes_driver = ImpesDriver(self.z_matrix)
        impes_driver.update_settings(self.impes_dict)
        self.im_labels = impes_driver.read_labels()
        print('beginning labels', self.im_labels)
        if self.qm_data_points is None:
           self.qm_data_points = []
           self.qm_energies = []
           for label in self.im_labels:
               qm_data_point = ImpesCoordinates(self.z_matrix)
               qm_data_point.read_hdf5(self.qm_datafile, label)
               self.qm_energies.append(qm_data_point.energy)
               self.qm_data_points.append(qm_data_point)

        self.mover_along_path = 0
        self.adjuster = 0
        self.last_point_added = 0
        self.cycle_iteration = self.unadded_cycles
        
        print('current datapoints around the given starting structure', self.desired_datpoint_density, self.density_around_data_point[0], self.density_around_data_point[1], '\n allowed derivation from the given structure', self.allowed_molecule_deviation, '\n ---------------------------------------')
        self.allowed_molecules = []
        self.coordinates = []
        self.coordinates_xyz = []
        self.gradients = []
        self.velocities = []
        self.velocities_np = []
        self.velocities_np.append(self.simulation.context.getState(getVelocities=True).getVelocities(True))
        # if excited states higher roots, gradient, energy and the current IM-Driver for excited states is initialized
        self.impes_drivers = []
        if self.roots > 1:

            for root in range(self.roots):
                # Dynamically create an attribute name
                attribute_name = f'impes_driver_{root}'

                # Initialize the object
                driver_object = ImpesDriver(self.z_matrix)

                driver_object.update_settings(self.impes_dict)
                # Set the object as an attribute of the instance
                setattr(self, attribute_name, driver_object)

                # Append the object to the list
                self.impes_drivers.append(driver_object)


            self.current_im_choice = self.impes_drivers[1]
            self.current_state = 1
            self.current_gradient = 0
            self.current_energy = 0
        else:
            self.current_state = 0
            driver_object = ImpesDriver(self.z_matrix)
            driver_object.non_core_symmetry_group = self.non_core_symmetry_groups
            driver_object.non_core_structure = self.non_core_symmetry_groups

            driver_object.update_settings(self.impes_dict)
            self.impes_drivers.append(driver_object)
            self.current_im_choice = self.impes_drivers[-1]
    
        self.FFlabels, self.qm_data_points = self.sort_points_with_association(self.im_labels, self.qm_data_points)

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
            R_kjmol = 0.00831446261815324
            particles = self.system.getNumParticles() * 1.5
            temp = kinetic.value_in_unit(unit.kilojoules_per_mole) / (particles * R_kjmol)
            self.temperatures.append(temp)

            # Total energy
            total = pot + kinetic
            self.total_energies.append(total.value_in_unit(unit.kilojoules_per_mole))
            print('scaling factor', self.scaling_factor)
            
            #if step % 10 == 0 and self.scaling_factor > 0.1:
            #    self.scaling_factor -= 0.1
                
            # Information output
            if step % save_freq == 0:
                
                print(f"Step: {step} / {self.nsteps} Time: {round((step * timestep) / 1000, 2)} ps")
                print('Potential Energy QM region:', qm, 'kJ/mol')
                print('Potential Energy MM region:', mm)
                print('QM/MM Interaction Energy:', qm_mm)
                print('Total Potential Energy:', pot)
                print('Kinetic Energy:', kinetic)
                print('Temperature:', temp, 'K')
                print('Total Energy:', total)
                print('-' * 60)
                print('Current Density', self.density_around_data_point[0] / self.roots, '-->', self.desired_datpoint_density, self.unadded_cycles)   

            self.output_file_writer(self.summary_output)
            self.step += 1
            if self.step == 5650:
                print('both')
                # exit()
            self.simulation.step(1)
            if step % 100 == 0 and step != 0:
                self.simulation.saveCheckpoint('checkpoint')
                #self.output_file_writer(self.summary_output)
                #print('cooridnates', simulation.context.getState(getPositions=True).getPositions())

            if step == self.nsteps and self.density_around_data_point[0] / self.roots != self.desired_datpoint_density:
                step = 0
            
            if self.density_around_data_point[0] / self.roots >= self.desired_datpoint_density:
                step = self.nsteps
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

        print(f'QM region parameters written to {filename}.xml')

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
        tree = ET.ElementTree(ForceField)
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
            ForceFieldGenerator object from VeloxChem.
        """

        from openmm import NonbondedForce

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
                if isinstance(force, NonbondedForce):
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

    def calculate_translation_coordinates(self, coordinates, rotation_point=None):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(coordinates, axis=0)
        translated_coordinates = coordinates - center

        return translated_coordinates, center

    def cartesian_just_distance(self, coordinate_1, coordinate_2):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.
           :param data_point:
                ImpesCoordinates object
        """
        # First, translate the cartesian coordinates to zero
        target_coordinates, center_target = self.calculate_translation_coordinates(coordinate_1)
        reference_coordinates, center_reference = (
            self.calculate_translation_coordinates(coordinate_2))
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
            ForceFieldGenerator object from VeloxChem.
        """

        print('scaling factor', self.scaling_factor)
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
        
        for root in range(self.roots):
            self.impes_drivers[root].compute(new_molecule, self.qm_data_points[root::self.roots], None, self.im_labels[root::self.roots], self.calc_NAC)
            # qm_energy, scf_tensors = self.compute_energy(new_molecule, self.basis)
            # print('\n\n ENergy difference ', qm_energy[0], self.impes_drivers[root].impes_coordinate.energy, abs(qm_energy[0] - self.impes_drivers[root].impes_coordinate.energy) * hartree_in_kcalpermol(), '\n\n')

        # determine the energy difference between difference states and if NAC determine non-adiabatic coupling matrix
        transitions = []
        if self.roots > 1:
            for root_1 in range(0, self.roots - 1):
                for root_2 in range(root_1 + 1, self.roots):

                    potential_kjmol = self.impes_drivers[root_1].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
                    potential_kjmol_2 = self.impes_drivers[root_2].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
                    transitions.append((root_1, root_2, potential_kjmol - potential_kjmol_2))
                    print(f'compare the energies between roots: {root_1} -> {root_2}', potential_kjmol - potential_kjmol_2)

                    if self.calc_NAC == True and np.linalg.norm(self.velocities_np[-1]) > 0.0:
                        current_NAC = self.impes_drivers[root_1].impes_coordinate.NAC.flatten()
                        current_velocitites = self.velocities_np[-1].flatten() * 4.566180e-4

                        hopping_potential = np.exp(-abs((np.pi/(4)) * (( self.impes_drivers[root_2].impes_coordinate.energy - self.impes_drivers[root_1].impes_coordinate.energy ) / np.linalg.multi_dot([current_NAC, current_velocitites]))))
                        print('#######################', '\n\n', hopping_potential, potential_kjmol_2 - potential_kjmol, '\n\n', '#######################')

                    if abs(potential_kjmol_2 - potential_kjmol) < 20:
                        # Choose a random integer between 0 and 1
                        random_integer = random.randint(0, 1)
                        if random_integer == 1 and self.current_im_choice == self.impes_drivers[root_1]:
                            self.current_im_choice = self.impes_drivers[root_2]
                            self.current_gradient = self.impes_drivers[root_2].impes_coordinate.gradient
                            self.current_energy = potential_kjmol_2
                            self.current_state = root_2
                        elif random_integer == 1 and self.current_im_choice == self.impes_drivers[root_2]:
                            self.current_im_choice = self.impes_drivers[root_1]
                            self.current_gradient = self.impes_drivers[root_1].impes_coordinate.gradient
                            self.current_energy = potential_kjmol
                            self.current_state = root_1
                        break

                    elif self.current_im_choice == self.impes_drivers[root_1]:
                        self.current_gradient = self.impes_drivers[root_1].impes_coordinate.gradient
                        self.current_energy = potential_kjmol

                    elif self.current_im_choice == self.impes_drivers[root_2]:
                        self.current_gradient = self.impes_drivers[root_2].impes_coordinate.gradient
                        self.current_energy = potential_kjmol_2

        else:
            potential_kjmol = self.impes_drivers[-1].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
            self.current_gradient = self.impes_drivers[-1].impes_coordinate.gradient
            self.current_energy = potential_kjmol
                        
            # Potential energy is in Hartree, convert to kJ/mol
            
        for initial_state, final_state, energy_gap in transitions:
            self.energy_gabs[self.step] = ((initial_state, final_state), energy_gap)
        print('Following the current state:', self.current_energy)
        if self.step == self.excitation_step and self.roots > 1 and self.current_im_choice == self.impes_drivers[0]:
            self.current_im_choice = self.impes_drivers[1]
            self.current_gradient = self.impes_drivers[1].impes_coordinate.gradient
            self.current_energy = self.impes_drivers[1].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
            print('the Interpolation driver has been switched')

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
        if self.density_around_data_point[1] is not None:
            current_dihedral = new_molecule.get_dihedral_in_degrees(self.density_around_data_point[1])
            lower, upper = self.allowed_molecule_deviation

            # Case 1: If boundaries do not wrap (e.g., [-60, 60])
            if lower < upper:
                allowed = lower <= current_dihedral <= upper
            else:
                # Case 2: If boundaries wrap around (e.g., [120, -120])
                allowed = current_dihedral >= lower or current_dihedral <= upper
        
        if allowed:
            if self.skipping_value == 0:
                min_sum_of_sqaures = np.inf
                for i, dp in enumerate(self.qm_data_points):
                        
                    _, _, distance_core_vector = self.calculate_distance_to_ref(new_molecule.get_coordinates_in_bohr(), dp.cartesian_coordinates)
                    sum_of_sqaures = 0.0
                    modes = dp.normal_modes_in_cartesian
                    for k, mode in enumerate(modes):
                        coeffcicent = np.dot(mode.ravel(), distance_core_vector.ravel())
                        weight = 1.0 / np.sqrt(abs(dp.normal_mode_frequencies[k])) * np.sqrt(abs(dp.normal_mode_frequencies[0]))
                        ratio = (coeffcicent / dp.normal_mode_displacement_coefficients[k]) * weight
                            
                        sum_of_sqaures += ratio**2
                        print('coeffcicent: ', k, sum_of_sqaures, min_sum_of_sqaures, weight)
                        
                    if abs(sum_of_sqaures) < min_sum_of_sqaures:
                        min_sum_of_sqaures = abs(sum_of_sqaures)
                
                if min_sum_of_sqaures > 1.0:
                    self.add_a_point = True
            
            else:
                self.skipping_value -= 1
            
            openmm_coordinate = context.getState(getPositions=True).getPositions()
            self.coordinates.append(openmm_coordinate)
            self.velocities.append(context.getState(getVelocities=True).getVelocities())
            self.gradients.append(gradient)
            self.point_checker += 1 

            if self.add_a_point == True and self.step > self.collect_qm_points or self.check_a_point == True and self.step > self.collect_qm_points:
                print('no point correlation ')
                self.point_correlation_check(new_molecule, self.basis)
            if self.point_checker == 0:            
                
                self.point_checker += 1

                context.setPositions(self.coordinates[0])
                self.coordinates = self.coordinates[:1]
                context.setVelocities(self.velocities[0])
                self.velocities = self.velocities[:1]
                
                if self.roots > 1:
                    self.current_im_choice = self.impes_drivers[1]

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
            self.molecules.append(new_molecule)
            for i, atom_idx in enumerate(self.qm_atoms):
                self.system.getForce(self.qm_force_index).setParticleParameters(i, atom_idx, force[i])
            self.system.getForce(self.qm_force_index).updateParametersInContext(context)
        
        else:
            context.setPositions(self.coordinates[0])
            self.coordinates = self.coordinates[:1]
            context.setVelocities(self.velocities[0])
            self.velocities = self.velocities[:1]
            new_positions = context.getState(getPositions=True).getPositions()
            qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])
            positions_ang = (qm_positions) * 10
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            gradient_2 = self.update_gradient(qm_positions)
            force = -np.array(gradient_2) * conversion_factor
            
            self.molecules.append(new_molecule)
            self.coordinates_xyz.append(qm_positions * 10)
            for i, atom_idx in enumerate(self.qm_atoms):
                self.system.getForce(self.qm_force_index).setParticleParameters(i, atom_idx, force[i])
            self.system.getForce(self.qm_force_index).updateParametersInContext(context)
    
    
    def get_qm_potential_energy(self):
        """
        Returns the potential energy of the QM region.

        Args:
            context: The OpenMM context object.
        Returns:
            The potential energy of the QM region.
        """

        #positions = context.getState(getPositions=True).getPositions()
        #qm_positions = np.array([positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])
        potential_energy = self.current_energy

        return potential_energy
    
    
    
        ####################################################################
    ################ Functions to expand the database ##################
    ####################################################################

    def point_correlation_check(self, molecule, basis=None):
        """ Takes the current point on the PES and checks if a QM-calculation
            is nessecary to perform in order to add a new point into the
            QM-databse.

            :param molecule:
                The molecule corresponding to the current conformation
                in the IM dynamics.
        """

        # Check energies and conformations and determine if a point should
        # be added or if a QM energy calculation should be performed
        if not(self.add_a_point):
            # Atom labels for the QM region
            for root in range(self.roots):
                    energy_diff_list = [abs(self.qm_energies[root::self.roots]
                        - self.current_im_choice.impes_coordinate.energy)]

                    if any(np.any(energy_diff * hartree_in_kcalpermol() < self.energy_threshold) for energy_diff in energy_diff_list):

                        for qm_data_point in self.qm_data_points[root::self.roots]:
                            length_vectors = (
                                self.current_im_choice.impes_coordinate.cartesian_distance_vector(
                            qm_data_point))
                        #print('here si the length ve vector', np.linalg.norm(length_vectors))
                        #TODO: Check if this is 1 Angstrom
                        if np.linalg.norm(length_vectors) < 1.0:
                            self.add_a_point = False
                            break
                        else:
                            print('######## the structure is to far away from the references ###########')
                            self.add_a_point = True
                    else:
                        print('############ The Energy difference is to big #############')
                        self.add_a_point = True
        # TODO: this is actually always True... because it is a condition for
        # calling this routine -> no because self.add_a_point can be True and then turn False.
        qm_energy = 0
        if self.add_a_point:
            print('############# Energy is QM claculated ############')
            qm_energy, scf_tensors = self.compute_energy(molecule, basis)
            energy_difference_list = []
            for root in range(self.roots):
                energy_difference_list.append(abs(qm_energy[root] - self.impes_drivers[root].impes_coordinate.energy))
                print('energy differences', energy_difference_list[-1] * hartree_in_kcalpermol())
                self.skipping_value = min(round(abs(self.energy_threshold / (energy_difference_list[-1] * hartree_in_kcalpermol())**2)), 10)
                print('skipping value is here ', self.skipping_value) 

            if any(energy_difference * hartree_in_kcalpermol() > self.energy_threshold for energy_difference in energy_difference_list):
                self.add_a_point = True
            else:
                self.allowed_molecules.append(molecule.get_coordinates_in_bohr())
                self.add_a_point = False
        if self.add_a_point:
            print('✨ A point is added! ✨', self.point_checker)
            labels = []
            for root in range(self.roots):
                labels.append("point_{0}".format((len(self.im_labels) + root)))
            self.add_point(molecule, labels, qm_energy, self.qm_datafile, basis, scf_results=scf_tensors)
            self.last_point_added = self.point_checker - 1
            self.point_checker = 0
            self.point_adding_molecule[self.step] = (molecule, qm_energy, labels)


    def add_point(self, molecule, labels, energy, filename, basis=None, scf_results=None):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param label:
                the label for the new point to be added.
            :param basis:
                the basis set (if required)
        """

        if self.qm_driver is None:
            raise ValueError("No energy driver defined.")
        if self.grad_driver is None:
            raise ValueError("No gradient driver defined.")
        if self.hess_driver is None:
            raise ValueError("No Hessian driver defined.")

        energy = energy
        gradient = None
        hessian = None

       
        gradient = self.compute_gradient(molecule, basis, scf_results)
        hessian = self.compute_hessian(molecule, basis)

        impes_coordinate = ImpesCoordinates(self.z_matrix)

        impes_coordinate.update_settings(self.impes_dict)
        
        qm_points_list = []
        for root in range(self.roots):
            impes_coordinate = ImpesCoordinates(self.z_matrix)
            impes_coordinate.update_settings(self.impes_dict)
            impes_coordinate.cartesian_coordinates = molecule.get_coordinates_in_bohr()

            impes_coordinate.energy = energy[root]
            impes_coordinate.gradient = gradient[root]
            impes_coordinate.hessian = hessian[root]
            impes_coordinate.transform_gradient_and_hessian()
            # impes_coordinate.normal_modes = normal_modes_relevant
            # impes_coordinate.coefficient_displacement_thresholds = 

            if self.impes_drivers is not None:
                self.impes_drivers[root].impes_coordinate.gradient = gradient[root]

                if self.impes_drivers[root] == self.current_im_choice:
                    self.current_im_choice.impes_coordinate.gradient = gradient[root]

            qm_points_list.append(impes_coordinate)
        

        ############################ Implement the Mode analysis for each data_point #########################################
        
        natoms = molecule.number_of_atoms()
        elem = molecule.get_labels()
        coords = molecule.get_coordinates_in_bohr().reshape(natoms * 3)
        

        vib_frequencies, normal_modes_vec, gibbs_energy = (
            geometric.normal_modes.frequency_analysis(
                coords,
                hessian[0],
                elem,
                energy=energy[0],
                temperature=self.starting_temperature,
                pressure=self.pressure,
                outfnm=f'vibrational_point_{energy[0]}',
                normalized=False))
        print('normal modes', normal_modes_vec[0], vib_frequencies)        
        
        # eigen_val, eigen_vec = np.linalg.eigh(internal_hessians)
        # reciprocal_sqrt_mass_vector = 1 / np.sqrt(mass_vector)

        # # Create a diagonal matrix with the reciprocal values
        # reciprocal_sqrt_mass_diagonal = np.diag(reciprocal_sqrt_mass_vector)
        
        # # matrix = np.zeros_like(internal_B_matrix.T)
        # # for i in range(internal_B_matrix.shape[1]):
        # #     for j in range(eigen_vec.shape[1]):
        # #         for n in range(internal_B_matrix.shape[0]):
                
        # #           matrix[i, j] += reciprocal_sqrt_mass_vector[i] * internal_B_matrix.T[i, n] * eigen_vec[n, j]

        # l_cart = np.linalg.multi_dot([reciprocal_sqrt_mass_diagonal, internal_B_matrix.T, eigen_vec])
        
        # print('shapes', internal_B_matrix.shape, eigen_val.shape, l_cart.shape, labels)
        
        # # l_cart = l_cart.reshape(-1, 3)
        # print('\n\n', l_cart, l_cart[:, 0])
        # exit()
        # frequencies_cm1 = 130.3 * np.sqrt(eigen_val)

        # print('frequencies', frequencies_cm1, self.non_core_symmetry_groups)
        original_coordinates = molecule.get_coordinates_in_bohr().copy()

        # Now 'mask' tells you which combos of (qi, qj) keep the energy diff <= 0.7
        # You might visualize it or store it.
        
        normal_modes_relevant = [(freq, normal_modes_vec[i]) for i, freq in enumerate(vib_frequencies) if abs(freq) <= 1000]
        normal_modes_h5_format = [normal_modes_vec[i] for i, freq in enumerate(vib_frequencies) if abs(freq) <= 1000]
        relevant_freqs = vib_frequencies[:len(normal_modes_relevant)]
        if len(normal_modes_relevant) == 0:
            normal_modes_relevant = [(vib_frequencies[0], normal_modes_vec[0])]
            normal_modes_h5_format = [normal_modes_vec[0]]
        self.normal_modes_vec[labels[0]] = normal_modes_relevant
        self.distplacement_threshold[labels[0]] = []
        self.distance_threshold_dict[labels[0]] = 0
        single_amplitude_threshold = []

        sign_swapp = [1.0, -1.0]
        for i, (freq, normal_mode) in enumerate(normal_modes_relevant):
            amplitudes_sign = []
            for sign in sign_swapp:
                print(f'\n Normal Mode {i} with Current frequency {freq} \n\n')
                current_energy_diff = 0
                displaced_coordinates = original_coordinates.copy()
                displaced_molecules = []
                step = sign * 0.05
                amplitude = 0.0
                max_threshold = 1.0
                while True:

                    displcement_matrix = normal_mode
                    amplitude += step
                    displaced_coordinates = (
                        original_coordinates.reshape(-1)
                        + amplitude * displcement_matrix.ravel()
                    )
                    current_displaced_coordinates_diff = displaced_coordinates - original_coordinates.reshape(natoms * 3)

                    coefficicent = np.dot(normal_mode.ravel(), current_displaced_coordinates_diff.ravel())

                    displaced_molecule = Molecule(molecule.get_labels(), displaced_coordinates.reshape(-1, 3), units='bohr')

                    displaced_molecules.append(displaced_molecule)
                    impes_driver = ImpesDriver(self.z_matrix)
                    impes_driver.update_settings(self.impes_dict)

                    impes_driver.distance_thrsh = 1.0
                    impes_driver.compute(displaced_molecule, qm_points_list, None, [labels[0]])
                    qm_energy, _ = self.compute_energy(displaced_molecule, basis)
                    
                    current_energy_diff = abs(qm_energy[0] - impes_driver.impes_coordinate.energy) * 627.509474

                    print('Here is the energy_diff', current_energy_diff, amplitude)

                    if current_energy_diff >= max_threshold:
                        # we overshot, so backtrack by one step
                        amplitude -= step
                        # reduce step
                        step *= 0.5
                        if abs(step) < 1e-4:
                            # we can't refine any further
                            break
                    else:
                        # if we are well below the threshold, maybe we let the step grow a bit
                        if current_energy_diff < 0.5 * max_threshold:
                            step = min(step * 1.2, 0.2)  # don't exceed some maximum step

                amplitudes_sign.append(coefficicent)
            print("Amplitude for mode i is about", amplitudes_sign, current_energy_diff)

            self.distplacement_threshold[labels[0]].append(min(amplitudes_sign))
            single_amplitude_threshold.append(min(amplitudes_sign))
        print('DISPLACEMENT Threshold', self.distplacement_threshold, single_amplitude_threshold)
        old_single_amplitude_thrsh = single_amplitude_threshold.copy()
        print('Here is the old thresholds', old_single_amplitude_thrsh)

        if len(normal_modes_relevant) != 1:
            def two_mode_scan(molecule, original_coordinates, normal_modes_vec, 
                    i, j, qi_range, qj_range, qm_datapoints, labels):
                """
                Scans in a 2D grid along mode i and j. Returns a 2D array of energy diffs
                and a boolean mask of which points are <= threshold.
                """

                def boundary_to_middle_order(range_values):
                    """Generate indices starting from both ends and moving towards the center."""
                    n = len(range_values)
                    left, right = 0, n - 1
                    ordered_indices = []

                    while left <= right:
                        if left != right:
                            ordered_indices.extend([left, right])
                        else:
                            ordered_indices.append(left)
                        left += 1
                        right -= 1

                    return ordered_indices

                def compute_skip(diff, threshold, default_step=1, max_skip=10):
                    """
                    Compute a skip factor based on the difference between the computed energy diff and the threshold.
                    """
                    if diff <= threshold:
                        return default_step  # already acceptable; just move one step
                    # The further diff is above the threshold, the larger the skip.
                    skip = int((diff - threshold) * 3.0)

                    print('compute skipping function:', min(max(default_step, skip), max_skip))
                    return min(max(default_step, skip), max_skip)
                
                def check_boundary(molecule, original_coordinates, mode_i, mode_j, qi, qj, qm_datapoints, labels, threshold=0.7):
                    """
                    Compute the energy difference for the displacement defined by qi_range[idx_i]
                    and qj_range[idx_j] and return a tuple:
                    (satisfied, energy_diff)
                    where satisfied is True if the difference is below the threshold.
                    """

                    # Build the combined displacement:
                    combined_displacement = qi * mode_i + qj * mode_j
                    displaced_coords = original_coordinates.ravel() + combined_displacement
                    # Create a new molecule with the displaced coordinates
                    displaced_molecule = Molecule(
                        molecule.get_labels(), 
                        displaced_coords.reshape(-1, 3), 
                        units='bohr'
                    )

                    impes_driver = ImpesDriver(self.z_matrix)
                    impes_driver.update_settings(self.impes_dict)
                    impes_driver.distance_thrsh = 1.0
                    # Compute energies (using your driver and energy function)
                    impes_driver.compute(displaced_molecule, qm_datapoints, None, [labels[0]])
                    qm_energy, scf_tensors = self.compute_energy(displaced_molecule, basis)
                    current_energy_diff = abs(qm_energy[0] - impes_driver.impes_coordinate.energy) * 627.509474
                    # Store the energy difference
                    energy_diffs[idx_i, idx_j] = current_energy_diff
                    # Return True if within threshold, else False
                    print('Check of energy diff', current_energy_diff < threshold)
                    return current_energy_diff < threshold, current_energy_diff
                
                
                # Prepare storage
                energy_diffs = np.zeros((len(qi_range), len(qj_range)))
                within_threshold = np.ones((len(qi_range), len(qj_range)), dtype=bool)

                mode_i = normal_modes_vec[i][1].ravel()
                mode_j = normal_modes_vec[j][1].ravel()


                qi_indices = boundary_to_middle_order(qi_range)
                qj_indices = boundary_to_middle_order(qj_range)

                # --- First, try the outer boundaries. ---
                #
                # Our idea is to test the “edge” values:
                # for qi, try the first and last element; for qj, try the first and last element.
                # If one edge (say, qi = first element) yields a value below threshold for one of the qj edges,
                # then we “accept” that edge. Then we check the opposite edge (qi = last element).
                #
                # If both opposite edges (for some qj choice) are within the threshold,
                # we can break out early and skip the full grid evaluation.

                # --- If the outer boundaries did not both pass, try “inward” boundaries. ---
                #
                # For example, if one of the outer qi values failed, you might want to try the next candidate.
                # In the code below, we iterate over all qi_indices and qj_indices,
                # but we “skip” calculations if we have already found that the boundaries (on both ends)
                # are OK. (You can modify this logic to “start” at the outer edge and move inward.)
                
                number_of_iterations = 0

                for idx_i in range(len(qi_indices) - 1):
                    print(idx_i, qi_indices[idx_i])
                    qi = qi_range[qi_indices[idx_i]]
                    qi_1 = qi_range[qi_indices[idx_i + 1]]
                    number_of_iterations += 1
                    
                    idx_j = 0
                    while idx_j < len(qj_indices):
                        
                        qj = qj_range[qj_indices[idx_j]]
                        qj_1 = qj_range[qj_indices[idx_j + 1]]
                        number_of_iterations += 1
                        boundaries_satisfied = False

                        # Define indices for the outer boundaries in qi and qj:
                        qi_outer = [qi, qi_1]
                        qj_outer = [qj, qj_1]
                        
                        print('Here is qi_outer and qj_outer', qi_outer, qj_outer)
                        # We will store whether each qi boundary passes (for at least one of the two qj values)
                        qi_boundary_ok = {}
                        qj_boundary_ok = {}
                        counter = 0
                        diff_list = []
                        for i, value_i in enumerate(qi_outer):
                            # For this qi boundary, check the two outer qj values.
                            qi_ok = False
                            qj_ok_list = []
                            for j, value_j in enumerate(qj_outer):
                                qj_ok = False
                                ok, diff = check_boundary(molecule, original_coordinates, mode_i, mode_j, value_i, value_j, qm_datapoints, labels, threshold=0.8)
                                if ok:
                                    # Once one boundary (for this qi) is within threshold, record it.
                                    qi_ok = True
                                    # Optionally, print or log:
                                    print(f"Boundary at qi index {value_i} and qj index {value_j} is OK (energy diff = {diff:.3f}).")
                                    print('Value of qi and qj',np.where(qi_range == value_i), np.where(qj_range == value_j))
                                    qj_ok = True
                                    qj_ok_list.append(qj_ok)

                                else:
                                    # Mark the cell as failing (if you use the within_threshold array):
                                    print(f"Boundary at qi index {value_i} and qj index {value_j} is NOT OK (energy diff = {diff:.3f}).")
                                    print('Value of qi and qj',np.where(qi_range == value_i), np.where(qj_range == value_j))
                                    within_threshold[np.where(qi_range == value_i), np.where(qj_range == value_j)] = 0

                                    qj_ok_list.append(qj_ok)
                                diff_list.append(diff)
                                qj_boundary_ok[counter] = qj_ok
                                counter += 1
                            qi_boundary_ok[i] = qi_ok
                        

                        if all(qi_boundary_ok.values()) and all(qj_boundary_ok.values()) and idx_i == idx_j:
                            print('Number of iterations', number_of_iterations)
                            return energy_diffs, within_threshold

                        # If both outer qi boundaries are OK then we are done.
                        if all(qi_boundary_ok.values()) and all(qj_boundary_ok.values()):
                            boundaries_satisfied = True
                        
                        if all(qj_ok_list):
                            boundaries_satisfied = True
                        if boundaries_satisfied:
                            print("Early termination: outer boundaries satisfied. Skipping full grid scan.")
                            print('Number of iterations', number_of_iterations)
                            break

                        step = compute_skip(min(diff_list), threshold=0.8)
                        step_pairs = step * 2
                        # For all the indices that will be skipped, mark them as False.
                        # (You may want to store a placeholder energy diff; here we simply use the current diff.)
                        for skip_offset in range(2, step_pairs):
                            print('HERE IS THE OFFSET:', idx_j + skip_offset)
                            if idx_j + skip_offset < len(qj_indices):
                                skipped_qj_pos = qj_indices[idx_j + skip_offset]
                                within_threshold[idx_i, skipped_qj_pos] = False
                                energy_diffs[idx_i, skipped_qj_pos] = diff  # or set to a flag value if desired
                                print(f"Skipping qi={qi} and qj={qj_range[skipped_qj_pos]} (marked as False).")
                        
                        # Jump ahead by the computed step size.
                        old_idx = idx_j
                        idx_j = idx_j + step_pairs

                        print('Here is step as is should be', step, idx_j, old_idx)
                        
                return energy_diffs, within_threshold

                                
            print('\n\n\n HERE THE CALCULATION SECTION IS STARTING \n\n\n')

            # Example usage in your code:
            
            def generate_dynamic_step_grid(max_range, min_step_fraction=0.5, boundary_fraction=0.5):
                """
                Generate an adaptive grid with step size based on the boundaries.

                Args:
                - max_range (float): Absolute boundary for the range.
                - min_step_fraction (float): Minimum step size as a fraction of the boundary step.
                - boundary_fraction (float): Region beyond which finer steps should apply.

                Returns:
                - np.array: Array of grid points.
                """
                base_step = 0.1 * max_range  # Base step proportional to the range size
                boundary_region = boundary_fraction * max_range  # Define where finer steps should apply
                min_step = min_step_fraction * base_step  # Minimum step size near boundaries

                grid_points = []
                x = -max_range

                while x <= max_range:
                    grid_points.append(x)

                    # Adapt step size based on position
                    if abs(x) < boundary_region:
                        step_size = base_step  # Coarser steps in the central region
                    else:
                        step_size = min_step  # Finer steps near the boundary

                    x += step_size

                return np.array(grid_points)
            
            def select_threshold_from_grid(success_rates, grid, cutoff=0.6):
                """
                Given an array of success rates (one per row or column) and the corresponding grid,
                select the grid value from the candidate index that:
                (a) has a success rate above the cutoff, and 
                (b) is as close as possible to the center of the grid.
                
                If no candidate exceeds the cutoff, return the middle value of the grid.
                
                Args:
                success_rates (np.array): 1D array of success rates.
                grid (np.array): 1D array of grid values.
                cutoff (float): The minimum success rate required.
                
                Returns:
                float: The selected grid value.
                """
                # Find indices where the success rate is at or above the cutoff
                candidate_indices = np.where(success_rates >= cutoff)[0]
                print('GRID IS HERE', grid, grid[-1], candidate_indices)
                if candidate_indices.size == 0:
                    # No candidate meets the criterion; return the grid’s midpoint.
                    return grid[-1]
                else:
                    # Define the "middle" index of the grid.
                    center_index = len(grid) // 2
                    # Pick the candidate that is closest to this center.
                    best_candidate = candidate_indices[np.argmax(np.abs(candidate_indices - center_index))]
                    if best_candidate == grid.shape[0]:
                        best_candidate -= 1

                    return grid[best_candidate]
                
            for i in range(len(normal_modes_relevant)):
                for j in range(len(normal_modes_relevant)):
                    
                    if j <= i:

                        continue

                    qi_range = generate_dynamic_step_grid(abs(self.distplacement_threshold[labels[0]][i])) # for instance
                    qj_range = generate_dynamic_step_grid(abs(self.distplacement_threshold[labels[0]][j]))

                    print('QI RANGE', qi_range, qj_range)

                    energy_diffs, mask = two_mode_scan(
                        molecule, 
                        original_coordinates, 
                        normal_modes_relevant, 
                        i, j, 
                        qi_range, 
                        qj_range,
                        qm_points_list,
                        labels, 
                    )
                    
                    print('MASK: \n\n', mask)
                    row_success = mask.mean(axis=0)  # fraction of Trues in each row

                    # Similarly, compute for each column (q_j candidates)
                    col_success = mask.mean(axis=1)

                    new_threshold_qi = select_threshold_from_grid(row_success, qi_range, cutoff=0.6)
                    # And similarly for q_j:
                    new_threshold_qj = select_threshold_from_grid(col_success, qj_range, cutoff=0.6)

                    print(f'\n\n SHOW THRESHOLDS {new_threshold_qi, new_threshold_qj} \n\n {self.distplacement_threshold[labels[0]][i], self.distplacement_threshold[labels[0]][j]} \n\n {col_success, row_success} \n\n')
                    
                    if abs(new_threshold_qi) < abs(self.distplacement_threshold[labels[0]][i]):    
                        self.distplacement_threshold[labels[0]][i] = new_threshold_qi
                        single_amplitude_threshold[i] = new_threshold_qi
                    if abs(new_threshold_qj) < abs(self.distplacement_threshold[labels[0]][j]):
                        self.distplacement_threshold[labels[0]][j] = new_threshold_qj
                        single_amplitude_threshold[j] = new_threshold_qj
            
            print(f'\n\n SHOW FINAL THRESHOLDS {new_threshold_qi, new_threshold_qj} \n\n {self.distplacement_threshold[labels[0]]}')

        for root, qm_datapoint in enumerate(qm_points_list):
            
            qm_datapoint.normal_modes_in_cartesian = normal_modes_h5_format
            qm_datapoint.normal_mode_displacement_coefficients = self.distplacement_threshold[labels[root]]
            qm_datapoint.normal_mode_frequencies = relevant_freqs

            qm_datapoint.write_hdf5(filename, labels[root])
            self.im_labels.append(labels[root])
            self.qm_energies.append(energy[root])
            self.qm_energies.append(qm_datapoint.energy)
            self.qm_data_points.append(qm_datapoint)

        self.density_around_data_point[0] += self.roots

    def compute_energy(self, molecule, basis=None):
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
        if isinstance(self.qm_driver, XtbDriver):

            self.qm_driver.compute(molecule)
            qm_energy = self.qm_driver.get_energy()
            qm_energy = np.array([qm_energy])
            print('qm_energy', qm_energy, qm_energy[0])

        # restricted SCF
        elif isinstance(self.qm_driver, ScfRestrictedDriver):
            self.qm_driver.ostream.mute()
            scf_tensors = self.qm_driver.compute(molecule, basis)
            qm_energy = self.qm_driver.scf_energy
            qm_energy = np.array([qm_energy])
            print('qm_energy', qm_energy)
            self.qm_driver.ostream.unmute()

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

    def compute_gradient(self, molecule, basis=None, scf_results=None):
        """ Computes the QM gradient using self.grad_driver.

            :param molecule:
                The molecule.
            :param basis:
                The basis set.

            :returns the QM gradient.
        """

        qm_gradient = None

        if isinstance(self.grad_driver, XtbGradientDriver):
            self.grad_driver.ostream.mute()
            self.grad_driver.compute(molecule)
            self.grad_driver.ostream.unmute()
            qm_gradient = self.grad_driver.gradient
            qm_gradient = np.array([qm_gradient])

        elif isinstance(self.grad_driver, ScfGradientDriver):
            self.grad_driver.ostream.mute()
            self.grad_driver.compute(molecule, basis, scf_results)
            qm_gradient = self.grad_driver.gradient
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
        return qm_gradient

    # TODO: mute outside to save time?
    def compute_hessian(self, molecule, basis=None):
        """ Computes the QM Hessian using self.hess_driver.

            :param molecule:
                The molecule.
            :param basis:
                The basis set.

            :returns the QM Hessian matrix.
        """

        qm_hessian = None

        if isinstance(self.hess_driver, XtbHessianDriver):
            self.hess_driver.ostream.mute()
            self.hess_driver.compute(molecule)
            qm_hessian = self.hess_driver.hessian
            self.hess_driver.ostream.unmute()
            qm_hessian = np.array([qm_hessian])

        elif isinstance(self.hess_driver, ScfHessianDriver):
            # self.hess_driver.ostream.mute()
            self.hess_driver.compute(molecule, basis)
            qm_hessian = self.hess_driver.hessian
            qm_hessian = np.array([qm_hessian])
            # self.hess_driver.ostream.unmute()

        elif isinstance(self.hesian_driver, TddftHessianDriver):
            self.hess_driver.ostream.mute()
            self.hess_driver.compute(molecule, basis)
            qm_hessian = self.hess_driver.hessian
            self.hess_driver.ostream.unmute()

        if qm_hessian is None:
            error_txt = "Could not compute the QM Hessian. "
            error_txt += "Please define a QM Hessian driver."
            raise ValueError(error_txt)

        print('hessina', qm_hessian)
        return qm_hessian

    def adiabatic_transformation(self, energies, gradients, hessians, NACs, NACs_deriv, ADTmatrix):

        # section for the energy in order to transform adiabatic energies into the diabatic basis

        diabatic_energies = np.linalg.multi_dot([ADTmatrix.T, energies, ADTmatrix])

        # section for the gradient in order to transform the adiabatic gradients into the dabatic basis
        # NACs structure is a list probably F = list([0, 1], [0, 2], ..., [1, 0], [1, 2]) assuming that for the transition densities the NAC of [0, 1] == [1, 0]

        nstates = len(energies)
        natoms = len(NACs[0].shape[0])

        adiabatic_Hamilton = np.diag(energies)
        NAC_Matrix = np.empty(adiabatic_Hamilton.shape, dtype=object)
        gradients_matrix = np.empty(adiabatic_Hamilton.shape, dtype=object)
        diabatic_gradients = np.zeros((nstates, nstates, natoms, 3))
        zero_matrix = np.zeros_like(NACs[0])

        for root_1 in range(len(energies)):
            for root_2 in range(len(energies)):
                if root_1 != root_2:
                    gradients_matrix[root_1, root_2] = zero_matrix
                    NAC_Matrix[root_1, root_2] = NACs[root_1 * len(energies) + root_2]

                else:
                    gradients_matrix[root_1, root_2] = gradients[root_2]
                    NAC_Matrix[root_1, root_2] = zero_matrix


        for atom in range(len(self.molecule.get_labels())):
            for direction in range(3):

                diab_sub_gradient = gradients_matrix[:, :, atom, direction] + NAC_Matrix[:, :, atom, direction] * adiabatic_Hamilton - adiabatic_Hamilton * NACs[:, :, atom, direction]

                diabatic_gradients[:, :, atom, direction] = np.linalg.multi_dot([ADTmatrix.T, diab_sub_gradient, ADTmatrix])


        # TODO Hessian tranformation from adiabatic to diabatic basis oriented at the gradient implementatioin above

        NAC_derivative_matrix = np.empty(adiabatic_Hamilton.shape, dtype=object)
        hessians_matrix = np.empty(adiabatic_Hamilton.shape, dtype=object)
        diabatic_hessians = np.zeros((nstates, nstates, natoms, 3, natoms, 3))
        zero_matrix_hessian = zero_matrix = np.zeros_like(hessians[0])

        for root_1 in range(len(energies)):
            for root_2 in range(len(energies)):

                if root_1 != root_2:
                    hessians_matrix[root_1, root_2] = zero_matrix_hessian
                    NAC_derivative_matrix[root_1, root_2] = NACs_deriv[root_1 * len(energies) + root_2]

                else:
                    hessians_matrix[root_1, root_2] = hessians[root_2]
                    NAC_derivative_matrix[root_1, root_2] = zero_matrix_hessian

        for atom_1 in range(len(self.molecule.get_labels())):
            for direction_1 in range(3):
                for atom_2 in range(len(self.molecule.get_labels())):
                    for direction_2 in range(3):

                        diab_sub_hessian = (0.5 * hessians_matrix[:, :, atom_1, direction_1, atom_2, direction_2] +
                                            0.5 * hessians_matrix[:, :, atom_2, direction_2, atom_1, direction_1] +
                                            NAC_Matrix[:, :, atom_1, direction_1] * gradients_matrix[:, :, atom_2, direction_2] +
                                            NAC_Matrix[:, :, atom_2, direction_2] * gradients_matrix[:, :, atom_1, direction_1] -
                                            gradients_matrix[:, :, atom_2, direction_2] * NAC_Matrix[:, :, atom_1, direction_1] -
                                            gradients_matrix[:, :, atom_1, direction_1] * NAC_Matrix[:, :, atom_2, direction_2] -
                                            NAC_Matrix[:, :, atom_1, direction_1] * adiabatic_Hamilton * NAC_Matrix[:, :, atom_2, direction_2] -
                                            NAC_Matrix[:, :, atom_2, direction_2] * adiabatic_Hamilton * NAC_Matrix[:, :, atom_1, direction_1] +
                                            0.5 * adiabatic_Hamilton * (NAC_Matrix[:, :, atom_1, direction_1] * NAC_Matrix[:, :, atom_2, direction_2] +
                                                                        NAC_Matrix[:, :, atom_2, direction_2] * NAC_Matrix[:, :, atom_1, direction_1] -
                                                                        NAC_derivative_matrix[:, :, atom_1, direction_1, atom_2, direction_2] -
                                                                        NAC_derivative_matrix[:, :, atom_2, direction_2, atom_1, direction_1]) +
                                            0.5 * (NAC_Matrix[:, :, atom_2, direction_2] * NAC_Matrix[:, :, atom_1, direction_1] +
                                                                        NAC_Matrix[:, :, atom_1, direction_1] * NAC_Matrix[:, :, atom_2, direction_2] +
                                                                        NAC_derivative_matrix[:, :, atom_1, direction_1, atom_2, direction_2] +
                                                                        NAC_derivative_matrix[:, :, atom_2, direction_2, atom_1, direction_1]) * adiabatic_Hamilton)

                        diabatic_hessians[:, :, atom_1, direction_1, atom_2, direction_2] = np.linalg.multi_dot([ADTmatrix.T, diab_sub_hessian, ADTmatrix])

        return diabatic_energies, diabatic_gradients, diabatic_hessians


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

        # Open the file in write mode ('w')
        with open(outputfile, 'a') as file:
            # Write the section header
            
            file.write(f"\n######################################\n")
            file.write(f"############## Step {self.step} ################\n")
            file.write(f"######################################\n")
            file.write("\n########## Coordinates (Angstrom) ##########\n\n")
            
            for i, coord in enumerate(self.coordinates_xyz[self.step]):

                file.write(f'{self.molecule.get_labels()[i]}   {coord[0]:.4f}    {coord[1]:.4f}     {coord[2]:.4f}\n')

            file.write("\n########## Gradient (hatree/bohr) ##########\n\n")
            
            for i, grad in enumerate(self.all_gradients[self.step]):

                file.write(f' {grad[0]:.4f}    {grad[1]:.4f}     {grad[2]:.4f}\n')
            
            file.write("########## Kinetic Energy (kJ mol^-1) ##########\n\n")
            file.write(f"kin E = {self.kinetic_energies[self.step]:.8f}\n\n")
            file.write("########## Potential Energy (kJ mol^-1) ##########\n\n")
            file.write(f"pot E = {self.total_potentials[self.step]:.8f}\n\n")
            file.write("########## Total Energy (kJ mol^-1) ##########\n\n")
            file.write(f"tot E = {self.total_energies[self.step]:.8f}\n\n")
            file.write("########## Temperature K ##########\n\n")
            file.write(f"T = {self.temperatures[self.step]:.8f}\n\n")
            file.write("########## ENERGY GAP (kJ mol^-1) ##########\n\n")
            
            # Write the column headers
            file.write("STATE | ENERGY GAP\n\n")
            
            # Loop through the energy gaps and write each entry in the desired format
            if self.roots > 1:
                states, energy_gap = self.energy_gabs[self.step]
                states = list(states)
                    
                file.write(f"{states[0]} --> {states[1]}     {energy_gap:.4f}\n")


    def calculate_translation_coordinates_analysis(self, given_coordinates):
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
                ImpesCoordinates object
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

        distance_vector_core = (ref_structure_check - rotated_coordinates_core)
        distance_vector_core_norm = np.zeros(reference_coordinates.shape[0])

        for i in range(len(distance_vector_core_norm)):
            distance_vector_core_norm[i] += np.linalg.norm(distance_vector_core[i])

        return ref_structure_check, distance_core, distance_vector_core
        