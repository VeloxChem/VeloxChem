#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from mpi4py import MPI
from pathlib import Path
from time import time
from xml.dom import minidom
import xml.etree.ElementTree as ET
import numpy as np
import json
import h5py
import sys
from sys import stdout
import math
import random
from copy import deepcopy

from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kjpermol, bohr_in_angstrom
from .molecule import Molecule
from .outputstream import OutputStream
from .mmforcefieldgenerator import MMForceFieldGenerator
from .solvationbuilder import SolvationBuilder

from .interpolationdatapoint import InterpolationDatapoint
from .interpolationdriver import InterpolationDriver

from .imrotorbuilder import (
    RotorDefinition,
    RotorClusterDefinition,
    RotorClusterInformation,
    RotorClusterStateDefinition,
    RotorClusterAngleLibrary,
)

# from .atommapper import AtomMapper
from .errorhandler import assert_msg_critical
from contextlib import redirect_stderr
from io import StringIO

with redirect_stderr(StringIO()) as fg_err:
    import geometric

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass


class OpenMMIMDynamics:
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
        - im_driver: The VeloxChem driver object. Options are XtbDriver or an ScfDriver.
        - basis: The basis set for the QM region if an SCF driver is used.
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

        assert_msg_critical(
            'openmm' in sys.modules,
            'OpenMM is required for OpenMMDynamics.')

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        try:
            self._mpi_ctrl_comm = self.comm.Dup()
        except Exception:
            self._mpi_ctrl_comm = self.comm

        try:
            self._mpi_interp_comm = self.comm.Dup()
        except Exception:
            self._mpi_interp_comm = self._mpi_ctrl_comm

        # MPI information (based on control communicator)
        self.rank = self._mpi_ctrl_comm.Get_rank()
        self.nodes = self._mpi_ctrl_comm.Get_size()

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
        self.unique_residues = []
        self.unique_molecules = []

        # filtering
        self.rmsd_threshold = 1.5
        self.energy_threshold = 1.5 # kj/mol

        # QM Region parameters
        self.im_drivers = {}
        self.use_symmetry = True
        self.symmetry_information = None
        self.excitation_pulse = (-1, 0)
        self.current_state = None
        self.current_molecule = None

        self.roots_z_matrix = None
        self.roots_to_follow = None
        self.qm_symmetry_datapoint_dict = None
        self.qm_data_point_dict = None
        self.qm_rotor_cluster_banks = None
        self.sorted_state_spec_im_labels = None
        self.root_spec_molecules = None
        self.interpolation_settings = None
        self.dihedrals_dict = None
        self.basis_set_label = 'def2-svp'
        self.ab_init_drivers = None

        self.step = None
        self.basis = None
        self.grad_driver = None
        self.qm_atoms = None
        self.mm_subregion = None
        self.linking_atoms = None
        self.qm_force_index = None
        self.driver_flag = None

        # QM stabilizer parameters
        self.scaling_factor = 0.01

        # Default value for the C-H linker distance
        self.linking_atom_distance = 1.0705 

    # Loading methods
    # # TODO: Integrate the guess with the read_pdb_file in Molecule.
    # def load_system_PDB(self, filename):
    #     """
    #     Loads a system from a PDB file, extracting unique residues and their details.

    #     :param filename: 
    #         Path to the PDB file. Example: 'system.pdb'.
    #     """
        
    #     pdb_path = Path(filename)
    #     if not pdb_path.is_file():
    #         raise FileNotFoundError(f"{filename} does not exist.")
        
    #     pdbstr = pdb_path.read_text()
    #     residues = {}

    #     # Formatting issues flags
    #     label_guess_warning = False
    #     conect_warning = True

    #     self.unique_residues = []  # Initialize unique_residues earlier

    #     for line in pdbstr.strip().splitlines():
    #         # Skip the header lines
    #         if line.startswith(('REMARK', 'HEADER', 'TITLE', 'COMPND', 'SOURCE')):
    #             continue
    #         # Skip the CRYST1 record it will be extracted from the System object
    #         if line.startswith('CRYST1'):
    #             continue
    #         # Add the ATOM and HETATM records to the residues dictionary
    #         if line.startswith(('ATOM', 'HETATM')):
    #             atom_label = line[76:78].strip()
    #             if not atom_label:
    #                 label_guess_warning = True
    #                 atom_name = line[12:16].strip()
    #                 atom_label = atom_name[0]

    #             residue = line[17:20].strip()
    #             residue_number = int(line[22:26])
    #             coordinates = [float(line[i:i+8]) for i in range(30, 54, 8)]

    #             residue_identifier = (residue, residue_number)

    #             if residue_identifier not in residues:
    #                 residues[residue_identifier] = {
    #                     'labels': [],
    #                     'coordinates': []
    #                 }
                
    #             residues[residue_identifier]['labels'].append(atom_label)
    #             residues[residue_identifier]['coordinates'].append(coordinates)
            
    #         # If a line starts with TER, skip it
    #         if line.startswith('TER'):
    #             continue
    #         # If any line starts with CONECT, set the conect_warning flag to False
    #         if line.startswith('CONECT'):
    #             conect_warning = False

    #         # If a line starts with END, break the loop
    #         if line.startswith('END'):
    #             break
            

    #     # Overwrite the PDB file with the guessed atom labels in columns 77-78
    #     # if the labels are missing
    #     if label_guess_warning:
    #         with open(filename, 'w') as f:
    #             for line in pdbstr.strip().splitlines():
    #                 if line.startswith(('ATOM', 'HETATM')):
    #                     atom_name = line[12:16].strip()
    #                     atom_label = atom_name[0]
    #                     new_line = line[:76] + f"{atom_label:>2}" + line[78:] + '\n'
    #                     f.write(new_line)
    #                 else:
    #                     f.write(line + '\n')
    #         msg = 'Atom labels were guessed based on atom names (first character).'
    #         msg += f'Please verify the atom labels in the {filename} PDB file.'
    #         self.ostream.print_warning(msg)
    #         self.ostream.flush()

    #     if conect_warning:
    #         msg = 'CONECT records not found in the PDB file.'
    #         msg += 'The connectivity matrix will be used to determine the bonds.'
    #         self.ostream.print_warning(msg)
    #         self.ostream.flush()
    #         # Create a molecule from the PDB file
    #         molecule = Molecule.read_pdb_file(filename)
    #         connectivity_matrix = molecule.get_connectivity_matrix()
    #         # Determine all the bonds in the molecule
    #         with open(filename, 'a') as f:
    #             for i in range(connectivity_matrix.shape[0]):
    #                 for j in range(i + 1, connectivity_matrix.shape[1]):
    #                     if connectivity_matrix[i, j] == 1:
    #                         # Convert indices to 1-based index for PDB format and ensure proper column alignment
    #                         i_index = i + 1
    #                         j_index = j + 1
    #                         # Align to the right 
    #                         con_string = "{:6s}{:>5d}{:>5d}".format('CONECT', i_index, j_index)
    #                         f.write(con_string + '\n')

    #     # Create VeloxChem Molecule objects for each unique residue

    #     molecules = []
    #     unq_residues = []

    #     for (residue, number), data in residues.items():
    #         coordinates_array = np.array(data['coordinates'])
    #         mol = Molecule(data['labels'], coordinates_array, "angstrom")
    #         molecules.append(mol)
    #         unq_residues.append((residue, number))


    #     # Initialize a set to track the first occurrence of each residue name
    #     seen_residue_names = set()

    #     # Lists to store the filtered residues and molecules
    #     self.unique_residues = []
    #     self.unique_molecules = []

    #     for index, (residue_name, number) in enumerate(unq_residues):
    #         if residue_name not in seen_residue_names:
    #             seen_residue_names.add(residue_name)
    #             self.unique_residues.append((residue_name, number))
    #             self.unique_molecules.append(molecules[index])

    #     # Print results
    #     info_msg = f"Unique Residues: {self.unique_residues}, saved as molecules."
    #     self.ostream.print_info(info_msg)
    #     self.ostream.flush()
    
    # def load_system_from_files(self, system_xml, system_pdb):
    #     """
    #     Loads a preexisting system from an XML file and a PDB file.

    #     :param system_xml: 
    #         XML file containing the system parameters.
    #     :param system_pdb: 
    #         PDB file containing the system coordinates.
    #     """
    #     self.pdb = app.PDBFile(system_pdb)

    #     with open(system_xml, 'r') as file:
    #         self.system = mm.XmlSerializer.deserialize(file.read())

    #     self.phase = 'gas'
    #     self.qm_atoms = []
    #     self.qm_force_index = -1  # Initialize with an invalid index

    #     # Iterate over all forces in the system to find the CustomExternalForce and detect QM atoms
    #     for i, force in enumerate(self.system.getForces()):
    #         if isinstance(force, mm.CustomExternalForce):
    #             self.qm_force_index = i  # Set the index of the QM force
    #             # Assuming that the particle indices are directly stored in the force
    #             for j in range(force.getNumParticles()):
    #                 # Get the actual particle index (usually force index and particle index are the same)
    #                 index, parameters = force.getParticleParameters(j)
    #                 self.qm_atoms.append(index)
    #             msg = f"QM Atoms detected: {self.qm_atoms}"
    #             self.ostream.print_info(msg)
    #             self.ostream.flush()
    #             break

    #     if self.qm_force_index == -1:
    #         msg = "No CustomExternalForce found, no QM atoms detected."
    #         self.ostream.print_info(msg)
    #         self.ostream.flush()
    #     else:
    #         msg = f"CustomExternalForce found at index {self.qm_force_index}."
    #         self.ostream.print_info(msg)
    #         self.ostream.flush()
            
    # Method to generate OpenMM system from VeloxChem objects
    def create_system_from_molecule(self, 
                             molecule,
                             interpol_settings, 
                             ff_gen, 
                             solvent='gas', 
                             qm_atoms=None, 
                             filename='residue', 
                             residue_name='MOL'):
        """
        Creates an OpenMM system from a VeloxChem molecule and a forcefield generator.
        
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
        self.interpolation_settings = interpol_settings

        self.set_up_the_system(molecule, interpol_settings)

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

        if solvent == 'gas':
            phase = 'gas'
            self.pdb = app.PDBFile(f'{filename}.pdb')     
            # Common forcefield loading, modified according to phase specifics
            forcefield_files = [f'{filename}.xml']

        if solvent != 'gas':
            # Solvate the molecule using the SolvationBuilder
            phase = 'periodic'
            sol_builder = SolvationBuilder()
            sol_builder.steps = 10000
            sol_builder.solvate(solute=molecule, 
                                solvent=solvent,
                                padding=self.padding,
                                equilibrate=True,
                                )
            
            sol_builder.write_openmm_files(solute_ff=ff_gen)
            
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

            # Load the PDB from the SolvationBuilder
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

    # # Methods to build a custom system from a PDB file and custom XML files
    # def create_md_system_from_files(self, 
    #                         pdb_file, 
    #                         xml_file,
    #                         filename='custom'):
    #     """
    #     Creates a system from a PDB file containing multiple residues and custom XML files.
        
    #     :param system_pdb:
    #         PDB file containing the system or a list of PDB files.
    #     :param xml_file:
    #         XML file containing the forcefield parameters or a list of XML files.
    #     :param filename:
    #         Base filename for output and intermediate files. Default is 'custom_system'.
    #     """
        
    #     # Load the PDB file and format it correctly
    #     # This method already prints the unique residues info!
    #     self.load_system_PDB(pdb_file)

    #     # Load the generated QM region topology and system PDB
    #     self.pdb = app.PDBFile(pdb_file)
    #     # Check if the system is periodic from the PDB file
    #     if self.pdb.topology.getUnitCellDimensions() is not None:
    #         periodic = True
    #         msg = 'PBC detected from the PDB file.'
    #         self.ostream.print_info(msg)
    #         self.ostream.flush()
    #     else:
    #         periodic = False

    #     # Combine XML files into a forcefield
    #     xml_files = [xml_file] if isinstance(xml_file, str) else xml_file
    #     msg = f'Added XML Files: {xml_files}'
    #     self.ostream.print_info(msg)
    #     self.ostream.flush()

    #     # Combine XML files into a forcefield
    #     forcefield = app.ForceField(*xml_files, self.parent_ff, self.water_ff)

    #     # Create the system with appropriate nonbonded settings
    #     if periodic:
    #         nonbondedMethod = app.PME  # Particle Mesh Ewald method for periodic systems
    #         nonbondedCutoff = self.cutoff * unit.nanometer  # Assuming self.cutoff is defined elsewhere appropriately
    #         constraints = app.HBonds
    #     else:
    #         nonbondedMethod = app.NoCutoff  # No cutoff for non-periodic systems
    #         nonbondedCutoff = None
    #         constraints = app.HBonds

    #     # Check if nonbondedCutoff is defined before using it in createSystem
    #     system_arguments = {
    #         'nonbondedMethod': nonbondedMethod,
    #         'constraints': constraints
    #     }   
        
    #     if nonbondedCutoff is not None:
    #         system_arguments['nonbondedCutoff'] = nonbondedCutoff

    #     self.system = forcefield.createSystem(self.pdb.topology, 
    #                                           **system_arguments)

    #     # Save system to XML and PDB for inspection and reuse
    #     with open(f'{filename}_system.xml', 'w') as f:
    #         f.write(mm.XmlSerializer.serialize(self.system))
    #         msg = f'System parameters written to {filename}_system.xml'
    #         self.ostream.print_info(msg)
    #         self.ostream.flush()

    #     with open(f'{filename}_system.pdb', 'w') as pdb_file:
    #         app.PDBFile.writeFile(self.pdb.topology, 
    #                               self.pdb.positions, 
    #                               pdb_file)
    #         msg = f'System coordinates written to {filename}_system.pdb'
    #         self.ostream.print_info(msg)
    #         self.ostream.flush()

    #     self.phase = 'gas'

    # def create_qmmm_system_from_files(self, 
    #                         pdb_file, 
    #                         xml_file, 
    #                         qm_residue, 
    #                         ff_gen_qm=None, 
    #                         qm_atoms='all', 
    #                         filename='custom'):
    #     """
    #     Builds a QM/MM system from a PDB file containing multiple residues and custom XML files.
        
    #     :param pdb_file: 
    #         PDB file containing the system.
    #     :param xml_file: 
    #         XML file containing the forcefield parameters. Can be a list of XML files.
    #     :param qm_residue: 
    #         Tuple containing the name and number of the QM residue. (e.g. ('MOL', 1))
    #     :param ff_gen_qm: 
    #         Optional custom forcefield generator for the QM region. If None, a new one will be created.
    #     :param qm_atoms:
    #         Options: 'all' or a list of atom indices for QM region.
    #     :param filename: 
    #         Base filename for output and intermediate files.
    #     """

    #     # Load the PDB file and format it correctly
    #     self.load_system_PDB(pdb_file)

    #     # Extracting the residue name and number for clarity and correct usage
    #     qm_residue_name, qm_residue_number = qm_residue
    #     qm_residue_index = self.unique_residues.index((qm_residue_name, qm_residue_number))
    #     msg = f'QM Residue: {qm_residue_name}'
    #     self.ostream.print_info(msg)
    #     self.ostream.flush()

    #     qm_molecule = self.unique_molecules[qm_residue_index]

    #     # Generate or use an existing forcefield generator for the QM region
    #     if ff_gen_qm is None:
    #         ff_gen_qm = MMForceFieldGenerator()
    #         ff_gen_qm.create_topology(qm_molecule)

    #     # Determine qm_atoms based on the QM residue
    #     if qm_atoms == 'all':
    #         msg = 'Full molecule as QM region'
    #         self.ostream.print_info(msg)
    #         qm_atoms = list(range(qm_molecule.number_of_atoms()))
    #         self._create_QM_residue(ff_gen_qm,
    #                     qm_atoms, 
    #                     filename=filename, 
    #                     residue_name=qm_residue_name)
    #     elif isinstance(qm_atoms, list):
    #         if qm_atoms == list(range(qm_molecule.number_of_atoms())): 
    #             msg = 'Full molecule as QM region'
    #             self.ostream.print_info(msg)
    #             self._create_QM_residue(ff_gen_qm,
    #                     qm_atoms, 
    #                     filename=filename, 
    #                     residue_name=qm_residue_name)
    #         msg = 'QM/MM partition inside the molecule'
    #         self.ostream.print_info(msg)
    #         qm_atoms = qm_atoms
    #         self._create_QM_subregion(ff_gen_qm,
    #                                   qm_atoms, 
    #                                   qm_molecule, 
    #                                   filename,
    #                                   residue_name=qm_residue_name)
    #     else:
    #         raise ValueError('Invalid value for qm_atoms. Please use "all" or a list of atom indices.')

    #     # Load the generated QM region topology and system PDB
    #     self.pdb = app.PDBFile(pdb_file)

    #     # Check if the system is periodic from the PDB file
    #     if self.pdb.topology.getUnitCellDimensions() is not None:
    #         periodic = True
    #         msg = 'PBC detected from the PDB file.'
    #         self.ostream.print_info(msg)
    #         self.ostream.flush()

    #     # Get self.positions from the topology
    #     self.positions = self.pdb.positions

    #     # Combine XML files into a forcefield
    #     xml_files = [f'{filename}.xml'] + ([xml_file] if isinstance(xml_file, str) else xml_file)
    #     msg = f'Added XML Files: {xml_files}'
    #     self.ostream.print_info(msg)
    #     self.ostream.flush()

    #     # Combine XML files into a forcefield
    #     forcefield = app.ForceField(*xml_files)

    #     # Create the system with appropriate nonbonded settings
    #     if periodic:
    #         nonbondedMethod = app.PME  # Particle Mesh Ewald method for periodic systems
    #         nonbondedCutoff = self.cutoff * unit.nanometer  # Assuming self.cutoff is defined elsewhere appropriately
    #         constraints = app.HBonds
    #     else:
    #         nonbondedMethod = app.NoCutoff  # No cutoff for non-periodic systems
    #         nonbondedCutoff = None
    #         constraints = app.HBonds

    #     # Check if nonbondedCutoff is defined before using it in createSystem
    #     system_arguments = {
    #         'nonbondedMethod': nonbondedMethod,
    #         'constraints': constraints
    #     }   

    #     if nonbondedCutoff is not None:
    #         system_arguments['nonbondedCutoff'] = nonbondedCutoff

    #     self.system = forcefield.createSystem(self.pdb.topology, **system_arguments)

    #     # Set QM atoms based on the QM residue
    #     #self.qm_atoms = [atom.index for atom in self.pdb.topology.atoms() if (atom.residue.name == qm_residue_name and atom.residue.id == str(qm_residue_number))]
    #     self.qm_atoms = qm_atoms

    #     # Setting QM/MM system
    #     self.set_qm_mm_system('periodic' if periodic else 'gas', ff_gen_qm)

    #     # Adding the stabilizer force to the QM region
    #     # self.qm_stabilizer(ff_gen_qm)

    #     # Save system to XML and PDB for inspection and reuse
    #     with open(f'{filename}_system.xml', 'w') as f:
    #         f.write(mm.XmlSerializer.serialize(self.system))
    #         msg = f'System parameters written to {filename}_system.xml'
    #         self.ostream.print_info(msg)
    #         self.ostream.flush()

    #     with open(f'{filename}_system.pdb', 'w') as pdb_file:
    #         app.PDBFile.writeFile(self.pdb.topology, self.pdb.positions, pdb_file)
    #         msg = f'System coordinates written to {filename}_system.pdb'
    #         self.ostream.print_info(msg)
    #         self.ostream.flush()

    #     # Correct phase setting
    #     self.phase = 'gas'
    
    def _reload_interpolation_root_from_hdf5(self, root, inv_sqrt_masses):
        driver_object = self.im_drivers[root]

        im_labels, _ = driver_object.read_labels()

        self.qm_data_point_dict[root] = []
        self.qm_symmetry_datapoint_dict[root] = {}
        self.sorted_state_spec_im_labels[root] = []

        old_label = None

        for label in im_labels:
            qm_data_point = InterpolationDatapoint(self.roots_z_matrix[root])
            qm_data_point.update_settings(self.interpolation_settings[root])
            qm_data_point.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)
            qm_data_point.inv_sqrt_masses = inv_sqrt_masses
            
            if qm_data_point.bank_role == "core" or "cluster" not in qm_data_point.point_label and "symmetry" not in qm_data_point.point_label:


                self.qm_data_point_dict[root].append(qm_data_point)
                self.sorted_state_spec_im_labels[root].append(label)

                old_label = qm_data_point.point_label
                self.qm_symmetry_datapoint_dict[root][old_label] = [qm_data_point]
            elif qm_data_point.bank_role == "symmetry":

                self.qm_symmetry_datapoint_dict[root][old_label].append(qm_data_point)

        driver_object.qm_symmetry_data_points = self.qm_symmetry_datapoint_dict[root]
        driver_object.qm_data_points = self.qm_data_point_dict[root]
        driver_object.labels = self.sorted_state_spec_im_labels[root]
        
        print(len(self.qm_data_point_dict[root]), driver_object.impes_coordinate.eq_bond_lengths, self.qm_data_point_dict[root][0].eq_bond_lengths)
        if len(self.qm_data_point_dict[root]) > 0:
            driver_object.impes_coordinate.eq_bond_lengths = self.qm_data_point_dict[root][0].eq_bond_lengths

        driver_object.mark_runtime_data_cache_dirty()
        driver_object.prepare_runtime_data_cache(force=True)

    def _load_rotor_cluster_bank_for_root(self, root):
        out = {}
        imff_file = self.interpolation_settings[root]["imforcefield_file"]
        families = self._list_rotor_cluster_families_from_registry(root)

        for family in families:
            cluster_info, angle_library, point_index = self._read_cluster_registry_for_family(
                imff_file, root, family
            )
            fam = {
                "cluster_info": cluster_info,
                "cluster_angle_library": angle_library,
                "point_index": point_index,
                "core": None,
                "clusters": {
                    cid: {
                        "cluster_type": cluster.cluster_type,
                        "rotor_ids": tuple(cluster.rotor_ids),
                        "expected_states": {state.state_id: None for state in angle_library.state_banks[cid]},
                    }
                    for cid, cluster in cluster_info.clusters.items()
                },
            }

            core_dp = InterpolationDatapoint(self.roots_z_matrix[root])
            core_dp.update_settings(self.interpolation_settings[root])
            core_dp.read_hdf5(imff_file, point_index["core_label"])
            fam["core"] = core_dp

            for cid, state_map in point_index["cluster_state_labels"].items():
                for sid, label in state_map.items():
                    dp = InterpolationDatapoint(self.roots_z_matrix[root])
                    dp.update_settings(self.interpolation_settings[root])
                    dp.read_hdf5(imff_file, label)
                    fam["clusters"][cid]["expected_states"][sid] = dp

            for cid, cbank in fam["clusters"].items():
                expected = cbank["expected_states"]
                if 0 in expected and expected[0] is None:
                    expected[0] = fam["core"]

            out[family] = fam
        return out
    
    def _list_rotor_cluster_families_from_registry(self, root):
        imff_file = self.interpolation_settings[root]["imforcefield_file"]
        prefix = f"rotor_cluster_registry/root_{root}"
        with h5py.File(imff_file, "r") as h5f:
            if prefix not in h5f:
                return []
            return sorted(name.replace("family_", "", 1) for name in h5f[prefix].keys())


    def _read_cluster_registry_for_family(self, imff_file, root, family_label):
        def _read_scalar_string(ds):
            val = ds[()]
            return val.decode("utf-8") if isinstance(val, bytes) else str(val)

        with h5py.File(imff_file, "r") as h5f:
            prefix = f"rotor_cluster_registry/root_{root}/family_{family_label}"
            info_json = _read_scalar_string(h5f[prefix + "/cluster_info_json"])
            library_json = _read_scalar_string(h5f[prefix + "/cluster_angle_library_json"])
            index_json = _read_scalar_string(h5f[prefix + "/point_index_json"])

        info_payload = json.loads(info_json)
        library_payload = json.loads(library_json)
        index_payload = json.loads(index_json)
        
        rotor_to_cluster_map = {
            int(rid): int(cid) for rid, cid in info_payload["rotor_to_cluster"].items()
        }
        cluster_info = RotorClusterInformation(
            dihedral_start=int(info_payload["dihedral_start"]),
            dihedral_end=int(info_payload["dihedral_end"]),
            rotor_to_cluster=rotor_to_cluster_map,
        )
        for rotor_id, rotor_data in info_payload["rotors"].items():
            cluster_info.rotors[int(rotor_id)] = RotorDefinition(
                rotor_id=int(rotor_id),
                center=tuple(rotor_data["center"]),
                torsion_rows=tuple(rotor_data["torsion_rows"]),
                torsion_coords=tuple(tuple(x) for x in rotor_data["torsion_coords"]),
                symmetry_order=int(rotor_data["symmetry_order"]),
                atom_group=tuple(rotor_data["atom_group"]),
            )
        for cid, cdata in info_payload["clusters"].items():
            cluster_info.clusters[int(cid)] = RotorClusterDefinition(
                cluster_id=int(cid),
                rotor_ids=tuple(cdata["rotor_ids"]),
                cluster_type=cdata["cluster_type"],
                torsion_rows=tuple(cdata["torsion_rows"]),
            )

        angle_library = RotorClusterAngleLibrary()
        for cid, state_list in library_payload.items():
            bank = []
            for state in state_list:
                aa = {}
                for key, val in state["angle_assignment"].items():
                    aa[tuple(int(x) for x in key.split(","))] = float(val)
                bank.append(
                    RotorClusterStateDefinition(
                        cluster_id=int(cid),
                        state_id=int(state["state_id"]),
                        cluster_type=state["cluster_type"],
                        rotor_ids=tuple(state["rotor_ids"]),
                        angle_assignment=aa,
                        dihedrals_to_rotate=(
                            None
                            if state["dihedrals_to_rotate"] is None
                            else tuple(tuple(int(x) for x in row) for row in state["dihedrals_to_rotate"])
                        ),
                        phase_signature=(
                            None if state["phase_signature"] is None
                            else np.asarray(state["phase_signature"], dtype=np.float64)
                        ),
                        is_anchor=bool(state.get("is_anchor", False)),
                        label_suffix=str(state.get("label_suffix", "")),
                    )
                )
            angle_library.state_banks[int(cid)] = tuple(bank)

        point_index = {
            "family_label": index_payload["family_label"],
            "core_label": index_payload["core_label"],
            "cluster_state_labels": {
                int(cid): {int(sid): lbl for sid, lbl in smap.items()}
                for cid, smap in index_payload["cluster_state_labels"].items()
            },
        }
        return cluster_info, angle_library, point_index

    def run_immm(self, 
                 roots_to_follow=[0],
                 restart_file = None, 
                 ensemble='NVE', 
                 temperature=298.15, 
                 pressure=1.0, 
                 friction=1.0,
                 timestep=0.5,
                 nsteps=1000, 
                 snapshots=100, 
                 traj_file='trajectory.pdb',
                 state_file='output.xml',
                 output_file='output'):
        """
        Runs a QM/MM simulation using OpenMM, storing the trajectory and simulation data.

        :param im_driver:
            QM driver object from VeloxChem.
        :param grad_driver:
            Gradient driver object from VeloxChem.
        :param basis:
            Molecular basis string.
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
        if any(roots_to_follow) > 0:
            raise RuntimeError('Currently Interpolation is restricted to ground state dynamics!')
        self.ensemble = ensemble
        self.temperature = temperature * unit.kelvin
        self.friction = friction / unit.picosecond
        self.timestep = timestep * unit.femtoseconds
        self.nsteps = nsteps
      
        interpolation_settings_dict = {}
        if len(self.interpolation_settings) != len(roots_to_follow):
            for i, root in enumerate(roots_to_follow):
                interpolation_settings_dict[root] = self.interpolation_settings[i]
        else:
            for i, root in enumerate(roots_to_follow):
                interpolation_settings_dict[root] = self.interpolation_settings[i]
        self.roots_to_follow = roots_to_follow
        self.qm_symmetry_datapoint_dict = {root: {} for root in self.roots_to_follow}
        self.qm_data_point_dict = {root: [] for root in self.roots_to_follow}
        self.sorted_state_spec_im_labels = {root: [] for root in self.roots_to_follow}
        self.root_spec_molecules = {root: [] for root in self.roots_to_follow}
        self.qm_rotor_cluster_banks = {root: {} for root in self.roots_to_follow}
        self.im_drivers = {root: None for root in self.roots_to_follow}
        self.interpolation_settings = interpolation_settings_dict
        masses = self.molecule.get_masses().copy()
    
        masses_cart = np.repeat(masses, 3)
        inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
        self.inv_sqrt_masses = inv_sqrt_masses
        for root in self.roots_to_follow:
            # Dynamically create an attribute name
            attribute_name = f'impes_driver_{root}'
            # Initialize the object
            driver_object = InterpolationDriver(self.roots_z_matrix[root])
            driver_object.update_settings(self.interpolation_settings[root])
            # print('Interpolation driver settings updated for root', root, self.eq_bond_force_constants)
            driver_object.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
            if root == 0:
                driver_object.symmetry_information = self.symmetry_information['gs']
            else:
                driver_object.symmetry_information = self.symmetry_information['es']
           
            self.im_drivers[root] = driver_object
            self._reload_interpolation_root_from_hdf5(root, inv_sqrt_masses)
            self.qm_rotor_cluster_banks[root] = self._load_rotor_cluster_bank_for_root(root)
            driver_object.qm_rotor_cluster_banks = self.qm_rotor_cluster_banks[root]
            if driver_object.qm_rotor_cluster_banks:
                first_family = next(iter(driver_object.qm_rotor_cluster_banks.values()))
                driver_object.rotor_cluster_information = first_family.get("cluster_info")
            else:
                driver_object.rotor_cluster_information = None
            print('In set up', driver_object.rotor_cluster_information)
            # Set the object as an attribute of the instance
            setattr(self, attribute_name, driver_object)
            # Append the object to the list
            self.im_drivers[root] = driver_object
            self.current_state = self.roots_to_follow[0]


        runtime_cache = self._build_run_qmmm_runtime_cache()
        self.qm_potentials = []
        self.qm_mm_interaction_energies = []
        self.mm_potentials = []
        self.total_potentials = []
        self.kinetic_energies = []
        self.temperatures = []
        self.total_energies = []
        self.dynamic_molecules = []
        
        save_freq = nsteps // snapshots if snapshots else nsteps

        # Create or update the integrator
        if self.integrator is None:
            new_integrator = self._create_integrator()
        else:
            new_integrator = self.integrator

        self.topology = self.pdb.topology

        self.positions = self.pdb.positions
        
        self.simulation = app.Simulation(self.topology, self.system, new_integrator)

        # Load the state if a restart file is provided
        if restart_file is not None:
            self.simulation.loadState(restart_file)
        
        else:
            self.simulation.context.setPositions(self.positions)

            # Set initial velocities if the ensemble is NVT or NPT
            if self.ensemble in ['NVT', 'NPT']:
                self.simulation.context.setVelocitiesToTemperature(self.temperature)
            # else:
            #     self.simulation.context.setVelocitiesToTemperature(10 * unit.kelvin)

        # There is no minimization step for QM/MM simulations
        # It causes instabilities in the MM region!
        
        # Set up reporting
        self.simulation.reporters.clear()
        self.simulation.reporters.append(app.PDBReporter(traj_file, save_freq))
        
        # Print header
        print('QM/MM Simulation Parameters')
        print('=' * 60)
        print('QM Driver:', self.driver_flag)
        print('Ensemble:', ensemble)
        print('Integration method:', new_integrator.__class__.__name__)
        if ensemble in ['NVT', 'NPT']:
            print('Temperature:', temperature, 'K')
            if ensemble == 'NPT':
                print('Pressure:', pressure, 'atm')
        print('Friction:', friction, '1/ps')
        print('Timestep:', timestep, 'fs')
        print('Total simulation time in ns:', nsteps * timestep / 1e6)
        print('=' * 60)

        start_time = time()
        self.step = 0
        for step in range(nsteps):

            # if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root():
            # qm_positions_nm = np.array([
            #     p.value_in_unit(unit.nanometer)
            #     for p in self.simulation.context.getState(getPositions=True).getPositions()
            # ])[self.qm_atoms]

                # sync_roots = sorted(self._mpi_pending_sync_roots)
                # self._mpi_pending_sync_roots.clear()

                # if self.mpi_reload_from_hdf5 and len(sync_roots) > 0:
                #     for root_sync in sync_roots:
                #         self._reload_interpolation_root_from_hdf5(root_sync, self.inv_sqrt_masses)

                # self._mpi_bcast_control({
                #     'cmd': self._MPI_CMD_STEP,
                #     'step_idx': int(step),
                #     'qm_positions_nm': np.asarray(qm_positions_nm, dtype=np.float64),
                #     'sync_roots': sync_roots,
                # })

            self.update_forces(self.simulation.context)
            observables = self._collect_step_observables_optimized(
                    self.simulation.context,
                    runtime_cache)
            
            self._append_step_observables(observables)

            self._print_run_qmmm_step(step, timestep, observables, save_freq)
            # Potential energies
            # QM region
            # qm = self.get_qm_potential_energy()
            # self.qm_potentials.append(qm)

            # # QM/MM interactions
            # qm_mm = self.simulation.context.getState(getEnergy=True, groups={1,2}).getPotentialEnergy()
            # self.qm_mm_interaction_energies.append(qm_mm.value_in_unit(unit.kilojoules_per_mole))

            # # MM region 
            # mm = self.simulation.context.getState(getEnergy=True, groups={3,4,5,6,7}).getPotentialEnergy()
            # self.mm_potentials.append(mm.value_in_unit(unit.kilojoules_per_mole))

            # # Total potential energy
            # pot = qm * unit.kilojoules_per_mole + qm_mm + mm
            # self.total_potentials.append(pot.value_in_unit(unit.kilojoules_per_mole))

            # # Kinetic energyex
            
            # kinetic = self.simulation.context.getState(getEnergy=True).getKineticEnergy()
            # self.kinetic_energies.append(kinetic.value_in_unit(unit.kilojoules_per_mole))

            # # Temperature
            # dof = self.system.getNumParticles() * 3 - self.system.getNumConstraints()
            # if hasattr(self.integrator, 'computeSystemTemperature'):
            #     temp = self.integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
            # else:
            #     temp = (2 * kinetic / (dof * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)
            # self.temperatures.append(temp)

            # # Total energy
            # total = pot + kinetic
            # self.total_energies.append(total.value_in_unit(unit.kilojoules_per_mole))


            # # Information output
            # if step % save_freq == 0:
                
            #     print(f"Step: {step} / {nsteps} Time: {round((step * timestep) / 1000, 2)} ps")
            #     print('Potential Energy QM region:', qm, 'kJ/mol')
            #     print('Potential Energy MM region:', mm)
            #     print('QM/MM Interaction Energy:', qm_mm)
            #     print('Total Potential Energy:', pot)
            #     print('Kinetic Energy:', kinetic)
            #     print('Temperature:', temp, 'K')
            #     print('Total Energy:', total)
            #     print('Current State (PES):', self.current_state)  
            #     print('-' * 60)   

            #     self.dynamic_molecules.append(self.current_molecule)    

            self.simulation.step(1)
            self.step += 1


        self.structures_to_xyz_file(self.dynamic_molecules, 'snapshots_traj.xyz')
        end_time = time()
        elapsed_time = end_time - start_time
        elapsed_time_days = elapsed_time / (24 * 3600)
        performance = (nsteps * timestep / 1e6) / elapsed_time_days

        print('QM/MM simulation completed!')
        print(f'Number of steps: {nsteps}')
        print(f'Trajectory saved as {traj_file}')
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
        print(f'Trajectory saved as {traj_file}')

        # Save the simulation data
        self.simulation.saveState(state_file)
        print(f'Simulation state saved as {state_file}')

        # Write the output to a file
        self._save_output(output_file)
        self.ostream.print_info(f'Simulation report saved as {output_file}.out')
        self.ostream.flush()

    def _collect_step_observables_optimized(self, context, runtime_cache):
        """
        Optimized observable collection avoiding repeated Quantity arithmetic.
        """

        energy_unit = runtime_cache['energy_unit']

        qm = float(self.current_energy)
        qm_mm = context.getState(
            getEnergy=True,
            groups=runtime_cache['qm_mm_groups']).getPotentialEnergy().value_in_unit(energy_unit)
        mm = context.getState(
            getEnergy=True,
            groups=runtime_cache['mm_groups']).getPotentialEnergy().value_in_unit(energy_unit)
        kinetic = context.getState(getEnergy=True).getKineticEnergy().value_in_unit(energy_unit)

        potential = qm + qm_mm + mm

        if runtime_cache['has_temp_method']:
            temperature = self.integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
        else:
            temperature = (2.0 * kinetic) / (runtime_cache['dof'] * runtime_cache['gas_constant'])

        total = potential + kinetic

        return {
            'qm': qm,
            'qm_mm': qm_mm,
            'mm': mm,
            'potential': potential,
            'kinetic': kinetic,
            'temperature': float(temperature),
            'total': total,
        }

    def _append_step_observables(self, observables):
        """
        Appends per-step observables to trajectory history.
        """

        self.qm_potentials.append(observables['qm'])
        self.qm_mm_interaction_energies.append(observables['qm_mm'])
        self.mm_potentials.append(observables['mm'])
        self.total_potentials.append(observables['potential'])
        self.kinetic_energies.append(observables['kinetic'])
        self.temperatures.append(observables['temperature'])
        self.total_energies.append(observables['total'])
    
    def _print_run_qmmm_step(self, step, timestep, observables, save_freq):
        """
        Prints per-step simulation diagnostics.
        """
        
        if step % save_freq != 0:
            return
 
        self.dynamic_molecules.append(self.current_molecule)
        print(f"Step: {step} / {self.nsteps} Time: {round((step * timestep) / 1000, 2)} ps")
        print('Potential Energy QM region:', observables['qm'], 'kJ/mol')
        print('Potential Energy MM region:', observables['mm'])
        print('QM/MM Interaction Energy:', observables['qm_mm'])
        print('Total Potential Energy:', observables['potential'])
        print('Kinetic Energy:', observables['kinetic'])
        print('Temperature:', observables['temperature'], 'K')
        print('Total Energy:', observables['total'], '±', np.std(self.total_energies), 'kJ/mol') 
        print('Current State (PES):', self.current_state)
        print('-' * 60)
    
    def _build_run_qmmm_runtime_cache(self):
        """
        Builds cached constants used during the QM/MM integration loop.
        """

        energy_unit = unit.kilojoules_per_mole
        gas_constant = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(energy_unit / unit.kelvin)

        return {
            'energy_unit': energy_unit,
            'qm_mm_groups': {1, 2},
            'mm_groups': {3, 4, 5, 6, 7},
            'dof': self.system.getNumParticles() * 3 - self.system.getNumConstraints(),
            'has_temp_method': hasattr(self.integrator, 'computeSystemTemperature'),
            'gas_constant': gas_constant,
        }

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
            Phase of the system ('gas', 'periodic').
        :param ff_gen:
            MMForceFieldGenerator object from VeloxChem.
        """

        from openmm import NonbondedForce

        # Set the QM/MM Interaction Groups
        total_atoms = self.system.getNumParticles()
        
        # The MM subregion is counted as regular MM atoms
        qm_group = list(self.qm_atoms)
        msg = f'QM region: {qm_group[0]} ... {qm_group[-1]}'
        self.ostream.print_info(msg)
        mm_group = list(set(range(total_atoms)) - set(qm_group))
       
        if mm_group != []:
            msg = f'MM region: {mm_group[0]} ... {mm_group[-1]}'
            self.ostream.print_info(msg)
        if not mm_group:
            self.ostream.print_info('No external MM atoms found in the system')
        self.ostream.flush()

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
            if phase == 'periodic':
                vdw.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
                vdw.setCutoffDistance(self.cutoff)
            vdw.addPerParticleParameter("sigma")
            vdw.addPerParticleParameter("epsilon")


            if phase == 'periodic':
                # OpenMM uses Reaction Field method for CustomNB with PBC.
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
    
    def qm_stabilizer(self, ff_gen_qm):
        
        """
        Implements a MM potential to stabilize the QM region.
        The forces are 1% of the regular MM forces.

        :param qm_atoms: 
            List of atom indices to be included in the QM region.
        :param ff_gen_qm: 
            MMForceFieldGenerator object from VeloxChem.
        """


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

        new_molecule = None
        positions_ang = (new_positions) * 10 
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
            self.current_molecule = new_molecule
 

        else:
            # Atom labels for the QM region
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            self.current_molecule = new_molecule
        for root in self.roots_to_follow:
            self.im_drivers[root].qm_data_points = self.qm_data_point_dict[root]
            self.im_drivers[root].compute(new_molecule)
        
        transitions = []
        if len(self.roots_to_follow) > 1:
            for root_1 in range(0, len(self.roots_to_follow)):
                for root_2 in range(root_1 + 1, len(self.roots_to_follow)):

                    potential_kjmol = self.im_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy * hartree_in_kjpermol()
                    potential_kjmol_2 = self.im_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy * hartree_in_kjpermol()
                    transitions.append((self.roots_to_follow[root_1], self.roots_to_follow[root_2], potential_kjmol - potential_kjmol_2))
                    print(f'compare the energies between roots: {self.roots_to_follow[root_1]} -> {self.roots_to_follow[root_2]}', potential_kjmol_2 - potential_kjmol)

                    if 1==2 and np.linalg.norm(self.velocities_np[-1]) > 0.0:
                        current_NAC = self.im_drivers[self.roots_to_follow[root_1]].impes_coordinate.NAC.flatten()
                        current_velocitites = self.velocities_np[-1].flatten() * 4.566180e-4

                        hopping_potential = np.exp(-abs((np.pi/(4)) * (( self.im_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy - self.im_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy ) / np.linalg.multi_dot([current_NAC, current_velocitites]))))
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
    
        self.root_spec_molecules[self.current_state].append(new_molecule)

        potential_kjmol = self.im_drivers[self.current_state].impes_coordinate.energy * hartree_in_kjpermol()
        self.current_gradient = self.im_drivers[self.current_state].impes_coordinate.gradient

        self.current_energy = potential_kjmol


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
    
    def update_potential_energy(self, new_positions):
        """
        Updates and returns the potential energy of the QM region.

        :param new_positions:
            The new positions of the atoms in the QM region.
        :return:
            The potential energy of the QM region.
        """
        _, potential_energy = self.update_gradient_and_energy(new_positions)

        return potential_energy

    def update_forces(self, context):
        """
        Updates the forces in the system based on a new gradient.

        Args:
            context: The OpenMM context object.
        """

        conversion_factor = (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom()) * unit.kilojoule_per_mole / unit.nanometer
        new_positions = context.getState(getPositions=True).getPositions()

        # Update the forces of the QM region
        qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])
        

        gradient = self.update_gradient(qm_positions)
        force = -np.array(gradient) * conversion_factor
        
        # Construct a set from the list of tuples self.broken_bonds
        #broken_bond_atoms = set([atom for bond in self.broken_bonds for atom in bond])

        for i, atom_idx in enumerate(self.qm_atoms):

            self.system.getForce(self.qm_force_index).setParticleParameters(i, atom_idx, force[i])
        self.system.getForce(self.qm_force_index).updateParametersInContext(context)
    
    def get_qm_potential_energy(self):
        """
        Returns the potential energy of the QM region.

        Returns:
            The potential energy of the QM region.
        """

        potential_energy = self.im_drivers[self.current_state].get_energy() * hartree_in_kjpermol()

        return potential_energy
    
    def _calculate_rmsd(self, coord1, coord2):
        """
        Calculate the root mean square deviation (RMSD) between two sets of atomic coordinates
        after optimal alignment using the Kabsch algorithm. This method accounts for rotational
        and translational differences between the molecular conformations.

        :param coord1: 
            The coordinates of the first molecule.
        :param coord2:
            The coordinates of the second molecule.


        :return:
        - rmsd_value: float
            The RMSD between the two aligned molecules.
        """

        # Center the coordinates to eliminate translational differences
        centroid1 = np.mean(coord1, axis=0)  
        centroid2 = np.mean(coord2, axis=0) 
        coord1_centered = coord1 - centroid1
        coord2_centered = coord2 - centroid2

        # The covariance matrix captures how the centered coordinates of the two molecules relate
        covariance_matrix = np.dot(coord1_centered.T, coord2_centered)
        
        # Perform Singular Value Decomposition (SVD) on the covariance matrix
        V, S, Wt = np.linalg.svd(covariance_matrix)

        # Correct for improper rotation (reflection)
        # Calculate the determinant of the rotation matrix to check for reflection
        # A determinant of -1 indicates a reflection, which we need to correct
        d = np.sign(np.linalg.det(np.dot(V, Wt)))

        # If the determinant is negative, adjust the matrices to ensure a proper rotation
        if d < 0:
            # Reflect the last column of V and invert the last singular value
            V[:, -1] *= -1
            S[-1] *= -1

        # Compute the optimal rotation matrix
        # Multiply V and Wt to get the rotation matrix that best aligns coord1 to coord2
        rotation_matrix = np.dot(V, Wt)

        # Apply the rotation to the centered coordinates of the first molecule
        # Rotate coord1_centered using the rotation matrix
        coord1_rotated = np.dot(coord1_centered, rotation_matrix)

        #Calculate the RMSD between the rotated coordinates and the centered coordinates of the second molecule
        diff = coord1_rotated - coord2_centered

        # Square the differences
        diff_squared = diff ** 2

        # Sum the squared differences over all atoms and coordinate axes (x, y, z)
        sum_diff_squared = np.sum(diff_squared)

        # Calculate the mean squared deviation
        mean_diff_squared = sum_diff_squared / len(coord1)

        # Take the square root of the mean squared deviation to obtain the RMSD
        rmsd_value = np.sqrt(mean_diff_squared)

        return rmsd_value

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
        datapoints = []

        for label in sorted_labels:
            impes_coordinate = InterpolationDatapoint(z_matrix)
            impes_coordinate.read_hdf5(datafile, label) # -> read in function from the ImpesDriver object
            coordinates_in_angstrom = impes_coordinate.cartesian_coordinates * bohr_in_angstrom()
            current_molecule = Molecule(mol_labels, coordinates_in_angstrom, 'angstrom') # -> creates a VeloxChem Molecule object
            
            datapoints.append(impes_coordinate)
            data_point_molecules.append(current_molecule)

        return data_point_molecules, datapoints
    
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

           :param current_coordinates:
                current molecular coordinates.
           
           :param data_point:
                InterpolationDatapoint object.

           :returns:
              Norm of the distance between 2 structures.
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

        return distance_core

    
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

    def set_up_the_system(self, molecule, interpolation_settings):

        """
        Assign the neccessary variables with respected values. 

        :param molecule: original molecule

        :param target_dihedrals: is a list of dihedrals that should be scanned during the dynamics

        :param sampling_structures: devides the searchspace around given rotatbale dihedrals
            
        """
        def _determine_cX_symmetry_groups(molecule):
            """
            Build symmetry groups for methyl rotors only.

            Returns the same tuple shape used from AtomMapper:
                (atom_map, symmetry_groups, rc_groups)
            where:
                atom_map: all atom indices
                symmetry_groups: list of equivalent H triplets in CH3 groups
                rc_groups: empty (not used in this CH3-only path)
            """
            labels = [str(x).strip().capitalize() for x in molecule.get_labels()]
            conn = molecule.get_connectivity_matrix()
            n_atoms = len(labels)

            atom_map = list(range(n_atoms))
            rc_groups = []

            ch3_groups = []
            seen = set()

            allowed_neighbours = ["H", "F", "Cl", "Br"]

            for c_idx, symbol in enumerate(labels):
                if symbol != "C":
                    continue

                neighbors = [j for j, bonded in enumerate(conn[c_idx]) if bonded]
                h_neighbors = [j for j in neighbors if labels[j] in allowed_neighbours]
                heavy_neighbors = [j for j in neighbors if labels[j] not in allowed_neighbours]

                # Rotor-relevant methyl carbon: C(H)3-X
                if len(neighbors) == 4 and len(h_neighbors) == 3 and len(heavy_neighbors) == 1:
                    group = tuple(sorted(h_neighbors))
                    if group not in seen:
                        seen.add(group)
                        ch3_groups.append(list(group))

            return atom_map, ch3_groups, rc_groups


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

            def determine_similar_groups(groups_to_process):
                separated_groups = []
                for group in groups_to_process:
                    # Dictionary to group atoms by their exact neighborhood
                    neighbor_map = {}
                    for atom in group:
                        # Find the indices of all atoms connected to this atom
                        neighbors = tuple(sorted([i for i, is_bonded in enumerate(conn[atom]) if is_bonded]))
                        neighbor_map.setdefault(neighbors, []).append(atom)
                    
                    # Extract the automatically separated groups
                    for subgroup in neighbor_map.values():
                        separated_groups.append(sorted(subgroup))
                        
                return separated_groups

            # Apply the separation logic to fix lumped groups before iterating
            groups = determine_similar_groups(groups)


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
        
        def _promote_nonrotatable_ring_torsions_to_impropers(zmat, ff_gen, rotatable_bonds_zero_based):
            """
            Promote proper ring torsions around non-rotatable bonds to impropers.

            Matching by full 4-atom set is too strict (and often fails), so we match
            by the middle bond of the proper torsion and collect impropers touching
            the same bond.
            """
            rot_set = {frozenset((int(i), int(j))) for i, j in rotatable_bonds_zero_based}

            # Existing impropers from z-matrix (if any)
            impropers = [tuple(int(x) for x in imp) for imp in zmat.get("impropers", [])]
            impropers_seen = set(impropers)

            # Build non-rotatable ring bond -> impropers map from MM impropers.
            impropers_ff = [tuple(int(x) for x in t) for t in ff_gen.impropers.keys()]
            bond_to_impropers = {}
            for imp in impropers_ff:
                center = int(imp[0])
                for neigh in imp[1:]:
                    neigh = int(neigh)
                    bond = frozenset((center, neigh))
                    if bond in rot_set:
                        continue
                    if not ff_gen.is_bond_in_ring(center, neigh):
                        continue
                    bond_to_impropers.setdefault(bond, []).append(imp)

            kept_dihedrals = []
            for dih_raw in zmat["dihedrals"]:
                dih = tuple(int(x) for x in dih_raw)
                mid = frozenset((dih[1], dih[2]))

                # Keep proper terms for rotatable or non-ring bonds.
                if mid in rot_set or ff_gen.is_bond_in_ring(dih[0], dih[1]) and ff_gen.is_bond_in_ring(dih[1], dih[2]) and ff_gen.is_bond_in_ring(dih[2], dih[3]):
                    kept_dihedrals.append(dih)
                    continue

                # Non-rotatable ring bond: attach any corresponding impropers.
                promoted = False
                for imp in bond_to_impropers.get(mid, []):
                    if imp not in impropers_seen:
                        impropers.append(imp)
                        impropers_seen.add(imp)
                    promoted = True

                # If no improper was found for this bond, keep the dihedral.
                if not promoted:
                    kept_dihedrals.append(dih)

            zmat["dihedrals"] = kept_dihedrals
            zmat["impropers"] = impropers
            return zmat

        # define global symmetry information object for gs and es
        # used in all other classes for the construction
        # 1. list of all atoms
        # 2. all rotatable CH3 groups for ground-state
        # 3. all groups that are assigned for rotation in the ground-state
        # 4. all atoms that are not in the symmetry groups (CH3-groups excluded)
        # 5. all CH3 symmetry groups
        # 6. indices 0-based of all rotatable bonds
        # 7. list of indices grouped for all dihedrals attached to a rotabale bond (z_matrix index)
        # 8. actual dihedrals with the indices assigned to each atom
        # 9. start and end of the dihedral section within the z_matrix
        self.symmetry_information = {'gs': (), 'es': ()}
        
        self.molecule = molecule
        self.roots_z_matrix = {}
        for root in interpolation_settings.keys():
            
            # generate the z-matrix based for the interpolation database provided
            int_driver = InterpolationDriver()
            int_driver.update_settings(interpolation_settings[root])
    
            _, z_matrix = int_driver.read_labels()
            self.roots_z_matrix[root] = z_matrix

            # The AtomMapper class is based on Turtlemap program which is used to
            # determine equivalent atoms within a molecular structure
            # atom_mapper = AtomMapper(molecule, molecule)
            # symmetry_groups = atom_mapper.determine_symmetry_group()

            symmetry_groups = _determine_cX_symmetry_groups(molecule)

            if not self.use_symmetry:
                symmetry_groups = (symmetry_groups[0], [], symmetry_groups[2])
            
            ff_gen = MMForceFieldGenerator()
            ff_gen.ostream.mute()
            ff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
            ff_gen.create_topology(molecule)

            rotatable_bonds = deepcopy(ff_gen.rotatable_bonds)

            # Work in zero-based indexing (same convention as z-matrix dihedrals)
            # and remove all symmetry-related rotatable bonds from the scan list.
            rotatable_bonds_zero_based = [tuple(sorted((i - 1, j - 1))) for (i, j) in rotatable_bonds]

            self.roots_z_matrix[root] = _promote_nonrotatable_ring_torsions_to_impropers(self.roots_z_matrix[root], ff_gen, rotatable_bonds_zero_based)

            dihedral_start = len(self.roots_z_matrix[root]['bonds']) + len(self.roots_z_matrix[root]['angles'])
            dihedral_end = dihedral_start + len(self.roots_z_matrix[root]['dihedrals'])

            self.all_rotatable_bonds = rotatable_bonds_zero_based
            all_exclision = [element for rot_bond in rotatable_bonds_zero_based for element in rot_bond]
            
            symmetry_groups_ref = [groups for groups in symmetry_groups[1] if not any(item in all_exclision for item in groups)]

            # reduce the symmetry to only CH3 or CH2 symmetry groups for the time being
            regrouped, rot_groups = regroup_by_rotatable_connection(molecule, symmetry_groups_ref, rotatable_bonds_zero_based, molecule.get_connectivity_matrix())
            dih_list = []
            if root == 0 or root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
            
                non_core_atoms = [element for group in regrouped['gs'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                angles_to_set, _, _, self.symmetry_dihedral_lists, dih_list = self._adjust_symmetry_dihedrals(molecule, rot_groups['gs'], rotatable_bonds_zero_based, self.roots_z_matrix[root])
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.roots_z_matrix[root]['dihedrals']):
                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)
                self.symmetry_information['gs'] = [symmetry_groups[0], rot_groups['gs'], regrouped['gs'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [dihedral_start, dihedral_end]]
            if root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip is False:
                non_core_atoms = [element for group in regrouped['gs'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                angles_to_set, _, _, self.symmetry_dihedral_lists, dih_list = self._adjust_symmetry_dihedrals(molecule, rot_groups['gs'], rotatable_bonds_zero_based,  self.roots_z_matrix[root])
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.roots_z_matrix[root]['dihedrals']):
                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)
                self.symmetry_information['es'] = [symmetry_groups[0], rot_groups['es'], regrouped['es'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [dihedral_start, dihedral_end]]
            if root >= 2 and len(self.symmetry_information['es']) == 0 and self.drivers['es']:
                non_core_atoms = [element for group in regrouped['es'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                print(rot_groups['es'], rotatable_bonds_zero_based)
                angles_to_set, _, _, self.symmetry_dihedral_lists, dih_list = self._adjust_symmetry_dihedrals(molecule, rot_groups['es'], rotatable_bonds_zero_based,  self.roots_z_matrix[root])
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.roots_z_matrix[root]['dihedrals']):
                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)
                self.symmetry_information['es'] = [symmetry_groups[0], rot_groups['es'], regrouped['es'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [dihedral_start, dihedral_end]]


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

    def _adjust_symmetry_dihedrals(self, molecule, symmetry_groups, rot_bonds, z_matrix):
        
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
        
        all_dihedrals = [element for element in z_matrix['dihedrals']]

        symmetry_group_dihedral_dict = {} 

        sym_groups_dih = {3: {}}

        angles_to_set = {}
        periodicities = {}
        dihedral_groups = {2: {}, 3: {}}

        for symmetry_group in symmetry_groups:
             
            symmetry_group_dihedral_list = symmetry_group_dihedral(symmetry_group, all_dihedrals, rot_bonds)
            symmetry_group_dihedral_dict[tuple(symmetry_group)] = symmetry_group_dihedral_list
            
            if len(symmetry_group) == 3:
                center = tuple(int(x) for x in symmetry_group_dihedral_list[0][1:3])
                atom_group = tuple(int(x) for x in symmetry_group)

                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = (
                    [0.0, np.pi / 6.0, np.pi / 3.0, np.pi / 2.0]
                )
                periodicities[tuple(symmetry_group_dihedral_list[0])] = 3

                dihedral_groups[3][center] = [
                    tuple(sorted(element, reverse=False))
                    for element in symmetry_group_dihedral_list
                ]

                center_key = tuple(int(x) for x in symmetry_group_dihedral_list[0][1:3])
                atom_group_key = tuple(int(x) for x in symmetry_group)

                sym_groups_dih[3][(center_key, atom_group_key)] = [
                    tuple(int(x) for x in element)
                    for element in symmetry_group_dihedral_list
                ]

            elif len(symmetry_group) == 2:

                # angles_to_set[symmetry_group_dihedral_list[0]] = ([0.0, np.pi/3.0])
                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = ([0.0, np.pi/2.0])

                periodicities[tuple(symmetry_group_dihedral_list[0])] = 2
                dihedral_groups[2].extend([tuple(sorted(element, reverse=False)) for element in symmetry_group_dihedral_list])

        return angles_to_set, periodicities, symmetry_group_dihedral_dict, dihedral_groups, sym_groups_dih