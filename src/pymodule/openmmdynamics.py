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
from pathlib import Path
from sys import stdout
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom

# TODO: added try import openmm to ensure that the
# python tests are run even if openmm is not installed.
# Are there better solutions?
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass
from .molecule import Molecule
from .veloxchemlib import mpi_master
from. veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .forcefieldgenerator import ForceFieldGenerator
from .systembuilder import SystemBuilder

# Drivers
#from .xtbdriver import XtbDriver
#from .impesdriver import ImpesDriver
from .optimizationdriver import OptimizationDriver

class OpenMMDynamics:
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
        self.padding = 2.0
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

        # QM Region parameters
        self.qm_driver = None
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
    # TODO: Integrate the guess with the read_pdb_file in Molecule.
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
            sys_builder.steps = 10000
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
            self.qm_stabilizer(ff_gen)
        
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
        else:
            periodic = False

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
        qm_residue_index = self.unique_residues.index((qm_residue_name, qm_residue_number))
        msg = f'QM Residue: {qm_residue_name}'
        self.ostream.print_info(msg)
        self.ostream.flush()

        qm_molecule = self.unique_molecules[qm_residue_index]

        # Generate or use an existing forcefield generator for the QM region
        if ff_gen_qm is None:
            ff_gen_qm = ForceFieldGenerator()
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

    def qm_minimization(self, molecule, scf_drv, constraints=None):
        """
        Minimizes the energy using a QM driver.

        :param molecule:
            VeloxChem Molecule object.
        :param scf_drv:
            VeloxChem SCF driver object.
        :param constraints:
            Constraints for the optimization. Default is None. Format example: ["set dihedral 6 1 2 3 0.0", "freeze distance 6 1", ...]
        """

        opt_drv = OptimizationDriver(scf_drv)
        opt_drv.ostream.mute()

        if constraints:
            opt_drv.constraints(constraints)

        results = opt_drv.compute(molecule)

        min_energy = results['opt_energies'][-1]
        min_coords = results['opt_geometries'][-1]

        return min_energy, min_coords
     
    def conformational_sampling(self, 
                                ensemble='NVT', 
                                temperature=700, 
                                timestep=2.0, 
                                nsteps=10000, 
                                snapshots=10,
                                lowest_conformations=None,
                                qm_driver=None,
                                constraints=None,
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
        :param lowest_conformations:
            Number of lowest energy conformations to save. Default is None.
        :param qm_minimization:
            QM driver object for energy minimization. Default is None. 
        :param minimize:
            If True, the energy of the conformations will be minimized. Default is True.

        :return:
            energies: List of potential energies of the conformations.
            opt_coordinates: List of XYZ format strings of the minimized coordinates.
        """

        if self.system is None:
            raise RuntimeError('System has not been created! First create a system using the appropriate method.')

        self.ensemble = ensemble
        self.temperature = temperature * unit.kelvin
        self.timestep = timestep * unit.femtosecond
        self.nsteps = nsteps

        self.integrator = self._create_integrator()
        topology = self.pdb.topology
        self.positions = self.pdb.positions

        self.simulation = app.Simulation(topology, self.system, self.integrator)
        self.simulation.context.setPositions(self.positions)
        self.simulation.context.setVelocitiesToTemperature(self.temperature)

        if self.molecule is None:
            if len(self.unique_molecules) > 1:
                # There are several residues in the system, print a warning because the conformational sampling may not be meaningful
                warning_msg = 'You are using a PDB file to run the simulation. The coordinates saved will correspond to the whole PDB file.'
                warning_msg += 'Make sure you have one residue in your PDB file or the results may not be meaningful.'
                self.ostream.print_warning(warning_msg)
                self.ostream.flush()
            # This is the case where the user has not created a molecule object and used files to create the system.
            self.labels = [atom.element.symbol for atom in self.pdb.topology.atoms()]

        # Determine the frequency of saving depending on the number of snapshots.
        save_freq = max(nsteps // snapshots, 1)
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

                    # Create a xyz with the minimized coordinates
                    energies.append(minimized_energy)
                    print(f'Minimized energy: {minimized_energy}')
                    xyz = f"{len(self.labels)}\n\n"
                    for label, coord in zip(self.labels, minimized_coordinates):
                        xyz += f"{label} {coord.x * 10} {coord.y * 10} {coord.z * 10}\n"  
                    print('Saved coordinates for step', step)
                    opt_coordinates.append(xyz)

                else:
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

        if lowest_conformations:
            msg = f'Looking for the {lowest_conformations} lowest energy conformations.'
            self.ostream.print_info(msg)
            self.ostream.flush()
            # Get the lowest energy conformations based on the number of conformations requested
            lowest_energies = sorted(energies)[:lowest_conformations]
            lowest_coords = opt_coordinates[:lowest_conformations]
            
            energies = lowest_energies
            opt_coordinates = lowest_coords

            msg = 'Lowest energy conformations saved'
            self.ostream.print_info(msg)
            self.ostream.flush()

        if qm_driver:
            # Use qm_miniization to minimize the energy of the conformations
            # Name flag based on if is instance of the XtbDriver or ScfDriver
            if isinstance(qm_driver, XtbDriver):
                drv_name = 'XTB Driver'
            else:
                drv_name = 'SCF Driver'

            msg = f'Requested QM minimization with {drv_name}'
            self.ostream.print_info(msg)
            self.ostream.flush()
            qm_energies = []
            qm_opt_coordinates = []
            for conf , coords in enumerate(opt_coordinates):
                print(f'Conformation {conf}')
                molecule = Molecule.from_xyz_string(coords)
                qm_energy, qm_coords = self.qm_minimization(molecule, qm_driver, constraints)
                qm_energies.append(qm_energy)
                qm_opt_coordinates.append(qm_coords)
                msg = f'Energy of the conformation: {qm_energy}'
                self.ostream.print_info(msg)
                self.ostream.flush()

            energies = qm_energies
            opt_coordinates = qm_opt_coordinates

        # Convert the opt_coordinates to molecule objects
        opt_molecules = [Molecule.from_xyz_string(coords) for coords in opt_coordinates]

        return energies, opt_molecules

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
        self.simulation.minimizeEnergy()

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

    def run_qmmm(self, 
                 qm_driver, 
                 grad_driver,
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
        
        self.ensemble = ensemble
        self.temperature = temperature * unit.kelvin
        self.friction = friction / unit.picosecond
        self.timestep = timestep * unit.femtoseconds
        self.nsteps = nsteps

        self.qm_driver = qm_driver
        self.grad_driver = grad_driver
        self.qm_driver.ostream.mute()

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
            else:
                self.simulation.context.setVelocitiesToTemperature(10 * unit.kelvin)

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

        for step in range(nsteps):

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


            # Information output
            if step % save_freq == 0:
                
                print(f"Step: {step} / {nsteps} Time: {round((step * timestep) / 1000, 2)} ps")
                print('Potential Energy QM region:', qm, 'kJ/mol')
                print('Potential Energy MM region:', mm)
                print('QM/MM Interaction Energy:', qm_mm)
                print('Total Potential Energy:', pot)
                print('Kinetic Energy:', kinetic)
                print('Temperature:', temp, 'K')
                print('Total Energy:', total)
                print('-' * 60)   

            self.simulation.step(1)

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
            ForceFieldGenerator object from VeloxChem.
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
            ForceFieldGenerator object from VeloxChem.
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

        else:
            # Atom labels for the QM region
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")

        if not isinstance(self.qm_driver, ImpesDriver):
            # Compute explicitly the energy and gradient
            self.qm_driver.compute(new_molecule)
            self.grad_driver.compute(new_molecule)
            potential_kjmol = self.qm_driver.get_energy() * hartree_in_kcalpermol() * 4.184
            gradient = self.grad_driver.get_gradient()
        else:
            # IM driver has the energy and gradient in the same driver
            self.qm_driver.compute(new_molecule)
            potential_kjmol = self.qm_driver.get_energy() * hartree_in_kcalpermol() * 4.184
            gradient = self.qm_driver.get_gradient()

        return gradient, potential_kjmol

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

        conversion_factor = (4.184 * hartree_in_kcalpermol() * 10.0 / bohr_in_angstrom()) * unit.kilojoule_per_mole / unit.nanometer
        new_positions = context.getState(getPositions=True).getPositions()

        # Update the forces of the QM region
        qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])

        gradient = self.update_gradient(qm_positions)
        force = -np.array(gradient) * conversion_factor

        custom_force = self.system.getForce(self.qm_force_index)

        # Construct a set from the list of tuples self.broken_bonds
        #broken_bond_atoms = set([atom for bond in self.broken_bonds for atom in bond])

        for i, atom_idx in enumerate(self.qm_atoms):
            custom_force.setParticleParameters(i, atom_idx, force[i])
            
        custom_force.updateParametersInContext(context)
    
    def get_qm_potential_energy(self):
        """
        Returns the potential energy of the QM region.

        Returns:
            The potential energy of the QM region.
        """

        potential_energy = self.qm_driver.get_energy() * hartree_in_kcalpermol() * 4.184

        return potential_energy
    




