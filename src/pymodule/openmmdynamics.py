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
from itertools import combinations, accumulate
import string
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kjpermol, bohr_in_angstrom
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream
from .mmforcefieldgenerator import MMForceFieldGenerator
from .solvationbuilder import SolvationBuilder
from .xtbdriver import XtbDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .optimizationdriver import OptimizationDriver
from .interpolationdriver import InterpolationDriver
from .errorhandler import assert_msg_critical
from .mofutils import svd_superimpose

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass


class OpenMMDynamics:
    """
    Implements the OpenMMDynamics.
    Performs classical MD and QM/MM simulations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - openmm_platform: The platform for OpenMM.
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
        - qm_driver: The VeloxChem driver object. Options are XtbDriver or an ScfDriver.
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

        # MPI and output stream
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # openmm platform
        self.openmm_platform = None

        # Instance variables
        # Simulation parameters
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
        self.qm_driver = None
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
    def create_system_from_molecule(self, 
                             molecule, 
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

    # Methods to build a custom system from a PDB file and custom XML files
    def create_md_system_from_files(self, 
                            pdb_file, 
                            xml_file,
                            filename='custom',
                            write_files = True):
        """
        Creates a system from a PDB file containing multiple residues and custom XML files.
        
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

        self.system = forcefield.createSystem(self.pdb.topology, 
                                              **system_arguments)

        if write_files:
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

    def create_qmmm_system_from_files(self, 
                            pdb_file, 
                            xml_file, 
                            qm_residue, 
                            ff_gen_qm=None, 
                            qm_atoms='all', 
                            filename='custom'):
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
            ff_gen_qm = MMForceFieldGenerator()
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
        # self.qm_stabilizer(ff_gen_qm)

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

    # Simulation methods
    def conformational_sampling_multiple(self,
                                molecules = None,
                                pdb_files=None,
                                xml_files=None,
                                partial_charges = None,
                                position_restraint = None,
                                temperature=700, 
                                timestep=2.0, 
                                nsteps=500000,
                                snapshots=100,
                                lowest_conformations=10,
                                save_pdb = False,
                                k = 1000,
                                r0 = 0.5,
                                water_model = 'ctip3p'
                                ):
        """
        Runs high-temperature conformational sampling for multiple residues in the system.

        :param molecules:
            List of VeloxChem Molecule objects to be used in the simulation.
        :param pdb_files:
            List of PDB files for each molecule to be used in the simulation.
        :param xml_files:
            List of XML files for the molecule force fields if pdb_files are used.
        :param partial_charges:
            List (of lists) of partial charges for all molecules. If None, RESP charges will be computed.
        :param position_restraint:
            List of 1-based indices of the molecule(s)/PDB file(s) in the input list whose positions will be constrained.
        :param temperature:
            Temperature of the system in Kelvin. Default is 700 K.
        :param timestep:
            Timestep of the simulation in femtoseconds. Default is 2.0 fs.
        :param nsteps:
            Number of steps in the simulation. Default is 500 000.
        :param snapshots:
            The number of snapshots to save. Default is 100.
        :param lowest_conformations:
            Number of lowest energy conformations to save. Default is 10.
        :param save_pdb:
            If true, a PDB file of the lowest energy conformers is written. Default is False.
        :param k:
            Force constant of the Centroid bond force. Default is 1000.
        :param r0:
            Equilibrium distance of the Centroid bond force. Default i 0.5 nm.
        :param water_model:
            The water model to be used if water molecules are present.

        :return:
            conformers_dict: Dictionary with lists of potential energies of the conformations, the minimized molecule objects, 
            and their corresponding coordinates in XYZ format.
        """
        if molecules:
            self.ostream.print_info("Generating system...")
            self.ostream.flush()
            self.atom_dict = {}
            xml_files = []
            water_msg = False

            for i, mol in enumerate(molecules):
                ff_gen = MMForceFieldGenerator()
                ff_gen.ostream.mute()
                duplicate = False
                is_water = False
                
                if mol.is_water_molecule():
                    is_water = True
                    water_msg = True

                for j in range(i):
                    if (mol.get_labels() == molecules[j].get_labels() and
                        mol.get_connectivity_matrix().shape == molecules[j].get_connectivity_matrix().shape and
                        (mol.get_connectivity_matrix() == molecules[j].get_connectivity_matrix()).all()):
                        ff_gen.create_topology(mol, resp=False, water_model=water_model)
                        self.atom_dict[f'{i}'] = ff_gen.atoms
                        duplicate = True
                        break

                if not duplicate:
                    if partial_charges:
                        ff_gen.partial_charges = partial_charges[i]
                    ff_gen.create_topology(mol, water_model=water_model)
                    ff_gen.generate_residue_xml(f'molecule_{i+1}.xml', f'M{i+1:02d}')
                    xml_files.append(f'molecule_{i+1}.xml')
                
                self.atom_dict[f'{i}'] = ff_gen.atoms
            
            if water_msg:    
                self.ostream.print_info(f'Water molecule detected. Using {water_model} water model.')
                self.ostream.flush()
            
            pdb_file = 'system.pdb'
            use_chain = False
            self._create_system_from_multiple_molecules(molecules, pdb_file, write_pdb = True)
            self.create_md_system_from_files(pdb_file, xml_files, write_files=False)
            self.topology = self.pdb.topology
            self.positions = self.pdb.positions
              
        elif pdb_files:
            # creating a system from merged xml files
            assert_msg_critical(
                xml_files is not None,
                "No XML files provided. Please provide either a list of VeloxChem Molecule objects or the XML files.")   
            
            self.ostream.print_info("Merging pdb files...")
            self.ostream.flush()
            
            # to ensure that coordinates are adjusted if identical files are provided
            molecules = [Molecule.read_pdb_file(pdb) for pdb in pdb_files]
            self._create_system_from_multiple_molecules(molecules, write_pdb = False)
            
            pdbs = [app.PDBFile(pdb) for pdb in pdb_files]
            coords = self.all_coords # list of the new coordinates
            atom_counter = [len(mol.get_labels()) for mol in molecules]

            modeller = None
            for pdb, (start,end) in zip(pdbs, zip([0,*accumulate(atom_counter[:-1])], accumulate(atom_counter))):
                positions = coords[start:end] * unit.angstrom
                if modeller is None:
                    modeller = app.Modeller(pdb.topology, positions)  
                else:
                    modeller.add(pdb.topology, positions)
            
            pdb_file = 'merged_system.pdb'
            use_chain = True

            app.PDBFile.writeFile(modeller.topology, modeller.positions, pdb_file)
            
            # create system
            self.topology = modeller.topology
            self.positions = modeller.positions
            forcefield = app.ForceField(*xml_files)
            self.ostream.print_info(f'Added XML Files: {xml_files}')
            self.ostream.flush()
            
            system_arguments = {
                'nonbondedMethod': app.NoCutoff,
                'constraints': app.HBonds
            }   
            
            self.system = forcefield.createSystem(modeller.topology, **system_arguments)

        # prepare system for simulation
        real_system = self.system # used for recalculation
        self.ensemble = 'NVT'
        self.temperature = temperature * unit.kelvin
        self.timestep = timestep * unit.femtosecond
        self.nsteps = nsteps

        self.integrator = self._create_integrator()
        
        if pdb_files:
            self.labels = [atom.element.symbol for atom in self.topology.atoms()]

        # add customcentroidbondforce to the system
        self._get_centroid_bond_force(k, r0, use_chain, position_restraint)
    
        self.simulation = app.Simulation(self.topology, self.system, self.integrator, platform=self._create_platform())
        self.simulation.context.setPositions(self.positions)
        self.simulation.context.setVelocitiesToTemperature(self.temperature)
        self.simulation.minimizeEnergy()

        # Determine the frequency of saving depending on the number of snapshots.
        save_freq = max(nsteps // snapshots, 1)
        conformations = []
        minimized_energies = []
        opt_coordinates = []
        pdb_coords = []

        # Run MD 
        self.ostream.print_info(f"Running high-temperature MD for {self.nsteps * self.timestep.value_in_unit(unit.nanoseconds):.2f} ns...")
        self.ostream.flush()
        for i in range(snapshots):
            self.simulation.step(save_freq)
            state = self.simulation.context.getState(getEnergy=True, getPositions=True)
            positions = state.getPositions().value_in_unit(unit.nanometers)
            self.ostream.print_info(f'Saved coordinates for step {save_freq * (i + 1)}')
            self.ostream.print_info(f"Energy for conformer {i+1}/{snapshots}: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.4f} kJ/mol")
            self.ostream.flush()
            conformations.append(positions)
        
        # Recalculate energies
        self.ostream.print_info("Recalculating energies for the conformations...")
        self.ostream.flush()
        simulation = app.Simulation(self.topology, real_system, self._create_integrator(), platform=self._create_platform())

        for idx, conformation in enumerate(conformations):    
            simulation.context.setPositions(conformation)
            simulation.minimizeEnergy()
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            
            minimized_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            minimized_energies.append(minimized_energy)
            minimized_coordinates = state.getPositions()
            pdb_coords.append(minimized_coordinates)
            self.ostream.print_info(f'Minimized energy of conformer {idx}: {minimized_energy:.4f} kJ/mol')
            self.ostream.flush()

            xyz = f"{len(self.labels)}\n\n"
            for label, coord in zip(self.labels, minimized_coordinates):
                xyz += f"{label} {coord.x * 10} {coord.y * 10} {coord.z * 10}\n"  

            opt_coordinates.append(xyz)

        # Sort the conformations by (increasing) energy
        index_ẹnergy = [(i, minimized_energies[i]) for i in range(len(minimized_energies))]
        sorted_index_energy = sorted(index_ẹnergy, key=lambda x: x[1]) 
        sorted_indices = [i[0] for i in sorted_index_energy]

        minimized_energies = [minimized_energies[i] for i in sorted_indices]
        opt_coordinates = [opt_coordinates[i] for i in sorted_indices]

        if lowest_conformations:
            self.ostream.print_info(f"Saving the {lowest_conformations} lowest energy conformations")
            self.ostream.flush()
            minimized_energies = minimized_energies[:lowest_conformations]
            opt_coordinates = opt_coordinates[:lowest_conformations]
        
        self.ostream.print_info('Conformational sampling completed!')
        self.ostream.flush()

        if save_pdb:
            filename = 'lowest_conformer.pdb'
            self.ostream.print_info(f"Saving PDB file of the lowest energy conformations as {filename}")
            self.ostream.flush()
            app.PDBFile.writeFile(self.topology, pdb_coords[sorted_indices[0]], filename)


        # Save final molecules, coordinates and corresponding energies to a dictionary
        conformers_dict = {
                'energies': minimized_energies,
                'molecules': [Molecule.from_xyz_string(coords) for coords in opt_coordinates],
                'geometries': opt_coordinates
        }

        self.conformer_dict = conformers_dict

        return conformers_dict

    def conformational_sampling(self, 
                                ensemble='NVT', 
                                temperature=700, 
                                timestep=2.0, 
                                nsteps=10000, 
                                snapshots=10,
                                unique_conformers=True,
                                qm_driver=None,
                                basis=None,
                                constraints=None):

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
        :param unique_conformers:
            If True, the method will remove any identical conformers. Default is True.
        :param qm_minimization:
            QM driver object for energy minimization. Default is None.
        :param basis:
            Basis set for the SCF driver if qm_driver is not None.
        :param constraints:
            Constraints for the system. Default is None.
        :param minimize:
            If True, the energy of the conformations will be minimized. Default is True.

        :return:
            conformers_dict: Dictionary with lists of potential energies of the conformations, the minimized molecule objects, 
            and their corresponding coordinates in XYZ format.
            
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

        self.simulation = app.Simulation(topology, self.system, self.integrator, platform=self._create_platform())
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

        # Run MD
        for i in range(snapshots):
            self.simulation.step(save_freq)
        
            minimized_system = app.Simulation(topology, self.system, self._create_integrator(), platform=self._create_platform())
            minimized_system.context.setPositions(self.simulation.context.getState(getPositions=True).getPositions())
            
            minimized_system.minimizeEnergy()
            minimized_state = minimized_system.context.getState(getPositions=True, getEnergy=True)
            
            minimized_coordinates = minimized_state.getPositions()
            minimized_energy = minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            energies.append(minimized_energy)
            
            self.ostream.print_info(f'Minimized energy: {minimized_energy}')
            self.ostream.flush()
            xyz = f"{len(self.labels)}\n\n"
            for label, coord in zip(self.labels, minimized_coordinates):
                xyz += f"{label} {coord.x * 10} {coord.y * 10} {coord.z * 10}\n"  
            self.ostream.print_info(f'Saved coordinates for step {save_freq * (i + 1)}')
            self.ostream.flush()
            opt_coordinates.append(xyz)

        self.ostream.print_info('Conformational sampling completed!')
        self.ostream.print_info(f'Number of conformations: {len(opt_coordinates)}')
        self.ostream.flush()

        equiv_conformer_pairs= []
        if unique_conformers:
            msg = f'Filtering for unique conformers'
            self.ostream.print_info(msg)
            self.ostream.flush()
            
            # reorder coordinates and energies by increasing energy
            index_ẹnergy = [(i, energies[i]) for i in range(len(energies))]
            sorted_index_energy = sorted(index_ẹnergy, key=lambda x: x[1]) # sorted by increasing energy
            sorted_indices = [i[0] for i in sorted_index_energy]


            minimized_energy = [energies[i] for i in sorted_indices]
            opt_coordinates = [opt_coordinates[i] for i in sorted_indices]


            # Filter out the conformations with RMSD < 0.1 and keep the lowest energy one.
            for i, coord in enumerate(opt_coordinates):
                mol= Molecule.read_xyz_string(coord)
                xyz_i = mol.get_coordinates_in_angstrom()
                ene_i = minimized_energy[i]

                for j in range(i + 1, len(opt_coordinates)):
                    mol_j= Molecule.read_xyz_string(opt_coordinates[j])
                    xyz_j = mol_j.get_coordinates_in_angstrom()
                    ene_j = minimized_energy[j]
                    if abs(ene_i - ene_j) < self.energy_threshold:
                        rmsd, rot, trans = svd_superimpose(xyz_j, xyz_i)
                        if rmsd < self.rmsd_threshold:
                            equiv_conformer_pairs.append((i, j))

            duplicate_conformers = [j for i, j in equiv_conformer_pairs]
            duplicate_conformers = sorted(list(set(duplicate_conformers)))

            filtered_energies = [
                e
                for i, e in enumerate(minimized_energy)
                if i not in duplicate_conformers
            ]

            filtered_geometries = [
                g
                for i, g in enumerate(opt_coordinates)
                if i not in duplicate_conformers
            ]
            
            opt_coordinates = filtered_geometries
            energies = filtered_energies


        
        
            # print table of results
            msg = f'\nNumber of unique conformers: {len(opt_coordinates)}'
            self.ostream.print_info(msg)
            self.ostream.flush()

        if qm_driver:
            # Use qm_miniization to minimize the energy of the conformations
            # Name flag based on if is instance of the XtbDriver or ScfDriver
            if isinstance(qm_driver, XtbDriver):
                drv_name = 'XTB Driver'
            elif isinstance(qm_driver, ScfRestrictedDriver):
                drv_name = 'RSCF Driver'
            elif isinstance(qm_driver, ScfUnrestrictedDriver):
                drv_name = 'USCF Driver'
            elif isinstance(qm_driver, ScfRestrictedOpenDriver):
                drv_name = 'ROSCF Driver'

            if drv_name != 'XTB Driver' and basis is None:
                raise ValueError('Basis set is required for the SCF driver.')

            msg = f'Requested QM minimization with {drv_name}'
            self.ostream.print_info(msg)
            self.ostream.flush()
            qm_energies = []
            qm_opt_coordinates = []
            for conf , coords in enumerate(opt_coordinates):
                self.ostream.print_info(f'Conformation {conf}')
                self.ostream.flush()
                molecule = Molecule.from_xyz_string(coords)
                qm_energy, qm_coords = self.qm_minimization(molecule, basis, qm_driver, constraints)
                qm_energies.append(qm_energy)
                qm_opt_coordinates.append(qm_coords)
                msg = f'Energy of the conformation: {qm_energy}'
                self.ostream.print_info(msg)
                self.ostream.flush()

            energies = qm_energies
            opt_coordinates = qm_opt_coordinates

        # Save final molecules, coordinates and corresponding energies to a dictionary
        conformers_dict = {
                'energies': energies,
                'molecules': [Molecule.from_xyz_string(coords) for coords in opt_coordinates],
                'geometries': opt_coordinates
        }

        self.conformer_dict = conformers_dict

        return conformers_dict
    
    def calculate_boltzmann_distribution(self, energies=None, T=300, unit='kj/mol'):
        """ 
        Calculate the Boltzmann distribution of conformers based on their energies.
        
        :param energies: 
            list of energies. If None, it will use the energies from the conformers_dict from the conformational sampling.
        :param T: 
            temperature in Kelvin
        :param unit: 
            unit of energy, options are 'kj/mol', 'kcal/mol', 'hartree'

        :return: 
            list of probabilities for each conformer
        
        """
        if unit=='kj/mol':
            R = 0.008314462618
        elif unit=='kcal/mol':
            R = 1.98720425864083e-3
        elif unit=='hartree':
            R = 3.166811563671e-6
        else:
            raise ValueError('Invalid unit')
        # calculate relative energies
        if energies is None:
            energies = self.conformer_dict['energies']
        else:
            energies = energies

        relative_energies = [energy - min(energies) for energy in energies]

        # calculate the boltzmann factors
        boltzmann_factors = [np.exp(-energy/(R*T)) for energy in relative_energies]
        # calculate the partition function
        partition_function = sum(boltzmann_factors)
        # calculate the probabilities
        probabilities = [factor/partition_function for factor in boltzmann_factors]
        
        return probabilities
    
    def show_conformers(self, number=5, atom_indices=False, atom_labels=False, boltzmann_distribution=True):
        weights = None
        if boltzmann_distribution:
            weights = self.calculate_boltzmann_distribution(T=300, unit='kj/mol')

        if number > len(self.conformer_dict["energies"]):
            number = len(self.conformer_dict["energies"])
            print(f"Only {number} conformers available, showing all.")
        for i in range(number):
            if weights is not None:
                msg = f'\nConformation {i+1}: Energy: {self.conformer_dict["energies"][i]:.3f} kJ/mol, Weight: {weights[i]:.4f}'   
            else:
                msg = f'\nConformation {i+1}: Energy: {self.conformer_dict["energies"][i]:.3f} kJ/mol'
            self.ostream.print_info(msg)
            self.ostream.flush()
            self.conformer_dict["molecules"][i].show(
                atom_indices=atom_indices, atom_labels=atom_labels)
                
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
        
        self.simulation = app.Simulation(self.topology, self.system, new_integrator, platform=self._create_platform())

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
                dof = self.system.getNumParticles() * 3 - self.system.getNumConstraints()
                if hasattr(self.integrator, 'computeSystemTemperature'):
                    temp = self.integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
                else:
                    temp = (2*state.getKineticEnergy()/(dof * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin) 

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
                 basis = None,
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
        
        if basis is not None:
            self.basis = basis

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
        elif isinstance(self.qm_driver, ScfRestrictedDriver):
            self.driver_flag = 'RSCF Driver'
        elif isinstance(self.qm_driver, ScfUnrestrictedDriver):
            self.driver_flag = 'USCF Driver'
        elif isinstance(self.qm_driver, ScfRestrictedOpenDriver):
            self.driver_flag = 'ROSCF Driver'
        elif isinstance(self.qm_driver, InterpolationDriver):
            self.driver_flag = 'IM Driver'
            _ = self.qm_driver.read_qm_data_points()

        
        else:
            raise ValueError('Invalid QM driver. Please use a valid VeloxChem driver.')

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
        
        self.simulation = app.Simulation(self.topology, self.system, new_integrator, platform=self._create_platform())

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

            # Kinetic energyex
            
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
    
    # Other methods:
    
    def qm_minimization(self, molecule, basis, scf_drv, constraints=None):
        """
        Minimizes the energy using a QM driver.

        :param molecule:
            VeloxChem Molecule object.
        :param basis:
            Basis set for the SCF driver.
        :param scf_drv:
            VeloxChem SCF driver object.
        :param constraints:
            Constraints for the optimization. Default is None. Format example: ["set dihedral 6 1 2 3 0.0", "freeze distance 6 1", ...]
        """

        opt_drv = OptimizationDriver(scf_drv)
        opt_drv.ostream.mute()

        if constraints:
            opt_drv.constraints(constraints)

        basis = MolecularBasis.read(molecule, basis)

        scf_results = scf_drv.compute(molecule, basis)

        results = opt_drv.compute(molecule, basis, scf_results)

        min_energy = results['opt_energies'][-1]
        min_coords = results['opt_geometries'][-1]

        return min_energy, min_coords
    
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

    def _create_platform(self):
        """
        Creates an OpenMM platform.

        Returns:
            OpenMM Platform.
        """

        if self.openmm_platform is None:
            return None
        else:
            platform = mm.Platform.getPlatformByName(self.openmm_platform)
            if self.openmm_platform == "CPU":
                platform.setPropertyDefaultValue("Threads", "1")
            return platform

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
            self.dynamic_molecules.append(new_molecule)

        else:
            # Atom labels for the QM region
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            self.dynamic_molecules.append(new_molecule)

        if self.basis is not None:
            basis = MolecularBasis.read(new_molecule, self.basis)
            scf_results = self.qm_driver.compute(new_molecule, basis)
            gradient = self.grad_driver.compute(new_molecule, basis, scf_results)
            potential_kjmol = self.qm_driver.get_scf_energy() * hartree_in_kjpermol()
            gradient = self.grad_driver.get_gradient()
        else:
            self.qm_driver.compute(new_molecule)
            if self.driver_flag != 'IM Driver':
                self.grad_driver.compute(new_molecule)
            potential_kjmol = self.qm_driver.get_energy() * hartree_in_kjpermol()
            gradient = self.grad_driver.get_gradient()

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

        conversion_factor = (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom()) * unit.kilojoule_per_mole / unit.nanometer
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
        if self.basis is not None:
            potential_energy = self.qm_driver.get_scf_energy() * hartree_in_kjpermol()
        else:
            potential_energy = self.qm_driver.get_energy() * hartree_in_kjpermol()

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

    def _get_centroid(self,coords):
        centroid = np.mean(coords, axis=0)
        rg_sq = np.mean(np.sum((coords - centroid) ** 2, axis=1))
        
        return np.sqrt(rg_sq)

    def _create_system_from_multiple_molecules(
        self, molecules, pdb_file='system.pdb', spacing_factor=1.5, max_attempts=1000, write_pdb = True
    ):
        placed_coords = []
        placed_rgs = []
        all_labels = []
        atom_mol_ids = []

        for mol_index, mol in enumerate(molecules):
            coords = mol.get_coordinates_in_angstrom()
            labels = mol.get_labels()
            centroid = np.mean(coords, axis=0)
            coords_centered = coords - centroid
            
            rg = self._get_centroid(coords)

            if mol_index == 0:
                offset = np.zeros(3)
                coords_shifted = coords_centered
            else:
                # Try finding a non-overlapping position
                for attempt in range(max_attempts):
                    vec = np.random.normal(size=3)
                    direction = vec / np.linalg.norm(vec)
                    distance = 0.0

                    # Calculate min distance needed from all placed molecules
                    for prev_rg, prev_coords in zip(placed_rgs, placed_coords):
                        dist_needed = spacing_factor * (rg + prev_rg)
                        distance = max(distance, dist_needed * (1 + 0.1 * attempt))

                    offset = direction * distance
                    coords_shifted = coords_centered + offset

                    # Check distance from all placed molecule centroids
                    overlaps = False
                    for prev_coords in placed_coords:
                        prev_centroid = np.mean(prev_coords, axis=0)
                        new_centroid = np.mean(coords_shifted, axis=0)
                        if np.linalg.norm(new_centroid - prev_centroid) < spacing_factor * (rg + self._get_centroid(prev_coords)):
                            overlaps = True
                            break

                    if not overlaps:
                        break
                else:
                    raise RuntimeError(f"Failed to place molecule {mol_index} without overlaps after {max_attempts} attempts.")

            placed_coords.append(coords_shifted)
            placed_rgs.append(rg)
            all_labels.extend(labels)
            atom_mol_ids.extend([mol_index] * len(labels))
        
        self.labels = all_labels
        all_coords = np.vstack(placed_coords)
        self.all_coords = all_coords
        
        if write_pdb:
            chain_id_list = list(string.ascii_uppercase + string.ascii_lowercase + string.digits)
            atom_names = [atom_info['name'] for mol in self.atom_dict.values() for atom_info in mol.values()]
            # Determine box size
            min_coords = np.min(all_coords, axis=0)
            max_coords = np.max(all_coords, axis=0)
            extent = max_coords - min_coords
            
            # Add the padding to the extent
            self.padding = 20 # (Å)
            box_size = np.max(extent) + 2 * self.padding

            with open(pdb_file, 'w') as f:
                # Write the box size (Ångstrom) to the PDB file
                f.write("HEADER    Generated by VeloxChem\n")
                f.write(f"CRYST1{box_size:9.3f}{box_size:9.3f}{box_size:9.3f}  90.00  90.00  90.00 P 1           1\n")

                for i, (coord, label, mol_id) in enumerate(zip(all_coords, all_labels, atom_mol_ids), 1):
                    resname = f"M{mol_id+1:02d}"
                    chain_id = chain_id_list[mol_id % len(chain_id_list)]
                    resseq = mol_id + 1
                    atom_name = atom_names[i - 1][:4]  # Ensure name is max 4 chars
                    occupancy = 1.00
                    temp_factor = 0.00
                    element_symbol = label[:2].rjust(2).capitalize()

                    # Format line per PDB spec (https://cupnet.net/pdb-format/)
                    line_str = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   ".format(
                        'HETATM', i, atom_name, '', resname[:3], chain_id, resseq, '')

                    line_str += "{:8.3f}{:8.3f}{:8.3f}".format(
                        coord[0], coord[1], coord[2])

                    line_str += "{:6.2f}{:6.2f}          {:>2s}".format(
                        occupancy, temp_factor, element_symbol)

                    f.write(line_str + '\n')

                f.write("END\n")

    def _get_centroid_bond_force(self, k, r0, use_chain = False, position_restraint = None):
        """
        Creates a CustomCentroidBondForce to restrain the molecules to eachother during conformational sampling
        of more than one molecule.
        """
        self.ostream.print_info("Applying centroid bond force to the system...")
        self.ostream.flush()

        residues = list(self.topology.chains()) if use_chain else list(self.topology.residues())
        if len(residues) < 2:
            raise ValueError("Expected at least two residues in the PDB file")

        restraint = mm.CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
        restraint.setUsesPeriodicBoundaryConditions(True)
        restraint.addGlobalParameter('r0', r0 * unit.nanometer)
        restraint.addGlobalParameter('k',k)
        
        groups = [[atom.index for atom in residue.atoms()] for residue in residues]

        group_indices = []
        for group in groups:
            group_idx = restraint.addGroup(group)
            group_indices.append(group_idx)

        # Add a bond between each pair of groups
        for g1, g2 in combinations(group_indices, 2):
            restraint.addBond([g1, g2], [])

        self.system.addForce(restraint)
        
        if position_restraint:
            assert isinstance(position_restraint, (list, tuple)), (
            f"'position_restraint' must be a list or tuple of 1-based molecule indices, "
            f"got {type(position_restraint).__name__}."
            )
            self.ostream.print_info(f"Setting position restraint on molecule(s): {', '.join(map(str, position_restraint))}")
            self.ostream.flush()

            pos_restraint = mm.CustomExternalForce('k_restraint*periodicdistance(x, y, z, x0, y0, z0)^2')
            self.system.addForce(pos_restraint)
            pos_restraint.addGlobalParameter('k_restraint', 10000 * unit.kilojoules_per_mole/unit.nanometer)
            pos_restraint.addPerParticleParameter('x0')
            pos_restraint.addPerParticleParameter('y0')
            pos_restraint.addPerParticleParameter('z0')
            
            for mol_idx in position_restraint:
                for atom_idx in groups[mol_idx - 1]:
                    pos_restraint.addParticle(atom_idx, self.positions[atom_idx])
                    
        

