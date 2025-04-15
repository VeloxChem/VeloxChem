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
import numpy as np
import sys
import time

from .veloxchemlib import mpi_master, bohr_in_angstrom
from .mmforcefieldgenerator import MMForceFieldGenerator
from .molecule import Molecule
from .outputstream import OutputStream


class SolvationBuilder:
    """
    Builds systems for molecular dynamics simulations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - threshold: The threshold for overlap between molecules in angstrom.
        - acceleration_threshold: The threshold of number of molecules to increase the batch size to speed up the process.
        - number_of_attempts: The number of attempts to insert a molecule.
        - random_rotation: Boolean flag to indicate if the solvent molecules will be randomly rotated. False for linear molecules.
        - solute: The VeloxChem molecule object of the solute.
        - solute_labels: The list of atom labels of the solute.
        - solute_ff: The ForceField object of the solute.
        - solvents: The list of VeloxChem molecule objects of the solvent molecules.
        - solvent_labels: The list of atom labels of the solvent molecules.
        - quantities: The list of quantities of each solvent molecule.
        - added_solvent_counts: The list of the number of solvent molecules added to the system.
        - solvent_name: The name of the solvent molecule.
        - system_molecule: The VeloxChem molecule object of the solvated system.
        - box: The dimensions of the box (x, y, z).
        - centered_solute: The list of tuples with the molecule id, atom labels and coordinates of the centered solute.
        - system: The list of tuples with the molecule id, atom labels and coordinates of the system.
        - equilibration_flag: Boolean flag to indicate if an has been requested.
        - temperature: The temperature of the for the equilibration in Kelvin.
        - pressure: The pressure for the equilibration in bar.
        - steps: The number of steps for the equilibration.
        - pcharge: The positive charge for the counterion. Default is 'Na'.
        - ncharge: The negative charge for the counterion. Default is 'Cl'.
        - counterion: The VeloxChem molecule object of the counterion.
        - parent_forcefield: The name of the parent forcefield. Default is 'amber03'.
    """

    def __init__(self, comm=None, ostream=None):
        '''
        Initialize the System Builder class.
        '''

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # MPI information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        # Output stream
        self.ostream = ostream

        # Packing configuration
        self.threshold = 1.8
        self.acceleration_threshold = 1000
        self.number_of_attempts = 100
        self.random_rotation = True
        self.failures_factor = 0.7

        # Molecules
        # Solute
        self.solute = None
        self.solute_labels = []
        self.solute_ff = None
        # Solvents
        self.solvents = []
        self.solvent_labels = []
        self.quantities = []
        self.added_solvent_counts = []
        self.solvent_name = None

        # System
        self.system_molecule = None
        self.box = None
        self.centered_solute = None
        # This is a list of tuples with the molecule id, atom labels and coordinates of the system
        self.system = []

        # NPT Equilibration options
        self.equilibration_flag = False
        self.temperature = 300
        self.pressure = 1
        self.steps = 5000

        # Neutralization options
        self.pcharge = 'Na'
        self.ncharge = 'Cl'
        self.counterion = None

        # Standard forcefield
        self.parent_forcefield = 'amber03'

    def solvate(self,
                solute,
                solvent='spce',
                solvent_molecule=None,
                padding=1,
                target_density=None,
                neutralize=True,
                equilibrate=False):
        """
        Create a solvated system with the most typical solvent molecules.

        :param solute:
            The VeloxChem molecule object of the solute.
        :param solvent:
            The name of the solvent molecule. The default is 'water'.
            Available options: 'spce', 'tip3p', 'ethanol', 'methanol', 'acetone', 
            'chloroform', 'hexane', 'toluene', 'dcm', 'benzene', 'dmso', 'thf', 
            'acetonitrile', 'other' or 'itself'.
                * 'other': The solvent molecule must be provided.
                * 'itself': The solute molecule is used as the solvent as well.
        :param solvent_molecule:
            The VeloxChem molecule object of the solvent molecule. Required if solvent is 'other'.
        :param padding:
            The padding to be added to the bounding box of the molecule in nm.
        :param target_density:
            The target density of the solvent in kg/m^3. If None, the experimental density is used.
            If the solvent is 'other' the target density must be provided or will be estimated from the solvent molecule.
        :param neutralize:
            If True, neutralizes the total charge of the molecule with a counterion.
        :param equilibrate:
            If True, perform an equilibration of the system.
        """

        from scipy.spatial import cKDTree

        # Save the solvent name
        self.solvent_name = solvent

        header_msg = "VeloxChem System Builder"
        self.ostream.print_header(header_msg)
        self.ostream.print_header("=" * len(header_msg))
        self.ostream.print_blank()
        self.ostream.flush()

        # Print solvent information
        self.ostream.print_info(
            f"Solvating the solute with {solvent} molecules")
        self.ostream.print_info(f"Padding: {padding} nm")
        self.ostream.print_blank()
        self.ostream.flush()
        if equilibrate:
            self.ostream.print_info("NPT Equilibration of the box requested")
            self.ostream.print_blank()
            self.ostream.flush()

        # Add the solute to the system
        self._load_solute_molecule(solute)

        # Determine the size of the box
        box_size = self._determine_cubic_box_size(
            solute.get_coordinates_in_angstrom(), padding)

        # Determine the solute volume and convert it to nm^3
        solute_volume = self._get_volume(solute) * 1e-3
        self.ostream.print_info(
            "The volume of the solute is: {:.2f} nm^3".format(solute_volume))
        self.ostream.flush()

        # Define the box
        self._define_box(box_size, box_size, box_size)
        self.ostream.print_info(
            "The box size is: {:.2f} x {:.2f} x {:.2f} nm^3".format(
                box_size * 0.1, box_size * 0.1, box_size * 0.1))
        self.ostream.flush()

        # Accesible volume for the solvent
        volume_nm3 = box_size**3 * 0.001 - solute_volume
        self.ostream.print_info(
            "The volume available for the solvent is: {:.2f} nm^3".format(
                volume_nm3))
        self.ostream.flush()

        if solvent == 'other':

            if solvent_molecule is None:
                raise ValueError(
                    f"The solvent molecule must be provided if the solvent is 'other'"
                )
            if target_density is None:
                raise ValueError(
                    f"The target density must be provided if the solvent is 'other'"
                )
            self.ostream.print_info(
                f"The target density of the solvent is: {target_density} kg/m^3"
            )
            self.ostream.flush()
            # Register the solvent molecule and its quantity to be added to the system
            mols_per_nm3 = self._density_to_mols_per_nm3(
                solvent_molecule, target_density)

            # Calculate the number of solvent molecules to be added (rounded to the nearest integer)
            number_of_solvents = int(mols_per_nm3 * volume_nm3)
            self.ostream.print_info(
                f"The number of solvent molecules to be added to match the target density is: {number_of_solvents}"
            )
            self.ostream.flush()

            # Register the solvent molecule and its quantity to be added to the system
            self._load_solvent_molecule(solvent_molecule, number_of_solvents)

        elif solvent == 'itself':

            if target_density is None:
                raise ValueError(
                    f"The target density must be provided if the solvent is 'itself'"
                )
            self.ostream.print_info(
                f"The target density of the solvent is: {target_density} kg/m^3"
            )
            self.ostream.flush()

            solvent_molecule = solute
            # Extract the properties of the solute
            mols_per_nm3 = self._density_to_mols_per_nm3(
                solvent_molecule, target_density)

            # Calculate the number of solvent molecules to be added (rounded to the nearest integer)
            number_of_solvents = int(mols_per_nm3 * volume_nm3)
            self.ostream.print_info(
                f"The number of solvent molecules to be added to match the target density is: {number_of_solvents}"
            )
            self.ostream.flush()

            # Register the solvent molecule and its quantity to be added to the system
            self._load_solvent_molecule(solute, number_of_solvents)

        else:
            # Extract the properties of the solvent
            mols_per_nm3, density, smiles_code = self._solvent_properties(
                solvent)

            # Calculate the number of solvent molecules to be added (rounded to the nearest integer)
            number_of_solvents = int(mols_per_nm3 * volume_nm3)
            self.ostream.print_info(
                f"The experimental density of {solvent} is: {density} kg/m^3")
            self.ostream.flush()

            # Register the solvent molecule and its quantity to be added to the system
            solvent_molecule = Molecule.read_smiles(smiles_code)
            self._load_solvent_molecule(solvent_molecule, number_of_solvents)

        # Process the solute
        solute_xyz = self.solute.get_coordinates_in_angstrom()
        centroid = np.mean(solute_xyz, axis=0)

        # Create a box with the origin at the centroid
        box_center = box_size / 2

        # Define the molecule id
        molecule_id = 0

        # Translate the solute to the center of the box
        translation = box_center - centroid
        self.centered_solute = [
            (molecule_id, label, coord + translation)
            for label, coord in zip(self.solute_labels, solute_xyz)
        ]

        # Add the centered solute to the system
        self.system.extend(self.centered_solute)

        # Initial build of the tree
        existing_coords = np.array([atom[-1] for atom in self.system])
        tree = cKDTree(existing_coords) if existing_coords.size > 0 else None

        if neutralize:

            charge = int(self.solute.get_charge())
            self.added_counterions = abs(charge)

            if charge == 0:
                self.ostream.print_info(
                    "The solute is neutral, no counterions will be added")
                self.ostream.flush()

            elif charge > 0:
                self.ostream.print_info(
                    f"The solute has a charge of {charge}, adding {abs(charge)}{self.ncharge} counterions"
                )
                self.ostream.flush()

                self.ion_name = self.ncharge
                self.counterion = self._counterion_molecules()
                # Update the required number of solvents
                number_of_solvents -= abs(charge)

            elif charge < 0:
                self.ostream.print_info(
                    f"The solute has a charge of {charge}, adding {abs(charge)} {self.pcharge} counterions"
                )
                self.ostream.flush()

                self.ion_name = self.pcharge
                self.counterion = self._counterion_molecules()
                number_of_solvents -= abs(charge)

        # Insert the counterions
        if self.counterion:
            for _ in range(abs(charge)):
                result = self._insert_molecule(self.counterion, tree)
                if result:
                    self.system.extend(result)
                    existing_coords = np.array(
                        [atom[-1] for atom in self.system])
                    # Update the KDTree with the counterion
                    tree = cKDTree(existing_coords)

        # Solvate the solute with the solvent molecules
        # This dynamic batch size is used to avoid building too often the KDTree.

        max_batch_size = int(number_of_solvents / 2)
        min_batch_size = 10

        # Check linearity of the solvent molecule
        if self.solvent_name == 'itself':
            solvent_molecule = solute
        else:
            solvent_molecule = solvent_molecule

        if self._molecule_linearity(solvent_molecule):
            self.ostream.print_info(
                f"The solvent molecule is linear, random rotation is disabled")
            self.ostream.flush()
            self.random_rotation = False

        # Solvate the solute with the solvent molecules
        start = time.time()
        for i, (solvent,
                quantity) in enumerate(zip(self.solvents, self.quantities)):

            added_count = 0
            # If the quantity is small, set the batch size to 1
            if quantity < self.acceleration_threshold:
                batch_size = 1
            else:
                # Speed up the process by increasing the batch size
                batch_size = self._compute_batch_size(
                    added_count=added_count,
                    max_batch_size=max_batch_size,
                    min_batch_size=min_batch_size,
                    total_quantity=quantity)

            failure_count = 0
            max_failures = self.failures_factor * quantity

            while added_count < quantity:
                attempts = min(batch_size, quantity - added_count)
                new_molecules = []
                for _ in range(attempts):
                    result = self._insert_molecule(solvent, tree)
                    if result:
                        new_molecules.extend(result)
                        added_count += 1
                        failure_count = 0
                    else:
                        failure_count += 1

                if failure_count >= max_failures:
                    self.ostream.print_info(
                        f"Failed to pack {quantity - added_count} out of {quantity} molecules after {failure_count} attempts"
                    )
                    self.ostream.flush()
                    break

                if new_molecules:
                    self.system.extend(new_molecules)
                    new_coords = np.array([atom[-1] for atom in new_molecules])
                    existing_coords = np.vstack((existing_coords, new_coords))
                    tree = cKDTree(existing_coords)

            self.added_solvent_counts.append(added_count)
            msg = f"Solvated system with {added_count} solvent molecules out of {quantity} requested"
            self.ostream.print_info(msg)
            self.ostream.flush()

        end = time.time()
        self.ostream.print_info(
            f"Time to solvate the system: {end - start:.2f} s")
        self.ostream.flush()
        # Print results
        self.ostream.print_info(
            f'The density of the solvent after packing is:{self._check_density(solvent_molecule, self.added_solvent_counts[0], volume_nm3)} kg/m^3'
        )
        self.ostream.flush()

        self.system_molecule = self._save_molecule()

        if equilibrate:
            self.equilibration_flag = True
            self.ostream.print_blank()
            self.ostream.print_info("Equilibrating the system")
            self.ostream.print_blank()
            self.ostream.flush()
            self.ostream.print_info(f"Duration: {self.steps/1000} ps")
            self.ostream.flush()
            self.ostream.print_info(f"Temperature: {self.temperature} K")
            self.ostream.flush()
            self.ostream.print_info(f"Pressure: {self.pressure} bar")
            self.ostream.flush()
            self.ostream.print_blank()
            start = time.time()
            self.perform_equilibration()
            self.ostream.print_info("Equilibration completed, system saved")
            self.ostream.flush()
            end = time.time()
            self.ostream.print_info(
                f"Elapsed time to equilibrate the system: {end - start:.2f} s")
            self.ostream.flush()

    def custom_solvate(self, solute, solvents, quantities, box):
        '''
        Solvate the solute

        :param solute:
            The VeloxChem molecule object of the solute
        :param solvents:
            The list of VeloxChem molecule objects of the solvent molecules
        :param quantities:
            The list of quantities of each solvent molecule. The order must match the order of the solvent molecules.
        :param box:
            The array with the dimensions of the box (x, y, z)
        '''

        from scipy.spatial import cKDTree

        # Add the solute to the system
        self._load_solute_molecule(solute)

        # Define the box
        self._define_box(*box)

        # Initialize solvent counts
        self.added_solvent_counts = [0] * len(solvents)

        # Load the solvent molecules
        for solvent, quantity in zip(solvents, quantities):
            self._load_solvent_molecule(solvent, quantity)

        # Extract the coordinates of the solute
        solute_xyz = self.solute.get_coordinates_in_angstrom()

        # Calculate the centroid of the solute
        centroid = np.mean(solute_xyz, axis=0)

        # Create a box with the origin at the centroid
        box_size = np.array(self.box)
        box_center = box_size / 2

        # Translate the solute to the center of the box
        molecule_id = 0
        translation = box_center - centroid
        # Define the centered solute with the molecule id
        self.centered_solute = [
            (molecule_id, label, coord + translation)
            for label, coord in zip(self.solute_labels, solute_xyz)
        ]

        # Add the centered solute to the system
        self.system.extend(self.centered_solute)

        # Initial build of the tree
        existing_coords = np.array([atom[-1] for atom in self.system])
        tree = cKDTree(existing_coords) if existing_coords.size > 0 else None

        # Solvate the solute with the solvent molecules
        for i, (solvent,
                quantity) in enumerate(zip(self.solvents, self.quantities)):
            added_count = 0
            while added_count < quantity:
                result = self._insert_molecule(solvent, tree)
                if result:
                    self.system.extend(result)
                    existing_coords = np.array(
                        [atom[-1] for atom in self.system])
                    # Update the KDTree with the new solvent
                    tree = cKDTree(existing_coords)
                    added_count += 1
            self.added_solvent_counts[i] = added_count
            msg = f"Solvated system with {added_count} solvent molecules out of {quantity} requested"
            self.ostream.print_info(msg)
            self.ostream.flush()

        self.system_molecule = self._save_molecule()

    def write_gromacs_files(self,
                            solute_ff=None,
                            solvent_ffs=None,
                            equilibration=False):
        '''
        Generates the ForceField for the system

        :param solute_ff:
            The ForceField object of the solute
        :param solvent_ffs:
            The list of ForceField objects of the solvent molecules
        :param equilibration:
            Boolean flag to indicate if the gromacs files will be used for equilibration.
            If True, printouts will not be displayed.
        '''

        self._generate_forcefields(solute_ff, solvent_ffs, equilibration)

        # Special case for 'itself' solvent
        if self.solvent_name == 'itself':
            # Write the itp and top files
            self.solute_ff.write_itp('liquid.itp', 'MOL')
            self.solute_ff.write_top('liquid.top', 'liquid.itp', 'MOL')
            # Change the number of molecules in solute.top to the total number of molecules
            with open('liquid.top', 'r') as f:
                lines = f.readlines()
            with open('liquid.top', 'w') as f:
                for line in lines:
                    if line.startswith('; Compound '):
                        f.write(line)
                        f.write(
                            f'MOL               {self.added_solvent_counts[0] + 1}\n'
                        )
                        break
                    else:
                        f.write(line)

            if not equilibration:
                self.ostream.print_info(
                    "liquid.itp, liquid.top, and solute.top files written")
                # Write the system GRO file
                self.ostream.flush()

            self._write_system_gro(filename='liquid.gro')

        else:
            # Write the itp files
            self.solute_ff.write_itp('solute.itp', 'MOL')

            if not equilibration:
                self.ostream.print_info("solute.itp file written")
                self.ostream.flush()

            # Only write the solvent itp files if they are not SPCE or TIP3P
            if self.solvent_ffs:
                for i, solvent_ff in enumerate(self.solvent_ffs):
                    solvent_ff.write_itp(f'solvent_{i+1}.itp', f'SOL{i+1}')
                    if not equilibration:
                        self.ostream.print_info(
                            f"solvent_{i+1}.itp file written")
                        self.ostream.flush()

            # Post-treatment of the files to ensure compatibility with GROMACS parsing.
            # Extract the atom types from the itp files
            atomtypes = []
            atomtypes.extend(self._extract_atomtypes('solute.itp'))
            if self.solvent_ffs:
                for i in range(len(self.solvent_ffs)):
                    atomtypes.extend(
                        self._extract_atomtypes(f'solvent_{i+1}.itp'))

            # Remove duplicated atom types
            atomtypes = list(set(atomtypes))

            # Remove the atomtypes section from the itp files
            self._remove_atomtypes_section('solute.itp')
            if self.solvent_ffs:
                for i in range(len(self.solvent_ffs)):
                    self._remove_atomtypes_section(f'solvent_{i+1}.itp')

            # If standard forcefields are used, include the parent forcefield and the water model
            directives_list = []
            standard_ff_flag = False

            if self.solvent_name in ['spce', 'tip3p'] or self.counterion:
                standard_ff_flag = True
                parent_directive = f'#include "{self.parent_forcefield}.ff/forcefield.itp"'
                directives_list.append(parent_directive)
                if self.solvent_name in ['spce', 'tip3p']:
                    solvent_directive = f'#include "{self.parent_forcefield}.ff/{self.solvent_name}.itp"'
                    directives_list.append(solvent_directive)
                if self.counterion:
                    counterion_directive = f'#include "{self.parent_forcefield}.ff/ions.itp"'
                    directives_list.append(counterion_directive)

            if standard_ff_flag:
                with open('system.top', 'w') as f:
                    f.write(parent_directive + '\n')
                    f.write('[ atomtypes ]\n')
                    for atomtype in atomtypes:
                        f.write(atomtype + '\n')
                    f.write('\n')
                    f.write(';Residue topologies\n')
                    f.write('#include "solute.itp"\n')
                    # Water model
                    if self.solvent_name in ['spce', 'tip3p']:
                        f.write(solvent_directive + '\n')
                    # Case for counterion with non-standard solvent
                    if self.solvent_ffs:
                        for i in range(len(self.solvent_ffs)):
                            f.write(f'#include "solvent_{i+1}.itp"\n')
                    # Counterion
                    if self.counterion:
                        f.write(counterion_directive + '\n')
                    f.write('\n[ system ]\n')
                    f.write('System\n\n')
                    f.write('[ molecules ]\n')
                    f.write('MOL 1\n')
                    for count in self.added_solvent_counts:
                        f.write(f'SOL {count}\n')
                    if self.counterion:
                        residue_name = self.ion_name.upper()
                        f.write(
                            f'{residue_name} {abs(self.added_counterions)}\n')
                if not equilibration:
                    self.ostream.print_info("system.top file written")
                    self.ostream.flush()

            else:
                # Write the top file based on the forcefields generated
                with open('system.top', 'w') as f:
                    f.write('[ defaults ]\n')
                    f.write(
                        '; nbfunc        comb-rule       gen-pairs        fudgeLJ   fudgeQQ\n'
                    )
                    f.write(
                        '1               2               yes             0.500000  0.833333\n\n'
                    )
                    f.write('[ atomtypes ]\n')
                    for atomtype in atomtypes:
                        f.write(atomtype + '\n')
                    f.write('\n')
                    f.write(';Residue topologies\n')
                    f.write('#include "solute.itp"\n')
                    for i in range(len(self.solvent_ffs)):
                        f.write(f'#include "solvent_{i+1}.itp"\n')
                    if self.counterion:
                        f.write('#include "counterion.itp"\n')
                    f.write('\n[ system ]\n')
                    f.write('System\n\n')
                    f.write('[ molecules ]\n')
                    f.write('MOL 1\n')
                    for i, count in enumerate(self.added_solvent_counts):
                        f.write(f'SOL{i+1} {count}\n')
                    if self.counterion:
                        f.write('ION 1\n')
                if not equilibration:
                    self.ostream.print_info("system.top file written")
                    self.ostream.flush()

            # Write the system GRO file
            self._write_system_gro()

            if not equilibration:
                self.ostream.print_info("system.gro file written")
                self.ostream.flush()

    def write_openmm_files(self, solute_ff=None, solvent_ffs=None):
        '''
        Generates the PDB and XML files for OpenMM.

        :param solute_ff:
            The ForceField object of the solute (optional).
        :param solvent_ffs:
            The list of ForceField objects of the solvent molecules (optional).
        '''

        # Generate the forcefields
        self._generate_forcefields(solute_ff, solvent_ffs)

        # Special case for 'itself' solvent
        if self.solvent_name == 'itself':
            # Write the XML and PDB files for the solute
            self.solute_ff.write_openmm_files('liquid', 'MOL')
            self.ostream.print_info("liquid.xml file written")
            self.ostream.flush()
            # Write the system PDB file
            filename = 'liquid.pdb'

        else:
            # Solute
            self.solute_ff.write_openmm_files('solute', 'MOL')
            self.ostream.print_info(
                "system.pdb, solute.pdb, and solute.xml files written")
            self.ostream.flush()

            # Not standard solvents
            if self.solvent_name not in ['spce', 'tip3p']:
                for i, solvent_ff in enumerate(self.solvent_ffs):
                    solvent_ff.write_openmm_files(f'solvent_{i+1}',
                                                  f'S{i+1:02d}')
                    self.ostream.print_info(
                        f"solvent_{i+1}.pdb and solvent_{i+1}.xml files written"
                    )
                    self.ostream.flush()
            else:
                self.ostream.print_info(
                    f"Using standard AMBER {self.solvent_name} forcefield, no solvent xml files will be written"
                )
                self.ostream.print_info(
                    f'Remember to include amber03.xml and the {self.solvent_name}.xml file while creating the OpenMM system'
                )
                self.ostream.flush()

            filename = 'system.pdb'

        # Write the system PDB file
        if self.equilibration_flag:
            # If the system was equilibrated, OpenMM has already written the PDB file
            Path('equilibrated_system.pdb').rename(filename)
        else:
            self._write_system_pdb(filename=filename)
        # Print information
        self.ostream.print_info(f"{filename} file written")
        self.ostream.flush()

    def perform_equilibration(self):
        """
        Performs an equilibration using OpenMM.
        """

        try:
            import openmm as mm
            import openmm.app as app
            import openmm.unit as unit

        except ImportError:
            raise ImportError("OpenMM is required for this functionality")

        # Generate the forcefields using semiempirical charges.

        solute_ff = MMForceFieldGenerator()
        solute_ff.ostream.mute()
        solute_ff.partial_charges = self.solute.get_partial_charges(
            self.solute.get_charge())
        solute_ff.create_topology(self.solute)

        if self.solvent_name in ['spce', 'tip3p', 'itself']:
            solvent_ffs = None

        else:
            solvent_ffs = []
            for i in range(len(self.solvents)):
                solvent_ff = MMForceFieldGenerator()
                solvent_ff.ostream.mute()
                solvent_ff.partial_charges = self.solvents[
                    i].get_partial_charges(self.solvents[i].get_charge())
                solvent_ff.create_topology(self.solvents[i])
                solvent_ffs.append(solvent_ff)

        self.write_gromacs_files(solute_ff, solvent_ffs, equilibration=True)

        # Load the system
        if self.solvent_name == 'itself':
            gro = app.GromacsGroFile('liquid.gro')
        else:
            gro = app.GromacsGroFile('system.gro')

        # # Create the force field
        if self.solvent_name == 'itself':
            forcefield = app.GromacsTopFile(
                'liquid.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
        else:
            forcefield = app.GromacsTopFile(
                'system.top', periodicBoxVectors=gro.getPeriodicBoxVectors())

        topology = forcefield.topology
        positions = gro.positions

        # Create the OpenMM system
        system = forcefield.createSystem(
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometers,
            constraints=app.HBonds,
        )

        # Set the temperature and pressure
        integrator = mm.LangevinIntegrator(self.temperature * unit.kelvin,
                                           1.0 / unit.picosecond,
                                           1 * unit.femtosecond)
        barostat = mm.MonteCarloBarostat(self.pressure * unit.bar,
                                         self.temperature * unit.kelvin, 25)
        system.addForce(barostat)

        # Create the simulation
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)

        # Minimize the energy
        simulation.minimizeEnergy()

        # Equilibrate
        simulation.reporters.append(
            app.StateDataReporter('equilibration.log',
                                  1000,
                                  step=True,
                                  potentialEnergy=True,
                                  temperature=True,
                                  volume=True))
        simulation.step(self.steps)

        # Get the final positions
        positions = simulation.context.getState(
            getPositions=True).getPositions()

        # Exctract the new box size from the simulation context
        box_vectors = simulation.context.getState().getPeriodicBoxVectors()
        # Update the box size in angstrom
        self.box = [
            box_vectors[0][0].value_in_unit(unit.angstroms),
            box_vectors[1][1].value_in_unit(unit.angstroms),
            box_vectors[2][2].value_in_unit(unit.angstroms)
        ]
        self.ostream.print_info(
            f'The box size after equilibration is: {self.box[0] * 0.1:.2f} x {self.box[1] * 0.1:.2f} x {self.box[2] * 0.1:.2f} nm^3'
        )
        self.ostream.flush()
        # Recalculate the available volume for the solvent
        volume_nm3 = (self.box[0] * self.box[1] * self.box[2] -
                      self._get_volume(self.solute)) * 1e-3
        # Recalculate the density of the solvent
        self.ostream.print_info(
            f'The density of the solvent after equilibration is: {self._check_density(self.solvents[0], self.added_solvent_counts[0], volume_nm3)} kg/m^3'
        )
        self.ostream.flush()
        # Write the PDB file
        with open('equilibrated_system.pdb', 'w') as f:
            app.PDBFile.writeFile(simulation.topology, positions, f)

        # Delete the produced gro and top files with Path
        if self.solvent_name == 'itself':
            Path('liquid.gro').unlink()
            Path('liquid.top').unlink()
            Path('liquid.itp').unlink()
        elif self.solvent_name in ['spce', 'tip3p']:
            Path('system.gro').unlink()
            Path('system.top').unlink()
            Path('solute.itp').unlink()
        else:
            Path('system.gro').unlink()
            Path('system.top').unlink()
            Path('solute.itp').unlink()
            for i in range(len(self.solvent_ffs)):
                Path(f'solvent_{i+1}.itp').unlink()

        # Update the system molecule
        self.system_molecule = Molecule.read_pdb_file('equilibrated_system.pdb')

    # Auxiliary functions

    def _determine_cubic_box_size(self, coordinates, padding):
        """
        Determines the size of a cubic box based on the molecular geometry and a padding parameter.
        
        :param coordinates: 
            The array of shape (n, 3) where n is the number of atoms, and each row contains the (x, y, z) coordinates of an atom.
        :param padding: 
            The padding to be added to the bounding box of the molecule in nm.
        :return: 
            The size of the cubic box in the same units as the input coordinates.
        """
        # Calculate the minimum and maximum coordinates along each axis
        min_coords = np.min(coordinates, axis=0)
        max_coords = np.max(coordinates, axis=0)

        # Determine the extent of the bounding box
        extent = max_coords - min_coords

        # Add the padding to the extent
        box_size = np.max(extent) + 2 * padding * 10

        return box_size

    def _load_solute_molecule(self, solute):
        '''
        Register the solute molecule to be added to the system
        '''

        self.solute = solute
        self.solute_labels = solute.get_labels()

    def _define_box(self, dim_x, dim_y, dim_z):
        '''
        Create a box around the solute
        '''

        self.box = [dim_x, dim_y, dim_z]

    def _load_solvent_molecule(self, solvent, quantity):
        '''
        Register the solvent molecule and its quantity to be added to the system

        :param solvent:
            The VeloxChem molecule object of the solvent
        :param quantity:
            The quantity of the solvent molecules to be added
        '''

        self.solvents.append(solvent)
        self.quantities.append(quantity)
        self.solvent_labels.append(solvent.get_labels())

    def _insert_molecule(self, new_molecule, tree):
        """
        Insert a molecule with a random rotation without rebuilding the KDTree each time.
        """

        from scipy.spatial.transform import Rotation

        new_molecule_xyz = new_molecule.get_coordinates_in_angstrom()
        new_molecule_labels = new_molecule.get_labels()

        total_attempts = self.number_of_attempts

        for attempt_num in range(1, total_attempts + 1):

            if self.random_rotation:
                # Generate a random rotation
                rotation_matrix = Rotation.random().as_matrix()
                # Rotate the molecule coordinates
                rotated_coords = new_molecule_xyz @ rotation_matrix.T

            else:
                rotated_coords = new_molecule_xyz

            # Generate a random position for the solvent molecule
            position = np.random.rand(3) * self.box

            # Translate the rotated solvent molecule to the new position
            translated_solvent_coords = rotated_coords + position

            # Check if the solvent molecule is within the box
            if not np.all((translated_solvent_coords >= 0)
                          & (translated_solvent_coords <= self.box)):
                continue

            # Check for overlap using KD-Tree
            if tree:
                # Query the KD-Tree for atoms within the threshold distance
                indices = tree.query_ball_point(translated_solvent_coords,
                                                self.threshold)
                if any(indices):
                    continue

            # No overlap detected; return the translated molecule
            molecule_id = len(self.system)
            translated_solvent = [
                (molecule_id, label, coord) for label, coord in zip(
                    new_molecule_labels, translated_solvent_coords)
            ]
            # If number of attempts approaches the total attempts print a warning
            if attempt_num == int(0.9 * total_attempts):
                self.ostream.print_info(
                    f"Warning: {attempt_num} attempts have been made to insert the solvent molecule"
                )
                self.ostream.print_info("Consider reducing the target density")
                self.ostream.flush()

            return translated_solvent

        return None

    def _check_overlap(self, new_coords, existing_tree):
        '''
        Check for overlap using a KD-Tree.
        '''
        indices = existing_tree.query_ball_point(new_coords, self.threshold)
        # If any indices are returned, overlap exists
        return any(indices)

    def _save_molecule(self):
        '''
        Creates a VeloxChem molecule object from the system
        '''

        # Extract atom labels and coordinates from the system
        labels = [atom[1] for atom in self.system]
        coordinates = np.array([atom[-1] for atom in self.system])

        # Create the VeloxChem molecule object
        molecule = Molecule(symbols=labels,
                            coordinates=coordinates,
                            units='angstrom')

        return molecule

    def _extract_atomtypes(self, itp_filename):
        """
        Extracts atom types from an ITP file, removes duplicates, and returns a list of unique atom types.
        
        :param itp_filename: The name of the ITP file from which to extract atom types.
        :return: A list of unique atom types.
        """
        atomtypes = []
        inside_atomtypes_block = False

        with open(itp_filename, 'r') as file:
            for line in file:
                if line.strip().startswith('[ atomtypes ]'):
                    inside_atomtypes_block = True
                    continue  # Skip the [ atomtypes ] header line

                if inside_atomtypes_block:
                    if line.strip() == '' or line.startswith('['):
                        # End of the atomtypes block
                        inside_atomtypes_block = False
                    elif not line.strip().startswith(';'):
                        # Add the atom type line to the list if it's not a comment
                        atomtypes.append(line.strip())

        # Remove duplicates by converting to a set and back to a list
        unique_atomtypes = list(set(atomtypes))

        return unique_atomtypes

    def _remove_atomtypes_section(self, itp_filename):
        """
        Removes the [ atomtypes ] section from an ITP file.
        
        :param itp_filename: The name of the ITP file from which to remove the atom types section.
        """
        new_lines = []
        inside_atomtypes_block = False

        with open(itp_filename, 'r') as file:
            for line in file:
                if line.strip().startswith('[ atomtypes ]'):
                    inside_atomtypes_block = True
                    continue  # Skip the [ atomtypes ] header line

                if inside_atomtypes_block:
                    if line.strip() == '' or line.startswith('['):
                        # End of the atomtypes block, resume normal line processing
                        inside_atomtypes_block = False
                        new_lines.append(line)
                else:
                    new_lines.append(line)

        # Rewrite the itp file without the atomtypes block
        with open(itp_filename, 'w') as file:
            file.writelines(new_lines)

    def _solvent_properties(self, solvent):
        '''
        Contains the density of the solvent

        :param solvent:
            The name of the solvent
        '''

        if solvent == 'spce':
            mols_per_nm3 = 33.3
            density = 1000
            smiles_code = 'O'

        elif solvent == 'tip3p':
            mols_per_nm3 = 33.3
            density = 1000
            smiles_code = 'O'

        elif solvent == 'ethanol':
            mols_per_nm3 = 10.3
            density = 789
            smiles_code = 'CCO'

        elif solvent == 'methanol':
            mols_per_nm3 = 14.9
            density = 791
            smiles_code = 'CO'

        elif solvent == 'acetone':
            mols_per_nm3 = 8.1
            density = 784
            smiles_code = 'CC(=O)C'

        elif solvent == 'chloroform':
            mols_per_nm3 = 7.6
            density = 1494
            smiles_code = 'C(Cl)(Cl)Cl'

        elif solvent == 'hexane':
            mols_per_nm3 = 4.7
            density = 655
            smiles_code = 'CCCCCC'

        elif solvent == 'toluene':
            mols_per_nm3 = 5.7
            density = 866
            smiles_code = 'CC1=CC=CC=C1'

        elif solvent == 'dcm':
            mols_per_nm3 = 8.3
            density = 1335
            smiles_code = 'ClCCl'

        elif solvent == 'benzene':
            mols_per_nm3 = 6.8
            density = 876
            smiles_code = 'c1ccccc1'

        elif solvent == 'dmso':
            mols_per_nm3 = 8.5
            density = 1101
            smiles_code = 'CS(=O)C'

        elif solvent == 'thf':
            mols_per_nm3 = 7.4
            density = 888
            smiles_code = 'C1CCOC1'

        elif solvent == 'acetonitrile':
            mols_per_nm3 = 11.53
            density = 786
            smiles_code = 'CC#N'

        elif solvent == 'dmf':
            mols_per_nm3 = 12.96
            density = 944
            smiles_code = 'CN(C)C=O'

        else:
            return None

        return mols_per_nm3, density, smiles_code

    def _counterion_molecules(self):
        """
        Generates a molecule object for the counterion.
        """
        coords = np.array([[0.0, 0.0, 0.0]])

        if self.ion_name == 'Na':
            labels = ['Na']
            charge = 1.0
        elif self.ion_name == 'K':
            labels = ['K']
            charge = 1.0
        elif self.ion_name == 'Li':
            labels = ['Li']
            charge = 1.0
        elif self.ion_name == 'Cl':
            labels = ['Cl']
            charge = -1.0

        ion_mol = Molecule(labels, coords, 'angstrom')
        ion_mol.set_charge(charge)

        return ion_mol

    def _check_density(self, molecule, number_of_molecules, volume):
        """
        Given the moecule mass in gr/mol and the volume in cubic nm,
        return the density in kg/m^3.

        :param molecule:
            The VeloxChem molecule object.
        :param number_of_molecules:
            The number of molecules.
        :param volume:
            The volume in cubic nm.
        """
        mass = sum(molecule.get_masses()) * 1e-3
        volume = volume * 1e-27
        moles = number_of_molecules / 6.022e23
        density = mass * moles / volume

        return int(density)

    def _density_to_mols_per_nm3(self, molecule, density):
        """
        Given the density in kg/m^3, return the number of moles per nm^3.

        :param density:
            The density in kg/m^3.

        :return:
            The number of moles per nm^3.
        """
        # Get the mass in kg/mol
        mass = sum(molecule.get_masses()) * 1e-3 / 6.022e23

        # Get the mols per m^3
        mols_per_m3 = density / mass

        # Convert to mols per nm^3
        mols_per_nm3 = mols_per_m3 * 1e-27

        return mols_per_nm3

    def _generate_forcefields(self,
                              solute_ff,
                              solvent_ffs,
                              equilibration=False):
        """
        Generate the force fields for the solute and solvent molecules.
        The forcefields get stored in the solute_ff and solvent_ffs attributes (lists).

        :param solute_ff:
            The ForceField object of the solute.
        :param solvent_ffs:
            The list of ForceField objects of the solvent molecules.
        :param equilibration:
            Boolean flag to indicate if the gromacs files will be used for equilibration.
            If True, printouts will not be displayed.
        """

        # Solute
        if not solute_ff:
            self.solute_ff = MMForceFieldGenerator()
            self.solute_ff.ostream.mute()
            if not equilibration:
                self.ostream.print_info(
                    'Generating the ForceField for the solute')
                self.ostream.flush()
            self.solute_ff.create_topology(self.solute)
            if not equilibration:
                self.ostream.print_info(
                    'Generated the ForceField for the solute')
                self.ostream.flush()
        else:
            self.solute_ff = solute_ff

        # Solvents
        # Special case for itself
        if not solvent_ffs:
            if self.solvent_name == 'itself':
                self.solvent_ffs = []
                self.solvent_ffs.append(self.solute_ff)
            # Non-water solvents
            elif self.solvent_name not in ['spce', 'tip3p']:
                self.solvent_ffs = []
                for solvent in self.solvents:
                    solvent_ff = MMForceFieldGenerator()
                    solvent_ff.ostream.mute()
                    if not equilibration:
                        self.ostream.print_info(
                            f'Generating the ForceField for the solvent')
                        self.ostream.flush()
                    solvent_ff.create_topology(solvent)
                    if not equilibration:
                        self.ostream.print_info(
                            f'Generated the ForceField for the solvent')
                        self.ostream.flush()
                    self.solvent_ffs.append(solvent_ff)
            else:
                self.solvent_ffs = None
        else:
            self.solvent_ffs = solvent_ffs

    def _write_system_gro(self, filename='system.gro'):
        """
        Write the system's molecule into a GRO file.
        """

        coords_in_nm = self.system_molecule.get_coordinates_in_angstrom() * 0.1

        if self.solvent_name == 'itself':
            # All the molecules are the same based on the solute
            # Write the GRO file
            with open(filename, 'w') as f:
                # Header
                f.write('Generated by VeloxChem\n')
                f.write(f'{self.system_molecule.number_of_atoms()}\n')

                # Initialize residue_offset before the molecule loop
                residue_offset = 0
                atom_counter = 1
                # Atoms
                for mols in range(self.added_solvent_counts[0] + 1):
                    residue_number = mols + 1
                    for i, atom in self.solute_ff.atoms.items():
                        atom_name = atom['name']
                        line_str = f'{residue_number:>5d}{"MOL":<5s}{atom_name:<5s}{atom_counter:>5d}'
                        if residue_number > 9999:
                            residue_number -= 9999
                        for d in range(3):
                            line_str += f'{coords_in_nm[residue_offset + i][d]:{8}.{3}f}'
                        line_str += '\n'
                        f.write(line_str)
                        atom_counter += 1
                        # GRO has a maximum of 5 digits for the atom index
                        if atom_counter > 99999:
                            atom_counter -= 99999
                    # Increment residue_offset after each molecule
                    residue_offset += len(self.solute_ff.atoms)

                # Box
                for d in range(3):
                    f.write(f'{self.box[d]*0.1:{10}.{5}f}')

        else:
            with open(filename, 'w') as f:
                # Header
                f.write('Generated by VeloxChem\n')
                f.write(f'{self.system_molecule.number_of_atoms()}\n')

                # Atoms
                # Solute
                # It will be added a counter to keep track of the atom index
                # in order to write the correct coordinates in the GRO file.
                # This first counter is for the size of the solute
                counter_1 = 0
                for i, atom in self.solute_ff.atoms.items():
                    atom_name = atom['name']
                    line_str = f'{1:>5d}{"MOL":<5s}{atom_name:<5s}{i + 1:>5d}'
                    for d in range(3):
                        line_str += f'{coords_in_nm[i][d]:{8}.{3}f}'
                    line_str += '\n'
                    f.write(line_str)
                    counter_1 += 1

                # Solvents
                # This second counter is for the size of the solvents
                counter_2 = counter_1

                # If the solvent is SPCE or TIP3P, the force field is standardized
                final_residx = 0
                if self.solvent_name in ['spce', 'tip3p']:
                    for i in range(self.added_solvent_counts[0]):
                        for j in range(3):
                            atom_name = ['OW', 'HW1', 'HW2'][j]
                            # The resdue index shall be reset to 1 when it reaches 9999
                            if i > 9999:
                                i -= 9999
                            line_str = f'{i + 1:>5d}{"SOL":<5s}{atom_name:<5s}{j + 1:>5d}'
                            for d in range(3):
                                line_str += f'{coords_in_nm[counter_2][d]:{8}.{3}f}'
                            line_str += '\n'
                            f.write(line_str)
                            counter_2 += 1
                            if counter_2 > 99999:
                                counter_2 -= 99999
                        final_residx = i

                # If the solvent is not SPCE or TIP3P, the information is coming from the force field
                else:
                    for i, solvent in enumerate(self.solvent_ffs):
                        for k in range(self.added_solvent_counts[i]):
                            for j, atom in solvent.atoms.items():
                                atom_name = atom['name']
                                if k > 9999:
                                    k -= 9999
                                line_str = f'{k:>5d}{f"SOL{i+1}":<5s}{atom_name:<5s}{j + 1:>5d}'
                                for d in range(3):
                                    line_str += f'{coords_in_nm[counter_2][d]:{8}.{3}f}'
                                line_str += '\n'
                                f.write(line_str)
                                counter_2 += 1
                                if counter_2 > 99999:
                                    counter_2 -= 99999
                            final_residx = k

                # Counterions
                if self.counterion:
                    # The counterion force field is standardized
                    # In GROMACS the ion residues are named as the ion name in capital letters
                    for i in range(self.added_counterions):
                        atom_name = self.ion_name.upper()
                        line_str = f'{final_residx + i + 2:>5d}{atom_name:<5s}{atom_name:<5s}{i + 1:>5d}'
                        for d in range(3):
                            line_str += f'{coords_in_nm[counter_2][d]:{8}.{3}f}'
                        line_str += '\n'
                        f.write(line_str)
                        counter_2 += 1

                # Box
                for d in range(3):
                    f.write(f'{self.box[d]*0.1:{10}.{5}f}')

    def _write_system_pdb(self, filename='system.pdb'):
        """
        Write the system's molecule into a PDB file.
        """

        # Write the system PDB file
        with open(filename, 'w') as f:
            f.write("HEADER    Generated by VeloxChem\n")
            f.write("TITLE     Solvated system\n")
            # Write the box size to the PDB file
            # Box size is in nm, convert to Angstrom
            f.write(
                f"CRYST1{self.box[0]:9.3f}{self.box[1]:9.3f}{self.box[2]:9.3f}  90.00  90.00  90.00 P 1           1\n"
            )

            atom_counter = 1
            residue_counter = 1
            chain_ids = ['A', 'B', 'C']
            pdb_atom_numbers = {}
            coordinates = self.system_molecule.get_coordinates_in_angstrom()

            # Solute
            residue_name = 'MOL'  # Residue name for the solute
            for i, atom in self.solute_ff.atoms.items():
                atom_name = atom['name']
                element = self.solute.get_labels()[i]
                x, y, z = coordinates[i]
                # PDB format string adhering to column specifications
                f.write(
                    "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                    .format('HETATM', atom_counter, atom_name, '', residue_name,
                            chain_ids[0], residue_counter, '', x, y, z, 1.00,
                            0.00, element))
                pdb_atom_numbers[('solute', i)] = atom_counter
                atom_counter += 1
            residue_counter += 1

            # Initialize coordinate counter after solute atoms
            coordinate_counter = len(self.solute_ff.atoms)

            # Special case for 'itself' solvent
            if self.solvent_name == 'itself':
                residue_name = 'MOL'
                num_atoms_per_molecule = len(self.solute_ff.atoms)
                for mols in range(self.added_solvent_counts[0]):
                    if residue_counter > 9999:
                        residue_counter -= 9999
                    for i, atom in self.solute_ff.atoms.items():
                        atom_name = atom['name']
                        element = self.solute.get_labels()[i]
                        x, y, z = coordinates[coordinate_counter]
                        f.write(
                            "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   "
                            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                            .format('HETATM', atom_counter, atom_name, '',
                                    residue_name, chain_ids[1], residue_counter,
                                    '', x, y, z, 1.00, 0.00, element))
                        pdb_atom_numbers[('solvent', mols, i)] = atom_counter
                        atom_counter += 1
                        coordinate_counter += 1
                    residue_counter += 1

            # Solvents
            elif self.solvent_name in ['spce', 'tip3p']:
                # The force field for SPCE and TIP3P water models are standardized
                # and do not require a separate force field object.
                elements = self.solvents[0].get_labels()
                residue_name = 'HOH'
                if self.added_solvent_counts[0] * len(elements) > 99999:
                    raise ValueError(
                        "The number of solvent atoms exceeds 99999. The PDB format does not support more than 99999 atoms. Write GROMACS files instead."
                    )
                # Atom names are O, H1, and H2
                for i in range(self.added_solvent_counts[0]):
                    if residue_counter > 9999:
                        residue_counter -= 9999
                    for j in range(3):
                        atom_name = ['O', 'H1', 'H2'][j]
                        element = elements[j]
                        x, y, z = coordinates[coordinate_counter]
                        f.write(
                            "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                            .format('ATOM', atom_counter, atom_name, '',
                                    residue_name, chain_ids[1], residue_counter,
                                    '', x, y, z, 1.00, 0.00, element))
                        pdb_atom_numbers[('solvent', 0, i, j)] = atom_counter
                        atom_counter += 1
                        coordinate_counter += 1
                    residue_counter += 1

            else:
                # The force field for other solvents are not standardized
                # and require a separate force field object.
                for i, solvent_ff in enumerate(self.solvent_ffs):
                    elements = self.solvents[i].get_labels()
                    if self.added_solvent_counts[i] * len(elements) > 99999:
                        raise ValueError(
                            "The number of solvent atoms exceeds 99999. The PDB format does not support more than 99999 atoms. Write GROMACS files instead."
                        )
                    for j in range(self.added_solvent_counts[i]):
                        if residue_counter > 9999:
                            residue_counter -= 9999
                        for k in range(len(elements)):
                            residue_name = f'S{i+1:02d}'
                            atom_name = solvent_ff.atoms[k]['name']
                            element = elements[k]
                            x, y, z = coordinates[coordinate_counter]
                            f.write(
                                "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                                .format('HETATM', atom_counter, atom_name, '',
                                        residue_name, chain_ids[1],
                                        residue_counter, '', x, y, z, 1.00,
                                        0.00, element))
                            pdb_atom_numbers[('solvent', i, j,
                                              k)] = atom_counter
                            atom_counter += 1
                            coordinate_counter += 1
                        residue_counter += 1

            # Counterions
            if self.counterion:
                for i in range(self.added_counterions):
                    if residue_counter > 9999:
                        residue_counter -= 9999
                    if self.ion_name in ['Na', 'K', 'Li']:
                        atom_name = self.ion_name + '+'
                        residue_name = self.ion_name + '+'
                    elif self.ion_name == 'Cl':
                        atom_name = 'Cl-'
                        residue_name = 'Cl-'
                    element = self.counterion.get_labels()[0]
                    x, y, z = coordinates[coordinate_counter]
                    f.write(
                        "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                        .format('ATOM', atom_counter, atom_name, '',
                                residue_name, chain_ids[2], residue_counter, '',
                                x, y, z, 1.00, 0.00, element))
                    pdb_atom_numbers[('counterion', i)] = atom_counter
                    atom_counter += 1
                    coordinate_counter += 1
                residue_counter += 1

            # Write CONECT records for bonds
            # Solute bonds
            for (i_atom, j_atom) in self.solute_ff.bonds:
                pdb_i = pdb_atom_numbers[('solute', i_atom)]
                pdb_j = pdb_atom_numbers[('solute', j_atom)]
                f.write(f"CONECT{pdb_i:>5d}{pdb_j:>5d}\n")

            # Solvent bonds
            # If the solvent is SPCE or TIP3P, the bonds are predefined as (0, 1) and (0, 2)
            # and do not require CONECT records.

            if self.solvent_ffs:
                # Solvent bonds for 'itself' solvent
                if self.solvent_name == 'itself':
                    num_atoms_per_molecule = len(self.solute_ff.atoms)
                    for mols in range(self.added_solvent_counts[0]):
                        for (i_atom, j_atom) in self.solute_ff.bonds:
                            pdb_i = pdb_atom_numbers[('solvent', mols, i_atom)]
                            pdb_j = pdb_atom_numbers[('solvent', mols, j_atom)]
                            f.write(f"CONECT{pdb_i:>5d}{pdb_j:>5d}\n")

                else:
                    # Use the bonds defined in the solvent force field
                    for i, solvent_ff in enumerate(self.solvent_ffs):
                        for j in range(self.added_solvent_counts[i]):
                            for (i_atom, j_atom) in solvent_ff.bonds:
                                pdb_i = pdb_atom_numbers[('solvent', i, j,
                                                          i_atom)]
                                pdb_j = pdb_atom_numbers[('solvent', i, j,
                                                          j_atom)]
                                f.write(f"CONECT{pdb_i:>5d}{pdb_j:>5d}\n")

            # Write END line
            f.write("END\n")

    def _get_volume(self, molecule):
        """
        Determine the volume of the molecule based on the vdW radii, with careful intersection correction.
        
        :return:
            The volume of the molecule in cubic Angstrom.
        """
        # Get the atomic van der Waals radii in Angstrom
        atomic_radii = molecule.vdw_radii_to_numpy() * bohr_in_angstrom()

        # Get the coordinates
        coords = molecule.get_coordinates_in_angstrom()

        # Get the number of atoms
        natoms = molecule.number_of_atoms()

        # Initialize the volume to be the sum of individual atomic volumes
        volume = 0.0

        # The initial volume is the sum of the volumes of the vdW spheres
        for i in range(natoms):
            volume += (4.0 / 3.0) * np.pi * atomic_radii[i]**3

        # Correct the volume by subtracting the intersection volumes
        for i in range(natoms):
            for j in range(i + 1, natoms):
                rij = np.linalg.norm(coords[i] - coords[j])
                if rij < (atomic_radii[i] + atomic_radii[j]):
                    # Calculate the volume of the intersection of the spheres
                    r1 = atomic_radii[i]
                    r2 = atomic_radii[j]
                    d = rij

                    # Calculate the intersection volume using a more precise formula
                    # Formula from: https://mathworld.wolfram.com/Sphere-SphereIntersection.html
                    if d > 0:
                        h1 = (r1 - r2 + d) * (r1 + r2 - d) / (2 * d)
                        h2 = (r2 - r1 + d) * (r2 + r1 - d) / (2 * d)

                        intersection_volume = (np.pi * h1**2 *
                                               (3 * r1 - h1) + np.pi * h2**2 *
                                               (3 * r2 - h2)) / 6
                        volume -= intersection_volume

        # Round the volume to 2 decimal places
        volume = round(volume, 2)

        return volume

    def _molecule_linearity(self, molecule):
        """
        Determines if a molecule is linear based on its moment of inertia tensor.

        Parameters:
        - molecule: An object with methods `get_coordinates_in_angstrom()` and
        `number_of_atoms()`. Each atom should have a mass accessible.

        Returns:
        - bool: True if the molecule is linear, False otherwise.
        """

        # Get the coordinates and masses
        coords = molecule.get_coordinates_in_angstrom()
        masses = molecule.get_masses()

        # Calculate the center of mass (COM)
        com = np.average(coords, axis=0, weights=masses)

        # Calculate the moment of inertia tensor as:
        # I = sum(m_i * (ri^2 * 1 - r_i x r_i))
        inertia_tensor = np.zeros((3, 3))
        for i in range(molecule.number_of_atoms()):
            r = coords[i] - com
            inertia_tensor += masses[i] * np.outer(r, r)

        # Diagonalize the inertia tensor
        eigvals, _ = np.linalg.eigh(inertia_tensor)

        # Sort eigenvalues in ascending order
        eigvals = np.sort(eigvals)

        # In a linear molecule one of the eigenvalues is much larger than the other two
        l_1 = eigvals[0]
        l_2 = eigvals[1]
        l_3 = eigvals[2]

        # Linearity coefficient is the ratio of the sum of the two smallest eigenvalues
        # to the largest eigenvalue
        linearity_coefficient = 0.5 * (l_1 + l_2) / l_3

        if linearity_coefficient < 0.05:
            is_linear = True
        else:
            is_linear = False

        return is_linear

    def _compute_batch_size(self, added_count, max_batch_size, min_batch_size,
                            total_quantity):
        """
        Compute the batch size dynamically based on the number of molecules added.
        
        :param added_count: Number of molecules already added to the system.
        :param max_batch_size: Maximum batch size to start with.
        :param min_batch_size: Minimum batch size to avoid being too small.
        :param total_quantity: Total number of molecules to be added.
        :return: Adjusted batch size.
        """
        fraction_completed = added_count / total_quantity
        batch_size = int(max_batch_size * (1 - fraction_completed))
        batch_size = max(batch_size, min_batch_size)
        return batch_size
