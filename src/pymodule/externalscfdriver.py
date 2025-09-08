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


### This class is designed to plug in an external qm_driver which is providing the energy, gradient, hessian, non-adiabatic couplings for excited states ###

from mpi4py import MPI
import numpy as np
from pathlib import Path
from sys import stdout
from time import time
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import re
import shutil
from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .molecule import Molecule
from .clustermanager import ClusterManager
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


class ExternalScfDriver:

    def __init__(self, program, charge, multiplicity, nprocs=8, comm=None, ostream=None, hostname=False, path_on_cluster=None):
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
        self.nprocs = nprocs
        # output stream
        self.ostream = ostream
        self.charge = charge
        self.spin = multiplicity

        self.xc_func = 'b3lyp'
        self.dispersion = 'D3'
        self.basis_set_label = 'def2-svp'
        self.method = 'RKS'

        self.tracked_roots = None
        self.energies = None
        self.spin_flip = False
        self.NAC = False
        self.cluster = False
        self.solvation = (False, 'CPCM', 'water')
        self.add_ghost_atom = False
        self.path_on_cluster = path_on_cluster
        
        self.xyz_filename = None
        self.input_files = []
        self.output_files = []
        self.xyz_structures = 'xyz_structures'
        self.input_filename = 'current_input.inp'
        self.output_filename = 'current_output.out'
        self.input_files.append(self.input_filename)
        self.output_files.append(self.output_filename)
        
        self.jobscript = 'jobscript'
        self.program = program
        self.cluster_manager = None
        self.job_ids = None
        if hostname == True:
            self.cluster_manager = ClusterManager('lyra', self.path_on_cluster)
            self.job_ids = []
    
    def prepare_job_script(self, template_path, input_file, output_file):
        with open(template_path, 'r') as template_file:
            script_content = template_file.read()

        # Replace placeholders with actual file names
        script_content = script_content.replace('INPUT', input_file)
        script_content = script_content.replace('OUTPUT', output_file)

        # Write the modified script to a temporary file
        temp_script_path = 'temp_job_script.sh'
        with open(temp_script_path, 'w') as temp_file:
            temp_file.write(script_content)

        os.chmod(temp_script_path, 0o755)  # Make the script executable
        return temp_script_path
    
    def create_shell_script(self, script_path, cluster=False):
        """
        Create a shell script to deactivate Conda environment, run pymolcas, and reactivate Conda environment.
        """

        orca_path = shutil.which('orca')
        if self.program == 'MOLCAS':
            script_content = f"""#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
# Deactivate the current Conda environment
conda deactivate

# Run the pymolcas command
pymolcas -nt 32 $1 > $2

# Reactivate the original Conda environment
conda activate vlxenv
export LD_LIBRARY_PATH=~/minconda3/envs/vlxenv/lib/
"""
        elif self.program == 'ORCA' and orca_path is not None:
            script_content = f"""#!/bin/bash

# Run the orca command
source /home/vlind06/miniconda3/etc/profile.d/conda.sh
conda deactivate
{orca_path} $1 > $2
conda activate vlxenv_simd_master
"""

        elif self.program == 'ORCA' and orca_path is not None and cluster == True:
            script_content = self.jobscript
        elif self.program == 'QCHEM':
            script_content = self.jobscript
        else:
            raise ValueError("Unsupported program or ORCA executable not found")
            
        with open(script_path, 'w') as script_file:
            script_file.write(script_content)
 
        # Make the script executable
        print('SKRIPT IS BEING MADE EXE', script_path)
        os.chmod(script_path, 0o755)
     
    
    def create_ghost_atom_line(self, coordinates, atom_label="Ne"):
        """
        Create a formatted string for a ghost atom at the centroid.

        Parameters:
        - coordinates: np.ndarray of shape (N_atoms, 3)
        - atom_label: str, default "Ne"

        Returns:
        - insert_line: str, e.g., "Ne: 1.234 2.345 3.456\n"
        """
        centroid = np.mean(coordinates, axis=0)
        insert_line = f"{atom_label}: {centroid[0]:.6f} {centroid[1]:.6f} {centroid[2]:.6f}\n"
        return insert_line
    
    def compute_energy(self, molecule, basis):
        
        self.basis_set_label = basis
        molecule.write_xyz_file('current_geometry.xyz')
        
        if self.add_ghost_atom:

            
            with open('current_geometry.xyz', 'r') as f:
                
                

                insert_line = self.create_ghost_atom_line(molecule.get_coordinates_in_angstrom())
                lines = f.readlines()

                # Insert the line at the 3rd index (line 4)
                lines.insert(2, insert_line)  # Index 3 = after line 3 (0-based indexing)

                # Write back the modified content
                with open('current_geometry.xyz', 'w') as f:
                    f.writelines(lines)

        
        lines = None
        with open('current_geometry.xyz', 'r') as f1:
            lines = f1.readlines()
        with open(self.xyz_structures, 'a') as f2:
            f2.writelines(lines)

        self.xyz_filename = 'current_geometry.xyz'

        """
        Run the pymolcas command with the input file and redirect output to the output file.
        """
        if self.program in ['MOLCAS', 'ORCA'] and self.cluster_manager is None:
            if self.program == 'MOLCAS':
                self.create_input_file_energy(self.roots)
                script_path = os.path.join(os.getcwd(), "run_pymolcas.sh")
                self.create_shell_script(script_path)

                bash_command = f"bash -c {script_path} {self.input_filename} {self.output_filename}"
                print(f"Executing command: {bash_command}")

                current_path = os.getcwd()
                workdir_path = os.path.join(current_path, "current_data")
                os.makedirs(workdir_path, exist_ok=True)
                print(f"MOLCAS_WORKDIR set to {workdir_path}")
                os.environ['MOLCAS_WORKDIR'] = workdir_path

            elif self.program == 'ORCA':
                
                self.create_input_file_energy(self.input_files[0])
                script_path = os.path.join(os.getcwd(), "run_orca.sh")
                self.create_shell_script(script_path)
                bash_command = f"bash -c '{script_path} {self.input_files[0]} {self.output_files[0]}'"

            try:
                clean_env = os.environ.copy()

                # Remove Conda-specific environment variables
                conda_vars = ["CONDA_PREFIX", "CONDA_SHLVL", "CONDA_DEFAULT_ENV"]
                for var in conda_vars:
                    clean_env.pop(var, None)  # Remove if exists

                # Adjust PATH to remove Conda-specific paths
                if "PATH" in clean_env:
                    clean_env["PATH"] = ":".join(
                        p for p in clean_env["PATH"].split(":") if "miniconda3" not in p
                    )
                result = subprocess.run(
                bash_command,
                shell=True,
                executable="/bin/bash",
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=clean_env,  # Inherit the current environment
                )

            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
        
        elif self.program in ['MOLCAS', 'ORCA'] and self.cluster_manager is not None:
            if self.program == 'ORCA':
                print('Here is the files list', self.input_files)
                self.create_input_file_energy(self.input_files[0])
                template_path = 'jobscript_orca'
                script_path = self.prepare_job_script(template_path, self.input_files[0], self.output_files[0])
                print('Jobscript is being created')
                current_energy_commit_id = self.cluster_manager.run_cluster_job_and_wait(
                script_path, input=self.input_files[0], geometry=self.xyz_filename, queue=False)
                self.cluster_manager.connect_copy(self.output_files[0])

                    
        elif self.program == 'QCHEM':
            
            self.create_input_file_energy(self.input_files[0])
            template_path = 'jobscript_qchem'
            script_path = self.prepare_job_script(template_path, self.input_files[0], self.output_files[0])
            print('Jobscript is being created')
            current_energy_commit_id = self.cluster_manager.run_cluster_job_and_wait(
            script_path, queue=False) 

        energies = self.extract_energies()
        self.energies = energies

        return energies

    
    def extract_energies(self):
        """
        Extract energy values for the given number of roots from the output file.
        """
        energies = []
        try:
            if self.program == 'MOLCAS':
                with open(self.output_filename, 'r') as file:
                    content = file.read()

                for root in range(1, num_roots + 1):
                    match = re.search(f'printout of CI-coefficients larger than  0.05 for root  {root}.*?energy=\\s*([-0-9.]+)', content, re.DOTALL)
                    if match:
                        energy = float(match.group(1))
                        energies.append(energy)
                    else:
                        print(f"Energy for root {root} not found.")
                        energies.append(None)
            
            elif self.program == 'ORCA':
                print('Here is the file that is being opend', self.output_files[0])
                with open(self.output_files[0], 'r', encoding="utf-8", errors="ignore") as file:
                    content = file.read()
                # Extract ground state energy
                ground_state_match = re.search(r'FINAL SINGLE POINT ENERGY\s*([-0-9.]+)', content)
                if ground_state_match:
                    energies.append(float(ground_state_match.group(1)))
                else:
                    print("Ground state energy not found.")
        
            elif self.program == 'QCHEM':
                energies = []
                spin_multplicity = None
                sf_dft_section = False

                with open(self.output_filename, 'r') as file:
                    content = file.read()

                # Check if this is an SF-DFT calculation
                if 'SF-DFT Excitation Energies' in content:
                    sf_dft_section = True

                if sf_dft_section:
                    # Extract ground state energy (first "excited" state in SF-DFT)
                    
                    ground_state_match = re.search(r'Total energy for state\s+1:\s+([-0-9.]+)\s+au', content)
                    if ground_state_match:
                        energies.append(float(ground_state_match.group(1)))
                    else:
                        print("Ground state energy not found.")
                else:
                    # If not SF-DFT, extract SCF energy as ground state
                    scf_energy_match = re.search(r'SCF\s+energy in the final basis set\s*=\s*([-0-9.]+)', content)
                    if scf_energy_match:
                        energies.append(float(scf_energy_match.group(1)))
                    else:
                        print("SCF energy not found.")
            
                print('Here is the spin_multi', spin_multplicity)

            reordered_energies = energies
            print('Here are the enegies', energies, reordered_energies, self.tracked_roots)

            return reordered_energies
        except Exception as e:
            print(f"Error extracting energies: {e}")
            return None
    
    
    def create_input_file_energy(self, input_file):
        """
        Create the specific input file with the given content.
        """
        try:
            if self.program == 'MOLCAS':
                with open(self.input_filename, 'w') as file:
                    file.write("&GATEWAY\n")
                    file.write(f" coord={self.xyz_filename}\n")
                    file.write(" group=NoSym\n")
                    file.write(" basis=6-31G*\n")
                    file.write("End of input\n\n")

                    file.write("&SEWARD\n")
                    file.write("Cholesky\n")
                    file.write("End of input\n\n")

                    file.write("&SCF\n")
                    file.write(" spin=1\n")
                    file.write(" charge=0\n")
                    file.write("End of input\n\n")

                    file.write("&RASSCF\n")
                    file.write(" Spin=1\n")
                    file.write(" Charge=0\n")
                    file.write(" NActEl=2 0 0\n")
                    file.write(" RAS2=2\n")
                    file.write(f" CIRoot={roots} {roots} 1\n")
                    file.write(" SDAV=1000\n")
                    file.write(" THRPr=0.0001\n")
                    file.write("End of input\n\n")

                    file.write("&LOCALISATION\n")
                    file.write("Boys\n")
                    file.write("End of input\n\n")

            elif self.program == 'ORCA':
                full_path = os.path.abspath(self.xyz_filename)
                if self.path_on_cluster is not None:
                    full_path = f'{self.path_on_cluster}/{self.xyz_filename}'
                with open(input_file, 'w') as file:
                    if self.solvation[0] is True:
                        file.write(f'!{self.method} {self.xc_func} {self.dispersion} {self.basis_set_label} {self.solvation[1]}({self.solvation[2]})\n')
                    else:
                        file.write(f'!{self.method} {self.xc_func} {self.dispersion} {self.basis_set_label}\n')
                    file.write(f'%maxcore 3000\n')
                    file.write(f'%PAL\n')
                    file.write(f'nprocs {self.nprocs}\n')
                    file.write('END\n')
                    file.write(f'* xyz {self.charge} {self.spin}\n')
                    with open(full_path, 'r') as geometry_lines:
                        for i, line in enumerate(geometry_lines):
                            if i < 2:
                                continue
                            file.write(line)
                    file.write('*\n')

            elif self.program == 'QCHEM':
                functional = 'b3lyp'
                if self.spin_flip is True:
                    self.spin = 3
                    functional = 'bhhlyp'
                full_path = os.path.abspath(self.xyz_filename)
                with open(input_file, 'w') as file:
                    file.write('$molecule\n')
                    file.write(f'{self.charge} {self.spin}\n')
                    with open(self.xyz_filename, 'r') as xyz_file:
                        # Skip the first two lines (atom count and comment)
                        next(xyz_file)
                        next(xyz_file)
                        # Write the remaining lines (atom coordinates)
                        for line in xyz_file:
                            file.write(line)
                    file.write('$end\n\n')
                    
                    file.write('$rem\n')
                    file.write(f'EXCHANGE        {functional}\n')
                    file.write('BASIS           6-311G(d,p)\n')
#                    file.write('SCF_GUESS       core\n')
                    file.write('UNRESTRICTED     true\n')
                    file.write(f'SPIN_FLIP      {self.spin_flip}\n')
                    file.write(f'CIS_N_ROOTS    {root}\n')
                    file.write('CIS_CONVERGENCE 0\n')
                    file.write(f'$end\n')



            print(f"Input file '{self.input_filename}' created successfully.")
        except Exception as e:
            print(f"Error creating input file: {e}")


