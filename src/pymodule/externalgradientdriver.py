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
from .externalscfdriver import ExternalScfDriver
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


class ExternalGradientDriver:

    def __init__(self, qm_driver, comm=None, ostream=None):
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

        self.do_four_point = False
        self.delta_h = 1e-3
        self.natoms = None

        self.xyz_filename = None
        self.input_files = []
        self.output_files = []
        self.input_filename = 'current_input_gradient.inp'
        self.output_filename = 'current_output_gradient.out'
        self.input_files.append(self.input_filename)
        self.output_files.append(self.output_filename)   
    
        self.program = qm_driver.program
        self.cluster_manager = qm_driver.cluster_manager
        self.NAC = qm_driver.NAC
        print('method in the gradient driver', qm_driver.method)
        self.method = qm_driver.method
        self.qm_driver = qm_driver    
    
    def prepare_job_script(self, template_path, input_file, output_file):
        with open(template_path, 'r') as template_file:
            script_content = template_file.read()

        # Replace placeholders with actual file names
        script_content = script_content.replace('INPUT', input_file)
        script_content = script_content.replace('OUTPUT', output_file)

        # Write the modified script to a temporary file
        temp_script_path = 'temp_jobscript'
        with open(temp_script_path, 'w') as temp_file:
            temp_file.write(script_content)

        #os.chmod(temp_script_path, 0o755)  # Make the script executable
        return temp_script_path
    
    def create_shell_script(self, script_path):
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

source $(conda info --base)/etc/profile.d/conda.sh

conda deactivate
# Run the orca command
{orca_path} $1 > $2
conda activate vlxenv_new_compile

"""
        else:
            raise ValueError("Unsupported program or ORCA executable not found")
            
        with open(script_path, 'w') as script_file:
            script_file.write(script_content)

        # Make the script executable
        os.chmod(script_path, 0o755)
        

    def extract_gradients(self):
        """
        Extract gradients from the output file.
        """
        gradients = []
        try:
            if self.program == 'MOLCAS':
                with open(self.full_output_filename, 'r', encoding='utf-8') as file:
                    lines = iter(file.readlines())
                    in_gradient_section = False
                    current_gradients = []
                    for line in lines:
                        if re.search(r'Numerical gradient, root', line):
                            if in_gradient_section:
                                gradients.append(np.array(current_gradients))
                                current_gradients = []
                            in_gradient_section = True
                            # Skip the next three lines (separator, header, separator)
                            next(lines)
                            next(lines)
                            next(lines)
                            continue
                        if in_gradient_section:
                            if re.match(r'^  -+', line):
                                gradients.append(np.array(current_gradients))
                                current_gradients = []
                                in_gradient_section = False
                            else:
                                parts = line.split()
                                if len(parts) == 4:
                                    _, x, y, z = parts
                                    current_gradients.append([float(x), float(y), float(z)])

                # Ensure the last set of gradients is added
                if current_gradients:
                    gradients.append(np.array(current_gradients))
            
            if self.program == 'ORCA':
                print('I am going in the ORCA part')
                gradients_orca = {'ground_state': None}

                print('ROOT', self.output_files[0])
                with open(self.output_files[0], 'r', encoding='utf-8') as file:
                    lines = iter(file.readlines())
                    in_ground_state_section = True
                    current_gradients = []
                    
                    for line in lines:
                        
                        # Check for the start of the ground state section
                        #else:
                        #    in_ground_state_section = True
                        # Check for the end of the ground state section
                        if in_ground_state_section and re.search(r'CARTESIAN GRADIENT', line):
                            # Skip the next line (header separator)
                            next(lines)
                            next(lines)
                            print(line)
                            in_ground_state_section = False
                            for gradient_line in lines:
                                if re.match(r'^\s*$', gradient_line):  # End of the gradient section
                                    break
                                parts = gradient_line.split()
                                if len(parts) == 6:
                                    _, _, _, x, y, z = parts
                                    current_gradients.append([float(x), float(y), float(z)])
                            gradients_orca['ground_state'] = np.array(current_gradients)
                            continue


                gradients.append(gradients_orca['ground_state'])
                
            if self.program == 'QCHEM':
                gradients_qchem = {'ground_state': None}
                
                for root in range(1):
                    current_state = None
                    current_gradient = []
                    with open(self.output_files[root], 'r', encoding='utf-8') as file:
                        lines = iter(file.readlines())
                        for line in lines:
                            # Check for the CIS state line
                            if self.qm_driver.spin_flip == True:
                                match = re.search(r'CIS (\d+) State Energy is', line)
                                if match:
                                    state_number = int(match.group(1))
                                    
                                    # Determine if this is ground state or excited state based on self.method
                                    if self.qm_driver.spin_flip is True:
                                        current_state = 'ground_state' if state_number == 1 else f'excited_state_{root + 1}'
                                    else:
                                        current_state = f'excited_state_{root + 1}'  # Assume excited state for non-SF methods
                                    
                                    # Skip the next line (header)
                                    next(lines)
                                    next(lines)
                                    
                                    # Read the gradient
                                    for _ in range(3):  # Assuming 3 lines of gradient data
                                        
                                        gradient_line = next(lines)
                                        values = gradient_line.split()[1:]  # Skip the first column (index)
                                        current_gradient.append([float(val) for val in values])
                                    
                                    # Store the gradient
                                    transposed_data = list(zip(*current_gradient))

                                    # Flatten the transposed data

                                    print('Current gradient', transposed_data)
                                    gradients_qchem[current_state] = np.array(transposed_data)
                            else:
                                match = re.search(r'Gradient of SCF Energy', line)
                                if match:
                                    next(lines)

                                    # Read the gradient
                                    for _ in range(3):  # Assuming 3 lines of gradient data

                                        gradient_line = next(lines)
                                        values = gradient_line.split()[1:]  # Skip the first column (index)
                                        current_gradient.append([float(val) for val in values])

                                    # Store the gradient
                                    transposed_data = list(zip(*current_gradient))

                                    # Flatten the transposed data

                                    print('Current gradient', transposed_data)
                                    gradients_qchem['ground_state'] = np.array(transposed_data)

                    gradients.append(gradients_qchem['ground_state'])

            print(gradients)
            return gradients
        except Exception as e:
            print(f"Error extracting gradients: {e}")
            return None
        
    
    
    def compute_gradient(self, current_molecule, molecule_geometry=None, paralell=False):
        
        """
        Run the pymolcas command with the input file and redirect output to the output file.
        """

       
        self.full_output_filename = self.output_filename
        if self.program in ['MOLCAS', 'ORCA'] and self.cluster_manager is None:
            if self.program == 'MOLCAS':
                self.create_input_file_gradient(self.roots)

                script_path = os.path.join(os.getcwd(), "run_pymolcas.sh")
                self.create_shell_script(script_path)

                bash_command = f"bash {script_path} {self.input_filename} {self.full_output_filename}"
                print(f"Executing command: {bash_command}")

                current_path = os.getcwd()
                workdir_path = os.path.join(current_path, "current_data_gradient")
                os.makedirs(workdir_path, exist_ok=True)
                print(f"MOLCAS_WORKDIR set to {workdir_path}")
                os.environ['MOLCAS_WORKDIR'] = workdir_path
            
            elif self.program == 'ORCA':
                self.create_input_file_gradient(self.input_files[0])
                script_path = os.path.join(os.getcwd(), "run_orca.sh")
                self.create_shell_script(script_path)
                bash_command = f"bash -c '{script_path} {self.input_files[0]} {self.output_files[0]}'"
                print(f"Executing command: {bash_command}")
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
                    #return e.stderr.decode('utf-8')
        
        elif self.program in ['MOLCAS', 'ORCA'] and self.cluster_manager is not None:

            if self.program == 'ORCA':
                print('Here is the files list', self.input_files)

                self.create_input_file_gradient(self.input_files[0])
                template_path = 'jobscript_orca'
                script_path = self.prepare_job_script(template_path, self.input_files[0], self.output_files[0])
                print('Jobscript is being created')
                current_gradient_commit_id = self.cluster_manager.run_cluster_job_and_wait(
                script_path, queue=False)
            
                self.qm_driver.job_ids.append(current_gradient_commit_id)

        elif self.program == 'QCHEM':
            
            for root in range(1):
                    self.create_input_file_gradient(self.input_files[root], root + 1)
                    template_path = 'jobscript_qchem'
                    script_path = self.prepare_job_script(template_path, self.input_files[root], self.output_files[root])
                    print('Jobscript is being created')

                    current_gradient_commit_id = self.cluster_manager.run_cluster_job_and_wait(
                    script_path, queue=False)

            self.qm_driver.job_ids.append(current_gradient_commit_id)
            print('commited jobs are finished')
        
        #if not paralell:
        #    print('################ extracting')
        #    gradients = self.extract_gradients()
        #    print('in not parallel', gradients)
        #NAC = None
        #if self.NAC and not molecule_geometry:
        #    NAC = self.extract_NACs()
        #print('NAC', NAC, gradients)
        #return gradients, NAC, current_job_id

    def create_input_file_gradient(self, input_file):
        """
        Create the specific input file with the given content.
        """
        print('Giving of the current input file', self.input_filename)
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

                    file.write('&ALASKA\n')
                    #file.write('PNEW\n\n')


                    file.write("&LOCALISATION\n")
                    file.write("Boys\n")
                    file.write("End of input\n\n")
            
            elif self.program == 'ORCA':

                full_path = os.path.abspath(self.qm_driver.xyz_filename)
                with open(input_file, 'w') as file:
                    if self.qm_driver.solvation[0] is True:
                        file.write(f'!{self.qm_driver.method} {self.qm_driver.xc_func} {self.qm_driver.dispersion} {self.qm_driver.basis_set_label} ENGRAD {self.qm_driver.solvation[1]}({self.qm_driver.solvation[2]})\n')
                    else:
                        file.write(f'!{self.qm_driver.method} {self.qm_driver.xc_func} {self.qm_driver.dispersion} {self.qm_driver.basis_set_label} ENGRAD\n')
                    
                    file.write(f'%maxcore 3000\n')
                    file.write(f'%PAL\n')
                    file.write(f'nprocs {self.qm_driver.nprocs * 2}\n')
                    file.write('END\n')
                    file.write(f'* xyz {self.qm_driver.charge} {self.qm_driver.spin}\n')
                    with open(full_path, 'r') as geometry_lines:
                        for i, line in enumerate(geometry_lines):
                            if i < 2:
                                continue
                            file.write(line)
                    file.write('*\n')
                
            elif self.program == 'QCHEM':
                    
                full_path = os.path.abspath(self.xyz_filename)
                with open(input_file, 'w') as file:
                    file.write('$molecule\n')
                    file.write(f'{self.qm_driver.charge} {self.qm_driver.spin}\n')
                    with open(self.xyz_filename, 'r') as xyz_file:
                        # Skip the first two lines (atom count and comment)
                        next(xyz_file)
                        next(xyz_file)
                        # Write the remaining lines (atom coordinates)
                        for line in xyz_file:
                            file.write(line)
                    file.write('$end\n\n')
                    if root == 1 and self.qm_driver.spin_flip == False:
                        
                        file.write('$rem\n')
                        file.write('JOBTYPE         force\n')
                        file.write('EXCHANGE        b3lyp\n')
                        file.write('BASIS           6-311G(d,p)\n')
                        file.write('SCF_GUESS       core\n')
                        file.write('$end\n')
                   
            print(f"Input file '{self.input_filename}' created successfully.")
        except Exception as e:
            print(f"Error creating input file: {e}")


