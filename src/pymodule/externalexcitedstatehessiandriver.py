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
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import re
import json
import shutil
from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .molecule import Molecule
from .externalscfdriver import ExternalScfDriver
from .externalgradientdriver import ExternalGradientDriver
from .externalexcitedstatedriver import ExternalExcitedStatesScfDriver
from .externalexcitedstategradientdriver import ExternalExcitedStatesGradientDriver
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


class ExternalExcitedStatesHessianDriver:

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

        self.roots_to_follow = qm_driver.roots
        self.do_four_point = False
        self.delta_h = 1e-3

        self.xyz_filename = None
        self.input_files = []
        self.output_files = []
        self.hessian_files = []
        
        for root in self.roots_to_follow:
            self.input_filename_excited = f'current_input_hessian_excited_{root}.inp'
            self.output_filename_excited = f'current_output_hessian_excited_{root}.out'
            self.hessian_file_excited = f'current_input_hessian_excited_{root}.hess'
            self.input_files.append(self.input_filename_excited)
            self.output_files.append(self.output_filename_excited)
            self.hessian_files.append(self.hessian_file_excited)

        self.analytical = False
        self.qm_driver = qm_driver
        self.cluster_manager = qm_driver.cluster_manager
        self.program = qm_driver.program
        self.qm_driver = qm_driver
        self.NAC = False
    

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

"""
        elif self.program == 'OpenQP':
            script_content = f"""#!/bin/bash

# Run the orca command
source $(conda info --base)/etc/profile.d/conda.sh
conda deactivate
export OMP_NUM_THREADS={self.qm_driver.nprocs * 3}
export OPENQP_ROOT={self.qm_driver.open_qp_path}
export MKL_INTERFACE_LAYER="@_MKL_INTERFACE_LAYER@"
export MKL_THREADING_LAYER=SEQUENTIAL
openqp $1 > $2
"""
        else:
            raise ValueError("Unsupported program or ORCA executable not found")
            
        with open(script_path, 'w') as script_file:
            script_file.write(script_content)

        # Make the script executable
        os.chmod(script_path, 0o755)

    
    def extract_hessians(self):
        """
        Extract gradients from the output file.
        """
        hessians = []
        try:
            if self.program == 'ORCA':
                    
                hessians_orca = {f'excited_state_{i+1}': None for i in self.roots_to_follow}
                for root in range(len(self.roots_to_follow)):
                    print('I am going in the ORCA part', self.hessian_files[root])
                    with open(self.hessian_files[root], 'r', encoding='utf-8') as file:
                        lines = file.readlines()
                        in_hessian_section = False
                        matrix_size = 0
                        matrix = None
                        column_indices = []

                        line_iter = iter(lines)
                        for line in line_iter:
                            line = line.strip()
                            if not in_hessian_section:
                                if "orca_hessian_file" in line:
                                    continue  # Skip this specific line
                                if "hessian" in line.lower():
                                    in_hessian_section = True
                                continue

                            if not line:
                                continue

                            if matrix_size == 0:
                                if line.isdigit():
                                    matrix_size = int(line)
                                    matrix = np.zeros((matrix_size, matrix_size))
                                    continue
                                else:
                                    continue  # Skip any lines until we get the matrix size

                            # Now we have matrix_size and we are in the hessian section

                            if '$' in line:
                                # End of hessian section
                                break

                            # Function to check if the line contains only integers (column indices)
                            def is_column_indices_line(line):
                                parts = line.strip().split()
                                return all(part.isdigit() for part in parts)

                            if is_column_indices_line(line):
                                column_indices = [int(x) for x in line.strip().split()]
                                continue

                            parts = line.strip().split()
                            if parts[0].isdigit():
                                row_index = int(parts[0])
                                values = [float(x) for x in parts[1:]]
                                for col_idx, value in zip(column_indices, values):
                                    matrix[row_index, col_idx] = value
                            else:
                                # Not a data line, skip
                                continue

                        hessians_orca[f'excited_state_{self.roots_to_follow[root]}'] = matrix
                        hessians.append(hessians_orca[f'excited_state_{self.roots_to_follow[root]}'])
                        print('Hessian matrix extracted successfully.')
                        if matrix is None:
                            print('No hessian data found')
                            exit()
                
            elif self.program == 'OpenQP':

                hessians_orca = {f'excited_state_{i+1}': None for i in self.roots_to_follow}

                for root_idx, root in enumerate(self.roots_to_follow):
                    p_in = Path(self.input_files[root_idx])    
                    new_output_file = str(Path(p_in).with_suffix(".hess.json"))
                    with open(new_output_file, "r", encoding="utf-8") as f:
                        data = json.load(f)

                    H = np.array(data['hessian'], dtype=float)

                    # If it was saved flattened (just in case), reshape to square
                    if H.ndim == 1:
                        n = int(round(H.size ** 0.5))
                        if n * n != H.size:
                            raise ValueError(f"{'hessian'} has length {H.size}, not a perfect square.")
                        H = H.reshape(n, n)

                    # Optional: enforce symmetry (useful if minor numeric noise exists)
                    if not np.allclose(H, H.T, atol=1e-10):
                        H = 0.5 * (H + H.T)
                    hessians_orca[f'excited_state_{self.roots_to_follow[root_idx]}'] = H
                    hessians.append(hessians_orca[f'excited_state_{self.roots_to_follow[root_idx]}'])
            
            elif self.program == 'QCHEM':

                
                hessians_qchem = {f'excited_state_{i+1}': None for i in range(self.roots)}
                for root in range(self.roots + 1):
                    print('I am going in the ORCA part', self.output_files[root])
                    with open(self.output_files[root], 'r', encoding='utf-8') as file:
                        lines = file.readlines()
                        in_hessian_section = False
                        in_orientation_section = False
                        matrix_size = 0
                        matrix = None
                        column_indices = []

                        line_iter = iter(lines)
                        for line in line_iter:
                            
                            if matrix_size == 0 and 'Standard Nuclear Orientation (Angstroms)' in line:
                                in_orientation_section = True
                                # Skip the next two header lines
                                next(line_iter)
                                next(line_iter)
                                continue
                            if in_orientation_section:
                                if '---' in line:
                                    # End of the orientation section
                                    in_orientation_section = False
                                    continue
                                else:
                                    fields = line.strip().split()
                                    if fields and fields[0].isdigit():
                                        number_of_atoms = int(fields[0]) 
                                        matrix_size = number_of_atoms * 3
                                        matrix = np.zeros((matrix_size, matrix_size))
                                        continue


                            line = line.strip()
                            if not in_hessian_section:
                                if "Final Hessian." in line:
                                    in_hessian_section = True
                                continue
            
                            # Now we have matrix_size and we are in the hessian section

                            if '----' in line and in_hessian_section is True:
                                # End of hessian section
                                break
                                    
                            # Function to check if the line contains only integers (column indices)
                            def is_column_indices_line(line):
                                parts = line.strip().split()
                                return all(part.isdigit() for part in parts)
                            #print('LINE', line, in_hessian_section)
                            if is_column_indices_line(line):
                                column_indices = [(int(x) -1) for x in line.strip().split()]
                                continue
                            
                            parts = line.strip().split()
                            if parts[0].isdigit():
                                row_index = int(parts[0]) - 1
                                values = [float(x) for x in parts[1:]]
                                for col_idx, value in zip(column_indices, values):
                                    matrix[row_index, col_idx] = value
                            else:
                                # Not a data line, skip
                                continue

                    hessians_qchem[f'excited_state_{root + 1}'] = matrix
                    hessians.append(hessians_qchem[f'excited_state_{root + 1}'])
                    print('Hessian matrix extracted successfully.')
                    if matrix is None:
                        print('No hessian data found')
                        exit()




            current_path = os.getcwd()
        
            pattern = os.path.join(current_path, 'current_*')
            files_to_remove = glob.glob(pattern)
            
            for path in files_to_remove:
                try:
                    if os.path.isfile(path):
                        os.remove(path)
                        print(f"Removed file: {path}")
                    elif os.path.isdir(path):
                        shutil.rmtree(path)   # removes folder and all its contents
                        print(f"Removed folder: {path}")
                    else:
                        print(f"Path not found: {path}")
                except Exception as e:
                    print(f"Could not remove {path}: {e}")

            return hessians
        except Exception as e:
            print(f"Error extracting gradients: {e}")
            return None
    
    def compute_analytical_hessian(self, current_molecule, paralell=False, molecule_geometry=None):
        
        """
        Run the pymolcas command with the input file and redirect output to the output file.
        """

        current_molecule.write_xyz_file('current_geometry.xyz')
        self.xyz_filename = 'current_geometry.xyz'
        current_job_id = None

        if molecule_geometry:

            self.full_output_filename = self.output_fils[0] + molecule_geometry
            current_molecule.write_xyz_file(f'current_geometry_{molecule_geometry}.xyz')
            self.xyz_filename = f'current_geometry_{molecule_geometry}.xyz'
        
        else:
            self.full_output_filename = self.output_files[0]
        if self.program in ['MOLCAS', 'ORCA', 'OpenQP'] and self.cluster_manager is None:
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
                
                for root in range(len(self.roots_to_follow)):

                    self.create_input_file_hessian(self.input_files[root], self.roots_to_follow[root])

                    script_path = os.path.join(os.getcwd(), "run_orca.sh")
                    self.create_shell_script(script_path)

                    bash_command = f"bash -c '{script_path} {self.input_files[root]} {self.output_files[root]}'"
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
                jobs_finished = 0

            elif self.program == 'OpenQP':

                script_path = None
                for root_idx, root in enumerate(self.roots_to_follow):
                    # p_in = Path(self.input_files[0])
                    # p_out = Path(self.output_files[0])
                    # new_input_file = p_in.with_suffix(".inp").with_stem(f"{p_in.stem}_{root}")
                    # new_output_file = p_out.with_suffix(".out").with_stem(f"{p_out.stem}_{root}")
                    inp0 = Path(self.input_files[root_idx]).resolve()
                    out0 = Path(self.output_files[root_idx]).resolve()
                    self.create_input_file_hessian(inp0, root)
                    script_path = os.path.join(os.getcwd(), "run_openqp.sh")
                    script_path = Path(script_path).resolve()
                    self.create_shell_script(script_path)
                    bash_command = f"'{script_path} {inp0} {out0}'"

                    cmd = f"source /etc/profile; source ~/.bashrc >/dev/null 2>&1 || true; /bin/bash '{script_path}' '{inp0}' '{out0}'"

                    # Start from the current environment; only strip conda if you *must*
                    env = os.environ.copy()

                    try:
                        result = subprocess.run(
                            ["/bin/bash", "-lc", cmd],
                            cwd=str(script_path.parent),
                            env=env,
                            check=True,
                            text=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                        )
                        print(result.stdout)
                    except subprocess.CalledProcessError as e:
                        print("Command failed.")
                        print("STDOUT:\n", e.stdout)
                        print("STDERR:\n", e.stderr)
                        raise
                jobs_finished = 0
                
        elif self.program in ['MOLCAS', 'ORCA'] and self.cluster_manager is not None:
            for root in range(self.roots + 1):

                self.create_input_file_hessian(self.input_files[root], self.qm_driver.tracked_roots[root])
                # Cluster execution for QCHEM
                template_path = 'jobscript_orca'  # Path to your template file
                script_path = self.prepare_job_script(template_path, self.input_files[root], self.output_files[root])
                
                # Run the job on the cluster and wait for completion
                current_hessian_commit_id = self.cluster_manager.run_cluster_job_and_wait(
                    script_path, queue=False)
               
                self.qm_driver.job_ids.append(current_hessian_commit_id)
                print('Here is the jobid list', self.qm_driver.job_ids)

            jobs_finished = None
            if self.cluster_manager is not None and 1 == 0:
                jobs_finished = self.cluster_manager.check_queue(self.qm_driver.job_ids)
            

        elif self.program == 'QCHEM' and self.cluster_manager is not None:
            
            for root in range(self.roots + 1):
                    self.create_input_file_hessian(self.input_files[root], self.qm_driver.tracked_roots[root])
                    template_path = 'jobscript_qchem'
                    script_path = self.prepare_job_script(template_path, self.input_files[root], self.output_files[root])
                    print('Jobscript is being created')

                    current_gradient_commit_id = self.cluster_manager.run_cluster_job_and_wait(
                    script_path, queue=False)

            self.qm_driver.job_ids.append(current_gradient_commit_id)

            jobs_finished = None
            if self.cluster_manager is not None and 1 == 0:
                jobs_finished = self.cluster_manager.check_queue(self.qm_driver.job_ids)

        print('calcualtion is over')
 
        return jobs_finished 
    
    # def compute_numerical_hessian(self, molecule):
    #     """
    #     Performs calculation of numerical Hessian.

    #     :param molecule:
    #         The molecule.
    #     """

    #     self.ostream.mute()

    #     # atom labels
    #     labels = molecule.get_labels()
        
    #     # number of atoms
    #     natm = molecule.number_of_atoms()

    #     # atom coordinates (nx3)
    #     coords = molecule.get_coordinates_in_bohr()

    #     # charge and spin multiplicity
    #     charge = molecule.get_charge()
    #     multiplicity = molecule.get_multiplicity()

    #     if self.rank == mpi_master():
    #         # Hessian
    #         hessian = np.zeros((self.roots, natm, 3, natm, 3))

    #         # numerical dipole gradient (3 dipole components,
    #         # no. atoms x 3 atom coords)
    #         self.dipole_gradient = np.zeros((3, 3 * natm))

    #     if not self.do_four_point:
    #         for i in range(natm):
    #             for d in range(3):
    #                 gradient_driver = ExternalQMGradientDriver(self.qm_driver)
    #                 coords[i, d] += self.delta_h
    #                 new_mol = Molecule(labels, coords, units='au')
    #                 new_mol.set_charge(charge)
    #                 new_mol.set_multiplicity(multiplicity)
    #                 current_string = f'{i}_{d}_plus'
    #                 # create a new XTB driver object;
    #                 # without this the energy is always zero...;
    #                 gradient, _, _ = gradient_driver.compute_gradient(new_mol, current_string)
    #                 grad_plus = gradient

    #                 coords[i, d] -= 2.0 * self.delta_h
    #                 new_mol = Molecule(labels, coords, units='au')
    #                 new_mol.set_charge(charge)
    #                 new_mol.set_multiplicity(multiplicity)
    #                 current_string = f'{i}_{d}_minus'

    #                 gradient, _, _ = gradient_driver.compute_gradient(new_mol, current_string)

    #                 grad_minus = gradient

    #                 print('Gradient_plus', grad_plus)
    #                 coords[i, d] += self.delta_h
    #                 if self.rank == mpi_master():
    #                     for root in range(self.roots):
    #                         hessian[root, i, d, :, :] = (grad_plus[root] - grad_minus[root]) / (2.0 * self.delta_h)

    #     reshaped_hessian = np.zeros((self.roots, 3 * natm, 3 * natm))

    #     # reshaped Hessian as member variable
    #     if self.rank == mpi_master():
    #         for root in range(self.roots):

    #             flat_hessian = hessian[root].flatten()

    #             # Ensure the total number of elements matches the desired shape
    #             expected_elements = (natm * 3)**2
    #             assert flat_hessian.size == expected_elements, "Total elements do not match the expected"
    #             reshaped_hessian[root] = flat_hessian.reshape(3 * natm, 3 * natm)
    #             reshaped_hessian[root] = ( reshaped_hessian[root] + reshaped_hessian[root].T ) / 2.0

    #         self.hessian = reshaped_hessian
    #     else:
    #         self.hessian = None

    #     self.ostream.unmute()
         
    #     print(self.hessian, self.hessian.shape, np.allclose(self.hessian[0], self.hessian[0].T))
    #     return reshaped_hessian
    
    #TODO Try to serialize that code in order to run multiple gradient calculations at the same time
    #in order to do so I need to do:
    #   - write the clustermanager so it keeps track of finished jobs
    #   - put the right elements into the correct elements in the matrix
    #   - merge everything together and let the program continue locally

    # def compute_numerical_hessian_parallel(self, molecule):
    #     """
    #     Performs calculation of numerical Hessian.

    #     :param molecule:
    #         The molecule.
    #     """

    #     self.ostream.mute()

    #     # atom labels
    #     labels = molecule.get_labels()
        
    #     # number of atoms
    #     natm = molecule.number_of_atoms()

    #     # atom coordinates (nx3)
    #     coords = molecule.get_coordinates_in_bohr()

    #     # charge and spin multiplicity
    #     charge = molecule.get_charge()
    #     multiplicity = molecule.get_multiplicity()

    #     if self.rank == mpi_master():
    #         # Hessian
    #         hessian = np.zeros((self.roots, natm, 3, natm, 3))

    #         # numerical dipole gradient (3 dipole components,
    #         # no. atoms x 3 atom coords)
    #         self.dipole_gradient = np.zeros((3, 3 * natm))
    #     job_ids = []
    #     if not self.do_four_point:
    #         for i in range(natm):
    #             for d in range(3):
    #                 gradient_driver = ExternalQMGradientDriver(self.qm_driver)
    #                 coords[i, d] += self.delta_h
    #                 new_mol = Molecule(labels, coords, units='au')
    #                 new_mol.set_charge(charge)
    #                 new_mol.set_multiplicity(multiplicity)
    #                 current_string = f'{i}_{d}_plus'
    #                 # create a new XTB driver object;
    #                 # without this the energy is always zero...;
    #                 _, _, job_id_plus = gradient_driver.compute_gradient(new_mol, current_string)

    #                 coords[i, d] -= 2.0 * self.delta_h
    #                 new_mol = Molecule(labels, coords, units='au')
    #                 new_mol.set_charge(charge)
    #                 new_mol.set_multiplicity(multiplicity)
    #                 current_string = f'{i}_{d}_minus'

    #                 _, _, job_id_minus = gradient_driver.compute_gradient(new_mol, current_string)

    #                 # coords[i, d] += self.delta_h
    #                 # if self.rank == mpi_master():
    #                 #     for root in range(0, self.roots):
    #                 #         hessian[root, i, d, :, :] = (grad_plus[root] - grad_minus[root]) / (2.0 * self.delta_h)
    #                 job_ids.append((job_id_plus, job_id_minus))
    #         print('JOB IDS that needs to be checked', job_ids)
    #         i = 0
    #         while len(job_ids) != 0:
                
    #             job_is_running = self.cluster_manager.is_job_running(job_ids[i])
    #             if job_is_running == False:
    #                 job_ids[i] = 0
    #             if all(job_id == 0 for job_id in job_ids):
    #                 job_ids = []
    #             i += 1
    #             if i + 1 == len(job_ids):
    #                 i = 0

    #         exit()

        
    #     reshaped_hessian = np.zeros((self.roots, 3 * natm, 3 * natm))
    #     gradient_driver = ExternalQMGradientDriver(self.qm_driver)
    #     all_gradients = gradient_driver.extract_gradients(job_ids)

    #     # reshaped Hessian as member variable
    #     if self.rank == mpi_master():
    #         for root in range(0, self.roots):

    #             flat_hessian = hessian[root].flatten()

    #             # Ensure the total number of elements matches the desired shape
    #             expected_elements = (natm * 3)**2
    #             assert flat_hessian.size == expected_elements, "Total elements do not match the expected"
    #             reshaped_hessian[root] = flat_hessian.reshape(3 * natm, 3 * natm)
    #             reshaped_hessian[root] = ( reshaped_hessian[root] + reshaped_hessian[root].T ) / 2.0

    #         self.hessian = reshaped_hessian
    #     else:
    #         self.hessian = None

    #     self.ostream.unmute()
         
    #     print(self.hessian, self.hessian.shape, np.allclose(self.hessian[0], self.hessian[0].T))
    #     return reshaped_hessian

    
    def create_input_file_hessian(self, input_file, root):
        """
        Create the specific input file with the given content.
        """
        print('Giving of the current input file', input_file)
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

                current_root = root 
                full_path = os.path.abspath(self.xyz_filename)
                with open(input_file, 'w') as file:
                    basis_set = self.qm_driver.basis_set_label
                    if "RIJK" in basis_set:
                        basis_set = basis_set.replace("RIJK", "RIJCOSX")
                    if not self.qm_driver.libxc:
                        if self.qm_driver.solvation[0] is True:
                            file.write(f'!{self.qm_driver.xc_func} {self.qm_driver.dispersion} {basis_set} freq {self.qm_driver.solvation[1]}({self.qm_driver.solvation[2]})\n')
                        else:
                            file.write(f'!{self.qm_driver.xc_func} {self.qm_driver.dispersion} {basis_set} freq\n')
                    else:
                        if self.qm_driver.solvation[0] is True:
                            file.write(f'! {basis_set} freq {self.qm_driver.solvation[1]}({self.qm_driver.solvation[2]})\n')
                        else:
                            file.write(f'! {basis_set} freq\n')
                    file.write(f'%TDDFT NROOTS {self.qm_driver.roots_to_check}\n')
                    file.write(f'sf {self.qm_driver.spin_flip}\n')
                    file.write(f'IROOT {current_root}\n')
                    file.write('END\n')
                    if self.qm_driver.libxc:
                        file.write(f'%method \n')
                        file.write(f'method dft\n')
                        file.write(f'functional {self.qm_driver.xc_func}\n')
                        file.write('END\n')
                    file.write(f'%maxcore 3000\n')
                    file.write(f'%PAL\n')
                    file.write(f'nprocs {self.qm_driver.nprocs * 3}\n')
                    file.write('END\n')
                    file.write(f'* xyzfile {self.qm_driver.charge} {self.qm_driver.spin} {full_path}\n')
            
            elif self.program == 'OpenQP':
                '''[input]
                    system = opt.xyz
                    charge = 0
                    functional = cam-b3lyp
                    basis = def2-svp
                    method = tdhf
                    runtype = hess

                    [scf]
                    type = rohf
                    maxit = 300
                    multiplicity = 3

                    [tdhf]
                    type = mrsf
                    nstate = 10

                    [hess]
                    state = 2'''
                
                full_path = os.path.abspath(self.xyz_filename)
                root_str = ', '.join(map(str, [root]))
                with open(input_file, "w") as f:
                    # [input]
                    f.write("[input]\n")
                    # Many tools want just the filename; if OpenQP expects a path, we pass the resolved path.
                    # If you prefer only the basename, switch to os.path.basename(system_path).
                    f.write(f"system={full_path}\n")
                    f.write(f"charge={self.qm_driver.charge}\n")
                    f.write(f"runtype=hess\n")
                    f.write(f"basis={self.qm_driver.basis_set_label}\n")
                    f.write(f"functional={self.qm_driver.xc_func}\n")
                    f.write(f"method={self.qm_driver.method}\n")
                    f.write("\n")

                    # [guess]
                    f.write("[guess]\n")
                    f.write(f"type=huckel\n")
                    f.write("\n")

                    # [scf]
                    f.write("[scf]\n")
                    f.write(f"multiplicity={self.qm_driver.spin}\n")
                    f.write(f"type=rohf\n")
                    f.write("\n")

                    # Method-specific block: use the method name as the section header (e.g., [tdhf] or [tddft])
                    f.write(f"[tdhf]\n")
                    f.write(f"type=mrsf\n")
                    f.write(f"nstate={root_str}\n")

                    f.write(f"[hess]\n")
                    f.write(f"state={root_str}\n")
            
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
                    if root == 0 and self.qm_driver.spin_flip is not True and self.roots == 0:

                        file.write('$rem\n')
                        file.write('JOBTYPE         freq\n')
                        file.write('EXCHANGE        b3lyp\n')
                        file.write('BASIS           6-31G*\n')
                        file.write('SCF_GUESS       core\n')
                        file.write('$end\n')
                    else:
                        file.write('$rem\n')
                        file.write('JOBTYPE          freq\n')
                        file.write('EXCHANGE         bhhlyp\n')
                        file.write('BASIS            6-311G(d,p)\n')
                        file.write('UNRESTRICTED     true\n')
#                        file.write('SCF_GUESS        core\n')
                        file.write(f'SPIN_FLIP       {self.qm_driver.spin_flip}\n')
                        file.write(f'CIS_N_ROOTS     {root}\n')
                        file.write(f'CIS_STATE_DERIV {root}\n')
                        file.write('CIS_CONVERGENCE 0\n')
                        file.write(f'$end\n')
                        #file.write('\n\n@@@\n\n')

                    if root + 1 == roots and self.qm_driver.spin_flip is True and self.NAC:
                        sequence = " ".join(str(i) for i in range(1, roots + 1))

                        file.write('$molecule\n')
                        file.write(f'0 {state}\n')
                        with open(self.xyz_filename, 'r') as xyz_file:
                            # Skip the first two lines (atom count and comment)
                            next(xyz_file)
                            next(xyz_file)
                            # Write the remaining lines (atom coordinates)
                            for line in xyz_file:
                                file.write(line)
                        file.write('$end\n\n')

                        file.write('$rem\n')
                        file.write('EXCHANGE        bhhlyp\n')
                        file.write('BASIS           6-31G*\n')
                        file.write('SCF_GUESS       core\n')
                        file.write(f'CALC_NAC       {self.NAC}\n')
                        file.write(f'SPIN_FLIP      {self.qm_driver.spin_flip}\n')
                        file.write(f'CIS_DER_NUMSTATE   {self.roots + 1}\n')
                        file.write('CIS_CONVERGENCE 0\n')
                        file.write(f'$end\n')
                        file.write('\n$derivative_coupling')
                        file.write('comment')
                        file.write(f'{sequence}')
                        file.write(f'$end\n')

                    elif root == roots and self.qm_driver.spin_flip is True and self.NAC:
                        sequence = " ".join(str(i) for i in range(1, roots + 1))

                        file.write('$molecule\n')
                        file.write(f'0 {state}\n')
                        with open(self.xyz_filename, 'r') as xyz_file:
                            # Skip the first two lines (atom count and comment)
                            next(xyz_file)
                            next(xyz_file)
                            # Write the remaining lines (atom coordinates)
                            for line in xyz_file:
                                file.write(line)
                        file.write('$end\n\n')
                        file.write('$rem\n')
                        file.write('EXCHANGE        bhhlyp\n')
                        file.write('BASIS           6-31G*\n')
                        file.write('SCF_GUESS       core\n')
                        file.write(f'CALC_NAC       {self.NAC}\n')
                        file.write(f'SPIN_FLIP      {self.qm_driver.spin_flip}\n')
                        file.write(f'CIS_DER_NUMSTATE   {self.roots + 1}\n')
                        file.write('CIS_CONVERGENCE 0\n')
                        file.write(f'$end\n')
                        file.write('\n$derivative_coupling\n')
                        file.write('comment\n')
                        file.write(f'{sequence}\n')
                        file.write(f'$end\n')

            print(f"Input file '{input_file}' created successfully.")
        except Exception as e:
            print(f"Error creating input file: {e}")
