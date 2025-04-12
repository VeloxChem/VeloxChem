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
import sys

from .veloxchemlib import mpi_master
from .inputparser import InputParser
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .errorhandler import assert_msg_critical


class MpiTask:
    """
    Implements the MPI task.

    :param fname_list:
        List of the input/output filenames.
    :param mpi_comm:
        The MPI communicator.

    Instance variable
        - mpi_comm: The MPI communicator.
        - mpi_rank: The MPI rank.
        - mpi_size: Number of MPI processes.
        - molecule: The molecule.
        - ao_basis: The AO basis set.
        - min_basis: The minimal AO basis set for generating initial guess.
        - ostream: The output stream.
        - input_dict: The input dictionary.
        - start_time: The start time of the task.
    """

    def __init__(self, fname_list, mpi_comm=None):
        """
        Initializes the MPI task.
        """

        # mpi settings
        if mpi_comm is None:
            self.mpi_comm = MPI.COMM_WORLD
        else:
            self.mpi_comm = mpi_comm
        self.mpi_rank = self.mpi_comm.Get_rank()
        self.mpi_size = self.mpi_comm.Get_size()

        input_fname = None
        output_fname = None

        if self.mpi_rank == mpi_master():
            assert_msg_critical(
                len(fname_list) > 0, 'MpiTask: Need input file name')

            input_fname = fname_list[0]
            output_fname = sys.stdout

            if len(fname_list) > 1:
                if fname_list[1] is None:
                    output_fname = None
                elif isinstance(fname_list[1], str):
                    if fname_list[1].strip().split():
                        output_fname = fname_list[1].strip()
                        if ('\0' in output_fname or output_fname == '-' or
                                output_fname == '.'):
                            output_fname = sys.stdout

            assert_msg_critical(
                Path(input_fname).is_file(),
                f'MpiTask: input file {input_fname} does not exist')

            assert_msg_critical(
                input_fname != output_fname,
                'MpiTask: input/output file cannot be the same')

        # initialize molecule, basis set and output stream

        self.molecule = Molecule()
        self.ao_basis = MolecularBasis()
        self.min_basis = MolecularBasis()
        self.ostream = OutputStream(output_fname)

        # process input file on master node

        self.input_dict = {}

        if self.mpi_rank == mpi_master():

            self.start_time = self.ostream.print_start_header(self.mpi_size)

            self.ostream.print_info(
                'Reading input file {}...'.format(input_fname))

            # read input file

            self.input_dict = InputParser(input_fname, output_fname).get_dict()
            self.ostream.print_blank()

            for group_key, group_val in self.input_dict.items():
                if group_key == 'molecule':
                    continue

                if not group_val:
                    continue

                if isinstance(group_val, dict):
                    self.ostream.print_info(f'@{group_key}')
                    for param_key, param_val in group_val.items():
                        if isinstance(param_val, str):
                            self.ostream.print_info(f'{param_key}: {param_val}')
                        elif isinstance(param_val, list):
                            self.ostream.print_info(f'{param_key}:')
                            for line in param_val:
                                self.ostream.print_info(line)
                    self.ostream.print_info('@end')
                    self.ostream.print_blank()

                elif isinstance(group_val, list):
                    self.ostream.print_info(f'@{group_key}')
                    for line in group_val:
                        self.ostream.print_info(line)
                    self.ostream.print_info('@end')
                    self.ostream.print_blank()

        # broadcast input dictionary

        self.input_dict = self.mpi_comm.bcast(self.input_dict,
                                              root=mpi_master())

        # create molecule

        if ('molecule' in self.input_dict and
            ('xyz' in self.input_dict['molecule'] or
             'xyzfile' in self.input_dict['molecule'])):

            # TODO: read xyz on master node and then broadcast xyz_string
            self.molecule = Molecule.from_dict(self.input_dict['molecule'])

            self.ostream.print_block(self.molecule.get_string())
            self.ostream.print_block(self.molecule.more_info())
            self.ostream.print_blank()

        # create basis set

        if 'xtb' not in self.input_dict['method_settings']:
            basis_path = '.'
            if 'basis_path' in self.input_dict['method_settings']:
                basis_path = self.input_dict['method_settings']['basis_path']

            basis_name = self.input_dict['method_settings']['basis'].upper()

            if ('molecule' in self.input_dict and
                ('xyz' in self.input_dict['molecule'] or
                 'xyzfile' in self.input_dict['molecule'])):

                # TODO: read basis on master node and then broadcast basis_dict
                self.ao_basis = MolecularBasis.read(self.molecule, basis_name,
                                                    basis_path, self.ostream)

                self.min_basis = MolecularBasis.read(self.molecule,
                                                     'AO-START-GUESS',
                                                     basis_path,
                                                     ostream=None)

    def finish(self):
        """
        Finalizes the MPI task.
        """

        if self.mpi_rank == mpi_master():
            self.ostream.print_finish_header(self.start_time)
