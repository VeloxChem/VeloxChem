#
#                           VELOXCHEM 1.0-RC
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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
                if isinstance(fname_list[1], str):
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

            for key in self.input_dict:
                if key == 'molecule':
                    continue
                if not isinstance(self.input_dict[key], dict):
                    continue
                if list(self.input_dict[key].keys()) == ['checkpoint_file']:
                    continue
                self.ostream.print_info('@{:s}'.format(key))
                for key_2 in self.input_dict[key]:
                    if isinstance(self.input_dict[key][key_2], str):
                        self.ostream.print_info('{:s}: {:s}'.format(
                            key_2, self.input_dict[key][key_2]))
                    elif isinstance(self.input_dict[key][key_2], list):
                        self.ostream.print_info('{:s}:'.format(key_2))
                        for line in self.input_dict[key][key_2]:
                            self.ostream.print_info(line)
                self.ostream.print_info('@end')
                self.ostream.print_blank()

            # create molecule

            self.ostream.print_info('Parsing @molecule group...')
            self.ostream.print_info('...done.')
            self.ostream.print_blank()

            self.molecule = Molecule.from_dict(self.input_dict['molecule'])

            self.ostream.print_block(self.molecule.get_string())
            self.ostream.print_block(self.molecule.more_info())
            self.ostream.print_blank()

            # create basis set

            self.ostream.print_info('Parsing @method settings group...')
            self.ostream.print_info('...done.')
            self.ostream.print_blank()

            if 'xtb' not in self.input_dict['method_settings']:
                basis_path = '.'
                if 'basis_path' in self.input_dict['method_settings']:
                    basis_path = self.input_dict['method_settings'][
                        'basis_path']

                basis_name = self.input_dict['method_settings']['basis'].upper()

                self.ao_basis = MolecularBasis.read(self.molecule, basis_name,
                                                    basis_path, self.ostream)

                self.min_basis = MolecularBasis.read(self.molecule,
                                                     'MIN-CC-PVDZ', basis_path)

                self.ostream.print_block(
                    self.ao_basis.get_string('Atomic Basis', self.molecule))

                self.ostream.flush()

        # broadcast input dictionary

        self.input_dict = self.mpi_comm.bcast(self.input_dict,
                                              root=mpi_master())

        # broadcast molecule and basis set

        self.molecule.broadcast(self.mpi_rank, self.mpi_comm)
        self.ao_basis.broadcast(self.mpi_rank, self.mpi_comm)
        self.min_basis.broadcast(self.mpi_rank, self.mpi_comm)

    def finish(self):
        """
        Finalizes the MPI task.
        """

        if self.mpi_rank == mpi_master():
            self.ostream.print_finish_header(self.start_time)
