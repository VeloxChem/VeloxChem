#
#                           VELOXCHEM 1.0-RC2
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

from .cppsolver import ComplexResponse
from .lrsolver import LinearResponseSolver
from .lreigensolver import LinearResponseEigenSolver
from .c6solver import C6Solver
from .tdaexcidriver import TDAExciDriver
from .tpafulldriver import TpaFullDriver
from .tpareddriver import TpaReducedDriver
from .shgdriver import SHGDriver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input


class ResponseDriver:
    """
    Implements response driver for molecular property calculations using
    conventional Hartree-Fock/Kohn-Sham response theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - prop_type: The type of the property to be calculated.
        - tamm_dancoff: The flag for using Tamm-Dancoff approximation.
        - rsp_dict: The dictionary of response input.
        - method_dict: The dictionary of method settings.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
        - is_converged: The flag for convergence.
    """

    def __init__(self, comm, ostream):
        """
        Initializes response driver to default setup.
        """

        # default calculation type
        self.prop_type = 'generic'
        self.tamm_dancoff = False
        self.rsp_dict = {}
        self.method_dict = {}

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # solver and convergence flag
        self.solver = None
        self.is_converged = False

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates settings in response solver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            self.method_dict = None
        else:
            self.method_dict = dict(method_dict)

        self.rsp_dict = dict(rsp_dict)
        self.rsp_dict['prop_type'] = self.rsp_dict['property']

        rsp_keywords = {
            'prop_type': 'str_lower',
            'tamm_dancoff': 'bool',
        }

        parse_input(self, rsp_keywords, self.rsp_dict)

        # Linear response eigensolver
        if (self.rsp_dict['order'] == 'linear' and
                self.rsp_dict['residue'] == 'single' and
                self.rsp_dict['complex'] == 'no'):
            if self.tamm_dancoff:
                self.solver = TDAExciDriver(self.comm, self.ostream)
            else:
                self.solver = LinearResponseEigenSolver(self.comm, self.ostream)
            self.solver.input_keywords['response'].update({
                'tamm_dancoff': ('bool', 'use Tamm-Dancoff approximation'),
            })

        # Linear response solver
        elif (self.rsp_dict['order'] == 'linear' and
              self.rsp_dict['residue'] == 'none' and
              self.rsp_dict['complex'] == 'no'):
            self.solver = LinearResponseSolver(self.comm, self.ostream)

        # Complex linear response solver
        elif (self.rsp_dict['order'] == 'linear' and
              self.rsp_dict['residue'] == 'none' and
              self.rsp_dict['onlystatic'] == 'no' and
              self.rsp_dict['complex'] == 'yes'):
            self.solver = ComplexResponse(self.comm, self.ostream)

        # C6 linear response solver
        elif (self.rsp_dict['order'] == 'linear' and
              self.rsp_dict['residue'] == 'none' and
              self.rsp_dict['onlystatic'] == 'yes' and
              self.rsp_dict['complex'] == 'yes'):
            self.solver = C6Solver(self.comm, self.ostream)

        # SHG
        if (self.rsp_dict['order'] == 'quadratic' and
                self.rsp_dict['complex'] == 'yes'):
            self.solver = SHGDriver(self.comm, self.ostream)

        # TPA
        elif (self.rsp_dict['order'] == 'cubic' and
              self.rsp_dict['complex'] == 'yes'):
            if ('tpa_type' not in self.rsp_dict or
                    self.rsp_dict['tpa_type'].lower() == 'full'):
                self.solver = TpaFullDriver(self.comm, self.ostream)
            elif ('tpa_type' in self.rsp_dict and
                  self.rsp_dict['tpa_type'].lower() == 'reduced'):
                self.solver = TpaReducedDriver(self.comm, self.ostream)
            self.solver.input_keywords['response'].update({
                'tpa_type': ('str_lower', 'full or reduced TPA calculation'),
            })

        self.solver.update_settings(self.rsp_dict, self.method_dict)

    def compute(self, molecule, ao_basis, scf_tensors):
        """
        Performs molecular property calculation using molecular data

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The results from the actual response solver.
        """

        assert_msg_critical(self.solver is not None,
                            'ResponseDriver: solver not initialized')

        result = self.solver.compute(molecule, ao_basis, scf_tensors)

        self.is_converged = self.solver.is_converged

        return result

    def prop_str(self):
        """
        Gets string with type of molecular property calculation (Excited
        states, linear and non-linear spectroscopies).

        :return:
            The string with type of molecular property calculation.
        """

        if self.prop_type == 'polarizability':
            return 'Polarizability'

        if self.prop_type == 'absorption':
            return 'Singlet Excited States'

        if self.prop_type == 'linear absorption cross-section':
            return 'Linear Absorption Cross-Section'

        if self.prop_type == 'circular dichroism spectrum':
            return 'Circular Dichroism Spectrum'

        if self.prop_type == 'c6':
            return 'C6 Dispersion Coefficient'

        if self.prop_type == 'custom':
            return 'Custom Response Property'

        return 'Undefined'
