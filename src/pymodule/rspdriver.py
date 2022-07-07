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
from .tpafulldriver import TPAFullDriver
from .tpareddriver import TPAReducedDriver
from .shgdriver import SHGDriver
from .quadraticresponsedriver import QuadraticResponseDriver
from .cubicresponsedriver import CubicResponseDriver
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
        self._rsp_dict = {}
        self._method_dict = {}

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # solver and convergence flag
        self._solver = None
        self._is_converged = False

    @property
    def is_converged(self):
        """
        Returns whether response driver is converged.
        """

        return self._is_converged

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates settings in response solver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            self._method_dict = None
        else:
            self._method_dict = dict(method_dict)

        self._rsp_dict = dict(rsp_dict)
        self._rsp_dict['prop_type'] = self._rsp_dict['property']

        rsp_keywords = {
            'prop_type': 'str_lower',
            'tamm_dancoff': 'bool',
        }

        parse_input(self, rsp_keywords, self._rsp_dict)

        # Linear response eigensolver
        if (self._rsp_dict['order'] == 'linear' and
                self._rsp_dict['residue'] == 'single' and
                self._rsp_dict['complex'] == 'no'):
            if self.tamm_dancoff:
                self._solver = TDAExciDriver(self.comm, self.ostream)
            else:
                self._solver = LinearResponseEigenSolver(
                    self.comm, self.ostream)
            self._solver._input_keywords['response'].update({
                'tamm_dancoff': ('bool', 'use Tamm-Dancoff approximation'),
            })

        # Linear response solver
        elif (self._rsp_dict['order'] == 'linear' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['complex'] == 'no'):
            self._solver = LinearResponseSolver(self.comm, self.ostream)

        # Complex linear response solver
        elif (self._rsp_dict['order'] == 'linear' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['onlystatic'] == 'no' and
              self._rsp_dict['complex'] == 'yes'):
            self._solver = ComplexResponse(self.comm, self.ostream)

        # C6 linear response solver
        elif (self._rsp_dict['order'] == 'linear' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['onlystatic'] == 'yes' and
              self._rsp_dict['complex'] == 'yes'):
            self._solver = C6Solver(self.comm, self.ostream)

        # SHG
        if (self._rsp_dict['order'] == 'quadratic' and
                self._rsp_dict['complex'] == 'yes'):
            self._solver = SHGDriver(self.comm, self.ostream)

        # TPA
        elif (self._rsp_dict['order'] == 'cubic' and
              self._rsp_dict['complex'] == 'yes'):
            if ('tpa_type' not in self._rsp_dict or
                    self._rsp_dict['tpa_type'].lower() == 'full'):
                self._solver = TPAFullDriver(self.comm, self.ostream)
            elif ('tpa_type' in self._rsp_dict and
                  self._rsp_dict['tpa_type'].lower() == 'reduced'):
                self._solver = TPAReducedDriver(self.comm, self.ostream)
            self._solver._input_keywords['response'].update({
                'tpa_type': ('str_lower', 'full or reduced TPA calculation'),
            })

        # Quadratic response driver
        if (self.prop_type == 'custom' and
                self._rsp_dict['order'] == 'quadratic' and
                self._rsp_dict['complex'] == 'yes'):
            self._solver = QuadraticResponseDriver(self.comm, self.ostream)

        # Cubic response driver
        if (self.prop_type == 'custom' and
                self._rsp_dict['order'] == 'cubic' and
                self._rsp_dict['complex'] == 'yes'):
            self._solver = CubicResponseDriver(self.comm, self.ostream)

        self._solver.update_settings(self._rsp_dict, self._method_dict)

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

        assert_msg_critical(self._solver is not None,
                            'ResponseDriver: solver not initialized')

        result = self._solver.compute(molecule, ao_basis, scf_tensors)

        self._is_converged = self._solver.is_converged

        return result

    def get_prop_str(self):
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
