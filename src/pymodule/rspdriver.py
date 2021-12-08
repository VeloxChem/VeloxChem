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

        # convergence
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

        # Linear response eigensolver

        if (self.rsp_dict['response'] == 'linear' and
                self.rsp_dict['residue'] == 'single' and
                self.rsp_dict['complex'] == 'no'):

            if self.tamm_dancoff:
                solver = TDAExciDriver(self.comm, self.ostream)
            else:
                solver = LinearResponseEigenSolver(self.comm, self.ostream)

            solver.update_settings(self.rsp_dict, self.method_dict)

            result = solver.compute(molecule, ao_basis, scf_tensors)

            self.is_converged = solver.is_converged

            return result

        # Linear response solver

        if (self.rsp_dict['response'] == 'linear' and
                self.rsp_dict['residue'] == 'none' and
                self.rsp_dict['complex'] == 'no'):

            solver = LinearResponseSolver(self.comm, self.ostream)

            solver.update_settings(self.rsp_dict, self.method_dict)

            result = solver.compute(molecule, ao_basis, scf_tensors)

            self.is_converged = solver.is_converged

            return result

        # Complex linear response solver

        if (self.rsp_dict['response'] == 'linear' and
                self.rsp_dict['residue'] == 'none' and
                self.rsp_dict['onlystatic'] == 'no' and
                self.rsp_dict['complex'] == 'yes'):

            clr_solver = ComplexResponse(self.comm, self.ostream)

            clr_solver.update_settings(self.rsp_dict, self.method_dict)

            clr_result = clr_solver.compute(molecule, ao_basis, scf_tensors)

            self.is_converged = clr_solver.is_converged

            return clr_result

        # C6 linear response solver

        if (self.rsp_dict['response'] == 'linear' and
                self.rsp_dict['residue'] == 'none' and
                self.rsp_dict['onlystatic'] == 'yes' and
                self.rsp_dict['complex'] == 'yes'):

            c6_solver = C6Solver(self.comm, self.ostream)

            c6_solver.update_settings(self.rsp_dict, self.method_dict)

            c6_result = c6_solver.compute(molecule, ao_basis, scf_tensors)

            self.is_converged = c6_solver.is_converged

            return c6_result
            
        # SHG

        if (self.rsp_dict['response'] == 'quadratic' and
                self.rsp_dict['complex'] == 'yes'):

            shg_solver = SHGDriver(self.comm, self.ostream)

            shg_solver.update_settings(self.rsp_dict, self.method_dict)

            shg_results = shg_solver.compute(molecule, ao_basis, scf_tensors,self.method_dict)

            self.is_converged = shg_solver.is_converged

            return shg_results


        # TPA

        if (self.rsp_dict['response'] == 'cubic' and
                self.rsp_dict['complex'] == 'yes'):

            if ('tpa_type' not in self.rsp_dict or
                    self.rsp_dict['tpa_type'].lower() == 'full'):
                tpa_solver = TpaFullDriver(self.comm, self.ostream)

            elif ('tpa_type' in self.rsp_dict and
                  self.rsp_dict['tpa_type'].lower() == 'reduced'):
                tpa_solver = TpaReducedDriver(self.comm, self.ostream)

            tpa_solver.update_settings(self.rsp_dict, self.method_dict)

            tpa_result = tpa_solver.compute(molecule, ao_basis, scf_tensors)

            self.is_converged = tpa_solver.is_converged

            return tpa_result

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
