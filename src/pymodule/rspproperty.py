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
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .cppsolver import ComplexResponse
from .tdacppsolver import ComplexResponseTDA
from .lrsolver import LinearResponseSolver
from .lreigensolver import LinearResponseEigenSolver
from .c6driver import C6Driver
from .tdaeigensolver import TdaEigenSolver
from .shgdriver import ShgDriver
from .tpatransitiondriver import TpaTransitionDriver
from .threepatransitiondriver import ThreePATransitionDriver
from .tpafulldriver import TpaFullDriver
from .tpareddriver import TpaReducedDriver
from .quadraticresponsedriver import QuadraticResponseDriver
from .cubicresponsedriver import CubicResponseDriver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input


class ResponseProperty:
    """
    Implements the base class for response property/spectroscopy.

    :param rsp_dict:
        The input dictionary that defines the property/spectroscopy.
    :param method_dict:
        The dictionary of method settings.

    Instance variables
        - rsp_dict: The dictionary of response input.
        - method_dict: The dictionary of method settings.
        - rsp_driver: The response driver.
        - rsp_property: The dictionary of response property.
    """

    def __init__(self, rsp_dict=None, method_dict=None):
        """
        Initializes response property/spectroscopy.
        """

        self._rsp_dict = rsp_dict
        self._method_dict = method_dict

        self._rsp_driver = None
        self._rsp_property = None

    def init_driver(self, comm=None, ostream=None):
        """
        Initializes response driver.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.prop_type = 'generic'
        self.tamm_dancoff = False

        rsp_keywords = {
            'prop_type': 'str_lower',
            'tamm_dancoff': 'bool',
        }
        self._rsp_dict.update({
            'prop_type': self._rsp_dict['property'],
        })
        parse_input(self, rsp_keywords, self._rsp_dict)

        self._rsp_driver = None
        self._is_converged = False

        # Custom linear response
        if (self.prop_type == 'custom' and
                self._rsp_dict['order'] == 'linear' and
                self._rsp_dict['residue'] == 'none' and
                self._rsp_dict['onlystatic'] == 'no'):

            if self._rsp_dict['is_complex'] == 'no':
                self._rsp_driver = LinearResponseSolver(self.comm, self.ostream)

            elif self._rsp_dict['is_complex'] == 'yes':
                self._rsp_driver = ComplexResponse(self.comm, self.ostream)

        # Linear response real solver
        elif (self._rsp_dict['order'] == 'linear' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['is_complex'] == 'no'):

            self._rsp_driver = LinearResponseSolver(self.comm, self.ostream)

        # Linear response complex solver
        elif (self._rsp_dict['order'] == 'linear' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['onlystatic'] == 'no' and
              self._rsp_dict['is_complex'] == 'yes'):

            if self.tamm_dancoff:
                self._rsp_driver = ComplexResponseTDA(self.comm, self.ostream)
            else:
                self._rsp_driver = ComplexResponse(self.comm, self.ostream)

            self._rsp_driver._input_keywords['response'].update({
                'tamm_dancoff': ('bool', 'use Tamm-Dancoff approximation'),
            })

            if self.prop_type in [
                    'linear absorption cross-section',
                    'linear absorption (cpp)',
                    'linear absorption(cpp)',
                    'absorption (cpp)',
                    'absorption(cpp)',
            ]:
                self._rsp_driver.set_cpp_flag('absorption')

            elif self.prop_type in [
                    'circular dichroism spectrum',
                    'circular dichroism (cpp)',
                    'circular dichroism(cpp)',
                    'ecd (cpp)',
                    'ecd(cpp)',
            ]:
                self._rsp_driver.set_cpp_flag('ecd')

        # Linear response C6 solver
        elif (self._rsp_dict['order'] == 'linear' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['onlystatic'] == 'yes' and
              self._rsp_dict['is_complex'] == 'yes'):

            self._rsp_driver = C6Driver(self.comm, self.ostream)

        # Linear response eigensolver (RPA/TDA)
        elif (self._rsp_dict['order'] == 'linear' and
              self._rsp_dict['residue'] == 'single' and
              self._rsp_dict['is_complex'] == 'no'):

            if self.tamm_dancoff:
                self._rsp_driver = TdaEigenSolver(self.comm, self.ostream)
            else:
                self._rsp_driver = LinearResponseEigenSolver(
                    self.comm, self.ostream)

            self._rsp_driver._input_keywords['response'].update({
                'tamm_dancoff': ('bool', 'use Tamm-Dancoff approximation'),
            })

        # Quadratic response driver
        elif (self.prop_type == 'custom' and
              self._rsp_dict['order'] == 'quadratic' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['is_complex'] == 'yes'):

            self._rsp_driver = QuadraticResponseDriver(self.comm, self.ostream)

        # SHG (quadratic response) driver
        elif (self._rsp_dict['order'] == 'quadratic' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['is_complex'] == 'yes'):

            self._rsp_driver = ShgDriver(self.comm, self.ostream)

        # TPA transition (quadratic response) driver
        elif (self._rsp_dict['order'] == 'quadratic' and
              self._rsp_dict['residue'] == 'single' and
              self._rsp_dict['is_complex'] == 'yes'):

            self._rsp_driver = TpaTransitionDriver(self.comm, self.ostream)

        # Cubic response driver
        elif (self.prop_type == 'custom' and
              self._rsp_dict['order'] == 'cubic' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['is_complex'] == 'yes'):

            self._rsp_driver = CubicResponseDriver(self.comm, self.ostream)

        # TPA (cubic response) driver
        elif (self._rsp_dict['order'] == 'cubic' and
              self._rsp_dict['residue'] == 'none' and
              self._rsp_dict['is_complex'] == 'yes'):

            if ('tpa_type' not in self._rsp_dict or
                    self._rsp_dict['tpa_type'].lower() == 'full'):
                self._rsp_driver = TpaFullDriver(self.comm, self.ostream)

            elif ('tpa_type' in self._rsp_dict and
                  self._rsp_dict['tpa_type'].lower() == 'reduced'):
                self._rsp_driver = TpaReducedDriver(self.comm, self.ostream)

            self._rsp_driver._input_keywords['response'].update({
                'tpa_type': ('str_lower', 'full or reduced TPA calculation'),
            })

        # 3PA (cubic response) driver
        elif (self._rsp_dict['order'] == 'cubic' and
              self._rsp_dict['residue'] == 'single'):

            self._rsp_driver = ThreePATransitionDriver(self.comm, self.ostream)

        # Update driver settings
        self._rsp_driver.update_settings(self._rsp_dict, self._method_dict)

    def print_keywords(self):
        """
        Prints input keywords for response property.
        """

        assert_msg_critical(
            self._rsp_driver is not None,
            'ResponseProperty: response driver not initialized')

        self._rsp_driver.print_keywords()

    def compute(self, molecule, basis, scf_tensors):
        """
        Computes response property/spectroscopy.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        """

        self._rsp_property = self._rsp_driver.compute(molecule, basis,
                                                      scf_tensors)

    @property
    def rsp_driver(self):
        """
        Returns the response driver.
        """

        return self._rsp_driver

    @property
    def is_converged(self):
        """
        Returns whether the response calculation is converged.
        """

        return self._rsp_driver.is_converged

    @property
    def rsp_property(self):
        """
        Returns the response property dictionary.
        """

        return self._rsp_property

    def get_property(self, key):
        """
        Gets response property/spectroscopy.

        :param key:
            The keyword for the property.

        :return:
            The property.
        """

        return self._rsp_property[key]

    def get_full_solution_vector(self, key):
        """
        Gets response solution vector for a given key.

        :param key:
            The key for the solution vector.

        :return:
            The solution vector.
        """

        return self._rsp_driver.get_full_solution_vector(
            self._rsp_property['solutions'][key])
