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

from mpi4py import MPI
import sys

from .veloxchemlib import mpi_master
from .rspdriver import ResponseDriver
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


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

        self.rsp_dict = rsp_dict
        self.method_dict = method_dict

        self.rsp_driver = None

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

        self.rsp_driver = ResponseDriver(comm, ostream)
        self.rsp_driver.update_settings(self.rsp_dict, self.method_dict)

    def print_keywords(self):
        """
        Prints input keywords for response property.
        """

        assert_msg_critical(
            self.rsp_driver is not None,
            'ResponseProperty: response driver not initialized')

        assert_msg_critical(
            self.rsp_driver.solver is not None,
            'ResponseProperty: response solver not initialized')

        self.rsp_driver.solver.print_keywords()

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

        self.rsp_property = self.rsp_driver.compute(molecule, basis,
                                                    scf_tensors)

        if not self.rsp_driver.is_converged:
            return

        if self.rsp_driver.rank == mpi_master():
            self.print_property(self.rsp_driver.ostream)

    def converged(self):
        """
        Checks if the response calculation is converged.

        :return:
            True if the response calculation is converged, False otherwise.
        """

        return self.rsp_driver.is_converged

    def get_property(self, key):
        """
        Gets response property/spectroscopy.

        :param key:
            The keyword for the property.

        :return:
            The property.
        """

        return None

    def print_property(self, ostream):
        """
        Prints response property/spectroscopy to output stream.

        :param ostream:
            The output stream.
        """

        return
