#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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
import numpy as np
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

try:
    from xtb.interface import Calculator as XtbCalculator
    from xtb.utils import get_method as xtb_get_method
except ImportError:
    pass


class XtbDriver:
    """
    Implements XTB driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - flag: The driver flag.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes XTB driver.
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

        self.xtb_method = 'gfn2'

        self._xtb_calc = None
        self._xtb_res = None

    @staticmethod
    def is_available():
        """
        Returns if XTB driver is available.
        """

        return ('xtb' in sys.modules)

    def set_method(self, xtb_method):
        """
        Sets XTB method.
        """

        if self.rank == mpi_master():
            self.xtb_method = xtb_method

    def get_method(self):
        """
        Gets XTB method.
        """

        if self.rank == mpi_master():
            return self.xtb_method
        else:
            return None

    def get_energy(self):
        """
        Gets XTB energy.
        """

        if self.rank == mpi_master() and self._xtb_res is not None:
            return self._xtb_res.get_energy()
        else:
            return None

    def get_gradient(self):
        """
        Gets XTB gradient.
        """

        if self.rank == mpi_master() and self._xtb_res is not None:
            return self._xtb_res.get_gradient()
        else:
            return None

    def get_dipole(self):
        """
        Gets XTB dipole.
        """

        if self.rank == mpi_master() and self._xtb_res is not None:
            return self._xtb_res.get_dipole()
        else:
            return None

    def compute(self, molecule):
        """
        Performs XTB calculation.

        :param molecule:
            The molecule.

        :return:
            The results from XTB calculation.
        """

        # sanity check

        errmsg = 'XtbDriver: xtb-python not available. Please install xtb-python.'

        assert_msg_critical(self.is_available(), errmsg)

        identifiers = np.array(molecule.get_identifiers())
        coords_in_au = molecule.get_coordinates_in_bohr()

        self._xtb_calc = XtbCalculator(
            xtb_get_method(self.xtb_method.upper() + '-xTB'), identifiers,
            coords_in_au)

        if self.rank == mpi_master():

            # set verbosity

            if self.ostream.is_muted:
                self._xtb_calc.set_verbosity('muted')
            else:
                self._xtb_calc.set_verbosity('full')

            # run XTB calculation

            self.print_title()

            self._xtb_res = self._xtb_calc.singlepoint()

            # process results

            energy = self._xtb_res.get_energy()
            gradient = self._xtb_res.get_gradient()
            dipole = self._xtb_res.get_dipole()
            partial_charges = self._xtb_res.get_charges()
            bond_orders = self._xtb_res.get_bond_orders()
            orbital_energies = self._xtb_res.get_orbital_eigenvalues()
            orbital_occupations = self._xtb_res.get_orbital_occupations()

            grad2 = np.sum(gradient**2, axis=1)
            rms_grad = np.sqrt(np.mean(grad2))
            max_grad = np.max(np.sqrt(grad2))

            xtb_results = {
                'energy': energy,
                'gradient': gradient,
                'max_gradient': max_grad,
                'rms_gradient': rms_grad,
                'dipole': dipole,
                'partial_charges': partial_charges,
                'bond_orders': bond_orders,
                'orbital_energies': orbital_energies,
                'orbital_occupations': orbital_occupations,
            }

            return xtb_results

        else:
            return None

    def print_title(self):
        """
        Prints title for XTB calculation.
        """

        self.ostream.print_blank()
        self.ostream.print_header('XTB Driver')
        self.ostream.print_header(12 * '=')
        self.ostream.print_blank()

        self.ostream.print_reference('Reference:')
        self.ostream.print_reference(self.get_reference())
        self.ostream.flush()

    def get_reference(self):
        """
        Gets reference string for XTB.
        """

        ref_str = 'C. Bannwarth, E. Caldeweyher, S. Ehlert, '
        ref_str += 'A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme, '
        ref_str += 'WIREs Comput. Mol. Sci., 2020, 11, e01493'

        return ref_str
