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

import numpy as np
import sys

try:
    from dftd4.interface import DampingParam as D4Param
    from dftd4.interface import DispersionModel as D4Model
except ImportError:
    pass


class DispersionModel:
    """
    Uses the dftd4-python interface.

    Instance variables
        - _energy: The dispersion correction to energy.
        - _gradient: The dispersion correction to gradient.
    """

    def __init__(self):
        """
        Initializes the model.
        """

        self._energy = None
        self._gradient = None

    @staticmethod
    def is_available():
        """
        Checks if dftd4-python is available.

        :return:
            True if dftd4-python is available, False otherwise.
        """

        return ('dftd4' in sys.modules)

    def compute(self, molecule, xc_label):
        """
        Uses the dftd4-python interface to compute dispersion correction.

        :param molecule:
            The molecule.
        :param xc_label:
            The label of XC functional.
        """

        identifiers_np = np.array(molecule.get_identifiers())
        coords_in_au = molecule.get_coordinates_in_bohr()
        net_charge = molecule.get_charge()

        disp_model = D4Model(numbers=identifiers_np,
                             positions=coords_in_au,
                             charge=net_charge,
                             model='d4')

        disp_res = disp_model.get_dispersion(D4Param(method=xc_label),
                                             grad=True)

        self._energy = disp_res.get("energy")
        self._gradient = disp_res.get("gradient")

    def get_energy(self):
        """
        Gets dispersion correction to energy.

        :return:
            The dispersion correction to energy.
        """

        return self._energy

    def get_gradient(self):
        """
        Gets dispersion correction to gradient.

        :return:
            The dispersion correction to gradient.
        """

        return self._gradient
