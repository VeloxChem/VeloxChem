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

import numpy as np
import sys

from .errorhandler import assert_msg_critical

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

    @staticmethod
    def _xc_label_to_dftd4(xc_label):
        """
        Determines dftd4 functional name from vlx functional name.

        :param xc_label:
            The vlx functional name.

        :return:
            The dftd4 functional name.
        """

        known_aliases = {
            'lrc-wpbeh': 'lc-wpbeh',
            'm06-l': 'm06l',
        }

        if xc_label.lower() in known_aliases:
            return known_aliases[xc_label.lower()]
        else:
            return xc_label.lower()

    def compute(self, molecule, xc_label):
        """
        Uses the dftd4-python interface to compute dispersion correction.

        :param molecule:
            The molecule.
        :param xc_label:
            The label of XC functional.
        """

        # sanity check
        errmsg = 'DispersionModel: dftd4-python is not available. '
        errmsg += 'Please install dftd4-python.'
        assert_msg_critical(self.is_available(), errmsg)

        identifiers_np = np.array(molecule.get_identifiers())
        coords_in_au = molecule.get_coordinates_in_bohr()
        net_charge = molecule.get_charge()

        disp_model = D4Model(numbers=identifiers_np,
                             positions=coords_in_au,
                             charge=net_charge,
                             model='d4')

        disp_xc_label = self._xc_label_to_dftd4(xc_label)
        disp_res = disp_model.get_dispersion(D4Param(method=disp_xc_label),
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

        return np.copy(self._gradient)

    def get_references(self):
        """
        Gets references.

        :return:
            The references.
        """

        dftd4_ref_1 = 'E. Caldeweyher, C. Bannwarth, S. Grimme,'
        dftd4_ref_1 += ' J. Chem. Phys., 2017, 147, 034112.'

        dftd4_ref_2 = 'E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer,'
        dftd4_ref_2 += ' S. Spicher, C. Bannwarth, S. Grimme,'
        dftd4_ref_2 += ' J. Chem Phys, 2019, 150, 154122.'

        return (dftd4_ref_1, dftd4_ref_2)
