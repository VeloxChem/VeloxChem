#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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

from .rspproperty import ResponseProperty


class LinearAbsorptionCrossSection(ResponseProperty):
    """
    Implements the linear absorption cross-section property.

    :param rsp_dict:
        The dictionary of response input.
    :param method_dict:
        The dictionary of method settings.

    Instance variables
        - rsp_dict: The dictionary of response input.
        - method_dict: The dictionary of method settings.
        - rsp_property: The dictionary of response property.
    """

    def __init__(self, rsp_dict=None, method_dict=None):
        """
        Initializes the linear absorption cross-section property.
        """

        if rsp_dict is None:
            rsp_dict = {}
        else:
            rsp_dict = dict(rsp_dict)

        if method_dict is None:
            method_dict = {}
        else:
            method_dict = dict(method_dict)

        rsp_dict['property'] = 'linear absorption cross-section'
        rsp_dict['order'] = 'linear'
        rsp_dict['residue'] = 'none'
        rsp_dict['onlystatic'] = 'no'
        rsp_dict['complex'] = 'yes'

        rsp_dict['a_operator'] = 'electric dipole'
        rsp_dict['a_components'] = 'xyz'

        rsp_dict['b_operator'] = 'electric dipole'
        rsp_dict['b_components'] = 'xyz'

        if 'frequencies' not in rsp_dict:
            rsp_dict['frequencies'] = '0'

        super().__init__(rsp_dict, method_dict)

    def get_property(self, key):
        """
        Gets response functions or solutions.

        :param key:
            The keyword 'response_functions' or 'solutions'.

        :return:
            The response functions or solutions.
        """

        return self._rsp_property[key]

    def print_property(self, ostream):
        """
        Prints response property to output stream.

        :param ostream:
            The output stream.
        """

        pass
