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

import numpy as np

from .rspproperty import ResponseProperty


class C6(ResponseProperty):
    """
    Implements the C6 dispersion coefficient  property.

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
        Initializes the C6 dispersion coefficient property.
        """

        if rsp_dict is None:
            rsp_dict = {}
        else:
            rsp_dict = dict(rsp_dict)

        if method_dict is None:
            method_dict = {}
        else:
            method_dict = dict(method_dict)

        rsp_dict['property'] = 'c6'
        rsp_dict['order'] = 'linear'
        rsp_dict['residue'] = 'none'
        rsp_dict['onlystatic'] = 'yes'
        rsp_dict['complex'] = 'yes'

        rsp_dict['a_operator'] = 'electric dipole'
        rsp_dict['a_components'] = 'xyz'

        rsp_dict['b_operator'] = 'electric dipole'
        rsp_dict['b_components'] = 'xyz'

        if 'n_points' not in rsp_dict:
            rsp_dict['n_points'] = '9'
        if 'w0' not in rsp_dict:
            rsp_dict['w0'] = '0.3'

        super().__init__(rsp_dict, method_dict)

    def get_property(self, key):
        """
        Gets excitation energies, CI vectors, or oscillator stengths.

        :param key:
            The keyword to the C6 property.

        :return:
            The C6 property.
        """

        return self._rsp_property[key]

    def print_property(self, ostream):
        """
        Prints response property to output stream.

        :param ostream:
            The output stream.
        """

        width = 92

        title = 'Response Functions at Given Imaginary Frequencies'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        w0 = float(self._rsp_dict['w0'])
        n_points = int(self._rsp_dict['n_points'])
        points, weights = np.polynomial.legendre.leggauss(n_points)
        imagfreqs = [w0 * (1 - t) / (1 + t) for t in points]
        printfreqs = np.append(imagfreqs, 0.0)

        for iw in printfreqs:
            title = '{:<7s} {:<7s} {:>10s} {:>15s} {:>16s}'.format(
                'Dipole', 'Dipole', 'Frequency', 'Real', 'Imaginary')
            ostream.print_header(title.ljust(width))
            ostream.print_header(('-' * len(title)).ljust(width))

            for a in self._rsp_dict['a_components']:
                for b in self._rsp_dict['b_components']:
                    prop = self._rsp_property['response_functions'][(a, b, iw)]
                    ops_label = '<<{:>3s}  ;  {:<3s}>> {:10.4f}'.format(
                        a.lower(), b.lower(), iw)
                    output = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, prop.real, prop.imag)
                    ostream.print_header(output.ljust(width))
            ostream.print_blank()

        title = self._rsp_driver.get_prop_str()
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        title = 'Reference: '
        title += 'Amos et al., '
        title += 'J. Chem. Phys. 89, 2186 (1985).'
        ostream.print_header(title.ljust(width))
        ostream.print_blank()

        c6 = self._rsp_property['c6']

        Gxx_i0 = self._rsp_property['response_functions'][('x', 'x', 0.0)].real
        Gyy_i0 = self._rsp_property['response_functions'][('y', 'y', 0.0)].real
        Gzz_i0 = self._rsp_property['response_functions'][('z', 'z', 0.0)].real

        alpha_i0 = -(Gxx_i0 + Gyy_i0 + Gzz_i0) / 3.0

        output = 'Homomolecular C_6 value        :    {:10.6f} a.u.'.format(c6)
        ostream.print_header(output.ljust(width))
        ostream.print_blank()
        output = 'Static polarizability alpha(0) :    {:10.6f} a.u.'.format(
            alpha_i0)
        ostream.print_header(output.ljust(width))

        ostream.print_blank()
