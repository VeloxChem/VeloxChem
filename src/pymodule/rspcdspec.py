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

from .veloxchemlib import hartree_in_ev
from .veloxchemlib import hartree_in_wavenumbers
from .veloxchemlib import molar_ellipticity_from_beta
from .rspproperty import ResponseProperty
from .inputparser import parse_seq_range


class CircularDichroismSpectrum(ResponseProperty):
    """
    Implements the circular dichroism spectrum property.

    :param rsp_dict:
        The dictionary of response input.
    :param method_dict:
        The dictionary of method settings.

    Instance variables
        - rsp_dict: The dictionary of response input.
        - method_dict: The dictionary of method settings.
        - rsp_property: The dictionary of response property.
    """

    def __init__(self, rsp_dict, method_dict=None):
        """
        Initialized the circular dichroism spectrum property.
        """

        rsp_dict = dict(rsp_dict)

        if method_dict is None:
            method_dict = {}
        else:
            method_dict = dict(method_dict)

        rsp_dict['property'] = 'circular dichroism spectrum'
        rsp_dict['response'] = 'linear'
        rsp_dict['residue'] = 'none'
        rsp_dict['onlystatic'] = 'no'
        rsp_dict['complex'] = 'yes'

        rsp_dict['a_operator'] = 'magnetic dipole'
        rsp_dict['a_components'] = 'xyz'

        rsp_dict['b_operator'] = 'linear momentum'
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

        return self.rsp_property[key]

    def print_property(self, ostream):
        """
        Prints response property to output stream.

        :param ostream:
            The output stream.
        """

        width = 92

        title = 'Response Functions at Given Frequencies'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        freqs = parse_seq_range(self.rsp_dict['frequencies'])

        for w in freqs:
            title = '{:<7s} {:<7s} {:>10s} {:>15s} {:>16s}'.format(
                'MagDip', 'LinMom', 'Frequency', 'Real', 'Imaginary')
            ostream.print_header(title.ljust(width))
            ostream.print_header(('-' * len(title)).ljust(width))

            for a in self.rsp_dict['a_components']:
                for b in self.rsp_dict['b_components']:
                    prop = self.rsp_property['response_functions'][(a, b, w)]
                    ops_label = '<<{:>3s}  ;  {:<3s}>> {:10.4f}'.format(
                        a.lower(), b.lower(), w)
                    output = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, prop.real, prop.imag)
                    ostream.print_header(output.ljust(width))
            ostream.print_blank()

        title = self.rsp_driver.prop_str()
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        if len(freqs) == 1 and freqs[0] == 0.0:
            text = '*** No circular dichroism spectrum at zero frequency.'
            ostream.print_header(text.ljust(width))
            ostream.print_blank()
            return

        title = 'Reference: '
        title += 'A. Jiemchooroj and P. Norman, '
        title += 'J. Chem. Phys. 126, 134102 (2007).'
        ostream.print_header(title.ljust(width))
        ostream.print_blank()

        title = '{:<20s}{:<20s}{:>28s}'.format('Frequency[a.u.]',
                                               'Frequency[eV]',
                                               'Delta_epsilon[L mol^-1 cm^-1]')
        ostream.print_header(title.ljust(width))
        ostream.print_header(('-' * len(title)).ljust(width))

        for w in freqs:
            if w == 0.0:
                continue

            Gxx = self.rsp_property['response_functions'][('x', 'x', w)].imag
            Gyy = self.rsp_property['response_functions'][('y', 'y', w)].imag
            Gzz = self.rsp_property['response_functions'][('z', 'z', w)].imag

            beta = -(Gxx + Gyy + Gzz) / (3.0 * w**2)
            wavenumber = w * hartree_in_wavenumbers()
            Delta_epsilon = (beta * wavenumber**2
                             * molar_ellipticity_from_beta() / (100.0 * 3298.8))

            output = '{:<20.4f}{:<20.5f}{:>18.8f}'.format(
                w, w * hartree_in_ev(), Delta_epsilon)
            ostream.print_header(output.ljust(width))

        ostream.print_blank()
