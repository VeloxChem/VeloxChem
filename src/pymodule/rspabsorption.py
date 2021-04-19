#
#                           VELOXCHEM 1.0-RC
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
from .rspproperty import ResponseProperty


class Absorption(ResponseProperty):
    """
    Implements the absorption property.

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
        Initializes the absorption property.
        """

        rsp_dict = dict(rsp_dict)

        if method_dict is None:
            method_dict = {}
        else:
            method_dict = dict(method_dict)

        rsp_dict['property'] = 'absorption'
        rsp_dict['response'] = 'linear'
        rsp_dict['residue'] = 'single'
        rsp_dict['complex'] = 'no'

        if 'nstates' not in rsp_dict:
            rsp_dict['nstates'] = '3'
        if 'tamm_dancoff' not in rsp_dict:
            rsp_dict['tamm_dancoff'] = 'no'

        super().__init__(rsp_dict, method_dict)

    def get_property(self, key):
        """
        Gets excitation energies, CI vectors, or oscillator stengths.

        :param key:
            The keyword to the absorption property.

        :return:
            The absorption property.
        """

        return self.rsp_property[key]

    def print_property(self, ostream):
        """
        Prints absorption to output stream.

        :param ostream:
            The output stream.
        """

        spin_str = 'S'

        self.print_transition_dipoles(
            ostream, spin_str,
            'Electric Transition Dipole Moments (dipole length, a.u.)',
            self.rsp_property['electric_transition_dipoles'])

        self.print_transition_dipoles(
            ostream, spin_str,
            'Electric Transition Dipole Moments (dipole velocity, a.u.)',
            self.rsp_property['velocity_transition_dipoles'])

        self.print_transition_dipoles(
            ostream, spin_str, 'Magnetic Transition Dipole Moments (a.u.)',
            self.rsp_property['magnetic_transition_dipoles'])

        self.print_absorption(ostream, spin_str, 'One-Photon Absorption')
        self.print_ecd(ostream, spin_str, 'Electronic Circular Dichroism')

    def print_transition_dipoles(self, ostream, spin_str, title, trans_dipoles):
        """
        Prints transition dipole moments to output stream.

        :param ostream:
            The output stream.
        :param spin_str:
            The string representation of spin.
        :param title:
            The title to be printed to the output stream.
        :param trans_dipoles:
            The transition dipole moments.
        """

        valstr = title
        ostream.print_header(valstr.ljust(92))
        ostream.print_header(('-' * len(valstr)).ljust(92))
        valstr = '                     '
        valstr += '{:>13s}{:>13s}{:>13s}'.format('X', 'Y', 'Z')
        ostream.print_header(valstr.ljust(92))
        for s, r in enumerate(trans_dipoles):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '{:13.6f}{:13.6f}{:13.6f}'.format(r[0], r[1], r[2])
            ostream.print_header(valstr.ljust(92))
        ostream.print_blank()

    def print_absorption(self, ostream, spin_str, title):
        """
        Prints absorption to output stream.

        :param ostream:
            The output stream.
        :param spin_str:
            The string representation of spin.
        :param title:
            The title to be printed to the output stream.
        """

        valstr = title
        ostream.print_header(valstr.ljust(92))
        ostream.print_header(('-' * len(valstr)).ljust(92))
        for s, e in enumerate(self.rsp_property['eigenvalues']):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '{:15.8f} a.u. '.format(e)
            valstr += '{:12.5f} eV'.format(e * hartree_in_ev())
            f = self.rsp_property['oscillator_strengths'][s]
            valstr += '    Osc.Str. {:9.4f}'.format(f)
            ostream.print_header(valstr.ljust(92))
        ostream.print_blank()

    def print_ecd(self, ostream, spin_str, title):
        """
        Prints electronic circular dichroism to output stream.

        :param ostream:
            The output stream.
        :param spin_str:
            The string representation of spin.
        :param title:
            The title to be printed to the output stream.
        """

        valstr = title
        ostream.print_header(valstr.ljust(92))
        ostream.print_header(('-' * len(valstr)).ljust(92))
        for s, R in enumerate(self.rsp_property['rotatory_strengths']):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '    Rot.Str. {:11.4f}'.format(R)
            valstr += '    [10**(-40) (esu**2)*(cm**2)]'
            ostream.print_header(valstr.ljust(92))
        ostream.print_blank()
