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

from .veloxchemlib import hartree_in_ev
from .veloxchemlib import rotatory_strength_in_cgs
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

    def __init__(self, rsp_dict=None, method_dict=None):
        """
        Initializes the absorption property.
        """

        if rsp_dict is None:
            rsp_dict = {}
        else:
            rsp_dict = dict(rsp_dict)

        if method_dict is None:
            method_dict = {}
        else:
            method_dict = dict(method_dict)

        rsp_dict['property'] = 'absorption'
        rsp_dict['order'] = 'linear'
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

        return self._rsp_property[key]

    def print_property(self, ostream):
        """
        Prints absorption to output stream.

        :param ostream:
            The output stream.
        """

        spin_str = 'S'

        self._print_transition_dipoles(
            ostream, spin_str,
            'Electric Transition Dipole Moments (dipole length, a.u.)',
            self._rsp_property['electric_transition_dipoles'])

        self._print_transition_dipoles(
            ostream, spin_str,
            'Electric Transition Dipole Moments (dipole velocity, a.u.)',
            self._rsp_property['velocity_transition_dipoles'])

        self._print_transition_dipoles(
            ostream, spin_str, 'Magnetic Transition Dipole Moments (a.u.)',
            self._rsp_property['magnetic_transition_dipoles'])

        self._print_absorption(ostream, spin_str, 'One-Photon Absorption')
        self._print_ecd(ostream, spin_str, 'Electronic Circular Dichroism')
        self._print_excitation_details(ostream, 'Character of excitations:')

    def _print_transition_dipoles(self, ostream, spin_str, title,
                                  trans_dipoles):
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

    def _print_absorption(self, ostream, spin_str, title):
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
        for s, e in enumerate(self._rsp_property['eigenvalues']):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '{:15.8f} a.u. '.format(e)
            valstr += '{:12.5f} eV'.format(e * hartree_in_ev())
            f = self._rsp_property['oscillator_strengths'][s]
            valstr += '    Osc.Str. {:9.4f}'.format(f)
            ostream.print_header(valstr.ljust(92))
        ostream.print_blank()

    def _print_ecd(self, ostream, spin_str, title):
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
        for s, R in enumerate(self._rsp_property['rotatory_strengths']):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '    Rot.Str. '
            valstr += f'{(R / rotatory_strength_in_cgs()):13.6f} a.u.'
            valstr += f'{R:11.4f} [10**(-40) cgs]'
            ostream.print_header(valstr.ljust(92))
        ostream.print_blank()

    def _print_excitation_details(self, ostream, title):

        ostream.print_header(title.ljust(92))
        ostream.print_blank()

        nstates = self._rsp_property['eigenvalues'].size
        excitation_details = self._rsp_property['excitation_details']

        for s in range(nstates):
            valstr = 'Excited state {}'.format(s + 1)
            ostream.print_header(valstr.ljust(92))
            ostream.print_header(('-' * len(valstr)).ljust(92))

            for exc_str in excitation_details[s]:
                ostream.print_header(exc_str.ljust(92))
            ostream.print_blank()

        ostream.flush()
