from .veloxchemlib import hartree_in_ev
from .rspproperty import ResponseProperty


class Absorption(ResponseProperty):
    """
    Implements the absorption property.

    :param rsp_input:
        The dictionary of response input.
    :param rsp_property:
        The dictionary of response property.
    """

    def __init__(self, rsp_dict):
        """
        Initializes the absorption property.

        :param rsp_dict:
            The dictionary of response input.
        """

        rsp_input = dict(rsp_dict)

        rsp_input['property'] = 'absorption'
        rsp_input['response'] = 'linear'
        rsp_input['residue'] = 'single'

        if 'nstates' not in rsp_dict:
            rsp_input['nstates'] = '3'
        if 'tamm_dancoff' not in rsp_dict:
            rsp_input['tamm_dancoff'] = 'no'
        if 'spin' not in rsp_dict:
            rsp_input['spin'] = 'singlet'

        super().__init__(rsp_input)

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

        ostream.print_blank()
        spin_str = 'T' if self.rsp_input['spin'][0].upper() == 'T' else 'S'

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
