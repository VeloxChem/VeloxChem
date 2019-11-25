import math

from .veloxchemlib import hartree_in_ev
from .rspproperty import ResponseProperty
from .inputparser import parse_frequencies


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

    def __init__(self, rsp_dict, method_dict={}):
        """
        Initializes the linear absorption cross-section property.
        """

        rsp_dict = dict(rsp_dict)
        method_dict = dict(method_dict)

        rsp_dict['property'] = 'linear absorption cross-section'
        rsp_dict['response'] = 'linear'
        rsp_dict['residue'] = 'none'
        rsp_dict['complex'] = 'yes'

        rsp_dict['a_operator'] = 'dipole'
        rsp_dict['a_components'] = 'xyz'

        rsp_dict['b_operator'] = 'dipole'
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

        for w in parse_frequencies(self.rsp_dict['frequencies']):
            title = '{:<7s} {:<7s} {:>10s} {:>15s} {:>16s}'.format(
                'Dipole', 'Dipole', 'Frequency', 'Real', 'Imaginary')
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

        freqs = parse_frequencies(self.rsp_dict['frequencies'])

        if len(freqs) == 1 and freqs[0] == 0.0:
            text = '*** No linear absorption spectrum at zero frequency.'
            ostream.print_header(text.ljust(width))
            ostream.print_blank()
            return

        title = 'Reference: '
        title += 'J. Kauczor and P. Norman, '
        title += 'J. Chem. Theory Comput. 2014, 10, 2449-2455.'
        ostream.print_header(title.ljust(width))
        ostream.print_blank()

        title = '{:<20s}{:<20s}{:>15s}'.format('Frequency[a.u.]',
                                               'Frequency[eV]',
                                               'sigma(w)[a.u.]')
        ostream.print_header(title.ljust(width))
        ostream.print_header(('-' * len(title)).ljust(width))

        for w in freqs:
            if w == 0.0:
                continue

            axx = -self.rsp_property['response_functions'][('x', 'x', w)].imag
            ayy = -self.rsp_property['response_functions'][('y', 'y', w)].imag
            azz = -self.rsp_property['response_functions'][('z', 'z', w)].imag

            alpha_bar = (axx + ayy + azz) / 3.0
            sigma = 4.0 * math.pi * w * alpha_bar / 137.035999

            output = '{:<20.4f}{:<20.5f}{:>13.8f}'.format(
                w, w * hartree_in_ev(), sigma)
            ostream.print_header(output.ljust(width))

        ostream.print_blank()
