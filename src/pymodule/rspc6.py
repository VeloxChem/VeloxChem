import numpy as np
import math

from .veloxchemlib import hartree_in_ev
from .rspproperty import ResponseProperty


class C6(ResponseProperty):
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

    def __init__(self, rsp_dict, method_dict={}):
        """
        Initializes the absorption property.
        """

        rsp_dict = dict(rsp_dict)
        method_dict = dict(method_dict)

        rsp_dict['property'] = 'c6'
        rsp_dict['response'] = 'linear'
        rsp_dict['residue'] = 'none'
        rsp_dict['onlystatic'] = 'yes'
        rsp_dict['complex'] = 'yes'

        rsp_dict['a_operator'] = 'dipole'
        rsp_dict['a_components'] = 'xyz'

        rsp_dict['b_operator'] = 'dipole'
        rsp_dict['b_components'] = 'xyz'

        if 'n_points' not in rsp_dict:
            rsp_dict['n_points'] = '9'

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
        Prints response property to output stream.

        :param ostream:
            The output stream.
        """

        width = 92

        title = 'Response Functions at Given Imaginary Frequencies'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        for iw in [0.3*(1-t)/(1+t) for t in np.polynomial.legendre.leggauss(int(self.rsp_dict['n_points']))][0]:
            title = '{:<7s} {:<7s} {:>10s} {:>15s} {:>16s}'.format(
                'Dipole', 'Dipole', 'Frequency', 'Real', 'Imaginary')
            ostream.print_header(title.ljust(width))
            ostream.print_header(('-' * len(title)).ljust(width))

            for a in self.rsp_dict['a_components']:
                for b in self.rsp_dict['b_components']:
                    prop = self.rsp_property['response_functions'][(a, b, iw)]
                    ops_label = '<<{:>3s}  ;  {:<3s}>> {:10.4f}'.format(
                        a.lower(), b.lower(), iw)
                    output = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, prop.real, prop.imag)
                    ostream.print_header(output.ljust(width))
            ostream.print_blank()

        title = self.rsp_driver.prop_str()
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        imagfreqs = [0.3*(1-t)/(1+t) for t in np.polynomial.legendre.leggauss(int(self.rsp_dict['n_points']))][0]
        weights = np.polynomial.legendre.leggauss(int(self.rsp_dict['n_points']))[1]
        integral = 0

        for iw in range(len(imagfreqs)):

            Gxx = self.rsp_property['response_functions'][('x', 'x', imagfreqs[iw])].real
            Gyy = self.rsp_property['response_functions'][('y', 'y', imagfreqs[iw])].real
            Gzz = self.rsp_property['response_functions'][('z', 'z', imagfreqs[iw])].real

            alpha = (Gxx + Gyy + Gzz) / 3.0
            weight = weights[iw]
            integral += alpha**2 * weight

        c6 = 3 * integral / math.pi

        output = 'C_6 value :    {:10.6f}'.format(c6)
        ostream.print_header(output.ljust(width))

        ostream.print_blank()
