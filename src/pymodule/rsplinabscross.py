from .rspproperty import ResponseProperty
from .inputparser import parse_frequencies


class LinearAbsorptionCrossSection(ResponseProperty):
    """
    Implements the linear absorption cross-section property.

    :param rsp_dict:
        The dictionary of response input.
    :param method_dict:
        The dictionary of method settings.
    :param rsp_property:
        The dictionary of response property.
    """

    def __init__(self, rsp_dict, method_dict={}):
        """
        Initializes the linear absorption cross-section property.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        rsp_dict = dict(rsp_dict)
        method_dict = dict(method_dict)

        rsp_dict['property'] = 'linear absorption cross-section'
        rsp_dict['response'] = 'linear'
        rsp_dict['complex'] = 'yes'
        rsp_dict['residue'] = 'none'

        rsp_dict['a_operator'] = 'dipole'
        rsp_dict['a_components'] = 'xyz'

        rsp_dict['b_operator'] = 'dipole'
        rsp_dict['b_components'] = 'xyz'

        if 'frequencies' not in rsp_dict:
            rsp_dict['frequencies'] = '0'

        super().__init__(rsp_dict, method_dict)

    def get_property(self, key):
        """
        Gets component of polarizability.

        :param key:
            The tuple of A component, B component, and frequency.

        :return:
            The component of polarizability.
        """

        # key example: ('x', 'y', 0.1)
        return self.rsp_property[key]

    def print_property(self, ostream):
        """
        Prints polarizability to output stream.

        :param ostream:
            The output stream.
        """

        width = 92

        title = self.rsp_driver.prop_str()
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        for w in parse_frequencies(self.rsp_dict['frequencies']):
            title = '{:<8s} {:<8s} {:>10s} {:>15s} {:>16s}'.format(
                'Dipole', 'Dipole', 'Frequency', 'Real', 'Imaginary')
            ostream.print_header(title.ljust(width))
            ostream.print_header(('-' * len(title)).ljust(width))

            for a in self.rsp_dict['a_components']:
                for b in self.rsp_dict['b_components']:
                    prop = -self.rsp_property['properties'][(a, b, w)]
                    ops_label = '{:<8s} {:<8s} {:10.4f}'.format(
                        a.upper(), b.upper(), w)
                    output = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, prop.real, prop.imag)
                    ostream.print_header(output.ljust(width))
            ostream.print_blank()
