from .rspproperty import ResponseProperty
from .inputparser import parse_frequencies


class Polarizability(ResponseProperty):
    """
    Implements the polarizability property.

    :param rsp_dict:
        The dictionary of response input.
    :param method_dict:
        The dictionary of method settings.
    :param rsp_property:
        The dictionary of response property.
    """

    def __init__(self, rsp_dict, method_dict={}):
        """
        Initializes the polarizability property.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        rsp_dict = dict(rsp_dict)
        method_dict = dict(method_dict)

        rsp_dict['property'] = 'polarizability'
        rsp_dict['response'] = 'linear'
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

        for w in parse_frequencies(self.rsp_dict['frequencies']):
            w_str = 'Polarizability (w={:.4f})'.format(w)
            ostream.print_header(w_str.ljust(width))
            ostream.print_header(('-' * len(w_str)).ljust(width))

            valstr = '{:<5s}'.format('')
            for b in self.rsp_dict['b_components']:
                valstr += '{:>15s}'.format(b.upper())
            ostream.print_header(valstr.ljust(width))

            for a in self.rsp_dict['a_components']:
                valstr = '{:<5s}'.format(a.upper())
                for b in self.rsp_dict['b_components']:
                    prop = -self.rsp_property[(a, b, w)]
                    valstr += '{:15.8f}'.format(prop)
                ostream.print_header(valstr.ljust(width))

            ostream.print_blank()
