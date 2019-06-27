from .rspproperty import ResponseProperty
from .inputparser import parse_frequencies


class Polarizability(ResponseProperty):
    """
    Implements the polarizability property.

    :param rsp_input:
        The dictionary of response input.
    :param rsp_property:
        The dictionary of response property.
    """

    def __init__(self, rsp_dict):
        """
        Initializes the polarizability property.

        :param rsp_dict:
            The dictionary of response input.
        """

        rsp_input = dict(rsp_dict)

        rsp_input['property'] = 'polarizability'
        rsp_input['response'] = 'linear'
        rsp_input['residue'] = 'none'

        rsp_input['a_operator'] = 'dipole'
        rsp_input['a_components'] = 'xyz'

        rsp_input['b_operator'] = 'dipole'
        rsp_input['b_components'] = 'xyz'

        if 'frequencies' not in rsp_dict:
            rsp_input['frequencies'] = '0'

        super().__init__(rsp_input)

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

        for w in parse_frequencies(self.rsp_input['frequencies']):
            w_str = 'Polarizability (w={:.4f})'.format(w)
            ostream.print_header(w_str.ljust(width))
            ostream.print_header(('-' * len(w_str)).ljust(width))

            valstr = '{:<5s}'.format('')
            for b in self.rsp_input['b_components']:
                valstr += '{:>15s}'.format(b.upper())
            ostream.print_header(valstr.ljust(width))

            for a in self.rsp_input['a_components']:
                valstr = '{:<5s}'.format(a.upper())
                for b in self.rsp_input['b_components']:
                    prop = -self.rsp_property[(a, b, w)]
                    valstr += '{:15.8f}'.format(prop)
                ostream.print_header(valstr.ljust(width))

            ostream.print_blank()
