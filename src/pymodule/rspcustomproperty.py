from .rspproperty import ResponseProperty
from .inputparser import parse_frequencies


class CustomProperty(ResponseProperty):
    """
    Implements a custom, user-defined  property.

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
        Initializes a custom, user-defined  property.
        """

        rsp_dict = dict(rsp_dict)

        if method_dict is None:
            method_dict = {}
        else:
            method_dict = dict(method_dict)

        rsp_dict['property'] = 'custom'
        if 'response' not in rsp_dict:
            rsp_dict['response'] = 'linear'
        rsp_dict['residue'] = 'none'
        rsp_dict['onlystatic'] = 'no'
        if 'complex' not in rsp_dict:
            rsp_dict['complex'] = 'no'

        if 'a_operator' not in rsp_dict:
            rsp_dict['a_operator'] = 'dipole'
        if 'a_components' not in rsp_dict:
            rsp_dict['a_components'] = 'xyz'

        if 'b_operator' not in rsp_dict:
            rsp_dict['b_operator'] = 'dipole'
        if 'b_components' not in rsp_dict:
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
        Prints custom property to output stream.

        :param ostream:
            The output stream.
        """

        width = 92

        for w in parse_frequencies(self.rsp_dict['frequencies']):
            w_str = 'Response Function (w={:.4f})'.format(w)
            ostream.print_header(w_str.ljust(width))
            ostream.print_header(('-' * len(w_str)).ljust(width))

            valstr = '{:<5s}'.format('')
            for b in self.rsp_dict['b_components']:
                if self.rsp_dict['complex'] == 'no':
                    valstr += '{:>15s}'.format(b.upper())
                else:
                    valstr += '{:>29s}'.format(b.upper())
            ostream.print_header(valstr.ljust(width))

            for a in self.rsp_dict['a_components']:
                valstr = '{:<5s}'.format(a.upper())
                for b in self.rsp_dict['b_components']:
                    prop = -self.rsp_property['response_functions'][(a, b, w)]
                    if self.rsp_dict['complex'] == 'no':
                        valstr += '{:15.8f}'.format(prop)
                    else:
                        valstr += '{:29.8f}'.format(prop)
                ostream.print_header(valstr.ljust(width))

            ostream.print_blank()
