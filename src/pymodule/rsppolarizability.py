from .rspproperty import ResponseProperty
from .rspdriver import ResponseDriver


class Polarizability(ResponseProperty):
    """Polarizability class"""

    def __init__(self, frequencies):
        """Initializes polarizability"""

        rsp_input = {
            'property': 'polarizability',
            'response': 'linear',
            'residue': 'none',
            'operators': ('xyz', 'xyz'),
            'frequencies': frequencies,
        }

        rsp_driver = ResponseDriver(rsp_input)

        super().__init__(rsp_input, rsp_driver)

    def get_property(self, key):
        """Gets polarizability component"""

        # key example: ('x', 'y', 0.1)

        return self.rsp_property[key]

    def print_property(self, ostream):
        """Prints polarizability to output stream"""

        for w in self.rsp_input['frequencies']:
            w_str = 'Polarizability (w={})'.format(w)
            ostream.print_header(w_str.ljust(68))
            ostream.print_header(('-' * len(w_str)).ljust(68))
            for a in self.rsp_input['operators'][0]:
                for b in self.rsp_input['operators'][1]:
                    ops_label = '<<{};{}>>_{}'.format(a, b, w)
                    output_alpha = '{:<15s} {:15.8f}'.format(
                        ops_label, self.rsp_property[(a, b, w)])
                    ostream.print_header(output_alpha.ljust(68))
            ostream.print_blank()
