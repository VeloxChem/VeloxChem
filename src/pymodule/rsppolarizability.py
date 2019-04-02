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

        ostream.print_header('Polarizability')
        ostream.print_header('--------------')
        for (a, b, w), alpha_abw in self.rsp_property.items():
            ops_label = '<<{};{}>>_{}'.format(a, b, w)
            output_alpha = '{:<15s} {:15.8f}'.format(ops_label, alpha_abw)
            ostream.print_header(output_alpha)
        ostream.print_blank()
