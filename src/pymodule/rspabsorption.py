from .rspproperty import ResponseProperty
from .rspdriver import ResponseDriver


class Absorption(ResponseProperty):
    """Absorption class"""

    def __init__(self, nstates):
        """Initializes absorption"""

        rsp_input = {
            'property': 'absorption',
            'response': 'linear',
            'residue': 'single',
            'operators': ('xyz',),
            'nstates': nstates,
            'spin': 'singlet',
        }

        rsp_driver = ResponseDriver(rsp_input)

        super().__init__(rsp_input, rsp_driver)

    def get_property(self, state):
        """Gets absorption component"""

        return self.rsp_property[state]

    def print_property(self, ostream):
        """Prints absorption to output stream"""

        ostream.print_header('Absorption')
        ostream.print_header('----------')
        for s, e in enumerate(self.rsp_property):
            output_abs = 'State {:3d}: {:15.8f}'.format(s + 1, e)
            ostream.print_header(output_abs)
        ostream.print_blank()
