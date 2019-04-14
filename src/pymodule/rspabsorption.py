from .rspproperty import ResponseProperty
from .veloxchemlib import hartree_in_ev


class Absorption(ResponseProperty):
    """Absorption class"""

    def __init__(self, rsp_dict):
        """Initializes absorption"""

        rsp_input = dict(rsp_dict)

        rsp_input['property'] = 'absorption'
        rsp_input['response'] = 'linear'
        rsp_input['residue'] = 'single'
        # rsp_input['operators'] = ('xyz',)

        # if 'spin' in rsp_dict:
        #    rsp_input['spin'] = rsp_dict['spin']
        # else:
        #    rsp_input['spin'] = 'singlet'

        super().__init__(rsp_input)

    def get_property(self, state):
        """Gets absorption component"""

        return self.rsp_property[state]

    def print_property(self, ostream):
        """Prints absorption to output stream"""

        ostream.print_blank()
        ostream.print_header('One-Photon Absorption'.ljust(92))
        ostream.print_header('---------------------'.ljust(92))
        for s, e in enumerate(self.rsp_property):
            output_abs = 'Excited State {:5d}: '.format(s + 1)
            output_abs += '{:15.8f} a.u. '.format(e)
            output_abs += '{:12.5f} eV'.format(e * hartree_in_ev())
            ostream.print_header(output_abs.ljust(92))
        ostream.print_blank()
