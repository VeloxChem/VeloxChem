from .veloxchemlib import hartree_in_ev
from .rspproperty import ResponseProperty


class Absorption(ResponseProperty):
    """Absorption class"""

    def __init__(self, rsp_dict):
        """Initializes absorption"""

        rsp_input = dict(rsp_dict)

        rsp_input['property'] = 'absorption'
        rsp_input['response'] = 'linear'
        rsp_input['residue'] = 'single'

        if 'nstates' not in rsp_dict:
            rsp_input['nstates'] = '3'
        if 'tamm_dancoff' not in rsp_dict:
            rsp_input['tamm_dancoff'] = 'no'
        if 'spin' not in rsp_dict:
            rsp_input['spin'] = 'singlet'

        super().__init__(rsp_input)

    def get_property(self, key):
        """Gets excitation energies, CI vectors, or oscillator stengths"""

        return self.rsp_property[key]

    def print_property(self, ostream):
        """Prints absorption to output stream"""

        ostream.print_blank()
        ostream.print_header('One-Photon Absorption'.ljust(92))
        ostream.print_header('---------------------'.ljust(92))
        spin_str = 'T' if self.rsp_input['spin'][0].upper() == 'T' else 'S'
        for s, e in enumerate(self.rsp_property['eigenvalues']):
            output_abs = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            output_abs += '{:15.8f} a.u. '.format(e)
            output_abs += '{:12.5f} eV'.format(e * hartree_in_ev())
            f = self.rsp_property['oscillator_strengths'][s]
            output_abs += '    osc.str.{:12.5f}'.format(f)
            ostream.print_header(output_abs.ljust(92))
        ostream.print_blank()

        ostream.print_header('Electronic Circular Dichroism'.ljust(92))
        ostream.print_header('-----------------------------'.ljust(92))
        spin_str = 'T' if self.rsp_input['spin'][0].upper() == 'T' else 'S'
        for s, e in enumerate(self.rsp_property['eigenvalues']):
            output_abs = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            output_abs += '{:15.8f} a.u. '.format(e)
            output_abs += '{:12.5f} eV'.format(e * hartree_in_ev())
            r = self.rsp_property['rotatory_strengths'][s]
            output_abs += '    rot.str.{:12.5f}'.format(r)
            ostream.print_header(output_abs.ljust(92))
        ostream.print_blank()
