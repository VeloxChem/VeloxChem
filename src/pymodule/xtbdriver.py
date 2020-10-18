import numpy as np

from .veloxchemlib import XTBDriver

def _XTBDriver_compute_energy(self, molecule, scf_dict, method_dict, ostream):
    """
    Computes DFT-B energy using XTB package.

    :param molecule:
        The molecule.
    :param scf_dict:
        The input dictionary of scf group.
    :param method_dict:
        The input dicitonary of method settings group.
    :param ostream:
        The output stream.
    """

    self.compute(molecule, method_dict['xtb'].lower())
    
    with open("xtb.scf.tempfile",'r') as f:
        for line in f.readlines(): 
            ostream.print_line(line.rstrip('\n'))
    
    return

XTBDriver.compute_energy = _XTBDriver_compute_energy
