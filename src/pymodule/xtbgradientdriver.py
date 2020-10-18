import numpy as np
import time as tm

from .molecule import Molecule
from .outputstream import OutputStream
from .scfrestdriver import ScfRestrictedDriver


class XTBGradientDriver:
    """
    Implements XTB gradient driver.

    :param comm:
        The MPI communicator. 
    :param xtb_drv: 
        The XTB driver.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm, xtb_drv, ostream):
        """
        Initializes XTB gradient driver.
        """

        self.comm = comm
        self.ostream = ostream
        self.xtb_drv = xtb_drv

    def compute(self, molecule): 
        """
        Performs calculation of XTB analytical gradient.

        :param molecule:
            The molecule.
        """

        self.print_header()

        self.gradient = self.xtb_drv.get_gradient()
        
        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule, molecule.get_labels())

        self.ostream.print_blank()
        self.ostream.flush()

    def get_gradient(self):
        """
        Gets the gradient.

        :return:
            The gradient.
        """

        return self.gradient

    def print_geometry(self, molecule):
        """
        Prints the gradient.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())

    def print_gradient(self, molecule, labels):
        """
        Prints the gradient.

        :param molecule:
            The molecule.
        :param labels:
            The atom labels.
        """

        title = 'Gradient (Hartree/Bohr)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))
        self.ostream.print_blank()

        valstr = '  Atom '
        valstr += '{:>20s}  '.format('Gradient X')
        valstr += '{:>20s}  '.format('Gradient Y')
        valstr += '{:>20s}  '.format('Gradient Z')
        self.ostream.print_header(valstr)
        self.ostream.print_blank()

        for i in range(molecule.number_of_atoms()):
            valstr = '  {:<4s}'.format(labels[i])
            for d in range(3):
                valstr += '{:22.12f}'.format(self.gradient[i, d])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self):
        """
        Prints XTB gradient calculation setup details to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Analytical XTB Gradient Driver")
        self.ostream.print_header("==============================")
        self.ostream.print_blank()
