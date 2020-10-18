import numpy as np
import time as tm

from .molecule import Molecule
from .outputstream import OutputStream
from .scfrestdriver import ScfRestrictedDriver


class GradientDriver:
    """
    Implements gradient driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - delta_h: The displacement for finite difference.
        - gradient: The gradient.
    """

    def __init__(self, comm, ostream):
        """
        Initializes gradient driver.
        """

        self.comm = comm
        self.ostream = ostream

        self.scf_drv = ScfRestrictedDriver(self.comm, OutputStream())

        self.delta_h = 0.001
        self.gradient = None

    def update_settings(self, scf_dict, method_dict=None):
        """
        Updates settings in gradient driver.

        :param scf_dict:
            The input dictionary of scf group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        self.scf_drv.update_settings(scf_dict, method_dict)

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Performs calculation of numerical gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        self.print_header()
        start_time = tm.time()

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # numerical gradient
        self.gradient = np.zeros((molecule.number_of_atoms(), 3))

        for i in range(molecule.number_of_atoms()):
            self.ostream.print_info(
                'Performing finite difference on atom {:d}...'.format(i + 1))
            self.ostream.flush()

            for d in range(3):
                coords[i, d] += self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                self.scf_drv.compute(new_mol, ao_basis, min_basis)
                e_plus = self.scf_drv.get_scf_energy()

                coords[i, d] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                self.scf_drv.compute(new_mol, ao_basis, min_basis)
                e_minus = self.scf_drv.get_scf_energy()

                coords[i, d] += self.delta_h
                self.gradient[i, d] = (e_plus - e_minus) / (2.0 * self.delta_h)

        self.ostream.print_blank()

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule, labels)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
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
        Prints gradient calculation setup details to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Numerical Gradient Driver")
        self.ostream.print_header(27 * "=")
        self.ostream.print_blank()

