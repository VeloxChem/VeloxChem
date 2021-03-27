from mpi4py import MPI
import numpy as np
import time as tm
import sys

from .molecule import Molecule
from .gradientdriver import GradientDriver
from .outputstream import OutputStream


class ScfGradientDriver(GradientDriver):
    """
    Implements SCF gradient driver.

    :param scf_drv:
        The SCF driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - delta_h: The displacement for finite difference.
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        super().__init__(comm, ostream)

        self.flag = 'SCF Gradient Driver'
        self.scf_drv = scf_drv
        self.delta_h = 0.001

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

        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # numerical gradient
        self.gradient = np.zeros((molecule.number_of_atoms(), 3))

        for i in range(molecule.number_of_atoms()):
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

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        self.scf_drv.ostream.state = scf_ostream_state

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule, labels)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()
