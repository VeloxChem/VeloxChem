from mpi4py import MPI
import sys

from .veloxchemlib import mpi_master
from .gradientdriver import GradientDriver
from .outputstream import OutputStream


class XTBGradientDriver(GradientDriver):
    """
    Implements XTB gradient driver.

    :param xtb_drv:
        The XTB driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - xtb_drv: The XTB driver.
    """

    def __init__(self, xtb_drv, comm=None, ostream=None):
        """
        Initializes XTB gradient driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        super().__init__(comm, ostream)

        self.flag = 'XTB Gradient Driver'
        self.xtb_drv = xtb_drv

    def compute(self, molecule):
        """
        Performs calculation of XTB analytical gradient.

        :param molecule:
            The molecule.
        """

        self.print_header()

        self.gradient = self.xtb_drv.get_gradient()
        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())

        self.print_geometry(molecule)
        self.print_gradient(molecule, molecule.get_labels())

        self.ostream.print_blank()
        self.ostream.flush()