import sys

from .rspdriver import ResponseDriver
from .outputstream import OutputStream


class ResponseProperty:
    """Implements response property/spectroscopy class.

    Implements the base class for response property/spectroscopy.
    """

    def __init__(self, rsp_input):
        """Initializes response property/spectroscopy.

        Parameters
        ----------
        rsp_input
            The input dictionary that defines the property/spectroscopy.
        """

        self.rsp_input = rsp_input

    def compute_task(self, mol_orbs, task):
        """Performs response property/spectroscopy calculation.

        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        task
            The MPI task.
        """

        rsp_driver = ResponseDriver(self.rsp_input, task.mpi_comm, task.ostream)
        self.rsp_property = rsp_driver.compute(mol_orbs, task.molecule,
                                               task.ao_basis)

    def compute(self,
                mol_orbs,
                molecule,
                basis,
                comm,
                ostream=OutputStream(sys.stdout)):
        """Performs response property/spectroscopy calculation.

        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        molecule
            The molecule.
        basis
            The AO basis set.
        comm
            The MPI communicator.
        ostream
            The output stream.
        """

        rsp_driver = ResponseDriver(self.rsp_input, comm, ostream)
        self.rsp_property = rsp_driver.compute(mol_orbs, molecule, basis)

    def get_property(self, req_list):
        """Gets response property/spectroscopy.

        Parameters
        ----------
        req_list
            The requested component of the response property.
        """

        return None

    def print_property(self, ostream):
        """Prints response property/spectroscopy to output stream."""

        return
