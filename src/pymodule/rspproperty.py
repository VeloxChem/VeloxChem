from .rspdriver import ResponseDriver
from .outputstream import OutputStream

import sys


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

        if 'conv_thresh' in rsp_input:
            rsp_input['conv_thresh'] = float(rsp_input['conv_thresh'])
        if 'max_iter' in rsp_input:
            rsp_input['max_iter'] = int(rsp_input['max_iter'])
        if 'eri_thresh' in rsp_input:
            rsp_input['eri_thresh'] = float(rsp_input['eri_thresh'])
        if 'qq_type' in rsp_input:
            rsp_input['qq_type'] = rsp_input['qq_type'].upper()

        self.rsp_input = rsp_input
        self.rsp_driver = ResponseDriver(rsp_input)

    def compute_task(self, mol_orbs, task):
        """Performs response property/spectroscopy calculation.

        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        task
            The MPI task.
        """

        self.rsp_property = self.rsp_driver.compute(
            mol_orbs, task.molecule, task.ao_basis, task.mpi_comm, task.ostream)

    def compute(self, mol_orbs, molecule, basis, comm,
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

        self.rsp_property = self.rsp_driver.compute(mol_orbs, molecule, basis,
                                                    comm, ostream)

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
