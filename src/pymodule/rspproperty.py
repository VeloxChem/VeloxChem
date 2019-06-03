from .rspdriver import ResponseDriver


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

    def init_driver(self, comm, ostream):

        self.rsp_driver = ResponseDriver(comm, ostream)
        self.rsp_driver.update_settings(self.rsp_input)

    def compute(self, molecule, basis, scf_tensors):
        """Performs response property/spectroscopy calculation.

        Parameters
        ----------
        molecule
            The molecule.
        basis
            The AO basis set.
        scf_tensors
            The tensors from converged SCF wavefunction.
        """

        self.rsp_property = self.rsp_driver.compute(molecule, basis,
                                                    scf_tensors)

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
