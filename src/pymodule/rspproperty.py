from .rspdriver import ResponseDriver


class ResponseProperty:
    """Implements response property/spectroscopy class.

    Implements the base class for response property/spectroscopy.

    Attributes
    ----------
    rsp_input
        The dictionary of response input.
    rsp_driver
        The response driver.
    rsp_property
        The dictionary of response property.
    """

    def __init__(self, rsp_input):
        """Initializes response property/spectroscopy.

        Initializes response property/spectroscopy.

        Parameters
        ----------
        rsp_input
            The input dictionary that defines the property/spectroscopy.
        """

        self.rsp_input = rsp_input

    def init_driver(self, comm, ostream):
        """Initializes response driver.

        Initializes response driver.

        Parameters
        ----------
        comm
            The MPI communicator.
        ostream
            The output stream.
        """

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
            The dictionary of tensors from converged SCF wavefunction.
        """

        self.rsp_property = self.rsp_driver.compute(molecule, basis,
                                                    scf_tensors)

    def get_property(self, key):
        """Gets response property/spectroscopy.

        Gets response property/spectroscopy.

        Parameters
        ----------
        key
            The keyword for the property.

        Returns
        -------
        dict_value
            The property.
        """

        return None

    def print_property(self, ostream):
        """Prints response property/spectroscopy to output stream.

        Prints response property/spectroscopy to output stream.

        Parameters
        ----------
        ostream
            The output stream.
        """

        return
