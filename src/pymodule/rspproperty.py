from .rspdriver import ResponseDriver


class ResponseProperty:
    """
    Implements the base class for response property/spectroscopy.

    :param rsp_input:
        The dictionary of response input.
    :param rsp_driver:
        The response driver.
    :param rsp_property:
        The dictionary of response property.
    """

    def __init__(self, rsp_input, method_input={}):
        """
        Initializes response property/spectroscopy.

        :param rsp_input:
            The input dictionary that defines the property/spectroscopy.
        """

        self.rsp_input = rsp_input
        self.method_input = method_input

    def init_driver(self, comm, ostream):
        """
        Initializes response driver.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        self.rsp_driver = ResponseDriver(comm, ostream)
        self.rsp_driver.update_settings(self.rsp_input, self.method_input)

    def compute(self, molecule, basis, scf_tensors):
        """
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        """

        self.rsp_property = self.rsp_driver.compute(molecule, basis,
                                                    scf_tensors)

    def get_property(self, key):
        """
        Gets response property/spectroscopy.

        :param key:
            The keyword for the property.

        :return:
            The property.
        """

        return None

    def print_property(self, ostream):
        """
        Prints response property/spectroscopy to output stream.

        :param ostream:
            The output stream.
        """

        return
