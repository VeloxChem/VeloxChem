from .rspdriver import ResponseDriver


class ResponseProperty:
    """
    Implements the base class for response property/spectroscopy.

    :param rsp_dict:
        The input dictionary that defines the property/spectroscopy.
    :param method_dict:
        The dictionary of method settings.

    Instance variables
        - rsp_dict: The dictionary of response input.
        - method_dict: The dictionary of method settings.
        - rsp_driver: The response driver.
        - rsp_property: The dictionary of response property.
    """

    def __init__(self, rsp_dict, method_dict=None):
        """
        Initializes response property/spectroscopy.
        """

        self.rsp_dict = rsp_dict
        self.method_dict = method_dict

    def init_driver(self, comm, ostream):
        """
        Initializes response driver.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        self.rsp_driver = ResponseDriver(comm, ostream)
        self.rsp_driver.update_settings(self.rsp_dict, self.method_dict)

    def compute(self, molecule, basis, scf_tensors):
        """
        Computes response property/spectroscopy.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        """

        self.rsp_property = self.rsp_driver.compute(molecule, basis,
                                                    scf_tensors)

    def converged(self):
        """
        Checks if the response calculation is converged.

        :return:
            True if the response calculation is converged, False otherwise.
        """

        return self.rsp_driver.is_converged

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
