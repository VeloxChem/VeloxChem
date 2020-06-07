from .rspproperty import ResponseProperty
from .inputparser import parse_frequencies


class TPA(ResponseProperty):
    """
    Implements the cubic Two-photon abs property.

    :param rsp_dict:
        The dictionary of response input.
    :param method_dict:
        The dictionary of method settings.

    Instance variables
        - rsp_dict: The dictionary of response input.
        - method_dict: The dictionary of method settings.
        - rsp_property: The dictionary of response property.
    """

    def __init__(self, rsp_dict, method_dict={}):
        """
        Initializes the polarizability property.
        """

        rsp_dict = dict(rsp_dict)
        method_dict = dict(method_dict)

        rsp_dict['property'] = 'tpa'
        rsp_dict['response'] = 'cubic'
        rsp_dict['residue'] = 'none'
        rsp_dict['complex'] = 'yes'


        if 'frequencies' not in rsp_dict:
            rsp_dict['frequencies'] = '0'

        super().__init__(rsp_dict, method_dict)

    def get_property(self, key):
        """
        Gets response functions or solutions.

        :param key:
            The keyword 'response_functions' or 'solutions'.

        :return:
            The response functions or solutions.
        """

        return self.rsp_property[key]
