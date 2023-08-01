import numpy as np

from .veloxchemlib import SubMatrix

def _SubMatrix_deepcopy(self, memo):
    """
    Implements deepcopy.

    :param memo:
        The memo dictionary for deepcopy.

    :return:
        A deepcopy of self.
    """

    return SubMatrix(self)

SubMatrix.__deepcopy__ = _SubMatrix_deepcopy
