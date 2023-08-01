from .veloxchemlib import Matrix

def _Matrix_deepcopy(self, memo):
    """
    Implements deepcopy.

    :param memo:
        The memo dictionary for deepcopy.

    :return:
        A deepcopy of self.
    """

    return Matrix(self)

Matrix.__deepcopy__ = _Matrix_deepcopy
