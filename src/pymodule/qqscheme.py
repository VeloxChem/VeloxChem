from .veloxchemlib import ericut


def get_qq_type(qq_type):
    """
    Gets string with type of electron repulsion integrals screening scheme
    (Cauchy Schwarz and it's variations).

    :param qq_type:
        The label of electron repulsion integrals screening scheme.

    :return:
        The string with type of electron repulsion integrals screening
        scheme.
    """

    if qq_type == "QQ":
        return "Cauchy Schwarz"

    if qq_type == "QQR":
        return "Distance Dependent Cauchy Schwarz"

    if qq_type == "QQ_DEN":
        return "Cauchy Schwarz + Density"

    if qq_type == "QQR_DEN":
        return "Distance Dependent Cauchy Schwarz + Density"

    return "Undefined"


def get_qq_scheme(qq_type):
    """
    Converts screening scheme string to C++ enum.

    :param qq_type:
        The label of electron repulsion integrals screening scheme.

    :return:
        The C++ enum with screening scheme.
    """

    if qq_type == "QQ":
        return ericut.qq

    if qq_type == "QQR":
        return ericut.qqr

    if qq_type == "QQ_DEN":
        return ericut.qqden

    if qq_type == "QQR_DEN":
        return ericut.qqrden

    return None
