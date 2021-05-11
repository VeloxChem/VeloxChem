#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
