#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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

from .veloxchemlib import xcfun
from .veloxchemlib import parse_xc_func
from .veloxchemlib import XCFunctional
from .errorhandler import assert_msg_critical


def get_default_grid_level(xc_func):
    """
    Gets default grid level for an exchange-correlation functional.

    :param xc_func:
        The exchange-correlation functional.

    :return:
        The default grid level.
    """

    if isinstance(xc_func, str):
        xc_func_obj = parse_xc_func(xc_func.upper())
        func_name = xc_func_obj.get_func_label()
        func_type = xc_func_obj.get_func_type()
    else:
        func_name = xc_func.get_func_label()
        func_type = xc_func.get_func_type()

    # LDA

    if func_type == xcfun.lda:

        return 4

    # GGA

    elif func_type == xcfun.gga:

        if func_name.upper() in [
                'B97',
                'B97-1',
                'B97-2',
                'B97-3',
        ]:
            return 5

        else:
            return 4

    # meta-GGA

    elif func_type == xcfun.mgga:

        if func_name.upper() in [
                'SCAN',
        ]:
            return 7

        elif func_name.upper() in [
                'RSCAN',
                'R2SCAN',
        ]:
            return 6

        elif func_name.upper() in [
                'M05',
                'M05-2X',
                'M06',
                'M06-2X',
                'M06-HF',
                'M06-L',
                'M11-L',
        ]:
            return 6

        else:
            return 5

    else:
        assert_msg_critical(
            False, 'get_default_grid_level: Invalid XC functional type')


def print_libxc_reference(xcfun, ostream):
    """
    Prints libxc reference.

    :param xcfun:
        The XC functional.
    :param ostream:
        The output stream.
    """

    if isinstance(xcfun, XCFunctional):
        valstr = 'Using the Libxc library '
        valstr += f'(version {xcfun.get_libxc_version()}).'
        ostream.print_info(valstr)
        ostream.print_blank()
        valstr = xcfun.get_libxc_reference()
        ostream.print_reference(valstr)
        ostream.print_blank()

        valstr = f'Using the {xcfun.get_func_label()} functional.'
        ostream.print_info(valstr)
        ostream.print_blank()
        printed_refs = []
        for ref in xcfun.get_functional_reference():
            if ref not in printed_refs:
                ostream.print_reference(ref)
                printed_refs.append(ref)
        ostream.print_blank()
