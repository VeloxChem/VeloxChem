#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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


def print_xc_reference(xcfun, ostream):
    """
    Prints libxc reference.

    :param xcfun:
        The XC functional.
    :param ostream:
        The output stream.
    """

    if isinstance(xcfun, XCFunctional):
        valstr = f'Using the {xcfun.get_func_label()} functional.'
        ostream.print_info(valstr)
        ostream.print_blank()
        printed_refs = []
        for ref in xcfun.get_functional_reference():
            if ref not in printed_refs:
                ostream.print_reference(ref)
                printed_refs.append(ref)
        ostream.print_blank()

        valstr = 'Using the Libxc library '
        valstr += f'(v{xcfun.get_libxc_version()}).'
        ostream.print_info(valstr)
        ostream.print_blank()
        ostream.print_reference(xcfun.get_libxc_reference())
        ostream.print_blank()

        valstr = 'Using the following algorithm for XC numerical integration.'
        ostream.print_info(valstr)
        ostream.print_blank()
        valstr = 'J. Kussmann, H. Laqua and C. Ochsenfeld, '
        valstr += 'J. Chem. Theory Comput. 2021, 17, 1512-1521'
        ostream.print_reference(valstr)
        ostream.print_blank()
