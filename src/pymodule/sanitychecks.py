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

from .veloxchemlib import parse_xc_func
from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical


def molecule_sanity_check(mol):
    """
    Checks molecule for charge/multiplicity combination and geometry.

    :param mol:
        The molecule.
    """

    mol.check_multiplicity()
    mol.check_proximity(0.1)


def dft_sanity_check(obj, method_flag='compute', response_flag='none'):
    """
    Checks DFT settings and updates relevant attributes.

    :param obj:
        The object (SCF or response driver) that is being checked.
    :param method_flag:
        The flag indicating the method in which the sanity check is
        called.
    :param response_flag:
        The flag indicating the type of response calculation in which the
        sanity check is called.
    """

    xcfun_is_none = (obj.xcfun is None)
    xcfun_is_hf = (isinstance(obj.xcfun, str) and obj.xcfun.lower() == 'hf')

    # Hartree-Fock: xcfun is None or 'hf'
    if xcfun_is_none or xcfun_is_hf:
        obj._dft = False

    # DFT: xcfun is functional object or string (other than 'hf')
    else:
        if isinstance(obj.xcfun, str):
            obj.xcfun = parse_xc_func(obj.xcfun.upper())
        assert_msg_critical(not obj.xcfun.is_undefined(),
                            f'{type(obj).__name__}: Undefined XC functional')
        obj._dft = True

    # check grid level
    if obj._dft and obj.grid_level is not None:
        if (obj.grid_level < 1 or obj.grid_level > 8):
            warn_msg = f'Invalid DFT grid level {obj.grid_level}. '
            warn_msg += 'Using default value.'
            obj.ostream.print_warning(warn_msg)
            obj.grid_level = None
        elif (method_flag.lower() == 'compute' and
              obj.grid_level < get_default_grid_level(obj.xcfun)):
            warn_msg = 'DFT grid level is below the recommended value. '
            warn_msg += 'Please double check.'
            obj.ostream.print_warning(warn_msg)
        obj.ostream.flush()
