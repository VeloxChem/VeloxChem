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

import numpy as np
import time as tm
import sys

from .veloxchemlib import XtbDriver
from .errorhandler import assert_msg_critical


def _XtbDriver_compute(self, molecule, ostream=None):
    """
    Computes DFT-B energy using XTB driver.

    :param molecule:
        The molecule.
    :param ostream:
        The output stream.
    """

    errmsg = 'XtbDriver: XTB not available. Please download and install XTB '
    errmsg += 'from https://github.com/grimme-lab/xtb, set XTBHOME environment '
    errmsg += 'variable, and reinstall VeloxChem.'
    assert_msg_critical(self.is_available(), errmsg)

    start_time = tm.time()

    if self.is_master_node() and ostream is not None:
        ostream.print_blank()
        ostream.print_header('XTB Driver')
        ostream.print_header(12 * '=')
        ostream.print_blank()
        ostream.flush()

    self._compute(molecule)
    sys.stdout.flush()

    if self.is_master_node():
        energy = self.get_energy()
        gradient = self.get_gradient()

        grad2 = np.sum(gradient**2, axis=1)
        rms_grad = np.sqrt(np.mean(grad2))
        max_grad = np.max(np.sqrt(grad2))

        if ostream is not None:
            valstr = '  Energy   : {:.10f} a.u.'.format(energy)
            ostream.print_info(valstr)

            valstr = '  Gradient : {:.6e} a.u. (RMS)'.format(rms_grad)
            ostream.print_info(valstr)
            valstr = '             {:.6e} a.u. (Max)'.format(max_grad)
            ostream.print_info(valstr)

            valstr = '  Time     : {:.2f} sec'.format(tm.time() - start_time)
            ostream.print_info(valstr)
            ostream.print_blank()

            valstr = 'Reference: C. Bannwarth, E. Caldeweyher, S. Ehlert,'
            ostream.print_info(valstr)
            valstr = 'A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme,'
            ostream.print_info(valstr)
            valstr = 'WIREs Comput. Mol. Sci., 2020, 11, e01493'
            ostream.print_info(valstr)
            ostream.print_blank()

            ostream.flush()


def _XtbDriver_get_reference(self):

    ref_str = 'C. Bannwarth, E. Caldeweyher, S. Ehlert, '
    ref_str += 'A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme, '
    ref_str += 'WIREs Comput. Mol. Sci., 2020, 11, e01493'
    return ref_str


XtbDriver.compute = _XtbDriver_compute
XtbDriver.get_reference = _XtbDriver_get_reference
