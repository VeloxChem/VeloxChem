#
#                           VELOXCHEM 1.0-RC
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

from pathlib import Path

from .veloxchemlib import XTBDriver


def _XTBDriver_compute(self, molecule, ostream):
    """
    Computes DFT-B energy using XTB driver.

    :param molecule:
        The molecule.
    :param ostream:
        The output stream.
    """

    ostream.print_blank()
    ostream.print_header('XTB Driver')
    ostream.print_header(12 * '=')
    ostream.flush()

    self._compute(molecule)

    if self.is_master_node():
        for line in self.get_output():
            ostream.print_line(line)
        ostream.flush()
        output = Path(self.get_output_filename())
        if output.is_file():
            output.unlink()


XTBDriver.compute = _XTBDriver_compute
