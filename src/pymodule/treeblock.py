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

import numpy as np

class TreeBlock:
    """
    Implements tree growth unit data storage for global optimization
    algorithm based on tree growth scheme.

    Instance variables
        - atoms: The list of atoms in tree growth unit.
        - bond_lengths: The list of bond lengths.
        - bond_angles: The list of bond angles.
        - ref_atoms: The list of reference atoms.
        - max_cnums: The list of maximum coordination numbers of
          reference atoms.
    """

    def __init__(self, data = None):
        """
        Initializes tree growth unit data storage.
        
        :param data:
            The list of stringd with tree block data.
        """
        
        self.atoms = []
        self.bond_legths = []
        self.bond_angles = []
        self.ref_atoms = []
        self.max_cnums = []

    def name(self):
        """
        Gets abbreviated chemical name of tree block.

        :return:
            The string with abbreviated chemical name.
        """
        
        name = ''
        for atom in self.atoms:
            name += atom
        return name
