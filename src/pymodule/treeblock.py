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
        
        # initialize lists
        self.atoms = []
        self.ref_atoms = []
        self.bond_lengths = []
        self.bond_angles = []
        self.max_cnums = []
        
        if data is not None:
            for line in data:
                keys = line.split()
                self.atoms.append(keys[0].upper())
                self.ref_atoms.append(keys[3].upper())
                self.bond_lengths.append(float(keys[1]))
                self.bond_angles.append(float(keys[2]))
                self.max_cnums.append(int(keys[4]))
                
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
        
    def number_of_atoms(self):
        """
        Gets number of atoms in tree block.

        :return:
            The number of atoms in tree block.
        """
        
        return len(self.atoms)
        
    def get_atom(self, iatom):
        """
        Gets capitalized label of specific atom in tree block.

        :param iatom:
            The index of atom in tree block.
        :return:
            The capatilized label of specific atom.
        """
        
        return self.atoms[iatom].upper()
        
    def get_ref_atom(self, iatom):
        """
        Gets capitalized label of specific reference atom in tree block.

        :param iatom:
            The index of reference atom in tree block.
        :return:
            The capatilized label of specific reference atom.
        """
    
        return self.ref_atoms[iatom].upper()
        
    def get_bond_length(self, ibond):
        """
        Gets specific bond length in tree block.

        :param ibond:
            The index of bond length in tree block.
        :return:
            The specific bond length.
        """
    
        return self.bond_lengths[ibond]
        
    def get_bond_angle(self, iangle):
        """
        Gets specific bond angle in tree block.

        :param iangle:
            The index of bond angle in tree block.
        :return:
            The specific bond angle.
        """
    
        return self.bond_angles[iangle]
        
    def get_coord_number(self, iatom):
        """
        Gets coordination number of specific reference atom in tree block.

        :param iatom:
            The index of reference atom in tree block.
        :return:
            The coordination number of specific reference atom.
        """
        
        return self.max_cnums[iatom]
    
    def geom_string(self):
        """
        Gets geometry string of tree block.

        :return:
            The geometry string of tree block.
        """
        
        width = 70
        geom_info = []

        geom_info.append('Geometry of Tree Growth Block:'.ljust(width))
        geom_info.append('=============================='.ljust(width))
        geom_info.append('                              '.ljust(width))
        geom_info.append('* Atom * Bond Length * Bond Angle * Reference * Max. Coordination *'.ljust(width))
        for i in range(self.number_of_atoms()):
            curr_str = '*' + self.get_atom(i).center(6)
            curr_str += '*' + f'{self.bond_lengths[i]:.2f}'.center(13)
            curr_str += '*' + f'{self.bond_angles[i]:.2f}'.center(12)
            curr_str += '*' + self.get_ref_atom(i).center(11)
            curr_str += '*' + f'{self.max_cnums[i]}'.center(19) + '*'
            geom_info.append(curr_str.ljust(width))
    
        return '\n'.join(geom_info)


