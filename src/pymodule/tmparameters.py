#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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


def get_tm_parameters():
    """
    Initialize transition metal (TM) parameters. Reference: F. Šebesta, V. Sláma, J. Melcr, 
    Z. Futera, and J. V. Burda. J. Chem. Theory Comput. 2016 12 (8), 3681-3688.

    :return:
        A dictionary containing TM parameters. When multiple oxidation states available in the reference, 
        parameters are chosen for the lower oxidation state. 
    """

    tm_elements = [
        'Sc','Ti','V','Cr','Mn',
        'Fe','Co','Ni','Cu','Zn',
        'Ru','Rh','Pt','Hg'
    ]   
    
    sigma_values = [
        3.330, 2.869, 2.767, 2.734, 2.904,
        2.969, 2.805, 2.720, 2.668, 2.847, 
        2.939, 2.715, 2.670, 2.801
    ]
    
    epsilon_values = [
        0.110, 0.917, 1.904, 1.518, 1.089,
        0.698, 1.196, 2.650, 2.148, 1.095, 
        0.418, 3.946, 5.073, 1.955
    ]

    tm_parameters = {}

    for idx, (s, e) in enumerate(zip(sigma_values, epsilon_values)):
        key = tm_elements[idx]

        sigma = s * 2.0 # Å
        epsilon = e * 4.184  # kJ/mol

        tm_parameters[key] = {'sigma': sigma, 'epsilon': epsilon}


    return tm_parameters