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


def get_uff_parameters():
    """
    Initialize UFF parameters. Reference: A. K. Rappé, C. J. Casewit, K. S.
    Colwell, W. A. Goddard III, W. M. Skiff, J. Am. Chem. Soc. 1992, 114,
    10024-10035.

    :return:
        A dictionary containing UFF parameters.
    """

    periodic_table = [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
        'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
        'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
        'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
        'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
    ]

    # H-Zn
    uff_x_values = [
        2.886, 2.362, 2.451, 2.745, 4.083, 3.851, 3.660, 3.500, 3.364, 3.243,
        2.983, 3.021, 4.499, 4.295, 4.147, 4.035, 3.947, 3.868, 3.812, 3.399,
        3.295, 3.175, 3.144, 3.023, 2.961, 2.912, 2.872, 2.834, 3.495, 2.763
    ]

    # Ga-Nd
    uff_x_values += [
        4.383, 4.280, 4.230, 4.205, 4.189, 4.141, 4.114, 3.641, 3.345, 3.124,
        3.165, 3.052, 2.998, 2.963, 2.929, 2.899, 3.148, 2.848, 4.463, 4.392,
        4.420, 4.470, 4.500, 4.404, 4.517, 3.703, 3.522, 3.556, 3.606, 3.575
    ]

    # Pm-Th
    uff_x_values += [
        3.547, 3.520, 3.493, 3.368, 3.451, 3.428, 3.409, 3.391, 3.374, 3.355,
        3.640, 3.141, 3.170, 3.069, 2.954, 3.120, 2.840, 2.754, 3.293, 2.705,
        4.347, 4.297, 4.370, 4.709, 4.750, 4.765, 4.900, 3.677, 3.478, 3.396
    ]

    # Pa-No
    uff_x_values += [
        3.424, 3.395, 3.424, 3.424, 3.381, 3.326, 3.339, 3.313, 3.299, 3.286,
        3.274, 3.248
    ]

    # H-Zn
    uff_D_values = [
        0.044, 0.056, 0.025, 0.085, 0.180, 0.105, 0.069, 0.060, 0.050, 0.042,
        0.030, 0.111, 0.505, 0.402, 0.305, 0.274, 0.227, 0.185, 0.035, 0.238,
        0.019, 0.017, 0.016, 0.015, 0.013, 0.013, 0.014, 0.015, 0.005, 0.124
    ]

    # Ga-Nd
    uff_D_values += [
        0.415, 0.379, 0.309, 0.291, 0.251, 0.220, 0.040, 0.235, 0.072, 0.069,
        0.059, 0.056, 0.048, 0.056, 0.053, 0.048, 0.036, 0.228, 0.599, 0.567,
        0.449, 0.398, 0.339, 0.332, 0.045, 0.364, 0.017, 0.013, 0.010, 0.010
    ]

    # Pm-Th
    uff_D_values += [
        0.009, 0.008, 0.008, 0.009, 0.007, 0.007, 0.007, 0.007, 0.006, 0.228,
        0.041, 0.072, 0.081, 0.067, 0.066, 0.037, 0.073, 0.080, 0.039, 0.385,
        0.680, 0.663, 0.518, 0.325, 0.284, 0.248, 0.050, 0.404, 0.033, 0.026
    ]

    # Pa-No
    uff_D_values += [
        0.022, 0.022, 0.019, 0.016, 0.014, 0.013, 0.013, 0.013, 0.012, 0.012,
        0.011, 0.011
    ]

    uff_parameters = {}

    for idx, (x, D) in enumerate(zip(uff_x_values, uff_D_values)):
        key = periodic_table[idx]

        # rename C to cx, H to hx, etc.
        if key in ['C', 'H', 'O', 'N', 'S', 'P']:
            key = key.lower() + 'x'

        sigma = 0.1 * x / 2.0**(1.0 / 6.0)
        epsilon = D * 4.184

        uff_parameters[key] = {'sigma': sigma, 'epsilon': epsilon}

    return uff_parameters
