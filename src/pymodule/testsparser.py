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

from pathlib import Path


def print_features():
    """
    Parses the tests and prints the features.
    """

    tests_path = Path(__file__).parent / 'tests'

    tests_files = sorted((f for f in tests_path.iterdir() if f.is_file()))

    for f in tests_files:
        with f.open('r') as f_test:
            for line in f_test:
                if line.strip().startswith('# vlxtag:'):
                    tags = line.strip().split(':')[1]
                    tags = tags.replace(',', ' ').split()
                    for tag in tags:
                        print(tag.upper() + ' ', end='')
                    print()
