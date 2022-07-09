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


def _get_list_of_features():
    """
    Gets the list of features.

    :return:
        The list of features.
    """

    list_of_features = []

    tests_path = Path(__file__).parent / 'tests'
    tests_files = sorted((f for f in tests_path.iterdir() if f.is_file()))

    for f in tests_files:
        with f.open('r') as f_test:
            for line in f_test:
                if line.strip().startswith('# vlxtag:'):
                    tags = line.strip().split(':')[1]
                    tags = tags.replace(',', ' ').split()
                    list_of_features.append(tuple(tags))

    return list_of_features


def _print_list_of_features(list_of_features):
    """
    Prints a list of features.
    """

    for feature in list_of_features:
        for tag in feature:
            print(f'  {tag.upper()}', end='')
        print()


def print_features(keyword=None):
    """
    Prints the features related to keyword.
    """

    all_features = _get_list_of_features()

    if keyword is None:
        _print_list_of_features(all_features)

    else:
        available_features = []

        for feature in all_features:
            if keyword.lower() == feature[1].lower():
                available_features.append(feature)

        if available_features:
            print(f'Available features for keyword \"{keyword}\":')
            _print_list_of_features(available_features)

        else:
            print(f'No available features for keyword \"{keyword}\".\n')
            print('Please check the features with other keywords:')
            _print_list_of_features(all_features)
