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

from pathlib import Path


def _get_list_of_features():
    """
    Gets the list of features.

    :return:
        The list of features.
    """

    list_of_features = []

    vlx_py_path = Path(__file__).parent / 'tests'
    vlx_py_files = sorted((f for f in vlx_py_path.iterdir() if f.is_file()))

    for f in vlx_py_files:
        with f.open('r') as fh:
            for line in fh:
                if line.strip().startswith('# vlxtag:'):
                    tags = line.strip().split(':')[1]
                    tags = tuple(tags.replace(',', ' ').split())
                    if tags not in list_of_features:
                        list_of_features.append(tags)

    return list_of_features


def _print_list_of_features(list_of_features):
    """
    Prints a list of features.
    """

    sorted_features = sorted(list_of_features, key=lambda f: f[1])

    for feature in sorted_features:
        for tag in feature:
            print(f'  {tag.upper()}', end='')
        print()


def _known_aliases_for_keywords():
    """
    Gets the aliased keyword.

    :param keyword:
        The input keyword.

    :return:
        The mapped keyword.
    """

    return {
        'hf': ('rhf', 'uhf', 'rohf'),
        'dft': ('rks', 'uks', 'roks'),
        'rsp': ('lr', 'qr', 'cr'),
        'response': ('lr', 'qr', 'cr'),
        'linear response': 'lr',
        'linear-response': 'lr',
        'quadratic response': 'qr',
        'quadratic-response': 'qr',
        'cubic response': 'cr',
        'cubic-response': 'cr',
        'td': ('cis', 'tdhf', 'tda', 'tddft'),
        'rpa': 'tdhf',
        'uv-vis': 'absorption',
        'cd': ('ecd',),
        'spectrum': ('absorption', 'ecd'),
        'mm force field': 'mm_force_field_generation',
        'resp': 'resp_charges',
        'esp': ('esp_charges', 'esp_on_points'),
    }


def _found_tag_in_keyword(tag, keyword):
    """
    Checks if a tag is associated with a keyword.
    """

    known_aliases = _known_aliases_for_keywords()

    if keyword.lower() in known_aliases:
        aliased_keywords = known_aliases[keyword.lower()]
    else:
        aliased_keywords = keyword

    if (isinstance(aliased_keywords, str) and
            tag.lower() == aliased_keywords.lower()):
        return True

    if (isinstance(aliased_keywords, tuple) and
            tag.lower() in aliased_keywords):
        return True

    return False


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
            for tag in feature:
                if (_found_tag_in_keyword(tag, keyword) and
                        feature not in available_features):
                    available_features.append(feature)
                    break

        if available_features:
            print(f'Available features for keyword \"{keyword}\":')
            _print_list_of_features(available_features)

        else:
            print(f'No available features for keyword \"{keyword}\".\n')
            print('Please check the features with other keywords:')
            _print_list_of_features(all_features)
