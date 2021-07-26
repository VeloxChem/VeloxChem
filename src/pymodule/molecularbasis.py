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

from os import environ
from pathlib import Path

from .errorhandler import assert_msg_critical
from .inputparser import InputParser
from .outputstream import OutputStream
from .veloxchemlib import (AtomBasis, BasisFunction, ChemicalElement,
                           MolecularBasis, to_angular_momentum)


def _known_aliases_for_basis_sets():
    """
    Maps from common basis name to basis filename.

    :return:
        A dictionary that maps from common basis name to basis filename.
    """

    return {
        '6-31G*': '6-31G_D_',
        '6-31G**': '6-31G_D,P_',
        '6-31+G*': '6-31+G_D_',
        '6-31+G**': '6-31+G_D,P_',
        '6-31++G*': '6-31++G_D_',
        '6-31++G**': '6-31++G_D,P_',
        '6-311G*': '6-311G_D_',
        '6-311G**': '6-311G_D,P_',
        '6-311+G*': '6-311+G_D_',
        '6-311+G**': '6-311+G_D,P_',
        '6-311++G*': '6-311++G_D_',
        '6-311++G**': '6-311++G_D,P_',
        '6-31G(2DF,P)': '6-31G_2DF,P_',
        '6-31G(3DF,3PD)': '6-31G_3DF,3PD_',
        '6-311G(2DF,2PD)': '6-311G_2DF,2PD_',
        '6-311+G(2D,P)': '6-311+G_2D,P_',
        '6-311++G(2D,2P)': '6-311++G_2D,2P_',
        '6-311++G(3DF,3PD)': '6-311++G_3DF,3PD_',
        'DEF2-SV(P)': 'DEF2-SV_P_',
    }


def _basis_name_to_file(name):
    """
    Determines basis set file name from common name.

    :param name:
        Common name of the basis set.

    :return:
        The file name of the basis set.
    """

    known_aliases = _known_aliases_for_basis_sets()

    if name in known_aliases:
        return known_aliases[name]
    else:
        return name


def _basis_file_to_name(fname):
    """
    Determines basis set common name from file name.

    :param fname:
        Name of the basis set file.

    :return:
        The common name of the basis set.
    """

    common_names = []

    for key, val in _known_aliases_for_basis_sets().items():
        if str(fname) == val and key not in common_names:
            common_names.append(key)

    if common_names:
        return sorted(common_names)[0]
    else:
        return str(fname)


@staticmethod
def _MolecularBasis_read(mol, basis_name, basis_path='.', ostream=None):
    """
    Reads AO basis set from file.

    :param mol:
        The molecule.
    :param basis_name:
        Name of the basis set.
    :param basis_path:
        Path to the basis set.
    :param ostream:
        The output stream.

    :return:
        The AO basis set.
    """

    if ostream is None:
        ostream = OutputStream(None)

    err_gc = 'MolcularBasis.read: '
    err_gc += 'General contraction currently is not supported'

    # de-alias basis set name to basis set file
    fname = _basis_name_to_file(basis_name.upper())

    # searching order:
    # 1. given basis_path
    # 2. current directory
    # 3. VLXBASISPATH

    basis_file = Path(basis_path, fname)

    if not basis_file.is_file() and basis_path != '.':
        basis_file = Path('.', fname)

    if not basis_file.is_file() and 'VLXBASISPATH' in environ:
        basis_file = Path(environ['VLXBASISPATH'], fname)

    basis_info = 'Reading basis set from file: ' + str(basis_file)
    ostream.print_info(basis_info)
    ostream.print_blank()

    basis_dict = InputParser(str(basis_file)).input_dict

    assert_msg_critical(
        fname.upper() == basis_dict['basis_set_name'].upper(),
        'MolecularBasis.read: Inconsistent basis set name',
    )

    mol_basis = MolecularBasis()

    elem_comp = mol.get_elemental_composition()

    for elem_id in elem_comp:

        elem = ChemicalElement()
        err = elem.set_atom_type(elem_id)
        assert_msg_critical(err, 'ChemicalElement.set_atom_type')

        basis_key = 'atombasis_{}'.format(elem.get_name().lower())
        basis_list = [entry for entry in basis_dict[basis_key]]

        atom_basis = AtomBasis()

        while basis_list:
            shell_title = basis_list.pop(0).split()
            assert_msg_critical(
                len(shell_title) == 3,
                'Basis set parser (shell): {}'.format(' '.join(shell_title)),
            )

            angl = to_angular_momentum(shell_title[0])
            npgto = int(shell_title[1])
            ncgto = int(shell_title[2])

            assert_msg_critical(ncgto == 1, err_gc)

            expons = [0.0] * npgto
            coeffs = [0.0] * npgto * ncgto

            for i in range(npgto):
                prims = basis_list.pop(0).split()
                assert_msg_critical(
                    len(prims) == ncgto + 1,
                    'Basis set parser (primitive): {}'.format(' '.join(prims)),
                )

                expons[i] = float(prims[0])
                for k in range(ncgto):
                    coeffs[k * npgto + i] = float(prims[k + 1])

            bf = BasisFunction(expons, coeffs, angl)
            bf.normalize()

            atom_basis.add_basis_function(bf)

        atom_basis.set_elemental_id(elem_id)

        mol_basis.add_atom_basis(atom_basis)

    mol_basis.set_label(basis_name.upper())

    return mol_basis


@staticmethod
def _MolecularBasis_get_avail_basis(element_label):
    """
    Gets the names of available basis sets for an element.

    :param element_label:
        The label of the chemical element.

    :return:
        The tuple of basis sets.
    """

    avail_basis = set()

    basis_path = Path(environ['VLXBASISPATH'])
    basis_files = sorted((x for x in basis_path.iterdir() if x.is_file()))

    for x in basis_files:
        name = _basis_file_to_name(x.name)
        basis = InputParser(str(x)).input_dict
        # check that the given element appears as key
        # and that its value is a non-empty list
        elem = f'atombasis_{element_label.lower()}'
        if elem in basis.keys():
            if basis[elem]:
                avail_basis.add(name)

    return sorted(list(avail_basis))


MolecularBasis.read = _MolecularBasis_read
MolecularBasis.get_avail_basis = _MolecularBasis_get_avail_basis
