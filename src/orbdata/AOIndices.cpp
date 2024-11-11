//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "AOIndices.hpp"

#include <algorithm>

namespace aoindices {  // aoindices namespace

std::vector<std::vector<int>>
getDimerAOIndices(const CMolecule& mol_1, const CMolecule& mol_2, const CMolecularBasis& basis_1, const CMolecularBasis& basis_2)
{
    // get indices of atomic orbitals located on each molecule

    // angular momentum based AO ordering
    // S, P-1, P0, P+1, D-2, D-1, D0, D+1, D+2, ...

    // AO_type:  S S S S S S S S P-1 P-1 P0 P0 P+1 P+1 ...
    // Molecule: A A A A B B B B A   B   A  B  A   B   ...

    int max_angl_1 = basis_1.max_angular_momentum();

    int max_angl_2 = basis_2.max_angular_momentum();

    int max_angl = std::max(max_angl_1, max_angl_2);

    std::vector<std::vector<int>> aoinds_mol;

    aoinds_mol.push_back(std::vector<int>());

    aoinds_mol.push_back(std::vector<int>());

    for (int aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        int nao_1 = basis_1.number_of_basis_functions(angl);

        int nao_2 = basis_2.number_of_basis_functions(angl);

        for (int s = -angl; s <= angl; s++)
        {
            for (int i = 0; i < nao_1; i++, aoidx++)
            {
                aoinds_mol[0].push_back(aoidx);
            }

            for (int i = 0; i < nao_2; i++, aoidx++)
            {
                aoinds_mol[1].push_back(aoidx);
            }
        }
    }

    return aoinds_mol;
}

}  // namespace aoindices
