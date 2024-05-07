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

#include "GtoFunc.hpp"

namespace gtofunc {  // gtofunc namespace

auto
makeGtoBlocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoBlock>
{
    if (const auto mang = basis.getMaxAngularMomentum(); mang >= 0)
    {
        std::vector<CGtoBlock> gto_blocks;

        for (int64_t i = 0; i <= mang; i++)
        {
            for (const auto npgtos : basis.getContractionDepths(i))
            {
                gto_blocks.push_back(CGtoBlock(basis, molecule, i, npgtos));
            }
        }

        return gto_blocks;
    }
    else
    {
        return std::vector<CGtoBlock>();
    }
}

auto
makeGtoBlocks(const CMolecularBasis& basis, const CMolecule& molecule, const std::vector<int64_t>& atoms) -> std::vector<CGtoBlock>
{
    if (const auto mang = basis.getMaxAngularMomentum(atoms); mang >= 0)
    {
        std::vector<CGtoBlock> gto_blocks;

        for (int64_t i = 0; i <= mang; i++)
        {
            for (const auto npgtos : basis.getContractionDepths(atoms, i))
            {
                gto_blocks.push_back(CGtoBlock(basis, molecule, atoms, i, npgtos));
            }
        }

        return gto_blocks;
    }
    else
    {
        return std::vector<CGtoBlock>();
    }
}

auto
getNumberOfAtomicOrbitals(const std::vector<CGtoBlock>& gto_blocks) -> int64_t
{
    int64_t naos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();

        const auto ang = gto_block.getAngularMomentum();

        naos += ncgtos * (ang * 2 + 1);
    }

    return naos;
}

}  // namespace gtofunc
