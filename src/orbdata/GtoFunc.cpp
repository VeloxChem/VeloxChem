//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
