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

#include "GtoPairBlockFunc.hpp"

#include <ranges>

#include "CustomViews.hpp"
#include "GtoFunc.hpp"

namespace gtofunc {  // gtofunc namespace

auto
make_gto_pair_blocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoPairBlock>
{
    return gtofunc::make_gto_pair_blocks(gtofunc::make_gto_blocks(basis, molecule));
}

auto
make_gto_pair_blocks(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<CGtoPairBlock>
{
    if (const auto nblocks = gto_blocks.size(); nblocks > 0)
    {
        std::vector<CGtoPairBlock> gpblocks;

        gpblocks.reserve(nblocks * (nblocks + 1) / 2);

        std::ranges::for_each(views::triangular(nblocks), [&](const auto& index) {
            const auto [i, j] = index;
            if (i == j)
            {
                gpblocks.push_back(CGtoPairBlock(gto_blocks[i]));
            }
            else
            {
                gpblocks.push_back(CGtoPairBlock(gto_blocks[i], gto_blocks[j]));
            }
        });

        return gpblocks;
    }
    else
    {
        return std::vector<CGtoPairBlock>();
    }
}

auto
make_gto_pair_blocks(const std::vector<CGtoBlock>& bra_gto_blocks,
                     const std::vector<CGtoBlock>& ket_gto_blocks) -> std::vector<CGtoPairBlock>
{
    const auto nbra_blocks = bra_gto_blocks.size();
    
    const auto nket_blocks = ket_gto_blocks.size();
    
    if  ((nbra_blocks > 0) && (nket_blocks > 0))
    {
        std::vector<CGtoPairBlock> gpblocks;
            
        std::ranges::for_each(views::rectangular(nbra_blocks, nket_blocks), [&](const auto& index) {
            const auto [i, j] = index;
            gpblocks.push_back(CGtoPairBlock(bra_gto_blocks[i], ket_gto_blocks[j]));
        });
            
        return gpblocks;
    }
    else
    {
        return std::vector<CGtoPairBlock>();
    }
}

}  // namespace gtofunc
