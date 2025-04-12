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

#include "OpenMPFunc.hpp"

#include "BatchFunc.hpp"

namespace omp {  // omp namespace

auto
makeWorkGroup(const std::vector<CGtoBlock>& gto_blocks) -> TWorkGroup
{
    const auto ntasks = static_cast<int64_t>(omp::getNumberOfThreads());

    auto wgroups = TWorkGroup(ntasks, TGraph());

    if (const auto nblocks = static_cast<int64_t>(gto_blocks.size()); nblocks > 0)
    {
        for (int64_t i = 0; i < nblocks; i++)
        {
            const auto bra_size = gto_blocks[i].getNumberOfBasisFunctions();

            for (int64_t itask = 0; itask < ntasks; itask++)
            {
                const auto first = batch::getBatchIndex(itask, bra_size, ntasks);

                const auto last = batch::getBatchIndex(itask + 1, bra_size, ntasks);

                if (first != last) wgroups[itask].push_back(T4Index({i, i, first, last}));
            }

            for (int64_t j = i + 1; j < nblocks; j++)
            {
                const auto ket_size = gto_blocks[j].getNumberOfBasisFunctions();

                if (bra_size > ket_size)
                {
                    for (int64_t itask = 0; itask < ntasks; itask++)
                    {
                        const auto first = batch::getBatchIndex(itask, bra_size, ntasks);

                        const auto last = batch::getBatchIndex(itask + 1, bra_size, ntasks);

                        if (first != last) wgroups[itask].push_back(T4Index({i, j, first, last}));
                    }
                }
                else
                {
                    for (int64_t itask = 0; itask < ntasks; itask++)
                    {
                        const auto first = batch::getBatchIndex(itask, ket_size, ntasks);

                        const auto last = batch::getBatchIndex(itask + 1, ket_size, ntasks);

                        if (first != last) wgroups[itask].push_back(T4Index({j, i, first, last}));
                    }
                }
            }
        }
    }

    return wgroups;
}

auto
makeWorkGroup(const std::vector<CGtoBlock>& bra_gto_blocks, const std::vector<CGtoBlock>& ket_gto_blocks) -> TWorkGroup
{
    const auto ntasks = static_cast<int64_t>(omp::getNumberOfThreads());

    auto wgroups = TWorkGroup(ntasks, TGraph());

    if (const auto bra_nblocks = static_cast<int64_t>(bra_gto_blocks.size()); bra_nblocks > 0)
    {
        if (const auto ket_nblocks = static_cast<int64_t>(ket_gto_blocks.size()); ket_nblocks > 0)
        {
            for (int64_t i = 0; i < bra_nblocks; i++)
            {
                const auto bra_size = bra_gto_blocks[i].getNumberOfBasisFunctions();

                for (int64_t j = 0; j < ket_nblocks; j++)
                {
                    for (int64_t itask = 0; itask < ntasks; itask++)
                    {
                        const auto first = batch::getBatchIndex(itask, bra_size, ntasks);

                        const auto last = batch::getBatchIndex(itask + 1, bra_size, ntasks);

                        if (first != last) wgroups[itask].push_back(T4Index({i, j, first, last}));
                    }
                }
            }
        }
    }

    return wgroups;
}

}  // namespace omp
