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
