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

#ifndef BatchFunc_hpp
#define BatchFunc_hpp

#include <cstdint>
#include <utility>

namespace batch {  // batch namespace

/**
 Gets starting index of requested batch.

 @param ibatch the index of batch.
 @param nelements the number of elements in vector.
 @param nbatches the number of batches to partition vector.
 @return the starting index of batch.
 */
inline auto
getBatchIndex(const int64_t ibatch, const int64_t nelements, const int64_t nbatches) -> int64_t
{
    const auto bdim = nelements / nbatches;

    const auto brem = nelements % nbatches;

    if (const auto boff = bdim * ibatch; boff > nelements)
    {
        return nelements;
    }
    else
    {
        return boff + ((ibatch <= brem) ? ibatch : brem);
    }
};

/**
 Gets number of batches to partition vector.

 @param nelements the number of elements in vector.
 @param bsize the batch size used to partition vector.
 @return the number of batches.
 */
inline auto
getNumberOfBatches(const int64_t nelements, const int64_t bsize) -> int64_t
{
    const auto nbatches = nelements / bsize;

    if ((nelements % bsize) != 0)
    {
        return nbatches + 1;
    }
    else
    {
        return nbatches;
    }
};

/**
 Gets range of specific batch in partitioned vector.

 @param ibatch the index of batch.
 @param nelements the number of elements in vector.
 @param bsize the batch size used to partition vector.
 @return the number of loop passes.
 */
inline auto
getBatchRange(const int64_t ibatch, const int64_t nelements, const int64_t bsize) -> std::pair<int64_t, int64_t>
{
    const auto first = ibatch * bsize;

    if (first > nelements)
    {
        return {nelements, nelements};
    }

    const auto last = first + bsize;

    if (last > nelements)
    {
        return {first, nelements};
    }
    else
    {
        return {first, last};
    }
};

}  // namespace batch

#endif /* BatchFunc_hpp */
