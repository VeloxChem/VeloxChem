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

#ifndef BatchFunc_hpp
#define BatchFunc_hpp

#include <utility>

#include "CustomConstrains.hpp"

namespace batch {  // batch

/// @brief Gets number of batches to partition vector.
/// @param nelements The number of elements in vector.
/// @param bsize The batch size used to partition vector.
/// @return The number of batches.
template <Integral T>
inline auto
number_of_batches(const T nelements, const T bsize) -> T
{
    const auto nbatches = nelements / bsize;

    return ((nelements % bsize) != 0) ? nbatches + 1 : nbatches;
}

/// @brief Gets range of specific batch in partitioned vector.
/// @param ibatch The index of batch.
/// @param nelements The number of elements in vector.
/// @param bsize The batch size used to partition vector.
/// @return The  [first, last) pair  for requested batch.
template <Integral T>
inline auto
batch_range(const T ibatch, const T nelements, const T bsize) -> std::pair<T, T>
{
    const auto first = ibatch * bsize;

    if (first > nelements) return {nelements, nelements};

    const auto last = first + bsize;

    if (last > nelements)
    {
        return {first, nelements};
    }
    else
    {
        return {first, last};
    }
}

/// @brief Gets range of specific batch in partitioned vector.
/// @param ibatch The index of batch.
/// @param nelements The number of elements in vector.
/// @param bsize The batch size used to partition vector.
/// @param position The shift of pair positions.
/// @return The  [first, last) pair  for requested batch.
template <Integral T>
inline auto
batch_range(const T ibatch, const T nelements, const T bsize, const T position) -> std::pair<T, T>
{
    auto range = batch::batch_range(ibatch, nelements, bsize);

    return {range.first + position, range.second + position};
}

}  // namespace batch

#endif /* BatchFunc_hpp */
