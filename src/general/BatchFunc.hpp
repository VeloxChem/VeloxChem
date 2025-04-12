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
