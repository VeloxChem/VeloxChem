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


#include "SimdTransformPP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_pp(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t pp,
             const size_t nmax) -> void
{
    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_1 = buffer.data(pp + 1);
    const auto *pp_2 = buffer.data(pp + 2);
    const auto *pp_3 = buffer.data(pp + 3);
    const auto *pp_4 = buffer.data(pp + 4);
    const auto *pp_5 = buffer.data(pp + 5);
    const auto *pp_6 = buffer.data(pp + 6);
    const auto *pp_7 = buffer.data(pp + 7);
    const auto *pp_8 = buffer.data(pp + 8);

#pragma omp simd aligned(pp_1, pp_2, pp_3, pp_4, pp_5, pp_6, pp_7, \
                         pp_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = pp_4[k];

        g_1[k] = pp_5[k];

        g_2[k] = pp_3[k];

        g_3[k] = pp_7[k];

        g_4[k] = pp_8[k];

        g_5[k] = pp_6[k];

        g_6[k] = pp_1[k];

        g_7[k] = pp_2[k];
    }

#pragma omp simd aligned(pp_0 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = pp_0[k];
    }
}

auto
transform_pp_tri(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t pp,
                 const size_t nmax) -> void
{
    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_3 = buffer.data(pp + 3);
    const auto *pp_4 = buffer.data(pp + 4);
    const auto *pp_5 = buffer.data(pp + 5);
    const auto *pp_6 = buffer.data(pp + 6);
    const auto *pp_8 = buffer.data(pp + 8);

#pragma omp simd aligned(pp_0, pp_3, pp_4, pp_5, pp_6, pp_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = pp_4[k];

        g_1[k] = pp_5[k];
        g_3[k] = g_1[k];

        g_2[k] = pp_3[k];
        g_6[k] = g_2[k];

        g_4[k] = pp_8[k];

        g_5[k] = pp_6[k];
        g_7[k] = g_5[k];

        g_8[k] = pp_0[k];
    }
}

}  // namespace simdtrf
