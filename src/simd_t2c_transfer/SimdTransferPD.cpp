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


#include "SimdTransferPD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_pd(double *values, const size_t nvalues, CSimdMatrix &buffer,
               const CSimdMatrix &coordinates, const size_t sd, const size_t sf,
               const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = std::sqrt(3.0);
    const auto f_1 = 0.5 * std::sqrt(3.0);

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
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);
    const auto *sd_3 = buffer.data(sd + 3);
    const auto *sd_4 = buffer.data(sd + 4);
    const auto *sd_5 = buffer.data(sd + 5);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);
    const auto *sf_9 = buffer.data(sf + 9);

#pragma omp simd aligned(ab_y, sd_0, sd_1, sd_3, sd_4, sd_5, sf_1, sf_3, sf_6, sf_7, \
                         sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -f_0 * ab_y[k] * sd_1[k]
                 + f_0 * sf_3[k];

        g_1[k] = -f_0 * ab_y[k] * sd_4[k]
                 + f_0 * sf_7[k];

        g_2[k] = 0.5 * ab_y[k] * sd_0[k]
                 + 0.5 * ab_y[k] * sd_3[k]
                 - ab_y[k] * sd_5[k]
                 - 0.5 * sf_1[k]
                 - 0.5 * sf_6[k]
                 + sf_8[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sd_0, sd_1, sd_2, sd_3, sd_4, sf_1, sf_4, sf_6, \
                         sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_0 * ab_y[k] * sd_2[k]
                 + f_0 * sf_4[k];

        g_4[k] = -f_1 * ab_y[k] * sd_0[k]
                 + f_1 * ab_y[k] * sd_3[k]
                 + f_1 * sf_1[k]
                 - f_1 * sf_6[k];

        g_5[k] = -f_0 * ab_z[k] * sd_1[k]
                 + f_0 * sf_4[k];

        g_6[k] = -f_0 * ab_z[k] * sd_4[k]
                 + f_0 * sf_8[k];
    }

#pragma omp simd aligned(ab_x, ab_z, sd_0, sd_1, sd_2, sd_3, sd_5, sf_1, sf_2, sf_5, sf_7, \
                         sf_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = 0.5 * ab_z[k] * sd_0[k]
                 + 0.5 * ab_z[k] * sd_3[k]
                 - ab_z[k] * sd_5[k]
                 - 0.5 * sf_2[k]
                 - 0.5 * sf_7[k]
                 + sf_9[k];

        g_8[k] = -f_0 * ab_z[k] * sd_2[k]
                 + f_0 * sf_5[k];

        g_9[k] = -f_1 * ab_z[k] * sd_0[k]
                 + f_1 * ab_z[k] * sd_3[k]
                 + f_1 * sf_2[k]
                 - f_1 * sf_7[k];

        g_10[k] = -f_0 * ab_x[k] * sd_1[k]
                  + f_0 * sf_1[k];
    }

#pragma omp simd aligned(ab_x, sd_0, sd_2, sd_3, sd_4, sd_5, sf_0, sf_2, sf_3, sf_4, \
                         sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_0 * ab_x[k] * sd_4[k]
                  + f_0 * sf_4[k];

        g_12[k] = 0.5 * ab_x[k] * sd_0[k]
                  + 0.5 * ab_x[k] * sd_3[k]
                  - ab_x[k] * sd_5[k]
                  - 0.5 * sf_0[k]
                  - 0.5 * sf_3[k]
                  + sf_5[k];

        g_13[k] = -f_0 * ab_x[k] * sd_2[k]
                  + f_0 * sf_2[k];

        g_14[k] = -f_1 * ab_x[k] * sd_0[k]
                  + f_1 * ab_x[k] * sd_3[k]
                  + f_1 * sf_0[k]
                  - f_1 * sf_3[k];
    }
}

}  // namespace simdovl
