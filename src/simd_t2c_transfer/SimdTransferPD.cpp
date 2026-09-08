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

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_pd_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t sd, const size_t sf,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

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

auto
compute_hrr_pd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t sd, const size_t sf, const size_t nmax) -> void
{
    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, sd_0, sd_1, sd_2, sd_3, sd_4, sf_0, \
                         sf_1, sf_2, sf_3, sf_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * sd_0[k]
                 + sf_0[k];

        t_1[k] = -ab_x[k] * sd_1[k]
                 + sf_1[k];

        t_2[k] = -ab_x[k] * sd_2[k]
                 + sf_2[k];

        t_3[k] = -ab_x[k] * sd_3[k]
                 + sf_3[k];

        t_4[k] = -ab_x[k] * sd_4[k]
                 + sf_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_y, sd_0, sd_1, sd_2, sd_5, sf_1, sf_3, \
                         sf_4, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * sd_5[k]
                 + sf_5[k];

        t_6[k] = -ab_y[k] * sd_0[k]
                 + sf_1[k];

        t_7[k] = -ab_y[k] * sd_1[k]
                 + sf_3[k];

        t_8[k] = -ab_y[k] * sd_2[k]
                 + sf_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_y, ab_z, sd_0, sd_3, sd_4, sd_5, sf_2, \
                         sf_6, sf_7, sf_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_9[k] = -ab_y[k] * sd_3[k]
                 + sf_6[k];

        t_10[k] = -ab_y[k] * sd_4[k]
                  + sf_7[k];

        t_11[k] = -ab_y[k] * sd_5[k]
                  + sf_8[k];

        t_12[k] = -ab_z[k] * sd_0[k]
                  + sf_2[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_z, sd_1, sd_2, sd_3, sd_4, sd_5, \
                         sf_4, sf_5, sf_7, sf_8, sf_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_13[k] = -ab_z[k] * sd_1[k]
                  + sf_4[k];

        t_14[k] = -ab_z[k] * sd_2[k]
                  + sf_5[k];

        t_15[k] = -ab_z[k] * sd_3[k]
                  + sf_7[k];

        t_16[k] = -ab_z[k] * sd_4[k]
                  + sf_8[k];

        t_17[k] = -ab_z[k] * sd_5[k]
                  + sf_9[k];
    }
}

}  // namespace simdtrf
