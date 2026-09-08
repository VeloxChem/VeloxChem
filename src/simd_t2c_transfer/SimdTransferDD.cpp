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


#include "SimdTransferDD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_dd(double *values, const size_t nvalues, CSimdMatrix &buffer,
               const CSimdMatrix &coordinates, const size_t sd, const size_t sf, const size_t sg,
               const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 * std::sqrt(3.0);
    const auto f_1 = std::sqrt(3.0);
    const auto f_2 = 2.0 * std::sqrt(3.0);
    const auto f_3 = 0.25 * std::sqrt(3.0);

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
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;

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

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

#pragma omp simd aligned(ab_x, ab_y, sd_1, sd_4, sf_1, sf_3, sf_4, sf_7, sg_3, \
                         sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 3.0 * ab_x[k] * ab_y[k] * sd_1[k]
                 - 3.0 * ab_y[k] * sf_1[k]
                 - 3.0 * ab_x[k] * sf_3[k]
                 + 3.0 * sg_3[k];

        g_1[k] = 3.0 * ab_x[k] * ab_y[k] * sd_4[k]
                 - 3.0 * ab_y[k] * sf_4[k]
                 - 3.0 * ab_x[k] * sf_7[k]
                 + 3.0 * sg_7[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_0, sd_3, sd_5, sf_0, sf_1, sf_3, sf_5, sf_6, sf_8, \
                         sg_1, sg_6, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_0 * ab_x[k] * ab_y[k] * sd_0[k]
                 - f_0 * ab_x[k] * ab_y[k] * sd_3[k]
                 + f_1 * ab_x[k] * ab_y[k] * sd_5[k]
                 + f_0 * ab_y[k] * sf_0[k]
                 + f_0 * ab_x[k] * sf_1[k]
                 + f_0 * ab_y[k] * sf_3[k]
                 - f_1 * ab_y[k] * sf_5[k]
                 + f_0 * ab_x[k] * sf_6[k]
                 - f_1 * ab_x[k] * sf_8[k]
                 - f_0 * sg_1[k]
                 - f_0 * sg_6[k]
                 + f_1 * sg_8[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_0, sd_2, sd_3, sf_0, sf_1, sf_2, sf_3, sf_4, sf_6, \
                         sg_1, sg_4, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = 3.0 * ab_x[k] * ab_y[k] * sd_2[k]
                 - 3.0 * ab_y[k] * sf_2[k]
                 - 3.0 * ab_x[k] * sf_4[k]
                 + 3.0 * sg_4[k];

        g_4[k] = 1.5 * ab_x[k] * ab_y[k] * sd_0[k]
                 - 1.5 * ab_x[k] * ab_y[k] * sd_3[k]
                 - 1.5 * ab_y[k] * sf_0[k]
                 - 1.5 * ab_x[k] * sf_1[k]
                 + 1.5 * ab_y[k] * sf_3[k]
                 + 1.5 * ab_x[k] * sf_6[k]
                 + 1.5 * sg_1[k]
                 - 1.5 * sg_6[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sd_1, sd_4, sf_3, sf_4, sf_7, sf_8, sg_7, \
                         sg_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = 3.0 * ab_y[k] * ab_z[k] * sd_1[k]
                 - 3.0 * ab_z[k] * sf_3[k]
                 - 3.0 * ab_y[k] * sf_4[k]
                 + 3.0 * sg_7[k];

        g_6[k] = 3.0 * ab_y[k] * ab_z[k] * sd_4[k]
                 - 3.0 * ab_z[k] * sf_7[k]
                 - 3.0 * ab_y[k] * sf_8[k]
                 + 3.0 * sg_12[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sd_0, sd_3, sd_5, sf_1, sf_2, sf_6, sf_7, sf_8, sf_9, \
                         sg_4, sg_11, sg_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_0 * ab_y[k] * ab_z[k] * sd_0[k]
                 - f_0 * ab_y[k] * ab_z[k] * sd_3[k]
                 + f_1 * ab_y[k] * ab_z[k] * sd_5[k]
                 + f_0 * ab_z[k] * sf_1[k]
                 + f_0 * ab_y[k] * sf_2[k]
                 + f_0 * ab_z[k] * sf_6[k]
                 + f_0 * ab_y[k] * sf_7[k]
                 - f_1 * ab_z[k] * sf_8[k]
                 - f_1 * ab_y[k] * sf_9[k]
                 - f_0 * sg_4[k]
                 - f_0 * sg_11[k]
                 + f_1 * sg_13[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sd_0, sd_2, sd_3, sf_1, sf_2, sf_4, sf_5, sf_6, sf_7, \
                         sg_4, sg_8, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = 3.0 * ab_y[k] * ab_z[k] * sd_2[k]
                 - 3.0 * ab_z[k] * sf_4[k]
                 - 3.0 * ab_y[k] * sf_5[k]
                 + 3.0 * sg_8[k];

        g_9[k] = 1.5 * ab_y[k] * ab_z[k] * sd_0[k]
                 - 1.5 * ab_y[k] * ab_z[k] * sd_3[k]
                 - 1.5 * ab_z[k] * sf_1[k]
                 - 1.5 * ab_y[k] * sf_2[k]
                 + 1.5 * ab_z[k] * sf_6[k]
                 + 1.5 * ab_y[k] * sf_7[k]
                 + 1.5 * sg_4[k]
                 - 1.5 * sg_11[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, sd_1, sf_1, sf_3, sf_4, sg_1, sg_6, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_0 * ab_x[k] * ab_x[k] * sd_1[k]
                  - f_0 * ab_y[k] * ab_y[k] * sd_1[k]
                  + f_1 * ab_z[k] * ab_z[k] * sd_1[k]
                  + f_1 * ab_x[k] * sf_1[k]
                  + f_1 * ab_y[k] * sf_3[k]
                  - f_2 * ab_z[k] * sf_4[k]
                  - f_0 * sg_1[k]
                  - f_0 * sg_6[k]
                  + f_1 * sg_8[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, sd_4, sf_4, sf_7, sf_8, sg_4, sg_11, \
                         sg_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_0 * ab_x[k] * ab_x[k] * sd_4[k]
                  - f_0 * ab_y[k] * ab_y[k] * sd_4[k]
                  + f_1 * ab_z[k] * ab_z[k] * sd_4[k]
                  + f_1 * ab_x[k] * sf_4[k]
                  + f_1 * ab_y[k] * sf_7[k]
                  - f_2 * ab_z[k] * sf_8[k]
                  - f_0 * sg_4[k]
                  - f_0 * sg_11[k]
                  + f_1 * sg_13[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, sd_0, sd_3, sd_5, sf_0, sf_1, sf_2, sf_3, sf_5, \
                         sf_6, sf_7, sf_8, sf_9, sg_0, sg_3, sg_5, sg_10, sg_12, \
                         sg_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = 0.25 * ab_x[k] * ab_x[k] * sd_0[k]
                  + 0.25 * ab_y[k] * ab_y[k] * sd_0[k]
                  - 0.5 * ab_z[k] * ab_z[k] * sd_0[k]
                  + 0.25 * ab_x[k] * ab_x[k] * sd_3[k]
                  + 0.25 * ab_y[k] * ab_y[k] * sd_3[k]
                  - 0.5 * ab_z[k] * ab_z[k] * sd_3[k]
                  - 0.5 * ab_x[k] * ab_x[k] * sd_5[k]
                  - 0.5 * ab_y[k] * ab_y[k] * sd_5[k]
                  + ab_z[k] * ab_z[k] * sd_5[k]
                  - 0.5 * ab_x[k] * sf_0[k]
                  - 0.5 * ab_y[k] * sf_1[k]
                  + ab_z[k] * sf_2[k]
                  - 0.5 * ab_x[k] * sf_3[k]
                  + ab_x[k] * sf_5[k]
                  - 0.5 * ab_y[k] * sf_6[k]
                  + ab_z[k] * sf_7[k]
                  + ab_y[k] * sf_8[k]
                  - 2.0 * ab_z[k] * sf_9[k]
                  + 0.25 * sg_0[k]
                  + 0.5 * sg_3[k]
                  - sg_5[k]
                  + 0.25 * sg_10[k]
                  - sg_12[k]
                  + sg_14[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, sd_2, sf_2, sf_4, sf_5, sg_2, sg_7, \
                         sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_0 * ab_x[k] * ab_x[k] * sd_2[k]
                  - f_0 * ab_y[k] * ab_y[k] * sd_2[k]
                  + f_1 * ab_z[k] * ab_z[k] * sd_2[k]
                  + f_1 * ab_x[k] * sf_2[k]
                  + f_1 * ab_y[k] * sf_4[k]
                  - f_2 * ab_z[k] * sf_5[k]
                  - f_0 * sg_2[k]
                  - f_0 * sg_7[k]
                  + f_1 * sg_9[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, sd_0, sd_3, sf_0, sf_1, sf_2, sf_3, sf_6, sf_7, \
                         sg_0, sg_5, sg_10, sg_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_3 * ab_x[k] * ab_x[k] * sd_0[k]
                  - f_3 * ab_y[k] * ab_y[k] * sd_0[k]
                  + f_0 * ab_z[k] * ab_z[k] * sd_0[k]
                  + f_3 * ab_x[k] * ab_x[k] * sd_3[k]
                  + f_3 * ab_y[k] * ab_y[k] * sd_3[k]
                  - f_0 * ab_z[k] * ab_z[k] * sd_3[k]
                  + f_0 * ab_x[k] * sf_0[k]
                  + f_0 * ab_y[k] * sf_1[k]
                  - f_1 * ab_z[k] * sf_2[k]
                  - f_0 * ab_x[k] * sf_3[k]
                  - f_0 * ab_y[k] * sf_6[k]
                  + f_1 * ab_z[k] * sf_7[k]
                  - f_3 * sg_0[k]
                  + f_0 * sg_5[k]
                  + f_3 * sg_10[k]
                  - f_0 * sg_12[k];
    }

#pragma omp simd aligned(ab_x, ab_z, sd_1, sd_4, sf_1, sf_4, sf_8, sg_4, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = 3.0 * ab_x[k] * ab_z[k] * sd_1[k]
                  - 3.0 * ab_z[k] * sf_1[k]
                  - 3.0 * ab_x[k] * sf_4[k]
                  + 3.0 * sg_4[k];

        g_16[k] = 3.0 * ab_x[k] * ab_z[k] * sd_4[k]
                  - 3.0 * ab_z[k] * sf_4[k]
                  - 3.0 * ab_x[k] * sf_8[k]
                  + 3.0 * sg_8[k];
    }

#pragma omp simd aligned(ab_x, ab_z, sd_0, sd_3, sd_5, sf_0, sf_2, sf_3, sf_5, sf_7, sf_9, \
                         sg_2, sg_7, sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_0 * ab_x[k] * ab_z[k] * sd_0[k]
                  - f_0 * ab_x[k] * ab_z[k] * sd_3[k]
                  + f_1 * ab_x[k] * ab_z[k] * sd_5[k]
                  + f_0 * ab_z[k] * sf_0[k]
                  + f_0 * ab_x[k] * sf_2[k]
                  + f_0 * ab_z[k] * sf_3[k]
                  - f_1 * ab_z[k] * sf_5[k]
                  + f_0 * ab_x[k] * sf_7[k]
                  - f_1 * ab_x[k] * sf_9[k]
                  - f_0 * sg_2[k]
                  - f_0 * sg_7[k]
                  + f_1 * sg_9[k];
    }

#pragma omp simd aligned(ab_x, ab_z, sd_0, sd_2, sd_3, sf_0, sf_2, sf_3, sf_5, sf_7, sg_2, \
                         sg_5, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = 3.0 * ab_x[k] * ab_z[k] * sd_2[k]
                  - 3.0 * ab_z[k] * sf_2[k]
                  - 3.0 * ab_x[k] * sf_5[k]
                  + 3.0 * sg_5[k];

        g_19[k] = 1.5 * ab_x[k] * ab_z[k] * sd_0[k]
                  - 1.5 * ab_x[k] * ab_z[k] * sd_3[k]
                  - 1.5 * ab_z[k] * sf_0[k]
                  - 1.5 * ab_x[k] * sf_2[k]
                  + 1.5 * ab_z[k] * sf_3[k]
                  + 1.5 * ab_x[k] * sf_7[k]
                  + 1.5 * sg_2[k]
                  - 1.5 * sg_7[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_1, sd_4, sf_1, sf_3, sf_4, sf_7, sg_1, sg_4, sg_6, \
                         sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = 1.5 * ab_x[k] * ab_x[k] * sd_1[k]
                  - 1.5 * ab_y[k] * ab_y[k] * sd_1[k]
                  - 3.0 * ab_x[k] * sf_1[k]
                  + 3.0 * ab_y[k] * sf_3[k]
                  + 1.5 * sg_1[k]
                  - 1.5 * sg_6[k];

        g_21[k] = 1.5 * ab_x[k] * ab_x[k] * sd_4[k]
                  - 1.5 * ab_y[k] * ab_y[k] * sd_4[k]
                  - 3.0 * ab_x[k] * sf_4[k]
                  + 3.0 * ab_y[k] * sf_7[k]
                  + 1.5 * sg_4[k]
                  - 1.5 * sg_11[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_0, sd_3, sd_5, sf_0, sf_1, sf_3, sf_5, sf_6, sf_8, \
                         sg_0, sg_5, sg_10, sg_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_3 * ab_x[k] * ab_x[k] * sd_0[k]
                  + f_3 * ab_y[k] * ab_y[k] * sd_0[k]
                  - f_3 * ab_x[k] * ab_x[k] * sd_3[k]
                  + f_3 * ab_y[k] * ab_y[k] * sd_3[k]
                  + f_0 * ab_x[k] * ab_x[k] * sd_5[k]
                  - f_0 * ab_y[k] * ab_y[k] * sd_5[k]
                  + f_0 * ab_x[k] * sf_0[k]
                  - f_0 * ab_y[k] * sf_1[k]
                  + f_0 * ab_x[k] * sf_3[k]
                  - f_1 * ab_x[k] * sf_5[k]
                  - f_0 * ab_y[k] * sf_6[k]
                  + f_1 * ab_y[k] * sf_8[k]
                  - f_3 * sg_0[k]
                  + f_0 * sg_5[k]
                  + f_3 * sg_10[k]
                  - f_0 * sg_12[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_2, sf_2, sf_4, sg_2, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = 1.5 * ab_x[k] * ab_x[k] * sd_2[k]
                  - 1.5 * ab_y[k] * ab_y[k] * sd_2[k]
                  - 3.0 * ab_x[k] * sf_2[k]
                  + 3.0 * ab_y[k] * sf_4[k]
                  + 1.5 * sg_2[k]
                  - 1.5 * sg_7[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_0, sd_3, sf_0, sf_1, sf_3, sf_6, sg_0, sg_3, \
                         sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 0.75 * ab_x[k] * ab_x[k] * sd_0[k]
                  - 0.75 * ab_y[k] * ab_y[k] * sd_0[k]
                  - 0.75 * ab_x[k] * ab_x[k] * sd_3[k]
                  + 0.75 * ab_y[k] * ab_y[k] * sd_3[k]
                  - 1.5 * ab_x[k] * sf_0[k]
                  + 1.5 * ab_y[k] * sf_1[k]
                  + 1.5 * ab_x[k] * sf_3[k]
                  - 1.5 * ab_y[k] * sf_6[k]
                  + 0.75 * sg_0[k]
                  - 1.5 * sg_3[k]
                  + 0.75 * sg_10[k];
    }
}

auto
compute_hrr_dd_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t sd, const size_t sf,
                   const size_t sg, const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 * std::sqrt(3.0);
    const auto f_1 = std::sqrt(3.0);
    const auto f_2 = 2.0 * std::sqrt(3.0);
    const auto f_3 = 0.25 * std::sqrt(3.0);

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
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;

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

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

#pragma omp simd aligned(ab_x, ab_y, sd_1, sd_4, sf_1, sf_3, sf_4, sf_7, sg_3, \
                         sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 3.0 * ab_x[k] * ab_y[k] * sd_1[k]
                 - 3.0 * ab_y[k] * sf_1[k]
                 - 3.0 * ab_x[k] * sf_3[k]
                 + 3.0 * sg_3[k];

        g_1[k] = 3.0 * ab_x[k] * ab_y[k] * sd_4[k]
                 - 3.0 * ab_y[k] * sf_4[k]
                 - 3.0 * ab_x[k] * sf_7[k]
                 + 3.0 * sg_7[k];
        g_5[k] = g_1[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_0, sd_3, sd_5, sf_0, sf_1, sf_3, sf_5, sf_6, sf_8, \
                         sg_1, sg_6, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_0 * ab_x[k] * ab_y[k] * sd_0[k]
                 - f_0 * ab_x[k] * ab_y[k] * sd_3[k]
                 + f_1 * ab_x[k] * ab_y[k] * sd_5[k]
                 + f_0 * ab_y[k] * sf_0[k]
                 + f_0 * ab_x[k] * sf_1[k]
                 + f_0 * ab_y[k] * sf_3[k]
                 - f_1 * ab_y[k] * sf_5[k]
                 + f_0 * ab_x[k] * sf_6[k]
                 - f_1 * ab_x[k] * sf_8[k]
                 - f_0 * sg_1[k]
                 - f_0 * sg_6[k]
                 + f_1 * sg_8[k];
        g_10[k] = g_2[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_0, sd_2, sd_3, sf_0, sf_1, sf_2, sf_3, sf_4, sf_6, \
                         sg_1, sg_4, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = 3.0 * ab_x[k] * ab_y[k] * sd_2[k]
                 - 3.0 * ab_y[k] * sf_2[k]
                 - 3.0 * ab_x[k] * sf_4[k]
                 + 3.0 * sg_4[k];
        g_15[k] = g_3[k];

        g_4[k] = 1.5 * ab_x[k] * ab_y[k] * sd_0[k]
                 - 1.5 * ab_x[k] * ab_y[k] * sd_3[k]
                 - 1.5 * ab_y[k] * sf_0[k]
                 - 1.5 * ab_x[k] * sf_1[k]
                 + 1.5 * ab_y[k] * sf_3[k]
                 + 1.5 * ab_x[k] * sf_6[k]
                 + 1.5 * sg_1[k]
                 - 1.5 * sg_6[k];
        g_20[k] = g_4[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sd_4, sf_7, sf_8, sg_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = 3.0 * ab_y[k] * ab_z[k] * sd_4[k]
                 - 3.0 * ab_z[k] * sf_7[k]
                 - 3.0 * ab_y[k] * sf_8[k]
                 + 3.0 * sg_12[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sd_0, sd_3, sd_5, sf_1, sf_2, sf_6, sf_7, sf_8, sf_9, \
                         sg_4, sg_11, sg_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_0 * ab_y[k] * ab_z[k] * sd_0[k]
                 - f_0 * ab_y[k] * ab_z[k] * sd_3[k]
                 + f_1 * ab_y[k] * ab_z[k] * sd_5[k]
                 + f_0 * ab_z[k] * sf_1[k]
                 + f_0 * ab_y[k] * sf_2[k]
                 + f_0 * ab_z[k] * sf_6[k]
                 + f_0 * ab_y[k] * sf_7[k]
                 - f_1 * ab_z[k] * sf_8[k]
                 - f_1 * ab_y[k] * sf_9[k]
                 - f_0 * sg_4[k]
                 - f_0 * sg_11[k]
                 + f_1 * sg_13[k];
        g_11[k] = g_7[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sd_0, sd_2, sd_3, sf_1, sf_2, sf_4, sf_5, sf_6, sf_7, \
                         sg_4, sg_8, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = 3.0 * ab_y[k] * ab_z[k] * sd_2[k]
                 - 3.0 * ab_z[k] * sf_4[k]
                 - 3.0 * ab_y[k] * sf_5[k]
                 + 3.0 * sg_8[k];
        g_16[k] = g_8[k];

        g_9[k] = 1.5 * ab_y[k] * ab_z[k] * sd_0[k]
                 - 1.5 * ab_y[k] * ab_z[k] * sd_3[k]
                 - 1.5 * ab_z[k] * sf_1[k]
                 - 1.5 * ab_y[k] * sf_2[k]
                 + 1.5 * ab_z[k] * sf_6[k]
                 + 1.5 * ab_y[k] * sf_7[k]
                 + 1.5 * sg_4[k]
                 - 1.5 * sg_11[k];
        g_21[k] = g_9[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, sd_0, sd_3, sd_5, sf_0, sf_1, sf_2, sf_3, sf_5, \
                         sf_6, sf_7, sf_8, sf_9, sg_0, sg_3, sg_5, sg_10, sg_12, \
                         sg_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = 0.25 * ab_x[k] * ab_x[k] * sd_0[k]
                  + 0.25 * ab_y[k] * ab_y[k] * sd_0[k]
                  - 0.5 * ab_z[k] * ab_z[k] * sd_0[k]
                  + 0.25 * ab_x[k] * ab_x[k] * sd_3[k]
                  + 0.25 * ab_y[k] * ab_y[k] * sd_3[k]
                  - 0.5 * ab_z[k] * ab_z[k] * sd_3[k]
                  - 0.5 * ab_x[k] * ab_x[k] * sd_5[k]
                  - 0.5 * ab_y[k] * ab_y[k] * sd_5[k]
                  + ab_z[k] * ab_z[k] * sd_5[k]
                  - 0.5 * ab_x[k] * sf_0[k]
                  - 0.5 * ab_y[k] * sf_1[k]
                  + ab_z[k] * sf_2[k]
                  - 0.5 * ab_x[k] * sf_3[k]
                  + ab_x[k] * sf_5[k]
                  - 0.5 * ab_y[k] * sf_6[k]
                  + ab_z[k] * sf_7[k]
                  + ab_y[k] * sf_8[k]
                  - 2.0 * ab_z[k] * sf_9[k]
                  + 0.25 * sg_0[k]
                  + 0.5 * sg_3[k]
                  - sg_5[k]
                  + 0.25 * sg_10[k]
                  - sg_12[k]
                  + sg_14[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, sd_2, sf_2, sf_4, sf_5, sg_2, sg_7, \
                         sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_0 * ab_x[k] * ab_x[k] * sd_2[k]
                  - f_0 * ab_y[k] * ab_y[k] * sd_2[k]
                  + f_1 * ab_z[k] * ab_z[k] * sd_2[k]
                  + f_1 * ab_x[k] * sf_2[k]
                  + f_1 * ab_y[k] * sf_4[k]
                  - f_2 * ab_z[k] * sf_5[k]
                  - f_0 * sg_2[k]
                  - f_0 * sg_7[k]
                  + f_1 * sg_9[k];
        g_17[k] = g_13[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, sd_0, sd_3, sf_0, sf_1, sf_2, sf_3, sf_6, sf_7, \
                         sg_0, sg_5, sg_10, sg_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_3 * ab_x[k] * ab_x[k] * sd_0[k]
                  - f_3 * ab_y[k] * ab_y[k] * sd_0[k]
                  + f_0 * ab_z[k] * ab_z[k] * sd_0[k]
                  + f_3 * ab_x[k] * ab_x[k] * sd_3[k]
                  + f_3 * ab_y[k] * ab_y[k] * sd_3[k]
                  - f_0 * ab_z[k] * ab_z[k] * sd_3[k]
                  + f_0 * ab_x[k] * sf_0[k]
                  + f_0 * ab_y[k] * sf_1[k]
                  - f_1 * ab_z[k] * sf_2[k]
                  - f_0 * ab_x[k] * sf_3[k]
                  - f_0 * ab_y[k] * sf_6[k]
                  + f_1 * ab_z[k] * sf_7[k]
                  - f_3 * sg_0[k]
                  + f_0 * sg_5[k]
                  + f_3 * sg_10[k]
                  - f_0 * sg_12[k];
        g_22[k] = g_14[k];
    }

#pragma omp simd aligned(ab_x, ab_z, sd_0, sd_2, sd_3, sf_0, sf_2, sf_3, sf_5, sf_7, sg_2, \
                         sg_5, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = 3.0 * ab_x[k] * ab_z[k] * sd_2[k]
                  - 3.0 * ab_z[k] * sf_2[k]
                  - 3.0 * ab_x[k] * sf_5[k]
                  + 3.0 * sg_5[k];

        g_19[k] = 1.5 * ab_x[k] * ab_z[k] * sd_0[k]
                  - 1.5 * ab_x[k] * ab_z[k] * sd_3[k]
                  - 1.5 * ab_z[k] * sf_0[k]
                  - 1.5 * ab_x[k] * sf_2[k]
                  + 1.5 * ab_z[k] * sf_3[k]
                  + 1.5 * ab_x[k] * sf_7[k]
                  + 1.5 * sg_2[k]
                  - 1.5 * sg_7[k];
        g_23[k] = g_19[k];
    }

#pragma omp simd aligned(ab_x, ab_y, sd_0, sd_3, sf_0, sf_1, sf_3, sf_6, sg_0, sg_3, \
                         sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 0.75 * ab_x[k] * ab_x[k] * sd_0[k]
                  - 0.75 * ab_y[k] * ab_y[k] * sd_0[k]
                  - 0.75 * ab_x[k] * ab_x[k] * sd_3[k]
                  + 0.75 * ab_y[k] * ab_y[k] * sd_3[k]
                  - 1.5 * ab_x[k] * sf_0[k]
                  + 1.5 * ab_y[k] * sf_1[k]
                  + 1.5 * ab_x[k] * sf_3[k]
                  - 1.5 * ab_y[k] * sf_6[k]
                  + 0.75 * sg_0[k]
                  - 1.5 * sg_3[k]
                  + 0.75 * sg_10[k];
    }
}

}  // namespace simdovl
