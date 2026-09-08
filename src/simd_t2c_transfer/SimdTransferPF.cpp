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


#include "SimdTransferPF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_pf_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t sf, const size_t sg,
                   const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.75 * std::sqrt(10.0);
    const auto f_1 = 0.25 * std::sqrt(10.0);
    const auto f_2 = std::sqrt(15.0);
    const auto f_3 = 0.25 * std::sqrt(6.0);
    const auto f_4 = std::sqrt(6.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

#pragma omp simd aligned(ab_y, sf_1, sf_4, sf_6, sf_8, sg_3, sg_7, sg_10, \
                         sg_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -f_0 * ab_y[k] * sf_1[k]
                 + f_1 * ab_y[k] * sf_6[k]
                 + f_0 * sg_3[k]
                 - f_1 * sg_10[k];

        g_1[k] = -f_2 * ab_y[k] * sf_4[k]
                 + f_2 * sg_7[k];

        g_2[k] = f_3 * ab_y[k] * sf_1[k]
                 + f_3 * ab_y[k] * sf_6[k]
                 - f_4 * ab_y[k] * sf_8[k]
                 - f_3 * sg_3[k]
                 - f_3 * sg_10[k]
                 + f_4 * sg_12[k];
    }

#pragma omp simd aligned(ab_y, sf_0, sf_2, sf_3, sf_5, sf_7, sf_9, sg_1, sg_4, sg_6, sg_8, \
                         sg_11, sg_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = 1.5 * ab_y[k] * sf_2[k]
                 + 1.5 * ab_y[k] * sf_7[k]
                 - ab_y[k] * sf_9[k]
                 - 1.5 * sg_4[k]
                 - 1.5 * sg_11[k]
                 + sg_13[k];

        g_4[k] = f_3 * ab_y[k] * sf_0[k]
                 + f_3 * ab_y[k] * sf_3[k]
                 - f_4 * ab_y[k] * sf_5[k]
                 - f_3 * sg_1[k]
                 - f_3 * sg_6[k]
                 + f_4 * sg_8[k];

        g_5[k] = -f_5 * ab_y[k] * sf_2[k]
                 + f_5 * ab_y[k] * sf_7[k]
                 + f_5 * sg_4[k]
                 - f_5 * sg_11[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sf_0, sf_1, sf_3, sf_4, sf_6, sg_1, sg_4, sg_6, sg_8, \
                         sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_1 * ab_y[k] * sf_0[k]
                 + f_0 * ab_y[k] * sf_3[k]
                 + f_1 * sg_1[k]
                 - f_0 * sg_6[k];

        g_7[k] = -f_0 * ab_z[k] * sf_1[k]
                 + f_1 * ab_z[k] * sf_6[k]
                 + f_0 * sg_4[k]
                 - f_1 * sg_11[k];

        g_8[k] = -f_2 * ab_z[k] * sf_4[k]
                 + f_2 * sg_8[k];
    }

#pragma omp simd aligned(ab_z, sf_1, sf_2, sf_6, sf_7, sf_8, sf_9, sg_4, sg_5, sg_11, sg_12, \
                         sg_13, sg_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_3 * ab_z[k] * sf_1[k]
                 + f_3 * ab_z[k] * sf_6[k]
                 - f_4 * ab_z[k] * sf_8[k]
                 - f_3 * sg_4[k]
                 - f_3 * sg_11[k]
                 + f_4 * sg_13[k];

        g_10[k] = 1.5 * ab_z[k] * sf_2[k]
                  + 1.5 * ab_z[k] * sf_7[k]
                  - ab_z[k] * sf_9[k]
                  - 1.5 * sg_5[k]
                  - 1.5 * sg_12[k]
                  + sg_14[k];
    }

#pragma omp simd aligned(ab_z, sf_0, sf_2, sf_3, sf_5, sf_7, sg_2, sg_5, sg_7, sg_9, \
                         sg_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_3 * ab_z[k] * sf_0[k]
                  + f_3 * ab_z[k] * sf_3[k]
                  - f_4 * ab_z[k] * sf_5[k]
                  - f_3 * sg_2[k]
                  - f_3 * sg_7[k]
                  + f_4 * sg_9[k];

        g_12[k] = -f_5 * ab_z[k] * sf_2[k]
                  + f_5 * ab_z[k] * sf_7[k]
                  + f_5 * sg_5[k]
                  - f_5 * sg_12[k];

        g_13[k] = -f_1 * ab_z[k] * sf_0[k]
                  + f_0 * ab_z[k] * sf_3[k]
                  + f_1 * sg_2[k]
                  - f_0 * sg_7[k];
    }

#pragma omp simd aligned(ab_x, sf_1, sf_4, sf_6, sf_8, sg_1, sg_4, sg_6, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_0 * ab_x[k] * sf_1[k]
                  + f_1 * ab_x[k] * sf_6[k]
                  + f_0 * sg_1[k]
                  - f_1 * sg_6[k];

        g_15[k] = -f_2 * ab_x[k] * sf_4[k]
                  + f_2 * sg_4[k];

        g_16[k] = f_3 * ab_x[k] * sf_1[k]
                  + f_3 * ab_x[k] * sf_6[k]
                  - f_4 * ab_x[k] * sf_8[k]
                  - f_3 * sg_1[k]
                  - f_3 * sg_6[k]
                  + f_4 * sg_8[k];
    }

#pragma omp simd aligned(ab_x, sf_0, sf_2, sf_3, sf_5, sf_7, sf_9, sg_0, sg_2, sg_3, sg_5, \
                         sg_7, sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = 1.5 * ab_x[k] * sf_2[k]
                  + 1.5 * ab_x[k] * sf_7[k]
                  - ab_x[k] * sf_9[k]
                  - 1.5 * sg_2[k]
                  - 1.5 * sg_7[k]
                  + sg_9[k];

        g_18[k] = f_3 * ab_x[k] * sf_0[k]
                  + f_3 * ab_x[k] * sf_3[k]
                  - f_4 * ab_x[k] * sf_5[k]
                  - f_3 * sg_0[k]
                  - f_3 * sg_3[k]
                  + f_4 * sg_5[k];

        g_19[k] = -f_5 * ab_x[k] * sf_2[k]
                  + f_5 * ab_x[k] * sf_7[k]
                  + f_5 * sg_2[k]
                  - f_5 * sg_7[k];
    }

#pragma omp simd aligned(ab_x, sf_0, sf_3, sg_0, sg_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_1 * ab_x[k] * sf_0[k]
                  + f_0 * ab_x[k] * sf_3[k]
                  + f_1 * sg_0[k]
                  - f_0 * sg_3[k];
    }
}

auto
compute_hrr_pf(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t sf, const size_t sg, const size_t nmax) -> void
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
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, sf_0, sf_1, sf_2, sf_3, sf_4, sg_0, \
                         sg_1, sg_2, sg_3, sg_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * sf_0[k]
                 + sg_0[k];

        t_1[k] = -ab_x[k] * sf_1[k]
                 + sg_1[k];

        t_2[k] = -ab_x[k] * sf_2[k]
                 + sg_2[k];

        t_3[k] = -ab_x[k] * sf_3[k]
                 + sg_3[k];

        t_4[k] = -ab_x[k] * sf_4[k]
                 + sg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, sf_5, sf_6, sf_7, sf_8, sf_9, sg_5, \
                         sg_6, sg_7, sg_8, sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * sf_5[k]
                 + sg_5[k];

        t_6[k] = -ab_x[k] * sf_6[k]
                 + sg_6[k];

        t_7[k] = -ab_x[k] * sf_7[k]
                 + sg_7[k];

        t_8[k] = -ab_x[k] * sf_8[k]
                 + sg_8[k];

        t_9[k] = -ab_x[k] * sf_9[k]
                 + sg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, sf_0, sf_1, sf_2, sf_3, sf_4, \
                         sg_1, sg_3, sg_4, sg_6, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_y[k] * sf_0[k]
                  + sg_1[k];

        t_11[k] = -ab_y[k] * sf_1[k]
                  + sg_3[k];

        t_12[k] = -ab_y[k] * sf_2[k]
                  + sg_4[k];

        t_13[k] = -ab_y[k] * sf_3[k]
                  + sg_6[k];

        t_14[k] = -ab_y[k] * sf_4[k]
                  + sg_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_y, sf_5, sf_6, sf_7, sf_8, sf_9, \
                         sg_8, sg_10, sg_11, sg_12, sg_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_y[k] * sf_5[k]
                  + sg_8[k];

        t_16[k] = -ab_y[k] * sf_6[k]
                  + sg_10[k];

        t_17[k] = -ab_y[k] * sf_7[k]
                  + sg_11[k];

        t_18[k] = -ab_y[k] * sf_8[k]
                  + sg_12[k];

        t_19[k] = -ab_y[k] * sf_9[k]
                  + sg_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_z, sf_0, sf_1, sf_2, sf_3, sf_4, \
                         sg_2, sg_4, sg_5, sg_7, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_z[k] * sf_0[k]
                  + sg_2[k];

        t_21[k] = -ab_z[k] * sf_1[k]
                  + sg_4[k];

        t_22[k] = -ab_z[k] * sf_2[k]
                  + sg_5[k];

        t_23[k] = -ab_z[k] * sf_3[k]
                  + sg_7[k];

        t_24[k] = -ab_z[k] * sf_4[k]
                  + sg_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_z, sf_5, sf_6, sf_7, sf_8, sf_9, \
                         sg_9, sg_11, sg_12, sg_13, sg_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_z[k] * sf_5[k]
                  + sg_9[k];

        t_26[k] = -ab_z[k] * sf_6[k]
                  + sg_11[k];

        t_27[k] = -ab_z[k] * sf_7[k]
                  + sg_12[k];

        t_28[k] = -ab_z[k] * sf_8[k]
                  + sg_13[k];

        t_29[k] = -ab_z[k] * sf_9[k]
                  + sg_14[k];
    }
}

}  // namespace simdovl
