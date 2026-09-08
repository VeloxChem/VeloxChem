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


#include "SimdTransferFD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_fd_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t fp, const size_t gp,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.75 * std::sqrt(30.0);
    const auto f_1 = 0.25 * std::sqrt(30.0);
    const auto f_2 = 0.375 * std::sqrt(10.0);
    const auto f_3 = 0.75 * std::sqrt(10.0);
    const auto f_4 = 0.125 * std::sqrt(10.0);
    const auto f_5 = 0.25 * std::sqrt(10.0);
    const auto f_6 = 0.375 * std::sqrt(30.0);
    const auto f_7 = 0.125 * std::sqrt(30.0);
    const auto f_8 = 3.0 * std::sqrt(5.0);
    const auto f_9 = 0.5 * std::sqrt(15.0);
    const auto f_10 = std::sqrt(15.0);
    const auto f_11 = 1.5 * std::sqrt(5.0);
    const auto f_12 = 0.75 * std::sqrt(2.0);
    const auto f_13 = 3.0 * std::sqrt(2.0);
    const auto f_14 = 0.125 * std::sqrt(6.0);
    const auto f_15 = 0.25 * std::sqrt(6.0);
    const auto f_16 = 0.5 * std::sqrt(6.0);
    const auto f_17 = std::sqrt(6.0);
    const auto f_18 = 0.375 * std::sqrt(2.0);
    const auto f_19 = 1.5 * std::sqrt(2.0);
    const auto f_20 = 1.5 * std::sqrt(3.0);
    const auto f_21 = std::sqrt(3.0);
    const auto f_22 = 0.75 * std::sqrt(3.0);
    const auto f_23 = 0.5 * std::sqrt(3.0);
    const auto f_24 = 0.25 * std::sqrt(15.0);
    const auto f_25 = 0.75 * std::sqrt(5.0);

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
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);
    const auto *fp_21 = buffer.data(fp + 21);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_24 = buffer.data(fp + 24);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_44 = buffer.data(gp + 44);

#pragma omp simd aligned(ab_x, ab_y, fp_4, fp_5, fp_19, fp_20, gp_4, gp_11, gp_19, \
                         gp_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_x[k] * fp_4[k]
                 - f_1 * ab_x[k] * fp_19[k]
                 + f_0 * gp_4[k]
                 - f_1 * gp_19[k];

        g_1[k] = f_0 * ab_y[k] * fp_5[k]
                 - f_1 * ab_y[k] * fp_20[k]
                 + f_0 * gp_11[k]
                 - f_1 * gp_32[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fp_3, fp_4, fp_5, fp_18, fp_19, fp_20, gp_3, gp_10, \
                         gp_14, gp_18, gp_31, gp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_2 * ab_x[k] * fp_3[k]
                 - f_2 * ab_y[k] * fp_4[k]
                 + f_3 * ab_z[k] * fp_5[k]
                 + f_4 * ab_x[k] * fp_18[k]
                 + f_4 * ab_y[k] * fp_19[k]
                 - f_5 * ab_z[k] * fp_20[k]
                 - f_2 * gp_3[k]
                 - f_2 * gp_10[k]
                 + f_3 * gp_14[k]
                 + f_4 * gp_18[k]
                 + f_4 * gp_31[k]
                 - f_5 * gp_35[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_3, fp_4, fp_5, fp_18, fp_19, fp_20, gp_3, gp_5, gp_10, \
                         gp_18, gp_20, gp_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_0 * ab_x[k] * fp_5[k]
                 - f_1 * ab_x[k] * fp_20[k]
                 + f_0 * gp_5[k]
                 - f_1 * gp_20[k];

        g_4[k] = f_6 * ab_x[k] * fp_3[k]
                 - f_6 * ab_y[k] * fp_4[k]
                 - f_7 * ab_x[k] * fp_18[k]
                 + f_7 * ab_y[k] * fp_19[k]
                 + f_6 * gp_3[k]
                 - f_6 * gp_10[k]
                 - f_7 * gp_18[k]
                 + f_7 * gp_31[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fp_12, fp_13, fp_14, gp_12, gp_13, gp_14, gp_22, \
                         gp_23, gp_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_8 * ab_x[k] * fp_13[k]
                 + f_8 * gp_13[k];

        g_6[k] = f_8 * ab_y[k] * fp_14[k]
                 + f_8 * gp_23[k];

        g_7[k] = -f_9 * ab_x[k] * fp_12[k]
                 - f_9 * ab_y[k] * fp_13[k]
                 + f_10 * ab_z[k] * fp_14[k]
                 - f_9 * gp_12[k]
                 - f_9 * gp_22[k]
                 + f_10 * gp_26[k];

        g_8[k] = f_8 * ab_x[k] * fp_14[k]
                 + f_8 * gp_14[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_4, fp_12, fp_13, fp_19, fp_25, gp_4, gp_12, gp_19, \
                         gp_22, gp_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_11 * ab_x[k] * fp_12[k]
                 - f_11 * ab_y[k] * fp_13[k]
                 + f_11 * gp_12[k]
                 - f_11 * gp_22[k];

        g_10[k] = -f_12 * ab_x[k] * fp_4[k]
                  - f_12 * ab_x[k] * fp_19[k]
                  + f_13 * ab_x[k] * fp_25[k]
                  - f_12 * gp_4[k]
                  - f_12 * gp_19[k]
                  + f_13 * gp_25[k];
    }

#pragma omp simd aligned(ab_y, fp_5, fp_20, fp_26, gp_11, gp_32, \
                         gp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_12 * ab_y[k] * fp_5[k]
                  - f_12 * ab_y[k] * fp_20[k]
                  + f_13 * ab_y[k] * fp_26[k]
                  - f_12 * gp_11[k]
                  - f_12 * gp_32[k]
                  + f_13 * gp_38[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fp_3, fp_4, fp_5, fp_18, fp_19, fp_20, fp_24, \
                         fp_25, fp_26, gp_3, gp_10, gp_14, gp_18, gp_24, gp_31, gp_35, gp_37, \
                         gp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_14 * ab_x[k] * fp_3[k]
                  + f_14 * ab_y[k] * fp_4[k]
                  - f_15 * ab_z[k] * fp_5[k]
                  + f_14 * ab_x[k] * fp_18[k]
                  + f_14 * ab_y[k] * fp_19[k]
                  - f_15 * ab_z[k] * fp_20[k]
                  - f_16 * ab_x[k] * fp_24[k]
                  - f_16 * ab_y[k] * fp_25[k]
                  + f_17 * ab_z[k] * fp_26[k]
                  + f_14 * gp_3[k]
                  + f_14 * gp_10[k]
                  - f_15 * gp_14[k]
                  + f_14 * gp_18[k]
                  - f_16 * gp_24[k]
                  + f_14 * gp_31[k]
                  - f_15 * gp_35[k]
                  - f_16 * gp_37[k]
                  + f_17 * gp_41[k];
    }

#pragma omp simd aligned(ab_x, fp_5, fp_20, fp_26, gp_5, gp_20, gp_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_12 * ab_x[k] * fp_5[k]
                  - f_12 * ab_x[k] * fp_20[k]
                  + f_13 * ab_x[k] * fp_26[k]
                  - f_12 * gp_5[k]
                  - f_12 * gp_20[k]
                  + f_13 * gp_26[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_3, fp_4, fp_18, fp_19, fp_24, fp_25, gp_3, gp_10, \
                         gp_18, gp_24, gp_31, gp_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_18 * ab_x[k] * fp_3[k]
                  + f_18 * ab_y[k] * fp_4[k]
                  - f_18 * ab_x[k] * fp_18[k]
                  + f_18 * ab_y[k] * fp_19[k]
                  + f_19 * ab_x[k] * fp_24[k]
                  - f_19 * ab_y[k] * fp_25[k]
                  - f_18 * gp_3[k]
                  + f_18 * gp_10[k]
                  - f_18 * gp_18[k]
                  + f_19 * gp_24[k]
                  + f_18 * gp_31[k]
                  - f_19 * gp_37[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_7, fp_8, fp_22, fp_23, fp_28, fp_29, gp_7, gp_14, \
                         gp_22, gp_28, gp_35, gp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_20 * ab_x[k] * fp_7[k]
                  - f_20 * ab_x[k] * fp_22[k]
                  + f_21 * ab_x[k] * fp_28[k]
                  - f_20 * gp_7[k]
                  - f_20 * gp_22[k]
                  + f_21 * gp_28[k];

        g_16[k] = -f_20 * ab_y[k] * fp_8[k]
                  - f_20 * ab_y[k] * fp_23[k]
                  + f_21 * ab_y[k] * fp_29[k]
                  - f_20 * gp_14[k]
                  - f_20 * gp_35[k]
                  + f_21 * gp_41[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fp_6, fp_7, fp_8, fp_21, fp_22, fp_23, fp_27, \
                         fp_28, fp_29, gp_6, gp_13, gp_17, gp_21, gp_27, gp_34, gp_38, gp_40, \
                         gp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = 0.75 * ab_x[k] * fp_6[k]
                  + 0.75 * ab_y[k] * fp_7[k]
                  - 1.5 * ab_z[k] * fp_8[k]
                  + 0.75 * ab_x[k] * fp_21[k]
                  + 0.75 * ab_y[k] * fp_22[k]
                  - 1.5 * ab_z[k] * fp_23[k]
                  - 0.5 * ab_x[k] * fp_27[k]
                  - 0.5 * ab_y[k] * fp_28[k]
                  + ab_z[k] * fp_29[k]
                  + 0.75 * gp_6[k]
                  + 0.75 * gp_13[k]
                  - 1.5 * gp_17[k]
                  + 0.75 * gp_21[k]
                  - 0.5 * gp_27[k]
                  + 0.75 * gp_34[k]
                  - 1.5 * gp_38[k]
                  - 0.5 * gp_40[k]
                  + gp_44[k];
    }

#pragma omp simd aligned(ab_x, fp_8, fp_23, fp_29, gp_8, gp_23, gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_20 * ab_x[k] * fp_8[k]
                  - f_20 * ab_x[k] * fp_23[k]
                  + f_21 * ab_x[k] * fp_29[k]
                  - f_20 * gp_8[k]
                  - f_20 * gp_23[k]
                  + f_21 * gp_29[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_6, fp_7, fp_21, fp_22, fp_27, fp_28, gp_6, gp_13, \
                         gp_21, gp_27, gp_34, gp_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_22 * ab_x[k] * fp_6[k]
                  + f_22 * ab_y[k] * fp_7[k]
                  - f_22 * ab_x[k] * fp_21[k]
                  + f_22 * ab_y[k] * fp_22[k]
                  + f_23 * ab_x[k] * fp_27[k]
                  - f_23 * ab_y[k] * fp_28[k]
                  - f_22 * gp_6[k]
                  + f_22 * gp_13[k]
                  - f_22 * gp_21[k]
                  + f_23 * gp_27[k]
                  + f_22 * gp_34[k]
                  - f_23 * gp_40[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_1, fp_2, fp_10, fp_11, fp_16, fp_17, gp_1, gp_5, \
                         gp_10, gp_16, gp_20, gp_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_12 * ab_x[k] * fp_1[k]
                  - f_12 * ab_x[k] * fp_10[k]
                  + f_13 * ab_x[k] * fp_16[k]
                  - f_12 * gp_1[k]
                  - f_12 * gp_10[k]
                  + f_13 * gp_16[k];

        g_21[k] = -f_12 * ab_y[k] * fp_2[k]
                  - f_12 * ab_y[k] * fp_11[k]
                  + f_13 * ab_y[k] * fp_17[k]
                  - f_12 * gp_5[k]
                  - f_12 * gp_20[k]
                  + f_13 * gp_26[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fp_0, fp_1, fp_2, fp_9, fp_10, fp_11, fp_15, fp_16, \
                         fp_17, gp_0, gp_4, gp_8, gp_9, gp_15, gp_19, gp_23, gp_25, \
                         gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_14 * ab_x[k] * fp_0[k]
                  + f_14 * ab_y[k] * fp_1[k]
                  - f_15 * ab_z[k] * fp_2[k]
                  + f_14 * ab_x[k] * fp_9[k]
                  + f_14 * ab_y[k] * fp_10[k]
                  - f_15 * ab_z[k] * fp_11[k]
                  - f_16 * ab_x[k] * fp_15[k]
                  - f_16 * ab_y[k] * fp_16[k]
                  + f_17 * ab_z[k] * fp_17[k]
                  + f_14 * gp_0[k]
                  + f_14 * gp_4[k]
                  - f_15 * gp_8[k]
                  + f_14 * gp_9[k]
                  - f_16 * gp_15[k]
                  + f_14 * gp_19[k]
                  - f_15 * gp_23[k]
                  - f_16 * gp_25[k]
                  + f_17 * gp_29[k];
    }

#pragma omp simd aligned(ab_x, fp_2, fp_11, fp_17, gp_2, gp_11, gp_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_12 * ab_x[k] * fp_2[k]
                  - f_12 * ab_x[k] * fp_11[k]
                  + f_13 * ab_x[k] * fp_17[k]
                  - f_12 * gp_2[k]
                  - f_12 * gp_11[k]
                  + f_13 * gp_17[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_0, fp_1, fp_9, fp_10, fp_15, fp_16, gp_0, gp_4, gp_9, \
                         gp_15, gp_19, gp_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_18 * ab_x[k] * fp_0[k]
                  + f_18 * ab_y[k] * fp_1[k]
                  - f_18 * ab_x[k] * fp_9[k]
                  + f_18 * ab_y[k] * fp_10[k]
                  + f_19 * ab_x[k] * fp_15[k]
                  - f_19 * ab_y[k] * fp_16[k]
                  - f_18 * gp_0[k]
                  + f_18 * gp_4[k]
                  - f_18 * gp_9[k]
                  + f_19 * gp_15[k]
                  + f_18 * gp_19[k]
                  - f_19 * gp_25[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_7, fp_8, fp_22, fp_23, gp_7, gp_14, gp_22, \
                         gp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_11 * ab_x[k] * fp_7[k]
                  - f_11 * ab_x[k] * fp_22[k]
                  + f_11 * gp_7[k]
                  - f_11 * gp_22[k];

        g_26[k] = f_11 * ab_y[k] * fp_8[k]
                  - f_11 * ab_y[k] * fp_23[k]
                  + f_11 * gp_14[k]
                  - f_11 * gp_35[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fp_6, fp_7, fp_8, fp_21, fp_22, fp_23, gp_6, gp_13, \
                         gp_17, gp_21, gp_34, gp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_24 * ab_x[k] * fp_6[k]
                  - f_24 * ab_y[k] * fp_7[k]
                  + f_9 * ab_z[k] * fp_8[k]
                  + f_24 * ab_x[k] * fp_21[k]
                  + f_24 * ab_y[k] * fp_22[k]
                  - f_9 * ab_z[k] * fp_23[k]
                  - f_24 * gp_6[k]
                  - f_24 * gp_13[k]
                  + f_9 * gp_17[k]
                  + f_24 * gp_21[k]
                  + f_24 * gp_34[k]
                  - f_9 * gp_38[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_6, fp_7, fp_8, fp_21, fp_22, fp_23, gp_6, gp_8, gp_13, \
                         gp_21, gp_23, gp_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_11 * ab_x[k] * fp_8[k]
                  - f_11 * ab_x[k] * fp_23[k]
                  + f_11 * gp_8[k]
                  - f_11 * gp_23[k];

        g_29[k] = f_25 * ab_x[k] * fp_6[k]
                  - f_25 * ab_y[k] * fp_7[k]
                  - f_25 * ab_x[k] * fp_21[k]
                  + f_25 * ab_y[k] * fp_22[k]
                  + f_25 * gp_6[k]
                  - f_25 * gp_13[k]
                  - f_25 * gp_21[k]
                  + f_25 * gp_34[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_1, fp_2, fp_10, fp_11, gp_1, gp_5, gp_10, \
                         gp_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_1 * ab_x[k] * fp_1[k]
                  - f_0 * ab_x[k] * fp_10[k]
                  + f_1 * gp_1[k]
                  - f_0 * gp_10[k];

        g_31[k] = f_1 * ab_y[k] * fp_2[k]
                  - f_0 * ab_y[k] * fp_11[k]
                  + f_1 * gp_5[k]
                  - f_0 * gp_20[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fp_0, fp_1, fp_2, fp_9, fp_10, fp_11, gp_0, gp_4, \
                         gp_8, gp_9, gp_19, gp_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_4 * ab_x[k] * fp_0[k]
                  - f_4 * ab_y[k] * fp_1[k]
                  + f_5 * ab_z[k] * fp_2[k]
                  + f_2 * ab_x[k] * fp_9[k]
                  + f_2 * ab_y[k] * fp_10[k]
                  - f_3 * ab_z[k] * fp_11[k]
                  - f_4 * gp_0[k]
                  - f_4 * gp_4[k]
                  + f_5 * gp_8[k]
                  + f_2 * gp_9[k]
                  + f_2 * gp_19[k]
                  - f_3 * gp_23[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fp_0, fp_1, fp_2, fp_9, fp_10, fp_11, gp_0, gp_2, gp_4, \
                         gp_9, gp_11, gp_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_1 * ab_x[k] * fp_2[k]
                  - f_0 * ab_x[k] * fp_11[k]
                  + f_1 * gp_2[k]
                  - f_0 * gp_11[k];

        g_34[k] = f_7 * ab_x[k] * fp_0[k]
                  - f_7 * ab_y[k] * fp_1[k]
                  - f_6 * ab_x[k] * fp_9[k]
                  + f_6 * ab_y[k] * fp_10[k]
                  + f_7 * gp_0[k]
                  - f_7 * gp_4[k]
                  - f_6 * gp_9[k]
                  + f_6 * gp_19[k];
    }
}

}  // namespace simdtrf
