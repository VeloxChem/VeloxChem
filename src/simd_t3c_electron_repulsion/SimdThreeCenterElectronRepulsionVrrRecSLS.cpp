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


#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sls_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sis0, const size_t sis1,
                                                   const size_t sks0, const size_t sks1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 3.5 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.5 * gamma / (p * q);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 2.0 * gamma / (p * q);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.5 * gamma / (p * q);
    const auto f_9 = 1.0 / p;
    const auto f_10 = gamma / (p * q);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 0.5 * gamma / (p * q);

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
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sis0_0 = buffer.data(sis0 + 0);
    const auto *sis0_3 = buffer.data(sis0 + 3);
    const auto *sis0_5 = buffer.data(sis0 + 5);
    const auto *sis0_6 = buffer.data(sis0 + 6);
    const auto *sis0_9 = buffer.data(sis0 + 9);
    const auto *sis0_10 = buffer.data(sis0 + 10);
    const auto *sis0_12 = buffer.data(sis0 + 12);
    const auto *sis0_14 = buffer.data(sis0 + 14);
    const auto *sis0_15 = buffer.data(sis0 + 15);
    const auto *sis0_17 = buffer.data(sis0 + 17);
    const auto *sis0_18 = buffer.data(sis0 + 18);
    const auto *sis0_20 = buffer.data(sis0 + 20);
    const auto *sis0_21 = buffer.data(sis0 + 21);
    const auto *sis0_23 = buffer.data(sis0 + 23);
    const auto *sis0_24 = buffer.data(sis0 + 24);
    const auto *sis0_25 = buffer.data(sis0 + 25);
    const auto *sis0_26 = buffer.data(sis0 + 26);
    const auto *sis0_27 = buffer.data(sis0 + 27);

    const auto *sis1_0 = buffer.data(sis1 + 0);
    const auto *sis1_3 = buffer.data(sis1 + 3);
    const auto *sis1_5 = buffer.data(sis1 + 5);
    const auto *sis1_6 = buffer.data(sis1 + 6);
    const auto *sis1_9 = buffer.data(sis1 + 9);
    const auto *sis1_10 = buffer.data(sis1 + 10);
    const auto *sis1_12 = buffer.data(sis1 + 12);
    const auto *sis1_14 = buffer.data(sis1 + 14);
    const auto *sis1_15 = buffer.data(sis1 + 15);
    const auto *sis1_17 = buffer.data(sis1 + 17);
    const auto *sis1_18 = buffer.data(sis1 + 18);
    const auto *sis1_20 = buffer.data(sis1 + 20);
    const auto *sis1_21 = buffer.data(sis1 + 21);
    const auto *sis1_23 = buffer.data(sis1 + 23);
    const auto *sis1_24 = buffer.data(sis1 + 24);
    const auto *sis1_25 = buffer.data(sis1 + 25);
    const auto *sis1_26 = buffer.data(sis1 + 26);
    const auto *sis1_27 = buffer.data(sis1 + 27);

    const auto *sks0_0 = buffer.data(sks0 + 0);
    const auto *sks0_2 = buffer.data(sks0 + 2);
    const auto *sks0_3 = buffer.data(sks0 + 3);
    const auto *sks0_5 = buffer.data(sks0 + 5);
    const auto *sks0_6 = buffer.data(sks0 + 6);
    const auto *sks0_9 = buffer.data(sks0 + 9);
    const auto *sks0_10 = buffer.data(sks0 + 10);
    const auto *sks0_12 = buffer.data(sks0 + 12);
    const auto *sks0_14 = buffer.data(sks0 + 14);
    const auto *sks0_15 = buffer.data(sks0 + 15);
    const auto *sks0_17 = buffer.data(sks0 + 17);
    const auto *sks0_18 = buffer.data(sks0 + 18);
    const auto *sks0_20 = buffer.data(sks0 + 20);
    const auto *sks0_21 = buffer.data(sks0 + 21);
    const auto *sks0_23 = buffer.data(sks0 + 23);
    const auto *sks0_24 = buffer.data(sks0 + 24);
    const auto *sks0_25 = buffer.data(sks0 + 25);
    const auto *sks0_27 = buffer.data(sks0 + 27);
    const auto *sks0_28 = buffer.data(sks0 + 28);
    const auto *sks0_29 = buffer.data(sks0 + 29);
    const auto *sks0_30 = buffer.data(sks0 + 30);
    const auto *sks0_31 = buffer.data(sks0 + 31);
    const auto *sks0_32 = buffer.data(sks0 + 32);
    const auto *sks0_33 = buffer.data(sks0 + 33);
    const auto *sks0_34 = buffer.data(sks0 + 34);
    const auto *sks0_35 = buffer.data(sks0 + 35);

    const auto *sks1_0 = buffer.data(sks1 + 0);
    const auto *sks1_2 = buffer.data(sks1 + 2);
    const auto *sks1_3 = buffer.data(sks1 + 3);
    const auto *sks1_5 = buffer.data(sks1 + 5);
    const auto *sks1_6 = buffer.data(sks1 + 6);
    const auto *sks1_9 = buffer.data(sks1 + 9);
    const auto *sks1_10 = buffer.data(sks1 + 10);
    const auto *sks1_12 = buffer.data(sks1 + 12);
    const auto *sks1_14 = buffer.data(sks1 + 14);
    const auto *sks1_15 = buffer.data(sks1 + 15);
    const auto *sks1_17 = buffer.data(sks1 + 17);
    const auto *sks1_18 = buffer.data(sks1 + 18);
    const auto *sks1_20 = buffer.data(sks1 + 20);
    const auto *sks1_21 = buffer.data(sks1 + 21);
    const auto *sks1_23 = buffer.data(sks1 + 23);
    const auto *sks1_24 = buffer.data(sks1 + 24);
    const auto *sks1_25 = buffer.data(sks1 + 25);
    const auto *sks1_27 = buffer.data(sks1 + 27);
    const auto *sks1_28 = buffer.data(sks1 + 28);
    const auto *sks1_29 = buffer.data(sks1 + 29);
    const auto *sks1_30 = buffer.data(sks1 + 30);
    const auto *sks1_31 = buffer.data(sks1 + 31);
    const auto *sks1_32 = buffer.data(sks1 + 32);
    const auto *sks1_33 = buffer.data(sks1 + 33);
    const auto *sks1_34 = buffer.data(sks1 + 34);
    const auto *sks1_35 = buffer.data(sks1 + 35);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sis0_0, sis1_0, \
                         sks0_0, sks1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sis0_0[k]
                 - f_1 * sis1_0[k]
                 + pb_x[k] * sks0_0[k]
                 - f_2 * pc_x[k] * sks1_0[k];

        t_1[k] = pb_y[k] * sks0_0[k]
                 - f_2 * pc_y[k] * sks1_0[k];

        t_2[k] = pb_z[k] * sks0_0[k]
                 - f_2 * pc_z[k] * sks1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sis0_3, sis1_3, sks0_2, sks0_3, \
                         sks1_2, sks1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * sis0_3[k]
                 - f_4 * sis1_3[k]
                 + pb_x[k] * sks0_3[k]
                 - f_2 * pc_x[k] * sks1_3[k];

        t_4[k] = pb_y[k] * sks0_2[k]
                 - f_2 * pc_y[k] * sks1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pb_x, pc_x, sis0_5, sis0_6, sis1_5, sis1_6, sks0_5, sks0_6, \
                         sks1_5, sks1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * sis0_5[k]
                 - f_4 * sis1_5[k]
                 + pb_x[k] * sks0_5[k]
                 - f_2 * pc_x[k] * sks1_5[k];

        t_6[k] = f_5 * sis0_6[k]
                 - f_6 * sis1_6[k]
                 + pb_x[k] * sks0_6[k]
                 - f_2 * pc_x[k] * sks1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, pc_y, pc_z, sks0_3, sks0_5, sks1_3, \
                         sks1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * sks0_3[k]
                 - f_2 * pc_z[k] * sks1_3[k];

        t_8[k] = pb_y[k] * sks0_5[k]
                 - f_2 * pc_y[k] * sks1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_x, pc_x, sis0_9, sis0_10, sis1_9, sis1_10, sks0_9, \
                         sks0_10, sks1_9, sks1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * sis0_9[k]
                 - f_6 * sis1_9[k]
                 + pb_x[k] * sks0_9[k]
                 - f_2 * pc_x[k] * sks1_9[k];

        t_10[k] = f_7 * sis0_10[k]
                  - f_8 * sis1_10[k]
                  + pb_x[k] * sks0_10[k]
                  - f_2 * pc_x[k] * sks1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pb_x, pb_z, pc_x, pc_z, sis0_12, sis1_12, sks0_6, \
                         sks0_12, sks1_6, sks1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * sks0_6[k]
                  - f_2 * pc_z[k] * sks1_6[k];

        t_12[k] = f_7 * sis0_12[k]
                  - f_8 * sis1_12[k]
                  + pb_x[k] * sks0_12[k]
                  - f_2 * pc_x[k] * sks1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_x, pb_y, pc_x, pc_y, sis0_14, sis1_14, sks0_9, \
                         sks0_14, sks1_9, sks1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * sks0_9[k]
                  - f_2 * pc_y[k] * sks1_9[k];

        t_14[k] = f_7 * sis0_14[k]
                  - f_8 * sis1_14[k]
                  + pb_x[k] * sks0_14[k]
                  - f_2 * pc_x[k] * sks1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_x, pb_z, pc_x, pc_z, sis0_15, sis1_15, sks0_10, \
                         sks0_15, sks1_10, sks1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * sis0_15[k]
                  - f_10 * sis1_15[k]
                  + pb_x[k] * sks0_15[k]
                  - f_2 * pc_x[k] * sks1_15[k];

        t_16[k] = pb_z[k] * sks0_10[k]
                  - f_2 * pc_z[k] * sks1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pb_x, pc_x, sis0_17, sis0_18, sis1_17, sis1_18, sks0_17, \
                         sks0_18, sks1_17, sks1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * sis0_17[k]
                  - f_10 * sis1_17[k]
                  + pb_x[k] * sks0_17[k]
                  - f_2 * pc_x[k] * sks1_17[k];

        t_18[k] = f_9 * sis0_18[k]
                  - f_10 * sis1_18[k]
                  + pb_x[k] * sks0_18[k]
                  - f_2 * pc_x[k] * sks1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_x, pb_y, pc_x, pc_y, sis0_20, sis1_20, sks0_14, \
                         sks0_20, sks1_14, sks1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * sks0_14[k]
                  - f_2 * pc_y[k] * sks1_14[k];

        t_20[k] = f_9 * sis0_20[k]
                  - f_10 * sis1_20[k]
                  + pb_x[k] * sks0_20[k]
                  - f_2 * pc_x[k] * sks1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, pb_z, pc_x, pc_z, sis0_21, sis1_21, sks0_15, \
                         sks0_21, sks1_15, sks1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * sis0_21[k]
                  - f_12 * sis1_21[k]
                  + pb_x[k] * sks0_21[k]
                  - f_2 * pc_x[k] * sks1_21[k];

        t_22[k] = pb_z[k] * sks0_15[k]
                  - f_2 * pc_z[k] * sks1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pb_x, pc_x, sis0_23, sis0_24, sis1_23, sis1_24, sks0_23, \
                         sks0_24, sks1_23, sks1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * sis0_23[k]
                  - f_12 * sis1_23[k]
                  + pb_x[k] * sks0_23[k]
                  - f_2 * pc_x[k] * sks1_23[k];

        t_24[k] = f_11 * sis0_24[k]
                  - f_12 * sis1_24[k]
                  + pb_x[k] * sks0_24[k]
                  - f_2 * pc_x[k] * sks1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_x, pb_y, pc_x, pc_y, sis0_25, sis1_25, sks0_20, \
                         sks0_25, sks1_20, sks1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * sis0_25[k]
                  - f_12 * sis1_25[k]
                  + pb_x[k] * sks0_25[k]
                  - f_2 * pc_x[k] * sks1_25[k];

        t_26[k] = pb_y[k] * sks0_20[k]
                  - f_2 * pc_y[k] * sks1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pc_x, sis0_27, sis1_27, sks0_27, \
                         sks0_28, sks0_29, sks0_30, sks1_27, sks1_28, sks1_29, \
                         sks1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * sis0_27[k]
                  - f_12 * sis1_27[k]
                  + pb_x[k] * sks0_27[k]
                  - f_2 * pc_x[k] * sks1_27[k];

        t_28[k] = pb_x[k] * sks0_28[k]
                  - f_2 * pc_x[k] * sks1_28[k];

        t_29[k] = pb_x[k] * sks0_29[k]
                  - f_2 * pc_x[k] * sks1_29[k];

        t_30[k] = pb_x[k] * sks0_30[k]
                  - f_2 * pc_x[k] * sks1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_x, pc_x, sks0_31, sks0_32, sks0_33, \
                         sks0_34, sks1_31, sks1_32, sks1_33, sks1_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_x[k] * sks0_31[k]
                  - f_2 * pc_x[k] * sks1_31[k];

        t_32[k] = pb_x[k] * sks0_32[k]
                  - f_2 * pc_x[k] * sks1_32[k];

        t_33[k] = pb_x[k] * sks0_33[k]
                  - f_2 * pc_x[k] * sks1_33[k];

        t_34[k] = pb_x[k] * sks0_34[k]
                  - f_2 * pc_x[k] * sks1_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sis0_21, \
                         sis1_21, sks0_28, sks0_35, sks1_28, sks1_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_x[k] * sks0_35[k]
                  - f_2 * pc_x[k] * sks1_35[k];

        t_36[k] = f_0 * sis0_21[k]
                  - f_1 * sis1_21[k]
                  + pb_y[k] * sks0_28[k]
                  - f_2 * pc_y[k] * sks1_28[k];

        t_37[k] = pb_z[k] * sks0_28[k]
                  - f_2 * pc_z[k] * sks1_28[k];
    }

#pragma omp simd aligned(t_38, t_39, pb_y, pc_y, sis0_23, sis0_24, sis1_23, sis1_24, sks0_30, \
                         sks0_31, sks1_30, sks1_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * sis0_23[k]
                  - f_4 * sis1_23[k]
                  + pb_y[k] * sks0_30[k]
                  - f_2 * pc_y[k] * sks1_30[k];

        t_39[k] = f_5 * sis0_24[k]
                  - f_6 * sis1_24[k]
                  + pb_y[k] * sks0_31[k]
                  - f_2 * pc_y[k] * sks1_31[k];
    }

#pragma omp simd aligned(t_40, t_41, pb_y, pc_y, sis0_25, sis0_26, sis1_25, sis1_26, sks0_32, \
                         sks0_33, sks1_32, sks1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * sis0_25[k]
                  - f_8 * sis1_25[k]
                  + pb_y[k] * sks0_32[k]
                  - f_2 * pc_y[k] * sks1_32[k];

        t_41[k] = f_9 * sis0_26[k]
                  - f_10 * sis1_26[k]
                  + pb_y[k] * sks0_33[k]
                  - f_2 * pc_y[k] * sks1_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_y, pb_z, pc_y, pc_z, sis0_27, sis1_27, sks0_34, \
                         sks0_35, sks1_34, sks1_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_11 * sis0_27[k]
                  - f_12 * sis1_27[k]
                  + pb_y[k] * sks0_34[k]
                  - f_2 * pc_y[k] * sks1_34[k];

        t_43[k] = pb_y[k] * sks0_35[k]
                  - f_2 * pc_y[k] * sks1_35[k];

        t_44[k] = f_0 * sis0_27[k]
                  - f_1 * sis1_27[k]
                  + pb_z[k] * sks0_35[k]
                  - f_2 * pc_z[k] * sks1_35[k];
    }
}

}  // namespace simdt3ceri
