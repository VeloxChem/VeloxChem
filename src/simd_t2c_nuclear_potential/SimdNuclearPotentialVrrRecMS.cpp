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


#include "SimdNuclearPotentialVrrRecMS.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_ms_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                    const size_t pc, const size_t ks0, const size_t ks1,
                                    const size_t ls0, const size_t ls1, const size_t ncols,
                                    const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 3.0 / p;
    const auto f_2 = 2.5 / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;

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
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);
    const auto *ks0_18 = buffer.data(ks0 + 18);
    const auto *ks0_20 = buffer.data(ks0 + 20);
    const auto *ks0_21 = buffer.data(ks0 + 21);
    const auto *ks0_23 = buffer.data(ks0 + 23);
    const auto *ks0_24 = buffer.data(ks0 + 24);
    const auto *ks0_25 = buffer.data(ks0 + 25);
    const auto *ks0_27 = buffer.data(ks0 + 27);
    const auto *ks0_28 = buffer.data(ks0 + 28);
    const auto *ks0_30 = buffer.data(ks0 + 30);
    const auto *ks0_31 = buffer.data(ks0 + 31);
    const auto *ks0_32 = buffer.data(ks0 + 32);
    const auto *ks0_33 = buffer.data(ks0 + 33);
    const auto *ks0_34 = buffer.data(ks0 + 34);
    const auto *ks0_35 = buffer.data(ks0 + 35);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);
    const auto *ks1_18 = buffer.data(ks1 + 18);
    const auto *ks1_20 = buffer.data(ks1 + 20);
    const auto *ks1_21 = buffer.data(ks1 + 21);
    const auto *ks1_23 = buffer.data(ks1 + 23);
    const auto *ks1_24 = buffer.data(ks1 + 24);
    const auto *ks1_25 = buffer.data(ks1 + 25);
    const auto *ks1_27 = buffer.data(ks1 + 27);
    const auto *ks1_28 = buffer.data(ks1 + 28);
    const auto *ks1_30 = buffer.data(ks1 + 30);
    const auto *ks1_31 = buffer.data(ks1 + 31);
    const auto *ks1_32 = buffer.data(ks1 + 32);
    const auto *ks1_33 = buffer.data(ks1 + 33);
    const auto *ks1_34 = buffer.data(ks1 + 34);
    const auto *ks1_35 = buffer.data(ks1 + 35);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_17 = buffer.data(ls0 + 17);
    const auto *ls0_18 = buffer.data(ls0 + 18);
    const auto *ls0_20 = buffer.data(ls0 + 20);
    const auto *ls0_21 = buffer.data(ls0 + 21);
    const auto *ls0_23 = buffer.data(ls0 + 23);
    const auto *ls0_24 = buffer.data(ls0 + 24);
    const auto *ls0_25 = buffer.data(ls0 + 25);
    const auto *ls0_27 = buffer.data(ls0 + 27);
    const auto *ls0_28 = buffer.data(ls0 + 28);
    const auto *ls0_30 = buffer.data(ls0 + 30);
    const auto *ls0_31 = buffer.data(ls0 + 31);
    const auto *ls0_32 = buffer.data(ls0 + 32);
    const auto *ls0_33 = buffer.data(ls0 + 33);
    const auto *ls0_35 = buffer.data(ls0 + 35);
    const auto *ls0_36 = buffer.data(ls0 + 36);
    const auto *ls0_37 = buffer.data(ls0 + 37);
    const auto *ls0_38 = buffer.data(ls0 + 38);
    const auto *ls0_39 = buffer.data(ls0 + 39);
    const auto *ls0_40 = buffer.data(ls0 + 40);
    const auto *ls0_41 = buffer.data(ls0 + 41);
    const auto *ls0_42 = buffer.data(ls0 + 42);
    const auto *ls0_43 = buffer.data(ls0 + 43);
    const auto *ls0_44 = buffer.data(ls0 + 44);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_17 = buffer.data(ls1 + 17);
    const auto *ls1_18 = buffer.data(ls1 + 18);
    const auto *ls1_20 = buffer.data(ls1 + 20);
    const auto *ls1_21 = buffer.data(ls1 + 21);
    const auto *ls1_23 = buffer.data(ls1 + 23);
    const auto *ls1_24 = buffer.data(ls1 + 24);
    const auto *ls1_25 = buffer.data(ls1 + 25);
    const auto *ls1_27 = buffer.data(ls1 + 27);
    const auto *ls1_28 = buffer.data(ls1 + 28);
    const auto *ls1_30 = buffer.data(ls1 + 30);
    const auto *ls1_31 = buffer.data(ls1 + 31);
    const auto *ls1_32 = buffer.data(ls1 + 32);
    const auto *ls1_33 = buffer.data(ls1 + 33);
    const auto *ls1_35 = buffer.data(ls1 + 35);
    const auto *ls1_36 = buffer.data(ls1 + 36);
    const auto *ls1_37 = buffer.data(ls1 + 37);
    const auto *ls1_38 = buffer.data(ls1 + 38);
    const auto *ls1_39 = buffer.data(ls1 + 39);
    const auto *ls1_40 = buffer.data(ls1 + 40);
    const auto *ls1_41 = buffer.data(ls1 + 41);
    const auto *ls1_42 = buffer.data(ls1 + 42);
    const auto *ls1_43 = buffer.data(ls1 + 43);
    const auto *ls1_44 = buffer.data(ls1 + 44);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ks0_0, ks1_0, \
                         ls0_0, ls1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ks0_0[k]
                 - f_0 * ks1_0[k]
                 + pa_x[k] * ls0_0[k]
                 - pc_x[k] * ls1_0[k];

        t_1[k] = pa_y[k] * ls0_0[k]
                 - pc_y[k] * ls1_0[k];

        t_2[k] = pa_z[k] * ls0_0[k]
                 - pc_z[k] * ls1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, ks0_3, ks1_3, ls0_2, ls0_3, ls1_2, \
                         ls1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * ks0_3[k]
                 - f_1 * ks1_3[k]
                 + pa_x[k] * ls0_3[k]
                 - pc_x[k] * ls1_3[k];

        t_4[k] = pa_y[k] * ls0_2[k]
                 - pc_y[k] * ls1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, ks0_5, ks0_6, ks1_5, ks1_6, ls0_5, ls0_6, \
                         ls1_5, ls1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ks0_5[k]
                 - f_1 * ks1_5[k]
                 + pa_x[k] * ls0_5[k]
                 - pc_x[k] * ls1_5[k];

        t_6[k] = f_2 * ks0_6[k]
                 - f_2 * ks1_6[k]
                 + pa_x[k] * ls0_6[k]
                 - pc_x[k] * ls1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, ls0_3, ls0_5, ls1_3, \
                         ls1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * ls0_3[k]
                 - pc_z[k] * ls1_3[k];

        t_8[k] = pa_y[k] * ls0_5[k]
                 - pc_y[k] * ls1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pc_x, ks0_9, ks0_10, ks1_9, ks1_10, ls0_9, ls0_10, \
                         ls1_9, ls1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * ks0_9[k]
                 - f_2 * ks1_9[k]
                 + pa_x[k] * ls0_9[k]
                 - pc_x[k] * ls1_9[k];

        t_10[k] = f_3 * ks0_10[k]
                  - f_3 * ks1_10[k]
                  + pa_x[k] * ls0_10[k]
                  - pc_x[k] * ls1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, pc_x, pc_z, ks0_12, ks1_12, ls0_6, ls0_12, \
                         ls1_6, ls1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * ls0_6[k]
                  - pc_z[k] * ls1_6[k];

        t_12[k] = f_3 * ks0_12[k]
                  - f_3 * ks1_12[k]
                  + pa_x[k] * ls0_12[k]
                  - pc_x[k] * ls1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, pc_x, pc_y, ks0_14, ks1_14, ls0_9, ls0_14, \
                         ls1_9, ls1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * ls0_9[k]
                  - pc_y[k] * ls1_9[k];

        t_14[k] = f_3 * ks0_14[k]
                  - f_3 * ks1_14[k]
                  + pa_x[k] * ls0_14[k]
                  - pc_x[k] * ls1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_x, pa_z, pc_x, pc_z, ks0_15, ks1_15, ls0_10, ls0_15, \
                         ls1_10, ls1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * ks0_15[k]
                  - f_4 * ks1_15[k]
                  + pa_x[k] * ls0_15[k]
                  - pc_x[k] * ls1_15[k];

        t_16[k] = pa_z[k] * ls0_10[k]
                  - pc_z[k] * ls1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pc_x, ks0_17, ks0_18, ks1_17, ks1_18, ls0_17, \
                         ls0_18, ls1_17, ls1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * ks0_17[k]
                  - f_4 * ks1_17[k]
                  + pa_x[k] * ls0_17[k]
                  - pc_x[k] * ls1_17[k];

        t_18[k] = f_4 * ks0_18[k]
                  - f_4 * ks1_18[k]
                  + pa_x[k] * ls0_18[k]
                  - pc_x[k] * ls1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_x, pa_y, pc_x, pc_y, ks0_20, ks1_20, ls0_14, ls0_20, \
                         ls1_14, ls1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * ls0_14[k]
                  - pc_y[k] * ls1_14[k];

        t_20[k] = f_4 * ks0_20[k]
                  - f_4 * ks1_20[k]
                  + pa_x[k] * ls0_20[k]
                  - pc_x[k] * ls1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pa_x, pa_z, pc_x, pc_z, ks0_21, ks1_21, ls0_15, ls0_21, \
                         ls1_15, ls1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * ks0_21[k]
                  - f_5 * ks1_21[k]
                  + pa_x[k] * ls0_21[k]
                  - pc_x[k] * ls1_21[k];

        t_22[k] = pa_z[k] * ls0_15[k]
                  - pc_z[k] * ls1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_x, pc_x, ks0_23, ks0_24, ks1_23, ks1_24, ls0_23, \
                         ls0_24, ls1_23, ls1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * ks0_23[k]
                  - f_5 * ks1_23[k]
                  + pa_x[k] * ls0_23[k]
                  - pc_x[k] * ls1_23[k];

        t_24[k] = f_5 * ks0_24[k]
                  - f_5 * ks1_24[k]
                  + pa_x[k] * ls0_24[k]
                  - pc_x[k] * ls1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pa_x, pa_y, pc_x, pc_y, ks0_25, ks1_25, ls0_20, ls0_25, \
                         ls1_20, ls1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * ks0_25[k]
                  - f_5 * ks1_25[k]
                  + pa_x[k] * ls0_25[k]
                  - pc_x[k] * ls1_25[k];

        t_26[k] = pa_y[k] * ls0_20[k]
                  - pc_y[k] * ls1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pa_x, pc_x, ks0_27, ks0_28, ks1_27, ks1_28, ls0_27, \
                         ls0_28, ls1_27, ls1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * ks0_27[k]
                  - f_5 * ks1_27[k]
                  + pa_x[k] * ls0_27[k]
                  - pc_x[k] * ls1_27[k];

        t_28[k] = f_6 * ks0_28[k]
                  - f_6 * ks1_28[k]
                  + pa_x[k] * ls0_28[k]
                  - pc_x[k] * ls1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_x, pa_z, pc_x, pc_z, ks0_30, ks1_30, ls0_21, ls0_30, \
                         ls1_21, ls1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * ls0_21[k]
                  - pc_z[k] * ls1_21[k];

        t_30[k] = f_6 * ks0_30[k]
                  - f_6 * ks1_30[k]
                  + pa_x[k] * ls0_30[k]
                  - pc_x[k] * ls1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pa_x, pc_x, ks0_31, ks0_32, ks1_31, ks1_32, ls0_31, \
                         ls0_32, ls1_31, ls1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * ks0_31[k]
                  - f_6 * ks1_31[k]
                  + pa_x[k] * ls0_31[k]
                  - pc_x[k] * ls1_31[k];

        t_32[k] = f_6 * ks0_32[k]
                  - f_6 * ks1_32[k]
                  + pa_x[k] * ls0_32[k]
                  - pc_x[k] * ls1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pa_x, pa_y, pc_x, pc_y, ks0_33, ks1_33, ls0_27, ls0_33, \
                         ls1_27, ls1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_6 * ks0_33[k]
                  - f_6 * ks1_33[k]
                  + pa_x[k] * ls0_33[k]
                  - pc_x[k] * ls1_33[k];

        t_34[k] = pa_y[k] * ls0_27[k]
                  - pc_y[k] * ls1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pc_x, ks0_35, ks1_35, ls0_35, ls0_36, \
                         ls0_37, ls0_38, ls1_35, ls1_36, ls1_37, \
                         ls1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * ks0_35[k]
                  - f_6 * ks1_35[k]
                  + pa_x[k] * ls0_35[k]
                  - pc_x[k] * ls1_35[k];

        t_36[k] = pa_x[k] * ls0_36[k]
                  - pc_x[k] * ls1_36[k];

        t_37[k] = pa_x[k] * ls0_37[k]
                  - pc_x[k] * ls1_37[k];

        t_38[k] = pa_x[k] * ls0_38[k]
                  - pc_x[k] * ls1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pc_x, ls0_39, ls0_40, ls0_41, ls0_42, \
                         ls1_39, ls1_40, ls1_41, ls1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_x[k] * ls0_39[k]
                  - pc_x[k] * ls1_39[k];

        t_40[k] = pa_x[k] * ls0_40[k]
                  - pc_x[k] * ls1_40[k];

        t_41[k] = pa_x[k] * ls0_41[k]
                  - pc_x[k] * ls1_41[k];

        t_42[k] = pa_x[k] * ls0_42[k]
                  - pc_x[k] * ls1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_y, pc_x, pc_y, ks0_28, ks1_28, ls0_36, \
                         ls0_43, ls0_44, ls1_36, ls1_43, ls1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_x[k] * ls0_43[k]
                  - pc_x[k] * ls1_43[k];

        t_44[k] = pa_x[k] * ls0_44[k]
                  - pc_x[k] * ls1_44[k];

        t_45[k] = f_0 * ks0_28[k]
                  - f_0 * ks1_28[k]
                  + pa_y[k] * ls0_36[k]
                  - pc_y[k] * ls1_36[k];
    }

#pragma omp simd aligned(t_46, t_47, pa_y, pa_z, pc_y, pc_z, ks0_30, ks1_30, ls0_36, ls0_38, \
                         ls1_36, ls1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * ls0_36[k]
                  - pc_z[k] * ls1_36[k];

        t_47[k] = f_1 * ks0_30[k]
                  - f_1 * ks1_30[k]
                  + pa_y[k] * ls0_38[k]
                  - pc_y[k] * ls1_38[k];
    }

#pragma omp simd aligned(t_48, t_49, pa_y, pc_y, ks0_31, ks0_32, ks1_31, ks1_32, ls0_39, \
                         ls0_40, ls1_39, ls1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_2 * ks0_31[k]
                  - f_2 * ks1_31[k]
                  + pa_y[k] * ls0_39[k]
                  - pc_y[k] * ls1_39[k];

        t_49[k] = f_3 * ks0_32[k]
                  - f_3 * ks1_32[k]
                  + pa_y[k] * ls0_40[k]
                  - pc_y[k] * ls1_40[k];
    }

#pragma omp simd aligned(t_50, t_51, pa_y, pc_y, ks0_33, ks0_34, ks1_33, ks1_34, ls0_41, \
                         ls0_42, ls1_41, ls1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_4 * ks0_33[k]
                  - f_4 * ks1_33[k]
                  + pa_y[k] * ls0_41[k]
                  - pc_y[k] * ls1_41[k];

        t_51[k] = f_5 * ks0_34[k]
                  - f_5 * ks1_34[k]
                  + pa_y[k] * ls0_42[k]
                  - pc_y[k] * ls1_42[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_y, pa_z, pc_y, pc_z, ks0_35, ks1_35, ls0_43, \
                         ls0_44, ls1_43, ls1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_6 * ks0_35[k]
                  - f_6 * ks1_35[k]
                  + pa_y[k] * ls0_43[k]
                  - pc_y[k] * ls1_43[k];

        t_53[k] = pa_y[k] * ls0_44[k]
                  - pc_y[k] * ls1_44[k];

        t_54[k] = f_0 * ks0_35[k]
                  - f_0 * ks1_35[k]
                  + pa_z[k] * ls0_44[k]
                  - pc_z[k] * ls1_44[k];
    }
}

}  // namespace simdnpot
