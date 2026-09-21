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


#include "SimdNuclearPotentialVrrRecNS.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_ns_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                    const size_t pc, const size_t ls0, const size_t ls1,
                                    const size_t ms0, const size_t ms1, const size_t ncols,
                                    const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 3.5 / p;
    const auto f_2 = 3.0 / p;
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 0.5 / p;

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
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ls0_0 = buffer.data(ls0 + 0);
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
    const auto *ls0_38 = buffer.data(ls0 + 38);
    const auto *ls0_39 = buffer.data(ls0 + 39);
    const auto *ls0_40 = buffer.data(ls0 + 40);
    const auto *ls0_41 = buffer.data(ls0 + 41);
    const auto *ls0_42 = buffer.data(ls0 + 42);
    const auto *ls0_43 = buffer.data(ls0 + 43);
    const auto *ls0_44 = buffer.data(ls0 + 44);

    const auto *ls1_0 = buffer.data(ls1 + 0);
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
    const auto *ls1_38 = buffer.data(ls1 + 38);
    const auto *ls1_39 = buffer.data(ls1 + 39);
    const auto *ls1_40 = buffer.data(ls1 + 40);
    const auto *ls1_41 = buffer.data(ls1 + 41);
    const auto *ls1_42 = buffer.data(ls1 + 42);
    const auto *ls1_43 = buffer.data(ls1 + 43);
    const auto *ls1_44 = buffer.data(ls1 + 44);

    const auto *ms0_0 = buffer.data(ms0 + 0);
    const auto *ms0_2 = buffer.data(ms0 + 2);
    const auto *ms0_3 = buffer.data(ms0 + 3);
    const auto *ms0_5 = buffer.data(ms0 + 5);
    const auto *ms0_6 = buffer.data(ms0 + 6);
    const auto *ms0_9 = buffer.data(ms0 + 9);
    const auto *ms0_10 = buffer.data(ms0 + 10);
    const auto *ms0_12 = buffer.data(ms0 + 12);
    const auto *ms0_14 = buffer.data(ms0 + 14);
    const auto *ms0_15 = buffer.data(ms0 + 15);
    const auto *ms0_17 = buffer.data(ms0 + 17);
    const auto *ms0_18 = buffer.data(ms0 + 18);
    const auto *ms0_20 = buffer.data(ms0 + 20);
    const auto *ms0_21 = buffer.data(ms0 + 21);
    const auto *ms0_23 = buffer.data(ms0 + 23);
    const auto *ms0_24 = buffer.data(ms0 + 24);
    const auto *ms0_25 = buffer.data(ms0 + 25);
    const auto *ms0_27 = buffer.data(ms0 + 27);
    const auto *ms0_28 = buffer.data(ms0 + 28);
    const auto *ms0_30 = buffer.data(ms0 + 30);
    const auto *ms0_31 = buffer.data(ms0 + 31);
    const auto *ms0_32 = buffer.data(ms0 + 32);
    const auto *ms0_33 = buffer.data(ms0 + 33);
    const auto *ms0_35 = buffer.data(ms0 + 35);
    const auto *ms0_36 = buffer.data(ms0 + 36);
    const auto *ms0_38 = buffer.data(ms0 + 38);
    const auto *ms0_39 = buffer.data(ms0 + 39);
    const auto *ms0_40 = buffer.data(ms0 + 40);
    const auto *ms0_41 = buffer.data(ms0 + 41);
    const auto *ms0_42 = buffer.data(ms0 + 42);
    const auto *ms0_44 = buffer.data(ms0 + 44);
    const auto *ms0_45 = buffer.data(ms0 + 45);
    const auto *ms0_46 = buffer.data(ms0 + 46);
    const auto *ms0_47 = buffer.data(ms0 + 47);
    const auto *ms0_48 = buffer.data(ms0 + 48);
    const auto *ms0_49 = buffer.data(ms0 + 49);
    const auto *ms0_50 = buffer.data(ms0 + 50);
    const auto *ms0_51 = buffer.data(ms0 + 51);
    const auto *ms0_52 = buffer.data(ms0 + 52);
    const auto *ms0_53 = buffer.data(ms0 + 53);
    const auto *ms0_54 = buffer.data(ms0 + 54);

    const auto *ms1_0 = buffer.data(ms1 + 0);
    const auto *ms1_2 = buffer.data(ms1 + 2);
    const auto *ms1_3 = buffer.data(ms1 + 3);
    const auto *ms1_5 = buffer.data(ms1 + 5);
    const auto *ms1_6 = buffer.data(ms1 + 6);
    const auto *ms1_9 = buffer.data(ms1 + 9);
    const auto *ms1_10 = buffer.data(ms1 + 10);
    const auto *ms1_12 = buffer.data(ms1 + 12);
    const auto *ms1_14 = buffer.data(ms1 + 14);
    const auto *ms1_15 = buffer.data(ms1 + 15);
    const auto *ms1_17 = buffer.data(ms1 + 17);
    const auto *ms1_18 = buffer.data(ms1 + 18);
    const auto *ms1_20 = buffer.data(ms1 + 20);
    const auto *ms1_21 = buffer.data(ms1 + 21);
    const auto *ms1_23 = buffer.data(ms1 + 23);
    const auto *ms1_24 = buffer.data(ms1 + 24);
    const auto *ms1_25 = buffer.data(ms1 + 25);
    const auto *ms1_27 = buffer.data(ms1 + 27);
    const auto *ms1_28 = buffer.data(ms1 + 28);
    const auto *ms1_30 = buffer.data(ms1 + 30);
    const auto *ms1_31 = buffer.data(ms1 + 31);
    const auto *ms1_32 = buffer.data(ms1 + 32);
    const auto *ms1_33 = buffer.data(ms1 + 33);
    const auto *ms1_35 = buffer.data(ms1 + 35);
    const auto *ms1_36 = buffer.data(ms1 + 36);
    const auto *ms1_38 = buffer.data(ms1 + 38);
    const auto *ms1_39 = buffer.data(ms1 + 39);
    const auto *ms1_40 = buffer.data(ms1 + 40);
    const auto *ms1_41 = buffer.data(ms1 + 41);
    const auto *ms1_42 = buffer.data(ms1 + 42);
    const auto *ms1_44 = buffer.data(ms1 + 44);
    const auto *ms1_45 = buffer.data(ms1 + 45);
    const auto *ms1_46 = buffer.data(ms1 + 46);
    const auto *ms1_47 = buffer.data(ms1 + 47);
    const auto *ms1_48 = buffer.data(ms1 + 48);
    const auto *ms1_49 = buffer.data(ms1 + 49);
    const auto *ms1_50 = buffer.data(ms1 + 50);
    const auto *ms1_51 = buffer.data(ms1 + 51);
    const auto *ms1_52 = buffer.data(ms1 + 52);
    const auto *ms1_53 = buffer.data(ms1 + 53);
    const auto *ms1_54 = buffer.data(ms1 + 54);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ls0_0, ls1_0, \
                         ms0_0, ms1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ls0_0[k]
                 - f_0 * ls1_0[k]
                 + pa_x[k] * ms0_0[k]
                 - pc_x[k] * ms1_0[k];

        t_1[k] = pa_y[k] * ms0_0[k]
                 - pc_y[k] * ms1_0[k];

        t_2[k] = pa_z[k] * ms0_0[k]
                 - pc_z[k] * ms1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, ls0_3, ls1_3, ms0_2, ms0_3, ms1_2, \
                         ms1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * ls0_3[k]
                 - f_1 * ls1_3[k]
                 + pa_x[k] * ms0_3[k]
                 - pc_x[k] * ms1_3[k];

        t_4[k] = pa_y[k] * ms0_2[k]
                 - pc_y[k] * ms1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, ls0_5, ls0_6, ls1_5, ls1_6, ms0_5, ms0_6, \
                         ms1_5, ms1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ls0_5[k]
                 - f_1 * ls1_5[k]
                 + pa_x[k] * ms0_5[k]
                 - pc_x[k] * ms1_5[k];

        t_6[k] = f_2 * ls0_6[k]
                 - f_2 * ls1_6[k]
                 + pa_x[k] * ms0_6[k]
                 - pc_x[k] * ms1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, ms0_3, ms0_5, ms1_3, \
                         ms1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * ms0_3[k]
                 - pc_z[k] * ms1_3[k];

        t_8[k] = pa_y[k] * ms0_5[k]
                 - pc_y[k] * ms1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pc_x, ls0_9, ls0_10, ls1_9, ls1_10, ms0_9, ms0_10, \
                         ms1_9, ms1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * ls0_9[k]
                 - f_2 * ls1_9[k]
                 + pa_x[k] * ms0_9[k]
                 - pc_x[k] * ms1_9[k];

        t_10[k] = f_3 * ls0_10[k]
                  - f_3 * ls1_10[k]
                  + pa_x[k] * ms0_10[k]
                  - pc_x[k] * ms1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, pc_x, pc_z, ls0_12, ls1_12, ms0_6, ms0_12, \
                         ms1_6, ms1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * ms0_6[k]
                  - pc_z[k] * ms1_6[k];

        t_12[k] = f_3 * ls0_12[k]
                  - f_3 * ls1_12[k]
                  + pa_x[k] * ms0_12[k]
                  - pc_x[k] * ms1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, pc_x, pc_y, ls0_14, ls1_14, ms0_9, ms0_14, \
                         ms1_9, ms1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * ms0_9[k]
                  - pc_y[k] * ms1_9[k];

        t_14[k] = f_3 * ls0_14[k]
                  - f_3 * ls1_14[k]
                  + pa_x[k] * ms0_14[k]
                  - pc_x[k] * ms1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_x, pa_z, pc_x, pc_z, ls0_15, ls1_15, ms0_10, ms0_15, \
                         ms1_10, ms1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * ls0_15[k]
                  - f_4 * ls1_15[k]
                  + pa_x[k] * ms0_15[k]
                  - pc_x[k] * ms1_15[k];

        t_16[k] = pa_z[k] * ms0_10[k]
                  - pc_z[k] * ms1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pc_x, ls0_17, ls0_18, ls1_17, ls1_18, ms0_17, \
                         ms0_18, ms1_17, ms1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * ls0_17[k]
                  - f_4 * ls1_17[k]
                  + pa_x[k] * ms0_17[k]
                  - pc_x[k] * ms1_17[k];

        t_18[k] = f_4 * ls0_18[k]
                  - f_4 * ls1_18[k]
                  + pa_x[k] * ms0_18[k]
                  - pc_x[k] * ms1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_x, pa_y, pc_x, pc_y, ls0_20, ls1_20, ms0_14, ms0_20, \
                         ms1_14, ms1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * ms0_14[k]
                  - pc_y[k] * ms1_14[k];

        t_20[k] = f_4 * ls0_20[k]
                  - f_4 * ls1_20[k]
                  + pa_x[k] * ms0_20[k]
                  - pc_x[k] * ms1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pa_x, pa_z, pc_x, pc_z, ls0_21, ls1_21, ms0_15, ms0_21, \
                         ms1_15, ms1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * ls0_21[k]
                  - f_5 * ls1_21[k]
                  + pa_x[k] * ms0_21[k]
                  - pc_x[k] * ms1_21[k];

        t_22[k] = pa_z[k] * ms0_15[k]
                  - pc_z[k] * ms1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_x, pc_x, ls0_23, ls0_24, ls1_23, ls1_24, ms0_23, \
                         ms0_24, ms1_23, ms1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * ls0_23[k]
                  - f_5 * ls1_23[k]
                  + pa_x[k] * ms0_23[k]
                  - pc_x[k] * ms1_23[k];

        t_24[k] = f_5 * ls0_24[k]
                  - f_5 * ls1_24[k]
                  + pa_x[k] * ms0_24[k]
                  - pc_x[k] * ms1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pa_x, pa_y, pc_x, pc_y, ls0_25, ls1_25, ms0_20, ms0_25, \
                         ms1_20, ms1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * ls0_25[k]
                  - f_5 * ls1_25[k]
                  + pa_x[k] * ms0_25[k]
                  - pc_x[k] * ms1_25[k];

        t_26[k] = pa_y[k] * ms0_20[k]
                  - pc_y[k] * ms1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pa_x, pc_x, ls0_27, ls0_28, ls1_27, ls1_28, ms0_27, \
                         ms0_28, ms1_27, ms1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * ls0_27[k]
                  - f_5 * ls1_27[k]
                  + pa_x[k] * ms0_27[k]
                  - pc_x[k] * ms1_27[k];

        t_28[k] = f_6 * ls0_28[k]
                  - f_6 * ls1_28[k]
                  + pa_x[k] * ms0_28[k]
                  - pc_x[k] * ms1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_x, pa_z, pc_x, pc_z, ls0_30, ls1_30, ms0_21, ms0_30, \
                         ms1_21, ms1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * ms0_21[k]
                  - pc_z[k] * ms1_21[k];

        t_30[k] = f_6 * ls0_30[k]
                  - f_6 * ls1_30[k]
                  + pa_x[k] * ms0_30[k]
                  - pc_x[k] * ms1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pa_x, pc_x, ls0_31, ls0_32, ls1_31, ls1_32, ms0_31, \
                         ms0_32, ms1_31, ms1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * ls0_31[k]
                  - f_6 * ls1_31[k]
                  + pa_x[k] * ms0_31[k]
                  - pc_x[k] * ms1_31[k];

        t_32[k] = f_6 * ls0_32[k]
                  - f_6 * ls1_32[k]
                  + pa_x[k] * ms0_32[k]
                  - pc_x[k] * ms1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pa_x, pa_y, pc_x, pc_y, ls0_33, ls1_33, ms0_27, ms0_33, \
                         ms1_27, ms1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_6 * ls0_33[k]
                  - f_6 * ls1_33[k]
                  + pa_x[k] * ms0_33[k]
                  - pc_x[k] * ms1_33[k];

        t_34[k] = pa_y[k] * ms0_27[k]
                  - pc_y[k] * ms1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, pa_x, pc_x, ls0_35, ls0_36, ls1_35, ls1_36, ms0_35, \
                         ms0_36, ms1_35, ms1_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * ls0_35[k]
                  - f_6 * ls1_35[k]
                  + pa_x[k] * ms0_35[k]
                  - pc_x[k] * ms1_35[k];

        t_36[k] = f_7 * ls0_36[k]
                  - f_7 * ls1_36[k]
                  + pa_x[k] * ms0_36[k]
                  - pc_x[k] * ms1_36[k];
    }

#pragma omp simd aligned(t_37, t_38, pa_x, pa_z, pc_x, pc_z, ls0_38, ls1_38, ms0_28, ms0_38, \
                         ms1_28, ms1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_z[k] * ms0_28[k]
                  - pc_z[k] * ms1_28[k];

        t_38[k] = f_7 * ls0_38[k]
                  - f_7 * ls1_38[k]
                  + pa_x[k] * ms0_38[k]
                  - pc_x[k] * ms1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, pa_x, pc_x, ls0_39, ls0_40, ls1_39, ls1_40, ms0_39, \
                         ms0_40, ms1_39, ms1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * ls0_39[k]
                  - f_7 * ls1_39[k]
                  + pa_x[k] * ms0_39[k]
                  - pc_x[k] * ms1_39[k];

        t_40[k] = f_7 * ls0_40[k]
                  - f_7 * ls1_40[k]
                  + pa_x[k] * ms0_40[k]
                  - pc_x[k] * ms1_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pa_x, pc_x, ls0_41, ls0_42, ls1_41, ls1_42, ms0_41, \
                         ms0_42, ms1_41, ms1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * ls0_41[k]
                  - f_7 * ls1_41[k]
                  + pa_x[k] * ms0_41[k]
                  - pc_x[k] * ms1_41[k];

        t_42[k] = f_7 * ls0_42[k]
                  - f_7 * ls1_42[k]
                  + pa_x[k] * ms0_42[k]
                  - pc_x[k] * ms1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_y, pc_x, pc_y, ls0_44, ls1_44, ms0_35, \
                         ms0_44, ms0_45, ms1_35, ms1_44, ms1_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_y[k] * ms0_35[k]
                  - pc_y[k] * ms1_35[k];

        t_44[k] = f_7 * ls0_44[k]
                  - f_7 * ls1_44[k]
                  + pa_x[k] * ms0_44[k]
                  - pc_x[k] * ms1_44[k];

        t_45[k] = pa_x[k] * ms0_45[k]
                  - pc_x[k] * ms1_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pc_x, ms0_46, ms0_47, ms0_48, ms0_49, \
                         ms1_46, ms1_47, ms1_48, ms1_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * ms0_46[k]
                  - pc_x[k] * ms1_46[k];

        t_47[k] = pa_x[k] * ms0_47[k]
                  - pc_x[k] * ms1_47[k];

        t_48[k] = pa_x[k] * ms0_48[k]
                  - pc_x[k] * ms1_48[k];

        t_49[k] = pa_x[k] * ms0_49[k]
                  - pc_x[k] * ms1_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pc_x, ms0_50, ms0_51, ms0_52, ms0_53, \
                         ms1_50, ms1_51, ms1_52, ms1_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_x[k] * ms0_50[k]
                  - pc_x[k] * ms1_50[k];

        t_51[k] = pa_x[k] * ms0_51[k]
                  - pc_x[k] * ms1_51[k];

        t_52[k] = pa_x[k] * ms0_52[k]
                  - pc_x[k] * ms1_52[k];

        t_53[k] = pa_x[k] * ms0_53[k]
                  - pc_x[k] * ms1_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ls0_36, ls1_36, \
                         ms0_45, ms0_54, ms1_45, ms1_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_x[k] * ms0_54[k]
                  - pc_x[k] * ms1_54[k];

        t_55[k] = f_0 * ls0_36[k]
                  - f_0 * ls1_36[k]
                  + pa_y[k] * ms0_45[k]
                  - pc_y[k] * ms1_45[k];

        t_56[k] = pa_z[k] * ms0_45[k]
                  - pc_z[k] * ms1_45[k];
    }

#pragma omp simd aligned(t_57, t_58, pa_y, pc_y, ls0_38, ls0_39, ls1_38, ls1_39, ms0_47, \
                         ms0_48, ms1_47, ms1_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_1 * ls0_38[k]
                  - f_1 * ls1_38[k]
                  + pa_y[k] * ms0_47[k]
                  - pc_y[k] * ms1_47[k];

        t_58[k] = f_2 * ls0_39[k]
                  - f_2 * ls1_39[k]
                  + pa_y[k] * ms0_48[k]
                  - pc_y[k] * ms1_48[k];
    }

#pragma omp simd aligned(t_59, t_60, pa_y, pc_y, ls0_40, ls0_41, ls1_40, ls1_41, ms0_49, \
                         ms0_50, ms1_49, ms1_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_3 * ls0_40[k]
                  - f_3 * ls1_40[k]
                  + pa_y[k] * ms0_49[k]
                  - pc_y[k] * ms1_49[k];

        t_60[k] = f_4 * ls0_41[k]
                  - f_4 * ls1_41[k]
                  + pa_y[k] * ms0_50[k]
                  - pc_y[k] * ms1_50[k];
    }

#pragma omp simd aligned(t_61, t_62, pa_y, pc_y, ls0_42, ls0_43, ls1_42, ls1_43, ms0_51, \
                         ms0_52, ms1_51, ms1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * ls0_42[k]
                  - f_5 * ls1_42[k]
                  + pa_y[k] * ms0_51[k]
                  - pc_y[k] * ms1_51[k];

        t_62[k] = f_6 * ls0_43[k]
                  - f_6 * ls1_43[k]
                  + pa_y[k] * ms0_52[k]
                  - pc_y[k] * ms1_52[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, ls0_44, ls1_44, ms0_53, \
                         ms0_54, ms1_53, ms1_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_7 * ls0_44[k]
                  - f_7 * ls1_44[k]
                  + pa_y[k] * ms0_53[k]
                  - pc_y[k] * ms1_53[k];

        t_64[k] = pa_y[k] * ms0_54[k]
                  - pc_y[k] * ms1_54[k];

        t_65[k] = f_0 * ls0_44[k]
                  - f_0 * ls1_44[k]
                  + pa_z[k] * ms0_54[k]
                  - pc_z[k] * ms1_54[k];
    }
}

}  // namespace simdnpot
