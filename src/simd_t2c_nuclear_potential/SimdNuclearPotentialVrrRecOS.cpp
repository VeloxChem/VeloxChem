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


#include "SimdNuclearPotentialVrrRecOS.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_os_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                    const size_t pc, const size_t ms0, const size_t ms1,
                                    const size_t ns0, const size_t ns1, const size_t ncols,
                                    const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / p;
    const auto f_1 = 4.0 / p;
    const auto f_2 = 3.5 / p;
    const auto f_3 = 3.0 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 0.5 / p;

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
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ms0_0 = buffer.data(ms0 + 0);
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
    const auto *ms0_47 = buffer.data(ms0 + 47);
    const auto *ms0_48 = buffer.data(ms0 + 48);
    const auto *ms0_49 = buffer.data(ms0 + 49);
    const auto *ms0_50 = buffer.data(ms0 + 50);
    const auto *ms0_51 = buffer.data(ms0 + 51);
    const auto *ms0_52 = buffer.data(ms0 + 52);
    const auto *ms0_53 = buffer.data(ms0 + 53);
    const auto *ms0_54 = buffer.data(ms0 + 54);

    const auto *ms1_0 = buffer.data(ms1 + 0);
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
    const auto *ms1_47 = buffer.data(ms1 + 47);
    const auto *ms1_48 = buffer.data(ms1 + 48);
    const auto *ms1_49 = buffer.data(ms1 + 49);
    const auto *ms1_50 = buffer.data(ms1 + 50);
    const auto *ms1_51 = buffer.data(ms1 + 51);
    const auto *ms1_52 = buffer.data(ms1 + 52);
    const auto *ms1_53 = buffer.data(ms1 + 53);
    const auto *ms1_54 = buffer.data(ms1 + 54);

    const auto *ns0_0 = buffer.data(ns0 + 0);
    const auto *ns0_2 = buffer.data(ns0 + 2);
    const auto *ns0_3 = buffer.data(ns0 + 3);
    const auto *ns0_5 = buffer.data(ns0 + 5);
    const auto *ns0_6 = buffer.data(ns0 + 6);
    const auto *ns0_9 = buffer.data(ns0 + 9);
    const auto *ns0_10 = buffer.data(ns0 + 10);
    const auto *ns0_12 = buffer.data(ns0 + 12);
    const auto *ns0_14 = buffer.data(ns0 + 14);
    const auto *ns0_15 = buffer.data(ns0 + 15);
    const auto *ns0_17 = buffer.data(ns0 + 17);
    const auto *ns0_18 = buffer.data(ns0 + 18);
    const auto *ns0_20 = buffer.data(ns0 + 20);
    const auto *ns0_21 = buffer.data(ns0 + 21);
    const auto *ns0_23 = buffer.data(ns0 + 23);
    const auto *ns0_24 = buffer.data(ns0 + 24);
    const auto *ns0_25 = buffer.data(ns0 + 25);
    const auto *ns0_27 = buffer.data(ns0 + 27);
    const auto *ns0_28 = buffer.data(ns0 + 28);
    const auto *ns0_30 = buffer.data(ns0 + 30);
    const auto *ns0_31 = buffer.data(ns0 + 31);
    const auto *ns0_32 = buffer.data(ns0 + 32);
    const auto *ns0_33 = buffer.data(ns0 + 33);
    const auto *ns0_35 = buffer.data(ns0 + 35);
    const auto *ns0_36 = buffer.data(ns0 + 36);
    const auto *ns0_38 = buffer.data(ns0 + 38);
    const auto *ns0_39 = buffer.data(ns0 + 39);
    const auto *ns0_40 = buffer.data(ns0 + 40);
    const auto *ns0_41 = buffer.data(ns0 + 41);
    const auto *ns0_42 = buffer.data(ns0 + 42);
    const auto *ns0_44 = buffer.data(ns0 + 44);
    const auto *ns0_45 = buffer.data(ns0 + 45);
    const auto *ns0_47 = buffer.data(ns0 + 47);
    const auto *ns0_48 = buffer.data(ns0 + 48);
    const auto *ns0_49 = buffer.data(ns0 + 49);
    const auto *ns0_50 = buffer.data(ns0 + 50);
    const auto *ns0_51 = buffer.data(ns0 + 51);
    const auto *ns0_52 = buffer.data(ns0 + 52);
    const auto *ns0_54 = buffer.data(ns0 + 54);
    const auto *ns0_55 = buffer.data(ns0 + 55);
    const auto *ns0_56 = buffer.data(ns0 + 56);
    const auto *ns0_57 = buffer.data(ns0 + 57);
    const auto *ns0_58 = buffer.data(ns0 + 58);
    const auto *ns0_59 = buffer.data(ns0 + 59);
    const auto *ns0_60 = buffer.data(ns0 + 60);
    const auto *ns0_61 = buffer.data(ns0 + 61);
    const auto *ns0_62 = buffer.data(ns0 + 62);
    const auto *ns0_63 = buffer.data(ns0 + 63);
    const auto *ns0_64 = buffer.data(ns0 + 64);
    const auto *ns0_65 = buffer.data(ns0 + 65);

    const auto *ns1_0 = buffer.data(ns1 + 0);
    const auto *ns1_2 = buffer.data(ns1 + 2);
    const auto *ns1_3 = buffer.data(ns1 + 3);
    const auto *ns1_5 = buffer.data(ns1 + 5);
    const auto *ns1_6 = buffer.data(ns1 + 6);
    const auto *ns1_9 = buffer.data(ns1 + 9);
    const auto *ns1_10 = buffer.data(ns1 + 10);
    const auto *ns1_12 = buffer.data(ns1 + 12);
    const auto *ns1_14 = buffer.data(ns1 + 14);
    const auto *ns1_15 = buffer.data(ns1 + 15);
    const auto *ns1_17 = buffer.data(ns1 + 17);
    const auto *ns1_18 = buffer.data(ns1 + 18);
    const auto *ns1_20 = buffer.data(ns1 + 20);
    const auto *ns1_21 = buffer.data(ns1 + 21);
    const auto *ns1_23 = buffer.data(ns1 + 23);
    const auto *ns1_24 = buffer.data(ns1 + 24);
    const auto *ns1_25 = buffer.data(ns1 + 25);
    const auto *ns1_27 = buffer.data(ns1 + 27);
    const auto *ns1_28 = buffer.data(ns1 + 28);
    const auto *ns1_30 = buffer.data(ns1 + 30);
    const auto *ns1_31 = buffer.data(ns1 + 31);
    const auto *ns1_32 = buffer.data(ns1 + 32);
    const auto *ns1_33 = buffer.data(ns1 + 33);
    const auto *ns1_35 = buffer.data(ns1 + 35);
    const auto *ns1_36 = buffer.data(ns1 + 36);
    const auto *ns1_38 = buffer.data(ns1 + 38);
    const auto *ns1_39 = buffer.data(ns1 + 39);
    const auto *ns1_40 = buffer.data(ns1 + 40);
    const auto *ns1_41 = buffer.data(ns1 + 41);
    const auto *ns1_42 = buffer.data(ns1 + 42);
    const auto *ns1_44 = buffer.data(ns1 + 44);
    const auto *ns1_45 = buffer.data(ns1 + 45);
    const auto *ns1_47 = buffer.data(ns1 + 47);
    const auto *ns1_48 = buffer.data(ns1 + 48);
    const auto *ns1_49 = buffer.data(ns1 + 49);
    const auto *ns1_50 = buffer.data(ns1 + 50);
    const auto *ns1_51 = buffer.data(ns1 + 51);
    const auto *ns1_52 = buffer.data(ns1 + 52);
    const auto *ns1_54 = buffer.data(ns1 + 54);
    const auto *ns1_55 = buffer.data(ns1 + 55);
    const auto *ns1_56 = buffer.data(ns1 + 56);
    const auto *ns1_57 = buffer.data(ns1 + 57);
    const auto *ns1_58 = buffer.data(ns1 + 58);
    const auto *ns1_59 = buffer.data(ns1 + 59);
    const auto *ns1_60 = buffer.data(ns1 + 60);
    const auto *ns1_61 = buffer.data(ns1 + 61);
    const auto *ns1_62 = buffer.data(ns1 + 62);
    const auto *ns1_63 = buffer.data(ns1 + 63);
    const auto *ns1_64 = buffer.data(ns1 + 64);
    const auto *ns1_65 = buffer.data(ns1 + 65);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ms0_0, ms1_0, \
                         ns0_0, ns1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ms0_0[k]
                 - f_0 * ms1_0[k]
                 + pa_x[k] * ns0_0[k]
                 - pc_x[k] * ns1_0[k];

        t_1[k] = pa_y[k] * ns0_0[k]
                 - pc_y[k] * ns1_0[k];

        t_2[k] = pa_z[k] * ns0_0[k]
                 - pc_z[k] * ns1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, ms0_3, ms1_3, ns0_2, ns0_3, ns1_2, \
                         ns1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * ms0_3[k]
                 - f_1 * ms1_3[k]
                 + pa_x[k] * ns0_3[k]
                 - pc_x[k] * ns1_3[k];

        t_4[k] = pa_y[k] * ns0_2[k]
                 - pc_y[k] * ns1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, ms0_5, ms0_6, ms1_5, ms1_6, ns0_5, ns0_6, \
                         ns1_5, ns1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ms0_5[k]
                 - f_1 * ms1_5[k]
                 + pa_x[k] * ns0_5[k]
                 - pc_x[k] * ns1_5[k];

        t_6[k] = f_2 * ms0_6[k]
                 - f_2 * ms1_6[k]
                 + pa_x[k] * ns0_6[k]
                 - pc_x[k] * ns1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, ns0_3, ns0_5, ns1_3, \
                         ns1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * ns0_3[k]
                 - pc_z[k] * ns1_3[k];

        t_8[k] = pa_y[k] * ns0_5[k]
                 - pc_y[k] * ns1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pc_x, ms0_9, ms0_10, ms1_9, ms1_10, ns0_9, ns0_10, \
                         ns1_9, ns1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * ms0_9[k]
                 - f_2 * ms1_9[k]
                 + pa_x[k] * ns0_9[k]
                 - pc_x[k] * ns1_9[k];

        t_10[k] = f_3 * ms0_10[k]
                  - f_3 * ms1_10[k]
                  + pa_x[k] * ns0_10[k]
                  - pc_x[k] * ns1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, pc_x, pc_z, ms0_12, ms1_12, ns0_6, ns0_12, \
                         ns1_6, ns1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * ns0_6[k]
                  - pc_z[k] * ns1_6[k];

        t_12[k] = f_3 * ms0_12[k]
                  - f_3 * ms1_12[k]
                  + pa_x[k] * ns0_12[k]
                  - pc_x[k] * ns1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, pc_x, pc_y, ms0_14, ms1_14, ns0_9, ns0_14, \
                         ns1_9, ns1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * ns0_9[k]
                  - pc_y[k] * ns1_9[k];

        t_14[k] = f_3 * ms0_14[k]
                  - f_3 * ms1_14[k]
                  + pa_x[k] * ns0_14[k]
                  - pc_x[k] * ns1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_x, pa_z, pc_x, pc_z, ms0_15, ms1_15, ns0_10, ns0_15, \
                         ns1_10, ns1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * ms0_15[k]
                  - f_4 * ms1_15[k]
                  + pa_x[k] * ns0_15[k]
                  - pc_x[k] * ns1_15[k];

        t_16[k] = pa_z[k] * ns0_10[k]
                  - pc_z[k] * ns1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pc_x, ms0_17, ms0_18, ms1_17, ms1_18, ns0_17, \
                         ns0_18, ns1_17, ns1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * ms0_17[k]
                  - f_4 * ms1_17[k]
                  + pa_x[k] * ns0_17[k]
                  - pc_x[k] * ns1_17[k];

        t_18[k] = f_4 * ms0_18[k]
                  - f_4 * ms1_18[k]
                  + pa_x[k] * ns0_18[k]
                  - pc_x[k] * ns1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_x, pa_y, pc_x, pc_y, ms0_20, ms1_20, ns0_14, ns0_20, \
                         ns1_14, ns1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * ns0_14[k]
                  - pc_y[k] * ns1_14[k];

        t_20[k] = f_4 * ms0_20[k]
                  - f_4 * ms1_20[k]
                  + pa_x[k] * ns0_20[k]
                  - pc_x[k] * ns1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pa_x, pa_z, pc_x, pc_z, ms0_21, ms1_21, ns0_15, ns0_21, \
                         ns1_15, ns1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * ms0_21[k]
                  - f_5 * ms1_21[k]
                  + pa_x[k] * ns0_21[k]
                  - pc_x[k] * ns1_21[k];

        t_22[k] = pa_z[k] * ns0_15[k]
                  - pc_z[k] * ns1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_x, pc_x, ms0_23, ms0_24, ms1_23, ms1_24, ns0_23, \
                         ns0_24, ns1_23, ns1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * ms0_23[k]
                  - f_5 * ms1_23[k]
                  + pa_x[k] * ns0_23[k]
                  - pc_x[k] * ns1_23[k];

        t_24[k] = f_5 * ms0_24[k]
                  - f_5 * ms1_24[k]
                  + pa_x[k] * ns0_24[k]
                  - pc_x[k] * ns1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pa_x, pa_y, pc_x, pc_y, ms0_25, ms1_25, ns0_20, ns0_25, \
                         ns1_20, ns1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * ms0_25[k]
                  - f_5 * ms1_25[k]
                  + pa_x[k] * ns0_25[k]
                  - pc_x[k] * ns1_25[k];

        t_26[k] = pa_y[k] * ns0_20[k]
                  - pc_y[k] * ns1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pa_x, pc_x, ms0_27, ms0_28, ms1_27, ms1_28, ns0_27, \
                         ns0_28, ns1_27, ns1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * ms0_27[k]
                  - f_5 * ms1_27[k]
                  + pa_x[k] * ns0_27[k]
                  - pc_x[k] * ns1_27[k];

        t_28[k] = f_6 * ms0_28[k]
                  - f_6 * ms1_28[k]
                  + pa_x[k] * ns0_28[k]
                  - pc_x[k] * ns1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_x, pa_z, pc_x, pc_z, ms0_30, ms1_30, ns0_21, ns0_30, \
                         ns1_21, ns1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * ns0_21[k]
                  - pc_z[k] * ns1_21[k];

        t_30[k] = f_6 * ms0_30[k]
                  - f_6 * ms1_30[k]
                  + pa_x[k] * ns0_30[k]
                  - pc_x[k] * ns1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pa_x, pc_x, ms0_31, ms0_32, ms1_31, ms1_32, ns0_31, \
                         ns0_32, ns1_31, ns1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * ms0_31[k]
                  - f_6 * ms1_31[k]
                  + pa_x[k] * ns0_31[k]
                  - pc_x[k] * ns1_31[k];

        t_32[k] = f_6 * ms0_32[k]
                  - f_6 * ms1_32[k]
                  + pa_x[k] * ns0_32[k]
                  - pc_x[k] * ns1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pa_x, pa_y, pc_x, pc_y, ms0_33, ms1_33, ns0_27, ns0_33, \
                         ns1_27, ns1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_6 * ms0_33[k]
                  - f_6 * ms1_33[k]
                  + pa_x[k] * ns0_33[k]
                  - pc_x[k] * ns1_33[k];

        t_34[k] = pa_y[k] * ns0_27[k]
                  - pc_y[k] * ns1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, pa_x, pc_x, ms0_35, ms0_36, ms1_35, ms1_36, ns0_35, \
                         ns0_36, ns1_35, ns1_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * ms0_35[k]
                  - f_6 * ms1_35[k]
                  + pa_x[k] * ns0_35[k]
                  - pc_x[k] * ns1_35[k];

        t_36[k] = f_7 * ms0_36[k]
                  - f_7 * ms1_36[k]
                  + pa_x[k] * ns0_36[k]
                  - pc_x[k] * ns1_36[k];
    }

#pragma omp simd aligned(t_37, t_38, pa_x, pa_z, pc_x, pc_z, ms0_38, ms1_38, ns0_28, ns0_38, \
                         ns1_28, ns1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_z[k] * ns0_28[k]
                  - pc_z[k] * ns1_28[k];

        t_38[k] = f_7 * ms0_38[k]
                  - f_7 * ms1_38[k]
                  + pa_x[k] * ns0_38[k]
                  - pc_x[k] * ns1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, pa_x, pc_x, ms0_39, ms0_40, ms1_39, ms1_40, ns0_39, \
                         ns0_40, ns1_39, ns1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * ms0_39[k]
                  - f_7 * ms1_39[k]
                  + pa_x[k] * ns0_39[k]
                  - pc_x[k] * ns1_39[k];

        t_40[k] = f_7 * ms0_40[k]
                  - f_7 * ms1_40[k]
                  + pa_x[k] * ns0_40[k]
                  - pc_x[k] * ns1_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pa_x, pc_x, ms0_41, ms0_42, ms1_41, ms1_42, ns0_41, \
                         ns0_42, ns1_41, ns1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * ms0_41[k]
                  - f_7 * ms1_41[k]
                  + pa_x[k] * ns0_41[k]
                  - pc_x[k] * ns1_41[k];

        t_42[k] = f_7 * ms0_42[k]
                  - f_7 * ms1_42[k]
                  + pa_x[k] * ns0_42[k]
                  - pc_x[k] * ns1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, pa_x, pa_y, pc_x, pc_y, ms0_44, ms1_44, ns0_35, ns0_44, \
                         ns1_35, ns1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_y[k] * ns0_35[k]
                  - pc_y[k] * ns1_35[k];

        t_44[k] = f_7 * ms0_44[k]
                  - f_7 * ms1_44[k]
                  + pa_x[k] * ns0_44[k]
                  - pc_x[k] * ns1_44[k];
    }

#pragma omp simd aligned(t_45, t_46, pa_x, pa_z, pc_x, pc_z, ms0_45, ms1_45, ns0_36, ns0_45, \
                         ns1_36, ns1_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * ms0_45[k]
                  - f_8 * ms1_45[k]
                  + pa_x[k] * ns0_45[k]
                  - pc_x[k] * ns1_45[k];

        t_46[k] = pa_z[k] * ns0_36[k]
                  - pc_z[k] * ns1_36[k];
    }

#pragma omp simd aligned(t_47, t_48, pa_x, pc_x, ms0_47, ms0_48, ms1_47, ms1_48, ns0_47, \
                         ns0_48, ns1_47, ns1_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_8 * ms0_47[k]
                  - f_8 * ms1_47[k]
                  + pa_x[k] * ns0_47[k]
                  - pc_x[k] * ns1_47[k];

        t_48[k] = f_8 * ms0_48[k]
                  - f_8 * ms1_48[k]
                  + pa_x[k] * ns0_48[k]
                  - pc_x[k] * ns1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, pa_x, pc_x, ms0_49, ms0_50, ms1_49, ms1_50, ns0_49, \
                         ns0_50, ns1_49, ns1_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * ms0_49[k]
                  - f_8 * ms1_49[k]
                  + pa_x[k] * ns0_49[k]
                  - pc_x[k] * ns1_49[k];

        t_50[k] = f_8 * ms0_50[k]
                  - f_8 * ms1_50[k]
                  + pa_x[k] * ns0_50[k]
                  - pc_x[k] * ns1_50[k];
    }

#pragma omp simd aligned(t_51, t_52, pa_x, pc_x, ms0_51, ms0_52, ms1_51, ms1_52, ns0_51, \
                         ns0_52, ns1_51, ns1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * ms0_51[k]
                  - f_8 * ms1_51[k]
                  + pa_x[k] * ns0_51[k]
                  - pc_x[k] * ns1_51[k];

        t_52[k] = f_8 * ms0_52[k]
                  - f_8 * ms1_52[k]
                  + pa_x[k] * ns0_52[k]
                  - pc_x[k] * ns1_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_x, pa_y, pc_x, pc_y, ms0_54, ms1_54, ns0_44, \
                         ns0_54, ns0_55, ns1_44, ns1_54, ns1_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_y[k] * ns0_44[k]
                  - pc_y[k] * ns1_44[k];

        t_54[k] = f_8 * ms0_54[k]
                  - f_8 * ms1_54[k]
                  + pa_x[k] * ns0_54[k]
                  - pc_x[k] * ns1_54[k];

        t_55[k] = pa_x[k] * ns0_55[k]
                  - pc_x[k] * ns1_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pc_x, ns0_56, ns0_57, ns0_58, ns0_59, \
                         ns1_56, ns1_57, ns1_58, ns1_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_x[k] * ns0_56[k]
                  - pc_x[k] * ns1_56[k];

        t_57[k] = pa_x[k] * ns0_57[k]
                  - pc_x[k] * ns1_57[k];

        t_58[k] = pa_x[k] * ns0_58[k]
                  - pc_x[k] * ns1_58[k];

        t_59[k] = pa_x[k] * ns0_59[k]
                  - pc_x[k] * ns1_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_x, pc_x, ns0_60, ns0_61, ns0_62, ns0_63, \
                         ns1_60, ns1_61, ns1_62, ns1_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * ns0_60[k]
                  - pc_x[k] * ns1_60[k];

        t_61[k] = pa_x[k] * ns0_61[k]
                  - pc_x[k] * ns1_61[k];

        t_62[k] = pa_x[k] * ns0_62[k]
                  - pc_x[k] * ns1_62[k];

        t_63[k] = pa_x[k] * ns0_63[k]
                  - pc_x[k] * ns1_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_x, pa_y, pc_x, pc_y, ms0_45, ms1_45, ns0_55, \
                         ns0_64, ns0_65, ns1_55, ns1_64, ns1_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pa_x[k] * ns0_64[k]
                  - pc_x[k] * ns1_64[k];

        t_65[k] = pa_x[k] * ns0_65[k]
                  - pc_x[k] * ns1_65[k];

        t_66[k] = f_0 * ms0_45[k]
                  - f_0 * ms1_45[k]
                  + pa_y[k] * ns0_55[k]
                  - pc_y[k] * ns1_55[k];
    }

#pragma omp simd aligned(t_67, t_68, pa_y, pa_z, pc_y, pc_z, ms0_47, ms1_47, ns0_55, ns0_57, \
                         ns1_55, ns1_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pa_z[k] * ns0_55[k]
                  - pc_z[k] * ns1_55[k];

        t_68[k] = f_1 * ms0_47[k]
                  - f_1 * ms1_47[k]
                  + pa_y[k] * ns0_57[k]
                  - pc_y[k] * ns1_57[k];
    }

#pragma omp simd aligned(t_69, t_70, pa_y, pc_y, ms0_48, ms0_49, ms1_48, ms1_49, ns0_58, \
                         ns0_59, ns1_58, ns1_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_2 * ms0_48[k]
                  - f_2 * ms1_48[k]
                  + pa_y[k] * ns0_58[k]
                  - pc_y[k] * ns1_58[k];

        t_70[k] = f_3 * ms0_49[k]
                  - f_3 * ms1_49[k]
                  + pa_y[k] * ns0_59[k]
                  - pc_y[k] * ns1_59[k];
    }

#pragma omp simd aligned(t_71, t_72, pa_y, pc_y, ms0_50, ms0_51, ms1_50, ms1_51, ns0_60, \
                         ns0_61, ns1_60, ns1_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_4 * ms0_50[k]
                  - f_4 * ms1_50[k]
                  + pa_y[k] * ns0_60[k]
                  - pc_y[k] * ns1_60[k];

        t_72[k] = f_5 * ms0_51[k]
                  - f_5 * ms1_51[k]
                  + pa_y[k] * ns0_61[k]
                  - pc_y[k] * ns1_61[k];
    }

#pragma omp simd aligned(t_73, t_74, pa_y, pc_y, ms0_52, ms0_53, ms1_52, ms1_53, ns0_62, \
                         ns0_63, ns1_62, ns1_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * ms0_52[k]
                  - f_6 * ms1_52[k]
                  + pa_y[k] * ns0_62[k]
                  - pc_y[k] * ns1_62[k];

        t_74[k] = f_7 * ms0_53[k]
                  - f_7 * ms1_53[k]
                  + pa_y[k] * ns0_63[k]
                  - pc_y[k] * ns1_63[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_y, pa_z, pc_y, pc_z, ms0_54, ms1_54, ns0_64, \
                         ns0_65, ns1_64, ns1_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * ms0_54[k]
                  - f_8 * ms1_54[k]
                  + pa_y[k] * ns0_64[k]
                  - pc_y[k] * ns1_64[k];

        t_76[k] = pa_y[k] * ns0_65[k]
                  - pc_y[k] * ns1_65[k];

        t_77[k] = f_0 * ms0_54[k]
                  - f_0 * ms1_54[k]
                  + pa_z[k] * ns0_65[k]
                  - pc_z[k] * ns1_65[k];
    }
}

}  // namespace simdnpot
