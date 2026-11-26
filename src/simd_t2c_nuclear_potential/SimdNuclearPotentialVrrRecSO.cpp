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


#include "SimdNuclearPotentialVrrRecSO.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_so_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                    const size_t pc, const size_t sm0, const size_t sn0,
                                    const size_t sm1, const size_t sn1, const size_t ncols,
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sm0_0 = buffer.data(sm0 + 0);
    const auto *sm0_3 = buffer.data(sm0 + 3);
    const auto *sm0_5 = buffer.data(sm0 + 5);
    const auto *sm0_6 = buffer.data(sm0 + 6);
    const auto *sm0_9 = buffer.data(sm0 + 9);
    const auto *sm0_10 = buffer.data(sm0 + 10);
    const auto *sm0_12 = buffer.data(sm0 + 12);
    const auto *sm0_14 = buffer.data(sm0 + 14);
    const auto *sm0_15 = buffer.data(sm0 + 15);
    const auto *sm0_17 = buffer.data(sm0 + 17);
    const auto *sm0_18 = buffer.data(sm0 + 18);
    const auto *sm0_20 = buffer.data(sm0 + 20);
    const auto *sm0_21 = buffer.data(sm0 + 21);
    const auto *sm0_23 = buffer.data(sm0 + 23);
    const auto *sm0_24 = buffer.data(sm0 + 24);
    const auto *sm0_25 = buffer.data(sm0 + 25);
    const auto *sm0_27 = buffer.data(sm0 + 27);
    const auto *sm0_28 = buffer.data(sm0 + 28);
    const auto *sm0_30 = buffer.data(sm0 + 30);
    const auto *sm0_31 = buffer.data(sm0 + 31);
    const auto *sm0_32 = buffer.data(sm0 + 32);
    const auto *sm0_33 = buffer.data(sm0 + 33);
    const auto *sm0_35 = buffer.data(sm0 + 35);
    const auto *sm0_36 = buffer.data(sm0 + 36);
    const auto *sm0_38 = buffer.data(sm0 + 38);
    const auto *sm0_39 = buffer.data(sm0 + 39);
    const auto *sm0_40 = buffer.data(sm0 + 40);
    const auto *sm0_41 = buffer.data(sm0 + 41);
    const auto *sm0_42 = buffer.data(sm0 + 42);
    const auto *sm0_44 = buffer.data(sm0 + 44);
    const auto *sm0_45 = buffer.data(sm0 + 45);
    const auto *sm0_47 = buffer.data(sm0 + 47);
    const auto *sm0_48 = buffer.data(sm0 + 48);
    const auto *sm0_49 = buffer.data(sm0 + 49);
    const auto *sm0_50 = buffer.data(sm0 + 50);
    const auto *sm0_51 = buffer.data(sm0 + 51);
    const auto *sm0_52 = buffer.data(sm0 + 52);
    const auto *sm0_53 = buffer.data(sm0 + 53);
    const auto *sm0_54 = buffer.data(sm0 + 54);

    const auto *sn0_0 = buffer.data(sn0 + 0);
    const auto *sn0_2 = buffer.data(sn0 + 2);
    const auto *sn0_3 = buffer.data(sn0 + 3);
    const auto *sn0_5 = buffer.data(sn0 + 5);
    const auto *sn0_6 = buffer.data(sn0 + 6);
    const auto *sn0_9 = buffer.data(sn0 + 9);
    const auto *sn0_10 = buffer.data(sn0 + 10);
    const auto *sn0_12 = buffer.data(sn0 + 12);
    const auto *sn0_14 = buffer.data(sn0 + 14);
    const auto *sn0_15 = buffer.data(sn0 + 15);
    const auto *sn0_17 = buffer.data(sn0 + 17);
    const auto *sn0_18 = buffer.data(sn0 + 18);
    const auto *sn0_20 = buffer.data(sn0 + 20);
    const auto *sn0_21 = buffer.data(sn0 + 21);
    const auto *sn0_23 = buffer.data(sn0 + 23);
    const auto *sn0_24 = buffer.data(sn0 + 24);
    const auto *sn0_25 = buffer.data(sn0 + 25);
    const auto *sn0_27 = buffer.data(sn0 + 27);
    const auto *sn0_28 = buffer.data(sn0 + 28);
    const auto *sn0_30 = buffer.data(sn0 + 30);
    const auto *sn0_31 = buffer.data(sn0 + 31);
    const auto *sn0_32 = buffer.data(sn0 + 32);
    const auto *sn0_33 = buffer.data(sn0 + 33);
    const auto *sn0_35 = buffer.data(sn0 + 35);
    const auto *sn0_36 = buffer.data(sn0 + 36);
    const auto *sn0_38 = buffer.data(sn0 + 38);
    const auto *sn0_39 = buffer.data(sn0 + 39);
    const auto *sn0_40 = buffer.data(sn0 + 40);
    const auto *sn0_41 = buffer.data(sn0 + 41);
    const auto *sn0_42 = buffer.data(sn0 + 42);
    const auto *sn0_44 = buffer.data(sn0 + 44);
    const auto *sn0_45 = buffer.data(sn0 + 45);
    const auto *sn0_47 = buffer.data(sn0 + 47);
    const auto *sn0_48 = buffer.data(sn0 + 48);
    const auto *sn0_49 = buffer.data(sn0 + 49);
    const auto *sn0_50 = buffer.data(sn0 + 50);
    const auto *sn0_51 = buffer.data(sn0 + 51);
    const auto *sn0_52 = buffer.data(sn0 + 52);
    const auto *sn0_54 = buffer.data(sn0 + 54);
    const auto *sn0_55 = buffer.data(sn0 + 55);
    const auto *sn0_56 = buffer.data(sn0 + 56);
    const auto *sn0_57 = buffer.data(sn0 + 57);
    const auto *sn0_58 = buffer.data(sn0 + 58);
    const auto *sn0_59 = buffer.data(sn0 + 59);
    const auto *sn0_60 = buffer.data(sn0 + 60);
    const auto *sn0_61 = buffer.data(sn0 + 61);
    const auto *sn0_62 = buffer.data(sn0 + 62);
    const auto *sn0_63 = buffer.data(sn0 + 63);
    const auto *sn0_64 = buffer.data(sn0 + 64);
    const auto *sn0_65 = buffer.data(sn0 + 65);

    const auto *sm1_0 = buffer.data(sm1 + 0);
    const auto *sm1_3 = buffer.data(sm1 + 3);
    const auto *sm1_5 = buffer.data(sm1 + 5);
    const auto *sm1_6 = buffer.data(sm1 + 6);
    const auto *sm1_9 = buffer.data(sm1 + 9);
    const auto *sm1_10 = buffer.data(sm1 + 10);
    const auto *sm1_12 = buffer.data(sm1 + 12);
    const auto *sm1_14 = buffer.data(sm1 + 14);
    const auto *sm1_15 = buffer.data(sm1 + 15);
    const auto *sm1_17 = buffer.data(sm1 + 17);
    const auto *sm1_18 = buffer.data(sm1 + 18);
    const auto *sm1_20 = buffer.data(sm1 + 20);
    const auto *sm1_21 = buffer.data(sm1 + 21);
    const auto *sm1_23 = buffer.data(sm1 + 23);
    const auto *sm1_24 = buffer.data(sm1 + 24);
    const auto *sm1_25 = buffer.data(sm1 + 25);
    const auto *sm1_27 = buffer.data(sm1 + 27);
    const auto *sm1_28 = buffer.data(sm1 + 28);
    const auto *sm1_30 = buffer.data(sm1 + 30);
    const auto *sm1_31 = buffer.data(sm1 + 31);
    const auto *sm1_32 = buffer.data(sm1 + 32);
    const auto *sm1_33 = buffer.data(sm1 + 33);
    const auto *sm1_35 = buffer.data(sm1 + 35);
    const auto *sm1_36 = buffer.data(sm1 + 36);
    const auto *sm1_38 = buffer.data(sm1 + 38);
    const auto *sm1_39 = buffer.data(sm1 + 39);
    const auto *sm1_40 = buffer.data(sm1 + 40);
    const auto *sm1_41 = buffer.data(sm1 + 41);
    const auto *sm1_42 = buffer.data(sm1 + 42);
    const auto *sm1_44 = buffer.data(sm1 + 44);
    const auto *sm1_45 = buffer.data(sm1 + 45);
    const auto *sm1_47 = buffer.data(sm1 + 47);
    const auto *sm1_48 = buffer.data(sm1 + 48);
    const auto *sm1_49 = buffer.data(sm1 + 49);
    const auto *sm1_50 = buffer.data(sm1 + 50);
    const auto *sm1_51 = buffer.data(sm1 + 51);
    const auto *sm1_52 = buffer.data(sm1 + 52);
    const auto *sm1_53 = buffer.data(sm1 + 53);
    const auto *sm1_54 = buffer.data(sm1 + 54);

    const auto *sn1_0 = buffer.data(sn1 + 0);
    const auto *sn1_2 = buffer.data(sn1 + 2);
    const auto *sn1_3 = buffer.data(sn1 + 3);
    const auto *sn1_5 = buffer.data(sn1 + 5);
    const auto *sn1_6 = buffer.data(sn1 + 6);
    const auto *sn1_9 = buffer.data(sn1 + 9);
    const auto *sn1_10 = buffer.data(sn1 + 10);
    const auto *sn1_12 = buffer.data(sn1 + 12);
    const auto *sn1_14 = buffer.data(sn1 + 14);
    const auto *sn1_15 = buffer.data(sn1 + 15);
    const auto *sn1_17 = buffer.data(sn1 + 17);
    const auto *sn1_18 = buffer.data(sn1 + 18);
    const auto *sn1_20 = buffer.data(sn1 + 20);
    const auto *sn1_21 = buffer.data(sn1 + 21);
    const auto *sn1_23 = buffer.data(sn1 + 23);
    const auto *sn1_24 = buffer.data(sn1 + 24);
    const auto *sn1_25 = buffer.data(sn1 + 25);
    const auto *sn1_27 = buffer.data(sn1 + 27);
    const auto *sn1_28 = buffer.data(sn1 + 28);
    const auto *sn1_30 = buffer.data(sn1 + 30);
    const auto *sn1_31 = buffer.data(sn1 + 31);
    const auto *sn1_32 = buffer.data(sn1 + 32);
    const auto *sn1_33 = buffer.data(sn1 + 33);
    const auto *sn1_35 = buffer.data(sn1 + 35);
    const auto *sn1_36 = buffer.data(sn1 + 36);
    const auto *sn1_38 = buffer.data(sn1 + 38);
    const auto *sn1_39 = buffer.data(sn1 + 39);
    const auto *sn1_40 = buffer.data(sn1 + 40);
    const auto *sn1_41 = buffer.data(sn1 + 41);
    const auto *sn1_42 = buffer.data(sn1 + 42);
    const auto *sn1_44 = buffer.data(sn1 + 44);
    const auto *sn1_45 = buffer.data(sn1 + 45);
    const auto *sn1_47 = buffer.data(sn1 + 47);
    const auto *sn1_48 = buffer.data(sn1 + 48);
    const auto *sn1_49 = buffer.data(sn1 + 49);
    const auto *sn1_50 = buffer.data(sn1 + 50);
    const auto *sn1_51 = buffer.data(sn1 + 51);
    const auto *sn1_52 = buffer.data(sn1 + 52);
    const auto *sn1_54 = buffer.data(sn1 + 54);
    const auto *sn1_55 = buffer.data(sn1 + 55);
    const auto *sn1_56 = buffer.data(sn1 + 56);
    const auto *sn1_57 = buffer.data(sn1 + 57);
    const auto *sn1_58 = buffer.data(sn1 + 58);
    const auto *sn1_59 = buffer.data(sn1 + 59);
    const auto *sn1_60 = buffer.data(sn1 + 60);
    const auto *sn1_61 = buffer.data(sn1 + 61);
    const auto *sn1_62 = buffer.data(sn1 + 62);
    const auto *sn1_63 = buffer.data(sn1 + 63);
    const auto *sn1_64 = buffer.data(sn1 + 64);
    const auto *sn1_65 = buffer.data(sn1 + 65);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sm0_0, sn0_0, \
                         sm1_0, sn1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sm0_0[k]
                 + pb_x[k] * sn0_0[k]
                 - f_0 * sm1_0[k]
                 - pc_x[k] * sn1_0[k];

        t_1[k] = pb_y[k] * sn0_0[k]
                 - pc_y[k] * sn1_0[k];

        t_2[k] = pb_z[k] * sn0_0[k]
                 - pc_z[k] * sn1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sm0_3, sn0_2, sn0_3, sm1_3, sn1_2, \
                         sn1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * sm0_3[k]
                 + pb_x[k] * sn0_3[k]
                 - f_1 * sm1_3[k]
                 - pc_x[k] * sn1_3[k];

        t_4[k] = pb_y[k] * sn0_2[k]
                 - pc_y[k] * sn1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pb_x, pc_x, sm0_5, sm0_6, sn0_5, sn0_6, sm1_5, sm1_6, \
                         sn1_5, sn1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sm0_5[k]
                 + pb_x[k] * sn0_5[k]
                 - f_1 * sm1_5[k]
                 - pc_x[k] * sn1_5[k];

        t_6[k] = f_2 * sm0_6[k]
                 + pb_x[k] * sn0_6[k]
                 - f_2 * sm1_6[k]
                 - pc_x[k] * sn1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, pc_y, pc_z, sn0_3, sn0_5, sn1_3, \
                         sn1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * sn0_3[k]
                 - pc_z[k] * sn1_3[k];

        t_8[k] = pb_y[k] * sn0_5[k]
                 - pc_y[k] * sn1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_x, pc_x, sm0_9, sm0_10, sn0_9, sn0_10, sm1_9, sm1_10, \
                         sn1_9, sn1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * sm0_9[k]
                 + pb_x[k] * sn0_9[k]
                 - f_2 * sm1_9[k]
                 - pc_x[k] * sn1_9[k];

        t_10[k] = f_3 * sm0_10[k]
                  + pb_x[k] * sn0_10[k]
                  - f_3 * sm1_10[k]
                  - pc_x[k] * sn1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pb_x, pb_z, pc_x, pc_z, sm0_12, sn0_6, sn0_12, sm1_12, \
                         sn1_6, sn1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * sn0_6[k]
                  - pc_z[k] * sn1_6[k];

        t_12[k] = f_3 * sm0_12[k]
                  + pb_x[k] * sn0_12[k]
                  - f_3 * sm1_12[k]
                  - pc_x[k] * sn1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_x, pb_y, pc_x, pc_y, sm0_14, sn0_9, sn0_14, sm1_14, \
                         sn1_9, sn1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * sn0_9[k]
                  - pc_y[k] * sn1_9[k];

        t_14[k] = f_3 * sm0_14[k]
                  + pb_x[k] * sn0_14[k]
                  - f_3 * sm1_14[k]
                  - pc_x[k] * sn1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_x, pb_z, pc_x, pc_z, sm0_15, sn0_10, sn0_15, sm1_15, \
                         sn1_10, sn1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * sm0_15[k]
                  + pb_x[k] * sn0_15[k]
                  - f_4 * sm1_15[k]
                  - pc_x[k] * sn1_15[k];

        t_16[k] = pb_z[k] * sn0_10[k]
                  - pc_z[k] * sn1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pb_x, pc_x, sm0_17, sm0_18, sn0_17, sn0_18, sm1_17, \
                         sm1_18, sn1_17, sn1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * sm0_17[k]
                  + pb_x[k] * sn0_17[k]
                  - f_4 * sm1_17[k]
                  - pc_x[k] * sn1_17[k];

        t_18[k] = f_4 * sm0_18[k]
                  + pb_x[k] * sn0_18[k]
                  - f_4 * sm1_18[k]
                  - pc_x[k] * sn1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_x, pb_y, pc_x, pc_y, sm0_20, sn0_14, sn0_20, sm1_20, \
                         sn1_14, sn1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * sn0_14[k]
                  - pc_y[k] * sn1_14[k];

        t_20[k] = f_4 * sm0_20[k]
                  + pb_x[k] * sn0_20[k]
                  - f_4 * sm1_20[k]
                  - pc_x[k] * sn1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, pb_z, pc_x, pc_z, sm0_21, sn0_15, sn0_21, sm1_21, \
                         sn1_15, sn1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * sm0_21[k]
                  + pb_x[k] * sn0_21[k]
                  - f_5 * sm1_21[k]
                  - pc_x[k] * sn1_21[k];

        t_22[k] = pb_z[k] * sn0_15[k]
                  - pc_z[k] * sn1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pb_x, pc_x, sm0_23, sm0_24, sn0_23, sn0_24, sm1_23, \
                         sm1_24, sn1_23, sn1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * sm0_23[k]
                  + pb_x[k] * sn0_23[k]
                  - f_5 * sm1_23[k]
                  - pc_x[k] * sn1_23[k];

        t_24[k] = f_5 * sm0_24[k]
                  + pb_x[k] * sn0_24[k]
                  - f_5 * sm1_24[k]
                  - pc_x[k] * sn1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_x, pb_y, pc_x, pc_y, sm0_25, sn0_20, sn0_25, sm1_25, \
                         sn1_20, sn1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * sm0_25[k]
                  + pb_x[k] * sn0_25[k]
                  - f_5 * sm1_25[k]
                  - pc_x[k] * sn1_25[k];

        t_26[k] = pb_y[k] * sn0_20[k]
                  - pc_y[k] * sn1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pb_x, pc_x, sm0_27, sm0_28, sn0_27, sn0_28, sm1_27, \
                         sm1_28, sn1_27, sn1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * sm0_27[k]
                  + pb_x[k] * sn0_27[k]
                  - f_5 * sm1_27[k]
                  - pc_x[k] * sn1_27[k];

        t_28[k] = f_6 * sm0_28[k]
                  + pb_x[k] * sn0_28[k]
                  - f_6 * sm1_28[k]
                  - pc_x[k] * sn1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pb_x, pb_z, pc_x, pc_z, sm0_30, sn0_21, sn0_30, sm1_30, \
                         sn1_21, sn1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * sn0_21[k]
                  - pc_z[k] * sn1_21[k];

        t_30[k] = f_6 * sm0_30[k]
                  + pb_x[k] * sn0_30[k]
                  - f_6 * sm1_30[k]
                  - pc_x[k] * sn1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, pc_x, sm0_31, sm0_32, sn0_31, sn0_32, sm1_31, \
                         sm1_32, sn1_31, sn1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * sm0_31[k]
                  + pb_x[k] * sn0_31[k]
                  - f_6 * sm1_31[k]
                  - pc_x[k] * sn1_31[k];

        t_32[k] = f_6 * sm0_32[k]
                  + pb_x[k] * sn0_32[k]
                  - f_6 * sm1_32[k]
                  - pc_x[k] * sn1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_y, pc_x, pc_y, sm0_33, sn0_27, sn0_33, sm1_33, \
                         sn1_27, sn1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_6 * sm0_33[k]
                  + pb_x[k] * sn0_33[k]
                  - f_6 * sm1_33[k]
                  - pc_x[k] * sn1_33[k];

        t_34[k] = pb_y[k] * sn0_27[k]
                  - pc_y[k] * sn1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, pb_x, pc_x, sm0_35, sm0_36, sn0_35, sn0_36, sm1_35, \
                         sm1_36, sn1_35, sn1_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * sm0_35[k]
                  + pb_x[k] * sn0_35[k]
                  - f_6 * sm1_35[k]
                  - pc_x[k] * sn1_35[k];

        t_36[k] = f_7 * sm0_36[k]
                  + pb_x[k] * sn0_36[k]
                  - f_7 * sm1_36[k]
                  - pc_x[k] * sn1_36[k];
    }

#pragma omp simd aligned(t_37, t_38, pb_x, pb_z, pc_x, pc_z, sm0_38, sn0_28, sn0_38, sm1_38, \
                         sn1_28, sn1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_z[k] * sn0_28[k]
                  - pc_z[k] * sn1_28[k];

        t_38[k] = f_7 * sm0_38[k]
                  + pb_x[k] * sn0_38[k]
                  - f_7 * sm1_38[k]
                  - pc_x[k] * sn1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, pb_x, pc_x, sm0_39, sm0_40, sn0_39, sn0_40, sm1_39, \
                         sm1_40, sn1_39, sn1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * sm0_39[k]
                  + pb_x[k] * sn0_39[k]
                  - f_7 * sm1_39[k]
                  - pc_x[k] * sn1_39[k];

        t_40[k] = f_7 * sm0_40[k]
                  + pb_x[k] * sn0_40[k]
                  - f_7 * sm1_40[k]
                  - pc_x[k] * sn1_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pb_x, pc_x, sm0_41, sm0_42, sn0_41, sn0_42, sm1_41, \
                         sm1_42, sn1_41, sn1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * sm0_41[k]
                  + pb_x[k] * sn0_41[k]
                  - f_7 * sm1_41[k]
                  - pc_x[k] * sn1_41[k];

        t_42[k] = f_7 * sm0_42[k]
                  + pb_x[k] * sn0_42[k]
                  - f_7 * sm1_42[k]
                  - pc_x[k] * sn1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_x, pb_y, pc_x, pc_y, sm0_44, sn0_35, sn0_44, sm1_44, \
                         sn1_35, sn1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * sn0_35[k]
                  - pc_y[k] * sn1_35[k];

        t_44[k] = f_7 * sm0_44[k]
                  + pb_x[k] * sn0_44[k]
                  - f_7 * sm1_44[k]
                  - pc_x[k] * sn1_44[k];
    }

#pragma omp simd aligned(t_45, t_46, pb_x, pb_z, pc_x, pc_z, sm0_45, sn0_36, sn0_45, sm1_45, \
                         sn1_36, sn1_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * sm0_45[k]
                  + pb_x[k] * sn0_45[k]
                  - f_8 * sm1_45[k]
                  - pc_x[k] * sn1_45[k];

        t_46[k] = pb_z[k] * sn0_36[k]
                  - pc_z[k] * sn1_36[k];
    }

#pragma omp simd aligned(t_47, t_48, pb_x, pc_x, sm0_47, sm0_48, sn0_47, sn0_48, sm1_47, \
                         sm1_48, sn1_47, sn1_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_8 * sm0_47[k]
                  + pb_x[k] * sn0_47[k]
                  - f_8 * sm1_47[k]
                  - pc_x[k] * sn1_47[k];

        t_48[k] = f_8 * sm0_48[k]
                  + pb_x[k] * sn0_48[k]
                  - f_8 * sm1_48[k]
                  - pc_x[k] * sn1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, pb_x, pc_x, sm0_49, sm0_50, sn0_49, sn0_50, sm1_49, \
                         sm1_50, sn1_49, sn1_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * sm0_49[k]
                  + pb_x[k] * sn0_49[k]
                  - f_8 * sm1_49[k]
                  - pc_x[k] * sn1_49[k];

        t_50[k] = f_8 * sm0_50[k]
                  + pb_x[k] * sn0_50[k]
                  - f_8 * sm1_50[k]
                  - pc_x[k] * sn1_50[k];
    }

#pragma omp simd aligned(t_51, t_52, pb_x, pc_x, sm0_51, sm0_52, sn0_51, sn0_52, sm1_51, \
                         sm1_52, sn1_51, sn1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * sm0_51[k]
                  + pb_x[k] * sn0_51[k]
                  - f_8 * sm1_51[k]
                  - pc_x[k] * sn1_51[k];

        t_52[k] = f_8 * sm0_52[k]
                  + pb_x[k] * sn0_52[k]
                  - f_8 * sm1_52[k]
                  - pc_x[k] * sn1_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, pc_x, pc_y, sm0_54, sn0_44, sn0_54, \
                         sn0_55, sm1_54, sn1_44, sn1_54, sn1_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * sn0_44[k]
                  - pc_y[k] * sn1_44[k];

        t_54[k] = f_8 * sm0_54[k]
                  + pb_x[k] * sn0_54[k]
                  - f_8 * sm1_54[k]
                  - pc_x[k] * sn1_54[k];

        t_55[k] = pb_x[k] * sn0_55[k]
                  - pc_x[k] * sn1_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pc_x, sn0_56, sn0_57, sn0_58, sn0_59, \
                         sn1_56, sn1_57, sn1_58, sn1_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * sn0_56[k]
                  - pc_x[k] * sn1_56[k];

        t_57[k] = pb_x[k] * sn0_57[k]
                  - pc_x[k] * sn1_57[k];

        t_58[k] = pb_x[k] * sn0_58[k]
                  - pc_x[k] * sn1_58[k];

        t_59[k] = pb_x[k] * sn0_59[k]
                  - pc_x[k] * sn1_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pb_x, pc_x, sn0_60, sn0_61, sn0_62, sn0_63, \
                         sn1_60, sn1_61, sn1_62, sn1_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pb_x[k] * sn0_60[k]
                  - pc_x[k] * sn1_60[k];

        t_61[k] = pb_x[k] * sn0_61[k]
                  - pc_x[k] * sn1_61[k];

        t_62[k] = pb_x[k] * sn0_62[k]
                  - pc_x[k] * sn1_62[k];

        t_63[k] = pb_x[k] * sn0_63[k]
                  - pc_x[k] * sn1_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_y, pc_x, pc_y, sm0_45, sn0_55, sn0_64, \
                         sn0_65, sm1_45, sn1_55, sn1_64, sn1_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_x[k] * sn0_64[k]
                  - pc_x[k] * sn1_64[k];

        t_65[k] = pb_x[k] * sn0_65[k]
                  - pc_x[k] * sn1_65[k];

        t_66[k] = f_0 * sm0_45[k]
                  + pb_y[k] * sn0_55[k]
                  - f_0 * sm1_45[k]
                  - pc_y[k] * sn1_55[k];
    }

#pragma omp simd aligned(t_67, t_68, pb_y, pb_z, pc_y, pc_z, sm0_47, sn0_55, sn0_57, sm1_47, \
                         sn1_55, sn1_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_z[k] * sn0_55[k]
                  - pc_z[k] * sn1_55[k];

        t_68[k] = f_1 * sm0_47[k]
                  + pb_y[k] * sn0_57[k]
                  - f_1 * sm1_47[k]
                  - pc_y[k] * sn1_57[k];
    }

#pragma omp simd aligned(t_69, t_70, pb_y, pc_y, sm0_48, sm0_49, sn0_58, sn0_59, sm1_48, \
                         sm1_49, sn1_58, sn1_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_2 * sm0_48[k]
                  + pb_y[k] * sn0_58[k]
                  - f_2 * sm1_48[k]
                  - pc_y[k] * sn1_58[k];

        t_70[k] = f_3 * sm0_49[k]
                  + pb_y[k] * sn0_59[k]
                  - f_3 * sm1_49[k]
                  - pc_y[k] * sn1_59[k];
    }

#pragma omp simd aligned(t_71, t_72, pb_y, pc_y, sm0_50, sm0_51, sn0_60, sn0_61, sm1_50, \
                         sm1_51, sn1_60, sn1_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_4 * sm0_50[k]
                  + pb_y[k] * sn0_60[k]
                  - f_4 * sm1_50[k]
                  - pc_y[k] * sn1_60[k];

        t_72[k] = f_5 * sm0_51[k]
                  + pb_y[k] * sn0_61[k]
                  - f_5 * sm1_51[k]
                  - pc_y[k] * sn1_61[k];
    }

#pragma omp simd aligned(t_73, t_74, pb_y, pc_y, sm0_52, sm0_53, sn0_62, sn0_63, sm1_52, \
                         sm1_53, sn1_62, sn1_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * sm0_52[k]
                  + pb_y[k] * sn0_62[k]
                  - f_6 * sm1_52[k]
                  - pc_y[k] * sn1_62[k];

        t_74[k] = f_7 * sm0_53[k]
                  + pb_y[k] * sn0_63[k]
                  - f_7 * sm1_53[k]
                  - pc_y[k] * sn1_63[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_y, pb_z, pc_y, pc_z, sm0_54, sn0_64, sn0_65, \
                         sm1_54, sn1_64, sn1_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * sm0_54[k]
                  + pb_y[k] * sn0_64[k]
                  - f_8 * sm1_54[k]
                  - pc_y[k] * sn1_64[k];

        t_76[k] = pb_y[k] * sn0_65[k]
                  - pc_y[k] * sn1_65[k];

        t_77[k] = f_0 * sm0_54[k]
                  + pb_z[k] * sn0_65[k]
                  - f_0 * sm1_54[k]
                  - pc_z[k] * sn1_65[k];
    }
}

}  // namespace simdnpot
