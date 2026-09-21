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


#include "SimdNuclearPotentialVrrRecSQ.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

static auto
compute_prim_sq_nuclear_potential_0_piece0(CSimdMatrix &buffer, const size_t target,
                                           const size_t pb, const size_t pc, const size_t sn0,
                                           const size_t so0, const size_t sn1, const size_t so1,
                                           const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / p;
    const auto f_1 = 4.5 / p;
    const auto f_2 = 4.0 / p;
    const auto f_3 = 3.5 / p;
    const auto f_4 = 3.0 / p;
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 0.5 / p;

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
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sn0_0 = buffer.data(sn0 + 0);
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
    const auto *sn0_57 = buffer.data(sn0 + 57);
    const auto *sn0_58 = buffer.data(sn0 + 58);
    const auto *sn0_59 = buffer.data(sn0 + 59);
    const auto *sn0_60 = buffer.data(sn0 + 60);
    const auto *sn0_61 = buffer.data(sn0 + 61);
    const auto *sn0_62 = buffer.data(sn0 + 62);
    const auto *sn0_63 = buffer.data(sn0 + 63);
    const auto *sn0_64 = buffer.data(sn0 + 64);
    const auto *sn0_65 = buffer.data(sn0 + 65);

    const auto *so0_0 = buffer.data(so0 + 0);
    const auto *so0_2 = buffer.data(so0 + 2);
    const auto *so0_3 = buffer.data(so0 + 3);
    const auto *so0_5 = buffer.data(so0 + 5);
    const auto *so0_6 = buffer.data(so0 + 6);
    const auto *so0_9 = buffer.data(so0 + 9);
    const auto *so0_10 = buffer.data(so0 + 10);
    const auto *so0_12 = buffer.data(so0 + 12);
    const auto *so0_14 = buffer.data(so0 + 14);
    const auto *so0_15 = buffer.data(so0 + 15);
    const auto *so0_17 = buffer.data(so0 + 17);
    const auto *so0_18 = buffer.data(so0 + 18);
    const auto *so0_20 = buffer.data(so0 + 20);
    const auto *so0_21 = buffer.data(so0 + 21);
    const auto *so0_23 = buffer.data(so0 + 23);
    const auto *so0_24 = buffer.data(so0 + 24);
    const auto *so0_25 = buffer.data(so0 + 25);
    const auto *so0_27 = buffer.data(so0 + 27);
    const auto *so0_28 = buffer.data(so0 + 28);
    const auto *so0_30 = buffer.data(so0 + 30);
    const auto *so0_31 = buffer.data(so0 + 31);
    const auto *so0_32 = buffer.data(so0 + 32);
    const auto *so0_33 = buffer.data(so0 + 33);
    const auto *so0_35 = buffer.data(so0 + 35);
    const auto *so0_36 = buffer.data(so0 + 36);
    const auto *so0_38 = buffer.data(so0 + 38);
    const auto *so0_39 = buffer.data(so0 + 39);
    const auto *so0_40 = buffer.data(so0 + 40);
    const auto *so0_41 = buffer.data(so0 + 41);
    const auto *so0_42 = buffer.data(so0 + 42);
    const auto *so0_44 = buffer.data(so0 + 44);
    const auto *so0_45 = buffer.data(so0 + 45);
    const auto *so0_47 = buffer.data(so0 + 47);
    const auto *so0_48 = buffer.data(so0 + 48);
    const auto *so0_49 = buffer.data(so0 + 49);
    const auto *so0_50 = buffer.data(so0 + 50);
    const auto *so0_51 = buffer.data(so0 + 51);
    const auto *so0_52 = buffer.data(so0 + 52);
    const auto *so0_54 = buffer.data(so0 + 54);
    const auto *so0_55 = buffer.data(so0 + 55);
    const auto *so0_57 = buffer.data(so0 + 57);
    const auto *so0_58 = buffer.data(so0 + 58);
    const auto *so0_59 = buffer.data(so0 + 59);
    const auto *so0_60 = buffer.data(so0 + 60);
    const auto *so0_61 = buffer.data(so0 + 61);
    const auto *so0_62 = buffer.data(so0 + 62);
    const auto *so0_63 = buffer.data(so0 + 63);
    const auto *so0_65 = buffer.data(so0 + 65);
    const auto *so0_66 = buffer.data(so0 + 66);
    const auto *so0_67 = buffer.data(so0 + 67);
    const auto *so0_68 = buffer.data(so0 + 68);
    const auto *so0_69 = buffer.data(so0 + 69);
    const auto *so0_70 = buffer.data(so0 + 70);
    const auto *so0_71 = buffer.data(so0 + 71);
    const auto *so0_72 = buffer.data(so0 + 72);
    const auto *so0_73 = buffer.data(so0 + 73);
    const auto *so0_74 = buffer.data(so0 + 74);
    const auto *so0_75 = buffer.data(so0 + 75);
    const auto *so0_76 = buffer.data(so0 + 76);
    const auto *so0_77 = buffer.data(so0 + 77);

    const auto *sn1_0 = buffer.data(sn1 + 0);
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
    const auto *sn1_57 = buffer.data(sn1 + 57);
    const auto *sn1_58 = buffer.data(sn1 + 58);
    const auto *sn1_59 = buffer.data(sn1 + 59);
    const auto *sn1_60 = buffer.data(sn1 + 60);
    const auto *sn1_61 = buffer.data(sn1 + 61);
    const auto *sn1_62 = buffer.data(sn1 + 62);
    const auto *sn1_63 = buffer.data(sn1 + 63);
    const auto *sn1_64 = buffer.data(sn1 + 64);
    const auto *sn1_65 = buffer.data(sn1 + 65);

    const auto *so1_0 = buffer.data(so1 + 0);
    const auto *so1_2 = buffer.data(so1 + 2);
    const auto *so1_3 = buffer.data(so1 + 3);
    const auto *so1_5 = buffer.data(so1 + 5);
    const auto *so1_6 = buffer.data(so1 + 6);
    const auto *so1_9 = buffer.data(so1 + 9);
    const auto *so1_10 = buffer.data(so1 + 10);
    const auto *so1_12 = buffer.data(so1 + 12);
    const auto *so1_14 = buffer.data(so1 + 14);
    const auto *so1_15 = buffer.data(so1 + 15);
    const auto *so1_17 = buffer.data(so1 + 17);
    const auto *so1_18 = buffer.data(so1 + 18);
    const auto *so1_20 = buffer.data(so1 + 20);
    const auto *so1_21 = buffer.data(so1 + 21);
    const auto *so1_23 = buffer.data(so1 + 23);
    const auto *so1_24 = buffer.data(so1 + 24);
    const auto *so1_25 = buffer.data(so1 + 25);
    const auto *so1_27 = buffer.data(so1 + 27);
    const auto *so1_28 = buffer.data(so1 + 28);
    const auto *so1_30 = buffer.data(so1 + 30);
    const auto *so1_31 = buffer.data(so1 + 31);
    const auto *so1_32 = buffer.data(so1 + 32);
    const auto *so1_33 = buffer.data(so1 + 33);
    const auto *so1_35 = buffer.data(so1 + 35);
    const auto *so1_36 = buffer.data(so1 + 36);
    const auto *so1_38 = buffer.data(so1 + 38);
    const auto *so1_39 = buffer.data(so1 + 39);
    const auto *so1_40 = buffer.data(so1 + 40);
    const auto *so1_41 = buffer.data(so1 + 41);
    const auto *so1_42 = buffer.data(so1 + 42);
    const auto *so1_44 = buffer.data(so1 + 44);
    const auto *so1_45 = buffer.data(so1 + 45);
    const auto *so1_47 = buffer.data(so1 + 47);
    const auto *so1_48 = buffer.data(so1 + 48);
    const auto *so1_49 = buffer.data(so1 + 49);
    const auto *so1_50 = buffer.data(so1 + 50);
    const auto *so1_51 = buffer.data(so1 + 51);
    const auto *so1_52 = buffer.data(so1 + 52);
    const auto *so1_54 = buffer.data(so1 + 54);
    const auto *so1_55 = buffer.data(so1 + 55);
    const auto *so1_57 = buffer.data(so1 + 57);
    const auto *so1_58 = buffer.data(so1 + 58);
    const auto *so1_59 = buffer.data(so1 + 59);
    const auto *so1_60 = buffer.data(so1 + 60);
    const auto *so1_61 = buffer.data(so1 + 61);
    const auto *so1_62 = buffer.data(so1 + 62);
    const auto *so1_63 = buffer.data(so1 + 63);
    const auto *so1_65 = buffer.data(so1 + 65);
    const auto *so1_66 = buffer.data(so1 + 66);
    const auto *so1_67 = buffer.data(so1 + 67);
    const auto *so1_68 = buffer.data(so1 + 68);
    const auto *so1_69 = buffer.data(so1 + 69);
    const auto *so1_70 = buffer.data(so1 + 70);
    const auto *so1_71 = buffer.data(so1 + 71);
    const auto *so1_72 = buffer.data(so1 + 72);
    const auto *so1_73 = buffer.data(so1 + 73);
    const auto *so1_74 = buffer.data(so1 + 74);
    const auto *so1_75 = buffer.data(so1 + 75);
    const auto *so1_76 = buffer.data(so1 + 76);
    const auto *so1_77 = buffer.data(so1 + 77);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sn0_0, so0_0, \
                         sn1_0, so1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sn0_0[k]
                 + pb_x[k] * so0_0[k]
                 - f_0 * sn1_0[k]
                 - pc_x[k] * so1_0[k];

        t_1[k] = pb_y[k] * so0_0[k]
                 - pc_y[k] * so1_0[k];

        t_2[k] = pb_z[k] * so0_0[k]
                 - pc_z[k] * so1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sn0_3, so0_2, so0_3, sn1_3, so1_2, \
                         so1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * sn0_3[k]
                 + pb_x[k] * so0_3[k]
                 - f_1 * sn1_3[k]
                 - pc_x[k] * so1_3[k];

        t_4[k] = pb_y[k] * so0_2[k]
                 - pc_y[k] * so1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pb_x, pc_x, sn0_5, sn0_6, so0_5, so0_6, sn1_5, sn1_6, \
                         so1_5, so1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sn0_5[k]
                 + pb_x[k] * so0_5[k]
                 - f_1 * sn1_5[k]
                 - pc_x[k] * so1_5[k];

        t_6[k] = f_2 * sn0_6[k]
                 + pb_x[k] * so0_6[k]
                 - f_2 * sn1_6[k]
                 - pc_x[k] * so1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, pc_y, pc_z, so0_3, so0_5, so1_3, \
                         so1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * so0_3[k]
                 - pc_z[k] * so1_3[k];

        t_8[k] = pb_y[k] * so0_5[k]
                 - pc_y[k] * so1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_x, pc_x, sn0_9, sn0_10, so0_9, so0_10, sn1_9, sn1_10, \
                         so1_9, so1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * sn0_9[k]
                 + pb_x[k] * so0_9[k]
                 - f_2 * sn1_9[k]
                 - pc_x[k] * so1_9[k];

        t_10[k] = f_3 * sn0_10[k]
                  + pb_x[k] * so0_10[k]
                  - f_3 * sn1_10[k]
                  - pc_x[k] * so1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pb_x, pb_z, pc_x, pc_z, sn0_12, so0_6, so0_12, sn1_12, \
                         so1_6, so1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * so0_6[k]
                  - pc_z[k] * so1_6[k];

        t_12[k] = f_3 * sn0_12[k]
                  + pb_x[k] * so0_12[k]
                  - f_3 * sn1_12[k]
                  - pc_x[k] * so1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_x, pb_y, pc_x, pc_y, sn0_14, so0_9, so0_14, sn1_14, \
                         so1_9, so1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * so0_9[k]
                  - pc_y[k] * so1_9[k];

        t_14[k] = f_3 * sn0_14[k]
                  + pb_x[k] * so0_14[k]
                  - f_3 * sn1_14[k]
                  - pc_x[k] * so1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_x, pb_z, pc_x, pc_z, sn0_15, so0_10, so0_15, sn1_15, \
                         so1_10, so1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * sn0_15[k]
                  + pb_x[k] * so0_15[k]
                  - f_4 * sn1_15[k]
                  - pc_x[k] * so1_15[k];

        t_16[k] = pb_z[k] * so0_10[k]
                  - pc_z[k] * so1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pb_x, pc_x, sn0_17, sn0_18, so0_17, so0_18, sn1_17, \
                         sn1_18, so1_17, so1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * sn0_17[k]
                  + pb_x[k] * so0_17[k]
                  - f_4 * sn1_17[k]
                  - pc_x[k] * so1_17[k];

        t_18[k] = f_4 * sn0_18[k]
                  + pb_x[k] * so0_18[k]
                  - f_4 * sn1_18[k]
                  - pc_x[k] * so1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_x, pb_y, pc_x, pc_y, sn0_20, so0_14, so0_20, sn1_20, \
                         so1_14, so1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * so0_14[k]
                  - pc_y[k] * so1_14[k];

        t_20[k] = f_4 * sn0_20[k]
                  + pb_x[k] * so0_20[k]
                  - f_4 * sn1_20[k]
                  - pc_x[k] * so1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, pb_z, pc_x, pc_z, sn0_21, so0_15, so0_21, sn1_21, \
                         so1_15, so1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * sn0_21[k]
                  + pb_x[k] * so0_21[k]
                  - f_5 * sn1_21[k]
                  - pc_x[k] * so1_21[k];

        t_22[k] = pb_z[k] * so0_15[k]
                  - pc_z[k] * so1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pb_x, pc_x, sn0_23, sn0_24, so0_23, so0_24, sn1_23, \
                         sn1_24, so1_23, so1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * sn0_23[k]
                  + pb_x[k] * so0_23[k]
                  - f_5 * sn1_23[k]
                  - pc_x[k] * so1_23[k];

        t_24[k] = f_5 * sn0_24[k]
                  + pb_x[k] * so0_24[k]
                  - f_5 * sn1_24[k]
                  - pc_x[k] * so1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_x, pb_y, pc_x, pc_y, sn0_25, so0_20, so0_25, sn1_25, \
                         so1_20, so1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * sn0_25[k]
                  + pb_x[k] * so0_25[k]
                  - f_5 * sn1_25[k]
                  - pc_x[k] * so1_25[k];

        t_26[k] = pb_y[k] * so0_20[k]
                  - pc_y[k] * so1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pb_x, pc_x, sn0_27, sn0_28, so0_27, so0_28, sn1_27, \
                         sn1_28, so1_27, so1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * sn0_27[k]
                  + pb_x[k] * so0_27[k]
                  - f_5 * sn1_27[k]
                  - pc_x[k] * so1_27[k];

        t_28[k] = f_6 * sn0_28[k]
                  + pb_x[k] * so0_28[k]
                  - f_6 * sn1_28[k]
                  - pc_x[k] * so1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pb_x, pb_z, pc_x, pc_z, sn0_30, so0_21, so0_30, sn1_30, \
                         so1_21, so1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * so0_21[k]
                  - pc_z[k] * so1_21[k];

        t_30[k] = f_6 * sn0_30[k]
                  + pb_x[k] * so0_30[k]
                  - f_6 * sn1_30[k]
                  - pc_x[k] * so1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, pc_x, sn0_31, sn0_32, so0_31, so0_32, sn1_31, \
                         sn1_32, so1_31, so1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * sn0_31[k]
                  + pb_x[k] * so0_31[k]
                  - f_6 * sn1_31[k]
                  - pc_x[k] * so1_31[k];

        t_32[k] = f_6 * sn0_32[k]
                  + pb_x[k] * so0_32[k]
                  - f_6 * sn1_32[k]
                  - pc_x[k] * so1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_y, pc_x, pc_y, sn0_33, so0_27, so0_33, sn1_33, \
                         so1_27, so1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_6 * sn0_33[k]
                  + pb_x[k] * so0_33[k]
                  - f_6 * sn1_33[k]
                  - pc_x[k] * so1_33[k];

        t_34[k] = pb_y[k] * so0_27[k]
                  - pc_y[k] * so1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, pb_x, pc_x, sn0_35, sn0_36, so0_35, so0_36, sn1_35, \
                         sn1_36, so1_35, so1_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * sn0_35[k]
                  + pb_x[k] * so0_35[k]
                  - f_6 * sn1_35[k]
                  - pc_x[k] * so1_35[k];

        t_36[k] = f_7 * sn0_36[k]
                  + pb_x[k] * so0_36[k]
                  - f_7 * sn1_36[k]
                  - pc_x[k] * so1_36[k];
    }

#pragma omp simd aligned(t_37, t_38, pb_x, pb_z, pc_x, pc_z, sn0_38, so0_28, so0_38, sn1_38, \
                         so1_28, so1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_z[k] * so0_28[k]
                  - pc_z[k] * so1_28[k];

        t_38[k] = f_7 * sn0_38[k]
                  + pb_x[k] * so0_38[k]
                  - f_7 * sn1_38[k]
                  - pc_x[k] * so1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, pb_x, pc_x, sn0_39, sn0_40, so0_39, so0_40, sn1_39, \
                         sn1_40, so1_39, so1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * sn0_39[k]
                  + pb_x[k] * so0_39[k]
                  - f_7 * sn1_39[k]
                  - pc_x[k] * so1_39[k];

        t_40[k] = f_7 * sn0_40[k]
                  + pb_x[k] * so0_40[k]
                  - f_7 * sn1_40[k]
                  - pc_x[k] * so1_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pb_x, pc_x, sn0_41, sn0_42, so0_41, so0_42, sn1_41, \
                         sn1_42, so1_41, so1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * sn0_41[k]
                  + pb_x[k] * so0_41[k]
                  - f_7 * sn1_41[k]
                  - pc_x[k] * so1_41[k];

        t_42[k] = f_7 * sn0_42[k]
                  + pb_x[k] * so0_42[k]
                  - f_7 * sn1_42[k]
                  - pc_x[k] * so1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_x, pb_y, pc_x, pc_y, sn0_44, so0_35, so0_44, sn1_44, \
                         so1_35, so1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * so0_35[k]
                  - pc_y[k] * so1_35[k];

        t_44[k] = f_7 * sn0_44[k]
                  + pb_x[k] * so0_44[k]
                  - f_7 * sn1_44[k]
                  - pc_x[k] * so1_44[k];
    }

#pragma omp simd aligned(t_45, t_46, pb_x, pb_z, pc_x, pc_z, sn0_45, so0_36, so0_45, sn1_45, \
                         so1_36, so1_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * sn0_45[k]
                  + pb_x[k] * so0_45[k]
                  - f_8 * sn1_45[k]
                  - pc_x[k] * so1_45[k];

        t_46[k] = pb_z[k] * so0_36[k]
                  - pc_z[k] * so1_36[k];
    }

#pragma omp simd aligned(t_47, t_48, pb_x, pc_x, sn0_47, sn0_48, so0_47, so0_48, sn1_47, \
                         sn1_48, so1_47, so1_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_8 * sn0_47[k]
                  + pb_x[k] * so0_47[k]
                  - f_8 * sn1_47[k]
                  - pc_x[k] * so1_47[k];

        t_48[k] = f_8 * sn0_48[k]
                  + pb_x[k] * so0_48[k]
                  - f_8 * sn1_48[k]
                  - pc_x[k] * so1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, pb_x, pc_x, sn0_49, sn0_50, so0_49, so0_50, sn1_49, \
                         sn1_50, so1_49, so1_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * sn0_49[k]
                  + pb_x[k] * so0_49[k]
                  - f_8 * sn1_49[k]
                  - pc_x[k] * so1_49[k];

        t_50[k] = f_8 * sn0_50[k]
                  + pb_x[k] * so0_50[k]
                  - f_8 * sn1_50[k]
                  - pc_x[k] * so1_50[k];
    }

#pragma omp simd aligned(t_51, t_52, pb_x, pc_x, sn0_51, sn0_52, so0_51, so0_52, sn1_51, \
                         sn1_52, so1_51, so1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * sn0_51[k]
                  + pb_x[k] * so0_51[k]
                  - f_8 * sn1_51[k]
                  - pc_x[k] * so1_51[k];

        t_52[k] = f_8 * sn0_52[k]
                  + pb_x[k] * so0_52[k]
                  - f_8 * sn1_52[k]
                  - pc_x[k] * so1_52[k];
    }

#pragma omp simd aligned(t_53, t_54, pb_x, pb_y, pc_x, pc_y, sn0_54, so0_44, so0_54, sn1_54, \
                         so1_44, so1_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * so0_44[k]
                  - pc_y[k] * so1_44[k];

        t_54[k] = f_8 * sn0_54[k]
                  + pb_x[k] * so0_54[k]
                  - f_8 * sn1_54[k]
                  - pc_x[k] * so1_54[k];
    }

#pragma omp simd aligned(t_55, t_56, pb_x, pb_z, pc_x, pc_z, sn0_55, so0_45, so0_55, sn1_55, \
                         so1_45, so1_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_9 * sn0_55[k]
                  + pb_x[k] * so0_55[k]
                  - f_9 * sn1_55[k]
                  - pc_x[k] * so1_55[k];

        t_56[k] = pb_z[k] * so0_45[k]
                  - pc_z[k] * so1_45[k];
    }

#pragma omp simd aligned(t_57, t_58, pb_x, pc_x, sn0_57, sn0_58, so0_57, so0_58, sn1_57, \
                         sn1_58, so1_57, so1_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_9 * sn0_57[k]
                  + pb_x[k] * so0_57[k]
                  - f_9 * sn1_57[k]
                  - pc_x[k] * so1_57[k];

        t_58[k] = f_9 * sn0_58[k]
                  + pb_x[k] * so0_58[k]
                  - f_9 * sn1_58[k]
                  - pc_x[k] * so1_58[k];
    }

#pragma omp simd aligned(t_59, t_60, pb_x, pc_x, sn0_59, sn0_60, so0_59, so0_60, sn1_59, \
                         sn1_60, so1_59, so1_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_9 * sn0_59[k]
                  + pb_x[k] * so0_59[k]
                  - f_9 * sn1_59[k]
                  - pc_x[k] * so1_59[k];

        t_60[k] = f_9 * sn0_60[k]
                  + pb_x[k] * so0_60[k]
                  - f_9 * sn1_60[k]
                  - pc_x[k] * so1_60[k];
    }

#pragma omp simd aligned(t_61, t_62, pb_x, pc_x, sn0_61, sn0_62, so0_61, so0_62, sn1_61, \
                         sn1_62, so1_61, so1_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_9 * sn0_61[k]
                  + pb_x[k] * so0_61[k]
                  - f_9 * sn1_61[k]
                  - pc_x[k] * so1_61[k];

        t_62[k] = f_9 * sn0_62[k]
                  + pb_x[k] * so0_62[k]
                  - f_9 * sn1_62[k]
                  - pc_x[k] * so1_62[k];
    }

#pragma omp simd aligned(t_63, t_64, pb_x, pb_y, pc_x, pc_y, sn0_63, so0_54, so0_63, sn1_63, \
                         so1_54, so1_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * sn0_63[k]
                  + pb_x[k] * so0_63[k]
                  - f_9 * sn1_63[k]
                  - pc_x[k] * so1_63[k];

        t_64[k] = pb_y[k] * so0_54[k]
                  - pc_y[k] * so1_54[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_x, pc_x, sn0_65, so0_65, so0_66, so0_67, \
                         so0_68, sn1_65, so1_65, so1_66, so1_67, \
                         so1_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_9 * sn0_65[k]
                  + pb_x[k] * so0_65[k]
                  - f_9 * sn1_65[k]
                  - pc_x[k] * so1_65[k];

        t_66[k] = pb_x[k] * so0_66[k]
                  - pc_x[k] * so1_66[k];

        t_67[k] = pb_x[k] * so0_67[k]
                  - pc_x[k] * so1_67[k];

        t_68[k] = pb_x[k] * so0_68[k]
                  - pc_x[k] * so1_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pc_x, so0_69, so0_70, so0_71, so0_72, \
                         so1_69, so1_70, so1_71, so1_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_x[k] * so0_69[k]
                  - pc_x[k] * so1_69[k];

        t_70[k] = pb_x[k] * so0_70[k]
                  - pc_x[k] * so1_70[k];

        t_71[k] = pb_x[k] * so0_71[k]
                  - pc_x[k] * so1_71[k];

        t_72[k] = pb_x[k] * so0_72[k]
                  - pc_x[k] * so1_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pb_x, pc_x, so0_73, so0_74, so0_75, so0_76, \
                         so1_73, so1_74, so1_75, so1_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_x[k] * so0_73[k]
                  - pc_x[k] * so1_73[k];

        t_74[k] = pb_x[k] * so0_74[k]
                  - pc_x[k] * so1_74[k];

        t_75[k] = pb_x[k] * so0_75[k]
                  - pc_x[k] * so1_75[k];

        t_76[k] = pb_x[k] * so0_76[k]
                  - pc_x[k] * so1_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sn0_55, so0_66, \
                         so0_77, sn1_55, so1_66, so1_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_x[k] * so0_77[k]
                  - pc_x[k] * so1_77[k];

        t_78[k] = f_0 * sn0_55[k]
                  + pb_y[k] * so0_66[k]
                  - f_0 * sn1_55[k]
                  - pc_y[k] * so1_66[k];

        t_79[k] = pb_z[k] * so0_66[k]
                  - pc_z[k] * so1_66[k];
    }

#pragma omp simd aligned(t_80, t_81, pb_y, pc_y, sn0_57, sn0_58, so0_68, so0_69, sn1_57, \
                         sn1_58, so1_68, so1_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * sn0_57[k]
                  + pb_y[k] * so0_68[k]
                  - f_1 * sn1_57[k]
                  - pc_y[k] * so1_68[k];

        t_81[k] = f_2 * sn0_58[k]
                  + pb_y[k] * so0_69[k]
                  - f_2 * sn1_58[k]
                  - pc_y[k] * so1_69[k];
    }

#pragma omp simd aligned(t_82, t_83, pb_y, pc_y, sn0_59, sn0_60, so0_70, so0_71, sn1_59, \
                         sn1_60, so1_70, so1_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * sn0_59[k]
                  + pb_y[k] * so0_70[k]
                  - f_3 * sn1_59[k]
                  - pc_y[k] * so1_70[k];

        t_83[k] = f_4 * sn0_60[k]
                  + pb_y[k] * so0_71[k]
                  - f_4 * sn1_60[k]
                  - pc_y[k] * so1_71[k];
    }

#pragma omp simd aligned(t_84, t_85, pb_y, pc_y, sn0_61, sn0_62, so0_72, so0_73, sn1_61, \
                         sn1_62, so1_72, so1_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_5 * sn0_61[k]
                  + pb_y[k] * so0_72[k]
                  - f_5 * sn1_61[k]
                  - pc_y[k] * so1_72[k];

        t_85[k] = f_6 * sn0_62[k]
                  + pb_y[k] * so0_73[k]
                  - f_6 * sn1_62[k]
                  - pc_y[k] * so1_73[k];
    }

#pragma omp simd aligned(t_86, t_87, pb_y, pc_y, sn0_63, sn0_64, so0_74, so0_75, sn1_63, \
                         sn1_64, so1_74, so1_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_7 * sn0_63[k]
                  + pb_y[k] * so0_74[k]
                  - f_7 * sn1_63[k]
                  - pc_y[k] * so1_74[k];

        t_87[k] = f_8 * sn0_64[k]
                  + pb_y[k] * so0_75[k]
                  - f_8 * sn1_64[k]
                  - pc_y[k] * so1_75[k];
    }
}

static auto
compute_prim_sq_nuclear_potential_0_piece1(CSimdMatrix &buffer, const size_t target,
                                           const size_t pb, const size_t pc, const size_t sn0,
                                           const size_t so0, const size_t sn1, const size_t so1,
                                           const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / p;
    const auto f_9 = 0.5 / p;

    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sn0_65 = buffer.data(sn0 + 65);

    const auto *so0_76 = buffer.data(so0 + 76);
    const auto *so0_77 = buffer.data(so0 + 77);

    const auto *sn1_65 = buffer.data(sn1 + 65);

    const auto *so1_76 = buffer.data(so1 + 76);
    const auto *so1_77 = buffer.data(so1 + 77);

#pragma omp simd aligned(t_88, t_89, t_90, pb_y, pb_z, pc_y, pc_z, sn0_65, so0_76, so0_77, \
                         sn1_65, so1_76, so1_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_9 * sn0_65[k]
                  + pb_y[k] * so0_76[k]
                  - f_9 * sn1_65[k]
                  - pc_y[k] * so1_76[k];

        t_89[k] = pb_y[k] * so0_77[k]
                  - pc_y[k] * so1_77[k];

        t_90[k] = f_0 * sn0_65[k]
                  + pb_z[k] * so0_77[k]
                  - f_0 * sn1_65[k]
                  - pc_z[k] * so1_77[k];
    }
}

auto
compute_prim_sq_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                    const size_t pc, const size_t sn0, const size_t so0,
                                    const size_t sn1, const size_t so1, const size_t ncols,
                                    const double p) -> void
{
    compute_prim_sq_nuclear_potential_0_piece0(buffer, target, pb, pc, sn0, so0, sn1, so1,
                                               ncols, p);

    compute_prim_sq_nuclear_potential_0_piece1(buffer, target, pb, pc, sn0, so0, sn1, so1,
                                               ncols, p);
}

}  // namespace simdnpot
