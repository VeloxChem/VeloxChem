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


#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_0 = buffer.data(sii0 + 0);
    const auto *sii0_3 = buffer.data(sii0 + 3);
    const auto *sii0_5 = buffer.data(sii0 + 5);
    const auto *sii0_6 = buffer.data(sii0 + 6);
    const auto *sii0_9 = buffer.data(sii0 + 9);
    const auto *sii0_10 = buffer.data(sii0 + 10);
    const auto *sii0_12 = buffer.data(sii0 + 12);
    const auto *sii0_14 = buffer.data(sii0 + 14);
    const auto *sii0_21 = buffer.data(sii0 + 21);
    const auto *sii0_27 = buffer.data(sii0 + 27);
    const auto *sii0_31 = buffer.data(sii0 + 31);
    const auto *sii0_34 = buffer.data(sii0 + 34);
    const auto *sii0_38 = buffer.data(sii0 + 38);
    const auto *sii0_56 = buffer.data(sii0 + 56);
    const auto *sii0_61 = buffer.data(sii0 + 61);
    const auto *sii0_65 = buffer.data(sii0 + 65);

    const auto *sih_0 = buffer.data(sih + 0);
    const auto *sih_1 = buffer.data(sih + 1);
    const auto *sih_2 = buffer.data(sih + 2);
    const auto *sih_3 = buffer.data(sih + 3);
    const auto *sih_5 = buffer.data(sih + 5);
    const auto *sih_6 = buffer.data(sih + 6);
    const auto *sih_7 = buffer.data(sih + 7);
    const auto *sih_8 = buffer.data(sih + 8);
    const auto *sih_9 = buffer.data(sih + 9);
    const auto *sih_10 = buffer.data(sih + 10);
    const auto *sih_12 = buffer.data(sih + 12);
    const auto *sih_14 = buffer.data(sih + 14);
    const auto *sih_15 = buffer.data(sih + 15);
    const auto *sih_16 = buffer.data(sih + 16);
    const auto *sih_17 = buffer.data(sih + 17);
    const auto *sih_18 = buffer.data(sih + 18);
    const auto *sih_19 = buffer.data(sih + 19);
    const auto *sih_20 = buffer.data(sih + 20);
    const auto *sih_21 = buffer.data(sih + 21);
    const auto *sih_23 = buffer.data(sih + 23);
    const auto *sih_24 = buffer.data(sih + 24);
    const auto *sih_26 = buffer.data(sih + 26);
    const auto *sih_27 = buffer.data(sih + 27);
    const auto *sih_30 = buffer.data(sih + 30);
    const auto *sih_36 = buffer.data(sih + 36);
    const auto *sih_37 = buffer.data(sih + 37);
    const auto *sih_38 = buffer.data(sih + 38);
    const auto *sih_39 = buffer.data(sih + 39);
    const auto *sih_40 = buffer.data(sih + 40);
    const auto *sih_41 = buffer.data(sih + 41);
    const auto *sih_42 = buffer.data(sih + 42);
    const auto *sih_44 = buffer.data(sih + 44);
    const auto *sih_47 = buffer.data(sih + 47);
    const auto *sih_57 = buffer.data(sih + 57);
    const auto *sih_58 = buffer.data(sih + 58);
    const auto *sih_59 = buffer.data(sih + 59);
    const auto *sih_60 = buffer.data(sih + 60);
    const auto *sih_61 = buffer.data(sih + 61);
    const auto *sih_62 = buffer.data(sih + 62);
    const auto *sih_63 = buffer.data(sih + 63);
    const auto *sih_66 = buffer.data(sih + 66);
    const auto *sih_68 = buffer.data(sih + 68);
    const auto *sih_69 = buffer.data(sih + 69);
    const auto *sih_72 = buffer.data(sih + 72);
    const auto *sih_73 = buffer.data(sih + 73);
    const auto *sih_75 = buffer.data(sih + 75);
    const auto *sih_77 = buffer.data(sih + 77);
    const auto *sih_78 = buffer.data(sih + 78);
    const auto *sih_79 = buffer.data(sih + 79);
    const auto *sih_80 = buffer.data(sih + 80);
    const auto *sih_81 = buffer.data(sih + 81);
    const auto *sih_82 = buffer.data(sih + 82);
    const auto *sih_83 = buffer.data(sih + 83);

    const auto *sii1_0 = buffer.data(sii1 + 0);
    const auto *sii1_3 = buffer.data(sii1 + 3);
    const auto *sii1_5 = buffer.data(sii1 + 5);
    const auto *sii1_6 = buffer.data(sii1 + 6);
    const auto *sii1_9 = buffer.data(sii1 + 9);
    const auto *sii1_10 = buffer.data(sii1 + 10);
    const auto *sii1_12 = buffer.data(sii1 + 12);
    const auto *sii1_14 = buffer.data(sii1 + 14);
    const auto *sii1_21 = buffer.data(sii1 + 21);
    const auto *sii1_27 = buffer.data(sii1 + 27);
    const auto *sii1_31 = buffer.data(sii1 + 31);
    const auto *sii1_34 = buffer.data(sii1 + 34);
    const auto *sii1_38 = buffer.data(sii1 + 38);
    const auto *sii1_56 = buffer.data(sii1 + 56);
    const auto *sii1_61 = buffer.data(sii1 + 61);
    const auto *sii1_65 = buffer.data(sii1 + 65);

    const auto *skg0_0 = buffer.data(skg0 + 0);
    const auto *skg0_3 = buffer.data(skg0 + 3);
    const auto *skg0_5 = buffer.data(skg0 + 5);
    const auto *skg0_6 = buffer.data(skg0 + 6);
    const auto *skg0_9 = buffer.data(skg0 + 9);
    const auto *skg0_10 = buffer.data(skg0 + 10);
    const auto *skg0_12 = buffer.data(skg0 + 12);
    const auto *skg0_13 = buffer.data(skg0 + 13);
    const auto *skg0_14 = buffer.data(skg0 + 14);
    const auto *skg0_25 = buffer.data(skg0 + 25);
    const auto *skg0_27 = buffer.data(skg0 + 27);
    const auto *skg0_28 = buffer.data(skg0 + 28);
    const auto *skg0_29 = buffer.data(skg0 + 29);
    const auto *skg0_42 = buffer.data(skg0 + 42);
    const auto *skg0_43 = buffer.data(skg0 + 43);
    const auto *skg0_44 = buffer.data(skg0 + 44);
    const auto *skg0_45 = buffer.data(skg0 + 45);
    const auto *skg0_48 = buffer.data(skg0 + 48);
    const auto *skg0_50 = buffer.data(skg0 + 50);
    const auto *skg0_51 = buffer.data(skg0 + 51);
    const auto *skg0_54 = buffer.data(skg0 + 54);
    const auto *skg0_55 = buffer.data(skg0 + 55);
    const auto *skg0_57 = buffer.data(skg0 + 57);
    const auto *skg0_58 = buffer.data(skg0 + 58);
    const auto *skg0_59 = buffer.data(skg0 + 59);

    const auto *skg1_0 = buffer.data(skg1 + 0);
    const auto *skg1_3 = buffer.data(skg1 + 3);
    const auto *skg1_5 = buffer.data(skg1 + 5);
    const auto *skg1_6 = buffer.data(skg1 + 6);
    const auto *skg1_9 = buffer.data(skg1 + 9);
    const auto *skg1_10 = buffer.data(skg1 + 10);
    const auto *skg1_12 = buffer.data(skg1 + 12);
    const auto *skg1_13 = buffer.data(skg1 + 13);
    const auto *skg1_14 = buffer.data(skg1 + 14);
    const auto *skg1_25 = buffer.data(skg1 + 25);
    const auto *skg1_27 = buffer.data(skg1 + 27);
    const auto *skg1_28 = buffer.data(skg1 + 28);
    const auto *skg1_29 = buffer.data(skg1 + 29);
    const auto *skg1_42 = buffer.data(skg1 + 42);
    const auto *skg1_43 = buffer.data(skg1 + 43);
    const auto *skg1_44 = buffer.data(skg1 + 44);
    const auto *skg1_45 = buffer.data(skg1 + 45);
    const auto *skg1_48 = buffer.data(skg1 + 48);
    const auto *skg1_50 = buffer.data(skg1 + 50);
    const auto *skg1_51 = buffer.data(skg1 + 51);
    const auto *skg1_54 = buffer.data(skg1 + 54);
    const auto *skg1_55 = buffer.data(skg1 + 55);
    const auto *skg1_57 = buffer.data(skg1 + 57);
    const auto *skg1_58 = buffer.data(skg1 + 58);
    const auto *skg1_59 = buffer.data(skg1 + 59);

    const auto *skh_0 = buffer.data(skh + 0);
    const auto *skh_2 = buffer.data(skh + 2);
    const auto *skh_3 = buffer.data(skh + 3);
    const auto *skh_5 = buffer.data(skh + 5);
    const auto *skh_6 = buffer.data(skh + 6);
    const auto *skh_9 = buffer.data(skh + 9);
    const auto *skh_10 = buffer.data(skh + 10);
    const auto *skh_12 = buffer.data(skh + 12);
    const auto *skh_14 = buffer.data(skh + 14);
    const auto *skh_15 = buffer.data(skh + 15);
    const auto *skh_16 = buffer.data(skh + 16);
    const auto *skh_17 = buffer.data(skh + 17);
    const auto *skh_18 = buffer.data(skh + 18);
    const auto *skh_19 = buffer.data(skh + 19);
    const auto *skh_20 = buffer.data(skh + 20);
    const auto *skh_21 = buffer.data(skh + 21);
    const auto *skh_23 = buffer.data(skh + 23);
    const auto *skh_24 = buffer.data(skh + 24);
    const auto *skh_26 = buffer.data(skh + 26);
    const auto *skh_27 = buffer.data(skh + 27);
    const auto *skh_30 = buffer.data(skh + 30);
    const auto *skh_36 = buffer.data(skh + 36);
    const auto *skh_37 = buffer.data(skh + 37);
    const auto *skh_38 = buffer.data(skh + 38);
    const auto *skh_39 = buffer.data(skh + 39);
    const auto *skh_40 = buffer.data(skh + 40);
    const auto *skh_41 = buffer.data(skh + 41);
    const auto *skh_42 = buffer.data(skh + 42);
    const auto *skh_44 = buffer.data(skh + 44);
    const auto *skh_45 = buffer.data(skh + 45);
    const auto *skh_47 = buffer.data(skh + 47);
    const auto *skh_48 = buffer.data(skh + 48);
    const auto *skh_51 = buffer.data(skh + 51);
    const auto *skh_57 = buffer.data(skh + 57);
    const auto *skh_58 = buffer.data(skh + 58);
    const auto *skh_59 = buffer.data(skh + 59);
    const auto *skh_60 = buffer.data(skh + 60);
    const auto *skh_61 = buffer.data(skh + 61);
    const auto *skh_62 = buffer.data(skh + 62);
    const auto *skh_63 = buffer.data(skh + 63);
    const auto *skh_65 = buffer.data(skh + 65);
    const auto *skh_66 = buffer.data(skh + 66);
    const auto *skh_68 = buffer.data(skh + 68);
    const auto *skh_69 = buffer.data(skh + 69);
    const auto *skh_72 = buffer.data(skh + 72);
    const auto *skh_73 = buffer.data(skh + 73);
    const auto *skh_75 = buffer.data(skh + 75);
    const auto *skh_77 = buffer.data(skh + 77);
    const auto *skh_78 = buffer.data(skh + 78);
    const auto *skh_79 = buffer.data(skh + 79);
    const auto *skh_80 = buffer.data(skh + 80);
    const auto *skh_81 = buffer.data(skh + 81);
    const auto *skh_82 = buffer.data(skh + 82);
    const auto *skh_83 = buffer.data(skh + 83);
    const auto *skh_84 = buffer.data(skh + 84);
    const auto *skh_86 = buffer.data(skh + 86);
    const auto *skh_87 = buffer.data(skh + 87);
    const auto *skh_89 = buffer.data(skh + 89);
    const auto *skh_90 = buffer.data(skh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sih_0, sih_3, skg0_0, skg0_3, \
                         skg1_0, skg1_3, skh_0, skh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sih_0[k]
                 + f_1 * skg0_0[k]
                 - f_2 * skg1_0[k]
                 + f_3 * pc_x[k] * skh_0[k];

        t_1[k] = f_3 * pc_y[k] * skh_0[k];

        t_2[k] = f_3 * pc_z[k] * skh_0[k];

        t_3[k] = f_0 * sih_3[k]
                 + f_4 * skg0_3[k]
                 - f_5 * skg1_3[k]
                 + f_3 * pc_x[k] * skh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sih_5, sih_6, skg0_5, skg0_6, skg1_5, \
                         skg1_6, skh_2, skh_5, skh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * skh_2[k];

        t_5[k] = f_0 * sih_5[k]
                 + f_4 * skg0_5[k]
                 - f_5 * skg1_5[k]
                 + f_3 * pc_x[k] * skh_5[k];

        t_6[k] = f_0 * sih_6[k]
                 + f_6 * skg0_6[k]
                 - f_7 * skg1_6[k]
                 + f_3 * pc_x[k] * skh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sih_9, skg0_9, skg1_9, skh_3, skh_5, \
                         skh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * skh_3[k];

        t_8[k] = f_3 * pc_y[k] * skh_5[k];

        t_9[k] = f_0 * sih_9[k]
                 + f_6 * skg0_9[k]
                 - f_7 * skg1_9[k]
                 + f_3 * pc_x[k] * skh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sih_10, sih_12, skg0_10, skg0_12, \
                         skg1_10, skg1_12, skh_6, skh_10, skh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sih_10[k]
                  + f_8 * skg0_10[k]
                  - f_9 * skg1_10[k]
                  + f_3 * pc_x[k] * skh_10[k];

        t_11[k] = f_3 * pc_z[k] * skh_6[k];

        t_12[k] = f_0 * sih_12[k]
                  + f_8 * skg0_12[k]
                  - f_9 * skg1_12[k]
                  + f_3 * pc_x[k] * skh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, sih_14, sih_15, sih_16, skg0_14, \
                         skg1_14, skh_9, skh_14, skh_15, skh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * skh_9[k];

        t_14[k] = f_0 * sih_14[k]
                  + f_8 * skg0_14[k]
                  - f_9 * skg1_14[k]
                  + f_3 * pc_x[k] * skh_14[k];

        t_15[k] = f_0 * sih_15[k]
                  + f_3 * pc_x[k] * skh_15[k];

        t_16[k] = f_0 * sih_16[k]
                  + f_3 * pc_x[k] * skh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, sih_17, sih_18, sih_19, sih_20, skh_17, \
                         skh_18, skh_19, skh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * sih_17[k]
                  + f_3 * pc_x[k] * skh_17[k];

        t_18[k] = f_0 * sih_18[k]
                  + f_3 * pc_x[k] * skh_18[k];

        t_19[k] = f_0 * sih_19[k]
                  + f_3 * pc_x[k] * skh_19[k];

        t_20[k] = f_0 * sih_20[k]
                  + f_3 * pc_x[k] * skh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, skg0_10, skg0_12, skg0_13, \
                         skg1_10, skg1_12, skg1_13, skh_15, skh_17, \
                         skh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * skg0_10[k]
                  - f_2 * skg1_10[k]
                  + f_3 * pc_y[k] * skh_15[k];

        t_22[k] = f_3 * pc_z[k] * skh_15[k];

        t_23[k] = f_4 * skg0_12[k]
                  - f_5 * skg1_12[k]
                  + f_3 * pc_y[k] * skh_17[k];

        t_24[k] = f_6 * skg0_13[k]
                  - f_7 * skg1_13[k]
                  + f_3 * pc_y[k] * skh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, sii0_0, sih_0, \
                         sii1_0, skg0_14, skg1_14, skh_19, skh_20, \
                         skh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * skg0_14[k]
                  - f_9 * skg1_14[k]
                  + f_3 * pc_y[k] * skh_19[k];

        t_26[k] = f_3 * pc_y[k] * skh_20[k];

        t_27[k] = f_1 * skg0_14[k]
                  - f_2 * skg1_14[k]
                  + f_3 * pc_z[k] * skh_20[k];

        t_28[k] = pb_y[k] * sii0_0[k]
                  - f_10 * pc_y[k] * sii1_0[k];

        t_29[k] = f_11 * sih_0[k]
                  + f_3 * pc_y[k] * skh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, sii0_3, sii0_5, sih_1, \
                         sih_2, sii1_3, sii1_5, skh_21, skh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * skh_21[k];

        t_31[k] = pb_y[k] * sii0_3[k]
                  + f_12 * sih_1[k]
                  - f_10 * pc_y[k] * sii1_3[k];

        t_32[k] = f_11 * sih_2[k]
                  + f_3 * pc_y[k] * skh_23[k];

        t_33[k] = pb_y[k] * sii0_5[k]
                  - f_10 * pc_y[k] * sii1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, sii0_6, sii0_9, sih_3, \
                         sih_5, sii1_6, sii1_9, skh_24, skh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * sii0_6[k]
                  + f_13 * sih_3[k]
                  - f_10 * pc_y[k] * sii1_6[k];

        t_35[k] = f_3 * pc_z[k] * skh_24[k];

        t_36[k] = f_11 * sih_5[k]
                  + f_3 * pc_y[k] * skh_26[k];

        t_37[k] = pb_y[k] * sii0_9[k]
                  - f_10 * pc_y[k] * sii1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sii0_10, sii0_12, sih_6, \
                         sih_8, sih_9, sii1_10, sii1_12, skh_27, \
                         skh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * sii0_10[k]
                  + f_14 * sih_6[k]
                  - f_10 * pc_y[k] * sii1_10[k];

        t_39[k] = f_3 * pc_z[k] * skh_27[k];

        t_40[k] = pb_y[k] * sii0_12[k]
                  + f_12 * sih_8[k]
                  - f_10 * pc_y[k] * sii1_12[k];

        t_41[k] = f_11 * sih_9[k]
                  + f_3 * pc_y[k] * skh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, sii0_14, sih_36, sih_37, \
                         sih_38, sii1_14, skh_36, skh_37, skh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * sii0_14[k]
                  - f_10 * pc_y[k] * sii1_14[k];

        t_43[k] = f_15 * sih_36[k]
                  + f_3 * pc_x[k] * skh_36[k];

        t_44[k] = f_15 * sih_37[k]
                  + f_3 * pc_x[k] * skh_37[k];

        t_45[k] = f_15 * sih_38[k]
                  + f_3 * pc_x[k] * skh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, sih_15, sih_39, sih_40, sih_41, \
                         skg0_25, skg1_25, skh_36, skh_39, skh_40, \
                         skh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * sih_39[k]
                  + f_3 * pc_x[k] * skh_39[k];

        t_47[k] = f_15 * sih_40[k]
                  + f_3 * pc_x[k] * skh_40[k];

        t_48[k] = f_15 * sih_41[k]
                  + f_3 * pc_x[k] * skh_41[k];

        t_49[k] = f_11 * sih_15[k]
                  + f_1 * skg0_25[k]
                  - f_2 * skg1_25[k]
                  + f_3 * pc_y[k] * skh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, sih_17, sih_18, skg0_27, skg0_28, \
                         skg1_27, skg1_28, skh_36, skh_38, skh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * skh_36[k];

        t_51[k] = f_11 * sih_17[k]
                  + f_4 * skg0_27[k]
                  - f_5 * skg1_27[k]
                  + f_3 * pc_y[k] * skh_38[k];

        t_52[k] = f_11 * sih_18[k]
                  + f_6 * skg0_28[k]
                  - f_7 * skg1_28[k]
                  + f_3 * pc_y[k] * skh_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, sii0_27, sih_19, sih_20, sii1_27, \
                         skg0_29, skg1_29, skh_40, skh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * sih_19[k]
                  + f_8 * skg0_29[k]
                  - f_9 * skg1_29[k]
                  + f_3 * pc_y[k] * skh_40[k];

        t_54[k] = f_11 * sih_20[k]
                  + f_3 * pc_y[k] * skh_41[k];

        t_55[k] = pb_y[k] * sii0_27[k]
                  - f_10 * pc_y[k] * sii1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, sii0_0, sii0_3, \
                         sih_0, sii1_0, sii1_3, skh_42, skh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * sii0_0[k]
                  - f_10 * pc_z[k] * sii1_0[k];

        t_57[k] = f_3 * pc_y[k] * skh_42[k];

        t_58[k] = f_11 * sih_0[k]
                  + f_3 * pc_z[k] * skh_42[k];

        t_59[k] = pb_z[k] * sii0_3[k]
                  - f_10 * pc_z[k] * sii1_3[k];

        t_60[k] = f_3 * pc_y[k] * skh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sii0_5, sii0_6, sih_2, \
                         sih_3, sii1_5, sii1_6, skh_45, skh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * sii0_5[k]
                  + f_12 * sih_2[k]
                  - f_10 * pc_z[k] * sii1_5[k];

        t_62[k] = pb_z[k] * sii0_6[k]
                  - f_10 * pc_z[k] * sii1_6[k];

        t_63[k] = f_11 * sih_3[k]
                  + f_3 * pc_z[k] * skh_45[k];

        t_64[k] = f_3 * pc_y[k] * skh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, sii0_9, sii0_10, sii0_12, sih_5, \
                         sih_6, sih_7, sii1_9, sii1_10, sii1_12, \
                         skh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * sii0_9[k]
                  + f_13 * sih_5[k]
                  - f_10 * pc_z[k] * sii1_9[k];

        t_66[k] = pb_z[k] * sii0_10[k]
                  - f_10 * pc_z[k] * sii1_10[k];

        t_67[k] = f_11 * sih_6[k]
                  + f_3 * pc_z[k] * skh_48[k];

        t_68[k] = pb_z[k] * sii0_12[k]
                  + f_12 * sih_7[k]
                  - f_10 * pc_z[k] * sii1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, sii0_14, sih_9, \
                         sih_57, sih_58, sii1_14, skh_51, skh_57, \
                         skh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * skh_51[k];

        t_70[k] = pb_z[k] * sii0_14[k]
                  + f_14 * sih_9[k]
                  - f_10 * pc_z[k] * sii1_14[k];

        t_71[k] = f_15 * sih_57[k]
                  + f_3 * pc_x[k] * skh_57[k];

        t_72[k] = f_15 * sih_58[k]
                  + f_3 * pc_x[k] * skh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, sih_59, sih_60, sih_61, sih_62, skh_59, \
                         skh_60, skh_61, skh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * sih_59[k]
                  + f_3 * pc_x[k] * skh_59[k];

        t_74[k] = f_15 * sih_60[k]
                  + f_3 * pc_x[k] * skh_60[k];

        t_75[k] = f_15 * sih_61[k]
                  + f_3 * pc_x[k] * skh_61[k];

        t_76[k] = f_15 * sih_62[k]
                  + f_3 * pc_x[k] * skh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, sii0_21, sih_15, sii1_21, \
                         skg0_42, skg1_42, skh_57, skh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * sii0_21[k]
                  - f_10 * pc_z[k] * sii1_21[k];

        t_78[k] = f_11 * sih_15[k]
                  + f_3 * pc_z[k] * skh_57[k];

        t_79[k] = f_4 * skg0_42[k]
                  - f_5 * skg1_42[k]
                  + f_3 * pc_y[k] * skh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, sih_20, skg0_43, skg0_44, \
                         skg1_43, skg1_44, skh_60, skh_61, skh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * skg0_43[k]
                  - f_7 * skg1_43[k]
                  + f_3 * pc_y[k] * skh_60[k];

        t_81[k] = f_8 * skg0_44[k]
                  - f_9 * skg1_44[k]
                  + f_3 * pc_y[k] * skh_61[k];

        t_82[k] = f_3 * pc_y[k] * skh_62[k];

        t_83[k] = f_11 * sih_20[k]
                  + f_1 * skg0_44[k]
                  - f_2 * skg1_44[k]
                  + f_3 * pc_z[k] * skh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, sih_21, sih_63, sih_66, \
                         skg0_45, skg0_48, skg1_45, skg1_48, skh_63, \
                         skh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_16 * sih_63[k]
                  + f_1 * skg0_45[k]
                  - f_2 * skg1_45[k]
                  + f_3 * pc_x[k] * skh_63[k];

        t_85[k] = f_12 * sih_21[k]
                  + f_3 * pc_y[k] * skh_63[k];

        t_86[k] = f_3 * pc_z[k] * skh_63[k];

        t_87[k] = f_16 * sih_66[k]
                  + f_4 * skg0_48[k]
                  - f_5 * skg1_48[k]
                  + f_3 * pc_x[k] * skh_66[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pc_x, pc_y, sih_23, sih_68, sih_69, skg0_50, \
                         skg0_51, skg1_50, skg1_51, skh_65, skh_68, \
                         skh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * sih_23[k]
                  + f_3 * pc_y[k] * skh_65[k];

        t_89[k] = f_16 * sih_68[k]
                  + f_4 * skg0_50[k]
                  - f_5 * skg1_50[k]
                  + f_3 * pc_x[k] * skh_68[k];

        t_90[k] = f_16 * sih_69[k]
                  + f_6 * skg0_51[k]
                  - f_7 * skg1_51[k]
                  + f_3 * pc_x[k] * skh_69[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pc_x, pc_y, pc_z, sih_26, sih_72, skg0_54, skg1_54, \
                         skh_66, skh_68, skh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * skh_66[k];

        t_92[k] = f_12 * sih_26[k]
                  + f_3 * pc_y[k] * skh_68[k];

        t_93[k] = f_16 * sih_72[k]
                  + f_6 * skg0_54[k]
                  - f_7 * skg1_54[k]
                  + f_3 * pc_x[k] * skh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pc_x, pc_z, sih_73, sih_75, skg0_55, skg0_57, \
                         skg1_55, skg1_57, skh_69, skh_73, skh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_16 * sih_73[k]
                  + f_8 * skg0_55[k]
                  - f_9 * skg1_55[k]
                  + f_3 * pc_x[k] * skh_73[k];

        t_95[k] = f_3 * pc_z[k] * skh_69[k];

        t_96[k] = f_16 * sih_75[k]
                  + f_8 * skg0_57[k]
                  - f_9 * skg1_57[k]
                  + f_3 * pc_x[k] * skh_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, sih_30, sih_77, sih_78, sih_79, \
                         skg0_59, skg1_59, skh_72, skh_77, skh_78, \
                         skh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * sih_30[k]
                  + f_3 * pc_y[k] * skh_72[k];

        t_98[k] = f_16 * sih_77[k]
                  + f_8 * skg0_59[k]
                  - f_9 * skg1_59[k]
                  + f_3 * pc_x[k] * skh_77[k];

        t_99[k] = f_16 * sih_78[k]
                  + f_3 * pc_x[k] * skh_78[k];

        t_100[k] = f_16 * sih_79[k]
                   + f_3 * pc_x[k] * skh_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, sih_80, sih_81, sih_82, sih_83, \
                         skh_80, skh_81, skh_82, skh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_16 * sih_80[k]
                   + f_3 * pc_x[k] * skh_80[k];

        t_102[k] = f_16 * sih_81[k]
                   + f_3 * pc_x[k] * skh_81[k];

        t_103[k] = f_16 * sih_82[k]
                   + f_3 * pc_x[k] * skh_82[k];

        t_104[k] = f_16 * sih_83[k]
                   + f_3 * pc_x[k] * skh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, sih_36, sih_38, skg0_55, skg0_57, \
                         skg1_55, skg1_57, skh_78, skh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * sih_36[k]
                   + f_1 * skg0_55[k]
                   - f_2 * skg1_55[k]
                   + f_3 * pc_y[k] * skh_78[k];

        t_106[k] = f_3 * pc_z[k] * skh_78[k];

        t_107[k] = f_12 * sih_38[k]
                   + f_4 * skg0_57[k]
                   - f_5 * skg1_57[k]
                   + f_3 * pc_y[k] * skh_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, sih_39, sih_40, sih_41, \
                         skg0_58, skg0_59, skg1_58, skg1_59, skh_81, skh_82, \
                         skh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_12 * sih_39[k]
                   + f_6 * skg0_58[k]
                   - f_7 * skg1_58[k]
                   + f_3 * pc_y[k] * skh_81[k];

        t_109[k] = f_12 * sih_40[k]
                   + f_8 * skg0_59[k]
                   - f_9 * skg1_59[k]
                   + f_3 * pc_y[k] * skh_82[k];

        t_110[k] = f_12 * sih_41[k]
                   + f_3 * pc_y[k] * skh_83[k];

        t_111[k] = f_1 * skg0_59[k]
                   - f_2 * skg1_59[k]
                   + f_3 * pc_z[k] * skh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, sii0_31, sii0_56, \
                         sih_21, sih_42, sii1_31, sii1_56, skh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * sii0_56[k]
                   - f_10 * pc_y[k] * sii1_56[k];

        t_113[k] = f_11 * sih_42[k]
                   + f_3 * pc_y[k] * skh_84[k];

        t_114[k] = f_11 * sih_21[k]
                   + f_3 * pc_z[k] * skh_84[k];

        t_115[k] = pb_z[k] * sii0_31[k]
                   - f_10 * pc_z[k] * sii1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, sii0_34, sii0_61, \
                         sih_24, sih_44, sii1_34, sii1_61, skh_86, \
                         skh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * sih_44[k]
                   + f_3 * pc_y[k] * skh_86[k];

        t_117[k] = pb_y[k] * sii0_61[k]
                   - f_10 * pc_y[k] * sii1_61[k];

        t_118[k] = pb_z[k] * sii0_34[k]
                   - f_10 * pc_z[k] * sii1_34[k];

        t_119[k] = f_11 * sih_24[k]
                   + f_3 * pc_z[k] * skh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, sii0_38, sii0_65, \
                         sih_27, sih_47, sii1_38, sii1_65, skh_89, \
                         skh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * sih_47[k]
                   + f_3 * pc_y[k] * skh_89[k];

        t_121[k] = pb_y[k] * sii0_65[k]
                   - f_10 * pc_y[k] * sii1_65[k];

        t_122[k] = pb_z[k] * sii0_38[k]
                   - f_10 * pc_z[k] * sii1_38[k];

        t_123[k] = f_11 * sih_27[k]
                   + f_3 * pc_z[k] * skh_90[k];
    }
}

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.5 / q;

    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_49 = buffer.data(sii0 + 49);
    const auto *sii0_68 = buffer.data(sii0 + 68);
    const auto *sii0_70 = buffer.data(sii0 + 70);
    const auto *sii0_83 = buffer.data(sii0 + 83);
    const auto *sii0_84 = buffer.data(sii0 + 84);
    const auto *sii0_87 = buffer.data(sii0 + 87);
    const auto *sii0_90 = buffer.data(sii0 + 90);
    const auto *sii0_94 = buffer.data(sii0 + 94);
    const auto *sii0_96 = buffer.data(sii0 + 96);
    const auto *sii0_105 = buffer.data(sii0 + 105);
    const auto *sii0_140 = buffer.data(sii0 + 140);
    const auto *sii0_143 = buffer.data(sii0 + 143);
    const auto *sii0_145 = buffer.data(sii0 + 145);
    const auto *sii0_146 = buffer.data(sii0 + 146);
    const auto *sii0_149 = buffer.data(sii0 + 149);
    const auto *sii0_150 = buffer.data(sii0 + 150);
    const auto *sii0_152 = buffer.data(sii0 + 152);
    const auto *sii0_154 = buffer.data(sii0 + 154);

    const auto *sih_36 = buffer.data(sih + 36);
    const auto *sih_42 = buffer.data(sih + 42);
    const auto *sih_45 = buffer.data(sih + 45);
    const auto *sih_48 = buffer.data(sih + 48);
    const auto *sih_50 = buffer.data(sih + 50);
    const auto *sih_51 = buffer.data(sih + 51);
    const auto *sih_57 = buffer.data(sih + 57);
    const auto *sih_59 = buffer.data(sih + 59);
    const auto *sih_60 = buffer.data(sih + 60);
    const auto *sih_61 = buffer.data(sih + 61);
    const auto *sih_62 = buffer.data(sih + 62);
    const auto *sih_63 = buffer.data(sih + 63);
    const auto *sih_65 = buffer.data(sih + 65);
    const auto *sih_66 = buffer.data(sih + 66);
    const auto *sih_68 = buffer.data(sih + 68);
    const auto *sih_69 = buffer.data(sih + 69);
    const auto *sih_70 = buffer.data(sih + 70);
    const auto *sih_72 = buffer.data(sih + 72);
    const auto *sih_78 = buffer.data(sih + 78);
    const auto *sih_80 = buffer.data(sih + 80);
    const auto *sih_81 = buffer.data(sih + 81);
    const auto *sih_82 = buffer.data(sih + 82);
    const auto *sih_83 = buffer.data(sih + 83);
    const auto *sih_84 = buffer.data(sih + 84);
    const auto *sih_86 = buffer.data(sih + 86);
    const auto *sih_87 = buffer.data(sih + 87);
    const auto *sih_89 = buffer.data(sih + 89);
    const auto *sih_90 = buffer.data(sih + 90);
    const auto *sih_93 = buffer.data(sih + 93);
    const auto *sih_99 = buffer.data(sih + 99);
    const auto *sih_100 = buffer.data(sih + 100);
    const auto *sih_101 = buffer.data(sih + 101);
    const auto *sih_102 = buffer.data(sih + 102);
    const auto *sih_103 = buffer.data(sih + 103);
    const auto *sih_104 = buffer.data(sih + 104);
    const auto *sih_105 = buffer.data(sih + 105);
    const auto *sih_106 = buffer.data(sih + 106);
    const auto *sih_107 = buffer.data(sih + 107);
    const auto *sih_108 = buffer.data(sih + 108);
    const auto *sih_110 = buffer.data(sih + 110);
    const auto *sih_111 = buffer.data(sih + 111);
    const auto *sih_113 = buffer.data(sih + 113);
    const auto *sih_114 = buffer.data(sih + 114);
    const auto *sih_115 = buffer.data(sih + 115);
    const auto *sih_117 = buffer.data(sih + 117);
    const auto *sih_119 = buffer.data(sih + 119);
    const auto *sih_120 = buffer.data(sih + 120);
    const auto *sih_121 = buffer.data(sih + 121);
    const auto *sih_122 = buffer.data(sih + 122);
    const auto *sih_123 = buffer.data(sih + 123);
    const auto *sih_124 = buffer.data(sih + 124);
    const auto *sih_125 = buffer.data(sih + 125);
    const auto *sih_126 = buffer.data(sih + 126);
    const auto *sih_129 = buffer.data(sih + 129);
    const auto *sih_131 = buffer.data(sih + 131);
    const auto *sih_132 = buffer.data(sih + 132);
    const auto *sih_135 = buffer.data(sih + 135);
    const auto *sih_136 = buffer.data(sih + 136);
    const auto *sih_138 = buffer.data(sih + 138);
    const auto *sih_140 = buffer.data(sih + 140);
    const auto *sih_141 = buffer.data(sih + 141);
    const auto *sih_142 = buffer.data(sih + 142);
    const auto *sih_143 = buffer.data(sih + 143);
    const auto *sih_144 = buffer.data(sih + 144);
    const auto *sih_145 = buffer.data(sih + 145);
    const auto *sih_146 = buffer.data(sih + 146);
    const auto *sih_152 = buffer.data(sih + 152);
    const auto *sih_156 = buffer.data(sih + 156);
    const auto *sih_161 = buffer.data(sih + 161);
    const auto *sih_162 = buffer.data(sih + 162);
    const auto *sih_163 = buffer.data(sih + 163);
    const auto *sih_164 = buffer.data(sih + 164);
    const auto *sih_165 = buffer.data(sih + 165);
    const auto *sih_166 = buffer.data(sih + 166);
    const auto *sih_167 = buffer.data(sih + 167);
    const auto *sih_183 = buffer.data(sih + 183);

    const auto *sii1_49 = buffer.data(sii1 + 49);
    const auto *sii1_68 = buffer.data(sii1 + 68);
    const auto *sii1_70 = buffer.data(sii1 + 70);
    const auto *sii1_83 = buffer.data(sii1 + 83);
    const auto *sii1_84 = buffer.data(sii1 + 84);
    const auto *sii1_87 = buffer.data(sii1 + 87);
    const auto *sii1_90 = buffer.data(sii1 + 90);
    const auto *sii1_94 = buffer.data(sii1 + 94);
    const auto *sii1_96 = buffer.data(sii1 + 96);
    const auto *sii1_105 = buffer.data(sii1 + 105);
    const auto *sii1_140 = buffer.data(sii1 + 140);
    const auto *sii1_143 = buffer.data(sii1 + 143);
    const auto *sii1_145 = buffer.data(sii1 + 145);
    const auto *sii1_146 = buffer.data(sii1 + 146);
    const auto *sii1_149 = buffer.data(sii1 + 149);
    const auto *sii1_150 = buffer.data(sii1 + 150);
    const auto *sii1_152 = buffer.data(sii1 + 152);
    const auto *sii1_154 = buffer.data(sii1 + 154);

    const auto *skg0_72 = buffer.data(skg0 + 72);
    const auto *skg0_73 = buffer.data(skg0 + 73);
    const auto *skg0_74 = buffer.data(skg0 + 74);
    const auto *skg0_75 = buffer.data(skg0 + 75);
    const auto *skg0_78 = buffer.data(skg0 + 78);
    const auto *skg0_80 = buffer.data(skg0 + 80);
    const auto *skg0_81 = buffer.data(skg0 + 81);
    const auto *skg0_84 = buffer.data(skg0 + 84);
    const auto *skg0_85 = buffer.data(skg0 + 85);
    const auto *skg0_87 = buffer.data(skg0 + 87);
    const auto *skg0_88 = buffer.data(skg0 + 88);
    const auto *skg0_89 = buffer.data(skg0 + 89);
    const auto *skg0_90 = buffer.data(skg0 + 90);
    const auto *skg0_93 = buffer.data(skg0 + 93);
    const auto *skg0_95 = buffer.data(skg0 + 95);
    const auto *skg0_96 = buffer.data(skg0 + 96);
    const auto *skg0_99 = buffer.data(skg0 + 99);
    const auto *skg0_100 = buffer.data(skg0 + 100);
    const auto *skg0_102 = buffer.data(skg0 + 102);
    const auto *skg0_103 = buffer.data(skg0 + 103);
    const auto *skg0_104 = buffer.data(skg0 + 104);
    const auto *skg0_110 = buffer.data(skg0 + 110);
    const auto *skg0_114 = buffer.data(skg0 + 114);
    const auto *skg0_117 = buffer.data(skg0 + 117);
    const auto *skg0_118 = buffer.data(skg0 + 118);
    const auto *skg0_119 = buffer.data(skg0 + 119);

    const auto *skg1_72 = buffer.data(skg1 + 72);
    const auto *skg1_73 = buffer.data(skg1 + 73);
    const auto *skg1_74 = buffer.data(skg1 + 74);
    const auto *skg1_75 = buffer.data(skg1 + 75);
    const auto *skg1_78 = buffer.data(skg1 + 78);
    const auto *skg1_80 = buffer.data(skg1 + 80);
    const auto *skg1_81 = buffer.data(skg1 + 81);
    const auto *skg1_84 = buffer.data(skg1 + 84);
    const auto *skg1_85 = buffer.data(skg1 + 85);
    const auto *skg1_87 = buffer.data(skg1 + 87);
    const auto *skg1_88 = buffer.data(skg1 + 88);
    const auto *skg1_89 = buffer.data(skg1 + 89);
    const auto *skg1_90 = buffer.data(skg1 + 90);
    const auto *skg1_93 = buffer.data(skg1 + 93);
    const auto *skg1_95 = buffer.data(skg1 + 95);
    const auto *skg1_96 = buffer.data(skg1 + 96);
    const auto *skg1_99 = buffer.data(skg1 + 99);
    const auto *skg1_100 = buffer.data(skg1 + 100);
    const auto *skg1_102 = buffer.data(skg1 + 102);
    const auto *skg1_103 = buffer.data(skg1 + 103);
    const auto *skg1_104 = buffer.data(skg1 + 104);
    const auto *skg1_110 = buffer.data(skg1 + 110);
    const auto *skg1_114 = buffer.data(skg1 + 114);
    const auto *skg1_117 = buffer.data(skg1 + 117);
    const auto *skg1_118 = buffer.data(skg1 + 118);
    const auto *skg1_119 = buffer.data(skg1 + 119);

    const auto *skh_93 = buffer.data(skh + 93);
    const auto *skh_99 = buffer.data(skh + 99);
    const auto *skh_100 = buffer.data(skh + 100);
    const auto *skh_101 = buffer.data(skh + 101);
    const auto *skh_102 = buffer.data(skh + 102);
    const auto *skh_103 = buffer.data(skh + 103);
    const auto *skh_104 = buffer.data(skh + 104);
    const auto *skh_105 = buffer.data(skh + 105);
    const auto *skh_107 = buffer.data(skh + 107);
    const auto *skh_108 = buffer.data(skh + 108);
    const auto *skh_110 = buffer.data(skh + 110);
    const auto *skh_111 = buffer.data(skh + 111);
    const auto *skh_114 = buffer.data(skh + 114);
    const auto *skh_115 = buffer.data(skh + 115);
    const auto *skh_117 = buffer.data(skh + 117);
    const auto *skh_119 = buffer.data(skh + 119);
    const auto *skh_120 = buffer.data(skh + 120);
    const auto *skh_121 = buffer.data(skh + 121);
    const auto *skh_122 = buffer.data(skh + 122);
    const auto *skh_123 = buffer.data(skh + 123);
    const auto *skh_124 = buffer.data(skh + 124);
    const auto *skh_125 = buffer.data(skh + 125);
    const auto *skh_126 = buffer.data(skh + 126);
    const auto *skh_128 = buffer.data(skh + 128);
    const auto *skh_129 = buffer.data(skh + 129);
    const auto *skh_131 = buffer.data(skh + 131);
    const auto *skh_132 = buffer.data(skh + 132);
    const auto *skh_135 = buffer.data(skh + 135);
    const auto *skh_136 = buffer.data(skh + 136);
    const auto *skh_138 = buffer.data(skh + 138);
    const auto *skh_140 = buffer.data(skh + 140);
    const auto *skh_141 = buffer.data(skh + 141);
    const auto *skh_142 = buffer.data(skh + 142);
    const auto *skh_143 = buffer.data(skh + 143);
    const auto *skh_144 = buffer.data(skh + 144);
    const auto *skh_145 = buffer.data(skh + 145);
    const auto *skh_146 = buffer.data(skh + 146);
    const auto *skh_147 = buffer.data(skh + 147);
    const auto *skh_149 = buffer.data(skh + 149);
    const auto *skh_150 = buffer.data(skh + 150);
    const auto *skh_152 = buffer.data(skh + 152);
    const auto *skh_153 = buffer.data(skh + 153);
    const auto *skh_156 = buffer.data(skh + 156);
    const auto *skh_161 = buffer.data(skh + 161);
    const auto *skh_162 = buffer.data(skh + 162);
    const auto *skh_163 = buffer.data(skh + 163);
    const auto *skh_164 = buffer.data(skh + 164);
    const auto *skh_165 = buffer.data(skh + 165);
    const auto *skh_166 = buffer.data(skh + 166);
    const auto *skh_167 = buffer.data(skh + 167);
    const auto *skh_168 = buffer.data(skh + 168);
    const auto *skh_170 = buffer.data(skh + 170);
    const auto *skh_171 = buffer.data(skh + 171);
    const auto *skh_173 = buffer.data(skh + 173);
    const auto *skh_174 = buffer.data(skh + 174);
    const auto *skh_177 = buffer.data(skh + 177);
    const auto *skh_183 = buffer.data(skh + 183);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, sii0_68, sii0_70, \
                         sih_50, sih_51, sih_99, sii1_68, sii1_70, skh_93, \
                         skh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * sii0_68[k]
                   + f_12 * sih_50[k]
                   - f_10 * pc_y[k] * sii1_68[k];

        t_125[k] = f_11 * sih_51[k]
                   + f_3 * pc_y[k] * skh_93[k];

        t_126[k] = pb_y[k] * sii0_70[k]
                   - f_10 * pc_y[k] * sii1_70[k];

        t_127[k] = f_16 * sih_99[k]
                   + f_3 * pc_x[k] * skh_99[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, sih_100, sih_101, sih_102, \
                         sih_103, sih_104, skh_100, skh_101, skh_102, skh_103, \
                         skh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_16 * sih_100[k]
                   + f_3 * pc_x[k] * skh_100[k];

        t_129[k] = f_16 * sih_101[k]
                   + f_3 * pc_x[k] * skh_101[k];

        t_130[k] = f_16 * sih_102[k]
                   + f_3 * pc_x[k] * skh_102[k];

        t_131[k] = f_16 * sih_103[k]
                   + f_3 * pc_x[k] * skh_103[k];

        t_132[k] = f_16 * sih_104[k]
                   + f_3 * pc_x[k] * skh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, pc_y, pc_z, sii0_49, sih_36, sih_59, \
                         sii1_49, skg0_72, skg1_72, skh_99, skh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_z[k] * sii0_49[k]
                   - f_10 * pc_z[k] * sii1_49[k];

        t_134[k] = f_11 * sih_36[k]
                   + f_3 * pc_z[k] * skh_99[k];

        t_135[k] = f_11 * sih_59[k]
                   + f_4 * skg0_72[k]
                   - f_5 * skg1_72[k]
                   + f_3 * pc_y[k] * skh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, sih_60, sih_61, sih_62, skg0_73, skg0_74, \
                         skg1_73, skg1_74, skh_102, skh_103, skh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * sih_60[k]
                   + f_6 * skg0_73[k]
                   - f_7 * skg1_73[k]
                   + f_3 * pc_y[k] * skh_102[k];

        t_137[k] = f_11 * sih_61[k]
                   + f_8 * skg0_74[k]
                   - f_9 * skg1_74[k]
                   + f_3 * pc_y[k] * skh_103[k];

        t_138[k] = f_11 * sih_62[k]
                   + f_3 * pc_y[k] * skh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_y, pc_x, pc_y, pc_z, sii0_83, sih_42, \
                         sih_105, sii1_83, skg0_75, skg1_75, skh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * sii0_83[k]
                   - f_10 * pc_y[k] * sii1_83[k];

        t_140[k] = f_16 * sih_105[k]
                   + f_1 * skg0_75[k]
                   - f_2 * skg1_75[k]
                   + f_3 * pc_x[k] * skh_105[k];

        t_141[k] = f_3 * pc_y[k] * skh_105[k];

        t_142[k] = f_12 * sih_42[k]
                   + f_3 * pc_z[k] * skh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, sih_108, sih_110, skg0_78, skg0_80, \
                         skg1_78, skg1_80, skh_107, skh_108, skh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * sih_108[k]
                   + f_4 * skg0_78[k]
                   - f_5 * skg1_78[k]
                   + f_3 * pc_x[k] * skh_108[k];

        t_144[k] = f_3 * pc_y[k] * skh_107[k];

        t_145[k] = f_16 * sih_110[k]
                   + f_4 * skg0_80[k]
                   - f_5 * skg1_80[k]
                   + f_3 * pc_x[k] * skh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_x, pc_y, pc_z, sih_45, sih_111, skg0_81, \
                         skg1_81, skh_108, skh_110, skh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_16 * sih_111[k]
                   + f_6 * skg0_81[k]
                   - f_7 * skg1_81[k]
                   + f_3 * pc_x[k] * skh_111[k];

        t_147[k] = f_12 * sih_45[k]
                   + f_3 * pc_z[k] * skh_108[k];

        t_148[k] = f_3 * pc_y[k] * skh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, sih_48, sih_114, sih_115, skg0_84, \
                         skg0_85, skg1_84, skg1_85, skh_111, skh_114, \
                         skh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_16 * sih_114[k]
                   + f_6 * skg0_84[k]
                   - f_7 * skg1_84[k]
                   + f_3 * pc_x[k] * skh_114[k];

        t_150[k] = f_16 * sih_115[k]
                   + f_8 * skg0_85[k]
                   - f_9 * skg1_85[k]
                   + f_3 * pc_x[k] * skh_115[k];

        t_151[k] = f_12 * sih_48[k]
                   + f_3 * pc_z[k] * skh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pc_x, pc_y, sih_117, sih_119, skg0_87, skg0_89, \
                         skg1_87, skg1_89, skh_114, skh_117, skh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_16 * sih_117[k]
                   + f_8 * skg0_87[k]
                   - f_9 * skg1_87[k]
                   + f_3 * pc_x[k] * skh_117[k];

        t_153[k] = f_3 * pc_y[k] * skh_114[k];

        t_154[k] = f_16 * sih_119[k]
                   + f_8 * skg0_89[k]
                   - f_9 * skg1_89[k]
                   + f_3 * pc_x[k] * skh_119[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, sih_120, sih_121, sih_122, \
                         sih_123, sih_124, skh_120, skh_121, skh_122, skh_123, \
                         skh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_16 * sih_120[k]
                   + f_3 * pc_x[k] * skh_120[k];

        t_156[k] = f_16 * sih_121[k]
                   + f_3 * pc_x[k] * skh_121[k];

        t_157[k] = f_16 * sih_122[k]
                   + f_3 * pc_x[k] * skh_122[k];

        t_158[k] = f_16 * sih_123[k]
                   + f_3 * pc_x[k] * skh_123[k];

        t_159[k] = f_16 * sih_124[k]
                   + f_3 * pc_x[k] * skh_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, sih_57, sih_125, \
                         skg0_85, skg0_87, skg1_85, skg1_87, skh_120, skh_122, \
                         skh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * sih_125[k]
                   + f_3 * pc_x[k] * skh_125[k];

        t_161[k] = f_1 * skg0_85[k]
                   - f_2 * skg1_85[k]
                   + f_3 * pc_y[k] * skh_120[k];

        t_162[k] = f_12 * sih_57[k]
                   + f_3 * pc_z[k] * skh_120[k];

        t_163[k] = f_4 * skg0_87[k]
                   - f_5 * skg1_87[k]
                   + f_3 * pc_y[k] * skh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, sih_62, skg0_88, skg0_89, \
                         skg1_88, skg1_89, skh_123, skh_124, skh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * skg0_88[k]
                   - f_7 * skg1_88[k]
                   + f_3 * pc_y[k] * skh_123[k];

        t_165[k] = f_8 * skg0_89[k]
                   - f_9 * skg1_89[k]
                   + f_3 * pc_y[k] * skh_124[k];

        t_166[k] = f_3 * pc_y[k] * skh_125[k];

        t_167[k] = f_12 * sih_62[k]
                   + f_1 * skg0_89[k]
                   - f_2 * skg1_89[k]
                   + f_3 * pc_z[k] * skh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, sih_63, sih_126, \
                         sih_129, skg0_90, skg0_93, skg1_90, skg1_93, skh_126, \
                         skh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_14 * sih_126[k]
                   + f_1 * skg0_90[k]
                   - f_2 * skg1_90[k]
                   + f_3 * pc_x[k] * skh_126[k];

        t_169[k] = f_13 * sih_63[k]
                   + f_3 * pc_y[k] * skh_126[k];

        t_170[k] = f_3 * pc_z[k] * skh_126[k];

        t_171[k] = f_14 * sih_129[k]
                   + f_4 * skg0_93[k]
                   - f_5 * skg1_93[k]
                   + f_3 * pc_x[k] * skh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pc_x, pc_y, sih_65, sih_131, sih_132, skg0_95, \
                         skg0_96, skg1_95, skg1_96, skh_128, skh_131, \
                         skh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * sih_65[k]
                   + f_3 * pc_y[k] * skh_128[k];

        t_173[k] = f_14 * sih_131[k]
                   + f_4 * skg0_95[k]
                   - f_5 * skg1_95[k]
                   + f_3 * pc_x[k] * skh_131[k];

        t_174[k] = f_14 * sih_132[k]
                   + f_6 * skg0_96[k]
                   - f_7 * skg1_96[k]
                   + f_3 * pc_x[k] * skh_132[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pc_x, pc_y, pc_z, sih_68, sih_135, skg0_99, \
                         skg1_99, skh_129, skh_131, skh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * pc_z[k] * skh_129[k];

        t_176[k] = f_13 * sih_68[k]
                   + f_3 * pc_y[k] * skh_131[k];

        t_177[k] = f_14 * sih_135[k]
                   + f_6 * skg0_99[k]
                   - f_7 * skg1_99[k]
                   + f_3 * pc_x[k] * skh_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_z, sih_136, sih_138, skg0_100, \
                         skg0_102, skg1_100, skg1_102, skh_132, skh_136, \
                         skh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * sih_136[k]
                   + f_8 * skg0_100[k]
                   - f_9 * skg1_100[k]
                   + f_3 * pc_x[k] * skh_136[k];

        t_179[k] = f_3 * pc_z[k] * skh_132[k];

        t_180[k] = f_14 * sih_138[k]
                   + f_8 * skg0_102[k]
                   - f_9 * skg1_102[k]
                   + f_3 * pc_x[k] * skh_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, sih_72, sih_140, sih_141, \
                         sih_142, skg0_104, skg1_104, skh_135, skh_140, skh_141, \
                         skh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * sih_72[k]
                   + f_3 * pc_y[k] * skh_135[k];

        t_182[k] = f_14 * sih_140[k]
                   + f_8 * skg0_104[k]
                   - f_9 * skg1_104[k]
                   + f_3 * pc_x[k] * skh_140[k];

        t_183[k] = f_14 * sih_141[k]
                   + f_3 * pc_x[k] * skh_141[k];

        t_184[k] = f_14 * sih_142[k]
                   + f_3 * pc_x[k] * skh_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, sih_143, sih_144, sih_145, sih_146, \
                         skh_143, skh_144, skh_145, skh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_14 * sih_143[k]
                   + f_3 * pc_x[k] * skh_143[k];

        t_186[k] = f_14 * sih_144[k]
                   + f_3 * pc_x[k] * skh_144[k];

        t_187[k] = f_14 * sih_145[k]
                   + f_3 * pc_x[k] * skh_145[k];

        t_188[k] = f_14 * sih_146[k]
                   + f_3 * pc_x[k] * skh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_y, pc_z, sih_78, sih_80, skg0_100, skg0_102, \
                         skg1_100, skg1_102, skh_141, skh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * sih_78[k]
                   + f_1 * skg0_100[k]
                   - f_2 * skg1_100[k]
                   + f_3 * pc_y[k] * skh_141[k];

        t_190[k] = f_3 * pc_z[k] * skh_141[k];

        t_191[k] = f_13 * sih_80[k]
                   + f_4 * skg0_102[k]
                   - f_5 * skg1_102[k]
                   + f_3 * pc_y[k] * skh_143[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, sih_81, sih_82, sih_83, \
                         skg0_103, skg0_104, skg1_103, skg1_104, skh_144, skh_145, \
                         skh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_13 * sih_81[k]
                   + f_6 * skg0_103[k]
                   - f_7 * skg1_103[k]
                   + f_3 * pc_y[k] * skh_144[k];

        t_193[k] = f_13 * sih_82[k]
                   + f_8 * skg0_104[k]
                   - f_9 * skg1_104[k]
                   + f_3 * pc_y[k] * skh_145[k];

        t_194[k] = f_13 * sih_83[k]
                   + f_3 * pc_y[k] * skh_146[k];

        t_195[k] = f_1 * skg0_104[k]
                   - f_2 * skg1_104[k]
                   + f_3 * pc_z[k] * skh_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_z, pc_y, pc_z, sii0_84, sii0_87, \
                         sih_63, sih_84, sii1_84, sii1_87, skh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * sii0_84[k]
                   - f_10 * pc_z[k] * sii1_84[k];

        t_197[k] = f_12 * sih_84[k]
                   + f_3 * pc_y[k] * skh_147[k];

        t_198[k] = f_11 * sih_63[k]
                   + f_3 * pc_z[k] * skh_147[k];

        t_199[k] = pb_z[k] * sii0_87[k]
                   - f_10 * pc_z[k] * sii1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_z, pc_x, pc_y, pc_z, sii0_90, sih_86, \
                         sih_152, sii1_90, skg0_110, skg1_110, skh_149, \
                         skh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * sih_86[k]
                   + f_3 * pc_y[k] * skh_149[k];

        t_201[k] = f_14 * sih_152[k]
                   + f_4 * skg0_110[k]
                   - f_5 * skg1_110[k]
                   + f_3 * pc_x[k] * skh_152[k];

        t_202[k] = pb_z[k] * sii0_90[k]
                   - f_10 * pc_z[k] * sii1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, pc_y, pc_z, sih_66, sih_89, sih_156, \
                         skg0_114, skg1_114, skh_150, skh_152, \
                         skh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * sih_66[k]
                   + f_3 * pc_z[k] * skh_150[k];

        t_204[k] = f_12 * sih_89[k]
                   + f_3 * pc_y[k] * skh_152[k];

        t_205[k] = f_14 * sih_156[k]
                   + f_6 * skg0_114[k]
                   - f_7 * skg1_114[k]
                   + f_3 * pc_x[k] * skh_156[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_z, pc_y, pc_z, sii0_94, sii0_96, \
                         sih_69, sih_70, sih_93, sii1_94, sii1_96, skh_153, \
                         skh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * sii0_94[k]
                   - f_10 * pc_z[k] * sii1_94[k];

        t_207[k] = f_11 * sih_69[k]
                   + f_3 * pc_z[k] * skh_153[k];

        t_208[k] = pb_z[k] * sii0_96[k]
                   + f_12 * sih_70[k]
                   - f_10 * pc_z[k] * sii1_96[k];

        t_209[k] = f_12 * sih_93[k]
                   + f_3 * pc_y[k] * skh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, sih_161, sih_162, sih_163, sih_164, \
                         skg0_119, skg1_119, skh_161, skh_162, skh_163, \
                         skh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_14 * sih_161[k]
                   + f_8 * skg0_119[k]
                   - f_9 * skg1_119[k]
                   + f_3 * pc_x[k] * skh_161[k];

        t_211[k] = f_14 * sih_162[k]
                   + f_3 * pc_x[k] * skh_162[k];

        t_212[k] = f_14 * sih_163[k]
                   + f_3 * pc_x[k] * skh_163[k];

        t_213[k] = f_14 * sih_164[k]
                   + f_3 * pc_x[k] * skh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, sii0_105, sih_165, \
                         sih_166, sih_167, sii1_105, skh_165, skh_166, \
                         skh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_14 * sih_165[k]
                   + f_3 * pc_x[k] * skh_165[k];

        t_215[k] = f_14 * sih_166[k]
                   + f_3 * pc_x[k] * skh_166[k];

        t_216[k] = f_14 * sih_167[k]
                   + f_3 * pc_x[k] * skh_167[k];

        t_217[k] = pb_z[k] * sii0_105[k]
                   - f_10 * pc_z[k] * sii1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, sih_78, sih_101, sih_102, skg0_117, \
                         skg0_118, skg1_117, skg1_118, skh_162, skh_164, \
                         skh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * sih_78[k]
                   + f_3 * pc_z[k] * skh_162[k];

        t_219[k] = f_12 * sih_101[k]
                   + f_4 * skg0_117[k]
                   - f_5 * skg1_117[k]
                   + f_3 * pc_y[k] * skh_164[k];

        t_220[k] = f_12 * sih_102[k]
                   + f_6 * skg0_118[k]
                   - f_7 * skg1_118[k]
                   + f_3 * pc_y[k] * skh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_y, pc_y, pc_z, sii0_140, sih_83, \
                         sih_103, sih_104, sii1_140, skg0_119, skg1_119, skh_166, \
                         skh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * sih_103[k]
                   + f_8 * skg0_119[k]
                   - f_9 * skg1_119[k]
                   + f_3 * pc_y[k] * skh_166[k];

        t_222[k] = f_12 * sih_104[k]
                   + f_3 * pc_y[k] * skh_167[k];

        t_223[k] = f_11 * sih_83[k]
                   + f_1 * skg0_119[k]
                   - f_2 * skg1_119[k]
                   + f_3 * pc_z[k] * skh_167[k];

        t_224[k] = pb_y[k] * sii0_140[k]
                   - f_10 * pc_y[k] * sii1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_y, pc_y, pc_z, sii0_143, sih_84, \
                         sih_105, sih_106, sih_107, sii1_143, skh_168, \
                         skh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * sih_105[k]
                   + f_3 * pc_y[k] * skh_168[k];

        t_226[k] = f_12 * sih_84[k]
                   + f_3 * pc_z[k] * skh_168[k];

        t_227[k] = pb_y[k] * sii0_143[k]
                   + f_12 * sih_106[k]
                   - f_10 * pc_y[k] * sii1_143[k];

        t_228[k] = f_11 * sih_107[k]
                   + f_3 * pc_y[k] * skh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_y, pc_y, pc_z, sii0_145, sii0_146, \
                         sih_87, sih_108, sih_110, sii1_145, sii1_146, skh_171, \
                         skh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * sii0_145[k]
                   - f_10 * pc_y[k] * sii1_145[k];

        t_230[k] = pb_y[k] * sii0_146[k]
                   + f_13 * sih_108[k]
                   - f_10 * pc_y[k] * sii1_146[k];

        t_231[k] = f_12 * sih_87[k]
                   + f_3 * pc_z[k] * skh_171[k];

        t_232[k] = f_11 * sih_110[k]
                   + f_3 * pc_y[k] * skh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pc_y, pc_z, sii0_149, sii0_150, sih_90, \
                         sih_111, sii1_149, sii1_150, skh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * sii0_149[k]
                   - f_10 * pc_y[k] * sii1_149[k];

        t_234[k] = pb_y[k] * sii0_150[k]
                   + f_14 * sih_111[k]
                   - f_10 * pc_y[k] * sii1_150[k];

        t_235[k] = f_12 * sih_90[k]
                   + f_3 * pc_z[k] * skh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_y, pc_x, pc_y, sii0_152, sii0_154, \
                         sih_113, sih_114, sih_183, sii1_152, sii1_154, skh_177, \
                         skh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_y[k] * sii0_152[k]
                   + f_12 * sih_113[k]
                   - f_10 * pc_y[k] * sii1_152[k];

        t_237[k] = f_11 * sih_114[k]
                   + f_3 * pc_y[k] * skh_177[k];

        t_238[k] = pb_y[k] * sii0_154[k]
                   - f_10 * pc_y[k] * sii1_154[k];

        t_239[k] = f_14 * sih_183[k]
                   + f_3 * pc_x[k] * skh_183[k];
    }
}

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;

    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);
    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_167 = buffer.data(sii0 + 167);
    const auto *sii0_168 = buffer.data(sii0 + 168);
    const auto *sii0_171 = buffer.data(sii0 + 171);
    const auto *sii0_174 = buffer.data(sii0 + 174);
    const auto *sii0_178 = buffer.data(sii0 + 178);
    const auto *sii0_180 = buffer.data(sii0 + 180);
    const auto *sii0_189 = buffer.data(sii0 + 189);

    const auto *sih_99 = buffer.data(sih + 99);
    const auto *sih_105 = buffer.data(sih + 105);
    const auto *sih_108 = buffer.data(sih + 108);
    const auto *sih_111 = buffer.data(sih + 111);
    const auto *sih_120 = buffer.data(sih + 120);
    const auto *sih_122 = buffer.data(sih + 122);
    const auto *sih_123 = buffer.data(sih + 123);
    const auto *sih_124 = buffer.data(sih + 124);
    const auto *sih_125 = buffer.data(sih + 125);
    const auto *sih_126 = buffer.data(sih + 126);
    const auto *sih_128 = buffer.data(sih + 128);
    const auto *sih_129 = buffer.data(sih + 129);
    const auto *sih_131 = buffer.data(sih + 131);
    const auto *sih_132 = buffer.data(sih + 132);
    const auto *sih_133 = buffer.data(sih + 133);
    const auto *sih_135 = buffer.data(sih + 135);
    const auto *sih_141 = buffer.data(sih + 141);
    const auto *sih_143 = buffer.data(sih + 143);
    const auto *sih_144 = buffer.data(sih + 144);
    const auto *sih_145 = buffer.data(sih + 145);
    const auto *sih_146 = buffer.data(sih + 146);
    const auto *sih_147 = buffer.data(sih + 147);
    const auto *sih_149 = buffer.data(sih + 149);
    const auto *sih_150 = buffer.data(sih + 150);
    const auto *sih_152 = buffer.data(sih + 152);
    const auto *sih_153 = buffer.data(sih + 153);
    const auto *sih_156 = buffer.data(sih + 156);
    const auto *sih_164 = buffer.data(sih + 164);
    const auto *sih_165 = buffer.data(sih + 165);
    const auto *sih_166 = buffer.data(sih + 166);
    const auto *sih_167 = buffer.data(sih + 167);
    const auto *sih_168 = buffer.data(sih + 168);
    const auto *sih_170 = buffer.data(sih + 170);
    const auto *sih_173 = buffer.data(sih + 173);
    const auto *sih_177 = buffer.data(sih + 177);
    const auto *sih_184 = buffer.data(sih + 184);
    const auto *sih_185 = buffer.data(sih + 185);
    const auto *sih_186 = buffer.data(sih + 186);
    const auto *sih_187 = buffer.data(sih + 187);
    const auto *sih_188 = buffer.data(sih + 188);
    const auto *sih_189 = buffer.data(sih + 189);
    const auto *sih_192 = buffer.data(sih + 192);
    const auto *sih_194 = buffer.data(sih + 194);
    const auto *sih_195 = buffer.data(sih + 195);
    const auto *sih_198 = buffer.data(sih + 198);
    const auto *sih_199 = buffer.data(sih + 199);
    const auto *sih_201 = buffer.data(sih + 201);
    const auto *sih_203 = buffer.data(sih + 203);
    const auto *sih_204 = buffer.data(sih + 204);
    const auto *sih_205 = buffer.data(sih + 205);
    const auto *sih_206 = buffer.data(sih + 206);
    const auto *sih_207 = buffer.data(sih + 207);
    const auto *sih_208 = buffer.data(sih + 208);
    const auto *sih_209 = buffer.data(sih + 209);
    const auto *sih_210 = buffer.data(sih + 210);
    const auto *sih_213 = buffer.data(sih + 213);
    const auto *sih_215 = buffer.data(sih + 215);
    const auto *sih_216 = buffer.data(sih + 216);
    const auto *sih_219 = buffer.data(sih + 219);
    const auto *sih_220 = buffer.data(sih + 220);
    const auto *sih_222 = buffer.data(sih + 222);
    const auto *sih_224 = buffer.data(sih + 224);
    const auto *sih_225 = buffer.data(sih + 225);
    const auto *sih_226 = buffer.data(sih + 226);
    const auto *sih_227 = buffer.data(sih + 227);
    const auto *sih_228 = buffer.data(sih + 228);
    const auto *sih_229 = buffer.data(sih + 229);
    const auto *sih_230 = buffer.data(sih + 230);
    const auto *sih_236 = buffer.data(sih + 236);
    const auto *sih_240 = buffer.data(sih + 240);
    const auto *sih_245 = buffer.data(sih + 245);
    const auto *sih_246 = buffer.data(sih + 246);
    const auto *sih_247 = buffer.data(sih + 247);
    const auto *sih_248 = buffer.data(sih + 248);
    const auto *sih_249 = buffer.data(sih + 249);
    const auto *sih_250 = buffer.data(sih + 250);
    const auto *sih_251 = buffer.data(sih + 251);
    const auto *sih_252 = buffer.data(sih + 252);
    const auto *sih_255 = buffer.data(sih + 255);
    const auto *sih_257 = buffer.data(sih + 257);
    const auto *sih_258 = buffer.data(sih + 258);
    const auto *sih_261 = buffer.data(sih + 261);
    const auto *sih_262 = buffer.data(sih + 262);
    const auto *sih_264 = buffer.data(sih + 264);
    const auto *sih_266 = buffer.data(sih + 266);

    const auto *sii1_167 = buffer.data(sii1 + 167);
    const auto *sii1_168 = buffer.data(sii1 + 168);
    const auto *sii1_171 = buffer.data(sii1 + 171);
    const auto *sii1_174 = buffer.data(sii1 + 174);
    const auto *sii1_178 = buffer.data(sii1 + 178);
    const auto *sii1_180 = buffer.data(sii1 + 180);
    const auto *sii1_189 = buffer.data(sii1 + 189);

    const auto *skg0_130 = buffer.data(skg0 + 130);
    const auto *skg0_132 = buffer.data(skg0 + 132);
    const auto *skg0_133 = buffer.data(skg0 + 133);
    const auto *skg0_134 = buffer.data(skg0 + 134);
    const auto *skg0_135 = buffer.data(skg0 + 135);
    const auto *skg0_138 = buffer.data(skg0 + 138);
    const auto *skg0_140 = buffer.data(skg0 + 140);
    const auto *skg0_141 = buffer.data(skg0 + 141);
    const auto *skg0_144 = buffer.data(skg0 + 144);
    const auto *skg0_145 = buffer.data(skg0 + 145);
    const auto *skg0_147 = buffer.data(skg0 + 147);
    const auto *skg0_148 = buffer.data(skg0 + 148);
    const auto *skg0_149 = buffer.data(skg0 + 149);
    const auto *skg0_150 = buffer.data(skg0 + 150);
    const auto *skg0_153 = buffer.data(skg0 + 153);
    const auto *skg0_155 = buffer.data(skg0 + 155);
    const auto *skg0_156 = buffer.data(skg0 + 156);
    const auto *skg0_159 = buffer.data(skg0 + 159);
    const auto *skg0_160 = buffer.data(skg0 + 160);
    const auto *skg0_162 = buffer.data(skg0 + 162);
    const auto *skg0_163 = buffer.data(skg0 + 163);
    const auto *skg0_164 = buffer.data(skg0 + 164);
    const auto *skg0_170 = buffer.data(skg0 + 170);
    const auto *skg0_174 = buffer.data(skg0 + 174);
    const auto *skg0_177 = buffer.data(skg0 + 177);
    const auto *skg0_178 = buffer.data(skg0 + 178);
    const auto *skg0_179 = buffer.data(skg0 + 179);
    const auto *skg0_180 = buffer.data(skg0 + 180);
    const auto *skg0_183 = buffer.data(skg0 + 183);
    const auto *skg0_185 = buffer.data(skg0 + 185);
    const auto *skg0_186 = buffer.data(skg0 + 186);
    const auto *skg0_189 = buffer.data(skg0 + 189);
    const auto *skg0_190 = buffer.data(skg0 + 190);
    const auto *skg0_192 = buffer.data(skg0 + 192);
    const auto *skg0_194 = buffer.data(skg0 + 194);

    const auto *skg1_130 = buffer.data(skg1 + 130);
    const auto *skg1_132 = buffer.data(skg1 + 132);
    const auto *skg1_133 = buffer.data(skg1 + 133);
    const auto *skg1_134 = buffer.data(skg1 + 134);
    const auto *skg1_135 = buffer.data(skg1 + 135);
    const auto *skg1_138 = buffer.data(skg1 + 138);
    const auto *skg1_140 = buffer.data(skg1 + 140);
    const auto *skg1_141 = buffer.data(skg1 + 141);
    const auto *skg1_144 = buffer.data(skg1 + 144);
    const auto *skg1_145 = buffer.data(skg1 + 145);
    const auto *skg1_147 = buffer.data(skg1 + 147);
    const auto *skg1_148 = buffer.data(skg1 + 148);
    const auto *skg1_149 = buffer.data(skg1 + 149);
    const auto *skg1_150 = buffer.data(skg1 + 150);
    const auto *skg1_153 = buffer.data(skg1 + 153);
    const auto *skg1_155 = buffer.data(skg1 + 155);
    const auto *skg1_156 = buffer.data(skg1 + 156);
    const auto *skg1_159 = buffer.data(skg1 + 159);
    const auto *skg1_160 = buffer.data(skg1 + 160);
    const auto *skg1_162 = buffer.data(skg1 + 162);
    const auto *skg1_163 = buffer.data(skg1 + 163);
    const auto *skg1_164 = buffer.data(skg1 + 164);
    const auto *skg1_170 = buffer.data(skg1 + 170);
    const auto *skg1_174 = buffer.data(skg1 + 174);
    const auto *skg1_177 = buffer.data(skg1 + 177);
    const auto *skg1_178 = buffer.data(skg1 + 178);
    const auto *skg1_179 = buffer.data(skg1 + 179);
    const auto *skg1_180 = buffer.data(skg1 + 180);
    const auto *skg1_183 = buffer.data(skg1 + 183);
    const auto *skg1_185 = buffer.data(skg1 + 185);
    const auto *skg1_186 = buffer.data(skg1 + 186);
    const auto *skg1_189 = buffer.data(skg1 + 189);
    const auto *skg1_190 = buffer.data(skg1 + 190);
    const auto *skg1_192 = buffer.data(skg1 + 192);
    const auto *skg1_194 = buffer.data(skg1 + 194);

    const auto *skh_183 = buffer.data(skh + 183);
    const auto *skh_184 = buffer.data(skh + 184);
    const auto *skh_185 = buffer.data(skh + 185);
    const auto *skh_186 = buffer.data(skh + 186);
    const auto *skh_187 = buffer.data(skh + 187);
    const auto *skh_188 = buffer.data(skh + 188);
    const auto *skh_189 = buffer.data(skh + 189);
    const auto *skh_191 = buffer.data(skh + 191);
    const auto *skh_192 = buffer.data(skh + 192);
    const auto *skh_194 = buffer.data(skh + 194);
    const auto *skh_195 = buffer.data(skh + 195);
    const auto *skh_198 = buffer.data(skh + 198);
    const auto *skh_199 = buffer.data(skh + 199);
    const auto *skh_201 = buffer.data(skh + 201);
    const auto *skh_203 = buffer.data(skh + 203);
    const auto *skh_204 = buffer.data(skh + 204);
    const auto *skh_205 = buffer.data(skh + 205);
    const auto *skh_206 = buffer.data(skh + 206);
    const auto *skh_207 = buffer.data(skh + 207);
    const auto *skh_208 = buffer.data(skh + 208);
    const auto *skh_209 = buffer.data(skh + 209);
    const auto *skh_210 = buffer.data(skh + 210);
    const auto *skh_212 = buffer.data(skh + 212);
    const auto *skh_213 = buffer.data(skh + 213);
    const auto *skh_215 = buffer.data(skh + 215);
    const auto *skh_216 = buffer.data(skh + 216);
    const auto *skh_219 = buffer.data(skh + 219);
    const auto *skh_220 = buffer.data(skh + 220);
    const auto *skh_222 = buffer.data(skh + 222);
    const auto *skh_224 = buffer.data(skh + 224);
    const auto *skh_225 = buffer.data(skh + 225);
    const auto *skh_226 = buffer.data(skh + 226);
    const auto *skh_227 = buffer.data(skh + 227);
    const auto *skh_228 = buffer.data(skh + 228);
    const auto *skh_229 = buffer.data(skh + 229);
    const auto *skh_230 = buffer.data(skh + 230);
    const auto *skh_231 = buffer.data(skh + 231);
    const auto *skh_233 = buffer.data(skh + 233);
    const auto *skh_234 = buffer.data(skh + 234);
    const auto *skh_236 = buffer.data(skh + 236);
    const auto *skh_237 = buffer.data(skh + 237);
    const auto *skh_240 = buffer.data(skh + 240);
    const auto *skh_245 = buffer.data(skh + 245);
    const auto *skh_246 = buffer.data(skh + 246);
    const auto *skh_247 = buffer.data(skh + 247);
    const auto *skh_248 = buffer.data(skh + 248);
    const auto *skh_249 = buffer.data(skh + 249);
    const auto *skh_250 = buffer.data(skh + 250);
    const auto *skh_251 = buffer.data(skh + 251);
    const auto *skh_252 = buffer.data(skh + 252);
    const auto *skh_254 = buffer.data(skh + 254);
    const auto *skh_255 = buffer.data(skh + 255);
    const auto *skh_257 = buffer.data(skh + 257);
    const auto *skh_258 = buffer.data(skh + 258);
    const auto *skh_261 = buffer.data(skh + 261);
    const auto *skh_262 = buffer.data(skh + 262);
    const auto *skh_264 = buffer.data(skh + 264);
    const auto *skh_266 = buffer.data(skh + 266);

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, sih_184, sih_185, sih_186, \
                         sih_187, sih_188, skh_184, skh_185, skh_186, skh_187, \
                         skh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_14 * sih_184[k]
                   + f_3 * pc_x[k] * skh_184[k];

        t_241[k] = f_14 * sih_185[k]
                   + f_3 * pc_x[k] * skh_185[k];

        t_242[k] = f_14 * sih_186[k]
                   + f_3 * pc_x[k] * skh_186[k];

        t_243[k] = f_14 * sih_187[k]
                   + f_3 * pc_x[k] * skh_187[k];

        t_244[k] = f_14 * sih_188[k]
                   + f_3 * pc_x[k] * skh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, sih_99, sih_120, sih_122, skg0_130, \
                         skg0_132, skg1_130, skg1_132, skh_183, \
                         skh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * sih_120[k]
                   + f_1 * skg0_130[k]
                   - f_2 * skg1_130[k]
                   + f_3 * pc_y[k] * skh_183[k];

        t_246[k] = f_12 * sih_99[k]
                   + f_3 * pc_z[k] * skh_183[k];

        t_247[k] = f_11 * sih_122[k]
                   + f_4 * skg0_132[k]
                   - f_5 * skg1_132[k]
                   + f_3 * pc_y[k] * skh_185[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, sih_123, sih_124, sih_125, skg0_133, \
                         skg0_134, skg1_133, skg1_134, skh_186, skh_187, \
                         skh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * sih_123[k]
                   + f_6 * skg0_133[k]
                   - f_7 * skg1_133[k]
                   + f_3 * pc_y[k] * skh_186[k];

        t_249[k] = f_11 * sih_124[k]
                   + f_8 * skg0_134[k]
                   - f_9 * skg1_134[k]
                   + f_3 * pc_y[k] * skh_187[k];

        t_250[k] = f_11 * sih_125[k]
                   + f_3 * pc_y[k] * skh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_x, pc_y, pc_z, sii0_167, \
                         sih_105, sih_189, sii1_167, skg0_135, skg1_135, \
                         skh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * sii0_167[k]
                   - f_10 * pc_y[k] * sii1_167[k];

        t_252[k] = f_14 * sih_189[k]
                   + f_1 * skg0_135[k]
                   - f_2 * skg1_135[k]
                   + f_3 * pc_x[k] * skh_189[k];

        t_253[k] = f_3 * pc_y[k] * skh_189[k];

        t_254[k] = f_13 * sih_105[k]
                   + f_3 * pc_z[k] * skh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, sih_192, sih_194, skg0_138, \
                         skg0_140, skg1_138, skg1_140, skh_191, skh_192, \
                         skh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_14 * sih_192[k]
                   + f_4 * skg0_138[k]
                   - f_5 * skg1_138[k]
                   + f_3 * pc_x[k] * skh_192[k];

        t_256[k] = f_3 * pc_y[k] * skh_191[k];

        t_257[k] = f_14 * sih_194[k]
                   + f_4 * skg0_140[k]
                   - f_5 * skg1_140[k]
                   + f_3 * pc_x[k] * skh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, pc_z, sih_108, sih_195, skg0_141, \
                         skg1_141, skh_192, skh_194, skh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_14 * sih_195[k]
                   + f_6 * skg0_141[k]
                   - f_7 * skg1_141[k]
                   + f_3 * pc_x[k] * skh_195[k];

        t_259[k] = f_13 * sih_108[k]
                   + f_3 * pc_z[k] * skh_192[k];

        t_260[k] = f_3 * pc_y[k] * skh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, sih_111, sih_198, sih_199, skg0_144, \
                         skg0_145, skg1_144, skg1_145, skh_195, skh_198, \
                         skh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_14 * sih_198[k]
                   + f_6 * skg0_144[k]
                   - f_7 * skg1_144[k]
                   + f_3 * pc_x[k] * skh_198[k];

        t_262[k] = f_14 * sih_199[k]
                   + f_8 * skg0_145[k]
                   - f_9 * skg1_145[k]
                   + f_3 * pc_x[k] * skh_199[k];

        t_263[k] = f_13 * sih_111[k]
                   + f_3 * pc_z[k] * skh_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, sih_201, sih_203, skg0_147, \
                         skg0_149, skg1_147, skg1_149, skh_198, skh_201, \
                         skh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_14 * sih_201[k]
                   + f_8 * skg0_147[k]
                   - f_9 * skg1_147[k]
                   + f_3 * pc_x[k] * skh_201[k];

        t_265[k] = f_3 * pc_y[k] * skh_198[k];

        t_266[k] = f_14 * sih_203[k]
                   + f_8 * skg0_149[k]
                   - f_9 * skg1_149[k]
                   + f_3 * pc_x[k] * skh_203[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pc_x, sih_204, sih_205, sih_206, \
                         sih_207, sih_208, skh_204, skh_205, skh_206, skh_207, \
                         skh_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * sih_204[k]
                   + f_3 * pc_x[k] * skh_204[k];

        t_268[k] = f_14 * sih_205[k]
                   + f_3 * pc_x[k] * skh_205[k];

        t_269[k] = f_14 * sih_206[k]
                   + f_3 * pc_x[k] * skh_206[k];

        t_270[k] = f_14 * sih_207[k]
                   + f_3 * pc_x[k] * skh_207[k];

        t_271[k] = f_14 * sih_208[k]
                   + f_3 * pc_x[k] * skh_208[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, sih_120, sih_209, \
                         skg0_145, skg0_147, skg1_145, skg1_147, skh_204, skh_206, \
                         skh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_14 * sih_209[k]
                   + f_3 * pc_x[k] * skh_209[k];

        t_273[k] = f_1 * skg0_145[k]
                   - f_2 * skg1_145[k]
                   + f_3 * pc_y[k] * skh_204[k];

        t_274[k] = f_13 * sih_120[k]
                   + f_3 * pc_z[k] * skh_204[k];

        t_275[k] = f_4 * skg0_147[k]
                   - f_5 * skg1_147[k]
                   + f_3 * pc_y[k] * skh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, sih_125, skg0_148, skg0_149, \
                         skg1_148, skg1_149, skh_207, skh_208, \
                         skh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * skg0_148[k]
                   - f_7 * skg1_148[k]
                   + f_3 * pc_y[k] * skh_207[k];

        t_277[k] = f_8 * skg0_149[k]
                   - f_9 * skg1_149[k]
                   + f_3 * pc_y[k] * skh_208[k];

        t_278[k] = f_3 * pc_y[k] * skh_209[k];

        t_279[k] = f_13 * sih_125[k]
                   + f_1 * skg0_149[k]
                   - f_2 * skg1_149[k]
                   + f_3 * pc_z[k] * skh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, sih_126, sih_210, \
                         sih_213, skg0_150, skg0_153, skg1_150, skg1_153, skh_210, \
                         skh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_13 * sih_210[k]
                   + f_1 * skg0_150[k]
                   - f_2 * skg1_150[k]
                   + f_3 * pc_x[k] * skh_210[k];

        t_281[k] = f_14 * sih_126[k]
                   + f_3 * pc_y[k] * skh_210[k];

        t_282[k] = f_3 * pc_z[k] * skh_210[k];

        t_283[k] = f_13 * sih_213[k]
                   + f_4 * skg0_153[k]
                   - f_5 * skg1_153[k]
                   + f_3 * pc_x[k] * skh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_y, sih_128, sih_215, sih_216, skg0_155, \
                         skg0_156, skg1_155, skg1_156, skh_212, skh_215, \
                         skh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_14 * sih_128[k]
                   + f_3 * pc_y[k] * skh_212[k];

        t_285[k] = f_13 * sih_215[k]
                   + f_4 * skg0_155[k]
                   - f_5 * skg1_155[k]
                   + f_3 * pc_x[k] * skh_215[k];

        t_286[k] = f_13 * sih_216[k]
                   + f_6 * skg0_156[k]
                   - f_7 * skg1_156[k]
                   + f_3 * pc_x[k] * skh_216[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, pc_z, sih_131, sih_219, skg0_159, \
                         skg1_159, skh_213, skh_215, skh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * skh_213[k];

        t_288[k] = f_14 * sih_131[k]
                   + f_3 * pc_y[k] * skh_215[k];

        t_289[k] = f_13 * sih_219[k]
                   + f_6 * skg0_159[k]
                   - f_7 * skg1_159[k]
                   + f_3 * pc_x[k] * skh_219[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, sih_220, sih_222, skg0_160, \
                         skg0_162, skg1_160, skg1_162, skh_216, skh_220, \
                         skh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_13 * sih_220[k]
                   + f_8 * skg0_160[k]
                   - f_9 * skg1_160[k]
                   + f_3 * pc_x[k] * skh_220[k];

        t_291[k] = f_3 * pc_z[k] * skh_216[k];

        t_292[k] = f_13 * sih_222[k]
                   + f_8 * skg0_162[k]
                   - f_9 * skg1_162[k]
                   + f_3 * pc_x[k] * skh_222[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pc_x, pc_y, sih_135, sih_224, sih_225, \
                         sih_226, skg0_164, skg1_164, skh_219, skh_224, skh_225, \
                         skh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_14 * sih_135[k]
                   + f_3 * pc_y[k] * skh_219[k];

        t_294[k] = f_13 * sih_224[k]
                   + f_8 * skg0_164[k]
                   - f_9 * skg1_164[k]
                   + f_3 * pc_x[k] * skh_224[k];

        t_295[k] = f_13 * sih_225[k]
                   + f_3 * pc_x[k] * skh_225[k];

        t_296[k] = f_13 * sih_226[k]
                   + f_3 * pc_x[k] * skh_226[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pc_x, sih_227, sih_228, sih_229, sih_230, \
                         skh_227, skh_228, skh_229, skh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_13 * sih_227[k]
                   + f_3 * pc_x[k] * skh_227[k];

        t_298[k] = f_13 * sih_228[k]
                   + f_3 * pc_x[k] * skh_228[k];

        t_299[k] = f_13 * sih_229[k]
                   + f_3 * pc_x[k] * skh_229[k];

        t_300[k] = f_13 * sih_230[k]
                   + f_3 * pc_x[k] * skh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, pc_y, pc_z, sih_141, sih_143, skg0_160, \
                         skg0_162, skg1_160, skg1_162, skh_225, \
                         skh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * sih_141[k]
                   + f_1 * skg0_160[k]
                   - f_2 * skg1_160[k]
                   + f_3 * pc_y[k] * skh_225[k];

        t_302[k] = f_3 * pc_z[k] * skh_225[k];

        t_303[k] = f_14 * sih_143[k]
                   + f_4 * skg0_162[k]
                   - f_5 * skg1_162[k]
                   + f_3 * pc_y[k] * skh_227[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_y, pc_z, sih_144, sih_145, sih_146, \
                         skg0_163, skg0_164, skg1_163, skg1_164, skh_228, skh_229, \
                         skh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_14 * sih_144[k]
                   + f_6 * skg0_163[k]
                   - f_7 * skg1_163[k]
                   + f_3 * pc_y[k] * skh_228[k];

        t_305[k] = f_14 * sih_145[k]
                   + f_8 * skg0_164[k]
                   - f_9 * skg1_164[k]
                   + f_3 * pc_y[k] * skh_229[k];

        t_306[k] = f_14 * sih_146[k]
                   + f_3 * pc_y[k] * skh_230[k];

        t_307[k] = f_1 * skg0_164[k]
                   - f_2 * skg1_164[k]
                   + f_3 * pc_z[k] * skh_230[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, sii0_168, sii0_171, \
                         sih_126, sih_147, sii1_168, sii1_171, \
                         skh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * sii0_168[k]
                   - f_10 * pc_z[k] * sii1_168[k];

        t_309[k] = f_13 * sih_147[k]
                   + f_3 * pc_y[k] * skh_231[k];

        t_310[k] = f_11 * sih_126[k]
                   + f_3 * pc_z[k] * skh_231[k];

        t_311[k] = pb_z[k] * sii0_171[k]
                   - f_10 * pc_z[k] * sii1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_z, pc_x, pc_y, pc_z, sii0_174, sih_149, \
                         sih_236, sii1_174, skg0_170, skg1_170, skh_233, \
                         skh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * sih_149[k]
                   + f_3 * pc_y[k] * skh_233[k];

        t_313[k] = f_13 * sih_236[k]
                   + f_4 * skg0_170[k]
                   - f_5 * skg1_170[k]
                   + f_3 * pc_x[k] * skh_236[k];

        t_314[k] = pb_z[k] * sii0_174[k]
                   - f_10 * pc_z[k] * sii1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pc_x, pc_y, pc_z, sih_129, sih_152, sih_240, \
                         skg0_174, skg1_174, skh_234, skh_236, \
                         skh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * sih_129[k]
                   + f_3 * pc_z[k] * skh_234[k];

        t_316[k] = f_13 * sih_152[k]
                   + f_3 * pc_y[k] * skh_236[k];

        t_317[k] = f_13 * sih_240[k]
                   + f_6 * skg0_174[k]
                   - f_7 * skg1_174[k]
                   + f_3 * pc_x[k] * skh_240[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_z, pc_y, pc_z, sii0_178, sii0_180, \
                         sih_132, sih_133, sih_156, sii1_178, sii1_180, skh_237, \
                         skh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * sii0_178[k]
                   - f_10 * pc_z[k] * sii1_178[k];

        t_319[k] = f_11 * sih_132[k]
                   + f_3 * pc_z[k] * skh_237[k];

        t_320[k] = pb_z[k] * sii0_180[k]
                   + f_12 * sih_133[k]
                   - f_10 * pc_z[k] * sii1_180[k];

        t_321[k] = f_13 * sih_156[k]
                   + f_3 * pc_y[k] * skh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, sih_245, sih_246, sih_247, sih_248, \
                         skg0_179, skg1_179, skh_245, skh_246, skh_247, \
                         skh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_13 * sih_245[k]
                   + f_8 * skg0_179[k]
                   - f_9 * skg1_179[k]
                   + f_3 * pc_x[k] * skh_245[k];

        t_323[k] = f_13 * sih_246[k]
                   + f_3 * pc_x[k] * skh_246[k];

        t_324[k] = f_13 * sih_247[k]
                   + f_3 * pc_x[k] * skh_247[k];

        t_325[k] = f_13 * sih_248[k]
                   + f_3 * pc_x[k] * skh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pb_z, pc_x, pc_z, sii0_189, sih_249, \
                         sih_250, sih_251, sii1_189, skh_249, skh_250, \
                         skh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_13 * sih_249[k]
                   + f_3 * pc_x[k] * skh_249[k];

        t_327[k] = f_13 * sih_250[k]
                   + f_3 * pc_x[k] * skh_250[k];

        t_328[k] = f_13 * sih_251[k]
                   + f_3 * pc_x[k] * skh_251[k];

        t_329[k] = pb_z[k] * sii0_189[k]
                   - f_10 * pc_z[k] * sii1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, sih_141, sih_164, sih_165, skg0_177, \
                         skg0_178, skg1_177, skg1_178, skh_246, skh_248, \
                         skh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * sih_141[k]
                   + f_3 * pc_z[k] * skh_246[k];

        t_331[k] = f_13 * sih_164[k]
                   + f_4 * skg0_177[k]
                   - f_5 * skg1_177[k]
                   + f_3 * pc_y[k] * skh_248[k];

        t_332[k] = f_13 * sih_165[k]
                   + f_6 * skg0_178[k]
                   - f_7 * skg1_178[k]
                   + f_3 * pc_y[k] * skh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, sih_146, sih_166, sih_167, skg0_179, \
                         skg1_179, skh_250, skh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * sih_166[k]
                   + f_8 * skg0_179[k]
                   - f_9 * skg1_179[k]
                   + f_3 * pc_y[k] * skh_250[k];

        t_334[k] = f_13 * sih_167[k]
                   + f_3 * pc_y[k] * skh_251[k];

        t_335[k] = f_11 * sih_146[k]
                   + f_1 * skg0_179[k]
                   - f_2 * skg1_179[k]
                   + f_3 * pc_z[k] * skh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, sih_147, sih_168, sih_252, \
                         skg0_180, skg1_180, skh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_13 * sih_252[k]
                   + f_1 * skg0_180[k]
                   - f_2 * skg1_180[k]
                   + f_3 * pc_x[k] * skh_252[k];

        t_337[k] = f_12 * sih_168[k]
                   + f_3 * pc_y[k] * skh_252[k];

        t_338[k] = f_12 * sih_147[k]
                   + f_3 * pc_z[k] * skh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, sih_170, sih_255, sih_257, skg0_183, \
                         skg0_185, skg1_183, skg1_185, skh_254, skh_255, \
                         skh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_13 * sih_255[k]
                   + f_4 * skg0_183[k]
                   - f_5 * skg1_183[k]
                   + f_3 * pc_x[k] * skh_255[k];

        t_340[k] = f_12 * sih_170[k]
                   + f_3 * pc_y[k] * skh_254[k];

        t_341[k] = f_13 * sih_257[k]
                   + f_4 * skg0_185[k]
                   - f_5 * skg1_185[k]
                   + f_3 * pc_x[k] * skh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, sih_150, sih_173, sih_258, \
                         skg0_186, skg1_186, skh_255, skh_257, \
                         skh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_13 * sih_258[k]
                   + f_6 * skg0_186[k]
                   - f_7 * skg1_186[k]
                   + f_3 * pc_x[k] * skh_258[k];

        t_343[k] = f_12 * sih_150[k]
                   + f_3 * pc_z[k] * skh_255[k];

        t_344[k] = f_12 * sih_173[k]
                   + f_3 * pc_y[k] * skh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, sih_153, sih_261, sih_262, skg0_189, \
                         skg0_190, skg1_189, skg1_190, skh_258, skh_261, \
                         skh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_13 * sih_261[k]
                   + f_6 * skg0_189[k]
                   - f_7 * skg1_189[k]
                   + f_3 * pc_x[k] * skh_261[k];

        t_346[k] = f_13 * sih_262[k]
                   + f_8 * skg0_190[k]
                   - f_9 * skg1_190[k]
                   + f_3 * pc_x[k] * skh_262[k];

        t_347[k] = f_12 * sih_153[k]
                   + f_3 * pc_z[k] * skh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, sih_177, sih_264, sih_266, skg0_192, \
                         skg0_194, skg1_192, skg1_194, skh_261, skh_264, \
                         skh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_13 * sih_264[k]
                   + f_8 * skg0_192[k]
                   - f_9 * skg1_192[k]
                   + f_3 * pc_x[k] * skh_264[k];

        t_349[k] = f_12 * sih_177[k]
                   + f_3 * pc_y[k] * skh_261[k];

        t_350[k] = f_13 * sih_266[k]
                   + f_8 * skg0_194[k]
                   - f_9 * skg1_194[k]
                   + f_3 * pc_x[k] * skh_266[k];
    }
}

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.5 / q;

    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_252 = buffer.data(sii0 + 252);
    const auto *sii0_255 = buffer.data(sii0 + 255);
    const auto *sii0_257 = buffer.data(sii0 + 257);
    const auto *sii0_258 = buffer.data(sii0 + 258);
    const auto *sii0_261 = buffer.data(sii0 + 261);
    const auto *sii0_262 = buffer.data(sii0 + 262);
    const auto *sii0_264 = buffer.data(sii0 + 264);
    const auto *sii0_266 = buffer.data(sii0 + 266);
    const auto *sii0_279 = buffer.data(sii0 + 279);
    const auto *sii0_280 = buffer.data(sii0 + 280);
    const auto *sii0_283 = buffer.data(sii0 + 283);
    const auto *sii0_286 = buffer.data(sii0 + 286);
    const auto *sii0_290 = buffer.data(sii0 + 290);
    const auto *sii0_292 = buffer.data(sii0 + 292);

    const auto *sih_162 = buffer.data(sih + 162);
    const auto *sih_167 = buffer.data(sih + 167);
    const auto *sih_168 = buffer.data(sih + 168);
    const auto *sih_171 = buffer.data(sih + 171);
    const auto *sih_174 = buffer.data(sih + 174);
    const auto *sih_183 = buffer.data(sih + 183);
    const auto *sih_185 = buffer.data(sih + 185);
    const auto *sih_186 = buffer.data(sih + 186);
    const auto *sih_187 = buffer.data(sih + 187);
    const auto *sih_188 = buffer.data(sih + 188);
    const auto *sih_189 = buffer.data(sih + 189);
    const auto *sih_190 = buffer.data(sih + 190);
    const auto *sih_191 = buffer.data(sih + 191);
    const auto *sih_192 = buffer.data(sih + 192);
    const auto *sih_194 = buffer.data(sih + 194);
    const auto *sih_195 = buffer.data(sih + 195);
    const auto *sih_197 = buffer.data(sih + 197);
    const auto *sih_198 = buffer.data(sih + 198);
    const auto *sih_204 = buffer.data(sih + 204);
    const auto *sih_206 = buffer.data(sih + 206);
    const auto *sih_207 = buffer.data(sih + 207);
    const auto *sih_208 = buffer.data(sih + 208);
    const auto *sih_209 = buffer.data(sih + 209);
    const auto *sih_210 = buffer.data(sih + 210);
    const auto *sih_212 = buffer.data(sih + 212);
    const auto *sih_213 = buffer.data(sih + 213);
    const auto *sih_215 = buffer.data(sih + 215);
    const auto *sih_216 = buffer.data(sih + 216);
    const auto *sih_217 = buffer.data(sih + 217);
    const auto *sih_219 = buffer.data(sih + 219);
    const auto *sih_225 = buffer.data(sih + 225);
    const auto *sih_227 = buffer.data(sih + 227);
    const auto *sih_228 = buffer.data(sih + 228);
    const auto *sih_229 = buffer.data(sih + 229);
    const auto *sih_230 = buffer.data(sih + 230);
    const auto *sih_231 = buffer.data(sih + 231);
    const auto *sih_233 = buffer.data(sih + 233);
    const auto *sih_236 = buffer.data(sih + 236);
    const auto *sih_240 = buffer.data(sih + 240);
    const auto *sih_267 = buffer.data(sih + 267);
    const auto *sih_268 = buffer.data(sih + 268);
    const auto *sih_269 = buffer.data(sih + 269);
    const auto *sih_270 = buffer.data(sih + 270);
    const auto *sih_271 = buffer.data(sih + 271);
    const auto *sih_272 = buffer.data(sih + 272);
    const auto *sih_288 = buffer.data(sih + 288);
    const auto *sih_289 = buffer.data(sih + 289);
    const auto *sih_290 = buffer.data(sih + 290);
    const auto *sih_291 = buffer.data(sih + 291);
    const auto *sih_292 = buffer.data(sih + 292);
    const auto *sih_293 = buffer.data(sih + 293);
    const auto *sih_294 = buffer.data(sih + 294);
    const auto *sih_297 = buffer.data(sih + 297);
    const auto *sih_299 = buffer.data(sih + 299);
    const auto *sih_300 = buffer.data(sih + 300);
    const auto *sih_303 = buffer.data(sih + 303);
    const auto *sih_304 = buffer.data(sih + 304);
    const auto *sih_306 = buffer.data(sih + 306);
    const auto *sih_308 = buffer.data(sih + 308);
    const auto *sih_309 = buffer.data(sih + 309);
    const auto *sih_310 = buffer.data(sih + 310);
    const auto *sih_311 = buffer.data(sih + 311);
    const auto *sih_312 = buffer.data(sih + 312);
    const auto *sih_313 = buffer.data(sih + 313);
    const auto *sih_314 = buffer.data(sih + 314);
    const auto *sih_315 = buffer.data(sih + 315);
    const auto *sih_318 = buffer.data(sih + 318);
    const auto *sih_320 = buffer.data(sih + 320);
    const auto *sih_321 = buffer.data(sih + 321);
    const auto *sih_324 = buffer.data(sih + 324);
    const auto *sih_325 = buffer.data(sih + 325);
    const auto *sih_327 = buffer.data(sih + 327);
    const auto *sih_329 = buffer.data(sih + 329);
    const auto *sih_330 = buffer.data(sih + 330);
    const auto *sih_331 = buffer.data(sih + 331);
    const auto *sih_332 = buffer.data(sih + 332);
    const auto *sih_333 = buffer.data(sih + 333);
    const auto *sih_334 = buffer.data(sih + 334);
    const auto *sih_335 = buffer.data(sih + 335);
    const auto *sih_341 = buffer.data(sih + 341);
    const auto *sih_345 = buffer.data(sih + 345);
    const auto *sih_350 = buffer.data(sih + 350);
    const auto *sih_351 = buffer.data(sih + 351);
    const auto *sih_352 = buffer.data(sih + 352);
    const auto *sih_353 = buffer.data(sih + 353);

    const auto *sii1_252 = buffer.data(sii1 + 252);
    const auto *sii1_255 = buffer.data(sii1 + 255);
    const auto *sii1_257 = buffer.data(sii1 + 257);
    const auto *sii1_258 = buffer.data(sii1 + 258);
    const auto *sii1_261 = buffer.data(sii1 + 261);
    const auto *sii1_262 = buffer.data(sii1 + 262);
    const auto *sii1_264 = buffer.data(sii1 + 264);
    const auto *sii1_266 = buffer.data(sii1 + 266);
    const auto *sii1_279 = buffer.data(sii1 + 279);
    const auto *sii1_280 = buffer.data(sii1 + 280);
    const auto *sii1_283 = buffer.data(sii1 + 283);
    const auto *sii1_286 = buffer.data(sii1 + 286);
    const auto *sii1_290 = buffer.data(sii1 + 290);
    const auto *sii1_292 = buffer.data(sii1 + 292);

    const auto *skg0_190 = buffer.data(skg0 + 190);
    const auto *skg0_192 = buffer.data(skg0 + 192);
    const auto *skg0_193 = buffer.data(skg0 + 193);
    const auto *skg0_194 = buffer.data(skg0 + 194);
    const auto *skg0_205 = buffer.data(skg0 + 205);
    const auto *skg0_207 = buffer.data(skg0 + 207);
    const auto *skg0_208 = buffer.data(skg0 + 208);
    const auto *skg0_209 = buffer.data(skg0 + 209);
    const auto *skg0_210 = buffer.data(skg0 + 210);
    const auto *skg0_213 = buffer.data(skg0 + 213);
    const auto *skg0_215 = buffer.data(skg0 + 215);
    const auto *skg0_216 = buffer.data(skg0 + 216);
    const auto *skg0_219 = buffer.data(skg0 + 219);
    const auto *skg0_220 = buffer.data(skg0 + 220);
    const auto *skg0_222 = buffer.data(skg0 + 222);
    const auto *skg0_223 = buffer.data(skg0 + 223);
    const auto *skg0_224 = buffer.data(skg0 + 224);
    const auto *skg0_225 = buffer.data(skg0 + 225);
    const auto *skg0_228 = buffer.data(skg0 + 228);
    const auto *skg0_230 = buffer.data(skg0 + 230);
    const auto *skg0_231 = buffer.data(skg0 + 231);
    const auto *skg0_234 = buffer.data(skg0 + 234);
    const auto *skg0_235 = buffer.data(skg0 + 235);
    const auto *skg0_237 = buffer.data(skg0 + 237);
    const auto *skg0_238 = buffer.data(skg0 + 238);
    const auto *skg0_239 = buffer.data(skg0 + 239);
    const auto *skg0_245 = buffer.data(skg0 + 245);
    const auto *skg0_249 = buffer.data(skg0 + 249);
    const auto *skg0_254 = buffer.data(skg0 + 254);

    const auto *skg1_190 = buffer.data(skg1 + 190);
    const auto *skg1_192 = buffer.data(skg1 + 192);
    const auto *skg1_193 = buffer.data(skg1 + 193);
    const auto *skg1_194 = buffer.data(skg1 + 194);
    const auto *skg1_205 = buffer.data(skg1 + 205);
    const auto *skg1_207 = buffer.data(skg1 + 207);
    const auto *skg1_208 = buffer.data(skg1 + 208);
    const auto *skg1_209 = buffer.data(skg1 + 209);
    const auto *skg1_210 = buffer.data(skg1 + 210);
    const auto *skg1_213 = buffer.data(skg1 + 213);
    const auto *skg1_215 = buffer.data(skg1 + 215);
    const auto *skg1_216 = buffer.data(skg1 + 216);
    const auto *skg1_219 = buffer.data(skg1 + 219);
    const auto *skg1_220 = buffer.data(skg1 + 220);
    const auto *skg1_222 = buffer.data(skg1 + 222);
    const auto *skg1_223 = buffer.data(skg1 + 223);
    const auto *skg1_224 = buffer.data(skg1 + 224);
    const auto *skg1_225 = buffer.data(skg1 + 225);
    const auto *skg1_228 = buffer.data(skg1 + 228);
    const auto *skg1_230 = buffer.data(skg1 + 230);
    const auto *skg1_231 = buffer.data(skg1 + 231);
    const auto *skg1_234 = buffer.data(skg1 + 234);
    const auto *skg1_235 = buffer.data(skg1 + 235);
    const auto *skg1_237 = buffer.data(skg1 + 237);
    const auto *skg1_238 = buffer.data(skg1 + 238);
    const auto *skg1_239 = buffer.data(skg1 + 239);
    const auto *skg1_245 = buffer.data(skg1 + 245);
    const auto *skg1_249 = buffer.data(skg1 + 249);
    const auto *skg1_254 = buffer.data(skg1 + 254);

    const auto *skh_267 = buffer.data(skh + 267);
    const auto *skh_268 = buffer.data(skh + 268);
    const auto *skh_269 = buffer.data(skh + 269);
    const auto *skh_270 = buffer.data(skh + 270);
    const auto *skh_271 = buffer.data(skh + 271);
    const auto *skh_272 = buffer.data(skh + 272);
    const auto *skh_273 = buffer.data(skh + 273);
    const auto *skh_275 = buffer.data(skh + 275);
    const auto *skh_276 = buffer.data(skh + 276);
    const auto *skh_278 = buffer.data(skh + 278);
    const auto *skh_279 = buffer.data(skh + 279);
    const auto *skh_282 = buffer.data(skh + 282);
    const auto *skh_288 = buffer.data(skh + 288);
    const auto *skh_289 = buffer.data(skh + 289);
    const auto *skh_290 = buffer.data(skh + 290);
    const auto *skh_291 = buffer.data(skh + 291);
    const auto *skh_292 = buffer.data(skh + 292);
    const auto *skh_293 = buffer.data(skh + 293);
    const auto *skh_294 = buffer.data(skh + 294);
    const auto *skh_296 = buffer.data(skh + 296);
    const auto *skh_297 = buffer.data(skh + 297);
    const auto *skh_299 = buffer.data(skh + 299);
    const auto *skh_300 = buffer.data(skh + 300);
    const auto *skh_303 = buffer.data(skh + 303);
    const auto *skh_304 = buffer.data(skh + 304);
    const auto *skh_306 = buffer.data(skh + 306);
    const auto *skh_308 = buffer.data(skh + 308);
    const auto *skh_309 = buffer.data(skh + 309);
    const auto *skh_310 = buffer.data(skh + 310);
    const auto *skh_311 = buffer.data(skh + 311);
    const auto *skh_312 = buffer.data(skh + 312);
    const auto *skh_313 = buffer.data(skh + 313);
    const auto *skh_314 = buffer.data(skh + 314);
    const auto *skh_315 = buffer.data(skh + 315);
    const auto *skh_317 = buffer.data(skh + 317);
    const auto *skh_318 = buffer.data(skh + 318);
    const auto *skh_320 = buffer.data(skh + 320);
    const auto *skh_321 = buffer.data(skh + 321);
    const auto *skh_324 = buffer.data(skh + 324);
    const auto *skh_325 = buffer.data(skh + 325);
    const auto *skh_327 = buffer.data(skh + 327);
    const auto *skh_329 = buffer.data(skh + 329);
    const auto *skh_330 = buffer.data(skh + 330);
    const auto *skh_331 = buffer.data(skh + 331);
    const auto *skh_332 = buffer.data(skh + 332);
    const auto *skh_333 = buffer.data(skh + 333);
    const auto *skh_334 = buffer.data(skh + 334);
    const auto *skh_335 = buffer.data(skh + 335);
    const auto *skh_336 = buffer.data(skh + 336);
    const auto *skh_338 = buffer.data(skh + 338);
    const auto *skh_339 = buffer.data(skh + 339);
    const auto *skh_341 = buffer.data(skh + 341);
    const auto *skh_342 = buffer.data(skh + 342);
    const auto *skh_345 = buffer.data(skh + 345);
    const auto *skh_350 = buffer.data(skh + 350);
    const auto *skh_351 = buffer.data(skh + 351);
    const auto *skh_352 = buffer.data(skh + 352);
    const auto *skh_353 = buffer.data(skh + 353);

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, sih_267, sih_268, sih_269, \
                         sih_270, sih_271, skh_267, skh_268, skh_269, skh_270, \
                         skh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_13 * sih_267[k]
                   + f_3 * pc_x[k] * skh_267[k];

        t_352[k] = f_13 * sih_268[k]
                   + f_3 * pc_x[k] * skh_268[k];

        t_353[k] = f_13 * sih_269[k]
                   + f_3 * pc_x[k] * skh_269[k];

        t_354[k] = f_13 * sih_270[k]
                   + f_3 * pc_x[k] * skh_270[k];

        t_355[k] = f_13 * sih_271[k]
                   + f_3 * pc_x[k] * skh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, sih_162, sih_183, sih_272, \
                         skg0_190, skg1_190, skh_267, skh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_13 * sih_272[k]
                   + f_3 * pc_x[k] * skh_272[k];

        t_357[k] = f_12 * sih_183[k]
                   + f_1 * skg0_190[k]
                   - f_2 * skg1_190[k]
                   + f_3 * pc_y[k] * skh_267[k];

        t_358[k] = f_12 * sih_162[k]
                   + f_3 * pc_z[k] * skh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, sih_185, sih_186, sih_187, skg0_192, \
                         skg0_193, skg0_194, skg1_192, skg1_193, skg1_194, skh_269, skh_270, \
                         skh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * sih_185[k]
                   + f_4 * skg0_192[k]
                   - f_5 * skg1_192[k]
                   + f_3 * pc_y[k] * skh_269[k];

        t_360[k] = f_12 * sih_186[k]
                   + f_6 * skg0_193[k]
                   - f_7 * skg1_193[k]
                   + f_3 * pc_y[k] * skh_270[k];

        t_361[k] = f_12 * sih_187[k]
                   + f_8 * skg0_194[k]
                   - f_9 * skg1_194[k]
                   + f_3 * pc_y[k] * skh_271[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pb_y, pc_y, pc_z, sii0_252, sih_167, \
                         sih_188, sih_189, sii1_252, skg0_194, skg1_194, skh_272, \
                         skh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * sih_188[k]
                   + f_3 * pc_y[k] * skh_272[k];

        t_363[k] = f_12 * sih_167[k]
                   + f_1 * skg0_194[k]
                   - f_2 * skg1_194[k]
                   + f_3 * pc_z[k] * skh_272[k];

        t_364[k] = pb_y[k] * sii0_252[k]
                   - f_10 * pc_y[k] * sii1_252[k];

        t_365[k] = f_11 * sih_189[k]
                   + f_3 * pc_y[k] * skh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_y, pc_y, pc_z, sii0_255, sii0_257, \
                         sih_168, sih_190, sih_191, sii1_255, sii1_257, skh_273, \
                         skh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * sih_168[k]
                   + f_3 * pc_z[k] * skh_273[k];

        t_367[k] = pb_y[k] * sii0_255[k]
                   + f_12 * sih_190[k]
                   - f_10 * pc_y[k] * sii1_255[k];

        t_368[k] = f_11 * sih_191[k]
                   + f_3 * pc_y[k] * skh_275[k];

        t_369[k] = pb_y[k] * sii0_257[k]
                   - f_10 * pc_y[k] * sii1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_y, pc_y, pc_z, sii0_258, sii0_261, \
                         sih_171, sih_192, sih_194, sii1_258, sii1_261, skh_276, \
                         skh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_y[k] * sii0_258[k]
                   + f_13 * sih_192[k]
                   - f_10 * pc_y[k] * sii1_258[k];

        t_371[k] = f_13 * sih_171[k]
                   + f_3 * pc_z[k] * skh_276[k];

        t_372[k] = f_11 * sih_194[k]
                   + f_3 * pc_y[k] * skh_278[k];

        t_373[k] = pb_y[k] * sii0_261[k]
                   - f_10 * pc_y[k] * sii1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pb_y, pc_y, pc_z, sii0_262, sii0_264, sih_174, \
                         sih_195, sih_197, sii1_262, sii1_264, \
                         skh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_y[k] * sii0_262[k]
                   + f_14 * sih_195[k]
                   - f_10 * pc_y[k] * sii1_262[k];

        t_375[k] = f_13 * sih_174[k]
                   + f_3 * pc_z[k] * skh_279[k];

        t_376[k] = pb_y[k] * sii0_264[k]
                   + f_12 * sih_197[k]
                   - f_10 * pc_y[k] * sii1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_y, pc_x, pc_y, sii0_266, sih_198, \
                         sih_288, sih_289, sii1_266, skh_282, skh_288, \
                         skh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * sih_198[k]
                   + f_3 * pc_y[k] * skh_282[k];

        t_378[k] = pb_y[k] * sii0_266[k]
                   - f_10 * pc_y[k] * sii1_266[k];

        t_379[k] = f_13 * sih_288[k]
                   + f_3 * pc_x[k] * skh_288[k];

        t_380[k] = f_13 * sih_289[k]
                   + f_3 * pc_x[k] * skh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, sih_290, sih_291, sih_292, sih_293, \
                         skh_290, skh_291, skh_292, skh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_13 * sih_290[k]
                   + f_3 * pc_x[k] * skh_290[k];

        t_382[k] = f_13 * sih_291[k]
                   + f_3 * pc_x[k] * skh_291[k];

        t_383[k] = f_13 * sih_292[k]
                   + f_3 * pc_x[k] * skh_292[k];

        t_384[k] = f_13 * sih_293[k]
                   + f_3 * pc_x[k] * skh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, sih_183, sih_204, sih_206, skg0_205, \
                         skg0_207, skg1_205, skg1_207, skh_288, \
                         skh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * sih_204[k]
                   + f_1 * skg0_205[k]
                   - f_2 * skg1_205[k]
                   + f_3 * pc_y[k] * skh_288[k];

        t_386[k] = f_13 * sih_183[k]
                   + f_3 * pc_z[k] * skh_288[k];

        t_387[k] = f_11 * sih_206[k]
                   + f_4 * skg0_207[k]
                   - f_5 * skg1_207[k]
                   + f_3 * pc_y[k] * skh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, sih_207, sih_208, sih_209, skg0_208, \
                         skg0_209, skg1_208, skg1_209, skh_291, skh_292, \
                         skh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * sih_207[k]
                   + f_6 * skg0_208[k]
                   - f_7 * skg1_208[k]
                   + f_3 * pc_y[k] * skh_291[k];

        t_389[k] = f_11 * sih_208[k]
                   + f_8 * skg0_209[k]
                   - f_9 * skg1_209[k]
                   + f_3 * pc_y[k] * skh_292[k];

        t_390[k] = f_11 * sih_209[k]
                   + f_3 * pc_y[k] * skh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_y, pc_x, pc_y, pc_z, sii0_279, \
                         sih_189, sih_294, sii1_279, skg0_210, skg1_210, \
                         skh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pb_y[k] * sii0_279[k]
                   - f_10 * pc_y[k] * sii1_279[k];

        t_392[k] = f_13 * sih_294[k]
                   + f_1 * skg0_210[k]
                   - f_2 * skg1_210[k]
                   + f_3 * pc_x[k] * skh_294[k];

        t_393[k] = f_3 * pc_y[k] * skh_294[k];

        t_394[k] = f_14 * sih_189[k]
                   + f_3 * pc_z[k] * skh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, sih_297, sih_299, skg0_213, \
                         skg0_215, skg1_213, skg1_215, skh_296, skh_297, \
                         skh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_13 * sih_297[k]
                   + f_4 * skg0_213[k]
                   - f_5 * skg1_213[k]
                   + f_3 * pc_x[k] * skh_297[k];

        t_396[k] = f_3 * pc_y[k] * skh_296[k];

        t_397[k] = f_13 * sih_299[k]
                   + f_4 * skg0_215[k]
                   - f_5 * skg1_215[k]
                   + f_3 * pc_x[k] * skh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_x, pc_y, pc_z, sih_192, sih_300, skg0_216, \
                         skg1_216, skh_297, skh_299, skh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * sih_300[k]
                   + f_6 * skg0_216[k]
                   - f_7 * skg1_216[k]
                   + f_3 * pc_x[k] * skh_300[k];

        t_399[k] = f_14 * sih_192[k]
                   + f_3 * pc_z[k] * skh_297[k];

        t_400[k] = f_3 * pc_y[k] * skh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_z, sih_195, sih_303, sih_304, skg0_219, \
                         skg0_220, skg1_219, skg1_220, skh_300, skh_303, \
                         skh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_13 * sih_303[k]
                   + f_6 * skg0_219[k]
                   - f_7 * skg1_219[k]
                   + f_3 * pc_x[k] * skh_303[k];

        t_402[k] = f_13 * sih_304[k]
                   + f_8 * skg0_220[k]
                   - f_9 * skg1_220[k]
                   + f_3 * pc_x[k] * skh_304[k];

        t_403[k] = f_14 * sih_195[k]
                   + f_3 * pc_z[k] * skh_300[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pc_x, pc_y, sih_306, sih_308, skg0_222, \
                         skg0_224, skg1_222, skg1_224, skh_303, skh_306, \
                         skh_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_13 * sih_306[k]
                   + f_8 * skg0_222[k]
                   - f_9 * skg1_222[k]
                   + f_3 * pc_x[k] * skh_306[k];

        t_405[k] = f_3 * pc_y[k] * skh_303[k];

        t_406[k] = f_13 * sih_308[k]
                   + f_8 * skg0_224[k]
                   - f_9 * skg1_224[k]
                   + f_3 * pc_x[k] * skh_308[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pc_x, sih_309, sih_310, sih_311, \
                         sih_312, sih_313, skh_309, skh_310, skh_311, skh_312, \
                         skh_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * sih_309[k]
                   + f_3 * pc_x[k] * skh_309[k];

        t_408[k] = f_13 * sih_310[k]
                   + f_3 * pc_x[k] * skh_310[k];

        t_409[k] = f_13 * sih_311[k]
                   + f_3 * pc_x[k] * skh_311[k];

        t_410[k] = f_13 * sih_312[k]
                   + f_3 * pc_x[k] * skh_312[k];

        t_411[k] = f_13 * sih_313[k]
                   + f_3 * pc_x[k] * skh_313[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pc_x, pc_y, pc_z, sih_204, sih_314, \
                         skg0_220, skg0_222, skg1_220, skg1_222, skh_309, skh_311, \
                         skh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_13 * sih_314[k]
                   + f_3 * pc_x[k] * skh_314[k];

        t_413[k] = f_1 * skg0_220[k]
                   - f_2 * skg1_220[k]
                   + f_3 * pc_y[k] * skh_309[k];

        t_414[k] = f_14 * sih_204[k]
                   + f_3 * pc_z[k] * skh_309[k];

        t_415[k] = f_4 * skg0_222[k]
                   - f_5 * skg1_222[k]
                   + f_3 * pc_y[k] * skh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, sih_209, skg0_223, skg0_224, \
                         skg1_223, skg1_224, skh_312, skh_313, \
                         skh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * skg0_223[k]
                   - f_7 * skg1_223[k]
                   + f_3 * pc_y[k] * skh_312[k];

        t_417[k] = f_8 * skg0_224[k]
                   - f_9 * skg1_224[k]
                   + f_3 * pc_y[k] * skh_313[k];

        t_418[k] = f_3 * pc_y[k] * skh_314[k];

        t_419[k] = f_14 * sih_209[k]
                   + f_1 * skg0_224[k]
                   - f_2 * skg1_224[k]
                   + f_3 * pc_z[k] * skh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, sih_210, sih_315, \
                         sih_318, skg0_225, skg0_228, skg1_225, skg1_228, skh_315, \
                         skh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_12 * sih_315[k]
                   + f_1 * skg0_225[k]
                   - f_2 * skg1_225[k]
                   + f_3 * pc_x[k] * skh_315[k];

        t_421[k] = f_16 * sih_210[k]
                   + f_3 * pc_y[k] * skh_315[k];

        t_422[k] = f_3 * pc_z[k] * skh_315[k];

        t_423[k] = f_12 * sih_318[k]
                   + f_4 * skg0_228[k]
                   - f_5 * skg1_228[k]
                   + f_3 * pc_x[k] * skh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, sih_212, sih_320, sih_321, skg0_230, \
                         skg0_231, skg1_230, skg1_231, skh_317, skh_320, \
                         skh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_16 * sih_212[k]
                   + f_3 * pc_y[k] * skh_317[k];

        t_425[k] = f_12 * sih_320[k]
                   + f_4 * skg0_230[k]
                   - f_5 * skg1_230[k]
                   + f_3 * pc_x[k] * skh_320[k];

        t_426[k] = f_12 * sih_321[k]
                   + f_6 * skg0_231[k]
                   - f_7 * skg1_231[k]
                   + f_3 * pc_x[k] * skh_321[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_x, pc_y, pc_z, sih_215, sih_324, skg0_234, \
                         skg1_234, skh_318, skh_320, skh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_3 * pc_z[k] * skh_318[k];

        t_428[k] = f_16 * sih_215[k]
                   + f_3 * pc_y[k] * skh_320[k];

        t_429[k] = f_12 * sih_324[k]
                   + f_6 * skg0_234[k]
                   - f_7 * skg1_234[k]
                   + f_3 * pc_x[k] * skh_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_z, sih_325, sih_327, skg0_235, \
                         skg0_237, skg1_235, skg1_237, skh_321, skh_325, \
                         skh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_12 * sih_325[k]
                   + f_8 * skg0_235[k]
                   - f_9 * skg1_235[k]
                   + f_3 * pc_x[k] * skh_325[k];

        t_431[k] = f_3 * pc_z[k] * skh_321[k];

        t_432[k] = f_12 * sih_327[k]
                   + f_8 * skg0_237[k]
                   - f_9 * skg1_237[k]
                   + f_3 * pc_x[k] * skh_327[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, sih_219, sih_329, sih_330, \
                         sih_331, skg0_239, skg1_239, skh_324, skh_329, skh_330, \
                         skh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_16 * sih_219[k]
                   + f_3 * pc_y[k] * skh_324[k];

        t_434[k] = f_12 * sih_329[k]
                   + f_8 * skg0_239[k]
                   - f_9 * skg1_239[k]
                   + f_3 * pc_x[k] * skh_329[k];

        t_435[k] = f_12 * sih_330[k]
                   + f_3 * pc_x[k] * skh_330[k];

        t_436[k] = f_12 * sih_331[k]
                   + f_3 * pc_x[k] * skh_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, sih_332, sih_333, sih_334, sih_335, \
                         skh_332, skh_333, skh_334, skh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_12 * sih_332[k]
                   + f_3 * pc_x[k] * skh_332[k];

        t_438[k] = f_12 * sih_333[k]
                   + f_3 * pc_x[k] * skh_333[k];

        t_439[k] = f_12 * sih_334[k]
                   + f_3 * pc_x[k] * skh_334[k];

        t_440[k] = f_12 * sih_335[k]
                   + f_3 * pc_x[k] * skh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, pc_z, sih_225, sih_227, skg0_235, \
                         skg0_237, skg1_235, skg1_237, skh_330, \
                         skh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_16 * sih_225[k]
                   + f_1 * skg0_235[k]
                   - f_2 * skg1_235[k]
                   + f_3 * pc_y[k] * skh_330[k];

        t_442[k] = f_3 * pc_z[k] * skh_330[k];

        t_443[k] = f_16 * sih_227[k]
                   + f_4 * skg0_237[k]
                   - f_5 * skg1_237[k]
                   + f_3 * pc_y[k] * skh_332[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, sih_228, sih_229, sih_230, \
                         skg0_238, skg0_239, skg1_238, skg1_239, skh_333, skh_334, \
                         skh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_16 * sih_228[k]
                   + f_6 * skg0_238[k]
                   - f_7 * skg1_238[k]
                   + f_3 * pc_y[k] * skh_333[k];

        t_445[k] = f_16 * sih_229[k]
                   + f_8 * skg0_239[k]
                   - f_9 * skg1_239[k]
                   + f_3 * pc_y[k] * skh_334[k];

        t_446[k] = f_16 * sih_230[k]
                   + f_3 * pc_y[k] * skh_335[k];

        t_447[k] = f_1 * skg0_239[k]
                   - f_2 * skg1_239[k]
                   + f_3 * pc_z[k] * skh_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_z, pc_y, pc_z, sii0_280, sii0_283, \
                         sih_210, sih_231, sii1_280, sii1_283, \
                         skh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_z[k] * sii0_280[k]
                   - f_10 * pc_z[k] * sii1_280[k];

        t_449[k] = f_14 * sih_231[k]
                   + f_3 * pc_y[k] * skh_336[k];

        t_450[k] = f_11 * sih_210[k]
                   + f_3 * pc_z[k] * skh_336[k];

        t_451[k] = pb_z[k] * sii0_283[k]
                   - f_10 * pc_z[k] * sii1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_z, pc_x, pc_y, pc_z, sii0_286, sih_233, \
                         sih_341, sii1_286, skg0_245, skg1_245, skh_338, \
                         skh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * sih_233[k]
                   + f_3 * pc_y[k] * skh_338[k];

        t_453[k] = f_12 * sih_341[k]
                   + f_4 * skg0_245[k]
                   - f_5 * skg1_245[k]
                   + f_3 * pc_x[k] * skh_341[k];

        t_454[k] = pb_z[k] * sii0_286[k]
                   - f_10 * pc_z[k] * sii1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_y, pc_z, sih_213, sih_236, sih_345, \
                         skg0_249, skg1_249, skh_339, skh_341, \
                         skh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * sih_213[k]
                   + f_3 * pc_z[k] * skh_339[k];

        t_456[k] = f_14 * sih_236[k]
                   + f_3 * pc_y[k] * skh_341[k];

        t_457[k] = f_12 * sih_345[k]
                   + f_6 * skg0_249[k]
                   - f_7 * skg1_249[k]
                   + f_3 * pc_x[k] * skh_345[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pb_z, pc_y, pc_z, sii0_290, sii0_292, \
                         sih_216, sih_217, sih_240, sii1_290, sii1_292, skh_342, \
                         skh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pb_z[k] * sii0_290[k]
                   - f_10 * pc_z[k] * sii1_290[k];

        t_459[k] = f_11 * sih_216[k]
                   + f_3 * pc_z[k] * skh_342[k];

        t_460[k] = pb_z[k] * sii0_292[k]
                   + f_12 * sih_217[k]
                   - f_10 * pc_z[k] * sii1_292[k];

        t_461[k] = f_14 * sih_240[k]
                   + f_3 * pc_y[k] * skh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, sih_350, sih_351, sih_352, sih_353, \
                         skg0_254, skg1_254, skh_350, skh_351, skh_352, \
                         skh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_12 * sih_350[k]
                   + f_8 * skg0_254[k]
                   - f_9 * skg1_254[k]
                   + f_3 * pc_x[k] * skh_350[k];

        t_463[k] = f_12 * sih_351[k]
                   + f_3 * pc_x[k] * skh_351[k];

        t_464[k] = f_12 * sih_352[k]
                   + f_3 * pc_x[k] * skh_352[k];

        t_465[k] = f_12 * sih_353[k]
                   + f_3 * pc_x[k] * skh_353[k];
    }
}

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.5 / q;

    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_301 = buffer.data(sii0 + 301);
    const auto *sii0_392 = buffer.data(sii0 + 392);
    const auto *sii0_395 = buffer.data(sii0 + 395);
    const auto *sii0_397 = buffer.data(sii0 + 397);
    const auto *sii0_398 = buffer.data(sii0 + 398);
    const auto *sii0_401 = buffer.data(sii0 + 401);
    const auto *sii0_402 = buffer.data(sii0 + 402);
    const auto *sii0_404 = buffer.data(sii0 + 404);
    const auto *sii0_406 = buffer.data(sii0 + 406);
    const auto *sii0_419 = buffer.data(sii0 + 419);

    const auto *sih_225 = buffer.data(sih + 225);
    const auto *sih_230 = buffer.data(sih + 230);
    const auto *sih_231 = buffer.data(sih + 231);
    const auto *sih_234 = buffer.data(sih + 234);
    const auto *sih_237 = buffer.data(sih + 237);
    const auto *sih_246 = buffer.data(sih + 246);
    const auto *sih_248 = buffer.data(sih + 248);
    const auto *sih_249 = buffer.data(sih + 249);
    const auto *sih_250 = buffer.data(sih + 250);
    const auto *sih_251 = buffer.data(sih + 251);
    const auto *sih_252 = buffer.data(sih + 252);
    const auto *sih_254 = buffer.data(sih + 254);
    const auto *sih_255 = buffer.data(sih + 255);
    const auto *sih_257 = buffer.data(sih + 257);
    const auto *sih_258 = buffer.data(sih + 258);
    const auto *sih_261 = buffer.data(sih + 261);
    const auto *sih_267 = buffer.data(sih + 267);
    const auto *sih_269 = buffer.data(sih + 269);
    const auto *sih_270 = buffer.data(sih + 270);
    const auto *sih_271 = buffer.data(sih + 271);
    const auto *sih_272 = buffer.data(sih + 272);
    const auto *sih_273 = buffer.data(sih + 273);
    const auto *sih_275 = buffer.data(sih + 275);
    const auto *sih_276 = buffer.data(sih + 276);
    const auto *sih_278 = buffer.data(sih + 278);
    const auto *sih_279 = buffer.data(sih + 279);
    const auto *sih_282 = buffer.data(sih + 282);
    const auto *sih_288 = buffer.data(sih + 288);
    const auto *sih_290 = buffer.data(sih + 290);
    const auto *sih_291 = buffer.data(sih + 291);
    const auto *sih_292 = buffer.data(sih + 292);
    const auto *sih_293 = buffer.data(sih + 293);
    const auto *sih_294 = buffer.data(sih + 294);
    const auto *sih_295 = buffer.data(sih + 295);
    const auto *sih_296 = buffer.data(sih + 296);
    const auto *sih_297 = buffer.data(sih + 297);
    const auto *sih_299 = buffer.data(sih + 299);
    const auto *sih_300 = buffer.data(sih + 300);
    const auto *sih_302 = buffer.data(sih + 302);
    const auto *sih_303 = buffer.data(sih + 303);
    const auto *sih_309 = buffer.data(sih + 309);
    const auto *sih_311 = buffer.data(sih + 311);
    const auto *sih_312 = buffer.data(sih + 312);
    const auto *sih_313 = buffer.data(sih + 313);
    const auto *sih_314 = buffer.data(sih + 314);
    const auto *sih_354 = buffer.data(sih + 354);
    const auto *sih_355 = buffer.data(sih + 355);
    const auto *sih_356 = buffer.data(sih + 356);
    const auto *sih_357 = buffer.data(sih + 357);
    const auto *sih_360 = buffer.data(sih + 360);
    const auto *sih_362 = buffer.data(sih + 362);
    const auto *sih_363 = buffer.data(sih + 363);
    const auto *sih_366 = buffer.data(sih + 366);
    const auto *sih_367 = buffer.data(sih + 367);
    const auto *sih_369 = buffer.data(sih + 369);
    const auto *sih_371 = buffer.data(sih + 371);
    const auto *sih_372 = buffer.data(sih + 372);
    const auto *sih_373 = buffer.data(sih + 373);
    const auto *sih_374 = buffer.data(sih + 374);
    const auto *sih_375 = buffer.data(sih + 375);
    const auto *sih_376 = buffer.data(sih + 376);
    const auto *sih_377 = buffer.data(sih + 377);
    const auto *sih_378 = buffer.data(sih + 378);
    const auto *sih_381 = buffer.data(sih + 381);
    const auto *sih_383 = buffer.data(sih + 383);
    const auto *sih_384 = buffer.data(sih + 384);
    const auto *sih_387 = buffer.data(sih + 387);
    const auto *sih_388 = buffer.data(sih + 388);
    const auto *sih_390 = buffer.data(sih + 390);
    const auto *sih_392 = buffer.data(sih + 392);
    const auto *sih_393 = buffer.data(sih + 393);
    const auto *sih_394 = buffer.data(sih + 394);
    const auto *sih_395 = buffer.data(sih + 395);
    const auto *sih_396 = buffer.data(sih + 396);
    const auto *sih_397 = buffer.data(sih + 397);
    const auto *sih_398 = buffer.data(sih + 398);
    const auto *sih_414 = buffer.data(sih + 414);
    const auto *sih_415 = buffer.data(sih + 415);
    const auto *sih_416 = buffer.data(sih + 416);
    const auto *sih_417 = buffer.data(sih + 417);
    const auto *sih_418 = buffer.data(sih + 418);
    const auto *sih_419 = buffer.data(sih + 419);
    const auto *sih_420 = buffer.data(sih + 420);
    const auto *sih_423 = buffer.data(sih + 423);
    const auto *sih_425 = buffer.data(sih + 425);
    const auto *sih_426 = buffer.data(sih + 426);
    const auto *sih_429 = buffer.data(sih + 429);
    const auto *sih_430 = buffer.data(sih + 430);
    const auto *sih_432 = buffer.data(sih + 432);
    const auto *sih_434 = buffer.data(sih + 434);

    const auto *sii1_301 = buffer.data(sii1 + 301);
    const auto *sii1_392 = buffer.data(sii1 + 392);
    const auto *sii1_395 = buffer.data(sii1 + 395);
    const auto *sii1_397 = buffer.data(sii1 + 397);
    const auto *sii1_398 = buffer.data(sii1 + 398);
    const auto *sii1_401 = buffer.data(sii1 + 401);
    const auto *sii1_402 = buffer.data(sii1 + 402);
    const auto *sii1_404 = buffer.data(sii1 + 404);
    const auto *sii1_406 = buffer.data(sii1 + 406);
    const auto *sii1_419 = buffer.data(sii1 + 419);

    const auto *skg0_252 = buffer.data(skg0 + 252);
    const auto *skg0_253 = buffer.data(skg0 + 253);
    const auto *skg0_254 = buffer.data(skg0 + 254);
    const auto *skg0_255 = buffer.data(skg0 + 255);
    const auto *skg0_258 = buffer.data(skg0 + 258);
    const auto *skg0_260 = buffer.data(skg0 + 260);
    const auto *skg0_261 = buffer.data(skg0 + 261);
    const auto *skg0_264 = buffer.data(skg0 + 264);
    const auto *skg0_265 = buffer.data(skg0 + 265);
    const auto *skg0_267 = buffer.data(skg0 + 267);
    const auto *skg0_268 = buffer.data(skg0 + 268);
    const auto *skg0_269 = buffer.data(skg0 + 269);
    const auto *skg0_270 = buffer.data(skg0 + 270);
    const auto *skg0_273 = buffer.data(skg0 + 273);
    const auto *skg0_275 = buffer.data(skg0 + 275);
    const auto *skg0_276 = buffer.data(skg0 + 276);
    const auto *skg0_279 = buffer.data(skg0 + 279);
    const auto *skg0_280 = buffer.data(skg0 + 280);
    const auto *skg0_282 = buffer.data(skg0 + 282);
    const auto *skg0_283 = buffer.data(skg0 + 283);
    const auto *skg0_284 = buffer.data(skg0 + 284);
    const auto *skg0_295 = buffer.data(skg0 + 295);
    const auto *skg0_297 = buffer.data(skg0 + 297);
    const auto *skg0_298 = buffer.data(skg0 + 298);
    const auto *skg0_299 = buffer.data(skg0 + 299);
    const auto *skg0_300 = buffer.data(skg0 + 300);
    const auto *skg0_303 = buffer.data(skg0 + 303);
    const auto *skg0_305 = buffer.data(skg0 + 305);
    const auto *skg0_306 = buffer.data(skg0 + 306);
    const auto *skg0_309 = buffer.data(skg0 + 309);
    const auto *skg0_310 = buffer.data(skg0 + 310);
    const auto *skg0_312 = buffer.data(skg0 + 312);
    const auto *skg0_314 = buffer.data(skg0 + 314);

    const auto *skg1_252 = buffer.data(skg1 + 252);
    const auto *skg1_253 = buffer.data(skg1 + 253);
    const auto *skg1_254 = buffer.data(skg1 + 254);
    const auto *skg1_255 = buffer.data(skg1 + 255);
    const auto *skg1_258 = buffer.data(skg1 + 258);
    const auto *skg1_260 = buffer.data(skg1 + 260);
    const auto *skg1_261 = buffer.data(skg1 + 261);
    const auto *skg1_264 = buffer.data(skg1 + 264);
    const auto *skg1_265 = buffer.data(skg1 + 265);
    const auto *skg1_267 = buffer.data(skg1 + 267);
    const auto *skg1_268 = buffer.data(skg1 + 268);
    const auto *skg1_269 = buffer.data(skg1 + 269);
    const auto *skg1_270 = buffer.data(skg1 + 270);
    const auto *skg1_273 = buffer.data(skg1 + 273);
    const auto *skg1_275 = buffer.data(skg1 + 275);
    const auto *skg1_276 = buffer.data(skg1 + 276);
    const auto *skg1_279 = buffer.data(skg1 + 279);
    const auto *skg1_280 = buffer.data(skg1 + 280);
    const auto *skg1_282 = buffer.data(skg1 + 282);
    const auto *skg1_283 = buffer.data(skg1 + 283);
    const auto *skg1_284 = buffer.data(skg1 + 284);
    const auto *skg1_295 = buffer.data(skg1 + 295);
    const auto *skg1_297 = buffer.data(skg1 + 297);
    const auto *skg1_298 = buffer.data(skg1 + 298);
    const auto *skg1_299 = buffer.data(skg1 + 299);
    const auto *skg1_300 = buffer.data(skg1 + 300);
    const auto *skg1_303 = buffer.data(skg1 + 303);
    const auto *skg1_305 = buffer.data(skg1 + 305);
    const auto *skg1_306 = buffer.data(skg1 + 306);
    const auto *skg1_309 = buffer.data(skg1 + 309);
    const auto *skg1_310 = buffer.data(skg1 + 310);
    const auto *skg1_312 = buffer.data(skg1 + 312);
    const auto *skg1_314 = buffer.data(skg1 + 314);

    const auto *skh_351 = buffer.data(skh + 351);
    const auto *skh_353 = buffer.data(skh + 353);
    const auto *skh_354 = buffer.data(skh + 354);
    const auto *skh_355 = buffer.data(skh + 355);
    const auto *skh_356 = buffer.data(skh + 356);
    const auto *skh_357 = buffer.data(skh + 357);
    const auto *skh_359 = buffer.data(skh + 359);
    const auto *skh_360 = buffer.data(skh + 360);
    const auto *skh_362 = buffer.data(skh + 362);
    const auto *skh_363 = buffer.data(skh + 363);
    const auto *skh_366 = buffer.data(skh + 366);
    const auto *skh_367 = buffer.data(skh + 367);
    const auto *skh_369 = buffer.data(skh + 369);
    const auto *skh_371 = buffer.data(skh + 371);
    const auto *skh_372 = buffer.data(skh + 372);
    const auto *skh_373 = buffer.data(skh + 373);
    const auto *skh_374 = buffer.data(skh + 374);
    const auto *skh_375 = buffer.data(skh + 375);
    const auto *skh_376 = buffer.data(skh + 376);
    const auto *skh_377 = buffer.data(skh + 377);
    const auto *skh_378 = buffer.data(skh + 378);
    const auto *skh_380 = buffer.data(skh + 380);
    const auto *skh_381 = buffer.data(skh + 381);
    const auto *skh_383 = buffer.data(skh + 383);
    const auto *skh_384 = buffer.data(skh + 384);
    const auto *skh_387 = buffer.data(skh + 387);
    const auto *skh_388 = buffer.data(skh + 388);
    const auto *skh_390 = buffer.data(skh + 390);
    const auto *skh_392 = buffer.data(skh + 392);
    const auto *skh_393 = buffer.data(skh + 393);
    const auto *skh_394 = buffer.data(skh + 394);
    const auto *skh_395 = buffer.data(skh + 395);
    const auto *skh_396 = buffer.data(skh + 396);
    const auto *skh_397 = buffer.data(skh + 397);
    const auto *skh_398 = buffer.data(skh + 398);
    const auto *skh_399 = buffer.data(skh + 399);
    const auto *skh_401 = buffer.data(skh + 401);
    const auto *skh_402 = buffer.data(skh + 402);
    const auto *skh_404 = buffer.data(skh + 404);
    const auto *skh_405 = buffer.data(skh + 405);
    const auto *skh_408 = buffer.data(skh + 408);
    const auto *skh_414 = buffer.data(skh + 414);
    const auto *skh_415 = buffer.data(skh + 415);
    const auto *skh_416 = buffer.data(skh + 416);
    const auto *skh_417 = buffer.data(skh + 417);
    const auto *skh_418 = buffer.data(skh + 418);
    const auto *skh_419 = buffer.data(skh + 419);
    const auto *skh_420 = buffer.data(skh + 420);
    const auto *skh_422 = buffer.data(skh + 422);
    const auto *skh_423 = buffer.data(skh + 423);
    const auto *skh_425 = buffer.data(skh + 425);
    const auto *skh_426 = buffer.data(skh + 426);
    const auto *skh_429 = buffer.data(skh + 429);
    const auto *skh_430 = buffer.data(skh + 430);
    const auto *skh_432 = buffer.data(skh + 432);
    const auto *skh_434 = buffer.data(skh + 434);

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_z, pc_x, pc_z, sii0_301, sih_354, \
                         sih_355, sih_356, sii1_301, skh_354, skh_355, \
                         skh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_12 * sih_354[k]
                   + f_3 * pc_x[k] * skh_354[k];

        t_467[k] = f_12 * sih_355[k]
                   + f_3 * pc_x[k] * skh_355[k];

        t_468[k] = f_12 * sih_356[k]
                   + f_3 * pc_x[k] * skh_356[k];

        t_469[k] = pb_z[k] * sii0_301[k]
                   - f_10 * pc_z[k] * sii1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, sih_225, sih_248, sih_249, skg0_252, \
                         skg0_253, skg1_252, skg1_253, skh_351, skh_353, \
                         skh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * sih_225[k]
                   + f_3 * pc_z[k] * skh_351[k];

        t_471[k] = f_14 * sih_248[k]
                   + f_4 * skg0_252[k]
                   - f_5 * skg1_252[k]
                   + f_3 * pc_y[k] * skh_353[k];

        t_472[k] = f_14 * sih_249[k]
                   + f_6 * skg0_253[k]
                   - f_7 * skg1_253[k]
                   + f_3 * pc_y[k] * skh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, sih_230, sih_250, sih_251, skg0_254, \
                         skg1_254, skh_355, skh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * sih_250[k]
                   + f_8 * skg0_254[k]
                   - f_9 * skg1_254[k]
                   + f_3 * pc_y[k] * skh_355[k];

        t_474[k] = f_14 * sih_251[k]
                   + f_3 * pc_y[k] * skh_356[k];

        t_475[k] = f_11 * sih_230[k]
                   + f_1 * skg0_254[k]
                   - f_2 * skg1_254[k]
                   + f_3 * pc_z[k] * skh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, sih_231, sih_252, sih_357, \
                         skg0_255, skg1_255, skh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_12 * sih_357[k]
                   + f_1 * skg0_255[k]
                   - f_2 * skg1_255[k]
                   + f_3 * pc_x[k] * skh_357[k];

        t_477[k] = f_13 * sih_252[k]
                   + f_3 * pc_y[k] * skh_357[k];

        t_478[k] = f_12 * sih_231[k]
                   + f_3 * pc_z[k] * skh_357[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, sih_254, sih_360, sih_362, skg0_258, \
                         skg0_260, skg1_258, skg1_260, skh_359, skh_360, \
                         skh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_12 * sih_360[k]
                   + f_4 * skg0_258[k]
                   - f_5 * skg1_258[k]
                   + f_3 * pc_x[k] * skh_360[k];

        t_480[k] = f_13 * sih_254[k]
                   + f_3 * pc_y[k] * skh_359[k];

        t_481[k] = f_12 * sih_362[k]
                   + f_4 * skg0_260[k]
                   - f_5 * skg1_260[k]
                   + f_3 * pc_x[k] * skh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, sih_234, sih_257, sih_363, \
                         skg0_261, skg1_261, skh_360, skh_362, \
                         skh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_12 * sih_363[k]
                   + f_6 * skg0_261[k]
                   - f_7 * skg1_261[k]
                   + f_3 * pc_x[k] * skh_363[k];

        t_483[k] = f_12 * sih_234[k]
                   + f_3 * pc_z[k] * skh_360[k];

        t_484[k] = f_13 * sih_257[k]
                   + f_3 * pc_y[k] * skh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, sih_237, sih_366, sih_367, skg0_264, \
                         skg0_265, skg1_264, skg1_265, skh_363, skh_366, \
                         skh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_12 * sih_366[k]
                   + f_6 * skg0_264[k]
                   - f_7 * skg1_264[k]
                   + f_3 * pc_x[k] * skh_366[k];

        t_486[k] = f_12 * sih_367[k]
                   + f_8 * skg0_265[k]
                   - f_9 * skg1_265[k]
                   + f_3 * pc_x[k] * skh_367[k];

        t_487[k] = f_12 * sih_237[k]
                   + f_3 * pc_z[k] * skh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, sih_261, sih_369, sih_371, skg0_267, \
                         skg0_269, skg1_267, skg1_269, skh_366, skh_369, \
                         skh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_12 * sih_369[k]
                   + f_8 * skg0_267[k]
                   - f_9 * skg1_267[k]
                   + f_3 * pc_x[k] * skh_369[k];

        t_489[k] = f_13 * sih_261[k]
                   + f_3 * pc_y[k] * skh_366[k];

        t_490[k] = f_12 * sih_371[k]
                   + f_8 * skg0_269[k]
                   - f_9 * skg1_269[k]
                   + f_3 * pc_x[k] * skh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, sih_372, sih_373, sih_374, \
                         sih_375, sih_376, skh_372, skh_373, skh_374, skh_375, \
                         skh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_12 * sih_372[k]
                   + f_3 * pc_x[k] * skh_372[k];

        t_492[k] = f_12 * sih_373[k]
                   + f_3 * pc_x[k] * skh_373[k];

        t_493[k] = f_12 * sih_374[k]
                   + f_3 * pc_x[k] * skh_374[k];

        t_494[k] = f_12 * sih_375[k]
                   + f_3 * pc_x[k] * skh_375[k];

        t_495[k] = f_12 * sih_376[k]
                   + f_3 * pc_x[k] * skh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, sih_246, sih_267, sih_377, \
                         skg0_265, skg1_265, skh_372, skh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_12 * sih_377[k]
                   + f_3 * pc_x[k] * skh_377[k];

        t_497[k] = f_13 * sih_267[k]
                   + f_1 * skg0_265[k]
                   - f_2 * skg1_265[k]
                   + f_3 * pc_y[k] * skh_372[k];

        t_498[k] = f_12 * sih_246[k]
                   + f_3 * pc_z[k] * skh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, sih_269, sih_270, sih_271, skg0_267, \
                         skg0_268, skg0_269, skg1_267, skg1_268, skg1_269, skh_374, skh_375, \
                         skh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * sih_269[k]
                   + f_4 * skg0_267[k]
                   - f_5 * skg1_267[k]
                   + f_3 * pc_y[k] * skh_374[k];

        t_500[k] = f_13 * sih_270[k]
                   + f_6 * skg0_268[k]
                   - f_7 * skg1_268[k]
                   + f_3 * pc_y[k] * skh_375[k];

        t_501[k] = f_13 * sih_271[k]
                   + f_8 * skg0_269[k]
                   - f_9 * skg1_269[k]
                   + f_3 * pc_y[k] * skh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, sih_251, sih_272, sih_378, \
                         skg0_269, skg0_270, skg1_269, skg1_270, skh_377, \
                         skh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * sih_272[k]
                   + f_3 * pc_y[k] * skh_377[k];

        t_503[k] = f_12 * sih_251[k]
                   + f_1 * skg0_269[k]
                   - f_2 * skg1_269[k]
                   + f_3 * pc_z[k] * skh_377[k];

        t_504[k] = f_12 * sih_378[k]
                   + f_1 * skg0_270[k]
                   - f_2 * skg1_270[k]
                   + f_3 * pc_x[k] * skh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, sih_252, sih_273, \
                         sih_275, sih_381, skg0_273, skg1_273, skh_378, skh_380, \
                         skh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * sih_273[k]
                   + f_3 * pc_y[k] * skh_378[k];

        t_506[k] = f_13 * sih_252[k]
                   + f_3 * pc_z[k] * skh_378[k];

        t_507[k] = f_12 * sih_381[k]
                   + f_4 * skg0_273[k]
                   - f_5 * skg1_273[k]
                   + f_3 * pc_x[k] * skh_381[k];

        t_508[k] = f_12 * sih_275[k]
                   + f_3 * pc_y[k] * skh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, sih_255, sih_383, sih_384, skg0_275, \
                         skg0_276, skg1_275, skg1_276, skh_381, skh_383, \
                         skh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_12 * sih_383[k]
                   + f_4 * skg0_275[k]
                   - f_5 * skg1_275[k]
                   + f_3 * pc_x[k] * skh_383[k];

        t_510[k] = f_12 * sih_384[k]
                   + f_6 * skg0_276[k]
                   - f_7 * skg1_276[k]
                   + f_3 * pc_x[k] * skh_384[k];

        t_511[k] = f_13 * sih_255[k]
                   + f_3 * pc_z[k] * skh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, sih_278, sih_387, sih_388, skg0_279, \
                         skg0_280, skg1_279, skg1_280, skh_383, skh_387, \
                         skh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * sih_278[k]
                   + f_3 * pc_y[k] * skh_383[k];

        t_513[k] = f_12 * sih_387[k]
                   + f_6 * skg0_279[k]
                   - f_7 * skg1_279[k]
                   + f_3 * pc_x[k] * skh_387[k];

        t_514[k] = f_12 * sih_388[k]
                   + f_8 * skg0_280[k]
                   - f_9 * skg1_280[k]
                   + f_3 * pc_x[k] * skh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, sih_258, sih_282, sih_390, \
                         skg0_282, skg1_282, skh_384, skh_387, \
                         skh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * sih_258[k]
                   + f_3 * pc_z[k] * skh_384[k];

        t_516[k] = f_12 * sih_390[k]
                   + f_8 * skg0_282[k]
                   - f_9 * skg1_282[k]
                   + f_3 * pc_x[k] * skh_390[k];

        t_517[k] = f_12 * sih_282[k]
                   + f_3 * pc_y[k] * skh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, sih_392, sih_393, sih_394, sih_395, \
                         skg0_284, skg1_284, skh_392, skh_393, skh_394, \
                         skh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_12 * sih_392[k]
                   + f_8 * skg0_284[k]
                   - f_9 * skg1_284[k]
                   + f_3 * pc_x[k] * skh_392[k];

        t_519[k] = f_12 * sih_393[k]
                   + f_3 * pc_x[k] * skh_393[k];

        t_520[k] = f_12 * sih_394[k]
                   + f_3 * pc_x[k] * skh_394[k];

        t_521[k] = f_12 * sih_395[k]
                   + f_3 * pc_x[k] * skh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, sih_288, sih_396, sih_397, \
                         sih_398, skg0_280, skg1_280, skh_393, skh_396, skh_397, \
                         skh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_12 * sih_396[k]
                   + f_3 * pc_x[k] * skh_396[k];

        t_523[k] = f_12 * sih_397[k]
                   + f_3 * pc_x[k] * skh_397[k];

        t_524[k] = f_12 * sih_398[k]
                   + f_3 * pc_x[k] * skh_398[k];

        t_525[k] = f_12 * sih_288[k]
                   + f_1 * skg0_280[k]
                   - f_2 * skg1_280[k]
                   + f_3 * pc_y[k] * skh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, sih_267, sih_290, sih_291, skg0_282, \
                         skg0_283, skg1_282, skg1_283, skh_393, skh_395, \
                         skh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * sih_267[k]
                   + f_3 * pc_z[k] * skh_393[k];

        t_527[k] = f_12 * sih_290[k]
                   + f_4 * skg0_282[k]
                   - f_5 * skg1_282[k]
                   + f_3 * pc_y[k] * skh_395[k];

        t_528[k] = f_12 * sih_291[k]
                   + f_6 * skg0_283[k]
                   - f_7 * skg1_283[k]
                   + f_3 * pc_y[k] * skh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pb_y, pc_y, pc_z, sii0_392, sih_272, \
                         sih_292, sih_293, sii1_392, skg0_284, skg1_284, skh_397, \
                         skh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * sih_292[k]
                   + f_8 * skg0_284[k]
                   - f_9 * skg1_284[k]
                   + f_3 * pc_y[k] * skh_397[k];

        t_530[k] = f_12 * sih_293[k]
                   + f_3 * pc_y[k] * skh_398[k];

        t_531[k] = f_13 * sih_272[k]
                   + f_1 * skg0_284[k]
                   - f_2 * skg1_284[k]
                   + f_3 * pc_z[k] * skh_398[k];

        t_532[k] = pb_y[k] * sii0_392[k]
                   - f_10 * pc_y[k] * sii1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pb_y, pc_y, pc_z, sii0_395, sih_273, \
                         sih_294, sih_295, sih_296, sii1_395, skh_399, \
                         skh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * sih_294[k]
                   + f_3 * pc_y[k] * skh_399[k];

        t_534[k] = f_14 * sih_273[k]
                   + f_3 * pc_z[k] * skh_399[k];

        t_535[k] = pb_y[k] * sii0_395[k]
                   + f_12 * sih_295[k]
                   - f_10 * pc_y[k] * sii1_395[k];

        t_536[k] = f_11 * sih_296[k]
                   + f_3 * pc_y[k] * skh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pb_y, pc_y, pc_z, sii0_397, sii0_398, \
                         sih_276, sih_297, sih_299, sii1_397, sii1_398, skh_402, \
                         skh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pb_y[k] * sii0_397[k]
                   - f_10 * pc_y[k] * sii1_397[k];

        t_538[k] = pb_y[k] * sii0_398[k]
                   + f_13 * sih_297[k]
                   - f_10 * pc_y[k] * sii1_398[k];

        t_539[k] = f_14 * sih_276[k]
                   + f_3 * pc_z[k] * skh_402[k];

        t_540[k] = f_11 * sih_299[k]
                   + f_3 * pc_y[k] * skh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pb_y, pc_y, pc_z, sii0_401, sii0_402, sih_279, \
                         sih_300, sii1_401, sii1_402, skh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pb_y[k] * sii0_401[k]
                   - f_10 * pc_y[k] * sii1_401[k];

        t_542[k] = pb_y[k] * sii0_402[k]
                   + f_14 * sih_300[k]
                   - f_10 * pc_y[k] * sii1_402[k];

        t_543[k] = f_14 * sih_279[k]
                   + f_3 * pc_z[k] * skh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_x, pc_y, sii0_404, sii0_406, \
                         sih_302, sih_303, sih_414, sii1_404, sii1_406, skh_408, \
                         skh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pb_y[k] * sii0_404[k]
                   + f_12 * sih_302[k]
                   - f_10 * pc_y[k] * sii1_404[k];

        t_545[k] = f_11 * sih_303[k]
                   + f_3 * pc_y[k] * skh_408[k];

        t_546[k] = pb_y[k] * sii0_406[k]
                   - f_10 * pc_y[k] * sii1_406[k];

        t_547[k] = f_12 * sih_414[k]
                   + f_3 * pc_x[k] * skh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, sih_415, sih_416, sih_417, \
                         sih_418, sih_419, skh_415, skh_416, skh_417, skh_418, \
                         skh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_12 * sih_415[k]
                   + f_3 * pc_x[k] * skh_415[k];

        t_549[k] = f_12 * sih_416[k]
                   + f_3 * pc_x[k] * skh_416[k];

        t_550[k] = f_12 * sih_417[k]
                   + f_3 * pc_x[k] * skh_417[k];

        t_551[k] = f_12 * sih_418[k]
                   + f_3 * pc_x[k] * skh_418[k];

        t_552[k] = f_12 * sih_419[k]
                   + f_3 * pc_x[k] * skh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, sih_288, sih_309, sih_311, skg0_295, \
                         skg0_297, skg1_295, skg1_297, skh_414, \
                         skh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * sih_309[k]
                   + f_1 * skg0_295[k]
                   - f_2 * skg1_295[k]
                   + f_3 * pc_y[k] * skh_414[k];

        t_554[k] = f_14 * sih_288[k]
                   + f_3 * pc_z[k] * skh_414[k];

        t_555[k] = f_11 * sih_311[k]
                   + f_4 * skg0_297[k]
                   - f_5 * skg1_297[k]
                   + f_3 * pc_y[k] * skh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, sih_312, sih_313, sih_314, skg0_298, \
                         skg0_299, skg1_298, skg1_299, skh_417, skh_418, \
                         skh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * sih_312[k]
                   + f_6 * skg0_298[k]
                   - f_7 * skg1_298[k]
                   + f_3 * pc_y[k] * skh_417[k];

        t_557[k] = f_11 * sih_313[k]
                   + f_8 * skg0_299[k]
                   - f_9 * skg1_299[k]
                   + f_3 * pc_y[k] * skh_418[k];

        t_558[k] = f_11 * sih_314[k]
                   + f_3 * pc_y[k] * skh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_y, pc_x, pc_y, pc_z, sii0_419, \
                         sih_294, sih_420, sii1_419, skg0_300, skg1_300, \
                         skh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pb_y[k] * sii0_419[k]
                   - f_10 * pc_y[k] * sii1_419[k];

        t_560[k] = f_12 * sih_420[k]
                   + f_1 * skg0_300[k]
                   - f_2 * skg1_300[k]
                   + f_3 * pc_x[k] * skh_420[k];

        t_561[k] = f_3 * pc_y[k] * skh_420[k];

        t_562[k] = f_16 * sih_294[k]
                   + f_3 * pc_z[k] * skh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, sih_423, sih_425, skg0_303, \
                         skg0_305, skg1_303, skg1_305, skh_422, skh_423, \
                         skh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_12 * sih_423[k]
                   + f_4 * skg0_303[k]
                   - f_5 * skg1_303[k]
                   + f_3 * pc_x[k] * skh_423[k];

        t_564[k] = f_3 * pc_y[k] * skh_422[k];

        t_565[k] = f_12 * sih_425[k]
                   + f_4 * skg0_305[k]
                   - f_5 * skg1_305[k]
                   + f_3 * pc_x[k] * skh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_x, pc_y, pc_z, sih_297, sih_426, skg0_306, \
                         skg1_306, skh_423, skh_425, skh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_12 * sih_426[k]
                   + f_6 * skg0_306[k]
                   - f_7 * skg1_306[k]
                   + f_3 * pc_x[k] * skh_426[k];

        t_567[k] = f_16 * sih_297[k]
                   + f_3 * pc_z[k] * skh_423[k];

        t_568[k] = f_3 * pc_y[k] * skh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_z, sih_300, sih_429, sih_430, skg0_309, \
                         skg0_310, skg1_309, skg1_310, skh_426, skh_429, \
                         skh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_12 * sih_429[k]
                   + f_6 * skg0_309[k]
                   - f_7 * skg1_309[k]
                   + f_3 * pc_x[k] * skh_429[k];

        t_570[k] = f_12 * sih_430[k]
                   + f_8 * skg0_310[k]
                   - f_9 * skg1_310[k]
                   + f_3 * pc_x[k] * skh_430[k];

        t_571[k] = f_16 * sih_300[k]
                   + f_3 * pc_z[k] * skh_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pc_x, pc_y, sih_432, sih_434, skg0_312, \
                         skg0_314, skg1_312, skg1_314, skh_429, skh_432, \
                         skh_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_12 * sih_432[k]
                   + f_8 * skg0_312[k]
                   - f_9 * skg1_312[k]
                   + f_3 * pc_x[k] * skh_432[k];

        t_573[k] = f_3 * pc_y[k] * skh_429[k];

        t_574[k] = f_12 * sih_434[k]
                   + f_8 * skg0_314[k]
                   - f_9 * skg1_314[k]
                   + f_3 * pc_x[k] * skh_434[k];
    }
}

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_420 = buffer.data(sii0 + 420);
    const auto *sii0_423 = buffer.data(sii0 + 423);
    const auto *sii0_426 = buffer.data(sii0 + 426);
    const auto *sii0_430 = buffer.data(sii0 + 430);
    const auto *sii0_588 = buffer.data(sii0 + 588);
    const auto *sii0_591 = buffer.data(sii0 + 591);
    const auto *sii0_593 = buffer.data(sii0 + 593);
    const auto *sii0_594 = buffer.data(sii0 + 594);
    const auto *sii0_597 = buffer.data(sii0 + 597);
    const auto *sii0_598 = buffer.data(sii0 + 598);
    const auto *sii0_600 = buffer.data(sii0 + 600);
    const auto *sii0_602 = buffer.data(sii0 + 602);
    const auto *sii0_609 = buffer.data(sii0 + 609);
    const auto *sii0_611 = buffer.data(sii0 + 611);
    const auto *sii0_612 = buffer.data(sii0 + 612);
    const auto *sii0_613 = buffer.data(sii0 + 613);
    const auto *sii0_615 = buffer.data(sii0 + 615);
    const auto *sii0_621 = buffer.data(sii0 + 621);
    const auto *sii0_625 = buffer.data(sii0 + 625);
    const auto *sii0_628 = buffer.data(sii0 + 628);
    const auto *sii0_630 = buffer.data(sii0 + 630);
    const auto *sii0_637 = buffer.data(sii0 + 637);
    const auto *sii0_639 = buffer.data(sii0 + 639);
    const auto *sii0_640 = buffer.data(sii0 + 640);
    const auto *sii0_641 = buffer.data(sii0 + 641);
    const auto *sii0_643 = buffer.data(sii0 + 643);
    const auto *sii0_644 = buffer.data(sii0 + 644);
    const auto *sii0_647 = buffer.data(sii0 + 647);
    const auto *sii0_649 = buffer.data(sii0 + 649);
    const auto *sii0_650 = buffer.data(sii0 + 650);
    const auto *sii0_653 = buffer.data(sii0 + 653);
    const auto *sii0_654 = buffer.data(sii0 + 654);
    const auto *sii0_656 = buffer.data(sii0 + 656);
    const auto *sii0_658 = buffer.data(sii0 + 658);
    const auto *sii0_665 = buffer.data(sii0 + 665);
    const auto *sii0_667 = buffer.data(sii0 + 667);
    const auto *sii0_668 = buffer.data(sii0 + 668);
    const auto *sii0_669 = buffer.data(sii0 + 669);
    const auto *sii0_671 = buffer.data(sii0 + 671);
    const auto *sii0_672 = buffer.data(sii0 + 672);
    const auto *sii0_675 = buffer.data(sii0 + 675);
    const auto *sii0_677 = buffer.data(sii0 + 677);
    const auto *sii0_678 = buffer.data(sii0 + 678);
    const auto *sii0_681 = buffer.data(sii0 + 681);
    const auto *sii0_682 = buffer.data(sii0 + 682);
    const auto *sii0_684 = buffer.data(sii0 + 684);
    const auto *sii0_686 = buffer.data(sii0 + 686);
    const auto *sii0_693 = buffer.data(sii0 + 693);
    const auto *sii0_695 = buffer.data(sii0 + 695);
    const auto *sii0_696 = buffer.data(sii0 + 696);
    const auto *sii0_697 = buffer.data(sii0 + 697);

    const auto *sih_309 = buffer.data(sih + 309);
    const auto *sih_314 = buffer.data(sih + 314);
    const auto *sih_315 = buffer.data(sih + 315);
    const auto *sih_317 = buffer.data(sih + 317);
    const auto *sih_318 = buffer.data(sih + 318);
    const auto *sih_320 = buffer.data(sih + 320);
    const auto *sih_321 = buffer.data(sih + 321);
    const auto *sih_324 = buffer.data(sih + 324);
    const auto *sih_330 = buffer.data(sih + 330);
    const auto *sih_335 = buffer.data(sih + 335);
    const auto *sih_336 = buffer.data(sih + 336);
    const auto *sih_338 = buffer.data(sih + 338);
    const auto *sih_339 = buffer.data(sih + 339);
    const auto *sih_341 = buffer.data(sih + 341);
    const auto *sih_342 = buffer.data(sih + 342);
    const auto *sih_345 = buffer.data(sih + 345);
    const auto *sih_351 = buffer.data(sih + 351);
    const auto *sih_356 = buffer.data(sih + 356);
    const auto *sih_357 = buffer.data(sih + 357);
    const auto *sih_359 = buffer.data(sih + 359);
    const auto *sih_360 = buffer.data(sih + 360);
    const auto *sih_362 = buffer.data(sih + 362);
    const auto *sih_363 = buffer.data(sih + 363);
    const auto *sih_366 = buffer.data(sih + 366);
    const auto *sih_372 = buffer.data(sih + 372);
    const auto *sih_377 = buffer.data(sih + 377);
    const auto *sih_378 = buffer.data(sih + 378);
    const auto *sih_380 = buffer.data(sih + 380);
    const auto *sih_383 = buffer.data(sih + 383);
    const auto *sih_387 = buffer.data(sih + 387);
    const auto *sih_435 = buffer.data(sih + 435);
    const auto *sih_436 = buffer.data(sih + 436);
    const auto *sih_437 = buffer.data(sih + 437);
    const auto *sih_438 = buffer.data(sih + 438);
    const auto *sih_439 = buffer.data(sih + 439);
    const auto *sih_440 = buffer.data(sih + 440);
    const auto *sih_441 = buffer.data(sih + 441);
    const auto *sih_444 = buffer.data(sih + 444);
    const auto *sih_446 = buffer.data(sih + 446);
    const auto *sih_447 = buffer.data(sih + 447);
    const auto *sih_450 = buffer.data(sih + 450);
    const auto *sih_451 = buffer.data(sih + 451);
    const auto *sih_453 = buffer.data(sih + 453);
    const auto *sih_455 = buffer.data(sih + 455);
    const auto *sih_456 = buffer.data(sih + 456);
    const auto *sih_457 = buffer.data(sih + 457);
    const auto *sih_458 = buffer.data(sih + 458);
    const auto *sih_459 = buffer.data(sih + 459);
    const auto *sih_460 = buffer.data(sih + 460);
    const auto *sih_461 = buffer.data(sih + 461);
    const auto *sih_467 = buffer.data(sih + 467);
    const auto *sih_471 = buffer.data(sih + 471);
    const auto *sih_474 = buffer.data(sih + 474);
    const auto *sih_476 = buffer.data(sih + 476);
    const auto *sih_477 = buffer.data(sih + 477);
    const auto *sih_478 = buffer.data(sih + 478);
    const auto *sih_479 = buffer.data(sih + 479);
    const auto *sih_480 = buffer.data(sih + 480);
    const auto *sih_481 = buffer.data(sih + 481);
    const auto *sih_482 = buffer.data(sih + 482);
    const auto *sih_483 = buffer.data(sih + 483);
    const auto *sih_486 = buffer.data(sih + 486);
    const auto *sih_488 = buffer.data(sih + 488);
    const auto *sih_489 = buffer.data(sih + 489);
    const auto *sih_492 = buffer.data(sih + 492);
    const auto *sih_493 = buffer.data(sih + 493);
    const auto *sih_495 = buffer.data(sih + 495);
    const auto *sih_497 = buffer.data(sih + 497);
    const auto *sih_498 = buffer.data(sih + 498);
    const auto *sih_499 = buffer.data(sih + 499);
    const auto *sih_500 = buffer.data(sih + 500);
    const auto *sih_501 = buffer.data(sih + 501);
    const auto *sih_502 = buffer.data(sih + 502);
    const auto *sih_503 = buffer.data(sih + 503);
    const auto *sih_504 = buffer.data(sih + 504);
    const auto *sih_507 = buffer.data(sih + 507);
    const auto *sih_509 = buffer.data(sih + 509);
    const auto *sih_510 = buffer.data(sih + 510);
    const auto *sih_513 = buffer.data(sih + 513);
    const auto *sih_514 = buffer.data(sih + 514);
    const auto *sih_516 = buffer.data(sih + 516);
    const auto *sih_518 = buffer.data(sih + 518);
    const auto *sih_519 = buffer.data(sih + 519);
    const auto *sih_520 = buffer.data(sih + 520);
    const auto *sih_521 = buffer.data(sih + 521);
    const auto *sih_522 = buffer.data(sih + 522);
    const auto *sih_523 = buffer.data(sih + 523);
    const auto *sih_524 = buffer.data(sih + 524);

    const auto *sii1_420 = buffer.data(sii1 + 420);
    const auto *sii1_423 = buffer.data(sii1 + 423);
    const auto *sii1_426 = buffer.data(sii1 + 426);
    const auto *sii1_430 = buffer.data(sii1 + 430);
    const auto *sii1_588 = buffer.data(sii1 + 588);
    const auto *sii1_591 = buffer.data(sii1 + 591);
    const auto *sii1_593 = buffer.data(sii1 + 593);
    const auto *sii1_594 = buffer.data(sii1 + 594);
    const auto *sii1_597 = buffer.data(sii1 + 597);
    const auto *sii1_598 = buffer.data(sii1 + 598);
    const auto *sii1_600 = buffer.data(sii1 + 600);
    const auto *sii1_602 = buffer.data(sii1 + 602);
    const auto *sii1_609 = buffer.data(sii1 + 609);
    const auto *sii1_611 = buffer.data(sii1 + 611);
    const auto *sii1_612 = buffer.data(sii1 + 612);
    const auto *sii1_613 = buffer.data(sii1 + 613);
    const auto *sii1_615 = buffer.data(sii1 + 615);
    const auto *sii1_621 = buffer.data(sii1 + 621);
    const auto *sii1_625 = buffer.data(sii1 + 625);
    const auto *sii1_628 = buffer.data(sii1 + 628);
    const auto *sii1_630 = buffer.data(sii1 + 630);
    const auto *sii1_637 = buffer.data(sii1 + 637);
    const auto *sii1_639 = buffer.data(sii1 + 639);
    const auto *sii1_640 = buffer.data(sii1 + 640);
    const auto *sii1_641 = buffer.data(sii1 + 641);
    const auto *sii1_643 = buffer.data(sii1 + 643);
    const auto *sii1_644 = buffer.data(sii1 + 644);
    const auto *sii1_647 = buffer.data(sii1 + 647);
    const auto *sii1_649 = buffer.data(sii1 + 649);
    const auto *sii1_650 = buffer.data(sii1 + 650);
    const auto *sii1_653 = buffer.data(sii1 + 653);
    const auto *sii1_654 = buffer.data(sii1 + 654);
    const auto *sii1_656 = buffer.data(sii1 + 656);
    const auto *sii1_658 = buffer.data(sii1 + 658);
    const auto *sii1_665 = buffer.data(sii1 + 665);
    const auto *sii1_667 = buffer.data(sii1 + 667);
    const auto *sii1_668 = buffer.data(sii1 + 668);
    const auto *sii1_669 = buffer.data(sii1 + 669);
    const auto *sii1_671 = buffer.data(sii1 + 671);
    const auto *sii1_672 = buffer.data(sii1 + 672);
    const auto *sii1_675 = buffer.data(sii1 + 675);
    const auto *sii1_677 = buffer.data(sii1 + 677);
    const auto *sii1_678 = buffer.data(sii1 + 678);
    const auto *sii1_681 = buffer.data(sii1 + 681);
    const auto *sii1_682 = buffer.data(sii1 + 682);
    const auto *sii1_684 = buffer.data(sii1 + 684);
    const auto *sii1_686 = buffer.data(sii1 + 686);
    const auto *sii1_693 = buffer.data(sii1 + 693);
    const auto *sii1_695 = buffer.data(sii1 + 695);
    const auto *sii1_696 = buffer.data(sii1 + 696);
    const auto *sii1_697 = buffer.data(sii1 + 697);

    const auto *skg0_310 = buffer.data(skg0 + 310);
    const auto *skg0_312 = buffer.data(skg0 + 312);
    const auto *skg0_313 = buffer.data(skg0 + 313);
    const auto *skg0_314 = buffer.data(skg0 + 314);

    const auto *skg1_310 = buffer.data(skg1 + 310);
    const auto *skg1_312 = buffer.data(skg1 + 312);
    const auto *skg1_313 = buffer.data(skg1 + 313);
    const auto *skg1_314 = buffer.data(skg1 + 314);

    const auto *skh_435 = buffer.data(skh + 435);
    const auto *skh_436 = buffer.data(skh + 436);
    const auto *skh_437 = buffer.data(skh + 437);
    const auto *skh_438 = buffer.data(skh + 438);
    const auto *skh_439 = buffer.data(skh + 439);
    const auto *skh_440 = buffer.data(skh + 440);
    const auto *skh_441 = buffer.data(skh + 441);
    const auto *skh_443 = buffer.data(skh + 443);
    const auto *skh_444 = buffer.data(skh + 444);
    const auto *skh_446 = buffer.data(skh + 446);
    const auto *skh_447 = buffer.data(skh + 447);
    const auto *skh_450 = buffer.data(skh + 450);
    const auto *skh_456 = buffer.data(skh + 456);
    const auto *skh_457 = buffer.data(skh + 457);
    const auto *skh_458 = buffer.data(skh + 458);
    const auto *skh_459 = buffer.data(skh + 459);
    const auto *skh_460 = buffer.data(skh + 460);
    const auto *skh_461 = buffer.data(skh + 461);
    const auto *skh_462 = buffer.data(skh + 462);
    const auto *skh_464 = buffer.data(skh + 464);
    const auto *skh_465 = buffer.data(skh + 465);
    const auto *skh_467 = buffer.data(skh + 467);
    const auto *skh_468 = buffer.data(skh + 468);
    const auto *skh_471 = buffer.data(skh + 471);
    const auto *skh_477 = buffer.data(skh + 477);
    const auto *skh_478 = buffer.data(skh + 478);
    const auto *skh_479 = buffer.data(skh + 479);
    const auto *skh_480 = buffer.data(skh + 480);
    const auto *skh_481 = buffer.data(skh + 481);
    const auto *skh_482 = buffer.data(skh + 482);
    const auto *skh_483 = buffer.data(skh + 483);
    const auto *skh_485 = buffer.data(skh + 485);
    const auto *skh_486 = buffer.data(skh + 486);
    const auto *skh_488 = buffer.data(skh + 488);
    const auto *skh_489 = buffer.data(skh + 489);
    const auto *skh_492 = buffer.data(skh + 492);
    const auto *skh_498 = buffer.data(skh + 498);
    const auto *skh_499 = buffer.data(skh + 499);
    const auto *skh_500 = buffer.data(skh + 500);
    const auto *skh_501 = buffer.data(skh + 501);
    const auto *skh_502 = buffer.data(skh + 502);
    const auto *skh_503 = buffer.data(skh + 503);
    const auto *skh_504 = buffer.data(skh + 504);
    const auto *skh_506 = buffer.data(skh + 506);
    const auto *skh_507 = buffer.data(skh + 507);
    const auto *skh_509 = buffer.data(skh + 509);
    const auto *skh_510 = buffer.data(skh + 510);
    const auto *skh_513 = buffer.data(skh + 513);
    const auto *skh_519 = buffer.data(skh + 519);
    const auto *skh_520 = buffer.data(skh + 520);
    const auto *skh_521 = buffer.data(skh + 521);
    const auto *skh_522 = buffer.data(skh + 522);
    const auto *skh_523 = buffer.data(skh + 523);
    const auto *skh_524 = buffer.data(skh + 524);

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, pc_x, sih_435, sih_436, sih_437, \
                         sih_438, sih_439, skh_435, skh_436, skh_437, skh_438, \
                         skh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_12 * sih_435[k]
                   + f_3 * pc_x[k] * skh_435[k];

        t_576[k] = f_12 * sih_436[k]
                   + f_3 * pc_x[k] * skh_436[k];

        t_577[k] = f_12 * sih_437[k]
                   + f_3 * pc_x[k] * skh_437[k];

        t_578[k] = f_12 * sih_438[k]
                   + f_3 * pc_x[k] * skh_438[k];

        t_579[k] = f_12 * sih_439[k]
                   + f_3 * pc_x[k] * skh_439[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, pc_z, sih_309, sih_440, \
                         skg0_310, skg0_312, skg1_310, skg1_312, skh_435, skh_437, \
                         skh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_12 * sih_440[k]
                   + f_3 * pc_x[k] * skh_440[k];

        t_581[k] = f_1 * skg0_310[k]
                   - f_2 * skg1_310[k]
                   + f_3 * pc_y[k] * skh_435[k];

        t_582[k] = f_16 * sih_309[k]
                   + f_3 * pc_z[k] * skh_435[k];

        t_583[k] = f_4 * skg0_312[k]
                   - f_5 * skg1_312[k]
                   + f_3 * pc_y[k] * skh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, sih_314, skg0_313, skg0_314, \
                         skg1_313, skg1_314, skh_438, skh_439, \
                         skh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * skg0_313[k]
                   - f_7 * skg1_313[k]
                   + f_3 * pc_y[k] * skh_438[k];

        t_585[k] = f_8 * skg0_314[k]
                   - f_9 * skg1_314[k]
                   + f_3 * pc_y[k] * skh_439[k];

        t_586[k] = f_3 * pc_y[k] * skh_440[k];

        t_587[k] = f_16 * sih_314[k]
                   + f_1 * skg0_314[k]
                   - f_2 * skg1_314[k]
                   + f_3 * pc_z[k] * skh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pb_x, pc_x, pc_y, pc_z, sii0_588, \
                         sii0_591, sih_315, sih_441, sih_444, sii1_588, sii1_591, \
                         skh_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pb_x[k] * sii0_588[k]
                   + f_15 * sih_441[k]
                   - f_10 * pc_x[k] * sii1_588[k];

        t_589[k] = f_15 * sih_315[k]
                   + f_3 * pc_y[k] * skh_441[k];

        t_590[k] = f_3 * pc_z[k] * skh_441[k];

        t_591[k] = pb_x[k] * sii0_591[k]
                   + f_14 * sih_444[k]
                   - f_10 * pc_x[k] * sii1_591[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, pb_x, pc_x, pc_y, sii0_593, sii0_594, sih_317, \
                         sih_446, sih_447, sii1_593, sii1_594, \
                         skh_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_15 * sih_317[k]
                   + f_3 * pc_y[k] * skh_443[k];

        t_593[k] = pb_x[k] * sii0_593[k]
                   + f_14 * sih_446[k]
                   - f_10 * pc_x[k] * sii1_593[k];

        t_594[k] = pb_x[k] * sii0_594[k]
                   + f_13 * sih_447[k]
                   - f_10 * pc_x[k] * sii1_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pb_x, pc_x, pc_y, pc_z, sii0_597, sih_320, \
                         sih_450, sii1_597, skh_444, skh_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_3 * pc_z[k] * skh_444[k];

        t_596[k] = f_15 * sih_320[k]
                   + f_3 * pc_y[k] * skh_446[k];

        t_597[k] = pb_x[k] * sii0_597[k]
                   + f_13 * sih_450[k]
                   - f_10 * pc_x[k] * sii1_597[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pb_x, pc_x, pc_z, sii0_598, sii0_600, sih_451, \
                         sih_453, sii1_598, sii1_600, skh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = pb_x[k] * sii0_598[k]
                   + f_12 * sih_451[k]
                   - f_10 * pc_x[k] * sii1_598[k];

        t_599[k] = f_3 * pc_z[k] * skh_447[k];

        t_600[k] = pb_x[k] * sii0_600[k]
                   + f_12 * sih_453[k]
                   - f_10 * pc_x[k] * sii1_600[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pb_x, pc_x, pc_y, sii0_602, sih_324, \
                         sih_455, sih_456, sih_457, sii1_602, skh_450, skh_456, \
                         skh_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_15 * sih_324[k]
                   + f_3 * pc_y[k] * skh_450[k];

        t_602[k] = pb_x[k] * sii0_602[k]
                   + f_12 * sih_455[k]
                   - f_10 * pc_x[k] * sii1_602[k];

        t_603[k] = f_11 * sih_456[k]
                   + f_3 * pc_x[k] * skh_456[k];

        t_604[k] = f_11 * sih_457[k]
                   + f_3 * pc_x[k] * skh_457[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pc_x, sih_458, sih_459, sih_460, sih_461, \
                         skh_458, skh_459, skh_460, skh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_11 * sih_458[k]
                   + f_3 * pc_x[k] * skh_458[k];

        t_606[k] = f_11 * sih_459[k]
                   + f_3 * pc_x[k] * skh_459[k];

        t_607[k] = f_11 * sih_460[k]
                   + f_3 * pc_x[k] * skh_460[k];

        t_608[k] = f_11 * sih_461[k]
                   + f_3 * pc_x[k] * skh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pb_x, pc_x, pc_z, sii0_609, sii0_611, \
                         sii0_612, sii1_609, sii1_611, sii1_612, \
                         skh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = pb_x[k] * sii0_609[k]
                   - f_10 * pc_x[k] * sii1_609[k];

        t_610[k] = f_3 * pc_z[k] * skh_456[k];

        t_611[k] = pb_x[k] * sii0_611[k]
                   - f_10 * pc_x[k] * sii1_611[k];

        t_612[k] = pb_x[k] * sii0_612[k]
                   - f_10 * pc_x[k] * sii1_612[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, pb_x, pc_x, pc_y, sii0_613, sii0_615, sih_335, \
                         sii1_613, sii1_615, skh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = pb_x[k] * sii0_613[k]
                   - f_10 * pc_x[k] * sii1_613[k];

        t_614[k] = f_15 * sih_335[k]
                   + f_3 * pc_y[k] * skh_461[k];

        t_615[k] = pb_x[k] * sii0_615[k]
                   - f_10 * pc_x[k] * sii1_615[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pb_z, pc_y, pc_z, sii0_420, sii0_423, \
                         sih_315, sih_336, sii1_420, sii1_423, \
                         skh_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pb_z[k] * sii0_420[k]
                   - f_10 * pc_z[k] * sii1_420[k];

        t_617[k] = f_16 * sih_336[k]
                   + f_3 * pc_y[k] * skh_462[k];

        t_618[k] = f_11 * sih_315[k]
                   + f_3 * pc_z[k] * skh_462[k];

        t_619[k] = pb_z[k] * sii0_423[k]
                   - f_10 * pc_z[k] * sii1_423[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pb_x, pb_z, pc_x, pc_y, pc_z, sii0_426, \
                         sii0_621, sih_338, sih_467, sii1_426, sii1_621, \
                         skh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_16 * sih_338[k]
                   + f_3 * pc_y[k] * skh_464[k];

        t_621[k] = pb_x[k] * sii0_621[k]
                   + f_14 * sih_467[k]
                   - f_10 * pc_x[k] * sii1_621[k];

        t_622[k] = pb_z[k] * sii0_426[k]
                   - f_10 * pc_z[k] * sii1_426[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pb_x, pc_x, pc_y, pc_z, sii0_625, sih_318, \
                         sih_341, sih_471, sii1_625, skh_465, skh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_11 * sih_318[k]
                   + f_3 * pc_z[k] * skh_465[k];

        t_624[k] = f_16 * sih_341[k]
                   + f_3 * pc_y[k] * skh_467[k];

        t_625[k] = pb_x[k] * sii0_625[k]
                   + f_13 * sih_471[k]
                   - f_10 * pc_x[k] * sii1_625[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pb_x, pb_z, pc_x, pc_z, sii0_430, sii0_628, \
                         sih_321, sih_474, sii1_430, sii1_628, \
                         skh_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pb_z[k] * sii0_430[k]
                   - f_10 * pc_z[k] * sii1_430[k];

        t_627[k] = f_11 * sih_321[k]
                   + f_3 * pc_z[k] * skh_468[k];

        t_628[k] = pb_x[k] * sii0_628[k]
                   + f_12 * sih_474[k]
                   - f_10 * pc_x[k] * sii1_628[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pb_x, pc_x, pc_y, sii0_630, sih_345, \
                         sih_476, sih_477, sih_478, sii1_630, skh_471, skh_477, \
                         skh_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_16 * sih_345[k]
                   + f_3 * pc_y[k] * skh_471[k];

        t_630[k] = pb_x[k] * sii0_630[k]
                   + f_12 * sih_476[k]
                   - f_10 * pc_x[k] * sii1_630[k];

        t_631[k] = f_11 * sih_477[k]
                   + f_3 * pc_x[k] * skh_477[k];

        t_632[k] = f_11 * sih_478[k]
                   + f_3 * pc_x[k] * skh_478[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, pc_x, sih_479, sih_480, sih_481, sih_482, \
                         skh_479, skh_480, skh_481, skh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_11 * sih_479[k]
                   + f_3 * pc_x[k] * skh_479[k];

        t_634[k] = f_11 * sih_480[k]
                   + f_3 * pc_x[k] * skh_480[k];

        t_635[k] = f_11 * sih_481[k]
                   + f_3 * pc_x[k] * skh_481[k];

        t_636[k] = f_11 * sih_482[k]
                   + f_3 * pc_x[k] * skh_482[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, pb_x, pc_x, pc_z, sii0_637, sii0_639, \
                         sii0_640, sih_330, sii1_637, sii1_639, sii1_640, \
                         skh_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = pb_x[k] * sii0_637[k]
                   - f_10 * pc_x[k] * sii1_637[k];

        t_638[k] = f_11 * sih_330[k]
                   + f_3 * pc_z[k] * skh_477[k];

        t_639[k] = pb_x[k] * sii0_639[k]
                   - f_10 * pc_x[k] * sii1_639[k];

        t_640[k] = pb_x[k] * sii0_640[k]
                   - f_10 * pc_x[k] * sii1_640[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, pb_x, pc_x, pc_y, sii0_641, sii0_643, \
                         sii0_644, sih_356, sih_483, sii1_641, sii1_643, sii1_644, \
                         skh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = pb_x[k] * sii0_641[k]
                   - f_10 * pc_x[k] * sii1_641[k];

        t_642[k] = f_16 * sih_356[k]
                   + f_3 * pc_y[k] * skh_482[k];

        t_643[k] = pb_x[k] * sii0_643[k]
                   - f_10 * pc_x[k] * sii1_643[k];

        t_644[k] = pb_x[k] * sii0_644[k]
                   + f_15 * sih_483[k]
                   - f_10 * pc_x[k] * sii1_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, pb_x, pc_x, pc_y, pc_z, sii0_647, \
                         sih_336, sih_357, sih_359, sih_486, sii1_647, skh_483, \
                         skh_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_14 * sih_357[k]
                   + f_3 * pc_y[k] * skh_483[k];

        t_646[k] = f_12 * sih_336[k]
                   + f_3 * pc_z[k] * skh_483[k];

        t_647[k] = pb_x[k] * sii0_647[k]
                   + f_14 * sih_486[k]
                   - f_10 * pc_x[k] * sii1_647[k];

        t_648[k] = f_14 * sih_359[k]
                   + f_3 * pc_y[k] * skh_485[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pb_x, pc_x, pc_z, sii0_649, sii0_650, sih_339, \
                         sih_488, sih_489, sii1_649, sii1_650, \
                         skh_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = pb_x[k] * sii0_649[k]
                   + f_14 * sih_488[k]
                   - f_10 * pc_x[k] * sii1_649[k];

        t_650[k] = pb_x[k] * sii0_650[k]
                   + f_13 * sih_489[k]
                   - f_10 * pc_x[k] * sii1_650[k];

        t_651[k] = f_12 * sih_339[k]
                   + f_3 * pc_z[k] * skh_486[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, pb_x, pc_x, pc_y, sii0_653, sii0_654, sih_362, \
                         sih_492, sih_493, sii1_653, sii1_654, \
                         skh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_14 * sih_362[k]
                   + f_3 * pc_y[k] * skh_488[k];

        t_653[k] = pb_x[k] * sii0_653[k]
                   + f_13 * sih_492[k]
                   - f_10 * pc_x[k] * sii1_653[k];

        t_654[k] = pb_x[k] * sii0_654[k]
                   + f_12 * sih_493[k]
                   - f_10 * pc_x[k] * sii1_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pb_x, pc_x, pc_y, pc_z, sii0_656, sih_342, \
                         sih_366, sih_495, sii1_656, skh_489, skh_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_12 * sih_342[k]
                   + f_3 * pc_z[k] * skh_489[k];

        t_656[k] = pb_x[k] * sii0_656[k]
                   + f_12 * sih_495[k]
                   - f_10 * pc_x[k] * sii1_656[k];

        t_657[k] = f_14 * sih_366[k]
                   + f_3 * pc_y[k] * skh_492[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pb_x, pc_x, sii0_658, sih_497, sih_498, \
                         sih_499, sih_500, sii1_658, skh_498, skh_499, \
                         skh_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = pb_x[k] * sii0_658[k]
                   + f_12 * sih_497[k]
                   - f_10 * pc_x[k] * sii1_658[k];

        t_659[k] = f_11 * sih_498[k]
                   + f_3 * pc_x[k] * skh_498[k];

        t_660[k] = f_11 * sih_499[k]
                   + f_3 * pc_x[k] * skh_499[k];

        t_661[k] = f_11 * sih_500[k]
                   + f_3 * pc_x[k] * skh_500[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, pb_x, pc_x, sii0_665, sih_501, sih_502, \
                         sih_503, sii1_665, skh_501, skh_502, skh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_11 * sih_501[k]
                   + f_3 * pc_x[k] * skh_501[k];

        t_663[k] = f_11 * sih_502[k]
                   + f_3 * pc_x[k] * skh_502[k];

        t_664[k] = f_11 * sih_503[k]
                   + f_3 * pc_x[k] * skh_503[k];

        t_665[k] = pb_x[k] * sii0_665[k]
                   - f_10 * pc_x[k] * sii1_665[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pb_x, pc_x, pc_z, sii0_667, sii0_668, \
                         sii0_669, sih_351, sii1_667, sii1_668, sii1_669, \
                         skh_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_12 * sih_351[k]
                   + f_3 * pc_z[k] * skh_498[k];

        t_667[k] = pb_x[k] * sii0_667[k]
                   - f_10 * pc_x[k] * sii1_667[k];

        t_668[k] = pb_x[k] * sii0_668[k]
                   - f_10 * pc_x[k] * sii1_668[k];

        t_669[k] = pb_x[k] * sii0_669[k]
                   - f_10 * pc_x[k] * sii1_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pb_x, pc_x, pc_y, sii0_671, sii0_672, \
                         sih_377, sih_378, sih_504, sii1_671, sii1_672, skh_503, \
                         skh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * sih_377[k]
                   + f_3 * pc_y[k] * skh_503[k];

        t_671[k] = pb_x[k] * sii0_671[k]
                   - f_10 * pc_x[k] * sii1_671[k];

        t_672[k] = pb_x[k] * sii0_672[k]
                   + f_15 * sih_504[k]
                   - f_10 * pc_x[k] * sii1_672[k];

        t_673[k] = f_13 * sih_378[k]
                   + f_3 * pc_y[k] * skh_504[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, pb_x, pc_x, pc_y, pc_z, sii0_675, sih_357, \
                         sih_380, sih_507, sii1_675, skh_504, skh_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_13 * sih_357[k]
                   + f_3 * pc_z[k] * skh_504[k];

        t_675[k] = pb_x[k] * sii0_675[k]
                   + f_14 * sih_507[k]
                   - f_10 * pc_x[k] * sii1_675[k];

        t_676[k] = f_13 * sih_380[k]
                   + f_3 * pc_y[k] * skh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pb_x, pc_x, pc_z, sii0_677, sii0_678, sih_360, \
                         sih_509, sih_510, sii1_677, sii1_678, \
                         skh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = pb_x[k] * sii0_677[k]
                   + f_14 * sih_509[k]
                   - f_10 * pc_x[k] * sii1_677[k];

        t_678[k] = pb_x[k] * sii0_678[k]
                   + f_13 * sih_510[k]
                   - f_10 * pc_x[k] * sii1_678[k];

        t_679[k] = f_13 * sih_360[k]
                   + f_3 * pc_z[k] * skh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pb_x, pc_x, pc_y, sii0_681, sii0_682, sih_383, \
                         sih_513, sih_514, sii1_681, sii1_682, \
                         skh_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * sih_383[k]
                   + f_3 * pc_y[k] * skh_509[k];

        t_681[k] = pb_x[k] * sii0_681[k]
                   + f_13 * sih_513[k]
                   - f_10 * pc_x[k] * sii1_681[k];

        t_682[k] = pb_x[k] * sii0_682[k]
                   + f_12 * sih_514[k]
                   - f_10 * pc_x[k] * sii1_682[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pb_x, pc_x, pc_y, pc_z, sii0_684, sih_363, \
                         sih_387, sih_516, sii1_684, skh_510, skh_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * sih_363[k]
                   + f_3 * pc_z[k] * skh_510[k];

        t_684[k] = pb_x[k] * sii0_684[k]
                   + f_12 * sih_516[k]
                   - f_10 * pc_x[k] * sii1_684[k];

        t_685[k] = f_13 * sih_387[k]
                   + f_3 * pc_y[k] * skh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pb_x, pc_x, sii0_686, sih_518, sih_519, \
                         sih_520, sih_521, sii1_686, skh_519, skh_520, \
                         skh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = pb_x[k] * sii0_686[k]
                   + f_12 * sih_518[k]
                   - f_10 * pc_x[k] * sii1_686[k];

        t_687[k] = f_11 * sih_519[k]
                   + f_3 * pc_x[k] * skh_519[k];

        t_688[k] = f_11 * sih_520[k]
                   + f_3 * pc_x[k] * skh_520[k];

        t_689[k] = f_11 * sih_521[k]
                   + f_3 * pc_x[k] * skh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pb_x, pc_x, sii0_693, sih_522, sih_523, \
                         sih_524, sii1_693, skh_522, skh_523, skh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_11 * sih_522[k]
                   + f_3 * pc_x[k] * skh_522[k];

        t_691[k] = f_11 * sih_523[k]
                   + f_3 * pc_x[k] * skh_523[k];

        t_692[k] = f_11 * sih_524[k]
                   + f_3 * pc_x[k] * skh_524[k];

        t_693[k] = pb_x[k] * sii0_693[k]
                   - f_10 * pc_x[k] * sii1_693[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pb_x, pc_x, pc_z, sii0_695, sii0_696, \
                         sii0_697, sih_372, sii1_695, sii1_696, sii1_697, \
                         skh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * sih_372[k]
                   + f_3 * pc_z[k] * skh_519[k];

        t_695[k] = pb_x[k] * sii0_695[k]
                   - f_10 * pc_x[k] * sii1_695[k];

        t_696[k] = pb_x[k] * sii0_696[k]
                   - f_10 * pc_x[k] * sii1_696[k];

        t_697[k] = pb_x[k] * sii0_697[k]
                   - f_10 * pc_x[k] * sii1_697[k];
    }
}

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_560 = buffer.data(sii0 + 560);
    const auto *sii0_565 = buffer.data(sii0 + 565);
    const auto *sii0_569 = buffer.data(sii0 + 569);
    const auto *sii0_574 = buffer.data(sii0 + 574);
    const auto *sii0_588 = buffer.data(sii0 + 588);
    const auto *sii0_591 = buffer.data(sii0 + 591);
    const auto *sii0_594 = buffer.data(sii0 + 594);
    const auto *sii0_699 = buffer.data(sii0 + 699);
    const auto *sii0_700 = buffer.data(sii0 + 700);
    const auto *sii0_703 = buffer.data(sii0 + 703);
    const auto *sii0_705 = buffer.data(sii0 + 705);
    const auto *sii0_706 = buffer.data(sii0 + 706);
    const auto *sii0_709 = buffer.data(sii0 + 709);
    const auto *sii0_710 = buffer.data(sii0 + 710);
    const auto *sii0_712 = buffer.data(sii0 + 712);
    const auto *sii0_714 = buffer.data(sii0 + 714);
    const auto *sii0_721 = buffer.data(sii0 + 721);
    const auto *sii0_723 = buffer.data(sii0 + 723);
    const auto *sii0_724 = buffer.data(sii0 + 724);
    const auto *sii0_725 = buffer.data(sii0 + 725);
    const auto *sii0_727 = buffer.data(sii0 + 727);
    const auto *sii0_731 = buffer.data(sii0 + 731);
    const auto *sii0_734 = buffer.data(sii0 + 734);
    const auto *sii0_738 = buffer.data(sii0 + 738);
    const auto *sii0_740 = buffer.data(sii0 + 740);
    const auto *sii0_749 = buffer.data(sii0 + 749);
    const auto *sii0_751 = buffer.data(sii0 + 751);
    const auto *sii0_752 = buffer.data(sii0 + 752);
    const auto *sii0_753 = buffer.data(sii0 + 753);
    const auto *sii0_755 = buffer.data(sii0 + 755);
    const auto *sii0_756 = buffer.data(sii0 + 756);
    const auto *sii0_759 = buffer.data(sii0 + 759);
    const auto *sii0_761 = buffer.data(sii0 + 761);
    const auto *sii0_762 = buffer.data(sii0 + 762);
    const auto *sii0_765 = buffer.data(sii0 + 765);
    const auto *sii0_766 = buffer.data(sii0 + 766);
    const auto *sii0_768 = buffer.data(sii0 + 768);
    const auto *sii0_770 = buffer.data(sii0 + 770);
    const auto *sii0_777 = buffer.data(sii0 + 777);
    const auto *sii0_779 = buffer.data(sii0 + 779);
    const auto *sii0_780 = buffer.data(sii0 + 780);
    const auto *sii0_781 = buffer.data(sii0 + 781);
    const auto *sii0_783 = buffer.data(sii0 + 783);

    const auto *sih_378 = buffer.data(sih + 378);
    const auto *sih_381 = buffer.data(sih + 381);
    const auto *sih_384 = buffer.data(sih + 384);
    const auto *sih_393 = buffer.data(sih + 393);
    const auto *sih_398 = buffer.data(sih + 398);
    const auto *sih_399 = buffer.data(sih + 399);
    const auto *sih_401 = buffer.data(sih + 401);
    const auto *sih_402 = buffer.data(sih + 402);
    const auto *sih_404 = buffer.data(sih + 404);
    const auto *sih_405 = buffer.data(sih + 405);
    const auto *sih_408 = buffer.data(sih + 408);
    const auto *sih_414 = buffer.data(sih + 414);
    const auto *sih_419 = buffer.data(sih + 419);
    const auto *sih_420 = buffer.data(sih + 420);
    const auto *sih_422 = buffer.data(sih + 422);
    const auto *sih_423 = buffer.data(sih + 423);
    const auto *sih_425 = buffer.data(sih + 425);
    const auto *sih_426 = buffer.data(sih + 426);
    const auto *sih_429 = buffer.data(sih + 429);
    const auto *sih_435 = buffer.data(sih + 435);
    const auto *sih_440 = buffer.data(sih + 440);
    const auto *sih_441 = buffer.data(sih + 441);
    const auto *sih_443 = buffer.data(sih + 443);
    const auto *sih_444 = buffer.data(sih + 444);
    const auto *sih_446 = buffer.data(sih + 446);
    const auto *sih_450 = buffer.data(sih + 450);
    const auto *sih_456 = buffer.data(sih + 456);
    const auto *sih_458 = buffer.data(sih + 458);
    const auto *sih_459 = buffer.data(sih + 459);
    const auto *sih_460 = buffer.data(sih + 460);
    const auto *sih_461 = buffer.data(sih + 461);
    const auto *sih_462 = buffer.data(sih + 462);
    const auto *sih_464 = buffer.data(sih + 464);
    const auto *sih_467 = buffer.data(sih + 467);
    const auto *sih_525 = buffer.data(sih + 525);
    const auto *sih_528 = buffer.data(sih + 528);
    const auto *sih_530 = buffer.data(sih + 530);
    const auto *sih_531 = buffer.data(sih + 531);
    const auto *sih_534 = buffer.data(sih + 534);
    const auto *sih_535 = buffer.data(sih + 535);
    const auto *sih_537 = buffer.data(sih + 537);
    const auto *sih_539 = buffer.data(sih + 539);
    const auto *sih_540 = buffer.data(sih + 540);
    const auto *sih_541 = buffer.data(sih + 541);
    const auto *sih_542 = buffer.data(sih + 542);
    const auto *sih_543 = buffer.data(sih + 543);
    const auto *sih_544 = buffer.data(sih + 544);
    const auto *sih_545 = buffer.data(sih + 545);
    const auto *sih_549 = buffer.data(sih + 549);
    const auto *sih_552 = buffer.data(sih + 552);
    const auto *sih_556 = buffer.data(sih + 556);
    const auto *sih_558 = buffer.data(sih + 558);
    const auto *sih_561 = buffer.data(sih + 561);
    const auto *sih_562 = buffer.data(sih + 562);
    const auto *sih_563 = buffer.data(sih + 563);
    const auto *sih_564 = buffer.data(sih + 564);
    const auto *sih_565 = buffer.data(sih + 565);
    const auto *sih_566 = buffer.data(sih + 566);
    const auto *sih_567 = buffer.data(sih + 567);
    const auto *sih_570 = buffer.data(sih + 570);
    const auto *sih_572 = buffer.data(sih + 572);
    const auto *sih_573 = buffer.data(sih + 573);
    const auto *sih_576 = buffer.data(sih + 576);
    const auto *sih_577 = buffer.data(sih + 577);
    const auto *sih_579 = buffer.data(sih + 579);
    const auto *sih_581 = buffer.data(sih + 581);
    const auto *sih_582 = buffer.data(sih + 582);
    const auto *sih_583 = buffer.data(sih + 583);
    const auto *sih_584 = buffer.data(sih + 584);
    const auto *sih_585 = buffer.data(sih + 585);
    const auto *sih_586 = buffer.data(sih + 586);
    const auto *sih_587 = buffer.data(sih + 587);

    const auto *sii1_560 = buffer.data(sii1 + 560);
    const auto *sii1_565 = buffer.data(sii1 + 565);
    const auto *sii1_569 = buffer.data(sii1 + 569);
    const auto *sii1_574 = buffer.data(sii1 + 574);
    const auto *sii1_588 = buffer.data(sii1 + 588);
    const auto *sii1_591 = buffer.data(sii1 + 591);
    const auto *sii1_594 = buffer.data(sii1 + 594);
    const auto *sii1_699 = buffer.data(sii1 + 699);
    const auto *sii1_700 = buffer.data(sii1 + 700);
    const auto *sii1_703 = buffer.data(sii1 + 703);
    const auto *sii1_705 = buffer.data(sii1 + 705);
    const auto *sii1_706 = buffer.data(sii1 + 706);
    const auto *sii1_709 = buffer.data(sii1 + 709);
    const auto *sii1_710 = buffer.data(sii1 + 710);
    const auto *sii1_712 = buffer.data(sii1 + 712);
    const auto *sii1_714 = buffer.data(sii1 + 714);
    const auto *sii1_721 = buffer.data(sii1 + 721);
    const auto *sii1_723 = buffer.data(sii1 + 723);
    const auto *sii1_724 = buffer.data(sii1 + 724);
    const auto *sii1_725 = buffer.data(sii1 + 725);
    const auto *sii1_727 = buffer.data(sii1 + 727);
    const auto *sii1_731 = buffer.data(sii1 + 731);
    const auto *sii1_734 = buffer.data(sii1 + 734);
    const auto *sii1_738 = buffer.data(sii1 + 738);
    const auto *sii1_740 = buffer.data(sii1 + 740);
    const auto *sii1_749 = buffer.data(sii1 + 749);
    const auto *sii1_751 = buffer.data(sii1 + 751);
    const auto *sii1_752 = buffer.data(sii1 + 752);
    const auto *sii1_753 = buffer.data(sii1 + 753);
    const auto *sii1_755 = buffer.data(sii1 + 755);
    const auto *sii1_756 = buffer.data(sii1 + 756);
    const auto *sii1_759 = buffer.data(sii1 + 759);
    const auto *sii1_761 = buffer.data(sii1 + 761);
    const auto *sii1_762 = buffer.data(sii1 + 762);
    const auto *sii1_765 = buffer.data(sii1 + 765);
    const auto *sii1_766 = buffer.data(sii1 + 766);
    const auto *sii1_768 = buffer.data(sii1 + 768);
    const auto *sii1_770 = buffer.data(sii1 + 770);
    const auto *sii1_777 = buffer.data(sii1 + 777);
    const auto *sii1_779 = buffer.data(sii1 + 779);
    const auto *sii1_780 = buffer.data(sii1 + 780);
    const auto *sii1_781 = buffer.data(sii1 + 781);
    const auto *sii1_783 = buffer.data(sii1 + 783);

    const auto *skg0_420 = buffer.data(skg0 + 420);
    const auto *skg0_423 = buffer.data(skg0 + 423);
    const auto *skg0_425 = buffer.data(skg0 + 425);
    const auto *skg0_426 = buffer.data(skg0 + 426);
    const auto *skg0_429 = buffer.data(skg0 + 429);
    const auto *skg0_430 = buffer.data(skg0 + 430);
    const auto *skg0_432 = buffer.data(skg0 + 432);
    const auto *skg0_433 = buffer.data(skg0 + 433);
    const auto *skg0_434 = buffer.data(skg0 + 434);
    const auto *skg0_440 = buffer.data(skg0 + 440);
    const auto *skg0_444 = buffer.data(skg0 + 444);

    const auto *skg1_420 = buffer.data(skg1 + 420);
    const auto *skg1_423 = buffer.data(skg1 + 423);
    const auto *skg1_425 = buffer.data(skg1 + 425);
    const auto *skg1_426 = buffer.data(skg1 + 426);
    const auto *skg1_429 = buffer.data(skg1 + 429);
    const auto *skg1_430 = buffer.data(skg1 + 430);
    const auto *skg1_432 = buffer.data(skg1 + 432);
    const auto *skg1_433 = buffer.data(skg1 + 433);
    const auto *skg1_434 = buffer.data(skg1 + 434);
    const auto *skg1_440 = buffer.data(skg1 + 440);
    const auto *skg1_444 = buffer.data(skg1 + 444);

    const auto *skh_524 = buffer.data(skh + 524);
    const auto *skh_525 = buffer.data(skh + 525);
    const auto *skh_527 = buffer.data(skh + 527);
    const auto *skh_528 = buffer.data(skh + 528);
    const auto *skh_530 = buffer.data(skh + 530);
    const auto *skh_531 = buffer.data(skh + 531);
    const auto *skh_534 = buffer.data(skh + 534);
    const auto *skh_540 = buffer.data(skh + 540);
    const auto *skh_541 = buffer.data(skh + 541);
    const auto *skh_542 = buffer.data(skh + 542);
    const auto *skh_543 = buffer.data(skh + 543);
    const auto *skh_544 = buffer.data(skh + 544);
    const auto *skh_545 = buffer.data(skh + 545);
    const auto *skh_546 = buffer.data(skh + 546);
    const auto *skh_548 = buffer.data(skh + 548);
    const auto *skh_549 = buffer.data(skh + 549);
    const auto *skh_551 = buffer.data(skh + 551);
    const auto *skh_552 = buffer.data(skh + 552);
    const auto *skh_555 = buffer.data(skh + 555);
    const auto *skh_561 = buffer.data(skh + 561);
    const auto *skh_562 = buffer.data(skh + 562);
    const auto *skh_563 = buffer.data(skh + 563);
    const auto *skh_564 = buffer.data(skh + 564);
    const auto *skh_565 = buffer.data(skh + 565);
    const auto *skh_566 = buffer.data(skh + 566);
    const auto *skh_567 = buffer.data(skh + 567);
    const auto *skh_569 = buffer.data(skh + 569);
    const auto *skh_570 = buffer.data(skh + 570);
    const auto *skh_572 = buffer.data(skh + 572);
    const auto *skh_573 = buffer.data(skh + 573);
    const auto *skh_576 = buffer.data(skh + 576);
    const auto *skh_582 = buffer.data(skh + 582);
    const auto *skh_583 = buffer.data(skh + 583);
    const auto *skh_584 = buffer.data(skh + 584);
    const auto *skh_585 = buffer.data(skh + 585);
    const auto *skh_586 = buffer.data(skh + 586);
    const auto *skh_587 = buffer.data(skh + 587);
    const auto *skh_588 = buffer.data(skh + 588);
    const auto *skh_590 = buffer.data(skh + 590);
    const auto *skh_591 = buffer.data(skh + 591);
    const auto *skh_593 = buffer.data(skh + 593);
    const auto *skh_594 = buffer.data(skh + 594);
    const auto *skh_597 = buffer.data(skh + 597);
    const auto *skh_598 = buffer.data(skh + 598);
    const auto *skh_600 = buffer.data(skh + 600);
    const auto *skh_602 = buffer.data(skh + 602);
    const auto *skh_603 = buffer.data(skh + 603);
    const auto *skh_604 = buffer.data(skh + 604);
    const auto *skh_605 = buffer.data(skh + 605);
    const auto *skh_606 = buffer.data(skh + 606);
    const auto *skh_607 = buffer.data(skh + 607);
    const auto *skh_608 = buffer.data(skh + 608);
    const auto *skh_609 = buffer.data(skh + 609);
    const auto *skh_611 = buffer.data(skh + 611);
    const auto *skh_612 = buffer.data(skh + 612);
    const auto *skh_614 = buffer.data(skh + 614);
    const auto *skh_618 = buffer.data(skh + 618);

#pragma omp simd aligned(t_698, t_699, t_700, t_701, pb_x, pc_x, pc_y, sii0_699, sii0_700, \
                         sih_398, sih_399, sih_525, sii1_699, sii1_700, skh_524, \
                         skh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_13 * sih_398[k]
                   + f_3 * pc_y[k] * skh_524[k];

        t_699[k] = pb_x[k] * sii0_699[k]
                   - f_10 * pc_x[k] * sii1_699[k];

        t_700[k] = pb_x[k] * sii0_700[k]
                   + f_15 * sih_525[k]
                   - f_10 * pc_x[k] * sii1_700[k];

        t_701[k] = f_12 * sih_399[k]
                   + f_3 * pc_y[k] * skh_525[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, pb_x, pc_x, pc_y, pc_z, sii0_703, sih_378, \
                         sih_401, sih_528, sii1_703, skh_525, skh_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = f_14 * sih_378[k]
                   + f_3 * pc_z[k] * skh_525[k];

        t_703[k] = pb_x[k] * sii0_703[k]
                   + f_14 * sih_528[k]
                   - f_10 * pc_x[k] * sii1_703[k];

        t_704[k] = f_12 * sih_401[k]
                   + f_3 * pc_y[k] * skh_527[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, pb_x, pc_x, pc_z, sii0_705, sii0_706, sih_381, \
                         sih_530, sih_531, sii1_705, sii1_706, \
                         skh_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = pb_x[k] * sii0_705[k]
                   + f_14 * sih_530[k]
                   - f_10 * pc_x[k] * sii1_705[k];

        t_706[k] = pb_x[k] * sii0_706[k]
                   + f_13 * sih_531[k]
                   - f_10 * pc_x[k] * sii1_706[k];

        t_707[k] = f_14 * sih_381[k]
                   + f_3 * pc_z[k] * skh_528[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, pb_x, pc_x, pc_y, sii0_709, sii0_710, sih_404, \
                         sih_534, sih_535, sii1_709, sii1_710, \
                         skh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_12 * sih_404[k]
                   + f_3 * pc_y[k] * skh_530[k];

        t_709[k] = pb_x[k] * sii0_709[k]
                   + f_13 * sih_534[k]
                   - f_10 * pc_x[k] * sii1_709[k];

        t_710[k] = pb_x[k] * sii0_710[k]
                   + f_12 * sih_535[k]
                   - f_10 * pc_x[k] * sii1_710[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, pb_x, pc_x, pc_y, pc_z, sii0_712, sih_384, \
                         sih_408, sih_537, sii1_712, skh_531, skh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_14 * sih_384[k]
                   + f_3 * pc_z[k] * skh_531[k];

        t_712[k] = pb_x[k] * sii0_712[k]
                   + f_12 * sih_537[k]
                   - f_10 * pc_x[k] * sii1_712[k];

        t_713[k] = f_12 * sih_408[k]
                   + f_3 * pc_y[k] * skh_534[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pb_x, pc_x, sii0_714, sih_539, sih_540, \
                         sih_541, sih_542, sii1_714, skh_540, skh_541, \
                         skh_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pb_x[k] * sii0_714[k]
                   + f_12 * sih_539[k]
                   - f_10 * pc_x[k] * sii1_714[k];

        t_715[k] = f_11 * sih_540[k]
                   + f_3 * pc_x[k] * skh_540[k];

        t_716[k] = f_11 * sih_541[k]
                   + f_3 * pc_x[k] * skh_541[k];

        t_717[k] = f_11 * sih_542[k]
                   + f_3 * pc_x[k] * skh_542[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pb_x, pc_x, sii0_721, sih_543, sih_544, \
                         sih_545, sii1_721, skh_543, skh_544, skh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_11 * sih_543[k]
                   + f_3 * pc_x[k] * skh_543[k];

        t_719[k] = f_11 * sih_544[k]
                   + f_3 * pc_x[k] * skh_544[k];

        t_720[k] = f_11 * sih_545[k]
                   + f_3 * pc_x[k] * skh_545[k];

        t_721[k] = pb_x[k] * sii0_721[k]
                   - f_10 * pc_x[k] * sii1_721[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, pb_x, pc_x, pc_z, sii0_723, sii0_724, \
                         sii0_725, sih_393, sii1_723, sii1_724, sii1_725, \
                         skh_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_14 * sih_393[k]
                   + f_3 * pc_z[k] * skh_540[k];

        t_723[k] = pb_x[k] * sii0_723[k]
                   - f_10 * pc_x[k] * sii1_723[k];

        t_724[k] = pb_x[k] * sii0_724[k]
                   - f_10 * pc_x[k] * sii1_724[k];

        t_725[k] = pb_x[k] * sii0_725[k]
                   - f_10 * pc_x[k] * sii1_725[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pb_x, pb_y, pc_x, pc_y, sii0_560, \
                         sii0_727, sih_419, sih_420, sii1_560, sii1_727, skh_545, \
                         skh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * sih_419[k]
                   + f_3 * pc_y[k] * skh_545[k];

        t_727[k] = pb_x[k] * sii0_727[k]
                   - f_10 * pc_x[k] * sii1_727[k];

        t_728[k] = pb_y[k] * sii0_560[k]
                   - f_10 * pc_y[k] * sii1_560[k];

        t_729[k] = f_11 * sih_420[k]
                   + f_3 * pc_y[k] * skh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pb_x, pc_x, pc_y, pc_z, sii0_731, sih_399, \
                         sih_422, sih_549, sii1_731, skh_546, skh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_16 * sih_399[k]
                   + f_3 * pc_z[k] * skh_546[k];

        t_731[k] = pb_x[k] * sii0_731[k]
                   + f_14 * sih_549[k]
                   - f_10 * pc_x[k] * sii1_731[k];

        t_732[k] = f_11 * sih_422[k]
                   + f_3 * pc_y[k] * skh_548[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pb_x, pb_y, pc_x, pc_y, pc_z, sii0_565, \
                         sii0_734, sih_402, sih_552, sii1_565, sii1_734, \
                         skh_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = pb_y[k] * sii0_565[k]
                   - f_10 * pc_y[k] * sii1_565[k];

        t_734[k] = pb_x[k] * sii0_734[k]
                   + f_13 * sih_552[k]
                   - f_10 * pc_x[k] * sii1_734[k];

        t_735[k] = f_16 * sih_402[k]
                   + f_3 * pc_z[k] * skh_549[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pb_x, pb_y, pc_x, pc_y, sii0_569, sii0_738, \
                         sih_425, sih_556, sii1_569, sii1_738, \
                         skh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_11 * sih_425[k]
                   + f_3 * pc_y[k] * skh_551[k];

        t_737[k] = pb_y[k] * sii0_569[k]
                   - f_10 * pc_y[k] * sii1_569[k];

        t_738[k] = pb_x[k] * sii0_738[k]
                   + f_12 * sih_556[k]
                   - f_10 * pc_x[k] * sii1_738[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pb_x, pc_x, pc_y, pc_z, sii0_740, sih_405, \
                         sih_429, sih_558, sii1_740, skh_552, skh_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_16 * sih_405[k]
                   + f_3 * pc_z[k] * skh_552[k];

        t_740[k] = pb_x[k] * sii0_740[k]
                   + f_12 * sih_558[k]
                   - f_10 * pc_x[k] * sii1_740[k];

        t_741[k] = f_11 * sih_429[k]
                   + f_3 * pc_y[k] * skh_555[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pb_y, pc_x, pc_y, sii0_574, sih_561, \
                         sih_562, sih_563, sii1_574, skh_561, skh_562, \
                         skh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = pb_y[k] * sii0_574[k]
                   - f_10 * pc_y[k] * sii1_574[k];

        t_743[k] = f_11 * sih_561[k]
                   + f_3 * pc_x[k] * skh_561[k];

        t_744[k] = f_11 * sih_562[k]
                   + f_3 * pc_x[k] * skh_562[k];

        t_745[k] = f_11 * sih_563[k]
                   + f_3 * pc_x[k] * skh_563[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pb_x, pc_x, sii0_749, sih_564, sih_565, \
                         sih_566, sii1_749, skh_564, skh_565, skh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_11 * sih_564[k]
                   + f_3 * pc_x[k] * skh_564[k];

        t_747[k] = f_11 * sih_565[k]
                   + f_3 * pc_x[k] * skh_565[k];

        t_748[k] = f_11 * sih_566[k]
                   + f_3 * pc_x[k] * skh_566[k];

        t_749[k] = pb_x[k] * sii0_749[k]
                   - f_10 * pc_x[k] * sii1_749[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pb_x, pc_x, pc_z, sii0_751, sii0_752, \
                         sii0_753, sih_414, sii1_751, sii1_752, sii1_753, \
                         skh_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_16 * sih_414[k]
                   + f_3 * pc_z[k] * skh_561[k];

        t_751[k] = pb_x[k] * sii0_751[k]
                   - f_10 * pc_x[k] * sii1_751[k];

        t_752[k] = pb_x[k] * sii0_752[k]
                   - f_10 * pc_x[k] * sii1_752[k];

        t_753[k] = pb_x[k] * sii0_753[k]
                   - f_10 * pc_x[k] * sii1_753[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, pb_x, pc_x, pc_y, sii0_755, sii0_756, \
                         sih_440, sih_567, sii1_755, sii1_756, skh_566, \
                         skh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_11 * sih_440[k]
                   + f_3 * pc_y[k] * skh_566[k];

        t_755[k] = pb_x[k] * sii0_755[k]
                   - f_10 * pc_x[k] * sii1_755[k];

        t_756[k] = pb_x[k] * sii0_756[k]
                   + f_15 * sih_567[k]
                   - f_10 * pc_x[k] * sii1_756[k];

        t_757[k] = f_3 * pc_y[k] * skh_567[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pb_x, pc_x, pc_y, pc_z, sii0_759, sih_420, \
                         sih_570, sii1_759, skh_567, skh_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_15 * sih_420[k]
                   + f_3 * pc_z[k] * skh_567[k];

        t_759[k] = pb_x[k] * sii0_759[k]
                   + f_14 * sih_570[k]
                   - f_10 * pc_x[k] * sii1_759[k];

        t_760[k] = f_3 * pc_y[k] * skh_569[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pb_x, pc_x, pc_z, sii0_761, sii0_762, sih_423, \
                         sih_572, sih_573, sii1_761, sii1_762, \
                         skh_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = pb_x[k] * sii0_761[k]
                   + f_14 * sih_572[k]
                   - f_10 * pc_x[k] * sii1_761[k];

        t_762[k] = pb_x[k] * sii0_762[k]
                   + f_13 * sih_573[k]
                   - f_10 * pc_x[k] * sii1_762[k];

        t_763[k] = f_15 * sih_423[k]
                   + f_3 * pc_z[k] * skh_570[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pb_x, pc_x, pc_y, sii0_765, sii0_766, sih_576, \
                         sih_577, sii1_765, sii1_766, skh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_3 * pc_y[k] * skh_572[k];

        t_765[k] = pb_x[k] * sii0_765[k]
                   + f_13 * sih_576[k]
                   - f_10 * pc_x[k] * sii1_765[k];

        t_766[k] = pb_x[k] * sii0_766[k]
                   + f_12 * sih_577[k]
                   - f_10 * pc_x[k] * sii1_766[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pb_x, pc_x, pc_y, pc_z, sii0_768, sih_426, \
                         sih_579, sii1_768, skh_573, skh_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_15 * sih_426[k]
                   + f_3 * pc_z[k] * skh_573[k];

        t_768[k] = pb_x[k] * sii0_768[k]
                   + f_12 * sih_579[k]
                   - f_10 * pc_x[k] * sii1_768[k];

        t_769[k] = f_3 * pc_y[k] * skh_576[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, pb_x, pc_x, sii0_770, sih_581, sih_582, \
                         sih_583, sih_584, sii1_770, skh_582, skh_583, \
                         skh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = pb_x[k] * sii0_770[k]
                   + f_12 * sih_581[k]
                   - f_10 * pc_x[k] * sii1_770[k];

        t_771[k] = f_11 * sih_582[k]
                   + f_3 * pc_x[k] * skh_582[k];

        t_772[k] = f_11 * sih_583[k]
                   + f_3 * pc_x[k] * skh_583[k];

        t_773[k] = f_11 * sih_584[k]
                   + f_3 * pc_x[k] * skh_584[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pb_x, pc_x, sii0_777, sih_585, sih_586, \
                         sih_587, sii1_777, skh_585, skh_586, skh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_11 * sih_585[k]
                   + f_3 * pc_x[k] * skh_585[k];

        t_775[k] = f_11 * sih_586[k]
                   + f_3 * pc_x[k] * skh_586[k];

        t_776[k] = f_11 * sih_587[k]
                   + f_3 * pc_x[k] * skh_587[k];

        t_777[k] = pb_x[k] * sii0_777[k]
                   - f_10 * pc_x[k] * sii1_777[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, pb_x, pc_x, pc_z, sii0_779, sii0_780, \
                         sii0_781, sih_435, sii1_779, sii1_780, sii1_781, \
                         skh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_15 * sih_435[k]
                   + f_3 * pc_z[k] * skh_582[k];

        t_779[k] = pb_x[k] * sii0_779[k]
                   - f_10 * pc_x[k] * sii1_779[k];

        t_780[k] = pb_x[k] * sii0_780[k]
                   - f_10 * pc_x[k] * sii1_780[k];

        t_781[k] = pb_x[k] * sii0_781[k]
                   - f_10 * pc_x[k] * sii1_781[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, t_786, pb_x, pc_x, pc_y, pc_z, sii0_783, \
                         sih_441, sii1_783, skg0_420, skg1_420, skh_587, \
                         skh_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_3 * pc_y[k] * skh_587[k];

        t_783[k] = pb_x[k] * sii0_783[k]
                   - f_10 * pc_x[k] * sii1_783[k];

        t_784[k] = f_1 * skg0_420[k]
                   - f_2 * skg1_420[k]
                   + f_3 * pc_x[k] * skh_588[k];

        t_785[k] = f_0 * sih_441[k]
                   + f_3 * pc_y[k] * skh_588[k];

        t_786[k] = f_3 * pc_z[k] * skh_588[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_x, pc_y, sih_443, skg0_423, skg0_425, \
                         skg1_423, skg1_425, skh_590, skh_591, \
                         skh_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_4 * skg0_423[k]
                   - f_5 * skg1_423[k]
                   + f_3 * pc_x[k] * skh_591[k];

        t_788[k] = f_0 * sih_443[k]
                   + f_3 * pc_y[k] * skh_590[k];

        t_789[k] = f_4 * skg0_425[k]
                   - f_5 * skg1_425[k]
                   + f_3 * pc_x[k] * skh_593[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pc_x, pc_y, pc_z, sih_446, skg0_426, \
                         skg0_429, skg1_426, skg1_429, skh_591, skh_593, skh_594, \
                         skh_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_6 * skg0_426[k]
                   - f_7 * skg1_426[k]
                   + f_3 * pc_x[k] * skh_594[k];

        t_791[k] = f_3 * pc_z[k] * skh_591[k];

        t_792[k] = f_0 * sih_446[k]
                   + f_3 * pc_y[k] * skh_593[k];

        t_793[k] = f_6 * skg0_429[k]
                   - f_7 * skg1_429[k]
                   + f_3 * pc_x[k] * skh_597[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, pc_x, pc_y, pc_z, sih_450, skg0_430, \
                         skg0_432, skg1_430, skg1_432, skh_594, skh_597, skh_598, \
                         skh_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_8 * skg0_430[k]
                   - f_9 * skg1_430[k]
                   + f_3 * pc_x[k] * skh_598[k];

        t_795[k] = f_3 * pc_z[k] * skh_594[k];

        t_796[k] = f_8 * skg0_432[k]
                   - f_9 * skg1_432[k]
                   + f_3 * pc_x[k] * skh_600[k];

        t_797[k] = f_0 * sih_450[k]
                   + f_3 * pc_y[k] * skh_597[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, t_802, t_803, pc_x, skg0_434, skg1_434, \
                         skh_602, skh_603, skh_604, skh_605, skh_606, \
                         skh_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = f_8 * skg0_434[k]
                   - f_9 * skg1_434[k]
                   + f_3 * pc_x[k] * skh_602[k];

        t_799[k] = f_3 * pc_x[k] * skh_603[k];

        t_800[k] = f_3 * pc_x[k] * skh_604[k];

        t_801[k] = f_3 * pc_x[k] * skh_605[k];

        t_802[k] = f_3 * pc_x[k] * skh_606[k];

        t_803[k] = f_3 * pc_x[k] * skh_607[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, pc_x, pc_y, pc_z, sih_456, sih_458, \
                         skg0_430, skg0_432, skg1_430, skg1_432, skh_603, skh_605, \
                         skh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_3 * pc_x[k] * skh_608[k];

        t_805[k] = f_0 * sih_456[k]
                   + f_1 * skg0_430[k]
                   - f_2 * skg1_430[k]
                   + f_3 * pc_y[k] * skh_603[k];

        t_806[k] = f_3 * pc_z[k] * skh_603[k];

        t_807[k] = f_0 * sih_458[k]
                   + f_4 * skg0_432[k]
                   - f_5 * skg1_432[k]
                   + f_3 * pc_y[k] * skh_605[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pc_y, pc_z, sih_459, sih_460, sih_461, \
                         skg0_433, skg0_434, skg1_433, skg1_434, skh_606, skh_607, \
                         skh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_0 * sih_459[k]
                   + f_6 * skg0_433[k]
                   - f_7 * skg1_433[k]
                   + f_3 * pc_y[k] * skh_606[k];

        t_809[k] = f_0 * sih_460[k]
                   + f_8 * skg0_434[k]
                   - f_9 * skg1_434[k]
                   + f_3 * pc_y[k] * skh_607[k];

        t_810[k] = f_0 * sih_461[k]
                   + f_3 * pc_y[k] * skh_608[k];

        t_811[k] = f_1 * skg0_434[k]
                   - f_2 * skg1_434[k]
                   + f_3 * pc_z[k] * skh_608[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pb_z, pc_y, pc_z, sii0_588, sii0_591, \
                         sih_441, sih_462, sii1_588, sii1_591, \
                         skh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = pb_z[k] * sii0_588[k]
                   - f_10 * pc_z[k] * sii1_588[k];

        t_813[k] = f_15 * sih_462[k]
                   + f_3 * pc_y[k] * skh_609[k];

        t_814[k] = f_11 * sih_441[k]
                   + f_3 * pc_z[k] * skh_609[k];

        t_815[k] = pb_z[k] * sii0_591[k]
                   - f_10 * pc_z[k] * sii1_591[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pb_z, pc_x, pc_y, pc_z, sii0_594, sih_464, \
                         sii1_594, skg0_440, skg1_440, skh_611, \
                         skh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_15 * sih_464[k]
                   + f_3 * pc_y[k] * skh_611[k];

        t_817[k] = f_4 * skg0_440[k]
                   - f_5 * skg1_440[k]
                   + f_3 * pc_x[k] * skh_614[k];

        t_818[k] = pb_z[k] * sii0_594[k]
                   - f_10 * pc_z[k] * sii1_594[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, pc_y, pc_z, sih_444, sih_467, skg0_444, \
                         skg1_444, skh_612, skh_614, skh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_11 * sih_444[k]
                   + f_3 * pc_z[k] * skh_612[k];

        t_820[k] = f_15 * sih_467[k]
                   + f_3 * pc_y[k] * skh_614[k];

        t_821[k] = f_6 * skg0_444[k]
                   - f_7 * skg1_444[k]
                   + f_3 * pc_x[k] * skh_618[k];
    }
}

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_598 = buffer.data(sii0 + 598);
    const auto *sii0_609 = buffer.data(sii0 + 609);
    const auto *sii0_611 = buffer.data(sii0 + 611);
    const auto *sii0_612 = buffer.data(sii0 + 612);
    const auto *sii0_613 = buffer.data(sii0 + 613);

    const auto *sih_447 = buffer.data(sih + 447);
    const auto *sih_456 = buffer.data(sih + 456);
    const auto *sih_457 = buffer.data(sih + 457);
    const auto *sih_458 = buffer.data(sih + 458);
    const auto *sih_459 = buffer.data(sih + 459);
    const auto *sih_461 = buffer.data(sih + 461);
    const auto *sih_462 = buffer.data(sih + 462);
    const auto *sih_465 = buffer.data(sih + 465);
    const auto *sih_468 = buffer.data(sih + 468);
    const auto *sih_471 = buffer.data(sih + 471);
    const auto *sih_477 = buffer.data(sih + 477);
    const auto *sih_482 = buffer.data(sih + 482);
    const auto *sih_483 = buffer.data(sih + 483);
    const auto *sih_485 = buffer.data(sih + 485);
    const auto *sih_486 = buffer.data(sih + 486);
    const auto *sih_488 = buffer.data(sih + 488);
    const auto *sih_489 = buffer.data(sih + 489);
    const auto *sih_492 = buffer.data(sih + 492);
    const auto *sih_498 = buffer.data(sih + 498);
    const auto *sih_500 = buffer.data(sih + 500);
    const auto *sih_501 = buffer.data(sih + 501);
    const auto *sih_502 = buffer.data(sih + 502);
    const auto *sih_503 = buffer.data(sih + 503);
    const auto *sih_504 = buffer.data(sih + 504);
    const auto *sih_506 = buffer.data(sih + 506);
    const auto *sih_507 = buffer.data(sih + 507);
    const auto *sih_509 = buffer.data(sih + 509);
    const auto *sih_510 = buffer.data(sih + 510);
    const auto *sih_513 = buffer.data(sih + 513);
    const auto *sih_519 = buffer.data(sih + 519);
    const auto *sih_521 = buffer.data(sih + 521);
    const auto *sih_522 = buffer.data(sih + 522);
    const auto *sih_523 = buffer.data(sih + 523);
    const auto *sih_524 = buffer.data(sih + 524);
    const auto *sih_525 = buffer.data(sih + 525);
    const auto *sih_527 = buffer.data(sih + 527);
    const auto *sih_528 = buffer.data(sih + 528);
    const auto *sih_530 = buffer.data(sih + 530);
    const auto *sih_531 = buffer.data(sih + 531);
    const auto *sih_534 = buffer.data(sih + 534);
    const auto *sih_540 = buffer.data(sih + 540);
    const auto *sih_542 = buffer.data(sih + 542);
    const auto *sih_543 = buffer.data(sih + 543);
    const auto *sih_544 = buffer.data(sih + 544);
    const auto *sih_545 = buffer.data(sih + 545);
    const auto *sih_546 = buffer.data(sih + 546);
    const auto *sih_548 = buffer.data(sih + 548);
    const auto *sih_551 = buffer.data(sih + 551);
    const auto *sih_555 = buffer.data(sih + 555);

    const auto *sii1_598 = buffer.data(sii1 + 598);
    const auto *sii1_609 = buffer.data(sii1 + 609);
    const auto *sii1_611 = buffer.data(sii1 + 611);
    const auto *sii1_612 = buffer.data(sii1 + 612);
    const auto *sii1_613 = buffer.data(sii1 + 613);

    const auto *skg0_447 = buffer.data(skg0 + 447);
    const auto *skg0_449 = buffer.data(skg0 + 449);
    const auto *skg0_450 = buffer.data(skg0 + 450);
    const auto *skg0_453 = buffer.data(skg0 + 453);
    const auto *skg0_455 = buffer.data(skg0 + 455);
    const auto *skg0_456 = buffer.data(skg0 + 456);
    const auto *skg0_459 = buffer.data(skg0 + 459);
    const auto *skg0_460 = buffer.data(skg0 + 460);
    const auto *skg0_462 = buffer.data(skg0 + 462);
    const auto *skg0_463 = buffer.data(skg0 + 463);
    const auto *skg0_464 = buffer.data(skg0 + 464);
    const auto *skg0_465 = buffer.data(skg0 + 465);
    const auto *skg0_468 = buffer.data(skg0 + 468);
    const auto *skg0_470 = buffer.data(skg0 + 470);
    const auto *skg0_471 = buffer.data(skg0 + 471);
    const auto *skg0_474 = buffer.data(skg0 + 474);
    const auto *skg0_475 = buffer.data(skg0 + 475);
    const auto *skg0_477 = buffer.data(skg0 + 477);
    const auto *skg0_478 = buffer.data(skg0 + 478);
    const auto *skg0_479 = buffer.data(skg0 + 479);
    const auto *skg0_480 = buffer.data(skg0 + 480);
    const auto *skg0_483 = buffer.data(skg0 + 483);
    const auto *skg0_485 = buffer.data(skg0 + 485);
    const auto *skg0_486 = buffer.data(skg0 + 486);
    const auto *skg0_489 = buffer.data(skg0 + 489);
    const auto *skg0_490 = buffer.data(skg0 + 490);
    const auto *skg0_492 = buffer.data(skg0 + 492);
    const auto *skg0_493 = buffer.data(skg0 + 493);
    const auto *skg0_494 = buffer.data(skg0 + 494);
    const auto *skg0_495 = buffer.data(skg0 + 495);
    const auto *skg0_498 = buffer.data(skg0 + 498);
    const auto *skg0_500 = buffer.data(skg0 + 500);
    const auto *skg0_501 = buffer.data(skg0 + 501);
    const auto *skg0_504 = buffer.data(skg0 + 504);
    const auto *skg0_505 = buffer.data(skg0 + 505);
    const auto *skg0_507 = buffer.data(skg0 + 507);
    const auto *skg0_509 = buffer.data(skg0 + 509);

    const auto *skg1_447 = buffer.data(skg1 + 447);
    const auto *skg1_449 = buffer.data(skg1 + 449);
    const auto *skg1_450 = buffer.data(skg1 + 450);
    const auto *skg1_453 = buffer.data(skg1 + 453);
    const auto *skg1_455 = buffer.data(skg1 + 455);
    const auto *skg1_456 = buffer.data(skg1 + 456);
    const auto *skg1_459 = buffer.data(skg1 + 459);
    const auto *skg1_460 = buffer.data(skg1 + 460);
    const auto *skg1_462 = buffer.data(skg1 + 462);
    const auto *skg1_463 = buffer.data(skg1 + 463);
    const auto *skg1_464 = buffer.data(skg1 + 464);
    const auto *skg1_465 = buffer.data(skg1 + 465);
    const auto *skg1_468 = buffer.data(skg1 + 468);
    const auto *skg1_470 = buffer.data(skg1 + 470);
    const auto *skg1_471 = buffer.data(skg1 + 471);
    const auto *skg1_474 = buffer.data(skg1 + 474);
    const auto *skg1_475 = buffer.data(skg1 + 475);
    const auto *skg1_477 = buffer.data(skg1 + 477);
    const auto *skg1_478 = buffer.data(skg1 + 478);
    const auto *skg1_479 = buffer.data(skg1 + 479);
    const auto *skg1_480 = buffer.data(skg1 + 480);
    const auto *skg1_483 = buffer.data(skg1 + 483);
    const auto *skg1_485 = buffer.data(skg1 + 485);
    const auto *skg1_486 = buffer.data(skg1 + 486);
    const auto *skg1_489 = buffer.data(skg1 + 489);
    const auto *skg1_490 = buffer.data(skg1 + 490);
    const auto *skg1_492 = buffer.data(skg1 + 492);
    const auto *skg1_493 = buffer.data(skg1 + 493);
    const auto *skg1_494 = buffer.data(skg1 + 494);
    const auto *skg1_495 = buffer.data(skg1 + 495);
    const auto *skg1_498 = buffer.data(skg1 + 498);
    const auto *skg1_500 = buffer.data(skg1 + 500);
    const auto *skg1_501 = buffer.data(skg1 + 501);
    const auto *skg1_504 = buffer.data(skg1 + 504);
    const auto *skg1_505 = buffer.data(skg1 + 505);
    const auto *skg1_507 = buffer.data(skg1 + 507);
    const auto *skg1_509 = buffer.data(skg1 + 509);

    const auto *skh_615 = buffer.data(skh + 615);
    const auto *skh_618 = buffer.data(skh + 618);
    const auto *skh_621 = buffer.data(skh + 621);
    const auto *skh_623 = buffer.data(skh + 623);
    const auto *skh_624 = buffer.data(skh + 624);
    const auto *skh_625 = buffer.data(skh + 625);
    const auto *skh_626 = buffer.data(skh + 626);
    const auto *skh_627 = buffer.data(skh + 627);
    const auto *skh_628 = buffer.data(skh + 628);
    const auto *skh_629 = buffer.data(skh + 629);
    const auto *skh_630 = buffer.data(skh + 630);
    const auto *skh_632 = buffer.data(skh + 632);
    const auto *skh_633 = buffer.data(skh + 633);
    const auto *skh_635 = buffer.data(skh + 635);
    const auto *skh_636 = buffer.data(skh + 636);
    const auto *skh_639 = buffer.data(skh + 639);
    const auto *skh_640 = buffer.data(skh + 640);
    const auto *skh_642 = buffer.data(skh + 642);
    const auto *skh_644 = buffer.data(skh + 644);
    const auto *skh_645 = buffer.data(skh + 645);
    const auto *skh_646 = buffer.data(skh + 646);
    const auto *skh_647 = buffer.data(skh + 647);
    const auto *skh_648 = buffer.data(skh + 648);
    const auto *skh_649 = buffer.data(skh + 649);
    const auto *skh_650 = buffer.data(skh + 650);
    const auto *skh_651 = buffer.data(skh + 651);
    const auto *skh_653 = buffer.data(skh + 653);
    const auto *skh_654 = buffer.data(skh + 654);
    const auto *skh_656 = buffer.data(skh + 656);
    const auto *skh_657 = buffer.data(skh + 657);
    const auto *skh_660 = buffer.data(skh + 660);
    const auto *skh_661 = buffer.data(skh + 661);
    const auto *skh_663 = buffer.data(skh + 663);
    const auto *skh_665 = buffer.data(skh + 665);
    const auto *skh_666 = buffer.data(skh + 666);
    const auto *skh_667 = buffer.data(skh + 667);
    const auto *skh_668 = buffer.data(skh + 668);
    const auto *skh_669 = buffer.data(skh + 669);
    const auto *skh_670 = buffer.data(skh + 670);
    const auto *skh_671 = buffer.data(skh + 671);
    const auto *skh_672 = buffer.data(skh + 672);
    const auto *skh_674 = buffer.data(skh + 674);
    const auto *skh_675 = buffer.data(skh + 675);
    const auto *skh_677 = buffer.data(skh + 677);
    const auto *skh_678 = buffer.data(skh + 678);
    const auto *skh_681 = buffer.data(skh + 681);
    const auto *skh_682 = buffer.data(skh + 682);
    const auto *skh_684 = buffer.data(skh + 684);
    const auto *skh_686 = buffer.data(skh + 686);
    const auto *skh_687 = buffer.data(skh + 687);
    const auto *skh_688 = buffer.data(skh + 688);
    const auto *skh_689 = buffer.data(skh + 689);
    const auto *skh_690 = buffer.data(skh + 690);
    const auto *skh_691 = buffer.data(skh + 691);
    const auto *skh_692 = buffer.data(skh + 692);
    const auto *skh_693 = buffer.data(skh + 693);
    const auto *skh_695 = buffer.data(skh + 695);
    const auto *skh_696 = buffer.data(skh + 696);
    const auto *skh_698 = buffer.data(skh + 698);
    const auto *skh_699 = buffer.data(skh + 699);
    const auto *skh_702 = buffer.data(skh + 702);
    const auto *skh_703 = buffer.data(skh + 703);
    const auto *skh_705 = buffer.data(skh + 705);
    const auto *skh_707 = buffer.data(skh + 707);
    const auto *skh_708 = buffer.data(skh + 708);
    const auto *skh_709 = buffer.data(skh + 709);
    const auto *skh_710 = buffer.data(skh + 710);
    const auto *skh_711 = buffer.data(skh + 711);
    const auto *skh_712 = buffer.data(skh + 712);
    const auto *skh_713 = buffer.data(skh + 713);

#pragma omp simd aligned(t_822, t_823, t_824, pb_z, pc_x, pc_z, sii0_598, sih_447, sii1_598, \
                         skg0_447, skg1_447, skh_615, skh_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pb_z[k] * sii0_598[k]
                   - f_10 * pc_z[k] * sii1_598[k];

        t_823[k] = f_11 * sih_447[k]
                   + f_3 * pc_z[k] * skh_615[k];

        t_824[k] = f_8 * skg0_447[k]
                   - f_9 * skg1_447[k]
                   + f_3 * pc_x[k] * skh_621[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, pc_x, pc_y, sih_471, skg0_449, \
                         skg1_449, skh_618, skh_623, skh_624, skh_625, \
                         skh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = f_15 * sih_471[k]
                   + f_3 * pc_y[k] * skh_618[k];

        t_826[k] = f_8 * skg0_449[k]
                   - f_9 * skg1_449[k]
                   + f_3 * pc_x[k] * skh_623[k];

        t_827[k] = f_3 * pc_x[k] * skh_624[k];

        t_828[k] = f_3 * pc_x[k] * skh_625[k];

        t_829[k] = f_3 * pc_x[k] * skh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, pb_z, pc_x, pc_z, sii0_609, \
                         sih_456, sii1_609, skh_624, skh_627, skh_628, \
                         skh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_3 * pc_x[k] * skh_627[k];

        t_831[k] = f_3 * pc_x[k] * skh_628[k];

        t_832[k] = f_3 * pc_x[k] * skh_629[k];

        t_833[k] = pb_z[k] * sii0_609[k]
                   - f_10 * pc_z[k] * sii1_609[k];

        t_834[k] = f_11 * sih_456[k]
                   + f_3 * pc_z[k] * skh_624[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, pb_z, pc_z, sii0_611, sii0_612, sii0_613, \
                         sih_457, sih_458, sih_459, sii1_611, sii1_612, \
                         sii1_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = pb_z[k] * sii0_611[k]
                   + f_12 * sih_457[k]
                   - f_10 * pc_z[k] * sii1_611[k];

        t_836[k] = pb_z[k] * sii0_612[k]
                   + f_13 * sih_458[k]
                   - f_10 * pc_z[k] * sii1_612[k];

        t_837[k] = pb_z[k] * sii0_613[k]
                   + f_14 * sih_459[k]
                   - f_10 * pc_z[k] * sii1_613[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, t_841, pc_x, pc_y, pc_z, sih_461, sih_482, \
                         sih_483, skg0_449, skg0_450, skg1_449, skg1_450, skh_629, \
                         skh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_15 * sih_482[k]
                   + f_3 * pc_y[k] * skh_629[k];

        t_839[k] = f_11 * sih_461[k]
                   + f_1 * skg0_449[k]
                   - f_2 * skg1_449[k]
                   + f_3 * pc_z[k] * skh_629[k];

        t_840[k] = f_1 * skg0_450[k]
                   - f_2 * skg1_450[k]
                   + f_3 * pc_x[k] * skh_630[k];

        t_841[k] = f_16 * sih_483[k]
                   + f_3 * pc_y[k] * skh_630[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, pc_y, pc_z, sih_462, sih_485, skg0_453, \
                         skg1_453, skh_630, skh_632, skh_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_12 * sih_462[k]
                   + f_3 * pc_z[k] * skh_630[k];

        t_843[k] = f_4 * skg0_453[k]
                   - f_5 * skg1_453[k]
                   + f_3 * pc_x[k] * skh_633[k];

        t_844[k] = f_16 * sih_485[k]
                   + f_3 * pc_y[k] * skh_632[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, pc_x, pc_y, pc_z, sih_465, sih_488, \
                         skg0_455, skg0_456, skg1_455, skg1_456, skh_633, skh_635, \
                         skh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_4 * skg0_455[k]
                   - f_5 * skg1_455[k]
                   + f_3 * pc_x[k] * skh_635[k];

        t_846[k] = f_6 * skg0_456[k]
                   - f_7 * skg1_456[k]
                   + f_3 * pc_x[k] * skh_636[k];

        t_847[k] = f_12 * sih_465[k]
                   + f_3 * pc_z[k] * skh_633[k];

        t_848[k] = f_16 * sih_488[k]
                   + f_3 * pc_y[k] * skh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, sih_468, skg0_459, skg0_460, \
                         skg1_459, skg1_460, skh_636, skh_639, \
                         skh_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_6 * skg0_459[k]
                   - f_7 * skg1_459[k]
                   + f_3 * pc_x[k] * skh_639[k];

        t_850[k] = f_8 * skg0_460[k]
                   - f_9 * skg1_460[k]
                   + f_3 * pc_x[k] * skh_640[k];

        t_851[k] = f_12 * sih_468[k]
                   + f_3 * pc_z[k] * skh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, pc_y, sih_492, skg0_462, skg0_464, \
                         skg1_462, skg1_464, skh_639, skh_642, skh_644, \
                         skh_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_8 * skg0_462[k]
                   - f_9 * skg1_462[k]
                   + f_3 * pc_x[k] * skh_642[k];

        t_853[k] = f_16 * sih_492[k]
                   + f_3 * pc_y[k] * skh_639[k];

        t_854[k] = f_8 * skg0_464[k]
                   - f_9 * skg1_464[k]
                   + f_3 * pc_x[k] * skh_644[k];

        t_855[k] = f_3 * pc_x[k] * skh_645[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, t_859, t_860, pc_x, skh_646, skh_647, skh_648, \
                         skh_649, skh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_3 * pc_x[k] * skh_646[k];

        t_857[k] = f_3 * pc_x[k] * skh_647[k];

        t_858[k] = f_3 * pc_x[k] * skh_648[k];

        t_859[k] = f_3 * pc_x[k] * skh_649[k];

        t_860[k] = f_3 * pc_x[k] * skh_650[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, pc_y, pc_z, sih_477, sih_498, sih_500, skg0_460, \
                         skg0_462, skg1_460, skg1_462, skh_645, \
                         skh_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = f_16 * sih_498[k]
                   + f_1 * skg0_460[k]
                   - f_2 * skg1_460[k]
                   + f_3 * pc_y[k] * skh_645[k];

        t_862[k] = f_12 * sih_477[k]
                   + f_3 * pc_z[k] * skh_645[k];

        t_863[k] = f_16 * sih_500[k]
                   + f_4 * skg0_462[k]
                   - f_5 * skg1_462[k]
                   + f_3 * pc_y[k] * skh_647[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, pc_y, sih_501, sih_502, sih_503, skg0_463, \
                         skg0_464, skg1_463, skg1_464, skh_648, skh_649, \
                         skh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = f_16 * sih_501[k]
                   + f_6 * skg0_463[k]
                   - f_7 * skg1_463[k]
                   + f_3 * pc_y[k] * skh_648[k];

        t_865[k] = f_16 * sih_502[k]
                   + f_8 * skg0_464[k]
                   - f_9 * skg1_464[k]
                   + f_3 * pc_y[k] * skh_649[k];

        t_866[k] = f_16 * sih_503[k]
                   + f_3 * pc_y[k] * skh_650[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, t_870, pc_x, pc_y, pc_z, sih_482, sih_483, \
                         sih_504, skg0_464, skg0_465, skg1_464, skg1_465, skh_650, \
                         skh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = f_12 * sih_482[k]
                   + f_1 * skg0_464[k]
                   - f_2 * skg1_464[k]
                   + f_3 * pc_z[k] * skh_650[k];

        t_868[k] = f_1 * skg0_465[k]
                   - f_2 * skg1_465[k]
                   + f_3 * pc_x[k] * skh_651[k];

        t_869[k] = f_14 * sih_504[k]
                   + f_3 * pc_y[k] * skh_651[k];

        t_870[k] = f_13 * sih_483[k]
                   + f_3 * pc_z[k] * skh_651[k];
    }

#pragma omp simd aligned(t_871, t_872, t_873, pc_x, pc_y, sih_506, skg0_468, skg0_470, \
                         skg1_468, skg1_470, skh_653, skh_654, \
                         skh_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_871[k] = f_4 * skg0_468[k]
                   - f_5 * skg1_468[k]
                   + f_3 * pc_x[k] * skh_654[k];

        t_872[k] = f_14 * sih_506[k]
                   + f_3 * pc_y[k] * skh_653[k];

        t_873[k] = f_4 * skg0_470[k]
                   - f_5 * skg1_470[k]
                   + f_3 * pc_x[k] * skh_656[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, pc_x, pc_y, pc_z, sih_486, sih_509, skg0_471, \
                         skg1_471, skh_654, skh_656, skh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_6 * skg0_471[k]
                   - f_7 * skg1_471[k]
                   + f_3 * pc_x[k] * skh_657[k];

        t_875[k] = f_13 * sih_486[k]
                   + f_3 * pc_z[k] * skh_654[k];

        t_876[k] = f_14 * sih_509[k]
                   + f_3 * pc_y[k] * skh_656[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_x, pc_z, sih_489, skg0_474, skg0_475, \
                         skg1_474, skg1_475, skh_657, skh_660, \
                         skh_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_6 * skg0_474[k]
                   - f_7 * skg1_474[k]
                   + f_3 * pc_x[k] * skh_660[k];

        t_878[k] = f_8 * skg0_475[k]
                   - f_9 * skg1_475[k]
                   + f_3 * pc_x[k] * skh_661[k];

        t_879[k] = f_13 * sih_489[k]
                   + f_3 * pc_z[k] * skh_657[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pc_x, pc_y, sih_513, skg0_477, skg0_479, \
                         skg1_477, skg1_479, skh_660, skh_663, skh_665, \
                         skh_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_8 * skg0_477[k]
                   - f_9 * skg1_477[k]
                   + f_3 * pc_x[k] * skh_663[k];

        t_881[k] = f_14 * sih_513[k]
                   + f_3 * pc_y[k] * skh_660[k];

        t_882[k] = f_8 * skg0_479[k]
                   - f_9 * skg1_479[k]
                   + f_3 * pc_x[k] * skh_665[k];

        t_883[k] = f_3 * pc_x[k] * skh_666[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, pc_x, skh_667, skh_668, skh_669, \
                         skh_670, skh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_3 * pc_x[k] * skh_667[k];

        t_885[k] = f_3 * pc_x[k] * skh_668[k];

        t_886[k] = f_3 * pc_x[k] * skh_669[k];

        t_887[k] = f_3 * pc_x[k] * skh_670[k];

        t_888[k] = f_3 * pc_x[k] * skh_671[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, pc_y, pc_z, sih_498, sih_519, sih_521, skg0_475, \
                         skg0_477, skg1_475, skg1_477, skh_666, \
                         skh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_14 * sih_519[k]
                   + f_1 * skg0_475[k]
                   - f_2 * skg1_475[k]
                   + f_3 * pc_y[k] * skh_666[k];

        t_890[k] = f_13 * sih_498[k]
                   + f_3 * pc_z[k] * skh_666[k];

        t_891[k] = f_14 * sih_521[k]
                   + f_4 * skg0_477[k]
                   - f_5 * skg1_477[k]
                   + f_3 * pc_y[k] * skh_668[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pc_y, sih_522, sih_523, sih_524, skg0_478, \
                         skg0_479, skg1_478, skg1_479, skh_669, skh_670, \
                         skh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_14 * sih_522[k]
                   + f_6 * skg0_478[k]
                   - f_7 * skg1_478[k]
                   + f_3 * pc_y[k] * skh_669[k];

        t_893[k] = f_14 * sih_523[k]
                   + f_8 * skg0_479[k]
                   - f_9 * skg1_479[k]
                   + f_3 * pc_y[k] * skh_670[k];

        t_894[k] = f_14 * sih_524[k]
                   + f_3 * pc_y[k] * skh_671[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, pc_x, pc_y, pc_z, sih_503, sih_504, \
                         sih_525, skg0_479, skg0_480, skg1_479, skg1_480, skh_671, \
                         skh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_13 * sih_503[k]
                   + f_1 * skg0_479[k]
                   - f_2 * skg1_479[k]
                   + f_3 * pc_z[k] * skh_671[k];

        t_896[k] = f_1 * skg0_480[k]
                   - f_2 * skg1_480[k]
                   + f_3 * pc_x[k] * skh_672[k];

        t_897[k] = f_13 * sih_525[k]
                   + f_3 * pc_y[k] * skh_672[k];

        t_898[k] = f_14 * sih_504[k]
                   + f_3 * pc_z[k] * skh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, sih_527, skg0_483, skg0_485, \
                         skg1_483, skg1_485, skh_674, skh_675, \
                         skh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_4 * skg0_483[k]
                   - f_5 * skg1_483[k]
                   + f_3 * pc_x[k] * skh_675[k];

        t_900[k] = f_13 * sih_527[k]
                   + f_3 * pc_y[k] * skh_674[k];

        t_901[k] = f_4 * skg0_485[k]
                   - f_5 * skg1_485[k]
                   + f_3 * pc_x[k] * skh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, sih_507, sih_530, skg0_486, \
                         skg1_486, skh_675, skh_677, skh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_6 * skg0_486[k]
                   - f_7 * skg1_486[k]
                   + f_3 * pc_x[k] * skh_678[k];

        t_903[k] = f_14 * sih_507[k]
                   + f_3 * pc_z[k] * skh_675[k];

        t_904[k] = f_13 * sih_530[k]
                   + f_3 * pc_y[k] * skh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, sih_510, skg0_489, skg0_490, \
                         skg1_489, skg1_490, skh_678, skh_681, \
                         skh_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_6 * skg0_489[k]
                   - f_7 * skg1_489[k]
                   + f_3 * pc_x[k] * skh_681[k];

        t_906[k] = f_8 * skg0_490[k]
                   - f_9 * skg1_490[k]
                   + f_3 * pc_x[k] * skh_682[k];

        t_907[k] = f_14 * sih_510[k]
                   + f_3 * pc_z[k] * skh_678[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, t_911, pc_x, pc_y, sih_534, skg0_492, skg0_494, \
                         skg1_492, skg1_494, skh_681, skh_684, skh_686, \
                         skh_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_8 * skg0_492[k]
                   - f_9 * skg1_492[k]
                   + f_3 * pc_x[k] * skh_684[k];

        t_909[k] = f_13 * sih_534[k]
                   + f_3 * pc_y[k] * skh_681[k];

        t_910[k] = f_8 * skg0_494[k]
                   - f_9 * skg1_494[k]
                   + f_3 * pc_x[k] * skh_686[k];

        t_911[k] = f_3 * pc_x[k] * skh_687[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, t_915, t_916, pc_x, skh_688, skh_689, skh_690, \
                         skh_691, skh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_3 * pc_x[k] * skh_688[k];

        t_913[k] = f_3 * pc_x[k] * skh_689[k];

        t_914[k] = f_3 * pc_x[k] * skh_690[k];

        t_915[k] = f_3 * pc_x[k] * skh_691[k];

        t_916[k] = f_3 * pc_x[k] * skh_692[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pc_y, pc_z, sih_519, sih_540, sih_542, skg0_490, \
                         skg0_492, skg1_490, skg1_492, skh_687, \
                         skh_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_13 * sih_540[k]
                   + f_1 * skg0_490[k]
                   - f_2 * skg1_490[k]
                   + f_3 * pc_y[k] * skh_687[k];

        t_918[k] = f_14 * sih_519[k]
                   + f_3 * pc_z[k] * skh_687[k];

        t_919[k] = f_13 * sih_542[k]
                   + f_4 * skg0_492[k]
                   - f_5 * skg1_492[k]
                   + f_3 * pc_y[k] * skh_689[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, pc_y, sih_543, sih_544, sih_545, skg0_493, \
                         skg0_494, skg1_493, skg1_494, skh_690, skh_691, \
                         skh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_13 * sih_543[k]
                   + f_6 * skg0_493[k]
                   - f_7 * skg1_493[k]
                   + f_3 * pc_y[k] * skh_690[k];

        t_921[k] = f_13 * sih_544[k]
                   + f_8 * skg0_494[k]
                   - f_9 * skg1_494[k]
                   + f_3 * pc_y[k] * skh_691[k];

        t_922[k] = f_13 * sih_545[k]
                   + f_3 * pc_y[k] * skh_692[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, t_926, pc_x, pc_y, pc_z, sih_524, sih_525, \
                         sih_546, skg0_494, skg0_495, skg1_494, skg1_495, skh_692, \
                         skh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = f_14 * sih_524[k]
                   + f_1 * skg0_494[k]
                   - f_2 * skg1_494[k]
                   + f_3 * pc_z[k] * skh_692[k];

        t_924[k] = f_1 * skg0_495[k]
                   - f_2 * skg1_495[k]
                   + f_3 * pc_x[k] * skh_693[k];

        t_925[k] = f_12 * sih_546[k]
                   + f_3 * pc_y[k] * skh_693[k];

        t_926[k] = f_16 * sih_525[k]
                   + f_3 * pc_z[k] * skh_693[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, pc_x, pc_y, sih_548, skg0_498, skg0_500, \
                         skg1_498, skg1_500, skh_695, skh_696, \
                         skh_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_4 * skg0_498[k]
                   - f_5 * skg1_498[k]
                   + f_3 * pc_x[k] * skh_696[k];

        t_928[k] = f_12 * sih_548[k]
                   + f_3 * pc_y[k] * skh_695[k];

        t_929[k] = f_4 * skg0_500[k]
                   - f_5 * skg1_500[k]
                   + f_3 * pc_x[k] * skh_698[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, pc_x, pc_y, pc_z, sih_528, sih_551, skg0_501, \
                         skg1_501, skh_696, skh_698, skh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_6 * skg0_501[k]
                   - f_7 * skg1_501[k]
                   + f_3 * pc_x[k] * skh_699[k];

        t_931[k] = f_16 * sih_528[k]
                   + f_3 * pc_z[k] * skh_696[k];

        t_932[k] = f_12 * sih_551[k]
                   + f_3 * pc_y[k] * skh_698[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, pc_x, pc_z, sih_531, skg0_504, skg0_505, \
                         skg1_504, skg1_505, skh_699, skh_702, \
                         skh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_6 * skg0_504[k]
                   - f_7 * skg1_504[k]
                   + f_3 * pc_x[k] * skh_702[k];

        t_934[k] = f_8 * skg0_505[k]
                   - f_9 * skg1_505[k]
                   + f_3 * pc_x[k] * skh_703[k];

        t_935[k] = f_16 * sih_531[k]
                   + f_3 * pc_z[k] * skh_699[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, pc_x, pc_y, sih_555, skg0_507, skg0_509, \
                         skg1_507, skg1_509, skh_702, skh_705, skh_707, \
                         skh_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_8 * skg0_507[k]
                   - f_9 * skg1_507[k]
                   + f_3 * pc_x[k] * skh_705[k];

        t_937[k] = f_12 * sih_555[k]
                   + f_3 * pc_y[k] * skh_702[k];

        t_938[k] = f_8 * skg0_509[k]
                   - f_9 * skg1_509[k]
                   + f_3 * pc_x[k] * skh_707[k];

        t_939[k] = f_3 * pc_x[k] * skh_708[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, pc_x, skh_709, skh_710, skh_711, \
                         skh_712, skh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_3 * pc_x[k] * skh_709[k];

        t_941[k] = f_3 * pc_x[k] * skh_710[k];

        t_942[k] = f_3 * pc_x[k] * skh_711[k];

        t_943[k] = f_3 * pc_x[k] * skh_712[k];

        t_944[k] = f_3 * pc_x[k] * skh_713[k];
    }
}

static auto
compute_prim_ski_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sii0,
                                                          const size_t sih, const size_t sii1,
                                                          const size_t skg0, const size_t skg1,
                                                          const size_t skh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sii0_756 = buffer.data(sii0 + 756);
    const auto *sii0_761 = buffer.data(sii0 + 761);
    const auto *sii0_765 = buffer.data(sii0 + 765);
    const auto *sii0_770 = buffer.data(sii0 + 770);
    const auto *sii0_777 = buffer.data(sii0 + 777);
    const auto *sii0_779 = buffer.data(sii0 + 779);
    const auto *sii0_780 = buffer.data(sii0 + 780);
    const auto *sii0_781 = buffer.data(sii0 + 781);
    const auto *sii0_783 = buffer.data(sii0 + 783);

    const auto *sih_540 = buffer.data(sih + 540);
    const auto *sih_545 = buffer.data(sih + 545);
    const auto *sih_546 = buffer.data(sih + 546);
    const auto *sih_549 = buffer.data(sih + 549);
    const auto *sih_552 = buffer.data(sih + 552);
    const auto *sih_561 = buffer.data(sih + 561);
    const auto *sih_563 = buffer.data(sih + 563);
    const auto *sih_564 = buffer.data(sih + 564);
    const auto *sih_565 = buffer.data(sih + 565);
    const auto *sih_566 = buffer.data(sih + 566);
    const auto *sih_567 = buffer.data(sih + 567);
    const auto *sih_569 = buffer.data(sih + 569);
    const auto *sih_570 = buffer.data(sih + 570);
    const auto *sih_572 = buffer.data(sih + 572);
    const auto *sih_573 = buffer.data(sih + 573);
    const auto *sih_576 = buffer.data(sih + 576);
    const auto *sih_582 = buffer.data(sih + 582);
    const auto *sih_584 = buffer.data(sih + 584);
    const auto *sih_585 = buffer.data(sih + 585);
    const auto *sih_586 = buffer.data(sih + 586);
    const auto *sih_587 = buffer.data(sih + 587);

    const auto *sii1_756 = buffer.data(sii1 + 756);
    const auto *sii1_761 = buffer.data(sii1 + 761);
    const auto *sii1_765 = buffer.data(sii1 + 765);
    const auto *sii1_770 = buffer.data(sii1 + 770);
    const auto *sii1_777 = buffer.data(sii1 + 777);
    const auto *sii1_779 = buffer.data(sii1 + 779);
    const auto *sii1_780 = buffer.data(sii1 + 780);
    const auto *sii1_781 = buffer.data(sii1 + 781);
    const auto *sii1_783 = buffer.data(sii1 + 783);

    const auto *skg0_505 = buffer.data(skg0 + 505);
    const auto *skg0_507 = buffer.data(skg0 + 507);
    const auto *skg0_508 = buffer.data(skg0 + 508);
    const auto *skg0_509 = buffer.data(skg0 + 509);
    const auto *skg0_513 = buffer.data(skg0 + 513);
    const auto *skg0_516 = buffer.data(skg0 + 516);
    const auto *skg0_520 = buffer.data(skg0 + 520);
    const auto *skg0_522 = buffer.data(skg0 + 522);
    const auto *skg0_525 = buffer.data(skg0 + 525);
    const auto *skg0_528 = buffer.data(skg0 + 528);
    const auto *skg0_530 = buffer.data(skg0 + 530);
    const auto *skg0_531 = buffer.data(skg0 + 531);
    const auto *skg0_534 = buffer.data(skg0 + 534);
    const auto *skg0_535 = buffer.data(skg0 + 535);
    const auto *skg0_537 = buffer.data(skg0 + 537);
    const auto *skg0_538 = buffer.data(skg0 + 538);
    const auto *skg0_539 = buffer.data(skg0 + 539);

    const auto *skg1_505 = buffer.data(skg1 + 505);
    const auto *skg1_507 = buffer.data(skg1 + 507);
    const auto *skg1_508 = buffer.data(skg1 + 508);
    const auto *skg1_509 = buffer.data(skg1 + 509);
    const auto *skg1_513 = buffer.data(skg1 + 513);
    const auto *skg1_516 = buffer.data(skg1 + 516);
    const auto *skg1_520 = buffer.data(skg1 + 520);
    const auto *skg1_522 = buffer.data(skg1 + 522);
    const auto *skg1_525 = buffer.data(skg1 + 525);
    const auto *skg1_528 = buffer.data(skg1 + 528);
    const auto *skg1_530 = buffer.data(skg1 + 530);
    const auto *skg1_531 = buffer.data(skg1 + 531);
    const auto *skg1_534 = buffer.data(skg1 + 534);
    const auto *skg1_535 = buffer.data(skg1 + 535);
    const auto *skg1_537 = buffer.data(skg1 + 537);
    const auto *skg1_538 = buffer.data(skg1 + 538);
    const auto *skg1_539 = buffer.data(skg1 + 539);

    const auto *skh_708 = buffer.data(skh + 708);
    const auto *skh_710 = buffer.data(skh + 710);
    const auto *skh_711 = buffer.data(skh + 711);
    const auto *skh_712 = buffer.data(skh + 712);
    const auto *skh_713 = buffer.data(skh + 713);
    const auto *skh_714 = buffer.data(skh + 714);
    const auto *skh_716 = buffer.data(skh + 716);
    const auto *skh_717 = buffer.data(skh + 717);
    const auto *skh_719 = buffer.data(skh + 719);
    const auto *skh_720 = buffer.data(skh + 720);
    const auto *skh_723 = buffer.data(skh + 723);
    const auto *skh_724 = buffer.data(skh + 724);
    const auto *skh_726 = buffer.data(skh + 726);
    const auto *skh_729 = buffer.data(skh + 729);
    const auto *skh_730 = buffer.data(skh + 730);
    const auto *skh_731 = buffer.data(skh + 731);
    const auto *skh_732 = buffer.data(skh + 732);
    const auto *skh_733 = buffer.data(skh + 733);
    const auto *skh_734 = buffer.data(skh + 734);
    const auto *skh_735 = buffer.data(skh + 735);
    const auto *skh_737 = buffer.data(skh + 737);
    const auto *skh_738 = buffer.data(skh + 738);
    const auto *skh_740 = buffer.data(skh + 740);
    const auto *skh_741 = buffer.data(skh + 741);
    const auto *skh_744 = buffer.data(skh + 744);
    const auto *skh_745 = buffer.data(skh + 745);
    const auto *skh_747 = buffer.data(skh + 747);
    const auto *skh_749 = buffer.data(skh + 749);
    const auto *skh_750 = buffer.data(skh + 750);
    const auto *skh_751 = buffer.data(skh + 751);
    const auto *skh_752 = buffer.data(skh + 752);
    const auto *skh_753 = buffer.data(skh + 753);
    const auto *skh_754 = buffer.data(skh + 754);
    const auto *skh_755 = buffer.data(skh + 755);

#pragma omp simd aligned(t_945, t_946, t_947, pc_y, pc_z, sih_540, sih_561, sih_563, skg0_505, \
                         skg0_507, skg1_505, skg1_507, skh_708, \
                         skh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_12 * sih_561[k]
                   + f_1 * skg0_505[k]
                   - f_2 * skg1_505[k]
                   + f_3 * pc_y[k] * skh_708[k];

        t_946[k] = f_16 * sih_540[k]
                   + f_3 * pc_z[k] * skh_708[k];

        t_947[k] = f_12 * sih_563[k]
                   + f_4 * skg0_507[k]
                   - f_5 * skg1_507[k]
                   + f_3 * pc_y[k] * skh_710[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, pc_y, sih_564, sih_565, sih_566, skg0_508, \
                         skg0_509, skg1_508, skg1_509, skh_711, skh_712, \
                         skh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = f_12 * sih_564[k]
                   + f_6 * skg0_508[k]
                   - f_7 * skg1_508[k]
                   + f_3 * pc_y[k] * skh_711[k];

        t_949[k] = f_12 * sih_565[k]
                   + f_8 * skg0_509[k]
                   - f_9 * skg1_509[k]
                   + f_3 * pc_y[k] * skh_712[k];

        t_950[k] = f_12 * sih_566[k]
                   + f_3 * pc_y[k] * skh_713[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, t_954, pb_y, pc_y, pc_z, sii0_756, sih_545, \
                         sih_546, sih_567, sii1_756, skg0_509, skg1_509, skh_713, \
                         skh_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_16 * sih_545[k]
                   + f_1 * skg0_509[k]
                   - f_2 * skg1_509[k]
                   + f_3 * pc_z[k] * skh_713[k];

        t_952[k] = pb_y[k] * sii0_756[k]
                   - f_10 * pc_y[k] * sii1_756[k];

        t_953[k] = f_11 * sih_567[k]
                   + f_3 * pc_y[k] * skh_714[k];

        t_954[k] = f_15 * sih_546[k]
                   + f_3 * pc_z[k] * skh_714[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, pb_y, pc_x, pc_y, sii0_761, sih_569, sii1_761, \
                         skg0_513, skg1_513, skh_716, skh_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_4 * skg0_513[k]
                   - f_5 * skg1_513[k]
                   + f_3 * pc_x[k] * skh_717[k];

        t_956[k] = f_11 * sih_569[k]
                   + f_3 * pc_y[k] * skh_716[k];

        t_957[k] = pb_y[k] * sii0_761[k]
                   - f_10 * pc_y[k] * sii1_761[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, pc_x, pc_y, pc_z, sih_549, sih_572, skg0_516, \
                         skg1_516, skh_717, skh_719, skh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = f_6 * skg0_516[k]
                   - f_7 * skg1_516[k]
                   + f_3 * pc_x[k] * skh_720[k];

        t_959[k] = f_15 * sih_549[k]
                   + f_3 * pc_z[k] * skh_717[k];

        t_960[k] = f_11 * sih_572[k]
                   + f_3 * pc_y[k] * skh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pb_y, pc_x, pc_y, pc_z, sii0_765, sih_552, \
                         sii1_765, skg0_520, skg1_520, skh_720, \
                         skh_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pb_y[k] * sii0_765[k]
                   - f_10 * pc_y[k] * sii1_765[k];

        t_962[k] = f_8 * skg0_520[k]
                   - f_9 * skg1_520[k]
                   + f_3 * pc_x[k] * skh_724[k];

        t_963[k] = f_15 * sih_552[k]
                   + f_3 * pc_z[k] * skh_720[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pb_y, pc_x, pc_y, sii0_770, sih_576, \
                         sii1_770, skg0_522, skg1_522, skh_723, skh_726, \
                         skh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = f_8 * skg0_522[k]
                   - f_9 * skg1_522[k]
                   + f_3 * pc_x[k] * skh_726[k];

        t_965[k] = f_11 * sih_576[k]
                   + f_3 * pc_y[k] * skh_723[k];

        t_966[k] = pb_y[k] * sii0_770[k]
                   - f_10 * pc_y[k] * sii1_770[k];

        t_967[k] = f_3 * pc_x[k] * skh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, skh_730, skh_731, skh_732, \
                         skh_733, skh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_3 * pc_x[k] * skh_730[k];

        t_969[k] = f_3 * pc_x[k] * skh_731[k];

        t_970[k] = f_3 * pc_x[k] * skh_732[k];

        t_971[k] = f_3 * pc_x[k] * skh_733[k];

        t_972[k] = f_3 * pc_x[k] * skh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pb_y, pc_y, pc_z, sii0_777, sii0_779, sih_561, \
                         sih_582, sih_584, sii1_777, sii1_779, \
                         skh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = pb_y[k] * sii0_777[k]
                   + f_15 * sih_582[k]
                   - f_10 * pc_y[k] * sii1_777[k];

        t_974[k] = f_15 * sih_561[k]
                   + f_3 * pc_z[k] * skh_729[k];

        t_975[k] = pb_y[k] * sii0_779[k]
                   + f_14 * sih_584[k]
                   - f_10 * pc_y[k] * sii1_779[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, t_979, pb_y, pc_y, sii0_780, sii0_781, sii0_783, \
                         sih_585, sih_586, sih_587, sii1_780, sii1_781, sii1_783, \
                         skh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = pb_y[k] * sii0_780[k]
                   + f_13 * sih_585[k]
                   - f_10 * pc_y[k] * sii1_780[k];

        t_977[k] = pb_y[k] * sii0_781[k]
                   + f_12 * sih_586[k]
                   - f_10 * pc_y[k] * sii1_781[k];

        t_978[k] = f_11 * sih_587[k]
                   + f_3 * pc_y[k] * skh_734[k];

        t_979[k] = pb_y[k] * sii0_783[k]
                   - f_10 * pc_y[k] * sii1_783[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, t_984, pc_x, pc_y, pc_z, sih_567, \
                         skg0_525, skg0_528, skg1_525, skg1_528, skh_735, skh_737, \
                         skh_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_1 * skg0_525[k]
                   - f_2 * skg1_525[k]
                   + f_3 * pc_x[k] * skh_735[k];

        t_981[k] = f_3 * pc_y[k] * skh_735[k];

        t_982[k] = f_0 * sih_567[k]
                   + f_3 * pc_z[k] * skh_735[k];

        t_983[k] = f_4 * skg0_528[k]
                   - f_5 * skg1_528[k]
                   + f_3 * pc_x[k] * skh_738[k];

        t_984[k] = f_3 * pc_y[k] * skh_737[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, pc_x, pc_y, pc_z, sih_570, skg0_530, \
                         skg0_531, skg1_530, skg1_531, skh_738, skh_740, \
                         skh_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_4 * skg0_530[k]
                   - f_5 * skg1_530[k]
                   + f_3 * pc_x[k] * skh_740[k];

        t_986[k] = f_6 * skg0_531[k]
                   - f_7 * skg1_531[k]
                   + f_3 * pc_x[k] * skh_741[k];

        t_987[k] = f_0 * sih_570[k]
                   + f_3 * pc_z[k] * skh_738[k];

        t_988[k] = f_3 * pc_y[k] * skh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_z, sih_573, skg0_534, skg0_535, \
                         skg1_534, skg1_535, skh_741, skh_744, \
                         skh_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_6 * skg0_534[k]
                   - f_7 * skg1_534[k]
                   + f_3 * pc_x[k] * skh_744[k];

        t_990[k] = f_8 * skg0_535[k]
                   - f_9 * skg1_535[k]
                   + f_3 * pc_x[k] * skh_745[k];

        t_991[k] = f_0 * sih_573[k]
                   + f_3 * pc_z[k] * skh_741[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, t_996, pc_x, pc_y, skg0_537, skg0_539, \
                         skg1_537, skg1_539, skh_744, skh_747, skh_749, skh_750, \
                         skh_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_8 * skg0_537[k]
                   - f_9 * skg1_537[k]
                   + f_3 * pc_x[k] * skh_747[k];

        t_993[k] = f_3 * pc_y[k] * skh_744[k];

        t_994[k] = f_8 * skg0_539[k]
                   - f_9 * skg1_539[k]
                   + f_3 * pc_x[k] * skh_749[k];

        t_995[k] = f_3 * pc_x[k] * skh_750[k];

        t_996[k] = f_3 * pc_x[k] * skh_751[k];
    }

#pragma omp simd aligned(t_997, t_998, t_999, t_1000, t_1001, pc_x, pc_y, skg0_535, skg1_535, \
                         skh_750, skh_752, skh_753, skh_754, skh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = f_3 * pc_x[k] * skh_752[k];

        t_998[k] = f_3 * pc_x[k] * skh_753[k];

        t_999[k] = f_3 * pc_x[k] * skh_754[k];

        t_1000[k] = f_3 * pc_x[k] * skh_755[k];

        t_1001[k] = f_1 * skg0_535[k]
                    - f_2 * skg1_535[k]
                    + f_3 * pc_y[k] * skh_750[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, pc_z, sih_582, skg0_537, skg0_538, \
                         skg1_537, skg1_538, skh_750, skh_752, \
                         skh_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_0 * sih_582[k]
                    + f_3 * pc_z[k] * skh_750[k];

        t_1003[k] = f_4 * skg0_537[k]
                    - f_5 * skg1_537[k]
                    + f_3 * pc_y[k] * skh_752[k];

        t_1004[k] = f_6 * skg0_538[k]
                    - f_7 * skg1_538[k]
                    + f_3 * pc_y[k] * skh_753[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pc_y, pc_z, sih_587, skg0_539, skg1_539, \
                         skh_754, skh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_8 * skg0_539[k]
                    - f_9 * skg1_539[k]
                    + f_3 * pc_y[k] * skh_754[k];

        t_1006[k] = f_3 * pc_y[k] * skh_755[k];

        t_1007[k] = f_0 * sih_587[k]
                    + f_1 * skg0_539[k]
                    - f_2 * skg1_539[k]
                    + f_3 * pc_z[k] * skh_755[k];
    }
}

auto
compute_prim_ski_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sii0, const size_t sih,
                                                   const size_t sii1, const size_t skg0,
                                                   const size_t skg1, const size_t skh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ski_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);

    compute_prim_ski_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);

    compute_prim_ski_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);

    compute_prim_ski_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);

    compute_prim_ski_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);

    compute_prim_ski_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);

    compute_prim_ski_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);

    compute_prim_ski_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);

    compute_prim_ski_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, sii0, sih,
                                                              sii1, skg0, skg1, skh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
