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


#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_skh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sih0,
                                                          const size_t sig, const size_t sih1,
                                                          const size_t skf0, const size_t skf1,
                                                          const size_t skg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;

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
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sih0_0 = buffer.data(sih0 + 0);
    const auto *sih0_3 = buffer.data(sih0 + 3);
    const auto *sih0_5 = buffer.data(sih0 + 5);
    const auto *sih0_6 = buffer.data(sih0 + 6);
    const auto *sih0_9 = buffer.data(sih0 + 9);
    const auto *sih0_15 = buffer.data(sih0 + 15);
    const auto *sih0_20 = buffer.data(sih0 + 20);
    const auto *sih0_24 = buffer.data(sih0 + 24);
    const auto *sih0_27 = buffer.data(sih0 + 27);
    const auto *sih0_36 = buffer.data(sih0 + 36);
    const auto *sih0_42 = buffer.data(sih0 + 42);
    const auto *sih0_47 = buffer.data(sih0 + 47);
    const auto *sih0_51 = buffer.data(sih0 + 51);
    const auto *sih0_62 = buffer.data(sih0 + 62);

    const auto *sig_0 = buffer.data(sig + 0);
    const auto *sig_1 = buffer.data(sig + 1);
    const auto *sig_2 = buffer.data(sig + 2);
    const auto *sig_3 = buffer.data(sig + 3);
    const auto *sig_5 = buffer.data(sig + 5);
    const auto *sig_6 = buffer.data(sig + 6);
    const auto *sig_9 = buffer.data(sig + 9);
    const auto *sig_10 = buffer.data(sig + 10);
    const auto *sig_11 = buffer.data(sig + 11);
    const auto *sig_12 = buffer.data(sig + 12);
    const auto *sig_13 = buffer.data(sig + 13);
    const auto *sig_14 = buffer.data(sig + 14);
    const auto *sig_15 = buffer.data(sig + 15);
    const auto *sig_17 = buffer.data(sig + 17);
    const auto *sig_18 = buffer.data(sig + 18);
    const auto *sig_20 = buffer.data(sig + 20);
    const auto *sig_25 = buffer.data(sig + 25);
    const auto *sig_26 = buffer.data(sig + 26);
    const auto *sig_27 = buffer.data(sig + 27);
    const auto *sig_28 = buffer.data(sig + 28);
    const auto *sig_29 = buffer.data(sig + 29);
    const auto *sig_30 = buffer.data(sig + 30);
    const auto *sig_32 = buffer.data(sig + 32);
    const auto *sig_33 = buffer.data(sig + 33);
    const auto *sig_35 = buffer.data(sig + 35);
    const auto *sig_40 = buffer.data(sig + 40);
    const auto *sig_41 = buffer.data(sig + 41);
    const auto *sig_42 = buffer.data(sig + 42);
    const auto *sig_43 = buffer.data(sig + 43);
    const auto *sig_44 = buffer.data(sig + 44);
    const auto *sig_45 = buffer.data(sig + 45);
    const auto *sig_48 = buffer.data(sig + 48);
    const auto *sig_50 = buffer.data(sig + 50);
    const auto *sig_51 = buffer.data(sig + 51);
    const auto *sig_54 = buffer.data(sig + 54);
    const auto *sig_55 = buffer.data(sig + 55);
    const auto *sig_56 = buffer.data(sig + 56);
    const auto *sig_57 = buffer.data(sig + 57);
    const auto *sig_58 = buffer.data(sig + 58);
    const auto *sig_59 = buffer.data(sig + 59);
    const auto *sig_70 = buffer.data(sig + 70);
    const auto *sig_71 = buffer.data(sig + 71);
    const auto *sig_72 = buffer.data(sig + 72);
    const auto *sig_73 = buffer.data(sig + 73);
    const auto *sig_74 = buffer.data(sig + 74);
    const auto *sig_75 = buffer.data(sig + 75);
    const auto *sig_78 = buffer.data(sig + 78);
    const auto *sig_80 = buffer.data(sig + 80);
    const auto *sig_81 = buffer.data(sig + 81);
    const auto *sig_84 = buffer.data(sig + 84);
    const auto *sig_85 = buffer.data(sig + 85);
    const auto *sig_86 = buffer.data(sig + 86);
    const auto *sig_87 = buffer.data(sig + 87);
    const auto *sig_88 = buffer.data(sig + 88);
    const auto *sig_89 = buffer.data(sig + 89);

    const auto *sih1_0 = buffer.data(sih1 + 0);
    const auto *sih1_3 = buffer.data(sih1 + 3);
    const auto *sih1_5 = buffer.data(sih1 + 5);
    const auto *sih1_6 = buffer.data(sih1 + 6);
    const auto *sih1_9 = buffer.data(sih1 + 9);
    const auto *sih1_15 = buffer.data(sih1 + 15);
    const auto *sih1_20 = buffer.data(sih1 + 20);
    const auto *sih1_24 = buffer.data(sih1 + 24);
    const auto *sih1_27 = buffer.data(sih1 + 27);
    const auto *sih1_36 = buffer.data(sih1 + 36);
    const auto *sih1_42 = buffer.data(sih1 + 42);
    const auto *sih1_47 = buffer.data(sih1 + 47);
    const auto *sih1_51 = buffer.data(sih1 + 51);
    const auto *sih1_62 = buffer.data(sih1 + 62);

    const auto *skf0_0 = buffer.data(skf0 + 0);
    const auto *skf0_3 = buffer.data(skf0 + 3);
    const auto *skf0_5 = buffer.data(skf0 + 5);
    const auto *skf0_6 = buffer.data(skf0 + 6);
    const auto *skf0_8 = buffer.data(skf0 + 8);
    const auto *skf0_9 = buffer.data(skf0 + 9);
    const auto *skf0_16 = buffer.data(skf0 + 16);
    const auto *skf0_18 = buffer.data(skf0 + 18);
    const auto *skf0_19 = buffer.data(skf0 + 19);
    const auto *skf0_28 = buffer.data(skf0 + 28);
    const auto *skf0_29 = buffer.data(skf0 + 29);
    const auto *skf0_30 = buffer.data(skf0 + 30);
    const auto *skf0_33 = buffer.data(skf0 + 33);
    const auto *skf0_35 = buffer.data(skf0 + 35);
    const auto *skf0_36 = buffer.data(skf0 + 36);
    const auto *skf0_38 = buffer.data(skf0 + 38);
    const auto *skf0_39 = buffer.data(skf0 + 39);
    const auto *skf0_48 = buffer.data(skf0 + 48);
    const auto *skf0_49 = buffer.data(skf0 + 49);
    const auto *skf0_50 = buffer.data(skf0 + 50);
    const auto *skf0_53 = buffer.data(skf0 + 53);
    const auto *skf0_55 = buffer.data(skf0 + 55);
    const auto *skf0_56 = buffer.data(skf0 + 56);
    const auto *skf0_58 = buffer.data(skf0 + 58);
    const auto *skf0_59 = buffer.data(skf0 + 59);

    const auto *skf1_0 = buffer.data(skf1 + 0);
    const auto *skf1_3 = buffer.data(skf1 + 3);
    const auto *skf1_5 = buffer.data(skf1 + 5);
    const auto *skf1_6 = buffer.data(skf1 + 6);
    const auto *skf1_8 = buffer.data(skf1 + 8);
    const auto *skf1_9 = buffer.data(skf1 + 9);
    const auto *skf1_16 = buffer.data(skf1 + 16);
    const auto *skf1_18 = buffer.data(skf1 + 18);
    const auto *skf1_19 = buffer.data(skf1 + 19);
    const auto *skf1_28 = buffer.data(skf1 + 28);
    const auto *skf1_29 = buffer.data(skf1 + 29);
    const auto *skf1_30 = buffer.data(skf1 + 30);
    const auto *skf1_33 = buffer.data(skf1 + 33);
    const auto *skf1_35 = buffer.data(skf1 + 35);
    const auto *skf1_36 = buffer.data(skf1 + 36);
    const auto *skf1_38 = buffer.data(skf1 + 38);
    const auto *skf1_39 = buffer.data(skf1 + 39);
    const auto *skf1_48 = buffer.data(skf1 + 48);
    const auto *skf1_49 = buffer.data(skf1 + 49);
    const auto *skf1_50 = buffer.data(skf1 + 50);
    const auto *skf1_53 = buffer.data(skf1 + 53);
    const auto *skf1_55 = buffer.data(skf1 + 55);
    const auto *skf1_56 = buffer.data(skf1 + 56);
    const auto *skf1_58 = buffer.data(skf1 + 58);
    const auto *skf1_59 = buffer.data(skf1 + 59);

    const auto *skg_0 = buffer.data(skg + 0);
    const auto *skg_2 = buffer.data(skg + 2);
    const auto *skg_3 = buffer.data(skg + 3);
    const auto *skg_5 = buffer.data(skg + 5);
    const auto *skg_6 = buffer.data(skg + 6);
    const auto *skg_9 = buffer.data(skg + 9);
    const auto *skg_10 = buffer.data(skg + 10);
    const auto *skg_11 = buffer.data(skg + 11);
    const auto *skg_12 = buffer.data(skg + 12);
    const auto *skg_13 = buffer.data(skg + 13);
    const auto *skg_14 = buffer.data(skg + 14);
    const auto *skg_15 = buffer.data(skg + 15);
    const auto *skg_17 = buffer.data(skg + 17);
    const auto *skg_18 = buffer.data(skg + 18);
    const auto *skg_20 = buffer.data(skg + 20);
    const auto *skg_25 = buffer.data(skg + 25);
    const auto *skg_26 = buffer.data(skg + 26);
    const auto *skg_27 = buffer.data(skg + 27);
    const auto *skg_28 = buffer.data(skg + 28);
    const auto *skg_29 = buffer.data(skg + 29);
    const auto *skg_30 = buffer.data(skg + 30);
    const auto *skg_32 = buffer.data(skg + 32);
    const auto *skg_33 = buffer.data(skg + 33);
    const auto *skg_35 = buffer.data(skg + 35);
    const auto *skg_40 = buffer.data(skg + 40);
    const auto *skg_41 = buffer.data(skg + 41);
    const auto *skg_42 = buffer.data(skg + 42);
    const auto *skg_43 = buffer.data(skg + 43);
    const auto *skg_44 = buffer.data(skg + 44);
    const auto *skg_45 = buffer.data(skg + 45);
    const auto *skg_47 = buffer.data(skg + 47);
    const auto *skg_48 = buffer.data(skg + 48);
    const auto *skg_50 = buffer.data(skg + 50);
    const auto *skg_51 = buffer.data(skg + 51);
    const auto *skg_54 = buffer.data(skg + 54);
    const auto *skg_55 = buffer.data(skg + 55);
    const auto *skg_56 = buffer.data(skg + 56);
    const auto *skg_57 = buffer.data(skg + 57);
    const auto *skg_58 = buffer.data(skg + 58);
    const auto *skg_59 = buffer.data(skg + 59);
    const auto *skg_60 = buffer.data(skg + 60);
    const auto *skg_62 = buffer.data(skg + 62);
    const auto *skg_63 = buffer.data(skg + 63);
    const auto *skg_65 = buffer.data(skg + 65);
    const auto *skg_70 = buffer.data(skg + 70);
    const auto *skg_71 = buffer.data(skg + 71);
    const auto *skg_72 = buffer.data(skg + 72);
    const auto *skg_73 = buffer.data(skg + 73);
    const auto *skg_74 = buffer.data(skg + 74);
    const auto *skg_75 = buffer.data(skg + 75);
    const auto *skg_77 = buffer.data(skg + 77);
    const auto *skg_78 = buffer.data(skg + 78);
    const auto *skg_80 = buffer.data(skg + 80);
    const auto *skg_81 = buffer.data(skg + 81);
    const auto *skg_84 = buffer.data(skg + 84);
    const auto *skg_85 = buffer.data(skg + 85);
    const auto *skg_86 = buffer.data(skg + 86);
    const auto *skg_87 = buffer.data(skg + 87);
    const auto *skg_88 = buffer.data(skg + 88);
    const auto *skg_89 = buffer.data(skg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sig_0, sig_3, skf0_0, skf0_3, \
                         skf1_0, skf1_3, skg_0, skg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sig_0[k]
                 + f_1 * skf0_0[k]
                 - f_2 * skf1_0[k]
                 + f_3 * pc_x[k] * skg_0[k];

        t_1[k] = f_3 * pc_y[k] * skg_0[k];

        t_2[k] = f_3 * pc_z[k] * skg_0[k];

        t_3[k] = f_0 * sig_3[k]
                 + f_4 * skf0_3[k]
                 - f_5 * skf1_3[k]
                 + f_3 * pc_x[k] * skg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sig_5, sig_6, skf0_5, skf0_6, skf1_5, \
                         skf1_6, skg_2, skg_5, skg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * skg_2[k];

        t_5[k] = f_0 * sig_5[k]
                 + f_4 * skf0_5[k]
                 - f_5 * skf1_5[k]
                 + f_3 * pc_x[k] * skg_5[k];

        t_6[k] = f_0 * sig_6[k]
                 + f_6 * skf0_6[k]
                 - f_7 * skf1_6[k]
                 + f_3 * pc_x[k] * skg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, sig_9, sig_10, skf0_9, skf1_9, \
                         skg_3, skg_5, skg_9, skg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * skg_3[k];

        t_8[k] = f_3 * pc_y[k] * skg_5[k];

        t_9[k] = f_0 * sig_9[k]
                 + f_6 * skf0_9[k]
                 - f_7 * skf1_9[k]
                 + f_3 * pc_x[k] * skg_9[k];

        t_10[k] = f_0 * sig_10[k]
                  + f_3 * pc_x[k] * skg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, sig_11, sig_12, sig_13, sig_14, skg_11, \
                         skg_12, skg_13, skg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * sig_11[k]
                  + f_3 * pc_x[k] * skg_11[k];

        t_12[k] = f_0 * sig_12[k]
                  + f_3 * pc_x[k] * skg_12[k];

        t_13[k] = f_0 * sig_13[k]
                  + f_3 * pc_x[k] * skg_13[k];

        t_14[k] = f_0 * sig_14[k]
                  + f_3 * pc_x[k] * skg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, skf0_6, skf0_8, skf0_9, skf1_6, \
                         skf1_8, skf1_9, skg_10, skg_12, skg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * skf0_6[k]
                  - f_2 * skf1_6[k]
                  + f_3 * pc_y[k] * skg_10[k];

        t_16[k] = f_3 * pc_z[k] * skg_10[k];

        t_17[k] = f_4 * skf0_8[k]
                  - f_5 * skf1_8[k]
                  + f_3 * pc_y[k] * skg_12[k];

        t_18[k] = f_6 * skf0_9[k]
                  - f_7 * skf1_9[k]
                  + f_3 * pc_y[k] * skg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, sih0_0, sig_0, \
                         sih1_0, skf0_9, skf1_9, skg_14, skg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * skg_14[k];

        t_20[k] = f_1 * skf0_9[k]
                  - f_2 * skf1_9[k]
                  + f_3 * pc_z[k] * skg_14[k];

        t_21[k] = pb_y[k] * sih0_0[k]
                  - f_8 * pc_y[k] * sih1_0[k];

        t_22[k] = f_9 * sig_0[k]
                  + f_3 * pc_y[k] * skg_15[k];

        t_23[k] = f_3 * pc_z[k] * skg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, sih0_3, sih0_5, sih0_6, sig_1, \
                         sig_2, sig_3, sih1_3, sih1_5, sih1_6, skg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * sih0_3[k]
                  + f_10 * sig_1[k]
                  - f_8 * pc_y[k] * sih1_3[k];

        t_25[k] = f_9 * sig_2[k]
                  + f_3 * pc_y[k] * skg_17[k];

        t_26[k] = pb_y[k] * sih0_5[k]
                  - f_8 * pc_y[k] * sih1_5[k];

        t_27[k] = pb_y[k] * sih0_6[k]
                  + f_11 * sig_3[k]
                  - f_8 * pc_y[k] * sih1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, sih0_9, sig_5, \
                         sig_25, sih1_9, skg_18, skg_20, skg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * skg_18[k];

        t_29[k] = f_9 * sig_5[k]
                  + f_3 * pc_y[k] * skg_20[k];

        t_30[k] = pb_y[k] * sih0_9[k]
                  - f_8 * pc_y[k] * sih1_9[k];

        t_31[k] = f_12 * sig_25[k]
                  + f_3 * pc_x[k] * skg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, sig_26, sig_27, sig_28, sig_29, skg_26, \
                         skg_27, skg_28, skg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * sig_26[k]
                  + f_3 * pc_x[k] * skg_26[k];

        t_33[k] = f_12 * sig_27[k]
                  + f_3 * pc_x[k] * skg_27[k];

        t_34[k] = f_12 * sig_28[k]
                  + f_3 * pc_x[k] * skg_28[k];

        t_35[k] = f_12 * sig_29[k]
                  + f_3 * pc_x[k] * skg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, sig_10, sig_12, skf0_16, skf0_18, \
                         skf1_16, skf1_18, skg_25, skg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * sig_10[k]
                  + f_1 * skf0_16[k]
                  - f_2 * skf1_16[k]
                  + f_3 * pc_y[k] * skg_25[k];

        t_37[k] = f_3 * pc_z[k] * skg_25[k];

        t_38[k] = f_9 * sig_12[k]
                  + f_4 * skf0_18[k]
                  - f_5 * skf1_18[k]
                  + f_3 * pc_y[k] * skg_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, sih0_20, sig_13, sig_14, sih1_20, \
                         skf0_19, skf1_19, skg_28, skg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * sig_13[k]
                  + f_6 * skf0_19[k]
                  - f_7 * skf1_19[k]
                  + f_3 * pc_y[k] * skg_28[k];

        t_40[k] = f_9 * sig_14[k]
                  + f_3 * pc_y[k] * skg_29[k];

        t_41[k] = pb_y[k] * sih0_20[k]
                  - f_8 * pc_y[k] * sih1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, sih0_0, sih0_3, \
                         sig_0, sih1_0, sih1_3, skg_30, skg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sih0_0[k]
                  - f_8 * pc_z[k] * sih1_0[k];

        t_43[k] = f_3 * pc_y[k] * skg_30[k];

        t_44[k] = f_9 * sig_0[k]
                  + f_3 * pc_z[k] * skg_30[k];

        t_45[k] = pb_z[k] * sih0_3[k]
                  - f_8 * pc_z[k] * sih1_3[k];

        t_46[k] = f_3 * pc_y[k] * skg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, sih0_5, sih0_6, sig_2, \
                         sig_3, sih1_5, sih1_6, skg_33, skg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * sih0_5[k]
                  + f_10 * sig_2[k]
                  - f_8 * pc_z[k] * sih1_5[k];

        t_48[k] = pb_z[k] * sih0_6[k]
                  - f_8 * pc_z[k] * sih1_6[k];

        t_49[k] = f_9 * sig_3[k]
                  + f_3 * pc_z[k] * skg_33[k];

        t_50[k] = f_3 * pc_y[k] * skg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, sih0_9, sig_5, sig_40, \
                         sig_41, sig_42, sih1_9, skg_40, skg_41, \
                         skg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * sih0_9[k]
                  + f_11 * sig_5[k]
                  - f_8 * pc_z[k] * sih1_9[k];

        t_52[k] = f_12 * sig_40[k]
                  + f_3 * pc_x[k] * skg_40[k];

        t_53[k] = f_12 * sig_41[k]
                  + f_3 * pc_x[k] * skg_41[k];

        t_54[k] = f_12 * sig_42[k]
                  + f_3 * pc_x[k] * skg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, sih0_15, sig_10, sig_43, \
                         sig_44, sih1_15, skg_40, skg_43, skg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_12 * sig_43[k]
                  + f_3 * pc_x[k] * skg_43[k];

        t_56[k] = f_12 * sig_44[k]
                  + f_3 * pc_x[k] * skg_44[k];

        t_57[k] = pb_z[k] * sih0_15[k]
                  - f_8 * pc_z[k] * sih1_15[k];

        t_58[k] = f_9 * sig_10[k]
                  + f_3 * pc_z[k] * skg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, sig_14, skf0_28, skf0_29, \
                         skf1_28, skf1_29, skg_42, skg_43, skg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * skf0_28[k]
                  - f_5 * skf1_28[k]
                  + f_3 * pc_y[k] * skg_42[k];

        t_60[k] = f_6 * skf0_29[k]
                  - f_7 * skf1_29[k]
                  + f_3 * pc_y[k] * skg_43[k];

        t_61[k] = f_3 * pc_y[k] * skg_44[k];

        t_62[k] = f_9 * sig_14[k]
                  + f_1 * skf0_29[k]
                  - f_2 * skf1_29[k]
                  + f_3 * pc_z[k] * skg_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, sig_15, sig_45, sig_48, \
                         skf0_30, skf0_33, skf1_30, skf1_33, skg_45, \
                         skg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_13 * sig_45[k]
                  + f_1 * skf0_30[k]
                  - f_2 * skf1_30[k]
                  + f_3 * pc_x[k] * skg_45[k];

        t_64[k] = f_10 * sig_15[k]
                  + f_3 * pc_y[k] * skg_45[k];

        t_65[k] = f_3 * pc_z[k] * skg_45[k];

        t_66[k] = f_13 * sig_48[k]
                  + f_4 * skf0_33[k]
                  - f_5 * skf1_33[k]
                  + f_3 * pc_x[k] * skg_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_x, pc_y, sig_17, sig_50, sig_51, skf0_35, \
                         skf0_36, skf1_35, skf1_36, skg_47, skg_50, \
                         skg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * sig_17[k]
                  + f_3 * pc_y[k] * skg_47[k];

        t_68[k] = f_13 * sig_50[k]
                  + f_4 * skf0_35[k]
                  - f_5 * skf1_35[k]
                  + f_3 * pc_x[k] * skg_50[k];

        t_69[k] = f_13 * sig_51[k]
                  + f_6 * skf0_36[k]
                  - f_7 * skf1_36[k]
                  + f_3 * pc_x[k] * skg_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, sig_20, sig_54, sig_55, \
                         skf0_39, skf1_39, skg_48, skg_50, skg_54, \
                         skg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * skg_48[k];

        t_71[k] = f_10 * sig_20[k]
                  + f_3 * pc_y[k] * skg_50[k];

        t_72[k] = f_13 * sig_54[k]
                  + f_6 * skf0_39[k]
                  - f_7 * skf1_39[k]
                  + f_3 * pc_x[k] * skg_54[k];

        t_73[k] = f_13 * sig_55[k]
                  + f_3 * pc_x[k] * skg_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, sig_56, sig_57, sig_58, sig_59, skg_56, \
                         skg_57, skg_58, skg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * sig_56[k]
                  + f_3 * pc_x[k] * skg_56[k];

        t_75[k] = f_13 * sig_57[k]
                  + f_3 * pc_x[k] * skg_57[k];

        t_76[k] = f_13 * sig_58[k]
                  + f_3 * pc_x[k] * skg_58[k];

        t_77[k] = f_13 * sig_59[k]
                  + f_3 * pc_x[k] * skg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, sig_25, sig_27, skf0_36, skf0_38, \
                         skf1_36, skf1_38, skg_55, skg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * sig_25[k]
                  + f_1 * skf0_36[k]
                  - f_2 * skf1_36[k]
                  + f_3 * pc_y[k] * skg_55[k];

        t_79[k] = f_3 * pc_z[k] * skg_55[k];

        t_80[k] = f_10 * sig_27[k]
                  + f_4 * skf0_38[k]
                  - f_5 * skf1_38[k]
                  + f_3 * pc_y[k] * skg_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, sih0_42, sig_28, sig_29, \
                         sih1_42, skf0_39, skf1_39, skg_58, skg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * sig_28[k]
                  + f_6 * skf0_39[k]
                  - f_7 * skf1_39[k]
                  + f_3 * pc_y[k] * skg_58[k];

        t_82[k] = f_10 * sig_29[k]
                  + f_3 * pc_y[k] * skg_59[k];

        t_83[k] = f_1 * skf0_39[k]
                  - f_2 * skf1_39[k]
                  + f_3 * pc_z[k] * skg_59[k];

        t_84[k] = pb_y[k] * sih0_42[k]
                  - f_8 * pc_y[k] * sih1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, sih0_24, sig_15, sig_30, \
                         sig_32, sih1_24, skg_60, skg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * sig_30[k]
                  + f_3 * pc_y[k] * skg_60[k];

        t_86[k] = f_9 * sig_15[k]
                  + f_3 * pc_z[k] * skg_60[k];

        t_87[k] = pb_z[k] * sih0_24[k]
                  - f_8 * pc_z[k] * sih1_24[k];

        t_88[k] = f_9 * sig_32[k]
                  + f_3 * pc_y[k] * skg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, sih0_27, sih0_47, \
                         sig_18, sig_35, sih1_27, sih1_47, skg_63, \
                         skg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * sih0_47[k]
                  - f_8 * pc_y[k] * sih1_47[k];

        t_90[k] = pb_z[k] * sih0_27[k]
                  - f_8 * pc_z[k] * sih1_27[k];

        t_91[k] = f_9 * sig_18[k]
                  + f_3 * pc_z[k] * skg_63[k];

        t_92[k] = f_9 * sig_35[k]
                  + f_3 * pc_y[k] * skg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, sih0_51, sig_70, sig_71, \
                         sig_72, sih1_51, skg_70, skg_71, skg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * sih0_51[k]
                  - f_8 * pc_y[k] * sih1_51[k];

        t_94[k] = f_13 * sig_70[k]
                  + f_3 * pc_x[k] * skg_70[k];

        t_95[k] = f_13 * sig_71[k]
                  + f_3 * pc_x[k] * skg_71[k];

        t_96[k] = f_13 * sig_72[k]
                  + f_3 * pc_x[k] * skg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, sih0_36, sig_25, sig_73, \
                         sig_74, sih1_36, skg_70, skg_73, skg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_13 * sig_73[k]
                  + f_3 * pc_x[k] * skg_73[k];

        t_98[k] = f_13 * sig_74[k]
                  + f_3 * pc_x[k] * skg_74[k];

        t_99[k] = pb_z[k] * sih0_36[k]
                  - f_8 * pc_z[k] * sih1_36[k];

        t_100[k] = f_9 * sig_25[k]
                   + f_3 * pc_z[k] * skg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, sig_42, sig_43, sig_44, skf0_48, skf0_49, \
                         skf1_48, skf1_49, skg_72, skg_73, skg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * sig_42[k]
                   + f_4 * skf0_48[k]
                   - f_5 * skf1_48[k]
                   + f_3 * pc_y[k] * skg_72[k];

        t_102[k] = f_9 * sig_43[k]
                   + f_6 * skf0_49[k]
                   - f_7 * skf1_49[k]
                   + f_3 * pc_y[k] * skg_73[k];

        t_103[k] = f_9 * sig_44[k]
                   + f_3 * pc_y[k] * skg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, sih0_62, sig_30, \
                         sig_75, sih1_62, skf0_50, skf1_50, skg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * sih0_62[k]
                   - f_8 * pc_y[k] * sih1_62[k];

        t_105[k] = f_13 * sig_75[k]
                   + f_1 * skf0_50[k]
                   - f_2 * skf1_50[k]
                   + f_3 * pc_x[k] * skg_75[k];

        t_106[k] = f_3 * pc_y[k] * skg_75[k];

        t_107[k] = f_10 * sig_30[k]
                   + f_3 * pc_z[k] * skg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, sig_78, sig_80, skf0_53, skf0_55, \
                         skf1_53, skf1_55, skg_77, skg_78, skg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_13 * sig_78[k]
                   + f_4 * skf0_53[k]
                   - f_5 * skf1_53[k]
                   + f_3 * pc_x[k] * skg_78[k];

        t_109[k] = f_3 * pc_y[k] * skg_77[k];

        t_110[k] = f_13 * sig_80[k]
                   + f_4 * skf0_55[k]
                   - f_5 * skf1_55[k]
                   + f_3 * pc_x[k] * skg_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_y, pc_z, sig_33, sig_81, skf0_56, \
                         skf1_56, skg_78, skg_80, skg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * sig_81[k]
                   + f_6 * skf0_56[k]
                   - f_7 * skf1_56[k]
                   + f_3 * pc_x[k] * skg_81[k];

        t_112[k] = f_10 * sig_33[k]
                   + f_3 * pc_z[k] * skg_78[k];

        t_113[k] = f_3 * pc_y[k] * skg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, sig_84, sig_85, sig_86, sig_87, \
                         skf0_59, skf1_59, skg_84, skg_85, skg_86, \
                         skg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * sig_84[k]
                   + f_6 * skf0_59[k]
                   - f_7 * skf1_59[k]
                   + f_3 * pc_x[k] * skg_84[k];

        t_115[k] = f_13 * sig_85[k]
                   + f_3 * pc_x[k] * skg_85[k];

        t_116[k] = f_13 * sig_86[k]
                   + f_3 * pc_x[k] * skg_86[k];

        t_117[k] = f_13 * sig_87[k]
                   + f_3 * pc_x[k] * skg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, sig_40, sig_88, sig_89, \
                         skf0_56, skf1_56, skg_85, skg_88, skg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_13 * sig_88[k]
                   + f_3 * pc_x[k] * skg_88[k];

        t_119[k] = f_13 * sig_89[k]
                   + f_3 * pc_x[k] * skg_89[k];

        t_120[k] = f_1 * skf0_56[k]
                   - f_2 * skf1_56[k]
                   + f_3 * pc_y[k] * skg_85[k];

        t_121[k] = f_10 * sig_40[k]
                   + f_3 * pc_z[k] * skg_85[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, sig_44, skf0_58, skf0_59, \
                         skf1_58, skf1_59, skg_87, skg_88, skg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * skf0_58[k]
                   - f_5 * skf1_58[k]
                   + f_3 * pc_y[k] * skg_87[k];

        t_123[k] = f_6 * skf0_59[k]
                   - f_7 * skf1_59[k]
                   + f_3 * pc_y[k] * skg_88[k];

        t_124[k] = f_3 * pc_y[k] * skg_89[k];

        t_125[k] = f_10 * sig_44[k]
                   + f_1 * skf0_59[k]
                   - f_2 * skf1_59[k]
                   + f_3 * pc_z[k] * skg_89[k];
    }
}

static auto
compute_prim_skh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sih0,
                                                          const size_t sig, const size_t sih1,
                                                          const size_t skf0, const size_t skf1,
                                                          const size_t skg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 2.0 / q;

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
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sih0_63 = buffer.data(sih0 + 63);
    const auto *sih0_66 = buffer.data(sih0 + 66);
    const auto *sih0_69 = buffer.data(sih0 + 69);
    const auto *sih0_78 = buffer.data(sih0 + 78);
    const auto *sih0_105 = buffer.data(sih0 + 105);
    const auto *sih0_108 = buffer.data(sih0 + 108);
    const auto *sih0_110 = buffer.data(sih0 + 110);
    const auto *sih0_111 = buffer.data(sih0 + 111);
    const auto *sih0_114 = buffer.data(sih0 + 114);
    const auto *sih0_125 = buffer.data(sih0 + 125);
    const auto *sih0_126 = buffer.data(sih0 + 126);
    const auto *sih0_129 = buffer.data(sih0 + 129);
    const auto *sih0_132 = buffer.data(sih0 + 132);

    const auto *sig_45 = buffer.data(sig + 45);
    const auto *sig_47 = buffer.data(sig + 47);
    const auto *sig_48 = buffer.data(sig + 48);
    const auto *sig_50 = buffer.data(sig + 50);
    const auto *sig_55 = buffer.data(sig + 55);
    const auto *sig_57 = buffer.data(sig + 57);
    const auto *sig_58 = buffer.data(sig + 58);
    const auto *sig_59 = buffer.data(sig + 59);
    const auto *sig_60 = buffer.data(sig + 60);
    const auto *sig_62 = buffer.data(sig + 62);
    const auto *sig_63 = buffer.data(sig + 63);
    const auto *sig_65 = buffer.data(sig + 65);
    const auto *sig_70 = buffer.data(sig + 70);
    const auto *sig_72 = buffer.data(sig + 72);
    const auto *sig_73 = buffer.data(sig + 73);
    const auto *sig_74 = buffer.data(sig + 74);
    const auto *sig_75 = buffer.data(sig + 75);
    const auto *sig_76 = buffer.data(sig + 76);
    const auto *sig_77 = buffer.data(sig + 77);
    const auto *sig_78 = buffer.data(sig + 78);
    const auto *sig_80 = buffer.data(sig + 80);
    const auto *sig_85 = buffer.data(sig + 85);
    const auto *sig_87 = buffer.data(sig + 87);
    const auto *sig_88 = buffer.data(sig + 88);
    const auto *sig_89 = buffer.data(sig + 89);
    const auto *sig_90 = buffer.data(sig + 90);
    const auto *sig_92 = buffer.data(sig + 92);
    const auto *sig_93 = buffer.data(sig + 93);
    const auto *sig_95 = buffer.data(sig + 95);
    const auto *sig_96 = buffer.data(sig + 96);
    const auto *sig_99 = buffer.data(sig + 99);
    const auto *sig_100 = buffer.data(sig + 100);
    const auto *sig_101 = buffer.data(sig + 101);
    const auto *sig_102 = buffer.data(sig + 102);
    const auto *sig_103 = buffer.data(sig + 103);
    const auto *sig_104 = buffer.data(sig + 104);
    const auto *sig_105 = buffer.data(sig + 105);
    const auto *sig_107 = buffer.data(sig + 107);
    const auto *sig_110 = buffer.data(sig + 110);
    const auto *sig_114 = buffer.data(sig + 114);
    const auto *sig_115 = buffer.data(sig + 115);
    const auto *sig_116 = buffer.data(sig + 116);
    const auto *sig_117 = buffer.data(sig + 117);
    const auto *sig_118 = buffer.data(sig + 118);
    const auto *sig_119 = buffer.data(sig + 119);
    const auto *sig_130 = buffer.data(sig + 130);
    const auto *sig_131 = buffer.data(sig + 131);
    const auto *sig_132 = buffer.data(sig + 132);
    const auto *sig_133 = buffer.data(sig + 133);
    const auto *sig_134 = buffer.data(sig + 134);
    const auto *sig_135 = buffer.data(sig + 135);
    const auto *sig_138 = buffer.data(sig + 138);
    const auto *sig_140 = buffer.data(sig + 140);
    const auto *sig_141 = buffer.data(sig + 141);
    const auto *sig_144 = buffer.data(sig + 144);
    const auto *sig_145 = buffer.data(sig + 145);
    const auto *sig_146 = buffer.data(sig + 146);
    const auto *sig_147 = buffer.data(sig + 147);
    const auto *sig_148 = buffer.data(sig + 148);
    const auto *sig_149 = buffer.data(sig + 149);
    const auto *sig_150 = buffer.data(sig + 150);
    const auto *sig_153 = buffer.data(sig + 153);
    const auto *sig_155 = buffer.data(sig + 155);
    const auto *sig_156 = buffer.data(sig + 156);
    const auto *sig_159 = buffer.data(sig + 159);
    const auto *sig_160 = buffer.data(sig + 160);
    const auto *sig_161 = buffer.data(sig + 161);
    const auto *sig_162 = buffer.data(sig + 162);
    const auto *sig_163 = buffer.data(sig + 163);
    const auto *sig_164 = buffer.data(sig + 164);
    const auto *sig_170 = buffer.data(sig + 170);
    const auto *sig_174 = buffer.data(sig + 174);
    const auto *sig_175 = buffer.data(sig + 175);
    const auto *sig_176 = buffer.data(sig + 176);

    const auto *sih1_63 = buffer.data(sih1 + 63);
    const auto *sih1_66 = buffer.data(sih1 + 66);
    const auto *sih1_69 = buffer.data(sih1 + 69);
    const auto *sih1_78 = buffer.data(sih1 + 78);
    const auto *sih1_105 = buffer.data(sih1 + 105);
    const auto *sih1_108 = buffer.data(sih1 + 108);
    const auto *sih1_110 = buffer.data(sih1 + 110);
    const auto *sih1_111 = buffer.data(sih1 + 111);
    const auto *sih1_114 = buffer.data(sih1 + 114);
    const auto *sih1_125 = buffer.data(sih1 + 125);
    const auto *sih1_126 = buffer.data(sih1 + 126);
    const auto *sih1_129 = buffer.data(sih1 + 129);
    const auto *sih1_132 = buffer.data(sih1 + 132);

    const auto *skf0_60 = buffer.data(skf0 + 60);
    const auto *skf0_63 = buffer.data(skf0 + 63);
    const auto *skf0_65 = buffer.data(skf0 + 65);
    const auto *skf0_66 = buffer.data(skf0 + 66);
    const auto *skf0_68 = buffer.data(skf0 + 68);
    const auto *skf0_69 = buffer.data(skf0 + 69);
    const auto *skf0_75 = buffer.data(skf0 + 75);
    const auto *skf0_78 = buffer.data(skf0 + 78);
    const auto *skf0_79 = buffer.data(skf0 + 79);
    const auto *skf0_86 = buffer.data(skf0 + 86);
    const auto *skf0_88 = buffer.data(skf0 + 88);
    const auto *skf0_89 = buffer.data(skf0 + 89);
    const auto *skf0_90 = buffer.data(skf0 + 90);
    const auto *skf0_93 = buffer.data(skf0 + 93);
    const auto *skf0_95 = buffer.data(skf0 + 95);
    const auto *skf0_96 = buffer.data(skf0 + 96);
    const auto *skf0_98 = buffer.data(skf0 + 98);
    const auto *skf0_99 = buffer.data(skf0 + 99);
    const auto *skf0_100 = buffer.data(skf0 + 100);
    const auto *skf0_103 = buffer.data(skf0 + 103);
    const auto *skf0_105 = buffer.data(skf0 + 105);
    const auto *skf0_106 = buffer.data(skf0 + 106);
    const auto *skf0_108 = buffer.data(skf0 + 108);
    const auto *skf0_109 = buffer.data(skf0 + 109);
    const auto *skf0_115 = buffer.data(skf0 + 115);
    const auto *skf0_119 = buffer.data(skf0 + 119);

    const auto *skf1_60 = buffer.data(skf1 + 60);
    const auto *skf1_63 = buffer.data(skf1 + 63);
    const auto *skf1_65 = buffer.data(skf1 + 65);
    const auto *skf1_66 = buffer.data(skf1 + 66);
    const auto *skf1_68 = buffer.data(skf1 + 68);
    const auto *skf1_69 = buffer.data(skf1 + 69);
    const auto *skf1_75 = buffer.data(skf1 + 75);
    const auto *skf1_78 = buffer.data(skf1 + 78);
    const auto *skf1_79 = buffer.data(skf1 + 79);
    const auto *skf1_86 = buffer.data(skf1 + 86);
    const auto *skf1_88 = buffer.data(skf1 + 88);
    const auto *skf1_89 = buffer.data(skf1 + 89);
    const auto *skf1_90 = buffer.data(skf1 + 90);
    const auto *skf1_93 = buffer.data(skf1 + 93);
    const auto *skf1_95 = buffer.data(skf1 + 95);
    const auto *skf1_96 = buffer.data(skf1 + 96);
    const auto *skf1_98 = buffer.data(skf1 + 98);
    const auto *skf1_99 = buffer.data(skf1 + 99);
    const auto *skf1_100 = buffer.data(skf1 + 100);
    const auto *skf1_103 = buffer.data(skf1 + 103);
    const auto *skf1_105 = buffer.data(skf1 + 105);
    const auto *skf1_106 = buffer.data(skf1 + 106);
    const auto *skf1_108 = buffer.data(skf1 + 108);
    const auto *skf1_109 = buffer.data(skf1 + 109);
    const auto *skf1_115 = buffer.data(skf1 + 115);
    const auto *skf1_119 = buffer.data(skf1 + 119);

    const auto *skg_90 = buffer.data(skg + 90);
    const auto *skg_92 = buffer.data(skg + 92);
    const auto *skg_93 = buffer.data(skg + 93);
    const auto *skg_95 = buffer.data(skg + 95);
    const auto *skg_96 = buffer.data(skg + 96);
    const auto *skg_99 = buffer.data(skg + 99);
    const auto *skg_100 = buffer.data(skg + 100);
    const auto *skg_101 = buffer.data(skg + 101);
    const auto *skg_102 = buffer.data(skg + 102);
    const auto *skg_103 = buffer.data(skg + 103);
    const auto *skg_104 = buffer.data(skg + 104);
    const auto *skg_105 = buffer.data(skg + 105);
    const auto *skg_107 = buffer.data(skg + 107);
    const auto *skg_108 = buffer.data(skg + 108);
    const auto *skg_110 = buffer.data(skg + 110);
    const auto *skg_114 = buffer.data(skg + 114);
    const auto *skg_115 = buffer.data(skg + 115);
    const auto *skg_116 = buffer.data(skg + 116);
    const auto *skg_117 = buffer.data(skg + 117);
    const auto *skg_118 = buffer.data(skg + 118);
    const auto *skg_119 = buffer.data(skg + 119);
    const auto *skg_120 = buffer.data(skg + 120);
    const auto *skg_122 = buffer.data(skg + 122);
    const auto *skg_123 = buffer.data(skg + 123);
    const auto *skg_125 = buffer.data(skg + 125);
    const auto *skg_130 = buffer.data(skg + 130);
    const auto *skg_131 = buffer.data(skg + 131);
    const auto *skg_132 = buffer.data(skg + 132);
    const auto *skg_133 = buffer.data(skg + 133);
    const auto *skg_134 = buffer.data(skg + 134);
    const auto *skg_135 = buffer.data(skg + 135);
    const auto *skg_137 = buffer.data(skg + 137);
    const auto *skg_138 = buffer.data(skg + 138);
    const auto *skg_140 = buffer.data(skg + 140);
    const auto *skg_141 = buffer.data(skg + 141);
    const auto *skg_144 = buffer.data(skg + 144);
    const auto *skg_145 = buffer.data(skg + 145);
    const auto *skg_146 = buffer.data(skg + 146);
    const auto *skg_147 = buffer.data(skg + 147);
    const auto *skg_148 = buffer.data(skg + 148);
    const auto *skg_149 = buffer.data(skg + 149);
    const auto *skg_150 = buffer.data(skg + 150);
    const auto *skg_152 = buffer.data(skg + 152);
    const auto *skg_153 = buffer.data(skg + 153);
    const auto *skg_155 = buffer.data(skg + 155);
    const auto *skg_156 = buffer.data(skg + 156);
    const auto *skg_159 = buffer.data(skg + 159);
    const auto *skg_160 = buffer.data(skg + 160);
    const auto *skg_161 = buffer.data(skg + 161);
    const auto *skg_162 = buffer.data(skg + 162);
    const auto *skg_163 = buffer.data(skg + 163);
    const auto *skg_164 = buffer.data(skg + 164);
    const auto *skg_165 = buffer.data(skg + 165);
    const auto *skg_167 = buffer.data(skg + 167);
    const auto *skg_168 = buffer.data(skg + 168);
    const auto *skg_170 = buffer.data(skg + 170);
    const auto *skg_174 = buffer.data(skg + 174);
    const auto *skg_175 = buffer.data(skg + 175);
    const auto *skg_176 = buffer.data(skg + 176);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, sig_45, sig_90, sig_93, \
                         skf0_60, skf0_63, skf1_60, skf1_63, skg_90, \
                         skg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_14 * sig_90[k]
                   + f_1 * skf0_60[k]
                   - f_2 * skf1_60[k]
                   + f_3 * pc_x[k] * skg_90[k];

        t_127[k] = f_11 * sig_45[k]
                   + f_3 * pc_y[k] * skg_90[k];

        t_128[k] = f_3 * pc_z[k] * skg_90[k];

        t_129[k] = f_14 * sig_93[k]
                   + f_4 * skf0_63[k]
                   - f_5 * skf1_63[k]
                   + f_3 * pc_x[k] * skg_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pc_x, pc_y, sig_47, sig_95, sig_96, skf0_65, \
                         skf0_66, skf1_65, skf1_66, skg_92, skg_95, \
                         skg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * sig_47[k]
                   + f_3 * pc_y[k] * skg_92[k];

        t_131[k] = f_14 * sig_95[k]
                   + f_4 * skf0_65[k]
                   - f_5 * skf1_65[k]
                   + f_3 * pc_x[k] * skg_95[k];

        t_132[k] = f_14 * sig_96[k]
                   + f_6 * skf0_66[k]
                   - f_7 * skf1_66[k]
                   + f_3 * pc_x[k] * skg_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, pc_z, sig_50, sig_99, \
                         sig_100, skf0_69, skf1_69, skg_93, skg_95, skg_99, \
                         skg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * skg_93[k];

        t_134[k] = f_11 * sig_50[k]
                   + f_3 * pc_y[k] * skg_95[k];

        t_135[k] = f_14 * sig_99[k]
                   + f_6 * skf0_69[k]
                   - f_7 * skf1_69[k]
                   + f_3 * pc_x[k] * skg_99[k];

        t_136[k] = f_14 * sig_100[k]
                   + f_3 * pc_x[k] * skg_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, sig_101, sig_102, sig_103, sig_104, \
                         skg_101, skg_102, skg_103, skg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * sig_101[k]
                   + f_3 * pc_x[k] * skg_101[k];

        t_138[k] = f_14 * sig_102[k]
                   + f_3 * pc_x[k] * skg_102[k];

        t_139[k] = f_14 * sig_103[k]
                   + f_3 * pc_x[k] * skg_103[k];

        t_140[k] = f_14 * sig_104[k]
                   + f_3 * pc_x[k] * skg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, sig_55, sig_57, skf0_66, skf0_68, \
                         skf1_66, skf1_68, skg_100, skg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * sig_55[k]
                   + f_1 * skf0_66[k]
                   - f_2 * skf1_66[k]
                   + f_3 * pc_y[k] * skg_100[k];

        t_142[k] = f_3 * pc_z[k] * skg_100[k];

        t_143[k] = f_11 * sig_57[k]
                   + f_4 * skf0_68[k]
                   - f_5 * skf1_68[k]
                   + f_3 * pc_y[k] * skg_102[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_z, pc_y, pc_z, sih0_63, sig_58, \
                         sig_59, sih1_63, skf0_69, skf1_69, skg_103, \
                         skg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_11 * sig_58[k]
                   + f_6 * skf0_69[k]
                   - f_7 * skf1_69[k]
                   + f_3 * pc_y[k] * skg_103[k];

        t_145[k] = f_11 * sig_59[k]
                   + f_3 * pc_y[k] * skg_104[k];

        t_146[k] = f_1 * skf0_69[k]
                   - f_2 * skf1_69[k]
                   + f_3 * pc_z[k] * skg_104[k];

        t_147[k] = pb_z[k] * sih0_63[k]
                   - f_8 * pc_z[k] * sih1_63[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_z, pc_y, pc_z, sih0_66, sig_45, \
                         sig_60, sig_62, sih1_66, skg_105, skg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * sig_60[k]
                   + f_3 * pc_y[k] * skg_105[k];

        t_149[k] = f_9 * sig_45[k]
                   + f_3 * pc_z[k] * skg_105[k];

        t_150[k] = pb_z[k] * sih0_66[k]
                   - f_8 * pc_z[k] * sih1_66[k];

        t_151[k] = f_10 * sig_62[k]
                   + f_3 * pc_y[k] * skg_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_z, pc_x, pc_z, sih0_69, sig_48, sig_110, \
                         sih1_69, skf0_75, skf1_75, skg_108, skg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_14 * sig_110[k]
                   + f_4 * skf0_75[k]
                   - f_5 * skf1_75[k]
                   + f_3 * pc_x[k] * skg_110[k];

        t_153[k] = pb_z[k] * sih0_69[k]
                   - f_8 * pc_z[k] * sih1_69[k];

        t_154[k] = f_9 * sig_48[k]
                   + f_3 * pc_z[k] * skg_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_y, sig_65, sig_114, sig_115, \
                         sig_116, skf0_79, skf1_79, skg_110, skg_114, skg_115, \
                         skg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * sig_65[k]
                   + f_3 * pc_y[k] * skg_110[k];

        t_156[k] = f_14 * sig_114[k]
                   + f_6 * skf0_79[k]
                   - f_7 * skf1_79[k]
                   + f_3 * pc_x[k] * skg_114[k];

        t_157[k] = f_14 * sig_115[k]
                   + f_3 * pc_x[k] * skg_115[k];

        t_158[k] = f_14 * sig_116[k]
                   + f_3 * pc_x[k] * skg_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_z, pc_x, pc_z, sih0_78, sig_117, \
                         sig_118, sig_119, sih1_78, skg_117, skg_118, \
                         skg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_14 * sig_117[k]
                   + f_3 * pc_x[k] * skg_117[k];

        t_160[k] = f_14 * sig_118[k]
                   + f_3 * pc_x[k] * skg_118[k];

        t_161[k] = f_14 * sig_119[k]
                   + f_3 * pc_x[k] * skg_119[k];

        t_162[k] = pb_z[k] * sih0_78[k]
                   - f_8 * pc_z[k] * sih1_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pc_y, pc_z, sig_55, sig_72, sig_73, skf0_78, \
                         skf0_79, skf1_78, skf1_79, skg_115, skg_117, \
                         skg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * sig_55[k]
                   + f_3 * pc_z[k] * skg_115[k];

        t_164[k] = f_10 * sig_72[k]
                   + f_4 * skf0_78[k]
                   - f_5 * skf1_78[k]
                   + f_3 * pc_y[k] * skg_117[k];

        t_165[k] = f_10 * sig_73[k]
                   + f_6 * skf0_79[k]
                   - f_7 * skf1_79[k]
                   + f_3 * pc_y[k] * skg_118[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_y, pc_y, pc_z, sih0_105, sig_59, \
                         sig_74, sig_75, sih1_105, skf0_79, skf1_79, skg_119, \
                         skg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * sig_74[k]
                   + f_3 * pc_y[k] * skg_119[k];

        t_167[k] = f_9 * sig_59[k]
                   + f_1 * skf0_79[k]
                   - f_2 * skf1_79[k]
                   + f_3 * pc_z[k] * skg_119[k];

        t_168[k] = pb_y[k] * sih0_105[k]
                   - f_8 * pc_y[k] * sih1_105[k];

        t_169[k] = f_9 * sig_75[k]
                   + f_3 * pc_y[k] * skg_120[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_y, pc_y, pc_z, sih0_108, sih0_110, \
                         sig_60, sig_76, sig_77, sih1_108, sih1_110, skg_120, \
                         skg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * sig_60[k]
                   + f_3 * pc_z[k] * skg_120[k];

        t_171[k] = pb_y[k] * sih0_108[k]
                   + f_10 * sig_76[k]
                   - f_8 * pc_y[k] * sih1_108[k];

        t_172[k] = f_9 * sig_77[k]
                   + f_3 * pc_y[k] * skg_122[k];

        t_173[k] = pb_y[k] * sih0_110[k]
                   - f_8 * pc_y[k] * sih1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_y, pc_y, pc_z, sih0_111, sih0_114, \
                         sig_63, sig_78, sig_80, sih1_111, sih1_114, skg_123, \
                         skg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_y[k] * sih0_111[k]
                   + f_11 * sig_78[k]
                   - f_8 * pc_y[k] * sih1_111[k];

        t_175[k] = f_10 * sig_63[k]
                   + f_3 * pc_z[k] * skg_123[k];

        t_176[k] = f_9 * sig_80[k]
                   + f_3 * pc_y[k] * skg_125[k];

        t_177[k] = pb_y[k] * sih0_114[k]
                   - f_8 * pc_y[k] * sih1_114[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pc_x, sig_130, sig_131, sig_132, \
                         sig_133, sig_134, skg_130, skg_131, skg_132, skg_133, \
                         skg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * sig_130[k]
                   + f_3 * pc_x[k] * skg_130[k];

        t_179[k] = f_14 * sig_131[k]
                   + f_3 * pc_x[k] * skg_131[k];

        t_180[k] = f_14 * sig_132[k]
                   + f_3 * pc_x[k] * skg_132[k];

        t_181[k] = f_14 * sig_133[k]
                   + f_3 * pc_x[k] * skg_133[k];

        t_182[k] = f_14 * sig_134[k]
                   + f_3 * pc_x[k] * skg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, sig_70, sig_85, sig_87, skf0_86, \
                         skf0_88, skf1_86, skf1_88, skg_130, skg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * sig_85[k]
                   + f_1 * skf0_86[k]
                   - f_2 * skf1_86[k]
                   + f_3 * pc_y[k] * skg_130[k];

        t_184[k] = f_10 * sig_70[k]
                   + f_3 * pc_z[k] * skg_130[k];

        t_185[k] = f_9 * sig_87[k]
                   + f_4 * skf0_88[k]
                   - f_5 * skf1_88[k]
                   + f_3 * pc_y[k] * skg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_y, pc_y, sih0_125, sig_88, sig_89, sih1_125, \
                         skf0_89, skf1_89, skg_133, skg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * sig_88[k]
                   + f_6 * skf0_89[k]
                   - f_7 * skf1_89[k]
                   + f_3 * pc_y[k] * skg_133[k];

        t_187[k] = f_9 * sig_89[k]
                   + f_3 * pc_y[k] * skg_134[k];

        t_188[k] = pb_y[k] * sih0_125[k]
                   - f_8 * pc_y[k] * sih1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_x, pc_y, pc_z, sig_75, sig_135, \
                         sig_138, skf0_90, skf0_93, skf1_90, skf1_93, skg_135, \
                         skg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * sig_135[k]
                   + f_1 * skf0_90[k]
                   - f_2 * skf1_90[k]
                   + f_3 * pc_x[k] * skg_135[k];

        t_190[k] = f_3 * pc_y[k] * skg_135[k];

        t_191[k] = f_11 * sig_75[k]
                   + f_3 * pc_z[k] * skg_135[k];

        t_192[k] = f_14 * sig_138[k]
                   + f_4 * skf0_93[k]
                   - f_5 * skf1_93[k]
                   + f_3 * pc_x[k] * skg_138[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, sig_140, sig_141, skf0_95, skf0_96, \
                         skf1_95, skf1_96, skg_137, skg_140, skg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * skg_137[k];

        t_194[k] = f_14 * sig_140[k]
                   + f_4 * skf0_95[k]
                   - f_5 * skf1_95[k]
                   + f_3 * pc_x[k] * skg_140[k];

        t_195[k] = f_14 * sig_141[k]
                   + f_6 * skf0_96[k]
                   - f_7 * skf1_96[k]
                   + f_3 * pc_x[k] * skg_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, pc_z, sig_78, sig_144, \
                         sig_145, skf0_99, skf1_99, skg_138, skg_140, skg_144, \
                         skg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_11 * sig_78[k]
                   + f_3 * pc_z[k] * skg_138[k];

        t_197[k] = f_3 * pc_y[k] * skg_140[k];

        t_198[k] = f_14 * sig_144[k]
                   + f_6 * skf0_99[k]
                   - f_7 * skf1_99[k]
                   + f_3 * pc_x[k] * skg_144[k];

        t_199[k] = f_14 * sig_145[k]
                   + f_3 * pc_x[k] * skg_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, sig_146, sig_147, sig_148, sig_149, \
                         skg_146, skg_147, skg_148, skg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_14 * sig_146[k]
                   + f_3 * pc_x[k] * skg_146[k];

        t_201[k] = f_14 * sig_147[k]
                   + f_3 * pc_x[k] * skg_147[k];

        t_202[k] = f_14 * sig_148[k]
                   + f_3 * pc_x[k] * skg_148[k];

        t_203[k] = f_14 * sig_149[k]
                   + f_3 * pc_x[k] * skg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, sig_85, skf0_96, skf0_98, \
                         skf0_99, skf1_96, skf1_98, skf1_99, skg_145, skg_147, \
                         skg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * skf0_96[k]
                   - f_2 * skf1_96[k]
                   + f_3 * pc_y[k] * skg_145[k];

        t_205[k] = f_11 * sig_85[k]
                   + f_3 * pc_z[k] * skg_145[k];

        t_206[k] = f_4 * skf0_98[k]
                   - f_5 * skf1_98[k]
                   + f_3 * pc_y[k] * skg_147[k];

        t_207[k] = f_6 * skf0_99[k]
                   - f_7 * skf1_99[k]
                   + f_3 * pc_y[k] * skg_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_x, pc_y, pc_z, sig_89, sig_90, \
                         sig_150, skf0_99, skf0_100, skf1_99, skf1_100, skg_149, \
                         skg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * skg_149[k];

        t_209[k] = f_11 * sig_89[k]
                   + f_1 * skf0_99[k]
                   - f_2 * skf1_99[k]
                   + f_3 * pc_z[k] * skg_149[k];

        t_210[k] = f_11 * sig_150[k]
                   + f_1 * skf0_100[k]
                   - f_2 * skf1_100[k]
                   + f_3 * pc_x[k] * skg_150[k];

        t_211[k] = f_14 * sig_90[k]
                   + f_3 * pc_y[k] * skg_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, sig_92, sig_153, skf0_103, \
                         skf1_103, skg_150, skg_152, skg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * skg_150[k];

        t_213[k] = f_11 * sig_153[k]
                   + f_4 * skf0_103[k]
                   - f_5 * skf1_103[k]
                   + f_3 * pc_x[k] * skg_153[k];

        t_214[k] = f_14 * sig_92[k]
                   + f_3 * pc_y[k] * skg_152[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pc_x, pc_z, sig_155, sig_156, skf0_105, \
                         skf0_106, skf1_105, skf1_106, skg_153, skg_155, \
                         skg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_11 * sig_155[k]
                   + f_4 * skf0_105[k]
                   - f_5 * skf1_105[k]
                   + f_3 * pc_x[k] * skg_155[k];

        t_216[k] = f_11 * sig_156[k]
                   + f_6 * skf0_106[k]
                   - f_7 * skf1_106[k]
                   + f_3 * pc_x[k] * skg_156[k];

        t_217[k] = f_3 * pc_z[k] * skg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, sig_95, sig_159, sig_160, \
                         sig_161, skf0_109, skf1_109, skg_155, skg_159, skg_160, \
                         skg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_14 * sig_95[k]
                   + f_3 * pc_y[k] * skg_155[k];

        t_219[k] = f_11 * sig_159[k]
                   + f_6 * skf0_109[k]
                   - f_7 * skf1_109[k]
                   + f_3 * pc_x[k] * skg_159[k];

        t_220[k] = f_11 * sig_160[k]
                   + f_3 * pc_x[k] * skg_160[k];

        t_221[k] = f_11 * sig_161[k]
                   + f_3 * pc_x[k] * skg_161[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, sig_100, sig_162, sig_163, \
                         sig_164, skf0_106, skf1_106, skg_160, skg_162, skg_163, \
                         skg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_11 * sig_162[k]
                   + f_3 * pc_x[k] * skg_162[k];

        t_223[k] = f_11 * sig_163[k]
                   + f_3 * pc_x[k] * skg_163[k];

        t_224[k] = f_11 * sig_164[k]
                   + f_3 * pc_x[k] * skg_164[k];

        t_225[k] = f_14 * sig_100[k]
                   + f_1 * skf0_106[k]
                   - f_2 * skf1_106[k]
                   + f_3 * pc_y[k] * skg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pc_z, sig_102, sig_103, skf0_108, \
                         skf0_109, skf1_108, skf1_109, skg_160, skg_162, \
                         skg_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * skg_160[k];

        t_227[k] = f_14 * sig_102[k]
                   + f_4 * skf0_108[k]
                   - f_5 * skf1_108[k]
                   + f_3 * pc_y[k] * skg_162[k];

        t_228[k] = f_14 * sig_103[k]
                   + f_6 * skf0_109[k]
                   - f_7 * skf1_109[k]
                   + f_3 * pc_y[k] * skg_163[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_z, pc_y, pc_z, sih0_126, sig_104, \
                         sig_105, sih1_126, skf0_109, skf1_109, skg_164, \
                         skg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_14 * sig_104[k]
                   + f_3 * pc_y[k] * skg_164[k];

        t_230[k] = f_1 * skf0_109[k]
                   - f_2 * skf1_109[k]
                   + f_3 * pc_z[k] * skg_164[k];

        t_231[k] = pb_z[k] * sih0_126[k]
                   - f_8 * pc_z[k] * sih1_126[k];

        t_232[k] = f_11 * sig_105[k]
                   + f_3 * pc_y[k] * skg_165[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_z, pc_y, pc_z, sih0_129, sig_90, sig_107, \
                         sih1_129, skg_165, skg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_9 * sig_90[k]
                   + f_3 * pc_z[k] * skg_165[k];

        t_234[k] = pb_z[k] * sih0_129[k]
                   - f_8 * pc_z[k] * sih1_129[k];

        t_235[k] = f_11 * sig_107[k]
                   + f_3 * pc_y[k] * skg_167[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_z, pc_x, pc_z, sih0_132, sig_93, sig_170, \
                         sih1_132, skf0_115, skf1_115, skg_168, \
                         skg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_11 * sig_170[k]
                   + f_4 * skf0_115[k]
                   - f_5 * skf1_115[k]
                   + f_3 * pc_x[k] * skg_170[k];

        t_237[k] = pb_z[k] * sih0_132[k]
                   - f_8 * pc_z[k] * sih1_132[k];

        t_238[k] = f_9 * sig_93[k]
                   + f_3 * pc_z[k] * skg_168[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pc_x, pc_y, sig_110, sig_174, sig_175, \
                         sig_176, skf0_119, skf1_119, skg_170, skg_174, skg_175, \
                         skg_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_11 * sig_110[k]
                   + f_3 * pc_y[k] * skg_170[k];

        t_240[k] = f_11 * sig_174[k]
                   + f_6 * skf0_119[k]
                   - f_7 * skf1_119[k]
                   + f_3 * pc_x[k] * skg_174[k];

        t_241[k] = f_11 * sig_175[k]
                   + f_3 * pc_x[k] * skg_175[k];

        t_242[k] = f_11 * sig_176[k]
                   + f_3 * pc_x[k] * skg_176[k];
    }
}

static auto
compute_prim_skh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sih0,
                                                          const size_t sig, const size_t sih1,
                                                          const size_t skf0, const size_t skf1,
                                                          const size_t skg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sih0_141 = buffer.data(sih0 + 141);
    const auto *sih0_189 = buffer.data(sih0 + 189);
    const auto *sih0_192 = buffer.data(sih0 + 192);
    const auto *sih0_194 = buffer.data(sih0 + 194);
    const auto *sih0_195 = buffer.data(sih0 + 195);
    const auto *sih0_198 = buffer.data(sih0 + 198);
    const auto *sih0_209 = buffer.data(sih0 + 209);
    const auto *sih0_210 = buffer.data(sih0 + 210);
    const auto *sih0_213 = buffer.data(sih0 + 213);
    const auto *sih0_216 = buffer.data(sih0 + 216);
    const auto *sih0_225 = buffer.data(sih0 + 225);

    const auto *sig_100 = buffer.data(sig + 100);
    const auto *sig_104 = buffer.data(sig + 104);
    const auto *sig_105 = buffer.data(sig + 105);
    const auto *sig_108 = buffer.data(sig + 108);
    const auto *sig_115 = buffer.data(sig + 115);
    const auto *sig_117 = buffer.data(sig + 117);
    const auto *sig_118 = buffer.data(sig + 118);
    const auto *sig_119 = buffer.data(sig + 119);
    const auto *sig_120 = buffer.data(sig + 120);
    const auto *sig_122 = buffer.data(sig + 122);
    const auto *sig_123 = buffer.data(sig + 123);
    const auto *sig_125 = buffer.data(sig + 125);
    const auto *sig_130 = buffer.data(sig + 130);
    const auto *sig_132 = buffer.data(sig + 132);
    const auto *sig_133 = buffer.data(sig + 133);
    const auto *sig_134 = buffer.data(sig + 134);
    const auto *sig_135 = buffer.data(sig + 135);
    const auto *sig_136 = buffer.data(sig + 136);
    const auto *sig_137 = buffer.data(sig + 137);
    const auto *sig_138 = buffer.data(sig + 138);
    const auto *sig_140 = buffer.data(sig + 140);
    const auto *sig_145 = buffer.data(sig + 145);
    const auto *sig_147 = buffer.data(sig + 147);
    const auto *sig_148 = buffer.data(sig + 148);
    const auto *sig_149 = buffer.data(sig + 149);
    const auto *sig_150 = buffer.data(sig + 150);
    const auto *sig_152 = buffer.data(sig + 152);
    const auto *sig_153 = buffer.data(sig + 153);
    const auto *sig_155 = buffer.data(sig + 155);
    const auto *sig_160 = buffer.data(sig + 160);
    const auto *sig_162 = buffer.data(sig + 162);
    const auto *sig_163 = buffer.data(sig + 163);
    const auto *sig_164 = buffer.data(sig + 164);
    const auto *sig_165 = buffer.data(sig + 165);
    const auto *sig_167 = buffer.data(sig + 167);
    const auto *sig_170 = buffer.data(sig + 170);
    const auto *sig_177 = buffer.data(sig + 177);
    const auto *sig_178 = buffer.data(sig + 178);
    const auto *sig_179 = buffer.data(sig + 179);
    const auto *sig_180 = buffer.data(sig + 180);
    const auto *sig_183 = buffer.data(sig + 183);
    const auto *sig_185 = buffer.data(sig + 185);
    const auto *sig_186 = buffer.data(sig + 186);
    const auto *sig_189 = buffer.data(sig + 189);
    const auto *sig_190 = buffer.data(sig + 190);
    const auto *sig_191 = buffer.data(sig + 191);
    const auto *sig_192 = buffer.data(sig + 192);
    const auto *sig_193 = buffer.data(sig + 193);
    const auto *sig_194 = buffer.data(sig + 194);
    const auto *sig_205 = buffer.data(sig + 205);
    const auto *sig_206 = buffer.data(sig + 206);
    const auto *sig_207 = buffer.data(sig + 207);
    const auto *sig_208 = buffer.data(sig + 208);
    const auto *sig_209 = buffer.data(sig + 209);
    const auto *sig_210 = buffer.data(sig + 210);
    const auto *sig_213 = buffer.data(sig + 213);
    const auto *sig_215 = buffer.data(sig + 215);
    const auto *sig_216 = buffer.data(sig + 216);
    const auto *sig_219 = buffer.data(sig + 219);
    const auto *sig_220 = buffer.data(sig + 220);
    const auto *sig_221 = buffer.data(sig + 221);
    const auto *sig_222 = buffer.data(sig + 222);
    const auto *sig_223 = buffer.data(sig + 223);
    const auto *sig_224 = buffer.data(sig + 224);
    const auto *sig_225 = buffer.data(sig + 225);
    const auto *sig_228 = buffer.data(sig + 228);
    const auto *sig_230 = buffer.data(sig + 230);
    const auto *sig_231 = buffer.data(sig + 231);
    const auto *sig_234 = buffer.data(sig + 234);
    const auto *sig_235 = buffer.data(sig + 235);
    const auto *sig_236 = buffer.data(sig + 236);
    const auto *sig_237 = buffer.data(sig + 237);
    const auto *sig_238 = buffer.data(sig + 238);
    const auto *sig_239 = buffer.data(sig + 239);
    const auto *sig_245 = buffer.data(sig + 245);
    const auto *sig_249 = buffer.data(sig + 249);
    const auto *sig_250 = buffer.data(sig + 250);
    const auto *sig_251 = buffer.data(sig + 251);
    const auto *sig_252 = buffer.data(sig + 252);
    const auto *sig_253 = buffer.data(sig + 253);
    const auto *sig_254 = buffer.data(sig + 254);
    const auto *sig_255 = buffer.data(sig + 255);

    const auto *sih1_141 = buffer.data(sih1 + 141);
    const auto *sih1_189 = buffer.data(sih1 + 189);
    const auto *sih1_192 = buffer.data(sih1 + 192);
    const auto *sih1_194 = buffer.data(sih1 + 194);
    const auto *sih1_195 = buffer.data(sih1 + 195);
    const auto *sih1_198 = buffer.data(sih1 + 198);
    const auto *sih1_209 = buffer.data(sih1 + 209);
    const auto *sih1_210 = buffer.data(sih1 + 210);
    const auto *sih1_213 = buffer.data(sih1 + 213);
    const auto *sih1_216 = buffer.data(sih1 + 216);
    const auto *sih1_225 = buffer.data(sih1 + 225);

    const auto *skf0_118 = buffer.data(skf0 + 118);
    const auto *skf0_119 = buffer.data(skf0 + 119);
    const auto *skf0_120 = buffer.data(skf0 + 120);
    const auto *skf0_123 = buffer.data(skf0 + 123);
    const auto *skf0_125 = buffer.data(skf0 + 125);
    const auto *skf0_126 = buffer.data(skf0 + 126);
    const auto *skf0_128 = buffer.data(skf0 + 128);
    const auto *skf0_129 = buffer.data(skf0 + 129);
    const auto *skf0_136 = buffer.data(skf0 + 136);
    const auto *skf0_138 = buffer.data(skf0 + 138);
    const auto *skf0_139 = buffer.data(skf0 + 139);
    const auto *skf0_140 = buffer.data(skf0 + 140);
    const auto *skf0_143 = buffer.data(skf0 + 143);
    const auto *skf0_145 = buffer.data(skf0 + 145);
    const auto *skf0_146 = buffer.data(skf0 + 146);
    const auto *skf0_148 = buffer.data(skf0 + 148);
    const auto *skf0_149 = buffer.data(skf0 + 149);
    const auto *skf0_150 = buffer.data(skf0 + 150);
    const auto *skf0_153 = buffer.data(skf0 + 153);
    const auto *skf0_155 = buffer.data(skf0 + 155);
    const auto *skf0_156 = buffer.data(skf0 + 156);
    const auto *skf0_158 = buffer.data(skf0 + 158);
    const auto *skf0_159 = buffer.data(skf0 + 159);
    const auto *skf0_165 = buffer.data(skf0 + 165);
    const auto *skf0_168 = buffer.data(skf0 + 168);
    const auto *skf0_169 = buffer.data(skf0 + 169);
    const auto *skf0_170 = buffer.data(skf0 + 170);

    const auto *skf1_118 = buffer.data(skf1 + 118);
    const auto *skf1_119 = buffer.data(skf1 + 119);
    const auto *skf1_120 = buffer.data(skf1 + 120);
    const auto *skf1_123 = buffer.data(skf1 + 123);
    const auto *skf1_125 = buffer.data(skf1 + 125);
    const auto *skf1_126 = buffer.data(skf1 + 126);
    const auto *skf1_128 = buffer.data(skf1 + 128);
    const auto *skf1_129 = buffer.data(skf1 + 129);
    const auto *skf1_136 = buffer.data(skf1 + 136);
    const auto *skf1_138 = buffer.data(skf1 + 138);
    const auto *skf1_139 = buffer.data(skf1 + 139);
    const auto *skf1_140 = buffer.data(skf1 + 140);
    const auto *skf1_143 = buffer.data(skf1 + 143);
    const auto *skf1_145 = buffer.data(skf1 + 145);
    const auto *skf1_146 = buffer.data(skf1 + 146);
    const auto *skf1_148 = buffer.data(skf1 + 148);
    const auto *skf1_149 = buffer.data(skf1 + 149);
    const auto *skf1_150 = buffer.data(skf1 + 150);
    const auto *skf1_153 = buffer.data(skf1 + 153);
    const auto *skf1_155 = buffer.data(skf1 + 155);
    const auto *skf1_156 = buffer.data(skf1 + 156);
    const auto *skf1_158 = buffer.data(skf1 + 158);
    const auto *skf1_159 = buffer.data(skf1 + 159);
    const auto *skf1_165 = buffer.data(skf1 + 165);
    const auto *skf1_168 = buffer.data(skf1 + 168);
    const auto *skf1_169 = buffer.data(skf1 + 169);
    const auto *skf1_170 = buffer.data(skf1 + 170);

    const auto *skg_175 = buffer.data(skg + 175);
    const auto *skg_177 = buffer.data(skg + 177);
    const auto *skg_178 = buffer.data(skg + 178);
    const auto *skg_179 = buffer.data(skg + 179);
    const auto *skg_180 = buffer.data(skg + 180);
    const auto *skg_182 = buffer.data(skg + 182);
    const auto *skg_183 = buffer.data(skg + 183);
    const auto *skg_185 = buffer.data(skg + 185);
    const auto *skg_186 = buffer.data(skg + 186);
    const auto *skg_189 = buffer.data(skg + 189);
    const auto *skg_190 = buffer.data(skg + 190);
    const auto *skg_191 = buffer.data(skg + 191);
    const auto *skg_192 = buffer.data(skg + 192);
    const auto *skg_193 = buffer.data(skg + 193);
    const auto *skg_194 = buffer.data(skg + 194);
    const auto *skg_195 = buffer.data(skg + 195);
    const auto *skg_197 = buffer.data(skg + 197);
    const auto *skg_198 = buffer.data(skg + 198);
    const auto *skg_200 = buffer.data(skg + 200);
    const auto *skg_205 = buffer.data(skg + 205);
    const auto *skg_206 = buffer.data(skg + 206);
    const auto *skg_207 = buffer.data(skg + 207);
    const auto *skg_208 = buffer.data(skg + 208);
    const auto *skg_209 = buffer.data(skg + 209);
    const auto *skg_210 = buffer.data(skg + 210);
    const auto *skg_212 = buffer.data(skg + 212);
    const auto *skg_213 = buffer.data(skg + 213);
    const auto *skg_215 = buffer.data(skg + 215);
    const auto *skg_216 = buffer.data(skg + 216);
    const auto *skg_219 = buffer.data(skg + 219);
    const auto *skg_220 = buffer.data(skg + 220);
    const auto *skg_221 = buffer.data(skg + 221);
    const auto *skg_222 = buffer.data(skg + 222);
    const auto *skg_223 = buffer.data(skg + 223);
    const auto *skg_224 = buffer.data(skg + 224);
    const auto *skg_225 = buffer.data(skg + 225);
    const auto *skg_227 = buffer.data(skg + 227);
    const auto *skg_228 = buffer.data(skg + 228);
    const auto *skg_230 = buffer.data(skg + 230);
    const auto *skg_231 = buffer.data(skg + 231);
    const auto *skg_234 = buffer.data(skg + 234);
    const auto *skg_235 = buffer.data(skg + 235);
    const auto *skg_236 = buffer.data(skg + 236);
    const auto *skg_237 = buffer.data(skg + 237);
    const auto *skg_238 = buffer.data(skg + 238);
    const auto *skg_239 = buffer.data(skg + 239);
    const auto *skg_240 = buffer.data(skg + 240);
    const auto *skg_242 = buffer.data(skg + 242);
    const auto *skg_243 = buffer.data(skg + 243);
    const auto *skg_245 = buffer.data(skg + 245);
    const auto *skg_249 = buffer.data(skg + 249);
    const auto *skg_250 = buffer.data(skg + 250);
    const auto *skg_251 = buffer.data(skg + 251);
    const auto *skg_252 = buffer.data(skg + 252);
    const auto *skg_253 = buffer.data(skg + 253);
    const auto *skg_254 = buffer.data(skg + 254);
    const auto *skg_255 = buffer.data(skg + 255);

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_z, pc_x, pc_z, sih0_141, sig_177, \
                         sig_178, sig_179, sih1_141, skg_177, skg_178, \
                         skg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_11 * sig_177[k]
                   + f_3 * pc_x[k] * skg_177[k];

        t_244[k] = f_11 * sig_178[k]
                   + f_3 * pc_x[k] * skg_178[k];

        t_245[k] = f_11 * sig_179[k]
                   + f_3 * pc_x[k] * skg_179[k];

        t_246[k] = pb_z[k] * sih0_141[k]
                   - f_8 * pc_z[k] * sih1_141[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, pc_z, sig_100, sig_117, sig_118, skf0_118, \
                         skf0_119, skf1_118, skf1_119, skg_175, skg_177, \
                         skg_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_9 * sig_100[k]
                   + f_3 * pc_z[k] * skg_175[k];

        t_248[k] = f_11 * sig_117[k]
                   + f_4 * skf0_118[k]
                   - f_5 * skf1_118[k]
                   + f_3 * pc_y[k] * skg_177[k];

        t_249[k] = f_11 * sig_118[k]
                   + f_6 * skf0_119[k]
                   - f_7 * skf1_119[k]
                   + f_3 * pc_y[k] * skg_178[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, sig_104, sig_119, sig_180, \
                         skf0_119, skf0_120, skf1_119, skf1_120, skg_179, \
                         skg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * sig_119[k]
                   + f_3 * pc_y[k] * skg_179[k];

        t_251[k] = f_9 * sig_104[k]
                   + f_1 * skf0_119[k]
                   - f_2 * skf1_119[k]
                   + f_3 * pc_z[k] * skg_179[k];

        t_252[k] = f_11 * sig_180[k]
                   + f_1 * skf0_120[k]
                   - f_2 * skf1_120[k]
                   + f_3 * pc_x[k] * skg_180[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, pc_z, sig_105, sig_120, \
                         sig_122, sig_183, skf0_123, skf1_123, skg_180, skg_182, \
                         skg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * sig_120[k]
                   + f_3 * pc_y[k] * skg_180[k];

        t_254[k] = f_10 * sig_105[k]
                   + f_3 * pc_z[k] * skg_180[k];

        t_255[k] = f_11 * sig_183[k]
                   + f_4 * skf0_123[k]
                   - f_5 * skf1_123[k]
                   + f_3 * pc_x[k] * skg_183[k];

        t_256[k] = f_10 * sig_122[k]
                   + f_3 * pc_y[k] * skg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, sig_108, sig_185, sig_186, skf0_125, \
                         skf0_126, skf1_125, skf1_126, skg_183, skg_185, \
                         skg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_11 * sig_185[k]
                   + f_4 * skf0_125[k]
                   - f_5 * skf1_125[k]
                   + f_3 * pc_x[k] * skg_185[k];

        t_258[k] = f_11 * sig_186[k]
                   + f_6 * skf0_126[k]
                   - f_7 * skf1_126[k]
                   + f_3 * pc_x[k] * skg_186[k];

        t_259[k] = f_10 * sig_108[k]
                   + f_3 * pc_z[k] * skg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, sig_125, sig_189, sig_190, \
                         sig_191, skf0_129, skf1_129, skg_185, skg_189, skg_190, \
                         skg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * sig_125[k]
                   + f_3 * pc_y[k] * skg_185[k];

        t_261[k] = f_11 * sig_189[k]
                   + f_6 * skf0_129[k]
                   - f_7 * skf1_129[k]
                   + f_3 * pc_x[k] * skg_189[k];

        t_262[k] = f_11 * sig_190[k]
                   + f_3 * pc_x[k] * skg_190[k];

        t_263[k] = f_11 * sig_191[k]
                   + f_3 * pc_x[k] * skg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, sig_130, sig_192, sig_193, \
                         sig_194, skf0_126, skf1_126, skg_190, skg_192, skg_193, \
                         skg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_11 * sig_192[k]
                   + f_3 * pc_x[k] * skg_192[k];

        t_265[k] = f_11 * sig_193[k]
                   + f_3 * pc_x[k] * skg_193[k];

        t_266[k] = f_11 * sig_194[k]
                   + f_3 * pc_x[k] * skg_194[k];

        t_267[k] = f_10 * sig_130[k]
                   + f_1 * skf0_126[k]
                   - f_2 * skf1_126[k]
                   + f_3 * pc_y[k] * skg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, sig_115, sig_132, sig_133, skf0_128, \
                         skf0_129, skf1_128, skf1_129, skg_190, skg_192, \
                         skg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * sig_115[k]
                   + f_3 * pc_z[k] * skg_190[k];

        t_269[k] = f_10 * sig_132[k]
                   + f_4 * skf0_128[k]
                   - f_5 * skf1_128[k]
                   + f_3 * pc_y[k] * skg_192[k];

        t_270[k] = f_10 * sig_133[k]
                   + f_6 * skf0_129[k]
                   - f_7 * skf1_129[k]
                   + f_3 * pc_y[k] * skg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_y, pc_y, pc_z, sih0_189, sig_119, \
                         sig_134, sig_135, sih1_189, skf0_129, skf1_129, skg_194, \
                         skg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * sig_134[k]
                   + f_3 * pc_y[k] * skg_194[k];

        t_272[k] = f_10 * sig_119[k]
                   + f_1 * skf0_129[k]
                   - f_2 * skf1_129[k]
                   + f_3 * pc_z[k] * skg_194[k];

        t_273[k] = pb_y[k] * sih0_189[k]
                   - f_8 * pc_y[k] * sih1_189[k];

        t_274[k] = f_9 * sig_135[k]
                   + f_3 * pc_y[k] * skg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_y, pc_y, pc_z, sih0_192, sih0_194, \
                         sig_120, sig_136, sig_137, sih1_192, sih1_194, skg_195, \
                         skg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * sig_120[k]
                   + f_3 * pc_z[k] * skg_195[k];

        t_276[k] = pb_y[k] * sih0_192[k]
                   + f_10 * sig_136[k]
                   - f_8 * pc_y[k] * sih1_192[k];

        t_277[k] = f_9 * sig_137[k]
                   + f_3 * pc_y[k] * skg_197[k];

        t_278[k] = pb_y[k] * sih0_194[k]
                   - f_8 * pc_y[k] * sih1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_y, pc_y, pc_z, sih0_195, sih0_198, \
                         sig_123, sig_138, sig_140, sih1_195, sih1_198, skg_198, \
                         skg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_y[k] * sih0_195[k]
                   + f_11 * sig_138[k]
                   - f_8 * pc_y[k] * sih1_195[k];

        t_280[k] = f_11 * sig_123[k]
                   + f_3 * pc_z[k] * skg_198[k];

        t_281[k] = f_9 * sig_140[k]
                   + f_3 * pc_y[k] * skg_200[k];

        t_282[k] = pb_y[k] * sih0_198[k]
                   - f_8 * pc_y[k] * sih1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, sig_205, sig_206, sig_207, \
                         sig_208, sig_209, skg_205, skg_206, skg_207, skg_208, \
                         skg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_11 * sig_205[k]
                   + f_3 * pc_x[k] * skg_205[k];

        t_284[k] = f_11 * sig_206[k]
                   + f_3 * pc_x[k] * skg_206[k];

        t_285[k] = f_11 * sig_207[k]
                   + f_3 * pc_x[k] * skg_207[k];

        t_286[k] = f_11 * sig_208[k]
                   + f_3 * pc_x[k] * skg_208[k];

        t_287[k] = f_11 * sig_209[k]
                   + f_3 * pc_x[k] * skg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, sig_130, sig_145, sig_147, skf0_136, \
                         skf0_138, skf1_136, skf1_138, skg_205, \
                         skg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * sig_145[k]
                   + f_1 * skf0_136[k]
                   - f_2 * skf1_136[k]
                   + f_3 * pc_y[k] * skg_205[k];

        t_289[k] = f_11 * sig_130[k]
                   + f_3 * pc_z[k] * skg_205[k];

        t_290[k] = f_9 * sig_147[k]
                   + f_4 * skf0_138[k]
                   - f_5 * skf1_138[k]
                   + f_3 * pc_y[k] * skg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_y, pc_y, sih0_209, sig_148, sig_149, \
                         sih1_209, skf0_139, skf1_139, skg_208, \
                         skg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * sig_148[k]
                   + f_6 * skf0_139[k]
                   - f_7 * skf1_139[k]
                   + f_3 * pc_y[k] * skg_208[k];

        t_292[k] = f_9 * sig_149[k]
                   + f_3 * pc_y[k] * skg_209[k];

        t_293[k] = pb_y[k] * sih0_209[k]
                   - f_8 * pc_y[k] * sih1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_x, pc_y, pc_z, sig_135, sig_210, \
                         sig_213, skf0_140, skf0_143, skf1_140, skf1_143, skg_210, \
                         skg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_11 * sig_210[k]
                   + f_1 * skf0_140[k]
                   - f_2 * skf1_140[k]
                   + f_3 * pc_x[k] * skg_210[k];

        t_295[k] = f_3 * pc_y[k] * skg_210[k];

        t_296[k] = f_14 * sig_135[k]
                   + f_3 * pc_z[k] * skg_210[k];

        t_297[k] = f_11 * sig_213[k]
                   + f_4 * skf0_143[k]
                   - f_5 * skf1_143[k]
                   + f_3 * pc_x[k] * skg_213[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pc_x, pc_y, sig_215, sig_216, skf0_145, \
                         skf0_146, skf1_145, skf1_146, skg_212, skg_215, \
                         skg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_3 * pc_y[k] * skg_212[k];

        t_299[k] = f_11 * sig_215[k]
                   + f_4 * skf0_145[k]
                   - f_5 * skf1_145[k]
                   + f_3 * pc_x[k] * skg_215[k];

        t_300[k] = f_11 * sig_216[k]
                   + f_6 * skf0_146[k]
                   - f_7 * skf1_146[k]
                   + f_3 * pc_x[k] * skg_216[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_x, pc_y, pc_z, sig_138, sig_219, \
                         sig_220, skf0_149, skf1_149, skg_213, skg_215, skg_219, \
                         skg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * sig_138[k]
                   + f_3 * pc_z[k] * skg_213[k];

        t_302[k] = f_3 * pc_y[k] * skg_215[k];

        t_303[k] = f_11 * sig_219[k]
                   + f_6 * skf0_149[k]
                   - f_7 * skf1_149[k]
                   + f_3 * pc_x[k] * skg_219[k];

        t_304[k] = f_11 * sig_220[k]
                   + f_3 * pc_x[k] * skg_220[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, sig_221, sig_222, sig_223, sig_224, \
                         skg_221, skg_222, skg_223, skg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_11 * sig_221[k]
                   + f_3 * pc_x[k] * skg_221[k];

        t_306[k] = f_11 * sig_222[k]
                   + f_3 * pc_x[k] * skg_222[k];

        t_307[k] = f_11 * sig_223[k]
                   + f_3 * pc_x[k] * skg_223[k];

        t_308[k] = f_11 * sig_224[k]
                   + f_3 * pc_x[k] * skg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pc_y, pc_z, sig_145, skf0_146, skf0_148, \
                         skf0_149, skf1_146, skf1_148, skf1_149, skg_220, skg_222, \
                         skg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * skf0_146[k]
                   - f_2 * skf1_146[k]
                   + f_3 * pc_y[k] * skg_220[k];

        t_310[k] = f_14 * sig_145[k]
                   + f_3 * pc_z[k] * skg_220[k];

        t_311[k] = f_4 * skf0_148[k]
                   - f_5 * skf1_148[k]
                   + f_3 * pc_y[k] * skg_222[k];

        t_312[k] = f_6 * skf0_149[k]
                   - f_7 * skf1_149[k]
                   + f_3 * pc_y[k] * skg_223[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pc_x, pc_y, pc_z, sig_149, sig_150, \
                         sig_225, skf0_149, skf0_150, skf1_149, skf1_150, skg_224, \
                         skg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * skg_224[k];

        t_314[k] = f_14 * sig_149[k]
                   + f_1 * skf0_149[k]
                   - f_2 * skf1_149[k]
                   + f_3 * pc_z[k] * skg_224[k];

        t_315[k] = f_10 * sig_225[k]
                   + f_1 * skf0_150[k]
                   - f_2 * skf1_150[k]
                   + f_3 * pc_x[k] * skg_225[k];

        t_316[k] = f_13 * sig_150[k]
                   + f_3 * pc_y[k] * skg_225[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pc_x, pc_y, pc_z, sig_152, sig_228, skf0_153, \
                         skf1_153, skg_225, skg_227, skg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_3 * pc_z[k] * skg_225[k];

        t_318[k] = f_10 * sig_228[k]
                   + f_4 * skf0_153[k]
                   - f_5 * skf1_153[k]
                   + f_3 * pc_x[k] * skg_228[k];

        t_319[k] = f_13 * sig_152[k]
                   + f_3 * pc_y[k] * skg_227[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pc_x, pc_z, sig_230, sig_231, skf0_155, \
                         skf0_156, skf1_155, skf1_156, skg_228, skg_230, \
                         skg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_10 * sig_230[k]
                   + f_4 * skf0_155[k]
                   - f_5 * skf1_155[k]
                   + f_3 * pc_x[k] * skg_230[k];

        t_321[k] = f_10 * sig_231[k]
                   + f_6 * skf0_156[k]
                   - f_7 * skf1_156[k]
                   + f_3 * pc_x[k] * skg_231[k];

        t_322[k] = f_3 * pc_z[k] * skg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, sig_155, sig_234, sig_235, \
                         sig_236, skf0_159, skf1_159, skg_230, skg_234, skg_235, \
                         skg_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_13 * sig_155[k]
                   + f_3 * pc_y[k] * skg_230[k];

        t_324[k] = f_10 * sig_234[k]
                   + f_6 * skf0_159[k]
                   - f_7 * skf1_159[k]
                   + f_3 * pc_x[k] * skg_234[k];

        t_325[k] = f_10 * sig_235[k]
                   + f_3 * pc_x[k] * skg_235[k];

        t_326[k] = f_10 * sig_236[k]
                   + f_3 * pc_x[k] * skg_236[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, sig_160, sig_237, sig_238, \
                         sig_239, skf0_156, skf1_156, skg_235, skg_237, skg_238, \
                         skg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_10 * sig_237[k]
                   + f_3 * pc_x[k] * skg_237[k];

        t_328[k] = f_10 * sig_238[k]
                   + f_3 * pc_x[k] * skg_238[k];

        t_329[k] = f_10 * sig_239[k]
                   + f_3 * pc_x[k] * skg_239[k];

        t_330[k] = f_13 * sig_160[k]
                   + f_1 * skf0_156[k]
                   - f_2 * skf1_156[k]
                   + f_3 * pc_y[k] * skg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_y, pc_z, sig_162, sig_163, skf0_158, \
                         skf0_159, skf1_158, skf1_159, skg_235, skg_237, \
                         skg_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * skg_235[k];

        t_332[k] = f_13 * sig_162[k]
                   + f_4 * skf0_158[k]
                   - f_5 * skf1_158[k]
                   + f_3 * pc_y[k] * skg_237[k];

        t_333[k] = f_13 * sig_163[k]
                   + f_6 * skf0_159[k]
                   - f_7 * skf1_159[k]
                   + f_3 * pc_y[k] * skg_238[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pb_z, pc_y, pc_z, sih0_210, sig_164, \
                         sig_165, sih1_210, skf0_159, skf1_159, skg_239, \
                         skg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_13 * sig_164[k]
                   + f_3 * pc_y[k] * skg_239[k];

        t_335[k] = f_1 * skf0_159[k]
                   - f_2 * skf1_159[k]
                   + f_3 * pc_z[k] * skg_239[k];

        t_336[k] = pb_z[k] * sih0_210[k]
                   - f_8 * pc_z[k] * sih1_210[k];

        t_337[k] = f_14 * sig_165[k]
                   + f_3 * pc_y[k] * skg_240[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_y, pc_z, sih0_213, sig_150, sig_167, \
                         sih1_213, skg_240, skg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_9 * sig_150[k]
                   + f_3 * pc_z[k] * skg_240[k];

        t_339[k] = pb_z[k] * sih0_213[k]
                   - f_8 * pc_z[k] * sih1_213[k];

        t_340[k] = f_14 * sig_167[k]
                   + f_3 * pc_y[k] * skg_242[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pb_z, pc_x, pc_z, sih0_216, sig_153, sig_245, \
                         sih1_216, skf0_165, skf1_165, skg_243, \
                         skg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_10 * sig_245[k]
                   + f_4 * skf0_165[k]
                   - f_5 * skf1_165[k]
                   + f_3 * pc_x[k] * skg_245[k];

        t_342[k] = pb_z[k] * sih0_216[k]
                   - f_8 * pc_z[k] * sih1_216[k];

        t_343[k] = f_9 * sig_153[k]
                   + f_3 * pc_z[k] * skg_243[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, sig_170, sig_249, sig_250, \
                         sig_251, skf0_169, skf1_169, skg_245, skg_249, skg_250, \
                         skg_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_14 * sig_170[k]
                   + f_3 * pc_y[k] * skg_245[k];

        t_345[k] = f_10 * sig_249[k]
                   + f_6 * skf0_169[k]
                   - f_7 * skf1_169[k]
                   + f_3 * pc_x[k] * skg_249[k];

        t_346[k] = f_10 * sig_250[k]
                   + f_3 * pc_x[k] * skg_250[k];

        t_347[k] = f_10 * sig_251[k]
                   + f_3 * pc_x[k] * skg_251[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pb_z, pc_x, pc_z, sih0_225, sig_252, \
                         sig_253, sig_254, sih1_225, skg_252, skg_253, \
                         skg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_10 * sig_252[k]
                   + f_3 * pc_x[k] * skg_252[k];

        t_349[k] = f_10 * sig_253[k]
                   + f_3 * pc_x[k] * skg_253[k];

        t_350[k] = f_10 * sig_254[k]
                   + f_3 * pc_x[k] * skg_254[k];

        t_351[k] = pb_z[k] * sih0_225[k]
                   - f_8 * pc_z[k] * sih1_225[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pc_y, pc_z, sig_160, sig_177, sig_178, skf0_168, \
                         skf0_169, skf1_168, skf1_169, skg_250, skg_252, \
                         skg_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_9 * sig_160[k]
                   + f_3 * pc_z[k] * skg_250[k];

        t_353[k] = f_14 * sig_177[k]
                   + f_4 * skf0_168[k]
                   - f_5 * skf1_168[k]
                   + f_3 * pc_y[k] * skg_252[k];

        t_354[k] = f_14 * sig_178[k]
                   + f_6 * skf0_169[k]
                   - f_7 * skf1_169[k]
                   + f_3 * pc_y[k] * skg_253[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_x, pc_y, pc_z, sig_164, sig_179, sig_255, \
                         skf0_169, skf0_170, skf1_169, skf1_170, skg_254, \
                         skg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_14 * sig_179[k]
                   + f_3 * pc_y[k] * skg_254[k];

        t_356[k] = f_9 * sig_164[k]
                   + f_1 * skf0_169[k]
                   - f_2 * skf1_169[k]
                   + f_3 * pc_z[k] * skg_254[k];

        t_357[k] = f_10 * sig_255[k]
                   + f_1 * skf0_170[k]
                   - f_2 * skf1_170[k]
                   + f_3 * pc_x[k] * skg_255[k];
    }
}

static auto
compute_prim_skh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sih0,
                                                          const size_t sig, const size_t sih1,
                                                          const size_t skf0, const size_t skf1,
                                                          const size_t skg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sih0_294 = buffer.data(sih0 + 294);
    const auto *sih0_297 = buffer.data(sih0 + 297);
    const auto *sih0_299 = buffer.data(sih0 + 299);
    const auto *sih0_300 = buffer.data(sih0 + 300);
    const auto *sih0_303 = buffer.data(sih0 + 303);
    const auto *sih0_314 = buffer.data(sih0 + 314);
    const auto *sih0_315 = buffer.data(sih0 + 315);
    const auto *sih0_318 = buffer.data(sih0 + 318);
    const auto *sih0_321 = buffer.data(sih0 + 321);
    const auto *sih0_441 = buffer.data(sih0 + 441);
    const auto *sih0_444 = buffer.data(sih0 + 444);
    const auto *sih0_446 = buffer.data(sih0 + 446);
    const auto *sih0_447 = buffer.data(sih0 + 447);
    const auto *sih0_450 = buffer.data(sih0 + 450);
    const auto *sih0_456 = buffer.data(sih0 + 456);
    const auto *sih0_458 = buffer.data(sih0 + 458);
    const auto *sih0_459 = buffer.data(sih0 + 459);
    const auto *sih0_461 = buffer.data(sih0 + 461);
    const auto *sih0_467 = buffer.data(sih0 + 467);
    const auto *sih0_471 = buffer.data(sih0 + 471);

    const auto *sig_165 = buffer.data(sig + 165);
    const auto *sig_168 = buffer.data(sig + 168);
    const auto *sig_175 = buffer.data(sig + 175);
    const auto *sig_179 = buffer.data(sig + 179);
    const auto *sig_180 = buffer.data(sig + 180);
    const auto *sig_182 = buffer.data(sig + 182);
    const auto *sig_183 = buffer.data(sig + 183);
    const auto *sig_185 = buffer.data(sig + 185);
    const auto *sig_190 = buffer.data(sig + 190);
    const auto *sig_192 = buffer.data(sig + 192);
    const auto *sig_193 = buffer.data(sig + 193);
    const auto *sig_194 = buffer.data(sig + 194);
    const auto *sig_195 = buffer.data(sig + 195);
    const auto *sig_197 = buffer.data(sig + 197);
    const auto *sig_198 = buffer.data(sig + 198);
    const auto *sig_200 = buffer.data(sig + 200);
    const auto *sig_205 = buffer.data(sig + 205);
    const auto *sig_207 = buffer.data(sig + 207);
    const auto *sig_208 = buffer.data(sig + 208);
    const auto *sig_209 = buffer.data(sig + 209);
    const auto *sig_210 = buffer.data(sig + 210);
    const auto *sig_211 = buffer.data(sig + 211);
    const auto *sig_212 = buffer.data(sig + 212);
    const auto *sig_213 = buffer.data(sig + 213);
    const auto *sig_215 = buffer.data(sig + 215);
    const auto *sig_220 = buffer.data(sig + 220);
    const auto *sig_222 = buffer.data(sig + 222);
    const auto *sig_223 = buffer.data(sig + 223);
    const auto *sig_224 = buffer.data(sig + 224);
    const auto *sig_225 = buffer.data(sig + 225);
    const auto *sig_227 = buffer.data(sig + 227);
    const auto *sig_228 = buffer.data(sig + 228);
    const auto *sig_230 = buffer.data(sig + 230);
    const auto *sig_239 = buffer.data(sig + 239);
    const auto *sig_240 = buffer.data(sig + 240);
    const auto *sig_242 = buffer.data(sig + 242);
    const auto *sig_245 = buffer.data(sig + 245);
    const auto *sig_258 = buffer.data(sig + 258);
    const auto *sig_260 = buffer.data(sig + 260);
    const auto *sig_261 = buffer.data(sig + 261);
    const auto *sig_264 = buffer.data(sig + 264);
    const auto *sig_265 = buffer.data(sig + 265);
    const auto *sig_266 = buffer.data(sig + 266);
    const auto *sig_267 = buffer.data(sig + 267);
    const auto *sig_268 = buffer.data(sig + 268);
    const auto *sig_269 = buffer.data(sig + 269);
    const auto *sig_270 = buffer.data(sig + 270);
    const auto *sig_273 = buffer.data(sig + 273);
    const auto *sig_275 = buffer.data(sig + 275);
    const auto *sig_276 = buffer.data(sig + 276);
    const auto *sig_279 = buffer.data(sig + 279);
    const auto *sig_280 = buffer.data(sig + 280);
    const auto *sig_281 = buffer.data(sig + 281);
    const auto *sig_282 = buffer.data(sig + 282);
    const auto *sig_283 = buffer.data(sig + 283);
    const auto *sig_284 = buffer.data(sig + 284);
    const auto *sig_295 = buffer.data(sig + 295);
    const auto *sig_296 = buffer.data(sig + 296);
    const auto *sig_297 = buffer.data(sig + 297);
    const auto *sig_298 = buffer.data(sig + 298);
    const auto *sig_299 = buffer.data(sig + 299);
    const auto *sig_300 = buffer.data(sig + 300);
    const auto *sig_303 = buffer.data(sig + 303);
    const auto *sig_305 = buffer.data(sig + 305);
    const auto *sig_306 = buffer.data(sig + 306);
    const auto *sig_309 = buffer.data(sig + 309);
    const auto *sig_310 = buffer.data(sig + 310);
    const auto *sig_311 = buffer.data(sig + 311);
    const auto *sig_312 = buffer.data(sig + 312);
    const auto *sig_313 = buffer.data(sig + 313);
    const auto *sig_314 = buffer.data(sig + 314);
    const auto *sig_315 = buffer.data(sig + 315);
    const auto *sig_318 = buffer.data(sig + 318);
    const auto *sig_320 = buffer.data(sig + 320);
    const auto *sig_321 = buffer.data(sig + 321);
    const auto *sig_324 = buffer.data(sig + 324);
    const auto *sig_325 = buffer.data(sig + 325);
    const auto *sig_326 = buffer.data(sig + 326);
    const auto *sig_327 = buffer.data(sig + 327);
    const auto *sig_328 = buffer.data(sig + 328);
    const auto *sig_329 = buffer.data(sig + 329);
    const auto *sig_335 = buffer.data(sig + 335);
    const auto *sig_339 = buffer.data(sig + 339);
    const auto *sig_340 = buffer.data(sig + 340);
    const auto *sig_341 = buffer.data(sig + 341);
    const auto *sig_342 = buffer.data(sig + 342);

    const auto *sih1_294 = buffer.data(sih1 + 294);
    const auto *sih1_297 = buffer.data(sih1 + 297);
    const auto *sih1_299 = buffer.data(sih1 + 299);
    const auto *sih1_300 = buffer.data(sih1 + 300);
    const auto *sih1_303 = buffer.data(sih1 + 303);
    const auto *sih1_314 = buffer.data(sih1 + 314);
    const auto *sih1_315 = buffer.data(sih1 + 315);
    const auto *sih1_318 = buffer.data(sih1 + 318);
    const auto *sih1_321 = buffer.data(sih1 + 321);
    const auto *sih1_441 = buffer.data(sih1 + 441);
    const auto *sih1_444 = buffer.data(sih1 + 444);
    const auto *sih1_446 = buffer.data(sih1 + 446);
    const auto *sih1_447 = buffer.data(sih1 + 447);
    const auto *sih1_450 = buffer.data(sih1 + 450);
    const auto *sih1_456 = buffer.data(sih1 + 456);
    const auto *sih1_458 = buffer.data(sih1 + 458);
    const auto *sih1_459 = buffer.data(sih1 + 459);
    const auto *sih1_461 = buffer.data(sih1 + 461);
    const auto *sih1_467 = buffer.data(sih1 + 467);
    const auto *sih1_471 = buffer.data(sih1 + 471);

    const auto *skf0_173 = buffer.data(skf0 + 173);
    const auto *skf0_175 = buffer.data(skf0 + 175);
    const auto *skf0_176 = buffer.data(skf0 + 176);
    const auto *skf0_178 = buffer.data(skf0 + 178);
    const auto *skf0_179 = buffer.data(skf0 + 179);
    const auto *skf0_180 = buffer.data(skf0 + 180);
    const auto *skf0_183 = buffer.data(skf0 + 183);
    const auto *skf0_185 = buffer.data(skf0 + 185);
    const auto *skf0_186 = buffer.data(skf0 + 186);
    const auto *skf0_188 = buffer.data(skf0 + 188);
    const auto *skf0_189 = buffer.data(skf0 + 189);
    const auto *skf0_196 = buffer.data(skf0 + 196);
    const auto *skf0_198 = buffer.data(skf0 + 198);
    const auto *skf0_199 = buffer.data(skf0 + 199);
    const auto *skf0_200 = buffer.data(skf0 + 200);
    const auto *skf0_203 = buffer.data(skf0 + 203);
    const auto *skf0_205 = buffer.data(skf0 + 205);
    const auto *skf0_206 = buffer.data(skf0 + 206);
    const auto *skf0_208 = buffer.data(skf0 + 208);
    const auto *skf0_209 = buffer.data(skf0 + 209);

    const auto *skf1_173 = buffer.data(skf1 + 173);
    const auto *skf1_175 = buffer.data(skf1 + 175);
    const auto *skf1_176 = buffer.data(skf1 + 176);
    const auto *skf1_178 = buffer.data(skf1 + 178);
    const auto *skf1_179 = buffer.data(skf1 + 179);
    const auto *skf1_180 = buffer.data(skf1 + 180);
    const auto *skf1_183 = buffer.data(skf1 + 183);
    const auto *skf1_185 = buffer.data(skf1 + 185);
    const auto *skf1_186 = buffer.data(skf1 + 186);
    const auto *skf1_188 = buffer.data(skf1 + 188);
    const auto *skf1_189 = buffer.data(skf1 + 189);
    const auto *skf1_196 = buffer.data(skf1 + 196);
    const auto *skf1_198 = buffer.data(skf1 + 198);
    const auto *skf1_199 = buffer.data(skf1 + 199);
    const auto *skf1_200 = buffer.data(skf1 + 200);
    const auto *skf1_203 = buffer.data(skf1 + 203);
    const auto *skf1_205 = buffer.data(skf1 + 205);
    const auto *skf1_206 = buffer.data(skf1 + 206);
    const auto *skf1_208 = buffer.data(skf1 + 208);
    const auto *skf1_209 = buffer.data(skf1 + 209);

    const auto *skg_255 = buffer.data(skg + 255);
    const auto *skg_257 = buffer.data(skg + 257);
    const auto *skg_258 = buffer.data(skg + 258);
    const auto *skg_260 = buffer.data(skg + 260);
    const auto *skg_261 = buffer.data(skg + 261);
    const auto *skg_264 = buffer.data(skg + 264);
    const auto *skg_265 = buffer.data(skg + 265);
    const auto *skg_266 = buffer.data(skg + 266);
    const auto *skg_267 = buffer.data(skg + 267);
    const auto *skg_268 = buffer.data(skg + 268);
    const auto *skg_269 = buffer.data(skg + 269);
    const auto *skg_270 = buffer.data(skg + 270);
    const auto *skg_272 = buffer.data(skg + 272);
    const auto *skg_273 = buffer.data(skg + 273);
    const auto *skg_275 = buffer.data(skg + 275);
    const auto *skg_276 = buffer.data(skg + 276);
    const auto *skg_279 = buffer.data(skg + 279);
    const auto *skg_280 = buffer.data(skg + 280);
    const auto *skg_281 = buffer.data(skg + 281);
    const auto *skg_282 = buffer.data(skg + 282);
    const auto *skg_283 = buffer.data(skg + 283);
    const auto *skg_284 = buffer.data(skg + 284);
    const auto *skg_285 = buffer.data(skg + 285);
    const auto *skg_287 = buffer.data(skg + 287);
    const auto *skg_288 = buffer.data(skg + 288);
    const auto *skg_290 = buffer.data(skg + 290);
    const auto *skg_295 = buffer.data(skg + 295);
    const auto *skg_296 = buffer.data(skg + 296);
    const auto *skg_297 = buffer.data(skg + 297);
    const auto *skg_298 = buffer.data(skg + 298);
    const auto *skg_299 = buffer.data(skg + 299);
    const auto *skg_300 = buffer.data(skg + 300);
    const auto *skg_302 = buffer.data(skg + 302);
    const auto *skg_303 = buffer.data(skg + 303);
    const auto *skg_305 = buffer.data(skg + 305);
    const auto *skg_306 = buffer.data(skg + 306);
    const auto *skg_309 = buffer.data(skg + 309);
    const auto *skg_310 = buffer.data(skg + 310);
    const auto *skg_311 = buffer.data(skg + 311);
    const auto *skg_312 = buffer.data(skg + 312);
    const auto *skg_313 = buffer.data(skg + 313);
    const auto *skg_314 = buffer.data(skg + 314);
    const auto *skg_315 = buffer.data(skg + 315);
    const auto *skg_317 = buffer.data(skg + 317);
    const auto *skg_318 = buffer.data(skg + 318);
    const auto *skg_320 = buffer.data(skg + 320);
    const auto *skg_325 = buffer.data(skg + 325);
    const auto *skg_326 = buffer.data(skg + 326);
    const auto *skg_327 = buffer.data(skg + 327);
    const auto *skg_328 = buffer.data(skg + 328);
    const auto *skg_329 = buffer.data(skg + 329);
    const auto *skg_330 = buffer.data(skg + 330);
    const auto *skg_332 = buffer.data(skg + 332);
    const auto *skg_333 = buffer.data(skg + 333);
    const auto *skg_335 = buffer.data(skg + 335);
    const auto *skg_340 = buffer.data(skg + 340);
    const auto *skg_341 = buffer.data(skg + 341);
    const auto *skg_342 = buffer.data(skg + 342);

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pc_x, pc_y, pc_z, sig_165, sig_180, \
                         sig_182, sig_258, skf0_173, skf1_173, skg_255, skg_257, \
                         skg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_11 * sig_180[k]
                   + f_3 * pc_y[k] * skg_255[k];

        t_359[k] = f_10 * sig_165[k]
                   + f_3 * pc_z[k] * skg_255[k];

        t_360[k] = f_10 * sig_258[k]
                   + f_4 * skf0_173[k]
                   - f_5 * skf1_173[k]
                   + f_3 * pc_x[k] * skg_258[k];

        t_361[k] = f_11 * sig_182[k]
                   + f_3 * pc_y[k] * skg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, sig_168, sig_260, sig_261, skf0_175, \
                         skf0_176, skf1_175, skf1_176, skg_258, skg_260, \
                         skg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_10 * sig_260[k]
                   + f_4 * skf0_175[k]
                   - f_5 * skf1_175[k]
                   + f_3 * pc_x[k] * skg_260[k];

        t_363[k] = f_10 * sig_261[k]
                   + f_6 * skf0_176[k]
                   - f_7 * skf1_176[k]
                   + f_3 * pc_x[k] * skg_261[k];

        t_364[k] = f_10 * sig_168[k]
                   + f_3 * pc_z[k] * skg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, sig_185, sig_264, sig_265, \
                         sig_266, skf0_179, skf1_179, skg_260, skg_264, skg_265, \
                         skg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * sig_185[k]
                   + f_3 * pc_y[k] * skg_260[k];

        t_366[k] = f_10 * sig_264[k]
                   + f_6 * skf0_179[k]
                   - f_7 * skf1_179[k]
                   + f_3 * pc_x[k] * skg_264[k];

        t_367[k] = f_10 * sig_265[k]
                   + f_3 * pc_x[k] * skg_265[k];

        t_368[k] = f_10 * sig_266[k]
                   + f_3 * pc_x[k] * skg_266[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, sig_190, sig_267, sig_268, \
                         sig_269, skf0_176, skf1_176, skg_265, skg_267, skg_268, \
                         skg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_10 * sig_267[k]
                   + f_3 * pc_x[k] * skg_267[k];

        t_370[k] = f_10 * sig_268[k]
                   + f_3 * pc_x[k] * skg_268[k];

        t_371[k] = f_10 * sig_269[k]
                   + f_3 * pc_x[k] * skg_269[k];

        t_372[k] = f_11 * sig_190[k]
                   + f_1 * skf0_176[k]
                   - f_2 * skf1_176[k]
                   + f_3 * pc_y[k] * skg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, sig_175, sig_192, sig_193, skf0_178, \
                         skf0_179, skf1_178, skf1_179, skg_265, skg_267, \
                         skg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * sig_175[k]
                   + f_3 * pc_z[k] * skg_265[k];

        t_374[k] = f_11 * sig_192[k]
                   + f_4 * skf0_178[k]
                   - f_5 * skf1_178[k]
                   + f_3 * pc_y[k] * skg_267[k];

        t_375[k] = f_11 * sig_193[k]
                   + f_6 * skf0_179[k]
                   - f_7 * skf1_179[k]
                   + f_3 * pc_y[k] * skg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, sig_179, sig_194, sig_270, \
                         skf0_179, skf0_180, skf1_179, skf1_180, skg_269, \
                         skg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * sig_194[k]
                   + f_3 * pc_y[k] * skg_269[k];

        t_377[k] = f_10 * sig_179[k]
                   + f_1 * skf0_179[k]
                   - f_2 * skf1_179[k]
                   + f_3 * pc_z[k] * skg_269[k];

        t_378[k] = f_10 * sig_270[k]
                   + f_1 * skf0_180[k]
                   - f_2 * skf1_180[k]
                   + f_3 * pc_x[k] * skg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, sig_180, sig_195, \
                         sig_197, sig_273, skf0_183, skf1_183, skg_270, skg_272, \
                         skg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * sig_195[k]
                   + f_3 * pc_y[k] * skg_270[k];

        t_380[k] = f_11 * sig_180[k]
                   + f_3 * pc_z[k] * skg_270[k];

        t_381[k] = f_10 * sig_273[k]
                   + f_4 * skf0_183[k]
                   - f_5 * skf1_183[k]
                   + f_3 * pc_x[k] * skg_273[k];

        t_382[k] = f_10 * sig_197[k]
                   + f_3 * pc_y[k] * skg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, sig_183, sig_275, sig_276, skf0_185, \
                         skf0_186, skf1_185, skf1_186, skg_273, skg_275, \
                         skg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_10 * sig_275[k]
                   + f_4 * skf0_185[k]
                   - f_5 * skf1_185[k]
                   + f_3 * pc_x[k] * skg_275[k];

        t_384[k] = f_10 * sig_276[k]
                   + f_6 * skf0_186[k]
                   - f_7 * skf1_186[k]
                   + f_3 * pc_x[k] * skg_276[k];

        t_385[k] = f_11 * sig_183[k]
                   + f_3 * pc_z[k] * skg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, sig_200, sig_279, sig_280, \
                         sig_281, skf0_189, skf1_189, skg_275, skg_279, skg_280, \
                         skg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * sig_200[k]
                   + f_3 * pc_y[k] * skg_275[k];

        t_387[k] = f_10 * sig_279[k]
                   + f_6 * skf0_189[k]
                   - f_7 * skf1_189[k]
                   + f_3 * pc_x[k] * skg_279[k];

        t_388[k] = f_10 * sig_280[k]
                   + f_3 * pc_x[k] * skg_280[k];

        t_389[k] = f_10 * sig_281[k]
                   + f_3 * pc_x[k] * skg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, sig_205, sig_282, sig_283, \
                         sig_284, skf0_186, skf1_186, skg_280, skg_282, skg_283, \
                         skg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_10 * sig_282[k]
                   + f_3 * pc_x[k] * skg_282[k];

        t_391[k] = f_10 * sig_283[k]
                   + f_3 * pc_x[k] * skg_283[k];

        t_392[k] = f_10 * sig_284[k]
                   + f_3 * pc_x[k] * skg_284[k];

        t_393[k] = f_10 * sig_205[k]
                   + f_1 * skf0_186[k]
                   - f_2 * skf1_186[k]
                   + f_3 * pc_y[k] * skg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, sig_190, sig_207, sig_208, skf0_188, \
                         skf0_189, skf1_188, skf1_189, skg_280, skg_282, \
                         skg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * sig_190[k]
                   + f_3 * pc_z[k] * skg_280[k];

        t_395[k] = f_10 * sig_207[k]
                   + f_4 * skf0_188[k]
                   - f_5 * skf1_188[k]
                   + f_3 * pc_y[k] * skg_282[k];

        t_396[k] = f_10 * sig_208[k]
                   + f_6 * skf0_189[k]
                   - f_7 * skf1_189[k]
                   + f_3 * pc_y[k] * skg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pb_y, pc_y, pc_z, sih0_294, sig_194, \
                         sig_209, sig_210, sih1_294, skf0_189, skf1_189, skg_284, \
                         skg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * sig_209[k]
                   + f_3 * pc_y[k] * skg_284[k];

        t_398[k] = f_11 * sig_194[k]
                   + f_1 * skf0_189[k]
                   - f_2 * skf1_189[k]
                   + f_3 * pc_z[k] * skg_284[k];

        t_399[k] = pb_y[k] * sih0_294[k]
                   - f_8 * pc_y[k] * sih1_294[k];

        t_400[k] = f_9 * sig_210[k]
                   + f_3 * pc_y[k] * skg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pb_y, pc_y, pc_z, sih0_297, sih0_299, \
                         sig_195, sig_211, sig_212, sih1_297, sih1_299, skg_285, \
                         skg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_14 * sig_195[k]
                   + f_3 * pc_z[k] * skg_285[k];

        t_402[k] = pb_y[k] * sih0_297[k]
                   + f_10 * sig_211[k]
                   - f_8 * pc_y[k] * sih1_297[k];

        t_403[k] = f_9 * sig_212[k]
                   + f_3 * pc_y[k] * skg_287[k];

        t_404[k] = pb_y[k] * sih0_299[k]
                   - f_8 * pc_y[k] * sih1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pb_y, pc_y, pc_z, sih0_300, sih0_303, \
                         sig_198, sig_213, sig_215, sih1_300, sih1_303, skg_288, \
                         skg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pb_y[k] * sih0_300[k]
                   + f_11 * sig_213[k]
                   - f_8 * pc_y[k] * sih1_300[k];

        t_406[k] = f_14 * sig_198[k]
                   + f_3 * pc_z[k] * skg_288[k];

        t_407[k] = f_9 * sig_215[k]
                   + f_3 * pc_y[k] * skg_290[k];

        t_408[k] = pb_y[k] * sih0_303[k]
                   - f_8 * pc_y[k] * sih1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, sig_295, sig_296, sig_297, \
                         sig_298, sig_299, skg_295, skg_296, skg_297, skg_298, \
                         skg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_10 * sig_295[k]
                   + f_3 * pc_x[k] * skg_295[k];

        t_410[k] = f_10 * sig_296[k]
                   + f_3 * pc_x[k] * skg_296[k];

        t_411[k] = f_10 * sig_297[k]
                   + f_3 * pc_x[k] * skg_297[k];

        t_412[k] = f_10 * sig_298[k]
                   + f_3 * pc_x[k] * skg_298[k];

        t_413[k] = f_10 * sig_299[k]
                   + f_3 * pc_x[k] * skg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, sig_205, sig_220, sig_222, skf0_196, \
                         skf0_198, skf1_196, skf1_198, skg_295, \
                         skg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * sig_220[k]
                   + f_1 * skf0_196[k]
                   - f_2 * skf1_196[k]
                   + f_3 * pc_y[k] * skg_295[k];

        t_415[k] = f_14 * sig_205[k]
                   + f_3 * pc_z[k] * skg_295[k];

        t_416[k] = f_9 * sig_222[k]
                   + f_4 * skf0_198[k]
                   - f_5 * skf1_198[k]
                   + f_3 * pc_y[k] * skg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pc_y, sih0_314, sig_223, sig_224, \
                         sih1_314, skf0_199, skf1_199, skg_298, \
                         skg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * sig_223[k]
                   + f_6 * skf0_199[k]
                   - f_7 * skf1_199[k]
                   + f_3 * pc_y[k] * skg_298[k];

        t_418[k] = f_9 * sig_224[k]
                   + f_3 * pc_y[k] * skg_299[k];

        t_419[k] = pb_y[k] * sih0_314[k]
                   - f_8 * pc_y[k] * sih1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, sig_210, sig_300, \
                         sig_303, skf0_200, skf0_203, skf1_200, skf1_203, skg_300, \
                         skg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_10 * sig_300[k]
                   + f_1 * skf0_200[k]
                   - f_2 * skf1_200[k]
                   + f_3 * pc_x[k] * skg_300[k];

        t_421[k] = f_3 * pc_y[k] * skg_300[k];

        t_422[k] = f_13 * sig_210[k]
                   + f_3 * pc_z[k] * skg_300[k];

        t_423[k] = f_10 * sig_303[k]
                   + f_4 * skf0_203[k]
                   - f_5 * skf1_203[k]
                   + f_3 * pc_x[k] * skg_303[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, sig_305, sig_306, skf0_205, \
                         skf0_206, skf1_205, skf1_206, skg_302, skg_305, \
                         skg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_y[k] * skg_302[k];

        t_425[k] = f_10 * sig_305[k]
                   + f_4 * skf0_205[k]
                   - f_5 * skf1_205[k]
                   + f_3 * pc_x[k] * skg_305[k];

        t_426[k] = f_10 * sig_306[k]
                   + f_6 * skf0_206[k]
                   - f_7 * skf1_206[k]
                   + f_3 * pc_x[k] * skg_306[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pc_x, pc_y, pc_z, sig_213, sig_309, \
                         sig_310, skf0_209, skf1_209, skg_303, skg_305, skg_309, \
                         skg_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_13 * sig_213[k]
                   + f_3 * pc_z[k] * skg_303[k];

        t_428[k] = f_3 * pc_y[k] * skg_305[k];

        t_429[k] = f_10 * sig_309[k]
                   + f_6 * skf0_209[k]
                   - f_7 * skf1_209[k]
                   + f_3 * pc_x[k] * skg_309[k];

        t_430[k] = f_10 * sig_310[k]
                   + f_3 * pc_x[k] * skg_310[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, sig_311, sig_312, sig_313, sig_314, \
                         skg_311, skg_312, skg_313, skg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_10 * sig_311[k]
                   + f_3 * pc_x[k] * skg_311[k];

        t_432[k] = f_10 * sig_312[k]
                   + f_3 * pc_x[k] * skg_312[k];

        t_433[k] = f_10 * sig_313[k]
                   + f_3 * pc_x[k] * skg_313[k];

        t_434[k] = f_10 * sig_314[k]
                   + f_3 * pc_x[k] * skg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pc_y, pc_z, sig_220, skf0_206, skf0_208, \
                         skf0_209, skf1_206, skf1_208, skf1_209, skg_310, skg_312, \
                         skg_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * skf0_206[k]
                   - f_2 * skf1_206[k]
                   + f_3 * pc_y[k] * skg_310[k];

        t_436[k] = f_13 * sig_220[k]
                   + f_3 * pc_z[k] * skg_310[k];

        t_437[k] = f_4 * skf0_208[k]
                   - f_5 * skf1_208[k]
                   + f_3 * pc_y[k] * skg_312[k];

        t_438[k] = f_6 * skf0_209[k]
                   - f_7 * skf1_209[k]
                   + f_3 * pc_y[k] * skg_313[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, pb_x, pc_x, pc_y, pc_z, sih0_441, sig_224, \
                         sig_315, sih1_441, skf0_209, skf1_209, \
                         skg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_3 * pc_y[k] * skg_314[k];

        t_440[k] = f_13 * sig_224[k]
                   + f_1 * skf0_209[k]
                   - f_2 * skf1_209[k]
                   + f_3 * pc_z[k] * skg_314[k];

        t_441[k] = pb_x[k] * sih0_441[k]
                   + f_13 * sig_315[k]
                   - f_8 * pc_x[k] * sih1_441[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, pb_x, pc_x, pc_y, pc_z, sih0_444, \
                         sig_225, sig_227, sig_318, sih1_444, skg_315, \
                         skg_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = f_12 * sig_225[k]
                   + f_3 * pc_y[k] * skg_315[k];

        t_443[k] = f_3 * pc_z[k] * skg_315[k];

        t_444[k] = pb_x[k] * sih0_444[k]
                   + f_11 * sig_318[k]
                   - f_8 * pc_x[k] * sih1_444[k];

        t_445[k] = f_12 * sig_227[k]
                   + f_3 * pc_y[k] * skg_317[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pb_x, pc_x, pc_z, sih0_446, sih0_447, sig_320, \
                         sig_321, sih1_446, sih1_447, skg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = pb_x[k] * sih0_446[k]
                   + f_11 * sig_320[k]
                   - f_8 * pc_x[k] * sih1_446[k];

        t_447[k] = pb_x[k] * sih0_447[k]
                   + f_10 * sig_321[k]
                   - f_8 * pc_x[k] * sih1_447[k];

        t_448[k] = f_3 * pc_z[k] * skg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pb_x, pc_x, pc_y, sih0_450, sig_230, \
                         sig_324, sig_325, sig_326, sih1_450, skg_320, skg_325, \
                         skg_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_12 * sig_230[k]
                   + f_3 * pc_y[k] * skg_320[k];

        t_450[k] = pb_x[k] * sih0_450[k]
                   + f_10 * sig_324[k]
                   - f_8 * pc_x[k] * sih1_450[k];

        t_451[k] = f_9 * sig_325[k]
                   + f_3 * pc_x[k] * skg_325[k];

        t_452[k] = f_9 * sig_326[k]
                   + f_3 * pc_x[k] * skg_326[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pb_x, pc_x, sih0_456, sig_327, sig_328, \
                         sig_329, sih1_456, skg_327, skg_328, skg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_9 * sig_327[k]
                   + f_3 * pc_x[k] * skg_327[k];

        t_454[k] = f_9 * sig_328[k]
                   + f_3 * pc_x[k] * skg_328[k];

        t_455[k] = f_9 * sig_329[k]
                   + f_3 * pc_x[k] * skg_329[k];

        t_456[k] = pb_x[k] * sih0_456[k]
                   - f_8 * pc_x[k] * sih1_456[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pb_x, pc_x, pc_y, pc_z, sih0_458, \
                         sih0_459, sig_239, sih1_458, sih1_459, skg_325, \
                         skg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * skg_325[k];

        t_458[k] = pb_x[k] * sih0_458[k]
                   - f_8 * pc_x[k] * sih1_458[k];

        t_459[k] = pb_x[k] * sih0_459[k]
                   - f_8 * pc_x[k] * sih1_459[k];

        t_460[k] = f_12 * sig_239[k]
                   + f_3 * pc_y[k] * skg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pb_x, pb_z, pc_x, pc_y, pc_z, sih0_315, \
                         sih0_461, sig_225, sig_240, sih1_315, sih1_461, \
                         skg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = pb_x[k] * sih0_461[k]
                   - f_8 * pc_x[k] * sih1_461[k];

        t_462[k] = pb_z[k] * sih0_315[k]
                   - f_8 * pc_z[k] * sih1_315[k];

        t_463[k] = f_13 * sig_240[k]
                   + f_3 * pc_y[k] * skg_330[k];

        t_464[k] = f_9 * sig_225[k]
                   + f_3 * pc_z[k] * skg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pb_x, pb_z, pc_x, pc_y, pc_z, sih0_318, \
                         sih0_467, sig_242, sig_335, sih1_318, sih1_467, \
                         skg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pb_z[k] * sih0_318[k]
                   - f_8 * pc_z[k] * sih1_318[k];

        t_466[k] = f_13 * sig_242[k]
                   + f_3 * pc_y[k] * skg_332[k];

        t_467[k] = pb_x[k] * sih0_467[k]
                   + f_11 * sig_335[k]
                   - f_8 * pc_x[k] * sih1_467[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pb_z, pc_y, pc_z, sih0_321, sig_228, sig_245, \
                         sih1_321, skg_333, skg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pb_z[k] * sih0_321[k]
                   - f_8 * pc_z[k] * sih1_321[k];

        t_469[k] = f_9 * sig_228[k]
                   + f_3 * pc_z[k] * skg_333[k];

        t_470[k] = f_13 * sig_245[k]
                   + f_3 * pc_y[k] * skg_335[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pb_x, pc_x, sih0_471, sig_339, sig_340, \
                         sig_341, sig_342, sih1_471, skg_340, skg_341, \
                         skg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pb_x[k] * sih0_471[k]
                   + f_10 * sig_339[k]
                   - f_8 * pc_x[k] * sih1_471[k];

        t_472[k] = f_9 * sig_340[k]
                   + f_3 * pc_x[k] * skg_340[k];

        t_473[k] = f_9 * sig_341[k]
                   + f_3 * pc_x[k] * skg_341[k];

        t_474[k] = f_9 * sig_342[k]
                   + f_3 * pc_x[k] * skg_342[k];
    }
}

static auto
compute_prim_skh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sih0,
                                                          const size_t sig, const size_t sih1,
                                                          const size_t skf0, const size_t skf1,
                                                          const size_t skg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sih0_420 = buffer.data(sih0 + 420);
    const auto *sih0_425 = buffer.data(sih0 + 425);
    const auto *sih0_429 = buffer.data(sih0 + 429);
    const auto *sih0_477 = buffer.data(sih0 + 477);
    const auto *sih0_479 = buffer.data(sih0 + 479);
    const auto *sih0_480 = buffer.data(sih0 + 480);
    const auto *sih0_482 = buffer.data(sih0 + 482);
    const auto *sih0_483 = buffer.data(sih0 + 483);
    const auto *sih0_486 = buffer.data(sih0 + 486);
    const auto *sih0_488 = buffer.data(sih0 + 488);
    const auto *sih0_489 = buffer.data(sih0 + 489);
    const auto *sih0_492 = buffer.data(sih0 + 492);
    const auto *sih0_498 = buffer.data(sih0 + 498);
    const auto *sih0_500 = buffer.data(sih0 + 500);
    const auto *sih0_501 = buffer.data(sih0 + 501);
    const auto *sih0_503 = buffer.data(sih0 + 503);
    const auto *sih0_504 = buffer.data(sih0 + 504);
    const auto *sih0_507 = buffer.data(sih0 + 507);
    const auto *sih0_509 = buffer.data(sih0 + 509);
    const auto *sih0_510 = buffer.data(sih0 + 510);
    const auto *sih0_513 = buffer.data(sih0 + 513);
    const auto *sih0_519 = buffer.data(sih0 + 519);
    const auto *sih0_521 = buffer.data(sih0 + 521);
    const auto *sih0_522 = buffer.data(sih0 + 522);
    const auto *sih0_524 = buffer.data(sih0 + 524);
    const auto *sih0_525 = buffer.data(sih0 + 525);
    const auto *sih0_528 = buffer.data(sih0 + 528);
    const auto *sih0_530 = buffer.data(sih0 + 530);
    const auto *sih0_531 = buffer.data(sih0 + 531);
    const auto *sih0_534 = buffer.data(sih0 + 534);
    const auto *sih0_540 = buffer.data(sih0 + 540);
    const auto *sih0_542 = buffer.data(sih0 + 542);
    const auto *sih0_543 = buffer.data(sih0 + 543);
    const auto *sih0_545 = buffer.data(sih0 + 545);
    const auto *sih0_549 = buffer.data(sih0 + 549);
    const auto *sih0_552 = buffer.data(sih0 + 552);
    const auto *sih0_561 = buffer.data(sih0 + 561);
    const auto *sih0_563 = buffer.data(sih0 + 563);
    const auto *sih0_564 = buffer.data(sih0 + 564);
    const auto *sih0_566 = buffer.data(sih0 + 566);
    const auto *sih0_567 = buffer.data(sih0 + 567);
    const auto *sih0_570 = buffer.data(sih0 + 570);
    const auto *sih0_572 = buffer.data(sih0 + 572);
    const auto *sih0_573 = buffer.data(sih0 + 573);
    const auto *sih0_576 = buffer.data(sih0 + 576);
    const auto *sih0_582 = buffer.data(sih0 + 582);
    const auto *sih0_584 = buffer.data(sih0 + 584);
    const auto *sih0_585 = buffer.data(sih0 + 585);
    const auto *sih0_587 = buffer.data(sih0 + 587);

    const auto *sig_235 = buffer.data(sig + 235);
    const auto *sig_240 = buffer.data(sig + 240);
    const auto *sig_243 = buffer.data(sig + 243);
    const auto *sig_250 = buffer.data(sig + 250);
    const auto *sig_254 = buffer.data(sig + 254);
    const auto *sig_255 = buffer.data(sig + 255);
    const auto *sig_257 = buffer.data(sig + 257);
    const auto *sig_258 = buffer.data(sig + 258);
    const auto *sig_260 = buffer.data(sig + 260);
    const auto *sig_265 = buffer.data(sig + 265);
    const auto *sig_269 = buffer.data(sig + 269);
    const auto *sig_270 = buffer.data(sig + 270);
    const auto *sig_272 = buffer.data(sig + 272);
    const auto *sig_273 = buffer.data(sig + 273);
    const auto *sig_275 = buffer.data(sig + 275);
    const auto *sig_280 = buffer.data(sig + 280);
    const auto *sig_284 = buffer.data(sig + 284);
    const auto *sig_285 = buffer.data(sig + 285);
    const auto *sig_287 = buffer.data(sig + 287);
    const auto *sig_288 = buffer.data(sig + 288);
    const auto *sig_290 = buffer.data(sig + 290);
    const auto *sig_295 = buffer.data(sig + 295);
    const auto *sig_299 = buffer.data(sig + 299);
    const auto *sig_300 = buffer.data(sig + 300);
    const auto *sig_302 = buffer.data(sig + 302);
    const auto *sig_303 = buffer.data(sig + 303);
    const auto *sig_305 = buffer.data(sig + 305);
    const auto *sig_310 = buffer.data(sig + 310);
    const auto *sig_314 = buffer.data(sig + 314);
    const auto *sig_315 = buffer.data(sig + 315);
    const auto *sig_317 = buffer.data(sig + 317);
    const auto *sig_320 = buffer.data(sig + 320);
    const auto *sig_343 = buffer.data(sig + 343);
    const auto *sig_344 = buffer.data(sig + 344);
    const auto *sig_345 = buffer.data(sig + 345);
    const auto *sig_348 = buffer.data(sig + 348);
    const auto *sig_350 = buffer.data(sig + 350);
    const auto *sig_351 = buffer.data(sig + 351);
    const auto *sig_354 = buffer.data(sig + 354);
    const auto *sig_355 = buffer.data(sig + 355);
    const auto *sig_356 = buffer.data(sig + 356);
    const auto *sig_357 = buffer.data(sig + 357);
    const auto *sig_358 = buffer.data(sig + 358);
    const auto *sig_359 = buffer.data(sig + 359);
    const auto *sig_360 = buffer.data(sig + 360);
    const auto *sig_363 = buffer.data(sig + 363);
    const auto *sig_365 = buffer.data(sig + 365);
    const auto *sig_366 = buffer.data(sig + 366);
    const auto *sig_369 = buffer.data(sig + 369);
    const auto *sig_370 = buffer.data(sig + 370);
    const auto *sig_371 = buffer.data(sig + 371);
    const auto *sig_372 = buffer.data(sig + 372);
    const auto *sig_373 = buffer.data(sig + 373);
    const auto *sig_374 = buffer.data(sig + 374);
    const auto *sig_375 = buffer.data(sig + 375);
    const auto *sig_378 = buffer.data(sig + 378);
    const auto *sig_380 = buffer.data(sig + 380);
    const auto *sig_381 = buffer.data(sig + 381);
    const auto *sig_384 = buffer.data(sig + 384);
    const auto *sig_385 = buffer.data(sig + 385);
    const auto *sig_386 = buffer.data(sig + 386);
    const auto *sig_387 = buffer.data(sig + 387);
    const auto *sig_388 = buffer.data(sig + 388);
    const auto *sig_389 = buffer.data(sig + 389);
    const auto *sig_393 = buffer.data(sig + 393);
    const auto *sig_396 = buffer.data(sig + 396);
    const auto *sig_400 = buffer.data(sig + 400);
    const auto *sig_401 = buffer.data(sig + 401);
    const auto *sig_402 = buffer.data(sig + 402);
    const auto *sig_403 = buffer.data(sig + 403);
    const auto *sig_404 = buffer.data(sig + 404);
    const auto *sig_405 = buffer.data(sig + 405);
    const auto *sig_408 = buffer.data(sig + 408);
    const auto *sig_410 = buffer.data(sig + 410);
    const auto *sig_411 = buffer.data(sig + 411);
    const auto *sig_414 = buffer.data(sig + 414);
    const auto *sig_415 = buffer.data(sig + 415);
    const auto *sig_416 = buffer.data(sig + 416);
    const auto *sig_417 = buffer.data(sig + 417);
    const auto *sig_418 = buffer.data(sig + 418);
    const auto *sig_419 = buffer.data(sig + 419);

    const auto *sih1_420 = buffer.data(sih1 + 420);
    const auto *sih1_425 = buffer.data(sih1 + 425);
    const auto *sih1_429 = buffer.data(sih1 + 429);
    const auto *sih1_477 = buffer.data(sih1 + 477);
    const auto *sih1_479 = buffer.data(sih1 + 479);
    const auto *sih1_480 = buffer.data(sih1 + 480);
    const auto *sih1_482 = buffer.data(sih1 + 482);
    const auto *sih1_483 = buffer.data(sih1 + 483);
    const auto *sih1_486 = buffer.data(sih1 + 486);
    const auto *sih1_488 = buffer.data(sih1 + 488);
    const auto *sih1_489 = buffer.data(sih1 + 489);
    const auto *sih1_492 = buffer.data(sih1 + 492);
    const auto *sih1_498 = buffer.data(sih1 + 498);
    const auto *sih1_500 = buffer.data(sih1 + 500);
    const auto *sih1_501 = buffer.data(sih1 + 501);
    const auto *sih1_503 = buffer.data(sih1 + 503);
    const auto *sih1_504 = buffer.data(sih1 + 504);
    const auto *sih1_507 = buffer.data(sih1 + 507);
    const auto *sih1_509 = buffer.data(sih1 + 509);
    const auto *sih1_510 = buffer.data(sih1 + 510);
    const auto *sih1_513 = buffer.data(sih1 + 513);
    const auto *sih1_519 = buffer.data(sih1 + 519);
    const auto *sih1_521 = buffer.data(sih1 + 521);
    const auto *sih1_522 = buffer.data(sih1 + 522);
    const auto *sih1_524 = buffer.data(sih1 + 524);
    const auto *sih1_525 = buffer.data(sih1 + 525);
    const auto *sih1_528 = buffer.data(sih1 + 528);
    const auto *sih1_530 = buffer.data(sih1 + 530);
    const auto *sih1_531 = buffer.data(sih1 + 531);
    const auto *sih1_534 = buffer.data(sih1 + 534);
    const auto *sih1_540 = buffer.data(sih1 + 540);
    const auto *sih1_542 = buffer.data(sih1 + 542);
    const auto *sih1_543 = buffer.data(sih1 + 543);
    const auto *sih1_545 = buffer.data(sih1 + 545);
    const auto *sih1_549 = buffer.data(sih1 + 549);
    const auto *sih1_552 = buffer.data(sih1 + 552);
    const auto *sih1_561 = buffer.data(sih1 + 561);
    const auto *sih1_563 = buffer.data(sih1 + 563);
    const auto *sih1_564 = buffer.data(sih1 + 564);
    const auto *sih1_566 = buffer.data(sih1 + 566);
    const auto *sih1_567 = buffer.data(sih1 + 567);
    const auto *sih1_570 = buffer.data(sih1 + 570);
    const auto *sih1_572 = buffer.data(sih1 + 572);
    const auto *sih1_573 = buffer.data(sih1 + 573);
    const auto *sih1_576 = buffer.data(sih1 + 576);
    const auto *sih1_582 = buffer.data(sih1 + 582);
    const auto *sih1_584 = buffer.data(sih1 + 584);
    const auto *sih1_585 = buffer.data(sih1 + 585);
    const auto *sih1_587 = buffer.data(sih1 + 587);

    const auto *skf0_280 = buffer.data(skf0 + 280);
    const auto *skf0_283 = buffer.data(skf0 + 283);
    const auto *skf0_285 = buffer.data(skf0 + 285);
    const auto *skf0_286 = buffer.data(skf0 + 286);
    const auto *skf0_289 = buffer.data(skf0 + 289);

    const auto *skf1_280 = buffer.data(skf1 + 280);
    const auto *skf1_283 = buffer.data(skf1 + 283);
    const auto *skf1_285 = buffer.data(skf1 + 285);
    const auto *skf1_286 = buffer.data(skf1 + 286);
    const auto *skf1_289 = buffer.data(skf1 + 289);

    const auto *skg_340 = buffer.data(skg + 340);
    const auto *skg_343 = buffer.data(skg + 343);
    const auto *skg_344 = buffer.data(skg + 344);
    const auto *skg_345 = buffer.data(skg + 345);
    const auto *skg_347 = buffer.data(skg + 347);
    const auto *skg_348 = buffer.data(skg + 348);
    const auto *skg_350 = buffer.data(skg + 350);
    const auto *skg_355 = buffer.data(skg + 355);
    const auto *skg_356 = buffer.data(skg + 356);
    const auto *skg_357 = buffer.data(skg + 357);
    const auto *skg_358 = buffer.data(skg + 358);
    const auto *skg_359 = buffer.data(skg + 359);
    const auto *skg_360 = buffer.data(skg + 360);
    const auto *skg_362 = buffer.data(skg + 362);
    const auto *skg_363 = buffer.data(skg + 363);
    const auto *skg_365 = buffer.data(skg + 365);
    const auto *skg_370 = buffer.data(skg + 370);
    const auto *skg_371 = buffer.data(skg + 371);
    const auto *skg_372 = buffer.data(skg + 372);
    const auto *skg_373 = buffer.data(skg + 373);
    const auto *skg_374 = buffer.data(skg + 374);
    const auto *skg_375 = buffer.data(skg + 375);
    const auto *skg_377 = buffer.data(skg + 377);
    const auto *skg_378 = buffer.data(skg + 378);
    const auto *skg_380 = buffer.data(skg + 380);
    const auto *skg_385 = buffer.data(skg + 385);
    const auto *skg_386 = buffer.data(skg + 386);
    const auto *skg_387 = buffer.data(skg + 387);
    const auto *skg_388 = buffer.data(skg + 388);
    const auto *skg_389 = buffer.data(skg + 389);
    const auto *skg_390 = buffer.data(skg + 390);
    const auto *skg_392 = buffer.data(skg + 392);
    const auto *skg_393 = buffer.data(skg + 393);
    const auto *skg_395 = buffer.data(skg + 395);
    const auto *skg_400 = buffer.data(skg + 400);
    const auto *skg_401 = buffer.data(skg + 401);
    const auto *skg_402 = buffer.data(skg + 402);
    const auto *skg_403 = buffer.data(skg + 403);
    const auto *skg_404 = buffer.data(skg + 404);
    const auto *skg_405 = buffer.data(skg + 405);
    const auto *skg_407 = buffer.data(skg + 407);
    const auto *skg_408 = buffer.data(skg + 408);
    const auto *skg_410 = buffer.data(skg + 410);
    const auto *skg_415 = buffer.data(skg + 415);
    const auto *skg_416 = buffer.data(skg + 416);
    const auto *skg_417 = buffer.data(skg + 417);
    const auto *skg_418 = buffer.data(skg + 418);
    const auto *skg_419 = buffer.data(skg + 419);
    const auto *skg_420 = buffer.data(skg + 420);
    const auto *skg_422 = buffer.data(skg + 422);
    const auto *skg_423 = buffer.data(skg + 423);
    const auto *skg_425 = buffer.data(skg + 425);
    const auto *skg_426 = buffer.data(skg + 426);
    const auto *skg_429 = buffer.data(skg + 429);

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pb_x, pc_x, pc_z, sih0_477, sig_235, \
                         sig_343, sig_344, sih1_477, skg_340, skg_343, \
                         skg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_9 * sig_343[k]
                   + f_3 * pc_x[k] * skg_343[k];

        t_476[k] = f_9 * sig_344[k]
                   + f_3 * pc_x[k] * skg_344[k];

        t_477[k] = pb_x[k] * sih0_477[k]
                   - f_8 * pc_x[k] * sih1_477[k];

        t_478[k] = f_9 * sig_235[k]
                   + f_3 * pc_z[k] * skg_340[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pb_x, pc_x, pc_y, sih0_479, sih0_480, \
                         sih0_482, sig_254, sih1_479, sih1_480, sih1_482, \
                         skg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = pb_x[k] * sih0_479[k]
                   - f_8 * pc_x[k] * sih1_479[k];

        t_480[k] = pb_x[k] * sih0_480[k]
                   - f_8 * pc_x[k] * sih1_480[k];

        t_481[k] = f_13 * sig_254[k]
                   + f_3 * pc_y[k] * skg_344[k];

        t_482[k] = pb_x[k] * sih0_482[k]
                   - f_8 * pc_x[k] * sih1_482[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pb_x, pc_x, pc_y, pc_z, sih0_483, sig_240, \
                         sig_255, sig_345, sih1_483, skg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = pb_x[k] * sih0_483[k]
                   + f_13 * sig_345[k]
                   - f_8 * pc_x[k] * sih1_483[k];

        t_484[k] = f_14 * sig_255[k]
                   + f_3 * pc_y[k] * skg_345[k];

        t_485[k] = f_10 * sig_240[k]
                   + f_3 * pc_z[k] * skg_345[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pb_x, pc_x, pc_y, sih0_486, sih0_488, sig_257, \
                         sig_348, sig_350, sih1_486, sih1_488, \
                         skg_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = pb_x[k] * sih0_486[k]
                   + f_11 * sig_348[k]
                   - f_8 * pc_x[k] * sih1_486[k];

        t_487[k] = f_14 * sig_257[k]
                   + f_3 * pc_y[k] * skg_347[k];

        t_488[k] = pb_x[k] * sih0_488[k]
                   + f_11 * sig_350[k]
                   - f_8 * pc_x[k] * sih1_488[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pb_x, pc_x, pc_y, pc_z, sih0_489, sig_243, \
                         sig_260, sig_351, sih1_489, skg_348, skg_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = pb_x[k] * sih0_489[k]
                   + f_10 * sig_351[k]
                   - f_8 * pc_x[k] * sih1_489[k];

        t_490[k] = f_10 * sig_243[k]
                   + f_3 * pc_z[k] * skg_348[k];

        t_491[k] = f_14 * sig_260[k]
                   + f_3 * pc_y[k] * skg_350[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pb_x, pc_x, sih0_492, sig_354, sig_355, \
                         sig_356, sig_357, sih1_492, skg_355, skg_356, \
                         skg_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = pb_x[k] * sih0_492[k]
                   + f_10 * sig_354[k]
                   - f_8 * pc_x[k] * sih1_492[k];

        t_493[k] = f_9 * sig_355[k]
                   + f_3 * pc_x[k] * skg_355[k];

        t_494[k] = f_9 * sig_356[k]
                   + f_3 * pc_x[k] * skg_356[k];

        t_495[k] = f_9 * sig_357[k]
                   + f_3 * pc_x[k] * skg_357[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_x, pc_x, pc_z, sih0_498, sig_250, \
                         sig_358, sig_359, sih1_498, skg_355, skg_358, \
                         skg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_9 * sig_358[k]
                   + f_3 * pc_x[k] * skg_358[k];

        t_497[k] = f_9 * sig_359[k]
                   + f_3 * pc_x[k] * skg_359[k];

        t_498[k] = pb_x[k] * sih0_498[k]
                   - f_8 * pc_x[k] * sih1_498[k];

        t_499[k] = f_10 * sig_250[k]
                   + f_3 * pc_z[k] * skg_355[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pb_x, pc_x, pc_y, sih0_500, sih0_501, \
                         sih0_503, sig_269, sih1_500, sih1_501, sih1_503, \
                         skg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pb_x[k] * sih0_500[k]
                   - f_8 * pc_x[k] * sih1_500[k];

        t_501[k] = pb_x[k] * sih0_501[k]
                   - f_8 * pc_x[k] * sih1_501[k];

        t_502[k] = f_14 * sig_269[k]
                   + f_3 * pc_y[k] * skg_359[k];

        t_503[k] = pb_x[k] * sih0_503[k]
                   - f_8 * pc_x[k] * sih1_503[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pb_x, pc_x, pc_y, pc_z, sih0_504, sig_255, \
                         sig_270, sig_360, sih1_504, skg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = pb_x[k] * sih0_504[k]
                   + f_13 * sig_360[k]
                   - f_8 * pc_x[k] * sih1_504[k];

        t_505[k] = f_11 * sig_270[k]
                   + f_3 * pc_y[k] * skg_360[k];

        t_506[k] = f_11 * sig_255[k]
                   + f_3 * pc_z[k] * skg_360[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pb_x, pc_x, pc_y, sih0_507, sih0_509, sig_272, \
                         sig_363, sig_365, sih1_507, sih1_509, \
                         skg_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = pb_x[k] * sih0_507[k]
                   + f_11 * sig_363[k]
                   - f_8 * pc_x[k] * sih1_507[k];

        t_508[k] = f_11 * sig_272[k]
                   + f_3 * pc_y[k] * skg_362[k];

        t_509[k] = pb_x[k] * sih0_509[k]
                   + f_11 * sig_365[k]
                   - f_8 * pc_x[k] * sih1_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, pb_x, pc_x, pc_y, pc_z, sih0_510, sig_258, \
                         sig_275, sig_366, sih1_510, skg_363, skg_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = pb_x[k] * sih0_510[k]
                   + f_10 * sig_366[k]
                   - f_8 * pc_x[k] * sih1_510[k];

        t_511[k] = f_11 * sig_258[k]
                   + f_3 * pc_z[k] * skg_363[k];

        t_512[k] = f_11 * sig_275[k]
                   + f_3 * pc_y[k] * skg_365[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, t_516, pb_x, pc_x, sih0_513, sig_369, sig_370, \
                         sig_371, sig_372, sih1_513, skg_370, skg_371, \
                         skg_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = pb_x[k] * sih0_513[k]
                   + f_10 * sig_369[k]
                   - f_8 * pc_x[k] * sih1_513[k];

        t_514[k] = f_9 * sig_370[k]
                   + f_3 * pc_x[k] * skg_370[k];

        t_515[k] = f_9 * sig_371[k]
                   + f_3 * pc_x[k] * skg_371[k];

        t_516[k] = f_9 * sig_372[k]
                   + f_3 * pc_x[k] * skg_372[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pb_x, pc_x, pc_z, sih0_519, sig_265, \
                         sig_373, sig_374, sih1_519, skg_370, skg_373, \
                         skg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_9 * sig_373[k]
                   + f_3 * pc_x[k] * skg_373[k];

        t_518[k] = f_9 * sig_374[k]
                   + f_3 * pc_x[k] * skg_374[k];

        t_519[k] = pb_x[k] * sih0_519[k]
                   - f_8 * pc_x[k] * sih1_519[k];

        t_520[k] = f_11 * sig_265[k]
                   + f_3 * pc_z[k] * skg_370[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pb_x, pc_x, pc_y, sih0_521, sih0_522, \
                         sih0_524, sig_284, sih1_521, sih1_522, sih1_524, \
                         skg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = pb_x[k] * sih0_521[k]
                   - f_8 * pc_x[k] * sih1_521[k];

        t_522[k] = pb_x[k] * sih0_522[k]
                   - f_8 * pc_x[k] * sih1_522[k];

        t_523[k] = f_11 * sig_284[k]
                   + f_3 * pc_y[k] * skg_374[k];

        t_524[k] = pb_x[k] * sih0_524[k]
                   - f_8 * pc_x[k] * sih1_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, pb_x, pc_x, pc_y, pc_z, sih0_525, sig_270, \
                         sig_285, sig_375, sih1_525, skg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = pb_x[k] * sih0_525[k]
                   + f_13 * sig_375[k]
                   - f_8 * pc_x[k] * sih1_525[k];

        t_526[k] = f_10 * sig_285[k]
                   + f_3 * pc_y[k] * skg_375[k];

        t_527[k] = f_14 * sig_270[k]
                   + f_3 * pc_z[k] * skg_375[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, pb_x, pc_x, pc_y, sih0_528, sih0_530, sig_287, \
                         sig_378, sig_380, sih1_528, sih1_530, \
                         skg_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = pb_x[k] * sih0_528[k]
                   + f_11 * sig_378[k]
                   - f_8 * pc_x[k] * sih1_528[k];

        t_529[k] = f_10 * sig_287[k]
                   + f_3 * pc_y[k] * skg_377[k];

        t_530[k] = pb_x[k] * sih0_530[k]
                   + f_11 * sig_380[k]
                   - f_8 * pc_x[k] * sih1_530[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pb_x, pc_x, pc_y, pc_z, sih0_531, sig_273, \
                         sig_290, sig_381, sih1_531, skg_378, skg_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = pb_x[k] * sih0_531[k]
                   + f_10 * sig_381[k]
                   - f_8 * pc_x[k] * sih1_531[k];

        t_532[k] = f_14 * sig_273[k]
                   + f_3 * pc_z[k] * skg_378[k];

        t_533[k] = f_10 * sig_290[k]
                   + f_3 * pc_y[k] * skg_380[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pb_x, pc_x, sih0_534, sig_384, sig_385, \
                         sig_386, sig_387, sih1_534, skg_385, skg_386, \
                         skg_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = pb_x[k] * sih0_534[k]
                   + f_10 * sig_384[k]
                   - f_8 * pc_x[k] * sih1_534[k];

        t_535[k] = f_9 * sig_385[k]
                   + f_3 * pc_x[k] * skg_385[k];

        t_536[k] = f_9 * sig_386[k]
                   + f_3 * pc_x[k] * skg_386[k];

        t_537[k] = f_9 * sig_387[k]
                   + f_3 * pc_x[k] * skg_387[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pb_x, pc_x, pc_z, sih0_540, sig_280, \
                         sig_388, sig_389, sih1_540, skg_385, skg_388, \
                         skg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_9 * sig_388[k]
                   + f_3 * pc_x[k] * skg_388[k];

        t_539[k] = f_9 * sig_389[k]
                   + f_3 * pc_x[k] * skg_389[k];

        t_540[k] = pb_x[k] * sih0_540[k]
                   - f_8 * pc_x[k] * sih1_540[k];

        t_541[k] = f_14 * sig_280[k]
                   + f_3 * pc_z[k] * skg_385[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, pb_x, pc_x, pc_y, sih0_542, sih0_543, \
                         sih0_545, sig_299, sih1_542, sih1_543, sih1_545, \
                         skg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = pb_x[k] * sih0_542[k]
                   - f_8 * pc_x[k] * sih1_542[k];

        t_543[k] = pb_x[k] * sih0_543[k]
                   - f_8 * pc_x[k] * sih1_543[k];

        t_544[k] = f_10 * sig_299[k]
                   + f_3 * pc_y[k] * skg_389[k];

        t_545[k] = pb_x[k] * sih0_545[k]
                   - f_8 * pc_x[k] * sih1_545[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, pb_y, pc_y, pc_z, sih0_420, sig_285, sig_300, \
                         sih1_420, skg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pb_y[k] * sih0_420[k]
                   - f_8 * pc_y[k] * sih1_420[k];

        t_547[k] = f_9 * sig_300[k]
                   + f_3 * pc_y[k] * skg_390[k];

        t_548[k] = f_13 * sig_285[k]
                   + f_3 * pc_z[k] * skg_390[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pb_x, pb_y, pc_x, pc_y, sih0_425, sih0_549, \
                         sig_302, sig_393, sih1_425, sih1_549, \
                         skg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = pb_x[k] * sih0_549[k]
                   + f_11 * sig_393[k]
                   - f_8 * pc_x[k] * sih1_549[k];

        t_550[k] = f_9 * sig_302[k]
                   + f_3 * pc_y[k] * skg_392[k];

        t_551[k] = pb_y[k] * sih0_425[k]
                   - f_8 * pc_y[k] * sih1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pb_x, pc_x, pc_y, pc_z, sih0_552, sig_288, \
                         sig_305, sig_396, sih1_552, skg_393, skg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pb_x[k] * sih0_552[k]
                   + f_10 * sig_396[k]
                   - f_8 * pc_x[k] * sih1_552[k];

        t_553[k] = f_13 * sig_288[k]
                   + f_3 * pc_z[k] * skg_393[k];

        t_554[k] = f_9 * sig_305[k]
                   + f_3 * pc_y[k] * skg_395[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, pb_y, pc_x, pc_y, sih0_429, sig_400, \
                         sig_401, sig_402, sih1_429, skg_400, skg_401, \
                         skg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pb_y[k] * sih0_429[k]
                   - f_8 * pc_y[k] * sih1_429[k];

        t_556[k] = f_9 * sig_400[k]
                   + f_3 * pc_x[k] * skg_400[k];

        t_557[k] = f_9 * sig_401[k]
                   + f_3 * pc_x[k] * skg_401[k];

        t_558[k] = f_9 * sig_402[k]
                   + f_3 * pc_x[k] * skg_402[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_x, pc_x, pc_z, sih0_561, sig_295, \
                         sig_403, sig_404, sih1_561, skg_400, skg_403, \
                         skg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_9 * sig_403[k]
                   + f_3 * pc_x[k] * skg_403[k];

        t_560[k] = f_9 * sig_404[k]
                   + f_3 * pc_x[k] * skg_404[k];

        t_561[k] = pb_x[k] * sih0_561[k]
                   - f_8 * pc_x[k] * sih1_561[k];

        t_562[k] = f_13 * sig_295[k]
                   + f_3 * pc_z[k] * skg_400[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, pb_x, pc_x, pc_y, sih0_563, sih0_564, \
                         sih0_566, sig_314, sih1_563, sih1_564, sih1_566, \
                         skg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = pb_x[k] * sih0_563[k]
                   - f_8 * pc_x[k] * sih1_563[k];

        t_564[k] = pb_x[k] * sih0_564[k]
                   - f_8 * pc_x[k] * sih1_564[k];

        t_565[k] = f_9 * sig_314[k]
                   + f_3 * pc_y[k] * skg_404[k];

        t_566[k] = pb_x[k] * sih0_566[k]
                   - f_8 * pc_x[k] * sih1_566[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pb_x, pc_x, pc_y, pc_z, sih0_567, \
                         sih0_570, sig_300, sig_405, sig_408, sih1_567, sih1_570, \
                         skg_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = pb_x[k] * sih0_567[k]
                   + f_13 * sig_405[k]
                   - f_8 * pc_x[k] * sih1_567[k];

        t_568[k] = f_3 * pc_y[k] * skg_405[k];

        t_569[k] = f_12 * sig_300[k]
                   + f_3 * pc_z[k] * skg_405[k];

        t_570[k] = pb_x[k] * sih0_570[k]
                   + f_11 * sig_408[k]
                   - f_8 * pc_x[k] * sih1_570[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pb_x, pc_x, pc_y, sih0_572, sih0_573, sig_410, \
                         sig_411, sih1_572, sih1_573, skg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_3 * pc_y[k] * skg_407[k];

        t_572[k] = pb_x[k] * sih0_572[k]
                   + f_11 * sig_410[k]
                   - f_8 * pc_x[k] * sih1_572[k];

        t_573[k] = pb_x[k] * sih0_573[k]
                   + f_10 * sig_411[k]
                   - f_8 * pc_x[k] * sih1_573[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pb_x, pc_x, pc_y, pc_z, sih0_576, \
                         sig_303, sig_414, sig_415, sih1_576, skg_408, skg_410, \
                         skg_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_12 * sig_303[k]
                   + f_3 * pc_z[k] * skg_408[k];

        t_575[k] = f_3 * pc_y[k] * skg_410[k];

        t_576[k] = pb_x[k] * sih0_576[k]
                   + f_10 * sig_414[k]
                   - f_8 * pc_x[k] * sih1_576[k];

        t_577[k] = f_9 * sig_415[k]
                   + f_3 * pc_x[k] * skg_415[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pc_x, sig_416, sig_417, sig_418, sig_419, \
                         skg_416, skg_417, skg_418, skg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_9 * sig_416[k]
                   + f_3 * pc_x[k] * skg_416[k];

        t_579[k] = f_9 * sig_417[k]
                   + f_3 * pc_x[k] * skg_417[k];

        t_580[k] = f_9 * sig_418[k]
                   + f_3 * pc_x[k] * skg_418[k];

        t_581[k] = f_9 * sig_419[k]
                   + f_3 * pc_x[k] * skg_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pb_x, pc_x, pc_z, sih0_582, sih0_584, \
                         sih0_585, sig_310, sih1_582, sih1_584, sih1_585, \
                         skg_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = pb_x[k] * sih0_582[k]
                   - f_8 * pc_x[k] * sih1_582[k];

        t_583[k] = f_12 * sig_310[k]
                   + f_3 * pc_z[k] * skg_415[k];

        t_584[k] = pb_x[k] * sih0_584[k]
                   - f_8 * pc_x[k] * sih1_584[k];

        t_585[k] = pb_x[k] * sih0_585[k]
                   - f_8 * pc_x[k] * sih1_585[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, pb_x, pc_x, pc_y, pc_z, sih0_587, \
                         sig_315, sih1_587, skf0_280, skf1_280, skg_419, \
                         skg_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * skg_419[k];

        t_587[k] = pb_x[k] * sih0_587[k]
                   - f_8 * pc_x[k] * sih1_587[k];

        t_588[k] = f_1 * skf0_280[k]
                   - f_2 * skf1_280[k]
                   + f_3 * pc_x[k] * skg_420[k];

        t_589[k] = f_0 * sig_315[k]
                   + f_3 * pc_y[k] * skg_420[k];

        t_590[k] = f_3 * pc_z[k] * skg_420[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, pc_x, pc_y, sig_317, skf0_283, skf0_285, \
                         skf1_283, skf1_285, skg_422, skg_423, \
                         skg_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_4 * skf0_283[k]
                   - f_5 * skf1_283[k]
                   + f_3 * pc_x[k] * skg_423[k];

        t_592[k] = f_0 * sig_317[k]
                   + f_3 * pc_y[k] * skg_422[k];

        t_593[k] = f_4 * skf0_285[k]
                   - f_5 * skf1_285[k]
                   + f_3 * pc_x[k] * skg_425[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pc_x, pc_y, pc_z, sig_320, skf0_286, \
                         skf0_289, skf1_286, skf1_289, skg_423, skg_425, skg_426, \
                         skg_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_6 * skf0_286[k]
                   - f_7 * skf1_286[k]
                   + f_3 * pc_x[k] * skg_426[k];

        t_595[k] = f_3 * pc_z[k] * skg_423[k];

        t_596[k] = f_0 * sig_320[k]
                   + f_3 * pc_y[k] * skg_425[k];

        t_597[k] = f_6 * skf0_289[k]
                   - f_7 * skf1_289[k]
                   + f_3 * pc_x[k] * skg_429[k];
    }
}

static auto
compute_prim_skh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sih0,
                                                          const size_t sig, const size_t sih1,
                                                          const size_t skf0, const size_t skf1,
                                                          const size_t skg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sih0_441 = buffer.data(sih0 + 441);
    const auto *sih0_444 = buffer.data(sih0 + 444);
    const auto *sih0_447 = buffer.data(sih0 + 447);
    const auto *sih0_456 = buffer.data(sih0 + 456);
    const auto *sih0_458 = buffer.data(sih0 + 458);
    const auto *sih0_459 = buffer.data(sih0 + 459);
    const auto *sih0_567 = buffer.data(sih0 + 567);
    const auto *sih0_572 = buffer.data(sih0 + 572);
    const auto *sih0_576 = buffer.data(sih0 + 576);

    const auto *sig_315 = buffer.data(sig + 315);
    const auto *sig_318 = buffer.data(sig + 318);
    const auto *sig_325 = buffer.data(sig + 325);
    const auto *sig_326 = buffer.data(sig + 326);
    const auto *sig_327 = buffer.data(sig + 327);
    const auto *sig_328 = buffer.data(sig + 328);
    const auto *sig_329 = buffer.data(sig + 329);
    const auto *sig_330 = buffer.data(sig + 330);
    const auto *sig_332 = buffer.data(sig + 332);
    const auto *sig_333 = buffer.data(sig + 333);
    const auto *sig_335 = buffer.data(sig + 335);
    const auto *sig_340 = buffer.data(sig + 340);
    const auto *sig_344 = buffer.data(sig + 344);
    const auto *sig_345 = buffer.data(sig + 345);
    const auto *sig_347 = buffer.data(sig + 347);
    const auto *sig_348 = buffer.data(sig + 348);
    const auto *sig_350 = buffer.data(sig + 350);
    const auto *sig_355 = buffer.data(sig + 355);
    const auto *sig_357 = buffer.data(sig + 357);
    const auto *sig_358 = buffer.data(sig + 358);
    const auto *sig_359 = buffer.data(sig + 359);
    const auto *sig_360 = buffer.data(sig + 360);
    const auto *sig_362 = buffer.data(sig + 362);
    const auto *sig_363 = buffer.data(sig + 363);
    const auto *sig_365 = buffer.data(sig + 365);
    const auto *sig_370 = buffer.data(sig + 370);
    const auto *sig_372 = buffer.data(sig + 372);
    const auto *sig_373 = buffer.data(sig + 373);
    const auto *sig_374 = buffer.data(sig + 374);
    const auto *sig_375 = buffer.data(sig + 375);
    const auto *sig_377 = buffer.data(sig + 377);
    const auto *sig_378 = buffer.data(sig + 378);
    const auto *sig_380 = buffer.data(sig + 380);
    const auto *sig_385 = buffer.data(sig + 385);
    const auto *sig_387 = buffer.data(sig + 387);
    const auto *sig_388 = buffer.data(sig + 388);
    const auto *sig_389 = buffer.data(sig + 389);
    const auto *sig_390 = buffer.data(sig + 390);
    const auto *sig_392 = buffer.data(sig + 392);
    const auto *sig_393 = buffer.data(sig + 393);
    const auto *sig_395 = buffer.data(sig + 395);
    const auto *sig_400 = buffer.data(sig + 400);
    const auto *sig_402 = buffer.data(sig + 402);
    const auto *sig_403 = buffer.data(sig + 403);
    const auto *sig_404 = buffer.data(sig + 404);
    const auto *sig_405 = buffer.data(sig + 405);
    const auto *sig_407 = buffer.data(sig + 407);
    const auto *sig_410 = buffer.data(sig + 410);

    const auto *sih1_441 = buffer.data(sih1 + 441);
    const auto *sih1_444 = buffer.data(sih1 + 444);
    const auto *sih1_447 = buffer.data(sih1 + 447);
    const auto *sih1_456 = buffer.data(sih1 + 456);
    const auto *sih1_458 = buffer.data(sih1 + 458);
    const auto *sih1_459 = buffer.data(sih1 + 459);
    const auto *sih1_567 = buffer.data(sih1 + 567);
    const auto *sih1_572 = buffer.data(sih1 + 572);
    const auto *sih1_576 = buffer.data(sih1 + 576);

    const auto *skf0_286 = buffer.data(skf0 + 286);
    const auto *skf0_288 = buffer.data(skf0 + 288);
    const auto *skf0_289 = buffer.data(skf0 + 289);
    const auto *skf0_295 = buffer.data(skf0 + 295);
    const auto *skf0_299 = buffer.data(skf0 + 299);
    const auto *skf0_300 = buffer.data(skf0 + 300);
    const auto *skf0_303 = buffer.data(skf0 + 303);
    const auto *skf0_305 = buffer.data(skf0 + 305);
    const auto *skf0_306 = buffer.data(skf0 + 306);
    const auto *skf0_308 = buffer.data(skf0 + 308);
    const auto *skf0_309 = buffer.data(skf0 + 309);
    const auto *skf0_310 = buffer.data(skf0 + 310);
    const auto *skf0_313 = buffer.data(skf0 + 313);
    const auto *skf0_315 = buffer.data(skf0 + 315);
    const auto *skf0_316 = buffer.data(skf0 + 316);
    const auto *skf0_318 = buffer.data(skf0 + 318);
    const auto *skf0_319 = buffer.data(skf0 + 319);
    const auto *skf0_320 = buffer.data(skf0 + 320);
    const auto *skf0_323 = buffer.data(skf0 + 323);
    const auto *skf0_325 = buffer.data(skf0 + 325);
    const auto *skf0_326 = buffer.data(skf0 + 326);
    const auto *skf0_328 = buffer.data(skf0 + 328);
    const auto *skf0_329 = buffer.data(skf0 + 329);
    const auto *skf0_330 = buffer.data(skf0 + 330);
    const auto *skf0_333 = buffer.data(skf0 + 333);
    const auto *skf0_335 = buffer.data(skf0 + 335);
    const auto *skf0_336 = buffer.data(skf0 + 336);
    const auto *skf0_338 = buffer.data(skf0 + 338);
    const auto *skf0_339 = buffer.data(skf0 + 339);
    const auto *skf0_343 = buffer.data(skf0 + 343);
    const auto *skf0_346 = buffer.data(skf0 + 346);

    const auto *skf1_286 = buffer.data(skf1 + 286);
    const auto *skf1_288 = buffer.data(skf1 + 288);
    const auto *skf1_289 = buffer.data(skf1 + 289);
    const auto *skf1_295 = buffer.data(skf1 + 295);
    const auto *skf1_299 = buffer.data(skf1 + 299);
    const auto *skf1_300 = buffer.data(skf1 + 300);
    const auto *skf1_303 = buffer.data(skf1 + 303);
    const auto *skf1_305 = buffer.data(skf1 + 305);
    const auto *skf1_306 = buffer.data(skf1 + 306);
    const auto *skf1_308 = buffer.data(skf1 + 308);
    const auto *skf1_309 = buffer.data(skf1 + 309);
    const auto *skf1_310 = buffer.data(skf1 + 310);
    const auto *skf1_313 = buffer.data(skf1 + 313);
    const auto *skf1_315 = buffer.data(skf1 + 315);
    const auto *skf1_316 = buffer.data(skf1 + 316);
    const auto *skf1_318 = buffer.data(skf1 + 318);
    const auto *skf1_319 = buffer.data(skf1 + 319);
    const auto *skf1_320 = buffer.data(skf1 + 320);
    const auto *skf1_323 = buffer.data(skf1 + 323);
    const auto *skf1_325 = buffer.data(skf1 + 325);
    const auto *skf1_326 = buffer.data(skf1 + 326);
    const auto *skf1_328 = buffer.data(skf1 + 328);
    const auto *skf1_329 = buffer.data(skf1 + 329);
    const auto *skf1_330 = buffer.data(skf1 + 330);
    const auto *skf1_333 = buffer.data(skf1 + 333);
    const auto *skf1_335 = buffer.data(skf1 + 335);
    const auto *skf1_336 = buffer.data(skf1 + 336);
    const auto *skf1_338 = buffer.data(skf1 + 338);
    const auto *skf1_339 = buffer.data(skf1 + 339);
    const auto *skf1_343 = buffer.data(skf1 + 343);
    const auto *skf1_346 = buffer.data(skf1 + 346);

    const auto *skg_430 = buffer.data(skg + 430);
    const auto *skg_431 = buffer.data(skg + 431);
    const auto *skg_432 = buffer.data(skg + 432);
    const auto *skg_433 = buffer.data(skg + 433);
    const auto *skg_434 = buffer.data(skg + 434);
    const auto *skg_435 = buffer.data(skg + 435);
    const auto *skg_437 = buffer.data(skg + 437);
    const auto *skg_438 = buffer.data(skg + 438);
    const auto *skg_440 = buffer.data(skg + 440);
    const auto *skg_444 = buffer.data(skg + 444);
    const auto *skg_445 = buffer.data(skg + 445);
    const auto *skg_446 = buffer.data(skg + 446);
    const auto *skg_447 = buffer.data(skg + 447);
    const auto *skg_448 = buffer.data(skg + 448);
    const auto *skg_449 = buffer.data(skg + 449);
    const auto *skg_450 = buffer.data(skg + 450);
    const auto *skg_452 = buffer.data(skg + 452);
    const auto *skg_453 = buffer.data(skg + 453);
    const auto *skg_455 = buffer.data(skg + 455);
    const auto *skg_456 = buffer.data(skg + 456);
    const auto *skg_459 = buffer.data(skg + 459);
    const auto *skg_460 = buffer.data(skg + 460);
    const auto *skg_461 = buffer.data(skg + 461);
    const auto *skg_462 = buffer.data(skg + 462);
    const auto *skg_463 = buffer.data(skg + 463);
    const auto *skg_464 = buffer.data(skg + 464);
    const auto *skg_465 = buffer.data(skg + 465);
    const auto *skg_467 = buffer.data(skg + 467);
    const auto *skg_468 = buffer.data(skg + 468);
    const auto *skg_470 = buffer.data(skg + 470);
    const auto *skg_471 = buffer.data(skg + 471);
    const auto *skg_474 = buffer.data(skg + 474);
    const auto *skg_475 = buffer.data(skg + 475);
    const auto *skg_476 = buffer.data(skg + 476);
    const auto *skg_477 = buffer.data(skg + 477);
    const auto *skg_478 = buffer.data(skg + 478);
    const auto *skg_479 = buffer.data(skg + 479);
    const auto *skg_480 = buffer.data(skg + 480);
    const auto *skg_482 = buffer.data(skg + 482);
    const auto *skg_483 = buffer.data(skg + 483);
    const auto *skg_485 = buffer.data(skg + 485);
    const auto *skg_486 = buffer.data(skg + 486);
    const auto *skg_489 = buffer.data(skg + 489);
    const auto *skg_490 = buffer.data(skg + 490);
    const auto *skg_491 = buffer.data(skg + 491);
    const auto *skg_492 = buffer.data(skg + 492);
    const auto *skg_493 = buffer.data(skg + 493);
    const auto *skg_494 = buffer.data(skg + 494);
    const auto *skg_495 = buffer.data(skg + 495);
    const auto *skg_497 = buffer.data(skg + 497);
    const auto *skg_498 = buffer.data(skg + 498);
    const auto *skg_500 = buffer.data(skg + 500);
    const auto *skg_501 = buffer.data(skg + 501);
    const auto *skg_504 = buffer.data(skg + 504);
    const auto *skg_505 = buffer.data(skg + 505);
    const auto *skg_506 = buffer.data(skg + 506);
    const auto *skg_507 = buffer.data(skg + 507);
    const auto *skg_508 = buffer.data(skg + 508);
    const auto *skg_509 = buffer.data(skg + 509);
    const auto *skg_510 = buffer.data(skg + 510);
    const auto *skg_512 = buffer.data(skg + 512);
    const auto *skg_513 = buffer.data(skg + 513);
    const auto *skg_515 = buffer.data(skg + 515);
    const auto *skg_516 = buffer.data(skg + 516);
    const auto *skg_520 = buffer.data(skg + 520);
    const auto *skg_521 = buffer.data(skg + 521);
    const auto *skg_522 = buffer.data(skg + 522);

#pragma omp simd aligned(t_598, t_599, t_600, t_601, t_602, t_603, pc_x, pc_y, sig_325, \
                         skf0_286, skf1_286, skg_430, skg_431, skg_432, skg_433, \
                         skg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_3 * pc_x[k] * skg_430[k];

        t_599[k] = f_3 * pc_x[k] * skg_431[k];

        t_600[k] = f_3 * pc_x[k] * skg_432[k];

        t_601[k] = f_3 * pc_x[k] * skg_433[k];

        t_602[k] = f_3 * pc_x[k] * skg_434[k];

        t_603[k] = f_0 * sig_325[k]
                   + f_1 * skf0_286[k]
                   - f_2 * skf1_286[k]
                   + f_3 * pc_y[k] * skg_430[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pc_y, pc_z, sig_327, sig_328, skf0_288, \
                         skf0_289, skf1_288, skf1_289, skg_430, skg_432, \
                         skg_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * skg_430[k];

        t_605[k] = f_0 * sig_327[k]
                   + f_4 * skf0_288[k]
                   - f_5 * skf1_288[k]
                   + f_3 * pc_y[k] * skg_432[k];

        t_606[k] = f_0 * sig_328[k]
                   + f_6 * skf0_289[k]
                   - f_7 * skf1_289[k]
                   + f_3 * pc_y[k] * skg_433[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, pb_z, pc_y, pc_z, sih0_441, sig_329, \
                         sig_330, sih1_441, skf0_289, skf1_289, skg_434, \
                         skg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_0 * sig_329[k]
                   + f_3 * pc_y[k] * skg_434[k];

        t_608[k] = f_1 * skf0_289[k]
                   - f_2 * skf1_289[k]
                   + f_3 * pc_z[k] * skg_434[k];

        t_609[k] = pb_z[k] * sih0_441[k]
                   - f_8 * pc_z[k] * sih1_441[k];

        t_610[k] = f_12 * sig_330[k]
                   + f_3 * pc_y[k] * skg_435[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, pb_z, pc_y, pc_z, sih0_444, sig_315, sig_332, \
                         sih1_444, skg_435, skg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_9 * sig_315[k]
                   + f_3 * pc_z[k] * skg_435[k];

        t_612[k] = pb_z[k] * sih0_444[k]
                   - f_8 * pc_z[k] * sih1_444[k];

        t_613[k] = f_12 * sig_332[k]
                   + f_3 * pc_y[k] * skg_437[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, pb_z, pc_x, pc_y, pc_z, sih0_447, \
                         sig_318, sig_335, sih1_447, skf0_295, skf1_295, skg_438, \
                         skg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_4 * skf0_295[k]
                   - f_5 * skf1_295[k]
                   + f_3 * pc_x[k] * skg_440[k];

        t_615[k] = pb_z[k] * sih0_447[k]
                   - f_8 * pc_z[k] * sih1_447[k];

        t_616[k] = f_9 * sig_318[k]
                   + f_3 * pc_z[k] * skg_438[k];

        t_617[k] = f_12 * sig_335[k]
                   + f_3 * pc_y[k] * skg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, t_622, t_623, pc_x, skf0_299, skf1_299, \
                         skg_444, skg_445, skg_446, skg_447, skg_448, \
                         skg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_6 * skf0_299[k]
                   - f_7 * skf1_299[k]
                   + f_3 * pc_x[k] * skg_444[k];

        t_619[k] = f_3 * pc_x[k] * skg_445[k];

        t_620[k] = f_3 * pc_x[k] * skg_446[k];

        t_621[k] = f_3 * pc_x[k] * skg_447[k];

        t_622[k] = f_3 * pc_x[k] * skg_448[k];

        t_623[k] = f_3 * pc_x[k] * skg_449[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, pb_z, pc_z, sih0_456, sih0_458, sih0_459, \
                         sig_325, sig_326, sig_327, sih1_456, sih1_458, sih1_459, \
                         skg_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pb_z[k] * sih0_456[k]
                   - f_8 * pc_z[k] * sih1_456[k];

        t_625[k] = f_9 * sig_325[k]
                   + f_3 * pc_z[k] * skg_445[k];

        t_626[k] = pb_z[k] * sih0_458[k]
                   + f_10 * sig_326[k]
                   - f_8 * pc_z[k] * sih1_458[k];

        t_627[k] = pb_z[k] * sih0_459[k]
                   + f_11 * sig_327[k]
                   - f_8 * pc_z[k] * sih1_459[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, t_631, pc_x, pc_y, pc_z, sig_329, sig_344, \
                         sig_345, skf0_299, skf0_300, skf1_299, skf1_300, skg_449, \
                         skg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_12 * sig_344[k]
                   + f_3 * pc_y[k] * skg_449[k];

        t_629[k] = f_9 * sig_329[k]
                   + f_1 * skf0_299[k]
                   - f_2 * skf1_299[k]
                   + f_3 * pc_z[k] * skg_449[k];

        t_630[k] = f_1 * skf0_300[k]
                   - f_2 * skf1_300[k]
                   + f_3 * pc_x[k] * skg_450[k];

        t_631[k] = f_13 * sig_345[k]
                   + f_3 * pc_y[k] * skg_450[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pc_x, pc_y, pc_z, sig_330, sig_347, skf0_303, \
                         skf1_303, skg_450, skg_452, skg_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_10 * sig_330[k]
                   + f_3 * pc_z[k] * skg_450[k];

        t_633[k] = f_4 * skf0_303[k]
                   - f_5 * skf1_303[k]
                   + f_3 * pc_x[k] * skg_453[k];

        t_634[k] = f_13 * sig_347[k]
                   + f_3 * pc_y[k] * skg_452[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, pc_x, pc_y, pc_z, sig_333, sig_350, \
                         skf0_305, skf0_306, skf1_305, skf1_306, skg_453, skg_455, \
                         skg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_4 * skf0_305[k]
                   - f_5 * skf1_305[k]
                   + f_3 * pc_x[k] * skg_455[k];

        t_636[k] = f_6 * skf0_306[k]
                   - f_7 * skf1_306[k]
                   + f_3 * pc_x[k] * skg_456[k];

        t_637[k] = f_10 * sig_333[k]
                   + f_3 * pc_z[k] * skg_453[k];

        t_638[k] = f_13 * sig_350[k]
                   + f_3 * pc_y[k] * skg_455[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, t_644, pc_x, skf0_309, skf1_309, \
                         skg_459, skg_460, skg_461, skg_462, skg_463, \
                         skg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_6 * skf0_309[k]
                   - f_7 * skf1_309[k]
                   + f_3 * pc_x[k] * skg_459[k];

        t_640[k] = f_3 * pc_x[k] * skg_460[k];

        t_641[k] = f_3 * pc_x[k] * skg_461[k];

        t_642[k] = f_3 * pc_x[k] * skg_462[k];

        t_643[k] = f_3 * pc_x[k] * skg_463[k];

        t_644[k] = f_3 * pc_x[k] * skg_464[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_y, pc_z, sig_340, sig_355, sig_357, skf0_306, \
                         skf0_308, skf1_306, skf1_308, skg_460, \
                         skg_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_13 * sig_355[k]
                   + f_1 * skf0_306[k]
                   - f_2 * skf1_306[k]
                   + f_3 * pc_y[k] * skg_460[k];

        t_646[k] = f_10 * sig_340[k]
                   + f_3 * pc_z[k] * skg_460[k];

        t_647[k] = f_13 * sig_357[k]
                   + f_4 * skf0_308[k]
                   - f_5 * skf1_308[k]
                   + f_3 * pc_y[k] * skg_462[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_y, pc_z, sig_344, sig_358, sig_359, skf0_309, \
                         skf1_309, skg_463, skg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_13 * sig_358[k]
                   + f_6 * skf0_309[k]
                   - f_7 * skf1_309[k]
                   + f_3 * pc_y[k] * skg_463[k];

        t_649[k] = f_13 * sig_359[k]
                   + f_3 * pc_y[k] * skg_464[k];

        t_650[k] = f_10 * sig_344[k]
                   + f_1 * skf0_309[k]
                   - f_2 * skf1_309[k]
                   + f_3 * pc_z[k] * skg_464[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, pc_x, pc_y, pc_z, sig_345, sig_360, \
                         skf0_310, skf0_313, skf1_310, skf1_313, skg_465, \
                         skg_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_1 * skf0_310[k]
                   - f_2 * skf1_310[k]
                   + f_3 * pc_x[k] * skg_465[k];

        t_652[k] = f_14 * sig_360[k]
                   + f_3 * pc_y[k] * skg_465[k];

        t_653[k] = f_11 * sig_345[k]
                   + f_3 * pc_z[k] * skg_465[k];

        t_654[k] = f_4 * skf0_313[k]
                   - f_5 * skf1_313[k]
                   + f_3 * pc_x[k] * skg_468[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pc_x, pc_y, sig_362, skf0_315, skf0_316, \
                         skf1_315, skf1_316, skg_467, skg_470, \
                         skg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_14 * sig_362[k]
                   + f_3 * pc_y[k] * skg_467[k];

        t_656[k] = f_4 * skf0_315[k]
                   - f_5 * skf1_315[k]
                   + f_3 * pc_x[k] * skg_470[k];

        t_657[k] = f_6 * skf0_316[k]
                   - f_7 * skf1_316[k]
                   + f_3 * pc_x[k] * skg_471[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pc_x, pc_y, pc_z, sig_348, sig_365, \
                         skf0_319, skf1_319, skg_468, skg_470, skg_474, \
                         skg_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_11 * sig_348[k]
                   + f_3 * pc_z[k] * skg_468[k];

        t_659[k] = f_14 * sig_365[k]
                   + f_3 * pc_y[k] * skg_470[k];

        t_660[k] = f_6 * skf0_319[k]
                   - f_7 * skf1_319[k]
                   + f_3 * pc_x[k] * skg_474[k];

        t_661[k] = f_3 * pc_x[k] * skg_475[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, pc_x, pc_y, sig_370, skf0_316, \
                         skf1_316, skg_475, skg_476, skg_477, skg_478, \
                         skg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_3 * pc_x[k] * skg_476[k];

        t_663[k] = f_3 * pc_x[k] * skg_477[k];

        t_664[k] = f_3 * pc_x[k] * skg_478[k];

        t_665[k] = f_3 * pc_x[k] * skg_479[k];

        t_666[k] = f_14 * sig_370[k]
                   + f_1 * skf0_316[k]
                   - f_2 * skf1_316[k]
                   + f_3 * pc_y[k] * skg_475[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, pc_z, sig_355, sig_372, sig_373, skf0_318, \
                         skf0_319, skf1_318, skf1_319, skg_475, skg_477, \
                         skg_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_11 * sig_355[k]
                   + f_3 * pc_z[k] * skg_475[k];

        t_668[k] = f_14 * sig_372[k]
                   + f_4 * skf0_318[k]
                   - f_5 * skf1_318[k]
                   + f_3 * pc_y[k] * skg_477[k];

        t_669[k] = f_14 * sig_373[k]
                   + f_6 * skf0_319[k]
                   - f_7 * skf1_319[k]
                   + f_3 * pc_y[k] * skg_478[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pc_x, pc_y, pc_z, sig_359, sig_374, \
                         sig_375, skf0_319, skf0_320, skf1_319, skf1_320, skg_479, \
                         skg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * sig_374[k]
                   + f_3 * pc_y[k] * skg_479[k];

        t_671[k] = f_11 * sig_359[k]
                   + f_1 * skf0_319[k]
                   - f_2 * skf1_319[k]
                   + f_3 * pc_z[k] * skg_479[k];

        t_672[k] = f_1 * skf0_320[k]
                   - f_2 * skf1_320[k]
                   + f_3 * pc_x[k] * skg_480[k];

        t_673[k] = f_11 * sig_375[k]
                   + f_3 * pc_y[k] * skg_480[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, pc_x, pc_y, pc_z, sig_360, sig_377, skf0_323, \
                         skf1_323, skg_480, skg_482, skg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_14 * sig_360[k]
                   + f_3 * pc_z[k] * skg_480[k];

        t_675[k] = f_4 * skf0_323[k]
                   - f_5 * skf1_323[k]
                   + f_3 * pc_x[k] * skg_483[k];

        t_676[k] = f_11 * sig_377[k]
                   + f_3 * pc_y[k] * skg_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, pc_x, pc_y, pc_z, sig_363, sig_380, \
                         skf0_325, skf0_326, skf1_325, skf1_326, skg_483, skg_485, \
                         skg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_4 * skf0_325[k]
                   - f_5 * skf1_325[k]
                   + f_3 * pc_x[k] * skg_485[k];

        t_678[k] = f_6 * skf0_326[k]
                   - f_7 * skf1_326[k]
                   + f_3 * pc_x[k] * skg_486[k];

        t_679[k] = f_14 * sig_363[k]
                   + f_3 * pc_z[k] * skg_483[k];

        t_680[k] = f_11 * sig_380[k]
                   + f_3 * pc_y[k] * skg_485[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, t_685, t_686, pc_x, skf0_329, skf1_329, \
                         skg_489, skg_490, skg_491, skg_492, skg_493, \
                         skg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_6 * skf0_329[k]
                   - f_7 * skf1_329[k]
                   + f_3 * pc_x[k] * skg_489[k];

        t_682[k] = f_3 * pc_x[k] * skg_490[k];

        t_683[k] = f_3 * pc_x[k] * skg_491[k];

        t_684[k] = f_3 * pc_x[k] * skg_492[k];

        t_685[k] = f_3 * pc_x[k] * skg_493[k];

        t_686[k] = f_3 * pc_x[k] * skg_494[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pc_y, pc_z, sig_370, sig_385, sig_387, skf0_326, \
                         skf0_328, skf1_326, skf1_328, skg_490, \
                         skg_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_11 * sig_385[k]
                   + f_1 * skf0_326[k]
                   - f_2 * skf1_326[k]
                   + f_3 * pc_y[k] * skg_490[k];

        t_688[k] = f_14 * sig_370[k]
                   + f_3 * pc_z[k] * skg_490[k];

        t_689[k] = f_11 * sig_387[k]
                   + f_4 * skf0_328[k]
                   - f_5 * skf1_328[k]
                   + f_3 * pc_y[k] * skg_492[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, pc_y, pc_z, sig_374, sig_388, sig_389, skf0_329, \
                         skf1_329, skg_493, skg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_11 * sig_388[k]
                   + f_6 * skf0_329[k]
                   - f_7 * skf1_329[k]
                   + f_3 * pc_y[k] * skg_493[k];

        t_691[k] = f_11 * sig_389[k]
                   + f_3 * pc_y[k] * skg_494[k];

        t_692[k] = f_14 * sig_374[k]
                   + f_1 * skf0_329[k]
                   - f_2 * skf1_329[k]
                   + f_3 * pc_z[k] * skg_494[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, t_696, pc_x, pc_y, pc_z, sig_375, sig_390, \
                         skf0_330, skf0_333, skf1_330, skf1_333, skg_495, \
                         skg_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_1 * skf0_330[k]
                   - f_2 * skf1_330[k]
                   + f_3 * pc_x[k] * skg_495[k];

        t_694[k] = f_10 * sig_390[k]
                   + f_3 * pc_y[k] * skg_495[k];

        t_695[k] = f_13 * sig_375[k]
                   + f_3 * pc_z[k] * skg_495[k];

        t_696[k] = f_4 * skf0_333[k]
                   - f_5 * skf1_333[k]
                   + f_3 * pc_x[k] * skg_498[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_x, pc_y, sig_392, skf0_335, skf0_336, \
                         skf1_335, skf1_336, skg_497, skg_500, \
                         skg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_10 * sig_392[k]
                   + f_3 * pc_y[k] * skg_497[k];

        t_698[k] = f_4 * skf0_335[k]
                   - f_5 * skf1_335[k]
                   + f_3 * pc_x[k] * skg_500[k];

        t_699[k] = f_6 * skf0_336[k]
                   - f_7 * skf1_336[k]
                   + f_3 * pc_x[k] * skg_501[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, pc_x, pc_y, pc_z, sig_378, sig_395, \
                         skf0_339, skf1_339, skg_498, skg_500, skg_504, \
                         skg_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_13 * sig_378[k]
                   + f_3 * pc_z[k] * skg_498[k];

        t_701[k] = f_10 * sig_395[k]
                   + f_3 * pc_y[k] * skg_500[k];

        t_702[k] = f_6 * skf0_339[k]
                   - f_7 * skf1_339[k]
                   + f_3 * pc_x[k] * skg_504[k];

        t_703[k] = f_3 * pc_x[k] * skg_505[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, pc_x, pc_y, sig_400, skf0_336, \
                         skf1_336, skg_505, skg_506, skg_507, skg_508, \
                         skg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_3 * pc_x[k] * skg_506[k];

        t_705[k] = f_3 * pc_x[k] * skg_507[k];

        t_706[k] = f_3 * pc_x[k] * skg_508[k];

        t_707[k] = f_3 * pc_x[k] * skg_509[k];

        t_708[k] = f_10 * sig_400[k]
                   + f_1 * skf0_336[k]
                   - f_2 * skf1_336[k]
                   + f_3 * pc_y[k] * skg_505[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_y, pc_z, sig_385, sig_402, sig_403, skf0_338, \
                         skf0_339, skf1_338, skf1_339, skg_505, skg_507, \
                         skg_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_13 * sig_385[k]
                   + f_3 * pc_z[k] * skg_505[k];

        t_710[k] = f_10 * sig_402[k]
                   + f_4 * skf0_338[k]
                   - f_5 * skf1_338[k]
                   + f_3 * pc_y[k] * skg_507[k];

        t_711[k] = f_10 * sig_403[k]
                   + f_6 * skf0_339[k]
                   - f_7 * skf1_339[k]
                   + f_3 * pc_y[k] * skg_508[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pb_y, pc_y, pc_z, sih0_567, sig_389, \
                         sig_404, sig_405, sih1_567, skf0_339, skf1_339, skg_509, \
                         skg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_10 * sig_404[k]
                   + f_3 * pc_y[k] * skg_509[k];

        t_713[k] = f_13 * sig_389[k]
                   + f_1 * skf0_339[k]
                   - f_2 * skf1_339[k]
                   + f_3 * pc_z[k] * skg_509[k];

        t_714[k] = pb_y[k] * sih0_567[k]
                   - f_8 * pc_y[k] * sih1_567[k];

        t_715[k] = f_9 * sig_405[k]
                   + f_3 * pc_y[k] * skg_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, pc_x, pc_y, pc_z, sig_390, sig_407, skf0_343, \
                         skf1_343, skg_510, skg_512, skg_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_12 * sig_390[k]
                   + f_3 * pc_z[k] * skg_510[k];

        t_717[k] = f_4 * skf0_343[k]
                   - f_5 * skf1_343[k]
                   + f_3 * pc_x[k] * skg_513[k];

        t_718[k] = f_9 * sig_407[k]
                   + f_3 * pc_y[k] * skg_512[k];
    }

#pragma omp simd aligned(t_719, t_720, t_721, pb_y, pc_x, pc_y, pc_z, sih0_572, sig_393, \
                         sih1_572, skf0_346, skf1_346, skg_513, \
                         skg_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_719[k] = pb_y[k] * sih0_572[k]
                   - f_8 * pc_y[k] * sih1_572[k];

        t_720[k] = f_6 * skf0_346[k]
                   - f_7 * skf1_346[k]
                   + f_3 * pc_x[k] * skg_516[k];

        t_721[k] = f_12 * sig_393[k]
                   + f_3 * pc_z[k] * skg_513[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, t_726, pb_y, pc_x, pc_y, sih0_576, \
                         sig_410, sih1_576, skg_515, skg_520, skg_521, \
                         skg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_9 * sig_410[k]
                   + f_3 * pc_y[k] * skg_515[k];

        t_723[k] = pb_y[k] * sih0_576[k]
                   - f_8 * pc_y[k] * sih1_576[k];

        t_724[k] = f_3 * pc_x[k] * skg_520[k];

        t_725[k] = f_3 * pc_x[k] * skg_521[k];

        t_726[k] = f_3 * pc_x[k] * skg_522[k];
    }
}

static auto
compute_prim_skh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sih0,
                                                          const size_t sig, const size_t sih1,
                                                          const size_t skf0, const size_t skf1,
                                                          const size_t skg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sih0_582 = buffer.data(sih0 + 582);
    const auto *sih0_584 = buffer.data(sih0 + 584);
    const auto *sih0_585 = buffer.data(sih0 + 585);
    const auto *sih0_587 = buffer.data(sih0 + 587);

    const auto *sig_400 = buffer.data(sig + 400);
    const auto *sig_405 = buffer.data(sig + 405);
    const auto *sig_408 = buffer.data(sig + 408);
    const auto *sig_415 = buffer.data(sig + 415);
    const auto *sig_417 = buffer.data(sig + 417);
    const auto *sig_418 = buffer.data(sig + 418);
    const auto *sig_419 = buffer.data(sig + 419);

    const auto *sih1_582 = buffer.data(sih1 + 582);
    const auto *sih1_584 = buffer.data(sih1 + 584);
    const auto *sih1_585 = buffer.data(sih1 + 585);
    const auto *sih1_587 = buffer.data(sih1 + 587);

    const auto *skf0_350 = buffer.data(skf0 + 350);
    const auto *skf0_353 = buffer.data(skf0 + 353);
    const auto *skf0_355 = buffer.data(skf0 + 355);
    const auto *skf0_356 = buffer.data(skf0 + 356);
    const auto *skf0_358 = buffer.data(skf0 + 358);
    const auto *skf0_359 = buffer.data(skf0 + 359);

    const auto *skf1_350 = buffer.data(skf1 + 350);
    const auto *skf1_353 = buffer.data(skf1 + 353);
    const auto *skf1_355 = buffer.data(skf1 + 355);
    const auto *skf1_356 = buffer.data(skf1 + 356);
    const auto *skf1_358 = buffer.data(skf1 + 358);
    const auto *skf1_359 = buffer.data(skf1 + 359);

    const auto *skg_520 = buffer.data(skg + 520);
    const auto *skg_523 = buffer.data(skg + 523);
    const auto *skg_524 = buffer.data(skg + 524);
    const auto *skg_525 = buffer.data(skg + 525);
    const auto *skg_527 = buffer.data(skg + 527);
    const auto *skg_528 = buffer.data(skg + 528);
    const auto *skg_530 = buffer.data(skg + 530);
    const auto *skg_531 = buffer.data(skg + 531);
    const auto *skg_534 = buffer.data(skg + 534);
    const auto *skg_535 = buffer.data(skg + 535);
    const auto *skg_536 = buffer.data(skg + 536);
    const auto *skg_537 = buffer.data(skg + 537);
    const auto *skg_538 = buffer.data(skg + 538);
    const auto *skg_539 = buffer.data(skg + 539);

#pragma omp simd aligned(t_727, t_728, t_729, t_730, pb_y, pc_x, pc_y, pc_z, sih0_582, \
                         sig_400, sig_415, sih1_582, skg_520, skg_523, \
                         skg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_3 * pc_x[k] * skg_523[k];

        t_728[k] = f_3 * pc_x[k] * skg_524[k];

        t_729[k] = pb_y[k] * sih0_582[k]
                   + f_13 * sig_415[k]
                   - f_8 * pc_y[k] * sih1_582[k];

        t_730[k] = f_12 * sig_400[k]
                   + f_3 * pc_z[k] * skg_520[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, pb_y, pc_y, sih0_584, sih0_585, sih0_587, \
                         sig_417, sig_418, sig_419, sih1_584, sih1_585, sih1_587, \
                         skg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = pb_y[k] * sih0_584[k]
                   + f_11 * sig_417[k]
                   - f_8 * pc_y[k] * sih1_584[k];

        t_732[k] = pb_y[k] * sih0_585[k]
                   + f_10 * sig_418[k]
                   - f_8 * pc_y[k] * sih1_585[k];

        t_733[k] = f_9 * sig_419[k]
                   + f_3 * pc_y[k] * skg_524[k];

        t_734[k] = pb_y[k] * sih0_587[k]
                   - f_8 * pc_y[k] * sih1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, pc_x, pc_y, pc_z, sig_405, \
                         skf0_350, skf0_353, skf1_350, skf1_353, skg_525, skg_527, \
                         skg_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_1 * skf0_350[k]
                   - f_2 * skf1_350[k]
                   + f_3 * pc_x[k] * skg_525[k];

        t_736[k] = f_3 * pc_y[k] * skg_525[k];

        t_737[k] = f_0 * sig_405[k]
                   + f_3 * pc_z[k] * skg_525[k];

        t_738[k] = f_4 * skf0_353[k]
                   - f_5 * skf1_353[k]
                   + f_3 * pc_x[k] * skg_528[k];

        t_739[k] = f_3 * pc_y[k] * skg_527[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, pc_y, pc_z, sig_408, skf0_355, \
                         skf0_356, skf1_355, skf1_356, skg_528, skg_530, \
                         skg_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_4 * skf0_355[k]
                   - f_5 * skf1_355[k]
                   + f_3 * pc_x[k] * skg_530[k];

        t_741[k] = f_6 * skf0_356[k]
                   - f_7 * skf1_356[k]
                   + f_3 * pc_x[k] * skg_531[k];

        t_742[k] = f_0 * sig_408[k]
                   + f_3 * pc_z[k] * skg_528[k];

        t_743[k] = f_3 * pc_y[k] * skg_530[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, t_749, pc_x, skf0_359, skf1_359, \
                         skg_534, skg_535, skg_536, skg_537, skg_538, \
                         skg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_6 * skf0_359[k]
                   - f_7 * skf1_359[k]
                   + f_3 * pc_x[k] * skg_534[k];

        t_745[k] = f_3 * pc_x[k] * skg_535[k];

        t_746[k] = f_3 * pc_x[k] * skg_536[k];

        t_747[k] = f_3 * pc_x[k] * skg_537[k];

        t_748[k] = f_3 * pc_x[k] * skg_538[k];

        t_749[k] = f_3 * pc_x[k] * skg_539[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pc_y, pc_z, sig_415, skf0_356, skf0_358, \
                         skf0_359, skf1_356, skf1_358, skf1_359, skg_535, skg_537, \
                         skg_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * skf0_356[k]
                   - f_2 * skf1_356[k]
                   + f_3 * pc_y[k] * skg_535[k];

        t_751[k] = f_0 * sig_415[k]
                   + f_3 * pc_z[k] * skg_535[k];

        t_752[k] = f_4 * skf0_358[k]
                   - f_5 * skf1_358[k]
                   + f_3 * pc_y[k] * skg_537[k];

        t_753[k] = f_6 * skf0_359[k]
                   - f_7 * skf1_359[k]
                   + f_3 * pc_y[k] * skg_538[k];
    }

#pragma omp simd aligned(t_754, t_755, pc_y, pc_z, sig_419, skf0_359, skf1_359, \
                         skg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_3 * pc_y[k] * skg_539[k];

        t_755[k] = f_0 * sig_419[k]
                   + f_1 * skf0_359[k]
                   - f_2 * skf1_359[k]
                   + f_3 * pc_z[k] * skg_539[k];
    }
}

auto
compute_prim_skh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sih0, const size_t sig,
                                                   const size_t sih1, const size_t skf0,
                                                   const size_t skf1, const size_t skg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_skh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sih0, sig,
                                                              sih1, skf0, skf1, skg, ncols,
                                                              gamma, p, q);

    compute_prim_skh_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sih0, sig,
                                                              sih1, skf0, skf1, skg, ncols,
                                                              gamma, p, q);

    compute_prim_skh_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sih0, sig,
                                                              sih1, skf0, skf1, skg, ncols,
                                                              gamma, p, q);

    compute_prim_skh_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sih0, sig,
                                                              sih1, skf0, skf1, skg, ncols,
                                                              gamma, p, q);

    compute_prim_skh_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sih0, sig,
                                                              sih1, skf0, skf1, skg, ncols,
                                                              gamma, p, q);

    compute_prim_skh_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sih0, sig,
                                                              sih1, skf0, skf1, skg, ncols,
                                                              gamma, p, q);

    compute_prim_skh_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sih0, sig,
                                                              sih1, skf0, skf1, skg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
