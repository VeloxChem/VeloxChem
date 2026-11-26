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


#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gsh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsh0,
                                                          const size_t fsg, const size_t fsh1,
                                                          const size_t gsf0, const size_t gsf1,
                                                          const size_t gsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_14 = 2.5 / q;

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
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsh0_0 = buffer.data(fsh0 + 0);
    const auto *fsh0_3 = buffer.data(fsh0 + 3);
    const auto *fsh0_5 = buffer.data(fsh0 + 5);
    const auto *fsh0_6 = buffer.data(fsh0 + 6);
    const auto *fsh0_9 = buffer.data(fsh0 + 9);
    const auto *fsh0_15 = buffer.data(fsh0 + 15);
    const auto *fsh0_20 = buffer.data(fsh0 + 20);
    const auto *fsh0_24 = buffer.data(fsh0 + 24);
    const auto *fsh0_27 = buffer.data(fsh0 + 27);
    const auto *fsh0_36 = buffer.data(fsh0 + 36);
    const auto *fsh0_42 = buffer.data(fsh0 + 42);
    const auto *fsh0_47 = buffer.data(fsh0 + 47);
    const auto *fsh0_51 = buffer.data(fsh0 + 51);
    const auto *fsh0_62 = buffer.data(fsh0 + 62);
    const auto *fsh0_126 = buffer.data(fsh0 + 126);
    const auto *fsh0_129 = buffer.data(fsh0 + 129);

    const auto *fsg_0 = buffer.data(fsg + 0);
    const auto *fsg_1 = buffer.data(fsg + 1);
    const auto *fsg_2 = buffer.data(fsg + 2);
    const auto *fsg_3 = buffer.data(fsg + 3);
    const auto *fsg_5 = buffer.data(fsg + 5);
    const auto *fsg_10 = buffer.data(fsg + 10);
    const auto *fsg_12 = buffer.data(fsg + 12);
    const auto *fsg_14 = buffer.data(fsg + 14);
    const auto *fsg_15 = buffer.data(fsg + 15);
    const auto *fsg_18 = buffer.data(fsg + 18);
    const auto *fsg_20 = buffer.data(fsg + 20);
    const auto *fsg_25 = buffer.data(fsg + 25);
    const auto *fsg_27 = buffer.data(fsg + 27);
    const auto *fsg_28 = buffer.data(fsg + 28);
    const auto *fsg_29 = buffer.data(fsg + 29);
    const auto *fsg_30 = buffer.data(fsg + 30);
    const auto *fsg_32 = buffer.data(fsg + 32);
    const auto *fsg_35 = buffer.data(fsg + 35);
    const auto *fsg_40 = buffer.data(fsg + 40);
    const auto *fsg_41 = buffer.data(fsg + 41);
    const auto *fsg_42 = buffer.data(fsg + 42);
    const auto *fsg_43 = buffer.data(fsg + 43);
    const auto *fsg_44 = buffer.data(fsg + 44);
    const auto *fsg_45 = buffer.data(fsg + 45);
    const auto *fsg_48 = buffer.data(fsg + 48);
    const auto *fsg_51 = buffer.data(fsg + 51);
    const auto *fsg_55 = buffer.data(fsg + 55);
    const auto *fsg_57 = buffer.data(fsg + 57);
    const auto *fsg_58 = buffer.data(fsg + 58);
    const auto *fsg_59 = buffer.data(fsg + 59);
    const auto *fsg_70 = buffer.data(fsg + 70);
    const auto *fsg_71 = buffer.data(fsg + 71);
    const auto *fsg_72 = buffer.data(fsg + 72);
    const auto *fsg_73 = buffer.data(fsg + 73);
    const auto *fsg_74 = buffer.data(fsg + 74);
    const auto *fsg_75 = buffer.data(fsg + 75);
    const auto *fsg_80 = buffer.data(fsg + 80);
    const auto *fsg_84 = buffer.data(fsg + 84);
    const auto *fsg_85 = buffer.data(fsg + 85);
    const auto *fsg_86 = buffer.data(fsg + 86);
    const auto *fsg_87 = buffer.data(fsg + 87);
    const auto *fsg_89 = buffer.data(fsg + 89);
    const auto *fsg_90 = buffer.data(fsg + 90);
    const auto *fsg_93 = buffer.data(fsg + 93);

    const auto *fsh1_0 = buffer.data(fsh1 + 0);
    const auto *fsh1_3 = buffer.data(fsh1 + 3);
    const auto *fsh1_5 = buffer.data(fsh1 + 5);
    const auto *fsh1_6 = buffer.data(fsh1 + 6);
    const auto *fsh1_9 = buffer.data(fsh1 + 9);
    const auto *fsh1_15 = buffer.data(fsh1 + 15);
    const auto *fsh1_20 = buffer.data(fsh1 + 20);
    const auto *fsh1_24 = buffer.data(fsh1 + 24);
    const auto *fsh1_27 = buffer.data(fsh1 + 27);
    const auto *fsh1_36 = buffer.data(fsh1 + 36);
    const auto *fsh1_42 = buffer.data(fsh1 + 42);
    const auto *fsh1_47 = buffer.data(fsh1 + 47);
    const auto *fsh1_51 = buffer.data(fsh1 + 51);
    const auto *fsh1_62 = buffer.data(fsh1 + 62);
    const auto *fsh1_126 = buffer.data(fsh1 + 126);
    const auto *fsh1_129 = buffer.data(fsh1 + 129);

    const auto *gsf0_0 = buffer.data(gsf0 + 0);
    const auto *gsf0_1 = buffer.data(gsf0 + 1);
    const auto *gsf0_2 = buffer.data(gsf0 + 2);
    const auto *gsf0_6 = buffer.data(gsf0 + 6);
    const auto *gsf0_8 = buffer.data(gsf0 + 8);
    const auto *gsf0_9 = buffer.data(gsf0 + 9);
    const auto *gsf0_16 = buffer.data(gsf0 + 16);
    const auto *gsf0_17 = buffer.data(gsf0 + 17);
    const auto *gsf0_22 = buffer.data(gsf0 + 22);
    const auto *gsf0_27 = buffer.data(gsf0 + 27);
    const auto *gsf0_28 = buffer.data(gsf0 + 28);
    const auto *gsf0_29 = buffer.data(gsf0 + 29);
    const auto *gsf0_30 = buffer.data(gsf0 + 30);
    const auto *gsf0_32 = buffer.data(gsf0 + 32);
    const auto *gsf0_33 = buffer.data(gsf0 + 33);
    const auto *gsf0_36 = buffer.data(gsf0 + 36);
    const auto *gsf0_37 = buffer.data(gsf0 + 37);
    const auto *gsf0_39 = buffer.data(gsf0 + 39);
    const auto *gsf0_48 = buffer.data(gsf0 + 48);
    const auto *gsf0_49 = buffer.data(gsf0 + 49);
    const auto *gsf0_50 = buffer.data(gsf0 + 50);
    const auto *gsf0_51 = buffer.data(gsf0 + 51);
    const auto *gsf0_52 = buffer.data(gsf0 + 52);
    const auto *gsf0_55 = buffer.data(gsf0 + 55);
    const auto *gsf0_56 = buffer.data(gsf0 + 56);
    const auto *gsf0_57 = buffer.data(gsf0 + 57);
    const auto *gsf0_58 = buffer.data(gsf0 + 58);
    const auto *gsf0_59 = buffer.data(gsf0 + 59);

    const auto *gsf1_0 = buffer.data(gsf1 + 0);
    const auto *gsf1_1 = buffer.data(gsf1 + 1);
    const auto *gsf1_2 = buffer.data(gsf1 + 2);
    const auto *gsf1_6 = buffer.data(gsf1 + 6);
    const auto *gsf1_8 = buffer.data(gsf1 + 8);
    const auto *gsf1_9 = buffer.data(gsf1 + 9);
    const auto *gsf1_16 = buffer.data(gsf1 + 16);
    const auto *gsf1_17 = buffer.data(gsf1 + 17);
    const auto *gsf1_22 = buffer.data(gsf1 + 22);
    const auto *gsf1_27 = buffer.data(gsf1 + 27);
    const auto *gsf1_28 = buffer.data(gsf1 + 28);
    const auto *gsf1_29 = buffer.data(gsf1 + 29);
    const auto *gsf1_30 = buffer.data(gsf1 + 30);
    const auto *gsf1_32 = buffer.data(gsf1 + 32);
    const auto *gsf1_33 = buffer.data(gsf1 + 33);
    const auto *gsf1_36 = buffer.data(gsf1 + 36);
    const auto *gsf1_37 = buffer.data(gsf1 + 37);
    const auto *gsf1_39 = buffer.data(gsf1 + 39);
    const auto *gsf1_48 = buffer.data(gsf1 + 48);
    const auto *gsf1_49 = buffer.data(gsf1 + 49);
    const auto *gsf1_50 = buffer.data(gsf1 + 50);
    const auto *gsf1_51 = buffer.data(gsf1 + 51);
    const auto *gsf1_52 = buffer.data(gsf1 + 52);
    const auto *gsf1_55 = buffer.data(gsf1 + 55);
    const auto *gsf1_56 = buffer.data(gsf1 + 56);
    const auto *gsf1_57 = buffer.data(gsf1 + 57);
    const auto *gsf1_58 = buffer.data(gsf1 + 58);
    const auto *gsf1_59 = buffer.data(gsf1 + 59);

    const auto *gsg_0 = buffer.data(gsg + 0);
    const auto *gsg_1 = buffer.data(gsg + 1);
    const auto *gsg_2 = buffer.data(gsg + 2);
    const auto *gsg_3 = buffer.data(gsg + 3);
    const auto *gsg_5 = buffer.data(gsg + 5);
    const auto *gsg_6 = buffer.data(gsg + 6);
    const auto *gsg_9 = buffer.data(gsg + 9);
    const auto *gsg_10 = buffer.data(gsg + 10);
    const auto *gsg_12 = buffer.data(gsg + 12);
    const auto *gsg_13 = buffer.data(gsg + 13);
    const auto *gsg_14 = buffer.data(gsg + 14);
    const auto *gsg_15 = buffer.data(gsg + 15);
    const auto *gsg_16 = buffer.data(gsg + 16);
    const auto *gsg_18 = buffer.data(gsg + 18);
    const auto *gsg_20 = buffer.data(gsg + 20);
    const auto *gsg_21 = buffer.data(gsg + 21);
    const auto *gsg_25 = buffer.data(gsg + 25);
    const auto *gsg_26 = buffer.data(gsg + 26);
    const auto *gsg_27 = buffer.data(gsg + 27);
    const auto *gsg_28 = buffer.data(gsg + 28);
    const auto *gsg_29 = buffer.data(gsg + 29);
    const auto *gsg_30 = buffer.data(gsg + 30);
    const auto *gsg_32 = buffer.data(gsg + 32);
    const auto *gsg_34 = buffer.data(gsg + 34);
    const auto *gsg_35 = buffer.data(gsg + 35);
    const auto *gsg_39 = buffer.data(gsg + 39);
    const auto *gsg_40 = buffer.data(gsg + 40);
    const auto *gsg_41 = buffer.data(gsg + 41);
    const auto *gsg_42 = buffer.data(gsg + 42);
    const auto *gsg_43 = buffer.data(gsg + 43);
    const auto *gsg_44 = buffer.data(gsg + 44);
    const auto *gsg_45 = buffer.data(gsg + 45);
    const auto *gsg_46 = buffer.data(gsg + 46);
    const auto *gsg_47 = buffer.data(gsg + 47);
    const auto *gsg_48 = buffer.data(gsg + 48);
    const auto *gsg_50 = buffer.data(gsg + 50);
    const auto *gsg_51 = buffer.data(gsg + 51);
    const auto *gsg_55 = buffer.data(gsg + 55);
    const auto *gsg_56 = buffer.data(gsg + 56);
    const auto *gsg_57 = buffer.data(gsg + 57);
    const auto *gsg_58 = buffer.data(gsg + 58);
    const auto *gsg_59 = buffer.data(gsg + 59);
    const auto *gsg_60 = buffer.data(gsg + 60);
    const auto *gsg_62 = buffer.data(gsg + 62);
    const auto *gsg_63 = buffer.data(gsg + 63);
    const auto *gsg_65 = buffer.data(gsg + 65);
    const auto *gsg_70 = buffer.data(gsg + 70);
    const auto *gsg_71 = buffer.data(gsg + 71);
    const auto *gsg_72 = buffer.data(gsg + 72);
    const auto *gsg_73 = buffer.data(gsg + 73);
    const auto *gsg_74 = buffer.data(gsg + 74);
    const auto *gsg_75 = buffer.data(gsg + 75);
    const auto *gsg_76 = buffer.data(gsg + 76);
    const auto *gsg_77 = buffer.data(gsg + 77);
    const auto *gsg_78 = buffer.data(gsg + 78);
    const auto *gsg_79 = buffer.data(gsg + 79);
    const auto *gsg_80 = buffer.data(gsg + 80);
    const auto *gsg_84 = buffer.data(gsg + 84);
    const auto *gsg_85 = buffer.data(gsg + 85);
    const auto *gsg_86 = buffer.data(gsg + 86);
    const auto *gsg_87 = buffer.data(gsg + 87);
    const auto *gsg_88 = buffer.data(gsg + 88);
    const auto *gsg_89 = buffer.data(gsg + 89);
    const auto *gsg_90 = buffer.data(gsg + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fsg_0, gsf0_0, \
                         gsf1_0, gsg_0, gsg_1, gsg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fsg_0[k]
                 + f_1 * gsf0_0[k]
                 - f_2 * gsf1_0[k]
                 + f_3 * pc_x[k] * gsg_0[k];

        t_1[k] = f_3 * pc_y[k] * gsg_0[k];

        t_2[k] = f_3 * pc_z[k] * gsg_0[k];

        t_3[k] = f_4 * gsf0_0[k]
                 - f_5 * gsf1_0[k]
                 + f_3 * pc_y[k] * gsg_1[k];

        t_4[k] = f_3 * pc_y[k] * gsg_2[k];

        t_5[k] = f_4 * gsf0_0[k]
                 - f_5 * gsf1_0[k]
                 + f_3 * pc_z[k] * gsg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, fsg_10, gsf0_1, gsf0_2, \
                         gsf1_1, gsf1_2, gsg_3, gsg_5, gsg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * gsf0_1[k]
                 - f_7 * gsf1_1[k]
                 + f_3 * pc_y[k] * gsg_3[k];

        t_7[k] = f_3 * pc_z[k] * gsg_3[k];

        t_8[k] = f_3 * pc_y[k] * gsg_5[k];

        t_9[k] = f_6 * gsf0_2[k]
                 - f_7 * gsf1_2[k]
                 + f_3 * pc_z[k] * gsg_5[k];

        t_10[k] = f_0 * fsg_10[k]
                  + f_3 * pc_x[k] * gsg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, fsg_12, fsg_14, gsg_6, \
                         gsg_9, gsg_12, gsg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * gsg_6[k];

        t_12[k] = f_0 * fsg_12[k]
                  + f_3 * pc_x[k] * gsg_12[k];

        t_13[k] = f_3 * pc_y[k] * gsg_9[k];

        t_14[k] = f_0 * fsg_14[k]
                  + f_3 * pc_x[k] * gsg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, gsf0_6, gsf0_8, gsf0_9, gsf1_6, \
                         gsf1_8, gsf1_9, gsg_10, gsg_12, gsg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * gsf0_6[k]
                  - f_2 * gsf1_6[k]
                  + f_3 * pc_y[k] * gsg_10[k];

        t_16[k] = f_3 * pc_z[k] * gsg_10[k];

        t_17[k] = f_6 * gsf0_8[k]
                  - f_7 * gsf1_8[k]
                  + f_3 * pc_y[k] * gsg_12[k];

        t_18[k] = f_4 * gsf0_9[k]
                  - f_5 * gsf1_9[k]
                  + f_3 * pc_y[k] * gsg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, fsh0_0, fsg_0, \
                         fsh1_0, gsf0_9, gsf1_9, gsg_14, gsg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * gsg_14[k];

        t_20[k] = f_1 * gsf0_9[k]
                  - f_2 * gsf1_9[k]
                  + f_3 * pc_z[k] * gsg_14[k];

        t_21[k] = pa_y[k] * fsh0_0[k]
                  - f_8 * pc_y[k] * fsh1_0[k];

        t_22[k] = f_9 * fsg_0[k]
                  + f_3 * pc_y[k] * gsg_15[k];

        t_23[k] = f_3 * pc_z[k] * gsg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, fsh0_3, fsh0_5, fsh0_6, \
                         fsg_1, fsg_3, fsh1_3, fsh1_5, fsh1_6, gsg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * fsh0_3[k]
                  + f_10 * fsg_1[k]
                  - f_8 * pc_y[k] * fsh1_3[k];

        t_25[k] = f_3 * pc_z[k] * gsg_16[k];

        t_26[k] = pa_y[k] * fsh0_5[k]
                  - f_8 * pc_y[k] * fsh1_5[k];

        t_27[k] = pa_y[k] * fsh0_6[k]
                  + f_11 * fsg_3[k]
                  - f_8 * pc_y[k] * fsh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, fsh0_9, fsg_5, \
                         fsg_25, fsh1_9, gsg_18, gsg_20, gsg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * gsg_18[k];

        t_29[k] = f_9 * fsg_5[k]
                  + f_3 * pc_y[k] * gsg_20[k];

        t_30[k] = pa_y[k] * fsh0_9[k]
                  - f_8 * pc_y[k] * fsh1_9[k];

        t_31[k] = f_11 * fsg_25[k]
                  + f_3 * pc_x[k] * gsg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, fsg_27, fsg_28, fsg_29, gsg_21, \
                         gsg_27, gsg_28, gsg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * gsg_21[k];

        t_33[k] = f_11 * fsg_27[k]
                  + f_3 * pc_x[k] * gsg_27[k];

        t_34[k] = f_11 * fsg_28[k]
                  + f_3 * pc_x[k] * gsg_28[k];

        t_35[k] = f_11 * fsg_29[k]
                  + f_3 * pc_x[k] * gsg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, fsg_10, gsf0_16, gsf0_17, \
                         gsf1_16, gsf1_17, gsg_25, gsg_26, gsg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * fsg_10[k]
                  + f_1 * gsf0_16[k]
                  - f_2 * gsf1_16[k]
                  + f_3 * pc_y[k] * gsg_25[k];

        t_37[k] = f_3 * pc_z[k] * gsg_25[k];

        t_38[k] = f_4 * gsf0_16[k]
                  - f_5 * gsf1_16[k]
                  + f_3 * pc_z[k] * gsg_26[k];

        t_39[k] = f_6 * gsf0_17[k]
                  - f_7 * gsf1_17[k]
                  + f_3 * pc_z[k] * gsg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, fsh0_0, fsh0_20, \
                         fsg_14, fsh1_0, fsh1_20, gsg_29, gsg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * fsg_14[k]
                  + f_3 * pc_y[k] * gsg_29[k];

        t_41[k] = pa_y[k] * fsh0_20[k]
                  - f_8 * pc_y[k] * fsh1_20[k];

        t_42[k] = pa_z[k] * fsh0_0[k]
                  - f_8 * pc_z[k] * fsh1_0[k];

        t_43[k] = f_3 * pc_y[k] * gsg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, fsh0_3, fsh0_5, fsg_0, \
                         fsg_2, fsh1_3, fsh1_5, gsg_30, gsg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * fsg_0[k]
                  + f_3 * pc_z[k] * gsg_30[k];

        t_45[k] = pa_z[k] * fsh0_3[k]
                  - f_8 * pc_z[k] * fsh1_3[k];

        t_46[k] = f_3 * pc_y[k] * gsg_32[k];

        t_47[k] = pa_z[k] * fsh0_5[k]
                  + f_10 * fsg_2[k]
                  - f_8 * pc_z[k] * fsh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, fsh0_6, fsh0_9, fsg_5, \
                         fsh1_6, fsh1_9, gsf0_22, gsf1_22, gsg_34, \
                         gsg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * fsh0_6[k]
                  - f_8 * pc_z[k] * fsh1_6[k];

        t_49[k] = f_4 * gsf0_22[k]
                  - f_5 * gsf1_22[k]
                  + f_3 * pc_y[k] * gsg_34[k];

        t_50[k] = f_3 * pc_y[k] * gsg_35[k];

        t_51[k] = pa_z[k] * fsh0_9[k]
                  + f_11 * fsg_5[k]
                  - f_8 * pc_z[k] * fsh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, fsg_40, fsg_41, fsg_42, \
                         fsg_44, gsg_39, gsg_40, gsg_41, gsg_42, \
                         gsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_11 * fsg_40[k]
                  + f_3 * pc_x[k] * gsg_40[k];

        t_53[k] = f_11 * fsg_41[k]
                  + f_3 * pc_x[k] * gsg_41[k];

        t_54[k] = f_11 * fsg_42[k]
                  + f_3 * pc_x[k] * gsg_42[k];

        t_55[k] = f_3 * pc_y[k] * gsg_39[k];

        t_56[k] = f_11 * fsg_44[k]
                  + f_3 * pc_x[k] * gsg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, fsh0_15, fsh1_15, gsf0_27, \
                         gsf0_28, gsf1_27, gsf1_28, gsg_41, gsg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * fsh0_15[k]
                  - f_8 * pc_z[k] * fsh1_15[k];

        t_58[k] = f_12 * gsf0_27[k]
                  - f_13 * gsf1_27[k]
                  + f_3 * pc_y[k] * gsg_41[k];

        t_59[k] = f_6 * gsf0_28[k]
                  - f_7 * gsf1_28[k]
                  + f_3 * pc_y[k] * gsg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, fsg_14, fsg_45, gsf0_29, \
                         gsf0_30, gsf1_29, gsf1_30, gsg_43, gsg_44, \
                         gsg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * gsf0_29[k]
                  - f_5 * gsf1_29[k]
                  + f_3 * pc_y[k] * gsg_43[k];

        t_61[k] = f_3 * pc_y[k] * gsg_44[k];

        t_62[k] = f_9 * fsg_14[k]
                  + f_1 * gsf0_29[k]
                  - f_2 * gsf1_29[k]
                  + f_3 * pc_z[k] * gsg_44[k];

        t_63[k] = f_10 * fsg_45[k]
                  + f_1 * gsf0_30[k]
                  - f_2 * gsf1_30[k]
                  + f_3 * pc_x[k] * gsg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, fsg_15, fsg_48, gsf0_33, \
                         gsf1_33, gsg_45, gsg_46, gsg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * fsg_15[k]
                  + f_3 * pc_y[k] * gsg_45[k];

        t_65[k] = f_3 * pc_z[k] * gsg_45[k];

        t_66[k] = f_10 * fsg_48[k]
                  + f_6 * gsf0_33[k]
                  - f_7 * gsf1_33[k]
                  + f_3 * pc_x[k] * gsg_48[k];

        t_67[k] = f_3 * pc_z[k] * gsg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, fsg_51, gsf0_30, gsf0_36, gsf1_30, \
                         gsf1_36, gsg_47, gsg_48, gsg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * gsf0_30[k]
                  - f_5 * gsf1_30[k]
                  + f_3 * pc_z[k] * gsg_47[k];

        t_69[k] = f_10 * fsg_51[k]
                  + f_4 * gsf0_36[k]
                  - f_5 * gsf1_36[k]
                  + f_3 * pc_x[k] * gsg_51[k];

        t_70[k] = f_3 * pc_z[k] * gsg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, fsg_20, fsg_55, gsf0_32, \
                         gsf1_32, gsg_50, gsg_51, gsg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * fsg_20[k]
                  + f_3 * pc_y[k] * gsg_50[k];

        t_72[k] = f_6 * gsf0_32[k]
                  - f_7 * gsf1_32[k]
                  + f_3 * pc_z[k] * gsg_50[k];

        t_73[k] = f_10 * fsg_55[k]
                  + f_3 * pc_x[k] * gsg_55[k];

        t_74[k] = f_3 * pc_z[k] * gsg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, fsg_25, fsg_57, fsg_58, fsg_59, \
                         gsf0_36, gsf1_36, gsg_55, gsg_57, gsg_58, \
                         gsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * fsg_57[k]
                  + f_3 * pc_x[k] * gsg_57[k];

        t_76[k] = f_10 * fsg_58[k]
                  + f_3 * pc_x[k] * gsg_58[k];

        t_77[k] = f_10 * fsg_59[k]
                  + f_3 * pc_x[k] * gsg_59[k];

        t_78[k] = f_10 * fsg_25[k]
                  + f_1 * gsf0_36[k]
                  - f_2 * gsf1_36[k]
                  + f_3 * pc_y[k] * gsg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, fsg_29, gsf0_36, gsf0_37, \
                         gsf1_36, gsf1_37, gsg_55, gsg_56, gsg_57, \
                         gsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * gsg_55[k];

        t_80[k] = f_4 * gsf0_36[k]
                  - f_5 * gsf1_36[k]
                  + f_3 * pc_z[k] * gsg_56[k];

        t_81[k] = f_6 * gsf0_37[k]
                  - f_7 * gsf1_37[k]
                  + f_3 * pc_z[k] * gsg_57[k];

        t_82[k] = f_10 * fsg_29[k]
                  + f_3 * pc_y[k] * gsg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, fsh0_42, fsg_15, fsg_30, \
                         fsh1_42, gsf0_39, gsf1_39, gsg_59, gsg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * gsf0_39[k]
                  - f_2 * gsf1_39[k]
                  + f_3 * pc_z[k] * gsg_59[k];

        t_84[k] = pa_y[k] * fsh0_42[k]
                  - f_8 * pc_y[k] * fsh1_42[k];

        t_85[k] = f_9 * fsg_30[k]
                  + f_3 * pc_y[k] * gsg_60[k];

        t_86[k] = f_9 * fsg_15[k]
                  + f_3 * pc_z[k] * gsg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, fsh0_24, fsh0_27, \
                         fsh0_47, fsg_32, fsh1_24, fsh1_27, fsh1_47, \
                         gsg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * fsh0_24[k]
                  - f_8 * pc_z[k] * fsh1_24[k];

        t_88[k] = f_9 * fsg_32[k]
                  + f_3 * pc_y[k] * gsg_62[k];

        t_89[k] = pa_y[k] * fsh0_47[k]
                  - f_8 * pc_y[k] * fsh1_47[k];

        t_90[k] = pa_z[k] * fsh0_27[k]
                  - f_8 * pc_z[k] * fsh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, fsh0_51, fsg_18, \
                         fsg_35, fsg_70, fsh1_51, gsg_63, gsg_65, \
                         gsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * fsg_18[k]
                  + f_3 * pc_z[k] * gsg_63[k];

        t_92[k] = f_9 * fsg_35[k]
                  + f_3 * pc_y[k] * gsg_65[k];

        t_93[k] = pa_y[k] * fsh0_51[k]
                  - f_8 * pc_y[k] * fsh1_51[k];

        t_94[k] = f_10 * fsg_70[k]
                  + f_3 * pc_x[k] * gsg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, fsg_71, fsg_72, fsg_73, fsg_74, gsg_71, \
                         gsg_72, gsg_73, gsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_10 * fsg_71[k]
                  + f_3 * pc_x[k] * gsg_71[k];

        t_96[k] = f_10 * fsg_72[k]
                  + f_3 * pc_x[k] * gsg_72[k];

        t_97[k] = f_10 * fsg_73[k]
                  + f_3 * pc_x[k] * gsg_73[k];

        t_98[k] = f_10 * fsg_74[k]
                  + f_3 * pc_x[k] * gsg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, fsh0_36, fsg_25, fsg_42, \
                         fsh1_36, gsf0_48, gsf1_48, gsg_70, gsg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * fsh0_36[k]
                  - f_8 * pc_z[k] * fsh1_36[k];

        t_100[k] = f_9 * fsg_25[k]
                   + f_3 * pc_z[k] * gsg_70[k];

        t_101[k] = f_9 * fsg_42[k]
                   + f_6 * gsf0_48[k]
                   - f_7 * gsf1_48[k]
                   + f_3 * pc_y[k] * gsg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, fsh0_62, fsg_43, fsg_44, fsh1_62, \
                         gsf0_49, gsf1_49, gsg_73, gsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * fsg_43[k]
                   + f_4 * gsf0_49[k]
                   - f_5 * gsf1_49[k]
                   + f_3 * pc_y[k] * gsg_73[k];

        t_103[k] = f_9 * fsg_44[k]
                   + f_3 * pc_y[k] * gsg_74[k];

        t_104[k] = pa_y[k] * fsh0_62[k]
                   - f_8 * pc_y[k] * fsh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, fsg_30, fsg_75, \
                         gsf0_50, gsf1_50, gsg_75, gsg_76, gsg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * fsg_75[k]
                   + f_1 * gsf0_50[k]
                   - f_2 * gsf1_50[k]
                   + f_3 * pc_x[k] * gsg_75[k];

        t_106[k] = f_3 * pc_y[k] * gsg_75[k];

        t_107[k] = f_10 * fsg_30[k]
                   + f_3 * pc_z[k] * gsg_75[k];

        t_108[k] = f_4 * gsf0_50[k]
                   - f_5 * gsf1_50[k]
                   + f_3 * pc_y[k] * gsg_76[k];

        t_109[k] = f_3 * pc_y[k] * gsg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, fsg_80, gsf0_51, gsf0_52, \
                         gsf0_55, gsf1_51, gsf1_52, gsf1_55, gsg_78, gsg_79, \
                         gsg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_10 * fsg_80[k]
                   + f_6 * gsf0_55[k]
                   - f_7 * gsf1_55[k]
                   + f_3 * pc_x[k] * gsg_80[k];

        t_111[k] = f_6 * gsf0_51[k]
                   - f_7 * gsf1_51[k]
                   + f_3 * pc_y[k] * gsg_78[k];

        t_112[k] = f_4 * gsf0_52[k]
                   - f_5 * gsf1_52[k]
                   + f_3 * pc_y[k] * gsg_79[k];

        t_113[k] = f_3 * pc_y[k] * gsg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, fsg_84, fsg_85, fsg_86, fsg_87, \
                         gsf0_59, gsf1_59, gsg_84, gsg_85, gsg_86, \
                         gsg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_10 * fsg_84[k]
                   + f_4 * gsf0_59[k]
                   - f_5 * gsf1_59[k]
                   + f_3 * pc_x[k] * gsg_84[k];

        t_115[k] = f_10 * fsg_85[k]
                   + f_3 * pc_x[k] * gsg_85[k];

        t_116[k] = f_10 * fsg_86[k]
                   + f_3 * pc_x[k] * gsg_86[k];

        t_117[k] = f_10 * fsg_87[k]
                   + f_3 * pc_x[k] * gsg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, fsg_89, gsf0_56, gsf0_57, \
                         gsf1_56, gsf1_57, gsg_84, gsg_85, gsg_86, \
                         gsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * gsg_84[k];

        t_119[k] = f_10 * fsg_89[k]
                   + f_3 * pc_x[k] * gsg_89[k];

        t_120[k] = f_1 * gsf0_56[k]
                   - f_2 * gsf1_56[k]
                   + f_3 * pc_y[k] * gsg_85[k];

        t_121[k] = f_12 * gsf0_57[k]
                   - f_13 * gsf1_57[k]
                   + f_3 * pc_y[k] * gsg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, fsg_44, gsf0_58, gsf0_59, \
                         gsf1_58, gsf1_59, gsg_87, gsg_88, gsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * gsf0_58[k]
                   - f_7 * gsf1_58[k]
                   + f_3 * pc_y[k] * gsg_87[k];

        t_123[k] = f_4 * gsf0_59[k]
                   - f_5 * gsf1_59[k]
                   + f_3 * pc_y[k] * gsg_88[k];

        t_124[k] = f_3 * pc_y[k] * gsg_89[k];

        t_125[k] = f_10 * fsg_44[k]
                   + f_1 * gsf0_59[k]
                   - f_2 * gsf1_59[k]
                   + f_3 * pc_z[k] * gsg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_x, pc_x, pc_y, pc_z, fsh0_126, \
                         fsh0_129, fsg_45, fsg_90, fsg_93, fsh1_126, fsh1_129, \
                         gsg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_x[k] * fsh0_126[k]
                   + f_14 * fsg_90[k]
                   - f_8 * pc_x[k] * fsh1_126[k];

        t_127[k] = f_11 * fsg_45[k]
                   + f_3 * pc_y[k] * gsg_90[k];

        t_128[k] = f_3 * pc_z[k] * gsg_90[k];

        t_129[k] = pa_x[k] * fsh0_129[k]
                   + f_11 * fsg_93[k]
                   - f_8 * pc_x[k] * fsh1_129[k];
    }
}

static auto
compute_prim_gsh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsh0,
                                                          const size_t fsg, const size_t fsh1,
                                                          const size_t gsf0, const size_t gsf1,
                                                          const size_t gsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_14 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsh0_63 = buffer.data(fsh0 + 63);
    const auto *fsh0_66 = buffer.data(fsh0 + 66);
    const auto *fsh0_69 = buffer.data(fsh0 + 69);
    const auto *fsh0_105 = buffer.data(fsh0 + 105);
    const auto *fsh0_110 = buffer.data(fsh0 + 110);
    const auto *fsh0_114 = buffer.data(fsh0 + 114);
    const auto *fsh0_126 = buffer.data(fsh0 + 126);
    const auto *fsh0_127 = buffer.data(fsh0 + 127);
    const auto *fsh0_129 = buffer.data(fsh0 + 129);
    const auto *fsh0_132 = buffer.data(fsh0 + 132);
    const auto *fsh0_141 = buffer.data(fsh0 + 141);
    const auto *fsh0_143 = buffer.data(fsh0 + 143);
    const auto *fsh0_144 = buffer.data(fsh0 + 144);
    const auto *fsh0_146 = buffer.data(fsh0 + 146);
    const auto *fsh0_152 = buffer.data(fsh0 + 152);
    const auto *fsh0_156 = buffer.data(fsh0 + 156);
    const auto *fsh0_162 = buffer.data(fsh0 + 162);
    const auto *fsh0_164 = buffer.data(fsh0 + 164);
    const auto *fsh0_165 = buffer.data(fsh0 + 165);
    const auto *fsh0_167 = buffer.data(fsh0 + 167);
    const auto *fsh0_171 = buffer.data(fsh0 + 171);
    const auto *fsh0_174 = buffer.data(fsh0 + 174);
    const auto *fsh0_183 = buffer.data(fsh0 + 183);
    const auto *fsh0_185 = buffer.data(fsh0 + 185);
    const auto *fsh0_186 = buffer.data(fsh0 + 186);
    const auto *fsh0_188 = buffer.data(fsh0 + 188);
    const auto *fsh0_189 = buffer.data(fsh0 + 189);
    const auto *fsh0_194 = buffer.data(fsh0 + 194);
    const auto *fsh0_198 = buffer.data(fsh0 + 198);
    const auto *fsh0_204 = buffer.data(fsh0 + 204);
    const auto *fsh0_205 = buffer.data(fsh0 + 205);
    const auto *fsh0_206 = buffer.data(fsh0 + 206);
    const auto *fsh0_207 = buffer.data(fsh0 + 207);
    const auto *fsh0_209 = buffer.data(fsh0 + 209);

    const auto *fsg_45 = buffer.data(fsg + 45);
    const auto *fsg_48 = buffer.data(fsg + 48);
    const auto *fsg_50 = buffer.data(fsg + 50);
    const auto *fsg_55 = buffer.data(fsg + 55);
    const auto *fsg_59 = buffer.data(fsg + 59);
    const auto *fsg_60 = buffer.data(fsg + 60);
    const auto *fsg_62 = buffer.data(fsg + 62);
    const auto *fsg_63 = buffer.data(fsg + 63);
    const auto *fsg_65 = buffer.data(fsg + 65);
    const auto *fsg_70 = buffer.data(fsg + 70);
    const auto *fsg_74 = buffer.data(fsg + 74);
    const auto *fsg_75 = buffer.data(fsg + 75);
    const auto *fsg_77 = buffer.data(fsg + 77);
    const auto *fsg_80 = buffer.data(fsg + 80);
    const auto *fsg_89 = buffer.data(fsg + 89);
    const auto *fsg_96 = buffer.data(fsg + 96);
    const auto *fsg_100 = buffer.data(fsg + 100);
    const auto *fsg_101 = buffer.data(fsg + 101);
    const auto *fsg_102 = buffer.data(fsg + 102);
    const auto *fsg_103 = buffer.data(fsg + 103);
    const auto *fsg_104 = buffer.data(fsg + 104);
    const auto *fsg_110 = buffer.data(fsg + 110);
    const auto *fsg_114 = buffer.data(fsg + 114);
    const auto *fsg_115 = buffer.data(fsg + 115);
    const auto *fsg_116 = buffer.data(fsg + 116);
    const auto *fsg_117 = buffer.data(fsg + 117);
    const auto *fsg_118 = buffer.data(fsg + 118);
    const auto *fsg_119 = buffer.data(fsg + 119);
    const auto *fsg_123 = buffer.data(fsg + 123);
    const auto *fsg_126 = buffer.data(fsg + 126);
    const auto *fsg_130 = buffer.data(fsg + 130);
    const auto *fsg_131 = buffer.data(fsg + 131);
    const auto *fsg_132 = buffer.data(fsg + 132);
    const auto *fsg_133 = buffer.data(fsg + 133);
    const auto *fsg_134 = buffer.data(fsg + 134);
    const auto *fsg_135 = buffer.data(fsg + 135);
    const auto *fsg_140 = buffer.data(fsg + 140);
    const auto *fsg_144 = buffer.data(fsg + 144);
    const auto *fsg_145 = buffer.data(fsg + 145);
    const auto *fsg_146 = buffer.data(fsg + 146);
    const auto *fsg_147 = buffer.data(fsg + 147);
    const auto *fsg_149 = buffer.data(fsg + 149);

    const auto *fsh1_63 = buffer.data(fsh1 + 63);
    const auto *fsh1_66 = buffer.data(fsh1 + 66);
    const auto *fsh1_69 = buffer.data(fsh1 + 69);
    const auto *fsh1_105 = buffer.data(fsh1 + 105);
    const auto *fsh1_110 = buffer.data(fsh1 + 110);
    const auto *fsh1_114 = buffer.data(fsh1 + 114);
    const auto *fsh1_126 = buffer.data(fsh1 + 126);
    const auto *fsh1_127 = buffer.data(fsh1 + 127);
    const auto *fsh1_129 = buffer.data(fsh1 + 129);
    const auto *fsh1_132 = buffer.data(fsh1 + 132);
    const auto *fsh1_141 = buffer.data(fsh1 + 141);
    const auto *fsh1_143 = buffer.data(fsh1 + 143);
    const auto *fsh1_144 = buffer.data(fsh1 + 144);
    const auto *fsh1_146 = buffer.data(fsh1 + 146);
    const auto *fsh1_152 = buffer.data(fsh1 + 152);
    const auto *fsh1_156 = buffer.data(fsh1 + 156);
    const auto *fsh1_162 = buffer.data(fsh1 + 162);
    const auto *fsh1_164 = buffer.data(fsh1 + 164);
    const auto *fsh1_165 = buffer.data(fsh1 + 165);
    const auto *fsh1_167 = buffer.data(fsh1 + 167);
    const auto *fsh1_171 = buffer.data(fsh1 + 171);
    const auto *fsh1_174 = buffer.data(fsh1 + 174);
    const auto *fsh1_183 = buffer.data(fsh1 + 183);
    const auto *fsh1_185 = buffer.data(fsh1 + 185);
    const auto *fsh1_186 = buffer.data(fsh1 + 186);
    const auto *fsh1_188 = buffer.data(fsh1 + 188);
    const auto *fsh1_189 = buffer.data(fsh1 + 189);
    const auto *fsh1_194 = buffer.data(fsh1 + 194);
    const auto *fsh1_198 = buffer.data(fsh1 + 198);
    const auto *fsh1_204 = buffer.data(fsh1 + 204);
    const auto *fsh1_205 = buffer.data(fsh1 + 205);
    const auto *fsh1_206 = buffer.data(fsh1 + 206);
    const auto *fsh1_207 = buffer.data(fsh1 + 207);
    const auto *fsh1_209 = buffer.data(fsh1 + 209);

    const auto *gsf0_60 = buffer.data(gsf0 + 60);
    const auto *gsf0_62 = buffer.data(gsf0 + 62);
    const auto *gsf0_90 = buffer.data(gsf0 + 90);
    const auto *gsf0_91 = buffer.data(gsf0 + 91);
    const auto *gsf0_92 = buffer.data(gsf0 + 92);
    const auto *gsf0_100 = buffer.data(gsf0 + 100);
    const auto *gsf0_101 = buffer.data(gsf0 + 101);
    const auto *gsf0_103 = buffer.data(gsf0 + 103);
    const auto *gsf0_105 = buffer.data(gsf0 + 105);
    const auto *gsf0_106 = buffer.data(gsf0 + 106);
    const auto *gsf0_107 = buffer.data(gsf0 + 107);
    const auto *gsf0_108 = buffer.data(gsf0 + 108);
    const auto *gsf0_109 = buffer.data(gsf0 + 109);
    const auto *gsf0_112 = buffer.data(gsf0 + 112);
    const auto *gsf0_114 = buffer.data(gsf0 + 114);
    const auto *gsf0_115 = buffer.data(gsf0 + 115);
    const auto *gsf0_117 = buffer.data(gsf0 + 117);
    const auto *gsf0_118 = buffer.data(gsf0 + 118);
    const auto *gsf0_119 = buffer.data(gsf0 + 119);
    const auto *gsf0_120 = buffer.data(gsf0 + 120);
    const auto *gsf0_121 = buffer.data(gsf0 + 121);
    const auto *gsf0_122 = buffer.data(gsf0 + 122);
    const auto *gsf0_123 = buffer.data(gsf0 + 123);
    const auto *gsf0_124 = buffer.data(gsf0 + 124);

    const auto *gsf1_60 = buffer.data(gsf1 + 60);
    const auto *gsf1_62 = buffer.data(gsf1 + 62);
    const auto *gsf1_90 = buffer.data(gsf1 + 90);
    const auto *gsf1_91 = buffer.data(gsf1 + 91);
    const auto *gsf1_92 = buffer.data(gsf1 + 92);
    const auto *gsf1_100 = buffer.data(gsf1 + 100);
    const auto *gsf1_101 = buffer.data(gsf1 + 101);
    const auto *gsf1_103 = buffer.data(gsf1 + 103);
    const auto *gsf1_105 = buffer.data(gsf1 + 105);
    const auto *gsf1_106 = buffer.data(gsf1 + 106);
    const auto *gsf1_107 = buffer.data(gsf1 + 107);
    const auto *gsf1_108 = buffer.data(gsf1 + 108);
    const auto *gsf1_109 = buffer.data(gsf1 + 109);
    const auto *gsf1_112 = buffer.data(gsf1 + 112);
    const auto *gsf1_114 = buffer.data(gsf1 + 114);
    const auto *gsf1_115 = buffer.data(gsf1 + 115);
    const auto *gsf1_117 = buffer.data(gsf1 + 117);
    const auto *gsf1_118 = buffer.data(gsf1 + 118);
    const auto *gsf1_119 = buffer.data(gsf1 + 119);
    const auto *gsf1_120 = buffer.data(gsf1 + 120);
    const auto *gsf1_121 = buffer.data(gsf1 + 121);
    const auto *gsf1_122 = buffer.data(gsf1 + 122);
    const auto *gsf1_123 = buffer.data(gsf1 + 123);
    const auto *gsf1_124 = buffer.data(gsf1 + 124);

    const auto *gsg_91 = buffer.data(gsg + 91);
    const auto *gsg_92 = buffer.data(gsg + 92);
    const auto *gsg_93 = buffer.data(gsg + 93);
    const auto *gsg_95 = buffer.data(gsg + 95);
    const auto *gsg_96 = buffer.data(gsg + 96);
    const auto *gsg_100 = buffer.data(gsg + 100);
    const auto *gsg_102 = buffer.data(gsg + 102);
    const auto *gsg_103 = buffer.data(gsg + 103);
    const auto *gsg_104 = buffer.data(gsg + 104);
    const auto *gsg_105 = buffer.data(gsg + 105);
    const auto *gsg_107 = buffer.data(gsg + 107);
    const auto *gsg_108 = buffer.data(gsg + 108);
    const auto *gsg_110 = buffer.data(gsg + 110);
    const auto *gsg_115 = buffer.data(gsg + 115);
    const auto *gsg_116 = buffer.data(gsg + 116);
    const auto *gsg_117 = buffer.data(gsg + 117);
    const auto *gsg_118 = buffer.data(gsg + 118);
    const auto *gsg_119 = buffer.data(gsg + 119);
    const auto *gsg_120 = buffer.data(gsg + 120);
    const auto *gsg_122 = buffer.data(gsg + 122);
    const auto *gsg_123 = buffer.data(gsg + 123);
    const auto *gsg_125 = buffer.data(gsg + 125);
    const auto *gsg_130 = buffer.data(gsg + 130);
    const auto *gsg_131 = buffer.data(gsg + 131);
    const auto *gsg_132 = buffer.data(gsg + 132);
    const auto *gsg_133 = buffer.data(gsg + 133);
    const auto *gsg_134 = buffer.data(gsg + 134);
    const auto *gsg_135 = buffer.data(gsg + 135);
    const auto *gsg_136 = buffer.data(gsg + 136);
    const auto *gsg_137 = buffer.data(gsg + 137);
    const auto *gsg_138 = buffer.data(gsg + 138);
    const auto *gsg_139 = buffer.data(gsg + 139);
    const auto *gsg_140 = buffer.data(gsg + 140);
    const auto *gsg_144 = buffer.data(gsg + 144);
    const auto *gsg_145 = buffer.data(gsg + 145);
    const auto *gsg_146 = buffer.data(gsg + 146);
    const auto *gsg_147 = buffer.data(gsg + 147);
    const auto *gsg_149 = buffer.data(gsg + 149);
    const auto *gsg_150 = buffer.data(gsg + 150);
    const auto *gsg_151 = buffer.data(gsg + 151);
    const auto *gsg_153 = buffer.data(gsg + 153);
    const auto *gsg_155 = buffer.data(gsg + 155);
    const auto *gsg_156 = buffer.data(gsg + 156);
    const auto *gsg_158 = buffer.data(gsg + 158);
    const auto *gsg_159 = buffer.data(gsg + 159);
    const auto *gsg_160 = buffer.data(gsg + 160);
    const auto *gsg_161 = buffer.data(gsg + 161);
    const auto *gsg_162 = buffer.data(gsg + 162);
    const auto *gsg_163 = buffer.data(gsg + 163);
    const auto *gsg_164 = buffer.data(gsg + 164);
    const auto *gsg_167 = buffer.data(gsg + 167);
    const auto *gsg_169 = buffer.data(gsg + 169);
    const auto *gsg_170 = buffer.data(gsg + 170);
    const auto *gsg_172 = buffer.data(gsg + 172);
    const auto *gsg_173 = buffer.data(gsg + 173);
    const auto *gsg_174 = buffer.data(gsg + 174);
    const auto *gsg_175 = buffer.data(gsg + 175);
    const auto *gsg_176 = buffer.data(gsg + 176);
    const auto *gsg_177 = buffer.data(gsg + 177);
    const auto *gsg_178 = buffer.data(gsg + 178);
    const auto *gsg_179 = buffer.data(gsg + 179);
    const auto *gsg_180 = buffer.data(gsg + 180);
    const auto *gsg_181 = buffer.data(gsg + 181);
    const auto *gsg_182 = buffer.data(gsg + 182);
    const auto *gsg_183 = buffer.data(gsg + 183);
    const auto *gsg_184 = buffer.data(gsg + 184);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_x, pc_x, pc_z, fsh0_132, fsg_96, \
                         fsh1_132, gsf0_60, gsf1_60, gsg_91, gsg_92, \
                         gsg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * gsg_91[k];

        t_131[k] = f_4 * gsf0_60[k]
                   - f_5 * gsf1_60[k]
                   + f_3 * pc_z[k] * gsg_92[k];

        t_132[k] = pa_x[k] * fsh0_132[k]
                   + f_10 * fsg_96[k]
                   - f_8 * pc_x[k] * fsh1_132[k];

        t_133[k] = f_3 * pc_z[k] * gsg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, fsg_50, fsg_100, \
                         gsf0_62, gsf1_62, gsg_95, gsg_96, gsg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * fsg_50[k]
                   + f_3 * pc_y[k] * gsg_95[k];

        t_135[k] = f_6 * gsf0_62[k]
                   - f_7 * gsf1_62[k]
                   + f_3 * pc_z[k] * gsg_95[k];

        t_136[k] = f_9 * fsg_100[k]
                   + f_3 * pc_x[k] * gsg_100[k];

        t_137[k] = f_3 * pc_z[k] * gsg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_x, pc_x, fsh0_141, fsg_102, fsg_103, \
                         fsg_104, fsh1_141, gsg_102, gsg_103, gsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_9 * fsg_102[k]
                   + f_3 * pc_x[k] * gsg_102[k];

        t_139[k] = f_9 * fsg_103[k]
                   + f_3 * pc_x[k] * gsg_103[k];

        t_140[k] = f_9 * fsg_104[k]
                   + f_3 * pc_x[k] * gsg_104[k];

        t_141[k] = pa_x[k] * fsh0_141[k]
                   - f_8 * pc_x[k] * fsh1_141[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pa_x, pc_x, pc_y, pc_z, fsh0_143, \
                         fsh0_144, fsg_59, fsh1_143, fsh1_144, gsg_100, \
                         gsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * gsg_100[k];

        t_143[k] = pa_x[k] * fsh0_143[k]
                   - f_8 * pc_x[k] * fsh1_143[k];

        t_144[k] = pa_x[k] * fsh0_144[k]
                   - f_8 * pc_x[k] * fsh1_144[k];

        t_145[k] = f_11 * fsg_59[k]
                   + f_3 * pc_y[k] * gsg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pa_z, pc_x, pc_y, pc_z, fsh0_63, \
                         fsh0_146, fsg_45, fsg_60, fsh1_63, fsh1_146, \
                         gsg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_x[k] * fsh0_146[k]
                   - f_8 * pc_x[k] * fsh1_146[k];

        t_147[k] = pa_z[k] * fsh0_63[k]
                   - f_8 * pc_z[k] * fsh1_63[k];

        t_148[k] = f_10 * fsg_60[k]
                   + f_3 * pc_y[k] * gsg_105[k];

        t_149[k] = f_9 * fsg_45[k]
                   + f_3 * pc_z[k] * gsg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_x, pa_z, pc_x, pc_y, pc_z, fsh0_66, fsh0_152, \
                         fsg_62, fsg_110, fsh1_66, fsh1_152, gsg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * fsh0_66[k]
                   - f_8 * pc_z[k] * fsh1_66[k];

        t_151[k] = f_10 * fsg_62[k]
                   + f_3 * pc_y[k] * gsg_107[k];

        t_152[k] = pa_x[k] * fsh0_152[k]
                   + f_11 * fsg_110[k]
                   - f_8 * pc_x[k] * fsh1_152[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, fsh0_69, fsg_48, fsg_65, \
                         fsh1_69, gsg_108, gsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * fsh0_69[k]
                   - f_8 * pc_z[k] * fsh1_69[k];

        t_154[k] = f_9 * fsg_48[k]
                   + f_3 * pc_z[k] * gsg_108[k];

        t_155[k] = f_10 * fsg_65[k]
                   + f_3 * pc_y[k] * gsg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pa_x, pc_x, fsh0_156, fsg_114, fsg_115, \
                         fsg_116, fsg_117, fsh1_156, gsg_115, gsg_116, \
                         gsg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_x[k] * fsh0_156[k]
                   + f_10 * fsg_114[k]
                   - f_8 * pc_x[k] * fsh1_156[k];

        t_157[k] = f_9 * fsg_115[k]
                   + f_3 * pc_x[k] * gsg_115[k];

        t_158[k] = f_9 * fsg_116[k]
                   + f_3 * pc_x[k] * gsg_116[k];

        t_159[k] = f_9 * fsg_117[k]
                   + f_3 * pc_x[k] * gsg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pc_x, pc_z, fsh0_162, fsg_55, \
                         fsg_118, fsg_119, fsh1_162, gsg_115, gsg_118, \
                         gsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_9 * fsg_118[k]
                   + f_3 * pc_x[k] * gsg_118[k];

        t_161[k] = f_9 * fsg_119[k]
                   + f_3 * pc_x[k] * gsg_119[k];

        t_162[k] = pa_x[k] * fsh0_162[k]
                   - f_8 * pc_x[k] * fsh1_162[k];

        t_163[k] = f_9 * fsg_55[k]
                   + f_3 * pc_z[k] * gsg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_x, pc_x, pc_y, fsh0_164, fsh0_165, \
                         fsh0_167, fsg_74, fsh1_164, fsh1_165, fsh1_167, \
                         gsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pa_x[k] * fsh0_164[k]
                   - f_8 * pc_x[k] * fsh1_164[k];

        t_165[k] = pa_x[k] * fsh0_165[k]
                   - f_8 * pc_x[k] * fsh1_165[k];

        t_166[k] = f_10 * fsg_74[k]
                   + f_3 * pc_y[k] * gsg_119[k];

        t_167[k] = pa_x[k] * fsh0_167[k]
                   - f_8 * pc_x[k] * fsh1_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_y, pc_y, pc_z, fsh0_105, fsg_60, fsg_75, \
                         fsh1_105, gsg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_y[k] * fsh0_105[k]
                   - f_8 * pc_y[k] * fsh1_105[k];

        t_169[k] = f_9 * fsg_75[k]
                   + f_3 * pc_y[k] * gsg_120[k];

        t_170[k] = f_10 * fsg_60[k]
                   + f_3 * pc_z[k] * gsg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_x, pa_y, pc_x, pc_y, fsh0_110, fsh0_171, \
                         fsg_77, fsg_123, fsh1_110, fsh1_171, gsg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_x[k] * fsh0_171[k]
                   + f_11 * fsg_123[k]
                   - f_8 * pc_x[k] * fsh1_171[k];

        t_172[k] = f_9 * fsg_77[k]
                   + f_3 * pc_y[k] * gsg_122[k];

        t_173[k] = pa_y[k] * fsh0_110[k]
                   - f_8 * pc_y[k] * fsh1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pa_x, pc_x, pc_y, pc_z, fsh0_174, fsg_63, \
                         fsg_80, fsg_126, fsh1_174, gsg_123, gsg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_x[k] * fsh0_174[k]
                   + f_10 * fsg_126[k]
                   - f_8 * pc_x[k] * fsh1_174[k];

        t_175[k] = f_10 * fsg_63[k]
                   + f_3 * pc_z[k] * gsg_123[k];

        t_176[k] = f_9 * fsg_80[k]
                   + f_3 * pc_y[k] * gsg_125[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_y, pc_x, pc_y, fsh0_114, fsg_130, \
                         fsg_131, fsg_132, fsh1_114, gsg_130, gsg_131, \
                         gsg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pa_y[k] * fsh0_114[k]
                   - f_8 * pc_y[k] * fsh1_114[k];

        t_178[k] = f_9 * fsg_130[k]
                   + f_3 * pc_x[k] * gsg_130[k];

        t_179[k] = f_9 * fsg_131[k]
                   + f_3 * pc_x[k] * gsg_131[k];

        t_180[k] = f_9 * fsg_132[k]
                   + f_3 * pc_x[k] * gsg_132[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_x, pc_x, pc_z, fsh0_183, fsg_70, \
                         fsg_133, fsg_134, fsh1_183, gsg_130, gsg_133, \
                         gsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_9 * fsg_133[k]
                   + f_3 * pc_x[k] * gsg_133[k];

        t_182[k] = f_9 * fsg_134[k]
                   + f_3 * pc_x[k] * gsg_134[k];

        t_183[k] = pa_x[k] * fsh0_183[k]
                   - f_8 * pc_x[k] * fsh1_183[k];

        t_184[k] = f_10 * fsg_70[k]
                   + f_3 * pc_z[k] * gsg_130[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_x, pc_x, pc_y, fsh0_185, fsh0_186, \
                         fsh0_188, fsg_89, fsh1_185, fsh1_186, fsh1_188, \
                         gsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pa_x[k] * fsh0_185[k]
                   - f_8 * pc_x[k] * fsh1_185[k];

        t_186[k] = pa_x[k] * fsh0_186[k]
                   - f_8 * pc_x[k] * fsh1_186[k];

        t_187[k] = f_9 * fsg_89[k]
                   + f_3 * pc_y[k] * gsg_134[k];

        t_188[k] = pa_x[k] * fsh0_188[k]
                   - f_8 * pc_x[k] * fsh1_188[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pc_x, pc_y, pc_z, fsh0_189, fsg_75, \
                         fsg_135, fsh1_189, gsf0_90, gsf1_90, gsg_135, \
                         gsg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pa_x[k] * fsh0_189[k]
                   + f_14 * fsg_135[k]
                   - f_8 * pc_x[k] * fsh1_189[k];

        t_190[k] = f_3 * pc_y[k] * gsg_135[k];

        t_191[k] = f_11 * fsg_75[k]
                   + f_3 * pc_z[k] * gsg_135[k];

        t_192[k] = f_4 * gsf0_90[k]
                   - f_5 * gsf1_90[k]
                   + f_3 * pc_y[k] * gsg_136[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pa_x, pc_x, pc_y, fsh0_194, fsg_140, fsh1_194, \
                         gsf0_91, gsf1_91, gsg_137, gsg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * gsg_137[k];

        t_194[k] = pa_x[k] * fsh0_194[k]
                   + f_11 * fsg_140[k]
                   - f_8 * pc_x[k] * fsh1_194[k];

        t_195[k] = f_6 * gsf0_91[k]
                   - f_7 * gsf1_91[k]
                   + f_3 * pc_y[k] * gsg_138[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_x, pc_x, pc_y, fsh0_198, fsg_144, \
                         fsg_145, fsh1_198, gsf0_92, gsf1_92, gsg_139, gsg_140, \
                         gsg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_4 * gsf0_92[k]
                   - f_5 * gsf1_92[k]
                   + f_3 * pc_y[k] * gsg_139[k];

        t_197[k] = f_3 * pc_y[k] * gsg_140[k];

        t_198[k] = pa_x[k] * fsh0_198[k]
                   + f_10 * fsg_144[k]
                   - f_8 * pc_x[k] * fsh1_198[k];

        t_199[k] = f_9 * fsg_145[k]
                   + f_3 * pc_x[k] * gsg_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, fsg_146, fsg_147, fsg_149, \
                         gsg_144, gsg_146, gsg_147, gsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_9 * fsg_146[k]
                   + f_3 * pc_x[k] * gsg_146[k];

        t_201[k] = f_9 * fsg_147[k]
                   + f_3 * pc_x[k] * gsg_147[k];

        t_202[k] = f_3 * pc_y[k] * gsg_144[k];

        t_203[k] = f_9 * fsg_149[k]
                   + f_3 * pc_x[k] * gsg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_x, pc_x, fsh0_204, fsh0_205, fsh0_206, \
                         fsh0_207, fsh1_204, fsh1_205, fsh1_206, \
                         fsh1_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_x[k] * fsh0_204[k]
                   - f_8 * pc_x[k] * fsh1_204[k];

        t_205[k] = pa_x[k] * fsh0_205[k]
                   - f_8 * pc_x[k] * fsh1_205[k];

        t_206[k] = pa_x[k] * fsh0_206[k]
                   - f_8 * pc_x[k] * fsh1_206[k];

        t_207[k] = pa_x[k] * fsh0_207[k]
                   - f_8 * pc_x[k] * fsh1_207[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pc_x, pc_y, fsh0_209, fsh1_209, \
                         gsf0_100, gsf0_101, gsf1_100, gsf1_101, gsg_149, gsg_150, \
                         gsg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * gsg_149[k];

        t_209[k] = pa_x[k] * fsh0_209[k]
                   - f_8 * pc_x[k] * fsh1_209[k];

        t_210[k] = f_1 * gsf0_100[k]
                   - f_2 * gsf1_100[k]
                   + f_3 * pc_x[k] * gsg_150[k];

        t_211[k] = f_12 * gsf0_101[k]
                   - f_13 * gsf1_101[k]
                   + f_3 * pc_x[k] * gsg_151[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_z, gsf0_103, gsf0_105, gsf1_103, \
                         gsf1_105, gsg_150, gsg_151, gsg_153, gsg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * gsg_150[k];

        t_213[k] = f_6 * gsf0_103[k]
                   - f_7 * gsf1_103[k]
                   + f_3 * pc_x[k] * gsg_153[k];

        t_214[k] = f_3 * pc_z[k] * gsg_151[k];

        t_215[k] = f_6 * gsf0_105[k]
                   - f_7 * gsf1_105[k]
                   + f_3 * pc_x[k] * gsg_155[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_z, gsf0_106, gsf0_108, gsf0_109, \
                         gsf1_106, gsf1_108, gsf1_109, gsg_153, gsg_156, gsg_158, \
                         gsg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_4 * gsf0_106[k]
                   - f_5 * gsf1_106[k]
                   + f_3 * pc_x[k] * gsg_156[k];

        t_217[k] = f_3 * pc_z[k] * gsg_153[k];

        t_218[k] = f_4 * gsf0_108[k]
                   - f_5 * gsf1_108[k]
                   + f_3 * pc_x[k] * gsg_158[k];

        t_219[k] = f_4 * gsf0_109[k]
                   - f_5 * gsf1_109[k]
                   + f_3 * pc_x[k] * gsg_159[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, t_225, pc_x, pc_y, fsg_100, \
                         gsf0_106, gsf1_106, gsg_160, gsg_161, gsg_162, gsg_163, \
                         gsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_3 * pc_x[k] * gsg_160[k];

        t_221[k] = f_3 * pc_x[k] * gsg_161[k];

        t_222[k] = f_3 * pc_x[k] * gsg_162[k];

        t_223[k] = f_3 * pc_x[k] * gsg_163[k];

        t_224[k] = f_3 * pc_x[k] * gsg_164[k];

        t_225[k] = f_0 * fsg_100[k]
                   + f_1 * gsf0_106[k]
                   - f_2 * gsf1_106[k]
                   + f_3 * pc_y[k] * gsg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, fsg_104, gsf0_106, gsf0_107, \
                         gsf1_106, gsf1_107, gsg_160, gsg_161, gsg_162, \
                         gsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * gsg_160[k];

        t_227[k] = f_4 * gsf0_106[k]
                   - f_5 * gsf1_106[k]
                   + f_3 * pc_z[k] * gsg_161[k];

        t_228[k] = f_6 * gsf0_107[k]
                   - f_7 * gsf1_107[k]
                   + f_3 * pc_z[k] * gsg_162[k];

        t_229[k] = f_0 * fsg_104[k]
                   + f_3 * pc_y[k] * gsg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pa_z, pc_z, fsh0_126, fsh0_127, fsh1_126, \
                         fsh1_127, gsf0_109, gsf1_109, gsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * gsf0_109[k]
                   - f_2 * gsf1_109[k]
                   + f_3 * pc_z[k] * gsg_164[k];

        t_231[k] = pa_z[k] * fsh0_126[k]
                   - f_8 * pc_z[k] * fsh1_126[k];

        t_232[k] = pa_z[k] * fsh0_127[k]
                   - f_8 * pc_z[k] * fsh1_127[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_z, pc_x, pc_z, fsh0_129, fsh1_129, gsf0_112, \
                         gsf0_114, gsf1_112, gsf1_114, gsg_167, \
                         gsg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_12 * gsf0_112[k]
                   - f_13 * gsf1_112[k]
                   + f_3 * pc_x[k] * gsg_167[k];

        t_234[k] = pa_z[k] * fsh0_129[k]
                   - f_8 * pc_z[k] * fsh1_129[k];

        t_235[k] = f_6 * gsf0_114[k]
                   - f_7 * gsf1_114[k]
                   + f_3 * pc_x[k] * gsg_169[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_z, pc_x, pc_z, fsh0_132, fsh1_132, gsf0_115, \
                         gsf0_117, gsf1_115, gsf1_117, gsg_170, \
                         gsg_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_6 * gsf0_115[k]
                   - f_7 * gsf1_115[k]
                   + f_3 * pc_x[k] * gsg_170[k];

        t_237[k] = pa_z[k] * fsh0_132[k]
                   - f_8 * pc_z[k] * fsh1_132[k];

        t_238[k] = f_4 * gsf0_117[k]
                   - f_5 * gsf1_117[k]
                   + f_3 * pc_x[k] * gsg_172[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pc_x, gsf0_118, gsf0_119, \
                         gsf1_118, gsf1_119, gsg_173, gsg_174, gsg_175, gsg_176, \
                         gsg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_4 * gsf0_118[k]
                   - f_5 * gsf1_118[k]
                   + f_3 * pc_x[k] * gsg_173[k];

        t_240[k] = f_4 * gsf0_119[k]
                   - f_5 * gsf1_119[k]
                   + f_3 * pc_x[k] * gsg_174[k];

        t_241[k] = f_3 * pc_x[k] * gsg_175[k];

        t_242[k] = f_3 * pc_x[k] * gsg_176[k];

        t_243[k] = f_3 * pc_x[k] * gsg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pc_x, pc_z, fsh0_141, fsg_100, \
                         fsh1_141, gsg_175, gsg_178, gsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_3 * pc_x[k] * gsg_178[k];

        t_245[k] = f_3 * pc_x[k] * gsg_179[k];

        t_246[k] = pa_z[k] * fsh0_141[k]
                   - f_8 * pc_z[k] * fsh1_141[k];

        t_247[k] = f_9 * fsg_100[k]
                   + f_3 * pc_z[k] * gsg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pa_z, pc_y, pc_z, fsh0_143, fsh0_144, fsg_101, \
                         fsg_102, fsg_119, fsh1_143, fsh1_144, \
                         gsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pa_z[k] * fsh0_143[k]
                   + f_10 * fsg_101[k]
                   - f_8 * pc_z[k] * fsh1_143[k];

        t_249[k] = pa_z[k] * fsh0_144[k]
                   + f_11 * fsg_102[k]
                   - f_8 * pc_z[k] * fsh1_144[k];

        t_250[k] = f_11 * fsg_119[k]
                   + f_3 * pc_y[k] * gsg_179[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_z, fsg_104, gsf0_119, gsf0_120, \
                         gsf0_121, gsf1_119, gsf1_120, gsf1_121, gsg_179, gsg_180, \
                         gsg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * fsg_104[k]
                   + f_1 * gsf0_119[k]
                   - f_2 * gsf1_119[k]
                   + f_3 * pc_z[k] * gsg_179[k];

        t_252[k] = f_1 * gsf0_120[k]
                   - f_2 * gsf1_120[k]
                   + f_3 * pc_x[k] * gsg_180[k];

        t_253[k] = f_12 * gsf0_121[k]
                   - f_13 * gsf1_121[k]
                   + f_3 * pc_x[k] * gsg_181[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, gsf0_122, gsf0_123, gsf0_124, gsf1_122, \
                         gsf1_123, gsf1_124, gsg_182, gsg_183, \
                         gsg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_12 * gsf0_122[k]
                   - f_13 * gsf1_122[k]
                   + f_3 * pc_x[k] * gsg_182[k];

        t_255[k] = f_6 * gsf0_123[k]
                   - f_7 * gsf1_123[k]
                   + f_3 * pc_x[k] * gsg_183[k];

        t_256[k] = f_6 * gsf0_124[k]
                   - f_7 * gsf1_124[k]
                   + f_3 * pc_x[k] * gsg_184[k];
    }
}

static auto
compute_prim_gsh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsh0,
                                                          const size_t fsg, const size_t fsh1,
                                                          const size_t gsf0, const size_t gsf1,
                                                          const size_t gsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_14 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsh0_189 = buffer.data(fsh0 + 189);
    const auto *fsh0_191 = buffer.data(fsh0 + 191);
    const auto *fsh0_194 = buffer.data(fsh0 + 194);
    const auto *fsh0_198 = buffer.data(fsh0 + 198);
    const auto *fsh0_204 = buffer.data(fsh0 + 204);
    const auto *fsh0_206 = buffer.data(fsh0 + 206);
    const auto *fsh0_207 = buffer.data(fsh0 + 207);
    const auto *fsh0_209 = buffer.data(fsh0 + 209);

    const auto *fsg_115 = buffer.data(fsg + 115);
    const auto *fsg_119 = buffer.data(fsg + 119);
    const auto *fsg_130 = buffer.data(fsg + 130);
    const auto *fsg_132 = buffer.data(fsg + 132);
    const auto *fsg_133 = buffer.data(fsg + 133);
    const auto *fsg_134 = buffer.data(fsg + 134);
    const auto *fsg_145 = buffer.data(fsg + 145);
    const auto *fsg_147 = buffer.data(fsg + 147);
    const auto *fsg_148 = buffer.data(fsg + 148);
    const auto *fsg_149 = buffer.data(fsg + 149);

    const auto *fsh1_189 = buffer.data(fsh1 + 189);
    const auto *fsh1_191 = buffer.data(fsh1 + 191);
    const auto *fsh1_194 = buffer.data(fsh1 + 194);
    const auto *fsh1_198 = buffer.data(fsh1 + 198);
    const auto *fsh1_204 = buffer.data(fsh1 + 204);
    const auto *fsh1_206 = buffer.data(fsh1 + 206);
    const auto *fsh1_207 = buffer.data(fsh1 + 207);
    const auto *fsh1_209 = buffer.data(fsh1 + 209);

    const auto *gsf0_125 = buffer.data(gsf0 + 125);
    const auto *gsf0_126 = buffer.data(gsf0 + 126);
    const auto *gsf0_127 = buffer.data(gsf0 + 127);
    const auto *gsf0_128 = buffer.data(gsf0 + 128);
    const auto *gsf0_129 = buffer.data(gsf0 + 129);
    const auto *gsf0_131 = buffer.data(gsf0 + 131);
    const auto *gsf0_133 = buffer.data(gsf0 + 133);
    const auto *gsf0_134 = buffer.data(gsf0 + 134);
    const auto *gsf0_136 = buffer.data(gsf0 + 136);
    const auto *gsf0_137 = buffer.data(gsf0 + 137);
    const auto *gsf0_138 = buffer.data(gsf0 + 138);
    const auto *gsf0_140 = buffer.data(gsf0 + 140);
    const auto *gsf0_142 = buffer.data(gsf0 + 142);
    const auto *gsf0_143 = buffer.data(gsf0 + 143);
    const auto *gsf0_145 = buffer.data(gsf0 + 145);
    const auto *gsf0_146 = buffer.data(gsf0 + 146);
    const auto *gsf0_147 = buffer.data(gsf0 + 147);
    const auto *gsf0_148 = buffer.data(gsf0 + 148);
    const auto *gsf0_149 = buffer.data(gsf0 + 149);

    const auto *gsf1_125 = buffer.data(gsf1 + 125);
    const auto *gsf1_126 = buffer.data(gsf1 + 126);
    const auto *gsf1_127 = buffer.data(gsf1 + 127);
    const auto *gsf1_128 = buffer.data(gsf1 + 128);
    const auto *gsf1_129 = buffer.data(gsf1 + 129);
    const auto *gsf1_131 = buffer.data(gsf1 + 131);
    const auto *gsf1_133 = buffer.data(gsf1 + 133);
    const auto *gsf1_134 = buffer.data(gsf1 + 134);
    const auto *gsf1_136 = buffer.data(gsf1 + 136);
    const auto *gsf1_137 = buffer.data(gsf1 + 137);
    const auto *gsf1_138 = buffer.data(gsf1 + 138);
    const auto *gsf1_140 = buffer.data(gsf1 + 140);
    const auto *gsf1_142 = buffer.data(gsf1 + 142);
    const auto *gsf1_143 = buffer.data(gsf1 + 143);
    const auto *gsf1_145 = buffer.data(gsf1 + 145);
    const auto *gsf1_146 = buffer.data(gsf1 + 146);
    const auto *gsf1_147 = buffer.data(gsf1 + 147);
    const auto *gsf1_148 = buffer.data(gsf1 + 148);
    const auto *gsf1_149 = buffer.data(gsf1 + 149);

    const auto *gsg_185 = buffer.data(gsg + 185);
    const auto *gsg_186 = buffer.data(gsg + 186);
    const auto *gsg_187 = buffer.data(gsg + 187);
    const auto *gsg_188 = buffer.data(gsg + 188);
    const auto *gsg_189 = buffer.data(gsg + 189);
    const auto *gsg_190 = buffer.data(gsg + 190);
    const auto *gsg_191 = buffer.data(gsg + 191);
    const auto *gsg_192 = buffer.data(gsg + 192);
    const auto *gsg_193 = buffer.data(gsg + 193);
    const auto *gsg_194 = buffer.data(gsg + 194);
    const auto *gsg_196 = buffer.data(gsg + 196);
    const auto *gsg_198 = buffer.data(gsg + 198);
    const auto *gsg_199 = buffer.data(gsg + 199);
    const auto *gsg_201 = buffer.data(gsg + 201);
    const auto *gsg_202 = buffer.data(gsg + 202);
    const auto *gsg_203 = buffer.data(gsg + 203);
    const auto *gsg_205 = buffer.data(gsg + 205);
    const auto *gsg_206 = buffer.data(gsg + 206);
    const auto *gsg_207 = buffer.data(gsg + 207);
    const auto *gsg_208 = buffer.data(gsg + 208);
    const auto *gsg_209 = buffer.data(gsg + 209);
    const auto *gsg_210 = buffer.data(gsg + 210);
    const auto *gsg_212 = buffer.data(gsg + 212);
    const auto *gsg_213 = buffer.data(gsg + 213);
    const auto *gsg_215 = buffer.data(gsg + 215);
    const auto *gsg_216 = buffer.data(gsg + 216);
    const auto *gsg_217 = buffer.data(gsg + 217);
    const auto *gsg_219 = buffer.data(gsg + 219);
    const auto *gsg_220 = buffer.data(gsg + 220);
    const auto *gsg_221 = buffer.data(gsg + 221);
    const auto *gsg_222 = buffer.data(gsg + 222);
    const auto *gsg_223 = buffer.data(gsg + 223);
    const auto *gsg_224 = buffer.data(gsg + 224);

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, gsf0_125, gsf0_126, gsf0_127, gsf1_125, \
                         gsf1_126, gsf1_127, gsg_185, gsg_186, \
                         gsg_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_6 * gsf0_125[k]
                   - f_7 * gsf1_125[k]
                   + f_3 * pc_x[k] * gsg_185[k];

        t_258[k] = f_4 * gsf0_126[k]
                   - f_5 * gsf1_126[k]
                   + f_3 * pc_x[k] * gsg_186[k];

        t_259[k] = f_4 * gsf0_127[k]
                   - f_5 * gsf1_127[k]
                   + f_3 * pc_x[k] * gsg_187[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, pc_x, gsf0_128, gsf0_129, \
                         gsf1_128, gsf1_129, gsg_188, gsg_189, gsg_190, gsg_191, \
                         gsg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_4 * gsf0_128[k]
                   - f_5 * gsf1_128[k]
                   + f_3 * pc_x[k] * gsg_188[k];

        t_261[k] = f_4 * gsf0_129[k]
                   - f_5 * gsf1_129[k]
                   + f_3 * pc_x[k] * gsg_189[k];

        t_262[k] = f_3 * pc_x[k] * gsg_190[k];

        t_263[k] = f_3 * pc_x[k] * gsg_191[k];

        t_264[k] = f_3 * pc_x[k] * gsg_192[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pc_x, pc_y, pc_z, fsg_115, fsg_130, \
                         gsf0_126, gsf1_126, gsg_190, gsg_193, \
                         gsg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_3 * pc_x[k] * gsg_193[k];

        t_266[k] = f_3 * pc_x[k] * gsg_194[k];

        t_267[k] = f_10 * fsg_130[k]
                   + f_1 * gsf0_126[k]
                   - f_2 * gsf1_126[k]
                   + f_3 * pc_y[k] * gsg_190[k];

        t_268[k] = f_10 * fsg_115[k]
                   + f_3 * pc_z[k] * gsg_190[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pc_y, fsg_132, fsg_133, fsg_134, gsf0_128, \
                         gsf0_129, gsf1_128, gsf1_129, gsg_192, gsg_193, \
                         gsg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_10 * fsg_132[k]
                   + f_6 * gsf0_128[k]
                   - f_7 * gsf1_128[k]
                   + f_3 * pc_y[k] * gsg_192[k];

        t_270[k] = f_10 * fsg_133[k]
                   + f_4 * gsf0_129[k]
                   - f_5 * gsf1_129[k]
                   + f_3 * pc_y[k] * gsg_193[k];

        t_271[k] = f_10 * fsg_134[k]
                   + f_3 * pc_y[k] * gsg_194[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pa_y, pc_x, pc_y, pc_z, fsh0_189, fsg_119, \
                         fsh1_189, gsf0_129, gsf0_131, gsf1_129, gsf1_131, gsg_194, \
                         gsg_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_10 * fsg_119[k]
                   + f_1 * gsf0_129[k]
                   - f_2 * gsf1_129[k]
                   + f_3 * pc_z[k] * gsg_194[k];

        t_273[k] = pa_y[k] * fsh0_189[k]
                   - f_8 * pc_y[k] * fsh1_189[k];

        t_274[k] = f_12 * gsf0_131[k]
                   - f_13 * gsf1_131[k]
                   + f_3 * pc_x[k] * gsg_196[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pa_y, pc_x, pc_y, fsh0_191, fsh1_191, gsf0_133, \
                         gsf0_134, gsf1_133, gsf1_134, gsg_198, \
                         gsg_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = pa_y[k] * fsh0_191[k]
                   - f_8 * pc_y[k] * fsh1_191[k];

        t_276[k] = f_6 * gsf0_133[k]
                   - f_7 * gsf1_133[k]
                   + f_3 * pc_x[k] * gsg_198[k];

        t_277[k] = f_6 * gsf0_134[k]
                   - f_7 * gsf1_134[k]
                   + f_3 * pc_x[k] * gsg_199[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pa_y, pc_x, pc_y, fsh0_194, fsh1_194, gsf0_136, \
                         gsf0_137, gsf1_136, gsf1_137, gsg_201, \
                         gsg_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = pa_y[k] * fsh0_194[k]
                   - f_8 * pc_y[k] * fsh1_194[k];

        t_279[k] = f_4 * gsf0_136[k]
                   - f_5 * gsf1_136[k]
                   + f_3 * pc_x[k] * gsg_201[k];

        t_280[k] = f_4 * gsf0_137[k]
                   - f_5 * gsf1_137[k]
                   + f_3 * pc_x[k] * gsg_202[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, pa_y, pc_x, pc_y, fsh0_198, \
                         fsh1_198, gsf0_138, gsf1_138, gsg_203, gsg_205, gsg_206, \
                         gsg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_4 * gsf0_138[k]
                   - f_5 * gsf1_138[k]
                   + f_3 * pc_x[k] * gsg_203[k];

        t_282[k] = pa_y[k] * fsh0_198[k]
                   - f_8 * pc_y[k] * fsh1_198[k];

        t_283[k] = f_3 * pc_x[k] * gsg_205[k];

        t_284[k] = f_3 * pc_x[k] * gsg_206[k];

        t_285[k] = f_3 * pc_x[k] * gsg_207[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pc_x, pc_y, pc_z, fsh0_204, \
                         fsg_130, fsg_145, fsh1_204, gsg_205, gsg_208, \
                         gsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_3 * pc_x[k] * gsg_208[k];

        t_287[k] = f_3 * pc_x[k] * gsg_209[k];

        t_288[k] = pa_y[k] * fsh0_204[k]
                   + f_14 * fsg_145[k]
                   - f_8 * pc_y[k] * fsh1_204[k];

        t_289[k] = f_11 * fsg_130[k]
                   + f_3 * pc_z[k] * gsg_205[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pc_y, fsh0_206, fsh0_207, fsh0_209, \
                         fsg_147, fsg_148, fsg_149, fsh1_206, fsh1_207, fsh1_209, \
                         gsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_y[k] * fsh0_206[k]
                   + f_11 * fsg_147[k]
                   - f_8 * pc_y[k] * fsh1_206[k];

        t_291[k] = pa_y[k] * fsh0_207[k]
                   + f_10 * fsg_148[k]
                   - f_8 * pc_y[k] * fsh1_207[k];

        t_292[k] = f_9 * fsg_149[k]
                   + f_3 * pc_y[k] * gsg_209[k];

        t_293[k] = pa_y[k] * fsh0_209[k]
                   - f_8 * pc_y[k] * fsh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, gsf0_140, gsf0_142, \
                         gsf0_143, gsf1_140, gsf1_142, gsf1_143, gsg_210, gsg_212, \
                         gsg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_1 * gsf0_140[k]
                   - f_2 * gsf1_140[k]
                   + f_3 * pc_x[k] * gsg_210[k];

        t_295[k] = f_3 * pc_y[k] * gsg_210[k];

        t_296[k] = f_12 * gsf0_142[k]
                   - f_13 * gsf1_142[k]
                   + f_3 * pc_x[k] * gsg_212[k];

        t_297[k] = f_6 * gsf0_143[k]
                   - f_7 * gsf1_143[k]
                   + f_3 * pc_x[k] * gsg_213[k];

        t_298[k] = f_3 * pc_y[k] * gsg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, gsf0_145, gsf0_146, gsf0_147, \
                         gsf1_145, gsf1_146, gsf1_147, gsg_215, gsg_216, \
                         gsg_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_6 * gsf0_145[k]
                   - f_7 * gsf1_145[k]
                   + f_3 * pc_x[k] * gsg_215[k];

        t_300[k] = f_4 * gsf0_146[k]
                   - f_5 * gsf1_146[k]
                   + f_3 * pc_x[k] * gsg_216[k];

        t_301[k] = f_4 * gsf0_147[k]
                   - f_5 * gsf1_147[k]
                   + f_3 * pc_x[k] * gsg_217[k];

        t_302[k] = f_3 * pc_y[k] * gsg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, t_308, pc_x, gsf0_149, gsf1_149, \
                         gsg_219, gsg_220, gsg_221, gsg_222, gsg_223, \
                         gsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_4 * gsf0_149[k]
                   - f_5 * gsf1_149[k]
                   + f_3 * pc_x[k] * gsg_219[k];

        t_304[k] = f_3 * pc_x[k] * gsg_220[k];

        t_305[k] = f_3 * pc_x[k] * gsg_221[k];

        t_306[k] = f_3 * pc_x[k] * gsg_222[k];

        t_307[k] = f_3 * pc_x[k] * gsg_223[k];

        t_308[k] = f_3 * pc_x[k] * gsg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_y, gsf0_146, gsf0_147, gsf0_148, gsf1_146, \
                         gsf1_147, gsf1_148, gsg_220, gsg_221, \
                         gsg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * gsf0_146[k]
                   - f_2 * gsf1_146[k]
                   + f_3 * pc_y[k] * gsg_220[k];

        t_310[k] = f_12 * gsf0_147[k]
                   - f_13 * gsf1_147[k]
                   + f_3 * pc_y[k] * gsg_221[k];

        t_311[k] = f_6 * gsf0_148[k]
                   - f_7 * gsf1_148[k]
                   + f_3 * pc_y[k] * gsg_222[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pc_y, pc_z, fsg_149, gsf0_149, gsf1_149, \
                         gsg_223, gsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * gsf0_149[k]
                   - f_5 * gsf1_149[k]
                   + f_3 * pc_y[k] * gsg_223[k];

        t_313[k] = f_3 * pc_y[k] * gsg_224[k];

        t_314[k] = f_0 * fsg_149[k]
                   + f_1 * gsf0_149[k]
                   - f_2 * gsf1_149[k]
                   + f_3 * pc_z[k] * gsg_224[k];
    }
}

auto
compute_prim_gsh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fsh0, const size_t fsg,
                                                   const size_t fsh1, const size_t gsf0,
                                                   const size_t gsf1, const size_t gsg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gsh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, fsh0, fsg,
                                                              fsh1, gsf0, gsf1, gsg, ncols,
                                                              gamma, p, q);

    compute_prim_gsh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, fsh0, fsg,
                                                              fsh1, gsf0, gsf1, gsg, ncols,
                                                              gamma, p, q);

    compute_prim_gsh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, fsh0, fsg,
                                                              fsh1, gsf0, gsf1, gsg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
