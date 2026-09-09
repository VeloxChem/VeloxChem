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


#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_12 = 4.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_0 = buffer.data(lsh0 + 0);
    const auto *lsh0_3 = buffer.data(lsh0 + 3);
    const auto *lsh0_5 = buffer.data(lsh0 + 5);
    const auto *lsh0_6 = buffer.data(lsh0 + 6);
    const auto *lsh0_9 = buffer.data(lsh0 + 9);
    const auto *lsh0_15 = buffer.data(lsh0 + 15);
    const auto *lsh0_20 = buffer.data(lsh0 + 20);
    const auto *lsh0_24 = buffer.data(lsh0 + 24);
    const auto *lsh0_27 = buffer.data(lsh0 + 27);
    const auto *lsh0_36 = buffer.data(lsh0 + 36);
    const auto *lsh0_42 = buffer.data(lsh0 + 42);
    const auto *lsh0_47 = buffer.data(lsh0 + 47);
    const auto *lsh0_51 = buffer.data(lsh0 + 51);
    const auto *lsh0_62 = buffer.data(lsh0 + 62);

    const auto *lsg_0 = buffer.data(lsg + 0);
    const auto *lsg_1 = buffer.data(lsg + 1);
    const auto *lsg_2 = buffer.data(lsg + 2);
    const auto *lsg_3 = buffer.data(lsg + 3);
    const auto *lsg_5 = buffer.data(lsg + 5);
    const auto *lsg_10 = buffer.data(lsg + 10);
    const auto *lsg_12 = buffer.data(lsg + 12);
    const auto *lsg_14 = buffer.data(lsg + 14);
    const auto *lsg_15 = buffer.data(lsg + 15);
    const auto *lsg_18 = buffer.data(lsg + 18);
    const auto *lsg_20 = buffer.data(lsg + 20);
    const auto *lsg_25 = buffer.data(lsg + 25);
    const auto *lsg_27 = buffer.data(lsg + 27);
    const auto *lsg_28 = buffer.data(lsg + 28);
    const auto *lsg_29 = buffer.data(lsg + 29);
    const auto *lsg_30 = buffer.data(lsg + 30);
    const auto *lsg_32 = buffer.data(lsg + 32);
    const auto *lsg_35 = buffer.data(lsg + 35);
    const auto *lsg_40 = buffer.data(lsg + 40);
    const auto *lsg_41 = buffer.data(lsg + 41);
    const auto *lsg_42 = buffer.data(lsg + 42);
    const auto *lsg_43 = buffer.data(lsg + 43);
    const auto *lsg_44 = buffer.data(lsg + 44);
    const auto *lsg_45 = buffer.data(lsg + 45);
    const auto *lsg_48 = buffer.data(lsg + 48);
    const auto *lsg_51 = buffer.data(lsg + 51);
    const auto *lsg_55 = buffer.data(lsg + 55);
    const auto *lsg_57 = buffer.data(lsg + 57);
    const auto *lsg_58 = buffer.data(lsg + 58);
    const auto *lsg_59 = buffer.data(lsg + 59);
    const auto *lsg_70 = buffer.data(lsg + 70);
    const auto *lsg_71 = buffer.data(lsg + 71);
    const auto *lsg_72 = buffer.data(lsg + 72);
    const auto *lsg_73 = buffer.data(lsg + 73);
    const auto *lsg_74 = buffer.data(lsg + 74);
    const auto *lsg_75 = buffer.data(lsg + 75);
    const auto *lsg_80 = buffer.data(lsg + 80);
    const auto *lsg_84 = buffer.data(lsg + 84);
    const auto *lsg_85 = buffer.data(lsg + 85);
    const auto *lsg_86 = buffer.data(lsg + 86);
    const auto *lsg_87 = buffer.data(lsg + 87);
    const auto *lsg_89 = buffer.data(lsg + 89);
    const auto *lsg_90 = buffer.data(lsg + 90);
    const auto *lsg_93 = buffer.data(lsg + 93);

    const auto *lsh1_0 = buffer.data(lsh1 + 0);
    const auto *lsh1_3 = buffer.data(lsh1 + 3);
    const auto *lsh1_5 = buffer.data(lsh1 + 5);
    const auto *lsh1_6 = buffer.data(lsh1 + 6);
    const auto *lsh1_9 = buffer.data(lsh1 + 9);
    const auto *lsh1_15 = buffer.data(lsh1 + 15);
    const auto *lsh1_20 = buffer.data(lsh1 + 20);
    const auto *lsh1_24 = buffer.data(lsh1 + 24);
    const auto *lsh1_27 = buffer.data(lsh1 + 27);
    const auto *lsh1_36 = buffer.data(lsh1 + 36);
    const auto *lsh1_42 = buffer.data(lsh1 + 42);
    const auto *lsh1_47 = buffer.data(lsh1 + 47);
    const auto *lsh1_51 = buffer.data(lsh1 + 51);
    const auto *lsh1_62 = buffer.data(lsh1 + 62);

    const auto *msf0_0 = buffer.data(msf0 + 0);
    const auto *msf0_1 = buffer.data(msf0 + 1);
    const auto *msf0_2 = buffer.data(msf0 + 2);
    const auto *msf0_6 = buffer.data(msf0 + 6);
    const auto *msf0_8 = buffer.data(msf0 + 8);
    const auto *msf0_9 = buffer.data(msf0 + 9);
    const auto *msf0_16 = buffer.data(msf0 + 16);
    const auto *msf0_17 = buffer.data(msf0 + 17);
    const auto *msf0_22 = buffer.data(msf0 + 22);
    const auto *msf0_27 = buffer.data(msf0 + 27);
    const auto *msf0_28 = buffer.data(msf0 + 28);
    const auto *msf0_29 = buffer.data(msf0 + 29);
    const auto *msf0_30 = buffer.data(msf0 + 30);
    const auto *msf0_32 = buffer.data(msf0 + 32);
    const auto *msf0_33 = buffer.data(msf0 + 33);
    const auto *msf0_36 = buffer.data(msf0 + 36);
    const auto *msf0_37 = buffer.data(msf0 + 37);
    const auto *msf0_39 = buffer.data(msf0 + 39);
    const auto *msf0_48 = buffer.data(msf0 + 48);
    const auto *msf0_49 = buffer.data(msf0 + 49);
    const auto *msf0_50 = buffer.data(msf0 + 50);
    const auto *msf0_51 = buffer.data(msf0 + 51);
    const auto *msf0_52 = buffer.data(msf0 + 52);
    const auto *msf0_55 = buffer.data(msf0 + 55);
    const auto *msf0_56 = buffer.data(msf0 + 56);
    const auto *msf0_57 = buffer.data(msf0 + 57);
    const auto *msf0_58 = buffer.data(msf0 + 58);
    const auto *msf0_59 = buffer.data(msf0 + 59);
    const auto *msf0_60 = buffer.data(msf0 + 60);
    const auto *msf0_63 = buffer.data(msf0 + 63);

    const auto *msf1_0 = buffer.data(msf1 + 0);
    const auto *msf1_1 = buffer.data(msf1 + 1);
    const auto *msf1_2 = buffer.data(msf1 + 2);
    const auto *msf1_6 = buffer.data(msf1 + 6);
    const auto *msf1_8 = buffer.data(msf1 + 8);
    const auto *msf1_9 = buffer.data(msf1 + 9);
    const auto *msf1_16 = buffer.data(msf1 + 16);
    const auto *msf1_17 = buffer.data(msf1 + 17);
    const auto *msf1_22 = buffer.data(msf1 + 22);
    const auto *msf1_27 = buffer.data(msf1 + 27);
    const auto *msf1_28 = buffer.data(msf1 + 28);
    const auto *msf1_29 = buffer.data(msf1 + 29);
    const auto *msf1_30 = buffer.data(msf1 + 30);
    const auto *msf1_32 = buffer.data(msf1 + 32);
    const auto *msf1_33 = buffer.data(msf1 + 33);
    const auto *msf1_36 = buffer.data(msf1 + 36);
    const auto *msf1_37 = buffer.data(msf1 + 37);
    const auto *msf1_39 = buffer.data(msf1 + 39);
    const auto *msf1_48 = buffer.data(msf1 + 48);
    const auto *msf1_49 = buffer.data(msf1 + 49);
    const auto *msf1_50 = buffer.data(msf1 + 50);
    const auto *msf1_51 = buffer.data(msf1 + 51);
    const auto *msf1_52 = buffer.data(msf1 + 52);
    const auto *msf1_55 = buffer.data(msf1 + 55);
    const auto *msf1_56 = buffer.data(msf1 + 56);
    const auto *msf1_57 = buffer.data(msf1 + 57);
    const auto *msf1_58 = buffer.data(msf1 + 58);
    const auto *msf1_59 = buffer.data(msf1 + 59);
    const auto *msf1_60 = buffer.data(msf1 + 60);
    const auto *msf1_63 = buffer.data(msf1 + 63);

    const auto *msg_0 = buffer.data(msg + 0);
    const auto *msg_1 = buffer.data(msg + 1);
    const auto *msg_2 = buffer.data(msg + 2);
    const auto *msg_3 = buffer.data(msg + 3);
    const auto *msg_5 = buffer.data(msg + 5);
    const auto *msg_6 = buffer.data(msg + 6);
    const auto *msg_9 = buffer.data(msg + 9);
    const auto *msg_10 = buffer.data(msg + 10);
    const auto *msg_12 = buffer.data(msg + 12);
    const auto *msg_13 = buffer.data(msg + 13);
    const auto *msg_14 = buffer.data(msg + 14);
    const auto *msg_15 = buffer.data(msg + 15);
    const auto *msg_16 = buffer.data(msg + 16);
    const auto *msg_18 = buffer.data(msg + 18);
    const auto *msg_20 = buffer.data(msg + 20);
    const auto *msg_21 = buffer.data(msg + 21);
    const auto *msg_25 = buffer.data(msg + 25);
    const auto *msg_26 = buffer.data(msg + 26);
    const auto *msg_27 = buffer.data(msg + 27);
    const auto *msg_28 = buffer.data(msg + 28);
    const auto *msg_29 = buffer.data(msg + 29);
    const auto *msg_30 = buffer.data(msg + 30);
    const auto *msg_32 = buffer.data(msg + 32);
    const auto *msg_34 = buffer.data(msg + 34);
    const auto *msg_35 = buffer.data(msg + 35);
    const auto *msg_39 = buffer.data(msg + 39);
    const auto *msg_40 = buffer.data(msg + 40);
    const auto *msg_41 = buffer.data(msg + 41);
    const auto *msg_42 = buffer.data(msg + 42);
    const auto *msg_43 = buffer.data(msg + 43);
    const auto *msg_44 = buffer.data(msg + 44);
    const auto *msg_45 = buffer.data(msg + 45);
    const auto *msg_46 = buffer.data(msg + 46);
    const auto *msg_47 = buffer.data(msg + 47);
    const auto *msg_48 = buffer.data(msg + 48);
    const auto *msg_50 = buffer.data(msg + 50);
    const auto *msg_51 = buffer.data(msg + 51);
    const auto *msg_55 = buffer.data(msg + 55);
    const auto *msg_56 = buffer.data(msg + 56);
    const auto *msg_57 = buffer.data(msg + 57);
    const auto *msg_58 = buffer.data(msg + 58);
    const auto *msg_59 = buffer.data(msg + 59);
    const auto *msg_60 = buffer.data(msg + 60);
    const auto *msg_62 = buffer.data(msg + 62);
    const auto *msg_63 = buffer.data(msg + 63);
    const auto *msg_65 = buffer.data(msg + 65);
    const auto *msg_70 = buffer.data(msg + 70);
    const auto *msg_71 = buffer.data(msg + 71);
    const auto *msg_72 = buffer.data(msg + 72);
    const auto *msg_73 = buffer.data(msg + 73);
    const auto *msg_74 = buffer.data(msg + 74);
    const auto *msg_75 = buffer.data(msg + 75);
    const auto *msg_76 = buffer.data(msg + 76);
    const auto *msg_77 = buffer.data(msg + 77);
    const auto *msg_78 = buffer.data(msg + 78);
    const auto *msg_79 = buffer.data(msg + 79);
    const auto *msg_80 = buffer.data(msg + 80);
    const auto *msg_84 = buffer.data(msg + 84);
    const auto *msg_85 = buffer.data(msg + 85);
    const auto *msg_86 = buffer.data(msg + 86);
    const auto *msg_87 = buffer.data(msg + 87);
    const auto *msg_88 = buffer.data(msg + 88);
    const auto *msg_89 = buffer.data(msg + 89);
    const auto *msg_90 = buffer.data(msg + 90);
    const auto *msg_93 = buffer.data(msg + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, lsg_0, msf0_0, \
                         msf1_0, msg_0, msg_1, msg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lsg_0[k]
                 + f_1 * msf0_0[k]
                 - f_2 * msf1_0[k]
                 + f_3 * pc_x[k] * msg_0[k];

        t_1[k] = f_3 * pc_y[k] * msg_0[k];

        t_2[k] = f_3 * pc_z[k] * msg_0[k];

        t_3[k] = f_4 * msf0_0[k]
                 - f_5 * msf1_0[k]
                 + f_3 * pc_y[k] * msg_1[k];

        t_4[k] = f_3 * pc_y[k] * msg_2[k];

        t_5[k] = f_4 * msf0_0[k]
                 - f_5 * msf1_0[k]
                 + f_3 * pc_z[k] * msg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, lsg_10, msf0_1, msf0_2, \
                         msf1_1, msf1_2, msg_3, msg_5, msg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * msf0_1[k]
                 - f_7 * msf1_1[k]
                 + f_3 * pc_y[k] * msg_3[k];

        t_7[k] = f_3 * pc_z[k] * msg_3[k];

        t_8[k] = f_3 * pc_y[k] * msg_5[k];

        t_9[k] = f_6 * msf0_2[k]
                 - f_7 * msf1_2[k]
                 + f_3 * pc_z[k] * msg_5[k];

        t_10[k] = f_0 * lsg_10[k]
                  + f_3 * pc_x[k] * msg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, lsg_12, lsg_14, msg_6, \
                         msg_9, msg_12, msg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * msg_6[k];

        t_12[k] = f_0 * lsg_12[k]
                  + f_3 * pc_x[k] * msg_12[k];

        t_13[k] = f_3 * pc_y[k] * msg_9[k];

        t_14[k] = f_0 * lsg_14[k]
                  + f_3 * pc_x[k] * msg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, msf0_6, msf0_8, msf0_9, msf1_6, \
                         msf1_8, msf1_9, msg_10, msg_12, msg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * msf0_6[k]
                  - f_2 * msf1_6[k]
                  + f_3 * pc_y[k] * msg_10[k];

        t_16[k] = f_3 * pc_z[k] * msg_10[k];

        t_17[k] = f_6 * msf0_8[k]
                  - f_7 * msf1_8[k]
                  + f_3 * pc_y[k] * msg_12[k];

        t_18[k] = f_4 * msf0_9[k]
                  - f_5 * msf1_9[k]
                  + f_3 * pc_y[k] * msg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, lsh0_0, lsg_0, \
                         lsh1_0, msf0_9, msf1_9, msg_14, msg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * msg_14[k];

        t_20[k] = f_1 * msf0_9[k]
                  - f_2 * msf1_9[k]
                  + f_3 * pc_z[k] * msg_14[k];

        t_21[k] = pa_y[k] * lsh0_0[k]
                  - f_8 * pc_y[k] * lsh1_0[k];

        t_22[k] = f_9 * lsg_0[k]
                  + f_3 * pc_y[k] * msg_15[k];

        t_23[k] = f_3 * pc_z[k] * msg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, lsh0_3, lsh0_5, lsh0_6, \
                         lsg_1, lsg_3, lsh1_3, lsh1_5, lsh1_6, msg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * lsh0_3[k]
                  + f_10 * lsg_1[k]
                  - f_8 * pc_y[k] * lsh1_3[k];

        t_25[k] = f_3 * pc_z[k] * msg_16[k];

        t_26[k] = pa_y[k] * lsh0_5[k]
                  - f_8 * pc_y[k] * lsh1_5[k];

        t_27[k] = pa_y[k] * lsh0_6[k]
                  + f_11 * lsg_3[k]
                  - f_8 * pc_y[k] * lsh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, lsh0_9, lsg_5, \
                         lsg_25, lsh1_9, msg_18, msg_20, msg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * msg_18[k];

        t_29[k] = f_9 * lsg_5[k]
                  + f_3 * pc_y[k] * msg_20[k];

        t_30[k] = pa_y[k] * lsh0_9[k]
                  - f_8 * pc_y[k] * lsh1_9[k];

        t_31[k] = f_12 * lsg_25[k]
                  + f_3 * pc_x[k] * msg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, lsg_27, lsg_28, lsg_29, msg_21, \
                         msg_27, msg_28, msg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * msg_21[k];

        t_33[k] = f_12 * lsg_27[k]
                  + f_3 * pc_x[k] * msg_27[k];

        t_34[k] = f_12 * lsg_28[k]
                  + f_3 * pc_x[k] * msg_28[k];

        t_35[k] = f_12 * lsg_29[k]
                  + f_3 * pc_x[k] * msg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, lsg_10, msf0_16, msf0_17, \
                         msf1_16, msf1_17, msg_25, msg_26, msg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * lsg_10[k]
                  + f_1 * msf0_16[k]
                  - f_2 * msf1_16[k]
                  + f_3 * pc_y[k] * msg_25[k];

        t_37[k] = f_3 * pc_z[k] * msg_25[k];

        t_38[k] = f_4 * msf0_16[k]
                  - f_5 * msf1_16[k]
                  + f_3 * pc_z[k] * msg_26[k];

        t_39[k] = f_6 * msf0_17[k]
                  - f_7 * msf1_17[k]
                  + f_3 * pc_z[k] * msg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, lsh0_0, lsh0_20, \
                         lsg_14, lsh1_0, lsh1_20, msg_29, msg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * lsg_14[k]
                  + f_3 * pc_y[k] * msg_29[k];

        t_41[k] = pa_y[k] * lsh0_20[k]
                  - f_8 * pc_y[k] * lsh1_20[k];

        t_42[k] = pa_z[k] * lsh0_0[k]
                  - f_8 * pc_z[k] * lsh1_0[k];

        t_43[k] = f_3 * pc_y[k] * msg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, lsh0_3, lsh0_5, lsg_0, \
                         lsg_2, lsh1_3, lsh1_5, msg_30, msg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * lsg_0[k]
                  + f_3 * pc_z[k] * msg_30[k];

        t_45[k] = pa_z[k] * lsh0_3[k]
                  - f_8 * pc_z[k] * lsh1_3[k];

        t_46[k] = f_3 * pc_y[k] * msg_32[k];

        t_47[k] = pa_z[k] * lsh0_5[k]
                  + f_10 * lsg_2[k]
                  - f_8 * pc_z[k] * lsh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, lsh0_6, lsh0_9, lsg_5, \
                         lsh1_6, lsh1_9, msf0_22, msf1_22, msg_34, \
                         msg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * lsh0_6[k]
                  - f_8 * pc_z[k] * lsh1_6[k];

        t_49[k] = f_4 * msf0_22[k]
                  - f_5 * msf1_22[k]
                  + f_3 * pc_y[k] * msg_34[k];

        t_50[k] = f_3 * pc_y[k] * msg_35[k];

        t_51[k] = pa_z[k] * lsh0_9[k]
                  + f_11 * lsg_5[k]
                  - f_8 * pc_z[k] * lsh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, lsg_40, lsg_41, lsg_42, \
                         lsg_44, msg_39, msg_40, msg_41, msg_42, \
                         msg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * lsg_40[k]
                  + f_3 * pc_x[k] * msg_40[k];

        t_53[k] = f_12 * lsg_41[k]
                  + f_3 * pc_x[k] * msg_41[k];

        t_54[k] = f_12 * lsg_42[k]
                  + f_3 * pc_x[k] * msg_42[k];

        t_55[k] = f_3 * pc_y[k] * msg_39[k];

        t_56[k] = f_12 * lsg_44[k]
                  + f_3 * pc_x[k] * msg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, lsh0_15, lsh1_15, msf0_27, \
                         msf0_28, msf1_27, msf1_28, msg_41, msg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * lsh0_15[k]
                  - f_8 * pc_z[k] * lsh1_15[k];

        t_58[k] = f_13 * msf0_27[k]
                  - f_14 * msf1_27[k]
                  + f_3 * pc_y[k] * msg_41[k];

        t_59[k] = f_6 * msf0_28[k]
                  - f_7 * msf1_28[k]
                  + f_3 * pc_y[k] * msg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, lsg_14, lsg_45, msf0_29, \
                         msf0_30, msf1_29, msf1_30, msg_43, msg_44, \
                         msg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * msf0_29[k]
                  - f_5 * msf1_29[k]
                  + f_3 * pc_y[k] * msg_43[k];

        t_61[k] = f_3 * pc_y[k] * msg_44[k];

        t_62[k] = f_9 * lsg_14[k]
                  + f_1 * msf0_29[k]
                  - f_2 * msf1_29[k]
                  + f_3 * pc_z[k] * msg_44[k];

        t_63[k] = f_15 * lsg_45[k]
                  + f_1 * msf0_30[k]
                  - f_2 * msf1_30[k]
                  + f_3 * pc_x[k] * msg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, lsg_15, lsg_48, msf0_33, \
                         msf1_33, msg_45, msg_46, msg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * lsg_15[k]
                  + f_3 * pc_y[k] * msg_45[k];

        t_65[k] = f_3 * pc_z[k] * msg_45[k];

        t_66[k] = f_15 * lsg_48[k]
                  + f_6 * msf0_33[k]
                  - f_7 * msf1_33[k]
                  + f_3 * pc_x[k] * msg_48[k];

        t_67[k] = f_3 * pc_z[k] * msg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, lsg_51, msf0_30, msf0_36, msf1_30, \
                         msf1_36, msg_47, msg_48, msg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * msf0_30[k]
                  - f_5 * msf1_30[k]
                  + f_3 * pc_z[k] * msg_47[k];

        t_69[k] = f_15 * lsg_51[k]
                  + f_4 * msf0_36[k]
                  - f_5 * msf1_36[k]
                  + f_3 * pc_x[k] * msg_51[k];

        t_70[k] = f_3 * pc_z[k] * msg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, lsg_20, lsg_55, msf0_32, \
                         msf1_32, msg_50, msg_51, msg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * lsg_20[k]
                  + f_3 * pc_y[k] * msg_50[k];

        t_72[k] = f_6 * msf0_32[k]
                  - f_7 * msf1_32[k]
                  + f_3 * pc_z[k] * msg_50[k];

        t_73[k] = f_15 * lsg_55[k]
                  + f_3 * pc_x[k] * msg_55[k];

        t_74[k] = f_3 * pc_z[k] * msg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, lsg_25, lsg_57, lsg_58, lsg_59, \
                         msf0_36, msf1_36, msg_55, msg_57, msg_58, \
                         msg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_15 * lsg_57[k]
                  + f_3 * pc_x[k] * msg_57[k];

        t_76[k] = f_15 * lsg_58[k]
                  + f_3 * pc_x[k] * msg_58[k];

        t_77[k] = f_15 * lsg_59[k]
                  + f_3 * pc_x[k] * msg_59[k];

        t_78[k] = f_10 * lsg_25[k]
                  + f_1 * msf0_36[k]
                  - f_2 * msf1_36[k]
                  + f_3 * pc_y[k] * msg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, lsg_29, msf0_36, msf0_37, \
                         msf1_36, msf1_37, msg_55, msg_56, msg_57, \
                         msg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * msg_55[k];

        t_80[k] = f_4 * msf0_36[k]
                  - f_5 * msf1_36[k]
                  + f_3 * pc_z[k] * msg_56[k];

        t_81[k] = f_6 * msf0_37[k]
                  - f_7 * msf1_37[k]
                  + f_3 * pc_z[k] * msg_57[k];

        t_82[k] = f_10 * lsg_29[k]
                  + f_3 * pc_y[k] * msg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, lsh0_42, lsg_15, lsg_30, \
                         lsh1_42, msf0_39, msf1_39, msg_59, msg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * msf0_39[k]
                  - f_2 * msf1_39[k]
                  + f_3 * pc_z[k] * msg_59[k];

        t_84[k] = pa_y[k] * lsh0_42[k]
                  - f_8 * pc_y[k] * lsh1_42[k];

        t_85[k] = f_9 * lsg_30[k]
                  + f_3 * pc_y[k] * msg_60[k];

        t_86[k] = f_9 * lsg_15[k]
                  + f_3 * pc_z[k] * msg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, lsh0_24, lsh0_27, \
                         lsh0_47, lsg_32, lsh1_24, lsh1_27, lsh1_47, \
                         msg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * lsh0_24[k]
                  - f_8 * pc_z[k] * lsh1_24[k];

        t_88[k] = f_9 * lsg_32[k]
                  + f_3 * pc_y[k] * msg_62[k];

        t_89[k] = pa_y[k] * lsh0_47[k]
                  - f_8 * pc_y[k] * lsh1_47[k];

        t_90[k] = pa_z[k] * lsh0_27[k]
                  - f_8 * pc_z[k] * lsh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, lsh0_51, lsg_18, \
                         lsg_35, lsg_70, lsh1_51, msg_63, msg_65, \
                         msg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * lsg_18[k]
                  + f_3 * pc_z[k] * msg_63[k];

        t_92[k] = f_9 * lsg_35[k]
                  + f_3 * pc_y[k] * msg_65[k];

        t_93[k] = pa_y[k] * lsh0_51[k]
                  - f_8 * pc_y[k] * lsh1_51[k];

        t_94[k] = f_15 * lsg_70[k]
                  + f_3 * pc_x[k] * msg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, lsg_71, lsg_72, lsg_73, lsg_74, msg_71, \
                         msg_72, msg_73, msg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_15 * lsg_71[k]
                  + f_3 * pc_x[k] * msg_71[k];

        t_96[k] = f_15 * lsg_72[k]
                  + f_3 * pc_x[k] * msg_72[k];

        t_97[k] = f_15 * lsg_73[k]
                  + f_3 * pc_x[k] * msg_73[k];

        t_98[k] = f_15 * lsg_74[k]
                  + f_3 * pc_x[k] * msg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, lsh0_36, lsg_25, lsg_42, \
                         lsh1_36, msf0_48, msf1_48, msg_70, msg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * lsh0_36[k]
                  - f_8 * pc_z[k] * lsh1_36[k];

        t_100[k] = f_9 * lsg_25[k]
                   + f_3 * pc_z[k] * msg_70[k];

        t_101[k] = f_9 * lsg_42[k]
                   + f_6 * msf0_48[k]
                   - f_7 * msf1_48[k]
                   + f_3 * pc_y[k] * msg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, lsh0_62, lsg_43, lsg_44, lsh1_62, \
                         msf0_49, msf1_49, msg_73, msg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * lsg_43[k]
                   + f_4 * msf0_49[k]
                   - f_5 * msf1_49[k]
                   + f_3 * pc_y[k] * msg_73[k];

        t_103[k] = f_9 * lsg_44[k]
                   + f_3 * pc_y[k] * msg_74[k];

        t_104[k] = pa_y[k] * lsh0_62[k]
                   - f_8 * pc_y[k] * lsh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, lsg_30, lsg_75, \
                         msf0_50, msf1_50, msg_75, msg_76, msg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_15 * lsg_75[k]
                   + f_1 * msf0_50[k]
                   - f_2 * msf1_50[k]
                   + f_3 * pc_x[k] * msg_75[k];

        t_106[k] = f_3 * pc_y[k] * msg_75[k];

        t_107[k] = f_10 * lsg_30[k]
                   + f_3 * pc_z[k] * msg_75[k];

        t_108[k] = f_4 * msf0_50[k]
                   - f_5 * msf1_50[k]
                   + f_3 * pc_y[k] * msg_76[k];

        t_109[k] = f_3 * pc_y[k] * msg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, lsg_80, msf0_51, msf0_52, \
                         msf0_55, msf1_51, msf1_52, msf1_55, msg_78, msg_79, \
                         msg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_15 * lsg_80[k]
                   + f_6 * msf0_55[k]
                   - f_7 * msf1_55[k]
                   + f_3 * pc_x[k] * msg_80[k];

        t_111[k] = f_6 * msf0_51[k]
                   - f_7 * msf1_51[k]
                   + f_3 * pc_y[k] * msg_78[k];

        t_112[k] = f_4 * msf0_52[k]
                   - f_5 * msf1_52[k]
                   + f_3 * pc_y[k] * msg_79[k];

        t_113[k] = f_3 * pc_y[k] * msg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, lsg_84, lsg_85, lsg_86, lsg_87, \
                         msf0_59, msf1_59, msg_84, msg_85, msg_86, \
                         msg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * lsg_84[k]
                   + f_4 * msf0_59[k]
                   - f_5 * msf1_59[k]
                   + f_3 * pc_x[k] * msg_84[k];

        t_115[k] = f_15 * lsg_85[k]
                   + f_3 * pc_x[k] * msg_85[k];

        t_116[k] = f_15 * lsg_86[k]
                   + f_3 * pc_x[k] * msg_86[k];

        t_117[k] = f_15 * lsg_87[k]
                   + f_3 * pc_x[k] * msg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, lsg_89, msf0_56, msf0_57, \
                         msf1_56, msf1_57, msg_84, msg_85, msg_86, \
                         msg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * msg_84[k];

        t_119[k] = f_15 * lsg_89[k]
                   + f_3 * pc_x[k] * msg_89[k];

        t_120[k] = f_1 * msf0_56[k]
                   - f_2 * msf1_56[k]
                   + f_3 * pc_y[k] * msg_85[k];

        t_121[k] = f_13 * msf0_57[k]
                   - f_14 * msf1_57[k]
                   + f_3 * pc_y[k] * msg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, lsg_44, msf0_58, msf0_59, \
                         msf1_58, msf1_59, msg_87, msg_88, msg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * msf0_58[k]
                   - f_7 * msf1_58[k]
                   + f_3 * pc_y[k] * msg_87[k];

        t_123[k] = f_4 * msf0_59[k]
                   - f_5 * msf1_59[k]
                   + f_3 * pc_y[k] * msg_88[k];

        t_124[k] = f_3 * pc_y[k] * msg_89[k];

        t_125[k] = f_10 * lsg_44[k]
                   + f_1 * msf0_59[k]
                   - f_2 * msf1_59[k]
                   + f_3 * pc_z[k] * msg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, lsg_45, lsg_90, lsg_93, \
                         msf0_60, msf0_63, msf1_60, msf1_63, msg_90, \
                         msg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_16 * lsg_90[k]
                   + f_1 * msf0_60[k]
                   - f_2 * msf1_60[k]
                   + f_3 * pc_x[k] * msg_90[k];

        t_127[k] = f_11 * lsg_45[k]
                   + f_3 * pc_y[k] * msg_90[k];

        t_128[k] = f_3 * pc_z[k] * msg_90[k];

        t_129[k] = f_16 * lsg_93[k]
                   + f_6 * msf0_63[k]
                   - f_7 * msf1_63[k]
                   + f_3 * pc_x[k] * msg_93[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_63 = buffer.data(lsh0 + 63);
    const auto *lsh0_66 = buffer.data(lsh0 + 66);
    const auto *lsh0_69 = buffer.data(lsh0 + 69);
    const auto *lsh0_78 = buffer.data(lsh0 + 78);
    const auto *lsh0_105 = buffer.data(lsh0 + 105);
    const auto *lsh0_108 = buffer.data(lsh0 + 108);
    const auto *lsh0_110 = buffer.data(lsh0 + 110);
    const auto *lsh0_111 = buffer.data(lsh0 + 111);
    const auto *lsh0_114 = buffer.data(lsh0 + 114);
    const auto *lsh0_125 = buffer.data(lsh0 + 125);
    const auto *lsh0_126 = buffer.data(lsh0 + 126);
    const auto *lsh0_129 = buffer.data(lsh0 + 129);
    const auto *lsh0_132 = buffer.data(lsh0 + 132);
    const auto *lsh0_141 = buffer.data(lsh0 + 141);

    const auto *lsg_45 = buffer.data(lsg + 45);
    const auto *lsg_48 = buffer.data(lsg + 48);
    const auto *lsg_50 = buffer.data(lsg + 50);
    const auto *lsg_55 = buffer.data(lsg + 55);
    const auto *lsg_59 = buffer.data(lsg + 59);
    const auto *lsg_60 = buffer.data(lsg + 60);
    const auto *lsg_62 = buffer.data(lsg + 62);
    const auto *lsg_63 = buffer.data(lsg + 63);
    const auto *lsg_65 = buffer.data(lsg + 65);
    const auto *lsg_70 = buffer.data(lsg + 70);
    const auto *lsg_72 = buffer.data(lsg + 72);
    const auto *lsg_73 = buffer.data(lsg + 73);
    const auto *lsg_74 = buffer.data(lsg + 74);
    const auto *lsg_75 = buffer.data(lsg + 75);
    const auto *lsg_76 = buffer.data(lsg + 76);
    const auto *lsg_77 = buffer.data(lsg + 77);
    const auto *lsg_78 = buffer.data(lsg + 78);
    const auto *lsg_80 = buffer.data(lsg + 80);
    const auto *lsg_85 = buffer.data(lsg + 85);
    const auto *lsg_87 = buffer.data(lsg + 87);
    const auto *lsg_88 = buffer.data(lsg + 88);
    const auto *lsg_89 = buffer.data(lsg + 89);
    const auto *lsg_90 = buffer.data(lsg + 90);
    const auto *lsg_93 = buffer.data(lsg + 93);
    const auto *lsg_95 = buffer.data(lsg + 95);
    const auto *lsg_96 = buffer.data(lsg + 96);
    const auto *lsg_100 = buffer.data(lsg + 100);
    const auto *lsg_102 = buffer.data(lsg + 102);
    const auto *lsg_103 = buffer.data(lsg + 103);
    const auto *lsg_104 = buffer.data(lsg + 104);
    const auto *lsg_105 = buffer.data(lsg + 105);
    const auto *lsg_107 = buffer.data(lsg + 107);
    const auto *lsg_110 = buffer.data(lsg + 110);
    const auto *lsg_114 = buffer.data(lsg + 114);
    const auto *lsg_115 = buffer.data(lsg + 115);
    const auto *lsg_116 = buffer.data(lsg + 116);
    const auto *lsg_117 = buffer.data(lsg + 117);
    const auto *lsg_118 = buffer.data(lsg + 118);
    const auto *lsg_119 = buffer.data(lsg + 119);
    const auto *lsg_130 = buffer.data(lsg + 130);
    const auto *lsg_131 = buffer.data(lsg + 131);
    const auto *lsg_132 = buffer.data(lsg + 132);
    const auto *lsg_133 = buffer.data(lsg + 133);
    const auto *lsg_134 = buffer.data(lsg + 134);
    const auto *lsg_135 = buffer.data(lsg + 135);
    const auto *lsg_140 = buffer.data(lsg + 140);
    const auto *lsg_144 = buffer.data(lsg + 144);
    const auto *lsg_145 = buffer.data(lsg + 145);
    const auto *lsg_146 = buffer.data(lsg + 146);
    const auto *lsg_147 = buffer.data(lsg + 147);
    const auto *lsg_149 = buffer.data(lsg + 149);
    const auto *lsg_150 = buffer.data(lsg + 150);
    const auto *lsg_153 = buffer.data(lsg + 153);
    const auto *lsg_156 = buffer.data(lsg + 156);
    const auto *lsg_160 = buffer.data(lsg + 160);
    const auto *lsg_162 = buffer.data(lsg + 162);
    const auto *lsg_163 = buffer.data(lsg + 163);
    const auto *lsg_164 = buffer.data(lsg + 164);
    const auto *lsg_170 = buffer.data(lsg + 170);
    const auto *lsg_174 = buffer.data(lsg + 174);
    const auto *lsg_175 = buffer.data(lsg + 175);
    const auto *lsg_176 = buffer.data(lsg + 176);
    const auto *lsg_177 = buffer.data(lsg + 177);
    const auto *lsg_178 = buffer.data(lsg + 178);
    const auto *lsg_179 = buffer.data(lsg + 179);

    const auto *lsh1_63 = buffer.data(lsh1 + 63);
    const auto *lsh1_66 = buffer.data(lsh1 + 66);
    const auto *lsh1_69 = buffer.data(lsh1 + 69);
    const auto *lsh1_78 = buffer.data(lsh1 + 78);
    const auto *lsh1_105 = buffer.data(lsh1 + 105);
    const auto *lsh1_108 = buffer.data(lsh1 + 108);
    const auto *lsh1_110 = buffer.data(lsh1 + 110);
    const auto *lsh1_111 = buffer.data(lsh1 + 111);
    const auto *lsh1_114 = buffer.data(lsh1 + 114);
    const auto *lsh1_125 = buffer.data(lsh1 + 125);
    const auto *lsh1_126 = buffer.data(lsh1 + 126);
    const auto *lsh1_129 = buffer.data(lsh1 + 129);
    const auto *lsh1_132 = buffer.data(lsh1 + 132);
    const auto *lsh1_141 = buffer.data(lsh1 + 141);

    const auto *msf0_60 = buffer.data(msf0 + 60);
    const auto *msf0_62 = buffer.data(msf0 + 62);
    const auto *msf0_66 = buffer.data(msf0 + 66);
    const auto *msf0_67 = buffer.data(msf0 + 67);
    const auto *msf0_69 = buffer.data(msf0 + 69);
    const auto *msf0_75 = buffer.data(msf0 + 75);
    const auto *msf0_78 = buffer.data(msf0 + 78);
    const auto *msf0_79 = buffer.data(msf0 + 79);
    const auto *msf0_86 = buffer.data(msf0 + 86);
    const auto *msf0_88 = buffer.data(msf0 + 88);
    const auto *msf0_89 = buffer.data(msf0 + 89);
    const auto *msf0_90 = buffer.data(msf0 + 90);
    const auto *msf0_91 = buffer.data(msf0 + 91);
    const auto *msf0_92 = buffer.data(msf0 + 92);
    const auto *msf0_95 = buffer.data(msf0 + 95);
    const auto *msf0_96 = buffer.data(msf0 + 96);
    const auto *msf0_97 = buffer.data(msf0 + 97);
    const auto *msf0_98 = buffer.data(msf0 + 98);
    const auto *msf0_99 = buffer.data(msf0 + 99);
    const auto *msf0_100 = buffer.data(msf0 + 100);
    const auto *msf0_102 = buffer.data(msf0 + 102);
    const auto *msf0_103 = buffer.data(msf0 + 103);
    const auto *msf0_106 = buffer.data(msf0 + 106);
    const auto *msf0_107 = buffer.data(msf0 + 107);
    const auto *msf0_109 = buffer.data(msf0 + 109);
    const auto *msf0_115 = buffer.data(msf0 + 115);
    const auto *msf0_118 = buffer.data(msf0 + 118);
    const auto *msf0_119 = buffer.data(msf0 + 119);

    const auto *msf1_60 = buffer.data(msf1 + 60);
    const auto *msf1_62 = buffer.data(msf1 + 62);
    const auto *msf1_66 = buffer.data(msf1 + 66);
    const auto *msf1_67 = buffer.data(msf1 + 67);
    const auto *msf1_69 = buffer.data(msf1 + 69);
    const auto *msf1_75 = buffer.data(msf1 + 75);
    const auto *msf1_78 = buffer.data(msf1 + 78);
    const auto *msf1_79 = buffer.data(msf1 + 79);
    const auto *msf1_86 = buffer.data(msf1 + 86);
    const auto *msf1_88 = buffer.data(msf1 + 88);
    const auto *msf1_89 = buffer.data(msf1 + 89);
    const auto *msf1_90 = buffer.data(msf1 + 90);
    const auto *msf1_91 = buffer.data(msf1 + 91);
    const auto *msf1_92 = buffer.data(msf1 + 92);
    const auto *msf1_95 = buffer.data(msf1 + 95);
    const auto *msf1_96 = buffer.data(msf1 + 96);
    const auto *msf1_97 = buffer.data(msf1 + 97);
    const auto *msf1_98 = buffer.data(msf1 + 98);
    const auto *msf1_99 = buffer.data(msf1 + 99);
    const auto *msf1_100 = buffer.data(msf1 + 100);
    const auto *msf1_102 = buffer.data(msf1 + 102);
    const auto *msf1_103 = buffer.data(msf1 + 103);
    const auto *msf1_106 = buffer.data(msf1 + 106);
    const auto *msf1_107 = buffer.data(msf1 + 107);
    const auto *msf1_109 = buffer.data(msf1 + 109);
    const auto *msf1_115 = buffer.data(msf1 + 115);
    const auto *msf1_118 = buffer.data(msf1 + 118);
    const auto *msf1_119 = buffer.data(msf1 + 119);

    const auto *msg_91 = buffer.data(msg + 91);
    const auto *msg_92 = buffer.data(msg + 92);
    const auto *msg_93 = buffer.data(msg + 93);
    const auto *msg_95 = buffer.data(msg + 95);
    const auto *msg_96 = buffer.data(msg + 96);
    const auto *msg_100 = buffer.data(msg + 100);
    const auto *msg_101 = buffer.data(msg + 101);
    const auto *msg_102 = buffer.data(msg + 102);
    const auto *msg_103 = buffer.data(msg + 103);
    const auto *msg_104 = buffer.data(msg + 104);
    const auto *msg_105 = buffer.data(msg + 105);
    const auto *msg_107 = buffer.data(msg + 107);
    const auto *msg_108 = buffer.data(msg + 108);
    const auto *msg_110 = buffer.data(msg + 110);
    const auto *msg_114 = buffer.data(msg + 114);
    const auto *msg_115 = buffer.data(msg + 115);
    const auto *msg_116 = buffer.data(msg + 116);
    const auto *msg_117 = buffer.data(msg + 117);
    const auto *msg_118 = buffer.data(msg + 118);
    const auto *msg_119 = buffer.data(msg + 119);
    const auto *msg_120 = buffer.data(msg + 120);
    const auto *msg_122 = buffer.data(msg + 122);
    const auto *msg_123 = buffer.data(msg + 123);
    const auto *msg_125 = buffer.data(msg + 125);
    const auto *msg_130 = buffer.data(msg + 130);
    const auto *msg_131 = buffer.data(msg + 131);
    const auto *msg_132 = buffer.data(msg + 132);
    const auto *msg_133 = buffer.data(msg + 133);
    const auto *msg_134 = buffer.data(msg + 134);
    const auto *msg_135 = buffer.data(msg + 135);
    const auto *msg_136 = buffer.data(msg + 136);
    const auto *msg_137 = buffer.data(msg + 137);
    const auto *msg_138 = buffer.data(msg + 138);
    const auto *msg_139 = buffer.data(msg + 139);
    const auto *msg_140 = buffer.data(msg + 140);
    const auto *msg_144 = buffer.data(msg + 144);
    const auto *msg_145 = buffer.data(msg + 145);
    const auto *msg_146 = buffer.data(msg + 146);
    const auto *msg_147 = buffer.data(msg + 147);
    const auto *msg_148 = buffer.data(msg + 148);
    const auto *msg_149 = buffer.data(msg + 149);
    const auto *msg_150 = buffer.data(msg + 150);
    const auto *msg_151 = buffer.data(msg + 151);
    const auto *msg_152 = buffer.data(msg + 152);
    const auto *msg_153 = buffer.data(msg + 153);
    const auto *msg_155 = buffer.data(msg + 155);
    const auto *msg_156 = buffer.data(msg + 156);
    const auto *msg_160 = buffer.data(msg + 160);
    const auto *msg_161 = buffer.data(msg + 161);
    const auto *msg_162 = buffer.data(msg + 162);
    const auto *msg_163 = buffer.data(msg + 163);
    const auto *msg_164 = buffer.data(msg + 164);
    const auto *msg_165 = buffer.data(msg + 165);
    const auto *msg_167 = buffer.data(msg + 167);
    const auto *msg_168 = buffer.data(msg + 168);
    const auto *msg_170 = buffer.data(msg + 170);
    const auto *msg_174 = buffer.data(msg + 174);
    const auto *msg_175 = buffer.data(msg + 175);
    const auto *msg_176 = buffer.data(msg + 176);
    const auto *msg_177 = buffer.data(msg + 177);
    const auto *msg_178 = buffer.data(msg + 178);
    const auto *msg_179 = buffer.data(msg + 179);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_x, pc_z, lsg_96, msf0_60, msf0_66, \
                         msf1_60, msf1_66, msg_91, msg_92, msg_93, \
                         msg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * msg_91[k];

        t_131[k] = f_4 * msf0_60[k]
                   - f_5 * msf1_60[k]
                   + f_3 * pc_z[k] * msg_92[k];

        t_132[k] = f_16 * lsg_96[k]
                   + f_4 * msf0_66[k]
                   - f_5 * msf1_66[k]
                   + f_3 * pc_x[k] * msg_96[k];

        t_133[k] = f_3 * pc_z[k] * msg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, lsg_50, lsg_100, \
                         msf0_62, msf1_62, msg_95, msg_96, msg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * lsg_50[k]
                   + f_3 * pc_y[k] * msg_95[k];

        t_135[k] = f_6 * msf0_62[k]
                   - f_7 * msf1_62[k]
                   + f_3 * pc_z[k] * msg_95[k];

        t_136[k] = f_16 * lsg_100[k]
                   + f_3 * pc_x[k] * msg_100[k];

        t_137[k] = f_3 * pc_z[k] * msg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_y, lsg_55, lsg_102, lsg_103, \
                         lsg_104, msf0_66, msf1_66, msg_100, msg_102, msg_103, \
                         msg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_16 * lsg_102[k]
                   + f_3 * pc_x[k] * msg_102[k];

        t_139[k] = f_16 * lsg_103[k]
                   + f_3 * pc_x[k] * msg_103[k];

        t_140[k] = f_16 * lsg_104[k]
                   + f_3 * pc_x[k] * msg_104[k];

        t_141[k] = f_11 * lsg_55[k]
                   + f_1 * msf0_66[k]
                   - f_2 * msf1_66[k]
                   + f_3 * pc_y[k] * msg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, pc_z, lsg_59, msf0_66, msf0_67, \
                         msf1_66, msf1_67, msg_100, msg_101, msg_102, \
                         msg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * msg_100[k];

        t_143[k] = f_4 * msf0_66[k]
                   - f_5 * msf1_66[k]
                   + f_3 * pc_z[k] * msg_101[k];

        t_144[k] = f_6 * msf0_67[k]
                   - f_7 * msf1_67[k]
                   + f_3 * pc_z[k] * msg_102[k];

        t_145[k] = f_11 * lsg_59[k]
                   + f_3 * pc_y[k] * msg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, lsh0_63, lsg_45, \
                         lsg_60, lsh1_63, msf0_69, msf1_69, msg_104, \
                         msg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * msf0_69[k]
                   - f_2 * msf1_69[k]
                   + f_3 * pc_z[k] * msg_104[k];

        t_147[k] = pa_z[k] * lsh0_63[k]
                   - f_8 * pc_z[k] * lsh1_63[k];

        t_148[k] = f_10 * lsg_60[k]
                   + f_3 * pc_y[k] * msg_105[k];

        t_149[k] = f_9 * lsg_45[k]
                   + f_3 * pc_z[k] * msg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_z, pc_x, pc_y, pc_z, lsh0_66, lsg_62, \
                         lsg_110, lsh1_66, msf0_75, msf1_75, msg_107, \
                         msg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * lsh0_66[k]
                   - f_8 * pc_z[k] * lsh1_66[k];

        t_151[k] = f_10 * lsg_62[k]
                   + f_3 * pc_y[k] * msg_107[k];

        t_152[k] = f_16 * lsg_110[k]
                   + f_6 * msf0_75[k]
                   - f_7 * msf1_75[k]
                   + f_3 * pc_x[k] * msg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, lsh0_69, lsg_48, lsg_65, \
                         lsh1_69, msg_108, msg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * lsh0_69[k]
                   - f_8 * pc_z[k] * lsh1_69[k];

        t_154[k] = f_9 * lsg_48[k]
                   + f_3 * pc_z[k] * msg_108[k];

        t_155[k] = f_10 * lsg_65[k]
                   + f_3 * pc_y[k] * msg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, lsg_114, lsg_115, lsg_116, lsg_117, \
                         msf0_79, msf1_79, msg_114, msg_115, msg_116, \
                         msg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * lsg_114[k]
                   + f_4 * msf0_79[k]
                   - f_5 * msf1_79[k]
                   + f_3 * pc_x[k] * msg_114[k];

        t_157[k] = f_16 * lsg_115[k]
                   + f_3 * pc_x[k] * msg_115[k];

        t_158[k] = f_16 * lsg_116[k]
                   + f_3 * pc_x[k] * msg_116[k];

        t_159[k] = f_16 * lsg_117[k]
                   + f_3 * pc_x[k] * msg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, lsh0_78, lsg_55, \
                         lsg_118, lsg_119, lsh1_78, msg_115, msg_118, \
                         msg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * lsg_118[k]
                   + f_3 * pc_x[k] * msg_118[k];

        t_161[k] = f_16 * lsg_119[k]
                   + f_3 * pc_x[k] * msg_119[k];

        t_162[k] = pa_z[k] * lsh0_78[k]
                   - f_8 * pc_z[k] * lsh1_78[k];

        t_163[k] = f_9 * lsg_55[k]
                   + f_3 * pc_z[k] * msg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pc_y, lsg_72, lsg_73, lsg_74, msf0_78, msf0_79, \
                         msf1_78, msf1_79, msg_117, msg_118, msg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * lsg_72[k]
                   + f_6 * msf0_78[k]
                   - f_7 * msf1_78[k]
                   + f_3 * pc_y[k] * msg_117[k];

        t_165[k] = f_10 * lsg_73[k]
                   + f_4 * msf0_79[k]
                   - f_5 * msf1_79[k]
                   + f_3 * pc_y[k] * msg_118[k];

        t_166[k] = f_10 * lsg_74[k]
                   + f_3 * pc_y[k] * msg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pc_y, pc_z, lsh0_105, lsg_59, \
                         lsg_60, lsg_75, lsh1_105, msf0_79, msf1_79, msg_119, \
                         msg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * lsg_59[k]
                   + f_1 * msf0_79[k]
                   - f_2 * msf1_79[k]
                   + f_3 * pc_z[k] * msg_119[k];

        t_168[k] = pa_y[k] * lsh0_105[k]
                   - f_8 * pc_y[k] * lsh1_105[k];

        t_169[k] = f_9 * lsg_75[k]
                   + f_3 * pc_y[k] * msg_120[k];

        t_170[k] = f_10 * lsg_60[k]
                   + f_3 * pc_z[k] * msg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, pc_y, lsh0_108, lsh0_110, lsh0_111, \
                         lsg_76, lsg_77, lsg_78, lsh1_108, lsh1_110, lsh1_111, \
                         msg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * lsh0_108[k]
                   + f_10 * lsg_76[k]
                   - f_8 * pc_y[k] * lsh1_108[k];

        t_172[k] = f_9 * lsg_77[k]
                   + f_3 * pc_y[k] * msg_122[k];

        t_173[k] = pa_y[k] * lsh0_110[k]
                   - f_8 * pc_y[k] * lsh1_110[k];

        t_174[k] = pa_y[k] * lsh0_111[k]
                   + f_11 * lsg_78[k]
                   - f_8 * pc_y[k] * lsh1_111[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, pc_z, lsh0_114, lsg_63, \
                         lsg_80, lsg_130, lsh1_114, msg_123, msg_125, \
                         msg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * lsg_63[k]
                   + f_3 * pc_z[k] * msg_123[k];

        t_176[k] = f_9 * lsg_80[k]
                   + f_3 * pc_y[k] * msg_125[k];

        t_177[k] = pa_y[k] * lsh0_114[k]
                   - f_8 * pc_y[k] * lsh1_114[k];

        t_178[k] = f_16 * lsg_130[k]
                   + f_3 * pc_x[k] * msg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, lsg_131, lsg_132, lsg_133, lsg_134, \
                         msg_131, msg_132, msg_133, msg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_16 * lsg_131[k]
                   + f_3 * pc_x[k] * msg_131[k];

        t_180[k] = f_16 * lsg_132[k]
                   + f_3 * pc_x[k] * msg_132[k];

        t_181[k] = f_16 * lsg_133[k]
                   + f_3 * pc_x[k] * msg_133[k];

        t_182[k] = f_16 * lsg_134[k]
                   + f_3 * pc_x[k] * msg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, lsg_70, lsg_85, lsg_87, msf0_86, \
                         msf0_88, msf1_86, msf1_88, msg_130, msg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * lsg_85[k]
                   + f_1 * msf0_86[k]
                   - f_2 * msf1_86[k]
                   + f_3 * pc_y[k] * msg_130[k];

        t_184[k] = f_10 * lsg_70[k]
                   + f_3 * pc_z[k] * msg_130[k];

        t_185[k] = f_9 * lsg_87[k]
                   + f_6 * msf0_88[k]
                   - f_7 * msf1_88[k]
                   + f_3 * pc_y[k] * msg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pc_y, lsh0_125, lsg_88, lsg_89, lsh1_125, \
                         msf0_89, msf1_89, msg_133, msg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * lsg_88[k]
                   + f_4 * msf0_89[k]
                   - f_5 * msf1_89[k]
                   + f_3 * pc_y[k] * msg_133[k];

        t_187[k] = f_9 * lsg_89[k]
                   + f_3 * pc_y[k] * msg_134[k];

        t_188[k] = pa_y[k] * lsh0_125[k]
                   - f_8 * pc_y[k] * lsh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, lsg_75, lsg_135, \
                         msf0_90, msf1_90, msg_135, msg_136, msg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_16 * lsg_135[k]
                   + f_1 * msf0_90[k]
                   - f_2 * msf1_90[k]
                   + f_3 * pc_x[k] * msg_135[k];

        t_190[k] = f_3 * pc_y[k] * msg_135[k];

        t_191[k] = f_11 * lsg_75[k]
                   + f_3 * pc_z[k] * msg_135[k];

        t_192[k] = f_4 * msf0_90[k]
                   - f_5 * msf1_90[k]
                   + f_3 * pc_y[k] * msg_136[k];

        t_193[k] = f_3 * pc_y[k] * msg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, lsg_140, msf0_91, msf0_92, \
                         msf0_95, msf1_91, msf1_92, msf1_95, msg_138, msg_139, \
                         msg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_16 * lsg_140[k]
                   + f_6 * msf0_95[k]
                   - f_7 * msf1_95[k]
                   + f_3 * pc_x[k] * msg_140[k];

        t_195[k] = f_6 * msf0_91[k]
                   - f_7 * msf1_91[k]
                   + f_3 * pc_y[k] * msg_138[k];

        t_196[k] = f_4 * msf0_92[k]
                   - f_5 * msf1_92[k]
                   + f_3 * pc_y[k] * msg_139[k];

        t_197[k] = f_3 * pc_y[k] * msg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, lsg_144, lsg_145, lsg_146, lsg_147, \
                         msf0_99, msf1_99, msg_144, msg_145, msg_146, \
                         msg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_16 * lsg_144[k]
                   + f_4 * msf0_99[k]
                   - f_5 * msf1_99[k]
                   + f_3 * pc_x[k] * msg_144[k];

        t_199[k] = f_16 * lsg_145[k]
                   + f_3 * pc_x[k] * msg_145[k];

        t_200[k] = f_16 * lsg_146[k]
                   + f_3 * pc_x[k] * msg_146[k];

        t_201[k] = f_16 * lsg_147[k]
                   + f_3 * pc_x[k] * msg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_y, lsg_149, msf0_96, msf0_97, \
                         msf1_96, msf1_97, msg_144, msg_145, msg_146, \
                         msg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * pc_y[k] * msg_144[k];

        t_203[k] = f_16 * lsg_149[k]
                   + f_3 * pc_x[k] * msg_149[k];

        t_204[k] = f_1 * msf0_96[k]
                   - f_2 * msf1_96[k]
                   + f_3 * pc_y[k] * msg_145[k];

        t_205[k] = f_13 * msf0_97[k]
                   - f_14 * msf1_97[k]
                   + f_3 * pc_y[k] * msg_146[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, lsg_89, msf0_98, msf0_99, \
                         msf1_98, msf1_99, msg_147, msg_148, msg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * msf0_98[k]
                   - f_7 * msf1_98[k]
                   + f_3 * pc_y[k] * msg_147[k];

        t_207[k] = f_4 * msf0_99[k]
                   - f_5 * msf1_99[k]
                   + f_3 * pc_y[k] * msg_148[k];

        t_208[k] = f_3 * pc_y[k] * msg_149[k];

        t_209[k] = f_11 * lsg_89[k]
                   + f_1 * msf0_99[k]
                   - f_2 * msf1_99[k]
                   + f_3 * pc_z[k] * msg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, pc_y, pc_z, lsg_90, lsg_150, \
                         lsg_153, msf0_100, msf0_103, msf1_100, msf1_103, msg_150, \
                         msg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * lsg_150[k]
                   + f_1 * msf0_100[k]
                   - f_2 * msf1_100[k]
                   + f_3 * pc_x[k] * msg_150[k];

        t_211[k] = f_18 * lsg_90[k]
                   + f_3 * pc_y[k] * msg_150[k];

        t_212[k] = f_3 * pc_z[k] * msg_150[k];

        t_213[k] = f_17 * lsg_153[k]
                   + f_6 * msf0_103[k]
                   - f_7 * msf1_103[k]
                   + f_3 * pc_x[k] * msg_153[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pc_x, pc_z, lsg_156, msf0_100, msf0_106, \
                         msf1_100, msf1_106, msg_151, msg_152, msg_153, \
                         msg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_z[k] * msg_151[k];

        t_215[k] = f_4 * msf0_100[k]
                   - f_5 * msf1_100[k]
                   + f_3 * pc_z[k] * msg_152[k];

        t_216[k] = f_17 * lsg_156[k]
                   + f_4 * msf0_106[k]
                   - f_5 * msf1_106[k]
                   + f_3 * pc_x[k] * msg_156[k];

        t_217[k] = f_3 * pc_z[k] * msg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, lsg_95, lsg_160, \
                         msf0_102, msf1_102, msg_155, msg_156, \
                         msg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_18 * lsg_95[k]
                   + f_3 * pc_y[k] * msg_155[k];

        t_219[k] = f_6 * msf0_102[k]
                   - f_7 * msf1_102[k]
                   + f_3 * pc_z[k] * msg_155[k];

        t_220[k] = f_17 * lsg_160[k]
                   + f_3 * pc_x[k] * msg_160[k];

        t_221[k] = f_3 * pc_z[k] * msg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, lsg_100, lsg_162, lsg_163, \
                         lsg_164, msf0_106, msf1_106, msg_160, msg_162, msg_163, \
                         msg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_17 * lsg_162[k]
                   + f_3 * pc_x[k] * msg_162[k];

        t_223[k] = f_17 * lsg_163[k]
                   + f_3 * pc_x[k] * msg_163[k];

        t_224[k] = f_17 * lsg_164[k]
                   + f_3 * pc_x[k] * msg_164[k];

        t_225[k] = f_18 * lsg_100[k]
                   + f_1 * msf0_106[k]
                   - f_2 * msf1_106[k]
                   + f_3 * pc_y[k] * msg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, lsg_104, msf0_106, msf0_107, \
                         msf1_106, msf1_107, msg_160, msg_161, msg_162, \
                         msg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * msg_160[k];

        t_227[k] = f_4 * msf0_106[k]
                   - f_5 * msf1_106[k]
                   + f_3 * pc_z[k] * msg_161[k];

        t_228[k] = f_6 * msf0_107[k]
                   - f_7 * msf1_107[k]
                   + f_3 * pc_z[k] * msg_162[k];

        t_229[k] = f_18 * lsg_104[k]
                   + f_3 * pc_y[k] * msg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pa_z, pc_y, pc_z, lsh0_126, lsg_90, \
                         lsg_105, lsh1_126, msf0_109, msf1_109, msg_164, \
                         msg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * msf0_109[k]
                   - f_2 * msf1_109[k]
                   + f_3 * pc_z[k] * msg_164[k];

        t_231[k] = pa_z[k] * lsh0_126[k]
                   - f_8 * pc_z[k] * lsh1_126[k];

        t_232[k] = f_11 * lsg_105[k]
                   + f_3 * pc_y[k] * msg_165[k];

        t_233[k] = f_9 * lsg_90[k]
                   + f_3 * pc_z[k] * msg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_z, pc_x, pc_y, pc_z, lsh0_129, lsg_107, \
                         lsg_170, lsh1_129, msf0_115, msf1_115, msg_167, \
                         msg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * lsh0_129[k]
                   - f_8 * pc_z[k] * lsh1_129[k];

        t_235[k] = f_11 * lsg_107[k]
                   + f_3 * pc_y[k] * msg_167[k];

        t_236[k] = f_17 * lsg_170[k]
                   + f_6 * msf0_115[k]
                   - f_7 * msf1_115[k]
                   + f_3 * pc_x[k] * msg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_z, pc_y, pc_z, lsh0_132, lsg_93, lsg_110, \
                         lsh1_132, msg_168, msg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * lsh0_132[k]
                   - f_8 * pc_z[k] * lsh1_132[k];

        t_238[k] = f_9 * lsg_93[k]
                   + f_3 * pc_z[k] * msg_168[k];

        t_239[k] = f_11 * lsg_110[k]
                   + f_3 * pc_y[k] * msg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, lsg_174, lsg_175, lsg_176, lsg_177, \
                         msf0_119, msf1_119, msg_174, msg_175, msg_176, \
                         msg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * lsg_174[k]
                   + f_4 * msf0_119[k]
                   - f_5 * msf1_119[k]
                   + f_3 * pc_x[k] * msg_174[k];

        t_241[k] = f_17 * lsg_175[k]
                   + f_3 * pc_x[k] * msg_175[k];

        t_242[k] = f_17 * lsg_176[k]
                   + f_3 * pc_x[k] * msg_176[k];

        t_243[k] = f_17 * lsg_177[k]
                   + f_3 * pc_x[k] * msg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pc_x, pc_z, lsh0_141, lsg_100, \
                         lsg_178, lsg_179, lsh1_141, msg_175, msg_178, \
                         msg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_17 * lsg_178[k]
                   + f_3 * pc_x[k] * msg_178[k];

        t_245[k] = f_17 * lsg_179[k]
                   + f_3 * pc_x[k] * msg_179[k];

        t_246[k] = pa_z[k] * lsh0_141[k]
                   - f_8 * pc_z[k] * lsh1_141[k];

        t_247[k] = f_9 * lsg_100[k]
                   + f_3 * pc_z[k] * msg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, lsg_117, lsg_118, lsg_119, msf0_118, \
                         msf0_119, msf1_118, msf1_119, msg_177, msg_178, \
                         msg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * lsg_117[k]
                   + f_6 * msf0_118[k]
                   - f_7 * msf1_118[k]
                   + f_3 * pc_y[k] * msg_177[k];

        t_249[k] = f_11 * lsg_118[k]
                   + f_4 * msf0_119[k]
                   - f_5 * msf1_119[k]
                   + f_3 * pc_y[k] * msg_178[k];

        t_250[k] = f_11 * lsg_119[k]
                   + f_3 * pc_y[k] * msg_179[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_189 = buffer.data(lsh0 + 189);
    const auto *lsh0_192 = buffer.data(lsh0 + 192);
    const auto *lsh0_194 = buffer.data(lsh0 + 194);
    const auto *lsh0_195 = buffer.data(lsh0 + 195);
    const auto *lsh0_198 = buffer.data(lsh0 + 198);
    const auto *lsh0_209 = buffer.data(lsh0 + 209);
    const auto *lsh0_210 = buffer.data(lsh0 + 210);
    const auto *lsh0_213 = buffer.data(lsh0 + 213);
    const auto *lsh0_216 = buffer.data(lsh0 + 216);
    const auto *lsh0_225 = buffer.data(lsh0 + 225);

    const auto *lsg_104 = buffer.data(lsg + 104);
    const auto *lsg_105 = buffer.data(lsg + 105);
    const auto *lsg_108 = buffer.data(lsg + 108);
    const auto *lsg_115 = buffer.data(lsg + 115);
    const auto *lsg_119 = buffer.data(lsg + 119);
    const auto *lsg_120 = buffer.data(lsg + 120);
    const auto *lsg_122 = buffer.data(lsg + 122);
    const auto *lsg_123 = buffer.data(lsg + 123);
    const auto *lsg_125 = buffer.data(lsg + 125);
    const auto *lsg_130 = buffer.data(lsg + 130);
    const auto *lsg_132 = buffer.data(lsg + 132);
    const auto *lsg_133 = buffer.data(lsg + 133);
    const auto *lsg_134 = buffer.data(lsg + 134);
    const auto *lsg_135 = buffer.data(lsg + 135);
    const auto *lsg_136 = buffer.data(lsg + 136);
    const auto *lsg_137 = buffer.data(lsg + 137);
    const auto *lsg_138 = buffer.data(lsg + 138);
    const auto *lsg_140 = buffer.data(lsg + 140);
    const auto *lsg_145 = buffer.data(lsg + 145);
    const auto *lsg_147 = buffer.data(lsg + 147);
    const auto *lsg_148 = buffer.data(lsg + 148);
    const auto *lsg_149 = buffer.data(lsg + 149);
    const auto *lsg_150 = buffer.data(lsg + 150);
    const auto *lsg_153 = buffer.data(lsg + 153);
    const auto *lsg_155 = buffer.data(lsg + 155);
    const auto *lsg_160 = buffer.data(lsg + 160);
    const auto *lsg_164 = buffer.data(lsg + 164);
    const auto *lsg_165 = buffer.data(lsg + 165);
    const auto *lsg_167 = buffer.data(lsg + 167);
    const auto *lsg_168 = buffer.data(lsg + 168);
    const auto *lsg_170 = buffer.data(lsg + 170);
    const auto *lsg_177 = buffer.data(lsg + 177);
    const auto *lsg_178 = buffer.data(lsg + 178);
    const auto *lsg_179 = buffer.data(lsg + 179);
    const auto *lsg_180 = buffer.data(lsg + 180);
    const auto *lsg_182 = buffer.data(lsg + 182);
    const auto *lsg_183 = buffer.data(lsg + 183);
    const auto *lsg_185 = buffer.data(lsg + 185);
    const auto *lsg_186 = buffer.data(lsg + 186);
    const auto *lsg_189 = buffer.data(lsg + 189);
    const auto *lsg_190 = buffer.data(lsg + 190);
    const auto *lsg_191 = buffer.data(lsg + 191);
    const auto *lsg_192 = buffer.data(lsg + 192);
    const auto *lsg_193 = buffer.data(lsg + 193);
    const auto *lsg_194 = buffer.data(lsg + 194);
    const auto *lsg_205 = buffer.data(lsg + 205);
    const auto *lsg_206 = buffer.data(lsg + 206);
    const auto *lsg_207 = buffer.data(lsg + 207);
    const auto *lsg_208 = buffer.data(lsg + 208);
    const auto *lsg_209 = buffer.data(lsg + 209);
    const auto *lsg_210 = buffer.data(lsg + 210);
    const auto *lsg_215 = buffer.data(lsg + 215);
    const auto *lsg_219 = buffer.data(lsg + 219);
    const auto *lsg_220 = buffer.data(lsg + 220);
    const auto *lsg_221 = buffer.data(lsg + 221);
    const auto *lsg_222 = buffer.data(lsg + 222);
    const auto *lsg_224 = buffer.data(lsg + 224);
    const auto *lsg_225 = buffer.data(lsg + 225);
    const auto *lsg_228 = buffer.data(lsg + 228);
    const auto *lsg_231 = buffer.data(lsg + 231);
    const auto *lsg_235 = buffer.data(lsg + 235);
    const auto *lsg_237 = buffer.data(lsg + 237);
    const auto *lsg_238 = buffer.data(lsg + 238);
    const auto *lsg_239 = buffer.data(lsg + 239);
    const auto *lsg_245 = buffer.data(lsg + 245);
    const auto *lsg_249 = buffer.data(lsg + 249);
    const auto *lsg_250 = buffer.data(lsg + 250);
    const auto *lsg_251 = buffer.data(lsg + 251);
    const auto *lsg_252 = buffer.data(lsg + 252);
    const auto *lsg_253 = buffer.data(lsg + 253);
    const auto *lsg_254 = buffer.data(lsg + 254);
    const auto *lsg_255 = buffer.data(lsg + 255);
    const auto *lsg_258 = buffer.data(lsg + 258);
    const auto *lsg_260 = buffer.data(lsg + 260);
    const auto *lsg_261 = buffer.data(lsg + 261);
    const auto *lsg_264 = buffer.data(lsg + 264);
    const auto *lsg_265 = buffer.data(lsg + 265);
    const auto *lsg_266 = buffer.data(lsg + 266);

    const auto *lsh1_189 = buffer.data(lsh1 + 189);
    const auto *lsh1_192 = buffer.data(lsh1 + 192);
    const auto *lsh1_194 = buffer.data(lsh1 + 194);
    const auto *lsh1_195 = buffer.data(lsh1 + 195);
    const auto *lsh1_198 = buffer.data(lsh1 + 198);
    const auto *lsh1_209 = buffer.data(lsh1 + 209);
    const auto *lsh1_210 = buffer.data(lsh1 + 210);
    const auto *lsh1_213 = buffer.data(lsh1 + 213);
    const auto *lsh1_216 = buffer.data(lsh1 + 216);
    const auto *lsh1_225 = buffer.data(lsh1 + 225);

    const auto *msf0_119 = buffer.data(msf0 + 119);
    const auto *msf0_120 = buffer.data(msf0 + 120);
    const auto *msf0_123 = buffer.data(msf0 + 123);
    const auto *msf0_125 = buffer.data(msf0 + 125);
    const auto *msf0_126 = buffer.data(msf0 + 126);
    const auto *msf0_128 = buffer.data(msf0 + 128);
    const auto *msf0_129 = buffer.data(msf0 + 129);
    const auto *msf0_136 = buffer.data(msf0 + 136);
    const auto *msf0_138 = buffer.data(msf0 + 138);
    const auto *msf0_139 = buffer.data(msf0 + 139);
    const auto *msf0_140 = buffer.data(msf0 + 140);
    const auto *msf0_141 = buffer.data(msf0 + 141);
    const auto *msf0_142 = buffer.data(msf0 + 142);
    const auto *msf0_145 = buffer.data(msf0 + 145);
    const auto *msf0_146 = buffer.data(msf0 + 146);
    const auto *msf0_147 = buffer.data(msf0 + 147);
    const auto *msf0_148 = buffer.data(msf0 + 148);
    const auto *msf0_149 = buffer.data(msf0 + 149);
    const auto *msf0_150 = buffer.data(msf0 + 150);
    const auto *msf0_152 = buffer.data(msf0 + 152);
    const auto *msf0_153 = buffer.data(msf0 + 153);
    const auto *msf0_156 = buffer.data(msf0 + 156);
    const auto *msf0_157 = buffer.data(msf0 + 157);
    const auto *msf0_159 = buffer.data(msf0 + 159);
    const auto *msf0_165 = buffer.data(msf0 + 165);
    const auto *msf0_168 = buffer.data(msf0 + 168);
    const auto *msf0_169 = buffer.data(msf0 + 169);
    const auto *msf0_170 = buffer.data(msf0 + 170);
    const auto *msf0_173 = buffer.data(msf0 + 173);
    const auto *msf0_175 = buffer.data(msf0 + 175);
    const auto *msf0_176 = buffer.data(msf0 + 176);
    const auto *msf0_179 = buffer.data(msf0 + 179);

    const auto *msf1_119 = buffer.data(msf1 + 119);
    const auto *msf1_120 = buffer.data(msf1 + 120);
    const auto *msf1_123 = buffer.data(msf1 + 123);
    const auto *msf1_125 = buffer.data(msf1 + 125);
    const auto *msf1_126 = buffer.data(msf1 + 126);
    const auto *msf1_128 = buffer.data(msf1 + 128);
    const auto *msf1_129 = buffer.data(msf1 + 129);
    const auto *msf1_136 = buffer.data(msf1 + 136);
    const auto *msf1_138 = buffer.data(msf1 + 138);
    const auto *msf1_139 = buffer.data(msf1 + 139);
    const auto *msf1_140 = buffer.data(msf1 + 140);
    const auto *msf1_141 = buffer.data(msf1 + 141);
    const auto *msf1_142 = buffer.data(msf1 + 142);
    const auto *msf1_145 = buffer.data(msf1 + 145);
    const auto *msf1_146 = buffer.data(msf1 + 146);
    const auto *msf1_147 = buffer.data(msf1 + 147);
    const auto *msf1_148 = buffer.data(msf1 + 148);
    const auto *msf1_149 = buffer.data(msf1 + 149);
    const auto *msf1_150 = buffer.data(msf1 + 150);
    const auto *msf1_152 = buffer.data(msf1 + 152);
    const auto *msf1_153 = buffer.data(msf1 + 153);
    const auto *msf1_156 = buffer.data(msf1 + 156);
    const auto *msf1_157 = buffer.data(msf1 + 157);
    const auto *msf1_159 = buffer.data(msf1 + 159);
    const auto *msf1_165 = buffer.data(msf1 + 165);
    const auto *msf1_168 = buffer.data(msf1 + 168);
    const auto *msf1_169 = buffer.data(msf1 + 169);
    const auto *msf1_170 = buffer.data(msf1 + 170);
    const auto *msf1_173 = buffer.data(msf1 + 173);
    const auto *msf1_175 = buffer.data(msf1 + 175);
    const auto *msf1_176 = buffer.data(msf1 + 176);
    const auto *msf1_179 = buffer.data(msf1 + 179);

    const auto *msg_179 = buffer.data(msg + 179);
    const auto *msg_180 = buffer.data(msg + 180);
    const auto *msg_182 = buffer.data(msg + 182);
    const auto *msg_183 = buffer.data(msg + 183);
    const auto *msg_185 = buffer.data(msg + 185);
    const auto *msg_186 = buffer.data(msg + 186);
    const auto *msg_189 = buffer.data(msg + 189);
    const auto *msg_190 = buffer.data(msg + 190);
    const auto *msg_191 = buffer.data(msg + 191);
    const auto *msg_192 = buffer.data(msg + 192);
    const auto *msg_193 = buffer.data(msg + 193);
    const auto *msg_194 = buffer.data(msg + 194);
    const auto *msg_195 = buffer.data(msg + 195);
    const auto *msg_197 = buffer.data(msg + 197);
    const auto *msg_198 = buffer.data(msg + 198);
    const auto *msg_200 = buffer.data(msg + 200);
    const auto *msg_205 = buffer.data(msg + 205);
    const auto *msg_206 = buffer.data(msg + 206);
    const auto *msg_207 = buffer.data(msg + 207);
    const auto *msg_208 = buffer.data(msg + 208);
    const auto *msg_209 = buffer.data(msg + 209);
    const auto *msg_210 = buffer.data(msg + 210);
    const auto *msg_211 = buffer.data(msg + 211);
    const auto *msg_212 = buffer.data(msg + 212);
    const auto *msg_213 = buffer.data(msg + 213);
    const auto *msg_214 = buffer.data(msg + 214);
    const auto *msg_215 = buffer.data(msg + 215);
    const auto *msg_219 = buffer.data(msg + 219);
    const auto *msg_220 = buffer.data(msg + 220);
    const auto *msg_221 = buffer.data(msg + 221);
    const auto *msg_222 = buffer.data(msg + 222);
    const auto *msg_223 = buffer.data(msg + 223);
    const auto *msg_224 = buffer.data(msg + 224);
    const auto *msg_225 = buffer.data(msg + 225);
    const auto *msg_226 = buffer.data(msg + 226);
    const auto *msg_227 = buffer.data(msg + 227);
    const auto *msg_228 = buffer.data(msg + 228);
    const auto *msg_230 = buffer.data(msg + 230);
    const auto *msg_231 = buffer.data(msg + 231);
    const auto *msg_235 = buffer.data(msg + 235);
    const auto *msg_236 = buffer.data(msg + 236);
    const auto *msg_237 = buffer.data(msg + 237);
    const auto *msg_238 = buffer.data(msg + 238);
    const auto *msg_239 = buffer.data(msg + 239);
    const auto *msg_240 = buffer.data(msg + 240);
    const auto *msg_242 = buffer.data(msg + 242);
    const auto *msg_243 = buffer.data(msg + 243);
    const auto *msg_245 = buffer.data(msg + 245);
    const auto *msg_249 = buffer.data(msg + 249);
    const auto *msg_250 = buffer.data(msg + 250);
    const auto *msg_251 = buffer.data(msg + 251);
    const auto *msg_252 = buffer.data(msg + 252);
    const auto *msg_253 = buffer.data(msg + 253);
    const auto *msg_254 = buffer.data(msg + 254);
    const auto *msg_255 = buffer.data(msg + 255);
    const auto *msg_257 = buffer.data(msg + 257);
    const auto *msg_258 = buffer.data(msg + 258);
    const auto *msg_260 = buffer.data(msg + 260);
    const auto *msg_261 = buffer.data(msg + 261);
    const auto *msg_264 = buffer.data(msg + 264);
    const auto *msg_265 = buffer.data(msg + 265);
    const auto *msg_266 = buffer.data(msg + 266);

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_y, pc_z, lsg_104, lsg_120, lsg_180, \
                         msf0_119, msf0_120, msf1_119, msf1_120, msg_179, \
                         msg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * lsg_104[k]
                   + f_1 * msf0_119[k]
                   - f_2 * msf1_119[k]
                   + f_3 * pc_z[k] * msg_179[k];

        t_252[k] = f_17 * lsg_180[k]
                   + f_1 * msf0_120[k]
                   - f_2 * msf1_120[k]
                   + f_3 * pc_x[k] * msg_180[k];

        t_253[k] = f_10 * lsg_120[k]
                   + f_3 * pc_y[k] * msg_180[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_y, pc_z, lsg_105, lsg_122, lsg_183, \
                         msf0_123, msf1_123, msg_180, msg_182, \
                         msg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * lsg_105[k]
                   + f_3 * pc_z[k] * msg_180[k];

        t_255[k] = f_17 * lsg_183[k]
                   + f_6 * msf0_123[k]
                   - f_7 * msf1_123[k]
                   + f_3 * pc_x[k] * msg_183[k];

        t_256[k] = f_10 * lsg_122[k]
                   + f_3 * pc_y[k] * msg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, lsg_108, lsg_185, lsg_186, msf0_125, \
                         msf0_126, msf1_125, msf1_126, msg_183, msg_185, \
                         msg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_17 * lsg_185[k]
                   + f_6 * msf0_125[k]
                   - f_7 * msf1_125[k]
                   + f_3 * pc_x[k] * msg_185[k];

        t_258[k] = f_17 * lsg_186[k]
                   + f_4 * msf0_126[k]
                   - f_5 * msf1_126[k]
                   + f_3 * pc_x[k] * msg_186[k];

        t_259[k] = f_10 * lsg_108[k]
                   + f_3 * pc_z[k] * msg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, lsg_125, lsg_189, lsg_190, \
                         lsg_191, msf0_129, msf1_129, msg_185, msg_189, msg_190, \
                         msg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * lsg_125[k]
                   + f_3 * pc_y[k] * msg_185[k];

        t_261[k] = f_17 * lsg_189[k]
                   + f_4 * msf0_129[k]
                   - f_5 * msf1_129[k]
                   + f_3 * pc_x[k] * msg_189[k];

        t_262[k] = f_17 * lsg_190[k]
                   + f_3 * pc_x[k] * msg_190[k];

        t_263[k] = f_17 * lsg_191[k]
                   + f_3 * pc_x[k] * msg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, lsg_130, lsg_192, lsg_193, \
                         lsg_194, msf0_126, msf1_126, msg_190, msg_192, msg_193, \
                         msg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * lsg_192[k]
                   + f_3 * pc_x[k] * msg_192[k];

        t_265[k] = f_17 * lsg_193[k]
                   + f_3 * pc_x[k] * msg_193[k];

        t_266[k] = f_17 * lsg_194[k]
                   + f_3 * pc_x[k] * msg_194[k];

        t_267[k] = f_10 * lsg_130[k]
                   + f_1 * msf0_126[k]
                   - f_2 * msf1_126[k]
                   + f_3 * pc_y[k] * msg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, lsg_115, lsg_132, lsg_133, msf0_128, \
                         msf0_129, msf1_128, msf1_129, msg_190, msg_192, \
                         msg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * lsg_115[k]
                   + f_3 * pc_z[k] * msg_190[k];

        t_269[k] = f_10 * lsg_132[k]
                   + f_6 * msf0_128[k]
                   - f_7 * msf1_128[k]
                   + f_3 * pc_y[k] * msg_192[k];

        t_270[k] = f_10 * lsg_133[k]
                   + f_4 * msf0_129[k]
                   - f_5 * msf1_129[k]
                   + f_3 * pc_y[k] * msg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_y, pc_y, pc_z, lsh0_189, lsg_119, \
                         lsg_134, lsg_135, lsh1_189, msf0_129, msf1_129, msg_194, \
                         msg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * lsg_134[k]
                   + f_3 * pc_y[k] * msg_194[k];

        t_272[k] = f_10 * lsg_119[k]
                   + f_1 * msf0_129[k]
                   - f_2 * msf1_129[k]
                   + f_3 * pc_z[k] * msg_194[k];

        t_273[k] = pa_y[k] * lsh0_189[k]
                   - f_8 * pc_y[k] * lsh1_189[k];

        t_274[k] = f_9 * lsg_135[k]
                   + f_3 * pc_y[k] * msg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_y, pc_y, pc_z, lsh0_192, lsh0_194, \
                         lsg_120, lsg_136, lsg_137, lsh1_192, lsh1_194, msg_195, \
                         msg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * lsg_120[k]
                   + f_3 * pc_z[k] * msg_195[k];

        t_276[k] = pa_y[k] * lsh0_192[k]
                   + f_10 * lsg_136[k]
                   - f_8 * pc_y[k] * lsh1_192[k];

        t_277[k] = f_9 * lsg_137[k]
                   + f_3 * pc_y[k] * msg_197[k];

        t_278[k] = pa_y[k] * lsh0_194[k]
                   - f_8 * pc_y[k] * lsh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_y, pc_y, pc_z, lsh0_195, lsh0_198, \
                         lsg_123, lsg_138, lsg_140, lsh1_195, lsh1_198, msg_198, \
                         msg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_y[k] * lsh0_195[k]
                   + f_11 * lsg_138[k]
                   - f_8 * pc_y[k] * lsh1_195[k];

        t_280[k] = f_11 * lsg_123[k]
                   + f_3 * pc_z[k] * msg_198[k];

        t_281[k] = f_9 * lsg_140[k]
                   + f_3 * pc_y[k] * msg_200[k];

        t_282[k] = pa_y[k] * lsh0_198[k]
                   - f_8 * pc_y[k] * lsh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, lsg_205, lsg_206, lsg_207, \
                         lsg_208, lsg_209, msg_205, msg_206, msg_207, msg_208, \
                         msg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_17 * lsg_205[k]
                   + f_3 * pc_x[k] * msg_205[k];

        t_284[k] = f_17 * lsg_206[k]
                   + f_3 * pc_x[k] * msg_206[k];

        t_285[k] = f_17 * lsg_207[k]
                   + f_3 * pc_x[k] * msg_207[k];

        t_286[k] = f_17 * lsg_208[k]
                   + f_3 * pc_x[k] * msg_208[k];

        t_287[k] = f_17 * lsg_209[k]
                   + f_3 * pc_x[k] * msg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, lsg_130, lsg_145, lsg_147, msf0_136, \
                         msf0_138, msf1_136, msf1_138, msg_205, \
                         msg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * lsg_145[k]
                   + f_1 * msf0_136[k]
                   - f_2 * msf1_136[k]
                   + f_3 * pc_y[k] * msg_205[k];

        t_289[k] = f_11 * lsg_130[k]
                   + f_3 * pc_z[k] * msg_205[k];

        t_290[k] = f_9 * lsg_147[k]
                   + f_6 * msf0_138[k]
                   - f_7 * msf1_138[k]
                   + f_3 * pc_y[k] * msg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pc_y, lsh0_209, lsg_148, lsg_149, \
                         lsh1_209, msf0_139, msf1_139, msg_208, \
                         msg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * lsg_148[k]
                   + f_4 * msf0_139[k]
                   - f_5 * msf1_139[k]
                   + f_3 * pc_y[k] * msg_208[k];

        t_292[k] = f_9 * lsg_149[k]
                   + f_3 * pc_y[k] * msg_209[k];

        t_293[k] = pa_y[k] * lsh0_209[k]
                   - f_8 * pc_y[k] * lsh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, pc_z, lsg_135, \
                         lsg_210, msf0_140, msf1_140, msg_210, msg_211, \
                         msg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_17 * lsg_210[k]
                   + f_1 * msf0_140[k]
                   - f_2 * msf1_140[k]
                   + f_3 * pc_x[k] * msg_210[k];

        t_295[k] = f_3 * pc_y[k] * msg_210[k];

        t_296[k] = f_18 * lsg_135[k]
                   + f_3 * pc_z[k] * msg_210[k];

        t_297[k] = f_4 * msf0_140[k]
                   - f_5 * msf1_140[k]
                   + f_3 * pc_y[k] * msg_211[k];

        t_298[k] = f_3 * pc_y[k] * msg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, lsg_215, msf0_141, msf0_142, \
                         msf0_145, msf1_141, msf1_142, msf1_145, msg_213, msg_214, \
                         msg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_17 * lsg_215[k]
                   + f_6 * msf0_145[k]
                   - f_7 * msf1_145[k]
                   + f_3 * pc_x[k] * msg_215[k];

        t_300[k] = f_6 * msf0_141[k]
                   - f_7 * msf1_141[k]
                   + f_3 * pc_y[k] * msg_213[k];

        t_301[k] = f_4 * msf0_142[k]
                   - f_5 * msf1_142[k]
                   + f_3 * pc_y[k] * msg_214[k];

        t_302[k] = f_3 * pc_y[k] * msg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, lsg_219, lsg_220, lsg_221, lsg_222, \
                         msf0_149, msf1_149, msg_219, msg_220, msg_221, \
                         msg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_17 * lsg_219[k]
                   + f_4 * msf0_149[k]
                   - f_5 * msf1_149[k]
                   + f_3 * pc_x[k] * msg_219[k];

        t_304[k] = f_17 * lsg_220[k]
                   + f_3 * pc_x[k] * msg_220[k];

        t_305[k] = f_17 * lsg_221[k]
                   + f_3 * pc_x[k] * msg_221[k];

        t_306[k] = f_17 * lsg_222[k]
                   + f_3 * pc_x[k] * msg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_x, pc_y, lsg_224, msf0_146, msf0_147, \
                         msf1_146, msf1_147, msg_219, msg_220, msg_221, \
                         msg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_y[k] * msg_219[k];

        t_308[k] = f_17 * lsg_224[k]
                   + f_3 * pc_x[k] * msg_224[k];

        t_309[k] = f_1 * msf0_146[k]
                   - f_2 * msf1_146[k]
                   + f_3 * pc_y[k] * msg_220[k];

        t_310[k] = f_13 * msf0_147[k]
                   - f_14 * msf1_147[k]
                   + f_3 * pc_y[k] * msg_221[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, lsg_149, msf0_148, msf0_149, \
                         msf1_148, msf1_149, msg_222, msg_223, \
                         msg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_6 * msf0_148[k]
                   - f_7 * msf1_148[k]
                   + f_3 * pc_y[k] * msg_222[k];

        t_312[k] = f_4 * msf0_149[k]
                   - f_5 * msf1_149[k]
                   + f_3 * pc_y[k] * msg_223[k];

        t_313[k] = f_3 * pc_y[k] * msg_224[k];

        t_314[k] = f_18 * lsg_149[k]
                   + f_1 * msf0_149[k]
                   - f_2 * msf1_149[k]
                   + f_3 * pc_z[k] * msg_224[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, lsg_150, lsg_225, \
                         lsg_228, msf0_150, msf0_153, msf1_150, msf1_153, msg_225, \
                         msg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_18 * lsg_225[k]
                   + f_1 * msf0_150[k]
                   - f_2 * msf1_150[k]
                   + f_3 * pc_x[k] * msg_225[k];

        t_316[k] = f_17 * lsg_150[k]
                   + f_3 * pc_y[k] * msg_225[k];

        t_317[k] = f_3 * pc_z[k] * msg_225[k];

        t_318[k] = f_18 * lsg_228[k]
                   + f_6 * msf0_153[k]
                   - f_7 * msf1_153[k]
                   + f_3 * pc_x[k] * msg_228[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pc_x, pc_z, lsg_231, msf0_150, msf0_156, \
                         msf1_150, msf1_156, msg_226, msg_227, msg_228, \
                         msg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_3 * pc_z[k] * msg_226[k];

        t_320[k] = f_4 * msf0_150[k]
                   - f_5 * msf1_150[k]
                   + f_3 * pc_z[k] * msg_227[k];

        t_321[k] = f_18 * lsg_231[k]
                   + f_4 * msf0_156[k]
                   - f_5 * msf1_156[k]
                   + f_3 * pc_x[k] * msg_231[k];

        t_322[k] = f_3 * pc_z[k] * msg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, pc_z, lsg_155, lsg_235, \
                         msf0_152, msf1_152, msg_230, msg_231, \
                         msg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_17 * lsg_155[k]
                   + f_3 * pc_y[k] * msg_230[k];

        t_324[k] = f_6 * msf0_152[k]
                   - f_7 * msf1_152[k]
                   + f_3 * pc_z[k] * msg_230[k];

        t_325[k] = f_18 * lsg_235[k]
                   + f_3 * pc_x[k] * msg_235[k];

        t_326[k] = f_3 * pc_z[k] * msg_231[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, lsg_160, lsg_237, lsg_238, \
                         lsg_239, msf0_156, msf1_156, msg_235, msg_237, msg_238, \
                         msg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_18 * lsg_237[k]
                   + f_3 * pc_x[k] * msg_237[k];

        t_328[k] = f_18 * lsg_238[k]
                   + f_3 * pc_x[k] * msg_238[k];

        t_329[k] = f_18 * lsg_239[k]
                   + f_3 * pc_x[k] * msg_239[k];

        t_330[k] = f_17 * lsg_160[k]
                   + f_1 * msf0_156[k]
                   - f_2 * msf1_156[k]
                   + f_3 * pc_y[k] * msg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pc_y, pc_z, lsg_164, msf0_156, msf0_157, \
                         msf1_156, msf1_157, msg_235, msg_236, msg_237, \
                         msg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * msg_235[k];

        t_332[k] = f_4 * msf0_156[k]
                   - f_5 * msf1_156[k]
                   + f_3 * pc_z[k] * msg_236[k];

        t_333[k] = f_6 * msf0_157[k]
                   - f_7 * msf1_157[k]
                   + f_3 * pc_z[k] * msg_237[k];

        t_334[k] = f_17 * lsg_164[k]
                   + f_3 * pc_y[k] * msg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_z, pc_y, pc_z, lsh0_210, lsg_150, \
                         lsg_165, lsh1_210, msf0_159, msf1_159, msg_239, \
                         msg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_1 * msf0_159[k]
                   - f_2 * msf1_159[k]
                   + f_3 * pc_z[k] * msg_239[k];

        t_336[k] = pa_z[k] * lsh0_210[k]
                   - f_8 * pc_z[k] * lsh1_210[k];

        t_337[k] = f_18 * lsg_165[k]
                   + f_3 * pc_y[k] * msg_240[k];

        t_338[k] = f_9 * lsg_150[k]
                   + f_3 * pc_z[k] * msg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pc_x, pc_y, pc_z, lsh0_213, lsg_167, \
                         lsg_245, lsh1_213, msf0_165, msf1_165, msg_242, \
                         msg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * lsh0_213[k]
                   - f_8 * pc_z[k] * lsh1_213[k];

        t_340[k] = f_18 * lsg_167[k]
                   + f_3 * pc_y[k] * msg_242[k];

        t_341[k] = f_18 * lsg_245[k]
                   + f_6 * msf0_165[k]
                   - f_7 * msf1_165[k]
                   + f_3 * pc_x[k] * msg_245[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_z, pc_y, pc_z, lsh0_216, lsg_153, lsg_170, \
                         lsh1_216, msg_243, msg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * lsh0_216[k]
                   - f_8 * pc_z[k] * lsh1_216[k];

        t_343[k] = f_9 * lsg_153[k]
                   + f_3 * pc_z[k] * msg_243[k];

        t_344[k] = f_18 * lsg_170[k]
                   + f_3 * pc_y[k] * msg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pc_x, lsg_249, lsg_250, lsg_251, lsg_252, \
                         msf0_169, msf1_169, msg_249, msg_250, msg_251, \
                         msg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_18 * lsg_249[k]
                   + f_4 * msf0_169[k]
                   - f_5 * msf1_169[k]
                   + f_3 * pc_x[k] * msg_249[k];

        t_346[k] = f_18 * lsg_250[k]
                   + f_3 * pc_x[k] * msg_250[k];

        t_347[k] = f_18 * lsg_251[k]
                   + f_3 * pc_x[k] * msg_251[k];

        t_348[k] = f_18 * lsg_252[k]
                   + f_3 * pc_x[k] * msg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_z, pc_x, pc_z, lsh0_225, lsg_160, \
                         lsg_253, lsg_254, lsh1_225, msg_250, msg_253, \
                         msg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_18 * lsg_253[k]
                   + f_3 * pc_x[k] * msg_253[k];

        t_350[k] = f_18 * lsg_254[k]
                   + f_3 * pc_x[k] * msg_254[k];

        t_351[k] = pa_z[k] * lsh0_225[k]
                   - f_8 * pc_z[k] * lsh1_225[k];

        t_352[k] = f_9 * lsg_160[k]
                   + f_3 * pc_z[k] * msg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, lsg_177, lsg_178, lsg_179, msf0_168, \
                         msf0_169, msf1_168, msf1_169, msg_252, msg_253, \
                         msg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_18 * lsg_177[k]
                   + f_6 * msf0_168[k]
                   - f_7 * msf1_168[k]
                   + f_3 * pc_y[k] * msg_252[k];

        t_354[k] = f_18 * lsg_178[k]
                   + f_4 * msf0_169[k]
                   - f_5 * msf1_169[k]
                   + f_3 * pc_y[k] * msg_253[k];

        t_355[k] = f_18 * lsg_179[k]
                   + f_3 * pc_y[k] * msg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, lsg_164, lsg_180, lsg_255, \
                         msf0_169, msf0_170, msf1_169, msf1_170, msg_254, \
                         msg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_9 * lsg_164[k]
                   + f_1 * msf0_169[k]
                   - f_2 * msf1_169[k]
                   + f_3 * pc_z[k] * msg_254[k];

        t_357[k] = f_18 * lsg_255[k]
                   + f_1 * msf0_170[k]
                   - f_2 * msf1_170[k]
                   + f_3 * pc_x[k] * msg_255[k];

        t_358[k] = f_11 * lsg_180[k]
                   + f_3 * pc_y[k] * msg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, pc_y, pc_z, lsg_165, lsg_182, lsg_258, \
                         msf0_173, msf1_173, msg_255, msg_257, \
                         msg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_10 * lsg_165[k]
                   + f_3 * pc_z[k] * msg_255[k];

        t_360[k] = f_18 * lsg_258[k]
                   + f_6 * msf0_173[k]
                   - f_7 * msf1_173[k]
                   + f_3 * pc_x[k] * msg_258[k];

        t_361[k] = f_11 * lsg_182[k]
                   + f_3 * pc_y[k] * msg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, lsg_168, lsg_260, lsg_261, msf0_175, \
                         msf0_176, msf1_175, msf1_176, msg_258, msg_260, \
                         msg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_18 * lsg_260[k]
                   + f_6 * msf0_175[k]
                   - f_7 * msf1_175[k]
                   + f_3 * pc_x[k] * msg_260[k];

        t_363[k] = f_18 * lsg_261[k]
                   + f_4 * msf0_176[k]
                   - f_5 * msf1_176[k]
                   + f_3 * pc_x[k] * msg_261[k];

        t_364[k] = f_10 * lsg_168[k]
                   + f_3 * pc_z[k] * msg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, lsg_185, lsg_264, lsg_265, \
                         lsg_266, msf0_179, msf1_179, msg_260, msg_264, msg_265, \
                         msg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * lsg_185[k]
                   + f_3 * pc_y[k] * msg_260[k];

        t_366[k] = f_18 * lsg_264[k]
                   + f_4 * msf0_179[k]
                   - f_5 * msf1_179[k]
                   + f_3 * pc_x[k] * msg_264[k];

        t_367[k] = f_18 * lsg_265[k]
                   + f_3 * pc_x[k] * msg_265[k];

        t_368[k] = f_18 * lsg_266[k]
                   + f_3 * pc_x[k] * msg_266[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_294 = buffer.data(lsh0 + 294);
    const auto *lsh0_297 = buffer.data(lsh0 + 297);
    const auto *lsh0_299 = buffer.data(lsh0 + 299);
    const auto *lsh0_300 = buffer.data(lsh0 + 300);
    const auto *lsh0_303 = buffer.data(lsh0 + 303);
    const auto *lsh0_314 = buffer.data(lsh0 + 314);
    const auto *lsh0_315 = buffer.data(lsh0 + 315);
    const auto *lsh0_318 = buffer.data(lsh0 + 318);
    const auto *lsh0_321 = buffer.data(lsh0 + 321);
    const auto *lsh0_330 = buffer.data(lsh0 + 330);

    const auto *lsg_175 = buffer.data(lsg + 175);
    const auto *lsg_179 = buffer.data(lsg + 179);
    const auto *lsg_180 = buffer.data(lsg + 180);
    const auto *lsg_183 = buffer.data(lsg + 183);
    const auto *lsg_190 = buffer.data(lsg + 190);
    const auto *lsg_192 = buffer.data(lsg + 192);
    const auto *lsg_193 = buffer.data(lsg + 193);
    const auto *lsg_194 = buffer.data(lsg + 194);
    const auto *lsg_195 = buffer.data(lsg + 195);
    const auto *lsg_197 = buffer.data(lsg + 197);
    const auto *lsg_198 = buffer.data(lsg + 198);
    const auto *lsg_200 = buffer.data(lsg + 200);
    const auto *lsg_205 = buffer.data(lsg + 205);
    const auto *lsg_207 = buffer.data(lsg + 207);
    const auto *lsg_208 = buffer.data(lsg + 208);
    const auto *lsg_209 = buffer.data(lsg + 209);
    const auto *lsg_210 = buffer.data(lsg + 210);
    const auto *lsg_211 = buffer.data(lsg + 211);
    const auto *lsg_212 = buffer.data(lsg + 212);
    const auto *lsg_213 = buffer.data(lsg + 213);
    const auto *lsg_215 = buffer.data(lsg + 215);
    const auto *lsg_220 = buffer.data(lsg + 220);
    const auto *lsg_222 = buffer.data(lsg + 222);
    const auto *lsg_223 = buffer.data(lsg + 223);
    const auto *lsg_224 = buffer.data(lsg + 224);
    const auto *lsg_225 = buffer.data(lsg + 225);
    const auto *lsg_228 = buffer.data(lsg + 228);
    const auto *lsg_230 = buffer.data(lsg + 230);
    const auto *lsg_235 = buffer.data(lsg + 235);
    const auto *lsg_239 = buffer.data(lsg + 239);
    const auto *lsg_240 = buffer.data(lsg + 240);
    const auto *lsg_242 = buffer.data(lsg + 242);
    const auto *lsg_245 = buffer.data(lsg + 245);
    const auto *lsg_252 = buffer.data(lsg + 252);
    const auto *lsg_253 = buffer.data(lsg + 253);
    const auto *lsg_254 = buffer.data(lsg + 254);
    const auto *lsg_255 = buffer.data(lsg + 255);
    const auto *lsg_257 = buffer.data(lsg + 257);
    const auto *lsg_267 = buffer.data(lsg + 267);
    const auto *lsg_268 = buffer.data(lsg + 268);
    const auto *lsg_269 = buffer.data(lsg + 269);
    const auto *lsg_270 = buffer.data(lsg + 270);
    const auto *lsg_273 = buffer.data(lsg + 273);
    const auto *lsg_275 = buffer.data(lsg + 275);
    const auto *lsg_276 = buffer.data(lsg + 276);
    const auto *lsg_279 = buffer.data(lsg + 279);
    const auto *lsg_280 = buffer.data(lsg + 280);
    const auto *lsg_281 = buffer.data(lsg + 281);
    const auto *lsg_282 = buffer.data(lsg + 282);
    const auto *lsg_283 = buffer.data(lsg + 283);
    const auto *lsg_284 = buffer.data(lsg + 284);
    const auto *lsg_295 = buffer.data(lsg + 295);
    const auto *lsg_296 = buffer.data(lsg + 296);
    const auto *lsg_297 = buffer.data(lsg + 297);
    const auto *lsg_298 = buffer.data(lsg + 298);
    const auto *lsg_299 = buffer.data(lsg + 299);
    const auto *lsg_300 = buffer.data(lsg + 300);
    const auto *lsg_305 = buffer.data(lsg + 305);
    const auto *lsg_309 = buffer.data(lsg + 309);
    const auto *lsg_310 = buffer.data(lsg + 310);
    const auto *lsg_311 = buffer.data(lsg + 311);
    const auto *lsg_312 = buffer.data(lsg + 312);
    const auto *lsg_314 = buffer.data(lsg + 314);
    const auto *lsg_315 = buffer.data(lsg + 315);
    const auto *lsg_318 = buffer.data(lsg + 318);
    const auto *lsg_321 = buffer.data(lsg + 321);
    const auto *lsg_325 = buffer.data(lsg + 325);
    const auto *lsg_327 = buffer.data(lsg + 327);
    const auto *lsg_328 = buffer.data(lsg + 328);
    const auto *lsg_329 = buffer.data(lsg + 329);
    const auto *lsg_335 = buffer.data(lsg + 335);
    const auto *lsg_339 = buffer.data(lsg + 339);
    const auto *lsg_340 = buffer.data(lsg + 340);
    const auto *lsg_341 = buffer.data(lsg + 341);
    const auto *lsg_342 = buffer.data(lsg + 342);
    const auto *lsg_343 = buffer.data(lsg + 343);
    const auto *lsg_344 = buffer.data(lsg + 344);
    const auto *lsg_345 = buffer.data(lsg + 345);
    const auto *lsg_348 = buffer.data(lsg + 348);

    const auto *lsh1_294 = buffer.data(lsh1 + 294);
    const auto *lsh1_297 = buffer.data(lsh1 + 297);
    const auto *lsh1_299 = buffer.data(lsh1 + 299);
    const auto *lsh1_300 = buffer.data(lsh1 + 300);
    const auto *lsh1_303 = buffer.data(lsh1 + 303);
    const auto *lsh1_314 = buffer.data(lsh1 + 314);
    const auto *lsh1_315 = buffer.data(lsh1 + 315);
    const auto *lsh1_318 = buffer.data(lsh1 + 318);
    const auto *lsh1_321 = buffer.data(lsh1 + 321);
    const auto *lsh1_330 = buffer.data(lsh1 + 330);

    const auto *msf0_176 = buffer.data(msf0 + 176);
    const auto *msf0_178 = buffer.data(msf0 + 178);
    const auto *msf0_179 = buffer.data(msf0 + 179);
    const auto *msf0_180 = buffer.data(msf0 + 180);
    const auto *msf0_183 = buffer.data(msf0 + 183);
    const auto *msf0_185 = buffer.data(msf0 + 185);
    const auto *msf0_186 = buffer.data(msf0 + 186);
    const auto *msf0_188 = buffer.data(msf0 + 188);
    const auto *msf0_189 = buffer.data(msf0 + 189);
    const auto *msf0_196 = buffer.data(msf0 + 196);
    const auto *msf0_198 = buffer.data(msf0 + 198);
    const auto *msf0_199 = buffer.data(msf0 + 199);
    const auto *msf0_200 = buffer.data(msf0 + 200);
    const auto *msf0_201 = buffer.data(msf0 + 201);
    const auto *msf0_202 = buffer.data(msf0 + 202);
    const auto *msf0_205 = buffer.data(msf0 + 205);
    const auto *msf0_206 = buffer.data(msf0 + 206);
    const auto *msf0_207 = buffer.data(msf0 + 207);
    const auto *msf0_208 = buffer.data(msf0 + 208);
    const auto *msf0_209 = buffer.data(msf0 + 209);
    const auto *msf0_210 = buffer.data(msf0 + 210);
    const auto *msf0_212 = buffer.data(msf0 + 212);
    const auto *msf0_213 = buffer.data(msf0 + 213);
    const auto *msf0_216 = buffer.data(msf0 + 216);
    const auto *msf0_217 = buffer.data(msf0 + 217);
    const auto *msf0_219 = buffer.data(msf0 + 219);
    const auto *msf0_225 = buffer.data(msf0 + 225);
    const auto *msf0_228 = buffer.data(msf0 + 228);
    const auto *msf0_229 = buffer.data(msf0 + 229);
    const auto *msf0_230 = buffer.data(msf0 + 230);
    const auto *msf0_233 = buffer.data(msf0 + 233);

    const auto *msf1_176 = buffer.data(msf1 + 176);
    const auto *msf1_178 = buffer.data(msf1 + 178);
    const auto *msf1_179 = buffer.data(msf1 + 179);
    const auto *msf1_180 = buffer.data(msf1 + 180);
    const auto *msf1_183 = buffer.data(msf1 + 183);
    const auto *msf1_185 = buffer.data(msf1 + 185);
    const auto *msf1_186 = buffer.data(msf1 + 186);
    const auto *msf1_188 = buffer.data(msf1 + 188);
    const auto *msf1_189 = buffer.data(msf1 + 189);
    const auto *msf1_196 = buffer.data(msf1 + 196);
    const auto *msf1_198 = buffer.data(msf1 + 198);
    const auto *msf1_199 = buffer.data(msf1 + 199);
    const auto *msf1_200 = buffer.data(msf1 + 200);
    const auto *msf1_201 = buffer.data(msf1 + 201);
    const auto *msf1_202 = buffer.data(msf1 + 202);
    const auto *msf1_205 = buffer.data(msf1 + 205);
    const auto *msf1_206 = buffer.data(msf1 + 206);
    const auto *msf1_207 = buffer.data(msf1 + 207);
    const auto *msf1_208 = buffer.data(msf1 + 208);
    const auto *msf1_209 = buffer.data(msf1 + 209);
    const auto *msf1_210 = buffer.data(msf1 + 210);
    const auto *msf1_212 = buffer.data(msf1 + 212);
    const auto *msf1_213 = buffer.data(msf1 + 213);
    const auto *msf1_216 = buffer.data(msf1 + 216);
    const auto *msf1_217 = buffer.data(msf1 + 217);
    const auto *msf1_219 = buffer.data(msf1 + 219);
    const auto *msf1_225 = buffer.data(msf1 + 225);
    const auto *msf1_228 = buffer.data(msf1 + 228);
    const auto *msf1_229 = buffer.data(msf1 + 229);
    const auto *msf1_230 = buffer.data(msf1 + 230);
    const auto *msf1_233 = buffer.data(msf1 + 233);

    const auto *msg_265 = buffer.data(msg + 265);
    const auto *msg_267 = buffer.data(msg + 267);
    const auto *msg_268 = buffer.data(msg + 268);
    const auto *msg_269 = buffer.data(msg + 269);
    const auto *msg_270 = buffer.data(msg + 270);
    const auto *msg_272 = buffer.data(msg + 272);
    const auto *msg_273 = buffer.data(msg + 273);
    const auto *msg_275 = buffer.data(msg + 275);
    const auto *msg_276 = buffer.data(msg + 276);
    const auto *msg_279 = buffer.data(msg + 279);
    const auto *msg_280 = buffer.data(msg + 280);
    const auto *msg_281 = buffer.data(msg + 281);
    const auto *msg_282 = buffer.data(msg + 282);
    const auto *msg_283 = buffer.data(msg + 283);
    const auto *msg_284 = buffer.data(msg + 284);
    const auto *msg_285 = buffer.data(msg + 285);
    const auto *msg_287 = buffer.data(msg + 287);
    const auto *msg_288 = buffer.data(msg + 288);
    const auto *msg_290 = buffer.data(msg + 290);
    const auto *msg_295 = buffer.data(msg + 295);
    const auto *msg_296 = buffer.data(msg + 296);
    const auto *msg_297 = buffer.data(msg + 297);
    const auto *msg_298 = buffer.data(msg + 298);
    const auto *msg_299 = buffer.data(msg + 299);
    const auto *msg_300 = buffer.data(msg + 300);
    const auto *msg_301 = buffer.data(msg + 301);
    const auto *msg_302 = buffer.data(msg + 302);
    const auto *msg_303 = buffer.data(msg + 303);
    const auto *msg_304 = buffer.data(msg + 304);
    const auto *msg_305 = buffer.data(msg + 305);
    const auto *msg_309 = buffer.data(msg + 309);
    const auto *msg_310 = buffer.data(msg + 310);
    const auto *msg_311 = buffer.data(msg + 311);
    const auto *msg_312 = buffer.data(msg + 312);
    const auto *msg_313 = buffer.data(msg + 313);
    const auto *msg_314 = buffer.data(msg + 314);
    const auto *msg_315 = buffer.data(msg + 315);
    const auto *msg_316 = buffer.data(msg + 316);
    const auto *msg_317 = buffer.data(msg + 317);
    const auto *msg_318 = buffer.data(msg + 318);
    const auto *msg_320 = buffer.data(msg + 320);
    const auto *msg_321 = buffer.data(msg + 321);
    const auto *msg_325 = buffer.data(msg + 325);
    const auto *msg_326 = buffer.data(msg + 326);
    const auto *msg_327 = buffer.data(msg + 327);
    const auto *msg_328 = buffer.data(msg + 328);
    const auto *msg_329 = buffer.data(msg + 329);
    const auto *msg_330 = buffer.data(msg + 330);
    const auto *msg_332 = buffer.data(msg + 332);
    const auto *msg_333 = buffer.data(msg + 333);
    const auto *msg_335 = buffer.data(msg + 335);
    const auto *msg_339 = buffer.data(msg + 339);
    const auto *msg_340 = buffer.data(msg + 340);
    const auto *msg_341 = buffer.data(msg + 341);
    const auto *msg_342 = buffer.data(msg + 342);
    const auto *msg_343 = buffer.data(msg + 343);
    const auto *msg_344 = buffer.data(msg + 344);
    const auto *msg_345 = buffer.data(msg + 345);
    const auto *msg_347 = buffer.data(msg + 347);
    const auto *msg_348 = buffer.data(msg + 348);

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, lsg_190, lsg_267, lsg_268, \
                         lsg_269, msf0_176, msf1_176, msg_265, msg_267, msg_268, \
                         msg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_18 * lsg_267[k]
                   + f_3 * pc_x[k] * msg_267[k];

        t_370[k] = f_18 * lsg_268[k]
                   + f_3 * pc_x[k] * msg_268[k];

        t_371[k] = f_18 * lsg_269[k]
                   + f_3 * pc_x[k] * msg_269[k];

        t_372[k] = f_11 * lsg_190[k]
                   + f_1 * msf0_176[k]
                   - f_2 * msf1_176[k]
                   + f_3 * pc_y[k] * msg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, lsg_175, lsg_192, lsg_193, msf0_178, \
                         msf0_179, msf1_178, msf1_179, msg_265, msg_267, \
                         msg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * lsg_175[k]
                   + f_3 * pc_z[k] * msg_265[k];

        t_374[k] = f_11 * lsg_192[k]
                   + f_6 * msf0_178[k]
                   - f_7 * msf1_178[k]
                   + f_3 * pc_y[k] * msg_267[k];

        t_375[k] = f_11 * lsg_193[k]
                   + f_4 * msf0_179[k]
                   - f_5 * msf1_179[k]
                   + f_3 * pc_y[k] * msg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, lsg_179, lsg_194, lsg_270, \
                         msf0_179, msf0_180, msf1_179, msf1_180, msg_269, \
                         msg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * lsg_194[k]
                   + f_3 * pc_y[k] * msg_269[k];

        t_377[k] = f_10 * lsg_179[k]
                   + f_1 * msf0_179[k]
                   - f_2 * msf1_179[k]
                   + f_3 * pc_z[k] * msg_269[k];

        t_378[k] = f_18 * lsg_270[k]
                   + f_1 * msf0_180[k]
                   - f_2 * msf1_180[k]
                   + f_3 * pc_x[k] * msg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, lsg_180, lsg_195, \
                         lsg_197, lsg_273, msf0_183, msf1_183, msg_270, msg_272, \
                         msg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * lsg_195[k]
                   + f_3 * pc_y[k] * msg_270[k];

        t_380[k] = f_11 * lsg_180[k]
                   + f_3 * pc_z[k] * msg_270[k];

        t_381[k] = f_18 * lsg_273[k]
                   + f_6 * msf0_183[k]
                   - f_7 * msf1_183[k]
                   + f_3 * pc_x[k] * msg_273[k];

        t_382[k] = f_10 * lsg_197[k]
                   + f_3 * pc_y[k] * msg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, lsg_183, lsg_275, lsg_276, msf0_185, \
                         msf0_186, msf1_185, msf1_186, msg_273, msg_275, \
                         msg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_18 * lsg_275[k]
                   + f_6 * msf0_185[k]
                   - f_7 * msf1_185[k]
                   + f_3 * pc_x[k] * msg_275[k];

        t_384[k] = f_18 * lsg_276[k]
                   + f_4 * msf0_186[k]
                   - f_5 * msf1_186[k]
                   + f_3 * pc_x[k] * msg_276[k];

        t_385[k] = f_11 * lsg_183[k]
                   + f_3 * pc_z[k] * msg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, lsg_200, lsg_279, lsg_280, \
                         lsg_281, msf0_189, msf1_189, msg_275, msg_279, msg_280, \
                         msg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * lsg_200[k]
                   + f_3 * pc_y[k] * msg_275[k];

        t_387[k] = f_18 * lsg_279[k]
                   + f_4 * msf0_189[k]
                   - f_5 * msf1_189[k]
                   + f_3 * pc_x[k] * msg_279[k];

        t_388[k] = f_18 * lsg_280[k]
                   + f_3 * pc_x[k] * msg_280[k];

        t_389[k] = f_18 * lsg_281[k]
                   + f_3 * pc_x[k] * msg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, lsg_205, lsg_282, lsg_283, \
                         lsg_284, msf0_186, msf1_186, msg_280, msg_282, msg_283, \
                         msg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_18 * lsg_282[k]
                   + f_3 * pc_x[k] * msg_282[k];

        t_391[k] = f_18 * lsg_283[k]
                   + f_3 * pc_x[k] * msg_283[k];

        t_392[k] = f_18 * lsg_284[k]
                   + f_3 * pc_x[k] * msg_284[k];

        t_393[k] = f_10 * lsg_205[k]
                   + f_1 * msf0_186[k]
                   - f_2 * msf1_186[k]
                   + f_3 * pc_y[k] * msg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, lsg_190, lsg_207, lsg_208, msf0_188, \
                         msf0_189, msf1_188, msf1_189, msg_280, msg_282, \
                         msg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * lsg_190[k]
                   + f_3 * pc_z[k] * msg_280[k];

        t_395[k] = f_10 * lsg_207[k]
                   + f_6 * msf0_188[k]
                   - f_7 * msf1_188[k]
                   + f_3 * pc_y[k] * msg_282[k];

        t_396[k] = f_10 * lsg_208[k]
                   + f_4 * msf0_189[k]
                   - f_5 * msf1_189[k]
                   + f_3 * pc_y[k] * msg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pc_y, pc_z, lsh0_294, lsg_194, \
                         lsg_209, lsg_210, lsh1_294, msf0_189, msf1_189, msg_284, \
                         msg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * lsg_209[k]
                   + f_3 * pc_y[k] * msg_284[k];

        t_398[k] = f_11 * lsg_194[k]
                   + f_1 * msf0_189[k]
                   - f_2 * msf1_189[k]
                   + f_3 * pc_z[k] * msg_284[k];

        t_399[k] = pa_y[k] * lsh0_294[k]
                   - f_8 * pc_y[k] * lsh1_294[k];

        t_400[k] = f_9 * lsg_210[k]
                   + f_3 * pc_y[k] * msg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pc_y, pc_z, lsh0_297, lsh0_299, \
                         lsg_195, lsg_211, lsg_212, lsh1_297, lsh1_299, msg_285, \
                         msg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_18 * lsg_195[k]
                   + f_3 * pc_z[k] * msg_285[k];

        t_402[k] = pa_y[k] * lsh0_297[k]
                   + f_10 * lsg_211[k]
                   - f_8 * pc_y[k] * lsh1_297[k];

        t_403[k] = f_9 * lsg_212[k]
                   + f_3 * pc_y[k] * msg_287[k];

        t_404[k] = pa_y[k] * lsh0_299[k]
                   - f_8 * pc_y[k] * lsh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_y, pc_y, pc_z, lsh0_300, lsh0_303, \
                         lsg_198, lsg_213, lsg_215, lsh1_300, lsh1_303, msg_288, \
                         msg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_y[k] * lsh0_300[k]
                   + f_11 * lsg_213[k]
                   - f_8 * pc_y[k] * lsh1_300[k];

        t_406[k] = f_18 * lsg_198[k]
                   + f_3 * pc_z[k] * msg_288[k];

        t_407[k] = f_9 * lsg_215[k]
                   + f_3 * pc_y[k] * msg_290[k];

        t_408[k] = pa_y[k] * lsh0_303[k]
                   - f_8 * pc_y[k] * lsh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, lsg_295, lsg_296, lsg_297, \
                         lsg_298, lsg_299, msg_295, msg_296, msg_297, msg_298, \
                         msg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_18 * lsg_295[k]
                   + f_3 * pc_x[k] * msg_295[k];

        t_410[k] = f_18 * lsg_296[k]
                   + f_3 * pc_x[k] * msg_296[k];

        t_411[k] = f_18 * lsg_297[k]
                   + f_3 * pc_x[k] * msg_297[k];

        t_412[k] = f_18 * lsg_298[k]
                   + f_3 * pc_x[k] * msg_298[k];

        t_413[k] = f_18 * lsg_299[k]
                   + f_3 * pc_x[k] * msg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, lsg_205, lsg_220, lsg_222, msf0_196, \
                         msf0_198, msf1_196, msf1_198, msg_295, \
                         msg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * lsg_220[k]
                   + f_1 * msf0_196[k]
                   - f_2 * msf1_196[k]
                   + f_3 * pc_y[k] * msg_295[k];

        t_415[k] = f_18 * lsg_205[k]
                   + f_3 * pc_z[k] * msg_295[k];

        t_416[k] = f_9 * lsg_222[k]
                   + f_6 * msf0_198[k]
                   - f_7 * msf1_198[k]
                   + f_3 * pc_y[k] * msg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, lsh0_314, lsg_223, lsg_224, \
                         lsh1_314, msf0_199, msf1_199, msg_298, \
                         msg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * lsg_223[k]
                   + f_4 * msf0_199[k]
                   - f_5 * msf1_199[k]
                   + f_3 * pc_y[k] * msg_298[k];

        t_418[k] = f_9 * lsg_224[k]
                   + f_3 * pc_y[k] * msg_299[k];

        t_419[k] = pa_y[k] * lsh0_314[k]
                   - f_8 * pc_y[k] * lsh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, lsg_210, \
                         lsg_300, msf0_200, msf1_200, msg_300, msg_301, \
                         msg_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_18 * lsg_300[k]
                   + f_1 * msf0_200[k]
                   - f_2 * msf1_200[k]
                   + f_3 * pc_x[k] * msg_300[k];

        t_421[k] = f_3 * pc_y[k] * msg_300[k];

        t_422[k] = f_17 * lsg_210[k]
                   + f_3 * pc_z[k] * msg_300[k];

        t_423[k] = f_4 * msf0_200[k]
                   - f_5 * msf1_200[k]
                   + f_3 * pc_y[k] * msg_301[k];

        t_424[k] = f_3 * pc_y[k] * msg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, lsg_305, msf0_201, msf0_202, \
                         msf0_205, msf1_201, msf1_202, msf1_205, msg_303, msg_304, \
                         msg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_18 * lsg_305[k]
                   + f_6 * msf0_205[k]
                   - f_7 * msf1_205[k]
                   + f_3 * pc_x[k] * msg_305[k];

        t_426[k] = f_6 * msf0_201[k]
                   - f_7 * msf1_201[k]
                   + f_3 * pc_y[k] * msg_303[k];

        t_427[k] = f_4 * msf0_202[k]
                   - f_5 * msf1_202[k]
                   + f_3 * pc_y[k] * msg_304[k];

        t_428[k] = f_3 * pc_y[k] * msg_305[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pc_x, lsg_309, lsg_310, lsg_311, lsg_312, \
                         msf0_209, msf1_209, msg_309, msg_310, msg_311, \
                         msg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_18 * lsg_309[k]
                   + f_4 * msf0_209[k]
                   - f_5 * msf1_209[k]
                   + f_3 * pc_x[k] * msg_309[k];

        t_430[k] = f_18 * lsg_310[k]
                   + f_3 * pc_x[k] * msg_310[k];

        t_431[k] = f_18 * lsg_311[k]
                   + f_3 * pc_x[k] * msg_311[k];

        t_432[k] = f_18 * lsg_312[k]
                   + f_3 * pc_x[k] * msg_312[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, lsg_314, msf0_206, msf0_207, \
                         msf1_206, msf1_207, msg_309, msg_310, msg_311, \
                         msg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_3 * pc_y[k] * msg_309[k];

        t_434[k] = f_18 * lsg_314[k]
                   + f_3 * pc_x[k] * msg_314[k];

        t_435[k] = f_1 * msf0_206[k]
                   - f_2 * msf1_206[k]
                   + f_3 * pc_y[k] * msg_310[k];

        t_436[k] = f_13 * msf0_207[k]
                   - f_14 * msf1_207[k]
                   + f_3 * pc_y[k] * msg_311[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_y, pc_z, lsg_224, msf0_208, msf0_209, \
                         msf1_208, msf1_209, msg_312, msg_313, \
                         msg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_6 * msf0_208[k]
                   - f_7 * msf1_208[k]
                   + f_3 * pc_y[k] * msg_312[k];

        t_438[k] = f_4 * msf0_209[k]
                   - f_5 * msf1_209[k]
                   + f_3 * pc_y[k] * msg_313[k];

        t_439[k] = f_3 * pc_y[k] * msg_314[k];

        t_440[k] = f_17 * lsg_224[k]
                   + f_1 * msf0_209[k]
                   - f_2 * msf1_209[k]
                   + f_3 * pc_z[k] * msg_314[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, pc_y, pc_z, lsg_225, lsg_315, \
                         lsg_318, msf0_210, msf0_213, msf1_210, msf1_213, msg_315, \
                         msg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_11 * lsg_315[k]
                   + f_1 * msf0_210[k]
                   - f_2 * msf1_210[k]
                   + f_3 * pc_x[k] * msg_315[k];

        t_442[k] = f_16 * lsg_225[k]
                   + f_3 * pc_y[k] * msg_315[k];

        t_443[k] = f_3 * pc_z[k] * msg_315[k];

        t_444[k] = f_11 * lsg_318[k]
                   + f_6 * msf0_213[k]
                   - f_7 * msf1_213[k]
                   + f_3 * pc_x[k] * msg_318[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_z, lsg_321, msf0_210, msf0_216, \
                         msf1_210, msf1_216, msg_316, msg_317, msg_318, \
                         msg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_3 * pc_z[k] * msg_316[k];

        t_446[k] = f_4 * msf0_210[k]
                   - f_5 * msf1_210[k]
                   + f_3 * pc_z[k] * msg_317[k];

        t_447[k] = f_11 * lsg_321[k]
                   + f_4 * msf0_216[k]
                   - f_5 * msf1_216[k]
                   + f_3 * pc_x[k] * msg_321[k];

        t_448[k] = f_3 * pc_z[k] * msg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, lsg_230, lsg_325, \
                         msf0_212, msf1_212, msg_320, msg_321, \
                         msg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_16 * lsg_230[k]
                   + f_3 * pc_y[k] * msg_320[k];

        t_450[k] = f_6 * msf0_212[k]
                   - f_7 * msf1_212[k]
                   + f_3 * pc_z[k] * msg_320[k];

        t_451[k] = f_11 * lsg_325[k]
                   + f_3 * pc_x[k] * msg_325[k];

        t_452[k] = f_3 * pc_z[k] * msg_321[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, lsg_235, lsg_327, lsg_328, \
                         lsg_329, msf0_216, msf1_216, msg_325, msg_327, msg_328, \
                         msg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_11 * lsg_327[k]
                   + f_3 * pc_x[k] * msg_327[k];

        t_454[k] = f_11 * lsg_328[k]
                   + f_3 * pc_x[k] * msg_328[k];

        t_455[k] = f_11 * lsg_329[k]
                   + f_3 * pc_x[k] * msg_329[k];

        t_456[k] = f_16 * lsg_235[k]
                   + f_1 * msf0_216[k]
                   - f_2 * msf1_216[k]
                   + f_3 * pc_y[k] * msg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pc_y, pc_z, lsg_239, msf0_216, msf0_217, \
                         msf1_216, msf1_217, msg_325, msg_326, msg_327, \
                         msg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * msg_325[k];

        t_458[k] = f_4 * msf0_216[k]
                   - f_5 * msf1_216[k]
                   + f_3 * pc_z[k] * msg_326[k];

        t_459[k] = f_6 * msf0_217[k]
                   - f_7 * msf1_217[k]
                   + f_3 * pc_z[k] * msg_327[k];

        t_460[k] = f_16 * lsg_239[k]
                   + f_3 * pc_y[k] * msg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pa_z, pc_y, pc_z, lsh0_315, lsg_225, \
                         lsg_240, lsh1_315, msf0_219, msf1_219, msg_329, \
                         msg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_1 * msf0_219[k]
                   - f_2 * msf1_219[k]
                   + f_3 * pc_z[k] * msg_329[k];

        t_462[k] = pa_z[k] * lsh0_315[k]
                   - f_8 * pc_z[k] * lsh1_315[k];

        t_463[k] = f_17 * lsg_240[k]
                   + f_3 * pc_y[k] * msg_330[k];

        t_464[k] = f_9 * lsg_225[k]
                   + f_3 * pc_z[k] * msg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pa_z, pc_x, pc_y, pc_z, lsh0_318, lsg_242, \
                         lsg_335, lsh1_318, msf0_225, msf1_225, msg_332, \
                         msg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * lsh0_318[k]
                   - f_8 * pc_z[k] * lsh1_318[k];

        t_466[k] = f_17 * lsg_242[k]
                   + f_3 * pc_y[k] * msg_332[k];

        t_467[k] = f_11 * lsg_335[k]
                   + f_6 * msf0_225[k]
                   - f_7 * msf1_225[k]
                   + f_3 * pc_x[k] * msg_335[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pa_z, pc_y, pc_z, lsh0_321, lsg_228, lsg_245, \
                         lsh1_321, msg_333, msg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pa_z[k] * lsh0_321[k]
                   - f_8 * pc_z[k] * lsh1_321[k];

        t_469[k] = f_9 * lsg_228[k]
                   + f_3 * pc_z[k] * msg_333[k];

        t_470[k] = f_17 * lsg_245[k]
                   + f_3 * pc_y[k] * msg_335[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, lsg_339, lsg_340, lsg_341, lsg_342, \
                         msf0_229, msf1_229, msg_339, msg_340, msg_341, \
                         msg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_11 * lsg_339[k]
                   + f_4 * msf0_229[k]
                   - f_5 * msf1_229[k]
                   + f_3 * pc_x[k] * msg_339[k];

        t_472[k] = f_11 * lsg_340[k]
                   + f_3 * pc_x[k] * msg_340[k];

        t_473[k] = f_11 * lsg_341[k]
                   + f_3 * pc_x[k] * msg_341[k];

        t_474[k] = f_11 * lsg_342[k]
                   + f_3 * pc_x[k] * msg_342[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_z, pc_x, pc_z, lsh0_330, lsg_235, \
                         lsg_343, lsg_344, lsh1_330, msg_340, msg_343, \
                         msg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_11 * lsg_343[k]
                   + f_3 * pc_x[k] * msg_343[k];

        t_476[k] = f_11 * lsg_344[k]
                   + f_3 * pc_x[k] * msg_344[k];

        t_477[k] = pa_z[k] * lsh0_330[k]
                   - f_8 * pc_z[k] * lsh1_330[k];

        t_478[k] = f_9 * lsg_235[k]
                   + f_3 * pc_z[k] * msg_340[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_y, lsg_252, lsg_253, lsg_254, msf0_228, \
                         msf0_229, msf1_228, msf1_229, msg_342, msg_343, \
                         msg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_17 * lsg_252[k]
                   + f_6 * msf0_228[k]
                   - f_7 * msf1_228[k]
                   + f_3 * pc_y[k] * msg_342[k];

        t_480[k] = f_17 * lsg_253[k]
                   + f_4 * msf0_229[k]
                   - f_5 * msf1_229[k]
                   + f_3 * pc_y[k] * msg_343[k];

        t_481[k] = f_17 * lsg_254[k]
                   + f_3 * pc_y[k] * msg_344[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, lsg_239, lsg_255, lsg_345, \
                         msf0_229, msf0_230, msf1_229, msf1_230, msg_344, \
                         msg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_9 * lsg_239[k]
                   + f_1 * msf0_229[k]
                   - f_2 * msf1_229[k]
                   + f_3 * pc_z[k] * msg_344[k];

        t_483[k] = f_11 * lsg_345[k]
                   + f_1 * msf0_230[k]
                   - f_2 * msf1_230[k]
                   + f_3 * pc_x[k] * msg_345[k];

        t_484[k] = f_18 * lsg_255[k]
                   + f_3 * pc_y[k] * msg_345[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, pc_z, lsg_240, lsg_257, lsg_348, \
                         msf0_233, msf1_233, msg_345, msg_347, \
                         msg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_10 * lsg_240[k]
                   + f_3 * pc_z[k] * msg_345[k];

        t_486[k] = f_11 * lsg_348[k]
                   + f_6 * msf0_233[k]
                   - f_7 * msf1_233[k]
                   + f_3 * pc_x[k] * msg_348[k];

        t_487[k] = f_18 * lsg_257[k]
                   + f_3 * pc_y[k] * msg_347[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;

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
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_420 = buffer.data(lsh0 + 420);
    const auto *lsh0_423 = buffer.data(lsh0 + 423);
    const auto *lsh0_425 = buffer.data(lsh0 + 425);
    const auto *lsh0_426 = buffer.data(lsh0 + 426);
    const auto *lsh0_429 = buffer.data(lsh0 + 429);
    const auto *lsh0_440 = buffer.data(lsh0 + 440);

    const auto *lsg_243 = buffer.data(lsg + 243);
    const auto *lsg_250 = buffer.data(lsg + 250);
    const auto *lsg_254 = buffer.data(lsg + 254);
    const auto *lsg_255 = buffer.data(lsg + 255);
    const auto *lsg_258 = buffer.data(lsg + 258);
    const auto *lsg_260 = buffer.data(lsg + 260);
    const auto *lsg_265 = buffer.data(lsg + 265);
    const auto *lsg_267 = buffer.data(lsg + 267);
    const auto *lsg_268 = buffer.data(lsg + 268);
    const auto *lsg_269 = buffer.data(lsg + 269);
    const auto *lsg_270 = buffer.data(lsg + 270);
    const auto *lsg_272 = buffer.data(lsg + 272);
    const auto *lsg_273 = buffer.data(lsg + 273);
    const auto *lsg_275 = buffer.data(lsg + 275);
    const auto *lsg_280 = buffer.data(lsg + 280);
    const auto *lsg_282 = buffer.data(lsg + 282);
    const auto *lsg_283 = buffer.data(lsg + 283);
    const auto *lsg_284 = buffer.data(lsg + 284);
    const auto *lsg_285 = buffer.data(lsg + 285);
    const auto *lsg_287 = buffer.data(lsg + 287);
    const auto *lsg_288 = buffer.data(lsg + 288);
    const auto *lsg_290 = buffer.data(lsg + 290);
    const auto *lsg_295 = buffer.data(lsg + 295);
    const auto *lsg_297 = buffer.data(lsg + 297);
    const auto *lsg_298 = buffer.data(lsg + 298);
    const auto *lsg_299 = buffer.data(lsg + 299);
    const auto *lsg_300 = buffer.data(lsg + 300);
    const auto *lsg_301 = buffer.data(lsg + 301);
    const auto *lsg_302 = buffer.data(lsg + 302);
    const auto *lsg_303 = buffer.data(lsg + 303);
    const auto *lsg_305 = buffer.data(lsg + 305);
    const auto *lsg_310 = buffer.data(lsg + 310);
    const auto *lsg_312 = buffer.data(lsg + 312);
    const auto *lsg_313 = buffer.data(lsg + 313);
    const auto *lsg_314 = buffer.data(lsg + 314);
    const auto *lsg_315 = buffer.data(lsg + 315);
    const auto *lsg_320 = buffer.data(lsg + 320);
    const auto *lsg_325 = buffer.data(lsg + 325);
    const auto *lsg_350 = buffer.data(lsg + 350);
    const auto *lsg_351 = buffer.data(lsg + 351);
    const auto *lsg_354 = buffer.data(lsg + 354);
    const auto *lsg_355 = buffer.data(lsg + 355);
    const auto *lsg_356 = buffer.data(lsg + 356);
    const auto *lsg_357 = buffer.data(lsg + 357);
    const auto *lsg_358 = buffer.data(lsg + 358);
    const auto *lsg_359 = buffer.data(lsg + 359);
    const auto *lsg_360 = buffer.data(lsg + 360);
    const auto *lsg_363 = buffer.data(lsg + 363);
    const auto *lsg_365 = buffer.data(lsg + 365);
    const auto *lsg_366 = buffer.data(lsg + 366);
    const auto *lsg_369 = buffer.data(lsg + 369);
    const auto *lsg_370 = buffer.data(lsg + 370);
    const auto *lsg_371 = buffer.data(lsg + 371);
    const auto *lsg_372 = buffer.data(lsg + 372);
    const auto *lsg_373 = buffer.data(lsg + 373);
    const auto *lsg_374 = buffer.data(lsg + 374);
    const auto *lsg_375 = buffer.data(lsg + 375);
    const auto *lsg_378 = buffer.data(lsg + 378);
    const auto *lsg_380 = buffer.data(lsg + 380);
    const auto *lsg_381 = buffer.data(lsg + 381);
    const auto *lsg_384 = buffer.data(lsg + 384);
    const auto *lsg_385 = buffer.data(lsg + 385);
    const auto *lsg_386 = buffer.data(lsg + 386);
    const auto *lsg_387 = buffer.data(lsg + 387);
    const auto *lsg_388 = buffer.data(lsg + 388);
    const auto *lsg_389 = buffer.data(lsg + 389);
    const auto *lsg_400 = buffer.data(lsg + 400);
    const auto *lsg_401 = buffer.data(lsg + 401);
    const auto *lsg_402 = buffer.data(lsg + 402);
    const auto *lsg_403 = buffer.data(lsg + 403);
    const auto *lsg_404 = buffer.data(lsg + 404);
    const auto *lsg_405 = buffer.data(lsg + 405);
    const auto *lsg_410 = buffer.data(lsg + 410);
    const auto *lsg_414 = buffer.data(lsg + 414);
    const auto *lsg_415 = buffer.data(lsg + 415);
    const auto *lsg_416 = buffer.data(lsg + 416);
    const auto *lsg_417 = buffer.data(lsg + 417);
    const auto *lsg_419 = buffer.data(lsg + 419);
    const auto *lsg_420 = buffer.data(lsg + 420);
    const auto *lsg_423 = buffer.data(lsg + 423);
    const auto *lsg_426 = buffer.data(lsg + 426);
    const auto *lsg_430 = buffer.data(lsg + 430);
    const auto *lsg_432 = buffer.data(lsg + 432);
    const auto *lsg_433 = buffer.data(lsg + 433);
    const auto *lsg_434 = buffer.data(lsg + 434);

    const auto *lsh1_420 = buffer.data(lsh1 + 420);
    const auto *lsh1_423 = buffer.data(lsh1 + 423);
    const auto *lsh1_425 = buffer.data(lsh1 + 425);
    const auto *lsh1_426 = buffer.data(lsh1 + 426);
    const auto *lsh1_429 = buffer.data(lsh1 + 429);
    const auto *lsh1_440 = buffer.data(lsh1 + 440);

    const auto *msf0_235 = buffer.data(msf0 + 235);
    const auto *msf0_236 = buffer.data(msf0 + 236);
    const auto *msf0_238 = buffer.data(msf0 + 238);
    const auto *msf0_239 = buffer.data(msf0 + 239);
    const auto *msf0_240 = buffer.data(msf0 + 240);
    const auto *msf0_243 = buffer.data(msf0 + 243);
    const auto *msf0_245 = buffer.data(msf0 + 245);
    const auto *msf0_246 = buffer.data(msf0 + 246);
    const auto *msf0_248 = buffer.data(msf0 + 248);
    const auto *msf0_249 = buffer.data(msf0 + 249);
    const auto *msf0_250 = buffer.data(msf0 + 250);
    const auto *msf0_253 = buffer.data(msf0 + 253);
    const auto *msf0_255 = buffer.data(msf0 + 255);
    const auto *msf0_256 = buffer.data(msf0 + 256);
    const auto *msf0_258 = buffer.data(msf0 + 258);
    const auto *msf0_259 = buffer.data(msf0 + 259);
    const auto *msf0_266 = buffer.data(msf0 + 266);
    const auto *msf0_268 = buffer.data(msf0 + 268);
    const auto *msf0_269 = buffer.data(msf0 + 269);
    const auto *msf0_270 = buffer.data(msf0 + 270);
    const auto *msf0_271 = buffer.data(msf0 + 271);
    const auto *msf0_272 = buffer.data(msf0 + 272);
    const auto *msf0_275 = buffer.data(msf0 + 275);
    const auto *msf0_276 = buffer.data(msf0 + 276);
    const auto *msf0_277 = buffer.data(msf0 + 277);
    const auto *msf0_278 = buffer.data(msf0 + 278);
    const auto *msf0_279 = buffer.data(msf0 + 279);
    const auto *msf0_280 = buffer.data(msf0 + 280);
    const auto *msf0_282 = buffer.data(msf0 + 282);
    const auto *msf0_283 = buffer.data(msf0 + 283);
    const auto *msf0_286 = buffer.data(msf0 + 286);

    const auto *msf1_235 = buffer.data(msf1 + 235);
    const auto *msf1_236 = buffer.data(msf1 + 236);
    const auto *msf1_238 = buffer.data(msf1 + 238);
    const auto *msf1_239 = buffer.data(msf1 + 239);
    const auto *msf1_240 = buffer.data(msf1 + 240);
    const auto *msf1_243 = buffer.data(msf1 + 243);
    const auto *msf1_245 = buffer.data(msf1 + 245);
    const auto *msf1_246 = buffer.data(msf1 + 246);
    const auto *msf1_248 = buffer.data(msf1 + 248);
    const auto *msf1_249 = buffer.data(msf1 + 249);
    const auto *msf1_250 = buffer.data(msf1 + 250);
    const auto *msf1_253 = buffer.data(msf1 + 253);
    const auto *msf1_255 = buffer.data(msf1 + 255);
    const auto *msf1_256 = buffer.data(msf1 + 256);
    const auto *msf1_258 = buffer.data(msf1 + 258);
    const auto *msf1_259 = buffer.data(msf1 + 259);
    const auto *msf1_266 = buffer.data(msf1 + 266);
    const auto *msf1_268 = buffer.data(msf1 + 268);
    const auto *msf1_269 = buffer.data(msf1 + 269);
    const auto *msf1_270 = buffer.data(msf1 + 270);
    const auto *msf1_271 = buffer.data(msf1 + 271);
    const auto *msf1_272 = buffer.data(msf1 + 272);
    const auto *msf1_275 = buffer.data(msf1 + 275);
    const auto *msf1_276 = buffer.data(msf1 + 276);
    const auto *msf1_277 = buffer.data(msf1 + 277);
    const auto *msf1_278 = buffer.data(msf1 + 278);
    const auto *msf1_279 = buffer.data(msf1 + 279);
    const auto *msf1_280 = buffer.data(msf1 + 280);
    const auto *msf1_282 = buffer.data(msf1 + 282);
    const auto *msf1_283 = buffer.data(msf1 + 283);
    const auto *msf1_286 = buffer.data(msf1 + 286);

    const auto *msg_348 = buffer.data(msg + 348);
    const auto *msg_350 = buffer.data(msg + 350);
    const auto *msg_351 = buffer.data(msg + 351);
    const auto *msg_354 = buffer.data(msg + 354);
    const auto *msg_355 = buffer.data(msg + 355);
    const auto *msg_356 = buffer.data(msg + 356);
    const auto *msg_357 = buffer.data(msg + 357);
    const auto *msg_358 = buffer.data(msg + 358);
    const auto *msg_359 = buffer.data(msg + 359);
    const auto *msg_360 = buffer.data(msg + 360);
    const auto *msg_362 = buffer.data(msg + 362);
    const auto *msg_363 = buffer.data(msg + 363);
    const auto *msg_365 = buffer.data(msg + 365);
    const auto *msg_366 = buffer.data(msg + 366);
    const auto *msg_369 = buffer.data(msg + 369);
    const auto *msg_370 = buffer.data(msg + 370);
    const auto *msg_371 = buffer.data(msg + 371);
    const auto *msg_372 = buffer.data(msg + 372);
    const auto *msg_373 = buffer.data(msg + 373);
    const auto *msg_374 = buffer.data(msg + 374);
    const auto *msg_375 = buffer.data(msg + 375);
    const auto *msg_377 = buffer.data(msg + 377);
    const auto *msg_378 = buffer.data(msg + 378);
    const auto *msg_380 = buffer.data(msg + 380);
    const auto *msg_381 = buffer.data(msg + 381);
    const auto *msg_384 = buffer.data(msg + 384);
    const auto *msg_385 = buffer.data(msg + 385);
    const auto *msg_386 = buffer.data(msg + 386);
    const auto *msg_387 = buffer.data(msg + 387);
    const auto *msg_388 = buffer.data(msg + 388);
    const auto *msg_389 = buffer.data(msg + 389);
    const auto *msg_390 = buffer.data(msg + 390);
    const auto *msg_392 = buffer.data(msg + 392);
    const auto *msg_393 = buffer.data(msg + 393);
    const auto *msg_395 = buffer.data(msg + 395);
    const auto *msg_400 = buffer.data(msg + 400);
    const auto *msg_401 = buffer.data(msg + 401);
    const auto *msg_402 = buffer.data(msg + 402);
    const auto *msg_403 = buffer.data(msg + 403);
    const auto *msg_404 = buffer.data(msg + 404);
    const auto *msg_405 = buffer.data(msg + 405);
    const auto *msg_406 = buffer.data(msg + 406);
    const auto *msg_407 = buffer.data(msg + 407);
    const auto *msg_408 = buffer.data(msg + 408);
    const auto *msg_409 = buffer.data(msg + 409);
    const auto *msg_410 = buffer.data(msg + 410);
    const auto *msg_414 = buffer.data(msg + 414);
    const auto *msg_415 = buffer.data(msg + 415);
    const auto *msg_416 = buffer.data(msg + 416);
    const auto *msg_417 = buffer.data(msg + 417);
    const auto *msg_418 = buffer.data(msg + 418);
    const auto *msg_419 = buffer.data(msg + 419);
    const auto *msg_420 = buffer.data(msg + 420);
    const auto *msg_421 = buffer.data(msg + 421);
    const auto *msg_422 = buffer.data(msg + 422);
    const auto *msg_423 = buffer.data(msg + 423);
    const auto *msg_425 = buffer.data(msg + 425);
    const auto *msg_426 = buffer.data(msg + 426);
    const auto *msg_430 = buffer.data(msg + 430);
    const auto *msg_432 = buffer.data(msg + 432);
    const auto *msg_433 = buffer.data(msg + 433);
    const auto *msg_434 = buffer.data(msg + 434);

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, lsg_243, lsg_350, lsg_351, msf0_235, \
                         msf0_236, msf1_235, msf1_236, msg_348, msg_350, \
                         msg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_11 * lsg_350[k]
                   + f_6 * msf0_235[k]
                   - f_7 * msf1_235[k]
                   + f_3 * pc_x[k] * msg_350[k];

        t_489[k] = f_11 * lsg_351[k]
                   + f_4 * msf0_236[k]
                   - f_5 * msf1_236[k]
                   + f_3 * pc_x[k] * msg_351[k];

        t_490[k] = f_10 * lsg_243[k]
                   + f_3 * pc_z[k] * msg_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, lsg_260, lsg_354, lsg_355, \
                         lsg_356, msf0_239, msf1_239, msg_350, msg_354, msg_355, \
                         msg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_18 * lsg_260[k]
                   + f_3 * pc_y[k] * msg_350[k];

        t_492[k] = f_11 * lsg_354[k]
                   + f_4 * msf0_239[k]
                   - f_5 * msf1_239[k]
                   + f_3 * pc_x[k] * msg_354[k];

        t_493[k] = f_11 * lsg_355[k]
                   + f_3 * pc_x[k] * msg_355[k];

        t_494[k] = f_11 * lsg_356[k]
                   + f_3 * pc_x[k] * msg_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, lsg_265, lsg_357, lsg_358, \
                         lsg_359, msf0_236, msf1_236, msg_355, msg_357, msg_358, \
                         msg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_11 * lsg_357[k]
                   + f_3 * pc_x[k] * msg_357[k];

        t_496[k] = f_11 * lsg_358[k]
                   + f_3 * pc_x[k] * msg_358[k];

        t_497[k] = f_11 * lsg_359[k]
                   + f_3 * pc_x[k] * msg_359[k];

        t_498[k] = f_18 * lsg_265[k]
                   + f_1 * msf0_236[k]
                   - f_2 * msf1_236[k]
                   + f_3 * pc_y[k] * msg_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, lsg_250, lsg_267, lsg_268, msf0_238, \
                         msf0_239, msf1_238, msf1_239, msg_355, msg_357, \
                         msg_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * lsg_250[k]
                   + f_3 * pc_z[k] * msg_355[k];

        t_500[k] = f_18 * lsg_267[k]
                   + f_6 * msf0_238[k]
                   - f_7 * msf1_238[k]
                   + f_3 * pc_y[k] * msg_357[k];

        t_501[k] = f_18 * lsg_268[k]
                   + f_4 * msf0_239[k]
                   - f_5 * msf1_239[k]
                   + f_3 * pc_y[k] * msg_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, lsg_254, lsg_269, lsg_360, \
                         msf0_239, msf0_240, msf1_239, msf1_240, msg_359, \
                         msg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_18 * lsg_269[k]
                   + f_3 * pc_y[k] * msg_359[k];

        t_503[k] = f_10 * lsg_254[k]
                   + f_1 * msf0_239[k]
                   - f_2 * msf1_239[k]
                   + f_3 * pc_z[k] * msg_359[k];

        t_504[k] = f_11 * lsg_360[k]
                   + f_1 * msf0_240[k]
                   - f_2 * msf1_240[k]
                   + f_3 * pc_x[k] * msg_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, lsg_255, lsg_270, \
                         lsg_272, lsg_363, msf0_243, msf1_243, msg_360, msg_362, \
                         msg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * lsg_270[k]
                   + f_3 * pc_y[k] * msg_360[k];

        t_506[k] = f_11 * lsg_255[k]
                   + f_3 * pc_z[k] * msg_360[k];

        t_507[k] = f_11 * lsg_363[k]
                   + f_6 * msf0_243[k]
                   - f_7 * msf1_243[k]
                   + f_3 * pc_x[k] * msg_363[k];

        t_508[k] = f_11 * lsg_272[k]
                   + f_3 * pc_y[k] * msg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, lsg_258, lsg_365, lsg_366, msf0_245, \
                         msf0_246, msf1_245, msf1_246, msg_363, msg_365, \
                         msg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_11 * lsg_365[k]
                   + f_6 * msf0_245[k]
                   - f_7 * msf1_245[k]
                   + f_3 * pc_x[k] * msg_365[k];

        t_510[k] = f_11 * lsg_366[k]
                   + f_4 * msf0_246[k]
                   - f_5 * msf1_246[k]
                   + f_3 * pc_x[k] * msg_366[k];

        t_511[k] = f_11 * lsg_258[k]
                   + f_3 * pc_z[k] * msg_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, lsg_275, lsg_369, lsg_370, \
                         lsg_371, msf0_249, msf1_249, msg_365, msg_369, msg_370, \
                         msg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * lsg_275[k]
                   + f_3 * pc_y[k] * msg_365[k];

        t_513[k] = f_11 * lsg_369[k]
                   + f_4 * msf0_249[k]
                   - f_5 * msf1_249[k]
                   + f_3 * pc_x[k] * msg_369[k];

        t_514[k] = f_11 * lsg_370[k]
                   + f_3 * pc_x[k] * msg_370[k];

        t_515[k] = f_11 * lsg_371[k]
                   + f_3 * pc_x[k] * msg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, lsg_280, lsg_372, lsg_373, \
                         lsg_374, msf0_246, msf1_246, msg_370, msg_372, msg_373, \
                         msg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_11 * lsg_372[k]
                   + f_3 * pc_x[k] * msg_372[k];

        t_517[k] = f_11 * lsg_373[k]
                   + f_3 * pc_x[k] * msg_373[k];

        t_518[k] = f_11 * lsg_374[k]
                   + f_3 * pc_x[k] * msg_374[k];

        t_519[k] = f_11 * lsg_280[k]
                   + f_1 * msf0_246[k]
                   - f_2 * msf1_246[k]
                   + f_3 * pc_y[k] * msg_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, lsg_265, lsg_282, lsg_283, msf0_248, \
                         msf0_249, msf1_248, msf1_249, msg_370, msg_372, \
                         msg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * lsg_265[k]
                   + f_3 * pc_z[k] * msg_370[k];

        t_521[k] = f_11 * lsg_282[k]
                   + f_6 * msf0_248[k]
                   - f_7 * msf1_248[k]
                   + f_3 * pc_y[k] * msg_372[k];

        t_522[k] = f_11 * lsg_283[k]
                   + f_4 * msf0_249[k]
                   - f_5 * msf1_249[k]
                   + f_3 * pc_y[k] * msg_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, lsg_269, lsg_284, lsg_375, \
                         msf0_249, msf0_250, msf1_249, msf1_250, msg_374, \
                         msg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * lsg_284[k]
                   + f_3 * pc_y[k] * msg_374[k];

        t_524[k] = f_11 * lsg_269[k]
                   + f_1 * msf0_249[k]
                   - f_2 * msf1_249[k]
                   + f_3 * pc_z[k] * msg_374[k];

        t_525[k] = f_11 * lsg_375[k]
                   + f_1 * msf0_250[k]
                   - f_2 * msf1_250[k]
                   + f_3 * pc_x[k] * msg_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, lsg_270, lsg_285, \
                         lsg_287, lsg_378, msf0_253, msf1_253, msg_375, msg_377, \
                         msg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * lsg_285[k]
                   + f_3 * pc_y[k] * msg_375[k];

        t_527[k] = f_18 * lsg_270[k]
                   + f_3 * pc_z[k] * msg_375[k];

        t_528[k] = f_11 * lsg_378[k]
                   + f_6 * msf0_253[k]
                   - f_7 * msf1_253[k]
                   + f_3 * pc_x[k] * msg_378[k];

        t_529[k] = f_10 * lsg_287[k]
                   + f_3 * pc_y[k] * msg_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, lsg_273, lsg_380, lsg_381, msf0_255, \
                         msf0_256, msf1_255, msf1_256, msg_378, msg_380, \
                         msg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_11 * lsg_380[k]
                   + f_6 * msf0_255[k]
                   - f_7 * msf1_255[k]
                   + f_3 * pc_x[k] * msg_380[k];

        t_531[k] = f_11 * lsg_381[k]
                   + f_4 * msf0_256[k]
                   - f_5 * msf1_256[k]
                   + f_3 * pc_x[k] * msg_381[k];

        t_532[k] = f_18 * lsg_273[k]
                   + f_3 * pc_z[k] * msg_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, lsg_290, lsg_384, lsg_385, \
                         lsg_386, msf0_259, msf1_259, msg_380, msg_384, msg_385, \
                         msg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * lsg_290[k]
                   + f_3 * pc_y[k] * msg_380[k];

        t_534[k] = f_11 * lsg_384[k]
                   + f_4 * msf0_259[k]
                   - f_5 * msf1_259[k]
                   + f_3 * pc_x[k] * msg_384[k];

        t_535[k] = f_11 * lsg_385[k]
                   + f_3 * pc_x[k] * msg_385[k];

        t_536[k] = f_11 * lsg_386[k]
                   + f_3 * pc_x[k] * msg_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, lsg_295, lsg_387, lsg_388, \
                         lsg_389, msf0_256, msf1_256, msg_385, msg_387, msg_388, \
                         msg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_11 * lsg_387[k]
                   + f_3 * pc_x[k] * msg_387[k];

        t_538[k] = f_11 * lsg_388[k]
                   + f_3 * pc_x[k] * msg_388[k];

        t_539[k] = f_11 * lsg_389[k]
                   + f_3 * pc_x[k] * msg_389[k];

        t_540[k] = f_10 * lsg_295[k]
                   + f_1 * msf0_256[k]
                   - f_2 * msf1_256[k]
                   + f_3 * pc_y[k] * msg_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, lsg_280, lsg_297, lsg_298, msf0_258, \
                         msf0_259, msf1_258, msf1_259, msg_385, msg_387, \
                         msg_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_18 * lsg_280[k]
                   + f_3 * pc_z[k] * msg_385[k];

        t_542[k] = f_10 * lsg_297[k]
                   + f_6 * msf0_258[k]
                   - f_7 * msf1_258[k]
                   + f_3 * pc_y[k] * msg_387[k];

        t_543[k] = f_10 * lsg_298[k]
                   + f_4 * msf0_259[k]
                   - f_5 * msf1_259[k]
                   + f_3 * pc_y[k] * msg_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_y, pc_z, lsh0_420, lsg_284, \
                         lsg_299, lsg_300, lsh1_420, msf0_259, msf1_259, msg_389, \
                         msg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * lsg_299[k]
                   + f_3 * pc_y[k] * msg_389[k];

        t_545[k] = f_18 * lsg_284[k]
                   + f_1 * msf0_259[k]
                   - f_2 * msf1_259[k]
                   + f_3 * pc_z[k] * msg_389[k];

        t_546[k] = pa_y[k] * lsh0_420[k]
                   - f_8 * pc_y[k] * lsh1_420[k];

        t_547[k] = f_9 * lsg_300[k]
                   + f_3 * pc_y[k] * msg_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pa_y, pc_y, pc_z, lsh0_423, lsh0_425, \
                         lsg_285, lsg_301, lsg_302, lsh1_423, lsh1_425, msg_390, \
                         msg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_17 * lsg_285[k]
                   + f_3 * pc_z[k] * msg_390[k];

        t_549[k] = pa_y[k] * lsh0_423[k]
                   + f_10 * lsg_301[k]
                   - f_8 * pc_y[k] * lsh1_423[k];

        t_550[k] = f_9 * lsg_302[k]
                   + f_3 * pc_y[k] * msg_392[k];

        t_551[k] = pa_y[k] * lsh0_425[k]
                   - f_8 * pc_y[k] * lsh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pa_y, pc_y, pc_z, lsh0_426, lsh0_429, \
                         lsg_288, lsg_303, lsg_305, lsh1_426, lsh1_429, msg_393, \
                         msg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_y[k] * lsh0_426[k]
                   + f_11 * lsg_303[k]
                   - f_8 * pc_y[k] * lsh1_426[k];

        t_553[k] = f_17 * lsg_288[k]
                   + f_3 * pc_z[k] * msg_393[k];

        t_554[k] = f_9 * lsg_305[k]
                   + f_3 * pc_y[k] * msg_395[k];

        t_555[k] = pa_y[k] * lsh0_429[k]
                   - f_8 * pc_y[k] * lsh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, lsg_400, lsg_401, lsg_402, \
                         lsg_403, lsg_404, msg_400, msg_401, msg_402, msg_403, \
                         msg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * lsg_400[k]
                   + f_3 * pc_x[k] * msg_400[k];

        t_557[k] = f_11 * lsg_401[k]
                   + f_3 * pc_x[k] * msg_401[k];

        t_558[k] = f_11 * lsg_402[k]
                   + f_3 * pc_x[k] * msg_402[k];

        t_559[k] = f_11 * lsg_403[k]
                   + f_3 * pc_x[k] * msg_403[k];

        t_560[k] = f_11 * lsg_404[k]
                   + f_3 * pc_x[k] * msg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, lsg_295, lsg_310, lsg_312, msf0_266, \
                         msf0_268, msf1_266, msf1_268, msg_400, \
                         msg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * lsg_310[k]
                   + f_1 * msf0_266[k]
                   - f_2 * msf1_266[k]
                   + f_3 * pc_y[k] * msg_400[k];

        t_562[k] = f_17 * lsg_295[k]
                   + f_3 * pc_z[k] * msg_400[k];

        t_563[k] = f_9 * lsg_312[k]
                   + f_6 * msf0_268[k]
                   - f_7 * msf1_268[k]
                   + f_3 * pc_y[k] * msg_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_y, pc_y, lsh0_440, lsg_313, lsg_314, \
                         lsh1_440, msf0_269, msf1_269, msg_403, \
                         msg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * lsg_313[k]
                   + f_4 * msf0_269[k]
                   - f_5 * msf1_269[k]
                   + f_3 * pc_y[k] * msg_403[k];

        t_565[k] = f_9 * lsg_314[k]
                   + f_3 * pc_y[k] * msg_404[k];

        t_566[k] = pa_y[k] * lsh0_440[k]
                   - f_8 * pc_y[k] * lsh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pc_x, pc_y, pc_z, lsg_300, \
                         lsg_405, msf0_270, msf1_270, msg_405, msg_406, \
                         msg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_11 * lsg_405[k]
                   + f_1 * msf0_270[k]
                   - f_2 * msf1_270[k]
                   + f_3 * pc_x[k] * msg_405[k];

        t_568[k] = f_3 * pc_y[k] * msg_405[k];

        t_569[k] = f_16 * lsg_300[k]
                   + f_3 * pc_z[k] * msg_405[k];

        t_570[k] = f_4 * msf0_270[k]
                   - f_5 * msf1_270[k]
                   + f_3 * pc_y[k] * msg_406[k];

        t_571[k] = f_3 * pc_y[k] * msg_407[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, lsg_410, msf0_271, msf0_272, \
                         msf0_275, msf1_271, msf1_272, msf1_275, msg_408, msg_409, \
                         msg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_11 * lsg_410[k]
                   + f_6 * msf0_275[k]
                   - f_7 * msf1_275[k]
                   + f_3 * pc_x[k] * msg_410[k];

        t_573[k] = f_6 * msf0_271[k]
                   - f_7 * msf1_271[k]
                   + f_3 * pc_y[k] * msg_408[k];

        t_574[k] = f_4 * msf0_272[k]
                   - f_5 * msf1_272[k]
                   + f_3 * pc_y[k] * msg_409[k];

        t_575[k] = f_3 * pc_y[k] * msg_410[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pc_x, lsg_414, lsg_415, lsg_416, lsg_417, \
                         msf0_279, msf1_279, msg_414, msg_415, msg_416, \
                         msg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_11 * lsg_414[k]
                   + f_4 * msf0_279[k]
                   - f_5 * msf1_279[k]
                   + f_3 * pc_x[k] * msg_414[k];

        t_577[k] = f_11 * lsg_415[k]
                   + f_3 * pc_x[k] * msg_415[k];

        t_578[k] = f_11 * lsg_416[k]
                   + f_3 * pc_x[k] * msg_416[k];

        t_579[k] = f_11 * lsg_417[k]
                   + f_3 * pc_x[k] * msg_417[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, lsg_419, msf0_276, msf0_277, \
                         msf1_276, msf1_277, msg_414, msg_415, msg_416, \
                         msg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_3 * pc_y[k] * msg_414[k];

        t_581[k] = f_11 * lsg_419[k]
                   + f_3 * pc_x[k] * msg_419[k];

        t_582[k] = f_1 * msf0_276[k]
                   - f_2 * msf1_276[k]
                   + f_3 * pc_y[k] * msg_415[k];

        t_583[k] = f_13 * msf0_277[k]
                   - f_14 * msf1_277[k]
                   + f_3 * pc_y[k] * msg_416[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, lsg_314, msf0_278, msf0_279, \
                         msf1_278, msf1_279, msg_417, msg_418, \
                         msg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * msf0_278[k]
                   - f_7 * msf1_278[k]
                   + f_3 * pc_y[k] * msg_417[k];

        t_585[k] = f_4 * msf0_279[k]
                   - f_5 * msf1_279[k]
                   + f_3 * pc_y[k] * msg_418[k];

        t_586[k] = f_3 * pc_y[k] * msg_419[k];

        t_587[k] = f_16 * lsg_314[k]
                   + f_1 * msf0_279[k]
                   - f_2 * msf1_279[k]
                   + f_3 * pc_z[k] * msg_419[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, lsg_315, lsg_420, \
                         lsg_423, msf0_280, msf0_283, msf1_280, msf1_283, msg_420, \
                         msg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_10 * lsg_420[k]
                   + f_1 * msf0_280[k]
                   - f_2 * msf1_280[k]
                   + f_3 * pc_x[k] * msg_420[k];

        t_589[k] = f_15 * lsg_315[k]
                   + f_3 * pc_y[k] * msg_420[k];

        t_590[k] = f_3 * pc_z[k] * msg_420[k];

        t_591[k] = f_10 * lsg_423[k]
                   + f_6 * msf0_283[k]
                   - f_7 * msf1_283[k]
                   + f_3 * pc_x[k] * msg_423[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, lsg_426, msf0_280, msf0_286, \
                         msf1_280, msf1_286, msg_421, msg_422, msg_423, \
                         msg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * msg_421[k];

        t_593[k] = f_4 * msf0_280[k]
                   - f_5 * msf1_280[k]
                   + f_3 * pc_z[k] * msg_422[k];

        t_594[k] = f_10 * lsg_426[k]
                   + f_4 * msf0_286[k]
                   - f_5 * msf1_286[k]
                   + f_3 * pc_x[k] * msg_426[k];

        t_595[k] = f_3 * pc_z[k] * msg_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, lsg_320, lsg_430, \
                         msf0_282, msf1_282, msg_425, msg_426, \
                         msg_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_15 * lsg_320[k]
                   + f_3 * pc_y[k] * msg_425[k];

        t_597[k] = f_6 * msf0_282[k]
                   - f_7 * msf1_282[k]
                   + f_3 * pc_z[k] * msg_425[k];

        t_598[k] = f_10 * lsg_430[k]
                   + f_3 * pc_x[k] * msg_430[k];

        t_599[k] = f_3 * pc_z[k] * msg_426[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, lsg_325, lsg_432, lsg_433, \
                         lsg_434, msf0_286, msf1_286, msg_430, msg_432, msg_433, \
                         msg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_10 * lsg_432[k]
                   + f_3 * pc_x[k] * msg_432[k];

        t_601[k] = f_10 * lsg_433[k]
                   + f_3 * pc_x[k] * msg_433[k];

        t_602[k] = f_10 * lsg_434[k]
                   + f_3 * pc_x[k] * msg_434[k];

        t_603[k] = f_15 * lsg_325[k]
                   + f_1 * msf0_286[k]
                   - f_2 * msf1_286[k]
                   + f_3 * pc_y[k] * msg_430[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_441 = buffer.data(lsh0 + 441);
    const auto *lsh0_444 = buffer.data(lsh0 + 444);
    const auto *lsh0_447 = buffer.data(lsh0 + 447);
    const auto *lsh0_456 = buffer.data(lsh0 + 456);

    const auto *lsg_315 = buffer.data(lsg + 315);
    const auto *lsg_318 = buffer.data(lsg + 318);
    const auto *lsg_325 = buffer.data(lsg + 325);
    const auto *lsg_329 = buffer.data(lsg + 329);
    const auto *lsg_330 = buffer.data(lsg + 330);
    const auto *lsg_332 = buffer.data(lsg + 332);
    const auto *lsg_333 = buffer.data(lsg + 333);
    const auto *lsg_335 = buffer.data(lsg + 335);
    const auto *lsg_340 = buffer.data(lsg + 340);
    const auto *lsg_342 = buffer.data(lsg + 342);
    const auto *lsg_343 = buffer.data(lsg + 343);
    const auto *lsg_344 = buffer.data(lsg + 344);
    const auto *lsg_345 = buffer.data(lsg + 345);
    const auto *lsg_347 = buffer.data(lsg + 347);
    const auto *lsg_348 = buffer.data(lsg + 348);
    const auto *lsg_350 = buffer.data(lsg + 350);
    const auto *lsg_355 = buffer.data(lsg + 355);
    const auto *lsg_357 = buffer.data(lsg + 357);
    const auto *lsg_358 = buffer.data(lsg + 358);
    const auto *lsg_359 = buffer.data(lsg + 359);
    const auto *lsg_360 = buffer.data(lsg + 360);
    const auto *lsg_362 = buffer.data(lsg + 362);
    const auto *lsg_363 = buffer.data(lsg + 363);
    const auto *lsg_365 = buffer.data(lsg + 365);
    const auto *lsg_370 = buffer.data(lsg + 370);
    const auto *lsg_372 = buffer.data(lsg + 372);
    const auto *lsg_373 = buffer.data(lsg + 373);
    const auto *lsg_374 = buffer.data(lsg + 374);
    const auto *lsg_375 = buffer.data(lsg + 375);
    const auto *lsg_377 = buffer.data(lsg + 377);
    const auto *lsg_378 = buffer.data(lsg + 378);
    const auto *lsg_380 = buffer.data(lsg + 380);
    const auto *lsg_385 = buffer.data(lsg + 385);
    const auto *lsg_387 = buffer.data(lsg + 387);
    const auto *lsg_388 = buffer.data(lsg + 388);
    const auto *lsg_389 = buffer.data(lsg + 389);
    const auto *lsg_390 = buffer.data(lsg + 390);
    const auto *lsg_392 = buffer.data(lsg + 392);
    const auto *lsg_395 = buffer.data(lsg + 395);
    const auto *lsg_400 = buffer.data(lsg + 400);
    const auto *lsg_402 = buffer.data(lsg + 402);
    const auto *lsg_403 = buffer.data(lsg + 403);
    const auto *lsg_440 = buffer.data(lsg + 440);
    const auto *lsg_444 = buffer.data(lsg + 444);
    const auto *lsg_445 = buffer.data(lsg + 445);
    const auto *lsg_446 = buffer.data(lsg + 446);
    const auto *lsg_447 = buffer.data(lsg + 447);
    const auto *lsg_448 = buffer.data(lsg + 448);
    const auto *lsg_449 = buffer.data(lsg + 449);
    const auto *lsg_450 = buffer.data(lsg + 450);
    const auto *lsg_453 = buffer.data(lsg + 453);
    const auto *lsg_455 = buffer.data(lsg + 455);
    const auto *lsg_456 = buffer.data(lsg + 456);
    const auto *lsg_459 = buffer.data(lsg + 459);
    const auto *lsg_460 = buffer.data(lsg + 460);
    const auto *lsg_461 = buffer.data(lsg + 461);
    const auto *lsg_462 = buffer.data(lsg + 462);
    const auto *lsg_463 = buffer.data(lsg + 463);
    const auto *lsg_464 = buffer.data(lsg + 464);
    const auto *lsg_465 = buffer.data(lsg + 465);
    const auto *lsg_468 = buffer.data(lsg + 468);
    const auto *lsg_470 = buffer.data(lsg + 470);
    const auto *lsg_471 = buffer.data(lsg + 471);
    const auto *lsg_474 = buffer.data(lsg + 474);
    const auto *lsg_475 = buffer.data(lsg + 475);
    const auto *lsg_476 = buffer.data(lsg + 476);
    const auto *lsg_477 = buffer.data(lsg + 477);
    const auto *lsg_478 = buffer.data(lsg + 478);
    const auto *lsg_479 = buffer.data(lsg + 479);
    const auto *lsg_480 = buffer.data(lsg + 480);
    const auto *lsg_483 = buffer.data(lsg + 483);
    const auto *lsg_485 = buffer.data(lsg + 485);
    const auto *lsg_486 = buffer.data(lsg + 486);
    const auto *lsg_489 = buffer.data(lsg + 489);
    const auto *lsg_490 = buffer.data(lsg + 490);
    const auto *lsg_491 = buffer.data(lsg + 491);
    const auto *lsg_492 = buffer.data(lsg + 492);
    const auto *lsg_493 = buffer.data(lsg + 493);
    const auto *lsg_494 = buffer.data(lsg + 494);
    const auto *lsg_495 = buffer.data(lsg + 495);
    const auto *lsg_498 = buffer.data(lsg + 498);
    const auto *lsg_500 = buffer.data(lsg + 500);
    const auto *lsg_501 = buffer.data(lsg + 501);
    const auto *lsg_504 = buffer.data(lsg + 504);
    const auto *lsg_505 = buffer.data(lsg + 505);
    const auto *lsg_506 = buffer.data(lsg + 506);
    const auto *lsg_507 = buffer.data(lsg + 507);
    const auto *lsg_508 = buffer.data(lsg + 508);
    const auto *lsg_509 = buffer.data(lsg + 509);

    const auto *lsh1_441 = buffer.data(lsh1 + 441);
    const auto *lsh1_444 = buffer.data(lsh1 + 444);
    const auto *lsh1_447 = buffer.data(lsh1 + 447);
    const auto *lsh1_456 = buffer.data(lsh1 + 456);

    const auto *msf0_286 = buffer.data(msf0 + 286);
    const auto *msf0_287 = buffer.data(msf0 + 287);
    const auto *msf0_289 = buffer.data(msf0 + 289);
    const auto *msf0_295 = buffer.data(msf0 + 295);
    const auto *msf0_298 = buffer.data(msf0 + 298);
    const auto *msf0_299 = buffer.data(msf0 + 299);
    const auto *msf0_300 = buffer.data(msf0 + 300);
    const auto *msf0_303 = buffer.data(msf0 + 303);
    const auto *msf0_305 = buffer.data(msf0 + 305);
    const auto *msf0_306 = buffer.data(msf0 + 306);
    const auto *msf0_308 = buffer.data(msf0 + 308);
    const auto *msf0_309 = buffer.data(msf0 + 309);
    const auto *msf0_310 = buffer.data(msf0 + 310);
    const auto *msf0_313 = buffer.data(msf0 + 313);
    const auto *msf0_315 = buffer.data(msf0 + 315);
    const auto *msf0_316 = buffer.data(msf0 + 316);
    const auto *msf0_318 = buffer.data(msf0 + 318);
    const auto *msf0_319 = buffer.data(msf0 + 319);
    const auto *msf0_320 = buffer.data(msf0 + 320);
    const auto *msf0_323 = buffer.data(msf0 + 323);
    const auto *msf0_325 = buffer.data(msf0 + 325);
    const auto *msf0_326 = buffer.data(msf0 + 326);
    const auto *msf0_328 = buffer.data(msf0 + 328);
    const auto *msf0_329 = buffer.data(msf0 + 329);
    const auto *msf0_330 = buffer.data(msf0 + 330);
    const auto *msf0_333 = buffer.data(msf0 + 333);
    const auto *msf0_335 = buffer.data(msf0 + 335);
    const auto *msf0_336 = buffer.data(msf0 + 336);
    const auto *msf0_338 = buffer.data(msf0 + 338);
    const auto *msf0_339 = buffer.data(msf0 + 339);

    const auto *msf1_286 = buffer.data(msf1 + 286);
    const auto *msf1_287 = buffer.data(msf1 + 287);
    const auto *msf1_289 = buffer.data(msf1 + 289);
    const auto *msf1_295 = buffer.data(msf1 + 295);
    const auto *msf1_298 = buffer.data(msf1 + 298);
    const auto *msf1_299 = buffer.data(msf1 + 299);
    const auto *msf1_300 = buffer.data(msf1 + 300);
    const auto *msf1_303 = buffer.data(msf1 + 303);
    const auto *msf1_305 = buffer.data(msf1 + 305);
    const auto *msf1_306 = buffer.data(msf1 + 306);
    const auto *msf1_308 = buffer.data(msf1 + 308);
    const auto *msf1_309 = buffer.data(msf1 + 309);
    const auto *msf1_310 = buffer.data(msf1 + 310);
    const auto *msf1_313 = buffer.data(msf1 + 313);
    const auto *msf1_315 = buffer.data(msf1 + 315);
    const auto *msf1_316 = buffer.data(msf1 + 316);
    const auto *msf1_318 = buffer.data(msf1 + 318);
    const auto *msf1_319 = buffer.data(msf1 + 319);
    const auto *msf1_320 = buffer.data(msf1 + 320);
    const auto *msf1_323 = buffer.data(msf1 + 323);
    const auto *msf1_325 = buffer.data(msf1 + 325);
    const auto *msf1_326 = buffer.data(msf1 + 326);
    const auto *msf1_328 = buffer.data(msf1 + 328);
    const auto *msf1_329 = buffer.data(msf1 + 329);
    const auto *msf1_330 = buffer.data(msf1 + 330);
    const auto *msf1_333 = buffer.data(msf1 + 333);
    const auto *msf1_335 = buffer.data(msf1 + 335);
    const auto *msf1_336 = buffer.data(msf1 + 336);
    const auto *msf1_338 = buffer.data(msf1 + 338);
    const auto *msf1_339 = buffer.data(msf1 + 339);

    const auto *msg_430 = buffer.data(msg + 430);
    const auto *msg_431 = buffer.data(msg + 431);
    const auto *msg_432 = buffer.data(msg + 432);
    const auto *msg_434 = buffer.data(msg + 434);
    const auto *msg_435 = buffer.data(msg + 435);
    const auto *msg_437 = buffer.data(msg + 437);
    const auto *msg_438 = buffer.data(msg + 438);
    const auto *msg_440 = buffer.data(msg + 440);
    const auto *msg_444 = buffer.data(msg + 444);
    const auto *msg_445 = buffer.data(msg + 445);
    const auto *msg_446 = buffer.data(msg + 446);
    const auto *msg_447 = buffer.data(msg + 447);
    const auto *msg_448 = buffer.data(msg + 448);
    const auto *msg_449 = buffer.data(msg + 449);
    const auto *msg_450 = buffer.data(msg + 450);
    const auto *msg_452 = buffer.data(msg + 452);
    const auto *msg_453 = buffer.data(msg + 453);
    const auto *msg_455 = buffer.data(msg + 455);
    const auto *msg_456 = buffer.data(msg + 456);
    const auto *msg_459 = buffer.data(msg + 459);
    const auto *msg_460 = buffer.data(msg + 460);
    const auto *msg_461 = buffer.data(msg + 461);
    const auto *msg_462 = buffer.data(msg + 462);
    const auto *msg_463 = buffer.data(msg + 463);
    const auto *msg_464 = buffer.data(msg + 464);
    const auto *msg_465 = buffer.data(msg + 465);
    const auto *msg_467 = buffer.data(msg + 467);
    const auto *msg_468 = buffer.data(msg + 468);
    const auto *msg_470 = buffer.data(msg + 470);
    const auto *msg_471 = buffer.data(msg + 471);
    const auto *msg_474 = buffer.data(msg + 474);
    const auto *msg_475 = buffer.data(msg + 475);
    const auto *msg_476 = buffer.data(msg + 476);
    const auto *msg_477 = buffer.data(msg + 477);
    const auto *msg_478 = buffer.data(msg + 478);
    const auto *msg_479 = buffer.data(msg + 479);
    const auto *msg_480 = buffer.data(msg + 480);
    const auto *msg_482 = buffer.data(msg + 482);
    const auto *msg_483 = buffer.data(msg + 483);
    const auto *msg_485 = buffer.data(msg + 485);
    const auto *msg_486 = buffer.data(msg + 486);
    const auto *msg_489 = buffer.data(msg + 489);
    const auto *msg_490 = buffer.data(msg + 490);
    const auto *msg_491 = buffer.data(msg + 491);
    const auto *msg_492 = buffer.data(msg + 492);
    const auto *msg_493 = buffer.data(msg + 493);
    const auto *msg_494 = buffer.data(msg + 494);
    const auto *msg_495 = buffer.data(msg + 495);
    const auto *msg_497 = buffer.data(msg + 497);
    const auto *msg_498 = buffer.data(msg + 498);
    const auto *msg_500 = buffer.data(msg + 500);
    const auto *msg_501 = buffer.data(msg + 501);
    const auto *msg_504 = buffer.data(msg + 504);
    const auto *msg_505 = buffer.data(msg + 505);
    const auto *msg_506 = buffer.data(msg + 506);
    const auto *msg_507 = buffer.data(msg + 507);
    const auto *msg_508 = buffer.data(msg + 508);
    const auto *msg_509 = buffer.data(msg + 509);

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pc_y, pc_z, lsg_329, msf0_286, msf0_287, \
                         msf1_286, msf1_287, msg_430, msg_431, msg_432, \
                         msg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * msg_430[k];

        t_605[k] = f_4 * msf0_286[k]
                   - f_5 * msf1_286[k]
                   + f_3 * pc_z[k] * msg_431[k];

        t_606[k] = f_6 * msf0_287[k]
                   - f_7 * msf1_287[k]
                   + f_3 * pc_z[k] * msg_432[k];

        t_607[k] = f_15 * lsg_329[k]
                   + f_3 * pc_y[k] * msg_434[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pa_z, pc_y, pc_z, lsh0_441, lsg_315, \
                         lsg_330, lsh1_441, msf0_289, msf1_289, msg_434, \
                         msg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_1 * msf0_289[k]
                   - f_2 * msf1_289[k]
                   + f_3 * pc_z[k] * msg_434[k];

        t_609[k] = pa_z[k] * lsh0_441[k]
                   - f_8 * pc_z[k] * lsh1_441[k];

        t_610[k] = f_16 * lsg_330[k]
                   + f_3 * pc_y[k] * msg_435[k];

        t_611[k] = f_9 * lsg_315[k]
                   + f_3 * pc_z[k] * msg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_z, pc_x, pc_y, pc_z, lsh0_444, lsg_332, \
                         lsg_440, lsh1_444, msf0_295, msf1_295, msg_437, \
                         msg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pa_z[k] * lsh0_444[k]
                   - f_8 * pc_z[k] * lsh1_444[k];

        t_613[k] = f_16 * lsg_332[k]
                   + f_3 * pc_y[k] * msg_437[k];

        t_614[k] = f_10 * lsg_440[k]
                   + f_6 * msf0_295[k]
                   - f_7 * msf1_295[k]
                   + f_3 * pc_x[k] * msg_440[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pa_z, pc_y, pc_z, lsh0_447, lsg_318, lsg_335, \
                         lsh1_447, msg_438, msg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = pa_z[k] * lsh0_447[k]
                   - f_8 * pc_z[k] * lsh1_447[k];

        t_616[k] = f_9 * lsg_318[k]
                   + f_3 * pc_z[k] * msg_438[k];

        t_617[k] = f_16 * lsg_335[k]
                   + f_3 * pc_y[k] * msg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pc_x, lsg_444, lsg_445, lsg_446, lsg_447, \
                         msf0_299, msf1_299, msg_444, msg_445, msg_446, \
                         msg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_10 * lsg_444[k]
                   + f_4 * msf0_299[k]
                   - f_5 * msf1_299[k]
                   + f_3 * pc_x[k] * msg_444[k];

        t_619[k] = f_10 * lsg_445[k]
                   + f_3 * pc_x[k] * msg_445[k];

        t_620[k] = f_10 * lsg_446[k]
                   + f_3 * pc_x[k] * msg_446[k];

        t_621[k] = f_10 * lsg_447[k]
                   + f_3 * pc_x[k] * msg_447[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_z, pc_x, pc_z, lsh0_456, lsg_325, \
                         lsg_448, lsg_449, lsh1_456, msg_445, msg_448, \
                         msg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_10 * lsg_448[k]
                   + f_3 * pc_x[k] * msg_448[k];

        t_623[k] = f_10 * lsg_449[k]
                   + f_3 * pc_x[k] * msg_449[k];

        t_624[k] = pa_z[k] * lsh0_456[k]
                   - f_8 * pc_z[k] * lsh1_456[k];

        t_625[k] = f_9 * lsg_325[k]
                   + f_3 * pc_z[k] * msg_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, lsg_342, lsg_343, lsg_344, msf0_298, \
                         msf0_299, msf1_298, msf1_299, msg_447, msg_448, \
                         msg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_16 * lsg_342[k]
                   + f_6 * msf0_298[k]
                   - f_7 * msf1_298[k]
                   + f_3 * pc_y[k] * msg_447[k];

        t_627[k] = f_16 * lsg_343[k]
                   + f_4 * msf0_299[k]
                   - f_5 * msf1_299[k]
                   + f_3 * pc_y[k] * msg_448[k];

        t_628[k] = f_16 * lsg_344[k]
                   + f_3 * pc_y[k] * msg_449[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, pc_z, lsg_329, lsg_345, lsg_450, \
                         msf0_299, msf0_300, msf1_299, msf1_300, msg_449, \
                         msg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_9 * lsg_329[k]
                   + f_1 * msf0_299[k]
                   - f_2 * msf1_299[k]
                   + f_3 * pc_z[k] * msg_449[k];

        t_630[k] = f_10 * lsg_450[k]
                   + f_1 * msf0_300[k]
                   - f_2 * msf1_300[k]
                   + f_3 * pc_x[k] * msg_450[k];

        t_631[k] = f_17 * lsg_345[k]
                   + f_3 * pc_y[k] * msg_450[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pc_x, pc_y, pc_z, lsg_330, lsg_347, lsg_453, \
                         msf0_303, msf1_303, msg_450, msg_452, \
                         msg_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_10 * lsg_330[k]
                   + f_3 * pc_z[k] * msg_450[k];

        t_633[k] = f_10 * lsg_453[k]
                   + f_6 * msf0_303[k]
                   - f_7 * msf1_303[k]
                   + f_3 * pc_x[k] * msg_453[k];

        t_634[k] = f_17 * lsg_347[k]
                   + f_3 * pc_y[k] * msg_452[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pc_x, pc_z, lsg_333, lsg_455, lsg_456, msf0_305, \
                         msf0_306, msf1_305, msf1_306, msg_453, msg_455, \
                         msg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_10 * lsg_455[k]
                   + f_6 * msf0_305[k]
                   - f_7 * msf1_305[k]
                   + f_3 * pc_x[k] * msg_455[k];

        t_636[k] = f_10 * lsg_456[k]
                   + f_4 * msf0_306[k]
                   - f_5 * msf1_306[k]
                   + f_3 * pc_x[k] * msg_456[k];

        t_637[k] = f_10 * lsg_333[k]
                   + f_3 * pc_z[k] * msg_453[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pc_x, pc_y, lsg_350, lsg_459, lsg_460, \
                         lsg_461, msf0_309, msf1_309, msg_455, msg_459, msg_460, \
                         msg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_17 * lsg_350[k]
                   + f_3 * pc_y[k] * msg_455[k];

        t_639[k] = f_10 * lsg_459[k]
                   + f_4 * msf0_309[k]
                   - f_5 * msf1_309[k]
                   + f_3 * pc_x[k] * msg_459[k];

        t_640[k] = f_10 * lsg_460[k]
                   + f_3 * pc_x[k] * msg_460[k];

        t_641[k] = f_10 * lsg_461[k]
                   + f_3 * pc_x[k] * msg_461[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, lsg_355, lsg_462, lsg_463, \
                         lsg_464, msf0_306, msf1_306, msg_460, msg_462, msg_463, \
                         msg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_10 * lsg_462[k]
                   + f_3 * pc_x[k] * msg_462[k];

        t_643[k] = f_10 * lsg_463[k]
                   + f_3 * pc_x[k] * msg_463[k];

        t_644[k] = f_10 * lsg_464[k]
                   + f_3 * pc_x[k] * msg_464[k];

        t_645[k] = f_17 * lsg_355[k]
                   + f_1 * msf0_306[k]
                   - f_2 * msf1_306[k]
                   + f_3 * pc_y[k] * msg_460[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_y, pc_z, lsg_340, lsg_357, lsg_358, msf0_308, \
                         msf0_309, msf1_308, msf1_309, msg_460, msg_462, \
                         msg_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_10 * lsg_340[k]
                   + f_3 * pc_z[k] * msg_460[k];

        t_647[k] = f_17 * lsg_357[k]
                   + f_6 * msf0_308[k]
                   - f_7 * msf1_308[k]
                   + f_3 * pc_y[k] * msg_462[k];

        t_648[k] = f_17 * lsg_358[k]
                   + f_4 * msf0_309[k]
                   - f_5 * msf1_309[k]
                   + f_3 * pc_y[k] * msg_463[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pc_x, pc_y, pc_z, lsg_344, lsg_359, lsg_465, \
                         msf0_309, msf0_310, msf1_309, msf1_310, msg_464, \
                         msg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_17 * lsg_359[k]
                   + f_3 * pc_y[k] * msg_464[k];

        t_650[k] = f_10 * lsg_344[k]
                   + f_1 * msf0_309[k]
                   - f_2 * msf1_309[k]
                   + f_3 * pc_z[k] * msg_464[k];

        t_651[k] = f_10 * lsg_465[k]
                   + f_1 * msf0_310[k]
                   - f_2 * msf1_310[k]
                   + f_3 * pc_x[k] * msg_465[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, pc_x, pc_y, pc_z, lsg_345, lsg_360, \
                         lsg_362, lsg_468, msf0_313, msf1_313, msg_465, msg_467, \
                         msg_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_18 * lsg_360[k]
                   + f_3 * pc_y[k] * msg_465[k];

        t_653[k] = f_11 * lsg_345[k]
                   + f_3 * pc_z[k] * msg_465[k];

        t_654[k] = f_10 * lsg_468[k]
                   + f_6 * msf0_313[k]
                   - f_7 * msf1_313[k]
                   + f_3 * pc_x[k] * msg_468[k];

        t_655[k] = f_18 * lsg_362[k]
                   + f_3 * pc_y[k] * msg_467[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_z, lsg_348, lsg_470, lsg_471, msf0_315, \
                         msf0_316, msf1_315, msf1_316, msg_468, msg_470, \
                         msg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_10 * lsg_470[k]
                   + f_6 * msf0_315[k]
                   - f_7 * msf1_315[k]
                   + f_3 * pc_x[k] * msg_470[k];

        t_657[k] = f_10 * lsg_471[k]
                   + f_4 * msf0_316[k]
                   - f_5 * msf1_316[k]
                   + f_3 * pc_x[k] * msg_471[k];

        t_658[k] = f_11 * lsg_348[k]
                   + f_3 * pc_z[k] * msg_468[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pc_x, pc_y, lsg_365, lsg_474, lsg_475, \
                         lsg_476, msf0_319, msf1_319, msg_470, msg_474, msg_475, \
                         msg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_18 * lsg_365[k]
                   + f_3 * pc_y[k] * msg_470[k];

        t_660[k] = f_10 * lsg_474[k]
                   + f_4 * msf0_319[k]
                   - f_5 * msf1_319[k]
                   + f_3 * pc_x[k] * msg_474[k];

        t_661[k] = f_10 * lsg_475[k]
                   + f_3 * pc_x[k] * msg_475[k];

        t_662[k] = f_10 * lsg_476[k]
                   + f_3 * pc_x[k] * msg_476[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, pc_x, pc_y, lsg_370, lsg_477, lsg_478, \
                         lsg_479, msf0_316, msf1_316, msg_475, msg_477, msg_478, \
                         msg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_10 * lsg_477[k]
                   + f_3 * pc_x[k] * msg_477[k];

        t_664[k] = f_10 * lsg_478[k]
                   + f_3 * pc_x[k] * msg_478[k];

        t_665[k] = f_10 * lsg_479[k]
                   + f_3 * pc_x[k] * msg_479[k];

        t_666[k] = f_18 * lsg_370[k]
                   + f_1 * msf0_316[k]
                   - f_2 * msf1_316[k]
                   + f_3 * pc_y[k] * msg_475[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, pc_z, lsg_355, lsg_372, lsg_373, msf0_318, \
                         msf0_319, msf1_318, msf1_319, msg_475, msg_477, \
                         msg_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_11 * lsg_355[k]
                   + f_3 * pc_z[k] * msg_475[k];

        t_668[k] = f_18 * lsg_372[k]
                   + f_6 * msf0_318[k]
                   - f_7 * msf1_318[k]
                   + f_3 * pc_y[k] * msg_477[k];

        t_669[k] = f_18 * lsg_373[k]
                   + f_4 * msf0_319[k]
                   - f_5 * msf1_319[k]
                   + f_3 * pc_y[k] * msg_478[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, lsg_359, lsg_374, lsg_480, \
                         msf0_319, msf0_320, msf1_319, msf1_320, msg_479, \
                         msg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_18 * lsg_374[k]
                   + f_3 * pc_y[k] * msg_479[k];

        t_671[k] = f_11 * lsg_359[k]
                   + f_1 * msf0_319[k]
                   - f_2 * msf1_319[k]
                   + f_3 * pc_z[k] * msg_479[k];

        t_672[k] = f_10 * lsg_480[k]
                   + f_1 * msf0_320[k]
                   - f_2 * msf1_320[k]
                   + f_3 * pc_x[k] * msg_480[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, lsg_360, lsg_375, \
                         lsg_377, lsg_483, msf0_323, msf1_323, msg_480, msg_482, \
                         msg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_11 * lsg_375[k]
                   + f_3 * pc_y[k] * msg_480[k];

        t_674[k] = f_18 * lsg_360[k]
                   + f_3 * pc_z[k] * msg_480[k];

        t_675[k] = f_10 * lsg_483[k]
                   + f_6 * msf0_323[k]
                   - f_7 * msf1_323[k]
                   + f_3 * pc_x[k] * msg_483[k];

        t_676[k] = f_11 * lsg_377[k]
                   + f_3 * pc_y[k] * msg_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, lsg_363, lsg_485, lsg_486, msf0_325, \
                         msf0_326, msf1_325, msf1_326, msg_483, msg_485, \
                         msg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_10 * lsg_485[k]
                   + f_6 * msf0_325[k]
                   - f_7 * msf1_325[k]
                   + f_3 * pc_x[k] * msg_485[k];

        t_678[k] = f_10 * lsg_486[k]
                   + f_4 * msf0_326[k]
                   - f_5 * msf1_326[k]
                   + f_3 * pc_x[k] * msg_486[k];

        t_679[k] = f_18 * lsg_363[k]
                   + f_3 * pc_z[k] * msg_483[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_y, lsg_380, lsg_489, lsg_490, \
                         lsg_491, msf0_329, msf1_329, msg_485, msg_489, msg_490, \
                         msg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * lsg_380[k]
                   + f_3 * pc_y[k] * msg_485[k];

        t_681[k] = f_10 * lsg_489[k]
                   + f_4 * msf0_329[k]
                   - f_5 * msf1_329[k]
                   + f_3 * pc_x[k] * msg_489[k];

        t_682[k] = f_10 * lsg_490[k]
                   + f_3 * pc_x[k] * msg_490[k];

        t_683[k] = f_10 * lsg_491[k]
                   + f_3 * pc_x[k] * msg_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, lsg_385, lsg_492, lsg_493, \
                         lsg_494, msf0_326, msf1_326, msg_490, msg_492, msg_493, \
                         msg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_10 * lsg_492[k]
                   + f_3 * pc_x[k] * msg_492[k];

        t_685[k] = f_10 * lsg_493[k]
                   + f_3 * pc_x[k] * msg_493[k];

        t_686[k] = f_10 * lsg_494[k]
                   + f_3 * pc_x[k] * msg_494[k];

        t_687[k] = f_11 * lsg_385[k]
                   + f_1 * msf0_326[k]
                   - f_2 * msf1_326[k]
                   + f_3 * pc_y[k] * msg_490[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pc_y, pc_z, lsg_370, lsg_387, lsg_388, msf0_328, \
                         msf0_329, msf1_328, msf1_329, msg_490, msg_492, \
                         msg_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_18 * lsg_370[k]
                   + f_3 * pc_z[k] * msg_490[k];

        t_689[k] = f_11 * lsg_387[k]
                   + f_6 * msf0_328[k]
                   - f_7 * msf1_328[k]
                   + f_3 * pc_y[k] * msg_492[k];

        t_690[k] = f_11 * lsg_388[k]
                   + f_4 * msf0_329[k]
                   - f_5 * msf1_329[k]
                   + f_3 * pc_y[k] * msg_493[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pc_x, pc_y, pc_z, lsg_374, lsg_389, lsg_495, \
                         msf0_329, msf0_330, msf1_329, msf1_330, msg_494, \
                         msg_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_11 * lsg_389[k]
                   + f_3 * pc_y[k] * msg_494[k];

        t_692[k] = f_18 * lsg_374[k]
                   + f_1 * msf0_329[k]
                   - f_2 * msf1_329[k]
                   + f_3 * pc_z[k] * msg_494[k];

        t_693[k] = f_10 * lsg_495[k]
                   + f_1 * msf0_330[k]
                   - f_2 * msf1_330[k]
                   + f_3 * pc_x[k] * msg_495[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pc_x, pc_y, pc_z, lsg_375, lsg_390, \
                         lsg_392, lsg_498, msf0_333, msf1_333, msg_495, msg_497, \
                         msg_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_10 * lsg_390[k]
                   + f_3 * pc_y[k] * msg_495[k];

        t_695[k] = f_17 * lsg_375[k]
                   + f_3 * pc_z[k] * msg_495[k];

        t_696[k] = f_10 * lsg_498[k]
                   + f_6 * msf0_333[k]
                   - f_7 * msf1_333[k]
                   + f_3 * pc_x[k] * msg_498[k];

        t_697[k] = f_10 * lsg_392[k]
                   + f_3 * pc_y[k] * msg_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_z, lsg_378, lsg_500, lsg_501, msf0_335, \
                         msf0_336, msf1_335, msf1_336, msg_498, msg_500, \
                         msg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_10 * lsg_500[k]
                   + f_6 * msf0_335[k]
                   - f_7 * msf1_335[k]
                   + f_3 * pc_x[k] * msg_500[k];

        t_699[k] = f_10 * lsg_501[k]
                   + f_4 * msf0_336[k]
                   - f_5 * msf1_336[k]
                   + f_3 * pc_x[k] * msg_501[k];

        t_700[k] = f_17 * lsg_378[k]
                   + f_3 * pc_z[k] * msg_498[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pc_x, pc_y, lsg_395, lsg_504, lsg_505, \
                         lsg_506, msf0_339, msf1_339, msg_500, msg_504, msg_505, \
                         msg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * lsg_395[k]
                   + f_3 * pc_y[k] * msg_500[k];

        t_702[k] = f_10 * lsg_504[k]
                   + f_4 * msf0_339[k]
                   - f_5 * msf1_339[k]
                   + f_3 * pc_x[k] * msg_504[k];

        t_703[k] = f_10 * lsg_505[k]
                   + f_3 * pc_x[k] * msg_505[k];

        t_704[k] = f_10 * lsg_506[k]
                   + f_3 * pc_x[k] * msg_506[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_y, lsg_400, lsg_507, lsg_508, \
                         lsg_509, msf0_336, msf1_336, msg_505, msg_507, msg_508, \
                         msg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_10 * lsg_507[k]
                   + f_3 * pc_x[k] * msg_507[k];

        t_706[k] = f_10 * lsg_508[k]
                   + f_3 * pc_x[k] * msg_508[k];

        t_707[k] = f_10 * lsg_509[k]
                   + f_3 * pc_x[k] * msg_509[k];

        t_708[k] = f_10 * lsg_400[k]
                   + f_1 * msf0_336[k]
                   - f_2 * msf1_336[k]
                   + f_3 * pc_y[k] * msg_505[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_y, pc_z, lsg_385, lsg_402, lsg_403, msf0_338, \
                         msf0_339, msf1_338, msf1_339, msg_505, msg_507, \
                         msg_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_17 * lsg_385[k]
                   + f_3 * pc_z[k] * msg_505[k];

        t_710[k] = f_10 * lsg_402[k]
                   + f_6 * msf0_338[k]
                   - f_7 * msf1_338[k]
                   + f_3 * pc_y[k] * msg_507[k];

        t_711[k] = f_10 * lsg_403[k]
                   + f_4 * msf0_339[k]
                   - f_5 * msf1_339[k]
                   + f_3 * pc_y[k] * msg_508[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_12 = 4.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_567 = buffer.data(lsh0 + 567);
    const auto *lsh0_570 = buffer.data(lsh0 + 570);
    const auto *lsh0_572 = buffer.data(lsh0 + 572);
    const auto *lsh0_573 = buffer.data(lsh0 + 573);
    const auto *lsh0_576 = buffer.data(lsh0 + 576);
    const auto *lsh0_587 = buffer.data(lsh0 + 587);
    const auto *lsh0_588 = buffer.data(lsh0 + 588);
    const auto *lsh0_591 = buffer.data(lsh0 + 591);
    const auto *lsh0_594 = buffer.data(lsh0 + 594);
    const auto *lsh0_756 = buffer.data(lsh0 + 756);
    const auto *lsh0_759 = buffer.data(lsh0 + 759);
    const auto *lsh0_762 = buffer.data(lsh0 + 762);
    const auto *lsh0_771 = buffer.data(lsh0 + 771);
    const auto *lsh0_773 = buffer.data(lsh0 + 773);
    const auto *lsh0_774 = buffer.data(lsh0 + 774);
    const auto *lsh0_776 = buffer.data(lsh0 + 776);
    const auto *lsh0_782 = buffer.data(lsh0 + 782);
    const auto *lsh0_786 = buffer.data(lsh0 + 786);
    const auto *lsh0_792 = buffer.data(lsh0 + 792);
    const auto *lsh0_794 = buffer.data(lsh0 + 794);
    const auto *lsh0_795 = buffer.data(lsh0 + 795);
    const auto *lsh0_797 = buffer.data(lsh0 + 797);
    const auto *lsh0_798 = buffer.data(lsh0 + 798);
    const auto *lsh0_801 = buffer.data(lsh0 + 801);
    const auto *lsh0_803 = buffer.data(lsh0 + 803);
    const auto *lsh0_804 = buffer.data(lsh0 + 804);
    const auto *lsh0_807 = buffer.data(lsh0 + 807);
    const auto *lsh0_813 = buffer.data(lsh0 + 813);
    const auto *lsh0_815 = buffer.data(lsh0 + 815);
    const auto *lsh0_816 = buffer.data(lsh0 + 816);
    const auto *lsh0_818 = buffer.data(lsh0 + 818);
    const auto *lsh0_819 = buffer.data(lsh0 + 819);
    const auto *lsh0_822 = buffer.data(lsh0 + 822);
    const auto *lsh0_824 = buffer.data(lsh0 + 824);
    const auto *lsh0_825 = buffer.data(lsh0 + 825);
    const auto *lsh0_828 = buffer.data(lsh0 + 828);
    const auto *lsh0_834 = buffer.data(lsh0 + 834);

    const auto *lsg_389 = buffer.data(lsg + 389);
    const auto *lsg_390 = buffer.data(lsg + 390);
    const auto *lsg_393 = buffer.data(lsg + 393);
    const auto *lsg_400 = buffer.data(lsg + 400);
    const auto *lsg_404 = buffer.data(lsg + 404);
    const auto *lsg_405 = buffer.data(lsg + 405);
    const auto *lsg_406 = buffer.data(lsg + 406);
    const auto *lsg_407 = buffer.data(lsg + 407);
    const auto *lsg_408 = buffer.data(lsg + 408);
    const auto *lsg_410 = buffer.data(lsg + 410);
    const auto *lsg_415 = buffer.data(lsg + 415);
    const auto *lsg_417 = buffer.data(lsg + 417);
    const auto *lsg_418 = buffer.data(lsg + 418);
    const auto *lsg_419 = buffer.data(lsg + 419);
    const auto *lsg_420 = buffer.data(lsg + 420);
    const auto *lsg_423 = buffer.data(lsg + 423);
    const auto *lsg_425 = buffer.data(lsg + 425);
    const auto *lsg_430 = buffer.data(lsg + 430);
    const auto *lsg_434 = buffer.data(lsg + 434);
    const auto *lsg_435 = buffer.data(lsg + 435);
    const auto *lsg_437 = buffer.data(lsg + 437);
    const auto *lsg_438 = buffer.data(lsg + 438);
    const auto *lsg_440 = buffer.data(lsg + 440);
    const auto *lsg_445 = buffer.data(lsg + 445);
    const auto *lsg_449 = buffer.data(lsg + 449);
    const auto *lsg_450 = buffer.data(lsg + 450);
    const auto *lsg_452 = buffer.data(lsg + 452);
    const auto *lsg_453 = buffer.data(lsg + 453);
    const auto *lsg_455 = buffer.data(lsg + 455);
    const auto *lsg_460 = buffer.data(lsg + 460);
    const auto *lsg_464 = buffer.data(lsg + 464);
    const auto *lsg_465 = buffer.data(lsg + 465);
    const auto *lsg_467 = buffer.data(lsg + 467);
    const auto *lsg_470 = buffer.data(lsg + 470);
    const auto *lsg_520 = buffer.data(lsg + 520);
    const auto *lsg_521 = buffer.data(lsg + 521);
    const auto *lsg_522 = buffer.data(lsg + 522);
    const auto *lsg_523 = buffer.data(lsg + 523);
    const auto *lsg_524 = buffer.data(lsg + 524);
    const auto *lsg_525 = buffer.data(lsg + 525);
    const auto *lsg_530 = buffer.data(lsg + 530);
    const auto *lsg_534 = buffer.data(lsg + 534);
    const auto *lsg_535 = buffer.data(lsg + 535);
    const auto *lsg_536 = buffer.data(lsg + 536);
    const auto *lsg_537 = buffer.data(lsg + 537);
    const auto *lsg_539 = buffer.data(lsg + 539);
    const auto *lsg_540 = buffer.data(lsg + 540);
    const auto *lsg_543 = buffer.data(lsg + 543);
    const auto *lsg_546 = buffer.data(lsg + 546);
    const auto *lsg_550 = buffer.data(lsg + 550);
    const auto *lsg_552 = buffer.data(lsg + 552);
    const auto *lsg_553 = buffer.data(lsg + 553);
    const auto *lsg_554 = buffer.data(lsg + 554);
    const auto *lsg_560 = buffer.data(lsg + 560);
    const auto *lsg_564 = buffer.data(lsg + 564);
    const auto *lsg_565 = buffer.data(lsg + 565);
    const auto *lsg_566 = buffer.data(lsg + 566);
    const auto *lsg_567 = buffer.data(lsg + 567);
    const auto *lsg_568 = buffer.data(lsg + 568);
    const auto *lsg_569 = buffer.data(lsg + 569);
    const auto *lsg_570 = buffer.data(lsg + 570);
    const auto *lsg_573 = buffer.data(lsg + 573);
    const auto *lsg_575 = buffer.data(lsg + 575);
    const auto *lsg_576 = buffer.data(lsg + 576);
    const auto *lsg_579 = buffer.data(lsg + 579);
    const auto *lsg_580 = buffer.data(lsg + 580);
    const auto *lsg_581 = buffer.data(lsg + 581);
    const auto *lsg_582 = buffer.data(lsg + 582);
    const auto *lsg_583 = buffer.data(lsg + 583);
    const auto *lsg_584 = buffer.data(lsg + 584);
    const auto *lsg_585 = buffer.data(lsg + 585);
    const auto *lsg_588 = buffer.data(lsg + 588);
    const auto *lsg_590 = buffer.data(lsg + 590);
    const auto *lsg_591 = buffer.data(lsg + 591);
    const auto *lsg_594 = buffer.data(lsg + 594);
    const auto *lsg_595 = buffer.data(lsg + 595);
    const auto *lsg_596 = buffer.data(lsg + 596);
    const auto *lsg_597 = buffer.data(lsg + 597);
    const auto *lsg_598 = buffer.data(lsg + 598);
    const auto *lsg_599 = buffer.data(lsg + 599);

    const auto *lsh1_567 = buffer.data(lsh1 + 567);
    const auto *lsh1_570 = buffer.data(lsh1 + 570);
    const auto *lsh1_572 = buffer.data(lsh1 + 572);
    const auto *lsh1_573 = buffer.data(lsh1 + 573);
    const auto *lsh1_576 = buffer.data(lsh1 + 576);
    const auto *lsh1_587 = buffer.data(lsh1 + 587);
    const auto *lsh1_588 = buffer.data(lsh1 + 588);
    const auto *lsh1_591 = buffer.data(lsh1 + 591);
    const auto *lsh1_594 = buffer.data(lsh1 + 594);
    const auto *lsh1_756 = buffer.data(lsh1 + 756);
    const auto *lsh1_759 = buffer.data(lsh1 + 759);
    const auto *lsh1_762 = buffer.data(lsh1 + 762);
    const auto *lsh1_771 = buffer.data(lsh1 + 771);
    const auto *lsh1_773 = buffer.data(lsh1 + 773);
    const auto *lsh1_774 = buffer.data(lsh1 + 774);
    const auto *lsh1_776 = buffer.data(lsh1 + 776);
    const auto *lsh1_782 = buffer.data(lsh1 + 782);
    const auto *lsh1_786 = buffer.data(lsh1 + 786);
    const auto *lsh1_792 = buffer.data(lsh1 + 792);
    const auto *lsh1_794 = buffer.data(lsh1 + 794);
    const auto *lsh1_795 = buffer.data(lsh1 + 795);
    const auto *lsh1_797 = buffer.data(lsh1 + 797);
    const auto *lsh1_798 = buffer.data(lsh1 + 798);
    const auto *lsh1_801 = buffer.data(lsh1 + 801);
    const auto *lsh1_803 = buffer.data(lsh1 + 803);
    const auto *lsh1_804 = buffer.data(lsh1 + 804);
    const auto *lsh1_807 = buffer.data(lsh1 + 807);
    const auto *lsh1_813 = buffer.data(lsh1 + 813);
    const auto *lsh1_815 = buffer.data(lsh1 + 815);
    const auto *lsh1_816 = buffer.data(lsh1 + 816);
    const auto *lsh1_818 = buffer.data(lsh1 + 818);
    const auto *lsh1_819 = buffer.data(lsh1 + 819);
    const auto *lsh1_822 = buffer.data(lsh1 + 822);
    const auto *lsh1_824 = buffer.data(lsh1 + 824);
    const auto *lsh1_825 = buffer.data(lsh1 + 825);
    const auto *lsh1_828 = buffer.data(lsh1 + 828);
    const auto *lsh1_834 = buffer.data(lsh1 + 834);

    const auto *msf0_339 = buffer.data(msf0 + 339);
    const auto *msf0_346 = buffer.data(msf0 + 346);
    const auto *msf0_348 = buffer.data(msf0 + 348);
    const auto *msf0_349 = buffer.data(msf0 + 349);
    const auto *msf0_350 = buffer.data(msf0 + 350);
    const auto *msf0_351 = buffer.data(msf0 + 351);
    const auto *msf0_352 = buffer.data(msf0 + 352);
    const auto *msf0_355 = buffer.data(msf0 + 355);
    const auto *msf0_356 = buffer.data(msf0 + 356);
    const auto *msf0_357 = buffer.data(msf0 + 357);
    const auto *msf0_358 = buffer.data(msf0 + 358);
    const auto *msf0_359 = buffer.data(msf0 + 359);
    const auto *msf0_360 = buffer.data(msf0 + 360);
    const auto *msf0_362 = buffer.data(msf0 + 362);

    const auto *msf1_339 = buffer.data(msf1 + 339);
    const auto *msf1_346 = buffer.data(msf1 + 346);
    const auto *msf1_348 = buffer.data(msf1 + 348);
    const auto *msf1_349 = buffer.data(msf1 + 349);
    const auto *msf1_350 = buffer.data(msf1 + 350);
    const auto *msf1_351 = buffer.data(msf1 + 351);
    const auto *msf1_352 = buffer.data(msf1 + 352);
    const auto *msf1_355 = buffer.data(msf1 + 355);
    const auto *msf1_356 = buffer.data(msf1 + 356);
    const auto *msf1_357 = buffer.data(msf1 + 357);
    const auto *msf1_358 = buffer.data(msf1 + 358);
    const auto *msf1_359 = buffer.data(msf1 + 359);
    const auto *msf1_360 = buffer.data(msf1 + 360);
    const auto *msf1_362 = buffer.data(msf1 + 362);

    const auto *msg_509 = buffer.data(msg + 509);
    const auto *msg_510 = buffer.data(msg + 510);
    const auto *msg_512 = buffer.data(msg + 512);
    const auto *msg_513 = buffer.data(msg + 513);
    const auto *msg_515 = buffer.data(msg + 515);
    const auto *msg_520 = buffer.data(msg + 520);
    const auto *msg_521 = buffer.data(msg + 521);
    const auto *msg_522 = buffer.data(msg + 522);
    const auto *msg_523 = buffer.data(msg + 523);
    const auto *msg_524 = buffer.data(msg + 524);
    const auto *msg_525 = buffer.data(msg + 525);
    const auto *msg_526 = buffer.data(msg + 526);
    const auto *msg_527 = buffer.data(msg + 527);
    const auto *msg_528 = buffer.data(msg + 528);
    const auto *msg_529 = buffer.data(msg + 529);
    const auto *msg_530 = buffer.data(msg + 530);
    const auto *msg_534 = buffer.data(msg + 534);
    const auto *msg_535 = buffer.data(msg + 535);
    const auto *msg_536 = buffer.data(msg + 536);
    const auto *msg_537 = buffer.data(msg + 537);
    const auto *msg_538 = buffer.data(msg + 538);
    const auto *msg_539 = buffer.data(msg + 539);
    const auto *msg_540 = buffer.data(msg + 540);
    const auto *msg_541 = buffer.data(msg + 541);
    const auto *msg_542 = buffer.data(msg + 542);
    const auto *msg_543 = buffer.data(msg + 543);
    const auto *msg_545 = buffer.data(msg + 545);
    const auto *msg_546 = buffer.data(msg + 546);
    const auto *msg_550 = buffer.data(msg + 550);
    const auto *msg_552 = buffer.data(msg + 552);
    const auto *msg_553 = buffer.data(msg + 553);
    const auto *msg_554 = buffer.data(msg + 554);
    const auto *msg_555 = buffer.data(msg + 555);
    const auto *msg_557 = buffer.data(msg + 557);
    const auto *msg_558 = buffer.data(msg + 558);
    const auto *msg_560 = buffer.data(msg + 560);
    const auto *msg_565 = buffer.data(msg + 565);
    const auto *msg_566 = buffer.data(msg + 566);
    const auto *msg_567 = buffer.data(msg + 567);
    const auto *msg_568 = buffer.data(msg + 568);
    const auto *msg_569 = buffer.data(msg + 569);
    const auto *msg_570 = buffer.data(msg + 570);
    const auto *msg_572 = buffer.data(msg + 572);
    const auto *msg_573 = buffer.data(msg + 573);
    const auto *msg_575 = buffer.data(msg + 575);
    const auto *msg_580 = buffer.data(msg + 580);
    const auto *msg_581 = buffer.data(msg + 581);
    const auto *msg_582 = buffer.data(msg + 582);
    const auto *msg_583 = buffer.data(msg + 583);
    const auto *msg_584 = buffer.data(msg + 584);
    const auto *msg_585 = buffer.data(msg + 585);
    const auto *msg_587 = buffer.data(msg + 587);
    const auto *msg_588 = buffer.data(msg + 588);
    const auto *msg_590 = buffer.data(msg + 590);
    const auto *msg_595 = buffer.data(msg + 595);
    const auto *msg_596 = buffer.data(msg + 596);
    const auto *msg_597 = buffer.data(msg + 597);
    const auto *msg_598 = buffer.data(msg + 598);
    const auto *msg_599 = buffer.data(msg + 599);

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_y, pc_y, pc_z, lsh0_567, lsg_389, \
                         lsg_404, lsg_405, lsh1_567, msf0_339, msf1_339, msg_509, \
                         msg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_10 * lsg_404[k]
                   + f_3 * pc_y[k] * msg_509[k];

        t_713[k] = f_17 * lsg_389[k]
                   + f_1 * msf0_339[k]
                   - f_2 * msf1_339[k]
                   + f_3 * pc_z[k] * msg_509[k];

        t_714[k] = pa_y[k] * lsh0_567[k]
                   - f_8 * pc_y[k] * lsh1_567[k];

        t_715[k] = f_9 * lsg_405[k]
                   + f_3 * pc_y[k] * msg_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pa_y, pc_y, pc_z, lsh0_570, lsh0_572, \
                         lsg_390, lsg_406, lsg_407, lsh1_570, lsh1_572, msg_510, \
                         msg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_16 * lsg_390[k]
                   + f_3 * pc_z[k] * msg_510[k];

        t_717[k] = pa_y[k] * lsh0_570[k]
                   + f_10 * lsg_406[k]
                   - f_8 * pc_y[k] * lsh1_570[k];

        t_718[k] = f_9 * lsg_407[k]
                   + f_3 * pc_y[k] * msg_512[k];

        t_719[k] = pa_y[k] * lsh0_572[k]
                   - f_8 * pc_y[k] * lsh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_y, pc_y, pc_z, lsh0_573, lsh0_576, \
                         lsg_393, lsg_408, lsg_410, lsh1_573, lsh1_576, msg_513, \
                         msg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_y[k] * lsh0_573[k]
                   + f_11 * lsg_408[k]
                   - f_8 * pc_y[k] * lsh1_573[k];

        t_721[k] = f_16 * lsg_393[k]
                   + f_3 * pc_z[k] * msg_513[k];

        t_722[k] = f_9 * lsg_410[k]
                   + f_3 * pc_y[k] * msg_515[k];

        t_723[k] = pa_y[k] * lsh0_576[k]
                   - f_8 * pc_y[k] * lsh1_576[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, lsg_520, lsg_521, lsg_522, \
                         lsg_523, lsg_524, msg_520, msg_521, msg_522, msg_523, \
                         msg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_10 * lsg_520[k]
                   + f_3 * pc_x[k] * msg_520[k];

        t_725[k] = f_10 * lsg_521[k]
                   + f_3 * pc_x[k] * msg_521[k];

        t_726[k] = f_10 * lsg_522[k]
                   + f_3 * pc_x[k] * msg_522[k];

        t_727[k] = f_10 * lsg_523[k]
                   + f_3 * pc_x[k] * msg_523[k];

        t_728[k] = f_10 * lsg_524[k]
                   + f_3 * pc_x[k] * msg_524[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_y, pc_z, lsg_400, lsg_415, lsg_417, msf0_346, \
                         msf0_348, msf1_346, msf1_348, msg_520, \
                         msg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_9 * lsg_415[k]
                   + f_1 * msf0_346[k]
                   - f_2 * msf1_346[k]
                   + f_3 * pc_y[k] * msg_520[k];

        t_730[k] = f_16 * lsg_400[k]
                   + f_3 * pc_z[k] * msg_520[k];

        t_731[k] = f_9 * lsg_417[k]
                   + f_6 * msf0_348[k]
                   - f_7 * msf1_348[k]
                   + f_3 * pc_y[k] * msg_522[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pa_y, pc_y, lsh0_587, lsg_418, lsg_419, \
                         lsh1_587, msf0_349, msf1_349, msg_523, \
                         msg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_9 * lsg_418[k]
                   + f_4 * msf0_349[k]
                   - f_5 * msf1_349[k]
                   + f_3 * pc_y[k] * msg_523[k];

        t_733[k] = f_9 * lsg_419[k]
                   + f_3 * pc_y[k] * msg_524[k];

        t_734[k] = pa_y[k] * lsh0_587[k]
                   - f_8 * pc_y[k] * lsh1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, pc_x, pc_y, pc_z, lsg_405, \
                         lsg_525, msf0_350, msf1_350, msg_525, msg_526, \
                         msg_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_10 * lsg_525[k]
                   + f_1 * msf0_350[k]
                   - f_2 * msf1_350[k]
                   + f_3 * pc_x[k] * msg_525[k];

        t_736[k] = f_3 * pc_y[k] * msg_525[k];

        t_737[k] = f_15 * lsg_405[k]
                   + f_3 * pc_z[k] * msg_525[k];

        t_738[k] = f_4 * msf0_350[k]
                   - f_5 * msf1_350[k]
                   + f_3 * pc_y[k] * msg_526[k];

        t_739[k] = f_3 * pc_y[k] * msg_527[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, pc_y, lsg_530, msf0_351, msf0_352, \
                         msf0_355, msf1_351, msf1_352, msf1_355, msg_528, msg_529, \
                         msg_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_10 * lsg_530[k]
                   + f_6 * msf0_355[k]
                   - f_7 * msf1_355[k]
                   + f_3 * pc_x[k] * msg_530[k];

        t_741[k] = f_6 * msf0_351[k]
                   - f_7 * msf1_351[k]
                   + f_3 * pc_y[k] * msg_528[k];

        t_742[k] = f_4 * msf0_352[k]
                   - f_5 * msf1_352[k]
                   + f_3 * pc_y[k] * msg_529[k];

        t_743[k] = f_3 * pc_y[k] * msg_530[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, pc_x, lsg_534, lsg_535, lsg_536, lsg_537, \
                         msf0_359, msf1_359, msg_534, msg_535, msg_536, \
                         msg_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_10 * lsg_534[k]
                   + f_4 * msf0_359[k]
                   - f_5 * msf1_359[k]
                   + f_3 * pc_x[k] * msg_534[k];

        t_745[k] = f_10 * lsg_535[k]
                   + f_3 * pc_x[k] * msg_535[k];

        t_746[k] = f_10 * lsg_536[k]
                   + f_3 * pc_x[k] * msg_536[k];

        t_747[k] = f_10 * lsg_537[k]
                   + f_3 * pc_x[k] * msg_537[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_x, pc_y, lsg_539, msf0_356, msf0_357, \
                         msf1_356, msf1_357, msg_534, msg_535, msg_536, \
                         msg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_3 * pc_y[k] * msg_534[k];

        t_749[k] = f_10 * lsg_539[k]
                   + f_3 * pc_x[k] * msg_539[k];

        t_750[k] = f_1 * msf0_356[k]
                   - f_2 * msf1_356[k]
                   + f_3 * pc_y[k] * msg_535[k];

        t_751[k] = f_13 * msf0_357[k]
                   - f_14 * msf1_357[k]
                   + f_3 * pc_y[k] * msg_536[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pc_y, pc_z, lsg_419, msf0_358, msf0_359, \
                         msf1_358, msf1_359, msg_537, msg_538, \
                         msg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_6 * msf0_358[k]
                   - f_7 * msf1_358[k]
                   + f_3 * pc_y[k] * msg_537[k];

        t_753[k] = f_4 * msf0_359[k]
                   - f_5 * msf1_359[k]
                   + f_3 * pc_y[k] * msg_538[k];

        t_754[k] = f_3 * pc_y[k] * msg_539[k];

        t_755[k] = f_15 * lsg_419[k]
                   + f_1 * msf0_359[k]
                   - f_2 * msf1_359[k]
                   + f_3 * pc_z[k] * msg_539[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pa_x, pc_x, pc_y, pc_z, lsh0_756, \
                         lsh0_759, lsg_420, lsg_540, lsg_543, lsh1_756, lsh1_759, \
                         msg_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = pa_x[k] * lsh0_756[k]
                   + f_17 * lsg_540[k]
                   - f_8 * pc_x[k] * lsh1_756[k];

        t_757[k] = f_12 * lsg_420[k]
                   + f_3 * pc_y[k] * msg_540[k];

        t_758[k] = f_3 * pc_z[k] * msg_540[k];

        t_759[k] = pa_x[k] * lsh0_759[k]
                   + f_11 * lsg_543[k]
                   - f_8 * pc_x[k] * lsh1_759[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, pa_x, pc_x, pc_z, lsh0_762, lsg_546, \
                         lsh1_762, msf0_360, msf1_360, msg_541, msg_542, \
                         msg_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_3 * pc_z[k] * msg_541[k];

        t_761[k] = f_4 * msf0_360[k]
                   - f_5 * msf1_360[k]
                   + f_3 * pc_z[k] * msg_542[k];

        t_762[k] = pa_x[k] * lsh0_762[k]
                   + f_10 * lsg_546[k]
                   - f_8 * pc_x[k] * lsh1_762[k];

        t_763[k] = f_3 * pc_z[k] * msg_543[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, lsg_425, lsg_550, \
                         msf0_362, msf1_362, msg_545, msg_546, \
                         msg_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_12 * lsg_425[k]
                   + f_3 * pc_y[k] * msg_545[k];

        t_765[k] = f_6 * msf0_362[k]
                   - f_7 * msf1_362[k]
                   + f_3 * pc_z[k] * msg_545[k];

        t_766[k] = f_9 * lsg_550[k]
                   + f_3 * pc_x[k] * msg_550[k];

        t_767[k] = f_3 * pc_z[k] * msg_546[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pa_x, pc_x, lsh0_771, lsg_552, lsg_553, \
                         lsg_554, lsh1_771, msg_552, msg_553, msg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_9 * lsg_552[k]
                   + f_3 * pc_x[k] * msg_552[k];

        t_769[k] = f_9 * lsg_553[k]
                   + f_3 * pc_x[k] * msg_553[k];

        t_770[k] = f_9 * lsg_554[k]
                   + f_3 * pc_x[k] * msg_554[k];

        t_771[k] = pa_x[k] * lsh0_771[k]
                   - f_8 * pc_x[k] * lsh1_771[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, pa_x, pc_x, pc_y, pc_z, lsh0_773, \
                         lsh0_774, lsg_434, lsh1_773, lsh1_774, msg_550, \
                         msg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * msg_550[k];

        t_773[k] = pa_x[k] * lsh0_773[k]
                   - f_8 * pc_x[k] * lsh1_773[k];

        t_774[k] = pa_x[k] * lsh0_774[k]
                   - f_8 * pc_x[k] * lsh1_774[k];

        t_775[k] = f_12 * lsg_434[k]
                   + f_3 * pc_y[k] * msg_554[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pa_x, pa_z, pc_x, pc_y, pc_z, lsh0_588, \
                         lsh0_776, lsg_420, lsg_435, lsh1_588, lsh1_776, \
                         msg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = pa_x[k] * lsh0_776[k]
                   - f_8 * pc_x[k] * lsh1_776[k];

        t_777[k] = pa_z[k] * lsh0_588[k]
                   - f_8 * pc_z[k] * lsh1_588[k];

        t_778[k] = f_15 * lsg_435[k]
                   + f_3 * pc_y[k] * msg_555[k];

        t_779[k] = f_9 * lsg_420[k]
                   + f_3 * pc_z[k] * msg_555[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pa_x, pa_z, pc_x, pc_y, pc_z, lsh0_591, \
                         lsh0_782, lsg_437, lsg_560, lsh1_591, lsh1_782, \
                         msg_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = pa_z[k] * lsh0_591[k]
                   - f_8 * pc_z[k] * lsh1_591[k];

        t_781[k] = f_15 * lsg_437[k]
                   + f_3 * pc_y[k] * msg_557[k];

        t_782[k] = pa_x[k] * lsh0_782[k]
                   + f_11 * lsg_560[k]
                   - f_8 * pc_x[k] * lsh1_782[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pa_z, pc_y, pc_z, lsh0_594, lsg_423, lsg_440, \
                         lsh1_594, msg_558, msg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = pa_z[k] * lsh0_594[k]
                   - f_8 * pc_z[k] * lsh1_594[k];

        t_784[k] = f_9 * lsg_423[k]
                   + f_3 * pc_z[k] * msg_558[k];

        t_785[k] = f_15 * lsg_440[k]
                   + f_3 * pc_y[k] * msg_560[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pa_x, pc_x, lsh0_786, lsg_564, lsg_565, \
                         lsg_566, lsg_567, lsh1_786, msg_565, msg_566, \
                         msg_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = pa_x[k] * lsh0_786[k]
                   + f_10 * lsg_564[k]
                   - f_8 * pc_x[k] * lsh1_786[k];

        t_787[k] = f_9 * lsg_565[k]
                   + f_3 * pc_x[k] * msg_565[k];

        t_788[k] = f_9 * lsg_566[k]
                   + f_3 * pc_x[k] * msg_566[k];

        t_789[k] = f_9 * lsg_567[k]
                   + f_3 * pc_x[k] * msg_567[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_x, pc_x, pc_z, lsh0_792, lsg_430, \
                         lsg_568, lsg_569, lsh1_792, msg_565, msg_568, \
                         msg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_9 * lsg_568[k]
                   + f_3 * pc_x[k] * msg_568[k];

        t_791[k] = f_9 * lsg_569[k]
                   + f_3 * pc_x[k] * msg_569[k];

        t_792[k] = pa_x[k] * lsh0_792[k]
                   - f_8 * pc_x[k] * lsh1_792[k];

        t_793[k] = f_9 * lsg_430[k]
                   + f_3 * pc_z[k] * msg_565[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, pa_x, pc_x, pc_y, lsh0_794, lsh0_795, \
                         lsh0_797, lsg_449, lsh1_794, lsh1_795, lsh1_797, \
                         msg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = pa_x[k] * lsh0_794[k]
                   - f_8 * pc_x[k] * lsh1_794[k];

        t_795[k] = pa_x[k] * lsh0_795[k]
                   - f_8 * pc_x[k] * lsh1_795[k];

        t_796[k] = f_15 * lsg_449[k]
                   + f_3 * pc_y[k] * msg_569[k];

        t_797[k] = pa_x[k] * lsh0_797[k]
                   - f_8 * pc_x[k] * lsh1_797[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, pa_x, pc_x, pc_y, pc_z, lsh0_798, lsg_435, \
                         lsg_450, lsg_570, lsh1_798, msg_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = pa_x[k] * lsh0_798[k]
                   + f_17 * lsg_570[k]
                   - f_8 * pc_x[k] * lsh1_798[k];

        t_799[k] = f_16 * lsg_450[k]
                   + f_3 * pc_y[k] * msg_570[k];

        t_800[k] = f_10 * lsg_435[k]
                   + f_3 * pc_z[k] * msg_570[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, pa_x, pc_x, pc_y, lsh0_801, lsh0_803, lsg_452, \
                         lsg_573, lsg_575, lsh1_801, lsh1_803, \
                         msg_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = pa_x[k] * lsh0_801[k]
                   + f_11 * lsg_573[k]
                   - f_8 * pc_x[k] * lsh1_801[k];

        t_802[k] = f_16 * lsg_452[k]
                   + f_3 * pc_y[k] * msg_572[k];

        t_803[k] = pa_x[k] * lsh0_803[k]
                   + f_11 * lsg_575[k]
                   - f_8 * pc_x[k] * lsh1_803[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pa_x, pc_x, pc_y, pc_z, lsh0_804, lsg_438, \
                         lsg_455, lsg_576, lsh1_804, msg_573, msg_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = pa_x[k] * lsh0_804[k]
                   + f_10 * lsg_576[k]
                   - f_8 * pc_x[k] * lsh1_804[k];

        t_805[k] = f_10 * lsg_438[k]
                   + f_3 * pc_z[k] * msg_573[k];

        t_806[k] = f_16 * lsg_455[k]
                   + f_3 * pc_y[k] * msg_575[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, pa_x, pc_x, lsh0_807, lsg_579, lsg_580, \
                         lsg_581, lsg_582, lsh1_807, msg_580, msg_581, \
                         msg_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = pa_x[k] * lsh0_807[k]
                   + f_10 * lsg_579[k]
                   - f_8 * pc_x[k] * lsh1_807[k];

        t_808[k] = f_9 * lsg_580[k]
                   + f_3 * pc_x[k] * msg_580[k];

        t_809[k] = f_9 * lsg_581[k]
                   + f_3 * pc_x[k] * msg_581[k];

        t_810[k] = f_9 * lsg_582[k]
                   + f_3 * pc_x[k] * msg_582[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, pa_x, pc_x, pc_z, lsh0_813, lsg_445, \
                         lsg_583, lsg_584, lsh1_813, msg_580, msg_583, \
                         msg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_9 * lsg_583[k]
                   + f_3 * pc_x[k] * msg_583[k];

        t_812[k] = f_9 * lsg_584[k]
                   + f_3 * pc_x[k] * msg_584[k];

        t_813[k] = pa_x[k] * lsh0_813[k]
                   - f_8 * pc_x[k] * lsh1_813[k];

        t_814[k] = f_10 * lsg_445[k]
                   + f_3 * pc_z[k] * msg_580[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, pa_x, pc_x, pc_y, lsh0_815, lsh0_816, \
                         lsh0_818, lsg_464, lsh1_815, lsh1_816, lsh1_818, \
                         msg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = pa_x[k] * lsh0_815[k]
                   - f_8 * pc_x[k] * lsh1_815[k];

        t_816[k] = pa_x[k] * lsh0_816[k]
                   - f_8 * pc_x[k] * lsh1_816[k];

        t_817[k] = f_16 * lsg_464[k]
                   + f_3 * pc_y[k] * msg_584[k];

        t_818[k] = pa_x[k] * lsh0_818[k]
                   - f_8 * pc_x[k] * lsh1_818[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pa_x, pc_x, pc_y, pc_z, lsh0_819, lsg_450, \
                         lsg_465, lsg_585, lsh1_819, msg_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = pa_x[k] * lsh0_819[k]
                   + f_17 * lsg_585[k]
                   - f_8 * pc_x[k] * lsh1_819[k];

        t_820[k] = f_17 * lsg_465[k]
                   + f_3 * pc_y[k] * msg_585[k];

        t_821[k] = f_11 * lsg_450[k]
                   + f_3 * pc_z[k] * msg_585[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pa_x, pc_x, pc_y, lsh0_822, lsh0_824, lsg_467, \
                         lsg_588, lsg_590, lsh1_822, lsh1_824, \
                         msg_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pa_x[k] * lsh0_822[k]
                   + f_11 * lsg_588[k]
                   - f_8 * pc_x[k] * lsh1_822[k];

        t_823[k] = f_17 * lsg_467[k]
                   + f_3 * pc_y[k] * msg_587[k];

        t_824[k] = pa_x[k] * lsh0_824[k]
                   + f_11 * lsg_590[k]
                   - f_8 * pc_x[k] * lsh1_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, pa_x, pc_x, pc_y, pc_z, lsh0_825, lsg_453, \
                         lsg_470, lsg_591, lsh1_825, msg_588, msg_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = pa_x[k] * lsh0_825[k]
                   + f_10 * lsg_591[k]
                   - f_8 * pc_x[k] * lsh1_825[k];

        t_826[k] = f_11 * lsg_453[k]
                   + f_3 * pc_z[k] * msg_588[k];

        t_827[k] = f_17 * lsg_470[k]
                   + f_3 * pc_y[k] * msg_590[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, t_831, pa_x, pc_x, lsh0_828, lsg_594, lsg_595, \
                         lsg_596, lsg_597, lsh1_828, msg_595, msg_596, \
                         msg_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = pa_x[k] * lsh0_828[k]
                   + f_10 * lsg_594[k]
                   - f_8 * pc_x[k] * lsh1_828[k];

        t_829[k] = f_9 * lsg_595[k]
                   + f_3 * pc_x[k] * msg_595[k];

        t_830[k] = f_9 * lsg_596[k]
                   + f_3 * pc_x[k] * msg_596[k];

        t_831[k] = f_9 * lsg_597[k]
                   + f_3 * pc_x[k] * msg_597[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, pa_x, pc_x, pc_z, lsh0_834, lsg_460, \
                         lsg_598, lsg_599, lsh1_834, msg_595, msg_598, \
                         msg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_9 * lsg_598[k]
                   + f_3 * pc_x[k] * msg_598[k];

        t_833[k] = f_9 * lsg_599[k]
                   + f_3 * pc_x[k] * msg_599[k];

        t_834[k] = pa_x[k] * lsh0_834[k]
                   - f_8 * pc_x[k] * lsh1_834[k];

        t_835[k] = f_11 * lsg_460[k]
                   + f_3 * pc_z[k] * msg_595[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_12 = 4.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_735 = buffer.data(lsh0 + 735);
    const auto *lsh0_740 = buffer.data(lsh0 + 740);
    const auto *lsh0_744 = buffer.data(lsh0 + 744);
    const auto *lsh0_836 = buffer.data(lsh0 + 836);
    const auto *lsh0_837 = buffer.data(lsh0 + 837);
    const auto *lsh0_839 = buffer.data(lsh0 + 839);
    const auto *lsh0_840 = buffer.data(lsh0 + 840);
    const auto *lsh0_843 = buffer.data(lsh0 + 843);
    const auto *lsh0_845 = buffer.data(lsh0 + 845);
    const auto *lsh0_846 = buffer.data(lsh0 + 846);
    const auto *lsh0_849 = buffer.data(lsh0 + 849);
    const auto *lsh0_855 = buffer.data(lsh0 + 855);
    const auto *lsh0_857 = buffer.data(lsh0 + 857);
    const auto *lsh0_858 = buffer.data(lsh0 + 858);
    const auto *lsh0_860 = buffer.data(lsh0 + 860);
    const auto *lsh0_861 = buffer.data(lsh0 + 861);
    const auto *lsh0_864 = buffer.data(lsh0 + 864);
    const auto *lsh0_866 = buffer.data(lsh0 + 866);
    const auto *lsh0_867 = buffer.data(lsh0 + 867);
    const auto *lsh0_870 = buffer.data(lsh0 + 870);
    const auto *lsh0_876 = buffer.data(lsh0 + 876);
    const auto *lsh0_878 = buffer.data(lsh0 + 878);
    const auto *lsh0_879 = buffer.data(lsh0 + 879);
    const auto *lsh0_881 = buffer.data(lsh0 + 881);
    const auto *lsh0_882 = buffer.data(lsh0 + 882);
    const auto *lsh0_885 = buffer.data(lsh0 + 885);
    const auto *lsh0_887 = buffer.data(lsh0 + 887);
    const auto *lsh0_888 = buffer.data(lsh0 + 888);
    const auto *lsh0_891 = buffer.data(lsh0 + 891);
    const auto *lsh0_897 = buffer.data(lsh0 + 897);
    const auto *lsh0_899 = buffer.data(lsh0 + 899);
    const auto *lsh0_900 = buffer.data(lsh0 + 900);
    const auto *lsh0_902 = buffer.data(lsh0 + 902);
    const auto *lsh0_906 = buffer.data(lsh0 + 906);
    const auto *lsh0_909 = buffer.data(lsh0 + 909);
    const auto *lsh0_918 = buffer.data(lsh0 + 918);
    const auto *lsh0_920 = buffer.data(lsh0 + 920);
    const auto *lsh0_921 = buffer.data(lsh0 + 921);
    const auto *lsh0_923 = buffer.data(lsh0 + 923);
    const auto *lsh0_924 = buffer.data(lsh0 + 924);
    const auto *lsh0_929 = buffer.data(lsh0 + 929);
    const auto *lsh0_933 = buffer.data(lsh0 + 933);
    const auto *lsh0_939 = buffer.data(lsh0 + 939);
    const auto *lsh0_940 = buffer.data(lsh0 + 940);
    const auto *lsh0_941 = buffer.data(lsh0 + 941);
    const auto *lsh0_942 = buffer.data(lsh0 + 942);
    const auto *lsh0_944 = buffer.data(lsh0 + 944);

    const auto *lsg_465 = buffer.data(lsg + 465);
    const auto *lsg_468 = buffer.data(lsg + 468);
    const auto *lsg_475 = buffer.data(lsg + 475);
    const auto *lsg_479 = buffer.data(lsg + 479);
    const auto *lsg_480 = buffer.data(lsg + 480);
    const auto *lsg_482 = buffer.data(lsg + 482);
    const auto *lsg_483 = buffer.data(lsg + 483);
    const auto *lsg_485 = buffer.data(lsg + 485);
    const auto *lsg_490 = buffer.data(lsg + 490);
    const auto *lsg_494 = buffer.data(lsg + 494);
    const auto *lsg_495 = buffer.data(lsg + 495);
    const auto *lsg_497 = buffer.data(lsg + 497);
    const auto *lsg_498 = buffer.data(lsg + 498);
    const auto *lsg_500 = buffer.data(lsg + 500);
    const auto *lsg_505 = buffer.data(lsg + 505);
    const auto *lsg_509 = buffer.data(lsg + 509);
    const auto *lsg_510 = buffer.data(lsg + 510);
    const auto *lsg_512 = buffer.data(lsg + 512);
    const auto *lsg_513 = buffer.data(lsg + 513);
    const auto *lsg_515 = buffer.data(lsg + 515);
    const auto *lsg_520 = buffer.data(lsg + 520);
    const auto *lsg_524 = buffer.data(lsg + 524);
    const auto *lsg_525 = buffer.data(lsg + 525);
    const auto *lsg_527 = buffer.data(lsg + 527);
    const auto *lsg_530 = buffer.data(lsg + 530);
    const auto *lsg_539 = buffer.data(lsg + 539);
    const auto *lsg_550 = buffer.data(lsg + 550);
    const auto *lsg_600 = buffer.data(lsg + 600);
    const auto *lsg_603 = buffer.data(lsg + 603);
    const auto *lsg_605 = buffer.data(lsg + 605);
    const auto *lsg_606 = buffer.data(lsg + 606);
    const auto *lsg_609 = buffer.data(lsg + 609);
    const auto *lsg_610 = buffer.data(lsg + 610);
    const auto *lsg_611 = buffer.data(lsg + 611);
    const auto *lsg_612 = buffer.data(lsg + 612);
    const auto *lsg_613 = buffer.data(lsg + 613);
    const auto *lsg_614 = buffer.data(lsg + 614);
    const auto *lsg_615 = buffer.data(lsg + 615);
    const auto *lsg_618 = buffer.data(lsg + 618);
    const auto *lsg_620 = buffer.data(lsg + 620);
    const auto *lsg_621 = buffer.data(lsg + 621);
    const auto *lsg_624 = buffer.data(lsg + 624);
    const auto *lsg_625 = buffer.data(lsg + 625);
    const auto *lsg_626 = buffer.data(lsg + 626);
    const auto *lsg_627 = buffer.data(lsg + 627);
    const auto *lsg_628 = buffer.data(lsg + 628);
    const auto *lsg_629 = buffer.data(lsg + 629);
    const auto *lsg_630 = buffer.data(lsg + 630);
    const auto *lsg_633 = buffer.data(lsg + 633);
    const auto *lsg_635 = buffer.data(lsg + 635);
    const auto *lsg_636 = buffer.data(lsg + 636);
    const auto *lsg_639 = buffer.data(lsg + 639);
    const auto *lsg_640 = buffer.data(lsg + 640);
    const auto *lsg_641 = buffer.data(lsg + 641);
    const auto *lsg_642 = buffer.data(lsg + 642);
    const auto *lsg_643 = buffer.data(lsg + 643);
    const auto *lsg_644 = buffer.data(lsg + 644);
    const auto *lsg_648 = buffer.data(lsg + 648);
    const auto *lsg_651 = buffer.data(lsg + 651);
    const auto *lsg_655 = buffer.data(lsg + 655);
    const auto *lsg_656 = buffer.data(lsg + 656);
    const auto *lsg_657 = buffer.data(lsg + 657);
    const auto *lsg_658 = buffer.data(lsg + 658);
    const auto *lsg_659 = buffer.data(lsg + 659);
    const auto *lsg_660 = buffer.data(lsg + 660);
    const auto *lsg_665 = buffer.data(lsg + 665);
    const auto *lsg_669 = buffer.data(lsg + 669);
    const auto *lsg_670 = buffer.data(lsg + 670);
    const auto *lsg_671 = buffer.data(lsg + 671);
    const auto *lsg_672 = buffer.data(lsg + 672);
    const auto *lsg_674 = buffer.data(lsg + 674);

    const auto *lsh1_735 = buffer.data(lsh1 + 735);
    const auto *lsh1_740 = buffer.data(lsh1 + 740);
    const auto *lsh1_744 = buffer.data(lsh1 + 744);
    const auto *lsh1_836 = buffer.data(lsh1 + 836);
    const auto *lsh1_837 = buffer.data(lsh1 + 837);
    const auto *lsh1_839 = buffer.data(lsh1 + 839);
    const auto *lsh1_840 = buffer.data(lsh1 + 840);
    const auto *lsh1_843 = buffer.data(lsh1 + 843);
    const auto *lsh1_845 = buffer.data(lsh1 + 845);
    const auto *lsh1_846 = buffer.data(lsh1 + 846);
    const auto *lsh1_849 = buffer.data(lsh1 + 849);
    const auto *lsh1_855 = buffer.data(lsh1 + 855);
    const auto *lsh1_857 = buffer.data(lsh1 + 857);
    const auto *lsh1_858 = buffer.data(lsh1 + 858);
    const auto *lsh1_860 = buffer.data(lsh1 + 860);
    const auto *lsh1_861 = buffer.data(lsh1 + 861);
    const auto *lsh1_864 = buffer.data(lsh1 + 864);
    const auto *lsh1_866 = buffer.data(lsh1 + 866);
    const auto *lsh1_867 = buffer.data(lsh1 + 867);
    const auto *lsh1_870 = buffer.data(lsh1 + 870);
    const auto *lsh1_876 = buffer.data(lsh1 + 876);
    const auto *lsh1_878 = buffer.data(lsh1 + 878);
    const auto *lsh1_879 = buffer.data(lsh1 + 879);
    const auto *lsh1_881 = buffer.data(lsh1 + 881);
    const auto *lsh1_882 = buffer.data(lsh1 + 882);
    const auto *lsh1_885 = buffer.data(lsh1 + 885);
    const auto *lsh1_887 = buffer.data(lsh1 + 887);
    const auto *lsh1_888 = buffer.data(lsh1 + 888);
    const auto *lsh1_891 = buffer.data(lsh1 + 891);
    const auto *lsh1_897 = buffer.data(lsh1 + 897);
    const auto *lsh1_899 = buffer.data(lsh1 + 899);
    const auto *lsh1_900 = buffer.data(lsh1 + 900);
    const auto *lsh1_902 = buffer.data(lsh1 + 902);
    const auto *lsh1_906 = buffer.data(lsh1 + 906);
    const auto *lsh1_909 = buffer.data(lsh1 + 909);
    const auto *lsh1_918 = buffer.data(lsh1 + 918);
    const auto *lsh1_920 = buffer.data(lsh1 + 920);
    const auto *lsh1_921 = buffer.data(lsh1 + 921);
    const auto *lsh1_923 = buffer.data(lsh1 + 923);
    const auto *lsh1_924 = buffer.data(lsh1 + 924);
    const auto *lsh1_929 = buffer.data(lsh1 + 929);
    const auto *lsh1_933 = buffer.data(lsh1 + 933);
    const auto *lsh1_939 = buffer.data(lsh1 + 939);
    const auto *lsh1_940 = buffer.data(lsh1 + 940);
    const auto *lsh1_941 = buffer.data(lsh1 + 941);
    const auto *lsh1_942 = buffer.data(lsh1 + 942);
    const auto *lsh1_944 = buffer.data(lsh1 + 944);

    const auto *msf0_440 = buffer.data(msf0 + 440);
    const auto *msf0_441 = buffer.data(msf0 + 441);
    const auto *msf0_442 = buffer.data(msf0 + 442);
    const auto *msf0_450 = buffer.data(msf0 + 450);
    const auto *msf0_451 = buffer.data(msf0 + 451);
    const auto *msf0_453 = buffer.data(msf0 + 453);
    const auto *msf0_455 = buffer.data(msf0 + 455);
    const auto *msf0_456 = buffer.data(msf0 + 456);
    const auto *msf0_458 = buffer.data(msf0 + 458);
    const auto *msf0_459 = buffer.data(msf0 + 459);

    const auto *msf1_440 = buffer.data(msf1 + 440);
    const auto *msf1_441 = buffer.data(msf1 + 441);
    const auto *msf1_442 = buffer.data(msf1 + 442);
    const auto *msf1_450 = buffer.data(msf1 + 450);
    const auto *msf1_451 = buffer.data(msf1 + 451);
    const auto *msf1_453 = buffer.data(msf1 + 453);
    const auto *msf1_455 = buffer.data(msf1 + 455);
    const auto *msf1_456 = buffer.data(msf1 + 456);
    const auto *msf1_458 = buffer.data(msf1 + 458);
    const auto *msf1_459 = buffer.data(msf1 + 459);

    const auto *msg_599 = buffer.data(msg + 599);
    const auto *msg_600 = buffer.data(msg + 600);
    const auto *msg_602 = buffer.data(msg + 602);
    const auto *msg_603 = buffer.data(msg + 603);
    const auto *msg_605 = buffer.data(msg + 605);
    const auto *msg_610 = buffer.data(msg + 610);
    const auto *msg_611 = buffer.data(msg + 611);
    const auto *msg_612 = buffer.data(msg + 612);
    const auto *msg_613 = buffer.data(msg + 613);
    const auto *msg_614 = buffer.data(msg + 614);
    const auto *msg_615 = buffer.data(msg + 615);
    const auto *msg_617 = buffer.data(msg + 617);
    const auto *msg_618 = buffer.data(msg + 618);
    const auto *msg_620 = buffer.data(msg + 620);
    const auto *msg_625 = buffer.data(msg + 625);
    const auto *msg_626 = buffer.data(msg + 626);
    const auto *msg_627 = buffer.data(msg + 627);
    const auto *msg_628 = buffer.data(msg + 628);
    const auto *msg_629 = buffer.data(msg + 629);
    const auto *msg_630 = buffer.data(msg + 630);
    const auto *msg_632 = buffer.data(msg + 632);
    const auto *msg_633 = buffer.data(msg + 633);
    const auto *msg_635 = buffer.data(msg + 635);
    const auto *msg_640 = buffer.data(msg + 640);
    const auto *msg_641 = buffer.data(msg + 641);
    const auto *msg_642 = buffer.data(msg + 642);
    const auto *msg_643 = buffer.data(msg + 643);
    const auto *msg_644 = buffer.data(msg + 644);
    const auto *msg_645 = buffer.data(msg + 645);
    const auto *msg_647 = buffer.data(msg + 647);
    const auto *msg_648 = buffer.data(msg + 648);
    const auto *msg_650 = buffer.data(msg + 650);
    const auto *msg_655 = buffer.data(msg + 655);
    const auto *msg_656 = buffer.data(msg + 656);
    const auto *msg_657 = buffer.data(msg + 657);
    const auto *msg_658 = buffer.data(msg + 658);
    const auto *msg_659 = buffer.data(msg + 659);
    const auto *msg_660 = buffer.data(msg + 660);
    const auto *msg_661 = buffer.data(msg + 661);
    const auto *msg_662 = buffer.data(msg + 662);
    const auto *msg_663 = buffer.data(msg + 663);
    const auto *msg_664 = buffer.data(msg + 664);
    const auto *msg_665 = buffer.data(msg + 665);
    const auto *msg_669 = buffer.data(msg + 669);
    const auto *msg_670 = buffer.data(msg + 670);
    const auto *msg_671 = buffer.data(msg + 671);
    const auto *msg_672 = buffer.data(msg + 672);
    const auto *msg_674 = buffer.data(msg + 674);
    const auto *msg_675 = buffer.data(msg + 675);
    const auto *msg_676 = buffer.data(msg + 676);
    const auto *msg_678 = buffer.data(msg + 678);
    const auto *msg_680 = buffer.data(msg + 680);
    const auto *msg_681 = buffer.data(msg + 681);
    const auto *msg_683 = buffer.data(msg + 683);
    const auto *msg_684 = buffer.data(msg + 684);
    const auto *msg_685 = buffer.data(msg + 685);
    const auto *msg_686 = buffer.data(msg + 686);
    const auto *msg_687 = buffer.data(msg + 687);
    const auto *msg_688 = buffer.data(msg + 688);
    const auto *msg_689 = buffer.data(msg + 689);

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pa_x, pc_x, pc_y, lsh0_836, lsh0_837, \
                         lsh0_839, lsg_479, lsh1_836, lsh1_837, lsh1_839, \
                         msg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = pa_x[k] * lsh0_836[k]
                   - f_8 * pc_x[k] * lsh1_836[k];

        t_837[k] = pa_x[k] * lsh0_837[k]
                   - f_8 * pc_x[k] * lsh1_837[k];

        t_838[k] = f_17 * lsg_479[k]
                   + f_3 * pc_y[k] * msg_599[k];

        t_839[k] = pa_x[k] * lsh0_839[k]
                   - f_8 * pc_x[k] * lsh1_839[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pa_x, pc_x, pc_y, pc_z, lsh0_840, lsg_465, \
                         lsg_480, lsg_600, lsh1_840, msg_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = pa_x[k] * lsh0_840[k]
                   + f_17 * lsg_600[k]
                   - f_8 * pc_x[k] * lsh1_840[k];

        t_841[k] = f_18 * lsg_480[k]
                   + f_3 * pc_y[k] * msg_600[k];

        t_842[k] = f_18 * lsg_465[k]
                   + f_3 * pc_z[k] * msg_600[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pa_x, pc_x, pc_y, lsh0_843, lsh0_845, lsg_482, \
                         lsg_603, lsg_605, lsh1_843, lsh1_845, \
                         msg_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = pa_x[k] * lsh0_843[k]
                   + f_11 * lsg_603[k]
                   - f_8 * pc_x[k] * lsh1_843[k];

        t_844[k] = f_18 * lsg_482[k]
                   + f_3 * pc_y[k] * msg_602[k];

        t_845[k] = pa_x[k] * lsh0_845[k]
                   + f_11 * lsg_605[k]
                   - f_8 * pc_x[k] * lsh1_845[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pa_x, pc_x, pc_y, pc_z, lsh0_846, lsg_468, \
                         lsg_485, lsg_606, lsh1_846, msg_603, msg_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = pa_x[k] * lsh0_846[k]
                   + f_10 * lsg_606[k]
                   - f_8 * pc_x[k] * lsh1_846[k];

        t_847[k] = f_18 * lsg_468[k]
                   + f_3 * pc_z[k] * msg_603[k];

        t_848[k] = f_18 * lsg_485[k]
                   + f_3 * pc_y[k] * msg_605[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, pa_x, pc_x, lsh0_849, lsg_609, lsg_610, \
                         lsg_611, lsg_612, lsh1_849, msg_610, msg_611, \
                         msg_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = pa_x[k] * lsh0_849[k]
                   + f_10 * lsg_609[k]
                   - f_8 * pc_x[k] * lsh1_849[k];

        t_850[k] = f_9 * lsg_610[k]
                   + f_3 * pc_x[k] * msg_610[k];

        t_851[k] = f_9 * lsg_611[k]
                   + f_3 * pc_x[k] * msg_611[k];

        t_852[k] = f_9 * lsg_612[k]
                   + f_3 * pc_x[k] * msg_612[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, pa_x, pc_x, pc_z, lsh0_855, lsg_475, \
                         lsg_613, lsg_614, lsh1_855, msg_610, msg_613, \
                         msg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_9 * lsg_613[k]
                   + f_3 * pc_x[k] * msg_613[k];

        t_854[k] = f_9 * lsg_614[k]
                   + f_3 * pc_x[k] * msg_614[k];

        t_855[k] = pa_x[k] * lsh0_855[k]
                   - f_8 * pc_x[k] * lsh1_855[k];

        t_856[k] = f_18 * lsg_475[k]
                   + f_3 * pc_z[k] * msg_610[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, t_860, pa_x, pc_x, pc_y, lsh0_857, lsh0_858, \
                         lsh0_860, lsg_494, lsh1_857, lsh1_858, lsh1_860, \
                         msg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = pa_x[k] * lsh0_857[k]
                   - f_8 * pc_x[k] * lsh1_857[k];

        t_858[k] = pa_x[k] * lsh0_858[k]
                   - f_8 * pc_x[k] * lsh1_858[k];

        t_859[k] = f_18 * lsg_494[k]
                   + f_3 * pc_y[k] * msg_614[k];

        t_860[k] = pa_x[k] * lsh0_860[k]
                   - f_8 * pc_x[k] * lsh1_860[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, pa_x, pc_x, pc_y, pc_z, lsh0_861, lsg_480, \
                         lsg_495, lsg_615, lsh1_861, msg_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = pa_x[k] * lsh0_861[k]
                   + f_17 * lsg_615[k]
                   - f_8 * pc_x[k] * lsh1_861[k];

        t_862[k] = f_11 * lsg_495[k]
                   + f_3 * pc_y[k] * msg_615[k];

        t_863[k] = f_17 * lsg_480[k]
                   + f_3 * pc_z[k] * msg_615[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, pa_x, pc_x, pc_y, lsh0_864, lsh0_866, lsg_497, \
                         lsg_618, lsg_620, lsh1_864, lsh1_866, \
                         msg_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = pa_x[k] * lsh0_864[k]
                   + f_11 * lsg_618[k]
                   - f_8 * pc_x[k] * lsh1_864[k];

        t_865[k] = f_11 * lsg_497[k]
                   + f_3 * pc_y[k] * msg_617[k];

        t_866[k] = pa_x[k] * lsh0_866[k]
                   + f_11 * lsg_620[k]
                   - f_8 * pc_x[k] * lsh1_866[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, pa_x, pc_x, pc_y, pc_z, lsh0_867, lsg_483, \
                         lsg_500, lsg_621, lsh1_867, msg_618, msg_620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = pa_x[k] * lsh0_867[k]
                   + f_10 * lsg_621[k]
                   - f_8 * pc_x[k] * lsh1_867[k];

        t_868[k] = f_17 * lsg_483[k]
                   + f_3 * pc_z[k] * msg_618[k];

        t_869[k] = f_11 * lsg_500[k]
                   + f_3 * pc_y[k] * msg_620[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, pa_x, pc_x, lsh0_870, lsg_624, lsg_625, \
                         lsg_626, lsg_627, lsh1_870, msg_625, msg_626, \
                         msg_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pa_x[k] * lsh0_870[k]
                   + f_10 * lsg_624[k]
                   - f_8 * pc_x[k] * lsh1_870[k];

        t_871[k] = f_9 * lsg_625[k]
                   + f_3 * pc_x[k] * msg_625[k];

        t_872[k] = f_9 * lsg_626[k]
                   + f_3 * pc_x[k] * msg_626[k];

        t_873[k] = f_9 * lsg_627[k]
                   + f_3 * pc_x[k] * msg_627[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, pa_x, pc_x, pc_z, lsh0_876, lsg_490, \
                         lsg_628, lsg_629, lsh1_876, msg_625, msg_628, \
                         msg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_9 * lsg_628[k]
                   + f_3 * pc_x[k] * msg_628[k];

        t_875[k] = f_9 * lsg_629[k]
                   + f_3 * pc_x[k] * msg_629[k];

        t_876[k] = pa_x[k] * lsh0_876[k]
                   - f_8 * pc_x[k] * lsh1_876[k];

        t_877[k] = f_17 * lsg_490[k]
                   + f_3 * pc_z[k] * msg_625[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, t_881, pa_x, pc_x, pc_y, lsh0_878, lsh0_879, \
                         lsh0_881, lsg_509, lsh1_878, lsh1_879, lsh1_881, \
                         msg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = pa_x[k] * lsh0_878[k]
                   - f_8 * pc_x[k] * lsh1_878[k];

        t_879[k] = pa_x[k] * lsh0_879[k]
                   - f_8 * pc_x[k] * lsh1_879[k];

        t_880[k] = f_11 * lsg_509[k]
                   + f_3 * pc_y[k] * msg_629[k];

        t_881[k] = pa_x[k] * lsh0_881[k]
                   - f_8 * pc_x[k] * lsh1_881[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, pa_x, pc_x, pc_y, pc_z, lsh0_882, lsg_495, \
                         lsg_510, lsg_630, lsh1_882, msg_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pa_x[k] * lsh0_882[k]
                   + f_17 * lsg_630[k]
                   - f_8 * pc_x[k] * lsh1_882[k];

        t_883[k] = f_10 * lsg_510[k]
                   + f_3 * pc_y[k] * msg_630[k];

        t_884[k] = f_16 * lsg_495[k]
                   + f_3 * pc_z[k] * msg_630[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, pa_x, pc_x, pc_y, lsh0_885, lsh0_887, lsg_512, \
                         lsg_633, lsg_635, lsh1_885, lsh1_887, \
                         msg_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = pa_x[k] * lsh0_885[k]
                   + f_11 * lsg_633[k]
                   - f_8 * pc_x[k] * lsh1_885[k];

        t_886[k] = f_10 * lsg_512[k]
                   + f_3 * pc_y[k] * msg_632[k];

        t_887[k] = pa_x[k] * lsh0_887[k]
                   + f_11 * lsg_635[k]
                   - f_8 * pc_x[k] * lsh1_887[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, pa_x, pc_x, pc_y, pc_z, lsh0_888, lsg_498, \
                         lsg_515, lsg_636, lsh1_888, msg_633, msg_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = pa_x[k] * lsh0_888[k]
                   + f_10 * lsg_636[k]
                   - f_8 * pc_x[k] * lsh1_888[k];

        t_889[k] = f_16 * lsg_498[k]
                   + f_3 * pc_z[k] * msg_633[k];

        t_890[k] = f_10 * lsg_515[k]
                   + f_3 * pc_y[k] * msg_635[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, pa_x, pc_x, lsh0_891, lsg_639, lsg_640, \
                         lsg_641, lsg_642, lsh1_891, msg_640, msg_641, \
                         msg_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = pa_x[k] * lsh0_891[k]
                   + f_10 * lsg_639[k]
                   - f_8 * pc_x[k] * lsh1_891[k];

        t_892[k] = f_9 * lsg_640[k]
                   + f_3 * pc_x[k] * msg_640[k];

        t_893[k] = f_9 * lsg_641[k]
                   + f_3 * pc_x[k] * msg_641[k];

        t_894[k] = f_9 * lsg_642[k]
                   + f_3 * pc_x[k] * msg_642[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, pa_x, pc_x, pc_z, lsh0_897, lsg_505, \
                         lsg_643, lsg_644, lsh1_897, msg_640, msg_643, \
                         msg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_9 * lsg_643[k]
                   + f_3 * pc_x[k] * msg_643[k];

        t_896[k] = f_9 * lsg_644[k]
                   + f_3 * pc_x[k] * msg_644[k];

        t_897[k] = pa_x[k] * lsh0_897[k]
                   - f_8 * pc_x[k] * lsh1_897[k];

        t_898[k] = f_16 * lsg_505[k]
                   + f_3 * pc_z[k] * msg_640[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pa_x, pc_x, pc_y, lsh0_899, lsh0_900, \
                         lsh0_902, lsg_524, lsh1_899, lsh1_900, lsh1_902, \
                         msg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pa_x[k] * lsh0_899[k]
                   - f_8 * pc_x[k] * lsh1_899[k];

        t_900[k] = pa_x[k] * lsh0_900[k]
                   - f_8 * pc_x[k] * lsh1_900[k];

        t_901[k] = f_10 * lsg_524[k]
                   + f_3 * pc_y[k] * msg_644[k];

        t_902[k] = pa_x[k] * lsh0_902[k]
                   - f_8 * pc_x[k] * lsh1_902[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pa_y, pc_y, pc_z, lsh0_735, lsg_510, lsg_525, \
                         lsh1_735, msg_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = pa_y[k] * lsh0_735[k]
                   - f_8 * pc_y[k] * lsh1_735[k];

        t_904[k] = f_9 * lsg_525[k]
                   + f_3 * pc_y[k] * msg_645[k];

        t_905[k] = f_15 * lsg_510[k]
                   + f_3 * pc_z[k] * msg_645[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pa_x, pa_y, pc_x, pc_y, lsh0_740, lsh0_906, \
                         lsg_527, lsg_648, lsh1_740, lsh1_906, \
                         msg_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = pa_x[k] * lsh0_906[k]
                   + f_11 * lsg_648[k]
                   - f_8 * pc_x[k] * lsh1_906[k];

        t_907[k] = f_9 * lsg_527[k]
                   + f_3 * pc_y[k] * msg_647[k];

        t_908[k] = pa_y[k] * lsh0_740[k]
                   - f_8 * pc_y[k] * lsh1_740[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pa_x, pc_x, pc_y, pc_z, lsh0_909, lsg_513, \
                         lsg_530, lsg_651, lsh1_909, msg_648, msg_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pa_x[k] * lsh0_909[k]
                   + f_10 * lsg_651[k]
                   - f_8 * pc_x[k] * lsh1_909[k];

        t_910[k] = f_15 * lsg_513[k]
                   + f_3 * pc_z[k] * msg_648[k];

        t_911[k] = f_9 * lsg_530[k]
                   + f_3 * pc_y[k] * msg_650[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, t_915, pa_y, pc_x, pc_y, lsh0_744, lsg_655, \
                         lsg_656, lsg_657, lsh1_744, msg_655, msg_656, \
                         msg_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = pa_y[k] * lsh0_744[k]
                   - f_8 * pc_y[k] * lsh1_744[k];

        t_913[k] = f_9 * lsg_655[k]
                   + f_3 * pc_x[k] * msg_655[k];

        t_914[k] = f_9 * lsg_656[k]
                   + f_3 * pc_x[k] * msg_656[k];

        t_915[k] = f_9 * lsg_657[k]
                   + f_3 * pc_x[k] * msg_657[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, t_919, pa_x, pc_x, pc_z, lsh0_918, lsg_520, \
                         lsg_658, lsg_659, lsh1_918, msg_655, msg_658, \
                         msg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_9 * lsg_658[k]
                   + f_3 * pc_x[k] * msg_658[k];

        t_917[k] = f_9 * lsg_659[k]
                   + f_3 * pc_x[k] * msg_659[k];

        t_918[k] = pa_x[k] * lsh0_918[k]
                   - f_8 * pc_x[k] * lsh1_918[k];

        t_919[k] = f_15 * lsg_520[k]
                   + f_3 * pc_z[k] * msg_655[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pa_x, pc_x, pc_y, lsh0_920, lsh0_921, \
                         lsh0_923, lsg_539, lsh1_920, lsh1_921, lsh1_923, \
                         msg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = pa_x[k] * lsh0_920[k]
                   - f_8 * pc_x[k] * lsh1_920[k];

        t_921[k] = pa_x[k] * lsh0_921[k]
                   - f_8 * pc_x[k] * lsh1_921[k];

        t_922[k] = f_9 * lsg_539[k]
                   + f_3 * pc_y[k] * msg_659[k];

        t_923[k] = pa_x[k] * lsh0_923[k]
                   - f_8 * pc_x[k] * lsh1_923[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pa_x, pc_x, pc_y, pc_z, lsh0_924, \
                         lsg_525, lsg_660, lsh1_924, msf0_440, msf1_440, msg_660, \
                         msg_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = pa_x[k] * lsh0_924[k]
                   + f_17 * lsg_660[k]
                   - f_8 * pc_x[k] * lsh1_924[k];

        t_925[k] = f_3 * pc_y[k] * msg_660[k];

        t_926[k] = f_12 * lsg_525[k]
                   + f_3 * pc_z[k] * msg_660[k];

        t_927[k] = f_4 * msf0_440[k]
                   - f_5 * msf1_440[k]
                   + f_3 * pc_y[k] * msg_661[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pa_x, pc_x, pc_y, lsh0_929, lsg_665, lsh1_929, \
                         msf0_441, msf1_441, msg_662, msg_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_3 * pc_y[k] * msg_662[k];

        t_929[k] = pa_x[k] * lsh0_929[k]
                   + f_11 * lsg_665[k]
                   - f_8 * pc_x[k] * lsh1_929[k];

        t_930[k] = f_6 * msf0_441[k]
                   - f_7 * msf1_441[k]
                   + f_3 * pc_y[k] * msg_663[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, pa_x, pc_x, pc_y, lsh0_933, lsg_669, \
                         lsg_670, lsh1_933, msf0_442, msf1_442, msg_664, msg_665, \
                         msg_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_4 * msf0_442[k]
                   - f_5 * msf1_442[k]
                   + f_3 * pc_y[k] * msg_664[k];

        t_932[k] = f_3 * pc_y[k] * msg_665[k];

        t_933[k] = pa_x[k] * lsh0_933[k]
                   + f_10 * lsg_669[k]
                   - f_8 * pc_x[k] * lsh1_933[k];

        t_934[k] = f_9 * lsg_670[k]
                   + f_3 * pc_x[k] * msg_670[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, pc_x, pc_y, lsg_671, lsg_672, lsg_674, \
                         msg_669, msg_671, msg_672, msg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_9 * lsg_671[k]
                   + f_3 * pc_x[k] * msg_671[k];

        t_936[k] = f_9 * lsg_672[k]
                   + f_3 * pc_x[k] * msg_672[k];

        t_937[k] = f_3 * pc_y[k] * msg_669[k];

        t_938[k] = f_9 * lsg_674[k]
                   + f_3 * pc_x[k] * msg_674[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pa_x, pc_x, lsh0_939, lsh0_940, lsh0_941, \
                         lsh0_942, lsh1_939, lsh1_940, lsh1_941, \
                         lsh1_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = pa_x[k] * lsh0_939[k]
                   - f_8 * pc_x[k] * lsh1_939[k];

        t_940[k] = pa_x[k] * lsh0_940[k]
                   - f_8 * pc_x[k] * lsh1_940[k];

        t_941[k] = pa_x[k] * lsh0_941[k]
                   - f_8 * pc_x[k] * lsh1_941[k];

        t_942[k] = pa_x[k] * lsh0_942[k]
                   - f_8 * pc_x[k] * lsh1_942[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, pa_x, pc_x, pc_y, lsh0_944, lsh1_944, \
                         msf0_450, msf0_451, msf1_450, msf1_451, msg_674, msg_675, \
                         msg_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_3 * pc_y[k] * msg_674[k];

        t_944[k] = pa_x[k] * lsh0_944[k]
                   - f_8 * pc_x[k] * lsh1_944[k];

        t_945[k] = f_1 * msf0_450[k]
                   - f_2 * msf1_450[k]
                   + f_3 * pc_x[k] * msg_675[k];

        t_946[k] = f_13 * msf0_451[k]
                   - f_14 * msf1_451[k]
                   + f_3 * pc_x[k] * msg_676[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, t_950, pc_x, pc_z, msf0_453, msf0_455, msf1_453, \
                         msf1_455, msg_675, msg_676, msg_678, msg_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_3 * pc_z[k] * msg_675[k];

        t_948[k] = f_6 * msf0_453[k]
                   - f_7 * msf1_453[k]
                   + f_3 * pc_x[k] * msg_678[k];

        t_949[k] = f_3 * pc_z[k] * msg_676[k];

        t_950[k] = f_6 * msf0_455[k]
                   - f_7 * msf1_455[k]
                   + f_3 * pc_x[k] * msg_680[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, t_954, pc_x, pc_z, msf0_456, msf0_458, msf0_459, \
                         msf1_456, msf1_458, msf1_459, msg_678, msg_681, msg_683, \
                         msg_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_4 * msf0_456[k]
                   - f_5 * msf1_456[k]
                   + f_3 * pc_x[k] * msg_681[k];

        t_952[k] = f_3 * pc_z[k] * msg_678[k];

        t_953[k] = f_4 * msf0_458[k]
                   - f_5 * msf1_458[k]
                   + f_3 * pc_x[k] * msg_683[k];

        t_954[k] = f_4 * msf0_459[k]
                   - f_5 * msf1_459[k]
                   + f_3 * pc_x[k] * msg_684[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, t_960, pc_x, pc_y, lsg_550, \
                         msf0_456, msf1_456, msg_685, msg_686, msg_687, msg_688, \
                         msg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_3 * pc_x[k] * msg_685[k];

        t_956[k] = f_3 * pc_x[k] * msg_686[k];

        t_957[k] = f_3 * pc_x[k] * msg_687[k];

        t_958[k] = f_3 * pc_x[k] * msg_688[k];

        t_959[k] = f_3 * pc_x[k] * msg_689[k];

        t_960[k] = f_0 * lsg_550[k]
                   + f_1 * msf0_456[k]
                   - f_2 * msf1_456[k]
                   + f_3 * pc_y[k] * msg_685[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_12 = 4.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;

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
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_756 = buffer.data(lsh0 + 756);
    const auto *lsh0_757 = buffer.data(lsh0 + 757);
    const auto *lsh0_759 = buffer.data(lsh0 + 759);
    const auto *lsh0_762 = buffer.data(lsh0 + 762);
    const auto *lsh0_771 = buffer.data(lsh0 + 771);
    const auto *lsh0_773 = buffer.data(lsh0 + 773);
    const auto *lsh0_774 = buffer.data(lsh0 + 774);

    const auto *lsg_550 = buffer.data(lsg + 550);
    const auto *lsg_551 = buffer.data(lsg + 551);
    const auto *lsg_552 = buffer.data(lsg + 552);
    const auto *lsg_554 = buffer.data(lsg + 554);
    const auto *lsg_565 = buffer.data(lsg + 565);
    const auto *lsg_569 = buffer.data(lsg + 569);
    const auto *lsg_580 = buffer.data(lsg + 580);
    const auto *lsg_582 = buffer.data(lsg + 582);
    const auto *lsg_583 = buffer.data(lsg + 583);
    const auto *lsg_584 = buffer.data(lsg + 584);
    const auto *lsg_595 = buffer.data(lsg + 595);
    const auto *lsg_597 = buffer.data(lsg + 597);
    const auto *lsg_598 = buffer.data(lsg + 598);
    const auto *lsg_599 = buffer.data(lsg + 599);
    const auto *lsg_610 = buffer.data(lsg + 610);
    const auto *lsg_612 = buffer.data(lsg + 612);
    const auto *lsg_613 = buffer.data(lsg + 613);
    const auto *lsg_614 = buffer.data(lsg + 614);
    const auto *lsg_625 = buffer.data(lsg + 625);
    const auto *lsg_627 = buffer.data(lsg + 627);
    const auto *lsg_628 = buffer.data(lsg + 628);
    const auto *lsg_629 = buffer.data(lsg + 629);

    const auto *lsh1_756 = buffer.data(lsh1 + 756);
    const auto *lsh1_757 = buffer.data(lsh1 + 757);
    const auto *lsh1_759 = buffer.data(lsh1 + 759);
    const auto *lsh1_762 = buffer.data(lsh1 + 762);
    const auto *lsh1_771 = buffer.data(lsh1 + 771);
    const auto *lsh1_773 = buffer.data(lsh1 + 773);
    const auto *lsh1_774 = buffer.data(lsh1 + 774);

    const auto *msf0_456 = buffer.data(msf0 + 456);
    const auto *msf0_457 = buffer.data(msf0 + 457);
    const auto *msf0_459 = buffer.data(msf0 + 459);
    const auto *msf0_462 = buffer.data(msf0 + 462);
    const auto *msf0_464 = buffer.data(msf0 + 464);
    const auto *msf0_465 = buffer.data(msf0 + 465);
    const auto *msf0_467 = buffer.data(msf0 + 467);
    const auto *msf0_468 = buffer.data(msf0 + 468);
    const auto *msf0_469 = buffer.data(msf0 + 469);
    const auto *msf0_470 = buffer.data(msf0 + 470);
    const auto *msf0_471 = buffer.data(msf0 + 471);
    const auto *msf0_472 = buffer.data(msf0 + 472);
    const auto *msf0_473 = buffer.data(msf0 + 473);
    const auto *msf0_474 = buffer.data(msf0 + 474);
    const auto *msf0_475 = buffer.data(msf0 + 475);
    const auto *msf0_476 = buffer.data(msf0 + 476);
    const auto *msf0_477 = buffer.data(msf0 + 477);
    const auto *msf0_478 = buffer.data(msf0 + 478);
    const auto *msf0_479 = buffer.data(msf0 + 479);
    const auto *msf0_480 = buffer.data(msf0 + 480);
    const auto *msf0_481 = buffer.data(msf0 + 481);
    const auto *msf0_482 = buffer.data(msf0 + 482);
    const auto *msf0_483 = buffer.data(msf0 + 483);
    const auto *msf0_484 = buffer.data(msf0 + 484);
    const auto *msf0_485 = buffer.data(msf0 + 485);
    const auto *msf0_486 = buffer.data(msf0 + 486);
    const auto *msf0_487 = buffer.data(msf0 + 487);
    const auto *msf0_488 = buffer.data(msf0 + 488);
    const auto *msf0_489 = buffer.data(msf0 + 489);
    const auto *msf0_490 = buffer.data(msf0 + 490);
    const auto *msf0_491 = buffer.data(msf0 + 491);
    const auto *msf0_492 = buffer.data(msf0 + 492);
    const auto *msf0_493 = buffer.data(msf0 + 493);
    const auto *msf0_494 = buffer.data(msf0 + 494);
    const auto *msf0_495 = buffer.data(msf0 + 495);
    const auto *msf0_496 = buffer.data(msf0 + 496);
    const auto *msf0_497 = buffer.data(msf0 + 497);
    const auto *msf0_498 = buffer.data(msf0 + 498);
    const auto *msf0_499 = buffer.data(msf0 + 499);
    const auto *msf0_500 = buffer.data(msf0 + 500);
    const auto *msf0_501 = buffer.data(msf0 + 501);
    const auto *msf0_502 = buffer.data(msf0 + 502);
    const auto *msf0_503 = buffer.data(msf0 + 503);
    const auto *msf0_504 = buffer.data(msf0 + 504);
    const auto *msf0_505 = buffer.data(msf0 + 505);
    const auto *msf0_506 = buffer.data(msf0 + 506);
    const auto *msf0_507 = buffer.data(msf0 + 507);
    const auto *msf0_508 = buffer.data(msf0 + 508);
    const auto *msf0_509 = buffer.data(msf0 + 509);
    const auto *msf0_510 = buffer.data(msf0 + 510);
    const auto *msf0_511 = buffer.data(msf0 + 511);
    const auto *msf0_512 = buffer.data(msf0 + 512);
    const auto *msf0_513 = buffer.data(msf0 + 513);
    const auto *msf0_514 = buffer.data(msf0 + 514);
    const auto *msf0_515 = buffer.data(msf0 + 515);
    const auto *msf0_516 = buffer.data(msf0 + 516);
    const auto *msf0_517 = buffer.data(msf0 + 517);

    const auto *msf1_456 = buffer.data(msf1 + 456);
    const auto *msf1_457 = buffer.data(msf1 + 457);
    const auto *msf1_459 = buffer.data(msf1 + 459);
    const auto *msf1_462 = buffer.data(msf1 + 462);
    const auto *msf1_464 = buffer.data(msf1 + 464);
    const auto *msf1_465 = buffer.data(msf1 + 465);
    const auto *msf1_467 = buffer.data(msf1 + 467);
    const auto *msf1_468 = buffer.data(msf1 + 468);
    const auto *msf1_469 = buffer.data(msf1 + 469);
    const auto *msf1_470 = buffer.data(msf1 + 470);
    const auto *msf1_471 = buffer.data(msf1 + 471);
    const auto *msf1_472 = buffer.data(msf1 + 472);
    const auto *msf1_473 = buffer.data(msf1 + 473);
    const auto *msf1_474 = buffer.data(msf1 + 474);
    const auto *msf1_475 = buffer.data(msf1 + 475);
    const auto *msf1_476 = buffer.data(msf1 + 476);
    const auto *msf1_477 = buffer.data(msf1 + 477);
    const auto *msf1_478 = buffer.data(msf1 + 478);
    const auto *msf1_479 = buffer.data(msf1 + 479);
    const auto *msf1_480 = buffer.data(msf1 + 480);
    const auto *msf1_481 = buffer.data(msf1 + 481);
    const auto *msf1_482 = buffer.data(msf1 + 482);
    const auto *msf1_483 = buffer.data(msf1 + 483);
    const auto *msf1_484 = buffer.data(msf1 + 484);
    const auto *msf1_485 = buffer.data(msf1 + 485);
    const auto *msf1_486 = buffer.data(msf1 + 486);
    const auto *msf1_487 = buffer.data(msf1 + 487);
    const auto *msf1_488 = buffer.data(msf1 + 488);
    const auto *msf1_489 = buffer.data(msf1 + 489);
    const auto *msf1_490 = buffer.data(msf1 + 490);
    const auto *msf1_491 = buffer.data(msf1 + 491);
    const auto *msf1_492 = buffer.data(msf1 + 492);
    const auto *msf1_493 = buffer.data(msf1 + 493);
    const auto *msf1_494 = buffer.data(msf1 + 494);
    const auto *msf1_495 = buffer.data(msf1 + 495);
    const auto *msf1_496 = buffer.data(msf1 + 496);
    const auto *msf1_497 = buffer.data(msf1 + 497);
    const auto *msf1_498 = buffer.data(msf1 + 498);
    const auto *msf1_499 = buffer.data(msf1 + 499);
    const auto *msf1_500 = buffer.data(msf1 + 500);
    const auto *msf1_501 = buffer.data(msf1 + 501);
    const auto *msf1_502 = buffer.data(msf1 + 502);
    const auto *msf1_503 = buffer.data(msf1 + 503);
    const auto *msf1_504 = buffer.data(msf1 + 504);
    const auto *msf1_505 = buffer.data(msf1 + 505);
    const auto *msf1_506 = buffer.data(msf1 + 506);
    const auto *msf1_507 = buffer.data(msf1 + 507);
    const auto *msf1_508 = buffer.data(msf1 + 508);
    const auto *msf1_509 = buffer.data(msf1 + 509);
    const auto *msf1_510 = buffer.data(msf1 + 510);
    const auto *msf1_511 = buffer.data(msf1 + 511);
    const auto *msf1_512 = buffer.data(msf1 + 512);
    const auto *msf1_513 = buffer.data(msf1 + 513);
    const auto *msf1_514 = buffer.data(msf1 + 514);
    const auto *msf1_515 = buffer.data(msf1 + 515);
    const auto *msf1_516 = buffer.data(msf1 + 516);
    const auto *msf1_517 = buffer.data(msf1 + 517);

    const auto *msg_685 = buffer.data(msg + 685);
    const auto *msg_686 = buffer.data(msg + 686);
    const auto *msg_687 = buffer.data(msg + 687);
    const auto *msg_689 = buffer.data(msg + 689);
    const auto *msg_692 = buffer.data(msg + 692);
    const auto *msg_694 = buffer.data(msg + 694);
    const auto *msg_695 = buffer.data(msg + 695);
    const auto *msg_697 = buffer.data(msg + 697);
    const auto *msg_698 = buffer.data(msg + 698);
    const auto *msg_699 = buffer.data(msg + 699);
    const auto *msg_700 = buffer.data(msg + 700);
    const auto *msg_701 = buffer.data(msg + 701);
    const auto *msg_702 = buffer.data(msg + 702);
    const auto *msg_703 = buffer.data(msg + 703);
    const auto *msg_704 = buffer.data(msg + 704);
    const auto *msg_705 = buffer.data(msg + 705);
    const auto *msg_706 = buffer.data(msg + 706);
    const auto *msg_707 = buffer.data(msg + 707);
    const auto *msg_708 = buffer.data(msg + 708);
    const auto *msg_709 = buffer.data(msg + 709);
    const auto *msg_710 = buffer.data(msg + 710);
    const auto *msg_711 = buffer.data(msg + 711);
    const auto *msg_712 = buffer.data(msg + 712);
    const auto *msg_713 = buffer.data(msg + 713);
    const auto *msg_714 = buffer.data(msg + 714);
    const auto *msg_715 = buffer.data(msg + 715);
    const auto *msg_716 = buffer.data(msg + 716);
    const auto *msg_717 = buffer.data(msg + 717);
    const auto *msg_718 = buffer.data(msg + 718);
    const auto *msg_719 = buffer.data(msg + 719);
    const auto *msg_720 = buffer.data(msg + 720);
    const auto *msg_721 = buffer.data(msg + 721);
    const auto *msg_722 = buffer.data(msg + 722);
    const auto *msg_723 = buffer.data(msg + 723);
    const auto *msg_724 = buffer.data(msg + 724);
    const auto *msg_725 = buffer.data(msg + 725);
    const auto *msg_726 = buffer.data(msg + 726);
    const auto *msg_727 = buffer.data(msg + 727);
    const auto *msg_728 = buffer.data(msg + 728);
    const auto *msg_729 = buffer.data(msg + 729);
    const auto *msg_730 = buffer.data(msg + 730);
    const auto *msg_731 = buffer.data(msg + 731);
    const auto *msg_732 = buffer.data(msg + 732);
    const auto *msg_733 = buffer.data(msg + 733);
    const auto *msg_734 = buffer.data(msg + 734);
    const auto *msg_735 = buffer.data(msg + 735);
    const auto *msg_736 = buffer.data(msg + 736);
    const auto *msg_737 = buffer.data(msg + 737);
    const auto *msg_738 = buffer.data(msg + 738);
    const auto *msg_739 = buffer.data(msg + 739);
    const auto *msg_740 = buffer.data(msg + 740);
    const auto *msg_741 = buffer.data(msg + 741);
    const auto *msg_742 = buffer.data(msg + 742);
    const auto *msg_743 = buffer.data(msg + 743);
    const auto *msg_744 = buffer.data(msg + 744);
    const auto *msg_745 = buffer.data(msg + 745);
    const auto *msg_746 = buffer.data(msg + 746);
    const auto *msg_747 = buffer.data(msg + 747);
    const auto *msg_748 = buffer.data(msg + 748);
    const auto *msg_749 = buffer.data(msg + 749);
    const auto *msg_750 = buffer.data(msg + 750);
    const auto *msg_751 = buffer.data(msg + 751);
    const auto *msg_752 = buffer.data(msg + 752);
    const auto *msg_753 = buffer.data(msg + 753);
    const auto *msg_754 = buffer.data(msg + 754);
    const auto *msg_755 = buffer.data(msg + 755);
    const auto *msg_756 = buffer.data(msg + 756);
    const auto *msg_757 = buffer.data(msg + 757);
    const auto *msg_758 = buffer.data(msg + 758);
    const auto *msg_759 = buffer.data(msg + 759);
    const auto *msg_760 = buffer.data(msg + 760);
    const auto *msg_761 = buffer.data(msg + 761);
    const auto *msg_762 = buffer.data(msg + 762);
    const auto *msg_763 = buffer.data(msg + 763);
    const auto *msg_764 = buffer.data(msg + 764);
    const auto *msg_765 = buffer.data(msg + 765);
    const auto *msg_766 = buffer.data(msg + 766);
    const auto *msg_767 = buffer.data(msg + 767);
    const auto *msg_768 = buffer.data(msg + 768);
    const auto *msg_769 = buffer.data(msg + 769);
    const auto *msg_770 = buffer.data(msg + 770);
    const auto *msg_771 = buffer.data(msg + 771);
    const auto *msg_772 = buffer.data(msg + 772);

#pragma omp simd aligned(t_961, t_962, t_963, t_964, pc_y, pc_z, lsg_554, msf0_456, msf0_457, \
                         msf1_456, msf1_457, msg_685, msg_686, msg_687, \
                         msg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_3 * pc_z[k] * msg_685[k];

        t_962[k] = f_4 * msf0_456[k]
                   - f_5 * msf1_456[k]
                   + f_3 * pc_z[k] * msg_686[k];

        t_963[k] = f_6 * msf0_457[k]
                   - f_7 * msf1_457[k]
                   + f_3 * pc_z[k] * msg_687[k];

        t_964[k] = f_0 * lsg_554[k]
                   + f_3 * pc_y[k] * msg_689[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, pa_z, pc_z, lsh0_756, lsh0_757, lsh1_756, \
                         lsh1_757, msf0_459, msf1_459, msg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = f_1 * msf0_459[k]
                   - f_2 * msf1_459[k]
                   + f_3 * pc_z[k] * msg_689[k];

        t_966[k] = pa_z[k] * lsh0_756[k]
                   - f_8 * pc_z[k] * lsh1_756[k];

        t_967[k] = pa_z[k] * lsh0_757[k]
                   - f_8 * pc_z[k] * lsh1_757[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, pa_z, pc_x, pc_z, lsh0_759, lsh1_759, msf0_462, \
                         msf0_464, msf1_462, msf1_464, msg_692, \
                         msg_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_13 * msf0_462[k]
                   - f_14 * msf1_462[k]
                   + f_3 * pc_x[k] * msg_692[k];

        t_969[k] = pa_z[k] * lsh0_759[k]
                   - f_8 * pc_z[k] * lsh1_759[k];

        t_970[k] = f_6 * msf0_464[k]
                   - f_7 * msf1_464[k]
                   + f_3 * pc_x[k] * msg_694[k];
    }

#pragma omp simd aligned(t_971, t_972, t_973, pa_z, pc_x, pc_z, lsh0_762, lsh1_762, msf0_465, \
                         msf0_467, msf1_465, msf1_467, msg_695, \
                         msg_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_971[k] = f_6 * msf0_465[k]
                   - f_7 * msf1_465[k]
                   + f_3 * pc_x[k] * msg_695[k];

        t_972[k] = pa_z[k] * lsh0_762[k]
                   - f_8 * pc_z[k] * lsh1_762[k];

        t_973[k] = f_4 * msf0_467[k]
                   - f_5 * msf1_467[k]
                   + f_3 * pc_x[k] * msg_697[k];
    }

#pragma omp simd aligned(t_974, t_975, t_976, t_977, t_978, pc_x, msf0_468, msf0_469, \
                         msf1_468, msf1_469, msg_698, msg_699, msg_700, msg_701, \
                         msg_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_974[k] = f_4 * msf0_468[k]
                   - f_5 * msf1_468[k]
                   + f_3 * pc_x[k] * msg_698[k];

        t_975[k] = f_4 * msf0_469[k]
                   - f_5 * msf1_469[k]
                   + f_3 * pc_x[k] * msg_699[k];

        t_976[k] = f_3 * pc_x[k] * msg_700[k];

        t_977[k] = f_3 * pc_x[k] * msg_701[k];

        t_978[k] = f_3 * pc_x[k] * msg_702[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pa_z, pc_x, pc_z, lsh0_771, lsg_550, \
                         lsh1_771, msg_700, msg_703, msg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_3 * pc_x[k] * msg_703[k];

        t_980[k] = f_3 * pc_x[k] * msg_704[k];

        t_981[k] = pa_z[k] * lsh0_771[k]
                   - f_8 * pc_z[k] * lsh1_771[k];

        t_982[k] = f_9 * lsg_550[k]
                   + f_3 * pc_z[k] * msg_700[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pa_z, pc_y, pc_z, lsh0_773, lsh0_774, lsg_551, \
                         lsg_552, lsg_569, lsh1_773, lsh1_774, \
                         msg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = pa_z[k] * lsh0_773[k]
                   + f_10 * lsg_551[k]
                   - f_8 * pc_z[k] * lsh1_773[k];

        t_984[k] = pa_z[k] * lsh0_774[k]
                   + f_11 * lsg_552[k]
                   - f_8 * pc_z[k] * lsh1_774[k];

        t_985[k] = f_12 * lsg_569[k]
                   + f_3 * pc_y[k] * msg_704[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_x, pc_z, lsg_554, msf0_469, msf0_470, \
                         msf0_471, msf1_469, msf1_470, msf1_471, msg_704, msg_705, \
                         msg_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_9 * lsg_554[k]
                   + f_1 * msf0_469[k]
                   - f_2 * msf1_469[k]
                   + f_3 * pc_z[k] * msg_704[k];

        t_987[k] = f_1 * msf0_470[k]
                   - f_2 * msf1_470[k]
                   + f_3 * pc_x[k] * msg_705[k];

        t_988[k] = f_13 * msf0_471[k]
                   - f_14 * msf1_471[k]
                   + f_3 * pc_x[k] * msg_706[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, msf0_472, msf0_473, msf0_474, msf1_472, \
                         msf1_473, msf1_474, msg_707, msg_708, \
                         msg_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_13 * msf0_472[k]
                   - f_14 * msf1_472[k]
                   + f_3 * pc_x[k] * msg_707[k];

        t_990[k] = f_6 * msf0_473[k]
                   - f_7 * msf1_473[k]
                   + f_3 * pc_x[k] * msg_708[k];

        t_991[k] = f_6 * msf0_474[k]
                   - f_7 * msf1_474[k]
                   + f_3 * pc_x[k] * msg_709[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, pc_x, msf0_475, msf0_476, msf0_477, msf1_475, \
                         msf1_476, msf1_477, msg_710, msg_711, \
                         msg_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_6 * msf0_475[k]
                   - f_7 * msf1_475[k]
                   + f_3 * pc_x[k] * msg_710[k];

        t_993[k] = f_4 * msf0_476[k]
                   - f_5 * msf1_476[k]
                   + f_3 * pc_x[k] * msg_711[k];

        t_994[k] = f_4 * msf0_477[k]
                   - f_5 * msf1_477[k]
                   + f_3 * pc_x[k] * msg_712[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pc_x, msf0_478, msf0_479, \
                         msf1_478, msf1_479, msg_713, msg_714, msg_715, msg_716, \
                         msg_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_4 * msf0_478[k]
                   - f_5 * msf1_478[k]
                   + f_3 * pc_x[k] * msg_713[k];

        t_996[k] = f_4 * msf0_479[k]
                   - f_5 * msf1_479[k]
                   + f_3 * pc_x[k] * msg_714[k];

        t_997[k] = f_3 * pc_x[k] * msg_715[k];

        t_998[k] = f_3 * pc_x[k] * msg_716[k];

        t_999[k] = f_3 * pc_x[k] * msg_717[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pc_x, pc_y, pc_z, lsg_565, lsg_580, \
                         msf0_476, msf1_476, msg_715, msg_718, \
                         msg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_3 * pc_x[k] * msg_718[k];

        t_1001[k] = f_3 * pc_x[k] * msg_719[k];

        t_1002[k] = f_15 * lsg_580[k]
                    + f_1 * msf0_476[k]
                    - f_2 * msf1_476[k]
                    + f_3 * pc_y[k] * msg_715[k];

        t_1003[k] = f_10 * lsg_565[k]
                    + f_3 * pc_z[k] * msg_715[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, pc_y, lsg_582, lsg_583, lsg_584, msf0_478, \
                         msf0_479, msf1_478, msf1_479, msg_717, msg_718, \
                         msg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_15 * lsg_582[k]
                    + f_6 * msf0_478[k]
                    - f_7 * msf1_478[k]
                    + f_3 * pc_y[k] * msg_717[k];

        t_1005[k] = f_15 * lsg_583[k]
                    + f_4 * msf0_479[k]
                    - f_5 * msf1_479[k]
                    + f_3 * pc_y[k] * msg_718[k];

        t_1006[k] = f_15 * lsg_584[k]
                    + f_3 * pc_y[k] * msg_719[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, pc_x, pc_z, lsg_569, msf0_479, msf0_480, \
                         msf0_481, msf1_479, msf1_480, msf1_481, msg_719, msg_720, \
                         msg_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = f_10 * lsg_569[k]
                    + f_1 * msf0_479[k]
                    - f_2 * msf1_479[k]
                    + f_3 * pc_z[k] * msg_719[k];

        t_1008[k] = f_1 * msf0_480[k]
                    - f_2 * msf1_480[k]
                    + f_3 * pc_x[k] * msg_720[k];

        t_1009[k] = f_13 * msf0_481[k]
                    - f_14 * msf1_481[k]
                    + f_3 * pc_x[k] * msg_721[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, pc_x, msf0_482, msf0_483, msf0_484, msf1_482, \
                         msf1_483, msf1_484, msg_722, msg_723, \
                         msg_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_13 * msf0_482[k]
                    - f_14 * msf1_482[k]
                    + f_3 * pc_x[k] * msg_722[k];

        t_1011[k] = f_6 * msf0_483[k]
                    - f_7 * msf1_483[k]
                    + f_3 * pc_x[k] * msg_723[k];

        t_1012[k] = f_6 * msf0_484[k]
                    - f_7 * msf1_484[k]
                    + f_3 * pc_x[k] * msg_724[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pc_x, msf0_485, msf0_486, msf0_487, msf1_485, \
                         msf1_486, msf1_487, msg_725, msg_726, \
                         msg_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_6 * msf0_485[k]
                    - f_7 * msf1_485[k]
                    + f_3 * pc_x[k] * msg_725[k];

        t_1014[k] = f_4 * msf0_486[k]
                    - f_5 * msf1_486[k]
                    + f_3 * pc_x[k] * msg_726[k];

        t_1015[k] = f_4 * msf0_487[k]
                    - f_5 * msf1_487[k]
                    + f_3 * pc_x[k] * msg_727[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, t_1020, pc_x, msf0_488, msf0_489, \
                         msf1_488, msf1_489, msg_728, msg_729, msg_730, msg_731, \
                         msg_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_4 * msf0_488[k]
                    - f_5 * msf1_488[k]
                    + f_3 * pc_x[k] * msg_728[k];

        t_1017[k] = f_4 * msf0_489[k]
                    - f_5 * msf1_489[k]
                    + f_3 * pc_x[k] * msg_729[k];

        t_1018[k] = f_3 * pc_x[k] * msg_730[k];

        t_1019[k] = f_3 * pc_x[k] * msg_731[k];

        t_1020[k] = f_3 * pc_x[k] * msg_732[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pc_x, pc_y, pc_z, lsg_580, lsg_595, \
                         msf0_486, msf1_486, msg_730, msg_733, \
                         msg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_3 * pc_x[k] * msg_733[k];

        t_1022[k] = f_3 * pc_x[k] * msg_734[k];

        t_1023[k] = f_16 * lsg_595[k]
                    + f_1 * msf0_486[k]
                    - f_2 * msf1_486[k]
                    + f_3 * pc_y[k] * msg_730[k];

        t_1024[k] = f_11 * lsg_580[k]
                    + f_3 * pc_z[k] * msg_730[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, pc_y, lsg_597, lsg_598, lsg_599, msf0_488, \
                         msf0_489, msf1_488, msf1_489, msg_732, msg_733, \
                         msg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_16 * lsg_597[k]
                    + f_6 * msf0_488[k]
                    - f_7 * msf1_488[k]
                    + f_3 * pc_y[k] * msg_732[k];

        t_1026[k] = f_16 * lsg_598[k]
                    + f_4 * msf0_489[k]
                    - f_5 * msf1_489[k]
                    + f_3 * pc_y[k] * msg_733[k];

        t_1027[k] = f_16 * lsg_599[k]
                    + f_3 * pc_y[k] * msg_734[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, pc_x, pc_z, lsg_584, msf0_489, msf0_490, \
                         msf0_491, msf1_489, msf1_490, msf1_491, msg_734, msg_735, \
                         msg_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_11 * lsg_584[k]
                    + f_1 * msf0_489[k]
                    - f_2 * msf1_489[k]
                    + f_3 * pc_z[k] * msg_734[k];

        t_1029[k] = f_1 * msf0_490[k]
                    - f_2 * msf1_490[k]
                    + f_3 * pc_x[k] * msg_735[k];

        t_1030[k] = f_13 * msf0_491[k]
                    - f_14 * msf1_491[k]
                    + f_3 * pc_x[k] * msg_736[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, pc_x, msf0_492, msf0_493, msf0_494, msf1_492, \
                         msf1_493, msf1_494, msg_737, msg_738, \
                         msg_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_13 * msf0_492[k]
                    - f_14 * msf1_492[k]
                    + f_3 * pc_x[k] * msg_737[k];

        t_1032[k] = f_6 * msf0_493[k]
                    - f_7 * msf1_493[k]
                    + f_3 * pc_x[k] * msg_738[k];

        t_1033[k] = f_6 * msf0_494[k]
                    - f_7 * msf1_494[k]
                    + f_3 * pc_x[k] * msg_739[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, pc_x, msf0_495, msf0_496, msf0_497, msf1_495, \
                         msf1_496, msf1_497, msg_740, msg_741, \
                         msg_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_6 * msf0_495[k]
                    - f_7 * msf1_495[k]
                    + f_3 * pc_x[k] * msg_740[k];

        t_1035[k] = f_4 * msf0_496[k]
                    - f_5 * msf1_496[k]
                    + f_3 * pc_x[k] * msg_741[k];

        t_1036[k] = f_4 * msf0_497[k]
                    - f_5 * msf1_497[k]
                    + f_3 * pc_x[k] * msg_742[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, t_1041, pc_x, msf0_498, msf0_499, \
                         msf1_498, msf1_499, msg_743, msg_744, msg_745, msg_746, \
                         msg_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_4 * msf0_498[k]
                    - f_5 * msf1_498[k]
                    + f_3 * pc_x[k] * msg_743[k];

        t_1038[k] = f_4 * msf0_499[k]
                    - f_5 * msf1_499[k]
                    + f_3 * pc_x[k] * msg_744[k];

        t_1039[k] = f_3 * pc_x[k] * msg_745[k];

        t_1040[k] = f_3 * pc_x[k] * msg_746[k];

        t_1041[k] = f_3 * pc_x[k] * msg_747[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pc_x, pc_y, pc_z, lsg_595, lsg_610, \
                         msf0_496, msf1_496, msg_745, msg_748, \
                         msg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_3 * pc_x[k] * msg_748[k];

        t_1043[k] = f_3 * pc_x[k] * msg_749[k];

        t_1044[k] = f_17 * lsg_610[k]
                    + f_1 * msf0_496[k]
                    - f_2 * msf1_496[k]
                    + f_3 * pc_y[k] * msg_745[k];

        t_1045[k] = f_18 * lsg_595[k]
                    + f_3 * pc_z[k] * msg_745[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pc_y, lsg_612, lsg_613, lsg_614, msf0_498, \
                         msf0_499, msf1_498, msf1_499, msg_747, msg_748, \
                         msg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_17 * lsg_612[k]
                    + f_6 * msf0_498[k]
                    - f_7 * msf1_498[k]
                    + f_3 * pc_y[k] * msg_747[k];

        t_1047[k] = f_17 * lsg_613[k]
                    + f_4 * msf0_499[k]
                    - f_5 * msf1_499[k]
                    + f_3 * pc_y[k] * msg_748[k];

        t_1048[k] = f_17 * lsg_614[k]
                    + f_3 * pc_y[k] * msg_749[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pc_x, pc_z, lsg_599, msf0_499, msf0_500, \
                         msf0_501, msf1_499, msf1_500, msf1_501, msg_749, msg_750, \
                         msg_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_18 * lsg_599[k]
                    + f_1 * msf0_499[k]
                    - f_2 * msf1_499[k]
                    + f_3 * pc_z[k] * msg_749[k];

        t_1050[k] = f_1 * msf0_500[k]
                    - f_2 * msf1_500[k]
                    + f_3 * pc_x[k] * msg_750[k];

        t_1051[k] = f_13 * msf0_501[k]
                    - f_14 * msf1_501[k]
                    + f_3 * pc_x[k] * msg_751[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pc_x, msf0_502, msf0_503, msf0_504, msf1_502, \
                         msf1_503, msf1_504, msg_752, msg_753, \
                         msg_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_13 * msf0_502[k]
                    - f_14 * msf1_502[k]
                    + f_3 * pc_x[k] * msg_752[k];

        t_1053[k] = f_6 * msf0_503[k]
                    - f_7 * msf1_503[k]
                    + f_3 * pc_x[k] * msg_753[k];

        t_1054[k] = f_6 * msf0_504[k]
                    - f_7 * msf1_504[k]
                    + f_3 * pc_x[k] * msg_754[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, msf0_505, msf0_506, msf0_507, msf1_505, \
                         msf1_506, msf1_507, msg_755, msg_756, \
                         msg_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_6 * msf0_505[k]
                    - f_7 * msf1_505[k]
                    + f_3 * pc_x[k] * msg_755[k];

        t_1056[k] = f_4 * msf0_506[k]
                    - f_5 * msf1_506[k]
                    + f_3 * pc_x[k] * msg_756[k];

        t_1057[k] = f_4 * msf0_507[k]
                    - f_5 * msf1_507[k]
                    + f_3 * pc_x[k] * msg_757[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, t_1061, t_1062, pc_x, msf0_508, msf0_509, \
                         msf1_508, msf1_509, msg_758, msg_759, msg_760, msg_761, \
                         msg_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_4 * msf0_508[k]
                    - f_5 * msf1_508[k]
                    + f_3 * pc_x[k] * msg_758[k];

        t_1059[k] = f_4 * msf0_509[k]
                    - f_5 * msf1_509[k]
                    + f_3 * pc_x[k] * msg_759[k];

        t_1060[k] = f_3 * pc_x[k] * msg_760[k];

        t_1061[k] = f_3 * pc_x[k] * msg_761[k];

        t_1062[k] = f_3 * pc_x[k] * msg_762[k];
    }

#pragma omp simd aligned(t_1063, t_1064, t_1065, t_1066, pc_x, pc_y, pc_z, lsg_610, lsg_625, \
                         msf0_506, msf1_506, msg_760, msg_763, \
                         msg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1063[k] = f_3 * pc_x[k] * msg_763[k];

        t_1064[k] = f_3 * pc_x[k] * msg_764[k];

        t_1065[k] = f_18 * lsg_625[k]
                    + f_1 * msf0_506[k]
                    - f_2 * msf1_506[k]
                    + f_3 * pc_y[k] * msg_760[k];

        t_1066[k] = f_17 * lsg_610[k]
                    + f_3 * pc_z[k] * msg_760[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pc_y, lsg_627, lsg_628, lsg_629, msf0_508, \
                         msf0_509, msf1_508, msf1_509, msg_762, msg_763, \
                         msg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_18 * lsg_627[k]
                    + f_6 * msf0_508[k]
                    - f_7 * msf1_508[k]
                    + f_3 * pc_y[k] * msg_762[k];

        t_1068[k] = f_18 * lsg_628[k]
                    + f_4 * msf0_509[k]
                    - f_5 * msf1_509[k]
                    + f_3 * pc_y[k] * msg_763[k];

        t_1069[k] = f_18 * lsg_629[k]
                    + f_3 * pc_y[k] * msg_764[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_z, lsg_614, msf0_509, msf0_510, \
                         msf0_511, msf1_509, msf1_510, msf1_511, msg_764, msg_765, \
                         msg_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_17 * lsg_614[k]
                    + f_1 * msf0_509[k]
                    - f_2 * msf1_509[k]
                    + f_3 * pc_z[k] * msg_764[k];

        t_1071[k] = f_1 * msf0_510[k]
                    - f_2 * msf1_510[k]
                    + f_3 * pc_x[k] * msg_765[k];

        t_1072[k] = f_13 * msf0_511[k]
                    - f_14 * msf1_511[k]
                    + f_3 * pc_x[k] * msg_766[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, msf0_512, msf0_513, msf0_514, msf1_512, \
                         msf1_513, msf1_514, msg_767, msg_768, \
                         msg_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_13 * msf0_512[k]
                    - f_14 * msf1_512[k]
                    + f_3 * pc_x[k] * msg_767[k];

        t_1074[k] = f_6 * msf0_513[k]
                    - f_7 * msf1_513[k]
                    + f_3 * pc_x[k] * msg_768[k];

        t_1075[k] = f_6 * msf0_514[k]
                    - f_7 * msf1_514[k]
                    + f_3 * pc_x[k] * msg_769[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, msf0_515, msf0_516, msf0_517, msf1_515, \
                         msf1_516, msf1_517, msg_770, msg_771, \
                         msg_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_6 * msf0_515[k]
                    - f_7 * msf1_515[k]
                    + f_3 * pc_x[k] * msg_770[k];

        t_1077[k] = f_4 * msf0_516[k]
                    - f_5 * msf1_516[k]
                    + f_3 * pc_x[k] * msg_771[k];

        t_1078[k] = f_4 * msf0_517[k]
                    - f_5 * msf1_517[k]
                    + f_3 * pc_x[k] * msg_772[k];
    }
}

static auto
compute_prim_msh_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsh0,
                                                          const size_t lsg, const size_t lsh1,
                                                          const size_t msf0, const size_t msf1,
                                                          const size_t msg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_12 = 4.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;

    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh0_924 = buffer.data(lsh0 + 924);
    const auto *lsh0_926 = buffer.data(lsh0 + 926);
    const auto *lsh0_929 = buffer.data(lsh0 + 929);
    const auto *lsh0_933 = buffer.data(lsh0 + 933);
    const auto *lsh0_939 = buffer.data(lsh0 + 939);
    const auto *lsh0_941 = buffer.data(lsh0 + 941);
    const auto *lsh0_942 = buffer.data(lsh0 + 942);
    const auto *lsh0_944 = buffer.data(lsh0 + 944);

    const auto *lsg_625 = buffer.data(lsg + 625);
    const auto *lsg_629 = buffer.data(lsg + 629);
    const auto *lsg_640 = buffer.data(lsg + 640);
    const auto *lsg_642 = buffer.data(lsg + 642);
    const auto *lsg_643 = buffer.data(lsg + 643);
    const auto *lsg_644 = buffer.data(lsg + 644);
    const auto *lsg_655 = buffer.data(lsg + 655);
    const auto *lsg_657 = buffer.data(lsg + 657);
    const auto *lsg_658 = buffer.data(lsg + 658);
    const auto *lsg_659 = buffer.data(lsg + 659);
    const auto *lsg_670 = buffer.data(lsg + 670);
    const auto *lsg_672 = buffer.data(lsg + 672);
    const auto *lsg_673 = buffer.data(lsg + 673);
    const auto *lsg_674 = buffer.data(lsg + 674);

    const auto *lsh1_924 = buffer.data(lsh1 + 924);
    const auto *lsh1_926 = buffer.data(lsh1 + 926);
    const auto *lsh1_929 = buffer.data(lsh1 + 929);
    const auto *lsh1_933 = buffer.data(lsh1 + 933);
    const auto *lsh1_939 = buffer.data(lsh1 + 939);
    const auto *lsh1_941 = buffer.data(lsh1 + 941);
    const auto *lsh1_942 = buffer.data(lsh1 + 942);
    const auto *lsh1_944 = buffer.data(lsh1 + 944);

    const auto *msf0_516 = buffer.data(msf0 + 516);
    const auto *msf0_518 = buffer.data(msf0 + 518);
    const auto *msf0_519 = buffer.data(msf0 + 519);
    const auto *msf0_520 = buffer.data(msf0 + 520);
    const auto *msf0_521 = buffer.data(msf0 + 521);
    const auto *msf0_522 = buffer.data(msf0 + 522);
    const auto *msf0_523 = buffer.data(msf0 + 523);
    const auto *msf0_524 = buffer.data(msf0 + 524);
    const auto *msf0_525 = buffer.data(msf0 + 525);
    const auto *msf0_526 = buffer.data(msf0 + 526);
    const auto *msf0_527 = buffer.data(msf0 + 527);
    const auto *msf0_528 = buffer.data(msf0 + 528);
    const auto *msf0_529 = buffer.data(msf0 + 529);
    const auto *msf0_531 = buffer.data(msf0 + 531);
    const auto *msf0_533 = buffer.data(msf0 + 533);
    const auto *msf0_534 = buffer.data(msf0 + 534);
    const auto *msf0_536 = buffer.data(msf0 + 536);
    const auto *msf0_537 = buffer.data(msf0 + 537);
    const auto *msf0_538 = buffer.data(msf0 + 538);
    const auto *msf0_540 = buffer.data(msf0 + 540);
    const auto *msf0_542 = buffer.data(msf0 + 542);
    const auto *msf0_543 = buffer.data(msf0 + 543);
    const auto *msf0_545 = buffer.data(msf0 + 545);
    const auto *msf0_546 = buffer.data(msf0 + 546);
    const auto *msf0_547 = buffer.data(msf0 + 547);
    const auto *msf0_548 = buffer.data(msf0 + 548);
    const auto *msf0_549 = buffer.data(msf0 + 549);

    const auto *msf1_516 = buffer.data(msf1 + 516);
    const auto *msf1_518 = buffer.data(msf1 + 518);
    const auto *msf1_519 = buffer.data(msf1 + 519);
    const auto *msf1_520 = buffer.data(msf1 + 520);
    const auto *msf1_521 = buffer.data(msf1 + 521);
    const auto *msf1_522 = buffer.data(msf1 + 522);
    const auto *msf1_523 = buffer.data(msf1 + 523);
    const auto *msf1_524 = buffer.data(msf1 + 524);
    const auto *msf1_525 = buffer.data(msf1 + 525);
    const auto *msf1_526 = buffer.data(msf1 + 526);
    const auto *msf1_527 = buffer.data(msf1 + 527);
    const auto *msf1_528 = buffer.data(msf1 + 528);
    const auto *msf1_529 = buffer.data(msf1 + 529);
    const auto *msf1_531 = buffer.data(msf1 + 531);
    const auto *msf1_533 = buffer.data(msf1 + 533);
    const auto *msf1_534 = buffer.data(msf1 + 534);
    const auto *msf1_536 = buffer.data(msf1 + 536);
    const auto *msf1_537 = buffer.data(msf1 + 537);
    const auto *msf1_538 = buffer.data(msf1 + 538);
    const auto *msf1_540 = buffer.data(msf1 + 540);
    const auto *msf1_542 = buffer.data(msf1 + 542);
    const auto *msf1_543 = buffer.data(msf1 + 543);
    const auto *msf1_545 = buffer.data(msf1 + 545);
    const auto *msf1_546 = buffer.data(msf1 + 546);
    const auto *msf1_547 = buffer.data(msf1 + 547);
    const auto *msf1_548 = buffer.data(msf1 + 548);
    const auto *msf1_549 = buffer.data(msf1 + 549);

    const auto *msg_773 = buffer.data(msg + 773);
    const auto *msg_774 = buffer.data(msg + 774);
    const auto *msg_775 = buffer.data(msg + 775);
    const auto *msg_776 = buffer.data(msg + 776);
    const auto *msg_777 = buffer.data(msg + 777);
    const auto *msg_778 = buffer.data(msg + 778);
    const auto *msg_779 = buffer.data(msg + 779);
    const auto *msg_780 = buffer.data(msg + 780);
    const auto *msg_781 = buffer.data(msg + 781);
    const auto *msg_782 = buffer.data(msg + 782);
    const auto *msg_783 = buffer.data(msg + 783);
    const auto *msg_784 = buffer.data(msg + 784);
    const auto *msg_785 = buffer.data(msg + 785);
    const auto *msg_786 = buffer.data(msg + 786);
    const auto *msg_787 = buffer.data(msg + 787);
    const auto *msg_788 = buffer.data(msg + 788);
    const auto *msg_789 = buffer.data(msg + 789);
    const auto *msg_790 = buffer.data(msg + 790);
    const auto *msg_791 = buffer.data(msg + 791);
    const auto *msg_792 = buffer.data(msg + 792);
    const auto *msg_793 = buffer.data(msg + 793);
    const auto *msg_794 = buffer.data(msg + 794);
    const auto *msg_796 = buffer.data(msg + 796);
    const auto *msg_798 = buffer.data(msg + 798);
    const auto *msg_799 = buffer.data(msg + 799);
    const auto *msg_801 = buffer.data(msg + 801);
    const auto *msg_802 = buffer.data(msg + 802);
    const auto *msg_803 = buffer.data(msg + 803);
    const auto *msg_805 = buffer.data(msg + 805);
    const auto *msg_806 = buffer.data(msg + 806);
    const auto *msg_807 = buffer.data(msg + 807);
    const auto *msg_808 = buffer.data(msg + 808);
    const auto *msg_809 = buffer.data(msg + 809);
    const auto *msg_810 = buffer.data(msg + 810);
    const auto *msg_812 = buffer.data(msg + 812);
    const auto *msg_813 = buffer.data(msg + 813);
    const auto *msg_815 = buffer.data(msg + 815);
    const auto *msg_816 = buffer.data(msg + 816);
    const auto *msg_817 = buffer.data(msg + 817);
    const auto *msg_819 = buffer.data(msg + 819);
    const auto *msg_820 = buffer.data(msg + 820);
    const auto *msg_821 = buffer.data(msg + 821);
    const auto *msg_822 = buffer.data(msg + 822);
    const auto *msg_823 = buffer.data(msg + 823);
    const auto *msg_824 = buffer.data(msg + 824);

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, t_1083, pc_x, msf0_518, msf0_519, \
                         msf1_518, msf1_519, msg_773, msg_774, msg_775, msg_776, \
                         msg_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_4 * msf0_518[k]
                    - f_5 * msf1_518[k]
                    + f_3 * pc_x[k] * msg_773[k];

        t_1080[k] = f_4 * msf0_519[k]
                    - f_5 * msf1_519[k]
                    + f_3 * pc_x[k] * msg_774[k];

        t_1081[k] = f_3 * pc_x[k] * msg_775[k];

        t_1082[k] = f_3 * pc_x[k] * msg_776[k];

        t_1083[k] = f_3 * pc_x[k] * msg_777[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, t_1087, pc_x, pc_y, pc_z, lsg_625, lsg_640, \
                         msf0_516, msf1_516, msg_775, msg_778, \
                         msg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_3 * pc_x[k] * msg_778[k];

        t_1085[k] = f_3 * pc_x[k] * msg_779[k];

        t_1086[k] = f_11 * lsg_640[k]
                    + f_1 * msf0_516[k]
                    - f_2 * msf1_516[k]
                    + f_3 * pc_y[k] * msg_775[k];

        t_1087[k] = f_16 * lsg_625[k]
                    + f_3 * pc_z[k] * msg_775[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_y, lsg_642, lsg_643, lsg_644, msf0_518, \
                         msf0_519, msf1_518, msf1_519, msg_777, msg_778, \
                         msg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_11 * lsg_642[k]
                    + f_6 * msf0_518[k]
                    - f_7 * msf1_518[k]
                    + f_3 * pc_y[k] * msg_777[k];

        t_1089[k] = f_11 * lsg_643[k]
                    + f_4 * msf0_519[k]
                    - f_5 * msf1_519[k]
                    + f_3 * pc_y[k] * msg_778[k];

        t_1090[k] = f_11 * lsg_644[k]
                    + f_3 * pc_y[k] * msg_779[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_z, lsg_629, msf0_519, msf0_520, \
                         msf0_521, msf1_519, msf1_520, msf1_521, msg_779, msg_780, \
                         msg_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_16 * lsg_629[k]
                    + f_1 * msf0_519[k]
                    - f_2 * msf1_519[k]
                    + f_3 * pc_z[k] * msg_779[k];

        t_1092[k] = f_1 * msf0_520[k]
                    - f_2 * msf1_520[k]
                    + f_3 * pc_x[k] * msg_780[k];

        t_1093[k] = f_13 * msf0_521[k]
                    - f_14 * msf1_521[k]
                    + f_3 * pc_x[k] * msg_781[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, msf0_522, msf0_523, msf0_524, msf1_522, \
                         msf1_523, msf1_524, msg_782, msg_783, \
                         msg_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_13 * msf0_522[k]
                    - f_14 * msf1_522[k]
                    + f_3 * pc_x[k] * msg_782[k];

        t_1095[k] = f_6 * msf0_523[k]
                    - f_7 * msf1_523[k]
                    + f_3 * pc_x[k] * msg_783[k];

        t_1096[k] = f_6 * msf0_524[k]
                    - f_7 * msf1_524[k]
                    + f_3 * pc_x[k] * msg_784[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, msf0_525, msf0_526, msf0_527, msf1_525, \
                         msf1_526, msf1_527, msg_785, msg_786, \
                         msg_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_6 * msf0_525[k]
                    - f_7 * msf1_525[k]
                    + f_3 * pc_x[k] * msg_785[k];

        t_1098[k] = f_4 * msf0_526[k]
                    - f_5 * msf1_526[k]
                    + f_3 * pc_x[k] * msg_786[k];

        t_1099[k] = f_4 * msf0_527[k]
                    - f_5 * msf1_527[k]
                    + f_3 * pc_x[k] * msg_787[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, t_1104, pc_x, msf0_528, msf0_529, \
                         msf1_528, msf1_529, msg_788, msg_789, msg_790, msg_791, \
                         msg_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_4 * msf0_528[k]
                    - f_5 * msf1_528[k]
                    + f_3 * pc_x[k] * msg_788[k];

        t_1101[k] = f_4 * msf0_529[k]
                    - f_5 * msf1_529[k]
                    + f_3 * pc_x[k] * msg_789[k];

        t_1102[k] = f_3 * pc_x[k] * msg_790[k];

        t_1103[k] = f_3 * pc_x[k] * msg_791[k];

        t_1104[k] = f_3 * pc_x[k] * msg_792[k];
    }

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, pc_x, pc_y, pc_z, lsg_640, lsg_655, \
                         msf0_526, msf1_526, msg_790, msg_793, \
                         msg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = f_3 * pc_x[k] * msg_793[k];

        t_1106[k] = f_3 * pc_x[k] * msg_794[k];

        t_1107[k] = f_10 * lsg_655[k]
                    + f_1 * msf0_526[k]
                    - f_2 * msf1_526[k]
                    + f_3 * pc_y[k] * msg_790[k];

        t_1108[k] = f_15 * lsg_640[k]
                    + f_3 * pc_z[k] * msg_790[k];
    }

#pragma omp simd aligned(t_1109, t_1110, t_1111, pc_y, lsg_657, lsg_658, lsg_659, msf0_528, \
                         msf0_529, msf1_528, msf1_529, msg_792, msg_793, \
                         msg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1109[k] = f_10 * lsg_657[k]
                    + f_6 * msf0_528[k]
                    - f_7 * msf1_528[k]
                    + f_3 * pc_y[k] * msg_792[k];

        t_1110[k] = f_10 * lsg_658[k]
                    + f_4 * msf0_529[k]
                    - f_5 * msf1_529[k]
                    + f_3 * pc_y[k] * msg_793[k];

        t_1111[k] = f_10 * lsg_659[k]
                    + f_3 * pc_y[k] * msg_794[k];
    }

#pragma omp simd aligned(t_1112, t_1113, t_1114, pa_y, pc_x, pc_y, pc_z, lsh0_924, lsg_644, \
                         lsh1_924, msf0_529, msf0_531, msf1_529, msf1_531, msg_794, \
                         msg_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1112[k] = f_15 * lsg_644[k]
                    + f_1 * msf0_529[k]
                    - f_2 * msf1_529[k]
                    + f_3 * pc_z[k] * msg_794[k];

        t_1113[k] = pa_y[k] * lsh0_924[k]
                    - f_8 * pc_y[k] * lsh1_924[k];

        t_1114[k] = f_13 * msf0_531[k]
                    - f_14 * msf1_531[k]
                    + f_3 * pc_x[k] * msg_796[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, pa_y, pc_x, pc_y, lsh0_926, lsh1_926, \
                         msf0_533, msf0_534, msf1_533, msf1_534, msg_798, \
                         msg_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = pa_y[k] * lsh0_926[k]
                    - f_8 * pc_y[k] * lsh1_926[k];

        t_1116[k] = f_6 * msf0_533[k]
                    - f_7 * msf1_533[k]
                    + f_3 * pc_x[k] * msg_798[k];

        t_1117[k] = f_6 * msf0_534[k]
                    - f_7 * msf1_534[k]
                    + f_3 * pc_x[k] * msg_799[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pa_y, pc_x, pc_y, lsh0_929, lsh1_929, \
                         msf0_536, msf0_537, msf1_536, msf1_537, msg_801, \
                         msg_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = pa_y[k] * lsh0_929[k]
                    - f_8 * pc_y[k] * lsh1_929[k];

        t_1119[k] = f_4 * msf0_536[k]
                    - f_5 * msf1_536[k]
                    + f_3 * pc_x[k] * msg_801[k];

        t_1120[k] = f_4 * msf0_537[k]
                    - f_5 * msf1_537[k]
                    + f_3 * pc_x[k] * msg_802[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, t_1124, t_1125, pa_y, pc_x, pc_y, lsh0_933, \
                         lsh1_933, msf0_538, msf1_538, msg_803, msg_805, msg_806, \
                         msg_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_4 * msf0_538[k]
                    - f_5 * msf1_538[k]
                    + f_3 * pc_x[k] * msg_803[k];

        t_1122[k] = pa_y[k] * lsh0_933[k]
                    - f_8 * pc_y[k] * lsh1_933[k];

        t_1123[k] = f_3 * pc_x[k] * msg_805[k];

        t_1124[k] = f_3 * pc_x[k] * msg_806[k];

        t_1125[k] = f_3 * pc_x[k] * msg_807[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, t_1129, pa_y, pc_x, pc_y, pc_z, lsh0_939, \
                         lsg_655, lsg_670, lsh1_939, msg_805, msg_808, \
                         msg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_3 * pc_x[k] * msg_808[k];

        t_1127[k] = f_3 * pc_x[k] * msg_809[k];

        t_1128[k] = pa_y[k] * lsh0_939[k]
                    + f_17 * lsg_670[k]
                    - f_8 * pc_y[k] * lsh1_939[k];

        t_1129[k] = f_12 * lsg_655[k]
                    + f_3 * pc_z[k] * msg_805[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, t_1133, pa_y, pc_y, lsh0_941, lsh0_942, \
                         lsh0_944, lsg_672, lsg_673, lsg_674, lsh1_941, lsh1_942, lsh1_944, \
                         msg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = pa_y[k] * lsh0_941[k]
                    + f_11 * lsg_672[k]
                    - f_8 * pc_y[k] * lsh1_941[k];

        t_1131[k] = pa_y[k] * lsh0_942[k]
                    + f_10 * lsg_673[k]
                    - f_8 * pc_y[k] * lsh1_942[k];

        t_1132[k] = f_9 * lsg_674[k]
                    + f_3 * pc_y[k] * msg_809[k];

        t_1133[k] = pa_y[k] * lsh0_944[k]
                    - f_8 * pc_y[k] * lsh1_944[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, t_1138, pc_x, pc_y, msf0_540, \
                         msf0_542, msf0_543, msf1_540, msf1_542, msf1_543, msg_810, msg_812, \
                         msg_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_1 * msf0_540[k]
                    - f_2 * msf1_540[k]
                    + f_3 * pc_x[k] * msg_810[k];

        t_1135[k] = f_3 * pc_y[k] * msg_810[k];

        t_1136[k] = f_13 * msf0_542[k]
                    - f_14 * msf1_542[k]
                    + f_3 * pc_x[k] * msg_812[k];

        t_1137[k] = f_6 * msf0_543[k]
                    - f_7 * msf1_543[k]
                    + f_3 * pc_x[k] * msg_813[k];

        t_1138[k] = f_3 * pc_y[k] * msg_812[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, t_1142, pc_x, pc_y, msf0_545, msf0_546, \
                         msf0_547, msf1_545, msf1_546, msf1_547, msg_815, msg_816, \
                         msg_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_6 * msf0_545[k]
                    - f_7 * msf1_545[k]
                    + f_3 * pc_x[k] * msg_815[k];

        t_1140[k] = f_4 * msf0_546[k]
                    - f_5 * msf1_546[k]
                    + f_3 * pc_x[k] * msg_816[k];

        t_1141[k] = f_4 * msf0_547[k]
                    - f_5 * msf1_547[k]
                    + f_3 * pc_x[k] * msg_817[k];

        t_1142[k] = f_3 * pc_y[k] * msg_815[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, t_1146, t_1147, t_1148, pc_x, msf0_549, \
                         msf1_549, msg_819, msg_820, msg_821, msg_822, msg_823, \
                         msg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_4 * msf0_549[k]
                    - f_5 * msf1_549[k]
                    + f_3 * pc_x[k] * msg_819[k];

        t_1144[k] = f_3 * pc_x[k] * msg_820[k];

        t_1145[k] = f_3 * pc_x[k] * msg_821[k];

        t_1146[k] = f_3 * pc_x[k] * msg_822[k];

        t_1147[k] = f_3 * pc_x[k] * msg_823[k];

        t_1148[k] = f_3 * pc_x[k] * msg_824[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, pc_y, msf0_546, msf0_547, msf0_548, msf1_546, \
                         msf1_547, msf1_548, msg_820, msg_821, \
                         msg_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_1 * msf0_546[k]
                    - f_2 * msf1_546[k]
                    + f_3 * pc_y[k] * msg_820[k];

        t_1150[k] = f_13 * msf0_547[k]
                    - f_14 * msf1_547[k]
                    + f_3 * pc_y[k] * msg_821[k];

        t_1151[k] = f_6 * msf0_548[k]
                    - f_7 * msf1_548[k]
                    + f_3 * pc_y[k] * msg_822[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, pc_y, pc_z, lsg_674, msf0_549, msf1_549, \
                         msg_823, msg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = f_4 * msf0_549[k]
                    - f_5 * msf1_549[k]
                    + f_3 * pc_y[k] * msg_823[k];

        t_1153[k] = f_3 * pc_y[k] * msg_824[k];

        t_1154[k] = f_0 * lsg_674[k]
                    + f_1 * msf0_549[k]
                    - f_2 * msf1_549[k]
                    + f_3 * pc_z[k] * msg_824[k];
    }
}

auto
compute_prim_msh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t lsh0, const size_t lsg,
                                                   const size_t lsh1, const size_t msf0,
                                                   const size_t msf1, const size_t msg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_msh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);

    compute_prim_msh_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, lsh0, lsg,
                                                              lsh1, msf0, msf1, msg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
