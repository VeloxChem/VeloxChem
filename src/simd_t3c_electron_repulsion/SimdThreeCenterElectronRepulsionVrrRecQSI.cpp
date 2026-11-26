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


#include "SimdThreeCenterElectronRepulsionVrrRecQSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 5.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_0 = buffer.data(osi0 + 0);
    const auto *osi0_3 = buffer.data(osi0 + 3);
    const auto *osi0_5 = buffer.data(osi0 + 5);
    const auto *osi0_6 = buffer.data(osi0 + 6);
    const auto *osi0_9 = buffer.data(osi0 + 9);
    const auto *osi0_10 = buffer.data(osi0 + 10);
    const auto *osi0_14 = buffer.data(osi0 + 14);
    const auto *osi0_21 = buffer.data(osi0 + 21);
    const auto *osi0_27 = buffer.data(osi0 + 27);
    const auto *osi0_31 = buffer.data(osi0 + 31);
    const auto *osi0_34 = buffer.data(osi0 + 34);
    const auto *osi0_38 = buffer.data(osi0 + 38);
    const auto *osi0_56 = buffer.data(osi0 + 56);
    const auto *osi0_61 = buffer.data(osi0 + 61);
    const auto *osi0_65 = buffer.data(osi0 + 65);
    const auto *osi0_68 = buffer.data(osi0 + 68);
    const auto *osi0_70 = buffer.data(osi0 + 70);

    const auto *osh_0 = buffer.data(osh + 0);
    const auto *osh_1 = buffer.data(osh + 1);
    const auto *osh_2 = buffer.data(osh + 2);
    const auto *osh_3 = buffer.data(osh + 3);
    const auto *osh_5 = buffer.data(osh + 5);
    const auto *osh_6 = buffer.data(osh + 6);
    const auto *osh_9 = buffer.data(osh + 9);
    const auto *osh_15 = buffer.data(osh + 15);
    const auto *osh_17 = buffer.data(osh + 17);
    const auto *osh_18 = buffer.data(osh + 18);
    const auto *osh_20 = buffer.data(osh + 20);
    const auto *osh_21 = buffer.data(osh + 21);
    const auto *osh_24 = buffer.data(osh + 24);
    const auto *osh_26 = buffer.data(osh + 26);
    const auto *osh_27 = buffer.data(osh + 27);
    const auto *osh_30 = buffer.data(osh + 30);
    const auto *osh_36 = buffer.data(osh + 36);
    const auto *osh_38 = buffer.data(osh + 38);
    const auto *osh_39 = buffer.data(osh + 39);
    const auto *osh_40 = buffer.data(osh + 40);
    const auto *osh_41 = buffer.data(osh + 41);
    const auto *osh_42 = buffer.data(osh + 42);
    const auto *osh_44 = buffer.data(osh + 44);
    const auto *osh_47 = buffer.data(osh + 47);
    const auto *osh_50 = buffer.data(osh + 50);
    const auto *osh_51 = buffer.data(osh + 51);
    const auto *osh_57 = buffer.data(osh + 57);
    const auto *osh_58 = buffer.data(osh + 58);
    const auto *osh_59 = buffer.data(osh + 59);
    const auto *osh_60 = buffer.data(osh + 60);
    const auto *osh_62 = buffer.data(osh + 62);
    const auto *osh_63 = buffer.data(osh + 63);
    const auto *osh_66 = buffer.data(osh + 66);
    const auto *osh_69 = buffer.data(osh + 69);
    const auto *osh_73 = buffer.data(osh + 73);
    const auto *osh_78 = buffer.data(osh + 78);
    const auto *osh_80 = buffer.data(osh + 80);
    const auto *osh_81 = buffer.data(osh + 81);
    const auto *osh_82 = buffer.data(osh + 82);
    const auto *osh_83 = buffer.data(osh + 83);
    const auto *osh_99 = buffer.data(osh + 99);

    const auto *osi1_0 = buffer.data(osi1 + 0);
    const auto *osi1_3 = buffer.data(osi1 + 3);
    const auto *osi1_5 = buffer.data(osi1 + 5);
    const auto *osi1_6 = buffer.data(osi1 + 6);
    const auto *osi1_9 = buffer.data(osi1 + 9);
    const auto *osi1_10 = buffer.data(osi1 + 10);
    const auto *osi1_14 = buffer.data(osi1 + 14);
    const auto *osi1_21 = buffer.data(osi1 + 21);
    const auto *osi1_27 = buffer.data(osi1 + 27);
    const auto *osi1_31 = buffer.data(osi1 + 31);
    const auto *osi1_34 = buffer.data(osi1 + 34);
    const auto *osi1_38 = buffer.data(osi1 + 38);
    const auto *osi1_56 = buffer.data(osi1 + 56);
    const auto *osi1_61 = buffer.data(osi1 + 61);
    const auto *osi1_65 = buffer.data(osi1 + 65);
    const auto *osi1_68 = buffer.data(osi1 + 68);
    const auto *osi1_70 = buffer.data(osi1 + 70);

    const auto *qsg0_0 = buffer.data(qsg0 + 0);
    const auto *qsg0_1 = buffer.data(qsg0 + 1);
    const auto *qsg0_2 = buffer.data(qsg0 + 2);
    const auto *qsg0_3 = buffer.data(qsg0 + 3);
    const auto *qsg0_5 = buffer.data(qsg0 + 5);
    const auto *qsg0_10 = buffer.data(qsg0 + 10);
    const auto *qsg0_12 = buffer.data(qsg0 + 12);
    const auto *qsg0_13 = buffer.data(qsg0 + 13);
    const auto *qsg0_14 = buffer.data(qsg0 + 14);
    const auto *qsg0_18 = buffer.data(qsg0 + 18);
    const auto *qsg0_25 = buffer.data(qsg0 + 25);
    const auto *qsg0_26 = buffer.data(qsg0 + 26);
    const auto *qsg0_27 = buffer.data(qsg0 + 27);
    const auto *qsg0_32 = buffer.data(qsg0 + 32);
    const auto *qsg0_34 = buffer.data(qsg0 + 34);
    const auto *qsg0_35 = buffer.data(qsg0 + 35);
    const auto *qsg0_41 = buffer.data(qsg0 + 41);
    const auto *qsg0_42 = buffer.data(qsg0 + 42);
    const auto *qsg0_43 = buffer.data(qsg0 + 43);
    const auto *qsg0_44 = buffer.data(qsg0 + 44);
    const auto *qsg0_45 = buffer.data(qsg0 + 45);
    const auto *qsg0_47 = buffer.data(qsg0 + 47);
    const auto *qsg0_48 = buffer.data(qsg0 + 48);
    const auto *qsg0_50 = buffer.data(qsg0 + 50);
    const auto *qsg0_51 = buffer.data(qsg0 + 51);
    const auto *qsg0_55 = buffer.data(qsg0 + 55);
    const auto *qsg0_56 = buffer.data(qsg0 + 56);
    const auto *qsg0_57 = buffer.data(qsg0 + 57);
    const auto *qsg0_59 = buffer.data(qsg0 + 59);

    const auto *qsg1_0 = buffer.data(qsg1 + 0);
    const auto *qsg1_1 = buffer.data(qsg1 + 1);
    const auto *qsg1_2 = buffer.data(qsg1 + 2);
    const auto *qsg1_3 = buffer.data(qsg1 + 3);
    const auto *qsg1_5 = buffer.data(qsg1 + 5);
    const auto *qsg1_10 = buffer.data(qsg1 + 10);
    const auto *qsg1_12 = buffer.data(qsg1 + 12);
    const auto *qsg1_13 = buffer.data(qsg1 + 13);
    const auto *qsg1_14 = buffer.data(qsg1 + 14);
    const auto *qsg1_18 = buffer.data(qsg1 + 18);
    const auto *qsg1_25 = buffer.data(qsg1 + 25);
    const auto *qsg1_26 = buffer.data(qsg1 + 26);
    const auto *qsg1_27 = buffer.data(qsg1 + 27);
    const auto *qsg1_32 = buffer.data(qsg1 + 32);
    const auto *qsg1_34 = buffer.data(qsg1 + 34);
    const auto *qsg1_35 = buffer.data(qsg1 + 35);
    const auto *qsg1_41 = buffer.data(qsg1 + 41);
    const auto *qsg1_42 = buffer.data(qsg1 + 42);
    const auto *qsg1_43 = buffer.data(qsg1 + 43);
    const auto *qsg1_44 = buffer.data(qsg1 + 44);
    const auto *qsg1_45 = buffer.data(qsg1 + 45);
    const auto *qsg1_47 = buffer.data(qsg1 + 47);
    const auto *qsg1_48 = buffer.data(qsg1 + 48);
    const auto *qsg1_50 = buffer.data(qsg1 + 50);
    const auto *qsg1_51 = buffer.data(qsg1 + 51);
    const auto *qsg1_55 = buffer.data(qsg1 + 55);
    const auto *qsg1_56 = buffer.data(qsg1 + 56);
    const auto *qsg1_57 = buffer.data(qsg1 + 57);
    const auto *qsg1_59 = buffer.data(qsg1 + 59);

    const auto *qsh_0 = buffer.data(qsh + 0);
    const auto *qsh_1 = buffer.data(qsh + 1);
    const auto *qsh_2 = buffer.data(qsh + 2);
    const auto *qsh_3 = buffer.data(qsh + 3);
    const auto *qsh_5 = buffer.data(qsh + 5);
    const auto *qsh_6 = buffer.data(qsh + 6);
    const auto *qsh_8 = buffer.data(qsh + 8);
    const auto *qsh_9 = buffer.data(qsh + 9);
    const auto *qsh_10 = buffer.data(qsh + 10);
    const auto *qsh_14 = buffer.data(qsh + 14);
    const auto *qsh_15 = buffer.data(qsh + 15);
    const auto *qsh_17 = buffer.data(qsh + 17);
    const auto *qsh_18 = buffer.data(qsh + 18);
    const auto *qsh_19 = buffer.data(qsh + 19);
    const auto *qsh_20 = buffer.data(qsh + 20);
    const auto *qsh_21 = buffer.data(qsh + 21);
    const auto *qsh_22 = buffer.data(qsh + 22);
    const auto *qsh_24 = buffer.data(qsh + 24);
    const auto *qsh_26 = buffer.data(qsh + 26);
    const auto *qsh_27 = buffer.data(qsh + 27);
    const auto *qsh_28 = buffer.data(qsh + 28);
    const auto *qsh_30 = buffer.data(qsh + 30);
    const auto *qsh_31 = buffer.data(qsh + 31);
    const auto *qsh_36 = buffer.data(qsh + 36);
    const auto *qsh_37 = buffer.data(qsh + 37);
    const auto *qsh_38 = buffer.data(qsh + 38);
    const auto *qsh_39 = buffer.data(qsh + 39);
    const auto *qsh_40 = buffer.data(qsh + 40);
    const auto *qsh_41 = buffer.data(qsh + 41);
    const auto *qsh_42 = buffer.data(qsh + 42);
    const auto *qsh_44 = buffer.data(qsh + 44);
    const auto *qsh_46 = buffer.data(qsh + 46);
    const auto *qsh_47 = buffer.data(qsh + 47);
    const auto *qsh_49 = buffer.data(qsh + 49);
    const auto *qsh_50 = buffer.data(qsh + 50);
    const auto *qsh_51 = buffer.data(qsh + 51);
    const auto *qsh_56 = buffer.data(qsh + 56);
    const auto *qsh_57 = buffer.data(qsh + 57);
    const auto *qsh_58 = buffer.data(qsh + 58);
    const auto *qsh_59 = buffer.data(qsh + 59);
    const auto *qsh_60 = buffer.data(qsh + 60);
    const auto *qsh_61 = buffer.data(qsh + 61);
    const auto *qsh_62 = buffer.data(qsh + 62);
    const auto *qsh_63 = buffer.data(qsh + 63);
    const auto *qsh_64 = buffer.data(qsh + 64);
    const auto *qsh_65 = buffer.data(qsh + 65);
    const auto *qsh_66 = buffer.data(qsh + 66);
    const auto *qsh_68 = buffer.data(qsh + 68);
    const auto *qsh_69 = buffer.data(qsh + 69);
    const auto *qsh_70 = buffer.data(qsh + 70);
    const auto *qsh_72 = buffer.data(qsh + 72);
    const auto *qsh_73 = buffer.data(qsh + 73);
    const auto *qsh_78 = buffer.data(qsh + 78);
    const auto *qsh_79 = buffer.data(qsh + 79);
    const auto *qsh_80 = buffer.data(qsh + 80);
    const auto *qsh_81 = buffer.data(qsh + 81);
    const auto *qsh_82 = buffer.data(qsh + 82);
    const auto *qsh_83 = buffer.data(qsh + 83);
    const auto *qsh_84 = buffer.data(qsh + 84);
    const auto *qsh_86 = buffer.data(qsh + 86);
    const auto *qsh_87 = buffer.data(qsh + 87);
    const auto *qsh_89 = buffer.data(qsh + 89);
    const auto *qsh_90 = buffer.data(qsh + 90);
    const auto *qsh_93 = buffer.data(qsh + 93);
    const auto *qsh_99 = buffer.data(qsh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, osh_0, qsg0_0, \
                         qsg1_0, qsh_0, qsh_1, qsh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * osh_0[k]
                 + f_1 * qsg0_0[k]
                 - f_2 * qsg1_0[k]
                 + f_3 * pc_x[k] * qsh_0[k];

        t_1[k] = f_3 * pc_y[k] * qsh_0[k];

        t_2[k] = f_3 * pc_z[k] * qsh_0[k];

        t_3[k] = f_4 * qsg0_0[k]
                 - f_5 * qsg1_0[k]
                 + f_3 * pc_y[k] * qsh_1[k];

        t_4[k] = f_3 * pc_y[k] * qsh_2[k];

        t_5[k] = f_4 * qsg0_0[k]
                 - f_5 * qsg1_0[k]
                 + f_3 * pc_z[k] * qsh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, qsg0_1, qsg0_2, qsg0_3, qsg1_1, \
                         qsg1_2, qsg1_3, qsh_3, qsh_5, qsh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * qsg0_1[k]
                 - f_7 * qsg1_1[k]
                 + f_3 * pc_y[k] * qsh_3[k];

        t_7[k] = f_3 * pc_z[k] * qsh_3[k];

        t_8[k] = f_3 * pc_y[k] * qsh_5[k];

        t_9[k] = f_6 * qsg0_2[k]
                 - f_7 * qsg1_2[k]
                 + f_3 * pc_z[k] * qsh_5[k];

        t_10[k] = f_8 * qsg0_3[k]
                  - f_9 * qsg1_3[k]
                  + f_3 * pc_y[k] * qsh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, osh_15, qsg0_5, \
                         qsg1_5, qsh_6, qsh_8, qsh_9, qsh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * qsh_6[k];

        t_12[k] = f_4 * qsg0_5[k]
                  - f_5 * qsg1_5[k]
                  + f_3 * pc_y[k] * qsh_8[k];

        t_13[k] = f_3 * pc_y[k] * qsh_9[k];

        t_14[k] = f_8 * qsg0_5[k]
                  - f_9 * qsg1_5[k]
                  + f_3 * pc_z[k] * qsh_9[k];

        t_15[k] = f_0 * osh_15[k]
                  + f_3 * pc_x[k] * qsh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, osh_17, osh_18, \
                         osh_20, qsh_10, qsh_14, qsh_17, qsh_18, \
                         qsh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * qsh_10[k];

        t_17[k] = f_0 * osh_17[k]
                  + f_3 * pc_x[k] * qsh_17[k];

        t_18[k] = f_0 * osh_18[k]
                  + f_3 * pc_x[k] * qsh_18[k];

        t_19[k] = f_3 * pc_y[k] * qsh_14[k];

        t_20[k] = f_0 * osh_20[k]
                  + f_3 * pc_x[k] * qsh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, qsg0_10, qsg0_12, qsg0_13, \
                         qsg1_10, qsg1_12, qsg1_13, qsh_15, qsh_17, \
                         qsh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * qsg0_10[k]
                  - f_2 * qsg1_10[k]
                  + f_3 * pc_y[k] * qsh_15[k];

        t_22[k] = f_3 * pc_z[k] * qsh_15[k];

        t_23[k] = f_8 * qsg0_12[k]
                  - f_9 * qsg1_12[k]
                  + f_3 * pc_y[k] * qsh_17[k];

        t_24[k] = f_6 * qsg0_13[k]
                  - f_7 * qsg1_13[k]
                  + f_3 * pc_y[k] * qsh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, osi0_0, osh_0, \
                         osi1_0, qsg0_14, qsg1_14, qsh_19, qsh_20, \
                         qsh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * qsg0_14[k]
                  - f_5 * qsg1_14[k]
                  + f_3 * pc_y[k] * qsh_19[k];

        t_26[k] = f_3 * pc_y[k] * qsh_20[k];

        t_27[k] = f_1 * qsg0_14[k]
                  - f_2 * qsg1_14[k]
                  + f_3 * pc_z[k] * qsh_20[k];

        t_28[k] = pa_y[k] * osi0_0[k]
                  - f_10 * pc_y[k] * osi1_0[k];

        t_29[k] = f_11 * osh_0[k]
                  + f_3 * pc_y[k] * qsh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, osi0_3, osi0_5, osh_1, \
                         osi1_3, osi1_5, qsh_21, qsh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * qsh_21[k];

        t_31[k] = pa_y[k] * osi0_3[k]
                  + f_12 * osh_1[k]
                  - f_10 * pc_y[k] * osi1_3[k];

        t_32[k] = f_3 * pc_z[k] * qsh_22[k];

        t_33[k] = pa_y[k] * osi0_5[k]
                  - f_10 * pc_y[k] * osi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, osi0_6, osi0_9, osh_3, \
                         osh_5, osi1_6, osi1_9, qsh_24, qsh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * osi0_6[k]
                  + f_13 * osh_3[k]
                  - f_10 * pc_y[k] * osi1_6[k];

        t_35[k] = f_3 * pc_z[k] * qsh_24[k];

        t_36[k] = f_11 * osh_5[k]
                  + f_3 * pc_y[k] * qsh_26[k];

        t_37[k] = pa_y[k] * osi0_9[k]
                  - f_10 * pc_y[k] * osi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, osi0_10, osh_6, osh_9, \
                         osi1_10, qsg0_18, qsg1_18, qsh_27, qsh_28, \
                         qsh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * osi0_10[k]
                  + f_14 * osh_6[k]
                  - f_10 * pc_y[k] * osi1_10[k];

        t_39[k] = f_3 * pc_z[k] * qsh_27[k];

        t_40[k] = f_4 * qsg0_18[k]
                  - f_5 * qsg1_18[k]
                  + f_3 * pc_z[k] * qsh_28[k];

        t_41[k] = f_11 * osh_9[k]
                  + f_3 * pc_y[k] * qsh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, osi0_14, osh_36, \
                         osh_38, osi1_14, qsh_31, qsh_36, qsh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * osi0_14[k]
                  - f_10 * pc_y[k] * osi1_14[k];

        t_43[k] = f_15 * osh_36[k]
                  + f_3 * pc_x[k] * qsh_36[k];

        t_44[k] = f_3 * pc_z[k] * qsh_31[k];

        t_45[k] = f_15 * osh_38[k]
                  + f_3 * pc_x[k] * qsh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, osh_15, osh_39, osh_40, osh_41, \
                         qsg0_25, qsg1_25, qsh_36, qsh_39, qsh_40, \
                         qsh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * osh_39[k]
                  + f_3 * pc_x[k] * qsh_39[k];

        t_47[k] = f_15 * osh_40[k]
                  + f_3 * pc_x[k] * qsh_40[k];

        t_48[k] = f_15 * osh_41[k]
                  + f_3 * pc_x[k] * qsh_41[k];

        t_49[k] = f_11 * osh_15[k]
                  + f_1 * qsg0_25[k]
                  - f_2 * qsg1_25[k]
                  + f_3 * pc_y[k] * qsh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, qsg0_25, qsg0_26, qsg0_27, qsg1_25, \
                         qsg1_26, qsg1_27, qsh_36, qsh_37, qsh_38, \
                         qsh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * qsh_36[k];

        t_51[k] = f_4 * qsg0_25[k]
                  - f_5 * qsg1_25[k]
                  + f_3 * pc_z[k] * qsh_37[k];

        t_52[k] = f_6 * qsg0_26[k]
                  - f_7 * qsg1_26[k]
                  + f_3 * pc_z[k] * qsh_38[k];

        t_53[k] = f_8 * qsg0_27[k]
                  - f_9 * qsg1_27[k]
                  + f_3 * pc_z[k] * qsh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, osi0_0, osi0_27, \
                         osh_20, osi1_0, osi1_27, qsh_41, qsh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * osh_20[k]
                  + f_3 * pc_y[k] * qsh_41[k];

        t_55[k] = pa_y[k] * osi0_27[k]
                  - f_10 * pc_y[k] * osi1_27[k];

        t_56[k] = pa_z[k] * osi0_0[k]
                  - f_10 * pc_z[k] * osi1_0[k];

        t_57[k] = f_3 * pc_y[k] * qsh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, osi0_3, osi0_5, osh_0, \
                         osh_2, osi1_3, osi1_5, qsh_42, qsh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * osh_0[k]
                  + f_3 * pc_z[k] * qsh_42[k];

        t_59[k] = pa_z[k] * osi0_3[k]
                  - f_10 * pc_z[k] * osi1_3[k];

        t_60[k] = f_3 * pc_y[k] * qsh_44[k];

        t_61[k] = pa_z[k] * osi0_5[k]
                  + f_12 * osh_2[k]
                  - f_10 * pc_z[k] * osi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, osi0_6, osi0_9, osh_5, \
                         osi1_6, osi1_9, qsg0_32, qsg1_32, qsh_46, \
                         qsh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * osi0_6[k]
                  - f_10 * pc_z[k] * osi1_6[k];

        t_63[k] = f_4 * qsg0_32[k]
                  - f_5 * qsg1_32[k]
                  + f_3 * pc_y[k] * qsh_46[k];

        t_64[k] = f_3 * pc_y[k] * qsh_47[k];

        t_65[k] = pa_z[k] * osi0_9[k]
                  + f_13 * osh_5[k]
                  - f_10 * pc_z[k] * osi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, osi0_10, osi1_10, qsg0_34, \
                         qsg0_35, qsg1_34, qsg1_35, qsh_49, qsh_50, \
                         qsh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * osi0_10[k]
                  - f_10 * pc_z[k] * osi1_10[k];

        t_67[k] = f_6 * qsg0_34[k]
                  - f_7 * qsg1_34[k]
                  + f_3 * pc_y[k] * qsh_49[k];

        t_68[k] = f_4 * qsg0_35[k]
                  - f_5 * qsg1_35[k]
                  + f_3 * pc_y[k] * qsh_50[k];

        t_69[k] = f_3 * pc_y[k] * qsh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, osi0_14, osh_9, osh_57, \
                         osh_58, osh_59, osi1_14, qsh_57, qsh_58, \
                         qsh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * osi0_14[k]
                  + f_14 * osh_9[k]
                  - f_10 * pc_z[k] * osi1_14[k];

        t_71[k] = f_15 * osh_57[k]
                  + f_3 * pc_x[k] * qsh_57[k];

        t_72[k] = f_15 * osh_58[k]
                  + f_3 * pc_x[k] * qsh_58[k];

        t_73[k] = f_15 * osh_59[k]
                  + f_3 * pc_x[k] * qsh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, osi0_21, osh_60, \
                         osh_62, osi1_21, qsh_56, qsh_60, qsh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_15 * osh_60[k]
                  + f_3 * pc_x[k] * qsh_60[k];

        t_75[k] = f_3 * pc_y[k] * qsh_56[k];

        t_76[k] = f_15 * osh_62[k]
                  + f_3 * pc_x[k] * qsh_62[k];

        t_77[k] = pa_z[k] * osi0_21[k]
                  - f_10 * pc_z[k] * osi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, qsg0_41, qsg0_42, qsg0_43, qsg1_41, qsg1_42, \
                         qsg1_43, qsh_58, qsh_59, qsh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * qsg0_41[k]
                  - f_17 * qsg1_41[k]
                  + f_3 * pc_y[k] * qsh_58[k];

        t_79[k] = f_8 * qsg0_42[k]
                  - f_9 * qsg1_42[k]
                  + f_3 * pc_y[k] * qsh_59[k];

        t_80[k] = f_6 * qsg0_43[k]
                  - f_7 * qsg1_43[k]
                  + f_3 * pc_y[k] * qsh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, osh_20, osh_63, qsg0_44, \
                         qsg0_45, qsg1_44, qsg1_45, qsh_61, qsh_62, \
                         qsh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * qsg0_44[k]
                  - f_5 * qsg1_44[k]
                  + f_3 * pc_y[k] * qsh_61[k];

        t_82[k] = f_3 * pc_y[k] * qsh_62[k];

        t_83[k] = f_11 * osh_20[k]
                  + f_1 * qsg0_44[k]
                  - f_2 * qsg1_44[k]
                  + f_3 * pc_z[k] * qsh_62[k];

        t_84[k] = f_18 * osh_63[k]
                  + f_1 * qsg0_45[k]
                  - f_2 * qsg1_45[k]
                  + f_3 * pc_x[k] * qsh_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, osh_21, osh_66, qsg0_48, \
                         qsg1_48, qsh_63, qsh_64, qsh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * osh_21[k]
                  + f_3 * pc_y[k] * qsh_63[k];

        t_86[k] = f_3 * pc_z[k] * qsh_63[k];

        t_87[k] = f_18 * osh_66[k]
                  + f_8 * qsg0_48[k]
                  - f_9 * qsg1_48[k]
                  + f_3 * pc_x[k] * qsh_66[k];

        t_88[k] = f_3 * pc_z[k] * qsh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, osh_69, qsg0_45, qsg0_51, qsg1_45, \
                         qsg1_51, qsh_65, qsh_66, qsh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * qsg0_45[k]
                  - f_5 * qsg1_45[k]
                  + f_3 * pc_z[k] * qsh_65[k];

        t_90[k] = f_18 * osh_69[k]
                  + f_6 * qsg0_51[k]
                  - f_7 * qsg1_51[k]
                  + f_3 * pc_x[k] * qsh_69[k];

        t_91[k] = f_3 * pc_z[k] * qsh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, osh_26, osh_73, qsg0_47, \
                         qsg0_55, qsg1_47, qsg1_55, qsh_68, qsh_69, \
                         qsh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * osh_26[k]
                  + f_3 * pc_y[k] * qsh_68[k];

        t_93[k] = f_6 * qsg0_47[k]
                  - f_7 * qsg1_47[k]
                  + f_3 * pc_z[k] * qsh_68[k];

        t_94[k] = f_18 * osh_73[k]
                  + f_4 * qsg0_55[k]
                  - f_5 * qsg1_55[k]
                  + f_3 * pc_x[k] * qsh_73[k];

        t_95[k] = f_3 * pc_z[k] * qsh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, osh_30, osh_78, qsg0_48, \
                         qsg0_50, qsg1_48, qsg1_50, qsh_70, qsh_72, \
                         qsh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * qsg0_48[k]
                  - f_5 * qsg1_48[k]
                  + f_3 * pc_z[k] * qsh_70[k];

        t_97[k] = f_12 * osh_30[k]
                  + f_3 * pc_y[k] * qsh_72[k];

        t_98[k] = f_8 * qsg0_50[k]
                  - f_9 * qsg1_50[k]
                  + f_3 * pc_z[k] * qsh_72[k];

        t_99[k] = f_18 * osh_78[k]
                  + f_3 * pc_x[k] * qsh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, osh_80, osh_81, \
                         osh_82, osh_83, qsh_73, qsh_80, qsh_81, qsh_82, \
                         qsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * qsh_73[k];

        t_101[k] = f_18 * osh_80[k]
                   + f_3 * pc_x[k] * qsh_80[k];

        t_102[k] = f_18 * osh_81[k]
                   + f_3 * pc_x[k] * qsh_81[k];

        t_103[k] = f_18 * osh_82[k]
                   + f_3 * pc_x[k] * qsh_82[k];

        t_104[k] = f_18 * osh_83[k]
                   + f_3 * pc_x[k] * qsh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, osh_36, qsg0_55, qsg0_56, \
                         qsg1_55, qsg1_56, qsh_78, qsh_79, qsh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * osh_36[k]
                   + f_1 * qsg0_55[k]
                   - f_2 * qsg1_55[k]
                   + f_3 * pc_y[k] * qsh_78[k];

        t_106[k] = f_3 * pc_z[k] * qsh_78[k];

        t_107[k] = f_4 * qsg0_55[k]
                   - f_5 * qsg1_55[k]
                   + f_3 * pc_z[k] * qsh_79[k];

        t_108[k] = f_6 * qsg0_56[k]
                   - f_7 * qsg1_56[k]
                   + f_3 * pc_z[k] * qsh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, osi0_56, osh_41, \
                         osi1_56, qsg0_57, qsg0_59, qsg1_57, qsg1_59, qsh_81, \
                         qsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * qsg0_57[k]
                   - f_9 * qsg1_57[k]
                   + f_3 * pc_z[k] * qsh_81[k];

        t_110[k] = f_12 * osh_41[k]
                   + f_3 * pc_y[k] * qsh_83[k];

        t_111[k] = f_1 * qsg0_59[k]
                   - f_2 * qsg1_59[k]
                   + f_3 * pc_z[k] * qsh_83[k];

        t_112[k] = pa_y[k] * osi0_56[k]
                   - f_10 * pc_y[k] * osi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, osi0_31, osh_21, \
                         osh_42, osh_44, osi1_31, qsh_84, qsh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * osh_42[k]
                   + f_3 * pc_y[k] * qsh_84[k];

        t_114[k] = f_11 * osh_21[k]
                   + f_3 * pc_z[k] * qsh_84[k];

        t_115[k] = pa_z[k] * osi0_31[k]
                   - f_10 * pc_z[k] * osi1_31[k];

        t_116[k] = f_11 * osh_44[k]
                   + f_3 * pc_y[k] * qsh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, osi0_34, osi0_61, \
                         osh_24, osh_47, osi1_34, osi1_61, qsh_87, \
                         qsh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * osi0_61[k]
                   - f_10 * pc_y[k] * osi1_61[k];

        t_118[k] = pa_z[k] * osi0_34[k]
                   - f_10 * pc_z[k] * osi1_34[k];

        t_119[k] = f_11 * osh_24[k]
                   + f_3 * pc_z[k] * qsh_87[k];

        t_120[k] = f_11 * osh_47[k]
                   + f_3 * pc_y[k] * qsh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, osi0_38, osi0_65, \
                         osh_27, osi1_38, osi1_65, qsh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * osi0_65[k]
                   - f_10 * pc_y[k] * osi1_65[k];

        t_122[k] = pa_z[k] * osi0_38[k]
                   - f_10 * pc_z[k] * osi1_38[k];

        t_123[k] = f_11 * osh_27[k]
                   + f_3 * pc_z[k] * qsh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, osi0_68, osi0_70, \
                         osh_50, osh_51, osh_99, osi1_68, osi1_70, qsh_93, \
                         qsh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * osi0_68[k]
                   + f_12 * osh_50[k]
                   - f_10 * pc_y[k] * osi1_68[k];

        t_125[k] = f_11 * osh_51[k]
                   + f_3 * pc_y[k] * qsh_93[k];

        t_126[k] = pa_y[k] * osi0_70[k]
                   - f_10 * pc_y[k] * osi1_70[k];

        t_127[k] = f_18 * osh_99[k]
                   + f_3 * pc_x[k] * qsh_99[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 5.0 / q;
    const auto f_19 = 4.5 / q;

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
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_49 = buffer.data(osi0 + 49);
    const auto *osi0_83 = buffer.data(osi0 + 83);
    const auto *osi0_84 = buffer.data(osi0 + 84);
    const auto *osi0_87 = buffer.data(osi0 + 87);
    const auto *osi0_90 = buffer.data(osi0 + 90);
    const auto *osi0_94 = buffer.data(osi0 + 94);
    const auto *osi0_96 = buffer.data(osi0 + 96);
    const auto *osi0_105 = buffer.data(osi0 + 105);
    const auto *osi0_140 = buffer.data(osi0 + 140);
    const auto *osi0_143 = buffer.data(osi0 + 143);
    const auto *osi0_145 = buffer.data(osi0 + 145);
    const auto *osi0_146 = buffer.data(osi0 + 146);
    const auto *osi0_149 = buffer.data(osi0 + 149);
    const auto *osi0_150 = buffer.data(osi0 + 150);
    const auto *osi0_152 = buffer.data(osi0 + 152);
    const auto *osi0_154 = buffer.data(osi0 + 154);

    const auto *osh_36 = buffer.data(osh + 36);
    const auto *osh_42 = buffer.data(osh + 42);
    const auto *osh_59 = buffer.data(osh + 59);
    const auto *osh_60 = buffer.data(osh + 60);
    const auto *osh_61 = buffer.data(osh + 61);
    const auto *osh_62 = buffer.data(osh + 62);
    const auto *osh_63 = buffer.data(osh + 63);
    const auto *osh_66 = buffer.data(osh + 66);
    const auto *osh_68 = buffer.data(osh + 68);
    const auto *osh_69 = buffer.data(osh + 69);
    const auto *osh_70 = buffer.data(osh + 70);
    const auto *osh_72 = buffer.data(osh + 72);
    const auto *osh_78 = buffer.data(osh + 78);
    const auto *osh_83 = buffer.data(osh + 83);
    const auto *osh_84 = buffer.data(osh + 84);
    const auto *osh_86 = buffer.data(osh + 86);
    const auto *osh_87 = buffer.data(osh + 87);
    const auto *osh_89 = buffer.data(osh + 89);
    const auto *osh_90 = buffer.data(osh + 90);
    const auto *osh_93 = buffer.data(osh + 93);
    const auto *osh_99 = buffer.data(osh + 99);
    const auto *osh_100 = buffer.data(osh + 100);
    const auto *osh_101 = buffer.data(osh + 101);
    const auto *osh_102 = buffer.data(osh + 102);
    const auto *osh_103 = buffer.data(osh + 103);
    const auto *osh_104 = buffer.data(osh + 104);
    const auto *osh_105 = buffer.data(osh + 105);
    const auto *osh_106 = buffer.data(osh + 106);
    const auto *osh_107 = buffer.data(osh + 107);
    const auto *osh_108 = buffer.data(osh + 108);
    const auto *osh_110 = buffer.data(osh + 110);
    const auto *osh_111 = buffer.data(osh + 111);
    const auto *osh_113 = buffer.data(osh + 113);
    const auto *osh_114 = buffer.data(osh + 114);
    const auto *osh_119 = buffer.data(osh + 119);
    const auto *osh_120 = buffer.data(osh + 120);
    const auto *osh_121 = buffer.data(osh + 121);
    const auto *osh_122 = buffer.data(osh + 122);
    const auto *osh_123 = buffer.data(osh + 123);
    const auto *osh_125 = buffer.data(osh + 125);
    const auto *osh_126 = buffer.data(osh + 126);
    const auto *osh_129 = buffer.data(osh + 129);
    const auto *osh_132 = buffer.data(osh + 132);
    const auto *osh_136 = buffer.data(osh + 136);
    const auto *osh_141 = buffer.data(osh + 141);
    const auto *osh_143 = buffer.data(osh + 143);
    const auto *osh_144 = buffer.data(osh + 144);
    const auto *osh_145 = buffer.data(osh + 145);
    const auto *osh_146 = buffer.data(osh + 146);
    const auto *osh_152 = buffer.data(osh + 152);
    const auto *osh_156 = buffer.data(osh + 156);
    const auto *osh_161 = buffer.data(osh + 161);
    const auto *osh_162 = buffer.data(osh + 162);
    const auto *osh_163 = buffer.data(osh + 163);
    const auto *osh_164 = buffer.data(osh + 164);
    const auto *osh_165 = buffer.data(osh + 165);
    const auto *osh_166 = buffer.data(osh + 166);
    const auto *osh_167 = buffer.data(osh + 167);
    const auto *osh_183 = buffer.data(osh + 183);
    const auto *osh_184 = buffer.data(osh + 184);
    const auto *osh_185 = buffer.data(osh + 185);
    const auto *osh_186 = buffer.data(osh + 186);
    const auto *osh_187 = buffer.data(osh + 187);
    const auto *osh_188 = buffer.data(osh + 188);

    const auto *osi1_49 = buffer.data(osi1 + 49);
    const auto *osi1_83 = buffer.data(osi1 + 83);
    const auto *osi1_84 = buffer.data(osi1 + 84);
    const auto *osi1_87 = buffer.data(osi1 + 87);
    const auto *osi1_90 = buffer.data(osi1 + 90);
    const auto *osi1_94 = buffer.data(osi1 + 94);
    const auto *osi1_96 = buffer.data(osi1 + 96);
    const auto *osi1_105 = buffer.data(osi1 + 105);
    const auto *osi1_140 = buffer.data(osi1 + 140);
    const auto *osi1_143 = buffer.data(osi1 + 143);
    const auto *osi1_145 = buffer.data(osi1 + 145);
    const auto *osi1_146 = buffer.data(osi1 + 146);
    const auto *osi1_149 = buffer.data(osi1 + 149);
    const auto *osi1_150 = buffer.data(osi1 + 150);
    const auto *osi1_152 = buffer.data(osi1 + 152);
    const auto *osi1_154 = buffer.data(osi1 + 154);

    const auto *qsg0_72 = buffer.data(qsg0 + 72);
    const auto *qsg0_73 = buffer.data(qsg0 + 73);
    const auto *qsg0_74 = buffer.data(qsg0 + 74);
    const auto *qsg0_75 = buffer.data(qsg0 + 75);
    const auto *qsg0_76 = buffer.data(qsg0 + 76);
    const auto *qsg0_77 = buffer.data(qsg0 + 77);
    const auto *qsg0_78 = buffer.data(qsg0 + 78);
    const auto *qsg0_79 = buffer.data(qsg0 + 79);
    const auto *qsg0_80 = buffer.data(qsg0 + 80);
    const auto *qsg0_84 = buffer.data(qsg0 + 84);
    const auto *qsg0_85 = buffer.data(qsg0 + 85);
    const auto *qsg0_86 = buffer.data(qsg0 + 86);
    const auto *qsg0_87 = buffer.data(qsg0 + 87);
    const auto *qsg0_88 = buffer.data(qsg0 + 88);
    const auto *qsg0_89 = buffer.data(qsg0 + 89);
    const auto *qsg0_90 = buffer.data(qsg0 + 90);
    const auto *qsg0_92 = buffer.data(qsg0 + 92);
    const auto *qsg0_93 = buffer.data(qsg0 + 93);
    const auto *qsg0_95 = buffer.data(qsg0 + 95);
    const auto *qsg0_96 = buffer.data(qsg0 + 96);
    const auto *qsg0_100 = buffer.data(qsg0 + 100);
    const auto *qsg0_101 = buffer.data(qsg0 + 101);
    const auto *qsg0_102 = buffer.data(qsg0 + 102);
    const auto *qsg0_104 = buffer.data(qsg0 + 104);
    const auto *qsg0_110 = buffer.data(qsg0 + 110);
    const auto *qsg0_114 = buffer.data(qsg0 + 114);
    const auto *qsg0_117 = buffer.data(qsg0 + 117);
    const auto *qsg0_118 = buffer.data(qsg0 + 118);
    const auto *qsg0_119 = buffer.data(qsg0 + 119);
    const auto *qsg0_130 = buffer.data(qsg0 + 130);
    const auto *qsg0_132 = buffer.data(qsg0 + 132);

    const auto *qsg1_72 = buffer.data(qsg1 + 72);
    const auto *qsg1_73 = buffer.data(qsg1 + 73);
    const auto *qsg1_74 = buffer.data(qsg1 + 74);
    const auto *qsg1_75 = buffer.data(qsg1 + 75);
    const auto *qsg1_76 = buffer.data(qsg1 + 76);
    const auto *qsg1_77 = buffer.data(qsg1 + 77);
    const auto *qsg1_78 = buffer.data(qsg1 + 78);
    const auto *qsg1_79 = buffer.data(qsg1 + 79);
    const auto *qsg1_80 = buffer.data(qsg1 + 80);
    const auto *qsg1_84 = buffer.data(qsg1 + 84);
    const auto *qsg1_85 = buffer.data(qsg1 + 85);
    const auto *qsg1_86 = buffer.data(qsg1 + 86);
    const auto *qsg1_87 = buffer.data(qsg1 + 87);
    const auto *qsg1_88 = buffer.data(qsg1 + 88);
    const auto *qsg1_89 = buffer.data(qsg1 + 89);
    const auto *qsg1_90 = buffer.data(qsg1 + 90);
    const auto *qsg1_92 = buffer.data(qsg1 + 92);
    const auto *qsg1_93 = buffer.data(qsg1 + 93);
    const auto *qsg1_95 = buffer.data(qsg1 + 95);
    const auto *qsg1_96 = buffer.data(qsg1 + 96);
    const auto *qsg1_100 = buffer.data(qsg1 + 100);
    const auto *qsg1_101 = buffer.data(qsg1 + 101);
    const auto *qsg1_102 = buffer.data(qsg1 + 102);
    const auto *qsg1_104 = buffer.data(qsg1 + 104);
    const auto *qsg1_110 = buffer.data(qsg1 + 110);
    const auto *qsg1_114 = buffer.data(qsg1 + 114);
    const auto *qsg1_117 = buffer.data(qsg1 + 117);
    const auto *qsg1_118 = buffer.data(qsg1 + 118);
    const auto *qsg1_119 = buffer.data(qsg1 + 119);
    const auto *qsg1_130 = buffer.data(qsg1 + 130);
    const auto *qsg1_132 = buffer.data(qsg1 + 132);

    const auto *qsh_99 = buffer.data(qsh + 99);
    const auto *qsh_100 = buffer.data(qsh + 100);
    const auto *qsh_101 = buffer.data(qsh + 101);
    const auto *qsh_102 = buffer.data(qsh + 102);
    const auto *qsh_103 = buffer.data(qsh + 103);
    const auto *qsh_104 = buffer.data(qsh + 104);
    const auto *qsh_105 = buffer.data(qsh + 105);
    const auto *qsh_106 = buffer.data(qsh + 106);
    const auto *qsh_107 = buffer.data(qsh + 107);
    const auto *qsh_108 = buffer.data(qsh + 108);
    const auto *qsh_109 = buffer.data(qsh + 109);
    const auto *qsh_110 = buffer.data(qsh + 110);
    const auto *qsh_111 = buffer.data(qsh + 111);
    const auto *qsh_112 = buffer.data(qsh + 112);
    const auto *qsh_113 = buffer.data(qsh + 113);
    const auto *qsh_114 = buffer.data(qsh + 114);
    const auto *qsh_119 = buffer.data(qsh + 119);
    const auto *qsh_120 = buffer.data(qsh + 120);
    const auto *qsh_121 = buffer.data(qsh + 121);
    const auto *qsh_122 = buffer.data(qsh + 122);
    const auto *qsh_123 = buffer.data(qsh + 123);
    const auto *qsh_124 = buffer.data(qsh + 124);
    const auto *qsh_125 = buffer.data(qsh + 125);
    const auto *qsh_126 = buffer.data(qsh + 126);
    const auto *qsh_127 = buffer.data(qsh + 127);
    const auto *qsh_128 = buffer.data(qsh + 128);
    const auto *qsh_129 = buffer.data(qsh + 129);
    const auto *qsh_131 = buffer.data(qsh + 131);
    const auto *qsh_132 = buffer.data(qsh + 132);
    const auto *qsh_133 = buffer.data(qsh + 133);
    const auto *qsh_135 = buffer.data(qsh + 135);
    const auto *qsh_136 = buffer.data(qsh + 136);
    const auto *qsh_141 = buffer.data(qsh + 141);
    const auto *qsh_142 = buffer.data(qsh + 142);
    const auto *qsh_143 = buffer.data(qsh + 143);
    const auto *qsh_144 = buffer.data(qsh + 144);
    const auto *qsh_145 = buffer.data(qsh + 145);
    const auto *qsh_146 = buffer.data(qsh + 146);
    const auto *qsh_147 = buffer.data(qsh + 147);
    const auto *qsh_149 = buffer.data(qsh + 149);
    const auto *qsh_150 = buffer.data(qsh + 150);
    const auto *qsh_152 = buffer.data(qsh + 152);
    const auto *qsh_153 = buffer.data(qsh + 153);
    const auto *qsh_156 = buffer.data(qsh + 156);
    const auto *qsh_161 = buffer.data(qsh + 161);
    const auto *qsh_162 = buffer.data(qsh + 162);
    const auto *qsh_163 = buffer.data(qsh + 163);
    const auto *qsh_164 = buffer.data(qsh + 164);
    const auto *qsh_165 = buffer.data(qsh + 165);
    const auto *qsh_166 = buffer.data(qsh + 166);
    const auto *qsh_167 = buffer.data(qsh + 167);
    const auto *qsh_168 = buffer.data(qsh + 168);
    const auto *qsh_170 = buffer.data(qsh + 170);
    const auto *qsh_171 = buffer.data(qsh + 171);
    const auto *qsh_173 = buffer.data(qsh + 173);
    const auto *qsh_174 = buffer.data(qsh + 174);
    const auto *qsh_177 = buffer.data(qsh + 177);
    const auto *qsh_183 = buffer.data(qsh + 183);
    const auto *qsh_184 = buffer.data(qsh + 184);
    const auto *qsh_185 = buffer.data(qsh + 185);
    const auto *qsh_186 = buffer.data(qsh + 186);
    const auto *qsh_187 = buffer.data(qsh + 187);
    const auto *qsh_188 = buffer.data(qsh + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, osh_100, osh_101, osh_102, \
                         osh_103, osh_104, qsh_100, qsh_101, qsh_102, qsh_103, \
                         qsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_18 * osh_100[k]
                   + f_3 * pc_x[k] * qsh_100[k];

        t_129[k] = f_18 * osh_101[k]
                   + f_3 * pc_x[k] * qsh_101[k];

        t_130[k] = f_18 * osh_102[k]
                   + f_3 * pc_x[k] * qsh_102[k];

        t_131[k] = f_18 * osh_103[k]
                   + f_3 * pc_x[k] * qsh_103[k];

        t_132[k] = f_18 * osh_104[k]
                   + f_3 * pc_x[k] * qsh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, osi0_49, osh_36, osh_59, \
                         osi1_49, qsg0_72, qsg1_72, qsh_99, qsh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * osi0_49[k]
                   - f_10 * pc_z[k] * osi1_49[k];

        t_134[k] = f_11 * osh_36[k]
                   + f_3 * pc_z[k] * qsh_99[k];

        t_135[k] = f_11 * osh_59[k]
                   + f_8 * qsg0_72[k]
                   - f_9 * qsg1_72[k]
                   + f_3 * pc_y[k] * qsh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, osh_60, osh_61, osh_62, qsg0_73, qsg0_74, \
                         qsg1_73, qsg1_74, qsh_102, qsh_103, qsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * osh_60[k]
                   + f_6 * qsg0_73[k]
                   - f_7 * qsg1_73[k]
                   + f_3 * pc_y[k] * qsh_102[k];

        t_137[k] = f_11 * osh_61[k]
                   + f_4 * qsg0_74[k]
                   - f_5 * qsg1_74[k]
                   + f_3 * pc_y[k] * qsh_103[k];

        t_138[k] = f_11 * osh_62[k]
                   + f_3 * pc_y[k] * qsh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, osi0_83, osh_42, \
                         osh_105, osi1_83, qsg0_75, qsg1_75, qsh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * osi0_83[k]
                   - f_10 * pc_y[k] * osi1_83[k];

        t_140[k] = f_18 * osh_105[k]
                   + f_1 * qsg0_75[k]
                   - f_2 * qsg1_75[k]
                   + f_3 * pc_x[k] * qsh_105[k];

        t_141[k] = f_3 * pc_y[k] * qsh_105[k];

        t_142[k] = f_12 * osh_42[k]
                   + f_3 * pc_z[k] * qsh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, osh_110, qsg0_75, qsg0_80, qsg1_75, \
                         qsg1_80, qsh_106, qsh_107, qsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * qsg0_75[k]
                   - f_5 * qsg1_75[k]
                   + f_3 * pc_y[k] * qsh_106[k];

        t_144[k] = f_3 * pc_y[k] * qsh_107[k];

        t_145[k] = f_18 * osh_110[k]
                   + f_8 * qsg0_80[k]
                   - f_9 * qsg1_80[k]
                   + f_3 * pc_x[k] * qsh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, qsg0_76, qsg0_77, qsg1_76, qsg1_77, \
                         qsh_108, qsh_109, qsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * qsg0_76[k]
                   - f_7 * qsg1_76[k]
                   + f_3 * pc_y[k] * qsh_108[k];

        t_147[k] = f_4 * qsg0_77[k]
                   - f_5 * qsg1_77[k]
                   + f_3 * pc_y[k] * qsh_109[k];

        t_148[k] = f_3 * pc_y[k] * qsh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, osh_114, qsg0_78, qsg0_79, qsg0_84, \
                         qsg1_78, qsg1_79, qsg1_84, qsh_111, qsh_112, \
                         qsh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_18 * osh_114[k]
                   + f_6 * qsg0_84[k]
                   - f_7 * qsg1_84[k]
                   + f_3 * pc_x[k] * qsh_114[k];

        t_150[k] = f_8 * qsg0_78[k]
                   - f_9 * qsg1_78[k]
                   + f_3 * pc_y[k] * qsh_111[k];

        t_151[k] = f_6 * qsg0_79[k]
                   - f_7 * qsg1_79[k]
                   + f_3 * pc_y[k] * qsh_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, osh_119, osh_120, qsg0_80, \
                         qsg0_89, qsg1_80, qsg1_89, qsh_113, qsh_114, qsh_119, \
                         qsh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * qsg0_80[k]
                   - f_5 * qsg1_80[k]
                   + f_3 * pc_y[k] * qsh_113[k];

        t_153[k] = f_3 * pc_y[k] * qsh_114[k];

        t_154[k] = f_18 * osh_119[k]
                   + f_4 * qsg0_89[k]
                   - f_5 * qsg1_89[k]
                   + f_3 * pc_x[k] * qsh_119[k];

        t_155[k] = f_18 * osh_120[k]
                   + f_3 * pc_x[k] * qsh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, osh_121, osh_122, \
                         osh_123, osh_125, qsh_119, qsh_121, qsh_122, qsh_123, \
                         qsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_18 * osh_121[k]
                   + f_3 * pc_x[k] * qsh_121[k];

        t_157[k] = f_18 * osh_122[k]
                   + f_3 * pc_x[k] * qsh_122[k];

        t_158[k] = f_18 * osh_123[k]
                   + f_3 * pc_x[k] * qsh_123[k];

        t_159[k] = f_3 * pc_y[k] * qsh_119[k];

        t_160[k] = f_18 * osh_125[k]
                   + f_3 * pc_x[k] * qsh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, qsg0_85, qsg0_86, qsg0_87, qsg1_85, \
                         qsg1_86, qsg1_87, qsh_120, qsh_121, qsh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * qsg0_85[k]
                   - f_2 * qsg1_85[k]
                   + f_3 * pc_y[k] * qsh_120[k];

        t_162[k] = f_16 * qsg0_86[k]
                   - f_17 * qsg1_86[k]
                   + f_3 * pc_y[k] * qsh_121[k];

        t_163[k] = f_8 * qsg0_87[k]
                   - f_9 * qsg1_87[k]
                   + f_3 * pc_y[k] * qsh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, osh_62, qsg0_88, qsg0_89, \
                         qsg1_88, qsg1_89, qsh_123, qsh_124, qsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * qsg0_88[k]
                   - f_7 * qsg1_88[k]
                   + f_3 * pc_y[k] * qsh_123[k];

        t_165[k] = f_4 * qsg0_89[k]
                   - f_5 * qsg1_89[k]
                   + f_3 * pc_y[k] * qsh_124[k];

        t_166[k] = f_3 * pc_y[k] * qsh_125[k];

        t_167[k] = f_12 * osh_62[k]
                   + f_1 * qsg0_89[k]
                   - f_2 * qsg1_89[k]
                   + f_3 * pc_z[k] * qsh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, osh_63, osh_126, \
                         osh_129, qsg0_90, qsg0_93, qsg1_90, qsg1_93, qsh_126, \
                         qsh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * osh_126[k]
                   + f_1 * qsg0_90[k]
                   - f_2 * qsg1_90[k]
                   + f_3 * pc_x[k] * qsh_126[k];

        t_169[k] = f_13 * osh_63[k]
                   + f_3 * pc_y[k] * qsh_126[k];

        t_170[k] = f_3 * pc_z[k] * qsh_126[k];

        t_171[k] = f_19 * osh_129[k]
                   + f_8 * qsg0_93[k]
                   - f_9 * qsg1_93[k]
                   + f_3 * pc_x[k] * qsh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_x, pc_z, osh_132, qsg0_90, qsg0_96, \
                         qsg1_90, qsg1_96, qsh_127, qsh_128, qsh_129, \
                         qsh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * qsh_127[k];

        t_173[k] = f_4 * qsg0_90[k]
                   - f_5 * qsg1_90[k]
                   + f_3 * pc_z[k] * qsh_128[k];

        t_174[k] = f_19 * osh_132[k]
                   + f_6 * qsg0_96[k]
                   - f_7 * qsg1_96[k]
                   + f_3 * pc_x[k] * qsh_132[k];

        t_175[k] = f_3 * pc_z[k] * qsh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, osh_68, osh_136, \
                         qsg0_92, qsg0_100, qsg1_92, qsg1_100, qsh_131, qsh_132, \
                         qsh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * osh_68[k]
                   + f_3 * pc_y[k] * qsh_131[k];

        t_177[k] = f_6 * qsg0_92[k]
                   - f_7 * qsg1_92[k]
                   + f_3 * pc_z[k] * qsh_131[k];

        t_178[k] = f_19 * osh_136[k]
                   + f_4 * qsg0_100[k]
                   - f_5 * qsg1_100[k]
                   + f_3 * pc_x[k] * qsh_136[k];

        t_179[k] = f_3 * pc_z[k] * qsh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, osh_72, osh_141, \
                         qsg0_93, qsg0_95, qsg1_93, qsg1_95, qsh_133, qsh_135, \
                         qsh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * qsg0_93[k]
                   - f_5 * qsg1_93[k]
                   + f_3 * pc_z[k] * qsh_133[k];

        t_181[k] = f_13 * osh_72[k]
                   + f_3 * pc_y[k] * qsh_135[k];

        t_182[k] = f_8 * qsg0_95[k]
                   - f_9 * qsg1_95[k]
                   + f_3 * pc_z[k] * qsh_135[k];

        t_183[k] = f_19 * osh_141[k]
                   + f_3 * pc_x[k] * qsh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, osh_143, osh_144, \
                         osh_145, osh_146, qsh_136, qsh_143, qsh_144, qsh_145, \
                         qsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * qsh_136[k];

        t_185[k] = f_19 * osh_143[k]
                   + f_3 * pc_x[k] * qsh_143[k];

        t_186[k] = f_19 * osh_144[k]
                   + f_3 * pc_x[k] * qsh_144[k];

        t_187[k] = f_19 * osh_145[k]
                   + f_3 * pc_x[k] * qsh_145[k];

        t_188[k] = f_19 * osh_146[k]
                   + f_3 * pc_x[k] * qsh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_y, pc_z, osh_78, qsg0_100, qsg0_101, \
                         qsg1_100, qsg1_101, qsh_141, qsh_142, \
                         qsh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * osh_78[k]
                   + f_1 * qsg0_100[k]
                   - f_2 * qsg1_100[k]
                   + f_3 * pc_y[k] * qsh_141[k];

        t_190[k] = f_3 * pc_z[k] * qsh_141[k];

        t_191[k] = f_4 * qsg0_100[k]
                   - f_5 * qsg1_100[k]
                   + f_3 * pc_z[k] * qsh_142[k];

        t_192[k] = f_6 * qsg0_101[k]
                   - f_7 * qsg1_101[k]
                   + f_3 * pc_z[k] * qsh_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pc_y, pc_z, osi0_84, osh_83, \
                         osi1_84, qsg0_102, qsg0_104, qsg1_102, qsg1_104, qsh_144, \
                         qsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * qsg0_102[k]
                   - f_9 * qsg1_102[k]
                   + f_3 * pc_z[k] * qsh_144[k];

        t_194[k] = f_13 * osh_83[k]
                   + f_3 * pc_y[k] * qsh_146[k];

        t_195[k] = f_1 * qsg0_104[k]
                   - f_2 * qsg1_104[k]
                   + f_3 * pc_z[k] * qsh_146[k];

        t_196[k] = pa_z[k] * osi0_84[k]
                   - f_10 * pc_z[k] * osi1_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pc_y, pc_z, osi0_87, osh_63, \
                         osh_84, osh_86, osi1_87, qsh_147, qsh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * osh_84[k]
                   + f_3 * pc_y[k] * qsh_147[k];

        t_198[k] = f_11 * osh_63[k]
                   + f_3 * pc_z[k] * qsh_147[k];

        t_199[k] = pa_z[k] * osi0_87[k]
                   - f_10 * pc_z[k] * osi1_87[k];

        t_200[k] = f_12 * osh_86[k]
                   + f_3 * pc_y[k] * qsh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pc_x, pc_z, osi0_90, osh_66, osh_152, \
                         osi1_90, qsg0_110, qsg1_110, qsh_150, \
                         qsh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_19 * osh_152[k]
                   + f_8 * qsg0_110[k]
                   - f_9 * qsg1_110[k]
                   + f_3 * pc_x[k] * qsh_152[k];

        t_202[k] = pa_z[k] * osi0_90[k]
                   - f_10 * pc_z[k] * osi1_90[k];

        t_203[k] = f_11 * osh_66[k]
                   + f_3 * pc_z[k] * qsh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_y, pc_z, osi0_94, osh_89, \
                         osh_156, osi1_94, qsg0_114, qsg1_114, qsh_152, \
                         qsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * osh_89[k]
                   + f_3 * pc_y[k] * qsh_152[k];

        t_205[k] = f_19 * osh_156[k]
                   + f_6 * qsg0_114[k]
                   - f_7 * qsg1_114[k]
                   + f_3 * pc_x[k] * qsh_156[k];

        t_206[k] = pa_z[k] * osi0_94[k]
                   - f_10 * pc_z[k] * osi1_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pc_y, pc_z, osi0_96, osh_69, osh_70, \
                         osh_93, osi1_96, qsh_153, qsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_11 * osh_69[k]
                   + f_3 * pc_z[k] * qsh_153[k];

        t_208[k] = pa_z[k] * osi0_96[k]
                   + f_12 * osh_70[k]
                   - f_10 * pc_z[k] * osi1_96[k];

        t_209[k] = f_12 * osh_93[k]
                   + f_3 * pc_y[k] * qsh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, osh_161, osh_162, osh_163, osh_164, \
                         qsg0_119, qsg1_119, qsh_161, qsh_162, qsh_163, \
                         qsh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_19 * osh_161[k]
                   + f_4 * qsg0_119[k]
                   - f_5 * qsg1_119[k]
                   + f_3 * pc_x[k] * qsh_161[k];

        t_211[k] = f_19 * osh_162[k]
                   + f_3 * pc_x[k] * qsh_162[k];

        t_212[k] = f_19 * osh_163[k]
                   + f_3 * pc_x[k] * qsh_163[k];

        t_213[k] = f_19 * osh_164[k]
                   + f_3 * pc_x[k] * qsh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, osi0_105, osh_165, \
                         osh_166, osh_167, osi1_105, qsh_165, qsh_166, \
                         qsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_19 * osh_165[k]
                   + f_3 * pc_x[k] * qsh_165[k];

        t_215[k] = f_19 * osh_166[k]
                   + f_3 * pc_x[k] * qsh_166[k];

        t_216[k] = f_19 * osh_167[k]
                   + f_3 * pc_x[k] * qsh_167[k];

        t_217[k] = pa_z[k] * osi0_105[k]
                   - f_10 * pc_z[k] * osi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, osh_78, osh_101, osh_102, qsg0_117, \
                         qsg0_118, qsg1_117, qsg1_118, qsh_162, qsh_164, \
                         qsh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * osh_78[k]
                   + f_3 * pc_z[k] * qsh_162[k];

        t_219[k] = f_12 * osh_101[k]
                   + f_8 * qsg0_117[k]
                   - f_9 * qsg1_117[k]
                   + f_3 * pc_y[k] * qsh_164[k];

        t_220[k] = f_12 * osh_102[k]
                   + f_6 * qsg0_118[k]
                   - f_7 * qsg1_118[k]
                   + f_3 * pc_y[k] * qsh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pc_y, pc_z, osi0_140, osh_83, \
                         osh_103, osh_104, osi1_140, qsg0_119, qsg1_119, qsh_166, \
                         qsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * osh_103[k]
                   + f_4 * qsg0_119[k]
                   - f_5 * qsg1_119[k]
                   + f_3 * pc_y[k] * qsh_166[k];

        t_222[k] = f_12 * osh_104[k]
                   + f_3 * pc_y[k] * qsh_167[k];

        t_223[k] = f_11 * osh_83[k]
                   + f_1 * qsg0_119[k]
                   - f_2 * qsg1_119[k]
                   + f_3 * pc_z[k] * qsh_167[k];

        t_224[k] = pa_y[k] * osi0_140[k]
                   - f_10 * pc_y[k] * osi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, osi0_143, osh_84, \
                         osh_105, osh_106, osh_107, osi1_143, qsh_168, \
                         qsh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * osh_105[k]
                   + f_3 * pc_y[k] * qsh_168[k];

        t_226[k] = f_12 * osh_84[k]
                   + f_3 * pc_z[k] * qsh_168[k];

        t_227[k] = pa_y[k] * osi0_143[k]
                   + f_12 * osh_106[k]
                   - f_10 * pc_y[k] * osi1_143[k];

        t_228[k] = f_11 * osh_107[k]
                   + f_3 * pc_y[k] * qsh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_y, pc_z, osi0_145, osi0_146, \
                         osh_87, osh_108, osh_110, osi1_145, osi1_146, qsh_171, \
                         qsh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * osi0_145[k]
                   - f_10 * pc_y[k] * osi1_145[k];

        t_230[k] = pa_y[k] * osi0_146[k]
                   + f_13 * osh_108[k]
                   - f_10 * pc_y[k] * osi1_146[k];

        t_231[k] = f_12 * osh_87[k]
                   + f_3 * pc_z[k] * qsh_171[k];

        t_232[k] = f_11 * osh_110[k]
                   + f_3 * pc_y[k] * qsh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_y, pc_z, osi0_149, osi0_150, osh_90, \
                         osh_111, osi1_149, osi1_150, qsh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * osi0_149[k]
                   - f_10 * pc_y[k] * osi1_149[k];

        t_234[k] = pa_y[k] * osi0_150[k]
                   + f_14 * osh_111[k]
                   - f_10 * pc_y[k] * osi1_150[k];

        t_235[k] = f_12 * osh_90[k]
                   + f_3 * pc_z[k] * qsh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pc_x, pc_y, osi0_152, osi0_154, \
                         osh_113, osh_114, osh_183, osi1_152, osi1_154, qsh_177, \
                         qsh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * osi0_152[k]
                   + f_12 * osh_113[k]
                   - f_10 * pc_y[k] * osi1_152[k];

        t_237[k] = f_11 * osh_114[k]
                   + f_3 * pc_y[k] * qsh_177[k];

        t_238[k] = pa_y[k] * osi0_154[k]
                   - f_10 * pc_y[k] * osi1_154[k];

        t_239[k] = f_19 * osh_183[k]
                   + f_3 * pc_x[k] * qsh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, osh_184, osh_185, osh_186, \
                         osh_187, osh_188, qsh_184, qsh_185, qsh_186, qsh_187, \
                         qsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_19 * osh_184[k]
                   + f_3 * pc_x[k] * qsh_184[k];

        t_241[k] = f_19 * osh_185[k]
                   + f_3 * pc_x[k] * qsh_185[k];

        t_242[k] = f_19 * osh_186[k]
                   + f_3 * pc_x[k] * qsh_186[k];

        t_243[k] = f_19 * osh_187[k]
                   + f_3 * pc_x[k] * qsh_187[k];

        t_244[k] = f_19 * osh_188[k]
                   + f_3 * pc_x[k] * qsh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, osh_99, osh_120, osh_122, qsg0_130, \
                         qsg0_132, qsg1_130, qsg1_132, qsh_183, \
                         qsh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * osh_120[k]
                   + f_1 * qsg0_130[k]
                   - f_2 * qsg1_130[k]
                   + f_3 * pc_y[k] * qsh_183[k];

        t_246[k] = f_12 * osh_99[k]
                   + f_3 * pc_z[k] * qsh_183[k];

        t_247[k] = f_11 * osh_122[k]
                   + f_8 * qsg0_132[k]
                   - f_9 * qsg1_132[k]
                   + f_3 * pc_y[k] * qsh_185[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 4.5 / q;
    const auto f_20 = 4.0 / q;

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
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_167 = buffer.data(osi0 + 167);
    const auto *osi0_168 = buffer.data(osi0 + 168);
    const auto *osi0_171 = buffer.data(osi0 + 171);
    const auto *osi0_174 = buffer.data(osi0 + 174);
    const auto *osi0_178 = buffer.data(osi0 + 178);
    const auto *osi0_180 = buffer.data(osi0 + 180);
    const auto *osi0_189 = buffer.data(osi0 + 189);

    const auto *osh_105 = buffer.data(osh + 105);
    const auto *osh_123 = buffer.data(osh + 123);
    const auto *osh_124 = buffer.data(osh + 124);
    const auto *osh_125 = buffer.data(osh + 125);
    const auto *osh_126 = buffer.data(osh + 126);
    const auto *osh_129 = buffer.data(osh + 129);
    const auto *osh_131 = buffer.data(osh + 131);
    const auto *osh_132 = buffer.data(osh + 132);
    const auto *osh_133 = buffer.data(osh + 133);
    const auto *osh_135 = buffer.data(osh + 135);
    const auto *osh_141 = buffer.data(osh + 141);
    const auto *osh_146 = buffer.data(osh + 146);
    const auto *osh_147 = buffer.data(osh + 147);
    const auto *osh_149 = buffer.data(osh + 149);
    const auto *osh_150 = buffer.data(osh + 150);
    const auto *osh_152 = buffer.data(osh + 152);
    const auto *osh_153 = buffer.data(osh + 153);
    const auto *osh_156 = buffer.data(osh + 156);
    const auto *osh_162 = buffer.data(osh + 162);
    const auto *osh_164 = buffer.data(osh + 164);
    const auto *osh_165 = buffer.data(osh + 165);
    const auto *osh_166 = buffer.data(osh + 166);
    const auto *osh_167 = buffer.data(osh + 167);
    const auto *osh_168 = buffer.data(osh + 168);
    const auto *osh_170 = buffer.data(osh + 170);
    const auto *osh_173 = buffer.data(osh + 173);
    const auto *osh_177 = buffer.data(osh + 177);
    const auto *osh_183 = buffer.data(osh + 183);
    const auto *osh_185 = buffer.data(osh + 185);
    const auto *osh_186 = buffer.data(osh + 186);
    const auto *osh_187 = buffer.data(osh + 187);
    const auto *osh_189 = buffer.data(osh + 189);
    const auto *osh_194 = buffer.data(osh + 194);
    const auto *osh_198 = buffer.data(osh + 198);
    const auto *osh_203 = buffer.data(osh + 203);
    const auto *osh_204 = buffer.data(osh + 204);
    const auto *osh_205 = buffer.data(osh + 205);
    const auto *osh_206 = buffer.data(osh + 206);
    const auto *osh_207 = buffer.data(osh + 207);
    const auto *osh_209 = buffer.data(osh + 209);
    const auto *osh_210 = buffer.data(osh + 210);
    const auto *osh_213 = buffer.data(osh + 213);
    const auto *osh_216 = buffer.data(osh + 216);
    const auto *osh_220 = buffer.data(osh + 220);
    const auto *osh_225 = buffer.data(osh + 225);
    const auto *osh_227 = buffer.data(osh + 227);
    const auto *osh_228 = buffer.data(osh + 228);
    const auto *osh_229 = buffer.data(osh + 229);
    const auto *osh_230 = buffer.data(osh + 230);
    const auto *osh_236 = buffer.data(osh + 236);
    const auto *osh_240 = buffer.data(osh + 240);
    const auto *osh_245 = buffer.data(osh + 245);
    const auto *osh_246 = buffer.data(osh + 246);
    const auto *osh_247 = buffer.data(osh + 247);
    const auto *osh_248 = buffer.data(osh + 248);
    const auto *osh_249 = buffer.data(osh + 249);
    const auto *osh_250 = buffer.data(osh + 250);
    const auto *osh_251 = buffer.data(osh + 251);
    const auto *osh_252 = buffer.data(osh + 252);
    const auto *osh_255 = buffer.data(osh + 255);
    const auto *osh_257 = buffer.data(osh + 257);
    const auto *osh_258 = buffer.data(osh + 258);
    const auto *osh_261 = buffer.data(osh + 261);
    const auto *osh_262 = buffer.data(osh + 262);
    const auto *osh_264 = buffer.data(osh + 264);
    const auto *osh_266 = buffer.data(osh + 266);
    const auto *osh_267 = buffer.data(osh + 267);
    const auto *osh_268 = buffer.data(osh + 268);
    const auto *osh_269 = buffer.data(osh + 269);
    const auto *osh_270 = buffer.data(osh + 270);
    const auto *osh_271 = buffer.data(osh + 271);
    const auto *osh_272 = buffer.data(osh + 272);

    const auto *osi1_167 = buffer.data(osi1 + 167);
    const auto *osi1_168 = buffer.data(osi1 + 168);
    const auto *osi1_171 = buffer.data(osi1 + 171);
    const auto *osi1_174 = buffer.data(osi1 + 174);
    const auto *osi1_178 = buffer.data(osi1 + 178);
    const auto *osi1_180 = buffer.data(osi1 + 180);
    const auto *osi1_189 = buffer.data(osi1 + 189);

    const auto *qsg0_133 = buffer.data(qsg0 + 133);
    const auto *qsg0_134 = buffer.data(qsg0 + 134);
    const auto *qsg0_135 = buffer.data(qsg0 + 135);
    const auto *qsg0_136 = buffer.data(qsg0 + 136);
    const auto *qsg0_137 = buffer.data(qsg0 + 137);
    const auto *qsg0_138 = buffer.data(qsg0 + 138);
    const auto *qsg0_139 = buffer.data(qsg0 + 139);
    const auto *qsg0_140 = buffer.data(qsg0 + 140);
    const auto *qsg0_144 = buffer.data(qsg0 + 144);
    const auto *qsg0_145 = buffer.data(qsg0 + 145);
    const auto *qsg0_146 = buffer.data(qsg0 + 146);
    const auto *qsg0_147 = buffer.data(qsg0 + 147);
    const auto *qsg0_148 = buffer.data(qsg0 + 148);
    const auto *qsg0_149 = buffer.data(qsg0 + 149);
    const auto *qsg0_150 = buffer.data(qsg0 + 150);
    const auto *qsg0_152 = buffer.data(qsg0 + 152);
    const auto *qsg0_153 = buffer.data(qsg0 + 153);
    const auto *qsg0_155 = buffer.data(qsg0 + 155);
    const auto *qsg0_156 = buffer.data(qsg0 + 156);
    const auto *qsg0_160 = buffer.data(qsg0 + 160);
    const auto *qsg0_161 = buffer.data(qsg0 + 161);
    const auto *qsg0_162 = buffer.data(qsg0 + 162);
    const auto *qsg0_164 = buffer.data(qsg0 + 164);
    const auto *qsg0_170 = buffer.data(qsg0 + 170);
    const auto *qsg0_174 = buffer.data(qsg0 + 174);
    const auto *qsg0_177 = buffer.data(qsg0 + 177);
    const auto *qsg0_178 = buffer.data(qsg0 + 178);
    const auto *qsg0_179 = buffer.data(qsg0 + 179);
    const auto *qsg0_180 = buffer.data(qsg0 + 180);
    const auto *qsg0_183 = buffer.data(qsg0 + 183);
    const auto *qsg0_185 = buffer.data(qsg0 + 185);
    const auto *qsg0_186 = buffer.data(qsg0 + 186);
    const auto *qsg0_189 = buffer.data(qsg0 + 189);
    const auto *qsg0_190 = buffer.data(qsg0 + 190);
    const auto *qsg0_192 = buffer.data(qsg0 + 192);
    const auto *qsg0_193 = buffer.data(qsg0 + 193);
    const auto *qsg0_194 = buffer.data(qsg0 + 194);

    const auto *qsg1_133 = buffer.data(qsg1 + 133);
    const auto *qsg1_134 = buffer.data(qsg1 + 134);
    const auto *qsg1_135 = buffer.data(qsg1 + 135);
    const auto *qsg1_136 = buffer.data(qsg1 + 136);
    const auto *qsg1_137 = buffer.data(qsg1 + 137);
    const auto *qsg1_138 = buffer.data(qsg1 + 138);
    const auto *qsg1_139 = buffer.data(qsg1 + 139);
    const auto *qsg1_140 = buffer.data(qsg1 + 140);
    const auto *qsg1_144 = buffer.data(qsg1 + 144);
    const auto *qsg1_145 = buffer.data(qsg1 + 145);
    const auto *qsg1_146 = buffer.data(qsg1 + 146);
    const auto *qsg1_147 = buffer.data(qsg1 + 147);
    const auto *qsg1_148 = buffer.data(qsg1 + 148);
    const auto *qsg1_149 = buffer.data(qsg1 + 149);
    const auto *qsg1_150 = buffer.data(qsg1 + 150);
    const auto *qsg1_152 = buffer.data(qsg1 + 152);
    const auto *qsg1_153 = buffer.data(qsg1 + 153);
    const auto *qsg1_155 = buffer.data(qsg1 + 155);
    const auto *qsg1_156 = buffer.data(qsg1 + 156);
    const auto *qsg1_160 = buffer.data(qsg1 + 160);
    const auto *qsg1_161 = buffer.data(qsg1 + 161);
    const auto *qsg1_162 = buffer.data(qsg1 + 162);
    const auto *qsg1_164 = buffer.data(qsg1 + 164);
    const auto *qsg1_170 = buffer.data(qsg1 + 170);
    const auto *qsg1_174 = buffer.data(qsg1 + 174);
    const auto *qsg1_177 = buffer.data(qsg1 + 177);
    const auto *qsg1_178 = buffer.data(qsg1 + 178);
    const auto *qsg1_179 = buffer.data(qsg1 + 179);
    const auto *qsg1_180 = buffer.data(qsg1 + 180);
    const auto *qsg1_183 = buffer.data(qsg1 + 183);
    const auto *qsg1_185 = buffer.data(qsg1 + 185);
    const auto *qsg1_186 = buffer.data(qsg1 + 186);
    const auto *qsg1_189 = buffer.data(qsg1 + 189);
    const auto *qsg1_190 = buffer.data(qsg1 + 190);
    const auto *qsg1_192 = buffer.data(qsg1 + 192);
    const auto *qsg1_193 = buffer.data(qsg1 + 193);
    const auto *qsg1_194 = buffer.data(qsg1 + 194);

    const auto *qsh_186 = buffer.data(qsh + 186);
    const auto *qsh_187 = buffer.data(qsh + 187);
    const auto *qsh_188 = buffer.data(qsh + 188);
    const auto *qsh_189 = buffer.data(qsh + 189);
    const auto *qsh_190 = buffer.data(qsh + 190);
    const auto *qsh_191 = buffer.data(qsh + 191);
    const auto *qsh_192 = buffer.data(qsh + 192);
    const auto *qsh_193 = buffer.data(qsh + 193);
    const auto *qsh_194 = buffer.data(qsh + 194);
    const auto *qsh_195 = buffer.data(qsh + 195);
    const auto *qsh_196 = buffer.data(qsh + 196);
    const auto *qsh_197 = buffer.data(qsh + 197);
    const auto *qsh_198 = buffer.data(qsh + 198);
    const auto *qsh_203 = buffer.data(qsh + 203);
    const auto *qsh_204 = buffer.data(qsh + 204);
    const auto *qsh_205 = buffer.data(qsh + 205);
    const auto *qsh_206 = buffer.data(qsh + 206);
    const auto *qsh_207 = buffer.data(qsh + 207);
    const auto *qsh_208 = buffer.data(qsh + 208);
    const auto *qsh_209 = buffer.data(qsh + 209);
    const auto *qsh_210 = buffer.data(qsh + 210);
    const auto *qsh_211 = buffer.data(qsh + 211);
    const auto *qsh_212 = buffer.data(qsh + 212);
    const auto *qsh_213 = buffer.data(qsh + 213);
    const auto *qsh_215 = buffer.data(qsh + 215);
    const auto *qsh_216 = buffer.data(qsh + 216);
    const auto *qsh_217 = buffer.data(qsh + 217);
    const auto *qsh_219 = buffer.data(qsh + 219);
    const auto *qsh_220 = buffer.data(qsh + 220);
    const auto *qsh_225 = buffer.data(qsh + 225);
    const auto *qsh_226 = buffer.data(qsh + 226);
    const auto *qsh_227 = buffer.data(qsh + 227);
    const auto *qsh_228 = buffer.data(qsh + 228);
    const auto *qsh_229 = buffer.data(qsh + 229);
    const auto *qsh_230 = buffer.data(qsh + 230);
    const auto *qsh_231 = buffer.data(qsh + 231);
    const auto *qsh_233 = buffer.data(qsh + 233);
    const auto *qsh_234 = buffer.data(qsh + 234);
    const auto *qsh_236 = buffer.data(qsh + 236);
    const auto *qsh_237 = buffer.data(qsh + 237);
    const auto *qsh_240 = buffer.data(qsh + 240);
    const auto *qsh_245 = buffer.data(qsh + 245);
    const auto *qsh_246 = buffer.data(qsh + 246);
    const auto *qsh_247 = buffer.data(qsh + 247);
    const auto *qsh_248 = buffer.data(qsh + 248);
    const auto *qsh_249 = buffer.data(qsh + 249);
    const auto *qsh_250 = buffer.data(qsh + 250);
    const auto *qsh_251 = buffer.data(qsh + 251);
    const auto *qsh_252 = buffer.data(qsh + 252);
    const auto *qsh_254 = buffer.data(qsh + 254);
    const auto *qsh_255 = buffer.data(qsh + 255);
    const auto *qsh_257 = buffer.data(qsh + 257);
    const auto *qsh_258 = buffer.data(qsh + 258);
    const auto *qsh_261 = buffer.data(qsh + 261);
    const auto *qsh_262 = buffer.data(qsh + 262);
    const auto *qsh_264 = buffer.data(qsh + 264);
    const auto *qsh_266 = buffer.data(qsh + 266);
    const auto *qsh_267 = buffer.data(qsh + 267);
    const auto *qsh_268 = buffer.data(qsh + 268);
    const auto *qsh_269 = buffer.data(qsh + 269);
    const auto *qsh_270 = buffer.data(qsh + 270);
    const auto *qsh_271 = buffer.data(qsh + 271);
    const auto *qsh_272 = buffer.data(qsh + 272);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, osh_123, osh_124, osh_125, qsg0_133, \
                         qsg0_134, qsg1_133, qsg1_134, qsh_186, qsh_187, \
                         qsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * osh_123[k]
                   + f_6 * qsg0_133[k]
                   - f_7 * qsg1_133[k]
                   + f_3 * pc_y[k] * qsh_186[k];

        t_249[k] = f_11 * osh_124[k]
                   + f_4 * qsg0_134[k]
                   - f_5 * qsg1_134[k]
                   + f_3 * pc_y[k] * qsh_187[k];

        t_250[k] = f_11 * osh_125[k]
                   + f_3 * pc_y[k] * qsh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pc_x, pc_y, pc_z, osi0_167, \
                         osh_105, osh_189, osi1_167, qsg0_135, qsg1_135, \
                         qsh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * osi0_167[k]
                   - f_10 * pc_y[k] * osi1_167[k];

        t_252[k] = f_19 * osh_189[k]
                   + f_1 * qsg0_135[k]
                   - f_2 * qsg1_135[k]
                   + f_3 * pc_x[k] * qsh_189[k];

        t_253[k] = f_3 * pc_y[k] * qsh_189[k];

        t_254[k] = f_13 * osh_105[k]
                   + f_3 * pc_z[k] * qsh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, osh_194, qsg0_135, qsg0_140, \
                         qsg1_135, qsg1_140, qsh_190, qsh_191, \
                         qsh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * qsg0_135[k]
                   - f_5 * qsg1_135[k]
                   + f_3 * pc_y[k] * qsh_190[k];

        t_256[k] = f_3 * pc_y[k] * qsh_191[k];

        t_257[k] = f_19 * osh_194[k]
                   + f_8 * qsg0_140[k]
                   - f_9 * qsg1_140[k]
                   + f_3 * pc_x[k] * qsh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_y, qsg0_136, qsg0_137, qsg1_136, qsg1_137, \
                         qsh_192, qsh_193, qsh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_6 * qsg0_136[k]
                   - f_7 * qsg1_136[k]
                   + f_3 * pc_y[k] * qsh_192[k];

        t_259[k] = f_4 * qsg0_137[k]
                   - f_5 * qsg1_137[k]
                   + f_3 * pc_y[k] * qsh_193[k];

        t_260[k] = f_3 * pc_y[k] * qsh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_y, osh_198, qsg0_138, qsg0_139, \
                         qsg0_144, qsg1_138, qsg1_139, qsg1_144, qsh_195, qsh_196, \
                         qsh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_19 * osh_198[k]
                   + f_6 * qsg0_144[k]
                   - f_7 * qsg1_144[k]
                   + f_3 * pc_x[k] * qsh_198[k];

        t_262[k] = f_8 * qsg0_138[k]
                   - f_9 * qsg1_138[k]
                   + f_3 * pc_y[k] * qsh_195[k];

        t_263[k] = f_6 * qsg0_139[k]
                   - f_7 * qsg1_139[k]
                   + f_3 * pc_y[k] * qsh_196[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, osh_203, osh_204, qsg0_140, \
                         qsg0_149, qsg1_140, qsg1_149, qsh_197, qsh_198, qsh_203, \
                         qsh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * qsg0_140[k]
                   - f_5 * qsg1_140[k]
                   + f_3 * pc_y[k] * qsh_197[k];

        t_265[k] = f_3 * pc_y[k] * qsh_198[k];

        t_266[k] = f_19 * osh_203[k]
                   + f_4 * qsg0_149[k]
                   - f_5 * qsg1_149[k]
                   + f_3 * pc_x[k] * qsh_203[k];

        t_267[k] = f_19 * osh_204[k]
                   + f_3 * pc_x[k] * qsh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pc_x, pc_y, osh_205, osh_206, \
                         osh_207, osh_209, qsh_203, qsh_205, qsh_206, qsh_207, \
                         qsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_19 * osh_205[k]
                   + f_3 * pc_x[k] * qsh_205[k];

        t_269[k] = f_19 * osh_206[k]
                   + f_3 * pc_x[k] * qsh_206[k];

        t_270[k] = f_19 * osh_207[k]
                   + f_3 * pc_x[k] * qsh_207[k];

        t_271[k] = f_3 * pc_y[k] * qsh_203[k];

        t_272[k] = f_19 * osh_209[k]
                   + f_3 * pc_x[k] * qsh_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_y, qsg0_145, qsg0_146, qsg0_147, qsg1_145, \
                         qsg1_146, qsg1_147, qsh_204, qsh_205, \
                         qsh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * qsg0_145[k]
                   - f_2 * qsg1_145[k]
                   + f_3 * pc_y[k] * qsh_204[k];

        t_274[k] = f_16 * qsg0_146[k]
                   - f_17 * qsg1_146[k]
                   + f_3 * pc_y[k] * qsh_205[k];

        t_275[k] = f_8 * qsg0_147[k]
                   - f_9 * qsg1_147[k]
                   + f_3 * pc_y[k] * qsh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, osh_125, qsg0_148, qsg0_149, \
                         qsg1_148, qsg1_149, qsh_207, qsh_208, \
                         qsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * qsg0_148[k]
                   - f_7 * qsg1_148[k]
                   + f_3 * pc_y[k] * qsh_207[k];

        t_277[k] = f_4 * qsg0_149[k]
                   - f_5 * qsg1_149[k]
                   + f_3 * pc_y[k] * qsh_208[k];

        t_278[k] = f_3 * pc_y[k] * qsh_209[k];

        t_279[k] = f_13 * osh_125[k]
                   + f_1 * qsg0_149[k]
                   - f_2 * qsg1_149[k]
                   + f_3 * pc_z[k] * qsh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, osh_126, osh_210, \
                         osh_213, qsg0_150, qsg0_153, qsg1_150, qsg1_153, qsh_210, \
                         qsh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_20 * osh_210[k]
                   + f_1 * qsg0_150[k]
                   - f_2 * qsg1_150[k]
                   + f_3 * pc_x[k] * qsh_210[k];

        t_281[k] = f_14 * osh_126[k]
                   + f_3 * pc_y[k] * qsh_210[k];

        t_282[k] = f_3 * pc_z[k] * qsh_210[k];

        t_283[k] = f_20 * osh_213[k]
                   + f_8 * qsg0_153[k]
                   - f_9 * qsg1_153[k]
                   + f_3 * pc_x[k] * qsh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_z, osh_216, qsg0_150, qsg0_156, \
                         qsg1_150, qsg1_156, qsh_211, qsh_212, qsh_213, \
                         qsh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * qsh_211[k];

        t_285[k] = f_4 * qsg0_150[k]
                   - f_5 * qsg1_150[k]
                   + f_3 * pc_z[k] * qsh_212[k];

        t_286[k] = f_20 * osh_216[k]
                   + f_6 * qsg0_156[k]
                   - f_7 * qsg1_156[k]
                   + f_3 * pc_x[k] * qsh_216[k];

        t_287[k] = f_3 * pc_z[k] * qsh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, pc_z, osh_131, osh_220, \
                         qsg0_152, qsg0_160, qsg1_152, qsg1_160, qsh_215, qsh_216, \
                         qsh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_14 * osh_131[k]
                   + f_3 * pc_y[k] * qsh_215[k];

        t_289[k] = f_6 * qsg0_152[k]
                   - f_7 * qsg1_152[k]
                   + f_3 * pc_z[k] * qsh_215[k];

        t_290[k] = f_20 * osh_220[k]
                   + f_4 * qsg0_160[k]
                   - f_5 * qsg1_160[k]
                   + f_3 * pc_x[k] * qsh_220[k];

        t_291[k] = f_3 * pc_z[k] * qsh_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, osh_135, osh_225, \
                         qsg0_153, qsg0_155, qsg1_153, qsg1_155, qsh_217, qsh_219, \
                         qsh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * qsg0_153[k]
                   - f_5 * qsg1_153[k]
                   + f_3 * pc_z[k] * qsh_217[k];

        t_293[k] = f_14 * osh_135[k]
                   + f_3 * pc_y[k] * qsh_219[k];

        t_294[k] = f_8 * qsg0_155[k]
                   - f_9 * qsg1_155[k]
                   + f_3 * pc_z[k] * qsh_219[k];

        t_295[k] = f_20 * osh_225[k]
                   + f_3 * pc_x[k] * qsh_225[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, pc_z, osh_227, osh_228, \
                         osh_229, osh_230, qsh_220, qsh_227, qsh_228, qsh_229, \
                         qsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_3 * pc_z[k] * qsh_220[k];

        t_297[k] = f_20 * osh_227[k]
                   + f_3 * pc_x[k] * qsh_227[k];

        t_298[k] = f_20 * osh_228[k]
                   + f_3 * pc_x[k] * qsh_228[k];

        t_299[k] = f_20 * osh_229[k]
                   + f_3 * pc_x[k] * qsh_229[k];

        t_300[k] = f_20 * osh_230[k]
                   + f_3 * pc_x[k] * qsh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_y, pc_z, osh_141, qsg0_160, qsg0_161, \
                         qsg1_160, qsg1_161, qsh_225, qsh_226, \
                         qsh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * osh_141[k]
                   + f_1 * qsg0_160[k]
                   - f_2 * qsg1_160[k]
                   + f_3 * pc_y[k] * qsh_225[k];

        t_302[k] = f_3 * pc_z[k] * qsh_225[k];

        t_303[k] = f_4 * qsg0_160[k]
                   - f_5 * qsg1_160[k]
                   + f_3 * pc_z[k] * qsh_226[k];

        t_304[k] = f_6 * qsg0_161[k]
                   - f_7 * qsg1_161[k]
                   + f_3 * pc_z[k] * qsh_227[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pc_y, pc_z, osi0_168, osh_146, \
                         osi1_168, qsg0_162, qsg0_164, qsg1_162, qsg1_164, qsh_228, \
                         qsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_8 * qsg0_162[k]
                   - f_9 * qsg1_162[k]
                   + f_3 * pc_z[k] * qsh_228[k];

        t_306[k] = f_14 * osh_146[k]
                   + f_3 * pc_y[k] * qsh_230[k];

        t_307[k] = f_1 * qsg0_164[k]
                   - f_2 * qsg1_164[k]
                   + f_3 * pc_z[k] * qsh_230[k];

        t_308[k] = pa_z[k] * osi0_168[k]
                   - f_10 * pc_z[k] * osi1_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_z, pc_y, pc_z, osi0_171, osh_126, \
                         osh_147, osh_149, osi1_171, qsh_231, qsh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_13 * osh_147[k]
                   + f_3 * pc_y[k] * qsh_231[k];

        t_310[k] = f_11 * osh_126[k]
                   + f_3 * pc_z[k] * qsh_231[k];

        t_311[k] = pa_z[k] * osi0_171[k]
                   - f_10 * pc_z[k] * osi1_171[k];

        t_312[k] = f_13 * osh_149[k]
                   + f_3 * pc_y[k] * qsh_233[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_z, pc_x, pc_z, osi0_174, osh_129, osh_236, \
                         osi1_174, qsg0_170, qsg1_170, qsh_234, \
                         qsh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_20 * osh_236[k]
                   + f_8 * qsg0_170[k]
                   - f_9 * qsg1_170[k]
                   + f_3 * pc_x[k] * qsh_236[k];

        t_314[k] = pa_z[k] * osi0_174[k]
                   - f_10 * pc_z[k] * osi1_174[k];

        t_315[k] = f_11 * osh_129[k]
                   + f_3 * pc_z[k] * qsh_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_z, pc_x, pc_y, pc_z, osi0_178, osh_152, \
                         osh_240, osi1_178, qsg0_174, qsg1_174, qsh_236, \
                         qsh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * osh_152[k]
                   + f_3 * pc_y[k] * qsh_236[k];

        t_317[k] = f_20 * osh_240[k]
                   + f_6 * qsg0_174[k]
                   - f_7 * qsg1_174[k]
                   + f_3 * pc_x[k] * qsh_240[k];

        t_318[k] = pa_z[k] * osi0_178[k]
                   - f_10 * pc_z[k] * osi1_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_y, pc_z, osi0_180, osh_132, osh_133, \
                         osh_156, osi1_180, qsh_237, qsh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_11 * osh_132[k]
                   + f_3 * pc_z[k] * qsh_237[k];

        t_320[k] = pa_z[k] * osi0_180[k]
                   + f_12 * osh_133[k]
                   - f_10 * pc_z[k] * osi1_180[k];

        t_321[k] = f_13 * osh_156[k]
                   + f_3 * pc_y[k] * qsh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, osh_245, osh_246, osh_247, osh_248, \
                         qsg0_179, qsg1_179, qsh_245, qsh_246, qsh_247, \
                         qsh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_20 * osh_245[k]
                   + f_4 * qsg0_179[k]
                   - f_5 * qsg1_179[k]
                   + f_3 * pc_x[k] * qsh_245[k];

        t_323[k] = f_20 * osh_246[k]
                   + f_3 * pc_x[k] * qsh_246[k];

        t_324[k] = f_20 * osh_247[k]
                   + f_3 * pc_x[k] * qsh_247[k];

        t_325[k] = f_20 * osh_248[k]
                   + f_3 * pc_x[k] * qsh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pa_z, pc_x, pc_z, osi0_189, osh_249, \
                         osh_250, osh_251, osi1_189, qsh_249, qsh_250, \
                         qsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_20 * osh_249[k]
                   + f_3 * pc_x[k] * qsh_249[k];

        t_327[k] = f_20 * osh_250[k]
                   + f_3 * pc_x[k] * qsh_250[k];

        t_328[k] = f_20 * osh_251[k]
                   + f_3 * pc_x[k] * qsh_251[k];

        t_329[k] = pa_z[k] * osi0_189[k]
                   - f_10 * pc_z[k] * osi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, osh_141, osh_164, osh_165, qsg0_177, \
                         qsg0_178, qsg1_177, qsg1_178, qsh_246, qsh_248, \
                         qsh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * osh_141[k]
                   + f_3 * pc_z[k] * qsh_246[k];

        t_331[k] = f_13 * osh_164[k]
                   + f_8 * qsg0_177[k]
                   - f_9 * qsg1_177[k]
                   + f_3 * pc_y[k] * qsh_248[k];

        t_332[k] = f_13 * osh_165[k]
                   + f_6 * qsg0_178[k]
                   - f_7 * qsg1_178[k]
                   + f_3 * pc_y[k] * qsh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, osh_146, osh_166, osh_167, qsg0_179, \
                         qsg1_179, qsh_250, qsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * osh_166[k]
                   + f_4 * qsg0_179[k]
                   - f_5 * qsg1_179[k]
                   + f_3 * pc_y[k] * qsh_250[k];

        t_334[k] = f_13 * osh_167[k]
                   + f_3 * pc_y[k] * qsh_251[k];

        t_335[k] = f_11 * osh_146[k]
                   + f_1 * qsg0_179[k]
                   - f_2 * qsg1_179[k]
                   + f_3 * pc_z[k] * qsh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, osh_147, osh_168, osh_252, \
                         qsg0_180, qsg1_180, qsh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_20 * osh_252[k]
                   + f_1 * qsg0_180[k]
                   - f_2 * qsg1_180[k]
                   + f_3 * pc_x[k] * qsh_252[k];

        t_337[k] = f_12 * osh_168[k]
                   + f_3 * pc_y[k] * qsh_252[k];

        t_338[k] = f_12 * osh_147[k]
                   + f_3 * pc_z[k] * qsh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, osh_170, osh_255, osh_257, qsg0_183, \
                         qsg0_185, qsg1_183, qsg1_185, qsh_254, qsh_255, \
                         qsh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_20 * osh_255[k]
                   + f_8 * qsg0_183[k]
                   - f_9 * qsg1_183[k]
                   + f_3 * pc_x[k] * qsh_255[k];

        t_340[k] = f_12 * osh_170[k]
                   + f_3 * pc_y[k] * qsh_254[k];

        t_341[k] = f_20 * osh_257[k]
                   + f_8 * qsg0_185[k]
                   - f_9 * qsg1_185[k]
                   + f_3 * pc_x[k] * qsh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, osh_150, osh_173, osh_258, \
                         qsg0_186, qsg1_186, qsh_255, qsh_257, \
                         qsh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_20 * osh_258[k]
                   + f_6 * qsg0_186[k]
                   - f_7 * qsg1_186[k]
                   + f_3 * pc_x[k] * qsh_258[k];

        t_343[k] = f_12 * osh_150[k]
                   + f_3 * pc_z[k] * qsh_255[k];

        t_344[k] = f_12 * osh_173[k]
                   + f_3 * pc_y[k] * qsh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, osh_153, osh_261, osh_262, qsg0_189, \
                         qsg0_190, qsg1_189, qsg1_190, qsh_258, qsh_261, \
                         qsh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_20 * osh_261[k]
                   + f_6 * qsg0_189[k]
                   - f_7 * qsg1_189[k]
                   + f_3 * pc_x[k] * qsh_261[k];

        t_346[k] = f_20 * osh_262[k]
                   + f_4 * qsg0_190[k]
                   - f_5 * qsg1_190[k]
                   + f_3 * pc_x[k] * qsh_262[k];

        t_347[k] = f_12 * osh_153[k]
                   + f_3 * pc_z[k] * qsh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, osh_177, osh_264, osh_266, qsg0_192, \
                         qsg0_194, qsg1_192, qsg1_194, qsh_261, qsh_264, \
                         qsh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_20 * osh_264[k]
                   + f_4 * qsg0_192[k]
                   - f_5 * qsg1_192[k]
                   + f_3 * pc_x[k] * qsh_264[k];

        t_349[k] = f_12 * osh_177[k]
                   + f_3 * pc_y[k] * qsh_261[k];

        t_350[k] = f_20 * osh_266[k]
                   + f_4 * qsg0_194[k]
                   - f_5 * qsg1_194[k]
                   + f_3 * pc_x[k] * qsh_266[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, osh_267, osh_268, osh_269, \
                         osh_270, osh_271, qsh_267, qsh_268, qsh_269, qsh_270, \
                         qsh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_20 * osh_267[k]
                   + f_3 * pc_x[k] * qsh_267[k];

        t_352[k] = f_20 * osh_268[k]
                   + f_3 * pc_x[k] * qsh_268[k];

        t_353[k] = f_20 * osh_269[k]
                   + f_3 * pc_x[k] * qsh_269[k];

        t_354[k] = f_20 * osh_270[k]
                   + f_3 * pc_x[k] * qsh_270[k];

        t_355[k] = f_20 * osh_271[k]
                   + f_3 * pc_x[k] * qsh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, osh_162, osh_183, osh_272, \
                         qsg0_190, qsg1_190, qsh_267, qsh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_20 * osh_272[k]
                   + f_3 * pc_x[k] * qsh_272[k];

        t_357[k] = f_12 * osh_183[k]
                   + f_1 * qsg0_190[k]
                   - f_2 * qsg1_190[k]
                   + f_3 * pc_y[k] * qsh_267[k];

        t_358[k] = f_12 * osh_162[k]
                   + f_3 * pc_z[k] * qsh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, osh_185, osh_186, osh_187, qsg0_192, \
                         qsg0_193, qsg0_194, qsg1_192, qsg1_193, qsg1_194, qsh_269, qsh_270, \
                         qsh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * osh_185[k]
                   + f_8 * qsg0_192[k]
                   - f_9 * qsg1_192[k]
                   + f_3 * pc_y[k] * qsh_269[k];

        t_360[k] = f_12 * osh_186[k]
                   + f_6 * qsg0_193[k]
                   - f_7 * qsg1_193[k]
                   + f_3 * pc_y[k] * qsh_270[k];

        t_361[k] = f_12 * osh_187[k]
                   + f_4 * qsg0_194[k]
                   - f_5 * qsg1_194[k]
                   + f_3 * pc_y[k] * qsh_271[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;

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
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_252 = buffer.data(osi0 + 252);
    const auto *osi0_255 = buffer.data(osi0 + 255);
    const auto *osi0_257 = buffer.data(osi0 + 257);
    const auto *osi0_258 = buffer.data(osi0 + 258);
    const auto *osi0_261 = buffer.data(osi0 + 261);
    const auto *osi0_262 = buffer.data(osi0 + 262);
    const auto *osi0_264 = buffer.data(osi0 + 264);
    const auto *osi0_266 = buffer.data(osi0 + 266);
    const auto *osi0_279 = buffer.data(osi0 + 279);
    const auto *osi0_280 = buffer.data(osi0 + 280);
    const auto *osi0_283 = buffer.data(osi0 + 283);
    const auto *osi0_286 = buffer.data(osi0 + 286);
    const auto *osi0_290 = buffer.data(osi0 + 290);
    const auto *osi0_292 = buffer.data(osi0 + 292);
    const auto *osi0_301 = buffer.data(osi0 + 301);

    const auto *osh_167 = buffer.data(osh + 167);
    const auto *osh_168 = buffer.data(osh + 168);
    const auto *osh_171 = buffer.data(osh + 171);
    const auto *osh_174 = buffer.data(osh + 174);
    const auto *osh_183 = buffer.data(osh + 183);
    const auto *osh_188 = buffer.data(osh + 188);
    const auto *osh_189 = buffer.data(osh + 189);
    const auto *osh_190 = buffer.data(osh + 190);
    const auto *osh_191 = buffer.data(osh + 191);
    const auto *osh_192 = buffer.data(osh + 192);
    const auto *osh_194 = buffer.data(osh + 194);
    const auto *osh_195 = buffer.data(osh + 195);
    const auto *osh_197 = buffer.data(osh + 197);
    const auto *osh_198 = buffer.data(osh + 198);
    const auto *osh_204 = buffer.data(osh + 204);
    const auto *osh_206 = buffer.data(osh + 206);
    const auto *osh_207 = buffer.data(osh + 207);
    const auto *osh_208 = buffer.data(osh + 208);
    const auto *osh_209 = buffer.data(osh + 209);
    const auto *osh_210 = buffer.data(osh + 210);
    const auto *osh_213 = buffer.data(osh + 213);
    const auto *osh_215 = buffer.data(osh + 215);
    const auto *osh_216 = buffer.data(osh + 216);
    const auto *osh_217 = buffer.data(osh + 217);
    const auto *osh_219 = buffer.data(osh + 219);
    const auto *osh_225 = buffer.data(osh + 225);
    const auto *osh_230 = buffer.data(osh + 230);
    const auto *osh_231 = buffer.data(osh + 231);
    const auto *osh_233 = buffer.data(osh + 233);
    const auto *osh_236 = buffer.data(osh + 236);
    const auto *osh_240 = buffer.data(osh + 240);
    const auto *osh_248 = buffer.data(osh + 248);
    const auto *osh_249 = buffer.data(osh + 249);
    const auto *osh_250 = buffer.data(osh + 250);
    const auto *osh_251 = buffer.data(osh + 251);
    const auto *osh_252 = buffer.data(osh + 252);
    const auto *osh_288 = buffer.data(osh + 288);
    const auto *osh_289 = buffer.data(osh + 289);
    const auto *osh_290 = buffer.data(osh + 290);
    const auto *osh_291 = buffer.data(osh + 291);
    const auto *osh_292 = buffer.data(osh + 292);
    const auto *osh_293 = buffer.data(osh + 293);
    const auto *osh_294 = buffer.data(osh + 294);
    const auto *osh_299 = buffer.data(osh + 299);
    const auto *osh_303 = buffer.data(osh + 303);
    const auto *osh_308 = buffer.data(osh + 308);
    const auto *osh_309 = buffer.data(osh + 309);
    const auto *osh_310 = buffer.data(osh + 310);
    const auto *osh_311 = buffer.data(osh + 311);
    const auto *osh_312 = buffer.data(osh + 312);
    const auto *osh_314 = buffer.data(osh + 314);
    const auto *osh_315 = buffer.data(osh + 315);
    const auto *osh_318 = buffer.data(osh + 318);
    const auto *osh_321 = buffer.data(osh + 321);
    const auto *osh_325 = buffer.data(osh + 325);
    const auto *osh_330 = buffer.data(osh + 330);
    const auto *osh_332 = buffer.data(osh + 332);
    const auto *osh_333 = buffer.data(osh + 333);
    const auto *osh_334 = buffer.data(osh + 334);
    const auto *osh_335 = buffer.data(osh + 335);
    const auto *osh_341 = buffer.data(osh + 341);
    const auto *osh_345 = buffer.data(osh + 345);
    const auto *osh_350 = buffer.data(osh + 350);
    const auto *osh_351 = buffer.data(osh + 351);
    const auto *osh_352 = buffer.data(osh + 352);
    const auto *osh_353 = buffer.data(osh + 353);
    const auto *osh_354 = buffer.data(osh + 354);
    const auto *osh_355 = buffer.data(osh + 355);
    const auto *osh_356 = buffer.data(osh + 356);
    const auto *osh_357 = buffer.data(osh + 357);

    const auto *osi1_252 = buffer.data(osi1 + 252);
    const auto *osi1_255 = buffer.data(osi1 + 255);
    const auto *osi1_257 = buffer.data(osi1 + 257);
    const auto *osi1_258 = buffer.data(osi1 + 258);
    const auto *osi1_261 = buffer.data(osi1 + 261);
    const auto *osi1_262 = buffer.data(osi1 + 262);
    const auto *osi1_264 = buffer.data(osi1 + 264);
    const auto *osi1_266 = buffer.data(osi1 + 266);
    const auto *osi1_279 = buffer.data(osi1 + 279);
    const auto *osi1_280 = buffer.data(osi1 + 280);
    const auto *osi1_283 = buffer.data(osi1 + 283);
    const auto *osi1_286 = buffer.data(osi1 + 286);
    const auto *osi1_290 = buffer.data(osi1 + 290);
    const auto *osi1_292 = buffer.data(osi1 + 292);
    const auto *osi1_301 = buffer.data(osi1 + 301);

    const auto *qsg0_194 = buffer.data(qsg0 + 194);
    const auto *qsg0_205 = buffer.data(qsg0 + 205);
    const auto *qsg0_207 = buffer.data(qsg0 + 207);
    const auto *qsg0_208 = buffer.data(qsg0 + 208);
    const auto *qsg0_209 = buffer.data(qsg0 + 209);
    const auto *qsg0_210 = buffer.data(qsg0 + 210);
    const auto *qsg0_211 = buffer.data(qsg0 + 211);
    const auto *qsg0_212 = buffer.data(qsg0 + 212);
    const auto *qsg0_213 = buffer.data(qsg0 + 213);
    const auto *qsg0_214 = buffer.data(qsg0 + 214);
    const auto *qsg0_215 = buffer.data(qsg0 + 215);
    const auto *qsg0_219 = buffer.data(qsg0 + 219);
    const auto *qsg0_220 = buffer.data(qsg0 + 220);
    const auto *qsg0_221 = buffer.data(qsg0 + 221);
    const auto *qsg0_222 = buffer.data(qsg0 + 222);
    const auto *qsg0_223 = buffer.data(qsg0 + 223);
    const auto *qsg0_224 = buffer.data(qsg0 + 224);
    const auto *qsg0_225 = buffer.data(qsg0 + 225);
    const auto *qsg0_227 = buffer.data(qsg0 + 227);
    const auto *qsg0_228 = buffer.data(qsg0 + 228);
    const auto *qsg0_230 = buffer.data(qsg0 + 230);
    const auto *qsg0_231 = buffer.data(qsg0 + 231);
    const auto *qsg0_235 = buffer.data(qsg0 + 235);
    const auto *qsg0_236 = buffer.data(qsg0 + 236);
    const auto *qsg0_237 = buffer.data(qsg0 + 237);
    const auto *qsg0_239 = buffer.data(qsg0 + 239);
    const auto *qsg0_245 = buffer.data(qsg0 + 245);
    const auto *qsg0_249 = buffer.data(qsg0 + 249);
    const auto *qsg0_252 = buffer.data(qsg0 + 252);
    const auto *qsg0_253 = buffer.data(qsg0 + 253);
    const auto *qsg0_254 = buffer.data(qsg0 + 254);
    const auto *qsg0_255 = buffer.data(qsg0 + 255);

    const auto *qsg1_194 = buffer.data(qsg1 + 194);
    const auto *qsg1_205 = buffer.data(qsg1 + 205);
    const auto *qsg1_207 = buffer.data(qsg1 + 207);
    const auto *qsg1_208 = buffer.data(qsg1 + 208);
    const auto *qsg1_209 = buffer.data(qsg1 + 209);
    const auto *qsg1_210 = buffer.data(qsg1 + 210);
    const auto *qsg1_211 = buffer.data(qsg1 + 211);
    const auto *qsg1_212 = buffer.data(qsg1 + 212);
    const auto *qsg1_213 = buffer.data(qsg1 + 213);
    const auto *qsg1_214 = buffer.data(qsg1 + 214);
    const auto *qsg1_215 = buffer.data(qsg1 + 215);
    const auto *qsg1_219 = buffer.data(qsg1 + 219);
    const auto *qsg1_220 = buffer.data(qsg1 + 220);
    const auto *qsg1_221 = buffer.data(qsg1 + 221);
    const auto *qsg1_222 = buffer.data(qsg1 + 222);
    const auto *qsg1_223 = buffer.data(qsg1 + 223);
    const auto *qsg1_224 = buffer.data(qsg1 + 224);
    const auto *qsg1_225 = buffer.data(qsg1 + 225);
    const auto *qsg1_227 = buffer.data(qsg1 + 227);
    const auto *qsg1_228 = buffer.data(qsg1 + 228);
    const auto *qsg1_230 = buffer.data(qsg1 + 230);
    const auto *qsg1_231 = buffer.data(qsg1 + 231);
    const auto *qsg1_235 = buffer.data(qsg1 + 235);
    const auto *qsg1_236 = buffer.data(qsg1 + 236);
    const auto *qsg1_237 = buffer.data(qsg1 + 237);
    const auto *qsg1_239 = buffer.data(qsg1 + 239);
    const auto *qsg1_245 = buffer.data(qsg1 + 245);
    const auto *qsg1_249 = buffer.data(qsg1 + 249);
    const auto *qsg1_252 = buffer.data(qsg1 + 252);
    const auto *qsg1_253 = buffer.data(qsg1 + 253);
    const auto *qsg1_254 = buffer.data(qsg1 + 254);
    const auto *qsg1_255 = buffer.data(qsg1 + 255);

    const auto *qsh_272 = buffer.data(qsh + 272);
    const auto *qsh_273 = buffer.data(qsh + 273);
    const auto *qsh_275 = buffer.data(qsh + 275);
    const auto *qsh_276 = buffer.data(qsh + 276);
    const auto *qsh_278 = buffer.data(qsh + 278);
    const auto *qsh_279 = buffer.data(qsh + 279);
    const auto *qsh_282 = buffer.data(qsh + 282);
    const auto *qsh_288 = buffer.data(qsh + 288);
    const auto *qsh_289 = buffer.data(qsh + 289);
    const auto *qsh_290 = buffer.data(qsh + 290);
    const auto *qsh_291 = buffer.data(qsh + 291);
    const auto *qsh_292 = buffer.data(qsh + 292);
    const auto *qsh_293 = buffer.data(qsh + 293);
    const auto *qsh_294 = buffer.data(qsh + 294);
    const auto *qsh_295 = buffer.data(qsh + 295);
    const auto *qsh_296 = buffer.data(qsh + 296);
    const auto *qsh_297 = buffer.data(qsh + 297);
    const auto *qsh_298 = buffer.data(qsh + 298);
    const auto *qsh_299 = buffer.data(qsh + 299);
    const auto *qsh_300 = buffer.data(qsh + 300);
    const auto *qsh_301 = buffer.data(qsh + 301);
    const auto *qsh_302 = buffer.data(qsh + 302);
    const auto *qsh_303 = buffer.data(qsh + 303);
    const auto *qsh_308 = buffer.data(qsh + 308);
    const auto *qsh_309 = buffer.data(qsh + 309);
    const auto *qsh_310 = buffer.data(qsh + 310);
    const auto *qsh_311 = buffer.data(qsh + 311);
    const auto *qsh_312 = buffer.data(qsh + 312);
    const auto *qsh_313 = buffer.data(qsh + 313);
    const auto *qsh_314 = buffer.data(qsh + 314);
    const auto *qsh_315 = buffer.data(qsh + 315);
    const auto *qsh_316 = buffer.data(qsh + 316);
    const auto *qsh_317 = buffer.data(qsh + 317);
    const auto *qsh_318 = buffer.data(qsh + 318);
    const auto *qsh_320 = buffer.data(qsh + 320);
    const auto *qsh_321 = buffer.data(qsh + 321);
    const auto *qsh_322 = buffer.data(qsh + 322);
    const auto *qsh_324 = buffer.data(qsh + 324);
    const auto *qsh_325 = buffer.data(qsh + 325);
    const auto *qsh_330 = buffer.data(qsh + 330);
    const auto *qsh_331 = buffer.data(qsh + 331);
    const auto *qsh_332 = buffer.data(qsh + 332);
    const auto *qsh_333 = buffer.data(qsh + 333);
    const auto *qsh_334 = buffer.data(qsh + 334);
    const auto *qsh_335 = buffer.data(qsh + 335);
    const auto *qsh_336 = buffer.data(qsh + 336);
    const auto *qsh_338 = buffer.data(qsh + 338);
    const auto *qsh_339 = buffer.data(qsh + 339);
    const auto *qsh_341 = buffer.data(qsh + 341);
    const auto *qsh_342 = buffer.data(qsh + 342);
    const auto *qsh_345 = buffer.data(qsh + 345);
    const auto *qsh_350 = buffer.data(qsh + 350);
    const auto *qsh_351 = buffer.data(qsh + 351);
    const auto *qsh_352 = buffer.data(qsh + 352);
    const auto *qsh_353 = buffer.data(qsh + 353);
    const auto *qsh_354 = buffer.data(qsh + 354);
    const auto *qsh_355 = buffer.data(qsh + 355);
    const auto *qsh_356 = buffer.data(qsh + 356);
    const auto *qsh_357 = buffer.data(qsh + 357);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_y, pc_y, pc_z, osi0_252, osh_167, \
                         osh_188, osh_189, osi1_252, qsg0_194, qsg1_194, qsh_272, \
                         qsh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * osh_188[k]
                   + f_3 * pc_y[k] * qsh_272[k];

        t_363[k] = f_12 * osh_167[k]
                   + f_1 * qsg0_194[k]
                   - f_2 * qsg1_194[k]
                   + f_3 * pc_z[k] * qsh_272[k];

        t_364[k] = pa_y[k] * osi0_252[k]
                   - f_10 * pc_y[k] * osi1_252[k];

        t_365[k] = f_11 * osh_189[k]
                   + f_3 * pc_y[k] * qsh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pc_y, pc_z, osi0_255, osi0_257, \
                         osh_168, osh_190, osh_191, osi1_255, osi1_257, qsh_273, \
                         qsh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * osh_168[k]
                   + f_3 * pc_z[k] * qsh_273[k];

        t_367[k] = pa_y[k] * osi0_255[k]
                   + f_12 * osh_190[k]
                   - f_10 * pc_y[k] * osi1_255[k];

        t_368[k] = f_11 * osh_191[k]
                   + f_3 * pc_y[k] * qsh_275[k];

        t_369[k] = pa_y[k] * osi0_257[k]
                   - f_10 * pc_y[k] * osi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_y, pc_y, pc_z, osi0_258, osi0_261, \
                         osh_171, osh_192, osh_194, osi1_258, osi1_261, qsh_276, \
                         qsh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_y[k] * osi0_258[k]
                   + f_13 * osh_192[k]
                   - f_10 * pc_y[k] * osi1_258[k];

        t_371[k] = f_13 * osh_171[k]
                   + f_3 * pc_z[k] * qsh_276[k];

        t_372[k] = f_11 * osh_194[k]
                   + f_3 * pc_y[k] * qsh_278[k];

        t_373[k] = pa_y[k] * osi0_261[k]
                   - f_10 * pc_y[k] * osi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_y, pc_y, pc_z, osi0_262, osi0_264, osh_174, \
                         osh_195, osh_197, osi1_262, osi1_264, \
                         qsh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * osi0_262[k]
                   + f_14 * osh_195[k]
                   - f_10 * pc_y[k] * osi1_262[k];

        t_375[k] = f_13 * osh_174[k]
                   + f_3 * pc_z[k] * qsh_279[k];

        t_376[k] = pa_y[k] * osi0_264[k]
                   + f_12 * osh_197[k]
                   - f_10 * pc_y[k] * osi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pc_x, pc_y, osi0_266, osh_198, \
                         osh_288, osh_289, osi1_266, qsh_282, qsh_288, \
                         qsh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * osh_198[k]
                   + f_3 * pc_y[k] * qsh_282[k];

        t_378[k] = pa_y[k] * osi0_266[k]
                   - f_10 * pc_y[k] * osi1_266[k];

        t_379[k] = f_20 * osh_288[k]
                   + f_3 * pc_x[k] * qsh_288[k];

        t_380[k] = f_20 * osh_289[k]
                   + f_3 * pc_x[k] * qsh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, osh_290, osh_291, osh_292, osh_293, \
                         qsh_290, qsh_291, qsh_292, qsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_20 * osh_290[k]
                   + f_3 * pc_x[k] * qsh_290[k];

        t_382[k] = f_20 * osh_291[k]
                   + f_3 * pc_x[k] * qsh_291[k];

        t_383[k] = f_20 * osh_292[k]
                   + f_3 * pc_x[k] * qsh_292[k];

        t_384[k] = f_20 * osh_293[k]
                   + f_3 * pc_x[k] * qsh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, osh_183, osh_204, osh_206, qsg0_205, \
                         qsg0_207, qsg1_205, qsg1_207, qsh_288, \
                         qsh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * osh_204[k]
                   + f_1 * qsg0_205[k]
                   - f_2 * qsg1_205[k]
                   + f_3 * pc_y[k] * qsh_288[k];

        t_386[k] = f_13 * osh_183[k]
                   + f_3 * pc_z[k] * qsh_288[k];

        t_387[k] = f_11 * osh_206[k]
                   + f_8 * qsg0_207[k]
                   - f_9 * qsg1_207[k]
                   + f_3 * pc_y[k] * qsh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, osh_207, osh_208, osh_209, qsg0_208, \
                         qsg0_209, qsg1_208, qsg1_209, qsh_291, qsh_292, \
                         qsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * osh_207[k]
                   + f_6 * qsg0_208[k]
                   - f_7 * qsg1_208[k]
                   + f_3 * pc_y[k] * qsh_291[k];

        t_389[k] = f_11 * osh_208[k]
                   + f_4 * qsg0_209[k]
                   - f_5 * qsg1_209[k]
                   + f_3 * pc_y[k] * qsh_292[k];

        t_390[k] = f_11 * osh_209[k]
                   + f_3 * pc_y[k] * qsh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_x, pc_y, pc_z, osi0_279, \
                         osh_189, osh_294, osi1_279, qsg0_210, qsg1_210, \
                         qsh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_y[k] * osi0_279[k]
                   - f_10 * pc_y[k] * osi1_279[k];

        t_392[k] = f_20 * osh_294[k]
                   + f_1 * qsg0_210[k]
                   - f_2 * qsg1_210[k]
                   + f_3 * pc_x[k] * qsh_294[k];

        t_393[k] = f_3 * pc_y[k] * qsh_294[k];

        t_394[k] = f_14 * osh_189[k]
                   + f_3 * pc_z[k] * qsh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, osh_299, qsg0_210, qsg0_215, \
                         qsg1_210, qsg1_215, qsh_295, qsh_296, \
                         qsh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_4 * qsg0_210[k]
                   - f_5 * qsg1_210[k]
                   + f_3 * pc_y[k] * qsh_295[k];

        t_396[k] = f_3 * pc_y[k] * qsh_296[k];

        t_397[k] = f_20 * osh_299[k]
                   + f_8 * qsg0_215[k]
                   - f_9 * qsg1_215[k]
                   + f_3 * pc_x[k] * qsh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, qsg0_211, qsg0_212, qsg1_211, qsg1_212, \
                         qsh_297, qsh_298, qsh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_6 * qsg0_211[k]
                   - f_7 * qsg1_211[k]
                   + f_3 * pc_y[k] * qsh_297[k];

        t_399[k] = f_4 * qsg0_212[k]
                   - f_5 * qsg1_212[k]
                   + f_3 * pc_y[k] * qsh_298[k];

        t_400[k] = f_3 * pc_y[k] * qsh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_y, osh_303, qsg0_213, qsg0_214, \
                         qsg0_219, qsg1_213, qsg1_214, qsg1_219, qsh_300, qsh_301, \
                         qsh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_20 * osh_303[k]
                   + f_6 * qsg0_219[k]
                   - f_7 * qsg1_219[k]
                   + f_3 * pc_x[k] * qsh_303[k];

        t_402[k] = f_8 * qsg0_213[k]
                   - f_9 * qsg1_213[k]
                   + f_3 * pc_y[k] * qsh_300[k];

        t_403[k] = f_6 * qsg0_214[k]
                   - f_7 * qsg1_214[k]
                   + f_3 * pc_y[k] * qsh_301[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pc_x, pc_y, osh_308, osh_309, qsg0_215, \
                         qsg0_224, qsg1_215, qsg1_224, qsh_302, qsh_303, qsh_308, \
                         qsh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * qsg0_215[k]
                   - f_5 * qsg1_215[k]
                   + f_3 * pc_y[k] * qsh_302[k];

        t_405[k] = f_3 * pc_y[k] * qsh_303[k];

        t_406[k] = f_20 * osh_308[k]
                   + f_4 * qsg0_224[k]
                   - f_5 * qsg1_224[k]
                   + f_3 * pc_x[k] * qsh_308[k];

        t_407[k] = f_20 * osh_309[k]
                   + f_3 * pc_x[k] * qsh_309[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pc_x, pc_y, osh_310, osh_311, \
                         osh_312, osh_314, qsh_308, qsh_310, qsh_311, qsh_312, \
                         qsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_20 * osh_310[k]
                   + f_3 * pc_x[k] * qsh_310[k];

        t_409[k] = f_20 * osh_311[k]
                   + f_3 * pc_x[k] * qsh_311[k];

        t_410[k] = f_20 * osh_312[k]
                   + f_3 * pc_x[k] * qsh_312[k];

        t_411[k] = f_3 * pc_y[k] * qsh_308[k];

        t_412[k] = f_20 * osh_314[k]
                   + f_3 * pc_x[k] * qsh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pc_y, qsg0_220, qsg0_221, qsg0_222, qsg1_220, \
                         qsg1_221, qsg1_222, qsh_309, qsh_310, \
                         qsh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_1 * qsg0_220[k]
                   - f_2 * qsg1_220[k]
                   + f_3 * pc_y[k] * qsh_309[k];

        t_414[k] = f_16 * qsg0_221[k]
                   - f_17 * qsg1_221[k]
                   + f_3 * pc_y[k] * qsh_310[k];

        t_415[k] = f_8 * qsg0_222[k]
                   - f_9 * qsg1_222[k]
                   + f_3 * pc_y[k] * qsh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, osh_209, qsg0_223, qsg0_224, \
                         qsg1_223, qsg1_224, qsh_312, qsh_313, \
                         qsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * qsg0_223[k]
                   - f_7 * qsg1_223[k]
                   + f_3 * pc_y[k] * qsh_312[k];

        t_417[k] = f_4 * qsg0_224[k]
                   - f_5 * qsg1_224[k]
                   + f_3 * pc_y[k] * qsh_313[k];

        t_418[k] = f_3 * pc_y[k] * qsh_314[k];

        t_419[k] = f_14 * osh_209[k]
                   + f_1 * qsg0_224[k]
                   - f_2 * qsg1_224[k]
                   + f_3 * pc_z[k] * qsh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, osh_210, osh_315, \
                         osh_318, qsg0_225, qsg0_228, qsg1_225, qsg1_228, qsh_315, \
                         qsh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_21 * osh_315[k]
                   + f_1 * qsg0_225[k]
                   - f_2 * qsg1_225[k]
                   + f_3 * pc_x[k] * qsh_315[k];

        t_421[k] = f_22 * osh_210[k]
                   + f_3 * pc_y[k] * qsh_315[k];

        t_422[k] = f_3 * pc_z[k] * qsh_315[k];

        t_423[k] = f_21 * osh_318[k]
                   + f_8 * qsg0_228[k]
                   - f_9 * qsg1_228[k]
                   + f_3 * pc_x[k] * qsh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pc_x, pc_z, osh_321, qsg0_225, qsg0_231, \
                         qsg1_225, qsg1_231, qsh_316, qsh_317, qsh_318, \
                         qsh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_z[k] * qsh_316[k];

        t_425[k] = f_4 * qsg0_225[k]
                   - f_5 * qsg1_225[k]
                   + f_3 * pc_z[k] * qsh_317[k];

        t_426[k] = f_21 * osh_321[k]
                   + f_6 * qsg0_231[k]
                   - f_7 * qsg1_231[k]
                   + f_3 * pc_x[k] * qsh_321[k];

        t_427[k] = f_3 * pc_z[k] * qsh_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pc_x, pc_y, pc_z, osh_215, osh_325, \
                         qsg0_227, qsg0_235, qsg1_227, qsg1_235, qsh_320, qsh_321, \
                         qsh_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_22 * osh_215[k]
                   + f_3 * pc_y[k] * qsh_320[k];

        t_429[k] = f_6 * qsg0_227[k]
                   - f_7 * qsg1_227[k]
                   + f_3 * pc_z[k] * qsh_320[k];

        t_430[k] = f_21 * osh_325[k]
                   + f_4 * qsg0_235[k]
                   - f_5 * qsg1_235[k]
                   + f_3 * pc_x[k] * qsh_325[k];

        t_431[k] = f_3 * pc_z[k] * qsh_321[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, pc_y, pc_z, osh_219, osh_330, \
                         qsg0_228, qsg0_230, qsg1_228, qsg1_230, qsh_322, qsh_324, \
                         qsh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_4 * qsg0_228[k]
                   - f_5 * qsg1_228[k]
                   + f_3 * pc_z[k] * qsh_322[k];

        t_433[k] = f_22 * osh_219[k]
                   + f_3 * pc_y[k] * qsh_324[k];

        t_434[k] = f_8 * qsg0_230[k]
                   - f_9 * qsg1_230[k]
                   + f_3 * pc_z[k] * qsh_324[k];

        t_435[k] = f_21 * osh_330[k]
                   + f_3 * pc_x[k] * qsh_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_z, osh_332, osh_333, \
                         osh_334, osh_335, qsh_325, qsh_332, qsh_333, qsh_334, \
                         qsh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_3 * pc_z[k] * qsh_325[k];

        t_437[k] = f_21 * osh_332[k]
                   + f_3 * pc_x[k] * qsh_332[k];

        t_438[k] = f_21 * osh_333[k]
                   + f_3 * pc_x[k] * qsh_333[k];

        t_439[k] = f_21 * osh_334[k]
                   + f_3 * pc_x[k] * qsh_334[k];

        t_440[k] = f_21 * osh_335[k]
                   + f_3 * pc_x[k] * qsh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, osh_225, qsg0_235, qsg0_236, \
                         qsg1_235, qsg1_236, qsh_330, qsh_331, \
                         qsh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_22 * osh_225[k]
                   + f_1 * qsg0_235[k]
                   - f_2 * qsg1_235[k]
                   + f_3 * pc_y[k] * qsh_330[k];

        t_442[k] = f_3 * pc_z[k] * qsh_330[k];

        t_443[k] = f_4 * qsg0_235[k]
                   - f_5 * qsg1_235[k]
                   + f_3 * pc_z[k] * qsh_331[k];

        t_444[k] = f_6 * qsg0_236[k]
                   - f_7 * qsg1_236[k]
                   + f_3 * pc_z[k] * qsh_332[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pc_y, pc_z, osi0_280, osh_230, \
                         osi1_280, qsg0_237, qsg0_239, qsg1_237, qsg1_239, qsh_333, \
                         qsh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * qsg0_237[k]
                   - f_9 * qsg1_237[k]
                   + f_3 * pc_z[k] * qsh_333[k];

        t_446[k] = f_22 * osh_230[k]
                   + f_3 * pc_y[k] * qsh_335[k];

        t_447[k] = f_1 * qsg0_239[k]
                   - f_2 * qsg1_239[k]
                   + f_3 * pc_z[k] * qsh_335[k];

        t_448[k] = pa_z[k] * osi0_280[k]
                   - f_10 * pc_z[k] * osi1_280[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_z, pc_y, pc_z, osi0_283, osh_210, \
                         osh_231, osh_233, osi1_283, qsh_336, qsh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * osh_231[k]
                   + f_3 * pc_y[k] * qsh_336[k];

        t_450[k] = f_11 * osh_210[k]
                   + f_3 * pc_z[k] * qsh_336[k];

        t_451[k] = pa_z[k] * osi0_283[k]
                   - f_10 * pc_z[k] * osi1_283[k];

        t_452[k] = f_14 * osh_233[k]
                   + f_3 * pc_y[k] * qsh_338[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pa_z, pc_x, pc_z, osi0_286, osh_213, osh_341, \
                         osi1_286, qsg0_245, qsg1_245, qsh_339, \
                         qsh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_21 * osh_341[k]
                   + f_8 * qsg0_245[k]
                   - f_9 * qsg1_245[k]
                   + f_3 * pc_x[k] * qsh_341[k];

        t_454[k] = pa_z[k] * osi0_286[k]
                   - f_10 * pc_z[k] * osi1_286[k];

        t_455[k] = f_11 * osh_213[k]
                   + f_3 * pc_z[k] * qsh_339[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pa_z, pc_x, pc_y, pc_z, osi0_290, osh_236, \
                         osh_345, osi1_290, qsg0_249, qsg1_249, qsh_341, \
                         qsh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * osh_236[k]
                   + f_3 * pc_y[k] * qsh_341[k];

        t_457[k] = f_21 * osh_345[k]
                   + f_6 * qsg0_249[k]
                   - f_7 * qsg1_249[k]
                   + f_3 * pc_x[k] * qsh_345[k];

        t_458[k] = pa_z[k] * osi0_290[k]
                   - f_10 * pc_z[k] * osi1_290[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, pc_y, pc_z, osi0_292, osh_216, osh_217, \
                         osh_240, osi1_292, qsh_342, qsh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_11 * osh_216[k]
                   + f_3 * pc_z[k] * qsh_342[k];

        t_460[k] = pa_z[k] * osi0_292[k]
                   + f_12 * osh_217[k]
                   - f_10 * pc_z[k] * osi1_292[k];

        t_461[k] = f_14 * osh_240[k]
                   + f_3 * pc_y[k] * qsh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, osh_350, osh_351, osh_352, osh_353, \
                         qsg0_254, qsg1_254, qsh_350, qsh_351, qsh_352, \
                         qsh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_21 * osh_350[k]
                   + f_4 * qsg0_254[k]
                   - f_5 * qsg1_254[k]
                   + f_3 * pc_x[k] * qsh_350[k];

        t_463[k] = f_21 * osh_351[k]
                   + f_3 * pc_x[k] * qsh_351[k];

        t_464[k] = f_21 * osh_352[k]
                   + f_3 * pc_x[k] * qsh_352[k];

        t_465[k] = f_21 * osh_353[k]
                   + f_3 * pc_x[k] * qsh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pc_x, pc_z, osi0_301, osh_354, \
                         osh_355, osh_356, osi1_301, qsh_354, qsh_355, \
                         qsh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_21 * osh_354[k]
                   + f_3 * pc_x[k] * qsh_354[k];

        t_467[k] = f_21 * osh_355[k]
                   + f_3 * pc_x[k] * qsh_355[k];

        t_468[k] = f_21 * osh_356[k]
                   + f_3 * pc_x[k] * qsh_356[k];

        t_469[k] = pa_z[k] * osi0_301[k]
                   - f_10 * pc_z[k] * osi1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, osh_225, osh_248, osh_249, qsg0_252, \
                         qsg0_253, qsg1_252, qsg1_253, qsh_351, qsh_353, \
                         qsh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * osh_225[k]
                   + f_3 * pc_z[k] * qsh_351[k];

        t_471[k] = f_14 * osh_248[k]
                   + f_8 * qsg0_252[k]
                   - f_9 * qsg1_252[k]
                   + f_3 * pc_y[k] * qsh_353[k];

        t_472[k] = f_14 * osh_249[k]
                   + f_6 * qsg0_253[k]
                   - f_7 * qsg1_253[k]
                   + f_3 * pc_y[k] * qsh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, osh_230, osh_250, osh_251, qsg0_254, \
                         qsg1_254, qsh_355, qsh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * osh_250[k]
                   + f_4 * qsg0_254[k]
                   - f_5 * qsg1_254[k]
                   + f_3 * pc_y[k] * qsh_355[k];

        t_474[k] = f_14 * osh_251[k]
                   + f_3 * pc_y[k] * qsh_356[k];

        t_475[k] = f_11 * osh_230[k]
                   + f_1 * qsg0_254[k]
                   - f_2 * qsg1_254[k]
                   + f_3 * pc_z[k] * qsh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, osh_231, osh_252, osh_357, \
                         qsg0_255, qsg1_255, qsh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_21 * osh_357[k]
                   + f_1 * qsg0_255[k]
                   - f_2 * qsg1_255[k]
                   + f_3 * pc_x[k] * qsh_357[k];

        t_477[k] = f_13 * osh_252[k]
                   + f_3 * pc_y[k] * qsh_357[k];

        t_478[k] = f_12 * osh_231[k]
                   + f_3 * pc_z[k] * qsh_357[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_392 = buffer.data(osi0 + 392);
    const auto *osi0_395 = buffer.data(osi0 + 395);
    const auto *osi0_397 = buffer.data(osi0 + 397);
    const auto *osi0_398 = buffer.data(osi0 + 398);
    const auto *osi0_401 = buffer.data(osi0 + 401);
    const auto *osi0_402 = buffer.data(osi0 + 402);
    const auto *osi0_404 = buffer.data(osi0 + 404);
    const auto *osi0_406 = buffer.data(osi0 + 406);
    const auto *osi0_419 = buffer.data(osi0 + 419);

    const auto *osh_234 = buffer.data(osh + 234);
    const auto *osh_237 = buffer.data(osh + 237);
    const auto *osh_246 = buffer.data(osh + 246);
    const auto *osh_251 = buffer.data(osh + 251);
    const auto *osh_252 = buffer.data(osh + 252);
    const auto *osh_254 = buffer.data(osh + 254);
    const auto *osh_255 = buffer.data(osh + 255);
    const auto *osh_257 = buffer.data(osh + 257);
    const auto *osh_258 = buffer.data(osh + 258);
    const auto *osh_261 = buffer.data(osh + 261);
    const auto *osh_267 = buffer.data(osh + 267);
    const auto *osh_269 = buffer.data(osh + 269);
    const auto *osh_270 = buffer.data(osh + 270);
    const auto *osh_271 = buffer.data(osh + 271);
    const auto *osh_272 = buffer.data(osh + 272);
    const auto *osh_273 = buffer.data(osh + 273);
    const auto *osh_275 = buffer.data(osh + 275);
    const auto *osh_276 = buffer.data(osh + 276);
    const auto *osh_278 = buffer.data(osh + 278);
    const auto *osh_279 = buffer.data(osh + 279);
    const auto *osh_282 = buffer.data(osh + 282);
    const auto *osh_288 = buffer.data(osh + 288);
    const auto *osh_290 = buffer.data(osh + 290);
    const auto *osh_291 = buffer.data(osh + 291);
    const auto *osh_292 = buffer.data(osh + 292);
    const auto *osh_293 = buffer.data(osh + 293);
    const auto *osh_294 = buffer.data(osh + 294);
    const auto *osh_295 = buffer.data(osh + 295);
    const auto *osh_296 = buffer.data(osh + 296);
    const auto *osh_297 = buffer.data(osh + 297);
    const auto *osh_299 = buffer.data(osh + 299);
    const auto *osh_300 = buffer.data(osh + 300);
    const auto *osh_302 = buffer.data(osh + 302);
    const auto *osh_303 = buffer.data(osh + 303);
    const auto *osh_309 = buffer.data(osh + 309);
    const auto *osh_311 = buffer.data(osh + 311);
    const auto *osh_312 = buffer.data(osh + 312);
    const auto *osh_313 = buffer.data(osh + 313);
    const auto *osh_314 = buffer.data(osh + 314);
    const auto *osh_315 = buffer.data(osh + 315);
    const auto *osh_360 = buffer.data(osh + 360);
    const auto *osh_362 = buffer.data(osh + 362);
    const auto *osh_363 = buffer.data(osh + 363);
    const auto *osh_366 = buffer.data(osh + 366);
    const auto *osh_367 = buffer.data(osh + 367);
    const auto *osh_369 = buffer.data(osh + 369);
    const auto *osh_371 = buffer.data(osh + 371);
    const auto *osh_372 = buffer.data(osh + 372);
    const auto *osh_373 = buffer.data(osh + 373);
    const auto *osh_374 = buffer.data(osh + 374);
    const auto *osh_375 = buffer.data(osh + 375);
    const auto *osh_376 = buffer.data(osh + 376);
    const auto *osh_377 = buffer.data(osh + 377);
    const auto *osh_378 = buffer.data(osh + 378);
    const auto *osh_381 = buffer.data(osh + 381);
    const auto *osh_383 = buffer.data(osh + 383);
    const auto *osh_384 = buffer.data(osh + 384);
    const auto *osh_387 = buffer.data(osh + 387);
    const auto *osh_388 = buffer.data(osh + 388);
    const auto *osh_390 = buffer.data(osh + 390);
    const auto *osh_392 = buffer.data(osh + 392);
    const auto *osh_393 = buffer.data(osh + 393);
    const auto *osh_394 = buffer.data(osh + 394);
    const auto *osh_395 = buffer.data(osh + 395);
    const auto *osh_396 = buffer.data(osh + 396);
    const auto *osh_397 = buffer.data(osh + 397);
    const auto *osh_398 = buffer.data(osh + 398);
    const auto *osh_414 = buffer.data(osh + 414);
    const auto *osh_415 = buffer.data(osh + 415);
    const auto *osh_416 = buffer.data(osh + 416);
    const auto *osh_417 = buffer.data(osh + 417);
    const auto *osh_418 = buffer.data(osh + 418);
    const auto *osh_419 = buffer.data(osh + 419);
    const auto *osh_420 = buffer.data(osh + 420);
    const auto *osh_425 = buffer.data(osh + 425);
    const auto *osh_429 = buffer.data(osh + 429);
    const auto *osh_434 = buffer.data(osh + 434);
    const auto *osh_435 = buffer.data(osh + 435);
    const auto *osh_436 = buffer.data(osh + 436);
    const auto *osh_437 = buffer.data(osh + 437);
    const auto *osh_438 = buffer.data(osh + 438);
    const auto *osh_440 = buffer.data(osh + 440);
    const auto *osh_441 = buffer.data(osh + 441);
    const auto *osh_444 = buffer.data(osh + 444);

    const auto *osi1_392 = buffer.data(osi1 + 392);
    const auto *osi1_395 = buffer.data(osi1 + 395);
    const auto *osi1_397 = buffer.data(osi1 + 397);
    const auto *osi1_398 = buffer.data(osi1 + 398);
    const auto *osi1_401 = buffer.data(osi1 + 401);
    const auto *osi1_402 = buffer.data(osi1 + 402);
    const auto *osi1_404 = buffer.data(osi1 + 404);
    const auto *osi1_406 = buffer.data(osi1 + 406);
    const auto *osi1_419 = buffer.data(osi1 + 419);

    const auto *qsg0_258 = buffer.data(qsg0 + 258);
    const auto *qsg0_260 = buffer.data(qsg0 + 260);
    const auto *qsg0_261 = buffer.data(qsg0 + 261);
    const auto *qsg0_264 = buffer.data(qsg0 + 264);
    const auto *qsg0_265 = buffer.data(qsg0 + 265);
    const auto *qsg0_267 = buffer.data(qsg0 + 267);
    const auto *qsg0_268 = buffer.data(qsg0 + 268);
    const auto *qsg0_269 = buffer.data(qsg0 + 269);
    const auto *qsg0_270 = buffer.data(qsg0 + 270);
    const auto *qsg0_273 = buffer.data(qsg0 + 273);
    const auto *qsg0_275 = buffer.data(qsg0 + 275);
    const auto *qsg0_276 = buffer.data(qsg0 + 276);
    const auto *qsg0_279 = buffer.data(qsg0 + 279);
    const auto *qsg0_280 = buffer.data(qsg0 + 280);
    const auto *qsg0_282 = buffer.data(qsg0 + 282);
    const auto *qsg0_283 = buffer.data(qsg0 + 283);
    const auto *qsg0_284 = buffer.data(qsg0 + 284);
    const auto *qsg0_295 = buffer.data(qsg0 + 295);
    const auto *qsg0_297 = buffer.data(qsg0 + 297);
    const auto *qsg0_298 = buffer.data(qsg0 + 298);
    const auto *qsg0_299 = buffer.data(qsg0 + 299);
    const auto *qsg0_300 = buffer.data(qsg0 + 300);
    const auto *qsg0_301 = buffer.data(qsg0 + 301);
    const auto *qsg0_302 = buffer.data(qsg0 + 302);
    const auto *qsg0_303 = buffer.data(qsg0 + 303);
    const auto *qsg0_304 = buffer.data(qsg0 + 304);
    const auto *qsg0_305 = buffer.data(qsg0 + 305);
    const auto *qsg0_309 = buffer.data(qsg0 + 309);
    const auto *qsg0_310 = buffer.data(qsg0 + 310);
    const auto *qsg0_311 = buffer.data(qsg0 + 311);
    const auto *qsg0_312 = buffer.data(qsg0 + 312);
    const auto *qsg0_313 = buffer.data(qsg0 + 313);
    const auto *qsg0_314 = buffer.data(qsg0 + 314);
    const auto *qsg0_315 = buffer.data(qsg0 + 315);
    const auto *qsg0_318 = buffer.data(qsg0 + 318);

    const auto *qsg1_258 = buffer.data(qsg1 + 258);
    const auto *qsg1_260 = buffer.data(qsg1 + 260);
    const auto *qsg1_261 = buffer.data(qsg1 + 261);
    const auto *qsg1_264 = buffer.data(qsg1 + 264);
    const auto *qsg1_265 = buffer.data(qsg1 + 265);
    const auto *qsg1_267 = buffer.data(qsg1 + 267);
    const auto *qsg1_268 = buffer.data(qsg1 + 268);
    const auto *qsg1_269 = buffer.data(qsg1 + 269);
    const auto *qsg1_270 = buffer.data(qsg1 + 270);
    const auto *qsg1_273 = buffer.data(qsg1 + 273);
    const auto *qsg1_275 = buffer.data(qsg1 + 275);
    const auto *qsg1_276 = buffer.data(qsg1 + 276);
    const auto *qsg1_279 = buffer.data(qsg1 + 279);
    const auto *qsg1_280 = buffer.data(qsg1 + 280);
    const auto *qsg1_282 = buffer.data(qsg1 + 282);
    const auto *qsg1_283 = buffer.data(qsg1 + 283);
    const auto *qsg1_284 = buffer.data(qsg1 + 284);
    const auto *qsg1_295 = buffer.data(qsg1 + 295);
    const auto *qsg1_297 = buffer.data(qsg1 + 297);
    const auto *qsg1_298 = buffer.data(qsg1 + 298);
    const auto *qsg1_299 = buffer.data(qsg1 + 299);
    const auto *qsg1_300 = buffer.data(qsg1 + 300);
    const auto *qsg1_301 = buffer.data(qsg1 + 301);
    const auto *qsg1_302 = buffer.data(qsg1 + 302);
    const auto *qsg1_303 = buffer.data(qsg1 + 303);
    const auto *qsg1_304 = buffer.data(qsg1 + 304);
    const auto *qsg1_305 = buffer.data(qsg1 + 305);
    const auto *qsg1_309 = buffer.data(qsg1 + 309);
    const auto *qsg1_310 = buffer.data(qsg1 + 310);
    const auto *qsg1_311 = buffer.data(qsg1 + 311);
    const auto *qsg1_312 = buffer.data(qsg1 + 312);
    const auto *qsg1_313 = buffer.data(qsg1 + 313);
    const auto *qsg1_314 = buffer.data(qsg1 + 314);
    const auto *qsg1_315 = buffer.data(qsg1 + 315);
    const auto *qsg1_318 = buffer.data(qsg1 + 318);

    const auto *qsh_359 = buffer.data(qsh + 359);
    const auto *qsh_360 = buffer.data(qsh + 360);
    const auto *qsh_362 = buffer.data(qsh + 362);
    const auto *qsh_363 = buffer.data(qsh + 363);
    const auto *qsh_366 = buffer.data(qsh + 366);
    const auto *qsh_367 = buffer.data(qsh + 367);
    const auto *qsh_369 = buffer.data(qsh + 369);
    const auto *qsh_371 = buffer.data(qsh + 371);
    const auto *qsh_372 = buffer.data(qsh + 372);
    const auto *qsh_373 = buffer.data(qsh + 373);
    const auto *qsh_374 = buffer.data(qsh + 374);
    const auto *qsh_375 = buffer.data(qsh + 375);
    const auto *qsh_376 = buffer.data(qsh + 376);
    const auto *qsh_377 = buffer.data(qsh + 377);
    const auto *qsh_378 = buffer.data(qsh + 378);
    const auto *qsh_380 = buffer.data(qsh + 380);
    const auto *qsh_381 = buffer.data(qsh + 381);
    const auto *qsh_383 = buffer.data(qsh + 383);
    const auto *qsh_384 = buffer.data(qsh + 384);
    const auto *qsh_387 = buffer.data(qsh + 387);
    const auto *qsh_388 = buffer.data(qsh + 388);
    const auto *qsh_390 = buffer.data(qsh + 390);
    const auto *qsh_392 = buffer.data(qsh + 392);
    const auto *qsh_393 = buffer.data(qsh + 393);
    const auto *qsh_394 = buffer.data(qsh + 394);
    const auto *qsh_395 = buffer.data(qsh + 395);
    const auto *qsh_396 = buffer.data(qsh + 396);
    const auto *qsh_397 = buffer.data(qsh + 397);
    const auto *qsh_398 = buffer.data(qsh + 398);
    const auto *qsh_399 = buffer.data(qsh + 399);
    const auto *qsh_401 = buffer.data(qsh + 401);
    const auto *qsh_402 = buffer.data(qsh + 402);
    const auto *qsh_404 = buffer.data(qsh + 404);
    const auto *qsh_405 = buffer.data(qsh + 405);
    const auto *qsh_408 = buffer.data(qsh + 408);
    const auto *qsh_414 = buffer.data(qsh + 414);
    const auto *qsh_415 = buffer.data(qsh + 415);
    const auto *qsh_416 = buffer.data(qsh + 416);
    const auto *qsh_417 = buffer.data(qsh + 417);
    const auto *qsh_418 = buffer.data(qsh + 418);
    const auto *qsh_419 = buffer.data(qsh + 419);
    const auto *qsh_420 = buffer.data(qsh + 420);
    const auto *qsh_421 = buffer.data(qsh + 421);
    const auto *qsh_422 = buffer.data(qsh + 422);
    const auto *qsh_423 = buffer.data(qsh + 423);
    const auto *qsh_424 = buffer.data(qsh + 424);
    const auto *qsh_425 = buffer.data(qsh + 425);
    const auto *qsh_426 = buffer.data(qsh + 426);
    const auto *qsh_427 = buffer.data(qsh + 427);
    const auto *qsh_428 = buffer.data(qsh + 428);
    const auto *qsh_429 = buffer.data(qsh + 429);
    const auto *qsh_434 = buffer.data(qsh + 434);
    const auto *qsh_435 = buffer.data(qsh + 435);
    const auto *qsh_436 = buffer.data(qsh + 436);
    const auto *qsh_437 = buffer.data(qsh + 437);
    const auto *qsh_438 = buffer.data(qsh + 438);
    const auto *qsh_439 = buffer.data(qsh + 439);
    const auto *qsh_440 = buffer.data(qsh + 440);
    const auto *qsh_441 = buffer.data(qsh + 441);
    const auto *qsh_444 = buffer.data(qsh + 444);

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, osh_254, osh_360, osh_362, qsg0_258, \
                         qsg0_260, qsg1_258, qsg1_260, qsh_359, qsh_360, \
                         qsh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_21 * osh_360[k]
                   + f_8 * qsg0_258[k]
                   - f_9 * qsg1_258[k]
                   + f_3 * pc_x[k] * qsh_360[k];

        t_480[k] = f_13 * osh_254[k]
                   + f_3 * pc_y[k] * qsh_359[k];

        t_481[k] = f_21 * osh_362[k]
                   + f_8 * qsg0_260[k]
                   - f_9 * qsg1_260[k]
                   + f_3 * pc_x[k] * qsh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, osh_234, osh_257, osh_363, \
                         qsg0_261, qsg1_261, qsh_360, qsh_362, \
                         qsh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_21 * osh_363[k]
                   + f_6 * qsg0_261[k]
                   - f_7 * qsg1_261[k]
                   + f_3 * pc_x[k] * qsh_363[k];

        t_483[k] = f_12 * osh_234[k]
                   + f_3 * pc_z[k] * qsh_360[k];

        t_484[k] = f_13 * osh_257[k]
                   + f_3 * pc_y[k] * qsh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, osh_237, osh_366, osh_367, qsg0_264, \
                         qsg0_265, qsg1_264, qsg1_265, qsh_363, qsh_366, \
                         qsh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_21 * osh_366[k]
                   + f_6 * qsg0_264[k]
                   - f_7 * qsg1_264[k]
                   + f_3 * pc_x[k] * qsh_366[k];

        t_486[k] = f_21 * osh_367[k]
                   + f_4 * qsg0_265[k]
                   - f_5 * qsg1_265[k]
                   + f_3 * pc_x[k] * qsh_367[k];

        t_487[k] = f_12 * osh_237[k]
                   + f_3 * pc_z[k] * qsh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, osh_261, osh_369, osh_371, qsg0_267, \
                         qsg0_269, qsg1_267, qsg1_269, qsh_366, qsh_369, \
                         qsh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_21 * osh_369[k]
                   + f_4 * qsg0_267[k]
                   - f_5 * qsg1_267[k]
                   + f_3 * pc_x[k] * qsh_369[k];

        t_489[k] = f_13 * osh_261[k]
                   + f_3 * pc_y[k] * qsh_366[k];

        t_490[k] = f_21 * osh_371[k]
                   + f_4 * qsg0_269[k]
                   - f_5 * qsg1_269[k]
                   + f_3 * pc_x[k] * qsh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, osh_372, osh_373, osh_374, \
                         osh_375, osh_376, qsh_372, qsh_373, qsh_374, qsh_375, \
                         qsh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_21 * osh_372[k]
                   + f_3 * pc_x[k] * qsh_372[k];

        t_492[k] = f_21 * osh_373[k]
                   + f_3 * pc_x[k] * qsh_373[k];

        t_493[k] = f_21 * osh_374[k]
                   + f_3 * pc_x[k] * qsh_374[k];

        t_494[k] = f_21 * osh_375[k]
                   + f_3 * pc_x[k] * qsh_375[k];

        t_495[k] = f_21 * osh_376[k]
                   + f_3 * pc_x[k] * qsh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, osh_246, osh_267, osh_377, \
                         qsg0_265, qsg1_265, qsh_372, qsh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_21 * osh_377[k]
                   + f_3 * pc_x[k] * qsh_377[k];

        t_497[k] = f_13 * osh_267[k]
                   + f_1 * qsg0_265[k]
                   - f_2 * qsg1_265[k]
                   + f_3 * pc_y[k] * qsh_372[k];

        t_498[k] = f_12 * osh_246[k]
                   + f_3 * pc_z[k] * qsh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, osh_269, osh_270, osh_271, qsg0_267, \
                         qsg0_268, qsg0_269, qsg1_267, qsg1_268, qsg1_269, qsh_374, qsh_375, \
                         qsh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * osh_269[k]
                   + f_8 * qsg0_267[k]
                   - f_9 * qsg1_267[k]
                   + f_3 * pc_y[k] * qsh_374[k];

        t_500[k] = f_13 * osh_270[k]
                   + f_6 * qsg0_268[k]
                   - f_7 * qsg1_268[k]
                   + f_3 * pc_y[k] * qsh_375[k];

        t_501[k] = f_13 * osh_271[k]
                   + f_4 * qsg0_269[k]
                   - f_5 * qsg1_269[k]
                   + f_3 * pc_y[k] * qsh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, osh_251, osh_272, osh_378, \
                         qsg0_269, qsg0_270, qsg1_269, qsg1_270, qsh_377, \
                         qsh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * osh_272[k]
                   + f_3 * pc_y[k] * qsh_377[k];

        t_503[k] = f_12 * osh_251[k]
                   + f_1 * qsg0_269[k]
                   - f_2 * qsg1_269[k]
                   + f_3 * pc_z[k] * qsh_377[k];

        t_504[k] = f_21 * osh_378[k]
                   + f_1 * qsg0_270[k]
                   - f_2 * qsg1_270[k]
                   + f_3 * pc_x[k] * qsh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, osh_252, osh_273, \
                         osh_275, osh_381, qsg0_273, qsg1_273, qsh_378, qsh_380, \
                         qsh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * osh_273[k]
                   + f_3 * pc_y[k] * qsh_378[k];

        t_506[k] = f_13 * osh_252[k]
                   + f_3 * pc_z[k] * qsh_378[k];

        t_507[k] = f_21 * osh_381[k]
                   + f_8 * qsg0_273[k]
                   - f_9 * qsg1_273[k]
                   + f_3 * pc_x[k] * qsh_381[k];

        t_508[k] = f_12 * osh_275[k]
                   + f_3 * pc_y[k] * qsh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, osh_255, osh_383, osh_384, qsg0_275, \
                         qsg0_276, qsg1_275, qsg1_276, qsh_381, qsh_383, \
                         qsh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_21 * osh_383[k]
                   + f_8 * qsg0_275[k]
                   - f_9 * qsg1_275[k]
                   + f_3 * pc_x[k] * qsh_383[k];

        t_510[k] = f_21 * osh_384[k]
                   + f_6 * qsg0_276[k]
                   - f_7 * qsg1_276[k]
                   + f_3 * pc_x[k] * qsh_384[k];

        t_511[k] = f_13 * osh_255[k]
                   + f_3 * pc_z[k] * qsh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, osh_278, osh_387, osh_388, qsg0_279, \
                         qsg0_280, qsg1_279, qsg1_280, qsh_383, qsh_387, \
                         qsh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * osh_278[k]
                   + f_3 * pc_y[k] * qsh_383[k];

        t_513[k] = f_21 * osh_387[k]
                   + f_6 * qsg0_279[k]
                   - f_7 * qsg1_279[k]
                   + f_3 * pc_x[k] * qsh_387[k];

        t_514[k] = f_21 * osh_388[k]
                   + f_4 * qsg0_280[k]
                   - f_5 * qsg1_280[k]
                   + f_3 * pc_x[k] * qsh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, osh_258, osh_282, osh_390, \
                         qsg0_282, qsg1_282, qsh_384, qsh_387, \
                         qsh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * osh_258[k]
                   + f_3 * pc_z[k] * qsh_384[k];

        t_516[k] = f_21 * osh_390[k]
                   + f_4 * qsg0_282[k]
                   - f_5 * qsg1_282[k]
                   + f_3 * pc_x[k] * qsh_390[k];

        t_517[k] = f_12 * osh_282[k]
                   + f_3 * pc_y[k] * qsh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, osh_392, osh_393, osh_394, osh_395, \
                         qsg0_284, qsg1_284, qsh_392, qsh_393, qsh_394, \
                         qsh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_21 * osh_392[k]
                   + f_4 * qsg0_284[k]
                   - f_5 * qsg1_284[k]
                   + f_3 * pc_x[k] * qsh_392[k];

        t_519[k] = f_21 * osh_393[k]
                   + f_3 * pc_x[k] * qsh_393[k];

        t_520[k] = f_21 * osh_394[k]
                   + f_3 * pc_x[k] * qsh_394[k];

        t_521[k] = f_21 * osh_395[k]
                   + f_3 * pc_x[k] * qsh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, osh_288, osh_396, osh_397, \
                         osh_398, qsg0_280, qsg1_280, qsh_393, qsh_396, qsh_397, \
                         qsh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_21 * osh_396[k]
                   + f_3 * pc_x[k] * qsh_396[k];

        t_523[k] = f_21 * osh_397[k]
                   + f_3 * pc_x[k] * qsh_397[k];

        t_524[k] = f_21 * osh_398[k]
                   + f_3 * pc_x[k] * qsh_398[k];

        t_525[k] = f_12 * osh_288[k]
                   + f_1 * qsg0_280[k]
                   - f_2 * qsg1_280[k]
                   + f_3 * pc_y[k] * qsh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, osh_267, osh_290, osh_291, qsg0_282, \
                         qsg0_283, qsg1_282, qsg1_283, qsh_393, qsh_395, \
                         qsh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * osh_267[k]
                   + f_3 * pc_z[k] * qsh_393[k];

        t_527[k] = f_12 * osh_290[k]
                   + f_8 * qsg0_282[k]
                   - f_9 * qsg1_282[k]
                   + f_3 * pc_y[k] * qsh_395[k];

        t_528[k] = f_12 * osh_291[k]
                   + f_6 * qsg0_283[k]
                   - f_7 * qsg1_283[k]
                   + f_3 * pc_y[k] * qsh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pa_y, pc_y, pc_z, osi0_392, osh_272, \
                         osh_292, osh_293, osi1_392, qsg0_284, qsg1_284, qsh_397, \
                         qsh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * osh_292[k]
                   + f_4 * qsg0_284[k]
                   - f_5 * qsg1_284[k]
                   + f_3 * pc_y[k] * qsh_397[k];

        t_530[k] = f_12 * osh_293[k]
                   + f_3 * pc_y[k] * qsh_398[k];

        t_531[k] = f_13 * osh_272[k]
                   + f_1 * qsg0_284[k]
                   - f_2 * qsg1_284[k]
                   + f_3 * pc_z[k] * qsh_398[k];

        t_532[k] = pa_y[k] * osi0_392[k]
                   - f_10 * pc_y[k] * osi1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_y, pc_y, pc_z, osi0_395, osh_273, \
                         osh_294, osh_295, osh_296, osi1_395, qsh_399, \
                         qsh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * osh_294[k]
                   + f_3 * pc_y[k] * qsh_399[k];

        t_534[k] = f_14 * osh_273[k]
                   + f_3 * pc_z[k] * qsh_399[k];

        t_535[k] = pa_y[k] * osi0_395[k]
                   + f_12 * osh_295[k]
                   - f_10 * pc_y[k] * osi1_395[k];

        t_536[k] = f_11 * osh_296[k]
                   + f_3 * pc_y[k] * qsh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_y, pc_y, pc_z, osi0_397, osi0_398, \
                         osh_276, osh_297, osh_299, osi1_397, osi1_398, qsh_402, \
                         qsh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pa_y[k] * osi0_397[k]
                   - f_10 * pc_y[k] * osi1_397[k];

        t_538[k] = pa_y[k] * osi0_398[k]
                   + f_13 * osh_297[k]
                   - f_10 * pc_y[k] * osi1_398[k];

        t_539[k] = f_14 * osh_276[k]
                   + f_3 * pc_z[k] * qsh_402[k];

        t_540[k] = f_11 * osh_299[k]
                   + f_3 * pc_y[k] * qsh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pa_y, pc_y, pc_z, osi0_401, osi0_402, osh_279, \
                         osh_300, osi1_401, osi1_402, qsh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pa_y[k] * osi0_401[k]
                   - f_10 * pc_y[k] * osi1_401[k];

        t_542[k] = pa_y[k] * osi0_402[k]
                   + f_14 * osh_300[k]
                   - f_10 * pc_y[k] * osi1_402[k];

        t_543[k] = f_14 * osh_279[k]
                   + f_3 * pc_z[k] * qsh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_x, pc_y, osi0_404, osi0_406, \
                         osh_302, osh_303, osh_414, osi1_404, osi1_406, qsh_408, \
                         qsh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pa_y[k] * osi0_404[k]
                   + f_12 * osh_302[k]
                   - f_10 * pc_y[k] * osi1_404[k];

        t_545[k] = f_11 * osh_303[k]
                   + f_3 * pc_y[k] * qsh_408[k];

        t_546[k] = pa_y[k] * osi0_406[k]
                   - f_10 * pc_y[k] * osi1_406[k];

        t_547[k] = f_21 * osh_414[k]
                   + f_3 * pc_x[k] * qsh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, osh_415, osh_416, osh_417, \
                         osh_418, osh_419, qsh_415, qsh_416, qsh_417, qsh_418, \
                         qsh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_21 * osh_415[k]
                   + f_3 * pc_x[k] * qsh_415[k];

        t_549[k] = f_21 * osh_416[k]
                   + f_3 * pc_x[k] * qsh_416[k];

        t_550[k] = f_21 * osh_417[k]
                   + f_3 * pc_x[k] * qsh_417[k];

        t_551[k] = f_21 * osh_418[k]
                   + f_3 * pc_x[k] * qsh_418[k];

        t_552[k] = f_21 * osh_419[k]
                   + f_3 * pc_x[k] * qsh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, osh_288, osh_309, osh_311, qsg0_295, \
                         qsg0_297, qsg1_295, qsg1_297, qsh_414, \
                         qsh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * osh_309[k]
                   + f_1 * qsg0_295[k]
                   - f_2 * qsg1_295[k]
                   + f_3 * pc_y[k] * qsh_414[k];

        t_554[k] = f_14 * osh_288[k]
                   + f_3 * pc_z[k] * qsh_414[k];

        t_555[k] = f_11 * osh_311[k]
                   + f_8 * qsg0_297[k]
                   - f_9 * qsg1_297[k]
                   + f_3 * pc_y[k] * qsh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, osh_312, osh_313, osh_314, qsg0_298, \
                         qsg0_299, qsg1_298, qsg1_299, qsh_417, qsh_418, \
                         qsh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * osh_312[k]
                   + f_6 * qsg0_298[k]
                   - f_7 * qsg1_298[k]
                   + f_3 * pc_y[k] * qsh_417[k];

        t_557[k] = f_11 * osh_313[k]
                   + f_4 * qsg0_299[k]
                   - f_5 * qsg1_299[k]
                   + f_3 * pc_y[k] * qsh_418[k];

        t_558[k] = f_11 * osh_314[k]
                   + f_3 * pc_y[k] * qsh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_y, pc_x, pc_y, pc_z, osi0_419, \
                         osh_294, osh_420, osi1_419, qsg0_300, qsg1_300, \
                         qsh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pa_y[k] * osi0_419[k]
                   - f_10 * pc_y[k] * osi1_419[k];

        t_560[k] = f_21 * osh_420[k]
                   + f_1 * qsg0_300[k]
                   - f_2 * qsg1_300[k]
                   + f_3 * pc_x[k] * qsh_420[k];

        t_561[k] = f_3 * pc_y[k] * qsh_420[k];

        t_562[k] = f_22 * osh_294[k]
                   + f_3 * pc_z[k] * qsh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, osh_425, qsg0_300, qsg0_305, \
                         qsg1_300, qsg1_305, qsh_421, qsh_422, \
                         qsh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_4 * qsg0_300[k]
                   - f_5 * qsg1_300[k]
                   + f_3 * pc_y[k] * qsh_421[k];

        t_564[k] = f_3 * pc_y[k] * qsh_422[k];

        t_565[k] = f_21 * osh_425[k]
                   + f_8 * qsg0_305[k]
                   - f_9 * qsg1_305[k]
                   + f_3 * pc_x[k] * qsh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_y, qsg0_301, qsg0_302, qsg1_301, qsg1_302, \
                         qsh_423, qsh_424, qsh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_6 * qsg0_301[k]
                   - f_7 * qsg1_301[k]
                   + f_3 * pc_y[k] * qsh_423[k];

        t_567[k] = f_4 * qsg0_302[k]
                   - f_5 * qsg1_302[k]
                   + f_3 * pc_y[k] * qsh_424[k];

        t_568[k] = f_3 * pc_y[k] * qsh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_y, osh_429, qsg0_303, qsg0_304, \
                         qsg0_309, qsg1_303, qsg1_304, qsg1_309, qsh_426, qsh_427, \
                         qsh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_21 * osh_429[k]
                   + f_6 * qsg0_309[k]
                   - f_7 * qsg1_309[k]
                   + f_3 * pc_x[k] * qsh_429[k];

        t_570[k] = f_8 * qsg0_303[k]
                   - f_9 * qsg1_303[k]
                   + f_3 * pc_y[k] * qsh_426[k];

        t_571[k] = f_6 * qsg0_304[k]
                   - f_7 * qsg1_304[k]
                   + f_3 * pc_y[k] * qsh_427[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, osh_434, osh_435, qsg0_305, \
                         qsg0_314, qsg1_305, qsg1_314, qsh_428, qsh_429, qsh_434, \
                         qsh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_4 * qsg0_305[k]
                   - f_5 * qsg1_305[k]
                   + f_3 * pc_y[k] * qsh_428[k];

        t_573[k] = f_3 * pc_y[k] * qsh_429[k];

        t_574[k] = f_21 * osh_434[k]
                   + f_4 * qsg0_314[k]
                   - f_5 * qsg1_314[k]
                   + f_3 * pc_x[k] * qsh_434[k];

        t_575[k] = f_21 * osh_435[k]
                   + f_3 * pc_x[k] * qsh_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, pc_x, pc_y, osh_436, osh_437, \
                         osh_438, osh_440, qsh_434, qsh_436, qsh_437, qsh_438, \
                         qsh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_21 * osh_436[k]
                   + f_3 * pc_x[k] * qsh_436[k];

        t_577[k] = f_21 * osh_437[k]
                   + f_3 * pc_x[k] * qsh_437[k];

        t_578[k] = f_21 * osh_438[k]
                   + f_3 * pc_x[k] * qsh_438[k];

        t_579[k] = f_3 * pc_y[k] * qsh_434[k];

        t_580[k] = f_21 * osh_440[k]
                   + f_3 * pc_x[k] * qsh_440[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, qsg0_310, qsg0_311, qsg0_312, qsg1_310, \
                         qsg1_311, qsg1_312, qsh_435, qsh_436, \
                         qsh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_1 * qsg0_310[k]
                   - f_2 * qsg1_310[k]
                   + f_3 * pc_y[k] * qsh_435[k];

        t_582[k] = f_16 * qsg0_311[k]
                   - f_17 * qsg1_311[k]
                   + f_3 * pc_y[k] * qsh_436[k];

        t_583[k] = f_8 * qsg0_312[k]
                   - f_9 * qsg1_312[k]
                   + f_3 * pc_y[k] * qsh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, osh_314, qsg0_313, qsg0_314, \
                         qsg1_313, qsg1_314, qsh_438, qsh_439, \
                         qsh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * qsg0_313[k]
                   - f_7 * qsg1_313[k]
                   + f_3 * pc_y[k] * qsh_438[k];

        t_585[k] = f_4 * qsg0_314[k]
                   - f_5 * qsg1_314[k]
                   + f_3 * pc_y[k] * qsh_439[k];

        t_586[k] = f_3 * pc_y[k] * qsh_440[k];

        t_587[k] = f_22 * osh_314[k]
                   + f_1 * qsg0_314[k]
                   - f_2 * qsg1_314[k]
                   + f_3 * pc_z[k] * qsh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, osh_315, osh_441, \
                         osh_444, qsg0_315, qsg0_318, qsg1_315, qsg1_318, qsh_441, \
                         qsh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_23 * osh_441[k]
                   + f_1 * qsg0_315[k]
                   - f_2 * qsg1_315[k]
                   + f_3 * pc_x[k] * qsh_441[k];

        t_589[k] = f_23 * osh_315[k]
                   + f_3 * pc_y[k] * qsh_441[k];

        t_590[k] = f_3 * pc_z[k] * qsh_441[k];

        t_591[k] = f_23 * osh_444[k]
                   + f_8 * qsg0_318[k]
                   - f_9 * qsg1_318[k]
                   + f_3 * pc_x[k] * qsh_444[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

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
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_420 = buffer.data(osi0 + 420);
    const auto *osi0_423 = buffer.data(osi0 + 423);
    const auto *osi0_426 = buffer.data(osi0 + 426);
    const auto *osi0_430 = buffer.data(osi0 + 430);
    const auto *osi0_432 = buffer.data(osi0 + 432);
    const auto *osi0_441 = buffer.data(osi0 + 441);

    const auto *osh_315 = buffer.data(osh + 315);
    const auto *osh_318 = buffer.data(osh + 318);
    const auto *osh_320 = buffer.data(osh + 320);
    const auto *osh_321 = buffer.data(osh + 321);
    const auto *osh_322 = buffer.data(osh + 322);
    const auto *osh_324 = buffer.data(osh + 324);
    const auto *osh_330 = buffer.data(osh + 330);
    const auto *osh_335 = buffer.data(osh + 335);
    const auto *osh_336 = buffer.data(osh + 336);
    const auto *osh_338 = buffer.data(osh + 338);
    const auto *osh_339 = buffer.data(osh + 339);
    const auto *osh_341 = buffer.data(osh + 341);
    const auto *osh_342 = buffer.data(osh + 342);
    const auto *osh_345 = buffer.data(osh + 345);
    const auto *osh_351 = buffer.data(osh + 351);
    const auto *osh_353 = buffer.data(osh + 353);
    const auto *osh_354 = buffer.data(osh + 354);
    const auto *osh_355 = buffer.data(osh + 355);
    const auto *osh_356 = buffer.data(osh + 356);
    const auto *osh_357 = buffer.data(osh + 357);
    const auto *osh_359 = buffer.data(osh + 359);
    const auto *osh_360 = buffer.data(osh + 360);
    const auto *osh_362 = buffer.data(osh + 362);
    const auto *osh_363 = buffer.data(osh + 363);
    const auto *osh_366 = buffer.data(osh + 366);
    const auto *osh_372 = buffer.data(osh + 372);
    const auto *osh_374 = buffer.data(osh + 374);
    const auto *osh_375 = buffer.data(osh + 375);
    const auto *osh_376 = buffer.data(osh + 376);
    const auto *osh_377 = buffer.data(osh + 377);
    const auto *osh_378 = buffer.data(osh + 378);
    const auto *osh_380 = buffer.data(osh + 380);
    const auto *osh_383 = buffer.data(osh + 383);
    const auto *osh_387 = buffer.data(osh + 387);
    const auto *osh_393 = buffer.data(osh + 393);
    const auto *osh_395 = buffer.data(osh + 395);
    const auto *osh_396 = buffer.data(osh + 396);
    const auto *osh_397 = buffer.data(osh + 397);
    const auto *osh_398 = buffer.data(osh + 398);
    const auto *osh_399 = buffer.data(osh + 399);
    const auto *osh_447 = buffer.data(osh + 447);
    const auto *osh_451 = buffer.data(osh + 451);
    const auto *osh_456 = buffer.data(osh + 456);
    const auto *osh_458 = buffer.data(osh + 458);
    const auto *osh_459 = buffer.data(osh + 459);
    const auto *osh_460 = buffer.data(osh + 460);
    const auto *osh_461 = buffer.data(osh + 461);
    const auto *osh_467 = buffer.data(osh + 467);
    const auto *osh_471 = buffer.data(osh + 471);
    const auto *osh_476 = buffer.data(osh + 476);
    const auto *osh_477 = buffer.data(osh + 477);
    const auto *osh_478 = buffer.data(osh + 478);
    const auto *osh_479 = buffer.data(osh + 479);
    const auto *osh_480 = buffer.data(osh + 480);
    const auto *osh_481 = buffer.data(osh + 481);
    const auto *osh_482 = buffer.data(osh + 482);
    const auto *osh_483 = buffer.data(osh + 483);
    const auto *osh_486 = buffer.data(osh + 486);
    const auto *osh_488 = buffer.data(osh + 488);
    const auto *osh_489 = buffer.data(osh + 489);
    const auto *osh_492 = buffer.data(osh + 492);
    const auto *osh_493 = buffer.data(osh + 493);
    const auto *osh_495 = buffer.data(osh + 495);
    const auto *osh_497 = buffer.data(osh + 497);
    const auto *osh_498 = buffer.data(osh + 498);
    const auto *osh_499 = buffer.data(osh + 499);
    const auto *osh_500 = buffer.data(osh + 500);
    const auto *osh_501 = buffer.data(osh + 501);
    const auto *osh_502 = buffer.data(osh + 502);
    const auto *osh_503 = buffer.data(osh + 503);
    const auto *osh_504 = buffer.data(osh + 504);
    const auto *osh_507 = buffer.data(osh + 507);
    const auto *osh_509 = buffer.data(osh + 509);
    const auto *osh_510 = buffer.data(osh + 510);
    const auto *osh_513 = buffer.data(osh + 513);
    const auto *osh_514 = buffer.data(osh + 514);
    const auto *osh_516 = buffer.data(osh + 516);
    const auto *osh_518 = buffer.data(osh + 518);
    const auto *osh_519 = buffer.data(osh + 519);
    const auto *osh_520 = buffer.data(osh + 520);
    const auto *osh_521 = buffer.data(osh + 521);
    const auto *osh_522 = buffer.data(osh + 522);
    const auto *osh_523 = buffer.data(osh + 523);
    const auto *osh_524 = buffer.data(osh + 524);
    const auto *osh_525 = buffer.data(osh + 525);

    const auto *osi1_420 = buffer.data(osi1 + 420);
    const auto *osi1_423 = buffer.data(osi1 + 423);
    const auto *osi1_426 = buffer.data(osi1 + 426);
    const auto *osi1_430 = buffer.data(osi1 + 430);
    const auto *osi1_432 = buffer.data(osi1 + 432);
    const auto *osi1_441 = buffer.data(osi1 + 441);

    const auto *qsg0_315 = buffer.data(qsg0 + 315);
    const auto *qsg0_317 = buffer.data(qsg0 + 317);
    const auto *qsg0_318 = buffer.data(qsg0 + 318);
    const auto *qsg0_320 = buffer.data(qsg0 + 320);
    const auto *qsg0_321 = buffer.data(qsg0 + 321);
    const auto *qsg0_325 = buffer.data(qsg0 + 325);
    const auto *qsg0_326 = buffer.data(qsg0 + 326);
    const auto *qsg0_327 = buffer.data(qsg0 + 327);
    const auto *qsg0_329 = buffer.data(qsg0 + 329);
    const auto *qsg0_335 = buffer.data(qsg0 + 335);
    const auto *qsg0_339 = buffer.data(qsg0 + 339);
    const auto *qsg0_342 = buffer.data(qsg0 + 342);
    const auto *qsg0_343 = buffer.data(qsg0 + 343);
    const auto *qsg0_344 = buffer.data(qsg0 + 344);
    const auto *qsg0_345 = buffer.data(qsg0 + 345);
    const auto *qsg0_348 = buffer.data(qsg0 + 348);
    const auto *qsg0_350 = buffer.data(qsg0 + 350);
    const auto *qsg0_351 = buffer.data(qsg0 + 351);
    const auto *qsg0_354 = buffer.data(qsg0 + 354);
    const auto *qsg0_355 = buffer.data(qsg0 + 355);
    const auto *qsg0_357 = buffer.data(qsg0 + 357);
    const auto *qsg0_358 = buffer.data(qsg0 + 358);
    const auto *qsg0_359 = buffer.data(qsg0 + 359);
    const auto *qsg0_360 = buffer.data(qsg0 + 360);
    const auto *qsg0_363 = buffer.data(qsg0 + 363);
    const auto *qsg0_365 = buffer.data(qsg0 + 365);
    const auto *qsg0_366 = buffer.data(qsg0 + 366);
    const auto *qsg0_369 = buffer.data(qsg0 + 369);
    const auto *qsg0_370 = buffer.data(qsg0 + 370);
    const auto *qsg0_372 = buffer.data(qsg0 + 372);
    const auto *qsg0_373 = buffer.data(qsg0 + 373);
    const auto *qsg0_374 = buffer.data(qsg0 + 374);
    const auto *qsg0_375 = buffer.data(qsg0 + 375);

    const auto *qsg1_315 = buffer.data(qsg1 + 315);
    const auto *qsg1_317 = buffer.data(qsg1 + 317);
    const auto *qsg1_318 = buffer.data(qsg1 + 318);
    const auto *qsg1_320 = buffer.data(qsg1 + 320);
    const auto *qsg1_321 = buffer.data(qsg1 + 321);
    const auto *qsg1_325 = buffer.data(qsg1 + 325);
    const auto *qsg1_326 = buffer.data(qsg1 + 326);
    const auto *qsg1_327 = buffer.data(qsg1 + 327);
    const auto *qsg1_329 = buffer.data(qsg1 + 329);
    const auto *qsg1_335 = buffer.data(qsg1 + 335);
    const auto *qsg1_339 = buffer.data(qsg1 + 339);
    const auto *qsg1_342 = buffer.data(qsg1 + 342);
    const auto *qsg1_343 = buffer.data(qsg1 + 343);
    const auto *qsg1_344 = buffer.data(qsg1 + 344);
    const auto *qsg1_345 = buffer.data(qsg1 + 345);
    const auto *qsg1_348 = buffer.data(qsg1 + 348);
    const auto *qsg1_350 = buffer.data(qsg1 + 350);
    const auto *qsg1_351 = buffer.data(qsg1 + 351);
    const auto *qsg1_354 = buffer.data(qsg1 + 354);
    const auto *qsg1_355 = buffer.data(qsg1 + 355);
    const auto *qsg1_357 = buffer.data(qsg1 + 357);
    const auto *qsg1_358 = buffer.data(qsg1 + 358);
    const auto *qsg1_359 = buffer.data(qsg1 + 359);
    const auto *qsg1_360 = buffer.data(qsg1 + 360);
    const auto *qsg1_363 = buffer.data(qsg1 + 363);
    const auto *qsg1_365 = buffer.data(qsg1 + 365);
    const auto *qsg1_366 = buffer.data(qsg1 + 366);
    const auto *qsg1_369 = buffer.data(qsg1 + 369);
    const auto *qsg1_370 = buffer.data(qsg1 + 370);
    const auto *qsg1_372 = buffer.data(qsg1 + 372);
    const auto *qsg1_373 = buffer.data(qsg1 + 373);
    const auto *qsg1_374 = buffer.data(qsg1 + 374);
    const auto *qsg1_375 = buffer.data(qsg1 + 375);

    const auto *qsh_442 = buffer.data(qsh + 442);
    const auto *qsh_443 = buffer.data(qsh + 443);
    const auto *qsh_444 = buffer.data(qsh + 444);
    const auto *qsh_446 = buffer.data(qsh + 446);
    const auto *qsh_447 = buffer.data(qsh + 447);
    const auto *qsh_448 = buffer.data(qsh + 448);
    const auto *qsh_450 = buffer.data(qsh + 450);
    const auto *qsh_451 = buffer.data(qsh + 451);
    const auto *qsh_456 = buffer.data(qsh + 456);
    const auto *qsh_457 = buffer.data(qsh + 457);
    const auto *qsh_458 = buffer.data(qsh + 458);
    const auto *qsh_459 = buffer.data(qsh + 459);
    const auto *qsh_460 = buffer.data(qsh + 460);
    const auto *qsh_461 = buffer.data(qsh + 461);
    const auto *qsh_462 = buffer.data(qsh + 462);
    const auto *qsh_464 = buffer.data(qsh + 464);
    const auto *qsh_465 = buffer.data(qsh + 465);
    const auto *qsh_467 = buffer.data(qsh + 467);
    const auto *qsh_468 = buffer.data(qsh + 468);
    const auto *qsh_471 = buffer.data(qsh + 471);
    const auto *qsh_476 = buffer.data(qsh + 476);
    const auto *qsh_477 = buffer.data(qsh + 477);
    const auto *qsh_478 = buffer.data(qsh + 478);
    const auto *qsh_479 = buffer.data(qsh + 479);
    const auto *qsh_480 = buffer.data(qsh + 480);
    const auto *qsh_481 = buffer.data(qsh + 481);
    const auto *qsh_482 = buffer.data(qsh + 482);
    const auto *qsh_483 = buffer.data(qsh + 483);
    const auto *qsh_485 = buffer.data(qsh + 485);
    const auto *qsh_486 = buffer.data(qsh + 486);
    const auto *qsh_488 = buffer.data(qsh + 488);
    const auto *qsh_489 = buffer.data(qsh + 489);
    const auto *qsh_492 = buffer.data(qsh + 492);
    const auto *qsh_493 = buffer.data(qsh + 493);
    const auto *qsh_495 = buffer.data(qsh + 495);
    const auto *qsh_497 = buffer.data(qsh + 497);
    const auto *qsh_498 = buffer.data(qsh + 498);
    const auto *qsh_499 = buffer.data(qsh + 499);
    const auto *qsh_500 = buffer.data(qsh + 500);
    const auto *qsh_501 = buffer.data(qsh + 501);
    const auto *qsh_502 = buffer.data(qsh + 502);
    const auto *qsh_503 = buffer.data(qsh + 503);
    const auto *qsh_504 = buffer.data(qsh + 504);
    const auto *qsh_506 = buffer.data(qsh + 506);
    const auto *qsh_507 = buffer.data(qsh + 507);
    const auto *qsh_509 = buffer.data(qsh + 509);
    const auto *qsh_510 = buffer.data(qsh + 510);
    const auto *qsh_513 = buffer.data(qsh + 513);
    const auto *qsh_514 = buffer.data(qsh + 514);
    const auto *qsh_516 = buffer.data(qsh + 516);
    const auto *qsh_518 = buffer.data(qsh + 518);
    const auto *qsh_519 = buffer.data(qsh + 519);
    const auto *qsh_520 = buffer.data(qsh + 520);
    const auto *qsh_521 = buffer.data(qsh + 521);
    const auto *qsh_522 = buffer.data(qsh + 522);
    const auto *qsh_523 = buffer.data(qsh + 523);
    const auto *qsh_524 = buffer.data(qsh + 524);
    const auto *qsh_525 = buffer.data(qsh + 525);

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, osh_447, qsg0_315, qsg0_321, \
                         qsg1_315, qsg1_321, qsh_442, qsh_443, qsh_444, \
                         qsh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * qsh_442[k];

        t_593[k] = f_4 * qsg0_315[k]
                   - f_5 * qsg1_315[k]
                   + f_3 * pc_z[k] * qsh_443[k];

        t_594[k] = f_23 * osh_447[k]
                   + f_6 * qsg0_321[k]
                   - f_7 * qsg1_321[k]
                   + f_3 * pc_x[k] * qsh_447[k];

        t_595[k] = f_3 * pc_z[k] * qsh_444[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, osh_320, osh_451, \
                         qsg0_317, qsg0_325, qsg1_317, qsg1_325, qsh_446, qsh_447, \
                         qsh_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_23 * osh_320[k]
                   + f_3 * pc_y[k] * qsh_446[k];

        t_597[k] = f_6 * qsg0_317[k]
                   - f_7 * qsg1_317[k]
                   + f_3 * pc_z[k] * qsh_446[k];

        t_598[k] = f_23 * osh_451[k]
                   + f_4 * qsg0_325[k]
                   - f_5 * qsg1_325[k]
                   + f_3 * pc_x[k] * qsh_451[k];

        t_599[k] = f_3 * pc_z[k] * qsh_447[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, pc_z, osh_324, osh_456, \
                         qsg0_318, qsg0_320, qsg1_318, qsg1_320, qsh_448, qsh_450, \
                         qsh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_4 * qsg0_318[k]
                   - f_5 * qsg1_318[k]
                   + f_3 * pc_z[k] * qsh_448[k];

        t_601[k] = f_23 * osh_324[k]
                   + f_3 * pc_y[k] * qsh_450[k];

        t_602[k] = f_8 * qsg0_320[k]
                   - f_9 * qsg1_320[k]
                   + f_3 * pc_z[k] * qsh_450[k];

        t_603[k] = f_23 * osh_456[k]
                   + f_3 * pc_x[k] * qsh_456[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pc_x, pc_z, osh_458, osh_459, \
                         osh_460, osh_461, qsh_451, qsh_458, qsh_459, qsh_460, \
                         qsh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * qsh_451[k];

        t_605[k] = f_23 * osh_458[k]
                   + f_3 * pc_x[k] * qsh_458[k];

        t_606[k] = f_23 * osh_459[k]
                   + f_3 * pc_x[k] * qsh_459[k];

        t_607[k] = f_23 * osh_460[k]
                   + f_3 * pc_x[k] * qsh_460[k];

        t_608[k] = f_23 * osh_461[k]
                   + f_3 * pc_x[k] * qsh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pc_y, pc_z, osh_330, qsg0_325, qsg0_326, \
                         qsg1_325, qsg1_326, qsh_456, qsh_457, \
                         qsh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_23 * osh_330[k]
                   + f_1 * qsg0_325[k]
                   - f_2 * qsg1_325[k]
                   + f_3 * pc_y[k] * qsh_456[k];

        t_610[k] = f_3 * pc_z[k] * qsh_456[k];

        t_611[k] = f_4 * qsg0_325[k]
                   - f_5 * qsg1_325[k]
                   + f_3 * pc_z[k] * qsh_457[k];

        t_612[k] = f_6 * qsg0_326[k]
                   - f_7 * qsg1_326[k]
                   + f_3 * pc_z[k] * qsh_458[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pc_y, pc_z, osi0_420, osh_335, \
                         osi1_420, qsg0_327, qsg0_329, qsg1_327, qsg1_329, qsh_459, \
                         qsh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_8 * qsg0_327[k]
                   - f_9 * qsg1_327[k]
                   + f_3 * pc_z[k] * qsh_459[k];

        t_614[k] = f_23 * osh_335[k]
                   + f_3 * pc_y[k] * qsh_461[k];

        t_615[k] = f_1 * qsg0_329[k]
                   - f_2 * qsg1_329[k]
                   + f_3 * pc_z[k] * qsh_461[k];

        t_616[k] = pa_z[k] * osi0_420[k]
                   - f_10 * pc_z[k] * osi1_420[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pa_z, pc_y, pc_z, osi0_423, osh_315, \
                         osh_336, osh_338, osi1_423, qsh_462, qsh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_22 * osh_336[k]
                   + f_3 * pc_y[k] * qsh_462[k];

        t_618[k] = f_11 * osh_315[k]
                   + f_3 * pc_z[k] * qsh_462[k];

        t_619[k] = pa_z[k] * osi0_423[k]
                   - f_10 * pc_z[k] * osi1_423[k];

        t_620[k] = f_22 * osh_338[k]
                   + f_3 * pc_y[k] * qsh_464[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pa_z, pc_x, pc_z, osi0_426, osh_318, osh_467, \
                         osi1_426, qsg0_335, qsg1_335, qsh_465, \
                         qsh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_23 * osh_467[k]
                   + f_8 * qsg0_335[k]
                   - f_9 * qsg1_335[k]
                   + f_3 * pc_x[k] * qsh_467[k];

        t_622[k] = pa_z[k] * osi0_426[k]
                   - f_10 * pc_z[k] * osi1_426[k];

        t_623[k] = f_11 * osh_318[k]
                   + f_3 * pc_z[k] * qsh_465[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_z, pc_x, pc_y, pc_z, osi0_430, osh_341, \
                         osh_471, osi1_430, qsg0_339, qsg1_339, qsh_467, \
                         qsh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_22 * osh_341[k]
                   + f_3 * pc_y[k] * qsh_467[k];

        t_625[k] = f_23 * osh_471[k]
                   + f_6 * qsg0_339[k]
                   - f_7 * qsg1_339[k]
                   + f_3 * pc_x[k] * qsh_471[k];

        t_626[k] = pa_z[k] * osi0_430[k]
                   - f_10 * pc_z[k] * osi1_430[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_z, pc_y, pc_z, osi0_432, osh_321, osh_322, \
                         osh_345, osi1_432, qsh_468, qsh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_11 * osh_321[k]
                   + f_3 * pc_z[k] * qsh_468[k];

        t_628[k] = pa_z[k] * osi0_432[k]
                   + f_12 * osh_322[k]
                   - f_10 * pc_z[k] * osi1_432[k];

        t_629[k] = f_22 * osh_345[k]
                   + f_3 * pc_y[k] * qsh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, osh_476, osh_477, osh_478, osh_479, \
                         qsg0_344, qsg1_344, qsh_476, qsh_477, qsh_478, \
                         qsh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_23 * osh_476[k]
                   + f_4 * qsg0_344[k]
                   - f_5 * qsg1_344[k]
                   + f_3 * pc_x[k] * qsh_476[k];

        t_631[k] = f_23 * osh_477[k]
                   + f_3 * pc_x[k] * qsh_477[k];

        t_632[k] = f_23 * osh_478[k]
                   + f_3 * pc_x[k] * qsh_478[k];

        t_633[k] = f_23 * osh_479[k]
                   + f_3 * pc_x[k] * qsh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_z, pc_x, pc_z, osi0_441, osh_480, \
                         osh_481, osh_482, osi1_441, qsh_480, qsh_481, \
                         qsh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_23 * osh_480[k]
                   + f_3 * pc_x[k] * qsh_480[k];

        t_635[k] = f_23 * osh_481[k]
                   + f_3 * pc_x[k] * qsh_481[k];

        t_636[k] = f_23 * osh_482[k]
                   + f_3 * pc_x[k] * qsh_482[k];

        t_637[k] = pa_z[k] * osi0_441[k]
                   - f_10 * pc_z[k] * osi1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, osh_330, osh_353, osh_354, qsg0_342, \
                         qsg0_343, qsg1_342, qsg1_343, qsh_477, qsh_479, \
                         qsh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * osh_330[k]
                   + f_3 * pc_z[k] * qsh_477[k];

        t_639[k] = f_22 * osh_353[k]
                   + f_8 * qsg0_342[k]
                   - f_9 * qsg1_342[k]
                   + f_3 * pc_y[k] * qsh_479[k];

        t_640[k] = f_22 * osh_354[k]
                   + f_6 * qsg0_343[k]
                   - f_7 * qsg1_343[k]
                   + f_3 * pc_y[k] * qsh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, osh_335, osh_355, osh_356, qsg0_344, \
                         qsg1_344, qsh_481, qsh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_22 * osh_355[k]
                   + f_4 * qsg0_344[k]
                   - f_5 * qsg1_344[k]
                   + f_3 * pc_y[k] * qsh_481[k];

        t_642[k] = f_22 * osh_356[k]
                   + f_3 * pc_y[k] * qsh_482[k];

        t_643[k] = f_11 * osh_335[k]
                   + f_1 * qsg0_344[k]
                   - f_2 * qsg1_344[k]
                   + f_3 * pc_z[k] * qsh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, osh_336, osh_357, osh_483, \
                         qsg0_345, qsg1_345, qsh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_23 * osh_483[k]
                   + f_1 * qsg0_345[k]
                   - f_2 * qsg1_345[k]
                   + f_3 * pc_x[k] * qsh_483[k];

        t_645[k] = f_14 * osh_357[k]
                   + f_3 * pc_y[k] * qsh_483[k];

        t_646[k] = f_12 * osh_336[k]
                   + f_3 * pc_z[k] * qsh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, osh_359, osh_486, osh_488, qsg0_348, \
                         qsg0_350, qsg1_348, qsg1_350, qsh_485, qsh_486, \
                         qsh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_23 * osh_486[k]
                   + f_8 * qsg0_348[k]
                   - f_9 * qsg1_348[k]
                   + f_3 * pc_x[k] * qsh_486[k];

        t_648[k] = f_14 * osh_359[k]
                   + f_3 * pc_y[k] * qsh_485[k];

        t_649[k] = f_23 * osh_488[k]
                   + f_8 * qsg0_350[k]
                   - f_9 * qsg1_350[k]
                   + f_3 * pc_x[k] * qsh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, osh_339, osh_362, osh_489, \
                         qsg0_351, qsg1_351, qsh_486, qsh_488, \
                         qsh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_23 * osh_489[k]
                   + f_6 * qsg0_351[k]
                   - f_7 * qsg1_351[k]
                   + f_3 * pc_x[k] * qsh_489[k];

        t_651[k] = f_12 * osh_339[k]
                   + f_3 * pc_z[k] * qsh_486[k];

        t_652[k] = f_14 * osh_362[k]
                   + f_3 * pc_y[k] * qsh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, osh_342, osh_492, osh_493, qsg0_354, \
                         qsg0_355, qsg1_354, qsg1_355, qsh_489, qsh_492, \
                         qsh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_23 * osh_492[k]
                   + f_6 * qsg0_354[k]
                   - f_7 * qsg1_354[k]
                   + f_3 * pc_x[k] * qsh_492[k];

        t_654[k] = f_23 * osh_493[k]
                   + f_4 * qsg0_355[k]
                   - f_5 * qsg1_355[k]
                   + f_3 * pc_x[k] * qsh_493[k];

        t_655[k] = f_12 * osh_342[k]
                   + f_3 * pc_z[k] * qsh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, osh_366, osh_495, osh_497, qsg0_357, \
                         qsg0_359, qsg1_357, qsg1_359, qsh_492, qsh_495, \
                         qsh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_23 * osh_495[k]
                   + f_4 * qsg0_357[k]
                   - f_5 * qsg1_357[k]
                   + f_3 * pc_x[k] * qsh_495[k];

        t_657[k] = f_14 * osh_366[k]
                   + f_3 * pc_y[k] * qsh_492[k];

        t_658[k] = f_23 * osh_497[k]
                   + f_4 * qsg0_359[k]
                   - f_5 * qsg1_359[k]
                   + f_3 * pc_x[k] * qsh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, osh_498, osh_499, osh_500, \
                         osh_501, osh_502, qsh_498, qsh_499, qsh_500, qsh_501, \
                         qsh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_23 * osh_498[k]
                   + f_3 * pc_x[k] * qsh_498[k];

        t_660[k] = f_23 * osh_499[k]
                   + f_3 * pc_x[k] * qsh_499[k];

        t_661[k] = f_23 * osh_500[k]
                   + f_3 * pc_x[k] * qsh_500[k];

        t_662[k] = f_23 * osh_501[k]
                   + f_3 * pc_x[k] * qsh_501[k];

        t_663[k] = f_23 * osh_502[k]
                   + f_3 * pc_x[k] * qsh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, osh_351, osh_372, osh_503, \
                         qsg0_355, qsg1_355, qsh_498, qsh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_23 * osh_503[k]
                   + f_3 * pc_x[k] * qsh_503[k];

        t_665[k] = f_14 * osh_372[k]
                   + f_1 * qsg0_355[k]
                   - f_2 * qsg1_355[k]
                   + f_3 * pc_y[k] * qsh_498[k];

        t_666[k] = f_12 * osh_351[k]
                   + f_3 * pc_z[k] * qsh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, osh_374, osh_375, osh_376, qsg0_357, \
                         qsg0_358, qsg0_359, qsg1_357, qsg1_358, qsg1_359, qsh_500, qsh_501, \
                         qsh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * osh_374[k]
                   + f_8 * qsg0_357[k]
                   - f_9 * qsg1_357[k]
                   + f_3 * pc_y[k] * qsh_500[k];

        t_668[k] = f_14 * osh_375[k]
                   + f_6 * qsg0_358[k]
                   - f_7 * qsg1_358[k]
                   + f_3 * pc_y[k] * qsh_501[k];

        t_669[k] = f_14 * osh_376[k]
                   + f_4 * qsg0_359[k]
                   - f_5 * qsg1_359[k]
                   + f_3 * pc_y[k] * qsh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, osh_356, osh_377, osh_504, \
                         qsg0_359, qsg0_360, qsg1_359, qsg1_360, qsh_503, \
                         qsh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * osh_377[k]
                   + f_3 * pc_y[k] * qsh_503[k];

        t_671[k] = f_12 * osh_356[k]
                   + f_1 * qsg0_359[k]
                   - f_2 * qsg1_359[k]
                   + f_3 * pc_z[k] * qsh_503[k];

        t_672[k] = f_23 * osh_504[k]
                   + f_1 * qsg0_360[k]
                   - f_2 * qsg1_360[k]
                   + f_3 * pc_x[k] * qsh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, osh_357, osh_378, \
                         osh_380, osh_507, qsg0_363, qsg1_363, qsh_504, qsh_506, \
                         qsh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * osh_378[k]
                   + f_3 * pc_y[k] * qsh_504[k];

        t_674[k] = f_13 * osh_357[k]
                   + f_3 * pc_z[k] * qsh_504[k];

        t_675[k] = f_23 * osh_507[k]
                   + f_8 * qsg0_363[k]
                   - f_9 * qsg1_363[k]
                   + f_3 * pc_x[k] * qsh_507[k];

        t_676[k] = f_13 * osh_380[k]
                   + f_3 * pc_y[k] * qsh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, osh_360, osh_509, osh_510, qsg0_365, \
                         qsg0_366, qsg1_365, qsg1_366, qsh_507, qsh_509, \
                         qsh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_23 * osh_509[k]
                   + f_8 * qsg0_365[k]
                   - f_9 * qsg1_365[k]
                   + f_3 * pc_x[k] * qsh_509[k];

        t_678[k] = f_23 * osh_510[k]
                   + f_6 * qsg0_366[k]
                   - f_7 * qsg1_366[k]
                   + f_3 * pc_x[k] * qsh_510[k];

        t_679[k] = f_13 * osh_360[k]
                   + f_3 * pc_z[k] * qsh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, osh_383, osh_513, osh_514, qsg0_369, \
                         qsg0_370, qsg1_369, qsg1_370, qsh_509, qsh_513, \
                         qsh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * osh_383[k]
                   + f_3 * pc_y[k] * qsh_509[k];

        t_681[k] = f_23 * osh_513[k]
                   + f_6 * qsg0_369[k]
                   - f_7 * qsg1_369[k]
                   + f_3 * pc_x[k] * qsh_513[k];

        t_682[k] = f_23 * osh_514[k]
                   + f_4 * qsg0_370[k]
                   - f_5 * qsg1_370[k]
                   + f_3 * pc_x[k] * qsh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, osh_363, osh_387, osh_516, \
                         qsg0_372, qsg1_372, qsh_510, qsh_513, \
                         qsh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * osh_363[k]
                   + f_3 * pc_z[k] * qsh_510[k];

        t_684[k] = f_23 * osh_516[k]
                   + f_4 * qsg0_372[k]
                   - f_5 * qsg1_372[k]
                   + f_3 * pc_x[k] * qsh_516[k];

        t_685[k] = f_13 * osh_387[k]
                   + f_3 * pc_y[k] * qsh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, osh_518, osh_519, osh_520, osh_521, \
                         qsg0_374, qsg1_374, qsh_518, qsh_519, qsh_520, \
                         qsh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_23 * osh_518[k]
                   + f_4 * qsg0_374[k]
                   - f_5 * qsg1_374[k]
                   + f_3 * pc_x[k] * qsh_518[k];

        t_687[k] = f_23 * osh_519[k]
                   + f_3 * pc_x[k] * qsh_519[k];

        t_688[k] = f_23 * osh_520[k]
                   + f_3 * pc_x[k] * qsh_520[k];

        t_689[k] = f_23 * osh_521[k]
                   + f_3 * pc_x[k] * qsh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, osh_393, osh_522, osh_523, \
                         osh_524, qsg0_370, qsg1_370, qsh_519, qsh_522, qsh_523, \
                         qsh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_23 * osh_522[k]
                   + f_3 * pc_x[k] * qsh_522[k];

        t_691[k] = f_23 * osh_523[k]
                   + f_3 * pc_x[k] * qsh_523[k];

        t_692[k] = f_23 * osh_524[k]
                   + f_3 * pc_x[k] * qsh_524[k];

        t_693[k] = f_13 * osh_393[k]
                   + f_1 * qsg0_370[k]
                   - f_2 * qsg1_370[k]
                   + f_3 * pc_y[k] * qsh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, osh_372, osh_395, osh_396, qsg0_372, \
                         qsg0_373, qsg1_372, qsg1_373, qsh_519, qsh_521, \
                         qsh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * osh_372[k]
                   + f_3 * pc_z[k] * qsh_519[k];

        t_695[k] = f_13 * osh_395[k]
                   + f_8 * qsg0_372[k]
                   - f_9 * qsg1_372[k]
                   + f_3 * pc_y[k] * qsh_521[k];

        t_696[k] = f_13 * osh_396[k]
                   + f_6 * qsg0_373[k]
                   - f_7 * qsg1_373[k]
                   + f_3 * pc_y[k] * qsh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, osh_377, osh_397, osh_398, qsg0_374, \
                         qsg1_374, qsh_523, qsh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * osh_397[k]
                   + f_4 * qsg0_374[k]
                   - f_5 * qsg1_374[k]
                   + f_3 * pc_y[k] * qsh_523[k];

        t_698[k] = f_13 * osh_398[k]
                   + f_3 * pc_y[k] * qsh_524[k];

        t_699[k] = f_13 * osh_377[k]
                   + f_1 * qsg0_374[k]
                   - f_2 * qsg1_374[k]
                   + f_3 * pc_z[k] * qsh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, osh_378, osh_399, osh_525, \
                         qsg0_375, qsg1_375, qsh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_23 * osh_525[k]
                   + f_1 * qsg0_375[k]
                   - f_2 * qsg1_375[k]
                   + f_3 * pc_x[k] * qsh_525[k];

        t_701[k] = f_12 * osh_399[k]
                   + f_3 * pc_y[k] * qsh_525[k];

        t_702[k] = f_14 * osh_378[k]
                   + f_3 * pc_z[k] * qsh_525[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_560 = buffer.data(osi0 + 560);
    const auto *osi0_563 = buffer.data(osi0 + 563);
    const auto *osi0_565 = buffer.data(osi0 + 565);
    const auto *osi0_566 = buffer.data(osi0 + 566);
    const auto *osi0_569 = buffer.data(osi0 + 569);
    const auto *osi0_570 = buffer.data(osi0 + 570);
    const auto *osi0_572 = buffer.data(osi0 + 572);
    const auto *osi0_574 = buffer.data(osi0 + 574);
    const auto *osi0_587 = buffer.data(osi0 + 587);
    const auto *osi0_588 = buffer.data(osi0 + 588);
    const auto *osi0_591 = buffer.data(osi0 + 591);
    const auto *osi0_594 = buffer.data(osi0 + 594);

    const auto *osh_381 = buffer.data(osh + 381);
    const auto *osh_384 = buffer.data(osh + 384);
    const auto *osh_393 = buffer.data(osh + 393);
    const auto *osh_398 = buffer.data(osh + 398);
    const auto *osh_399 = buffer.data(osh + 399);
    const auto *osh_401 = buffer.data(osh + 401);
    const auto *osh_402 = buffer.data(osh + 402);
    const auto *osh_404 = buffer.data(osh + 404);
    const auto *osh_405 = buffer.data(osh + 405);
    const auto *osh_408 = buffer.data(osh + 408);
    const auto *osh_414 = buffer.data(osh + 414);
    const auto *osh_416 = buffer.data(osh + 416);
    const auto *osh_417 = buffer.data(osh + 417);
    const auto *osh_418 = buffer.data(osh + 418);
    const auto *osh_419 = buffer.data(osh + 419);
    const auto *osh_420 = buffer.data(osh + 420);
    const auto *osh_421 = buffer.data(osh + 421);
    const auto *osh_422 = buffer.data(osh + 422);
    const auto *osh_423 = buffer.data(osh + 423);
    const auto *osh_425 = buffer.data(osh + 425);
    const auto *osh_426 = buffer.data(osh + 426);
    const auto *osh_428 = buffer.data(osh + 428);
    const auto *osh_429 = buffer.data(osh + 429);
    const auto *osh_435 = buffer.data(osh + 435);
    const auto *osh_437 = buffer.data(osh + 437);
    const auto *osh_438 = buffer.data(osh + 438);
    const auto *osh_439 = buffer.data(osh + 439);
    const auto *osh_440 = buffer.data(osh + 440);
    const auto *osh_441 = buffer.data(osh + 441);
    const auto *osh_444 = buffer.data(osh + 444);
    const auto *osh_446 = buffer.data(osh + 446);
    const auto *osh_450 = buffer.data(osh + 450);
    const auto *osh_456 = buffer.data(osh + 456);
    const auto *osh_461 = buffer.data(osh + 461);
    const auto *osh_462 = buffer.data(osh + 462);
    const auto *osh_464 = buffer.data(osh + 464);
    const auto *osh_528 = buffer.data(osh + 528);
    const auto *osh_530 = buffer.data(osh + 530);
    const auto *osh_531 = buffer.data(osh + 531);
    const auto *osh_534 = buffer.data(osh + 534);
    const auto *osh_535 = buffer.data(osh + 535);
    const auto *osh_537 = buffer.data(osh + 537);
    const auto *osh_539 = buffer.data(osh + 539);
    const auto *osh_540 = buffer.data(osh + 540);
    const auto *osh_541 = buffer.data(osh + 541);
    const auto *osh_542 = buffer.data(osh + 542);
    const auto *osh_543 = buffer.data(osh + 543);
    const auto *osh_544 = buffer.data(osh + 544);
    const auto *osh_545 = buffer.data(osh + 545);
    const auto *osh_561 = buffer.data(osh + 561);
    const auto *osh_562 = buffer.data(osh + 562);
    const auto *osh_563 = buffer.data(osh + 563);
    const auto *osh_564 = buffer.data(osh + 564);
    const auto *osh_565 = buffer.data(osh + 565);
    const auto *osh_566 = buffer.data(osh + 566);
    const auto *osh_567 = buffer.data(osh + 567);
    const auto *osh_572 = buffer.data(osh + 572);
    const auto *osh_576 = buffer.data(osh + 576);
    const auto *osh_581 = buffer.data(osh + 581);
    const auto *osh_582 = buffer.data(osh + 582);
    const auto *osh_583 = buffer.data(osh + 583);
    const auto *osh_584 = buffer.data(osh + 584);
    const auto *osh_585 = buffer.data(osh + 585);
    const auto *osh_587 = buffer.data(osh + 587);
    const auto *osh_588 = buffer.data(osh + 588);
    const auto *osh_591 = buffer.data(osh + 591);
    const auto *osh_594 = buffer.data(osh + 594);
    const auto *osh_598 = buffer.data(osh + 598);
    const auto *osh_603 = buffer.data(osh + 603);
    const auto *osh_605 = buffer.data(osh + 605);
    const auto *osh_606 = buffer.data(osh + 606);
    const auto *osh_607 = buffer.data(osh + 607);
    const auto *osh_608 = buffer.data(osh + 608);
    const auto *osh_614 = buffer.data(osh + 614);

    const auto *osi1_560 = buffer.data(osi1 + 560);
    const auto *osi1_563 = buffer.data(osi1 + 563);
    const auto *osi1_565 = buffer.data(osi1 + 565);
    const auto *osi1_566 = buffer.data(osi1 + 566);
    const auto *osi1_569 = buffer.data(osi1 + 569);
    const auto *osi1_570 = buffer.data(osi1 + 570);
    const auto *osi1_572 = buffer.data(osi1 + 572);
    const auto *osi1_574 = buffer.data(osi1 + 574);
    const auto *osi1_587 = buffer.data(osi1 + 587);
    const auto *osi1_588 = buffer.data(osi1 + 588);
    const auto *osi1_591 = buffer.data(osi1 + 591);
    const auto *osi1_594 = buffer.data(osi1 + 594);

    const auto *qsg0_378 = buffer.data(qsg0 + 378);
    const auto *qsg0_380 = buffer.data(qsg0 + 380);
    const auto *qsg0_381 = buffer.data(qsg0 + 381);
    const auto *qsg0_384 = buffer.data(qsg0 + 384);
    const auto *qsg0_385 = buffer.data(qsg0 + 385);
    const auto *qsg0_387 = buffer.data(qsg0 + 387);
    const auto *qsg0_388 = buffer.data(qsg0 + 388);
    const auto *qsg0_389 = buffer.data(qsg0 + 389);
    const auto *qsg0_400 = buffer.data(qsg0 + 400);
    const auto *qsg0_402 = buffer.data(qsg0 + 402);
    const auto *qsg0_403 = buffer.data(qsg0 + 403);
    const auto *qsg0_404 = buffer.data(qsg0 + 404);
    const auto *qsg0_405 = buffer.data(qsg0 + 405);
    const auto *qsg0_406 = buffer.data(qsg0 + 406);
    const auto *qsg0_407 = buffer.data(qsg0 + 407);
    const auto *qsg0_408 = buffer.data(qsg0 + 408);
    const auto *qsg0_409 = buffer.data(qsg0 + 409);
    const auto *qsg0_410 = buffer.data(qsg0 + 410);
    const auto *qsg0_414 = buffer.data(qsg0 + 414);
    const auto *qsg0_415 = buffer.data(qsg0 + 415);
    const auto *qsg0_416 = buffer.data(qsg0 + 416);
    const auto *qsg0_417 = buffer.data(qsg0 + 417);
    const auto *qsg0_418 = buffer.data(qsg0 + 418);
    const auto *qsg0_419 = buffer.data(qsg0 + 419);
    const auto *qsg0_420 = buffer.data(qsg0 + 420);
    const auto *qsg0_422 = buffer.data(qsg0 + 422);
    const auto *qsg0_423 = buffer.data(qsg0 + 423);
    const auto *qsg0_425 = buffer.data(qsg0 + 425);
    const auto *qsg0_426 = buffer.data(qsg0 + 426);
    const auto *qsg0_430 = buffer.data(qsg0 + 430);
    const auto *qsg0_431 = buffer.data(qsg0 + 431);
    const auto *qsg0_432 = buffer.data(qsg0 + 432);
    const auto *qsg0_434 = buffer.data(qsg0 + 434);
    const auto *qsg0_440 = buffer.data(qsg0 + 440);

    const auto *qsg1_378 = buffer.data(qsg1 + 378);
    const auto *qsg1_380 = buffer.data(qsg1 + 380);
    const auto *qsg1_381 = buffer.data(qsg1 + 381);
    const auto *qsg1_384 = buffer.data(qsg1 + 384);
    const auto *qsg1_385 = buffer.data(qsg1 + 385);
    const auto *qsg1_387 = buffer.data(qsg1 + 387);
    const auto *qsg1_388 = buffer.data(qsg1 + 388);
    const auto *qsg1_389 = buffer.data(qsg1 + 389);
    const auto *qsg1_400 = buffer.data(qsg1 + 400);
    const auto *qsg1_402 = buffer.data(qsg1 + 402);
    const auto *qsg1_403 = buffer.data(qsg1 + 403);
    const auto *qsg1_404 = buffer.data(qsg1 + 404);
    const auto *qsg1_405 = buffer.data(qsg1 + 405);
    const auto *qsg1_406 = buffer.data(qsg1 + 406);
    const auto *qsg1_407 = buffer.data(qsg1 + 407);
    const auto *qsg1_408 = buffer.data(qsg1 + 408);
    const auto *qsg1_409 = buffer.data(qsg1 + 409);
    const auto *qsg1_410 = buffer.data(qsg1 + 410);
    const auto *qsg1_414 = buffer.data(qsg1 + 414);
    const auto *qsg1_415 = buffer.data(qsg1 + 415);
    const auto *qsg1_416 = buffer.data(qsg1 + 416);
    const auto *qsg1_417 = buffer.data(qsg1 + 417);
    const auto *qsg1_418 = buffer.data(qsg1 + 418);
    const auto *qsg1_419 = buffer.data(qsg1 + 419);
    const auto *qsg1_420 = buffer.data(qsg1 + 420);
    const auto *qsg1_422 = buffer.data(qsg1 + 422);
    const auto *qsg1_423 = buffer.data(qsg1 + 423);
    const auto *qsg1_425 = buffer.data(qsg1 + 425);
    const auto *qsg1_426 = buffer.data(qsg1 + 426);
    const auto *qsg1_430 = buffer.data(qsg1 + 430);
    const auto *qsg1_431 = buffer.data(qsg1 + 431);
    const auto *qsg1_432 = buffer.data(qsg1 + 432);
    const auto *qsg1_434 = buffer.data(qsg1 + 434);
    const auto *qsg1_440 = buffer.data(qsg1 + 440);

    const auto *qsh_527 = buffer.data(qsh + 527);
    const auto *qsh_528 = buffer.data(qsh + 528);
    const auto *qsh_530 = buffer.data(qsh + 530);
    const auto *qsh_531 = buffer.data(qsh + 531);
    const auto *qsh_534 = buffer.data(qsh + 534);
    const auto *qsh_535 = buffer.data(qsh + 535);
    const auto *qsh_537 = buffer.data(qsh + 537);
    const auto *qsh_539 = buffer.data(qsh + 539);
    const auto *qsh_540 = buffer.data(qsh + 540);
    const auto *qsh_541 = buffer.data(qsh + 541);
    const auto *qsh_542 = buffer.data(qsh + 542);
    const auto *qsh_543 = buffer.data(qsh + 543);
    const auto *qsh_544 = buffer.data(qsh + 544);
    const auto *qsh_545 = buffer.data(qsh + 545);
    const auto *qsh_546 = buffer.data(qsh + 546);
    const auto *qsh_548 = buffer.data(qsh + 548);
    const auto *qsh_549 = buffer.data(qsh + 549);
    const auto *qsh_551 = buffer.data(qsh + 551);
    const auto *qsh_552 = buffer.data(qsh + 552);
    const auto *qsh_555 = buffer.data(qsh + 555);
    const auto *qsh_561 = buffer.data(qsh + 561);
    const auto *qsh_562 = buffer.data(qsh + 562);
    const auto *qsh_563 = buffer.data(qsh + 563);
    const auto *qsh_564 = buffer.data(qsh + 564);
    const auto *qsh_565 = buffer.data(qsh + 565);
    const auto *qsh_566 = buffer.data(qsh + 566);
    const auto *qsh_567 = buffer.data(qsh + 567);
    const auto *qsh_568 = buffer.data(qsh + 568);
    const auto *qsh_569 = buffer.data(qsh + 569);
    const auto *qsh_570 = buffer.data(qsh + 570);
    const auto *qsh_571 = buffer.data(qsh + 571);
    const auto *qsh_572 = buffer.data(qsh + 572);
    const auto *qsh_573 = buffer.data(qsh + 573);
    const auto *qsh_574 = buffer.data(qsh + 574);
    const auto *qsh_575 = buffer.data(qsh + 575);
    const auto *qsh_576 = buffer.data(qsh + 576);
    const auto *qsh_581 = buffer.data(qsh + 581);
    const auto *qsh_582 = buffer.data(qsh + 582);
    const auto *qsh_583 = buffer.data(qsh + 583);
    const auto *qsh_584 = buffer.data(qsh + 584);
    const auto *qsh_585 = buffer.data(qsh + 585);
    const auto *qsh_586 = buffer.data(qsh + 586);
    const auto *qsh_587 = buffer.data(qsh + 587);
    const auto *qsh_588 = buffer.data(qsh + 588);
    const auto *qsh_589 = buffer.data(qsh + 589);
    const auto *qsh_590 = buffer.data(qsh + 590);
    const auto *qsh_591 = buffer.data(qsh + 591);
    const auto *qsh_593 = buffer.data(qsh + 593);
    const auto *qsh_594 = buffer.data(qsh + 594);
    const auto *qsh_595 = buffer.data(qsh + 595);
    const auto *qsh_597 = buffer.data(qsh + 597);
    const auto *qsh_598 = buffer.data(qsh + 598);
    const auto *qsh_603 = buffer.data(qsh + 603);
    const auto *qsh_604 = buffer.data(qsh + 604);
    const auto *qsh_605 = buffer.data(qsh + 605);
    const auto *qsh_606 = buffer.data(qsh + 606);
    const auto *qsh_607 = buffer.data(qsh + 607);
    const auto *qsh_608 = buffer.data(qsh + 608);
    const auto *qsh_609 = buffer.data(qsh + 609);
    const auto *qsh_611 = buffer.data(qsh + 611);
    const auto *qsh_612 = buffer.data(qsh + 612);
    const auto *qsh_614 = buffer.data(qsh + 614);

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, osh_401, osh_528, osh_530, qsg0_378, \
                         qsg0_380, qsg1_378, qsg1_380, qsh_527, qsh_528, \
                         qsh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_23 * osh_528[k]
                   + f_8 * qsg0_378[k]
                   - f_9 * qsg1_378[k]
                   + f_3 * pc_x[k] * qsh_528[k];

        t_704[k] = f_12 * osh_401[k]
                   + f_3 * pc_y[k] * qsh_527[k];

        t_705[k] = f_23 * osh_530[k]
                   + f_8 * qsg0_380[k]
                   - f_9 * qsg1_380[k]
                   + f_3 * pc_x[k] * qsh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, osh_381, osh_404, osh_531, \
                         qsg0_381, qsg1_381, qsh_528, qsh_530, \
                         qsh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_23 * osh_531[k]
                   + f_6 * qsg0_381[k]
                   - f_7 * qsg1_381[k]
                   + f_3 * pc_x[k] * qsh_531[k];

        t_707[k] = f_14 * osh_381[k]
                   + f_3 * pc_z[k] * qsh_528[k];

        t_708[k] = f_12 * osh_404[k]
                   + f_3 * pc_y[k] * qsh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, osh_384, osh_534, osh_535, qsg0_384, \
                         qsg0_385, qsg1_384, qsg1_385, qsh_531, qsh_534, \
                         qsh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_23 * osh_534[k]
                   + f_6 * qsg0_384[k]
                   - f_7 * qsg1_384[k]
                   + f_3 * pc_x[k] * qsh_534[k];

        t_710[k] = f_23 * osh_535[k]
                   + f_4 * qsg0_385[k]
                   - f_5 * qsg1_385[k]
                   + f_3 * pc_x[k] * qsh_535[k];

        t_711[k] = f_14 * osh_384[k]
                   + f_3 * pc_z[k] * qsh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, osh_408, osh_537, osh_539, qsg0_387, \
                         qsg0_389, qsg1_387, qsg1_389, qsh_534, qsh_537, \
                         qsh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_23 * osh_537[k]
                   + f_4 * qsg0_387[k]
                   - f_5 * qsg1_387[k]
                   + f_3 * pc_x[k] * qsh_537[k];

        t_713[k] = f_12 * osh_408[k]
                   + f_3 * pc_y[k] * qsh_534[k];

        t_714[k] = f_23 * osh_539[k]
                   + f_4 * qsg0_389[k]
                   - f_5 * qsg1_389[k]
                   + f_3 * pc_x[k] * qsh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, osh_540, osh_541, osh_542, \
                         osh_543, osh_544, qsh_540, qsh_541, qsh_542, qsh_543, \
                         qsh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_23 * osh_540[k]
                   + f_3 * pc_x[k] * qsh_540[k];

        t_716[k] = f_23 * osh_541[k]
                   + f_3 * pc_x[k] * qsh_541[k];

        t_717[k] = f_23 * osh_542[k]
                   + f_3 * pc_x[k] * qsh_542[k];

        t_718[k] = f_23 * osh_543[k]
                   + f_3 * pc_x[k] * qsh_543[k];

        t_719[k] = f_23 * osh_544[k]
                   + f_3 * pc_x[k] * qsh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, osh_393, osh_414, osh_545, \
                         qsg0_385, qsg1_385, qsh_540, qsh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_23 * osh_545[k]
                   + f_3 * pc_x[k] * qsh_545[k];

        t_721[k] = f_12 * osh_414[k]
                   + f_1 * qsg0_385[k]
                   - f_2 * qsg1_385[k]
                   + f_3 * pc_y[k] * qsh_540[k];

        t_722[k] = f_14 * osh_393[k]
                   + f_3 * pc_z[k] * qsh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, osh_416, osh_417, osh_418, qsg0_387, \
                         qsg0_388, qsg0_389, qsg1_387, qsg1_388, qsg1_389, qsh_542, qsh_543, \
                         qsh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * osh_416[k]
                   + f_8 * qsg0_387[k]
                   - f_9 * qsg1_387[k]
                   + f_3 * pc_y[k] * qsh_542[k];

        t_724[k] = f_12 * osh_417[k]
                   + f_6 * qsg0_388[k]
                   - f_7 * qsg1_388[k]
                   + f_3 * pc_y[k] * qsh_543[k];

        t_725[k] = f_12 * osh_418[k]
                   + f_4 * qsg0_389[k]
                   - f_5 * qsg1_389[k]
                   + f_3 * pc_y[k] * qsh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_y, pc_y, pc_z, osi0_560, osh_398, \
                         osh_419, osh_420, osi1_560, qsg0_389, qsg1_389, qsh_545, \
                         qsh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * osh_419[k]
                   + f_3 * pc_y[k] * qsh_545[k];

        t_727[k] = f_14 * osh_398[k]
                   + f_1 * qsg0_389[k]
                   - f_2 * qsg1_389[k]
                   + f_3 * pc_z[k] * qsh_545[k];

        t_728[k] = pa_y[k] * osi0_560[k]
                   - f_10 * pc_y[k] * osi1_560[k];

        t_729[k] = f_11 * osh_420[k]
                   + f_3 * pc_y[k] * qsh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_y, pc_y, pc_z, osi0_563, osi0_565, \
                         osh_399, osh_421, osh_422, osi1_563, osi1_565, qsh_546, \
                         qsh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_22 * osh_399[k]
                   + f_3 * pc_z[k] * qsh_546[k];

        t_731[k] = pa_y[k] * osi0_563[k]
                   + f_12 * osh_421[k]
                   - f_10 * pc_y[k] * osi1_563[k];

        t_732[k] = f_11 * osh_422[k]
                   + f_3 * pc_y[k] * qsh_548[k];

        t_733[k] = pa_y[k] * osi0_565[k]
                   - f_10 * pc_y[k] * osi1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pa_y, pc_y, pc_z, osi0_566, osi0_569, \
                         osh_402, osh_423, osh_425, osi1_566, osi1_569, qsh_549, \
                         qsh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pa_y[k] * osi0_566[k]
                   + f_13 * osh_423[k]
                   - f_10 * pc_y[k] * osi1_566[k];

        t_735[k] = f_22 * osh_402[k]
                   + f_3 * pc_z[k] * qsh_549[k];

        t_736[k] = f_11 * osh_425[k]
                   + f_3 * pc_y[k] * qsh_551[k];

        t_737[k] = pa_y[k] * osi0_569[k]
                   - f_10 * pc_y[k] * osi1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pa_y, pc_y, pc_z, osi0_570, osi0_572, osh_405, \
                         osh_426, osh_428, osi1_570, osi1_572, \
                         qsh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pa_y[k] * osi0_570[k]
                   + f_14 * osh_426[k]
                   - f_10 * pc_y[k] * osi1_570[k];

        t_739[k] = f_22 * osh_405[k]
                   + f_3 * pc_z[k] * qsh_552[k];

        t_740[k] = pa_y[k] * osi0_572[k]
                   + f_12 * osh_428[k]
                   - f_10 * pc_y[k] * osi1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pa_y, pc_x, pc_y, osi0_574, osh_429, \
                         osh_561, osh_562, osi1_574, qsh_555, qsh_561, \
                         qsh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * osh_429[k]
                   + f_3 * pc_y[k] * qsh_555[k];

        t_742[k] = pa_y[k] * osi0_574[k]
                   - f_10 * pc_y[k] * osi1_574[k];

        t_743[k] = f_23 * osh_561[k]
                   + f_3 * pc_x[k] * qsh_561[k];

        t_744[k] = f_23 * osh_562[k]
                   + f_3 * pc_x[k] * qsh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, osh_563, osh_564, osh_565, osh_566, \
                         qsh_563, qsh_564, qsh_565, qsh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_23 * osh_563[k]
                   + f_3 * pc_x[k] * qsh_563[k];

        t_746[k] = f_23 * osh_564[k]
                   + f_3 * pc_x[k] * qsh_564[k];

        t_747[k] = f_23 * osh_565[k]
                   + f_3 * pc_x[k] * qsh_565[k];

        t_748[k] = f_23 * osh_566[k]
                   + f_3 * pc_x[k] * qsh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, osh_414, osh_435, osh_437, qsg0_400, \
                         qsg0_402, qsg1_400, qsg1_402, qsh_561, \
                         qsh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * osh_435[k]
                   + f_1 * qsg0_400[k]
                   - f_2 * qsg1_400[k]
                   + f_3 * pc_y[k] * qsh_561[k];

        t_750[k] = f_22 * osh_414[k]
                   + f_3 * pc_z[k] * qsh_561[k];

        t_751[k] = f_11 * osh_437[k]
                   + f_8 * qsg0_402[k]
                   - f_9 * qsg1_402[k]
                   + f_3 * pc_y[k] * qsh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, osh_438, osh_439, osh_440, qsg0_403, \
                         qsg0_404, qsg1_403, qsg1_404, qsh_564, qsh_565, \
                         qsh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * osh_438[k]
                   + f_6 * qsg0_403[k]
                   - f_7 * qsg1_403[k]
                   + f_3 * pc_y[k] * qsh_564[k];

        t_753[k] = f_11 * osh_439[k]
                   + f_4 * qsg0_404[k]
                   - f_5 * qsg1_404[k]
                   + f_3 * pc_y[k] * qsh_565[k];

        t_754[k] = f_11 * osh_440[k]
                   + f_3 * pc_y[k] * qsh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pa_y, pc_x, pc_y, pc_z, osi0_587, \
                         osh_420, osh_567, osi1_587, qsg0_405, qsg1_405, \
                         qsh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pa_y[k] * osi0_587[k]
                   - f_10 * pc_y[k] * osi1_587[k];

        t_756[k] = f_23 * osh_567[k]
                   + f_1 * qsg0_405[k]
                   - f_2 * qsg1_405[k]
                   + f_3 * pc_x[k] * qsh_567[k];

        t_757[k] = f_3 * pc_y[k] * qsh_567[k];

        t_758[k] = f_23 * osh_420[k]
                   + f_3 * pc_z[k] * qsh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, osh_572, qsg0_405, qsg0_410, \
                         qsg1_405, qsg1_410, qsh_568, qsh_569, \
                         qsh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_4 * qsg0_405[k]
                   - f_5 * qsg1_405[k]
                   + f_3 * pc_y[k] * qsh_568[k];

        t_760[k] = f_3 * pc_y[k] * qsh_569[k];

        t_761[k] = f_23 * osh_572[k]
                   + f_8 * qsg0_410[k]
                   - f_9 * qsg1_410[k]
                   + f_3 * pc_x[k] * qsh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_y, qsg0_406, qsg0_407, qsg1_406, qsg1_407, \
                         qsh_570, qsh_571, qsh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_6 * qsg0_406[k]
                   - f_7 * qsg1_406[k]
                   + f_3 * pc_y[k] * qsh_570[k];

        t_763[k] = f_4 * qsg0_407[k]
                   - f_5 * qsg1_407[k]
                   + f_3 * pc_y[k] * qsh_571[k];

        t_764[k] = f_3 * pc_y[k] * qsh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_y, osh_576, qsg0_408, qsg0_409, \
                         qsg0_414, qsg1_408, qsg1_409, qsg1_414, qsh_573, qsh_574, \
                         qsh_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_23 * osh_576[k]
                   + f_6 * qsg0_414[k]
                   - f_7 * qsg1_414[k]
                   + f_3 * pc_x[k] * qsh_576[k];

        t_766[k] = f_8 * qsg0_408[k]
                   - f_9 * qsg1_408[k]
                   + f_3 * pc_y[k] * qsh_573[k];

        t_767[k] = f_6 * qsg0_409[k]
                   - f_7 * qsg1_409[k]
                   + f_3 * pc_y[k] * qsh_574[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, osh_581, osh_582, qsg0_410, \
                         qsg0_419, qsg1_410, qsg1_419, qsh_575, qsh_576, qsh_581, \
                         qsh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * qsg0_410[k]
                   - f_5 * qsg1_410[k]
                   + f_3 * pc_y[k] * qsh_575[k];

        t_769[k] = f_3 * pc_y[k] * qsh_576[k];

        t_770[k] = f_23 * osh_581[k]
                   + f_4 * qsg0_419[k]
                   - f_5 * qsg1_419[k]
                   + f_3 * pc_x[k] * qsh_581[k];

        t_771[k] = f_23 * osh_582[k]
                   + f_3 * pc_x[k] * qsh_582[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pc_x, pc_y, osh_583, osh_584, \
                         osh_585, osh_587, qsh_581, qsh_583, qsh_584, qsh_585, \
                         qsh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_23 * osh_583[k]
                   + f_3 * pc_x[k] * qsh_583[k];

        t_773[k] = f_23 * osh_584[k]
                   + f_3 * pc_x[k] * qsh_584[k];

        t_774[k] = f_23 * osh_585[k]
                   + f_3 * pc_x[k] * qsh_585[k];

        t_775[k] = f_3 * pc_y[k] * qsh_581[k];

        t_776[k] = f_23 * osh_587[k]
                   + f_3 * pc_x[k] * qsh_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, qsg0_415, qsg0_416, qsg0_417, qsg1_415, \
                         qsg1_416, qsg1_417, qsh_582, qsh_583, \
                         qsh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_1 * qsg0_415[k]
                   - f_2 * qsg1_415[k]
                   + f_3 * pc_y[k] * qsh_582[k];

        t_778[k] = f_16 * qsg0_416[k]
                   - f_17 * qsg1_416[k]
                   + f_3 * pc_y[k] * qsh_583[k];

        t_779[k] = f_8 * qsg0_417[k]
                   - f_9 * qsg1_417[k]
                   + f_3 * pc_y[k] * qsh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, osh_440, qsg0_418, qsg0_419, \
                         qsg1_418, qsg1_419, qsh_585, qsh_586, \
                         qsh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * qsg0_418[k]
                   - f_7 * qsg1_418[k]
                   + f_3 * pc_y[k] * qsh_585[k];

        t_781[k] = f_4 * qsg0_419[k]
                   - f_5 * qsg1_419[k]
                   + f_3 * pc_y[k] * qsh_586[k];

        t_782[k] = f_3 * pc_y[k] * qsh_587[k];

        t_783[k] = f_23 * osh_440[k]
                   + f_1 * qsg0_419[k]
                   - f_2 * qsg1_419[k]
                   + f_3 * pc_z[k] * qsh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pc_x, pc_y, pc_z, osh_441, osh_588, \
                         osh_591, qsg0_420, qsg0_423, qsg1_420, qsg1_423, qsh_588, \
                         qsh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_22 * osh_588[k]
                   + f_1 * qsg0_420[k]
                   - f_2 * qsg1_420[k]
                   + f_3 * pc_x[k] * qsh_588[k];

        t_785[k] = f_21 * osh_441[k]
                   + f_3 * pc_y[k] * qsh_588[k];

        t_786[k] = f_3 * pc_z[k] * qsh_588[k];

        t_787[k] = f_22 * osh_591[k]
                   + f_8 * qsg0_423[k]
                   - f_9 * qsg1_423[k]
                   + f_3 * pc_x[k] * qsh_591[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, t_791, pc_x, pc_z, osh_594, qsg0_420, qsg0_426, \
                         qsg1_420, qsg1_426, qsh_589, qsh_590, qsh_591, \
                         qsh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_3 * pc_z[k] * qsh_589[k];

        t_789[k] = f_4 * qsg0_420[k]
                   - f_5 * qsg1_420[k]
                   + f_3 * pc_z[k] * qsh_590[k];

        t_790[k] = f_22 * osh_594[k]
                   + f_6 * qsg0_426[k]
                   - f_7 * qsg1_426[k]
                   + f_3 * pc_x[k] * qsh_594[k];

        t_791[k] = f_3 * pc_z[k] * qsh_591[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pc_x, pc_y, pc_z, osh_446, osh_598, \
                         qsg0_422, qsg0_430, qsg1_422, qsg1_430, qsh_593, qsh_594, \
                         qsh_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_21 * osh_446[k]
                   + f_3 * pc_y[k] * qsh_593[k];

        t_793[k] = f_6 * qsg0_422[k]
                   - f_7 * qsg1_422[k]
                   + f_3 * pc_z[k] * qsh_593[k];

        t_794[k] = f_22 * osh_598[k]
                   + f_4 * qsg0_430[k]
                   - f_5 * qsg1_430[k]
                   + f_3 * pc_x[k] * qsh_598[k];

        t_795[k] = f_3 * pc_z[k] * qsh_594[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pc_x, pc_y, pc_z, osh_450, osh_603, \
                         qsg0_423, qsg0_425, qsg1_423, qsg1_425, qsh_595, qsh_597, \
                         qsh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_4 * qsg0_423[k]
                   - f_5 * qsg1_423[k]
                   + f_3 * pc_z[k] * qsh_595[k];

        t_797[k] = f_21 * osh_450[k]
                   + f_3 * pc_y[k] * qsh_597[k];

        t_798[k] = f_8 * qsg0_425[k]
                   - f_9 * qsg1_425[k]
                   + f_3 * pc_z[k] * qsh_597[k];

        t_799[k] = f_22 * osh_603[k]
                   + f_3 * pc_x[k] * qsh_603[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, pc_x, pc_z, osh_605, osh_606, \
                         osh_607, osh_608, qsh_598, qsh_605, qsh_606, qsh_607, \
                         qsh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_3 * pc_z[k] * qsh_598[k];

        t_801[k] = f_22 * osh_605[k]
                   + f_3 * pc_x[k] * qsh_605[k];

        t_802[k] = f_22 * osh_606[k]
                   + f_3 * pc_x[k] * qsh_606[k];

        t_803[k] = f_22 * osh_607[k]
                   + f_3 * pc_x[k] * qsh_607[k];

        t_804[k] = f_22 * osh_608[k]
                   + f_3 * pc_x[k] * qsh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, pc_y, pc_z, osh_456, qsg0_430, qsg0_431, \
                         qsg1_430, qsg1_431, qsh_603, qsh_604, \
                         qsh_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_21 * osh_456[k]
                   + f_1 * qsg0_430[k]
                   - f_2 * qsg1_430[k]
                   + f_3 * pc_y[k] * qsh_603[k];

        t_806[k] = f_3 * pc_z[k] * qsh_603[k];

        t_807[k] = f_4 * qsg0_430[k]
                   - f_5 * qsg1_430[k]
                   + f_3 * pc_z[k] * qsh_604[k];

        t_808[k] = f_6 * qsg0_431[k]
                   - f_7 * qsg1_431[k]
                   + f_3 * pc_z[k] * qsh_605[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, pa_z, pc_y, pc_z, osi0_588, osh_461, \
                         osi1_588, qsg0_432, qsg0_434, qsg1_432, qsg1_434, qsh_606, \
                         qsh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_8 * qsg0_432[k]
                   - f_9 * qsg1_432[k]
                   + f_3 * pc_z[k] * qsh_606[k];

        t_810[k] = f_21 * osh_461[k]
                   + f_3 * pc_y[k] * qsh_608[k];

        t_811[k] = f_1 * qsg0_434[k]
                   - f_2 * qsg1_434[k]
                   + f_3 * pc_z[k] * qsh_608[k];

        t_812[k] = pa_z[k] * osi0_588[k]
                   - f_10 * pc_z[k] * osi1_588[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, t_816, pa_z, pc_y, pc_z, osi0_591, osh_441, \
                         osh_462, osh_464, osi1_591, qsh_609, qsh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_23 * osh_462[k]
                   + f_3 * pc_y[k] * qsh_609[k];

        t_814[k] = f_11 * osh_441[k]
                   + f_3 * pc_z[k] * qsh_609[k];

        t_815[k] = pa_z[k] * osi0_591[k]
                   - f_10 * pc_z[k] * osi1_591[k];

        t_816[k] = f_23 * osh_464[k]
                   + f_3 * pc_y[k] * qsh_611[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pa_z, pc_x, pc_z, osi0_594, osh_444, osh_614, \
                         osi1_594, qsg0_440, qsg1_440, qsh_612, \
                         qsh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_22 * osh_614[k]
                   + f_8 * qsg0_440[k]
                   - f_9 * qsg1_440[k]
                   + f_3 * pc_x[k] * qsh_614[k];

        t_818[k] = pa_z[k] * osi0_594[k]
                   - f_10 * pc_z[k] * osi1_594[k];

        t_819[k] = f_11 * osh_444[k]
                   + f_3 * pc_z[k] * qsh_612[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_598 = buffer.data(osi0 + 598);
    const auto *osi0_600 = buffer.data(osi0 + 600);
    const auto *osi0_609 = buffer.data(osi0 + 609);

    const auto *osh_447 = buffer.data(osh + 447);
    const auto *osh_448 = buffer.data(osh + 448);
    const auto *osh_456 = buffer.data(osh + 456);
    const auto *osh_461 = buffer.data(osh + 461);
    const auto *osh_462 = buffer.data(osh + 462);
    const auto *osh_465 = buffer.data(osh + 465);
    const auto *osh_467 = buffer.data(osh + 467);
    const auto *osh_468 = buffer.data(osh + 468);
    const auto *osh_471 = buffer.data(osh + 471);
    const auto *osh_477 = buffer.data(osh + 477);
    const auto *osh_479 = buffer.data(osh + 479);
    const auto *osh_480 = buffer.data(osh + 480);
    const auto *osh_481 = buffer.data(osh + 481);
    const auto *osh_482 = buffer.data(osh + 482);
    const auto *osh_483 = buffer.data(osh + 483);
    const auto *osh_485 = buffer.data(osh + 485);
    const auto *osh_486 = buffer.data(osh + 486);
    const auto *osh_488 = buffer.data(osh + 488);
    const auto *osh_489 = buffer.data(osh + 489);
    const auto *osh_492 = buffer.data(osh + 492);
    const auto *osh_498 = buffer.data(osh + 498);
    const auto *osh_500 = buffer.data(osh + 500);
    const auto *osh_501 = buffer.data(osh + 501);
    const auto *osh_502 = buffer.data(osh + 502);
    const auto *osh_503 = buffer.data(osh + 503);
    const auto *osh_504 = buffer.data(osh + 504);
    const auto *osh_506 = buffer.data(osh + 506);
    const auto *osh_507 = buffer.data(osh + 507);
    const auto *osh_509 = buffer.data(osh + 509);
    const auto *osh_510 = buffer.data(osh + 510);
    const auto *osh_513 = buffer.data(osh + 513);
    const auto *osh_519 = buffer.data(osh + 519);
    const auto *osh_521 = buffer.data(osh + 521);
    const auto *osh_522 = buffer.data(osh + 522);
    const auto *osh_523 = buffer.data(osh + 523);
    const auto *osh_524 = buffer.data(osh + 524);
    const auto *osh_525 = buffer.data(osh + 525);
    const auto *osh_527 = buffer.data(osh + 527);
    const auto *osh_530 = buffer.data(osh + 530);
    const auto *osh_534 = buffer.data(osh + 534);
    const auto *osh_540 = buffer.data(osh + 540);
    const auto *osh_542 = buffer.data(osh + 542);
    const auto *osh_543 = buffer.data(osh + 543);
    const auto *osh_544 = buffer.data(osh + 544);
    const auto *osh_545 = buffer.data(osh + 545);
    const auto *osh_618 = buffer.data(osh + 618);
    const auto *osh_623 = buffer.data(osh + 623);
    const auto *osh_624 = buffer.data(osh + 624);
    const auto *osh_625 = buffer.data(osh + 625);
    const auto *osh_626 = buffer.data(osh + 626);
    const auto *osh_627 = buffer.data(osh + 627);
    const auto *osh_628 = buffer.data(osh + 628);
    const auto *osh_629 = buffer.data(osh + 629);
    const auto *osh_630 = buffer.data(osh + 630);
    const auto *osh_633 = buffer.data(osh + 633);
    const auto *osh_635 = buffer.data(osh + 635);
    const auto *osh_636 = buffer.data(osh + 636);
    const auto *osh_639 = buffer.data(osh + 639);
    const auto *osh_640 = buffer.data(osh + 640);
    const auto *osh_642 = buffer.data(osh + 642);
    const auto *osh_644 = buffer.data(osh + 644);
    const auto *osh_645 = buffer.data(osh + 645);
    const auto *osh_646 = buffer.data(osh + 646);
    const auto *osh_647 = buffer.data(osh + 647);
    const auto *osh_648 = buffer.data(osh + 648);
    const auto *osh_649 = buffer.data(osh + 649);
    const auto *osh_650 = buffer.data(osh + 650);
    const auto *osh_651 = buffer.data(osh + 651);
    const auto *osh_654 = buffer.data(osh + 654);
    const auto *osh_656 = buffer.data(osh + 656);
    const auto *osh_657 = buffer.data(osh + 657);
    const auto *osh_660 = buffer.data(osh + 660);
    const auto *osh_661 = buffer.data(osh + 661);
    const auto *osh_663 = buffer.data(osh + 663);
    const auto *osh_665 = buffer.data(osh + 665);
    const auto *osh_666 = buffer.data(osh + 666);
    const auto *osh_667 = buffer.data(osh + 667);
    const auto *osh_668 = buffer.data(osh + 668);
    const auto *osh_669 = buffer.data(osh + 669);
    const auto *osh_670 = buffer.data(osh + 670);
    const auto *osh_671 = buffer.data(osh + 671);
    const auto *osh_672 = buffer.data(osh + 672);
    const auto *osh_675 = buffer.data(osh + 675);
    const auto *osh_677 = buffer.data(osh + 677);
    const auto *osh_678 = buffer.data(osh + 678);
    const auto *osh_681 = buffer.data(osh + 681);
    const auto *osh_682 = buffer.data(osh + 682);
    const auto *osh_684 = buffer.data(osh + 684);
    const auto *osh_686 = buffer.data(osh + 686);
    const auto *osh_687 = buffer.data(osh + 687);
    const auto *osh_688 = buffer.data(osh + 688);
    const auto *osh_689 = buffer.data(osh + 689);
    const auto *osh_690 = buffer.data(osh + 690);
    const auto *osh_691 = buffer.data(osh + 691);
    const auto *osh_692 = buffer.data(osh + 692);
    const auto *osh_693 = buffer.data(osh + 693);

    const auto *osi1_598 = buffer.data(osi1 + 598);
    const auto *osi1_600 = buffer.data(osi1 + 600);
    const auto *osi1_609 = buffer.data(osi1 + 609);

    const auto *qsg0_444 = buffer.data(qsg0 + 444);
    const auto *qsg0_447 = buffer.data(qsg0 + 447);
    const auto *qsg0_448 = buffer.data(qsg0 + 448);
    const auto *qsg0_449 = buffer.data(qsg0 + 449);
    const auto *qsg0_450 = buffer.data(qsg0 + 450);
    const auto *qsg0_453 = buffer.data(qsg0 + 453);
    const auto *qsg0_455 = buffer.data(qsg0 + 455);
    const auto *qsg0_456 = buffer.data(qsg0 + 456);
    const auto *qsg0_459 = buffer.data(qsg0 + 459);
    const auto *qsg0_460 = buffer.data(qsg0 + 460);
    const auto *qsg0_462 = buffer.data(qsg0 + 462);
    const auto *qsg0_463 = buffer.data(qsg0 + 463);
    const auto *qsg0_464 = buffer.data(qsg0 + 464);
    const auto *qsg0_465 = buffer.data(qsg0 + 465);
    const auto *qsg0_468 = buffer.data(qsg0 + 468);
    const auto *qsg0_470 = buffer.data(qsg0 + 470);
    const auto *qsg0_471 = buffer.data(qsg0 + 471);
    const auto *qsg0_474 = buffer.data(qsg0 + 474);
    const auto *qsg0_475 = buffer.data(qsg0 + 475);
    const auto *qsg0_477 = buffer.data(qsg0 + 477);
    const auto *qsg0_478 = buffer.data(qsg0 + 478);
    const auto *qsg0_479 = buffer.data(qsg0 + 479);
    const auto *qsg0_480 = buffer.data(qsg0 + 480);
    const auto *qsg0_483 = buffer.data(qsg0 + 483);
    const auto *qsg0_485 = buffer.data(qsg0 + 485);
    const auto *qsg0_486 = buffer.data(qsg0 + 486);
    const auto *qsg0_489 = buffer.data(qsg0 + 489);
    const auto *qsg0_490 = buffer.data(qsg0 + 490);
    const auto *qsg0_492 = buffer.data(qsg0 + 492);
    const auto *qsg0_493 = buffer.data(qsg0 + 493);
    const auto *qsg0_494 = buffer.data(qsg0 + 494);
    const auto *qsg0_495 = buffer.data(qsg0 + 495);

    const auto *qsg1_444 = buffer.data(qsg1 + 444);
    const auto *qsg1_447 = buffer.data(qsg1 + 447);
    const auto *qsg1_448 = buffer.data(qsg1 + 448);
    const auto *qsg1_449 = buffer.data(qsg1 + 449);
    const auto *qsg1_450 = buffer.data(qsg1 + 450);
    const auto *qsg1_453 = buffer.data(qsg1 + 453);
    const auto *qsg1_455 = buffer.data(qsg1 + 455);
    const auto *qsg1_456 = buffer.data(qsg1 + 456);
    const auto *qsg1_459 = buffer.data(qsg1 + 459);
    const auto *qsg1_460 = buffer.data(qsg1 + 460);
    const auto *qsg1_462 = buffer.data(qsg1 + 462);
    const auto *qsg1_463 = buffer.data(qsg1 + 463);
    const auto *qsg1_464 = buffer.data(qsg1 + 464);
    const auto *qsg1_465 = buffer.data(qsg1 + 465);
    const auto *qsg1_468 = buffer.data(qsg1 + 468);
    const auto *qsg1_470 = buffer.data(qsg1 + 470);
    const auto *qsg1_471 = buffer.data(qsg1 + 471);
    const auto *qsg1_474 = buffer.data(qsg1 + 474);
    const auto *qsg1_475 = buffer.data(qsg1 + 475);
    const auto *qsg1_477 = buffer.data(qsg1 + 477);
    const auto *qsg1_478 = buffer.data(qsg1 + 478);
    const auto *qsg1_479 = buffer.data(qsg1 + 479);
    const auto *qsg1_480 = buffer.data(qsg1 + 480);
    const auto *qsg1_483 = buffer.data(qsg1 + 483);
    const auto *qsg1_485 = buffer.data(qsg1 + 485);
    const auto *qsg1_486 = buffer.data(qsg1 + 486);
    const auto *qsg1_489 = buffer.data(qsg1 + 489);
    const auto *qsg1_490 = buffer.data(qsg1 + 490);
    const auto *qsg1_492 = buffer.data(qsg1 + 492);
    const auto *qsg1_493 = buffer.data(qsg1 + 493);
    const auto *qsg1_494 = buffer.data(qsg1 + 494);
    const auto *qsg1_495 = buffer.data(qsg1 + 495);

    const auto *qsh_614 = buffer.data(qsh + 614);
    const auto *qsh_615 = buffer.data(qsh + 615);
    const auto *qsh_618 = buffer.data(qsh + 618);
    const auto *qsh_623 = buffer.data(qsh + 623);
    const auto *qsh_624 = buffer.data(qsh + 624);
    const auto *qsh_625 = buffer.data(qsh + 625);
    const auto *qsh_626 = buffer.data(qsh + 626);
    const auto *qsh_627 = buffer.data(qsh + 627);
    const auto *qsh_628 = buffer.data(qsh + 628);
    const auto *qsh_629 = buffer.data(qsh + 629);
    const auto *qsh_630 = buffer.data(qsh + 630);
    const auto *qsh_632 = buffer.data(qsh + 632);
    const auto *qsh_633 = buffer.data(qsh + 633);
    const auto *qsh_635 = buffer.data(qsh + 635);
    const auto *qsh_636 = buffer.data(qsh + 636);
    const auto *qsh_639 = buffer.data(qsh + 639);
    const auto *qsh_640 = buffer.data(qsh + 640);
    const auto *qsh_642 = buffer.data(qsh + 642);
    const auto *qsh_644 = buffer.data(qsh + 644);
    const auto *qsh_645 = buffer.data(qsh + 645);
    const auto *qsh_646 = buffer.data(qsh + 646);
    const auto *qsh_647 = buffer.data(qsh + 647);
    const auto *qsh_648 = buffer.data(qsh + 648);
    const auto *qsh_649 = buffer.data(qsh + 649);
    const auto *qsh_650 = buffer.data(qsh + 650);
    const auto *qsh_651 = buffer.data(qsh + 651);
    const auto *qsh_653 = buffer.data(qsh + 653);
    const auto *qsh_654 = buffer.data(qsh + 654);
    const auto *qsh_656 = buffer.data(qsh + 656);
    const auto *qsh_657 = buffer.data(qsh + 657);
    const auto *qsh_660 = buffer.data(qsh + 660);
    const auto *qsh_661 = buffer.data(qsh + 661);
    const auto *qsh_663 = buffer.data(qsh + 663);
    const auto *qsh_665 = buffer.data(qsh + 665);
    const auto *qsh_666 = buffer.data(qsh + 666);
    const auto *qsh_667 = buffer.data(qsh + 667);
    const auto *qsh_668 = buffer.data(qsh + 668);
    const auto *qsh_669 = buffer.data(qsh + 669);
    const auto *qsh_670 = buffer.data(qsh + 670);
    const auto *qsh_671 = buffer.data(qsh + 671);
    const auto *qsh_672 = buffer.data(qsh + 672);
    const auto *qsh_674 = buffer.data(qsh + 674);
    const auto *qsh_675 = buffer.data(qsh + 675);
    const auto *qsh_677 = buffer.data(qsh + 677);
    const auto *qsh_678 = buffer.data(qsh + 678);
    const auto *qsh_681 = buffer.data(qsh + 681);
    const auto *qsh_682 = buffer.data(qsh + 682);
    const auto *qsh_684 = buffer.data(qsh + 684);
    const auto *qsh_686 = buffer.data(qsh + 686);
    const auto *qsh_687 = buffer.data(qsh + 687);
    const auto *qsh_688 = buffer.data(qsh + 688);
    const auto *qsh_689 = buffer.data(qsh + 689);
    const auto *qsh_690 = buffer.data(qsh + 690);
    const auto *qsh_691 = buffer.data(qsh + 691);
    const auto *qsh_692 = buffer.data(qsh + 692);
    const auto *qsh_693 = buffer.data(qsh + 693);

#pragma omp simd aligned(t_820, t_821, t_822, pa_z, pc_x, pc_y, pc_z, osi0_598, osh_467, \
                         osh_618, osi1_598, qsg0_444, qsg1_444, qsh_614, \
                         qsh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_23 * osh_467[k]
                   + f_3 * pc_y[k] * qsh_614[k];

        t_821[k] = f_22 * osh_618[k]
                   + f_6 * qsg0_444[k]
                   - f_7 * qsg1_444[k]
                   + f_3 * pc_x[k] * qsh_618[k];

        t_822[k] = pa_z[k] * osi0_598[k]
                   - f_10 * pc_z[k] * osi1_598[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pa_z, pc_y, pc_z, osi0_600, osh_447, osh_448, \
                         osh_471, osi1_600, qsh_615, qsh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_11 * osh_447[k]
                   + f_3 * pc_z[k] * qsh_615[k];

        t_824[k] = pa_z[k] * osi0_600[k]
                   + f_12 * osh_448[k]
                   - f_10 * pc_z[k] * osi1_600[k];

        t_825[k] = f_23 * osh_471[k]
                   + f_3 * pc_y[k] * qsh_618[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, osh_623, osh_624, osh_625, osh_626, \
                         qsg0_449, qsg1_449, qsh_623, qsh_624, qsh_625, \
                         qsh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_22 * osh_623[k]
                   + f_4 * qsg0_449[k]
                   - f_5 * qsg1_449[k]
                   + f_3 * pc_x[k] * qsh_623[k];

        t_827[k] = f_22 * osh_624[k]
                   + f_3 * pc_x[k] * qsh_624[k];

        t_828[k] = f_22 * osh_625[k]
                   + f_3 * pc_x[k] * qsh_625[k];

        t_829[k] = f_22 * osh_626[k]
                   + f_3 * pc_x[k] * qsh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pa_z, pc_x, pc_z, osi0_609, osh_627, \
                         osh_628, osh_629, osi1_609, qsh_627, qsh_628, \
                         qsh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_22 * osh_627[k]
                   + f_3 * pc_x[k] * qsh_627[k];

        t_831[k] = f_22 * osh_628[k]
                   + f_3 * pc_x[k] * qsh_628[k];

        t_832[k] = f_22 * osh_629[k]
                   + f_3 * pc_x[k] * qsh_629[k];

        t_833[k] = pa_z[k] * osi0_609[k]
                   - f_10 * pc_z[k] * osi1_609[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_y, pc_z, osh_456, osh_479, osh_480, qsg0_447, \
                         qsg0_448, qsg1_447, qsg1_448, qsh_624, qsh_626, \
                         qsh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_11 * osh_456[k]
                   + f_3 * pc_z[k] * qsh_624[k];

        t_835[k] = f_23 * osh_479[k]
                   + f_8 * qsg0_447[k]
                   - f_9 * qsg1_447[k]
                   + f_3 * pc_y[k] * qsh_626[k];

        t_836[k] = f_23 * osh_480[k]
                   + f_6 * qsg0_448[k]
                   - f_7 * qsg1_448[k]
                   + f_3 * pc_y[k] * qsh_627[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, osh_461, osh_481, osh_482, qsg0_449, \
                         qsg1_449, qsh_628, qsh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_23 * osh_481[k]
                   + f_4 * qsg0_449[k]
                   - f_5 * qsg1_449[k]
                   + f_3 * pc_y[k] * qsh_628[k];

        t_838[k] = f_23 * osh_482[k]
                   + f_3 * pc_y[k] * qsh_629[k];

        t_839[k] = f_11 * osh_461[k]
                   + f_1 * qsg0_449[k]
                   - f_2 * qsg1_449[k]
                   + f_3 * pc_z[k] * qsh_629[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, pc_y, pc_z, osh_462, osh_483, osh_630, \
                         qsg0_450, qsg1_450, qsh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_22 * osh_630[k]
                   + f_1 * qsg0_450[k]
                   - f_2 * qsg1_450[k]
                   + f_3 * pc_x[k] * qsh_630[k];

        t_841[k] = f_22 * osh_483[k]
                   + f_3 * pc_y[k] * qsh_630[k];

        t_842[k] = f_12 * osh_462[k]
                   + f_3 * pc_z[k] * qsh_630[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, pc_y, osh_485, osh_633, osh_635, qsg0_453, \
                         qsg0_455, qsg1_453, qsg1_455, qsh_632, qsh_633, \
                         qsh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_22 * osh_633[k]
                   + f_8 * qsg0_453[k]
                   - f_9 * qsg1_453[k]
                   + f_3 * pc_x[k] * qsh_633[k];

        t_844[k] = f_22 * osh_485[k]
                   + f_3 * pc_y[k] * qsh_632[k];

        t_845[k] = f_22 * osh_635[k]
                   + f_8 * qsg0_455[k]
                   - f_9 * qsg1_455[k]
                   + f_3 * pc_x[k] * qsh_635[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pc_x, pc_y, pc_z, osh_465, osh_488, osh_636, \
                         qsg0_456, qsg1_456, qsh_633, qsh_635, \
                         qsh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_22 * osh_636[k]
                   + f_6 * qsg0_456[k]
                   - f_7 * qsg1_456[k]
                   + f_3 * pc_x[k] * qsh_636[k];

        t_847[k] = f_12 * osh_465[k]
                   + f_3 * pc_z[k] * qsh_633[k];

        t_848[k] = f_22 * osh_488[k]
                   + f_3 * pc_y[k] * qsh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, osh_468, osh_639, osh_640, qsg0_459, \
                         qsg0_460, qsg1_459, qsg1_460, qsh_636, qsh_639, \
                         qsh_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_22 * osh_639[k]
                   + f_6 * qsg0_459[k]
                   - f_7 * qsg1_459[k]
                   + f_3 * pc_x[k] * qsh_639[k];

        t_850[k] = f_22 * osh_640[k]
                   + f_4 * qsg0_460[k]
                   - f_5 * qsg1_460[k]
                   + f_3 * pc_x[k] * qsh_640[k];

        t_851[k] = f_12 * osh_468[k]
                   + f_3 * pc_z[k] * qsh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, pc_x, pc_y, osh_492, osh_642, osh_644, qsg0_462, \
                         qsg0_464, qsg1_462, qsg1_464, qsh_639, qsh_642, \
                         qsh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_22 * osh_642[k]
                   + f_4 * qsg0_462[k]
                   - f_5 * qsg1_462[k]
                   + f_3 * pc_x[k] * qsh_642[k];

        t_853[k] = f_22 * osh_492[k]
                   + f_3 * pc_y[k] * qsh_639[k];

        t_854[k] = f_22 * osh_644[k]
                   + f_4 * qsg0_464[k]
                   - f_5 * qsg1_464[k]
                   + f_3 * pc_x[k] * qsh_644[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, pc_x, osh_645, osh_646, osh_647, \
                         osh_648, osh_649, qsh_645, qsh_646, qsh_647, qsh_648, \
                         qsh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_22 * osh_645[k]
                   + f_3 * pc_x[k] * qsh_645[k];

        t_856[k] = f_22 * osh_646[k]
                   + f_3 * pc_x[k] * qsh_646[k];

        t_857[k] = f_22 * osh_647[k]
                   + f_3 * pc_x[k] * qsh_647[k];

        t_858[k] = f_22 * osh_648[k]
                   + f_3 * pc_x[k] * qsh_648[k];

        t_859[k] = f_22 * osh_649[k]
                   + f_3 * pc_x[k] * qsh_649[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pc_x, pc_y, pc_z, osh_477, osh_498, osh_650, \
                         qsg0_460, qsg1_460, qsh_645, qsh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_22 * osh_650[k]
                   + f_3 * pc_x[k] * qsh_650[k];

        t_861[k] = f_22 * osh_498[k]
                   + f_1 * qsg0_460[k]
                   - f_2 * qsg1_460[k]
                   + f_3 * pc_y[k] * qsh_645[k];

        t_862[k] = f_12 * osh_477[k]
                   + f_3 * pc_z[k] * qsh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_y, osh_500, osh_501, osh_502, qsg0_462, \
                         qsg0_463, qsg0_464, qsg1_462, qsg1_463, qsg1_464, qsh_647, qsh_648, \
                         qsh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_22 * osh_500[k]
                   + f_8 * qsg0_462[k]
                   - f_9 * qsg1_462[k]
                   + f_3 * pc_y[k] * qsh_647[k];

        t_864[k] = f_22 * osh_501[k]
                   + f_6 * qsg0_463[k]
                   - f_7 * qsg1_463[k]
                   + f_3 * pc_y[k] * qsh_648[k];

        t_865[k] = f_22 * osh_502[k]
                   + f_4 * qsg0_464[k]
                   - f_5 * qsg1_464[k]
                   + f_3 * pc_y[k] * qsh_649[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, osh_482, osh_503, osh_651, \
                         qsg0_464, qsg0_465, qsg1_464, qsg1_465, qsh_650, \
                         qsh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_22 * osh_503[k]
                   + f_3 * pc_y[k] * qsh_650[k];

        t_867[k] = f_12 * osh_482[k]
                   + f_1 * qsg0_464[k]
                   - f_2 * qsg1_464[k]
                   + f_3 * pc_z[k] * qsh_650[k];

        t_868[k] = f_22 * osh_651[k]
                   + f_1 * qsg0_465[k]
                   - f_2 * qsg1_465[k]
                   + f_3 * pc_x[k] * qsh_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, pc_z, osh_483, osh_504, \
                         osh_506, osh_654, qsg0_468, qsg1_468, qsh_651, qsh_653, \
                         qsh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * osh_504[k]
                   + f_3 * pc_y[k] * qsh_651[k];

        t_870[k] = f_13 * osh_483[k]
                   + f_3 * pc_z[k] * qsh_651[k];

        t_871[k] = f_22 * osh_654[k]
                   + f_8 * qsg0_468[k]
                   - f_9 * qsg1_468[k]
                   + f_3 * pc_x[k] * qsh_654[k];

        t_872[k] = f_14 * osh_506[k]
                   + f_3 * pc_y[k] * qsh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_z, osh_486, osh_656, osh_657, qsg0_470, \
                         qsg0_471, qsg1_470, qsg1_471, qsh_654, qsh_656, \
                         qsh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_22 * osh_656[k]
                   + f_8 * qsg0_470[k]
                   - f_9 * qsg1_470[k]
                   + f_3 * pc_x[k] * qsh_656[k];

        t_874[k] = f_22 * osh_657[k]
                   + f_6 * qsg0_471[k]
                   - f_7 * qsg1_471[k]
                   + f_3 * pc_x[k] * qsh_657[k];

        t_875[k] = f_13 * osh_486[k]
                   + f_3 * pc_z[k] * qsh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_y, osh_509, osh_660, osh_661, qsg0_474, \
                         qsg0_475, qsg1_474, qsg1_475, qsh_656, qsh_660, \
                         qsh_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * osh_509[k]
                   + f_3 * pc_y[k] * qsh_656[k];

        t_877[k] = f_22 * osh_660[k]
                   + f_6 * qsg0_474[k]
                   - f_7 * qsg1_474[k]
                   + f_3 * pc_x[k] * qsh_660[k];

        t_878[k] = f_22 * osh_661[k]
                   + f_4 * qsg0_475[k]
                   - f_5 * qsg1_475[k]
                   + f_3 * pc_x[k] * qsh_661[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, osh_489, osh_513, osh_663, \
                         qsg0_477, qsg1_477, qsh_657, qsh_660, \
                         qsh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * osh_489[k]
                   + f_3 * pc_z[k] * qsh_657[k];

        t_880[k] = f_22 * osh_663[k]
                   + f_4 * qsg0_477[k]
                   - f_5 * qsg1_477[k]
                   + f_3 * pc_x[k] * qsh_663[k];

        t_881[k] = f_14 * osh_513[k]
                   + f_3 * pc_y[k] * qsh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pc_x, osh_665, osh_666, osh_667, osh_668, \
                         qsg0_479, qsg1_479, qsh_665, qsh_666, qsh_667, \
                         qsh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_22 * osh_665[k]
                   + f_4 * qsg0_479[k]
                   - f_5 * qsg1_479[k]
                   + f_3 * pc_x[k] * qsh_665[k];

        t_883[k] = f_22 * osh_666[k]
                   + f_3 * pc_x[k] * qsh_666[k];

        t_884[k] = f_22 * osh_667[k]
                   + f_3 * pc_x[k] * qsh_667[k];

        t_885[k] = f_22 * osh_668[k]
                   + f_3 * pc_x[k] * qsh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pc_x, pc_y, osh_519, osh_669, osh_670, \
                         osh_671, qsg0_475, qsg1_475, qsh_666, qsh_669, qsh_670, \
                         qsh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_22 * osh_669[k]
                   + f_3 * pc_x[k] * qsh_669[k];

        t_887[k] = f_22 * osh_670[k]
                   + f_3 * pc_x[k] * qsh_670[k];

        t_888[k] = f_22 * osh_671[k]
                   + f_3 * pc_x[k] * qsh_671[k];

        t_889[k] = f_14 * osh_519[k]
                   + f_1 * qsg0_475[k]
                   - f_2 * qsg1_475[k]
                   + f_3 * pc_y[k] * qsh_666[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pc_y, pc_z, osh_498, osh_521, osh_522, qsg0_477, \
                         qsg0_478, qsg1_477, qsg1_478, qsh_666, qsh_668, \
                         qsh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * osh_498[k]
                   + f_3 * pc_z[k] * qsh_666[k];

        t_891[k] = f_14 * osh_521[k]
                   + f_8 * qsg0_477[k]
                   - f_9 * qsg1_477[k]
                   + f_3 * pc_y[k] * qsh_668[k];

        t_892[k] = f_14 * osh_522[k]
                   + f_6 * qsg0_478[k]
                   - f_7 * qsg1_478[k]
                   + f_3 * pc_y[k] * qsh_669[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, pc_z, osh_503, osh_523, osh_524, qsg0_479, \
                         qsg1_479, qsh_670, qsh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_14 * osh_523[k]
                   + f_4 * qsg0_479[k]
                   - f_5 * qsg1_479[k]
                   + f_3 * pc_y[k] * qsh_670[k];

        t_894[k] = f_14 * osh_524[k]
                   + f_3 * pc_y[k] * qsh_671[k];

        t_895[k] = f_13 * osh_503[k]
                   + f_1 * qsg0_479[k]
                   - f_2 * qsg1_479[k]
                   + f_3 * pc_z[k] * qsh_671[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_x, pc_y, pc_z, osh_504, osh_525, osh_672, \
                         qsg0_480, qsg1_480, qsh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_22 * osh_672[k]
                   + f_1 * qsg0_480[k]
                   - f_2 * qsg1_480[k]
                   + f_3 * pc_x[k] * qsh_672[k];

        t_897[k] = f_13 * osh_525[k]
                   + f_3 * pc_y[k] * qsh_672[k];

        t_898[k] = f_14 * osh_504[k]
                   + f_3 * pc_z[k] * qsh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, osh_527, osh_675, osh_677, qsg0_483, \
                         qsg0_485, qsg1_483, qsg1_485, qsh_674, qsh_675, \
                         qsh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_22 * osh_675[k]
                   + f_8 * qsg0_483[k]
                   - f_9 * qsg1_483[k]
                   + f_3 * pc_x[k] * qsh_675[k];

        t_900[k] = f_13 * osh_527[k]
                   + f_3 * pc_y[k] * qsh_674[k];

        t_901[k] = f_22 * osh_677[k]
                   + f_8 * qsg0_485[k]
                   - f_9 * qsg1_485[k]
                   + f_3 * pc_x[k] * qsh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, osh_507, osh_530, osh_678, \
                         qsg0_486, qsg1_486, qsh_675, qsh_677, \
                         qsh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_22 * osh_678[k]
                   + f_6 * qsg0_486[k]
                   - f_7 * qsg1_486[k]
                   + f_3 * pc_x[k] * qsh_678[k];

        t_903[k] = f_14 * osh_507[k]
                   + f_3 * pc_z[k] * qsh_675[k];

        t_904[k] = f_13 * osh_530[k]
                   + f_3 * pc_y[k] * qsh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, osh_510, osh_681, osh_682, qsg0_489, \
                         qsg0_490, qsg1_489, qsg1_490, qsh_678, qsh_681, \
                         qsh_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_22 * osh_681[k]
                   + f_6 * qsg0_489[k]
                   - f_7 * qsg1_489[k]
                   + f_3 * pc_x[k] * qsh_681[k];

        t_906[k] = f_22 * osh_682[k]
                   + f_4 * qsg0_490[k]
                   - f_5 * qsg1_490[k]
                   + f_3 * pc_x[k] * qsh_682[k];

        t_907[k] = f_14 * osh_510[k]
                   + f_3 * pc_z[k] * qsh_678[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, osh_534, osh_684, osh_686, qsg0_492, \
                         qsg0_494, qsg1_492, qsg1_494, qsh_681, qsh_684, \
                         qsh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_22 * osh_684[k]
                   + f_4 * qsg0_492[k]
                   - f_5 * qsg1_492[k]
                   + f_3 * pc_x[k] * qsh_684[k];

        t_909[k] = f_13 * osh_534[k]
                   + f_3 * pc_y[k] * qsh_681[k];

        t_910[k] = f_22 * osh_686[k]
                   + f_4 * qsg0_494[k]
                   - f_5 * qsg1_494[k]
                   + f_3 * pc_x[k] * qsh_686[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pc_x, osh_687, osh_688, osh_689, \
                         osh_690, osh_691, qsh_687, qsh_688, qsh_689, qsh_690, \
                         qsh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_22 * osh_687[k]
                   + f_3 * pc_x[k] * qsh_687[k];

        t_912[k] = f_22 * osh_688[k]
                   + f_3 * pc_x[k] * qsh_688[k];

        t_913[k] = f_22 * osh_689[k]
                   + f_3 * pc_x[k] * qsh_689[k];

        t_914[k] = f_22 * osh_690[k]
                   + f_3 * pc_x[k] * qsh_690[k];

        t_915[k] = f_22 * osh_691[k]
                   + f_3 * pc_x[k] * qsh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pc_x, pc_y, pc_z, osh_519, osh_540, osh_692, \
                         qsg0_490, qsg1_490, qsh_687, qsh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_22 * osh_692[k]
                   + f_3 * pc_x[k] * qsh_692[k];

        t_917[k] = f_13 * osh_540[k]
                   + f_1 * qsg0_490[k]
                   - f_2 * qsg1_490[k]
                   + f_3 * pc_y[k] * qsh_687[k];

        t_918[k] = f_14 * osh_519[k]
                   + f_3 * pc_z[k] * qsh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pc_y, osh_542, osh_543, osh_544, qsg0_492, \
                         qsg0_493, qsg0_494, qsg1_492, qsg1_493, qsg1_494, qsh_689, qsh_690, \
                         qsh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_13 * osh_542[k]
                   + f_8 * qsg0_492[k]
                   - f_9 * qsg1_492[k]
                   + f_3 * pc_y[k] * qsh_689[k];

        t_920[k] = f_13 * osh_543[k]
                   + f_6 * qsg0_493[k]
                   - f_7 * qsg1_493[k]
                   + f_3 * pc_y[k] * qsh_690[k];

        t_921[k] = f_13 * osh_544[k]
                   + f_4 * qsg0_494[k]
                   - f_5 * qsg1_494[k]
                   + f_3 * pc_y[k] * qsh_691[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pc_x, pc_y, pc_z, osh_524, osh_545, osh_693, \
                         qsg0_494, qsg0_495, qsg1_494, qsg1_495, qsh_692, \
                         qsh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * osh_545[k]
                   + f_3 * pc_y[k] * qsh_692[k];

        t_923[k] = f_14 * osh_524[k]
                   + f_1 * qsg0_494[k]
                   - f_2 * qsg1_494[k]
                   + f_3 * pc_z[k] * qsh_692[k];

        t_924[k] = f_22 * osh_693[k]
                   + f_1 * qsg0_495[k]
                   - f_2 * qsg1_495[k]
                   + f_3 * pc_x[k] * qsh_693[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_756 = buffer.data(osi0 + 756);
    const auto *osi0_759 = buffer.data(osi0 + 759);
    const auto *osi0_761 = buffer.data(osi0 + 761);
    const auto *osi0_762 = buffer.data(osi0 + 762);
    const auto *osi0_765 = buffer.data(osi0 + 765);
    const auto *osi0_766 = buffer.data(osi0 + 766);
    const auto *osi0_768 = buffer.data(osi0 + 768);
    const auto *osi0_770 = buffer.data(osi0 + 770);
    const auto *osi0_783 = buffer.data(osi0 + 783);
    const auto *osi0_784 = buffer.data(osi0 + 784);
    const auto *osi0_787 = buffer.data(osi0 + 787);
    const auto *osi0_790 = buffer.data(osi0 + 790);

    const auto *osh_525 = buffer.data(osh + 525);
    const auto *osh_528 = buffer.data(osh + 528);
    const auto *osh_531 = buffer.data(osh + 531);
    const auto *osh_540 = buffer.data(osh + 540);
    const auto *osh_545 = buffer.data(osh + 545);
    const auto *osh_546 = buffer.data(osh + 546);
    const auto *osh_548 = buffer.data(osh + 548);
    const auto *osh_549 = buffer.data(osh + 549);
    const auto *osh_551 = buffer.data(osh + 551);
    const auto *osh_552 = buffer.data(osh + 552);
    const auto *osh_555 = buffer.data(osh + 555);
    const auto *osh_561 = buffer.data(osh + 561);
    const auto *osh_563 = buffer.data(osh + 563);
    const auto *osh_564 = buffer.data(osh + 564);
    const auto *osh_565 = buffer.data(osh + 565);
    const auto *osh_566 = buffer.data(osh + 566);
    const auto *osh_567 = buffer.data(osh + 567);
    const auto *osh_568 = buffer.data(osh + 568);
    const auto *osh_569 = buffer.data(osh + 569);
    const auto *osh_570 = buffer.data(osh + 570);
    const auto *osh_572 = buffer.data(osh + 572);
    const auto *osh_573 = buffer.data(osh + 573);
    const auto *osh_575 = buffer.data(osh + 575);
    const auto *osh_576 = buffer.data(osh + 576);
    const auto *osh_582 = buffer.data(osh + 582);
    const auto *osh_584 = buffer.data(osh + 584);
    const auto *osh_585 = buffer.data(osh + 585);
    const auto *osh_586 = buffer.data(osh + 586);
    const auto *osh_587 = buffer.data(osh + 587);
    const auto *osh_588 = buffer.data(osh + 588);
    const auto *osh_591 = buffer.data(osh + 591);
    const auto *osh_593 = buffer.data(osh + 593);
    const auto *osh_597 = buffer.data(osh + 597);
    const auto *osh_603 = buffer.data(osh + 603);
    const auto *osh_608 = buffer.data(osh + 608);
    const auto *osh_609 = buffer.data(osh + 609);
    const auto *osh_611 = buffer.data(osh + 611);
    const auto *osh_696 = buffer.data(osh + 696);
    const auto *osh_698 = buffer.data(osh + 698);
    const auto *osh_699 = buffer.data(osh + 699);
    const auto *osh_702 = buffer.data(osh + 702);
    const auto *osh_703 = buffer.data(osh + 703);
    const auto *osh_705 = buffer.data(osh + 705);
    const auto *osh_707 = buffer.data(osh + 707);
    const auto *osh_708 = buffer.data(osh + 708);
    const auto *osh_709 = buffer.data(osh + 709);
    const auto *osh_710 = buffer.data(osh + 710);
    const auto *osh_711 = buffer.data(osh + 711);
    const auto *osh_712 = buffer.data(osh + 712);
    const auto *osh_713 = buffer.data(osh + 713);
    const auto *osh_729 = buffer.data(osh + 729);
    const auto *osh_730 = buffer.data(osh + 730);
    const auto *osh_731 = buffer.data(osh + 731);
    const auto *osh_732 = buffer.data(osh + 732);
    const auto *osh_733 = buffer.data(osh + 733);
    const auto *osh_734 = buffer.data(osh + 734);
    const auto *osh_735 = buffer.data(osh + 735);
    const auto *osh_740 = buffer.data(osh + 740);
    const auto *osh_744 = buffer.data(osh + 744);
    const auto *osh_749 = buffer.data(osh + 749);
    const auto *osh_750 = buffer.data(osh + 750);
    const auto *osh_751 = buffer.data(osh + 751);
    const auto *osh_752 = buffer.data(osh + 752);
    const auto *osh_753 = buffer.data(osh + 753);
    const auto *osh_755 = buffer.data(osh + 755);
    const auto *osh_756 = buffer.data(osh + 756);
    const auto *osh_759 = buffer.data(osh + 759);
    const auto *osh_762 = buffer.data(osh + 762);
    const auto *osh_766 = buffer.data(osh + 766);
    const auto *osh_771 = buffer.data(osh + 771);
    const auto *osh_773 = buffer.data(osh + 773);
    const auto *osh_774 = buffer.data(osh + 774);
    const auto *osh_775 = buffer.data(osh + 775);
    const auto *osh_776 = buffer.data(osh + 776);
    const auto *osh_782 = buffer.data(osh + 782);

    const auto *osi1_756 = buffer.data(osi1 + 756);
    const auto *osi1_759 = buffer.data(osi1 + 759);
    const auto *osi1_761 = buffer.data(osi1 + 761);
    const auto *osi1_762 = buffer.data(osi1 + 762);
    const auto *osi1_765 = buffer.data(osi1 + 765);
    const auto *osi1_766 = buffer.data(osi1 + 766);
    const auto *osi1_768 = buffer.data(osi1 + 768);
    const auto *osi1_770 = buffer.data(osi1 + 770);
    const auto *osi1_783 = buffer.data(osi1 + 783);
    const auto *osi1_784 = buffer.data(osi1 + 784);
    const auto *osi1_787 = buffer.data(osi1 + 787);
    const auto *osi1_790 = buffer.data(osi1 + 790);

    const auto *qsg0_498 = buffer.data(qsg0 + 498);
    const auto *qsg0_500 = buffer.data(qsg0 + 500);
    const auto *qsg0_501 = buffer.data(qsg0 + 501);
    const auto *qsg0_504 = buffer.data(qsg0 + 504);
    const auto *qsg0_505 = buffer.data(qsg0 + 505);
    const auto *qsg0_507 = buffer.data(qsg0 + 507);
    const auto *qsg0_508 = buffer.data(qsg0 + 508);
    const auto *qsg0_509 = buffer.data(qsg0 + 509);
    const auto *qsg0_520 = buffer.data(qsg0 + 520);
    const auto *qsg0_522 = buffer.data(qsg0 + 522);
    const auto *qsg0_523 = buffer.data(qsg0 + 523);
    const auto *qsg0_524 = buffer.data(qsg0 + 524);
    const auto *qsg0_525 = buffer.data(qsg0 + 525);
    const auto *qsg0_526 = buffer.data(qsg0 + 526);
    const auto *qsg0_527 = buffer.data(qsg0 + 527);
    const auto *qsg0_528 = buffer.data(qsg0 + 528);
    const auto *qsg0_529 = buffer.data(qsg0 + 529);
    const auto *qsg0_530 = buffer.data(qsg0 + 530);
    const auto *qsg0_534 = buffer.data(qsg0 + 534);
    const auto *qsg0_535 = buffer.data(qsg0 + 535);
    const auto *qsg0_536 = buffer.data(qsg0 + 536);
    const auto *qsg0_537 = buffer.data(qsg0 + 537);
    const auto *qsg0_538 = buffer.data(qsg0 + 538);
    const auto *qsg0_539 = buffer.data(qsg0 + 539);
    const auto *qsg0_540 = buffer.data(qsg0 + 540);
    const auto *qsg0_542 = buffer.data(qsg0 + 542);
    const auto *qsg0_543 = buffer.data(qsg0 + 543);
    const auto *qsg0_545 = buffer.data(qsg0 + 545);
    const auto *qsg0_546 = buffer.data(qsg0 + 546);
    const auto *qsg0_550 = buffer.data(qsg0 + 550);
    const auto *qsg0_551 = buffer.data(qsg0 + 551);
    const auto *qsg0_552 = buffer.data(qsg0 + 552);
    const auto *qsg0_554 = buffer.data(qsg0 + 554);
    const auto *qsg0_560 = buffer.data(qsg0 + 560);

    const auto *qsg1_498 = buffer.data(qsg1 + 498);
    const auto *qsg1_500 = buffer.data(qsg1 + 500);
    const auto *qsg1_501 = buffer.data(qsg1 + 501);
    const auto *qsg1_504 = buffer.data(qsg1 + 504);
    const auto *qsg1_505 = buffer.data(qsg1 + 505);
    const auto *qsg1_507 = buffer.data(qsg1 + 507);
    const auto *qsg1_508 = buffer.data(qsg1 + 508);
    const auto *qsg1_509 = buffer.data(qsg1 + 509);
    const auto *qsg1_520 = buffer.data(qsg1 + 520);
    const auto *qsg1_522 = buffer.data(qsg1 + 522);
    const auto *qsg1_523 = buffer.data(qsg1 + 523);
    const auto *qsg1_524 = buffer.data(qsg1 + 524);
    const auto *qsg1_525 = buffer.data(qsg1 + 525);
    const auto *qsg1_526 = buffer.data(qsg1 + 526);
    const auto *qsg1_527 = buffer.data(qsg1 + 527);
    const auto *qsg1_528 = buffer.data(qsg1 + 528);
    const auto *qsg1_529 = buffer.data(qsg1 + 529);
    const auto *qsg1_530 = buffer.data(qsg1 + 530);
    const auto *qsg1_534 = buffer.data(qsg1 + 534);
    const auto *qsg1_535 = buffer.data(qsg1 + 535);
    const auto *qsg1_536 = buffer.data(qsg1 + 536);
    const auto *qsg1_537 = buffer.data(qsg1 + 537);
    const auto *qsg1_538 = buffer.data(qsg1 + 538);
    const auto *qsg1_539 = buffer.data(qsg1 + 539);
    const auto *qsg1_540 = buffer.data(qsg1 + 540);
    const auto *qsg1_542 = buffer.data(qsg1 + 542);
    const auto *qsg1_543 = buffer.data(qsg1 + 543);
    const auto *qsg1_545 = buffer.data(qsg1 + 545);
    const auto *qsg1_546 = buffer.data(qsg1 + 546);
    const auto *qsg1_550 = buffer.data(qsg1 + 550);
    const auto *qsg1_551 = buffer.data(qsg1 + 551);
    const auto *qsg1_552 = buffer.data(qsg1 + 552);
    const auto *qsg1_554 = buffer.data(qsg1 + 554);
    const auto *qsg1_560 = buffer.data(qsg1 + 560);

    const auto *qsh_693 = buffer.data(qsh + 693);
    const auto *qsh_695 = buffer.data(qsh + 695);
    const auto *qsh_696 = buffer.data(qsh + 696);
    const auto *qsh_698 = buffer.data(qsh + 698);
    const auto *qsh_699 = buffer.data(qsh + 699);
    const auto *qsh_702 = buffer.data(qsh + 702);
    const auto *qsh_703 = buffer.data(qsh + 703);
    const auto *qsh_705 = buffer.data(qsh + 705);
    const auto *qsh_707 = buffer.data(qsh + 707);
    const auto *qsh_708 = buffer.data(qsh + 708);
    const auto *qsh_709 = buffer.data(qsh + 709);
    const auto *qsh_710 = buffer.data(qsh + 710);
    const auto *qsh_711 = buffer.data(qsh + 711);
    const auto *qsh_712 = buffer.data(qsh + 712);
    const auto *qsh_713 = buffer.data(qsh + 713);
    const auto *qsh_714 = buffer.data(qsh + 714);
    const auto *qsh_716 = buffer.data(qsh + 716);
    const auto *qsh_717 = buffer.data(qsh + 717);
    const auto *qsh_719 = buffer.data(qsh + 719);
    const auto *qsh_720 = buffer.data(qsh + 720);
    const auto *qsh_723 = buffer.data(qsh + 723);
    const auto *qsh_729 = buffer.data(qsh + 729);
    const auto *qsh_730 = buffer.data(qsh + 730);
    const auto *qsh_731 = buffer.data(qsh + 731);
    const auto *qsh_732 = buffer.data(qsh + 732);
    const auto *qsh_733 = buffer.data(qsh + 733);
    const auto *qsh_734 = buffer.data(qsh + 734);
    const auto *qsh_735 = buffer.data(qsh + 735);
    const auto *qsh_736 = buffer.data(qsh + 736);
    const auto *qsh_737 = buffer.data(qsh + 737);
    const auto *qsh_738 = buffer.data(qsh + 738);
    const auto *qsh_739 = buffer.data(qsh + 739);
    const auto *qsh_740 = buffer.data(qsh + 740);
    const auto *qsh_741 = buffer.data(qsh + 741);
    const auto *qsh_742 = buffer.data(qsh + 742);
    const auto *qsh_743 = buffer.data(qsh + 743);
    const auto *qsh_744 = buffer.data(qsh + 744);
    const auto *qsh_749 = buffer.data(qsh + 749);
    const auto *qsh_750 = buffer.data(qsh + 750);
    const auto *qsh_751 = buffer.data(qsh + 751);
    const auto *qsh_752 = buffer.data(qsh + 752);
    const auto *qsh_753 = buffer.data(qsh + 753);
    const auto *qsh_754 = buffer.data(qsh + 754);
    const auto *qsh_755 = buffer.data(qsh + 755);
    const auto *qsh_756 = buffer.data(qsh + 756);
    const auto *qsh_757 = buffer.data(qsh + 757);
    const auto *qsh_758 = buffer.data(qsh + 758);
    const auto *qsh_759 = buffer.data(qsh + 759);
    const auto *qsh_761 = buffer.data(qsh + 761);
    const auto *qsh_762 = buffer.data(qsh + 762);
    const auto *qsh_763 = buffer.data(qsh + 763);
    const auto *qsh_765 = buffer.data(qsh + 765);
    const auto *qsh_766 = buffer.data(qsh + 766);
    const auto *qsh_771 = buffer.data(qsh + 771);
    const auto *qsh_772 = buffer.data(qsh + 772);
    const auto *qsh_773 = buffer.data(qsh + 773);
    const auto *qsh_774 = buffer.data(qsh + 774);
    const auto *qsh_775 = buffer.data(qsh + 775);
    const auto *qsh_776 = buffer.data(qsh + 776);
    const auto *qsh_777 = buffer.data(qsh + 777);
    const auto *qsh_779 = buffer.data(qsh + 779);
    const auto *qsh_780 = buffer.data(qsh + 780);
    const auto *qsh_782 = buffer.data(qsh + 782);

#pragma omp simd aligned(t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, osh_525, osh_546, \
                         osh_548, osh_696, qsg0_498, qsg1_498, qsh_693, qsh_695, \
                         qsh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_12 * osh_546[k]
                   + f_3 * pc_y[k] * qsh_693[k];

        t_926[k] = f_22 * osh_525[k]
                   + f_3 * pc_z[k] * qsh_693[k];

        t_927[k] = f_22 * osh_696[k]
                   + f_8 * qsg0_498[k]
                   - f_9 * qsg1_498[k]
                   + f_3 * pc_x[k] * qsh_696[k];

        t_928[k] = f_12 * osh_548[k]
                   + f_3 * pc_y[k] * qsh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pc_x, pc_z, osh_528, osh_698, osh_699, qsg0_500, \
                         qsg0_501, qsg1_500, qsg1_501, qsh_696, qsh_698, \
                         qsh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_22 * osh_698[k]
                   + f_8 * qsg0_500[k]
                   - f_9 * qsg1_500[k]
                   + f_3 * pc_x[k] * qsh_698[k];

        t_930[k] = f_22 * osh_699[k]
                   + f_6 * qsg0_501[k]
                   - f_7 * qsg1_501[k]
                   + f_3 * pc_x[k] * qsh_699[k];

        t_931[k] = f_22 * osh_528[k]
                   + f_3 * pc_z[k] * qsh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pc_x, pc_y, osh_551, osh_702, osh_703, qsg0_504, \
                         qsg0_505, qsg1_504, qsg1_505, qsh_698, qsh_702, \
                         qsh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * osh_551[k]
                   + f_3 * pc_y[k] * qsh_698[k];

        t_933[k] = f_22 * osh_702[k]
                   + f_6 * qsg0_504[k]
                   - f_7 * qsg1_504[k]
                   + f_3 * pc_x[k] * qsh_702[k];

        t_934[k] = f_22 * osh_703[k]
                   + f_4 * qsg0_505[k]
                   - f_5 * qsg1_505[k]
                   + f_3 * pc_x[k] * qsh_703[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pc_x, pc_y, pc_z, osh_531, osh_555, osh_705, \
                         qsg0_507, qsg1_507, qsh_699, qsh_702, \
                         qsh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_22 * osh_531[k]
                   + f_3 * pc_z[k] * qsh_699[k];

        t_936[k] = f_22 * osh_705[k]
                   + f_4 * qsg0_507[k]
                   - f_5 * qsg1_507[k]
                   + f_3 * pc_x[k] * qsh_705[k];

        t_937[k] = f_12 * osh_555[k]
                   + f_3 * pc_y[k] * qsh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pc_x, osh_707, osh_708, osh_709, osh_710, \
                         qsg0_509, qsg1_509, qsh_707, qsh_708, qsh_709, \
                         qsh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_22 * osh_707[k]
                   + f_4 * qsg0_509[k]
                   - f_5 * qsg1_509[k]
                   + f_3 * pc_x[k] * qsh_707[k];

        t_939[k] = f_22 * osh_708[k]
                   + f_3 * pc_x[k] * qsh_708[k];

        t_940[k] = f_22 * osh_709[k]
                   + f_3 * pc_x[k] * qsh_709[k];

        t_941[k] = f_22 * osh_710[k]
                   + f_3 * pc_x[k] * qsh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, osh_561, osh_711, osh_712, \
                         osh_713, qsg0_505, qsg1_505, qsh_708, qsh_711, qsh_712, \
                         qsh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_22 * osh_711[k]
                   + f_3 * pc_x[k] * qsh_711[k];

        t_943[k] = f_22 * osh_712[k]
                   + f_3 * pc_x[k] * qsh_712[k];

        t_944[k] = f_22 * osh_713[k]
                   + f_3 * pc_x[k] * qsh_713[k];

        t_945[k] = f_12 * osh_561[k]
                   + f_1 * qsg0_505[k]
                   - f_2 * qsg1_505[k]
                   + f_3 * pc_y[k] * qsh_708[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pc_y, pc_z, osh_540, osh_563, osh_564, qsg0_507, \
                         qsg0_508, qsg1_507, qsg1_508, qsh_708, qsh_710, \
                         qsh_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_22 * osh_540[k]
                   + f_3 * pc_z[k] * qsh_708[k];

        t_947[k] = f_12 * osh_563[k]
                   + f_8 * qsg0_507[k]
                   - f_9 * qsg1_507[k]
                   + f_3 * pc_y[k] * qsh_710[k];

        t_948[k] = f_12 * osh_564[k]
                   + f_6 * qsg0_508[k]
                   - f_7 * qsg1_508[k]
                   + f_3 * pc_y[k] * qsh_711[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_y, pc_y, pc_z, osi0_756, osh_545, \
                         osh_565, osh_566, osi1_756, qsg0_509, qsg1_509, qsh_712, \
                         qsh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_12 * osh_565[k]
                   + f_4 * qsg0_509[k]
                   - f_5 * qsg1_509[k]
                   + f_3 * pc_y[k] * qsh_712[k];

        t_950[k] = f_12 * osh_566[k]
                   + f_3 * pc_y[k] * qsh_713[k];

        t_951[k] = f_22 * osh_545[k]
                   + f_1 * qsg0_509[k]
                   - f_2 * qsg1_509[k]
                   + f_3 * pc_z[k] * qsh_713[k];

        t_952[k] = pa_y[k] * osi0_756[k]
                   - f_10 * pc_y[k] * osi1_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pa_y, pc_y, pc_z, osi0_759, osh_546, \
                         osh_567, osh_568, osh_569, osi1_759, qsh_714, \
                         qsh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_11 * osh_567[k]
                   + f_3 * pc_y[k] * qsh_714[k];

        t_954[k] = f_23 * osh_546[k]
                   + f_3 * pc_z[k] * qsh_714[k];

        t_955[k] = pa_y[k] * osi0_759[k]
                   + f_12 * osh_568[k]
                   - f_10 * pc_y[k] * osi1_759[k];

        t_956[k] = f_11 * osh_569[k]
                   + f_3 * pc_y[k] * qsh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pa_y, pc_y, pc_z, osi0_761, osi0_762, \
                         osh_549, osh_570, osh_572, osi1_761, osi1_762, qsh_717, \
                         qsh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pa_y[k] * osi0_761[k]
                   - f_10 * pc_y[k] * osi1_761[k];

        t_958[k] = pa_y[k] * osi0_762[k]
                   + f_13 * osh_570[k]
                   - f_10 * pc_y[k] * osi1_762[k];

        t_959[k] = f_23 * osh_549[k]
                   + f_3 * pc_z[k] * qsh_717[k];

        t_960[k] = f_11 * osh_572[k]
                   + f_3 * pc_y[k] * qsh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pa_y, pc_y, pc_z, osi0_765, osi0_766, osh_552, \
                         osh_573, osi1_765, osi1_766, qsh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pa_y[k] * osi0_765[k]
                   - f_10 * pc_y[k] * osi1_765[k];

        t_962[k] = pa_y[k] * osi0_766[k]
                   + f_14 * osh_573[k]
                   - f_10 * pc_y[k] * osi1_766[k];

        t_963[k] = f_23 * osh_552[k]
                   + f_3 * pc_z[k] * qsh_720[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pa_y, pc_x, pc_y, osi0_768, osi0_770, \
                         osh_575, osh_576, osh_729, osi1_768, osi1_770, qsh_723, \
                         qsh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pa_y[k] * osi0_768[k]
                   + f_12 * osh_575[k]
                   - f_10 * pc_y[k] * osi1_768[k];

        t_965[k] = f_11 * osh_576[k]
                   + f_3 * pc_y[k] * qsh_723[k];

        t_966[k] = pa_y[k] * osi0_770[k]
                   - f_10 * pc_y[k] * osi1_770[k];

        t_967[k] = f_22 * osh_729[k]
                   + f_3 * pc_x[k] * qsh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, osh_730, osh_731, osh_732, \
                         osh_733, osh_734, qsh_730, qsh_731, qsh_732, qsh_733, \
                         qsh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_22 * osh_730[k]
                   + f_3 * pc_x[k] * qsh_730[k];

        t_969[k] = f_22 * osh_731[k]
                   + f_3 * pc_x[k] * qsh_731[k];

        t_970[k] = f_22 * osh_732[k]
                   + f_3 * pc_x[k] * qsh_732[k];

        t_971[k] = f_22 * osh_733[k]
                   + f_3 * pc_x[k] * qsh_733[k];

        t_972[k] = f_22 * osh_734[k]
                   + f_3 * pc_x[k] * qsh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pc_y, pc_z, osh_561, osh_582, osh_584, qsg0_520, \
                         qsg0_522, qsg1_520, qsg1_522, qsh_729, \
                         qsh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_11 * osh_582[k]
                   + f_1 * qsg0_520[k]
                   - f_2 * qsg1_520[k]
                   + f_3 * pc_y[k] * qsh_729[k];

        t_974[k] = f_23 * osh_561[k]
                   + f_3 * pc_z[k] * qsh_729[k];

        t_975[k] = f_11 * osh_584[k]
                   + f_8 * qsg0_522[k]
                   - f_9 * qsg1_522[k]
                   + f_3 * pc_y[k] * qsh_731[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_y, osh_585, osh_586, osh_587, qsg0_523, \
                         qsg0_524, qsg1_523, qsg1_524, qsh_732, qsh_733, \
                         qsh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_11 * osh_585[k]
                   + f_6 * qsg0_523[k]
                   - f_7 * qsg1_523[k]
                   + f_3 * pc_y[k] * qsh_732[k];

        t_977[k] = f_11 * osh_586[k]
                   + f_4 * qsg0_524[k]
                   - f_5 * qsg1_524[k]
                   + f_3 * pc_y[k] * qsh_733[k];

        t_978[k] = f_11 * osh_587[k]
                   + f_3 * pc_y[k] * qsh_734[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pa_y, pc_x, pc_y, pc_z, osi0_783, \
                         osh_567, osh_735, osi1_783, qsg0_525, qsg1_525, \
                         qsh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pa_y[k] * osi0_783[k]
                   - f_10 * pc_y[k] * osi1_783[k];

        t_980[k] = f_22 * osh_735[k]
                   + f_1 * qsg0_525[k]
                   - f_2 * qsg1_525[k]
                   + f_3 * pc_x[k] * qsh_735[k];

        t_981[k] = f_3 * pc_y[k] * qsh_735[k];

        t_982[k] = f_21 * osh_567[k]
                   + f_3 * pc_z[k] * qsh_735[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_x, pc_y, osh_740, qsg0_525, qsg0_530, \
                         qsg1_525, qsg1_530, qsh_736, qsh_737, \
                         qsh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_4 * qsg0_525[k]
                   - f_5 * qsg1_525[k]
                   + f_3 * pc_y[k] * qsh_736[k];

        t_984[k] = f_3 * pc_y[k] * qsh_737[k];

        t_985[k] = f_22 * osh_740[k]
                   + f_8 * qsg0_530[k]
                   - f_9 * qsg1_530[k]
                   + f_3 * pc_x[k] * qsh_740[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_y, qsg0_526, qsg0_527, qsg1_526, qsg1_527, \
                         qsh_738, qsh_739, qsh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_6 * qsg0_526[k]
                   - f_7 * qsg1_526[k]
                   + f_3 * pc_y[k] * qsh_738[k];

        t_987[k] = f_4 * qsg0_527[k]
                   - f_5 * qsg1_527[k]
                   + f_3 * pc_y[k] * qsh_739[k];

        t_988[k] = f_3 * pc_y[k] * qsh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_y, osh_744, qsg0_528, qsg0_529, \
                         qsg0_534, qsg1_528, qsg1_529, qsg1_534, qsh_741, qsh_742, \
                         qsh_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_22 * osh_744[k]
                   + f_6 * qsg0_534[k]
                   - f_7 * qsg1_534[k]
                   + f_3 * pc_x[k] * qsh_744[k];

        t_990[k] = f_8 * qsg0_528[k]
                   - f_9 * qsg1_528[k]
                   + f_3 * pc_y[k] * qsh_741[k];

        t_991[k] = f_6 * qsg0_529[k]
                   - f_7 * qsg1_529[k]
                   + f_3 * pc_y[k] * qsh_742[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, pc_x, pc_y, osh_749, osh_750, qsg0_530, \
                         qsg0_539, qsg1_530, qsg1_539, qsh_743, qsh_744, qsh_749, \
                         qsh_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_4 * qsg0_530[k]
                   - f_5 * qsg1_530[k]
                   + f_3 * pc_y[k] * qsh_743[k];

        t_993[k] = f_3 * pc_y[k] * qsh_744[k];

        t_994[k] = f_22 * osh_749[k]
                   + f_4 * qsg0_539[k]
                   - f_5 * qsg1_539[k]
                   + f_3 * pc_x[k] * qsh_749[k];

        t_995[k] = f_22 * osh_750[k]
                   + f_3 * pc_x[k] * qsh_750[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, pc_x, pc_y, osh_751, osh_752, \
                         osh_753, osh_755, qsh_749, qsh_751, qsh_752, qsh_753, \
                         qsh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_22 * osh_751[k]
                   + f_3 * pc_x[k] * qsh_751[k];

        t_997[k] = f_22 * osh_752[k]
                   + f_3 * pc_x[k] * qsh_752[k];

        t_998[k] = f_22 * osh_753[k]
                   + f_3 * pc_x[k] * qsh_753[k];

        t_999[k] = f_3 * pc_y[k] * qsh_749[k];

        t_1000[k] = f_22 * osh_755[k]
                    + f_3 * pc_x[k] * qsh_755[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, pc_y, qsg0_535, qsg0_536, qsg0_537, qsg1_535, \
                         qsg1_536, qsg1_537, qsh_750, qsh_751, \
                         qsh_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_1 * qsg0_535[k]
                    - f_2 * qsg1_535[k]
                    + f_3 * pc_y[k] * qsh_750[k];

        t_1002[k] = f_16 * qsg0_536[k]
                    - f_17 * qsg1_536[k]
                    + f_3 * pc_y[k] * qsh_751[k];

        t_1003[k] = f_8 * qsg0_537[k]
                    - f_9 * qsg1_537[k]
                    + f_3 * pc_y[k] * qsh_752[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, osh_587, qsg0_538, \
                         qsg0_539, qsg1_538, qsg1_539, qsh_753, qsh_754, \
                         qsh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_6 * qsg0_538[k]
                    - f_7 * qsg1_538[k]
                    + f_3 * pc_y[k] * qsh_753[k];

        t_1005[k] = f_4 * qsg0_539[k]
                    - f_5 * qsg1_539[k]
                    + f_3 * pc_y[k] * qsh_754[k];

        t_1006[k] = f_3 * pc_y[k] * qsh_755[k];

        t_1007[k] = f_21 * osh_587[k]
                    + f_1 * qsg0_539[k]
                    - f_2 * qsg1_539[k]
                    + f_3 * pc_z[k] * qsh_755[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, osh_588, osh_756, \
                         osh_759, qsg0_540, qsg0_543, qsg1_540, qsg1_543, qsh_756, \
                         qsh_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_14 * osh_756[k]
                    + f_1 * qsg0_540[k]
                    - f_2 * qsg1_540[k]
                    + f_3 * pc_x[k] * qsh_756[k];

        t_1009[k] = f_20 * osh_588[k]
                    + f_3 * pc_y[k] * qsh_756[k];

        t_1010[k] = f_3 * pc_z[k] * qsh_756[k];

        t_1011[k] = f_14 * osh_759[k]
                    + f_8 * qsg0_543[k]
                    - f_9 * qsg1_543[k]
                    + f_3 * pc_x[k] * qsh_759[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, t_1015, pc_x, pc_z, osh_762, qsg0_540, \
                         qsg0_546, qsg1_540, qsg1_546, qsh_757, qsh_758, qsh_759, \
                         qsh_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_3 * pc_z[k] * qsh_757[k];

        t_1013[k] = f_4 * qsg0_540[k]
                    - f_5 * qsg1_540[k]
                    + f_3 * pc_z[k] * qsh_758[k];

        t_1014[k] = f_14 * osh_762[k]
                    + f_6 * qsg0_546[k]
                    - f_7 * qsg1_546[k]
                    + f_3 * pc_x[k] * qsh_762[k];

        t_1015[k] = f_3 * pc_z[k] * qsh_759[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, pc_z, osh_593, osh_766, \
                         qsg0_542, qsg0_550, qsg1_542, qsg1_550, qsh_761, qsh_762, \
                         qsh_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_20 * osh_593[k]
                    + f_3 * pc_y[k] * qsh_761[k];

        t_1017[k] = f_6 * qsg0_542[k]
                    - f_7 * qsg1_542[k]
                    + f_3 * pc_z[k] * qsh_761[k];

        t_1018[k] = f_14 * osh_766[k]
                    + f_4 * qsg0_550[k]
                    - f_5 * qsg1_550[k]
                    + f_3 * pc_x[k] * qsh_766[k];

        t_1019[k] = f_3 * pc_z[k] * qsh_762[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, pc_x, pc_y, pc_z, osh_597, osh_771, \
                         qsg0_543, qsg0_545, qsg1_543, qsg1_545, qsh_763, qsh_765, \
                         qsh_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_4 * qsg0_543[k]
                    - f_5 * qsg1_543[k]
                    + f_3 * pc_z[k] * qsh_763[k];

        t_1021[k] = f_20 * osh_597[k]
                    + f_3 * pc_y[k] * qsh_765[k];

        t_1022[k] = f_8 * qsg0_545[k]
                    - f_9 * qsg1_545[k]
                    + f_3 * pc_z[k] * qsh_765[k];

        t_1023[k] = f_14 * osh_771[k]
                    + f_3 * pc_x[k] * qsh_771[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, t_1027, t_1028, pc_x, pc_z, osh_773, osh_774, \
                         osh_775, osh_776, qsh_766, qsh_773, qsh_774, qsh_775, \
                         qsh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_3 * pc_z[k] * qsh_766[k];

        t_1025[k] = f_14 * osh_773[k]
                    + f_3 * pc_x[k] * qsh_773[k];

        t_1026[k] = f_14 * osh_774[k]
                    + f_3 * pc_x[k] * qsh_774[k];

        t_1027[k] = f_14 * osh_775[k]
                    + f_3 * pc_x[k] * qsh_775[k];

        t_1028[k] = f_14 * osh_776[k]
                    + f_3 * pc_x[k] * qsh_776[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pc_y, pc_z, osh_603, qsg0_550, \
                         qsg0_551, qsg1_550, qsg1_551, qsh_771, qsh_772, \
                         qsh_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = f_20 * osh_603[k]
                    + f_1 * qsg0_550[k]
                    - f_2 * qsg1_550[k]
                    + f_3 * pc_y[k] * qsh_771[k];

        t_1030[k] = f_3 * pc_z[k] * qsh_771[k];

        t_1031[k] = f_4 * qsg0_550[k]
                    - f_5 * qsg1_550[k]
                    + f_3 * pc_z[k] * qsh_772[k];

        t_1032[k] = f_6 * qsg0_551[k]
                    - f_7 * qsg1_551[k]
                    + f_3 * pc_z[k] * qsh_773[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, t_1036, pa_z, pc_y, pc_z, osi0_784, osh_608, \
                         osi1_784, qsg0_552, qsg0_554, qsg1_552, qsg1_554, qsh_774, \
                         qsh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_8 * qsg0_552[k]
                    - f_9 * qsg1_552[k]
                    + f_3 * pc_z[k] * qsh_774[k];

        t_1034[k] = f_20 * osh_608[k]
                    + f_3 * pc_y[k] * qsh_776[k];

        t_1035[k] = f_1 * qsg0_554[k]
                    - f_2 * qsg1_554[k]
                    + f_3 * pc_z[k] * qsh_776[k];

        t_1036[k] = pa_z[k] * osi0_784[k]
                    - f_10 * pc_z[k] * osi1_784[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, pa_z, pc_y, pc_z, osi0_787, osh_588, \
                         osh_609, osh_611, osi1_787, qsh_777, qsh_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_21 * osh_609[k]
                    + f_3 * pc_y[k] * qsh_777[k];

        t_1038[k] = f_11 * osh_588[k]
                    + f_3 * pc_z[k] * qsh_777[k];

        t_1039[k] = pa_z[k] * osi0_787[k]
                    - f_10 * pc_z[k] * osi1_787[k];

        t_1040[k] = f_21 * osh_611[k]
                    + f_3 * pc_y[k] * qsh_779[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, pa_z, pc_x, pc_z, osi0_790, osh_591, osh_782, \
                         osi1_790, qsg0_560, qsg1_560, qsh_780, \
                         qsh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_14 * osh_782[k]
                    + f_8 * qsg0_560[k]
                    - f_9 * qsg1_560[k]
                    + f_3 * pc_x[k] * qsh_782[k];

        t_1042[k] = pa_z[k] * osi0_790[k]
                    - f_10 * pc_z[k] * osi1_790[k];

        t_1043[k] = f_11 * osh_591[k]
                    + f_3 * pc_z[k] * qsh_780[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t osi0,
                                                          const size_t osh, const size_t osi1,
                                                          const size_t qsg0, const size_t qsg1,
                                                          const size_t qsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_794 = buffer.data(osi0 + 794);
    const auto *osi0_796 = buffer.data(osi0 + 796);
    const auto *osi0_805 = buffer.data(osi0 + 805);

    const auto *osh_594 = buffer.data(osh + 594);
    const auto *osh_595 = buffer.data(osh + 595);
    const auto *osh_603 = buffer.data(osh + 603);
    const auto *osh_608 = buffer.data(osh + 608);
    const auto *osh_609 = buffer.data(osh + 609);
    const auto *osh_612 = buffer.data(osh + 612);
    const auto *osh_614 = buffer.data(osh + 614);
    const auto *osh_615 = buffer.data(osh + 615);
    const auto *osh_618 = buffer.data(osh + 618);
    const auto *osh_624 = buffer.data(osh + 624);
    const auto *osh_626 = buffer.data(osh + 626);
    const auto *osh_627 = buffer.data(osh + 627);
    const auto *osh_628 = buffer.data(osh + 628);
    const auto *osh_629 = buffer.data(osh + 629);
    const auto *osh_630 = buffer.data(osh + 630);
    const auto *osh_632 = buffer.data(osh + 632);
    const auto *osh_633 = buffer.data(osh + 633);
    const auto *osh_635 = buffer.data(osh + 635);
    const auto *osh_636 = buffer.data(osh + 636);
    const auto *osh_639 = buffer.data(osh + 639);
    const auto *osh_645 = buffer.data(osh + 645);
    const auto *osh_647 = buffer.data(osh + 647);
    const auto *osh_648 = buffer.data(osh + 648);
    const auto *osh_649 = buffer.data(osh + 649);
    const auto *osh_650 = buffer.data(osh + 650);
    const auto *osh_651 = buffer.data(osh + 651);
    const auto *osh_653 = buffer.data(osh + 653);
    const auto *osh_654 = buffer.data(osh + 654);
    const auto *osh_656 = buffer.data(osh + 656);
    const auto *osh_657 = buffer.data(osh + 657);
    const auto *osh_660 = buffer.data(osh + 660);
    const auto *osh_666 = buffer.data(osh + 666);
    const auto *osh_668 = buffer.data(osh + 668);
    const auto *osh_669 = buffer.data(osh + 669);
    const auto *osh_670 = buffer.data(osh + 670);
    const auto *osh_671 = buffer.data(osh + 671);
    const auto *osh_672 = buffer.data(osh + 672);
    const auto *osh_674 = buffer.data(osh + 674);
    const auto *osh_677 = buffer.data(osh + 677);
    const auto *osh_681 = buffer.data(osh + 681);
    const auto *osh_687 = buffer.data(osh + 687);
    const auto *osh_689 = buffer.data(osh + 689);
    const auto *osh_690 = buffer.data(osh + 690);
    const auto *osh_691 = buffer.data(osh + 691);
    const auto *osh_692 = buffer.data(osh + 692);
    const auto *osh_786 = buffer.data(osh + 786);
    const auto *osh_791 = buffer.data(osh + 791);
    const auto *osh_792 = buffer.data(osh + 792);
    const auto *osh_793 = buffer.data(osh + 793);
    const auto *osh_794 = buffer.data(osh + 794);
    const auto *osh_795 = buffer.data(osh + 795);
    const auto *osh_796 = buffer.data(osh + 796);
    const auto *osh_797 = buffer.data(osh + 797);
    const auto *osh_798 = buffer.data(osh + 798);
    const auto *osh_801 = buffer.data(osh + 801);
    const auto *osh_803 = buffer.data(osh + 803);
    const auto *osh_804 = buffer.data(osh + 804);
    const auto *osh_807 = buffer.data(osh + 807);
    const auto *osh_808 = buffer.data(osh + 808);
    const auto *osh_810 = buffer.data(osh + 810);
    const auto *osh_812 = buffer.data(osh + 812);
    const auto *osh_813 = buffer.data(osh + 813);
    const auto *osh_814 = buffer.data(osh + 814);
    const auto *osh_815 = buffer.data(osh + 815);
    const auto *osh_816 = buffer.data(osh + 816);
    const auto *osh_817 = buffer.data(osh + 817);
    const auto *osh_818 = buffer.data(osh + 818);
    const auto *osh_819 = buffer.data(osh + 819);
    const auto *osh_822 = buffer.data(osh + 822);
    const auto *osh_824 = buffer.data(osh + 824);
    const auto *osh_825 = buffer.data(osh + 825);
    const auto *osh_828 = buffer.data(osh + 828);
    const auto *osh_829 = buffer.data(osh + 829);
    const auto *osh_831 = buffer.data(osh + 831);
    const auto *osh_833 = buffer.data(osh + 833);
    const auto *osh_834 = buffer.data(osh + 834);
    const auto *osh_835 = buffer.data(osh + 835);
    const auto *osh_836 = buffer.data(osh + 836);
    const auto *osh_837 = buffer.data(osh + 837);
    const auto *osh_838 = buffer.data(osh + 838);
    const auto *osh_839 = buffer.data(osh + 839);
    const auto *osh_840 = buffer.data(osh + 840);
    const auto *osh_843 = buffer.data(osh + 843);
    const auto *osh_845 = buffer.data(osh + 845);
    const auto *osh_846 = buffer.data(osh + 846);
    const auto *osh_849 = buffer.data(osh + 849);
    const auto *osh_850 = buffer.data(osh + 850);
    const auto *osh_852 = buffer.data(osh + 852);
    const auto *osh_854 = buffer.data(osh + 854);
    const auto *osh_855 = buffer.data(osh + 855);
    const auto *osh_856 = buffer.data(osh + 856);
    const auto *osh_857 = buffer.data(osh + 857);
    const auto *osh_858 = buffer.data(osh + 858);
    const auto *osh_859 = buffer.data(osh + 859);
    const auto *osh_860 = buffer.data(osh + 860);
    const auto *osh_861 = buffer.data(osh + 861);

    const auto *osi1_794 = buffer.data(osi1 + 794);
    const auto *osi1_796 = buffer.data(osi1 + 796);
    const auto *osi1_805 = buffer.data(osi1 + 805);

    const auto *qsg0_564 = buffer.data(qsg0 + 564);
    const auto *qsg0_567 = buffer.data(qsg0 + 567);
    const auto *qsg0_568 = buffer.data(qsg0 + 568);
    const auto *qsg0_569 = buffer.data(qsg0 + 569);
    const auto *qsg0_570 = buffer.data(qsg0 + 570);
    const auto *qsg0_573 = buffer.data(qsg0 + 573);
    const auto *qsg0_575 = buffer.data(qsg0 + 575);
    const auto *qsg0_576 = buffer.data(qsg0 + 576);
    const auto *qsg0_579 = buffer.data(qsg0 + 579);
    const auto *qsg0_580 = buffer.data(qsg0 + 580);
    const auto *qsg0_582 = buffer.data(qsg0 + 582);
    const auto *qsg0_583 = buffer.data(qsg0 + 583);
    const auto *qsg0_584 = buffer.data(qsg0 + 584);
    const auto *qsg0_585 = buffer.data(qsg0 + 585);
    const auto *qsg0_588 = buffer.data(qsg0 + 588);
    const auto *qsg0_590 = buffer.data(qsg0 + 590);
    const auto *qsg0_591 = buffer.data(qsg0 + 591);
    const auto *qsg0_594 = buffer.data(qsg0 + 594);
    const auto *qsg0_595 = buffer.data(qsg0 + 595);
    const auto *qsg0_597 = buffer.data(qsg0 + 597);
    const auto *qsg0_598 = buffer.data(qsg0 + 598);
    const auto *qsg0_599 = buffer.data(qsg0 + 599);
    const auto *qsg0_600 = buffer.data(qsg0 + 600);
    const auto *qsg0_603 = buffer.data(qsg0 + 603);
    const auto *qsg0_605 = buffer.data(qsg0 + 605);
    const auto *qsg0_606 = buffer.data(qsg0 + 606);
    const auto *qsg0_609 = buffer.data(qsg0 + 609);
    const auto *qsg0_610 = buffer.data(qsg0 + 610);
    const auto *qsg0_612 = buffer.data(qsg0 + 612);
    const auto *qsg0_613 = buffer.data(qsg0 + 613);
    const auto *qsg0_614 = buffer.data(qsg0 + 614);
    const auto *qsg0_615 = buffer.data(qsg0 + 615);

    const auto *qsg1_564 = buffer.data(qsg1 + 564);
    const auto *qsg1_567 = buffer.data(qsg1 + 567);
    const auto *qsg1_568 = buffer.data(qsg1 + 568);
    const auto *qsg1_569 = buffer.data(qsg1 + 569);
    const auto *qsg1_570 = buffer.data(qsg1 + 570);
    const auto *qsg1_573 = buffer.data(qsg1 + 573);
    const auto *qsg1_575 = buffer.data(qsg1 + 575);
    const auto *qsg1_576 = buffer.data(qsg1 + 576);
    const auto *qsg1_579 = buffer.data(qsg1 + 579);
    const auto *qsg1_580 = buffer.data(qsg1 + 580);
    const auto *qsg1_582 = buffer.data(qsg1 + 582);
    const auto *qsg1_583 = buffer.data(qsg1 + 583);
    const auto *qsg1_584 = buffer.data(qsg1 + 584);
    const auto *qsg1_585 = buffer.data(qsg1 + 585);
    const auto *qsg1_588 = buffer.data(qsg1 + 588);
    const auto *qsg1_590 = buffer.data(qsg1 + 590);
    const auto *qsg1_591 = buffer.data(qsg1 + 591);
    const auto *qsg1_594 = buffer.data(qsg1 + 594);
    const auto *qsg1_595 = buffer.data(qsg1 + 595);
    const auto *qsg1_597 = buffer.data(qsg1 + 597);
    const auto *qsg1_598 = buffer.data(qsg1 + 598);
    const auto *qsg1_599 = buffer.data(qsg1 + 599);
    const auto *qsg1_600 = buffer.data(qsg1 + 600);
    const auto *qsg1_603 = buffer.data(qsg1 + 603);
    const auto *qsg1_605 = buffer.data(qsg1 + 605);
    const auto *qsg1_606 = buffer.data(qsg1 + 606);
    const auto *qsg1_609 = buffer.data(qsg1 + 609);
    const auto *qsg1_610 = buffer.data(qsg1 + 610);
    const auto *qsg1_612 = buffer.data(qsg1 + 612);
    const auto *qsg1_613 = buffer.data(qsg1 + 613);
    const auto *qsg1_614 = buffer.data(qsg1 + 614);
    const auto *qsg1_615 = buffer.data(qsg1 + 615);

    const auto *qsh_782 = buffer.data(qsh + 782);
    const auto *qsh_783 = buffer.data(qsh + 783);
    const auto *qsh_786 = buffer.data(qsh + 786);
    const auto *qsh_791 = buffer.data(qsh + 791);
    const auto *qsh_792 = buffer.data(qsh + 792);
    const auto *qsh_793 = buffer.data(qsh + 793);
    const auto *qsh_794 = buffer.data(qsh + 794);
    const auto *qsh_795 = buffer.data(qsh + 795);
    const auto *qsh_796 = buffer.data(qsh + 796);
    const auto *qsh_797 = buffer.data(qsh + 797);
    const auto *qsh_798 = buffer.data(qsh + 798);
    const auto *qsh_800 = buffer.data(qsh + 800);
    const auto *qsh_801 = buffer.data(qsh + 801);
    const auto *qsh_803 = buffer.data(qsh + 803);
    const auto *qsh_804 = buffer.data(qsh + 804);
    const auto *qsh_807 = buffer.data(qsh + 807);
    const auto *qsh_808 = buffer.data(qsh + 808);
    const auto *qsh_810 = buffer.data(qsh + 810);
    const auto *qsh_812 = buffer.data(qsh + 812);
    const auto *qsh_813 = buffer.data(qsh + 813);
    const auto *qsh_814 = buffer.data(qsh + 814);
    const auto *qsh_815 = buffer.data(qsh + 815);
    const auto *qsh_816 = buffer.data(qsh + 816);
    const auto *qsh_817 = buffer.data(qsh + 817);
    const auto *qsh_818 = buffer.data(qsh + 818);
    const auto *qsh_819 = buffer.data(qsh + 819);
    const auto *qsh_821 = buffer.data(qsh + 821);
    const auto *qsh_822 = buffer.data(qsh + 822);
    const auto *qsh_824 = buffer.data(qsh + 824);
    const auto *qsh_825 = buffer.data(qsh + 825);
    const auto *qsh_828 = buffer.data(qsh + 828);
    const auto *qsh_829 = buffer.data(qsh + 829);
    const auto *qsh_831 = buffer.data(qsh + 831);
    const auto *qsh_833 = buffer.data(qsh + 833);
    const auto *qsh_834 = buffer.data(qsh + 834);
    const auto *qsh_835 = buffer.data(qsh + 835);
    const auto *qsh_836 = buffer.data(qsh + 836);
    const auto *qsh_837 = buffer.data(qsh + 837);
    const auto *qsh_838 = buffer.data(qsh + 838);
    const auto *qsh_839 = buffer.data(qsh + 839);
    const auto *qsh_840 = buffer.data(qsh + 840);
    const auto *qsh_842 = buffer.data(qsh + 842);
    const auto *qsh_843 = buffer.data(qsh + 843);
    const auto *qsh_845 = buffer.data(qsh + 845);
    const auto *qsh_846 = buffer.data(qsh + 846);
    const auto *qsh_849 = buffer.data(qsh + 849);
    const auto *qsh_850 = buffer.data(qsh + 850);
    const auto *qsh_852 = buffer.data(qsh + 852);
    const auto *qsh_854 = buffer.data(qsh + 854);
    const auto *qsh_855 = buffer.data(qsh + 855);
    const auto *qsh_856 = buffer.data(qsh + 856);
    const auto *qsh_857 = buffer.data(qsh + 857);
    const auto *qsh_858 = buffer.data(qsh + 858);
    const auto *qsh_859 = buffer.data(qsh + 859);
    const auto *qsh_860 = buffer.data(qsh + 860);
    const auto *qsh_861 = buffer.data(qsh + 861);

#pragma omp simd aligned(t_1044, t_1045, t_1046, pa_z, pc_x, pc_y, pc_z, osi0_794, osh_614, \
                         osh_786, osi1_794, qsg0_564, qsg1_564, qsh_782, \
                         qsh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_21 * osh_614[k]
                    + f_3 * pc_y[k] * qsh_782[k];

        t_1045[k] = f_14 * osh_786[k]
                    + f_6 * qsg0_564[k]
                    - f_7 * qsg1_564[k]
                    + f_3 * pc_x[k] * qsh_786[k];

        t_1046[k] = pa_z[k] * osi0_794[k]
                    - f_10 * pc_z[k] * osi1_794[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pa_z, pc_y, pc_z, osi0_796, osh_594, osh_595, \
                         osh_618, osi1_796, qsh_783, qsh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = f_11 * osh_594[k]
                    + f_3 * pc_z[k] * qsh_783[k];

        t_1048[k] = pa_z[k] * osi0_796[k]
                    + f_12 * osh_595[k]
                    - f_10 * pc_z[k] * osi1_796[k];

        t_1049[k] = f_21 * osh_618[k]
                    + f_3 * pc_y[k] * qsh_786[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pc_x, osh_791, osh_792, osh_793, \
                         osh_794, qsg0_569, qsg1_569, qsh_791, qsh_792, qsh_793, \
                         qsh_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_14 * osh_791[k]
                    + f_4 * qsg0_569[k]
                    - f_5 * qsg1_569[k]
                    + f_3 * pc_x[k] * qsh_791[k];

        t_1051[k] = f_14 * osh_792[k]
                    + f_3 * pc_x[k] * qsh_792[k];

        t_1052[k] = f_14 * osh_793[k]
                    + f_3 * pc_x[k] * qsh_793[k];

        t_1053[k] = f_14 * osh_794[k]
                    + f_3 * pc_x[k] * qsh_794[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, pa_z, pc_x, pc_z, osi0_805, osh_795, \
                         osh_796, osh_797, osi1_805, qsh_795, qsh_796, \
                         qsh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_14 * osh_795[k]
                    + f_3 * pc_x[k] * qsh_795[k];

        t_1055[k] = f_14 * osh_796[k]
                    + f_3 * pc_x[k] * qsh_796[k];

        t_1056[k] = f_14 * osh_797[k]
                    + f_3 * pc_x[k] * qsh_797[k];

        t_1057[k] = pa_z[k] * osi0_805[k]
                    - f_10 * pc_z[k] * osi1_805[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_y, pc_z, osh_603, osh_626, osh_627, \
                         qsg0_567, qsg0_568, qsg1_567, qsg1_568, qsh_792, qsh_794, \
                         qsh_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_11 * osh_603[k]
                    + f_3 * pc_z[k] * qsh_792[k];

        t_1059[k] = f_21 * osh_626[k]
                    + f_8 * qsg0_567[k]
                    - f_9 * qsg1_567[k]
                    + f_3 * pc_y[k] * qsh_794[k];

        t_1060[k] = f_21 * osh_627[k]
                    + f_6 * qsg0_568[k]
                    - f_7 * qsg1_568[k]
                    + f_3 * pc_y[k] * qsh_795[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pc_y, pc_z, osh_608, osh_628, osh_629, \
                         qsg0_569, qsg1_569, qsh_796, qsh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_21 * osh_628[k]
                    + f_4 * qsg0_569[k]
                    - f_5 * qsg1_569[k]
                    + f_3 * pc_y[k] * qsh_796[k];

        t_1062[k] = f_21 * osh_629[k]
                    + f_3 * pc_y[k] * qsh_797[k];

        t_1063[k] = f_11 * osh_608[k]
                    + f_1 * qsg0_569[k]
                    - f_2 * qsg1_569[k]
                    + f_3 * pc_z[k] * qsh_797[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, pc_x, pc_y, pc_z, osh_609, osh_630, osh_798, \
                         qsg0_570, qsg1_570, qsh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_14 * osh_798[k]
                    + f_1 * qsg0_570[k]
                    - f_2 * qsg1_570[k]
                    + f_3 * pc_x[k] * qsh_798[k];

        t_1065[k] = f_23 * osh_630[k]
                    + f_3 * pc_y[k] * qsh_798[k];

        t_1066[k] = f_12 * osh_609[k]
                    + f_3 * pc_z[k] * qsh_798[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pc_x, pc_y, osh_632, osh_801, osh_803, \
                         qsg0_573, qsg0_575, qsg1_573, qsg1_575, qsh_800, qsh_801, \
                         qsh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_14 * osh_801[k]
                    + f_8 * qsg0_573[k]
                    - f_9 * qsg1_573[k]
                    + f_3 * pc_x[k] * qsh_801[k];

        t_1068[k] = f_23 * osh_632[k]
                    + f_3 * pc_y[k] * qsh_800[k];

        t_1069[k] = f_14 * osh_803[k]
                    + f_8 * qsg0_575[k]
                    - f_9 * qsg1_575[k]
                    + f_3 * pc_x[k] * qsh_803[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, osh_612, osh_635, osh_804, \
                         qsg0_576, qsg1_576, qsh_801, qsh_803, \
                         qsh_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_14 * osh_804[k]
                    + f_6 * qsg0_576[k]
                    - f_7 * qsg1_576[k]
                    + f_3 * pc_x[k] * qsh_804[k];

        t_1071[k] = f_12 * osh_612[k]
                    + f_3 * pc_z[k] * qsh_801[k];

        t_1072[k] = f_23 * osh_635[k]
                    + f_3 * pc_y[k] * qsh_803[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, pc_z, osh_615, osh_807, osh_808, \
                         qsg0_579, qsg0_580, qsg1_579, qsg1_580, qsh_804, qsh_807, \
                         qsh_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_14 * osh_807[k]
                    + f_6 * qsg0_579[k]
                    - f_7 * qsg1_579[k]
                    + f_3 * pc_x[k] * qsh_807[k];

        t_1074[k] = f_14 * osh_808[k]
                    + f_4 * qsg0_580[k]
                    - f_5 * qsg1_580[k]
                    + f_3 * pc_x[k] * qsh_808[k];

        t_1075[k] = f_12 * osh_615[k]
                    + f_3 * pc_z[k] * qsh_804[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, pc_y, osh_639, osh_810, osh_812, \
                         qsg0_582, qsg0_584, qsg1_582, qsg1_584, qsh_807, qsh_810, \
                         qsh_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_14 * osh_810[k]
                    + f_4 * qsg0_582[k]
                    - f_5 * qsg1_582[k]
                    + f_3 * pc_x[k] * qsh_810[k];

        t_1077[k] = f_23 * osh_639[k]
                    + f_3 * pc_y[k] * qsh_807[k];

        t_1078[k] = f_14 * osh_812[k]
                    + f_4 * qsg0_584[k]
                    - f_5 * qsg1_584[k]
                    + f_3 * pc_x[k] * qsh_812[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, t_1083, pc_x, osh_813, osh_814, \
                         osh_815, osh_816, osh_817, qsh_813, qsh_814, qsh_815, qsh_816, \
                         qsh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_14 * osh_813[k]
                    + f_3 * pc_x[k] * qsh_813[k];

        t_1080[k] = f_14 * osh_814[k]
                    + f_3 * pc_x[k] * qsh_814[k];

        t_1081[k] = f_14 * osh_815[k]
                    + f_3 * pc_x[k] * qsh_815[k];

        t_1082[k] = f_14 * osh_816[k]
                    + f_3 * pc_x[k] * qsh_816[k];

        t_1083[k] = f_14 * osh_817[k]
                    + f_3 * pc_x[k] * qsh_817[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, pc_x, pc_y, pc_z, osh_624, osh_645, osh_818, \
                         qsg0_580, qsg1_580, qsh_813, qsh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_14 * osh_818[k]
                    + f_3 * pc_x[k] * qsh_818[k];

        t_1085[k] = f_23 * osh_645[k]
                    + f_1 * qsg0_580[k]
                    - f_2 * qsg1_580[k]
                    + f_3 * pc_y[k] * qsh_813[k];

        t_1086[k] = f_12 * osh_624[k]
                    + f_3 * pc_z[k] * qsh_813[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, osh_647, osh_648, osh_649, qsg0_582, \
                         qsg0_583, qsg0_584, qsg1_582, qsg1_583, qsg1_584, qsh_815, qsh_816, \
                         qsh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_23 * osh_647[k]
                    + f_8 * qsg0_582[k]
                    - f_9 * qsg1_582[k]
                    + f_3 * pc_y[k] * qsh_815[k];

        t_1088[k] = f_23 * osh_648[k]
                    + f_6 * qsg0_583[k]
                    - f_7 * qsg1_583[k]
                    + f_3 * pc_y[k] * qsh_816[k];

        t_1089[k] = f_23 * osh_649[k]
                    + f_4 * qsg0_584[k]
                    - f_5 * qsg1_584[k]
                    + f_3 * pc_y[k] * qsh_817[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, osh_629, osh_650, osh_819, \
                         qsg0_584, qsg0_585, qsg1_584, qsg1_585, qsh_818, \
                         qsh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_23 * osh_650[k]
                    + f_3 * pc_y[k] * qsh_818[k];

        t_1091[k] = f_12 * osh_629[k]
                    + f_1 * qsg0_584[k]
                    - f_2 * qsg1_584[k]
                    + f_3 * pc_z[k] * qsh_818[k];

        t_1092[k] = f_14 * osh_819[k]
                    + f_1 * qsg0_585[k]
                    - f_2 * qsg1_585[k]
                    + f_3 * pc_x[k] * qsh_819[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pc_x, pc_y, pc_z, osh_630, osh_651, \
                         osh_653, osh_822, qsg0_588, qsg1_588, qsh_819, qsh_821, \
                         qsh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_22 * osh_651[k]
                    + f_3 * pc_y[k] * qsh_819[k];

        t_1094[k] = f_13 * osh_630[k]
                    + f_3 * pc_z[k] * qsh_819[k];

        t_1095[k] = f_14 * osh_822[k]
                    + f_8 * qsg0_588[k]
                    - f_9 * qsg1_588[k]
                    + f_3 * pc_x[k] * qsh_822[k];

        t_1096[k] = f_22 * osh_653[k]
                    + f_3 * pc_y[k] * qsh_821[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_z, osh_633, osh_824, osh_825, \
                         qsg0_590, qsg0_591, qsg1_590, qsg1_591, qsh_822, qsh_824, \
                         qsh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_14 * osh_824[k]
                    + f_8 * qsg0_590[k]
                    - f_9 * qsg1_590[k]
                    + f_3 * pc_x[k] * qsh_824[k];

        t_1098[k] = f_14 * osh_825[k]
                    + f_6 * qsg0_591[k]
                    - f_7 * qsg1_591[k]
                    + f_3 * pc_x[k] * qsh_825[k];

        t_1099[k] = f_13 * osh_633[k]
                    + f_3 * pc_z[k] * qsh_822[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_y, osh_656, osh_828, osh_829, \
                         qsg0_594, qsg0_595, qsg1_594, qsg1_595, qsh_824, qsh_828, \
                         qsh_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_22 * osh_656[k]
                    + f_3 * pc_y[k] * qsh_824[k];

        t_1101[k] = f_14 * osh_828[k]
                    + f_6 * qsg0_594[k]
                    - f_7 * qsg1_594[k]
                    + f_3 * pc_x[k] * qsh_828[k];

        t_1102[k] = f_14 * osh_829[k]
                    + f_4 * qsg0_595[k]
                    - f_5 * qsg1_595[k]
                    + f_3 * pc_x[k] * qsh_829[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, pc_y, pc_z, osh_636, osh_660, osh_831, \
                         qsg0_597, qsg1_597, qsh_825, qsh_828, \
                         qsh_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_13 * osh_636[k]
                    + f_3 * pc_z[k] * qsh_825[k];

        t_1104[k] = f_14 * osh_831[k]
                    + f_4 * qsg0_597[k]
                    - f_5 * qsg1_597[k]
                    + f_3 * pc_x[k] * qsh_831[k];

        t_1105[k] = f_22 * osh_660[k]
                    + f_3 * pc_y[k] * qsh_828[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, osh_833, osh_834, osh_835, \
                         osh_836, qsg0_599, qsg1_599, qsh_833, qsh_834, qsh_835, \
                         qsh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_14 * osh_833[k]
                    + f_4 * qsg0_599[k]
                    - f_5 * qsg1_599[k]
                    + f_3 * pc_x[k] * qsh_833[k];

        t_1107[k] = f_14 * osh_834[k]
                    + f_3 * pc_x[k] * qsh_834[k];

        t_1108[k] = f_14 * osh_835[k]
                    + f_3 * pc_x[k] * qsh_835[k];

        t_1109[k] = f_14 * osh_836[k]
                    + f_3 * pc_x[k] * qsh_836[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pc_x, pc_y, osh_666, osh_837, \
                         osh_838, osh_839, qsg0_595, qsg1_595, qsh_834, qsh_837, qsh_838, \
                         qsh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_14 * osh_837[k]
                    + f_3 * pc_x[k] * qsh_837[k];

        t_1111[k] = f_14 * osh_838[k]
                    + f_3 * pc_x[k] * qsh_838[k];

        t_1112[k] = f_14 * osh_839[k]
                    + f_3 * pc_x[k] * qsh_839[k];

        t_1113[k] = f_22 * osh_666[k]
                    + f_1 * qsg0_595[k]
                    - f_2 * qsg1_595[k]
                    + f_3 * pc_y[k] * qsh_834[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_y, pc_z, osh_645, osh_668, osh_669, \
                         qsg0_597, qsg0_598, qsg1_597, qsg1_598, qsh_834, qsh_836, \
                         qsh_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_13 * osh_645[k]
                    + f_3 * pc_z[k] * qsh_834[k];

        t_1115[k] = f_22 * osh_668[k]
                    + f_8 * qsg0_597[k]
                    - f_9 * qsg1_597[k]
                    + f_3 * pc_y[k] * qsh_836[k];

        t_1116[k] = f_22 * osh_669[k]
                    + f_6 * qsg0_598[k]
                    - f_7 * qsg1_598[k]
                    + f_3 * pc_y[k] * qsh_837[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, pc_y, pc_z, osh_650, osh_670, osh_671, \
                         qsg0_599, qsg1_599, qsh_838, qsh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_22 * osh_670[k]
                    + f_4 * qsg0_599[k]
                    - f_5 * qsg1_599[k]
                    + f_3 * pc_y[k] * qsh_838[k];

        t_1118[k] = f_22 * osh_671[k]
                    + f_3 * pc_y[k] * qsh_839[k];

        t_1119[k] = f_13 * osh_650[k]
                    + f_1 * qsg0_599[k]
                    - f_2 * qsg1_599[k]
                    + f_3 * pc_z[k] * qsh_839[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, pc_x, pc_y, pc_z, osh_651, osh_672, osh_840, \
                         qsg0_600, qsg1_600, qsh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = f_14 * osh_840[k]
                    + f_1 * qsg0_600[k]
                    - f_2 * qsg1_600[k]
                    + f_3 * pc_x[k] * qsh_840[k];

        t_1121[k] = f_14 * osh_672[k]
                    + f_3 * pc_y[k] * qsh_840[k];

        t_1122[k] = f_14 * osh_651[k]
                    + f_3 * pc_z[k] * qsh_840[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, pc_x, pc_y, osh_674, osh_843, osh_845, \
                         qsg0_603, qsg0_605, qsg1_603, qsg1_605, qsh_842, qsh_843, \
                         qsh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_14 * osh_843[k]
                    + f_8 * qsg0_603[k]
                    - f_9 * qsg1_603[k]
                    + f_3 * pc_x[k] * qsh_843[k];

        t_1124[k] = f_14 * osh_674[k]
                    + f_3 * pc_y[k] * qsh_842[k];

        t_1125[k] = f_14 * osh_845[k]
                    + f_8 * qsg0_605[k]
                    - f_9 * qsg1_605[k]
                    + f_3 * pc_x[k] * qsh_845[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, pc_x, pc_y, pc_z, osh_654, osh_677, osh_846, \
                         qsg0_606, qsg1_606, qsh_843, qsh_845, \
                         qsh_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_14 * osh_846[k]
                    + f_6 * qsg0_606[k]
                    - f_7 * qsg1_606[k]
                    + f_3 * pc_x[k] * qsh_846[k];

        t_1127[k] = f_14 * osh_654[k]
                    + f_3 * pc_z[k] * qsh_843[k];

        t_1128[k] = f_14 * osh_677[k]
                    + f_3 * pc_y[k] * qsh_845[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pc_x, pc_z, osh_657, osh_849, osh_850, \
                         qsg0_609, qsg0_610, qsg1_609, qsg1_610, qsh_846, qsh_849, \
                         qsh_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_14 * osh_849[k]
                    + f_6 * qsg0_609[k]
                    - f_7 * qsg1_609[k]
                    + f_3 * pc_x[k] * qsh_849[k];

        t_1130[k] = f_14 * osh_850[k]
                    + f_4 * qsg0_610[k]
                    - f_5 * qsg1_610[k]
                    + f_3 * pc_x[k] * qsh_850[k];

        t_1131[k] = f_14 * osh_657[k]
                    + f_3 * pc_z[k] * qsh_846[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, pc_x, pc_y, osh_681, osh_852, osh_854, \
                         qsg0_612, qsg0_614, qsg1_612, qsg1_614, qsh_849, qsh_852, \
                         qsh_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_14 * osh_852[k]
                    + f_4 * qsg0_612[k]
                    - f_5 * qsg1_612[k]
                    + f_3 * pc_x[k] * qsh_852[k];

        t_1133[k] = f_14 * osh_681[k]
                    + f_3 * pc_y[k] * qsh_849[k];

        t_1134[k] = f_14 * osh_854[k]
                    + f_4 * qsg0_614[k]
                    - f_5 * qsg1_614[k]
                    + f_3 * pc_x[k] * qsh_854[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, t_1138, t_1139, pc_x, osh_855, osh_856, \
                         osh_857, osh_858, osh_859, qsh_855, qsh_856, qsh_857, qsh_858, \
                         qsh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = f_14 * osh_855[k]
                    + f_3 * pc_x[k] * qsh_855[k];

        t_1136[k] = f_14 * osh_856[k]
                    + f_3 * pc_x[k] * qsh_856[k];

        t_1137[k] = f_14 * osh_857[k]
                    + f_3 * pc_x[k] * qsh_857[k];

        t_1138[k] = f_14 * osh_858[k]
                    + f_3 * pc_x[k] * qsh_858[k];

        t_1139[k] = f_14 * osh_859[k]
                    + f_3 * pc_x[k] * qsh_859[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, pc_x, pc_y, pc_z, osh_666, osh_687, osh_860, \
                         qsg0_610, qsg1_610, qsh_855, qsh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_14 * osh_860[k]
                    + f_3 * pc_x[k] * qsh_860[k];

        t_1141[k] = f_14 * osh_687[k]
                    + f_1 * qsg0_610[k]
                    - f_2 * qsg1_610[k]
                    + f_3 * pc_y[k] * qsh_855[k];

        t_1142[k] = f_14 * osh_666[k]
                    + f_3 * pc_z[k] * qsh_855[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, pc_y, osh_689, osh_690, osh_691, qsg0_612, \
                         qsg0_613, qsg0_614, qsg1_612, qsg1_613, qsg1_614, qsh_857, qsh_858, \
                         qsh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_14 * osh_689[k]
                    + f_8 * qsg0_612[k]
                    - f_9 * qsg1_612[k]
                    + f_3 * pc_y[k] * qsh_857[k];

        t_1144[k] = f_14 * osh_690[k]
                    + f_6 * qsg0_613[k]
                    - f_7 * qsg1_613[k]
                    + f_3 * pc_y[k] * qsh_858[k];

        t_1145[k] = f_14 * osh_691[k]
                    + f_4 * qsg0_614[k]
                    - f_5 * qsg1_614[k]
                    + f_3 * pc_y[k] * qsh_859[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pc_x, pc_y, pc_z, osh_671, osh_692, osh_861, \
                         qsg0_614, qsg0_615, qsg1_614, qsg1_615, qsh_860, \
                         qsh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_14 * osh_692[k]
                    + f_3 * pc_y[k] * qsh_860[k];

        t_1147[k] = f_14 * osh_671[k]
                    + f_1 * qsg0_614[k]
                    - f_2 * qsg1_614[k]
                    + f_3 * pc_z[k] * qsh_860[k];

        t_1148[k] = f_14 * osh_861[k]
                    + f_1 * qsg0_615[k]
                    - f_2 * qsg1_615[k]
                    + f_3 * pc_x[k] * qsh_861[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsg0, const size_t qsg1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_980 = buffer.data(osi0 + 980);
    const auto *osi0_983 = buffer.data(osi0 + 983);
    const auto *osi0_985 = buffer.data(osi0 + 985);
    const auto *osi0_986 = buffer.data(osi0 + 986);
    const auto *osi0_989 = buffer.data(osi0 + 989);
    const auto *osi0_990 = buffer.data(osi0 + 990);
    const auto *osi0_992 = buffer.data(osi0 + 992);
    const auto *osi0_994 = buffer.data(osi0 + 994);
    const auto *osi0_1007 = buffer.data(osi0 + 1007);

    const auto *osh_672 = buffer.data(osh + 672);
    const auto *osh_675 = buffer.data(osh + 675);
    const auto *osh_678 = buffer.data(osh + 678);
    const auto *osh_687 = buffer.data(osh + 687);
    const auto *osh_692 = buffer.data(osh + 692);
    const auto *osh_693 = buffer.data(osh + 693);
    const auto *osh_695 = buffer.data(osh + 695);
    const auto *osh_696 = buffer.data(osh + 696);
    const auto *osh_698 = buffer.data(osh + 698);
    const auto *osh_699 = buffer.data(osh + 699);
    const auto *osh_702 = buffer.data(osh + 702);
    const auto *osh_708 = buffer.data(osh + 708);
    const auto *osh_710 = buffer.data(osh + 710);
    const auto *osh_711 = buffer.data(osh + 711);
    const auto *osh_712 = buffer.data(osh + 712);
    const auto *osh_713 = buffer.data(osh + 713);
    const auto *osh_714 = buffer.data(osh + 714);
    const auto *osh_716 = buffer.data(osh + 716);
    const auto *osh_717 = buffer.data(osh + 717);
    const auto *osh_719 = buffer.data(osh + 719);
    const auto *osh_720 = buffer.data(osh + 720);
    const auto *osh_723 = buffer.data(osh + 723);
    const auto *osh_729 = buffer.data(osh + 729);
    const auto *osh_731 = buffer.data(osh + 731);
    const auto *osh_732 = buffer.data(osh + 732);
    const auto *osh_733 = buffer.data(osh + 733);
    const auto *osh_734 = buffer.data(osh + 734);
    const auto *osh_735 = buffer.data(osh + 735);
    const auto *osh_736 = buffer.data(osh + 736);
    const auto *osh_737 = buffer.data(osh + 737);
    const auto *osh_738 = buffer.data(osh + 738);
    const auto *osh_740 = buffer.data(osh + 740);
    const auto *osh_741 = buffer.data(osh + 741);
    const auto *osh_743 = buffer.data(osh + 743);
    const auto *osh_744 = buffer.data(osh + 744);
    const auto *osh_750 = buffer.data(osh + 750);
    const auto *osh_752 = buffer.data(osh + 752);
    const auto *osh_753 = buffer.data(osh + 753);
    const auto *osh_754 = buffer.data(osh + 754);
    const auto *osh_755 = buffer.data(osh + 755);
    const auto *osh_864 = buffer.data(osh + 864);
    const auto *osh_866 = buffer.data(osh + 866);
    const auto *osh_867 = buffer.data(osh + 867);
    const auto *osh_870 = buffer.data(osh + 870);
    const auto *osh_871 = buffer.data(osh + 871);
    const auto *osh_873 = buffer.data(osh + 873);
    const auto *osh_875 = buffer.data(osh + 875);
    const auto *osh_876 = buffer.data(osh + 876);
    const auto *osh_877 = buffer.data(osh + 877);
    const auto *osh_878 = buffer.data(osh + 878);
    const auto *osh_879 = buffer.data(osh + 879);
    const auto *osh_880 = buffer.data(osh + 880);
    const auto *osh_881 = buffer.data(osh + 881);
    const auto *osh_882 = buffer.data(osh + 882);
    const auto *osh_885 = buffer.data(osh + 885);
    const auto *osh_887 = buffer.data(osh + 887);
    const auto *osh_888 = buffer.data(osh + 888);
    const auto *osh_891 = buffer.data(osh + 891);
    const auto *osh_892 = buffer.data(osh + 892);
    const auto *osh_894 = buffer.data(osh + 894);
    const auto *osh_896 = buffer.data(osh + 896);
    const auto *osh_897 = buffer.data(osh + 897);
    const auto *osh_898 = buffer.data(osh + 898);
    const auto *osh_899 = buffer.data(osh + 899);
    const auto *osh_900 = buffer.data(osh + 900);
    const auto *osh_901 = buffer.data(osh + 901);
    const auto *osh_902 = buffer.data(osh + 902);
    const auto *osh_918 = buffer.data(osh + 918);
    const auto *osh_919 = buffer.data(osh + 919);
    const auto *osh_920 = buffer.data(osh + 920);
    const auto *osh_921 = buffer.data(osh + 921);
    const auto *osh_922 = buffer.data(osh + 922);
    const auto *osh_923 = buffer.data(osh + 923);
    const auto *osh_924 = buffer.data(osh + 924);
    const auto *osh_929 = buffer.data(osh + 929);
    const auto *osh_933 = buffer.data(osh + 933);
    const auto *osh_938 = buffer.data(osh + 938);
    const auto *osh_939 = buffer.data(osh + 939);
    const auto *osh_940 = buffer.data(osh + 940);
    const auto *osh_941 = buffer.data(osh + 941);
    const auto *osh_942 = buffer.data(osh + 942);
    const auto *osh_944 = buffer.data(osh + 944);

    const auto *osi1_980 = buffer.data(osi1 + 980);
    const auto *osi1_983 = buffer.data(osi1 + 983);
    const auto *osi1_985 = buffer.data(osi1 + 985);
    const auto *osi1_986 = buffer.data(osi1 + 986);
    const auto *osi1_989 = buffer.data(osi1 + 989);
    const auto *osi1_990 = buffer.data(osi1 + 990);
    const auto *osi1_992 = buffer.data(osi1 + 992);
    const auto *osi1_994 = buffer.data(osi1 + 994);
    const auto *osi1_1007 = buffer.data(osi1 + 1007);

    const auto *qsg0_618 = buffer.data(qsg0 + 618);
    const auto *qsg0_620 = buffer.data(qsg0 + 620);
    const auto *qsg0_621 = buffer.data(qsg0 + 621);
    const auto *qsg0_624 = buffer.data(qsg0 + 624);
    const auto *qsg0_625 = buffer.data(qsg0 + 625);
    const auto *qsg0_627 = buffer.data(qsg0 + 627);
    const auto *qsg0_628 = buffer.data(qsg0 + 628);
    const auto *qsg0_629 = buffer.data(qsg0 + 629);
    const auto *qsg0_630 = buffer.data(qsg0 + 630);
    const auto *qsg0_633 = buffer.data(qsg0 + 633);
    const auto *qsg0_635 = buffer.data(qsg0 + 635);
    const auto *qsg0_636 = buffer.data(qsg0 + 636);
    const auto *qsg0_639 = buffer.data(qsg0 + 639);
    const auto *qsg0_640 = buffer.data(qsg0 + 640);
    const auto *qsg0_642 = buffer.data(qsg0 + 642);
    const auto *qsg0_643 = buffer.data(qsg0 + 643);
    const auto *qsg0_644 = buffer.data(qsg0 + 644);
    const auto *qsg0_655 = buffer.data(qsg0 + 655);
    const auto *qsg0_657 = buffer.data(qsg0 + 657);
    const auto *qsg0_658 = buffer.data(qsg0 + 658);
    const auto *qsg0_659 = buffer.data(qsg0 + 659);
    const auto *qsg0_660 = buffer.data(qsg0 + 660);
    const auto *qsg0_661 = buffer.data(qsg0 + 661);
    const auto *qsg0_662 = buffer.data(qsg0 + 662);
    const auto *qsg0_663 = buffer.data(qsg0 + 663);
    const auto *qsg0_664 = buffer.data(qsg0 + 664);
    const auto *qsg0_665 = buffer.data(qsg0 + 665);
    const auto *qsg0_669 = buffer.data(qsg0 + 669);
    const auto *qsg0_670 = buffer.data(qsg0 + 670);
    const auto *qsg0_671 = buffer.data(qsg0 + 671);
    const auto *qsg0_672 = buffer.data(qsg0 + 672);
    const auto *qsg0_673 = buffer.data(qsg0 + 673);
    const auto *qsg0_674 = buffer.data(qsg0 + 674);

    const auto *qsg1_618 = buffer.data(qsg1 + 618);
    const auto *qsg1_620 = buffer.data(qsg1 + 620);
    const auto *qsg1_621 = buffer.data(qsg1 + 621);
    const auto *qsg1_624 = buffer.data(qsg1 + 624);
    const auto *qsg1_625 = buffer.data(qsg1 + 625);
    const auto *qsg1_627 = buffer.data(qsg1 + 627);
    const auto *qsg1_628 = buffer.data(qsg1 + 628);
    const auto *qsg1_629 = buffer.data(qsg1 + 629);
    const auto *qsg1_630 = buffer.data(qsg1 + 630);
    const auto *qsg1_633 = buffer.data(qsg1 + 633);
    const auto *qsg1_635 = buffer.data(qsg1 + 635);
    const auto *qsg1_636 = buffer.data(qsg1 + 636);
    const auto *qsg1_639 = buffer.data(qsg1 + 639);
    const auto *qsg1_640 = buffer.data(qsg1 + 640);
    const auto *qsg1_642 = buffer.data(qsg1 + 642);
    const auto *qsg1_643 = buffer.data(qsg1 + 643);
    const auto *qsg1_644 = buffer.data(qsg1 + 644);
    const auto *qsg1_655 = buffer.data(qsg1 + 655);
    const auto *qsg1_657 = buffer.data(qsg1 + 657);
    const auto *qsg1_658 = buffer.data(qsg1 + 658);
    const auto *qsg1_659 = buffer.data(qsg1 + 659);
    const auto *qsg1_660 = buffer.data(qsg1 + 660);
    const auto *qsg1_661 = buffer.data(qsg1 + 661);
    const auto *qsg1_662 = buffer.data(qsg1 + 662);
    const auto *qsg1_663 = buffer.data(qsg1 + 663);
    const auto *qsg1_664 = buffer.data(qsg1 + 664);
    const auto *qsg1_665 = buffer.data(qsg1 + 665);
    const auto *qsg1_669 = buffer.data(qsg1 + 669);
    const auto *qsg1_670 = buffer.data(qsg1 + 670);
    const auto *qsg1_671 = buffer.data(qsg1 + 671);
    const auto *qsg1_672 = buffer.data(qsg1 + 672);
    const auto *qsg1_673 = buffer.data(qsg1 + 673);
    const auto *qsg1_674 = buffer.data(qsg1 + 674);

    const auto *qsh_861 = buffer.data(qsh + 861);
    const auto *qsh_863 = buffer.data(qsh + 863);
    const auto *qsh_864 = buffer.data(qsh + 864);
    const auto *qsh_866 = buffer.data(qsh + 866);
    const auto *qsh_867 = buffer.data(qsh + 867);
    const auto *qsh_870 = buffer.data(qsh + 870);
    const auto *qsh_871 = buffer.data(qsh + 871);
    const auto *qsh_873 = buffer.data(qsh + 873);
    const auto *qsh_875 = buffer.data(qsh + 875);
    const auto *qsh_876 = buffer.data(qsh + 876);
    const auto *qsh_877 = buffer.data(qsh + 877);
    const auto *qsh_878 = buffer.data(qsh + 878);
    const auto *qsh_879 = buffer.data(qsh + 879);
    const auto *qsh_880 = buffer.data(qsh + 880);
    const auto *qsh_881 = buffer.data(qsh + 881);
    const auto *qsh_882 = buffer.data(qsh + 882);
    const auto *qsh_884 = buffer.data(qsh + 884);
    const auto *qsh_885 = buffer.data(qsh + 885);
    const auto *qsh_887 = buffer.data(qsh + 887);
    const auto *qsh_888 = buffer.data(qsh + 888);
    const auto *qsh_891 = buffer.data(qsh + 891);
    const auto *qsh_892 = buffer.data(qsh + 892);
    const auto *qsh_894 = buffer.data(qsh + 894);
    const auto *qsh_896 = buffer.data(qsh + 896);
    const auto *qsh_897 = buffer.data(qsh + 897);
    const auto *qsh_898 = buffer.data(qsh + 898);
    const auto *qsh_899 = buffer.data(qsh + 899);
    const auto *qsh_900 = buffer.data(qsh + 900);
    const auto *qsh_901 = buffer.data(qsh + 901);
    const auto *qsh_902 = buffer.data(qsh + 902);
    const auto *qsh_903 = buffer.data(qsh + 903);
    const auto *qsh_905 = buffer.data(qsh + 905);
    const auto *qsh_906 = buffer.data(qsh + 906);
    const auto *qsh_908 = buffer.data(qsh + 908);
    const auto *qsh_909 = buffer.data(qsh + 909);
    const auto *qsh_912 = buffer.data(qsh + 912);
    const auto *qsh_918 = buffer.data(qsh + 918);
    const auto *qsh_919 = buffer.data(qsh + 919);
    const auto *qsh_920 = buffer.data(qsh + 920);
    const auto *qsh_921 = buffer.data(qsh + 921);
    const auto *qsh_922 = buffer.data(qsh + 922);
    const auto *qsh_923 = buffer.data(qsh + 923);
    const auto *qsh_924 = buffer.data(qsh + 924);
    const auto *qsh_925 = buffer.data(qsh + 925);
    const auto *qsh_926 = buffer.data(qsh + 926);
    const auto *qsh_927 = buffer.data(qsh + 927);
    const auto *qsh_928 = buffer.data(qsh + 928);
    const auto *qsh_929 = buffer.data(qsh + 929);
    const auto *qsh_930 = buffer.data(qsh + 930);
    const auto *qsh_931 = buffer.data(qsh + 931);
    const auto *qsh_932 = buffer.data(qsh + 932);
    const auto *qsh_933 = buffer.data(qsh + 933);
    const auto *qsh_938 = buffer.data(qsh + 938);
    const auto *qsh_939 = buffer.data(qsh + 939);
    const auto *qsh_940 = buffer.data(qsh + 940);
    const auto *qsh_941 = buffer.data(qsh + 941);
    const auto *qsh_942 = buffer.data(qsh + 942);
    const auto *qsh_943 = buffer.data(qsh + 943);
    const auto *qsh_944 = buffer.data(qsh + 944);

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, osh_672, osh_693, \
                         osh_695, osh_864, qsg0_618, qsg1_618, qsh_861, qsh_863, \
                         qsh_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_13 * osh_693[k]
                    + f_3 * pc_y[k] * qsh_861[k];

        t_1150[k] = f_22 * osh_672[k]
                    + f_3 * pc_z[k] * qsh_861[k];

        t_1151[k] = f_14 * osh_864[k]
                    + f_8 * qsg0_618[k]
                    - f_9 * qsg1_618[k]
                    + f_3 * pc_x[k] * qsh_864[k];

        t_1152[k] = f_13 * osh_695[k]
                    + f_3 * pc_y[k] * qsh_863[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, pc_x, pc_z, osh_675, osh_866, osh_867, \
                         qsg0_620, qsg0_621, qsg1_620, qsg1_621, qsh_864, qsh_866, \
                         qsh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_14 * osh_866[k]
                    + f_8 * qsg0_620[k]
                    - f_9 * qsg1_620[k]
                    + f_3 * pc_x[k] * qsh_866[k];

        t_1154[k] = f_14 * osh_867[k]
                    + f_6 * qsg0_621[k]
                    - f_7 * qsg1_621[k]
                    + f_3 * pc_x[k] * qsh_867[k];

        t_1155[k] = f_22 * osh_675[k]
                    + f_3 * pc_z[k] * qsh_864[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, pc_x, pc_y, osh_698, osh_870, osh_871, \
                         qsg0_624, qsg0_625, qsg1_624, qsg1_625, qsh_866, qsh_870, \
                         qsh_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_13 * osh_698[k]
                    + f_3 * pc_y[k] * qsh_866[k];

        t_1157[k] = f_14 * osh_870[k]
                    + f_6 * qsg0_624[k]
                    - f_7 * qsg1_624[k]
                    + f_3 * pc_x[k] * qsh_870[k];

        t_1158[k] = f_14 * osh_871[k]
                    + f_4 * qsg0_625[k]
                    - f_5 * qsg1_625[k]
                    + f_3 * pc_x[k] * qsh_871[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, pc_x, pc_y, pc_z, osh_678, osh_702, osh_873, \
                         qsg0_627, qsg1_627, qsh_867, qsh_870, \
                         qsh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_22 * osh_678[k]
                    + f_3 * pc_z[k] * qsh_867[k];

        t_1160[k] = f_14 * osh_873[k]
                    + f_4 * qsg0_627[k]
                    - f_5 * qsg1_627[k]
                    + f_3 * pc_x[k] * qsh_873[k];

        t_1161[k] = f_13 * osh_702[k]
                    + f_3 * pc_y[k] * qsh_870[k];
    }

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, pc_x, osh_875, osh_876, osh_877, \
                         osh_878, qsg0_629, qsg1_629, qsh_875, qsh_876, qsh_877, \
                         qsh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = f_14 * osh_875[k]
                    + f_4 * qsg0_629[k]
                    - f_5 * qsg1_629[k]
                    + f_3 * pc_x[k] * qsh_875[k];

        t_1163[k] = f_14 * osh_876[k]
                    + f_3 * pc_x[k] * qsh_876[k];

        t_1164[k] = f_14 * osh_877[k]
                    + f_3 * pc_x[k] * qsh_877[k];

        t_1165[k] = f_14 * osh_878[k]
                    + f_3 * pc_x[k] * qsh_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, pc_x, pc_y, osh_708, osh_879, \
                         osh_880, osh_881, qsg0_625, qsg1_625, qsh_876, qsh_879, qsh_880, \
                         qsh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_14 * osh_879[k]
                    + f_3 * pc_x[k] * qsh_879[k];

        t_1167[k] = f_14 * osh_880[k]
                    + f_3 * pc_x[k] * qsh_880[k];

        t_1168[k] = f_14 * osh_881[k]
                    + f_3 * pc_x[k] * qsh_881[k];

        t_1169[k] = f_13 * osh_708[k]
                    + f_1 * qsg0_625[k]
                    - f_2 * qsg1_625[k]
                    + f_3 * pc_y[k] * qsh_876[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, pc_y, pc_z, osh_687, osh_710, osh_711, \
                         qsg0_627, qsg0_628, qsg1_627, qsg1_628, qsh_876, qsh_878, \
                         qsh_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_22 * osh_687[k]
                    + f_3 * pc_z[k] * qsh_876[k];

        t_1171[k] = f_13 * osh_710[k]
                    + f_8 * qsg0_627[k]
                    - f_9 * qsg1_627[k]
                    + f_3 * pc_y[k] * qsh_878[k];

        t_1172[k] = f_13 * osh_711[k]
                    + f_6 * qsg0_628[k]
                    - f_7 * qsg1_628[k]
                    + f_3 * pc_y[k] * qsh_879[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, pc_y, pc_z, osh_692, osh_712, osh_713, \
                         qsg0_629, qsg1_629, qsh_880, qsh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_13 * osh_712[k]
                    + f_4 * qsg0_629[k]
                    - f_5 * qsg1_629[k]
                    + f_3 * pc_y[k] * qsh_880[k];

        t_1174[k] = f_13 * osh_713[k]
                    + f_3 * pc_y[k] * qsh_881[k];

        t_1175[k] = f_22 * osh_692[k]
                    + f_1 * qsg0_629[k]
                    - f_2 * qsg1_629[k]
                    + f_3 * pc_z[k] * qsh_881[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, pc_x, pc_y, pc_z, osh_693, osh_714, osh_882, \
                         qsg0_630, qsg1_630, qsh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_14 * osh_882[k]
                    + f_1 * qsg0_630[k]
                    - f_2 * qsg1_630[k]
                    + f_3 * pc_x[k] * qsh_882[k];

        t_1177[k] = f_12 * osh_714[k]
                    + f_3 * pc_y[k] * qsh_882[k];

        t_1178[k] = f_23 * osh_693[k]
                    + f_3 * pc_z[k] * qsh_882[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pc_x, pc_y, osh_716, osh_885, osh_887, \
                         qsg0_633, qsg0_635, qsg1_633, qsg1_635, qsh_884, qsh_885, \
                         qsh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = f_14 * osh_885[k]
                    + f_8 * qsg0_633[k]
                    - f_9 * qsg1_633[k]
                    + f_3 * pc_x[k] * qsh_885[k];

        t_1180[k] = f_12 * osh_716[k]
                    + f_3 * pc_y[k] * qsh_884[k];

        t_1181[k] = f_14 * osh_887[k]
                    + f_8 * qsg0_635[k]
                    - f_9 * qsg1_635[k]
                    + f_3 * pc_x[k] * qsh_887[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pc_x, pc_y, pc_z, osh_696, osh_719, osh_888, \
                         qsg0_636, qsg1_636, qsh_885, qsh_887, \
                         qsh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_14 * osh_888[k]
                    + f_6 * qsg0_636[k]
                    - f_7 * qsg1_636[k]
                    + f_3 * pc_x[k] * qsh_888[k];

        t_1183[k] = f_23 * osh_696[k]
                    + f_3 * pc_z[k] * qsh_885[k];

        t_1184[k] = f_12 * osh_719[k]
                    + f_3 * pc_y[k] * qsh_887[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_x, pc_z, osh_699, osh_891, osh_892, \
                         qsg0_639, qsg0_640, qsg1_639, qsg1_640, qsh_888, qsh_891, \
                         qsh_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_14 * osh_891[k]
                    + f_6 * qsg0_639[k]
                    - f_7 * qsg1_639[k]
                    + f_3 * pc_x[k] * qsh_891[k];

        t_1186[k] = f_14 * osh_892[k]
                    + f_4 * qsg0_640[k]
                    - f_5 * qsg1_640[k]
                    + f_3 * pc_x[k] * qsh_892[k];

        t_1187[k] = f_23 * osh_699[k]
                    + f_3 * pc_z[k] * qsh_888[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, pc_x, pc_y, osh_723, osh_894, osh_896, \
                         qsg0_642, qsg0_644, qsg1_642, qsg1_644, qsh_891, qsh_894, \
                         qsh_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_14 * osh_894[k]
                    + f_4 * qsg0_642[k]
                    - f_5 * qsg1_642[k]
                    + f_3 * pc_x[k] * qsh_894[k];

        t_1189[k] = f_12 * osh_723[k]
                    + f_3 * pc_y[k] * qsh_891[k];

        t_1190[k] = f_14 * osh_896[k]
                    + f_4 * qsg0_644[k]
                    - f_5 * qsg1_644[k]
                    + f_3 * pc_x[k] * qsh_896[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, t_1194, t_1195, pc_x, osh_897, osh_898, \
                         osh_899, osh_900, osh_901, qsh_897, qsh_898, qsh_899, qsh_900, \
                         qsh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = f_14 * osh_897[k]
                    + f_3 * pc_x[k] * qsh_897[k];

        t_1192[k] = f_14 * osh_898[k]
                    + f_3 * pc_x[k] * qsh_898[k];

        t_1193[k] = f_14 * osh_899[k]
                    + f_3 * pc_x[k] * qsh_899[k];

        t_1194[k] = f_14 * osh_900[k]
                    + f_3 * pc_x[k] * qsh_900[k];

        t_1195[k] = f_14 * osh_901[k]
                    + f_3 * pc_x[k] * qsh_901[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, pc_z, osh_708, osh_729, osh_902, \
                         qsg0_640, qsg1_640, qsh_897, qsh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_14 * osh_902[k]
                    + f_3 * pc_x[k] * qsh_902[k];

        t_1197[k] = f_12 * osh_729[k]
                    + f_1 * qsg0_640[k]
                    - f_2 * qsg1_640[k]
                    + f_3 * pc_y[k] * qsh_897[k];

        t_1198[k] = f_23 * osh_708[k]
                    + f_3 * pc_z[k] * qsh_897[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_y, osh_731, osh_732, osh_733, qsg0_642, \
                         qsg0_643, qsg0_644, qsg1_642, qsg1_643, qsg1_644, qsh_899, qsh_900, \
                         qsh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_12 * osh_731[k]
                    + f_8 * qsg0_642[k]
                    - f_9 * qsg1_642[k]
                    + f_3 * pc_y[k] * qsh_899[k];

        t_1200[k] = f_12 * osh_732[k]
                    + f_6 * qsg0_643[k]
                    - f_7 * qsg1_643[k]
                    + f_3 * pc_y[k] * qsh_900[k];

        t_1201[k] = f_12 * osh_733[k]
                    + f_4 * qsg0_644[k]
                    - f_5 * qsg1_644[k]
                    + f_3 * pc_y[k] * qsh_901[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, pa_y, pc_y, pc_z, osi0_980, osh_713, \
                         osh_734, osh_735, osi1_980, qsg0_644, qsg1_644, qsh_902, \
                         qsh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_12 * osh_734[k]
                    + f_3 * pc_y[k] * qsh_902[k];

        t_1203[k] = f_23 * osh_713[k]
                    + f_1 * qsg0_644[k]
                    - f_2 * qsg1_644[k]
                    + f_3 * pc_z[k] * qsh_902[k];

        t_1204[k] = pa_y[k] * osi0_980[k]
                    - f_10 * pc_y[k] * osi1_980[k];

        t_1205[k] = f_11 * osh_735[k]
                    + f_3 * pc_y[k] * qsh_903[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pa_y, pc_y, pc_z, osi0_983, osi0_985, \
                         osh_714, osh_736, osh_737, osi1_983, osi1_985, qsh_903, \
                         qsh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_21 * osh_714[k]
                    + f_3 * pc_z[k] * qsh_903[k];

        t_1207[k] = pa_y[k] * osi0_983[k]
                    + f_12 * osh_736[k]
                    - f_10 * pc_y[k] * osi1_983[k];

        t_1208[k] = f_11 * osh_737[k]
                    + f_3 * pc_y[k] * qsh_905[k];

        t_1209[k] = pa_y[k] * osi0_985[k]
                    - f_10 * pc_y[k] * osi1_985[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pa_y, pc_y, pc_z, osi0_986, osi0_989, \
                         osh_717, osh_738, osh_740, osi1_986, osi1_989, qsh_906, \
                         qsh_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = pa_y[k] * osi0_986[k]
                    + f_13 * osh_738[k]
                    - f_10 * pc_y[k] * osi1_986[k];

        t_1211[k] = f_21 * osh_717[k]
                    + f_3 * pc_z[k] * qsh_906[k];

        t_1212[k] = f_11 * osh_740[k]
                    + f_3 * pc_y[k] * qsh_908[k];

        t_1213[k] = pa_y[k] * osi0_989[k]
                    - f_10 * pc_y[k] * osi1_989[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, pa_y, pc_y, pc_z, osi0_990, osi0_992, \
                         osh_720, osh_741, osh_743, osi1_990, osi1_992, \
                         qsh_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pa_y[k] * osi0_990[k]
                    + f_14 * osh_741[k]
                    - f_10 * pc_y[k] * osi1_990[k];

        t_1215[k] = f_21 * osh_720[k]
                    + f_3 * pc_z[k] * qsh_909[k];

        t_1216[k] = pa_y[k] * osi0_992[k]
                    + f_12 * osh_743[k]
                    - f_10 * pc_y[k] * osi1_992[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, pa_y, pc_x, pc_y, osi0_994, osh_744, \
                         osh_918, osh_919, osi1_994, qsh_912, qsh_918, \
                         qsh_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_11 * osh_744[k]
                    + f_3 * pc_y[k] * qsh_912[k];

        t_1218[k] = pa_y[k] * osi0_994[k]
                    - f_10 * pc_y[k] * osi1_994[k];

        t_1219[k] = f_14 * osh_918[k]
                    + f_3 * pc_x[k] * qsh_918[k];

        t_1220[k] = f_14 * osh_919[k]
                    + f_3 * pc_x[k] * qsh_919[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, t_1224, pc_x, osh_920, osh_921, osh_922, \
                         osh_923, qsh_920, qsh_921, qsh_922, qsh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_14 * osh_920[k]
                    + f_3 * pc_x[k] * qsh_920[k];

        t_1222[k] = f_14 * osh_921[k]
                    + f_3 * pc_x[k] * qsh_921[k];

        t_1223[k] = f_14 * osh_922[k]
                    + f_3 * pc_x[k] * qsh_922[k];

        t_1224[k] = f_14 * osh_923[k]
                    + f_3 * pc_x[k] * qsh_923[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, pc_y, pc_z, osh_729, osh_750, osh_752, \
                         qsg0_655, qsg0_657, qsg1_655, qsg1_657, qsh_918, \
                         qsh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_11 * osh_750[k]
                    + f_1 * qsg0_655[k]
                    - f_2 * qsg1_655[k]
                    + f_3 * pc_y[k] * qsh_918[k];

        t_1226[k] = f_21 * osh_729[k]
                    + f_3 * pc_z[k] * qsh_918[k];

        t_1227[k] = f_11 * osh_752[k]
                    + f_8 * qsg0_657[k]
                    - f_9 * qsg1_657[k]
                    + f_3 * pc_y[k] * qsh_920[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, pc_y, osh_753, osh_754, osh_755, qsg0_658, \
                         qsg0_659, qsg1_658, qsg1_659, qsh_921, qsh_922, \
                         qsh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = f_11 * osh_753[k]
                    + f_6 * qsg0_658[k]
                    - f_7 * qsg1_658[k]
                    + f_3 * pc_y[k] * qsh_921[k];

        t_1229[k] = f_11 * osh_754[k]
                    + f_4 * qsg0_659[k]
                    - f_5 * qsg1_659[k]
                    + f_3 * pc_y[k] * qsh_922[k];

        t_1230[k] = f_11 * osh_755[k]
                    + f_3 * pc_y[k] * qsh_923[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pa_y, pc_x, pc_y, pc_z, osi0_1007, \
                         osh_735, osh_924, osi1_1007, qsg0_660, qsg1_660, \
                         qsh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = pa_y[k] * osi0_1007[k]
                    - f_10 * pc_y[k] * osi1_1007[k];

        t_1232[k] = f_14 * osh_924[k]
                    + f_1 * qsg0_660[k]
                    - f_2 * qsg1_660[k]
                    + f_3 * pc_x[k] * qsh_924[k];

        t_1233[k] = f_3 * pc_y[k] * qsh_924[k];

        t_1234[k] = f_20 * osh_735[k]
                    + f_3 * pc_z[k] * qsh_924[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pc_x, pc_y, osh_929, qsg0_660, qsg0_665, \
                         qsg1_660, qsg1_665, qsh_925, qsh_926, \
                         qsh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_4 * qsg0_660[k]
                    - f_5 * qsg1_660[k]
                    + f_3 * pc_y[k] * qsh_925[k];

        t_1236[k] = f_3 * pc_y[k] * qsh_926[k];

        t_1237[k] = f_14 * osh_929[k]
                    + f_8 * qsg0_665[k]
                    - f_9 * qsg1_665[k]
                    + f_3 * pc_x[k] * qsh_929[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, pc_y, qsg0_661, qsg0_662, qsg1_661, qsg1_662, \
                         qsh_927, qsh_928, qsh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_6 * qsg0_661[k]
                    - f_7 * qsg1_661[k]
                    + f_3 * pc_y[k] * qsh_927[k];

        t_1239[k] = f_4 * qsg0_662[k]
                    - f_5 * qsg1_662[k]
                    + f_3 * pc_y[k] * qsh_928[k];

        t_1240[k] = f_3 * pc_y[k] * qsh_929[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pc_y, osh_933, qsg0_663, qsg0_664, \
                         qsg0_669, qsg1_663, qsg1_664, qsg1_669, qsh_930, qsh_931, \
                         qsh_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_14 * osh_933[k]
                    + f_6 * qsg0_669[k]
                    - f_7 * qsg1_669[k]
                    + f_3 * pc_x[k] * qsh_933[k];

        t_1242[k] = f_8 * qsg0_663[k]
                    - f_9 * qsg1_663[k]
                    + f_3 * pc_y[k] * qsh_930[k];

        t_1243[k] = f_6 * qsg0_664[k]
                    - f_7 * qsg1_664[k]
                    + f_3 * pc_y[k] * qsh_931[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, pc_x, pc_y, osh_938, osh_939, \
                         qsg0_665, qsg0_674, qsg1_665, qsg1_674, qsh_932, qsh_933, qsh_938, \
                         qsh_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_4 * qsg0_665[k]
                    - f_5 * qsg1_665[k]
                    + f_3 * pc_y[k] * qsh_932[k];

        t_1245[k] = f_3 * pc_y[k] * qsh_933[k];

        t_1246[k] = f_14 * osh_938[k]
                    + f_4 * qsg0_674[k]
                    - f_5 * qsg1_674[k]
                    + f_3 * pc_x[k] * qsh_938[k];

        t_1247[k] = f_14 * osh_939[k]
                    + f_3 * pc_x[k] * qsh_939[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, t_1252, pc_x, pc_y, osh_940, osh_941, \
                         osh_942, osh_944, qsh_938, qsh_940, qsh_941, qsh_942, \
                         qsh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = f_14 * osh_940[k]
                    + f_3 * pc_x[k] * qsh_940[k];

        t_1249[k] = f_14 * osh_941[k]
                    + f_3 * pc_x[k] * qsh_941[k];

        t_1250[k] = f_14 * osh_942[k]
                    + f_3 * pc_x[k] * qsh_942[k];

        t_1251[k] = f_3 * pc_y[k] * qsh_938[k];

        t_1252[k] = f_14 * osh_944[k]
                    + f_3 * pc_x[k] * qsh_944[k];
    }

#pragma omp simd aligned(t_1253, t_1254, t_1255, pc_y, qsg0_670, qsg0_671, qsg0_672, qsg1_670, \
                         qsg1_671, qsg1_672, qsh_939, qsh_940, \
                         qsh_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1253[k] = f_1 * qsg0_670[k]
                    - f_2 * qsg1_670[k]
                    + f_3 * pc_y[k] * qsh_939[k];

        t_1254[k] = f_16 * qsg0_671[k]
                    - f_17 * qsg1_671[k]
                    + f_3 * pc_y[k] * qsh_940[k];

        t_1255[k] = f_8 * qsg0_672[k]
                    - f_9 * qsg1_672[k]
                    + f_3 * pc_y[k] * qsh_941[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pc_y, pc_z, osh_755, qsg0_673, \
                         qsg0_674, qsg1_673, qsg1_674, qsh_942, qsh_943, \
                         qsh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_6 * qsg0_673[k]
                    - f_7 * qsg1_673[k]
                    + f_3 * pc_y[k] * qsh_942[k];

        t_1257[k] = f_4 * qsg0_674[k]
                    - f_5 * qsg1_674[k]
                    + f_3 * pc_y[k] * qsh_943[k];

        t_1258[k] = f_3 * pc_y[k] * qsh_944[k];

        t_1259[k] = f_20 * osh_755[k]
                    + f_1 * qsg0_674[k]
                    - f_2 * qsg1_674[k]
                    + f_3 * pc_z[k] * qsh_944[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsg0, const size_t qsg1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_19 = 4.5 / q;
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);
    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);
    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);
    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);
    auto *t_1296 = buffer.data(target + 1296);
    auto *t_1297 = buffer.data(target + 1297);
    auto *t_1298 = buffer.data(target + 1298);
    auto *t_1299 = buffer.data(target + 1299);
    auto *t_1300 = buffer.data(target + 1300);
    auto *t_1301 = buffer.data(target + 1301);
    auto *t_1302 = buffer.data(target + 1302);
    auto *t_1303 = buffer.data(target + 1303);
    auto *t_1304 = buffer.data(target + 1304);
    auto *t_1305 = buffer.data(target + 1305);
    auto *t_1306 = buffer.data(target + 1306);
    auto *t_1307 = buffer.data(target + 1307);
    auto *t_1308 = buffer.data(target + 1308);
    auto *t_1309 = buffer.data(target + 1309);
    auto *t_1310 = buffer.data(target + 1310);
    auto *t_1311 = buffer.data(target + 1311);
    auto *t_1312 = buffer.data(target + 1312);
    auto *t_1313 = buffer.data(target + 1313);
    auto *t_1314 = buffer.data(target + 1314);
    auto *t_1315 = buffer.data(target + 1315);
    auto *t_1316 = buffer.data(target + 1316);
    auto *t_1317 = buffer.data(target + 1317);
    auto *t_1318 = buffer.data(target + 1318);
    auto *t_1319 = buffer.data(target + 1319);
    auto *t_1320 = buffer.data(target + 1320);
    auto *t_1321 = buffer.data(target + 1321);
    auto *t_1322 = buffer.data(target + 1322);
    auto *t_1323 = buffer.data(target + 1323);
    auto *t_1324 = buffer.data(target + 1324);
    auto *t_1325 = buffer.data(target + 1325);
    auto *t_1326 = buffer.data(target + 1326);
    auto *t_1327 = buffer.data(target + 1327);
    auto *t_1328 = buffer.data(target + 1328);
    auto *t_1329 = buffer.data(target + 1329);
    auto *t_1330 = buffer.data(target + 1330);
    auto *t_1331 = buffer.data(target + 1331);
    auto *t_1332 = buffer.data(target + 1332);
    auto *t_1333 = buffer.data(target + 1333);
    auto *t_1334 = buffer.data(target + 1334);
    auto *t_1335 = buffer.data(target + 1335);
    auto *t_1336 = buffer.data(target + 1336);
    auto *t_1337 = buffer.data(target + 1337);
    auto *t_1338 = buffer.data(target + 1338);
    auto *t_1339 = buffer.data(target + 1339);
    auto *t_1340 = buffer.data(target + 1340);
    auto *t_1341 = buffer.data(target + 1341);
    auto *t_1342 = buffer.data(target + 1342);
    auto *t_1343 = buffer.data(target + 1343);
    auto *t_1344 = buffer.data(target + 1344);
    auto *t_1345 = buffer.data(target + 1345);
    auto *t_1346 = buffer.data(target + 1346);
    auto *t_1347 = buffer.data(target + 1347);
    auto *t_1348 = buffer.data(target + 1348);
    auto *t_1349 = buffer.data(target + 1349);
    auto *t_1350 = buffer.data(target + 1350);
    auto *t_1351 = buffer.data(target + 1351);
    auto *t_1352 = buffer.data(target + 1352);
    auto *t_1353 = buffer.data(target + 1353);
    auto *t_1354 = buffer.data(target + 1354);
    auto *t_1355 = buffer.data(target + 1355);
    auto *t_1356 = buffer.data(target + 1356);
    auto *t_1357 = buffer.data(target + 1357);
    auto *t_1358 = buffer.data(target + 1358);
    auto *t_1359 = buffer.data(target + 1359);
    auto *t_1360 = buffer.data(target + 1360);
    auto *t_1361 = buffer.data(target + 1361);
    auto *t_1362 = buffer.data(target + 1362);
    auto *t_1363 = buffer.data(target + 1363);
    auto *t_1364 = buffer.data(target + 1364);
    auto *t_1365 = buffer.data(target + 1365);
    auto *t_1366 = buffer.data(target + 1366);
    auto *t_1367 = buffer.data(target + 1367);
    auto *t_1368 = buffer.data(target + 1368);
    auto *t_1369 = buffer.data(target + 1369);
    auto *t_1370 = buffer.data(target + 1370);
    auto *t_1371 = buffer.data(target + 1371);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_1008 = buffer.data(osi0 + 1008);
    const auto *osi0_1011 = buffer.data(osi0 + 1011);
    const auto *osi0_1014 = buffer.data(osi0 + 1014);
    const auto *osi0_1018 = buffer.data(osi0 + 1018);
    const auto *osi0_1020 = buffer.data(osi0 + 1020);
    const auto *osi0_1029 = buffer.data(osi0 + 1029);

    const auto *osh_756 = buffer.data(osh + 756);
    const auto *osh_759 = buffer.data(osh + 759);
    const auto *osh_761 = buffer.data(osh + 761);
    const auto *osh_762 = buffer.data(osh + 762);
    const auto *osh_763 = buffer.data(osh + 763);
    const auto *osh_765 = buffer.data(osh + 765);
    const auto *osh_771 = buffer.data(osh + 771);
    const auto *osh_776 = buffer.data(osh + 776);
    const auto *osh_777 = buffer.data(osh + 777);
    const auto *osh_779 = buffer.data(osh + 779);
    const auto *osh_780 = buffer.data(osh + 780);
    const auto *osh_782 = buffer.data(osh + 782);
    const auto *osh_783 = buffer.data(osh + 783);
    const auto *osh_786 = buffer.data(osh + 786);
    const auto *osh_792 = buffer.data(osh + 792);
    const auto *osh_794 = buffer.data(osh + 794);
    const auto *osh_795 = buffer.data(osh + 795);
    const auto *osh_796 = buffer.data(osh + 796);
    const auto *osh_797 = buffer.data(osh + 797);
    const auto *osh_798 = buffer.data(osh + 798);
    const auto *osh_800 = buffer.data(osh + 800);
    const auto *osh_801 = buffer.data(osh + 801);
    const auto *osh_803 = buffer.data(osh + 803);
    const auto *osh_804 = buffer.data(osh + 804);
    const auto *osh_807 = buffer.data(osh + 807);
    const auto *osh_813 = buffer.data(osh + 813);
    const auto *osh_815 = buffer.data(osh + 815);
    const auto *osh_816 = buffer.data(osh + 816);
    const auto *osh_817 = buffer.data(osh + 817);
    const auto *osh_818 = buffer.data(osh + 818);
    const auto *osh_819 = buffer.data(osh + 819);
    const auto *osh_821 = buffer.data(osh + 821);
    const auto *osh_824 = buffer.data(osh + 824);
    const auto *osh_828 = buffer.data(osh + 828);
    const auto *osh_834 = buffer.data(osh + 834);
    const auto *osh_836 = buffer.data(osh + 836);
    const auto *osh_837 = buffer.data(osh + 837);
    const auto *osh_838 = buffer.data(osh + 838);
    const auto *osh_839 = buffer.data(osh + 839);
    const auto *osh_945 = buffer.data(osh + 945);
    const auto *osh_948 = buffer.data(osh + 948);
    const auto *osh_951 = buffer.data(osh + 951);
    const auto *osh_955 = buffer.data(osh + 955);
    const auto *osh_960 = buffer.data(osh + 960);
    const auto *osh_962 = buffer.data(osh + 962);
    const auto *osh_963 = buffer.data(osh + 963);
    const auto *osh_964 = buffer.data(osh + 964);
    const auto *osh_965 = buffer.data(osh + 965);
    const auto *osh_971 = buffer.data(osh + 971);
    const auto *osh_975 = buffer.data(osh + 975);
    const auto *osh_980 = buffer.data(osh + 980);
    const auto *osh_981 = buffer.data(osh + 981);
    const auto *osh_982 = buffer.data(osh + 982);
    const auto *osh_983 = buffer.data(osh + 983);
    const auto *osh_984 = buffer.data(osh + 984);
    const auto *osh_985 = buffer.data(osh + 985);
    const auto *osh_986 = buffer.data(osh + 986);
    const auto *osh_987 = buffer.data(osh + 987);
    const auto *osh_990 = buffer.data(osh + 990);
    const auto *osh_992 = buffer.data(osh + 992);
    const auto *osh_993 = buffer.data(osh + 993);
    const auto *osh_996 = buffer.data(osh + 996);
    const auto *osh_997 = buffer.data(osh + 997);
    const auto *osh_999 = buffer.data(osh + 999);
    const auto *osh_1001 = buffer.data(osh + 1001);
    const auto *osh_1002 = buffer.data(osh + 1002);
    const auto *osh_1003 = buffer.data(osh + 1003);
    const auto *osh_1004 = buffer.data(osh + 1004);
    const auto *osh_1005 = buffer.data(osh + 1005);
    const auto *osh_1006 = buffer.data(osh + 1006);
    const auto *osh_1007 = buffer.data(osh + 1007);
    const auto *osh_1008 = buffer.data(osh + 1008);
    const auto *osh_1011 = buffer.data(osh + 1011);
    const auto *osh_1013 = buffer.data(osh + 1013);
    const auto *osh_1014 = buffer.data(osh + 1014);
    const auto *osh_1017 = buffer.data(osh + 1017);
    const auto *osh_1018 = buffer.data(osh + 1018);
    const auto *osh_1020 = buffer.data(osh + 1020);
    const auto *osh_1022 = buffer.data(osh + 1022);
    const auto *osh_1023 = buffer.data(osh + 1023);
    const auto *osh_1024 = buffer.data(osh + 1024);
    const auto *osh_1025 = buffer.data(osh + 1025);
    const auto *osh_1026 = buffer.data(osh + 1026);
    const auto *osh_1027 = buffer.data(osh + 1027);
    const auto *osh_1028 = buffer.data(osh + 1028);

    const auto *osi1_1008 = buffer.data(osi1 + 1008);
    const auto *osi1_1011 = buffer.data(osi1 + 1011);
    const auto *osi1_1014 = buffer.data(osi1 + 1014);
    const auto *osi1_1018 = buffer.data(osi1 + 1018);
    const auto *osi1_1020 = buffer.data(osi1 + 1020);
    const auto *osi1_1029 = buffer.data(osi1 + 1029);

    const auto *qsg0_675 = buffer.data(qsg0 + 675);
    const auto *qsg0_677 = buffer.data(qsg0 + 677);
    const auto *qsg0_678 = buffer.data(qsg0 + 678);
    const auto *qsg0_680 = buffer.data(qsg0 + 680);
    const auto *qsg0_681 = buffer.data(qsg0 + 681);
    const auto *qsg0_685 = buffer.data(qsg0 + 685);
    const auto *qsg0_686 = buffer.data(qsg0 + 686);
    const auto *qsg0_687 = buffer.data(qsg0 + 687);
    const auto *qsg0_689 = buffer.data(qsg0 + 689);
    const auto *qsg0_695 = buffer.data(qsg0 + 695);
    const auto *qsg0_699 = buffer.data(qsg0 + 699);
    const auto *qsg0_702 = buffer.data(qsg0 + 702);
    const auto *qsg0_703 = buffer.data(qsg0 + 703);
    const auto *qsg0_704 = buffer.data(qsg0 + 704);
    const auto *qsg0_705 = buffer.data(qsg0 + 705);
    const auto *qsg0_708 = buffer.data(qsg0 + 708);
    const auto *qsg0_710 = buffer.data(qsg0 + 710);
    const auto *qsg0_711 = buffer.data(qsg0 + 711);
    const auto *qsg0_714 = buffer.data(qsg0 + 714);
    const auto *qsg0_715 = buffer.data(qsg0 + 715);
    const auto *qsg0_717 = buffer.data(qsg0 + 717);
    const auto *qsg0_718 = buffer.data(qsg0 + 718);
    const auto *qsg0_719 = buffer.data(qsg0 + 719);
    const auto *qsg0_720 = buffer.data(qsg0 + 720);
    const auto *qsg0_723 = buffer.data(qsg0 + 723);
    const auto *qsg0_725 = buffer.data(qsg0 + 725);
    const auto *qsg0_726 = buffer.data(qsg0 + 726);
    const auto *qsg0_729 = buffer.data(qsg0 + 729);
    const auto *qsg0_730 = buffer.data(qsg0 + 730);
    const auto *qsg0_732 = buffer.data(qsg0 + 732);
    const auto *qsg0_733 = buffer.data(qsg0 + 733);
    const auto *qsg0_734 = buffer.data(qsg0 + 734);

    const auto *qsg1_675 = buffer.data(qsg1 + 675);
    const auto *qsg1_677 = buffer.data(qsg1 + 677);
    const auto *qsg1_678 = buffer.data(qsg1 + 678);
    const auto *qsg1_680 = buffer.data(qsg1 + 680);
    const auto *qsg1_681 = buffer.data(qsg1 + 681);
    const auto *qsg1_685 = buffer.data(qsg1 + 685);
    const auto *qsg1_686 = buffer.data(qsg1 + 686);
    const auto *qsg1_687 = buffer.data(qsg1 + 687);
    const auto *qsg1_689 = buffer.data(qsg1 + 689);
    const auto *qsg1_695 = buffer.data(qsg1 + 695);
    const auto *qsg1_699 = buffer.data(qsg1 + 699);
    const auto *qsg1_702 = buffer.data(qsg1 + 702);
    const auto *qsg1_703 = buffer.data(qsg1 + 703);
    const auto *qsg1_704 = buffer.data(qsg1 + 704);
    const auto *qsg1_705 = buffer.data(qsg1 + 705);
    const auto *qsg1_708 = buffer.data(qsg1 + 708);
    const auto *qsg1_710 = buffer.data(qsg1 + 710);
    const auto *qsg1_711 = buffer.data(qsg1 + 711);
    const auto *qsg1_714 = buffer.data(qsg1 + 714);
    const auto *qsg1_715 = buffer.data(qsg1 + 715);
    const auto *qsg1_717 = buffer.data(qsg1 + 717);
    const auto *qsg1_718 = buffer.data(qsg1 + 718);
    const auto *qsg1_719 = buffer.data(qsg1 + 719);
    const auto *qsg1_720 = buffer.data(qsg1 + 720);
    const auto *qsg1_723 = buffer.data(qsg1 + 723);
    const auto *qsg1_725 = buffer.data(qsg1 + 725);
    const auto *qsg1_726 = buffer.data(qsg1 + 726);
    const auto *qsg1_729 = buffer.data(qsg1 + 729);
    const auto *qsg1_730 = buffer.data(qsg1 + 730);
    const auto *qsg1_732 = buffer.data(qsg1 + 732);
    const auto *qsg1_733 = buffer.data(qsg1 + 733);
    const auto *qsg1_734 = buffer.data(qsg1 + 734);

    const auto *qsh_945 = buffer.data(qsh + 945);
    const auto *qsh_946 = buffer.data(qsh + 946);
    const auto *qsh_947 = buffer.data(qsh + 947);
    const auto *qsh_948 = buffer.data(qsh + 948);
    const auto *qsh_950 = buffer.data(qsh + 950);
    const auto *qsh_951 = buffer.data(qsh + 951);
    const auto *qsh_952 = buffer.data(qsh + 952);
    const auto *qsh_954 = buffer.data(qsh + 954);
    const auto *qsh_955 = buffer.data(qsh + 955);
    const auto *qsh_960 = buffer.data(qsh + 960);
    const auto *qsh_961 = buffer.data(qsh + 961);
    const auto *qsh_962 = buffer.data(qsh + 962);
    const auto *qsh_963 = buffer.data(qsh + 963);
    const auto *qsh_964 = buffer.data(qsh + 964);
    const auto *qsh_965 = buffer.data(qsh + 965);
    const auto *qsh_966 = buffer.data(qsh + 966);
    const auto *qsh_968 = buffer.data(qsh + 968);
    const auto *qsh_969 = buffer.data(qsh + 969);
    const auto *qsh_971 = buffer.data(qsh + 971);
    const auto *qsh_972 = buffer.data(qsh + 972);
    const auto *qsh_975 = buffer.data(qsh + 975);
    const auto *qsh_980 = buffer.data(qsh + 980);
    const auto *qsh_981 = buffer.data(qsh + 981);
    const auto *qsh_982 = buffer.data(qsh + 982);
    const auto *qsh_983 = buffer.data(qsh + 983);
    const auto *qsh_984 = buffer.data(qsh + 984);
    const auto *qsh_985 = buffer.data(qsh + 985);
    const auto *qsh_986 = buffer.data(qsh + 986);
    const auto *qsh_987 = buffer.data(qsh + 987);
    const auto *qsh_989 = buffer.data(qsh + 989);
    const auto *qsh_990 = buffer.data(qsh + 990);
    const auto *qsh_992 = buffer.data(qsh + 992);
    const auto *qsh_993 = buffer.data(qsh + 993);
    const auto *qsh_996 = buffer.data(qsh + 996);
    const auto *qsh_997 = buffer.data(qsh + 997);
    const auto *qsh_999 = buffer.data(qsh + 999);
    const auto *qsh_1001 = buffer.data(qsh + 1001);
    const auto *qsh_1002 = buffer.data(qsh + 1002);
    const auto *qsh_1003 = buffer.data(qsh + 1003);
    const auto *qsh_1004 = buffer.data(qsh + 1004);
    const auto *qsh_1005 = buffer.data(qsh + 1005);
    const auto *qsh_1006 = buffer.data(qsh + 1006);
    const auto *qsh_1007 = buffer.data(qsh + 1007);
    const auto *qsh_1008 = buffer.data(qsh + 1008);
    const auto *qsh_1010 = buffer.data(qsh + 1010);
    const auto *qsh_1011 = buffer.data(qsh + 1011);
    const auto *qsh_1013 = buffer.data(qsh + 1013);
    const auto *qsh_1014 = buffer.data(qsh + 1014);
    const auto *qsh_1017 = buffer.data(qsh + 1017);
    const auto *qsh_1018 = buffer.data(qsh + 1018);
    const auto *qsh_1020 = buffer.data(qsh + 1020);
    const auto *qsh_1022 = buffer.data(qsh + 1022);
    const auto *qsh_1023 = buffer.data(qsh + 1023);
    const auto *qsh_1024 = buffer.data(qsh + 1024);
    const auto *qsh_1025 = buffer.data(qsh + 1025);
    const auto *qsh_1026 = buffer.data(qsh + 1026);
    const auto *qsh_1027 = buffer.data(qsh + 1027);
    const auto *qsh_1028 = buffer.data(qsh + 1028);

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pc_x, pc_y, pc_z, osh_756, osh_945, \
                         osh_948, qsg0_675, qsg0_678, qsg1_675, qsg1_678, qsh_945, \
                         qsh_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_13 * osh_945[k]
                    + f_1 * qsg0_675[k]
                    - f_2 * qsg1_675[k]
                    + f_3 * pc_x[k] * qsh_945[k];

        t_1261[k] = f_19 * osh_756[k]
                    + f_3 * pc_y[k] * qsh_945[k];

        t_1262[k] = f_3 * pc_z[k] * qsh_945[k];

        t_1263[k] = f_13 * osh_948[k]
                    + f_8 * qsg0_678[k]
                    - f_9 * qsg1_678[k]
                    + f_3 * pc_x[k] * qsh_948[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, t_1267, pc_x, pc_z, osh_951, qsg0_675, \
                         qsg0_681, qsg1_675, qsg1_681, qsh_946, qsh_947, qsh_948, \
                         qsh_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_3 * pc_z[k] * qsh_946[k];

        t_1265[k] = f_4 * qsg0_675[k]
                    - f_5 * qsg1_675[k]
                    + f_3 * pc_z[k] * qsh_947[k];

        t_1266[k] = f_13 * osh_951[k]
                    + f_6 * qsg0_681[k]
                    - f_7 * qsg1_681[k]
                    + f_3 * pc_x[k] * qsh_951[k];

        t_1267[k] = f_3 * pc_z[k] * qsh_948[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, t_1271, pc_x, pc_y, pc_z, osh_761, osh_955, \
                         qsg0_677, qsg0_685, qsg1_677, qsg1_685, qsh_950, qsh_951, \
                         qsh_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_19 * osh_761[k]
                    + f_3 * pc_y[k] * qsh_950[k];

        t_1269[k] = f_6 * qsg0_677[k]
                    - f_7 * qsg1_677[k]
                    + f_3 * pc_z[k] * qsh_950[k];

        t_1270[k] = f_13 * osh_955[k]
                    + f_4 * qsg0_685[k]
                    - f_5 * qsg1_685[k]
                    + f_3 * pc_x[k] * qsh_955[k];

        t_1271[k] = f_3 * pc_z[k] * qsh_951[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, t_1275, pc_x, pc_y, pc_z, osh_765, osh_960, \
                         qsg0_678, qsg0_680, qsg1_678, qsg1_680, qsh_952, qsh_954, \
                         qsh_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_4 * qsg0_678[k]
                    - f_5 * qsg1_678[k]
                    + f_3 * pc_z[k] * qsh_952[k];

        t_1273[k] = f_19 * osh_765[k]
                    + f_3 * pc_y[k] * qsh_954[k];

        t_1274[k] = f_8 * qsg0_680[k]
                    - f_9 * qsg1_680[k]
                    + f_3 * pc_z[k] * qsh_954[k];

        t_1275[k] = f_13 * osh_960[k]
                    + f_3 * pc_x[k] * qsh_960[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, t_1279, t_1280, pc_x, pc_z, osh_962, osh_963, \
                         osh_964, osh_965, qsh_955, qsh_962, qsh_963, qsh_964, \
                         qsh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_3 * pc_z[k] * qsh_955[k];

        t_1277[k] = f_13 * osh_962[k]
                    + f_3 * pc_x[k] * qsh_962[k];

        t_1278[k] = f_13 * osh_963[k]
                    + f_3 * pc_x[k] * qsh_963[k];

        t_1279[k] = f_13 * osh_964[k]
                    + f_3 * pc_x[k] * qsh_964[k];

        t_1280[k] = f_13 * osh_965[k]
                    + f_3 * pc_x[k] * qsh_965[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, t_1284, pc_y, pc_z, osh_771, qsg0_685, \
                         qsg0_686, qsg1_685, qsg1_686, qsh_960, qsh_961, \
                         qsh_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_19 * osh_771[k]
                    + f_1 * qsg0_685[k]
                    - f_2 * qsg1_685[k]
                    + f_3 * pc_y[k] * qsh_960[k];

        t_1282[k] = f_3 * pc_z[k] * qsh_960[k];

        t_1283[k] = f_4 * qsg0_685[k]
                    - f_5 * qsg1_685[k]
                    + f_3 * pc_z[k] * qsh_961[k];

        t_1284[k] = f_6 * qsg0_686[k]
                    - f_7 * qsg1_686[k]
                    + f_3 * pc_z[k] * qsh_962[k];
    }

#pragma omp simd aligned(t_1285, t_1286, t_1287, t_1288, pa_z, pc_y, pc_z, osi0_1008, osh_776, \
                         osi1_1008, qsg0_687, qsg0_689, qsg1_687, qsg1_689, qsh_963, \
                         qsh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1285[k] = f_8 * qsg0_687[k]
                    - f_9 * qsg1_687[k]
                    + f_3 * pc_z[k] * qsh_963[k];

        t_1286[k] = f_19 * osh_776[k]
                    + f_3 * pc_y[k] * qsh_965[k];

        t_1287[k] = f_1 * qsg0_689[k]
                    - f_2 * qsg1_689[k]
                    + f_3 * pc_z[k] * qsh_965[k];

        t_1288[k] = pa_z[k] * osi0_1008[k]
                    - f_10 * pc_z[k] * osi1_1008[k];
    }

#pragma omp simd aligned(t_1289, t_1290, t_1291, t_1292, pa_z, pc_y, pc_z, osi0_1011, osh_756, \
                         osh_777, osh_779, osi1_1011, qsh_966, \
                         qsh_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1289[k] = f_20 * osh_777[k]
                    + f_3 * pc_y[k] * qsh_966[k];

        t_1290[k] = f_11 * osh_756[k]
                    + f_3 * pc_z[k] * qsh_966[k];

        t_1291[k] = pa_z[k] * osi0_1011[k]
                    - f_10 * pc_z[k] * osi1_1011[k];

        t_1292[k] = f_20 * osh_779[k]
                    + f_3 * pc_y[k] * qsh_968[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, pa_z, pc_x, pc_z, osi0_1014, osh_759, \
                         osh_971, osi1_1014, qsg0_695, qsg1_695, qsh_969, \
                         qsh_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_13 * osh_971[k]
                    + f_8 * qsg0_695[k]
                    - f_9 * qsg1_695[k]
                    + f_3 * pc_x[k] * qsh_971[k];

        t_1294[k] = pa_z[k] * osi0_1014[k]
                    - f_10 * pc_z[k] * osi1_1014[k];

        t_1295[k] = f_11 * osh_759[k]
                    + f_3 * pc_z[k] * qsh_969[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, pa_z, pc_x, pc_y, pc_z, osi0_1018, osh_782, \
                         osh_975, osi1_1018, qsg0_699, qsg1_699, qsh_971, \
                         qsh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_20 * osh_782[k]
                    + f_3 * pc_y[k] * qsh_971[k];

        t_1297[k] = f_13 * osh_975[k]
                    + f_6 * qsg0_699[k]
                    - f_7 * qsg1_699[k]
                    + f_3 * pc_x[k] * qsh_975[k];

        t_1298[k] = pa_z[k] * osi0_1018[k]
                    - f_10 * pc_z[k] * osi1_1018[k];
    }

#pragma omp simd aligned(t_1299, t_1300, t_1301, pa_z, pc_y, pc_z, osi0_1020, osh_762, \
                         osh_763, osh_786, osi1_1020, qsh_972, \
                         qsh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1299[k] = f_11 * osh_762[k]
                    + f_3 * pc_z[k] * qsh_972[k];

        t_1300[k] = pa_z[k] * osi0_1020[k]
                    + f_12 * osh_763[k]
                    - f_10 * pc_z[k] * osi1_1020[k];

        t_1301[k] = f_20 * osh_786[k]
                    + f_3 * pc_y[k] * qsh_975[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, t_1305, pc_x, osh_980, osh_981, osh_982, \
                         osh_983, qsg0_704, qsg1_704, qsh_980, qsh_981, qsh_982, \
                         qsh_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = f_13 * osh_980[k]
                    + f_4 * qsg0_704[k]
                    - f_5 * qsg1_704[k]
                    + f_3 * pc_x[k] * qsh_980[k];

        t_1303[k] = f_13 * osh_981[k]
                    + f_3 * pc_x[k] * qsh_981[k];

        t_1304[k] = f_13 * osh_982[k]
                    + f_3 * pc_x[k] * qsh_982[k];

        t_1305[k] = f_13 * osh_983[k]
                    + f_3 * pc_x[k] * qsh_983[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, t_1309, pa_z, pc_x, pc_z, osi0_1029, osh_984, \
                         osh_985, osh_986, osi1_1029, qsh_984, qsh_985, \
                         qsh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_13 * osh_984[k]
                    + f_3 * pc_x[k] * qsh_984[k];

        t_1307[k] = f_13 * osh_985[k]
                    + f_3 * pc_x[k] * qsh_985[k];

        t_1308[k] = f_13 * osh_986[k]
                    + f_3 * pc_x[k] * qsh_986[k];

        t_1309[k] = pa_z[k] * osi0_1029[k]
                    - f_10 * pc_z[k] * osi1_1029[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, pc_y, pc_z, osh_771, osh_794, osh_795, \
                         qsg0_702, qsg0_703, qsg1_702, qsg1_703, qsh_981, qsh_983, \
                         qsh_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_11 * osh_771[k]
                    + f_3 * pc_z[k] * qsh_981[k];

        t_1311[k] = f_20 * osh_794[k]
                    + f_8 * qsg0_702[k]
                    - f_9 * qsg1_702[k]
                    + f_3 * pc_y[k] * qsh_983[k];

        t_1312[k] = f_20 * osh_795[k]
                    + f_6 * qsg0_703[k]
                    - f_7 * qsg1_703[k]
                    + f_3 * pc_y[k] * qsh_984[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, pc_y, pc_z, osh_776, osh_796, osh_797, \
                         qsg0_704, qsg1_704, qsh_985, qsh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_20 * osh_796[k]
                    + f_4 * qsg0_704[k]
                    - f_5 * qsg1_704[k]
                    + f_3 * pc_y[k] * qsh_985[k];

        t_1314[k] = f_20 * osh_797[k]
                    + f_3 * pc_y[k] * qsh_986[k];

        t_1315[k] = f_11 * osh_776[k]
                    + f_1 * qsg0_704[k]
                    - f_2 * qsg1_704[k]
                    + f_3 * pc_z[k] * qsh_986[k];
    }

#pragma omp simd aligned(t_1316, t_1317, t_1318, pc_x, pc_y, pc_z, osh_777, osh_798, osh_987, \
                         qsg0_705, qsg1_705, qsh_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1316[k] = f_13 * osh_987[k]
                    + f_1 * qsg0_705[k]
                    - f_2 * qsg1_705[k]
                    + f_3 * pc_x[k] * qsh_987[k];

        t_1317[k] = f_21 * osh_798[k]
                    + f_3 * pc_y[k] * qsh_987[k];

        t_1318[k] = f_12 * osh_777[k]
                    + f_3 * pc_z[k] * qsh_987[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, pc_x, pc_y, osh_800, osh_990, osh_992, \
                         qsg0_708, qsg0_710, qsg1_708, qsg1_710, qsh_989, qsh_990, \
                         qsh_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_13 * osh_990[k]
                    + f_8 * qsg0_708[k]
                    - f_9 * qsg1_708[k]
                    + f_3 * pc_x[k] * qsh_990[k];

        t_1320[k] = f_21 * osh_800[k]
                    + f_3 * pc_y[k] * qsh_989[k];

        t_1321[k] = f_13 * osh_992[k]
                    + f_8 * qsg0_710[k]
                    - f_9 * qsg1_710[k]
                    + f_3 * pc_x[k] * qsh_992[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, pc_x, pc_y, pc_z, osh_780, osh_803, osh_993, \
                         qsg0_711, qsg1_711, qsh_990, qsh_992, \
                         qsh_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = f_13 * osh_993[k]
                    + f_6 * qsg0_711[k]
                    - f_7 * qsg1_711[k]
                    + f_3 * pc_x[k] * qsh_993[k];

        t_1323[k] = f_12 * osh_780[k]
                    + f_3 * pc_z[k] * qsh_990[k];

        t_1324[k] = f_21 * osh_803[k]
                    + f_3 * pc_y[k] * qsh_992[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pc_x, pc_z, osh_783, osh_996, osh_997, \
                         qsg0_714, qsg0_715, qsg1_714, qsg1_715, qsh_993, qsh_996, \
                         qsh_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_13 * osh_996[k]
                    + f_6 * qsg0_714[k]
                    - f_7 * qsg1_714[k]
                    + f_3 * pc_x[k] * qsh_996[k];

        t_1326[k] = f_13 * osh_997[k]
                    + f_4 * qsg0_715[k]
                    - f_5 * qsg1_715[k]
                    + f_3 * pc_x[k] * qsh_997[k];

        t_1327[k] = f_12 * osh_783[k]
                    + f_3 * pc_z[k] * qsh_993[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pc_x, pc_y, osh_807, osh_999, osh_1001, \
                         qsg0_717, qsg0_719, qsg1_717, qsg1_719, qsh_996, qsh_999, \
                         qsh_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_13 * osh_999[k]
                    + f_4 * qsg0_717[k]
                    - f_5 * qsg1_717[k]
                    + f_3 * pc_x[k] * qsh_999[k];

        t_1329[k] = f_21 * osh_807[k]
                    + f_3 * pc_y[k] * qsh_996[k];

        t_1330[k] = f_13 * osh_1001[k]
                    + f_4 * qsg0_719[k]
                    - f_5 * qsg1_719[k]
                    + f_3 * pc_x[k] * qsh_1001[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, t_1335, pc_x, osh_1002, osh_1003, \
                         osh_1004, osh_1005, osh_1006, qsh_1002, qsh_1003, qsh_1004, qsh_1005, \
                         qsh_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_13 * osh_1002[k]
                    + f_3 * pc_x[k] * qsh_1002[k];

        t_1332[k] = f_13 * osh_1003[k]
                    + f_3 * pc_x[k] * qsh_1003[k];

        t_1333[k] = f_13 * osh_1004[k]
                    + f_3 * pc_x[k] * qsh_1004[k];

        t_1334[k] = f_13 * osh_1005[k]
                    + f_3 * pc_x[k] * qsh_1005[k];

        t_1335[k] = f_13 * osh_1006[k]
                    + f_3 * pc_x[k] * qsh_1006[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, pc_x, pc_y, pc_z, osh_792, osh_813, osh_1007, \
                         qsg0_715, qsg1_715, qsh_1002, qsh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_13 * osh_1007[k]
                    + f_3 * pc_x[k] * qsh_1007[k];

        t_1337[k] = f_21 * osh_813[k]
                    + f_1 * qsg0_715[k]
                    - f_2 * qsg1_715[k]
                    + f_3 * pc_y[k] * qsh_1002[k];

        t_1338[k] = f_12 * osh_792[k]
                    + f_3 * pc_z[k] * qsh_1002[k];
    }

#pragma omp simd aligned(t_1339, t_1340, t_1341, pc_y, osh_815, osh_816, osh_817, qsg0_717, \
                         qsg0_718, qsg0_719, qsg1_717, qsg1_718, qsg1_719, qsh_1004, qsh_1005, \
                         qsh_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1339[k] = f_21 * osh_815[k]
                    + f_8 * qsg0_717[k]
                    - f_9 * qsg1_717[k]
                    + f_3 * pc_y[k] * qsh_1004[k];

        t_1340[k] = f_21 * osh_816[k]
                    + f_6 * qsg0_718[k]
                    - f_7 * qsg1_718[k]
                    + f_3 * pc_y[k] * qsh_1005[k];

        t_1341[k] = f_21 * osh_817[k]
                    + f_4 * qsg0_719[k]
                    - f_5 * qsg1_719[k]
                    + f_3 * pc_y[k] * qsh_1006[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, pc_x, pc_y, pc_z, osh_797, osh_818, osh_1008, \
                         qsg0_719, qsg0_720, qsg1_719, qsg1_720, qsh_1007, \
                         qsh_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_21 * osh_818[k]
                    + f_3 * pc_y[k] * qsh_1007[k];

        t_1343[k] = f_12 * osh_797[k]
                    + f_1 * qsg0_719[k]
                    - f_2 * qsg1_719[k]
                    + f_3 * pc_z[k] * qsh_1007[k];

        t_1344[k] = f_13 * osh_1008[k]
                    + f_1 * qsg0_720[k]
                    - f_2 * qsg1_720[k]
                    + f_3 * pc_x[k] * qsh_1008[k];
    }

#pragma omp simd aligned(t_1345, t_1346, t_1347, t_1348, pc_x, pc_y, pc_z, osh_798, osh_819, \
                         osh_821, osh_1011, qsg0_723, qsg1_723, qsh_1008, qsh_1010, \
                         qsh_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1345[k] = f_23 * osh_819[k]
                    + f_3 * pc_y[k] * qsh_1008[k];

        t_1346[k] = f_13 * osh_798[k]
                    + f_3 * pc_z[k] * qsh_1008[k];

        t_1347[k] = f_13 * osh_1011[k]
                    + f_8 * qsg0_723[k]
                    - f_9 * qsg1_723[k]
                    + f_3 * pc_x[k] * qsh_1011[k];

        t_1348[k] = f_23 * osh_821[k]
                    + f_3 * pc_y[k] * qsh_1010[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pc_x, pc_z, osh_801, osh_1013, osh_1014, \
                         qsg0_725, qsg0_726, qsg1_725, qsg1_726, qsh_1011, qsh_1013, \
                         qsh_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = f_13 * osh_1013[k]
                    + f_8 * qsg0_725[k]
                    - f_9 * qsg1_725[k]
                    + f_3 * pc_x[k] * qsh_1013[k];

        t_1350[k] = f_13 * osh_1014[k]
                    + f_6 * qsg0_726[k]
                    - f_7 * qsg1_726[k]
                    + f_3 * pc_x[k] * qsh_1014[k];

        t_1351[k] = f_13 * osh_801[k]
                    + f_3 * pc_z[k] * qsh_1011[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, pc_x, pc_y, osh_824, osh_1017, osh_1018, \
                         qsg0_729, qsg0_730, qsg1_729, qsg1_730, qsh_1013, qsh_1017, \
                         qsh_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_23 * osh_824[k]
                    + f_3 * pc_y[k] * qsh_1013[k];

        t_1353[k] = f_13 * osh_1017[k]
                    + f_6 * qsg0_729[k]
                    - f_7 * qsg1_729[k]
                    + f_3 * pc_x[k] * qsh_1017[k];

        t_1354[k] = f_13 * osh_1018[k]
                    + f_4 * qsg0_730[k]
                    - f_5 * qsg1_730[k]
                    + f_3 * pc_x[k] * qsh_1018[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pc_x, pc_y, pc_z, osh_804, osh_828, osh_1020, \
                         qsg0_732, qsg1_732, qsh_1014, qsh_1017, \
                         qsh_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_13 * osh_804[k]
                    + f_3 * pc_z[k] * qsh_1014[k];

        t_1356[k] = f_13 * osh_1020[k]
                    + f_4 * qsg0_732[k]
                    - f_5 * qsg1_732[k]
                    + f_3 * pc_x[k] * qsh_1020[k];

        t_1357[k] = f_23 * osh_828[k]
                    + f_3 * pc_y[k] * qsh_1017[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, t_1361, pc_x, osh_1022, osh_1023, osh_1024, \
                         osh_1025, qsg0_734, qsg1_734, qsh_1022, qsh_1023, qsh_1024, \
                         qsh_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_13 * osh_1022[k]
                    + f_4 * qsg0_734[k]
                    - f_5 * qsg1_734[k]
                    + f_3 * pc_x[k] * qsh_1022[k];

        t_1359[k] = f_13 * osh_1023[k]
                    + f_3 * pc_x[k] * qsh_1023[k];

        t_1360[k] = f_13 * osh_1024[k]
                    + f_3 * pc_x[k] * qsh_1024[k];

        t_1361[k] = f_13 * osh_1025[k]
                    + f_3 * pc_x[k] * qsh_1025[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, t_1365, pc_x, pc_y, osh_834, osh_1026, \
                         osh_1027, osh_1028, qsg0_730, qsg1_730, qsh_1023, qsh_1026, qsh_1027, \
                         qsh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = f_13 * osh_1026[k]
                    + f_3 * pc_x[k] * qsh_1026[k];

        t_1363[k] = f_13 * osh_1027[k]
                    + f_3 * pc_x[k] * qsh_1027[k];

        t_1364[k] = f_13 * osh_1028[k]
                    + f_3 * pc_x[k] * qsh_1028[k];

        t_1365[k] = f_23 * osh_834[k]
                    + f_1 * qsg0_730[k]
                    - f_2 * qsg1_730[k]
                    + f_3 * pc_y[k] * qsh_1023[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, pc_y, pc_z, osh_813, osh_836, osh_837, \
                         qsg0_732, qsg0_733, qsg1_732, qsg1_733, qsh_1023, qsh_1025, \
                         qsh_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_13 * osh_813[k]
                    + f_3 * pc_z[k] * qsh_1023[k];

        t_1367[k] = f_23 * osh_836[k]
                    + f_8 * qsg0_732[k]
                    - f_9 * qsg1_732[k]
                    + f_3 * pc_y[k] * qsh_1025[k];

        t_1368[k] = f_23 * osh_837[k]
                    + f_6 * qsg0_733[k]
                    - f_7 * qsg1_733[k]
                    + f_3 * pc_y[k] * qsh_1026[k];
    }

#pragma omp simd aligned(t_1369, t_1370, t_1371, pc_y, pc_z, osh_818, osh_838, osh_839, \
                         qsg0_734, qsg1_734, qsh_1027, qsh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1369[k] = f_23 * osh_838[k]
                    + f_4 * qsg0_734[k]
                    - f_5 * qsg1_734[k]
                    + f_3 * pc_y[k] * qsh_1027[k];

        t_1370[k] = f_23 * osh_839[k]
                    + f_3 * pc_y[k] * qsh_1028[k];

        t_1371[k] = f_13 * osh_818[k]
                    + f_1 * qsg0_734[k]
                    - f_2 * qsg1_734[k]
                    + f_3 * pc_z[k] * qsh_1028[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t osh, const size_t qsg0,
                                                           const size_t qsg1, const size_t qsh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_1372 = buffer.data(target + 1372);
    auto *t_1373 = buffer.data(target + 1373);
    auto *t_1374 = buffer.data(target + 1374);
    auto *t_1375 = buffer.data(target + 1375);
    auto *t_1376 = buffer.data(target + 1376);
    auto *t_1377 = buffer.data(target + 1377);
    auto *t_1378 = buffer.data(target + 1378);
    auto *t_1379 = buffer.data(target + 1379);
    auto *t_1380 = buffer.data(target + 1380);
    auto *t_1381 = buffer.data(target + 1381);
    auto *t_1382 = buffer.data(target + 1382);
    auto *t_1383 = buffer.data(target + 1383);
    auto *t_1384 = buffer.data(target + 1384);
    auto *t_1385 = buffer.data(target + 1385);
    auto *t_1386 = buffer.data(target + 1386);
    auto *t_1387 = buffer.data(target + 1387);
    auto *t_1388 = buffer.data(target + 1388);
    auto *t_1389 = buffer.data(target + 1389);
    auto *t_1390 = buffer.data(target + 1390);
    auto *t_1391 = buffer.data(target + 1391);
    auto *t_1392 = buffer.data(target + 1392);
    auto *t_1393 = buffer.data(target + 1393);
    auto *t_1394 = buffer.data(target + 1394);
    auto *t_1395 = buffer.data(target + 1395);
    auto *t_1396 = buffer.data(target + 1396);
    auto *t_1397 = buffer.data(target + 1397);
    auto *t_1398 = buffer.data(target + 1398);
    auto *t_1399 = buffer.data(target + 1399);
    auto *t_1400 = buffer.data(target + 1400);
    auto *t_1401 = buffer.data(target + 1401);
    auto *t_1402 = buffer.data(target + 1402);
    auto *t_1403 = buffer.data(target + 1403);
    auto *t_1404 = buffer.data(target + 1404);
    auto *t_1405 = buffer.data(target + 1405);
    auto *t_1406 = buffer.data(target + 1406);
    auto *t_1407 = buffer.data(target + 1407);
    auto *t_1408 = buffer.data(target + 1408);
    auto *t_1409 = buffer.data(target + 1409);
    auto *t_1410 = buffer.data(target + 1410);
    auto *t_1411 = buffer.data(target + 1411);
    auto *t_1412 = buffer.data(target + 1412);
    auto *t_1413 = buffer.data(target + 1413);
    auto *t_1414 = buffer.data(target + 1414);
    auto *t_1415 = buffer.data(target + 1415);
    auto *t_1416 = buffer.data(target + 1416);
    auto *t_1417 = buffer.data(target + 1417);
    auto *t_1418 = buffer.data(target + 1418);
    auto *t_1419 = buffer.data(target + 1419);
    auto *t_1420 = buffer.data(target + 1420);
    auto *t_1421 = buffer.data(target + 1421);
    auto *t_1422 = buffer.data(target + 1422);
    auto *t_1423 = buffer.data(target + 1423);
    auto *t_1424 = buffer.data(target + 1424);
    auto *t_1425 = buffer.data(target + 1425);
    auto *t_1426 = buffer.data(target + 1426);
    auto *t_1427 = buffer.data(target + 1427);
    auto *t_1428 = buffer.data(target + 1428);
    auto *t_1429 = buffer.data(target + 1429);
    auto *t_1430 = buffer.data(target + 1430);
    auto *t_1431 = buffer.data(target + 1431);
    auto *t_1432 = buffer.data(target + 1432);
    auto *t_1433 = buffer.data(target + 1433);
    auto *t_1434 = buffer.data(target + 1434);
    auto *t_1435 = buffer.data(target + 1435);
    auto *t_1436 = buffer.data(target + 1436);
    auto *t_1437 = buffer.data(target + 1437);
    auto *t_1438 = buffer.data(target + 1438);
    auto *t_1439 = buffer.data(target + 1439);
    auto *t_1440 = buffer.data(target + 1440);
    auto *t_1441 = buffer.data(target + 1441);
    auto *t_1442 = buffer.data(target + 1442);
    auto *t_1443 = buffer.data(target + 1443);
    auto *t_1444 = buffer.data(target + 1444);
    auto *t_1445 = buffer.data(target + 1445);
    auto *t_1446 = buffer.data(target + 1446);
    auto *t_1447 = buffer.data(target + 1447);
    auto *t_1448 = buffer.data(target + 1448);
    auto *t_1449 = buffer.data(target + 1449);
    auto *t_1450 = buffer.data(target + 1450);
    auto *t_1451 = buffer.data(target + 1451);
    auto *t_1452 = buffer.data(target + 1452);
    auto *t_1453 = buffer.data(target + 1453);
    auto *t_1454 = buffer.data(target + 1454);
    auto *t_1455 = buffer.data(target + 1455);
    auto *t_1456 = buffer.data(target + 1456);
    auto *t_1457 = buffer.data(target + 1457);
    auto *t_1458 = buffer.data(target + 1458);
    auto *t_1459 = buffer.data(target + 1459);
    auto *t_1460 = buffer.data(target + 1460);
    auto *t_1461 = buffer.data(target + 1461);
    auto *t_1462 = buffer.data(target + 1462);
    auto *t_1463 = buffer.data(target + 1463);
    auto *t_1464 = buffer.data(target + 1464);
    auto *t_1465 = buffer.data(target + 1465);
    auto *t_1466 = buffer.data(target + 1466);
    auto *t_1467 = buffer.data(target + 1467);
    auto *t_1468 = buffer.data(target + 1468);
    auto *t_1469 = buffer.data(target + 1469);
    auto *t_1470 = buffer.data(target + 1470);
    auto *t_1471 = buffer.data(target + 1471);
    auto *t_1472 = buffer.data(target + 1472);
    auto *t_1473 = buffer.data(target + 1473);
    auto *t_1474 = buffer.data(target + 1474);
    auto *t_1475 = buffer.data(target + 1475);
    auto *t_1476 = buffer.data(target + 1476);
    auto *t_1477 = buffer.data(target + 1477);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh_819 = buffer.data(osh + 819);
    const auto *osh_822 = buffer.data(osh + 822);
    const auto *osh_825 = buffer.data(osh + 825);
    const auto *osh_834 = buffer.data(osh + 834);
    const auto *osh_839 = buffer.data(osh + 839);
    const auto *osh_840 = buffer.data(osh + 840);
    const auto *osh_842 = buffer.data(osh + 842);
    const auto *osh_843 = buffer.data(osh + 843);
    const auto *osh_845 = buffer.data(osh + 845);
    const auto *osh_846 = buffer.data(osh + 846);
    const auto *osh_849 = buffer.data(osh + 849);
    const auto *osh_855 = buffer.data(osh + 855);
    const auto *osh_857 = buffer.data(osh + 857);
    const auto *osh_858 = buffer.data(osh + 858);
    const auto *osh_859 = buffer.data(osh + 859);
    const auto *osh_860 = buffer.data(osh + 860);
    const auto *osh_861 = buffer.data(osh + 861);
    const auto *osh_863 = buffer.data(osh + 863);
    const auto *osh_864 = buffer.data(osh + 864);
    const auto *osh_866 = buffer.data(osh + 866);
    const auto *osh_867 = buffer.data(osh + 867);
    const auto *osh_870 = buffer.data(osh + 870);
    const auto *osh_876 = buffer.data(osh + 876);
    const auto *osh_878 = buffer.data(osh + 878);
    const auto *osh_879 = buffer.data(osh + 879);
    const auto *osh_880 = buffer.data(osh + 880);
    const auto *osh_881 = buffer.data(osh + 881);
    const auto *osh_882 = buffer.data(osh + 882);
    const auto *osh_884 = buffer.data(osh + 884);
    const auto *osh_885 = buffer.data(osh + 885);
    const auto *osh_887 = buffer.data(osh + 887);
    const auto *osh_888 = buffer.data(osh + 888);
    const auto *osh_891 = buffer.data(osh + 891);
    const auto *osh_897 = buffer.data(osh + 897);
    const auto *osh_899 = buffer.data(osh + 899);
    const auto *osh_900 = buffer.data(osh + 900);
    const auto *osh_901 = buffer.data(osh + 901);
    const auto *osh_902 = buffer.data(osh + 902);
    const auto *osh_903 = buffer.data(osh + 903);
    const auto *osh_905 = buffer.data(osh + 905);
    const auto *osh_908 = buffer.data(osh + 908);
    const auto *osh_912 = buffer.data(osh + 912);
    const auto *osh_918 = buffer.data(osh + 918);
    const auto *osh_1029 = buffer.data(osh + 1029);
    const auto *osh_1032 = buffer.data(osh + 1032);
    const auto *osh_1034 = buffer.data(osh + 1034);
    const auto *osh_1035 = buffer.data(osh + 1035);
    const auto *osh_1038 = buffer.data(osh + 1038);
    const auto *osh_1039 = buffer.data(osh + 1039);
    const auto *osh_1041 = buffer.data(osh + 1041);
    const auto *osh_1043 = buffer.data(osh + 1043);
    const auto *osh_1044 = buffer.data(osh + 1044);
    const auto *osh_1045 = buffer.data(osh + 1045);
    const auto *osh_1046 = buffer.data(osh + 1046);
    const auto *osh_1047 = buffer.data(osh + 1047);
    const auto *osh_1048 = buffer.data(osh + 1048);
    const auto *osh_1049 = buffer.data(osh + 1049);
    const auto *osh_1050 = buffer.data(osh + 1050);
    const auto *osh_1053 = buffer.data(osh + 1053);
    const auto *osh_1055 = buffer.data(osh + 1055);
    const auto *osh_1056 = buffer.data(osh + 1056);
    const auto *osh_1059 = buffer.data(osh + 1059);
    const auto *osh_1060 = buffer.data(osh + 1060);
    const auto *osh_1062 = buffer.data(osh + 1062);
    const auto *osh_1064 = buffer.data(osh + 1064);
    const auto *osh_1065 = buffer.data(osh + 1065);
    const auto *osh_1066 = buffer.data(osh + 1066);
    const auto *osh_1067 = buffer.data(osh + 1067);
    const auto *osh_1068 = buffer.data(osh + 1068);
    const auto *osh_1069 = buffer.data(osh + 1069);
    const auto *osh_1070 = buffer.data(osh + 1070);
    const auto *osh_1071 = buffer.data(osh + 1071);
    const auto *osh_1074 = buffer.data(osh + 1074);
    const auto *osh_1076 = buffer.data(osh + 1076);
    const auto *osh_1077 = buffer.data(osh + 1077);
    const auto *osh_1080 = buffer.data(osh + 1080);
    const auto *osh_1081 = buffer.data(osh + 1081);
    const auto *osh_1083 = buffer.data(osh + 1083);
    const auto *osh_1085 = buffer.data(osh + 1085);
    const auto *osh_1086 = buffer.data(osh + 1086);
    const auto *osh_1087 = buffer.data(osh + 1087);
    const auto *osh_1088 = buffer.data(osh + 1088);
    const auto *osh_1089 = buffer.data(osh + 1089);
    const auto *osh_1090 = buffer.data(osh + 1090);
    const auto *osh_1091 = buffer.data(osh + 1091);
    const auto *osh_1092 = buffer.data(osh + 1092);
    const auto *osh_1095 = buffer.data(osh + 1095);
    const auto *osh_1097 = buffer.data(osh + 1097);
    const auto *osh_1098 = buffer.data(osh + 1098);
    const auto *osh_1101 = buffer.data(osh + 1101);
    const auto *osh_1102 = buffer.data(osh + 1102);
    const auto *osh_1104 = buffer.data(osh + 1104);
    const auto *osh_1106 = buffer.data(osh + 1106);
    const auto *osh_1107 = buffer.data(osh + 1107);
    const auto *osh_1108 = buffer.data(osh + 1108);
    const auto *osh_1109 = buffer.data(osh + 1109);
    const auto *osh_1110 = buffer.data(osh + 1110);
    const auto *osh_1111 = buffer.data(osh + 1111);
    const auto *osh_1112 = buffer.data(osh + 1112);

    const auto *qsg0_735 = buffer.data(qsg0 + 735);
    const auto *qsg0_738 = buffer.data(qsg0 + 738);
    const auto *qsg0_740 = buffer.data(qsg0 + 740);
    const auto *qsg0_741 = buffer.data(qsg0 + 741);
    const auto *qsg0_744 = buffer.data(qsg0 + 744);
    const auto *qsg0_745 = buffer.data(qsg0 + 745);
    const auto *qsg0_747 = buffer.data(qsg0 + 747);
    const auto *qsg0_748 = buffer.data(qsg0 + 748);
    const auto *qsg0_749 = buffer.data(qsg0 + 749);
    const auto *qsg0_750 = buffer.data(qsg0 + 750);
    const auto *qsg0_753 = buffer.data(qsg0 + 753);
    const auto *qsg0_755 = buffer.data(qsg0 + 755);
    const auto *qsg0_756 = buffer.data(qsg0 + 756);
    const auto *qsg0_759 = buffer.data(qsg0 + 759);
    const auto *qsg0_760 = buffer.data(qsg0 + 760);
    const auto *qsg0_762 = buffer.data(qsg0 + 762);
    const auto *qsg0_763 = buffer.data(qsg0 + 763);
    const auto *qsg0_764 = buffer.data(qsg0 + 764);
    const auto *qsg0_765 = buffer.data(qsg0 + 765);
    const auto *qsg0_768 = buffer.data(qsg0 + 768);
    const auto *qsg0_770 = buffer.data(qsg0 + 770);
    const auto *qsg0_771 = buffer.data(qsg0 + 771);
    const auto *qsg0_774 = buffer.data(qsg0 + 774);
    const auto *qsg0_775 = buffer.data(qsg0 + 775);
    const auto *qsg0_777 = buffer.data(qsg0 + 777);
    const auto *qsg0_778 = buffer.data(qsg0 + 778);
    const auto *qsg0_779 = buffer.data(qsg0 + 779);
    const auto *qsg0_780 = buffer.data(qsg0 + 780);
    const auto *qsg0_783 = buffer.data(qsg0 + 783);
    const auto *qsg0_785 = buffer.data(qsg0 + 785);
    const auto *qsg0_786 = buffer.data(qsg0 + 786);
    const auto *qsg0_789 = buffer.data(qsg0 + 789);
    const auto *qsg0_790 = buffer.data(qsg0 + 790);
    const auto *qsg0_792 = buffer.data(qsg0 + 792);
    const auto *qsg0_794 = buffer.data(qsg0 + 794);

    const auto *qsg1_735 = buffer.data(qsg1 + 735);
    const auto *qsg1_738 = buffer.data(qsg1 + 738);
    const auto *qsg1_740 = buffer.data(qsg1 + 740);
    const auto *qsg1_741 = buffer.data(qsg1 + 741);
    const auto *qsg1_744 = buffer.data(qsg1 + 744);
    const auto *qsg1_745 = buffer.data(qsg1 + 745);
    const auto *qsg1_747 = buffer.data(qsg1 + 747);
    const auto *qsg1_748 = buffer.data(qsg1 + 748);
    const auto *qsg1_749 = buffer.data(qsg1 + 749);
    const auto *qsg1_750 = buffer.data(qsg1 + 750);
    const auto *qsg1_753 = buffer.data(qsg1 + 753);
    const auto *qsg1_755 = buffer.data(qsg1 + 755);
    const auto *qsg1_756 = buffer.data(qsg1 + 756);
    const auto *qsg1_759 = buffer.data(qsg1 + 759);
    const auto *qsg1_760 = buffer.data(qsg1 + 760);
    const auto *qsg1_762 = buffer.data(qsg1 + 762);
    const auto *qsg1_763 = buffer.data(qsg1 + 763);
    const auto *qsg1_764 = buffer.data(qsg1 + 764);
    const auto *qsg1_765 = buffer.data(qsg1 + 765);
    const auto *qsg1_768 = buffer.data(qsg1 + 768);
    const auto *qsg1_770 = buffer.data(qsg1 + 770);
    const auto *qsg1_771 = buffer.data(qsg1 + 771);
    const auto *qsg1_774 = buffer.data(qsg1 + 774);
    const auto *qsg1_775 = buffer.data(qsg1 + 775);
    const auto *qsg1_777 = buffer.data(qsg1 + 777);
    const auto *qsg1_778 = buffer.data(qsg1 + 778);
    const auto *qsg1_779 = buffer.data(qsg1 + 779);
    const auto *qsg1_780 = buffer.data(qsg1 + 780);
    const auto *qsg1_783 = buffer.data(qsg1 + 783);
    const auto *qsg1_785 = buffer.data(qsg1 + 785);
    const auto *qsg1_786 = buffer.data(qsg1 + 786);
    const auto *qsg1_789 = buffer.data(qsg1 + 789);
    const auto *qsg1_790 = buffer.data(qsg1 + 790);
    const auto *qsg1_792 = buffer.data(qsg1 + 792);
    const auto *qsg1_794 = buffer.data(qsg1 + 794);

    const auto *qsh_1029 = buffer.data(qsh + 1029);
    const auto *qsh_1031 = buffer.data(qsh + 1031);
    const auto *qsh_1032 = buffer.data(qsh + 1032);
    const auto *qsh_1034 = buffer.data(qsh + 1034);
    const auto *qsh_1035 = buffer.data(qsh + 1035);
    const auto *qsh_1038 = buffer.data(qsh + 1038);
    const auto *qsh_1039 = buffer.data(qsh + 1039);
    const auto *qsh_1041 = buffer.data(qsh + 1041);
    const auto *qsh_1043 = buffer.data(qsh + 1043);
    const auto *qsh_1044 = buffer.data(qsh + 1044);
    const auto *qsh_1045 = buffer.data(qsh + 1045);
    const auto *qsh_1046 = buffer.data(qsh + 1046);
    const auto *qsh_1047 = buffer.data(qsh + 1047);
    const auto *qsh_1048 = buffer.data(qsh + 1048);
    const auto *qsh_1049 = buffer.data(qsh + 1049);
    const auto *qsh_1050 = buffer.data(qsh + 1050);
    const auto *qsh_1052 = buffer.data(qsh + 1052);
    const auto *qsh_1053 = buffer.data(qsh + 1053);
    const auto *qsh_1055 = buffer.data(qsh + 1055);
    const auto *qsh_1056 = buffer.data(qsh + 1056);
    const auto *qsh_1059 = buffer.data(qsh + 1059);
    const auto *qsh_1060 = buffer.data(qsh + 1060);
    const auto *qsh_1062 = buffer.data(qsh + 1062);
    const auto *qsh_1064 = buffer.data(qsh + 1064);
    const auto *qsh_1065 = buffer.data(qsh + 1065);
    const auto *qsh_1066 = buffer.data(qsh + 1066);
    const auto *qsh_1067 = buffer.data(qsh + 1067);
    const auto *qsh_1068 = buffer.data(qsh + 1068);
    const auto *qsh_1069 = buffer.data(qsh + 1069);
    const auto *qsh_1070 = buffer.data(qsh + 1070);
    const auto *qsh_1071 = buffer.data(qsh + 1071);
    const auto *qsh_1073 = buffer.data(qsh + 1073);
    const auto *qsh_1074 = buffer.data(qsh + 1074);
    const auto *qsh_1076 = buffer.data(qsh + 1076);
    const auto *qsh_1077 = buffer.data(qsh + 1077);
    const auto *qsh_1080 = buffer.data(qsh + 1080);
    const auto *qsh_1081 = buffer.data(qsh + 1081);
    const auto *qsh_1083 = buffer.data(qsh + 1083);
    const auto *qsh_1085 = buffer.data(qsh + 1085);
    const auto *qsh_1086 = buffer.data(qsh + 1086);
    const auto *qsh_1087 = buffer.data(qsh + 1087);
    const auto *qsh_1088 = buffer.data(qsh + 1088);
    const auto *qsh_1089 = buffer.data(qsh + 1089);
    const auto *qsh_1090 = buffer.data(qsh + 1090);
    const auto *qsh_1091 = buffer.data(qsh + 1091);
    const auto *qsh_1092 = buffer.data(qsh + 1092);
    const auto *qsh_1094 = buffer.data(qsh + 1094);
    const auto *qsh_1095 = buffer.data(qsh + 1095);
    const auto *qsh_1097 = buffer.data(qsh + 1097);
    const auto *qsh_1098 = buffer.data(qsh + 1098);
    const auto *qsh_1101 = buffer.data(qsh + 1101);
    const auto *qsh_1102 = buffer.data(qsh + 1102);
    const auto *qsh_1104 = buffer.data(qsh + 1104);
    const auto *qsh_1106 = buffer.data(qsh + 1106);
    const auto *qsh_1107 = buffer.data(qsh + 1107);
    const auto *qsh_1108 = buffer.data(qsh + 1108);
    const auto *qsh_1109 = buffer.data(qsh + 1109);
    const auto *qsh_1110 = buffer.data(qsh + 1110);
    const auto *qsh_1111 = buffer.data(qsh + 1111);
    const auto *qsh_1112 = buffer.data(qsh + 1112);

#pragma omp simd aligned(t_1372, t_1373, t_1374, pc_x, pc_y, pc_z, osh_819, osh_840, osh_1029, \
                         qsg0_735, qsg1_735, qsh_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1372[k] = f_13 * osh_1029[k]
                    + f_1 * qsg0_735[k]
                    - f_2 * qsg1_735[k]
                    + f_3 * pc_x[k] * qsh_1029[k];

        t_1373[k] = f_22 * osh_840[k]
                    + f_3 * pc_y[k] * qsh_1029[k];

        t_1374[k] = f_14 * osh_819[k]
                    + f_3 * pc_z[k] * qsh_1029[k];
    }

#pragma omp simd aligned(t_1375, t_1376, t_1377, pc_x, pc_y, osh_842, osh_1032, osh_1034, \
                         qsg0_738, qsg0_740, qsg1_738, qsg1_740, qsh_1031, qsh_1032, \
                         qsh_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1375[k] = f_13 * osh_1032[k]
                    + f_8 * qsg0_738[k]
                    - f_9 * qsg1_738[k]
                    + f_3 * pc_x[k] * qsh_1032[k];

        t_1376[k] = f_22 * osh_842[k]
                    + f_3 * pc_y[k] * qsh_1031[k];

        t_1377[k] = f_13 * osh_1034[k]
                    + f_8 * qsg0_740[k]
                    - f_9 * qsg1_740[k]
                    + f_3 * pc_x[k] * qsh_1034[k];
    }

#pragma omp simd aligned(t_1378, t_1379, t_1380, pc_x, pc_y, pc_z, osh_822, osh_845, osh_1035, \
                         qsg0_741, qsg1_741, qsh_1032, qsh_1034, \
                         qsh_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1378[k] = f_13 * osh_1035[k]
                    + f_6 * qsg0_741[k]
                    - f_7 * qsg1_741[k]
                    + f_3 * pc_x[k] * qsh_1035[k];

        t_1379[k] = f_14 * osh_822[k]
                    + f_3 * pc_z[k] * qsh_1032[k];

        t_1380[k] = f_22 * osh_845[k]
                    + f_3 * pc_y[k] * qsh_1034[k];
    }

#pragma omp simd aligned(t_1381, t_1382, t_1383, pc_x, pc_z, osh_825, osh_1038, osh_1039, \
                         qsg0_744, qsg0_745, qsg1_744, qsg1_745, qsh_1035, qsh_1038, \
                         qsh_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1381[k] = f_13 * osh_1038[k]
                    + f_6 * qsg0_744[k]
                    - f_7 * qsg1_744[k]
                    + f_3 * pc_x[k] * qsh_1038[k];

        t_1382[k] = f_13 * osh_1039[k]
                    + f_4 * qsg0_745[k]
                    - f_5 * qsg1_745[k]
                    + f_3 * pc_x[k] * qsh_1039[k];

        t_1383[k] = f_14 * osh_825[k]
                    + f_3 * pc_z[k] * qsh_1035[k];
    }

#pragma omp simd aligned(t_1384, t_1385, t_1386, pc_x, pc_y, osh_849, osh_1041, osh_1043, \
                         qsg0_747, qsg0_749, qsg1_747, qsg1_749, qsh_1038, qsh_1041, \
                         qsh_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1384[k] = f_13 * osh_1041[k]
                    + f_4 * qsg0_747[k]
                    - f_5 * qsg1_747[k]
                    + f_3 * pc_x[k] * qsh_1041[k];

        t_1385[k] = f_22 * osh_849[k]
                    + f_3 * pc_y[k] * qsh_1038[k];

        t_1386[k] = f_13 * osh_1043[k]
                    + f_4 * qsg0_749[k]
                    - f_5 * qsg1_749[k]
                    + f_3 * pc_x[k] * qsh_1043[k];
    }

#pragma omp simd aligned(t_1387, t_1388, t_1389, t_1390, t_1391, pc_x, osh_1044, osh_1045, \
                         osh_1046, osh_1047, osh_1048, qsh_1044, qsh_1045, qsh_1046, qsh_1047, \
                         qsh_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1387[k] = f_13 * osh_1044[k]
                    + f_3 * pc_x[k] * qsh_1044[k];

        t_1388[k] = f_13 * osh_1045[k]
                    + f_3 * pc_x[k] * qsh_1045[k];

        t_1389[k] = f_13 * osh_1046[k]
                    + f_3 * pc_x[k] * qsh_1046[k];

        t_1390[k] = f_13 * osh_1047[k]
                    + f_3 * pc_x[k] * qsh_1047[k];

        t_1391[k] = f_13 * osh_1048[k]
                    + f_3 * pc_x[k] * qsh_1048[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, pc_x, pc_y, pc_z, osh_834, osh_855, osh_1049, \
                         qsg0_745, qsg1_745, qsh_1044, qsh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_13 * osh_1049[k]
                    + f_3 * pc_x[k] * qsh_1049[k];

        t_1393[k] = f_22 * osh_855[k]
                    + f_1 * qsg0_745[k]
                    - f_2 * qsg1_745[k]
                    + f_3 * pc_y[k] * qsh_1044[k];

        t_1394[k] = f_14 * osh_834[k]
                    + f_3 * pc_z[k] * qsh_1044[k];
    }

#pragma omp simd aligned(t_1395, t_1396, t_1397, pc_y, osh_857, osh_858, osh_859, qsg0_747, \
                         qsg0_748, qsg0_749, qsg1_747, qsg1_748, qsg1_749, qsh_1046, qsh_1047, \
                         qsh_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1395[k] = f_22 * osh_857[k]
                    + f_8 * qsg0_747[k]
                    - f_9 * qsg1_747[k]
                    + f_3 * pc_y[k] * qsh_1046[k];

        t_1396[k] = f_22 * osh_858[k]
                    + f_6 * qsg0_748[k]
                    - f_7 * qsg1_748[k]
                    + f_3 * pc_y[k] * qsh_1047[k];

        t_1397[k] = f_22 * osh_859[k]
                    + f_4 * qsg0_749[k]
                    - f_5 * qsg1_749[k]
                    + f_3 * pc_y[k] * qsh_1048[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, pc_x, pc_y, pc_z, osh_839, osh_860, osh_1050, \
                         qsg0_749, qsg0_750, qsg1_749, qsg1_750, qsh_1049, \
                         qsh_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_22 * osh_860[k]
                    + f_3 * pc_y[k] * qsh_1049[k];

        t_1399[k] = f_14 * osh_839[k]
                    + f_1 * qsg0_749[k]
                    - f_2 * qsg1_749[k]
                    + f_3 * pc_z[k] * qsh_1049[k];

        t_1400[k] = f_13 * osh_1050[k]
                    + f_1 * qsg0_750[k]
                    - f_2 * qsg1_750[k]
                    + f_3 * pc_x[k] * qsh_1050[k];
    }

#pragma omp simd aligned(t_1401, t_1402, t_1403, t_1404, pc_x, pc_y, pc_z, osh_840, osh_861, \
                         osh_863, osh_1053, qsg0_753, qsg1_753, qsh_1050, qsh_1052, \
                         qsh_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1401[k] = f_14 * osh_861[k]
                    + f_3 * pc_y[k] * qsh_1050[k];

        t_1402[k] = f_22 * osh_840[k]
                    + f_3 * pc_z[k] * qsh_1050[k];

        t_1403[k] = f_13 * osh_1053[k]
                    + f_8 * qsg0_753[k]
                    - f_9 * qsg1_753[k]
                    + f_3 * pc_x[k] * qsh_1053[k];

        t_1404[k] = f_14 * osh_863[k]
                    + f_3 * pc_y[k] * qsh_1052[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, pc_x, pc_z, osh_843, osh_1055, osh_1056, \
                         qsg0_755, qsg0_756, qsg1_755, qsg1_756, qsh_1053, qsh_1055, \
                         qsh_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = f_13 * osh_1055[k]
                    + f_8 * qsg0_755[k]
                    - f_9 * qsg1_755[k]
                    + f_3 * pc_x[k] * qsh_1055[k];

        t_1406[k] = f_13 * osh_1056[k]
                    + f_6 * qsg0_756[k]
                    - f_7 * qsg1_756[k]
                    + f_3 * pc_x[k] * qsh_1056[k];

        t_1407[k] = f_22 * osh_843[k]
                    + f_3 * pc_z[k] * qsh_1053[k];
    }

#pragma omp simd aligned(t_1408, t_1409, t_1410, pc_x, pc_y, osh_866, osh_1059, osh_1060, \
                         qsg0_759, qsg0_760, qsg1_759, qsg1_760, qsh_1055, qsh_1059, \
                         qsh_1060 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1408[k] = f_14 * osh_866[k]
                    + f_3 * pc_y[k] * qsh_1055[k];

        t_1409[k] = f_13 * osh_1059[k]
                    + f_6 * qsg0_759[k]
                    - f_7 * qsg1_759[k]
                    + f_3 * pc_x[k] * qsh_1059[k];

        t_1410[k] = f_13 * osh_1060[k]
                    + f_4 * qsg0_760[k]
                    - f_5 * qsg1_760[k]
                    + f_3 * pc_x[k] * qsh_1060[k];
    }

#pragma omp simd aligned(t_1411, t_1412, t_1413, pc_x, pc_y, pc_z, osh_846, osh_870, osh_1062, \
                         qsg0_762, qsg1_762, qsh_1056, qsh_1059, \
                         qsh_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1411[k] = f_22 * osh_846[k]
                    + f_3 * pc_z[k] * qsh_1056[k];

        t_1412[k] = f_13 * osh_1062[k]
                    + f_4 * qsg0_762[k]
                    - f_5 * qsg1_762[k]
                    + f_3 * pc_x[k] * qsh_1062[k];

        t_1413[k] = f_14 * osh_870[k]
                    + f_3 * pc_y[k] * qsh_1059[k];
    }

#pragma omp simd aligned(t_1414, t_1415, t_1416, t_1417, pc_x, osh_1064, osh_1065, osh_1066, \
                         osh_1067, qsg0_764, qsg1_764, qsh_1064, qsh_1065, qsh_1066, \
                         qsh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1414[k] = f_13 * osh_1064[k]
                    + f_4 * qsg0_764[k]
                    - f_5 * qsg1_764[k]
                    + f_3 * pc_x[k] * qsh_1064[k];

        t_1415[k] = f_13 * osh_1065[k]
                    + f_3 * pc_x[k] * qsh_1065[k];

        t_1416[k] = f_13 * osh_1066[k]
                    + f_3 * pc_x[k] * qsh_1066[k];

        t_1417[k] = f_13 * osh_1067[k]
                    + f_3 * pc_x[k] * qsh_1067[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, t_1421, pc_x, pc_y, osh_876, osh_1068, \
                         osh_1069, osh_1070, qsg0_760, qsg1_760, qsh_1065, qsh_1068, qsh_1069, \
                         qsh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_13 * osh_1068[k]
                    + f_3 * pc_x[k] * qsh_1068[k];

        t_1419[k] = f_13 * osh_1069[k]
                    + f_3 * pc_x[k] * qsh_1069[k];

        t_1420[k] = f_13 * osh_1070[k]
                    + f_3 * pc_x[k] * qsh_1070[k];

        t_1421[k] = f_14 * osh_876[k]
                    + f_1 * qsg0_760[k]
                    - f_2 * qsg1_760[k]
                    + f_3 * pc_y[k] * qsh_1065[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, pc_y, pc_z, osh_855, osh_878, osh_879, \
                         qsg0_762, qsg0_763, qsg1_762, qsg1_763, qsh_1065, qsh_1067, \
                         qsh_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_22 * osh_855[k]
                    + f_3 * pc_z[k] * qsh_1065[k];

        t_1423[k] = f_14 * osh_878[k]
                    + f_8 * qsg0_762[k]
                    - f_9 * qsg1_762[k]
                    + f_3 * pc_y[k] * qsh_1067[k];

        t_1424[k] = f_14 * osh_879[k]
                    + f_6 * qsg0_763[k]
                    - f_7 * qsg1_763[k]
                    + f_3 * pc_y[k] * qsh_1068[k];
    }

#pragma omp simd aligned(t_1425, t_1426, t_1427, pc_y, pc_z, osh_860, osh_880, osh_881, \
                         qsg0_764, qsg1_764, qsh_1069, qsh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1425[k] = f_14 * osh_880[k]
                    + f_4 * qsg0_764[k]
                    - f_5 * qsg1_764[k]
                    + f_3 * pc_y[k] * qsh_1069[k];

        t_1426[k] = f_14 * osh_881[k]
                    + f_3 * pc_y[k] * qsh_1070[k];

        t_1427[k] = f_22 * osh_860[k]
                    + f_1 * qsg0_764[k]
                    - f_2 * qsg1_764[k]
                    + f_3 * pc_z[k] * qsh_1070[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, pc_x, pc_y, pc_z, osh_861, osh_882, osh_1071, \
                         qsg0_765, qsg1_765, qsh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_13 * osh_1071[k]
                    + f_1 * qsg0_765[k]
                    - f_2 * qsg1_765[k]
                    + f_3 * pc_x[k] * qsh_1071[k];

        t_1429[k] = f_13 * osh_882[k]
                    + f_3 * pc_y[k] * qsh_1071[k];

        t_1430[k] = f_23 * osh_861[k]
                    + f_3 * pc_z[k] * qsh_1071[k];
    }

#pragma omp simd aligned(t_1431, t_1432, t_1433, pc_x, pc_y, osh_884, osh_1074, osh_1076, \
                         qsg0_768, qsg0_770, qsg1_768, qsg1_770, qsh_1073, qsh_1074, \
                         qsh_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1431[k] = f_13 * osh_1074[k]
                    + f_8 * qsg0_768[k]
                    - f_9 * qsg1_768[k]
                    + f_3 * pc_x[k] * qsh_1074[k];

        t_1432[k] = f_13 * osh_884[k]
                    + f_3 * pc_y[k] * qsh_1073[k];

        t_1433[k] = f_13 * osh_1076[k]
                    + f_8 * qsg0_770[k]
                    - f_9 * qsg1_770[k]
                    + f_3 * pc_x[k] * qsh_1076[k];
    }

#pragma omp simd aligned(t_1434, t_1435, t_1436, pc_x, pc_y, pc_z, osh_864, osh_887, osh_1077, \
                         qsg0_771, qsg1_771, qsh_1074, qsh_1076, \
                         qsh_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1434[k] = f_13 * osh_1077[k]
                    + f_6 * qsg0_771[k]
                    - f_7 * qsg1_771[k]
                    + f_3 * pc_x[k] * qsh_1077[k];

        t_1435[k] = f_23 * osh_864[k]
                    + f_3 * pc_z[k] * qsh_1074[k];

        t_1436[k] = f_13 * osh_887[k]
                    + f_3 * pc_y[k] * qsh_1076[k];
    }

#pragma omp simd aligned(t_1437, t_1438, t_1439, pc_x, pc_z, osh_867, osh_1080, osh_1081, \
                         qsg0_774, qsg0_775, qsg1_774, qsg1_775, qsh_1077, qsh_1080, \
                         qsh_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1437[k] = f_13 * osh_1080[k]
                    + f_6 * qsg0_774[k]
                    - f_7 * qsg1_774[k]
                    + f_3 * pc_x[k] * qsh_1080[k];

        t_1438[k] = f_13 * osh_1081[k]
                    + f_4 * qsg0_775[k]
                    - f_5 * qsg1_775[k]
                    + f_3 * pc_x[k] * qsh_1081[k];

        t_1439[k] = f_23 * osh_867[k]
                    + f_3 * pc_z[k] * qsh_1077[k];
    }

#pragma omp simd aligned(t_1440, t_1441, t_1442, pc_x, pc_y, osh_891, osh_1083, osh_1085, \
                         qsg0_777, qsg0_779, qsg1_777, qsg1_779, qsh_1080, qsh_1083, \
                         qsh_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1440[k] = f_13 * osh_1083[k]
                    + f_4 * qsg0_777[k]
                    - f_5 * qsg1_777[k]
                    + f_3 * pc_x[k] * qsh_1083[k];

        t_1441[k] = f_13 * osh_891[k]
                    + f_3 * pc_y[k] * qsh_1080[k];

        t_1442[k] = f_13 * osh_1085[k]
                    + f_4 * qsg0_779[k]
                    - f_5 * qsg1_779[k]
                    + f_3 * pc_x[k] * qsh_1085[k];
    }

#pragma omp simd aligned(t_1443, t_1444, t_1445, t_1446, t_1447, pc_x, osh_1086, osh_1087, \
                         osh_1088, osh_1089, osh_1090, qsh_1086, qsh_1087, qsh_1088, qsh_1089, \
                         qsh_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1443[k] = f_13 * osh_1086[k]
                    + f_3 * pc_x[k] * qsh_1086[k];

        t_1444[k] = f_13 * osh_1087[k]
                    + f_3 * pc_x[k] * qsh_1087[k];

        t_1445[k] = f_13 * osh_1088[k]
                    + f_3 * pc_x[k] * qsh_1088[k];

        t_1446[k] = f_13 * osh_1089[k]
                    + f_3 * pc_x[k] * qsh_1089[k];

        t_1447[k] = f_13 * osh_1090[k]
                    + f_3 * pc_x[k] * qsh_1090[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, pc_y, pc_z, osh_876, osh_897, osh_1091, \
                         qsg0_775, qsg1_775, qsh_1086, qsh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_13 * osh_1091[k]
                    + f_3 * pc_x[k] * qsh_1091[k];

        t_1449[k] = f_13 * osh_897[k]
                    + f_1 * qsg0_775[k]
                    - f_2 * qsg1_775[k]
                    + f_3 * pc_y[k] * qsh_1086[k];

        t_1450[k] = f_23 * osh_876[k]
                    + f_3 * pc_z[k] * qsh_1086[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_y, osh_899, osh_900, osh_901, qsg0_777, \
                         qsg0_778, qsg0_779, qsg1_777, qsg1_778, qsg1_779, qsh_1088, qsh_1089, \
                         qsh_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_13 * osh_899[k]
                    + f_8 * qsg0_777[k]
                    - f_9 * qsg1_777[k]
                    + f_3 * pc_y[k] * qsh_1088[k];

        t_1452[k] = f_13 * osh_900[k]
                    + f_6 * qsg0_778[k]
                    - f_7 * qsg1_778[k]
                    + f_3 * pc_y[k] * qsh_1089[k];

        t_1453[k] = f_13 * osh_901[k]
                    + f_4 * qsg0_779[k]
                    - f_5 * qsg1_779[k]
                    + f_3 * pc_y[k] * qsh_1090[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, pc_y, pc_z, osh_881, osh_902, osh_1092, \
                         qsg0_779, qsg0_780, qsg1_779, qsg1_780, qsh_1091, \
                         qsh_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_13 * osh_902[k]
                    + f_3 * pc_y[k] * qsh_1091[k];

        t_1455[k] = f_23 * osh_881[k]
                    + f_1 * qsg0_779[k]
                    - f_2 * qsg1_779[k]
                    + f_3 * pc_z[k] * qsh_1091[k];

        t_1456[k] = f_13 * osh_1092[k]
                    + f_1 * qsg0_780[k]
                    - f_2 * qsg1_780[k]
                    + f_3 * pc_x[k] * qsh_1092[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, t_1460, pc_x, pc_y, pc_z, osh_882, osh_903, \
                         osh_905, osh_1095, qsg0_783, qsg1_783, qsh_1092, qsh_1094, \
                         qsh_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_12 * osh_903[k]
                    + f_3 * pc_y[k] * qsh_1092[k];

        t_1458[k] = f_21 * osh_882[k]
                    + f_3 * pc_z[k] * qsh_1092[k];

        t_1459[k] = f_13 * osh_1095[k]
                    + f_8 * qsg0_783[k]
                    - f_9 * qsg1_783[k]
                    + f_3 * pc_x[k] * qsh_1095[k];

        t_1460[k] = f_12 * osh_905[k]
                    + f_3 * pc_y[k] * qsh_1094[k];
    }

#pragma omp simd aligned(t_1461, t_1462, t_1463, pc_x, pc_z, osh_885, osh_1097, osh_1098, \
                         qsg0_785, qsg0_786, qsg1_785, qsg1_786, qsh_1095, qsh_1097, \
                         qsh_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1461[k] = f_13 * osh_1097[k]
                    + f_8 * qsg0_785[k]
                    - f_9 * qsg1_785[k]
                    + f_3 * pc_x[k] * qsh_1097[k];

        t_1462[k] = f_13 * osh_1098[k]
                    + f_6 * qsg0_786[k]
                    - f_7 * qsg1_786[k]
                    + f_3 * pc_x[k] * qsh_1098[k];

        t_1463[k] = f_21 * osh_885[k]
                    + f_3 * pc_z[k] * qsh_1095[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, pc_x, pc_y, osh_908, osh_1101, osh_1102, \
                         qsg0_789, qsg0_790, qsg1_789, qsg1_790, qsh_1097, qsh_1101, \
                         qsh_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_12 * osh_908[k]
                    + f_3 * pc_y[k] * qsh_1097[k];

        t_1465[k] = f_13 * osh_1101[k]
                    + f_6 * qsg0_789[k]
                    - f_7 * qsg1_789[k]
                    + f_3 * pc_x[k] * qsh_1101[k];

        t_1466[k] = f_13 * osh_1102[k]
                    + f_4 * qsg0_790[k]
                    - f_5 * qsg1_790[k]
                    + f_3 * pc_x[k] * qsh_1102[k];
    }

#pragma omp simd aligned(t_1467, t_1468, t_1469, pc_x, pc_y, pc_z, osh_888, osh_912, osh_1104, \
                         qsg0_792, qsg1_792, qsh_1098, qsh_1101, \
                         qsh_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1467[k] = f_21 * osh_888[k]
                    + f_3 * pc_z[k] * qsh_1098[k];

        t_1468[k] = f_13 * osh_1104[k]
                    + f_4 * qsg0_792[k]
                    - f_5 * qsg1_792[k]
                    + f_3 * pc_x[k] * qsh_1104[k];

        t_1469[k] = f_12 * osh_912[k]
                    + f_3 * pc_y[k] * qsh_1101[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, pc_x, osh_1106, osh_1107, osh_1108, \
                         osh_1109, qsg0_794, qsg1_794, qsh_1106, qsh_1107, qsh_1108, \
                         qsh_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = f_13 * osh_1106[k]
                    + f_4 * qsg0_794[k]
                    - f_5 * qsg1_794[k]
                    + f_3 * pc_x[k] * qsh_1106[k];

        t_1471[k] = f_13 * osh_1107[k]
                    + f_3 * pc_x[k] * qsh_1107[k];

        t_1472[k] = f_13 * osh_1108[k]
                    + f_3 * pc_x[k] * qsh_1108[k];

        t_1473[k] = f_13 * osh_1109[k]
                    + f_3 * pc_x[k] * qsh_1109[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, pc_x, pc_y, osh_918, osh_1110, \
                         osh_1111, osh_1112, qsg0_790, qsg1_790, qsh_1107, qsh_1110, qsh_1111, \
                         qsh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_13 * osh_1110[k]
                    + f_3 * pc_x[k] * qsh_1110[k];

        t_1475[k] = f_13 * osh_1111[k]
                    + f_3 * pc_x[k] * qsh_1111[k];

        t_1476[k] = f_13 * osh_1112[k]
                    + f_3 * pc_x[k] * qsh_1112[k];

        t_1477[k] = f_12 * osh_918[k]
                    + f_1 * qsg0_790[k]
                    - f_2 * qsg1_790[k]
                    + f_3 * pc_y[k] * qsh_1107[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsg0, const size_t qsg1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 5.0 / q;
    const auto f_19 = 4.5 / q;
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;

    auto *t_1478 = buffer.data(target + 1478);
    auto *t_1479 = buffer.data(target + 1479);
    auto *t_1480 = buffer.data(target + 1480);
    auto *t_1481 = buffer.data(target + 1481);
    auto *t_1482 = buffer.data(target + 1482);
    auto *t_1483 = buffer.data(target + 1483);
    auto *t_1484 = buffer.data(target + 1484);
    auto *t_1485 = buffer.data(target + 1485);
    auto *t_1486 = buffer.data(target + 1486);
    auto *t_1487 = buffer.data(target + 1487);
    auto *t_1488 = buffer.data(target + 1488);
    auto *t_1489 = buffer.data(target + 1489);
    auto *t_1490 = buffer.data(target + 1490);
    auto *t_1491 = buffer.data(target + 1491);
    auto *t_1492 = buffer.data(target + 1492);
    auto *t_1493 = buffer.data(target + 1493);
    auto *t_1494 = buffer.data(target + 1494);
    auto *t_1495 = buffer.data(target + 1495);
    auto *t_1496 = buffer.data(target + 1496);
    auto *t_1497 = buffer.data(target + 1497);
    auto *t_1498 = buffer.data(target + 1498);
    auto *t_1499 = buffer.data(target + 1499);
    auto *t_1500 = buffer.data(target + 1500);
    auto *t_1501 = buffer.data(target + 1501);
    auto *t_1502 = buffer.data(target + 1502);
    auto *t_1503 = buffer.data(target + 1503);
    auto *t_1504 = buffer.data(target + 1504);
    auto *t_1505 = buffer.data(target + 1505);
    auto *t_1506 = buffer.data(target + 1506);
    auto *t_1507 = buffer.data(target + 1507);
    auto *t_1508 = buffer.data(target + 1508);
    auto *t_1509 = buffer.data(target + 1509);
    auto *t_1510 = buffer.data(target + 1510);
    auto *t_1511 = buffer.data(target + 1511);
    auto *t_1512 = buffer.data(target + 1512);
    auto *t_1513 = buffer.data(target + 1513);
    auto *t_1514 = buffer.data(target + 1514);
    auto *t_1515 = buffer.data(target + 1515);
    auto *t_1516 = buffer.data(target + 1516);
    auto *t_1517 = buffer.data(target + 1517);
    auto *t_1518 = buffer.data(target + 1518);
    auto *t_1519 = buffer.data(target + 1519);
    auto *t_1520 = buffer.data(target + 1520);
    auto *t_1521 = buffer.data(target + 1521);
    auto *t_1522 = buffer.data(target + 1522);
    auto *t_1523 = buffer.data(target + 1523);
    auto *t_1524 = buffer.data(target + 1524);
    auto *t_1525 = buffer.data(target + 1525);
    auto *t_1526 = buffer.data(target + 1526);
    auto *t_1527 = buffer.data(target + 1527);
    auto *t_1528 = buffer.data(target + 1528);
    auto *t_1529 = buffer.data(target + 1529);
    auto *t_1530 = buffer.data(target + 1530);
    auto *t_1531 = buffer.data(target + 1531);
    auto *t_1532 = buffer.data(target + 1532);
    auto *t_1533 = buffer.data(target + 1533);
    auto *t_1534 = buffer.data(target + 1534);
    auto *t_1535 = buffer.data(target + 1535);
    auto *t_1536 = buffer.data(target + 1536);
    auto *t_1537 = buffer.data(target + 1537);
    auto *t_1538 = buffer.data(target + 1538);
    auto *t_1539 = buffer.data(target + 1539);
    auto *t_1540 = buffer.data(target + 1540);
    auto *t_1541 = buffer.data(target + 1541);
    auto *t_1542 = buffer.data(target + 1542);
    auto *t_1543 = buffer.data(target + 1543);
    auto *t_1544 = buffer.data(target + 1544);
    auto *t_1545 = buffer.data(target + 1545);
    auto *t_1546 = buffer.data(target + 1546);
    auto *t_1547 = buffer.data(target + 1547);
    auto *t_1548 = buffer.data(target + 1548);
    auto *t_1549 = buffer.data(target + 1549);
    auto *t_1550 = buffer.data(target + 1550);
    auto *t_1551 = buffer.data(target + 1551);
    auto *t_1552 = buffer.data(target + 1552);
    auto *t_1553 = buffer.data(target + 1553);
    auto *t_1554 = buffer.data(target + 1554);
    auto *t_1555 = buffer.data(target + 1555);
    auto *t_1556 = buffer.data(target + 1556);
    auto *t_1557 = buffer.data(target + 1557);
    auto *t_1558 = buffer.data(target + 1558);
    auto *t_1559 = buffer.data(target + 1559);
    auto *t_1560 = buffer.data(target + 1560);
    auto *t_1561 = buffer.data(target + 1561);
    auto *t_1562 = buffer.data(target + 1562);
    auto *t_1563 = buffer.data(target + 1563);
    auto *t_1564 = buffer.data(target + 1564);
    auto *t_1565 = buffer.data(target + 1565);
    auto *t_1566 = buffer.data(target + 1566);
    auto *t_1567 = buffer.data(target + 1567);
    auto *t_1568 = buffer.data(target + 1568);
    auto *t_1569 = buffer.data(target + 1569);
    auto *t_1570 = buffer.data(target + 1570);
    auto *t_1571 = buffer.data(target + 1571);
    auto *t_1572 = buffer.data(target + 1572);
    auto *t_1573 = buffer.data(target + 1573);
    auto *t_1574 = buffer.data(target + 1574);
    auto *t_1575 = buffer.data(target + 1575);
    auto *t_1576 = buffer.data(target + 1576);
    auto *t_1577 = buffer.data(target + 1577);
    auto *t_1578 = buffer.data(target + 1578);
    auto *t_1579 = buffer.data(target + 1579);
    auto *t_1580 = buffer.data(target + 1580);
    auto *t_1581 = buffer.data(target + 1581);
    auto *t_1582 = buffer.data(target + 1582);
    auto *t_1583 = buffer.data(target + 1583);
    auto *t_1584 = buffer.data(target + 1584);
    auto *t_1585 = buffer.data(target + 1585);
    auto *t_1586 = buffer.data(target + 1586);
    auto *t_1587 = buffer.data(target + 1587);
    auto *t_1588 = buffer.data(target + 1588);
    auto *t_1589 = buffer.data(target + 1589);
    auto *t_1590 = buffer.data(target + 1590);
    auto *t_1591 = buffer.data(target + 1591);
    auto *t_1592 = buffer.data(target + 1592);
    auto *t_1593 = buffer.data(target + 1593);
    auto *t_1594 = buffer.data(target + 1594);
    auto *t_1595 = buffer.data(target + 1595);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_1232 = buffer.data(osi0 + 1232);
    const auto *osi0_1235 = buffer.data(osi0 + 1235);
    const auto *osi0_1237 = buffer.data(osi0 + 1237);
    const auto *osi0_1238 = buffer.data(osi0 + 1238);
    const auto *osi0_1241 = buffer.data(osi0 + 1241);
    const auto *osi0_1242 = buffer.data(osi0 + 1242);
    const auto *osi0_1244 = buffer.data(osi0 + 1244);
    const auto *osi0_1246 = buffer.data(osi0 + 1246);
    const auto *osi0_1259 = buffer.data(osi0 + 1259);
    const auto *osi0_1260 = buffer.data(osi0 + 1260);
    const auto *osi0_1263 = buffer.data(osi0 + 1263);
    const auto *osi0_1266 = buffer.data(osi0 + 1266);
    const auto *osi0_1270 = buffer.data(osi0 + 1270);
    const auto *osi0_1272 = buffer.data(osi0 + 1272);
    const auto *osi0_1281 = buffer.data(osi0 + 1281);

    const auto *osh_897 = buffer.data(osh + 897);
    const auto *osh_902 = buffer.data(osh + 902);
    const auto *osh_903 = buffer.data(osh + 903);
    const auto *osh_906 = buffer.data(osh + 906);
    const auto *osh_909 = buffer.data(osh + 909);
    const auto *osh_918 = buffer.data(osh + 918);
    const auto *osh_920 = buffer.data(osh + 920);
    const auto *osh_921 = buffer.data(osh + 921);
    const auto *osh_922 = buffer.data(osh + 922);
    const auto *osh_923 = buffer.data(osh + 923);
    const auto *osh_924 = buffer.data(osh + 924);
    const auto *osh_925 = buffer.data(osh + 925);
    const auto *osh_926 = buffer.data(osh + 926);
    const auto *osh_927 = buffer.data(osh + 927);
    const auto *osh_929 = buffer.data(osh + 929);
    const auto *osh_930 = buffer.data(osh + 930);
    const auto *osh_932 = buffer.data(osh + 932);
    const auto *osh_933 = buffer.data(osh + 933);
    const auto *osh_939 = buffer.data(osh + 939);
    const auto *osh_941 = buffer.data(osh + 941);
    const auto *osh_942 = buffer.data(osh + 942);
    const auto *osh_943 = buffer.data(osh + 943);
    const auto *osh_944 = buffer.data(osh + 944);
    const auto *osh_945 = buffer.data(osh + 945);
    const auto *osh_948 = buffer.data(osh + 948);
    const auto *osh_950 = buffer.data(osh + 950);
    const auto *osh_951 = buffer.data(osh + 951);
    const auto *osh_952 = buffer.data(osh + 952);
    const auto *osh_954 = buffer.data(osh + 954);
    const auto *osh_960 = buffer.data(osh + 960);
    const auto *osh_965 = buffer.data(osh + 965);
    const auto *osh_966 = buffer.data(osh + 966);
    const auto *osh_968 = buffer.data(osh + 968);
    const auto *osh_971 = buffer.data(osh + 971);
    const auto *osh_975 = buffer.data(osh + 975);
    const auto *osh_983 = buffer.data(osh + 983);
    const auto *osh_984 = buffer.data(osh + 984);
    const auto *osh_985 = buffer.data(osh + 985);
    const auto *osh_986 = buffer.data(osh + 986);
    const auto *osh_1128 = buffer.data(osh + 1128);
    const auto *osh_1129 = buffer.data(osh + 1129);
    const auto *osh_1130 = buffer.data(osh + 1130);
    const auto *osh_1131 = buffer.data(osh + 1131);
    const auto *osh_1132 = buffer.data(osh + 1132);
    const auto *osh_1133 = buffer.data(osh + 1133);
    const auto *osh_1134 = buffer.data(osh + 1134);
    const auto *osh_1139 = buffer.data(osh + 1139);
    const auto *osh_1143 = buffer.data(osh + 1143);
    const auto *osh_1148 = buffer.data(osh + 1148);
    const auto *osh_1149 = buffer.data(osh + 1149);
    const auto *osh_1150 = buffer.data(osh + 1150);
    const auto *osh_1151 = buffer.data(osh + 1151);
    const auto *osh_1152 = buffer.data(osh + 1152);
    const auto *osh_1154 = buffer.data(osh + 1154);
    const auto *osh_1155 = buffer.data(osh + 1155);
    const auto *osh_1158 = buffer.data(osh + 1158);
    const auto *osh_1161 = buffer.data(osh + 1161);
    const auto *osh_1165 = buffer.data(osh + 1165);
    const auto *osh_1170 = buffer.data(osh + 1170);
    const auto *osh_1172 = buffer.data(osh + 1172);
    const auto *osh_1173 = buffer.data(osh + 1173);
    const auto *osh_1174 = buffer.data(osh + 1174);
    const auto *osh_1175 = buffer.data(osh + 1175);
    const auto *osh_1181 = buffer.data(osh + 1181);
    const auto *osh_1185 = buffer.data(osh + 1185);
    const auto *osh_1190 = buffer.data(osh + 1190);
    const auto *osh_1191 = buffer.data(osh + 1191);
    const auto *osh_1192 = buffer.data(osh + 1192);
    const auto *osh_1193 = buffer.data(osh + 1193);
    const auto *osh_1194 = buffer.data(osh + 1194);
    const auto *osh_1195 = buffer.data(osh + 1195);
    const auto *osh_1196 = buffer.data(osh + 1196);

    const auto *osi1_1232 = buffer.data(osi1 + 1232);
    const auto *osi1_1235 = buffer.data(osi1 + 1235);
    const auto *osi1_1237 = buffer.data(osi1 + 1237);
    const auto *osi1_1238 = buffer.data(osi1 + 1238);
    const auto *osi1_1241 = buffer.data(osi1 + 1241);
    const auto *osi1_1242 = buffer.data(osi1 + 1242);
    const auto *osi1_1244 = buffer.data(osi1 + 1244);
    const auto *osi1_1246 = buffer.data(osi1 + 1246);
    const auto *osi1_1259 = buffer.data(osi1 + 1259);
    const auto *osi1_1260 = buffer.data(osi1 + 1260);
    const auto *osi1_1263 = buffer.data(osi1 + 1263);
    const auto *osi1_1266 = buffer.data(osi1 + 1266);
    const auto *osi1_1270 = buffer.data(osi1 + 1270);
    const auto *osi1_1272 = buffer.data(osi1 + 1272);
    const auto *osi1_1281 = buffer.data(osi1 + 1281);

    const auto *qsg0_792 = buffer.data(qsg0 + 792);
    const auto *qsg0_793 = buffer.data(qsg0 + 793);
    const auto *qsg0_794 = buffer.data(qsg0 + 794);
    const auto *qsg0_805 = buffer.data(qsg0 + 805);
    const auto *qsg0_807 = buffer.data(qsg0 + 807);
    const auto *qsg0_808 = buffer.data(qsg0 + 808);
    const auto *qsg0_809 = buffer.data(qsg0 + 809);
    const auto *qsg0_810 = buffer.data(qsg0 + 810);
    const auto *qsg0_811 = buffer.data(qsg0 + 811);
    const auto *qsg0_812 = buffer.data(qsg0 + 812);
    const auto *qsg0_813 = buffer.data(qsg0 + 813);
    const auto *qsg0_814 = buffer.data(qsg0 + 814);
    const auto *qsg0_815 = buffer.data(qsg0 + 815);
    const auto *qsg0_819 = buffer.data(qsg0 + 819);
    const auto *qsg0_820 = buffer.data(qsg0 + 820);
    const auto *qsg0_821 = buffer.data(qsg0 + 821);
    const auto *qsg0_822 = buffer.data(qsg0 + 822);
    const auto *qsg0_823 = buffer.data(qsg0 + 823);
    const auto *qsg0_824 = buffer.data(qsg0 + 824);
    const auto *qsg0_825 = buffer.data(qsg0 + 825);
    const auto *qsg0_827 = buffer.data(qsg0 + 827);
    const auto *qsg0_828 = buffer.data(qsg0 + 828);
    const auto *qsg0_830 = buffer.data(qsg0 + 830);
    const auto *qsg0_831 = buffer.data(qsg0 + 831);
    const auto *qsg0_835 = buffer.data(qsg0 + 835);
    const auto *qsg0_836 = buffer.data(qsg0 + 836);
    const auto *qsg0_837 = buffer.data(qsg0 + 837);
    const auto *qsg0_839 = buffer.data(qsg0 + 839);
    const auto *qsg0_845 = buffer.data(qsg0 + 845);
    const auto *qsg0_849 = buffer.data(qsg0 + 849);
    const auto *qsg0_852 = buffer.data(qsg0 + 852);
    const auto *qsg0_853 = buffer.data(qsg0 + 853);
    const auto *qsg0_854 = buffer.data(qsg0 + 854);

    const auto *qsg1_792 = buffer.data(qsg1 + 792);
    const auto *qsg1_793 = buffer.data(qsg1 + 793);
    const auto *qsg1_794 = buffer.data(qsg1 + 794);
    const auto *qsg1_805 = buffer.data(qsg1 + 805);
    const auto *qsg1_807 = buffer.data(qsg1 + 807);
    const auto *qsg1_808 = buffer.data(qsg1 + 808);
    const auto *qsg1_809 = buffer.data(qsg1 + 809);
    const auto *qsg1_810 = buffer.data(qsg1 + 810);
    const auto *qsg1_811 = buffer.data(qsg1 + 811);
    const auto *qsg1_812 = buffer.data(qsg1 + 812);
    const auto *qsg1_813 = buffer.data(qsg1 + 813);
    const auto *qsg1_814 = buffer.data(qsg1 + 814);
    const auto *qsg1_815 = buffer.data(qsg1 + 815);
    const auto *qsg1_819 = buffer.data(qsg1 + 819);
    const auto *qsg1_820 = buffer.data(qsg1 + 820);
    const auto *qsg1_821 = buffer.data(qsg1 + 821);
    const auto *qsg1_822 = buffer.data(qsg1 + 822);
    const auto *qsg1_823 = buffer.data(qsg1 + 823);
    const auto *qsg1_824 = buffer.data(qsg1 + 824);
    const auto *qsg1_825 = buffer.data(qsg1 + 825);
    const auto *qsg1_827 = buffer.data(qsg1 + 827);
    const auto *qsg1_828 = buffer.data(qsg1 + 828);
    const auto *qsg1_830 = buffer.data(qsg1 + 830);
    const auto *qsg1_831 = buffer.data(qsg1 + 831);
    const auto *qsg1_835 = buffer.data(qsg1 + 835);
    const auto *qsg1_836 = buffer.data(qsg1 + 836);
    const auto *qsg1_837 = buffer.data(qsg1 + 837);
    const auto *qsg1_839 = buffer.data(qsg1 + 839);
    const auto *qsg1_845 = buffer.data(qsg1 + 845);
    const auto *qsg1_849 = buffer.data(qsg1 + 849);
    const auto *qsg1_852 = buffer.data(qsg1 + 852);
    const auto *qsg1_853 = buffer.data(qsg1 + 853);
    const auto *qsg1_854 = buffer.data(qsg1 + 854);

    const auto *qsh_1107 = buffer.data(qsh + 1107);
    const auto *qsh_1109 = buffer.data(qsh + 1109);
    const auto *qsh_1110 = buffer.data(qsh + 1110);
    const auto *qsh_1111 = buffer.data(qsh + 1111);
    const auto *qsh_1112 = buffer.data(qsh + 1112);
    const auto *qsh_1113 = buffer.data(qsh + 1113);
    const auto *qsh_1115 = buffer.data(qsh + 1115);
    const auto *qsh_1116 = buffer.data(qsh + 1116);
    const auto *qsh_1118 = buffer.data(qsh + 1118);
    const auto *qsh_1119 = buffer.data(qsh + 1119);
    const auto *qsh_1122 = buffer.data(qsh + 1122);
    const auto *qsh_1128 = buffer.data(qsh + 1128);
    const auto *qsh_1129 = buffer.data(qsh + 1129);
    const auto *qsh_1130 = buffer.data(qsh + 1130);
    const auto *qsh_1131 = buffer.data(qsh + 1131);
    const auto *qsh_1132 = buffer.data(qsh + 1132);
    const auto *qsh_1133 = buffer.data(qsh + 1133);
    const auto *qsh_1134 = buffer.data(qsh + 1134);
    const auto *qsh_1135 = buffer.data(qsh + 1135);
    const auto *qsh_1136 = buffer.data(qsh + 1136);
    const auto *qsh_1137 = buffer.data(qsh + 1137);
    const auto *qsh_1138 = buffer.data(qsh + 1138);
    const auto *qsh_1139 = buffer.data(qsh + 1139);
    const auto *qsh_1140 = buffer.data(qsh + 1140);
    const auto *qsh_1141 = buffer.data(qsh + 1141);
    const auto *qsh_1142 = buffer.data(qsh + 1142);
    const auto *qsh_1143 = buffer.data(qsh + 1143);
    const auto *qsh_1148 = buffer.data(qsh + 1148);
    const auto *qsh_1149 = buffer.data(qsh + 1149);
    const auto *qsh_1150 = buffer.data(qsh + 1150);
    const auto *qsh_1151 = buffer.data(qsh + 1151);
    const auto *qsh_1152 = buffer.data(qsh + 1152);
    const auto *qsh_1153 = buffer.data(qsh + 1153);
    const auto *qsh_1154 = buffer.data(qsh + 1154);
    const auto *qsh_1155 = buffer.data(qsh + 1155);
    const auto *qsh_1156 = buffer.data(qsh + 1156);
    const auto *qsh_1157 = buffer.data(qsh + 1157);
    const auto *qsh_1158 = buffer.data(qsh + 1158);
    const auto *qsh_1160 = buffer.data(qsh + 1160);
    const auto *qsh_1161 = buffer.data(qsh + 1161);
    const auto *qsh_1162 = buffer.data(qsh + 1162);
    const auto *qsh_1164 = buffer.data(qsh + 1164);
    const auto *qsh_1165 = buffer.data(qsh + 1165);
    const auto *qsh_1170 = buffer.data(qsh + 1170);
    const auto *qsh_1171 = buffer.data(qsh + 1171);
    const auto *qsh_1172 = buffer.data(qsh + 1172);
    const auto *qsh_1173 = buffer.data(qsh + 1173);
    const auto *qsh_1174 = buffer.data(qsh + 1174);
    const auto *qsh_1175 = buffer.data(qsh + 1175);
    const auto *qsh_1176 = buffer.data(qsh + 1176);
    const auto *qsh_1178 = buffer.data(qsh + 1178);
    const auto *qsh_1179 = buffer.data(qsh + 1179);
    const auto *qsh_1181 = buffer.data(qsh + 1181);
    const auto *qsh_1182 = buffer.data(qsh + 1182);
    const auto *qsh_1185 = buffer.data(qsh + 1185);
    const auto *qsh_1190 = buffer.data(qsh + 1190);
    const auto *qsh_1191 = buffer.data(qsh + 1191);
    const auto *qsh_1192 = buffer.data(qsh + 1192);
    const auto *qsh_1193 = buffer.data(qsh + 1193);
    const auto *qsh_1194 = buffer.data(qsh + 1194);
    const auto *qsh_1195 = buffer.data(qsh + 1195);
    const auto *qsh_1196 = buffer.data(qsh + 1196);

#pragma omp simd aligned(t_1478, t_1479, t_1480, pc_y, pc_z, osh_897, osh_920, osh_921, \
                         qsg0_792, qsg0_793, qsg1_792, qsg1_793, qsh_1107, qsh_1109, \
                         qsh_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_21 * osh_897[k]
                    + f_3 * pc_z[k] * qsh_1107[k];

        t_1479[k] = f_12 * osh_920[k]
                    + f_8 * qsg0_792[k]
                    - f_9 * qsg1_792[k]
                    + f_3 * pc_y[k] * qsh_1109[k];

        t_1480[k] = f_12 * osh_921[k]
                    + f_6 * qsg0_793[k]
                    - f_7 * qsg1_793[k]
                    + f_3 * pc_y[k] * qsh_1110[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, t_1484, pa_y, pc_y, pc_z, osi0_1232, osh_902, \
                         osh_922, osh_923, osi1_1232, qsg0_794, qsg1_794, qsh_1111, \
                         qsh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_12 * osh_922[k]
                    + f_4 * qsg0_794[k]
                    - f_5 * qsg1_794[k]
                    + f_3 * pc_y[k] * qsh_1111[k];

        t_1482[k] = f_12 * osh_923[k]
                    + f_3 * pc_y[k] * qsh_1112[k];

        t_1483[k] = f_21 * osh_902[k]
                    + f_1 * qsg0_794[k]
                    - f_2 * qsg1_794[k]
                    + f_3 * pc_z[k] * qsh_1112[k];

        t_1484[k] = pa_y[k] * osi0_1232[k]
                    - f_10 * pc_y[k] * osi1_1232[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, t_1488, pa_y, pc_y, pc_z, osi0_1235, osh_903, \
                         osh_924, osh_925, osh_926, osi1_1235, qsh_1113, \
                         qsh_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = f_11 * osh_924[k]
                    + f_3 * pc_y[k] * qsh_1113[k];

        t_1486[k] = f_20 * osh_903[k]
                    + f_3 * pc_z[k] * qsh_1113[k];

        t_1487[k] = pa_y[k] * osi0_1235[k]
                    + f_12 * osh_925[k]
                    - f_10 * pc_y[k] * osi1_1235[k];

        t_1488[k] = f_11 * osh_926[k]
                    + f_3 * pc_y[k] * qsh_1115[k];
    }

#pragma omp simd aligned(t_1489, t_1490, t_1491, t_1492, pa_y, pc_y, pc_z, osi0_1237, \
                         osi0_1238, osh_906, osh_927, osh_929, osi1_1237, osi1_1238, qsh_1116, \
                         qsh_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1489[k] = pa_y[k] * osi0_1237[k]
                    - f_10 * pc_y[k] * osi1_1237[k];

        t_1490[k] = pa_y[k] * osi0_1238[k]
                    + f_13 * osh_927[k]
                    - f_10 * pc_y[k] * osi1_1238[k];

        t_1491[k] = f_20 * osh_906[k]
                    + f_3 * pc_z[k] * qsh_1116[k];

        t_1492[k] = f_11 * osh_929[k]
                    + f_3 * pc_y[k] * qsh_1118[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pa_y, pc_y, pc_z, osi0_1241, osi0_1242, \
                         osh_909, osh_930, osi1_1241, osi1_1242, \
                         qsh_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = pa_y[k] * osi0_1241[k]
                    - f_10 * pc_y[k] * osi1_1241[k];

        t_1494[k] = pa_y[k] * osi0_1242[k]
                    + f_14 * osh_930[k]
                    - f_10 * pc_y[k] * osi1_1242[k];

        t_1495[k] = f_20 * osh_909[k]
                    + f_3 * pc_z[k] * qsh_1119[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, t_1499, pa_y, pc_x, pc_y, osi0_1244, \
                         osi0_1246, osh_932, osh_933, osh_1128, osi1_1244, osi1_1246, \
                         qsh_1122, qsh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = pa_y[k] * osi0_1244[k]
                    + f_12 * osh_932[k]
                    - f_10 * pc_y[k] * osi1_1244[k];

        t_1497[k] = f_11 * osh_933[k]
                    + f_3 * pc_y[k] * qsh_1122[k];

        t_1498[k] = pa_y[k] * osi0_1246[k]
                    - f_10 * pc_y[k] * osi1_1246[k];

        t_1499[k] = f_13 * osh_1128[k]
                    + f_3 * pc_x[k] * qsh_1128[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, t_1503, t_1504, pc_x, osh_1129, osh_1130, \
                         osh_1131, osh_1132, osh_1133, qsh_1129, qsh_1130, qsh_1131, qsh_1132, \
                         qsh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_13 * osh_1129[k]
                    + f_3 * pc_x[k] * qsh_1129[k];

        t_1501[k] = f_13 * osh_1130[k]
                    + f_3 * pc_x[k] * qsh_1130[k];

        t_1502[k] = f_13 * osh_1131[k]
                    + f_3 * pc_x[k] * qsh_1131[k];

        t_1503[k] = f_13 * osh_1132[k]
                    + f_3 * pc_x[k] * qsh_1132[k];

        t_1504[k] = f_13 * osh_1133[k]
                    + f_3 * pc_x[k] * qsh_1133[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, pc_y, pc_z, osh_918, osh_939, osh_941, \
                         qsg0_805, qsg0_807, qsg1_805, qsg1_807, qsh_1128, \
                         qsh_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = f_11 * osh_939[k]
                    + f_1 * qsg0_805[k]
                    - f_2 * qsg1_805[k]
                    + f_3 * pc_y[k] * qsh_1128[k];

        t_1506[k] = f_20 * osh_918[k]
                    + f_3 * pc_z[k] * qsh_1128[k];

        t_1507[k] = f_11 * osh_941[k]
                    + f_8 * qsg0_807[k]
                    - f_9 * qsg1_807[k]
                    + f_3 * pc_y[k] * qsh_1130[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, pc_y, osh_942, osh_943, osh_944, qsg0_808, \
                         qsg0_809, qsg1_808, qsg1_809, qsh_1131, qsh_1132, \
                         qsh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = f_11 * osh_942[k]
                    + f_6 * qsg0_808[k]
                    - f_7 * qsg1_808[k]
                    + f_3 * pc_y[k] * qsh_1131[k];

        t_1509[k] = f_11 * osh_943[k]
                    + f_4 * qsg0_809[k]
                    - f_5 * qsg1_809[k]
                    + f_3 * pc_y[k] * qsh_1132[k];

        t_1510[k] = f_11 * osh_944[k]
                    + f_3 * pc_y[k] * qsh_1133[k];
    }

#pragma omp simd aligned(t_1511, t_1512, t_1513, t_1514, pa_y, pc_x, pc_y, pc_z, osi0_1259, \
                         osh_924, osh_1134, osi1_1259, qsg0_810, qsg1_810, \
                         qsh_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1511[k] = pa_y[k] * osi0_1259[k]
                    - f_10 * pc_y[k] * osi1_1259[k];

        t_1512[k] = f_13 * osh_1134[k]
                    + f_1 * qsg0_810[k]
                    - f_2 * qsg1_810[k]
                    + f_3 * pc_x[k] * qsh_1134[k];

        t_1513[k] = f_3 * pc_y[k] * qsh_1134[k];

        t_1514[k] = f_19 * osh_924[k]
                    + f_3 * pc_z[k] * qsh_1134[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, pc_x, pc_y, osh_1139, qsg0_810, qsg0_815, \
                         qsg1_810, qsg1_815, qsh_1135, qsh_1136, \
                         qsh_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = f_4 * qsg0_810[k]
                    - f_5 * qsg1_810[k]
                    + f_3 * pc_y[k] * qsh_1135[k];

        t_1516[k] = f_3 * pc_y[k] * qsh_1136[k];

        t_1517[k] = f_13 * osh_1139[k]
                    + f_8 * qsg0_815[k]
                    - f_9 * qsg1_815[k]
                    + f_3 * pc_x[k] * qsh_1139[k];
    }

#pragma omp simd aligned(t_1518, t_1519, t_1520, pc_y, qsg0_811, qsg0_812, qsg1_811, qsg1_812, \
                         qsh_1137, qsh_1138, qsh_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = f_6 * qsg0_811[k]
                    - f_7 * qsg1_811[k]
                    + f_3 * pc_y[k] * qsh_1137[k];

        t_1519[k] = f_4 * qsg0_812[k]
                    - f_5 * qsg1_812[k]
                    + f_3 * pc_y[k] * qsh_1138[k];

        t_1520[k] = f_3 * pc_y[k] * qsh_1139[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, pc_x, pc_y, osh_1143, qsg0_813, qsg0_814, \
                         qsg0_819, qsg1_813, qsg1_814, qsg1_819, qsh_1140, qsh_1141, \
                         qsh_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = f_13 * osh_1143[k]
                    + f_6 * qsg0_819[k]
                    - f_7 * qsg1_819[k]
                    + f_3 * pc_x[k] * qsh_1143[k];

        t_1522[k] = f_8 * qsg0_813[k]
                    - f_9 * qsg1_813[k]
                    + f_3 * pc_y[k] * qsh_1140[k];

        t_1523[k] = f_6 * qsg0_814[k]
                    - f_7 * qsg1_814[k]
                    + f_3 * pc_y[k] * qsh_1141[k];
    }

#pragma omp simd aligned(t_1524, t_1525, t_1526, t_1527, pc_x, pc_y, osh_1148, osh_1149, \
                         qsg0_815, qsg0_824, qsg1_815, qsg1_824, qsh_1142, qsh_1143, qsh_1148, \
                         qsh_1149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1524[k] = f_4 * qsg0_815[k]
                    - f_5 * qsg1_815[k]
                    + f_3 * pc_y[k] * qsh_1142[k];

        t_1525[k] = f_3 * pc_y[k] * qsh_1143[k];

        t_1526[k] = f_13 * osh_1148[k]
                    + f_4 * qsg0_824[k]
                    - f_5 * qsg1_824[k]
                    + f_3 * pc_x[k] * qsh_1148[k];

        t_1527[k] = f_13 * osh_1149[k]
                    + f_3 * pc_x[k] * qsh_1149[k];
    }

#pragma omp simd aligned(t_1528, t_1529, t_1530, t_1531, t_1532, pc_x, pc_y, osh_1150, \
                         osh_1151, osh_1152, osh_1154, qsh_1148, qsh_1150, qsh_1151, qsh_1152, \
                         qsh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1528[k] = f_13 * osh_1150[k]
                    + f_3 * pc_x[k] * qsh_1150[k];

        t_1529[k] = f_13 * osh_1151[k]
                    + f_3 * pc_x[k] * qsh_1151[k];

        t_1530[k] = f_13 * osh_1152[k]
                    + f_3 * pc_x[k] * qsh_1152[k];

        t_1531[k] = f_3 * pc_y[k] * qsh_1148[k];

        t_1532[k] = f_13 * osh_1154[k]
                    + f_3 * pc_x[k] * qsh_1154[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, pc_y, qsg0_820, qsg0_821, qsg0_822, qsg1_820, \
                         qsg1_821, qsg1_822, qsh_1149, qsh_1150, \
                         qsh_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = f_1 * qsg0_820[k]
                    - f_2 * qsg1_820[k]
                    + f_3 * pc_y[k] * qsh_1149[k];

        t_1534[k] = f_16 * qsg0_821[k]
                    - f_17 * qsg1_821[k]
                    + f_3 * pc_y[k] * qsh_1150[k];

        t_1535[k] = f_8 * qsg0_822[k]
                    - f_9 * qsg1_822[k]
                    + f_3 * pc_y[k] * qsh_1151[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, t_1539, pc_y, pc_z, osh_944, qsg0_823, \
                         qsg0_824, qsg1_823, qsg1_824, qsh_1152, qsh_1153, \
                         qsh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_6 * qsg0_823[k]
                    - f_7 * qsg1_823[k]
                    + f_3 * pc_y[k] * qsh_1152[k];

        t_1537[k] = f_4 * qsg0_824[k]
                    - f_5 * qsg1_824[k]
                    + f_3 * pc_y[k] * qsh_1153[k];

        t_1538[k] = f_3 * pc_y[k] * qsh_1154[k];

        t_1539[k] = f_19 * osh_944[k]
                    + f_1 * qsg0_824[k]
                    - f_2 * qsg1_824[k]
                    + f_3 * pc_z[k] * qsh_1154[k];
    }

#pragma omp simd aligned(t_1540, t_1541, t_1542, t_1543, pc_x, pc_y, pc_z, osh_945, osh_1155, \
                         osh_1158, qsg0_825, qsg0_828, qsg1_825, qsg1_828, qsh_1155, \
                         qsh_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1540[k] = f_12 * osh_1155[k]
                    + f_1 * qsg0_825[k]
                    - f_2 * qsg1_825[k]
                    + f_3 * pc_x[k] * qsh_1155[k];

        t_1541[k] = f_18 * osh_945[k]
                    + f_3 * pc_y[k] * qsh_1155[k];

        t_1542[k] = f_3 * pc_z[k] * qsh_1155[k];

        t_1543[k] = f_12 * osh_1158[k]
                    + f_8 * qsg0_828[k]
                    - f_9 * qsg1_828[k]
                    + f_3 * pc_x[k] * qsh_1158[k];
    }

#pragma omp simd aligned(t_1544, t_1545, t_1546, t_1547, pc_x, pc_z, osh_1161, qsg0_825, \
                         qsg0_831, qsg1_825, qsg1_831, qsh_1156, qsh_1157, qsh_1158, \
                         qsh_1161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1544[k] = f_3 * pc_z[k] * qsh_1156[k];

        t_1545[k] = f_4 * qsg0_825[k]
                    - f_5 * qsg1_825[k]
                    + f_3 * pc_z[k] * qsh_1157[k];

        t_1546[k] = f_12 * osh_1161[k]
                    + f_6 * qsg0_831[k]
                    - f_7 * qsg1_831[k]
                    + f_3 * pc_x[k] * qsh_1161[k];

        t_1547[k] = f_3 * pc_z[k] * qsh_1158[k];
    }

#pragma omp simd aligned(t_1548, t_1549, t_1550, t_1551, pc_x, pc_y, pc_z, osh_950, osh_1165, \
                         qsg0_827, qsg0_835, qsg1_827, qsg1_835, qsh_1160, qsh_1161, \
                         qsh_1165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1548[k] = f_18 * osh_950[k]
                    + f_3 * pc_y[k] * qsh_1160[k];

        t_1549[k] = f_6 * qsg0_827[k]
                    - f_7 * qsg1_827[k]
                    + f_3 * pc_z[k] * qsh_1160[k];

        t_1550[k] = f_12 * osh_1165[k]
                    + f_4 * qsg0_835[k]
                    - f_5 * qsg1_835[k]
                    + f_3 * pc_x[k] * qsh_1165[k];

        t_1551[k] = f_3 * pc_z[k] * qsh_1161[k];
    }

#pragma omp simd aligned(t_1552, t_1553, t_1554, t_1555, pc_x, pc_y, pc_z, osh_954, osh_1170, \
                         qsg0_828, qsg0_830, qsg1_828, qsg1_830, qsh_1162, qsh_1164, \
                         qsh_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1552[k] = f_4 * qsg0_828[k]
                    - f_5 * qsg1_828[k]
                    + f_3 * pc_z[k] * qsh_1162[k];

        t_1553[k] = f_18 * osh_954[k]
                    + f_3 * pc_y[k] * qsh_1164[k];

        t_1554[k] = f_8 * qsg0_830[k]
                    - f_9 * qsg1_830[k]
                    + f_3 * pc_z[k] * qsh_1164[k];

        t_1555[k] = f_12 * osh_1170[k]
                    + f_3 * pc_x[k] * qsh_1170[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, t_1559, t_1560, pc_x, pc_z, osh_1172, \
                         osh_1173, osh_1174, osh_1175, qsh_1165, qsh_1172, qsh_1173, qsh_1174, \
                         qsh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_3 * pc_z[k] * qsh_1165[k];

        t_1557[k] = f_12 * osh_1172[k]
                    + f_3 * pc_x[k] * qsh_1172[k];

        t_1558[k] = f_12 * osh_1173[k]
                    + f_3 * pc_x[k] * qsh_1173[k];

        t_1559[k] = f_12 * osh_1174[k]
                    + f_3 * pc_x[k] * qsh_1174[k];

        t_1560[k] = f_12 * osh_1175[k]
                    + f_3 * pc_x[k] * qsh_1175[k];
    }

#pragma omp simd aligned(t_1561, t_1562, t_1563, t_1564, pc_y, pc_z, osh_960, qsg0_835, \
                         qsg0_836, qsg1_835, qsg1_836, qsh_1170, qsh_1171, \
                         qsh_1172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1561[k] = f_18 * osh_960[k]
                    + f_1 * qsg0_835[k]
                    - f_2 * qsg1_835[k]
                    + f_3 * pc_y[k] * qsh_1170[k];

        t_1562[k] = f_3 * pc_z[k] * qsh_1170[k];

        t_1563[k] = f_4 * qsg0_835[k]
                    - f_5 * qsg1_835[k]
                    + f_3 * pc_z[k] * qsh_1171[k];

        t_1564[k] = f_6 * qsg0_836[k]
                    - f_7 * qsg1_836[k]
                    + f_3 * pc_z[k] * qsh_1172[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, t_1568, pa_z, pc_y, pc_z, osi0_1260, osh_965, \
                         osi1_1260, qsg0_837, qsg0_839, qsg1_837, qsg1_839, qsh_1173, \
                         qsh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = f_8 * qsg0_837[k]
                    - f_9 * qsg1_837[k]
                    + f_3 * pc_z[k] * qsh_1173[k];

        t_1566[k] = f_18 * osh_965[k]
                    + f_3 * pc_y[k] * qsh_1175[k];

        t_1567[k] = f_1 * qsg0_839[k]
                    - f_2 * qsg1_839[k]
                    + f_3 * pc_z[k] * qsh_1175[k];

        t_1568[k] = pa_z[k] * osi0_1260[k]
                    - f_10 * pc_z[k] * osi1_1260[k];
    }

#pragma omp simd aligned(t_1569, t_1570, t_1571, t_1572, pa_z, pc_y, pc_z, osi0_1263, osh_945, \
                         osh_966, osh_968, osi1_1263, qsh_1176, \
                         qsh_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1569[k] = f_19 * osh_966[k]
                    + f_3 * pc_y[k] * qsh_1176[k];

        t_1570[k] = f_11 * osh_945[k]
                    + f_3 * pc_z[k] * qsh_1176[k];

        t_1571[k] = pa_z[k] * osi0_1263[k]
                    - f_10 * pc_z[k] * osi1_1263[k];

        t_1572[k] = f_19 * osh_968[k]
                    + f_3 * pc_y[k] * qsh_1178[k];
    }

#pragma omp simd aligned(t_1573, t_1574, t_1575, pa_z, pc_x, pc_z, osi0_1266, osh_948, \
                         osh_1181, osi1_1266, qsg0_845, qsg1_845, qsh_1179, \
                         qsh_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1573[k] = f_12 * osh_1181[k]
                    + f_8 * qsg0_845[k]
                    - f_9 * qsg1_845[k]
                    + f_3 * pc_x[k] * qsh_1181[k];

        t_1574[k] = pa_z[k] * osi0_1266[k]
                    - f_10 * pc_z[k] * osi1_1266[k];

        t_1575[k] = f_11 * osh_948[k]
                    + f_3 * pc_z[k] * qsh_1179[k];
    }

#pragma omp simd aligned(t_1576, t_1577, t_1578, pa_z, pc_x, pc_y, pc_z, osi0_1270, osh_971, \
                         osh_1185, osi1_1270, qsg0_849, qsg1_849, qsh_1181, \
                         qsh_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1576[k] = f_19 * osh_971[k]
                    + f_3 * pc_y[k] * qsh_1181[k];

        t_1577[k] = f_12 * osh_1185[k]
                    + f_6 * qsg0_849[k]
                    - f_7 * qsg1_849[k]
                    + f_3 * pc_x[k] * qsh_1185[k];

        t_1578[k] = pa_z[k] * osi0_1270[k]
                    - f_10 * pc_z[k] * osi1_1270[k];
    }

#pragma omp simd aligned(t_1579, t_1580, t_1581, pa_z, pc_y, pc_z, osi0_1272, osh_951, \
                         osh_952, osh_975, osi1_1272, qsh_1182, \
                         qsh_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1579[k] = f_11 * osh_951[k]
                    + f_3 * pc_z[k] * qsh_1182[k];

        t_1580[k] = pa_z[k] * osi0_1272[k]
                    + f_12 * osh_952[k]
                    - f_10 * pc_z[k] * osi1_1272[k];

        t_1581[k] = f_19 * osh_975[k]
                    + f_3 * pc_y[k] * qsh_1185[k];
    }

#pragma omp simd aligned(t_1582, t_1583, t_1584, t_1585, pc_x, osh_1190, osh_1191, osh_1192, \
                         osh_1193, qsg0_854, qsg1_854, qsh_1190, qsh_1191, qsh_1192, \
                         qsh_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1582[k] = f_12 * osh_1190[k]
                    + f_4 * qsg0_854[k]
                    - f_5 * qsg1_854[k]
                    + f_3 * pc_x[k] * qsh_1190[k];

        t_1583[k] = f_12 * osh_1191[k]
                    + f_3 * pc_x[k] * qsh_1191[k];

        t_1584[k] = f_12 * osh_1192[k]
                    + f_3 * pc_x[k] * qsh_1192[k];

        t_1585[k] = f_12 * osh_1193[k]
                    + f_3 * pc_x[k] * qsh_1193[k];
    }

#pragma omp simd aligned(t_1586, t_1587, t_1588, t_1589, pa_z, pc_x, pc_z, osi0_1281, \
                         osh_1194, osh_1195, osh_1196, osi1_1281, qsh_1194, qsh_1195, \
                         qsh_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1586[k] = f_12 * osh_1194[k]
                    + f_3 * pc_x[k] * qsh_1194[k];

        t_1587[k] = f_12 * osh_1195[k]
                    + f_3 * pc_x[k] * qsh_1195[k];

        t_1588[k] = f_12 * osh_1196[k]
                    + f_3 * pc_x[k] * qsh_1196[k];

        t_1589[k] = pa_z[k] * osi0_1281[k]
                    - f_10 * pc_z[k] * osi1_1281[k];
    }

#pragma omp simd aligned(t_1590, t_1591, t_1592, pc_y, pc_z, osh_960, osh_983, osh_984, \
                         qsg0_852, qsg0_853, qsg1_852, qsg1_853, qsh_1191, qsh_1193, \
                         qsh_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1590[k] = f_11 * osh_960[k]
                    + f_3 * pc_z[k] * qsh_1191[k];

        t_1591[k] = f_19 * osh_983[k]
                    + f_8 * qsg0_852[k]
                    - f_9 * qsg1_852[k]
                    + f_3 * pc_y[k] * qsh_1193[k];

        t_1592[k] = f_19 * osh_984[k]
                    + f_6 * qsg0_853[k]
                    - f_7 * qsg1_853[k]
                    + f_3 * pc_y[k] * qsh_1194[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, pc_y, pc_z, osh_965, osh_985, osh_986, \
                         qsg0_854, qsg1_854, qsh_1195, qsh_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = f_19 * osh_985[k]
                    + f_4 * qsg0_854[k]
                    - f_5 * qsg1_854[k]
                    + f_3 * pc_y[k] * qsh_1195[k];

        t_1594[k] = f_19 * osh_986[k]
                    + f_3 * pc_y[k] * qsh_1196[k];

        t_1595[k] = f_11 * osh_965[k]
                    + f_1 * qsg0_854[k]
                    - f_2 * qsg1_854[k]
                    + f_3 * pc_z[k] * qsh_1196[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t osh, const size_t qsg0,
                                                           const size_t qsg1, const size_t qsh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_1596 = buffer.data(target + 1596);
    auto *t_1597 = buffer.data(target + 1597);
    auto *t_1598 = buffer.data(target + 1598);
    auto *t_1599 = buffer.data(target + 1599);
    auto *t_1600 = buffer.data(target + 1600);
    auto *t_1601 = buffer.data(target + 1601);
    auto *t_1602 = buffer.data(target + 1602);
    auto *t_1603 = buffer.data(target + 1603);
    auto *t_1604 = buffer.data(target + 1604);
    auto *t_1605 = buffer.data(target + 1605);
    auto *t_1606 = buffer.data(target + 1606);
    auto *t_1607 = buffer.data(target + 1607);
    auto *t_1608 = buffer.data(target + 1608);
    auto *t_1609 = buffer.data(target + 1609);
    auto *t_1610 = buffer.data(target + 1610);
    auto *t_1611 = buffer.data(target + 1611);
    auto *t_1612 = buffer.data(target + 1612);
    auto *t_1613 = buffer.data(target + 1613);
    auto *t_1614 = buffer.data(target + 1614);
    auto *t_1615 = buffer.data(target + 1615);
    auto *t_1616 = buffer.data(target + 1616);
    auto *t_1617 = buffer.data(target + 1617);
    auto *t_1618 = buffer.data(target + 1618);
    auto *t_1619 = buffer.data(target + 1619);
    auto *t_1620 = buffer.data(target + 1620);
    auto *t_1621 = buffer.data(target + 1621);
    auto *t_1622 = buffer.data(target + 1622);
    auto *t_1623 = buffer.data(target + 1623);
    auto *t_1624 = buffer.data(target + 1624);
    auto *t_1625 = buffer.data(target + 1625);
    auto *t_1626 = buffer.data(target + 1626);
    auto *t_1627 = buffer.data(target + 1627);
    auto *t_1628 = buffer.data(target + 1628);
    auto *t_1629 = buffer.data(target + 1629);
    auto *t_1630 = buffer.data(target + 1630);
    auto *t_1631 = buffer.data(target + 1631);
    auto *t_1632 = buffer.data(target + 1632);
    auto *t_1633 = buffer.data(target + 1633);
    auto *t_1634 = buffer.data(target + 1634);
    auto *t_1635 = buffer.data(target + 1635);
    auto *t_1636 = buffer.data(target + 1636);
    auto *t_1637 = buffer.data(target + 1637);
    auto *t_1638 = buffer.data(target + 1638);
    auto *t_1639 = buffer.data(target + 1639);
    auto *t_1640 = buffer.data(target + 1640);
    auto *t_1641 = buffer.data(target + 1641);
    auto *t_1642 = buffer.data(target + 1642);
    auto *t_1643 = buffer.data(target + 1643);
    auto *t_1644 = buffer.data(target + 1644);
    auto *t_1645 = buffer.data(target + 1645);
    auto *t_1646 = buffer.data(target + 1646);
    auto *t_1647 = buffer.data(target + 1647);
    auto *t_1648 = buffer.data(target + 1648);
    auto *t_1649 = buffer.data(target + 1649);
    auto *t_1650 = buffer.data(target + 1650);
    auto *t_1651 = buffer.data(target + 1651);
    auto *t_1652 = buffer.data(target + 1652);
    auto *t_1653 = buffer.data(target + 1653);
    auto *t_1654 = buffer.data(target + 1654);
    auto *t_1655 = buffer.data(target + 1655);
    auto *t_1656 = buffer.data(target + 1656);
    auto *t_1657 = buffer.data(target + 1657);
    auto *t_1658 = buffer.data(target + 1658);
    auto *t_1659 = buffer.data(target + 1659);
    auto *t_1660 = buffer.data(target + 1660);
    auto *t_1661 = buffer.data(target + 1661);
    auto *t_1662 = buffer.data(target + 1662);
    auto *t_1663 = buffer.data(target + 1663);
    auto *t_1664 = buffer.data(target + 1664);
    auto *t_1665 = buffer.data(target + 1665);
    auto *t_1666 = buffer.data(target + 1666);
    auto *t_1667 = buffer.data(target + 1667);
    auto *t_1668 = buffer.data(target + 1668);
    auto *t_1669 = buffer.data(target + 1669);
    auto *t_1670 = buffer.data(target + 1670);
    auto *t_1671 = buffer.data(target + 1671);
    auto *t_1672 = buffer.data(target + 1672);
    auto *t_1673 = buffer.data(target + 1673);
    auto *t_1674 = buffer.data(target + 1674);
    auto *t_1675 = buffer.data(target + 1675);
    auto *t_1676 = buffer.data(target + 1676);
    auto *t_1677 = buffer.data(target + 1677);
    auto *t_1678 = buffer.data(target + 1678);
    auto *t_1679 = buffer.data(target + 1679);
    auto *t_1680 = buffer.data(target + 1680);
    auto *t_1681 = buffer.data(target + 1681);
    auto *t_1682 = buffer.data(target + 1682);
    auto *t_1683 = buffer.data(target + 1683);
    auto *t_1684 = buffer.data(target + 1684);
    auto *t_1685 = buffer.data(target + 1685);
    auto *t_1686 = buffer.data(target + 1686);
    auto *t_1687 = buffer.data(target + 1687);
    auto *t_1688 = buffer.data(target + 1688);
    auto *t_1689 = buffer.data(target + 1689);
    auto *t_1690 = buffer.data(target + 1690);
    auto *t_1691 = buffer.data(target + 1691);
    auto *t_1692 = buffer.data(target + 1692);
    auto *t_1693 = buffer.data(target + 1693);
    auto *t_1694 = buffer.data(target + 1694);
    auto *t_1695 = buffer.data(target + 1695);
    auto *t_1696 = buffer.data(target + 1696);
    auto *t_1697 = buffer.data(target + 1697);
    auto *t_1698 = buffer.data(target + 1698);
    auto *t_1699 = buffer.data(target + 1699);
    auto *t_1700 = buffer.data(target + 1700);
    auto *t_1701 = buffer.data(target + 1701);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh_966 = buffer.data(osh + 966);
    const auto *osh_969 = buffer.data(osh + 969);
    const auto *osh_972 = buffer.data(osh + 972);
    const auto *osh_981 = buffer.data(osh + 981);
    const auto *osh_986 = buffer.data(osh + 986);
    const auto *osh_987 = buffer.data(osh + 987);
    const auto *osh_989 = buffer.data(osh + 989);
    const auto *osh_990 = buffer.data(osh + 990);
    const auto *osh_992 = buffer.data(osh + 992);
    const auto *osh_993 = buffer.data(osh + 993);
    const auto *osh_996 = buffer.data(osh + 996);
    const auto *osh_1002 = buffer.data(osh + 1002);
    const auto *osh_1004 = buffer.data(osh + 1004);
    const auto *osh_1005 = buffer.data(osh + 1005);
    const auto *osh_1006 = buffer.data(osh + 1006);
    const auto *osh_1007 = buffer.data(osh + 1007);
    const auto *osh_1008 = buffer.data(osh + 1008);
    const auto *osh_1010 = buffer.data(osh + 1010);
    const auto *osh_1011 = buffer.data(osh + 1011);
    const auto *osh_1013 = buffer.data(osh + 1013);
    const auto *osh_1014 = buffer.data(osh + 1014);
    const auto *osh_1017 = buffer.data(osh + 1017);
    const auto *osh_1023 = buffer.data(osh + 1023);
    const auto *osh_1025 = buffer.data(osh + 1025);
    const auto *osh_1026 = buffer.data(osh + 1026);
    const auto *osh_1027 = buffer.data(osh + 1027);
    const auto *osh_1028 = buffer.data(osh + 1028);
    const auto *osh_1029 = buffer.data(osh + 1029);
    const auto *osh_1031 = buffer.data(osh + 1031);
    const auto *osh_1032 = buffer.data(osh + 1032);
    const auto *osh_1034 = buffer.data(osh + 1034);
    const auto *osh_1035 = buffer.data(osh + 1035);
    const auto *osh_1038 = buffer.data(osh + 1038);
    const auto *osh_1044 = buffer.data(osh + 1044);
    const auto *osh_1046 = buffer.data(osh + 1046);
    const auto *osh_1047 = buffer.data(osh + 1047);
    const auto *osh_1048 = buffer.data(osh + 1048);
    const auto *osh_1049 = buffer.data(osh + 1049);
    const auto *osh_1050 = buffer.data(osh + 1050);
    const auto *osh_1052 = buffer.data(osh + 1052);
    const auto *osh_1055 = buffer.data(osh + 1055);
    const auto *osh_1059 = buffer.data(osh + 1059);
    const auto *osh_1065 = buffer.data(osh + 1065);
    const auto *osh_1197 = buffer.data(osh + 1197);
    const auto *osh_1200 = buffer.data(osh + 1200);
    const auto *osh_1202 = buffer.data(osh + 1202);
    const auto *osh_1203 = buffer.data(osh + 1203);
    const auto *osh_1206 = buffer.data(osh + 1206);
    const auto *osh_1207 = buffer.data(osh + 1207);
    const auto *osh_1209 = buffer.data(osh + 1209);
    const auto *osh_1211 = buffer.data(osh + 1211);
    const auto *osh_1212 = buffer.data(osh + 1212);
    const auto *osh_1213 = buffer.data(osh + 1213);
    const auto *osh_1214 = buffer.data(osh + 1214);
    const auto *osh_1215 = buffer.data(osh + 1215);
    const auto *osh_1216 = buffer.data(osh + 1216);
    const auto *osh_1217 = buffer.data(osh + 1217);
    const auto *osh_1218 = buffer.data(osh + 1218);
    const auto *osh_1221 = buffer.data(osh + 1221);
    const auto *osh_1223 = buffer.data(osh + 1223);
    const auto *osh_1224 = buffer.data(osh + 1224);
    const auto *osh_1227 = buffer.data(osh + 1227);
    const auto *osh_1228 = buffer.data(osh + 1228);
    const auto *osh_1230 = buffer.data(osh + 1230);
    const auto *osh_1232 = buffer.data(osh + 1232);
    const auto *osh_1233 = buffer.data(osh + 1233);
    const auto *osh_1234 = buffer.data(osh + 1234);
    const auto *osh_1235 = buffer.data(osh + 1235);
    const auto *osh_1236 = buffer.data(osh + 1236);
    const auto *osh_1237 = buffer.data(osh + 1237);
    const auto *osh_1238 = buffer.data(osh + 1238);
    const auto *osh_1239 = buffer.data(osh + 1239);
    const auto *osh_1242 = buffer.data(osh + 1242);
    const auto *osh_1244 = buffer.data(osh + 1244);
    const auto *osh_1245 = buffer.data(osh + 1245);
    const auto *osh_1248 = buffer.data(osh + 1248);
    const auto *osh_1249 = buffer.data(osh + 1249);
    const auto *osh_1251 = buffer.data(osh + 1251);
    const auto *osh_1253 = buffer.data(osh + 1253);
    const auto *osh_1254 = buffer.data(osh + 1254);
    const auto *osh_1255 = buffer.data(osh + 1255);
    const auto *osh_1256 = buffer.data(osh + 1256);
    const auto *osh_1257 = buffer.data(osh + 1257);
    const auto *osh_1258 = buffer.data(osh + 1258);
    const auto *osh_1259 = buffer.data(osh + 1259);
    const auto *osh_1260 = buffer.data(osh + 1260);
    const auto *osh_1263 = buffer.data(osh + 1263);
    const auto *osh_1265 = buffer.data(osh + 1265);
    const auto *osh_1266 = buffer.data(osh + 1266);
    const auto *osh_1269 = buffer.data(osh + 1269);
    const auto *osh_1270 = buffer.data(osh + 1270);
    const auto *osh_1272 = buffer.data(osh + 1272);
    const auto *osh_1274 = buffer.data(osh + 1274);
    const auto *osh_1275 = buffer.data(osh + 1275);
    const auto *osh_1276 = buffer.data(osh + 1276);
    const auto *osh_1277 = buffer.data(osh + 1277);
    const auto *osh_1278 = buffer.data(osh + 1278);
    const auto *osh_1279 = buffer.data(osh + 1279);
    const auto *osh_1280 = buffer.data(osh + 1280);

    const auto *qsg0_855 = buffer.data(qsg0 + 855);
    const auto *qsg0_858 = buffer.data(qsg0 + 858);
    const auto *qsg0_860 = buffer.data(qsg0 + 860);
    const auto *qsg0_861 = buffer.data(qsg0 + 861);
    const auto *qsg0_864 = buffer.data(qsg0 + 864);
    const auto *qsg0_865 = buffer.data(qsg0 + 865);
    const auto *qsg0_867 = buffer.data(qsg0 + 867);
    const auto *qsg0_868 = buffer.data(qsg0 + 868);
    const auto *qsg0_869 = buffer.data(qsg0 + 869);
    const auto *qsg0_870 = buffer.data(qsg0 + 870);
    const auto *qsg0_873 = buffer.data(qsg0 + 873);
    const auto *qsg0_875 = buffer.data(qsg0 + 875);
    const auto *qsg0_876 = buffer.data(qsg0 + 876);
    const auto *qsg0_879 = buffer.data(qsg0 + 879);
    const auto *qsg0_880 = buffer.data(qsg0 + 880);
    const auto *qsg0_882 = buffer.data(qsg0 + 882);
    const auto *qsg0_883 = buffer.data(qsg0 + 883);
    const auto *qsg0_884 = buffer.data(qsg0 + 884);
    const auto *qsg0_885 = buffer.data(qsg0 + 885);
    const auto *qsg0_888 = buffer.data(qsg0 + 888);
    const auto *qsg0_890 = buffer.data(qsg0 + 890);
    const auto *qsg0_891 = buffer.data(qsg0 + 891);
    const auto *qsg0_894 = buffer.data(qsg0 + 894);
    const auto *qsg0_895 = buffer.data(qsg0 + 895);
    const auto *qsg0_897 = buffer.data(qsg0 + 897);
    const auto *qsg0_898 = buffer.data(qsg0 + 898);
    const auto *qsg0_899 = buffer.data(qsg0 + 899);
    const auto *qsg0_900 = buffer.data(qsg0 + 900);
    const auto *qsg0_903 = buffer.data(qsg0 + 903);
    const auto *qsg0_905 = buffer.data(qsg0 + 905);
    const auto *qsg0_906 = buffer.data(qsg0 + 906);
    const auto *qsg0_909 = buffer.data(qsg0 + 909);
    const auto *qsg0_910 = buffer.data(qsg0 + 910);
    const auto *qsg0_912 = buffer.data(qsg0 + 912);
    const auto *qsg0_914 = buffer.data(qsg0 + 914);

    const auto *qsg1_855 = buffer.data(qsg1 + 855);
    const auto *qsg1_858 = buffer.data(qsg1 + 858);
    const auto *qsg1_860 = buffer.data(qsg1 + 860);
    const auto *qsg1_861 = buffer.data(qsg1 + 861);
    const auto *qsg1_864 = buffer.data(qsg1 + 864);
    const auto *qsg1_865 = buffer.data(qsg1 + 865);
    const auto *qsg1_867 = buffer.data(qsg1 + 867);
    const auto *qsg1_868 = buffer.data(qsg1 + 868);
    const auto *qsg1_869 = buffer.data(qsg1 + 869);
    const auto *qsg1_870 = buffer.data(qsg1 + 870);
    const auto *qsg1_873 = buffer.data(qsg1 + 873);
    const auto *qsg1_875 = buffer.data(qsg1 + 875);
    const auto *qsg1_876 = buffer.data(qsg1 + 876);
    const auto *qsg1_879 = buffer.data(qsg1 + 879);
    const auto *qsg1_880 = buffer.data(qsg1 + 880);
    const auto *qsg1_882 = buffer.data(qsg1 + 882);
    const auto *qsg1_883 = buffer.data(qsg1 + 883);
    const auto *qsg1_884 = buffer.data(qsg1 + 884);
    const auto *qsg1_885 = buffer.data(qsg1 + 885);
    const auto *qsg1_888 = buffer.data(qsg1 + 888);
    const auto *qsg1_890 = buffer.data(qsg1 + 890);
    const auto *qsg1_891 = buffer.data(qsg1 + 891);
    const auto *qsg1_894 = buffer.data(qsg1 + 894);
    const auto *qsg1_895 = buffer.data(qsg1 + 895);
    const auto *qsg1_897 = buffer.data(qsg1 + 897);
    const auto *qsg1_898 = buffer.data(qsg1 + 898);
    const auto *qsg1_899 = buffer.data(qsg1 + 899);
    const auto *qsg1_900 = buffer.data(qsg1 + 900);
    const auto *qsg1_903 = buffer.data(qsg1 + 903);
    const auto *qsg1_905 = buffer.data(qsg1 + 905);
    const auto *qsg1_906 = buffer.data(qsg1 + 906);
    const auto *qsg1_909 = buffer.data(qsg1 + 909);
    const auto *qsg1_910 = buffer.data(qsg1 + 910);
    const auto *qsg1_912 = buffer.data(qsg1 + 912);
    const auto *qsg1_914 = buffer.data(qsg1 + 914);

    const auto *qsh_1197 = buffer.data(qsh + 1197);
    const auto *qsh_1199 = buffer.data(qsh + 1199);
    const auto *qsh_1200 = buffer.data(qsh + 1200);
    const auto *qsh_1202 = buffer.data(qsh + 1202);
    const auto *qsh_1203 = buffer.data(qsh + 1203);
    const auto *qsh_1206 = buffer.data(qsh + 1206);
    const auto *qsh_1207 = buffer.data(qsh + 1207);
    const auto *qsh_1209 = buffer.data(qsh + 1209);
    const auto *qsh_1211 = buffer.data(qsh + 1211);
    const auto *qsh_1212 = buffer.data(qsh + 1212);
    const auto *qsh_1213 = buffer.data(qsh + 1213);
    const auto *qsh_1214 = buffer.data(qsh + 1214);
    const auto *qsh_1215 = buffer.data(qsh + 1215);
    const auto *qsh_1216 = buffer.data(qsh + 1216);
    const auto *qsh_1217 = buffer.data(qsh + 1217);
    const auto *qsh_1218 = buffer.data(qsh + 1218);
    const auto *qsh_1220 = buffer.data(qsh + 1220);
    const auto *qsh_1221 = buffer.data(qsh + 1221);
    const auto *qsh_1223 = buffer.data(qsh + 1223);
    const auto *qsh_1224 = buffer.data(qsh + 1224);
    const auto *qsh_1227 = buffer.data(qsh + 1227);
    const auto *qsh_1228 = buffer.data(qsh + 1228);
    const auto *qsh_1230 = buffer.data(qsh + 1230);
    const auto *qsh_1232 = buffer.data(qsh + 1232);
    const auto *qsh_1233 = buffer.data(qsh + 1233);
    const auto *qsh_1234 = buffer.data(qsh + 1234);
    const auto *qsh_1235 = buffer.data(qsh + 1235);
    const auto *qsh_1236 = buffer.data(qsh + 1236);
    const auto *qsh_1237 = buffer.data(qsh + 1237);
    const auto *qsh_1238 = buffer.data(qsh + 1238);
    const auto *qsh_1239 = buffer.data(qsh + 1239);
    const auto *qsh_1241 = buffer.data(qsh + 1241);
    const auto *qsh_1242 = buffer.data(qsh + 1242);
    const auto *qsh_1244 = buffer.data(qsh + 1244);
    const auto *qsh_1245 = buffer.data(qsh + 1245);
    const auto *qsh_1248 = buffer.data(qsh + 1248);
    const auto *qsh_1249 = buffer.data(qsh + 1249);
    const auto *qsh_1251 = buffer.data(qsh + 1251);
    const auto *qsh_1253 = buffer.data(qsh + 1253);
    const auto *qsh_1254 = buffer.data(qsh + 1254);
    const auto *qsh_1255 = buffer.data(qsh + 1255);
    const auto *qsh_1256 = buffer.data(qsh + 1256);
    const auto *qsh_1257 = buffer.data(qsh + 1257);
    const auto *qsh_1258 = buffer.data(qsh + 1258);
    const auto *qsh_1259 = buffer.data(qsh + 1259);
    const auto *qsh_1260 = buffer.data(qsh + 1260);
    const auto *qsh_1262 = buffer.data(qsh + 1262);
    const auto *qsh_1263 = buffer.data(qsh + 1263);
    const auto *qsh_1265 = buffer.data(qsh + 1265);
    const auto *qsh_1266 = buffer.data(qsh + 1266);
    const auto *qsh_1269 = buffer.data(qsh + 1269);
    const auto *qsh_1270 = buffer.data(qsh + 1270);
    const auto *qsh_1272 = buffer.data(qsh + 1272);
    const auto *qsh_1274 = buffer.data(qsh + 1274);
    const auto *qsh_1275 = buffer.data(qsh + 1275);
    const auto *qsh_1276 = buffer.data(qsh + 1276);
    const auto *qsh_1277 = buffer.data(qsh + 1277);
    const auto *qsh_1278 = buffer.data(qsh + 1278);
    const auto *qsh_1279 = buffer.data(qsh + 1279);
    const auto *qsh_1280 = buffer.data(qsh + 1280);

#pragma omp simd aligned(t_1596, t_1597, t_1598, pc_x, pc_y, pc_z, osh_966, osh_987, osh_1197, \
                         qsg0_855, qsg1_855, qsh_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = f_12 * osh_1197[k]
                    + f_1 * qsg0_855[k]
                    - f_2 * qsg1_855[k]
                    + f_3 * pc_x[k] * qsh_1197[k];

        t_1597[k] = f_20 * osh_987[k]
                    + f_3 * pc_y[k] * qsh_1197[k];

        t_1598[k] = f_12 * osh_966[k]
                    + f_3 * pc_z[k] * qsh_1197[k];
    }

#pragma omp simd aligned(t_1599, t_1600, t_1601, pc_x, pc_y, osh_989, osh_1200, osh_1202, \
                         qsg0_858, qsg0_860, qsg1_858, qsg1_860, qsh_1199, qsh_1200, \
                         qsh_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1599[k] = f_12 * osh_1200[k]
                    + f_8 * qsg0_858[k]
                    - f_9 * qsg1_858[k]
                    + f_3 * pc_x[k] * qsh_1200[k];

        t_1600[k] = f_20 * osh_989[k]
                    + f_3 * pc_y[k] * qsh_1199[k];

        t_1601[k] = f_12 * osh_1202[k]
                    + f_8 * qsg0_860[k]
                    - f_9 * qsg1_860[k]
                    + f_3 * pc_x[k] * qsh_1202[k];
    }

#pragma omp simd aligned(t_1602, t_1603, t_1604, pc_x, pc_y, pc_z, osh_969, osh_992, osh_1203, \
                         qsg0_861, qsg1_861, qsh_1200, qsh_1202, \
                         qsh_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1602[k] = f_12 * osh_1203[k]
                    + f_6 * qsg0_861[k]
                    - f_7 * qsg1_861[k]
                    + f_3 * pc_x[k] * qsh_1203[k];

        t_1603[k] = f_12 * osh_969[k]
                    + f_3 * pc_z[k] * qsh_1200[k];

        t_1604[k] = f_20 * osh_992[k]
                    + f_3 * pc_y[k] * qsh_1202[k];
    }

#pragma omp simd aligned(t_1605, t_1606, t_1607, pc_x, pc_z, osh_972, osh_1206, osh_1207, \
                         qsg0_864, qsg0_865, qsg1_864, qsg1_865, qsh_1203, qsh_1206, \
                         qsh_1207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1605[k] = f_12 * osh_1206[k]
                    + f_6 * qsg0_864[k]
                    - f_7 * qsg1_864[k]
                    + f_3 * pc_x[k] * qsh_1206[k];

        t_1606[k] = f_12 * osh_1207[k]
                    + f_4 * qsg0_865[k]
                    - f_5 * qsg1_865[k]
                    + f_3 * pc_x[k] * qsh_1207[k];

        t_1607[k] = f_12 * osh_972[k]
                    + f_3 * pc_z[k] * qsh_1203[k];
    }

#pragma omp simd aligned(t_1608, t_1609, t_1610, pc_x, pc_y, osh_996, osh_1209, osh_1211, \
                         qsg0_867, qsg0_869, qsg1_867, qsg1_869, qsh_1206, qsh_1209, \
                         qsh_1211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1608[k] = f_12 * osh_1209[k]
                    + f_4 * qsg0_867[k]
                    - f_5 * qsg1_867[k]
                    + f_3 * pc_x[k] * qsh_1209[k];

        t_1609[k] = f_20 * osh_996[k]
                    + f_3 * pc_y[k] * qsh_1206[k];

        t_1610[k] = f_12 * osh_1211[k]
                    + f_4 * qsg0_869[k]
                    - f_5 * qsg1_869[k]
                    + f_3 * pc_x[k] * qsh_1211[k];
    }

#pragma omp simd aligned(t_1611, t_1612, t_1613, t_1614, t_1615, pc_x, osh_1212, osh_1213, \
                         osh_1214, osh_1215, osh_1216, qsh_1212, qsh_1213, qsh_1214, qsh_1215, \
                         qsh_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1611[k] = f_12 * osh_1212[k]
                    + f_3 * pc_x[k] * qsh_1212[k];

        t_1612[k] = f_12 * osh_1213[k]
                    + f_3 * pc_x[k] * qsh_1213[k];

        t_1613[k] = f_12 * osh_1214[k]
                    + f_3 * pc_x[k] * qsh_1214[k];

        t_1614[k] = f_12 * osh_1215[k]
                    + f_3 * pc_x[k] * qsh_1215[k];

        t_1615[k] = f_12 * osh_1216[k]
                    + f_3 * pc_x[k] * qsh_1216[k];
    }

#pragma omp simd aligned(t_1616, t_1617, t_1618, pc_x, pc_y, pc_z, osh_981, osh_1002, \
                         osh_1217, qsg0_865, qsg1_865, qsh_1212, \
                         qsh_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1616[k] = f_12 * osh_1217[k]
                    + f_3 * pc_x[k] * qsh_1217[k];

        t_1617[k] = f_20 * osh_1002[k]
                    + f_1 * qsg0_865[k]
                    - f_2 * qsg1_865[k]
                    + f_3 * pc_y[k] * qsh_1212[k];

        t_1618[k] = f_12 * osh_981[k]
                    + f_3 * pc_z[k] * qsh_1212[k];
    }

#pragma omp simd aligned(t_1619, t_1620, t_1621, pc_y, osh_1004, osh_1005, osh_1006, qsg0_867, \
                         qsg0_868, qsg0_869, qsg1_867, qsg1_868, qsg1_869, qsh_1214, qsh_1215, \
                         qsh_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1619[k] = f_20 * osh_1004[k]
                    + f_8 * qsg0_867[k]
                    - f_9 * qsg1_867[k]
                    + f_3 * pc_y[k] * qsh_1214[k];

        t_1620[k] = f_20 * osh_1005[k]
                    + f_6 * qsg0_868[k]
                    - f_7 * qsg1_868[k]
                    + f_3 * pc_y[k] * qsh_1215[k];

        t_1621[k] = f_20 * osh_1006[k]
                    + f_4 * qsg0_869[k]
                    - f_5 * qsg1_869[k]
                    + f_3 * pc_y[k] * qsh_1216[k];
    }

#pragma omp simd aligned(t_1622, t_1623, t_1624, pc_x, pc_y, pc_z, osh_986, osh_1007, \
                         osh_1218, qsg0_869, qsg0_870, qsg1_869, qsg1_870, qsh_1217, \
                         qsh_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1622[k] = f_20 * osh_1007[k]
                    + f_3 * pc_y[k] * qsh_1217[k];

        t_1623[k] = f_12 * osh_986[k]
                    + f_1 * qsg0_869[k]
                    - f_2 * qsg1_869[k]
                    + f_3 * pc_z[k] * qsh_1217[k];

        t_1624[k] = f_12 * osh_1218[k]
                    + f_1 * qsg0_870[k]
                    - f_2 * qsg1_870[k]
                    + f_3 * pc_x[k] * qsh_1218[k];
    }

#pragma omp simd aligned(t_1625, t_1626, t_1627, t_1628, pc_x, pc_y, pc_z, osh_987, osh_1008, \
                         osh_1010, osh_1221, qsg0_873, qsg1_873, qsh_1218, qsh_1220, \
                         qsh_1221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1625[k] = f_21 * osh_1008[k]
                    + f_3 * pc_y[k] * qsh_1218[k];

        t_1626[k] = f_13 * osh_987[k]
                    + f_3 * pc_z[k] * qsh_1218[k];

        t_1627[k] = f_12 * osh_1221[k]
                    + f_8 * qsg0_873[k]
                    - f_9 * qsg1_873[k]
                    + f_3 * pc_x[k] * qsh_1221[k];

        t_1628[k] = f_21 * osh_1010[k]
                    + f_3 * pc_y[k] * qsh_1220[k];
    }

#pragma omp simd aligned(t_1629, t_1630, t_1631, pc_x, pc_z, osh_990, osh_1223, osh_1224, \
                         qsg0_875, qsg0_876, qsg1_875, qsg1_876, qsh_1221, qsh_1223, \
                         qsh_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1629[k] = f_12 * osh_1223[k]
                    + f_8 * qsg0_875[k]
                    - f_9 * qsg1_875[k]
                    + f_3 * pc_x[k] * qsh_1223[k];

        t_1630[k] = f_12 * osh_1224[k]
                    + f_6 * qsg0_876[k]
                    - f_7 * qsg1_876[k]
                    + f_3 * pc_x[k] * qsh_1224[k];

        t_1631[k] = f_13 * osh_990[k]
                    + f_3 * pc_z[k] * qsh_1221[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, pc_x, pc_y, osh_1013, osh_1227, osh_1228, \
                         qsg0_879, qsg0_880, qsg1_879, qsg1_880, qsh_1223, qsh_1227, \
                         qsh_1228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_21 * osh_1013[k]
                    + f_3 * pc_y[k] * qsh_1223[k];

        t_1633[k] = f_12 * osh_1227[k]
                    + f_6 * qsg0_879[k]
                    - f_7 * qsg1_879[k]
                    + f_3 * pc_x[k] * qsh_1227[k];

        t_1634[k] = f_12 * osh_1228[k]
                    + f_4 * qsg0_880[k]
                    - f_5 * qsg1_880[k]
                    + f_3 * pc_x[k] * qsh_1228[k];
    }

#pragma omp simd aligned(t_1635, t_1636, t_1637, pc_x, pc_y, pc_z, osh_993, osh_1017, \
                         osh_1230, qsg0_882, qsg1_882, qsh_1224, qsh_1227, \
                         qsh_1230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1635[k] = f_13 * osh_993[k]
                    + f_3 * pc_z[k] * qsh_1224[k];

        t_1636[k] = f_12 * osh_1230[k]
                    + f_4 * qsg0_882[k]
                    - f_5 * qsg1_882[k]
                    + f_3 * pc_x[k] * qsh_1230[k];

        t_1637[k] = f_21 * osh_1017[k]
                    + f_3 * pc_y[k] * qsh_1227[k];
    }

#pragma omp simd aligned(t_1638, t_1639, t_1640, t_1641, pc_x, osh_1232, osh_1233, osh_1234, \
                         osh_1235, qsg0_884, qsg1_884, qsh_1232, qsh_1233, qsh_1234, \
                         qsh_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1638[k] = f_12 * osh_1232[k]
                    + f_4 * qsg0_884[k]
                    - f_5 * qsg1_884[k]
                    + f_3 * pc_x[k] * qsh_1232[k];

        t_1639[k] = f_12 * osh_1233[k]
                    + f_3 * pc_x[k] * qsh_1233[k];

        t_1640[k] = f_12 * osh_1234[k]
                    + f_3 * pc_x[k] * qsh_1234[k];

        t_1641[k] = f_12 * osh_1235[k]
                    + f_3 * pc_x[k] * qsh_1235[k];
    }

#pragma omp simd aligned(t_1642, t_1643, t_1644, t_1645, pc_x, pc_y, osh_1023, osh_1236, \
                         osh_1237, osh_1238, qsg0_880, qsg1_880, qsh_1233, qsh_1236, qsh_1237, \
                         qsh_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1642[k] = f_12 * osh_1236[k]
                    + f_3 * pc_x[k] * qsh_1236[k];

        t_1643[k] = f_12 * osh_1237[k]
                    + f_3 * pc_x[k] * qsh_1237[k];

        t_1644[k] = f_12 * osh_1238[k]
                    + f_3 * pc_x[k] * qsh_1238[k];

        t_1645[k] = f_21 * osh_1023[k]
                    + f_1 * qsg0_880[k]
                    - f_2 * qsg1_880[k]
                    + f_3 * pc_y[k] * qsh_1233[k];
    }

#pragma omp simd aligned(t_1646, t_1647, t_1648, pc_y, pc_z, osh_1002, osh_1025, osh_1026, \
                         qsg0_882, qsg0_883, qsg1_882, qsg1_883, qsh_1233, qsh_1235, \
                         qsh_1236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1646[k] = f_13 * osh_1002[k]
                    + f_3 * pc_z[k] * qsh_1233[k];

        t_1647[k] = f_21 * osh_1025[k]
                    + f_8 * qsg0_882[k]
                    - f_9 * qsg1_882[k]
                    + f_3 * pc_y[k] * qsh_1235[k];

        t_1648[k] = f_21 * osh_1026[k]
                    + f_6 * qsg0_883[k]
                    - f_7 * qsg1_883[k]
                    + f_3 * pc_y[k] * qsh_1236[k];
    }

#pragma omp simd aligned(t_1649, t_1650, t_1651, pc_y, pc_z, osh_1007, osh_1027, osh_1028, \
                         qsg0_884, qsg1_884, qsh_1237, qsh_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1649[k] = f_21 * osh_1027[k]
                    + f_4 * qsg0_884[k]
                    - f_5 * qsg1_884[k]
                    + f_3 * pc_y[k] * qsh_1237[k];

        t_1650[k] = f_21 * osh_1028[k]
                    + f_3 * pc_y[k] * qsh_1238[k];

        t_1651[k] = f_13 * osh_1007[k]
                    + f_1 * qsg0_884[k]
                    - f_2 * qsg1_884[k]
                    + f_3 * pc_z[k] * qsh_1238[k];
    }

#pragma omp simd aligned(t_1652, t_1653, t_1654, pc_x, pc_y, pc_z, osh_1008, osh_1029, \
                         osh_1239, qsg0_885, qsg1_885, qsh_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1652[k] = f_12 * osh_1239[k]
                    + f_1 * qsg0_885[k]
                    - f_2 * qsg1_885[k]
                    + f_3 * pc_x[k] * qsh_1239[k];

        t_1653[k] = f_23 * osh_1029[k]
                    + f_3 * pc_y[k] * qsh_1239[k];

        t_1654[k] = f_14 * osh_1008[k]
                    + f_3 * pc_z[k] * qsh_1239[k];
    }

#pragma omp simd aligned(t_1655, t_1656, t_1657, pc_x, pc_y, osh_1031, osh_1242, osh_1244, \
                         qsg0_888, qsg0_890, qsg1_888, qsg1_890, qsh_1241, qsh_1242, \
                         qsh_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1655[k] = f_12 * osh_1242[k]
                    + f_8 * qsg0_888[k]
                    - f_9 * qsg1_888[k]
                    + f_3 * pc_x[k] * qsh_1242[k];

        t_1656[k] = f_23 * osh_1031[k]
                    + f_3 * pc_y[k] * qsh_1241[k];

        t_1657[k] = f_12 * osh_1244[k]
                    + f_8 * qsg0_890[k]
                    - f_9 * qsg1_890[k]
                    + f_3 * pc_x[k] * qsh_1244[k];
    }

#pragma omp simd aligned(t_1658, t_1659, t_1660, pc_x, pc_y, pc_z, osh_1011, osh_1034, \
                         osh_1245, qsg0_891, qsg1_891, qsh_1242, qsh_1244, \
                         qsh_1245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1658[k] = f_12 * osh_1245[k]
                    + f_6 * qsg0_891[k]
                    - f_7 * qsg1_891[k]
                    + f_3 * pc_x[k] * qsh_1245[k];

        t_1659[k] = f_14 * osh_1011[k]
                    + f_3 * pc_z[k] * qsh_1242[k];

        t_1660[k] = f_23 * osh_1034[k]
                    + f_3 * pc_y[k] * qsh_1244[k];
    }

#pragma omp simd aligned(t_1661, t_1662, t_1663, pc_x, pc_z, osh_1014, osh_1248, osh_1249, \
                         qsg0_894, qsg0_895, qsg1_894, qsg1_895, qsh_1245, qsh_1248, \
                         qsh_1249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1661[k] = f_12 * osh_1248[k]
                    + f_6 * qsg0_894[k]
                    - f_7 * qsg1_894[k]
                    + f_3 * pc_x[k] * qsh_1248[k];

        t_1662[k] = f_12 * osh_1249[k]
                    + f_4 * qsg0_895[k]
                    - f_5 * qsg1_895[k]
                    + f_3 * pc_x[k] * qsh_1249[k];

        t_1663[k] = f_14 * osh_1014[k]
                    + f_3 * pc_z[k] * qsh_1245[k];
    }

#pragma omp simd aligned(t_1664, t_1665, t_1666, pc_x, pc_y, osh_1038, osh_1251, osh_1253, \
                         qsg0_897, qsg0_899, qsg1_897, qsg1_899, qsh_1248, qsh_1251, \
                         qsh_1253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1664[k] = f_12 * osh_1251[k]
                    + f_4 * qsg0_897[k]
                    - f_5 * qsg1_897[k]
                    + f_3 * pc_x[k] * qsh_1251[k];

        t_1665[k] = f_23 * osh_1038[k]
                    + f_3 * pc_y[k] * qsh_1248[k];

        t_1666[k] = f_12 * osh_1253[k]
                    + f_4 * qsg0_899[k]
                    - f_5 * qsg1_899[k]
                    + f_3 * pc_x[k] * qsh_1253[k];
    }

#pragma omp simd aligned(t_1667, t_1668, t_1669, t_1670, t_1671, pc_x, osh_1254, osh_1255, \
                         osh_1256, osh_1257, osh_1258, qsh_1254, qsh_1255, qsh_1256, qsh_1257, \
                         qsh_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1667[k] = f_12 * osh_1254[k]
                    + f_3 * pc_x[k] * qsh_1254[k];

        t_1668[k] = f_12 * osh_1255[k]
                    + f_3 * pc_x[k] * qsh_1255[k];

        t_1669[k] = f_12 * osh_1256[k]
                    + f_3 * pc_x[k] * qsh_1256[k];

        t_1670[k] = f_12 * osh_1257[k]
                    + f_3 * pc_x[k] * qsh_1257[k];

        t_1671[k] = f_12 * osh_1258[k]
                    + f_3 * pc_x[k] * qsh_1258[k];
    }

#pragma omp simd aligned(t_1672, t_1673, t_1674, pc_x, pc_y, pc_z, osh_1023, osh_1044, \
                         osh_1259, qsg0_895, qsg1_895, qsh_1254, \
                         qsh_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1672[k] = f_12 * osh_1259[k]
                    + f_3 * pc_x[k] * qsh_1259[k];

        t_1673[k] = f_23 * osh_1044[k]
                    + f_1 * qsg0_895[k]
                    - f_2 * qsg1_895[k]
                    + f_3 * pc_y[k] * qsh_1254[k];

        t_1674[k] = f_14 * osh_1023[k]
                    + f_3 * pc_z[k] * qsh_1254[k];
    }

#pragma omp simd aligned(t_1675, t_1676, t_1677, pc_y, osh_1046, osh_1047, osh_1048, qsg0_897, \
                         qsg0_898, qsg0_899, qsg1_897, qsg1_898, qsg1_899, qsh_1256, qsh_1257, \
                         qsh_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1675[k] = f_23 * osh_1046[k]
                    + f_8 * qsg0_897[k]
                    - f_9 * qsg1_897[k]
                    + f_3 * pc_y[k] * qsh_1256[k];

        t_1676[k] = f_23 * osh_1047[k]
                    + f_6 * qsg0_898[k]
                    - f_7 * qsg1_898[k]
                    + f_3 * pc_y[k] * qsh_1257[k];

        t_1677[k] = f_23 * osh_1048[k]
                    + f_4 * qsg0_899[k]
                    - f_5 * qsg1_899[k]
                    + f_3 * pc_y[k] * qsh_1258[k];
    }

#pragma omp simd aligned(t_1678, t_1679, t_1680, pc_x, pc_y, pc_z, osh_1028, osh_1049, \
                         osh_1260, qsg0_899, qsg0_900, qsg1_899, qsg1_900, qsh_1259, \
                         qsh_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1678[k] = f_23 * osh_1049[k]
                    + f_3 * pc_y[k] * qsh_1259[k];

        t_1679[k] = f_14 * osh_1028[k]
                    + f_1 * qsg0_899[k]
                    - f_2 * qsg1_899[k]
                    + f_3 * pc_z[k] * qsh_1259[k];

        t_1680[k] = f_12 * osh_1260[k]
                    + f_1 * qsg0_900[k]
                    - f_2 * qsg1_900[k]
                    + f_3 * pc_x[k] * qsh_1260[k];
    }

#pragma omp simd aligned(t_1681, t_1682, t_1683, t_1684, pc_x, pc_y, pc_z, osh_1029, osh_1050, \
                         osh_1052, osh_1263, qsg0_903, qsg1_903, qsh_1260, qsh_1262, \
                         qsh_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1681[k] = f_22 * osh_1050[k]
                    + f_3 * pc_y[k] * qsh_1260[k];

        t_1682[k] = f_22 * osh_1029[k]
                    + f_3 * pc_z[k] * qsh_1260[k];

        t_1683[k] = f_12 * osh_1263[k]
                    + f_8 * qsg0_903[k]
                    - f_9 * qsg1_903[k]
                    + f_3 * pc_x[k] * qsh_1263[k];

        t_1684[k] = f_22 * osh_1052[k]
                    + f_3 * pc_y[k] * qsh_1262[k];
    }

#pragma omp simd aligned(t_1685, t_1686, t_1687, pc_x, pc_z, osh_1032, osh_1265, osh_1266, \
                         qsg0_905, qsg0_906, qsg1_905, qsg1_906, qsh_1263, qsh_1265, \
                         qsh_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1685[k] = f_12 * osh_1265[k]
                    + f_8 * qsg0_905[k]
                    - f_9 * qsg1_905[k]
                    + f_3 * pc_x[k] * qsh_1265[k];

        t_1686[k] = f_12 * osh_1266[k]
                    + f_6 * qsg0_906[k]
                    - f_7 * qsg1_906[k]
                    + f_3 * pc_x[k] * qsh_1266[k];

        t_1687[k] = f_22 * osh_1032[k]
                    + f_3 * pc_z[k] * qsh_1263[k];
    }

#pragma omp simd aligned(t_1688, t_1689, t_1690, pc_x, pc_y, osh_1055, osh_1269, osh_1270, \
                         qsg0_909, qsg0_910, qsg1_909, qsg1_910, qsh_1265, qsh_1269, \
                         qsh_1270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1688[k] = f_22 * osh_1055[k]
                    + f_3 * pc_y[k] * qsh_1265[k];

        t_1689[k] = f_12 * osh_1269[k]
                    + f_6 * qsg0_909[k]
                    - f_7 * qsg1_909[k]
                    + f_3 * pc_x[k] * qsh_1269[k];

        t_1690[k] = f_12 * osh_1270[k]
                    + f_4 * qsg0_910[k]
                    - f_5 * qsg1_910[k]
                    + f_3 * pc_x[k] * qsh_1270[k];
    }

#pragma omp simd aligned(t_1691, t_1692, t_1693, pc_x, pc_y, pc_z, osh_1035, osh_1059, \
                         osh_1272, qsg0_912, qsg1_912, qsh_1266, qsh_1269, \
                         qsh_1272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1691[k] = f_22 * osh_1035[k]
                    + f_3 * pc_z[k] * qsh_1266[k];

        t_1692[k] = f_12 * osh_1272[k]
                    + f_4 * qsg0_912[k]
                    - f_5 * qsg1_912[k]
                    + f_3 * pc_x[k] * qsh_1272[k];

        t_1693[k] = f_22 * osh_1059[k]
                    + f_3 * pc_y[k] * qsh_1269[k];
    }

#pragma omp simd aligned(t_1694, t_1695, t_1696, t_1697, pc_x, osh_1274, osh_1275, osh_1276, \
                         osh_1277, qsg0_914, qsg1_914, qsh_1274, qsh_1275, qsh_1276, \
                         qsh_1277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1694[k] = f_12 * osh_1274[k]
                    + f_4 * qsg0_914[k]
                    - f_5 * qsg1_914[k]
                    + f_3 * pc_x[k] * qsh_1274[k];

        t_1695[k] = f_12 * osh_1275[k]
                    + f_3 * pc_x[k] * qsh_1275[k];

        t_1696[k] = f_12 * osh_1276[k]
                    + f_3 * pc_x[k] * qsh_1276[k];

        t_1697[k] = f_12 * osh_1277[k]
                    + f_3 * pc_x[k] * qsh_1277[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, t_1701, pc_x, pc_y, osh_1065, osh_1278, \
                         osh_1279, osh_1280, qsg0_910, qsg1_910, qsh_1275, qsh_1278, qsh_1279, \
                         qsh_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_12 * osh_1278[k]
                    + f_3 * pc_x[k] * qsh_1278[k];

        t_1699[k] = f_12 * osh_1279[k]
                    + f_3 * pc_x[k] * qsh_1279[k];

        t_1700[k] = f_12 * osh_1280[k]
                    + f_3 * pc_x[k] * qsh_1280[k];

        t_1701[k] = f_22 * osh_1065[k]
                    + f_1 * qsg0_910[k]
                    - f_2 * qsg1_910[k]
                    + f_3 * pc_y[k] * qsh_1275[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsg0, const size_t qsg1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_19 = 4.5 / q;
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_1702 = buffer.data(target + 1702);
    auto *t_1703 = buffer.data(target + 1703);
    auto *t_1704 = buffer.data(target + 1704);
    auto *t_1705 = buffer.data(target + 1705);
    auto *t_1706 = buffer.data(target + 1706);
    auto *t_1707 = buffer.data(target + 1707);
    auto *t_1708 = buffer.data(target + 1708);
    auto *t_1709 = buffer.data(target + 1709);
    auto *t_1710 = buffer.data(target + 1710);
    auto *t_1711 = buffer.data(target + 1711);
    auto *t_1712 = buffer.data(target + 1712);
    auto *t_1713 = buffer.data(target + 1713);
    auto *t_1714 = buffer.data(target + 1714);
    auto *t_1715 = buffer.data(target + 1715);
    auto *t_1716 = buffer.data(target + 1716);
    auto *t_1717 = buffer.data(target + 1717);
    auto *t_1718 = buffer.data(target + 1718);
    auto *t_1719 = buffer.data(target + 1719);
    auto *t_1720 = buffer.data(target + 1720);
    auto *t_1721 = buffer.data(target + 1721);
    auto *t_1722 = buffer.data(target + 1722);
    auto *t_1723 = buffer.data(target + 1723);
    auto *t_1724 = buffer.data(target + 1724);
    auto *t_1725 = buffer.data(target + 1725);
    auto *t_1726 = buffer.data(target + 1726);
    auto *t_1727 = buffer.data(target + 1727);
    auto *t_1728 = buffer.data(target + 1728);
    auto *t_1729 = buffer.data(target + 1729);
    auto *t_1730 = buffer.data(target + 1730);
    auto *t_1731 = buffer.data(target + 1731);
    auto *t_1732 = buffer.data(target + 1732);
    auto *t_1733 = buffer.data(target + 1733);
    auto *t_1734 = buffer.data(target + 1734);
    auto *t_1735 = buffer.data(target + 1735);
    auto *t_1736 = buffer.data(target + 1736);
    auto *t_1737 = buffer.data(target + 1737);
    auto *t_1738 = buffer.data(target + 1738);
    auto *t_1739 = buffer.data(target + 1739);
    auto *t_1740 = buffer.data(target + 1740);
    auto *t_1741 = buffer.data(target + 1741);
    auto *t_1742 = buffer.data(target + 1742);
    auto *t_1743 = buffer.data(target + 1743);
    auto *t_1744 = buffer.data(target + 1744);
    auto *t_1745 = buffer.data(target + 1745);
    auto *t_1746 = buffer.data(target + 1746);
    auto *t_1747 = buffer.data(target + 1747);
    auto *t_1748 = buffer.data(target + 1748);
    auto *t_1749 = buffer.data(target + 1749);
    auto *t_1750 = buffer.data(target + 1750);
    auto *t_1751 = buffer.data(target + 1751);
    auto *t_1752 = buffer.data(target + 1752);
    auto *t_1753 = buffer.data(target + 1753);
    auto *t_1754 = buffer.data(target + 1754);
    auto *t_1755 = buffer.data(target + 1755);
    auto *t_1756 = buffer.data(target + 1756);
    auto *t_1757 = buffer.data(target + 1757);
    auto *t_1758 = buffer.data(target + 1758);
    auto *t_1759 = buffer.data(target + 1759);
    auto *t_1760 = buffer.data(target + 1760);
    auto *t_1761 = buffer.data(target + 1761);
    auto *t_1762 = buffer.data(target + 1762);
    auto *t_1763 = buffer.data(target + 1763);
    auto *t_1764 = buffer.data(target + 1764);
    auto *t_1765 = buffer.data(target + 1765);
    auto *t_1766 = buffer.data(target + 1766);
    auto *t_1767 = buffer.data(target + 1767);
    auto *t_1768 = buffer.data(target + 1768);
    auto *t_1769 = buffer.data(target + 1769);
    auto *t_1770 = buffer.data(target + 1770);
    auto *t_1771 = buffer.data(target + 1771);
    auto *t_1772 = buffer.data(target + 1772);
    auto *t_1773 = buffer.data(target + 1773);
    auto *t_1774 = buffer.data(target + 1774);
    auto *t_1775 = buffer.data(target + 1775);
    auto *t_1776 = buffer.data(target + 1776);
    auto *t_1777 = buffer.data(target + 1777);
    auto *t_1778 = buffer.data(target + 1778);
    auto *t_1779 = buffer.data(target + 1779);
    auto *t_1780 = buffer.data(target + 1780);
    auto *t_1781 = buffer.data(target + 1781);
    auto *t_1782 = buffer.data(target + 1782);
    auto *t_1783 = buffer.data(target + 1783);
    auto *t_1784 = buffer.data(target + 1784);
    auto *t_1785 = buffer.data(target + 1785);
    auto *t_1786 = buffer.data(target + 1786);
    auto *t_1787 = buffer.data(target + 1787);
    auto *t_1788 = buffer.data(target + 1788);
    auto *t_1789 = buffer.data(target + 1789);
    auto *t_1790 = buffer.data(target + 1790);
    auto *t_1791 = buffer.data(target + 1791);
    auto *t_1792 = buffer.data(target + 1792);
    auto *t_1793 = buffer.data(target + 1793);
    auto *t_1794 = buffer.data(target + 1794);
    auto *t_1795 = buffer.data(target + 1795);
    auto *t_1796 = buffer.data(target + 1796);
    auto *t_1797 = buffer.data(target + 1797);
    auto *t_1798 = buffer.data(target + 1798);
    auto *t_1799 = buffer.data(target + 1799);
    auto *t_1800 = buffer.data(target + 1800);
    auto *t_1801 = buffer.data(target + 1801);
    auto *t_1802 = buffer.data(target + 1802);
    auto *t_1803 = buffer.data(target + 1803);
    auto *t_1804 = buffer.data(target + 1804);
    auto *t_1805 = buffer.data(target + 1805);
    auto *t_1806 = buffer.data(target + 1806);
    auto *t_1807 = buffer.data(target + 1807);
    auto *t_1808 = buffer.data(target + 1808);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_1512 = buffer.data(osi0 + 1512);
    const auto *osi0_1515 = buffer.data(osi0 + 1515);
    const auto *osi0_1517 = buffer.data(osi0 + 1517);
    const auto *osi0_1518 = buffer.data(osi0 + 1518);
    const auto *osi0_1521 = buffer.data(osi0 + 1521);
    const auto *osi0_1522 = buffer.data(osi0 + 1522);
    const auto *osi0_1524 = buffer.data(osi0 + 1524);
    const auto *osi0_1526 = buffer.data(osi0 + 1526);

    const auto *osh_1044 = buffer.data(osh + 1044);
    const auto *osh_1049 = buffer.data(osh + 1049);
    const auto *osh_1050 = buffer.data(osh + 1050);
    const auto *osh_1053 = buffer.data(osh + 1053);
    const auto *osh_1056 = buffer.data(osh + 1056);
    const auto *osh_1065 = buffer.data(osh + 1065);
    const auto *osh_1067 = buffer.data(osh + 1067);
    const auto *osh_1068 = buffer.data(osh + 1068);
    const auto *osh_1069 = buffer.data(osh + 1069);
    const auto *osh_1070 = buffer.data(osh + 1070);
    const auto *osh_1071 = buffer.data(osh + 1071);
    const auto *osh_1073 = buffer.data(osh + 1073);
    const auto *osh_1074 = buffer.data(osh + 1074);
    const auto *osh_1076 = buffer.data(osh + 1076);
    const auto *osh_1077 = buffer.data(osh + 1077);
    const auto *osh_1080 = buffer.data(osh + 1080);
    const auto *osh_1086 = buffer.data(osh + 1086);
    const auto *osh_1088 = buffer.data(osh + 1088);
    const auto *osh_1089 = buffer.data(osh + 1089);
    const auto *osh_1090 = buffer.data(osh + 1090);
    const auto *osh_1091 = buffer.data(osh + 1091);
    const auto *osh_1092 = buffer.data(osh + 1092);
    const auto *osh_1094 = buffer.data(osh + 1094);
    const auto *osh_1095 = buffer.data(osh + 1095);
    const auto *osh_1097 = buffer.data(osh + 1097);
    const auto *osh_1098 = buffer.data(osh + 1098);
    const auto *osh_1101 = buffer.data(osh + 1101);
    const auto *osh_1107 = buffer.data(osh + 1107);
    const auto *osh_1109 = buffer.data(osh + 1109);
    const auto *osh_1110 = buffer.data(osh + 1110);
    const auto *osh_1111 = buffer.data(osh + 1111);
    const auto *osh_1112 = buffer.data(osh + 1112);
    const auto *osh_1113 = buffer.data(osh + 1113);
    const auto *osh_1115 = buffer.data(osh + 1115);
    const auto *osh_1116 = buffer.data(osh + 1116);
    const auto *osh_1118 = buffer.data(osh + 1118);
    const auto *osh_1119 = buffer.data(osh + 1119);
    const auto *osh_1122 = buffer.data(osh + 1122);
    const auto *osh_1128 = buffer.data(osh + 1128);
    const auto *osh_1130 = buffer.data(osh + 1130);
    const auto *osh_1131 = buffer.data(osh + 1131);
    const auto *osh_1132 = buffer.data(osh + 1132);
    const auto *osh_1133 = buffer.data(osh + 1133);
    const auto *osh_1134 = buffer.data(osh + 1134);
    const auto *osh_1135 = buffer.data(osh + 1135);
    const auto *osh_1136 = buffer.data(osh + 1136);
    const auto *osh_1137 = buffer.data(osh + 1137);
    const auto *osh_1139 = buffer.data(osh + 1139);
    const auto *osh_1140 = buffer.data(osh + 1140);
    const auto *osh_1142 = buffer.data(osh + 1142);
    const auto *osh_1143 = buffer.data(osh + 1143);
    const auto *osh_1281 = buffer.data(osh + 1281);
    const auto *osh_1284 = buffer.data(osh + 1284);
    const auto *osh_1286 = buffer.data(osh + 1286);
    const auto *osh_1287 = buffer.data(osh + 1287);
    const auto *osh_1290 = buffer.data(osh + 1290);
    const auto *osh_1291 = buffer.data(osh + 1291);
    const auto *osh_1293 = buffer.data(osh + 1293);
    const auto *osh_1295 = buffer.data(osh + 1295);
    const auto *osh_1296 = buffer.data(osh + 1296);
    const auto *osh_1297 = buffer.data(osh + 1297);
    const auto *osh_1298 = buffer.data(osh + 1298);
    const auto *osh_1299 = buffer.data(osh + 1299);
    const auto *osh_1300 = buffer.data(osh + 1300);
    const auto *osh_1301 = buffer.data(osh + 1301);
    const auto *osh_1302 = buffer.data(osh + 1302);
    const auto *osh_1305 = buffer.data(osh + 1305);
    const auto *osh_1307 = buffer.data(osh + 1307);
    const auto *osh_1308 = buffer.data(osh + 1308);
    const auto *osh_1311 = buffer.data(osh + 1311);
    const auto *osh_1312 = buffer.data(osh + 1312);
    const auto *osh_1314 = buffer.data(osh + 1314);
    const auto *osh_1316 = buffer.data(osh + 1316);
    const auto *osh_1317 = buffer.data(osh + 1317);
    const auto *osh_1318 = buffer.data(osh + 1318);
    const auto *osh_1319 = buffer.data(osh + 1319);
    const auto *osh_1320 = buffer.data(osh + 1320);
    const auto *osh_1321 = buffer.data(osh + 1321);
    const auto *osh_1322 = buffer.data(osh + 1322);
    const auto *osh_1323 = buffer.data(osh + 1323);
    const auto *osh_1326 = buffer.data(osh + 1326);
    const auto *osh_1328 = buffer.data(osh + 1328);
    const auto *osh_1329 = buffer.data(osh + 1329);
    const auto *osh_1332 = buffer.data(osh + 1332);
    const auto *osh_1333 = buffer.data(osh + 1333);
    const auto *osh_1335 = buffer.data(osh + 1335);
    const auto *osh_1337 = buffer.data(osh + 1337);
    const auto *osh_1338 = buffer.data(osh + 1338);
    const auto *osh_1339 = buffer.data(osh + 1339);
    const auto *osh_1340 = buffer.data(osh + 1340);
    const auto *osh_1341 = buffer.data(osh + 1341);
    const auto *osh_1342 = buffer.data(osh + 1342);
    const auto *osh_1343 = buffer.data(osh + 1343);
    const auto *osh_1359 = buffer.data(osh + 1359);
    const auto *osh_1360 = buffer.data(osh + 1360);

    const auto *osi1_1512 = buffer.data(osi1 + 1512);
    const auto *osi1_1515 = buffer.data(osi1 + 1515);
    const auto *osi1_1517 = buffer.data(osi1 + 1517);
    const auto *osi1_1518 = buffer.data(osi1 + 1518);
    const auto *osi1_1521 = buffer.data(osi1 + 1521);
    const auto *osi1_1522 = buffer.data(osi1 + 1522);
    const auto *osi1_1524 = buffer.data(osi1 + 1524);
    const auto *osi1_1526 = buffer.data(osi1 + 1526);

    const auto *qsg0_912 = buffer.data(qsg0 + 912);
    const auto *qsg0_913 = buffer.data(qsg0 + 913);
    const auto *qsg0_914 = buffer.data(qsg0 + 914);
    const auto *qsg0_915 = buffer.data(qsg0 + 915);
    const auto *qsg0_918 = buffer.data(qsg0 + 918);
    const auto *qsg0_920 = buffer.data(qsg0 + 920);
    const auto *qsg0_921 = buffer.data(qsg0 + 921);
    const auto *qsg0_924 = buffer.data(qsg0 + 924);
    const auto *qsg0_925 = buffer.data(qsg0 + 925);
    const auto *qsg0_927 = buffer.data(qsg0 + 927);
    const auto *qsg0_928 = buffer.data(qsg0 + 928);
    const auto *qsg0_929 = buffer.data(qsg0 + 929);
    const auto *qsg0_930 = buffer.data(qsg0 + 930);
    const auto *qsg0_933 = buffer.data(qsg0 + 933);
    const auto *qsg0_935 = buffer.data(qsg0 + 935);
    const auto *qsg0_936 = buffer.data(qsg0 + 936);
    const auto *qsg0_939 = buffer.data(qsg0 + 939);
    const auto *qsg0_940 = buffer.data(qsg0 + 940);
    const auto *qsg0_942 = buffer.data(qsg0 + 942);
    const auto *qsg0_943 = buffer.data(qsg0 + 943);
    const auto *qsg0_944 = buffer.data(qsg0 + 944);
    const auto *qsg0_945 = buffer.data(qsg0 + 945);
    const auto *qsg0_948 = buffer.data(qsg0 + 948);
    const auto *qsg0_950 = buffer.data(qsg0 + 950);
    const auto *qsg0_951 = buffer.data(qsg0 + 951);
    const auto *qsg0_954 = buffer.data(qsg0 + 954);
    const auto *qsg0_955 = buffer.data(qsg0 + 955);
    const auto *qsg0_957 = buffer.data(qsg0 + 957);
    const auto *qsg0_958 = buffer.data(qsg0 + 958);
    const auto *qsg0_959 = buffer.data(qsg0 + 959);

    const auto *qsg1_912 = buffer.data(qsg1 + 912);
    const auto *qsg1_913 = buffer.data(qsg1 + 913);
    const auto *qsg1_914 = buffer.data(qsg1 + 914);
    const auto *qsg1_915 = buffer.data(qsg1 + 915);
    const auto *qsg1_918 = buffer.data(qsg1 + 918);
    const auto *qsg1_920 = buffer.data(qsg1 + 920);
    const auto *qsg1_921 = buffer.data(qsg1 + 921);
    const auto *qsg1_924 = buffer.data(qsg1 + 924);
    const auto *qsg1_925 = buffer.data(qsg1 + 925);
    const auto *qsg1_927 = buffer.data(qsg1 + 927);
    const auto *qsg1_928 = buffer.data(qsg1 + 928);
    const auto *qsg1_929 = buffer.data(qsg1 + 929);
    const auto *qsg1_930 = buffer.data(qsg1 + 930);
    const auto *qsg1_933 = buffer.data(qsg1 + 933);
    const auto *qsg1_935 = buffer.data(qsg1 + 935);
    const auto *qsg1_936 = buffer.data(qsg1 + 936);
    const auto *qsg1_939 = buffer.data(qsg1 + 939);
    const auto *qsg1_940 = buffer.data(qsg1 + 940);
    const auto *qsg1_942 = buffer.data(qsg1 + 942);
    const auto *qsg1_943 = buffer.data(qsg1 + 943);
    const auto *qsg1_944 = buffer.data(qsg1 + 944);
    const auto *qsg1_945 = buffer.data(qsg1 + 945);
    const auto *qsg1_948 = buffer.data(qsg1 + 948);
    const auto *qsg1_950 = buffer.data(qsg1 + 950);
    const auto *qsg1_951 = buffer.data(qsg1 + 951);
    const auto *qsg1_954 = buffer.data(qsg1 + 954);
    const auto *qsg1_955 = buffer.data(qsg1 + 955);
    const auto *qsg1_957 = buffer.data(qsg1 + 957);
    const auto *qsg1_958 = buffer.data(qsg1 + 958);
    const auto *qsg1_959 = buffer.data(qsg1 + 959);

    const auto *qsh_1275 = buffer.data(qsh + 1275);
    const auto *qsh_1277 = buffer.data(qsh + 1277);
    const auto *qsh_1278 = buffer.data(qsh + 1278);
    const auto *qsh_1279 = buffer.data(qsh + 1279);
    const auto *qsh_1280 = buffer.data(qsh + 1280);
    const auto *qsh_1281 = buffer.data(qsh + 1281);
    const auto *qsh_1283 = buffer.data(qsh + 1283);
    const auto *qsh_1284 = buffer.data(qsh + 1284);
    const auto *qsh_1286 = buffer.data(qsh + 1286);
    const auto *qsh_1287 = buffer.data(qsh + 1287);
    const auto *qsh_1290 = buffer.data(qsh + 1290);
    const auto *qsh_1291 = buffer.data(qsh + 1291);
    const auto *qsh_1293 = buffer.data(qsh + 1293);
    const auto *qsh_1295 = buffer.data(qsh + 1295);
    const auto *qsh_1296 = buffer.data(qsh + 1296);
    const auto *qsh_1297 = buffer.data(qsh + 1297);
    const auto *qsh_1298 = buffer.data(qsh + 1298);
    const auto *qsh_1299 = buffer.data(qsh + 1299);
    const auto *qsh_1300 = buffer.data(qsh + 1300);
    const auto *qsh_1301 = buffer.data(qsh + 1301);
    const auto *qsh_1302 = buffer.data(qsh + 1302);
    const auto *qsh_1304 = buffer.data(qsh + 1304);
    const auto *qsh_1305 = buffer.data(qsh + 1305);
    const auto *qsh_1307 = buffer.data(qsh + 1307);
    const auto *qsh_1308 = buffer.data(qsh + 1308);
    const auto *qsh_1311 = buffer.data(qsh + 1311);
    const auto *qsh_1312 = buffer.data(qsh + 1312);
    const auto *qsh_1314 = buffer.data(qsh + 1314);
    const auto *qsh_1316 = buffer.data(qsh + 1316);
    const auto *qsh_1317 = buffer.data(qsh + 1317);
    const auto *qsh_1318 = buffer.data(qsh + 1318);
    const auto *qsh_1319 = buffer.data(qsh + 1319);
    const auto *qsh_1320 = buffer.data(qsh + 1320);
    const auto *qsh_1321 = buffer.data(qsh + 1321);
    const auto *qsh_1322 = buffer.data(qsh + 1322);
    const auto *qsh_1323 = buffer.data(qsh + 1323);
    const auto *qsh_1325 = buffer.data(qsh + 1325);
    const auto *qsh_1326 = buffer.data(qsh + 1326);
    const auto *qsh_1328 = buffer.data(qsh + 1328);
    const auto *qsh_1329 = buffer.data(qsh + 1329);
    const auto *qsh_1332 = buffer.data(qsh + 1332);
    const auto *qsh_1333 = buffer.data(qsh + 1333);
    const auto *qsh_1335 = buffer.data(qsh + 1335);
    const auto *qsh_1337 = buffer.data(qsh + 1337);
    const auto *qsh_1338 = buffer.data(qsh + 1338);
    const auto *qsh_1339 = buffer.data(qsh + 1339);
    const auto *qsh_1340 = buffer.data(qsh + 1340);
    const auto *qsh_1341 = buffer.data(qsh + 1341);
    const auto *qsh_1342 = buffer.data(qsh + 1342);
    const auto *qsh_1343 = buffer.data(qsh + 1343);
    const auto *qsh_1344 = buffer.data(qsh + 1344);
    const auto *qsh_1346 = buffer.data(qsh + 1346);
    const auto *qsh_1347 = buffer.data(qsh + 1347);
    const auto *qsh_1349 = buffer.data(qsh + 1349);
    const auto *qsh_1350 = buffer.data(qsh + 1350);
    const auto *qsh_1353 = buffer.data(qsh + 1353);
    const auto *qsh_1359 = buffer.data(qsh + 1359);
    const auto *qsh_1360 = buffer.data(qsh + 1360);

#pragma omp simd aligned(t_1702, t_1703, t_1704, pc_y, pc_z, osh_1044, osh_1067, osh_1068, \
                         qsg0_912, qsg0_913, qsg1_912, qsg1_913, qsh_1275, qsh_1277, \
                         qsh_1278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1702[k] = f_22 * osh_1044[k]
                    + f_3 * pc_z[k] * qsh_1275[k];

        t_1703[k] = f_22 * osh_1067[k]
                    + f_8 * qsg0_912[k]
                    - f_9 * qsg1_912[k]
                    + f_3 * pc_y[k] * qsh_1277[k];

        t_1704[k] = f_22 * osh_1068[k]
                    + f_6 * qsg0_913[k]
                    - f_7 * qsg1_913[k]
                    + f_3 * pc_y[k] * qsh_1278[k];
    }

#pragma omp simd aligned(t_1705, t_1706, t_1707, pc_y, pc_z, osh_1049, osh_1069, osh_1070, \
                         qsg0_914, qsg1_914, qsh_1279, qsh_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1705[k] = f_22 * osh_1069[k]
                    + f_4 * qsg0_914[k]
                    - f_5 * qsg1_914[k]
                    + f_3 * pc_y[k] * qsh_1279[k];

        t_1706[k] = f_22 * osh_1070[k]
                    + f_3 * pc_y[k] * qsh_1280[k];

        t_1707[k] = f_22 * osh_1049[k]
                    + f_1 * qsg0_914[k]
                    - f_2 * qsg1_914[k]
                    + f_3 * pc_z[k] * qsh_1280[k];
    }

#pragma omp simd aligned(t_1708, t_1709, t_1710, pc_x, pc_y, pc_z, osh_1050, osh_1071, \
                         osh_1281, qsg0_915, qsg1_915, qsh_1281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1708[k] = f_12 * osh_1281[k]
                    + f_1 * qsg0_915[k]
                    - f_2 * qsg1_915[k]
                    + f_3 * pc_x[k] * qsh_1281[k];

        t_1709[k] = f_14 * osh_1071[k]
                    + f_3 * pc_y[k] * qsh_1281[k];

        t_1710[k] = f_23 * osh_1050[k]
                    + f_3 * pc_z[k] * qsh_1281[k];
    }

#pragma omp simd aligned(t_1711, t_1712, t_1713, pc_x, pc_y, osh_1073, osh_1284, osh_1286, \
                         qsg0_918, qsg0_920, qsg1_918, qsg1_920, qsh_1283, qsh_1284, \
                         qsh_1286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1711[k] = f_12 * osh_1284[k]
                    + f_8 * qsg0_918[k]
                    - f_9 * qsg1_918[k]
                    + f_3 * pc_x[k] * qsh_1284[k];

        t_1712[k] = f_14 * osh_1073[k]
                    + f_3 * pc_y[k] * qsh_1283[k];

        t_1713[k] = f_12 * osh_1286[k]
                    + f_8 * qsg0_920[k]
                    - f_9 * qsg1_920[k]
                    + f_3 * pc_x[k] * qsh_1286[k];
    }

#pragma omp simd aligned(t_1714, t_1715, t_1716, pc_x, pc_y, pc_z, osh_1053, osh_1076, \
                         osh_1287, qsg0_921, qsg1_921, qsh_1284, qsh_1286, \
                         qsh_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1714[k] = f_12 * osh_1287[k]
                    + f_6 * qsg0_921[k]
                    - f_7 * qsg1_921[k]
                    + f_3 * pc_x[k] * qsh_1287[k];

        t_1715[k] = f_23 * osh_1053[k]
                    + f_3 * pc_z[k] * qsh_1284[k];

        t_1716[k] = f_14 * osh_1076[k]
                    + f_3 * pc_y[k] * qsh_1286[k];
    }

#pragma omp simd aligned(t_1717, t_1718, t_1719, pc_x, pc_z, osh_1056, osh_1290, osh_1291, \
                         qsg0_924, qsg0_925, qsg1_924, qsg1_925, qsh_1287, qsh_1290, \
                         qsh_1291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1717[k] = f_12 * osh_1290[k]
                    + f_6 * qsg0_924[k]
                    - f_7 * qsg1_924[k]
                    + f_3 * pc_x[k] * qsh_1290[k];

        t_1718[k] = f_12 * osh_1291[k]
                    + f_4 * qsg0_925[k]
                    - f_5 * qsg1_925[k]
                    + f_3 * pc_x[k] * qsh_1291[k];

        t_1719[k] = f_23 * osh_1056[k]
                    + f_3 * pc_z[k] * qsh_1287[k];
    }

#pragma omp simd aligned(t_1720, t_1721, t_1722, pc_x, pc_y, osh_1080, osh_1293, osh_1295, \
                         qsg0_927, qsg0_929, qsg1_927, qsg1_929, qsh_1290, qsh_1293, \
                         qsh_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1720[k] = f_12 * osh_1293[k]
                    + f_4 * qsg0_927[k]
                    - f_5 * qsg1_927[k]
                    + f_3 * pc_x[k] * qsh_1293[k];

        t_1721[k] = f_14 * osh_1080[k]
                    + f_3 * pc_y[k] * qsh_1290[k];

        t_1722[k] = f_12 * osh_1295[k]
                    + f_4 * qsg0_929[k]
                    - f_5 * qsg1_929[k]
                    + f_3 * pc_x[k] * qsh_1295[k];
    }

#pragma omp simd aligned(t_1723, t_1724, t_1725, t_1726, t_1727, pc_x, osh_1296, osh_1297, \
                         osh_1298, osh_1299, osh_1300, qsh_1296, qsh_1297, qsh_1298, qsh_1299, \
                         qsh_1300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1723[k] = f_12 * osh_1296[k]
                    + f_3 * pc_x[k] * qsh_1296[k];

        t_1724[k] = f_12 * osh_1297[k]
                    + f_3 * pc_x[k] * qsh_1297[k];

        t_1725[k] = f_12 * osh_1298[k]
                    + f_3 * pc_x[k] * qsh_1298[k];

        t_1726[k] = f_12 * osh_1299[k]
                    + f_3 * pc_x[k] * qsh_1299[k];

        t_1727[k] = f_12 * osh_1300[k]
                    + f_3 * pc_x[k] * qsh_1300[k];
    }

#pragma omp simd aligned(t_1728, t_1729, t_1730, pc_x, pc_y, pc_z, osh_1065, osh_1086, \
                         osh_1301, qsg0_925, qsg1_925, qsh_1296, \
                         qsh_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1728[k] = f_12 * osh_1301[k]
                    + f_3 * pc_x[k] * qsh_1301[k];

        t_1729[k] = f_14 * osh_1086[k]
                    + f_1 * qsg0_925[k]
                    - f_2 * qsg1_925[k]
                    + f_3 * pc_y[k] * qsh_1296[k];

        t_1730[k] = f_23 * osh_1065[k]
                    + f_3 * pc_z[k] * qsh_1296[k];
    }

#pragma omp simd aligned(t_1731, t_1732, t_1733, pc_y, osh_1088, osh_1089, osh_1090, qsg0_927, \
                         qsg0_928, qsg0_929, qsg1_927, qsg1_928, qsg1_929, qsh_1298, qsh_1299, \
                         qsh_1300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1731[k] = f_14 * osh_1088[k]
                    + f_8 * qsg0_927[k]
                    - f_9 * qsg1_927[k]
                    + f_3 * pc_y[k] * qsh_1298[k];

        t_1732[k] = f_14 * osh_1089[k]
                    + f_6 * qsg0_928[k]
                    - f_7 * qsg1_928[k]
                    + f_3 * pc_y[k] * qsh_1299[k];

        t_1733[k] = f_14 * osh_1090[k]
                    + f_4 * qsg0_929[k]
                    - f_5 * qsg1_929[k]
                    + f_3 * pc_y[k] * qsh_1300[k];
    }

#pragma omp simd aligned(t_1734, t_1735, t_1736, pc_x, pc_y, pc_z, osh_1070, osh_1091, \
                         osh_1302, qsg0_929, qsg0_930, qsg1_929, qsg1_930, qsh_1301, \
                         qsh_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1734[k] = f_14 * osh_1091[k]
                    + f_3 * pc_y[k] * qsh_1301[k];

        t_1735[k] = f_23 * osh_1070[k]
                    + f_1 * qsg0_929[k]
                    - f_2 * qsg1_929[k]
                    + f_3 * pc_z[k] * qsh_1301[k];

        t_1736[k] = f_12 * osh_1302[k]
                    + f_1 * qsg0_930[k]
                    - f_2 * qsg1_930[k]
                    + f_3 * pc_x[k] * qsh_1302[k];
    }

#pragma omp simd aligned(t_1737, t_1738, t_1739, t_1740, pc_x, pc_y, pc_z, osh_1071, osh_1092, \
                         osh_1094, osh_1305, qsg0_933, qsg1_933, qsh_1302, qsh_1304, \
                         qsh_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1737[k] = f_13 * osh_1092[k]
                    + f_3 * pc_y[k] * qsh_1302[k];

        t_1738[k] = f_21 * osh_1071[k]
                    + f_3 * pc_z[k] * qsh_1302[k];

        t_1739[k] = f_12 * osh_1305[k]
                    + f_8 * qsg0_933[k]
                    - f_9 * qsg1_933[k]
                    + f_3 * pc_x[k] * qsh_1305[k];

        t_1740[k] = f_13 * osh_1094[k]
                    + f_3 * pc_y[k] * qsh_1304[k];
    }

#pragma omp simd aligned(t_1741, t_1742, t_1743, pc_x, pc_z, osh_1074, osh_1307, osh_1308, \
                         qsg0_935, qsg0_936, qsg1_935, qsg1_936, qsh_1305, qsh_1307, \
                         qsh_1308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1741[k] = f_12 * osh_1307[k]
                    + f_8 * qsg0_935[k]
                    - f_9 * qsg1_935[k]
                    + f_3 * pc_x[k] * qsh_1307[k];

        t_1742[k] = f_12 * osh_1308[k]
                    + f_6 * qsg0_936[k]
                    - f_7 * qsg1_936[k]
                    + f_3 * pc_x[k] * qsh_1308[k];

        t_1743[k] = f_21 * osh_1074[k]
                    + f_3 * pc_z[k] * qsh_1305[k];
    }

#pragma omp simd aligned(t_1744, t_1745, t_1746, pc_x, pc_y, osh_1097, osh_1311, osh_1312, \
                         qsg0_939, qsg0_940, qsg1_939, qsg1_940, qsh_1307, qsh_1311, \
                         qsh_1312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1744[k] = f_13 * osh_1097[k]
                    + f_3 * pc_y[k] * qsh_1307[k];

        t_1745[k] = f_12 * osh_1311[k]
                    + f_6 * qsg0_939[k]
                    - f_7 * qsg1_939[k]
                    + f_3 * pc_x[k] * qsh_1311[k];

        t_1746[k] = f_12 * osh_1312[k]
                    + f_4 * qsg0_940[k]
                    - f_5 * qsg1_940[k]
                    + f_3 * pc_x[k] * qsh_1312[k];
    }

#pragma omp simd aligned(t_1747, t_1748, t_1749, pc_x, pc_y, pc_z, osh_1077, osh_1101, \
                         osh_1314, qsg0_942, qsg1_942, qsh_1308, qsh_1311, \
                         qsh_1314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1747[k] = f_21 * osh_1077[k]
                    + f_3 * pc_z[k] * qsh_1308[k];

        t_1748[k] = f_12 * osh_1314[k]
                    + f_4 * qsg0_942[k]
                    - f_5 * qsg1_942[k]
                    + f_3 * pc_x[k] * qsh_1314[k];

        t_1749[k] = f_13 * osh_1101[k]
                    + f_3 * pc_y[k] * qsh_1311[k];
    }

#pragma omp simd aligned(t_1750, t_1751, t_1752, t_1753, pc_x, osh_1316, osh_1317, osh_1318, \
                         osh_1319, qsg0_944, qsg1_944, qsh_1316, qsh_1317, qsh_1318, \
                         qsh_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1750[k] = f_12 * osh_1316[k]
                    + f_4 * qsg0_944[k]
                    - f_5 * qsg1_944[k]
                    + f_3 * pc_x[k] * qsh_1316[k];

        t_1751[k] = f_12 * osh_1317[k]
                    + f_3 * pc_x[k] * qsh_1317[k];

        t_1752[k] = f_12 * osh_1318[k]
                    + f_3 * pc_x[k] * qsh_1318[k];

        t_1753[k] = f_12 * osh_1319[k]
                    + f_3 * pc_x[k] * qsh_1319[k];
    }

#pragma omp simd aligned(t_1754, t_1755, t_1756, t_1757, pc_x, pc_y, osh_1107, osh_1320, \
                         osh_1321, osh_1322, qsg0_940, qsg1_940, qsh_1317, qsh_1320, qsh_1321, \
                         qsh_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1754[k] = f_12 * osh_1320[k]
                    + f_3 * pc_x[k] * qsh_1320[k];

        t_1755[k] = f_12 * osh_1321[k]
                    + f_3 * pc_x[k] * qsh_1321[k];

        t_1756[k] = f_12 * osh_1322[k]
                    + f_3 * pc_x[k] * qsh_1322[k];

        t_1757[k] = f_13 * osh_1107[k]
                    + f_1 * qsg0_940[k]
                    - f_2 * qsg1_940[k]
                    + f_3 * pc_y[k] * qsh_1317[k];
    }

#pragma omp simd aligned(t_1758, t_1759, t_1760, pc_y, pc_z, osh_1086, osh_1109, osh_1110, \
                         qsg0_942, qsg0_943, qsg1_942, qsg1_943, qsh_1317, qsh_1319, \
                         qsh_1320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1758[k] = f_21 * osh_1086[k]
                    + f_3 * pc_z[k] * qsh_1317[k];

        t_1759[k] = f_13 * osh_1109[k]
                    + f_8 * qsg0_942[k]
                    - f_9 * qsg1_942[k]
                    + f_3 * pc_y[k] * qsh_1319[k];

        t_1760[k] = f_13 * osh_1110[k]
                    + f_6 * qsg0_943[k]
                    - f_7 * qsg1_943[k]
                    + f_3 * pc_y[k] * qsh_1320[k];
    }

#pragma omp simd aligned(t_1761, t_1762, t_1763, pc_y, pc_z, osh_1091, osh_1111, osh_1112, \
                         qsg0_944, qsg1_944, qsh_1321, qsh_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1761[k] = f_13 * osh_1111[k]
                    + f_4 * qsg0_944[k]
                    - f_5 * qsg1_944[k]
                    + f_3 * pc_y[k] * qsh_1321[k];

        t_1762[k] = f_13 * osh_1112[k]
                    + f_3 * pc_y[k] * qsh_1322[k];

        t_1763[k] = f_21 * osh_1091[k]
                    + f_1 * qsg0_944[k]
                    - f_2 * qsg1_944[k]
                    + f_3 * pc_z[k] * qsh_1322[k];
    }

#pragma omp simd aligned(t_1764, t_1765, t_1766, pc_x, pc_y, pc_z, osh_1092, osh_1113, \
                         osh_1323, qsg0_945, qsg1_945, qsh_1323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1764[k] = f_12 * osh_1323[k]
                    + f_1 * qsg0_945[k]
                    - f_2 * qsg1_945[k]
                    + f_3 * pc_x[k] * qsh_1323[k];

        t_1765[k] = f_12 * osh_1113[k]
                    + f_3 * pc_y[k] * qsh_1323[k];

        t_1766[k] = f_20 * osh_1092[k]
                    + f_3 * pc_z[k] * qsh_1323[k];
    }

#pragma omp simd aligned(t_1767, t_1768, t_1769, pc_x, pc_y, osh_1115, osh_1326, osh_1328, \
                         qsg0_948, qsg0_950, qsg1_948, qsg1_950, qsh_1325, qsh_1326, \
                         qsh_1328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1767[k] = f_12 * osh_1326[k]
                    + f_8 * qsg0_948[k]
                    - f_9 * qsg1_948[k]
                    + f_3 * pc_x[k] * qsh_1326[k];

        t_1768[k] = f_12 * osh_1115[k]
                    + f_3 * pc_y[k] * qsh_1325[k];

        t_1769[k] = f_12 * osh_1328[k]
                    + f_8 * qsg0_950[k]
                    - f_9 * qsg1_950[k]
                    + f_3 * pc_x[k] * qsh_1328[k];
    }

#pragma omp simd aligned(t_1770, t_1771, t_1772, pc_x, pc_y, pc_z, osh_1095, osh_1118, \
                         osh_1329, qsg0_951, qsg1_951, qsh_1326, qsh_1328, \
                         qsh_1329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1770[k] = f_12 * osh_1329[k]
                    + f_6 * qsg0_951[k]
                    - f_7 * qsg1_951[k]
                    + f_3 * pc_x[k] * qsh_1329[k];

        t_1771[k] = f_20 * osh_1095[k]
                    + f_3 * pc_z[k] * qsh_1326[k];

        t_1772[k] = f_12 * osh_1118[k]
                    + f_3 * pc_y[k] * qsh_1328[k];
    }

#pragma omp simd aligned(t_1773, t_1774, t_1775, pc_x, pc_z, osh_1098, osh_1332, osh_1333, \
                         qsg0_954, qsg0_955, qsg1_954, qsg1_955, qsh_1329, qsh_1332, \
                         qsh_1333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1773[k] = f_12 * osh_1332[k]
                    + f_6 * qsg0_954[k]
                    - f_7 * qsg1_954[k]
                    + f_3 * pc_x[k] * qsh_1332[k];

        t_1774[k] = f_12 * osh_1333[k]
                    + f_4 * qsg0_955[k]
                    - f_5 * qsg1_955[k]
                    + f_3 * pc_x[k] * qsh_1333[k];

        t_1775[k] = f_20 * osh_1098[k]
                    + f_3 * pc_z[k] * qsh_1329[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, pc_x, pc_y, osh_1122, osh_1335, osh_1337, \
                         qsg0_957, qsg0_959, qsg1_957, qsg1_959, qsh_1332, qsh_1335, \
                         qsh_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = f_12 * osh_1335[k]
                    + f_4 * qsg0_957[k]
                    - f_5 * qsg1_957[k]
                    + f_3 * pc_x[k] * qsh_1335[k];

        t_1777[k] = f_12 * osh_1122[k]
                    + f_3 * pc_y[k] * qsh_1332[k];

        t_1778[k] = f_12 * osh_1337[k]
                    + f_4 * qsg0_959[k]
                    - f_5 * qsg1_959[k]
                    + f_3 * pc_x[k] * qsh_1337[k];
    }

#pragma omp simd aligned(t_1779, t_1780, t_1781, t_1782, t_1783, pc_x, osh_1338, osh_1339, \
                         osh_1340, osh_1341, osh_1342, qsh_1338, qsh_1339, qsh_1340, qsh_1341, \
                         qsh_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1779[k] = f_12 * osh_1338[k]
                    + f_3 * pc_x[k] * qsh_1338[k];

        t_1780[k] = f_12 * osh_1339[k]
                    + f_3 * pc_x[k] * qsh_1339[k];

        t_1781[k] = f_12 * osh_1340[k]
                    + f_3 * pc_x[k] * qsh_1340[k];

        t_1782[k] = f_12 * osh_1341[k]
                    + f_3 * pc_x[k] * qsh_1341[k];

        t_1783[k] = f_12 * osh_1342[k]
                    + f_3 * pc_x[k] * qsh_1342[k];
    }

#pragma omp simd aligned(t_1784, t_1785, t_1786, pc_x, pc_y, pc_z, osh_1107, osh_1128, \
                         osh_1343, qsg0_955, qsg1_955, qsh_1338, \
                         qsh_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1784[k] = f_12 * osh_1343[k]
                    + f_3 * pc_x[k] * qsh_1343[k];

        t_1785[k] = f_12 * osh_1128[k]
                    + f_1 * qsg0_955[k]
                    - f_2 * qsg1_955[k]
                    + f_3 * pc_y[k] * qsh_1338[k];

        t_1786[k] = f_20 * osh_1107[k]
                    + f_3 * pc_z[k] * qsh_1338[k];
    }

#pragma omp simd aligned(t_1787, t_1788, t_1789, pc_y, osh_1130, osh_1131, osh_1132, qsg0_957, \
                         qsg0_958, qsg0_959, qsg1_957, qsg1_958, qsg1_959, qsh_1340, qsh_1341, \
                         qsh_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1787[k] = f_12 * osh_1130[k]
                    + f_8 * qsg0_957[k]
                    - f_9 * qsg1_957[k]
                    + f_3 * pc_y[k] * qsh_1340[k];

        t_1788[k] = f_12 * osh_1131[k]
                    + f_6 * qsg0_958[k]
                    - f_7 * qsg1_958[k]
                    + f_3 * pc_y[k] * qsh_1341[k];

        t_1789[k] = f_12 * osh_1132[k]
                    + f_4 * qsg0_959[k]
                    - f_5 * qsg1_959[k]
                    + f_3 * pc_y[k] * qsh_1342[k];
    }

#pragma omp simd aligned(t_1790, t_1791, t_1792, t_1793, pa_y, pc_y, pc_z, osi0_1512, \
                         osh_1112, osh_1133, osh_1134, osi1_1512, qsg0_959, qsg1_959, \
                         qsh_1343, qsh_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1790[k] = f_12 * osh_1133[k]
                    + f_3 * pc_y[k] * qsh_1343[k];

        t_1791[k] = f_20 * osh_1112[k]
                    + f_1 * qsg0_959[k]
                    - f_2 * qsg1_959[k]
                    + f_3 * pc_z[k] * qsh_1343[k];

        t_1792[k] = pa_y[k] * osi0_1512[k]
                    - f_10 * pc_y[k] * osi1_1512[k];

        t_1793[k] = f_11 * osh_1134[k]
                    + f_3 * pc_y[k] * qsh_1344[k];
    }

#pragma omp simd aligned(t_1794, t_1795, t_1796, t_1797, pa_y, pc_y, pc_z, osi0_1515, \
                         osi0_1517, osh_1113, osh_1135, osh_1136, osi1_1515, osi1_1517, \
                         qsh_1344, qsh_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1794[k] = f_19 * osh_1113[k]
                    + f_3 * pc_z[k] * qsh_1344[k];

        t_1795[k] = pa_y[k] * osi0_1515[k]
                    + f_12 * osh_1135[k]
                    - f_10 * pc_y[k] * osi1_1515[k];

        t_1796[k] = f_11 * osh_1136[k]
                    + f_3 * pc_y[k] * qsh_1346[k];

        t_1797[k] = pa_y[k] * osi0_1517[k]
                    - f_10 * pc_y[k] * osi1_1517[k];
    }

#pragma omp simd aligned(t_1798, t_1799, t_1800, t_1801, pa_y, pc_y, pc_z, osi0_1518, \
                         osi0_1521, osh_1116, osh_1137, osh_1139, osi1_1518, osi1_1521, \
                         qsh_1347, qsh_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1798[k] = pa_y[k] * osi0_1518[k]
                    + f_13 * osh_1137[k]
                    - f_10 * pc_y[k] * osi1_1518[k];

        t_1799[k] = f_19 * osh_1116[k]
                    + f_3 * pc_z[k] * qsh_1347[k];

        t_1800[k] = f_11 * osh_1139[k]
                    + f_3 * pc_y[k] * qsh_1349[k];

        t_1801[k] = pa_y[k] * osi0_1521[k]
                    - f_10 * pc_y[k] * osi1_1521[k];
    }

#pragma omp simd aligned(t_1802, t_1803, t_1804, pa_y, pc_y, pc_z, osi0_1522, osi0_1524, \
                         osh_1119, osh_1140, osh_1142, osi1_1522, osi1_1524, \
                         qsh_1350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1802[k] = pa_y[k] * osi0_1522[k]
                    + f_14 * osh_1140[k]
                    - f_10 * pc_y[k] * osi1_1522[k];

        t_1803[k] = f_19 * osh_1119[k]
                    + f_3 * pc_z[k] * qsh_1350[k];

        t_1804[k] = pa_y[k] * osi0_1524[k]
                    + f_12 * osh_1142[k]
                    - f_10 * pc_y[k] * osi1_1524[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, t_1808, pa_y, pc_x, pc_y, osi0_1526, \
                         osh_1143, osh_1359, osh_1360, osi1_1526, qsh_1353, qsh_1359, \
                         qsh_1360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = f_11 * osh_1143[k]
                    + f_3 * pc_y[k] * qsh_1353[k];

        t_1806[k] = pa_y[k] * osi0_1526[k]
                    - f_10 * pc_y[k] * osi1_1526[k];

        t_1807[k] = f_12 * osh_1359[k]
                    + f_3 * pc_x[k] * qsh_1359[k];

        t_1808[k] = f_12 * osh_1360[k]
                    + f_3 * pc_x[k] * qsh_1360[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsg0, const size_t qsg1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 5.0 / q;
    const auto f_19 = 4.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_1809 = buffer.data(target + 1809);
    auto *t_1810 = buffer.data(target + 1810);
    auto *t_1811 = buffer.data(target + 1811);
    auto *t_1812 = buffer.data(target + 1812);
    auto *t_1813 = buffer.data(target + 1813);
    auto *t_1814 = buffer.data(target + 1814);
    auto *t_1815 = buffer.data(target + 1815);
    auto *t_1816 = buffer.data(target + 1816);
    auto *t_1817 = buffer.data(target + 1817);
    auto *t_1818 = buffer.data(target + 1818);
    auto *t_1819 = buffer.data(target + 1819);
    auto *t_1820 = buffer.data(target + 1820);
    auto *t_1821 = buffer.data(target + 1821);
    auto *t_1822 = buffer.data(target + 1822);
    auto *t_1823 = buffer.data(target + 1823);
    auto *t_1824 = buffer.data(target + 1824);
    auto *t_1825 = buffer.data(target + 1825);
    auto *t_1826 = buffer.data(target + 1826);
    auto *t_1827 = buffer.data(target + 1827);
    auto *t_1828 = buffer.data(target + 1828);
    auto *t_1829 = buffer.data(target + 1829);
    auto *t_1830 = buffer.data(target + 1830);
    auto *t_1831 = buffer.data(target + 1831);
    auto *t_1832 = buffer.data(target + 1832);
    auto *t_1833 = buffer.data(target + 1833);
    auto *t_1834 = buffer.data(target + 1834);
    auto *t_1835 = buffer.data(target + 1835);
    auto *t_1836 = buffer.data(target + 1836);
    auto *t_1837 = buffer.data(target + 1837);
    auto *t_1838 = buffer.data(target + 1838);
    auto *t_1839 = buffer.data(target + 1839);
    auto *t_1840 = buffer.data(target + 1840);
    auto *t_1841 = buffer.data(target + 1841);
    auto *t_1842 = buffer.data(target + 1842);
    auto *t_1843 = buffer.data(target + 1843);
    auto *t_1844 = buffer.data(target + 1844);
    auto *t_1845 = buffer.data(target + 1845);
    auto *t_1846 = buffer.data(target + 1846);
    auto *t_1847 = buffer.data(target + 1847);
    auto *t_1848 = buffer.data(target + 1848);
    auto *t_1849 = buffer.data(target + 1849);
    auto *t_1850 = buffer.data(target + 1850);
    auto *t_1851 = buffer.data(target + 1851);
    auto *t_1852 = buffer.data(target + 1852);
    auto *t_1853 = buffer.data(target + 1853);
    auto *t_1854 = buffer.data(target + 1854);
    auto *t_1855 = buffer.data(target + 1855);
    auto *t_1856 = buffer.data(target + 1856);
    auto *t_1857 = buffer.data(target + 1857);
    auto *t_1858 = buffer.data(target + 1858);
    auto *t_1859 = buffer.data(target + 1859);
    auto *t_1860 = buffer.data(target + 1860);
    auto *t_1861 = buffer.data(target + 1861);
    auto *t_1862 = buffer.data(target + 1862);
    auto *t_1863 = buffer.data(target + 1863);
    auto *t_1864 = buffer.data(target + 1864);
    auto *t_1865 = buffer.data(target + 1865);
    auto *t_1866 = buffer.data(target + 1866);
    auto *t_1867 = buffer.data(target + 1867);
    auto *t_1868 = buffer.data(target + 1868);
    auto *t_1869 = buffer.data(target + 1869);
    auto *t_1870 = buffer.data(target + 1870);
    auto *t_1871 = buffer.data(target + 1871);
    auto *t_1872 = buffer.data(target + 1872);
    auto *t_1873 = buffer.data(target + 1873);
    auto *t_1874 = buffer.data(target + 1874);
    auto *t_1875 = buffer.data(target + 1875);
    auto *t_1876 = buffer.data(target + 1876);
    auto *t_1877 = buffer.data(target + 1877);
    auto *t_1878 = buffer.data(target + 1878);
    auto *t_1879 = buffer.data(target + 1879);
    auto *t_1880 = buffer.data(target + 1880);
    auto *t_1881 = buffer.data(target + 1881);
    auto *t_1882 = buffer.data(target + 1882);
    auto *t_1883 = buffer.data(target + 1883);
    auto *t_1884 = buffer.data(target + 1884);
    auto *t_1885 = buffer.data(target + 1885);
    auto *t_1886 = buffer.data(target + 1886);
    auto *t_1887 = buffer.data(target + 1887);
    auto *t_1888 = buffer.data(target + 1888);
    auto *t_1889 = buffer.data(target + 1889);
    auto *t_1890 = buffer.data(target + 1890);
    auto *t_1891 = buffer.data(target + 1891);
    auto *t_1892 = buffer.data(target + 1892);
    auto *t_1893 = buffer.data(target + 1893);
    auto *t_1894 = buffer.data(target + 1894);
    auto *t_1895 = buffer.data(target + 1895);
    auto *t_1896 = buffer.data(target + 1896);
    auto *t_1897 = buffer.data(target + 1897);
    auto *t_1898 = buffer.data(target + 1898);
    auto *t_1899 = buffer.data(target + 1899);
    auto *t_1900 = buffer.data(target + 1900);
    auto *t_1901 = buffer.data(target + 1901);
    auto *t_1902 = buffer.data(target + 1902);
    auto *t_1903 = buffer.data(target + 1903);
    auto *t_1904 = buffer.data(target + 1904);
    auto *t_1905 = buffer.data(target + 1905);
    auto *t_1906 = buffer.data(target + 1906);
    auto *t_1907 = buffer.data(target + 1907);
    auto *t_1908 = buffer.data(target + 1908);
    auto *t_1909 = buffer.data(target + 1909);
    auto *t_1910 = buffer.data(target + 1910);
    auto *t_1911 = buffer.data(target + 1911);
    auto *t_1912 = buffer.data(target + 1912);
    auto *t_1913 = buffer.data(target + 1913);
    auto *t_1914 = buffer.data(target + 1914);
    auto *t_1915 = buffer.data(target + 1915);
    auto *t_1916 = buffer.data(target + 1916);
    auto *t_1917 = buffer.data(target + 1917);
    auto *t_1918 = buffer.data(target + 1918);
    auto *t_1919 = buffer.data(target + 1919);
    auto *t_1920 = buffer.data(target + 1920);
    auto *t_1921 = buffer.data(target + 1921);
    auto *t_1922 = buffer.data(target + 1922);
    auto *t_1923 = buffer.data(target + 1923);
    auto *t_1924 = buffer.data(target + 1924);
    auto *t_1925 = buffer.data(target + 1925);
    auto *t_1926 = buffer.data(target + 1926);
    auto *t_1927 = buffer.data(target + 1927);
    auto *t_1928 = buffer.data(target + 1928);
    auto *t_1929 = buffer.data(target + 1929);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_1539 = buffer.data(osi0 + 1539);
    const auto *osi0_1540 = buffer.data(osi0 + 1540);
    const auto *osi0_1543 = buffer.data(osi0 + 1543);
    const auto *osi0_1546 = buffer.data(osi0 + 1546);
    const auto *osi0_1550 = buffer.data(osi0 + 1550);
    const auto *osi0_1848 = buffer.data(osi0 + 1848);
    const auto *osi0_1851 = buffer.data(osi0 + 1851);
    const auto *osi0_1854 = buffer.data(osi0 + 1854);
    const auto *osi0_1858 = buffer.data(osi0 + 1858);
    const auto *osi0_1869 = buffer.data(osi0 + 1869);
    const auto *osi0_1871 = buffer.data(osi0 + 1871);
    const auto *osi0_1872 = buffer.data(osi0 + 1872);
    const auto *osi0_1873 = buffer.data(osi0 + 1873);
    const auto *osi0_1875 = buffer.data(osi0 + 1875);
    const auto *osi0_1881 = buffer.data(osi0 + 1881);
    const auto *osi0_1885 = buffer.data(osi0 + 1885);
    const auto *osi0_1888 = buffer.data(osi0 + 1888);
    const auto *osi0_1890 = buffer.data(osi0 + 1890);
    const auto *osi0_1897 = buffer.data(osi0 + 1897);
    const auto *osi0_1899 = buffer.data(osi0 + 1899);
    const auto *osi0_1900 = buffer.data(osi0 + 1900);
    const auto *osi0_1901 = buffer.data(osi0 + 1901);
    const auto *osi0_1903 = buffer.data(osi0 + 1903);
    const auto *osi0_1904 = buffer.data(osi0 + 1904);
    const auto *osi0_1907 = buffer.data(osi0 + 1907);
    const auto *osi0_1909 = buffer.data(osi0 + 1909);
    const auto *osi0_1910 = buffer.data(osi0 + 1910);
    const auto *osi0_1913 = buffer.data(osi0 + 1913);
    const auto *osi0_1914 = buffer.data(osi0 + 1914);
    const auto *osi0_1916 = buffer.data(osi0 + 1916);
    const auto *osi0_1918 = buffer.data(osi0 + 1918);
    const auto *osi0_1925 = buffer.data(osi0 + 1925);
    const auto *osi0_1927 = buffer.data(osi0 + 1927);
    const auto *osi0_1928 = buffer.data(osi0 + 1928);
    const auto *osi0_1929 = buffer.data(osi0 + 1929);

    const auto *osh_1128 = buffer.data(osh + 1128);
    const auto *osh_1134 = buffer.data(osh + 1134);
    const auto *osh_1149 = buffer.data(osh + 1149);
    const auto *osh_1151 = buffer.data(osh + 1151);
    const auto *osh_1152 = buffer.data(osh + 1152);
    const auto *osh_1153 = buffer.data(osh + 1153);
    const auto *osh_1154 = buffer.data(osh + 1154);
    const auto *osh_1155 = buffer.data(osh + 1155);
    const auto *osh_1158 = buffer.data(osh + 1158);
    const auto *osh_1160 = buffer.data(osh + 1160);
    const auto *osh_1161 = buffer.data(osh + 1161);
    const auto *osh_1164 = buffer.data(osh + 1164);
    const auto *osh_1170 = buffer.data(osh + 1170);
    const auto *osh_1175 = buffer.data(osh + 1175);
    const auto *osh_1176 = buffer.data(osh + 1176);
    const auto *osh_1178 = buffer.data(osh + 1178);
    const auto *osh_1179 = buffer.data(osh + 1179);
    const auto *osh_1181 = buffer.data(osh + 1181);
    const auto *osh_1182 = buffer.data(osh + 1182);
    const auto *osh_1185 = buffer.data(osh + 1185);
    const auto *osh_1191 = buffer.data(osh + 1191);
    const auto *osh_1196 = buffer.data(osh + 1196);
    const auto *osh_1197 = buffer.data(osh + 1197);
    const auto *osh_1199 = buffer.data(osh + 1199);
    const auto *osh_1202 = buffer.data(osh + 1202);
    const auto *osh_1206 = buffer.data(osh + 1206);
    const auto *osh_1361 = buffer.data(osh + 1361);
    const auto *osh_1362 = buffer.data(osh + 1362);
    const auto *osh_1363 = buffer.data(osh + 1363);
    const auto *osh_1364 = buffer.data(osh + 1364);
    const auto *osh_1365 = buffer.data(osh + 1365);
    const auto *osh_1370 = buffer.data(osh + 1370);
    const auto *osh_1374 = buffer.data(osh + 1374);
    const auto *osh_1379 = buffer.data(osh + 1379);
    const auto *osh_1380 = buffer.data(osh + 1380);
    const auto *osh_1381 = buffer.data(osh + 1381);
    const auto *osh_1382 = buffer.data(osh + 1382);
    const auto *osh_1383 = buffer.data(osh + 1383);
    const auto *osh_1385 = buffer.data(osh + 1385);
    const auto *osh_1386 = buffer.data(osh + 1386);
    const auto *osh_1389 = buffer.data(osh + 1389);
    const auto *osh_1392 = buffer.data(osh + 1392);
    const auto *osh_1396 = buffer.data(osh + 1396);
    const auto *osh_1401 = buffer.data(osh + 1401);
    const auto *osh_1403 = buffer.data(osh + 1403);
    const auto *osh_1404 = buffer.data(osh + 1404);
    const auto *osh_1405 = buffer.data(osh + 1405);
    const auto *osh_1406 = buffer.data(osh + 1406);
    const auto *osh_1412 = buffer.data(osh + 1412);
    const auto *osh_1416 = buffer.data(osh + 1416);
    const auto *osh_1419 = buffer.data(osh + 1419);
    const auto *osh_1421 = buffer.data(osh + 1421);
    const auto *osh_1422 = buffer.data(osh + 1422);
    const auto *osh_1423 = buffer.data(osh + 1423);
    const auto *osh_1424 = buffer.data(osh + 1424);
    const auto *osh_1425 = buffer.data(osh + 1425);
    const auto *osh_1426 = buffer.data(osh + 1426);
    const auto *osh_1427 = buffer.data(osh + 1427);
    const auto *osh_1428 = buffer.data(osh + 1428);
    const auto *osh_1431 = buffer.data(osh + 1431);
    const auto *osh_1433 = buffer.data(osh + 1433);
    const auto *osh_1434 = buffer.data(osh + 1434);
    const auto *osh_1437 = buffer.data(osh + 1437);
    const auto *osh_1438 = buffer.data(osh + 1438);
    const auto *osh_1440 = buffer.data(osh + 1440);
    const auto *osh_1442 = buffer.data(osh + 1442);
    const auto *osh_1443 = buffer.data(osh + 1443);
    const auto *osh_1444 = buffer.data(osh + 1444);
    const auto *osh_1445 = buffer.data(osh + 1445);
    const auto *osh_1446 = buffer.data(osh + 1446);
    const auto *osh_1447 = buffer.data(osh + 1447);
    const auto *osh_1448 = buffer.data(osh + 1448);

    const auto *osi1_1539 = buffer.data(osi1 + 1539);
    const auto *osi1_1540 = buffer.data(osi1 + 1540);
    const auto *osi1_1543 = buffer.data(osi1 + 1543);
    const auto *osi1_1546 = buffer.data(osi1 + 1546);
    const auto *osi1_1550 = buffer.data(osi1 + 1550);
    const auto *osi1_1848 = buffer.data(osi1 + 1848);
    const auto *osi1_1851 = buffer.data(osi1 + 1851);
    const auto *osi1_1854 = buffer.data(osi1 + 1854);
    const auto *osi1_1858 = buffer.data(osi1 + 1858);
    const auto *osi1_1869 = buffer.data(osi1 + 1869);
    const auto *osi1_1871 = buffer.data(osi1 + 1871);
    const auto *osi1_1872 = buffer.data(osi1 + 1872);
    const auto *osi1_1873 = buffer.data(osi1 + 1873);
    const auto *osi1_1875 = buffer.data(osi1 + 1875);
    const auto *osi1_1881 = buffer.data(osi1 + 1881);
    const auto *osi1_1885 = buffer.data(osi1 + 1885);
    const auto *osi1_1888 = buffer.data(osi1 + 1888);
    const auto *osi1_1890 = buffer.data(osi1 + 1890);
    const auto *osi1_1897 = buffer.data(osi1 + 1897);
    const auto *osi1_1899 = buffer.data(osi1 + 1899);
    const auto *osi1_1900 = buffer.data(osi1 + 1900);
    const auto *osi1_1901 = buffer.data(osi1 + 1901);
    const auto *osi1_1903 = buffer.data(osi1 + 1903);
    const auto *osi1_1904 = buffer.data(osi1 + 1904);
    const auto *osi1_1907 = buffer.data(osi1 + 1907);
    const auto *osi1_1909 = buffer.data(osi1 + 1909);
    const auto *osi1_1910 = buffer.data(osi1 + 1910);
    const auto *osi1_1913 = buffer.data(osi1 + 1913);
    const auto *osi1_1914 = buffer.data(osi1 + 1914);
    const auto *osi1_1916 = buffer.data(osi1 + 1916);
    const auto *osi1_1918 = buffer.data(osi1 + 1918);
    const auto *osi1_1925 = buffer.data(osi1 + 1925);
    const auto *osi1_1927 = buffer.data(osi1 + 1927);
    const auto *osi1_1928 = buffer.data(osi1 + 1928);
    const auto *osi1_1929 = buffer.data(osi1 + 1929);

    const auto *qsg0_970 = buffer.data(qsg0 + 970);
    const auto *qsg0_972 = buffer.data(qsg0 + 972);
    const auto *qsg0_973 = buffer.data(qsg0 + 973);
    const auto *qsg0_974 = buffer.data(qsg0 + 974);
    const auto *qsg0_975 = buffer.data(qsg0 + 975);
    const auto *qsg0_976 = buffer.data(qsg0 + 976);
    const auto *qsg0_977 = buffer.data(qsg0 + 977);
    const auto *qsg0_978 = buffer.data(qsg0 + 978);
    const auto *qsg0_979 = buffer.data(qsg0 + 979);
    const auto *qsg0_980 = buffer.data(qsg0 + 980);
    const auto *qsg0_984 = buffer.data(qsg0 + 984);
    const auto *qsg0_985 = buffer.data(qsg0 + 985);
    const auto *qsg0_986 = buffer.data(qsg0 + 986);
    const auto *qsg0_987 = buffer.data(qsg0 + 987);
    const auto *qsg0_988 = buffer.data(qsg0 + 988);
    const auto *qsg0_989 = buffer.data(qsg0 + 989);
    const auto *qsg0_990 = buffer.data(qsg0 + 990);
    const auto *qsg0_992 = buffer.data(qsg0 + 992);
    const auto *qsg0_993 = buffer.data(qsg0 + 993);
    const auto *qsg0_995 = buffer.data(qsg0 + 995);

    const auto *qsg1_970 = buffer.data(qsg1 + 970);
    const auto *qsg1_972 = buffer.data(qsg1 + 972);
    const auto *qsg1_973 = buffer.data(qsg1 + 973);
    const auto *qsg1_974 = buffer.data(qsg1 + 974);
    const auto *qsg1_975 = buffer.data(qsg1 + 975);
    const auto *qsg1_976 = buffer.data(qsg1 + 976);
    const auto *qsg1_977 = buffer.data(qsg1 + 977);
    const auto *qsg1_978 = buffer.data(qsg1 + 978);
    const auto *qsg1_979 = buffer.data(qsg1 + 979);
    const auto *qsg1_980 = buffer.data(qsg1 + 980);
    const auto *qsg1_984 = buffer.data(qsg1 + 984);
    const auto *qsg1_985 = buffer.data(qsg1 + 985);
    const auto *qsg1_986 = buffer.data(qsg1 + 986);
    const auto *qsg1_987 = buffer.data(qsg1 + 987);
    const auto *qsg1_988 = buffer.data(qsg1 + 988);
    const auto *qsg1_989 = buffer.data(qsg1 + 989);
    const auto *qsg1_990 = buffer.data(qsg1 + 990);
    const auto *qsg1_992 = buffer.data(qsg1 + 992);
    const auto *qsg1_993 = buffer.data(qsg1 + 993);
    const auto *qsg1_995 = buffer.data(qsg1 + 995);

    const auto *qsh_1359 = buffer.data(qsh + 1359);
    const auto *qsh_1361 = buffer.data(qsh + 1361);
    const auto *qsh_1362 = buffer.data(qsh + 1362);
    const auto *qsh_1363 = buffer.data(qsh + 1363);
    const auto *qsh_1364 = buffer.data(qsh + 1364);
    const auto *qsh_1365 = buffer.data(qsh + 1365);
    const auto *qsh_1366 = buffer.data(qsh + 1366);
    const auto *qsh_1367 = buffer.data(qsh + 1367);
    const auto *qsh_1368 = buffer.data(qsh + 1368);
    const auto *qsh_1369 = buffer.data(qsh + 1369);
    const auto *qsh_1370 = buffer.data(qsh + 1370);
    const auto *qsh_1371 = buffer.data(qsh + 1371);
    const auto *qsh_1372 = buffer.data(qsh + 1372);
    const auto *qsh_1373 = buffer.data(qsh + 1373);
    const auto *qsh_1374 = buffer.data(qsh + 1374);
    const auto *qsh_1379 = buffer.data(qsh + 1379);
    const auto *qsh_1380 = buffer.data(qsh + 1380);
    const auto *qsh_1381 = buffer.data(qsh + 1381);
    const auto *qsh_1382 = buffer.data(qsh + 1382);
    const auto *qsh_1383 = buffer.data(qsh + 1383);
    const auto *qsh_1384 = buffer.data(qsh + 1384);
    const auto *qsh_1385 = buffer.data(qsh + 1385);
    const auto *qsh_1386 = buffer.data(qsh + 1386);
    const auto *qsh_1387 = buffer.data(qsh + 1387);
    const auto *qsh_1388 = buffer.data(qsh + 1388);
    const auto *qsh_1389 = buffer.data(qsh + 1389);
    const auto *qsh_1391 = buffer.data(qsh + 1391);
    const auto *qsh_1392 = buffer.data(qsh + 1392);
    const auto *qsh_1393 = buffer.data(qsh + 1393);
    const auto *qsh_1395 = buffer.data(qsh + 1395);
    const auto *qsh_1396 = buffer.data(qsh + 1396);
    const auto *qsh_1401 = buffer.data(qsh + 1401);
    const auto *qsh_1403 = buffer.data(qsh + 1403);
    const auto *qsh_1404 = buffer.data(qsh + 1404);
    const auto *qsh_1405 = buffer.data(qsh + 1405);
    const auto *qsh_1406 = buffer.data(qsh + 1406);
    const auto *qsh_1407 = buffer.data(qsh + 1407);
    const auto *qsh_1409 = buffer.data(qsh + 1409);
    const auto *qsh_1410 = buffer.data(qsh + 1410);
    const auto *qsh_1412 = buffer.data(qsh + 1412);
    const auto *qsh_1413 = buffer.data(qsh + 1413);
    const auto *qsh_1416 = buffer.data(qsh + 1416);
    const auto *qsh_1422 = buffer.data(qsh + 1422);
    const auto *qsh_1423 = buffer.data(qsh + 1423);
    const auto *qsh_1424 = buffer.data(qsh + 1424);
    const auto *qsh_1425 = buffer.data(qsh + 1425);
    const auto *qsh_1426 = buffer.data(qsh + 1426);
    const auto *qsh_1427 = buffer.data(qsh + 1427);
    const auto *qsh_1428 = buffer.data(qsh + 1428);
    const auto *qsh_1430 = buffer.data(qsh + 1430);
    const auto *qsh_1431 = buffer.data(qsh + 1431);
    const auto *qsh_1433 = buffer.data(qsh + 1433);
    const auto *qsh_1434 = buffer.data(qsh + 1434);
    const auto *qsh_1437 = buffer.data(qsh + 1437);
    const auto *qsh_1443 = buffer.data(qsh + 1443);
    const auto *qsh_1444 = buffer.data(qsh + 1444);
    const auto *qsh_1445 = buffer.data(qsh + 1445);
    const auto *qsh_1446 = buffer.data(qsh + 1446);
    const auto *qsh_1447 = buffer.data(qsh + 1447);
    const auto *qsh_1448 = buffer.data(qsh + 1448);

#pragma omp simd aligned(t_1809, t_1810, t_1811, t_1812, pc_x, osh_1361, osh_1362, osh_1363, \
                         osh_1364, qsh_1361, qsh_1362, qsh_1363, \
                         qsh_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1809[k] = f_12 * osh_1361[k]
                    + f_3 * pc_x[k] * qsh_1361[k];

        t_1810[k] = f_12 * osh_1362[k]
                    + f_3 * pc_x[k] * qsh_1362[k];

        t_1811[k] = f_12 * osh_1363[k]
                    + f_3 * pc_x[k] * qsh_1363[k];

        t_1812[k] = f_12 * osh_1364[k]
                    + f_3 * pc_x[k] * qsh_1364[k];
    }

#pragma omp simd aligned(t_1813, t_1814, t_1815, pc_y, pc_z, osh_1128, osh_1149, osh_1151, \
                         qsg0_970, qsg0_972, qsg1_970, qsg1_972, qsh_1359, \
                         qsh_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1813[k] = f_11 * osh_1149[k]
                    + f_1 * qsg0_970[k]
                    - f_2 * qsg1_970[k]
                    + f_3 * pc_y[k] * qsh_1359[k];

        t_1814[k] = f_19 * osh_1128[k]
                    + f_3 * pc_z[k] * qsh_1359[k];

        t_1815[k] = f_11 * osh_1151[k]
                    + f_8 * qsg0_972[k]
                    - f_9 * qsg1_972[k]
                    + f_3 * pc_y[k] * qsh_1361[k];
    }

#pragma omp simd aligned(t_1816, t_1817, t_1818, pc_y, osh_1152, osh_1153, osh_1154, qsg0_973, \
                         qsg0_974, qsg1_973, qsg1_974, qsh_1362, qsh_1363, \
                         qsh_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1816[k] = f_11 * osh_1152[k]
                    + f_6 * qsg0_973[k]
                    - f_7 * qsg1_973[k]
                    + f_3 * pc_y[k] * qsh_1362[k];

        t_1817[k] = f_11 * osh_1153[k]
                    + f_4 * qsg0_974[k]
                    - f_5 * qsg1_974[k]
                    + f_3 * pc_y[k] * qsh_1363[k];

        t_1818[k] = f_11 * osh_1154[k]
                    + f_3 * pc_y[k] * qsh_1364[k];
    }

#pragma omp simd aligned(t_1819, t_1820, t_1821, t_1822, pa_y, pc_x, pc_y, pc_z, osi0_1539, \
                         osh_1134, osh_1365, osi1_1539, qsg0_975, qsg1_975, \
                         qsh_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1819[k] = pa_y[k] * osi0_1539[k]
                    - f_10 * pc_y[k] * osi1_1539[k];

        t_1820[k] = f_12 * osh_1365[k]
                    + f_1 * qsg0_975[k]
                    - f_2 * qsg1_975[k]
                    + f_3 * pc_x[k] * qsh_1365[k];

        t_1821[k] = f_3 * pc_y[k] * qsh_1365[k];

        t_1822[k] = f_18 * osh_1134[k]
                    + f_3 * pc_z[k] * qsh_1365[k];
    }

#pragma omp simd aligned(t_1823, t_1824, t_1825, pc_x, pc_y, osh_1370, qsg0_975, qsg0_980, \
                         qsg1_975, qsg1_980, qsh_1366, qsh_1367, \
                         qsh_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1823[k] = f_4 * qsg0_975[k]
                    - f_5 * qsg1_975[k]
                    + f_3 * pc_y[k] * qsh_1366[k];

        t_1824[k] = f_3 * pc_y[k] * qsh_1367[k];

        t_1825[k] = f_12 * osh_1370[k]
                    + f_8 * qsg0_980[k]
                    - f_9 * qsg1_980[k]
                    + f_3 * pc_x[k] * qsh_1370[k];
    }

#pragma omp simd aligned(t_1826, t_1827, t_1828, pc_y, qsg0_976, qsg0_977, qsg1_976, qsg1_977, \
                         qsh_1368, qsh_1369, qsh_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1826[k] = f_6 * qsg0_976[k]
                    - f_7 * qsg1_976[k]
                    + f_3 * pc_y[k] * qsh_1368[k];

        t_1827[k] = f_4 * qsg0_977[k]
                    - f_5 * qsg1_977[k]
                    + f_3 * pc_y[k] * qsh_1369[k];

        t_1828[k] = f_3 * pc_y[k] * qsh_1370[k];
    }

#pragma omp simd aligned(t_1829, t_1830, t_1831, pc_x, pc_y, osh_1374, qsg0_978, qsg0_979, \
                         qsg0_984, qsg1_978, qsg1_979, qsg1_984, qsh_1371, qsh_1372, \
                         qsh_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1829[k] = f_12 * osh_1374[k]
                    + f_6 * qsg0_984[k]
                    - f_7 * qsg1_984[k]
                    + f_3 * pc_x[k] * qsh_1374[k];

        t_1830[k] = f_8 * qsg0_978[k]
                    - f_9 * qsg1_978[k]
                    + f_3 * pc_y[k] * qsh_1371[k];

        t_1831[k] = f_6 * qsg0_979[k]
                    - f_7 * qsg1_979[k]
                    + f_3 * pc_y[k] * qsh_1372[k];
    }

#pragma omp simd aligned(t_1832, t_1833, t_1834, t_1835, pc_x, pc_y, osh_1379, osh_1380, \
                         qsg0_980, qsg0_989, qsg1_980, qsg1_989, qsh_1373, qsh_1374, qsh_1379, \
                         qsh_1380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1832[k] = f_4 * qsg0_980[k]
                    - f_5 * qsg1_980[k]
                    + f_3 * pc_y[k] * qsh_1373[k];

        t_1833[k] = f_3 * pc_y[k] * qsh_1374[k];

        t_1834[k] = f_12 * osh_1379[k]
                    + f_4 * qsg0_989[k]
                    - f_5 * qsg1_989[k]
                    + f_3 * pc_x[k] * qsh_1379[k];

        t_1835[k] = f_12 * osh_1380[k]
                    + f_3 * pc_x[k] * qsh_1380[k];
    }

#pragma omp simd aligned(t_1836, t_1837, t_1838, t_1839, t_1840, pc_x, pc_y, osh_1381, \
                         osh_1382, osh_1383, osh_1385, qsh_1379, qsh_1381, qsh_1382, qsh_1383, \
                         qsh_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1836[k] = f_12 * osh_1381[k]
                    + f_3 * pc_x[k] * qsh_1381[k];

        t_1837[k] = f_12 * osh_1382[k]
                    + f_3 * pc_x[k] * qsh_1382[k];

        t_1838[k] = f_12 * osh_1383[k]
                    + f_3 * pc_x[k] * qsh_1383[k];

        t_1839[k] = f_3 * pc_y[k] * qsh_1379[k];

        t_1840[k] = f_12 * osh_1385[k]
                    + f_3 * pc_x[k] * qsh_1385[k];
    }

#pragma omp simd aligned(t_1841, t_1842, t_1843, pc_y, qsg0_985, qsg0_986, qsg0_987, qsg1_985, \
                         qsg1_986, qsg1_987, qsh_1380, qsh_1381, \
                         qsh_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1841[k] = f_1 * qsg0_985[k]
                    - f_2 * qsg1_985[k]
                    + f_3 * pc_y[k] * qsh_1380[k];

        t_1842[k] = f_16 * qsg0_986[k]
                    - f_17 * qsg1_986[k]
                    + f_3 * pc_y[k] * qsh_1381[k];

        t_1843[k] = f_8 * qsg0_987[k]
                    - f_9 * qsg1_987[k]
                    + f_3 * pc_y[k] * qsh_1382[k];
    }

#pragma omp simd aligned(t_1844, t_1845, t_1846, t_1847, pc_y, pc_z, osh_1154, qsg0_988, \
                         qsg0_989, qsg1_988, qsg1_989, qsh_1383, qsh_1384, \
                         qsh_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1844[k] = f_6 * qsg0_988[k]
                    - f_7 * qsg1_988[k]
                    + f_3 * pc_y[k] * qsh_1383[k];

        t_1845[k] = f_4 * qsg0_989[k]
                    - f_5 * qsg1_989[k]
                    + f_3 * pc_y[k] * qsh_1384[k];

        t_1846[k] = f_3 * pc_y[k] * qsh_1385[k];

        t_1847[k] = f_18 * osh_1154[k]
                    + f_1 * qsg0_989[k]
                    - f_2 * qsg1_989[k]
                    + f_3 * pc_z[k] * qsh_1385[k];
    }

#pragma omp simd aligned(t_1848, t_1849, t_1850, t_1851, pa_x, pc_x, pc_y, pc_z, osi0_1848, \
                         osi0_1851, osh_1155, osh_1386, osh_1389, osi1_1848, osi1_1851, \
                         qsh_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1848[k] = pa_x[k] * osi0_1848[k]
                    + f_23 * osh_1386[k]
                    - f_10 * pc_x[k] * osi1_1848[k];

        t_1849[k] = f_15 * osh_1155[k]
                    + f_3 * pc_y[k] * qsh_1386[k];

        t_1850[k] = f_3 * pc_z[k] * qsh_1386[k];

        t_1851[k] = pa_x[k] * osi0_1851[k]
                    + f_14 * osh_1389[k]
                    - f_10 * pc_x[k] * osi1_1851[k];
    }

#pragma omp simd aligned(t_1852, t_1853, t_1854, t_1855, pa_x, pc_x, pc_z, osi0_1854, \
                         osh_1392, osi1_1854, qsg0_990, qsg1_990, qsh_1387, qsh_1388, \
                         qsh_1389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1852[k] = f_3 * pc_z[k] * qsh_1387[k];

        t_1853[k] = f_4 * qsg0_990[k]
                    - f_5 * qsg1_990[k]
                    + f_3 * pc_z[k] * qsh_1388[k];

        t_1854[k] = pa_x[k] * osi0_1854[k]
                    + f_13 * osh_1392[k]
                    - f_10 * pc_x[k] * osi1_1854[k];

        t_1855[k] = f_3 * pc_z[k] * qsh_1389[k];
    }

#pragma omp simd aligned(t_1856, t_1857, t_1858, t_1859, pa_x, pc_x, pc_y, pc_z, osi0_1858, \
                         osh_1160, osh_1396, osi1_1858, qsg0_992, qsg1_992, qsh_1391, \
                         qsh_1392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1856[k] = f_15 * osh_1160[k]
                    + f_3 * pc_y[k] * qsh_1391[k];

        t_1857[k] = f_6 * qsg0_992[k]
                    - f_7 * qsg1_992[k]
                    + f_3 * pc_z[k] * qsh_1391[k];

        t_1858[k] = pa_x[k] * osi0_1858[k]
                    + f_12 * osh_1396[k]
                    - f_10 * pc_x[k] * osi1_1858[k];

        t_1859[k] = f_3 * pc_z[k] * qsh_1392[k];
    }

#pragma omp simd aligned(t_1860, t_1861, t_1862, t_1863, pc_x, pc_y, pc_z, osh_1164, osh_1401, \
                         qsg0_993, qsg0_995, qsg1_993, qsg1_995, qsh_1393, qsh_1395, \
                         qsh_1401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1860[k] = f_4 * qsg0_993[k]
                    - f_5 * qsg1_993[k]
                    + f_3 * pc_z[k] * qsh_1393[k];

        t_1861[k] = f_15 * osh_1164[k]
                    + f_3 * pc_y[k] * qsh_1395[k];

        t_1862[k] = f_8 * qsg0_995[k]
                    - f_9 * qsg1_995[k]
                    + f_3 * pc_z[k] * qsh_1395[k];

        t_1863[k] = f_11 * osh_1401[k]
                    + f_3 * pc_x[k] * qsh_1401[k];
    }

#pragma omp simd aligned(t_1864, t_1865, t_1866, t_1867, t_1868, pc_x, pc_z, osh_1403, \
                         osh_1404, osh_1405, osh_1406, qsh_1396, qsh_1403, qsh_1404, qsh_1405, \
                         qsh_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1864[k] = f_3 * pc_z[k] * qsh_1396[k];

        t_1865[k] = f_11 * osh_1403[k]
                    + f_3 * pc_x[k] * qsh_1403[k];

        t_1866[k] = f_11 * osh_1404[k]
                    + f_3 * pc_x[k] * qsh_1404[k];

        t_1867[k] = f_11 * osh_1405[k]
                    + f_3 * pc_x[k] * qsh_1405[k];

        t_1868[k] = f_11 * osh_1406[k]
                    + f_3 * pc_x[k] * qsh_1406[k];
    }

#pragma omp simd aligned(t_1869, t_1870, t_1871, t_1872, pa_x, pc_x, pc_z, osi0_1869, \
                         osi0_1871, osi0_1872, osi1_1869, osi1_1871, osi1_1872, \
                         qsh_1401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1869[k] = pa_x[k] * osi0_1869[k]
                    - f_10 * pc_x[k] * osi1_1869[k];

        t_1870[k] = f_3 * pc_z[k] * qsh_1401[k];

        t_1871[k] = pa_x[k] * osi0_1871[k]
                    - f_10 * pc_x[k] * osi1_1871[k];

        t_1872[k] = pa_x[k] * osi0_1872[k]
                    - f_10 * pc_x[k] * osi1_1872[k];
    }

#pragma omp simd aligned(t_1873, t_1874, t_1875, pa_x, pc_x, pc_y, osi0_1873, osi0_1875, \
                         osh_1175, osi1_1873, osi1_1875, qsh_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1873[k] = pa_x[k] * osi0_1873[k]
                    - f_10 * pc_x[k] * osi1_1873[k];

        t_1874[k] = f_15 * osh_1175[k]
                    + f_3 * pc_y[k] * qsh_1406[k];

        t_1875[k] = pa_x[k] * osi0_1875[k]
                    - f_10 * pc_x[k] * osi1_1875[k];
    }

#pragma omp simd aligned(t_1876, t_1877, t_1878, t_1879, pa_z, pc_y, pc_z, osi0_1540, \
                         osi0_1543, osh_1155, osh_1176, osi1_1540, osi1_1543, \
                         qsh_1407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1876[k] = pa_z[k] * osi0_1540[k]
                    - f_10 * pc_z[k] * osi1_1540[k];

        t_1877[k] = f_18 * osh_1176[k]
                    + f_3 * pc_y[k] * qsh_1407[k];

        t_1878[k] = f_11 * osh_1155[k]
                    + f_3 * pc_z[k] * qsh_1407[k];

        t_1879[k] = pa_z[k] * osi0_1543[k]
                    - f_10 * pc_z[k] * osi1_1543[k];
    }

#pragma omp simd aligned(t_1880, t_1881, t_1882, pa_x, pa_z, pc_x, pc_y, pc_z, osi0_1546, \
                         osi0_1881, osh_1178, osh_1412, osi1_1546, osi1_1881, \
                         qsh_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1880[k] = f_18 * osh_1178[k]
                    + f_3 * pc_y[k] * qsh_1409[k];

        t_1881[k] = pa_x[k] * osi0_1881[k]
                    + f_14 * osh_1412[k]
                    - f_10 * pc_x[k] * osi1_1881[k];

        t_1882[k] = pa_z[k] * osi0_1546[k]
                    - f_10 * pc_z[k] * osi1_1546[k];
    }

#pragma omp simd aligned(t_1883, t_1884, t_1885, pa_x, pc_x, pc_y, pc_z, osi0_1885, osh_1158, \
                         osh_1181, osh_1416, osi1_1885, qsh_1410, \
                         qsh_1412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1883[k] = f_11 * osh_1158[k]
                    + f_3 * pc_z[k] * qsh_1410[k];

        t_1884[k] = f_18 * osh_1181[k]
                    + f_3 * pc_y[k] * qsh_1412[k];

        t_1885[k] = pa_x[k] * osi0_1885[k]
                    + f_13 * osh_1416[k]
                    - f_10 * pc_x[k] * osi1_1885[k];
    }

#pragma omp simd aligned(t_1886, t_1887, t_1888, pa_x, pa_z, pc_x, pc_z, osi0_1550, osi0_1888, \
                         osh_1161, osh_1419, osi1_1550, osi1_1888, \
                         qsh_1413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1886[k] = pa_z[k] * osi0_1550[k]
                    - f_10 * pc_z[k] * osi1_1550[k];

        t_1887[k] = f_11 * osh_1161[k]
                    + f_3 * pc_z[k] * qsh_1413[k];

        t_1888[k] = pa_x[k] * osi0_1888[k]
                    + f_12 * osh_1419[k]
                    - f_10 * pc_x[k] * osi1_1888[k];
    }

#pragma omp simd aligned(t_1889, t_1890, t_1891, t_1892, pa_x, pc_x, pc_y, osi0_1890, \
                         osh_1185, osh_1421, osh_1422, osh_1423, osi1_1890, qsh_1416, \
                         qsh_1422, qsh_1423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1889[k] = f_18 * osh_1185[k]
                    + f_3 * pc_y[k] * qsh_1416[k];

        t_1890[k] = pa_x[k] * osi0_1890[k]
                    + f_12 * osh_1421[k]
                    - f_10 * pc_x[k] * osi1_1890[k];

        t_1891[k] = f_11 * osh_1422[k]
                    + f_3 * pc_x[k] * qsh_1422[k];

        t_1892[k] = f_11 * osh_1423[k]
                    + f_3 * pc_x[k] * qsh_1423[k];
    }

#pragma omp simd aligned(t_1893, t_1894, t_1895, t_1896, pc_x, osh_1424, osh_1425, osh_1426, \
                         osh_1427, qsh_1424, qsh_1425, qsh_1426, \
                         qsh_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1893[k] = f_11 * osh_1424[k]
                    + f_3 * pc_x[k] * qsh_1424[k];

        t_1894[k] = f_11 * osh_1425[k]
                    + f_3 * pc_x[k] * qsh_1425[k];

        t_1895[k] = f_11 * osh_1426[k]
                    + f_3 * pc_x[k] * qsh_1426[k];

        t_1896[k] = f_11 * osh_1427[k]
                    + f_3 * pc_x[k] * qsh_1427[k];
    }

#pragma omp simd aligned(t_1897, t_1898, t_1899, t_1900, pa_x, pc_x, pc_z, osi0_1897, \
                         osi0_1899, osi0_1900, osh_1170, osi1_1897, osi1_1899, osi1_1900, \
                         qsh_1422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1897[k] = pa_x[k] * osi0_1897[k]
                    - f_10 * pc_x[k] * osi1_1897[k];

        t_1898[k] = f_11 * osh_1170[k]
                    + f_3 * pc_z[k] * qsh_1422[k];

        t_1899[k] = pa_x[k] * osi0_1899[k]
                    - f_10 * pc_x[k] * osi1_1899[k];

        t_1900[k] = pa_x[k] * osi0_1900[k]
                    - f_10 * pc_x[k] * osi1_1900[k];
    }

#pragma omp simd aligned(t_1901, t_1902, t_1903, t_1904, pa_x, pc_x, pc_y, osi0_1901, \
                         osi0_1903, osi0_1904, osh_1196, osh_1428, osi1_1901, osi1_1903, \
                         osi1_1904, qsh_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1901[k] = pa_x[k] * osi0_1901[k]
                    - f_10 * pc_x[k] * osi1_1901[k];

        t_1902[k] = f_18 * osh_1196[k]
                    + f_3 * pc_y[k] * qsh_1427[k];

        t_1903[k] = pa_x[k] * osi0_1903[k]
                    - f_10 * pc_x[k] * osi1_1903[k];

        t_1904[k] = pa_x[k] * osi0_1904[k]
                    + f_23 * osh_1428[k]
                    - f_10 * pc_x[k] * osi1_1904[k];
    }

#pragma omp simd aligned(t_1905, t_1906, t_1907, t_1908, pa_x, pc_x, pc_y, pc_z, osi0_1907, \
                         osh_1176, osh_1197, osh_1199, osh_1431, osi1_1907, qsh_1428, \
                         qsh_1430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1905[k] = f_19 * osh_1197[k]
                    + f_3 * pc_y[k] * qsh_1428[k];

        t_1906[k] = f_12 * osh_1176[k]
                    + f_3 * pc_z[k] * qsh_1428[k];

        t_1907[k] = pa_x[k] * osi0_1907[k]
                    + f_14 * osh_1431[k]
                    - f_10 * pc_x[k] * osi1_1907[k];

        t_1908[k] = f_19 * osh_1199[k]
                    + f_3 * pc_y[k] * qsh_1430[k];
    }

#pragma omp simd aligned(t_1909, t_1910, t_1911, pa_x, pc_x, pc_z, osi0_1909, osi0_1910, \
                         osh_1179, osh_1433, osh_1434, osi1_1909, osi1_1910, \
                         qsh_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1909[k] = pa_x[k] * osi0_1909[k]
                    + f_14 * osh_1433[k]
                    - f_10 * pc_x[k] * osi1_1909[k];

        t_1910[k] = pa_x[k] * osi0_1910[k]
                    + f_13 * osh_1434[k]
                    - f_10 * pc_x[k] * osi1_1910[k];

        t_1911[k] = f_12 * osh_1179[k]
                    + f_3 * pc_z[k] * qsh_1431[k];
    }

#pragma omp simd aligned(t_1912, t_1913, t_1914, pa_x, pc_x, pc_y, osi0_1913, osi0_1914, \
                         osh_1202, osh_1437, osh_1438, osi1_1913, osi1_1914, \
                         qsh_1433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1912[k] = f_19 * osh_1202[k]
                    + f_3 * pc_y[k] * qsh_1433[k];

        t_1913[k] = pa_x[k] * osi0_1913[k]
                    + f_13 * osh_1437[k]
                    - f_10 * pc_x[k] * osi1_1913[k];

        t_1914[k] = pa_x[k] * osi0_1914[k]
                    + f_12 * osh_1438[k]
                    - f_10 * pc_x[k] * osi1_1914[k];
    }

#pragma omp simd aligned(t_1915, t_1916, t_1917, pa_x, pc_x, pc_y, pc_z, osi0_1916, osh_1182, \
                         osh_1206, osh_1440, osi1_1916, qsh_1434, \
                         qsh_1437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1915[k] = f_12 * osh_1182[k]
                    + f_3 * pc_z[k] * qsh_1434[k];

        t_1916[k] = pa_x[k] * osi0_1916[k]
                    + f_12 * osh_1440[k]
                    - f_10 * pc_x[k] * osi1_1916[k];

        t_1917[k] = f_19 * osh_1206[k]
                    + f_3 * pc_y[k] * qsh_1437[k];
    }

#pragma omp simd aligned(t_1918, t_1919, t_1920, t_1921, pa_x, pc_x, osi0_1918, osh_1442, \
                         osh_1443, osh_1444, osh_1445, osi1_1918, qsh_1443, qsh_1444, \
                         qsh_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1918[k] = pa_x[k] * osi0_1918[k]
                    + f_12 * osh_1442[k]
                    - f_10 * pc_x[k] * osi1_1918[k];

        t_1919[k] = f_11 * osh_1443[k]
                    + f_3 * pc_x[k] * qsh_1443[k];

        t_1920[k] = f_11 * osh_1444[k]
                    + f_3 * pc_x[k] * qsh_1444[k];

        t_1921[k] = f_11 * osh_1445[k]
                    + f_3 * pc_x[k] * qsh_1445[k];
    }

#pragma omp simd aligned(t_1922, t_1923, t_1924, t_1925, pa_x, pc_x, osi0_1925, osh_1446, \
                         osh_1447, osh_1448, osi1_1925, qsh_1446, qsh_1447, \
                         qsh_1448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1922[k] = f_11 * osh_1446[k]
                    + f_3 * pc_x[k] * qsh_1446[k];

        t_1923[k] = f_11 * osh_1447[k]
                    + f_3 * pc_x[k] * qsh_1447[k];

        t_1924[k] = f_11 * osh_1448[k]
                    + f_3 * pc_x[k] * qsh_1448[k];

        t_1925[k] = pa_x[k] * osi0_1925[k]
                    - f_10 * pc_x[k] * osi1_1925[k];
    }

#pragma omp simd aligned(t_1926, t_1927, t_1928, t_1929, pa_x, pc_x, pc_z, osi0_1927, \
                         osi0_1928, osi0_1929, osh_1191, osi1_1927, osi1_1928, osi1_1929, \
                         qsh_1443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1926[k] = f_12 * osh_1191[k]
                    + f_3 * pc_z[k] * qsh_1443[k];

        t_1927[k] = pa_x[k] * osi0_1927[k]
                    - f_10 * pc_x[k] * osi1_1927[k];

        t_1928[k] = pa_x[k] * osi0_1928[k]
                    - f_10 * pc_x[k] * osi1_1928[k];

        t_1929[k] = pa_x[k] * osi0_1929[k]
                    - f_10 * pc_x[k] * osi1_1929[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_19 = 4.5 / q;
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_1930 = buffer.data(target + 1930);
    auto *t_1931 = buffer.data(target + 1931);
    auto *t_1932 = buffer.data(target + 1932);
    auto *t_1933 = buffer.data(target + 1933);
    auto *t_1934 = buffer.data(target + 1934);
    auto *t_1935 = buffer.data(target + 1935);
    auto *t_1936 = buffer.data(target + 1936);
    auto *t_1937 = buffer.data(target + 1937);
    auto *t_1938 = buffer.data(target + 1938);
    auto *t_1939 = buffer.data(target + 1939);
    auto *t_1940 = buffer.data(target + 1940);
    auto *t_1941 = buffer.data(target + 1941);
    auto *t_1942 = buffer.data(target + 1942);
    auto *t_1943 = buffer.data(target + 1943);
    auto *t_1944 = buffer.data(target + 1944);
    auto *t_1945 = buffer.data(target + 1945);
    auto *t_1946 = buffer.data(target + 1946);
    auto *t_1947 = buffer.data(target + 1947);
    auto *t_1948 = buffer.data(target + 1948);
    auto *t_1949 = buffer.data(target + 1949);
    auto *t_1950 = buffer.data(target + 1950);
    auto *t_1951 = buffer.data(target + 1951);
    auto *t_1952 = buffer.data(target + 1952);
    auto *t_1953 = buffer.data(target + 1953);
    auto *t_1954 = buffer.data(target + 1954);
    auto *t_1955 = buffer.data(target + 1955);
    auto *t_1956 = buffer.data(target + 1956);
    auto *t_1957 = buffer.data(target + 1957);
    auto *t_1958 = buffer.data(target + 1958);
    auto *t_1959 = buffer.data(target + 1959);
    auto *t_1960 = buffer.data(target + 1960);
    auto *t_1961 = buffer.data(target + 1961);
    auto *t_1962 = buffer.data(target + 1962);
    auto *t_1963 = buffer.data(target + 1963);
    auto *t_1964 = buffer.data(target + 1964);
    auto *t_1965 = buffer.data(target + 1965);
    auto *t_1966 = buffer.data(target + 1966);
    auto *t_1967 = buffer.data(target + 1967);
    auto *t_1968 = buffer.data(target + 1968);
    auto *t_1969 = buffer.data(target + 1969);
    auto *t_1970 = buffer.data(target + 1970);
    auto *t_1971 = buffer.data(target + 1971);
    auto *t_1972 = buffer.data(target + 1972);
    auto *t_1973 = buffer.data(target + 1973);
    auto *t_1974 = buffer.data(target + 1974);
    auto *t_1975 = buffer.data(target + 1975);
    auto *t_1976 = buffer.data(target + 1976);
    auto *t_1977 = buffer.data(target + 1977);
    auto *t_1978 = buffer.data(target + 1978);
    auto *t_1979 = buffer.data(target + 1979);
    auto *t_1980 = buffer.data(target + 1980);
    auto *t_1981 = buffer.data(target + 1981);
    auto *t_1982 = buffer.data(target + 1982);
    auto *t_1983 = buffer.data(target + 1983);
    auto *t_1984 = buffer.data(target + 1984);
    auto *t_1985 = buffer.data(target + 1985);
    auto *t_1986 = buffer.data(target + 1986);
    auto *t_1987 = buffer.data(target + 1987);
    auto *t_1988 = buffer.data(target + 1988);
    auto *t_1989 = buffer.data(target + 1989);
    auto *t_1990 = buffer.data(target + 1990);
    auto *t_1991 = buffer.data(target + 1991);
    auto *t_1992 = buffer.data(target + 1992);
    auto *t_1993 = buffer.data(target + 1993);
    auto *t_1994 = buffer.data(target + 1994);
    auto *t_1995 = buffer.data(target + 1995);
    auto *t_1996 = buffer.data(target + 1996);
    auto *t_1997 = buffer.data(target + 1997);
    auto *t_1998 = buffer.data(target + 1998);
    auto *t_1999 = buffer.data(target + 1999);
    auto *t_2000 = buffer.data(target + 2000);
    auto *t_2001 = buffer.data(target + 2001);
    auto *t_2002 = buffer.data(target + 2002);
    auto *t_2003 = buffer.data(target + 2003);
    auto *t_2004 = buffer.data(target + 2004);
    auto *t_2005 = buffer.data(target + 2005);
    auto *t_2006 = buffer.data(target + 2006);
    auto *t_2007 = buffer.data(target + 2007);
    auto *t_2008 = buffer.data(target + 2008);
    auto *t_2009 = buffer.data(target + 2009);
    auto *t_2010 = buffer.data(target + 2010);
    auto *t_2011 = buffer.data(target + 2011);
    auto *t_2012 = buffer.data(target + 2012);
    auto *t_2013 = buffer.data(target + 2013);
    auto *t_2014 = buffer.data(target + 2014);
    auto *t_2015 = buffer.data(target + 2015);
    auto *t_2016 = buffer.data(target + 2016);
    auto *t_2017 = buffer.data(target + 2017);
    auto *t_2018 = buffer.data(target + 2018);
    auto *t_2019 = buffer.data(target + 2019);
    auto *t_2020 = buffer.data(target + 2020);
    auto *t_2021 = buffer.data(target + 2021);
    auto *t_2022 = buffer.data(target + 2022);
    auto *t_2023 = buffer.data(target + 2023);
    auto *t_2024 = buffer.data(target + 2024);
    auto *t_2025 = buffer.data(target + 2025);
    auto *t_2026 = buffer.data(target + 2026);
    auto *t_2027 = buffer.data(target + 2027);
    auto *t_2028 = buffer.data(target + 2028);
    auto *t_2029 = buffer.data(target + 2029);
    auto *t_2030 = buffer.data(target + 2030);
    auto *t_2031 = buffer.data(target + 2031);
    auto *t_2032 = buffer.data(target + 2032);
    auto *t_2033 = buffer.data(target + 2033);
    auto *t_2034 = buffer.data(target + 2034);
    auto *t_2035 = buffer.data(target + 2035);
    auto *t_2036 = buffer.data(target + 2036);
    auto *t_2037 = buffer.data(target + 2037);
    auto *t_2038 = buffer.data(target + 2038);
    auto *t_2039 = buffer.data(target + 2039);
    auto *t_2040 = buffer.data(target + 2040);
    auto *t_2041 = buffer.data(target + 2041);
    auto *t_2042 = buffer.data(target + 2042);
    auto *t_2043 = buffer.data(target + 2043);
    auto *t_2044 = buffer.data(target + 2044);
    auto *t_2045 = buffer.data(target + 2045);
    auto *t_2046 = buffer.data(target + 2046);
    auto *t_2047 = buffer.data(target + 2047);
    auto *t_2048 = buffer.data(target + 2048);

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_1931 = buffer.data(osi0 + 1931);
    const auto *osi0_1932 = buffer.data(osi0 + 1932);
    const auto *osi0_1935 = buffer.data(osi0 + 1935);
    const auto *osi0_1937 = buffer.data(osi0 + 1937);
    const auto *osi0_1938 = buffer.data(osi0 + 1938);
    const auto *osi0_1941 = buffer.data(osi0 + 1941);
    const auto *osi0_1942 = buffer.data(osi0 + 1942);
    const auto *osi0_1944 = buffer.data(osi0 + 1944);
    const auto *osi0_1946 = buffer.data(osi0 + 1946);
    const auto *osi0_1953 = buffer.data(osi0 + 1953);
    const auto *osi0_1955 = buffer.data(osi0 + 1955);
    const auto *osi0_1956 = buffer.data(osi0 + 1956);
    const auto *osi0_1957 = buffer.data(osi0 + 1957);
    const auto *osi0_1959 = buffer.data(osi0 + 1959);
    const auto *osi0_1960 = buffer.data(osi0 + 1960);
    const auto *osi0_1963 = buffer.data(osi0 + 1963);
    const auto *osi0_1965 = buffer.data(osi0 + 1965);
    const auto *osi0_1966 = buffer.data(osi0 + 1966);
    const auto *osi0_1969 = buffer.data(osi0 + 1969);
    const auto *osi0_1970 = buffer.data(osi0 + 1970);
    const auto *osi0_1972 = buffer.data(osi0 + 1972);
    const auto *osi0_1974 = buffer.data(osi0 + 1974);
    const auto *osi0_1981 = buffer.data(osi0 + 1981);
    const auto *osi0_1983 = buffer.data(osi0 + 1983);
    const auto *osi0_1984 = buffer.data(osi0 + 1984);
    const auto *osi0_1985 = buffer.data(osi0 + 1985);
    const auto *osi0_1987 = buffer.data(osi0 + 1987);
    const auto *osi0_1988 = buffer.data(osi0 + 1988);
    const auto *osi0_1991 = buffer.data(osi0 + 1991);
    const auto *osi0_1993 = buffer.data(osi0 + 1993);
    const auto *osi0_1994 = buffer.data(osi0 + 1994);
    const auto *osi0_1997 = buffer.data(osi0 + 1997);
    const auto *osi0_1998 = buffer.data(osi0 + 1998);
    const auto *osi0_2000 = buffer.data(osi0 + 2000);
    const auto *osi0_2002 = buffer.data(osi0 + 2002);
    const auto *osi0_2009 = buffer.data(osi0 + 2009);
    const auto *osi0_2011 = buffer.data(osi0 + 2011);
    const auto *osi0_2012 = buffer.data(osi0 + 2012);
    const auto *osi0_2013 = buffer.data(osi0 + 2013);
    const auto *osi0_2015 = buffer.data(osi0 + 2015);
    const auto *osi0_2016 = buffer.data(osi0 + 2016);
    const auto *osi0_2019 = buffer.data(osi0 + 2019);
    const auto *osi0_2021 = buffer.data(osi0 + 2021);
    const auto *osi0_2022 = buffer.data(osi0 + 2022);
    const auto *osi0_2025 = buffer.data(osi0 + 2025);
    const auto *osi0_2026 = buffer.data(osi0 + 2026);
    const auto *osi0_2028 = buffer.data(osi0 + 2028);
    const auto *osi0_2030 = buffer.data(osi0 + 2030);
    const auto *osi0_2037 = buffer.data(osi0 + 2037);
    const auto *osi0_2039 = buffer.data(osi0 + 2039);
    const auto *osi0_2040 = buffer.data(osi0 + 2040);
    const auto *osi0_2041 = buffer.data(osi0 + 2041);
    const auto *osi0_2043 = buffer.data(osi0 + 2043);
    const auto *osi0_2044 = buffer.data(osi0 + 2044);
    const auto *osi0_2047 = buffer.data(osi0 + 2047);

    const auto *osh_1197 = buffer.data(osh + 1197);
    const auto *osh_1200 = buffer.data(osh + 1200);
    const auto *osh_1203 = buffer.data(osh + 1203);
    const auto *osh_1212 = buffer.data(osh + 1212);
    const auto *osh_1217 = buffer.data(osh + 1217);
    const auto *osh_1218 = buffer.data(osh + 1218);
    const auto *osh_1220 = buffer.data(osh + 1220);
    const auto *osh_1221 = buffer.data(osh + 1221);
    const auto *osh_1223 = buffer.data(osh + 1223);
    const auto *osh_1224 = buffer.data(osh + 1224);
    const auto *osh_1227 = buffer.data(osh + 1227);
    const auto *osh_1233 = buffer.data(osh + 1233);
    const auto *osh_1238 = buffer.data(osh + 1238);
    const auto *osh_1239 = buffer.data(osh + 1239);
    const auto *osh_1241 = buffer.data(osh + 1241);
    const auto *osh_1242 = buffer.data(osh + 1242);
    const auto *osh_1244 = buffer.data(osh + 1244);
    const auto *osh_1245 = buffer.data(osh + 1245);
    const auto *osh_1248 = buffer.data(osh + 1248);
    const auto *osh_1254 = buffer.data(osh + 1254);
    const auto *osh_1259 = buffer.data(osh + 1259);
    const auto *osh_1260 = buffer.data(osh + 1260);
    const auto *osh_1262 = buffer.data(osh + 1262);
    const auto *osh_1263 = buffer.data(osh + 1263);
    const auto *osh_1265 = buffer.data(osh + 1265);
    const auto *osh_1266 = buffer.data(osh + 1266);
    const auto *osh_1269 = buffer.data(osh + 1269);
    const auto *osh_1275 = buffer.data(osh + 1275);
    const auto *osh_1280 = buffer.data(osh + 1280);
    const auto *osh_1281 = buffer.data(osh + 1281);
    const auto *osh_1283 = buffer.data(osh + 1283);
    const auto *osh_1286 = buffer.data(osh + 1286);
    const auto *osh_1290 = buffer.data(osh + 1290);
    const auto *osh_1301 = buffer.data(osh + 1301);
    const auto *osh_1302 = buffer.data(osh + 1302);
    const auto *osh_1304 = buffer.data(osh + 1304);
    const auto *osh_1449 = buffer.data(osh + 1449);
    const auto *osh_1452 = buffer.data(osh + 1452);
    const auto *osh_1454 = buffer.data(osh + 1454);
    const auto *osh_1455 = buffer.data(osh + 1455);
    const auto *osh_1458 = buffer.data(osh + 1458);
    const auto *osh_1459 = buffer.data(osh + 1459);
    const auto *osh_1461 = buffer.data(osh + 1461);
    const auto *osh_1463 = buffer.data(osh + 1463);
    const auto *osh_1464 = buffer.data(osh + 1464);
    const auto *osh_1465 = buffer.data(osh + 1465);
    const auto *osh_1466 = buffer.data(osh + 1466);
    const auto *osh_1467 = buffer.data(osh + 1467);
    const auto *osh_1468 = buffer.data(osh + 1468);
    const auto *osh_1469 = buffer.data(osh + 1469);
    const auto *osh_1470 = buffer.data(osh + 1470);
    const auto *osh_1473 = buffer.data(osh + 1473);
    const auto *osh_1475 = buffer.data(osh + 1475);
    const auto *osh_1476 = buffer.data(osh + 1476);
    const auto *osh_1479 = buffer.data(osh + 1479);
    const auto *osh_1480 = buffer.data(osh + 1480);
    const auto *osh_1482 = buffer.data(osh + 1482);
    const auto *osh_1484 = buffer.data(osh + 1484);
    const auto *osh_1485 = buffer.data(osh + 1485);
    const auto *osh_1486 = buffer.data(osh + 1486);
    const auto *osh_1487 = buffer.data(osh + 1487);
    const auto *osh_1488 = buffer.data(osh + 1488);
    const auto *osh_1489 = buffer.data(osh + 1489);
    const auto *osh_1490 = buffer.data(osh + 1490);
    const auto *osh_1491 = buffer.data(osh + 1491);
    const auto *osh_1494 = buffer.data(osh + 1494);
    const auto *osh_1496 = buffer.data(osh + 1496);
    const auto *osh_1497 = buffer.data(osh + 1497);
    const auto *osh_1500 = buffer.data(osh + 1500);
    const auto *osh_1501 = buffer.data(osh + 1501);
    const auto *osh_1503 = buffer.data(osh + 1503);
    const auto *osh_1505 = buffer.data(osh + 1505);
    const auto *osh_1506 = buffer.data(osh + 1506);
    const auto *osh_1507 = buffer.data(osh + 1507);
    const auto *osh_1508 = buffer.data(osh + 1508);
    const auto *osh_1509 = buffer.data(osh + 1509);
    const auto *osh_1510 = buffer.data(osh + 1510);
    const auto *osh_1511 = buffer.data(osh + 1511);
    const auto *osh_1512 = buffer.data(osh + 1512);
    const auto *osh_1515 = buffer.data(osh + 1515);
    const auto *osh_1517 = buffer.data(osh + 1517);
    const auto *osh_1518 = buffer.data(osh + 1518);
    const auto *osh_1521 = buffer.data(osh + 1521);
    const auto *osh_1522 = buffer.data(osh + 1522);
    const auto *osh_1524 = buffer.data(osh + 1524);
    const auto *osh_1526 = buffer.data(osh + 1526);
    const auto *osh_1527 = buffer.data(osh + 1527);
    const auto *osh_1528 = buffer.data(osh + 1528);
    const auto *osh_1529 = buffer.data(osh + 1529);
    const auto *osh_1530 = buffer.data(osh + 1530);
    const auto *osh_1531 = buffer.data(osh + 1531);
    const auto *osh_1532 = buffer.data(osh + 1532);
    const auto *osh_1533 = buffer.data(osh + 1533);
    const auto *osh_1536 = buffer.data(osh + 1536);

    const auto *osi1_1931 = buffer.data(osi1 + 1931);
    const auto *osi1_1932 = buffer.data(osi1 + 1932);
    const auto *osi1_1935 = buffer.data(osi1 + 1935);
    const auto *osi1_1937 = buffer.data(osi1 + 1937);
    const auto *osi1_1938 = buffer.data(osi1 + 1938);
    const auto *osi1_1941 = buffer.data(osi1 + 1941);
    const auto *osi1_1942 = buffer.data(osi1 + 1942);
    const auto *osi1_1944 = buffer.data(osi1 + 1944);
    const auto *osi1_1946 = buffer.data(osi1 + 1946);
    const auto *osi1_1953 = buffer.data(osi1 + 1953);
    const auto *osi1_1955 = buffer.data(osi1 + 1955);
    const auto *osi1_1956 = buffer.data(osi1 + 1956);
    const auto *osi1_1957 = buffer.data(osi1 + 1957);
    const auto *osi1_1959 = buffer.data(osi1 + 1959);
    const auto *osi1_1960 = buffer.data(osi1 + 1960);
    const auto *osi1_1963 = buffer.data(osi1 + 1963);
    const auto *osi1_1965 = buffer.data(osi1 + 1965);
    const auto *osi1_1966 = buffer.data(osi1 + 1966);
    const auto *osi1_1969 = buffer.data(osi1 + 1969);
    const auto *osi1_1970 = buffer.data(osi1 + 1970);
    const auto *osi1_1972 = buffer.data(osi1 + 1972);
    const auto *osi1_1974 = buffer.data(osi1 + 1974);
    const auto *osi1_1981 = buffer.data(osi1 + 1981);
    const auto *osi1_1983 = buffer.data(osi1 + 1983);
    const auto *osi1_1984 = buffer.data(osi1 + 1984);
    const auto *osi1_1985 = buffer.data(osi1 + 1985);
    const auto *osi1_1987 = buffer.data(osi1 + 1987);
    const auto *osi1_1988 = buffer.data(osi1 + 1988);
    const auto *osi1_1991 = buffer.data(osi1 + 1991);
    const auto *osi1_1993 = buffer.data(osi1 + 1993);
    const auto *osi1_1994 = buffer.data(osi1 + 1994);
    const auto *osi1_1997 = buffer.data(osi1 + 1997);
    const auto *osi1_1998 = buffer.data(osi1 + 1998);
    const auto *osi1_2000 = buffer.data(osi1 + 2000);
    const auto *osi1_2002 = buffer.data(osi1 + 2002);
    const auto *osi1_2009 = buffer.data(osi1 + 2009);
    const auto *osi1_2011 = buffer.data(osi1 + 2011);
    const auto *osi1_2012 = buffer.data(osi1 + 2012);
    const auto *osi1_2013 = buffer.data(osi1 + 2013);
    const auto *osi1_2015 = buffer.data(osi1 + 2015);
    const auto *osi1_2016 = buffer.data(osi1 + 2016);
    const auto *osi1_2019 = buffer.data(osi1 + 2019);
    const auto *osi1_2021 = buffer.data(osi1 + 2021);
    const auto *osi1_2022 = buffer.data(osi1 + 2022);
    const auto *osi1_2025 = buffer.data(osi1 + 2025);
    const auto *osi1_2026 = buffer.data(osi1 + 2026);
    const auto *osi1_2028 = buffer.data(osi1 + 2028);
    const auto *osi1_2030 = buffer.data(osi1 + 2030);
    const auto *osi1_2037 = buffer.data(osi1 + 2037);
    const auto *osi1_2039 = buffer.data(osi1 + 2039);
    const auto *osi1_2040 = buffer.data(osi1 + 2040);
    const auto *osi1_2041 = buffer.data(osi1 + 2041);
    const auto *osi1_2043 = buffer.data(osi1 + 2043);
    const auto *osi1_2044 = buffer.data(osi1 + 2044);
    const auto *osi1_2047 = buffer.data(osi1 + 2047);

    const auto *qsh_1448 = buffer.data(qsh + 1448);
    const auto *qsh_1449 = buffer.data(qsh + 1449);
    const auto *qsh_1451 = buffer.data(qsh + 1451);
    const auto *qsh_1452 = buffer.data(qsh + 1452);
    const auto *qsh_1454 = buffer.data(qsh + 1454);
    const auto *qsh_1455 = buffer.data(qsh + 1455);
    const auto *qsh_1458 = buffer.data(qsh + 1458);
    const auto *qsh_1464 = buffer.data(qsh + 1464);
    const auto *qsh_1465 = buffer.data(qsh + 1465);
    const auto *qsh_1466 = buffer.data(qsh + 1466);
    const auto *qsh_1467 = buffer.data(qsh + 1467);
    const auto *qsh_1468 = buffer.data(qsh + 1468);
    const auto *qsh_1469 = buffer.data(qsh + 1469);
    const auto *qsh_1470 = buffer.data(qsh + 1470);
    const auto *qsh_1472 = buffer.data(qsh + 1472);
    const auto *qsh_1473 = buffer.data(qsh + 1473);
    const auto *qsh_1475 = buffer.data(qsh + 1475);
    const auto *qsh_1476 = buffer.data(qsh + 1476);
    const auto *qsh_1479 = buffer.data(qsh + 1479);
    const auto *qsh_1485 = buffer.data(qsh + 1485);
    const auto *qsh_1486 = buffer.data(qsh + 1486);
    const auto *qsh_1487 = buffer.data(qsh + 1487);
    const auto *qsh_1488 = buffer.data(qsh + 1488);
    const auto *qsh_1489 = buffer.data(qsh + 1489);
    const auto *qsh_1490 = buffer.data(qsh + 1490);
    const auto *qsh_1491 = buffer.data(qsh + 1491);
    const auto *qsh_1493 = buffer.data(qsh + 1493);
    const auto *qsh_1494 = buffer.data(qsh + 1494);
    const auto *qsh_1496 = buffer.data(qsh + 1496);
    const auto *qsh_1497 = buffer.data(qsh + 1497);
    const auto *qsh_1500 = buffer.data(qsh + 1500);
    const auto *qsh_1506 = buffer.data(qsh + 1506);
    const auto *qsh_1507 = buffer.data(qsh + 1507);
    const auto *qsh_1508 = buffer.data(qsh + 1508);
    const auto *qsh_1509 = buffer.data(qsh + 1509);
    const auto *qsh_1510 = buffer.data(qsh + 1510);
    const auto *qsh_1511 = buffer.data(qsh + 1511);
    const auto *qsh_1512 = buffer.data(qsh + 1512);
    const auto *qsh_1514 = buffer.data(qsh + 1514);
    const auto *qsh_1515 = buffer.data(qsh + 1515);
    const auto *qsh_1517 = buffer.data(qsh + 1517);
    const auto *qsh_1518 = buffer.data(qsh + 1518);
    const auto *qsh_1521 = buffer.data(qsh + 1521);
    const auto *qsh_1527 = buffer.data(qsh + 1527);
    const auto *qsh_1528 = buffer.data(qsh + 1528);
    const auto *qsh_1529 = buffer.data(qsh + 1529);
    const auto *qsh_1530 = buffer.data(qsh + 1530);
    const auto *qsh_1531 = buffer.data(qsh + 1531);
    const auto *qsh_1532 = buffer.data(qsh + 1532);
    const auto *qsh_1533 = buffer.data(qsh + 1533);
    const auto *qsh_1535 = buffer.data(qsh + 1535);

#pragma omp simd aligned(t_1930, t_1931, t_1932, t_1933, pa_x, pc_x, pc_y, osi0_1931, \
                         osi0_1932, osh_1217, osh_1218, osh_1449, osi1_1931, osi1_1932, \
                         qsh_1448, qsh_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1930[k] = f_19 * osh_1217[k]
                    + f_3 * pc_y[k] * qsh_1448[k];

        t_1931[k] = pa_x[k] * osi0_1931[k]
                    - f_10 * pc_x[k] * osi1_1931[k];

        t_1932[k] = pa_x[k] * osi0_1932[k]
                    + f_23 * osh_1449[k]
                    - f_10 * pc_x[k] * osi1_1932[k];

        t_1933[k] = f_20 * osh_1218[k]
                    + f_3 * pc_y[k] * qsh_1449[k];
    }

#pragma omp simd aligned(t_1934, t_1935, t_1936, pa_x, pc_x, pc_y, pc_z, osi0_1935, osh_1197, \
                         osh_1220, osh_1452, osi1_1935, qsh_1449, \
                         qsh_1451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1934[k] = f_13 * osh_1197[k]
                    + f_3 * pc_z[k] * qsh_1449[k];

        t_1935[k] = pa_x[k] * osi0_1935[k]
                    + f_14 * osh_1452[k]
                    - f_10 * pc_x[k] * osi1_1935[k];

        t_1936[k] = f_20 * osh_1220[k]
                    + f_3 * pc_y[k] * qsh_1451[k];
    }

#pragma omp simd aligned(t_1937, t_1938, t_1939, pa_x, pc_x, pc_z, osi0_1937, osi0_1938, \
                         osh_1200, osh_1454, osh_1455, osi1_1937, osi1_1938, \
                         qsh_1452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1937[k] = pa_x[k] * osi0_1937[k]
                    + f_14 * osh_1454[k]
                    - f_10 * pc_x[k] * osi1_1937[k];

        t_1938[k] = pa_x[k] * osi0_1938[k]
                    + f_13 * osh_1455[k]
                    - f_10 * pc_x[k] * osi1_1938[k];

        t_1939[k] = f_13 * osh_1200[k]
                    + f_3 * pc_z[k] * qsh_1452[k];
    }

#pragma omp simd aligned(t_1940, t_1941, t_1942, pa_x, pc_x, pc_y, osi0_1941, osi0_1942, \
                         osh_1223, osh_1458, osh_1459, osi1_1941, osi1_1942, \
                         qsh_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1940[k] = f_20 * osh_1223[k]
                    + f_3 * pc_y[k] * qsh_1454[k];

        t_1941[k] = pa_x[k] * osi0_1941[k]
                    + f_13 * osh_1458[k]
                    - f_10 * pc_x[k] * osi1_1941[k];

        t_1942[k] = pa_x[k] * osi0_1942[k]
                    + f_12 * osh_1459[k]
                    - f_10 * pc_x[k] * osi1_1942[k];
    }

#pragma omp simd aligned(t_1943, t_1944, t_1945, pa_x, pc_x, pc_y, pc_z, osi0_1944, osh_1203, \
                         osh_1227, osh_1461, osi1_1944, qsh_1455, \
                         qsh_1458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1943[k] = f_13 * osh_1203[k]
                    + f_3 * pc_z[k] * qsh_1455[k];

        t_1944[k] = pa_x[k] * osi0_1944[k]
                    + f_12 * osh_1461[k]
                    - f_10 * pc_x[k] * osi1_1944[k];

        t_1945[k] = f_20 * osh_1227[k]
                    + f_3 * pc_y[k] * qsh_1458[k];
    }

#pragma omp simd aligned(t_1946, t_1947, t_1948, t_1949, pa_x, pc_x, osi0_1946, osh_1463, \
                         osh_1464, osh_1465, osh_1466, osi1_1946, qsh_1464, qsh_1465, \
                         qsh_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1946[k] = pa_x[k] * osi0_1946[k]
                    + f_12 * osh_1463[k]
                    - f_10 * pc_x[k] * osi1_1946[k];

        t_1947[k] = f_11 * osh_1464[k]
                    + f_3 * pc_x[k] * qsh_1464[k];

        t_1948[k] = f_11 * osh_1465[k]
                    + f_3 * pc_x[k] * qsh_1465[k];

        t_1949[k] = f_11 * osh_1466[k]
                    + f_3 * pc_x[k] * qsh_1466[k];
    }

#pragma omp simd aligned(t_1950, t_1951, t_1952, t_1953, pa_x, pc_x, osi0_1953, osh_1467, \
                         osh_1468, osh_1469, osi1_1953, qsh_1467, qsh_1468, \
                         qsh_1469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1950[k] = f_11 * osh_1467[k]
                    + f_3 * pc_x[k] * qsh_1467[k];

        t_1951[k] = f_11 * osh_1468[k]
                    + f_3 * pc_x[k] * qsh_1468[k];

        t_1952[k] = f_11 * osh_1469[k]
                    + f_3 * pc_x[k] * qsh_1469[k];

        t_1953[k] = pa_x[k] * osi0_1953[k]
                    - f_10 * pc_x[k] * osi1_1953[k];
    }

#pragma omp simd aligned(t_1954, t_1955, t_1956, t_1957, pa_x, pc_x, pc_z, osi0_1955, \
                         osi0_1956, osi0_1957, osh_1212, osi1_1955, osi1_1956, osi1_1957, \
                         qsh_1464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1954[k] = f_13 * osh_1212[k]
                    + f_3 * pc_z[k] * qsh_1464[k];

        t_1955[k] = pa_x[k] * osi0_1955[k]
                    - f_10 * pc_x[k] * osi1_1955[k];

        t_1956[k] = pa_x[k] * osi0_1956[k]
                    - f_10 * pc_x[k] * osi1_1956[k];

        t_1957[k] = pa_x[k] * osi0_1957[k]
                    - f_10 * pc_x[k] * osi1_1957[k];
    }

#pragma omp simd aligned(t_1958, t_1959, t_1960, t_1961, pa_x, pc_x, pc_y, osi0_1959, \
                         osi0_1960, osh_1238, osh_1239, osh_1470, osi1_1959, osi1_1960, \
                         qsh_1469, qsh_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1958[k] = f_20 * osh_1238[k]
                    + f_3 * pc_y[k] * qsh_1469[k];

        t_1959[k] = pa_x[k] * osi0_1959[k]
                    - f_10 * pc_x[k] * osi1_1959[k];

        t_1960[k] = pa_x[k] * osi0_1960[k]
                    + f_23 * osh_1470[k]
                    - f_10 * pc_x[k] * osi1_1960[k];

        t_1961[k] = f_21 * osh_1239[k]
                    + f_3 * pc_y[k] * qsh_1470[k];
    }

#pragma omp simd aligned(t_1962, t_1963, t_1964, pa_x, pc_x, pc_y, pc_z, osi0_1963, osh_1218, \
                         osh_1241, osh_1473, osi1_1963, qsh_1470, \
                         qsh_1472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1962[k] = f_14 * osh_1218[k]
                    + f_3 * pc_z[k] * qsh_1470[k];

        t_1963[k] = pa_x[k] * osi0_1963[k]
                    + f_14 * osh_1473[k]
                    - f_10 * pc_x[k] * osi1_1963[k];

        t_1964[k] = f_21 * osh_1241[k]
                    + f_3 * pc_y[k] * qsh_1472[k];
    }

#pragma omp simd aligned(t_1965, t_1966, t_1967, pa_x, pc_x, pc_z, osi0_1965, osi0_1966, \
                         osh_1221, osh_1475, osh_1476, osi1_1965, osi1_1966, \
                         qsh_1473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1965[k] = pa_x[k] * osi0_1965[k]
                    + f_14 * osh_1475[k]
                    - f_10 * pc_x[k] * osi1_1965[k];

        t_1966[k] = pa_x[k] * osi0_1966[k]
                    + f_13 * osh_1476[k]
                    - f_10 * pc_x[k] * osi1_1966[k];

        t_1967[k] = f_14 * osh_1221[k]
                    + f_3 * pc_z[k] * qsh_1473[k];
    }

#pragma omp simd aligned(t_1968, t_1969, t_1970, pa_x, pc_x, pc_y, osi0_1969, osi0_1970, \
                         osh_1244, osh_1479, osh_1480, osi1_1969, osi1_1970, \
                         qsh_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1968[k] = f_21 * osh_1244[k]
                    + f_3 * pc_y[k] * qsh_1475[k];

        t_1969[k] = pa_x[k] * osi0_1969[k]
                    + f_13 * osh_1479[k]
                    - f_10 * pc_x[k] * osi1_1969[k];

        t_1970[k] = pa_x[k] * osi0_1970[k]
                    + f_12 * osh_1480[k]
                    - f_10 * pc_x[k] * osi1_1970[k];
    }

#pragma omp simd aligned(t_1971, t_1972, t_1973, pa_x, pc_x, pc_y, pc_z, osi0_1972, osh_1224, \
                         osh_1248, osh_1482, osi1_1972, qsh_1476, \
                         qsh_1479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1971[k] = f_14 * osh_1224[k]
                    + f_3 * pc_z[k] * qsh_1476[k];

        t_1972[k] = pa_x[k] * osi0_1972[k]
                    + f_12 * osh_1482[k]
                    - f_10 * pc_x[k] * osi1_1972[k];

        t_1973[k] = f_21 * osh_1248[k]
                    + f_3 * pc_y[k] * qsh_1479[k];
    }

#pragma omp simd aligned(t_1974, t_1975, t_1976, t_1977, pa_x, pc_x, osi0_1974, osh_1484, \
                         osh_1485, osh_1486, osh_1487, osi1_1974, qsh_1485, qsh_1486, \
                         qsh_1487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1974[k] = pa_x[k] * osi0_1974[k]
                    + f_12 * osh_1484[k]
                    - f_10 * pc_x[k] * osi1_1974[k];

        t_1975[k] = f_11 * osh_1485[k]
                    + f_3 * pc_x[k] * qsh_1485[k];

        t_1976[k] = f_11 * osh_1486[k]
                    + f_3 * pc_x[k] * qsh_1486[k];

        t_1977[k] = f_11 * osh_1487[k]
                    + f_3 * pc_x[k] * qsh_1487[k];
    }

#pragma omp simd aligned(t_1978, t_1979, t_1980, t_1981, pa_x, pc_x, osi0_1981, osh_1488, \
                         osh_1489, osh_1490, osi1_1981, qsh_1488, qsh_1489, \
                         qsh_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1978[k] = f_11 * osh_1488[k]
                    + f_3 * pc_x[k] * qsh_1488[k];

        t_1979[k] = f_11 * osh_1489[k]
                    + f_3 * pc_x[k] * qsh_1489[k];

        t_1980[k] = f_11 * osh_1490[k]
                    + f_3 * pc_x[k] * qsh_1490[k];

        t_1981[k] = pa_x[k] * osi0_1981[k]
                    - f_10 * pc_x[k] * osi1_1981[k];
    }

#pragma omp simd aligned(t_1982, t_1983, t_1984, t_1985, pa_x, pc_x, pc_z, osi0_1983, \
                         osi0_1984, osi0_1985, osh_1233, osi1_1983, osi1_1984, osi1_1985, \
                         qsh_1485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1982[k] = f_14 * osh_1233[k]
                    + f_3 * pc_z[k] * qsh_1485[k];

        t_1983[k] = pa_x[k] * osi0_1983[k]
                    - f_10 * pc_x[k] * osi1_1983[k];

        t_1984[k] = pa_x[k] * osi0_1984[k]
                    - f_10 * pc_x[k] * osi1_1984[k];

        t_1985[k] = pa_x[k] * osi0_1985[k]
                    - f_10 * pc_x[k] * osi1_1985[k];
    }

#pragma omp simd aligned(t_1986, t_1987, t_1988, t_1989, pa_x, pc_x, pc_y, osi0_1987, \
                         osi0_1988, osh_1259, osh_1260, osh_1491, osi1_1987, osi1_1988, \
                         qsh_1490, qsh_1491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1986[k] = f_21 * osh_1259[k]
                    + f_3 * pc_y[k] * qsh_1490[k];

        t_1987[k] = pa_x[k] * osi0_1987[k]
                    - f_10 * pc_x[k] * osi1_1987[k];

        t_1988[k] = pa_x[k] * osi0_1988[k]
                    + f_23 * osh_1491[k]
                    - f_10 * pc_x[k] * osi1_1988[k];

        t_1989[k] = f_23 * osh_1260[k]
                    + f_3 * pc_y[k] * qsh_1491[k];
    }

#pragma omp simd aligned(t_1990, t_1991, t_1992, pa_x, pc_x, pc_y, pc_z, osi0_1991, osh_1239, \
                         osh_1262, osh_1494, osi1_1991, qsh_1491, \
                         qsh_1493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1990[k] = f_22 * osh_1239[k]
                    + f_3 * pc_z[k] * qsh_1491[k];

        t_1991[k] = pa_x[k] * osi0_1991[k]
                    + f_14 * osh_1494[k]
                    - f_10 * pc_x[k] * osi1_1991[k];

        t_1992[k] = f_23 * osh_1262[k]
                    + f_3 * pc_y[k] * qsh_1493[k];
    }

#pragma omp simd aligned(t_1993, t_1994, t_1995, pa_x, pc_x, pc_z, osi0_1993, osi0_1994, \
                         osh_1242, osh_1496, osh_1497, osi1_1993, osi1_1994, \
                         qsh_1494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1993[k] = pa_x[k] * osi0_1993[k]
                    + f_14 * osh_1496[k]
                    - f_10 * pc_x[k] * osi1_1993[k];

        t_1994[k] = pa_x[k] * osi0_1994[k]
                    + f_13 * osh_1497[k]
                    - f_10 * pc_x[k] * osi1_1994[k];

        t_1995[k] = f_22 * osh_1242[k]
                    + f_3 * pc_z[k] * qsh_1494[k];
    }

#pragma omp simd aligned(t_1996, t_1997, t_1998, pa_x, pc_x, pc_y, osi0_1997, osi0_1998, \
                         osh_1265, osh_1500, osh_1501, osi1_1997, osi1_1998, \
                         qsh_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1996[k] = f_23 * osh_1265[k]
                    + f_3 * pc_y[k] * qsh_1496[k];

        t_1997[k] = pa_x[k] * osi0_1997[k]
                    + f_13 * osh_1500[k]
                    - f_10 * pc_x[k] * osi1_1997[k];

        t_1998[k] = pa_x[k] * osi0_1998[k]
                    + f_12 * osh_1501[k]
                    - f_10 * pc_x[k] * osi1_1998[k];
    }

#pragma omp simd aligned(t_1999, t_2000, t_2001, pa_x, pc_x, pc_y, pc_z, osi0_2000, osh_1245, \
                         osh_1269, osh_1503, osi1_2000, qsh_1497, \
                         qsh_1500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1999[k] = f_22 * osh_1245[k]
                    + f_3 * pc_z[k] * qsh_1497[k];

        t_2000[k] = pa_x[k] * osi0_2000[k]
                    + f_12 * osh_1503[k]
                    - f_10 * pc_x[k] * osi1_2000[k];

        t_2001[k] = f_23 * osh_1269[k]
                    + f_3 * pc_y[k] * qsh_1500[k];
    }

#pragma omp simd aligned(t_2002, t_2003, t_2004, t_2005, pa_x, pc_x, osi0_2002, osh_1505, \
                         osh_1506, osh_1507, osh_1508, osi1_2002, qsh_1506, qsh_1507, \
                         qsh_1508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2002[k] = pa_x[k] * osi0_2002[k]
                    + f_12 * osh_1505[k]
                    - f_10 * pc_x[k] * osi1_2002[k];

        t_2003[k] = f_11 * osh_1506[k]
                    + f_3 * pc_x[k] * qsh_1506[k];

        t_2004[k] = f_11 * osh_1507[k]
                    + f_3 * pc_x[k] * qsh_1507[k];

        t_2005[k] = f_11 * osh_1508[k]
                    + f_3 * pc_x[k] * qsh_1508[k];
    }

#pragma omp simd aligned(t_2006, t_2007, t_2008, t_2009, pa_x, pc_x, osi0_2009, osh_1509, \
                         osh_1510, osh_1511, osi1_2009, qsh_1509, qsh_1510, \
                         qsh_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2006[k] = f_11 * osh_1509[k]
                    + f_3 * pc_x[k] * qsh_1509[k];

        t_2007[k] = f_11 * osh_1510[k]
                    + f_3 * pc_x[k] * qsh_1510[k];

        t_2008[k] = f_11 * osh_1511[k]
                    + f_3 * pc_x[k] * qsh_1511[k];

        t_2009[k] = pa_x[k] * osi0_2009[k]
                    - f_10 * pc_x[k] * osi1_2009[k];
    }

#pragma omp simd aligned(t_2010, t_2011, t_2012, t_2013, pa_x, pc_x, pc_z, osi0_2011, \
                         osi0_2012, osi0_2013, osh_1254, osi1_2011, osi1_2012, osi1_2013, \
                         qsh_1506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2010[k] = f_22 * osh_1254[k]
                    + f_3 * pc_z[k] * qsh_1506[k];

        t_2011[k] = pa_x[k] * osi0_2011[k]
                    - f_10 * pc_x[k] * osi1_2011[k];

        t_2012[k] = pa_x[k] * osi0_2012[k]
                    - f_10 * pc_x[k] * osi1_2012[k];

        t_2013[k] = pa_x[k] * osi0_2013[k]
                    - f_10 * pc_x[k] * osi1_2013[k];
    }

#pragma omp simd aligned(t_2014, t_2015, t_2016, t_2017, pa_x, pc_x, pc_y, osi0_2015, \
                         osi0_2016, osh_1280, osh_1281, osh_1512, osi1_2015, osi1_2016, \
                         qsh_1511, qsh_1512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2014[k] = f_23 * osh_1280[k]
                    + f_3 * pc_y[k] * qsh_1511[k];

        t_2015[k] = pa_x[k] * osi0_2015[k]
                    - f_10 * pc_x[k] * osi1_2015[k];

        t_2016[k] = pa_x[k] * osi0_2016[k]
                    + f_23 * osh_1512[k]
                    - f_10 * pc_x[k] * osi1_2016[k];

        t_2017[k] = f_22 * osh_1281[k]
                    + f_3 * pc_y[k] * qsh_1512[k];
    }

#pragma omp simd aligned(t_2018, t_2019, t_2020, pa_x, pc_x, pc_y, pc_z, osi0_2019, osh_1260, \
                         osh_1283, osh_1515, osi1_2019, qsh_1512, \
                         qsh_1514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2018[k] = f_23 * osh_1260[k]
                    + f_3 * pc_z[k] * qsh_1512[k];

        t_2019[k] = pa_x[k] * osi0_2019[k]
                    + f_14 * osh_1515[k]
                    - f_10 * pc_x[k] * osi1_2019[k];

        t_2020[k] = f_22 * osh_1283[k]
                    + f_3 * pc_y[k] * qsh_1514[k];
    }

#pragma omp simd aligned(t_2021, t_2022, t_2023, pa_x, pc_x, pc_z, osi0_2021, osi0_2022, \
                         osh_1263, osh_1517, osh_1518, osi1_2021, osi1_2022, \
                         qsh_1515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2021[k] = pa_x[k] * osi0_2021[k]
                    + f_14 * osh_1517[k]
                    - f_10 * pc_x[k] * osi1_2021[k];

        t_2022[k] = pa_x[k] * osi0_2022[k]
                    + f_13 * osh_1518[k]
                    - f_10 * pc_x[k] * osi1_2022[k];

        t_2023[k] = f_23 * osh_1263[k]
                    + f_3 * pc_z[k] * qsh_1515[k];
    }

#pragma omp simd aligned(t_2024, t_2025, t_2026, pa_x, pc_x, pc_y, osi0_2025, osi0_2026, \
                         osh_1286, osh_1521, osh_1522, osi1_2025, osi1_2026, \
                         qsh_1517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2024[k] = f_22 * osh_1286[k]
                    + f_3 * pc_y[k] * qsh_1517[k];

        t_2025[k] = pa_x[k] * osi0_2025[k]
                    + f_13 * osh_1521[k]
                    - f_10 * pc_x[k] * osi1_2025[k];

        t_2026[k] = pa_x[k] * osi0_2026[k]
                    + f_12 * osh_1522[k]
                    - f_10 * pc_x[k] * osi1_2026[k];
    }

#pragma omp simd aligned(t_2027, t_2028, t_2029, pa_x, pc_x, pc_y, pc_z, osi0_2028, osh_1266, \
                         osh_1290, osh_1524, osi1_2028, qsh_1518, \
                         qsh_1521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2027[k] = f_23 * osh_1266[k]
                    + f_3 * pc_z[k] * qsh_1518[k];

        t_2028[k] = pa_x[k] * osi0_2028[k]
                    + f_12 * osh_1524[k]
                    - f_10 * pc_x[k] * osi1_2028[k];

        t_2029[k] = f_22 * osh_1290[k]
                    + f_3 * pc_y[k] * qsh_1521[k];
    }

#pragma omp simd aligned(t_2030, t_2031, t_2032, t_2033, pa_x, pc_x, osi0_2030, osh_1526, \
                         osh_1527, osh_1528, osh_1529, osi1_2030, qsh_1527, qsh_1528, \
                         qsh_1529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2030[k] = pa_x[k] * osi0_2030[k]
                    + f_12 * osh_1526[k]
                    - f_10 * pc_x[k] * osi1_2030[k];

        t_2031[k] = f_11 * osh_1527[k]
                    + f_3 * pc_x[k] * qsh_1527[k];

        t_2032[k] = f_11 * osh_1528[k]
                    + f_3 * pc_x[k] * qsh_1528[k];

        t_2033[k] = f_11 * osh_1529[k]
                    + f_3 * pc_x[k] * qsh_1529[k];
    }

#pragma omp simd aligned(t_2034, t_2035, t_2036, t_2037, pa_x, pc_x, osi0_2037, osh_1530, \
                         osh_1531, osh_1532, osi1_2037, qsh_1530, qsh_1531, \
                         qsh_1532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2034[k] = f_11 * osh_1530[k]
                    + f_3 * pc_x[k] * qsh_1530[k];

        t_2035[k] = f_11 * osh_1531[k]
                    + f_3 * pc_x[k] * qsh_1531[k];

        t_2036[k] = f_11 * osh_1532[k]
                    + f_3 * pc_x[k] * qsh_1532[k];

        t_2037[k] = pa_x[k] * osi0_2037[k]
                    - f_10 * pc_x[k] * osi1_2037[k];
    }

#pragma omp simd aligned(t_2038, t_2039, t_2040, t_2041, pa_x, pc_x, pc_z, osi0_2039, \
                         osi0_2040, osi0_2041, osh_1275, osi1_2039, osi1_2040, osi1_2041, \
                         qsh_1527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2038[k] = f_23 * osh_1275[k]
                    + f_3 * pc_z[k] * qsh_1527[k];

        t_2039[k] = pa_x[k] * osi0_2039[k]
                    - f_10 * pc_x[k] * osi1_2039[k];

        t_2040[k] = pa_x[k] * osi0_2040[k]
                    - f_10 * pc_x[k] * osi1_2040[k];

        t_2041[k] = pa_x[k] * osi0_2041[k]
                    - f_10 * pc_x[k] * osi1_2041[k];
    }

#pragma omp simd aligned(t_2042, t_2043, t_2044, t_2045, pa_x, pc_x, pc_y, osi0_2043, \
                         osi0_2044, osh_1301, osh_1302, osh_1533, osi1_2043, osi1_2044, \
                         qsh_1532, qsh_1533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2042[k] = f_22 * osh_1301[k]
                    + f_3 * pc_y[k] * qsh_1532[k];

        t_2043[k] = pa_x[k] * osi0_2043[k]
                    - f_10 * pc_x[k] * osi1_2043[k];

        t_2044[k] = pa_x[k] * osi0_2044[k]
                    + f_23 * osh_1533[k]
                    - f_10 * pc_x[k] * osi1_2044[k];

        t_2045[k] = f_14 * osh_1302[k]
                    + f_3 * pc_y[k] * qsh_1533[k];
    }

#pragma omp simd aligned(t_2046, t_2047, t_2048, pa_x, pc_x, pc_y, pc_z, osi0_2047, osh_1281, \
                         osh_1304, osh_1536, osi1_2047, qsh_1533, \
                         qsh_1535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2046[k] = f_21 * osh_1281[k]
                    + f_3 * pc_z[k] * qsh_1533[k];

        t_2047[k] = pa_x[k] * osi0_2047[k]
                    + f_14 * osh_1536[k]
                    - f_10 * pc_x[k] * osi1_2047[k];

        t_2048[k] = f_14 * osh_1304[k]
                    + f_3 * pc_y[k] * qsh_1535[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece18(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsg0, const size_t qsg1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.5 / q;
    const auto f_18 = 5.0 / q;
    const auto f_19 = 4.5 / q;
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_2049 = buffer.data(target + 2049);
    auto *t_2050 = buffer.data(target + 2050);
    auto *t_2051 = buffer.data(target + 2051);
    auto *t_2052 = buffer.data(target + 2052);
    auto *t_2053 = buffer.data(target + 2053);
    auto *t_2054 = buffer.data(target + 2054);
    auto *t_2055 = buffer.data(target + 2055);
    auto *t_2056 = buffer.data(target + 2056);
    auto *t_2057 = buffer.data(target + 2057);
    auto *t_2058 = buffer.data(target + 2058);
    auto *t_2059 = buffer.data(target + 2059);
    auto *t_2060 = buffer.data(target + 2060);
    auto *t_2061 = buffer.data(target + 2061);
    auto *t_2062 = buffer.data(target + 2062);
    auto *t_2063 = buffer.data(target + 2063);
    auto *t_2064 = buffer.data(target + 2064);
    auto *t_2065 = buffer.data(target + 2065);
    auto *t_2066 = buffer.data(target + 2066);
    auto *t_2067 = buffer.data(target + 2067);
    auto *t_2068 = buffer.data(target + 2068);
    auto *t_2069 = buffer.data(target + 2069);
    auto *t_2070 = buffer.data(target + 2070);
    auto *t_2071 = buffer.data(target + 2071);
    auto *t_2072 = buffer.data(target + 2072);
    auto *t_2073 = buffer.data(target + 2073);
    auto *t_2074 = buffer.data(target + 2074);
    auto *t_2075 = buffer.data(target + 2075);
    auto *t_2076 = buffer.data(target + 2076);
    auto *t_2077 = buffer.data(target + 2077);
    auto *t_2078 = buffer.data(target + 2078);
    auto *t_2079 = buffer.data(target + 2079);
    auto *t_2080 = buffer.data(target + 2080);
    auto *t_2081 = buffer.data(target + 2081);
    auto *t_2082 = buffer.data(target + 2082);
    auto *t_2083 = buffer.data(target + 2083);
    auto *t_2084 = buffer.data(target + 2084);
    auto *t_2085 = buffer.data(target + 2085);
    auto *t_2086 = buffer.data(target + 2086);
    auto *t_2087 = buffer.data(target + 2087);
    auto *t_2088 = buffer.data(target + 2088);
    auto *t_2089 = buffer.data(target + 2089);
    auto *t_2090 = buffer.data(target + 2090);
    auto *t_2091 = buffer.data(target + 2091);
    auto *t_2092 = buffer.data(target + 2092);
    auto *t_2093 = buffer.data(target + 2093);
    auto *t_2094 = buffer.data(target + 2094);
    auto *t_2095 = buffer.data(target + 2095);
    auto *t_2096 = buffer.data(target + 2096);
    auto *t_2097 = buffer.data(target + 2097);
    auto *t_2098 = buffer.data(target + 2098);
    auto *t_2099 = buffer.data(target + 2099);
    auto *t_2100 = buffer.data(target + 2100);
    auto *t_2101 = buffer.data(target + 2101);
    auto *t_2102 = buffer.data(target + 2102);
    auto *t_2103 = buffer.data(target + 2103);
    auto *t_2104 = buffer.data(target + 2104);
    auto *t_2105 = buffer.data(target + 2105);
    auto *t_2106 = buffer.data(target + 2106);
    auto *t_2107 = buffer.data(target + 2107);
    auto *t_2108 = buffer.data(target + 2108);
    auto *t_2109 = buffer.data(target + 2109);
    auto *t_2110 = buffer.data(target + 2110);
    auto *t_2111 = buffer.data(target + 2111);
    auto *t_2112 = buffer.data(target + 2112);
    auto *t_2113 = buffer.data(target + 2113);
    auto *t_2114 = buffer.data(target + 2114);
    auto *t_2115 = buffer.data(target + 2115);
    auto *t_2116 = buffer.data(target + 2116);
    auto *t_2117 = buffer.data(target + 2117);
    auto *t_2118 = buffer.data(target + 2118);
    auto *t_2119 = buffer.data(target + 2119);
    auto *t_2120 = buffer.data(target + 2120);
    auto *t_2121 = buffer.data(target + 2121);
    auto *t_2122 = buffer.data(target + 2122);
    auto *t_2123 = buffer.data(target + 2123);
    auto *t_2124 = buffer.data(target + 2124);
    auto *t_2125 = buffer.data(target + 2125);
    auto *t_2126 = buffer.data(target + 2126);
    auto *t_2127 = buffer.data(target + 2127);
    auto *t_2128 = buffer.data(target + 2128);
    auto *t_2129 = buffer.data(target + 2129);
    auto *t_2130 = buffer.data(target + 2130);
    auto *t_2131 = buffer.data(target + 2131);
    auto *t_2132 = buffer.data(target + 2132);
    auto *t_2133 = buffer.data(target + 2133);
    auto *t_2134 = buffer.data(target + 2134);
    auto *t_2135 = buffer.data(target + 2135);
    auto *t_2136 = buffer.data(target + 2136);
    auto *t_2137 = buffer.data(target + 2137);
    auto *t_2138 = buffer.data(target + 2138);
    auto *t_2139 = buffer.data(target + 2139);
    auto *t_2140 = buffer.data(target + 2140);
    auto *t_2141 = buffer.data(target + 2141);
    auto *t_2142 = buffer.data(target + 2142);
    auto *t_2143 = buffer.data(target + 2143);
    auto *t_2144 = buffer.data(target + 2144);
    auto *t_2145 = buffer.data(target + 2145);
    auto *t_2146 = buffer.data(target + 2146);
    auto *t_2147 = buffer.data(target + 2147);
    auto *t_2148 = buffer.data(target + 2148);
    auto *t_2149 = buffer.data(target + 2149);
    auto *t_2150 = buffer.data(target + 2150);
    auto *t_2151 = buffer.data(target + 2151);
    auto *t_2152 = buffer.data(target + 2152);
    auto *t_2153 = buffer.data(target + 2153);
    auto *t_2154 = buffer.data(target + 2154);
    auto *t_2155 = buffer.data(target + 2155);
    auto *t_2156 = buffer.data(target + 2156);
    auto *t_2157 = buffer.data(target + 2157);
    auto *t_2158 = buffer.data(target + 2158);
    auto *t_2159 = buffer.data(target + 2159);
    auto *t_2160 = buffer.data(target + 2160);
    auto *t_2161 = buffer.data(target + 2161);
    auto *t_2162 = buffer.data(target + 2162);
    auto *t_2163 = buffer.data(target + 2163);
    auto *t_2164 = buffer.data(target + 2164);
    auto *t_2165 = buffer.data(target + 2165);
    auto *t_2166 = buffer.data(target + 2166);
    auto *t_2167 = buffer.data(target + 2167);
    auto *t_2168 = buffer.data(target + 2168);
    auto *t_2169 = buffer.data(target + 2169);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_1820 = buffer.data(osi0 + 1820);
    const auto *osi0_1825 = buffer.data(osi0 + 1825);
    const auto *osi0_1829 = buffer.data(osi0 + 1829);
    const auto *osi0_1834 = buffer.data(osi0 + 1834);
    const auto *osi0_2049 = buffer.data(osi0 + 2049);
    const auto *osi0_2050 = buffer.data(osi0 + 2050);
    const auto *osi0_2053 = buffer.data(osi0 + 2053);
    const auto *osi0_2054 = buffer.data(osi0 + 2054);
    const auto *osi0_2056 = buffer.data(osi0 + 2056);
    const auto *osi0_2058 = buffer.data(osi0 + 2058);
    const auto *osi0_2065 = buffer.data(osi0 + 2065);
    const auto *osi0_2067 = buffer.data(osi0 + 2067);
    const auto *osi0_2068 = buffer.data(osi0 + 2068);
    const auto *osi0_2069 = buffer.data(osi0 + 2069);
    const auto *osi0_2071 = buffer.data(osi0 + 2071);
    const auto *osi0_2072 = buffer.data(osi0 + 2072);
    const auto *osi0_2075 = buffer.data(osi0 + 2075);
    const auto *osi0_2077 = buffer.data(osi0 + 2077);
    const auto *osi0_2078 = buffer.data(osi0 + 2078);
    const auto *osi0_2081 = buffer.data(osi0 + 2081);
    const auto *osi0_2082 = buffer.data(osi0 + 2082);
    const auto *osi0_2084 = buffer.data(osi0 + 2084);
    const auto *osi0_2086 = buffer.data(osi0 + 2086);
    const auto *osi0_2093 = buffer.data(osi0 + 2093);
    const auto *osi0_2095 = buffer.data(osi0 + 2095);
    const auto *osi0_2096 = buffer.data(osi0 + 2096);
    const auto *osi0_2097 = buffer.data(osi0 + 2097);
    const auto *osi0_2099 = buffer.data(osi0 + 2099);
    const auto *osi0_2100 = buffer.data(osi0 + 2100);
    const auto *osi0_2103 = buffer.data(osi0 + 2103);
    const auto *osi0_2105 = buffer.data(osi0 + 2105);
    const auto *osi0_2106 = buffer.data(osi0 + 2106);
    const auto *osi0_2109 = buffer.data(osi0 + 2109);
    const auto *osi0_2110 = buffer.data(osi0 + 2110);
    const auto *osi0_2112 = buffer.data(osi0 + 2112);
    const auto *osi0_2114 = buffer.data(osi0 + 2114);
    const auto *osi0_2121 = buffer.data(osi0 + 2121);
    const auto *osi0_2123 = buffer.data(osi0 + 2123);
    const auto *osi0_2124 = buffer.data(osi0 + 2124);
    const auto *osi0_2125 = buffer.data(osi0 + 2125);
    const auto *osi0_2127 = buffer.data(osi0 + 2127);
    const auto *osi0_2131 = buffer.data(osi0 + 2131);
    const auto *osi0_2134 = buffer.data(osi0 + 2134);
    const auto *osi0_2138 = buffer.data(osi0 + 2138);
    const auto *osi0_2140 = buffer.data(osi0 + 2140);
    const auto *osi0_2149 = buffer.data(osi0 + 2149);
    const auto *osi0_2151 = buffer.data(osi0 + 2151);
    const auto *osi0_2152 = buffer.data(osi0 + 2152);
    const auto *osi0_2153 = buffer.data(osi0 + 2153);
    const auto *osi0_2155 = buffer.data(osi0 + 2155);
    const auto *osi0_2156 = buffer.data(osi0 + 2156);
    const auto *osi0_2161 = buffer.data(osi0 + 2161);
    const auto *osi0_2165 = buffer.data(osi0 + 2165);

    const auto *osh_1284 = buffer.data(osh + 1284);
    const auto *osh_1287 = buffer.data(osh + 1287);
    const auto *osh_1296 = buffer.data(osh + 1296);
    const auto *osh_1302 = buffer.data(osh + 1302);
    const auto *osh_1305 = buffer.data(osh + 1305);
    const auto *osh_1307 = buffer.data(osh + 1307);
    const auto *osh_1308 = buffer.data(osh + 1308);
    const auto *osh_1311 = buffer.data(osh + 1311);
    const auto *osh_1317 = buffer.data(osh + 1317);
    const auto *osh_1322 = buffer.data(osh + 1322);
    const auto *osh_1323 = buffer.data(osh + 1323);
    const auto *osh_1325 = buffer.data(osh + 1325);
    const auto *osh_1326 = buffer.data(osh + 1326);
    const auto *osh_1328 = buffer.data(osh + 1328);
    const auto *osh_1329 = buffer.data(osh + 1329);
    const auto *osh_1332 = buffer.data(osh + 1332);
    const auto *osh_1338 = buffer.data(osh + 1338);
    const auto *osh_1343 = buffer.data(osh + 1343);
    const auto *osh_1344 = buffer.data(osh + 1344);
    const auto *osh_1346 = buffer.data(osh + 1346);
    const auto *osh_1347 = buffer.data(osh + 1347);
    const auto *osh_1349 = buffer.data(osh + 1349);
    const auto *osh_1350 = buffer.data(osh + 1350);
    const auto *osh_1353 = buffer.data(osh + 1353);
    const auto *osh_1359 = buffer.data(osh + 1359);
    const auto *osh_1364 = buffer.data(osh + 1364);
    const auto *osh_1365 = buffer.data(osh + 1365);
    const auto *osh_1367 = buffer.data(osh + 1367);
    const auto *osh_1370 = buffer.data(osh + 1370);
    const auto *osh_1374 = buffer.data(osh + 1374);
    const auto *osh_1385 = buffer.data(osh + 1385);
    const auto *osh_1538 = buffer.data(osh + 1538);
    const auto *osh_1539 = buffer.data(osh + 1539);
    const auto *osh_1542 = buffer.data(osh + 1542);
    const auto *osh_1543 = buffer.data(osh + 1543);
    const auto *osh_1545 = buffer.data(osh + 1545);
    const auto *osh_1547 = buffer.data(osh + 1547);
    const auto *osh_1548 = buffer.data(osh + 1548);
    const auto *osh_1549 = buffer.data(osh + 1549);
    const auto *osh_1550 = buffer.data(osh + 1550);
    const auto *osh_1551 = buffer.data(osh + 1551);
    const auto *osh_1552 = buffer.data(osh + 1552);
    const auto *osh_1553 = buffer.data(osh + 1553);
    const auto *osh_1554 = buffer.data(osh + 1554);
    const auto *osh_1557 = buffer.data(osh + 1557);
    const auto *osh_1559 = buffer.data(osh + 1559);
    const auto *osh_1560 = buffer.data(osh + 1560);
    const auto *osh_1563 = buffer.data(osh + 1563);
    const auto *osh_1564 = buffer.data(osh + 1564);
    const auto *osh_1566 = buffer.data(osh + 1566);
    const auto *osh_1568 = buffer.data(osh + 1568);
    const auto *osh_1569 = buffer.data(osh + 1569);
    const auto *osh_1570 = buffer.data(osh + 1570);
    const auto *osh_1571 = buffer.data(osh + 1571);
    const auto *osh_1572 = buffer.data(osh + 1572);
    const auto *osh_1573 = buffer.data(osh + 1573);
    const auto *osh_1574 = buffer.data(osh + 1574);
    const auto *osh_1575 = buffer.data(osh + 1575);
    const auto *osh_1578 = buffer.data(osh + 1578);
    const auto *osh_1580 = buffer.data(osh + 1580);
    const auto *osh_1581 = buffer.data(osh + 1581);
    const auto *osh_1584 = buffer.data(osh + 1584);
    const auto *osh_1585 = buffer.data(osh + 1585);
    const auto *osh_1587 = buffer.data(osh + 1587);
    const auto *osh_1589 = buffer.data(osh + 1589);
    const auto *osh_1590 = buffer.data(osh + 1590);
    const auto *osh_1591 = buffer.data(osh + 1591);
    const auto *osh_1592 = buffer.data(osh + 1592);
    const auto *osh_1593 = buffer.data(osh + 1593);
    const auto *osh_1594 = buffer.data(osh + 1594);
    const auto *osh_1595 = buffer.data(osh + 1595);
    const auto *osh_1599 = buffer.data(osh + 1599);
    const auto *osh_1602 = buffer.data(osh + 1602);
    const auto *osh_1606 = buffer.data(osh + 1606);
    const auto *osh_1608 = buffer.data(osh + 1608);
    const auto *osh_1611 = buffer.data(osh + 1611);
    const auto *osh_1612 = buffer.data(osh + 1612);
    const auto *osh_1613 = buffer.data(osh + 1613);
    const auto *osh_1614 = buffer.data(osh + 1614);
    const auto *osh_1615 = buffer.data(osh + 1615);
    const auto *osh_1616 = buffer.data(osh + 1616);
    const auto *osh_1617 = buffer.data(osh + 1617);
    const auto *osh_1622 = buffer.data(osh + 1622);
    const auto *osh_1626 = buffer.data(osh + 1626);

    const auto *osi1_1820 = buffer.data(osi1 + 1820);
    const auto *osi1_1825 = buffer.data(osi1 + 1825);
    const auto *osi1_1829 = buffer.data(osi1 + 1829);
    const auto *osi1_1834 = buffer.data(osi1 + 1834);
    const auto *osi1_2049 = buffer.data(osi1 + 2049);
    const auto *osi1_2050 = buffer.data(osi1 + 2050);
    const auto *osi1_2053 = buffer.data(osi1 + 2053);
    const auto *osi1_2054 = buffer.data(osi1 + 2054);
    const auto *osi1_2056 = buffer.data(osi1 + 2056);
    const auto *osi1_2058 = buffer.data(osi1 + 2058);
    const auto *osi1_2065 = buffer.data(osi1 + 2065);
    const auto *osi1_2067 = buffer.data(osi1 + 2067);
    const auto *osi1_2068 = buffer.data(osi1 + 2068);
    const auto *osi1_2069 = buffer.data(osi1 + 2069);
    const auto *osi1_2071 = buffer.data(osi1 + 2071);
    const auto *osi1_2072 = buffer.data(osi1 + 2072);
    const auto *osi1_2075 = buffer.data(osi1 + 2075);
    const auto *osi1_2077 = buffer.data(osi1 + 2077);
    const auto *osi1_2078 = buffer.data(osi1 + 2078);
    const auto *osi1_2081 = buffer.data(osi1 + 2081);
    const auto *osi1_2082 = buffer.data(osi1 + 2082);
    const auto *osi1_2084 = buffer.data(osi1 + 2084);
    const auto *osi1_2086 = buffer.data(osi1 + 2086);
    const auto *osi1_2093 = buffer.data(osi1 + 2093);
    const auto *osi1_2095 = buffer.data(osi1 + 2095);
    const auto *osi1_2096 = buffer.data(osi1 + 2096);
    const auto *osi1_2097 = buffer.data(osi1 + 2097);
    const auto *osi1_2099 = buffer.data(osi1 + 2099);
    const auto *osi1_2100 = buffer.data(osi1 + 2100);
    const auto *osi1_2103 = buffer.data(osi1 + 2103);
    const auto *osi1_2105 = buffer.data(osi1 + 2105);
    const auto *osi1_2106 = buffer.data(osi1 + 2106);
    const auto *osi1_2109 = buffer.data(osi1 + 2109);
    const auto *osi1_2110 = buffer.data(osi1 + 2110);
    const auto *osi1_2112 = buffer.data(osi1 + 2112);
    const auto *osi1_2114 = buffer.data(osi1 + 2114);
    const auto *osi1_2121 = buffer.data(osi1 + 2121);
    const auto *osi1_2123 = buffer.data(osi1 + 2123);
    const auto *osi1_2124 = buffer.data(osi1 + 2124);
    const auto *osi1_2125 = buffer.data(osi1 + 2125);
    const auto *osi1_2127 = buffer.data(osi1 + 2127);
    const auto *osi1_2131 = buffer.data(osi1 + 2131);
    const auto *osi1_2134 = buffer.data(osi1 + 2134);
    const auto *osi1_2138 = buffer.data(osi1 + 2138);
    const auto *osi1_2140 = buffer.data(osi1 + 2140);
    const auto *osi1_2149 = buffer.data(osi1 + 2149);
    const auto *osi1_2151 = buffer.data(osi1 + 2151);
    const auto *osi1_2152 = buffer.data(osi1 + 2152);
    const auto *osi1_2153 = buffer.data(osi1 + 2153);
    const auto *osi1_2155 = buffer.data(osi1 + 2155);
    const auto *osi1_2156 = buffer.data(osi1 + 2156);
    const auto *osi1_2161 = buffer.data(osi1 + 2161);
    const auto *osi1_2165 = buffer.data(osi1 + 2165);

    const auto *qsg0_1155 = buffer.data(qsg0 + 1155);
    const auto *qsg0_1156 = buffer.data(qsg0 + 1156);
    const auto *qsg0_1157 = buffer.data(qsg0 + 1157);
    const auto *qsg0_1158 = buffer.data(qsg0 + 1158);
    const auto *qsg0_1159 = buffer.data(qsg0 + 1159);
    const auto *qsg0_1160 = buffer.data(qsg0 + 1160);

    const auto *qsg1_1155 = buffer.data(qsg1 + 1155);
    const auto *qsg1_1156 = buffer.data(qsg1 + 1156);
    const auto *qsg1_1157 = buffer.data(qsg1 + 1157);
    const auto *qsg1_1158 = buffer.data(qsg1 + 1158);
    const auto *qsg1_1159 = buffer.data(qsg1 + 1159);
    const auto *qsg1_1160 = buffer.data(qsg1 + 1160);

    const auto *qsh_1536 = buffer.data(qsh + 1536);
    const auto *qsh_1538 = buffer.data(qsh + 1538);
    const auto *qsh_1539 = buffer.data(qsh + 1539);
    const auto *qsh_1542 = buffer.data(qsh + 1542);
    const auto *qsh_1548 = buffer.data(qsh + 1548);
    const auto *qsh_1549 = buffer.data(qsh + 1549);
    const auto *qsh_1550 = buffer.data(qsh + 1550);
    const auto *qsh_1551 = buffer.data(qsh + 1551);
    const auto *qsh_1552 = buffer.data(qsh + 1552);
    const auto *qsh_1553 = buffer.data(qsh + 1553);
    const auto *qsh_1554 = buffer.data(qsh + 1554);
    const auto *qsh_1556 = buffer.data(qsh + 1556);
    const auto *qsh_1557 = buffer.data(qsh + 1557);
    const auto *qsh_1559 = buffer.data(qsh + 1559);
    const auto *qsh_1560 = buffer.data(qsh + 1560);
    const auto *qsh_1563 = buffer.data(qsh + 1563);
    const auto *qsh_1569 = buffer.data(qsh + 1569);
    const auto *qsh_1570 = buffer.data(qsh + 1570);
    const auto *qsh_1571 = buffer.data(qsh + 1571);
    const auto *qsh_1572 = buffer.data(qsh + 1572);
    const auto *qsh_1573 = buffer.data(qsh + 1573);
    const auto *qsh_1574 = buffer.data(qsh + 1574);
    const auto *qsh_1575 = buffer.data(qsh + 1575);
    const auto *qsh_1577 = buffer.data(qsh + 1577);
    const auto *qsh_1578 = buffer.data(qsh + 1578);
    const auto *qsh_1580 = buffer.data(qsh + 1580);
    const auto *qsh_1581 = buffer.data(qsh + 1581);
    const auto *qsh_1584 = buffer.data(qsh + 1584);
    const auto *qsh_1590 = buffer.data(qsh + 1590);
    const auto *qsh_1591 = buffer.data(qsh + 1591);
    const auto *qsh_1592 = buffer.data(qsh + 1592);
    const auto *qsh_1593 = buffer.data(qsh + 1593);
    const auto *qsh_1594 = buffer.data(qsh + 1594);
    const auto *qsh_1595 = buffer.data(qsh + 1595);
    const auto *qsh_1596 = buffer.data(qsh + 1596);
    const auto *qsh_1598 = buffer.data(qsh + 1598);
    const auto *qsh_1599 = buffer.data(qsh + 1599);
    const auto *qsh_1601 = buffer.data(qsh + 1601);
    const auto *qsh_1602 = buffer.data(qsh + 1602);
    const auto *qsh_1605 = buffer.data(qsh + 1605);
    const auto *qsh_1611 = buffer.data(qsh + 1611);
    const auto *qsh_1612 = buffer.data(qsh + 1612);
    const auto *qsh_1613 = buffer.data(qsh + 1613);
    const auto *qsh_1614 = buffer.data(qsh + 1614);
    const auto *qsh_1615 = buffer.data(qsh + 1615);
    const auto *qsh_1616 = buffer.data(qsh + 1616);
    const auto *qsh_1617 = buffer.data(qsh + 1617);
    const auto *qsh_1618 = buffer.data(qsh + 1618);
    const auto *qsh_1619 = buffer.data(qsh + 1619);
    const auto *qsh_1620 = buffer.data(qsh + 1620);
    const auto *qsh_1621 = buffer.data(qsh + 1621);
    const auto *qsh_1622 = buffer.data(qsh + 1622);
    const auto *qsh_1623 = buffer.data(qsh + 1623);
    const auto *qsh_1624 = buffer.data(qsh + 1624);
    const auto *qsh_1625 = buffer.data(qsh + 1625);
    const auto *qsh_1626 = buffer.data(qsh + 1626);

#pragma omp simd aligned(t_2049, t_2050, t_2051, pa_x, pc_x, pc_z, osi0_2049, osi0_2050, \
                         osh_1284, osh_1538, osh_1539, osi1_2049, osi1_2050, \
                         qsh_1536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2049[k] = pa_x[k] * osi0_2049[k]
                    + f_14 * osh_1538[k]
                    - f_10 * pc_x[k] * osi1_2049[k];

        t_2050[k] = pa_x[k] * osi0_2050[k]
                    + f_13 * osh_1539[k]
                    - f_10 * pc_x[k] * osi1_2050[k];

        t_2051[k] = f_21 * osh_1284[k]
                    + f_3 * pc_z[k] * qsh_1536[k];
    }

#pragma omp simd aligned(t_2052, t_2053, t_2054, pa_x, pc_x, pc_y, osi0_2053, osi0_2054, \
                         osh_1307, osh_1542, osh_1543, osi1_2053, osi1_2054, \
                         qsh_1538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2052[k] = f_14 * osh_1307[k]
                    + f_3 * pc_y[k] * qsh_1538[k];

        t_2053[k] = pa_x[k] * osi0_2053[k]
                    + f_13 * osh_1542[k]
                    - f_10 * pc_x[k] * osi1_2053[k];

        t_2054[k] = pa_x[k] * osi0_2054[k]
                    + f_12 * osh_1543[k]
                    - f_10 * pc_x[k] * osi1_2054[k];
    }

#pragma omp simd aligned(t_2055, t_2056, t_2057, pa_x, pc_x, pc_y, pc_z, osi0_2056, osh_1287, \
                         osh_1311, osh_1545, osi1_2056, qsh_1539, \
                         qsh_1542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2055[k] = f_21 * osh_1287[k]
                    + f_3 * pc_z[k] * qsh_1539[k];

        t_2056[k] = pa_x[k] * osi0_2056[k]
                    + f_12 * osh_1545[k]
                    - f_10 * pc_x[k] * osi1_2056[k];

        t_2057[k] = f_14 * osh_1311[k]
                    + f_3 * pc_y[k] * qsh_1542[k];
    }

#pragma omp simd aligned(t_2058, t_2059, t_2060, t_2061, pa_x, pc_x, osi0_2058, osh_1547, \
                         osh_1548, osh_1549, osh_1550, osi1_2058, qsh_1548, qsh_1549, \
                         qsh_1550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2058[k] = pa_x[k] * osi0_2058[k]
                    + f_12 * osh_1547[k]
                    - f_10 * pc_x[k] * osi1_2058[k];

        t_2059[k] = f_11 * osh_1548[k]
                    + f_3 * pc_x[k] * qsh_1548[k];

        t_2060[k] = f_11 * osh_1549[k]
                    + f_3 * pc_x[k] * qsh_1549[k];

        t_2061[k] = f_11 * osh_1550[k]
                    + f_3 * pc_x[k] * qsh_1550[k];
    }

#pragma omp simd aligned(t_2062, t_2063, t_2064, t_2065, pa_x, pc_x, osi0_2065, osh_1551, \
                         osh_1552, osh_1553, osi1_2065, qsh_1551, qsh_1552, \
                         qsh_1553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2062[k] = f_11 * osh_1551[k]
                    + f_3 * pc_x[k] * qsh_1551[k];

        t_2063[k] = f_11 * osh_1552[k]
                    + f_3 * pc_x[k] * qsh_1552[k];

        t_2064[k] = f_11 * osh_1553[k]
                    + f_3 * pc_x[k] * qsh_1553[k];

        t_2065[k] = pa_x[k] * osi0_2065[k]
                    - f_10 * pc_x[k] * osi1_2065[k];
    }

#pragma omp simd aligned(t_2066, t_2067, t_2068, t_2069, pa_x, pc_x, pc_z, osi0_2067, \
                         osi0_2068, osi0_2069, osh_1296, osi1_2067, osi1_2068, osi1_2069, \
                         qsh_1548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2066[k] = f_21 * osh_1296[k]
                    + f_3 * pc_z[k] * qsh_1548[k];

        t_2067[k] = pa_x[k] * osi0_2067[k]
                    - f_10 * pc_x[k] * osi1_2067[k];

        t_2068[k] = pa_x[k] * osi0_2068[k]
                    - f_10 * pc_x[k] * osi1_2068[k];

        t_2069[k] = pa_x[k] * osi0_2069[k]
                    - f_10 * pc_x[k] * osi1_2069[k];
    }

#pragma omp simd aligned(t_2070, t_2071, t_2072, t_2073, pa_x, pc_x, pc_y, osi0_2071, \
                         osi0_2072, osh_1322, osh_1323, osh_1554, osi1_2071, osi1_2072, \
                         qsh_1553, qsh_1554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2070[k] = f_14 * osh_1322[k]
                    + f_3 * pc_y[k] * qsh_1553[k];

        t_2071[k] = pa_x[k] * osi0_2071[k]
                    - f_10 * pc_x[k] * osi1_2071[k];

        t_2072[k] = pa_x[k] * osi0_2072[k]
                    + f_23 * osh_1554[k]
                    - f_10 * pc_x[k] * osi1_2072[k];

        t_2073[k] = f_13 * osh_1323[k]
                    + f_3 * pc_y[k] * qsh_1554[k];
    }

#pragma omp simd aligned(t_2074, t_2075, t_2076, pa_x, pc_x, pc_y, pc_z, osi0_2075, osh_1302, \
                         osh_1325, osh_1557, osi1_2075, qsh_1554, \
                         qsh_1556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2074[k] = f_20 * osh_1302[k]
                    + f_3 * pc_z[k] * qsh_1554[k];

        t_2075[k] = pa_x[k] * osi0_2075[k]
                    + f_14 * osh_1557[k]
                    - f_10 * pc_x[k] * osi1_2075[k];

        t_2076[k] = f_13 * osh_1325[k]
                    + f_3 * pc_y[k] * qsh_1556[k];
    }

#pragma omp simd aligned(t_2077, t_2078, t_2079, pa_x, pc_x, pc_z, osi0_2077, osi0_2078, \
                         osh_1305, osh_1559, osh_1560, osi1_2077, osi1_2078, \
                         qsh_1557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2077[k] = pa_x[k] * osi0_2077[k]
                    + f_14 * osh_1559[k]
                    - f_10 * pc_x[k] * osi1_2077[k];

        t_2078[k] = pa_x[k] * osi0_2078[k]
                    + f_13 * osh_1560[k]
                    - f_10 * pc_x[k] * osi1_2078[k];

        t_2079[k] = f_20 * osh_1305[k]
                    + f_3 * pc_z[k] * qsh_1557[k];
    }

#pragma omp simd aligned(t_2080, t_2081, t_2082, pa_x, pc_x, pc_y, osi0_2081, osi0_2082, \
                         osh_1328, osh_1563, osh_1564, osi1_2081, osi1_2082, \
                         qsh_1559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2080[k] = f_13 * osh_1328[k]
                    + f_3 * pc_y[k] * qsh_1559[k];

        t_2081[k] = pa_x[k] * osi0_2081[k]
                    + f_13 * osh_1563[k]
                    - f_10 * pc_x[k] * osi1_2081[k];

        t_2082[k] = pa_x[k] * osi0_2082[k]
                    + f_12 * osh_1564[k]
                    - f_10 * pc_x[k] * osi1_2082[k];
    }

#pragma omp simd aligned(t_2083, t_2084, t_2085, pa_x, pc_x, pc_y, pc_z, osi0_2084, osh_1308, \
                         osh_1332, osh_1566, osi1_2084, qsh_1560, \
                         qsh_1563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2083[k] = f_20 * osh_1308[k]
                    + f_3 * pc_z[k] * qsh_1560[k];

        t_2084[k] = pa_x[k] * osi0_2084[k]
                    + f_12 * osh_1566[k]
                    - f_10 * pc_x[k] * osi1_2084[k];

        t_2085[k] = f_13 * osh_1332[k]
                    + f_3 * pc_y[k] * qsh_1563[k];
    }

#pragma omp simd aligned(t_2086, t_2087, t_2088, t_2089, pa_x, pc_x, osi0_2086, osh_1568, \
                         osh_1569, osh_1570, osh_1571, osi1_2086, qsh_1569, qsh_1570, \
                         qsh_1571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2086[k] = pa_x[k] * osi0_2086[k]
                    + f_12 * osh_1568[k]
                    - f_10 * pc_x[k] * osi1_2086[k];

        t_2087[k] = f_11 * osh_1569[k]
                    + f_3 * pc_x[k] * qsh_1569[k];

        t_2088[k] = f_11 * osh_1570[k]
                    + f_3 * pc_x[k] * qsh_1570[k];

        t_2089[k] = f_11 * osh_1571[k]
                    + f_3 * pc_x[k] * qsh_1571[k];
    }

#pragma omp simd aligned(t_2090, t_2091, t_2092, t_2093, pa_x, pc_x, osi0_2093, osh_1572, \
                         osh_1573, osh_1574, osi1_2093, qsh_1572, qsh_1573, \
                         qsh_1574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2090[k] = f_11 * osh_1572[k]
                    + f_3 * pc_x[k] * qsh_1572[k];

        t_2091[k] = f_11 * osh_1573[k]
                    + f_3 * pc_x[k] * qsh_1573[k];

        t_2092[k] = f_11 * osh_1574[k]
                    + f_3 * pc_x[k] * qsh_1574[k];

        t_2093[k] = pa_x[k] * osi0_2093[k]
                    - f_10 * pc_x[k] * osi1_2093[k];
    }

#pragma omp simd aligned(t_2094, t_2095, t_2096, t_2097, pa_x, pc_x, pc_z, osi0_2095, \
                         osi0_2096, osi0_2097, osh_1317, osi1_2095, osi1_2096, osi1_2097, \
                         qsh_1569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2094[k] = f_20 * osh_1317[k]
                    + f_3 * pc_z[k] * qsh_1569[k];

        t_2095[k] = pa_x[k] * osi0_2095[k]
                    - f_10 * pc_x[k] * osi1_2095[k];

        t_2096[k] = pa_x[k] * osi0_2096[k]
                    - f_10 * pc_x[k] * osi1_2096[k];

        t_2097[k] = pa_x[k] * osi0_2097[k]
                    - f_10 * pc_x[k] * osi1_2097[k];
    }

#pragma omp simd aligned(t_2098, t_2099, t_2100, t_2101, pa_x, pc_x, pc_y, osi0_2099, \
                         osi0_2100, osh_1343, osh_1344, osh_1575, osi1_2099, osi1_2100, \
                         qsh_1574, qsh_1575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2098[k] = f_13 * osh_1343[k]
                    + f_3 * pc_y[k] * qsh_1574[k];

        t_2099[k] = pa_x[k] * osi0_2099[k]
                    - f_10 * pc_x[k] * osi1_2099[k];

        t_2100[k] = pa_x[k] * osi0_2100[k]
                    + f_23 * osh_1575[k]
                    - f_10 * pc_x[k] * osi1_2100[k];

        t_2101[k] = f_12 * osh_1344[k]
                    + f_3 * pc_y[k] * qsh_1575[k];
    }

#pragma omp simd aligned(t_2102, t_2103, t_2104, pa_x, pc_x, pc_y, pc_z, osi0_2103, osh_1323, \
                         osh_1346, osh_1578, osi1_2103, qsh_1575, \
                         qsh_1577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2102[k] = f_19 * osh_1323[k]
                    + f_3 * pc_z[k] * qsh_1575[k];

        t_2103[k] = pa_x[k] * osi0_2103[k]
                    + f_14 * osh_1578[k]
                    - f_10 * pc_x[k] * osi1_2103[k];

        t_2104[k] = f_12 * osh_1346[k]
                    + f_3 * pc_y[k] * qsh_1577[k];
    }

#pragma omp simd aligned(t_2105, t_2106, t_2107, pa_x, pc_x, pc_z, osi0_2105, osi0_2106, \
                         osh_1326, osh_1580, osh_1581, osi1_2105, osi1_2106, \
                         qsh_1578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2105[k] = pa_x[k] * osi0_2105[k]
                    + f_14 * osh_1580[k]
                    - f_10 * pc_x[k] * osi1_2105[k];

        t_2106[k] = pa_x[k] * osi0_2106[k]
                    + f_13 * osh_1581[k]
                    - f_10 * pc_x[k] * osi1_2106[k];

        t_2107[k] = f_19 * osh_1326[k]
                    + f_3 * pc_z[k] * qsh_1578[k];
    }

#pragma omp simd aligned(t_2108, t_2109, t_2110, pa_x, pc_x, pc_y, osi0_2109, osi0_2110, \
                         osh_1349, osh_1584, osh_1585, osi1_2109, osi1_2110, \
                         qsh_1580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2108[k] = f_12 * osh_1349[k]
                    + f_3 * pc_y[k] * qsh_1580[k];

        t_2109[k] = pa_x[k] * osi0_2109[k]
                    + f_13 * osh_1584[k]
                    - f_10 * pc_x[k] * osi1_2109[k];

        t_2110[k] = pa_x[k] * osi0_2110[k]
                    + f_12 * osh_1585[k]
                    - f_10 * pc_x[k] * osi1_2110[k];
    }

#pragma omp simd aligned(t_2111, t_2112, t_2113, pa_x, pc_x, pc_y, pc_z, osi0_2112, osh_1329, \
                         osh_1353, osh_1587, osi1_2112, qsh_1581, \
                         qsh_1584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2111[k] = f_19 * osh_1329[k]
                    + f_3 * pc_z[k] * qsh_1581[k];

        t_2112[k] = pa_x[k] * osi0_2112[k]
                    + f_12 * osh_1587[k]
                    - f_10 * pc_x[k] * osi1_2112[k];

        t_2113[k] = f_12 * osh_1353[k]
                    + f_3 * pc_y[k] * qsh_1584[k];
    }

#pragma omp simd aligned(t_2114, t_2115, t_2116, t_2117, pa_x, pc_x, osi0_2114, osh_1589, \
                         osh_1590, osh_1591, osh_1592, osi1_2114, qsh_1590, qsh_1591, \
                         qsh_1592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2114[k] = pa_x[k] * osi0_2114[k]
                    + f_12 * osh_1589[k]
                    - f_10 * pc_x[k] * osi1_2114[k];

        t_2115[k] = f_11 * osh_1590[k]
                    + f_3 * pc_x[k] * qsh_1590[k];

        t_2116[k] = f_11 * osh_1591[k]
                    + f_3 * pc_x[k] * qsh_1591[k];

        t_2117[k] = f_11 * osh_1592[k]
                    + f_3 * pc_x[k] * qsh_1592[k];
    }

#pragma omp simd aligned(t_2118, t_2119, t_2120, t_2121, pa_x, pc_x, osi0_2121, osh_1593, \
                         osh_1594, osh_1595, osi1_2121, qsh_1593, qsh_1594, \
                         qsh_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2118[k] = f_11 * osh_1593[k]
                    + f_3 * pc_x[k] * qsh_1593[k];

        t_2119[k] = f_11 * osh_1594[k]
                    + f_3 * pc_x[k] * qsh_1594[k];

        t_2120[k] = f_11 * osh_1595[k]
                    + f_3 * pc_x[k] * qsh_1595[k];

        t_2121[k] = pa_x[k] * osi0_2121[k]
                    - f_10 * pc_x[k] * osi1_2121[k];
    }

#pragma omp simd aligned(t_2122, t_2123, t_2124, t_2125, pa_x, pc_x, pc_z, osi0_2123, \
                         osi0_2124, osi0_2125, osh_1338, osi1_2123, osi1_2124, osi1_2125, \
                         qsh_1590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2122[k] = f_19 * osh_1338[k]
                    + f_3 * pc_z[k] * qsh_1590[k];

        t_2123[k] = pa_x[k] * osi0_2123[k]
                    - f_10 * pc_x[k] * osi1_2123[k];

        t_2124[k] = pa_x[k] * osi0_2124[k]
                    - f_10 * pc_x[k] * osi1_2124[k];

        t_2125[k] = pa_x[k] * osi0_2125[k]
                    - f_10 * pc_x[k] * osi1_2125[k];
    }

#pragma omp simd aligned(t_2126, t_2127, t_2128, t_2129, pa_x, pa_y, pc_x, pc_y, osi0_1820, \
                         osi0_2127, osh_1364, osh_1365, osi1_1820, osi1_2127, qsh_1595, \
                         qsh_1596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2126[k] = f_12 * osh_1364[k]
                    + f_3 * pc_y[k] * qsh_1595[k];

        t_2127[k] = pa_x[k] * osi0_2127[k]
                    - f_10 * pc_x[k] * osi1_2127[k];

        t_2128[k] = pa_y[k] * osi0_1820[k]
                    - f_10 * pc_y[k] * osi1_1820[k];

        t_2129[k] = f_11 * osh_1365[k]
                    + f_3 * pc_y[k] * qsh_1596[k];
    }

#pragma omp simd aligned(t_2130, t_2131, t_2132, pa_x, pc_x, pc_y, pc_z, osi0_2131, osh_1344, \
                         osh_1367, osh_1599, osi1_2131, qsh_1596, \
                         qsh_1598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2130[k] = f_18 * osh_1344[k]
                    + f_3 * pc_z[k] * qsh_1596[k];

        t_2131[k] = pa_x[k] * osi0_2131[k]
                    + f_14 * osh_1599[k]
                    - f_10 * pc_x[k] * osi1_2131[k];

        t_2132[k] = f_11 * osh_1367[k]
                    + f_3 * pc_y[k] * qsh_1598[k];
    }

#pragma omp simd aligned(t_2133, t_2134, t_2135, pa_x, pa_y, pc_x, pc_y, pc_z, osi0_1825, \
                         osi0_2134, osh_1347, osh_1602, osi1_1825, osi1_2134, \
                         qsh_1599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2133[k] = pa_y[k] * osi0_1825[k]
                    - f_10 * pc_y[k] * osi1_1825[k];

        t_2134[k] = pa_x[k] * osi0_2134[k]
                    + f_13 * osh_1602[k]
                    - f_10 * pc_x[k] * osi1_2134[k];

        t_2135[k] = f_18 * osh_1347[k]
                    + f_3 * pc_z[k] * qsh_1599[k];
    }

#pragma omp simd aligned(t_2136, t_2137, t_2138, pa_x, pa_y, pc_x, pc_y, osi0_1829, osi0_2138, \
                         osh_1370, osh_1606, osi1_1829, osi1_2138, \
                         qsh_1601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2136[k] = f_11 * osh_1370[k]
                    + f_3 * pc_y[k] * qsh_1601[k];

        t_2137[k] = pa_y[k] * osi0_1829[k]
                    - f_10 * pc_y[k] * osi1_1829[k];

        t_2138[k] = pa_x[k] * osi0_2138[k]
                    + f_12 * osh_1606[k]
                    - f_10 * pc_x[k] * osi1_2138[k];
    }

#pragma omp simd aligned(t_2139, t_2140, t_2141, pa_x, pc_x, pc_y, pc_z, osi0_2140, osh_1350, \
                         osh_1374, osh_1608, osi1_2140, qsh_1602, \
                         qsh_1605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2139[k] = f_18 * osh_1350[k]
                    + f_3 * pc_z[k] * qsh_1602[k];

        t_2140[k] = pa_x[k] * osi0_2140[k]
                    + f_12 * osh_1608[k]
                    - f_10 * pc_x[k] * osi1_2140[k];

        t_2141[k] = f_11 * osh_1374[k]
                    + f_3 * pc_y[k] * qsh_1605[k];
    }

#pragma omp simd aligned(t_2142, t_2143, t_2144, t_2145, pa_y, pc_x, pc_y, osi0_1834, \
                         osh_1611, osh_1612, osh_1613, osi1_1834, qsh_1611, qsh_1612, \
                         qsh_1613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2142[k] = pa_y[k] * osi0_1834[k]
                    - f_10 * pc_y[k] * osi1_1834[k];

        t_2143[k] = f_11 * osh_1611[k]
                    + f_3 * pc_x[k] * qsh_1611[k];

        t_2144[k] = f_11 * osh_1612[k]
                    + f_3 * pc_x[k] * qsh_1612[k];

        t_2145[k] = f_11 * osh_1613[k]
                    + f_3 * pc_x[k] * qsh_1613[k];
    }

#pragma omp simd aligned(t_2146, t_2147, t_2148, t_2149, pa_x, pc_x, osi0_2149, osh_1614, \
                         osh_1615, osh_1616, osi1_2149, qsh_1614, qsh_1615, \
                         qsh_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2146[k] = f_11 * osh_1614[k]
                    + f_3 * pc_x[k] * qsh_1614[k];

        t_2147[k] = f_11 * osh_1615[k]
                    + f_3 * pc_x[k] * qsh_1615[k];

        t_2148[k] = f_11 * osh_1616[k]
                    + f_3 * pc_x[k] * qsh_1616[k];

        t_2149[k] = pa_x[k] * osi0_2149[k]
                    - f_10 * pc_x[k] * osi1_2149[k];
    }

#pragma omp simd aligned(t_2150, t_2151, t_2152, t_2153, pa_x, pc_x, pc_z, osi0_2151, \
                         osi0_2152, osi0_2153, osh_1359, osi1_2151, osi1_2152, osi1_2153, \
                         qsh_1611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2150[k] = f_18 * osh_1359[k]
                    + f_3 * pc_z[k] * qsh_1611[k];

        t_2151[k] = pa_x[k] * osi0_2151[k]
                    - f_10 * pc_x[k] * osi1_2151[k];

        t_2152[k] = pa_x[k] * osi0_2152[k]
                    - f_10 * pc_x[k] * osi1_2152[k];

        t_2153[k] = pa_x[k] * osi0_2153[k]
                    - f_10 * pc_x[k] * osi1_2153[k];
    }

#pragma omp simd aligned(t_2154, t_2155, t_2156, t_2157, pa_x, pc_x, pc_y, osi0_2155, \
                         osi0_2156, osh_1385, osh_1617, osi1_2155, osi1_2156, qsh_1616, \
                         qsh_1617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2154[k] = f_11 * osh_1385[k]
                    + f_3 * pc_y[k] * qsh_1616[k];

        t_2155[k] = pa_x[k] * osi0_2155[k]
                    - f_10 * pc_x[k] * osi1_2155[k];

        t_2156[k] = pa_x[k] * osi0_2156[k]
                    + f_23 * osh_1617[k]
                    - f_10 * pc_x[k] * osi1_2156[k];

        t_2157[k] = f_3 * pc_y[k] * qsh_1617[k];
    }

#pragma omp simd aligned(t_2158, t_2159, t_2160, pc_y, pc_z, osh_1365, qsg0_1155, qsg1_1155, \
                         qsh_1617, qsh_1618, qsh_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2158[k] = f_15 * osh_1365[k]
                    + f_3 * pc_z[k] * qsh_1617[k];

        t_2159[k] = f_4 * qsg0_1155[k]
                    - f_5 * qsg1_1155[k]
                    + f_3 * pc_y[k] * qsh_1618[k];

        t_2160[k] = f_3 * pc_y[k] * qsh_1619[k];
    }

#pragma omp simd aligned(t_2161, t_2162, t_2163, pa_x, pc_x, pc_y, osi0_2161, osh_1622, \
                         osi1_2161, qsg0_1156, qsg0_1157, qsg1_1156, qsg1_1157, qsh_1620, \
                         qsh_1621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2161[k] = pa_x[k] * osi0_2161[k]
                    + f_14 * osh_1622[k]
                    - f_10 * pc_x[k] * osi1_2161[k];

        t_2162[k] = f_6 * qsg0_1156[k]
                    - f_7 * qsg1_1156[k]
                    + f_3 * pc_y[k] * qsh_1620[k];

        t_2163[k] = f_4 * qsg0_1157[k]
                    - f_5 * qsg1_1157[k]
                    + f_3 * pc_y[k] * qsh_1621[k];
    }

#pragma omp simd aligned(t_2164, t_2165, t_2166, pa_x, pc_x, pc_y, osi0_2165, osh_1626, \
                         osi1_2165, qsg0_1158, qsg1_1158, qsh_1622, \
                         qsh_1623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2164[k] = f_3 * pc_y[k] * qsh_1622[k];

        t_2165[k] = pa_x[k] * osi0_2165[k]
                    + f_13 * osh_1626[k]
                    - f_10 * pc_x[k] * osi1_2165[k];

        t_2166[k] = f_8 * qsg0_1158[k]
                    - f_9 * qsg1_1158[k]
                    + f_3 * pc_y[k] * qsh_1623[k];
    }

#pragma omp simd aligned(t_2167, t_2168, t_2169, pc_y, qsg0_1159, qsg0_1160, qsg1_1159, \
                         qsg1_1160, qsh_1624, qsh_1625, qsh_1626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2167[k] = f_6 * qsg0_1159[k]
                    - f_7 * qsg1_1159[k]
                    + f_3 * pc_y[k] * qsh_1624[k];

        t_2168[k] = f_4 * qsg0_1160[k]
                    - f_5 * qsg1_1160[k]
                    + f_3 * pc_y[k] * qsh_1625[k];

        t_2169[k] = f_3 * pc_y[k] * qsh_1626[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece19(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsg0, const size_t qsg1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 5.0 / q;
    const auto f_19 = 4.5 / q;

    auto *t_2170 = buffer.data(target + 2170);
    auto *t_2171 = buffer.data(target + 2171);
    auto *t_2172 = buffer.data(target + 2172);
    auto *t_2173 = buffer.data(target + 2173);
    auto *t_2174 = buffer.data(target + 2174);
    auto *t_2175 = buffer.data(target + 2175);
    auto *t_2176 = buffer.data(target + 2176);
    auto *t_2177 = buffer.data(target + 2177);
    auto *t_2178 = buffer.data(target + 2178);
    auto *t_2179 = buffer.data(target + 2179);
    auto *t_2180 = buffer.data(target + 2180);
    auto *t_2181 = buffer.data(target + 2181);
    auto *t_2182 = buffer.data(target + 2182);
    auto *t_2183 = buffer.data(target + 2183);
    auto *t_2184 = buffer.data(target + 2184);
    auto *t_2185 = buffer.data(target + 2185);
    auto *t_2186 = buffer.data(target + 2186);
    auto *t_2187 = buffer.data(target + 2187);
    auto *t_2188 = buffer.data(target + 2188);
    auto *t_2189 = buffer.data(target + 2189);
    auto *t_2190 = buffer.data(target + 2190);
    auto *t_2191 = buffer.data(target + 2191);
    auto *t_2192 = buffer.data(target + 2192);
    auto *t_2193 = buffer.data(target + 2193);
    auto *t_2194 = buffer.data(target + 2194);
    auto *t_2195 = buffer.data(target + 2195);
    auto *t_2196 = buffer.data(target + 2196);
    auto *t_2197 = buffer.data(target + 2197);
    auto *t_2198 = buffer.data(target + 2198);
    auto *t_2199 = buffer.data(target + 2199);
    auto *t_2200 = buffer.data(target + 2200);
    auto *t_2201 = buffer.data(target + 2201);
    auto *t_2202 = buffer.data(target + 2202);
    auto *t_2203 = buffer.data(target + 2203);
    auto *t_2204 = buffer.data(target + 2204);
    auto *t_2205 = buffer.data(target + 2205);
    auto *t_2206 = buffer.data(target + 2206);
    auto *t_2207 = buffer.data(target + 2207);
    auto *t_2208 = buffer.data(target + 2208);
    auto *t_2209 = buffer.data(target + 2209);
    auto *t_2210 = buffer.data(target + 2210);
    auto *t_2211 = buffer.data(target + 2211);
    auto *t_2212 = buffer.data(target + 2212);
    auto *t_2213 = buffer.data(target + 2213);
    auto *t_2214 = buffer.data(target + 2214);
    auto *t_2215 = buffer.data(target + 2215);
    auto *t_2216 = buffer.data(target + 2216);
    auto *t_2217 = buffer.data(target + 2217);
    auto *t_2218 = buffer.data(target + 2218);
    auto *t_2219 = buffer.data(target + 2219);
    auto *t_2220 = buffer.data(target + 2220);
    auto *t_2221 = buffer.data(target + 2221);
    auto *t_2222 = buffer.data(target + 2222);
    auto *t_2223 = buffer.data(target + 2223);
    auto *t_2224 = buffer.data(target + 2224);
    auto *t_2225 = buffer.data(target + 2225);
    auto *t_2226 = buffer.data(target + 2226);
    auto *t_2227 = buffer.data(target + 2227);
    auto *t_2228 = buffer.data(target + 2228);
    auto *t_2229 = buffer.data(target + 2229);
    auto *t_2230 = buffer.data(target + 2230);
    auto *t_2231 = buffer.data(target + 2231);
    auto *t_2232 = buffer.data(target + 2232);
    auto *t_2233 = buffer.data(target + 2233);
    auto *t_2234 = buffer.data(target + 2234);
    auto *t_2235 = buffer.data(target + 2235);
    auto *t_2236 = buffer.data(target + 2236);
    auto *t_2237 = buffer.data(target + 2237);
    auto *t_2238 = buffer.data(target + 2238);
    auto *t_2239 = buffer.data(target + 2239);
    auto *t_2240 = buffer.data(target + 2240);
    auto *t_2241 = buffer.data(target + 2241);
    auto *t_2242 = buffer.data(target + 2242);
    auto *t_2243 = buffer.data(target + 2243);
    auto *t_2244 = buffer.data(target + 2244);
    auto *t_2245 = buffer.data(target + 2245);
    auto *t_2246 = buffer.data(target + 2246);
    auto *t_2247 = buffer.data(target + 2247);
    auto *t_2248 = buffer.data(target + 2248);
    auto *t_2249 = buffer.data(target + 2249);
    auto *t_2250 = buffer.data(target + 2250);
    auto *t_2251 = buffer.data(target + 2251);
    auto *t_2252 = buffer.data(target + 2252);
    auto *t_2253 = buffer.data(target + 2253);
    auto *t_2254 = buffer.data(target + 2254);
    auto *t_2255 = buffer.data(target + 2255);
    auto *t_2256 = buffer.data(target + 2256);
    auto *t_2257 = buffer.data(target + 2257);
    auto *t_2258 = buffer.data(target + 2258);
    auto *t_2259 = buffer.data(target + 2259);
    auto *t_2260 = buffer.data(target + 2260);
    auto *t_2261 = buffer.data(target + 2261);
    auto *t_2262 = buffer.data(target + 2262);
    auto *t_2263 = buffer.data(target + 2263);
    auto *t_2264 = buffer.data(target + 2264);
    auto *t_2265 = buffer.data(target + 2265);
    auto *t_2266 = buffer.data(target + 2266);
    auto *t_2267 = buffer.data(target + 2267);
    auto *t_2268 = buffer.data(target + 2268);
    auto *t_2269 = buffer.data(target + 2269);
    auto *t_2270 = buffer.data(target + 2270);
    auto *t_2271 = buffer.data(target + 2271);
    auto *t_2272 = buffer.data(target + 2272);
    auto *t_2273 = buffer.data(target + 2273);
    auto *t_2274 = buffer.data(target + 2274);
    auto *t_2275 = buffer.data(target + 2275);
    auto *t_2276 = buffer.data(target + 2276);
    auto *t_2277 = buffer.data(target + 2277);
    auto *t_2278 = buffer.data(target + 2278);
    auto *t_2279 = buffer.data(target + 2279);
    auto *t_2280 = buffer.data(target + 2280);
    auto *t_2281 = buffer.data(target + 2281);
    auto *t_2282 = buffer.data(target + 2282);
    auto *t_2283 = buffer.data(target + 2283);
    auto *t_2284 = buffer.data(target + 2284);
    auto *t_2285 = buffer.data(target + 2285);
    auto *t_2286 = buffer.data(target + 2286);
    auto *t_2287 = buffer.data(target + 2287);
    auto *t_2288 = buffer.data(target + 2288);
    auto *t_2289 = buffer.data(target + 2289);
    auto *t_2290 = buffer.data(target + 2290);
    auto *t_2291 = buffer.data(target + 2291);
    auto *t_2292 = buffer.data(target + 2292);
    auto *t_2293 = buffer.data(target + 2293);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_1848 = buffer.data(osi0 + 1848);
    const auto *osi0_1849 = buffer.data(osi0 + 1849);
    const auto *osi0_1851 = buffer.data(osi0 + 1851);
    const auto *osi0_1854 = buffer.data(osi0 + 1854);
    const auto *osi0_1858 = buffer.data(osi0 + 1858);
    const auto *osi0_1869 = buffer.data(osi0 + 1869);
    const auto *osi0_1871 = buffer.data(osi0 + 1871);
    const auto *osi0_1872 = buffer.data(osi0 + 1872);
    const auto *osi0_1873 = buffer.data(osi0 + 1873);
    const auto *osi0_2170 = buffer.data(osi0 + 2170);
    const auto *osi0_2177 = buffer.data(osi0 + 2177);
    const auto *osi0_2178 = buffer.data(osi0 + 2178);
    const auto *osi0_2179 = buffer.data(osi0 + 2179);
    const auto *osi0_2180 = buffer.data(osi0 + 2180);
    const auto *osi0_2181 = buffer.data(osi0 + 2181);
    const auto *osi0_2183 = buffer.data(osi0 + 2183);

    const auto *osh_1401 = buffer.data(osh + 1401);
    const auto *osh_1402 = buffer.data(osh + 1402);
    const auto *osh_1403 = buffer.data(osh + 1403);
    const auto *osh_1404 = buffer.data(osh + 1404);
    const auto *osh_1406 = buffer.data(osh + 1406);
    const auto *osh_1422 = buffer.data(osh + 1422);
    const auto *osh_1427 = buffer.data(osh + 1427);
    const auto *osh_1443 = buffer.data(osh + 1443);
    const auto *osh_1445 = buffer.data(osh + 1445);
    const auto *osh_1446 = buffer.data(osh + 1446);
    const auto *osh_1447 = buffer.data(osh + 1447);
    const auto *osh_1448 = buffer.data(osh + 1448);
    const auto *osh_1464 = buffer.data(osh + 1464);
    const auto *osh_1466 = buffer.data(osh + 1466);
    const auto *osh_1467 = buffer.data(osh + 1467);
    const auto *osh_1468 = buffer.data(osh + 1468);
    const auto *osh_1631 = buffer.data(osh + 1631);
    const auto *osh_1632 = buffer.data(osh + 1632);
    const auto *osh_1633 = buffer.data(osh + 1633);
    const auto *osh_1634 = buffer.data(osh + 1634);
    const auto *osh_1635 = buffer.data(osh + 1635);
    const auto *osh_1637 = buffer.data(osh + 1637);

    const auto *osi1_1848 = buffer.data(osi1 + 1848);
    const auto *osi1_1849 = buffer.data(osi1 + 1849);
    const auto *osi1_1851 = buffer.data(osi1 + 1851);
    const auto *osi1_1854 = buffer.data(osi1 + 1854);
    const auto *osi1_1858 = buffer.data(osi1 + 1858);
    const auto *osi1_1869 = buffer.data(osi1 + 1869);
    const auto *osi1_1871 = buffer.data(osi1 + 1871);
    const auto *osi1_1872 = buffer.data(osi1 + 1872);
    const auto *osi1_1873 = buffer.data(osi1 + 1873);
    const auto *osi1_2170 = buffer.data(osi1 + 2170);
    const auto *osi1_2177 = buffer.data(osi1 + 2177);
    const auto *osi1_2178 = buffer.data(osi1 + 2178);
    const auto *osi1_2179 = buffer.data(osi1 + 2179);
    const auto *osi1_2180 = buffer.data(osi1 + 2180);
    const auto *osi1_2181 = buffer.data(osi1 + 2181);
    const auto *osi1_2183 = buffer.data(osi1 + 2183);

    const auto *qsg0_1170 = buffer.data(qsg0 + 1170);
    const auto *qsg0_1171 = buffer.data(qsg0 + 1171);
    const auto *qsg0_1173 = buffer.data(qsg0 + 1173);
    const auto *qsg0_1175 = buffer.data(qsg0 + 1175);
    const auto *qsg0_1176 = buffer.data(qsg0 + 1176);
    const auto *qsg0_1178 = buffer.data(qsg0 + 1178);
    const auto *qsg0_1179 = buffer.data(qsg0 + 1179);
    const auto *qsg0_1180 = buffer.data(qsg0 + 1180);
    const auto *qsg0_1181 = buffer.data(qsg0 + 1181);
    const auto *qsg0_1182 = buffer.data(qsg0 + 1182);
    const auto *qsg0_1183 = buffer.data(qsg0 + 1183);
    const auto *qsg0_1184 = buffer.data(qsg0 + 1184);
    const auto *qsg0_1187 = buffer.data(qsg0 + 1187);
    const auto *qsg0_1189 = buffer.data(qsg0 + 1189);
    const auto *qsg0_1190 = buffer.data(qsg0 + 1190);
    const auto *qsg0_1192 = buffer.data(qsg0 + 1192);
    const auto *qsg0_1193 = buffer.data(qsg0 + 1193);
    const auto *qsg0_1194 = buffer.data(qsg0 + 1194);
    const auto *qsg0_1196 = buffer.data(qsg0 + 1196);
    const auto *qsg0_1197 = buffer.data(qsg0 + 1197);
    const auto *qsg0_1198 = buffer.data(qsg0 + 1198);
    const auto *qsg0_1199 = buffer.data(qsg0 + 1199);
    const auto *qsg0_1200 = buffer.data(qsg0 + 1200);
    const auto *qsg0_1201 = buffer.data(qsg0 + 1201);
    const auto *qsg0_1202 = buffer.data(qsg0 + 1202);
    const auto *qsg0_1203 = buffer.data(qsg0 + 1203);
    const auto *qsg0_1204 = buffer.data(qsg0 + 1204);
    const auto *qsg0_1205 = buffer.data(qsg0 + 1205);
    const auto *qsg0_1206 = buffer.data(qsg0 + 1206);
    const auto *qsg0_1207 = buffer.data(qsg0 + 1207);
    const auto *qsg0_1208 = buffer.data(qsg0 + 1208);
    const auto *qsg0_1209 = buffer.data(qsg0 + 1209);
    const auto *qsg0_1210 = buffer.data(qsg0 + 1210);
    const auto *qsg0_1211 = buffer.data(qsg0 + 1211);
    const auto *qsg0_1212 = buffer.data(qsg0 + 1212);
    const auto *qsg0_1213 = buffer.data(qsg0 + 1213);
    const auto *qsg0_1214 = buffer.data(qsg0 + 1214);
    const auto *qsg0_1215 = buffer.data(qsg0 + 1215);
    const auto *qsg0_1216 = buffer.data(qsg0 + 1216);
    const auto *qsg0_1217 = buffer.data(qsg0 + 1217);
    const auto *qsg0_1218 = buffer.data(qsg0 + 1218);
    const auto *qsg0_1219 = buffer.data(qsg0 + 1219);
    const auto *qsg0_1220 = buffer.data(qsg0 + 1220);
    const auto *qsg0_1221 = buffer.data(qsg0 + 1221);
    const auto *qsg0_1222 = buffer.data(qsg0 + 1222);
    const auto *qsg0_1223 = buffer.data(qsg0 + 1223);
    const auto *qsg0_1224 = buffer.data(qsg0 + 1224);
    const auto *qsg0_1225 = buffer.data(qsg0 + 1225);
    const auto *qsg0_1226 = buffer.data(qsg0 + 1226);
    const auto *qsg0_1227 = buffer.data(qsg0 + 1227);
    const auto *qsg0_1228 = buffer.data(qsg0 + 1228);
    const auto *qsg0_1229 = buffer.data(qsg0 + 1229);

    const auto *qsg1_1170 = buffer.data(qsg1 + 1170);
    const auto *qsg1_1171 = buffer.data(qsg1 + 1171);
    const auto *qsg1_1173 = buffer.data(qsg1 + 1173);
    const auto *qsg1_1175 = buffer.data(qsg1 + 1175);
    const auto *qsg1_1176 = buffer.data(qsg1 + 1176);
    const auto *qsg1_1178 = buffer.data(qsg1 + 1178);
    const auto *qsg1_1179 = buffer.data(qsg1 + 1179);
    const auto *qsg1_1180 = buffer.data(qsg1 + 1180);
    const auto *qsg1_1181 = buffer.data(qsg1 + 1181);
    const auto *qsg1_1182 = buffer.data(qsg1 + 1182);
    const auto *qsg1_1183 = buffer.data(qsg1 + 1183);
    const auto *qsg1_1184 = buffer.data(qsg1 + 1184);
    const auto *qsg1_1187 = buffer.data(qsg1 + 1187);
    const auto *qsg1_1189 = buffer.data(qsg1 + 1189);
    const auto *qsg1_1190 = buffer.data(qsg1 + 1190);
    const auto *qsg1_1192 = buffer.data(qsg1 + 1192);
    const auto *qsg1_1193 = buffer.data(qsg1 + 1193);
    const auto *qsg1_1194 = buffer.data(qsg1 + 1194);
    const auto *qsg1_1196 = buffer.data(qsg1 + 1196);
    const auto *qsg1_1197 = buffer.data(qsg1 + 1197);
    const auto *qsg1_1198 = buffer.data(qsg1 + 1198);
    const auto *qsg1_1199 = buffer.data(qsg1 + 1199);
    const auto *qsg1_1200 = buffer.data(qsg1 + 1200);
    const auto *qsg1_1201 = buffer.data(qsg1 + 1201);
    const auto *qsg1_1202 = buffer.data(qsg1 + 1202);
    const auto *qsg1_1203 = buffer.data(qsg1 + 1203);
    const auto *qsg1_1204 = buffer.data(qsg1 + 1204);
    const auto *qsg1_1205 = buffer.data(qsg1 + 1205);
    const auto *qsg1_1206 = buffer.data(qsg1 + 1206);
    const auto *qsg1_1207 = buffer.data(qsg1 + 1207);
    const auto *qsg1_1208 = buffer.data(qsg1 + 1208);
    const auto *qsg1_1209 = buffer.data(qsg1 + 1209);
    const auto *qsg1_1210 = buffer.data(qsg1 + 1210);
    const auto *qsg1_1211 = buffer.data(qsg1 + 1211);
    const auto *qsg1_1212 = buffer.data(qsg1 + 1212);
    const auto *qsg1_1213 = buffer.data(qsg1 + 1213);
    const auto *qsg1_1214 = buffer.data(qsg1 + 1214);
    const auto *qsg1_1215 = buffer.data(qsg1 + 1215);
    const auto *qsg1_1216 = buffer.data(qsg1 + 1216);
    const auto *qsg1_1217 = buffer.data(qsg1 + 1217);
    const auto *qsg1_1218 = buffer.data(qsg1 + 1218);
    const auto *qsg1_1219 = buffer.data(qsg1 + 1219);
    const auto *qsg1_1220 = buffer.data(qsg1 + 1220);
    const auto *qsg1_1221 = buffer.data(qsg1 + 1221);
    const auto *qsg1_1222 = buffer.data(qsg1 + 1222);
    const auto *qsg1_1223 = buffer.data(qsg1 + 1223);
    const auto *qsg1_1224 = buffer.data(qsg1 + 1224);
    const auto *qsg1_1225 = buffer.data(qsg1 + 1225);
    const auto *qsg1_1226 = buffer.data(qsg1 + 1226);
    const auto *qsg1_1227 = buffer.data(qsg1 + 1227);
    const auto *qsg1_1228 = buffer.data(qsg1 + 1228);
    const auto *qsg1_1229 = buffer.data(qsg1 + 1229);

    const auto *qsh_1631 = buffer.data(qsh + 1631);
    const auto *qsh_1632 = buffer.data(qsh + 1632);
    const auto *qsh_1633 = buffer.data(qsh + 1633);
    const auto *qsh_1634 = buffer.data(qsh + 1634);
    const auto *qsh_1635 = buffer.data(qsh + 1635);
    const auto *qsh_1637 = buffer.data(qsh + 1637);
    const auto *qsh_1638 = buffer.data(qsh + 1638);
    const auto *qsh_1639 = buffer.data(qsh + 1639);
    const auto *qsh_1641 = buffer.data(qsh + 1641);
    const auto *qsh_1643 = buffer.data(qsh + 1643);
    const auto *qsh_1644 = buffer.data(qsh + 1644);
    const auto *qsh_1646 = buffer.data(qsh + 1646);
    const auto *qsh_1647 = buffer.data(qsh + 1647);
    const auto *qsh_1648 = buffer.data(qsh + 1648);
    const auto *qsh_1650 = buffer.data(qsh + 1650);
    const auto *qsh_1651 = buffer.data(qsh + 1651);
    const auto *qsh_1652 = buffer.data(qsh + 1652);
    const auto *qsh_1653 = buffer.data(qsh + 1653);
    const auto *qsh_1654 = buffer.data(qsh + 1654);
    const auto *qsh_1655 = buffer.data(qsh + 1655);
    const auto *qsh_1656 = buffer.data(qsh + 1656);
    const auto *qsh_1657 = buffer.data(qsh + 1657);
    const auto *qsh_1658 = buffer.data(qsh + 1658);
    const auto *qsh_1661 = buffer.data(qsh + 1661);
    const auto *qsh_1663 = buffer.data(qsh + 1663);
    const auto *qsh_1664 = buffer.data(qsh + 1664);
    const auto *qsh_1666 = buffer.data(qsh + 1666);
    const auto *qsh_1667 = buffer.data(qsh + 1667);
    const auto *qsh_1668 = buffer.data(qsh + 1668);
    const auto *qsh_1670 = buffer.data(qsh + 1670);
    const auto *qsh_1671 = buffer.data(qsh + 1671);
    const auto *qsh_1672 = buffer.data(qsh + 1672);
    const auto *qsh_1673 = buffer.data(qsh + 1673);
    const auto *qsh_1674 = buffer.data(qsh + 1674);
    const auto *qsh_1675 = buffer.data(qsh + 1675);
    const auto *qsh_1676 = buffer.data(qsh + 1676);
    const auto *qsh_1677 = buffer.data(qsh + 1677);
    const auto *qsh_1678 = buffer.data(qsh + 1678);
    const auto *qsh_1679 = buffer.data(qsh + 1679);
    const auto *qsh_1680 = buffer.data(qsh + 1680);
    const auto *qsh_1681 = buffer.data(qsh + 1681);
    const auto *qsh_1682 = buffer.data(qsh + 1682);
    const auto *qsh_1683 = buffer.data(qsh + 1683);
    const auto *qsh_1684 = buffer.data(qsh + 1684);
    const auto *qsh_1685 = buffer.data(qsh + 1685);
    const auto *qsh_1686 = buffer.data(qsh + 1686);
    const auto *qsh_1687 = buffer.data(qsh + 1687);
    const auto *qsh_1688 = buffer.data(qsh + 1688);
    const auto *qsh_1689 = buffer.data(qsh + 1689);
    const auto *qsh_1690 = buffer.data(qsh + 1690);
    const auto *qsh_1691 = buffer.data(qsh + 1691);
    const auto *qsh_1692 = buffer.data(qsh + 1692);
    const auto *qsh_1693 = buffer.data(qsh + 1693);
    const auto *qsh_1694 = buffer.data(qsh + 1694);
    const auto *qsh_1695 = buffer.data(qsh + 1695);
    const auto *qsh_1696 = buffer.data(qsh + 1696);
    const auto *qsh_1697 = buffer.data(qsh + 1697);
    const auto *qsh_1698 = buffer.data(qsh + 1698);
    const auto *qsh_1699 = buffer.data(qsh + 1699);
    const auto *qsh_1700 = buffer.data(qsh + 1700);
    const auto *qsh_1701 = buffer.data(qsh + 1701);
    const auto *qsh_1702 = buffer.data(qsh + 1702);
    const auto *qsh_1703 = buffer.data(qsh + 1703);
    const auto *qsh_1704 = buffer.data(qsh + 1704);
    const auto *qsh_1705 = buffer.data(qsh + 1705);
    const auto *qsh_1706 = buffer.data(qsh + 1706);
    const auto *qsh_1707 = buffer.data(qsh + 1707);
    const auto *qsh_1708 = buffer.data(qsh + 1708);
    const auto *qsh_1709 = buffer.data(qsh + 1709);
    const auto *qsh_1710 = buffer.data(qsh + 1710);
    const auto *qsh_1711 = buffer.data(qsh + 1711);
    const auto *qsh_1712 = buffer.data(qsh + 1712);
    const auto *qsh_1713 = buffer.data(qsh + 1713);
    const auto *qsh_1714 = buffer.data(qsh + 1714);
    const auto *qsh_1715 = buffer.data(qsh + 1715);
    const auto *qsh_1716 = buffer.data(qsh + 1716);
    const auto *qsh_1717 = buffer.data(qsh + 1717);
    const auto *qsh_1718 = buffer.data(qsh + 1718);
    const auto *qsh_1719 = buffer.data(qsh + 1719);
    const auto *qsh_1720 = buffer.data(qsh + 1720);
    const auto *qsh_1721 = buffer.data(qsh + 1721);

#pragma omp simd aligned(t_2170, t_2171, t_2172, t_2173, pa_x, pc_x, osi0_2170, osh_1631, \
                         osh_1632, osh_1633, osh_1634, osi1_2170, qsh_1632, qsh_1633, \
                         qsh_1634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2170[k] = pa_x[k] * osi0_2170[k]
                    + f_12 * osh_1631[k]
                    - f_10 * pc_x[k] * osi1_2170[k];

        t_2171[k] = f_11 * osh_1632[k]
                    + f_3 * pc_x[k] * qsh_1632[k];

        t_2172[k] = f_11 * osh_1633[k]
                    + f_3 * pc_x[k] * qsh_1633[k];

        t_2173[k] = f_11 * osh_1634[k]
                    + f_3 * pc_x[k] * qsh_1634[k];
    }

#pragma omp simd aligned(t_2174, t_2175, t_2176, t_2177, pa_x, pc_x, pc_y, osi0_2177, \
                         osh_1635, osh_1637, osi1_2177, qsh_1631, qsh_1635, \
                         qsh_1637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2174[k] = f_11 * osh_1635[k]
                    + f_3 * pc_x[k] * qsh_1635[k];

        t_2175[k] = f_3 * pc_y[k] * qsh_1631[k];

        t_2176[k] = f_11 * osh_1637[k]
                    + f_3 * pc_x[k] * qsh_1637[k];

        t_2177[k] = pa_x[k] * osi0_2177[k]
                    - f_10 * pc_x[k] * osi1_2177[k];
    }

#pragma omp simd aligned(t_2178, t_2179, t_2180, t_2181, pa_x, pc_x, osi0_2178, osi0_2179, \
                         osi0_2180, osi0_2181, osi1_2178, osi1_2179, osi1_2180, \
                         osi1_2181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2178[k] = pa_x[k] * osi0_2178[k]
                    - f_10 * pc_x[k] * osi1_2178[k];

        t_2179[k] = pa_x[k] * osi0_2179[k]
                    - f_10 * pc_x[k] * osi1_2179[k];

        t_2180[k] = pa_x[k] * osi0_2180[k]
                    - f_10 * pc_x[k] * osi1_2180[k];

        t_2181[k] = pa_x[k] * osi0_2181[k]
                    - f_10 * pc_x[k] * osi1_2181[k];
    }

#pragma omp simd aligned(t_2182, t_2183, t_2184, t_2185, pa_x, pc_x, pc_y, osi0_2183, \
                         osi1_2183, qsg0_1170, qsg0_1171, qsg1_1170, qsg1_1171, qsh_1637, \
                         qsh_1638, qsh_1639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2182[k] = f_3 * pc_y[k] * qsh_1637[k];

        t_2183[k] = pa_x[k] * osi0_2183[k]
                    - f_10 * pc_x[k] * osi1_2183[k];

        t_2184[k] = f_1 * qsg0_1170[k]
                    - f_2 * qsg1_1170[k]
                    + f_3 * pc_x[k] * qsh_1638[k];

        t_2185[k] = f_16 * qsg0_1171[k]
                    - f_17 * qsg1_1171[k]
                    + f_3 * pc_x[k] * qsh_1639[k];
    }

#pragma omp simd aligned(t_2186, t_2187, t_2188, t_2189, pc_x, pc_z, qsg0_1173, qsg0_1175, \
                         qsg1_1173, qsg1_1175, qsh_1638, qsh_1639, qsh_1641, \
                         qsh_1643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2186[k] = f_3 * pc_z[k] * qsh_1638[k];

        t_2187[k] = f_8 * qsg0_1173[k]
                    - f_9 * qsg1_1173[k]
                    + f_3 * pc_x[k] * qsh_1641[k];

        t_2188[k] = f_3 * pc_z[k] * qsh_1639[k];

        t_2189[k] = f_8 * qsg0_1175[k]
                    - f_9 * qsg1_1175[k]
                    + f_3 * pc_x[k] * qsh_1643[k];
    }

#pragma omp simd aligned(t_2190, t_2191, t_2192, t_2193, pc_x, pc_z, qsg0_1176, qsg0_1178, \
                         qsg0_1179, qsg1_1176, qsg1_1178, qsg1_1179, qsh_1641, qsh_1644, \
                         qsh_1646, qsh_1647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2190[k] = f_6 * qsg0_1176[k]
                    - f_7 * qsg1_1176[k]
                    + f_3 * pc_x[k] * qsh_1644[k];

        t_2191[k] = f_3 * pc_z[k] * qsh_1641[k];

        t_2192[k] = f_6 * qsg0_1178[k]
                    - f_7 * qsg1_1178[k]
                    + f_3 * pc_x[k] * qsh_1646[k];

        t_2193[k] = f_6 * qsg0_1179[k]
                    - f_7 * qsg1_1179[k]
                    + f_3 * pc_x[k] * qsh_1647[k];
    }

#pragma omp simd aligned(t_2194, t_2195, t_2196, t_2197, pc_x, pc_z, qsg0_1180, qsg0_1182, \
                         qsg0_1183, qsg1_1180, qsg1_1182, qsg1_1183, qsh_1644, qsh_1648, \
                         qsh_1650, qsh_1651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2194[k] = f_4 * qsg0_1180[k]
                    - f_5 * qsg1_1180[k]
                    + f_3 * pc_x[k] * qsh_1648[k];

        t_2195[k] = f_3 * pc_z[k] * qsh_1644[k];

        t_2196[k] = f_4 * qsg0_1182[k]
                    - f_5 * qsg1_1182[k]
                    + f_3 * pc_x[k] * qsh_1650[k];

        t_2197[k] = f_4 * qsg0_1183[k]
                    - f_5 * qsg1_1183[k]
                    + f_3 * pc_x[k] * qsh_1651[k];
    }

#pragma omp simd aligned(t_2198, t_2199, t_2200, t_2201, t_2202, t_2203, pc_x, qsg0_1184, \
                         qsg1_1184, qsh_1652, qsh_1653, qsh_1654, qsh_1655, qsh_1656, \
                         qsh_1657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2198[k] = f_4 * qsg0_1184[k]
                    - f_5 * qsg1_1184[k]
                    + f_3 * pc_x[k] * qsh_1652[k];

        t_2199[k] = f_3 * pc_x[k] * qsh_1653[k];

        t_2200[k] = f_3 * pc_x[k] * qsh_1654[k];

        t_2201[k] = f_3 * pc_x[k] * qsh_1655[k];

        t_2202[k] = f_3 * pc_x[k] * qsh_1656[k];

        t_2203[k] = f_3 * pc_x[k] * qsh_1657[k];
    }

#pragma omp simd aligned(t_2204, t_2205, t_2206, t_2207, pc_x, pc_y, pc_z, osh_1401, \
                         qsg0_1180, qsg1_1180, qsh_1653, qsh_1654, \
                         qsh_1658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2204[k] = f_3 * pc_x[k] * qsh_1658[k];

        t_2205[k] = f_0 * osh_1401[k]
                    + f_1 * qsg0_1180[k]
                    - f_2 * qsg1_1180[k]
                    + f_3 * pc_y[k] * qsh_1653[k];

        t_2206[k] = f_3 * pc_z[k] * qsh_1653[k];

        t_2207[k] = f_4 * qsg0_1180[k]
                    - f_5 * qsg1_1180[k]
                    + f_3 * pc_z[k] * qsh_1654[k];
    }

#pragma omp simd aligned(t_2208, t_2209, t_2210, t_2211, pc_y, pc_z, osh_1406, qsg0_1181, \
                         qsg0_1182, qsg0_1184, qsg1_1181, qsg1_1182, qsg1_1184, qsh_1655, \
                         qsh_1656, qsh_1658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2208[k] = f_6 * qsg0_1181[k]
                    - f_7 * qsg1_1181[k]
                    + f_3 * pc_z[k] * qsh_1655[k];

        t_2209[k] = f_8 * qsg0_1182[k]
                    - f_9 * qsg1_1182[k]
                    + f_3 * pc_z[k] * qsh_1656[k];

        t_2210[k] = f_0 * osh_1406[k]
                    + f_3 * pc_y[k] * qsh_1658[k];

        t_2211[k] = f_1 * qsg0_1184[k]
                    - f_2 * qsg1_1184[k]
                    + f_3 * pc_z[k] * qsh_1658[k];
    }

#pragma omp simd aligned(t_2212, t_2213, t_2214, t_2215, pa_z, pc_x, pc_z, osi0_1848, \
                         osi0_1849, osi0_1851, osi1_1848, osi1_1849, osi1_1851, qsg0_1187, \
                         qsg1_1187, qsh_1661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2212[k] = pa_z[k] * osi0_1848[k]
                    - f_10 * pc_z[k] * osi1_1848[k];

        t_2213[k] = pa_z[k] * osi0_1849[k]
                    - f_10 * pc_z[k] * osi1_1849[k];

        t_2214[k] = f_16 * qsg0_1187[k]
                    - f_17 * qsg1_1187[k]
                    + f_3 * pc_x[k] * qsh_1661[k];

        t_2215[k] = pa_z[k] * osi0_1851[k]
                    - f_10 * pc_z[k] * osi1_1851[k];
    }

#pragma omp simd aligned(t_2216, t_2217, t_2218, pa_z, pc_x, pc_z, osi0_1854, osi1_1854, \
                         qsg0_1189, qsg0_1190, qsg1_1189, qsg1_1190, qsh_1663, \
                         qsh_1664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2216[k] = f_8 * qsg0_1189[k]
                    - f_9 * qsg1_1189[k]
                    + f_3 * pc_x[k] * qsh_1663[k];

        t_2217[k] = f_8 * qsg0_1190[k]
                    - f_9 * qsg1_1190[k]
                    + f_3 * pc_x[k] * qsh_1664[k];

        t_2218[k] = pa_z[k] * osi0_1854[k]
                    - f_10 * pc_z[k] * osi1_1854[k];
    }

#pragma omp simd aligned(t_2219, t_2220, t_2221, pc_x, qsg0_1192, qsg0_1193, qsg0_1194, \
                         qsg1_1192, qsg1_1193, qsg1_1194, qsh_1666, qsh_1667, \
                         qsh_1668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2219[k] = f_6 * qsg0_1192[k]
                    - f_7 * qsg1_1192[k]
                    + f_3 * pc_x[k] * qsh_1666[k];

        t_2220[k] = f_6 * qsg0_1193[k]
                    - f_7 * qsg1_1193[k]
                    + f_3 * pc_x[k] * qsh_1667[k];

        t_2221[k] = f_6 * qsg0_1194[k]
                    - f_7 * qsg1_1194[k]
                    + f_3 * pc_x[k] * qsh_1668[k];
    }

#pragma omp simd aligned(t_2222, t_2223, t_2224, pa_z, pc_x, pc_z, osi0_1858, osi1_1858, \
                         qsg0_1196, qsg0_1197, qsg1_1196, qsg1_1197, qsh_1670, \
                         qsh_1671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2222[k] = pa_z[k] * osi0_1858[k]
                    - f_10 * pc_z[k] * osi1_1858[k];

        t_2223[k] = f_4 * qsg0_1196[k]
                    - f_5 * qsg1_1196[k]
                    + f_3 * pc_x[k] * qsh_1670[k];

        t_2224[k] = f_4 * qsg0_1197[k]
                    - f_5 * qsg1_1197[k]
                    + f_3 * pc_x[k] * qsh_1671[k];
    }

#pragma omp simd aligned(t_2225, t_2226, t_2227, t_2228, t_2229, pc_x, qsg0_1198, qsg0_1199, \
                         qsg1_1198, qsg1_1199, qsh_1672, qsh_1673, qsh_1674, qsh_1675, \
                         qsh_1676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2225[k] = f_4 * qsg0_1198[k]
                    - f_5 * qsg1_1198[k]
                    + f_3 * pc_x[k] * qsh_1672[k];

        t_2226[k] = f_4 * qsg0_1199[k]
                    - f_5 * qsg1_1199[k]
                    + f_3 * pc_x[k] * qsh_1673[k];

        t_2227[k] = f_3 * pc_x[k] * qsh_1674[k];

        t_2228[k] = f_3 * pc_x[k] * qsh_1675[k];

        t_2229[k] = f_3 * pc_x[k] * qsh_1676[k];
    }

#pragma omp simd aligned(t_2230, t_2231, t_2232, t_2233, t_2234, pa_z, pc_x, pc_z, osi0_1869, \
                         osh_1401, osi1_1869, qsh_1674, qsh_1677, qsh_1678, \
                         qsh_1679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2230[k] = f_3 * pc_x[k] * qsh_1677[k];

        t_2231[k] = f_3 * pc_x[k] * qsh_1678[k];

        t_2232[k] = f_3 * pc_x[k] * qsh_1679[k];

        t_2233[k] = pa_z[k] * osi0_1869[k]
                    - f_10 * pc_z[k] * osi1_1869[k];

        t_2234[k] = f_11 * osh_1401[k]
                    + f_3 * pc_z[k] * qsh_1674[k];
    }

#pragma omp simd aligned(t_2235, t_2236, t_2237, pa_z, pc_z, osi0_1871, osi0_1872, osi0_1873, \
                         osh_1402, osh_1403, osh_1404, osi1_1871, osi1_1872, \
                         osi1_1873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2235[k] = pa_z[k] * osi0_1871[k]
                    + f_12 * osh_1402[k]
                    - f_10 * pc_z[k] * osi1_1871[k];

        t_2236[k] = pa_z[k] * osi0_1872[k]
                    + f_13 * osh_1403[k]
                    - f_10 * pc_z[k] * osi1_1872[k];

        t_2237[k] = pa_z[k] * osi0_1873[k]
                    + f_14 * osh_1404[k]
                    - f_10 * pc_z[k] * osi1_1873[k];
    }

#pragma omp simd aligned(t_2238, t_2239, t_2240, pc_x, pc_y, pc_z, osh_1406, osh_1427, \
                         qsg0_1199, qsg0_1200, qsg1_1199, qsg1_1200, qsh_1679, \
                         qsh_1680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2238[k] = f_15 * osh_1427[k]
                    + f_3 * pc_y[k] * qsh_1679[k];

        t_2239[k] = f_11 * osh_1406[k]
                    + f_1 * qsg0_1199[k]
                    - f_2 * qsg1_1199[k]
                    + f_3 * pc_z[k] * qsh_1679[k];

        t_2240[k] = f_1 * qsg0_1200[k]
                    - f_2 * qsg1_1200[k]
                    + f_3 * pc_x[k] * qsh_1680[k];
    }

#pragma omp simd aligned(t_2241, t_2242, t_2243, pc_x, qsg0_1201, qsg0_1202, qsg0_1203, \
                         qsg1_1201, qsg1_1202, qsg1_1203, qsh_1681, qsh_1682, \
                         qsh_1683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2241[k] = f_16 * qsg0_1201[k]
                    - f_17 * qsg1_1201[k]
                    + f_3 * pc_x[k] * qsh_1681[k];

        t_2242[k] = f_16 * qsg0_1202[k]
                    - f_17 * qsg1_1202[k]
                    + f_3 * pc_x[k] * qsh_1682[k];

        t_2243[k] = f_8 * qsg0_1203[k]
                    - f_9 * qsg1_1203[k]
                    + f_3 * pc_x[k] * qsh_1683[k];
    }

#pragma omp simd aligned(t_2244, t_2245, t_2246, pc_x, qsg0_1204, qsg0_1205, qsg0_1206, \
                         qsg1_1204, qsg1_1205, qsg1_1206, qsh_1684, qsh_1685, \
                         qsh_1686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2244[k] = f_8 * qsg0_1204[k]
                    - f_9 * qsg1_1204[k]
                    + f_3 * pc_x[k] * qsh_1684[k];

        t_2245[k] = f_8 * qsg0_1205[k]
                    - f_9 * qsg1_1205[k]
                    + f_3 * pc_x[k] * qsh_1685[k];

        t_2246[k] = f_6 * qsg0_1206[k]
                    - f_7 * qsg1_1206[k]
                    + f_3 * pc_x[k] * qsh_1686[k];
    }

#pragma omp simd aligned(t_2247, t_2248, t_2249, pc_x, qsg0_1207, qsg0_1208, qsg0_1209, \
                         qsg1_1207, qsg1_1208, qsg1_1209, qsh_1687, qsh_1688, \
                         qsh_1689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2247[k] = f_6 * qsg0_1207[k]
                    - f_7 * qsg1_1207[k]
                    + f_3 * pc_x[k] * qsh_1687[k];

        t_2248[k] = f_6 * qsg0_1208[k]
                    - f_7 * qsg1_1208[k]
                    + f_3 * pc_x[k] * qsh_1688[k];

        t_2249[k] = f_6 * qsg0_1209[k]
                    - f_7 * qsg1_1209[k]
                    + f_3 * pc_x[k] * qsh_1689[k];
    }

#pragma omp simd aligned(t_2250, t_2251, t_2252, pc_x, qsg0_1210, qsg0_1211, qsg0_1212, \
                         qsg1_1210, qsg1_1211, qsg1_1212, qsh_1690, qsh_1691, \
                         qsh_1692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2250[k] = f_4 * qsg0_1210[k]
                    - f_5 * qsg1_1210[k]
                    + f_3 * pc_x[k] * qsh_1690[k];

        t_2251[k] = f_4 * qsg0_1211[k]
                    - f_5 * qsg1_1211[k]
                    + f_3 * pc_x[k] * qsh_1691[k];

        t_2252[k] = f_4 * qsg0_1212[k]
                    - f_5 * qsg1_1212[k]
                    + f_3 * pc_x[k] * qsh_1692[k];
    }

#pragma omp simd aligned(t_2253, t_2254, t_2255, t_2256, t_2257, pc_x, qsg0_1213, qsg0_1214, \
                         qsg1_1213, qsg1_1214, qsh_1693, qsh_1694, qsh_1695, qsh_1696, \
                         qsh_1697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2253[k] = f_4 * qsg0_1213[k]
                    - f_5 * qsg1_1213[k]
                    + f_3 * pc_x[k] * qsh_1693[k];

        t_2254[k] = f_4 * qsg0_1214[k]
                    - f_5 * qsg1_1214[k]
                    + f_3 * pc_x[k] * qsh_1694[k];

        t_2255[k] = f_3 * pc_x[k] * qsh_1695[k];

        t_2256[k] = f_3 * pc_x[k] * qsh_1696[k];

        t_2257[k] = f_3 * pc_x[k] * qsh_1697[k];
    }

#pragma omp simd aligned(t_2258, t_2259, t_2260, t_2261, t_2262, pc_x, pc_y, pc_z, osh_1422, \
                         osh_1443, qsg0_1210, qsg1_1210, qsh_1695, qsh_1698, qsh_1699, \
                         qsh_1700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2258[k] = f_3 * pc_x[k] * qsh_1698[k];

        t_2259[k] = f_3 * pc_x[k] * qsh_1699[k];

        t_2260[k] = f_3 * pc_x[k] * qsh_1700[k];

        t_2261[k] = f_18 * osh_1443[k]
                    + f_1 * qsg0_1210[k]
                    - f_2 * qsg1_1210[k]
                    + f_3 * pc_y[k] * qsh_1695[k];

        t_2262[k] = f_12 * osh_1422[k]
                    + f_3 * pc_z[k] * qsh_1695[k];
    }

#pragma omp simd aligned(t_2263, t_2264, t_2265, pc_y, osh_1445, osh_1446, osh_1447, \
                         qsg0_1212, qsg0_1213, qsg0_1214, qsg1_1212, qsg1_1213, qsg1_1214, \
                         qsh_1697, qsh_1698, qsh_1699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2263[k] = f_18 * osh_1445[k]
                    + f_8 * qsg0_1212[k]
                    - f_9 * qsg1_1212[k]
                    + f_3 * pc_y[k] * qsh_1697[k];

        t_2264[k] = f_18 * osh_1446[k]
                    + f_6 * qsg0_1213[k]
                    - f_7 * qsg1_1213[k]
                    + f_3 * pc_y[k] * qsh_1698[k];

        t_2265[k] = f_18 * osh_1447[k]
                    + f_4 * qsg0_1214[k]
                    - f_5 * qsg1_1214[k]
                    + f_3 * pc_y[k] * qsh_1699[k];
    }

#pragma omp simd aligned(t_2266, t_2267, t_2268, pc_x, pc_y, pc_z, osh_1427, osh_1448, \
                         qsg0_1214, qsg0_1215, qsg1_1214, qsg1_1215, qsh_1700, \
                         qsh_1701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2266[k] = f_18 * osh_1448[k]
                    + f_3 * pc_y[k] * qsh_1700[k];

        t_2267[k] = f_12 * osh_1427[k]
                    + f_1 * qsg0_1214[k]
                    - f_2 * qsg1_1214[k]
                    + f_3 * pc_z[k] * qsh_1700[k];

        t_2268[k] = f_1 * qsg0_1215[k]
                    - f_2 * qsg1_1215[k]
                    + f_3 * pc_x[k] * qsh_1701[k];
    }

#pragma omp simd aligned(t_2269, t_2270, t_2271, pc_x, qsg0_1216, qsg0_1217, qsg0_1218, \
                         qsg1_1216, qsg1_1217, qsg1_1218, qsh_1702, qsh_1703, \
                         qsh_1704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2269[k] = f_16 * qsg0_1216[k]
                    - f_17 * qsg1_1216[k]
                    + f_3 * pc_x[k] * qsh_1702[k];

        t_2270[k] = f_16 * qsg0_1217[k]
                    - f_17 * qsg1_1217[k]
                    + f_3 * pc_x[k] * qsh_1703[k];

        t_2271[k] = f_8 * qsg0_1218[k]
                    - f_9 * qsg1_1218[k]
                    + f_3 * pc_x[k] * qsh_1704[k];
    }

#pragma omp simd aligned(t_2272, t_2273, t_2274, pc_x, qsg0_1219, qsg0_1220, qsg0_1221, \
                         qsg1_1219, qsg1_1220, qsg1_1221, qsh_1705, qsh_1706, \
                         qsh_1707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2272[k] = f_8 * qsg0_1219[k]
                    - f_9 * qsg1_1219[k]
                    + f_3 * pc_x[k] * qsh_1705[k];

        t_2273[k] = f_8 * qsg0_1220[k]
                    - f_9 * qsg1_1220[k]
                    + f_3 * pc_x[k] * qsh_1706[k];

        t_2274[k] = f_6 * qsg0_1221[k]
                    - f_7 * qsg1_1221[k]
                    + f_3 * pc_x[k] * qsh_1707[k];
    }

#pragma omp simd aligned(t_2275, t_2276, t_2277, pc_x, qsg0_1222, qsg0_1223, qsg0_1224, \
                         qsg1_1222, qsg1_1223, qsg1_1224, qsh_1708, qsh_1709, \
                         qsh_1710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2275[k] = f_6 * qsg0_1222[k]
                    - f_7 * qsg1_1222[k]
                    + f_3 * pc_x[k] * qsh_1708[k];

        t_2276[k] = f_6 * qsg0_1223[k]
                    - f_7 * qsg1_1223[k]
                    + f_3 * pc_x[k] * qsh_1709[k];

        t_2277[k] = f_6 * qsg0_1224[k]
                    - f_7 * qsg1_1224[k]
                    + f_3 * pc_x[k] * qsh_1710[k];
    }

#pragma omp simd aligned(t_2278, t_2279, t_2280, pc_x, qsg0_1225, qsg0_1226, qsg0_1227, \
                         qsg1_1225, qsg1_1226, qsg1_1227, qsh_1711, qsh_1712, \
                         qsh_1713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2278[k] = f_4 * qsg0_1225[k]
                    - f_5 * qsg1_1225[k]
                    + f_3 * pc_x[k] * qsh_1711[k];

        t_2279[k] = f_4 * qsg0_1226[k]
                    - f_5 * qsg1_1226[k]
                    + f_3 * pc_x[k] * qsh_1712[k];

        t_2280[k] = f_4 * qsg0_1227[k]
                    - f_5 * qsg1_1227[k]
                    + f_3 * pc_x[k] * qsh_1713[k];
    }

#pragma omp simd aligned(t_2281, t_2282, t_2283, t_2284, t_2285, pc_x, qsg0_1228, qsg0_1229, \
                         qsg1_1228, qsg1_1229, qsh_1714, qsh_1715, qsh_1716, qsh_1717, \
                         qsh_1718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2281[k] = f_4 * qsg0_1228[k]
                    - f_5 * qsg1_1228[k]
                    + f_3 * pc_x[k] * qsh_1714[k];

        t_2282[k] = f_4 * qsg0_1229[k]
                    - f_5 * qsg1_1229[k]
                    + f_3 * pc_x[k] * qsh_1715[k];

        t_2283[k] = f_3 * pc_x[k] * qsh_1716[k];

        t_2284[k] = f_3 * pc_x[k] * qsh_1717[k];

        t_2285[k] = f_3 * pc_x[k] * qsh_1718[k];
    }

#pragma omp simd aligned(t_2286, t_2287, t_2288, t_2289, t_2290, pc_x, pc_y, pc_z, osh_1443, \
                         osh_1464, qsg0_1225, qsg1_1225, qsh_1716, qsh_1719, qsh_1720, \
                         qsh_1721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2286[k] = f_3 * pc_x[k] * qsh_1719[k];

        t_2287[k] = f_3 * pc_x[k] * qsh_1720[k];

        t_2288[k] = f_3 * pc_x[k] * qsh_1721[k];

        t_2289[k] = f_19 * osh_1464[k]
                    + f_1 * qsg0_1225[k]
                    - f_2 * qsg1_1225[k]
                    + f_3 * pc_y[k] * qsh_1716[k];

        t_2290[k] = f_13 * osh_1443[k]
                    + f_3 * pc_z[k] * qsh_1716[k];
    }

#pragma omp simd aligned(t_2291, t_2292, t_2293, pc_y, osh_1466, osh_1467, osh_1468, \
                         qsg0_1227, qsg0_1228, qsg0_1229, qsg1_1227, qsg1_1228, qsg1_1229, \
                         qsh_1718, qsh_1719, qsh_1720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2291[k] = f_19 * osh_1466[k]
                    + f_8 * qsg0_1227[k]
                    - f_9 * qsg1_1227[k]
                    + f_3 * pc_y[k] * qsh_1718[k];

        t_2292[k] = f_19 * osh_1467[k]
                    + f_6 * qsg0_1228[k]
                    - f_7 * qsg1_1228[k]
                    + f_3 * pc_y[k] * qsh_1719[k];

        t_2293[k] = f_19 * osh_1468[k]
                    + f_4 * qsg0_1229[k]
                    - f_5 * qsg1_1229[k]
                    + f_3 * pc_y[k] * qsh_1720[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece20(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t osh, const size_t qsg0,
                                                           const size_t qsg1, const size_t qsh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 4.5 / q;
    const auto f_20 = 4.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 2.5 / q;
    const auto f_23 = 3.0 / q;

    auto *t_2294 = buffer.data(target + 2294);
    auto *t_2295 = buffer.data(target + 2295);
    auto *t_2296 = buffer.data(target + 2296);
    auto *t_2297 = buffer.data(target + 2297);
    auto *t_2298 = buffer.data(target + 2298);
    auto *t_2299 = buffer.data(target + 2299);
    auto *t_2300 = buffer.data(target + 2300);
    auto *t_2301 = buffer.data(target + 2301);
    auto *t_2302 = buffer.data(target + 2302);
    auto *t_2303 = buffer.data(target + 2303);
    auto *t_2304 = buffer.data(target + 2304);
    auto *t_2305 = buffer.data(target + 2305);
    auto *t_2306 = buffer.data(target + 2306);
    auto *t_2307 = buffer.data(target + 2307);
    auto *t_2308 = buffer.data(target + 2308);
    auto *t_2309 = buffer.data(target + 2309);
    auto *t_2310 = buffer.data(target + 2310);
    auto *t_2311 = buffer.data(target + 2311);
    auto *t_2312 = buffer.data(target + 2312);
    auto *t_2313 = buffer.data(target + 2313);
    auto *t_2314 = buffer.data(target + 2314);
    auto *t_2315 = buffer.data(target + 2315);
    auto *t_2316 = buffer.data(target + 2316);
    auto *t_2317 = buffer.data(target + 2317);
    auto *t_2318 = buffer.data(target + 2318);
    auto *t_2319 = buffer.data(target + 2319);
    auto *t_2320 = buffer.data(target + 2320);
    auto *t_2321 = buffer.data(target + 2321);
    auto *t_2322 = buffer.data(target + 2322);
    auto *t_2323 = buffer.data(target + 2323);
    auto *t_2324 = buffer.data(target + 2324);
    auto *t_2325 = buffer.data(target + 2325);
    auto *t_2326 = buffer.data(target + 2326);
    auto *t_2327 = buffer.data(target + 2327);
    auto *t_2328 = buffer.data(target + 2328);
    auto *t_2329 = buffer.data(target + 2329);
    auto *t_2330 = buffer.data(target + 2330);
    auto *t_2331 = buffer.data(target + 2331);
    auto *t_2332 = buffer.data(target + 2332);
    auto *t_2333 = buffer.data(target + 2333);
    auto *t_2334 = buffer.data(target + 2334);
    auto *t_2335 = buffer.data(target + 2335);
    auto *t_2336 = buffer.data(target + 2336);
    auto *t_2337 = buffer.data(target + 2337);
    auto *t_2338 = buffer.data(target + 2338);
    auto *t_2339 = buffer.data(target + 2339);
    auto *t_2340 = buffer.data(target + 2340);
    auto *t_2341 = buffer.data(target + 2341);
    auto *t_2342 = buffer.data(target + 2342);
    auto *t_2343 = buffer.data(target + 2343);
    auto *t_2344 = buffer.data(target + 2344);
    auto *t_2345 = buffer.data(target + 2345);
    auto *t_2346 = buffer.data(target + 2346);
    auto *t_2347 = buffer.data(target + 2347);
    auto *t_2348 = buffer.data(target + 2348);
    auto *t_2349 = buffer.data(target + 2349);
    auto *t_2350 = buffer.data(target + 2350);
    auto *t_2351 = buffer.data(target + 2351);
    auto *t_2352 = buffer.data(target + 2352);
    auto *t_2353 = buffer.data(target + 2353);
    auto *t_2354 = buffer.data(target + 2354);
    auto *t_2355 = buffer.data(target + 2355);
    auto *t_2356 = buffer.data(target + 2356);
    auto *t_2357 = buffer.data(target + 2357);
    auto *t_2358 = buffer.data(target + 2358);
    auto *t_2359 = buffer.data(target + 2359);
    auto *t_2360 = buffer.data(target + 2360);
    auto *t_2361 = buffer.data(target + 2361);
    auto *t_2362 = buffer.data(target + 2362);
    auto *t_2363 = buffer.data(target + 2363);
    auto *t_2364 = buffer.data(target + 2364);
    auto *t_2365 = buffer.data(target + 2365);
    auto *t_2366 = buffer.data(target + 2366);
    auto *t_2367 = buffer.data(target + 2367);
    auto *t_2368 = buffer.data(target + 2368);
    auto *t_2369 = buffer.data(target + 2369);
    auto *t_2370 = buffer.data(target + 2370);
    auto *t_2371 = buffer.data(target + 2371);
    auto *t_2372 = buffer.data(target + 2372);
    auto *t_2373 = buffer.data(target + 2373);
    auto *t_2374 = buffer.data(target + 2374);
    auto *t_2375 = buffer.data(target + 2375);
    auto *t_2376 = buffer.data(target + 2376);
    auto *t_2377 = buffer.data(target + 2377);
    auto *t_2378 = buffer.data(target + 2378);
    auto *t_2379 = buffer.data(target + 2379);
    auto *t_2380 = buffer.data(target + 2380);
    auto *t_2381 = buffer.data(target + 2381);
    auto *t_2382 = buffer.data(target + 2382);
    auto *t_2383 = buffer.data(target + 2383);
    auto *t_2384 = buffer.data(target + 2384);
    auto *t_2385 = buffer.data(target + 2385);
    auto *t_2386 = buffer.data(target + 2386);
    auto *t_2387 = buffer.data(target + 2387);
    auto *t_2388 = buffer.data(target + 2388);
    auto *t_2389 = buffer.data(target + 2389);
    auto *t_2390 = buffer.data(target + 2390);
    auto *t_2391 = buffer.data(target + 2391);
    auto *t_2392 = buffer.data(target + 2392);
    auto *t_2393 = buffer.data(target + 2393);
    auto *t_2394 = buffer.data(target + 2394);
    auto *t_2395 = buffer.data(target + 2395);
    auto *t_2396 = buffer.data(target + 2396);
    auto *t_2397 = buffer.data(target + 2397);
    auto *t_2398 = buffer.data(target + 2398);
    auto *t_2399 = buffer.data(target + 2399);
    auto *t_2400 = buffer.data(target + 2400);
    auto *t_2401 = buffer.data(target + 2401);
    auto *t_2402 = buffer.data(target + 2402);
    auto *t_2403 = buffer.data(target + 2403);
    auto *t_2404 = buffer.data(target + 2404);
    auto *t_2405 = buffer.data(target + 2405);
    auto *t_2406 = buffer.data(target + 2406);
    auto *t_2407 = buffer.data(target + 2407);
    auto *t_2408 = buffer.data(target + 2408);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh_1448 = buffer.data(osh + 1448);
    const auto *osh_1464 = buffer.data(osh + 1464);
    const auto *osh_1469 = buffer.data(osh + 1469);
    const auto *osh_1485 = buffer.data(osh + 1485);
    const auto *osh_1487 = buffer.data(osh + 1487);
    const auto *osh_1488 = buffer.data(osh + 1488);
    const auto *osh_1489 = buffer.data(osh + 1489);
    const auto *osh_1490 = buffer.data(osh + 1490);
    const auto *osh_1506 = buffer.data(osh + 1506);
    const auto *osh_1508 = buffer.data(osh + 1508);
    const auto *osh_1509 = buffer.data(osh + 1509);
    const auto *osh_1510 = buffer.data(osh + 1510);
    const auto *osh_1511 = buffer.data(osh + 1511);
    const auto *osh_1527 = buffer.data(osh + 1527);
    const auto *osh_1529 = buffer.data(osh + 1529);
    const auto *osh_1530 = buffer.data(osh + 1530);
    const auto *osh_1531 = buffer.data(osh + 1531);
    const auto *osh_1532 = buffer.data(osh + 1532);
    const auto *osh_1548 = buffer.data(osh + 1548);
    const auto *osh_1550 = buffer.data(osh + 1550);
    const auto *osh_1551 = buffer.data(osh + 1551);
    const auto *osh_1552 = buffer.data(osh + 1552);
    const auto *osh_1553 = buffer.data(osh + 1553);

    const auto *qsg0_1229 = buffer.data(qsg0 + 1229);
    const auto *qsg0_1230 = buffer.data(qsg0 + 1230);
    const auto *qsg0_1231 = buffer.data(qsg0 + 1231);
    const auto *qsg0_1232 = buffer.data(qsg0 + 1232);
    const auto *qsg0_1233 = buffer.data(qsg0 + 1233);
    const auto *qsg0_1234 = buffer.data(qsg0 + 1234);
    const auto *qsg0_1235 = buffer.data(qsg0 + 1235);
    const auto *qsg0_1236 = buffer.data(qsg0 + 1236);
    const auto *qsg0_1237 = buffer.data(qsg0 + 1237);
    const auto *qsg0_1238 = buffer.data(qsg0 + 1238);
    const auto *qsg0_1239 = buffer.data(qsg0 + 1239);
    const auto *qsg0_1240 = buffer.data(qsg0 + 1240);
    const auto *qsg0_1241 = buffer.data(qsg0 + 1241);
    const auto *qsg0_1242 = buffer.data(qsg0 + 1242);
    const auto *qsg0_1243 = buffer.data(qsg0 + 1243);
    const auto *qsg0_1244 = buffer.data(qsg0 + 1244);
    const auto *qsg0_1245 = buffer.data(qsg0 + 1245);
    const auto *qsg0_1246 = buffer.data(qsg0 + 1246);
    const auto *qsg0_1247 = buffer.data(qsg0 + 1247);
    const auto *qsg0_1248 = buffer.data(qsg0 + 1248);
    const auto *qsg0_1249 = buffer.data(qsg0 + 1249);
    const auto *qsg0_1250 = buffer.data(qsg0 + 1250);
    const auto *qsg0_1251 = buffer.data(qsg0 + 1251);
    const auto *qsg0_1252 = buffer.data(qsg0 + 1252);
    const auto *qsg0_1253 = buffer.data(qsg0 + 1253);
    const auto *qsg0_1254 = buffer.data(qsg0 + 1254);
    const auto *qsg0_1255 = buffer.data(qsg0 + 1255);
    const auto *qsg0_1256 = buffer.data(qsg0 + 1256);
    const auto *qsg0_1257 = buffer.data(qsg0 + 1257);
    const auto *qsg0_1258 = buffer.data(qsg0 + 1258);
    const auto *qsg0_1259 = buffer.data(qsg0 + 1259);
    const auto *qsg0_1260 = buffer.data(qsg0 + 1260);
    const auto *qsg0_1261 = buffer.data(qsg0 + 1261);
    const auto *qsg0_1262 = buffer.data(qsg0 + 1262);
    const auto *qsg0_1263 = buffer.data(qsg0 + 1263);
    const auto *qsg0_1264 = buffer.data(qsg0 + 1264);
    const auto *qsg0_1265 = buffer.data(qsg0 + 1265);
    const auto *qsg0_1266 = buffer.data(qsg0 + 1266);
    const auto *qsg0_1267 = buffer.data(qsg0 + 1267);
    const auto *qsg0_1268 = buffer.data(qsg0 + 1268);
    const auto *qsg0_1269 = buffer.data(qsg0 + 1269);
    const auto *qsg0_1270 = buffer.data(qsg0 + 1270);
    const auto *qsg0_1271 = buffer.data(qsg0 + 1271);
    const auto *qsg0_1272 = buffer.data(qsg0 + 1272);
    const auto *qsg0_1273 = buffer.data(qsg0 + 1273);
    const auto *qsg0_1274 = buffer.data(qsg0 + 1274);
    const auto *qsg0_1275 = buffer.data(qsg0 + 1275);
    const auto *qsg0_1276 = buffer.data(qsg0 + 1276);
    const auto *qsg0_1277 = buffer.data(qsg0 + 1277);
    const auto *qsg0_1278 = buffer.data(qsg0 + 1278);
    const auto *qsg0_1279 = buffer.data(qsg0 + 1279);
    const auto *qsg0_1280 = buffer.data(qsg0 + 1280);
    const auto *qsg0_1281 = buffer.data(qsg0 + 1281);
    const auto *qsg0_1282 = buffer.data(qsg0 + 1282);
    const auto *qsg0_1283 = buffer.data(qsg0 + 1283);
    const auto *qsg0_1284 = buffer.data(qsg0 + 1284);
    const auto *qsg0_1285 = buffer.data(qsg0 + 1285);
    const auto *qsg0_1286 = buffer.data(qsg0 + 1286);
    const auto *qsg0_1287 = buffer.data(qsg0 + 1287);
    const auto *qsg0_1288 = buffer.data(qsg0 + 1288);
    const auto *qsg0_1289 = buffer.data(qsg0 + 1289);
    const auto *qsg0_1290 = buffer.data(qsg0 + 1290);

    const auto *qsg1_1229 = buffer.data(qsg1 + 1229);
    const auto *qsg1_1230 = buffer.data(qsg1 + 1230);
    const auto *qsg1_1231 = buffer.data(qsg1 + 1231);
    const auto *qsg1_1232 = buffer.data(qsg1 + 1232);
    const auto *qsg1_1233 = buffer.data(qsg1 + 1233);
    const auto *qsg1_1234 = buffer.data(qsg1 + 1234);
    const auto *qsg1_1235 = buffer.data(qsg1 + 1235);
    const auto *qsg1_1236 = buffer.data(qsg1 + 1236);
    const auto *qsg1_1237 = buffer.data(qsg1 + 1237);
    const auto *qsg1_1238 = buffer.data(qsg1 + 1238);
    const auto *qsg1_1239 = buffer.data(qsg1 + 1239);
    const auto *qsg1_1240 = buffer.data(qsg1 + 1240);
    const auto *qsg1_1241 = buffer.data(qsg1 + 1241);
    const auto *qsg1_1242 = buffer.data(qsg1 + 1242);
    const auto *qsg1_1243 = buffer.data(qsg1 + 1243);
    const auto *qsg1_1244 = buffer.data(qsg1 + 1244);
    const auto *qsg1_1245 = buffer.data(qsg1 + 1245);
    const auto *qsg1_1246 = buffer.data(qsg1 + 1246);
    const auto *qsg1_1247 = buffer.data(qsg1 + 1247);
    const auto *qsg1_1248 = buffer.data(qsg1 + 1248);
    const auto *qsg1_1249 = buffer.data(qsg1 + 1249);
    const auto *qsg1_1250 = buffer.data(qsg1 + 1250);
    const auto *qsg1_1251 = buffer.data(qsg1 + 1251);
    const auto *qsg1_1252 = buffer.data(qsg1 + 1252);
    const auto *qsg1_1253 = buffer.data(qsg1 + 1253);
    const auto *qsg1_1254 = buffer.data(qsg1 + 1254);
    const auto *qsg1_1255 = buffer.data(qsg1 + 1255);
    const auto *qsg1_1256 = buffer.data(qsg1 + 1256);
    const auto *qsg1_1257 = buffer.data(qsg1 + 1257);
    const auto *qsg1_1258 = buffer.data(qsg1 + 1258);
    const auto *qsg1_1259 = buffer.data(qsg1 + 1259);
    const auto *qsg1_1260 = buffer.data(qsg1 + 1260);
    const auto *qsg1_1261 = buffer.data(qsg1 + 1261);
    const auto *qsg1_1262 = buffer.data(qsg1 + 1262);
    const auto *qsg1_1263 = buffer.data(qsg1 + 1263);
    const auto *qsg1_1264 = buffer.data(qsg1 + 1264);
    const auto *qsg1_1265 = buffer.data(qsg1 + 1265);
    const auto *qsg1_1266 = buffer.data(qsg1 + 1266);
    const auto *qsg1_1267 = buffer.data(qsg1 + 1267);
    const auto *qsg1_1268 = buffer.data(qsg1 + 1268);
    const auto *qsg1_1269 = buffer.data(qsg1 + 1269);
    const auto *qsg1_1270 = buffer.data(qsg1 + 1270);
    const auto *qsg1_1271 = buffer.data(qsg1 + 1271);
    const auto *qsg1_1272 = buffer.data(qsg1 + 1272);
    const auto *qsg1_1273 = buffer.data(qsg1 + 1273);
    const auto *qsg1_1274 = buffer.data(qsg1 + 1274);
    const auto *qsg1_1275 = buffer.data(qsg1 + 1275);
    const auto *qsg1_1276 = buffer.data(qsg1 + 1276);
    const auto *qsg1_1277 = buffer.data(qsg1 + 1277);
    const auto *qsg1_1278 = buffer.data(qsg1 + 1278);
    const auto *qsg1_1279 = buffer.data(qsg1 + 1279);
    const auto *qsg1_1280 = buffer.data(qsg1 + 1280);
    const auto *qsg1_1281 = buffer.data(qsg1 + 1281);
    const auto *qsg1_1282 = buffer.data(qsg1 + 1282);
    const auto *qsg1_1283 = buffer.data(qsg1 + 1283);
    const auto *qsg1_1284 = buffer.data(qsg1 + 1284);
    const auto *qsg1_1285 = buffer.data(qsg1 + 1285);
    const auto *qsg1_1286 = buffer.data(qsg1 + 1286);
    const auto *qsg1_1287 = buffer.data(qsg1 + 1287);
    const auto *qsg1_1288 = buffer.data(qsg1 + 1288);
    const auto *qsg1_1289 = buffer.data(qsg1 + 1289);
    const auto *qsg1_1290 = buffer.data(qsg1 + 1290);

    const auto *qsh_1721 = buffer.data(qsh + 1721);
    const auto *qsh_1722 = buffer.data(qsh + 1722);
    const auto *qsh_1723 = buffer.data(qsh + 1723);
    const auto *qsh_1724 = buffer.data(qsh + 1724);
    const auto *qsh_1725 = buffer.data(qsh + 1725);
    const auto *qsh_1726 = buffer.data(qsh + 1726);
    const auto *qsh_1727 = buffer.data(qsh + 1727);
    const auto *qsh_1728 = buffer.data(qsh + 1728);
    const auto *qsh_1729 = buffer.data(qsh + 1729);
    const auto *qsh_1730 = buffer.data(qsh + 1730);
    const auto *qsh_1731 = buffer.data(qsh + 1731);
    const auto *qsh_1732 = buffer.data(qsh + 1732);
    const auto *qsh_1733 = buffer.data(qsh + 1733);
    const auto *qsh_1734 = buffer.data(qsh + 1734);
    const auto *qsh_1735 = buffer.data(qsh + 1735);
    const auto *qsh_1736 = buffer.data(qsh + 1736);
    const auto *qsh_1737 = buffer.data(qsh + 1737);
    const auto *qsh_1738 = buffer.data(qsh + 1738);
    const auto *qsh_1739 = buffer.data(qsh + 1739);
    const auto *qsh_1740 = buffer.data(qsh + 1740);
    const auto *qsh_1741 = buffer.data(qsh + 1741);
    const auto *qsh_1742 = buffer.data(qsh + 1742);
    const auto *qsh_1743 = buffer.data(qsh + 1743);
    const auto *qsh_1744 = buffer.data(qsh + 1744);
    const auto *qsh_1745 = buffer.data(qsh + 1745);
    const auto *qsh_1746 = buffer.data(qsh + 1746);
    const auto *qsh_1747 = buffer.data(qsh + 1747);
    const auto *qsh_1748 = buffer.data(qsh + 1748);
    const auto *qsh_1749 = buffer.data(qsh + 1749);
    const auto *qsh_1750 = buffer.data(qsh + 1750);
    const auto *qsh_1751 = buffer.data(qsh + 1751);
    const auto *qsh_1752 = buffer.data(qsh + 1752);
    const auto *qsh_1753 = buffer.data(qsh + 1753);
    const auto *qsh_1754 = buffer.data(qsh + 1754);
    const auto *qsh_1755 = buffer.data(qsh + 1755);
    const auto *qsh_1756 = buffer.data(qsh + 1756);
    const auto *qsh_1757 = buffer.data(qsh + 1757);
    const auto *qsh_1758 = buffer.data(qsh + 1758);
    const auto *qsh_1759 = buffer.data(qsh + 1759);
    const auto *qsh_1760 = buffer.data(qsh + 1760);
    const auto *qsh_1761 = buffer.data(qsh + 1761);
    const auto *qsh_1762 = buffer.data(qsh + 1762);
    const auto *qsh_1763 = buffer.data(qsh + 1763);
    const auto *qsh_1764 = buffer.data(qsh + 1764);
    const auto *qsh_1765 = buffer.data(qsh + 1765);
    const auto *qsh_1766 = buffer.data(qsh + 1766);
    const auto *qsh_1767 = buffer.data(qsh + 1767);
    const auto *qsh_1768 = buffer.data(qsh + 1768);
    const auto *qsh_1769 = buffer.data(qsh + 1769);
    const auto *qsh_1770 = buffer.data(qsh + 1770);
    const auto *qsh_1771 = buffer.data(qsh + 1771);
    const auto *qsh_1772 = buffer.data(qsh + 1772);
    const auto *qsh_1773 = buffer.data(qsh + 1773);
    const auto *qsh_1774 = buffer.data(qsh + 1774);
    const auto *qsh_1775 = buffer.data(qsh + 1775);
    const auto *qsh_1776 = buffer.data(qsh + 1776);
    const auto *qsh_1777 = buffer.data(qsh + 1777);
    const auto *qsh_1778 = buffer.data(qsh + 1778);
    const auto *qsh_1779 = buffer.data(qsh + 1779);
    const auto *qsh_1780 = buffer.data(qsh + 1780);
    const auto *qsh_1781 = buffer.data(qsh + 1781);
    const auto *qsh_1782 = buffer.data(qsh + 1782);
    const auto *qsh_1783 = buffer.data(qsh + 1783);
    const auto *qsh_1784 = buffer.data(qsh + 1784);
    const auto *qsh_1785 = buffer.data(qsh + 1785);
    const auto *qsh_1786 = buffer.data(qsh + 1786);
    const auto *qsh_1787 = buffer.data(qsh + 1787);
    const auto *qsh_1788 = buffer.data(qsh + 1788);
    const auto *qsh_1789 = buffer.data(qsh + 1789);
    const auto *qsh_1790 = buffer.data(qsh + 1790);
    const auto *qsh_1791 = buffer.data(qsh + 1791);
    const auto *qsh_1792 = buffer.data(qsh + 1792);
    const auto *qsh_1793 = buffer.data(qsh + 1793);
    const auto *qsh_1794 = buffer.data(qsh + 1794);
    const auto *qsh_1795 = buffer.data(qsh + 1795);
    const auto *qsh_1796 = buffer.data(qsh + 1796);
    const auto *qsh_1797 = buffer.data(qsh + 1797);
    const auto *qsh_1798 = buffer.data(qsh + 1798);
    const auto *qsh_1799 = buffer.data(qsh + 1799);
    const auto *qsh_1800 = buffer.data(qsh + 1800);
    const auto *qsh_1801 = buffer.data(qsh + 1801);
    const auto *qsh_1802 = buffer.data(qsh + 1802);
    const auto *qsh_1803 = buffer.data(qsh + 1803);
    const auto *qsh_1804 = buffer.data(qsh + 1804);
    const auto *qsh_1805 = buffer.data(qsh + 1805);
    const auto *qsh_1806 = buffer.data(qsh + 1806);

#pragma omp simd aligned(t_2294, t_2295, t_2296, pc_x, pc_y, pc_z, osh_1448, osh_1469, \
                         qsg0_1229, qsg0_1230, qsg1_1229, qsg1_1230, qsh_1721, \
                         qsh_1722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2294[k] = f_19 * osh_1469[k]
                    + f_3 * pc_y[k] * qsh_1721[k];

        t_2295[k] = f_13 * osh_1448[k]
                    + f_1 * qsg0_1229[k]
                    - f_2 * qsg1_1229[k]
                    + f_3 * pc_z[k] * qsh_1721[k];

        t_2296[k] = f_1 * qsg0_1230[k]
                    - f_2 * qsg1_1230[k]
                    + f_3 * pc_x[k] * qsh_1722[k];
    }

#pragma omp simd aligned(t_2297, t_2298, t_2299, pc_x, qsg0_1231, qsg0_1232, qsg0_1233, \
                         qsg1_1231, qsg1_1232, qsg1_1233, qsh_1723, qsh_1724, \
                         qsh_1725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2297[k] = f_16 * qsg0_1231[k]
                    - f_17 * qsg1_1231[k]
                    + f_3 * pc_x[k] * qsh_1723[k];

        t_2298[k] = f_16 * qsg0_1232[k]
                    - f_17 * qsg1_1232[k]
                    + f_3 * pc_x[k] * qsh_1724[k];

        t_2299[k] = f_8 * qsg0_1233[k]
                    - f_9 * qsg1_1233[k]
                    + f_3 * pc_x[k] * qsh_1725[k];
    }

#pragma omp simd aligned(t_2300, t_2301, t_2302, pc_x, qsg0_1234, qsg0_1235, qsg0_1236, \
                         qsg1_1234, qsg1_1235, qsg1_1236, qsh_1726, qsh_1727, \
                         qsh_1728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2300[k] = f_8 * qsg0_1234[k]
                    - f_9 * qsg1_1234[k]
                    + f_3 * pc_x[k] * qsh_1726[k];

        t_2301[k] = f_8 * qsg0_1235[k]
                    - f_9 * qsg1_1235[k]
                    + f_3 * pc_x[k] * qsh_1727[k];

        t_2302[k] = f_6 * qsg0_1236[k]
                    - f_7 * qsg1_1236[k]
                    + f_3 * pc_x[k] * qsh_1728[k];
    }

#pragma omp simd aligned(t_2303, t_2304, t_2305, pc_x, qsg0_1237, qsg0_1238, qsg0_1239, \
                         qsg1_1237, qsg1_1238, qsg1_1239, qsh_1729, qsh_1730, \
                         qsh_1731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2303[k] = f_6 * qsg0_1237[k]
                    - f_7 * qsg1_1237[k]
                    + f_3 * pc_x[k] * qsh_1729[k];

        t_2304[k] = f_6 * qsg0_1238[k]
                    - f_7 * qsg1_1238[k]
                    + f_3 * pc_x[k] * qsh_1730[k];

        t_2305[k] = f_6 * qsg0_1239[k]
                    - f_7 * qsg1_1239[k]
                    + f_3 * pc_x[k] * qsh_1731[k];
    }

#pragma omp simd aligned(t_2306, t_2307, t_2308, pc_x, qsg0_1240, qsg0_1241, qsg0_1242, \
                         qsg1_1240, qsg1_1241, qsg1_1242, qsh_1732, qsh_1733, \
                         qsh_1734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2306[k] = f_4 * qsg0_1240[k]
                    - f_5 * qsg1_1240[k]
                    + f_3 * pc_x[k] * qsh_1732[k];

        t_2307[k] = f_4 * qsg0_1241[k]
                    - f_5 * qsg1_1241[k]
                    + f_3 * pc_x[k] * qsh_1733[k];

        t_2308[k] = f_4 * qsg0_1242[k]
                    - f_5 * qsg1_1242[k]
                    + f_3 * pc_x[k] * qsh_1734[k];
    }

#pragma omp simd aligned(t_2309, t_2310, t_2311, t_2312, t_2313, pc_x, qsg0_1243, qsg0_1244, \
                         qsg1_1243, qsg1_1244, qsh_1735, qsh_1736, qsh_1737, qsh_1738, \
                         qsh_1739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2309[k] = f_4 * qsg0_1243[k]
                    - f_5 * qsg1_1243[k]
                    + f_3 * pc_x[k] * qsh_1735[k];

        t_2310[k] = f_4 * qsg0_1244[k]
                    - f_5 * qsg1_1244[k]
                    + f_3 * pc_x[k] * qsh_1736[k];

        t_2311[k] = f_3 * pc_x[k] * qsh_1737[k];

        t_2312[k] = f_3 * pc_x[k] * qsh_1738[k];

        t_2313[k] = f_3 * pc_x[k] * qsh_1739[k];
    }

#pragma omp simd aligned(t_2314, t_2315, t_2316, t_2317, t_2318, pc_x, pc_y, pc_z, osh_1464, \
                         osh_1485, qsg0_1240, qsg1_1240, qsh_1737, qsh_1740, qsh_1741, \
                         qsh_1742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2314[k] = f_3 * pc_x[k] * qsh_1740[k];

        t_2315[k] = f_3 * pc_x[k] * qsh_1741[k];

        t_2316[k] = f_3 * pc_x[k] * qsh_1742[k];

        t_2317[k] = f_20 * osh_1485[k]
                    + f_1 * qsg0_1240[k]
                    - f_2 * qsg1_1240[k]
                    + f_3 * pc_y[k] * qsh_1737[k];

        t_2318[k] = f_14 * osh_1464[k]
                    + f_3 * pc_z[k] * qsh_1737[k];
    }

#pragma omp simd aligned(t_2319, t_2320, t_2321, pc_y, osh_1487, osh_1488, osh_1489, \
                         qsg0_1242, qsg0_1243, qsg0_1244, qsg1_1242, qsg1_1243, qsg1_1244, \
                         qsh_1739, qsh_1740, qsh_1741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2319[k] = f_20 * osh_1487[k]
                    + f_8 * qsg0_1242[k]
                    - f_9 * qsg1_1242[k]
                    + f_3 * pc_y[k] * qsh_1739[k];

        t_2320[k] = f_20 * osh_1488[k]
                    + f_6 * qsg0_1243[k]
                    - f_7 * qsg1_1243[k]
                    + f_3 * pc_y[k] * qsh_1740[k];

        t_2321[k] = f_20 * osh_1489[k]
                    + f_4 * qsg0_1244[k]
                    - f_5 * qsg1_1244[k]
                    + f_3 * pc_y[k] * qsh_1741[k];
    }

#pragma omp simd aligned(t_2322, t_2323, t_2324, pc_x, pc_y, pc_z, osh_1469, osh_1490, \
                         qsg0_1244, qsg0_1245, qsg1_1244, qsg1_1245, qsh_1742, \
                         qsh_1743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2322[k] = f_20 * osh_1490[k]
                    + f_3 * pc_y[k] * qsh_1742[k];

        t_2323[k] = f_14 * osh_1469[k]
                    + f_1 * qsg0_1244[k]
                    - f_2 * qsg1_1244[k]
                    + f_3 * pc_z[k] * qsh_1742[k];

        t_2324[k] = f_1 * qsg0_1245[k]
                    - f_2 * qsg1_1245[k]
                    + f_3 * pc_x[k] * qsh_1743[k];
    }

#pragma omp simd aligned(t_2325, t_2326, t_2327, pc_x, qsg0_1246, qsg0_1247, qsg0_1248, \
                         qsg1_1246, qsg1_1247, qsg1_1248, qsh_1744, qsh_1745, \
                         qsh_1746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2325[k] = f_16 * qsg0_1246[k]
                    - f_17 * qsg1_1246[k]
                    + f_3 * pc_x[k] * qsh_1744[k];

        t_2326[k] = f_16 * qsg0_1247[k]
                    - f_17 * qsg1_1247[k]
                    + f_3 * pc_x[k] * qsh_1745[k];

        t_2327[k] = f_8 * qsg0_1248[k]
                    - f_9 * qsg1_1248[k]
                    + f_3 * pc_x[k] * qsh_1746[k];
    }

#pragma omp simd aligned(t_2328, t_2329, t_2330, pc_x, qsg0_1249, qsg0_1250, qsg0_1251, \
                         qsg1_1249, qsg1_1250, qsg1_1251, qsh_1747, qsh_1748, \
                         qsh_1749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2328[k] = f_8 * qsg0_1249[k]
                    - f_9 * qsg1_1249[k]
                    + f_3 * pc_x[k] * qsh_1747[k];

        t_2329[k] = f_8 * qsg0_1250[k]
                    - f_9 * qsg1_1250[k]
                    + f_3 * pc_x[k] * qsh_1748[k];

        t_2330[k] = f_6 * qsg0_1251[k]
                    - f_7 * qsg1_1251[k]
                    + f_3 * pc_x[k] * qsh_1749[k];
    }

#pragma omp simd aligned(t_2331, t_2332, t_2333, pc_x, qsg0_1252, qsg0_1253, qsg0_1254, \
                         qsg1_1252, qsg1_1253, qsg1_1254, qsh_1750, qsh_1751, \
                         qsh_1752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2331[k] = f_6 * qsg0_1252[k]
                    - f_7 * qsg1_1252[k]
                    + f_3 * pc_x[k] * qsh_1750[k];

        t_2332[k] = f_6 * qsg0_1253[k]
                    - f_7 * qsg1_1253[k]
                    + f_3 * pc_x[k] * qsh_1751[k];

        t_2333[k] = f_6 * qsg0_1254[k]
                    - f_7 * qsg1_1254[k]
                    + f_3 * pc_x[k] * qsh_1752[k];
    }

#pragma omp simd aligned(t_2334, t_2335, t_2336, pc_x, qsg0_1255, qsg0_1256, qsg0_1257, \
                         qsg1_1255, qsg1_1256, qsg1_1257, qsh_1753, qsh_1754, \
                         qsh_1755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2334[k] = f_4 * qsg0_1255[k]
                    - f_5 * qsg1_1255[k]
                    + f_3 * pc_x[k] * qsh_1753[k];

        t_2335[k] = f_4 * qsg0_1256[k]
                    - f_5 * qsg1_1256[k]
                    + f_3 * pc_x[k] * qsh_1754[k];

        t_2336[k] = f_4 * qsg0_1257[k]
                    - f_5 * qsg1_1257[k]
                    + f_3 * pc_x[k] * qsh_1755[k];
    }

#pragma omp simd aligned(t_2337, t_2338, t_2339, t_2340, t_2341, pc_x, qsg0_1258, qsg0_1259, \
                         qsg1_1258, qsg1_1259, qsh_1756, qsh_1757, qsh_1758, qsh_1759, \
                         qsh_1760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2337[k] = f_4 * qsg0_1258[k]
                    - f_5 * qsg1_1258[k]
                    + f_3 * pc_x[k] * qsh_1756[k];

        t_2338[k] = f_4 * qsg0_1259[k]
                    - f_5 * qsg1_1259[k]
                    + f_3 * pc_x[k] * qsh_1757[k];

        t_2339[k] = f_3 * pc_x[k] * qsh_1758[k];

        t_2340[k] = f_3 * pc_x[k] * qsh_1759[k];

        t_2341[k] = f_3 * pc_x[k] * qsh_1760[k];
    }

#pragma omp simd aligned(t_2342, t_2343, t_2344, t_2345, t_2346, pc_x, pc_y, pc_z, osh_1485, \
                         osh_1506, qsg0_1255, qsg1_1255, qsh_1758, qsh_1761, qsh_1762, \
                         qsh_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2342[k] = f_3 * pc_x[k] * qsh_1761[k];

        t_2343[k] = f_3 * pc_x[k] * qsh_1762[k];

        t_2344[k] = f_3 * pc_x[k] * qsh_1763[k];

        t_2345[k] = f_21 * osh_1506[k]
                    + f_1 * qsg0_1255[k]
                    - f_2 * qsg1_1255[k]
                    + f_3 * pc_y[k] * qsh_1758[k];

        t_2346[k] = f_22 * osh_1485[k]
                    + f_3 * pc_z[k] * qsh_1758[k];
    }

#pragma omp simd aligned(t_2347, t_2348, t_2349, pc_y, osh_1508, osh_1509, osh_1510, \
                         qsg0_1257, qsg0_1258, qsg0_1259, qsg1_1257, qsg1_1258, qsg1_1259, \
                         qsh_1760, qsh_1761, qsh_1762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2347[k] = f_21 * osh_1508[k]
                    + f_8 * qsg0_1257[k]
                    - f_9 * qsg1_1257[k]
                    + f_3 * pc_y[k] * qsh_1760[k];

        t_2348[k] = f_21 * osh_1509[k]
                    + f_6 * qsg0_1258[k]
                    - f_7 * qsg1_1258[k]
                    + f_3 * pc_y[k] * qsh_1761[k];

        t_2349[k] = f_21 * osh_1510[k]
                    + f_4 * qsg0_1259[k]
                    - f_5 * qsg1_1259[k]
                    + f_3 * pc_y[k] * qsh_1762[k];
    }

#pragma omp simd aligned(t_2350, t_2351, t_2352, pc_x, pc_y, pc_z, osh_1490, osh_1511, \
                         qsg0_1259, qsg0_1260, qsg1_1259, qsg1_1260, qsh_1763, \
                         qsh_1764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2350[k] = f_21 * osh_1511[k]
                    + f_3 * pc_y[k] * qsh_1763[k];

        t_2351[k] = f_22 * osh_1490[k]
                    + f_1 * qsg0_1259[k]
                    - f_2 * qsg1_1259[k]
                    + f_3 * pc_z[k] * qsh_1763[k];

        t_2352[k] = f_1 * qsg0_1260[k]
                    - f_2 * qsg1_1260[k]
                    + f_3 * pc_x[k] * qsh_1764[k];
    }

#pragma omp simd aligned(t_2353, t_2354, t_2355, pc_x, qsg0_1261, qsg0_1262, qsg0_1263, \
                         qsg1_1261, qsg1_1262, qsg1_1263, qsh_1765, qsh_1766, \
                         qsh_1767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2353[k] = f_16 * qsg0_1261[k]
                    - f_17 * qsg1_1261[k]
                    + f_3 * pc_x[k] * qsh_1765[k];

        t_2354[k] = f_16 * qsg0_1262[k]
                    - f_17 * qsg1_1262[k]
                    + f_3 * pc_x[k] * qsh_1766[k];

        t_2355[k] = f_8 * qsg0_1263[k]
                    - f_9 * qsg1_1263[k]
                    + f_3 * pc_x[k] * qsh_1767[k];
    }

#pragma omp simd aligned(t_2356, t_2357, t_2358, pc_x, qsg0_1264, qsg0_1265, qsg0_1266, \
                         qsg1_1264, qsg1_1265, qsg1_1266, qsh_1768, qsh_1769, \
                         qsh_1770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2356[k] = f_8 * qsg0_1264[k]
                    - f_9 * qsg1_1264[k]
                    + f_3 * pc_x[k] * qsh_1768[k];

        t_2357[k] = f_8 * qsg0_1265[k]
                    - f_9 * qsg1_1265[k]
                    + f_3 * pc_x[k] * qsh_1769[k];

        t_2358[k] = f_6 * qsg0_1266[k]
                    - f_7 * qsg1_1266[k]
                    + f_3 * pc_x[k] * qsh_1770[k];
    }

#pragma omp simd aligned(t_2359, t_2360, t_2361, pc_x, qsg0_1267, qsg0_1268, qsg0_1269, \
                         qsg1_1267, qsg1_1268, qsg1_1269, qsh_1771, qsh_1772, \
                         qsh_1773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2359[k] = f_6 * qsg0_1267[k]
                    - f_7 * qsg1_1267[k]
                    + f_3 * pc_x[k] * qsh_1771[k];

        t_2360[k] = f_6 * qsg0_1268[k]
                    - f_7 * qsg1_1268[k]
                    + f_3 * pc_x[k] * qsh_1772[k];

        t_2361[k] = f_6 * qsg0_1269[k]
                    - f_7 * qsg1_1269[k]
                    + f_3 * pc_x[k] * qsh_1773[k];
    }

#pragma omp simd aligned(t_2362, t_2363, t_2364, pc_x, qsg0_1270, qsg0_1271, qsg0_1272, \
                         qsg1_1270, qsg1_1271, qsg1_1272, qsh_1774, qsh_1775, \
                         qsh_1776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2362[k] = f_4 * qsg0_1270[k]
                    - f_5 * qsg1_1270[k]
                    + f_3 * pc_x[k] * qsh_1774[k];

        t_2363[k] = f_4 * qsg0_1271[k]
                    - f_5 * qsg1_1271[k]
                    + f_3 * pc_x[k] * qsh_1775[k];

        t_2364[k] = f_4 * qsg0_1272[k]
                    - f_5 * qsg1_1272[k]
                    + f_3 * pc_x[k] * qsh_1776[k];
    }

#pragma omp simd aligned(t_2365, t_2366, t_2367, t_2368, t_2369, pc_x, qsg0_1273, qsg0_1274, \
                         qsg1_1273, qsg1_1274, qsh_1777, qsh_1778, qsh_1779, qsh_1780, \
                         qsh_1781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2365[k] = f_4 * qsg0_1273[k]
                    - f_5 * qsg1_1273[k]
                    + f_3 * pc_x[k] * qsh_1777[k];

        t_2366[k] = f_4 * qsg0_1274[k]
                    - f_5 * qsg1_1274[k]
                    + f_3 * pc_x[k] * qsh_1778[k];

        t_2367[k] = f_3 * pc_x[k] * qsh_1779[k];

        t_2368[k] = f_3 * pc_x[k] * qsh_1780[k];

        t_2369[k] = f_3 * pc_x[k] * qsh_1781[k];
    }

#pragma omp simd aligned(t_2370, t_2371, t_2372, t_2373, t_2374, pc_x, pc_y, pc_z, osh_1506, \
                         osh_1527, qsg0_1270, qsg1_1270, qsh_1779, qsh_1782, qsh_1783, \
                         qsh_1784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2370[k] = f_3 * pc_x[k] * qsh_1782[k];

        t_2371[k] = f_3 * pc_x[k] * qsh_1783[k];

        t_2372[k] = f_3 * pc_x[k] * qsh_1784[k];

        t_2373[k] = f_23 * osh_1527[k]
                    + f_1 * qsg0_1270[k]
                    - f_2 * qsg1_1270[k]
                    + f_3 * pc_y[k] * qsh_1779[k];

        t_2374[k] = f_23 * osh_1506[k]
                    + f_3 * pc_z[k] * qsh_1779[k];
    }

#pragma omp simd aligned(t_2375, t_2376, t_2377, pc_y, osh_1529, osh_1530, osh_1531, \
                         qsg0_1272, qsg0_1273, qsg0_1274, qsg1_1272, qsg1_1273, qsg1_1274, \
                         qsh_1781, qsh_1782, qsh_1783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2375[k] = f_23 * osh_1529[k]
                    + f_8 * qsg0_1272[k]
                    - f_9 * qsg1_1272[k]
                    + f_3 * pc_y[k] * qsh_1781[k];

        t_2376[k] = f_23 * osh_1530[k]
                    + f_6 * qsg0_1273[k]
                    - f_7 * qsg1_1273[k]
                    + f_3 * pc_y[k] * qsh_1782[k];

        t_2377[k] = f_23 * osh_1531[k]
                    + f_4 * qsg0_1274[k]
                    - f_5 * qsg1_1274[k]
                    + f_3 * pc_y[k] * qsh_1783[k];
    }

#pragma omp simd aligned(t_2378, t_2379, t_2380, pc_x, pc_y, pc_z, osh_1511, osh_1532, \
                         qsg0_1274, qsg0_1275, qsg1_1274, qsg1_1275, qsh_1784, \
                         qsh_1785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2378[k] = f_23 * osh_1532[k]
                    + f_3 * pc_y[k] * qsh_1784[k];

        t_2379[k] = f_23 * osh_1511[k]
                    + f_1 * qsg0_1274[k]
                    - f_2 * qsg1_1274[k]
                    + f_3 * pc_z[k] * qsh_1784[k];

        t_2380[k] = f_1 * qsg0_1275[k]
                    - f_2 * qsg1_1275[k]
                    + f_3 * pc_x[k] * qsh_1785[k];
    }

#pragma omp simd aligned(t_2381, t_2382, t_2383, pc_x, qsg0_1276, qsg0_1277, qsg0_1278, \
                         qsg1_1276, qsg1_1277, qsg1_1278, qsh_1786, qsh_1787, \
                         qsh_1788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2381[k] = f_16 * qsg0_1276[k]
                    - f_17 * qsg1_1276[k]
                    + f_3 * pc_x[k] * qsh_1786[k];

        t_2382[k] = f_16 * qsg0_1277[k]
                    - f_17 * qsg1_1277[k]
                    + f_3 * pc_x[k] * qsh_1787[k];

        t_2383[k] = f_8 * qsg0_1278[k]
                    - f_9 * qsg1_1278[k]
                    + f_3 * pc_x[k] * qsh_1788[k];
    }

#pragma omp simd aligned(t_2384, t_2385, t_2386, pc_x, qsg0_1279, qsg0_1280, qsg0_1281, \
                         qsg1_1279, qsg1_1280, qsg1_1281, qsh_1789, qsh_1790, \
                         qsh_1791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2384[k] = f_8 * qsg0_1279[k]
                    - f_9 * qsg1_1279[k]
                    + f_3 * pc_x[k] * qsh_1789[k];

        t_2385[k] = f_8 * qsg0_1280[k]
                    - f_9 * qsg1_1280[k]
                    + f_3 * pc_x[k] * qsh_1790[k];

        t_2386[k] = f_6 * qsg0_1281[k]
                    - f_7 * qsg1_1281[k]
                    + f_3 * pc_x[k] * qsh_1791[k];
    }

#pragma omp simd aligned(t_2387, t_2388, t_2389, pc_x, qsg0_1282, qsg0_1283, qsg0_1284, \
                         qsg1_1282, qsg1_1283, qsg1_1284, qsh_1792, qsh_1793, \
                         qsh_1794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2387[k] = f_6 * qsg0_1282[k]
                    - f_7 * qsg1_1282[k]
                    + f_3 * pc_x[k] * qsh_1792[k];

        t_2388[k] = f_6 * qsg0_1283[k]
                    - f_7 * qsg1_1283[k]
                    + f_3 * pc_x[k] * qsh_1793[k];

        t_2389[k] = f_6 * qsg0_1284[k]
                    - f_7 * qsg1_1284[k]
                    + f_3 * pc_x[k] * qsh_1794[k];
    }

#pragma omp simd aligned(t_2390, t_2391, t_2392, pc_x, qsg0_1285, qsg0_1286, qsg0_1287, \
                         qsg1_1285, qsg1_1286, qsg1_1287, qsh_1795, qsh_1796, \
                         qsh_1797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2390[k] = f_4 * qsg0_1285[k]
                    - f_5 * qsg1_1285[k]
                    + f_3 * pc_x[k] * qsh_1795[k];

        t_2391[k] = f_4 * qsg0_1286[k]
                    - f_5 * qsg1_1286[k]
                    + f_3 * pc_x[k] * qsh_1796[k];

        t_2392[k] = f_4 * qsg0_1287[k]
                    - f_5 * qsg1_1287[k]
                    + f_3 * pc_x[k] * qsh_1797[k];
    }

#pragma omp simd aligned(t_2393, t_2394, t_2395, t_2396, t_2397, pc_x, qsg0_1288, qsg0_1289, \
                         qsg1_1288, qsg1_1289, qsh_1798, qsh_1799, qsh_1800, qsh_1801, \
                         qsh_1802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2393[k] = f_4 * qsg0_1288[k]
                    - f_5 * qsg1_1288[k]
                    + f_3 * pc_x[k] * qsh_1798[k];

        t_2394[k] = f_4 * qsg0_1289[k]
                    - f_5 * qsg1_1289[k]
                    + f_3 * pc_x[k] * qsh_1799[k];

        t_2395[k] = f_3 * pc_x[k] * qsh_1800[k];

        t_2396[k] = f_3 * pc_x[k] * qsh_1801[k];

        t_2397[k] = f_3 * pc_x[k] * qsh_1802[k];
    }

#pragma omp simd aligned(t_2398, t_2399, t_2400, t_2401, t_2402, pc_x, pc_y, pc_z, osh_1527, \
                         osh_1548, qsg0_1285, qsg1_1285, qsh_1800, qsh_1803, qsh_1804, \
                         qsh_1805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2398[k] = f_3 * pc_x[k] * qsh_1803[k];

        t_2399[k] = f_3 * pc_x[k] * qsh_1804[k];

        t_2400[k] = f_3 * pc_x[k] * qsh_1805[k];

        t_2401[k] = f_22 * osh_1548[k]
                    + f_1 * qsg0_1285[k]
                    - f_2 * qsg1_1285[k]
                    + f_3 * pc_y[k] * qsh_1800[k];

        t_2402[k] = f_21 * osh_1527[k]
                    + f_3 * pc_z[k] * qsh_1800[k];
    }

#pragma omp simd aligned(t_2403, t_2404, t_2405, pc_y, osh_1550, osh_1551, osh_1552, \
                         qsg0_1287, qsg0_1288, qsg0_1289, qsg1_1287, qsg1_1288, qsg1_1289, \
                         qsh_1802, qsh_1803, qsh_1804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2403[k] = f_22 * osh_1550[k]
                    + f_8 * qsg0_1287[k]
                    - f_9 * qsg1_1287[k]
                    + f_3 * pc_y[k] * qsh_1802[k];

        t_2404[k] = f_22 * osh_1551[k]
                    + f_6 * qsg0_1288[k]
                    - f_7 * qsg1_1288[k]
                    + f_3 * pc_y[k] * qsh_1803[k];

        t_2405[k] = f_22 * osh_1552[k]
                    + f_4 * qsg0_1289[k]
                    - f_5 * qsg1_1289[k]
                    + f_3 * pc_y[k] * qsh_1804[k];
    }

#pragma omp simd aligned(t_2406, t_2407, t_2408, pc_x, pc_y, pc_z, osh_1532, osh_1553, \
                         qsg0_1289, qsg0_1290, qsg1_1289, qsg1_1290, qsh_1805, \
                         qsh_1806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2406[k] = f_22 * osh_1553[k]
                    + f_3 * pc_y[k] * qsh_1805[k];

        t_2407[k] = f_21 * osh_1532[k]
                    + f_1 * qsg0_1289[k]
                    - f_2 * qsg1_1289[k]
                    + f_3 * pc_z[k] * qsh_1805[k];

        t_2408[k] = f_1 * qsg0_1290[k]
                    - f_2 * qsg1_1290[k]
                    + f_3 * pc_x[k] * qsh_1806[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece21(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t osi0,
                                                           const size_t osh, const size_t osi1,
                                                           const size_t qsg0, const size_t qsg1,
                                                           const size_t qsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 5.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 5.0 / q;
    const auto f_19 = 4.5 / q;
    const auto f_20 = 4.0 / q;
    const auto f_23 = 3.0 / q;

    auto *t_2409 = buffer.data(target + 2409);
    auto *t_2410 = buffer.data(target + 2410);
    auto *t_2411 = buffer.data(target + 2411);
    auto *t_2412 = buffer.data(target + 2412);
    auto *t_2413 = buffer.data(target + 2413);
    auto *t_2414 = buffer.data(target + 2414);
    auto *t_2415 = buffer.data(target + 2415);
    auto *t_2416 = buffer.data(target + 2416);
    auto *t_2417 = buffer.data(target + 2417);
    auto *t_2418 = buffer.data(target + 2418);
    auto *t_2419 = buffer.data(target + 2419);
    auto *t_2420 = buffer.data(target + 2420);
    auto *t_2421 = buffer.data(target + 2421);
    auto *t_2422 = buffer.data(target + 2422);
    auto *t_2423 = buffer.data(target + 2423);
    auto *t_2424 = buffer.data(target + 2424);
    auto *t_2425 = buffer.data(target + 2425);
    auto *t_2426 = buffer.data(target + 2426);
    auto *t_2427 = buffer.data(target + 2427);
    auto *t_2428 = buffer.data(target + 2428);
    auto *t_2429 = buffer.data(target + 2429);
    auto *t_2430 = buffer.data(target + 2430);
    auto *t_2431 = buffer.data(target + 2431);
    auto *t_2432 = buffer.data(target + 2432);
    auto *t_2433 = buffer.data(target + 2433);
    auto *t_2434 = buffer.data(target + 2434);
    auto *t_2435 = buffer.data(target + 2435);
    auto *t_2436 = buffer.data(target + 2436);
    auto *t_2437 = buffer.data(target + 2437);
    auto *t_2438 = buffer.data(target + 2438);
    auto *t_2439 = buffer.data(target + 2439);
    auto *t_2440 = buffer.data(target + 2440);
    auto *t_2441 = buffer.data(target + 2441);
    auto *t_2442 = buffer.data(target + 2442);
    auto *t_2443 = buffer.data(target + 2443);
    auto *t_2444 = buffer.data(target + 2444);
    auto *t_2445 = buffer.data(target + 2445);
    auto *t_2446 = buffer.data(target + 2446);
    auto *t_2447 = buffer.data(target + 2447);
    auto *t_2448 = buffer.data(target + 2448);
    auto *t_2449 = buffer.data(target + 2449);
    auto *t_2450 = buffer.data(target + 2450);
    auto *t_2451 = buffer.data(target + 2451);
    auto *t_2452 = buffer.data(target + 2452);
    auto *t_2453 = buffer.data(target + 2453);
    auto *t_2454 = buffer.data(target + 2454);
    auto *t_2455 = buffer.data(target + 2455);
    auto *t_2456 = buffer.data(target + 2456);
    auto *t_2457 = buffer.data(target + 2457);
    auto *t_2458 = buffer.data(target + 2458);
    auto *t_2459 = buffer.data(target + 2459);
    auto *t_2460 = buffer.data(target + 2460);
    auto *t_2461 = buffer.data(target + 2461);
    auto *t_2462 = buffer.data(target + 2462);
    auto *t_2463 = buffer.data(target + 2463);
    auto *t_2464 = buffer.data(target + 2464);
    auto *t_2465 = buffer.data(target + 2465);
    auto *t_2466 = buffer.data(target + 2466);
    auto *t_2467 = buffer.data(target + 2467);
    auto *t_2468 = buffer.data(target + 2468);
    auto *t_2469 = buffer.data(target + 2469);
    auto *t_2470 = buffer.data(target + 2470);
    auto *t_2471 = buffer.data(target + 2471);
    auto *t_2472 = buffer.data(target + 2472);
    auto *t_2473 = buffer.data(target + 2473);
    auto *t_2474 = buffer.data(target + 2474);
    auto *t_2475 = buffer.data(target + 2475);
    auto *t_2476 = buffer.data(target + 2476);
    auto *t_2477 = buffer.data(target + 2477);
    auto *t_2478 = buffer.data(target + 2478);
    auto *t_2479 = buffer.data(target + 2479);
    auto *t_2480 = buffer.data(target + 2480);
    auto *t_2481 = buffer.data(target + 2481);
    auto *t_2482 = buffer.data(target + 2482);
    auto *t_2483 = buffer.data(target + 2483);
    auto *t_2484 = buffer.data(target + 2484);
    auto *t_2485 = buffer.data(target + 2485);
    auto *t_2486 = buffer.data(target + 2486);
    auto *t_2487 = buffer.data(target + 2487);
    auto *t_2488 = buffer.data(target + 2488);
    auto *t_2489 = buffer.data(target + 2489);
    auto *t_2490 = buffer.data(target + 2490);
    auto *t_2491 = buffer.data(target + 2491);
    auto *t_2492 = buffer.data(target + 2492);
    auto *t_2493 = buffer.data(target + 2493);
    auto *t_2494 = buffer.data(target + 2494);
    auto *t_2495 = buffer.data(target + 2495);
    auto *t_2496 = buffer.data(target + 2496);
    auto *t_2497 = buffer.data(target + 2497);
    auto *t_2498 = buffer.data(target + 2498);
    auto *t_2499 = buffer.data(target + 2499);
    auto *t_2500 = buffer.data(target + 2500);
    auto *t_2501 = buffer.data(target + 2501);
    auto *t_2502 = buffer.data(target + 2502);
    auto *t_2503 = buffer.data(target + 2503);
    auto *t_2504 = buffer.data(target + 2504);
    auto *t_2505 = buffer.data(target + 2505);
    auto *t_2506 = buffer.data(target + 2506);
    auto *t_2507 = buffer.data(target + 2507);
    auto *t_2508 = buffer.data(target + 2508);
    auto *t_2509 = buffer.data(target + 2509);
    auto *t_2510 = buffer.data(target + 2510);
    auto *t_2511 = buffer.data(target + 2511);
    auto *t_2512 = buffer.data(target + 2512);
    auto *t_2513 = buffer.data(target + 2513);
    auto *t_2514 = buffer.data(target + 2514);
    auto *t_2515 = buffer.data(target + 2515);
    auto *t_2516 = buffer.data(target + 2516);
    auto *t_2517 = buffer.data(target + 2517);
    auto *t_2518 = buffer.data(target + 2518);
    auto *t_2519 = buffer.data(target + 2519);
    auto *t_2520 = buffer.data(target + 2520);
    auto *t_2521 = buffer.data(target + 2521);
    auto *t_2522 = buffer.data(target + 2522);
    auto *t_2523 = buffer.data(target + 2523);
    auto *t_2524 = buffer.data(target + 2524);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osi0_2156 = buffer.data(osi0 + 2156);
    const auto *osi0_2158 = buffer.data(osi0 + 2158);
    const auto *osi0_2161 = buffer.data(osi0 + 2161);
    const auto *osi0_2165 = buffer.data(osi0 + 2165);
    const auto *osi0_2170 = buffer.data(osi0 + 2170);
    const auto *osi0_2177 = buffer.data(osi0 + 2177);
    const auto *osi0_2179 = buffer.data(osi0 + 2179);
    const auto *osi0_2180 = buffer.data(osi0 + 2180);
    const auto *osi0_2181 = buffer.data(osi0 + 2181);
    const auto *osi0_2183 = buffer.data(osi0 + 2183);

    const auto *osh_1548 = buffer.data(osh + 1548);
    const auto *osh_1553 = buffer.data(osh + 1553);
    const auto *osh_1569 = buffer.data(osh + 1569);
    const auto *osh_1571 = buffer.data(osh + 1571);
    const auto *osh_1572 = buffer.data(osh + 1572);
    const auto *osh_1573 = buffer.data(osh + 1573);
    const auto *osh_1574 = buffer.data(osh + 1574);
    const auto *osh_1590 = buffer.data(osh + 1590);
    const auto *osh_1592 = buffer.data(osh + 1592);
    const auto *osh_1593 = buffer.data(osh + 1593);
    const auto *osh_1594 = buffer.data(osh + 1594);
    const auto *osh_1595 = buffer.data(osh + 1595);
    const auto *osh_1611 = buffer.data(osh + 1611);
    const auto *osh_1613 = buffer.data(osh + 1613);
    const auto *osh_1614 = buffer.data(osh + 1614);
    const auto *osh_1615 = buffer.data(osh + 1615);
    const auto *osh_1616 = buffer.data(osh + 1616);
    const auto *osh_1632 = buffer.data(osh + 1632);
    const auto *osh_1634 = buffer.data(osh + 1634);
    const auto *osh_1635 = buffer.data(osh + 1635);
    const auto *osh_1636 = buffer.data(osh + 1636);
    const auto *osh_1637 = buffer.data(osh + 1637);

    const auto *osi1_2156 = buffer.data(osi1 + 2156);
    const auto *osi1_2158 = buffer.data(osi1 + 2158);
    const auto *osi1_2161 = buffer.data(osi1 + 2161);
    const auto *osi1_2165 = buffer.data(osi1 + 2165);
    const auto *osi1_2170 = buffer.data(osi1 + 2170);
    const auto *osi1_2177 = buffer.data(osi1 + 2177);
    const auto *osi1_2179 = buffer.data(osi1 + 2179);
    const auto *osi1_2180 = buffer.data(osi1 + 2180);
    const auto *osi1_2181 = buffer.data(osi1 + 2181);
    const auto *osi1_2183 = buffer.data(osi1 + 2183);

    const auto *qsg0_1291 = buffer.data(qsg0 + 1291);
    const auto *qsg0_1292 = buffer.data(qsg0 + 1292);
    const auto *qsg0_1293 = buffer.data(qsg0 + 1293);
    const auto *qsg0_1294 = buffer.data(qsg0 + 1294);
    const auto *qsg0_1295 = buffer.data(qsg0 + 1295);
    const auto *qsg0_1296 = buffer.data(qsg0 + 1296);
    const auto *qsg0_1297 = buffer.data(qsg0 + 1297);
    const auto *qsg0_1298 = buffer.data(qsg0 + 1298);
    const auto *qsg0_1299 = buffer.data(qsg0 + 1299);
    const auto *qsg0_1300 = buffer.data(qsg0 + 1300);
    const auto *qsg0_1301 = buffer.data(qsg0 + 1301);
    const auto *qsg0_1302 = buffer.data(qsg0 + 1302);
    const auto *qsg0_1303 = buffer.data(qsg0 + 1303);
    const auto *qsg0_1304 = buffer.data(qsg0 + 1304);
    const auto *qsg0_1305 = buffer.data(qsg0 + 1305);
    const auto *qsg0_1306 = buffer.data(qsg0 + 1306);
    const auto *qsg0_1307 = buffer.data(qsg0 + 1307);
    const auto *qsg0_1308 = buffer.data(qsg0 + 1308);
    const auto *qsg0_1309 = buffer.data(qsg0 + 1309);
    const auto *qsg0_1310 = buffer.data(qsg0 + 1310);
    const auto *qsg0_1311 = buffer.data(qsg0 + 1311);
    const auto *qsg0_1312 = buffer.data(qsg0 + 1312);
    const auto *qsg0_1313 = buffer.data(qsg0 + 1313);
    const auto *qsg0_1314 = buffer.data(qsg0 + 1314);
    const auto *qsg0_1315 = buffer.data(qsg0 + 1315);
    const auto *qsg0_1316 = buffer.data(qsg0 + 1316);
    const auto *qsg0_1317 = buffer.data(qsg0 + 1317);
    const auto *qsg0_1318 = buffer.data(qsg0 + 1318);
    const auto *qsg0_1319 = buffer.data(qsg0 + 1319);
    const auto *qsg0_1320 = buffer.data(qsg0 + 1320);
    const auto *qsg0_1321 = buffer.data(qsg0 + 1321);
    const auto *qsg0_1322 = buffer.data(qsg0 + 1322);
    const auto *qsg0_1323 = buffer.data(qsg0 + 1323);
    const auto *qsg0_1324 = buffer.data(qsg0 + 1324);
    const auto *qsg0_1325 = buffer.data(qsg0 + 1325);
    const auto *qsg0_1326 = buffer.data(qsg0 + 1326);
    const auto *qsg0_1327 = buffer.data(qsg0 + 1327);
    const auto *qsg0_1328 = buffer.data(qsg0 + 1328);
    const auto *qsg0_1329 = buffer.data(qsg0 + 1329);
    const auto *qsg0_1330 = buffer.data(qsg0 + 1330);
    const auto *qsg0_1331 = buffer.data(qsg0 + 1331);
    const auto *qsg0_1332 = buffer.data(qsg0 + 1332);
    const auto *qsg0_1333 = buffer.data(qsg0 + 1333);
    const auto *qsg0_1334 = buffer.data(qsg0 + 1334);
    const auto *qsg0_1336 = buffer.data(qsg0 + 1336);
    const auto *qsg0_1338 = buffer.data(qsg0 + 1338);
    const auto *qsg0_1339 = buffer.data(qsg0 + 1339);
    const auto *qsg0_1341 = buffer.data(qsg0 + 1341);
    const auto *qsg0_1342 = buffer.data(qsg0 + 1342);
    const auto *qsg0_1343 = buffer.data(qsg0 + 1343);
    const auto *qsg0_1345 = buffer.data(qsg0 + 1345);
    const auto *qsg0_1346 = buffer.data(qsg0 + 1346);
    const auto *qsg0_1347 = buffer.data(qsg0 + 1347);
    const auto *qsg0_1348 = buffer.data(qsg0 + 1348);
    const auto *qsg0_1350 = buffer.data(qsg0 + 1350);
    const auto *qsg0_1352 = buffer.data(qsg0 + 1352);
    const auto *qsg0_1353 = buffer.data(qsg0 + 1353);

    const auto *qsg1_1291 = buffer.data(qsg1 + 1291);
    const auto *qsg1_1292 = buffer.data(qsg1 + 1292);
    const auto *qsg1_1293 = buffer.data(qsg1 + 1293);
    const auto *qsg1_1294 = buffer.data(qsg1 + 1294);
    const auto *qsg1_1295 = buffer.data(qsg1 + 1295);
    const auto *qsg1_1296 = buffer.data(qsg1 + 1296);
    const auto *qsg1_1297 = buffer.data(qsg1 + 1297);
    const auto *qsg1_1298 = buffer.data(qsg1 + 1298);
    const auto *qsg1_1299 = buffer.data(qsg1 + 1299);
    const auto *qsg1_1300 = buffer.data(qsg1 + 1300);
    const auto *qsg1_1301 = buffer.data(qsg1 + 1301);
    const auto *qsg1_1302 = buffer.data(qsg1 + 1302);
    const auto *qsg1_1303 = buffer.data(qsg1 + 1303);
    const auto *qsg1_1304 = buffer.data(qsg1 + 1304);
    const auto *qsg1_1305 = buffer.data(qsg1 + 1305);
    const auto *qsg1_1306 = buffer.data(qsg1 + 1306);
    const auto *qsg1_1307 = buffer.data(qsg1 + 1307);
    const auto *qsg1_1308 = buffer.data(qsg1 + 1308);
    const auto *qsg1_1309 = buffer.data(qsg1 + 1309);
    const auto *qsg1_1310 = buffer.data(qsg1 + 1310);
    const auto *qsg1_1311 = buffer.data(qsg1 + 1311);
    const auto *qsg1_1312 = buffer.data(qsg1 + 1312);
    const auto *qsg1_1313 = buffer.data(qsg1 + 1313);
    const auto *qsg1_1314 = buffer.data(qsg1 + 1314);
    const auto *qsg1_1315 = buffer.data(qsg1 + 1315);
    const auto *qsg1_1316 = buffer.data(qsg1 + 1316);
    const auto *qsg1_1317 = buffer.data(qsg1 + 1317);
    const auto *qsg1_1318 = buffer.data(qsg1 + 1318);
    const auto *qsg1_1319 = buffer.data(qsg1 + 1319);
    const auto *qsg1_1320 = buffer.data(qsg1 + 1320);
    const auto *qsg1_1321 = buffer.data(qsg1 + 1321);
    const auto *qsg1_1322 = buffer.data(qsg1 + 1322);
    const auto *qsg1_1323 = buffer.data(qsg1 + 1323);
    const auto *qsg1_1324 = buffer.data(qsg1 + 1324);
    const auto *qsg1_1325 = buffer.data(qsg1 + 1325);
    const auto *qsg1_1326 = buffer.data(qsg1 + 1326);
    const auto *qsg1_1327 = buffer.data(qsg1 + 1327);
    const auto *qsg1_1328 = buffer.data(qsg1 + 1328);
    const auto *qsg1_1329 = buffer.data(qsg1 + 1329);
    const auto *qsg1_1330 = buffer.data(qsg1 + 1330);
    const auto *qsg1_1331 = buffer.data(qsg1 + 1331);
    const auto *qsg1_1332 = buffer.data(qsg1 + 1332);
    const auto *qsg1_1333 = buffer.data(qsg1 + 1333);
    const auto *qsg1_1334 = buffer.data(qsg1 + 1334);
    const auto *qsg1_1336 = buffer.data(qsg1 + 1336);
    const auto *qsg1_1338 = buffer.data(qsg1 + 1338);
    const auto *qsg1_1339 = buffer.data(qsg1 + 1339);
    const auto *qsg1_1341 = buffer.data(qsg1 + 1341);
    const auto *qsg1_1342 = buffer.data(qsg1 + 1342);
    const auto *qsg1_1343 = buffer.data(qsg1 + 1343);
    const auto *qsg1_1345 = buffer.data(qsg1 + 1345);
    const auto *qsg1_1346 = buffer.data(qsg1 + 1346);
    const auto *qsg1_1347 = buffer.data(qsg1 + 1347);
    const auto *qsg1_1348 = buffer.data(qsg1 + 1348);
    const auto *qsg1_1350 = buffer.data(qsg1 + 1350);
    const auto *qsg1_1352 = buffer.data(qsg1 + 1352);
    const auto *qsg1_1353 = buffer.data(qsg1 + 1353);

    const auto *qsh_1807 = buffer.data(qsh + 1807);
    const auto *qsh_1808 = buffer.data(qsh + 1808);
    const auto *qsh_1809 = buffer.data(qsh + 1809);
    const auto *qsh_1810 = buffer.data(qsh + 1810);
    const auto *qsh_1811 = buffer.data(qsh + 1811);
    const auto *qsh_1812 = buffer.data(qsh + 1812);
    const auto *qsh_1813 = buffer.data(qsh + 1813);
    const auto *qsh_1814 = buffer.data(qsh + 1814);
    const auto *qsh_1815 = buffer.data(qsh + 1815);
    const auto *qsh_1816 = buffer.data(qsh + 1816);
    const auto *qsh_1817 = buffer.data(qsh + 1817);
    const auto *qsh_1818 = buffer.data(qsh + 1818);
    const auto *qsh_1819 = buffer.data(qsh + 1819);
    const auto *qsh_1820 = buffer.data(qsh + 1820);
    const auto *qsh_1821 = buffer.data(qsh + 1821);
    const auto *qsh_1822 = buffer.data(qsh + 1822);
    const auto *qsh_1823 = buffer.data(qsh + 1823);
    const auto *qsh_1824 = buffer.data(qsh + 1824);
    const auto *qsh_1825 = buffer.data(qsh + 1825);
    const auto *qsh_1826 = buffer.data(qsh + 1826);
    const auto *qsh_1827 = buffer.data(qsh + 1827);
    const auto *qsh_1828 = buffer.data(qsh + 1828);
    const auto *qsh_1829 = buffer.data(qsh + 1829);
    const auto *qsh_1830 = buffer.data(qsh + 1830);
    const auto *qsh_1831 = buffer.data(qsh + 1831);
    const auto *qsh_1832 = buffer.data(qsh + 1832);
    const auto *qsh_1833 = buffer.data(qsh + 1833);
    const auto *qsh_1834 = buffer.data(qsh + 1834);
    const auto *qsh_1835 = buffer.data(qsh + 1835);
    const auto *qsh_1836 = buffer.data(qsh + 1836);
    const auto *qsh_1837 = buffer.data(qsh + 1837);
    const auto *qsh_1838 = buffer.data(qsh + 1838);
    const auto *qsh_1839 = buffer.data(qsh + 1839);
    const auto *qsh_1840 = buffer.data(qsh + 1840);
    const auto *qsh_1841 = buffer.data(qsh + 1841);
    const auto *qsh_1842 = buffer.data(qsh + 1842);
    const auto *qsh_1843 = buffer.data(qsh + 1843);
    const auto *qsh_1844 = buffer.data(qsh + 1844);
    const auto *qsh_1845 = buffer.data(qsh + 1845);
    const auto *qsh_1846 = buffer.data(qsh + 1846);
    const auto *qsh_1847 = buffer.data(qsh + 1847);
    const auto *qsh_1848 = buffer.data(qsh + 1848);
    const auto *qsh_1849 = buffer.data(qsh + 1849);
    const auto *qsh_1850 = buffer.data(qsh + 1850);
    const auto *qsh_1851 = buffer.data(qsh + 1851);
    const auto *qsh_1852 = buffer.data(qsh + 1852);
    const auto *qsh_1853 = buffer.data(qsh + 1853);
    const auto *qsh_1854 = buffer.data(qsh + 1854);
    const auto *qsh_1855 = buffer.data(qsh + 1855);
    const auto *qsh_1856 = buffer.data(qsh + 1856);
    const auto *qsh_1857 = buffer.data(qsh + 1857);
    const auto *qsh_1858 = buffer.data(qsh + 1858);
    const auto *qsh_1859 = buffer.data(qsh + 1859);
    const auto *qsh_1860 = buffer.data(qsh + 1860);
    const auto *qsh_1861 = buffer.data(qsh + 1861);
    const auto *qsh_1862 = buffer.data(qsh + 1862);
    const auto *qsh_1863 = buffer.data(qsh + 1863);
    const auto *qsh_1864 = buffer.data(qsh + 1864);
    const auto *qsh_1865 = buffer.data(qsh + 1865);
    const auto *qsh_1866 = buffer.data(qsh + 1866);
    const auto *qsh_1867 = buffer.data(qsh + 1867);
    const auto *qsh_1868 = buffer.data(qsh + 1868);
    const auto *qsh_1870 = buffer.data(qsh + 1870);
    const auto *qsh_1872 = buffer.data(qsh + 1872);
    const auto *qsh_1873 = buffer.data(qsh + 1873);
    const auto *qsh_1875 = buffer.data(qsh + 1875);
    const auto *qsh_1876 = buffer.data(qsh + 1876);
    const auto *qsh_1877 = buffer.data(qsh + 1877);
    const auto *qsh_1879 = buffer.data(qsh + 1879);
    const auto *qsh_1880 = buffer.data(qsh + 1880);
    const auto *qsh_1881 = buffer.data(qsh + 1881);
    const auto *qsh_1882 = buffer.data(qsh + 1882);
    const auto *qsh_1884 = buffer.data(qsh + 1884);
    const auto *qsh_1885 = buffer.data(qsh + 1885);
    const auto *qsh_1886 = buffer.data(qsh + 1886);
    const auto *qsh_1887 = buffer.data(qsh + 1887);
    const auto *qsh_1888 = buffer.data(qsh + 1888);
    const auto *qsh_1889 = buffer.data(qsh + 1889);
    const auto *qsh_1890 = buffer.data(qsh + 1890);
    const auto *qsh_1892 = buffer.data(qsh + 1892);
    const auto *qsh_1893 = buffer.data(qsh + 1893);

#pragma omp simd aligned(t_2409, t_2410, t_2411, pc_x, qsg0_1291, qsg0_1292, qsg0_1293, \
                         qsg1_1291, qsg1_1292, qsg1_1293, qsh_1807, qsh_1808, \
                         qsh_1809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2409[k] = f_16 * qsg0_1291[k]
                    - f_17 * qsg1_1291[k]
                    + f_3 * pc_x[k] * qsh_1807[k];

        t_2410[k] = f_16 * qsg0_1292[k]
                    - f_17 * qsg1_1292[k]
                    + f_3 * pc_x[k] * qsh_1808[k];

        t_2411[k] = f_8 * qsg0_1293[k]
                    - f_9 * qsg1_1293[k]
                    + f_3 * pc_x[k] * qsh_1809[k];
    }

#pragma omp simd aligned(t_2412, t_2413, t_2414, pc_x, qsg0_1294, qsg0_1295, qsg0_1296, \
                         qsg1_1294, qsg1_1295, qsg1_1296, qsh_1810, qsh_1811, \
                         qsh_1812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2412[k] = f_8 * qsg0_1294[k]
                    - f_9 * qsg1_1294[k]
                    + f_3 * pc_x[k] * qsh_1810[k];

        t_2413[k] = f_8 * qsg0_1295[k]
                    - f_9 * qsg1_1295[k]
                    + f_3 * pc_x[k] * qsh_1811[k];

        t_2414[k] = f_6 * qsg0_1296[k]
                    - f_7 * qsg1_1296[k]
                    + f_3 * pc_x[k] * qsh_1812[k];
    }

#pragma omp simd aligned(t_2415, t_2416, t_2417, pc_x, qsg0_1297, qsg0_1298, qsg0_1299, \
                         qsg1_1297, qsg1_1298, qsg1_1299, qsh_1813, qsh_1814, \
                         qsh_1815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2415[k] = f_6 * qsg0_1297[k]
                    - f_7 * qsg1_1297[k]
                    + f_3 * pc_x[k] * qsh_1813[k];

        t_2416[k] = f_6 * qsg0_1298[k]
                    - f_7 * qsg1_1298[k]
                    + f_3 * pc_x[k] * qsh_1814[k];

        t_2417[k] = f_6 * qsg0_1299[k]
                    - f_7 * qsg1_1299[k]
                    + f_3 * pc_x[k] * qsh_1815[k];
    }

#pragma omp simd aligned(t_2418, t_2419, t_2420, pc_x, qsg0_1300, qsg0_1301, qsg0_1302, \
                         qsg1_1300, qsg1_1301, qsg1_1302, qsh_1816, qsh_1817, \
                         qsh_1818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2418[k] = f_4 * qsg0_1300[k]
                    - f_5 * qsg1_1300[k]
                    + f_3 * pc_x[k] * qsh_1816[k];

        t_2419[k] = f_4 * qsg0_1301[k]
                    - f_5 * qsg1_1301[k]
                    + f_3 * pc_x[k] * qsh_1817[k];

        t_2420[k] = f_4 * qsg0_1302[k]
                    - f_5 * qsg1_1302[k]
                    + f_3 * pc_x[k] * qsh_1818[k];
    }

#pragma omp simd aligned(t_2421, t_2422, t_2423, t_2424, t_2425, pc_x, qsg0_1303, qsg0_1304, \
                         qsg1_1303, qsg1_1304, qsh_1819, qsh_1820, qsh_1821, qsh_1822, \
                         qsh_1823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2421[k] = f_4 * qsg0_1303[k]
                    - f_5 * qsg1_1303[k]
                    + f_3 * pc_x[k] * qsh_1819[k];

        t_2422[k] = f_4 * qsg0_1304[k]
                    - f_5 * qsg1_1304[k]
                    + f_3 * pc_x[k] * qsh_1820[k];

        t_2423[k] = f_3 * pc_x[k] * qsh_1821[k];

        t_2424[k] = f_3 * pc_x[k] * qsh_1822[k];

        t_2425[k] = f_3 * pc_x[k] * qsh_1823[k];
    }

#pragma omp simd aligned(t_2426, t_2427, t_2428, t_2429, t_2430, pc_x, pc_y, pc_z, osh_1548, \
                         osh_1569, qsg0_1300, qsg1_1300, qsh_1821, qsh_1824, qsh_1825, \
                         qsh_1826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2426[k] = f_3 * pc_x[k] * qsh_1824[k];

        t_2427[k] = f_3 * pc_x[k] * qsh_1825[k];

        t_2428[k] = f_3 * pc_x[k] * qsh_1826[k];

        t_2429[k] = f_14 * osh_1569[k]
                    + f_1 * qsg0_1300[k]
                    - f_2 * qsg1_1300[k]
                    + f_3 * pc_y[k] * qsh_1821[k];

        t_2430[k] = f_20 * osh_1548[k]
                    + f_3 * pc_z[k] * qsh_1821[k];
    }

#pragma omp simd aligned(t_2431, t_2432, t_2433, pc_y, osh_1571, osh_1572, osh_1573, \
                         qsg0_1302, qsg0_1303, qsg0_1304, qsg1_1302, qsg1_1303, qsg1_1304, \
                         qsh_1823, qsh_1824, qsh_1825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2431[k] = f_14 * osh_1571[k]
                    + f_8 * qsg0_1302[k]
                    - f_9 * qsg1_1302[k]
                    + f_3 * pc_y[k] * qsh_1823[k];

        t_2432[k] = f_14 * osh_1572[k]
                    + f_6 * qsg0_1303[k]
                    - f_7 * qsg1_1303[k]
                    + f_3 * pc_y[k] * qsh_1824[k];

        t_2433[k] = f_14 * osh_1573[k]
                    + f_4 * qsg0_1304[k]
                    - f_5 * qsg1_1304[k]
                    + f_3 * pc_y[k] * qsh_1825[k];
    }

#pragma omp simd aligned(t_2434, t_2435, t_2436, pc_x, pc_y, pc_z, osh_1553, osh_1574, \
                         qsg0_1304, qsg0_1305, qsg1_1304, qsg1_1305, qsh_1826, \
                         qsh_1827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2434[k] = f_14 * osh_1574[k]
                    + f_3 * pc_y[k] * qsh_1826[k];

        t_2435[k] = f_20 * osh_1553[k]
                    + f_1 * qsg0_1304[k]
                    - f_2 * qsg1_1304[k]
                    + f_3 * pc_z[k] * qsh_1826[k];

        t_2436[k] = f_1 * qsg0_1305[k]
                    - f_2 * qsg1_1305[k]
                    + f_3 * pc_x[k] * qsh_1827[k];
    }

#pragma omp simd aligned(t_2437, t_2438, t_2439, pc_x, qsg0_1306, qsg0_1307, qsg0_1308, \
                         qsg1_1306, qsg1_1307, qsg1_1308, qsh_1828, qsh_1829, \
                         qsh_1830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2437[k] = f_16 * qsg0_1306[k]
                    - f_17 * qsg1_1306[k]
                    + f_3 * pc_x[k] * qsh_1828[k];

        t_2438[k] = f_16 * qsg0_1307[k]
                    - f_17 * qsg1_1307[k]
                    + f_3 * pc_x[k] * qsh_1829[k];

        t_2439[k] = f_8 * qsg0_1308[k]
                    - f_9 * qsg1_1308[k]
                    + f_3 * pc_x[k] * qsh_1830[k];
    }

#pragma omp simd aligned(t_2440, t_2441, t_2442, pc_x, qsg0_1309, qsg0_1310, qsg0_1311, \
                         qsg1_1309, qsg1_1310, qsg1_1311, qsh_1831, qsh_1832, \
                         qsh_1833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2440[k] = f_8 * qsg0_1309[k]
                    - f_9 * qsg1_1309[k]
                    + f_3 * pc_x[k] * qsh_1831[k];

        t_2441[k] = f_8 * qsg0_1310[k]
                    - f_9 * qsg1_1310[k]
                    + f_3 * pc_x[k] * qsh_1832[k];

        t_2442[k] = f_6 * qsg0_1311[k]
                    - f_7 * qsg1_1311[k]
                    + f_3 * pc_x[k] * qsh_1833[k];
    }

#pragma omp simd aligned(t_2443, t_2444, t_2445, pc_x, qsg0_1312, qsg0_1313, qsg0_1314, \
                         qsg1_1312, qsg1_1313, qsg1_1314, qsh_1834, qsh_1835, \
                         qsh_1836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2443[k] = f_6 * qsg0_1312[k]
                    - f_7 * qsg1_1312[k]
                    + f_3 * pc_x[k] * qsh_1834[k];

        t_2444[k] = f_6 * qsg0_1313[k]
                    - f_7 * qsg1_1313[k]
                    + f_3 * pc_x[k] * qsh_1835[k];

        t_2445[k] = f_6 * qsg0_1314[k]
                    - f_7 * qsg1_1314[k]
                    + f_3 * pc_x[k] * qsh_1836[k];
    }

#pragma omp simd aligned(t_2446, t_2447, t_2448, pc_x, qsg0_1315, qsg0_1316, qsg0_1317, \
                         qsg1_1315, qsg1_1316, qsg1_1317, qsh_1837, qsh_1838, \
                         qsh_1839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2446[k] = f_4 * qsg0_1315[k]
                    - f_5 * qsg1_1315[k]
                    + f_3 * pc_x[k] * qsh_1837[k];

        t_2447[k] = f_4 * qsg0_1316[k]
                    - f_5 * qsg1_1316[k]
                    + f_3 * pc_x[k] * qsh_1838[k];

        t_2448[k] = f_4 * qsg0_1317[k]
                    - f_5 * qsg1_1317[k]
                    + f_3 * pc_x[k] * qsh_1839[k];
    }

#pragma omp simd aligned(t_2449, t_2450, t_2451, t_2452, t_2453, pc_x, qsg0_1318, qsg0_1319, \
                         qsg1_1318, qsg1_1319, qsh_1840, qsh_1841, qsh_1842, qsh_1843, \
                         qsh_1844 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2449[k] = f_4 * qsg0_1318[k]
                    - f_5 * qsg1_1318[k]
                    + f_3 * pc_x[k] * qsh_1840[k];

        t_2450[k] = f_4 * qsg0_1319[k]
                    - f_5 * qsg1_1319[k]
                    + f_3 * pc_x[k] * qsh_1841[k];

        t_2451[k] = f_3 * pc_x[k] * qsh_1842[k];

        t_2452[k] = f_3 * pc_x[k] * qsh_1843[k];

        t_2453[k] = f_3 * pc_x[k] * qsh_1844[k];
    }

#pragma omp simd aligned(t_2454, t_2455, t_2456, t_2457, t_2458, pc_x, pc_y, pc_z, osh_1569, \
                         osh_1590, qsg0_1315, qsg1_1315, qsh_1842, qsh_1845, qsh_1846, \
                         qsh_1847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2454[k] = f_3 * pc_x[k] * qsh_1845[k];

        t_2455[k] = f_3 * pc_x[k] * qsh_1846[k];

        t_2456[k] = f_3 * pc_x[k] * qsh_1847[k];

        t_2457[k] = f_13 * osh_1590[k]
                    + f_1 * qsg0_1315[k]
                    - f_2 * qsg1_1315[k]
                    + f_3 * pc_y[k] * qsh_1842[k];

        t_2458[k] = f_19 * osh_1569[k]
                    + f_3 * pc_z[k] * qsh_1842[k];
    }

#pragma omp simd aligned(t_2459, t_2460, t_2461, pc_y, osh_1592, osh_1593, osh_1594, \
                         qsg0_1317, qsg0_1318, qsg0_1319, qsg1_1317, qsg1_1318, qsg1_1319, \
                         qsh_1844, qsh_1845, qsh_1846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2459[k] = f_13 * osh_1592[k]
                    + f_8 * qsg0_1317[k]
                    - f_9 * qsg1_1317[k]
                    + f_3 * pc_y[k] * qsh_1844[k];

        t_2460[k] = f_13 * osh_1593[k]
                    + f_6 * qsg0_1318[k]
                    - f_7 * qsg1_1318[k]
                    + f_3 * pc_y[k] * qsh_1845[k];

        t_2461[k] = f_13 * osh_1594[k]
                    + f_4 * qsg0_1319[k]
                    - f_5 * qsg1_1319[k]
                    + f_3 * pc_y[k] * qsh_1846[k];
    }

#pragma omp simd aligned(t_2462, t_2463, t_2464, pc_x, pc_y, pc_z, osh_1574, osh_1595, \
                         qsg0_1319, qsg0_1320, qsg1_1319, qsg1_1320, qsh_1847, \
                         qsh_1848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2462[k] = f_13 * osh_1595[k]
                    + f_3 * pc_y[k] * qsh_1847[k];

        t_2463[k] = f_19 * osh_1574[k]
                    + f_1 * qsg0_1319[k]
                    - f_2 * qsg1_1319[k]
                    + f_3 * pc_z[k] * qsh_1847[k];

        t_2464[k] = f_1 * qsg0_1320[k]
                    - f_2 * qsg1_1320[k]
                    + f_3 * pc_x[k] * qsh_1848[k];
    }

#pragma omp simd aligned(t_2465, t_2466, t_2467, pc_x, qsg0_1321, qsg0_1322, qsg0_1323, \
                         qsg1_1321, qsg1_1322, qsg1_1323, qsh_1849, qsh_1850, \
                         qsh_1851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2465[k] = f_16 * qsg0_1321[k]
                    - f_17 * qsg1_1321[k]
                    + f_3 * pc_x[k] * qsh_1849[k];

        t_2466[k] = f_16 * qsg0_1322[k]
                    - f_17 * qsg1_1322[k]
                    + f_3 * pc_x[k] * qsh_1850[k];

        t_2467[k] = f_8 * qsg0_1323[k]
                    - f_9 * qsg1_1323[k]
                    + f_3 * pc_x[k] * qsh_1851[k];
    }

#pragma omp simd aligned(t_2468, t_2469, t_2470, pc_x, qsg0_1324, qsg0_1325, qsg0_1326, \
                         qsg1_1324, qsg1_1325, qsg1_1326, qsh_1852, qsh_1853, \
                         qsh_1854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2468[k] = f_8 * qsg0_1324[k]
                    - f_9 * qsg1_1324[k]
                    + f_3 * pc_x[k] * qsh_1852[k];

        t_2469[k] = f_8 * qsg0_1325[k]
                    - f_9 * qsg1_1325[k]
                    + f_3 * pc_x[k] * qsh_1853[k];

        t_2470[k] = f_6 * qsg0_1326[k]
                    - f_7 * qsg1_1326[k]
                    + f_3 * pc_x[k] * qsh_1854[k];
    }

#pragma omp simd aligned(t_2471, t_2472, t_2473, pc_x, qsg0_1327, qsg0_1328, qsg0_1329, \
                         qsg1_1327, qsg1_1328, qsg1_1329, qsh_1855, qsh_1856, \
                         qsh_1857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2471[k] = f_6 * qsg0_1327[k]
                    - f_7 * qsg1_1327[k]
                    + f_3 * pc_x[k] * qsh_1855[k];

        t_2472[k] = f_6 * qsg0_1328[k]
                    - f_7 * qsg1_1328[k]
                    + f_3 * pc_x[k] * qsh_1856[k];

        t_2473[k] = f_6 * qsg0_1329[k]
                    - f_7 * qsg1_1329[k]
                    + f_3 * pc_x[k] * qsh_1857[k];
    }

#pragma omp simd aligned(t_2474, t_2475, t_2476, pc_x, qsg0_1330, qsg0_1331, qsg0_1332, \
                         qsg1_1330, qsg1_1331, qsg1_1332, qsh_1858, qsh_1859, \
                         qsh_1860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2474[k] = f_4 * qsg0_1330[k]
                    - f_5 * qsg1_1330[k]
                    + f_3 * pc_x[k] * qsh_1858[k];

        t_2475[k] = f_4 * qsg0_1331[k]
                    - f_5 * qsg1_1331[k]
                    + f_3 * pc_x[k] * qsh_1859[k];

        t_2476[k] = f_4 * qsg0_1332[k]
                    - f_5 * qsg1_1332[k]
                    + f_3 * pc_x[k] * qsh_1860[k];
    }

#pragma omp simd aligned(t_2477, t_2478, t_2479, t_2480, t_2481, pc_x, qsg0_1333, qsg0_1334, \
                         qsg1_1333, qsg1_1334, qsh_1861, qsh_1862, qsh_1863, qsh_1864, \
                         qsh_1865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2477[k] = f_4 * qsg0_1333[k]
                    - f_5 * qsg1_1333[k]
                    + f_3 * pc_x[k] * qsh_1861[k];

        t_2478[k] = f_4 * qsg0_1334[k]
                    - f_5 * qsg1_1334[k]
                    + f_3 * pc_x[k] * qsh_1862[k];

        t_2479[k] = f_3 * pc_x[k] * qsh_1863[k];

        t_2480[k] = f_3 * pc_x[k] * qsh_1864[k];

        t_2481[k] = f_3 * pc_x[k] * qsh_1865[k];
    }

#pragma omp simd aligned(t_2482, t_2483, t_2484, t_2485, t_2486, pc_x, pc_y, pc_z, osh_1590, \
                         osh_1611, qsg0_1330, qsg1_1330, qsh_1863, qsh_1866, qsh_1867, \
                         qsh_1868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2482[k] = f_3 * pc_x[k] * qsh_1866[k];

        t_2483[k] = f_3 * pc_x[k] * qsh_1867[k];

        t_2484[k] = f_3 * pc_x[k] * qsh_1868[k];

        t_2485[k] = f_12 * osh_1611[k]
                    + f_1 * qsg0_1330[k]
                    - f_2 * qsg1_1330[k]
                    + f_3 * pc_y[k] * qsh_1863[k];

        t_2486[k] = f_18 * osh_1590[k]
                    + f_3 * pc_z[k] * qsh_1863[k];
    }

#pragma omp simd aligned(t_2487, t_2488, t_2489, pc_y, osh_1613, osh_1614, osh_1615, \
                         qsg0_1332, qsg0_1333, qsg0_1334, qsg1_1332, qsg1_1333, qsg1_1334, \
                         qsh_1865, qsh_1866, qsh_1867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2487[k] = f_12 * osh_1613[k]
                    + f_8 * qsg0_1332[k]
                    - f_9 * qsg1_1332[k]
                    + f_3 * pc_y[k] * qsh_1865[k];

        t_2488[k] = f_12 * osh_1614[k]
                    + f_6 * qsg0_1333[k]
                    - f_7 * qsg1_1333[k]
                    + f_3 * pc_y[k] * qsh_1866[k];

        t_2489[k] = f_12 * osh_1615[k]
                    + f_4 * qsg0_1334[k]
                    - f_5 * qsg1_1334[k]
                    + f_3 * pc_y[k] * qsh_1867[k];
    }

#pragma omp simd aligned(t_2490, t_2491, t_2492, pa_y, pc_y, pc_z, osi0_2156, osh_1595, \
                         osh_1616, osi1_2156, qsg0_1334, qsg1_1334, \
                         qsh_1868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2490[k] = f_12 * osh_1616[k]
                    + f_3 * pc_y[k] * qsh_1868[k];

        t_2491[k] = f_18 * osh_1595[k]
                    + f_1 * qsg0_1334[k]
                    - f_2 * qsg1_1334[k]
                    + f_3 * pc_z[k] * qsh_1868[k];

        t_2492[k] = pa_y[k] * osi0_2156[k]
                    - f_10 * pc_y[k] * osi1_2156[k];
    }

#pragma omp simd aligned(t_2493, t_2494, t_2495, pa_y, pc_x, pc_y, osi0_2158, osi1_2158, \
                         qsg0_1336, qsg0_1338, qsg1_1336, qsg1_1338, qsh_1870, \
                         qsh_1872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2493[k] = f_16 * qsg0_1336[k]
                    - f_17 * qsg1_1336[k]
                    + f_3 * pc_x[k] * qsh_1870[k];

        t_2494[k] = pa_y[k] * osi0_2158[k]
                    - f_10 * pc_y[k] * osi1_2158[k];

        t_2495[k] = f_8 * qsg0_1338[k]
                    - f_9 * qsg1_1338[k]
                    + f_3 * pc_x[k] * qsh_1872[k];
    }

#pragma omp simd aligned(t_2496, t_2497, t_2498, pa_y, pc_x, pc_y, osi0_2161, osi1_2161, \
                         qsg0_1339, qsg0_1341, qsg1_1339, qsg1_1341, qsh_1873, \
                         qsh_1875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2496[k] = f_8 * qsg0_1339[k]
                    - f_9 * qsg1_1339[k]
                    + f_3 * pc_x[k] * qsh_1873[k];

        t_2497[k] = pa_y[k] * osi0_2161[k]
                    - f_10 * pc_y[k] * osi1_2161[k];

        t_2498[k] = f_6 * qsg0_1341[k]
                    - f_7 * qsg1_1341[k]
                    + f_3 * pc_x[k] * qsh_1875[k];
    }

#pragma omp simd aligned(t_2499, t_2500, t_2501, pa_y, pc_x, pc_y, osi0_2165, osi1_2165, \
                         qsg0_1342, qsg0_1343, qsg1_1342, qsg1_1343, qsh_1876, \
                         qsh_1877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2499[k] = f_6 * qsg0_1342[k]
                    - f_7 * qsg1_1342[k]
                    + f_3 * pc_x[k] * qsh_1876[k];

        t_2500[k] = f_6 * qsg0_1343[k]
                    - f_7 * qsg1_1343[k]
                    + f_3 * pc_x[k] * qsh_1877[k];

        t_2501[k] = pa_y[k] * osi0_2165[k]
                    - f_10 * pc_y[k] * osi1_2165[k];
    }

#pragma omp simd aligned(t_2502, t_2503, t_2504, pc_x, qsg0_1345, qsg0_1346, qsg0_1347, \
                         qsg1_1345, qsg1_1346, qsg1_1347, qsh_1879, qsh_1880, \
                         qsh_1881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2502[k] = f_4 * qsg0_1345[k]
                    - f_5 * qsg1_1345[k]
                    + f_3 * pc_x[k] * qsh_1879[k];

        t_2503[k] = f_4 * qsg0_1346[k]
                    - f_5 * qsg1_1346[k]
                    + f_3 * pc_x[k] * qsh_1880[k];

        t_2504[k] = f_4 * qsg0_1347[k]
                    - f_5 * qsg1_1347[k]
                    + f_3 * pc_x[k] * qsh_1881[k];
    }

#pragma omp simd aligned(t_2505, t_2506, t_2507, t_2508, t_2509, pa_y, pc_x, pc_y, osi0_2170, \
                         osi1_2170, qsg0_1348, qsg1_1348, qsh_1882, qsh_1884, qsh_1885, \
                         qsh_1886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2505[k] = f_4 * qsg0_1348[k]
                    - f_5 * qsg1_1348[k]
                    + f_3 * pc_x[k] * qsh_1882[k];

        t_2506[k] = pa_y[k] * osi0_2170[k]
                    - f_10 * pc_y[k] * osi1_2170[k];

        t_2507[k] = f_3 * pc_x[k] * qsh_1884[k];

        t_2508[k] = f_3 * pc_x[k] * qsh_1885[k];

        t_2509[k] = f_3 * pc_x[k] * qsh_1886[k];
    }

#pragma omp simd aligned(t_2510, t_2511, t_2512, t_2513, pa_y, pc_x, pc_y, osi0_2177, \
                         osh_1632, osi1_2177, qsh_1887, qsh_1888, \
                         qsh_1889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2510[k] = f_3 * pc_x[k] * qsh_1887[k];

        t_2511[k] = f_3 * pc_x[k] * qsh_1888[k];

        t_2512[k] = f_3 * pc_x[k] * qsh_1889[k];

        t_2513[k] = pa_y[k] * osi0_2177[k]
                    + f_23 * osh_1632[k]
                    - f_10 * pc_y[k] * osi1_2177[k];
    }

#pragma omp simd aligned(t_2514, t_2515, t_2516, pa_y, pc_y, pc_z, osi0_2179, osi0_2180, \
                         osh_1611, osh_1634, osh_1635, osi1_2179, osi1_2180, \
                         qsh_1884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2514[k] = f_15 * osh_1611[k]
                    + f_3 * pc_z[k] * qsh_1884[k];

        t_2515[k] = pa_y[k] * osi0_2179[k]
                    + f_14 * osh_1634[k]
                    - f_10 * pc_y[k] * osi1_2179[k];

        t_2516[k] = pa_y[k] * osi0_2180[k]
                    + f_13 * osh_1635[k]
                    - f_10 * pc_y[k] * osi1_2180[k];
    }

#pragma omp simd aligned(t_2517, t_2518, t_2519, pa_y, pc_y, osi0_2181, osi0_2183, osh_1636, \
                         osh_1637, osi1_2181, osi1_2183, qsh_1889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2517[k] = pa_y[k] * osi0_2181[k]
                    + f_12 * osh_1636[k]
                    - f_10 * pc_y[k] * osi1_2181[k];

        t_2518[k] = f_11 * osh_1637[k]
                    + f_3 * pc_y[k] * qsh_1889[k];

        t_2519[k] = pa_y[k] * osi0_2183[k]
                    - f_10 * pc_y[k] * osi1_2183[k];
    }

#pragma omp simd aligned(t_2520, t_2521, t_2522, t_2523, t_2524, pc_x, pc_y, qsg0_1350, \
                         qsg0_1352, qsg0_1353, qsg1_1350, qsg1_1352, qsg1_1353, qsh_1890, \
                         qsh_1892, qsh_1893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2520[k] = f_1 * qsg0_1350[k]
                    - f_2 * qsg1_1350[k]
                    + f_3 * pc_x[k] * qsh_1890[k];

        t_2521[k] = f_3 * pc_y[k] * qsh_1890[k];

        t_2522[k] = f_16 * qsg0_1352[k]
                    - f_17 * qsg1_1352[k]
                    + f_3 * pc_x[k] * qsh_1892[k];

        t_2523[k] = f_8 * qsg0_1353[k]
                    - f_9 * qsg1_1353[k]
                    + f_3 * pc_x[k] * qsh_1893[k];

        t_2524[k] = f_3 * pc_y[k] * qsh_1892[k];
    }
}

static auto
compute_prim_qsi_three_center_electron_repulsion_0_piece22(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t osh, const size_t qsg0,
                                                           const size_t qsg1, const size_t qsh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

    auto *t_2525 = buffer.data(target + 2525);
    auto *t_2526 = buffer.data(target + 2526);
    auto *t_2527 = buffer.data(target + 2527);
    auto *t_2528 = buffer.data(target + 2528);
    auto *t_2529 = buffer.data(target + 2529);
    auto *t_2530 = buffer.data(target + 2530);
    auto *t_2531 = buffer.data(target + 2531);
    auto *t_2532 = buffer.data(target + 2532);
    auto *t_2533 = buffer.data(target + 2533);
    auto *t_2534 = buffer.data(target + 2534);
    auto *t_2535 = buffer.data(target + 2535);
    auto *t_2536 = buffer.data(target + 2536);
    auto *t_2537 = buffer.data(target + 2537);
    auto *t_2538 = buffer.data(target + 2538);
    auto *t_2539 = buffer.data(target + 2539);
    auto *t_2540 = buffer.data(target + 2540);
    auto *t_2541 = buffer.data(target + 2541);
    auto *t_2542 = buffer.data(target + 2542);
    auto *t_2543 = buffer.data(target + 2543);
    auto *t_2544 = buffer.data(target + 2544);
    auto *t_2545 = buffer.data(target + 2545);
    auto *t_2546 = buffer.data(target + 2546);
    auto *t_2547 = buffer.data(target + 2547);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *osh_1637 = buffer.data(osh + 1637);

    const auto *qsg0_1355 = buffer.data(qsg0 + 1355);
    const auto *qsg0_1356 = buffer.data(qsg0 + 1356);
    const auto *qsg0_1357 = buffer.data(qsg0 + 1357);
    const auto *qsg0_1359 = buffer.data(qsg0 + 1359);
    const auto *qsg0_1360 = buffer.data(qsg0 + 1360);
    const auto *qsg0_1361 = buffer.data(qsg0 + 1361);
    const auto *qsg0_1362 = buffer.data(qsg0 + 1362);
    const auto *qsg0_1363 = buffer.data(qsg0 + 1363);
    const auto *qsg0_1364 = buffer.data(qsg0 + 1364);

    const auto *qsg1_1355 = buffer.data(qsg1 + 1355);
    const auto *qsg1_1356 = buffer.data(qsg1 + 1356);
    const auto *qsg1_1357 = buffer.data(qsg1 + 1357);
    const auto *qsg1_1359 = buffer.data(qsg1 + 1359);
    const auto *qsg1_1360 = buffer.data(qsg1 + 1360);
    const auto *qsg1_1361 = buffer.data(qsg1 + 1361);
    const auto *qsg1_1362 = buffer.data(qsg1 + 1362);
    const auto *qsg1_1363 = buffer.data(qsg1 + 1363);
    const auto *qsg1_1364 = buffer.data(qsg1 + 1364);

    const auto *qsh_1895 = buffer.data(qsh + 1895);
    const auto *qsh_1896 = buffer.data(qsh + 1896);
    const auto *qsh_1897 = buffer.data(qsh + 1897);
    const auto *qsh_1899 = buffer.data(qsh + 1899);
    const auto *qsh_1900 = buffer.data(qsh + 1900);
    const auto *qsh_1901 = buffer.data(qsh + 1901);
    const auto *qsh_1902 = buffer.data(qsh + 1902);
    const auto *qsh_1904 = buffer.data(qsh + 1904);
    const auto *qsh_1905 = buffer.data(qsh + 1905);
    const auto *qsh_1906 = buffer.data(qsh + 1906);
    const auto *qsh_1907 = buffer.data(qsh + 1907);
    const auto *qsh_1908 = buffer.data(qsh + 1908);
    const auto *qsh_1909 = buffer.data(qsh + 1909);
    const auto *qsh_1910 = buffer.data(qsh + 1910);

#pragma omp simd aligned(t_2525, t_2526, t_2527, t_2528, pc_x, pc_y, qsg0_1355, qsg0_1356, \
                         qsg0_1357, qsg1_1355, qsg1_1356, qsg1_1357, qsh_1895, qsh_1896, \
                         qsh_1897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2525[k] = f_8 * qsg0_1355[k]
                    - f_9 * qsg1_1355[k]
                    + f_3 * pc_x[k] * qsh_1895[k];

        t_2526[k] = f_6 * qsg0_1356[k]
                    - f_7 * qsg1_1356[k]
                    + f_3 * pc_x[k] * qsh_1896[k];

        t_2527[k] = f_6 * qsg0_1357[k]
                    - f_7 * qsg1_1357[k]
                    + f_3 * pc_x[k] * qsh_1897[k];

        t_2528[k] = f_3 * pc_y[k] * qsh_1895[k];
    }

#pragma omp simd aligned(t_2529, t_2530, t_2531, pc_x, qsg0_1359, qsg0_1360, qsg0_1361, \
                         qsg1_1359, qsg1_1360, qsg1_1361, qsh_1899, qsh_1900, \
                         qsh_1901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2529[k] = f_6 * qsg0_1359[k]
                    - f_7 * qsg1_1359[k]
                    + f_3 * pc_x[k] * qsh_1899[k];

        t_2530[k] = f_4 * qsg0_1360[k]
                    - f_5 * qsg1_1360[k]
                    + f_3 * pc_x[k] * qsh_1900[k];

        t_2531[k] = f_4 * qsg0_1361[k]
                    - f_5 * qsg1_1361[k]
                    + f_3 * pc_x[k] * qsh_1901[k];
    }

#pragma omp simd aligned(t_2532, t_2533, t_2534, t_2535, t_2536, pc_x, pc_y, qsg0_1362, \
                         qsg0_1364, qsg1_1362, qsg1_1364, qsh_1899, qsh_1902, qsh_1904, \
                         qsh_1905, qsh_1906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2532[k] = f_4 * qsg0_1362[k]
                    - f_5 * qsg1_1362[k]
                    + f_3 * pc_x[k] * qsh_1902[k];

        t_2533[k] = f_3 * pc_y[k] * qsh_1899[k];

        t_2534[k] = f_4 * qsg0_1364[k]
                    - f_5 * qsg1_1364[k]
                    + f_3 * pc_x[k] * qsh_1904[k];

        t_2535[k] = f_3 * pc_x[k] * qsh_1905[k];

        t_2536[k] = f_3 * pc_x[k] * qsh_1906[k];
    }

#pragma omp simd aligned(t_2537, t_2538, t_2539, t_2540, t_2541, pc_x, pc_y, qsg0_1360, \
                         qsg1_1360, qsh_1905, qsh_1907, qsh_1908, qsh_1909, \
                         qsh_1910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2537[k] = f_3 * pc_x[k] * qsh_1907[k];

        t_2538[k] = f_3 * pc_x[k] * qsh_1908[k];

        t_2539[k] = f_3 * pc_x[k] * qsh_1909[k];

        t_2540[k] = f_3 * pc_x[k] * qsh_1910[k];

        t_2541[k] = f_1 * qsg0_1360[k]
                    - f_2 * qsg1_1360[k]
                    + f_3 * pc_y[k] * qsh_1905[k];
    }

#pragma omp simd aligned(t_2542, t_2543, t_2544, pc_y, qsg0_1361, qsg0_1362, qsg0_1363, \
                         qsg1_1361, qsg1_1362, qsg1_1363, qsh_1906, qsh_1907, \
                         qsh_1908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2542[k] = f_16 * qsg0_1361[k]
                    - f_17 * qsg1_1361[k]
                    + f_3 * pc_y[k] * qsh_1906[k];

        t_2543[k] = f_8 * qsg0_1362[k]
                    - f_9 * qsg1_1362[k]
                    + f_3 * pc_y[k] * qsh_1907[k];

        t_2544[k] = f_6 * qsg0_1363[k]
                    - f_7 * qsg1_1363[k]
                    + f_3 * pc_y[k] * qsh_1908[k];
    }

#pragma omp simd aligned(t_2545, t_2546, t_2547, pc_y, pc_z, osh_1637, qsg0_1364, qsg1_1364, \
                         qsh_1909, qsh_1910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2545[k] = f_4 * qsg0_1364[k]
                    - f_5 * qsg1_1364[k]
                    + f_3 * pc_y[k] * qsh_1909[k];

        t_2546[k] = f_3 * pc_y[k] * qsh_1910[k];

        t_2547[k] = f_0 * osh_1637[k]
                    + f_1 * qsg0_1364[k]
                    - f_2 * qsg1_1364[k]
                    + f_3 * pc_z[k] * qsh_1910[k];
    }
}

auto
compute_prim_qsi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t osi0, const size_t osh,
                                                   const size_t osi1, const size_t qsg0,
                                                   const size_t qsg1, const size_t qsh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_qsi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, osi0, osh,
                                                              osi1, qsg0, qsg1, qsh, ncols,
                                                              gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsg0, qsg1, qsh,
                                                               ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsg0, qsg1, qsh,
                                                               ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece12(buffer, target, pc, osh, qsg0,
                                                               qsg1, qsh, ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsg0, qsg1, qsh,
                                                               ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece14(buffer, target, pc, osh, qsg0,
                                                               qsg1, qsh, ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece15(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsg0, qsg1, qsh,
                                                               ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece16(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsg0, qsg1, qsh,
                                                               ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece17(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsh, ncols, gamma, p,
                                                               q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece18(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsg0, qsg1, qsh,
                                                               ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece19(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsg0, qsg1, qsh,
                                                               ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece20(buffer, target, pc, osh, qsg0,
                                                               qsg1, qsh, ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece21(buffer, target, pa, pc, osi0,
                                                               osh, osi1, qsg0, qsg1, qsh,
                                                               ncols, gamma, p, q);

    compute_prim_qsi_three_center_electron_repulsion_0_piece22(buffer, target, pc, osh, qsg0,
                                                               qsg1, qsh, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
