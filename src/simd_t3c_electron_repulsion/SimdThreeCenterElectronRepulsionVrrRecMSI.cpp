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


#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.5 / q;

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

    const auto *lsi0_0 = buffer.data(lsi0 + 0);
    const auto *lsi0_3 = buffer.data(lsi0 + 3);
    const auto *lsi0_5 = buffer.data(lsi0 + 5);
    const auto *lsi0_6 = buffer.data(lsi0 + 6);
    const auto *lsi0_9 = buffer.data(lsi0 + 9);
    const auto *lsi0_10 = buffer.data(lsi0 + 10);
    const auto *lsi0_14 = buffer.data(lsi0 + 14);
    const auto *lsi0_21 = buffer.data(lsi0 + 21);
    const auto *lsi0_27 = buffer.data(lsi0 + 27);
    const auto *lsi0_31 = buffer.data(lsi0 + 31);
    const auto *lsi0_34 = buffer.data(lsi0 + 34);
    const auto *lsi0_38 = buffer.data(lsi0 + 38);
    const auto *lsi0_56 = buffer.data(lsi0 + 56);
    const auto *lsi0_61 = buffer.data(lsi0 + 61);
    const auto *lsi0_65 = buffer.data(lsi0 + 65);
    const auto *lsi0_68 = buffer.data(lsi0 + 68);
    const auto *lsi0_70 = buffer.data(lsi0 + 70);

    const auto *lsh_0 = buffer.data(lsh + 0);
    const auto *lsh_1 = buffer.data(lsh + 1);
    const auto *lsh_2 = buffer.data(lsh + 2);
    const auto *lsh_3 = buffer.data(lsh + 3);
    const auto *lsh_5 = buffer.data(lsh + 5);
    const auto *lsh_6 = buffer.data(lsh + 6);
    const auto *lsh_9 = buffer.data(lsh + 9);
    const auto *lsh_15 = buffer.data(lsh + 15);
    const auto *lsh_17 = buffer.data(lsh + 17);
    const auto *lsh_18 = buffer.data(lsh + 18);
    const auto *lsh_20 = buffer.data(lsh + 20);
    const auto *lsh_21 = buffer.data(lsh + 21);
    const auto *lsh_24 = buffer.data(lsh + 24);
    const auto *lsh_26 = buffer.data(lsh + 26);
    const auto *lsh_27 = buffer.data(lsh + 27);
    const auto *lsh_30 = buffer.data(lsh + 30);
    const auto *lsh_36 = buffer.data(lsh + 36);
    const auto *lsh_38 = buffer.data(lsh + 38);
    const auto *lsh_39 = buffer.data(lsh + 39);
    const auto *lsh_40 = buffer.data(lsh + 40);
    const auto *lsh_41 = buffer.data(lsh + 41);
    const auto *lsh_42 = buffer.data(lsh + 42);
    const auto *lsh_44 = buffer.data(lsh + 44);
    const auto *lsh_47 = buffer.data(lsh + 47);
    const auto *lsh_50 = buffer.data(lsh + 50);
    const auto *lsh_51 = buffer.data(lsh + 51);
    const auto *lsh_57 = buffer.data(lsh + 57);
    const auto *lsh_58 = buffer.data(lsh + 58);
    const auto *lsh_59 = buffer.data(lsh + 59);
    const auto *lsh_60 = buffer.data(lsh + 60);
    const auto *lsh_62 = buffer.data(lsh + 62);
    const auto *lsh_63 = buffer.data(lsh + 63);
    const auto *lsh_66 = buffer.data(lsh + 66);
    const auto *lsh_69 = buffer.data(lsh + 69);
    const auto *lsh_73 = buffer.data(lsh + 73);
    const auto *lsh_78 = buffer.data(lsh + 78);
    const auto *lsh_80 = buffer.data(lsh + 80);
    const auto *lsh_81 = buffer.data(lsh + 81);
    const auto *lsh_82 = buffer.data(lsh + 82);
    const auto *lsh_83 = buffer.data(lsh + 83);
    const auto *lsh_99 = buffer.data(lsh + 99);

    const auto *lsi1_0 = buffer.data(lsi1 + 0);
    const auto *lsi1_3 = buffer.data(lsi1 + 3);
    const auto *lsi1_5 = buffer.data(lsi1 + 5);
    const auto *lsi1_6 = buffer.data(lsi1 + 6);
    const auto *lsi1_9 = buffer.data(lsi1 + 9);
    const auto *lsi1_10 = buffer.data(lsi1 + 10);
    const auto *lsi1_14 = buffer.data(lsi1 + 14);
    const auto *lsi1_21 = buffer.data(lsi1 + 21);
    const auto *lsi1_27 = buffer.data(lsi1 + 27);
    const auto *lsi1_31 = buffer.data(lsi1 + 31);
    const auto *lsi1_34 = buffer.data(lsi1 + 34);
    const auto *lsi1_38 = buffer.data(lsi1 + 38);
    const auto *lsi1_56 = buffer.data(lsi1 + 56);
    const auto *lsi1_61 = buffer.data(lsi1 + 61);
    const auto *lsi1_65 = buffer.data(lsi1 + 65);
    const auto *lsi1_68 = buffer.data(lsi1 + 68);
    const auto *lsi1_70 = buffer.data(lsi1 + 70);

    const auto *msg0_0 = buffer.data(msg0 + 0);
    const auto *msg0_1 = buffer.data(msg0 + 1);
    const auto *msg0_2 = buffer.data(msg0 + 2);
    const auto *msg0_3 = buffer.data(msg0 + 3);
    const auto *msg0_5 = buffer.data(msg0 + 5);
    const auto *msg0_10 = buffer.data(msg0 + 10);
    const auto *msg0_12 = buffer.data(msg0 + 12);
    const auto *msg0_13 = buffer.data(msg0 + 13);
    const auto *msg0_14 = buffer.data(msg0 + 14);
    const auto *msg0_18 = buffer.data(msg0 + 18);
    const auto *msg0_25 = buffer.data(msg0 + 25);
    const auto *msg0_26 = buffer.data(msg0 + 26);
    const auto *msg0_27 = buffer.data(msg0 + 27);
    const auto *msg0_32 = buffer.data(msg0 + 32);
    const auto *msg0_34 = buffer.data(msg0 + 34);
    const auto *msg0_35 = buffer.data(msg0 + 35);
    const auto *msg0_41 = buffer.data(msg0 + 41);
    const auto *msg0_42 = buffer.data(msg0 + 42);
    const auto *msg0_43 = buffer.data(msg0 + 43);
    const auto *msg0_44 = buffer.data(msg0 + 44);
    const auto *msg0_45 = buffer.data(msg0 + 45);
    const auto *msg0_47 = buffer.data(msg0 + 47);
    const auto *msg0_48 = buffer.data(msg0 + 48);
    const auto *msg0_50 = buffer.data(msg0 + 50);
    const auto *msg0_51 = buffer.data(msg0 + 51);
    const auto *msg0_55 = buffer.data(msg0 + 55);
    const auto *msg0_56 = buffer.data(msg0 + 56);
    const auto *msg0_57 = buffer.data(msg0 + 57);
    const auto *msg0_59 = buffer.data(msg0 + 59);

    const auto *msg1_0 = buffer.data(msg1 + 0);
    const auto *msg1_1 = buffer.data(msg1 + 1);
    const auto *msg1_2 = buffer.data(msg1 + 2);
    const auto *msg1_3 = buffer.data(msg1 + 3);
    const auto *msg1_5 = buffer.data(msg1 + 5);
    const auto *msg1_10 = buffer.data(msg1 + 10);
    const auto *msg1_12 = buffer.data(msg1 + 12);
    const auto *msg1_13 = buffer.data(msg1 + 13);
    const auto *msg1_14 = buffer.data(msg1 + 14);
    const auto *msg1_18 = buffer.data(msg1 + 18);
    const auto *msg1_25 = buffer.data(msg1 + 25);
    const auto *msg1_26 = buffer.data(msg1 + 26);
    const auto *msg1_27 = buffer.data(msg1 + 27);
    const auto *msg1_32 = buffer.data(msg1 + 32);
    const auto *msg1_34 = buffer.data(msg1 + 34);
    const auto *msg1_35 = buffer.data(msg1 + 35);
    const auto *msg1_41 = buffer.data(msg1 + 41);
    const auto *msg1_42 = buffer.data(msg1 + 42);
    const auto *msg1_43 = buffer.data(msg1 + 43);
    const auto *msg1_44 = buffer.data(msg1 + 44);
    const auto *msg1_45 = buffer.data(msg1 + 45);
    const auto *msg1_47 = buffer.data(msg1 + 47);
    const auto *msg1_48 = buffer.data(msg1 + 48);
    const auto *msg1_50 = buffer.data(msg1 + 50);
    const auto *msg1_51 = buffer.data(msg1 + 51);
    const auto *msg1_55 = buffer.data(msg1 + 55);
    const auto *msg1_56 = buffer.data(msg1 + 56);
    const auto *msg1_57 = buffer.data(msg1 + 57);
    const auto *msg1_59 = buffer.data(msg1 + 59);

    const auto *msh_0 = buffer.data(msh + 0);
    const auto *msh_1 = buffer.data(msh + 1);
    const auto *msh_2 = buffer.data(msh + 2);
    const auto *msh_3 = buffer.data(msh + 3);
    const auto *msh_5 = buffer.data(msh + 5);
    const auto *msh_6 = buffer.data(msh + 6);
    const auto *msh_8 = buffer.data(msh + 8);
    const auto *msh_9 = buffer.data(msh + 9);
    const auto *msh_10 = buffer.data(msh + 10);
    const auto *msh_14 = buffer.data(msh + 14);
    const auto *msh_15 = buffer.data(msh + 15);
    const auto *msh_17 = buffer.data(msh + 17);
    const auto *msh_18 = buffer.data(msh + 18);
    const auto *msh_19 = buffer.data(msh + 19);
    const auto *msh_20 = buffer.data(msh + 20);
    const auto *msh_21 = buffer.data(msh + 21);
    const auto *msh_22 = buffer.data(msh + 22);
    const auto *msh_24 = buffer.data(msh + 24);
    const auto *msh_26 = buffer.data(msh + 26);
    const auto *msh_27 = buffer.data(msh + 27);
    const auto *msh_28 = buffer.data(msh + 28);
    const auto *msh_30 = buffer.data(msh + 30);
    const auto *msh_31 = buffer.data(msh + 31);
    const auto *msh_36 = buffer.data(msh + 36);
    const auto *msh_37 = buffer.data(msh + 37);
    const auto *msh_38 = buffer.data(msh + 38);
    const auto *msh_39 = buffer.data(msh + 39);
    const auto *msh_40 = buffer.data(msh + 40);
    const auto *msh_41 = buffer.data(msh + 41);
    const auto *msh_42 = buffer.data(msh + 42);
    const auto *msh_44 = buffer.data(msh + 44);
    const auto *msh_46 = buffer.data(msh + 46);
    const auto *msh_47 = buffer.data(msh + 47);
    const auto *msh_49 = buffer.data(msh + 49);
    const auto *msh_50 = buffer.data(msh + 50);
    const auto *msh_51 = buffer.data(msh + 51);
    const auto *msh_56 = buffer.data(msh + 56);
    const auto *msh_57 = buffer.data(msh + 57);
    const auto *msh_58 = buffer.data(msh + 58);
    const auto *msh_59 = buffer.data(msh + 59);
    const auto *msh_60 = buffer.data(msh + 60);
    const auto *msh_61 = buffer.data(msh + 61);
    const auto *msh_62 = buffer.data(msh + 62);
    const auto *msh_63 = buffer.data(msh + 63);
    const auto *msh_64 = buffer.data(msh + 64);
    const auto *msh_65 = buffer.data(msh + 65);
    const auto *msh_66 = buffer.data(msh + 66);
    const auto *msh_68 = buffer.data(msh + 68);
    const auto *msh_69 = buffer.data(msh + 69);
    const auto *msh_70 = buffer.data(msh + 70);
    const auto *msh_72 = buffer.data(msh + 72);
    const auto *msh_73 = buffer.data(msh + 73);
    const auto *msh_78 = buffer.data(msh + 78);
    const auto *msh_79 = buffer.data(msh + 79);
    const auto *msh_80 = buffer.data(msh + 80);
    const auto *msh_81 = buffer.data(msh + 81);
    const auto *msh_82 = buffer.data(msh + 82);
    const auto *msh_83 = buffer.data(msh + 83);
    const auto *msh_84 = buffer.data(msh + 84);
    const auto *msh_86 = buffer.data(msh + 86);
    const auto *msh_87 = buffer.data(msh + 87);
    const auto *msh_89 = buffer.data(msh + 89);
    const auto *msh_90 = buffer.data(msh + 90);
    const auto *msh_93 = buffer.data(msh + 93);
    const auto *msh_99 = buffer.data(msh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, lsh_0, msg0_0, \
                         msg1_0, msh_0, msh_1, msh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lsh_0[k]
                 + f_1 * msg0_0[k]
                 - f_2 * msg1_0[k]
                 + f_3 * pc_x[k] * msh_0[k];

        t_1[k] = f_3 * pc_y[k] * msh_0[k];

        t_2[k] = f_3 * pc_z[k] * msh_0[k];

        t_3[k] = f_4 * msg0_0[k]
                 - f_5 * msg1_0[k]
                 + f_3 * pc_y[k] * msh_1[k];

        t_4[k] = f_3 * pc_y[k] * msh_2[k];

        t_5[k] = f_4 * msg0_0[k]
                 - f_5 * msg1_0[k]
                 + f_3 * pc_z[k] * msh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, msg0_1, msg0_2, msg0_3, msg1_1, \
                         msg1_2, msg1_3, msh_3, msh_5, msh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * msg0_1[k]
                 - f_7 * msg1_1[k]
                 + f_3 * pc_y[k] * msh_3[k];

        t_7[k] = f_3 * pc_z[k] * msh_3[k];

        t_8[k] = f_3 * pc_y[k] * msh_5[k];

        t_9[k] = f_6 * msg0_2[k]
                 - f_7 * msg1_2[k]
                 + f_3 * pc_z[k] * msh_5[k];

        t_10[k] = f_8 * msg0_3[k]
                  - f_9 * msg1_3[k]
                  + f_3 * pc_y[k] * msh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, lsh_15, msg0_5, \
                         msg1_5, msh_6, msh_8, msh_9, msh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * msh_6[k];

        t_12[k] = f_4 * msg0_5[k]
                  - f_5 * msg1_5[k]
                  + f_3 * pc_y[k] * msh_8[k];

        t_13[k] = f_3 * pc_y[k] * msh_9[k];

        t_14[k] = f_8 * msg0_5[k]
                  - f_9 * msg1_5[k]
                  + f_3 * pc_z[k] * msh_9[k];

        t_15[k] = f_0 * lsh_15[k]
                  + f_3 * pc_x[k] * msh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, lsh_17, lsh_18, \
                         lsh_20, msh_10, msh_14, msh_17, msh_18, \
                         msh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * msh_10[k];

        t_17[k] = f_0 * lsh_17[k]
                  + f_3 * pc_x[k] * msh_17[k];

        t_18[k] = f_0 * lsh_18[k]
                  + f_3 * pc_x[k] * msh_18[k];

        t_19[k] = f_3 * pc_y[k] * msh_14[k];

        t_20[k] = f_0 * lsh_20[k]
                  + f_3 * pc_x[k] * msh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, msg0_10, msg0_12, msg0_13, \
                         msg1_10, msg1_12, msg1_13, msh_15, msh_17, \
                         msh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * msg0_10[k]
                  - f_2 * msg1_10[k]
                  + f_3 * pc_y[k] * msh_15[k];

        t_22[k] = f_3 * pc_z[k] * msh_15[k];

        t_23[k] = f_8 * msg0_12[k]
                  - f_9 * msg1_12[k]
                  + f_3 * pc_y[k] * msh_17[k];

        t_24[k] = f_6 * msg0_13[k]
                  - f_7 * msg1_13[k]
                  + f_3 * pc_y[k] * msh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, lsi0_0, lsh_0, \
                         lsi1_0, msg0_14, msg1_14, msh_19, msh_20, \
                         msh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * msg0_14[k]
                  - f_5 * msg1_14[k]
                  + f_3 * pc_y[k] * msh_19[k];

        t_26[k] = f_3 * pc_y[k] * msh_20[k];

        t_27[k] = f_1 * msg0_14[k]
                  - f_2 * msg1_14[k]
                  + f_3 * pc_z[k] * msh_20[k];

        t_28[k] = pa_y[k] * lsi0_0[k]
                  - f_10 * pc_y[k] * lsi1_0[k];

        t_29[k] = f_11 * lsh_0[k]
                  + f_3 * pc_y[k] * msh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, lsi0_3, lsi0_5, lsh_1, \
                         lsi1_3, lsi1_5, msh_21, msh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * msh_21[k];

        t_31[k] = pa_y[k] * lsi0_3[k]
                  + f_12 * lsh_1[k]
                  - f_10 * pc_y[k] * lsi1_3[k];

        t_32[k] = f_3 * pc_z[k] * msh_22[k];

        t_33[k] = pa_y[k] * lsi0_5[k]
                  - f_10 * pc_y[k] * lsi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, lsi0_6, lsi0_9, lsh_3, \
                         lsh_5, lsi1_6, lsi1_9, msh_24, msh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * lsi0_6[k]
                  + f_13 * lsh_3[k]
                  - f_10 * pc_y[k] * lsi1_6[k];

        t_35[k] = f_3 * pc_z[k] * msh_24[k];

        t_36[k] = f_11 * lsh_5[k]
                  + f_3 * pc_y[k] * msh_26[k];

        t_37[k] = pa_y[k] * lsi0_9[k]
                  - f_10 * pc_y[k] * lsi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, lsi0_10, lsh_6, lsh_9, \
                         lsi1_10, msg0_18, msg1_18, msh_27, msh_28, \
                         msh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * lsi0_10[k]
                  + f_14 * lsh_6[k]
                  - f_10 * pc_y[k] * lsi1_10[k];

        t_39[k] = f_3 * pc_z[k] * msh_27[k];

        t_40[k] = f_4 * msg0_18[k]
                  - f_5 * msg1_18[k]
                  + f_3 * pc_z[k] * msh_28[k];

        t_41[k] = f_11 * lsh_9[k]
                  + f_3 * pc_y[k] * msh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, lsi0_14, lsh_36, \
                         lsh_38, lsi1_14, msh_31, msh_36, msh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * lsi0_14[k]
                  - f_10 * pc_y[k] * lsi1_14[k];

        t_43[k] = f_15 * lsh_36[k]
                  + f_3 * pc_x[k] * msh_36[k];

        t_44[k] = f_3 * pc_z[k] * msh_31[k];

        t_45[k] = f_15 * lsh_38[k]
                  + f_3 * pc_x[k] * msh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, lsh_15, lsh_39, lsh_40, lsh_41, \
                         msg0_25, msg1_25, msh_36, msh_39, msh_40, \
                         msh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * lsh_39[k]
                  + f_3 * pc_x[k] * msh_39[k];

        t_47[k] = f_15 * lsh_40[k]
                  + f_3 * pc_x[k] * msh_40[k];

        t_48[k] = f_15 * lsh_41[k]
                  + f_3 * pc_x[k] * msh_41[k];

        t_49[k] = f_11 * lsh_15[k]
                  + f_1 * msg0_25[k]
                  - f_2 * msg1_25[k]
                  + f_3 * pc_y[k] * msh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, msg0_25, msg0_26, msg0_27, msg1_25, \
                         msg1_26, msg1_27, msh_36, msh_37, msh_38, \
                         msh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * msh_36[k];

        t_51[k] = f_4 * msg0_25[k]
                  - f_5 * msg1_25[k]
                  + f_3 * pc_z[k] * msh_37[k];

        t_52[k] = f_6 * msg0_26[k]
                  - f_7 * msg1_26[k]
                  + f_3 * pc_z[k] * msh_38[k];

        t_53[k] = f_8 * msg0_27[k]
                  - f_9 * msg1_27[k]
                  + f_3 * pc_z[k] * msh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, lsi0_0, lsi0_27, \
                         lsh_20, lsi1_0, lsi1_27, msh_41, msh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * lsh_20[k]
                  + f_3 * pc_y[k] * msh_41[k];

        t_55[k] = pa_y[k] * lsi0_27[k]
                  - f_10 * pc_y[k] * lsi1_27[k];

        t_56[k] = pa_z[k] * lsi0_0[k]
                  - f_10 * pc_z[k] * lsi1_0[k];

        t_57[k] = f_3 * pc_y[k] * msh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, lsi0_3, lsi0_5, lsh_0, \
                         lsh_2, lsi1_3, lsi1_5, msh_42, msh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * lsh_0[k]
                  + f_3 * pc_z[k] * msh_42[k];

        t_59[k] = pa_z[k] * lsi0_3[k]
                  - f_10 * pc_z[k] * lsi1_3[k];

        t_60[k] = f_3 * pc_y[k] * msh_44[k];

        t_61[k] = pa_z[k] * lsi0_5[k]
                  + f_12 * lsh_2[k]
                  - f_10 * pc_z[k] * lsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, lsi0_6, lsi0_9, lsh_5, \
                         lsi1_6, lsi1_9, msg0_32, msg1_32, msh_46, \
                         msh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * lsi0_6[k]
                  - f_10 * pc_z[k] * lsi1_6[k];

        t_63[k] = f_4 * msg0_32[k]
                  - f_5 * msg1_32[k]
                  + f_3 * pc_y[k] * msh_46[k];

        t_64[k] = f_3 * pc_y[k] * msh_47[k];

        t_65[k] = pa_z[k] * lsi0_9[k]
                  + f_13 * lsh_5[k]
                  - f_10 * pc_z[k] * lsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, lsi0_10, lsi1_10, msg0_34, \
                         msg0_35, msg1_34, msg1_35, msh_49, msh_50, \
                         msh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * lsi0_10[k]
                  - f_10 * pc_z[k] * lsi1_10[k];

        t_67[k] = f_6 * msg0_34[k]
                  - f_7 * msg1_34[k]
                  + f_3 * pc_y[k] * msh_49[k];

        t_68[k] = f_4 * msg0_35[k]
                  - f_5 * msg1_35[k]
                  + f_3 * pc_y[k] * msh_50[k];

        t_69[k] = f_3 * pc_y[k] * msh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, lsi0_14, lsh_9, lsh_57, \
                         lsh_58, lsh_59, lsi1_14, msh_57, msh_58, \
                         msh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * lsi0_14[k]
                  + f_14 * lsh_9[k]
                  - f_10 * pc_z[k] * lsi1_14[k];

        t_71[k] = f_15 * lsh_57[k]
                  + f_3 * pc_x[k] * msh_57[k];

        t_72[k] = f_15 * lsh_58[k]
                  + f_3 * pc_x[k] * msh_58[k];

        t_73[k] = f_15 * lsh_59[k]
                  + f_3 * pc_x[k] * msh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, lsi0_21, lsh_60, \
                         lsh_62, lsi1_21, msh_56, msh_60, msh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_15 * lsh_60[k]
                  + f_3 * pc_x[k] * msh_60[k];

        t_75[k] = f_3 * pc_y[k] * msh_56[k];

        t_76[k] = f_15 * lsh_62[k]
                  + f_3 * pc_x[k] * msh_62[k];

        t_77[k] = pa_z[k] * lsi0_21[k]
                  - f_10 * pc_z[k] * lsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, msg0_41, msg0_42, msg0_43, msg1_41, msg1_42, \
                         msg1_43, msh_58, msh_59, msh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * msg0_41[k]
                  - f_17 * msg1_41[k]
                  + f_3 * pc_y[k] * msh_58[k];

        t_79[k] = f_8 * msg0_42[k]
                  - f_9 * msg1_42[k]
                  + f_3 * pc_y[k] * msh_59[k];

        t_80[k] = f_6 * msg0_43[k]
                  - f_7 * msg1_43[k]
                  + f_3 * pc_y[k] * msh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, lsh_20, lsh_63, msg0_44, \
                         msg0_45, msg1_44, msg1_45, msh_61, msh_62, \
                         msh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * msg0_44[k]
                  - f_5 * msg1_44[k]
                  + f_3 * pc_y[k] * msh_61[k];

        t_82[k] = f_3 * pc_y[k] * msh_62[k];

        t_83[k] = f_11 * lsh_20[k]
                  + f_1 * msg0_44[k]
                  - f_2 * msg1_44[k]
                  + f_3 * pc_z[k] * msh_62[k];

        t_84[k] = f_18 * lsh_63[k]
                  + f_1 * msg0_45[k]
                  - f_2 * msg1_45[k]
                  + f_3 * pc_x[k] * msh_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, lsh_21, lsh_66, msg0_48, \
                         msg1_48, msh_63, msh_64, msh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * lsh_21[k]
                  + f_3 * pc_y[k] * msh_63[k];

        t_86[k] = f_3 * pc_z[k] * msh_63[k];

        t_87[k] = f_18 * lsh_66[k]
                  + f_8 * msg0_48[k]
                  - f_9 * msg1_48[k]
                  + f_3 * pc_x[k] * msh_66[k];

        t_88[k] = f_3 * pc_z[k] * msh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, lsh_69, msg0_45, msg0_51, msg1_45, \
                         msg1_51, msh_65, msh_66, msh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * msg0_45[k]
                  - f_5 * msg1_45[k]
                  + f_3 * pc_z[k] * msh_65[k];

        t_90[k] = f_18 * lsh_69[k]
                  + f_6 * msg0_51[k]
                  - f_7 * msg1_51[k]
                  + f_3 * pc_x[k] * msh_69[k];

        t_91[k] = f_3 * pc_z[k] * msh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, lsh_26, lsh_73, msg0_47, \
                         msg0_55, msg1_47, msg1_55, msh_68, msh_69, \
                         msh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * lsh_26[k]
                  + f_3 * pc_y[k] * msh_68[k];

        t_93[k] = f_6 * msg0_47[k]
                  - f_7 * msg1_47[k]
                  + f_3 * pc_z[k] * msh_68[k];

        t_94[k] = f_18 * lsh_73[k]
                  + f_4 * msg0_55[k]
                  - f_5 * msg1_55[k]
                  + f_3 * pc_x[k] * msh_73[k];

        t_95[k] = f_3 * pc_z[k] * msh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, lsh_30, lsh_78, msg0_48, \
                         msg0_50, msg1_48, msg1_50, msh_70, msh_72, \
                         msh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * msg0_48[k]
                  - f_5 * msg1_48[k]
                  + f_3 * pc_z[k] * msh_70[k];

        t_97[k] = f_12 * lsh_30[k]
                  + f_3 * pc_y[k] * msh_72[k];

        t_98[k] = f_8 * msg0_50[k]
                  - f_9 * msg1_50[k]
                  + f_3 * pc_z[k] * msh_72[k];

        t_99[k] = f_18 * lsh_78[k]
                  + f_3 * pc_x[k] * msh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, lsh_80, lsh_81, \
                         lsh_82, lsh_83, msh_73, msh_80, msh_81, msh_82, \
                         msh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * msh_73[k];

        t_101[k] = f_18 * lsh_80[k]
                   + f_3 * pc_x[k] * msh_80[k];

        t_102[k] = f_18 * lsh_81[k]
                   + f_3 * pc_x[k] * msh_81[k];

        t_103[k] = f_18 * lsh_82[k]
                   + f_3 * pc_x[k] * msh_82[k];

        t_104[k] = f_18 * lsh_83[k]
                   + f_3 * pc_x[k] * msh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, lsh_36, msg0_55, msg0_56, \
                         msg1_55, msg1_56, msh_78, msh_79, msh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * lsh_36[k]
                   + f_1 * msg0_55[k]
                   - f_2 * msg1_55[k]
                   + f_3 * pc_y[k] * msh_78[k];

        t_106[k] = f_3 * pc_z[k] * msh_78[k];

        t_107[k] = f_4 * msg0_55[k]
                   - f_5 * msg1_55[k]
                   + f_3 * pc_z[k] * msh_79[k];

        t_108[k] = f_6 * msg0_56[k]
                   - f_7 * msg1_56[k]
                   + f_3 * pc_z[k] * msh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, lsi0_56, lsh_41, \
                         lsi1_56, msg0_57, msg0_59, msg1_57, msg1_59, msh_81, \
                         msh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * msg0_57[k]
                   - f_9 * msg1_57[k]
                   + f_3 * pc_z[k] * msh_81[k];

        t_110[k] = f_12 * lsh_41[k]
                   + f_3 * pc_y[k] * msh_83[k];

        t_111[k] = f_1 * msg0_59[k]
                   - f_2 * msg1_59[k]
                   + f_3 * pc_z[k] * msh_83[k];

        t_112[k] = pa_y[k] * lsi0_56[k]
                   - f_10 * pc_y[k] * lsi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, lsi0_31, lsh_21, \
                         lsh_42, lsh_44, lsi1_31, msh_84, msh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * lsh_42[k]
                   + f_3 * pc_y[k] * msh_84[k];

        t_114[k] = f_11 * lsh_21[k]
                   + f_3 * pc_z[k] * msh_84[k];

        t_115[k] = pa_z[k] * lsi0_31[k]
                   - f_10 * pc_z[k] * lsi1_31[k];

        t_116[k] = f_11 * lsh_44[k]
                   + f_3 * pc_y[k] * msh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, lsi0_34, lsi0_61, \
                         lsh_24, lsh_47, lsi1_34, lsi1_61, msh_87, \
                         msh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * lsi0_61[k]
                   - f_10 * pc_y[k] * lsi1_61[k];

        t_118[k] = pa_z[k] * lsi0_34[k]
                   - f_10 * pc_z[k] * lsi1_34[k];

        t_119[k] = f_11 * lsh_24[k]
                   + f_3 * pc_z[k] * msh_87[k];

        t_120[k] = f_11 * lsh_47[k]
                   + f_3 * pc_y[k] * msh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, lsi0_38, lsi0_65, \
                         lsh_27, lsi1_38, lsi1_65, msh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * lsi0_65[k]
                   - f_10 * pc_y[k] * lsi1_65[k];

        t_122[k] = pa_z[k] * lsi0_38[k]
                   - f_10 * pc_z[k] * lsi1_38[k];

        t_123[k] = f_11 * lsh_27[k]
                   + f_3 * pc_z[k] * msh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, lsi0_68, lsi0_70, \
                         lsh_50, lsh_51, lsh_99, lsi1_68, lsi1_70, msh_93, \
                         msh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * lsi0_68[k]
                   + f_12 * lsh_50[k]
                   - f_10 * pc_y[k] * lsi1_68[k];

        t_125[k] = f_11 * lsh_51[k]
                   + f_3 * pc_y[k] * msh_93[k];

        t_126[k] = pa_y[k] * lsi0_70[k]
                   - f_10 * pc_y[k] * lsi1_70[k];

        t_127[k] = f_18 * lsh_99[k]
                   + f_3 * pc_x[k] * msh_99[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *lsi0_49 = buffer.data(lsi0 + 49);
    const auto *lsi0_83 = buffer.data(lsi0 + 83);
    const auto *lsi0_84 = buffer.data(lsi0 + 84);
    const auto *lsi0_87 = buffer.data(lsi0 + 87);
    const auto *lsi0_90 = buffer.data(lsi0 + 90);
    const auto *lsi0_94 = buffer.data(lsi0 + 94);
    const auto *lsi0_96 = buffer.data(lsi0 + 96);
    const auto *lsi0_105 = buffer.data(lsi0 + 105);
    const auto *lsi0_140 = buffer.data(lsi0 + 140);
    const auto *lsi0_143 = buffer.data(lsi0 + 143);
    const auto *lsi0_145 = buffer.data(lsi0 + 145);
    const auto *lsi0_146 = buffer.data(lsi0 + 146);
    const auto *lsi0_149 = buffer.data(lsi0 + 149);
    const auto *lsi0_150 = buffer.data(lsi0 + 150);
    const auto *lsi0_152 = buffer.data(lsi0 + 152);
    const auto *lsi0_154 = buffer.data(lsi0 + 154);

    const auto *lsh_36 = buffer.data(lsh + 36);
    const auto *lsh_42 = buffer.data(lsh + 42);
    const auto *lsh_59 = buffer.data(lsh + 59);
    const auto *lsh_60 = buffer.data(lsh + 60);
    const auto *lsh_61 = buffer.data(lsh + 61);
    const auto *lsh_62 = buffer.data(lsh + 62);
    const auto *lsh_63 = buffer.data(lsh + 63);
    const auto *lsh_66 = buffer.data(lsh + 66);
    const auto *lsh_68 = buffer.data(lsh + 68);
    const auto *lsh_69 = buffer.data(lsh + 69);
    const auto *lsh_70 = buffer.data(lsh + 70);
    const auto *lsh_72 = buffer.data(lsh + 72);
    const auto *lsh_78 = buffer.data(lsh + 78);
    const auto *lsh_83 = buffer.data(lsh + 83);
    const auto *lsh_84 = buffer.data(lsh + 84);
    const auto *lsh_86 = buffer.data(lsh + 86);
    const auto *lsh_87 = buffer.data(lsh + 87);
    const auto *lsh_89 = buffer.data(lsh + 89);
    const auto *lsh_90 = buffer.data(lsh + 90);
    const auto *lsh_93 = buffer.data(lsh + 93);
    const auto *lsh_99 = buffer.data(lsh + 99);
    const auto *lsh_100 = buffer.data(lsh + 100);
    const auto *lsh_101 = buffer.data(lsh + 101);
    const auto *lsh_102 = buffer.data(lsh + 102);
    const auto *lsh_103 = buffer.data(lsh + 103);
    const auto *lsh_104 = buffer.data(lsh + 104);
    const auto *lsh_105 = buffer.data(lsh + 105);
    const auto *lsh_106 = buffer.data(lsh + 106);
    const auto *lsh_107 = buffer.data(lsh + 107);
    const auto *lsh_108 = buffer.data(lsh + 108);
    const auto *lsh_110 = buffer.data(lsh + 110);
    const auto *lsh_111 = buffer.data(lsh + 111);
    const auto *lsh_113 = buffer.data(lsh + 113);
    const auto *lsh_114 = buffer.data(lsh + 114);
    const auto *lsh_119 = buffer.data(lsh + 119);
    const auto *lsh_120 = buffer.data(lsh + 120);
    const auto *lsh_121 = buffer.data(lsh + 121);
    const auto *lsh_122 = buffer.data(lsh + 122);
    const auto *lsh_123 = buffer.data(lsh + 123);
    const auto *lsh_125 = buffer.data(lsh + 125);
    const auto *lsh_126 = buffer.data(lsh + 126);
    const auto *lsh_129 = buffer.data(lsh + 129);
    const auto *lsh_132 = buffer.data(lsh + 132);
    const auto *lsh_136 = buffer.data(lsh + 136);
    const auto *lsh_141 = buffer.data(lsh + 141);
    const auto *lsh_143 = buffer.data(lsh + 143);
    const auto *lsh_144 = buffer.data(lsh + 144);
    const auto *lsh_145 = buffer.data(lsh + 145);
    const auto *lsh_146 = buffer.data(lsh + 146);
    const auto *lsh_152 = buffer.data(lsh + 152);
    const auto *lsh_156 = buffer.data(lsh + 156);
    const auto *lsh_161 = buffer.data(lsh + 161);
    const auto *lsh_162 = buffer.data(lsh + 162);
    const auto *lsh_163 = buffer.data(lsh + 163);
    const auto *lsh_164 = buffer.data(lsh + 164);
    const auto *lsh_165 = buffer.data(lsh + 165);
    const auto *lsh_166 = buffer.data(lsh + 166);
    const auto *lsh_167 = buffer.data(lsh + 167);
    const auto *lsh_183 = buffer.data(lsh + 183);
    const auto *lsh_184 = buffer.data(lsh + 184);
    const auto *lsh_185 = buffer.data(lsh + 185);
    const auto *lsh_186 = buffer.data(lsh + 186);
    const auto *lsh_187 = buffer.data(lsh + 187);
    const auto *lsh_188 = buffer.data(lsh + 188);

    const auto *lsi1_49 = buffer.data(lsi1 + 49);
    const auto *lsi1_83 = buffer.data(lsi1 + 83);
    const auto *lsi1_84 = buffer.data(lsi1 + 84);
    const auto *lsi1_87 = buffer.data(lsi1 + 87);
    const auto *lsi1_90 = buffer.data(lsi1 + 90);
    const auto *lsi1_94 = buffer.data(lsi1 + 94);
    const auto *lsi1_96 = buffer.data(lsi1 + 96);
    const auto *lsi1_105 = buffer.data(lsi1 + 105);
    const auto *lsi1_140 = buffer.data(lsi1 + 140);
    const auto *lsi1_143 = buffer.data(lsi1 + 143);
    const auto *lsi1_145 = buffer.data(lsi1 + 145);
    const auto *lsi1_146 = buffer.data(lsi1 + 146);
    const auto *lsi1_149 = buffer.data(lsi1 + 149);
    const auto *lsi1_150 = buffer.data(lsi1 + 150);
    const auto *lsi1_152 = buffer.data(lsi1 + 152);
    const auto *lsi1_154 = buffer.data(lsi1 + 154);

    const auto *msg0_72 = buffer.data(msg0 + 72);
    const auto *msg0_73 = buffer.data(msg0 + 73);
    const auto *msg0_74 = buffer.data(msg0 + 74);
    const auto *msg0_75 = buffer.data(msg0 + 75);
    const auto *msg0_76 = buffer.data(msg0 + 76);
    const auto *msg0_77 = buffer.data(msg0 + 77);
    const auto *msg0_78 = buffer.data(msg0 + 78);
    const auto *msg0_79 = buffer.data(msg0 + 79);
    const auto *msg0_80 = buffer.data(msg0 + 80);
    const auto *msg0_84 = buffer.data(msg0 + 84);
    const auto *msg0_85 = buffer.data(msg0 + 85);
    const auto *msg0_86 = buffer.data(msg0 + 86);
    const auto *msg0_87 = buffer.data(msg0 + 87);
    const auto *msg0_88 = buffer.data(msg0 + 88);
    const auto *msg0_89 = buffer.data(msg0 + 89);
    const auto *msg0_90 = buffer.data(msg0 + 90);
    const auto *msg0_92 = buffer.data(msg0 + 92);
    const auto *msg0_93 = buffer.data(msg0 + 93);
    const auto *msg0_95 = buffer.data(msg0 + 95);
    const auto *msg0_96 = buffer.data(msg0 + 96);
    const auto *msg0_100 = buffer.data(msg0 + 100);
    const auto *msg0_101 = buffer.data(msg0 + 101);
    const auto *msg0_102 = buffer.data(msg0 + 102);
    const auto *msg0_104 = buffer.data(msg0 + 104);
    const auto *msg0_110 = buffer.data(msg0 + 110);
    const auto *msg0_114 = buffer.data(msg0 + 114);
    const auto *msg0_117 = buffer.data(msg0 + 117);
    const auto *msg0_118 = buffer.data(msg0 + 118);
    const auto *msg0_119 = buffer.data(msg0 + 119);
    const auto *msg0_130 = buffer.data(msg0 + 130);
    const auto *msg0_132 = buffer.data(msg0 + 132);

    const auto *msg1_72 = buffer.data(msg1 + 72);
    const auto *msg1_73 = buffer.data(msg1 + 73);
    const auto *msg1_74 = buffer.data(msg1 + 74);
    const auto *msg1_75 = buffer.data(msg1 + 75);
    const auto *msg1_76 = buffer.data(msg1 + 76);
    const auto *msg1_77 = buffer.data(msg1 + 77);
    const auto *msg1_78 = buffer.data(msg1 + 78);
    const auto *msg1_79 = buffer.data(msg1 + 79);
    const auto *msg1_80 = buffer.data(msg1 + 80);
    const auto *msg1_84 = buffer.data(msg1 + 84);
    const auto *msg1_85 = buffer.data(msg1 + 85);
    const auto *msg1_86 = buffer.data(msg1 + 86);
    const auto *msg1_87 = buffer.data(msg1 + 87);
    const auto *msg1_88 = buffer.data(msg1 + 88);
    const auto *msg1_89 = buffer.data(msg1 + 89);
    const auto *msg1_90 = buffer.data(msg1 + 90);
    const auto *msg1_92 = buffer.data(msg1 + 92);
    const auto *msg1_93 = buffer.data(msg1 + 93);
    const auto *msg1_95 = buffer.data(msg1 + 95);
    const auto *msg1_96 = buffer.data(msg1 + 96);
    const auto *msg1_100 = buffer.data(msg1 + 100);
    const auto *msg1_101 = buffer.data(msg1 + 101);
    const auto *msg1_102 = buffer.data(msg1 + 102);
    const auto *msg1_104 = buffer.data(msg1 + 104);
    const auto *msg1_110 = buffer.data(msg1 + 110);
    const auto *msg1_114 = buffer.data(msg1 + 114);
    const auto *msg1_117 = buffer.data(msg1 + 117);
    const auto *msg1_118 = buffer.data(msg1 + 118);
    const auto *msg1_119 = buffer.data(msg1 + 119);
    const auto *msg1_130 = buffer.data(msg1 + 130);
    const auto *msg1_132 = buffer.data(msg1 + 132);

    const auto *msh_99 = buffer.data(msh + 99);
    const auto *msh_100 = buffer.data(msh + 100);
    const auto *msh_101 = buffer.data(msh + 101);
    const auto *msh_102 = buffer.data(msh + 102);
    const auto *msh_103 = buffer.data(msh + 103);
    const auto *msh_104 = buffer.data(msh + 104);
    const auto *msh_105 = buffer.data(msh + 105);
    const auto *msh_106 = buffer.data(msh + 106);
    const auto *msh_107 = buffer.data(msh + 107);
    const auto *msh_108 = buffer.data(msh + 108);
    const auto *msh_109 = buffer.data(msh + 109);
    const auto *msh_110 = buffer.data(msh + 110);
    const auto *msh_111 = buffer.data(msh + 111);
    const auto *msh_112 = buffer.data(msh + 112);
    const auto *msh_113 = buffer.data(msh + 113);
    const auto *msh_114 = buffer.data(msh + 114);
    const auto *msh_119 = buffer.data(msh + 119);
    const auto *msh_120 = buffer.data(msh + 120);
    const auto *msh_121 = buffer.data(msh + 121);
    const auto *msh_122 = buffer.data(msh + 122);
    const auto *msh_123 = buffer.data(msh + 123);
    const auto *msh_124 = buffer.data(msh + 124);
    const auto *msh_125 = buffer.data(msh + 125);
    const auto *msh_126 = buffer.data(msh + 126);
    const auto *msh_127 = buffer.data(msh + 127);
    const auto *msh_128 = buffer.data(msh + 128);
    const auto *msh_129 = buffer.data(msh + 129);
    const auto *msh_131 = buffer.data(msh + 131);
    const auto *msh_132 = buffer.data(msh + 132);
    const auto *msh_133 = buffer.data(msh + 133);
    const auto *msh_135 = buffer.data(msh + 135);
    const auto *msh_136 = buffer.data(msh + 136);
    const auto *msh_141 = buffer.data(msh + 141);
    const auto *msh_142 = buffer.data(msh + 142);
    const auto *msh_143 = buffer.data(msh + 143);
    const auto *msh_144 = buffer.data(msh + 144);
    const auto *msh_145 = buffer.data(msh + 145);
    const auto *msh_146 = buffer.data(msh + 146);
    const auto *msh_147 = buffer.data(msh + 147);
    const auto *msh_149 = buffer.data(msh + 149);
    const auto *msh_150 = buffer.data(msh + 150);
    const auto *msh_152 = buffer.data(msh + 152);
    const auto *msh_153 = buffer.data(msh + 153);
    const auto *msh_156 = buffer.data(msh + 156);
    const auto *msh_161 = buffer.data(msh + 161);
    const auto *msh_162 = buffer.data(msh + 162);
    const auto *msh_163 = buffer.data(msh + 163);
    const auto *msh_164 = buffer.data(msh + 164);
    const auto *msh_165 = buffer.data(msh + 165);
    const auto *msh_166 = buffer.data(msh + 166);
    const auto *msh_167 = buffer.data(msh + 167);
    const auto *msh_168 = buffer.data(msh + 168);
    const auto *msh_170 = buffer.data(msh + 170);
    const auto *msh_171 = buffer.data(msh + 171);
    const auto *msh_173 = buffer.data(msh + 173);
    const auto *msh_174 = buffer.data(msh + 174);
    const auto *msh_177 = buffer.data(msh + 177);
    const auto *msh_183 = buffer.data(msh + 183);
    const auto *msh_184 = buffer.data(msh + 184);
    const auto *msh_185 = buffer.data(msh + 185);
    const auto *msh_186 = buffer.data(msh + 186);
    const auto *msh_187 = buffer.data(msh + 187);
    const auto *msh_188 = buffer.data(msh + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, lsh_100, lsh_101, lsh_102, \
                         lsh_103, lsh_104, msh_100, msh_101, msh_102, msh_103, \
                         msh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_18 * lsh_100[k]
                   + f_3 * pc_x[k] * msh_100[k];

        t_129[k] = f_18 * lsh_101[k]
                   + f_3 * pc_x[k] * msh_101[k];

        t_130[k] = f_18 * lsh_102[k]
                   + f_3 * pc_x[k] * msh_102[k];

        t_131[k] = f_18 * lsh_103[k]
                   + f_3 * pc_x[k] * msh_103[k];

        t_132[k] = f_18 * lsh_104[k]
                   + f_3 * pc_x[k] * msh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, lsi0_49, lsh_36, lsh_59, \
                         lsi1_49, msg0_72, msg1_72, msh_99, msh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * lsi0_49[k]
                   - f_10 * pc_z[k] * lsi1_49[k];

        t_134[k] = f_11 * lsh_36[k]
                   + f_3 * pc_z[k] * msh_99[k];

        t_135[k] = f_11 * lsh_59[k]
                   + f_8 * msg0_72[k]
                   - f_9 * msg1_72[k]
                   + f_3 * pc_y[k] * msh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, lsh_60, lsh_61, lsh_62, msg0_73, msg0_74, \
                         msg1_73, msg1_74, msh_102, msh_103, msh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * lsh_60[k]
                   + f_6 * msg0_73[k]
                   - f_7 * msg1_73[k]
                   + f_3 * pc_y[k] * msh_102[k];

        t_137[k] = f_11 * lsh_61[k]
                   + f_4 * msg0_74[k]
                   - f_5 * msg1_74[k]
                   + f_3 * pc_y[k] * msh_103[k];

        t_138[k] = f_11 * lsh_62[k]
                   + f_3 * pc_y[k] * msh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, lsi0_83, lsh_42, \
                         lsh_105, lsi1_83, msg0_75, msg1_75, msh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * lsi0_83[k]
                   - f_10 * pc_y[k] * lsi1_83[k];

        t_140[k] = f_18 * lsh_105[k]
                   + f_1 * msg0_75[k]
                   - f_2 * msg1_75[k]
                   + f_3 * pc_x[k] * msh_105[k];

        t_141[k] = f_3 * pc_y[k] * msh_105[k];

        t_142[k] = f_12 * lsh_42[k]
                   + f_3 * pc_z[k] * msh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, lsh_110, msg0_75, msg0_80, msg1_75, \
                         msg1_80, msh_106, msh_107, msh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * msg0_75[k]
                   - f_5 * msg1_75[k]
                   + f_3 * pc_y[k] * msh_106[k];

        t_144[k] = f_3 * pc_y[k] * msh_107[k];

        t_145[k] = f_18 * lsh_110[k]
                   + f_8 * msg0_80[k]
                   - f_9 * msg1_80[k]
                   + f_3 * pc_x[k] * msh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, msg0_76, msg0_77, msg1_76, msg1_77, \
                         msh_108, msh_109, msh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * msg0_76[k]
                   - f_7 * msg1_76[k]
                   + f_3 * pc_y[k] * msh_108[k];

        t_147[k] = f_4 * msg0_77[k]
                   - f_5 * msg1_77[k]
                   + f_3 * pc_y[k] * msh_109[k];

        t_148[k] = f_3 * pc_y[k] * msh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, lsh_114, msg0_78, msg0_79, msg0_84, \
                         msg1_78, msg1_79, msg1_84, msh_111, msh_112, \
                         msh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_18 * lsh_114[k]
                   + f_6 * msg0_84[k]
                   - f_7 * msg1_84[k]
                   + f_3 * pc_x[k] * msh_114[k];

        t_150[k] = f_8 * msg0_78[k]
                   - f_9 * msg1_78[k]
                   + f_3 * pc_y[k] * msh_111[k];

        t_151[k] = f_6 * msg0_79[k]
                   - f_7 * msg1_79[k]
                   + f_3 * pc_y[k] * msh_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, lsh_119, lsh_120, msg0_80, \
                         msg0_89, msg1_80, msg1_89, msh_113, msh_114, msh_119, \
                         msh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * msg0_80[k]
                   - f_5 * msg1_80[k]
                   + f_3 * pc_y[k] * msh_113[k];

        t_153[k] = f_3 * pc_y[k] * msh_114[k];

        t_154[k] = f_18 * lsh_119[k]
                   + f_4 * msg0_89[k]
                   - f_5 * msg1_89[k]
                   + f_3 * pc_x[k] * msh_119[k];

        t_155[k] = f_18 * lsh_120[k]
                   + f_3 * pc_x[k] * msh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, lsh_121, lsh_122, \
                         lsh_123, lsh_125, msh_119, msh_121, msh_122, msh_123, \
                         msh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_18 * lsh_121[k]
                   + f_3 * pc_x[k] * msh_121[k];

        t_157[k] = f_18 * lsh_122[k]
                   + f_3 * pc_x[k] * msh_122[k];

        t_158[k] = f_18 * lsh_123[k]
                   + f_3 * pc_x[k] * msh_123[k];

        t_159[k] = f_3 * pc_y[k] * msh_119[k];

        t_160[k] = f_18 * lsh_125[k]
                   + f_3 * pc_x[k] * msh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, msg0_85, msg0_86, msg0_87, msg1_85, \
                         msg1_86, msg1_87, msh_120, msh_121, msh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * msg0_85[k]
                   - f_2 * msg1_85[k]
                   + f_3 * pc_y[k] * msh_120[k];

        t_162[k] = f_16 * msg0_86[k]
                   - f_17 * msg1_86[k]
                   + f_3 * pc_y[k] * msh_121[k];

        t_163[k] = f_8 * msg0_87[k]
                   - f_9 * msg1_87[k]
                   + f_3 * pc_y[k] * msh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, lsh_62, msg0_88, msg0_89, \
                         msg1_88, msg1_89, msh_123, msh_124, msh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * msg0_88[k]
                   - f_7 * msg1_88[k]
                   + f_3 * pc_y[k] * msh_123[k];

        t_165[k] = f_4 * msg0_89[k]
                   - f_5 * msg1_89[k]
                   + f_3 * pc_y[k] * msh_124[k];

        t_166[k] = f_3 * pc_y[k] * msh_125[k];

        t_167[k] = f_12 * lsh_62[k]
                   + f_1 * msg0_89[k]
                   - f_2 * msg1_89[k]
                   + f_3 * pc_z[k] * msh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, lsh_63, lsh_126, \
                         lsh_129, msg0_90, msg0_93, msg1_90, msg1_93, msh_126, \
                         msh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * lsh_126[k]
                   + f_1 * msg0_90[k]
                   - f_2 * msg1_90[k]
                   + f_3 * pc_x[k] * msh_126[k];

        t_169[k] = f_13 * lsh_63[k]
                   + f_3 * pc_y[k] * msh_126[k];

        t_170[k] = f_3 * pc_z[k] * msh_126[k];

        t_171[k] = f_19 * lsh_129[k]
                   + f_8 * msg0_93[k]
                   - f_9 * msg1_93[k]
                   + f_3 * pc_x[k] * msh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_x, pc_z, lsh_132, msg0_90, msg0_96, \
                         msg1_90, msg1_96, msh_127, msh_128, msh_129, \
                         msh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * msh_127[k];

        t_173[k] = f_4 * msg0_90[k]
                   - f_5 * msg1_90[k]
                   + f_3 * pc_z[k] * msh_128[k];

        t_174[k] = f_19 * lsh_132[k]
                   + f_6 * msg0_96[k]
                   - f_7 * msg1_96[k]
                   + f_3 * pc_x[k] * msh_132[k];

        t_175[k] = f_3 * pc_z[k] * msh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, lsh_68, lsh_136, \
                         msg0_92, msg0_100, msg1_92, msg1_100, msh_131, msh_132, \
                         msh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * lsh_68[k]
                   + f_3 * pc_y[k] * msh_131[k];

        t_177[k] = f_6 * msg0_92[k]
                   - f_7 * msg1_92[k]
                   + f_3 * pc_z[k] * msh_131[k];

        t_178[k] = f_19 * lsh_136[k]
                   + f_4 * msg0_100[k]
                   - f_5 * msg1_100[k]
                   + f_3 * pc_x[k] * msh_136[k];

        t_179[k] = f_3 * pc_z[k] * msh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, lsh_72, lsh_141, \
                         msg0_93, msg0_95, msg1_93, msg1_95, msh_133, msh_135, \
                         msh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * msg0_93[k]
                   - f_5 * msg1_93[k]
                   + f_3 * pc_z[k] * msh_133[k];

        t_181[k] = f_13 * lsh_72[k]
                   + f_3 * pc_y[k] * msh_135[k];

        t_182[k] = f_8 * msg0_95[k]
                   - f_9 * msg1_95[k]
                   + f_3 * pc_z[k] * msh_135[k];

        t_183[k] = f_19 * lsh_141[k]
                   + f_3 * pc_x[k] * msh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, lsh_143, lsh_144, \
                         lsh_145, lsh_146, msh_136, msh_143, msh_144, msh_145, \
                         msh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * msh_136[k];

        t_185[k] = f_19 * lsh_143[k]
                   + f_3 * pc_x[k] * msh_143[k];

        t_186[k] = f_19 * lsh_144[k]
                   + f_3 * pc_x[k] * msh_144[k];

        t_187[k] = f_19 * lsh_145[k]
                   + f_3 * pc_x[k] * msh_145[k];

        t_188[k] = f_19 * lsh_146[k]
                   + f_3 * pc_x[k] * msh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_y, pc_z, lsh_78, msg0_100, msg0_101, \
                         msg1_100, msg1_101, msh_141, msh_142, \
                         msh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * lsh_78[k]
                   + f_1 * msg0_100[k]
                   - f_2 * msg1_100[k]
                   + f_3 * pc_y[k] * msh_141[k];

        t_190[k] = f_3 * pc_z[k] * msh_141[k];

        t_191[k] = f_4 * msg0_100[k]
                   - f_5 * msg1_100[k]
                   + f_3 * pc_z[k] * msh_142[k];

        t_192[k] = f_6 * msg0_101[k]
                   - f_7 * msg1_101[k]
                   + f_3 * pc_z[k] * msh_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pc_y, pc_z, lsi0_84, lsh_83, \
                         lsi1_84, msg0_102, msg0_104, msg1_102, msg1_104, msh_144, \
                         msh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * msg0_102[k]
                   - f_9 * msg1_102[k]
                   + f_3 * pc_z[k] * msh_144[k];

        t_194[k] = f_13 * lsh_83[k]
                   + f_3 * pc_y[k] * msh_146[k];

        t_195[k] = f_1 * msg0_104[k]
                   - f_2 * msg1_104[k]
                   + f_3 * pc_z[k] * msh_146[k];

        t_196[k] = pa_z[k] * lsi0_84[k]
                   - f_10 * pc_z[k] * lsi1_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pc_y, pc_z, lsi0_87, lsh_63, \
                         lsh_84, lsh_86, lsi1_87, msh_147, msh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * lsh_84[k]
                   + f_3 * pc_y[k] * msh_147[k];

        t_198[k] = f_11 * lsh_63[k]
                   + f_3 * pc_z[k] * msh_147[k];

        t_199[k] = pa_z[k] * lsi0_87[k]
                   - f_10 * pc_z[k] * lsi1_87[k];

        t_200[k] = f_12 * lsh_86[k]
                   + f_3 * pc_y[k] * msh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pc_x, pc_z, lsi0_90, lsh_66, lsh_152, \
                         lsi1_90, msg0_110, msg1_110, msh_150, \
                         msh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_19 * lsh_152[k]
                   + f_8 * msg0_110[k]
                   - f_9 * msg1_110[k]
                   + f_3 * pc_x[k] * msh_152[k];

        t_202[k] = pa_z[k] * lsi0_90[k]
                   - f_10 * pc_z[k] * lsi1_90[k];

        t_203[k] = f_11 * lsh_66[k]
                   + f_3 * pc_z[k] * msh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_y, pc_z, lsi0_94, lsh_89, \
                         lsh_156, lsi1_94, msg0_114, msg1_114, msh_152, \
                         msh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * lsh_89[k]
                   + f_3 * pc_y[k] * msh_152[k];

        t_205[k] = f_19 * lsh_156[k]
                   + f_6 * msg0_114[k]
                   - f_7 * msg1_114[k]
                   + f_3 * pc_x[k] * msh_156[k];

        t_206[k] = pa_z[k] * lsi0_94[k]
                   - f_10 * pc_z[k] * lsi1_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pc_y, pc_z, lsi0_96, lsh_69, lsh_70, \
                         lsh_93, lsi1_96, msh_153, msh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_11 * lsh_69[k]
                   + f_3 * pc_z[k] * msh_153[k];

        t_208[k] = pa_z[k] * lsi0_96[k]
                   + f_12 * lsh_70[k]
                   - f_10 * pc_z[k] * lsi1_96[k];

        t_209[k] = f_12 * lsh_93[k]
                   + f_3 * pc_y[k] * msh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, lsh_161, lsh_162, lsh_163, lsh_164, \
                         msg0_119, msg1_119, msh_161, msh_162, msh_163, \
                         msh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_19 * lsh_161[k]
                   + f_4 * msg0_119[k]
                   - f_5 * msg1_119[k]
                   + f_3 * pc_x[k] * msh_161[k];

        t_211[k] = f_19 * lsh_162[k]
                   + f_3 * pc_x[k] * msh_162[k];

        t_212[k] = f_19 * lsh_163[k]
                   + f_3 * pc_x[k] * msh_163[k];

        t_213[k] = f_19 * lsh_164[k]
                   + f_3 * pc_x[k] * msh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, lsi0_105, lsh_165, \
                         lsh_166, lsh_167, lsi1_105, msh_165, msh_166, \
                         msh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_19 * lsh_165[k]
                   + f_3 * pc_x[k] * msh_165[k];

        t_215[k] = f_19 * lsh_166[k]
                   + f_3 * pc_x[k] * msh_166[k];

        t_216[k] = f_19 * lsh_167[k]
                   + f_3 * pc_x[k] * msh_167[k];

        t_217[k] = pa_z[k] * lsi0_105[k]
                   - f_10 * pc_z[k] * lsi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, lsh_78, lsh_101, lsh_102, msg0_117, \
                         msg0_118, msg1_117, msg1_118, msh_162, msh_164, \
                         msh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * lsh_78[k]
                   + f_3 * pc_z[k] * msh_162[k];

        t_219[k] = f_12 * lsh_101[k]
                   + f_8 * msg0_117[k]
                   - f_9 * msg1_117[k]
                   + f_3 * pc_y[k] * msh_164[k];

        t_220[k] = f_12 * lsh_102[k]
                   + f_6 * msg0_118[k]
                   - f_7 * msg1_118[k]
                   + f_3 * pc_y[k] * msh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pc_y, pc_z, lsi0_140, lsh_83, \
                         lsh_103, lsh_104, lsi1_140, msg0_119, msg1_119, msh_166, \
                         msh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * lsh_103[k]
                   + f_4 * msg0_119[k]
                   - f_5 * msg1_119[k]
                   + f_3 * pc_y[k] * msh_166[k];

        t_222[k] = f_12 * lsh_104[k]
                   + f_3 * pc_y[k] * msh_167[k];

        t_223[k] = f_11 * lsh_83[k]
                   + f_1 * msg0_119[k]
                   - f_2 * msg1_119[k]
                   + f_3 * pc_z[k] * msh_167[k];

        t_224[k] = pa_y[k] * lsi0_140[k]
                   - f_10 * pc_y[k] * lsi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, lsi0_143, lsh_84, \
                         lsh_105, lsh_106, lsh_107, lsi1_143, msh_168, \
                         msh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * lsh_105[k]
                   + f_3 * pc_y[k] * msh_168[k];

        t_226[k] = f_12 * lsh_84[k]
                   + f_3 * pc_z[k] * msh_168[k];

        t_227[k] = pa_y[k] * lsi0_143[k]
                   + f_12 * lsh_106[k]
                   - f_10 * pc_y[k] * lsi1_143[k];

        t_228[k] = f_11 * lsh_107[k]
                   + f_3 * pc_y[k] * msh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_y, pc_z, lsi0_145, lsi0_146, \
                         lsh_87, lsh_108, lsh_110, lsi1_145, lsi1_146, msh_171, \
                         msh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * lsi0_145[k]
                   - f_10 * pc_y[k] * lsi1_145[k];

        t_230[k] = pa_y[k] * lsi0_146[k]
                   + f_13 * lsh_108[k]
                   - f_10 * pc_y[k] * lsi1_146[k];

        t_231[k] = f_12 * lsh_87[k]
                   + f_3 * pc_z[k] * msh_171[k];

        t_232[k] = f_11 * lsh_110[k]
                   + f_3 * pc_y[k] * msh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_y, pc_z, lsi0_149, lsi0_150, lsh_90, \
                         lsh_111, lsi1_149, lsi1_150, msh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * lsi0_149[k]
                   - f_10 * pc_y[k] * lsi1_149[k];

        t_234[k] = pa_y[k] * lsi0_150[k]
                   + f_14 * lsh_111[k]
                   - f_10 * pc_y[k] * lsi1_150[k];

        t_235[k] = f_12 * lsh_90[k]
                   + f_3 * pc_z[k] * msh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pc_x, pc_y, lsi0_152, lsi0_154, \
                         lsh_113, lsh_114, lsh_183, lsi1_152, lsi1_154, msh_177, \
                         msh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * lsi0_152[k]
                   + f_12 * lsh_113[k]
                   - f_10 * pc_y[k] * lsi1_152[k];

        t_237[k] = f_11 * lsh_114[k]
                   + f_3 * pc_y[k] * msh_177[k];

        t_238[k] = pa_y[k] * lsi0_154[k]
                   - f_10 * pc_y[k] * lsi1_154[k];

        t_239[k] = f_19 * lsh_183[k]
                   + f_3 * pc_x[k] * msh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, lsh_184, lsh_185, lsh_186, \
                         lsh_187, lsh_188, msh_184, msh_185, msh_186, msh_187, \
                         msh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_19 * lsh_184[k]
                   + f_3 * pc_x[k] * msh_184[k];

        t_241[k] = f_19 * lsh_185[k]
                   + f_3 * pc_x[k] * msh_185[k];

        t_242[k] = f_19 * lsh_186[k]
                   + f_3 * pc_x[k] * msh_186[k];

        t_243[k] = f_19 * lsh_187[k]
                   + f_3 * pc_x[k] * msh_187[k];

        t_244[k] = f_19 * lsh_188[k]
                   + f_3 * pc_x[k] * msh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, lsh_99, lsh_120, lsh_122, msg0_130, \
                         msg0_132, msg1_130, msg1_132, msh_183, \
                         msh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * lsh_120[k]
                   + f_1 * msg0_130[k]
                   - f_2 * msg1_130[k]
                   + f_3 * pc_y[k] * msh_183[k];

        t_246[k] = f_12 * lsh_99[k]
                   + f_3 * pc_z[k] * msh_183[k];

        t_247[k] = f_11 * lsh_122[k]
                   + f_8 * msg0_132[k]
                   - f_9 * msg1_132[k]
                   + f_3 * pc_y[k] * msh_185[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *lsi0_167 = buffer.data(lsi0 + 167);
    const auto *lsi0_168 = buffer.data(lsi0 + 168);
    const auto *lsi0_171 = buffer.data(lsi0 + 171);
    const auto *lsi0_174 = buffer.data(lsi0 + 174);
    const auto *lsi0_178 = buffer.data(lsi0 + 178);
    const auto *lsi0_180 = buffer.data(lsi0 + 180);
    const auto *lsi0_189 = buffer.data(lsi0 + 189);

    const auto *lsh_105 = buffer.data(lsh + 105);
    const auto *lsh_123 = buffer.data(lsh + 123);
    const auto *lsh_124 = buffer.data(lsh + 124);
    const auto *lsh_125 = buffer.data(lsh + 125);
    const auto *lsh_126 = buffer.data(lsh + 126);
    const auto *lsh_129 = buffer.data(lsh + 129);
    const auto *lsh_131 = buffer.data(lsh + 131);
    const auto *lsh_132 = buffer.data(lsh + 132);
    const auto *lsh_133 = buffer.data(lsh + 133);
    const auto *lsh_135 = buffer.data(lsh + 135);
    const auto *lsh_141 = buffer.data(lsh + 141);
    const auto *lsh_146 = buffer.data(lsh + 146);
    const auto *lsh_147 = buffer.data(lsh + 147);
    const auto *lsh_149 = buffer.data(lsh + 149);
    const auto *lsh_150 = buffer.data(lsh + 150);
    const auto *lsh_152 = buffer.data(lsh + 152);
    const auto *lsh_153 = buffer.data(lsh + 153);
    const auto *lsh_156 = buffer.data(lsh + 156);
    const auto *lsh_162 = buffer.data(lsh + 162);
    const auto *lsh_164 = buffer.data(lsh + 164);
    const auto *lsh_165 = buffer.data(lsh + 165);
    const auto *lsh_166 = buffer.data(lsh + 166);
    const auto *lsh_167 = buffer.data(lsh + 167);
    const auto *lsh_168 = buffer.data(lsh + 168);
    const auto *lsh_170 = buffer.data(lsh + 170);
    const auto *lsh_173 = buffer.data(lsh + 173);
    const auto *lsh_177 = buffer.data(lsh + 177);
    const auto *lsh_183 = buffer.data(lsh + 183);
    const auto *lsh_185 = buffer.data(lsh + 185);
    const auto *lsh_186 = buffer.data(lsh + 186);
    const auto *lsh_187 = buffer.data(lsh + 187);
    const auto *lsh_189 = buffer.data(lsh + 189);
    const auto *lsh_194 = buffer.data(lsh + 194);
    const auto *lsh_198 = buffer.data(lsh + 198);
    const auto *lsh_203 = buffer.data(lsh + 203);
    const auto *lsh_204 = buffer.data(lsh + 204);
    const auto *lsh_205 = buffer.data(lsh + 205);
    const auto *lsh_206 = buffer.data(lsh + 206);
    const auto *lsh_207 = buffer.data(lsh + 207);
    const auto *lsh_209 = buffer.data(lsh + 209);
    const auto *lsh_210 = buffer.data(lsh + 210);
    const auto *lsh_213 = buffer.data(lsh + 213);
    const auto *lsh_216 = buffer.data(lsh + 216);
    const auto *lsh_220 = buffer.data(lsh + 220);
    const auto *lsh_225 = buffer.data(lsh + 225);
    const auto *lsh_227 = buffer.data(lsh + 227);
    const auto *lsh_228 = buffer.data(lsh + 228);
    const auto *lsh_229 = buffer.data(lsh + 229);
    const auto *lsh_230 = buffer.data(lsh + 230);
    const auto *lsh_236 = buffer.data(lsh + 236);
    const auto *lsh_240 = buffer.data(lsh + 240);
    const auto *lsh_245 = buffer.data(lsh + 245);
    const auto *lsh_246 = buffer.data(lsh + 246);
    const auto *lsh_247 = buffer.data(lsh + 247);
    const auto *lsh_248 = buffer.data(lsh + 248);
    const auto *lsh_249 = buffer.data(lsh + 249);
    const auto *lsh_250 = buffer.data(lsh + 250);
    const auto *lsh_251 = buffer.data(lsh + 251);
    const auto *lsh_252 = buffer.data(lsh + 252);
    const auto *lsh_255 = buffer.data(lsh + 255);
    const auto *lsh_257 = buffer.data(lsh + 257);
    const auto *lsh_258 = buffer.data(lsh + 258);
    const auto *lsh_261 = buffer.data(lsh + 261);
    const auto *lsh_262 = buffer.data(lsh + 262);
    const auto *lsh_264 = buffer.data(lsh + 264);
    const auto *lsh_266 = buffer.data(lsh + 266);
    const auto *lsh_267 = buffer.data(lsh + 267);
    const auto *lsh_268 = buffer.data(lsh + 268);
    const auto *lsh_269 = buffer.data(lsh + 269);
    const auto *lsh_270 = buffer.data(lsh + 270);
    const auto *lsh_271 = buffer.data(lsh + 271);
    const auto *lsh_272 = buffer.data(lsh + 272);

    const auto *lsi1_167 = buffer.data(lsi1 + 167);
    const auto *lsi1_168 = buffer.data(lsi1 + 168);
    const auto *lsi1_171 = buffer.data(lsi1 + 171);
    const auto *lsi1_174 = buffer.data(lsi1 + 174);
    const auto *lsi1_178 = buffer.data(lsi1 + 178);
    const auto *lsi1_180 = buffer.data(lsi1 + 180);
    const auto *lsi1_189 = buffer.data(lsi1 + 189);

    const auto *msg0_133 = buffer.data(msg0 + 133);
    const auto *msg0_134 = buffer.data(msg0 + 134);
    const auto *msg0_135 = buffer.data(msg0 + 135);
    const auto *msg0_136 = buffer.data(msg0 + 136);
    const auto *msg0_137 = buffer.data(msg0 + 137);
    const auto *msg0_138 = buffer.data(msg0 + 138);
    const auto *msg0_139 = buffer.data(msg0 + 139);
    const auto *msg0_140 = buffer.data(msg0 + 140);
    const auto *msg0_144 = buffer.data(msg0 + 144);
    const auto *msg0_145 = buffer.data(msg0 + 145);
    const auto *msg0_146 = buffer.data(msg0 + 146);
    const auto *msg0_147 = buffer.data(msg0 + 147);
    const auto *msg0_148 = buffer.data(msg0 + 148);
    const auto *msg0_149 = buffer.data(msg0 + 149);
    const auto *msg0_150 = buffer.data(msg0 + 150);
    const auto *msg0_152 = buffer.data(msg0 + 152);
    const auto *msg0_153 = buffer.data(msg0 + 153);
    const auto *msg0_155 = buffer.data(msg0 + 155);
    const auto *msg0_156 = buffer.data(msg0 + 156);
    const auto *msg0_160 = buffer.data(msg0 + 160);
    const auto *msg0_161 = buffer.data(msg0 + 161);
    const auto *msg0_162 = buffer.data(msg0 + 162);
    const auto *msg0_164 = buffer.data(msg0 + 164);
    const auto *msg0_170 = buffer.data(msg0 + 170);
    const auto *msg0_174 = buffer.data(msg0 + 174);
    const auto *msg0_177 = buffer.data(msg0 + 177);
    const auto *msg0_178 = buffer.data(msg0 + 178);
    const auto *msg0_179 = buffer.data(msg0 + 179);
    const auto *msg0_180 = buffer.data(msg0 + 180);
    const auto *msg0_183 = buffer.data(msg0 + 183);
    const auto *msg0_185 = buffer.data(msg0 + 185);
    const auto *msg0_186 = buffer.data(msg0 + 186);
    const auto *msg0_189 = buffer.data(msg0 + 189);
    const auto *msg0_190 = buffer.data(msg0 + 190);
    const auto *msg0_192 = buffer.data(msg0 + 192);
    const auto *msg0_193 = buffer.data(msg0 + 193);
    const auto *msg0_194 = buffer.data(msg0 + 194);

    const auto *msg1_133 = buffer.data(msg1 + 133);
    const auto *msg1_134 = buffer.data(msg1 + 134);
    const auto *msg1_135 = buffer.data(msg1 + 135);
    const auto *msg1_136 = buffer.data(msg1 + 136);
    const auto *msg1_137 = buffer.data(msg1 + 137);
    const auto *msg1_138 = buffer.data(msg1 + 138);
    const auto *msg1_139 = buffer.data(msg1 + 139);
    const auto *msg1_140 = buffer.data(msg1 + 140);
    const auto *msg1_144 = buffer.data(msg1 + 144);
    const auto *msg1_145 = buffer.data(msg1 + 145);
    const auto *msg1_146 = buffer.data(msg1 + 146);
    const auto *msg1_147 = buffer.data(msg1 + 147);
    const auto *msg1_148 = buffer.data(msg1 + 148);
    const auto *msg1_149 = buffer.data(msg1 + 149);
    const auto *msg1_150 = buffer.data(msg1 + 150);
    const auto *msg1_152 = buffer.data(msg1 + 152);
    const auto *msg1_153 = buffer.data(msg1 + 153);
    const auto *msg1_155 = buffer.data(msg1 + 155);
    const auto *msg1_156 = buffer.data(msg1 + 156);
    const auto *msg1_160 = buffer.data(msg1 + 160);
    const auto *msg1_161 = buffer.data(msg1 + 161);
    const auto *msg1_162 = buffer.data(msg1 + 162);
    const auto *msg1_164 = buffer.data(msg1 + 164);
    const auto *msg1_170 = buffer.data(msg1 + 170);
    const auto *msg1_174 = buffer.data(msg1 + 174);
    const auto *msg1_177 = buffer.data(msg1 + 177);
    const auto *msg1_178 = buffer.data(msg1 + 178);
    const auto *msg1_179 = buffer.data(msg1 + 179);
    const auto *msg1_180 = buffer.data(msg1 + 180);
    const auto *msg1_183 = buffer.data(msg1 + 183);
    const auto *msg1_185 = buffer.data(msg1 + 185);
    const auto *msg1_186 = buffer.data(msg1 + 186);
    const auto *msg1_189 = buffer.data(msg1 + 189);
    const auto *msg1_190 = buffer.data(msg1 + 190);
    const auto *msg1_192 = buffer.data(msg1 + 192);
    const auto *msg1_193 = buffer.data(msg1 + 193);
    const auto *msg1_194 = buffer.data(msg1 + 194);

    const auto *msh_186 = buffer.data(msh + 186);
    const auto *msh_187 = buffer.data(msh + 187);
    const auto *msh_188 = buffer.data(msh + 188);
    const auto *msh_189 = buffer.data(msh + 189);
    const auto *msh_190 = buffer.data(msh + 190);
    const auto *msh_191 = buffer.data(msh + 191);
    const auto *msh_192 = buffer.data(msh + 192);
    const auto *msh_193 = buffer.data(msh + 193);
    const auto *msh_194 = buffer.data(msh + 194);
    const auto *msh_195 = buffer.data(msh + 195);
    const auto *msh_196 = buffer.data(msh + 196);
    const auto *msh_197 = buffer.data(msh + 197);
    const auto *msh_198 = buffer.data(msh + 198);
    const auto *msh_203 = buffer.data(msh + 203);
    const auto *msh_204 = buffer.data(msh + 204);
    const auto *msh_205 = buffer.data(msh + 205);
    const auto *msh_206 = buffer.data(msh + 206);
    const auto *msh_207 = buffer.data(msh + 207);
    const auto *msh_208 = buffer.data(msh + 208);
    const auto *msh_209 = buffer.data(msh + 209);
    const auto *msh_210 = buffer.data(msh + 210);
    const auto *msh_211 = buffer.data(msh + 211);
    const auto *msh_212 = buffer.data(msh + 212);
    const auto *msh_213 = buffer.data(msh + 213);
    const auto *msh_215 = buffer.data(msh + 215);
    const auto *msh_216 = buffer.data(msh + 216);
    const auto *msh_217 = buffer.data(msh + 217);
    const auto *msh_219 = buffer.data(msh + 219);
    const auto *msh_220 = buffer.data(msh + 220);
    const auto *msh_225 = buffer.data(msh + 225);
    const auto *msh_226 = buffer.data(msh + 226);
    const auto *msh_227 = buffer.data(msh + 227);
    const auto *msh_228 = buffer.data(msh + 228);
    const auto *msh_229 = buffer.data(msh + 229);
    const auto *msh_230 = buffer.data(msh + 230);
    const auto *msh_231 = buffer.data(msh + 231);
    const auto *msh_233 = buffer.data(msh + 233);
    const auto *msh_234 = buffer.data(msh + 234);
    const auto *msh_236 = buffer.data(msh + 236);
    const auto *msh_237 = buffer.data(msh + 237);
    const auto *msh_240 = buffer.data(msh + 240);
    const auto *msh_245 = buffer.data(msh + 245);
    const auto *msh_246 = buffer.data(msh + 246);
    const auto *msh_247 = buffer.data(msh + 247);
    const auto *msh_248 = buffer.data(msh + 248);
    const auto *msh_249 = buffer.data(msh + 249);
    const auto *msh_250 = buffer.data(msh + 250);
    const auto *msh_251 = buffer.data(msh + 251);
    const auto *msh_252 = buffer.data(msh + 252);
    const auto *msh_254 = buffer.data(msh + 254);
    const auto *msh_255 = buffer.data(msh + 255);
    const auto *msh_257 = buffer.data(msh + 257);
    const auto *msh_258 = buffer.data(msh + 258);
    const auto *msh_261 = buffer.data(msh + 261);
    const auto *msh_262 = buffer.data(msh + 262);
    const auto *msh_264 = buffer.data(msh + 264);
    const auto *msh_266 = buffer.data(msh + 266);
    const auto *msh_267 = buffer.data(msh + 267);
    const auto *msh_268 = buffer.data(msh + 268);
    const auto *msh_269 = buffer.data(msh + 269);
    const auto *msh_270 = buffer.data(msh + 270);
    const auto *msh_271 = buffer.data(msh + 271);
    const auto *msh_272 = buffer.data(msh + 272);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, lsh_123, lsh_124, lsh_125, msg0_133, \
                         msg0_134, msg1_133, msg1_134, msh_186, msh_187, \
                         msh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * lsh_123[k]
                   + f_6 * msg0_133[k]
                   - f_7 * msg1_133[k]
                   + f_3 * pc_y[k] * msh_186[k];

        t_249[k] = f_11 * lsh_124[k]
                   + f_4 * msg0_134[k]
                   - f_5 * msg1_134[k]
                   + f_3 * pc_y[k] * msh_187[k];

        t_250[k] = f_11 * lsh_125[k]
                   + f_3 * pc_y[k] * msh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pc_x, pc_y, pc_z, lsi0_167, \
                         lsh_105, lsh_189, lsi1_167, msg0_135, msg1_135, \
                         msh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * lsi0_167[k]
                   - f_10 * pc_y[k] * lsi1_167[k];

        t_252[k] = f_19 * lsh_189[k]
                   + f_1 * msg0_135[k]
                   - f_2 * msg1_135[k]
                   + f_3 * pc_x[k] * msh_189[k];

        t_253[k] = f_3 * pc_y[k] * msh_189[k];

        t_254[k] = f_13 * lsh_105[k]
                   + f_3 * pc_z[k] * msh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, lsh_194, msg0_135, msg0_140, \
                         msg1_135, msg1_140, msh_190, msh_191, \
                         msh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * msg0_135[k]
                   - f_5 * msg1_135[k]
                   + f_3 * pc_y[k] * msh_190[k];

        t_256[k] = f_3 * pc_y[k] * msh_191[k];

        t_257[k] = f_19 * lsh_194[k]
                   + f_8 * msg0_140[k]
                   - f_9 * msg1_140[k]
                   + f_3 * pc_x[k] * msh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_y, msg0_136, msg0_137, msg1_136, msg1_137, \
                         msh_192, msh_193, msh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_6 * msg0_136[k]
                   - f_7 * msg1_136[k]
                   + f_3 * pc_y[k] * msh_192[k];

        t_259[k] = f_4 * msg0_137[k]
                   - f_5 * msg1_137[k]
                   + f_3 * pc_y[k] * msh_193[k];

        t_260[k] = f_3 * pc_y[k] * msh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_y, lsh_198, msg0_138, msg0_139, \
                         msg0_144, msg1_138, msg1_139, msg1_144, msh_195, msh_196, \
                         msh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_19 * lsh_198[k]
                   + f_6 * msg0_144[k]
                   - f_7 * msg1_144[k]
                   + f_3 * pc_x[k] * msh_198[k];

        t_262[k] = f_8 * msg0_138[k]
                   - f_9 * msg1_138[k]
                   + f_3 * pc_y[k] * msh_195[k];

        t_263[k] = f_6 * msg0_139[k]
                   - f_7 * msg1_139[k]
                   + f_3 * pc_y[k] * msh_196[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, lsh_203, lsh_204, msg0_140, \
                         msg0_149, msg1_140, msg1_149, msh_197, msh_198, msh_203, \
                         msh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * msg0_140[k]
                   - f_5 * msg1_140[k]
                   + f_3 * pc_y[k] * msh_197[k];

        t_265[k] = f_3 * pc_y[k] * msh_198[k];

        t_266[k] = f_19 * lsh_203[k]
                   + f_4 * msg0_149[k]
                   - f_5 * msg1_149[k]
                   + f_3 * pc_x[k] * msh_203[k];

        t_267[k] = f_19 * lsh_204[k]
                   + f_3 * pc_x[k] * msh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pc_x, pc_y, lsh_205, lsh_206, \
                         lsh_207, lsh_209, msh_203, msh_205, msh_206, msh_207, \
                         msh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_19 * lsh_205[k]
                   + f_3 * pc_x[k] * msh_205[k];

        t_269[k] = f_19 * lsh_206[k]
                   + f_3 * pc_x[k] * msh_206[k];

        t_270[k] = f_19 * lsh_207[k]
                   + f_3 * pc_x[k] * msh_207[k];

        t_271[k] = f_3 * pc_y[k] * msh_203[k];

        t_272[k] = f_19 * lsh_209[k]
                   + f_3 * pc_x[k] * msh_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_y, msg0_145, msg0_146, msg0_147, msg1_145, \
                         msg1_146, msg1_147, msh_204, msh_205, \
                         msh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * msg0_145[k]
                   - f_2 * msg1_145[k]
                   + f_3 * pc_y[k] * msh_204[k];

        t_274[k] = f_16 * msg0_146[k]
                   - f_17 * msg1_146[k]
                   + f_3 * pc_y[k] * msh_205[k];

        t_275[k] = f_8 * msg0_147[k]
                   - f_9 * msg1_147[k]
                   + f_3 * pc_y[k] * msh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, lsh_125, msg0_148, msg0_149, \
                         msg1_148, msg1_149, msh_207, msh_208, \
                         msh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * msg0_148[k]
                   - f_7 * msg1_148[k]
                   + f_3 * pc_y[k] * msh_207[k];

        t_277[k] = f_4 * msg0_149[k]
                   - f_5 * msg1_149[k]
                   + f_3 * pc_y[k] * msh_208[k];

        t_278[k] = f_3 * pc_y[k] * msh_209[k];

        t_279[k] = f_13 * lsh_125[k]
                   + f_1 * msg0_149[k]
                   - f_2 * msg1_149[k]
                   + f_3 * pc_z[k] * msh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, lsh_126, lsh_210, \
                         lsh_213, msg0_150, msg0_153, msg1_150, msg1_153, msh_210, \
                         msh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_20 * lsh_210[k]
                   + f_1 * msg0_150[k]
                   - f_2 * msg1_150[k]
                   + f_3 * pc_x[k] * msh_210[k];

        t_281[k] = f_14 * lsh_126[k]
                   + f_3 * pc_y[k] * msh_210[k];

        t_282[k] = f_3 * pc_z[k] * msh_210[k];

        t_283[k] = f_20 * lsh_213[k]
                   + f_8 * msg0_153[k]
                   - f_9 * msg1_153[k]
                   + f_3 * pc_x[k] * msh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_z, lsh_216, msg0_150, msg0_156, \
                         msg1_150, msg1_156, msh_211, msh_212, msh_213, \
                         msh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * msh_211[k];

        t_285[k] = f_4 * msg0_150[k]
                   - f_5 * msg1_150[k]
                   + f_3 * pc_z[k] * msh_212[k];

        t_286[k] = f_20 * lsh_216[k]
                   + f_6 * msg0_156[k]
                   - f_7 * msg1_156[k]
                   + f_3 * pc_x[k] * msh_216[k];

        t_287[k] = f_3 * pc_z[k] * msh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, pc_z, lsh_131, lsh_220, \
                         msg0_152, msg0_160, msg1_152, msg1_160, msh_215, msh_216, \
                         msh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_14 * lsh_131[k]
                   + f_3 * pc_y[k] * msh_215[k];

        t_289[k] = f_6 * msg0_152[k]
                   - f_7 * msg1_152[k]
                   + f_3 * pc_z[k] * msh_215[k];

        t_290[k] = f_20 * lsh_220[k]
                   + f_4 * msg0_160[k]
                   - f_5 * msg1_160[k]
                   + f_3 * pc_x[k] * msh_220[k];

        t_291[k] = f_3 * pc_z[k] * msh_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, lsh_135, lsh_225, \
                         msg0_153, msg0_155, msg1_153, msg1_155, msh_217, msh_219, \
                         msh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * msg0_153[k]
                   - f_5 * msg1_153[k]
                   + f_3 * pc_z[k] * msh_217[k];

        t_293[k] = f_14 * lsh_135[k]
                   + f_3 * pc_y[k] * msh_219[k];

        t_294[k] = f_8 * msg0_155[k]
                   - f_9 * msg1_155[k]
                   + f_3 * pc_z[k] * msh_219[k];

        t_295[k] = f_20 * lsh_225[k]
                   + f_3 * pc_x[k] * msh_225[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, pc_z, lsh_227, lsh_228, \
                         lsh_229, lsh_230, msh_220, msh_227, msh_228, msh_229, \
                         msh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_3 * pc_z[k] * msh_220[k];

        t_297[k] = f_20 * lsh_227[k]
                   + f_3 * pc_x[k] * msh_227[k];

        t_298[k] = f_20 * lsh_228[k]
                   + f_3 * pc_x[k] * msh_228[k];

        t_299[k] = f_20 * lsh_229[k]
                   + f_3 * pc_x[k] * msh_229[k];

        t_300[k] = f_20 * lsh_230[k]
                   + f_3 * pc_x[k] * msh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_y, pc_z, lsh_141, msg0_160, msg0_161, \
                         msg1_160, msg1_161, msh_225, msh_226, \
                         msh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * lsh_141[k]
                   + f_1 * msg0_160[k]
                   - f_2 * msg1_160[k]
                   + f_3 * pc_y[k] * msh_225[k];

        t_302[k] = f_3 * pc_z[k] * msh_225[k];

        t_303[k] = f_4 * msg0_160[k]
                   - f_5 * msg1_160[k]
                   + f_3 * pc_z[k] * msh_226[k];

        t_304[k] = f_6 * msg0_161[k]
                   - f_7 * msg1_161[k]
                   + f_3 * pc_z[k] * msh_227[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pc_y, pc_z, lsi0_168, lsh_146, \
                         lsi1_168, msg0_162, msg0_164, msg1_162, msg1_164, msh_228, \
                         msh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_8 * msg0_162[k]
                   - f_9 * msg1_162[k]
                   + f_3 * pc_z[k] * msh_228[k];

        t_306[k] = f_14 * lsh_146[k]
                   + f_3 * pc_y[k] * msh_230[k];

        t_307[k] = f_1 * msg0_164[k]
                   - f_2 * msg1_164[k]
                   + f_3 * pc_z[k] * msh_230[k];

        t_308[k] = pa_z[k] * lsi0_168[k]
                   - f_10 * pc_z[k] * lsi1_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_z, pc_y, pc_z, lsi0_171, lsh_126, \
                         lsh_147, lsh_149, lsi1_171, msh_231, msh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_13 * lsh_147[k]
                   + f_3 * pc_y[k] * msh_231[k];

        t_310[k] = f_11 * lsh_126[k]
                   + f_3 * pc_z[k] * msh_231[k];

        t_311[k] = pa_z[k] * lsi0_171[k]
                   - f_10 * pc_z[k] * lsi1_171[k];

        t_312[k] = f_13 * lsh_149[k]
                   + f_3 * pc_y[k] * msh_233[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_z, pc_x, pc_z, lsi0_174, lsh_129, lsh_236, \
                         lsi1_174, msg0_170, msg1_170, msh_234, \
                         msh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_20 * lsh_236[k]
                   + f_8 * msg0_170[k]
                   - f_9 * msg1_170[k]
                   + f_3 * pc_x[k] * msh_236[k];

        t_314[k] = pa_z[k] * lsi0_174[k]
                   - f_10 * pc_z[k] * lsi1_174[k];

        t_315[k] = f_11 * lsh_129[k]
                   + f_3 * pc_z[k] * msh_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_z, pc_x, pc_y, pc_z, lsi0_178, lsh_152, \
                         lsh_240, lsi1_178, msg0_174, msg1_174, msh_236, \
                         msh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * lsh_152[k]
                   + f_3 * pc_y[k] * msh_236[k];

        t_317[k] = f_20 * lsh_240[k]
                   + f_6 * msg0_174[k]
                   - f_7 * msg1_174[k]
                   + f_3 * pc_x[k] * msh_240[k];

        t_318[k] = pa_z[k] * lsi0_178[k]
                   - f_10 * pc_z[k] * lsi1_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_y, pc_z, lsi0_180, lsh_132, lsh_133, \
                         lsh_156, lsi1_180, msh_237, msh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_11 * lsh_132[k]
                   + f_3 * pc_z[k] * msh_237[k];

        t_320[k] = pa_z[k] * lsi0_180[k]
                   + f_12 * lsh_133[k]
                   - f_10 * pc_z[k] * lsi1_180[k];

        t_321[k] = f_13 * lsh_156[k]
                   + f_3 * pc_y[k] * msh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, lsh_245, lsh_246, lsh_247, lsh_248, \
                         msg0_179, msg1_179, msh_245, msh_246, msh_247, \
                         msh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_20 * lsh_245[k]
                   + f_4 * msg0_179[k]
                   - f_5 * msg1_179[k]
                   + f_3 * pc_x[k] * msh_245[k];

        t_323[k] = f_20 * lsh_246[k]
                   + f_3 * pc_x[k] * msh_246[k];

        t_324[k] = f_20 * lsh_247[k]
                   + f_3 * pc_x[k] * msh_247[k];

        t_325[k] = f_20 * lsh_248[k]
                   + f_3 * pc_x[k] * msh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pa_z, pc_x, pc_z, lsi0_189, lsh_249, \
                         lsh_250, lsh_251, lsi1_189, msh_249, msh_250, \
                         msh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_20 * lsh_249[k]
                   + f_3 * pc_x[k] * msh_249[k];

        t_327[k] = f_20 * lsh_250[k]
                   + f_3 * pc_x[k] * msh_250[k];

        t_328[k] = f_20 * lsh_251[k]
                   + f_3 * pc_x[k] * msh_251[k];

        t_329[k] = pa_z[k] * lsi0_189[k]
                   - f_10 * pc_z[k] * lsi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, lsh_141, lsh_164, lsh_165, msg0_177, \
                         msg0_178, msg1_177, msg1_178, msh_246, msh_248, \
                         msh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * lsh_141[k]
                   + f_3 * pc_z[k] * msh_246[k];

        t_331[k] = f_13 * lsh_164[k]
                   + f_8 * msg0_177[k]
                   - f_9 * msg1_177[k]
                   + f_3 * pc_y[k] * msh_248[k];

        t_332[k] = f_13 * lsh_165[k]
                   + f_6 * msg0_178[k]
                   - f_7 * msg1_178[k]
                   + f_3 * pc_y[k] * msh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, lsh_146, lsh_166, lsh_167, msg0_179, \
                         msg1_179, msh_250, msh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * lsh_166[k]
                   + f_4 * msg0_179[k]
                   - f_5 * msg1_179[k]
                   + f_3 * pc_y[k] * msh_250[k];

        t_334[k] = f_13 * lsh_167[k]
                   + f_3 * pc_y[k] * msh_251[k];

        t_335[k] = f_11 * lsh_146[k]
                   + f_1 * msg0_179[k]
                   - f_2 * msg1_179[k]
                   + f_3 * pc_z[k] * msh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, lsh_147, lsh_168, lsh_252, \
                         msg0_180, msg1_180, msh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_20 * lsh_252[k]
                   + f_1 * msg0_180[k]
                   - f_2 * msg1_180[k]
                   + f_3 * pc_x[k] * msh_252[k];

        t_337[k] = f_12 * lsh_168[k]
                   + f_3 * pc_y[k] * msh_252[k];

        t_338[k] = f_12 * lsh_147[k]
                   + f_3 * pc_z[k] * msh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, lsh_170, lsh_255, lsh_257, msg0_183, \
                         msg0_185, msg1_183, msg1_185, msh_254, msh_255, \
                         msh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_20 * lsh_255[k]
                   + f_8 * msg0_183[k]
                   - f_9 * msg1_183[k]
                   + f_3 * pc_x[k] * msh_255[k];

        t_340[k] = f_12 * lsh_170[k]
                   + f_3 * pc_y[k] * msh_254[k];

        t_341[k] = f_20 * lsh_257[k]
                   + f_8 * msg0_185[k]
                   - f_9 * msg1_185[k]
                   + f_3 * pc_x[k] * msh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, lsh_150, lsh_173, lsh_258, \
                         msg0_186, msg1_186, msh_255, msh_257, \
                         msh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_20 * lsh_258[k]
                   + f_6 * msg0_186[k]
                   - f_7 * msg1_186[k]
                   + f_3 * pc_x[k] * msh_258[k];

        t_343[k] = f_12 * lsh_150[k]
                   + f_3 * pc_z[k] * msh_255[k];

        t_344[k] = f_12 * lsh_173[k]
                   + f_3 * pc_y[k] * msh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, lsh_153, lsh_261, lsh_262, msg0_189, \
                         msg0_190, msg1_189, msg1_190, msh_258, msh_261, \
                         msh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_20 * lsh_261[k]
                   + f_6 * msg0_189[k]
                   - f_7 * msg1_189[k]
                   + f_3 * pc_x[k] * msh_261[k];

        t_346[k] = f_20 * lsh_262[k]
                   + f_4 * msg0_190[k]
                   - f_5 * msg1_190[k]
                   + f_3 * pc_x[k] * msh_262[k];

        t_347[k] = f_12 * lsh_153[k]
                   + f_3 * pc_z[k] * msh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, lsh_177, lsh_264, lsh_266, msg0_192, \
                         msg0_194, msg1_192, msg1_194, msh_261, msh_264, \
                         msh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_20 * lsh_264[k]
                   + f_4 * msg0_192[k]
                   - f_5 * msg1_192[k]
                   + f_3 * pc_x[k] * msh_264[k];

        t_349[k] = f_12 * lsh_177[k]
                   + f_3 * pc_y[k] * msh_261[k];

        t_350[k] = f_20 * lsh_266[k]
                   + f_4 * msg0_194[k]
                   - f_5 * msg1_194[k]
                   + f_3 * pc_x[k] * msh_266[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, lsh_267, lsh_268, lsh_269, \
                         lsh_270, lsh_271, msh_267, msh_268, msh_269, msh_270, \
                         msh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_20 * lsh_267[k]
                   + f_3 * pc_x[k] * msh_267[k];

        t_352[k] = f_20 * lsh_268[k]
                   + f_3 * pc_x[k] * msh_268[k];

        t_353[k] = f_20 * lsh_269[k]
                   + f_3 * pc_x[k] * msh_269[k];

        t_354[k] = f_20 * lsh_270[k]
                   + f_3 * pc_x[k] * msh_270[k];

        t_355[k] = f_20 * lsh_271[k]
                   + f_3 * pc_x[k] * msh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, lsh_162, lsh_183, lsh_272, \
                         msg0_190, msg1_190, msh_267, msh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_20 * lsh_272[k]
                   + f_3 * pc_x[k] * msh_272[k];

        t_357[k] = f_12 * lsh_183[k]
                   + f_1 * msg0_190[k]
                   - f_2 * msg1_190[k]
                   + f_3 * pc_y[k] * msh_267[k];

        t_358[k] = f_12 * lsh_162[k]
                   + f_3 * pc_z[k] * msh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, lsh_185, lsh_186, lsh_187, msg0_192, \
                         msg0_193, msg0_194, msg1_192, msg1_193, msg1_194, msh_269, msh_270, \
                         msh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * lsh_185[k]
                   + f_8 * msg0_192[k]
                   - f_9 * msg1_192[k]
                   + f_3 * pc_y[k] * msh_269[k];

        t_360[k] = f_12 * lsh_186[k]
                   + f_6 * msg0_193[k]
                   - f_7 * msg1_193[k]
                   + f_3 * pc_y[k] * msh_270[k];

        t_361[k] = f_12 * lsh_187[k]
                   + f_4 * msg0_194[k]
                   - f_5 * msg1_194[k]
                   + f_3 * pc_y[k] * msh_271[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_20 = 2.5 / q;

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

    const auto *lsi0_252 = buffer.data(lsi0 + 252);
    const auto *lsi0_255 = buffer.data(lsi0 + 255);
    const auto *lsi0_257 = buffer.data(lsi0 + 257);
    const auto *lsi0_258 = buffer.data(lsi0 + 258);
    const auto *lsi0_261 = buffer.data(lsi0 + 261);
    const auto *lsi0_262 = buffer.data(lsi0 + 262);
    const auto *lsi0_264 = buffer.data(lsi0 + 264);
    const auto *lsi0_266 = buffer.data(lsi0 + 266);
    const auto *lsi0_279 = buffer.data(lsi0 + 279);
    const auto *lsi0_280 = buffer.data(lsi0 + 280);
    const auto *lsi0_283 = buffer.data(lsi0 + 283);
    const auto *lsi0_286 = buffer.data(lsi0 + 286);
    const auto *lsi0_290 = buffer.data(lsi0 + 290);
    const auto *lsi0_292 = buffer.data(lsi0 + 292);
    const auto *lsi0_301 = buffer.data(lsi0 + 301);

    const auto *lsh_167 = buffer.data(lsh + 167);
    const auto *lsh_168 = buffer.data(lsh + 168);
    const auto *lsh_171 = buffer.data(lsh + 171);
    const auto *lsh_174 = buffer.data(lsh + 174);
    const auto *lsh_183 = buffer.data(lsh + 183);
    const auto *lsh_188 = buffer.data(lsh + 188);
    const auto *lsh_189 = buffer.data(lsh + 189);
    const auto *lsh_190 = buffer.data(lsh + 190);
    const auto *lsh_191 = buffer.data(lsh + 191);
    const auto *lsh_192 = buffer.data(lsh + 192);
    const auto *lsh_194 = buffer.data(lsh + 194);
    const auto *lsh_195 = buffer.data(lsh + 195);
    const auto *lsh_197 = buffer.data(lsh + 197);
    const auto *lsh_198 = buffer.data(lsh + 198);
    const auto *lsh_204 = buffer.data(lsh + 204);
    const auto *lsh_206 = buffer.data(lsh + 206);
    const auto *lsh_207 = buffer.data(lsh + 207);
    const auto *lsh_208 = buffer.data(lsh + 208);
    const auto *lsh_209 = buffer.data(lsh + 209);
    const auto *lsh_210 = buffer.data(lsh + 210);
    const auto *lsh_213 = buffer.data(lsh + 213);
    const auto *lsh_215 = buffer.data(lsh + 215);
    const auto *lsh_216 = buffer.data(lsh + 216);
    const auto *lsh_217 = buffer.data(lsh + 217);
    const auto *lsh_219 = buffer.data(lsh + 219);
    const auto *lsh_225 = buffer.data(lsh + 225);
    const auto *lsh_230 = buffer.data(lsh + 230);
    const auto *lsh_231 = buffer.data(lsh + 231);
    const auto *lsh_233 = buffer.data(lsh + 233);
    const auto *lsh_236 = buffer.data(lsh + 236);
    const auto *lsh_240 = buffer.data(lsh + 240);
    const auto *lsh_248 = buffer.data(lsh + 248);
    const auto *lsh_249 = buffer.data(lsh + 249);
    const auto *lsh_250 = buffer.data(lsh + 250);
    const auto *lsh_251 = buffer.data(lsh + 251);
    const auto *lsh_252 = buffer.data(lsh + 252);
    const auto *lsh_288 = buffer.data(lsh + 288);
    const auto *lsh_289 = buffer.data(lsh + 289);
    const auto *lsh_290 = buffer.data(lsh + 290);
    const auto *lsh_291 = buffer.data(lsh + 291);
    const auto *lsh_292 = buffer.data(lsh + 292);
    const auto *lsh_293 = buffer.data(lsh + 293);
    const auto *lsh_294 = buffer.data(lsh + 294);
    const auto *lsh_299 = buffer.data(lsh + 299);
    const auto *lsh_303 = buffer.data(lsh + 303);
    const auto *lsh_308 = buffer.data(lsh + 308);
    const auto *lsh_309 = buffer.data(lsh + 309);
    const auto *lsh_310 = buffer.data(lsh + 310);
    const auto *lsh_311 = buffer.data(lsh + 311);
    const auto *lsh_312 = buffer.data(lsh + 312);
    const auto *lsh_314 = buffer.data(lsh + 314);
    const auto *lsh_315 = buffer.data(lsh + 315);
    const auto *lsh_318 = buffer.data(lsh + 318);
    const auto *lsh_321 = buffer.data(lsh + 321);
    const auto *lsh_325 = buffer.data(lsh + 325);
    const auto *lsh_330 = buffer.data(lsh + 330);
    const auto *lsh_332 = buffer.data(lsh + 332);
    const auto *lsh_333 = buffer.data(lsh + 333);
    const auto *lsh_334 = buffer.data(lsh + 334);
    const auto *lsh_335 = buffer.data(lsh + 335);
    const auto *lsh_341 = buffer.data(lsh + 341);
    const auto *lsh_345 = buffer.data(lsh + 345);
    const auto *lsh_350 = buffer.data(lsh + 350);
    const auto *lsh_351 = buffer.data(lsh + 351);
    const auto *lsh_352 = buffer.data(lsh + 352);
    const auto *lsh_353 = buffer.data(lsh + 353);
    const auto *lsh_354 = buffer.data(lsh + 354);
    const auto *lsh_355 = buffer.data(lsh + 355);
    const auto *lsh_356 = buffer.data(lsh + 356);
    const auto *lsh_357 = buffer.data(lsh + 357);

    const auto *lsi1_252 = buffer.data(lsi1 + 252);
    const auto *lsi1_255 = buffer.data(lsi1 + 255);
    const auto *lsi1_257 = buffer.data(lsi1 + 257);
    const auto *lsi1_258 = buffer.data(lsi1 + 258);
    const auto *lsi1_261 = buffer.data(lsi1 + 261);
    const auto *lsi1_262 = buffer.data(lsi1 + 262);
    const auto *lsi1_264 = buffer.data(lsi1 + 264);
    const auto *lsi1_266 = buffer.data(lsi1 + 266);
    const auto *lsi1_279 = buffer.data(lsi1 + 279);
    const auto *lsi1_280 = buffer.data(lsi1 + 280);
    const auto *lsi1_283 = buffer.data(lsi1 + 283);
    const auto *lsi1_286 = buffer.data(lsi1 + 286);
    const auto *lsi1_290 = buffer.data(lsi1 + 290);
    const auto *lsi1_292 = buffer.data(lsi1 + 292);
    const auto *lsi1_301 = buffer.data(lsi1 + 301);

    const auto *msg0_194 = buffer.data(msg0 + 194);
    const auto *msg0_205 = buffer.data(msg0 + 205);
    const auto *msg0_207 = buffer.data(msg0 + 207);
    const auto *msg0_208 = buffer.data(msg0 + 208);
    const auto *msg0_209 = buffer.data(msg0 + 209);
    const auto *msg0_210 = buffer.data(msg0 + 210);
    const auto *msg0_211 = buffer.data(msg0 + 211);
    const auto *msg0_212 = buffer.data(msg0 + 212);
    const auto *msg0_213 = buffer.data(msg0 + 213);
    const auto *msg0_214 = buffer.data(msg0 + 214);
    const auto *msg0_215 = buffer.data(msg0 + 215);
    const auto *msg0_219 = buffer.data(msg0 + 219);
    const auto *msg0_220 = buffer.data(msg0 + 220);
    const auto *msg0_221 = buffer.data(msg0 + 221);
    const auto *msg0_222 = buffer.data(msg0 + 222);
    const auto *msg0_223 = buffer.data(msg0 + 223);
    const auto *msg0_224 = buffer.data(msg0 + 224);
    const auto *msg0_225 = buffer.data(msg0 + 225);
    const auto *msg0_227 = buffer.data(msg0 + 227);
    const auto *msg0_228 = buffer.data(msg0 + 228);
    const auto *msg0_230 = buffer.data(msg0 + 230);
    const auto *msg0_231 = buffer.data(msg0 + 231);
    const auto *msg0_235 = buffer.data(msg0 + 235);
    const auto *msg0_236 = buffer.data(msg0 + 236);
    const auto *msg0_237 = buffer.data(msg0 + 237);
    const auto *msg0_239 = buffer.data(msg0 + 239);
    const auto *msg0_245 = buffer.data(msg0 + 245);
    const auto *msg0_249 = buffer.data(msg0 + 249);
    const auto *msg0_252 = buffer.data(msg0 + 252);
    const auto *msg0_253 = buffer.data(msg0 + 253);
    const auto *msg0_254 = buffer.data(msg0 + 254);
    const auto *msg0_255 = buffer.data(msg0 + 255);

    const auto *msg1_194 = buffer.data(msg1 + 194);
    const auto *msg1_205 = buffer.data(msg1 + 205);
    const auto *msg1_207 = buffer.data(msg1 + 207);
    const auto *msg1_208 = buffer.data(msg1 + 208);
    const auto *msg1_209 = buffer.data(msg1 + 209);
    const auto *msg1_210 = buffer.data(msg1 + 210);
    const auto *msg1_211 = buffer.data(msg1 + 211);
    const auto *msg1_212 = buffer.data(msg1 + 212);
    const auto *msg1_213 = buffer.data(msg1 + 213);
    const auto *msg1_214 = buffer.data(msg1 + 214);
    const auto *msg1_215 = buffer.data(msg1 + 215);
    const auto *msg1_219 = buffer.data(msg1 + 219);
    const auto *msg1_220 = buffer.data(msg1 + 220);
    const auto *msg1_221 = buffer.data(msg1 + 221);
    const auto *msg1_222 = buffer.data(msg1 + 222);
    const auto *msg1_223 = buffer.data(msg1 + 223);
    const auto *msg1_224 = buffer.data(msg1 + 224);
    const auto *msg1_225 = buffer.data(msg1 + 225);
    const auto *msg1_227 = buffer.data(msg1 + 227);
    const auto *msg1_228 = buffer.data(msg1 + 228);
    const auto *msg1_230 = buffer.data(msg1 + 230);
    const auto *msg1_231 = buffer.data(msg1 + 231);
    const auto *msg1_235 = buffer.data(msg1 + 235);
    const auto *msg1_236 = buffer.data(msg1 + 236);
    const auto *msg1_237 = buffer.data(msg1 + 237);
    const auto *msg1_239 = buffer.data(msg1 + 239);
    const auto *msg1_245 = buffer.data(msg1 + 245);
    const auto *msg1_249 = buffer.data(msg1 + 249);
    const auto *msg1_252 = buffer.data(msg1 + 252);
    const auto *msg1_253 = buffer.data(msg1 + 253);
    const auto *msg1_254 = buffer.data(msg1 + 254);
    const auto *msg1_255 = buffer.data(msg1 + 255);

    const auto *msh_272 = buffer.data(msh + 272);
    const auto *msh_273 = buffer.data(msh + 273);
    const auto *msh_275 = buffer.data(msh + 275);
    const auto *msh_276 = buffer.data(msh + 276);
    const auto *msh_278 = buffer.data(msh + 278);
    const auto *msh_279 = buffer.data(msh + 279);
    const auto *msh_282 = buffer.data(msh + 282);
    const auto *msh_288 = buffer.data(msh + 288);
    const auto *msh_289 = buffer.data(msh + 289);
    const auto *msh_290 = buffer.data(msh + 290);
    const auto *msh_291 = buffer.data(msh + 291);
    const auto *msh_292 = buffer.data(msh + 292);
    const auto *msh_293 = buffer.data(msh + 293);
    const auto *msh_294 = buffer.data(msh + 294);
    const auto *msh_295 = buffer.data(msh + 295);
    const auto *msh_296 = buffer.data(msh + 296);
    const auto *msh_297 = buffer.data(msh + 297);
    const auto *msh_298 = buffer.data(msh + 298);
    const auto *msh_299 = buffer.data(msh + 299);
    const auto *msh_300 = buffer.data(msh + 300);
    const auto *msh_301 = buffer.data(msh + 301);
    const auto *msh_302 = buffer.data(msh + 302);
    const auto *msh_303 = buffer.data(msh + 303);
    const auto *msh_308 = buffer.data(msh + 308);
    const auto *msh_309 = buffer.data(msh + 309);
    const auto *msh_310 = buffer.data(msh + 310);
    const auto *msh_311 = buffer.data(msh + 311);
    const auto *msh_312 = buffer.data(msh + 312);
    const auto *msh_313 = buffer.data(msh + 313);
    const auto *msh_314 = buffer.data(msh + 314);
    const auto *msh_315 = buffer.data(msh + 315);
    const auto *msh_316 = buffer.data(msh + 316);
    const auto *msh_317 = buffer.data(msh + 317);
    const auto *msh_318 = buffer.data(msh + 318);
    const auto *msh_320 = buffer.data(msh + 320);
    const auto *msh_321 = buffer.data(msh + 321);
    const auto *msh_322 = buffer.data(msh + 322);
    const auto *msh_324 = buffer.data(msh + 324);
    const auto *msh_325 = buffer.data(msh + 325);
    const auto *msh_330 = buffer.data(msh + 330);
    const auto *msh_331 = buffer.data(msh + 331);
    const auto *msh_332 = buffer.data(msh + 332);
    const auto *msh_333 = buffer.data(msh + 333);
    const auto *msh_334 = buffer.data(msh + 334);
    const auto *msh_335 = buffer.data(msh + 335);
    const auto *msh_336 = buffer.data(msh + 336);
    const auto *msh_338 = buffer.data(msh + 338);
    const auto *msh_339 = buffer.data(msh + 339);
    const auto *msh_341 = buffer.data(msh + 341);
    const auto *msh_342 = buffer.data(msh + 342);
    const auto *msh_345 = buffer.data(msh + 345);
    const auto *msh_350 = buffer.data(msh + 350);
    const auto *msh_351 = buffer.data(msh + 351);
    const auto *msh_352 = buffer.data(msh + 352);
    const auto *msh_353 = buffer.data(msh + 353);
    const auto *msh_354 = buffer.data(msh + 354);
    const auto *msh_355 = buffer.data(msh + 355);
    const auto *msh_356 = buffer.data(msh + 356);
    const auto *msh_357 = buffer.data(msh + 357);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_y, pc_y, pc_z, lsi0_252, lsh_167, \
                         lsh_188, lsh_189, lsi1_252, msg0_194, msg1_194, msh_272, \
                         msh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * lsh_188[k]
                   + f_3 * pc_y[k] * msh_272[k];

        t_363[k] = f_12 * lsh_167[k]
                   + f_1 * msg0_194[k]
                   - f_2 * msg1_194[k]
                   + f_3 * pc_z[k] * msh_272[k];

        t_364[k] = pa_y[k] * lsi0_252[k]
                   - f_10 * pc_y[k] * lsi1_252[k];

        t_365[k] = f_11 * lsh_189[k]
                   + f_3 * pc_y[k] * msh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pc_y, pc_z, lsi0_255, lsi0_257, \
                         lsh_168, lsh_190, lsh_191, lsi1_255, lsi1_257, msh_273, \
                         msh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * lsh_168[k]
                   + f_3 * pc_z[k] * msh_273[k];

        t_367[k] = pa_y[k] * lsi0_255[k]
                   + f_12 * lsh_190[k]
                   - f_10 * pc_y[k] * lsi1_255[k];

        t_368[k] = f_11 * lsh_191[k]
                   + f_3 * pc_y[k] * msh_275[k];

        t_369[k] = pa_y[k] * lsi0_257[k]
                   - f_10 * pc_y[k] * lsi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_y, pc_y, pc_z, lsi0_258, lsi0_261, \
                         lsh_171, lsh_192, lsh_194, lsi1_258, lsi1_261, msh_276, \
                         msh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_y[k] * lsi0_258[k]
                   + f_13 * lsh_192[k]
                   - f_10 * pc_y[k] * lsi1_258[k];

        t_371[k] = f_13 * lsh_171[k]
                   + f_3 * pc_z[k] * msh_276[k];

        t_372[k] = f_11 * lsh_194[k]
                   + f_3 * pc_y[k] * msh_278[k];

        t_373[k] = pa_y[k] * lsi0_261[k]
                   - f_10 * pc_y[k] * lsi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_y, pc_y, pc_z, lsi0_262, lsi0_264, lsh_174, \
                         lsh_195, lsh_197, lsi1_262, lsi1_264, \
                         msh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * lsi0_262[k]
                   + f_14 * lsh_195[k]
                   - f_10 * pc_y[k] * lsi1_262[k];

        t_375[k] = f_13 * lsh_174[k]
                   + f_3 * pc_z[k] * msh_279[k];

        t_376[k] = pa_y[k] * lsi0_264[k]
                   + f_12 * lsh_197[k]
                   - f_10 * pc_y[k] * lsi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pc_x, pc_y, lsi0_266, lsh_198, \
                         lsh_288, lsh_289, lsi1_266, msh_282, msh_288, \
                         msh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * lsh_198[k]
                   + f_3 * pc_y[k] * msh_282[k];

        t_378[k] = pa_y[k] * lsi0_266[k]
                   - f_10 * pc_y[k] * lsi1_266[k];

        t_379[k] = f_20 * lsh_288[k]
                   + f_3 * pc_x[k] * msh_288[k];

        t_380[k] = f_20 * lsh_289[k]
                   + f_3 * pc_x[k] * msh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, lsh_290, lsh_291, lsh_292, lsh_293, \
                         msh_290, msh_291, msh_292, msh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_20 * lsh_290[k]
                   + f_3 * pc_x[k] * msh_290[k];

        t_382[k] = f_20 * lsh_291[k]
                   + f_3 * pc_x[k] * msh_291[k];

        t_383[k] = f_20 * lsh_292[k]
                   + f_3 * pc_x[k] * msh_292[k];

        t_384[k] = f_20 * lsh_293[k]
                   + f_3 * pc_x[k] * msh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, lsh_183, lsh_204, lsh_206, msg0_205, \
                         msg0_207, msg1_205, msg1_207, msh_288, \
                         msh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * lsh_204[k]
                   + f_1 * msg0_205[k]
                   - f_2 * msg1_205[k]
                   + f_3 * pc_y[k] * msh_288[k];

        t_386[k] = f_13 * lsh_183[k]
                   + f_3 * pc_z[k] * msh_288[k];

        t_387[k] = f_11 * lsh_206[k]
                   + f_8 * msg0_207[k]
                   - f_9 * msg1_207[k]
                   + f_3 * pc_y[k] * msh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, lsh_207, lsh_208, lsh_209, msg0_208, \
                         msg0_209, msg1_208, msg1_209, msh_291, msh_292, \
                         msh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * lsh_207[k]
                   + f_6 * msg0_208[k]
                   - f_7 * msg1_208[k]
                   + f_3 * pc_y[k] * msh_291[k];

        t_389[k] = f_11 * lsh_208[k]
                   + f_4 * msg0_209[k]
                   - f_5 * msg1_209[k]
                   + f_3 * pc_y[k] * msh_292[k];

        t_390[k] = f_11 * lsh_209[k]
                   + f_3 * pc_y[k] * msh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_x, pc_y, pc_z, lsi0_279, \
                         lsh_189, lsh_294, lsi1_279, msg0_210, msg1_210, \
                         msh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_y[k] * lsi0_279[k]
                   - f_10 * pc_y[k] * lsi1_279[k];

        t_392[k] = f_20 * lsh_294[k]
                   + f_1 * msg0_210[k]
                   - f_2 * msg1_210[k]
                   + f_3 * pc_x[k] * msh_294[k];

        t_393[k] = f_3 * pc_y[k] * msh_294[k];

        t_394[k] = f_14 * lsh_189[k]
                   + f_3 * pc_z[k] * msh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, lsh_299, msg0_210, msg0_215, \
                         msg1_210, msg1_215, msh_295, msh_296, \
                         msh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_4 * msg0_210[k]
                   - f_5 * msg1_210[k]
                   + f_3 * pc_y[k] * msh_295[k];

        t_396[k] = f_3 * pc_y[k] * msh_296[k];

        t_397[k] = f_20 * lsh_299[k]
                   + f_8 * msg0_215[k]
                   - f_9 * msg1_215[k]
                   + f_3 * pc_x[k] * msh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, msg0_211, msg0_212, msg1_211, msg1_212, \
                         msh_297, msh_298, msh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_6 * msg0_211[k]
                   - f_7 * msg1_211[k]
                   + f_3 * pc_y[k] * msh_297[k];

        t_399[k] = f_4 * msg0_212[k]
                   - f_5 * msg1_212[k]
                   + f_3 * pc_y[k] * msh_298[k];

        t_400[k] = f_3 * pc_y[k] * msh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_y, lsh_303, msg0_213, msg0_214, \
                         msg0_219, msg1_213, msg1_214, msg1_219, msh_300, msh_301, \
                         msh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_20 * lsh_303[k]
                   + f_6 * msg0_219[k]
                   - f_7 * msg1_219[k]
                   + f_3 * pc_x[k] * msh_303[k];

        t_402[k] = f_8 * msg0_213[k]
                   - f_9 * msg1_213[k]
                   + f_3 * pc_y[k] * msh_300[k];

        t_403[k] = f_6 * msg0_214[k]
                   - f_7 * msg1_214[k]
                   + f_3 * pc_y[k] * msh_301[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pc_x, pc_y, lsh_308, lsh_309, msg0_215, \
                         msg0_224, msg1_215, msg1_224, msh_302, msh_303, msh_308, \
                         msh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * msg0_215[k]
                   - f_5 * msg1_215[k]
                   + f_3 * pc_y[k] * msh_302[k];

        t_405[k] = f_3 * pc_y[k] * msh_303[k];

        t_406[k] = f_20 * lsh_308[k]
                   + f_4 * msg0_224[k]
                   - f_5 * msg1_224[k]
                   + f_3 * pc_x[k] * msh_308[k];

        t_407[k] = f_20 * lsh_309[k]
                   + f_3 * pc_x[k] * msh_309[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pc_x, pc_y, lsh_310, lsh_311, \
                         lsh_312, lsh_314, msh_308, msh_310, msh_311, msh_312, \
                         msh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_20 * lsh_310[k]
                   + f_3 * pc_x[k] * msh_310[k];

        t_409[k] = f_20 * lsh_311[k]
                   + f_3 * pc_x[k] * msh_311[k];

        t_410[k] = f_20 * lsh_312[k]
                   + f_3 * pc_x[k] * msh_312[k];

        t_411[k] = f_3 * pc_y[k] * msh_308[k];

        t_412[k] = f_20 * lsh_314[k]
                   + f_3 * pc_x[k] * msh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pc_y, msg0_220, msg0_221, msg0_222, msg1_220, \
                         msg1_221, msg1_222, msh_309, msh_310, \
                         msh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_1 * msg0_220[k]
                   - f_2 * msg1_220[k]
                   + f_3 * pc_y[k] * msh_309[k];

        t_414[k] = f_16 * msg0_221[k]
                   - f_17 * msg1_221[k]
                   + f_3 * pc_y[k] * msh_310[k];

        t_415[k] = f_8 * msg0_222[k]
                   - f_9 * msg1_222[k]
                   + f_3 * pc_y[k] * msh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, lsh_209, msg0_223, msg0_224, \
                         msg1_223, msg1_224, msh_312, msh_313, \
                         msh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * msg0_223[k]
                   - f_7 * msg1_223[k]
                   + f_3 * pc_y[k] * msh_312[k];

        t_417[k] = f_4 * msg0_224[k]
                   - f_5 * msg1_224[k]
                   + f_3 * pc_y[k] * msh_313[k];

        t_418[k] = f_3 * pc_y[k] * msh_314[k];

        t_419[k] = f_14 * lsh_209[k]
                   + f_1 * msg0_224[k]
                   - f_2 * msg1_224[k]
                   + f_3 * pc_z[k] * msh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, lsh_210, lsh_315, \
                         lsh_318, msg0_225, msg0_228, msg1_225, msg1_228, msh_315, \
                         msh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_14 * lsh_315[k]
                   + f_1 * msg0_225[k]
                   - f_2 * msg1_225[k]
                   + f_3 * pc_x[k] * msh_315[k];

        t_421[k] = f_20 * lsh_210[k]
                   + f_3 * pc_y[k] * msh_315[k];

        t_422[k] = f_3 * pc_z[k] * msh_315[k];

        t_423[k] = f_14 * lsh_318[k]
                   + f_8 * msg0_228[k]
                   - f_9 * msg1_228[k]
                   + f_3 * pc_x[k] * msh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pc_x, pc_z, lsh_321, msg0_225, msg0_231, \
                         msg1_225, msg1_231, msh_316, msh_317, msh_318, \
                         msh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_z[k] * msh_316[k];

        t_425[k] = f_4 * msg0_225[k]
                   - f_5 * msg1_225[k]
                   + f_3 * pc_z[k] * msh_317[k];

        t_426[k] = f_14 * lsh_321[k]
                   + f_6 * msg0_231[k]
                   - f_7 * msg1_231[k]
                   + f_3 * pc_x[k] * msh_321[k];

        t_427[k] = f_3 * pc_z[k] * msh_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pc_x, pc_y, pc_z, lsh_215, lsh_325, \
                         msg0_227, msg0_235, msg1_227, msg1_235, msh_320, msh_321, \
                         msh_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_20 * lsh_215[k]
                   + f_3 * pc_y[k] * msh_320[k];

        t_429[k] = f_6 * msg0_227[k]
                   - f_7 * msg1_227[k]
                   + f_3 * pc_z[k] * msh_320[k];

        t_430[k] = f_14 * lsh_325[k]
                   + f_4 * msg0_235[k]
                   - f_5 * msg1_235[k]
                   + f_3 * pc_x[k] * msh_325[k];

        t_431[k] = f_3 * pc_z[k] * msh_321[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, pc_y, pc_z, lsh_219, lsh_330, \
                         msg0_228, msg0_230, msg1_228, msg1_230, msh_322, msh_324, \
                         msh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_4 * msg0_228[k]
                   - f_5 * msg1_228[k]
                   + f_3 * pc_z[k] * msh_322[k];

        t_433[k] = f_20 * lsh_219[k]
                   + f_3 * pc_y[k] * msh_324[k];

        t_434[k] = f_8 * msg0_230[k]
                   - f_9 * msg1_230[k]
                   + f_3 * pc_z[k] * msh_324[k];

        t_435[k] = f_14 * lsh_330[k]
                   + f_3 * pc_x[k] * msh_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_z, lsh_332, lsh_333, \
                         lsh_334, lsh_335, msh_325, msh_332, msh_333, msh_334, \
                         msh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_3 * pc_z[k] * msh_325[k];

        t_437[k] = f_14 * lsh_332[k]
                   + f_3 * pc_x[k] * msh_332[k];

        t_438[k] = f_14 * lsh_333[k]
                   + f_3 * pc_x[k] * msh_333[k];

        t_439[k] = f_14 * lsh_334[k]
                   + f_3 * pc_x[k] * msh_334[k];

        t_440[k] = f_14 * lsh_335[k]
                   + f_3 * pc_x[k] * msh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, lsh_225, msg0_235, msg0_236, \
                         msg1_235, msg1_236, msh_330, msh_331, \
                         msh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_20 * lsh_225[k]
                   + f_1 * msg0_235[k]
                   - f_2 * msg1_235[k]
                   + f_3 * pc_y[k] * msh_330[k];

        t_442[k] = f_3 * pc_z[k] * msh_330[k];

        t_443[k] = f_4 * msg0_235[k]
                   - f_5 * msg1_235[k]
                   + f_3 * pc_z[k] * msh_331[k];

        t_444[k] = f_6 * msg0_236[k]
                   - f_7 * msg1_236[k]
                   + f_3 * pc_z[k] * msh_332[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pc_y, pc_z, lsi0_280, lsh_230, \
                         lsi1_280, msg0_237, msg0_239, msg1_237, msg1_239, msh_333, \
                         msh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * msg0_237[k]
                   - f_9 * msg1_237[k]
                   + f_3 * pc_z[k] * msh_333[k];

        t_446[k] = f_20 * lsh_230[k]
                   + f_3 * pc_y[k] * msh_335[k];

        t_447[k] = f_1 * msg0_239[k]
                   - f_2 * msg1_239[k]
                   + f_3 * pc_z[k] * msh_335[k];

        t_448[k] = pa_z[k] * lsi0_280[k]
                   - f_10 * pc_z[k] * lsi1_280[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_z, pc_y, pc_z, lsi0_283, lsh_210, \
                         lsh_231, lsh_233, lsi1_283, msh_336, msh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * lsh_231[k]
                   + f_3 * pc_y[k] * msh_336[k];

        t_450[k] = f_11 * lsh_210[k]
                   + f_3 * pc_z[k] * msh_336[k];

        t_451[k] = pa_z[k] * lsi0_283[k]
                   - f_10 * pc_z[k] * lsi1_283[k];

        t_452[k] = f_14 * lsh_233[k]
                   + f_3 * pc_y[k] * msh_338[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pa_z, pc_x, pc_z, lsi0_286, lsh_213, lsh_341, \
                         lsi1_286, msg0_245, msg1_245, msh_339, \
                         msh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_14 * lsh_341[k]
                   + f_8 * msg0_245[k]
                   - f_9 * msg1_245[k]
                   + f_3 * pc_x[k] * msh_341[k];

        t_454[k] = pa_z[k] * lsi0_286[k]
                   - f_10 * pc_z[k] * lsi1_286[k];

        t_455[k] = f_11 * lsh_213[k]
                   + f_3 * pc_z[k] * msh_339[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pa_z, pc_x, pc_y, pc_z, lsi0_290, lsh_236, \
                         lsh_345, lsi1_290, msg0_249, msg1_249, msh_341, \
                         msh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * lsh_236[k]
                   + f_3 * pc_y[k] * msh_341[k];

        t_457[k] = f_14 * lsh_345[k]
                   + f_6 * msg0_249[k]
                   - f_7 * msg1_249[k]
                   + f_3 * pc_x[k] * msh_345[k];

        t_458[k] = pa_z[k] * lsi0_290[k]
                   - f_10 * pc_z[k] * lsi1_290[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, pc_y, pc_z, lsi0_292, lsh_216, lsh_217, \
                         lsh_240, lsi1_292, msh_342, msh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_11 * lsh_216[k]
                   + f_3 * pc_z[k] * msh_342[k];

        t_460[k] = pa_z[k] * lsi0_292[k]
                   + f_12 * lsh_217[k]
                   - f_10 * pc_z[k] * lsi1_292[k];

        t_461[k] = f_14 * lsh_240[k]
                   + f_3 * pc_y[k] * msh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, lsh_350, lsh_351, lsh_352, lsh_353, \
                         msg0_254, msg1_254, msh_350, msh_351, msh_352, \
                         msh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_14 * lsh_350[k]
                   + f_4 * msg0_254[k]
                   - f_5 * msg1_254[k]
                   + f_3 * pc_x[k] * msh_350[k];

        t_463[k] = f_14 * lsh_351[k]
                   + f_3 * pc_x[k] * msh_351[k];

        t_464[k] = f_14 * lsh_352[k]
                   + f_3 * pc_x[k] * msh_352[k];

        t_465[k] = f_14 * lsh_353[k]
                   + f_3 * pc_x[k] * msh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pc_x, pc_z, lsi0_301, lsh_354, \
                         lsh_355, lsh_356, lsi1_301, msh_354, msh_355, \
                         msh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * lsh_354[k]
                   + f_3 * pc_x[k] * msh_354[k];

        t_467[k] = f_14 * lsh_355[k]
                   + f_3 * pc_x[k] * msh_355[k];

        t_468[k] = f_14 * lsh_356[k]
                   + f_3 * pc_x[k] * msh_356[k];

        t_469[k] = pa_z[k] * lsi0_301[k]
                   - f_10 * pc_z[k] * lsi1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, lsh_225, lsh_248, lsh_249, msg0_252, \
                         msg0_253, msg1_252, msg1_253, msh_351, msh_353, \
                         msh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * lsh_225[k]
                   + f_3 * pc_z[k] * msh_351[k];

        t_471[k] = f_14 * lsh_248[k]
                   + f_8 * msg0_252[k]
                   - f_9 * msg1_252[k]
                   + f_3 * pc_y[k] * msh_353[k];

        t_472[k] = f_14 * lsh_249[k]
                   + f_6 * msg0_253[k]
                   - f_7 * msg1_253[k]
                   + f_3 * pc_y[k] * msh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, lsh_230, lsh_250, lsh_251, msg0_254, \
                         msg1_254, msh_355, msh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * lsh_250[k]
                   + f_4 * msg0_254[k]
                   - f_5 * msg1_254[k]
                   + f_3 * pc_y[k] * msh_355[k];

        t_474[k] = f_14 * lsh_251[k]
                   + f_3 * pc_y[k] * msh_356[k];

        t_475[k] = f_11 * lsh_230[k]
                   + f_1 * msg0_254[k]
                   - f_2 * msg1_254[k]
                   + f_3 * pc_z[k] * msh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, lsh_231, lsh_252, lsh_357, \
                         msg0_255, msg1_255, msh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_14 * lsh_357[k]
                   + f_1 * msg0_255[k]
                   - f_2 * msg1_255[k]
                   + f_3 * pc_x[k] * msh_357[k];

        t_477[k] = f_13 * lsh_252[k]
                   + f_3 * pc_y[k] * msh_357[k];

        t_478[k] = f_12 * lsh_231[k]
                   + f_3 * pc_z[k] * msh_357[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *lsi0_392 = buffer.data(lsi0 + 392);
    const auto *lsi0_395 = buffer.data(lsi0 + 395);
    const auto *lsi0_397 = buffer.data(lsi0 + 397);
    const auto *lsi0_398 = buffer.data(lsi0 + 398);
    const auto *lsi0_401 = buffer.data(lsi0 + 401);
    const auto *lsi0_402 = buffer.data(lsi0 + 402);
    const auto *lsi0_404 = buffer.data(lsi0 + 404);
    const auto *lsi0_406 = buffer.data(lsi0 + 406);
    const auto *lsi0_419 = buffer.data(lsi0 + 419);

    const auto *lsh_234 = buffer.data(lsh + 234);
    const auto *lsh_237 = buffer.data(lsh + 237);
    const auto *lsh_246 = buffer.data(lsh + 246);
    const auto *lsh_251 = buffer.data(lsh + 251);
    const auto *lsh_252 = buffer.data(lsh + 252);
    const auto *lsh_254 = buffer.data(lsh + 254);
    const auto *lsh_255 = buffer.data(lsh + 255);
    const auto *lsh_257 = buffer.data(lsh + 257);
    const auto *lsh_258 = buffer.data(lsh + 258);
    const auto *lsh_261 = buffer.data(lsh + 261);
    const auto *lsh_267 = buffer.data(lsh + 267);
    const auto *lsh_269 = buffer.data(lsh + 269);
    const auto *lsh_270 = buffer.data(lsh + 270);
    const auto *lsh_271 = buffer.data(lsh + 271);
    const auto *lsh_272 = buffer.data(lsh + 272);
    const auto *lsh_273 = buffer.data(lsh + 273);
    const auto *lsh_275 = buffer.data(lsh + 275);
    const auto *lsh_276 = buffer.data(lsh + 276);
    const auto *lsh_278 = buffer.data(lsh + 278);
    const auto *lsh_279 = buffer.data(lsh + 279);
    const auto *lsh_282 = buffer.data(lsh + 282);
    const auto *lsh_288 = buffer.data(lsh + 288);
    const auto *lsh_290 = buffer.data(lsh + 290);
    const auto *lsh_291 = buffer.data(lsh + 291);
    const auto *lsh_292 = buffer.data(lsh + 292);
    const auto *lsh_293 = buffer.data(lsh + 293);
    const auto *lsh_294 = buffer.data(lsh + 294);
    const auto *lsh_295 = buffer.data(lsh + 295);
    const auto *lsh_296 = buffer.data(lsh + 296);
    const auto *lsh_297 = buffer.data(lsh + 297);
    const auto *lsh_299 = buffer.data(lsh + 299);
    const auto *lsh_300 = buffer.data(lsh + 300);
    const auto *lsh_302 = buffer.data(lsh + 302);
    const auto *lsh_303 = buffer.data(lsh + 303);
    const auto *lsh_309 = buffer.data(lsh + 309);
    const auto *lsh_311 = buffer.data(lsh + 311);
    const auto *lsh_312 = buffer.data(lsh + 312);
    const auto *lsh_313 = buffer.data(lsh + 313);
    const auto *lsh_314 = buffer.data(lsh + 314);
    const auto *lsh_315 = buffer.data(lsh + 315);
    const auto *lsh_360 = buffer.data(lsh + 360);
    const auto *lsh_362 = buffer.data(lsh + 362);
    const auto *lsh_363 = buffer.data(lsh + 363);
    const auto *lsh_366 = buffer.data(lsh + 366);
    const auto *lsh_367 = buffer.data(lsh + 367);
    const auto *lsh_369 = buffer.data(lsh + 369);
    const auto *lsh_371 = buffer.data(lsh + 371);
    const auto *lsh_372 = buffer.data(lsh + 372);
    const auto *lsh_373 = buffer.data(lsh + 373);
    const auto *lsh_374 = buffer.data(lsh + 374);
    const auto *lsh_375 = buffer.data(lsh + 375);
    const auto *lsh_376 = buffer.data(lsh + 376);
    const auto *lsh_377 = buffer.data(lsh + 377);
    const auto *lsh_378 = buffer.data(lsh + 378);
    const auto *lsh_381 = buffer.data(lsh + 381);
    const auto *lsh_383 = buffer.data(lsh + 383);
    const auto *lsh_384 = buffer.data(lsh + 384);
    const auto *lsh_387 = buffer.data(lsh + 387);
    const auto *lsh_388 = buffer.data(lsh + 388);
    const auto *lsh_390 = buffer.data(lsh + 390);
    const auto *lsh_392 = buffer.data(lsh + 392);
    const auto *lsh_393 = buffer.data(lsh + 393);
    const auto *lsh_394 = buffer.data(lsh + 394);
    const auto *lsh_395 = buffer.data(lsh + 395);
    const auto *lsh_396 = buffer.data(lsh + 396);
    const auto *lsh_397 = buffer.data(lsh + 397);
    const auto *lsh_398 = buffer.data(lsh + 398);
    const auto *lsh_414 = buffer.data(lsh + 414);
    const auto *lsh_415 = buffer.data(lsh + 415);
    const auto *lsh_416 = buffer.data(lsh + 416);
    const auto *lsh_417 = buffer.data(lsh + 417);
    const auto *lsh_418 = buffer.data(lsh + 418);
    const auto *lsh_419 = buffer.data(lsh + 419);
    const auto *lsh_420 = buffer.data(lsh + 420);
    const auto *lsh_425 = buffer.data(lsh + 425);
    const auto *lsh_429 = buffer.data(lsh + 429);
    const auto *lsh_434 = buffer.data(lsh + 434);
    const auto *lsh_435 = buffer.data(lsh + 435);
    const auto *lsh_436 = buffer.data(lsh + 436);
    const auto *lsh_437 = buffer.data(lsh + 437);
    const auto *lsh_438 = buffer.data(lsh + 438);
    const auto *lsh_440 = buffer.data(lsh + 440);
    const auto *lsh_441 = buffer.data(lsh + 441);
    const auto *lsh_444 = buffer.data(lsh + 444);

    const auto *lsi1_392 = buffer.data(lsi1 + 392);
    const auto *lsi1_395 = buffer.data(lsi1 + 395);
    const auto *lsi1_397 = buffer.data(lsi1 + 397);
    const auto *lsi1_398 = buffer.data(lsi1 + 398);
    const auto *lsi1_401 = buffer.data(lsi1 + 401);
    const auto *lsi1_402 = buffer.data(lsi1 + 402);
    const auto *lsi1_404 = buffer.data(lsi1 + 404);
    const auto *lsi1_406 = buffer.data(lsi1 + 406);
    const auto *lsi1_419 = buffer.data(lsi1 + 419);

    const auto *msg0_258 = buffer.data(msg0 + 258);
    const auto *msg0_260 = buffer.data(msg0 + 260);
    const auto *msg0_261 = buffer.data(msg0 + 261);
    const auto *msg0_264 = buffer.data(msg0 + 264);
    const auto *msg0_265 = buffer.data(msg0 + 265);
    const auto *msg0_267 = buffer.data(msg0 + 267);
    const auto *msg0_268 = buffer.data(msg0 + 268);
    const auto *msg0_269 = buffer.data(msg0 + 269);
    const auto *msg0_270 = buffer.data(msg0 + 270);
    const auto *msg0_273 = buffer.data(msg0 + 273);
    const auto *msg0_275 = buffer.data(msg0 + 275);
    const auto *msg0_276 = buffer.data(msg0 + 276);
    const auto *msg0_279 = buffer.data(msg0 + 279);
    const auto *msg0_280 = buffer.data(msg0 + 280);
    const auto *msg0_282 = buffer.data(msg0 + 282);
    const auto *msg0_283 = buffer.data(msg0 + 283);
    const auto *msg0_284 = buffer.data(msg0 + 284);
    const auto *msg0_295 = buffer.data(msg0 + 295);
    const auto *msg0_297 = buffer.data(msg0 + 297);
    const auto *msg0_298 = buffer.data(msg0 + 298);
    const auto *msg0_299 = buffer.data(msg0 + 299);
    const auto *msg0_300 = buffer.data(msg0 + 300);
    const auto *msg0_301 = buffer.data(msg0 + 301);
    const auto *msg0_302 = buffer.data(msg0 + 302);
    const auto *msg0_303 = buffer.data(msg0 + 303);
    const auto *msg0_304 = buffer.data(msg0 + 304);
    const auto *msg0_305 = buffer.data(msg0 + 305);
    const auto *msg0_309 = buffer.data(msg0 + 309);
    const auto *msg0_310 = buffer.data(msg0 + 310);
    const auto *msg0_311 = buffer.data(msg0 + 311);
    const auto *msg0_312 = buffer.data(msg0 + 312);
    const auto *msg0_313 = buffer.data(msg0 + 313);
    const auto *msg0_314 = buffer.data(msg0 + 314);
    const auto *msg0_315 = buffer.data(msg0 + 315);
    const auto *msg0_318 = buffer.data(msg0 + 318);

    const auto *msg1_258 = buffer.data(msg1 + 258);
    const auto *msg1_260 = buffer.data(msg1 + 260);
    const auto *msg1_261 = buffer.data(msg1 + 261);
    const auto *msg1_264 = buffer.data(msg1 + 264);
    const auto *msg1_265 = buffer.data(msg1 + 265);
    const auto *msg1_267 = buffer.data(msg1 + 267);
    const auto *msg1_268 = buffer.data(msg1 + 268);
    const auto *msg1_269 = buffer.data(msg1 + 269);
    const auto *msg1_270 = buffer.data(msg1 + 270);
    const auto *msg1_273 = buffer.data(msg1 + 273);
    const auto *msg1_275 = buffer.data(msg1 + 275);
    const auto *msg1_276 = buffer.data(msg1 + 276);
    const auto *msg1_279 = buffer.data(msg1 + 279);
    const auto *msg1_280 = buffer.data(msg1 + 280);
    const auto *msg1_282 = buffer.data(msg1 + 282);
    const auto *msg1_283 = buffer.data(msg1 + 283);
    const auto *msg1_284 = buffer.data(msg1 + 284);
    const auto *msg1_295 = buffer.data(msg1 + 295);
    const auto *msg1_297 = buffer.data(msg1 + 297);
    const auto *msg1_298 = buffer.data(msg1 + 298);
    const auto *msg1_299 = buffer.data(msg1 + 299);
    const auto *msg1_300 = buffer.data(msg1 + 300);
    const auto *msg1_301 = buffer.data(msg1 + 301);
    const auto *msg1_302 = buffer.data(msg1 + 302);
    const auto *msg1_303 = buffer.data(msg1 + 303);
    const auto *msg1_304 = buffer.data(msg1 + 304);
    const auto *msg1_305 = buffer.data(msg1 + 305);
    const auto *msg1_309 = buffer.data(msg1 + 309);
    const auto *msg1_310 = buffer.data(msg1 + 310);
    const auto *msg1_311 = buffer.data(msg1 + 311);
    const auto *msg1_312 = buffer.data(msg1 + 312);
    const auto *msg1_313 = buffer.data(msg1 + 313);
    const auto *msg1_314 = buffer.data(msg1 + 314);
    const auto *msg1_315 = buffer.data(msg1 + 315);
    const auto *msg1_318 = buffer.data(msg1 + 318);

    const auto *msh_359 = buffer.data(msh + 359);
    const auto *msh_360 = buffer.data(msh + 360);
    const auto *msh_362 = buffer.data(msh + 362);
    const auto *msh_363 = buffer.data(msh + 363);
    const auto *msh_366 = buffer.data(msh + 366);
    const auto *msh_367 = buffer.data(msh + 367);
    const auto *msh_369 = buffer.data(msh + 369);
    const auto *msh_371 = buffer.data(msh + 371);
    const auto *msh_372 = buffer.data(msh + 372);
    const auto *msh_373 = buffer.data(msh + 373);
    const auto *msh_374 = buffer.data(msh + 374);
    const auto *msh_375 = buffer.data(msh + 375);
    const auto *msh_376 = buffer.data(msh + 376);
    const auto *msh_377 = buffer.data(msh + 377);
    const auto *msh_378 = buffer.data(msh + 378);
    const auto *msh_380 = buffer.data(msh + 380);
    const auto *msh_381 = buffer.data(msh + 381);
    const auto *msh_383 = buffer.data(msh + 383);
    const auto *msh_384 = buffer.data(msh + 384);
    const auto *msh_387 = buffer.data(msh + 387);
    const auto *msh_388 = buffer.data(msh + 388);
    const auto *msh_390 = buffer.data(msh + 390);
    const auto *msh_392 = buffer.data(msh + 392);
    const auto *msh_393 = buffer.data(msh + 393);
    const auto *msh_394 = buffer.data(msh + 394);
    const auto *msh_395 = buffer.data(msh + 395);
    const auto *msh_396 = buffer.data(msh + 396);
    const auto *msh_397 = buffer.data(msh + 397);
    const auto *msh_398 = buffer.data(msh + 398);
    const auto *msh_399 = buffer.data(msh + 399);
    const auto *msh_401 = buffer.data(msh + 401);
    const auto *msh_402 = buffer.data(msh + 402);
    const auto *msh_404 = buffer.data(msh + 404);
    const auto *msh_405 = buffer.data(msh + 405);
    const auto *msh_408 = buffer.data(msh + 408);
    const auto *msh_414 = buffer.data(msh + 414);
    const auto *msh_415 = buffer.data(msh + 415);
    const auto *msh_416 = buffer.data(msh + 416);
    const auto *msh_417 = buffer.data(msh + 417);
    const auto *msh_418 = buffer.data(msh + 418);
    const auto *msh_419 = buffer.data(msh + 419);
    const auto *msh_420 = buffer.data(msh + 420);
    const auto *msh_421 = buffer.data(msh + 421);
    const auto *msh_422 = buffer.data(msh + 422);
    const auto *msh_423 = buffer.data(msh + 423);
    const auto *msh_424 = buffer.data(msh + 424);
    const auto *msh_425 = buffer.data(msh + 425);
    const auto *msh_426 = buffer.data(msh + 426);
    const auto *msh_427 = buffer.data(msh + 427);
    const auto *msh_428 = buffer.data(msh + 428);
    const auto *msh_429 = buffer.data(msh + 429);
    const auto *msh_434 = buffer.data(msh + 434);
    const auto *msh_435 = buffer.data(msh + 435);
    const auto *msh_436 = buffer.data(msh + 436);
    const auto *msh_437 = buffer.data(msh + 437);
    const auto *msh_438 = buffer.data(msh + 438);
    const auto *msh_439 = buffer.data(msh + 439);
    const auto *msh_440 = buffer.data(msh + 440);
    const auto *msh_441 = buffer.data(msh + 441);
    const auto *msh_444 = buffer.data(msh + 444);

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, lsh_254, lsh_360, lsh_362, msg0_258, \
                         msg0_260, msg1_258, msg1_260, msh_359, msh_360, \
                         msh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_14 * lsh_360[k]
                   + f_8 * msg0_258[k]
                   - f_9 * msg1_258[k]
                   + f_3 * pc_x[k] * msh_360[k];

        t_480[k] = f_13 * lsh_254[k]
                   + f_3 * pc_y[k] * msh_359[k];

        t_481[k] = f_14 * lsh_362[k]
                   + f_8 * msg0_260[k]
                   - f_9 * msg1_260[k]
                   + f_3 * pc_x[k] * msh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, lsh_234, lsh_257, lsh_363, \
                         msg0_261, msg1_261, msh_360, msh_362, \
                         msh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_14 * lsh_363[k]
                   + f_6 * msg0_261[k]
                   - f_7 * msg1_261[k]
                   + f_3 * pc_x[k] * msh_363[k];

        t_483[k] = f_12 * lsh_234[k]
                   + f_3 * pc_z[k] * msh_360[k];

        t_484[k] = f_13 * lsh_257[k]
                   + f_3 * pc_y[k] * msh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, lsh_237, lsh_366, lsh_367, msg0_264, \
                         msg0_265, msg1_264, msg1_265, msh_363, msh_366, \
                         msh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_14 * lsh_366[k]
                   + f_6 * msg0_264[k]
                   - f_7 * msg1_264[k]
                   + f_3 * pc_x[k] * msh_366[k];

        t_486[k] = f_14 * lsh_367[k]
                   + f_4 * msg0_265[k]
                   - f_5 * msg1_265[k]
                   + f_3 * pc_x[k] * msh_367[k];

        t_487[k] = f_12 * lsh_237[k]
                   + f_3 * pc_z[k] * msh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, lsh_261, lsh_369, lsh_371, msg0_267, \
                         msg0_269, msg1_267, msg1_269, msh_366, msh_369, \
                         msh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_14 * lsh_369[k]
                   + f_4 * msg0_267[k]
                   - f_5 * msg1_267[k]
                   + f_3 * pc_x[k] * msh_369[k];

        t_489[k] = f_13 * lsh_261[k]
                   + f_3 * pc_y[k] * msh_366[k];

        t_490[k] = f_14 * lsh_371[k]
                   + f_4 * msg0_269[k]
                   - f_5 * msg1_269[k]
                   + f_3 * pc_x[k] * msh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, lsh_372, lsh_373, lsh_374, \
                         lsh_375, lsh_376, msh_372, msh_373, msh_374, msh_375, \
                         msh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_14 * lsh_372[k]
                   + f_3 * pc_x[k] * msh_372[k];

        t_492[k] = f_14 * lsh_373[k]
                   + f_3 * pc_x[k] * msh_373[k];

        t_493[k] = f_14 * lsh_374[k]
                   + f_3 * pc_x[k] * msh_374[k];

        t_494[k] = f_14 * lsh_375[k]
                   + f_3 * pc_x[k] * msh_375[k];

        t_495[k] = f_14 * lsh_376[k]
                   + f_3 * pc_x[k] * msh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, lsh_246, lsh_267, lsh_377, \
                         msg0_265, msg1_265, msh_372, msh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_14 * lsh_377[k]
                   + f_3 * pc_x[k] * msh_377[k];

        t_497[k] = f_13 * lsh_267[k]
                   + f_1 * msg0_265[k]
                   - f_2 * msg1_265[k]
                   + f_3 * pc_y[k] * msh_372[k];

        t_498[k] = f_12 * lsh_246[k]
                   + f_3 * pc_z[k] * msh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, lsh_269, lsh_270, lsh_271, msg0_267, \
                         msg0_268, msg0_269, msg1_267, msg1_268, msg1_269, msh_374, msh_375, \
                         msh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * lsh_269[k]
                   + f_8 * msg0_267[k]
                   - f_9 * msg1_267[k]
                   + f_3 * pc_y[k] * msh_374[k];

        t_500[k] = f_13 * lsh_270[k]
                   + f_6 * msg0_268[k]
                   - f_7 * msg1_268[k]
                   + f_3 * pc_y[k] * msh_375[k];

        t_501[k] = f_13 * lsh_271[k]
                   + f_4 * msg0_269[k]
                   - f_5 * msg1_269[k]
                   + f_3 * pc_y[k] * msh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, lsh_251, lsh_272, lsh_378, \
                         msg0_269, msg0_270, msg1_269, msg1_270, msh_377, \
                         msh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * lsh_272[k]
                   + f_3 * pc_y[k] * msh_377[k];

        t_503[k] = f_12 * lsh_251[k]
                   + f_1 * msg0_269[k]
                   - f_2 * msg1_269[k]
                   + f_3 * pc_z[k] * msh_377[k];

        t_504[k] = f_14 * lsh_378[k]
                   + f_1 * msg0_270[k]
                   - f_2 * msg1_270[k]
                   + f_3 * pc_x[k] * msh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, lsh_252, lsh_273, \
                         lsh_275, lsh_381, msg0_273, msg1_273, msh_378, msh_380, \
                         msh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * lsh_273[k]
                   + f_3 * pc_y[k] * msh_378[k];

        t_506[k] = f_13 * lsh_252[k]
                   + f_3 * pc_z[k] * msh_378[k];

        t_507[k] = f_14 * lsh_381[k]
                   + f_8 * msg0_273[k]
                   - f_9 * msg1_273[k]
                   + f_3 * pc_x[k] * msh_381[k];

        t_508[k] = f_12 * lsh_275[k]
                   + f_3 * pc_y[k] * msh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, lsh_255, lsh_383, lsh_384, msg0_275, \
                         msg0_276, msg1_275, msg1_276, msh_381, msh_383, \
                         msh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_14 * lsh_383[k]
                   + f_8 * msg0_275[k]
                   - f_9 * msg1_275[k]
                   + f_3 * pc_x[k] * msh_383[k];

        t_510[k] = f_14 * lsh_384[k]
                   + f_6 * msg0_276[k]
                   - f_7 * msg1_276[k]
                   + f_3 * pc_x[k] * msh_384[k];

        t_511[k] = f_13 * lsh_255[k]
                   + f_3 * pc_z[k] * msh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, lsh_278, lsh_387, lsh_388, msg0_279, \
                         msg0_280, msg1_279, msg1_280, msh_383, msh_387, \
                         msh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * lsh_278[k]
                   + f_3 * pc_y[k] * msh_383[k];

        t_513[k] = f_14 * lsh_387[k]
                   + f_6 * msg0_279[k]
                   - f_7 * msg1_279[k]
                   + f_3 * pc_x[k] * msh_387[k];

        t_514[k] = f_14 * lsh_388[k]
                   + f_4 * msg0_280[k]
                   - f_5 * msg1_280[k]
                   + f_3 * pc_x[k] * msh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, lsh_258, lsh_282, lsh_390, \
                         msg0_282, msg1_282, msh_384, msh_387, \
                         msh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * lsh_258[k]
                   + f_3 * pc_z[k] * msh_384[k];

        t_516[k] = f_14 * lsh_390[k]
                   + f_4 * msg0_282[k]
                   - f_5 * msg1_282[k]
                   + f_3 * pc_x[k] * msh_390[k];

        t_517[k] = f_12 * lsh_282[k]
                   + f_3 * pc_y[k] * msh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, lsh_392, lsh_393, lsh_394, lsh_395, \
                         msg0_284, msg1_284, msh_392, msh_393, msh_394, \
                         msh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_14 * lsh_392[k]
                   + f_4 * msg0_284[k]
                   - f_5 * msg1_284[k]
                   + f_3 * pc_x[k] * msh_392[k];

        t_519[k] = f_14 * lsh_393[k]
                   + f_3 * pc_x[k] * msh_393[k];

        t_520[k] = f_14 * lsh_394[k]
                   + f_3 * pc_x[k] * msh_394[k];

        t_521[k] = f_14 * lsh_395[k]
                   + f_3 * pc_x[k] * msh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, lsh_288, lsh_396, lsh_397, \
                         lsh_398, msg0_280, msg1_280, msh_393, msh_396, msh_397, \
                         msh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_14 * lsh_396[k]
                   + f_3 * pc_x[k] * msh_396[k];

        t_523[k] = f_14 * lsh_397[k]
                   + f_3 * pc_x[k] * msh_397[k];

        t_524[k] = f_14 * lsh_398[k]
                   + f_3 * pc_x[k] * msh_398[k];

        t_525[k] = f_12 * lsh_288[k]
                   + f_1 * msg0_280[k]
                   - f_2 * msg1_280[k]
                   + f_3 * pc_y[k] * msh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, lsh_267, lsh_290, lsh_291, msg0_282, \
                         msg0_283, msg1_282, msg1_283, msh_393, msh_395, \
                         msh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * lsh_267[k]
                   + f_3 * pc_z[k] * msh_393[k];

        t_527[k] = f_12 * lsh_290[k]
                   + f_8 * msg0_282[k]
                   - f_9 * msg1_282[k]
                   + f_3 * pc_y[k] * msh_395[k];

        t_528[k] = f_12 * lsh_291[k]
                   + f_6 * msg0_283[k]
                   - f_7 * msg1_283[k]
                   + f_3 * pc_y[k] * msh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pa_y, pc_y, pc_z, lsi0_392, lsh_272, \
                         lsh_292, lsh_293, lsi1_392, msg0_284, msg1_284, msh_397, \
                         msh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * lsh_292[k]
                   + f_4 * msg0_284[k]
                   - f_5 * msg1_284[k]
                   + f_3 * pc_y[k] * msh_397[k];

        t_530[k] = f_12 * lsh_293[k]
                   + f_3 * pc_y[k] * msh_398[k];

        t_531[k] = f_13 * lsh_272[k]
                   + f_1 * msg0_284[k]
                   - f_2 * msg1_284[k]
                   + f_3 * pc_z[k] * msh_398[k];

        t_532[k] = pa_y[k] * lsi0_392[k]
                   - f_10 * pc_y[k] * lsi1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_y, pc_y, pc_z, lsi0_395, lsh_273, \
                         lsh_294, lsh_295, lsh_296, lsi1_395, msh_399, \
                         msh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * lsh_294[k]
                   + f_3 * pc_y[k] * msh_399[k];

        t_534[k] = f_14 * lsh_273[k]
                   + f_3 * pc_z[k] * msh_399[k];

        t_535[k] = pa_y[k] * lsi0_395[k]
                   + f_12 * lsh_295[k]
                   - f_10 * pc_y[k] * lsi1_395[k];

        t_536[k] = f_11 * lsh_296[k]
                   + f_3 * pc_y[k] * msh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_y, pc_y, pc_z, lsi0_397, lsi0_398, \
                         lsh_276, lsh_297, lsh_299, lsi1_397, lsi1_398, msh_402, \
                         msh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pa_y[k] * lsi0_397[k]
                   - f_10 * pc_y[k] * lsi1_397[k];

        t_538[k] = pa_y[k] * lsi0_398[k]
                   + f_13 * lsh_297[k]
                   - f_10 * pc_y[k] * lsi1_398[k];

        t_539[k] = f_14 * lsh_276[k]
                   + f_3 * pc_z[k] * msh_402[k];

        t_540[k] = f_11 * lsh_299[k]
                   + f_3 * pc_y[k] * msh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pa_y, pc_y, pc_z, lsi0_401, lsi0_402, lsh_279, \
                         lsh_300, lsi1_401, lsi1_402, msh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pa_y[k] * lsi0_401[k]
                   - f_10 * pc_y[k] * lsi1_401[k];

        t_542[k] = pa_y[k] * lsi0_402[k]
                   + f_14 * lsh_300[k]
                   - f_10 * pc_y[k] * lsi1_402[k];

        t_543[k] = f_14 * lsh_279[k]
                   + f_3 * pc_z[k] * msh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_x, pc_y, lsi0_404, lsi0_406, \
                         lsh_302, lsh_303, lsh_414, lsi1_404, lsi1_406, msh_408, \
                         msh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pa_y[k] * lsi0_404[k]
                   + f_12 * lsh_302[k]
                   - f_10 * pc_y[k] * lsi1_404[k];

        t_545[k] = f_11 * lsh_303[k]
                   + f_3 * pc_y[k] * msh_408[k];

        t_546[k] = pa_y[k] * lsi0_406[k]
                   - f_10 * pc_y[k] * lsi1_406[k];

        t_547[k] = f_14 * lsh_414[k]
                   + f_3 * pc_x[k] * msh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, lsh_415, lsh_416, lsh_417, \
                         lsh_418, lsh_419, msh_415, msh_416, msh_417, msh_418, \
                         msh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_14 * lsh_415[k]
                   + f_3 * pc_x[k] * msh_415[k];

        t_549[k] = f_14 * lsh_416[k]
                   + f_3 * pc_x[k] * msh_416[k];

        t_550[k] = f_14 * lsh_417[k]
                   + f_3 * pc_x[k] * msh_417[k];

        t_551[k] = f_14 * lsh_418[k]
                   + f_3 * pc_x[k] * msh_418[k];

        t_552[k] = f_14 * lsh_419[k]
                   + f_3 * pc_x[k] * msh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, lsh_288, lsh_309, lsh_311, msg0_295, \
                         msg0_297, msg1_295, msg1_297, msh_414, \
                         msh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * lsh_309[k]
                   + f_1 * msg0_295[k]
                   - f_2 * msg1_295[k]
                   + f_3 * pc_y[k] * msh_414[k];

        t_554[k] = f_14 * lsh_288[k]
                   + f_3 * pc_z[k] * msh_414[k];

        t_555[k] = f_11 * lsh_311[k]
                   + f_8 * msg0_297[k]
                   - f_9 * msg1_297[k]
                   + f_3 * pc_y[k] * msh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, lsh_312, lsh_313, lsh_314, msg0_298, \
                         msg0_299, msg1_298, msg1_299, msh_417, msh_418, \
                         msh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * lsh_312[k]
                   + f_6 * msg0_298[k]
                   - f_7 * msg1_298[k]
                   + f_3 * pc_y[k] * msh_417[k];

        t_557[k] = f_11 * lsh_313[k]
                   + f_4 * msg0_299[k]
                   - f_5 * msg1_299[k]
                   + f_3 * pc_y[k] * msh_418[k];

        t_558[k] = f_11 * lsh_314[k]
                   + f_3 * pc_y[k] * msh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_y, pc_x, pc_y, pc_z, lsi0_419, \
                         lsh_294, lsh_420, lsi1_419, msg0_300, msg1_300, \
                         msh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pa_y[k] * lsi0_419[k]
                   - f_10 * pc_y[k] * lsi1_419[k];

        t_560[k] = f_14 * lsh_420[k]
                   + f_1 * msg0_300[k]
                   - f_2 * msg1_300[k]
                   + f_3 * pc_x[k] * msh_420[k];

        t_561[k] = f_3 * pc_y[k] * msh_420[k];

        t_562[k] = f_20 * lsh_294[k]
                   + f_3 * pc_z[k] * msh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, lsh_425, msg0_300, msg0_305, \
                         msg1_300, msg1_305, msh_421, msh_422, \
                         msh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_4 * msg0_300[k]
                   - f_5 * msg1_300[k]
                   + f_3 * pc_y[k] * msh_421[k];

        t_564[k] = f_3 * pc_y[k] * msh_422[k];

        t_565[k] = f_14 * lsh_425[k]
                   + f_8 * msg0_305[k]
                   - f_9 * msg1_305[k]
                   + f_3 * pc_x[k] * msh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_y, msg0_301, msg0_302, msg1_301, msg1_302, \
                         msh_423, msh_424, msh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_6 * msg0_301[k]
                   - f_7 * msg1_301[k]
                   + f_3 * pc_y[k] * msh_423[k];

        t_567[k] = f_4 * msg0_302[k]
                   - f_5 * msg1_302[k]
                   + f_3 * pc_y[k] * msh_424[k];

        t_568[k] = f_3 * pc_y[k] * msh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_y, lsh_429, msg0_303, msg0_304, \
                         msg0_309, msg1_303, msg1_304, msg1_309, msh_426, msh_427, \
                         msh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_14 * lsh_429[k]
                   + f_6 * msg0_309[k]
                   - f_7 * msg1_309[k]
                   + f_3 * pc_x[k] * msh_429[k];

        t_570[k] = f_8 * msg0_303[k]
                   - f_9 * msg1_303[k]
                   + f_3 * pc_y[k] * msh_426[k];

        t_571[k] = f_6 * msg0_304[k]
                   - f_7 * msg1_304[k]
                   + f_3 * pc_y[k] * msh_427[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, lsh_434, lsh_435, msg0_305, \
                         msg0_314, msg1_305, msg1_314, msh_428, msh_429, msh_434, \
                         msh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_4 * msg0_305[k]
                   - f_5 * msg1_305[k]
                   + f_3 * pc_y[k] * msh_428[k];

        t_573[k] = f_3 * pc_y[k] * msh_429[k];

        t_574[k] = f_14 * lsh_434[k]
                   + f_4 * msg0_314[k]
                   - f_5 * msg1_314[k]
                   + f_3 * pc_x[k] * msh_434[k];

        t_575[k] = f_14 * lsh_435[k]
                   + f_3 * pc_x[k] * msh_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, pc_x, pc_y, lsh_436, lsh_437, \
                         lsh_438, lsh_440, msh_434, msh_436, msh_437, msh_438, \
                         msh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_14 * lsh_436[k]
                   + f_3 * pc_x[k] * msh_436[k];

        t_577[k] = f_14 * lsh_437[k]
                   + f_3 * pc_x[k] * msh_437[k];

        t_578[k] = f_14 * lsh_438[k]
                   + f_3 * pc_x[k] * msh_438[k];

        t_579[k] = f_3 * pc_y[k] * msh_434[k];

        t_580[k] = f_14 * lsh_440[k]
                   + f_3 * pc_x[k] * msh_440[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, msg0_310, msg0_311, msg0_312, msg1_310, \
                         msg1_311, msg1_312, msh_435, msh_436, \
                         msh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_1 * msg0_310[k]
                   - f_2 * msg1_310[k]
                   + f_3 * pc_y[k] * msh_435[k];

        t_582[k] = f_16 * msg0_311[k]
                   - f_17 * msg1_311[k]
                   + f_3 * pc_y[k] * msh_436[k];

        t_583[k] = f_8 * msg0_312[k]
                   - f_9 * msg1_312[k]
                   + f_3 * pc_y[k] * msh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, lsh_314, msg0_313, msg0_314, \
                         msg1_313, msg1_314, msh_438, msh_439, \
                         msh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * msg0_313[k]
                   - f_7 * msg1_313[k]
                   + f_3 * pc_y[k] * msh_438[k];

        t_585[k] = f_4 * msg0_314[k]
                   - f_5 * msg1_314[k]
                   + f_3 * pc_y[k] * msh_439[k];

        t_586[k] = f_3 * pc_y[k] * msh_440[k];

        t_587[k] = f_20 * lsh_314[k]
                   + f_1 * msg0_314[k]
                   - f_2 * msg1_314[k]
                   + f_3 * pc_z[k] * msh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, lsh_315, lsh_441, \
                         lsh_444, msg0_315, msg0_318, msg1_315, msg1_318, msh_441, \
                         msh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_13 * lsh_441[k]
                   + f_1 * msg0_315[k]
                   - f_2 * msg1_315[k]
                   + f_3 * pc_x[k] * msh_441[k];

        t_589[k] = f_19 * lsh_315[k]
                   + f_3 * pc_y[k] * msh_441[k];

        t_590[k] = f_3 * pc_z[k] * msh_441[k];

        t_591[k] = f_13 * lsh_444[k]
                   + f_8 * msg0_318[k]
                   - f_9 * msg1_318[k]
                   + f_3 * pc_x[k] * msh_444[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *lsi0_420 = buffer.data(lsi0 + 420);
    const auto *lsi0_423 = buffer.data(lsi0 + 423);
    const auto *lsi0_426 = buffer.data(lsi0 + 426);
    const auto *lsi0_430 = buffer.data(lsi0 + 430);
    const auto *lsi0_432 = buffer.data(lsi0 + 432);
    const auto *lsi0_441 = buffer.data(lsi0 + 441);

    const auto *lsh_315 = buffer.data(lsh + 315);
    const auto *lsh_318 = buffer.data(lsh + 318);
    const auto *lsh_320 = buffer.data(lsh + 320);
    const auto *lsh_321 = buffer.data(lsh + 321);
    const auto *lsh_322 = buffer.data(lsh + 322);
    const auto *lsh_324 = buffer.data(lsh + 324);
    const auto *lsh_330 = buffer.data(lsh + 330);
    const auto *lsh_335 = buffer.data(lsh + 335);
    const auto *lsh_336 = buffer.data(lsh + 336);
    const auto *lsh_338 = buffer.data(lsh + 338);
    const auto *lsh_339 = buffer.data(lsh + 339);
    const auto *lsh_341 = buffer.data(lsh + 341);
    const auto *lsh_342 = buffer.data(lsh + 342);
    const auto *lsh_345 = buffer.data(lsh + 345);
    const auto *lsh_351 = buffer.data(lsh + 351);
    const auto *lsh_353 = buffer.data(lsh + 353);
    const auto *lsh_354 = buffer.data(lsh + 354);
    const auto *lsh_355 = buffer.data(lsh + 355);
    const auto *lsh_356 = buffer.data(lsh + 356);
    const auto *lsh_357 = buffer.data(lsh + 357);
    const auto *lsh_359 = buffer.data(lsh + 359);
    const auto *lsh_360 = buffer.data(lsh + 360);
    const auto *lsh_362 = buffer.data(lsh + 362);
    const auto *lsh_363 = buffer.data(lsh + 363);
    const auto *lsh_366 = buffer.data(lsh + 366);
    const auto *lsh_372 = buffer.data(lsh + 372);
    const auto *lsh_374 = buffer.data(lsh + 374);
    const auto *lsh_375 = buffer.data(lsh + 375);
    const auto *lsh_376 = buffer.data(lsh + 376);
    const auto *lsh_377 = buffer.data(lsh + 377);
    const auto *lsh_378 = buffer.data(lsh + 378);
    const auto *lsh_380 = buffer.data(lsh + 380);
    const auto *lsh_383 = buffer.data(lsh + 383);
    const auto *lsh_387 = buffer.data(lsh + 387);
    const auto *lsh_393 = buffer.data(lsh + 393);
    const auto *lsh_395 = buffer.data(lsh + 395);
    const auto *lsh_396 = buffer.data(lsh + 396);
    const auto *lsh_397 = buffer.data(lsh + 397);
    const auto *lsh_398 = buffer.data(lsh + 398);
    const auto *lsh_399 = buffer.data(lsh + 399);
    const auto *lsh_447 = buffer.data(lsh + 447);
    const auto *lsh_451 = buffer.data(lsh + 451);
    const auto *lsh_456 = buffer.data(lsh + 456);
    const auto *lsh_458 = buffer.data(lsh + 458);
    const auto *lsh_459 = buffer.data(lsh + 459);
    const auto *lsh_460 = buffer.data(lsh + 460);
    const auto *lsh_461 = buffer.data(lsh + 461);
    const auto *lsh_467 = buffer.data(lsh + 467);
    const auto *lsh_471 = buffer.data(lsh + 471);
    const auto *lsh_476 = buffer.data(lsh + 476);
    const auto *lsh_477 = buffer.data(lsh + 477);
    const auto *lsh_478 = buffer.data(lsh + 478);
    const auto *lsh_479 = buffer.data(lsh + 479);
    const auto *lsh_480 = buffer.data(lsh + 480);
    const auto *lsh_481 = buffer.data(lsh + 481);
    const auto *lsh_482 = buffer.data(lsh + 482);
    const auto *lsh_483 = buffer.data(lsh + 483);
    const auto *lsh_486 = buffer.data(lsh + 486);
    const auto *lsh_488 = buffer.data(lsh + 488);
    const auto *lsh_489 = buffer.data(lsh + 489);
    const auto *lsh_492 = buffer.data(lsh + 492);
    const auto *lsh_493 = buffer.data(lsh + 493);
    const auto *lsh_495 = buffer.data(lsh + 495);
    const auto *lsh_497 = buffer.data(lsh + 497);
    const auto *lsh_498 = buffer.data(lsh + 498);
    const auto *lsh_499 = buffer.data(lsh + 499);
    const auto *lsh_500 = buffer.data(lsh + 500);
    const auto *lsh_501 = buffer.data(lsh + 501);
    const auto *lsh_502 = buffer.data(lsh + 502);
    const auto *lsh_503 = buffer.data(lsh + 503);
    const auto *lsh_504 = buffer.data(lsh + 504);
    const auto *lsh_507 = buffer.data(lsh + 507);
    const auto *lsh_509 = buffer.data(lsh + 509);
    const auto *lsh_510 = buffer.data(lsh + 510);
    const auto *lsh_513 = buffer.data(lsh + 513);
    const auto *lsh_514 = buffer.data(lsh + 514);
    const auto *lsh_516 = buffer.data(lsh + 516);
    const auto *lsh_518 = buffer.data(lsh + 518);
    const auto *lsh_519 = buffer.data(lsh + 519);
    const auto *lsh_520 = buffer.data(lsh + 520);
    const auto *lsh_521 = buffer.data(lsh + 521);
    const auto *lsh_522 = buffer.data(lsh + 522);
    const auto *lsh_523 = buffer.data(lsh + 523);
    const auto *lsh_524 = buffer.data(lsh + 524);
    const auto *lsh_525 = buffer.data(lsh + 525);

    const auto *lsi1_420 = buffer.data(lsi1 + 420);
    const auto *lsi1_423 = buffer.data(lsi1 + 423);
    const auto *lsi1_426 = buffer.data(lsi1 + 426);
    const auto *lsi1_430 = buffer.data(lsi1 + 430);
    const auto *lsi1_432 = buffer.data(lsi1 + 432);
    const auto *lsi1_441 = buffer.data(lsi1 + 441);

    const auto *msg0_315 = buffer.data(msg0 + 315);
    const auto *msg0_317 = buffer.data(msg0 + 317);
    const auto *msg0_318 = buffer.data(msg0 + 318);
    const auto *msg0_320 = buffer.data(msg0 + 320);
    const auto *msg0_321 = buffer.data(msg0 + 321);
    const auto *msg0_325 = buffer.data(msg0 + 325);
    const auto *msg0_326 = buffer.data(msg0 + 326);
    const auto *msg0_327 = buffer.data(msg0 + 327);
    const auto *msg0_329 = buffer.data(msg0 + 329);
    const auto *msg0_335 = buffer.data(msg0 + 335);
    const auto *msg0_339 = buffer.data(msg0 + 339);
    const auto *msg0_342 = buffer.data(msg0 + 342);
    const auto *msg0_343 = buffer.data(msg0 + 343);
    const auto *msg0_344 = buffer.data(msg0 + 344);
    const auto *msg0_345 = buffer.data(msg0 + 345);
    const auto *msg0_348 = buffer.data(msg0 + 348);
    const auto *msg0_350 = buffer.data(msg0 + 350);
    const auto *msg0_351 = buffer.data(msg0 + 351);
    const auto *msg0_354 = buffer.data(msg0 + 354);
    const auto *msg0_355 = buffer.data(msg0 + 355);
    const auto *msg0_357 = buffer.data(msg0 + 357);
    const auto *msg0_358 = buffer.data(msg0 + 358);
    const auto *msg0_359 = buffer.data(msg0 + 359);
    const auto *msg0_360 = buffer.data(msg0 + 360);
    const auto *msg0_363 = buffer.data(msg0 + 363);
    const auto *msg0_365 = buffer.data(msg0 + 365);
    const auto *msg0_366 = buffer.data(msg0 + 366);
    const auto *msg0_369 = buffer.data(msg0 + 369);
    const auto *msg0_370 = buffer.data(msg0 + 370);
    const auto *msg0_372 = buffer.data(msg0 + 372);
    const auto *msg0_373 = buffer.data(msg0 + 373);
    const auto *msg0_374 = buffer.data(msg0 + 374);
    const auto *msg0_375 = buffer.data(msg0 + 375);

    const auto *msg1_315 = buffer.data(msg1 + 315);
    const auto *msg1_317 = buffer.data(msg1 + 317);
    const auto *msg1_318 = buffer.data(msg1 + 318);
    const auto *msg1_320 = buffer.data(msg1 + 320);
    const auto *msg1_321 = buffer.data(msg1 + 321);
    const auto *msg1_325 = buffer.data(msg1 + 325);
    const auto *msg1_326 = buffer.data(msg1 + 326);
    const auto *msg1_327 = buffer.data(msg1 + 327);
    const auto *msg1_329 = buffer.data(msg1 + 329);
    const auto *msg1_335 = buffer.data(msg1 + 335);
    const auto *msg1_339 = buffer.data(msg1 + 339);
    const auto *msg1_342 = buffer.data(msg1 + 342);
    const auto *msg1_343 = buffer.data(msg1 + 343);
    const auto *msg1_344 = buffer.data(msg1 + 344);
    const auto *msg1_345 = buffer.data(msg1 + 345);
    const auto *msg1_348 = buffer.data(msg1 + 348);
    const auto *msg1_350 = buffer.data(msg1 + 350);
    const auto *msg1_351 = buffer.data(msg1 + 351);
    const auto *msg1_354 = buffer.data(msg1 + 354);
    const auto *msg1_355 = buffer.data(msg1 + 355);
    const auto *msg1_357 = buffer.data(msg1 + 357);
    const auto *msg1_358 = buffer.data(msg1 + 358);
    const auto *msg1_359 = buffer.data(msg1 + 359);
    const auto *msg1_360 = buffer.data(msg1 + 360);
    const auto *msg1_363 = buffer.data(msg1 + 363);
    const auto *msg1_365 = buffer.data(msg1 + 365);
    const auto *msg1_366 = buffer.data(msg1 + 366);
    const auto *msg1_369 = buffer.data(msg1 + 369);
    const auto *msg1_370 = buffer.data(msg1 + 370);
    const auto *msg1_372 = buffer.data(msg1 + 372);
    const auto *msg1_373 = buffer.data(msg1 + 373);
    const auto *msg1_374 = buffer.data(msg1 + 374);
    const auto *msg1_375 = buffer.data(msg1 + 375);

    const auto *msh_442 = buffer.data(msh + 442);
    const auto *msh_443 = buffer.data(msh + 443);
    const auto *msh_444 = buffer.data(msh + 444);
    const auto *msh_446 = buffer.data(msh + 446);
    const auto *msh_447 = buffer.data(msh + 447);
    const auto *msh_448 = buffer.data(msh + 448);
    const auto *msh_450 = buffer.data(msh + 450);
    const auto *msh_451 = buffer.data(msh + 451);
    const auto *msh_456 = buffer.data(msh + 456);
    const auto *msh_457 = buffer.data(msh + 457);
    const auto *msh_458 = buffer.data(msh + 458);
    const auto *msh_459 = buffer.data(msh + 459);
    const auto *msh_460 = buffer.data(msh + 460);
    const auto *msh_461 = buffer.data(msh + 461);
    const auto *msh_462 = buffer.data(msh + 462);
    const auto *msh_464 = buffer.data(msh + 464);
    const auto *msh_465 = buffer.data(msh + 465);
    const auto *msh_467 = buffer.data(msh + 467);
    const auto *msh_468 = buffer.data(msh + 468);
    const auto *msh_471 = buffer.data(msh + 471);
    const auto *msh_476 = buffer.data(msh + 476);
    const auto *msh_477 = buffer.data(msh + 477);
    const auto *msh_478 = buffer.data(msh + 478);
    const auto *msh_479 = buffer.data(msh + 479);
    const auto *msh_480 = buffer.data(msh + 480);
    const auto *msh_481 = buffer.data(msh + 481);
    const auto *msh_482 = buffer.data(msh + 482);
    const auto *msh_483 = buffer.data(msh + 483);
    const auto *msh_485 = buffer.data(msh + 485);
    const auto *msh_486 = buffer.data(msh + 486);
    const auto *msh_488 = buffer.data(msh + 488);
    const auto *msh_489 = buffer.data(msh + 489);
    const auto *msh_492 = buffer.data(msh + 492);
    const auto *msh_493 = buffer.data(msh + 493);
    const auto *msh_495 = buffer.data(msh + 495);
    const auto *msh_497 = buffer.data(msh + 497);
    const auto *msh_498 = buffer.data(msh + 498);
    const auto *msh_499 = buffer.data(msh + 499);
    const auto *msh_500 = buffer.data(msh + 500);
    const auto *msh_501 = buffer.data(msh + 501);
    const auto *msh_502 = buffer.data(msh + 502);
    const auto *msh_503 = buffer.data(msh + 503);
    const auto *msh_504 = buffer.data(msh + 504);
    const auto *msh_506 = buffer.data(msh + 506);
    const auto *msh_507 = buffer.data(msh + 507);
    const auto *msh_509 = buffer.data(msh + 509);
    const auto *msh_510 = buffer.data(msh + 510);
    const auto *msh_513 = buffer.data(msh + 513);
    const auto *msh_514 = buffer.data(msh + 514);
    const auto *msh_516 = buffer.data(msh + 516);
    const auto *msh_518 = buffer.data(msh + 518);
    const auto *msh_519 = buffer.data(msh + 519);
    const auto *msh_520 = buffer.data(msh + 520);
    const auto *msh_521 = buffer.data(msh + 521);
    const auto *msh_522 = buffer.data(msh + 522);
    const auto *msh_523 = buffer.data(msh + 523);
    const auto *msh_524 = buffer.data(msh + 524);
    const auto *msh_525 = buffer.data(msh + 525);

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, lsh_447, msg0_315, msg0_321, \
                         msg1_315, msg1_321, msh_442, msh_443, msh_444, \
                         msh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * msh_442[k];

        t_593[k] = f_4 * msg0_315[k]
                   - f_5 * msg1_315[k]
                   + f_3 * pc_z[k] * msh_443[k];

        t_594[k] = f_13 * lsh_447[k]
                   + f_6 * msg0_321[k]
                   - f_7 * msg1_321[k]
                   + f_3 * pc_x[k] * msh_447[k];

        t_595[k] = f_3 * pc_z[k] * msh_444[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, lsh_320, lsh_451, \
                         msg0_317, msg0_325, msg1_317, msg1_325, msh_446, msh_447, \
                         msh_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_19 * lsh_320[k]
                   + f_3 * pc_y[k] * msh_446[k];

        t_597[k] = f_6 * msg0_317[k]
                   - f_7 * msg1_317[k]
                   + f_3 * pc_z[k] * msh_446[k];

        t_598[k] = f_13 * lsh_451[k]
                   + f_4 * msg0_325[k]
                   - f_5 * msg1_325[k]
                   + f_3 * pc_x[k] * msh_451[k];

        t_599[k] = f_3 * pc_z[k] * msh_447[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, pc_z, lsh_324, lsh_456, \
                         msg0_318, msg0_320, msg1_318, msg1_320, msh_448, msh_450, \
                         msh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_4 * msg0_318[k]
                   - f_5 * msg1_318[k]
                   + f_3 * pc_z[k] * msh_448[k];

        t_601[k] = f_19 * lsh_324[k]
                   + f_3 * pc_y[k] * msh_450[k];

        t_602[k] = f_8 * msg0_320[k]
                   - f_9 * msg1_320[k]
                   + f_3 * pc_z[k] * msh_450[k];

        t_603[k] = f_13 * lsh_456[k]
                   + f_3 * pc_x[k] * msh_456[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pc_x, pc_z, lsh_458, lsh_459, \
                         lsh_460, lsh_461, msh_451, msh_458, msh_459, msh_460, \
                         msh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * msh_451[k];

        t_605[k] = f_13 * lsh_458[k]
                   + f_3 * pc_x[k] * msh_458[k];

        t_606[k] = f_13 * lsh_459[k]
                   + f_3 * pc_x[k] * msh_459[k];

        t_607[k] = f_13 * lsh_460[k]
                   + f_3 * pc_x[k] * msh_460[k];

        t_608[k] = f_13 * lsh_461[k]
                   + f_3 * pc_x[k] * msh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pc_y, pc_z, lsh_330, msg0_325, msg0_326, \
                         msg1_325, msg1_326, msh_456, msh_457, \
                         msh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_19 * lsh_330[k]
                   + f_1 * msg0_325[k]
                   - f_2 * msg1_325[k]
                   + f_3 * pc_y[k] * msh_456[k];

        t_610[k] = f_3 * pc_z[k] * msh_456[k];

        t_611[k] = f_4 * msg0_325[k]
                   - f_5 * msg1_325[k]
                   + f_3 * pc_z[k] * msh_457[k];

        t_612[k] = f_6 * msg0_326[k]
                   - f_7 * msg1_326[k]
                   + f_3 * pc_z[k] * msh_458[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pc_y, pc_z, lsi0_420, lsh_335, \
                         lsi1_420, msg0_327, msg0_329, msg1_327, msg1_329, msh_459, \
                         msh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_8 * msg0_327[k]
                   - f_9 * msg1_327[k]
                   + f_3 * pc_z[k] * msh_459[k];

        t_614[k] = f_19 * lsh_335[k]
                   + f_3 * pc_y[k] * msh_461[k];

        t_615[k] = f_1 * msg0_329[k]
                   - f_2 * msg1_329[k]
                   + f_3 * pc_z[k] * msh_461[k];

        t_616[k] = pa_z[k] * lsi0_420[k]
                   - f_10 * pc_z[k] * lsi1_420[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pa_z, pc_y, pc_z, lsi0_423, lsh_315, \
                         lsh_336, lsh_338, lsi1_423, msh_462, msh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_20 * lsh_336[k]
                   + f_3 * pc_y[k] * msh_462[k];

        t_618[k] = f_11 * lsh_315[k]
                   + f_3 * pc_z[k] * msh_462[k];

        t_619[k] = pa_z[k] * lsi0_423[k]
                   - f_10 * pc_z[k] * lsi1_423[k];

        t_620[k] = f_20 * lsh_338[k]
                   + f_3 * pc_y[k] * msh_464[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pa_z, pc_x, pc_z, lsi0_426, lsh_318, lsh_467, \
                         lsi1_426, msg0_335, msg1_335, msh_465, \
                         msh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_13 * lsh_467[k]
                   + f_8 * msg0_335[k]
                   - f_9 * msg1_335[k]
                   + f_3 * pc_x[k] * msh_467[k];

        t_622[k] = pa_z[k] * lsi0_426[k]
                   - f_10 * pc_z[k] * lsi1_426[k];

        t_623[k] = f_11 * lsh_318[k]
                   + f_3 * pc_z[k] * msh_465[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_z, pc_x, pc_y, pc_z, lsi0_430, lsh_341, \
                         lsh_471, lsi1_430, msg0_339, msg1_339, msh_467, \
                         msh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_20 * lsh_341[k]
                   + f_3 * pc_y[k] * msh_467[k];

        t_625[k] = f_13 * lsh_471[k]
                   + f_6 * msg0_339[k]
                   - f_7 * msg1_339[k]
                   + f_3 * pc_x[k] * msh_471[k];

        t_626[k] = pa_z[k] * lsi0_430[k]
                   - f_10 * pc_z[k] * lsi1_430[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_z, pc_y, pc_z, lsi0_432, lsh_321, lsh_322, \
                         lsh_345, lsi1_432, msh_468, msh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_11 * lsh_321[k]
                   + f_3 * pc_z[k] * msh_468[k];

        t_628[k] = pa_z[k] * lsi0_432[k]
                   + f_12 * lsh_322[k]
                   - f_10 * pc_z[k] * lsi1_432[k];

        t_629[k] = f_20 * lsh_345[k]
                   + f_3 * pc_y[k] * msh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, lsh_476, lsh_477, lsh_478, lsh_479, \
                         msg0_344, msg1_344, msh_476, msh_477, msh_478, \
                         msh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_13 * lsh_476[k]
                   + f_4 * msg0_344[k]
                   - f_5 * msg1_344[k]
                   + f_3 * pc_x[k] * msh_476[k];

        t_631[k] = f_13 * lsh_477[k]
                   + f_3 * pc_x[k] * msh_477[k];

        t_632[k] = f_13 * lsh_478[k]
                   + f_3 * pc_x[k] * msh_478[k];

        t_633[k] = f_13 * lsh_479[k]
                   + f_3 * pc_x[k] * msh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_z, pc_x, pc_z, lsi0_441, lsh_480, \
                         lsh_481, lsh_482, lsi1_441, msh_480, msh_481, \
                         msh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_13 * lsh_480[k]
                   + f_3 * pc_x[k] * msh_480[k];

        t_635[k] = f_13 * lsh_481[k]
                   + f_3 * pc_x[k] * msh_481[k];

        t_636[k] = f_13 * lsh_482[k]
                   + f_3 * pc_x[k] * msh_482[k];

        t_637[k] = pa_z[k] * lsi0_441[k]
                   - f_10 * pc_z[k] * lsi1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, lsh_330, lsh_353, lsh_354, msg0_342, \
                         msg0_343, msg1_342, msg1_343, msh_477, msh_479, \
                         msh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * lsh_330[k]
                   + f_3 * pc_z[k] * msh_477[k];

        t_639[k] = f_20 * lsh_353[k]
                   + f_8 * msg0_342[k]
                   - f_9 * msg1_342[k]
                   + f_3 * pc_y[k] * msh_479[k];

        t_640[k] = f_20 * lsh_354[k]
                   + f_6 * msg0_343[k]
                   - f_7 * msg1_343[k]
                   + f_3 * pc_y[k] * msh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, lsh_335, lsh_355, lsh_356, msg0_344, \
                         msg1_344, msh_481, msh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_20 * lsh_355[k]
                   + f_4 * msg0_344[k]
                   - f_5 * msg1_344[k]
                   + f_3 * pc_y[k] * msh_481[k];

        t_642[k] = f_20 * lsh_356[k]
                   + f_3 * pc_y[k] * msh_482[k];

        t_643[k] = f_11 * lsh_335[k]
                   + f_1 * msg0_344[k]
                   - f_2 * msg1_344[k]
                   + f_3 * pc_z[k] * msh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, lsh_336, lsh_357, lsh_483, \
                         msg0_345, msg1_345, msh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_13 * lsh_483[k]
                   + f_1 * msg0_345[k]
                   - f_2 * msg1_345[k]
                   + f_3 * pc_x[k] * msh_483[k];

        t_645[k] = f_14 * lsh_357[k]
                   + f_3 * pc_y[k] * msh_483[k];

        t_646[k] = f_12 * lsh_336[k]
                   + f_3 * pc_z[k] * msh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, lsh_359, lsh_486, lsh_488, msg0_348, \
                         msg0_350, msg1_348, msg1_350, msh_485, msh_486, \
                         msh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_13 * lsh_486[k]
                   + f_8 * msg0_348[k]
                   - f_9 * msg1_348[k]
                   + f_3 * pc_x[k] * msh_486[k];

        t_648[k] = f_14 * lsh_359[k]
                   + f_3 * pc_y[k] * msh_485[k];

        t_649[k] = f_13 * lsh_488[k]
                   + f_8 * msg0_350[k]
                   - f_9 * msg1_350[k]
                   + f_3 * pc_x[k] * msh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, lsh_339, lsh_362, lsh_489, \
                         msg0_351, msg1_351, msh_486, msh_488, \
                         msh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_13 * lsh_489[k]
                   + f_6 * msg0_351[k]
                   - f_7 * msg1_351[k]
                   + f_3 * pc_x[k] * msh_489[k];

        t_651[k] = f_12 * lsh_339[k]
                   + f_3 * pc_z[k] * msh_486[k];

        t_652[k] = f_14 * lsh_362[k]
                   + f_3 * pc_y[k] * msh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, lsh_342, lsh_492, lsh_493, msg0_354, \
                         msg0_355, msg1_354, msg1_355, msh_489, msh_492, \
                         msh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_13 * lsh_492[k]
                   + f_6 * msg0_354[k]
                   - f_7 * msg1_354[k]
                   + f_3 * pc_x[k] * msh_492[k];

        t_654[k] = f_13 * lsh_493[k]
                   + f_4 * msg0_355[k]
                   - f_5 * msg1_355[k]
                   + f_3 * pc_x[k] * msh_493[k];

        t_655[k] = f_12 * lsh_342[k]
                   + f_3 * pc_z[k] * msh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, lsh_366, lsh_495, lsh_497, msg0_357, \
                         msg0_359, msg1_357, msg1_359, msh_492, msh_495, \
                         msh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_13 * lsh_495[k]
                   + f_4 * msg0_357[k]
                   - f_5 * msg1_357[k]
                   + f_3 * pc_x[k] * msh_495[k];

        t_657[k] = f_14 * lsh_366[k]
                   + f_3 * pc_y[k] * msh_492[k];

        t_658[k] = f_13 * lsh_497[k]
                   + f_4 * msg0_359[k]
                   - f_5 * msg1_359[k]
                   + f_3 * pc_x[k] * msh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, lsh_498, lsh_499, lsh_500, \
                         lsh_501, lsh_502, msh_498, msh_499, msh_500, msh_501, \
                         msh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_13 * lsh_498[k]
                   + f_3 * pc_x[k] * msh_498[k];

        t_660[k] = f_13 * lsh_499[k]
                   + f_3 * pc_x[k] * msh_499[k];

        t_661[k] = f_13 * lsh_500[k]
                   + f_3 * pc_x[k] * msh_500[k];

        t_662[k] = f_13 * lsh_501[k]
                   + f_3 * pc_x[k] * msh_501[k];

        t_663[k] = f_13 * lsh_502[k]
                   + f_3 * pc_x[k] * msh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, lsh_351, lsh_372, lsh_503, \
                         msg0_355, msg1_355, msh_498, msh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_13 * lsh_503[k]
                   + f_3 * pc_x[k] * msh_503[k];

        t_665[k] = f_14 * lsh_372[k]
                   + f_1 * msg0_355[k]
                   - f_2 * msg1_355[k]
                   + f_3 * pc_y[k] * msh_498[k];

        t_666[k] = f_12 * lsh_351[k]
                   + f_3 * pc_z[k] * msh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, lsh_374, lsh_375, lsh_376, msg0_357, \
                         msg0_358, msg0_359, msg1_357, msg1_358, msg1_359, msh_500, msh_501, \
                         msh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * lsh_374[k]
                   + f_8 * msg0_357[k]
                   - f_9 * msg1_357[k]
                   + f_3 * pc_y[k] * msh_500[k];

        t_668[k] = f_14 * lsh_375[k]
                   + f_6 * msg0_358[k]
                   - f_7 * msg1_358[k]
                   + f_3 * pc_y[k] * msh_501[k];

        t_669[k] = f_14 * lsh_376[k]
                   + f_4 * msg0_359[k]
                   - f_5 * msg1_359[k]
                   + f_3 * pc_y[k] * msh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, lsh_356, lsh_377, lsh_504, \
                         msg0_359, msg0_360, msg1_359, msg1_360, msh_503, \
                         msh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * lsh_377[k]
                   + f_3 * pc_y[k] * msh_503[k];

        t_671[k] = f_12 * lsh_356[k]
                   + f_1 * msg0_359[k]
                   - f_2 * msg1_359[k]
                   + f_3 * pc_z[k] * msh_503[k];

        t_672[k] = f_13 * lsh_504[k]
                   + f_1 * msg0_360[k]
                   - f_2 * msg1_360[k]
                   + f_3 * pc_x[k] * msh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, lsh_357, lsh_378, \
                         lsh_380, lsh_507, msg0_363, msg1_363, msh_504, msh_506, \
                         msh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * lsh_378[k]
                   + f_3 * pc_y[k] * msh_504[k];

        t_674[k] = f_13 * lsh_357[k]
                   + f_3 * pc_z[k] * msh_504[k];

        t_675[k] = f_13 * lsh_507[k]
                   + f_8 * msg0_363[k]
                   - f_9 * msg1_363[k]
                   + f_3 * pc_x[k] * msh_507[k];

        t_676[k] = f_13 * lsh_380[k]
                   + f_3 * pc_y[k] * msh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, lsh_360, lsh_509, lsh_510, msg0_365, \
                         msg0_366, msg1_365, msg1_366, msh_507, msh_509, \
                         msh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_13 * lsh_509[k]
                   + f_8 * msg0_365[k]
                   - f_9 * msg1_365[k]
                   + f_3 * pc_x[k] * msh_509[k];

        t_678[k] = f_13 * lsh_510[k]
                   + f_6 * msg0_366[k]
                   - f_7 * msg1_366[k]
                   + f_3 * pc_x[k] * msh_510[k];

        t_679[k] = f_13 * lsh_360[k]
                   + f_3 * pc_z[k] * msh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, lsh_383, lsh_513, lsh_514, msg0_369, \
                         msg0_370, msg1_369, msg1_370, msh_509, msh_513, \
                         msh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * lsh_383[k]
                   + f_3 * pc_y[k] * msh_509[k];

        t_681[k] = f_13 * lsh_513[k]
                   + f_6 * msg0_369[k]
                   - f_7 * msg1_369[k]
                   + f_3 * pc_x[k] * msh_513[k];

        t_682[k] = f_13 * lsh_514[k]
                   + f_4 * msg0_370[k]
                   - f_5 * msg1_370[k]
                   + f_3 * pc_x[k] * msh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, lsh_363, lsh_387, lsh_516, \
                         msg0_372, msg1_372, msh_510, msh_513, \
                         msh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * lsh_363[k]
                   + f_3 * pc_z[k] * msh_510[k];

        t_684[k] = f_13 * lsh_516[k]
                   + f_4 * msg0_372[k]
                   - f_5 * msg1_372[k]
                   + f_3 * pc_x[k] * msh_516[k];

        t_685[k] = f_13 * lsh_387[k]
                   + f_3 * pc_y[k] * msh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, lsh_518, lsh_519, lsh_520, lsh_521, \
                         msg0_374, msg1_374, msh_518, msh_519, msh_520, \
                         msh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_13 * lsh_518[k]
                   + f_4 * msg0_374[k]
                   - f_5 * msg1_374[k]
                   + f_3 * pc_x[k] * msh_518[k];

        t_687[k] = f_13 * lsh_519[k]
                   + f_3 * pc_x[k] * msh_519[k];

        t_688[k] = f_13 * lsh_520[k]
                   + f_3 * pc_x[k] * msh_520[k];

        t_689[k] = f_13 * lsh_521[k]
                   + f_3 * pc_x[k] * msh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, lsh_393, lsh_522, lsh_523, \
                         lsh_524, msg0_370, msg1_370, msh_519, msh_522, msh_523, \
                         msh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_13 * lsh_522[k]
                   + f_3 * pc_x[k] * msh_522[k];

        t_691[k] = f_13 * lsh_523[k]
                   + f_3 * pc_x[k] * msh_523[k];

        t_692[k] = f_13 * lsh_524[k]
                   + f_3 * pc_x[k] * msh_524[k];

        t_693[k] = f_13 * lsh_393[k]
                   + f_1 * msg0_370[k]
                   - f_2 * msg1_370[k]
                   + f_3 * pc_y[k] * msh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, lsh_372, lsh_395, lsh_396, msg0_372, \
                         msg0_373, msg1_372, msg1_373, msh_519, msh_521, \
                         msh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * lsh_372[k]
                   + f_3 * pc_z[k] * msh_519[k];

        t_695[k] = f_13 * lsh_395[k]
                   + f_8 * msg0_372[k]
                   - f_9 * msg1_372[k]
                   + f_3 * pc_y[k] * msh_521[k];

        t_696[k] = f_13 * lsh_396[k]
                   + f_6 * msg0_373[k]
                   - f_7 * msg1_373[k]
                   + f_3 * pc_y[k] * msh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, lsh_377, lsh_397, lsh_398, msg0_374, \
                         msg1_374, msh_523, msh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * lsh_397[k]
                   + f_4 * msg0_374[k]
                   - f_5 * msg1_374[k]
                   + f_3 * pc_y[k] * msh_523[k];

        t_698[k] = f_13 * lsh_398[k]
                   + f_3 * pc_y[k] * msh_524[k];

        t_699[k] = f_13 * lsh_377[k]
                   + f_1 * msg0_374[k]
                   - f_2 * msg1_374[k]
                   + f_3 * pc_z[k] * msh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, lsh_378, lsh_399, lsh_525, \
                         msg0_375, msg1_375, msh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_13 * lsh_525[k]
                   + f_1 * msg0_375[k]
                   - f_2 * msg1_375[k]
                   + f_3 * pc_x[k] * msh_525[k];

        t_701[k] = f_12 * lsh_399[k]
                   + f_3 * pc_y[k] * msh_525[k];

        t_702[k] = f_14 * lsh_378[k]
                   + f_3 * pc_z[k] * msh_525[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *lsi0_560 = buffer.data(lsi0 + 560);
    const auto *lsi0_563 = buffer.data(lsi0 + 563);
    const auto *lsi0_565 = buffer.data(lsi0 + 565);
    const auto *lsi0_566 = buffer.data(lsi0 + 566);
    const auto *lsi0_569 = buffer.data(lsi0 + 569);
    const auto *lsi0_570 = buffer.data(lsi0 + 570);
    const auto *lsi0_572 = buffer.data(lsi0 + 572);
    const auto *lsi0_574 = buffer.data(lsi0 + 574);
    const auto *lsi0_587 = buffer.data(lsi0 + 587);
    const auto *lsi0_588 = buffer.data(lsi0 + 588);
    const auto *lsi0_591 = buffer.data(lsi0 + 591);
    const auto *lsi0_594 = buffer.data(lsi0 + 594);

    const auto *lsh_381 = buffer.data(lsh + 381);
    const auto *lsh_384 = buffer.data(lsh + 384);
    const auto *lsh_393 = buffer.data(lsh + 393);
    const auto *lsh_398 = buffer.data(lsh + 398);
    const auto *lsh_399 = buffer.data(lsh + 399);
    const auto *lsh_401 = buffer.data(lsh + 401);
    const auto *lsh_402 = buffer.data(lsh + 402);
    const auto *lsh_404 = buffer.data(lsh + 404);
    const auto *lsh_405 = buffer.data(lsh + 405);
    const auto *lsh_408 = buffer.data(lsh + 408);
    const auto *lsh_414 = buffer.data(lsh + 414);
    const auto *lsh_416 = buffer.data(lsh + 416);
    const auto *lsh_417 = buffer.data(lsh + 417);
    const auto *lsh_418 = buffer.data(lsh + 418);
    const auto *lsh_419 = buffer.data(lsh + 419);
    const auto *lsh_420 = buffer.data(lsh + 420);
    const auto *lsh_421 = buffer.data(lsh + 421);
    const auto *lsh_422 = buffer.data(lsh + 422);
    const auto *lsh_423 = buffer.data(lsh + 423);
    const auto *lsh_425 = buffer.data(lsh + 425);
    const auto *lsh_426 = buffer.data(lsh + 426);
    const auto *lsh_428 = buffer.data(lsh + 428);
    const auto *lsh_429 = buffer.data(lsh + 429);
    const auto *lsh_435 = buffer.data(lsh + 435);
    const auto *lsh_437 = buffer.data(lsh + 437);
    const auto *lsh_438 = buffer.data(lsh + 438);
    const auto *lsh_439 = buffer.data(lsh + 439);
    const auto *lsh_440 = buffer.data(lsh + 440);
    const auto *lsh_441 = buffer.data(lsh + 441);
    const auto *lsh_444 = buffer.data(lsh + 444);
    const auto *lsh_446 = buffer.data(lsh + 446);
    const auto *lsh_450 = buffer.data(lsh + 450);
    const auto *lsh_456 = buffer.data(lsh + 456);
    const auto *lsh_461 = buffer.data(lsh + 461);
    const auto *lsh_462 = buffer.data(lsh + 462);
    const auto *lsh_464 = buffer.data(lsh + 464);
    const auto *lsh_528 = buffer.data(lsh + 528);
    const auto *lsh_530 = buffer.data(lsh + 530);
    const auto *lsh_531 = buffer.data(lsh + 531);
    const auto *lsh_534 = buffer.data(lsh + 534);
    const auto *lsh_535 = buffer.data(lsh + 535);
    const auto *lsh_537 = buffer.data(lsh + 537);
    const auto *lsh_539 = buffer.data(lsh + 539);
    const auto *lsh_540 = buffer.data(lsh + 540);
    const auto *lsh_541 = buffer.data(lsh + 541);
    const auto *lsh_542 = buffer.data(lsh + 542);
    const auto *lsh_543 = buffer.data(lsh + 543);
    const auto *lsh_544 = buffer.data(lsh + 544);
    const auto *lsh_545 = buffer.data(lsh + 545);
    const auto *lsh_561 = buffer.data(lsh + 561);
    const auto *lsh_562 = buffer.data(lsh + 562);
    const auto *lsh_563 = buffer.data(lsh + 563);
    const auto *lsh_564 = buffer.data(lsh + 564);
    const auto *lsh_565 = buffer.data(lsh + 565);
    const auto *lsh_566 = buffer.data(lsh + 566);
    const auto *lsh_567 = buffer.data(lsh + 567);
    const auto *lsh_572 = buffer.data(lsh + 572);
    const auto *lsh_576 = buffer.data(lsh + 576);
    const auto *lsh_581 = buffer.data(lsh + 581);
    const auto *lsh_582 = buffer.data(lsh + 582);
    const auto *lsh_583 = buffer.data(lsh + 583);
    const auto *lsh_584 = buffer.data(lsh + 584);
    const auto *lsh_585 = buffer.data(lsh + 585);
    const auto *lsh_587 = buffer.data(lsh + 587);
    const auto *lsh_588 = buffer.data(lsh + 588);
    const auto *lsh_591 = buffer.data(lsh + 591);
    const auto *lsh_594 = buffer.data(lsh + 594);
    const auto *lsh_598 = buffer.data(lsh + 598);
    const auto *lsh_603 = buffer.data(lsh + 603);
    const auto *lsh_605 = buffer.data(lsh + 605);
    const auto *lsh_606 = buffer.data(lsh + 606);
    const auto *lsh_607 = buffer.data(lsh + 607);
    const auto *lsh_608 = buffer.data(lsh + 608);
    const auto *lsh_614 = buffer.data(lsh + 614);

    const auto *lsi1_560 = buffer.data(lsi1 + 560);
    const auto *lsi1_563 = buffer.data(lsi1 + 563);
    const auto *lsi1_565 = buffer.data(lsi1 + 565);
    const auto *lsi1_566 = buffer.data(lsi1 + 566);
    const auto *lsi1_569 = buffer.data(lsi1 + 569);
    const auto *lsi1_570 = buffer.data(lsi1 + 570);
    const auto *lsi1_572 = buffer.data(lsi1 + 572);
    const auto *lsi1_574 = buffer.data(lsi1 + 574);
    const auto *lsi1_587 = buffer.data(lsi1 + 587);
    const auto *lsi1_588 = buffer.data(lsi1 + 588);
    const auto *lsi1_591 = buffer.data(lsi1 + 591);
    const auto *lsi1_594 = buffer.data(lsi1 + 594);

    const auto *msg0_378 = buffer.data(msg0 + 378);
    const auto *msg0_380 = buffer.data(msg0 + 380);
    const auto *msg0_381 = buffer.data(msg0 + 381);
    const auto *msg0_384 = buffer.data(msg0 + 384);
    const auto *msg0_385 = buffer.data(msg0 + 385);
    const auto *msg0_387 = buffer.data(msg0 + 387);
    const auto *msg0_388 = buffer.data(msg0 + 388);
    const auto *msg0_389 = buffer.data(msg0 + 389);
    const auto *msg0_400 = buffer.data(msg0 + 400);
    const auto *msg0_402 = buffer.data(msg0 + 402);
    const auto *msg0_403 = buffer.data(msg0 + 403);
    const auto *msg0_404 = buffer.data(msg0 + 404);
    const auto *msg0_405 = buffer.data(msg0 + 405);
    const auto *msg0_406 = buffer.data(msg0 + 406);
    const auto *msg0_407 = buffer.data(msg0 + 407);
    const auto *msg0_408 = buffer.data(msg0 + 408);
    const auto *msg0_409 = buffer.data(msg0 + 409);
    const auto *msg0_410 = buffer.data(msg0 + 410);
    const auto *msg0_414 = buffer.data(msg0 + 414);
    const auto *msg0_415 = buffer.data(msg0 + 415);
    const auto *msg0_416 = buffer.data(msg0 + 416);
    const auto *msg0_417 = buffer.data(msg0 + 417);
    const auto *msg0_418 = buffer.data(msg0 + 418);
    const auto *msg0_419 = buffer.data(msg0 + 419);
    const auto *msg0_420 = buffer.data(msg0 + 420);
    const auto *msg0_422 = buffer.data(msg0 + 422);
    const auto *msg0_423 = buffer.data(msg0 + 423);
    const auto *msg0_425 = buffer.data(msg0 + 425);
    const auto *msg0_426 = buffer.data(msg0 + 426);
    const auto *msg0_430 = buffer.data(msg0 + 430);
    const auto *msg0_431 = buffer.data(msg0 + 431);
    const auto *msg0_432 = buffer.data(msg0 + 432);
    const auto *msg0_434 = buffer.data(msg0 + 434);
    const auto *msg0_440 = buffer.data(msg0 + 440);

    const auto *msg1_378 = buffer.data(msg1 + 378);
    const auto *msg1_380 = buffer.data(msg1 + 380);
    const auto *msg1_381 = buffer.data(msg1 + 381);
    const auto *msg1_384 = buffer.data(msg1 + 384);
    const auto *msg1_385 = buffer.data(msg1 + 385);
    const auto *msg1_387 = buffer.data(msg1 + 387);
    const auto *msg1_388 = buffer.data(msg1 + 388);
    const auto *msg1_389 = buffer.data(msg1 + 389);
    const auto *msg1_400 = buffer.data(msg1 + 400);
    const auto *msg1_402 = buffer.data(msg1 + 402);
    const auto *msg1_403 = buffer.data(msg1 + 403);
    const auto *msg1_404 = buffer.data(msg1 + 404);
    const auto *msg1_405 = buffer.data(msg1 + 405);
    const auto *msg1_406 = buffer.data(msg1 + 406);
    const auto *msg1_407 = buffer.data(msg1 + 407);
    const auto *msg1_408 = buffer.data(msg1 + 408);
    const auto *msg1_409 = buffer.data(msg1 + 409);
    const auto *msg1_410 = buffer.data(msg1 + 410);
    const auto *msg1_414 = buffer.data(msg1 + 414);
    const auto *msg1_415 = buffer.data(msg1 + 415);
    const auto *msg1_416 = buffer.data(msg1 + 416);
    const auto *msg1_417 = buffer.data(msg1 + 417);
    const auto *msg1_418 = buffer.data(msg1 + 418);
    const auto *msg1_419 = buffer.data(msg1 + 419);
    const auto *msg1_420 = buffer.data(msg1 + 420);
    const auto *msg1_422 = buffer.data(msg1 + 422);
    const auto *msg1_423 = buffer.data(msg1 + 423);
    const auto *msg1_425 = buffer.data(msg1 + 425);
    const auto *msg1_426 = buffer.data(msg1 + 426);
    const auto *msg1_430 = buffer.data(msg1 + 430);
    const auto *msg1_431 = buffer.data(msg1 + 431);
    const auto *msg1_432 = buffer.data(msg1 + 432);
    const auto *msg1_434 = buffer.data(msg1 + 434);
    const auto *msg1_440 = buffer.data(msg1 + 440);

    const auto *msh_527 = buffer.data(msh + 527);
    const auto *msh_528 = buffer.data(msh + 528);
    const auto *msh_530 = buffer.data(msh + 530);
    const auto *msh_531 = buffer.data(msh + 531);
    const auto *msh_534 = buffer.data(msh + 534);
    const auto *msh_535 = buffer.data(msh + 535);
    const auto *msh_537 = buffer.data(msh + 537);
    const auto *msh_539 = buffer.data(msh + 539);
    const auto *msh_540 = buffer.data(msh + 540);
    const auto *msh_541 = buffer.data(msh + 541);
    const auto *msh_542 = buffer.data(msh + 542);
    const auto *msh_543 = buffer.data(msh + 543);
    const auto *msh_544 = buffer.data(msh + 544);
    const auto *msh_545 = buffer.data(msh + 545);
    const auto *msh_546 = buffer.data(msh + 546);
    const auto *msh_548 = buffer.data(msh + 548);
    const auto *msh_549 = buffer.data(msh + 549);
    const auto *msh_551 = buffer.data(msh + 551);
    const auto *msh_552 = buffer.data(msh + 552);
    const auto *msh_555 = buffer.data(msh + 555);
    const auto *msh_561 = buffer.data(msh + 561);
    const auto *msh_562 = buffer.data(msh + 562);
    const auto *msh_563 = buffer.data(msh + 563);
    const auto *msh_564 = buffer.data(msh + 564);
    const auto *msh_565 = buffer.data(msh + 565);
    const auto *msh_566 = buffer.data(msh + 566);
    const auto *msh_567 = buffer.data(msh + 567);
    const auto *msh_568 = buffer.data(msh + 568);
    const auto *msh_569 = buffer.data(msh + 569);
    const auto *msh_570 = buffer.data(msh + 570);
    const auto *msh_571 = buffer.data(msh + 571);
    const auto *msh_572 = buffer.data(msh + 572);
    const auto *msh_573 = buffer.data(msh + 573);
    const auto *msh_574 = buffer.data(msh + 574);
    const auto *msh_575 = buffer.data(msh + 575);
    const auto *msh_576 = buffer.data(msh + 576);
    const auto *msh_581 = buffer.data(msh + 581);
    const auto *msh_582 = buffer.data(msh + 582);
    const auto *msh_583 = buffer.data(msh + 583);
    const auto *msh_584 = buffer.data(msh + 584);
    const auto *msh_585 = buffer.data(msh + 585);
    const auto *msh_586 = buffer.data(msh + 586);
    const auto *msh_587 = buffer.data(msh + 587);
    const auto *msh_588 = buffer.data(msh + 588);
    const auto *msh_589 = buffer.data(msh + 589);
    const auto *msh_590 = buffer.data(msh + 590);
    const auto *msh_591 = buffer.data(msh + 591);
    const auto *msh_593 = buffer.data(msh + 593);
    const auto *msh_594 = buffer.data(msh + 594);
    const auto *msh_595 = buffer.data(msh + 595);
    const auto *msh_597 = buffer.data(msh + 597);
    const auto *msh_598 = buffer.data(msh + 598);
    const auto *msh_603 = buffer.data(msh + 603);
    const auto *msh_604 = buffer.data(msh + 604);
    const auto *msh_605 = buffer.data(msh + 605);
    const auto *msh_606 = buffer.data(msh + 606);
    const auto *msh_607 = buffer.data(msh + 607);
    const auto *msh_608 = buffer.data(msh + 608);
    const auto *msh_609 = buffer.data(msh + 609);
    const auto *msh_611 = buffer.data(msh + 611);
    const auto *msh_612 = buffer.data(msh + 612);
    const auto *msh_614 = buffer.data(msh + 614);

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, lsh_401, lsh_528, lsh_530, msg0_378, \
                         msg0_380, msg1_378, msg1_380, msh_527, msh_528, \
                         msh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_13 * lsh_528[k]
                   + f_8 * msg0_378[k]
                   - f_9 * msg1_378[k]
                   + f_3 * pc_x[k] * msh_528[k];

        t_704[k] = f_12 * lsh_401[k]
                   + f_3 * pc_y[k] * msh_527[k];

        t_705[k] = f_13 * lsh_530[k]
                   + f_8 * msg0_380[k]
                   - f_9 * msg1_380[k]
                   + f_3 * pc_x[k] * msh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, lsh_381, lsh_404, lsh_531, \
                         msg0_381, msg1_381, msh_528, msh_530, \
                         msh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_13 * lsh_531[k]
                   + f_6 * msg0_381[k]
                   - f_7 * msg1_381[k]
                   + f_3 * pc_x[k] * msh_531[k];

        t_707[k] = f_14 * lsh_381[k]
                   + f_3 * pc_z[k] * msh_528[k];

        t_708[k] = f_12 * lsh_404[k]
                   + f_3 * pc_y[k] * msh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, lsh_384, lsh_534, lsh_535, msg0_384, \
                         msg0_385, msg1_384, msg1_385, msh_531, msh_534, \
                         msh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_13 * lsh_534[k]
                   + f_6 * msg0_384[k]
                   - f_7 * msg1_384[k]
                   + f_3 * pc_x[k] * msh_534[k];

        t_710[k] = f_13 * lsh_535[k]
                   + f_4 * msg0_385[k]
                   - f_5 * msg1_385[k]
                   + f_3 * pc_x[k] * msh_535[k];

        t_711[k] = f_14 * lsh_384[k]
                   + f_3 * pc_z[k] * msh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, lsh_408, lsh_537, lsh_539, msg0_387, \
                         msg0_389, msg1_387, msg1_389, msh_534, msh_537, \
                         msh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_13 * lsh_537[k]
                   + f_4 * msg0_387[k]
                   - f_5 * msg1_387[k]
                   + f_3 * pc_x[k] * msh_537[k];

        t_713[k] = f_12 * lsh_408[k]
                   + f_3 * pc_y[k] * msh_534[k];

        t_714[k] = f_13 * lsh_539[k]
                   + f_4 * msg0_389[k]
                   - f_5 * msg1_389[k]
                   + f_3 * pc_x[k] * msh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, lsh_540, lsh_541, lsh_542, \
                         lsh_543, lsh_544, msh_540, msh_541, msh_542, msh_543, \
                         msh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_13 * lsh_540[k]
                   + f_3 * pc_x[k] * msh_540[k];

        t_716[k] = f_13 * lsh_541[k]
                   + f_3 * pc_x[k] * msh_541[k];

        t_717[k] = f_13 * lsh_542[k]
                   + f_3 * pc_x[k] * msh_542[k];

        t_718[k] = f_13 * lsh_543[k]
                   + f_3 * pc_x[k] * msh_543[k];

        t_719[k] = f_13 * lsh_544[k]
                   + f_3 * pc_x[k] * msh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, lsh_393, lsh_414, lsh_545, \
                         msg0_385, msg1_385, msh_540, msh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_13 * lsh_545[k]
                   + f_3 * pc_x[k] * msh_545[k];

        t_721[k] = f_12 * lsh_414[k]
                   + f_1 * msg0_385[k]
                   - f_2 * msg1_385[k]
                   + f_3 * pc_y[k] * msh_540[k];

        t_722[k] = f_14 * lsh_393[k]
                   + f_3 * pc_z[k] * msh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, lsh_416, lsh_417, lsh_418, msg0_387, \
                         msg0_388, msg0_389, msg1_387, msg1_388, msg1_389, msh_542, msh_543, \
                         msh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * lsh_416[k]
                   + f_8 * msg0_387[k]
                   - f_9 * msg1_387[k]
                   + f_3 * pc_y[k] * msh_542[k];

        t_724[k] = f_12 * lsh_417[k]
                   + f_6 * msg0_388[k]
                   - f_7 * msg1_388[k]
                   + f_3 * pc_y[k] * msh_543[k];

        t_725[k] = f_12 * lsh_418[k]
                   + f_4 * msg0_389[k]
                   - f_5 * msg1_389[k]
                   + f_3 * pc_y[k] * msh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_y, pc_y, pc_z, lsi0_560, lsh_398, \
                         lsh_419, lsh_420, lsi1_560, msg0_389, msg1_389, msh_545, \
                         msh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * lsh_419[k]
                   + f_3 * pc_y[k] * msh_545[k];

        t_727[k] = f_14 * lsh_398[k]
                   + f_1 * msg0_389[k]
                   - f_2 * msg1_389[k]
                   + f_3 * pc_z[k] * msh_545[k];

        t_728[k] = pa_y[k] * lsi0_560[k]
                   - f_10 * pc_y[k] * lsi1_560[k];

        t_729[k] = f_11 * lsh_420[k]
                   + f_3 * pc_y[k] * msh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_y, pc_y, pc_z, lsi0_563, lsi0_565, \
                         lsh_399, lsh_421, lsh_422, lsi1_563, lsi1_565, msh_546, \
                         msh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_20 * lsh_399[k]
                   + f_3 * pc_z[k] * msh_546[k];

        t_731[k] = pa_y[k] * lsi0_563[k]
                   + f_12 * lsh_421[k]
                   - f_10 * pc_y[k] * lsi1_563[k];

        t_732[k] = f_11 * lsh_422[k]
                   + f_3 * pc_y[k] * msh_548[k];

        t_733[k] = pa_y[k] * lsi0_565[k]
                   - f_10 * pc_y[k] * lsi1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pa_y, pc_y, pc_z, lsi0_566, lsi0_569, \
                         lsh_402, lsh_423, lsh_425, lsi1_566, lsi1_569, msh_549, \
                         msh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pa_y[k] * lsi0_566[k]
                   + f_13 * lsh_423[k]
                   - f_10 * pc_y[k] * lsi1_566[k];

        t_735[k] = f_20 * lsh_402[k]
                   + f_3 * pc_z[k] * msh_549[k];

        t_736[k] = f_11 * lsh_425[k]
                   + f_3 * pc_y[k] * msh_551[k];

        t_737[k] = pa_y[k] * lsi0_569[k]
                   - f_10 * pc_y[k] * lsi1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pa_y, pc_y, pc_z, lsi0_570, lsi0_572, lsh_405, \
                         lsh_426, lsh_428, lsi1_570, lsi1_572, \
                         msh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pa_y[k] * lsi0_570[k]
                   + f_14 * lsh_426[k]
                   - f_10 * pc_y[k] * lsi1_570[k];

        t_739[k] = f_20 * lsh_405[k]
                   + f_3 * pc_z[k] * msh_552[k];

        t_740[k] = pa_y[k] * lsi0_572[k]
                   + f_12 * lsh_428[k]
                   - f_10 * pc_y[k] * lsi1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pa_y, pc_x, pc_y, lsi0_574, lsh_429, \
                         lsh_561, lsh_562, lsi1_574, msh_555, msh_561, \
                         msh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * lsh_429[k]
                   + f_3 * pc_y[k] * msh_555[k];

        t_742[k] = pa_y[k] * lsi0_574[k]
                   - f_10 * pc_y[k] * lsi1_574[k];

        t_743[k] = f_13 * lsh_561[k]
                   + f_3 * pc_x[k] * msh_561[k];

        t_744[k] = f_13 * lsh_562[k]
                   + f_3 * pc_x[k] * msh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, lsh_563, lsh_564, lsh_565, lsh_566, \
                         msh_563, msh_564, msh_565, msh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_13 * lsh_563[k]
                   + f_3 * pc_x[k] * msh_563[k];

        t_746[k] = f_13 * lsh_564[k]
                   + f_3 * pc_x[k] * msh_564[k];

        t_747[k] = f_13 * lsh_565[k]
                   + f_3 * pc_x[k] * msh_565[k];

        t_748[k] = f_13 * lsh_566[k]
                   + f_3 * pc_x[k] * msh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, lsh_414, lsh_435, lsh_437, msg0_400, \
                         msg0_402, msg1_400, msg1_402, msh_561, \
                         msh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * lsh_435[k]
                   + f_1 * msg0_400[k]
                   - f_2 * msg1_400[k]
                   + f_3 * pc_y[k] * msh_561[k];

        t_750[k] = f_20 * lsh_414[k]
                   + f_3 * pc_z[k] * msh_561[k];

        t_751[k] = f_11 * lsh_437[k]
                   + f_8 * msg0_402[k]
                   - f_9 * msg1_402[k]
                   + f_3 * pc_y[k] * msh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, lsh_438, lsh_439, lsh_440, msg0_403, \
                         msg0_404, msg1_403, msg1_404, msh_564, msh_565, \
                         msh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * lsh_438[k]
                   + f_6 * msg0_403[k]
                   - f_7 * msg1_403[k]
                   + f_3 * pc_y[k] * msh_564[k];

        t_753[k] = f_11 * lsh_439[k]
                   + f_4 * msg0_404[k]
                   - f_5 * msg1_404[k]
                   + f_3 * pc_y[k] * msh_565[k];

        t_754[k] = f_11 * lsh_440[k]
                   + f_3 * pc_y[k] * msh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pa_y, pc_x, pc_y, pc_z, lsi0_587, \
                         lsh_420, lsh_567, lsi1_587, msg0_405, msg1_405, \
                         msh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pa_y[k] * lsi0_587[k]
                   - f_10 * pc_y[k] * lsi1_587[k];

        t_756[k] = f_13 * lsh_567[k]
                   + f_1 * msg0_405[k]
                   - f_2 * msg1_405[k]
                   + f_3 * pc_x[k] * msh_567[k];

        t_757[k] = f_3 * pc_y[k] * msh_567[k];

        t_758[k] = f_19 * lsh_420[k]
                   + f_3 * pc_z[k] * msh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, lsh_572, msg0_405, msg0_410, \
                         msg1_405, msg1_410, msh_568, msh_569, \
                         msh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_4 * msg0_405[k]
                   - f_5 * msg1_405[k]
                   + f_3 * pc_y[k] * msh_568[k];

        t_760[k] = f_3 * pc_y[k] * msh_569[k];

        t_761[k] = f_13 * lsh_572[k]
                   + f_8 * msg0_410[k]
                   - f_9 * msg1_410[k]
                   + f_3 * pc_x[k] * msh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_y, msg0_406, msg0_407, msg1_406, msg1_407, \
                         msh_570, msh_571, msh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_6 * msg0_406[k]
                   - f_7 * msg1_406[k]
                   + f_3 * pc_y[k] * msh_570[k];

        t_763[k] = f_4 * msg0_407[k]
                   - f_5 * msg1_407[k]
                   + f_3 * pc_y[k] * msh_571[k];

        t_764[k] = f_3 * pc_y[k] * msh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_y, lsh_576, msg0_408, msg0_409, \
                         msg0_414, msg1_408, msg1_409, msg1_414, msh_573, msh_574, \
                         msh_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_13 * lsh_576[k]
                   + f_6 * msg0_414[k]
                   - f_7 * msg1_414[k]
                   + f_3 * pc_x[k] * msh_576[k];

        t_766[k] = f_8 * msg0_408[k]
                   - f_9 * msg1_408[k]
                   + f_3 * pc_y[k] * msh_573[k];

        t_767[k] = f_6 * msg0_409[k]
                   - f_7 * msg1_409[k]
                   + f_3 * pc_y[k] * msh_574[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, lsh_581, lsh_582, msg0_410, \
                         msg0_419, msg1_410, msg1_419, msh_575, msh_576, msh_581, \
                         msh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * msg0_410[k]
                   - f_5 * msg1_410[k]
                   + f_3 * pc_y[k] * msh_575[k];

        t_769[k] = f_3 * pc_y[k] * msh_576[k];

        t_770[k] = f_13 * lsh_581[k]
                   + f_4 * msg0_419[k]
                   - f_5 * msg1_419[k]
                   + f_3 * pc_x[k] * msh_581[k];

        t_771[k] = f_13 * lsh_582[k]
                   + f_3 * pc_x[k] * msh_582[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pc_x, pc_y, lsh_583, lsh_584, \
                         lsh_585, lsh_587, msh_581, msh_583, msh_584, msh_585, \
                         msh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_13 * lsh_583[k]
                   + f_3 * pc_x[k] * msh_583[k];

        t_773[k] = f_13 * lsh_584[k]
                   + f_3 * pc_x[k] * msh_584[k];

        t_774[k] = f_13 * lsh_585[k]
                   + f_3 * pc_x[k] * msh_585[k];

        t_775[k] = f_3 * pc_y[k] * msh_581[k];

        t_776[k] = f_13 * lsh_587[k]
                   + f_3 * pc_x[k] * msh_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, msg0_415, msg0_416, msg0_417, msg1_415, \
                         msg1_416, msg1_417, msh_582, msh_583, \
                         msh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_1 * msg0_415[k]
                   - f_2 * msg1_415[k]
                   + f_3 * pc_y[k] * msh_582[k];

        t_778[k] = f_16 * msg0_416[k]
                   - f_17 * msg1_416[k]
                   + f_3 * pc_y[k] * msh_583[k];

        t_779[k] = f_8 * msg0_417[k]
                   - f_9 * msg1_417[k]
                   + f_3 * pc_y[k] * msh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, lsh_440, msg0_418, msg0_419, \
                         msg1_418, msg1_419, msh_585, msh_586, \
                         msh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * msg0_418[k]
                   - f_7 * msg1_418[k]
                   + f_3 * pc_y[k] * msh_585[k];

        t_781[k] = f_4 * msg0_419[k]
                   - f_5 * msg1_419[k]
                   + f_3 * pc_y[k] * msh_586[k];

        t_782[k] = f_3 * pc_y[k] * msh_587[k];

        t_783[k] = f_19 * lsh_440[k]
                   + f_1 * msg0_419[k]
                   - f_2 * msg1_419[k]
                   + f_3 * pc_z[k] * msh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pc_x, pc_y, pc_z, lsh_441, lsh_588, \
                         lsh_591, msg0_420, msg0_423, msg1_420, msg1_423, msh_588, \
                         msh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_12 * lsh_588[k]
                   + f_1 * msg0_420[k]
                   - f_2 * msg1_420[k]
                   + f_3 * pc_x[k] * msh_588[k];

        t_785[k] = f_18 * lsh_441[k]
                   + f_3 * pc_y[k] * msh_588[k];

        t_786[k] = f_3 * pc_z[k] * msh_588[k];

        t_787[k] = f_12 * lsh_591[k]
                   + f_8 * msg0_423[k]
                   - f_9 * msg1_423[k]
                   + f_3 * pc_x[k] * msh_591[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, t_791, pc_x, pc_z, lsh_594, msg0_420, msg0_426, \
                         msg1_420, msg1_426, msh_589, msh_590, msh_591, \
                         msh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_3 * pc_z[k] * msh_589[k];

        t_789[k] = f_4 * msg0_420[k]
                   - f_5 * msg1_420[k]
                   + f_3 * pc_z[k] * msh_590[k];

        t_790[k] = f_12 * lsh_594[k]
                   + f_6 * msg0_426[k]
                   - f_7 * msg1_426[k]
                   + f_3 * pc_x[k] * msh_594[k];

        t_791[k] = f_3 * pc_z[k] * msh_591[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pc_x, pc_y, pc_z, lsh_446, lsh_598, \
                         msg0_422, msg0_430, msg1_422, msg1_430, msh_593, msh_594, \
                         msh_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_18 * lsh_446[k]
                   + f_3 * pc_y[k] * msh_593[k];

        t_793[k] = f_6 * msg0_422[k]
                   - f_7 * msg1_422[k]
                   + f_3 * pc_z[k] * msh_593[k];

        t_794[k] = f_12 * lsh_598[k]
                   + f_4 * msg0_430[k]
                   - f_5 * msg1_430[k]
                   + f_3 * pc_x[k] * msh_598[k];

        t_795[k] = f_3 * pc_z[k] * msh_594[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pc_x, pc_y, pc_z, lsh_450, lsh_603, \
                         msg0_423, msg0_425, msg1_423, msg1_425, msh_595, msh_597, \
                         msh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_4 * msg0_423[k]
                   - f_5 * msg1_423[k]
                   + f_3 * pc_z[k] * msh_595[k];

        t_797[k] = f_18 * lsh_450[k]
                   + f_3 * pc_y[k] * msh_597[k];

        t_798[k] = f_8 * msg0_425[k]
                   - f_9 * msg1_425[k]
                   + f_3 * pc_z[k] * msh_597[k];

        t_799[k] = f_12 * lsh_603[k]
                   + f_3 * pc_x[k] * msh_603[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, pc_x, pc_z, lsh_605, lsh_606, \
                         lsh_607, lsh_608, msh_598, msh_605, msh_606, msh_607, \
                         msh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_3 * pc_z[k] * msh_598[k];

        t_801[k] = f_12 * lsh_605[k]
                   + f_3 * pc_x[k] * msh_605[k];

        t_802[k] = f_12 * lsh_606[k]
                   + f_3 * pc_x[k] * msh_606[k];

        t_803[k] = f_12 * lsh_607[k]
                   + f_3 * pc_x[k] * msh_607[k];

        t_804[k] = f_12 * lsh_608[k]
                   + f_3 * pc_x[k] * msh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, pc_y, pc_z, lsh_456, msg0_430, msg0_431, \
                         msg1_430, msg1_431, msh_603, msh_604, \
                         msh_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_18 * lsh_456[k]
                   + f_1 * msg0_430[k]
                   - f_2 * msg1_430[k]
                   + f_3 * pc_y[k] * msh_603[k];

        t_806[k] = f_3 * pc_z[k] * msh_603[k];

        t_807[k] = f_4 * msg0_430[k]
                   - f_5 * msg1_430[k]
                   + f_3 * pc_z[k] * msh_604[k];

        t_808[k] = f_6 * msg0_431[k]
                   - f_7 * msg1_431[k]
                   + f_3 * pc_z[k] * msh_605[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, pa_z, pc_y, pc_z, lsi0_588, lsh_461, \
                         lsi1_588, msg0_432, msg0_434, msg1_432, msg1_434, msh_606, \
                         msh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_8 * msg0_432[k]
                   - f_9 * msg1_432[k]
                   + f_3 * pc_z[k] * msh_606[k];

        t_810[k] = f_18 * lsh_461[k]
                   + f_3 * pc_y[k] * msh_608[k];

        t_811[k] = f_1 * msg0_434[k]
                   - f_2 * msg1_434[k]
                   + f_3 * pc_z[k] * msh_608[k];

        t_812[k] = pa_z[k] * lsi0_588[k]
                   - f_10 * pc_z[k] * lsi1_588[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, t_816, pa_z, pc_y, pc_z, lsi0_591, lsh_441, \
                         lsh_462, lsh_464, lsi1_591, msh_609, msh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_19 * lsh_462[k]
                   + f_3 * pc_y[k] * msh_609[k];

        t_814[k] = f_11 * lsh_441[k]
                   + f_3 * pc_z[k] * msh_609[k];

        t_815[k] = pa_z[k] * lsi0_591[k]
                   - f_10 * pc_z[k] * lsi1_591[k];

        t_816[k] = f_19 * lsh_464[k]
                   + f_3 * pc_y[k] * msh_611[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pa_z, pc_x, pc_z, lsi0_594, lsh_444, lsh_614, \
                         lsi1_594, msg0_440, msg1_440, msh_612, \
                         msh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_12 * lsh_614[k]
                   + f_8 * msg0_440[k]
                   - f_9 * msg1_440[k]
                   + f_3 * pc_x[k] * msh_614[k];

        t_818[k] = pa_z[k] * lsi0_594[k]
                   - f_10 * pc_z[k] * lsi1_594[k];

        t_819[k] = f_11 * lsh_444[k]
                   + f_3 * pc_z[k] * msh_612[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *lsi0_598 = buffer.data(lsi0 + 598);
    const auto *lsi0_600 = buffer.data(lsi0 + 600);
    const auto *lsi0_609 = buffer.data(lsi0 + 609);

    const auto *lsh_447 = buffer.data(lsh + 447);
    const auto *lsh_448 = buffer.data(lsh + 448);
    const auto *lsh_456 = buffer.data(lsh + 456);
    const auto *lsh_461 = buffer.data(lsh + 461);
    const auto *lsh_462 = buffer.data(lsh + 462);
    const auto *lsh_465 = buffer.data(lsh + 465);
    const auto *lsh_467 = buffer.data(lsh + 467);
    const auto *lsh_468 = buffer.data(lsh + 468);
    const auto *lsh_471 = buffer.data(lsh + 471);
    const auto *lsh_477 = buffer.data(lsh + 477);
    const auto *lsh_479 = buffer.data(lsh + 479);
    const auto *lsh_480 = buffer.data(lsh + 480);
    const auto *lsh_481 = buffer.data(lsh + 481);
    const auto *lsh_482 = buffer.data(lsh + 482);
    const auto *lsh_483 = buffer.data(lsh + 483);
    const auto *lsh_485 = buffer.data(lsh + 485);
    const auto *lsh_486 = buffer.data(lsh + 486);
    const auto *lsh_488 = buffer.data(lsh + 488);
    const auto *lsh_489 = buffer.data(lsh + 489);
    const auto *lsh_492 = buffer.data(lsh + 492);
    const auto *lsh_498 = buffer.data(lsh + 498);
    const auto *lsh_500 = buffer.data(lsh + 500);
    const auto *lsh_501 = buffer.data(lsh + 501);
    const auto *lsh_502 = buffer.data(lsh + 502);
    const auto *lsh_503 = buffer.data(lsh + 503);
    const auto *lsh_504 = buffer.data(lsh + 504);
    const auto *lsh_506 = buffer.data(lsh + 506);
    const auto *lsh_507 = buffer.data(lsh + 507);
    const auto *lsh_509 = buffer.data(lsh + 509);
    const auto *lsh_510 = buffer.data(lsh + 510);
    const auto *lsh_513 = buffer.data(lsh + 513);
    const auto *lsh_519 = buffer.data(lsh + 519);
    const auto *lsh_521 = buffer.data(lsh + 521);
    const auto *lsh_522 = buffer.data(lsh + 522);
    const auto *lsh_523 = buffer.data(lsh + 523);
    const auto *lsh_524 = buffer.data(lsh + 524);
    const auto *lsh_525 = buffer.data(lsh + 525);
    const auto *lsh_527 = buffer.data(lsh + 527);
    const auto *lsh_530 = buffer.data(lsh + 530);
    const auto *lsh_534 = buffer.data(lsh + 534);
    const auto *lsh_540 = buffer.data(lsh + 540);
    const auto *lsh_542 = buffer.data(lsh + 542);
    const auto *lsh_543 = buffer.data(lsh + 543);
    const auto *lsh_544 = buffer.data(lsh + 544);
    const auto *lsh_545 = buffer.data(lsh + 545);
    const auto *lsh_618 = buffer.data(lsh + 618);
    const auto *lsh_623 = buffer.data(lsh + 623);
    const auto *lsh_624 = buffer.data(lsh + 624);
    const auto *lsh_625 = buffer.data(lsh + 625);
    const auto *lsh_626 = buffer.data(lsh + 626);
    const auto *lsh_627 = buffer.data(lsh + 627);
    const auto *lsh_628 = buffer.data(lsh + 628);
    const auto *lsh_629 = buffer.data(lsh + 629);
    const auto *lsh_630 = buffer.data(lsh + 630);
    const auto *lsh_633 = buffer.data(lsh + 633);
    const auto *lsh_635 = buffer.data(lsh + 635);
    const auto *lsh_636 = buffer.data(lsh + 636);
    const auto *lsh_639 = buffer.data(lsh + 639);
    const auto *lsh_640 = buffer.data(lsh + 640);
    const auto *lsh_642 = buffer.data(lsh + 642);
    const auto *lsh_644 = buffer.data(lsh + 644);
    const auto *lsh_645 = buffer.data(lsh + 645);
    const auto *lsh_646 = buffer.data(lsh + 646);
    const auto *lsh_647 = buffer.data(lsh + 647);
    const auto *lsh_648 = buffer.data(lsh + 648);
    const auto *lsh_649 = buffer.data(lsh + 649);
    const auto *lsh_650 = buffer.data(lsh + 650);
    const auto *lsh_651 = buffer.data(lsh + 651);
    const auto *lsh_654 = buffer.data(lsh + 654);
    const auto *lsh_656 = buffer.data(lsh + 656);
    const auto *lsh_657 = buffer.data(lsh + 657);
    const auto *lsh_660 = buffer.data(lsh + 660);
    const auto *lsh_661 = buffer.data(lsh + 661);
    const auto *lsh_663 = buffer.data(lsh + 663);
    const auto *lsh_665 = buffer.data(lsh + 665);
    const auto *lsh_666 = buffer.data(lsh + 666);
    const auto *lsh_667 = buffer.data(lsh + 667);
    const auto *lsh_668 = buffer.data(lsh + 668);
    const auto *lsh_669 = buffer.data(lsh + 669);
    const auto *lsh_670 = buffer.data(lsh + 670);
    const auto *lsh_671 = buffer.data(lsh + 671);
    const auto *lsh_672 = buffer.data(lsh + 672);
    const auto *lsh_675 = buffer.data(lsh + 675);
    const auto *lsh_677 = buffer.data(lsh + 677);
    const auto *lsh_678 = buffer.data(lsh + 678);
    const auto *lsh_681 = buffer.data(lsh + 681);
    const auto *lsh_682 = buffer.data(lsh + 682);
    const auto *lsh_684 = buffer.data(lsh + 684);
    const auto *lsh_686 = buffer.data(lsh + 686);
    const auto *lsh_687 = buffer.data(lsh + 687);
    const auto *lsh_688 = buffer.data(lsh + 688);
    const auto *lsh_689 = buffer.data(lsh + 689);
    const auto *lsh_690 = buffer.data(lsh + 690);
    const auto *lsh_691 = buffer.data(lsh + 691);
    const auto *lsh_692 = buffer.data(lsh + 692);
    const auto *lsh_693 = buffer.data(lsh + 693);

    const auto *lsi1_598 = buffer.data(lsi1 + 598);
    const auto *lsi1_600 = buffer.data(lsi1 + 600);
    const auto *lsi1_609 = buffer.data(lsi1 + 609);

    const auto *msg0_444 = buffer.data(msg0 + 444);
    const auto *msg0_447 = buffer.data(msg0 + 447);
    const auto *msg0_448 = buffer.data(msg0 + 448);
    const auto *msg0_449 = buffer.data(msg0 + 449);
    const auto *msg0_450 = buffer.data(msg0 + 450);
    const auto *msg0_453 = buffer.data(msg0 + 453);
    const auto *msg0_455 = buffer.data(msg0 + 455);
    const auto *msg0_456 = buffer.data(msg0 + 456);
    const auto *msg0_459 = buffer.data(msg0 + 459);
    const auto *msg0_460 = buffer.data(msg0 + 460);
    const auto *msg0_462 = buffer.data(msg0 + 462);
    const auto *msg0_463 = buffer.data(msg0 + 463);
    const auto *msg0_464 = buffer.data(msg0 + 464);
    const auto *msg0_465 = buffer.data(msg0 + 465);
    const auto *msg0_468 = buffer.data(msg0 + 468);
    const auto *msg0_470 = buffer.data(msg0 + 470);
    const auto *msg0_471 = buffer.data(msg0 + 471);
    const auto *msg0_474 = buffer.data(msg0 + 474);
    const auto *msg0_475 = buffer.data(msg0 + 475);
    const auto *msg0_477 = buffer.data(msg0 + 477);
    const auto *msg0_478 = buffer.data(msg0 + 478);
    const auto *msg0_479 = buffer.data(msg0 + 479);
    const auto *msg0_480 = buffer.data(msg0 + 480);
    const auto *msg0_483 = buffer.data(msg0 + 483);
    const auto *msg0_485 = buffer.data(msg0 + 485);
    const auto *msg0_486 = buffer.data(msg0 + 486);
    const auto *msg0_489 = buffer.data(msg0 + 489);
    const auto *msg0_490 = buffer.data(msg0 + 490);
    const auto *msg0_492 = buffer.data(msg0 + 492);
    const auto *msg0_493 = buffer.data(msg0 + 493);
    const auto *msg0_494 = buffer.data(msg0 + 494);
    const auto *msg0_495 = buffer.data(msg0 + 495);

    const auto *msg1_444 = buffer.data(msg1 + 444);
    const auto *msg1_447 = buffer.data(msg1 + 447);
    const auto *msg1_448 = buffer.data(msg1 + 448);
    const auto *msg1_449 = buffer.data(msg1 + 449);
    const auto *msg1_450 = buffer.data(msg1 + 450);
    const auto *msg1_453 = buffer.data(msg1 + 453);
    const auto *msg1_455 = buffer.data(msg1 + 455);
    const auto *msg1_456 = buffer.data(msg1 + 456);
    const auto *msg1_459 = buffer.data(msg1 + 459);
    const auto *msg1_460 = buffer.data(msg1 + 460);
    const auto *msg1_462 = buffer.data(msg1 + 462);
    const auto *msg1_463 = buffer.data(msg1 + 463);
    const auto *msg1_464 = buffer.data(msg1 + 464);
    const auto *msg1_465 = buffer.data(msg1 + 465);
    const auto *msg1_468 = buffer.data(msg1 + 468);
    const auto *msg1_470 = buffer.data(msg1 + 470);
    const auto *msg1_471 = buffer.data(msg1 + 471);
    const auto *msg1_474 = buffer.data(msg1 + 474);
    const auto *msg1_475 = buffer.data(msg1 + 475);
    const auto *msg1_477 = buffer.data(msg1 + 477);
    const auto *msg1_478 = buffer.data(msg1 + 478);
    const auto *msg1_479 = buffer.data(msg1 + 479);
    const auto *msg1_480 = buffer.data(msg1 + 480);
    const auto *msg1_483 = buffer.data(msg1 + 483);
    const auto *msg1_485 = buffer.data(msg1 + 485);
    const auto *msg1_486 = buffer.data(msg1 + 486);
    const auto *msg1_489 = buffer.data(msg1 + 489);
    const auto *msg1_490 = buffer.data(msg1 + 490);
    const auto *msg1_492 = buffer.data(msg1 + 492);
    const auto *msg1_493 = buffer.data(msg1 + 493);
    const auto *msg1_494 = buffer.data(msg1 + 494);
    const auto *msg1_495 = buffer.data(msg1 + 495);

    const auto *msh_614 = buffer.data(msh + 614);
    const auto *msh_615 = buffer.data(msh + 615);
    const auto *msh_618 = buffer.data(msh + 618);
    const auto *msh_623 = buffer.data(msh + 623);
    const auto *msh_624 = buffer.data(msh + 624);
    const auto *msh_625 = buffer.data(msh + 625);
    const auto *msh_626 = buffer.data(msh + 626);
    const auto *msh_627 = buffer.data(msh + 627);
    const auto *msh_628 = buffer.data(msh + 628);
    const auto *msh_629 = buffer.data(msh + 629);
    const auto *msh_630 = buffer.data(msh + 630);
    const auto *msh_632 = buffer.data(msh + 632);
    const auto *msh_633 = buffer.data(msh + 633);
    const auto *msh_635 = buffer.data(msh + 635);
    const auto *msh_636 = buffer.data(msh + 636);
    const auto *msh_639 = buffer.data(msh + 639);
    const auto *msh_640 = buffer.data(msh + 640);
    const auto *msh_642 = buffer.data(msh + 642);
    const auto *msh_644 = buffer.data(msh + 644);
    const auto *msh_645 = buffer.data(msh + 645);
    const auto *msh_646 = buffer.data(msh + 646);
    const auto *msh_647 = buffer.data(msh + 647);
    const auto *msh_648 = buffer.data(msh + 648);
    const auto *msh_649 = buffer.data(msh + 649);
    const auto *msh_650 = buffer.data(msh + 650);
    const auto *msh_651 = buffer.data(msh + 651);
    const auto *msh_653 = buffer.data(msh + 653);
    const auto *msh_654 = buffer.data(msh + 654);
    const auto *msh_656 = buffer.data(msh + 656);
    const auto *msh_657 = buffer.data(msh + 657);
    const auto *msh_660 = buffer.data(msh + 660);
    const auto *msh_661 = buffer.data(msh + 661);
    const auto *msh_663 = buffer.data(msh + 663);
    const auto *msh_665 = buffer.data(msh + 665);
    const auto *msh_666 = buffer.data(msh + 666);
    const auto *msh_667 = buffer.data(msh + 667);
    const auto *msh_668 = buffer.data(msh + 668);
    const auto *msh_669 = buffer.data(msh + 669);
    const auto *msh_670 = buffer.data(msh + 670);
    const auto *msh_671 = buffer.data(msh + 671);
    const auto *msh_672 = buffer.data(msh + 672);
    const auto *msh_674 = buffer.data(msh + 674);
    const auto *msh_675 = buffer.data(msh + 675);
    const auto *msh_677 = buffer.data(msh + 677);
    const auto *msh_678 = buffer.data(msh + 678);
    const auto *msh_681 = buffer.data(msh + 681);
    const auto *msh_682 = buffer.data(msh + 682);
    const auto *msh_684 = buffer.data(msh + 684);
    const auto *msh_686 = buffer.data(msh + 686);
    const auto *msh_687 = buffer.data(msh + 687);
    const auto *msh_688 = buffer.data(msh + 688);
    const auto *msh_689 = buffer.data(msh + 689);
    const auto *msh_690 = buffer.data(msh + 690);
    const auto *msh_691 = buffer.data(msh + 691);
    const auto *msh_692 = buffer.data(msh + 692);
    const auto *msh_693 = buffer.data(msh + 693);

#pragma omp simd aligned(t_820, t_821, t_822, pa_z, pc_x, pc_y, pc_z, lsi0_598, lsh_467, \
                         lsh_618, lsi1_598, msg0_444, msg1_444, msh_614, \
                         msh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_19 * lsh_467[k]
                   + f_3 * pc_y[k] * msh_614[k];

        t_821[k] = f_12 * lsh_618[k]
                   + f_6 * msg0_444[k]
                   - f_7 * msg1_444[k]
                   + f_3 * pc_x[k] * msh_618[k];

        t_822[k] = pa_z[k] * lsi0_598[k]
                   - f_10 * pc_z[k] * lsi1_598[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pa_z, pc_y, pc_z, lsi0_600, lsh_447, lsh_448, \
                         lsh_471, lsi1_600, msh_615, msh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_11 * lsh_447[k]
                   + f_3 * pc_z[k] * msh_615[k];

        t_824[k] = pa_z[k] * lsi0_600[k]
                   + f_12 * lsh_448[k]
                   - f_10 * pc_z[k] * lsi1_600[k];

        t_825[k] = f_19 * lsh_471[k]
                   + f_3 * pc_y[k] * msh_618[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, lsh_623, lsh_624, lsh_625, lsh_626, \
                         msg0_449, msg1_449, msh_623, msh_624, msh_625, \
                         msh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_12 * lsh_623[k]
                   + f_4 * msg0_449[k]
                   - f_5 * msg1_449[k]
                   + f_3 * pc_x[k] * msh_623[k];

        t_827[k] = f_12 * lsh_624[k]
                   + f_3 * pc_x[k] * msh_624[k];

        t_828[k] = f_12 * lsh_625[k]
                   + f_3 * pc_x[k] * msh_625[k];

        t_829[k] = f_12 * lsh_626[k]
                   + f_3 * pc_x[k] * msh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pa_z, pc_x, pc_z, lsi0_609, lsh_627, \
                         lsh_628, lsh_629, lsi1_609, msh_627, msh_628, \
                         msh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_12 * lsh_627[k]
                   + f_3 * pc_x[k] * msh_627[k];

        t_831[k] = f_12 * lsh_628[k]
                   + f_3 * pc_x[k] * msh_628[k];

        t_832[k] = f_12 * lsh_629[k]
                   + f_3 * pc_x[k] * msh_629[k];

        t_833[k] = pa_z[k] * lsi0_609[k]
                   - f_10 * pc_z[k] * lsi1_609[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_y, pc_z, lsh_456, lsh_479, lsh_480, msg0_447, \
                         msg0_448, msg1_447, msg1_448, msh_624, msh_626, \
                         msh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_11 * lsh_456[k]
                   + f_3 * pc_z[k] * msh_624[k];

        t_835[k] = f_19 * lsh_479[k]
                   + f_8 * msg0_447[k]
                   - f_9 * msg1_447[k]
                   + f_3 * pc_y[k] * msh_626[k];

        t_836[k] = f_19 * lsh_480[k]
                   + f_6 * msg0_448[k]
                   - f_7 * msg1_448[k]
                   + f_3 * pc_y[k] * msh_627[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, lsh_461, lsh_481, lsh_482, msg0_449, \
                         msg1_449, msh_628, msh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_19 * lsh_481[k]
                   + f_4 * msg0_449[k]
                   - f_5 * msg1_449[k]
                   + f_3 * pc_y[k] * msh_628[k];

        t_838[k] = f_19 * lsh_482[k]
                   + f_3 * pc_y[k] * msh_629[k];

        t_839[k] = f_11 * lsh_461[k]
                   + f_1 * msg0_449[k]
                   - f_2 * msg1_449[k]
                   + f_3 * pc_z[k] * msh_629[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, pc_y, pc_z, lsh_462, lsh_483, lsh_630, \
                         msg0_450, msg1_450, msh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_12 * lsh_630[k]
                   + f_1 * msg0_450[k]
                   - f_2 * msg1_450[k]
                   + f_3 * pc_x[k] * msh_630[k];

        t_841[k] = f_20 * lsh_483[k]
                   + f_3 * pc_y[k] * msh_630[k];

        t_842[k] = f_12 * lsh_462[k]
                   + f_3 * pc_z[k] * msh_630[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, pc_y, lsh_485, lsh_633, lsh_635, msg0_453, \
                         msg0_455, msg1_453, msg1_455, msh_632, msh_633, \
                         msh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_12 * lsh_633[k]
                   + f_8 * msg0_453[k]
                   - f_9 * msg1_453[k]
                   + f_3 * pc_x[k] * msh_633[k];

        t_844[k] = f_20 * lsh_485[k]
                   + f_3 * pc_y[k] * msh_632[k];

        t_845[k] = f_12 * lsh_635[k]
                   + f_8 * msg0_455[k]
                   - f_9 * msg1_455[k]
                   + f_3 * pc_x[k] * msh_635[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pc_x, pc_y, pc_z, lsh_465, lsh_488, lsh_636, \
                         msg0_456, msg1_456, msh_633, msh_635, \
                         msh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_12 * lsh_636[k]
                   + f_6 * msg0_456[k]
                   - f_7 * msg1_456[k]
                   + f_3 * pc_x[k] * msh_636[k];

        t_847[k] = f_12 * lsh_465[k]
                   + f_3 * pc_z[k] * msh_633[k];

        t_848[k] = f_20 * lsh_488[k]
                   + f_3 * pc_y[k] * msh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, lsh_468, lsh_639, lsh_640, msg0_459, \
                         msg0_460, msg1_459, msg1_460, msh_636, msh_639, \
                         msh_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_12 * lsh_639[k]
                   + f_6 * msg0_459[k]
                   - f_7 * msg1_459[k]
                   + f_3 * pc_x[k] * msh_639[k];

        t_850[k] = f_12 * lsh_640[k]
                   + f_4 * msg0_460[k]
                   - f_5 * msg1_460[k]
                   + f_3 * pc_x[k] * msh_640[k];

        t_851[k] = f_12 * lsh_468[k]
                   + f_3 * pc_z[k] * msh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, pc_x, pc_y, lsh_492, lsh_642, lsh_644, msg0_462, \
                         msg0_464, msg1_462, msg1_464, msh_639, msh_642, \
                         msh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_12 * lsh_642[k]
                   + f_4 * msg0_462[k]
                   - f_5 * msg1_462[k]
                   + f_3 * pc_x[k] * msh_642[k];

        t_853[k] = f_20 * lsh_492[k]
                   + f_3 * pc_y[k] * msh_639[k];

        t_854[k] = f_12 * lsh_644[k]
                   + f_4 * msg0_464[k]
                   - f_5 * msg1_464[k]
                   + f_3 * pc_x[k] * msh_644[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, pc_x, lsh_645, lsh_646, lsh_647, \
                         lsh_648, lsh_649, msh_645, msh_646, msh_647, msh_648, \
                         msh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_12 * lsh_645[k]
                   + f_3 * pc_x[k] * msh_645[k];

        t_856[k] = f_12 * lsh_646[k]
                   + f_3 * pc_x[k] * msh_646[k];

        t_857[k] = f_12 * lsh_647[k]
                   + f_3 * pc_x[k] * msh_647[k];

        t_858[k] = f_12 * lsh_648[k]
                   + f_3 * pc_x[k] * msh_648[k];

        t_859[k] = f_12 * lsh_649[k]
                   + f_3 * pc_x[k] * msh_649[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pc_x, pc_y, pc_z, lsh_477, lsh_498, lsh_650, \
                         msg0_460, msg1_460, msh_645, msh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_12 * lsh_650[k]
                   + f_3 * pc_x[k] * msh_650[k];

        t_861[k] = f_20 * lsh_498[k]
                   + f_1 * msg0_460[k]
                   - f_2 * msg1_460[k]
                   + f_3 * pc_y[k] * msh_645[k];

        t_862[k] = f_12 * lsh_477[k]
                   + f_3 * pc_z[k] * msh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_y, lsh_500, lsh_501, lsh_502, msg0_462, \
                         msg0_463, msg0_464, msg1_462, msg1_463, msg1_464, msh_647, msh_648, \
                         msh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_20 * lsh_500[k]
                   + f_8 * msg0_462[k]
                   - f_9 * msg1_462[k]
                   + f_3 * pc_y[k] * msh_647[k];

        t_864[k] = f_20 * lsh_501[k]
                   + f_6 * msg0_463[k]
                   - f_7 * msg1_463[k]
                   + f_3 * pc_y[k] * msh_648[k];

        t_865[k] = f_20 * lsh_502[k]
                   + f_4 * msg0_464[k]
                   - f_5 * msg1_464[k]
                   + f_3 * pc_y[k] * msh_649[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, lsh_482, lsh_503, lsh_651, \
                         msg0_464, msg0_465, msg1_464, msg1_465, msh_650, \
                         msh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_20 * lsh_503[k]
                   + f_3 * pc_y[k] * msh_650[k];

        t_867[k] = f_12 * lsh_482[k]
                   + f_1 * msg0_464[k]
                   - f_2 * msg1_464[k]
                   + f_3 * pc_z[k] * msh_650[k];

        t_868[k] = f_12 * lsh_651[k]
                   + f_1 * msg0_465[k]
                   - f_2 * msg1_465[k]
                   + f_3 * pc_x[k] * msh_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, pc_z, lsh_483, lsh_504, \
                         lsh_506, lsh_654, msg0_468, msg1_468, msh_651, msh_653, \
                         msh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * lsh_504[k]
                   + f_3 * pc_y[k] * msh_651[k];

        t_870[k] = f_13 * lsh_483[k]
                   + f_3 * pc_z[k] * msh_651[k];

        t_871[k] = f_12 * lsh_654[k]
                   + f_8 * msg0_468[k]
                   - f_9 * msg1_468[k]
                   + f_3 * pc_x[k] * msh_654[k];

        t_872[k] = f_14 * lsh_506[k]
                   + f_3 * pc_y[k] * msh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_z, lsh_486, lsh_656, lsh_657, msg0_470, \
                         msg0_471, msg1_470, msg1_471, msh_654, msh_656, \
                         msh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_12 * lsh_656[k]
                   + f_8 * msg0_470[k]
                   - f_9 * msg1_470[k]
                   + f_3 * pc_x[k] * msh_656[k];

        t_874[k] = f_12 * lsh_657[k]
                   + f_6 * msg0_471[k]
                   - f_7 * msg1_471[k]
                   + f_3 * pc_x[k] * msh_657[k];

        t_875[k] = f_13 * lsh_486[k]
                   + f_3 * pc_z[k] * msh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_y, lsh_509, lsh_660, lsh_661, msg0_474, \
                         msg0_475, msg1_474, msg1_475, msh_656, msh_660, \
                         msh_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * lsh_509[k]
                   + f_3 * pc_y[k] * msh_656[k];

        t_877[k] = f_12 * lsh_660[k]
                   + f_6 * msg0_474[k]
                   - f_7 * msg1_474[k]
                   + f_3 * pc_x[k] * msh_660[k];

        t_878[k] = f_12 * lsh_661[k]
                   + f_4 * msg0_475[k]
                   - f_5 * msg1_475[k]
                   + f_3 * pc_x[k] * msh_661[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, lsh_489, lsh_513, lsh_663, \
                         msg0_477, msg1_477, msh_657, msh_660, \
                         msh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * lsh_489[k]
                   + f_3 * pc_z[k] * msh_657[k];

        t_880[k] = f_12 * lsh_663[k]
                   + f_4 * msg0_477[k]
                   - f_5 * msg1_477[k]
                   + f_3 * pc_x[k] * msh_663[k];

        t_881[k] = f_14 * lsh_513[k]
                   + f_3 * pc_y[k] * msh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pc_x, lsh_665, lsh_666, lsh_667, lsh_668, \
                         msg0_479, msg1_479, msh_665, msh_666, msh_667, \
                         msh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_12 * lsh_665[k]
                   + f_4 * msg0_479[k]
                   - f_5 * msg1_479[k]
                   + f_3 * pc_x[k] * msh_665[k];

        t_883[k] = f_12 * lsh_666[k]
                   + f_3 * pc_x[k] * msh_666[k];

        t_884[k] = f_12 * lsh_667[k]
                   + f_3 * pc_x[k] * msh_667[k];

        t_885[k] = f_12 * lsh_668[k]
                   + f_3 * pc_x[k] * msh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pc_x, pc_y, lsh_519, lsh_669, lsh_670, \
                         lsh_671, msg0_475, msg1_475, msh_666, msh_669, msh_670, \
                         msh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_12 * lsh_669[k]
                   + f_3 * pc_x[k] * msh_669[k];

        t_887[k] = f_12 * lsh_670[k]
                   + f_3 * pc_x[k] * msh_670[k];

        t_888[k] = f_12 * lsh_671[k]
                   + f_3 * pc_x[k] * msh_671[k];

        t_889[k] = f_14 * lsh_519[k]
                   + f_1 * msg0_475[k]
                   - f_2 * msg1_475[k]
                   + f_3 * pc_y[k] * msh_666[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pc_y, pc_z, lsh_498, lsh_521, lsh_522, msg0_477, \
                         msg0_478, msg1_477, msg1_478, msh_666, msh_668, \
                         msh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * lsh_498[k]
                   + f_3 * pc_z[k] * msh_666[k];

        t_891[k] = f_14 * lsh_521[k]
                   + f_8 * msg0_477[k]
                   - f_9 * msg1_477[k]
                   + f_3 * pc_y[k] * msh_668[k];

        t_892[k] = f_14 * lsh_522[k]
                   + f_6 * msg0_478[k]
                   - f_7 * msg1_478[k]
                   + f_3 * pc_y[k] * msh_669[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, pc_z, lsh_503, lsh_523, lsh_524, msg0_479, \
                         msg1_479, msh_670, msh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_14 * lsh_523[k]
                   + f_4 * msg0_479[k]
                   - f_5 * msg1_479[k]
                   + f_3 * pc_y[k] * msh_670[k];

        t_894[k] = f_14 * lsh_524[k]
                   + f_3 * pc_y[k] * msh_671[k];

        t_895[k] = f_13 * lsh_503[k]
                   + f_1 * msg0_479[k]
                   - f_2 * msg1_479[k]
                   + f_3 * pc_z[k] * msh_671[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_x, pc_y, pc_z, lsh_504, lsh_525, lsh_672, \
                         msg0_480, msg1_480, msh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_12 * lsh_672[k]
                   + f_1 * msg0_480[k]
                   - f_2 * msg1_480[k]
                   + f_3 * pc_x[k] * msh_672[k];

        t_897[k] = f_13 * lsh_525[k]
                   + f_3 * pc_y[k] * msh_672[k];

        t_898[k] = f_14 * lsh_504[k]
                   + f_3 * pc_z[k] * msh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, lsh_527, lsh_675, lsh_677, msg0_483, \
                         msg0_485, msg1_483, msg1_485, msh_674, msh_675, \
                         msh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_12 * lsh_675[k]
                   + f_8 * msg0_483[k]
                   - f_9 * msg1_483[k]
                   + f_3 * pc_x[k] * msh_675[k];

        t_900[k] = f_13 * lsh_527[k]
                   + f_3 * pc_y[k] * msh_674[k];

        t_901[k] = f_12 * lsh_677[k]
                   + f_8 * msg0_485[k]
                   - f_9 * msg1_485[k]
                   + f_3 * pc_x[k] * msh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, lsh_507, lsh_530, lsh_678, \
                         msg0_486, msg1_486, msh_675, msh_677, \
                         msh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_12 * lsh_678[k]
                   + f_6 * msg0_486[k]
                   - f_7 * msg1_486[k]
                   + f_3 * pc_x[k] * msh_678[k];

        t_903[k] = f_14 * lsh_507[k]
                   + f_3 * pc_z[k] * msh_675[k];

        t_904[k] = f_13 * lsh_530[k]
                   + f_3 * pc_y[k] * msh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, lsh_510, lsh_681, lsh_682, msg0_489, \
                         msg0_490, msg1_489, msg1_490, msh_678, msh_681, \
                         msh_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_12 * lsh_681[k]
                   + f_6 * msg0_489[k]
                   - f_7 * msg1_489[k]
                   + f_3 * pc_x[k] * msh_681[k];

        t_906[k] = f_12 * lsh_682[k]
                   + f_4 * msg0_490[k]
                   - f_5 * msg1_490[k]
                   + f_3 * pc_x[k] * msh_682[k];

        t_907[k] = f_14 * lsh_510[k]
                   + f_3 * pc_z[k] * msh_678[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, lsh_534, lsh_684, lsh_686, msg0_492, \
                         msg0_494, msg1_492, msg1_494, msh_681, msh_684, \
                         msh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_12 * lsh_684[k]
                   + f_4 * msg0_492[k]
                   - f_5 * msg1_492[k]
                   + f_3 * pc_x[k] * msh_684[k];

        t_909[k] = f_13 * lsh_534[k]
                   + f_3 * pc_y[k] * msh_681[k];

        t_910[k] = f_12 * lsh_686[k]
                   + f_4 * msg0_494[k]
                   - f_5 * msg1_494[k]
                   + f_3 * pc_x[k] * msh_686[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pc_x, lsh_687, lsh_688, lsh_689, \
                         lsh_690, lsh_691, msh_687, msh_688, msh_689, msh_690, \
                         msh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_12 * lsh_687[k]
                   + f_3 * pc_x[k] * msh_687[k];

        t_912[k] = f_12 * lsh_688[k]
                   + f_3 * pc_x[k] * msh_688[k];

        t_913[k] = f_12 * lsh_689[k]
                   + f_3 * pc_x[k] * msh_689[k];

        t_914[k] = f_12 * lsh_690[k]
                   + f_3 * pc_x[k] * msh_690[k];

        t_915[k] = f_12 * lsh_691[k]
                   + f_3 * pc_x[k] * msh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pc_x, pc_y, pc_z, lsh_519, lsh_540, lsh_692, \
                         msg0_490, msg1_490, msh_687, msh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_12 * lsh_692[k]
                   + f_3 * pc_x[k] * msh_692[k];

        t_917[k] = f_13 * lsh_540[k]
                   + f_1 * msg0_490[k]
                   - f_2 * msg1_490[k]
                   + f_3 * pc_y[k] * msh_687[k];

        t_918[k] = f_14 * lsh_519[k]
                   + f_3 * pc_z[k] * msh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pc_y, lsh_542, lsh_543, lsh_544, msg0_492, \
                         msg0_493, msg0_494, msg1_492, msg1_493, msg1_494, msh_689, msh_690, \
                         msh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_13 * lsh_542[k]
                   + f_8 * msg0_492[k]
                   - f_9 * msg1_492[k]
                   + f_3 * pc_y[k] * msh_689[k];

        t_920[k] = f_13 * lsh_543[k]
                   + f_6 * msg0_493[k]
                   - f_7 * msg1_493[k]
                   + f_3 * pc_y[k] * msh_690[k];

        t_921[k] = f_13 * lsh_544[k]
                   + f_4 * msg0_494[k]
                   - f_5 * msg1_494[k]
                   + f_3 * pc_y[k] * msh_691[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pc_x, pc_y, pc_z, lsh_524, lsh_545, lsh_693, \
                         msg0_494, msg0_495, msg1_494, msg1_495, msh_692, \
                         msh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * lsh_545[k]
                   + f_3 * pc_y[k] * msh_692[k];

        t_923[k] = f_14 * lsh_524[k]
                   + f_1 * msg0_494[k]
                   - f_2 * msg1_494[k]
                   + f_3 * pc_z[k] * msh_692[k];

        t_924[k] = f_12 * lsh_693[k]
                   + f_1 * msg0_495[k]
                   - f_2 * msg1_495[k]
                   + f_3 * pc_x[k] * msh_693[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msg0, const size_t msg1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsi0_756 = buffer.data(lsi0 + 756);
    const auto *lsi0_759 = buffer.data(lsi0 + 759);
    const auto *lsi0_761 = buffer.data(lsi0 + 761);
    const auto *lsi0_762 = buffer.data(lsi0 + 762);
    const auto *lsi0_765 = buffer.data(lsi0 + 765);
    const auto *lsi0_766 = buffer.data(lsi0 + 766);
    const auto *lsi0_768 = buffer.data(lsi0 + 768);
    const auto *lsi0_770 = buffer.data(lsi0 + 770);
    const auto *lsi0_783 = buffer.data(lsi0 + 783);
    const auto *lsi0_784 = buffer.data(lsi0 + 784);
    const auto *lsi0_787 = buffer.data(lsi0 + 787);
    const auto *lsi0_790 = buffer.data(lsi0 + 790);
    const auto *lsi0_1008 = buffer.data(lsi0 + 1008);
    const auto *lsi0_1011 = buffer.data(lsi0 + 1011);
    const auto *lsi0_1014 = buffer.data(lsi0 + 1014);
    const auto *lsi0_1018 = buffer.data(lsi0 + 1018);
    const auto *lsi0_1029 = buffer.data(lsi0 + 1029);
    const auto *lsi0_1031 = buffer.data(lsi0 + 1031);
    const auto *lsi0_1032 = buffer.data(lsi0 + 1032);
    const auto *lsi0_1033 = buffer.data(lsi0 + 1033);
    const auto *lsi0_1035 = buffer.data(lsi0 + 1035);
    const auto *lsi0_1041 = buffer.data(lsi0 + 1041);

    const auto *lsh_525 = buffer.data(lsh + 525);
    const auto *lsh_528 = buffer.data(lsh + 528);
    const auto *lsh_531 = buffer.data(lsh + 531);
    const auto *lsh_540 = buffer.data(lsh + 540);
    const auto *lsh_545 = buffer.data(lsh + 545);
    const auto *lsh_546 = buffer.data(lsh + 546);
    const auto *lsh_548 = buffer.data(lsh + 548);
    const auto *lsh_549 = buffer.data(lsh + 549);
    const auto *lsh_551 = buffer.data(lsh + 551);
    const auto *lsh_552 = buffer.data(lsh + 552);
    const auto *lsh_555 = buffer.data(lsh + 555);
    const auto *lsh_561 = buffer.data(lsh + 561);
    const auto *lsh_563 = buffer.data(lsh + 563);
    const auto *lsh_564 = buffer.data(lsh + 564);
    const auto *lsh_565 = buffer.data(lsh + 565);
    const auto *lsh_566 = buffer.data(lsh + 566);
    const auto *lsh_567 = buffer.data(lsh + 567);
    const auto *lsh_568 = buffer.data(lsh + 568);
    const auto *lsh_569 = buffer.data(lsh + 569);
    const auto *lsh_570 = buffer.data(lsh + 570);
    const auto *lsh_572 = buffer.data(lsh + 572);
    const auto *lsh_573 = buffer.data(lsh + 573);
    const auto *lsh_575 = buffer.data(lsh + 575);
    const auto *lsh_576 = buffer.data(lsh + 576);
    const auto *lsh_582 = buffer.data(lsh + 582);
    const auto *lsh_584 = buffer.data(lsh + 584);
    const auto *lsh_585 = buffer.data(lsh + 585);
    const auto *lsh_586 = buffer.data(lsh + 586);
    const auto *lsh_587 = buffer.data(lsh + 587);
    const auto *lsh_588 = buffer.data(lsh + 588);
    const auto *lsh_593 = buffer.data(lsh + 593);
    const auto *lsh_597 = buffer.data(lsh + 597);
    const auto *lsh_608 = buffer.data(lsh + 608);
    const auto *lsh_609 = buffer.data(lsh + 609);
    const auto *lsh_611 = buffer.data(lsh + 611);
    const auto *lsh_696 = buffer.data(lsh + 696);
    const auto *lsh_698 = buffer.data(lsh + 698);
    const auto *lsh_699 = buffer.data(lsh + 699);
    const auto *lsh_702 = buffer.data(lsh + 702);
    const auto *lsh_703 = buffer.data(lsh + 703);
    const auto *lsh_705 = buffer.data(lsh + 705);
    const auto *lsh_707 = buffer.data(lsh + 707);
    const auto *lsh_708 = buffer.data(lsh + 708);
    const auto *lsh_709 = buffer.data(lsh + 709);
    const auto *lsh_710 = buffer.data(lsh + 710);
    const auto *lsh_711 = buffer.data(lsh + 711);
    const auto *lsh_712 = buffer.data(lsh + 712);
    const auto *lsh_713 = buffer.data(lsh + 713);
    const auto *lsh_729 = buffer.data(lsh + 729);
    const auto *lsh_730 = buffer.data(lsh + 730);
    const auto *lsh_731 = buffer.data(lsh + 731);
    const auto *lsh_732 = buffer.data(lsh + 732);
    const auto *lsh_733 = buffer.data(lsh + 733);
    const auto *lsh_734 = buffer.data(lsh + 734);
    const auto *lsh_735 = buffer.data(lsh + 735);
    const auto *lsh_740 = buffer.data(lsh + 740);
    const auto *lsh_744 = buffer.data(lsh + 744);
    const auto *lsh_749 = buffer.data(lsh + 749);
    const auto *lsh_750 = buffer.data(lsh + 750);
    const auto *lsh_751 = buffer.data(lsh + 751);
    const auto *lsh_752 = buffer.data(lsh + 752);
    const auto *lsh_753 = buffer.data(lsh + 753);
    const auto *lsh_755 = buffer.data(lsh + 755);
    const auto *lsh_756 = buffer.data(lsh + 756);
    const auto *lsh_759 = buffer.data(lsh + 759);
    const auto *lsh_762 = buffer.data(lsh + 762);
    const auto *lsh_766 = buffer.data(lsh + 766);
    const auto *lsh_771 = buffer.data(lsh + 771);
    const auto *lsh_773 = buffer.data(lsh + 773);
    const auto *lsh_774 = buffer.data(lsh + 774);
    const auto *lsh_775 = buffer.data(lsh + 775);
    const auto *lsh_776 = buffer.data(lsh + 776);
    const auto *lsh_782 = buffer.data(lsh + 782);

    const auto *lsi1_756 = buffer.data(lsi1 + 756);
    const auto *lsi1_759 = buffer.data(lsi1 + 759);
    const auto *lsi1_761 = buffer.data(lsi1 + 761);
    const auto *lsi1_762 = buffer.data(lsi1 + 762);
    const auto *lsi1_765 = buffer.data(lsi1 + 765);
    const auto *lsi1_766 = buffer.data(lsi1 + 766);
    const auto *lsi1_768 = buffer.data(lsi1 + 768);
    const auto *lsi1_770 = buffer.data(lsi1 + 770);
    const auto *lsi1_783 = buffer.data(lsi1 + 783);
    const auto *lsi1_784 = buffer.data(lsi1 + 784);
    const auto *lsi1_787 = buffer.data(lsi1 + 787);
    const auto *lsi1_790 = buffer.data(lsi1 + 790);
    const auto *lsi1_1008 = buffer.data(lsi1 + 1008);
    const auto *lsi1_1011 = buffer.data(lsi1 + 1011);
    const auto *lsi1_1014 = buffer.data(lsi1 + 1014);
    const auto *lsi1_1018 = buffer.data(lsi1 + 1018);
    const auto *lsi1_1029 = buffer.data(lsi1 + 1029);
    const auto *lsi1_1031 = buffer.data(lsi1 + 1031);
    const auto *lsi1_1032 = buffer.data(lsi1 + 1032);
    const auto *lsi1_1033 = buffer.data(lsi1 + 1033);
    const auto *lsi1_1035 = buffer.data(lsi1 + 1035);
    const auto *lsi1_1041 = buffer.data(lsi1 + 1041);

    const auto *msg0_498 = buffer.data(msg0 + 498);
    const auto *msg0_500 = buffer.data(msg0 + 500);
    const auto *msg0_501 = buffer.data(msg0 + 501);
    const auto *msg0_504 = buffer.data(msg0 + 504);
    const auto *msg0_505 = buffer.data(msg0 + 505);
    const auto *msg0_507 = buffer.data(msg0 + 507);
    const auto *msg0_508 = buffer.data(msg0 + 508);
    const auto *msg0_509 = buffer.data(msg0 + 509);
    const auto *msg0_520 = buffer.data(msg0 + 520);
    const auto *msg0_522 = buffer.data(msg0 + 522);
    const auto *msg0_523 = buffer.data(msg0 + 523);
    const auto *msg0_524 = buffer.data(msg0 + 524);
    const auto *msg0_525 = buffer.data(msg0 + 525);
    const auto *msg0_526 = buffer.data(msg0 + 526);
    const auto *msg0_527 = buffer.data(msg0 + 527);
    const auto *msg0_528 = buffer.data(msg0 + 528);
    const auto *msg0_529 = buffer.data(msg0 + 529);
    const auto *msg0_530 = buffer.data(msg0 + 530);
    const auto *msg0_534 = buffer.data(msg0 + 534);
    const auto *msg0_535 = buffer.data(msg0 + 535);
    const auto *msg0_536 = buffer.data(msg0 + 536);
    const auto *msg0_537 = buffer.data(msg0 + 537);
    const auto *msg0_538 = buffer.data(msg0 + 538);
    const auto *msg0_539 = buffer.data(msg0 + 539);
    const auto *msg0_540 = buffer.data(msg0 + 540);
    const auto *msg0_542 = buffer.data(msg0 + 542);
    const auto *msg0_543 = buffer.data(msg0 + 543);
    const auto *msg0_545 = buffer.data(msg0 + 545);

    const auto *msg1_498 = buffer.data(msg1 + 498);
    const auto *msg1_500 = buffer.data(msg1 + 500);
    const auto *msg1_501 = buffer.data(msg1 + 501);
    const auto *msg1_504 = buffer.data(msg1 + 504);
    const auto *msg1_505 = buffer.data(msg1 + 505);
    const auto *msg1_507 = buffer.data(msg1 + 507);
    const auto *msg1_508 = buffer.data(msg1 + 508);
    const auto *msg1_509 = buffer.data(msg1 + 509);
    const auto *msg1_520 = buffer.data(msg1 + 520);
    const auto *msg1_522 = buffer.data(msg1 + 522);
    const auto *msg1_523 = buffer.data(msg1 + 523);
    const auto *msg1_524 = buffer.data(msg1 + 524);
    const auto *msg1_525 = buffer.data(msg1 + 525);
    const auto *msg1_526 = buffer.data(msg1 + 526);
    const auto *msg1_527 = buffer.data(msg1 + 527);
    const auto *msg1_528 = buffer.data(msg1 + 528);
    const auto *msg1_529 = buffer.data(msg1 + 529);
    const auto *msg1_530 = buffer.data(msg1 + 530);
    const auto *msg1_534 = buffer.data(msg1 + 534);
    const auto *msg1_535 = buffer.data(msg1 + 535);
    const auto *msg1_536 = buffer.data(msg1 + 536);
    const auto *msg1_537 = buffer.data(msg1 + 537);
    const auto *msg1_538 = buffer.data(msg1 + 538);
    const auto *msg1_539 = buffer.data(msg1 + 539);
    const auto *msg1_540 = buffer.data(msg1 + 540);
    const auto *msg1_542 = buffer.data(msg1 + 542);
    const auto *msg1_543 = buffer.data(msg1 + 543);
    const auto *msg1_545 = buffer.data(msg1 + 545);

    const auto *msh_693 = buffer.data(msh + 693);
    const auto *msh_695 = buffer.data(msh + 695);
    const auto *msh_696 = buffer.data(msh + 696);
    const auto *msh_698 = buffer.data(msh + 698);
    const auto *msh_699 = buffer.data(msh + 699);
    const auto *msh_702 = buffer.data(msh + 702);
    const auto *msh_703 = buffer.data(msh + 703);
    const auto *msh_705 = buffer.data(msh + 705);
    const auto *msh_707 = buffer.data(msh + 707);
    const auto *msh_708 = buffer.data(msh + 708);
    const auto *msh_709 = buffer.data(msh + 709);
    const auto *msh_710 = buffer.data(msh + 710);
    const auto *msh_711 = buffer.data(msh + 711);
    const auto *msh_712 = buffer.data(msh + 712);
    const auto *msh_713 = buffer.data(msh + 713);
    const auto *msh_714 = buffer.data(msh + 714);
    const auto *msh_716 = buffer.data(msh + 716);
    const auto *msh_717 = buffer.data(msh + 717);
    const auto *msh_719 = buffer.data(msh + 719);
    const auto *msh_720 = buffer.data(msh + 720);
    const auto *msh_723 = buffer.data(msh + 723);
    const auto *msh_729 = buffer.data(msh + 729);
    const auto *msh_730 = buffer.data(msh + 730);
    const auto *msh_731 = buffer.data(msh + 731);
    const auto *msh_732 = buffer.data(msh + 732);
    const auto *msh_733 = buffer.data(msh + 733);
    const auto *msh_734 = buffer.data(msh + 734);
    const auto *msh_735 = buffer.data(msh + 735);
    const auto *msh_736 = buffer.data(msh + 736);
    const auto *msh_737 = buffer.data(msh + 737);
    const auto *msh_738 = buffer.data(msh + 738);
    const auto *msh_739 = buffer.data(msh + 739);
    const auto *msh_740 = buffer.data(msh + 740);
    const auto *msh_741 = buffer.data(msh + 741);
    const auto *msh_742 = buffer.data(msh + 742);
    const auto *msh_743 = buffer.data(msh + 743);
    const auto *msh_744 = buffer.data(msh + 744);
    const auto *msh_749 = buffer.data(msh + 749);
    const auto *msh_750 = buffer.data(msh + 750);
    const auto *msh_751 = buffer.data(msh + 751);
    const auto *msh_752 = buffer.data(msh + 752);
    const auto *msh_753 = buffer.data(msh + 753);
    const auto *msh_754 = buffer.data(msh + 754);
    const auto *msh_755 = buffer.data(msh + 755);
    const auto *msh_756 = buffer.data(msh + 756);
    const auto *msh_757 = buffer.data(msh + 757);
    const auto *msh_758 = buffer.data(msh + 758);
    const auto *msh_759 = buffer.data(msh + 759);
    const auto *msh_761 = buffer.data(msh + 761);
    const auto *msh_762 = buffer.data(msh + 762);
    const auto *msh_763 = buffer.data(msh + 763);
    const auto *msh_765 = buffer.data(msh + 765);
    const auto *msh_766 = buffer.data(msh + 766);
    const auto *msh_771 = buffer.data(msh + 771);
    const auto *msh_773 = buffer.data(msh + 773);
    const auto *msh_774 = buffer.data(msh + 774);
    const auto *msh_775 = buffer.data(msh + 775);
    const auto *msh_776 = buffer.data(msh + 776);
    const auto *msh_777 = buffer.data(msh + 777);
    const auto *msh_779 = buffer.data(msh + 779);

#pragma omp simd aligned(t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, lsh_525, lsh_546, \
                         lsh_548, lsh_696, msg0_498, msg1_498, msh_693, msh_695, \
                         msh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_12 * lsh_546[k]
                   + f_3 * pc_y[k] * msh_693[k];

        t_926[k] = f_20 * lsh_525[k]
                   + f_3 * pc_z[k] * msh_693[k];

        t_927[k] = f_12 * lsh_696[k]
                   + f_8 * msg0_498[k]
                   - f_9 * msg1_498[k]
                   + f_3 * pc_x[k] * msh_696[k];

        t_928[k] = f_12 * lsh_548[k]
                   + f_3 * pc_y[k] * msh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pc_x, pc_z, lsh_528, lsh_698, lsh_699, msg0_500, \
                         msg0_501, msg1_500, msg1_501, msh_696, msh_698, \
                         msh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_12 * lsh_698[k]
                   + f_8 * msg0_500[k]
                   - f_9 * msg1_500[k]
                   + f_3 * pc_x[k] * msh_698[k];

        t_930[k] = f_12 * lsh_699[k]
                   + f_6 * msg0_501[k]
                   - f_7 * msg1_501[k]
                   + f_3 * pc_x[k] * msh_699[k];

        t_931[k] = f_20 * lsh_528[k]
                   + f_3 * pc_z[k] * msh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pc_x, pc_y, lsh_551, lsh_702, lsh_703, msg0_504, \
                         msg0_505, msg1_504, msg1_505, msh_698, msh_702, \
                         msh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * lsh_551[k]
                   + f_3 * pc_y[k] * msh_698[k];

        t_933[k] = f_12 * lsh_702[k]
                   + f_6 * msg0_504[k]
                   - f_7 * msg1_504[k]
                   + f_3 * pc_x[k] * msh_702[k];

        t_934[k] = f_12 * lsh_703[k]
                   + f_4 * msg0_505[k]
                   - f_5 * msg1_505[k]
                   + f_3 * pc_x[k] * msh_703[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pc_x, pc_y, pc_z, lsh_531, lsh_555, lsh_705, \
                         msg0_507, msg1_507, msh_699, msh_702, \
                         msh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_20 * lsh_531[k]
                   + f_3 * pc_z[k] * msh_699[k];

        t_936[k] = f_12 * lsh_705[k]
                   + f_4 * msg0_507[k]
                   - f_5 * msg1_507[k]
                   + f_3 * pc_x[k] * msh_705[k];

        t_937[k] = f_12 * lsh_555[k]
                   + f_3 * pc_y[k] * msh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pc_x, lsh_707, lsh_708, lsh_709, lsh_710, \
                         msg0_509, msg1_509, msh_707, msh_708, msh_709, \
                         msh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_12 * lsh_707[k]
                   + f_4 * msg0_509[k]
                   - f_5 * msg1_509[k]
                   + f_3 * pc_x[k] * msh_707[k];

        t_939[k] = f_12 * lsh_708[k]
                   + f_3 * pc_x[k] * msh_708[k];

        t_940[k] = f_12 * lsh_709[k]
                   + f_3 * pc_x[k] * msh_709[k];

        t_941[k] = f_12 * lsh_710[k]
                   + f_3 * pc_x[k] * msh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, lsh_561, lsh_711, lsh_712, \
                         lsh_713, msg0_505, msg1_505, msh_708, msh_711, msh_712, \
                         msh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_12 * lsh_711[k]
                   + f_3 * pc_x[k] * msh_711[k];

        t_943[k] = f_12 * lsh_712[k]
                   + f_3 * pc_x[k] * msh_712[k];

        t_944[k] = f_12 * lsh_713[k]
                   + f_3 * pc_x[k] * msh_713[k];

        t_945[k] = f_12 * lsh_561[k]
                   + f_1 * msg0_505[k]
                   - f_2 * msg1_505[k]
                   + f_3 * pc_y[k] * msh_708[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pc_y, pc_z, lsh_540, lsh_563, lsh_564, msg0_507, \
                         msg0_508, msg1_507, msg1_508, msh_708, msh_710, \
                         msh_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_20 * lsh_540[k]
                   + f_3 * pc_z[k] * msh_708[k];

        t_947[k] = f_12 * lsh_563[k]
                   + f_8 * msg0_507[k]
                   - f_9 * msg1_507[k]
                   + f_3 * pc_y[k] * msh_710[k];

        t_948[k] = f_12 * lsh_564[k]
                   + f_6 * msg0_508[k]
                   - f_7 * msg1_508[k]
                   + f_3 * pc_y[k] * msh_711[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_y, pc_y, pc_z, lsi0_756, lsh_545, \
                         lsh_565, lsh_566, lsi1_756, msg0_509, msg1_509, msh_712, \
                         msh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_12 * lsh_565[k]
                   + f_4 * msg0_509[k]
                   - f_5 * msg1_509[k]
                   + f_3 * pc_y[k] * msh_712[k];

        t_950[k] = f_12 * lsh_566[k]
                   + f_3 * pc_y[k] * msh_713[k];

        t_951[k] = f_20 * lsh_545[k]
                   + f_1 * msg0_509[k]
                   - f_2 * msg1_509[k]
                   + f_3 * pc_z[k] * msh_713[k];

        t_952[k] = pa_y[k] * lsi0_756[k]
                   - f_10 * pc_y[k] * lsi1_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pa_y, pc_y, pc_z, lsi0_759, lsh_546, \
                         lsh_567, lsh_568, lsh_569, lsi1_759, msh_714, \
                         msh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_11 * lsh_567[k]
                   + f_3 * pc_y[k] * msh_714[k];

        t_954[k] = f_19 * lsh_546[k]
                   + f_3 * pc_z[k] * msh_714[k];

        t_955[k] = pa_y[k] * lsi0_759[k]
                   + f_12 * lsh_568[k]
                   - f_10 * pc_y[k] * lsi1_759[k];

        t_956[k] = f_11 * lsh_569[k]
                   + f_3 * pc_y[k] * msh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pa_y, pc_y, pc_z, lsi0_761, lsi0_762, \
                         lsh_549, lsh_570, lsh_572, lsi1_761, lsi1_762, msh_717, \
                         msh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pa_y[k] * lsi0_761[k]
                   - f_10 * pc_y[k] * lsi1_761[k];

        t_958[k] = pa_y[k] * lsi0_762[k]
                   + f_13 * lsh_570[k]
                   - f_10 * pc_y[k] * lsi1_762[k];

        t_959[k] = f_19 * lsh_549[k]
                   + f_3 * pc_z[k] * msh_717[k];

        t_960[k] = f_11 * lsh_572[k]
                   + f_3 * pc_y[k] * msh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pa_y, pc_y, pc_z, lsi0_765, lsi0_766, lsh_552, \
                         lsh_573, lsi1_765, lsi1_766, msh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pa_y[k] * lsi0_765[k]
                   - f_10 * pc_y[k] * lsi1_765[k];

        t_962[k] = pa_y[k] * lsi0_766[k]
                   + f_14 * lsh_573[k]
                   - f_10 * pc_y[k] * lsi1_766[k];

        t_963[k] = f_19 * lsh_552[k]
                   + f_3 * pc_z[k] * msh_720[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pa_y, pc_x, pc_y, lsi0_768, lsi0_770, \
                         lsh_575, lsh_576, lsh_729, lsi1_768, lsi1_770, msh_723, \
                         msh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pa_y[k] * lsi0_768[k]
                   + f_12 * lsh_575[k]
                   - f_10 * pc_y[k] * lsi1_768[k];

        t_965[k] = f_11 * lsh_576[k]
                   + f_3 * pc_y[k] * msh_723[k];

        t_966[k] = pa_y[k] * lsi0_770[k]
                   - f_10 * pc_y[k] * lsi1_770[k];

        t_967[k] = f_12 * lsh_729[k]
                   + f_3 * pc_x[k] * msh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, lsh_730, lsh_731, lsh_732, \
                         lsh_733, lsh_734, msh_730, msh_731, msh_732, msh_733, \
                         msh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_12 * lsh_730[k]
                   + f_3 * pc_x[k] * msh_730[k];

        t_969[k] = f_12 * lsh_731[k]
                   + f_3 * pc_x[k] * msh_731[k];

        t_970[k] = f_12 * lsh_732[k]
                   + f_3 * pc_x[k] * msh_732[k];

        t_971[k] = f_12 * lsh_733[k]
                   + f_3 * pc_x[k] * msh_733[k];

        t_972[k] = f_12 * lsh_734[k]
                   + f_3 * pc_x[k] * msh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pc_y, pc_z, lsh_561, lsh_582, lsh_584, msg0_520, \
                         msg0_522, msg1_520, msg1_522, msh_729, \
                         msh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_11 * lsh_582[k]
                   + f_1 * msg0_520[k]
                   - f_2 * msg1_520[k]
                   + f_3 * pc_y[k] * msh_729[k];

        t_974[k] = f_19 * lsh_561[k]
                   + f_3 * pc_z[k] * msh_729[k];

        t_975[k] = f_11 * lsh_584[k]
                   + f_8 * msg0_522[k]
                   - f_9 * msg1_522[k]
                   + f_3 * pc_y[k] * msh_731[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_y, lsh_585, lsh_586, lsh_587, msg0_523, \
                         msg0_524, msg1_523, msg1_524, msh_732, msh_733, \
                         msh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_11 * lsh_585[k]
                   + f_6 * msg0_523[k]
                   - f_7 * msg1_523[k]
                   + f_3 * pc_y[k] * msh_732[k];

        t_977[k] = f_11 * lsh_586[k]
                   + f_4 * msg0_524[k]
                   - f_5 * msg1_524[k]
                   + f_3 * pc_y[k] * msh_733[k];

        t_978[k] = f_11 * lsh_587[k]
                   + f_3 * pc_y[k] * msh_734[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pa_y, pc_x, pc_y, pc_z, lsi0_783, \
                         lsh_567, lsh_735, lsi1_783, msg0_525, msg1_525, \
                         msh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pa_y[k] * lsi0_783[k]
                   - f_10 * pc_y[k] * lsi1_783[k];

        t_980[k] = f_12 * lsh_735[k]
                   + f_1 * msg0_525[k]
                   - f_2 * msg1_525[k]
                   + f_3 * pc_x[k] * msh_735[k];

        t_981[k] = f_3 * pc_y[k] * msh_735[k];

        t_982[k] = f_18 * lsh_567[k]
                   + f_3 * pc_z[k] * msh_735[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_x, pc_y, lsh_740, msg0_525, msg0_530, \
                         msg1_525, msg1_530, msh_736, msh_737, \
                         msh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_4 * msg0_525[k]
                   - f_5 * msg1_525[k]
                   + f_3 * pc_y[k] * msh_736[k];

        t_984[k] = f_3 * pc_y[k] * msh_737[k];

        t_985[k] = f_12 * lsh_740[k]
                   + f_8 * msg0_530[k]
                   - f_9 * msg1_530[k]
                   + f_3 * pc_x[k] * msh_740[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_y, msg0_526, msg0_527, msg1_526, msg1_527, \
                         msh_738, msh_739, msh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_6 * msg0_526[k]
                   - f_7 * msg1_526[k]
                   + f_3 * pc_y[k] * msh_738[k];

        t_987[k] = f_4 * msg0_527[k]
                   - f_5 * msg1_527[k]
                   + f_3 * pc_y[k] * msh_739[k];

        t_988[k] = f_3 * pc_y[k] * msh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_y, lsh_744, msg0_528, msg0_529, \
                         msg0_534, msg1_528, msg1_529, msg1_534, msh_741, msh_742, \
                         msh_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_12 * lsh_744[k]
                   + f_6 * msg0_534[k]
                   - f_7 * msg1_534[k]
                   + f_3 * pc_x[k] * msh_744[k];

        t_990[k] = f_8 * msg0_528[k]
                   - f_9 * msg1_528[k]
                   + f_3 * pc_y[k] * msh_741[k];

        t_991[k] = f_6 * msg0_529[k]
                   - f_7 * msg1_529[k]
                   + f_3 * pc_y[k] * msh_742[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, pc_x, pc_y, lsh_749, lsh_750, msg0_530, \
                         msg0_539, msg1_530, msg1_539, msh_743, msh_744, msh_749, \
                         msh_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_4 * msg0_530[k]
                   - f_5 * msg1_530[k]
                   + f_3 * pc_y[k] * msh_743[k];

        t_993[k] = f_3 * pc_y[k] * msh_744[k];

        t_994[k] = f_12 * lsh_749[k]
                   + f_4 * msg0_539[k]
                   - f_5 * msg1_539[k]
                   + f_3 * pc_x[k] * msh_749[k];

        t_995[k] = f_12 * lsh_750[k]
                   + f_3 * pc_x[k] * msh_750[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, pc_x, pc_y, lsh_751, lsh_752, \
                         lsh_753, lsh_755, msh_749, msh_751, msh_752, msh_753, \
                         msh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_12 * lsh_751[k]
                   + f_3 * pc_x[k] * msh_751[k];

        t_997[k] = f_12 * lsh_752[k]
                   + f_3 * pc_x[k] * msh_752[k];

        t_998[k] = f_12 * lsh_753[k]
                   + f_3 * pc_x[k] * msh_753[k];

        t_999[k] = f_3 * pc_y[k] * msh_749[k];

        t_1000[k] = f_12 * lsh_755[k]
                    + f_3 * pc_x[k] * msh_755[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, pc_y, msg0_535, msg0_536, msg0_537, msg1_535, \
                         msg1_536, msg1_537, msh_750, msh_751, \
                         msh_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_1 * msg0_535[k]
                    - f_2 * msg1_535[k]
                    + f_3 * pc_y[k] * msh_750[k];

        t_1002[k] = f_16 * msg0_536[k]
                    - f_17 * msg1_536[k]
                    + f_3 * pc_y[k] * msh_751[k];

        t_1003[k] = f_8 * msg0_537[k]
                    - f_9 * msg1_537[k]
                    + f_3 * pc_y[k] * msh_752[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, lsh_587, msg0_538, \
                         msg0_539, msg1_538, msg1_539, msh_753, msh_754, \
                         msh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_6 * msg0_538[k]
                    - f_7 * msg1_538[k]
                    + f_3 * pc_y[k] * msh_753[k];

        t_1005[k] = f_4 * msg0_539[k]
                    - f_5 * msg1_539[k]
                    + f_3 * pc_y[k] * msh_754[k];

        t_1006[k] = f_3 * pc_y[k] * msh_755[k];

        t_1007[k] = f_18 * lsh_587[k]
                    + f_1 * msg0_539[k]
                    - f_2 * msg1_539[k]
                    + f_3 * pc_z[k] * msh_755[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pa_x, pc_x, pc_y, pc_z, lsi0_1008, \
                         lsi0_1011, lsh_588, lsh_756, lsh_759, lsi1_1008, lsi1_1011, \
                         msh_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = pa_x[k] * lsi0_1008[k]
                    + f_19 * lsh_756[k]
                    - f_10 * pc_x[k] * lsi1_1008[k];

        t_1009[k] = f_15 * lsh_588[k]
                    + f_3 * pc_y[k] * msh_756[k];

        t_1010[k] = f_3 * pc_z[k] * msh_756[k];

        t_1011[k] = pa_x[k] * lsi0_1011[k]
                    + f_14 * lsh_759[k]
                    - f_10 * pc_x[k] * lsi1_1011[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, t_1015, pa_x, pc_x, pc_z, lsi0_1014, lsh_762, \
                         lsi1_1014, msg0_540, msg1_540, msh_757, msh_758, \
                         msh_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_3 * pc_z[k] * msh_757[k];

        t_1013[k] = f_4 * msg0_540[k]
                    - f_5 * msg1_540[k]
                    + f_3 * pc_z[k] * msh_758[k];

        t_1014[k] = pa_x[k] * lsi0_1014[k]
                    + f_13 * lsh_762[k]
                    - f_10 * pc_x[k] * lsi1_1014[k];

        t_1015[k] = f_3 * pc_z[k] * msh_759[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pa_x, pc_x, pc_y, pc_z, lsi0_1018, \
                         lsh_593, lsh_766, lsi1_1018, msg0_542, msg1_542, msh_761, \
                         msh_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_15 * lsh_593[k]
                    + f_3 * pc_y[k] * msh_761[k];

        t_1017[k] = f_6 * msg0_542[k]
                    - f_7 * msg1_542[k]
                    + f_3 * pc_z[k] * msh_761[k];

        t_1018[k] = pa_x[k] * lsi0_1018[k]
                    + f_12 * lsh_766[k]
                    - f_10 * pc_x[k] * lsi1_1018[k];

        t_1019[k] = f_3 * pc_z[k] * msh_762[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, pc_x, pc_y, pc_z, lsh_597, lsh_771, \
                         msg0_543, msg0_545, msg1_543, msg1_545, msh_763, msh_765, \
                         msh_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_4 * msg0_543[k]
                    - f_5 * msg1_543[k]
                    + f_3 * pc_z[k] * msh_763[k];

        t_1021[k] = f_15 * lsh_597[k]
                    + f_3 * pc_y[k] * msh_765[k];

        t_1022[k] = f_8 * msg0_545[k]
                    - f_9 * msg1_545[k]
                    + f_3 * pc_z[k] * msh_765[k];

        t_1023[k] = f_11 * lsh_771[k]
                    + f_3 * pc_x[k] * msh_771[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, t_1027, t_1028, pc_x, pc_z, lsh_773, lsh_774, \
                         lsh_775, lsh_776, msh_766, msh_773, msh_774, msh_775, \
                         msh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_3 * pc_z[k] * msh_766[k];

        t_1025[k] = f_11 * lsh_773[k]
                    + f_3 * pc_x[k] * msh_773[k];

        t_1026[k] = f_11 * lsh_774[k]
                    + f_3 * pc_x[k] * msh_774[k];

        t_1027[k] = f_11 * lsh_775[k]
                    + f_3 * pc_x[k] * msh_775[k];

        t_1028[k] = f_11 * lsh_776[k]
                    + f_3 * pc_x[k] * msh_776[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pa_x, pc_x, pc_z, lsi0_1029, \
                         lsi0_1031, lsi0_1032, lsi1_1029, lsi1_1031, lsi1_1032, \
                         msh_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = pa_x[k] * lsi0_1029[k]
                    - f_10 * pc_x[k] * lsi1_1029[k];

        t_1030[k] = f_3 * pc_z[k] * msh_771[k];

        t_1031[k] = pa_x[k] * lsi0_1031[k]
                    - f_10 * pc_x[k] * lsi1_1031[k];

        t_1032[k] = pa_x[k] * lsi0_1032[k]
                    - f_10 * pc_x[k] * lsi1_1032[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, pa_x, pc_x, pc_y, lsi0_1033, lsi0_1035, \
                         lsh_608, lsi1_1033, lsi1_1035, msh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = pa_x[k] * lsi0_1033[k]
                    - f_10 * pc_x[k] * lsi1_1033[k];

        t_1034[k] = f_15 * lsh_608[k]
                    + f_3 * pc_y[k] * msh_776[k];

        t_1035[k] = pa_x[k] * lsi0_1035[k]
                    - f_10 * pc_x[k] * lsi1_1035[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, t_1039, pa_z, pc_y, pc_z, lsi0_784, lsi0_787, \
                         lsh_588, lsh_609, lsi1_784, lsi1_787, \
                         msh_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = pa_z[k] * lsi0_784[k]
                    - f_10 * pc_z[k] * lsi1_784[k];

        t_1037[k] = f_18 * lsh_609[k]
                    + f_3 * pc_y[k] * msh_777[k];

        t_1038[k] = f_11 * lsh_588[k]
                    + f_3 * pc_z[k] * msh_777[k];

        t_1039[k] = pa_z[k] * lsi0_787[k]
                    - f_10 * pc_z[k] * lsi1_787[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pa_x, pa_z, pc_x, pc_y, pc_z, lsi0_790, \
                         lsi0_1041, lsh_611, lsh_782, lsi1_790, lsi1_1041, \
                         msh_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_18 * lsh_611[k]
                    + f_3 * pc_y[k] * msh_779[k];

        t_1041[k] = pa_x[k] * lsi0_1041[k]
                    + f_14 * lsh_782[k]
                    - f_10 * pc_x[k] * lsi1_1041[k];

        t_1042[k] = pa_z[k] * lsi0_790[k]
                    - f_10 * pc_z[k] * lsi1_790[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsi0,
                                                          const size_t lsh, const size_t lsi1,
                                                          const size_t msh, const size_t ncols,
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
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsi0_794 = buffer.data(lsi0 + 794);
    const auto *lsi0_1045 = buffer.data(lsi0 + 1045);
    const auto *lsi0_1048 = buffer.data(lsi0 + 1048);
    const auto *lsi0_1050 = buffer.data(lsi0 + 1050);
    const auto *lsi0_1057 = buffer.data(lsi0 + 1057);
    const auto *lsi0_1059 = buffer.data(lsi0 + 1059);
    const auto *lsi0_1060 = buffer.data(lsi0 + 1060);
    const auto *lsi0_1061 = buffer.data(lsi0 + 1061);
    const auto *lsi0_1063 = buffer.data(lsi0 + 1063);
    const auto *lsi0_1064 = buffer.data(lsi0 + 1064);
    const auto *lsi0_1067 = buffer.data(lsi0 + 1067);
    const auto *lsi0_1069 = buffer.data(lsi0 + 1069);
    const auto *lsi0_1070 = buffer.data(lsi0 + 1070);
    const auto *lsi0_1073 = buffer.data(lsi0 + 1073);
    const auto *lsi0_1074 = buffer.data(lsi0 + 1074);
    const auto *lsi0_1076 = buffer.data(lsi0 + 1076);
    const auto *lsi0_1078 = buffer.data(lsi0 + 1078);
    const auto *lsi0_1085 = buffer.data(lsi0 + 1085);
    const auto *lsi0_1087 = buffer.data(lsi0 + 1087);
    const auto *lsi0_1088 = buffer.data(lsi0 + 1088);
    const auto *lsi0_1089 = buffer.data(lsi0 + 1089);
    const auto *lsi0_1091 = buffer.data(lsi0 + 1091);
    const auto *lsi0_1092 = buffer.data(lsi0 + 1092);
    const auto *lsi0_1095 = buffer.data(lsi0 + 1095);
    const auto *lsi0_1097 = buffer.data(lsi0 + 1097);
    const auto *lsi0_1098 = buffer.data(lsi0 + 1098);
    const auto *lsi0_1101 = buffer.data(lsi0 + 1101);
    const auto *lsi0_1102 = buffer.data(lsi0 + 1102);
    const auto *lsi0_1104 = buffer.data(lsi0 + 1104);
    const auto *lsi0_1106 = buffer.data(lsi0 + 1106);
    const auto *lsi0_1113 = buffer.data(lsi0 + 1113);
    const auto *lsi0_1115 = buffer.data(lsi0 + 1115);
    const auto *lsi0_1116 = buffer.data(lsi0 + 1116);
    const auto *lsi0_1117 = buffer.data(lsi0 + 1117);
    const auto *lsi0_1119 = buffer.data(lsi0 + 1119);
    const auto *lsi0_1120 = buffer.data(lsi0 + 1120);
    const auto *lsi0_1123 = buffer.data(lsi0 + 1123);
    const auto *lsi0_1125 = buffer.data(lsi0 + 1125);
    const auto *lsi0_1126 = buffer.data(lsi0 + 1126);
    const auto *lsi0_1129 = buffer.data(lsi0 + 1129);
    const auto *lsi0_1130 = buffer.data(lsi0 + 1130);
    const auto *lsi0_1132 = buffer.data(lsi0 + 1132);
    const auto *lsi0_1134 = buffer.data(lsi0 + 1134);
    const auto *lsi0_1141 = buffer.data(lsi0 + 1141);
    const auto *lsi0_1143 = buffer.data(lsi0 + 1143);
    const auto *lsi0_1144 = buffer.data(lsi0 + 1144);
    const auto *lsi0_1145 = buffer.data(lsi0 + 1145);
    const auto *lsi0_1147 = buffer.data(lsi0 + 1147);
    const auto *lsi0_1148 = buffer.data(lsi0 + 1148);
    const auto *lsi0_1151 = buffer.data(lsi0 + 1151);
    const auto *lsi0_1153 = buffer.data(lsi0 + 1153);
    const auto *lsi0_1154 = buffer.data(lsi0 + 1154);
    const auto *lsi0_1157 = buffer.data(lsi0 + 1157);
    const auto *lsi0_1158 = buffer.data(lsi0 + 1158);
    const auto *lsi0_1160 = buffer.data(lsi0 + 1160);

    const auto *lsh_591 = buffer.data(lsh + 591);
    const auto *lsh_594 = buffer.data(lsh + 594);
    const auto *lsh_603 = buffer.data(lsh + 603);
    const auto *lsh_609 = buffer.data(lsh + 609);
    const auto *lsh_612 = buffer.data(lsh + 612);
    const auto *lsh_614 = buffer.data(lsh + 614);
    const auto *lsh_615 = buffer.data(lsh + 615);
    const auto *lsh_618 = buffer.data(lsh + 618);
    const auto *lsh_624 = buffer.data(lsh + 624);
    const auto *lsh_629 = buffer.data(lsh + 629);
    const auto *lsh_630 = buffer.data(lsh + 630);
    const auto *lsh_632 = buffer.data(lsh + 632);
    const auto *lsh_633 = buffer.data(lsh + 633);
    const auto *lsh_635 = buffer.data(lsh + 635);
    const auto *lsh_636 = buffer.data(lsh + 636);
    const auto *lsh_639 = buffer.data(lsh + 639);
    const auto *lsh_645 = buffer.data(lsh + 645);
    const auto *lsh_650 = buffer.data(lsh + 650);
    const auto *lsh_651 = buffer.data(lsh + 651);
    const auto *lsh_653 = buffer.data(lsh + 653);
    const auto *lsh_654 = buffer.data(lsh + 654);
    const auto *lsh_656 = buffer.data(lsh + 656);
    const auto *lsh_657 = buffer.data(lsh + 657);
    const auto *lsh_660 = buffer.data(lsh + 660);
    const auto *lsh_666 = buffer.data(lsh + 666);
    const auto *lsh_671 = buffer.data(lsh + 671);
    const auto *lsh_672 = buffer.data(lsh + 672);
    const auto *lsh_674 = buffer.data(lsh + 674);
    const auto *lsh_675 = buffer.data(lsh + 675);
    const auto *lsh_677 = buffer.data(lsh + 677);
    const auto *lsh_678 = buffer.data(lsh + 678);
    const auto *lsh_681 = buffer.data(lsh + 681);
    const auto *lsh_692 = buffer.data(lsh + 692);
    const auto *lsh_693 = buffer.data(lsh + 693);
    const auto *lsh_695 = buffer.data(lsh + 695);
    const auto *lsh_698 = buffer.data(lsh + 698);
    const auto *lsh_702 = buffer.data(lsh + 702);
    const auto *lsh_786 = buffer.data(lsh + 786);
    const auto *lsh_789 = buffer.data(lsh + 789);
    const auto *lsh_791 = buffer.data(lsh + 791);
    const auto *lsh_792 = buffer.data(lsh + 792);
    const auto *lsh_793 = buffer.data(lsh + 793);
    const auto *lsh_794 = buffer.data(lsh + 794);
    const auto *lsh_795 = buffer.data(lsh + 795);
    const auto *lsh_796 = buffer.data(lsh + 796);
    const auto *lsh_797 = buffer.data(lsh + 797);
    const auto *lsh_798 = buffer.data(lsh + 798);
    const auto *lsh_801 = buffer.data(lsh + 801);
    const auto *lsh_803 = buffer.data(lsh + 803);
    const auto *lsh_804 = buffer.data(lsh + 804);
    const auto *lsh_807 = buffer.data(lsh + 807);
    const auto *lsh_808 = buffer.data(lsh + 808);
    const auto *lsh_810 = buffer.data(lsh + 810);
    const auto *lsh_812 = buffer.data(lsh + 812);
    const auto *lsh_813 = buffer.data(lsh + 813);
    const auto *lsh_814 = buffer.data(lsh + 814);
    const auto *lsh_815 = buffer.data(lsh + 815);
    const auto *lsh_816 = buffer.data(lsh + 816);
    const auto *lsh_817 = buffer.data(lsh + 817);
    const auto *lsh_818 = buffer.data(lsh + 818);
    const auto *lsh_819 = buffer.data(lsh + 819);
    const auto *lsh_822 = buffer.data(lsh + 822);
    const auto *lsh_824 = buffer.data(lsh + 824);
    const auto *lsh_825 = buffer.data(lsh + 825);
    const auto *lsh_828 = buffer.data(lsh + 828);
    const auto *lsh_829 = buffer.data(lsh + 829);
    const auto *lsh_831 = buffer.data(lsh + 831);
    const auto *lsh_833 = buffer.data(lsh + 833);
    const auto *lsh_834 = buffer.data(lsh + 834);
    const auto *lsh_835 = buffer.data(lsh + 835);
    const auto *lsh_836 = buffer.data(lsh + 836);
    const auto *lsh_837 = buffer.data(lsh + 837);
    const auto *lsh_838 = buffer.data(lsh + 838);
    const auto *lsh_839 = buffer.data(lsh + 839);
    const auto *lsh_840 = buffer.data(lsh + 840);
    const auto *lsh_843 = buffer.data(lsh + 843);
    const auto *lsh_845 = buffer.data(lsh + 845);
    const auto *lsh_846 = buffer.data(lsh + 846);
    const auto *lsh_849 = buffer.data(lsh + 849);
    const auto *lsh_850 = buffer.data(lsh + 850);
    const auto *lsh_852 = buffer.data(lsh + 852);
    const auto *lsh_854 = buffer.data(lsh + 854);
    const auto *lsh_855 = buffer.data(lsh + 855);
    const auto *lsh_856 = buffer.data(lsh + 856);
    const auto *lsh_857 = buffer.data(lsh + 857);
    const auto *lsh_858 = buffer.data(lsh + 858);
    const auto *lsh_859 = buffer.data(lsh + 859);
    const auto *lsh_860 = buffer.data(lsh + 860);
    const auto *lsh_861 = buffer.data(lsh + 861);
    const auto *lsh_864 = buffer.data(lsh + 864);
    const auto *lsh_866 = buffer.data(lsh + 866);
    const auto *lsh_867 = buffer.data(lsh + 867);
    const auto *lsh_870 = buffer.data(lsh + 870);
    const auto *lsh_871 = buffer.data(lsh + 871);
    const auto *lsh_873 = buffer.data(lsh + 873);

    const auto *lsi1_794 = buffer.data(lsi1 + 794);
    const auto *lsi1_1045 = buffer.data(lsi1 + 1045);
    const auto *lsi1_1048 = buffer.data(lsi1 + 1048);
    const auto *lsi1_1050 = buffer.data(lsi1 + 1050);
    const auto *lsi1_1057 = buffer.data(lsi1 + 1057);
    const auto *lsi1_1059 = buffer.data(lsi1 + 1059);
    const auto *lsi1_1060 = buffer.data(lsi1 + 1060);
    const auto *lsi1_1061 = buffer.data(lsi1 + 1061);
    const auto *lsi1_1063 = buffer.data(lsi1 + 1063);
    const auto *lsi1_1064 = buffer.data(lsi1 + 1064);
    const auto *lsi1_1067 = buffer.data(lsi1 + 1067);
    const auto *lsi1_1069 = buffer.data(lsi1 + 1069);
    const auto *lsi1_1070 = buffer.data(lsi1 + 1070);
    const auto *lsi1_1073 = buffer.data(lsi1 + 1073);
    const auto *lsi1_1074 = buffer.data(lsi1 + 1074);
    const auto *lsi1_1076 = buffer.data(lsi1 + 1076);
    const auto *lsi1_1078 = buffer.data(lsi1 + 1078);
    const auto *lsi1_1085 = buffer.data(lsi1 + 1085);
    const auto *lsi1_1087 = buffer.data(lsi1 + 1087);
    const auto *lsi1_1088 = buffer.data(lsi1 + 1088);
    const auto *lsi1_1089 = buffer.data(lsi1 + 1089);
    const auto *lsi1_1091 = buffer.data(lsi1 + 1091);
    const auto *lsi1_1092 = buffer.data(lsi1 + 1092);
    const auto *lsi1_1095 = buffer.data(lsi1 + 1095);
    const auto *lsi1_1097 = buffer.data(lsi1 + 1097);
    const auto *lsi1_1098 = buffer.data(lsi1 + 1098);
    const auto *lsi1_1101 = buffer.data(lsi1 + 1101);
    const auto *lsi1_1102 = buffer.data(lsi1 + 1102);
    const auto *lsi1_1104 = buffer.data(lsi1 + 1104);
    const auto *lsi1_1106 = buffer.data(lsi1 + 1106);
    const auto *lsi1_1113 = buffer.data(lsi1 + 1113);
    const auto *lsi1_1115 = buffer.data(lsi1 + 1115);
    const auto *lsi1_1116 = buffer.data(lsi1 + 1116);
    const auto *lsi1_1117 = buffer.data(lsi1 + 1117);
    const auto *lsi1_1119 = buffer.data(lsi1 + 1119);
    const auto *lsi1_1120 = buffer.data(lsi1 + 1120);
    const auto *lsi1_1123 = buffer.data(lsi1 + 1123);
    const auto *lsi1_1125 = buffer.data(lsi1 + 1125);
    const auto *lsi1_1126 = buffer.data(lsi1 + 1126);
    const auto *lsi1_1129 = buffer.data(lsi1 + 1129);
    const auto *lsi1_1130 = buffer.data(lsi1 + 1130);
    const auto *lsi1_1132 = buffer.data(lsi1 + 1132);
    const auto *lsi1_1134 = buffer.data(lsi1 + 1134);
    const auto *lsi1_1141 = buffer.data(lsi1 + 1141);
    const auto *lsi1_1143 = buffer.data(lsi1 + 1143);
    const auto *lsi1_1144 = buffer.data(lsi1 + 1144);
    const auto *lsi1_1145 = buffer.data(lsi1 + 1145);
    const auto *lsi1_1147 = buffer.data(lsi1 + 1147);
    const auto *lsi1_1148 = buffer.data(lsi1 + 1148);
    const auto *lsi1_1151 = buffer.data(lsi1 + 1151);
    const auto *lsi1_1153 = buffer.data(lsi1 + 1153);
    const auto *lsi1_1154 = buffer.data(lsi1 + 1154);
    const auto *lsi1_1157 = buffer.data(lsi1 + 1157);
    const auto *lsi1_1158 = buffer.data(lsi1 + 1158);
    const auto *lsi1_1160 = buffer.data(lsi1 + 1160);

    const auto *msh_780 = buffer.data(msh + 780);
    const auto *msh_782 = buffer.data(msh + 782);
    const auto *msh_783 = buffer.data(msh + 783);
    const auto *msh_786 = buffer.data(msh + 786);
    const auto *msh_792 = buffer.data(msh + 792);
    const auto *msh_793 = buffer.data(msh + 793);
    const auto *msh_794 = buffer.data(msh + 794);
    const auto *msh_795 = buffer.data(msh + 795);
    const auto *msh_796 = buffer.data(msh + 796);
    const auto *msh_797 = buffer.data(msh + 797);
    const auto *msh_798 = buffer.data(msh + 798);
    const auto *msh_800 = buffer.data(msh + 800);
    const auto *msh_801 = buffer.data(msh + 801);
    const auto *msh_803 = buffer.data(msh + 803);
    const auto *msh_804 = buffer.data(msh + 804);
    const auto *msh_807 = buffer.data(msh + 807);
    const auto *msh_813 = buffer.data(msh + 813);
    const auto *msh_814 = buffer.data(msh + 814);
    const auto *msh_815 = buffer.data(msh + 815);
    const auto *msh_816 = buffer.data(msh + 816);
    const auto *msh_817 = buffer.data(msh + 817);
    const auto *msh_818 = buffer.data(msh + 818);
    const auto *msh_819 = buffer.data(msh + 819);
    const auto *msh_821 = buffer.data(msh + 821);
    const auto *msh_822 = buffer.data(msh + 822);
    const auto *msh_824 = buffer.data(msh + 824);
    const auto *msh_825 = buffer.data(msh + 825);
    const auto *msh_828 = buffer.data(msh + 828);
    const auto *msh_834 = buffer.data(msh + 834);
    const auto *msh_835 = buffer.data(msh + 835);
    const auto *msh_836 = buffer.data(msh + 836);
    const auto *msh_837 = buffer.data(msh + 837);
    const auto *msh_838 = buffer.data(msh + 838);
    const auto *msh_839 = buffer.data(msh + 839);
    const auto *msh_840 = buffer.data(msh + 840);
    const auto *msh_842 = buffer.data(msh + 842);
    const auto *msh_843 = buffer.data(msh + 843);
    const auto *msh_845 = buffer.data(msh + 845);
    const auto *msh_846 = buffer.data(msh + 846);
    const auto *msh_849 = buffer.data(msh + 849);
    const auto *msh_855 = buffer.data(msh + 855);
    const auto *msh_856 = buffer.data(msh + 856);
    const auto *msh_857 = buffer.data(msh + 857);
    const auto *msh_858 = buffer.data(msh + 858);
    const auto *msh_859 = buffer.data(msh + 859);
    const auto *msh_860 = buffer.data(msh + 860);
    const auto *msh_861 = buffer.data(msh + 861);
    const auto *msh_863 = buffer.data(msh + 863);
    const auto *msh_864 = buffer.data(msh + 864);
    const auto *msh_866 = buffer.data(msh + 866);
    const auto *msh_867 = buffer.data(msh + 867);
    const auto *msh_870 = buffer.data(msh + 870);

#pragma omp simd aligned(t_1043, t_1044, t_1045, pa_x, pc_x, pc_y, pc_z, lsi0_1045, lsh_591, \
                         lsh_614, lsh_786, lsi1_1045, msh_780, \
                         msh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_11 * lsh_591[k]
                    + f_3 * pc_z[k] * msh_780[k];

        t_1044[k] = f_18 * lsh_614[k]
                    + f_3 * pc_y[k] * msh_782[k];

        t_1045[k] = pa_x[k] * lsi0_1045[k]
                    + f_13 * lsh_786[k]
                    - f_10 * pc_x[k] * lsi1_1045[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pa_x, pa_z, pc_x, pc_z, lsi0_794, lsi0_1048, \
                         lsh_594, lsh_789, lsi1_794, lsi1_1048, \
                         msh_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pa_z[k] * lsi0_794[k]
                    - f_10 * pc_z[k] * lsi1_794[k];

        t_1047[k] = f_11 * lsh_594[k]
                    + f_3 * pc_z[k] * msh_783[k];

        t_1048[k] = pa_x[k] * lsi0_1048[k]
                    + f_12 * lsh_789[k]
                    - f_10 * pc_x[k] * lsi1_1048[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, t_1052, pa_x, pc_x, pc_y, lsi0_1050, lsh_618, \
                         lsh_791, lsh_792, lsh_793, lsi1_1050, msh_786, msh_792, \
                         msh_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_18 * lsh_618[k]
                    + f_3 * pc_y[k] * msh_786[k];

        t_1050[k] = pa_x[k] * lsi0_1050[k]
                    + f_12 * lsh_791[k]
                    - f_10 * pc_x[k] * lsi1_1050[k];

        t_1051[k] = f_11 * lsh_792[k]
                    + f_3 * pc_x[k] * msh_792[k];

        t_1052[k] = f_11 * lsh_793[k]
                    + f_3 * pc_x[k] * msh_793[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, t_1056, pc_x, lsh_794, lsh_795, lsh_796, \
                         lsh_797, msh_794, msh_795, msh_796, msh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = f_11 * lsh_794[k]
                    + f_3 * pc_x[k] * msh_794[k];

        t_1054[k] = f_11 * lsh_795[k]
                    + f_3 * pc_x[k] * msh_795[k];

        t_1055[k] = f_11 * lsh_796[k]
                    + f_3 * pc_x[k] * msh_796[k];

        t_1056[k] = f_11 * lsh_797[k]
                    + f_3 * pc_x[k] * msh_797[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, t_1060, pa_x, pc_x, pc_z, lsi0_1057, \
                         lsi0_1059, lsi0_1060, lsh_603, lsi1_1057, lsi1_1059, lsi1_1060, \
                         msh_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = pa_x[k] * lsi0_1057[k]
                    - f_10 * pc_x[k] * lsi1_1057[k];

        t_1058[k] = f_11 * lsh_603[k]
                    + f_3 * pc_z[k] * msh_792[k];

        t_1059[k] = pa_x[k] * lsi0_1059[k]
                    - f_10 * pc_x[k] * lsi1_1059[k];

        t_1060[k] = pa_x[k] * lsi0_1060[k]
                    - f_10 * pc_x[k] * lsi1_1060[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, pa_x, pc_x, pc_y, lsi0_1061, \
                         lsi0_1063, lsi0_1064, lsh_629, lsh_798, lsi1_1061, lsi1_1063, \
                         lsi1_1064, msh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = pa_x[k] * lsi0_1061[k]
                    - f_10 * pc_x[k] * lsi1_1061[k];

        t_1062[k] = f_18 * lsh_629[k]
                    + f_3 * pc_y[k] * msh_797[k];

        t_1063[k] = pa_x[k] * lsi0_1063[k]
                    - f_10 * pc_x[k] * lsi1_1063[k];

        t_1064[k] = pa_x[k] * lsi0_1064[k]
                    + f_19 * lsh_798[k]
                    - f_10 * pc_x[k] * lsi1_1064[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, pa_x, pc_x, pc_y, pc_z, lsi0_1067, \
                         lsh_609, lsh_630, lsh_632, lsh_801, lsi1_1067, msh_798, \
                         msh_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_19 * lsh_630[k]
                    + f_3 * pc_y[k] * msh_798[k];

        t_1066[k] = f_12 * lsh_609[k]
                    + f_3 * pc_z[k] * msh_798[k];

        t_1067[k] = pa_x[k] * lsi0_1067[k]
                    + f_14 * lsh_801[k]
                    - f_10 * pc_x[k] * lsi1_1067[k];

        t_1068[k] = f_19 * lsh_632[k]
                    + f_3 * pc_y[k] * msh_800[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, pa_x, pc_x, pc_z, lsi0_1069, lsi0_1070, \
                         lsh_612, lsh_803, lsh_804, lsi1_1069, lsi1_1070, \
                         msh_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = pa_x[k] * lsi0_1069[k]
                    + f_14 * lsh_803[k]
                    - f_10 * pc_x[k] * lsi1_1069[k];

        t_1070[k] = pa_x[k] * lsi0_1070[k]
                    + f_13 * lsh_804[k]
                    - f_10 * pc_x[k] * lsi1_1070[k];

        t_1071[k] = f_12 * lsh_612[k]
                    + f_3 * pc_z[k] * msh_801[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, pa_x, pc_x, pc_y, lsi0_1073, lsi0_1074, \
                         lsh_635, lsh_807, lsh_808, lsi1_1073, lsi1_1074, \
                         msh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = f_19 * lsh_635[k]
                    + f_3 * pc_y[k] * msh_803[k];

        t_1073[k] = pa_x[k] * lsi0_1073[k]
                    + f_13 * lsh_807[k]
                    - f_10 * pc_x[k] * lsi1_1073[k];

        t_1074[k] = pa_x[k] * lsi0_1074[k]
                    + f_12 * lsh_808[k]
                    - f_10 * pc_x[k] * lsi1_1074[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, pa_x, pc_x, pc_y, pc_z, lsi0_1076, lsh_615, \
                         lsh_639, lsh_810, lsi1_1076, msh_804, \
                         msh_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_12 * lsh_615[k]
                    + f_3 * pc_z[k] * msh_804[k];

        t_1076[k] = pa_x[k] * lsi0_1076[k]
                    + f_12 * lsh_810[k]
                    - f_10 * pc_x[k] * lsi1_1076[k];

        t_1077[k] = f_19 * lsh_639[k]
                    + f_3 * pc_y[k] * msh_807[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, t_1081, pa_x, pc_x, lsi0_1078, lsh_812, \
                         lsh_813, lsh_814, lsh_815, lsi1_1078, msh_813, msh_814, \
                         msh_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = pa_x[k] * lsi0_1078[k]
                    + f_12 * lsh_812[k]
                    - f_10 * pc_x[k] * lsi1_1078[k];

        t_1079[k] = f_11 * lsh_813[k]
                    + f_3 * pc_x[k] * msh_813[k];

        t_1080[k] = f_11 * lsh_814[k]
                    + f_3 * pc_x[k] * msh_814[k];

        t_1081[k] = f_11 * lsh_815[k]
                    + f_3 * pc_x[k] * msh_815[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, t_1085, pa_x, pc_x, lsi0_1085, lsh_816, \
                         lsh_817, lsh_818, lsi1_1085, msh_816, msh_817, \
                         msh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_11 * lsh_816[k]
                    + f_3 * pc_x[k] * msh_816[k];

        t_1083[k] = f_11 * lsh_817[k]
                    + f_3 * pc_x[k] * msh_817[k];

        t_1084[k] = f_11 * lsh_818[k]
                    + f_3 * pc_x[k] * msh_818[k];

        t_1085[k] = pa_x[k] * lsi0_1085[k]
                    - f_10 * pc_x[k] * lsi1_1085[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, t_1089, pa_x, pc_x, pc_z, lsi0_1087, \
                         lsi0_1088, lsi0_1089, lsh_624, lsi1_1087, lsi1_1088, lsi1_1089, \
                         msh_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_12 * lsh_624[k]
                    + f_3 * pc_z[k] * msh_813[k];

        t_1087[k] = pa_x[k] * lsi0_1087[k]
                    - f_10 * pc_x[k] * lsi1_1087[k];

        t_1088[k] = pa_x[k] * lsi0_1088[k]
                    - f_10 * pc_x[k] * lsi1_1088[k];

        t_1089[k] = pa_x[k] * lsi0_1089[k]
                    - f_10 * pc_x[k] * lsi1_1089[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, t_1093, pa_x, pc_x, pc_y, lsi0_1091, \
                         lsi0_1092, lsh_650, lsh_651, lsh_819, lsi1_1091, lsi1_1092, msh_818, \
                         msh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_19 * lsh_650[k]
                    + f_3 * pc_y[k] * msh_818[k];

        t_1091[k] = pa_x[k] * lsi0_1091[k]
                    - f_10 * pc_x[k] * lsi1_1091[k];

        t_1092[k] = pa_x[k] * lsi0_1092[k]
                    + f_19 * lsh_819[k]
                    - f_10 * pc_x[k] * lsi1_1092[k];

        t_1093[k] = f_20 * lsh_651[k]
                    + f_3 * pc_y[k] * msh_819[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pa_x, pc_x, pc_y, pc_z, lsi0_1095, lsh_630, \
                         lsh_653, lsh_822, lsi1_1095, msh_819, \
                         msh_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_13 * lsh_630[k]
                    + f_3 * pc_z[k] * msh_819[k];

        t_1095[k] = pa_x[k] * lsi0_1095[k]
                    + f_14 * lsh_822[k]
                    - f_10 * pc_x[k] * lsi1_1095[k];

        t_1096[k] = f_20 * lsh_653[k]
                    + f_3 * pc_y[k] * msh_821[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pa_x, pc_x, pc_z, lsi0_1097, lsi0_1098, \
                         lsh_633, lsh_824, lsh_825, lsi1_1097, lsi1_1098, \
                         msh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = pa_x[k] * lsi0_1097[k]
                    + f_14 * lsh_824[k]
                    - f_10 * pc_x[k] * lsi1_1097[k];

        t_1098[k] = pa_x[k] * lsi0_1098[k]
                    + f_13 * lsh_825[k]
                    - f_10 * pc_x[k] * lsi1_1098[k];

        t_1099[k] = f_13 * lsh_633[k]
                    + f_3 * pc_z[k] * msh_822[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pa_x, pc_x, pc_y, lsi0_1101, lsi0_1102, \
                         lsh_656, lsh_828, lsh_829, lsi1_1101, lsi1_1102, \
                         msh_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_20 * lsh_656[k]
                    + f_3 * pc_y[k] * msh_824[k];

        t_1101[k] = pa_x[k] * lsi0_1101[k]
                    + f_13 * lsh_828[k]
                    - f_10 * pc_x[k] * lsi1_1101[k];

        t_1102[k] = pa_x[k] * lsi0_1102[k]
                    + f_12 * lsh_829[k]
                    - f_10 * pc_x[k] * lsi1_1102[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pa_x, pc_x, pc_y, pc_z, lsi0_1104, lsh_636, \
                         lsh_660, lsh_831, lsi1_1104, msh_825, \
                         msh_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_13 * lsh_636[k]
                    + f_3 * pc_z[k] * msh_825[k];

        t_1104[k] = pa_x[k] * lsi0_1104[k]
                    + f_12 * lsh_831[k]
                    - f_10 * pc_x[k] * lsi1_1104[k];

        t_1105[k] = f_20 * lsh_660[k]
                    + f_3 * pc_y[k] * msh_828[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pa_x, pc_x, lsi0_1106, lsh_833, \
                         lsh_834, lsh_835, lsh_836, lsi1_1106, msh_834, msh_835, \
                         msh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = pa_x[k] * lsi0_1106[k]
                    + f_12 * lsh_833[k]
                    - f_10 * pc_x[k] * lsi1_1106[k];

        t_1107[k] = f_11 * lsh_834[k]
                    + f_3 * pc_x[k] * msh_834[k];

        t_1108[k] = f_11 * lsh_835[k]
                    + f_3 * pc_x[k] * msh_835[k];

        t_1109[k] = f_11 * lsh_836[k]
                    + f_3 * pc_x[k] * msh_836[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pa_x, pc_x, lsi0_1113, lsh_837, \
                         lsh_838, lsh_839, lsi1_1113, msh_837, msh_838, \
                         msh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_11 * lsh_837[k]
                    + f_3 * pc_x[k] * msh_837[k];

        t_1111[k] = f_11 * lsh_838[k]
                    + f_3 * pc_x[k] * msh_838[k];

        t_1112[k] = f_11 * lsh_839[k]
                    + f_3 * pc_x[k] * msh_839[k];

        t_1113[k] = pa_x[k] * lsi0_1113[k]
                    - f_10 * pc_x[k] * lsi1_1113[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, t_1117, pa_x, pc_x, pc_z, lsi0_1115, \
                         lsi0_1116, lsi0_1117, lsh_645, lsi1_1115, lsi1_1116, lsi1_1117, \
                         msh_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_13 * lsh_645[k]
                    + f_3 * pc_z[k] * msh_834[k];

        t_1115[k] = pa_x[k] * lsi0_1115[k]
                    - f_10 * pc_x[k] * lsi1_1115[k];

        t_1116[k] = pa_x[k] * lsi0_1116[k]
                    - f_10 * pc_x[k] * lsi1_1116[k];

        t_1117[k] = pa_x[k] * lsi0_1117[k]
                    - f_10 * pc_x[k] * lsi1_1117[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, t_1121, pa_x, pc_x, pc_y, lsi0_1119, \
                         lsi0_1120, lsh_671, lsh_672, lsh_840, lsi1_1119, lsi1_1120, msh_839, \
                         msh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_20 * lsh_671[k]
                    + f_3 * pc_y[k] * msh_839[k];

        t_1119[k] = pa_x[k] * lsi0_1119[k]
                    - f_10 * pc_x[k] * lsi1_1119[k];

        t_1120[k] = pa_x[k] * lsi0_1120[k]
                    + f_19 * lsh_840[k]
                    - f_10 * pc_x[k] * lsi1_1120[k];

        t_1121[k] = f_14 * lsh_672[k]
                    + f_3 * pc_y[k] * msh_840[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, pa_x, pc_x, pc_y, pc_z, lsi0_1123, lsh_651, \
                         lsh_674, lsh_843, lsi1_1123, msh_840, \
                         msh_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = f_14 * lsh_651[k]
                    + f_3 * pc_z[k] * msh_840[k];

        t_1123[k] = pa_x[k] * lsi0_1123[k]
                    + f_14 * lsh_843[k]
                    - f_10 * pc_x[k] * lsi1_1123[k];

        t_1124[k] = f_14 * lsh_674[k]
                    + f_3 * pc_y[k] * msh_842[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, pa_x, pc_x, pc_z, lsi0_1125, lsi0_1126, \
                         lsh_654, lsh_845, lsh_846, lsi1_1125, lsi1_1126, \
                         msh_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = pa_x[k] * lsi0_1125[k]
                    + f_14 * lsh_845[k]
                    - f_10 * pc_x[k] * lsi1_1125[k];

        t_1126[k] = pa_x[k] * lsi0_1126[k]
                    + f_13 * lsh_846[k]
                    - f_10 * pc_x[k] * lsi1_1126[k];

        t_1127[k] = f_14 * lsh_654[k]
                    + f_3 * pc_z[k] * msh_843[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pa_x, pc_x, pc_y, lsi0_1129, lsi0_1130, \
                         lsh_677, lsh_849, lsh_850, lsi1_1129, lsi1_1130, \
                         msh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_14 * lsh_677[k]
                    + f_3 * pc_y[k] * msh_845[k];

        t_1129[k] = pa_x[k] * lsi0_1129[k]
                    + f_13 * lsh_849[k]
                    - f_10 * pc_x[k] * lsi1_1129[k];

        t_1130[k] = pa_x[k] * lsi0_1130[k]
                    + f_12 * lsh_850[k]
                    - f_10 * pc_x[k] * lsi1_1130[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pa_x, pc_x, pc_y, pc_z, lsi0_1132, lsh_657, \
                         lsh_681, lsh_852, lsi1_1132, msh_846, \
                         msh_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_14 * lsh_657[k]
                    + f_3 * pc_z[k] * msh_846[k];

        t_1132[k] = pa_x[k] * lsi0_1132[k]
                    + f_12 * lsh_852[k]
                    - f_10 * pc_x[k] * lsi1_1132[k];

        t_1133[k] = f_14 * lsh_681[k]
                    + f_3 * pc_y[k] * msh_849[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, pa_x, pc_x, lsi0_1134, lsh_854, \
                         lsh_855, lsh_856, lsh_857, lsi1_1134, msh_855, msh_856, \
                         msh_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = pa_x[k] * lsi0_1134[k]
                    + f_12 * lsh_854[k]
                    - f_10 * pc_x[k] * lsi1_1134[k];

        t_1135[k] = f_11 * lsh_855[k]
                    + f_3 * pc_x[k] * msh_855[k];

        t_1136[k] = f_11 * lsh_856[k]
                    + f_3 * pc_x[k] * msh_856[k];

        t_1137[k] = f_11 * lsh_857[k]
                    + f_3 * pc_x[k] * msh_857[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, t_1141, pa_x, pc_x, lsi0_1141, lsh_858, \
                         lsh_859, lsh_860, lsi1_1141, msh_858, msh_859, \
                         msh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_11 * lsh_858[k]
                    + f_3 * pc_x[k] * msh_858[k];

        t_1139[k] = f_11 * lsh_859[k]
                    + f_3 * pc_x[k] * msh_859[k];

        t_1140[k] = f_11 * lsh_860[k]
                    + f_3 * pc_x[k] * msh_860[k];

        t_1141[k] = pa_x[k] * lsi0_1141[k]
                    - f_10 * pc_x[k] * lsi1_1141[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, t_1145, pa_x, pc_x, pc_z, lsi0_1143, \
                         lsi0_1144, lsi0_1145, lsh_666, lsi1_1143, lsi1_1144, lsi1_1145, \
                         msh_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_14 * lsh_666[k]
                    + f_3 * pc_z[k] * msh_855[k];

        t_1143[k] = pa_x[k] * lsi0_1143[k]
                    - f_10 * pc_x[k] * lsi1_1143[k];

        t_1144[k] = pa_x[k] * lsi0_1144[k]
                    - f_10 * pc_x[k] * lsi1_1144[k];

        t_1145[k] = pa_x[k] * lsi0_1145[k]
                    - f_10 * pc_x[k] * lsi1_1145[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, t_1149, pa_x, pc_x, pc_y, lsi0_1147, \
                         lsi0_1148, lsh_692, lsh_693, lsh_861, lsi1_1147, lsi1_1148, msh_860, \
                         msh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_14 * lsh_692[k]
                    + f_3 * pc_y[k] * msh_860[k];

        t_1147[k] = pa_x[k] * lsi0_1147[k]
                    - f_10 * pc_x[k] * lsi1_1147[k];

        t_1148[k] = pa_x[k] * lsi0_1148[k]
                    + f_19 * lsh_861[k]
                    - f_10 * pc_x[k] * lsi1_1148[k];

        t_1149[k] = f_13 * lsh_693[k]
                    + f_3 * pc_y[k] * msh_861[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, pa_x, pc_x, pc_y, pc_z, lsi0_1151, lsh_672, \
                         lsh_695, lsh_864, lsi1_1151, msh_861, \
                         msh_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_20 * lsh_672[k]
                    + f_3 * pc_z[k] * msh_861[k];

        t_1151[k] = pa_x[k] * lsi0_1151[k]
                    + f_14 * lsh_864[k]
                    - f_10 * pc_x[k] * lsi1_1151[k];

        t_1152[k] = f_13 * lsh_695[k]
                    + f_3 * pc_y[k] * msh_863[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, pa_x, pc_x, pc_z, lsi0_1153, lsi0_1154, \
                         lsh_675, lsh_866, lsh_867, lsi1_1153, lsi1_1154, \
                         msh_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = pa_x[k] * lsi0_1153[k]
                    + f_14 * lsh_866[k]
                    - f_10 * pc_x[k] * lsi1_1153[k];

        t_1154[k] = pa_x[k] * lsi0_1154[k]
                    + f_13 * lsh_867[k]
                    - f_10 * pc_x[k] * lsi1_1154[k];

        t_1155[k] = f_20 * lsh_675[k]
                    + f_3 * pc_z[k] * msh_864[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, pa_x, pc_x, pc_y, lsi0_1157, lsi0_1158, \
                         lsh_698, lsh_870, lsh_871, lsi1_1157, lsi1_1158, \
                         msh_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_13 * lsh_698[k]
                    + f_3 * pc_y[k] * msh_866[k];

        t_1157[k] = pa_x[k] * lsi0_1157[k]
                    + f_13 * lsh_870[k]
                    - f_10 * pc_x[k] * lsi1_1157[k];

        t_1158[k] = pa_x[k] * lsi0_1158[k]
                    + f_12 * lsh_871[k]
                    - f_10 * pc_x[k] * lsi1_1158[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, pa_x, pc_x, pc_y, pc_z, lsi0_1160, lsh_678, \
                         lsh_702, lsh_873, lsi1_1160, msh_867, \
                         msh_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_20 * lsh_678[k]
                    + f_3 * pc_z[k] * msh_867[k];

        t_1160[k] = pa_x[k] * lsi0_1160[k]
                    + f_12 * lsh_873[k]
                    - f_10 * pc_x[k] * lsi1_1160[k];

        t_1161[k] = f_13 * lsh_702[k]
                    + f_3 * pc_y[k] * msh_870[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsi0,
                                                           const size_t lsh, const size_t lsi1,
                                                           const size_t msg0, const size_t msg1,
                                                           const size_t msh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsi0_980 = buffer.data(lsi0 + 980);
    const auto *lsi0_985 = buffer.data(lsi0 + 985);
    const auto *lsi0_989 = buffer.data(lsi0 + 989);
    const auto *lsi0_994 = buffer.data(lsi0 + 994);
    const auto *lsi0_1162 = buffer.data(lsi0 + 1162);
    const auto *lsi0_1169 = buffer.data(lsi0 + 1169);
    const auto *lsi0_1171 = buffer.data(lsi0 + 1171);
    const auto *lsi0_1172 = buffer.data(lsi0 + 1172);
    const auto *lsi0_1173 = buffer.data(lsi0 + 1173);
    const auto *lsi0_1175 = buffer.data(lsi0 + 1175);
    const auto *lsi0_1176 = buffer.data(lsi0 + 1176);
    const auto *lsi0_1179 = buffer.data(lsi0 + 1179);
    const auto *lsi0_1181 = buffer.data(lsi0 + 1181);
    const auto *lsi0_1182 = buffer.data(lsi0 + 1182);
    const auto *lsi0_1185 = buffer.data(lsi0 + 1185);
    const auto *lsi0_1186 = buffer.data(lsi0 + 1186);
    const auto *lsi0_1188 = buffer.data(lsi0 + 1188);
    const auto *lsi0_1190 = buffer.data(lsi0 + 1190);
    const auto *lsi0_1197 = buffer.data(lsi0 + 1197);
    const auto *lsi0_1199 = buffer.data(lsi0 + 1199);
    const auto *lsi0_1200 = buffer.data(lsi0 + 1200);
    const auto *lsi0_1201 = buffer.data(lsi0 + 1201);
    const auto *lsi0_1203 = buffer.data(lsi0 + 1203);
    const auto *lsi0_1207 = buffer.data(lsi0 + 1207);
    const auto *lsi0_1210 = buffer.data(lsi0 + 1210);
    const auto *lsi0_1214 = buffer.data(lsi0 + 1214);
    const auto *lsi0_1216 = buffer.data(lsi0 + 1216);
    const auto *lsi0_1225 = buffer.data(lsi0 + 1225);
    const auto *lsi0_1227 = buffer.data(lsi0 + 1227);
    const auto *lsi0_1228 = buffer.data(lsi0 + 1228);
    const auto *lsi0_1229 = buffer.data(lsi0 + 1229);
    const auto *lsi0_1231 = buffer.data(lsi0 + 1231);
    const auto *lsi0_1232 = buffer.data(lsi0 + 1232);
    const auto *lsi0_1237 = buffer.data(lsi0 + 1237);
    const auto *lsi0_1241 = buffer.data(lsi0 + 1241);
    const auto *lsi0_1246 = buffer.data(lsi0 + 1246);
    const auto *lsi0_1253 = buffer.data(lsi0 + 1253);
    const auto *lsi0_1254 = buffer.data(lsi0 + 1254);
    const auto *lsi0_1255 = buffer.data(lsi0 + 1255);
    const auto *lsi0_1256 = buffer.data(lsi0 + 1256);
    const auto *lsi0_1257 = buffer.data(lsi0 + 1257);
    const auto *lsi0_1259 = buffer.data(lsi0 + 1259);

    const auto *lsh_687 = buffer.data(lsh + 687);
    const auto *lsh_693 = buffer.data(lsh + 693);
    const auto *lsh_696 = buffer.data(lsh + 696);
    const auto *lsh_699 = buffer.data(lsh + 699);
    const auto *lsh_708 = buffer.data(lsh + 708);
    const auto *lsh_713 = buffer.data(lsh + 713);
    const auto *lsh_714 = buffer.data(lsh + 714);
    const auto *lsh_716 = buffer.data(lsh + 716);
    const auto *lsh_717 = buffer.data(lsh + 717);
    const auto *lsh_719 = buffer.data(lsh + 719);
    const auto *lsh_720 = buffer.data(lsh + 720);
    const auto *lsh_723 = buffer.data(lsh + 723);
    const auto *lsh_729 = buffer.data(lsh + 729);
    const auto *lsh_734 = buffer.data(lsh + 734);
    const auto *lsh_735 = buffer.data(lsh + 735);
    const auto *lsh_737 = buffer.data(lsh + 737);
    const auto *lsh_740 = buffer.data(lsh + 740);
    const auto *lsh_744 = buffer.data(lsh + 744);
    const auto *lsh_755 = buffer.data(lsh + 755);
    const auto *lsh_771 = buffer.data(lsh + 771);
    const auto *lsh_776 = buffer.data(lsh + 776);
    const auto *lsh_875 = buffer.data(lsh + 875);
    const auto *lsh_876 = buffer.data(lsh + 876);
    const auto *lsh_877 = buffer.data(lsh + 877);
    const auto *lsh_878 = buffer.data(lsh + 878);
    const auto *lsh_879 = buffer.data(lsh + 879);
    const auto *lsh_880 = buffer.data(lsh + 880);
    const auto *lsh_881 = buffer.data(lsh + 881);
    const auto *lsh_882 = buffer.data(lsh + 882);
    const auto *lsh_885 = buffer.data(lsh + 885);
    const auto *lsh_887 = buffer.data(lsh + 887);
    const auto *lsh_888 = buffer.data(lsh + 888);
    const auto *lsh_891 = buffer.data(lsh + 891);
    const auto *lsh_892 = buffer.data(lsh + 892);
    const auto *lsh_894 = buffer.data(lsh + 894);
    const auto *lsh_896 = buffer.data(lsh + 896);
    const auto *lsh_897 = buffer.data(lsh + 897);
    const auto *lsh_898 = buffer.data(lsh + 898);
    const auto *lsh_899 = buffer.data(lsh + 899);
    const auto *lsh_900 = buffer.data(lsh + 900);
    const auto *lsh_901 = buffer.data(lsh + 901);
    const auto *lsh_902 = buffer.data(lsh + 902);
    const auto *lsh_906 = buffer.data(lsh + 906);
    const auto *lsh_909 = buffer.data(lsh + 909);
    const auto *lsh_913 = buffer.data(lsh + 913);
    const auto *lsh_915 = buffer.data(lsh + 915);
    const auto *lsh_918 = buffer.data(lsh + 918);
    const auto *lsh_919 = buffer.data(lsh + 919);
    const auto *lsh_920 = buffer.data(lsh + 920);
    const auto *lsh_921 = buffer.data(lsh + 921);
    const auto *lsh_922 = buffer.data(lsh + 922);
    const auto *lsh_923 = buffer.data(lsh + 923);
    const auto *lsh_924 = buffer.data(lsh + 924);
    const auto *lsh_929 = buffer.data(lsh + 929);
    const auto *lsh_933 = buffer.data(lsh + 933);
    const auto *lsh_938 = buffer.data(lsh + 938);
    const auto *lsh_939 = buffer.data(lsh + 939);
    const auto *lsh_940 = buffer.data(lsh + 940);
    const auto *lsh_941 = buffer.data(lsh + 941);
    const auto *lsh_942 = buffer.data(lsh + 942);
    const auto *lsh_944 = buffer.data(lsh + 944);

    const auto *lsi1_980 = buffer.data(lsi1 + 980);
    const auto *lsi1_985 = buffer.data(lsi1 + 985);
    const auto *lsi1_989 = buffer.data(lsi1 + 989);
    const auto *lsi1_994 = buffer.data(lsi1 + 994);
    const auto *lsi1_1162 = buffer.data(lsi1 + 1162);
    const auto *lsi1_1169 = buffer.data(lsi1 + 1169);
    const auto *lsi1_1171 = buffer.data(lsi1 + 1171);
    const auto *lsi1_1172 = buffer.data(lsi1 + 1172);
    const auto *lsi1_1173 = buffer.data(lsi1 + 1173);
    const auto *lsi1_1175 = buffer.data(lsi1 + 1175);
    const auto *lsi1_1176 = buffer.data(lsi1 + 1176);
    const auto *lsi1_1179 = buffer.data(lsi1 + 1179);
    const auto *lsi1_1181 = buffer.data(lsi1 + 1181);
    const auto *lsi1_1182 = buffer.data(lsi1 + 1182);
    const auto *lsi1_1185 = buffer.data(lsi1 + 1185);
    const auto *lsi1_1186 = buffer.data(lsi1 + 1186);
    const auto *lsi1_1188 = buffer.data(lsi1 + 1188);
    const auto *lsi1_1190 = buffer.data(lsi1 + 1190);
    const auto *lsi1_1197 = buffer.data(lsi1 + 1197);
    const auto *lsi1_1199 = buffer.data(lsi1 + 1199);
    const auto *lsi1_1200 = buffer.data(lsi1 + 1200);
    const auto *lsi1_1201 = buffer.data(lsi1 + 1201);
    const auto *lsi1_1203 = buffer.data(lsi1 + 1203);
    const auto *lsi1_1207 = buffer.data(lsi1 + 1207);
    const auto *lsi1_1210 = buffer.data(lsi1 + 1210);
    const auto *lsi1_1214 = buffer.data(lsi1 + 1214);
    const auto *lsi1_1216 = buffer.data(lsi1 + 1216);
    const auto *lsi1_1225 = buffer.data(lsi1 + 1225);
    const auto *lsi1_1227 = buffer.data(lsi1 + 1227);
    const auto *lsi1_1228 = buffer.data(lsi1 + 1228);
    const auto *lsi1_1229 = buffer.data(lsi1 + 1229);
    const auto *lsi1_1231 = buffer.data(lsi1 + 1231);
    const auto *lsi1_1232 = buffer.data(lsi1 + 1232);
    const auto *lsi1_1237 = buffer.data(lsi1 + 1237);
    const auto *lsi1_1241 = buffer.data(lsi1 + 1241);
    const auto *lsi1_1246 = buffer.data(lsi1 + 1246);
    const auto *lsi1_1253 = buffer.data(lsi1 + 1253);
    const auto *lsi1_1254 = buffer.data(lsi1 + 1254);
    const auto *lsi1_1255 = buffer.data(lsi1 + 1255);
    const auto *lsi1_1256 = buffer.data(lsi1 + 1256);
    const auto *lsi1_1257 = buffer.data(lsi1 + 1257);
    const auto *lsi1_1259 = buffer.data(lsi1 + 1259);

    const auto *msg0_660 = buffer.data(msg0 + 660);
    const auto *msg0_661 = buffer.data(msg0 + 661);
    const auto *msg0_662 = buffer.data(msg0 + 662);
    const auto *msg0_663 = buffer.data(msg0 + 663);
    const auto *msg0_664 = buffer.data(msg0 + 664);
    const auto *msg0_665 = buffer.data(msg0 + 665);
    const auto *msg0_675 = buffer.data(msg0 + 675);
    const auto *msg0_676 = buffer.data(msg0 + 676);
    const auto *msg0_678 = buffer.data(msg0 + 678);
    const auto *msg0_680 = buffer.data(msg0 + 680);
    const auto *msg0_681 = buffer.data(msg0 + 681);
    const auto *msg0_683 = buffer.data(msg0 + 683);
    const auto *msg0_684 = buffer.data(msg0 + 684);
    const auto *msg0_685 = buffer.data(msg0 + 685);
    const auto *msg0_686 = buffer.data(msg0 + 686);
    const auto *msg0_687 = buffer.data(msg0 + 687);
    const auto *msg0_688 = buffer.data(msg0 + 688);
    const auto *msg0_689 = buffer.data(msg0 + 689);

    const auto *msg1_660 = buffer.data(msg1 + 660);
    const auto *msg1_661 = buffer.data(msg1 + 661);
    const auto *msg1_662 = buffer.data(msg1 + 662);
    const auto *msg1_663 = buffer.data(msg1 + 663);
    const auto *msg1_664 = buffer.data(msg1 + 664);
    const auto *msg1_665 = buffer.data(msg1 + 665);
    const auto *msg1_675 = buffer.data(msg1 + 675);
    const auto *msg1_676 = buffer.data(msg1 + 676);
    const auto *msg1_678 = buffer.data(msg1 + 678);
    const auto *msg1_680 = buffer.data(msg1 + 680);
    const auto *msg1_681 = buffer.data(msg1 + 681);
    const auto *msg1_683 = buffer.data(msg1 + 683);
    const auto *msg1_684 = buffer.data(msg1 + 684);
    const auto *msg1_685 = buffer.data(msg1 + 685);
    const auto *msg1_686 = buffer.data(msg1 + 686);
    const auto *msg1_687 = buffer.data(msg1 + 687);
    const auto *msg1_688 = buffer.data(msg1 + 688);
    const auto *msg1_689 = buffer.data(msg1 + 689);

    const auto *msh_876 = buffer.data(msh + 876);
    const auto *msh_877 = buffer.data(msh + 877);
    const auto *msh_878 = buffer.data(msh + 878);
    const auto *msh_879 = buffer.data(msh + 879);
    const auto *msh_880 = buffer.data(msh + 880);
    const auto *msh_881 = buffer.data(msh + 881);
    const auto *msh_882 = buffer.data(msh + 882);
    const auto *msh_884 = buffer.data(msh + 884);
    const auto *msh_885 = buffer.data(msh + 885);
    const auto *msh_887 = buffer.data(msh + 887);
    const auto *msh_888 = buffer.data(msh + 888);
    const auto *msh_891 = buffer.data(msh + 891);
    const auto *msh_897 = buffer.data(msh + 897);
    const auto *msh_898 = buffer.data(msh + 898);
    const auto *msh_899 = buffer.data(msh + 899);
    const auto *msh_900 = buffer.data(msh + 900);
    const auto *msh_901 = buffer.data(msh + 901);
    const auto *msh_902 = buffer.data(msh + 902);
    const auto *msh_903 = buffer.data(msh + 903);
    const auto *msh_905 = buffer.data(msh + 905);
    const auto *msh_906 = buffer.data(msh + 906);
    const auto *msh_908 = buffer.data(msh + 908);
    const auto *msh_909 = buffer.data(msh + 909);
    const auto *msh_912 = buffer.data(msh + 912);
    const auto *msh_918 = buffer.data(msh + 918);
    const auto *msh_919 = buffer.data(msh + 919);
    const auto *msh_920 = buffer.data(msh + 920);
    const auto *msh_921 = buffer.data(msh + 921);
    const auto *msh_922 = buffer.data(msh + 922);
    const auto *msh_923 = buffer.data(msh + 923);
    const auto *msh_924 = buffer.data(msh + 924);
    const auto *msh_925 = buffer.data(msh + 925);
    const auto *msh_926 = buffer.data(msh + 926);
    const auto *msh_927 = buffer.data(msh + 927);
    const auto *msh_928 = buffer.data(msh + 928);
    const auto *msh_929 = buffer.data(msh + 929);
    const auto *msh_930 = buffer.data(msh + 930);
    const auto *msh_931 = buffer.data(msh + 931);
    const auto *msh_932 = buffer.data(msh + 932);
    const auto *msh_933 = buffer.data(msh + 933);
    const auto *msh_938 = buffer.data(msh + 938);
    const auto *msh_939 = buffer.data(msh + 939);
    const auto *msh_940 = buffer.data(msh + 940);
    const auto *msh_941 = buffer.data(msh + 941);
    const auto *msh_942 = buffer.data(msh + 942);
    const auto *msh_944 = buffer.data(msh + 944);
    const auto *msh_945 = buffer.data(msh + 945);
    const auto *msh_946 = buffer.data(msh + 946);
    const auto *msh_948 = buffer.data(msh + 948);
    const auto *msh_950 = buffer.data(msh + 950);
    const auto *msh_951 = buffer.data(msh + 951);
    const auto *msh_953 = buffer.data(msh + 953);
    const auto *msh_954 = buffer.data(msh + 954);
    const auto *msh_955 = buffer.data(msh + 955);
    const auto *msh_957 = buffer.data(msh + 957);
    const auto *msh_958 = buffer.data(msh + 958);
    const auto *msh_959 = buffer.data(msh + 959);
    const auto *msh_960 = buffer.data(msh + 960);
    const auto *msh_961 = buffer.data(msh + 961);
    const auto *msh_962 = buffer.data(msh + 962);
    const auto *msh_963 = buffer.data(msh + 963);
    const auto *msh_964 = buffer.data(msh + 964);
    const auto *msh_965 = buffer.data(msh + 965);

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, pa_x, pc_x, lsi0_1162, lsh_875, \
                         lsh_876, lsh_877, lsh_878, lsi1_1162, msh_876, msh_877, \
                         msh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = pa_x[k] * lsi0_1162[k]
                    + f_12 * lsh_875[k]
                    - f_10 * pc_x[k] * lsi1_1162[k];

        t_1163[k] = f_11 * lsh_876[k]
                    + f_3 * pc_x[k] * msh_876[k];

        t_1164[k] = f_11 * lsh_877[k]
                    + f_3 * pc_x[k] * msh_877[k];

        t_1165[k] = f_11 * lsh_878[k]
                    + f_3 * pc_x[k] * msh_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, pa_x, pc_x, lsi0_1169, lsh_879, \
                         lsh_880, lsh_881, lsi1_1169, msh_879, msh_880, \
                         msh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_11 * lsh_879[k]
                    + f_3 * pc_x[k] * msh_879[k];

        t_1167[k] = f_11 * lsh_880[k]
                    + f_3 * pc_x[k] * msh_880[k];

        t_1168[k] = f_11 * lsh_881[k]
                    + f_3 * pc_x[k] * msh_881[k];

        t_1169[k] = pa_x[k] * lsi0_1169[k]
                    - f_10 * pc_x[k] * lsi1_1169[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, t_1173, pa_x, pc_x, pc_z, lsi0_1171, \
                         lsi0_1172, lsi0_1173, lsh_687, lsi1_1171, lsi1_1172, lsi1_1173, \
                         msh_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_20 * lsh_687[k]
                    + f_3 * pc_z[k] * msh_876[k];

        t_1171[k] = pa_x[k] * lsi0_1171[k]
                    - f_10 * pc_x[k] * lsi1_1171[k];

        t_1172[k] = pa_x[k] * lsi0_1172[k]
                    - f_10 * pc_x[k] * lsi1_1172[k];

        t_1173[k] = pa_x[k] * lsi0_1173[k]
                    - f_10 * pc_x[k] * lsi1_1173[k];
    }

#pragma omp simd aligned(t_1174, t_1175, t_1176, t_1177, pa_x, pc_x, pc_y, lsi0_1175, \
                         lsi0_1176, lsh_713, lsh_714, lsh_882, lsi1_1175, lsi1_1176, msh_881, \
                         msh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1174[k] = f_13 * lsh_713[k]
                    + f_3 * pc_y[k] * msh_881[k];

        t_1175[k] = pa_x[k] * lsi0_1175[k]
                    - f_10 * pc_x[k] * lsi1_1175[k];

        t_1176[k] = pa_x[k] * lsi0_1176[k]
                    + f_19 * lsh_882[k]
                    - f_10 * pc_x[k] * lsi1_1176[k];

        t_1177[k] = f_12 * lsh_714[k]
                    + f_3 * pc_y[k] * msh_882[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, pa_x, pc_x, pc_y, pc_z, lsi0_1179, lsh_693, \
                         lsh_716, lsh_885, lsi1_1179, msh_882, \
                         msh_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_19 * lsh_693[k]
                    + f_3 * pc_z[k] * msh_882[k];

        t_1179[k] = pa_x[k] * lsi0_1179[k]
                    + f_14 * lsh_885[k]
                    - f_10 * pc_x[k] * lsi1_1179[k];

        t_1180[k] = f_12 * lsh_716[k]
                    + f_3 * pc_y[k] * msh_884[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, pa_x, pc_x, pc_z, lsi0_1181, lsi0_1182, \
                         lsh_696, lsh_887, lsh_888, lsi1_1181, lsi1_1182, \
                         msh_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = pa_x[k] * lsi0_1181[k]
                    + f_14 * lsh_887[k]
                    - f_10 * pc_x[k] * lsi1_1181[k];

        t_1182[k] = pa_x[k] * lsi0_1182[k]
                    + f_13 * lsh_888[k]
                    - f_10 * pc_x[k] * lsi1_1182[k];

        t_1183[k] = f_19 * lsh_696[k]
                    + f_3 * pc_z[k] * msh_885[k];
    }

#pragma omp simd aligned(t_1184, t_1185, t_1186, pa_x, pc_x, pc_y, lsi0_1185, lsi0_1186, \
                         lsh_719, lsh_891, lsh_892, lsi1_1185, lsi1_1186, \
                         msh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1184[k] = f_12 * lsh_719[k]
                    + f_3 * pc_y[k] * msh_887[k];

        t_1185[k] = pa_x[k] * lsi0_1185[k]
                    + f_13 * lsh_891[k]
                    - f_10 * pc_x[k] * lsi1_1185[k];

        t_1186[k] = pa_x[k] * lsi0_1186[k]
                    + f_12 * lsh_892[k]
                    - f_10 * pc_x[k] * lsi1_1186[k];
    }

#pragma omp simd aligned(t_1187, t_1188, t_1189, pa_x, pc_x, pc_y, pc_z, lsi0_1188, lsh_699, \
                         lsh_723, lsh_894, lsi1_1188, msh_888, \
                         msh_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1187[k] = f_19 * lsh_699[k]
                    + f_3 * pc_z[k] * msh_888[k];

        t_1188[k] = pa_x[k] * lsi0_1188[k]
                    + f_12 * lsh_894[k]
                    - f_10 * pc_x[k] * lsi1_1188[k];

        t_1189[k] = f_12 * lsh_723[k]
                    + f_3 * pc_y[k] * msh_891[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, t_1193, pa_x, pc_x, lsi0_1190, lsh_896, \
                         lsh_897, lsh_898, lsh_899, lsi1_1190, msh_897, msh_898, \
                         msh_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = pa_x[k] * lsi0_1190[k]
                    + f_12 * lsh_896[k]
                    - f_10 * pc_x[k] * lsi1_1190[k];

        t_1191[k] = f_11 * lsh_897[k]
                    + f_3 * pc_x[k] * msh_897[k];

        t_1192[k] = f_11 * lsh_898[k]
                    + f_3 * pc_x[k] * msh_898[k];

        t_1193[k] = f_11 * lsh_899[k]
                    + f_3 * pc_x[k] * msh_899[k];
    }

#pragma omp simd aligned(t_1194, t_1195, t_1196, t_1197, pa_x, pc_x, lsi0_1197, lsh_900, \
                         lsh_901, lsh_902, lsi1_1197, msh_900, msh_901, \
                         msh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1194[k] = f_11 * lsh_900[k]
                    + f_3 * pc_x[k] * msh_900[k];

        t_1195[k] = f_11 * lsh_901[k]
                    + f_3 * pc_x[k] * msh_901[k];

        t_1196[k] = f_11 * lsh_902[k]
                    + f_3 * pc_x[k] * msh_902[k];

        t_1197[k] = pa_x[k] * lsi0_1197[k]
                    - f_10 * pc_x[k] * lsi1_1197[k];
    }

#pragma omp simd aligned(t_1198, t_1199, t_1200, t_1201, pa_x, pc_x, pc_z, lsi0_1199, \
                         lsi0_1200, lsi0_1201, lsh_708, lsi1_1199, lsi1_1200, lsi1_1201, \
                         msh_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1198[k] = f_19 * lsh_708[k]
                    + f_3 * pc_z[k] * msh_897[k];

        t_1199[k] = pa_x[k] * lsi0_1199[k]
                    - f_10 * pc_x[k] * lsi1_1199[k];

        t_1200[k] = pa_x[k] * lsi0_1200[k]
                    - f_10 * pc_x[k] * lsi1_1200[k];

        t_1201[k] = pa_x[k] * lsi0_1201[k]
                    - f_10 * pc_x[k] * lsi1_1201[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, pa_x, pa_y, pc_x, pc_y, lsi0_980, \
                         lsi0_1203, lsh_734, lsh_735, lsi1_980, lsi1_1203, msh_902, \
                         msh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_12 * lsh_734[k]
                    + f_3 * pc_y[k] * msh_902[k];

        t_1203[k] = pa_x[k] * lsi0_1203[k]
                    - f_10 * pc_x[k] * lsi1_1203[k];

        t_1204[k] = pa_y[k] * lsi0_980[k]
                    - f_10 * pc_y[k] * lsi1_980[k];

        t_1205[k] = f_11 * lsh_735[k]
                    + f_3 * pc_y[k] * msh_903[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, pa_x, pc_x, pc_y, pc_z, lsi0_1207, lsh_714, \
                         lsh_737, lsh_906, lsi1_1207, msh_903, \
                         msh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_18 * lsh_714[k]
                    + f_3 * pc_z[k] * msh_903[k];

        t_1207[k] = pa_x[k] * lsi0_1207[k]
                    + f_14 * lsh_906[k]
                    - f_10 * pc_x[k] * lsi1_1207[k];

        t_1208[k] = f_11 * lsh_737[k]
                    + f_3 * pc_y[k] * msh_905[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, pa_x, pa_y, pc_x, pc_y, pc_z, lsi0_985, \
                         lsi0_1210, lsh_717, lsh_909, lsi1_985, lsi1_1210, \
                         msh_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = pa_y[k] * lsi0_985[k]
                    - f_10 * pc_y[k] * lsi1_985[k];

        t_1210[k] = pa_x[k] * lsi0_1210[k]
                    + f_13 * lsh_909[k]
                    - f_10 * pc_x[k] * lsi1_1210[k];

        t_1211[k] = f_18 * lsh_717[k]
                    + f_3 * pc_z[k] * msh_906[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, pa_x, pa_y, pc_x, pc_y, lsi0_989, lsi0_1214, \
                         lsh_740, lsh_913, lsi1_989, lsi1_1214, \
                         msh_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_11 * lsh_740[k]
                    + f_3 * pc_y[k] * msh_908[k];

        t_1213[k] = pa_y[k] * lsi0_989[k]
                    - f_10 * pc_y[k] * lsi1_989[k];

        t_1214[k] = pa_x[k] * lsi0_1214[k]
                    + f_12 * lsh_913[k]
                    - f_10 * pc_x[k] * lsi1_1214[k];
    }

#pragma omp simd aligned(t_1215, t_1216, t_1217, pa_x, pc_x, pc_y, pc_z, lsi0_1216, lsh_720, \
                         lsh_744, lsh_915, lsi1_1216, msh_909, \
                         msh_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1215[k] = f_18 * lsh_720[k]
                    + f_3 * pc_z[k] * msh_909[k];

        t_1216[k] = pa_x[k] * lsi0_1216[k]
                    + f_12 * lsh_915[k]
                    - f_10 * pc_x[k] * lsi1_1216[k];

        t_1217[k] = f_11 * lsh_744[k]
                    + f_3 * pc_y[k] * msh_912[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, t_1221, pa_y, pc_x, pc_y, lsi0_994, lsh_918, \
                         lsh_919, lsh_920, lsi1_994, msh_918, msh_919, \
                         msh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = pa_y[k] * lsi0_994[k]
                    - f_10 * pc_y[k] * lsi1_994[k];

        t_1219[k] = f_11 * lsh_918[k]
                    + f_3 * pc_x[k] * msh_918[k];

        t_1220[k] = f_11 * lsh_919[k]
                    + f_3 * pc_x[k] * msh_919[k];

        t_1221[k] = f_11 * lsh_920[k]
                    + f_3 * pc_x[k] * msh_920[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pa_x, pc_x, lsi0_1225, lsh_921, \
                         lsh_922, lsh_923, lsi1_1225, msh_921, msh_922, \
                         msh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_11 * lsh_921[k]
                    + f_3 * pc_x[k] * msh_921[k];

        t_1223[k] = f_11 * lsh_922[k]
                    + f_3 * pc_x[k] * msh_922[k];

        t_1224[k] = f_11 * lsh_923[k]
                    + f_3 * pc_x[k] * msh_923[k];

        t_1225[k] = pa_x[k] * lsi0_1225[k]
                    - f_10 * pc_x[k] * lsi1_1225[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, pa_x, pc_x, pc_z, lsi0_1227, \
                         lsi0_1228, lsi0_1229, lsh_729, lsi1_1227, lsi1_1228, lsi1_1229, \
                         msh_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_18 * lsh_729[k]
                    + f_3 * pc_z[k] * msh_918[k];

        t_1227[k] = pa_x[k] * lsi0_1227[k]
                    - f_10 * pc_x[k] * lsi1_1227[k];

        t_1228[k] = pa_x[k] * lsi0_1228[k]
                    - f_10 * pc_x[k] * lsi1_1228[k];

        t_1229[k] = pa_x[k] * lsi0_1229[k]
                    - f_10 * pc_x[k] * lsi1_1229[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pa_x, pc_x, pc_y, lsi0_1231, \
                         lsi0_1232, lsh_755, lsh_924, lsi1_1231, lsi1_1232, msh_923, \
                         msh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_11 * lsh_755[k]
                    + f_3 * pc_y[k] * msh_923[k];

        t_1231[k] = pa_x[k] * lsi0_1231[k]
                    - f_10 * pc_x[k] * lsi1_1231[k];

        t_1232[k] = pa_x[k] * lsi0_1232[k]
                    + f_19 * lsh_924[k]
                    - f_10 * pc_x[k] * lsi1_1232[k];

        t_1233[k] = f_3 * pc_y[k] * msh_924[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pc_y, pc_z, lsh_735, msg0_660, msg1_660, \
                         msh_924, msh_925, msh_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = f_15 * lsh_735[k]
                    + f_3 * pc_z[k] * msh_924[k];

        t_1235[k] = f_4 * msg0_660[k]
                    - f_5 * msg1_660[k]
                    + f_3 * pc_y[k] * msh_925[k];

        t_1236[k] = f_3 * pc_y[k] * msh_926[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, pa_x, pc_x, pc_y, lsi0_1237, lsh_929, \
                         lsi1_1237, msg0_661, msg0_662, msg1_661, msg1_662, msh_927, \
                         msh_928 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = pa_x[k] * lsi0_1237[k]
                    + f_14 * lsh_929[k]
                    - f_10 * pc_x[k] * lsi1_1237[k];

        t_1238[k] = f_6 * msg0_661[k]
                    - f_7 * msg1_661[k]
                    + f_3 * pc_y[k] * msh_927[k];

        t_1239[k] = f_4 * msg0_662[k]
                    - f_5 * msg1_662[k]
                    + f_3 * pc_y[k] * msh_928[k];
    }

#pragma omp simd aligned(t_1240, t_1241, t_1242, pa_x, pc_x, pc_y, lsi0_1241, lsh_933, \
                         lsi1_1241, msg0_663, msg1_663, msh_929, \
                         msh_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1240[k] = f_3 * pc_y[k] * msh_929[k];

        t_1241[k] = pa_x[k] * lsi0_1241[k]
                    + f_13 * lsh_933[k]
                    - f_10 * pc_x[k] * lsi1_1241[k];

        t_1242[k] = f_8 * msg0_663[k]
                    - f_9 * msg1_663[k]
                    + f_3 * pc_y[k] * msh_930[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, pc_y, msg0_664, msg0_665, msg1_664, msg1_665, \
                         msh_931, msh_932, msh_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_6 * msg0_664[k]
                    - f_7 * msg1_664[k]
                    + f_3 * pc_y[k] * msh_931[k];

        t_1244[k] = f_4 * msg0_665[k]
                    - f_5 * msg1_665[k]
                    + f_3 * pc_y[k] * msh_932[k];

        t_1245[k] = f_3 * pc_y[k] * msh_933[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, pa_x, pc_x, lsi0_1246, lsh_938, \
                         lsh_939, lsh_940, lsh_941, lsi1_1246, msh_939, msh_940, \
                         msh_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = pa_x[k] * lsi0_1246[k]
                    + f_12 * lsh_938[k]
                    - f_10 * pc_x[k] * lsi1_1246[k];

        t_1247[k] = f_11 * lsh_939[k]
                    + f_3 * pc_x[k] * msh_939[k];

        t_1248[k] = f_11 * lsh_940[k]
                    + f_3 * pc_x[k] * msh_940[k];

        t_1249[k] = f_11 * lsh_941[k]
                    + f_3 * pc_x[k] * msh_941[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pa_x, pc_x, pc_y, lsi0_1253, lsh_942, \
                         lsh_944, lsi1_1253, msh_938, msh_942, \
                         msh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_11 * lsh_942[k]
                    + f_3 * pc_x[k] * msh_942[k];

        t_1251[k] = f_3 * pc_y[k] * msh_938[k];

        t_1252[k] = f_11 * lsh_944[k]
                    + f_3 * pc_x[k] * msh_944[k];

        t_1253[k] = pa_x[k] * lsi0_1253[k]
                    - f_10 * pc_x[k] * lsi1_1253[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, t_1257, pa_x, pc_x, lsi0_1254, lsi0_1255, \
                         lsi0_1256, lsi0_1257, lsi1_1254, lsi1_1255, lsi1_1256, \
                         lsi1_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = pa_x[k] * lsi0_1254[k]
                    - f_10 * pc_x[k] * lsi1_1254[k];

        t_1255[k] = pa_x[k] * lsi0_1255[k]
                    - f_10 * pc_x[k] * lsi1_1255[k];

        t_1256[k] = pa_x[k] * lsi0_1256[k]
                    - f_10 * pc_x[k] * lsi1_1256[k];

        t_1257[k] = pa_x[k] * lsi0_1257[k]
                    - f_10 * pc_x[k] * lsi1_1257[k];
    }

#pragma omp simd aligned(t_1258, t_1259, t_1260, t_1261, pa_x, pc_x, pc_y, lsi0_1259, \
                         lsi1_1259, msg0_675, msg0_676, msg1_675, msg1_676, msh_944, msh_945, \
                         msh_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1258[k] = f_3 * pc_y[k] * msh_944[k];

        t_1259[k] = pa_x[k] * lsi0_1259[k]
                    - f_10 * pc_x[k] * lsi1_1259[k];

        t_1260[k] = f_1 * msg0_675[k]
                    - f_2 * msg1_675[k]
                    + f_3 * pc_x[k] * msh_945[k];

        t_1261[k] = f_16 * msg0_676[k]
                    - f_17 * msg1_676[k]
                    + f_3 * pc_x[k] * msh_946[k];
    }

#pragma omp simd aligned(t_1262, t_1263, t_1264, t_1265, pc_x, pc_z, msg0_678, msg0_680, \
                         msg1_678, msg1_680, msh_945, msh_946, msh_948, \
                         msh_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1262[k] = f_3 * pc_z[k] * msh_945[k];

        t_1263[k] = f_8 * msg0_678[k]
                    - f_9 * msg1_678[k]
                    + f_3 * pc_x[k] * msh_948[k];

        t_1264[k] = f_3 * pc_z[k] * msh_946[k];

        t_1265[k] = f_8 * msg0_680[k]
                    - f_9 * msg1_680[k]
                    + f_3 * pc_x[k] * msh_950[k];
    }

#pragma omp simd aligned(t_1266, t_1267, t_1268, t_1269, pc_x, pc_z, msg0_681, msg0_683, \
                         msg0_684, msg1_681, msg1_683, msg1_684, msh_948, msh_951, msh_953, \
                         msh_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1266[k] = f_6 * msg0_681[k]
                    - f_7 * msg1_681[k]
                    + f_3 * pc_x[k] * msh_951[k];

        t_1267[k] = f_3 * pc_z[k] * msh_948[k];

        t_1268[k] = f_6 * msg0_683[k]
                    - f_7 * msg1_683[k]
                    + f_3 * pc_x[k] * msh_953[k];

        t_1269[k] = f_6 * msg0_684[k]
                    - f_7 * msg1_684[k]
                    + f_3 * pc_x[k] * msh_954[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, t_1273, pc_x, pc_z, msg0_685, msg0_687, \
                         msg0_688, msg1_685, msg1_687, msg1_688, msh_951, msh_955, msh_957, \
                         msh_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = f_4 * msg0_685[k]
                    - f_5 * msg1_685[k]
                    + f_3 * pc_x[k] * msh_955[k];

        t_1271[k] = f_3 * pc_z[k] * msh_951[k];

        t_1272[k] = f_4 * msg0_687[k]
                    - f_5 * msg1_687[k]
                    + f_3 * pc_x[k] * msh_957[k];

        t_1273[k] = f_4 * msg0_688[k]
                    - f_5 * msg1_688[k]
                    + f_3 * pc_x[k] * msh_958[k];
    }

#pragma omp simd aligned(t_1274, t_1275, t_1276, t_1277, t_1278, t_1279, pc_x, msg0_689, \
                         msg1_689, msh_959, msh_960, msh_961, msh_962, msh_963, \
                         msh_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1274[k] = f_4 * msg0_689[k]
                    - f_5 * msg1_689[k]
                    + f_3 * pc_x[k] * msh_959[k];

        t_1275[k] = f_3 * pc_x[k] * msh_960[k];

        t_1276[k] = f_3 * pc_x[k] * msh_961[k];

        t_1277[k] = f_3 * pc_x[k] * msh_962[k];

        t_1278[k] = f_3 * pc_x[k] * msh_963[k];

        t_1279[k] = f_3 * pc_x[k] * msh_964[k];
    }

#pragma omp simd aligned(t_1280, t_1281, t_1282, t_1283, pc_x, pc_y, pc_z, lsh_771, msg0_685, \
                         msg1_685, msh_960, msh_961, msh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1280[k] = f_3 * pc_x[k] * msh_965[k];

        t_1281[k] = f_0 * lsh_771[k]
                    + f_1 * msg0_685[k]
                    - f_2 * msg1_685[k]
                    + f_3 * pc_y[k] * msh_960[k];

        t_1282[k] = f_3 * pc_z[k] * msh_960[k];

        t_1283[k] = f_4 * msg0_685[k]
                    - f_5 * msg1_685[k]
                    + f_3 * pc_z[k] * msh_961[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, t_1287, pc_y, pc_z, lsh_776, msg0_686, \
                         msg0_687, msg0_689, msg1_686, msg1_687, msg1_689, msh_962, msh_963, \
                         msh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_6 * msg0_686[k]
                    - f_7 * msg1_686[k]
                    + f_3 * pc_z[k] * msh_962[k];

        t_1285[k] = f_8 * msg0_687[k]
                    - f_9 * msg1_687[k]
                    + f_3 * pc_z[k] * msh_963[k];

        t_1286[k] = f_0 * lsh_776[k]
                    + f_3 * pc_y[k] * msh_965[k];

        t_1287[k] = f_1 * msg0_689[k]
                    - f_2 * msg1_689[k]
                    + f_3 * pc_z[k] * msh_965[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsi0,
                                                           const size_t lsh, const size_t lsi1,
                                                           const size_t msg0, const size_t msg1,
                                                           const size_t msh, const size_t ncols,
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
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsi0_1008 = buffer.data(lsi0 + 1008);
    const auto *lsi0_1009 = buffer.data(lsi0 + 1009);
    const auto *lsi0_1011 = buffer.data(lsi0 + 1011);
    const auto *lsi0_1014 = buffer.data(lsi0 + 1014);
    const auto *lsi0_1018 = buffer.data(lsi0 + 1018);
    const auto *lsi0_1029 = buffer.data(lsi0 + 1029);
    const auto *lsi0_1031 = buffer.data(lsi0 + 1031);
    const auto *lsi0_1032 = buffer.data(lsi0 + 1032);
    const auto *lsi0_1033 = buffer.data(lsi0 + 1033);

    const auto *lsh_771 = buffer.data(lsh + 771);
    const auto *lsh_772 = buffer.data(lsh + 772);
    const auto *lsh_773 = buffer.data(lsh + 773);
    const auto *lsh_774 = buffer.data(lsh + 774);
    const auto *lsh_776 = buffer.data(lsh + 776);
    const auto *lsh_792 = buffer.data(lsh + 792);
    const auto *lsh_797 = buffer.data(lsh + 797);
    const auto *lsh_813 = buffer.data(lsh + 813);
    const auto *lsh_815 = buffer.data(lsh + 815);
    const auto *lsh_816 = buffer.data(lsh + 816);
    const auto *lsh_817 = buffer.data(lsh + 817);
    const auto *lsh_818 = buffer.data(lsh + 818);
    const auto *lsh_834 = buffer.data(lsh + 834);
    const auto *lsh_836 = buffer.data(lsh + 836);
    const auto *lsh_837 = buffer.data(lsh + 837);
    const auto *lsh_838 = buffer.data(lsh + 838);
    const auto *lsh_839 = buffer.data(lsh + 839);
    const auto *lsh_855 = buffer.data(lsh + 855);
    const auto *lsh_857 = buffer.data(lsh + 857);
    const auto *lsh_858 = buffer.data(lsh + 858);
    const auto *lsh_859 = buffer.data(lsh + 859);
    const auto *lsh_860 = buffer.data(lsh + 860);

    const auto *lsi1_1008 = buffer.data(lsi1 + 1008);
    const auto *lsi1_1009 = buffer.data(lsi1 + 1009);
    const auto *lsi1_1011 = buffer.data(lsi1 + 1011);
    const auto *lsi1_1014 = buffer.data(lsi1 + 1014);
    const auto *lsi1_1018 = buffer.data(lsi1 + 1018);
    const auto *lsi1_1029 = buffer.data(lsi1 + 1029);
    const auto *lsi1_1031 = buffer.data(lsi1 + 1031);
    const auto *lsi1_1032 = buffer.data(lsi1 + 1032);
    const auto *lsi1_1033 = buffer.data(lsi1 + 1033);

    const auto *msg0_692 = buffer.data(msg0 + 692);
    const auto *msg0_694 = buffer.data(msg0 + 694);
    const auto *msg0_695 = buffer.data(msg0 + 695);
    const auto *msg0_697 = buffer.data(msg0 + 697);
    const auto *msg0_698 = buffer.data(msg0 + 698);
    const auto *msg0_699 = buffer.data(msg0 + 699);
    const auto *msg0_701 = buffer.data(msg0 + 701);
    const auto *msg0_702 = buffer.data(msg0 + 702);
    const auto *msg0_703 = buffer.data(msg0 + 703);
    const auto *msg0_704 = buffer.data(msg0 + 704);
    const auto *msg0_705 = buffer.data(msg0 + 705);
    const auto *msg0_706 = buffer.data(msg0 + 706);
    const auto *msg0_707 = buffer.data(msg0 + 707);
    const auto *msg0_708 = buffer.data(msg0 + 708);
    const auto *msg0_709 = buffer.data(msg0 + 709);
    const auto *msg0_710 = buffer.data(msg0 + 710);
    const auto *msg0_711 = buffer.data(msg0 + 711);
    const auto *msg0_712 = buffer.data(msg0 + 712);
    const auto *msg0_713 = buffer.data(msg0 + 713);
    const auto *msg0_714 = buffer.data(msg0 + 714);
    const auto *msg0_715 = buffer.data(msg0 + 715);
    const auto *msg0_716 = buffer.data(msg0 + 716);
    const auto *msg0_717 = buffer.data(msg0 + 717);
    const auto *msg0_718 = buffer.data(msg0 + 718);
    const auto *msg0_719 = buffer.data(msg0 + 719);
    const auto *msg0_720 = buffer.data(msg0 + 720);
    const auto *msg0_721 = buffer.data(msg0 + 721);
    const auto *msg0_722 = buffer.data(msg0 + 722);
    const auto *msg0_723 = buffer.data(msg0 + 723);
    const auto *msg0_724 = buffer.data(msg0 + 724);
    const auto *msg0_725 = buffer.data(msg0 + 725);
    const auto *msg0_726 = buffer.data(msg0 + 726);
    const auto *msg0_727 = buffer.data(msg0 + 727);
    const auto *msg0_728 = buffer.data(msg0 + 728);
    const auto *msg0_729 = buffer.data(msg0 + 729);
    const auto *msg0_730 = buffer.data(msg0 + 730);
    const auto *msg0_731 = buffer.data(msg0 + 731);
    const auto *msg0_732 = buffer.data(msg0 + 732);
    const auto *msg0_733 = buffer.data(msg0 + 733);
    const auto *msg0_734 = buffer.data(msg0 + 734);
    const auto *msg0_735 = buffer.data(msg0 + 735);
    const auto *msg0_736 = buffer.data(msg0 + 736);
    const auto *msg0_737 = buffer.data(msg0 + 737);
    const auto *msg0_738 = buffer.data(msg0 + 738);
    const auto *msg0_739 = buffer.data(msg0 + 739);
    const auto *msg0_740 = buffer.data(msg0 + 740);
    const auto *msg0_741 = buffer.data(msg0 + 741);
    const auto *msg0_742 = buffer.data(msg0 + 742);
    const auto *msg0_743 = buffer.data(msg0 + 743);
    const auto *msg0_744 = buffer.data(msg0 + 744);
    const auto *msg0_745 = buffer.data(msg0 + 745);
    const auto *msg0_746 = buffer.data(msg0 + 746);
    const auto *msg0_747 = buffer.data(msg0 + 747);
    const auto *msg0_748 = buffer.data(msg0 + 748);
    const auto *msg0_749 = buffer.data(msg0 + 749);
    const auto *msg0_750 = buffer.data(msg0 + 750);
    const auto *msg0_751 = buffer.data(msg0 + 751);
    const auto *msg0_752 = buffer.data(msg0 + 752);
    const auto *msg0_753 = buffer.data(msg0 + 753);

    const auto *msg1_692 = buffer.data(msg1 + 692);
    const auto *msg1_694 = buffer.data(msg1 + 694);
    const auto *msg1_695 = buffer.data(msg1 + 695);
    const auto *msg1_697 = buffer.data(msg1 + 697);
    const auto *msg1_698 = buffer.data(msg1 + 698);
    const auto *msg1_699 = buffer.data(msg1 + 699);
    const auto *msg1_701 = buffer.data(msg1 + 701);
    const auto *msg1_702 = buffer.data(msg1 + 702);
    const auto *msg1_703 = buffer.data(msg1 + 703);
    const auto *msg1_704 = buffer.data(msg1 + 704);
    const auto *msg1_705 = buffer.data(msg1 + 705);
    const auto *msg1_706 = buffer.data(msg1 + 706);
    const auto *msg1_707 = buffer.data(msg1 + 707);
    const auto *msg1_708 = buffer.data(msg1 + 708);
    const auto *msg1_709 = buffer.data(msg1 + 709);
    const auto *msg1_710 = buffer.data(msg1 + 710);
    const auto *msg1_711 = buffer.data(msg1 + 711);
    const auto *msg1_712 = buffer.data(msg1 + 712);
    const auto *msg1_713 = buffer.data(msg1 + 713);
    const auto *msg1_714 = buffer.data(msg1 + 714);
    const auto *msg1_715 = buffer.data(msg1 + 715);
    const auto *msg1_716 = buffer.data(msg1 + 716);
    const auto *msg1_717 = buffer.data(msg1 + 717);
    const auto *msg1_718 = buffer.data(msg1 + 718);
    const auto *msg1_719 = buffer.data(msg1 + 719);
    const auto *msg1_720 = buffer.data(msg1 + 720);
    const auto *msg1_721 = buffer.data(msg1 + 721);
    const auto *msg1_722 = buffer.data(msg1 + 722);
    const auto *msg1_723 = buffer.data(msg1 + 723);
    const auto *msg1_724 = buffer.data(msg1 + 724);
    const auto *msg1_725 = buffer.data(msg1 + 725);
    const auto *msg1_726 = buffer.data(msg1 + 726);
    const auto *msg1_727 = buffer.data(msg1 + 727);
    const auto *msg1_728 = buffer.data(msg1 + 728);
    const auto *msg1_729 = buffer.data(msg1 + 729);
    const auto *msg1_730 = buffer.data(msg1 + 730);
    const auto *msg1_731 = buffer.data(msg1 + 731);
    const auto *msg1_732 = buffer.data(msg1 + 732);
    const auto *msg1_733 = buffer.data(msg1 + 733);
    const auto *msg1_734 = buffer.data(msg1 + 734);
    const auto *msg1_735 = buffer.data(msg1 + 735);
    const auto *msg1_736 = buffer.data(msg1 + 736);
    const auto *msg1_737 = buffer.data(msg1 + 737);
    const auto *msg1_738 = buffer.data(msg1 + 738);
    const auto *msg1_739 = buffer.data(msg1 + 739);
    const auto *msg1_740 = buffer.data(msg1 + 740);
    const auto *msg1_741 = buffer.data(msg1 + 741);
    const auto *msg1_742 = buffer.data(msg1 + 742);
    const auto *msg1_743 = buffer.data(msg1 + 743);
    const auto *msg1_744 = buffer.data(msg1 + 744);
    const auto *msg1_745 = buffer.data(msg1 + 745);
    const auto *msg1_746 = buffer.data(msg1 + 746);
    const auto *msg1_747 = buffer.data(msg1 + 747);
    const auto *msg1_748 = buffer.data(msg1 + 748);
    const auto *msg1_749 = buffer.data(msg1 + 749);
    const auto *msg1_750 = buffer.data(msg1 + 750);
    const auto *msg1_751 = buffer.data(msg1 + 751);
    const auto *msg1_752 = buffer.data(msg1 + 752);
    const auto *msg1_753 = buffer.data(msg1 + 753);

    const auto *msh_968 = buffer.data(msh + 968);
    const auto *msh_970 = buffer.data(msh + 970);
    const auto *msh_971 = buffer.data(msh + 971);
    const auto *msh_973 = buffer.data(msh + 973);
    const auto *msh_974 = buffer.data(msh + 974);
    const auto *msh_975 = buffer.data(msh + 975);
    const auto *msh_977 = buffer.data(msh + 977);
    const auto *msh_978 = buffer.data(msh + 978);
    const auto *msh_979 = buffer.data(msh + 979);
    const auto *msh_980 = buffer.data(msh + 980);
    const auto *msh_981 = buffer.data(msh + 981);
    const auto *msh_982 = buffer.data(msh + 982);
    const auto *msh_983 = buffer.data(msh + 983);
    const auto *msh_984 = buffer.data(msh + 984);
    const auto *msh_985 = buffer.data(msh + 985);
    const auto *msh_986 = buffer.data(msh + 986);
    const auto *msh_987 = buffer.data(msh + 987);
    const auto *msh_988 = buffer.data(msh + 988);
    const auto *msh_989 = buffer.data(msh + 989);
    const auto *msh_990 = buffer.data(msh + 990);
    const auto *msh_991 = buffer.data(msh + 991);
    const auto *msh_992 = buffer.data(msh + 992);
    const auto *msh_993 = buffer.data(msh + 993);
    const auto *msh_994 = buffer.data(msh + 994);
    const auto *msh_995 = buffer.data(msh + 995);
    const auto *msh_996 = buffer.data(msh + 996);
    const auto *msh_997 = buffer.data(msh + 997);
    const auto *msh_998 = buffer.data(msh + 998);
    const auto *msh_999 = buffer.data(msh + 999);
    const auto *msh_1000 = buffer.data(msh + 1000);
    const auto *msh_1001 = buffer.data(msh + 1001);
    const auto *msh_1002 = buffer.data(msh + 1002);
    const auto *msh_1003 = buffer.data(msh + 1003);
    const auto *msh_1004 = buffer.data(msh + 1004);
    const auto *msh_1005 = buffer.data(msh + 1005);
    const auto *msh_1006 = buffer.data(msh + 1006);
    const auto *msh_1007 = buffer.data(msh + 1007);
    const auto *msh_1008 = buffer.data(msh + 1008);
    const auto *msh_1009 = buffer.data(msh + 1009);
    const auto *msh_1010 = buffer.data(msh + 1010);
    const auto *msh_1011 = buffer.data(msh + 1011);
    const auto *msh_1012 = buffer.data(msh + 1012);
    const auto *msh_1013 = buffer.data(msh + 1013);
    const auto *msh_1014 = buffer.data(msh + 1014);
    const auto *msh_1015 = buffer.data(msh + 1015);
    const auto *msh_1016 = buffer.data(msh + 1016);
    const auto *msh_1017 = buffer.data(msh + 1017);
    const auto *msh_1018 = buffer.data(msh + 1018);
    const auto *msh_1019 = buffer.data(msh + 1019);
    const auto *msh_1020 = buffer.data(msh + 1020);
    const auto *msh_1021 = buffer.data(msh + 1021);
    const auto *msh_1022 = buffer.data(msh + 1022);
    const auto *msh_1023 = buffer.data(msh + 1023);
    const auto *msh_1024 = buffer.data(msh + 1024);
    const auto *msh_1025 = buffer.data(msh + 1025);
    const auto *msh_1026 = buffer.data(msh + 1026);
    const auto *msh_1027 = buffer.data(msh + 1027);
    const auto *msh_1028 = buffer.data(msh + 1028);
    const auto *msh_1029 = buffer.data(msh + 1029);
    const auto *msh_1030 = buffer.data(msh + 1030);
    const auto *msh_1031 = buffer.data(msh + 1031);
    const auto *msh_1032 = buffer.data(msh + 1032);
    const auto *msh_1033 = buffer.data(msh + 1033);
    const auto *msh_1034 = buffer.data(msh + 1034);
    const auto *msh_1035 = buffer.data(msh + 1035);
    const auto *msh_1036 = buffer.data(msh + 1036);
    const auto *msh_1037 = buffer.data(msh + 1037);
    const auto *msh_1038 = buffer.data(msh + 1038);
    const auto *msh_1039 = buffer.data(msh + 1039);
    const auto *msh_1040 = buffer.data(msh + 1040);
    const auto *msh_1041 = buffer.data(msh + 1041);
    const auto *msh_1042 = buffer.data(msh + 1042);
    const auto *msh_1043 = buffer.data(msh + 1043);
    const auto *msh_1044 = buffer.data(msh + 1044);
    const auto *msh_1045 = buffer.data(msh + 1045);
    const auto *msh_1046 = buffer.data(msh + 1046);
    const auto *msh_1047 = buffer.data(msh + 1047);
    const auto *msh_1048 = buffer.data(msh + 1048);
    const auto *msh_1049 = buffer.data(msh + 1049);
    const auto *msh_1050 = buffer.data(msh + 1050);
    const auto *msh_1051 = buffer.data(msh + 1051);
    const auto *msh_1052 = buffer.data(msh + 1052);
    const auto *msh_1053 = buffer.data(msh + 1053);

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pa_z, pc_x, pc_z, lsi0_1008, \
                         lsi0_1009, lsi0_1011, lsi1_1008, lsi1_1009, lsi1_1011, msg0_692, \
                         msg1_692, msh_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = pa_z[k] * lsi0_1008[k]
                    - f_10 * pc_z[k] * lsi1_1008[k];

        t_1289[k] = pa_z[k] * lsi0_1009[k]
                    - f_10 * pc_z[k] * lsi1_1009[k];

        t_1290[k] = f_16 * msg0_692[k]
                    - f_17 * msg1_692[k]
                    + f_3 * pc_x[k] * msh_968[k];

        t_1291[k] = pa_z[k] * lsi0_1011[k]
                    - f_10 * pc_z[k] * lsi1_1011[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, pa_z, pc_x, pc_z, lsi0_1014, lsi1_1014, \
                         msg0_694, msg0_695, msg1_694, msg1_695, msh_970, \
                         msh_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_8 * msg0_694[k]
                    - f_9 * msg1_694[k]
                    + f_3 * pc_x[k] * msh_970[k];

        t_1293[k] = f_8 * msg0_695[k]
                    - f_9 * msg1_695[k]
                    + f_3 * pc_x[k] * msh_971[k];

        t_1294[k] = pa_z[k] * lsi0_1014[k]
                    - f_10 * pc_z[k] * lsi1_1014[k];
    }

#pragma omp simd aligned(t_1295, t_1296, t_1297, pc_x, msg0_697, msg0_698, msg0_699, msg1_697, \
                         msg1_698, msg1_699, msh_973, msh_974, \
                         msh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1295[k] = f_6 * msg0_697[k]
                    - f_7 * msg1_697[k]
                    + f_3 * pc_x[k] * msh_973[k];

        t_1296[k] = f_6 * msg0_698[k]
                    - f_7 * msg1_698[k]
                    + f_3 * pc_x[k] * msh_974[k];

        t_1297[k] = f_6 * msg0_699[k]
                    - f_7 * msg1_699[k]
                    + f_3 * pc_x[k] * msh_975[k];
    }

#pragma omp simd aligned(t_1298, t_1299, t_1300, pa_z, pc_x, pc_z, lsi0_1018, lsi1_1018, \
                         msg0_701, msg0_702, msg1_701, msg1_702, msh_977, \
                         msh_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1298[k] = pa_z[k] * lsi0_1018[k]
                    - f_10 * pc_z[k] * lsi1_1018[k];

        t_1299[k] = f_4 * msg0_701[k]
                    - f_5 * msg1_701[k]
                    + f_3 * pc_x[k] * msh_977[k];

        t_1300[k] = f_4 * msg0_702[k]
                    - f_5 * msg1_702[k]
                    + f_3 * pc_x[k] * msh_978[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, t_1305, pc_x, msg0_703, msg0_704, \
                         msg1_703, msg1_704, msh_979, msh_980, msh_981, msh_982, \
                         msh_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_4 * msg0_703[k]
                    - f_5 * msg1_703[k]
                    + f_3 * pc_x[k] * msh_979[k];

        t_1302[k] = f_4 * msg0_704[k]
                    - f_5 * msg1_704[k]
                    + f_3 * pc_x[k] * msh_980[k];

        t_1303[k] = f_3 * pc_x[k] * msh_981[k];

        t_1304[k] = f_3 * pc_x[k] * msh_982[k];

        t_1305[k] = f_3 * pc_x[k] * msh_983[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, t_1309, t_1310, pa_z, pc_x, pc_z, lsi0_1029, \
                         lsh_771, lsi1_1029, msh_981, msh_984, msh_985, \
                         msh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_3 * pc_x[k] * msh_984[k];

        t_1307[k] = f_3 * pc_x[k] * msh_985[k];

        t_1308[k] = f_3 * pc_x[k] * msh_986[k];

        t_1309[k] = pa_z[k] * lsi0_1029[k]
                    - f_10 * pc_z[k] * lsi1_1029[k];

        t_1310[k] = f_11 * lsh_771[k]
                    + f_3 * pc_z[k] * msh_981[k];
    }

#pragma omp simd aligned(t_1311, t_1312, t_1313, pa_z, pc_z, lsi0_1031, lsi0_1032, lsi0_1033, \
                         lsh_772, lsh_773, lsh_774, lsi1_1031, lsi1_1032, \
                         lsi1_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1311[k] = pa_z[k] * lsi0_1031[k]
                    + f_12 * lsh_772[k]
                    - f_10 * pc_z[k] * lsi1_1031[k];

        t_1312[k] = pa_z[k] * lsi0_1032[k]
                    + f_13 * lsh_773[k]
                    - f_10 * pc_z[k] * lsi1_1032[k];

        t_1313[k] = pa_z[k] * lsi0_1033[k]
                    + f_14 * lsh_774[k]
                    - f_10 * pc_z[k] * lsi1_1033[k];
    }

#pragma omp simd aligned(t_1314, t_1315, t_1316, pc_x, pc_y, pc_z, lsh_776, lsh_797, msg0_704, \
                         msg0_705, msg1_704, msg1_705, msh_986, \
                         msh_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1314[k] = f_15 * lsh_797[k]
                    + f_3 * pc_y[k] * msh_986[k];

        t_1315[k] = f_11 * lsh_776[k]
                    + f_1 * msg0_704[k]
                    - f_2 * msg1_704[k]
                    + f_3 * pc_z[k] * msh_986[k];

        t_1316[k] = f_1 * msg0_705[k]
                    - f_2 * msg1_705[k]
                    + f_3 * pc_x[k] * msh_987[k];
    }

#pragma omp simd aligned(t_1317, t_1318, t_1319, pc_x, msg0_706, msg0_707, msg0_708, msg1_706, \
                         msg1_707, msg1_708, msh_988, msh_989, \
                         msh_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1317[k] = f_16 * msg0_706[k]
                    - f_17 * msg1_706[k]
                    + f_3 * pc_x[k] * msh_988[k];

        t_1318[k] = f_16 * msg0_707[k]
                    - f_17 * msg1_707[k]
                    + f_3 * pc_x[k] * msh_989[k];

        t_1319[k] = f_8 * msg0_708[k]
                    - f_9 * msg1_708[k]
                    + f_3 * pc_x[k] * msh_990[k];
    }

#pragma omp simd aligned(t_1320, t_1321, t_1322, pc_x, msg0_709, msg0_710, msg0_711, msg1_709, \
                         msg1_710, msg1_711, msh_991, msh_992, \
                         msh_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1320[k] = f_8 * msg0_709[k]
                    - f_9 * msg1_709[k]
                    + f_3 * pc_x[k] * msh_991[k];

        t_1321[k] = f_8 * msg0_710[k]
                    - f_9 * msg1_710[k]
                    + f_3 * pc_x[k] * msh_992[k];

        t_1322[k] = f_6 * msg0_711[k]
                    - f_7 * msg1_711[k]
                    + f_3 * pc_x[k] * msh_993[k];
    }

#pragma omp simd aligned(t_1323, t_1324, t_1325, pc_x, msg0_712, msg0_713, msg0_714, msg1_712, \
                         msg1_713, msg1_714, msh_994, msh_995, \
                         msh_996 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1323[k] = f_6 * msg0_712[k]
                    - f_7 * msg1_712[k]
                    + f_3 * pc_x[k] * msh_994[k];

        t_1324[k] = f_6 * msg0_713[k]
                    - f_7 * msg1_713[k]
                    + f_3 * pc_x[k] * msh_995[k];

        t_1325[k] = f_6 * msg0_714[k]
                    - f_7 * msg1_714[k]
                    + f_3 * pc_x[k] * msh_996[k];
    }

#pragma omp simd aligned(t_1326, t_1327, t_1328, pc_x, msg0_715, msg0_716, msg0_717, msg1_715, \
                         msg1_716, msg1_717, msh_997, msh_998, \
                         msh_999 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1326[k] = f_4 * msg0_715[k]
                    - f_5 * msg1_715[k]
                    + f_3 * pc_x[k] * msh_997[k];

        t_1327[k] = f_4 * msg0_716[k]
                    - f_5 * msg1_716[k]
                    + f_3 * pc_x[k] * msh_998[k];

        t_1328[k] = f_4 * msg0_717[k]
                    - f_5 * msg1_717[k]
                    + f_3 * pc_x[k] * msh_999[k];
    }

#pragma omp simd aligned(t_1329, t_1330, t_1331, t_1332, t_1333, pc_x, msg0_718, msg0_719, \
                         msg1_718, msg1_719, msh_1000, msh_1001, msh_1002, msh_1003, \
                         msh_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1329[k] = f_4 * msg0_718[k]
                    - f_5 * msg1_718[k]
                    + f_3 * pc_x[k] * msh_1000[k];

        t_1330[k] = f_4 * msg0_719[k]
                    - f_5 * msg1_719[k]
                    + f_3 * pc_x[k] * msh_1001[k];

        t_1331[k] = f_3 * pc_x[k] * msh_1002[k];

        t_1332[k] = f_3 * pc_x[k] * msh_1003[k];

        t_1333[k] = f_3 * pc_x[k] * msh_1004[k];
    }

#pragma omp simd aligned(t_1334, t_1335, t_1336, t_1337, t_1338, pc_x, pc_y, pc_z, lsh_792, \
                         lsh_813, msg0_715, msg1_715, msh_1002, msh_1005, msh_1006, \
                         msh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1334[k] = f_3 * pc_x[k] * msh_1005[k];

        t_1335[k] = f_3 * pc_x[k] * msh_1006[k];

        t_1336[k] = f_3 * pc_x[k] * msh_1007[k];

        t_1337[k] = f_18 * lsh_813[k]
                    + f_1 * msg0_715[k]
                    - f_2 * msg1_715[k]
                    + f_3 * pc_y[k] * msh_1002[k];

        t_1338[k] = f_12 * lsh_792[k]
                    + f_3 * pc_z[k] * msh_1002[k];
    }

#pragma omp simd aligned(t_1339, t_1340, t_1341, pc_y, lsh_815, lsh_816, lsh_817, msg0_717, \
                         msg0_718, msg0_719, msg1_717, msg1_718, msg1_719, msh_1004, msh_1005, \
                         msh_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1339[k] = f_18 * lsh_815[k]
                    + f_8 * msg0_717[k]
                    - f_9 * msg1_717[k]
                    + f_3 * pc_y[k] * msh_1004[k];

        t_1340[k] = f_18 * lsh_816[k]
                    + f_6 * msg0_718[k]
                    - f_7 * msg1_718[k]
                    + f_3 * pc_y[k] * msh_1005[k];

        t_1341[k] = f_18 * lsh_817[k]
                    + f_4 * msg0_719[k]
                    - f_5 * msg1_719[k]
                    + f_3 * pc_y[k] * msh_1006[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, pc_x, pc_y, pc_z, lsh_797, lsh_818, msg0_719, \
                         msg0_720, msg1_719, msg1_720, msh_1007, \
                         msh_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_18 * lsh_818[k]
                    + f_3 * pc_y[k] * msh_1007[k];

        t_1343[k] = f_12 * lsh_797[k]
                    + f_1 * msg0_719[k]
                    - f_2 * msg1_719[k]
                    + f_3 * pc_z[k] * msh_1007[k];

        t_1344[k] = f_1 * msg0_720[k]
                    - f_2 * msg1_720[k]
                    + f_3 * pc_x[k] * msh_1008[k];
    }

#pragma omp simd aligned(t_1345, t_1346, t_1347, pc_x, msg0_721, msg0_722, msg0_723, msg1_721, \
                         msg1_722, msg1_723, msh_1009, msh_1010, \
                         msh_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1345[k] = f_16 * msg0_721[k]
                    - f_17 * msg1_721[k]
                    + f_3 * pc_x[k] * msh_1009[k];

        t_1346[k] = f_16 * msg0_722[k]
                    - f_17 * msg1_722[k]
                    + f_3 * pc_x[k] * msh_1010[k];

        t_1347[k] = f_8 * msg0_723[k]
                    - f_9 * msg1_723[k]
                    + f_3 * pc_x[k] * msh_1011[k];
    }

#pragma omp simd aligned(t_1348, t_1349, t_1350, pc_x, msg0_724, msg0_725, msg0_726, msg1_724, \
                         msg1_725, msg1_726, msh_1012, msh_1013, \
                         msh_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1348[k] = f_8 * msg0_724[k]
                    - f_9 * msg1_724[k]
                    + f_3 * pc_x[k] * msh_1012[k];

        t_1349[k] = f_8 * msg0_725[k]
                    - f_9 * msg1_725[k]
                    + f_3 * pc_x[k] * msh_1013[k];

        t_1350[k] = f_6 * msg0_726[k]
                    - f_7 * msg1_726[k]
                    + f_3 * pc_x[k] * msh_1014[k];
    }

#pragma omp simd aligned(t_1351, t_1352, t_1353, pc_x, msg0_727, msg0_728, msg0_729, msg1_727, \
                         msg1_728, msg1_729, msh_1015, msh_1016, \
                         msh_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1351[k] = f_6 * msg0_727[k]
                    - f_7 * msg1_727[k]
                    + f_3 * pc_x[k] * msh_1015[k];

        t_1352[k] = f_6 * msg0_728[k]
                    - f_7 * msg1_728[k]
                    + f_3 * pc_x[k] * msh_1016[k];

        t_1353[k] = f_6 * msg0_729[k]
                    - f_7 * msg1_729[k]
                    + f_3 * pc_x[k] * msh_1017[k];
    }

#pragma omp simd aligned(t_1354, t_1355, t_1356, pc_x, msg0_730, msg0_731, msg0_732, msg1_730, \
                         msg1_731, msg1_732, msh_1018, msh_1019, \
                         msh_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1354[k] = f_4 * msg0_730[k]
                    - f_5 * msg1_730[k]
                    + f_3 * pc_x[k] * msh_1018[k];

        t_1355[k] = f_4 * msg0_731[k]
                    - f_5 * msg1_731[k]
                    + f_3 * pc_x[k] * msh_1019[k];

        t_1356[k] = f_4 * msg0_732[k]
                    - f_5 * msg1_732[k]
                    + f_3 * pc_x[k] * msh_1020[k];
    }

#pragma omp simd aligned(t_1357, t_1358, t_1359, t_1360, t_1361, pc_x, msg0_733, msg0_734, \
                         msg1_733, msg1_734, msh_1021, msh_1022, msh_1023, msh_1024, \
                         msh_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1357[k] = f_4 * msg0_733[k]
                    - f_5 * msg1_733[k]
                    + f_3 * pc_x[k] * msh_1021[k];

        t_1358[k] = f_4 * msg0_734[k]
                    - f_5 * msg1_734[k]
                    + f_3 * pc_x[k] * msh_1022[k];

        t_1359[k] = f_3 * pc_x[k] * msh_1023[k];

        t_1360[k] = f_3 * pc_x[k] * msh_1024[k];

        t_1361[k] = f_3 * pc_x[k] * msh_1025[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, t_1365, t_1366, pc_x, pc_y, pc_z, lsh_813, \
                         lsh_834, msg0_730, msg1_730, msh_1023, msh_1026, msh_1027, \
                         msh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = f_3 * pc_x[k] * msh_1026[k];

        t_1363[k] = f_3 * pc_x[k] * msh_1027[k];

        t_1364[k] = f_3 * pc_x[k] * msh_1028[k];

        t_1365[k] = f_19 * lsh_834[k]
                    + f_1 * msg0_730[k]
                    - f_2 * msg1_730[k]
                    + f_3 * pc_y[k] * msh_1023[k];

        t_1366[k] = f_13 * lsh_813[k]
                    + f_3 * pc_z[k] * msh_1023[k];
    }

#pragma omp simd aligned(t_1367, t_1368, t_1369, pc_y, lsh_836, lsh_837, lsh_838, msg0_732, \
                         msg0_733, msg0_734, msg1_732, msg1_733, msg1_734, msh_1025, msh_1026, \
                         msh_1027 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1367[k] = f_19 * lsh_836[k]
                    + f_8 * msg0_732[k]
                    - f_9 * msg1_732[k]
                    + f_3 * pc_y[k] * msh_1025[k];

        t_1368[k] = f_19 * lsh_837[k]
                    + f_6 * msg0_733[k]
                    - f_7 * msg1_733[k]
                    + f_3 * pc_y[k] * msh_1026[k];

        t_1369[k] = f_19 * lsh_838[k]
                    + f_4 * msg0_734[k]
                    - f_5 * msg1_734[k]
                    + f_3 * pc_y[k] * msh_1027[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, pc_x, pc_y, pc_z, lsh_818, lsh_839, msg0_734, \
                         msg0_735, msg1_734, msg1_735, msh_1028, \
                         msh_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_19 * lsh_839[k]
                    + f_3 * pc_y[k] * msh_1028[k];

        t_1371[k] = f_13 * lsh_818[k]
                    + f_1 * msg0_734[k]
                    - f_2 * msg1_734[k]
                    + f_3 * pc_z[k] * msh_1028[k];

        t_1372[k] = f_1 * msg0_735[k]
                    - f_2 * msg1_735[k]
                    + f_3 * pc_x[k] * msh_1029[k];
    }

#pragma omp simd aligned(t_1373, t_1374, t_1375, pc_x, msg0_736, msg0_737, msg0_738, msg1_736, \
                         msg1_737, msg1_738, msh_1030, msh_1031, \
                         msh_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1373[k] = f_16 * msg0_736[k]
                    - f_17 * msg1_736[k]
                    + f_3 * pc_x[k] * msh_1030[k];

        t_1374[k] = f_16 * msg0_737[k]
                    - f_17 * msg1_737[k]
                    + f_3 * pc_x[k] * msh_1031[k];

        t_1375[k] = f_8 * msg0_738[k]
                    - f_9 * msg1_738[k]
                    + f_3 * pc_x[k] * msh_1032[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, pc_x, msg0_739, msg0_740, msg0_741, msg1_739, \
                         msg1_740, msg1_741, msh_1033, msh_1034, \
                         msh_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_8 * msg0_739[k]
                    - f_9 * msg1_739[k]
                    + f_3 * pc_x[k] * msh_1033[k];

        t_1377[k] = f_8 * msg0_740[k]
                    - f_9 * msg1_740[k]
                    + f_3 * pc_x[k] * msh_1034[k];

        t_1378[k] = f_6 * msg0_741[k]
                    - f_7 * msg1_741[k]
                    + f_3 * pc_x[k] * msh_1035[k];
    }

#pragma omp simd aligned(t_1379, t_1380, t_1381, pc_x, msg0_742, msg0_743, msg0_744, msg1_742, \
                         msg1_743, msg1_744, msh_1036, msh_1037, \
                         msh_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1379[k] = f_6 * msg0_742[k]
                    - f_7 * msg1_742[k]
                    + f_3 * pc_x[k] * msh_1036[k];

        t_1380[k] = f_6 * msg0_743[k]
                    - f_7 * msg1_743[k]
                    + f_3 * pc_x[k] * msh_1037[k];

        t_1381[k] = f_6 * msg0_744[k]
                    - f_7 * msg1_744[k]
                    + f_3 * pc_x[k] * msh_1038[k];
    }

#pragma omp simd aligned(t_1382, t_1383, t_1384, pc_x, msg0_745, msg0_746, msg0_747, msg1_745, \
                         msg1_746, msg1_747, msh_1039, msh_1040, \
                         msh_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1382[k] = f_4 * msg0_745[k]
                    - f_5 * msg1_745[k]
                    + f_3 * pc_x[k] * msh_1039[k];

        t_1383[k] = f_4 * msg0_746[k]
                    - f_5 * msg1_746[k]
                    + f_3 * pc_x[k] * msh_1040[k];

        t_1384[k] = f_4 * msg0_747[k]
                    - f_5 * msg1_747[k]
                    + f_3 * pc_x[k] * msh_1041[k];
    }

#pragma omp simd aligned(t_1385, t_1386, t_1387, t_1388, t_1389, pc_x, msg0_748, msg0_749, \
                         msg1_748, msg1_749, msh_1042, msh_1043, msh_1044, msh_1045, \
                         msh_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1385[k] = f_4 * msg0_748[k]
                    - f_5 * msg1_748[k]
                    + f_3 * pc_x[k] * msh_1042[k];

        t_1386[k] = f_4 * msg0_749[k]
                    - f_5 * msg1_749[k]
                    + f_3 * pc_x[k] * msh_1043[k];

        t_1387[k] = f_3 * pc_x[k] * msh_1044[k];

        t_1388[k] = f_3 * pc_x[k] * msh_1045[k];

        t_1389[k] = f_3 * pc_x[k] * msh_1046[k];
    }

#pragma omp simd aligned(t_1390, t_1391, t_1392, t_1393, t_1394, pc_x, pc_y, pc_z, lsh_834, \
                         lsh_855, msg0_745, msg1_745, msh_1044, msh_1047, msh_1048, \
                         msh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1390[k] = f_3 * pc_x[k] * msh_1047[k];

        t_1391[k] = f_3 * pc_x[k] * msh_1048[k];

        t_1392[k] = f_3 * pc_x[k] * msh_1049[k];

        t_1393[k] = f_20 * lsh_855[k]
                    + f_1 * msg0_745[k]
                    - f_2 * msg1_745[k]
                    + f_3 * pc_y[k] * msh_1044[k];

        t_1394[k] = f_14 * lsh_834[k]
                    + f_3 * pc_z[k] * msh_1044[k];
    }

#pragma omp simd aligned(t_1395, t_1396, t_1397, pc_y, lsh_857, lsh_858, lsh_859, msg0_747, \
                         msg0_748, msg0_749, msg1_747, msg1_748, msg1_749, msh_1046, msh_1047, \
                         msh_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1395[k] = f_20 * lsh_857[k]
                    + f_8 * msg0_747[k]
                    - f_9 * msg1_747[k]
                    + f_3 * pc_y[k] * msh_1046[k];

        t_1396[k] = f_20 * lsh_858[k]
                    + f_6 * msg0_748[k]
                    - f_7 * msg1_748[k]
                    + f_3 * pc_y[k] * msh_1047[k];

        t_1397[k] = f_20 * lsh_859[k]
                    + f_4 * msg0_749[k]
                    - f_5 * msg1_749[k]
                    + f_3 * pc_y[k] * msh_1048[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, pc_x, pc_y, pc_z, lsh_839, lsh_860, msg0_749, \
                         msg0_750, msg1_749, msg1_750, msh_1049, \
                         msh_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_20 * lsh_860[k]
                    + f_3 * pc_y[k] * msh_1049[k];

        t_1399[k] = f_14 * lsh_839[k]
                    + f_1 * msg0_749[k]
                    - f_2 * msg1_749[k]
                    + f_3 * pc_z[k] * msh_1049[k];

        t_1400[k] = f_1 * msg0_750[k]
                    - f_2 * msg1_750[k]
                    + f_3 * pc_x[k] * msh_1050[k];
    }

#pragma omp simd aligned(t_1401, t_1402, t_1403, pc_x, msg0_751, msg0_752, msg0_753, msg1_751, \
                         msg1_752, msg1_753, msh_1051, msh_1052, \
                         msh_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1401[k] = f_16 * msg0_751[k]
                    - f_17 * msg1_751[k]
                    + f_3 * pc_x[k] * msh_1051[k];

        t_1402[k] = f_16 * msg0_752[k]
                    - f_17 * msg1_752[k]
                    + f_3 * pc_x[k] * msh_1052[k];

        t_1403[k] = f_8 * msg0_753[k]
                    - f_9 * msg1_753[k]
                    + f_3 * pc_x[k] * msh_1053[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t lsi0,
                                                           const size_t lsh, const size_t lsi1,
                                                           const size_t msg0, const size_t msg1,
                                                           const size_t msh, const size_t ncols,
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
    const auto f_15 = 4.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsi0_1232 = buffer.data(lsi0 + 1232);
    const auto *lsi0_1234 = buffer.data(lsi0 + 1234);
    const auto *lsi0_1237 = buffer.data(lsi0 + 1237);
    const auto *lsi0_1241 = buffer.data(lsi0 + 1241);
    const auto *lsi0_1246 = buffer.data(lsi0 + 1246);
    const auto *lsi0_1253 = buffer.data(lsi0 + 1253);
    const auto *lsi0_1255 = buffer.data(lsi0 + 1255);
    const auto *lsi0_1256 = buffer.data(lsi0 + 1256);
    const auto *lsi0_1257 = buffer.data(lsi0 + 1257);
    const auto *lsi0_1259 = buffer.data(lsi0 + 1259);

    const auto *lsh_855 = buffer.data(lsh + 855);
    const auto *lsh_860 = buffer.data(lsh + 860);
    const auto *lsh_876 = buffer.data(lsh + 876);
    const auto *lsh_878 = buffer.data(lsh + 878);
    const auto *lsh_879 = buffer.data(lsh + 879);
    const auto *lsh_880 = buffer.data(lsh + 880);
    const auto *lsh_881 = buffer.data(lsh + 881);
    const auto *lsh_897 = buffer.data(lsh + 897);
    const auto *lsh_899 = buffer.data(lsh + 899);
    const auto *lsh_900 = buffer.data(lsh + 900);
    const auto *lsh_901 = buffer.data(lsh + 901);
    const auto *lsh_902 = buffer.data(lsh + 902);
    const auto *lsh_918 = buffer.data(lsh + 918);
    const auto *lsh_920 = buffer.data(lsh + 920);
    const auto *lsh_921 = buffer.data(lsh + 921);
    const auto *lsh_922 = buffer.data(lsh + 922);
    const auto *lsh_923 = buffer.data(lsh + 923);
    const auto *lsh_939 = buffer.data(lsh + 939);
    const auto *lsh_941 = buffer.data(lsh + 941);
    const auto *lsh_942 = buffer.data(lsh + 942);
    const auto *lsh_943 = buffer.data(lsh + 943);
    const auto *lsh_944 = buffer.data(lsh + 944);

    const auto *lsi1_1232 = buffer.data(lsi1 + 1232);
    const auto *lsi1_1234 = buffer.data(lsi1 + 1234);
    const auto *lsi1_1237 = buffer.data(lsi1 + 1237);
    const auto *lsi1_1241 = buffer.data(lsi1 + 1241);
    const auto *lsi1_1246 = buffer.data(lsi1 + 1246);
    const auto *lsi1_1253 = buffer.data(lsi1 + 1253);
    const auto *lsi1_1255 = buffer.data(lsi1 + 1255);
    const auto *lsi1_1256 = buffer.data(lsi1 + 1256);
    const auto *lsi1_1257 = buffer.data(lsi1 + 1257);
    const auto *lsi1_1259 = buffer.data(lsi1 + 1259);

    const auto *msg0_754 = buffer.data(msg0 + 754);
    const auto *msg0_755 = buffer.data(msg0 + 755);
    const auto *msg0_756 = buffer.data(msg0 + 756);
    const auto *msg0_757 = buffer.data(msg0 + 757);
    const auto *msg0_758 = buffer.data(msg0 + 758);
    const auto *msg0_759 = buffer.data(msg0 + 759);
    const auto *msg0_760 = buffer.data(msg0 + 760);
    const auto *msg0_761 = buffer.data(msg0 + 761);
    const auto *msg0_762 = buffer.data(msg0 + 762);
    const auto *msg0_763 = buffer.data(msg0 + 763);
    const auto *msg0_764 = buffer.data(msg0 + 764);
    const auto *msg0_765 = buffer.data(msg0 + 765);
    const auto *msg0_766 = buffer.data(msg0 + 766);
    const auto *msg0_767 = buffer.data(msg0 + 767);
    const auto *msg0_768 = buffer.data(msg0 + 768);
    const auto *msg0_769 = buffer.data(msg0 + 769);
    const auto *msg0_770 = buffer.data(msg0 + 770);
    const auto *msg0_771 = buffer.data(msg0 + 771);
    const auto *msg0_772 = buffer.data(msg0 + 772);
    const auto *msg0_773 = buffer.data(msg0 + 773);
    const auto *msg0_774 = buffer.data(msg0 + 774);
    const auto *msg0_775 = buffer.data(msg0 + 775);
    const auto *msg0_776 = buffer.data(msg0 + 776);
    const auto *msg0_777 = buffer.data(msg0 + 777);
    const auto *msg0_778 = buffer.data(msg0 + 778);
    const auto *msg0_779 = buffer.data(msg0 + 779);
    const auto *msg0_780 = buffer.data(msg0 + 780);
    const auto *msg0_781 = buffer.data(msg0 + 781);
    const auto *msg0_782 = buffer.data(msg0 + 782);
    const auto *msg0_783 = buffer.data(msg0 + 783);
    const auto *msg0_784 = buffer.data(msg0 + 784);
    const auto *msg0_785 = buffer.data(msg0 + 785);
    const auto *msg0_786 = buffer.data(msg0 + 786);
    const auto *msg0_787 = buffer.data(msg0 + 787);
    const auto *msg0_788 = buffer.data(msg0 + 788);
    const auto *msg0_789 = buffer.data(msg0 + 789);
    const auto *msg0_790 = buffer.data(msg0 + 790);
    const auto *msg0_791 = buffer.data(msg0 + 791);
    const auto *msg0_792 = buffer.data(msg0 + 792);
    const auto *msg0_793 = buffer.data(msg0 + 793);
    const auto *msg0_794 = buffer.data(msg0 + 794);
    const auto *msg0_796 = buffer.data(msg0 + 796);
    const auto *msg0_798 = buffer.data(msg0 + 798);
    const auto *msg0_799 = buffer.data(msg0 + 799);
    const auto *msg0_801 = buffer.data(msg0 + 801);
    const auto *msg0_802 = buffer.data(msg0 + 802);
    const auto *msg0_803 = buffer.data(msg0 + 803);
    const auto *msg0_805 = buffer.data(msg0 + 805);
    const auto *msg0_806 = buffer.data(msg0 + 806);
    const auto *msg0_807 = buffer.data(msg0 + 807);
    const auto *msg0_808 = buffer.data(msg0 + 808);
    const auto *msg0_810 = buffer.data(msg0 + 810);
    const auto *msg0_812 = buffer.data(msg0 + 812);
    const auto *msg0_813 = buffer.data(msg0 + 813);
    const auto *msg0_815 = buffer.data(msg0 + 815);
    const auto *msg0_816 = buffer.data(msg0 + 816);
    const auto *msg0_817 = buffer.data(msg0 + 817);

    const auto *msg1_754 = buffer.data(msg1 + 754);
    const auto *msg1_755 = buffer.data(msg1 + 755);
    const auto *msg1_756 = buffer.data(msg1 + 756);
    const auto *msg1_757 = buffer.data(msg1 + 757);
    const auto *msg1_758 = buffer.data(msg1 + 758);
    const auto *msg1_759 = buffer.data(msg1 + 759);
    const auto *msg1_760 = buffer.data(msg1 + 760);
    const auto *msg1_761 = buffer.data(msg1 + 761);
    const auto *msg1_762 = buffer.data(msg1 + 762);
    const auto *msg1_763 = buffer.data(msg1 + 763);
    const auto *msg1_764 = buffer.data(msg1 + 764);
    const auto *msg1_765 = buffer.data(msg1 + 765);
    const auto *msg1_766 = buffer.data(msg1 + 766);
    const auto *msg1_767 = buffer.data(msg1 + 767);
    const auto *msg1_768 = buffer.data(msg1 + 768);
    const auto *msg1_769 = buffer.data(msg1 + 769);
    const auto *msg1_770 = buffer.data(msg1 + 770);
    const auto *msg1_771 = buffer.data(msg1 + 771);
    const auto *msg1_772 = buffer.data(msg1 + 772);
    const auto *msg1_773 = buffer.data(msg1 + 773);
    const auto *msg1_774 = buffer.data(msg1 + 774);
    const auto *msg1_775 = buffer.data(msg1 + 775);
    const auto *msg1_776 = buffer.data(msg1 + 776);
    const auto *msg1_777 = buffer.data(msg1 + 777);
    const auto *msg1_778 = buffer.data(msg1 + 778);
    const auto *msg1_779 = buffer.data(msg1 + 779);
    const auto *msg1_780 = buffer.data(msg1 + 780);
    const auto *msg1_781 = buffer.data(msg1 + 781);
    const auto *msg1_782 = buffer.data(msg1 + 782);
    const auto *msg1_783 = buffer.data(msg1 + 783);
    const auto *msg1_784 = buffer.data(msg1 + 784);
    const auto *msg1_785 = buffer.data(msg1 + 785);
    const auto *msg1_786 = buffer.data(msg1 + 786);
    const auto *msg1_787 = buffer.data(msg1 + 787);
    const auto *msg1_788 = buffer.data(msg1 + 788);
    const auto *msg1_789 = buffer.data(msg1 + 789);
    const auto *msg1_790 = buffer.data(msg1 + 790);
    const auto *msg1_791 = buffer.data(msg1 + 791);
    const auto *msg1_792 = buffer.data(msg1 + 792);
    const auto *msg1_793 = buffer.data(msg1 + 793);
    const auto *msg1_794 = buffer.data(msg1 + 794);
    const auto *msg1_796 = buffer.data(msg1 + 796);
    const auto *msg1_798 = buffer.data(msg1 + 798);
    const auto *msg1_799 = buffer.data(msg1 + 799);
    const auto *msg1_801 = buffer.data(msg1 + 801);
    const auto *msg1_802 = buffer.data(msg1 + 802);
    const auto *msg1_803 = buffer.data(msg1 + 803);
    const auto *msg1_805 = buffer.data(msg1 + 805);
    const auto *msg1_806 = buffer.data(msg1 + 806);
    const auto *msg1_807 = buffer.data(msg1 + 807);
    const auto *msg1_808 = buffer.data(msg1 + 808);
    const auto *msg1_810 = buffer.data(msg1 + 810);
    const auto *msg1_812 = buffer.data(msg1 + 812);
    const auto *msg1_813 = buffer.data(msg1 + 813);
    const auto *msg1_815 = buffer.data(msg1 + 815);
    const auto *msg1_816 = buffer.data(msg1 + 816);
    const auto *msg1_817 = buffer.data(msg1 + 817);

    const auto *msh_1054 = buffer.data(msh + 1054);
    const auto *msh_1055 = buffer.data(msh + 1055);
    const auto *msh_1056 = buffer.data(msh + 1056);
    const auto *msh_1057 = buffer.data(msh + 1057);
    const auto *msh_1058 = buffer.data(msh + 1058);
    const auto *msh_1059 = buffer.data(msh + 1059);
    const auto *msh_1060 = buffer.data(msh + 1060);
    const auto *msh_1061 = buffer.data(msh + 1061);
    const auto *msh_1062 = buffer.data(msh + 1062);
    const auto *msh_1063 = buffer.data(msh + 1063);
    const auto *msh_1064 = buffer.data(msh + 1064);
    const auto *msh_1065 = buffer.data(msh + 1065);
    const auto *msh_1066 = buffer.data(msh + 1066);
    const auto *msh_1067 = buffer.data(msh + 1067);
    const auto *msh_1068 = buffer.data(msh + 1068);
    const auto *msh_1069 = buffer.data(msh + 1069);
    const auto *msh_1070 = buffer.data(msh + 1070);
    const auto *msh_1071 = buffer.data(msh + 1071);
    const auto *msh_1072 = buffer.data(msh + 1072);
    const auto *msh_1073 = buffer.data(msh + 1073);
    const auto *msh_1074 = buffer.data(msh + 1074);
    const auto *msh_1075 = buffer.data(msh + 1075);
    const auto *msh_1076 = buffer.data(msh + 1076);
    const auto *msh_1077 = buffer.data(msh + 1077);
    const auto *msh_1078 = buffer.data(msh + 1078);
    const auto *msh_1079 = buffer.data(msh + 1079);
    const auto *msh_1080 = buffer.data(msh + 1080);
    const auto *msh_1081 = buffer.data(msh + 1081);
    const auto *msh_1082 = buffer.data(msh + 1082);
    const auto *msh_1083 = buffer.data(msh + 1083);
    const auto *msh_1084 = buffer.data(msh + 1084);
    const auto *msh_1085 = buffer.data(msh + 1085);
    const auto *msh_1086 = buffer.data(msh + 1086);
    const auto *msh_1087 = buffer.data(msh + 1087);
    const auto *msh_1088 = buffer.data(msh + 1088);
    const auto *msh_1089 = buffer.data(msh + 1089);
    const auto *msh_1090 = buffer.data(msh + 1090);
    const auto *msh_1091 = buffer.data(msh + 1091);
    const auto *msh_1092 = buffer.data(msh + 1092);
    const auto *msh_1093 = buffer.data(msh + 1093);
    const auto *msh_1094 = buffer.data(msh + 1094);
    const auto *msh_1095 = buffer.data(msh + 1095);
    const auto *msh_1096 = buffer.data(msh + 1096);
    const auto *msh_1097 = buffer.data(msh + 1097);
    const auto *msh_1098 = buffer.data(msh + 1098);
    const auto *msh_1099 = buffer.data(msh + 1099);
    const auto *msh_1100 = buffer.data(msh + 1100);
    const auto *msh_1101 = buffer.data(msh + 1101);
    const auto *msh_1102 = buffer.data(msh + 1102);
    const auto *msh_1103 = buffer.data(msh + 1103);
    const auto *msh_1104 = buffer.data(msh + 1104);
    const auto *msh_1105 = buffer.data(msh + 1105);
    const auto *msh_1106 = buffer.data(msh + 1106);
    const auto *msh_1107 = buffer.data(msh + 1107);
    const auto *msh_1108 = buffer.data(msh + 1108);
    const auto *msh_1109 = buffer.data(msh + 1109);
    const auto *msh_1110 = buffer.data(msh + 1110);
    const auto *msh_1111 = buffer.data(msh + 1111);
    const auto *msh_1112 = buffer.data(msh + 1112);
    const auto *msh_1114 = buffer.data(msh + 1114);
    const auto *msh_1116 = buffer.data(msh + 1116);
    const auto *msh_1117 = buffer.data(msh + 1117);
    const auto *msh_1119 = buffer.data(msh + 1119);
    const auto *msh_1120 = buffer.data(msh + 1120);
    const auto *msh_1121 = buffer.data(msh + 1121);
    const auto *msh_1123 = buffer.data(msh + 1123);
    const auto *msh_1124 = buffer.data(msh + 1124);
    const auto *msh_1125 = buffer.data(msh + 1125);
    const auto *msh_1126 = buffer.data(msh + 1126);
    const auto *msh_1128 = buffer.data(msh + 1128);
    const auto *msh_1129 = buffer.data(msh + 1129);
    const auto *msh_1130 = buffer.data(msh + 1130);
    const auto *msh_1131 = buffer.data(msh + 1131);
    const auto *msh_1132 = buffer.data(msh + 1132);
    const auto *msh_1133 = buffer.data(msh + 1133);
    const auto *msh_1134 = buffer.data(msh + 1134);
    const auto *msh_1136 = buffer.data(msh + 1136);
    const auto *msh_1137 = buffer.data(msh + 1137);
    const auto *msh_1139 = buffer.data(msh + 1139);
    const auto *msh_1140 = buffer.data(msh + 1140);
    const auto *msh_1141 = buffer.data(msh + 1141);

#pragma omp simd aligned(t_1404, t_1405, t_1406, pc_x, msg0_754, msg0_755, msg0_756, msg1_754, \
                         msg1_755, msg1_756, msh_1054, msh_1055, \
                         msh_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1404[k] = f_8 * msg0_754[k]
                    - f_9 * msg1_754[k]
                    + f_3 * pc_x[k] * msh_1054[k];

        t_1405[k] = f_8 * msg0_755[k]
                    - f_9 * msg1_755[k]
                    + f_3 * pc_x[k] * msh_1055[k];

        t_1406[k] = f_6 * msg0_756[k]
                    - f_7 * msg1_756[k]
                    + f_3 * pc_x[k] * msh_1056[k];
    }

#pragma omp simd aligned(t_1407, t_1408, t_1409, pc_x, msg0_757, msg0_758, msg0_759, msg1_757, \
                         msg1_758, msg1_759, msh_1057, msh_1058, \
                         msh_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1407[k] = f_6 * msg0_757[k]
                    - f_7 * msg1_757[k]
                    + f_3 * pc_x[k] * msh_1057[k];

        t_1408[k] = f_6 * msg0_758[k]
                    - f_7 * msg1_758[k]
                    + f_3 * pc_x[k] * msh_1058[k];

        t_1409[k] = f_6 * msg0_759[k]
                    - f_7 * msg1_759[k]
                    + f_3 * pc_x[k] * msh_1059[k];
    }

#pragma omp simd aligned(t_1410, t_1411, t_1412, pc_x, msg0_760, msg0_761, msg0_762, msg1_760, \
                         msg1_761, msg1_762, msh_1060, msh_1061, \
                         msh_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1410[k] = f_4 * msg0_760[k]
                    - f_5 * msg1_760[k]
                    + f_3 * pc_x[k] * msh_1060[k];

        t_1411[k] = f_4 * msg0_761[k]
                    - f_5 * msg1_761[k]
                    + f_3 * pc_x[k] * msh_1061[k];

        t_1412[k] = f_4 * msg0_762[k]
                    - f_5 * msg1_762[k]
                    + f_3 * pc_x[k] * msh_1062[k];
    }

#pragma omp simd aligned(t_1413, t_1414, t_1415, t_1416, t_1417, pc_x, msg0_763, msg0_764, \
                         msg1_763, msg1_764, msh_1063, msh_1064, msh_1065, msh_1066, \
                         msh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1413[k] = f_4 * msg0_763[k]
                    - f_5 * msg1_763[k]
                    + f_3 * pc_x[k] * msh_1063[k];

        t_1414[k] = f_4 * msg0_764[k]
                    - f_5 * msg1_764[k]
                    + f_3 * pc_x[k] * msh_1064[k];

        t_1415[k] = f_3 * pc_x[k] * msh_1065[k];

        t_1416[k] = f_3 * pc_x[k] * msh_1066[k];

        t_1417[k] = f_3 * pc_x[k] * msh_1067[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, t_1421, t_1422, pc_x, pc_y, pc_z, lsh_855, \
                         lsh_876, msg0_760, msg1_760, msh_1065, msh_1068, msh_1069, \
                         msh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_3 * pc_x[k] * msh_1068[k];

        t_1419[k] = f_3 * pc_x[k] * msh_1069[k];

        t_1420[k] = f_3 * pc_x[k] * msh_1070[k];

        t_1421[k] = f_14 * lsh_876[k]
                    + f_1 * msg0_760[k]
                    - f_2 * msg1_760[k]
                    + f_3 * pc_y[k] * msh_1065[k];

        t_1422[k] = f_20 * lsh_855[k]
                    + f_3 * pc_z[k] * msh_1065[k];
    }

#pragma omp simd aligned(t_1423, t_1424, t_1425, pc_y, lsh_878, lsh_879, lsh_880, msg0_762, \
                         msg0_763, msg0_764, msg1_762, msg1_763, msg1_764, msh_1067, msh_1068, \
                         msh_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1423[k] = f_14 * lsh_878[k]
                    + f_8 * msg0_762[k]
                    - f_9 * msg1_762[k]
                    + f_3 * pc_y[k] * msh_1067[k];

        t_1424[k] = f_14 * lsh_879[k]
                    + f_6 * msg0_763[k]
                    - f_7 * msg1_763[k]
                    + f_3 * pc_y[k] * msh_1068[k];

        t_1425[k] = f_14 * lsh_880[k]
                    + f_4 * msg0_764[k]
                    - f_5 * msg1_764[k]
                    + f_3 * pc_y[k] * msh_1069[k];
    }

#pragma omp simd aligned(t_1426, t_1427, t_1428, pc_x, pc_y, pc_z, lsh_860, lsh_881, msg0_764, \
                         msg0_765, msg1_764, msg1_765, msh_1070, \
                         msh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1426[k] = f_14 * lsh_881[k]
                    + f_3 * pc_y[k] * msh_1070[k];

        t_1427[k] = f_20 * lsh_860[k]
                    + f_1 * msg0_764[k]
                    - f_2 * msg1_764[k]
                    + f_3 * pc_z[k] * msh_1070[k];

        t_1428[k] = f_1 * msg0_765[k]
                    - f_2 * msg1_765[k]
                    + f_3 * pc_x[k] * msh_1071[k];
    }

#pragma omp simd aligned(t_1429, t_1430, t_1431, pc_x, msg0_766, msg0_767, msg0_768, msg1_766, \
                         msg1_767, msg1_768, msh_1072, msh_1073, \
                         msh_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1429[k] = f_16 * msg0_766[k]
                    - f_17 * msg1_766[k]
                    + f_3 * pc_x[k] * msh_1072[k];

        t_1430[k] = f_16 * msg0_767[k]
                    - f_17 * msg1_767[k]
                    + f_3 * pc_x[k] * msh_1073[k];

        t_1431[k] = f_8 * msg0_768[k]
                    - f_9 * msg1_768[k]
                    + f_3 * pc_x[k] * msh_1074[k];
    }

#pragma omp simd aligned(t_1432, t_1433, t_1434, pc_x, msg0_769, msg0_770, msg0_771, msg1_769, \
                         msg1_770, msg1_771, msh_1075, msh_1076, \
                         msh_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1432[k] = f_8 * msg0_769[k]
                    - f_9 * msg1_769[k]
                    + f_3 * pc_x[k] * msh_1075[k];

        t_1433[k] = f_8 * msg0_770[k]
                    - f_9 * msg1_770[k]
                    + f_3 * pc_x[k] * msh_1076[k];

        t_1434[k] = f_6 * msg0_771[k]
                    - f_7 * msg1_771[k]
                    + f_3 * pc_x[k] * msh_1077[k];
    }

#pragma omp simd aligned(t_1435, t_1436, t_1437, pc_x, msg0_772, msg0_773, msg0_774, msg1_772, \
                         msg1_773, msg1_774, msh_1078, msh_1079, \
                         msh_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1435[k] = f_6 * msg0_772[k]
                    - f_7 * msg1_772[k]
                    + f_3 * pc_x[k] * msh_1078[k];

        t_1436[k] = f_6 * msg0_773[k]
                    - f_7 * msg1_773[k]
                    + f_3 * pc_x[k] * msh_1079[k];

        t_1437[k] = f_6 * msg0_774[k]
                    - f_7 * msg1_774[k]
                    + f_3 * pc_x[k] * msh_1080[k];
    }

#pragma omp simd aligned(t_1438, t_1439, t_1440, pc_x, msg0_775, msg0_776, msg0_777, msg1_775, \
                         msg1_776, msg1_777, msh_1081, msh_1082, \
                         msh_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1438[k] = f_4 * msg0_775[k]
                    - f_5 * msg1_775[k]
                    + f_3 * pc_x[k] * msh_1081[k];

        t_1439[k] = f_4 * msg0_776[k]
                    - f_5 * msg1_776[k]
                    + f_3 * pc_x[k] * msh_1082[k];

        t_1440[k] = f_4 * msg0_777[k]
                    - f_5 * msg1_777[k]
                    + f_3 * pc_x[k] * msh_1083[k];
    }

#pragma omp simd aligned(t_1441, t_1442, t_1443, t_1444, t_1445, pc_x, msg0_778, msg0_779, \
                         msg1_778, msg1_779, msh_1084, msh_1085, msh_1086, msh_1087, \
                         msh_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1441[k] = f_4 * msg0_778[k]
                    - f_5 * msg1_778[k]
                    + f_3 * pc_x[k] * msh_1084[k];

        t_1442[k] = f_4 * msg0_779[k]
                    - f_5 * msg1_779[k]
                    + f_3 * pc_x[k] * msh_1085[k];

        t_1443[k] = f_3 * pc_x[k] * msh_1086[k];

        t_1444[k] = f_3 * pc_x[k] * msh_1087[k];

        t_1445[k] = f_3 * pc_x[k] * msh_1088[k];
    }

#pragma omp simd aligned(t_1446, t_1447, t_1448, t_1449, t_1450, pc_x, pc_y, pc_z, lsh_876, \
                         lsh_897, msg0_775, msg1_775, msh_1086, msh_1089, msh_1090, \
                         msh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1446[k] = f_3 * pc_x[k] * msh_1089[k];

        t_1447[k] = f_3 * pc_x[k] * msh_1090[k];

        t_1448[k] = f_3 * pc_x[k] * msh_1091[k];

        t_1449[k] = f_13 * lsh_897[k]
                    + f_1 * msg0_775[k]
                    - f_2 * msg1_775[k]
                    + f_3 * pc_y[k] * msh_1086[k];

        t_1450[k] = f_19 * lsh_876[k]
                    + f_3 * pc_z[k] * msh_1086[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_y, lsh_899, lsh_900, lsh_901, msg0_777, \
                         msg0_778, msg0_779, msg1_777, msg1_778, msg1_779, msh_1088, msh_1089, \
                         msh_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_13 * lsh_899[k]
                    + f_8 * msg0_777[k]
                    - f_9 * msg1_777[k]
                    + f_3 * pc_y[k] * msh_1088[k];

        t_1452[k] = f_13 * lsh_900[k]
                    + f_6 * msg0_778[k]
                    - f_7 * msg1_778[k]
                    + f_3 * pc_y[k] * msh_1089[k];

        t_1453[k] = f_13 * lsh_901[k]
                    + f_4 * msg0_779[k]
                    - f_5 * msg1_779[k]
                    + f_3 * pc_y[k] * msh_1090[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, pc_y, pc_z, lsh_881, lsh_902, msg0_779, \
                         msg0_780, msg1_779, msg1_780, msh_1091, \
                         msh_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_13 * lsh_902[k]
                    + f_3 * pc_y[k] * msh_1091[k];

        t_1455[k] = f_19 * lsh_881[k]
                    + f_1 * msg0_779[k]
                    - f_2 * msg1_779[k]
                    + f_3 * pc_z[k] * msh_1091[k];

        t_1456[k] = f_1 * msg0_780[k]
                    - f_2 * msg1_780[k]
                    + f_3 * pc_x[k] * msh_1092[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pc_x, msg0_781, msg0_782, msg0_783, msg1_781, \
                         msg1_782, msg1_783, msh_1093, msh_1094, \
                         msh_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_16 * msg0_781[k]
                    - f_17 * msg1_781[k]
                    + f_3 * pc_x[k] * msh_1093[k];

        t_1458[k] = f_16 * msg0_782[k]
                    - f_17 * msg1_782[k]
                    + f_3 * pc_x[k] * msh_1094[k];

        t_1459[k] = f_8 * msg0_783[k]
                    - f_9 * msg1_783[k]
                    + f_3 * pc_x[k] * msh_1095[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, pc_x, msg0_784, msg0_785, msg0_786, msg1_784, \
                         msg1_785, msg1_786, msh_1096, msh_1097, \
                         msh_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_8 * msg0_784[k]
                    - f_9 * msg1_784[k]
                    + f_3 * pc_x[k] * msh_1096[k];

        t_1461[k] = f_8 * msg0_785[k]
                    - f_9 * msg1_785[k]
                    + f_3 * pc_x[k] * msh_1097[k];

        t_1462[k] = f_6 * msg0_786[k]
                    - f_7 * msg1_786[k]
                    + f_3 * pc_x[k] * msh_1098[k];
    }

#pragma omp simd aligned(t_1463, t_1464, t_1465, pc_x, msg0_787, msg0_788, msg0_789, msg1_787, \
                         msg1_788, msg1_789, msh_1099, msh_1100, \
                         msh_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1463[k] = f_6 * msg0_787[k]
                    - f_7 * msg1_787[k]
                    + f_3 * pc_x[k] * msh_1099[k];

        t_1464[k] = f_6 * msg0_788[k]
                    - f_7 * msg1_788[k]
                    + f_3 * pc_x[k] * msh_1100[k];

        t_1465[k] = f_6 * msg0_789[k]
                    - f_7 * msg1_789[k]
                    + f_3 * pc_x[k] * msh_1101[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, pc_x, msg0_790, msg0_791, msg0_792, msg1_790, \
                         msg1_791, msg1_792, msh_1102, msh_1103, \
                         msh_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_4 * msg0_790[k]
                    - f_5 * msg1_790[k]
                    + f_3 * pc_x[k] * msh_1102[k];

        t_1467[k] = f_4 * msg0_791[k]
                    - f_5 * msg1_791[k]
                    + f_3 * pc_x[k] * msh_1103[k];

        t_1468[k] = f_4 * msg0_792[k]
                    - f_5 * msg1_792[k]
                    + f_3 * pc_x[k] * msh_1104[k];
    }

#pragma omp simd aligned(t_1469, t_1470, t_1471, t_1472, t_1473, pc_x, msg0_793, msg0_794, \
                         msg1_793, msg1_794, msh_1105, msh_1106, msh_1107, msh_1108, \
                         msh_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1469[k] = f_4 * msg0_793[k]
                    - f_5 * msg1_793[k]
                    + f_3 * pc_x[k] * msh_1105[k];

        t_1470[k] = f_4 * msg0_794[k]
                    - f_5 * msg1_794[k]
                    + f_3 * pc_x[k] * msh_1106[k];

        t_1471[k] = f_3 * pc_x[k] * msh_1107[k];

        t_1472[k] = f_3 * pc_x[k] * msh_1108[k];

        t_1473[k] = f_3 * pc_x[k] * msh_1109[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, t_1478, pc_x, pc_y, pc_z, lsh_897, \
                         lsh_918, msg0_790, msg1_790, msh_1107, msh_1110, msh_1111, \
                         msh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_3 * pc_x[k] * msh_1110[k];

        t_1475[k] = f_3 * pc_x[k] * msh_1111[k];

        t_1476[k] = f_3 * pc_x[k] * msh_1112[k];

        t_1477[k] = f_12 * lsh_918[k]
                    + f_1 * msg0_790[k]
                    - f_2 * msg1_790[k]
                    + f_3 * pc_y[k] * msh_1107[k];

        t_1478[k] = f_18 * lsh_897[k]
                    + f_3 * pc_z[k] * msh_1107[k];
    }

#pragma omp simd aligned(t_1479, t_1480, t_1481, pc_y, lsh_920, lsh_921, lsh_922, msg0_792, \
                         msg0_793, msg0_794, msg1_792, msg1_793, msg1_794, msh_1109, msh_1110, \
                         msh_1111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1479[k] = f_12 * lsh_920[k]
                    + f_8 * msg0_792[k]
                    - f_9 * msg1_792[k]
                    + f_3 * pc_y[k] * msh_1109[k];

        t_1480[k] = f_12 * lsh_921[k]
                    + f_6 * msg0_793[k]
                    - f_7 * msg1_793[k]
                    + f_3 * pc_y[k] * msh_1110[k];

        t_1481[k] = f_12 * lsh_922[k]
                    + f_4 * msg0_794[k]
                    - f_5 * msg1_794[k]
                    + f_3 * pc_y[k] * msh_1111[k];
    }

#pragma omp simd aligned(t_1482, t_1483, t_1484, pa_y, pc_y, pc_z, lsi0_1232, lsh_902, \
                         lsh_923, lsi1_1232, msg0_794, msg1_794, \
                         msh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1482[k] = f_12 * lsh_923[k]
                    + f_3 * pc_y[k] * msh_1112[k];

        t_1483[k] = f_18 * lsh_902[k]
                    + f_1 * msg0_794[k]
                    - f_2 * msg1_794[k]
                    + f_3 * pc_z[k] * msh_1112[k];

        t_1484[k] = pa_y[k] * lsi0_1232[k]
                    - f_10 * pc_y[k] * lsi1_1232[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, pa_y, pc_x, pc_y, lsi0_1234, lsi1_1234, \
                         msg0_796, msg0_798, msg1_796, msg1_798, msh_1114, \
                         msh_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = f_16 * msg0_796[k]
                    - f_17 * msg1_796[k]
                    + f_3 * pc_x[k] * msh_1114[k];

        t_1486[k] = pa_y[k] * lsi0_1234[k]
                    - f_10 * pc_y[k] * lsi1_1234[k];

        t_1487[k] = f_8 * msg0_798[k]
                    - f_9 * msg1_798[k]
                    + f_3 * pc_x[k] * msh_1116[k];
    }

#pragma omp simd aligned(t_1488, t_1489, t_1490, pa_y, pc_x, pc_y, lsi0_1237, lsi1_1237, \
                         msg0_799, msg0_801, msg1_799, msg1_801, msh_1117, \
                         msh_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1488[k] = f_8 * msg0_799[k]
                    - f_9 * msg1_799[k]
                    + f_3 * pc_x[k] * msh_1117[k];

        t_1489[k] = pa_y[k] * lsi0_1237[k]
                    - f_10 * pc_y[k] * lsi1_1237[k];

        t_1490[k] = f_6 * msg0_801[k]
                    - f_7 * msg1_801[k]
                    + f_3 * pc_x[k] * msh_1119[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, pa_y, pc_x, pc_y, lsi0_1241, lsi1_1241, \
                         msg0_802, msg0_803, msg1_802, msg1_803, msh_1120, \
                         msh_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = f_6 * msg0_802[k]
                    - f_7 * msg1_802[k]
                    + f_3 * pc_x[k] * msh_1120[k];

        t_1492[k] = f_6 * msg0_803[k]
                    - f_7 * msg1_803[k]
                    + f_3 * pc_x[k] * msh_1121[k];

        t_1493[k] = pa_y[k] * lsi0_1241[k]
                    - f_10 * pc_y[k] * lsi1_1241[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, pc_x, msg0_805, msg0_806, msg0_807, msg1_805, \
                         msg1_806, msg1_807, msh_1123, msh_1124, \
                         msh_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = f_4 * msg0_805[k]
                    - f_5 * msg1_805[k]
                    + f_3 * pc_x[k] * msh_1123[k];

        t_1495[k] = f_4 * msg0_806[k]
                    - f_5 * msg1_806[k]
                    + f_3 * pc_x[k] * msh_1124[k];

        t_1496[k] = f_4 * msg0_807[k]
                    - f_5 * msg1_807[k]
                    + f_3 * pc_x[k] * msh_1125[k];
    }

#pragma omp simd aligned(t_1497, t_1498, t_1499, t_1500, t_1501, pa_y, pc_x, pc_y, lsi0_1246, \
                         lsi1_1246, msg0_808, msg1_808, msh_1126, msh_1128, msh_1129, \
                         msh_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1497[k] = f_4 * msg0_808[k]
                    - f_5 * msg1_808[k]
                    + f_3 * pc_x[k] * msh_1126[k];

        t_1498[k] = pa_y[k] * lsi0_1246[k]
                    - f_10 * pc_y[k] * lsi1_1246[k];

        t_1499[k] = f_3 * pc_x[k] * msh_1128[k];

        t_1500[k] = f_3 * pc_x[k] * msh_1129[k];

        t_1501[k] = f_3 * pc_x[k] * msh_1130[k];
    }

#pragma omp simd aligned(t_1502, t_1503, t_1504, t_1505, pa_y, pc_x, pc_y, lsi0_1253, lsh_939, \
                         lsi1_1253, msh_1131, msh_1132, msh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = f_3 * pc_x[k] * msh_1131[k];

        t_1503[k] = f_3 * pc_x[k] * msh_1132[k];

        t_1504[k] = f_3 * pc_x[k] * msh_1133[k];

        t_1505[k] = pa_y[k] * lsi0_1253[k]
                    + f_19 * lsh_939[k]
                    - f_10 * pc_y[k] * lsi1_1253[k];
    }

#pragma omp simd aligned(t_1506, t_1507, t_1508, pa_y, pc_y, pc_z, lsi0_1255, lsi0_1256, \
                         lsh_918, lsh_941, lsh_942, lsi1_1255, lsi1_1256, \
                         msh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1506[k] = f_15 * lsh_918[k]
                    + f_3 * pc_z[k] * msh_1128[k];

        t_1507[k] = pa_y[k] * lsi0_1255[k]
                    + f_14 * lsh_941[k]
                    - f_10 * pc_y[k] * lsi1_1255[k];

        t_1508[k] = pa_y[k] * lsi0_1256[k]
                    + f_13 * lsh_942[k]
                    - f_10 * pc_y[k] * lsi1_1256[k];
    }

#pragma omp simd aligned(t_1509, t_1510, t_1511, pa_y, pc_y, lsi0_1257, lsi0_1259, lsh_943, \
                         lsh_944, lsi1_1257, lsi1_1259, msh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1509[k] = pa_y[k] * lsi0_1257[k]
                    + f_12 * lsh_943[k]
                    - f_10 * pc_y[k] * lsi1_1257[k];

        t_1510[k] = f_11 * lsh_944[k]
                    + f_3 * pc_y[k] * msh_1133[k];

        t_1511[k] = pa_y[k] * lsi0_1259[k]
                    - f_10 * pc_y[k] * lsi1_1259[k];
    }

#pragma omp simd aligned(t_1512, t_1513, t_1514, t_1515, t_1516, pc_x, pc_y, msg0_810, \
                         msg0_812, msg0_813, msg1_810, msg1_812, msg1_813, msh_1134, msh_1136, \
                         msh_1137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1512[k] = f_1 * msg0_810[k]
                    - f_2 * msg1_810[k]
                    + f_3 * pc_x[k] * msh_1134[k];

        t_1513[k] = f_3 * pc_y[k] * msh_1134[k];

        t_1514[k] = f_16 * msg0_812[k]
                    - f_17 * msg1_812[k]
                    + f_3 * pc_x[k] * msh_1136[k];

        t_1515[k] = f_8 * msg0_813[k]
                    - f_9 * msg1_813[k]
                    + f_3 * pc_x[k] * msh_1137[k];

        t_1516[k] = f_3 * pc_y[k] * msh_1136[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, t_1520, pc_x, pc_y, msg0_815, msg0_816, \
                         msg0_817, msg1_815, msg1_816, msg1_817, msh_1139, msh_1140, \
                         msh_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = f_8 * msg0_815[k]
                    - f_9 * msg1_815[k]
                    + f_3 * pc_x[k] * msh_1139[k];

        t_1518[k] = f_6 * msg0_816[k]
                    - f_7 * msg1_816[k]
                    + f_3 * pc_x[k] * msh_1140[k];

        t_1519[k] = f_6 * msg0_817[k]
                    - f_7 * msg1_817[k]
                    + f_3 * pc_x[k] * msh_1141[k];

        t_1520[k] = f_3 * pc_y[k] * msh_1139[k];
    }
}

static auto
compute_prim_msi_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t lsh, const size_t msg0,
                                                           const size_t msg1, const size_t msh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsh_944 = buffer.data(lsh + 944);

    const auto *msg0_819 = buffer.data(msg0 + 819);
    const auto *msg0_820 = buffer.data(msg0 + 820);
    const auto *msg0_821 = buffer.data(msg0 + 821);
    const auto *msg0_822 = buffer.data(msg0 + 822);
    const auto *msg0_823 = buffer.data(msg0 + 823);
    const auto *msg0_824 = buffer.data(msg0 + 824);

    const auto *msg1_819 = buffer.data(msg1 + 819);
    const auto *msg1_820 = buffer.data(msg1 + 820);
    const auto *msg1_821 = buffer.data(msg1 + 821);
    const auto *msg1_822 = buffer.data(msg1 + 822);
    const auto *msg1_823 = buffer.data(msg1 + 823);
    const auto *msg1_824 = buffer.data(msg1 + 824);

    const auto *msh_1143 = buffer.data(msh + 1143);
    const auto *msh_1144 = buffer.data(msh + 1144);
    const auto *msh_1145 = buffer.data(msh + 1145);
    const auto *msh_1146 = buffer.data(msh + 1146);
    const auto *msh_1148 = buffer.data(msh + 1148);
    const auto *msh_1149 = buffer.data(msh + 1149);
    const auto *msh_1150 = buffer.data(msh + 1150);
    const auto *msh_1151 = buffer.data(msh + 1151);
    const auto *msh_1152 = buffer.data(msh + 1152);
    const auto *msh_1153 = buffer.data(msh + 1153);
    const auto *msh_1154 = buffer.data(msh + 1154);

#pragma omp simd aligned(t_1521, t_1522, t_1523, pc_x, msg0_819, msg0_820, msg0_821, msg1_819, \
                         msg1_820, msg1_821, msh_1143, msh_1144, \
                         msh_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = f_6 * msg0_819[k]
                    - f_7 * msg1_819[k]
                    + f_3 * pc_x[k] * msh_1143[k];

        t_1522[k] = f_4 * msg0_820[k]
                    - f_5 * msg1_820[k]
                    + f_3 * pc_x[k] * msh_1144[k];

        t_1523[k] = f_4 * msg0_821[k]
                    - f_5 * msg1_821[k]
                    + f_3 * pc_x[k] * msh_1145[k];
    }

#pragma omp simd aligned(t_1524, t_1525, t_1526, t_1527, t_1528, pc_x, pc_y, msg0_822, \
                         msg0_824, msg1_822, msg1_824, msh_1143, msh_1146, msh_1148, msh_1149, \
                         msh_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1524[k] = f_4 * msg0_822[k]
                    - f_5 * msg1_822[k]
                    + f_3 * pc_x[k] * msh_1146[k];

        t_1525[k] = f_3 * pc_y[k] * msh_1143[k];

        t_1526[k] = f_4 * msg0_824[k]
                    - f_5 * msg1_824[k]
                    + f_3 * pc_x[k] * msh_1148[k];

        t_1527[k] = f_3 * pc_x[k] * msh_1149[k];

        t_1528[k] = f_3 * pc_x[k] * msh_1150[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, t_1532, t_1533, pc_x, pc_y, msg0_820, \
                         msg1_820, msh_1149, msh_1151, msh_1152, msh_1153, \
                         msh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_3 * pc_x[k] * msh_1151[k];

        t_1530[k] = f_3 * pc_x[k] * msh_1152[k];

        t_1531[k] = f_3 * pc_x[k] * msh_1153[k];

        t_1532[k] = f_3 * pc_x[k] * msh_1154[k];

        t_1533[k] = f_1 * msg0_820[k]
                    - f_2 * msg1_820[k]
                    + f_3 * pc_y[k] * msh_1149[k];
    }

#pragma omp simd aligned(t_1534, t_1535, t_1536, pc_y, msg0_821, msg0_822, msg0_823, msg1_821, \
                         msg1_822, msg1_823, msh_1150, msh_1151, \
                         msh_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1534[k] = f_16 * msg0_821[k]
                    - f_17 * msg1_821[k]
                    + f_3 * pc_y[k] * msh_1150[k];

        t_1535[k] = f_8 * msg0_822[k]
                    - f_9 * msg1_822[k]
                    + f_3 * pc_y[k] * msh_1151[k];

        t_1536[k] = f_6 * msg0_823[k]
                    - f_7 * msg1_823[k]
                    + f_3 * pc_y[k] * msh_1152[k];
    }

#pragma omp simd aligned(t_1537, t_1538, t_1539, pc_y, pc_z, lsh_944, msg0_824, msg1_824, \
                         msh_1153, msh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1537[k] = f_4 * msg0_824[k]
                    - f_5 * msg1_824[k]
                    + f_3 * pc_y[k] * msh_1153[k];

        t_1538[k] = f_3 * pc_y[k] * msh_1154[k];

        t_1539[k] = f_0 * lsh_944[k]
                    + f_1 * msg0_824[k]
                    - f_2 * msg1_824[k]
                    + f_3 * pc_z[k] * msh_1154[k];
    }
}

auto
compute_prim_msi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t lsi0, const size_t lsh,
                                                   const size_t lsi1, const size_t msg0,
                                                   const size_t msg1, const size_t msh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_msi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msg0, msg1, msh, ncols,
                                                              gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, lsi0, lsh,
                                                              lsi1, msh, ncols, gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, lsi0,
                                                               lsh, lsi1, msg0, msg1, msh,
                                                               ncols, gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, lsi0,
                                                               lsh, lsi1, msg0, msg1, msh,
                                                               ncols, gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece12(buffer, target, pa, pc, lsi0,
                                                               lsh, lsi1, msg0, msg1, msh,
                                                               ncols, gamma, p, q);

    compute_prim_msi_three_center_electron_repulsion_0_piece13(buffer, target, pc, lsh, msg0,
                                                               msg1, msh, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
