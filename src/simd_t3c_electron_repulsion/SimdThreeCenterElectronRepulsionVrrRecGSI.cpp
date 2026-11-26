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


#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gsi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t gsg0, const size_t gsg1,
                                                          const size_t gsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);

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

    const auto *fsi0_0 = buffer.data(fsi0 + 0);
    const auto *fsi0_3 = buffer.data(fsi0 + 3);
    const auto *fsi0_5 = buffer.data(fsi0 + 5);
    const auto *fsi0_6 = buffer.data(fsi0 + 6);
    const auto *fsi0_9 = buffer.data(fsi0 + 9);
    const auto *fsi0_10 = buffer.data(fsi0 + 10);
    const auto *fsi0_14 = buffer.data(fsi0 + 14);
    const auto *fsi0_21 = buffer.data(fsi0 + 21);
    const auto *fsi0_27 = buffer.data(fsi0 + 27);
    const auto *fsi0_31 = buffer.data(fsi0 + 31);
    const auto *fsi0_34 = buffer.data(fsi0 + 34);
    const auto *fsi0_38 = buffer.data(fsi0 + 38);
    const auto *fsi0_56 = buffer.data(fsi0 + 56);
    const auto *fsi0_61 = buffer.data(fsi0 + 61);
    const auto *fsi0_65 = buffer.data(fsi0 + 65);
    const auto *fsi0_68 = buffer.data(fsi0 + 68);
    const auto *fsi0_70 = buffer.data(fsi0 + 70);

    const auto *fsh_0 = buffer.data(fsh + 0);
    const auto *fsh_1 = buffer.data(fsh + 1);
    const auto *fsh_2 = buffer.data(fsh + 2);
    const auto *fsh_3 = buffer.data(fsh + 3);
    const auto *fsh_5 = buffer.data(fsh + 5);
    const auto *fsh_6 = buffer.data(fsh + 6);
    const auto *fsh_9 = buffer.data(fsh + 9);
    const auto *fsh_15 = buffer.data(fsh + 15);
    const auto *fsh_17 = buffer.data(fsh + 17);
    const auto *fsh_18 = buffer.data(fsh + 18);
    const auto *fsh_20 = buffer.data(fsh + 20);
    const auto *fsh_21 = buffer.data(fsh + 21);
    const auto *fsh_24 = buffer.data(fsh + 24);
    const auto *fsh_26 = buffer.data(fsh + 26);
    const auto *fsh_27 = buffer.data(fsh + 27);
    const auto *fsh_30 = buffer.data(fsh + 30);
    const auto *fsh_36 = buffer.data(fsh + 36);
    const auto *fsh_38 = buffer.data(fsh + 38);
    const auto *fsh_39 = buffer.data(fsh + 39);
    const auto *fsh_40 = buffer.data(fsh + 40);
    const auto *fsh_41 = buffer.data(fsh + 41);
    const auto *fsh_42 = buffer.data(fsh + 42);
    const auto *fsh_44 = buffer.data(fsh + 44);
    const auto *fsh_47 = buffer.data(fsh + 47);
    const auto *fsh_50 = buffer.data(fsh + 50);
    const auto *fsh_51 = buffer.data(fsh + 51);
    const auto *fsh_57 = buffer.data(fsh + 57);
    const auto *fsh_58 = buffer.data(fsh + 58);
    const auto *fsh_59 = buffer.data(fsh + 59);
    const auto *fsh_60 = buffer.data(fsh + 60);
    const auto *fsh_62 = buffer.data(fsh + 62);
    const auto *fsh_63 = buffer.data(fsh + 63);
    const auto *fsh_66 = buffer.data(fsh + 66);
    const auto *fsh_69 = buffer.data(fsh + 69);
    const auto *fsh_73 = buffer.data(fsh + 73);
    const auto *fsh_78 = buffer.data(fsh + 78);
    const auto *fsh_80 = buffer.data(fsh + 80);
    const auto *fsh_81 = buffer.data(fsh + 81);
    const auto *fsh_82 = buffer.data(fsh + 82);
    const auto *fsh_83 = buffer.data(fsh + 83);
    const auto *fsh_99 = buffer.data(fsh + 99);

    const auto *fsi1_0 = buffer.data(fsi1 + 0);
    const auto *fsi1_3 = buffer.data(fsi1 + 3);
    const auto *fsi1_5 = buffer.data(fsi1 + 5);
    const auto *fsi1_6 = buffer.data(fsi1 + 6);
    const auto *fsi1_9 = buffer.data(fsi1 + 9);
    const auto *fsi1_10 = buffer.data(fsi1 + 10);
    const auto *fsi1_14 = buffer.data(fsi1 + 14);
    const auto *fsi1_21 = buffer.data(fsi1 + 21);
    const auto *fsi1_27 = buffer.data(fsi1 + 27);
    const auto *fsi1_31 = buffer.data(fsi1 + 31);
    const auto *fsi1_34 = buffer.data(fsi1 + 34);
    const auto *fsi1_38 = buffer.data(fsi1 + 38);
    const auto *fsi1_56 = buffer.data(fsi1 + 56);
    const auto *fsi1_61 = buffer.data(fsi1 + 61);
    const auto *fsi1_65 = buffer.data(fsi1 + 65);
    const auto *fsi1_68 = buffer.data(fsi1 + 68);
    const auto *fsi1_70 = buffer.data(fsi1 + 70);

    const auto *gsg0_0 = buffer.data(gsg0 + 0);
    const auto *gsg0_1 = buffer.data(gsg0 + 1);
    const auto *gsg0_2 = buffer.data(gsg0 + 2);
    const auto *gsg0_3 = buffer.data(gsg0 + 3);
    const auto *gsg0_5 = buffer.data(gsg0 + 5);
    const auto *gsg0_10 = buffer.data(gsg0 + 10);
    const auto *gsg0_12 = buffer.data(gsg0 + 12);
    const auto *gsg0_13 = buffer.data(gsg0 + 13);
    const auto *gsg0_14 = buffer.data(gsg0 + 14);
    const auto *gsg0_18 = buffer.data(gsg0 + 18);
    const auto *gsg0_25 = buffer.data(gsg0 + 25);
    const auto *gsg0_26 = buffer.data(gsg0 + 26);
    const auto *gsg0_27 = buffer.data(gsg0 + 27);
    const auto *gsg0_32 = buffer.data(gsg0 + 32);
    const auto *gsg0_34 = buffer.data(gsg0 + 34);
    const auto *gsg0_35 = buffer.data(gsg0 + 35);
    const auto *gsg0_41 = buffer.data(gsg0 + 41);
    const auto *gsg0_42 = buffer.data(gsg0 + 42);
    const auto *gsg0_43 = buffer.data(gsg0 + 43);
    const auto *gsg0_44 = buffer.data(gsg0 + 44);
    const auto *gsg0_45 = buffer.data(gsg0 + 45);
    const auto *gsg0_47 = buffer.data(gsg0 + 47);
    const auto *gsg0_48 = buffer.data(gsg0 + 48);
    const auto *gsg0_50 = buffer.data(gsg0 + 50);
    const auto *gsg0_51 = buffer.data(gsg0 + 51);
    const auto *gsg0_55 = buffer.data(gsg0 + 55);
    const auto *gsg0_56 = buffer.data(gsg0 + 56);
    const auto *gsg0_57 = buffer.data(gsg0 + 57);
    const auto *gsg0_59 = buffer.data(gsg0 + 59);

    const auto *gsg1_0 = buffer.data(gsg1 + 0);
    const auto *gsg1_1 = buffer.data(gsg1 + 1);
    const auto *gsg1_2 = buffer.data(gsg1 + 2);
    const auto *gsg1_3 = buffer.data(gsg1 + 3);
    const auto *gsg1_5 = buffer.data(gsg1 + 5);
    const auto *gsg1_10 = buffer.data(gsg1 + 10);
    const auto *gsg1_12 = buffer.data(gsg1 + 12);
    const auto *gsg1_13 = buffer.data(gsg1 + 13);
    const auto *gsg1_14 = buffer.data(gsg1 + 14);
    const auto *gsg1_18 = buffer.data(gsg1 + 18);
    const auto *gsg1_25 = buffer.data(gsg1 + 25);
    const auto *gsg1_26 = buffer.data(gsg1 + 26);
    const auto *gsg1_27 = buffer.data(gsg1 + 27);
    const auto *gsg1_32 = buffer.data(gsg1 + 32);
    const auto *gsg1_34 = buffer.data(gsg1 + 34);
    const auto *gsg1_35 = buffer.data(gsg1 + 35);
    const auto *gsg1_41 = buffer.data(gsg1 + 41);
    const auto *gsg1_42 = buffer.data(gsg1 + 42);
    const auto *gsg1_43 = buffer.data(gsg1 + 43);
    const auto *gsg1_44 = buffer.data(gsg1 + 44);
    const auto *gsg1_45 = buffer.data(gsg1 + 45);
    const auto *gsg1_47 = buffer.data(gsg1 + 47);
    const auto *gsg1_48 = buffer.data(gsg1 + 48);
    const auto *gsg1_50 = buffer.data(gsg1 + 50);
    const auto *gsg1_51 = buffer.data(gsg1 + 51);
    const auto *gsg1_55 = buffer.data(gsg1 + 55);
    const auto *gsg1_56 = buffer.data(gsg1 + 56);
    const auto *gsg1_57 = buffer.data(gsg1 + 57);
    const auto *gsg1_59 = buffer.data(gsg1 + 59);

    const auto *gsh_0 = buffer.data(gsh + 0);
    const auto *gsh_1 = buffer.data(gsh + 1);
    const auto *gsh_2 = buffer.data(gsh + 2);
    const auto *gsh_3 = buffer.data(gsh + 3);
    const auto *gsh_5 = buffer.data(gsh + 5);
    const auto *gsh_6 = buffer.data(gsh + 6);
    const auto *gsh_8 = buffer.data(gsh + 8);
    const auto *gsh_9 = buffer.data(gsh + 9);
    const auto *gsh_10 = buffer.data(gsh + 10);
    const auto *gsh_14 = buffer.data(gsh + 14);
    const auto *gsh_15 = buffer.data(gsh + 15);
    const auto *gsh_17 = buffer.data(gsh + 17);
    const auto *gsh_18 = buffer.data(gsh + 18);
    const auto *gsh_19 = buffer.data(gsh + 19);
    const auto *gsh_20 = buffer.data(gsh + 20);
    const auto *gsh_21 = buffer.data(gsh + 21);
    const auto *gsh_22 = buffer.data(gsh + 22);
    const auto *gsh_24 = buffer.data(gsh + 24);
    const auto *gsh_26 = buffer.data(gsh + 26);
    const auto *gsh_27 = buffer.data(gsh + 27);
    const auto *gsh_28 = buffer.data(gsh + 28);
    const auto *gsh_30 = buffer.data(gsh + 30);
    const auto *gsh_31 = buffer.data(gsh + 31);
    const auto *gsh_36 = buffer.data(gsh + 36);
    const auto *gsh_37 = buffer.data(gsh + 37);
    const auto *gsh_38 = buffer.data(gsh + 38);
    const auto *gsh_39 = buffer.data(gsh + 39);
    const auto *gsh_40 = buffer.data(gsh + 40);
    const auto *gsh_41 = buffer.data(gsh + 41);
    const auto *gsh_42 = buffer.data(gsh + 42);
    const auto *gsh_44 = buffer.data(gsh + 44);
    const auto *gsh_46 = buffer.data(gsh + 46);
    const auto *gsh_47 = buffer.data(gsh + 47);
    const auto *gsh_49 = buffer.data(gsh + 49);
    const auto *gsh_50 = buffer.data(gsh + 50);
    const auto *gsh_51 = buffer.data(gsh + 51);
    const auto *gsh_56 = buffer.data(gsh + 56);
    const auto *gsh_57 = buffer.data(gsh + 57);
    const auto *gsh_58 = buffer.data(gsh + 58);
    const auto *gsh_59 = buffer.data(gsh + 59);
    const auto *gsh_60 = buffer.data(gsh + 60);
    const auto *gsh_61 = buffer.data(gsh + 61);
    const auto *gsh_62 = buffer.data(gsh + 62);
    const auto *gsh_63 = buffer.data(gsh + 63);
    const auto *gsh_64 = buffer.data(gsh + 64);
    const auto *gsh_65 = buffer.data(gsh + 65);
    const auto *gsh_66 = buffer.data(gsh + 66);
    const auto *gsh_68 = buffer.data(gsh + 68);
    const auto *gsh_69 = buffer.data(gsh + 69);
    const auto *gsh_70 = buffer.data(gsh + 70);
    const auto *gsh_72 = buffer.data(gsh + 72);
    const auto *gsh_73 = buffer.data(gsh + 73);
    const auto *gsh_78 = buffer.data(gsh + 78);
    const auto *gsh_79 = buffer.data(gsh + 79);
    const auto *gsh_80 = buffer.data(gsh + 80);
    const auto *gsh_81 = buffer.data(gsh + 81);
    const auto *gsh_82 = buffer.data(gsh + 82);
    const auto *gsh_83 = buffer.data(gsh + 83);
    const auto *gsh_84 = buffer.data(gsh + 84);
    const auto *gsh_86 = buffer.data(gsh + 86);
    const auto *gsh_87 = buffer.data(gsh + 87);
    const auto *gsh_89 = buffer.data(gsh + 89);
    const auto *gsh_90 = buffer.data(gsh + 90);
    const auto *gsh_93 = buffer.data(gsh + 93);
    const auto *gsh_99 = buffer.data(gsh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fsh_0, gsg0_0, \
                         gsg1_0, gsh_0, gsh_1, gsh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fsh_0[k]
                 + f_1 * gsg0_0[k]
                 - f_2 * gsg1_0[k]
                 + f_3 * pc_x[k] * gsh_0[k];

        t_1[k] = f_3 * pc_y[k] * gsh_0[k];

        t_2[k] = f_3 * pc_z[k] * gsh_0[k];

        t_3[k] = f_4 * gsg0_0[k]
                 - f_5 * gsg1_0[k]
                 + f_3 * pc_y[k] * gsh_1[k];

        t_4[k] = f_3 * pc_y[k] * gsh_2[k];

        t_5[k] = f_4 * gsg0_0[k]
                 - f_5 * gsg1_0[k]
                 + f_3 * pc_z[k] * gsh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, gsg0_1, gsg0_2, gsg0_3, gsg1_1, \
                         gsg1_2, gsg1_3, gsh_3, gsh_5, gsh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * gsg0_1[k]
                 - f_7 * gsg1_1[k]
                 + f_3 * pc_y[k] * gsh_3[k];

        t_7[k] = f_3 * pc_z[k] * gsh_3[k];

        t_8[k] = f_3 * pc_y[k] * gsh_5[k];

        t_9[k] = f_6 * gsg0_2[k]
                 - f_7 * gsg1_2[k]
                 + f_3 * pc_z[k] * gsh_5[k];

        t_10[k] = f_8 * gsg0_3[k]
                  - f_9 * gsg1_3[k]
                  + f_3 * pc_y[k] * gsh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, fsh_15, gsg0_5, \
                         gsg1_5, gsh_6, gsh_8, gsh_9, gsh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * gsh_6[k];

        t_12[k] = f_4 * gsg0_5[k]
                  - f_5 * gsg1_5[k]
                  + f_3 * pc_y[k] * gsh_8[k];

        t_13[k] = f_3 * pc_y[k] * gsh_9[k];

        t_14[k] = f_8 * gsg0_5[k]
                  - f_9 * gsg1_5[k]
                  + f_3 * pc_z[k] * gsh_9[k];

        t_15[k] = f_0 * fsh_15[k]
                  + f_3 * pc_x[k] * gsh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, fsh_17, fsh_18, \
                         fsh_20, gsh_10, gsh_14, gsh_17, gsh_18, \
                         gsh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * gsh_10[k];

        t_17[k] = f_0 * fsh_17[k]
                  + f_3 * pc_x[k] * gsh_17[k];

        t_18[k] = f_0 * fsh_18[k]
                  + f_3 * pc_x[k] * gsh_18[k];

        t_19[k] = f_3 * pc_y[k] * gsh_14[k];

        t_20[k] = f_0 * fsh_20[k]
                  + f_3 * pc_x[k] * gsh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, gsg0_10, gsg0_12, gsg0_13, \
                         gsg1_10, gsg1_12, gsg1_13, gsh_15, gsh_17, \
                         gsh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * gsg0_10[k]
                  - f_2 * gsg1_10[k]
                  + f_3 * pc_y[k] * gsh_15[k];

        t_22[k] = f_3 * pc_z[k] * gsh_15[k];

        t_23[k] = f_8 * gsg0_12[k]
                  - f_9 * gsg1_12[k]
                  + f_3 * pc_y[k] * gsh_17[k];

        t_24[k] = f_6 * gsg0_13[k]
                  - f_7 * gsg1_13[k]
                  + f_3 * pc_y[k] * gsh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, fsi0_0, fsh_0, \
                         fsi1_0, gsg0_14, gsg1_14, gsh_19, gsh_20, \
                         gsh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * gsg0_14[k]
                  - f_5 * gsg1_14[k]
                  + f_3 * pc_y[k] * gsh_19[k];

        t_26[k] = f_3 * pc_y[k] * gsh_20[k];

        t_27[k] = f_1 * gsg0_14[k]
                  - f_2 * gsg1_14[k]
                  + f_3 * pc_z[k] * gsh_20[k];

        t_28[k] = pa_y[k] * fsi0_0[k]
                  - f_10 * pc_y[k] * fsi1_0[k];

        t_29[k] = f_11 * fsh_0[k]
                  + f_3 * pc_y[k] * gsh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, fsi0_3, fsi0_5, fsh_1, \
                         fsi1_3, fsi1_5, gsh_21, gsh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * gsh_21[k];

        t_31[k] = pa_y[k] * fsi0_3[k]
                  + f_12 * fsh_1[k]
                  - f_10 * pc_y[k] * fsi1_3[k];

        t_32[k] = f_3 * pc_z[k] * gsh_22[k];

        t_33[k] = pa_y[k] * fsi0_5[k]
                  - f_10 * pc_y[k] * fsi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, fsi0_6, fsi0_9, fsh_3, \
                         fsh_5, fsi1_6, fsi1_9, gsh_24, gsh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * fsi0_6[k]
                  + f_13 * fsh_3[k]
                  - f_10 * pc_y[k] * fsi1_6[k];

        t_35[k] = f_3 * pc_z[k] * gsh_24[k];

        t_36[k] = f_11 * fsh_5[k]
                  + f_3 * pc_y[k] * gsh_26[k];

        t_37[k] = pa_y[k] * fsi0_9[k]
                  - f_10 * pc_y[k] * fsi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, fsi0_10, fsh_6, fsh_9, \
                         fsi1_10, gsg0_18, gsg1_18, gsh_27, gsh_28, \
                         gsh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * fsi0_10[k]
                  + f_0 * fsh_6[k]
                  - f_10 * pc_y[k] * fsi1_10[k];

        t_39[k] = f_3 * pc_z[k] * gsh_27[k];

        t_40[k] = f_4 * gsg0_18[k]
                  - f_5 * gsg1_18[k]
                  + f_3 * pc_z[k] * gsh_28[k];

        t_41[k] = f_11 * fsh_9[k]
                  + f_3 * pc_y[k] * gsh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, fsi0_14, fsh_36, \
                         fsh_38, fsi1_14, gsh_31, gsh_36, gsh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * fsi0_14[k]
                  - f_10 * pc_y[k] * fsi1_14[k];

        t_43[k] = f_13 * fsh_36[k]
                  + f_3 * pc_x[k] * gsh_36[k];

        t_44[k] = f_3 * pc_z[k] * gsh_31[k];

        t_45[k] = f_13 * fsh_38[k]
                  + f_3 * pc_x[k] * gsh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, fsh_15, fsh_39, fsh_40, fsh_41, \
                         gsg0_25, gsg1_25, gsh_36, gsh_39, gsh_40, \
                         gsh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_13 * fsh_39[k]
                  + f_3 * pc_x[k] * gsh_39[k];

        t_47[k] = f_13 * fsh_40[k]
                  + f_3 * pc_x[k] * gsh_40[k];

        t_48[k] = f_13 * fsh_41[k]
                  + f_3 * pc_x[k] * gsh_41[k];

        t_49[k] = f_11 * fsh_15[k]
                  + f_1 * gsg0_25[k]
                  - f_2 * gsg1_25[k]
                  + f_3 * pc_y[k] * gsh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, gsg0_25, gsg0_26, gsg0_27, gsg1_25, \
                         gsg1_26, gsg1_27, gsh_36, gsh_37, gsh_38, \
                         gsh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * gsh_36[k];

        t_51[k] = f_4 * gsg0_25[k]
                  - f_5 * gsg1_25[k]
                  + f_3 * pc_z[k] * gsh_37[k];

        t_52[k] = f_6 * gsg0_26[k]
                  - f_7 * gsg1_26[k]
                  + f_3 * pc_z[k] * gsh_38[k];

        t_53[k] = f_8 * gsg0_27[k]
                  - f_9 * gsg1_27[k]
                  + f_3 * pc_z[k] * gsh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, fsi0_0, fsi0_27, \
                         fsh_20, fsi1_0, fsi1_27, gsh_41, gsh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * fsh_20[k]
                  + f_3 * pc_y[k] * gsh_41[k];

        t_55[k] = pa_y[k] * fsi0_27[k]
                  - f_10 * pc_y[k] * fsi1_27[k];

        t_56[k] = pa_z[k] * fsi0_0[k]
                  - f_10 * pc_z[k] * fsi1_0[k];

        t_57[k] = f_3 * pc_y[k] * gsh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, fsi0_3, fsi0_5, fsh_0, \
                         fsh_2, fsi1_3, fsi1_5, gsh_42, gsh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * fsh_0[k]
                  + f_3 * pc_z[k] * gsh_42[k];

        t_59[k] = pa_z[k] * fsi0_3[k]
                  - f_10 * pc_z[k] * fsi1_3[k];

        t_60[k] = f_3 * pc_y[k] * gsh_44[k];

        t_61[k] = pa_z[k] * fsi0_5[k]
                  + f_12 * fsh_2[k]
                  - f_10 * pc_z[k] * fsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, fsi0_6, fsi0_9, fsh_5, \
                         fsi1_6, fsi1_9, gsg0_32, gsg1_32, gsh_46, \
                         gsh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * fsi0_6[k]
                  - f_10 * pc_z[k] * fsi1_6[k];

        t_63[k] = f_4 * gsg0_32[k]
                  - f_5 * gsg1_32[k]
                  + f_3 * pc_y[k] * gsh_46[k];

        t_64[k] = f_3 * pc_y[k] * gsh_47[k];

        t_65[k] = pa_z[k] * fsi0_9[k]
                  + f_13 * fsh_5[k]
                  - f_10 * pc_z[k] * fsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, fsi0_10, fsi1_10, gsg0_34, \
                         gsg0_35, gsg1_34, gsg1_35, gsh_49, gsh_50, \
                         gsh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * fsi0_10[k]
                  - f_10 * pc_z[k] * fsi1_10[k];

        t_67[k] = f_6 * gsg0_34[k]
                  - f_7 * gsg1_34[k]
                  + f_3 * pc_y[k] * gsh_49[k];

        t_68[k] = f_4 * gsg0_35[k]
                  - f_5 * gsg1_35[k]
                  + f_3 * pc_y[k] * gsh_50[k];

        t_69[k] = f_3 * pc_y[k] * gsh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, fsi0_14, fsh_9, fsh_57, \
                         fsh_58, fsh_59, fsi1_14, gsh_57, gsh_58, \
                         gsh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * fsi0_14[k]
                  + f_0 * fsh_9[k]
                  - f_10 * pc_z[k] * fsi1_14[k];

        t_71[k] = f_13 * fsh_57[k]
                  + f_3 * pc_x[k] * gsh_57[k];

        t_72[k] = f_13 * fsh_58[k]
                  + f_3 * pc_x[k] * gsh_58[k];

        t_73[k] = f_13 * fsh_59[k]
                  + f_3 * pc_x[k] * gsh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, fsi0_21, fsh_60, \
                         fsh_62, fsi1_21, gsh_56, gsh_60, gsh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * fsh_60[k]
                  + f_3 * pc_x[k] * gsh_60[k];

        t_75[k] = f_3 * pc_y[k] * gsh_56[k];

        t_76[k] = f_13 * fsh_62[k]
                  + f_3 * pc_x[k] * gsh_62[k];

        t_77[k] = pa_z[k] * fsi0_21[k]
                  - f_10 * pc_z[k] * fsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, gsg0_41, gsg0_42, gsg0_43, gsg1_41, gsg1_42, \
                         gsg1_43, gsh_58, gsh_59, gsh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_14 * gsg0_41[k]
                  - f_15 * gsg1_41[k]
                  + f_3 * pc_y[k] * gsh_58[k];

        t_79[k] = f_8 * gsg0_42[k]
                  - f_9 * gsg1_42[k]
                  + f_3 * pc_y[k] * gsh_59[k];

        t_80[k] = f_6 * gsg0_43[k]
                  - f_7 * gsg1_43[k]
                  + f_3 * pc_y[k] * gsh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, fsh_20, fsh_63, gsg0_44, \
                         gsg0_45, gsg1_44, gsg1_45, gsh_61, gsh_62, \
                         gsh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * gsg0_44[k]
                  - f_5 * gsg1_44[k]
                  + f_3 * pc_y[k] * gsh_61[k];

        t_82[k] = f_3 * pc_y[k] * gsh_62[k];

        t_83[k] = f_11 * fsh_20[k]
                  + f_1 * gsg0_44[k]
                  - f_2 * gsg1_44[k]
                  + f_3 * pc_z[k] * gsh_62[k];

        t_84[k] = f_12 * fsh_63[k]
                  + f_1 * gsg0_45[k]
                  - f_2 * gsg1_45[k]
                  + f_3 * pc_x[k] * gsh_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, fsh_21, fsh_66, gsg0_48, \
                         gsg1_48, gsh_63, gsh_64, gsh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * fsh_21[k]
                  + f_3 * pc_y[k] * gsh_63[k];

        t_86[k] = f_3 * pc_z[k] * gsh_63[k];

        t_87[k] = f_12 * fsh_66[k]
                  + f_8 * gsg0_48[k]
                  - f_9 * gsg1_48[k]
                  + f_3 * pc_x[k] * gsh_66[k];

        t_88[k] = f_3 * pc_z[k] * gsh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, fsh_69, gsg0_45, gsg0_51, gsg1_45, \
                         gsg1_51, gsh_65, gsh_66, gsh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * gsg0_45[k]
                  - f_5 * gsg1_45[k]
                  + f_3 * pc_z[k] * gsh_65[k];

        t_90[k] = f_12 * fsh_69[k]
                  + f_6 * gsg0_51[k]
                  - f_7 * gsg1_51[k]
                  + f_3 * pc_x[k] * gsh_69[k];

        t_91[k] = f_3 * pc_z[k] * gsh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, fsh_26, fsh_73, gsg0_47, \
                         gsg0_55, gsg1_47, gsg1_55, gsh_68, gsh_69, \
                         gsh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * fsh_26[k]
                  + f_3 * pc_y[k] * gsh_68[k];

        t_93[k] = f_6 * gsg0_47[k]
                  - f_7 * gsg1_47[k]
                  + f_3 * pc_z[k] * gsh_68[k];

        t_94[k] = f_12 * fsh_73[k]
                  + f_4 * gsg0_55[k]
                  - f_5 * gsg1_55[k]
                  + f_3 * pc_x[k] * gsh_73[k];

        t_95[k] = f_3 * pc_z[k] * gsh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, fsh_30, fsh_78, gsg0_48, \
                         gsg0_50, gsg1_48, gsg1_50, gsh_70, gsh_72, \
                         gsh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * gsg0_48[k]
                  - f_5 * gsg1_48[k]
                  + f_3 * pc_z[k] * gsh_70[k];

        t_97[k] = f_12 * fsh_30[k]
                  + f_3 * pc_y[k] * gsh_72[k];

        t_98[k] = f_8 * gsg0_50[k]
                  - f_9 * gsg1_50[k]
                  + f_3 * pc_z[k] * gsh_72[k];

        t_99[k] = f_12 * fsh_78[k]
                  + f_3 * pc_x[k] * gsh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, fsh_80, fsh_81, \
                         fsh_82, fsh_83, gsh_73, gsh_80, gsh_81, gsh_82, \
                         gsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * gsh_73[k];

        t_101[k] = f_12 * fsh_80[k]
                   + f_3 * pc_x[k] * gsh_80[k];

        t_102[k] = f_12 * fsh_81[k]
                   + f_3 * pc_x[k] * gsh_81[k];

        t_103[k] = f_12 * fsh_82[k]
                   + f_3 * pc_x[k] * gsh_82[k];

        t_104[k] = f_12 * fsh_83[k]
                   + f_3 * pc_x[k] * gsh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, fsh_36, gsg0_55, gsg0_56, \
                         gsg1_55, gsg1_56, gsh_78, gsh_79, gsh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * fsh_36[k]
                   + f_1 * gsg0_55[k]
                   - f_2 * gsg1_55[k]
                   + f_3 * pc_y[k] * gsh_78[k];

        t_106[k] = f_3 * pc_z[k] * gsh_78[k];

        t_107[k] = f_4 * gsg0_55[k]
                   - f_5 * gsg1_55[k]
                   + f_3 * pc_z[k] * gsh_79[k];

        t_108[k] = f_6 * gsg0_56[k]
                   - f_7 * gsg1_56[k]
                   + f_3 * pc_z[k] * gsh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, fsi0_56, fsh_41, \
                         fsi1_56, gsg0_57, gsg0_59, gsg1_57, gsg1_59, gsh_81, \
                         gsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * gsg0_57[k]
                   - f_9 * gsg1_57[k]
                   + f_3 * pc_z[k] * gsh_81[k];

        t_110[k] = f_12 * fsh_41[k]
                   + f_3 * pc_y[k] * gsh_83[k];

        t_111[k] = f_1 * gsg0_59[k]
                   - f_2 * gsg1_59[k]
                   + f_3 * pc_z[k] * gsh_83[k];

        t_112[k] = pa_y[k] * fsi0_56[k]
                   - f_10 * pc_y[k] * fsi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, fsi0_31, fsh_21, \
                         fsh_42, fsh_44, fsi1_31, gsh_84, gsh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * fsh_42[k]
                   + f_3 * pc_y[k] * gsh_84[k];

        t_114[k] = f_11 * fsh_21[k]
                   + f_3 * pc_z[k] * gsh_84[k];

        t_115[k] = pa_z[k] * fsi0_31[k]
                   - f_10 * pc_z[k] * fsi1_31[k];

        t_116[k] = f_11 * fsh_44[k]
                   + f_3 * pc_y[k] * gsh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, fsi0_34, fsi0_61, \
                         fsh_24, fsh_47, fsi1_34, fsi1_61, gsh_87, \
                         gsh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * fsi0_61[k]
                   - f_10 * pc_y[k] * fsi1_61[k];

        t_118[k] = pa_z[k] * fsi0_34[k]
                   - f_10 * pc_z[k] * fsi1_34[k];

        t_119[k] = f_11 * fsh_24[k]
                   + f_3 * pc_z[k] * gsh_87[k];

        t_120[k] = f_11 * fsh_47[k]
                   + f_3 * pc_y[k] * gsh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, fsi0_38, fsi0_65, \
                         fsh_27, fsi1_38, fsi1_65, gsh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * fsi0_65[k]
                   - f_10 * pc_y[k] * fsi1_65[k];

        t_122[k] = pa_z[k] * fsi0_38[k]
                   - f_10 * pc_z[k] * fsi1_38[k];

        t_123[k] = f_11 * fsh_27[k]
                   + f_3 * pc_z[k] * gsh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, fsi0_68, fsi0_70, \
                         fsh_50, fsh_51, fsh_99, fsi1_68, fsi1_70, gsh_93, \
                         gsh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * fsi0_68[k]
                   + f_12 * fsh_50[k]
                   - f_10 * pc_y[k] * fsi1_68[k];

        t_125[k] = f_11 * fsh_51[k]
                   + f_3 * pc_y[k] * gsh_93[k];

        t_126[k] = pa_y[k] * fsi0_70[k]
                   - f_10 * pc_y[k] * fsi1_70[k];

        t_127[k] = f_12 * fsh_99[k]
                   + f_3 * pc_x[k] * gsh_99[k];
    }
}

static auto
compute_prim_gsi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t gsg0, const size_t gsg1,
                                                          const size_t gsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);
    const auto f_16 = 3.0 / q;

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
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsi0_49 = buffer.data(fsi0 + 49);
    const auto *fsi0_83 = buffer.data(fsi0 + 83);
    const auto *fsi0_84 = buffer.data(fsi0 + 84);
    const auto *fsi0_87 = buffer.data(fsi0 + 87);
    const auto *fsi0_90 = buffer.data(fsi0 + 90);
    const auto *fsi0_94 = buffer.data(fsi0 + 94);
    const auto *fsi0_140 = buffer.data(fsi0 + 140);
    const auto *fsi0_145 = buffer.data(fsi0 + 145);
    const auto *fsi0_149 = buffer.data(fsi0 + 149);
    const auto *fsi0_154 = buffer.data(fsi0 + 154);
    const auto *fsi0_168 = buffer.data(fsi0 + 168);
    const auto *fsi0_171 = buffer.data(fsi0 + 171);
    const auto *fsi0_174 = buffer.data(fsi0 + 174);
    const auto *fsi0_178 = buffer.data(fsi0 + 178);
    const auto *fsi0_189 = buffer.data(fsi0 + 189);
    const auto *fsi0_191 = buffer.data(fsi0 + 191);
    const auto *fsi0_192 = buffer.data(fsi0 + 192);
    const auto *fsi0_193 = buffer.data(fsi0 + 193);
    const auto *fsi0_195 = buffer.data(fsi0 + 195);
    const auto *fsi0_201 = buffer.data(fsi0 + 201);
    const auto *fsi0_205 = buffer.data(fsi0 + 205);
    const auto *fsi0_208 = buffer.data(fsi0 + 208);
    const auto *fsi0_210 = buffer.data(fsi0 + 210);
    const auto *fsi0_217 = buffer.data(fsi0 + 217);
    const auto *fsi0_219 = buffer.data(fsi0 + 219);
    const auto *fsi0_220 = buffer.data(fsi0 + 220);
    const auto *fsi0_221 = buffer.data(fsi0 + 221);
    const auto *fsi0_223 = buffer.data(fsi0 + 223);
    const auto *fsi0_227 = buffer.data(fsi0 + 227);
    const auto *fsi0_230 = buffer.data(fsi0 + 230);
    const auto *fsi0_234 = buffer.data(fsi0 + 234);
    const auto *fsi0_236 = buffer.data(fsi0 + 236);
    const auto *fsi0_245 = buffer.data(fsi0 + 245);
    const auto *fsi0_247 = buffer.data(fsi0 + 247);
    const auto *fsi0_248 = buffer.data(fsi0 + 248);
    const auto *fsi0_249 = buffer.data(fsi0 + 249);

    const auto *fsh_36 = buffer.data(fsh + 36);
    const auto *fsh_42 = buffer.data(fsh + 42);
    const auto *fsh_59 = buffer.data(fsh + 59);
    const auto *fsh_60 = buffer.data(fsh + 60);
    const auto *fsh_61 = buffer.data(fsh + 61);
    const auto *fsh_62 = buffer.data(fsh + 62);
    const auto *fsh_63 = buffer.data(fsh + 63);
    const auto *fsh_66 = buffer.data(fsh + 66);
    const auto *fsh_68 = buffer.data(fsh + 68);
    const auto *fsh_69 = buffer.data(fsh + 69);
    const auto *fsh_72 = buffer.data(fsh + 72);
    const auto *fsh_78 = buffer.data(fsh + 78);
    const auto *fsh_83 = buffer.data(fsh + 83);
    const auto *fsh_84 = buffer.data(fsh + 84);
    const auto *fsh_86 = buffer.data(fsh + 86);
    const auto *fsh_87 = buffer.data(fsh + 87);
    const auto *fsh_89 = buffer.data(fsh + 89);
    const auto *fsh_90 = buffer.data(fsh + 90);
    const auto *fsh_93 = buffer.data(fsh + 93);
    const auto *fsh_99 = buffer.data(fsh + 99);
    const auto *fsh_100 = buffer.data(fsh + 100);
    const auto *fsh_101 = buffer.data(fsh + 101);
    const auto *fsh_102 = buffer.data(fsh + 102);
    const auto *fsh_103 = buffer.data(fsh + 103);
    const auto *fsh_104 = buffer.data(fsh + 104);
    const auto *fsh_105 = buffer.data(fsh + 105);
    const auto *fsh_107 = buffer.data(fsh + 107);
    const auto *fsh_110 = buffer.data(fsh + 110);
    const auto *fsh_114 = buffer.data(fsh + 114);
    const auto *fsh_119 = buffer.data(fsh + 119);
    const auto *fsh_120 = buffer.data(fsh + 120);
    const auto *fsh_121 = buffer.data(fsh + 121);
    const auto *fsh_122 = buffer.data(fsh + 122);
    const auto *fsh_123 = buffer.data(fsh + 123);
    const auto *fsh_125 = buffer.data(fsh + 125);
    const auto *fsh_126 = buffer.data(fsh + 126);
    const auto *fsh_129 = buffer.data(fsh + 129);
    const auto *fsh_132 = buffer.data(fsh + 132);
    const auto *fsh_136 = buffer.data(fsh + 136);
    const auto *fsh_141 = buffer.data(fsh + 141);
    const auto *fsh_143 = buffer.data(fsh + 143);
    const auto *fsh_144 = buffer.data(fsh + 144);
    const auto *fsh_145 = buffer.data(fsh + 145);
    const auto *fsh_146 = buffer.data(fsh + 146);
    const auto *fsh_152 = buffer.data(fsh + 152);
    const auto *fsh_156 = buffer.data(fsh + 156);
    const auto *fsh_159 = buffer.data(fsh + 159);
    const auto *fsh_161 = buffer.data(fsh + 161);
    const auto *fsh_162 = buffer.data(fsh + 162);
    const auto *fsh_163 = buffer.data(fsh + 163);
    const auto *fsh_164 = buffer.data(fsh + 164);
    const auto *fsh_165 = buffer.data(fsh + 165);
    const auto *fsh_166 = buffer.data(fsh + 166);
    const auto *fsh_167 = buffer.data(fsh + 167);
    const auto *fsh_171 = buffer.data(fsh + 171);
    const auto *fsh_174 = buffer.data(fsh + 174);
    const auto *fsh_178 = buffer.data(fsh + 178);
    const auto *fsh_180 = buffer.data(fsh + 180);
    const auto *fsh_183 = buffer.data(fsh + 183);
    const auto *fsh_184 = buffer.data(fsh + 184);
    const auto *fsh_185 = buffer.data(fsh + 185);
    const auto *fsh_186 = buffer.data(fsh + 186);
    const auto *fsh_187 = buffer.data(fsh + 187);
    const auto *fsh_188 = buffer.data(fsh + 188);

    const auto *fsi1_49 = buffer.data(fsi1 + 49);
    const auto *fsi1_83 = buffer.data(fsi1 + 83);
    const auto *fsi1_84 = buffer.data(fsi1 + 84);
    const auto *fsi1_87 = buffer.data(fsi1 + 87);
    const auto *fsi1_90 = buffer.data(fsi1 + 90);
    const auto *fsi1_94 = buffer.data(fsi1 + 94);
    const auto *fsi1_140 = buffer.data(fsi1 + 140);
    const auto *fsi1_145 = buffer.data(fsi1 + 145);
    const auto *fsi1_149 = buffer.data(fsi1 + 149);
    const auto *fsi1_154 = buffer.data(fsi1 + 154);
    const auto *fsi1_168 = buffer.data(fsi1 + 168);
    const auto *fsi1_171 = buffer.data(fsi1 + 171);
    const auto *fsi1_174 = buffer.data(fsi1 + 174);
    const auto *fsi1_178 = buffer.data(fsi1 + 178);
    const auto *fsi1_189 = buffer.data(fsi1 + 189);
    const auto *fsi1_191 = buffer.data(fsi1 + 191);
    const auto *fsi1_192 = buffer.data(fsi1 + 192);
    const auto *fsi1_193 = buffer.data(fsi1 + 193);
    const auto *fsi1_195 = buffer.data(fsi1 + 195);
    const auto *fsi1_201 = buffer.data(fsi1 + 201);
    const auto *fsi1_205 = buffer.data(fsi1 + 205);
    const auto *fsi1_208 = buffer.data(fsi1 + 208);
    const auto *fsi1_210 = buffer.data(fsi1 + 210);
    const auto *fsi1_217 = buffer.data(fsi1 + 217);
    const auto *fsi1_219 = buffer.data(fsi1 + 219);
    const auto *fsi1_220 = buffer.data(fsi1 + 220);
    const auto *fsi1_221 = buffer.data(fsi1 + 221);
    const auto *fsi1_223 = buffer.data(fsi1 + 223);
    const auto *fsi1_227 = buffer.data(fsi1 + 227);
    const auto *fsi1_230 = buffer.data(fsi1 + 230);
    const auto *fsi1_234 = buffer.data(fsi1 + 234);
    const auto *fsi1_236 = buffer.data(fsi1 + 236);
    const auto *fsi1_245 = buffer.data(fsi1 + 245);
    const auto *fsi1_247 = buffer.data(fsi1 + 247);
    const auto *fsi1_248 = buffer.data(fsi1 + 248);
    const auto *fsi1_249 = buffer.data(fsi1 + 249);

    const auto *gsg0_72 = buffer.data(gsg0 + 72);
    const auto *gsg0_73 = buffer.data(gsg0 + 73);
    const auto *gsg0_74 = buffer.data(gsg0 + 74);
    const auto *gsg0_75 = buffer.data(gsg0 + 75);
    const auto *gsg0_76 = buffer.data(gsg0 + 76);
    const auto *gsg0_77 = buffer.data(gsg0 + 77);
    const auto *gsg0_78 = buffer.data(gsg0 + 78);
    const auto *gsg0_79 = buffer.data(gsg0 + 79);
    const auto *gsg0_80 = buffer.data(gsg0 + 80);
    const auto *gsg0_84 = buffer.data(gsg0 + 84);
    const auto *gsg0_85 = buffer.data(gsg0 + 85);
    const auto *gsg0_86 = buffer.data(gsg0 + 86);
    const auto *gsg0_87 = buffer.data(gsg0 + 87);
    const auto *gsg0_88 = buffer.data(gsg0 + 88);
    const auto *gsg0_89 = buffer.data(gsg0 + 89);
    const auto *gsg0_90 = buffer.data(gsg0 + 90);
    const auto *gsg0_92 = buffer.data(gsg0 + 92);
    const auto *gsg0_93 = buffer.data(gsg0 + 93);
    const auto *gsg0_95 = buffer.data(gsg0 + 95);

    const auto *gsg1_72 = buffer.data(gsg1 + 72);
    const auto *gsg1_73 = buffer.data(gsg1 + 73);
    const auto *gsg1_74 = buffer.data(gsg1 + 74);
    const auto *gsg1_75 = buffer.data(gsg1 + 75);
    const auto *gsg1_76 = buffer.data(gsg1 + 76);
    const auto *gsg1_77 = buffer.data(gsg1 + 77);
    const auto *gsg1_78 = buffer.data(gsg1 + 78);
    const auto *gsg1_79 = buffer.data(gsg1 + 79);
    const auto *gsg1_80 = buffer.data(gsg1 + 80);
    const auto *gsg1_84 = buffer.data(gsg1 + 84);
    const auto *gsg1_85 = buffer.data(gsg1 + 85);
    const auto *gsg1_86 = buffer.data(gsg1 + 86);
    const auto *gsg1_87 = buffer.data(gsg1 + 87);
    const auto *gsg1_88 = buffer.data(gsg1 + 88);
    const auto *gsg1_89 = buffer.data(gsg1 + 89);
    const auto *gsg1_90 = buffer.data(gsg1 + 90);
    const auto *gsg1_92 = buffer.data(gsg1 + 92);
    const auto *gsg1_93 = buffer.data(gsg1 + 93);
    const auto *gsg1_95 = buffer.data(gsg1 + 95);

    const auto *gsh_99 = buffer.data(gsh + 99);
    const auto *gsh_100 = buffer.data(gsh + 100);
    const auto *gsh_101 = buffer.data(gsh + 101);
    const auto *gsh_102 = buffer.data(gsh + 102);
    const auto *gsh_103 = buffer.data(gsh + 103);
    const auto *gsh_104 = buffer.data(gsh + 104);
    const auto *gsh_105 = buffer.data(gsh + 105);
    const auto *gsh_106 = buffer.data(gsh + 106);
    const auto *gsh_107 = buffer.data(gsh + 107);
    const auto *gsh_108 = buffer.data(gsh + 108);
    const auto *gsh_109 = buffer.data(gsh + 109);
    const auto *gsh_110 = buffer.data(gsh + 110);
    const auto *gsh_111 = buffer.data(gsh + 111);
    const auto *gsh_112 = buffer.data(gsh + 112);
    const auto *gsh_113 = buffer.data(gsh + 113);
    const auto *gsh_114 = buffer.data(gsh + 114);
    const auto *gsh_119 = buffer.data(gsh + 119);
    const auto *gsh_120 = buffer.data(gsh + 120);
    const auto *gsh_121 = buffer.data(gsh + 121);
    const auto *gsh_122 = buffer.data(gsh + 122);
    const auto *gsh_123 = buffer.data(gsh + 123);
    const auto *gsh_124 = buffer.data(gsh + 124);
    const auto *gsh_125 = buffer.data(gsh + 125);
    const auto *gsh_126 = buffer.data(gsh + 126);
    const auto *gsh_127 = buffer.data(gsh + 127);
    const auto *gsh_128 = buffer.data(gsh + 128);
    const auto *gsh_129 = buffer.data(gsh + 129);
    const auto *gsh_131 = buffer.data(gsh + 131);
    const auto *gsh_132 = buffer.data(gsh + 132);
    const auto *gsh_133 = buffer.data(gsh + 133);
    const auto *gsh_135 = buffer.data(gsh + 135);
    const auto *gsh_136 = buffer.data(gsh + 136);
    const auto *gsh_141 = buffer.data(gsh + 141);
    const auto *gsh_143 = buffer.data(gsh + 143);
    const auto *gsh_144 = buffer.data(gsh + 144);
    const auto *gsh_145 = buffer.data(gsh + 145);
    const auto *gsh_146 = buffer.data(gsh + 146);
    const auto *gsh_147 = buffer.data(gsh + 147);
    const auto *gsh_149 = buffer.data(gsh + 149);
    const auto *gsh_150 = buffer.data(gsh + 150);
    const auto *gsh_152 = buffer.data(gsh + 152);
    const auto *gsh_153 = buffer.data(gsh + 153);
    const auto *gsh_156 = buffer.data(gsh + 156);
    const auto *gsh_162 = buffer.data(gsh + 162);
    const auto *gsh_163 = buffer.data(gsh + 163);
    const auto *gsh_164 = buffer.data(gsh + 164);
    const auto *gsh_165 = buffer.data(gsh + 165);
    const auto *gsh_166 = buffer.data(gsh + 166);
    const auto *gsh_167 = buffer.data(gsh + 167);
    const auto *gsh_168 = buffer.data(gsh + 168);
    const auto *gsh_170 = buffer.data(gsh + 170);
    const auto *gsh_171 = buffer.data(gsh + 171);
    const auto *gsh_173 = buffer.data(gsh + 173);
    const auto *gsh_174 = buffer.data(gsh + 174);
    const auto *gsh_177 = buffer.data(gsh + 177);
    const auto *gsh_183 = buffer.data(gsh + 183);
    const auto *gsh_184 = buffer.data(gsh + 184);
    const auto *gsh_185 = buffer.data(gsh + 185);
    const auto *gsh_186 = buffer.data(gsh + 186);
    const auto *gsh_187 = buffer.data(gsh + 187);
    const auto *gsh_188 = buffer.data(gsh + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, fsh_100, fsh_101, fsh_102, \
                         fsh_103, fsh_104, gsh_100, gsh_101, gsh_102, gsh_103, \
                         gsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_12 * fsh_100[k]
                   + f_3 * pc_x[k] * gsh_100[k];

        t_129[k] = f_12 * fsh_101[k]
                   + f_3 * pc_x[k] * gsh_101[k];

        t_130[k] = f_12 * fsh_102[k]
                   + f_3 * pc_x[k] * gsh_102[k];

        t_131[k] = f_12 * fsh_103[k]
                   + f_3 * pc_x[k] * gsh_103[k];

        t_132[k] = f_12 * fsh_104[k]
                   + f_3 * pc_x[k] * gsh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, fsi0_49, fsh_36, fsh_59, \
                         fsi1_49, gsg0_72, gsg1_72, gsh_99, gsh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * fsi0_49[k]
                   - f_10 * pc_z[k] * fsi1_49[k];

        t_134[k] = f_11 * fsh_36[k]
                   + f_3 * pc_z[k] * gsh_99[k];

        t_135[k] = f_11 * fsh_59[k]
                   + f_8 * gsg0_72[k]
                   - f_9 * gsg1_72[k]
                   + f_3 * pc_y[k] * gsh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, fsh_60, fsh_61, fsh_62, gsg0_73, gsg0_74, \
                         gsg1_73, gsg1_74, gsh_102, gsh_103, gsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * fsh_60[k]
                   + f_6 * gsg0_73[k]
                   - f_7 * gsg1_73[k]
                   + f_3 * pc_y[k] * gsh_102[k];

        t_137[k] = f_11 * fsh_61[k]
                   + f_4 * gsg0_74[k]
                   - f_5 * gsg1_74[k]
                   + f_3 * pc_y[k] * gsh_103[k];

        t_138[k] = f_11 * fsh_62[k]
                   + f_3 * pc_y[k] * gsh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, fsi0_83, fsh_42, \
                         fsh_105, fsi1_83, gsg0_75, gsg1_75, gsh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * fsi0_83[k]
                   - f_10 * pc_y[k] * fsi1_83[k];

        t_140[k] = f_12 * fsh_105[k]
                   + f_1 * gsg0_75[k]
                   - f_2 * gsg1_75[k]
                   + f_3 * pc_x[k] * gsh_105[k];

        t_141[k] = f_3 * pc_y[k] * gsh_105[k];

        t_142[k] = f_12 * fsh_42[k]
                   + f_3 * pc_z[k] * gsh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, fsh_110, gsg0_75, gsg0_80, gsg1_75, \
                         gsg1_80, gsh_106, gsh_107, gsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * gsg0_75[k]
                   - f_5 * gsg1_75[k]
                   + f_3 * pc_y[k] * gsh_106[k];

        t_144[k] = f_3 * pc_y[k] * gsh_107[k];

        t_145[k] = f_12 * fsh_110[k]
                   + f_8 * gsg0_80[k]
                   - f_9 * gsg1_80[k]
                   + f_3 * pc_x[k] * gsh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, gsg0_76, gsg0_77, gsg1_76, gsg1_77, \
                         gsh_108, gsh_109, gsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * gsg0_76[k]
                   - f_7 * gsg1_76[k]
                   + f_3 * pc_y[k] * gsh_108[k];

        t_147[k] = f_4 * gsg0_77[k]
                   - f_5 * gsg1_77[k]
                   + f_3 * pc_y[k] * gsh_109[k];

        t_148[k] = f_3 * pc_y[k] * gsh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, fsh_114, gsg0_78, gsg0_79, gsg0_84, \
                         gsg1_78, gsg1_79, gsg1_84, gsh_111, gsh_112, \
                         gsh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_12 * fsh_114[k]
                   + f_6 * gsg0_84[k]
                   - f_7 * gsg1_84[k]
                   + f_3 * pc_x[k] * gsh_114[k];

        t_150[k] = f_8 * gsg0_78[k]
                   - f_9 * gsg1_78[k]
                   + f_3 * pc_y[k] * gsh_111[k];

        t_151[k] = f_6 * gsg0_79[k]
                   - f_7 * gsg1_79[k]
                   + f_3 * pc_y[k] * gsh_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, fsh_119, fsh_120, gsg0_80, \
                         gsg0_89, gsg1_80, gsg1_89, gsh_113, gsh_114, gsh_119, \
                         gsh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * gsg0_80[k]
                   - f_5 * gsg1_80[k]
                   + f_3 * pc_y[k] * gsh_113[k];

        t_153[k] = f_3 * pc_y[k] * gsh_114[k];

        t_154[k] = f_12 * fsh_119[k]
                   + f_4 * gsg0_89[k]
                   - f_5 * gsg1_89[k]
                   + f_3 * pc_x[k] * gsh_119[k];

        t_155[k] = f_12 * fsh_120[k]
                   + f_3 * pc_x[k] * gsh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, fsh_121, fsh_122, \
                         fsh_123, fsh_125, gsh_119, gsh_121, gsh_122, gsh_123, \
                         gsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_12 * fsh_121[k]
                   + f_3 * pc_x[k] * gsh_121[k];

        t_157[k] = f_12 * fsh_122[k]
                   + f_3 * pc_x[k] * gsh_122[k];

        t_158[k] = f_12 * fsh_123[k]
                   + f_3 * pc_x[k] * gsh_123[k];

        t_159[k] = f_3 * pc_y[k] * gsh_119[k];

        t_160[k] = f_12 * fsh_125[k]
                   + f_3 * pc_x[k] * gsh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, gsg0_85, gsg0_86, gsg0_87, gsg1_85, \
                         gsg1_86, gsg1_87, gsh_120, gsh_121, gsh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * gsg0_85[k]
                   - f_2 * gsg1_85[k]
                   + f_3 * pc_y[k] * gsh_120[k];

        t_162[k] = f_14 * gsg0_86[k]
                   - f_15 * gsg1_86[k]
                   + f_3 * pc_y[k] * gsh_121[k];

        t_163[k] = f_8 * gsg0_87[k]
                   - f_9 * gsg1_87[k]
                   + f_3 * pc_y[k] * gsh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, fsh_62, gsg0_88, gsg0_89, \
                         gsg1_88, gsg1_89, gsh_123, gsh_124, gsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * gsg0_88[k]
                   - f_7 * gsg1_88[k]
                   + f_3 * pc_y[k] * gsh_123[k];

        t_165[k] = f_4 * gsg0_89[k]
                   - f_5 * gsg1_89[k]
                   + f_3 * pc_y[k] * gsh_124[k];

        t_166[k] = f_3 * pc_y[k] * gsh_125[k];

        t_167[k] = f_12 * fsh_62[k]
                   + f_1 * gsg0_89[k]
                   - f_2 * gsg1_89[k]
                   + f_3 * pc_z[k] * gsh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_x, pc_x, pc_y, pc_z, fsi0_168, \
                         fsi0_171, fsh_63, fsh_126, fsh_129, fsi1_168, fsi1_171, \
                         gsh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_x[k] * fsi0_168[k]
                   + f_16 * fsh_126[k]
                   - f_10 * pc_x[k] * fsi1_168[k];

        t_169[k] = f_13 * fsh_63[k]
                   + f_3 * pc_y[k] * gsh_126[k];

        t_170[k] = f_3 * pc_z[k] * gsh_126[k];

        t_171[k] = pa_x[k] * fsi0_171[k]
                   + f_0 * fsh_129[k]
                   - f_10 * pc_x[k] * fsi1_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_x, pc_x, pc_z, fsi0_174, fsh_132, \
                         fsi1_174, gsg0_90, gsg1_90, gsh_127, gsh_128, \
                         gsh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * gsh_127[k];

        t_173[k] = f_4 * gsg0_90[k]
                   - f_5 * gsg1_90[k]
                   + f_3 * pc_z[k] * gsh_128[k];

        t_174[k] = pa_x[k] * fsi0_174[k]
                   + f_13 * fsh_132[k]
                   - f_10 * pc_x[k] * fsi1_174[k];

        t_175[k] = f_3 * pc_z[k] * gsh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_x, pc_x, pc_y, pc_z, fsi0_178, fsh_68, \
                         fsh_136, fsi1_178, gsg0_92, gsg1_92, gsh_131, \
                         gsh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * fsh_68[k]
                   + f_3 * pc_y[k] * gsh_131[k];

        t_177[k] = f_6 * gsg0_92[k]
                   - f_7 * gsg1_92[k]
                   + f_3 * pc_z[k] * gsh_131[k];

        t_178[k] = pa_x[k] * fsi0_178[k]
                   + f_12 * fsh_136[k]
                   - f_10 * pc_x[k] * fsi1_178[k];

        t_179[k] = f_3 * pc_z[k] * gsh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, fsh_72, fsh_141, \
                         gsg0_93, gsg0_95, gsg1_93, gsg1_95, gsh_133, gsh_135, \
                         gsh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * gsg0_93[k]
                   - f_5 * gsg1_93[k]
                   + f_3 * pc_z[k] * gsh_133[k];

        t_181[k] = f_13 * fsh_72[k]
                   + f_3 * pc_y[k] * gsh_135[k];

        t_182[k] = f_8 * gsg0_95[k]
                   - f_9 * gsg1_95[k]
                   + f_3 * pc_z[k] * gsh_135[k];

        t_183[k] = f_11 * fsh_141[k]
                   + f_3 * pc_x[k] * gsh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, fsh_143, fsh_144, \
                         fsh_145, fsh_146, gsh_136, gsh_143, gsh_144, gsh_145, \
                         gsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * gsh_136[k];

        t_185[k] = f_11 * fsh_143[k]
                   + f_3 * pc_x[k] * gsh_143[k];

        t_186[k] = f_11 * fsh_144[k]
                   + f_3 * pc_x[k] * gsh_144[k];

        t_187[k] = f_11 * fsh_145[k]
                   + f_3 * pc_x[k] * gsh_145[k];

        t_188[k] = f_11 * fsh_146[k]
                   + f_3 * pc_x[k] * gsh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pc_x, pc_z, fsi0_189, fsi0_191, \
                         fsi0_192, fsi1_189, fsi1_191, fsi1_192, \
                         gsh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pa_x[k] * fsi0_189[k]
                   - f_10 * pc_x[k] * fsi1_189[k];

        t_190[k] = f_3 * pc_z[k] * gsh_141[k];

        t_191[k] = pa_x[k] * fsi0_191[k]
                   - f_10 * pc_x[k] * fsi1_191[k];

        t_192[k] = pa_x[k] * fsi0_192[k]
                   - f_10 * pc_x[k] * fsi1_192[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pa_x, pc_x, pc_y, fsi0_193, fsi0_195, fsh_83, \
                         fsi1_193, fsi1_195, gsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pa_x[k] * fsi0_193[k]
                   - f_10 * pc_x[k] * fsi1_193[k];

        t_194[k] = f_13 * fsh_83[k]
                   + f_3 * pc_y[k] * gsh_146[k];

        t_195[k] = pa_x[k] * fsi0_195[k]
                   - f_10 * pc_x[k] * fsi1_195[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_z, pc_y, pc_z, fsi0_84, fsi0_87, \
                         fsh_63, fsh_84, fsi1_84, fsi1_87, gsh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_z[k] * fsi0_84[k]
                   - f_10 * pc_z[k] * fsi1_84[k];

        t_197[k] = f_12 * fsh_84[k]
                   + f_3 * pc_y[k] * gsh_147[k];

        t_198[k] = f_11 * fsh_63[k]
                   + f_3 * pc_z[k] * gsh_147[k];

        t_199[k] = pa_z[k] * fsi0_87[k]
                   - f_10 * pc_z[k] * fsi1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pa_x, pa_z, pc_x, pc_y, pc_z, fsi0_90, fsi0_201, \
                         fsh_86, fsh_152, fsi1_90, fsi1_201, gsh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * fsh_86[k]
                   + f_3 * pc_y[k] * gsh_149[k];

        t_201[k] = pa_x[k] * fsi0_201[k]
                   + f_0 * fsh_152[k]
                   - f_10 * pc_x[k] * fsi1_201[k];

        t_202[k] = pa_z[k] * fsi0_90[k]
                   - f_10 * pc_z[k] * fsi1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pa_x, pc_x, pc_y, pc_z, fsi0_205, fsh_66, \
                         fsh_89, fsh_156, fsi1_205, gsh_150, gsh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * fsh_66[k]
                   + f_3 * pc_z[k] * gsh_150[k];

        t_204[k] = f_12 * fsh_89[k]
                   + f_3 * pc_y[k] * gsh_152[k];

        t_205[k] = pa_x[k] * fsi0_205[k]
                   + f_13 * fsh_156[k]
                   - f_10 * pc_x[k] * fsi1_205[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pa_x, pa_z, pc_x, pc_z, fsi0_94, fsi0_208, \
                         fsh_69, fsh_159, fsi1_94, fsi1_208, gsh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_z[k] * fsi0_94[k]
                   - f_10 * pc_z[k] * fsi1_94[k];

        t_207[k] = f_11 * fsh_69[k]
                   + f_3 * pc_z[k] * gsh_153[k];

        t_208[k] = pa_x[k] * fsi0_208[k]
                   + f_12 * fsh_159[k]
                   - f_10 * pc_x[k] * fsi1_208[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pa_x, pc_x, pc_y, fsi0_210, fsh_93, \
                         fsh_161, fsh_162, fsh_163, fsi1_210, gsh_156, gsh_162, \
                         gsh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_12 * fsh_93[k]
                   + f_3 * pc_y[k] * gsh_156[k];

        t_210[k] = pa_x[k] * fsi0_210[k]
                   + f_12 * fsh_161[k]
                   - f_10 * pc_x[k] * fsi1_210[k];

        t_211[k] = f_11 * fsh_162[k]
                   + f_3 * pc_x[k] * gsh_162[k];

        t_212[k] = f_11 * fsh_163[k]
                   + f_3 * pc_x[k] * gsh_163[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, fsh_164, fsh_165, fsh_166, fsh_167, \
                         gsh_164, gsh_165, gsh_166, gsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_11 * fsh_164[k]
                   + f_3 * pc_x[k] * gsh_164[k];

        t_214[k] = f_11 * fsh_165[k]
                   + f_3 * pc_x[k] * gsh_165[k];

        t_215[k] = f_11 * fsh_166[k]
                   + f_3 * pc_x[k] * gsh_166[k];

        t_216[k] = f_11 * fsh_167[k]
                   + f_3 * pc_x[k] * gsh_167[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_x, pc_x, pc_z, fsi0_217, fsi0_219, \
                         fsi0_220, fsh_78, fsi1_217, fsi1_219, fsi1_220, \
                         gsh_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_x[k] * fsi0_217[k]
                   - f_10 * pc_x[k] * fsi1_217[k];

        t_218[k] = f_11 * fsh_78[k]
                   + f_3 * pc_z[k] * gsh_162[k];

        t_219[k] = pa_x[k] * fsi0_219[k]
                   - f_10 * pc_x[k] * fsi1_219[k];

        t_220[k] = pa_x[k] * fsi0_220[k]
                   - f_10 * pc_x[k] * fsi1_220[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_x, pa_y, pc_x, pc_y, fsi0_140, \
                         fsi0_221, fsi0_223, fsh_104, fsi1_140, fsi1_221, fsi1_223, \
                         gsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pa_x[k] * fsi0_221[k]
                   - f_10 * pc_x[k] * fsi1_221[k];

        t_222[k] = f_12 * fsh_104[k]
                   + f_3 * pc_y[k] * gsh_167[k];

        t_223[k] = pa_x[k] * fsi0_223[k]
                   - f_10 * pc_x[k] * fsi1_223[k];

        t_224[k] = pa_y[k] * fsi0_140[k]
                   - f_10 * pc_y[k] * fsi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_x, pc_x, pc_y, pc_z, fsi0_227, fsh_84, \
                         fsh_105, fsh_107, fsh_171, fsi1_227, gsh_168, \
                         gsh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * fsh_105[k]
                   + f_3 * pc_y[k] * gsh_168[k];

        t_226[k] = f_12 * fsh_84[k]
                   + f_3 * pc_z[k] * gsh_168[k];

        t_227[k] = pa_x[k] * fsi0_227[k]
                   + f_0 * fsh_171[k]
                   - f_10 * pc_x[k] * fsi1_227[k];

        t_228[k] = f_11 * fsh_107[k]
                   + f_3 * pc_y[k] * gsh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pa_x, pa_y, pc_x, pc_y, pc_z, fsi0_145, \
                         fsi0_230, fsh_87, fsh_174, fsi1_145, fsi1_230, \
                         gsh_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * fsi0_145[k]
                   - f_10 * pc_y[k] * fsi1_145[k];

        t_230[k] = pa_x[k] * fsi0_230[k]
                   + f_13 * fsh_174[k]
                   - f_10 * pc_x[k] * fsi1_230[k];

        t_231[k] = f_12 * fsh_87[k]
                   + f_3 * pc_z[k] * gsh_171[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pa_x, pa_y, pc_x, pc_y, fsi0_149, fsi0_234, \
                         fsh_110, fsh_178, fsi1_149, fsi1_234, \
                         gsh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_11 * fsh_110[k]
                   + f_3 * pc_y[k] * gsh_173[k];

        t_233[k] = pa_y[k] * fsi0_149[k]
                   - f_10 * pc_y[k] * fsi1_149[k];

        t_234[k] = pa_x[k] * fsi0_234[k]
                   + f_12 * fsh_178[k]
                   - f_10 * pc_x[k] * fsi1_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pa_x, pc_x, pc_y, pc_z, fsi0_236, fsh_90, \
                         fsh_114, fsh_180, fsi1_236, gsh_174, gsh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_12 * fsh_90[k]
                   + f_3 * pc_z[k] * gsh_174[k];

        t_236[k] = pa_x[k] * fsi0_236[k]
                   + f_12 * fsh_180[k]
                   - f_10 * pc_x[k] * fsi1_236[k];

        t_237[k] = f_11 * fsh_114[k]
                   + f_3 * pc_y[k] * gsh_177[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_y, pc_x, pc_y, fsi0_154, fsh_183, \
                         fsh_184, fsh_185, fsi1_154, gsh_183, gsh_184, \
                         gsh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_y[k] * fsi0_154[k]
                   - f_10 * pc_y[k] * fsi1_154[k];

        t_239[k] = f_11 * fsh_183[k]
                   + f_3 * pc_x[k] * gsh_183[k];

        t_240[k] = f_11 * fsh_184[k]
                   + f_3 * pc_x[k] * gsh_184[k];

        t_241[k] = f_11 * fsh_185[k]
                   + f_3 * pc_x[k] * gsh_185[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_x, pc_x, fsi0_245, fsh_186, fsh_187, \
                         fsh_188, fsi1_245, gsh_186, gsh_187, gsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_11 * fsh_186[k]
                   + f_3 * pc_x[k] * gsh_186[k];

        t_243[k] = f_11 * fsh_187[k]
                   + f_3 * pc_x[k] * gsh_187[k];

        t_244[k] = f_11 * fsh_188[k]
                   + f_3 * pc_x[k] * gsh_188[k];

        t_245[k] = pa_x[k] * fsi0_245[k]
                   - f_10 * pc_x[k] * fsi1_245[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pa_x, pc_x, pc_z, fsi0_247, fsi0_248, \
                         fsi0_249, fsh_99, fsi1_247, fsi1_248, fsi1_249, \
                         gsh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * fsh_99[k]
                   + f_3 * pc_z[k] * gsh_183[k];

        t_247[k] = pa_x[k] * fsi0_247[k]
                   - f_10 * pc_x[k] * fsi1_247[k];

        t_248[k] = pa_x[k] * fsi0_248[k]
                   - f_10 * pc_x[k] * fsi1_248[k];

        t_249[k] = pa_x[k] * fsi0_249[k]
                   - f_10 * pc_x[k] * fsi1_249[k];
    }
}

static auto
compute_prim_gsi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t gsg0, const size_t gsg1,
                                                          const size_t gsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);
    const auto f_16 = 3.0 / q;

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
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsi0_168 = buffer.data(fsi0 + 168);
    const auto *fsi0_169 = buffer.data(fsi0 + 169);
    const auto *fsi0_171 = buffer.data(fsi0 + 171);
    const auto *fsi0_174 = buffer.data(fsi0 + 174);
    const auto *fsi0_178 = buffer.data(fsi0 + 178);
    const auto *fsi0_189 = buffer.data(fsi0 + 189);
    const auto *fsi0_191 = buffer.data(fsi0 + 191);
    const auto *fsi0_192 = buffer.data(fsi0 + 192);
    const auto *fsi0_193 = buffer.data(fsi0 + 193);
    const auto *fsi0_251 = buffer.data(fsi0 + 251);
    const auto *fsi0_252 = buffer.data(fsi0 + 252);
    const auto *fsi0_254 = buffer.data(fsi0 + 254);
    const auto *fsi0_257 = buffer.data(fsi0 + 257);
    const auto *fsi0_261 = buffer.data(fsi0 + 261);
    const auto *fsi0_266 = buffer.data(fsi0 + 266);
    const auto *fsi0_273 = buffer.data(fsi0 + 273);
    const auto *fsi0_274 = buffer.data(fsi0 + 274);
    const auto *fsi0_275 = buffer.data(fsi0 + 275);
    const auto *fsi0_276 = buffer.data(fsi0 + 276);
    const auto *fsi0_277 = buffer.data(fsi0 + 277);
    const auto *fsi0_279 = buffer.data(fsi0 + 279);

    const auto *fsh_105 = buffer.data(fsh + 105);
    const auto *fsh_125 = buffer.data(fsh + 125);
    const auto *fsh_141 = buffer.data(fsh + 141);
    const auto *fsh_142 = buffer.data(fsh + 142);
    const auto *fsh_143 = buffer.data(fsh + 143);
    const auto *fsh_144 = buffer.data(fsh + 144);
    const auto *fsh_146 = buffer.data(fsh + 146);
    const auto *fsh_162 = buffer.data(fsh + 162);
    const auto *fsh_167 = buffer.data(fsh + 167);
    const auto *fsh_183 = buffer.data(fsh + 183);
    const auto *fsh_185 = buffer.data(fsh + 185);
    const auto *fsh_186 = buffer.data(fsh + 186);
    const auto *fsh_187 = buffer.data(fsh + 187);
    const auto *fsh_188 = buffer.data(fsh + 188);
    const auto *fsh_189 = buffer.data(fsh + 189);
    const auto *fsh_194 = buffer.data(fsh + 194);
    const auto *fsh_198 = buffer.data(fsh + 198);
    const auto *fsh_203 = buffer.data(fsh + 203);
    const auto *fsh_204 = buffer.data(fsh + 204);
    const auto *fsh_205 = buffer.data(fsh + 205);
    const auto *fsh_206 = buffer.data(fsh + 206);
    const auto *fsh_207 = buffer.data(fsh + 207);
    const auto *fsh_209 = buffer.data(fsh + 209);

    const auto *fsi1_168 = buffer.data(fsi1 + 168);
    const auto *fsi1_169 = buffer.data(fsi1 + 169);
    const auto *fsi1_171 = buffer.data(fsi1 + 171);
    const auto *fsi1_174 = buffer.data(fsi1 + 174);
    const auto *fsi1_178 = buffer.data(fsi1 + 178);
    const auto *fsi1_189 = buffer.data(fsi1 + 189);
    const auto *fsi1_191 = buffer.data(fsi1 + 191);
    const auto *fsi1_192 = buffer.data(fsi1 + 192);
    const auto *fsi1_193 = buffer.data(fsi1 + 193);
    const auto *fsi1_251 = buffer.data(fsi1 + 251);
    const auto *fsi1_252 = buffer.data(fsi1 + 252);
    const auto *fsi1_254 = buffer.data(fsi1 + 254);
    const auto *fsi1_257 = buffer.data(fsi1 + 257);
    const auto *fsi1_261 = buffer.data(fsi1 + 261);
    const auto *fsi1_266 = buffer.data(fsi1 + 266);
    const auto *fsi1_273 = buffer.data(fsi1 + 273);
    const auto *fsi1_274 = buffer.data(fsi1 + 274);
    const auto *fsi1_275 = buffer.data(fsi1 + 275);
    const auto *fsi1_276 = buffer.data(fsi1 + 276);
    const auto *fsi1_277 = buffer.data(fsi1 + 277);
    const auto *fsi1_279 = buffer.data(fsi1 + 279);

    const auto *gsg0_135 = buffer.data(gsg0 + 135);
    const auto *gsg0_136 = buffer.data(gsg0 + 136);
    const auto *gsg0_137 = buffer.data(gsg0 + 137);
    const auto *gsg0_138 = buffer.data(gsg0 + 138);
    const auto *gsg0_139 = buffer.data(gsg0 + 139);
    const auto *gsg0_140 = buffer.data(gsg0 + 140);
    const auto *gsg0_150 = buffer.data(gsg0 + 150);
    const auto *gsg0_151 = buffer.data(gsg0 + 151);
    const auto *gsg0_153 = buffer.data(gsg0 + 153);
    const auto *gsg0_155 = buffer.data(gsg0 + 155);
    const auto *gsg0_156 = buffer.data(gsg0 + 156);
    const auto *gsg0_158 = buffer.data(gsg0 + 158);
    const auto *gsg0_159 = buffer.data(gsg0 + 159);
    const auto *gsg0_160 = buffer.data(gsg0 + 160);
    const auto *gsg0_161 = buffer.data(gsg0 + 161);
    const auto *gsg0_162 = buffer.data(gsg0 + 162);
    const auto *gsg0_163 = buffer.data(gsg0 + 163);
    const auto *gsg0_164 = buffer.data(gsg0 + 164);
    const auto *gsg0_167 = buffer.data(gsg0 + 167);
    const auto *gsg0_169 = buffer.data(gsg0 + 169);
    const auto *gsg0_170 = buffer.data(gsg0 + 170);
    const auto *gsg0_172 = buffer.data(gsg0 + 172);
    const auto *gsg0_173 = buffer.data(gsg0 + 173);
    const auto *gsg0_174 = buffer.data(gsg0 + 174);
    const auto *gsg0_176 = buffer.data(gsg0 + 176);
    const auto *gsg0_177 = buffer.data(gsg0 + 177);
    const auto *gsg0_178 = buffer.data(gsg0 + 178);
    const auto *gsg0_179 = buffer.data(gsg0 + 179);
    const auto *gsg0_180 = buffer.data(gsg0 + 180);
    const auto *gsg0_181 = buffer.data(gsg0 + 181);
    const auto *gsg0_182 = buffer.data(gsg0 + 182);
    const auto *gsg0_183 = buffer.data(gsg0 + 183);
    const auto *gsg0_184 = buffer.data(gsg0 + 184);
    const auto *gsg0_185 = buffer.data(gsg0 + 185);
    const auto *gsg0_186 = buffer.data(gsg0 + 186);
    const auto *gsg0_187 = buffer.data(gsg0 + 187);
    const auto *gsg0_188 = buffer.data(gsg0 + 188);
    const auto *gsg0_189 = buffer.data(gsg0 + 189);
    const auto *gsg0_190 = buffer.data(gsg0 + 190);
    const auto *gsg0_191 = buffer.data(gsg0 + 191);
    const auto *gsg0_192 = buffer.data(gsg0 + 192);
    const auto *gsg0_193 = buffer.data(gsg0 + 193);
    const auto *gsg0_194 = buffer.data(gsg0 + 194);
    const auto *gsg0_196 = buffer.data(gsg0 + 196);
    const auto *gsg0_198 = buffer.data(gsg0 + 198);
    const auto *gsg0_199 = buffer.data(gsg0 + 199);
    const auto *gsg0_201 = buffer.data(gsg0 + 201);

    const auto *gsg1_135 = buffer.data(gsg1 + 135);
    const auto *gsg1_136 = buffer.data(gsg1 + 136);
    const auto *gsg1_137 = buffer.data(gsg1 + 137);
    const auto *gsg1_138 = buffer.data(gsg1 + 138);
    const auto *gsg1_139 = buffer.data(gsg1 + 139);
    const auto *gsg1_140 = buffer.data(gsg1 + 140);
    const auto *gsg1_150 = buffer.data(gsg1 + 150);
    const auto *gsg1_151 = buffer.data(gsg1 + 151);
    const auto *gsg1_153 = buffer.data(gsg1 + 153);
    const auto *gsg1_155 = buffer.data(gsg1 + 155);
    const auto *gsg1_156 = buffer.data(gsg1 + 156);
    const auto *gsg1_158 = buffer.data(gsg1 + 158);
    const auto *gsg1_159 = buffer.data(gsg1 + 159);
    const auto *gsg1_160 = buffer.data(gsg1 + 160);
    const auto *gsg1_161 = buffer.data(gsg1 + 161);
    const auto *gsg1_162 = buffer.data(gsg1 + 162);
    const auto *gsg1_163 = buffer.data(gsg1 + 163);
    const auto *gsg1_164 = buffer.data(gsg1 + 164);
    const auto *gsg1_167 = buffer.data(gsg1 + 167);
    const auto *gsg1_169 = buffer.data(gsg1 + 169);
    const auto *gsg1_170 = buffer.data(gsg1 + 170);
    const auto *gsg1_172 = buffer.data(gsg1 + 172);
    const auto *gsg1_173 = buffer.data(gsg1 + 173);
    const auto *gsg1_174 = buffer.data(gsg1 + 174);
    const auto *gsg1_176 = buffer.data(gsg1 + 176);
    const auto *gsg1_177 = buffer.data(gsg1 + 177);
    const auto *gsg1_178 = buffer.data(gsg1 + 178);
    const auto *gsg1_179 = buffer.data(gsg1 + 179);
    const auto *gsg1_180 = buffer.data(gsg1 + 180);
    const auto *gsg1_181 = buffer.data(gsg1 + 181);
    const auto *gsg1_182 = buffer.data(gsg1 + 182);
    const auto *gsg1_183 = buffer.data(gsg1 + 183);
    const auto *gsg1_184 = buffer.data(gsg1 + 184);
    const auto *gsg1_185 = buffer.data(gsg1 + 185);
    const auto *gsg1_186 = buffer.data(gsg1 + 186);
    const auto *gsg1_187 = buffer.data(gsg1 + 187);
    const auto *gsg1_188 = buffer.data(gsg1 + 188);
    const auto *gsg1_189 = buffer.data(gsg1 + 189);
    const auto *gsg1_190 = buffer.data(gsg1 + 190);
    const auto *gsg1_191 = buffer.data(gsg1 + 191);
    const auto *gsg1_192 = buffer.data(gsg1 + 192);
    const auto *gsg1_193 = buffer.data(gsg1 + 193);
    const auto *gsg1_194 = buffer.data(gsg1 + 194);
    const auto *gsg1_196 = buffer.data(gsg1 + 196);
    const auto *gsg1_198 = buffer.data(gsg1 + 198);
    const auto *gsg1_199 = buffer.data(gsg1 + 199);
    const auto *gsg1_201 = buffer.data(gsg1 + 201);

    const auto *gsh_188 = buffer.data(gsh + 188);
    const auto *gsh_189 = buffer.data(gsh + 189);
    const auto *gsh_190 = buffer.data(gsh + 190);
    const auto *gsh_191 = buffer.data(gsh + 191);
    const auto *gsh_192 = buffer.data(gsh + 192);
    const auto *gsh_193 = buffer.data(gsh + 193);
    const auto *gsh_194 = buffer.data(gsh + 194);
    const auto *gsh_195 = buffer.data(gsh + 195);
    const auto *gsh_196 = buffer.data(gsh + 196);
    const auto *gsh_197 = buffer.data(gsh + 197);
    const auto *gsh_198 = buffer.data(gsh + 198);
    const auto *gsh_203 = buffer.data(gsh + 203);
    const auto *gsh_204 = buffer.data(gsh + 204);
    const auto *gsh_205 = buffer.data(gsh + 205);
    const auto *gsh_206 = buffer.data(gsh + 206);
    const auto *gsh_207 = buffer.data(gsh + 207);
    const auto *gsh_209 = buffer.data(gsh + 209);
    const auto *gsh_210 = buffer.data(gsh + 210);
    const auto *gsh_211 = buffer.data(gsh + 211);
    const auto *gsh_213 = buffer.data(gsh + 213);
    const auto *gsh_215 = buffer.data(gsh + 215);
    const auto *gsh_216 = buffer.data(gsh + 216);
    const auto *gsh_218 = buffer.data(gsh + 218);
    const auto *gsh_219 = buffer.data(gsh + 219);
    const auto *gsh_220 = buffer.data(gsh + 220);
    const auto *gsh_222 = buffer.data(gsh + 222);
    const auto *gsh_223 = buffer.data(gsh + 223);
    const auto *gsh_224 = buffer.data(gsh + 224);
    const auto *gsh_225 = buffer.data(gsh + 225);
    const auto *gsh_226 = buffer.data(gsh + 226);
    const auto *gsh_227 = buffer.data(gsh + 227);
    const auto *gsh_228 = buffer.data(gsh + 228);
    const auto *gsh_229 = buffer.data(gsh + 229);
    const auto *gsh_230 = buffer.data(gsh + 230);
    const auto *gsh_233 = buffer.data(gsh + 233);
    const auto *gsh_235 = buffer.data(gsh + 235);
    const auto *gsh_236 = buffer.data(gsh + 236);
    const auto *gsh_238 = buffer.data(gsh + 238);
    const auto *gsh_239 = buffer.data(gsh + 239);
    const auto *gsh_240 = buffer.data(gsh + 240);
    const auto *gsh_242 = buffer.data(gsh + 242);
    const auto *gsh_243 = buffer.data(gsh + 243);
    const auto *gsh_244 = buffer.data(gsh + 244);
    const auto *gsh_245 = buffer.data(gsh + 245);
    const auto *gsh_246 = buffer.data(gsh + 246);
    const auto *gsh_247 = buffer.data(gsh + 247);
    const auto *gsh_248 = buffer.data(gsh + 248);
    const auto *gsh_249 = buffer.data(gsh + 249);
    const auto *gsh_250 = buffer.data(gsh + 250);
    const auto *gsh_251 = buffer.data(gsh + 251);
    const auto *gsh_252 = buffer.data(gsh + 252);
    const auto *gsh_253 = buffer.data(gsh + 253);
    const auto *gsh_254 = buffer.data(gsh + 254);
    const auto *gsh_255 = buffer.data(gsh + 255);
    const auto *gsh_256 = buffer.data(gsh + 256);
    const auto *gsh_257 = buffer.data(gsh + 257);
    const auto *gsh_258 = buffer.data(gsh + 258);
    const auto *gsh_259 = buffer.data(gsh + 259);
    const auto *gsh_260 = buffer.data(gsh + 260);
    const auto *gsh_261 = buffer.data(gsh + 261);
    const auto *gsh_262 = buffer.data(gsh + 262);
    const auto *gsh_263 = buffer.data(gsh + 263);
    const auto *gsh_264 = buffer.data(gsh + 264);
    const auto *gsh_265 = buffer.data(gsh + 265);
    const auto *gsh_266 = buffer.data(gsh + 266);
    const auto *gsh_267 = buffer.data(gsh + 267);
    const auto *gsh_268 = buffer.data(gsh + 268);
    const auto *gsh_269 = buffer.data(gsh + 269);
    const auto *gsh_270 = buffer.data(gsh + 270);
    const auto *gsh_271 = buffer.data(gsh + 271);
    const auto *gsh_272 = buffer.data(gsh + 272);
    const auto *gsh_274 = buffer.data(gsh + 274);
    const auto *gsh_276 = buffer.data(gsh + 276);
    const auto *gsh_277 = buffer.data(gsh + 277);
    const auto *gsh_279 = buffer.data(gsh + 279);

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_x, pc_x, pc_y, fsi0_251, fsi0_252, \
                         fsh_125, fsh_189, fsi1_251, fsi1_252, gsh_188, \
                         gsh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * fsh_125[k]
                   + f_3 * pc_y[k] * gsh_188[k];

        t_251[k] = pa_x[k] * fsi0_251[k]
                   - f_10 * pc_x[k] * fsi1_251[k];

        t_252[k] = pa_x[k] * fsi0_252[k]
                   + f_16 * fsh_189[k]
                   - f_10 * pc_x[k] * fsi1_252[k];

        t_253[k] = f_3 * pc_y[k] * gsh_189[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_y, pc_z, fsh_105, gsg0_135, gsg1_135, \
                         gsh_189, gsh_190, gsh_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * fsh_105[k]
                   + f_3 * pc_z[k] * gsh_189[k];

        t_255[k] = f_4 * gsg0_135[k]
                   - f_5 * gsg1_135[k]
                   + f_3 * pc_y[k] * gsh_190[k];

        t_256[k] = f_3 * pc_y[k] * gsh_191[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_x, pc_x, pc_y, fsi0_257, fsh_194, fsi1_257, \
                         gsg0_136, gsg0_137, gsg1_136, gsg1_137, gsh_192, \
                         gsh_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = pa_x[k] * fsi0_257[k]
                   + f_0 * fsh_194[k]
                   - f_10 * pc_x[k] * fsi1_257[k];

        t_258[k] = f_6 * gsg0_136[k]
                   - f_7 * gsg1_136[k]
                   + f_3 * pc_y[k] * gsh_192[k];

        t_259[k] = f_4 * gsg0_137[k]
                   - f_5 * gsg1_137[k]
                   + f_3 * pc_y[k] * gsh_193[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_x, pc_x, pc_y, fsi0_261, fsh_198, fsi1_261, \
                         gsg0_138, gsg1_138, gsh_194, gsh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_3 * pc_y[k] * gsh_194[k];

        t_261[k] = pa_x[k] * fsi0_261[k]
                   + f_13 * fsh_198[k]
                   - f_10 * pc_x[k] * fsi1_261[k];

        t_262[k] = f_8 * gsg0_138[k]
                   - f_9 * gsg1_138[k]
                   + f_3 * pc_y[k] * gsh_195[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pc_y, gsg0_139, gsg0_140, gsg1_139, gsg1_140, \
                         gsh_196, gsh_197, gsh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_6 * gsg0_139[k]
                   - f_7 * gsg1_139[k]
                   + f_3 * pc_y[k] * gsh_196[k];

        t_264[k] = f_4 * gsg0_140[k]
                   - f_5 * gsg1_140[k]
                   + f_3 * pc_y[k] * gsh_197[k];

        t_265[k] = f_3 * pc_y[k] * gsh_198[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pa_x, pc_x, fsi0_266, fsh_203, fsh_204, \
                         fsh_205, fsh_206, fsi1_266, gsh_204, gsh_205, \
                         gsh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pa_x[k] * fsi0_266[k]
                   + f_12 * fsh_203[k]
                   - f_10 * pc_x[k] * fsi1_266[k];

        t_267[k] = f_11 * fsh_204[k]
                   + f_3 * pc_x[k] * gsh_204[k];

        t_268[k] = f_11 * fsh_205[k]
                   + f_3 * pc_x[k] * gsh_205[k];

        t_269[k] = f_11 * fsh_206[k]
                   + f_3 * pc_x[k] * gsh_206[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pa_x, pc_x, pc_y, fsi0_273, fsh_207, \
                         fsh_209, fsi1_273, gsh_203, gsh_207, gsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_11 * fsh_207[k]
                   + f_3 * pc_x[k] * gsh_207[k];

        t_271[k] = f_3 * pc_y[k] * gsh_203[k];

        t_272[k] = f_11 * fsh_209[k]
                   + f_3 * pc_x[k] * gsh_209[k];

        t_273[k] = pa_x[k] * fsi0_273[k]
                   - f_10 * pc_x[k] * fsi1_273[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pa_x, pc_x, fsi0_274, fsi0_275, fsi0_276, \
                         fsi0_277, fsi1_274, fsi1_275, fsi1_276, \
                         fsi1_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = pa_x[k] * fsi0_274[k]
                   - f_10 * pc_x[k] * fsi1_274[k];

        t_275[k] = pa_x[k] * fsi0_275[k]
                   - f_10 * pc_x[k] * fsi1_275[k];

        t_276[k] = pa_x[k] * fsi0_276[k]
                   - f_10 * pc_x[k] * fsi1_276[k];

        t_277[k] = pa_x[k] * fsi0_277[k]
                   - f_10 * pc_x[k] * fsi1_277[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pa_x, pc_x, pc_y, fsi0_279, fsi1_279, \
                         gsg0_150, gsg0_151, gsg1_150, gsg1_151, gsh_209, gsh_210, \
                         gsh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_3 * pc_y[k] * gsh_209[k];

        t_279[k] = pa_x[k] * fsi0_279[k]
                   - f_10 * pc_x[k] * fsi1_279[k];

        t_280[k] = f_1 * gsg0_150[k]
                   - f_2 * gsg1_150[k]
                   + f_3 * pc_x[k] * gsh_210[k];

        t_281[k] = f_14 * gsg0_151[k]
                   - f_15 * gsg1_151[k]
                   + f_3 * pc_x[k] * gsh_211[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pc_x, pc_z, gsg0_153, gsg0_155, gsg1_153, \
                         gsg1_155, gsh_210, gsh_211, gsh_213, gsh_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_3 * pc_z[k] * gsh_210[k];

        t_283[k] = f_8 * gsg0_153[k]
                   - f_9 * gsg1_153[k]
                   + f_3 * pc_x[k] * gsh_213[k];

        t_284[k] = f_3 * pc_z[k] * gsh_211[k];

        t_285[k] = f_8 * gsg0_155[k]
                   - f_9 * gsg1_155[k]
                   + f_3 * pc_x[k] * gsh_215[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pc_x, pc_z, gsg0_156, gsg0_158, gsg0_159, \
                         gsg1_156, gsg1_158, gsg1_159, gsh_213, gsh_216, gsh_218, \
                         gsh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_6 * gsg0_156[k]
                   - f_7 * gsg1_156[k]
                   + f_3 * pc_x[k] * gsh_216[k];

        t_287[k] = f_3 * pc_z[k] * gsh_213[k];

        t_288[k] = f_6 * gsg0_158[k]
                   - f_7 * gsg1_158[k]
                   + f_3 * pc_x[k] * gsh_218[k];

        t_289[k] = f_6 * gsg0_159[k]
                   - f_7 * gsg1_159[k]
                   + f_3 * pc_x[k] * gsh_219[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pc_x, pc_z, gsg0_160, gsg0_162, gsg0_163, \
                         gsg1_160, gsg1_162, gsg1_163, gsh_216, gsh_220, gsh_222, \
                         gsh_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_4 * gsg0_160[k]
                   - f_5 * gsg1_160[k]
                   + f_3 * pc_x[k] * gsh_220[k];

        t_291[k] = f_3 * pc_z[k] * gsh_216[k];

        t_292[k] = f_4 * gsg0_162[k]
                   - f_5 * gsg1_162[k]
                   + f_3 * pc_x[k] * gsh_222[k];

        t_293[k] = f_4 * gsg0_163[k]
                   - f_5 * gsg1_163[k]
                   + f_3 * pc_x[k] * gsh_223[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, pc_x, gsg0_164, gsg1_164, \
                         gsh_224, gsh_225, gsh_226, gsh_227, gsh_228, \
                         gsh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_4 * gsg0_164[k]
                   - f_5 * gsg1_164[k]
                   + f_3 * pc_x[k] * gsh_224[k];

        t_295[k] = f_3 * pc_x[k] * gsh_225[k];

        t_296[k] = f_3 * pc_x[k] * gsh_226[k];

        t_297[k] = f_3 * pc_x[k] * gsh_227[k];

        t_298[k] = f_3 * pc_x[k] * gsh_228[k];

        t_299[k] = f_3 * pc_x[k] * gsh_229[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, fsh_141, gsg0_160, \
                         gsg1_160, gsh_225, gsh_226, gsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_3 * pc_x[k] * gsh_230[k];

        t_301[k] = f_0 * fsh_141[k]
                   + f_1 * gsg0_160[k]
                   - f_2 * gsg1_160[k]
                   + f_3 * pc_y[k] * gsh_225[k];

        t_302[k] = f_3 * pc_z[k] * gsh_225[k];

        t_303[k] = f_4 * gsg0_160[k]
                   - f_5 * gsg1_160[k]
                   + f_3 * pc_z[k] * gsh_226[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_y, pc_z, fsh_146, gsg0_161, gsg0_162, \
                         gsg0_164, gsg1_161, gsg1_162, gsg1_164, gsh_227, gsh_228, \
                         gsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_6 * gsg0_161[k]
                   - f_7 * gsg1_161[k]
                   + f_3 * pc_z[k] * gsh_227[k];

        t_305[k] = f_8 * gsg0_162[k]
                   - f_9 * gsg1_162[k]
                   + f_3 * pc_z[k] * gsh_228[k];

        t_306[k] = f_0 * fsh_146[k]
                   + f_3 * pc_y[k] * gsh_230[k];

        t_307[k] = f_1 * gsg0_164[k]
                   - f_2 * gsg1_164[k]
                   + f_3 * pc_z[k] * gsh_230[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_z, pc_x, pc_z, fsi0_168, fsi0_169, \
                         fsi0_171, fsi1_168, fsi1_169, fsi1_171, gsg0_167, gsg1_167, \
                         gsh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pa_z[k] * fsi0_168[k]
                   - f_10 * pc_z[k] * fsi1_168[k];

        t_309[k] = pa_z[k] * fsi0_169[k]
                   - f_10 * pc_z[k] * fsi1_169[k];

        t_310[k] = f_14 * gsg0_167[k]
                   - f_15 * gsg1_167[k]
                   + f_3 * pc_x[k] * gsh_233[k];

        t_311[k] = pa_z[k] * fsi0_171[k]
                   - f_10 * pc_z[k] * fsi1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_z, pc_x, pc_z, fsi0_174, fsi1_174, gsg0_169, \
                         gsg0_170, gsg1_169, gsg1_170, gsh_235, \
                         gsh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_8 * gsg0_169[k]
                   - f_9 * gsg1_169[k]
                   + f_3 * pc_x[k] * gsh_235[k];

        t_313[k] = f_8 * gsg0_170[k]
                   - f_9 * gsg1_170[k]
                   + f_3 * pc_x[k] * gsh_236[k];

        t_314[k] = pa_z[k] * fsi0_174[k]
                   - f_10 * pc_z[k] * fsi1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pc_x, gsg0_172, gsg0_173, gsg0_174, gsg1_172, \
                         gsg1_173, gsg1_174, gsh_238, gsh_239, \
                         gsh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_6 * gsg0_172[k]
                   - f_7 * gsg1_172[k]
                   + f_3 * pc_x[k] * gsh_238[k];

        t_316[k] = f_6 * gsg0_173[k]
                   - f_7 * gsg1_173[k]
                   + f_3 * pc_x[k] * gsh_239[k];

        t_317[k] = f_6 * gsg0_174[k]
                   - f_7 * gsg1_174[k]
                   + f_3 * pc_x[k] * gsh_240[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pa_z, pc_x, pc_z, fsi0_178, fsi1_178, gsg0_176, \
                         gsg0_177, gsg1_176, gsg1_177, gsh_242, \
                         gsh_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pa_z[k] * fsi0_178[k]
                   - f_10 * pc_z[k] * fsi1_178[k];

        t_319[k] = f_4 * gsg0_176[k]
                   - f_5 * gsg1_176[k]
                   + f_3 * pc_x[k] * gsh_242[k];

        t_320[k] = f_4 * gsg0_177[k]
                   - f_5 * gsg1_177[k]
                   + f_3 * pc_x[k] * gsh_243[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, t_325, pc_x, gsg0_178, gsg0_179, \
                         gsg1_178, gsg1_179, gsh_244, gsh_245, gsh_246, gsh_247, \
                         gsh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_4 * gsg0_178[k]
                   - f_5 * gsg1_178[k]
                   + f_3 * pc_x[k] * gsh_244[k];

        t_322[k] = f_4 * gsg0_179[k]
                   - f_5 * gsg1_179[k]
                   + f_3 * pc_x[k] * gsh_245[k];

        t_323[k] = f_3 * pc_x[k] * gsh_246[k];

        t_324[k] = f_3 * pc_x[k] * gsh_247[k];

        t_325[k] = f_3 * pc_x[k] * gsh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, pa_z, pc_x, pc_z, fsi0_189, \
                         fsh_141, fsi1_189, gsh_246, gsh_249, gsh_250, \
                         gsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_3 * pc_x[k] * gsh_249[k];

        t_327[k] = f_3 * pc_x[k] * gsh_250[k];

        t_328[k] = f_3 * pc_x[k] * gsh_251[k];

        t_329[k] = pa_z[k] * fsi0_189[k]
                   - f_10 * pc_z[k] * fsi1_189[k];

        t_330[k] = f_11 * fsh_141[k]
                   + f_3 * pc_z[k] * gsh_246[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pa_z, pc_z, fsi0_191, fsi0_192, fsi0_193, \
                         fsh_142, fsh_143, fsh_144, fsi1_191, fsi1_192, \
                         fsi1_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_z[k] * fsi0_191[k]
                   + f_12 * fsh_142[k]
                   - f_10 * pc_z[k] * fsi1_191[k];

        t_332[k] = pa_z[k] * fsi0_192[k]
                   + f_13 * fsh_143[k]
                   - f_10 * pc_z[k] * fsi1_192[k];

        t_333[k] = pa_z[k] * fsi0_193[k]
                   + f_0 * fsh_144[k]
                   - f_10 * pc_z[k] * fsi1_193[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pc_x, pc_y, pc_z, fsh_146, fsh_167, gsg0_179, \
                         gsg0_180, gsg1_179, gsg1_180, gsh_251, \
                         gsh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_13 * fsh_167[k]
                   + f_3 * pc_y[k] * gsh_251[k];

        t_335[k] = f_11 * fsh_146[k]
                   + f_1 * gsg0_179[k]
                   - f_2 * gsg1_179[k]
                   + f_3 * pc_z[k] * gsh_251[k];

        t_336[k] = f_1 * gsg0_180[k]
                   - f_2 * gsg1_180[k]
                   + f_3 * pc_x[k] * gsh_252[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pc_x, gsg0_181, gsg0_182, gsg0_183, gsg1_181, \
                         gsg1_182, gsg1_183, gsh_253, gsh_254, \
                         gsh_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_14 * gsg0_181[k]
                   - f_15 * gsg1_181[k]
                   + f_3 * pc_x[k] * gsh_253[k];

        t_338[k] = f_14 * gsg0_182[k]
                   - f_15 * gsg1_182[k]
                   + f_3 * pc_x[k] * gsh_254[k];

        t_339[k] = f_8 * gsg0_183[k]
                   - f_9 * gsg1_183[k]
                   + f_3 * pc_x[k] * gsh_255[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pc_x, gsg0_184, gsg0_185, gsg0_186, gsg1_184, \
                         gsg1_185, gsg1_186, gsh_256, gsh_257, \
                         gsh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_8 * gsg0_184[k]
                   - f_9 * gsg1_184[k]
                   + f_3 * pc_x[k] * gsh_256[k];

        t_341[k] = f_8 * gsg0_185[k]
                   - f_9 * gsg1_185[k]
                   + f_3 * pc_x[k] * gsh_257[k];

        t_342[k] = f_6 * gsg0_186[k]
                   - f_7 * gsg1_186[k]
                   + f_3 * pc_x[k] * gsh_258[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pc_x, gsg0_187, gsg0_188, gsg0_189, gsg1_187, \
                         gsg1_188, gsg1_189, gsh_259, gsh_260, \
                         gsh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_6 * gsg0_187[k]
                   - f_7 * gsg1_187[k]
                   + f_3 * pc_x[k] * gsh_259[k];

        t_344[k] = f_6 * gsg0_188[k]
                   - f_7 * gsg1_188[k]
                   + f_3 * pc_x[k] * gsh_260[k];

        t_345[k] = f_6 * gsg0_189[k]
                   - f_7 * gsg1_189[k]
                   + f_3 * pc_x[k] * gsh_261[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, pc_x, gsg0_190, gsg0_191, gsg0_192, gsg1_190, \
                         gsg1_191, gsg1_192, gsh_262, gsh_263, \
                         gsh_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_4 * gsg0_190[k]
                   - f_5 * gsg1_190[k]
                   + f_3 * pc_x[k] * gsh_262[k];

        t_347[k] = f_4 * gsg0_191[k]
                   - f_5 * gsg1_191[k]
                   + f_3 * pc_x[k] * gsh_263[k];

        t_348[k] = f_4 * gsg0_192[k]
                   - f_5 * gsg1_192[k]
                   + f_3 * pc_x[k] * gsh_264[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, pc_x, gsg0_193, gsg0_194, \
                         gsg1_193, gsg1_194, gsh_265, gsh_266, gsh_267, gsh_268, \
                         gsh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_4 * gsg0_193[k]
                   - f_5 * gsg1_193[k]
                   + f_3 * pc_x[k] * gsh_265[k];

        t_350[k] = f_4 * gsg0_194[k]
                   - f_5 * gsg1_194[k]
                   + f_3 * pc_x[k] * gsh_266[k];

        t_351[k] = f_3 * pc_x[k] * gsh_267[k];

        t_352[k] = f_3 * pc_x[k] * gsh_268[k];

        t_353[k] = f_3 * pc_x[k] * gsh_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, t_358, pc_x, pc_y, pc_z, fsh_162, \
                         fsh_183, gsg0_190, gsg1_190, gsh_267, gsh_270, gsh_271, \
                         gsh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_3 * pc_x[k] * gsh_270[k];

        t_355[k] = f_3 * pc_x[k] * gsh_271[k];

        t_356[k] = f_3 * pc_x[k] * gsh_272[k];

        t_357[k] = f_12 * fsh_183[k]
                   + f_1 * gsg0_190[k]
                   - f_2 * gsg1_190[k]
                   + f_3 * pc_y[k] * gsh_267[k];

        t_358[k] = f_12 * fsh_162[k]
                   + f_3 * pc_z[k] * gsh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, fsh_185, fsh_186, fsh_187, gsg0_192, \
                         gsg0_193, gsg0_194, gsg1_192, gsg1_193, gsg1_194, gsh_269, gsh_270, \
                         gsh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * fsh_185[k]
                   + f_8 * gsg0_192[k]
                   - f_9 * gsg1_192[k]
                   + f_3 * pc_y[k] * gsh_269[k];

        t_360[k] = f_12 * fsh_186[k]
                   + f_6 * gsg0_193[k]
                   - f_7 * gsg1_193[k]
                   + f_3 * pc_y[k] * gsh_270[k];

        t_361[k] = f_12 * fsh_187[k]
                   + f_4 * gsg0_194[k]
                   - f_5 * gsg1_194[k]
                   + f_3 * pc_y[k] * gsh_271[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pa_y, pc_y, pc_z, fsi0_252, fsh_167, fsh_188, \
                         fsi1_252, gsg0_194, gsg1_194, gsh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * fsh_188[k]
                   + f_3 * pc_y[k] * gsh_272[k];

        t_363[k] = f_12 * fsh_167[k]
                   + f_1 * gsg0_194[k]
                   - f_2 * gsg1_194[k]
                   + f_3 * pc_z[k] * gsh_272[k];

        t_364[k] = pa_y[k] * fsi0_252[k]
                   - f_10 * pc_y[k] * fsi1_252[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, pa_y, pc_x, pc_y, fsi0_254, fsi1_254, gsg0_196, \
                         gsg0_198, gsg1_196, gsg1_198, gsh_274, \
                         gsh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_14 * gsg0_196[k]
                   - f_15 * gsg1_196[k]
                   + f_3 * pc_x[k] * gsh_274[k];

        t_366[k] = pa_y[k] * fsi0_254[k]
                   - f_10 * pc_y[k] * fsi1_254[k];

        t_367[k] = f_8 * gsg0_198[k]
                   - f_9 * gsg1_198[k]
                   + f_3 * pc_x[k] * gsh_276[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, pa_y, pc_x, pc_y, fsi0_257, fsi1_257, gsg0_199, \
                         gsg0_201, gsg1_199, gsg1_201, gsh_277, \
                         gsh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_8 * gsg0_199[k]
                   - f_9 * gsg1_199[k]
                   + f_3 * pc_x[k] * gsh_277[k];

        t_369[k] = pa_y[k] * fsi0_257[k]
                   - f_10 * pc_y[k] * fsi1_257[k];

        t_370[k] = f_6 * gsg0_201[k]
                   - f_7 * gsg1_201[k]
                   + f_3 * pc_x[k] * gsh_279[k];
    }
}

static auto
compute_prim_gsi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t gsg0, const size_t gsg1,
                                                          const size_t gsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);
    const auto f_16 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsi0_261 = buffer.data(fsi0 + 261);
    const auto *fsi0_266 = buffer.data(fsi0 + 266);
    const auto *fsi0_273 = buffer.data(fsi0 + 273);
    const auto *fsi0_275 = buffer.data(fsi0 + 275);
    const auto *fsi0_276 = buffer.data(fsi0 + 276);
    const auto *fsi0_277 = buffer.data(fsi0 + 277);
    const auto *fsi0_279 = buffer.data(fsi0 + 279);

    const auto *fsh_183 = buffer.data(fsh + 183);
    const auto *fsh_204 = buffer.data(fsh + 204);
    const auto *fsh_206 = buffer.data(fsh + 206);
    const auto *fsh_207 = buffer.data(fsh + 207);
    const auto *fsh_208 = buffer.data(fsh + 208);
    const auto *fsh_209 = buffer.data(fsh + 209);

    const auto *fsi1_261 = buffer.data(fsi1 + 261);
    const auto *fsi1_266 = buffer.data(fsi1 + 266);
    const auto *fsi1_273 = buffer.data(fsi1 + 273);
    const auto *fsi1_275 = buffer.data(fsi1 + 275);
    const auto *fsi1_276 = buffer.data(fsi1 + 276);
    const auto *fsi1_277 = buffer.data(fsi1 + 277);
    const auto *fsi1_279 = buffer.data(fsi1 + 279);

    const auto *gsg0_202 = buffer.data(gsg0 + 202);
    const auto *gsg0_203 = buffer.data(gsg0 + 203);
    const auto *gsg0_205 = buffer.data(gsg0 + 205);
    const auto *gsg0_206 = buffer.data(gsg0 + 206);
    const auto *gsg0_207 = buffer.data(gsg0 + 207);
    const auto *gsg0_208 = buffer.data(gsg0 + 208);
    const auto *gsg0_210 = buffer.data(gsg0 + 210);
    const auto *gsg0_212 = buffer.data(gsg0 + 212);
    const auto *gsg0_213 = buffer.data(gsg0 + 213);
    const auto *gsg0_215 = buffer.data(gsg0 + 215);
    const auto *gsg0_216 = buffer.data(gsg0 + 216);
    const auto *gsg0_217 = buffer.data(gsg0 + 217);
    const auto *gsg0_219 = buffer.data(gsg0 + 219);
    const auto *gsg0_220 = buffer.data(gsg0 + 220);
    const auto *gsg0_221 = buffer.data(gsg0 + 221);
    const auto *gsg0_222 = buffer.data(gsg0 + 222);
    const auto *gsg0_223 = buffer.data(gsg0 + 223);
    const auto *gsg0_224 = buffer.data(gsg0 + 224);

    const auto *gsg1_202 = buffer.data(gsg1 + 202);
    const auto *gsg1_203 = buffer.data(gsg1 + 203);
    const auto *gsg1_205 = buffer.data(gsg1 + 205);
    const auto *gsg1_206 = buffer.data(gsg1 + 206);
    const auto *gsg1_207 = buffer.data(gsg1 + 207);
    const auto *gsg1_208 = buffer.data(gsg1 + 208);
    const auto *gsg1_210 = buffer.data(gsg1 + 210);
    const auto *gsg1_212 = buffer.data(gsg1 + 212);
    const auto *gsg1_213 = buffer.data(gsg1 + 213);
    const auto *gsg1_215 = buffer.data(gsg1 + 215);
    const auto *gsg1_216 = buffer.data(gsg1 + 216);
    const auto *gsg1_217 = buffer.data(gsg1 + 217);
    const auto *gsg1_219 = buffer.data(gsg1 + 219);
    const auto *gsg1_220 = buffer.data(gsg1 + 220);
    const auto *gsg1_221 = buffer.data(gsg1 + 221);
    const auto *gsg1_222 = buffer.data(gsg1 + 222);
    const auto *gsg1_223 = buffer.data(gsg1 + 223);
    const auto *gsg1_224 = buffer.data(gsg1 + 224);

    const auto *gsh_280 = buffer.data(gsh + 280);
    const auto *gsh_281 = buffer.data(gsh + 281);
    const auto *gsh_283 = buffer.data(gsh + 283);
    const auto *gsh_284 = buffer.data(gsh + 284);
    const auto *gsh_285 = buffer.data(gsh + 285);
    const auto *gsh_286 = buffer.data(gsh + 286);
    const auto *gsh_288 = buffer.data(gsh + 288);
    const auto *gsh_289 = buffer.data(gsh + 289);
    const auto *gsh_290 = buffer.data(gsh + 290);
    const auto *gsh_291 = buffer.data(gsh + 291);
    const auto *gsh_292 = buffer.data(gsh + 292);
    const auto *gsh_293 = buffer.data(gsh + 293);
    const auto *gsh_294 = buffer.data(gsh + 294);
    const auto *gsh_296 = buffer.data(gsh + 296);
    const auto *gsh_297 = buffer.data(gsh + 297);
    const auto *gsh_299 = buffer.data(gsh + 299);
    const auto *gsh_300 = buffer.data(gsh + 300);
    const auto *gsh_301 = buffer.data(gsh + 301);
    const auto *gsh_303 = buffer.data(gsh + 303);
    const auto *gsh_304 = buffer.data(gsh + 304);
    const auto *gsh_305 = buffer.data(gsh + 305);
    const auto *gsh_306 = buffer.data(gsh + 306);
    const auto *gsh_308 = buffer.data(gsh + 308);
    const auto *gsh_309 = buffer.data(gsh + 309);
    const auto *gsh_310 = buffer.data(gsh + 310);
    const auto *gsh_311 = buffer.data(gsh + 311);
    const auto *gsh_312 = buffer.data(gsh + 312);
    const auto *gsh_313 = buffer.data(gsh + 313);
    const auto *gsh_314 = buffer.data(gsh + 314);

#pragma omp simd aligned(t_371, t_372, t_373, pa_y, pc_x, pc_y, fsi0_261, fsi1_261, gsg0_202, \
                         gsg0_203, gsg1_202, gsg1_203, gsh_280, \
                         gsh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_6 * gsg0_202[k]
                   - f_7 * gsg1_202[k]
                   + f_3 * pc_x[k] * gsh_280[k];

        t_372[k] = f_6 * gsg0_203[k]
                   - f_7 * gsg1_203[k]
                   + f_3 * pc_x[k] * gsh_281[k];

        t_373[k] = pa_y[k] * fsi0_261[k]
                   - f_10 * pc_y[k] * fsi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pc_x, gsg0_205, gsg0_206, gsg0_207, gsg1_205, \
                         gsg1_206, gsg1_207, gsh_283, gsh_284, \
                         gsh_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_4 * gsg0_205[k]
                   - f_5 * gsg1_205[k]
                   + f_3 * pc_x[k] * gsh_283[k];

        t_375[k] = f_4 * gsg0_206[k]
                   - f_5 * gsg1_206[k]
                   + f_3 * pc_x[k] * gsh_284[k];

        t_376[k] = f_4 * gsg0_207[k]
                   - f_5 * gsg1_207[k]
                   + f_3 * pc_x[k] * gsh_285[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, pa_y, pc_x, pc_y, fsi0_266, \
                         fsi1_266, gsg0_208, gsg1_208, gsh_286, gsh_288, gsh_289, \
                         gsh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_4 * gsg0_208[k]
                   - f_5 * gsg1_208[k]
                   + f_3 * pc_x[k] * gsh_286[k];

        t_378[k] = pa_y[k] * fsi0_266[k]
                   - f_10 * pc_y[k] * fsi1_266[k];

        t_379[k] = f_3 * pc_x[k] * gsh_288[k];

        t_380[k] = f_3 * pc_x[k] * gsh_289[k];

        t_381[k] = f_3 * pc_x[k] * gsh_290[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_y, pc_x, pc_y, fsi0_273, fsh_204, \
                         fsi1_273, gsh_291, gsh_292, gsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_3 * pc_x[k] * gsh_291[k];

        t_383[k] = f_3 * pc_x[k] * gsh_292[k];

        t_384[k] = f_3 * pc_x[k] * gsh_293[k];

        t_385[k] = pa_y[k] * fsi0_273[k]
                   + f_16 * fsh_204[k]
                   - f_10 * pc_y[k] * fsi1_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, pa_y, pc_y, pc_z, fsi0_275, fsi0_276, fsh_183, \
                         fsh_206, fsh_207, fsi1_275, fsi1_276, \
                         gsh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_13 * fsh_183[k]
                   + f_3 * pc_z[k] * gsh_288[k];

        t_387[k] = pa_y[k] * fsi0_275[k]
                   + f_0 * fsh_206[k]
                   - f_10 * pc_y[k] * fsi1_275[k];

        t_388[k] = pa_y[k] * fsi0_276[k]
                   + f_13 * fsh_207[k]
                   - f_10 * pc_y[k] * fsi1_276[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, pa_y, pc_y, fsi0_277, fsi0_279, fsh_208, \
                         fsh_209, fsi1_277, fsi1_279, gsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = pa_y[k] * fsi0_277[k]
                   + f_12 * fsh_208[k]
                   - f_10 * pc_y[k] * fsi1_277[k];

        t_390[k] = f_11 * fsh_209[k]
                   + f_3 * pc_y[k] * gsh_293[k];

        t_391[k] = pa_y[k] * fsi0_279[k]
                   - f_10 * pc_y[k] * fsi1_279[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, pc_x, pc_y, gsg0_210, gsg0_212, \
                         gsg0_213, gsg1_210, gsg1_212, gsg1_213, gsh_294, gsh_296, \
                         gsh_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_1 * gsg0_210[k]
                   - f_2 * gsg1_210[k]
                   + f_3 * pc_x[k] * gsh_294[k];

        t_393[k] = f_3 * pc_y[k] * gsh_294[k];

        t_394[k] = f_14 * gsg0_212[k]
                   - f_15 * gsg1_212[k]
                   + f_3 * pc_x[k] * gsh_296[k];

        t_395[k] = f_8 * gsg0_213[k]
                   - f_9 * gsg1_213[k]
                   + f_3 * pc_x[k] * gsh_297[k];

        t_396[k] = f_3 * pc_y[k] * gsh_296[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pc_x, pc_y, gsg0_215, gsg0_216, gsg0_217, \
                         gsg1_215, gsg1_216, gsg1_217, gsh_299, gsh_300, \
                         gsh_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_8 * gsg0_215[k]
                   - f_9 * gsg1_215[k]
                   + f_3 * pc_x[k] * gsh_299[k];

        t_398[k] = f_6 * gsg0_216[k]
                   - f_7 * gsg1_216[k]
                   + f_3 * pc_x[k] * gsh_300[k];

        t_399[k] = f_6 * gsg0_217[k]
                   - f_7 * gsg1_217[k]
                   + f_3 * pc_x[k] * gsh_301[k];

        t_400[k] = f_3 * pc_y[k] * gsh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, gsg0_219, gsg0_220, gsg0_221, gsg1_219, \
                         gsg1_220, gsg1_221, gsh_303, gsh_304, \
                         gsh_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_6 * gsg0_219[k]
                   - f_7 * gsg1_219[k]
                   + f_3 * pc_x[k] * gsh_303[k];

        t_402[k] = f_4 * gsg0_220[k]
                   - f_5 * gsg1_220[k]
                   + f_3 * pc_x[k] * gsh_304[k];

        t_403[k] = f_4 * gsg0_221[k]
                   - f_5 * gsg1_221[k]
                   + f_3 * pc_x[k] * gsh_305[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, pc_x, pc_y, gsg0_222, gsg0_224, \
                         gsg1_222, gsg1_224, gsh_303, gsh_306, gsh_308, gsh_309, \
                         gsh_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * gsg0_222[k]
                   - f_5 * gsg1_222[k]
                   + f_3 * pc_x[k] * gsh_306[k];

        t_405[k] = f_3 * pc_y[k] * gsh_303[k];

        t_406[k] = f_4 * gsg0_224[k]
                   - f_5 * gsg1_224[k]
                   + f_3 * pc_x[k] * gsh_308[k];

        t_407[k] = f_3 * pc_x[k] * gsh_309[k];

        t_408[k] = f_3 * pc_x[k] * gsh_310[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, pc_y, gsg0_220, gsg1_220, \
                         gsh_309, gsh_311, gsh_312, gsh_313, gsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_3 * pc_x[k] * gsh_311[k];

        t_410[k] = f_3 * pc_x[k] * gsh_312[k];

        t_411[k] = f_3 * pc_x[k] * gsh_313[k];

        t_412[k] = f_3 * pc_x[k] * gsh_314[k];

        t_413[k] = f_1 * gsg0_220[k]
                   - f_2 * gsg1_220[k]
                   + f_3 * pc_y[k] * gsh_309[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, gsg0_221, gsg0_222, gsg0_223, gsg1_221, \
                         gsg1_222, gsg1_223, gsh_310, gsh_311, \
                         gsh_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_14 * gsg0_221[k]
                   - f_15 * gsg1_221[k]
                   + f_3 * pc_y[k] * gsh_310[k];

        t_415[k] = f_8 * gsg0_222[k]
                   - f_9 * gsg1_222[k]
                   + f_3 * pc_y[k] * gsh_311[k];

        t_416[k] = f_6 * gsg0_223[k]
                   - f_7 * gsg1_223[k]
                   + f_3 * pc_y[k] * gsh_312[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_y, pc_z, fsh_209, gsg0_224, gsg1_224, \
                         gsh_313, gsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * gsg0_224[k]
                   - f_5 * gsg1_224[k]
                   + f_3 * pc_y[k] * gsh_313[k];

        t_418[k] = f_3 * pc_y[k] * gsh_314[k];

        t_419[k] = f_0 * fsh_209[k]
                   + f_1 * gsg0_224[k]
                   - f_2 * gsg1_224[k]
                   + f_3 * pc_z[k] * gsh_314[k];
    }
}

auto
compute_prim_gsi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fsi0, const size_t fsh,
                                                   const size_t fsi1, const size_t gsg0,
                                                   const size_t gsg1, const size_t gsh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gsi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, fsi0, fsh,
                                                              fsi1, gsg0, gsg1, gsh, ncols,
                                                              gamma, p, q);

    compute_prim_gsi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, fsi0, fsh,
                                                              fsi1, gsg0, gsg1, gsh, ncols,
                                                              gamma, p, q);

    compute_prim_gsi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, fsi0, fsh,
                                                              fsi1, gsg0, gsg1, gsh, ncols,
                                                              gamma, p, q);

    compute_prim_gsi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, fsi0, fsh,
                                                              fsi1, gsg0, gsg1, gsh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
