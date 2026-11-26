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


#include "SimdThreeCenterElectronRepulsionVrrRecGSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gsk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsk0,
                                                          const size_t fsi, const size_t fsk1,
                                                          const size_t gsh0, const size_t gsh1,
                                                          const size_t gsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsk0_0 = buffer.data(fsk0 + 0);
    const auto *fsk0_3 = buffer.data(fsk0 + 3);
    const auto *fsk0_5 = buffer.data(fsk0 + 5);
    const auto *fsk0_6 = buffer.data(fsk0 + 6);
    const auto *fsk0_9 = buffer.data(fsk0 + 9);
    const auto *fsk0_10 = buffer.data(fsk0 + 10);
    const auto *fsk0_14 = buffer.data(fsk0 + 14);
    const auto *fsk0_15 = buffer.data(fsk0 + 15);
    const auto *fsk0_20 = buffer.data(fsk0 + 20);
    const auto *fsk0_28 = buffer.data(fsk0 + 28);
    const auto *fsk0_35 = buffer.data(fsk0 + 35);

    const auto *fsi_0 = buffer.data(fsi + 0);
    const auto *fsi_1 = buffer.data(fsi + 1);
    const auto *fsi_2 = buffer.data(fsi + 2);
    const auto *fsi_3 = buffer.data(fsi + 3);
    const auto *fsi_5 = buffer.data(fsi + 5);
    const auto *fsi_6 = buffer.data(fsi + 6);
    const auto *fsi_9 = buffer.data(fsi + 9);
    const auto *fsi_10 = buffer.data(fsi + 10);
    const auto *fsi_14 = buffer.data(fsi + 14);
    const auto *fsi_21 = buffer.data(fsi + 21);
    const auto *fsi_23 = buffer.data(fsi + 23);
    const auto *fsi_24 = buffer.data(fsi + 24);
    const auto *fsi_25 = buffer.data(fsi + 25);
    const auto *fsi_27 = buffer.data(fsi + 27);
    const auto *fsi_28 = buffer.data(fsi + 28);
    const auto *fsi_33 = buffer.data(fsi + 33);
    const auto *fsi_37 = buffer.data(fsi + 37);
    const auto *fsi_42 = buffer.data(fsi + 42);
    const auto *fsi_49 = buffer.data(fsi + 49);
    const auto *fsi_51 = buffer.data(fsi + 51);
    const auto *fsi_52 = buffer.data(fsi + 52);
    const auto *fsi_53 = buffer.data(fsi + 53);
    const auto *fsi_54 = buffer.data(fsi + 54);
    const auto *fsi_55 = buffer.data(fsi + 55);
    const auto *fsi_77 = buffer.data(fsi + 77);
    const auto *fsi_78 = buffer.data(fsi + 78);
    const auto *fsi_79 = buffer.data(fsi + 79);
    const auto *fsi_80 = buffer.data(fsi + 80);
    const auto *fsi_81 = buffer.data(fsi + 81);
    const auto *fsi_83 = buffer.data(fsi + 83);
    const auto *fsi_84 = buffer.data(fsi + 84);
    const auto *fsi_87 = buffer.data(fsi + 87);
    const auto *fsi_90 = buffer.data(fsi + 90);
    const auto *fsi_94 = buffer.data(fsi + 94);
    const auto *fsi_99 = buffer.data(fsi + 99);

    const auto *fsk1_0 = buffer.data(fsk1 + 0);
    const auto *fsk1_3 = buffer.data(fsk1 + 3);
    const auto *fsk1_5 = buffer.data(fsk1 + 5);
    const auto *fsk1_6 = buffer.data(fsk1 + 6);
    const auto *fsk1_9 = buffer.data(fsk1 + 9);
    const auto *fsk1_10 = buffer.data(fsk1 + 10);
    const auto *fsk1_14 = buffer.data(fsk1 + 14);
    const auto *fsk1_15 = buffer.data(fsk1 + 15);
    const auto *fsk1_20 = buffer.data(fsk1 + 20);
    const auto *fsk1_28 = buffer.data(fsk1 + 28);
    const auto *fsk1_35 = buffer.data(fsk1 + 35);

    const auto *gsh0_0 = buffer.data(gsh0 + 0);
    const auto *gsh0_1 = buffer.data(gsh0 + 1);
    const auto *gsh0_2 = buffer.data(gsh0 + 2);
    const auto *gsh0_3 = buffer.data(gsh0 + 3);
    const auto *gsh0_5 = buffer.data(gsh0 + 5);
    const auto *gsh0_6 = buffer.data(gsh0 + 6);
    const auto *gsh0_8 = buffer.data(gsh0 + 8);
    const auto *gsh0_9 = buffer.data(gsh0 + 9);
    const auto *gsh0_15 = buffer.data(gsh0 + 15);
    const auto *gsh0_17 = buffer.data(gsh0 + 17);
    const auto *gsh0_18 = buffer.data(gsh0 + 18);
    const auto *gsh0_19 = buffer.data(gsh0 + 19);
    const auto *gsh0_20 = buffer.data(gsh0 + 20);
    const auto *gsh0_24 = buffer.data(gsh0 + 24);
    const auto *gsh0_27 = buffer.data(gsh0 + 27);
    const auto *gsh0_28 = buffer.data(gsh0 + 28);
    const auto *gsh0_36 = buffer.data(gsh0 + 36);
    const auto *gsh0_37 = buffer.data(gsh0 + 37);
    const auto *gsh0_38 = buffer.data(gsh0 + 38);
    const auto *gsh0_39 = buffer.data(gsh0 + 39);
    const auto *gsh0_44 = buffer.data(gsh0 + 44);
    const auto *gsh0_46 = buffer.data(gsh0 + 46);
    const auto *gsh0_47 = buffer.data(gsh0 + 47);
    const auto *gsh0_49 = buffer.data(gsh0 + 49);
    const auto *gsh0_50 = buffer.data(gsh0 + 50);
    const auto *gsh0_51 = buffer.data(gsh0 + 51);
    const auto *gsh0_58 = buffer.data(gsh0 + 58);
    const auto *gsh0_59 = buffer.data(gsh0 + 59);
    const auto *gsh0_60 = buffer.data(gsh0 + 60);
    const auto *gsh0_61 = buffer.data(gsh0 + 61);
    const auto *gsh0_62 = buffer.data(gsh0 + 62);
    const auto *gsh0_63 = buffer.data(gsh0 + 63);
    const auto *gsh0_65 = buffer.data(gsh0 + 65);
    const auto *gsh0_66 = buffer.data(gsh0 + 66);
    const auto *gsh0_68 = buffer.data(gsh0 + 68);
    const auto *gsh0_69 = buffer.data(gsh0 + 69);
    const auto *gsh0_70 = buffer.data(gsh0 + 70);
    const auto *gsh0_72 = buffer.data(gsh0 + 72);
    const auto *gsh0_73 = buffer.data(gsh0 + 73);
    const auto *gsh0_78 = buffer.data(gsh0 + 78);

    const auto *gsh1_0 = buffer.data(gsh1 + 0);
    const auto *gsh1_1 = buffer.data(gsh1 + 1);
    const auto *gsh1_2 = buffer.data(gsh1 + 2);
    const auto *gsh1_3 = buffer.data(gsh1 + 3);
    const auto *gsh1_5 = buffer.data(gsh1 + 5);
    const auto *gsh1_6 = buffer.data(gsh1 + 6);
    const auto *gsh1_8 = buffer.data(gsh1 + 8);
    const auto *gsh1_9 = buffer.data(gsh1 + 9);
    const auto *gsh1_15 = buffer.data(gsh1 + 15);
    const auto *gsh1_17 = buffer.data(gsh1 + 17);
    const auto *gsh1_18 = buffer.data(gsh1 + 18);
    const auto *gsh1_19 = buffer.data(gsh1 + 19);
    const auto *gsh1_20 = buffer.data(gsh1 + 20);
    const auto *gsh1_24 = buffer.data(gsh1 + 24);
    const auto *gsh1_27 = buffer.data(gsh1 + 27);
    const auto *gsh1_28 = buffer.data(gsh1 + 28);
    const auto *gsh1_36 = buffer.data(gsh1 + 36);
    const auto *gsh1_37 = buffer.data(gsh1 + 37);
    const auto *gsh1_38 = buffer.data(gsh1 + 38);
    const auto *gsh1_39 = buffer.data(gsh1 + 39);
    const auto *gsh1_44 = buffer.data(gsh1 + 44);
    const auto *gsh1_46 = buffer.data(gsh1 + 46);
    const auto *gsh1_47 = buffer.data(gsh1 + 47);
    const auto *gsh1_49 = buffer.data(gsh1 + 49);
    const auto *gsh1_50 = buffer.data(gsh1 + 50);
    const auto *gsh1_51 = buffer.data(gsh1 + 51);
    const auto *gsh1_58 = buffer.data(gsh1 + 58);
    const auto *gsh1_59 = buffer.data(gsh1 + 59);
    const auto *gsh1_60 = buffer.data(gsh1 + 60);
    const auto *gsh1_61 = buffer.data(gsh1 + 61);
    const auto *gsh1_62 = buffer.data(gsh1 + 62);
    const auto *gsh1_63 = buffer.data(gsh1 + 63);
    const auto *gsh1_65 = buffer.data(gsh1 + 65);
    const auto *gsh1_66 = buffer.data(gsh1 + 66);
    const auto *gsh1_68 = buffer.data(gsh1 + 68);
    const auto *gsh1_69 = buffer.data(gsh1 + 69);
    const auto *gsh1_70 = buffer.data(gsh1 + 70);
    const auto *gsh1_72 = buffer.data(gsh1 + 72);
    const auto *gsh1_73 = buffer.data(gsh1 + 73);
    const auto *gsh1_78 = buffer.data(gsh1 + 78);

    const auto *gsi_0 = buffer.data(gsi + 0);
    const auto *gsi_1 = buffer.data(gsi + 1);
    const auto *gsi_2 = buffer.data(gsi + 2);
    const auto *gsi_3 = buffer.data(gsi + 3);
    const auto *gsi_5 = buffer.data(gsi + 5);
    const auto *gsi_6 = buffer.data(gsi + 6);
    const auto *gsi_8 = buffer.data(gsi + 8);
    const auto *gsi_9 = buffer.data(gsi + 9);
    const auto *gsi_10 = buffer.data(gsi + 10);
    const auto *gsi_12 = buffer.data(gsi + 12);
    const auto *gsi_13 = buffer.data(gsi + 13);
    const auto *gsi_14 = buffer.data(gsi + 14);
    const auto *gsi_15 = buffer.data(gsi + 15);
    const auto *gsi_20 = buffer.data(gsi + 20);
    const auto *gsi_21 = buffer.data(gsi + 21);
    const auto *gsi_23 = buffer.data(gsi + 23);
    const auto *gsi_24 = buffer.data(gsi + 24);
    const auto *gsi_25 = buffer.data(gsi + 25);
    const auto *gsi_26 = buffer.data(gsi + 26);
    const auto *gsi_27 = buffer.data(gsi + 27);
    const auto *gsi_28 = buffer.data(gsi + 28);
    const auto *gsi_29 = buffer.data(gsi + 29);
    const auto *gsi_31 = buffer.data(gsi + 31);
    const auto *gsi_33 = buffer.data(gsi + 33);
    const auto *gsi_34 = buffer.data(gsi + 34);
    const auto *gsi_35 = buffer.data(gsi + 35);
    const auto *gsi_37 = buffer.data(gsi + 37);
    const auto *gsi_38 = buffer.data(gsi + 38);
    const auto *gsi_39 = buffer.data(gsi + 39);
    const auto *gsi_40 = buffer.data(gsi + 40);
    const auto *gsi_42 = buffer.data(gsi + 42);
    const auto *gsi_43 = buffer.data(gsi + 43);
    const auto *gsi_49 = buffer.data(gsi + 49);
    const auto *gsi_50 = buffer.data(gsi + 50);
    const auto *gsi_51 = buffer.data(gsi + 51);
    const auto *gsi_52 = buffer.data(gsi + 52);
    const auto *gsi_53 = buffer.data(gsi + 53);
    const auto *gsi_54 = buffer.data(gsi + 54);
    const auto *gsi_55 = buffer.data(gsi + 55);
    const auto *gsi_56 = buffer.data(gsi + 56);
    const auto *gsi_58 = buffer.data(gsi + 58);
    const auto *gsi_60 = buffer.data(gsi + 60);
    const auto *gsi_61 = buffer.data(gsi + 61);
    const auto *gsi_63 = buffer.data(gsi + 63);
    const auto *gsi_64 = buffer.data(gsi + 64);
    const auto *gsi_65 = buffer.data(gsi + 65);
    const auto *gsi_67 = buffer.data(gsi + 67);
    const auto *gsi_68 = buffer.data(gsi + 68);
    const auto *gsi_69 = buffer.data(gsi + 69);
    const auto *gsi_70 = buffer.data(gsi + 70);
    const auto *gsi_76 = buffer.data(gsi + 76);
    const auto *gsi_77 = buffer.data(gsi + 77);
    const auto *gsi_78 = buffer.data(gsi + 78);
    const auto *gsi_79 = buffer.data(gsi + 79);
    const auto *gsi_80 = buffer.data(gsi + 80);
    const auto *gsi_81 = buffer.data(gsi + 81);
    const auto *gsi_82 = buffer.data(gsi + 82);
    const auto *gsi_83 = buffer.data(gsi + 83);
    const auto *gsi_84 = buffer.data(gsi + 84);
    const auto *gsi_85 = buffer.data(gsi + 85);
    const auto *gsi_86 = buffer.data(gsi + 86);
    const auto *gsi_87 = buffer.data(gsi + 87);
    const auto *gsi_89 = buffer.data(gsi + 89);
    const auto *gsi_90 = buffer.data(gsi + 90);
    const auto *gsi_91 = buffer.data(gsi + 91);
    const auto *gsi_93 = buffer.data(gsi + 93);
    const auto *gsi_94 = buffer.data(gsi + 94);
    const auto *gsi_95 = buffer.data(gsi + 95);
    const auto *gsi_96 = buffer.data(gsi + 96);
    const auto *gsi_98 = buffer.data(gsi + 98);
    const auto *gsi_99 = buffer.data(gsi + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fsi_0, gsh0_0, \
                         gsh1_0, gsi_0, gsi_1, gsi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fsi_0[k]
                 + f_1 * gsh0_0[k]
                 - f_2 * gsh1_0[k]
                 + f_3 * pc_x[k] * gsi_0[k];

        t_1[k] = f_3 * pc_y[k] * gsi_0[k];

        t_2[k] = f_3 * pc_z[k] * gsi_0[k];

        t_3[k] = f_4 * gsh0_0[k]
                 - f_5 * gsh1_0[k]
                 + f_3 * pc_y[k] * gsi_1[k];

        t_4[k] = f_3 * pc_y[k] * gsi_2[k];

        t_5[k] = f_4 * gsh0_0[k]
                 - f_5 * gsh1_0[k]
                 + f_3 * pc_z[k] * gsi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, gsh0_1, gsh0_2, gsh0_3, gsh1_1, \
                         gsh1_2, gsh1_3, gsi_3, gsi_5, gsi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * gsh0_1[k]
                 - f_7 * gsh1_1[k]
                 + f_3 * pc_y[k] * gsi_3[k];

        t_7[k] = f_3 * pc_z[k] * gsi_3[k];

        t_8[k] = f_3 * pc_y[k] * gsi_5[k];

        t_9[k] = f_6 * gsh0_2[k]
                 - f_7 * gsh1_2[k]
                 + f_3 * pc_z[k] * gsi_5[k];

        t_10[k] = f_8 * gsh0_3[k]
                  - f_9 * gsh1_3[k]
                  + f_3 * pc_y[k] * gsi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, gsh0_5, gsh0_6, \
                         gsh1_5, gsh1_6, gsi_6, gsi_8, gsi_9, gsi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * gsi_6[k];

        t_12[k] = f_4 * gsh0_5[k]
                  - f_5 * gsh1_5[k]
                  + f_3 * pc_y[k] * gsi_8[k];

        t_13[k] = f_3 * pc_y[k] * gsi_9[k];

        t_14[k] = f_8 * gsh0_5[k]
                  - f_9 * gsh1_5[k]
                  + f_3 * pc_z[k] * gsi_9[k];

        t_15[k] = f_10 * gsh0_6[k]
                  - f_11 * gsh1_6[k]
                  + f_3 * pc_y[k] * gsi_10[k];

        t_16[k] = f_3 * pc_z[k] * gsi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, gsh0_8, gsh0_9, gsh1_8, gsh1_9, \
                         gsi_12, gsi_13, gsi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * gsh0_8[k]
                  - f_7 * gsh1_8[k]
                  + f_3 * pc_y[k] * gsi_12[k];

        t_18[k] = f_4 * gsh0_9[k]
                  - f_5 * gsh1_9[k]
                  + f_3 * pc_y[k] * gsi_13[k];

        t_19[k] = f_3 * pc_y[k] * gsi_14[k];

        t_20[k] = f_10 * gsh0_9[k]
                  - f_11 * gsh1_9[k]
                  + f_3 * pc_z[k] * gsi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, fsi_21, fsi_23, fsi_24, \
                         fsi_25, gsi_15, gsi_21, gsi_23, gsi_24, \
                         gsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * fsi_21[k]
                  + f_3 * pc_x[k] * gsi_21[k];

        t_22[k] = f_3 * pc_z[k] * gsi_15[k];

        t_23[k] = f_0 * fsi_23[k]
                  + f_3 * pc_x[k] * gsi_23[k];

        t_24[k] = f_0 * fsi_24[k]
                  + f_3 * pc_x[k] * gsi_24[k];

        t_25[k] = f_0 * fsi_25[k]
                  + f_3 * pc_x[k] * gsi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, fsi_27, gsh0_15, gsh1_15, \
                         gsi_20, gsi_21, gsi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * gsi_20[k];

        t_27[k] = f_0 * fsi_27[k]
                  + f_3 * pc_x[k] * gsi_27[k];

        t_28[k] = f_1 * gsh0_15[k]
                  - f_2 * gsh1_15[k]
                  + f_3 * pc_y[k] * gsi_21[k];

        t_29[k] = f_3 * pc_z[k] * gsi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, gsh0_17, gsh0_18, gsh0_19, gsh1_17, gsh1_18, \
                         gsh1_19, gsi_23, gsi_24, gsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * gsh0_17[k]
                  - f_11 * gsh1_17[k]
                  + f_3 * pc_y[k] * gsi_23[k];

        t_31[k] = f_8 * gsh0_18[k]
                  - f_9 * gsh1_18[k]
                  + f_3 * pc_y[k] * gsi_24[k];

        t_32[k] = f_6 * gsh0_19[k]
                  - f_7 * gsh1_19[k]
                  + f_3 * pc_y[k] * gsi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, fsk0_0, fsi_0, \
                         fsk1_0, gsh0_20, gsh1_20, gsi_26, gsi_27, \
                         gsi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * gsh0_20[k]
                  - f_5 * gsh1_20[k]
                  + f_3 * pc_y[k] * gsi_26[k];

        t_34[k] = f_3 * pc_y[k] * gsi_27[k];

        t_35[k] = f_1 * gsh0_20[k]
                  - f_2 * gsh1_20[k]
                  + f_3 * pc_z[k] * gsi_27[k];

        t_36[k] = pa_y[k] * fsk0_0[k]
                  - f_12 * pc_y[k] * fsk1_0[k];

        t_37[k] = f_13 * fsi_0[k]
                  + f_3 * pc_y[k] * gsi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, fsk0_3, fsk0_5, fsi_1, \
                         fsk1_3, fsk1_5, gsi_28, gsi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * gsi_28[k];

        t_39[k] = pa_y[k] * fsk0_3[k]
                  + f_14 * fsi_1[k]
                  - f_12 * pc_y[k] * fsk1_3[k];

        t_40[k] = f_3 * pc_z[k] * gsi_29[k];

        t_41[k] = pa_y[k] * fsk0_5[k]
                  - f_12 * pc_y[k] * fsk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_y, pc_z, fsk0_6, fsk0_9, fsi_3, \
                         fsi_5, fsk1_6, fsk1_9, gsi_31, gsi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * fsk0_6[k]
                  + f_15 * fsi_3[k]
                  - f_12 * pc_y[k] * fsk1_6[k];

        t_43[k] = f_3 * pc_z[k] * gsi_31[k];

        t_44[k] = f_13 * fsi_5[k]
                  + f_3 * pc_y[k] * gsi_33[k];

        t_45[k] = pa_y[k] * fsk0_9[k]
                  - f_12 * pc_y[k] * fsk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, fsk0_10, fsi_6, fsi_9, \
                         fsk1_10, gsh0_24, gsh1_24, gsi_34, gsi_35, \
                         gsi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * fsk0_10[k]
                  + f_0 * fsi_6[k]
                  - f_12 * pc_y[k] * fsk1_10[k];

        t_47[k] = f_3 * pc_z[k] * gsi_34[k];

        t_48[k] = f_4 * gsh0_24[k]
                  - f_5 * gsh1_24[k]
                  + f_3 * pc_z[k] * gsi_35[k];

        t_49[k] = f_13 * fsi_9[k]
                  + f_3 * pc_y[k] * gsi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_y, pc_z, fsk0_14, fsk0_15, fsi_10, \
                         fsk1_14, fsk1_15, gsh0_27, gsh1_27, gsi_38, \
                         gsi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * fsk0_14[k]
                  - f_12 * pc_y[k] * fsk1_14[k];

        t_51[k] = pa_y[k] * fsk0_15[k]
                  + f_16 * fsi_10[k]
                  - f_12 * pc_y[k] * fsk1_15[k];

        t_52[k] = f_3 * pc_z[k] * gsi_38[k];

        t_53[k] = f_4 * gsh0_27[k]
                  - f_5 * gsh1_27[k]
                  + f_3 * pc_z[k] * gsi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, fsk0_20, fsi_14, fsk1_20, \
                         gsh0_28, gsh1_28, gsi_40, gsi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * gsh0_28[k]
                  - f_7 * gsh1_28[k]
                  + f_3 * pc_z[k] * gsi_40[k];

        t_55[k] = f_13 * fsi_14[k]
                  + f_3 * pc_y[k] * gsi_42[k];

        t_56[k] = pa_y[k] * fsk0_20[k]
                  - f_12 * pc_y[k] * fsk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, fsi_49, fsi_51, fsi_52, \
                         fsi_53, gsi_43, gsi_49, gsi_51, gsi_52, \
                         gsi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_15 * fsi_49[k]
                  + f_3 * pc_x[k] * gsi_49[k];

        t_58[k] = f_3 * pc_z[k] * gsi_43[k];

        t_59[k] = f_15 * fsi_51[k]
                  + f_3 * pc_x[k] * gsi_51[k];

        t_60[k] = f_15 * fsi_52[k]
                  + f_3 * pc_x[k] * gsi_52[k];

        t_61[k] = f_15 * fsi_53[k]
                  + f_3 * pc_x[k] * gsi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, fsi_21, fsi_54, fsi_55, \
                         gsh0_36, gsh1_36, gsi_49, gsi_54, gsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_15 * fsi_54[k]
                  + f_3 * pc_x[k] * gsi_54[k];

        t_63[k] = f_15 * fsi_55[k]
                  + f_3 * pc_x[k] * gsi_55[k];

        t_64[k] = f_13 * fsi_21[k]
                  + f_1 * gsh0_36[k]
                  - f_2 * gsh1_36[k]
                  + f_3 * pc_y[k] * gsi_49[k];

        t_65[k] = f_3 * pc_z[k] * gsi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_z, gsh0_36, gsh0_37, gsh0_38, gsh1_36, gsh1_37, \
                         gsh1_38, gsi_50, gsi_51, gsi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * gsh0_36[k]
                  - f_5 * gsh1_36[k]
                  + f_3 * pc_z[k] * gsi_50[k];

        t_67[k] = f_6 * gsh0_37[k]
                  - f_7 * gsh1_37[k]
                  + f_3 * pc_z[k] * gsi_51[k];

        t_68[k] = f_8 * gsh0_38[k]
                  - f_9 * gsh1_38[k]
                  + f_3 * pc_z[k] * gsi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, fsk0_35, fsi_27, fsk1_35, \
                         gsh0_39, gsh1_39, gsi_53, gsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * gsh0_39[k]
                  - f_11 * gsh1_39[k]
                  + f_3 * pc_z[k] * gsi_53[k];

        t_70[k] = f_13 * fsi_27[k]
                  + f_3 * pc_y[k] * gsi_55[k];

        t_71[k] = pa_y[k] * fsk0_35[k]
                  - f_12 * pc_y[k] * fsk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pc_y, pc_z, fsk0_0, fsk0_3, \
                         fsi_0, fsk1_0, fsk1_3, gsi_56, gsi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * fsk0_0[k]
                  - f_12 * pc_z[k] * fsk1_0[k];

        t_73[k] = f_3 * pc_y[k] * gsi_56[k];

        t_74[k] = f_13 * fsi_0[k]
                  + f_3 * pc_z[k] * gsi_56[k];

        t_75[k] = pa_z[k] * fsk0_3[k]
                  - f_12 * pc_z[k] * fsk1_3[k];

        t_76[k] = f_3 * pc_y[k] * gsi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pc_y, pc_z, fsk0_5, fsk0_6, fsi_2, \
                         fsk1_5, fsk1_6, gsh0_44, gsh1_44, gsi_60, \
                         gsi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * fsk0_5[k]
                  + f_14 * fsi_2[k]
                  - f_12 * pc_z[k] * fsk1_5[k];

        t_78[k] = pa_z[k] * fsk0_6[k]
                  - f_12 * pc_z[k] * fsk1_6[k];

        t_79[k] = f_4 * gsh0_44[k]
                  - f_5 * gsh1_44[k]
                  + f_3 * pc_y[k] * gsi_60[k];

        t_80[k] = f_3 * pc_y[k] * gsi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, fsk0_9, fsk0_10, fsi_5, fsk1_9, \
                         fsk1_10, gsh0_46, gsh1_46, gsi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * fsk0_9[k]
                  + f_15 * fsi_5[k]
                  - f_12 * pc_z[k] * fsk1_9[k];

        t_82[k] = pa_z[k] * fsk0_10[k]
                  - f_12 * pc_z[k] * fsk1_10[k];

        t_83[k] = f_6 * gsh0_46[k]
                  - f_7 * gsh1_46[k]
                  + f_3 * pc_y[k] * gsi_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pc_y, pc_z, fsk0_14, fsk0_15, fsi_9, \
                         fsk1_14, fsk1_15, gsh0_47, gsh1_47, gsi_64, \
                         gsi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * gsh0_47[k]
                  - f_5 * gsh1_47[k]
                  + f_3 * pc_y[k] * gsi_64[k];

        t_85[k] = f_3 * pc_y[k] * gsi_65[k];

        t_86[k] = pa_z[k] * fsk0_14[k]
                  + f_0 * fsi_9[k]
                  - f_12 * pc_z[k] * fsk1_14[k];

        t_87[k] = pa_z[k] * fsk0_15[k]
                  - f_12 * pc_z[k] * fsk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, gsh0_49, gsh0_50, gsh0_51, gsh1_49, \
                         gsh1_50, gsh1_51, gsi_67, gsi_68, gsi_69, \
                         gsi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_8 * gsh0_49[k]
                  - f_9 * gsh1_49[k]
                  + f_3 * pc_y[k] * gsi_67[k];

        t_89[k] = f_6 * gsh0_50[k]
                  - f_7 * gsh1_50[k]
                  + f_3 * pc_y[k] * gsi_68[k];

        t_90[k] = f_4 * gsh0_51[k]
                  - f_5 * gsh1_51[k]
                  + f_3 * pc_y[k] * gsi_69[k];

        t_91[k] = f_3 * pc_y[k] * gsi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_z, fsk0_20, fsi_14, fsi_77, \
                         fsi_78, fsi_79, fsk1_20, gsi_77, gsi_78, \
                         gsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * fsk0_20[k]
                  + f_16 * fsi_14[k]
                  - f_12 * pc_z[k] * fsk1_20[k];

        t_93[k] = f_15 * fsi_77[k]
                  + f_3 * pc_x[k] * gsi_77[k];

        t_94[k] = f_15 * fsi_78[k]
                  + f_3 * pc_x[k] * gsi_78[k];

        t_95[k] = f_15 * fsi_79[k]
                  + f_3 * pc_x[k] * gsi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, fsi_80, fsi_81, fsi_83, gsi_76, \
                         gsi_80, gsi_81, gsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_15 * fsi_80[k]
                  + f_3 * pc_x[k] * gsi_80[k];

        t_97[k] = f_15 * fsi_81[k]
                  + f_3 * pc_x[k] * gsi_81[k];

        t_98[k] = f_3 * pc_y[k] * gsi_76[k];

        t_99[k] = f_15 * fsi_83[k]
                  + f_3 * pc_x[k] * gsi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_y, pc_z, fsk0_28, fsk1_28, gsh0_58, \
                         gsh0_59, gsh1_58, gsh1_59, gsi_78, gsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * fsk0_28[k]
                   - f_12 * pc_z[k] * fsk1_28[k];

        t_101[k] = f_17 * gsh0_58[k]
                   - f_18 * gsh1_58[k]
                   + f_3 * pc_y[k] * gsi_78[k];

        t_102[k] = f_10 * gsh0_59[k]
                   - f_11 * gsh1_59[k]
                   + f_3 * pc_y[k] * gsi_79[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_y, gsh0_60, gsh0_61, gsh0_62, gsh1_60, \
                         gsh1_61, gsh1_62, gsi_80, gsi_81, gsi_82, \
                         gsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * gsh0_60[k]
                   - f_9 * gsh1_60[k]
                   + f_3 * pc_y[k] * gsi_80[k];

        t_104[k] = f_6 * gsh0_61[k]
                   - f_7 * gsh1_61[k]
                   + f_3 * pc_y[k] * gsi_81[k];

        t_105[k] = f_4 * gsh0_62[k]
                   - f_5 * gsh1_62[k]
                   + f_3 * pc_y[k] * gsi_82[k];

        t_106[k] = f_3 * pc_y[k] * gsi_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, pc_z, fsi_27, fsi_28, fsi_84, \
                         gsh0_62, gsh0_63, gsh1_62, gsh1_63, gsi_83, \
                         gsi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_13 * fsi_27[k]
                   + f_1 * gsh0_62[k]
                   - f_2 * gsh1_62[k]
                   + f_3 * pc_z[k] * gsi_83[k];

        t_108[k] = f_14 * fsi_84[k]
                   + f_1 * gsh0_63[k]
                   - f_2 * gsh1_63[k]
                   + f_3 * pc_x[k] * gsi_84[k];

        t_109[k] = f_14 * fsi_28[k]
                   + f_3 * pc_y[k] * gsi_84[k];

        t_110[k] = f_3 * pc_z[k] * gsi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_z, fsi_87, gsh0_63, gsh0_66, gsh1_63, \
                         gsh1_66, gsi_85, gsi_86, gsi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_14 * fsi_87[k]
                   + f_10 * gsh0_66[k]
                   - f_11 * gsh1_66[k]
                   + f_3 * pc_x[k] * gsi_87[k];

        t_112[k] = f_3 * pc_z[k] * gsi_85[k];

        t_113[k] = f_4 * gsh0_63[k]
                   - f_5 * gsh1_63[k]
                   + f_3 * pc_z[k] * gsi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, pc_y, pc_z, fsi_33, fsi_90, \
                         gsh0_65, gsh0_69, gsh1_65, gsh1_69, gsi_87, gsi_89, \
                         gsi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_14 * fsi_90[k]
                   + f_8 * gsh0_69[k]
                   - f_9 * gsh1_69[k]
                   + f_3 * pc_x[k] * gsi_90[k];

        t_115[k] = f_3 * pc_z[k] * gsi_87[k];

        t_116[k] = f_14 * fsi_33[k]
                   + f_3 * pc_y[k] * gsi_89[k];

        t_117[k] = f_6 * gsh0_65[k]
                   - f_7 * gsh1_65[k]
                   + f_3 * pc_z[k] * gsi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_z, fsi_94, gsh0_66, gsh0_73, gsh1_66, \
                         gsh1_73, gsi_90, gsi_91, gsi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_14 * fsi_94[k]
                   + f_6 * gsh0_73[k]
                   - f_7 * gsh1_73[k]
                   + f_3 * pc_x[k] * gsi_94[k];

        t_119[k] = f_3 * pc_z[k] * gsi_90[k];

        t_120[k] = f_4 * gsh0_66[k]
                   - f_5 * gsh1_66[k]
                   + f_3 * pc_z[k] * gsi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, fsi_37, fsi_99, \
                         gsh0_68, gsh0_78, gsh1_68, gsh1_78, gsi_93, gsi_94, \
                         gsi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_14 * fsi_37[k]
                   + f_3 * pc_y[k] * gsi_93[k];

        t_122[k] = f_8 * gsh0_68[k]
                   - f_9 * gsh1_68[k]
                   + f_3 * pc_z[k] * gsi_93[k];

        t_123[k] = f_14 * fsi_99[k]
                   + f_4 * gsh0_78[k]
                   - f_5 * gsh1_78[k]
                   + f_3 * pc_x[k] * gsi_99[k];

        t_124[k] = f_3 * pc_z[k] * gsi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_y, pc_z, fsi_42, gsh0_69, gsh0_70, \
                         gsh0_72, gsh1_69, gsh1_70, gsh1_72, gsi_95, gsi_96, \
                         gsi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * gsh0_69[k]
                   - f_5 * gsh1_69[k]
                   + f_3 * pc_z[k] * gsi_95[k];

        t_126[k] = f_6 * gsh0_70[k]
                   - f_7 * gsh1_70[k]
                   + f_3 * pc_z[k] * gsi_96[k];

        t_127[k] = f_14 * fsi_42[k]
                   + f_3 * pc_y[k] * gsi_98[k];

        t_128[k] = f_10 * gsh0_72[k]
                   - f_11 * gsh1_72[k]
                   + f_3 * pc_z[k] * gsi_98[k];
    }
}

static auto
compute_prim_gsk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsk0,
                                                          const size_t fsi, const size_t fsk1,
                                                          const size_t gsh0, const size_t gsh1,
                                                          const size_t gsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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
    auto *t_250 = buffer.data(target + 250);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsk0_39 = buffer.data(fsk0 + 39);
    const auto *fsk0_42 = buffer.data(fsk0 + 42);
    const auto *fsk0_46 = buffer.data(fsk0 + 46);
    const auto *fsk0_51 = buffer.data(fsk0 + 51);
    const auto *fsk0_64 = buffer.data(fsk0 + 64);
    const auto *fsk0_72 = buffer.data(fsk0 + 72);
    const auto *fsk0_77 = buffer.data(fsk0 + 77);
    const auto *fsk0_81 = buffer.data(fsk0 + 81);
    const auto *fsk0_84 = buffer.data(fsk0 + 84);
    const auto *fsk0_86 = buffer.data(fsk0 + 86);
    const auto *fsk0_89 = buffer.data(fsk0 + 89);
    const auto *fsk0_90 = buffer.data(fsk0 + 90);
    const auto *fsk0_92 = buffer.data(fsk0 + 92);
    const auto *fsk0_107 = buffer.data(fsk0 + 107);
    const auto *fsk0_216 = buffer.data(fsk0 + 216);
    const auto *fsk0_219 = buffer.data(fsk0 + 219);
    const auto *fsk0_222 = buffer.data(fsk0 + 222);
    const auto *fsk0_226 = buffer.data(fsk0 + 226);
    const auto *fsk0_231 = buffer.data(fsk0 + 231);
    const auto *fsk0_244 = buffer.data(fsk0 + 244);
    const auto *fsk0_246 = buffer.data(fsk0 + 246);
    const auto *fsk0_247 = buffer.data(fsk0 + 247);
    const auto *fsk0_248 = buffer.data(fsk0 + 248);
    const auto *fsk0_249 = buffer.data(fsk0 + 249);

    const auto *fsi_28 = buffer.data(fsi + 28);
    const auto *fsi_31 = buffer.data(fsi + 31);
    const auto *fsi_34 = buffer.data(fsi + 34);
    const auto *fsi_38 = buffer.data(fsi + 38);
    const auto *fsi_49 = buffer.data(fsi + 49);
    const auto *fsi_55 = buffer.data(fsi + 55);
    const auto *fsi_56 = buffer.data(fsi + 56);
    const auto *fsi_58 = buffer.data(fsi + 58);
    const auto *fsi_61 = buffer.data(fsi + 61);
    const auto *fsi_64 = buffer.data(fsi + 64);
    const auto *fsi_65 = buffer.data(fsi + 65);
    const auto *fsi_68 = buffer.data(fsi + 68);
    const auto *fsi_69 = buffer.data(fsi + 69);
    const auto *fsi_70 = buffer.data(fsi + 70);
    const auto *fsi_79 = buffer.data(fsi + 79);
    const auto *fsi_80 = buffer.data(fsi + 80);
    const auto *fsi_81 = buffer.data(fsi + 81);
    const auto *fsi_82 = buffer.data(fsi + 82);
    const auto *fsi_83 = buffer.data(fsi + 83);
    const auto *fsi_84 = buffer.data(fsi + 84);
    const auto *fsi_89 = buffer.data(fsi + 89);
    const auto *fsi_93 = buffer.data(fsi + 93);
    const auto *fsi_98 = buffer.data(fsi + 98);
    const auto *fsi_105 = buffer.data(fsi + 105);
    const auto *fsi_107 = buffer.data(fsi + 107);
    const auto *fsi_108 = buffer.data(fsi + 108);
    const auto *fsi_109 = buffer.data(fsi + 109);
    const auto *fsi_110 = buffer.data(fsi + 110);
    const auto *fsi_111 = buffer.data(fsi + 111);
    const auto *fsi_133 = buffer.data(fsi + 133);
    const auto *fsi_134 = buffer.data(fsi + 134);
    const auto *fsi_135 = buffer.data(fsi + 135);
    const auto *fsi_136 = buffer.data(fsi + 136);
    const auto *fsi_137 = buffer.data(fsi + 137);
    const auto *fsi_138 = buffer.data(fsi + 138);
    const auto *fsi_139 = buffer.data(fsi + 139);
    const auto *fsi_140 = buffer.data(fsi + 140);
    const auto *fsi_145 = buffer.data(fsi + 145);
    const auto *fsi_149 = buffer.data(fsi + 149);
    const auto *fsi_154 = buffer.data(fsi + 154);
    const auto *fsi_160 = buffer.data(fsi + 160);
    const auto *fsi_161 = buffer.data(fsi + 161);
    const auto *fsi_162 = buffer.data(fsi + 162);
    const auto *fsi_163 = buffer.data(fsi + 163);
    const auto *fsi_164 = buffer.data(fsi + 164);
    const auto *fsi_165 = buffer.data(fsi + 165);
    const auto *fsi_167 = buffer.data(fsi + 167);
    const auto *fsi_168 = buffer.data(fsi + 168);
    const auto *fsi_171 = buffer.data(fsi + 171);
    const auto *fsi_174 = buffer.data(fsi + 174);
    const auto *fsi_178 = buffer.data(fsi + 178);
    const auto *fsi_183 = buffer.data(fsi + 183);
    const auto *fsi_189 = buffer.data(fsi + 189);
    const auto *fsi_191 = buffer.data(fsi + 191);
    const auto *fsi_192 = buffer.data(fsi + 192);
    const auto *fsi_193 = buffer.data(fsi + 193);
    const auto *fsi_194 = buffer.data(fsi + 194);
    const auto *fsi_195 = buffer.data(fsi + 195);

    const auto *fsk1_39 = buffer.data(fsk1 + 39);
    const auto *fsk1_42 = buffer.data(fsk1 + 42);
    const auto *fsk1_46 = buffer.data(fsk1 + 46);
    const auto *fsk1_51 = buffer.data(fsk1 + 51);
    const auto *fsk1_64 = buffer.data(fsk1 + 64);
    const auto *fsk1_72 = buffer.data(fsk1 + 72);
    const auto *fsk1_77 = buffer.data(fsk1 + 77);
    const auto *fsk1_81 = buffer.data(fsk1 + 81);
    const auto *fsk1_84 = buffer.data(fsk1 + 84);
    const auto *fsk1_86 = buffer.data(fsk1 + 86);
    const auto *fsk1_89 = buffer.data(fsk1 + 89);
    const auto *fsk1_90 = buffer.data(fsk1 + 90);
    const auto *fsk1_92 = buffer.data(fsk1 + 92);
    const auto *fsk1_107 = buffer.data(fsk1 + 107);
    const auto *fsk1_216 = buffer.data(fsk1 + 216);
    const auto *fsk1_219 = buffer.data(fsk1 + 219);
    const auto *fsk1_222 = buffer.data(fsk1 + 222);
    const auto *fsk1_226 = buffer.data(fsk1 + 226);
    const auto *fsk1_231 = buffer.data(fsk1 + 231);
    const auto *fsk1_244 = buffer.data(fsk1 + 244);
    const auto *fsk1_246 = buffer.data(fsk1 + 246);
    const auto *fsk1_247 = buffer.data(fsk1 + 247);
    const auto *fsk1_248 = buffer.data(fsk1 + 248);
    const auto *fsk1_249 = buffer.data(fsk1 + 249);

    const auto *gsh0_78 = buffer.data(gsh0 + 78);
    const auto *gsh0_79 = buffer.data(gsh0 + 79);
    const auto *gsh0_80 = buffer.data(gsh0 + 80);
    const auto *gsh0_81 = buffer.data(gsh0 + 81);
    const auto *gsh0_83 = buffer.data(gsh0 + 83);
    const auto *gsh0_101 = buffer.data(gsh0 + 101);
    const auto *gsh0_102 = buffer.data(gsh0 + 102);
    const auto *gsh0_103 = buffer.data(gsh0 + 103);
    const auto *gsh0_104 = buffer.data(gsh0 + 104);
    const auto *gsh0_105 = buffer.data(gsh0 + 105);
    const auto *gsh0_106 = buffer.data(gsh0 + 106);
    const auto *gsh0_107 = buffer.data(gsh0 + 107);
    const auto *gsh0_108 = buffer.data(gsh0 + 108);
    const auto *gsh0_109 = buffer.data(gsh0 + 109);
    const auto *gsh0_110 = buffer.data(gsh0 + 110);
    const auto *gsh0_111 = buffer.data(gsh0 + 111);
    const auto *gsh0_112 = buffer.data(gsh0 + 112);
    const auto *gsh0_113 = buffer.data(gsh0 + 113);
    const auto *gsh0_114 = buffer.data(gsh0 + 114);
    const auto *gsh0_119 = buffer.data(gsh0 + 119);
    const auto *gsh0_120 = buffer.data(gsh0 + 120);
    const auto *gsh0_121 = buffer.data(gsh0 + 121);
    const auto *gsh0_122 = buffer.data(gsh0 + 122);
    const auto *gsh0_123 = buffer.data(gsh0 + 123);
    const auto *gsh0_124 = buffer.data(gsh0 + 124);
    const auto *gsh0_125 = buffer.data(gsh0 + 125);
    const auto *gsh0_126 = buffer.data(gsh0 + 126);
    const auto *gsh0_128 = buffer.data(gsh0 + 128);
    const auto *gsh0_129 = buffer.data(gsh0 + 129);
    const auto *gsh0_131 = buffer.data(gsh0 + 131);
    const auto *gsh0_132 = buffer.data(gsh0 + 132);
    const auto *gsh0_133 = buffer.data(gsh0 + 133);
    const auto *gsh0_135 = buffer.data(gsh0 + 135);

    const auto *gsh1_78 = buffer.data(gsh1 + 78);
    const auto *gsh1_79 = buffer.data(gsh1 + 79);
    const auto *gsh1_80 = buffer.data(gsh1 + 80);
    const auto *gsh1_81 = buffer.data(gsh1 + 81);
    const auto *gsh1_83 = buffer.data(gsh1 + 83);
    const auto *gsh1_101 = buffer.data(gsh1 + 101);
    const auto *gsh1_102 = buffer.data(gsh1 + 102);
    const auto *gsh1_103 = buffer.data(gsh1 + 103);
    const auto *gsh1_104 = buffer.data(gsh1 + 104);
    const auto *gsh1_105 = buffer.data(gsh1 + 105);
    const auto *gsh1_106 = buffer.data(gsh1 + 106);
    const auto *gsh1_107 = buffer.data(gsh1 + 107);
    const auto *gsh1_108 = buffer.data(gsh1 + 108);
    const auto *gsh1_109 = buffer.data(gsh1 + 109);
    const auto *gsh1_110 = buffer.data(gsh1 + 110);
    const auto *gsh1_111 = buffer.data(gsh1 + 111);
    const auto *gsh1_112 = buffer.data(gsh1 + 112);
    const auto *gsh1_113 = buffer.data(gsh1 + 113);
    const auto *gsh1_114 = buffer.data(gsh1 + 114);
    const auto *gsh1_119 = buffer.data(gsh1 + 119);
    const auto *gsh1_120 = buffer.data(gsh1 + 120);
    const auto *gsh1_121 = buffer.data(gsh1 + 121);
    const auto *gsh1_122 = buffer.data(gsh1 + 122);
    const auto *gsh1_123 = buffer.data(gsh1 + 123);
    const auto *gsh1_124 = buffer.data(gsh1 + 124);
    const auto *gsh1_125 = buffer.data(gsh1 + 125);
    const auto *gsh1_126 = buffer.data(gsh1 + 126);
    const auto *gsh1_128 = buffer.data(gsh1 + 128);
    const auto *gsh1_129 = buffer.data(gsh1 + 129);
    const auto *gsh1_131 = buffer.data(gsh1 + 131);
    const auto *gsh1_132 = buffer.data(gsh1 + 132);
    const auto *gsh1_133 = buffer.data(gsh1 + 133);
    const auto *gsh1_135 = buffer.data(gsh1 + 135);

    const auto *gsi_99 = buffer.data(gsi + 99);
    const auto *gsi_105 = buffer.data(gsi + 105);
    const auto *gsi_106 = buffer.data(gsi + 106);
    const auto *gsi_107 = buffer.data(gsi + 107);
    const auto *gsi_108 = buffer.data(gsi + 108);
    const auto *gsi_109 = buffer.data(gsi + 109);
    const auto *gsi_110 = buffer.data(gsi + 110);
    const auto *gsi_111 = buffer.data(gsi + 111);
    const auto *gsi_112 = buffer.data(gsi + 112);
    const auto *gsi_114 = buffer.data(gsi + 114);
    const auto *gsi_115 = buffer.data(gsi + 115);
    const auto *gsi_117 = buffer.data(gsi + 117);
    const auto *gsi_118 = buffer.data(gsi + 118);
    const auto *gsi_121 = buffer.data(gsi + 121);
    const auto *gsi_122 = buffer.data(gsi + 122);
    const auto *gsi_126 = buffer.data(gsi + 126);
    const auto *gsi_133 = buffer.data(gsi + 133);
    const auto *gsi_134 = buffer.data(gsi + 134);
    const auto *gsi_135 = buffer.data(gsi + 135);
    const auto *gsi_136 = buffer.data(gsi + 136);
    const auto *gsi_137 = buffer.data(gsi + 137);
    const auto *gsi_138 = buffer.data(gsi + 138);
    const auto *gsi_139 = buffer.data(gsi + 139);
    const auto *gsi_140 = buffer.data(gsi + 140);
    const auto *gsi_141 = buffer.data(gsi + 141);
    const auto *gsi_142 = buffer.data(gsi + 142);
    const auto *gsi_143 = buffer.data(gsi + 143);
    const auto *gsi_144 = buffer.data(gsi + 144);
    const auto *gsi_145 = buffer.data(gsi + 145);
    const auto *gsi_146 = buffer.data(gsi + 146);
    const auto *gsi_147 = buffer.data(gsi + 147);
    const auto *gsi_148 = buffer.data(gsi + 148);
    const auto *gsi_149 = buffer.data(gsi + 149);
    const auto *gsi_150 = buffer.data(gsi + 150);
    const auto *gsi_151 = buffer.data(gsi + 151);
    const auto *gsi_152 = buffer.data(gsi + 152);
    const auto *gsi_153 = buffer.data(gsi + 153);
    const auto *gsi_154 = buffer.data(gsi + 154);
    const auto *gsi_160 = buffer.data(gsi + 160);
    const auto *gsi_161 = buffer.data(gsi + 161);
    const auto *gsi_162 = buffer.data(gsi + 162);
    const auto *gsi_163 = buffer.data(gsi + 163);
    const auto *gsi_164 = buffer.data(gsi + 164);
    const auto *gsi_165 = buffer.data(gsi + 165);
    const auto *gsi_166 = buffer.data(gsi + 166);
    const auto *gsi_167 = buffer.data(gsi + 167);
    const auto *gsi_168 = buffer.data(gsi + 168);
    const auto *gsi_169 = buffer.data(gsi + 169);
    const auto *gsi_170 = buffer.data(gsi + 170);
    const auto *gsi_171 = buffer.data(gsi + 171);
    const auto *gsi_173 = buffer.data(gsi + 173);
    const auto *gsi_174 = buffer.data(gsi + 174);
    const auto *gsi_175 = buffer.data(gsi + 175);
    const auto *gsi_177 = buffer.data(gsi + 177);
    const auto *gsi_178 = buffer.data(gsi + 178);
    const auto *gsi_179 = buffer.data(gsi + 179);
    const auto *gsi_180 = buffer.data(gsi + 180);
    const auto *gsi_182 = buffer.data(gsi + 182);
    const auto *gsi_183 = buffer.data(gsi + 183);
    const auto *gsi_189 = buffer.data(gsi + 189);
    const auto *gsi_191 = buffer.data(gsi + 191);
    const auto *gsi_192 = buffer.data(gsi + 192);
    const auto *gsi_193 = buffer.data(gsi + 193);
    const auto *gsi_194 = buffer.data(gsi + 194);
    const auto *gsi_195 = buffer.data(gsi + 195);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pc_x, pc_z, fsi_105, fsi_107, \
                         fsi_108, fsi_109, gsi_99, gsi_105, gsi_107, gsi_108, \
                         gsi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_14 * fsi_105[k]
                   + f_3 * pc_x[k] * gsi_105[k];

        t_130[k] = f_3 * pc_z[k] * gsi_99[k];

        t_131[k] = f_14 * fsi_107[k]
                   + f_3 * pc_x[k] * gsi_107[k];

        t_132[k] = f_14 * fsi_108[k]
                   + f_3 * pc_x[k] * gsi_108[k];

        t_133[k] = f_14 * fsi_109[k]
                   + f_3 * pc_x[k] * gsi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, fsi_49, fsi_110, \
                         fsi_111, gsh0_78, gsh1_78, gsi_105, gsi_110, \
                         gsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_14 * fsi_110[k]
                   + f_3 * pc_x[k] * gsi_110[k];

        t_135[k] = f_14 * fsi_111[k]
                   + f_3 * pc_x[k] * gsi_111[k];

        t_136[k] = f_14 * fsi_49[k]
                   + f_1 * gsh0_78[k]
                   - f_2 * gsh1_78[k]
                   + f_3 * pc_y[k] * gsi_105[k];

        t_137[k] = f_3 * pc_z[k] * gsi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_z, gsh0_78, gsh0_79, gsh0_80, gsh1_78, \
                         gsh1_79, gsh1_80, gsi_106, gsi_107, gsi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_4 * gsh0_78[k]
                   - f_5 * gsh1_78[k]
                   + f_3 * pc_z[k] * gsi_106[k];

        t_139[k] = f_6 * gsh0_79[k]
                   - f_7 * gsh1_79[k]
                   + f_3 * pc_z[k] * gsi_107[k];

        t_140[k] = f_8 * gsh0_80[k]
                   - f_9 * gsh1_80[k]
                   + f_3 * pc_z[k] * gsi_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pc_y, pc_z, fsk0_72, fsi_55, \
                         fsk1_72, gsh0_81, gsh0_83, gsh1_81, gsh1_83, gsi_109, \
                         gsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_10 * gsh0_81[k]
                   - f_11 * gsh1_81[k]
                   + f_3 * pc_z[k] * gsi_109[k];

        t_142[k] = f_14 * fsi_55[k]
                   + f_3 * pc_y[k] * gsi_111[k];

        t_143[k] = f_1 * gsh0_83[k]
                   - f_2 * gsh1_83[k]
                   + f_3 * pc_z[k] * gsi_111[k];

        t_144[k] = pa_y[k] * fsk0_72[k]
                   - f_12 * pc_y[k] * fsk1_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_z, pc_y, pc_z, fsk0_39, fsi_28, \
                         fsi_56, fsi_58, fsk1_39, gsi_112, gsi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * fsi_56[k]
                   + f_3 * pc_y[k] * gsi_112[k];

        t_146[k] = f_13 * fsi_28[k]
                   + f_3 * pc_z[k] * gsi_112[k];

        t_147[k] = pa_z[k] * fsk0_39[k]
                   - f_12 * pc_z[k] * fsk1_39[k];

        t_148[k] = f_13 * fsi_58[k]
                   + f_3 * pc_y[k] * gsi_114[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pa_z, pc_y, pc_z, fsk0_42, fsk0_77, \
                         fsi_31, fsi_61, fsk1_42, fsk1_77, gsi_115, \
                         gsi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * fsk0_77[k]
                   - f_12 * pc_y[k] * fsk1_77[k];

        t_150[k] = pa_z[k] * fsk0_42[k]
                   - f_12 * pc_z[k] * fsk1_42[k];

        t_151[k] = f_13 * fsi_31[k]
                   + f_3 * pc_z[k] * gsi_115[k];

        t_152[k] = f_13 * fsi_61[k]
                   + f_3 * pc_y[k] * gsi_117[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_y, pa_z, pc_y, pc_z, fsk0_46, fsk0_81, \
                         fsi_34, fsk1_46, fsk1_81, gsi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_y[k] * fsk0_81[k]
                   - f_12 * pc_y[k] * fsk1_81[k];

        t_154[k] = pa_z[k] * fsk0_46[k]
                   - f_12 * pc_z[k] * fsk1_46[k];

        t_155[k] = f_13 * fsi_34[k]
                   + f_3 * pc_z[k] * gsi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_y, pc_y, fsk0_84, fsk0_86, fsi_64, fsi_65, \
                         fsk1_84, fsk1_86, gsi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_y[k] * fsk0_84[k]
                   + f_14 * fsi_64[k]
                   - f_12 * pc_y[k] * fsk1_84[k];

        t_157[k] = f_13 * fsi_65[k]
                   + f_3 * pc_y[k] * gsi_121[k];

        t_158[k] = pa_y[k] * fsk0_86[k]
                   - f_12 * pc_y[k] * fsk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pa_z, pc_y, pc_z, fsk0_51, fsk0_89, \
                         fsi_38, fsi_68, fsk1_51, fsk1_89, gsi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * fsk0_51[k]
                   - f_12 * pc_z[k] * fsk1_51[k];

        t_160[k] = f_13 * fsi_38[k]
                   + f_3 * pc_z[k] * gsi_122[k];

        t_161[k] = pa_y[k] * fsk0_89[k]
                   + f_15 * fsi_68[k]
                   - f_12 * pc_y[k] * fsk1_89[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_y, pc_x, pc_y, fsk0_90, fsk0_92, \
                         fsi_69, fsi_70, fsi_133, fsk1_90, fsk1_92, gsi_126, \
                         gsi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_y[k] * fsk0_90[k]
                   + f_14 * fsi_69[k]
                   - f_12 * pc_y[k] * fsk1_90[k];

        t_163[k] = f_13 * fsi_70[k]
                   + f_3 * pc_y[k] * gsi_126[k];

        t_164[k] = pa_y[k] * fsk0_92[k]
                   - f_12 * pc_y[k] * fsk1_92[k];

        t_165[k] = f_14 * fsi_133[k]
                   + f_3 * pc_x[k] * gsi_133[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, fsi_134, fsi_135, fsi_136, \
                         fsi_137, fsi_138, gsi_134, gsi_135, gsi_136, gsi_137, \
                         gsi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_14 * fsi_134[k]
                   + f_3 * pc_x[k] * gsi_134[k];

        t_167[k] = f_14 * fsi_135[k]
                   + f_3 * pc_x[k] * gsi_135[k];

        t_168[k] = f_14 * fsi_136[k]
                   + f_3 * pc_x[k] * gsi_136[k];

        t_169[k] = f_14 * fsi_137[k]
                   + f_3 * pc_x[k] * gsi_137[k];

        t_170[k] = f_14 * fsi_138[k]
                   + f_3 * pc_x[k] * gsi_138[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_z, pc_x, pc_z, fsk0_64, fsi_49, fsi_139, \
                         fsk1_64, gsi_133, gsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_14 * fsi_139[k]
                   + f_3 * pc_x[k] * gsi_139[k];

        t_172[k] = pa_z[k] * fsk0_64[k]
                   - f_12 * pc_z[k] * fsk1_64[k];

        t_173[k] = f_13 * fsi_49[k]
                   + f_3 * pc_z[k] * gsi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, fsi_79, fsi_80, fsi_81, gsh0_101, \
                         gsh0_102, gsh0_103, gsh1_101, gsh1_102, gsh1_103, gsi_135, gsi_136, \
                         gsi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * fsi_79[k]
                   + f_10 * gsh0_101[k]
                   - f_11 * gsh1_101[k]
                   + f_3 * pc_y[k] * gsi_135[k];

        t_175[k] = f_13 * fsi_80[k]
                   + f_8 * gsh0_102[k]
                   - f_9 * gsh1_102[k]
                   + f_3 * pc_y[k] * gsi_136[k];

        t_176[k] = f_13 * fsi_81[k]
                   + f_6 * gsh0_103[k]
                   - f_7 * gsh1_103[k]
                   + f_3 * pc_y[k] * gsi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, fsk0_107, fsi_82, fsi_83, fsk1_107, \
                         gsh0_104, gsh1_104, gsi_138, gsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * fsi_82[k]
                   + f_4 * gsh0_104[k]
                   - f_5 * gsh1_104[k]
                   + f_3 * pc_y[k] * gsi_138[k];

        t_178[k] = f_13 * fsi_83[k]
                   + f_3 * pc_y[k] * gsi_139[k];

        t_179[k] = pa_y[k] * fsk0_107[k]
                   - f_12 * pc_y[k] * fsk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, fsi_56, fsi_140, \
                         gsh0_105, gsh1_105, gsi_140, gsi_141, \
                         gsi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_14 * fsi_140[k]
                   + f_1 * gsh0_105[k]
                   - f_2 * gsh1_105[k]
                   + f_3 * pc_x[k] * gsi_140[k];

        t_181[k] = f_3 * pc_y[k] * gsi_140[k];

        t_182[k] = f_14 * fsi_56[k]
                   + f_3 * pc_z[k] * gsi_140[k];

        t_183[k] = f_4 * gsh0_105[k]
                   - f_5 * gsh1_105[k]
                   + f_3 * pc_y[k] * gsi_141[k];

        t_184[k] = f_3 * pc_y[k] * gsi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, fsi_145, gsh0_106, gsh0_107, \
                         gsh0_110, gsh1_106, gsh1_107, gsh1_110, gsi_143, gsi_144, \
                         gsi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_14 * fsi_145[k]
                   + f_10 * gsh0_110[k]
                   - f_11 * gsh1_110[k]
                   + f_3 * pc_x[k] * gsi_145[k];

        t_186[k] = f_6 * gsh0_106[k]
                   - f_7 * gsh1_106[k]
                   + f_3 * pc_y[k] * gsi_143[k];

        t_187[k] = f_4 * gsh0_107[k]
                   - f_5 * gsh1_107[k]
                   + f_3 * pc_y[k] * gsi_144[k];

        t_188[k] = f_3 * pc_y[k] * gsi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, fsi_149, gsh0_108, gsh0_109, \
                         gsh0_114, gsh1_108, gsh1_109, gsh1_114, gsi_146, gsi_147, \
                         gsi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * fsi_149[k]
                   + f_8 * gsh0_114[k]
                   - f_9 * gsh1_114[k]
                   + f_3 * pc_x[k] * gsi_149[k];

        t_190[k] = f_8 * gsh0_108[k]
                   - f_9 * gsh1_108[k]
                   + f_3 * pc_y[k] * gsi_146[k];

        t_191[k] = f_6 * gsh0_109[k]
                   - f_7 * gsh1_109[k]
                   + f_3 * pc_y[k] * gsi_147[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pc_x, pc_y, fsi_154, gsh0_110, gsh0_119, \
                         gsh1_110, gsh1_119, gsi_148, gsi_149, \
                         gsi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_4 * gsh0_110[k]
                   - f_5 * gsh1_110[k]
                   + f_3 * pc_y[k] * gsi_148[k];

        t_193[k] = f_3 * pc_y[k] * gsi_149[k];

        t_194[k] = f_14 * fsi_154[k]
                   + f_6 * gsh0_119[k]
                   - f_7 * gsh1_119[k]
                   + f_3 * pc_x[k] * gsi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, gsh0_111, gsh0_112, gsh0_113, gsh1_111, \
                         gsh1_112, gsh1_113, gsi_150, gsi_151, \
                         gsi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * gsh0_111[k]
                   - f_11 * gsh1_111[k]
                   + f_3 * pc_y[k] * gsi_150[k];

        t_196[k] = f_8 * gsh0_112[k]
                   - f_9 * gsh1_112[k]
                   + f_3 * pc_y[k] * gsi_151[k];

        t_197[k] = f_6 * gsh0_113[k]
                   - f_7 * gsh1_113[k]
                   + f_3 * pc_y[k] * gsi_152[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, fsi_160, fsi_161, gsh0_114, \
                         gsh0_125, gsh1_114, gsh1_125, gsi_153, gsi_154, gsi_160, \
                         gsi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * gsh0_114[k]
                   - f_5 * gsh1_114[k]
                   + f_3 * pc_y[k] * gsi_153[k];

        t_199[k] = f_3 * pc_y[k] * gsi_154[k];

        t_200[k] = f_14 * fsi_160[k]
                   + f_4 * gsh0_125[k]
                   - f_5 * gsh1_125[k]
                   + f_3 * pc_x[k] * gsi_160[k];

        t_201[k] = f_14 * fsi_161[k]
                   + f_3 * pc_x[k] * gsi_161[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pc_x, pc_y, fsi_162, fsi_163, \
                         fsi_164, fsi_165, gsi_160, gsi_162, gsi_163, gsi_164, \
                         gsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_14 * fsi_162[k]
                   + f_3 * pc_x[k] * gsi_162[k];

        t_203[k] = f_14 * fsi_163[k]
                   + f_3 * pc_x[k] * gsi_163[k];

        t_204[k] = f_14 * fsi_164[k]
                   + f_3 * pc_x[k] * gsi_164[k];

        t_205[k] = f_14 * fsi_165[k]
                   + f_3 * pc_x[k] * gsi_165[k];

        t_206[k] = f_3 * pc_y[k] * gsi_160[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pc_x, pc_y, fsi_167, gsh0_120, gsh0_121, \
                         gsh1_120, gsh1_121, gsi_161, gsi_162, \
                         gsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_14 * fsi_167[k]
                   + f_3 * pc_x[k] * gsi_167[k];

        t_208[k] = f_1 * gsh0_120[k]
                   - f_2 * gsh1_120[k]
                   + f_3 * pc_y[k] * gsi_161[k];

        t_209[k] = f_17 * gsh0_121[k]
                   - f_18 * gsh1_121[k]
                   + f_3 * pc_y[k] * gsi_162[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pc_y, gsh0_122, gsh0_123, gsh0_124, gsh1_122, \
                         gsh1_123, gsh1_124, gsi_163, gsi_164, \
                         gsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * gsh0_122[k]
                   - f_11 * gsh1_122[k]
                   + f_3 * pc_y[k] * gsi_163[k];

        t_211[k] = f_8 * gsh0_123[k]
                   - f_9 * gsh1_123[k]
                   + f_3 * pc_y[k] * gsi_164[k];

        t_212[k] = f_6 * gsh0_124[k]
                   - f_7 * gsh1_124[k]
                   + f_3 * pc_y[k] * gsi_165[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_x, pc_x, pc_y, pc_z, fsk0_216, fsi_83, \
                         fsi_168, fsk1_216, gsh0_125, gsh1_125, gsi_166, \
                         gsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * gsh0_125[k]
                   - f_5 * gsh1_125[k]
                   + f_3 * pc_y[k] * gsi_166[k];

        t_214[k] = f_3 * pc_y[k] * gsi_167[k];

        t_215[k] = f_14 * fsi_83[k]
                   + f_1 * gsh0_125[k]
                   - f_2 * gsh1_125[k]
                   + f_3 * pc_z[k] * gsi_167[k];

        t_216[k] = pa_x[k] * fsk0_216[k]
                   + f_19 * fsi_168[k]
                   - f_12 * pc_x[k] * fsk1_216[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_x, pc_x, pc_y, pc_z, fsk0_219, fsi_84, \
                         fsi_171, fsk1_219, gsi_168, gsi_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * fsi_84[k]
                   + f_3 * pc_y[k] * gsi_168[k];

        t_218[k] = f_3 * pc_z[k] * gsi_168[k];

        t_219[k] = pa_x[k] * fsk0_219[k]
                   + f_16 * fsi_171[k]
                   - f_12 * pc_x[k] * fsk1_219[k];

        t_220[k] = f_3 * pc_z[k] * gsi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pa_x, pc_x, pc_z, fsk0_222, fsi_174, fsk1_222, \
                         gsh0_126, gsh1_126, gsi_170, gsi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_4 * gsh0_126[k]
                   - f_5 * gsh1_126[k]
                   + f_3 * pc_z[k] * gsi_170[k];

        t_222[k] = pa_x[k] * fsk0_222[k]
                   + f_0 * fsi_174[k]
                   - f_12 * pc_x[k] * fsk1_222[k];

        t_223[k] = f_3 * pc_z[k] * gsi_171[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_x, pc_x, pc_y, pc_z, fsk0_226, fsi_89, \
                         fsi_178, fsk1_226, gsh0_128, gsh1_128, gsi_173, \
                         gsi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_15 * fsi_89[k]
                   + f_3 * pc_y[k] * gsi_173[k];

        t_225[k] = f_6 * gsh0_128[k]
                   - f_7 * gsh1_128[k]
                   + f_3 * pc_z[k] * gsi_173[k];

        t_226[k] = pa_x[k] * fsk0_226[k]
                   + f_15 * fsi_178[k]
                   - f_12 * pc_x[k] * fsk1_226[k];

        t_227[k] = f_3 * pc_z[k] * gsi_174[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_y, pc_z, fsi_93, gsh0_129, gsh0_131, \
                         gsh1_129, gsh1_131, gsi_175, gsi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * gsh0_129[k]
                   - f_5 * gsh1_129[k]
                   + f_3 * pc_z[k] * gsi_175[k];

        t_229[k] = f_15 * fsi_93[k]
                   + f_3 * pc_y[k] * gsi_177[k];

        t_230[k] = f_8 * gsh0_131[k]
                   - f_9 * gsh1_131[k]
                   + f_3 * pc_z[k] * gsi_177[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pa_x, pc_x, pc_z, fsk0_231, fsi_183, fsk1_231, \
                         gsh0_132, gsh1_132, gsi_178, gsi_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = pa_x[k] * fsk0_231[k]
                   + f_14 * fsi_183[k]
                   - f_12 * pc_x[k] * fsk1_231[k];

        t_232[k] = f_3 * pc_z[k] * gsi_178[k];

        t_233[k] = f_4 * gsh0_132[k]
                   - f_5 * gsh1_132[k]
                   + f_3 * pc_z[k] * gsi_179[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, fsi_98, fsi_189, \
                         gsh0_133, gsh0_135, gsh1_133, gsh1_135, gsi_180, gsi_182, \
                         gsi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_6 * gsh0_133[k]
                   - f_7 * gsh1_133[k]
                   + f_3 * pc_z[k] * gsi_180[k];

        t_235[k] = f_15 * fsi_98[k]
                   + f_3 * pc_y[k] * gsi_182[k];

        t_236[k] = f_10 * gsh0_135[k]
                   - f_11 * gsh1_135[k]
                   + f_3 * pc_z[k] * gsi_182[k];

        t_237[k] = f_13 * fsi_189[k]
                   + f_3 * pc_x[k] * gsi_189[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pc_x, pc_z, fsi_191, fsi_192, \
                         fsi_193, fsi_194, gsi_183, gsi_191, gsi_192, gsi_193, \
                         gsi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_3 * pc_z[k] * gsi_183[k];

        t_239[k] = f_13 * fsi_191[k]
                   + f_3 * pc_x[k] * gsi_191[k];

        t_240[k] = f_13 * fsi_192[k]
                   + f_3 * pc_x[k] * gsi_192[k];

        t_241[k] = f_13 * fsi_193[k]
                   + f_3 * pc_x[k] * gsi_193[k];

        t_242[k] = f_13 * fsi_194[k]
                   + f_3 * pc_x[k] * gsi_194[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_x, pc_x, pc_z, fsk0_244, fsk0_246, \
                         fsi_195, fsk1_244, fsk1_246, gsi_189, \
                         gsi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_13 * fsi_195[k]
                   + f_3 * pc_x[k] * gsi_195[k];

        t_244[k] = pa_x[k] * fsk0_244[k]
                   - f_12 * pc_x[k] * fsk1_244[k];

        t_245[k] = f_3 * pc_z[k] * gsi_189[k];

        t_246[k] = pa_x[k] * fsk0_246[k]
                   - f_12 * pc_x[k] * fsk1_246[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pa_x, pc_x, pc_y, fsk0_247, fsk0_248, \
                         fsk0_249, fsi_111, fsk1_247, fsk1_248, fsk1_249, \
                         gsi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = pa_x[k] * fsk0_247[k]
                   - f_12 * pc_x[k] * fsk1_247[k];

        t_248[k] = pa_x[k] * fsk0_248[k]
                   - f_12 * pc_x[k] * fsk1_248[k];

        t_249[k] = pa_x[k] * fsk0_249[k]
                   - f_12 * pc_x[k] * fsk1_249[k];

        t_250[k] = f_15 * fsi_111[k]
                   + f_3 * pc_y[k] * gsi_195[k];
    }
}

static auto
compute_prim_gsk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsk0,
                                                          const size_t fsi, const size_t fsk1,
                                                          const size_t gsh0, const size_t gsh1,
                                                          const size_t gsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsk0_108 = buffer.data(fsk0 + 108);
    const auto *fsk0_111 = buffer.data(fsk0 + 111);
    const auto *fsk0_114 = buffer.data(fsk0 + 114);
    const auto *fsk0_118 = buffer.data(fsk0 + 118);
    const auto *fsk0_123 = buffer.data(fsk0 + 123);
    const auto *fsk0_180 = buffer.data(fsk0 + 180);
    const auto *fsk0_185 = buffer.data(fsk0 + 185);
    const auto *fsk0_189 = buffer.data(fsk0 + 189);
    const auto *fsk0_194 = buffer.data(fsk0 + 194);
    const auto *fsk0_200 = buffer.data(fsk0 + 200);
    const auto *fsk0_251 = buffer.data(fsk0 + 251);
    const auto *fsk0_257 = buffer.data(fsk0 + 257);
    const auto *fsk0_261 = buffer.data(fsk0 + 261);
    const auto *fsk0_264 = buffer.data(fsk0 + 264);
    const auto *fsk0_266 = buffer.data(fsk0 + 266);
    const auto *fsk0_269 = buffer.data(fsk0 + 269);
    const auto *fsk0_270 = buffer.data(fsk0 + 270);
    const auto *fsk0_272 = buffer.data(fsk0 + 272);
    const auto *fsk0_280 = buffer.data(fsk0 + 280);
    const auto *fsk0_282 = buffer.data(fsk0 + 282);
    const auto *fsk0_283 = buffer.data(fsk0 + 283);
    const auto *fsk0_284 = buffer.data(fsk0 + 284);
    const auto *fsk0_285 = buffer.data(fsk0 + 285);
    const auto *fsk0_287 = buffer.data(fsk0 + 287);
    const auto *fsk0_291 = buffer.data(fsk0 + 291);
    const auto *fsk0_294 = buffer.data(fsk0 + 294);
    const auto *fsk0_298 = buffer.data(fsk0 + 298);
    const auto *fsk0_300 = buffer.data(fsk0 + 300);
    const auto *fsk0_303 = buffer.data(fsk0 + 303);
    const auto *fsk0_305 = buffer.data(fsk0 + 305);
    const auto *fsk0_306 = buffer.data(fsk0 + 306);
    const auto *fsk0_316 = buffer.data(fsk0 + 316);
    const auto *fsk0_318 = buffer.data(fsk0 + 318);
    const auto *fsk0_319 = buffer.data(fsk0 + 319);
    const auto *fsk0_320 = buffer.data(fsk0 + 320);
    const auto *fsk0_321 = buffer.data(fsk0 + 321);
    const auto *fsk0_323 = buffer.data(fsk0 + 323);
    const auto *fsk0_324 = buffer.data(fsk0 + 324);
    const auto *fsk0_329 = buffer.data(fsk0 + 329);
    const auto *fsk0_333 = buffer.data(fsk0 + 333);
    const auto *fsk0_338 = buffer.data(fsk0 + 338);
    const auto *fsk0_344 = buffer.data(fsk0 + 344);
    const auto *fsk0_352 = buffer.data(fsk0 + 352);
    const auto *fsk0_353 = buffer.data(fsk0 + 353);
    const auto *fsk0_354 = buffer.data(fsk0 + 354);
    const auto *fsk0_355 = buffer.data(fsk0 + 355);
    const auto *fsk0_356 = buffer.data(fsk0 + 356);
    const auto *fsk0_357 = buffer.data(fsk0 + 357);
    const auto *fsk0_359 = buffer.data(fsk0 + 359);

    const auto *fsi_84 = buffer.data(fsi + 84);
    const auto *fsi_87 = buffer.data(fsi + 87);
    const auto *fsi_90 = buffer.data(fsi + 90);
    const auto *fsi_94 = buffer.data(fsi + 94);
    const auto *fsi_105 = buffer.data(fsi + 105);
    const auto *fsi_112 = buffer.data(fsi + 112);
    const auto *fsi_114 = buffer.data(fsi + 114);
    const auto *fsi_115 = buffer.data(fsi + 115);
    const auto *fsi_117 = buffer.data(fsi + 117);
    const auto *fsi_118 = buffer.data(fsi + 118);
    const auto *fsi_121 = buffer.data(fsi + 121);
    const auto *fsi_122 = buffer.data(fsi + 122);
    const auto *fsi_126 = buffer.data(fsi + 126);
    const auto *fsi_133 = buffer.data(fsi + 133);
    const auto *fsi_139 = buffer.data(fsi + 139);
    const auto *fsi_140 = buffer.data(fsi + 140);
    const auto *fsi_142 = buffer.data(fsi + 142);
    const auto *fsi_145 = buffer.data(fsi + 145);
    const auto *fsi_149 = buffer.data(fsi + 149);
    const auto *fsi_154 = buffer.data(fsi + 154);
    const auto *fsi_167 = buffer.data(fsi + 167);
    const auto *fsi_201 = buffer.data(fsi + 201);
    const auto *fsi_205 = buffer.data(fsi + 205);
    const auto *fsi_208 = buffer.data(fsi + 208);
    const auto *fsi_210 = buffer.data(fsi + 210);
    const auto *fsi_213 = buffer.data(fsi + 213);
    const auto *fsi_214 = buffer.data(fsi + 214);
    const auto *fsi_216 = buffer.data(fsi + 216);
    const auto *fsi_217 = buffer.data(fsi + 217);
    const auto *fsi_218 = buffer.data(fsi + 218);
    const auto *fsi_219 = buffer.data(fsi + 219);
    const auto *fsi_220 = buffer.data(fsi + 220);
    const auto *fsi_221 = buffer.data(fsi + 221);
    const auto *fsi_222 = buffer.data(fsi + 222);
    const auto *fsi_223 = buffer.data(fsi + 223);
    const auto *fsi_227 = buffer.data(fsi + 227);
    const auto *fsi_230 = buffer.data(fsi + 230);
    const auto *fsi_234 = buffer.data(fsi + 234);
    const auto *fsi_236 = buffer.data(fsi + 236);
    const auto *fsi_239 = buffer.data(fsi + 239);
    const auto *fsi_241 = buffer.data(fsi + 241);
    const auto *fsi_242 = buffer.data(fsi + 242);
    const auto *fsi_245 = buffer.data(fsi + 245);
    const auto *fsi_246 = buffer.data(fsi + 246);
    const auto *fsi_247 = buffer.data(fsi + 247);
    const auto *fsi_248 = buffer.data(fsi + 248);
    const auto *fsi_249 = buffer.data(fsi + 249);
    const auto *fsi_250 = buffer.data(fsi + 250);
    const auto *fsi_251 = buffer.data(fsi + 251);
    const auto *fsi_252 = buffer.data(fsi + 252);
    const auto *fsi_257 = buffer.data(fsi + 257);
    const auto *fsi_261 = buffer.data(fsi + 261);
    const auto *fsi_266 = buffer.data(fsi + 266);
    const auto *fsi_272 = buffer.data(fsi + 272);
    const auto *fsi_273 = buffer.data(fsi + 273);
    const auto *fsi_274 = buffer.data(fsi + 274);
    const auto *fsi_275 = buffer.data(fsi + 275);
    const auto *fsi_276 = buffer.data(fsi + 276);
    const auto *fsi_277 = buffer.data(fsi + 277);
    const auto *fsi_279 = buffer.data(fsi + 279);

    const auto *fsk1_108 = buffer.data(fsk1 + 108);
    const auto *fsk1_111 = buffer.data(fsk1 + 111);
    const auto *fsk1_114 = buffer.data(fsk1 + 114);
    const auto *fsk1_118 = buffer.data(fsk1 + 118);
    const auto *fsk1_123 = buffer.data(fsk1 + 123);
    const auto *fsk1_180 = buffer.data(fsk1 + 180);
    const auto *fsk1_185 = buffer.data(fsk1 + 185);
    const auto *fsk1_189 = buffer.data(fsk1 + 189);
    const auto *fsk1_194 = buffer.data(fsk1 + 194);
    const auto *fsk1_200 = buffer.data(fsk1 + 200);
    const auto *fsk1_251 = buffer.data(fsk1 + 251);
    const auto *fsk1_257 = buffer.data(fsk1 + 257);
    const auto *fsk1_261 = buffer.data(fsk1 + 261);
    const auto *fsk1_264 = buffer.data(fsk1 + 264);
    const auto *fsk1_266 = buffer.data(fsk1 + 266);
    const auto *fsk1_269 = buffer.data(fsk1 + 269);
    const auto *fsk1_270 = buffer.data(fsk1 + 270);
    const auto *fsk1_272 = buffer.data(fsk1 + 272);
    const auto *fsk1_280 = buffer.data(fsk1 + 280);
    const auto *fsk1_282 = buffer.data(fsk1 + 282);
    const auto *fsk1_283 = buffer.data(fsk1 + 283);
    const auto *fsk1_284 = buffer.data(fsk1 + 284);
    const auto *fsk1_285 = buffer.data(fsk1 + 285);
    const auto *fsk1_287 = buffer.data(fsk1 + 287);
    const auto *fsk1_291 = buffer.data(fsk1 + 291);
    const auto *fsk1_294 = buffer.data(fsk1 + 294);
    const auto *fsk1_298 = buffer.data(fsk1 + 298);
    const auto *fsk1_300 = buffer.data(fsk1 + 300);
    const auto *fsk1_303 = buffer.data(fsk1 + 303);
    const auto *fsk1_305 = buffer.data(fsk1 + 305);
    const auto *fsk1_306 = buffer.data(fsk1 + 306);
    const auto *fsk1_316 = buffer.data(fsk1 + 316);
    const auto *fsk1_318 = buffer.data(fsk1 + 318);
    const auto *fsk1_319 = buffer.data(fsk1 + 319);
    const auto *fsk1_320 = buffer.data(fsk1 + 320);
    const auto *fsk1_321 = buffer.data(fsk1 + 321);
    const auto *fsk1_323 = buffer.data(fsk1 + 323);
    const auto *fsk1_324 = buffer.data(fsk1 + 324);
    const auto *fsk1_329 = buffer.data(fsk1 + 329);
    const auto *fsk1_333 = buffer.data(fsk1 + 333);
    const auto *fsk1_338 = buffer.data(fsk1 + 338);
    const auto *fsk1_344 = buffer.data(fsk1 + 344);
    const auto *fsk1_352 = buffer.data(fsk1 + 352);
    const auto *fsk1_353 = buffer.data(fsk1 + 353);
    const auto *fsk1_354 = buffer.data(fsk1 + 354);
    const auto *fsk1_355 = buffer.data(fsk1 + 355);
    const auto *fsk1_356 = buffer.data(fsk1 + 356);
    const auto *fsk1_357 = buffer.data(fsk1 + 357);
    const auto *fsk1_359 = buffer.data(fsk1 + 359);

    const auto *gsh0_189 = buffer.data(gsh0 + 189);
    const auto *gsh0_190 = buffer.data(gsh0 + 190);
    const auto *gsh0_191 = buffer.data(gsh0 + 191);
    const auto *gsh0_192 = buffer.data(gsh0 + 192);
    const auto *gsh0_193 = buffer.data(gsh0 + 193);
    const auto *gsh0_194 = buffer.data(gsh0 + 194);
    const auto *gsh0_195 = buffer.data(gsh0 + 195);
    const auto *gsh0_196 = buffer.data(gsh0 + 196);
    const auto *gsh0_197 = buffer.data(gsh0 + 197);
    const auto *gsh0_198 = buffer.data(gsh0 + 198);
    const auto *gsh0_210 = buffer.data(gsh0 + 210);
    const auto *gsh0_211 = buffer.data(gsh0 + 211);
    const auto *gsh0_213 = buffer.data(gsh0 + 213);
    const auto *gsh0_215 = buffer.data(gsh0 + 215);
    const auto *gsh0_216 = buffer.data(gsh0 + 216);
    const auto *gsh0_218 = buffer.data(gsh0 + 218);
    const auto *gsh0_219 = buffer.data(gsh0 + 219);
    const auto *gsh0_220 = buffer.data(gsh0 + 220);
    const auto *gsh0_222 = buffer.data(gsh0 + 222);

    const auto *gsh1_189 = buffer.data(gsh1 + 189);
    const auto *gsh1_190 = buffer.data(gsh1 + 190);
    const auto *gsh1_191 = buffer.data(gsh1 + 191);
    const auto *gsh1_192 = buffer.data(gsh1 + 192);
    const auto *gsh1_193 = buffer.data(gsh1 + 193);
    const auto *gsh1_194 = buffer.data(gsh1 + 194);
    const auto *gsh1_195 = buffer.data(gsh1 + 195);
    const auto *gsh1_196 = buffer.data(gsh1 + 196);
    const auto *gsh1_197 = buffer.data(gsh1 + 197);
    const auto *gsh1_198 = buffer.data(gsh1 + 198);
    const auto *gsh1_210 = buffer.data(gsh1 + 210);
    const auto *gsh1_211 = buffer.data(gsh1 + 211);
    const auto *gsh1_213 = buffer.data(gsh1 + 213);
    const auto *gsh1_215 = buffer.data(gsh1 + 215);
    const auto *gsh1_216 = buffer.data(gsh1 + 216);
    const auto *gsh1_218 = buffer.data(gsh1 + 218);
    const auto *gsh1_219 = buffer.data(gsh1 + 219);
    const auto *gsh1_220 = buffer.data(gsh1 + 220);
    const auto *gsh1_222 = buffer.data(gsh1 + 222);

    const auto *gsi_196 = buffer.data(gsi + 196);
    const auto *gsi_198 = buffer.data(gsi + 198);
    const auto *gsi_199 = buffer.data(gsi + 199);
    const auto *gsi_201 = buffer.data(gsi + 201);
    const auto *gsi_202 = buffer.data(gsi + 202);
    const auto *gsi_205 = buffer.data(gsi + 205);
    const auto *gsi_206 = buffer.data(gsi + 206);
    const auto *gsi_210 = buffer.data(gsi + 210);
    const auto *gsi_217 = buffer.data(gsi + 217);
    const auto *gsi_218 = buffer.data(gsi + 218);
    const auto *gsi_219 = buffer.data(gsi + 219);
    const auto *gsi_220 = buffer.data(gsi + 220);
    const auto *gsi_221 = buffer.data(gsi + 221);
    const auto *gsi_222 = buffer.data(gsi + 222);
    const auto *gsi_223 = buffer.data(gsi + 223);
    const auto *gsi_224 = buffer.data(gsi + 224);
    const auto *gsi_226 = buffer.data(gsi + 226);
    const auto *gsi_227 = buffer.data(gsi + 227);
    const auto *gsi_229 = buffer.data(gsi + 229);
    const auto *gsi_230 = buffer.data(gsi + 230);
    const auto *gsi_233 = buffer.data(gsi + 233);
    const auto *gsi_234 = buffer.data(gsi + 234);
    const auto *gsi_238 = buffer.data(gsi + 238);
    const auto *gsi_245 = buffer.data(gsi + 245);
    const auto *gsi_246 = buffer.data(gsi + 246);
    const auto *gsi_247 = buffer.data(gsi + 247);
    const auto *gsi_248 = buffer.data(gsi + 248);
    const auto *gsi_249 = buffer.data(gsi + 249);
    const auto *gsi_250 = buffer.data(gsi + 250);
    const auto *gsi_251 = buffer.data(gsi + 251);
    const auto *gsi_252 = buffer.data(gsi + 252);
    const auto *gsi_253 = buffer.data(gsi + 253);
    const auto *gsi_254 = buffer.data(gsi + 254);
    const auto *gsi_255 = buffer.data(gsi + 255);
    const auto *gsi_256 = buffer.data(gsi + 256);
    const auto *gsi_257 = buffer.data(gsi + 257);
    const auto *gsi_258 = buffer.data(gsi + 258);
    const auto *gsi_259 = buffer.data(gsi + 259);
    const auto *gsi_260 = buffer.data(gsi + 260);
    const auto *gsi_261 = buffer.data(gsi + 261);
    const auto *gsi_262 = buffer.data(gsi + 262);
    const auto *gsi_263 = buffer.data(gsi + 263);
    const auto *gsi_264 = buffer.data(gsi + 264);
    const auto *gsi_265 = buffer.data(gsi + 265);
    const auto *gsi_266 = buffer.data(gsi + 266);
    const auto *gsi_272 = buffer.data(gsi + 272);
    const auto *gsi_273 = buffer.data(gsi + 273);
    const auto *gsi_274 = buffer.data(gsi + 274);
    const auto *gsi_275 = buffer.data(gsi + 275);
    const auto *gsi_276 = buffer.data(gsi + 276);
    const auto *gsi_277 = buffer.data(gsi + 277);
    const auto *gsi_279 = buffer.data(gsi + 279);
    const auto *gsi_280 = buffer.data(gsi + 280);
    const auto *gsi_281 = buffer.data(gsi + 281);
    const auto *gsi_283 = buffer.data(gsi + 283);
    const auto *gsi_285 = buffer.data(gsi + 285);
    const auto *gsi_286 = buffer.data(gsi + 286);
    const auto *gsi_288 = buffer.data(gsi + 288);
    const auto *gsi_289 = buffer.data(gsi + 289);
    const auto *gsi_290 = buffer.data(gsi + 290);
    const auto *gsi_292 = buffer.data(gsi + 292);

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_x, pa_z, pc_x, pc_y, pc_z, fsk0_108, \
                         fsk0_251, fsi_84, fsi_112, fsk1_108, fsk1_251, \
                         gsi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_x[k] * fsk0_251[k]
                   - f_12 * pc_x[k] * fsk1_251[k];

        t_252[k] = pa_z[k] * fsk0_108[k]
                   - f_12 * pc_z[k] * fsk1_108[k];

        t_253[k] = f_14 * fsi_112[k]
                   + f_3 * pc_y[k] * gsi_196[k];

        t_254[k] = f_13 * fsi_84[k]
                   + f_3 * pc_z[k] * gsi_196[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pa_x, pa_z, pc_x, pc_y, pc_z, fsk0_111, \
                         fsk0_257, fsi_114, fsi_201, fsk1_111, fsk1_257, \
                         gsi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = pa_z[k] * fsk0_111[k]
                   - f_12 * pc_z[k] * fsk1_111[k];

        t_256[k] = f_14 * fsi_114[k]
                   + f_3 * pc_y[k] * gsi_198[k];

        t_257[k] = pa_x[k] * fsk0_257[k]
                   + f_16 * fsi_201[k]
                   - f_12 * pc_x[k] * fsk1_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pa_z, pc_y, pc_z, fsk0_114, fsi_87, fsi_117, \
                         fsk1_114, gsi_199, gsi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pa_z[k] * fsk0_114[k]
                   - f_12 * pc_z[k] * fsk1_114[k];

        t_259[k] = f_13 * fsi_87[k]
                   + f_3 * pc_z[k] * gsi_199[k];

        t_260[k] = f_14 * fsi_117[k]
                   + f_3 * pc_y[k] * gsi_201[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pa_x, pa_z, pc_x, pc_z, fsk0_118, fsk0_261, \
                         fsi_90, fsi_205, fsk1_118, fsk1_261, gsi_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pa_x[k] * fsk0_261[k]
                   + f_0 * fsi_205[k]
                   - f_12 * pc_x[k] * fsk1_261[k];

        t_262[k] = pa_z[k] * fsk0_118[k]
                   - f_12 * pc_z[k] * fsk1_118[k];

        t_263[k] = f_13 * fsi_90[k]
                   + f_3 * pc_z[k] * gsi_202[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_x, pc_x, pc_y, fsk0_264, fsk0_266, fsi_121, \
                         fsi_208, fsi_210, fsk1_264, fsk1_266, \
                         gsi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = pa_x[k] * fsk0_264[k]
                   + f_15 * fsi_208[k]
                   - f_12 * pc_x[k] * fsk1_264[k];

        t_265[k] = f_14 * fsi_121[k]
                   + f_3 * pc_y[k] * gsi_205[k];

        t_266[k] = pa_x[k] * fsk0_266[k]
                   + f_15 * fsi_210[k]
                   - f_12 * pc_x[k] * fsk1_266[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_x, pa_z, pc_x, pc_z, fsk0_123, fsk0_269, \
                         fsi_94, fsi_213, fsk1_123, fsk1_269, gsi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_z[k] * fsk0_123[k]
                   - f_12 * pc_z[k] * fsk1_123[k];

        t_268[k] = f_13 * fsi_94[k]
                   + f_3 * pc_z[k] * gsi_206[k];

        t_269[k] = pa_x[k] * fsk0_269[k]
                   + f_14 * fsi_213[k]
                   - f_12 * pc_x[k] * fsk1_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pa_x, pc_x, pc_y, fsk0_270, fsk0_272, fsi_126, \
                         fsi_214, fsi_216, fsk1_270, fsk1_272, \
                         gsi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = pa_x[k] * fsk0_270[k]
                   + f_14 * fsi_214[k]
                   - f_12 * pc_x[k] * fsk1_270[k];

        t_271[k] = f_14 * fsi_126[k]
                   + f_3 * pc_y[k] * gsi_210[k];

        t_272[k] = pa_x[k] * fsk0_272[k]
                   + f_14 * fsi_216[k]
                   - f_12 * pc_x[k] * fsk1_272[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, pc_x, fsi_217, fsi_218, fsi_219, \
                         fsi_220, fsi_221, gsi_217, gsi_218, gsi_219, gsi_220, \
                         gsi_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_13 * fsi_217[k]
                   + f_3 * pc_x[k] * gsi_217[k];

        t_274[k] = f_13 * fsi_218[k]
                   + f_3 * pc_x[k] * gsi_218[k];

        t_275[k] = f_13 * fsi_219[k]
                   + f_3 * pc_x[k] * gsi_219[k];

        t_276[k] = f_13 * fsi_220[k]
                   + f_3 * pc_x[k] * gsi_220[k];

        t_277[k] = f_13 * fsi_221[k]
                   + f_3 * pc_x[k] * gsi_221[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pa_x, pc_x, pc_z, fsk0_280, fsi_105, \
                         fsi_222, fsi_223, fsk1_280, gsi_217, gsi_222, \
                         gsi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_13 * fsi_222[k]
                   + f_3 * pc_x[k] * gsi_222[k];

        t_279[k] = f_13 * fsi_223[k]
                   + f_3 * pc_x[k] * gsi_223[k];

        t_280[k] = pa_x[k] * fsk0_280[k]
                   - f_12 * pc_x[k] * fsk1_280[k];

        t_281[k] = f_13 * fsi_105[k]
                   + f_3 * pc_z[k] * gsi_217[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_x, pc_x, fsk0_282, fsk0_283, fsk0_284, \
                         fsk0_285, fsk1_282, fsk1_283, fsk1_284, \
                         fsk1_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = pa_x[k] * fsk0_282[k]
                   - f_12 * pc_x[k] * fsk1_282[k];

        t_283[k] = pa_x[k] * fsk0_283[k]
                   - f_12 * pc_x[k] * fsk1_283[k];

        t_284[k] = pa_x[k] * fsk0_284[k]
                   - f_12 * pc_x[k] * fsk1_284[k];

        t_285[k] = pa_x[k] * fsk0_285[k]
                   - f_12 * pc_x[k] * fsk1_285[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_x, pa_y, pc_x, pc_y, fsk0_180, \
                         fsk0_287, fsi_139, fsi_140, fsk1_180, fsk1_287, gsi_223, \
                         gsi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * fsi_139[k]
                   + f_3 * pc_y[k] * gsi_223[k];

        t_287[k] = pa_x[k] * fsk0_287[k]
                   - f_12 * pc_x[k] * fsk1_287[k];

        t_288[k] = pa_y[k] * fsk0_180[k]
                   - f_12 * pc_y[k] * fsk1_180[k];

        t_289[k] = f_13 * fsi_140[k]
                   + f_3 * pc_y[k] * gsi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pa_x, pc_x, pc_y, pc_z, fsk0_291, fsi_112, \
                         fsi_142, fsi_227, fsk1_291, gsi_224, gsi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * fsi_112[k]
                   + f_3 * pc_z[k] * gsi_224[k];

        t_291[k] = pa_x[k] * fsk0_291[k]
                   + f_16 * fsi_227[k]
                   - f_12 * pc_x[k] * fsk1_291[k];

        t_292[k] = f_13 * fsi_142[k]
                   + f_3 * pc_y[k] * gsi_226[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pa_x, pa_y, pc_x, pc_y, pc_z, fsk0_185, \
                         fsk0_294, fsi_115, fsi_230, fsk1_185, fsk1_294, \
                         gsi_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_y[k] * fsk0_185[k]
                   - f_12 * pc_y[k] * fsk1_185[k];

        t_294[k] = pa_x[k] * fsk0_294[k]
                   + f_0 * fsi_230[k]
                   - f_12 * pc_x[k] * fsk1_294[k];

        t_295[k] = f_14 * fsi_115[k]
                   + f_3 * pc_z[k] * gsi_227[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, pa_x, pa_y, pc_x, pc_y, fsk0_189, fsk0_298, \
                         fsi_145, fsi_234, fsk1_189, fsk1_298, \
                         gsi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_13 * fsi_145[k]
                   + f_3 * pc_y[k] * gsi_229[k];

        t_297[k] = pa_y[k] * fsk0_189[k]
                   - f_12 * pc_y[k] * fsk1_189[k];

        t_298[k] = pa_x[k] * fsk0_298[k]
                   + f_15 * fsi_234[k]
                   - f_12 * pc_x[k] * fsk1_298[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pa_x, pc_x, pc_y, pc_z, fsk0_300, fsi_118, \
                         fsi_149, fsi_236, fsk1_300, gsi_230, gsi_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_14 * fsi_118[k]
                   + f_3 * pc_z[k] * gsi_230[k];

        t_300[k] = pa_x[k] * fsk0_300[k]
                   + f_15 * fsi_236[k]
                   - f_12 * pc_x[k] * fsk1_300[k];

        t_301[k] = f_13 * fsi_149[k]
                   + f_3 * pc_y[k] * gsi_233[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pa_x, pa_y, pc_x, pc_y, pc_z, fsk0_194, \
                         fsk0_303, fsi_122, fsi_239, fsk1_194, fsk1_303, \
                         gsi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = pa_y[k] * fsk0_194[k]
                   - f_12 * pc_y[k] * fsk1_194[k];

        t_303[k] = pa_x[k] * fsk0_303[k]
                   + f_14 * fsi_239[k]
                   - f_12 * pc_x[k] * fsk1_303[k];

        t_304[k] = f_14 * fsi_122[k]
                   + f_3 * pc_z[k] * gsi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pa_x, pc_x, pc_y, fsk0_305, fsk0_306, fsi_154, \
                         fsi_241, fsi_242, fsk1_305, fsk1_306, \
                         gsi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_x[k] * fsk0_305[k]
                   + f_14 * fsi_241[k]
                   - f_12 * pc_x[k] * fsk1_305[k];

        t_306[k] = pa_x[k] * fsk0_306[k]
                   + f_14 * fsi_242[k]
                   - f_12 * pc_x[k] * fsk1_306[k];

        t_307[k] = f_13 * fsi_154[k]
                   + f_3 * pc_y[k] * gsi_238[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_y, pc_x, pc_y, fsk0_200, fsi_245, \
                         fsi_246, fsi_247, fsk1_200, gsi_245, gsi_246, \
                         gsi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pa_y[k] * fsk0_200[k]
                   - f_12 * pc_y[k] * fsk1_200[k];

        t_309[k] = f_13 * fsi_245[k]
                   + f_3 * pc_x[k] * gsi_245[k];

        t_310[k] = f_13 * fsi_246[k]
                   + f_3 * pc_x[k] * gsi_246[k];

        t_311[k] = f_13 * fsi_247[k]
                   + f_3 * pc_x[k] * gsi_247[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, fsi_248, fsi_249, fsi_250, fsi_251, \
                         gsi_248, gsi_249, gsi_250, gsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * fsi_248[k]
                   + f_3 * pc_x[k] * gsi_248[k];

        t_313[k] = f_13 * fsi_249[k]
                   + f_3 * pc_x[k] * gsi_249[k];

        t_314[k] = f_13 * fsi_250[k]
                   + f_3 * pc_x[k] * gsi_250[k];

        t_315[k] = f_13 * fsi_251[k]
                   + f_3 * pc_x[k] * gsi_251[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_x, pc_x, pc_z, fsk0_316, fsk0_318, \
                         fsk0_319, fsi_133, fsk1_316, fsk1_318, fsk1_319, \
                         gsi_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_x[k] * fsk0_316[k]
                   - f_12 * pc_x[k] * fsk1_316[k];

        t_317[k] = f_14 * fsi_133[k]
                   + f_3 * pc_z[k] * gsi_245[k];

        t_318[k] = pa_x[k] * fsk0_318[k]
                   - f_12 * pc_x[k] * fsk1_318[k];

        t_319[k] = pa_x[k] * fsk0_319[k]
                   - f_12 * pc_x[k] * fsk1_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_x, pc_x, pc_y, fsk0_320, fsk0_321, \
                         fsk0_323, fsi_167, fsk1_320, fsk1_321, fsk1_323, \
                         gsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pa_x[k] * fsk0_320[k]
                   - f_12 * pc_x[k] * fsk1_320[k];

        t_321[k] = pa_x[k] * fsk0_321[k]
                   - f_12 * pc_x[k] * fsk1_321[k];

        t_322[k] = f_13 * fsi_167[k]
                   + f_3 * pc_y[k] * gsi_251[k];

        t_323[k] = pa_x[k] * fsk0_323[k]
                   - f_12 * pc_x[k] * fsk1_323[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_x, pc_x, pc_y, pc_z, fsk0_324, \
                         fsi_140, fsi_252, fsk1_324, gsh0_189, gsh1_189, gsi_252, \
                         gsi_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = pa_x[k] * fsk0_324[k]
                   + f_19 * fsi_252[k]
                   - f_12 * pc_x[k] * fsk1_324[k];

        t_325[k] = f_3 * pc_y[k] * gsi_252[k];

        t_326[k] = f_15 * fsi_140[k]
                   + f_3 * pc_z[k] * gsi_252[k];

        t_327[k] = f_4 * gsh0_189[k]
                   - f_5 * gsh1_189[k]
                   + f_3 * pc_y[k] * gsi_253[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pa_x, pc_x, pc_y, fsk0_329, fsi_257, fsk1_329, \
                         gsh0_190, gsh1_190, gsi_254, gsi_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_3 * pc_y[k] * gsi_254[k];

        t_329[k] = pa_x[k] * fsk0_329[k]
                   + f_16 * fsi_257[k]
                   - f_12 * pc_x[k] * fsk1_329[k];

        t_330[k] = f_6 * gsh0_190[k]
                   - f_7 * gsh1_190[k]
                   + f_3 * pc_y[k] * gsi_255[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pa_x, pc_x, pc_y, fsk0_333, fsi_261, fsk1_333, \
                         gsh0_191, gsh1_191, gsi_256, gsi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_4 * gsh0_191[k]
                   - f_5 * gsh1_191[k]
                   + f_3 * pc_y[k] * gsi_256[k];

        t_332[k] = f_3 * pc_y[k] * gsi_257[k];

        t_333[k] = pa_x[k] * fsk0_333[k]
                   + f_0 * fsi_261[k]
                   - f_12 * pc_x[k] * fsk1_333[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pc_y, gsh0_192, gsh0_193, gsh0_194, \
                         gsh1_192, gsh1_193, gsh1_194, gsi_258, gsi_259, gsi_260, \
                         gsi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_8 * gsh0_192[k]
                   - f_9 * gsh1_192[k]
                   + f_3 * pc_y[k] * gsi_258[k];

        t_335[k] = f_6 * gsh0_193[k]
                   - f_7 * gsh1_193[k]
                   + f_3 * pc_y[k] * gsi_259[k];

        t_336[k] = f_4 * gsh0_194[k]
                   - f_5 * gsh1_194[k]
                   + f_3 * pc_y[k] * gsi_260[k];

        t_337[k] = f_3 * pc_y[k] * gsi_261[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pa_x, pc_x, pc_y, fsk0_338, fsi_266, fsk1_338, \
                         gsh0_195, gsh0_196, gsh1_195, gsh1_196, gsi_262, \
                         gsi_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pa_x[k] * fsk0_338[k]
                   + f_15 * fsi_266[k]
                   - f_12 * pc_x[k] * fsk1_338[k];

        t_339[k] = f_10 * gsh0_195[k]
                   - f_11 * gsh1_195[k]
                   + f_3 * pc_y[k] * gsi_262[k];

        t_340[k] = f_8 * gsh0_196[k]
                   - f_9 * gsh1_196[k]
                   + f_3 * pc_y[k] * gsi_263[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pc_y, gsh0_197, gsh0_198, gsh1_197, gsh1_198, \
                         gsi_264, gsi_265, gsi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_6 * gsh0_197[k]
                   - f_7 * gsh1_197[k]
                   + f_3 * pc_y[k] * gsi_264[k];

        t_342[k] = f_4 * gsh0_198[k]
                   - f_5 * gsh1_198[k]
                   + f_3 * pc_y[k] * gsi_265[k];

        t_343[k] = f_3 * pc_y[k] * gsi_266[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_x, pc_x, fsk0_344, fsi_272, fsi_273, \
                         fsi_274, fsi_275, fsk1_344, gsi_273, gsi_274, \
                         gsi_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = pa_x[k] * fsk0_344[k]
                   + f_14 * fsi_272[k]
                   - f_12 * pc_x[k] * fsk1_344[k];

        t_345[k] = f_13 * fsi_273[k]
                   + f_3 * pc_x[k] * gsi_273[k];

        t_346[k] = f_13 * fsi_274[k]
                   + f_3 * pc_x[k] * gsi_274[k];

        t_347[k] = f_13 * fsi_275[k]
                   + f_3 * pc_x[k] * gsi_275[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pc_x, pc_y, fsi_276, fsi_277, fsi_279, \
                         gsi_272, gsi_276, gsi_277, gsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_13 * fsi_276[k]
                   + f_3 * pc_x[k] * gsi_276[k];

        t_349[k] = f_13 * fsi_277[k]
                   + f_3 * pc_x[k] * gsi_277[k];

        t_350[k] = f_3 * pc_y[k] * gsi_272[k];

        t_351[k] = f_13 * fsi_279[k]
                   + f_3 * pc_x[k] * gsi_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_x, pc_x, fsk0_352, fsk0_353, fsk0_354, \
                         fsk0_355, fsk1_352, fsk1_353, fsk1_354, \
                         fsk1_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = pa_x[k] * fsk0_352[k]
                   - f_12 * pc_x[k] * fsk1_352[k];

        t_353[k] = pa_x[k] * fsk0_353[k]
                   - f_12 * pc_x[k] * fsk1_353[k];

        t_354[k] = pa_x[k] * fsk0_354[k]
                   - f_12 * pc_x[k] * fsk1_354[k];

        t_355[k] = pa_x[k] * fsk0_355[k]
                   - f_12 * pc_x[k] * fsk1_355[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_x, pc_x, pc_y, fsk0_356, fsk0_357, \
                         fsk0_359, fsk1_356, fsk1_357, fsk1_359, \
                         gsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pa_x[k] * fsk0_356[k]
                   - f_12 * pc_x[k] * fsk1_356[k];

        t_357[k] = pa_x[k] * fsk0_357[k]
                   - f_12 * pc_x[k] * fsk1_357[k];

        t_358[k] = f_3 * pc_y[k] * gsi_279[k];

        t_359[k] = pa_x[k] * fsk0_359[k]
                   - f_12 * pc_x[k] * fsk1_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, pc_x, pc_z, gsh0_210, gsh0_211, \
                         gsh0_213, gsh1_210, gsh1_211, gsh1_213, gsi_280, gsi_281, \
                         gsi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_1 * gsh0_210[k]
                   - f_2 * gsh1_210[k]
                   + f_3 * pc_x[k] * gsi_280[k];

        t_361[k] = f_17 * gsh0_211[k]
                   - f_18 * gsh1_211[k]
                   + f_3 * pc_x[k] * gsi_281[k];

        t_362[k] = f_3 * pc_z[k] * gsi_280[k];

        t_363[k] = f_10 * gsh0_213[k]
                   - f_11 * gsh1_213[k]
                   + f_3 * pc_x[k] * gsi_283[k];

        t_364[k] = f_3 * pc_z[k] * gsi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_z, gsh0_215, gsh0_216, gsh0_218, \
                         gsh1_215, gsh1_216, gsh1_218, gsi_283, gsi_285, gsi_286, \
                         gsi_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_10 * gsh0_215[k]
                   - f_11 * gsh1_215[k]
                   + f_3 * pc_x[k] * gsi_285[k];

        t_366[k] = f_8 * gsh0_216[k]
                   - f_9 * gsh1_216[k]
                   + f_3 * pc_x[k] * gsi_286[k];

        t_367[k] = f_3 * pc_z[k] * gsi_283[k];

        t_368[k] = f_8 * gsh0_218[k]
                   - f_9 * gsh1_218[k]
                   + f_3 * pc_x[k] * gsi_288[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_z, gsh0_219, gsh0_220, gsh0_222, \
                         gsh1_219, gsh1_220, gsh1_222, gsi_286, gsi_289, gsi_290, \
                         gsi_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_8 * gsh0_219[k]
                   - f_9 * gsh1_219[k]
                   + f_3 * pc_x[k] * gsi_289[k];

        t_370[k] = f_6 * gsh0_220[k]
                   - f_7 * gsh1_220[k]
                   + f_3 * pc_x[k] * gsi_290[k];

        t_371[k] = f_3 * pc_z[k] * gsi_286[k];

        t_372[k] = f_6 * gsh0_222[k]
                   - f_7 * gsh1_222[k]
                   + f_3 * pc_x[k] * gsi_292[k];
    }
}

static auto
compute_prim_gsk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsk0,
                                                          const size_t fsi, const size_t fsk1,
                                                          const size_t gsh0, const size_t gsh1,
                                                          const size_t gsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);

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
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsk0_216 = buffer.data(fsk0 + 216);
    const auto *fsk0_217 = buffer.data(fsk0 + 217);
    const auto *fsk0_219 = buffer.data(fsk0 + 219);
    const auto *fsk0_222 = buffer.data(fsk0 + 222);
    const auto *fsk0_226 = buffer.data(fsk0 + 226);
    const auto *fsk0_231 = buffer.data(fsk0 + 231);
    const auto *fsk0_244 = buffer.data(fsk0 + 244);
    const auto *fsk0_246 = buffer.data(fsk0 + 246);
    const auto *fsk0_247 = buffer.data(fsk0 + 247);
    const auto *fsk0_248 = buffer.data(fsk0 + 248);
    const auto *fsk0_249 = buffer.data(fsk0 + 249);
    const auto *fsk0_324 = buffer.data(fsk0 + 324);
    const auto *fsk0_326 = buffer.data(fsk0 + 326);
    const auto *fsk0_329 = buffer.data(fsk0 + 329);
    const auto *fsk0_333 = buffer.data(fsk0 + 333);
    const auto *fsk0_338 = buffer.data(fsk0 + 338);
    const auto *fsk0_344 = buffer.data(fsk0 + 344);

    const auto *fsi_189 = buffer.data(fsi + 189);
    const auto *fsi_190 = buffer.data(fsi + 190);
    const auto *fsi_191 = buffer.data(fsi + 191);
    const auto *fsi_192 = buffer.data(fsi + 192);
    const auto *fsi_193 = buffer.data(fsi + 193);
    const auto *fsi_195 = buffer.data(fsi + 195);
    const auto *fsi_217 = buffer.data(fsi + 217);
    const auto *fsi_223 = buffer.data(fsi + 223);
    const auto *fsi_245 = buffer.data(fsi + 245);
    const auto *fsi_247 = buffer.data(fsi + 247);
    const auto *fsi_248 = buffer.data(fsi + 248);
    const auto *fsi_249 = buffer.data(fsi + 249);
    const auto *fsi_250 = buffer.data(fsi + 250);
    const auto *fsi_251 = buffer.data(fsi + 251);

    const auto *fsk1_216 = buffer.data(fsk1 + 216);
    const auto *fsk1_217 = buffer.data(fsk1 + 217);
    const auto *fsk1_219 = buffer.data(fsk1 + 219);
    const auto *fsk1_222 = buffer.data(fsk1 + 222);
    const auto *fsk1_226 = buffer.data(fsk1 + 226);
    const auto *fsk1_231 = buffer.data(fsk1 + 231);
    const auto *fsk1_244 = buffer.data(fsk1 + 244);
    const auto *fsk1_246 = buffer.data(fsk1 + 246);
    const auto *fsk1_247 = buffer.data(fsk1 + 247);
    const auto *fsk1_248 = buffer.data(fsk1 + 248);
    const auto *fsk1_249 = buffer.data(fsk1 + 249);
    const auto *fsk1_324 = buffer.data(fsk1 + 324);
    const auto *fsk1_326 = buffer.data(fsk1 + 326);
    const auto *fsk1_329 = buffer.data(fsk1 + 329);
    const auto *fsk1_333 = buffer.data(fsk1 + 333);
    const auto *fsk1_338 = buffer.data(fsk1 + 338);
    const auto *fsk1_344 = buffer.data(fsk1 + 344);

    const auto *gsh0_223 = buffer.data(gsh0 + 223);
    const auto *gsh0_224 = buffer.data(gsh0 + 224);
    const auto *gsh0_225 = buffer.data(gsh0 + 225);
    const auto *gsh0_226 = buffer.data(gsh0 + 226);
    const auto *gsh0_227 = buffer.data(gsh0 + 227);
    const auto *gsh0_228 = buffer.data(gsh0 + 228);
    const auto *gsh0_229 = buffer.data(gsh0 + 229);
    const auto *gsh0_230 = buffer.data(gsh0 + 230);
    const auto *gsh0_233 = buffer.data(gsh0 + 233);
    const auto *gsh0_235 = buffer.data(gsh0 + 235);
    const auto *gsh0_236 = buffer.data(gsh0 + 236);
    const auto *gsh0_238 = buffer.data(gsh0 + 238);
    const auto *gsh0_239 = buffer.data(gsh0 + 239);
    const auto *gsh0_240 = buffer.data(gsh0 + 240);
    const auto *gsh0_242 = buffer.data(gsh0 + 242);
    const auto *gsh0_243 = buffer.data(gsh0 + 243);
    const auto *gsh0_244 = buffer.data(gsh0 + 244);
    const auto *gsh0_245 = buffer.data(gsh0 + 245);
    const auto *gsh0_247 = buffer.data(gsh0 + 247);
    const auto *gsh0_248 = buffer.data(gsh0 + 248);
    const auto *gsh0_249 = buffer.data(gsh0 + 249);
    const auto *gsh0_250 = buffer.data(gsh0 + 250);
    const auto *gsh0_251 = buffer.data(gsh0 + 251);
    const auto *gsh0_252 = buffer.data(gsh0 + 252);
    const auto *gsh0_253 = buffer.data(gsh0 + 253);
    const auto *gsh0_254 = buffer.data(gsh0 + 254);
    const auto *gsh0_255 = buffer.data(gsh0 + 255);
    const auto *gsh0_256 = buffer.data(gsh0 + 256);
    const auto *gsh0_257 = buffer.data(gsh0 + 257);
    const auto *gsh0_258 = buffer.data(gsh0 + 258);
    const auto *gsh0_259 = buffer.data(gsh0 + 259);
    const auto *gsh0_260 = buffer.data(gsh0 + 260);
    const auto *gsh0_261 = buffer.data(gsh0 + 261);
    const auto *gsh0_262 = buffer.data(gsh0 + 262);
    const auto *gsh0_263 = buffer.data(gsh0 + 263);
    const auto *gsh0_264 = buffer.data(gsh0 + 264);
    const auto *gsh0_265 = buffer.data(gsh0 + 265);
    const auto *gsh0_266 = buffer.data(gsh0 + 266);
    const auto *gsh0_267 = buffer.data(gsh0 + 267);
    const auto *gsh0_268 = buffer.data(gsh0 + 268);
    const auto *gsh0_269 = buffer.data(gsh0 + 269);
    const auto *gsh0_270 = buffer.data(gsh0 + 270);
    const auto *gsh0_271 = buffer.data(gsh0 + 271);
    const auto *gsh0_272 = buffer.data(gsh0 + 272);
    const auto *gsh0_274 = buffer.data(gsh0 + 274);
    const auto *gsh0_276 = buffer.data(gsh0 + 276);
    const auto *gsh0_277 = buffer.data(gsh0 + 277);
    const auto *gsh0_279 = buffer.data(gsh0 + 279);
    const auto *gsh0_280 = buffer.data(gsh0 + 280);
    const auto *gsh0_281 = buffer.data(gsh0 + 281);
    const auto *gsh0_283 = buffer.data(gsh0 + 283);
    const auto *gsh0_284 = buffer.data(gsh0 + 284);
    const auto *gsh0_285 = buffer.data(gsh0 + 285);
    const auto *gsh0_286 = buffer.data(gsh0 + 286);
    const auto *gsh0_288 = buffer.data(gsh0 + 288);
    const auto *gsh0_289 = buffer.data(gsh0 + 289);
    const auto *gsh0_290 = buffer.data(gsh0 + 290);
    const auto *gsh0_291 = buffer.data(gsh0 + 291);
    const auto *gsh0_292 = buffer.data(gsh0 + 292);

    const auto *gsh1_223 = buffer.data(gsh1 + 223);
    const auto *gsh1_224 = buffer.data(gsh1 + 224);
    const auto *gsh1_225 = buffer.data(gsh1 + 225);
    const auto *gsh1_226 = buffer.data(gsh1 + 226);
    const auto *gsh1_227 = buffer.data(gsh1 + 227);
    const auto *gsh1_228 = buffer.data(gsh1 + 228);
    const auto *gsh1_229 = buffer.data(gsh1 + 229);
    const auto *gsh1_230 = buffer.data(gsh1 + 230);
    const auto *gsh1_233 = buffer.data(gsh1 + 233);
    const auto *gsh1_235 = buffer.data(gsh1 + 235);
    const auto *gsh1_236 = buffer.data(gsh1 + 236);
    const auto *gsh1_238 = buffer.data(gsh1 + 238);
    const auto *gsh1_239 = buffer.data(gsh1 + 239);
    const auto *gsh1_240 = buffer.data(gsh1 + 240);
    const auto *gsh1_242 = buffer.data(gsh1 + 242);
    const auto *gsh1_243 = buffer.data(gsh1 + 243);
    const auto *gsh1_244 = buffer.data(gsh1 + 244);
    const auto *gsh1_245 = buffer.data(gsh1 + 245);
    const auto *gsh1_247 = buffer.data(gsh1 + 247);
    const auto *gsh1_248 = buffer.data(gsh1 + 248);
    const auto *gsh1_249 = buffer.data(gsh1 + 249);
    const auto *gsh1_250 = buffer.data(gsh1 + 250);
    const auto *gsh1_251 = buffer.data(gsh1 + 251);
    const auto *gsh1_252 = buffer.data(gsh1 + 252);
    const auto *gsh1_253 = buffer.data(gsh1 + 253);
    const auto *gsh1_254 = buffer.data(gsh1 + 254);
    const auto *gsh1_255 = buffer.data(gsh1 + 255);
    const auto *gsh1_256 = buffer.data(gsh1 + 256);
    const auto *gsh1_257 = buffer.data(gsh1 + 257);
    const auto *gsh1_258 = buffer.data(gsh1 + 258);
    const auto *gsh1_259 = buffer.data(gsh1 + 259);
    const auto *gsh1_260 = buffer.data(gsh1 + 260);
    const auto *gsh1_261 = buffer.data(gsh1 + 261);
    const auto *gsh1_262 = buffer.data(gsh1 + 262);
    const auto *gsh1_263 = buffer.data(gsh1 + 263);
    const auto *gsh1_264 = buffer.data(gsh1 + 264);
    const auto *gsh1_265 = buffer.data(gsh1 + 265);
    const auto *gsh1_266 = buffer.data(gsh1 + 266);
    const auto *gsh1_267 = buffer.data(gsh1 + 267);
    const auto *gsh1_268 = buffer.data(gsh1 + 268);
    const auto *gsh1_269 = buffer.data(gsh1 + 269);
    const auto *gsh1_270 = buffer.data(gsh1 + 270);
    const auto *gsh1_271 = buffer.data(gsh1 + 271);
    const auto *gsh1_272 = buffer.data(gsh1 + 272);
    const auto *gsh1_274 = buffer.data(gsh1 + 274);
    const auto *gsh1_276 = buffer.data(gsh1 + 276);
    const auto *gsh1_277 = buffer.data(gsh1 + 277);
    const auto *gsh1_279 = buffer.data(gsh1 + 279);
    const auto *gsh1_280 = buffer.data(gsh1 + 280);
    const auto *gsh1_281 = buffer.data(gsh1 + 281);
    const auto *gsh1_283 = buffer.data(gsh1 + 283);
    const auto *gsh1_284 = buffer.data(gsh1 + 284);
    const auto *gsh1_285 = buffer.data(gsh1 + 285);
    const auto *gsh1_286 = buffer.data(gsh1 + 286);
    const auto *gsh1_288 = buffer.data(gsh1 + 288);
    const auto *gsh1_289 = buffer.data(gsh1 + 289);
    const auto *gsh1_290 = buffer.data(gsh1 + 290);
    const auto *gsh1_291 = buffer.data(gsh1 + 291);
    const auto *gsh1_292 = buffer.data(gsh1 + 292);

    const auto *gsi_290 = buffer.data(gsi + 290);
    const auto *gsi_293 = buffer.data(gsi + 293);
    const auto *gsi_294 = buffer.data(gsi + 294);
    const auto *gsi_295 = buffer.data(gsi + 295);
    const auto *gsi_297 = buffer.data(gsi + 297);
    const auto *gsi_298 = buffer.data(gsi + 298);
    const auto *gsi_299 = buffer.data(gsi + 299);
    const auto *gsi_300 = buffer.data(gsi + 300);
    const auto *gsi_301 = buffer.data(gsi + 301);
    const auto *gsi_302 = buffer.data(gsi + 302);
    const auto *gsi_303 = buffer.data(gsi + 303);
    const auto *gsi_304 = buffer.data(gsi + 304);
    const auto *gsi_305 = buffer.data(gsi + 305);
    const auto *gsi_306 = buffer.data(gsi + 306);
    const auto *gsi_307 = buffer.data(gsi + 307);
    const auto *gsi_310 = buffer.data(gsi + 310);
    const auto *gsi_312 = buffer.data(gsi + 312);
    const auto *gsi_313 = buffer.data(gsi + 313);
    const auto *gsi_315 = buffer.data(gsi + 315);
    const auto *gsi_316 = buffer.data(gsi + 316);
    const auto *gsi_317 = buffer.data(gsi + 317);
    const auto *gsi_319 = buffer.data(gsi + 319);
    const auto *gsi_320 = buffer.data(gsi + 320);
    const auto *gsi_321 = buffer.data(gsi + 321);
    const auto *gsi_322 = buffer.data(gsi + 322);
    const auto *gsi_324 = buffer.data(gsi + 324);
    const auto *gsi_325 = buffer.data(gsi + 325);
    const auto *gsi_326 = buffer.data(gsi + 326);
    const auto *gsi_327 = buffer.data(gsi + 327);
    const auto *gsi_328 = buffer.data(gsi + 328);
    const auto *gsi_329 = buffer.data(gsi + 329);
    const auto *gsi_330 = buffer.data(gsi + 330);
    const auto *gsi_331 = buffer.data(gsi + 331);
    const auto *gsi_332 = buffer.data(gsi + 332);
    const auto *gsi_333 = buffer.data(gsi + 333);
    const auto *gsi_334 = buffer.data(gsi + 334);
    const auto *gsi_335 = buffer.data(gsi + 335);
    const auto *gsi_336 = buffer.data(gsi + 336);
    const auto *gsi_337 = buffer.data(gsi + 337);
    const auto *gsi_338 = buffer.data(gsi + 338);
    const auto *gsi_339 = buffer.data(gsi + 339);
    const auto *gsi_340 = buffer.data(gsi + 340);
    const auto *gsi_341 = buffer.data(gsi + 341);
    const auto *gsi_342 = buffer.data(gsi + 342);
    const auto *gsi_343 = buffer.data(gsi + 343);
    const auto *gsi_344 = buffer.data(gsi + 344);
    const auto *gsi_345 = buffer.data(gsi + 345);
    const auto *gsi_346 = buffer.data(gsi + 346);
    const auto *gsi_347 = buffer.data(gsi + 347);
    const auto *gsi_348 = buffer.data(gsi + 348);
    const auto *gsi_349 = buffer.data(gsi + 349);
    const auto *gsi_350 = buffer.data(gsi + 350);
    const auto *gsi_351 = buffer.data(gsi + 351);
    const auto *gsi_352 = buffer.data(gsi + 352);
    const auto *gsi_353 = buffer.data(gsi + 353);
    const auto *gsi_354 = buffer.data(gsi + 354);
    const auto *gsi_355 = buffer.data(gsi + 355);
    const auto *gsi_356 = buffer.data(gsi + 356);
    const auto *gsi_357 = buffer.data(gsi + 357);
    const auto *gsi_358 = buffer.data(gsi + 358);
    const auto *gsi_359 = buffer.data(gsi + 359);
    const auto *gsi_360 = buffer.data(gsi + 360);
    const auto *gsi_361 = buffer.data(gsi + 361);
    const auto *gsi_362 = buffer.data(gsi + 362);
    const auto *gsi_363 = buffer.data(gsi + 363);
    const auto *gsi_365 = buffer.data(gsi + 365);
    const auto *gsi_367 = buffer.data(gsi + 367);
    const auto *gsi_368 = buffer.data(gsi + 368);
    const auto *gsi_370 = buffer.data(gsi + 370);
    const auto *gsi_371 = buffer.data(gsi + 371);
    const auto *gsi_372 = buffer.data(gsi + 372);
    const auto *gsi_374 = buffer.data(gsi + 374);
    const auto *gsi_375 = buffer.data(gsi + 375);
    const auto *gsi_376 = buffer.data(gsi + 376);
    const auto *gsi_377 = buffer.data(gsi + 377);
    const auto *gsi_379 = buffer.data(gsi + 379);
    const auto *gsi_380 = buffer.data(gsi + 380);
    const auto *gsi_381 = buffer.data(gsi + 381);
    const auto *gsi_382 = buffer.data(gsi + 382);
    const auto *gsi_383 = buffer.data(gsi + 383);
    const auto *gsi_385 = buffer.data(gsi + 385);
    const auto *gsi_386 = buffer.data(gsi + 386);
    const auto *gsi_387 = buffer.data(gsi + 387);

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pc_x, pc_z, gsh0_223, gsh0_224, gsh0_225, \
                         gsh1_223, gsh1_224, gsh1_225, gsi_290, gsi_293, gsi_294, \
                         gsi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_6 * gsh0_223[k]
                   - f_7 * gsh1_223[k]
                   + f_3 * pc_x[k] * gsi_293[k];

        t_374[k] = f_6 * gsh0_224[k]
                   - f_7 * gsh1_224[k]
                   + f_3 * pc_x[k] * gsi_294[k];

        t_375[k] = f_4 * gsh0_225[k]
                   - f_5 * gsh1_225[k]
                   + f_3 * pc_x[k] * gsi_295[k];

        t_376[k] = f_3 * pc_z[k] * gsi_290[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pc_x, gsh0_227, gsh0_228, gsh0_229, gsh1_227, \
                         gsh1_228, gsh1_229, gsi_297, gsi_298, \
                         gsi_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_4 * gsh0_227[k]
                   - f_5 * gsh1_227[k]
                   + f_3 * pc_x[k] * gsi_297[k];

        t_378[k] = f_4 * gsh0_228[k]
                   - f_5 * gsh1_228[k]
                   + f_3 * pc_x[k] * gsi_298[k];

        t_379[k] = f_4 * gsh0_229[k]
                   - f_5 * gsh1_229[k]
                   + f_3 * pc_x[k] * gsi_299[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, t_385, pc_x, gsh0_230, gsh1_230, \
                         gsi_300, gsi_301, gsi_302, gsi_303, gsi_304, \
                         gsi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_4 * gsh0_230[k]
                   - f_5 * gsh1_230[k]
                   + f_3 * pc_x[k] * gsi_300[k];

        t_381[k] = f_3 * pc_x[k] * gsi_301[k];

        t_382[k] = f_3 * pc_x[k] * gsi_302[k];

        t_383[k] = f_3 * pc_x[k] * gsi_303[k];

        t_384[k] = f_3 * pc_x[k] * gsi_304[k];

        t_385[k] = f_3 * pc_x[k] * gsi_305[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, pc_x, pc_y, pc_z, fsi_189, \
                         gsh0_225, gsh1_225, gsi_301, gsi_302, gsi_306, \
                         gsi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_3 * pc_x[k] * gsi_306[k];

        t_387[k] = f_3 * pc_x[k] * gsi_307[k];

        t_388[k] = f_0 * fsi_189[k]
                   + f_1 * gsh0_225[k]
                   - f_2 * gsh1_225[k]
                   + f_3 * pc_y[k] * gsi_301[k];

        t_389[k] = f_3 * pc_z[k] * gsi_301[k];

        t_390[k] = f_4 * gsh0_225[k]
                   - f_5 * gsh1_225[k]
                   + f_3 * pc_z[k] * gsi_302[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_z, gsh0_226, gsh0_227, gsh0_228, gsh1_226, \
                         gsh1_227, gsh1_228, gsi_303, gsi_304, \
                         gsi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_6 * gsh0_226[k]
                   - f_7 * gsh1_226[k]
                   + f_3 * pc_z[k] * gsi_303[k];

        t_392[k] = f_8 * gsh0_227[k]
                   - f_9 * gsh1_227[k]
                   + f_3 * pc_z[k] * gsi_304[k];

        t_393[k] = f_10 * gsh0_228[k]
                   - f_11 * gsh1_228[k]
                   + f_3 * pc_z[k] * gsi_305[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pa_z, pc_y, pc_z, fsk0_216, fsk0_217, \
                         fsi_195, fsk1_216, fsk1_217, gsh0_230, gsh1_230, \
                         gsi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_0 * fsi_195[k]
                   + f_3 * pc_y[k] * gsi_307[k];

        t_395[k] = f_1 * gsh0_230[k]
                   - f_2 * gsh1_230[k]
                   + f_3 * pc_z[k] * gsi_307[k];

        t_396[k] = pa_z[k] * fsk0_216[k]
                   - f_12 * pc_z[k] * fsk1_216[k];

        t_397[k] = pa_z[k] * fsk0_217[k]
                   - f_12 * pc_z[k] * fsk1_217[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_z, pc_x, pc_z, fsk0_219, fsk1_219, gsh0_233, \
                         gsh0_235, gsh1_233, gsh1_235, gsi_310, \
                         gsi_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_17 * gsh0_233[k]
                   - f_18 * gsh1_233[k]
                   + f_3 * pc_x[k] * gsi_310[k];

        t_399[k] = pa_z[k] * fsk0_219[k]
                   - f_12 * pc_z[k] * fsk1_219[k];

        t_400[k] = f_10 * gsh0_235[k]
                   - f_11 * gsh1_235[k]
                   + f_3 * pc_x[k] * gsi_312[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_z, pc_x, pc_z, fsk0_222, fsk1_222, gsh0_236, \
                         gsh0_238, gsh1_236, gsh1_238, gsi_313, \
                         gsi_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_10 * gsh0_236[k]
                   - f_11 * gsh1_236[k]
                   + f_3 * pc_x[k] * gsi_313[k];

        t_402[k] = pa_z[k] * fsk0_222[k]
                   - f_12 * pc_z[k] * fsk1_222[k];

        t_403[k] = f_8 * gsh0_238[k]
                   - f_9 * gsh1_238[k]
                   + f_3 * pc_x[k] * gsi_315[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_z, pc_x, pc_z, fsk0_226, fsk1_226, gsh0_239, \
                         gsh0_240, gsh1_239, gsh1_240, gsi_316, \
                         gsi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_8 * gsh0_239[k]
                   - f_9 * gsh1_239[k]
                   + f_3 * pc_x[k] * gsi_316[k];

        t_405[k] = f_8 * gsh0_240[k]
                   - f_9 * gsh1_240[k]
                   + f_3 * pc_x[k] * gsi_317[k];

        t_406[k] = pa_z[k] * fsk0_226[k]
                   - f_12 * pc_z[k] * fsk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pc_x, gsh0_242, gsh0_243, gsh0_244, gsh1_242, \
                         gsh1_243, gsh1_244, gsi_319, gsi_320, \
                         gsi_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_6 * gsh0_242[k]
                   - f_7 * gsh1_242[k]
                   + f_3 * pc_x[k] * gsi_319[k];

        t_408[k] = f_6 * gsh0_243[k]
                   - f_7 * gsh1_243[k]
                   + f_3 * pc_x[k] * gsi_320[k];

        t_409[k] = f_6 * gsh0_244[k]
                   - f_7 * gsh1_244[k]
                   + f_3 * pc_x[k] * gsi_321[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pa_z, pc_x, pc_z, fsk0_231, fsk1_231, gsh0_245, \
                         gsh0_247, gsh1_245, gsh1_247, gsi_322, \
                         gsi_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_6 * gsh0_245[k]
                   - f_7 * gsh1_245[k]
                   + f_3 * pc_x[k] * gsi_322[k];

        t_411[k] = pa_z[k] * fsk0_231[k]
                   - f_12 * pc_z[k] * fsk1_231[k];

        t_412[k] = f_4 * gsh0_247[k]
                   - f_5 * gsh1_247[k]
                   + f_3 * pc_x[k] * gsi_324[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pc_x, gsh0_248, gsh0_249, gsh0_250, gsh1_248, \
                         gsh1_249, gsh1_250, gsi_325, gsi_326, \
                         gsi_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_4 * gsh0_248[k]
                   - f_5 * gsh1_248[k]
                   + f_3 * pc_x[k] * gsi_325[k];

        t_414[k] = f_4 * gsh0_249[k]
                   - f_5 * gsh1_249[k]
                   + f_3 * pc_x[k] * gsi_326[k];

        t_415[k] = f_4 * gsh0_250[k]
                   - f_5 * gsh1_250[k]
                   + f_3 * pc_x[k] * gsi_327[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, t_421, pc_x, gsh0_251, gsh1_251, \
                         gsi_328, gsi_329, gsi_330, gsi_331, gsi_332, \
                         gsi_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_4 * gsh0_251[k]
                   - f_5 * gsh1_251[k]
                   + f_3 * pc_x[k] * gsi_328[k];

        t_417[k] = f_3 * pc_x[k] * gsi_329[k];

        t_418[k] = f_3 * pc_x[k] * gsi_330[k];

        t_419[k] = f_3 * pc_x[k] * gsi_331[k];

        t_420[k] = f_3 * pc_x[k] * gsi_332[k];

        t_421[k] = f_3 * pc_x[k] * gsi_333[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pa_z, pc_x, pc_z, fsk0_244, fsi_189, \
                         fsk1_244, gsi_329, gsi_334, gsi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_3 * pc_x[k] * gsi_334[k];

        t_423[k] = f_3 * pc_x[k] * gsi_335[k];

        t_424[k] = pa_z[k] * fsk0_244[k]
                   - f_12 * pc_z[k] * fsk1_244[k];

        t_425[k] = f_13 * fsi_189[k]
                   + f_3 * pc_z[k] * gsi_329[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pa_z, pc_z, fsk0_246, fsk0_247, fsk0_248, \
                         fsi_190, fsi_191, fsi_192, fsk1_246, fsk1_247, \
                         fsk1_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = pa_z[k] * fsk0_246[k]
                   + f_14 * fsi_190[k]
                   - f_12 * pc_z[k] * fsk1_246[k];

        t_427[k] = pa_z[k] * fsk0_247[k]
                   + f_15 * fsi_191[k]
                   - f_12 * pc_z[k] * fsk1_247[k];

        t_428[k] = pa_z[k] * fsk0_248[k]
                   + f_0 * fsi_192[k]
                   - f_12 * pc_z[k] * fsk1_248[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pa_z, pc_y, pc_z, fsk0_249, fsi_193, fsi_195, \
                         fsi_223, fsk1_249, gsh0_251, gsh1_251, \
                         gsi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = pa_z[k] * fsk0_249[k]
                   + f_16 * fsi_193[k]
                   - f_12 * pc_z[k] * fsk1_249[k];

        t_430[k] = f_15 * fsi_223[k]
                   + f_3 * pc_y[k] * gsi_335[k];

        t_431[k] = f_13 * fsi_195[k]
                   + f_1 * gsh0_251[k]
                   - f_2 * gsh1_251[k]
                   + f_3 * pc_z[k] * gsi_335[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, pc_x, gsh0_252, gsh0_253, gsh0_254, gsh1_252, \
                         gsh1_253, gsh1_254, gsi_336, gsi_337, \
                         gsi_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_1 * gsh0_252[k]
                   - f_2 * gsh1_252[k]
                   + f_3 * pc_x[k] * gsi_336[k];

        t_433[k] = f_17 * gsh0_253[k]
                   - f_18 * gsh1_253[k]
                   + f_3 * pc_x[k] * gsi_337[k];

        t_434[k] = f_17 * gsh0_254[k]
                   - f_18 * gsh1_254[k]
                   + f_3 * pc_x[k] * gsi_338[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, pc_x, gsh0_255, gsh0_256, gsh0_257, gsh1_255, \
                         gsh1_256, gsh1_257, gsi_339, gsi_340, \
                         gsi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_10 * gsh0_255[k]
                   - f_11 * gsh1_255[k]
                   + f_3 * pc_x[k] * gsi_339[k];

        t_436[k] = f_10 * gsh0_256[k]
                   - f_11 * gsh1_256[k]
                   + f_3 * pc_x[k] * gsi_340[k];

        t_437[k] = f_10 * gsh0_257[k]
                   - f_11 * gsh1_257[k]
                   + f_3 * pc_x[k] * gsi_341[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pc_x, gsh0_258, gsh0_259, gsh0_260, gsh1_258, \
                         gsh1_259, gsh1_260, gsi_342, gsi_343, \
                         gsi_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_8 * gsh0_258[k]
                   - f_9 * gsh1_258[k]
                   + f_3 * pc_x[k] * gsi_342[k];

        t_439[k] = f_8 * gsh0_259[k]
                   - f_9 * gsh1_259[k]
                   + f_3 * pc_x[k] * gsi_343[k];

        t_440[k] = f_8 * gsh0_260[k]
                   - f_9 * gsh1_260[k]
                   + f_3 * pc_x[k] * gsi_344[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_x, gsh0_261, gsh0_262, gsh0_263, gsh1_261, \
                         gsh1_262, gsh1_263, gsi_345, gsi_346, \
                         gsi_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_8 * gsh0_261[k]
                   - f_9 * gsh1_261[k]
                   + f_3 * pc_x[k] * gsi_345[k];

        t_442[k] = f_6 * gsh0_262[k]
                   - f_7 * gsh1_262[k]
                   + f_3 * pc_x[k] * gsi_346[k];

        t_443[k] = f_6 * gsh0_263[k]
                   - f_7 * gsh1_263[k]
                   + f_3 * pc_x[k] * gsi_347[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_x, gsh0_264, gsh0_265, gsh0_266, gsh1_264, \
                         gsh1_265, gsh1_266, gsi_348, gsi_349, \
                         gsi_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_6 * gsh0_264[k]
                   - f_7 * gsh1_264[k]
                   + f_3 * pc_x[k] * gsi_348[k];

        t_445[k] = f_6 * gsh0_265[k]
                   - f_7 * gsh1_265[k]
                   + f_3 * pc_x[k] * gsi_349[k];

        t_446[k] = f_6 * gsh0_266[k]
                   - f_7 * gsh1_266[k]
                   + f_3 * pc_x[k] * gsi_350[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_x, gsh0_267, gsh0_268, gsh0_269, gsh1_267, \
                         gsh1_268, gsh1_269, gsi_351, gsi_352, \
                         gsi_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_4 * gsh0_267[k]
                   - f_5 * gsh1_267[k]
                   + f_3 * pc_x[k] * gsi_351[k];

        t_448[k] = f_4 * gsh0_268[k]
                   - f_5 * gsh1_268[k]
                   + f_3 * pc_x[k] * gsi_352[k];

        t_449[k] = f_4 * gsh0_269[k]
                   - f_5 * gsh1_269[k]
                   + f_3 * pc_x[k] * gsi_353[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pc_x, gsh0_270, gsh0_271, gsh0_272, \
                         gsh1_270, gsh1_271, gsh1_272, gsi_354, gsi_355, gsi_356, \
                         gsi_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_4 * gsh0_270[k]
                   - f_5 * gsh1_270[k]
                   + f_3 * pc_x[k] * gsi_354[k];

        t_451[k] = f_4 * gsh0_271[k]
                   - f_5 * gsh1_271[k]
                   + f_3 * pc_x[k] * gsi_355[k];

        t_452[k] = f_4 * gsh0_272[k]
                   - f_5 * gsh1_272[k]
                   + f_3 * pc_x[k] * gsi_356[k];

        t_453[k] = f_3 * pc_x[k] * gsi_357[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, t_459, pc_x, gsi_358, gsi_359, \
                         gsi_360, gsi_361, gsi_362, gsi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_3 * pc_x[k] * gsi_358[k];

        t_455[k] = f_3 * pc_x[k] * gsi_359[k];

        t_456[k] = f_3 * pc_x[k] * gsi_360[k];

        t_457[k] = f_3 * pc_x[k] * gsi_361[k];

        t_458[k] = f_3 * pc_x[k] * gsi_362[k];

        t_459[k] = f_3 * pc_x[k] * gsi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, fsi_217, fsi_245, fsi_247, gsh0_267, \
                         gsh0_269, gsh1_267, gsh1_269, gsi_357, \
                         gsi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * fsi_245[k]
                   + f_1 * gsh0_267[k]
                   - f_2 * gsh1_267[k]
                   + f_3 * pc_y[k] * gsi_357[k];

        t_461[k] = f_14 * fsi_217[k]
                   + f_3 * pc_z[k] * gsi_357[k];

        t_462[k] = f_14 * fsi_247[k]
                   + f_10 * gsh0_269[k]
                   - f_11 * gsh1_269[k]
                   + f_3 * pc_y[k] * gsi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, fsi_248, fsi_249, fsi_250, gsh0_270, \
                         gsh0_271, gsh0_272, gsh1_270, gsh1_271, gsh1_272, gsi_360, gsi_361, \
                         gsi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * fsi_248[k]
                   + f_8 * gsh0_270[k]
                   - f_9 * gsh1_270[k]
                   + f_3 * pc_y[k] * gsi_360[k];

        t_464[k] = f_14 * fsi_249[k]
                   + f_6 * gsh0_271[k]
                   - f_7 * gsh1_271[k]
                   + f_3 * pc_y[k] * gsi_361[k];

        t_465[k] = f_14 * fsi_250[k]
                   + f_4 * gsh0_272[k]
                   - f_5 * gsh1_272[k]
                   + f_3 * pc_y[k] * gsi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, pa_y, pc_y, pc_z, fsk0_324, fsi_223, fsi_251, \
                         fsk1_324, gsh0_272, gsh1_272, gsi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * fsi_251[k]
                   + f_3 * pc_y[k] * gsi_363[k];

        t_467[k] = f_14 * fsi_223[k]
                   + f_1 * gsh0_272[k]
                   - f_2 * gsh1_272[k]
                   + f_3 * pc_z[k] * gsi_363[k];

        t_468[k] = pa_y[k] * fsk0_324[k]
                   - f_12 * pc_y[k] * fsk1_324[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, pa_y, pc_x, pc_y, fsk0_326, fsk1_326, gsh0_274, \
                         gsh0_276, gsh1_274, gsh1_276, gsi_365, \
                         gsi_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_17 * gsh0_274[k]
                   - f_18 * gsh1_274[k]
                   + f_3 * pc_x[k] * gsi_365[k];

        t_470[k] = pa_y[k] * fsk0_326[k]
                   - f_12 * pc_y[k] * fsk1_326[k];

        t_471[k] = f_10 * gsh0_276[k]
                   - f_11 * gsh1_276[k]
                   + f_3 * pc_x[k] * gsi_367[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, pa_y, pc_x, pc_y, fsk0_329, fsk1_329, gsh0_277, \
                         gsh0_279, gsh1_277, gsh1_279, gsi_368, \
                         gsi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = f_10 * gsh0_277[k]
                   - f_11 * gsh1_277[k]
                   + f_3 * pc_x[k] * gsi_368[k];

        t_473[k] = pa_y[k] * fsk0_329[k]
                   - f_12 * pc_y[k] * fsk1_329[k];

        t_474[k] = f_8 * gsh0_279[k]
                   - f_9 * gsh1_279[k]
                   + f_3 * pc_x[k] * gsi_370[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, pa_y, pc_x, pc_y, fsk0_333, fsk1_333, gsh0_280, \
                         gsh0_281, gsh1_280, gsh1_281, gsi_371, \
                         gsi_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_8 * gsh0_280[k]
                   - f_9 * gsh1_280[k]
                   + f_3 * pc_x[k] * gsi_371[k];

        t_476[k] = f_8 * gsh0_281[k]
                   - f_9 * gsh1_281[k]
                   + f_3 * pc_x[k] * gsi_372[k];

        t_477[k] = pa_y[k] * fsk0_333[k]
                   - f_12 * pc_y[k] * fsk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_x, gsh0_283, gsh0_284, gsh0_285, gsh1_283, \
                         gsh1_284, gsh1_285, gsi_374, gsi_375, \
                         gsi_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_6 * gsh0_283[k]
                   - f_7 * gsh1_283[k]
                   + f_3 * pc_x[k] * gsi_374[k];

        t_479[k] = f_6 * gsh0_284[k]
                   - f_7 * gsh1_284[k]
                   + f_3 * pc_x[k] * gsi_375[k];

        t_480[k] = f_6 * gsh0_285[k]
                   - f_7 * gsh1_285[k]
                   + f_3 * pc_x[k] * gsi_376[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, pa_y, pc_x, pc_y, fsk0_338, fsk1_338, gsh0_286, \
                         gsh0_288, gsh1_286, gsh1_288, gsi_377, \
                         gsi_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_6 * gsh0_286[k]
                   - f_7 * gsh1_286[k]
                   + f_3 * pc_x[k] * gsi_377[k];

        t_482[k] = pa_y[k] * fsk0_338[k]
                   - f_12 * pc_y[k] * fsk1_338[k];

        t_483[k] = f_4 * gsh0_288[k]
                   - f_5 * gsh1_288[k]
                   + f_3 * pc_x[k] * gsi_379[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, pc_x, gsh0_289, gsh0_290, gsh0_291, gsh1_289, \
                         gsh1_290, gsh1_291, gsi_380, gsi_381, \
                         gsi_382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_4 * gsh0_289[k]
                   - f_5 * gsh1_289[k]
                   + f_3 * pc_x[k] * gsi_380[k];

        t_485[k] = f_4 * gsh0_290[k]
                   - f_5 * gsh1_290[k]
                   + f_3 * pc_x[k] * gsi_381[k];

        t_486[k] = f_4 * gsh0_291[k]
                   - f_5 * gsh1_291[k]
                   + f_3 * pc_x[k] * gsi_382[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, pa_y, pc_x, pc_y, fsk0_344, \
                         fsk1_344, gsh0_292, gsh1_292, gsi_383, gsi_385, gsi_386, \
                         gsi_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_4 * gsh0_292[k]
                   - f_5 * gsh1_292[k]
                   + f_3 * pc_x[k] * gsi_383[k];

        t_488[k] = pa_y[k] * fsk0_344[k]
                   - f_12 * pc_y[k] * fsk1_344[k];

        t_489[k] = f_3 * pc_x[k] * gsi_385[k];

        t_490[k] = f_3 * pc_x[k] * gsi_386[k];

        t_491[k] = f_3 * pc_x[k] * gsi_387[k];
    }
}

static auto
compute_prim_gsk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsk0,
                                                          const size_t fsi, const size_t fsk1,
                                                          const size_t gsh0, const size_t gsh1,
                                                          const size_t gsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsk0_352 = buffer.data(fsk0 + 352);
    const auto *fsk0_354 = buffer.data(fsk0 + 354);
    const auto *fsk0_355 = buffer.data(fsk0 + 355);
    const auto *fsk0_356 = buffer.data(fsk0 + 356);
    const auto *fsk0_357 = buffer.data(fsk0 + 357);
    const auto *fsk0_359 = buffer.data(fsk0 + 359);

    const auto *fsi_245 = buffer.data(fsi + 245);
    const auto *fsi_273 = buffer.data(fsi + 273);
    const auto *fsi_275 = buffer.data(fsi + 275);
    const auto *fsi_276 = buffer.data(fsi + 276);
    const auto *fsi_277 = buffer.data(fsi + 277);
    const auto *fsi_278 = buffer.data(fsi + 278);
    const auto *fsi_279 = buffer.data(fsi + 279);

    const auto *fsk1_352 = buffer.data(fsk1 + 352);
    const auto *fsk1_354 = buffer.data(fsk1 + 354);
    const auto *fsk1_355 = buffer.data(fsk1 + 355);
    const auto *fsk1_356 = buffer.data(fsk1 + 356);
    const auto *fsk1_357 = buffer.data(fsk1 + 357);
    const auto *fsk1_359 = buffer.data(fsk1 + 359);

    const auto *gsh0_294 = buffer.data(gsh0 + 294);
    const auto *gsh0_296 = buffer.data(gsh0 + 296);
    const auto *gsh0_297 = buffer.data(gsh0 + 297);
    const auto *gsh0_299 = buffer.data(gsh0 + 299);
    const auto *gsh0_300 = buffer.data(gsh0 + 300);
    const auto *gsh0_301 = buffer.data(gsh0 + 301);
    const auto *gsh0_303 = buffer.data(gsh0 + 303);
    const auto *gsh0_304 = buffer.data(gsh0 + 304);
    const auto *gsh0_305 = buffer.data(gsh0 + 305);
    const auto *gsh0_306 = buffer.data(gsh0 + 306);
    const auto *gsh0_308 = buffer.data(gsh0 + 308);
    const auto *gsh0_309 = buffer.data(gsh0 + 309);
    const auto *gsh0_310 = buffer.data(gsh0 + 310);
    const auto *gsh0_311 = buffer.data(gsh0 + 311);
    const auto *gsh0_312 = buffer.data(gsh0 + 312);
    const auto *gsh0_313 = buffer.data(gsh0 + 313);
    const auto *gsh0_314 = buffer.data(gsh0 + 314);

    const auto *gsh1_294 = buffer.data(gsh1 + 294);
    const auto *gsh1_296 = buffer.data(gsh1 + 296);
    const auto *gsh1_297 = buffer.data(gsh1 + 297);
    const auto *gsh1_299 = buffer.data(gsh1 + 299);
    const auto *gsh1_300 = buffer.data(gsh1 + 300);
    const auto *gsh1_301 = buffer.data(gsh1 + 301);
    const auto *gsh1_303 = buffer.data(gsh1 + 303);
    const auto *gsh1_304 = buffer.data(gsh1 + 304);
    const auto *gsh1_305 = buffer.data(gsh1 + 305);
    const auto *gsh1_306 = buffer.data(gsh1 + 306);
    const auto *gsh1_308 = buffer.data(gsh1 + 308);
    const auto *gsh1_309 = buffer.data(gsh1 + 309);
    const auto *gsh1_310 = buffer.data(gsh1 + 310);
    const auto *gsh1_311 = buffer.data(gsh1 + 311);
    const auto *gsh1_312 = buffer.data(gsh1 + 312);
    const auto *gsh1_313 = buffer.data(gsh1 + 313);
    const auto *gsh1_314 = buffer.data(gsh1 + 314);

    const auto *gsi_385 = buffer.data(gsi + 385);
    const auto *gsi_388 = buffer.data(gsi + 388);
    const auto *gsi_389 = buffer.data(gsi + 389);
    const auto *gsi_390 = buffer.data(gsi + 390);
    const auto *gsi_391 = buffer.data(gsi + 391);
    const auto *gsi_392 = buffer.data(gsi + 392);
    const auto *gsi_394 = buffer.data(gsi + 394);
    const auto *gsi_395 = buffer.data(gsi + 395);
    const auto *gsi_397 = buffer.data(gsi + 397);
    const auto *gsi_398 = buffer.data(gsi + 398);
    const auto *gsi_399 = buffer.data(gsi + 399);
    const auto *gsi_401 = buffer.data(gsi + 401);
    const auto *gsi_402 = buffer.data(gsi + 402);
    const auto *gsi_403 = buffer.data(gsi + 403);
    const auto *gsi_404 = buffer.data(gsi + 404);
    const auto *gsi_406 = buffer.data(gsi + 406);
    const auto *gsi_407 = buffer.data(gsi + 407);
    const auto *gsi_408 = buffer.data(gsi + 408);
    const auto *gsi_409 = buffer.data(gsi + 409);
    const auto *gsi_410 = buffer.data(gsi + 410);
    const auto *gsi_412 = buffer.data(gsi + 412);
    const auto *gsi_413 = buffer.data(gsi + 413);
    const auto *gsi_414 = buffer.data(gsi + 414);
    const auto *gsi_415 = buffer.data(gsi + 415);
    const auto *gsi_416 = buffer.data(gsi + 416);
    const auto *gsi_417 = buffer.data(gsi + 417);
    const auto *gsi_418 = buffer.data(gsi + 418);
    const auto *gsi_419 = buffer.data(gsi + 419);

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, pa_y, pc_x, pc_y, fsk0_352, \
                         fsi_273, fsk1_352, gsi_388, gsi_389, gsi_390, \
                         gsi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_3 * pc_x[k] * gsi_388[k];

        t_493[k] = f_3 * pc_x[k] * gsi_389[k];

        t_494[k] = f_3 * pc_x[k] * gsi_390[k];

        t_495[k] = f_3 * pc_x[k] * gsi_391[k];

        t_496[k] = pa_y[k] * fsk0_352[k]
                   + f_19 * fsi_273[k]
                   - f_12 * pc_y[k] * fsk1_352[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, pa_y, pc_y, pc_z, fsk0_354, fsk0_355, fsi_245, \
                         fsi_275, fsi_276, fsk1_354, fsk1_355, \
                         gsi_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = f_15 * fsi_245[k]
                   + f_3 * pc_z[k] * gsi_385[k];

        t_498[k] = pa_y[k] * fsk0_354[k]
                   + f_16 * fsi_275[k]
                   - f_12 * pc_y[k] * fsk1_354[k];

        t_499[k] = pa_y[k] * fsk0_355[k]
                   + f_0 * fsi_276[k]
                   - f_12 * pc_y[k] * fsk1_355[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_y, pc_y, fsk0_356, fsk0_357, fsk0_359, \
                         fsi_277, fsi_278, fsi_279, fsk1_356, fsk1_357, fsk1_359, \
                         gsi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pa_y[k] * fsk0_356[k]
                   + f_15 * fsi_277[k]
                   - f_12 * pc_y[k] * fsk1_356[k];

        t_501[k] = pa_y[k] * fsk0_357[k]
                   + f_14 * fsi_278[k]
                   - f_12 * pc_y[k] * fsk1_357[k];

        t_502[k] = f_13 * fsi_279[k]
                   + f_3 * pc_y[k] * gsi_391[k];

        t_503[k] = pa_y[k] * fsk0_359[k]
                   - f_12 * pc_y[k] * fsk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, pc_x, pc_y, gsh0_294, gsh0_296, \
                         gsh0_297, gsh1_294, gsh1_296, gsh1_297, gsi_392, gsi_394, \
                         gsi_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_1 * gsh0_294[k]
                   - f_2 * gsh1_294[k]
                   + f_3 * pc_x[k] * gsi_392[k];

        t_505[k] = f_3 * pc_y[k] * gsi_392[k];

        t_506[k] = f_17 * gsh0_296[k]
                   - f_18 * gsh1_296[k]
                   + f_3 * pc_x[k] * gsi_394[k];

        t_507[k] = f_10 * gsh0_297[k]
                   - f_11 * gsh1_297[k]
                   + f_3 * pc_x[k] * gsi_395[k];

        t_508[k] = f_3 * pc_y[k] * gsi_394[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, pc_x, pc_y, gsh0_299, gsh0_300, gsh0_301, \
                         gsh1_299, gsh1_300, gsh1_301, gsi_397, gsi_398, \
                         gsi_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_10 * gsh0_299[k]
                   - f_11 * gsh1_299[k]
                   + f_3 * pc_x[k] * gsi_397[k];

        t_510[k] = f_8 * gsh0_300[k]
                   - f_9 * gsh1_300[k]
                   + f_3 * pc_x[k] * gsi_398[k];

        t_511[k] = f_8 * gsh0_301[k]
                   - f_9 * gsh1_301[k]
                   + f_3 * pc_x[k] * gsi_399[k];

        t_512[k] = f_3 * pc_y[k] * gsi_397[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, gsh0_303, gsh0_304, gsh0_305, gsh1_303, \
                         gsh1_304, gsh1_305, gsi_401, gsi_402, \
                         gsi_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_8 * gsh0_303[k]
                   - f_9 * gsh1_303[k]
                   + f_3 * pc_x[k] * gsi_401[k];

        t_514[k] = f_6 * gsh0_304[k]
                   - f_7 * gsh1_304[k]
                   + f_3 * pc_x[k] * gsi_402[k];

        t_515[k] = f_6 * gsh0_305[k]
                   - f_7 * gsh1_305[k]
                   + f_3 * pc_x[k] * gsi_403[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, gsh0_306, gsh0_308, gsh0_309, \
                         gsh1_306, gsh1_308, gsh1_309, gsi_401, gsi_404, gsi_406, \
                         gsi_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_6 * gsh0_306[k]
                   - f_7 * gsh1_306[k]
                   + f_3 * pc_x[k] * gsi_404[k];

        t_517[k] = f_3 * pc_y[k] * gsi_401[k];

        t_518[k] = f_6 * gsh0_308[k]
                   - f_7 * gsh1_308[k]
                   + f_3 * pc_x[k] * gsi_406[k];

        t_519[k] = f_4 * gsh0_309[k]
                   - f_5 * gsh1_309[k]
                   + f_3 * pc_x[k] * gsi_407[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, pc_x, pc_y, gsh0_310, gsh0_311, gsh0_312, \
                         gsh1_310, gsh1_311, gsh1_312, gsi_406, gsi_408, gsi_409, \
                         gsi_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_4 * gsh0_310[k]
                   - f_5 * gsh1_310[k]
                   + f_3 * pc_x[k] * gsi_408[k];

        t_521[k] = f_4 * gsh0_311[k]
                   - f_5 * gsh1_311[k]
                   + f_3 * pc_x[k] * gsi_409[k];

        t_522[k] = f_4 * gsh0_312[k]
                   - f_5 * gsh1_312[k]
                   + f_3 * pc_x[k] * gsi_410[k];

        t_523[k] = f_3 * pc_y[k] * gsi_406[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, t_529, pc_x, gsh0_314, gsh1_314, \
                         gsi_412, gsi_413, gsi_414, gsi_415, gsi_416, \
                         gsi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_4 * gsh0_314[k]
                   - f_5 * gsh1_314[k]
                   + f_3 * pc_x[k] * gsi_412[k];

        t_525[k] = f_3 * pc_x[k] * gsi_413[k];

        t_526[k] = f_3 * pc_x[k] * gsi_414[k];

        t_527[k] = f_3 * pc_x[k] * gsi_415[k];

        t_528[k] = f_3 * pc_x[k] * gsi_416[k];

        t_529[k] = f_3 * pc_x[k] * gsi_417[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pc_x, pc_y, gsh0_309, gsh0_310, gsh1_309, \
                         gsh1_310, gsi_413, gsi_414, gsi_418, gsi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_3 * pc_x[k] * gsi_418[k];

        t_531[k] = f_3 * pc_x[k] * gsi_419[k];

        t_532[k] = f_1 * gsh0_309[k]
                   - f_2 * gsh1_309[k]
                   + f_3 * pc_y[k] * gsi_413[k];

        t_533[k] = f_17 * gsh0_310[k]
                   - f_18 * gsh1_310[k]
                   + f_3 * pc_y[k] * gsi_414[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_y, gsh0_311, gsh0_312, gsh0_313, gsh1_311, \
                         gsh1_312, gsh1_313, gsi_415, gsi_416, \
                         gsi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_10 * gsh0_311[k]
                   - f_11 * gsh1_311[k]
                   + f_3 * pc_y[k] * gsi_415[k];

        t_535[k] = f_8 * gsh0_312[k]
                   - f_9 * gsh1_312[k]
                   + f_3 * pc_y[k] * gsi_416[k];

        t_536[k] = f_6 * gsh0_313[k]
                   - f_7 * gsh1_313[k]
                   + f_3 * pc_y[k] * gsi_417[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pc_y, pc_z, fsi_279, gsh0_314, gsh1_314, \
                         gsi_418, gsi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * gsh0_314[k]
                   - f_5 * gsh1_314[k]
                   + f_3 * pc_y[k] * gsi_418[k];

        t_538[k] = f_3 * pc_y[k] * gsi_419[k];

        t_539[k] = f_0 * fsi_279[k]
                   + f_1 * gsh0_314[k]
                   - f_2 * gsh1_314[k]
                   + f_3 * pc_z[k] * gsi_419[k];
    }
}

auto
compute_prim_gsk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fsk0, const size_t fsi,
                                                   const size_t fsk1, const size_t gsh0,
                                                   const size_t gsh1, const size_t gsi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gsk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, fsk0, fsi,
                                                              fsk1, gsh0, gsh1, gsi, ncols,
                                                              gamma, p, q);

    compute_prim_gsk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, fsk0, fsi,
                                                              fsk1, gsh0, gsh1, gsi, ncols,
                                                              gamma, p, q);

    compute_prim_gsk_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, fsk0, fsi,
                                                              fsk1, gsh0, gsh1, gsi, ncols,
                                                              gamma, p, q);

    compute_prim_gsk_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, fsk0, fsi,
                                                              fsk1, gsh0, gsh1, gsi, ncols,
                                                              gamma, p, q);

    compute_prim_gsk_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, fsk0, fsi,
                                                              fsk1, gsh0, gsh1, gsi, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
