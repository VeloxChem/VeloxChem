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


#include "SimdThreeCenterElectronRepulsionVrrRecHSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_hsk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsk0,
                                                          const size_t gsi, const size_t gsk1,
                                                          const size_t hsh0, const size_t hsh1,
                                                          const size_t hsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_16 = 2.0 / q;
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

    const auto *gsk0_0 = buffer.data(gsk0 + 0);
    const auto *gsk0_3 = buffer.data(gsk0 + 3);
    const auto *gsk0_5 = buffer.data(gsk0 + 5);
    const auto *gsk0_6 = buffer.data(gsk0 + 6);
    const auto *gsk0_9 = buffer.data(gsk0 + 9);
    const auto *gsk0_10 = buffer.data(gsk0 + 10);
    const auto *gsk0_14 = buffer.data(gsk0 + 14);
    const auto *gsk0_15 = buffer.data(gsk0 + 15);
    const auto *gsk0_20 = buffer.data(gsk0 + 20);
    const auto *gsk0_28 = buffer.data(gsk0 + 28);
    const auto *gsk0_35 = buffer.data(gsk0 + 35);

    const auto *gsi_0 = buffer.data(gsi + 0);
    const auto *gsi_1 = buffer.data(gsi + 1);
    const auto *gsi_2 = buffer.data(gsi + 2);
    const auto *gsi_3 = buffer.data(gsi + 3);
    const auto *gsi_5 = buffer.data(gsi + 5);
    const auto *gsi_6 = buffer.data(gsi + 6);
    const auto *gsi_9 = buffer.data(gsi + 9);
    const auto *gsi_10 = buffer.data(gsi + 10);
    const auto *gsi_14 = buffer.data(gsi + 14);
    const auto *gsi_21 = buffer.data(gsi + 21);
    const auto *gsi_23 = buffer.data(gsi + 23);
    const auto *gsi_24 = buffer.data(gsi + 24);
    const auto *gsi_25 = buffer.data(gsi + 25);
    const auto *gsi_27 = buffer.data(gsi + 27);
    const auto *gsi_28 = buffer.data(gsi + 28);
    const auto *gsi_33 = buffer.data(gsi + 33);
    const auto *gsi_37 = buffer.data(gsi + 37);
    const auto *gsi_42 = buffer.data(gsi + 42);
    const auto *gsi_49 = buffer.data(gsi + 49);
    const auto *gsi_51 = buffer.data(gsi + 51);
    const auto *gsi_52 = buffer.data(gsi + 52);
    const auto *gsi_53 = buffer.data(gsi + 53);
    const auto *gsi_54 = buffer.data(gsi + 54);
    const auto *gsi_55 = buffer.data(gsi + 55);
    const auto *gsi_77 = buffer.data(gsi + 77);
    const auto *gsi_78 = buffer.data(gsi + 78);
    const auto *gsi_79 = buffer.data(gsi + 79);
    const auto *gsi_80 = buffer.data(gsi + 80);
    const auto *gsi_81 = buffer.data(gsi + 81);
    const auto *gsi_83 = buffer.data(gsi + 83);
    const auto *gsi_84 = buffer.data(gsi + 84);
    const auto *gsi_87 = buffer.data(gsi + 87);
    const auto *gsi_90 = buffer.data(gsi + 90);
    const auto *gsi_94 = buffer.data(gsi + 94);
    const auto *gsi_99 = buffer.data(gsi + 99);

    const auto *gsk1_0 = buffer.data(gsk1 + 0);
    const auto *gsk1_3 = buffer.data(gsk1 + 3);
    const auto *gsk1_5 = buffer.data(gsk1 + 5);
    const auto *gsk1_6 = buffer.data(gsk1 + 6);
    const auto *gsk1_9 = buffer.data(gsk1 + 9);
    const auto *gsk1_10 = buffer.data(gsk1 + 10);
    const auto *gsk1_14 = buffer.data(gsk1 + 14);
    const auto *gsk1_15 = buffer.data(gsk1 + 15);
    const auto *gsk1_20 = buffer.data(gsk1 + 20);
    const auto *gsk1_28 = buffer.data(gsk1 + 28);
    const auto *gsk1_35 = buffer.data(gsk1 + 35);

    const auto *hsh0_0 = buffer.data(hsh0 + 0);
    const auto *hsh0_1 = buffer.data(hsh0 + 1);
    const auto *hsh0_2 = buffer.data(hsh0 + 2);
    const auto *hsh0_3 = buffer.data(hsh0 + 3);
    const auto *hsh0_5 = buffer.data(hsh0 + 5);
    const auto *hsh0_6 = buffer.data(hsh0 + 6);
    const auto *hsh0_8 = buffer.data(hsh0 + 8);
    const auto *hsh0_9 = buffer.data(hsh0 + 9);
    const auto *hsh0_15 = buffer.data(hsh0 + 15);
    const auto *hsh0_17 = buffer.data(hsh0 + 17);
    const auto *hsh0_18 = buffer.data(hsh0 + 18);
    const auto *hsh0_19 = buffer.data(hsh0 + 19);
    const auto *hsh0_20 = buffer.data(hsh0 + 20);
    const auto *hsh0_24 = buffer.data(hsh0 + 24);
    const auto *hsh0_27 = buffer.data(hsh0 + 27);
    const auto *hsh0_28 = buffer.data(hsh0 + 28);
    const auto *hsh0_36 = buffer.data(hsh0 + 36);
    const auto *hsh0_37 = buffer.data(hsh0 + 37);
    const auto *hsh0_38 = buffer.data(hsh0 + 38);
    const auto *hsh0_39 = buffer.data(hsh0 + 39);
    const auto *hsh0_44 = buffer.data(hsh0 + 44);
    const auto *hsh0_46 = buffer.data(hsh0 + 46);
    const auto *hsh0_47 = buffer.data(hsh0 + 47);
    const auto *hsh0_49 = buffer.data(hsh0 + 49);
    const auto *hsh0_50 = buffer.data(hsh0 + 50);
    const auto *hsh0_51 = buffer.data(hsh0 + 51);
    const auto *hsh0_58 = buffer.data(hsh0 + 58);
    const auto *hsh0_59 = buffer.data(hsh0 + 59);
    const auto *hsh0_60 = buffer.data(hsh0 + 60);
    const auto *hsh0_61 = buffer.data(hsh0 + 61);
    const auto *hsh0_62 = buffer.data(hsh0 + 62);
    const auto *hsh0_63 = buffer.data(hsh0 + 63);
    const auto *hsh0_65 = buffer.data(hsh0 + 65);
    const auto *hsh0_66 = buffer.data(hsh0 + 66);
    const auto *hsh0_68 = buffer.data(hsh0 + 68);
    const auto *hsh0_69 = buffer.data(hsh0 + 69);
    const auto *hsh0_70 = buffer.data(hsh0 + 70);
    const auto *hsh0_72 = buffer.data(hsh0 + 72);
    const auto *hsh0_73 = buffer.data(hsh0 + 73);
    const auto *hsh0_78 = buffer.data(hsh0 + 78);

    const auto *hsh1_0 = buffer.data(hsh1 + 0);
    const auto *hsh1_1 = buffer.data(hsh1 + 1);
    const auto *hsh1_2 = buffer.data(hsh1 + 2);
    const auto *hsh1_3 = buffer.data(hsh1 + 3);
    const auto *hsh1_5 = buffer.data(hsh1 + 5);
    const auto *hsh1_6 = buffer.data(hsh1 + 6);
    const auto *hsh1_8 = buffer.data(hsh1 + 8);
    const auto *hsh1_9 = buffer.data(hsh1 + 9);
    const auto *hsh1_15 = buffer.data(hsh1 + 15);
    const auto *hsh1_17 = buffer.data(hsh1 + 17);
    const auto *hsh1_18 = buffer.data(hsh1 + 18);
    const auto *hsh1_19 = buffer.data(hsh1 + 19);
    const auto *hsh1_20 = buffer.data(hsh1 + 20);
    const auto *hsh1_24 = buffer.data(hsh1 + 24);
    const auto *hsh1_27 = buffer.data(hsh1 + 27);
    const auto *hsh1_28 = buffer.data(hsh1 + 28);
    const auto *hsh1_36 = buffer.data(hsh1 + 36);
    const auto *hsh1_37 = buffer.data(hsh1 + 37);
    const auto *hsh1_38 = buffer.data(hsh1 + 38);
    const auto *hsh1_39 = buffer.data(hsh1 + 39);
    const auto *hsh1_44 = buffer.data(hsh1 + 44);
    const auto *hsh1_46 = buffer.data(hsh1 + 46);
    const auto *hsh1_47 = buffer.data(hsh1 + 47);
    const auto *hsh1_49 = buffer.data(hsh1 + 49);
    const auto *hsh1_50 = buffer.data(hsh1 + 50);
    const auto *hsh1_51 = buffer.data(hsh1 + 51);
    const auto *hsh1_58 = buffer.data(hsh1 + 58);
    const auto *hsh1_59 = buffer.data(hsh1 + 59);
    const auto *hsh1_60 = buffer.data(hsh1 + 60);
    const auto *hsh1_61 = buffer.data(hsh1 + 61);
    const auto *hsh1_62 = buffer.data(hsh1 + 62);
    const auto *hsh1_63 = buffer.data(hsh1 + 63);
    const auto *hsh1_65 = buffer.data(hsh1 + 65);
    const auto *hsh1_66 = buffer.data(hsh1 + 66);
    const auto *hsh1_68 = buffer.data(hsh1 + 68);
    const auto *hsh1_69 = buffer.data(hsh1 + 69);
    const auto *hsh1_70 = buffer.data(hsh1 + 70);
    const auto *hsh1_72 = buffer.data(hsh1 + 72);
    const auto *hsh1_73 = buffer.data(hsh1 + 73);
    const auto *hsh1_78 = buffer.data(hsh1 + 78);

    const auto *hsi_0 = buffer.data(hsi + 0);
    const auto *hsi_1 = buffer.data(hsi + 1);
    const auto *hsi_2 = buffer.data(hsi + 2);
    const auto *hsi_3 = buffer.data(hsi + 3);
    const auto *hsi_5 = buffer.data(hsi + 5);
    const auto *hsi_6 = buffer.data(hsi + 6);
    const auto *hsi_8 = buffer.data(hsi + 8);
    const auto *hsi_9 = buffer.data(hsi + 9);
    const auto *hsi_10 = buffer.data(hsi + 10);
    const auto *hsi_12 = buffer.data(hsi + 12);
    const auto *hsi_13 = buffer.data(hsi + 13);
    const auto *hsi_14 = buffer.data(hsi + 14);
    const auto *hsi_15 = buffer.data(hsi + 15);
    const auto *hsi_20 = buffer.data(hsi + 20);
    const auto *hsi_21 = buffer.data(hsi + 21);
    const auto *hsi_23 = buffer.data(hsi + 23);
    const auto *hsi_24 = buffer.data(hsi + 24);
    const auto *hsi_25 = buffer.data(hsi + 25);
    const auto *hsi_26 = buffer.data(hsi + 26);
    const auto *hsi_27 = buffer.data(hsi + 27);
    const auto *hsi_28 = buffer.data(hsi + 28);
    const auto *hsi_29 = buffer.data(hsi + 29);
    const auto *hsi_31 = buffer.data(hsi + 31);
    const auto *hsi_33 = buffer.data(hsi + 33);
    const auto *hsi_34 = buffer.data(hsi + 34);
    const auto *hsi_35 = buffer.data(hsi + 35);
    const auto *hsi_37 = buffer.data(hsi + 37);
    const auto *hsi_38 = buffer.data(hsi + 38);
    const auto *hsi_39 = buffer.data(hsi + 39);
    const auto *hsi_40 = buffer.data(hsi + 40);
    const auto *hsi_42 = buffer.data(hsi + 42);
    const auto *hsi_43 = buffer.data(hsi + 43);
    const auto *hsi_49 = buffer.data(hsi + 49);
    const auto *hsi_50 = buffer.data(hsi + 50);
    const auto *hsi_51 = buffer.data(hsi + 51);
    const auto *hsi_52 = buffer.data(hsi + 52);
    const auto *hsi_53 = buffer.data(hsi + 53);
    const auto *hsi_54 = buffer.data(hsi + 54);
    const auto *hsi_55 = buffer.data(hsi + 55);
    const auto *hsi_56 = buffer.data(hsi + 56);
    const auto *hsi_58 = buffer.data(hsi + 58);
    const auto *hsi_60 = buffer.data(hsi + 60);
    const auto *hsi_61 = buffer.data(hsi + 61);
    const auto *hsi_63 = buffer.data(hsi + 63);
    const auto *hsi_64 = buffer.data(hsi + 64);
    const auto *hsi_65 = buffer.data(hsi + 65);
    const auto *hsi_67 = buffer.data(hsi + 67);
    const auto *hsi_68 = buffer.data(hsi + 68);
    const auto *hsi_69 = buffer.data(hsi + 69);
    const auto *hsi_70 = buffer.data(hsi + 70);
    const auto *hsi_76 = buffer.data(hsi + 76);
    const auto *hsi_77 = buffer.data(hsi + 77);
    const auto *hsi_78 = buffer.data(hsi + 78);
    const auto *hsi_79 = buffer.data(hsi + 79);
    const auto *hsi_80 = buffer.data(hsi + 80);
    const auto *hsi_81 = buffer.data(hsi + 81);
    const auto *hsi_82 = buffer.data(hsi + 82);
    const auto *hsi_83 = buffer.data(hsi + 83);
    const auto *hsi_84 = buffer.data(hsi + 84);
    const auto *hsi_85 = buffer.data(hsi + 85);
    const auto *hsi_86 = buffer.data(hsi + 86);
    const auto *hsi_87 = buffer.data(hsi + 87);
    const auto *hsi_89 = buffer.data(hsi + 89);
    const auto *hsi_90 = buffer.data(hsi + 90);
    const auto *hsi_91 = buffer.data(hsi + 91);
    const auto *hsi_93 = buffer.data(hsi + 93);
    const auto *hsi_94 = buffer.data(hsi + 94);
    const auto *hsi_95 = buffer.data(hsi + 95);
    const auto *hsi_96 = buffer.data(hsi + 96);
    const auto *hsi_98 = buffer.data(hsi + 98);
    const auto *hsi_99 = buffer.data(hsi + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, gsi_0, hsh0_0, \
                         hsh1_0, hsi_0, hsi_1, hsi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gsi_0[k]
                 + f_1 * hsh0_0[k]
                 - f_2 * hsh1_0[k]
                 + f_3 * pc_x[k] * hsi_0[k];

        t_1[k] = f_3 * pc_y[k] * hsi_0[k];

        t_2[k] = f_3 * pc_z[k] * hsi_0[k];

        t_3[k] = f_4 * hsh0_0[k]
                 - f_5 * hsh1_0[k]
                 + f_3 * pc_y[k] * hsi_1[k];

        t_4[k] = f_3 * pc_y[k] * hsi_2[k];

        t_5[k] = f_4 * hsh0_0[k]
                 - f_5 * hsh1_0[k]
                 + f_3 * pc_z[k] * hsi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, hsh0_1, hsh0_2, hsh0_3, hsh1_1, \
                         hsh1_2, hsh1_3, hsi_3, hsi_5, hsi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hsh0_1[k]
                 - f_7 * hsh1_1[k]
                 + f_3 * pc_y[k] * hsi_3[k];

        t_7[k] = f_3 * pc_z[k] * hsi_3[k];

        t_8[k] = f_3 * pc_y[k] * hsi_5[k];

        t_9[k] = f_6 * hsh0_2[k]
                 - f_7 * hsh1_2[k]
                 + f_3 * pc_z[k] * hsi_5[k];

        t_10[k] = f_8 * hsh0_3[k]
                  - f_9 * hsh1_3[k]
                  + f_3 * pc_y[k] * hsi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, hsh0_5, hsh0_6, \
                         hsh1_5, hsh1_6, hsi_6, hsi_8, hsi_9, hsi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * hsi_6[k];

        t_12[k] = f_4 * hsh0_5[k]
                  - f_5 * hsh1_5[k]
                  + f_3 * pc_y[k] * hsi_8[k];

        t_13[k] = f_3 * pc_y[k] * hsi_9[k];

        t_14[k] = f_8 * hsh0_5[k]
                  - f_9 * hsh1_5[k]
                  + f_3 * pc_z[k] * hsi_9[k];

        t_15[k] = f_10 * hsh0_6[k]
                  - f_11 * hsh1_6[k]
                  + f_3 * pc_y[k] * hsi_10[k];

        t_16[k] = f_3 * pc_z[k] * hsi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, hsh0_8, hsh0_9, hsh1_8, hsh1_9, \
                         hsi_12, hsi_13, hsi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * hsh0_8[k]
                  - f_7 * hsh1_8[k]
                  + f_3 * pc_y[k] * hsi_12[k];

        t_18[k] = f_4 * hsh0_9[k]
                  - f_5 * hsh1_9[k]
                  + f_3 * pc_y[k] * hsi_13[k];

        t_19[k] = f_3 * pc_y[k] * hsi_14[k];

        t_20[k] = f_10 * hsh0_9[k]
                  - f_11 * hsh1_9[k]
                  + f_3 * pc_z[k] * hsi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, gsi_21, gsi_23, gsi_24, \
                         gsi_25, hsi_15, hsi_21, hsi_23, hsi_24, \
                         hsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * gsi_21[k]
                  + f_3 * pc_x[k] * hsi_21[k];

        t_22[k] = f_3 * pc_z[k] * hsi_15[k];

        t_23[k] = f_0 * gsi_23[k]
                  + f_3 * pc_x[k] * hsi_23[k];

        t_24[k] = f_0 * gsi_24[k]
                  + f_3 * pc_x[k] * hsi_24[k];

        t_25[k] = f_0 * gsi_25[k]
                  + f_3 * pc_x[k] * hsi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, gsi_27, hsh0_15, hsh1_15, \
                         hsi_20, hsi_21, hsi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * hsi_20[k];

        t_27[k] = f_0 * gsi_27[k]
                  + f_3 * pc_x[k] * hsi_27[k];

        t_28[k] = f_1 * hsh0_15[k]
                  - f_2 * hsh1_15[k]
                  + f_3 * pc_y[k] * hsi_21[k];

        t_29[k] = f_3 * pc_z[k] * hsi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, hsh0_17, hsh0_18, hsh0_19, hsh1_17, hsh1_18, \
                         hsh1_19, hsi_23, hsi_24, hsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * hsh0_17[k]
                  - f_11 * hsh1_17[k]
                  + f_3 * pc_y[k] * hsi_23[k];

        t_31[k] = f_8 * hsh0_18[k]
                  - f_9 * hsh1_18[k]
                  + f_3 * pc_y[k] * hsi_24[k];

        t_32[k] = f_6 * hsh0_19[k]
                  - f_7 * hsh1_19[k]
                  + f_3 * pc_y[k] * hsi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, gsk0_0, gsi_0, \
                         gsk1_0, hsh0_20, hsh1_20, hsi_26, hsi_27, \
                         hsi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * hsh0_20[k]
                  - f_5 * hsh1_20[k]
                  + f_3 * pc_y[k] * hsi_26[k];

        t_34[k] = f_3 * pc_y[k] * hsi_27[k];

        t_35[k] = f_1 * hsh0_20[k]
                  - f_2 * hsh1_20[k]
                  + f_3 * pc_z[k] * hsi_27[k];

        t_36[k] = pa_y[k] * gsk0_0[k]
                  - f_12 * pc_y[k] * gsk1_0[k];

        t_37[k] = f_13 * gsi_0[k]
                  + f_3 * pc_y[k] * hsi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, gsk0_3, gsk0_5, gsi_1, \
                         gsk1_3, gsk1_5, hsi_28, hsi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * hsi_28[k];

        t_39[k] = pa_y[k] * gsk0_3[k]
                  + f_14 * gsi_1[k]
                  - f_12 * pc_y[k] * gsk1_3[k];

        t_40[k] = f_3 * pc_z[k] * hsi_29[k];

        t_41[k] = pa_y[k] * gsk0_5[k]
                  - f_12 * pc_y[k] * gsk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_y, pc_z, gsk0_6, gsk0_9, gsi_3, \
                         gsi_5, gsk1_6, gsk1_9, hsi_31, hsi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * gsk0_6[k]
                  + f_15 * gsi_3[k]
                  - f_12 * pc_y[k] * gsk1_6[k];

        t_43[k] = f_3 * pc_z[k] * hsi_31[k];

        t_44[k] = f_13 * gsi_5[k]
                  + f_3 * pc_y[k] * hsi_33[k];

        t_45[k] = pa_y[k] * gsk0_9[k]
                  - f_12 * pc_y[k] * gsk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, gsk0_10, gsi_6, gsi_9, \
                         gsk1_10, hsh0_24, hsh1_24, hsi_34, hsi_35, \
                         hsi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * gsk0_10[k]
                  + f_16 * gsi_6[k]
                  - f_12 * pc_y[k] * gsk1_10[k];

        t_47[k] = f_3 * pc_z[k] * hsi_34[k];

        t_48[k] = f_4 * hsh0_24[k]
                  - f_5 * hsh1_24[k]
                  + f_3 * pc_z[k] * hsi_35[k];

        t_49[k] = f_13 * gsi_9[k]
                  + f_3 * pc_y[k] * hsi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_y, pc_z, gsk0_14, gsk0_15, gsi_10, \
                         gsk1_14, gsk1_15, hsh0_27, hsh1_27, hsi_38, \
                         hsi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * gsk0_14[k]
                  - f_12 * pc_y[k] * gsk1_14[k];

        t_51[k] = pa_y[k] * gsk0_15[k]
                  + f_0 * gsi_10[k]
                  - f_12 * pc_y[k] * gsk1_15[k];

        t_52[k] = f_3 * pc_z[k] * hsi_38[k];

        t_53[k] = f_4 * hsh0_27[k]
                  - f_5 * hsh1_27[k]
                  + f_3 * pc_z[k] * hsi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, gsk0_20, gsi_14, gsk1_20, \
                         hsh0_28, hsh1_28, hsi_40, hsi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * hsh0_28[k]
                  - f_7 * hsh1_28[k]
                  + f_3 * pc_z[k] * hsi_40[k];

        t_55[k] = f_13 * gsi_14[k]
                  + f_3 * pc_y[k] * hsi_42[k];

        t_56[k] = pa_y[k] * gsk0_20[k]
                  - f_12 * pc_y[k] * gsk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, gsi_49, gsi_51, gsi_52, \
                         gsi_53, hsi_43, hsi_49, hsi_51, hsi_52, \
                         hsi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_16 * gsi_49[k]
                  + f_3 * pc_x[k] * hsi_49[k];

        t_58[k] = f_3 * pc_z[k] * hsi_43[k];

        t_59[k] = f_16 * gsi_51[k]
                  + f_3 * pc_x[k] * hsi_51[k];

        t_60[k] = f_16 * gsi_52[k]
                  + f_3 * pc_x[k] * hsi_52[k];

        t_61[k] = f_16 * gsi_53[k]
                  + f_3 * pc_x[k] * hsi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, gsi_21, gsi_54, gsi_55, \
                         hsh0_36, hsh1_36, hsi_49, hsi_54, hsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_16 * gsi_54[k]
                  + f_3 * pc_x[k] * hsi_54[k];

        t_63[k] = f_16 * gsi_55[k]
                  + f_3 * pc_x[k] * hsi_55[k];

        t_64[k] = f_13 * gsi_21[k]
                  + f_1 * hsh0_36[k]
                  - f_2 * hsh1_36[k]
                  + f_3 * pc_y[k] * hsi_49[k];

        t_65[k] = f_3 * pc_z[k] * hsi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_z, hsh0_36, hsh0_37, hsh0_38, hsh1_36, hsh1_37, \
                         hsh1_38, hsi_50, hsi_51, hsi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * hsh0_36[k]
                  - f_5 * hsh1_36[k]
                  + f_3 * pc_z[k] * hsi_50[k];

        t_67[k] = f_6 * hsh0_37[k]
                  - f_7 * hsh1_37[k]
                  + f_3 * pc_z[k] * hsi_51[k];

        t_68[k] = f_8 * hsh0_38[k]
                  - f_9 * hsh1_38[k]
                  + f_3 * pc_z[k] * hsi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, gsk0_35, gsi_27, gsk1_35, \
                         hsh0_39, hsh1_39, hsi_53, hsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * hsh0_39[k]
                  - f_11 * hsh1_39[k]
                  + f_3 * pc_z[k] * hsi_53[k];

        t_70[k] = f_13 * gsi_27[k]
                  + f_3 * pc_y[k] * hsi_55[k];

        t_71[k] = pa_y[k] * gsk0_35[k]
                  - f_12 * pc_y[k] * gsk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pc_y, pc_z, gsk0_0, gsk0_3, \
                         gsi_0, gsk1_0, gsk1_3, hsi_56, hsi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * gsk0_0[k]
                  - f_12 * pc_z[k] * gsk1_0[k];

        t_73[k] = f_3 * pc_y[k] * hsi_56[k];

        t_74[k] = f_13 * gsi_0[k]
                  + f_3 * pc_z[k] * hsi_56[k];

        t_75[k] = pa_z[k] * gsk0_3[k]
                  - f_12 * pc_z[k] * gsk1_3[k];

        t_76[k] = f_3 * pc_y[k] * hsi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pc_y, pc_z, gsk0_5, gsk0_6, gsi_2, \
                         gsk1_5, gsk1_6, hsh0_44, hsh1_44, hsi_60, \
                         hsi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * gsk0_5[k]
                  + f_14 * gsi_2[k]
                  - f_12 * pc_z[k] * gsk1_5[k];

        t_78[k] = pa_z[k] * gsk0_6[k]
                  - f_12 * pc_z[k] * gsk1_6[k];

        t_79[k] = f_4 * hsh0_44[k]
                  - f_5 * hsh1_44[k]
                  + f_3 * pc_y[k] * hsi_60[k];

        t_80[k] = f_3 * pc_y[k] * hsi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, gsk0_9, gsk0_10, gsi_5, gsk1_9, \
                         gsk1_10, hsh0_46, hsh1_46, hsi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * gsk0_9[k]
                  + f_15 * gsi_5[k]
                  - f_12 * pc_z[k] * gsk1_9[k];

        t_82[k] = pa_z[k] * gsk0_10[k]
                  - f_12 * pc_z[k] * gsk1_10[k];

        t_83[k] = f_6 * hsh0_46[k]
                  - f_7 * hsh1_46[k]
                  + f_3 * pc_y[k] * hsi_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pc_y, pc_z, gsk0_14, gsk0_15, gsi_9, \
                         gsk1_14, gsk1_15, hsh0_47, hsh1_47, hsi_64, \
                         hsi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * hsh0_47[k]
                  - f_5 * hsh1_47[k]
                  + f_3 * pc_y[k] * hsi_64[k];

        t_85[k] = f_3 * pc_y[k] * hsi_65[k];

        t_86[k] = pa_z[k] * gsk0_14[k]
                  + f_16 * gsi_9[k]
                  - f_12 * pc_z[k] * gsk1_14[k];

        t_87[k] = pa_z[k] * gsk0_15[k]
                  - f_12 * pc_z[k] * gsk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, hsh0_49, hsh0_50, hsh0_51, hsh1_49, \
                         hsh1_50, hsh1_51, hsi_67, hsi_68, hsi_69, \
                         hsi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_8 * hsh0_49[k]
                  - f_9 * hsh1_49[k]
                  + f_3 * pc_y[k] * hsi_67[k];

        t_89[k] = f_6 * hsh0_50[k]
                  - f_7 * hsh1_50[k]
                  + f_3 * pc_y[k] * hsi_68[k];

        t_90[k] = f_4 * hsh0_51[k]
                  - f_5 * hsh1_51[k]
                  + f_3 * pc_y[k] * hsi_69[k];

        t_91[k] = f_3 * pc_y[k] * hsi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_z, gsk0_20, gsi_14, gsi_77, \
                         gsi_78, gsi_79, gsk1_20, hsi_77, hsi_78, \
                         hsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * gsk0_20[k]
                  + f_0 * gsi_14[k]
                  - f_12 * pc_z[k] * gsk1_20[k];

        t_93[k] = f_16 * gsi_77[k]
                  + f_3 * pc_x[k] * hsi_77[k];

        t_94[k] = f_16 * gsi_78[k]
                  + f_3 * pc_x[k] * hsi_78[k];

        t_95[k] = f_16 * gsi_79[k]
                  + f_3 * pc_x[k] * hsi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, gsi_80, gsi_81, gsi_83, hsi_76, \
                         hsi_80, hsi_81, hsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_16 * gsi_80[k]
                  + f_3 * pc_x[k] * hsi_80[k];

        t_97[k] = f_16 * gsi_81[k]
                  + f_3 * pc_x[k] * hsi_81[k];

        t_98[k] = f_3 * pc_y[k] * hsi_76[k];

        t_99[k] = f_16 * gsi_83[k]
                  + f_3 * pc_x[k] * hsi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_y, pc_z, gsk0_28, gsk1_28, hsh0_58, \
                         hsh0_59, hsh1_58, hsh1_59, hsi_78, hsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * gsk0_28[k]
                   - f_12 * pc_z[k] * gsk1_28[k];

        t_101[k] = f_17 * hsh0_58[k]
                   - f_18 * hsh1_58[k]
                   + f_3 * pc_y[k] * hsi_78[k];

        t_102[k] = f_10 * hsh0_59[k]
                   - f_11 * hsh1_59[k]
                   + f_3 * pc_y[k] * hsi_79[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_y, hsh0_60, hsh0_61, hsh0_62, hsh1_60, \
                         hsh1_61, hsh1_62, hsi_80, hsi_81, hsi_82, \
                         hsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * hsh0_60[k]
                   - f_9 * hsh1_60[k]
                   + f_3 * pc_y[k] * hsi_80[k];

        t_104[k] = f_6 * hsh0_61[k]
                   - f_7 * hsh1_61[k]
                   + f_3 * pc_y[k] * hsi_81[k];

        t_105[k] = f_4 * hsh0_62[k]
                   - f_5 * hsh1_62[k]
                   + f_3 * pc_y[k] * hsi_82[k];

        t_106[k] = f_3 * pc_y[k] * hsi_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, pc_z, gsi_27, gsi_28, gsi_84, \
                         hsh0_62, hsh0_63, hsh1_62, hsh1_63, hsi_83, \
                         hsi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_13 * gsi_27[k]
                   + f_1 * hsh0_62[k]
                   - f_2 * hsh1_62[k]
                   + f_3 * pc_z[k] * hsi_83[k];

        t_108[k] = f_15 * gsi_84[k]
                   + f_1 * hsh0_63[k]
                   - f_2 * hsh1_63[k]
                   + f_3 * pc_x[k] * hsi_84[k];

        t_109[k] = f_14 * gsi_28[k]
                   + f_3 * pc_y[k] * hsi_84[k];

        t_110[k] = f_3 * pc_z[k] * hsi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_z, gsi_87, hsh0_63, hsh0_66, hsh1_63, \
                         hsh1_66, hsi_85, hsi_86, hsi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_15 * gsi_87[k]
                   + f_10 * hsh0_66[k]
                   - f_11 * hsh1_66[k]
                   + f_3 * pc_x[k] * hsi_87[k];

        t_112[k] = f_3 * pc_z[k] * hsi_85[k];

        t_113[k] = f_4 * hsh0_63[k]
                   - f_5 * hsh1_63[k]
                   + f_3 * pc_z[k] * hsi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, pc_y, pc_z, gsi_33, gsi_90, \
                         hsh0_65, hsh0_69, hsh1_65, hsh1_69, hsi_87, hsi_89, \
                         hsi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * gsi_90[k]
                   + f_8 * hsh0_69[k]
                   - f_9 * hsh1_69[k]
                   + f_3 * pc_x[k] * hsi_90[k];

        t_115[k] = f_3 * pc_z[k] * hsi_87[k];

        t_116[k] = f_14 * gsi_33[k]
                   + f_3 * pc_y[k] * hsi_89[k];

        t_117[k] = f_6 * hsh0_65[k]
                   - f_7 * hsh1_65[k]
                   + f_3 * pc_z[k] * hsi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_z, gsi_94, hsh0_66, hsh0_73, hsh1_66, \
                         hsh1_73, hsi_90, hsi_91, hsi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_15 * gsi_94[k]
                   + f_6 * hsh0_73[k]
                   - f_7 * hsh1_73[k]
                   + f_3 * pc_x[k] * hsi_94[k];

        t_119[k] = f_3 * pc_z[k] * hsi_90[k];

        t_120[k] = f_4 * hsh0_66[k]
                   - f_5 * hsh1_66[k]
                   + f_3 * pc_z[k] * hsi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, gsi_37, gsi_99, \
                         hsh0_68, hsh0_78, hsh1_68, hsh1_78, hsi_93, hsi_94, \
                         hsi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_14 * gsi_37[k]
                   + f_3 * pc_y[k] * hsi_93[k];

        t_122[k] = f_8 * hsh0_68[k]
                   - f_9 * hsh1_68[k]
                   + f_3 * pc_z[k] * hsi_93[k];

        t_123[k] = f_15 * gsi_99[k]
                   + f_4 * hsh0_78[k]
                   - f_5 * hsh1_78[k]
                   + f_3 * pc_x[k] * hsi_99[k];

        t_124[k] = f_3 * pc_z[k] * hsi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_y, pc_z, gsi_42, hsh0_69, hsh0_70, \
                         hsh0_72, hsh1_69, hsh1_70, hsh1_72, hsi_95, hsi_96, \
                         hsi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * hsh0_69[k]
                   - f_5 * hsh1_69[k]
                   + f_3 * pc_z[k] * hsi_95[k];

        t_126[k] = f_6 * hsh0_70[k]
                   - f_7 * hsh1_70[k]
                   + f_3 * pc_z[k] * hsi_96[k];

        t_127[k] = f_14 * gsi_42[k]
                   + f_3 * pc_y[k] * hsi_98[k];

        t_128[k] = f_10 * hsh0_72[k]
                   - f_11 * hsh1_72[k]
                   + f_3 * pc_z[k] * hsi_98[k];
    }
}

static auto
compute_prim_hsk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsk0,
                                                          const size_t gsi, const size_t gsk1,
                                                          const size_t hsh0, const size_t hsh1,
                                                          const size_t hsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsk0_39 = buffer.data(gsk0 + 39);
    const auto *gsk0_42 = buffer.data(gsk0 + 42);
    const auto *gsk0_46 = buffer.data(gsk0 + 46);
    const auto *gsk0_51 = buffer.data(gsk0 + 51);
    const auto *gsk0_64 = buffer.data(gsk0 + 64);
    const auto *gsk0_72 = buffer.data(gsk0 + 72);
    const auto *gsk0_77 = buffer.data(gsk0 + 77);
    const auto *gsk0_81 = buffer.data(gsk0 + 81);
    const auto *gsk0_84 = buffer.data(gsk0 + 84);
    const auto *gsk0_86 = buffer.data(gsk0 + 86);
    const auto *gsk0_89 = buffer.data(gsk0 + 89);
    const auto *gsk0_90 = buffer.data(gsk0 + 90);
    const auto *gsk0_92 = buffer.data(gsk0 + 92);
    const auto *gsk0_107 = buffer.data(gsk0 + 107);

    const auto *gsi_28 = buffer.data(gsi + 28);
    const auto *gsi_31 = buffer.data(gsi + 31);
    const auto *gsi_34 = buffer.data(gsi + 34);
    const auto *gsi_38 = buffer.data(gsi + 38);
    const auto *gsi_49 = buffer.data(gsi + 49);
    const auto *gsi_55 = buffer.data(gsi + 55);
    const auto *gsi_56 = buffer.data(gsi + 56);
    const auto *gsi_58 = buffer.data(gsi + 58);
    const auto *gsi_61 = buffer.data(gsi + 61);
    const auto *gsi_64 = buffer.data(gsi + 64);
    const auto *gsi_65 = buffer.data(gsi + 65);
    const auto *gsi_68 = buffer.data(gsi + 68);
    const auto *gsi_69 = buffer.data(gsi + 69);
    const auto *gsi_70 = buffer.data(gsi + 70);
    const auto *gsi_79 = buffer.data(gsi + 79);
    const auto *gsi_80 = buffer.data(gsi + 80);
    const auto *gsi_81 = buffer.data(gsi + 81);
    const auto *gsi_82 = buffer.data(gsi + 82);
    const auto *gsi_83 = buffer.data(gsi + 83);
    const auto *gsi_84 = buffer.data(gsi + 84);
    const auto *gsi_89 = buffer.data(gsi + 89);
    const auto *gsi_93 = buffer.data(gsi + 93);
    const auto *gsi_98 = buffer.data(gsi + 98);
    const auto *gsi_105 = buffer.data(gsi + 105);
    const auto *gsi_107 = buffer.data(gsi + 107);
    const auto *gsi_108 = buffer.data(gsi + 108);
    const auto *gsi_109 = buffer.data(gsi + 109);
    const auto *gsi_110 = buffer.data(gsi + 110);
    const auto *gsi_111 = buffer.data(gsi + 111);
    const auto *gsi_133 = buffer.data(gsi + 133);
    const auto *gsi_134 = buffer.data(gsi + 134);
    const auto *gsi_135 = buffer.data(gsi + 135);
    const auto *gsi_136 = buffer.data(gsi + 136);
    const auto *gsi_137 = buffer.data(gsi + 137);
    const auto *gsi_138 = buffer.data(gsi + 138);
    const auto *gsi_139 = buffer.data(gsi + 139);
    const auto *gsi_140 = buffer.data(gsi + 140);
    const auto *gsi_145 = buffer.data(gsi + 145);
    const auto *gsi_149 = buffer.data(gsi + 149);
    const auto *gsi_154 = buffer.data(gsi + 154);
    const auto *gsi_160 = buffer.data(gsi + 160);
    const auto *gsi_161 = buffer.data(gsi + 161);
    const auto *gsi_162 = buffer.data(gsi + 162);
    const auto *gsi_163 = buffer.data(gsi + 163);
    const auto *gsi_164 = buffer.data(gsi + 164);
    const auto *gsi_165 = buffer.data(gsi + 165);
    const auto *gsi_167 = buffer.data(gsi + 167);
    const auto *gsi_168 = buffer.data(gsi + 168);
    const auto *gsi_171 = buffer.data(gsi + 171);
    const auto *gsi_174 = buffer.data(gsi + 174);
    const auto *gsi_178 = buffer.data(gsi + 178);
    const auto *gsi_183 = buffer.data(gsi + 183);
    const auto *gsi_189 = buffer.data(gsi + 189);
    const auto *gsi_191 = buffer.data(gsi + 191);
    const auto *gsi_192 = buffer.data(gsi + 192);
    const auto *gsi_193 = buffer.data(gsi + 193);
    const auto *gsi_194 = buffer.data(gsi + 194);
    const auto *gsi_195 = buffer.data(gsi + 195);

    const auto *gsk1_39 = buffer.data(gsk1 + 39);
    const auto *gsk1_42 = buffer.data(gsk1 + 42);
    const auto *gsk1_46 = buffer.data(gsk1 + 46);
    const auto *gsk1_51 = buffer.data(gsk1 + 51);
    const auto *gsk1_64 = buffer.data(gsk1 + 64);
    const auto *gsk1_72 = buffer.data(gsk1 + 72);
    const auto *gsk1_77 = buffer.data(gsk1 + 77);
    const auto *gsk1_81 = buffer.data(gsk1 + 81);
    const auto *gsk1_84 = buffer.data(gsk1 + 84);
    const auto *gsk1_86 = buffer.data(gsk1 + 86);
    const auto *gsk1_89 = buffer.data(gsk1 + 89);
    const auto *gsk1_90 = buffer.data(gsk1 + 90);
    const auto *gsk1_92 = buffer.data(gsk1 + 92);
    const auto *gsk1_107 = buffer.data(gsk1 + 107);

    const auto *hsh0_78 = buffer.data(hsh0 + 78);
    const auto *hsh0_79 = buffer.data(hsh0 + 79);
    const auto *hsh0_80 = buffer.data(hsh0 + 80);
    const auto *hsh0_81 = buffer.data(hsh0 + 81);
    const auto *hsh0_83 = buffer.data(hsh0 + 83);
    const auto *hsh0_101 = buffer.data(hsh0 + 101);
    const auto *hsh0_102 = buffer.data(hsh0 + 102);
    const auto *hsh0_103 = buffer.data(hsh0 + 103);
    const auto *hsh0_104 = buffer.data(hsh0 + 104);
    const auto *hsh0_105 = buffer.data(hsh0 + 105);
    const auto *hsh0_106 = buffer.data(hsh0 + 106);
    const auto *hsh0_107 = buffer.data(hsh0 + 107);
    const auto *hsh0_108 = buffer.data(hsh0 + 108);
    const auto *hsh0_109 = buffer.data(hsh0 + 109);
    const auto *hsh0_110 = buffer.data(hsh0 + 110);
    const auto *hsh0_111 = buffer.data(hsh0 + 111);
    const auto *hsh0_112 = buffer.data(hsh0 + 112);
    const auto *hsh0_113 = buffer.data(hsh0 + 113);
    const auto *hsh0_114 = buffer.data(hsh0 + 114);
    const auto *hsh0_119 = buffer.data(hsh0 + 119);
    const auto *hsh0_120 = buffer.data(hsh0 + 120);
    const auto *hsh0_121 = buffer.data(hsh0 + 121);
    const auto *hsh0_122 = buffer.data(hsh0 + 122);
    const auto *hsh0_123 = buffer.data(hsh0 + 123);
    const auto *hsh0_124 = buffer.data(hsh0 + 124);
    const auto *hsh0_125 = buffer.data(hsh0 + 125);
    const auto *hsh0_126 = buffer.data(hsh0 + 126);
    const auto *hsh0_128 = buffer.data(hsh0 + 128);
    const auto *hsh0_129 = buffer.data(hsh0 + 129);
    const auto *hsh0_131 = buffer.data(hsh0 + 131);
    const auto *hsh0_132 = buffer.data(hsh0 + 132);
    const auto *hsh0_133 = buffer.data(hsh0 + 133);
    const auto *hsh0_135 = buffer.data(hsh0 + 135);
    const auto *hsh0_136 = buffer.data(hsh0 + 136);
    const auto *hsh0_141 = buffer.data(hsh0 + 141);
    const auto *hsh0_142 = buffer.data(hsh0 + 142);
    const auto *hsh0_143 = buffer.data(hsh0 + 143);
    const auto *hsh0_144 = buffer.data(hsh0 + 144);

    const auto *hsh1_78 = buffer.data(hsh1 + 78);
    const auto *hsh1_79 = buffer.data(hsh1 + 79);
    const auto *hsh1_80 = buffer.data(hsh1 + 80);
    const auto *hsh1_81 = buffer.data(hsh1 + 81);
    const auto *hsh1_83 = buffer.data(hsh1 + 83);
    const auto *hsh1_101 = buffer.data(hsh1 + 101);
    const auto *hsh1_102 = buffer.data(hsh1 + 102);
    const auto *hsh1_103 = buffer.data(hsh1 + 103);
    const auto *hsh1_104 = buffer.data(hsh1 + 104);
    const auto *hsh1_105 = buffer.data(hsh1 + 105);
    const auto *hsh1_106 = buffer.data(hsh1 + 106);
    const auto *hsh1_107 = buffer.data(hsh1 + 107);
    const auto *hsh1_108 = buffer.data(hsh1 + 108);
    const auto *hsh1_109 = buffer.data(hsh1 + 109);
    const auto *hsh1_110 = buffer.data(hsh1 + 110);
    const auto *hsh1_111 = buffer.data(hsh1 + 111);
    const auto *hsh1_112 = buffer.data(hsh1 + 112);
    const auto *hsh1_113 = buffer.data(hsh1 + 113);
    const auto *hsh1_114 = buffer.data(hsh1 + 114);
    const auto *hsh1_119 = buffer.data(hsh1 + 119);
    const auto *hsh1_120 = buffer.data(hsh1 + 120);
    const auto *hsh1_121 = buffer.data(hsh1 + 121);
    const auto *hsh1_122 = buffer.data(hsh1 + 122);
    const auto *hsh1_123 = buffer.data(hsh1 + 123);
    const auto *hsh1_124 = buffer.data(hsh1 + 124);
    const auto *hsh1_125 = buffer.data(hsh1 + 125);
    const auto *hsh1_126 = buffer.data(hsh1 + 126);
    const auto *hsh1_128 = buffer.data(hsh1 + 128);
    const auto *hsh1_129 = buffer.data(hsh1 + 129);
    const auto *hsh1_131 = buffer.data(hsh1 + 131);
    const auto *hsh1_132 = buffer.data(hsh1 + 132);
    const auto *hsh1_133 = buffer.data(hsh1 + 133);
    const auto *hsh1_135 = buffer.data(hsh1 + 135);
    const auto *hsh1_136 = buffer.data(hsh1 + 136);
    const auto *hsh1_141 = buffer.data(hsh1 + 141);
    const auto *hsh1_142 = buffer.data(hsh1 + 142);
    const auto *hsh1_143 = buffer.data(hsh1 + 143);
    const auto *hsh1_144 = buffer.data(hsh1 + 144);

    const auto *hsi_99 = buffer.data(hsi + 99);
    const auto *hsi_105 = buffer.data(hsi + 105);
    const auto *hsi_106 = buffer.data(hsi + 106);
    const auto *hsi_107 = buffer.data(hsi + 107);
    const auto *hsi_108 = buffer.data(hsi + 108);
    const auto *hsi_109 = buffer.data(hsi + 109);
    const auto *hsi_110 = buffer.data(hsi + 110);
    const auto *hsi_111 = buffer.data(hsi + 111);
    const auto *hsi_112 = buffer.data(hsi + 112);
    const auto *hsi_114 = buffer.data(hsi + 114);
    const auto *hsi_115 = buffer.data(hsi + 115);
    const auto *hsi_117 = buffer.data(hsi + 117);
    const auto *hsi_118 = buffer.data(hsi + 118);
    const auto *hsi_121 = buffer.data(hsi + 121);
    const auto *hsi_122 = buffer.data(hsi + 122);
    const auto *hsi_126 = buffer.data(hsi + 126);
    const auto *hsi_133 = buffer.data(hsi + 133);
    const auto *hsi_134 = buffer.data(hsi + 134);
    const auto *hsi_135 = buffer.data(hsi + 135);
    const auto *hsi_136 = buffer.data(hsi + 136);
    const auto *hsi_137 = buffer.data(hsi + 137);
    const auto *hsi_138 = buffer.data(hsi + 138);
    const auto *hsi_139 = buffer.data(hsi + 139);
    const auto *hsi_140 = buffer.data(hsi + 140);
    const auto *hsi_141 = buffer.data(hsi + 141);
    const auto *hsi_142 = buffer.data(hsi + 142);
    const auto *hsi_143 = buffer.data(hsi + 143);
    const auto *hsi_144 = buffer.data(hsi + 144);
    const auto *hsi_145 = buffer.data(hsi + 145);
    const auto *hsi_146 = buffer.data(hsi + 146);
    const auto *hsi_147 = buffer.data(hsi + 147);
    const auto *hsi_148 = buffer.data(hsi + 148);
    const auto *hsi_149 = buffer.data(hsi + 149);
    const auto *hsi_150 = buffer.data(hsi + 150);
    const auto *hsi_151 = buffer.data(hsi + 151);
    const auto *hsi_152 = buffer.data(hsi + 152);
    const auto *hsi_153 = buffer.data(hsi + 153);
    const auto *hsi_154 = buffer.data(hsi + 154);
    const auto *hsi_160 = buffer.data(hsi + 160);
    const auto *hsi_161 = buffer.data(hsi + 161);
    const auto *hsi_162 = buffer.data(hsi + 162);
    const auto *hsi_163 = buffer.data(hsi + 163);
    const auto *hsi_164 = buffer.data(hsi + 164);
    const auto *hsi_165 = buffer.data(hsi + 165);
    const auto *hsi_166 = buffer.data(hsi + 166);
    const auto *hsi_167 = buffer.data(hsi + 167);
    const auto *hsi_168 = buffer.data(hsi + 168);
    const auto *hsi_169 = buffer.data(hsi + 169);
    const auto *hsi_170 = buffer.data(hsi + 170);
    const auto *hsi_171 = buffer.data(hsi + 171);
    const auto *hsi_173 = buffer.data(hsi + 173);
    const auto *hsi_174 = buffer.data(hsi + 174);
    const auto *hsi_175 = buffer.data(hsi + 175);
    const auto *hsi_177 = buffer.data(hsi + 177);
    const auto *hsi_178 = buffer.data(hsi + 178);
    const auto *hsi_179 = buffer.data(hsi + 179);
    const auto *hsi_180 = buffer.data(hsi + 180);
    const auto *hsi_182 = buffer.data(hsi + 182);
    const auto *hsi_183 = buffer.data(hsi + 183);
    const auto *hsi_189 = buffer.data(hsi + 189);
    const auto *hsi_190 = buffer.data(hsi + 190);
    const auto *hsi_191 = buffer.data(hsi + 191);
    const auto *hsi_192 = buffer.data(hsi + 192);
    const auto *hsi_193 = buffer.data(hsi + 193);
    const auto *hsi_194 = buffer.data(hsi + 194);
    const auto *hsi_195 = buffer.data(hsi + 195);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pc_x, pc_z, gsi_105, gsi_107, \
                         gsi_108, gsi_109, hsi_99, hsi_105, hsi_107, hsi_108, \
                         hsi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_15 * gsi_105[k]
                   + f_3 * pc_x[k] * hsi_105[k];

        t_130[k] = f_3 * pc_z[k] * hsi_99[k];

        t_131[k] = f_15 * gsi_107[k]
                   + f_3 * pc_x[k] * hsi_107[k];

        t_132[k] = f_15 * gsi_108[k]
                   + f_3 * pc_x[k] * hsi_108[k];

        t_133[k] = f_15 * gsi_109[k]
                   + f_3 * pc_x[k] * hsi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, gsi_49, gsi_110, \
                         gsi_111, hsh0_78, hsh1_78, hsi_105, hsi_110, \
                         hsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * gsi_110[k]
                   + f_3 * pc_x[k] * hsi_110[k];

        t_135[k] = f_15 * gsi_111[k]
                   + f_3 * pc_x[k] * hsi_111[k];

        t_136[k] = f_14 * gsi_49[k]
                   + f_1 * hsh0_78[k]
                   - f_2 * hsh1_78[k]
                   + f_3 * pc_y[k] * hsi_105[k];

        t_137[k] = f_3 * pc_z[k] * hsi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_z, hsh0_78, hsh0_79, hsh0_80, hsh1_78, \
                         hsh1_79, hsh1_80, hsi_106, hsi_107, hsi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_4 * hsh0_78[k]
                   - f_5 * hsh1_78[k]
                   + f_3 * pc_z[k] * hsi_106[k];

        t_139[k] = f_6 * hsh0_79[k]
                   - f_7 * hsh1_79[k]
                   + f_3 * pc_z[k] * hsi_107[k];

        t_140[k] = f_8 * hsh0_80[k]
                   - f_9 * hsh1_80[k]
                   + f_3 * pc_z[k] * hsi_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pc_y, pc_z, gsk0_72, gsi_55, \
                         gsk1_72, hsh0_81, hsh0_83, hsh1_81, hsh1_83, hsi_109, \
                         hsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_10 * hsh0_81[k]
                   - f_11 * hsh1_81[k]
                   + f_3 * pc_z[k] * hsi_109[k];

        t_142[k] = f_14 * gsi_55[k]
                   + f_3 * pc_y[k] * hsi_111[k];

        t_143[k] = f_1 * hsh0_83[k]
                   - f_2 * hsh1_83[k]
                   + f_3 * pc_z[k] * hsi_111[k];

        t_144[k] = pa_y[k] * gsk0_72[k]
                   - f_12 * pc_y[k] * gsk1_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_z, pc_y, pc_z, gsk0_39, gsi_28, \
                         gsi_56, gsi_58, gsk1_39, hsi_112, hsi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * gsi_56[k]
                   + f_3 * pc_y[k] * hsi_112[k];

        t_146[k] = f_13 * gsi_28[k]
                   + f_3 * pc_z[k] * hsi_112[k];

        t_147[k] = pa_z[k] * gsk0_39[k]
                   - f_12 * pc_z[k] * gsk1_39[k];

        t_148[k] = f_13 * gsi_58[k]
                   + f_3 * pc_y[k] * hsi_114[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pa_z, pc_y, pc_z, gsk0_42, gsk0_77, \
                         gsi_31, gsi_61, gsk1_42, gsk1_77, hsi_115, \
                         hsi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * gsk0_77[k]
                   - f_12 * pc_y[k] * gsk1_77[k];

        t_150[k] = pa_z[k] * gsk0_42[k]
                   - f_12 * pc_z[k] * gsk1_42[k];

        t_151[k] = f_13 * gsi_31[k]
                   + f_3 * pc_z[k] * hsi_115[k];

        t_152[k] = f_13 * gsi_61[k]
                   + f_3 * pc_y[k] * hsi_117[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_y, pa_z, pc_y, pc_z, gsk0_46, gsk0_81, \
                         gsi_34, gsk1_46, gsk1_81, hsi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_y[k] * gsk0_81[k]
                   - f_12 * pc_y[k] * gsk1_81[k];

        t_154[k] = pa_z[k] * gsk0_46[k]
                   - f_12 * pc_z[k] * gsk1_46[k];

        t_155[k] = f_13 * gsi_34[k]
                   + f_3 * pc_z[k] * hsi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_y, pc_y, gsk0_84, gsk0_86, gsi_64, gsi_65, \
                         gsk1_84, gsk1_86, hsi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_y[k] * gsk0_84[k]
                   + f_14 * gsi_64[k]
                   - f_12 * pc_y[k] * gsk1_84[k];

        t_157[k] = f_13 * gsi_65[k]
                   + f_3 * pc_y[k] * hsi_121[k];

        t_158[k] = pa_y[k] * gsk0_86[k]
                   - f_12 * pc_y[k] * gsk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pa_z, pc_y, pc_z, gsk0_51, gsk0_89, \
                         gsi_38, gsi_68, gsk1_51, gsk1_89, hsi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * gsk0_51[k]
                   - f_12 * pc_z[k] * gsk1_51[k];

        t_160[k] = f_13 * gsi_38[k]
                   + f_3 * pc_z[k] * hsi_122[k];

        t_161[k] = pa_y[k] * gsk0_89[k]
                   + f_15 * gsi_68[k]
                   - f_12 * pc_y[k] * gsk1_89[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_y, pc_x, pc_y, gsk0_90, gsk0_92, \
                         gsi_69, gsi_70, gsi_133, gsk1_90, gsk1_92, hsi_126, \
                         hsi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_y[k] * gsk0_90[k]
                   + f_14 * gsi_69[k]
                   - f_12 * pc_y[k] * gsk1_90[k];

        t_163[k] = f_13 * gsi_70[k]
                   + f_3 * pc_y[k] * hsi_126[k];

        t_164[k] = pa_y[k] * gsk0_92[k]
                   - f_12 * pc_y[k] * gsk1_92[k];

        t_165[k] = f_15 * gsi_133[k]
                   + f_3 * pc_x[k] * hsi_133[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, gsi_134, gsi_135, gsi_136, \
                         gsi_137, gsi_138, hsi_134, hsi_135, hsi_136, hsi_137, \
                         hsi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_15 * gsi_134[k]
                   + f_3 * pc_x[k] * hsi_134[k];

        t_167[k] = f_15 * gsi_135[k]
                   + f_3 * pc_x[k] * hsi_135[k];

        t_168[k] = f_15 * gsi_136[k]
                   + f_3 * pc_x[k] * hsi_136[k];

        t_169[k] = f_15 * gsi_137[k]
                   + f_3 * pc_x[k] * hsi_137[k];

        t_170[k] = f_15 * gsi_138[k]
                   + f_3 * pc_x[k] * hsi_138[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_z, pc_x, pc_z, gsk0_64, gsi_49, gsi_139, \
                         gsk1_64, hsi_133, hsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_15 * gsi_139[k]
                   + f_3 * pc_x[k] * hsi_139[k];

        t_172[k] = pa_z[k] * gsk0_64[k]
                   - f_12 * pc_z[k] * gsk1_64[k];

        t_173[k] = f_13 * gsi_49[k]
                   + f_3 * pc_z[k] * hsi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, gsi_79, gsi_80, gsi_81, hsh0_101, \
                         hsh0_102, hsh0_103, hsh1_101, hsh1_102, hsh1_103, hsi_135, hsi_136, \
                         hsi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * gsi_79[k]
                   + f_10 * hsh0_101[k]
                   - f_11 * hsh1_101[k]
                   + f_3 * pc_y[k] * hsi_135[k];

        t_175[k] = f_13 * gsi_80[k]
                   + f_8 * hsh0_102[k]
                   - f_9 * hsh1_102[k]
                   + f_3 * pc_y[k] * hsi_136[k];

        t_176[k] = f_13 * gsi_81[k]
                   + f_6 * hsh0_103[k]
                   - f_7 * hsh1_103[k]
                   + f_3 * pc_y[k] * hsi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, gsk0_107, gsi_82, gsi_83, gsk1_107, \
                         hsh0_104, hsh1_104, hsi_138, hsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * gsi_82[k]
                   + f_4 * hsh0_104[k]
                   - f_5 * hsh1_104[k]
                   + f_3 * pc_y[k] * hsi_138[k];

        t_178[k] = f_13 * gsi_83[k]
                   + f_3 * pc_y[k] * hsi_139[k];

        t_179[k] = pa_y[k] * gsk0_107[k]
                   - f_12 * pc_y[k] * gsk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, gsi_56, gsi_140, \
                         hsh0_105, hsh1_105, hsi_140, hsi_141, \
                         hsi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_15 * gsi_140[k]
                   + f_1 * hsh0_105[k]
                   - f_2 * hsh1_105[k]
                   + f_3 * pc_x[k] * hsi_140[k];

        t_181[k] = f_3 * pc_y[k] * hsi_140[k];

        t_182[k] = f_14 * gsi_56[k]
                   + f_3 * pc_z[k] * hsi_140[k];

        t_183[k] = f_4 * hsh0_105[k]
                   - f_5 * hsh1_105[k]
                   + f_3 * pc_y[k] * hsi_141[k];

        t_184[k] = f_3 * pc_y[k] * hsi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, gsi_145, hsh0_106, hsh0_107, \
                         hsh0_110, hsh1_106, hsh1_107, hsh1_110, hsi_143, hsi_144, \
                         hsi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * gsi_145[k]
                   + f_10 * hsh0_110[k]
                   - f_11 * hsh1_110[k]
                   + f_3 * pc_x[k] * hsi_145[k];

        t_186[k] = f_6 * hsh0_106[k]
                   - f_7 * hsh1_106[k]
                   + f_3 * pc_y[k] * hsi_143[k];

        t_187[k] = f_4 * hsh0_107[k]
                   - f_5 * hsh1_107[k]
                   + f_3 * pc_y[k] * hsi_144[k];

        t_188[k] = f_3 * pc_y[k] * hsi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, gsi_149, hsh0_108, hsh0_109, \
                         hsh0_114, hsh1_108, hsh1_109, hsh1_114, hsi_146, hsi_147, \
                         hsi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_15 * gsi_149[k]
                   + f_8 * hsh0_114[k]
                   - f_9 * hsh1_114[k]
                   + f_3 * pc_x[k] * hsi_149[k];

        t_190[k] = f_8 * hsh0_108[k]
                   - f_9 * hsh1_108[k]
                   + f_3 * pc_y[k] * hsi_146[k];

        t_191[k] = f_6 * hsh0_109[k]
                   - f_7 * hsh1_109[k]
                   + f_3 * pc_y[k] * hsi_147[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pc_x, pc_y, gsi_154, hsh0_110, hsh0_119, \
                         hsh1_110, hsh1_119, hsi_148, hsi_149, \
                         hsi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_4 * hsh0_110[k]
                   - f_5 * hsh1_110[k]
                   + f_3 * pc_y[k] * hsi_148[k];

        t_193[k] = f_3 * pc_y[k] * hsi_149[k];

        t_194[k] = f_15 * gsi_154[k]
                   + f_6 * hsh0_119[k]
                   - f_7 * hsh1_119[k]
                   + f_3 * pc_x[k] * hsi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, hsh0_111, hsh0_112, hsh0_113, hsh1_111, \
                         hsh1_112, hsh1_113, hsi_150, hsi_151, \
                         hsi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * hsh0_111[k]
                   - f_11 * hsh1_111[k]
                   + f_3 * pc_y[k] * hsi_150[k];

        t_196[k] = f_8 * hsh0_112[k]
                   - f_9 * hsh1_112[k]
                   + f_3 * pc_y[k] * hsi_151[k];

        t_197[k] = f_6 * hsh0_113[k]
                   - f_7 * hsh1_113[k]
                   + f_3 * pc_y[k] * hsi_152[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, gsi_160, gsi_161, hsh0_114, \
                         hsh0_125, hsh1_114, hsh1_125, hsi_153, hsi_154, hsi_160, \
                         hsi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * hsh0_114[k]
                   - f_5 * hsh1_114[k]
                   + f_3 * pc_y[k] * hsi_153[k];

        t_199[k] = f_3 * pc_y[k] * hsi_154[k];

        t_200[k] = f_15 * gsi_160[k]
                   + f_4 * hsh0_125[k]
                   - f_5 * hsh1_125[k]
                   + f_3 * pc_x[k] * hsi_160[k];

        t_201[k] = f_15 * gsi_161[k]
                   + f_3 * pc_x[k] * hsi_161[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pc_x, pc_y, gsi_162, gsi_163, \
                         gsi_164, gsi_165, hsi_160, hsi_162, hsi_163, hsi_164, \
                         hsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_15 * gsi_162[k]
                   + f_3 * pc_x[k] * hsi_162[k];

        t_203[k] = f_15 * gsi_163[k]
                   + f_3 * pc_x[k] * hsi_163[k];

        t_204[k] = f_15 * gsi_164[k]
                   + f_3 * pc_x[k] * hsi_164[k];

        t_205[k] = f_15 * gsi_165[k]
                   + f_3 * pc_x[k] * hsi_165[k];

        t_206[k] = f_3 * pc_y[k] * hsi_160[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pc_x, pc_y, gsi_167, hsh0_120, hsh0_121, \
                         hsh1_120, hsh1_121, hsi_161, hsi_162, \
                         hsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_15 * gsi_167[k]
                   + f_3 * pc_x[k] * hsi_167[k];

        t_208[k] = f_1 * hsh0_120[k]
                   - f_2 * hsh1_120[k]
                   + f_3 * pc_y[k] * hsi_161[k];

        t_209[k] = f_17 * hsh0_121[k]
                   - f_18 * hsh1_121[k]
                   + f_3 * pc_y[k] * hsi_162[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pc_y, hsh0_122, hsh0_123, hsh0_124, hsh1_122, \
                         hsh1_123, hsh1_124, hsi_163, hsi_164, \
                         hsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * hsh0_122[k]
                   - f_11 * hsh1_122[k]
                   + f_3 * pc_y[k] * hsi_163[k];

        t_211[k] = f_8 * hsh0_123[k]
                   - f_9 * hsh1_123[k]
                   + f_3 * pc_y[k] * hsi_164[k];

        t_212[k] = f_6 * hsh0_124[k]
                   - f_7 * hsh1_124[k]
                   + f_3 * pc_y[k] * hsi_165[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, pc_z, gsi_83, gsi_168, \
                         hsh0_125, hsh0_126, hsh1_125, hsh1_126, hsi_166, hsi_167, \
                         hsi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * hsh0_125[k]
                   - f_5 * hsh1_125[k]
                   + f_3 * pc_y[k] * hsi_166[k];

        t_214[k] = f_3 * pc_y[k] * hsi_167[k];

        t_215[k] = f_14 * gsi_83[k]
                   + f_1 * hsh0_125[k]
                   - f_2 * hsh1_125[k]
                   + f_3 * pc_z[k] * hsi_167[k];

        t_216[k] = f_14 * gsi_168[k]
                   + f_1 * hsh0_126[k]
                   - f_2 * hsh1_126[k]
                   + f_3 * pc_x[k] * hsi_168[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, gsi_84, gsi_171, \
                         hsh0_129, hsh1_129, hsi_168, hsi_169, \
                         hsi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * gsi_84[k]
                   + f_3 * pc_y[k] * hsi_168[k];

        t_218[k] = f_3 * pc_z[k] * hsi_168[k];

        t_219[k] = f_14 * gsi_171[k]
                   + f_10 * hsh0_129[k]
                   - f_11 * hsh1_129[k]
                   + f_3 * pc_x[k] * hsi_171[k];

        t_220[k] = f_3 * pc_z[k] * hsi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pc_x, pc_z, gsi_174, hsh0_126, hsh0_132, \
                         hsh1_126, hsh1_132, hsi_170, hsi_171, \
                         hsi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_4 * hsh0_126[k]
                   - f_5 * hsh1_126[k]
                   + f_3 * pc_z[k] * hsi_170[k];

        t_222[k] = f_14 * gsi_174[k]
                   + f_8 * hsh0_132[k]
                   - f_9 * hsh1_132[k]
                   + f_3 * pc_x[k] * hsi_174[k];

        t_223[k] = f_3 * pc_z[k] * hsi_171[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pc_x, pc_y, pc_z, gsi_89, gsi_178, \
                         hsh0_128, hsh0_136, hsh1_128, hsh1_136, hsi_173, hsi_174, \
                         hsi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_15 * gsi_89[k]
                   + f_3 * pc_y[k] * hsi_173[k];

        t_225[k] = f_6 * hsh0_128[k]
                   - f_7 * hsh1_128[k]
                   + f_3 * pc_z[k] * hsi_173[k];

        t_226[k] = f_14 * gsi_178[k]
                   + f_6 * hsh0_136[k]
                   - f_7 * hsh1_136[k]
                   + f_3 * pc_x[k] * hsi_178[k];

        t_227[k] = f_3 * pc_z[k] * hsi_174[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_y, pc_z, gsi_93, hsh0_129, hsh0_131, \
                         hsh1_129, hsh1_131, hsi_175, hsi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * hsh0_129[k]
                   - f_5 * hsh1_129[k]
                   + f_3 * pc_z[k] * hsi_175[k];

        t_229[k] = f_15 * gsi_93[k]
                   + f_3 * pc_y[k] * hsi_177[k];

        t_230[k] = f_8 * hsh0_131[k]
                   - f_9 * hsh1_131[k]
                   + f_3 * pc_z[k] * hsi_177[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pc_x, pc_z, gsi_183, hsh0_132, hsh0_141, \
                         hsh1_132, hsh1_141, hsi_178, hsi_179, \
                         hsi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_14 * gsi_183[k]
                   + f_4 * hsh0_141[k]
                   - f_5 * hsh1_141[k]
                   + f_3 * pc_x[k] * hsi_183[k];

        t_232[k] = f_3 * pc_z[k] * hsi_178[k];

        t_233[k] = f_4 * hsh0_132[k]
                   - f_5 * hsh1_132[k]
                   + f_3 * pc_z[k] * hsi_179[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, gsi_98, gsi_189, \
                         hsh0_133, hsh0_135, hsh1_133, hsh1_135, hsi_180, hsi_182, \
                         hsi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_6 * hsh0_133[k]
                   - f_7 * hsh1_133[k]
                   + f_3 * pc_z[k] * hsi_180[k];

        t_235[k] = f_15 * gsi_98[k]
                   + f_3 * pc_y[k] * hsi_182[k];

        t_236[k] = f_10 * hsh0_135[k]
                   - f_11 * hsh1_135[k]
                   + f_3 * pc_z[k] * hsi_182[k];

        t_237[k] = f_14 * gsi_189[k]
                   + f_3 * pc_x[k] * hsi_189[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pc_x, pc_z, gsi_191, gsi_192, \
                         gsi_193, gsi_194, hsi_183, hsi_191, hsi_192, hsi_193, \
                         hsi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_3 * pc_z[k] * hsi_183[k];

        t_239[k] = f_14 * gsi_191[k]
                   + f_3 * pc_x[k] * hsi_191[k];

        t_240[k] = f_14 * gsi_192[k]
                   + f_3 * pc_x[k] * hsi_192[k];

        t_241[k] = f_14 * gsi_193[k]
                   + f_3 * pc_x[k] * hsi_193[k];

        t_242[k] = f_14 * gsi_194[k]
                   + f_3 * pc_x[k] * hsi_194[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pc_x, pc_y, pc_z, gsi_105, gsi_195, \
                         hsh0_141, hsh1_141, hsi_189, hsi_190, \
                         hsi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_14 * gsi_195[k]
                   + f_3 * pc_x[k] * hsi_195[k];

        t_244[k] = f_15 * gsi_105[k]
                   + f_1 * hsh0_141[k]
                   - f_2 * hsh1_141[k]
                   + f_3 * pc_y[k] * hsi_189[k];

        t_245[k] = f_3 * pc_z[k] * hsi_189[k];

        t_246[k] = f_4 * hsh0_141[k]
                   - f_5 * hsh1_141[k]
                   + f_3 * pc_z[k] * hsi_190[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_z, hsh0_142, hsh0_143, hsh0_144, hsh1_142, \
                         hsh1_143, hsh1_144, hsi_191, hsi_192, \
                         hsi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * hsh0_142[k]
                   - f_7 * hsh1_142[k]
                   + f_3 * pc_z[k] * hsi_191[k];

        t_248[k] = f_8 * hsh0_143[k]
                   - f_9 * hsh1_143[k]
                   + f_3 * pc_z[k] * hsi_192[k];

        t_249[k] = f_10 * hsh0_144[k]
                   - f_11 * hsh1_144[k]
                   + f_3 * pc_z[k] * hsi_193[k];
    }
}

static auto
compute_prim_hsk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsk0,
                                                          const size_t gsi, const size_t gsk1,
                                                          const size_t hsh0, const size_t hsh1,
                                                          const size_t hsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsk0_108 = buffer.data(gsk0 + 108);
    const auto *gsk0_111 = buffer.data(gsk0 + 111);
    const auto *gsk0_114 = buffer.data(gsk0 + 114);
    const auto *gsk0_118 = buffer.data(gsk0 + 118);
    const auto *gsk0_120 = buffer.data(gsk0 + 120);
    const auto *gsk0_123 = buffer.data(gsk0 + 123);
    const auto *gsk0_125 = buffer.data(gsk0 + 125);
    const auto *gsk0_126 = buffer.data(gsk0 + 126);
    const auto *gsk0_136 = buffer.data(gsk0 + 136);
    const auto *gsk0_180 = buffer.data(gsk0 + 180);
    const auto *gsk0_183 = buffer.data(gsk0 + 183);
    const auto *gsk0_185 = buffer.data(gsk0 + 185);
    const auto *gsk0_186 = buffer.data(gsk0 + 186);
    const auto *gsk0_189 = buffer.data(gsk0 + 189);
    const auto *gsk0_190 = buffer.data(gsk0 + 190);
    const auto *gsk0_192 = buffer.data(gsk0 + 192);
    const auto *gsk0_194 = buffer.data(gsk0 + 194);
    const auto *gsk0_195 = buffer.data(gsk0 + 195);
    const auto *gsk0_197 = buffer.data(gsk0 + 197);
    const auto *gsk0_198 = buffer.data(gsk0 + 198);
    const auto *gsk0_200 = buffer.data(gsk0 + 200);
    const auto *gsk0_215 = buffer.data(gsk0 + 215);
    const auto *gsk0_360 = buffer.data(gsk0 + 360);
    const auto *gsk0_363 = buffer.data(gsk0 + 363);
    const auto *gsk0_366 = buffer.data(gsk0 + 366);

    const auto *gsi_84 = buffer.data(gsi + 84);
    const auto *gsi_87 = buffer.data(gsi + 87);
    const auto *gsi_90 = buffer.data(gsi + 90);
    const auto *gsi_91 = buffer.data(gsi + 91);
    const auto *gsi_94 = buffer.data(gsi + 94);
    const auto *gsi_95 = buffer.data(gsi + 95);
    const auto *gsi_96 = buffer.data(gsi + 96);
    const auto *gsi_105 = buffer.data(gsi + 105);
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
    const auto *gsi_135 = buffer.data(gsi + 135);
    const auto *gsi_136 = buffer.data(gsi + 136);
    const auto *gsi_137 = buffer.data(gsi + 137);
    const auto *gsi_138 = buffer.data(gsi + 138);
    const auto *gsi_139 = buffer.data(gsi + 139);
    const auto *gsi_140 = buffer.data(gsi + 140);
    const auto *gsi_141 = buffer.data(gsi + 141);
    const auto *gsi_142 = buffer.data(gsi + 142);
    const auto *gsi_143 = buffer.data(gsi + 143);
    const auto *gsi_145 = buffer.data(gsi + 145);
    const auto *gsi_146 = buffer.data(gsi + 146);
    const auto *gsi_148 = buffer.data(gsi + 148);
    const auto *gsi_149 = buffer.data(gsi + 149);
    const auto *gsi_150 = buffer.data(gsi + 150);
    const auto *gsi_152 = buffer.data(gsi + 152);
    const auto *gsi_153 = buffer.data(gsi + 153);
    const auto *gsi_154 = buffer.data(gsi + 154);
    const auto *gsi_161 = buffer.data(gsi + 161);
    const auto *gsi_163 = buffer.data(gsi + 163);
    const auto *gsi_164 = buffer.data(gsi + 164);
    const auto *gsi_165 = buffer.data(gsi + 165);
    const auto *gsi_166 = buffer.data(gsi + 166);
    const auto *gsi_167 = buffer.data(gsi + 167);
    const auto *gsi_168 = buffer.data(gsi + 168);
    const auto *gsi_201 = buffer.data(gsi + 201);
    const auto *gsi_205 = buffer.data(gsi + 205);
    const auto *gsi_210 = buffer.data(gsi + 210);
    const auto *gsi_216 = buffer.data(gsi + 216);
    const auto *gsi_217 = buffer.data(gsi + 217);
    const auto *gsi_218 = buffer.data(gsi + 218);
    const auto *gsi_219 = buffer.data(gsi + 219);
    const auto *gsi_220 = buffer.data(gsi + 220);
    const auto *gsi_221 = buffer.data(gsi + 221);
    const auto *gsi_222 = buffer.data(gsi + 222);
    const auto *gsi_223 = buffer.data(gsi + 223);
    const auto *gsi_245 = buffer.data(gsi + 245);
    const auto *gsi_246 = buffer.data(gsi + 246);
    const auto *gsi_247 = buffer.data(gsi + 247);
    const auto *gsi_248 = buffer.data(gsi + 248);
    const auto *gsi_249 = buffer.data(gsi + 249);
    const auto *gsi_250 = buffer.data(gsi + 250);
    const auto *gsi_251 = buffer.data(gsi + 251);
    const auto *gsi_252 = buffer.data(gsi + 252);
    const auto *gsi_257 = buffer.data(gsi + 257);
    const auto *gsi_261 = buffer.data(gsi + 261);
    const auto *gsi_266 = buffer.data(gsi + 266);
    const auto *gsi_272 = buffer.data(gsi + 272);
    const auto *gsi_273 = buffer.data(gsi + 273);
    const auto *gsi_274 = buffer.data(gsi + 274);
    const auto *gsi_275 = buffer.data(gsi + 275);
    const auto *gsi_276 = buffer.data(gsi + 276);
    const auto *gsi_277 = buffer.data(gsi + 277);
    const auto *gsi_279 = buffer.data(gsi + 279);
    const auto *gsi_280 = buffer.data(gsi + 280);
    const auto *gsi_283 = buffer.data(gsi + 283);
    const auto *gsi_286 = buffer.data(gsi + 286);

    const auto *gsk1_108 = buffer.data(gsk1 + 108);
    const auto *gsk1_111 = buffer.data(gsk1 + 111);
    const auto *gsk1_114 = buffer.data(gsk1 + 114);
    const auto *gsk1_118 = buffer.data(gsk1 + 118);
    const auto *gsk1_120 = buffer.data(gsk1 + 120);
    const auto *gsk1_123 = buffer.data(gsk1 + 123);
    const auto *gsk1_125 = buffer.data(gsk1 + 125);
    const auto *gsk1_126 = buffer.data(gsk1 + 126);
    const auto *gsk1_136 = buffer.data(gsk1 + 136);
    const auto *gsk1_180 = buffer.data(gsk1 + 180);
    const auto *gsk1_183 = buffer.data(gsk1 + 183);
    const auto *gsk1_185 = buffer.data(gsk1 + 185);
    const auto *gsk1_186 = buffer.data(gsk1 + 186);
    const auto *gsk1_189 = buffer.data(gsk1 + 189);
    const auto *gsk1_190 = buffer.data(gsk1 + 190);
    const auto *gsk1_192 = buffer.data(gsk1 + 192);
    const auto *gsk1_194 = buffer.data(gsk1 + 194);
    const auto *gsk1_195 = buffer.data(gsk1 + 195);
    const auto *gsk1_197 = buffer.data(gsk1 + 197);
    const auto *gsk1_198 = buffer.data(gsk1 + 198);
    const auto *gsk1_200 = buffer.data(gsk1 + 200);
    const auto *gsk1_215 = buffer.data(gsk1 + 215);
    const auto *gsk1_360 = buffer.data(gsk1 + 360);
    const auto *gsk1_363 = buffer.data(gsk1 + 363);
    const auto *gsk1_366 = buffer.data(gsk1 + 366);

    const auto *hsh0_146 = buffer.data(hsh0 + 146);
    const auto *hsh0_152 = buffer.data(hsh0 + 152);
    const auto *hsh0_156 = buffer.data(hsh0 + 156);
    const auto *hsh0_161 = buffer.data(hsh0 + 161);
    const auto *hsh0_164 = buffer.data(hsh0 + 164);
    const auto *hsh0_165 = buffer.data(hsh0 + 165);
    const auto *hsh0_166 = buffer.data(hsh0 + 166);
    const auto *hsh0_167 = buffer.data(hsh0 + 167);
    const auto *hsh0_183 = buffer.data(hsh0 + 183);
    const auto *hsh0_185 = buffer.data(hsh0 + 185);
    const auto *hsh0_186 = buffer.data(hsh0 + 186);
    const auto *hsh0_187 = buffer.data(hsh0 + 187);
    const auto *hsh0_188 = buffer.data(hsh0 + 188);
    const auto *hsh0_189 = buffer.data(hsh0 + 189);
    const auto *hsh0_190 = buffer.data(hsh0 + 190);
    const auto *hsh0_191 = buffer.data(hsh0 + 191);
    const auto *hsh0_192 = buffer.data(hsh0 + 192);
    const auto *hsh0_193 = buffer.data(hsh0 + 193);
    const auto *hsh0_194 = buffer.data(hsh0 + 194);
    const auto *hsh0_195 = buffer.data(hsh0 + 195);
    const auto *hsh0_196 = buffer.data(hsh0 + 196);
    const auto *hsh0_197 = buffer.data(hsh0 + 197);
    const auto *hsh0_198 = buffer.data(hsh0 + 198);
    const auto *hsh0_203 = buffer.data(hsh0 + 203);
    const auto *hsh0_204 = buffer.data(hsh0 + 204);
    const auto *hsh0_205 = buffer.data(hsh0 + 205);
    const auto *hsh0_206 = buffer.data(hsh0 + 206);
    const auto *hsh0_207 = buffer.data(hsh0 + 207);
    const auto *hsh0_208 = buffer.data(hsh0 + 208);
    const auto *hsh0_209 = buffer.data(hsh0 + 209);
    const auto *hsh0_210 = buffer.data(hsh0 + 210);

    const auto *hsh1_146 = buffer.data(hsh1 + 146);
    const auto *hsh1_152 = buffer.data(hsh1 + 152);
    const auto *hsh1_156 = buffer.data(hsh1 + 156);
    const auto *hsh1_161 = buffer.data(hsh1 + 161);
    const auto *hsh1_164 = buffer.data(hsh1 + 164);
    const auto *hsh1_165 = buffer.data(hsh1 + 165);
    const auto *hsh1_166 = buffer.data(hsh1 + 166);
    const auto *hsh1_167 = buffer.data(hsh1 + 167);
    const auto *hsh1_183 = buffer.data(hsh1 + 183);
    const auto *hsh1_185 = buffer.data(hsh1 + 185);
    const auto *hsh1_186 = buffer.data(hsh1 + 186);
    const auto *hsh1_187 = buffer.data(hsh1 + 187);
    const auto *hsh1_188 = buffer.data(hsh1 + 188);
    const auto *hsh1_189 = buffer.data(hsh1 + 189);
    const auto *hsh1_190 = buffer.data(hsh1 + 190);
    const auto *hsh1_191 = buffer.data(hsh1 + 191);
    const auto *hsh1_192 = buffer.data(hsh1 + 192);
    const auto *hsh1_193 = buffer.data(hsh1 + 193);
    const auto *hsh1_194 = buffer.data(hsh1 + 194);
    const auto *hsh1_195 = buffer.data(hsh1 + 195);
    const auto *hsh1_196 = buffer.data(hsh1 + 196);
    const auto *hsh1_197 = buffer.data(hsh1 + 197);
    const auto *hsh1_198 = buffer.data(hsh1 + 198);
    const auto *hsh1_203 = buffer.data(hsh1 + 203);
    const auto *hsh1_204 = buffer.data(hsh1 + 204);
    const auto *hsh1_205 = buffer.data(hsh1 + 205);
    const auto *hsh1_206 = buffer.data(hsh1 + 206);
    const auto *hsh1_207 = buffer.data(hsh1 + 207);
    const auto *hsh1_208 = buffer.data(hsh1 + 208);
    const auto *hsh1_209 = buffer.data(hsh1 + 209);
    const auto *hsh1_210 = buffer.data(hsh1 + 210);

    const auto *hsi_195 = buffer.data(hsi + 195);
    const auto *hsi_196 = buffer.data(hsi + 196);
    const auto *hsi_198 = buffer.data(hsi + 198);
    const auto *hsi_199 = buffer.data(hsi + 199);
    const auto *hsi_201 = buffer.data(hsi + 201);
    const auto *hsi_202 = buffer.data(hsi + 202);
    const auto *hsi_205 = buffer.data(hsi + 205);
    const auto *hsi_206 = buffer.data(hsi + 206);
    const auto *hsi_210 = buffer.data(hsi + 210);
    const auto *hsi_216 = buffer.data(hsi + 216);
    const auto *hsi_217 = buffer.data(hsi + 217);
    const auto *hsi_218 = buffer.data(hsi + 218);
    const auto *hsi_219 = buffer.data(hsi + 219);
    const auto *hsi_220 = buffer.data(hsi + 220);
    const auto *hsi_221 = buffer.data(hsi + 221);
    const auto *hsi_222 = buffer.data(hsi + 222);
    const auto *hsi_223 = buffer.data(hsi + 223);
    const auto *hsi_224 = buffer.data(hsi + 224);
    const auto *hsi_226 = buffer.data(hsi + 226);
    const auto *hsi_227 = buffer.data(hsi + 227);
    const auto *hsi_229 = buffer.data(hsi + 229);
    const auto *hsi_230 = buffer.data(hsi + 230);
    const auto *hsi_233 = buffer.data(hsi + 233);
    const auto *hsi_234 = buffer.data(hsi + 234);
    const auto *hsi_238 = buffer.data(hsi + 238);
    const auto *hsi_245 = buffer.data(hsi + 245);
    const auto *hsi_246 = buffer.data(hsi + 246);
    const auto *hsi_247 = buffer.data(hsi + 247);
    const auto *hsi_248 = buffer.data(hsi + 248);
    const auto *hsi_249 = buffer.data(hsi + 249);
    const auto *hsi_250 = buffer.data(hsi + 250);
    const auto *hsi_251 = buffer.data(hsi + 251);
    const auto *hsi_252 = buffer.data(hsi + 252);
    const auto *hsi_253 = buffer.data(hsi + 253);
    const auto *hsi_254 = buffer.data(hsi + 254);
    const auto *hsi_255 = buffer.data(hsi + 255);
    const auto *hsi_256 = buffer.data(hsi + 256);
    const auto *hsi_257 = buffer.data(hsi + 257);
    const auto *hsi_258 = buffer.data(hsi + 258);
    const auto *hsi_259 = buffer.data(hsi + 259);
    const auto *hsi_260 = buffer.data(hsi + 260);
    const auto *hsi_261 = buffer.data(hsi + 261);
    const auto *hsi_262 = buffer.data(hsi + 262);
    const auto *hsi_263 = buffer.data(hsi + 263);
    const auto *hsi_264 = buffer.data(hsi + 264);
    const auto *hsi_265 = buffer.data(hsi + 265);
    const auto *hsi_266 = buffer.data(hsi + 266);
    const auto *hsi_272 = buffer.data(hsi + 272);
    const auto *hsi_273 = buffer.data(hsi + 273);
    const auto *hsi_274 = buffer.data(hsi + 274);
    const auto *hsi_275 = buffer.data(hsi + 275);
    const auto *hsi_276 = buffer.data(hsi + 276);
    const auto *hsi_277 = buffer.data(hsi + 277);
    const auto *hsi_278 = buffer.data(hsi + 278);
    const auto *hsi_279 = buffer.data(hsi + 279);
    const auto *hsi_280 = buffer.data(hsi + 280);
    const auto *hsi_281 = buffer.data(hsi + 281);
    const auto *hsi_282 = buffer.data(hsi + 282);
    const auto *hsi_283 = buffer.data(hsi + 283);

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_z, pc_y, pc_z, gsk0_108, gsi_111, \
                         gsi_112, gsk1_108, hsh0_146, hsh1_146, hsi_195, \
                         hsi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * gsi_111[k]
                   + f_3 * pc_y[k] * hsi_195[k];

        t_251[k] = f_1 * hsh0_146[k]
                   - f_2 * hsh1_146[k]
                   + f_3 * pc_z[k] * hsi_195[k];

        t_252[k] = pa_z[k] * gsk0_108[k]
                   - f_12 * pc_z[k] * gsk1_108[k];

        t_253[k] = f_14 * gsi_112[k]
                   + f_3 * pc_y[k] * hsi_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pa_z, pc_y, pc_z, gsk0_111, gsi_84, gsi_114, \
                         gsk1_111, hsi_196, hsi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * gsi_84[k]
                   + f_3 * pc_z[k] * hsi_196[k];

        t_255[k] = pa_z[k] * gsk0_111[k]
                   - f_12 * pc_z[k] * gsk1_111[k];

        t_256[k] = f_14 * gsi_114[k]
                   + f_3 * pc_y[k] * hsi_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_z, pc_x, pc_z, gsk0_114, gsi_87, gsi_201, \
                         gsk1_114, hsh0_152, hsh1_152, hsi_199, \
                         hsi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_14 * gsi_201[k]
                   + f_10 * hsh0_152[k]
                   - f_11 * hsh1_152[k]
                   + f_3 * pc_x[k] * hsi_201[k];

        t_258[k] = pa_z[k] * gsk0_114[k]
                   - f_12 * pc_z[k] * gsk1_114[k];

        t_259[k] = f_13 * gsi_87[k]
                   + f_3 * pc_z[k] * hsi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_z, pc_x, pc_y, pc_z, gsk0_118, gsi_117, \
                         gsi_205, gsk1_118, hsh0_156, hsh1_156, hsi_201, \
                         hsi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * gsi_117[k]
                   + f_3 * pc_y[k] * hsi_201[k];

        t_261[k] = f_14 * gsi_205[k]
                   + f_8 * hsh0_156[k]
                   - f_9 * hsh1_156[k]
                   + f_3 * pc_x[k] * hsi_205[k];

        t_262[k] = pa_z[k] * gsk0_118[k]
                   - f_12 * pc_z[k] * gsk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pa_z, pc_y, pc_z, gsk0_120, gsi_90, gsi_91, \
                         gsi_121, gsk1_120, hsi_202, hsi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * gsi_90[k]
                   + f_3 * pc_z[k] * hsi_202[k];

        t_264[k] = pa_z[k] * gsk0_120[k]
                   + f_14 * gsi_91[k]
                   - f_12 * pc_z[k] * gsk1_120[k];

        t_265[k] = f_14 * gsi_121[k]
                   + f_3 * pc_y[k] * hsi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_z, pc_x, pc_z, gsk0_123, gsi_94, gsi_210, \
                         gsk1_123, hsh0_161, hsh1_161, hsi_206, \
                         hsi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_14 * gsi_210[k]
                   + f_6 * hsh0_161[k]
                   - f_7 * hsh1_161[k]
                   + f_3 * pc_x[k] * hsi_210[k];

        t_267[k] = pa_z[k] * gsk0_123[k]
                   - f_12 * pc_z[k] * gsk1_123[k];

        t_268[k] = f_13 * gsi_94[k]
                   + f_3 * pc_z[k] * hsi_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_z, pc_y, pc_z, gsk0_125, gsk0_126, gsi_95, \
                         gsi_96, gsi_126, gsk1_125, gsk1_126, hsi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_z[k] * gsk0_125[k]
                   + f_14 * gsi_95[k]
                   - f_12 * pc_z[k] * gsk1_125[k];

        t_270[k] = pa_z[k] * gsk0_126[k]
                   + f_15 * gsi_96[k]
                   - f_12 * pc_z[k] * gsk1_126[k];

        t_271[k] = f_14 * gsi_126[k]
                   + f_3 * pc_y[k] * hsi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, gsi_216, gsi_217, gsi_218, gsi_219, \
                         hsh0_167, hsh1_167, hsi_216, hsi_217, hsi_218, \
                         hsi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_14 * gsi_216[k]
                   + f_4 * hsh0_167[k]
                   - f_5 * hsh1_167[k]
                   + f_3 * pc_x[k] * hsi_216[k];

        t_273[k] = f_14 * gsi_217[k]
                   + f_3 * pc_x[k] * hsi_217[k];

        t_274[k] = f_14 * gsi_218[k]
                   + f_3 * pc_x[k] * hsi_218[k];

        t_275[k] = f_14 * gsi_219[k]
                   + f_3 * pc_x[k] * hsi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, gsi_220, gsi_221, gsi_222, gsi_223, \
                         hsi_220, hsi_221, hsi_222, hsi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_14 * gsi_220[k]
                   + f_3 * pc_x[k] * hsi_220[k];

        t_277[k] = f_14 * gsi_221[k]
                   + f_3 * pc_x[k] * hsi_221[k];

        t_278[k] = f_14 * gsi_222[k]
                   + f_3 * pc_x[k] * hsi_222[k];

        t_279[k] = f_14 * gsi_223[k]
                   + f_3 * pc_x[k] * hsi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_z, pc_y, pc_z, gsk0_136, gsi_105, gsi_135, \
                         gsk1_136, hsh0_164, hsh1_164, hsi_217, \
                         hsi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * gsk0_136[k]
                   - f_12 * pc_z[k] * gsk1_136[k];

        t_281[k] = f_13 * gsi_105[k]
                   + f_3 * pc_z[k] * hsi_217[k];

        t_282[k] = f_14 * gsi_135[k]
                   + f_10 * hsh0_164[k]
                   - f_11 * hsh1_164[k]
                   + f_3 * pc_y[k] * hsi_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, gsi_136, gsi_137, gsi_138, hsh0_165, \
                         hsh0_166, hsh0_167, hsh1_165, hsh1_166, hsh1_167, hsi_220, hsi_221, \
                         hsi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * gsi_136[k]
                   + f_8 * hsh0_165[k]
                   - f_9 * hsh1_165[k]
                   + f_3 * pc_y[k] * hsi_220[k];

        t_284[k] = f_14 * gsi_137[k]
                   + f_6 * hsh0_166[k]
                   - f_7 * hsh1_166[k]
                   + f_3 * pc_y[k] * hsi_221[k];

        t_285[k] = f_14 * gsi_138[k]
                   + f_4 * hsh0_167[k]
                   - f_5 * hsh1_167[k]
                   + f_3 * pc_y[k] * hsi_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pc_y, pc_z, gsk0_180, gsi_111, \
                         gsi_139, gsi_140, gsk1_180, hsh0_167, hsh1_167, hsi_223, \
                         hsi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * gsi_139[k]
                   + f_3 * pc_y[k] * hsi_223[k];

        t_287[k] = f_13 * gsi_111[k]
                   + f_1 * hsh0_167[k]
                   - f_2 * hsh1_167[k]
                   + f_3 * pc_z[k] * hsi_223[k];

        t_288[k] = pa_y[k] * gsk0_180[k]
                   - f_12 * pc_y[k] * gsk1_180[k];

        t_289[k] = f_13 * gsi_140[k]
                   + f_3 * pc_y[k] * hsi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pc_y, pc_z, gsk0_183, gsk0_185, \
                         gsi_112, gsi_141, gsi_142, gsk1_183, gsk1_185, hsi_224, \
                         hsi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * gsi_112[k]
                   + f_3 * pc_z[k] * hsi_224[k];

        t_291[k] = pa_y[k] * gsk0_183[k]
                   + f_14 * gsi_141[k]
                   - f_12 * pc_y[k] * gsk1_183[k];

        t_292[k] = f_13 * gsi_142[k]
                   + f_3 * pc_y[k] * hsi_226[k];

        t_293[k] = pa_y[k] * gsk0_185[k]
                   - f_12 * pc_y[k] * gsk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pc_y, pc_z, gsk0_186, gsk0_189, \
                         gsi_115, gsi_143, gsi_145, gsk1_186, gsk1_189, hsi_227, \
                         hsi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pa_y[k] * gsk0_186[k]
                   + f_15 * gsi_143[k]
                   - f_12 * pc_y[k] * gsk1_186[k];

        t_295[k] = f_14 * gsi_115[k]
                   + f_3 * pc_z[k] * hsi_227[k];

        t_296[k] = f_13 * gsi_145[k]
                   + f_3 * pc_y[k] * hsi_229[k];

        t_297[k] = pa_y[k] * gsk0_189[k]
                   - f_12 * pc_y[k] * gsk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pa_y, pc_y, pc_z, gsk0_190, gsk0_192, gsi_118, \
                         gsi_146, gsi_148, gsk1_190, gsk1_192, \
                         hsi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pa_y[k] * gsk0_190[k]
                   + f_16 * gsi_146[k]
                   - f_12 * pc_y[k] * gsk1_190[k];

        t_299[k] = f_14 * gsi_118[k]
                   + f_3 * pc_z[k] * hsi_230[k];

        t_300[k] = pa_y[k] * gsk0_192[k]
                   + f_14 * gsi_148[k]
                   - f_12 * pc_y[k] * gsk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_y, pc_y, pc_z, gsk0_194, gsk0_195, \
                         gsi_122, gsi_149, gsi_150, gsk1_194, gsk1_195, hsi_233, \
                         hsi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * gsi_149[k]
                   + f_3 * pc_y[k] * hsi_233[k];

        t_302[k] = pa_y[k] * gsk0_194[k]
                   - f_12 * pc_y[k] * gsk1_194[k];

        t_303[k] = pa_y[k] * gsk0_195[k]
                   + f_0 * gsi_150[k]
                   - f_12 * pc_y[k] * gsk1_195[k];

        t_304[k] = f_14 * gsi_122[k]
                   + f_3 * pc_z[k] * hsi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_y, pc_y, gsk0_197, gsk0_198, gsk0_200, \
                         gsi_152, gsi_153, gsi_154, gsk1_197, gsk1_198, gsk1_200, \
                         hsi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_y[k] * gsk0_197[k]
                   + f_15 * gsi_152[k]
                   - f_12 * pc_y[k] * gsk1_197[k];

        t_306[k] = pa_y[k] * gsk0_198[k]
                   + f_14 * gsi_153[k]
                   - f_12 * pc_y[k] * gsk1_198[k];

        t_307[k] = f_13 * gsi_154[k]
                   + f_3 * pc_y[k] * hsi_238[k];

        t_308[k] = pa_y[k] * gsk0_200[k]
                   - f_12 * pc_y[k] * gsk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, gsi_245, gsi_246, gsi_247, \
                         gsi_248, gsi_249, hsi_245, hsi_246, hsi_247, hsi_248, \
                         hsi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_14 * gsi_245[k]
                   + f_3 * pc_x[k] * hsi_245[k];

        t_310[k] = f_14 * gsi_246[k]
                   + f_3 * pc_x[k] * hsi_246[k];

        t_311[k] = f_14 * gsi_247[k]
                   + f_3 * pc_x[k] * hsi_247[k];

        t_312[k] = f_14 * gsi_248[k]
                   + f_3 * pc_x[k] * hsi_248[k];

        t_313[k] = f_14 * gsi_249[k]
                   + f_3 * pc_x[k] * hsi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, gsi_133, gsi_161, \
                         gsi_250, gsi_251, hsh0_183, hsh1_183, hsi_245, hsi_250, \
                         hsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_14 * gsi_250[k]
                   + f_3 * pc_x[k] * hsi_250[k];

        t_315[k] = f_14 * gsi_251[k]
                   + f_3 * pc_x[k] * hsi_251[k];

        t_316[k] = f_13 * gsi_161[k]
                   + f_1 * hsh0_183[k]
                   - f_2 * hsh1_183[k]
                   + f_3 * pc_y[k] * hsi_245[k];

        t_317[k] = f_14 * gsi_133[k]
                   + f_3 * pc_z[k] * hsi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, gsi_163, gsi_164, gsi_165, hsh0_185, \
                         hsh0_186, hsh0_187, hsh1_185, hsh1_186, hsh1_187, hsi_247, hsi_248, \
                         hsi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * gsi_163[k]
                   + f_10 * hsh0_185[k]
                   - f_11 * hsh1_185[k]
                   + f_3 * pc_y[k] * hsi_247[k];

        t_319[k] = f_13 * gsi_164[k]
                   + f_8 * hsh0_186[k]
                   - f_9 * hsh1_186[k]
                   + f_3 * pc_y[k] * hsi_248[k];

        t_320[k] = f_13 * gsi_165[k]
                   + f_6 * hsh0_187[k]
                   - f_7 * hsh1_187[k]
                   + f_3 * pc_y[k] * hsi_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_y, pc_y, gsk0_215, gsi_166, gsi_167, \
                         gsk1_215, hsh0_188, hsh1_188, hsi_250, \
                         hsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * gsi_166[k]
                   + f_4 * hsh0_188[k]
                   - f_5 * hsh1_188[k]
                   + f_3 * pc_y[k] * hsi_250[k];

        t_322[k] = f_13 * gsi_167[k]
                   + f_3 * pc_y[k] * hsi_251[k];

        t_323[k] = pa_y[k] * gsk0_215[k]
                   - f_12 * pc_y[k] * gsk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, gsi_140, \
                         gsi_252, hsh0_189, hsh1_189, hsi_252, hsi_253, \
                         hsi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_14 * gsi_252[k]
                   + f_1 * hsh0_189[k]
                   - f_2 * hsh1_189[k]
                   + f_3 * pc_x[k] * hsi_252[k];

        t_325[k] = f_3 * pc_y[k] * hsi_252[k];

        t_326[k] = f_15 * gsi_140[k]
                   + f_3 * pc_z[k] * hsi_252[k];

        t_327[k] = f_4 * hsh0_189[k]
                   - f_5 * hsh1_189[k]
                   + f_3 * pc_y[k] * hsi_253[k];

        t_328[k] = f_3 * pc_y[k] * hsi_254[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pc_x, pc_y, gsi_257, hsh0_190, hsh0_191, \
                         hsh0_194, hsh1_190, hsh1_191, hsh1_194, hsi_255, hsi_256, \
                         hsi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_14 * gsi_257[k]
                   + f_10 * hsh0_194[k]
                   - f_11 * hsh1_194[k]
                   + f_3 * pc_x[k] * hsi_257[k];

        t_330[k] = f_6 * hsh0_190[k]
                   - f_7 * hsh1_190[k]
                   + f_3 * pc_y[k] * hsi_255[k];

        t_331[k] = f_4 * hsh0_191[k]
                   - f_5 * hsh1_191[k]
                   + f_3 * pc_y[k] * hsi_256[k];

        t_332[k] = f_3 * pc_y[k] * hsi_257[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_x, pc_y, gsi_261, hsh0_192, hsh0_193, \
                         hsh0_198, hsh1_192, hsh1_193, hsh1_198, hsi_258, hsi_259, \
                         hsi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_14 * gsi_261[k]
                   + f_8 * hsh0_198[k]
                   - f_9 * hsh1_198[k]
                   + f_3 * pc_x[k] * hsi_261[k];

        t_334[k] = f_8 * hsh0_192[k]
                   - f_9 * hsh1_192[k]
                   + f_3 * pc_y[k] * hsi_258[k];

        t_335[k] = f_6 * hsh0_193[k]
                   - f_7 * hsh1_193[k]
                   + f_3 * pc_y[k] * hsi_259[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, gsi_266, hsh0_194, hsh0_203, \
                         hsh1_194, hsh1_203, hsi_260, hsi_261, \
                         hsi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_4 * hsh0_194[k]
                   - f_5 * hsh1_194[k]
                   + f_3 * pc_y[k] * hsi_260[k];

        t_337[k] = f_3 * pc_y[k] * hsi_261[k];

        t_338[k] = f_14 * gsi_266[k]
                   + f_6 * hsh0_203[k]
                   - f_7 * hsh1_203[k]
                   + f_3 * pc_x[k] * hsi_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_y, hsh0_195, hsh0_196, hsh0_197, hsh1_195, \
                         hsh1_196, hsh1_197, hsi_262, hsi_263, \
                         hsi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_10 * hsh0_195[k]
                   - f_11 * hsh1_195[k]
                   + f_3 * pc_y[k] * hsi_262[k];

        t_340[k] = f_8 * hsh0_196[k]
                   - f_9 * hsh1_196[k]
                   + f_3 * pc_y[k] * hsi_263[k];

        t_341[k] = f_6 * hsh0_197[k]
                   - f_7 * hsh1_197[k]
                   + f_3 * pc_y[k] * hsi_264[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_y, gsi_272, gsi_273, hsh0_198, \
                         hsh0_209, hsh1_198, hsh1_209, hsi_265, hsi_266, hsi_272, \
                         hsi_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_4 * hsh0_198[k]
                   - f_5 * hsh1_198[k]
                   + f_3 * pc_y[k] * hsi_265[k];

        t_343[k] = f_3 * pc_y[k] * hsi_266[k];

        t_344[k] = f_14 * gsi_272[k]
                   + f_4 * hsh0_209[k]
                   - f_5 * hsh1_209[k]
                   + f_3 * pc_x[k] * hsi_272[k];

        t_345[k] = f_14 * gsi_273[k]
                   + f_3 * pc_x[k] * hsi_273[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, pc_x, pc_y, gsi_274, gsi_275, \
                         gsi_276, gsi_277, hsi_272, hsi_274, hsi_275, hsi_276, \
                         hsi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_14 * gsi_274[k]
                   + f_3 * pc_x[k] * hsi_274[k];

        t_347[k] = f_14 * gsi_275[k]
                   + f_3 * pc_x[k] * hsi_275[k];

        t_348[k] = f_14 * gsi_276[k]
                   + f_3 * pc_x[k] * hsi_276[k];

        t_349[k] = f_14 * gsi_277[k]
                   + f_3 * pc_x[k] * hsi_277[k];

        t_350[k] = f_3 * pc_y[k] * hsi_272[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pc_x, pc_y, gsi_279, hsh0_204, hsh0_205, \
                         hsh1_204, hsh1_205, hsi_273, hsi_274, \
                         hsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_14 * gsi_279[k]
                   + f_3 * pc_x[k] * hsi_279[k];

        t_352[k] = f_1 * hsh0_204[k]
                   - f_2 * hsh1_204[k]
                   + f_3 * pc_y[k] * hsi_273[k];

        t_353[k] = f_17 * hsh0_205[k]
                   - f_18 * hsh1_205[k]
                   + f_3 * pc_y[k] * hsi_274[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_y, hsh0_206, hsh0_207, hsh0_208, hsh1_206, \
                         hsh1_207, hsh1_208, hsi_275, hsi_276, \
                         hsi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_10 * hsh0_206[k]
                   - f_11 * hsh1_206[k]
                   + f_3 * pc_y[k] * hsi_275[k];

        t_355[k] = f_8 * hsh0_207[k]
                   - f_9 * hsh1_207[k]
                   + f_3 * pc_y[k] * hsi_276[k];

        t_356[k] = f_6 * hsh0_208[k]
                   - f_7 * hsh1_208[k]
                   + f_3 * pc_y[k] * hsi_277[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pa_x, pc_x, pc_y, pc_z, gsk0_360, \
                         gsi_167, gsi_280, gsk1_360, hsh0_209, hsh1_209, hsi_278, \
                         hsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_4 * hsh0_209[k]
                   - f_5 * hsh1_209[k]
                   + f_3 * pc_y[k] * hsi_278[k];

        t_358[k] = f_3 * pc_y[k] * hsi_279[k];

        t_359[k] = f_15 * gsi_167[k]
                   + f_1 * hsh0_209[k]
                   - f_2 * hsh1_209[k]
                   + f_3 * pc_z[k] * hsi_279[k];

        t_360[k] = pa_x[k] * gsk0_360[k]
                   + f_19 * gsi_280[k]
                   - f_12 * pc_x[k] * gsk1_360[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pa_x, pc_x, pc_y, pc_z, gsk0_363, \
                         gsi_168, gsi_283, gsk1_363, hsi_280, hsi_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_16 * gsi_168[k]
                   + f_3 * pc_y[k] * hsi_280[k];

        t_362[k] = f_3 * pc_z[k] * hsi_280[k];

        t_363[k] = pa_x[k] * gsk0_363[k]
                   + f_0 * gsi_283[k]
                   - f_12 * pc_x[k] * gsk1_363[k];

        t_364[k] = f_3 * pc_z[k] * hsi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, pa_x, pc_x, pc_z, gsk0_366, gsi_286, gsk1_366, \
                         hsh0_210, hsh1_210, hsi_282, hsi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_4 * hsh0_210[k]
                   - f_5 * hsh1_210[k]
                   + f_3 * pc_z[k] * hsi_282[k];

        t_366[k] = pa_x[k] * gsk0_366[k]
                   + f_16 * gsi_286[k]
                   - f_12 * pc_x[k] * gsk1_366[k];

        t_367[k] = f_3 * pc_z[k] * hsi_283[k];
    }
}

static auto
compute_prim_hsk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsk0,
                                                          const size_t gsi, const size_t gsk1,
                                                          const size_t hsh0, const size_t hsh1,
                                                          const size_t hsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_16 = 2.0 / q;
    const auto f_19 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsk0_216 = buffer.data(gsk0 + 216);
    const auto *gsk0_219 = buffer.data(gsk0 + 219);
    const auto *gsk0_222 = buffer.data(gsk0 + 222);
    const auto *gsk0_226 = buffer.data(gsk0 + 226);
    const auto *gsk0_231 = buffer.data(gsk0 + 231);
    const auto *gsk0_324 = buffer.data(gsk0 + 324);
    const auto *gsk0_329 = buffer.data(gsk0 + 329);
    const auto *gsk0_333 = buffer.data(gsk0 + 333);
    const auto *gsk0_338 = buffer.data(gsk0 + 338);
    const auto *gsk0_344 = buffer.data(gsk0 + 344);
    const auto *gsk0_370 = buffer.data(gsk0 + 370);
    const auto *gsk0_375 = buffer.data(gsk0 + 375);
    const auto *gsk0_388 = buffer.data(gsk0 + 388);
    const auto *gsk0_390 = buffer.data(gsk0 + 390);
    const auto *gsk0_391 = buffer.data(gsk0 + 391);
    const auto *gsk0_392 = buffer.data(gsk0 + 392);
    const auto *gsk0_393 = buffer.data(gsk0 + 393);
    const auto *gsk0_395 = buffer.data(gsk0 + 395);
    const auto *gsk0_401 = buffer.data(gsk0 + 401);
    const auto *gsk0_405 = buffer.data(gsk0 + 405);
    const auto *gsk0_408 = buffer.data(gsk0 + 408);
    const auto *gsk0_410 = buffer.data(gsk0 + 410);
    const auto *gsk0_413 = buffer.data(gsk0 + 413);
    const auto *gsk0_414 = buffer.data(gsk0 + 414);
    const auto *gsk0_416 = buffer.data(gsk0 + 416);
    const auto *gsk0_424 = buffer.data(gsk0 + 424);
    const auto *gsk0_426 = buffer.data(gsk0 + 426);
    const auto *gsk0_427 = buffer.data(gsk0 + 427);
    const auto *gsk0_428 = buffer.data(gsk0 + 428);
    const auto *gsk0_429 = buffer.data(gsk0 + 429);
    const auto *gsk0_431 = buffer.data(gsk0 + 431);
    const auto *gsk0_432 = buffer.data(gsk0 + 432);
    const auto *gsk0_435 = buffer.data(gsk0 + 435);
    const auto *gsk0_437 = buffer.data(gsk0 + 437);
    const auto *gsk0_438 = buffer.data(gsk0 + 438);
    const auto *gsk0_441 = buffer.data(gsk0 + 441);
    const auto *gsk0_442 = buffer.data(gsk0 + 442);
    const auto *gsk0_444 = buffer.data(gsk0 + 444);
    const auto *gsk0_446 = buffer.data(gsk0 + 446);
    const auto *gsk0_447 = buffer.data(gsk0 + 447);
    const auto *gsk0_449 = buffer.data(gsk0 + 449);
    const auto *gsk0_450 = buffer.data(gsk0 + 450);
    const auto *gsk0_452 = buffer.data(gsk0 + 452);
    const auto *gsk0_460 = buffer.data(gsk0 + 460);
    const auto *gsk0_462 = buffer.data(gsk0 + 462);
    const auto *gsk0_463 = buffer.data(gsk0 + 463);
    const auto *gsk0_464 = buffer.data(gsk0 + 464);
    const auto *gsk0_465 = buffer.data(gsk0 + 465);
    const auto *gsk0_467 = buffer.data(gsk0 + 467);
    const auto *gsk0_471 = buffer.data(gsk0 + 471);
    const auto *gsk0_474 = buffer.data(gsk0 + 474);
    const auto *gsk0_478 = buffer.data(gsk0 + 478);
    const auto *gsk0_480 = buffer.data(gsk0 + 480);
    const auto *gsk0_483 = buffer.data(gsk0 + 483);
    const auto *gsk0_485 = buffer.data(gsk0 + 485);
    const auto *gsk0_486 = buffer.data(gsk0 + 486);

    const auto *gsi_168 = buffer.data(gsi + 168);
    const auto *gsi_171 = buffer.data(gsi + 171);
    const auto *gsi_173 = buffer.data(gsi + 173);
    const auto *gsi_174 = buffer.data(gsi + 174);
    const auto *gsi_177 = buffer.data(gsi + 177);
    const auto *gsi_178 = buffer.data(gsi + 178);
    const auto *gsi_182 = buffer.data(gsi + 182);
    const auto *gsi_189 = buffer.data(gsi + 189);
    const auto *gsi_195 = buffer.data(gsi + 195);
    const auto *gsi_196 = buffer.data(gsi + 196);
    const auto *gsi_198 = buffer.data(gsi + 198);
    const auto *gsi_199 = buffer.data(gsi + 199);
    const auto *gsi_201 = buffer.data(gsi + 201);
    const auto *gsi_202 = buffer.data(gsi + 202);
    const auto *gsi_205 = buffer.data(gsi + 205);
    const auto *gsi_206 = buffer.data(gsi + 206);
    const auto *gsi_210 = buffer.data(gsi + 210);
    const auto *gsi_217 = buffer.data(gsi + 217);
    const auto *gsi_223 = buffer.data(gsi + 223);
    const auto *gsi_224 = buffer.data(gsi + 224);
    const auto *gsi_226 = buffer.data(gsi + 226);
    const auto *gsi_227 = buffer.data(gsi + 227);
    const auto *gsi_229 = buffer.data(gsi + 229);
    const auto *gsi_230 = buffer.data(gsi + 230);
    const auto *gsi_233 = buffer.data(gsi + 233);
    const auto *gsi_234 = buffer.data(gsi + 234);
    const auto *gsi_238 = buffer.data(gsi + 238);
    const auto *gsi_251 = buffer.data(gsi + 251);
    const auto *gsi_252 = buffer.data(gsi + 252);
    const auto *gsi_254 = buffer.data(gsi + 254);
    const auto *gsi_257 = buffer.data(gsi + 257);
    const auto *gsi_261 = buffer.data(gsi + 261);
    const auto *gsi_266 = buffer.data(gsi + 266);
    const auto *gsi_290 = buffer.data(gsi + 290);
    const auto *gsi_295 = buffer.data(gsi + 295);
    const auto *gsi_301 = buffer.data(gsi + 301);
    const auto *gsi_303 = buffer.data(gsi + 303);
    const auto *gsi_304 = buffer.data(gsi + 304);
    const auto *gsi_305 = buffer.data(gsi + 305);
    const auto *gsi_306 = buffer.data(gsi + 306);
    const auto *gsi_307 = buffer.data(gsi + 307);
    const auto *gsi_313 = buffer.data(gsi + 313);
    const auto *gsi_317 = buffer.data(gsi + 317);
    const auto *gsi_320 = buffer.data(gsi + 320);
    const auto *gsi_322 = buffer.data(gsi + 322);
    const auto *gsi_325 = buffer.data(gsi + 325);
    const auto *gsi_326 = buffer.data(gsi + 326);
    const auto *gsi_328 = buffer.data(gsi + 328);
    const auto *gsi_329 = buffer.data(gsi + 329);
    const auto *gsi_330 = buffer.data(gsi + 330);
    const auto *gsi_331 = buffer.data(gsi + 331);
    const auto *gsi_332 = buffer.data(gsi + 332);
    const auto *gsi_333 = buffer.data(gsi + 333);
    const auto *gsi_334 = buffer.data(gsi + 334);
    const auto *gsi_335 = buffer.data(gsi + 335);
    const auto *gsi_336 = buffer.data(gsi + 336);
    const auto *gsi_339 = buffer.data(gsi + 339);
    const auto *gsi_341 = buffer.data(gsi + 341);
    const auto *gsi_342 = buffer.data(gsi + 342);
    const auto *gsi_345 = buffer.data(gsi + 345);
    const auto *gsi_346 = buffer.data(gsi + 346);
    const auto *gsi_348 = buffer.data(gsi + 348);
    const auto *gsi_350 = buffer.data(gsi + 350);
    const auto *gsi_351 = buffer.data(gsi + 351);
    const auto *gsi_353 = buffer.data(gsi + 353);
    const auto *gsi_354 = buffer.data(gsi + 354);
    const auto *gsi_356 = buffer.data(gsi + 356);
    const auto *gsi_357 = buffer.data(gsi + 357);
    const auto *gsi_358 = buffer.data(gsi + 358);
    const auto *gsi_359 = buffer.data(gsi + 359);
    const auto *gsi_360 = buffer.data(gsi + 360);
    const auto *gsi_361 = buffer.data(gsi + 361);
    const auto *gsi_362 = buffer.data(gsi + 362);
    const auto *gsi_363 = buffer.data(gsi + 363);
    const auto *gsi_367 = buffer.data(gsi + 367);
    const auto *gsi_370 = buffer.data(gsi + 370);
    const auto *gsi_374 = buffer.data(gsi + 374);
    const auto *gsi_376 = buffer.data(gsi + 376);
    const auto *gsi_379 = buffer.data(gsi + 379);
    const auto *gsi_381 = buffer.data(gsi + 381);
    const auto *gsi_382 = buffer.data(gsi + 382);

    const auto *gsk1_216 = buffer.data(gsk1 + 216);
    const auto *gsk1_219 = buffer.data(gsk1 + 219);
    const auto *gsk1_222 = buffer.data(gsk1 + 222);
    const auto *gsk1_226 = buffer.data(gsk1 + 226);
    const auto *gsk1_231 = buffer.data(gsk1 + 231);
    const auto *gsk1_324 = buffer.data(gsk1 + 324);
    const auto *gsk1_329 = buffer.data(gsk1 + 329);
    const auto *gsk1_333 = buffer.data(gsk1 + 333);
    const auto *gsk1_338 = buffer.data(gsk1 + 338);
    const auto *gsk1_344 = buffer.data(gsk1 + 344);
    const auto *gsk1_370 = buffer.data(gsk1 + 370);
    const auto *gsk1_375 = buffer.data(gsk1 + 375);
    const auto *gsk1_388 = buffer.data(gsk1 + 388);
    const auto *gsk1_390 = buffer.data(gsk1 + 390);
    const auto *gsk1_391 = buffer.data(gsk1 + 391);
    const auto *gsk1_392 = buffer.data(gsk1 + 392);
    const auto *gsk1_393 = buffer.data(gsk1 + 393);
    const auto *gsk1_395 = buffer.data(gsk1 + 395);
    const auto *gsk1_401 = buffer.data(gsk1 + 401);
    const auto *gsk1_405 = buffer.data(gsk1 + 405);
    const auto *gsk1_408 = buffer.data(gsk1 + 408);
    const auto *gsk1_410 = buffer.data(gsk1 + 410);
    const auto *gsk1_413 = buffer.data(gsk1 + 413);
    const auto *gsk1_414 = buffer.data(gsk1 + 414);
    const auto *gsk1_416 = buffer.data(gsk1 + 416);
    const auto *gsk1_424 = buffer.data(gsk1 + 424);
    const auto *gsk1_426 = buffer.data(gsk1 + 426);
    const auto *gsk1_427 = buffer.data(gsk1 + 427);
    const auto *gsk1_428 = buffer.data(gsk1 + 428);
    const auto *gsk1_429 = buffer.data(gsk1 + 429);
    const auto *gsk1_431 = buffer.data(gsk1 + 431);
    const auto *gsk1_432 = buffer.data(gsk1 + 432);
    const auto *gsk1_435 = buffer.data(gsk1 + 435);
    const auto *gsk1_437 = buffer.data(gsk1 + 437);
    const auto *gsk1_438 = buffer.data(gsk1 + 438);
    const auto *gsk1_441 = buffer.data(gsk1 + 441);
    const auto *gsk1_442 = buffer.data(gsk1 + 442);
    const auto *gsk1_444 = buffer.data(gsk1 + 444);
    const auto *gsk1_446 = buffer.data(gsk1 + 446);
    const auto *gsk1_447 = buffer.data(gsk1 + 447);
    const auto *gsk1_449 = buffer.data(gsk1 + 449);
    const auto *gsk1_450 = buffer.data(gsk1 + 450);
    const auto *gsk1_452 = buffer.data(gsk1 + 452);
    const auto *gsk1_460 = buffer.data(gsk1 + 460);
    const auto *gsk1_462 = buffer.data(gsk1 + 462);
    const auto *gsk1_463 = buffer.data(gsk1 + 463);
    const auto *gsk1_464 = buffer.data(gsk1 + 464);
    const auto *gsk1_465 = buffer.data(gsk1 + 465);
    const auto *gsk1_467 = buffer.data(gsk1 + 467);
    const auto *gsk1_471 = buffer.data(gsk1 + 471);
    const auto *gsk1_474 = buffer.data(gsk1 + 474);
    const auto *gsk1_478 = buffer.data(gsk1 + 478);
    const auto *gsk1_480 = buffer.data(gsk1 + 480);
    const auto *gsk1_483 = buffer.data(gsk1 + 483);
    const auto *gsk1_485 = buffer.data(gsk1 + 485);
    const auto *gsk1_486 = buffer.data(gsk1 + 486);

    const auto *hsh0_212 = buffer.data(hsh0 + 212);
    const auto *hsh0_213 = buffer.data(hsh0 + 213);
    const auto *hsh0_215 = buffer.data(hsh0 + 215);
    const auto *hsh0_216 = buffer.data(hsh0 + 216);
    const auto *hsh0_217 = buffer.data(hsh0 + 217);
    const auto *hsh0_219 = buffer.data(hsh0 + 219);

    const auto *hsh1_212 = buffer.data(hsh1 + 212);
    const auto *hsh1_213 = buffer.data(hsh1 + 213);
    const auto *hsh1_215 = buffer.data(hsh1 + 215);
    const auto *hsh1_216 = buffer.data(hsh1 + 216);
    const auto *hsh1_217 = buffer.data(hsh1 + 217);
    const auto *hsh1_219 = buffer.data(hsh1 + 219);

    const auto *hsi_285 = buffer.data(hsi + 285);
    const auto *hsi_286 = buffer.data(hsi + 286);
    const auto *hsi_287 = buffer.data(hsi + 287);
    const auto *hsi_289 = buffer.data(hsi + 289);
    const auto *hsi_290 = buffer.data(hsi + 290);
    const auto *hsi_291 = buffer.data(hsi + 291);
    const auto *hsi_292 = buffer.data(hsi + 292);
    const auto *hsi_294 = buffer.data(hsi + 294);
    const auto *hsi_295 = buffer.data(hsi + 295);
    const auto *hsi_301 = buffer.data(hsi + 301);
    const auto *hsi_303 = buffer.data(hsi + 303);
    const auto *hsi_304 = buffer.data(hsi + 304);
    const auto *hsi_305 = buffer.data(hsi + 305);
    const auto *hsi_306 = buffer.data(hsi + 306);
    const auto *hsi_307 = buffer.data(hsi + 307);
    const auto *hsi_308 = buffer.data(hsi + 308);
    const auto *hsi_310 = buffer.data(hsi + 310);
    const auto *hsi_311 = buffer.data(hsi + 311);
    const auto *hsi_313 = buffer.data(hsi + 313);
    const auto *hsi_314 = buffer.data(hsi + 314);
    const auto *hsi_317 = buffer.data(hsi + 317);
    const auto *hsi_318 = buffer.data(hsi + 318);
    const auto *hsi_322 = buffer.data(hsi + 322);
    const auto *hsi_329 = buffer.data(hsi + 329);
    const auto *hsi_330 = buffer.data(hsi + 330);
    const auto *hsi_331 = buffer.data(hsi + 331);
    const auto *hsi_332 = buffer.data(hsi + 332);
    const auto *hsi_333 = buffer.data(hsi + 333);
    const auto *hsi_334 = buffer.data(hsi + 334);
    const auto *hsi_335 = buffer.data(hsi + 335);
    const auto *hsi_336 = buffer.data(hsi + 336);
    const auto *hsi_338 = buffer.data(hsi + 338);
    const auto *hsi_339 = buffer.data(hsi + 339);
    const auto *hsi_341 = buffer.data(hsi + 341);
    const auto *hsi_342 = buffer.data(hsi + 342);
    const auto *hsi_345 = buffer.data(hsi + 345);
    const auto *hsi_346 = buffer.data(hsi + 346);
    const auto *hsi_350 = buffer.data(hsi + 350);
    const auto *hsi_357 = buffer.data(hsi + 357);
    const auto *hsi_358 = buffer.data(hsi + 358);
    const auto *hsi_359 = buffer.data(hsi + 359);
    const auto *hsi_360 = buffer.data(hsi + 360);
    const auto *hsi_361 = buffer.data(hsi + 361);
    const auto *hsi_362 = buffer.data(hsi + 362);
    const auto *hsi_363 = buffer.data(hsi + 363);
    const auto *hsi_364 = buffer.data(hsi + 364);
    const auto *hsi_366 = buffer.data(hsi + 366);
    const auto *hsi_367 = buffer.data(hsi + 367);
    const auto *hsi_369 = buffer.data(hsi + 369);
    const auto *hsi_370 = buffer.data(hsi + 370);
    const auto *hsi_373 = buffer.data(hsi + 373);
    const auto *hsi_374 = buffer.data(hsi + 374);
    const auto *hsi_378 = buffer.data(hsi + 378);

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pa_x, pc_x, pc_y, pc_z, gsk0_370, \
                         gsi_173, gsi_290, gsk1_370, hsh0_212, hsh1_212, hsi_285, \
                         hsi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_16 * gsi_173[k]
                   + f_3 * pc_y[k] * hsi_285[k];

        t_369[k] = f_6 * hsh0_212[k]
                   - f_7 * hsh1_212[k]
                   + f_3 * pc_z[k] * hsi_285[k];

        t_370[k] = pa_x[k] * gsk0_370[k]
                   + f_15 * gsi_290[k]
                   - f_12 * pc_x[k] * gsk1_370[k];

        t_371[k] = f_3 * pc_z[k] * hsi_286[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, gsi_177, hsh0_213, hsh0_215, \
                         hsh1_213, hsh1_215, hsi_287, hsi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_4 * hsh0_213[k]
                   - f_5 * hsh1_213[k]
                   + f_3 * pc_z[k] * hsi_287[k];

        t_373[k] = f_16 * gsi_177[k]
                   + f_3 * pc_y[k] * hsi_289[k];

        t_374[k] = f_8 * hsh0_215[k]
                   - f_9 * hsh1_215[k]
                   + f_3 * pc_z[k] * hsi_289[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_x, pc_x, pc_z, gsk0_375, gsi_295, gsk1_375, \
                         hsh0_216, hsh1_216, hsi_290, hsi_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pa_x[k] * gsk0_375[k]
                   + f_14 * gsi_295[k]
                   - f_12 * pc_x[k] * gsk1_375[k];

        t_376[k] = f_3 * pc_z[k] * hsi_290[k];

        t_377[k] = f_4 * hsh0_216[k]
                   - f_5 * hsh1_216[k]
                   + f_3 * pc_z[k] * hsi_291[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, pc_y, pc_z, gsi_182, gsi_301, \
                         hsh0_217, hsh0_219, hsh1_217, hsh1_219, hsi_292, hsi_294, \
                         hsi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_6 * hsh0_217[k]
                   - f_7 * hsh1_217[k]
                   + f_3 * pc_z[k] * hsi_292[k];

        t_379[k] = f_16 * gsi_182[k]
                   + f_3 * pc_y[k] * hsi_294[k];

        t_380[k] = f_10 * hsh0_219[k]
                   - f_11 * hsh1_219[k]
                   + f_3 * pc_z[k] * hsi_294[k];

        t_381[k] = f_13 * gsi_301[k]
                   + f_3 * pc_x[k] * hsi_301[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_z, gsi_303, gsi_304, \
                         gsi_305, gsi_306, hsi_295, hsi_303, hsi_304, hsi_305, \
                         hsi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_3 * pc_z[k] * hsi_295[k];

        t_383[k] = f_13 * gsi_303[k]
                   + f_3 * pc_x[k] * hsi_303[k];

        t_384[k] = f_13 * gsi_304[k]
                   + f_3 * pc_x[k] * hsi_304[k];

        t_385[k] = f_13 * gsi_305[k]
                   + f_3 * pc_x[k] * hsi_305[k];

        t_386[k] = f_13 * gsi_306[k]
                   + f_3 * pc_x[k] * hsi_306[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_x, pc_x, pc_z, gsk0_388, gsk0_390, \
                         gsi_307, gsk1_388, gsk1_390, hsi_301, \
                         hsi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_13 * gsi_307[k]
                   + f_3 * pc_x[k] * hsi_307[k];

        t_388[k] = pa_x[k] * gsk0_388[k]
                   - f_12 * pc_x[k] * gsk1_388[k];

        t_389[k] = f_3 * pc_z[k] * hsi_301[k];

        t_390[k] = pa_x[k] * gsk0_390[k]
                   - f_12 * pc_x[k] * gsk1_390[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_x, pc_x, pc_y, gsk0_391, gsk0_392, \
                         gsk0_393, gsi_195, gsk1_391, gsk1_392, gsk1_393, \
                         hsi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_x[k] * gsk0_391[k]
                   - f_12 * pc_x[k] * gsk1_391[k];

        t_392[k] = pa_x[k] * gsk0_392[k]
                   - f_12 * pc_x[k] * gsk1_392[k];

        t_393[k] = pa_x[k] * gsk0_393[k]
                   - f_12 * pc_x[k] * gsk1_393[k];

        t_394[k] = f_16 * gsi_195[k]
                   + f_3 * pc_y[k] * hsi_307[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_x, pa_z, pc_x, pc_y, pc_z, gsk0_216, \
                         gsk0_395, gsi_168, gsi_196, gsk1_216, gsk1_395, \
                         hsi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_x[k] * gsk0_395[k]
                   - f_12 * pc_x[k] * gsk1_395[k];

        t_396[k] = pa_z[k] * gsk0_216[k]
                   - f_12 * pc_z[k] * gsk1_216[k];

        t_397[k] = f_15 * gsi_196[k]
                   + f_3 * pc_y[k] * hsi_308[k];

        t_398[k] = f_13 * gsi_168[k]
                   + f_3 * pc_z[k] * hsi_308[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pa_x, pa_z, pc_x, pc_y, pc_z, gsk0_219, \
                         gsk0_401, gsi_198, gsi_313, gsk1_219, gsk1_401, \
                         hsi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = pa_z[k] * gsk0_219[k]
                   - f_12 * pc_z[k] * gsk1_219[k];

        t_400[k] = f_15 * gsi_198[k]
                   + f_3 * pc_y[k] * hsi_310[k];

        t_401[k] = pa_x[k] * gsk0_401[k]
                   + f_0 * gsi_313[k]
                   - f_12 * pc_x[k] * gsk1_401[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_z, pc_y, pc_z, gsk0_222, gsi_171, gsi_201, \
                         gsk1_222, hsi_311, hsi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_z[k] * gsk0_222[k]
                   - f_12 * pc_z[k] * gsk1_222[k];

        t_403[k] = f_13 * gsi_171[k]
                   + f_3 * pc_z[k] * hsi_311[k];

        t_404[k] = f_15 * gsi_201[k]
                   + f_3 * pc_y[k] * hsi_313[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, pa_x, pa_z, pc_x, pc_z, gsk0_226, gsk0_405, \
                         gsi_174, gsi_317, gsk1_226, gsk1_405, \
                         hsi_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_x[k] * gsk0_405[k]
                   + f_16 * gsi_317[k]
                   - f_12 * pc_x[k] * gsk1_405[k];

        t_406[k] = pa_z[k] * gsk0_226[k]
                   - f_12 * pc_z[k] * gsk1_226[k];

        t_407[k] = f_13 * gsi_174[k]
                   + f_3 * pc_z[k] * hsi_314[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pa_x, pc_x, pc_y, gsk0_408, gsk0_410, gsi_205, \
                         gsi_320, gsi_322, gsk1_408, gsk1_410, \
                         hsi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = pa_x[k] * gsk0_408[k]
                   + f_15 * gsi_320[k]
                   - f_12 * pc_x[k] * gsk1_408[k];

        t_409[k] = f_15 * gsi_205[k]
                   + f_3 * pc_y[k] * hsi_317[k];

        t_410[k] = pa_x[k] * gsk0_410[k]
                   + f_15 * gsi_322[k]
                   - f_12 * pc_x[k] * gsk1_410[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pa_x, pa_z, pc_x, pc_z, gsk0_231, gsk0_413, \
                         gsi_178, gsi_325, gsk1_231, gsk1_413, \
                         hsi_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = pa_z[k] * gsk0_231[k]
                   - f_12 * pc_z[k] * gsk1_231[k];

        t_412[k] = f_13 * gsi_178[k]
                   + f_3 * pc_z[k] * hsi_318[k];

        t_413[k] = pa_x[k] * gsk0_413[k]
                   + f_14 * gsi_325[k]
                   - f_12 * pc_x[k] * gsk1_413[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pa_x, pc_x, pc_y, gsk0_414, gsk0_416, gsi_210, \
                         gsi_326, gsi_328, gsk1_414, gsk1_416, \
                         hsi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = pa_x[k] * gsk0_414[k]
                   + f_14 * gsi_326[k]
                   - f_12 * pc_x[k] * gsk1_414[k];

        t_415[k] = f_15 * gsi_210[k]
                   + f_3 * pc_y[k] * hsi_322[k];

        t_416[k] = pa_x[k] * gsk0_416[k]
                   + f_14 * gsi_328[k]
                   - f_12 * pc_x[k] * gsk1_416[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, pc_x, gsi_329, gsi_330, gsi_331, \
                         gsi_332, gsi_333, hsi_329, hsi_330, hsi_331, hsi_332, \
                         hsi_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_13 * gsi_329[k]
                   + f_3 * pc_x[k] * hsi_329[k];

        t_418[k] = f_13 * gsi_330[k]
                   + f_3 * pc_x[k] * hsi_330[k];

        t_419[k] = f_13 * gsi_331[k]
                   + f_3 * pc_x[k] * hsi_331[k];

        t_420[k] = f_13 * gsi_332[k]
                   + f_3 * pc_x[k] * hsi_332[k];

        t_421[k] = f_13 * gsi_333[k]
                   + f_3 * pc_x[k] * hsi_333[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pa_x, pc_x, pc_z, gsk0_424, gsi_189, \
                         gsi_334, gsi_335, gsk1_424, hsi_329, hsi_334, \
                         hsi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_13 * gsi_334[k]
                   + f_3 * pc_x[k] * hsi_334[k];

        t_423[k] = f_13 * gsi_335[k]
                   + f_3 * pc_x[k] * hsi_335[k];

        t_424[k] = pa_x[k] * gsk0_424[k]
                   - f_12 * pc_x[k] * gsk1_424[k];

        t_425[k] = f_13 * gsi_189[k]
                   + f_3 * pc_z[k] * hsi_329[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pa_x, pc_x, gsk0_426, gsk0_427, gsk0_428, \
                         gsk0_429, gsk1_426, gsk1_427, gsk1_428, \
                         gsk1_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = pa_x[k] * gsk0_426[k]
                   - f_12 * pc_x[k] * gsk1_426[k];

        t_427[k] = pa_x[k] * gsk0_427[k]
                   - f_12 * pc_x[k] * gsk1_427[k];

        t_428[k] = pa_x[k] * gsk0_428[k]
                   - f_12 * pc_x[k] * gsk1_428[k];

        t_429[k] = pa_x[k] * gsk0_429[k]
                   - f_12 * pc_x[k] * gsk1_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pa_x, pc_x, pc_y, gsk0_431, gsk0_432, \
                         gsi_223, gsi_224, gsi_336, gsk1_431, gsk1_432, hsi_335, \
                         hsi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * gsi_223[k]
                   + f_3 * pc_y[k] * hsi_335[k];

        t_431[k] = pa_x[k] * gsk0_431[k]
                   - f_12 * pc_x[k] * gsk1_431[k];

        t_432[k] = pa_x[k] * gsk0_432[k]
                   + f_19 * gsi_336[k]
                   - f_12 * pc_x[k] * gsk1_432[k];

        t_433[k] = f_14 * gsi_224[k]
                   + f_3 * pc_y[k] * hsi_336[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, pa_x, pc_x, pc_y, pc_z, gsk0_435, gsi_196, \
                         gsi_226, gsi_339, gsk1_435, hsi_336, hsi_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_14 * gsi_196[k]
                   + f_3 * pc_z[k] * hsi_336[k];

        t_435[k] = pa_x[k] * gsk0_435[k]
                   + f_0 * gsi_339[k]
                   - f_12 * pc_x[k] * gsk1_435[k];

        t_436[k] = f_14 * gsi_226[k]
                   + f_3 * pc_y[k] * hsi_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pa_x, pc_x, pc_z, gsk0_437, gsk0_438, gsi_199, \
                         gsi_341, gsi_342, gsk1_437, gsk1_438, \
                         hsi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = pa_x[k] * gsk0_437[k]
                   + f_0 * gsi_341[k]
                   - f_12 * pc_x[k] * gsk1_437[k];

        t_438[k] = pa_x[k] * gsk0_438[k]
                   + f_16 * gsi_342[k]
                   - f_12 * pc_x[k] * gsk1_438[k];

        t_439[k] = f_14 * gsi_199[k]
                   + f_3 * pc_z[k] * hsi_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pa_x, pc_x, pc_y, gsk0_441, gsk0_442, gsi_229, \
                         gsi_345, gsi_346, gsk1_441, gsk1_442, \
                         hsi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * gsi_229[k]
                   + f_3 * pc_y[k] * hsi_341[k];

        t_441[k] = pa_x[k] * gsk0_441[k]
                   + f_16 * gsi_345[k]
                   - f_12 * pc_x[k] * gsk1_441[k];

        t_442[k] = pa_x[k] * gsk0_442[k]
                   + f_15 * gsi_346[k]
                   - f_12 * pc_x[k] * gsk1_442[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pa_x, pc_x, pc_y, pc_z, gsk0_444, gsi_202, \
                         gsi_233, gsi_348, gsk1_444, hsi_342, hsi_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * gsi_202[k]
                   + f_3 * pc_z[k] * hsi_342[k];

        t_444[k] = pa_x[k] * gsk0_444[k]
                   + f_15 * gsi_348[k]
                   - f_12 * pc_x[k] * gsk1_444[k];

        t_445[k] = f_14 * gsi_233[k]
                   + f_3 * pc_y[k] * hsi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pa_x, pc_x, pc_z, gsk0_446, gsk0_447, gsi_206, \
                         gsi_350, gsi_351, gsk1_446, gsk1_447, \
                         hsi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = pa_x[k] * gsk0_446[k]
                   + f_15 * gsi_350[k]
                   - f_12 * pc_x[k] * gsk1_446[k];

        t_447[k] = pa_x[k] * gsk0_447[k]
                   + f_14 * gsi_351[k]
                   - f_12 * pc_x[k] * gsk1_447[k];

        t_448[k] = f_14 * gsi_206[k]
                   + f_3 * pc_z[k] * hsi_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pa_x, pc_x, pc_y, gsk0_449, gsk0_450, gsi_238, \
                         gsi_353, gsi_354, gsk1_449, gsk1_450, \
                         hsi_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pa_x[k] * gsk0_449[k]
                   + f_14 * gsi_353[k]
                   - f_12 * pc_x[k] * gsk1_449[k];

        t_450[k] = pa_x[k] * gsk0_450[k]
                   + f_14 * gsi_354[k]
                   - f_12 * pc_x[k] * gsk1_450[k];

        t_451[k] = f_14 * gsi_238[k]
                   + f_3 * pc_y[k] * hsi_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pa_x, pc_x, gsk0_452, gsi_356, gsi_357, \
                         gsi_358, gsi_359, gsk1_452, hsi_357, hsi_358, \
                         hsi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = pa_x[k] * gsk0_452[k]
                   + f_14 * gsi_356[k]
                   - f_12 * pc_x[k] * gsk1_452[k];

        t_453[k] = f_13 * gsi_357[k]
                   + f_3 * pc_x[k] * hsi_357[k];

        t_454[k] = f_13 * gsi_358[k]
                   + f_3 * pc_x[k] * hsi_358[k];

        t_455[k] = f_13 * gsi_359[k]
                   + f_3 * pc_x[k] * hsi_359[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, gsi_360, gsi_361, gsi_362, gsi_363, \
                         hsi_360, hsi_361, hsi_362, hsi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_13 * gsi_360[k]
                   + f_3 * pc_x[k] * hsi_360[k];

        t_457[k] = f_13 * gsi_361[k]
                   + f_3 * pc_x[k] * hsi_361[k];

        t_458[k] = f_13 * gsi_362[k]
                   + f_3 * pc_x[k] * hsi_362[k];

        t_459[k] = f_13 * gsi_363[k]
                   + f_3 * pc_x[k] * hsi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pa_x, pc_x, pc_z, gsk0_460, gsk0_462, \
                         gsk0_463, gsi_217, gsk1_460, gsk1_462, gsk1_463, \
                         hsi_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = pa_x[k] * gsk0_460[k]
                   - f_12 * pc_x[k] * gsk1_460[k];

        t_461[k] = f_14 * gsi_217[k]
                   + f_3 * pc_z[k] * hsi_357[k];

        t_462[k] = pa_x[k] * gsk0_462[k]
                   - f_12 * pc_x[k] * gsk1_462[k];

        t_463[k] = pa_x[k] * gsk0_463[k]
                   - f_12 * pc_x[k] * gsk1_463[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, pa_x, pc_x, pc_y, gsk0_464, gsk0_465, \
                         gsk0_467, gsi_251, gsk1_464, gsk1_465, gsk1_467, \
                         hsi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = pa_x[k] * gsk0_464[k]
                   - f_12 * pc_x[k] * gsk1_464[k];

        t_465[k] = pa_x[k] * gsk0_465[k]
                   - f_12 * pc_x[k] * gsk1_465[k];

        t_466[k] = f_14 * gsi_251[k]
                   + f_3 * pc_y[k] * hsi_363[k];

        t_467[k] = pa_x[k] * gsk0_467[k]
                   - f_12 * pc_x[k] * gsk1_467[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pa_y, pc_y, pc_z, gsk0_324, gsi_224, gsi_252, \
                         gsk1_324, hsi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pa_y[k] * gsk0_324[k]
                   - f_12 * pc_y[k] * gsk1_324[k];

        t_469[k] = f_13 * gsi_252[k]
                   + f_3 * pc_y[k] * hsi_364[k];

        t_470[k] = f_15 * gsi_224[k]
                   + f_3 * pc_z[k] * hsi_364[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pa_x, pa_y, pc_x, pc_y, gsk0_329, gsk0_471, \
                         gsi_254, gsi_367, gsk1_329, gsk1_471, \
                         hsi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pa_x[k] * gsk0_471[k]
                   + f_0 * gsi_367[k]
                   - f_12 * pc_x[k] * gsk1_471[k];

        t_472[k] = f_13 * gsi_254[k]
                   + f_3 * pc_y[k] * hsi_366[k];

        t_473[k] = pa_y[k] * gsk0_329[k]
                   - f_12 * pc_y[k] * gsk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pa_x, pc_x, pc_y, pc_z, gsk0_474, gsi_227, \
                         gsi_257, gsi_370, gsk1_474, hsi_367, hsi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pa_x[k] * gsk0_474[k]
                   + f_16 * gsi_370[k]
                   - f_12 * pc_x[k] * gsk1_474[k];

        t_475[k] = f_15 * gsi_227[k]
                   + f_3 * pc_z[k] * hsi_367[k];

        t_476[k] = f_13 * gsi_257[k]
                   + f_3 * pc_y[k] * hsi_369[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, pa_x, pa_y, pc_x, pc_y, pc_z, gsk0_333, \
                         gsk0_478, gsi_230, gsi_374, gsk1_333, gsk1_478, \
                         hsi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = pa_y[k] * gsk0_333[k]
                   - f_12 * pc_y[k] * gsk1_333[k];

        t_478[k] = pa_x[k] * gsk0_478[k]
                   + f_15 * gsi_374[k]
                   - f_12 * pc_x[k] * gsk1_478[k];

        t_479[k] = f_15 * gsi_230[k]
                   + f_3 * pc_z[k] * hsi_370[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, pa_x, pa_y, pc_x, pc_y, gsk0_338, gsk0_480, \
                         gsi_261, gsi_376, gsk1_338, gsk1_480, \
                         hsi_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = pa_x[k] * gsk0_480[k]
                   + f_15 * gsi_376[k]
                   - f_12 * pc_x[k] * gsk1_480[k];

        t_481[k] = f_13 * gsi_261[k]
                   + f_3 * pc_y[k] * hsi_373[k];

        t_482[k] = pa_y[k] * gsk0_338[k]
                   - f_12 * pc_y[k] * gsk1_338[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pa_x, pc_x, pc_z, gsk0_483, gsk0_485, gsi_234, \
                         gsi_379, gsi_381, gsk1_483, gsk1_485, \
                         hsi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = pa_x[k] * gsk0_483[k]
                   + f_14 * gsi_379[k]
                   - f_12 * pc_x[k] * gsk1_483[k];

        t_484[k] = f_15 * gsi_234[k]
                   + f_3 * pc_z[k] * hsi_374[k];

        t_485[k] = pa_x[k] * gsk0_485[k]
                   + f_14 * gsi_381[k]
                   - f_12 * pc_x[k] * gsk1_485[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pa_x, pa_y, pc_x, pc_y, gsk0_344, gsk0_486, \
                         gsi_266, gsi_382, gsk1_344, gsk1_486, \
                         hsi_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = pa_x[k] * gsk0_486[k]
                   + f_14 * gsi_382[k]
                   - f_12 * pc_x[k] * gsk1_486[k];

        t_487[k] = f_13 * gsi_266[k]
                   + f_3 * pc_y[k] * hsi_378[k];

        t_488[k] = pa_y[k] * gsk0_344[k]
                   - f_12 * pc_y[k] * gsk1_344[k];
    }
}

static auto
compute_prim_hsk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsk0,
                                                          const size_t gsi, const size_t gsk1,
                                                          const size_t hsh0, const size_t hsh1,
                                                          const size_t hsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsk0_360 = buffer.data(gsk0 + 360);
    const auto *gsk0_361 = buffer.data(gsk0 + 361);
    const auto *gsk0_363 = buffer.data(gsk0 + 363);
    const auto *gsk0_366 = buffer.data(gsk0 + 366);
    const auto *gsk0_370 = buffer.data(gsk0 + 370);
    const auto *gsk0_375 = buffer.data(gsk0 + 375);
    const auto *gsk0_388 = buffer.data(gsk0 + 388);
    const auto *gsk0_390 = buffer.data(gsk0 + 390);
    const auto *gsk0_391 = buffer.data(gsk0 + 391);
    const auto *gsk0_392 = buffer.data(gsk0 + 392);
    const auto *gsk0_393 = buffer.data(gsk0 + 393);
    const auto *gsk0_496 = buffer.data(gsk0 + 496);
    const auto *gsk0_498 = buffer.data(gsk0 + 498);
    const auto *gsk0_499 = buffer.data(gsk0 + 499);
    const auto *gsk0_500 = buffer.data(gsk0 + 500);
    const auto *gsk0_501 = buffer.data(gsk0 + 501);
    const auto *gsk0_503 = buffer.data(gsk0 + 503);
    const auto *gsk0_504 = buffer.data(gsk0 + 504);
    const auto *gsk0_509 = buffer.data(gsk0 + 509);
    const auto *gsk0_513 = buffer.data(gsk0 + 513);
    const auto *gsk0_518 = buffer.data(gsk0 + 518);
    const auto *gsk0_524 = buffer.data(gsk0 + 524);
    const auto *gsk0_532 = buffer.data(gsk0 + 532);
    const auto *gsk0_533 = buffer.data(gsk0 + 533);
    const auto *gsk0_534 = buffer.data(gsk0 + 534);
    const auto *gsk0_535 = buffer.data(gsk0 + 535);
    const auto *gsk0_536 = buffer.data(gsk0 + 536);
    const auto *gsk0_537 = buffer.data(gsk0 + 537);
    const auto *gsk0_539 = buffer.data(gsk0 + 539);

    const auto *gsi_245 = buffer.data(gsi + 245);
    const auto *gsi_252 = buffer.data(gsi + 252);
    const auto *gsi_279 = buffer.data(gsi + 279);
    const auto *gsi_301 = buffer.data(gsi + 301);
    const auto *gsi_302 = buffer.data(gsi + 302);
    const auto *gsi_303 = buffer.data(gsi + 303);
    const auto *gsi_304 = buffer.data(gsi + 304);
    const auto *gsi_305 = buffer.data(gsi + 305);
    const auto *gsi_307 = buffer.data(gsi + 307);
    const auto *gsi_335 = buffer.data(gsi + 335);
    const auto *gsi_385 = buffer.data(gsi + 385);
    const auto *gsi_386 = buffer.data(gsi + 386);
    const auto *gsi_387 = buffer.data(gsi + 387);
    const auto *gsi_388 = buffer.data(gsi + 388);
    const auto *gsi_389 = buffer.data(gsi + 389);
    const auto *gsi_390 = buffer.data(gsi + 390);
    const auto *gsi_391 = buffer.data(gsi + 391);
    const auto *gsi_392 = buffer.data(gsi + 392);
    const auto *gsi_397 = buffer.data(gsi + 397);
    const auto *gsi_401 = buffer.data(gsi + 401);
    const auto *gsi_406 = buffer.data(gsi + 406);
    const auto *gsi_412 = buffer.data(gsi + 412);
    const auto *gsi_413 = buffer.data(gsi + 413);
    const auto *gsi_414 = buffer.data(gsi + 414);
    const auto *gsi_415 = buffer.data(gsi + 415);
    const auto *gsi_416 = buffer.data(gsi + 416);
    const auto *gsi_417 = buffer.data(gsi + 417);
    const auto *gsi_419 = buffer.data(gsi + 419);

    const auto *gsk1_360 = buffer.data(gsk1 + 360);
    const auto *gsk1_361 = buffer.data(gsk1 + 361);
    const auto *gsk1_363 = buffer.data(gsk1 + 363);
    const auto *gsk1_366 = buffer.data(gsk1 + 366);
    const auto *gsk1_370 = buffer.data(gsk1 + 370);
    const auto *gsk1_375 = buffer.data(gsk1 + 375);
    const auto *gsk1_388 = buffer.data(gsk1 + 388);
    const auto *gsk1_390 = buffer.data(gsk1 + 390);
    const auto *gsk1_391 = buffer.data(gsk1 + 391);
    const auto *gsk1_392 = buffer.data(gsk1 + 392);
    const auto *gsk1_393 = buffer.data(gsk1 + 393);
    const auto *gsk1_496 = buffer.data(gsk1 + 496);
    const auto *gsk1_498 = buffer.data(gsk1 + 498);
    const auto *gsk1_499 = buffer.data(gsk1 + 499);
    const auto *gsk1_500 = buffer.data(gsk1 + 500);
    const auto *gsk1_501 = buffer.data(gsk1 + 501);
    const auto *gsk1_503 = buffer.data(gsk1 + 503);
    const auto *gsk1_504 = buffer.data(gsk1 + 504);
    const auto *gsk1_509 = buffer.data(gsk1 + 509);
    const auto *gsk1_513 = buffer.data(gsk1 + 513);
    const auto *gsk1_518 = buffer.data(gsk1 + 518);
    const auto *gsk1_524 = buffer.data(gsk1 + 524);
    const auto *gsk1_532 = buffer.data(gsk1 + 532);
    const auto *gsk1_533 = buffer.data(gsk1 + 533);
    const auto *gsk1_534 = buffer.data(gsk1 + 534);
    const auto *gsk1_535 = buffer.data(gsk1 + 535);
    const auto *gsk1_536 = buffer.data(gsk1 + 536);
    const auto *gsk1_537 = buffer.data(gsk1 + 537);
    const auto *gsk1_539 = buffer.data(gsk1 + 539);

    const auto *hsh0_294 = buffer.data(hsh0 + 294);
    const auto *hsh0_295 = buffer.data(hsh0 + 295);
    const auto *hsh0_296 = buffer.data(hsh0 + 296);
    const auto *hsh0_297 = buffer.data(hsh0 + 297);
    const auto *hsh0_298 = buffer.data(hsh0 + 298);
    const auto *hsh0_299 = buffer.data(hsh0 + 299);
    const auto *hsh0_300 = buffer.data(hsh0 + 300);
    const auto *hsh0_301 = buffer.data(hsh0 + 301);
    const auto *hsh0_302 = buffer.data(hsh0 + 302);
    const auto *hsh0_303 = buffer.data(hsh0 + 303);
    const auto *hsh0_315 = buffer.data(hsh0 + 315);
    const auto *hsh0_316 = buffer.data(hsh0 + 316);
    const auto *hsh0_318 = buffer.data(hsh0 + 318);
    const auto *hsh0_320 = buffer.data(hsh0 + 320);
    const auto *hsh0_321 = buffer.data(hsh0 + 321);
    const auto *hsh0_323 = buffer.data(hsh0 + 323);
    const auto *hsh0_324 = buffer.data(hsh0 + 324);
    const auto *hsh0_325 = buffer.data(hsh0 + 325);
    const auto *hsh0_327 = buffer.data(hsh0 + 327);
    const auto *hsh0_328 = buffer.data(hsh0 + 328);
    const auto *hsh0_329 = buffer.data(hsh0 + 329);
    const auto *hsh0_330 = buffer.data(hsh0 + 330);
    const auto *hsh0_331 = buffer.data(hsh0 + 331);
    const auto *hsh0_332 = buffer.data(hsh0 + 332);
    const auto *hsh0_333 = buffer.data(hsh0 + 333);
    const auto *hsh0_334 = buffer.data(hsh0 + 334);
    const auto *hsh0_335 = buffer.data(hsh0 + 335);
    const auto *hsh0_338 = buffer.data(hsh0 + 338);
    const auto *hsh0_340 = buffer.data(hsh0 + 340);
    const auto *hsh0_341 = buffer.data(hsh0 + 341);
    const auto *hsh0_343 = buffer.data(hsh0 + 343);
    const auto *hsh0_344 = buffer.data(hsh0 + 344);
    const auto *hsh0_345 = buffer.data(hsh0 + 345);
    const auto *hsh0_347 = buffer.data(hsh0 + 347);
    const auto *hsh0_348 = buffer.data(hsh0 + 348);
    const auto *hsh0_349 = buffer.data(hsh0 + 349);
    const auto *hsh0_350 = buffer.data(hsh0 + 350);
    const auto *hsh0_352 = buffer.data(hsh0 + 352);
    const auto *hsh0_353 = buffer.data(hsh0 + 353);
    const auto *hsh0_354 = buffer.data(hsh0 + 354);
    const auto *hsh0_355 = buffer.data(hsh0 + 355);
    const auto *hsh0_356 = buffer.data(hsh0 + 356);
    const auto *hsh0_357 = buffer.data(hsh0 + 357);
    const auto *hsh0_358 = buffer.data(hsh0 + 358);
    const auto *hsh0_359 = buffer.data(hsh0 + 359);

    const auto *hsh1_294 = buffer.data(hsh1 + 294);
    const auto *hsh1_295 = buffer.data(hsh1 + 295);
    const auto *hsh1_296 = buffer.data(hsh1 + 296);
    const auto *hsh1_297 = buffer.data(hsh1 + 297);
    const auto *hsh1_298 = buffer.data(hsh1 + 298);
    const auto *hsh1_299 = buffer.data(hsh1 + 299);
    const auto *hsh1_300 = buffer.data(hsh1 + 300);
    const auto *hsh1_301 = buffer.data(hsh1 + 301);
    const auto *hsh1_302 = buffer.data(hsh1 + 302);
    const auto *hsh1_303 = buffer.data(hsh1 + 303);
    const auto *hsh1_315 = buffer.data(hsh1 + 315);
    const auto *hsh1_316 = buffer.data(hsh1 + 316);
    const auto *hsh1_318 = buffer.data(hsh1 + 318);
    const auto *hsh1_320 = buffer.data(hsh1 + 320);
    const auto *hsh1_321 = buffer.data(hsh1 + 321);
    const auto *hsh1_323 = buffer.data(hsh1 + 323);
    const auto *hsh1_324 = buffer.data(hsh1 + 324);
    const auto *hsh1_325 = buffer.data(hsh1 + 325);
    const auto *hsh1_327 = buffer.data(hsh1 + 327);
    const auto *hsh1_328 = buffer.data(hsh1 + 328);
    const auto *hsh1_329 = buffer.data(hsh1 + 329);
    const auto *hsh1_330 = buffer.data(hsh1 + 330);
    const auto *hsh1_331 = buffer.data(hsh1 + 331);
    const auto *hsh1_332 = buffer.data(hsh1 + 332);
    const auto *hsh1_333 = buffer.data(hsh1 + 333);
    const auto *hsh1_334 = buffer.data(hsh1 + 334);
    const auto *hsh1_335 = buffer.data(hsh1 + 335);
    const auto *hsh1_338 = buffer.data(hsh1 + 338);
    const auto *hsh1_340 = buffer.data(hsh1 + 340);
    const auto *hsh1_341 = buffer.data(hsh1 + 341);
    const auto *hsh1_343 = buffer.data(hsh1 + 343);
    const auto *hsh1_344 = buffer.data(hsh1 + 344);
    const auto *hsh1_345 = buffer.data(hsh1 + 345);
    const auto *hsh1_347 = buffer.data(hsh1 + 347);
    const auto *hsh1_348 = buffer.data(hsh1 + 348);
    const auto *hsh1_349 = buffer.data(hsh1 + 349);
    const auto *hsh1_350 = buffer.data(hsh1 + 350);
    const auto *hsh1_352 = buffer.data(hsh1 + 352);
    const auto *hsh1_353 = buffer.data(hsh1 + 353);
    const auto *hsh1_354 = buffer.data(hsh1 + 354);
    const auto *hsh1_355 = buffer.data(hsh1 + 355);
    const auto *hsh1_356 = buffer.data(hsh1 + 356);
    const auto *hsh1_357 = buffer.data(hsh1 + 357);
    const auto *hsh1_358 = buffer.data(hsh1 + 358);
    const auto *hsh1_359 = buffer.data(hsh1 + 359);

    const auto *hsi_385 = buffer.data(hsi + 385);
    const auto *hsi_386 = buffer.data(hsi + 386);
    const auto *hsi_387 = buffer.data(hsi + 387);
    const auto *hsi_388 = buffer.data(hsi + 388);
    const auto *hsi_389 = buffer.data(hsi + 389);
    const auto *hsi_390 = buffer.data(hsi + 390);
    const auto *hsi_391 = buffer.data(hsi + 391);
    const auto *hsi_392 = buffer.data(hsi + 392);
    const auto *hsi_393 = buffer.data(hsi + 393);
    const auto *hsi_394 = buffer.data(hsi + 394);
    const auto *hsi_395 = buffer.data(hsi + 395);
    const auto *hsi_396 = buffer.data(hsi + 396);
    const auto *hsi_397 = buffer.data(hsi + 397);
    const auto *hsi_398 = buffer.data(hsi + 398);
    const auto *hsi_399 = buffer.data(hsi + 399);
    const auto *hsi_400 = buffer.data(hsi + 400);
    const auto *hsi_401 = buffer.data(hsi + 401);
    const auto *hsi_402 = buffer.data(hsi + 402);
    const auto *hsi_403 = buffer.data(hsi + 403);
    const auto *hsi_404 = buffer.data(hsi + 404);
    const auto *hsi_405 = buffer.data(hsi + 405);
    const auto *hsi_406 = buffer.data(hsi + 406);
    const auto *hsi_412 = buffer.data(hsi + 412);
    const auto *hsi_413 = buffer.data(hsi + 413);
    const auto *hsi_414 = buffer.data(hsi + 414);
    const auto *hsi_415 = buffer.data(hsi + 415);
    const auto *hsi_416 = buffer.data(hsi + 416);
    const auto *hsi_417 = buffer.data(hsi + 417);
    const auto *hsi_419 = buffer.data(hsi + 419);
    const auto *hsi_420 = buffer.data(hsi + 420);
    const auto *hsi_421 = buffer.data(hsi + 421);
    const auto *hsi_423 = buffer.data(hsi + 423);
    const auto *hsi_425 = buffer.data(hsi + 425);
    const auto *hsi_426 = buffer.data(hsi + 426);
    const auto *hsi_428 = buffer.data(hsi + 428);
    const auto *hsi_429 = buffer.data(hsi + 429);
    const auto *hsi_430 = buffer.data(hsi + 430);
    const auto *hsi_432 = buffer.data(hsi + 432);
    const auto *hsi_433 = buffer.data(hsi + 433);
    const auto *hsi_434 = buffer.data(hsi + 434);
    const auto *hsi_435 = buffer.data(hsi + 435);
    const auto *hsi_437 = buffer.data(hsi + 437);
    const auto *hsi_438 = buffer.data(hsi + 438);
    const auto *hsi_439 = buffer.data(hsi + 439);
    const auto *hsi_440 = buffer.data(hsi + 440);
    const auto *hsi_441 = buffer.data(hsi + 441);
    const auto *hsi_442 = buffer.data(hsi + 442);
    const auto *hsi_443 = buffer.data(hsi + 443);
    const auto *hsi_444 = buffer.data(hsi + 444);
    const auto *hsi_445 = buffer.data(hsi + 445);
    const auto *hsi_446 = buffer.data(hsi + 446);
    const auto *hsi_447 = buffer.data(hsi + 447);
    const auto *hsi_450 = buffer.data(hsi + 450);
    const auto *hsi_452 = buffer.data(hsi + 452);
    const auto *hsi_453 = buffer.data(hsi + 453);
    const auto *hsi_455 = buffer.data(hsi + 455);
    const auto *hsi_456 = buffer.data(hsi + 456);
    const auto *hsi_457 = buffer.data(hsi + 457);
    const auto *hsi_459 = buffer.data(hsi + 459);
    const auto *hsi_460 = buffer.data(hsi + 460);
    const auto *hsi_461 = buffer.data(hsi + 461);
    const auto *hsi_462 = buffer.data(hsi + 462);
    const auto *hsi_464 = buffer.data(hsi + 464);
    const auto *hsi_465 = buffer.data(hsi + 465);
    const auto *hsi_466 = buffer.data(hsi + 466);
    const auto *hsi_467 = buffer.data(hsi + 467);
    const auto *hsi_468 = buffer.data(hsi + 468);
    const auto *hsi_469 = buffer.data(hsi + 469);
    const auto *hsi_470 = buffer.data(hsi + 470);
    const auto *hsi_471 = buffer.data(hsi + 471);
    const auto *hsi_472 = buffer.data(hsi + 472);
    const auto *hsi_473 = buffer.data(hsi + 473);
    const auto *hsi_474 = buffer.data(hsi + 474);
    const auto *hsi_475 = buffer.data(hsi + 475);
    const auto *hsi_476 = buffer.data(hsi + 476);
    const auto *hsi_477 = buffer.data(hsi + 477);
    const auto *hsi_478 = buffer.data(hsi + 478);

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, gsi_385, gsi_386, gsi_387, \
                         gsi_388, gsi_389, hsi_385, hsi_386, hsi_387, hsi_388, \
                         hsi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_13 * gsi_385[k]
                   + f_3 * pc_x[k] * hsi_385[k];

        t_490[k] = f_13 * gsi_386[k]
                   + f_3 * pc_x[k] * hsi_386[k];

        t_491[k] = f_13 * gsi_387[k]
                   + f_3 * pc_x[k] * hsi_387[k];

        t_492[k] = f_13 * gsi_388[k]
                   + f_3 * pc_x[k] * hsi_388[k];

        t_493[k] = f_13 * gsi_389[k]
                   + f_3 * pc_x[k] * hsi_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pa_x, pc_x, pc_z, gsk0_496, gsi_245, \
                         gsi_390, gsi_391, gsk1_496, hsi_385, hsi_390, \
                         hsi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_13 * gsi_390[k]
                   + f_3 * pc_x[k] * hsi_390[k];

        t_495[k] = f_13 * gsi_391[k]
                   + f_3 * pc_x[k] * hsi_391[k];

        t_496[k] = pa_x[k] * gsk0_496[k]
                   - f_12 * pc_x[k] * gsk1_496[k];

        t_497[k] = f_15 * gsi_245[k]
                   + f_3 * pc_z[k] * hsi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, pa_x, pc_x, gsk0_498, gsk0_499, gsk0_500, \
                         gsk0_501, gsk1_498, gsk1_499, gsk1_500, \
                         gsk1_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = pa_x[k] * gsk0_498[k]
                   - f_12 * pc_x[k] * gsk1_498[k];

        t_499[k] = pa_x[k] * gsk0_499[k]
                   - f_12 * pc_x[k] * gsk1_499[k];

        t_500[k] = pa_x[k] * gsk0_500[k]
                   - f_12 * pc_x[k] * gsk1_500[k];

        t_501[k] = pa_x[k] * gsk0_501[k]
                   - f_12 * pc_x[k] * gsk1_501[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, pa_x, pc_x, pc_y, gsk0_503, gsk0_504, \
                         gsi_279, gsi_392, gsk1_503, gsk1_504, hsi_391, \
                         hsi_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * gsi_279[k]
                   + f_3 * pc_y[k] * hsi_391[k];

        t_503[k] = pa_x[k] * gsk0_503[k]
                   - f_12 * pc_x[k] * gsk1_503[k];

        t_504[k] = pa_x[k] * gsk0_504[k]
                   + f_19 * gsi_392[k]
                   - f_12 * pc_x[k] * gsk1_504[k];

        t_505[k] = f_3 * pc_y[k] * hsi_392[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pc_y, pc_z, gsi_252, hsh0_294, hsh1_294, \
                         hsi_392, hsi_393, hsi_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_16 * gsi_252[k]
                   + f_3 * pc_z[k] * hsi_392[k];

        t_507[k] = f_4 * hsh0_294[k]
                   - f_5 * hsh1_294[k]
                   + f_3 * pc_y[k] * hsi_393[k];

        t_508[k] = f_3 * pc_y[k] * hsi_394[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_x, pc_x, pc_y, gsk0_509, gsi_397, gsk1_509, \
                         hsh0_295, hsh0_296, hsh1_295, hsh1_296, hsi_395, \
                         hsi_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pa_x[k] * gsk0_509[k]
                   + f_0 * gsi_397[k]
                   - f_12 * pc_x[k] * gsk1_509[k];

        t_510[k] = f_6 * hsh0_295[k]
                   - f_7 * hsh1_295[k]
                   + f_3 * pc_y[k] * hsi_395[k];

        t_511[k] = f_4 * hsh0_296[k]
                   - f_5 * hsh1_296[k]
                   + f_3 * pc_y[k] * hsi_396[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pa_x, pc_x, pc_y, gsk0_513, gsi_401, gsk1_513, \
                         hsh0_297, hsh1_297, hsi_397, hsi_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_3 * pc_y[k] * hsi_397[k];

        t_513[k] = pa_x[k] * gsk0_513[k]
                   + f_16 * gsi_401[k]
                   - f_12 * pc_x[k] * gsk1_513[k];

        t_514[k] = f_8 * hsh0_297[k]
                   - f_9 * hsh1_297[k]
                   + f_3 * pc_y[k] * hsi_398[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_y, hsh0_298, hsh0_299, hsh1_298, hsh1_299, \
                         hsi_399, hsi_400, hsi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_6 * hsh0_298[k]
                   - f_7 * hsh1_298[k]
                   + f_3 * pc_y[k] * hsi_399[k];

        t_516[k] = f_4 * hsh0_299[k]
                   - f_5 * hsh1_299[k]
                   + f_3 * pc_y[k] * hsi_400[k];

        t_517[k] = f_3 * pc_y[k] * hsi_401[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pa_x, pc_x, pc_y, gsk0_518, gsi_406, gsk1_518, \
                         hsh0_300, hsh0_301, hsh1_300, hsh1_301, hsi_402, \
                         hsi_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pa_x[k] * gsk0_518[k]
                   + f_15 * gsi_406[k]
                   - f_12 * pc_x[k] * gsk1_518[k];

        t_519[k] = f_10 * hsh0_300[k]
                   - f_11 * hsh1_300[k]
                   + f_3 * pc_y[k] * hsi_402[k];

        t_520[k] = f_8 * hsh0_301[k]
                   - f_9 * hsh1_301[k]
                   + f_3 * pc_y[k] * hsi_403[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, pc_y, hsh0_302, hsh0_303, hsh1_302, hsh1_303, \
                         hsi_404, hsi_405, hsi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_6 * hsh0_302[k]
                   - f_7 * hsh1_302[k]
                   + f_3 * pc_y[k] * hsi_404[k];

        t_522[k] = f_4 * hsh0_303[k]
                   - f_5 * hsh1_303[k]
                   + f_3 * pc_y[k] * hsi_405[k];

        t_523[k] = f_3 * pc_y[k] * hsi_406[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, pa_x, pc_x, gsk0_524, gsi_412, gsi_413, \
                         gsi_414, gsi_415, gsk1_524, hsi_413, hsi_414, \
                         hsi_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = pa_x[k] * gsk0_524[k]
                   + f_14 * gsi_412[k]
                   - f_12 * pc_x[k] * gsk1_524[k];

        t_525[k] = f_13 * gsi_413[k]
                   + f_3 * pc_x[k] * hsi_413[k];

        t_526[k] = f_13 * gsi_414[k]
                   + f_3 * pc_x[k] * hsi_414[k];

        t_527[k] = f_13 * gsi_415[k]
                   + f_3 * pc_x[k] * hsi_415[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, pc_x, pc_y, gsi_416, gsi_417, gsi_419, \
                         hsi_412, hsi_416, hsi_417, hsi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_13 * gsi_416[k]
                   + f_3 * pc_x[k] * hsi_416[k];

        t_529[k] = f_13 * gsi_417[k]
                   + f_3 * pc_x[k] * hsi_417[k];

        t_530[k] = f_3 * pc_y[k] * hsi_412[k];

        t_531[k] = f_13 * gsi_419[k]
                   + f_3 * pc_x[k] * hsi_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pa_x, pc_x, gsk0_532, gsk0_533, gsk0_534, \
                         gsk0_535, gsk1_532, gsk1_533, gsk1_534, \
                         gsk1_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = pa_x[k] * gsk0_532[k]
                   - f_12 * pc_x[k] * gsk1_532[k];

        t_533[k] = pa_x[k] * gsk0_533[k]
                   - f_12 * pc_x[k] * gsk1_533[k];

        t_534[k] = pa_x[k] * gsk0_534[k]
                   - f_12 * pc_x[k] * gsk1_534[k];

        t_535[k] = pa_x[k] * gsk0_535[k]
                   - f_12 * pc_x[k] * gsk1_535[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pa_x, pc_x, pc_y, gsk0_536, gsk0_537, \
                         gsk0_539, gsk1_536, gsk1_537, gsk1_539, \
                         hsi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = pa_x[k] * gsk0_536[k]
                   - f_12 * pc_x[k] * gsk1_536[k];

        t_537[k] = pa_x[k] * gsk0_537[k]
                   - f_12 * pc_x[k] * gsk1_537[k];

        t_538[k] = f_3 * pc_y[k] * hsi_419[k];

        t_539[k] = pa_x[k] * gsk0_539[k]
                   - f_12 * pc_x[k] * gsk1_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, pc_x, pc_z, hsh0_315, hsh0_316, \
                         hsh0_318, hsh1_315, hsh1_316, hsh1_318, hsi_420, hsi_421, \
                         hsi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * hsh0_315[k]
                   - f_2 * hsh1_315[k]
                   + f_3 * pc_x[k] * hsi_420[k];

        t_541[k] = f_17 * hsh0_316[k]
                   - f_18 * hsh1_316[k]
                   + f_3 * pc_x[k] * hsi_421[k];

        t_542[k] = f_3 * pc_z[k] * hsi_420[k];

        t_543[k] = f_10 * hsh0_318[k]
                   - f_11 * hsh1_318[k]
                   + f_3 * pc_x[k] * hsi_423[k];

        t_544[k] = f_3 * pc_z[k] * hsi_421[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pc_x, pc_z, hsh0_320, hsh0_321, hsh0_323, \
                         hsh1_320, hsh1_321, hsh1_323, hsi_423, hsi_425, hsi_426, \
                         hsi_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_10 * hsh0_320[k]
                   - f_11 * hsh1_320[k]
                   + f_3 * pc_x[k] * hsi_425[k];

        t_546[k] = f_8 * hsh0_321[k]
                   - f_9 * hsh1_321[k]
                   + f_3 * pc_x[k] * hsi_426[k];

        t_547[k] = f_3 * pc_z[k] * hsi_423[k];

        t_548[k] = f_8 * hsh0_323[k]
                   - f_9 * hsh1_323[k]
                   + f_3 * pc_x[k] * hsi_428[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, pc_x, pc_z, hsh0_324, hsh0_325, hsh0_327, \
                         hsh1_324, hsh1_325, hsh1_327, hsi_426, hsi_429, hsi_430, \
                         hsi_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_8 * hsh0_324[k]
                   - f_9 * hsh1_324[k]
                   + f_3 * pc_x[k] * hsi_429[k];

        t_550[k] = f_6 * hsh0_325[k]
                   - f_7 * hsh1_325[k]
                   + f_3 * pc_x[k] * hsi_430[k];

        t_551[k] = f_3 * pc_z[k] * hsi_426[k];

        t_552[k] = f_6 * hsh0_327[k]
                   - f_7 * hsh1_327[k]
                   + f_3 * pc_x[k] * hsi_432[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pc_x, pc_z, hsh0_328, hsh0_329, hsh0_330, \
                         hsh1_328, hsh1_329, hsh1_330, hsi_430, hsi_433, hsi_434, \
                         hsi_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_6 * hsh0_328[k]
                   - f_7 * hsh1_328[k]
                   + f_3 * pc_x[k] * hsi_433[k];

        t_554[k] = f_6 * hsh0_329[k]
                   - f_7 * hsh1_329[k]
                   + f_3 * pc_x[k] * hsi_434[k];

        t_555[k] = f_4 * hsh0_330[k]
                   - f_5 * hsh1_330[k]
                   + f_3 * pc_x[k] * hsi_435[k];

        t_556[k] = f_3 * pc_z[k] * hsi_430[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, hsh0_332, hsh0_333, hsh0_334, hsh1_332, \
                         hsh1_333, hsh1_334, hsi_437, hsi_438, \
                         hsi_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_4 * hsh0_332[k]
                   - f_5 * hsh1_332[k]
                   + f_3 * pc_x[k] * hsi_437[k];

        t_558[k] = f_4 * hsh0_333[k]
                   - f_5 * hsh1_333[k]
                   + f_3 * pc_x[k] * hsi_438[k];

        t_559[k] = f_4 * hsh0_334[k]
                   - f_5 * hsh1_334[k]
                   + f_3 * pc_x[k] * hsi_439[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, t_565, pc_x, hsh0_335, hsh1_335, \
                         hsi_440, hsi_441, hsi_442, hsi_443, hsi_444, \
                         hsi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_4 * hsh0_335[k]
                   - f_5 * hsh1_335[k]
                   + f_3 * pc_x[k] * hsi_440[k];

        t_561[k] = f_3 * pc_x[k] * hsi_441[k];

        t_562[k] = f_3 * pc_x[k] * hsi_442[k];

        t_563[k] = f_3 * pc_x[k] * hsi_443[k];

        t_564[k] = f_3 * pc_x[k] * hsi_444[k];

        t_565[k] = f_3 * pc_x[k] * hsi_445[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, gsi_301, \
                         hsh0_330, hsh1_330, hsi_441, hsi_442, hsi_446, \
                         hsi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_3 * pc_x[k] * hsi_446[k];

        t_567[k] = f_3 * pc_x[k] * hsi_447[k];

        t_568[k] = f_0 * gsi_301[k]
                   + f_1 * hsh0_330[k]
                   - f_2 * hsh1_330[k]
                   + f_3 * pc_y[k] * hsi_441[k];

        t_569[k] = f_3 * pc_z[k] * hsi_441[k];

        t_570[k] = f_4 * hsh0_330[k]
                   - f_5 * hsh1_330[k]
                   + f_3 * pc_z[k] * hsi_442[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_z, hsh0_331, hsh0_332, hsh0_333, hsh1_331, \
                         hsh1_332, hsh1_333, hsi_443, hsi_444, \
                         hsi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_6 * hsh0_331[k]
                   - f_7 * hsh1_331[k]
                   + f_3 * pc_z[k] * hsi_443[k];

        t_572[k] = f_8 * hsh0_332[k]
                   - f_9 * hsh1_332[k]
                   + f_3 * pc_z[k] * hsi_444[k];

        t_573[k] = f_10 * hsh0_333[k]
                   - f_11 * hsh1_333[k]
                   + f_3 * pc_z[k] * hsi_445[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_z, pc_y, pc_z, gsk0_360, gsk0_361, \
                         gsi_307, gsk1_360, gsk1_361, hsh0_335, hsh1_335, \
                         hsi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_0 * gsi_307[k]
                   + f_3 * pc_y[k] * hsi_447[k];

        t_575[k] = f_1 * hsh0_335[k]
                   - f_2 * hsh1_335[k]
                   + f_3 * pc_z[k] * hsi_447[k];

        t_576[k] = pa_z[k] * gsk0_360[k]
                   - f_12 * pc_z[k] * gsk1_360[k];

        t_577[k] = pa_z[k] * gsk0_361[k]
                   - f_12 * pc_z[k] * gsk1_361[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pa_z, pc_x, pc_z, gsk0_363, gsk1_363, hsh0_338, \
                         hsh0_340, hsh1_338, hsh1_340, hsi_450, \
                         hsi_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_17 * hsh0_338[k]
                   - f_18 * hsh1_338[k]
                   + f_3 * pc_x[k] * hsi_450[k];

        t_579[k] = pa_z[k] * gsk0_363[k]
                   - f_12 * pc_z[k] * gsk1_363[k];

        t_580[k] = f_10 * hsh0_340[k]
                   - f_11 * hsh1_340[k]
                   + f_3 * pc_x[k] * hsi_452[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pa_z, pc_x, pc_z, gsk0_366, gsk1_366, hsh0_341, \
                         hsh0_343, hsh1_341, hsh1_343, hsi_453, \
                         hsi_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_10 * hsh0_341[k]
                   - f_11 * hsh1_341[k]
                   + f_3 * pc_x[k] * hsi_453[k];

        t_582[k] = pa_z[k] * gsk0_366[k]
                   - f_12 * pc_z[k] * gsk1_366[k];

        t_583[k] = f_8 * hsh0_343[k]
                   - f_9 * hsh1_343[k]
                   + f_3 * pc_x[k] * hsi_455[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pa_z, pc_x, pc_z, gsk0_370, gsk1_370, hsh0_344, \
                         hsh0_345, hsh1_344, hsh1_345, hsi_456, \
                         hsi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_8 * hsh0_344[k]
                   - f_9 * hsh1_344[k]
                   + f_3 * pc_x[k] * hsi_456[k];

        t_585[k] = f_8 * hsh0_345[k]
                   - f_9 * hsh1_345[k]
                   + f_3 * pc_x[k] * hsi_457[k];

        t_586[k] = pa_z[k] * gsk0_370[k]
                   - f_12 * pc_z[k] * gsk1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pc_x, hsh0_347, hsh0_348, hsh0_349, hsh1_347, \
                         hsh1_348, hsh1_349, hsi_459, hsi_460, \
                         hsi_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_6 * hsh0_347[k]
                   - f_7 * hsh1_347[k]
                   + f_3 * pc_x[k] * hsi_459[k];

        t_588[k] = f_6 * hsh0_348[k]
                   - f_7 * hsh1_348[k]
                   + f_3 * pc_x[k] * hsi_460[k];

        t_589[k] = f_6 * hsh0_349[k]
                   - f_7 * hsh1_349[k]
                   + f_3 * pc_x[k] * hsi_461[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pa_z, pc_x, pc_z, gsk0_375, gsk1_375, hsh0_350, \
                         hsh0_352, hsh1_350, hsh1_352, hsi_462, \
                         hsi_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_6 * hsh0_350[k]
                   - f_7 * hsh1_350[k]
                   + f_3 * pc_x[k] * hsi_462[k];

        t_591[k] = pa_z[k] * gsk0_375[k]
                   - f_12 * pc_z[k] * gsk1_375[k];

        t_592[k] = f_4 * hsh0_352[k]
                   - f_5 * hsh1_352[k]
                   + f_3 * pc_x[k] * hsi_464[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pc_x, hsh0_353, hsh0_354, hsh0_355, hsh1_353, \
                         hsh1_354, hsh1_355, hsi_465, hsi_466, \
                         hsi_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_4 * hsh0_353[k]
                   - f_5 * hsh1_353[k]
                   + f_3 * pc_x[k] * hsi_465[k];

        t_594[k] = f_4 * hsh0_354[k]
                   - f_5 * hsh1_354[k]
                   + f_3 * pc_x[k] * hsi_466[k];

        t_595[k] = f_4 * hsh0_355[k]
                   - f_5 * hsh1_355[k]
                   + f_3 * pc_x[k] * hsi_467[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, t_600, t_601, pc_x, hsh0_356, hsh1_356, \
                         hsi_468, hsi_469, hsi_470, hsi_471, hsi_472, \
                         hsi_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_4 * hsh0_356[k]
                   - f_5 * hsh1_356[k]
                   + f_3 * pc_x[k] * hsi_468[k];

        t_597[k] = f_3 * pc_x[k] * hsi_469[k];

        t_598[k] = f_3 * pc_x[k] * hsi_470[k];

        t_599[k] = f_3 * pc_x[k] * hsi_471[k];

        t_600[k] = f_3 * pc_x[k] * hsi_472[k];

        t_601[k] = f_3 * pc_x[k] * hsi_473[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, pa_z, pc_x, pc_z, gsk0_388, gsi_301, \
                         gsk1_388, hsi_469, hsi_474, hsi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_3 * pc_x[k] * hsi_474[k];

        t_603[k] = f_3 * pc_x[k] * hsi_475[k];

        t_604[k] = pa_z[k] * gsk0_388[k]
                   - f_12 * pc_z[k] * gsk1_388[k];

        t_605[k] = f_13 * gsi_301[k]
                   + f_3 * pc_z[k] * hsi_469[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, pa_z, pc_z, gsk0_390, gsk0_391, gsk0_392, \
                         gsi_302, gsi_303, gsi_304, gsk1_390, gsk1_391, \
                         gsk1_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = pa_z[k] * gsk0_390[k]
                   + f_14 * gsi_302[k]
                   - f_12 * pc_z[k] * gsk1_390[k];

        t_607[k] = pa_z[k] * gsk0_391[k]
                   + f_15 * gsi_303[k]
                   - f_12 * pc_z[k] * gsk1_391[k];

        t_608[k] = pa_z[k] * gsk0_392[k]
                   + f_16 * gsi_304[k]
                   - f_12 * pc_z[k] * gsk1_392[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pa_z, pc_y, pc_z, gsk0_393, gsi_305, gsi_307, \
                         gsi_335, gsk1_393, hsh0_356, hsh1_356, \
                         hsi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = pa_z[k] * gsk0_393[k]
                   + f_0 * gsi_305[k]
                   - f_12 * pc_z[k] * gsk1_393[k];

        t_610[k] = f_16 * gsi_335[k]
                   + f_3 * pc_y[k] * hsi_475[k];

        t_611[k] = f_13 * gsi_307[k]
                   + f_1 * hsh0_356[k]
                   - f_2 * hsh1_356[k]
                   + f_3 * pc_z[k] * hsi_475[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pc_x, hsh0_357, hsh0_358, hsh0_359, hsh1_357, \
                         hsh1_358, hsh1_359, hsi_476, hsi_477, \
                         hsi_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_1 * hsh0_357[k]
                   - f_2 * hsh1_357[k]
                   + f_3 * pc_x[k] * hsi_476[k];

        t_613[k] = f_17 * hsh0_358[k]
                   - f_18 * hsh1_358[k]
                   + f_3 * pc_x[k] * hsi_477[k];

        t_614[k] = f_17 * hsh0_359[k]
                   - f_18 * hsh1_359[k]
                   + f_3 * pc_x[k] * hsi_478[k];
    }
}

static auto
compute_prim_hsk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsk0,
                                                          const size_t gsi, const size_t gsk1,
                                                          const size_t hsh0, const size_t hsh1,
                                                          const size_t hsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);
    const auto f_19 = 3.5 / q;

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
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsk0_504 = buffer.data(gsk0 + 504);
    const auto *gsk0_506 = buffer.data(gsk0 + 506);
    const auto *gsk0_509 = buffer.data(gsk0 + 509);
    const auto *gsk0_513 = buffer.data(gsk0 + 513);
    const auto *gsk0_518 = buffer.data(gsk0 + 518);
    const auto *gsk0_524 = buffer.data(gsk0 + 524);
    const auto *gsk0_532 = buffer.data(gsk0 + 532);
    const auto *gsk0_534 = buffer.data(gsk0 + 534);
    const auto *gsk0_535 = buffer.data(gsk0 + 535);
    const auto *gsk0_536 = buffer.data(gsk0 + 536);
    const auto *gsk0_537 = buffer.data(gsk0 + 537);
    const auto *gsk0_539 = buffer.data(gsk0 + 539);

    const auto *gsi_329 = buffer.data(gsi + 329);
    const auto *gsi_335 = buffer.data(gsi + 335);
    const auto *gsi_357 = buffer.data(gsi + 357);
    const auto *gsi_359 = buffer.data(gsi + 359);
    const auto *gsi_360 = buffer.data(gsi + 360);
    const auto *gsi_361 = buffer.data(gsi + 361);
    const auto *gsi_362 = buffer.data(gsi + 362);
    const auto *gsi_363 = buffer.data(gsi + 363);
    const auto *gsi_385 = buffer.data(gsi + 385);
    const auto *gsi_387 = buffer.data(gsi + 387);
    const auto *gsi_388 = buffer.data(gsi + 388);
    const auto *gsi_389 = buffer.data(gsi + 389);
    const auto *gsi_390 = buffer.data(gsi + 390);
    const auto *gsi_391 = buffer.data(gsi + 391);
    const auto *gsi_413 = buffer.data(gsi + 413);
    const auto *gsi_415 = buffer.data(gsi + 415);
    const auto *gsi_416 = buffer.data(gsi + 416);
    const auto *gsi_417 = buffer.data(gsi + 417);
    const auto *gsi_418 = buffer.data(gsi + 418);
    const auto *gsi_419 = buffer.data(gsi + 419);

    const auto *gsk1_504 = buffer.data(gsk1 + 504);
    const auto *gsk1_506 = buffer.data(gsk1 + 506);
    const auto *gsk1_509 = buffer.data(gsk1 + 509);
    const auto *gsk1_513 = buffer.data(gsk1 + 513);
    const auto *gsk1_518 = buffer.data(gsk1 + 518);
    const auto *gsk1_524 = buffer.data(gsk1 + 524);
    const auto *gsk1_532 = buffer.data(gsk1 + 532);
    const auto *gsk1_534 = buffer.data(gsk1 + 534);
    const auto *gsk1_535 = buffer.data(gsk1 + 535);
    const auto *gsk1_536 = buffer.data(gsk1 + 536);
    const auto *gsk1_537 = buffer.data(gsk1 + 537);
    const auto *gsk1_539 = buffer.data(gsk1 + 539);

    const auto *hsh0_360 = buffer.data(hsh0 + 360);
    const auto *hsh0_361 = buffer.data(hsh0 + 361);
    const auto *hsh0_362 = buffer.data(hsh0 + 362);
    const auto *hsh0_363 = buffer.data(hsh0 + 363);
    const auto *hsh0_364 = buffer.data(hsh0 + 364);
    const auto *hsh0_365 = buffer.data(hsh0 + 365);
    const auto *hsh0_366 = buffer.data(hsh0 + 366);
    const auto *hsh0_367 = buffer.data(hsh0 + 367);
    const auto *hsh0_368 = buffer.data(hsh0 + 368);
    const auto *hsh0_369 = buffer.data(hsh0 + 369);
    const auto *hsh0_370 = buffer.data(hsh0 + 370);
    const auto *hsh0_371 = buffer.data(hsh0 + 371);
    const auto *hsh0_372 = buffer.data(hsh0 + 372);
    const auto *hsh0_373 = buffer.data(hsh0 + 373);
    const auto *hsh0_374 = buffer.data(hsh0 + 374);
    const auto *hsh0_375 = buffer.data(hsh0 + 375);
    const auto *hsh0_376 = buffer.data(hsh0 + 376);
    const auto *hsh0_377 = buffer.data(hsh0 + 377);
    const auto *hsh0_378 = buffer.data(hsh0 + 378);
    const auto *hsh0_379 = buffer.data(hsh0 + 379);
    const auto *hsh0_380 = buffer.data(hsh0 + 380);
    const auto *hsh0_381 = buffer.data(hsh0 + 381);
    const auto *hsh0_382 = buffer.data(hsh0 + 382);
    const auto *hsh0_383 = buffer.data(hsh0 + 383);
    const auto *hsh0_384 = buffer.data(hsh0 + 384);
    const auto *hsh0_385 = buffer.data(hsh0 + 385);
    const auto *hsh0_386 = buffer.data(hsh0 + 386);
    const auto *hsh0_387 = buffer.data(hsh0 + 387);
    const auto *hsh0_388 = buffer.data(hsh0 + 388);
    const auto *hsh0_389 = buffer.data(hsh0 + 389);
    const auto *hsh0_390 = buffer.data(hsh0 + 390);
    const auto *hsh0_391 = buffer.data(hsh0 + 391);
    const auto *hsh0_392 = buffer.data(hsh0 + 392);
    const auto *hsh0_393 = buffer.data(hsh0 + 393);
    const auto *hsh0_394 = buffer.data(hsh0 + 394);
    const auto *hsh0_395 = buffer.data(hsh0 + 395);
    const auto *hsh0_396 = buffer.data(hsh0 + 396);
    const auto *hsh0_397 = buffer.data(hsh0 + 397);
    const auto *hsh0_398 = buffer.data(hsh0 + 398);
    const auto *hsh0_400 = buffer.data(hsh0 + 400);
    const auto *hsh0_402 = buffer.data(hsh0 + 402);
    const auto *hsh0_403 = buffer.data(hsh0 + 403);
    const auto *hsh0_405 = buffer.data(hsh0 + 405);
    const auto *hsh0_406 = buffer.data(hsh0 + 406);
    const auto *hsh0_407 = buffer.data(hsh0 + 407);
    const auto *hsh0_409 = buffer.data(hsh0 + 409);
    const auto *hsh0_410 = buffer.data(hsh0 + 410);
    const auto *hsh0_411 = buffer.data(hsh0 + 411);
    const auto *hsh0_412 = buffer.data(hsh0 + 412);
    const auto *hsh0_414 = buffer.data(hsh0 + 414);
    const auto *hsh0_415 = buffer.data(hsh0 + 415);
    const auto *hsh0_416 = buffer.data(hsh0 + 416);
    const auto *hsh0_417 = buffer.data(hsh0 + 417);
    const auto *hsh0_418 = buffer.data(hsh0 + 418);
    const auto *hsh0_420 = buffer.data(hsh0 + 420);
    const auto *hsh0_422 = buffer.data(hsh0 + 422);
    const auto *hsh0_423 = buffer.data(hsh0 + 423);
    const auto *hsh0_425 = buffer.data(hsh0 + 425);
    const auto *hsh0_426 = buffer.data(hsh0 + 426);
    const auto *hsh0_427 = buffer.data(hsh0 + 427);
    const auto *hsh0_429 = buffer.data(hsh0 + 429);
    const auto *hsh0_430 = buffer.data(hsh0 + 430);
    const auto *hsh0_431 = buffer.data(hsh0 + 431);

    const auto *hsh1_360 = buffer.data(hsh1 + 360);
    const auto *hsh1_361 = buffer.data(hsh1 + 361);
    const auto *hsh1_362 = buffer.data(hsh1 + 362);
    const auto *hsh1_363 = buffer.data(hsh1 + 363);
    const auto *hsh1_364 = buffer.data(hsh1 + 364);
    const auto *hsh1_365 = buffer.data(hsh1 + 365);
    const auto *hsh1_366 = buffer.data(hsh1 + 366);
    const auto *hsh1_367 = buffer.data(hsh1 + 367);
    const auto *hsh1_368 = buffer.data(hsh1 + 368);
    const auto *hsh1_369 = buffer.data(hsh1 + 369);
    const auto *hsh1_370 = buffer.data(hsh1 + 370);
    const auto *hsh1_371 = buffer.data(hsh1 + 371);
    const auto *hsh1_372 = buffer.data(hsh1 + 372);
    const auto *hsh1_373 = buffer.data(hsh1 + 373);
    const auto *hsh1_374 = buffer.data(hsh1 + 374);
    const auto *hsh1_375 = buffer.data(hsh1 + 375);
    const auto *hsh1_376 = buffer.data(hsh1 + 376);
    const auto *hsh1_377 = buffer.data(hsh1 + 377);
    const auto *hsh1_378 = buffer.data(hsh1 + 378);
    const auto *hsh1_379 = buffer.data(hsh1 + 379);
    const auto *hsh1_380 = buffer.data(hsh1 + 380);
    const auto *hsh1_381 = buffer.data(hsh1 + 381);
    const auto *hsh1_382 = buffer.data(hsh1 + 382);
    const auto *hsh1_383 = buffer.data(hsh1 + 383);
    const auto *hsh1_384 = buffer.data(hsh1 + 384);
    const auto *hsh1_385 = buffer.data(hsh1 + 385);
    const auto *hsh1_386 = buffer.data(hsh1 + 386);
    const auto *hsh1_387 = buffer.data(hsh1 + 387);
    const auto *hsh1_388 = buffer.data(hsh1 + 388);
    const auto *hsh1_389 = buffer.data(hsh1 + 389);
    const auto *hsh1_390 = buffer.data(hsh1 + 390);
    const auto *hsh1_391 = buffer.data(hsh1 + 391);
    const auto *hsh1_392 = buffer.data(hsh1 + 392);
    const auto *hsh1_393 = buffer.data(hsh1 + 393);
    const auto *hsh1_394 = buffer.data(hsh1 + 394);
    const auto *hsh1_395 = buffer.data(hsh1 + 395);
    const auto *hsh1_396 = buffer.data(hsh1 + 396);
    const auto *hsh1_397 = buffer.data(hsh1 + 397);
    const auto *hsh1_398 = buffer.data(hsh1 + 398);
    const auto *hsh1_400 = buffer.data(hsh1 + 400);
    const auto *hsh1_402 = buffer.data(hsh1 + 402);
    const auto *hsh1_403 = buffer.data(hsh1 + 403);
    const auto *hsh1_405 = buffer.data(hsh1 + 405);
    const auto *hsh1_406 = buffer.data(hsh1 + 406);
    const auto *hsh1_407 = buffer.data(hsh1 + 407);
    const auto *hsh1_409 = buffer.data(hsh1 + 409);
    const auto *hsh1_410 = buffer.data(hsh1 + 410);
    const auto *hsh1_411 = buffer.data(hsh1 + 411);
    const auto *hsh1_412 = buffer.data(hsh1 + 412);
    const auto *hsh1_414 = buffer.data(hsh1 + 414);
    const auto *hsh1_415 = buffer.data(hsh1 + 415);
    const auto *hsh1_416 = buffer.data(hsh1 + 416);
    const auto *hsh1_417 = buffer.data(hsh1 + 417);
    const auto *hsh1_418 = buffer.data(hsh1 + 418);
    const auto *hsh1_420 = buffer.data(hsh1 + 420);
    const auto *hsh1_422 = buffer.data(hsh1 + 422);
    const auto *hsh1_423 = buffer.data(hsh1 + 423);
    const auto *hsh1_425 = buffer.data(hsh1 + 425);
    const auto *hsh1_426 = buffer.data(hsh1 + 426);
    const auto *hsh1_427 = buffer.data(hsh1 + 427);
    const auto *hsh1_429 = buffer.data(hsh1 + 429);
    const auto *hsh1_430 = buffer.data(hsh1 + 430);
    const auto *hsh1_431 = buffer.data(hsh1 + 431);

    const auto *hsi_479 = buffer.data(hsi + 479);
    const auto *hsi_480 = buffer.data(hsi + 480);
    const auto *hsi_481 = buffer.data(hsi + 481);
    const auto *hsi_482 = buffer.data(hsi + 482);
    const auto *hsi_483 = buffer.data(hsi + 483);
    const auto *hsi_484 = buffer.data(hsi + 484);
    const auto *hsi_485 = buffer.data(hsi + 485);
    const auto *hsi_486 = buffer.data(hsi + 486);
    const auto *hsi_487 = buffer.data(hsi + 487);
    const auto *hsi_488 = buffer.data(hsi + 488);
    const auto *hsi_489 = buffer.data(hsi + 489);
    const auto *hsi_490 = buffer.data(hsi + 490);
    const auto *hsi_491 = buffer.data(hsi + 491);
    const auto *hsi_492 = buffer.data(hsi + 492);
    const auto *hsi_493 = buffer.data(hsi + 493);
    const auto *hsi_494 = buffer.data(hsi + 494);
    const auto *hsi_495 = buffer.data(hsi + 495);
    const auto *hsi_496 = buffer.data(hsi + 496);
    const auto *hsi_497 = buffer.data(hsi + 497);
    const auto *hsi_498 = buffer.data(hsi + 498);
    const auto *hsi_499 = buffer.data(hsi + 499);
    const auto *hsi_500 = buffer.data(hsi + 500);
    const auto *hsi_501 = buffer.data(hsi + 501);
    const auto *hsi_502 = buffer.data(hsi + 502);
    const auto *hsi_503 = buffer.data(hsi + 503);
    const auto *hsi_504 = buffer.data(hsi + 504);
    const auto *hsi_505 = buffer.data(hsi + 505);
    const auto *hsi_506 = buffer.data(hsi + 506);
    const auto *hsi_507 = buffer.data(hsi + 507);
    const auto *hsi_508 = buffer.data(hsi + 508);
    const auto *hsi_509 = buffer.data(hsi + 509);
    const auto *hsi_510 = buffer.data(hsi + 510);
    const auto *hsi_511 = buffer.data(hsi + 511);
    const auto *hsi_512 = buffer.data(hsi + 512);
    const auto *hsi_513 = buffer.data(hsi + 513);
    const auto *hsi_514 = buffer.data(hsi + 514);
    const auto *hsi_515 = buffer.data(hsi + 515);
    const auto *hsi_516 = buffer.data(hsi + 516);
    const auto *hsi_517 = buffer.data(hsi + 517);
    const auto *hsi_518 = buffer.data(hsi + 518);
    const auto *hsi_519 = buffer.data(hsi + 519);
    const auto *hsi_520 = buffer.data(hsi + 520);
    const auto *hsi_521 = buffer.data(hsi + 521);
    const auto *hsi_522 = buffer.data(hsi + 522);
    const auto *hsi_523 = buffer.data(hsi + 523);
    const auto *hsi_524 = buffer.data(hsi + 524);
    const auto *hsi_525 = buffer.data(hsi + 525);
    const auto *hsi_526 = buffer.data(hsi + 526);
    const auto *hsi_527 = buffer.data(hsi + 527);
    const auto *hsi_528 = buffer.data(hsi + 528);
    const auto *hsi_529 = buffer.data(hsi + 529);
    const auto *hsi_530 = buffer.data(hsi + 530);
    const auto *hsi_531 = buffer.data(hsi + 531);
    const auto *hsi_533 = buffer.data(hsi + 533);
    const auto *hsi_535 = buffer.data(hsi + 535);
    const auto *hsi_536 = buffer.data(hsi + 536);
    const auto *hsi_538 = buffer.data(hsi + 538);
    const auto *hsi_539 = buffer.data(hsi + 539);
    const auto *hsi_540 = buffer.data(hsi + 540);
    const auto *hsi_542 = buffer.data(hsi + 542);
    const auto *hsi_543 = buffer.data(hsi + 543);
    const auto *hsi_544 = buffer.data(hsi + 544);
    const auto *hsi_545 = buffer.data(hsi + 545);
    const auto *hsi_547 = buffer.data(hsi + 547);
    const auto *hsi_548 = buffer.data(hsi + 548);
    const auto *hsi_549 = buffer.data(hsi + 549);
    const auto *hsi_550 = buffer.data(hsi + 550);
    const auto *hsi_551 = buffer.data(hsi + 551);
    const auto *hsi_553 = buffer.data(hsi + 553);
    const auto *hsi_554 = buffer.data(hsi + 554);
    const auto *hsi_555 = buffer.data(hsi + 555);
    const auto *hsi_556 = buffer.data(hsi + 556);
    const auto *hsi_557 = buffer.data(hsi + 557);
    const auto *hsi_558 = buffer.data(hsi + 558);
    const auto *hsi_559 = buffer.data(hsi + 559);
    const auto *hsi_560 = buffer.data(hsi + 560);
    const auto *hsi_562 = buffer.data(hsi + 562);
    const auto *hsi_563 = buffer.data(hsi + 563);
    const auto *hsi_565 = buffer.data(hsi + 565);
    const auto *hsi_566 = buffer.data(hsi + 566);
    const auto *hsi_567 = buffer.data(hsi + 567);
    const auto *hsi_569 = buffer.data(hsi + 569);
    const auto *hsi_570 = buffer.data(hsi + 570);
    const auto *hsi_571 = buffer.data(hsi + 571);

#pragma omp simd aligned(t_615, t_616, t_617, pc_x, hsh0_360, hsh0_361, hsh0_362, hsh1_360, \
                         hsh1_361, hsh1_362, hsi_479, hsi_480, \
                         hsi_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_10 * hsh0_360[k]
                   - f_11 * hsh1_360[k]
                   + f_3 * pc_x[k] * hsi_479[k];

        t_616[k] = f_10 * hsh0_361[k]
                   - f_11 * hsh1_361[k]
                   + f_3 * pc_x[k] * hsi_480[k];

        t_617[k] = f_10 * hsh0_362[k]
                   - f_11 * hsh1_362[k]
                   + f_3 * pc_x[k] * hsi_481[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, pc_x, hsh0_363, hsh0_364, hsh0_365, hsh1_363, \
                         hsh1_364, hsh1_365, hsi_482, hsi_483, \
                         hsi_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_8 * hsh0_363[k]
                   - f_9 * hsh1_363[k]
                   + f_3 * pc_x[k] * hsi_482[k];

        t_619[k] = f_8 * hsh0_364[k]
                   - f_9 * hsh1_364[k]
                   + f_3 * pc_x[k] * hsi_483[k];

        t_620[k] = f_8 * hsh0_365[k]
                   - f_9 * hsh1_365[k]
                   + f_3 * pc_x[k] * hsi_484[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pc_x, hsh0_366, hsh0_367, hsh0_368, hsh1_366, \
                         hsh1_367, hsh1_368, hsi_485, hsi_486, \
                         hsi_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_8 * hsh0_366[k]
                   - f_9 * hsh1_366[k]
                   + f_3 * pc_x[k] * hsi_485[k];

        t_622[k] = f_6 * hsh0_367[k]
                   - f_7 * hsh1_367[k]
                   + f_3 * pc_x[k] * hsi_486[k];

        t_623[k] = f_6 * hsh0_368[k]
                   - f_7 * hsh1_368[k]
                   + f_3 * pc_x[k] * hsi_487[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pc_x, hsh0_369, hsh0_370, hsh0_371, hsh1_369, \
                         hsh1_370, hsh1_371, hsi_488, hsi_489, \
                         hsi_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_6 * hsh0_369[k]
                   - f_7 * hsh1_369[k]
                   + f_3 * pc_x[k] * hsi_488[k];

        t_625[k] = f_6 * hsh0_370[k]
                   - f_7 * hsh1_370[k]
                   + f_3 * pc_x[k] * hsi_489[k];

        t_626[k] = f_6 * hsh0_371[k]
                   - f_7 * hsh1_371[k]
                   + f_3 * pc_x[k] * hsi_490[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pc_x, hsh0_372, hsh0_373, hsh0_374, hsh1_372, \
                         hsh1_373, hsh1_374, hsi_491, hsi_492, \
                         hsi_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_4 * hsh0_372[k]
                   - f_5 * hsh1_372[k]
                   + f_3 * pc_x[k] * hsi_491[k];

        t_628[k] = f_4 * hsh0_373[k]
                   - f_5 * hsh1_373[k]
                   + f_3 * pc_x[k] * hsi_492[k];

        t_629[k] = f_4 * hsh0_374[k]
                   - f_5 * hsh1_374[k]
                   + f_3 * pc_x[k] * hsi_493[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, hsh0_375, hsh0_376, hsh0_377, \
                         hsh1_375, hsh1_376, hsh1_377, hsi_494, hsi_495, hsi_496, \
                         hsi_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_4 * hsh0_375[k]
                   - f_5 * hsh1_375[k]
                   + f_3 * pc_x[k] * hsi_494[k];

        t_631[k] = f_4 * hsh0_376[k]
                   - f_5 * hsh1_376[k]
                   + f_3 * pc_x[k] * hsi_495[k];

        t_632[k] = f_4 * hsh0_377[k]
                   - f_5 * hsh1_377[k]
                   + f_3 * pc_x[k] * hsi_496[k];

        t_633[k] = f_3 * pc_x[k] * hsi_497[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, t_639, pc_x, hsi_498, hsi_499, \
                         hsi_500, hsi_501, hsi_502, hsi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_3 * pc_x[k] * hsi_498[k];

        t_635[k] = f_3 * pc_x[k] * hsi_499[k];

        t_636[k] = f_3 * pc_x[k] * hsi_500[k];

        t_637[k] = f_3 * pc_x[k] * hsi_501[k];

        t_638[k] = f_3 * pc_x[k] * hsi_502[k];

        t_639[k] = f_3 * pc_x[k] * hsi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, gsi_329, gsi_357, gsi_359, hsh0_372, \
                         hsh0_374, hsh1_372, hsh1_374, hsi_497, \
                         hsi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * gsi_357[k]
                   + f_1 * hsh0_372[k]
                   - f_2 * hsh1_372[k]
                   + f_3 * pc_y[k] * hsi_497[k];

        t_641[k] = f_14 * gsi_329[k]
                   + f_3 * pc_z[k] * hsi_497[k];

        t_642[k] = f_15 * gsi_359[k]
                   + f_10 * hsh0_374[k]
                   - f_11 * hsh1_374[k]
                   + f_3 * pc_y[k] * hsi_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, gsi_360, gsi_361, gsi_362, hsh0_375, \
                         hsh0_376, hsh0_377, hsh1_375, hsh1_376, hsh1_377, hsi_500, hsi_501, \
                         hsi_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * gsi_360[k]
                   + f_8 * hsh0_375[k]
                   - f_9 * hsh1_375[k]
                   + f_3 * pc_y[k] * hsi_500[k];

        t_644[k] = f_15 * gsi_361[k]
                   + f_6 * hsh0_376[k]
                   - f_7 * hsh1_376[k]
                   + f_3 * pc_y[k] * hsi_501[k];

        t_645[k] = f_15 * gsi_362[k]
                   + f_4 * hsh0_377[k]
                   - f_5 * hsh1_377[k]
                   + f_3 * pc_y[k] * hsi_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, gsi_335, gsi_363, hsh0_377, \
                         hsh0_378, hsh1_377, hsh1_378, hsi_503, \
                         hsi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * gsi_363[k]
                   + f_3 * pc_y[k] * hsi_503[k];

        t_647[k] = f_14 * gsi_335[k]
                   + f_1 * hsh0_377[k]
                   - f_2 * hsh1_377[k]
                   + f_3 * pc_z[k] * hsi_503[k];

        t_648[k] = f_1 * hsh0_378[k]
                   - f_2 * hsh1_378[k]
                   + f_3 * pc_x[k] * hsi_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pc_x, hsh0_379, hsh0_380, hsh0_381, hsh1_379, \
                         hsh1_380, hsh1_381, hsi_505, hsi_506, \
                         hsi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_17 * hsh0_379[k]
                   - f_18 * hsh1_379[k]
                   + f_3 * pc_x[k] * hsi_505[k];

        t_650[k] = f_17 * hsh0_380[k]
                   - f_18 * hsh1_380[k]
                   + f_3 * pc_x[k] * hsi_506[k];

        t_651[k] = f_10 * hsh0_381[k]
                   - f_11 * hsh1_381[k]
                   + f_3 * pc_x[k] * hsi_507[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, pc_x, hsh0_382, hsh0_383, hsh0_384, hsh1_382, \
                         hsh1_383, hsh1_384, hsi_508, hsi_509, \
                         hsi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_10 * hsh0_382[k]
                   - f_11 * hsh1_382[k]
                   + f_3 * pc_x[k] * hsi_508[k];

        t_653[k] = f_10 * hsh0_383[k]
                   - f_11 * hsh1_383[k]
                   + f_3 * pc_x[k] * hsi_509[k];

        t_654[k] = f_8 * hsh0_384[k]
                   - f_9 * hsh1_384[k]
                   + f_3 * pc_x[k] * hsi_510[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pc_x, hsh0_385, hsh0_386, hsh0_387, hsh1_385, \
                         hsh1_386, hsh1_387, hsi_511, hsi_512, \
                         hsi_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_8 * hsh0_385[k]
                   - f_9 * hsh1_385[k]
                   + f_3 * pc_x[k] * hsi_511[k];

        t_656[k] = f_8 * hsh0_386[k]
                   - f_9 * hsh1_386[k]
                   + f_3 * pc_x[k] * hsi_512[k];

        t_657[k] = f_8 * hsh0_387[k]
                   - f_9 * hsh1_387[k]
                   + f_3 * pc_x[k] * hsi_513[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, pc_x, hsh0_388, hsh0_389, hsh0_390, hsh1_388, \
                         hsh1_389, hsh1_390, hsi_514, hsi_515, \
                         hsi_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_6 * hsh0_388[k]
                   - f_7 * hsh1_388[k]
                   + f_3 * pc_x[k] * hsi_514[k];

        t_659[k] = f_6 * hsh0_389[k]
                   - f_7 * hsh1_389[k]
                   + f_3 * pc_x[k] * hsi_515[k];

        t_660[k] = f_6 * hsh0_390[k]
                   - f_7 * hsh1_390[k]
                   + f_3 * pc_x[k] * hsi_516[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, pc_x, hsh0_391, hsh0_392, hsh0_393, hsh1_391, \
                         hsh1_392, hsh1_393, hsi_517, hsi_518, \
                         hsi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_6 * hsh0_391[k]
                   - f_7 * hsh1_391[k]
                   + f_3 * pc_x[k] * hsi_517[k];

        t_662[k] = f_6 * hsh0_392[k]
                   - f_7 * hsh1_392[k]
                   + f_3 * pc_x[k] * hsi_518[k];

        t_663[k] = f_4 * hsh0_393[k]
                   - f_5 * hsh1_393[k]
                   + f_3 * pc_x[k] * hsi_519[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, hsh0_394, hsh0_395, hsh0_396, hsh1_394, \
                         hsh1_395, hsh1_396, hsi_520, hsi_521, \
                         hsi_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_4 * hsh0_394[k]
                   - f_5 * hsh1_394[k]
                   + f_3 * pc_x[k] * hsi_520[k];

        t_665[k] = f_4 * hsh0_395[k]
                   - f_5 * hsh1_395[k]
                   + f_3 * pc_x[k] * hsi_521[k];

        t_666[k] = f_4 * hsh0_396[k]
                   - f_5 * hsh1_396[k]
                   + f_3 * pc_x[k] * hsi_522[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, pc_x, hsh0_397, hsh0_398, \
                         hsh1_397, hsh1_398, hsi_523, hsi_524, hsi_525, hsi_526, \
                         hsi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_4 * hsh0_397[k]
                   - f_5 * hsh1_397[k]
                   + f_3 * pc_x[k] * hsi_523[k];

        t_668[k] = f_4 * hsh0_398[k]
                   - f_5 * hsh1_398[k]
                   + f_3 * pc_x[k] * hsi_524[k];

        t_669[k] = f_3 * pc_x[k] * hsi_525[k];

        t_670[k] = f_3 * pc_x[k] * hsi_526[k];

        t_671[k] = f_3 * pc_x[k] * hsi_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, pc_x, pc_y, gsi_385, hsh0_393, \
                         hsh1_393, hsi_525, hsi_528, hsi_529, hsi_530, \
                         hsi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_3 * pc_x[k] * hsi_528[k];

        t_673[k] = f_3 * pc_x[k] * hsi_529[k];

        t_674[k] = f_3 * pc_x[k] * hsi_530[k];

        t_675[k] = f_3 * pc_x[k] * hsi_531[k];

        t_676[k] = f_14 * gsi_385[k]
                   + f_1 * hsh0_393[k]
                   - f_2 * hsh1_393[k]
                   + f_3 * pc_y[k] * hsi_525[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_y, pc_z, gsi_357, gsi_387, gsi_388, hsh0_395, \
                         hsh0_396, hsh1_395, hsh1_396, hsi_525, hsi_527, \
                         hsi_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_15 * gsi_357[k]
                   + f_3 * pc_z[k] * hsi_525[k];

        t_678[k] = f_14 * gsi_387[k]
                   + f_10 * hsh0_395[k]
                   - f_11 * hsh1_395[k]
                   + f_3 * pc_y[k] * hsi_527[k];

        t_679[k] = f_14 * gsi_388[k]
                   + f_8 * hsh0_396[k]
                   - f_9 * hsh1_396[k]
                   + f_3 * pc_y[k] * hsi_528[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_y, gsi_389, gsi_390, gsi_391, hsh0_397, \
                         hsh0_398, hsh1_397, hsh1_398, hsi_529, hsi_530, \
                         hsi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_14 * gsi_389[k]
                   + f_6 * hsh0_397[k]
                   - f_7 * hsh1_397[k]
                   + f_3 * pc_y[k] * hsi_529[k];

        t_681[k] = f_14 * gsi_390[k]
                   + f_4 * hsh0_398[k]
                   - f_5 * hsh1_398[k]
                   + f_3 * pc_y[k] * hsi_530[k];

        t_682[k] = f_14 * gsi_391[k]
                   + f_3 * pc_y[k] * hsi_531[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pa_y, pc_x, pc_y, pc_z, gsk0_504, gsi_363, \
                         gsk1_504, hsh0_398, hsh0_400, hsh1_398, hsh1_400, hsi_531, \
                         hsi_533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_15 * gsi_363[k]
                   + f_1 * hsh0_398[k]
                   - f_2 * hsh1_398[k]
                   + f_3 * pc_z[k] * hsi_531[k];

        t_684[k] = pa_y[k] * gsk0_504[k]
                   - f_12 * pc_y[k] * gsk1_504[k];

        t_685[k] = f_17 * hsh0_400[k]
                   - f_18 * hsh1_400[k]
                   + f_3 * pc_x[k] * hsi_533[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, pa_y, pc_x, pc_y, gsk0_506, gsk1_506, hsh0_402, \
                         hsh0_403, hsh1_402, hsh1_403, hsi_535, \
                         hsi_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = pa_y[k] * gsk0_506[k]
                   - f_12 * pc_y[k] * gsk1_506[k];

        t_687[k] = f_10 * hsh0_402[k]
                   - f_11 * hsh1_402[k]
                   + f_3 * pc_x[k] * hsi_535[k];

        t_688[k] = f_10 * hsh0_403[k]
                   - f_11 * hsh1_403[k]
                   + f_3 * pc_x[k] * hsi_536[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, pa_y, pc_x, pc_y, gsk0_509, gsk1_509, hsh0_405, \
                         hsh0_406, hsh1_405, hsh1_406, hsi_538, \
                         hsi_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = pa_y[k] * gsk0_509[k]
                   - f_12 * pc_y[k] * gsk1_509[k];

        t_690[k] = f_8 * hsh0_405[k]
                   - f_9 * hsh1_405[k]
                   + f_3 * pc_x[k] * hsi_538[k];

        t_691[k] = f_8 * hsh0_406[k]
                   - f_9 * hsh1_406[k]
                   + f_3 * pc_x[k] * hsi_539[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pa_y, pc_x, pc_y, gsk0_513, gsk1_513, hsh0_407, \
                         hsh0_409, hsh1_407, hsh1_409, hsi_540, \
                         hsi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_8 * hsh0_407[k]
                   - f_9 * hsh1_407[k]
                   + f_3 * pc_x[k] * hsi_540[k];

        t_693[k] = pa_y[k] * gsk0_513[k]
                   - f_12 * pc_y[k] * gsk1_513[k];

        t_694[k] = f_6 * hsh0_409[k]
                   - f_7 * hsh1_409[k]
                   + f_3 * pc_x[k] * hsi_542[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, pc_x, hsh0_410, hsh0_411, hsh0_412, hsh1_410, \
                         hsh1_411, hsh1_412, hsi_543, hsi_544, \
                         hsi_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_6 * hsh0_410[k]
                   - f_7 * hsh1_410[k]
                   + f_3 * pc_x[k] * hsi_543[k];

        t_696[k] = f_6 * hsh0_411[k]
                   - f_7 * hsh1_411[k]
                   + f_3 * pc_x[k] * hsi_544[k];

        t_697[k] = f_6 * hsh0_412[k]
                   - f_7 * hsh1_412[k]
                   + f_3 * pc_x[k] * hsi_545[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pa_y, pc_x, pc_y, gsk0_518, gsk1_518, hsh0_414, \
                         hsh0_415, hsh1_414, hsh1_415, hsi_547, \
                         hsi_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = pa_y[k] * gsk0_518[k]
                   - f_12 * pc_y[k] * gsk1_518[k];

        t_699[k] = f_4 * hsh0_414[k]
                   - f_5 * hsh1_414[k]
                   + f_3 * pc_x[k] * hsi_547[k];

        t_700[k] = f_4 * hsh0_415[k]
                   - f_5 * hsh1_415[k]
                   + f_3 * pc_x[k] * hsi_548[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, pc_x, hsh0_416, hsh0_417, hsh0_418, hsh1_416, \
                         hsh1_417, hsh1_418, hsi_549, hsi_550, \
                         hsi_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_4 * hsh0_416[k]
                   - f_5 * hsh1_416[k]
                   + f_3 * pc_x[k] * hsi_549[k];

        t_702[k] = f_4 * hsh0_417[k]
                   - f_5 * hsh1_417[k]
                   + f_3 * pc_x[k] * hsi_550[k];

        t_703[k] = f_4 * hsh0_418[k]
                   - f_5 * hsh1_418[k]
                   + f_3 * pc_x[k] * hsi_551[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, t_709, pa_y, pc_x, pc_y, gsk0_524, \
                         gsk1_524, hsi_553, hsi_554, hsi_555, hsi_556, \
                         hsi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = pa_y[k] * gsk0_524[k]
                   - f_12 * pc_y[k] * gsk1_524[k];

        t_705[k] = f_3 * pc_x[k] * hsi_553[k];

        t_706[k] = f_3 * pc_x[k] * hsi_554[k];

        t_707[k] = f_3 * pc_x[k] * hsi_555[k];

        t_708[k] = f_3 * pc_x[k] * hsi_556[k];

        t_709[k] = f_3 * pc_x[k] * hsi_557[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pa_y, pc_x, pc_y, pc_z, gsk0_532, \
                         gsi_385, gsi_413, gsk1_532, hsi_553, hsi_558, \
                         hsi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_3 * pc_x[k] * hsi_558[k];

        t_711[k] = f_3 * pc_x[k] * hsi_559[k];

        t_712[k] = pa_y[k] * gsk0_532[k]
                   + f_19 * gsi_413[k]
                   - f_12 * pc_y[k] * gsk1_532[k];

        t_713[k] = f_16 * gsi_385[k]
                   + f_3 * pc_z[k] * hsi_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pa_y, pc_y, gsk0_534, gsk0_535, gsk0_536, \
                         gsi_415, gsi_416, gsi_417, gsk1_534, gsk1_535, \
                         gsk1_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pa_y[k] * gsk0_534[k]
                   + f_0 * gsi_415[k]
                   - f_12 * pc_y[k] * gsk1_534[k];

        t_715[k] = pa_y[k] * gsk0_535[k]
                   + f_16 * gsi_416[k]
                   - f_12 * pc_y[k] * gsk1_535[k];

        t_716[k] = pa_y[k] * gsk0_536[k]
                   + f_15 * gsi_417[k]
                   - f_12 * pc_y[k] * gsk1_536[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pa_y, pc_y, gsk0_537, gsk0_539, gsi_418, \
                         gsi_419, gsk1_537, gsk1_539, hsi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = pa_y[k] * gsk0_537[k]
                   + f_14 * gsi_418[k]
                   - f_12 * pc_y[k] * gsk1_537[k];

        t_718[k] = f_13 * gsi_419[k]
                   + f_3 * pc_y[k] * hsi_559[k];

        t_719[k] = pa_y[k] * gsk0_539[k]
                   - f_12 * pc_y[k] * gsk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, pc_x, pc_y, hsh0_420, hsh0_422, \
                         hsh0_423, hsh1_420, hsh1_422, hsh1_423, hsi_560, hsi_562, \
                         hsi_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_1 * hsh0_420[k]
                   - f_2 * hsh1_420[k]
                   + f_3 * pc_x[k] * hsi_560[k];

        t_721[k] = f_3 * pc_y[k] * hsi_560[k];

        t_722[k] = f_17 * hsh0_422[k]
                   - f_18 * hsh1_422[k]
                   + f_3 * pc_x[k] * hsi_562[k];

        t_723[k] = f_10 * hsh0_423[k]
                   - f_11 * hsh1_423[k]
                   + f_3 * pc_x[k] * hsi_563[k];

        t_724[k] = f_3 * pc_y[k] * hsi_562[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pc_x, pc_y, hsh0_425, hsh0_426, hsh0_427, \
                         hsh1_425, hsh1_426, hsh1_427, hsi_565, hsi_566, \
                         hsi_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_10 * hsh0_425[k]
                   - f_11 * hsh1_425[k]
                   + f_3 * pc_x[k] * hsi_565[k];

        t_726[k] = f_8 * hsh0_426[k]
                   - f_9 * hsh1_426[k]
                   + f_3 * pc_x[k] * hsi_566[k];

        t_727[k] = f_8 * hsh0_427[k]
                   - f_9 * hsh1_427[k]
                   + f_3 * pc_x[k] * hsi_567[k];

        t_728[k] = f_3 * pc_y[k] * hsi_565[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, hsh0_429, hsh0_430, hsh0_431, hsh1_429, \
                         hsh1_430, hsh1_431, hsi_569, hsi_570, \
                         hsi_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_8 * hsh0_429[k]
                   - f_9 * hsh1_429[k]
                   + f_3 * pc_x[k] * hsi_569[k];

        t_730[k] = f_6 * hsh0_430[k]
                   - f_7 * hsh1_430[k]
                   + f_3 * pc_x[k] * hsi_570[k];

        t_731[k] = f_6 * hsh0_431[k]
                   - f_7 * hsh1_431[k]
                   + f_3 * pc_x[k] * hsi_571[k];
    }
}

static auto
compute_prim_hsk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t gsi, const size_t hsh0,
                                                          const size_t hsh1, const size_t hsi,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_17 = 2.5 / gamma;
    const auto f_18 = 2.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsi_419 = buffer.data(gsi + 419);

    const auto *hsh0_432 = buffer.data(hsh0 + 432);
    const auto *hsh0_434 = buffer.data(hsh0 + 434);
    const auto *hsh0_435 = buffer.data(hsh0 + 435);
    const auto *hsh0_436 = buffer.data(hsh0 + 436);
    const auto *hsh0_437 = buffer.data(hsh0 + 437);
    const auto *hsh0_438 = buffer.data(hsh0 + 438);
    const auto *hsh0_439 = buffer.data(hsh0 + 439);
    const auto *hsh0_440 = buffer.data(hsh0 + 440);

    const auto *hsh1_432 = buffer.data(hsh1 + 432);
    const auto *hsh1_434 = buffer.data(hsh1 + 434);
    const auto *hsh1_435 = buffer.data(hsh1 + 435);
    const auto *hsh1_436 = buffer.data(hsh1 + 436);
    const auto *hsh1_437 = buffer.data(hsh1 + 437);
    const auto *hsh1_438 = buffer.data(hsh1 + 438);
    const auto *hsh1_439 = buffer.data(hsh1 + 439);
    const auto *hsh1_440 = buffer.data(hsh1 + 440);

    const auto *hsi_569 = buffer.data(hsi + 569);
    const auto *hsi_572 = buffer.data(hsi + 572);
    const auto *hsi_574 = buffer.data(hsi + 574);
    const auto *hsi_575 = buffer.data(hsi + 575);
    const auto *hsi_576 = buffer.data(hsi + 576);
    const auto *hsi_577 = buffer.data(hsi + 577);
    const auto *hsi_578 = buffer.data(hsi + 578);
    const auto *hsi_580 = buffer.data(hsi + 580);
    const auto *hsi_581 = buffer.data(hsi + 581);
    const auto *hsi_582 = buffer.data(hsi + 582);
    const auto *hsi_583 = buffer.data(hsi + 583);
    const auto *hsi_584 = buffer.data(hsi + 584);
    const auto *hsi_585 = buffer.data(hsi + 585);
    const auto *hsi_586 = buffer.data(hsi + 586);
    const auto *hsi_587 = buffer.data(hsi + 587);

#pragma omp simd aligned(t_732, t_733, t_734, t_735, pc_x, pc_y, hsh0_432, hsh0_434, hsh0_435, \
                         hsh1_432, hsh1_434, hsh1_435, hsi_569, hsi_572, hsi_574, \
                         hsi_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_6 * hsh0_432[k]
                   - f_7 * hsh1_432[k]
                   + f_3 * pc_x[k] * hsi_572[k];

        t_733[k] = f_3 * pc_y[k] * hsi_569[k];

        t_734[k] = f_6 * hsh0_434[k]
                   - f_7 * hsh1_434[k]
                   + f_3 * pc_x[k] * hsi_574[k];

        t_735[k] = f_4 * hsh0_435[k]
                   - f_5 * hsh1_435[k]
                   + f_3 * pc_x[k] * hsi_575[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, t_739, pc_x, pc_y, hsh0_436, hsh0_437, hsh0_438, \
                         hsh1_436, hsh1_437, hsh1_438, hsi_574, hsi_576, hsi_577, \
                         hsi_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_4 * hsh0_436[k]
                   - f_5 * hsh1_436[k]
                   + f_3 * pc_x[k] * hsi_576[k];

        t_737[k] = f_4 * hsh0_437[k]
                   - f_5 * hsh1_437[k]
                   + f_3 * pc_x[k] * hsi_577[k];

        t_738[k] = f_4 * hsh0_438[k]
                   - f_5 * hsh1_438[k]
                   + f_3 * pc_x[k] * hsi_578[k];

        t_739[k] = f_3 * pc_y[k] * hsi_574[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, t_745, pc_x, hsh0_440, hsh1_440, \
                         hsi_580, hsi_581, hsi_582, hsi_583, hsi_584, \
                         hsi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_4 * hsh0_440[k]
                   - f_5 * hsh1_440[k]
                   + f_3 * pc_x[k] * hsi_580[k];

        t_741[k] = f_3 * pc_x[k] * hsi_581[k];

        t_742[k] = f_3 * pc_x[k] * hsi_582[k];

        t_743[k] = f_3 * pc_x[k] * hsi_583[k];

        t_744[k] = f_3 * pc_x[k] * hsi_584[k];

        t_745[k] = f_3 * pc_x[k] * hsi_585[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, hsh0_435, hsh0_436, hsh1_435, \
                         hsh1_436, hsi_581, hsi_582, hsi_586, hsi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_3 * pc_x[k] * hsi_586[k];

        t_747[k] = f_3 * pc_x[k] * hsi_587[k];

        t_748[k] = f_1 * hsh0_435[k]
                   - f_2 * hsh1_435[k]
                   + f_3 * pc_y[k] * hsi_581[k];

        t_749[k] = f_17 * hsh0_436[k]
                   - f_18 * hsh1_436[k]
                   + f_3 * pc_y[k] * hsi_582[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_y, hsh0_437, hsh0_438, hsh0_439, hsh1_437, \
                         hsh1_438, hsh1_439, hsi_583, hsi_584, \
                         hsi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_10 * hsh0_437[k]
                   - f_11 * hsh1_437[k]
                   + f_3 * pc_y[k] * hsi_583[k];

        t_751[k] = f_8 * hsh0_438[k]
                   - f_9 * hsh1_438[k]
                   + f_3 * pc_y[k] * hsi_584[k];

        t_752[k] = f_6 * hsh0_439[k]
                   - f_7 * hsh1_439[k]
                   + f_3 * pc_y[k] * hsi_585[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pc_y, pc_z, gsi_419, hsh0_440, hsh1_440, \
                         hsi_586, hsi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_4 * hsh0_440[k]
                   - f_5 * hsh1_440[k]
                   + f_3 * pc_y[k] * hsi_586[k];

        t_754[k] = f_3 * pc_y[k] * hsi_587[k];

        t_755[k] = f_0 * gsi_419[k]
                   + f_1 * hsh0_440[k]
                   - f_2 * hsh1_440[k]
                   + f_3 * pc_z[k] * hsi_587[k];
    }
}

auto
compute_prim_hsk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t gsk0, const size_t gsi,
                                                   const size_t gsk1, const size_t hsh0,
                                                   const size_t hsh1, const size_t hsi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_hsk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, gsk0, gsi,
                                                              gsk1, hsh0, hsh1, hsi, ncols,
                                                              gamma, p, q);

    compute_prim_hsk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, gsk0, gsi,
                                                              gsk1, hsh0, hsh1, hsi, ncols,
                                                              gamma, p, q);

    compute_prim_hsk_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, gsk0, gsi,
                                                              gsk1, hsh0, hsh1, hsi, ncols,
                                                              gamma, p, q);

    compute_prim_hsk_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, gsk0, gsi,
                                                              gsk1, hsh0, hsh1, hsi, ncols,
                                                              gamma, p, q);

    compute_prim_hsk_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, gsk0, gsi,
                                                              gsk1, hsh0, hsh1, hsi, ncols,
                                                              gamma, p, q);

    compute_prim_hsk_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, gsk0, gsi,
                                                              gsk1, hsh0, hsh1, hsi, ncols,
                                                              gamma, p, q);

    compute_prim_hsk_three_center_electron_repulsion_0_piece6(buffer, target, pc, gsi, hsh0,
                                                              hsh1, hsi, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
