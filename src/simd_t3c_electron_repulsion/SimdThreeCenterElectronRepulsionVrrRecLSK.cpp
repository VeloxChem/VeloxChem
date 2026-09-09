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


#include "SimdThreeCenterElectronRepulsionVrrRecLSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsh0, const size_t lsh1,
                                                          const size_t lsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.0 / q;

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

    const auto *ksk0_0 = buffer.data(ksk0 + 0);
    const auto *ksk0_3 = buffer.data(ksk0 + 3);
    const auto *ksk0_5 = buffer.data(ksk0 + 5);
    const auto *ksk0_6 = buffer.data(ksk0 + 6);
    const auto *ksk0_9 = buffer.data(ksk0 + 9);
    const auto *ksk0_10 = buffer.data(ksk0 + 10);
    const auto *ksk0_14 = buffer.data(ksk0 + 14);
    const auto *ksk0_15 = buffer.data(ksk0 + 15);
    const auto *ksk0_20 = buffer.data(ksk0 + 20);
    const auto *ksk0_28 = buffer.data(ksk0 + 28);
    const auto *ksk0_35 = buffer.data(ksk0 + 35);

    const auto *ksi_0 = buffer.data(ksi + 0);
    const auto *ksi_1 = buffer.data(ksi + 1);
    const auto *ksi_2 = buffer.data(ksi + 2);
    const auto *ksi_3 = buffer.data(ksi + 3);
    const auto *ksi_5 = buffer.data(ksi + 5);
    const auto *ksi_6 = buffer.data(ksi + 6);
    const auto *ksi_9 = buffer.data(ksi + 9);
    const auto *ksi_10 = buffer.data(ksi + 10);
    const auto *ksi_14 = buffer.data(ksi + 14);
    const auto *ksi_21 = buffer.data(ksi + 21);
    const auto *ksi_23 = buffer.data(ksi + 23);
    const auto *ksi_24 = buffer.data(ksi + 24);
    const auto *ksi_25 = buffer.data(ksi + 25);
    const auto *ksi_27 = buffer.data(ksi + 27);
    const auto *ksi_28 = buffer.data(ksi + 28);
    const auto *ksi_33 = buffer.data(ksi + 33);
    const auto *ksi_37 = buffer.data(ksi + 37);
    const auto *ksi_42 = buffer.data(ksi + 42);
    const auto *ksi_49 = buffer.data(ksi + 49);
    const auto *ksi_51 = buffer.data(ksi + 51);
    const auto *ksi_52 = buffer.data(ksi + 52);
    const auto *ksi_53 = buffer.data(ksi + 53);
    const auto *ksi_54 = buffer.data(ksi + 54);
    const auto *ksi_55 = buffer.data(ksi + 55);
    const auto *ksi_77 = buffer.data(ksi + 77);
    const auto *ksi_78 = buffer.data(ksi + 78);
    const auto *ksi_79 = buffer.data(ksi + 79);
    const auto *ksi_80 = buffer.data(ksi + 80);
    const auto *ksi_81 = buffer.data(ksi + 81);
    const auto *ksi_83 = buffer.data(ksi + 83);
    const auto *ksi_84 = buffer.data(ksi + 84);
    const auto *ksi_87 = buffer.data(ksi + 87);
    const auto *ksi_90 = buffer.data(ksi + 90);
    const auto *ksi_94 = buffer.data(ksi + 94);
    const auto *ksi_99 = buffer.data(ksi + 99);

    const auto *ksk1_0 = buffer.data(ksk1 + 0);
    const auto *ksk1_3 = buffer.data(ksk1 + 3);
    const auto *ksk1_5 = buffer.data(ksk1 + 5);
    const auto *ksk1_6 = buffer.data(ksk1 + 6);
    const auto *ksk1_9 = buffer.data(ksk1 + 9);
    const auto *ksk1_10 = buffer.data(ksk1 + 10);
    const auto *ksk1_14 = buffer.data(ksk1 + 14);
    const auto *ksk1_15 = buffer.data(ksk1 + 15);
    const auto *ksk1_20 = buffer.data(ksk1 + 20);
    const auto *ksk1_28 = buffer.data(ksk1 + 28);
    const auto *ksk1_35 = buffer.data(ksk1 + 35);

    const auto *lsh0_0 = buffer.data(lsh0 + 0);
    const auto *lsh0_1 = buffer.data(lsh0 + 1);
    const auto *lsh0_2 = buffer.data(lsh0 + 2);
    const auto *lsh0_3 = buffer.data(lsh0 + 3);
    const auto *lsh0_5 = buffer.data(lsh0 + 5);
    const auto *lsh0_6 = buffer.data(lsh0 + 6);
    const auto *lsh0_8 = buffer.data(lsh0 + 8);
    const auto *lsh0_9 = buffer.data(lsh0 + 9);
    const auto *lsh0_15 = buffer.data(lsh0 + 15);
    const auto *lsh0_17 = buffer.data(lsh0 + 17);
    const auto *lsh0_18 = buffer.data(lsh0 + 18);
    const auto *lsh0_19 = buffer.data(lsh0 + 19);
    const auto *lsh0_20 = buffer.data(lsh0 + 20);
    const auto *lsh0_24 = buffer.data(lsh0 + 24);
    const auto *lsh0_27 = buffer.data(lsh0 + 27);
    const auto *lsh0_28 = buffer.data(lsh0 + 28);
    const auto *lsh0_36 = buffer.data(lsh0 + 36);
    const auto *lsh0_37 = buffer.data(lsh0 + 37);
    const auto *lsh0_38 = buffer.data(lsh0 + 38);
    const auto *lsh0_39 = buffer.data(lsh0 + 39);
    const auto *lsh0_44 = buffer.data(lsh0 + 44);
    const auto *lsh0_46 = buffer.data(lsh0 + 46);
    const auto *lsh0_47 = buffer.data(lsh0 + 47);
    const auto *lsh0_49 = buffer.data(lsh0 + 49);
    const auto *lsh0_50 = buffer.data(lsh0 + 50);
    const auto *lsh0_51 = buffer.data(lsh0 + 51);
    const auto *lsh0_58 = buffer.data(lsh0 + 58);
    const auto *lsh0_59 = buffer.data(lsh0 + 59);
    const auto *lsh0_60 = buffer.data(lsh0 + 60);
    const auto *lsh0_61 = buffer.data(lsh0 + 61);
    const auto *lsh0_62 = buffer.data(lsh0 + 62);
    const auto *lsh0_63 = buffer.data(lsh0 + 63);
    const auto *lsh0_65 = buffer.data(lsh0 + 65);
    const auto *lsh0_66 = buffer.data(lsh0 + 66);
    const auto *lsh0_68 = buffer.data(lsh0 + 68);
    const auto *lsh0_69 = buffer.data(lsh0 + 69);
    const auto *lsh0_70 = buffer.data(lsh0 + 70);
    const auto *lsh0_72 = buffer.data(lsh0 + 72);
    const auto *lsh0_73 = buffer.data(lsh0 + 73);
    const auto *lsh0_78 = buffer.data(lsh0 + 78);

    const auto *lsh1_0 = buffer.data(lsh1 + 0);
    const auto *lsh1_1 = buffer.data(lsh1 + 1);
    const auto *lsh1_2 = buffer.data(lsh1 + 2);
    const auto *lsh1_3 = buffer.data(lsh1 + 3);
    const auto *lsh1_5 = buffer.data(lsh1 + 5);
    const auto *lsh1_6 = buffer.data(lsh1 + 6);
    const auto *lsh1_8 = buffer.data(lsh1 + 8);
    const auto *lsh1_9 = buffer.data(lsh1 + 9);
    const auto *lsh1_15 = buffer.data(lsh1 + 15);
    const auto *lsh1_17 = buffer.data(lsh1 + 17);
    const auto *lsh1_18 = buffer.data(lsh1 + 18);
    const auto *lsh1_19 = buffer.data(lsh1 + 19);
    const auto *lsh1_20 = buffer.data(lsh1 + 20);
    const auto *lsh1_24 = buffer.data(lsh1 + 24);
    const auto *lsh1_27 = buffer.data(lsh1 + 27);
    const auto *lsh1_28 = buffer.data(lsh1 + 28);
    const auto *lsh1_36 = buffer.data(lsh1 + 36);
    const auto *lsh1_37 = buffer.data(lsh1 + 37);
    const auto *lsh1_38 = buffer.data(lsh1 + 38);
    const auto *lsh1_39 = buffer.data(lsh1 + 39);
    const auto *lsh1_44 = buffer.data(lsh1 + 44);
    const auto *lsh1_46 = buffer.data(lsh1 + 46);
    const auto *lsh1_47 = buffer.data(lsh1 + 47);
    const auto *lsh1_49 = buffer.data(lsh1 + 49);
    const auto *lsh1_50 = buffer.data(lsh1 + 50);
    const auto *lsh1_51 = buffer.data(lsh1 + 51);
    const auto *lsh1_58 = buffer.data(lsh1 + 58);
    const auto *lsh1_59 = buffer.data(lsh1 + 59);
    const auto *lsh1_60 = buffer.data(lsh1 + 60);
    const auto *lsh1_61 = buffer.data(lsh1 + 61);
    const auto *lsh1_62 = buffer.data(lsh1 + 62);
    const auto *lsh1_63 = buffer.data(lsh1 + 63);
    const auto *lsh1_65 = buffer.data(lsh1 + 65);
    const auto *lsh1_66 = buffer.data(lsh1 + 66);
    const auto *lsh1_68 = buffer.data(lsh1 + 68);
    const auto *lsh1_69 = buffer.data(lsh1 + 69);
    const auto *lsh1_70 = buffer.data(lsh1 + 70);
    const auto *lsh1_72 = buffer.data(lsh1 + 72);
    const auto *lsh1_73 = buffer.data(lsh1 + 73);
    const auto *lsh1_78 = buffer.data(lsh1 + 78);

    const auto *lsi_0 = buffer.data(lsi + 0);
    const auto *lsi_1 = buffer.data(lsi + 1);
    const auto *lsi_2 = buffer.data(lsi + 2);
    const auto *lsi_3 = buffer.data(lsi + 3);
    const auto *lsi_5 = buffer.data(lsi + 5);
    const auto *lsi_6 = buffer.data(lsi + 6);
    const auto *lsi_8 = buffer.data(lsi + 8);
    const auto *lsi_9 = buffer.data(lsi + 9);
    const auto *lsi_10 = buffer.data(lsi + 10);
    const auto *lsi_12 = buffer.data(lsi + 12);
    const auto *lsi_13 = buffer.data(lsi + 13);
    const auto *lsi_14 = buffer.data(lsi + 14);
    const auto *lsi_15 = buffer.data(lsi + 15);
    const auto *lsi_20 = buffer.data(lsi + 20);
    const auto *lsi_21 = buffer.data(lsi + 21);
    const auto *lsi_23 = buffer.data(lsi + 23);
    const auto *lsi_24 = buffer.data(lsi + 24);
    const auto *lsi_25 = buffer.data(lsi + 25);
    const auto *lsi_26 = buffer.data(lsi + 26);
    const auto *lsi_27 = buffer.data(lsi + 27);
    const auto *lsi_28 = buffer.data(lsi + 28);
    const auto *lsi_29 = buffer.data(lsi + 29);
    const auto *lsi_31 = buffer.data(lsi + 31);
    const auto *lsi_33 = buffer.data(lsi + 33);
    const auto *lsi_34 = buffer.data(lsi + 34);
    const auto *lsi_35 = buffer.data(lsi + 35);
    const auto *lsi_37 = buffer.data(lsi + 37);
    const auto *lsi_38 = buffer.data(lsi + 38);
    const auto *lsi_39 = buffer.data(lsi + 39);
    const auto *lsi_40 = buffer.data(lsi + 40);
    const auto *lsi_42 = buffer.data(lsi + 42);
    const auto *lsi_43 = buffer.data(lsi + 43);
    const auto *lsi_49 = buffer.data(lsi + 49);
    const auto *lsi_50 = buffer.data(lsi + 50);
    const auto *lsi_51 = buffer.data(lsi + 51);
    const auto *lsi_52 = buffer.data(lsi + 52);
    const auto *lsi_53 = buffer.data(lsi + 53);
    const auto *lsi_54 = buffer.data(lsi + 54);
    const auto *lsi_55 = buffer.data(lsi + 55);
    const auto *lsi_56 = buffer.data(lsi + 56);
    const auto *lsi_58 = buffer.data(lsi + 58);
    const auto *lsi_60 = buffer.data(lsi + 60);
    const auto *lsi_61 = buffer.data(lsi + 61);
    const auto *lsi_63 = buffer.data(lsi + 63);
    const auto *lsi_64 = buffer.data(lsi + 64);
    const auto *lsi_65 = buffer.data(lsi + 65);
    const auto *lsi_67 = buffer.data(lsi + 67);
    const auto *lsi_68 = buffer.data(lsi + 68);
    const auto *lsi_69 = buffer.data(lsi + 69);
    const auto *lsi_70 = buffer.data(lsi + 70);
    const auto *lsi_76 = buffer.data(lsi + 76);
    const auto *lsi_77 = buffer.data(lsi + 77);
    const auto *lsi_78 = buffer.data(lsi + 78);
    const auto *lsi_79 = buffer.data(lsi + 79);
    const auto *lsi_80 = buffer.data(lsi + 80);
    const auto *lsi_81 = buffer.data(lsi + 81);
    const auto *lsi_82 = buffer.data(lsi + 82);
    const auto *lsi_83 = buffer.data(lsi + 83);
    const auto *lsi_84 = buffer.data(lsi + 84);
    const auto *lsi_85 = buffer.data(lsi + 85);
    const auto *lsi_86 = buffer.data(lsi + 86);
    const auto *lsi_87 = buffer.data(lsi + 87);
    const auto *lsi_89 = buffer.data(lsi + 89);
    const auto *lsi_90 = buffer.data(lsi + 90);
    const auto *lsi_91 = buffer.data(lsi + 91);
    const auto *lsi_93 = buffer.data(lsi + 93);
    const auto *lsi_94 = buffer.data(lsi + 94);
    const auto *lsi_95 = buffer.data(lsi + 95);
    const auto *lsi_96 = buffer.data(lsi + 96);
    const auto *lsi_98 = buffer.data(lsi + 98);
    const auto *lsi_99 = buffer.data(lsi + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ksi_0, lsh0_0, \
                         lsh1_0, lsi_0, lsi_1, lsi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ksi_0[k]
                 + f_1 * lsh0_0[k]
                 - f_2 * lsh1_0[k]
                 + f_3 * pc_x[k] * lsi_0[k];

        t_1[k] = f_3 * pc_y[k] * lsi_0[k];

        t_2[k] = f_3 * pc_z[k] * lsi_0[k];

        t_3[k] = f_4 * lsh0_0[k]
                 - f_5 * lsh1_0[k]
                 + f_3 * pc_y[k] * lsi_1[k];

        t_4[k] = f_3 * pc_y[k] * lsi_2[k];

        t_5[k] = f_4 * lsh0_0[k]
                 - f_5 * lsh1_0[k]
                 + f_3 * pc_z[k] * lsi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, lsh0_1, lsh0_2, lsh0_3, lsh1_1, \
                         lsh1_2, lsh1_3, lsi_3, lsi_5, lsi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * lsh0_1[k]
                 - f_7 * lsh1_1[k]
                 + f_3 * pc_y[k] * lsi_3[k];

        t_7[k] = f_3 * pc_z[k] * lsi_3[k];

        t_8[k] = f_3 * pc_y[k] * lsi_5[k];

        t_9[k] = f_6 * lsh0_2[k]
                 - f_7 * lsh1_2[k]
                 + f_3 * pc_z[k] * lsi_5[k];

        t_10[k] = f_8 * lsh0_3[k]
                  - f_9 * lsh1_3[k]
                  + f_3 * pc_y[k] * lsi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, lsh0_5, lsh0_6, \
                         lsh1_5, lsh1_6, lsi_6, lsi_8, lsi_9, lsi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * lsi_6[k];

        t_12[k] = f_4 * lsh0_5[k]
                  - f_5 * lsh1_5[k]
                  + f_3 * pc_y[k] * lsi_8[k];

        t_13[k] = f_3 * pc_y[k] * lsi_9[k];

        t_14[k] = f_8 * lsh0_5[k]
                  - f_9 * lsh1_5[k]
                  + f_3 * pc_z[k] * lsi_9[k];

        t_15[k] = f_10 * lsh0_6[k]
                  - f_11 * lsh1_6[k]
                  + f_3 * pc_y[k] * lsi_10[k];

        t_16[k] = f_3 * pc_z[k] * lsi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, lsh0_8, lsh0_9, lsh1_8, lsh1_9, \
                         lsi_12, lsi_13, lsi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * lsh0_8[k]
                  - f_7 * lsh1_8[k]
                  + f_3 * pc_y[k] * lsi_12[k];

        t_18[k] = f_4 * lsh0_9[k]
                  - f_5 * lsh1_9[k]
                  + f_3 * pc_y[k] * lsi_13[k];

        t_19[k] = f_3 * pc_y[k] * lsi_14[k];

        t_20[k] = f_10 * lsh0_9[k]
                  - f_11 * lsh1_9[k]
                  + f_3 * pc_z[k] * lsi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pc_x, pc_z, ksi_21, ksi_23, ksi_24, \
                         ksi_25, lsi_15, lsi_21, lsi_23, lsi_24, \
                         lsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * ksi_21[k]
                  + f_3 * pc_x[k] * lsi_21[k];

        t_22[k] = f_3 * pc_z[k] * lsi_15[k];

        t_23[k] = f_0 * ksi_23[k]
                  + f_3 * pc_x[k] * lsi_23[k];

        t_24[k] = f_0 * ksi_24[k]
                  + f_3 * pc_x[k] * lsi_24[k];

        t_25[k] = f_0 * ksi_25[k]
                  + f_3 * pc_x[k] * lsi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, ksi_27, lsh0_15, lsh1_15, \
                         lsi_20, lsi_21, lsi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * pc_y[k] * lsi_20[k];

        t_27[k] = f_0 * ksi_27[k]
                  + f_3 * pc_x[k] * lsi_27[k];

        t_28[k] = f_1 * lsh0_15[k]
                  - f_2 * lsh1_15[k]
                  + f_3 * pc_y[k] * lsi_21[k];

        t_29[k] = f_3 * pc_z[k] * lsi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pc_y, lsh0_17, lsh0_18, lsh0_19, lsh1_17, lsh1_18, \
                         lsh1_19, lsi_23, lsi_24, lsi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * lsh0_17[k]
                  - f_11 * lsh1_17[k]
                  + f_3 * pc_y[k] * lsi_23[k];

        t_31[k] = f_8 * lsh0_18[k]
                  - f_9 * lsh1_18[k]
                  + f_3 * pc_y[k] * lsi_24[k];

        t_32[k] = f_6 * lsh0_19[k]
                  - f_7 * lsh1_19[k]
                  + f_3 * pc_y[k] * lsi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, ksk0_0, ksi_0, \
                         ksk1_0, lsh0_20, lsh1_20, lsi_26, lsi_27, \
                         lsi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * lsh0_20[k]
                  - f_5 * lsh1_20[k]
                  + f_3 * pc_y[k] * lsi_26[k];

        t_34[k] = f_3 * pc_y[k] * lsi_27[k];

        t_35[k] = f_1 * lsh0_20[k]
                  - f_2 * lsh1_20[k]
                  + f_3 * pc_z[k] * lsi_27[k];

        t_36[k] = pa_y[k] * ksk0_0[k]
                  - f_12 * pc_y[k] * ksk1_0[k];

        t_37[k] = f_13 * ksi_0[k]
                  + f_3 * pc_y[k] * lsi_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, ksk0_3, ksk0_5, ksi_1, \
                         ksk1_3, ksk1_5, lsi_28, lsi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * pc_z[k] * lsi_28[k];

        t_39[k] = pa_y[k] * ksk0_3[k]
                  + f_14 * ksi_1[k]
                  - f_12 * pc_y[k] * ksk1_3[k];

        t_40[k] = f_3 * pc_z[k] * lsi_29[k];

        t_41[k] = pa_y[k] * ksk0_5[k]
                  - f_12 * pc_y[k] * ksk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_y, pc_z, ksk0_6, ksk0_9, ksi_3, \
                         ksi_5, ksk1_6, ksk1_9, lsi_31, lsi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * ksk0_6[k]
                  + f_15 * ksi_3[k]
                  - f_12 * pc_y[k] * ksk1_6[k];

        t_43[k] = f_3 * pc_z[k] * lsi_31[k];

        t_44[k] = f_13 * ksi_5[k]
                  + f_3 * pc_y[k] * lsi_33[k];

        t_45[k] = pa_y[k] * ksk0_9[k]
                  - f_12 * pc_y[k] * ksk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, ksk0_10, ksi_6, ksi_9, \
                         ksk1_10, lsh0_24, lsh1_24, lsi_34, lsi_35, \
                         lsi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * ksk0_10[k]
                  + f_16 * ksi_6[k]
                  - f_12 * pc_y[k] * ksk1_10[k];

        t_47[k] = f_3 * pc_z[k] * lsi_34[k];

        t_48[k] = f_4 * lsh0_24[k]
                  - f_5 * lsh1_24[k]
                  + f_3 * pc_z[k] * lsi_35[k];

        t_49[k] = f_13 * ksi_9[k]
                  + f_3 * pc_y[k] * lsi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_y, pc_z, ksk0_14, ksk0_15, ksi_10, \
                         ksk1_14, ksk1_15, lsh0_27, lsh1_27, lsi_38, \
                         lsi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * ksk0_14[k]
                  - f_12 * pc_y[k] * ksk1_14[k];

        t_51[k] = pa_y[k] * ksk0_15[k]
                  + f_17 * ksi_10[k]
                  - f_12 * pc_y[k] * ksk1_15[k];

        t_52[k] = f_3 * pc_z[k] * lsi_38[k];

        t_53[k] = f_4 * lsh0_27[k]
                  - f_5 * lsh1_27[k]
                  + f_3 * pc_z[k] * lsi_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_y, pc_z, ksk0_20, ksi_14, ksk1_20, \
                         lsh0_28, lsh1_28, lsi_40, lsi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * lsh0_28[k]
                  - f_7 * lsh1_28[k]
                  + f_3 * pc_z[k] * lsi_40[k];

        t_55[k] = f_13 * ksi_14[k]
                  + f_3 * pc_y[k] * lsi_42[k];

        t_56[k] = pa_y[k] * ksk0_20[k]
                  - f_12 * pc_y[k] * ksk1_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pc_x, pc_z, ksi_49, ksi_51, ksi_52, \
                         ksi_53, lsi_43, lsi_49, lsi_51, lsi_52, \
                         lsi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_18 * ksi_49[k]
                  + f_3 * pc_x[k] * lsi_49[k];

        t_58[k] = f_3 * pc_z[k] * lsi_43[k];

        t_59[k] = f_18 * ksi_51[k]
                  + f_3 * pc_x[k] * lsi_51[k];

        t_60[k] = f_18 * ksi_52[k]
                  + f_3 * pc_x[k] * lsi_52[k];

        t_61[k] = f_18 * ksi_53[k]
                  + f_3 * pc_x[k] * lsi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, ksi_21, ksi_54, ksi_55, \
                         lsh0_36, lsh1_36, lsi_49, lsi_54, lsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_18 * ksi_54[k]
                  + f_3 * pc_x[k] * lsi_54[k];

        t_63[k] = f_18 * ksi_55[k]
                  + f_3 * pc_x[k] * lsi_55[k];

        t_64[k] = f_13 * ksi_21[k]
                  + f_1 * lsh0_36[k]
                  - f_2 * lsh1_36[k]
                  + f_3 * pc_y[k] * lsi_49[k];

        t_65[k] = f_3 * pc_z[k] * lsi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_z, lsh0_36, lsh0_37, lsh0_38, lsh1_36, lsh1_37, \
                         lsh1_38, lsi_50, lsi_51, lsi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * lsh0_36[k]
                  - f_5 * lsh1_36[k]
                  + f_3 * pc_z[k] * lsi_50[k];

        t_67[k] = f_6 * lsh0_37[k]
                  - f_7 * lsh1_37[k]
                  + f_3 * pc_z[k] * lsi_51[k];

        t_68[k] = f_8 * lsh0_38[k]
                  - f_9 * lsh1_38[k]
                  + f_3 * pc_z[k] * lsi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, ksk0_35, ksi_27, ksk1_35, \
                         lsh0_39, lsh1_39, lsi_53, lsi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * lsh0_39[k]
                  - f_11 * lsh1_39[k]
                  + f_3 * pc_z[k] * lsi_53[k];

        t_70[k] = f_13 * ksi_27[k]
                  + f_3 * pc_y[k] * lsi_55[k];

        t_71[k] = pa_y[k] * ksk0_35[k]
                  - f_12 * pc_y[k] * ksk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pc_y, pc_z, ksk0_0, ksk0_3, \
                         ksi_0, ksk1_0, ksk1_3, lsi_56, lsi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * ksk0_0[k]
                  - f_12 * pc_z[k] * ksk1_0[k];

        t_73[k] = f_3 * pc_y[k] * lsi_56[k];

        t_74[k] = f_13 * ksi_0[k]
                  + f_3 * pc_z[k] * lsi_56[k];

        t_75[k] = pa_z[k] * ksk0_3[k]
                  - f_12 * pc_z[k] * ksk1_3[k];

        t_76[k] = f_3 * pc_y[k] * lsi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pc_y, pc_z, ksk0_5, ksk0_6, ksi_2, \
                         ksk1_5, ksk1_6, lsh0_44, lsh1_44, lsi_60, \
                         lsi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * ksk0_5[k]
                  + f_14 * ksi_2[k]
                  - f_12 * pc_z[k] * ksk1_5[k];

        t_78[k] = pa_z[k] * ksk0_6[k]
                  - f_12 * pc_z[k] * ksk1_6[k];

        t_79[k] = f_4 * lsh0_44[k]
                  - f_5 * lsh1_44[k]
                  + f_3 * pc_y[k] * lsi_60[k];

        t_80[k] = f_3 * pc_y[k] * lsi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, ksk0_9, ksk0_10, ksi_5, ksk1_9, \
                         ksk1_10, lsh0_46, lsh1_46, lsi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * ksk0_9[k]
                  + f_15 * ksi_5[k]
                  - f_12 * pc_z[k] * ksk1_9[k];

        t_82[k] = pa_z[k] * ksk0_10[k]
                  - f_12 * pc_z[k] * ksk1_10[k];

        t_83[k] = f_6 * lsh0_46[k]
                  - f_7 * lsh1_46[k]
                  + f_3 * pc_y[k] * lsi_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pc_y, pc_z, ksk0_14, ksk0_15, ksi_9, \
                         ksk1_14, ksk1_15, lsh0_47, lsh1_47, lsi_64, \
                         lsi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * lsh0_47[k]
                  - f_5 * lsh1_47[k]
                  + f_3 * pc_y[k] * lsi_64[k];

        t_85[k] = f_3 * pc_y[k] * lsi_65[k];

        t_86[k] = pa_z[k] * ksk0_14[k]
                  + f_16 * ksi_9[k]
                  - f_12 * pc_z[k] * ksk1_14[k];

        t_87[k] = pa_z[k] * ksk0_15[k]
                  - f_12 * pc_z[k] * ksk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, lsh0_49, lsh0_50, lsh0_51, lsh1_49, \
                         lsh1_50, lsh1_51, lsi_67, lsi_68, lsi_69, \
                         lsi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_8 * lsh0_49[k]
                  - f_9 * lsh1_49[k]
                  + f_3 * pc_y[k] * lsi_67[k];

        t_89[k] = f_6 * lsh0_50[k]
                  - f_7 * lsh1_50[k]
                  + f_3 * pc_y[k] * lsi_68[k];

        t_90[k] = f_4 * lsh0_51[k]
                  - f_5 * lsh1_51[k]
                  + f_3 * pc_y[k] * lsi_69[k];

        t_91[k] = f_3 * pc_y[k] * lsi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_z, ksk0_20, ksi_14, ksi_77, \
                         ksi_78, ksi_79, ksk1_20, lsi_77, lsi_78, \
                         lsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * ksk0_20[k]
                  + f_17 * ksi_14[k]
                  - f_12 * pc_z[k] * ksk1_20[k];

        t_93[k] = f_18 * ksi_77[k]
                  + f_3 * pc_x[k] * lsi_77[k];

        t_94[k] = f_18 * ksi_78[k]
                  + f_3 * pc_x[k] * lsi_78[k];

        t_95[k] = f_18 * ksi_79[k]
                  + f_3 * pc_x[k] * lsi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, ksi_80, ksi_81, ksi_83, lsi_76, \
                         lsi_80, lsi_81, lsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_18 * ksi_80[k]
                  + f_3 * pc_x[k] * lsi_80[k];

        t_97[k] = f_18 * ksi_81[k]
                  + f_3 * pc_x[k] * lsi_81[k];

        t_98[k] = f_3 * pc_y[k] * lsi_76[k];

        t_99[k] = f_18 * ksi_83[k]
                  + f_3 * pc_x[k] * lsi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_y, pc_z, ksk0_28, ksk1_28, lsh0_58, \
                         lsh0_59, lsh1_58, lsh1_59, lsi_78, lsi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * ksk0_28[k]
                   - f_12 * pc_z[k] * ksk1_28[k];

        t_101[k] = f_19 * lsh0_58[k]
                   - f_20 * lsh1_58[k]
                   + f_3 * pc_y[k] * lsi_78[k];

        t_102[k] = f_10 * lsh0_59[k]
                   - f_11 * lsh1_59[k]
                   + f_3 * pc_y[k] * lsi_79[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_y, lsh0_60, lsh0_61, lsh0_62, lsh1_60, \
                         lsh1_61, lsh1_62, lsi_80, lsi_81, lsi_82, \
                         lsi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * lsh0_60[k]
                   - f_9 * lsh1_60[k]
                   + f_3 * pc_y[k] * lsi_80[k];

        t_104[k] = f_6 * lsh0_61[k]
                   - f_7 * lsh1_61[k]
                   + f_3 * pc_y[k] * lsi_81[k];

        t_105[k] = f_4 * lsh0_62[k]
                   - f_5 * lsh1_62[k]
                   + f_3 * pc_y[k] * lsi_82[k];

        t_106[k] = f_3 * pc_y[k] * lsi_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, pc_z, ksi_27, ksi_28, ksi_84, \
                         lsh0_62, lsh0_63, lsh1_62, lsh1_63, lsi_83, \
                         lsi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_13 * ksi_27[k]
                   + f_1 * lsh0_62[k]
                   - f_2 * lsh1_62[k]
                   + f_3 * pc_z[k] * lsi_83[k];

        t_108[k] = f_21 * ksi_84[k]
                   + f_1 * lsh0_63[k]
                   - f_2 * lsh1_63[k]
                   + f_3 * pc_x[k] * lsi_84[k];

        t_109[k] = f_14 * ksi_28[k]
                   + f_3 * pc_y[k] * lsi_84[k];

        t_110[k] = f_3 * pc_z[k] * lsi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_z, ksi_87, lsh0_63, lsh0_66, lsh1_63, \
                         lsh1_66, lsi_85, lsi_86, lsi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_21 * ksi_87[k]
                   + f_10 * lsh0_66[k]
                   - f_11 * lsh1_66[k]
                   + f_3 * pc_x[k] * lsi_87[k];

        t_112[k] = f_3 * pc_z[k] * lsi_85[k];

        t_113[k] = f_4 * lsh0_63[k]
                   - f_5 * lsh1_63[k]
                   + f_3 * pc_z[k] * lsi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, pc_y, pc_z, ksi_33, ksi_90, \
                         lsh0_65, lsh0_69, lsh1_65, lsh1_69, lsi_87, lsi_89, \
                         lsi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_21 * ksi_90[k]
                   + f_8 * lsh0_69[k]
                   - f_9 * lsh1_69[k]
                   + f_3 * pc_x[k] * lsi_90[k];

        t_115[k] = f_3 * pc_z[k] * lsi_87[k];

        t_116[k] = f_14 * ksi_33[k]
                   + f_3 * pc_y[k] * lsi_89[k];

        t_117[k] = f_6 * lsh0_65[k]
                   - f_7 * lsh1_65[k]
                   + f_3 * pc_z[k] * lsi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_z, ksi_94, lsh0_66, lsh0_73, lsh1_66, \
                         lsh1_73, lsi_90, lsi_91, lsi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_21 * ksi_94[k]
                   + f_6 * lsh0_73[k]
                   - f_7 * lsh1_73[k]
                   + f_3 * pc_x[k] * lsi_94[k];

        t_119[k] = f_3 * pc_z[k] * lsi_90[k];

        t_120[k] = f_4 * lsh0_66[k]
                   - f_5 * lsh1_66[k]
                   + f_3 * pc_z[k] * lsi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, ksi_37, ksi_99, \
                         lsh0_68, lsh0_78, lsh1_68, lsh1_78, lsi_93, lsi_94, \
                         lsi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_14 * ksi_37[k]
                   + f_3 * pc_y[k] * lsi_93[k];

        t_122[k] = f_8 * lsh0_68[k]
                   - f_9 * lsh1_68[k]
                   + f_3 * pc_z[k] * lsi_93[k];

        t_123[k] = f_21 * ksi_99[k]
                   + f_4 * lsh0_78[k]
                   - f_5 * lsh1_78[k]
                   + f_3 * pc_x[k] * lsi_99[k];

        t_124[k] = f_3 * pc_z[k] * lsi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_y, pc_z, ksi_42, lsh0_69, lsh0_70, \
                         lsh0_72, lsh1_69, lsh1_70, lsh1_72, lsi_95, lsi_96, \
                         lsi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * lsh0_69[k]
                   - f_5 * lsh1_69[k]
                   + f_3 * pc_z[k] * lsi_95[k];

        t_126[k] = f_6 * lsh0_70[k]
                   - f_7 * lsh1_70[k]
                   + f_3 * pc_z[k] * lsi_96[k];

        t_127[k] = f_14 * ksi_42[k]
                   + f_3 * pc_y[k] * lsi_98[k];

        t_128[k] = f_10 * lsh0_72[k]
                   - f_11 * lsh1_72[k]
                   + f_3 * pc_z[k] * lsi_98[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsh0, const size_t lsh1,
                                                          const size_t lsi, const size_t ncols,
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
    const auto f_17 = 2.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.0 / q;

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

    const auto *ksk0_39 = buffer.data(ksk0 + 39);
    const auto *ksk0_42 = buffer.data(ksk0 + 42);
    const auto *ksk0_46 = buffer.data(ksk0 + 46);
    const auto *ksk0_51 = buffer.data(ksk0 + 51);
    const auto *ksk0_64 = buffer.data(ksk0 + 64);
    const auto *ksk0_72 = buffer.data(ksk0 + 72);
    const auto *ksk0_77 = buffer.data(ksk0 + 77);
    const auto *ksk0_81 = buffer.data(ksk0 + 81);
    const auto *ksk0_84 = buffer.data(ksk0 + 84);
    const auto *ksk0_86 = buffer.data(ksk0 + 86);
    const auto *ksk0_89 = buffer.data(ksk0 + 89);
    const auto *ksk0_90 = buffer.data(ksk0 + 90);
    const auto *ksk0_92 = buffer.data(ksk0 + 92);
    const auto *ksk0_107 = buffer.data(ksk0 + 107);

    const auto *ksi_28 = buffer.data(ksi + 28);
    const auto *ksi_31 = buffer.data(ksi + 31);
    const auto *ksi_34 = buffer.data(ksi + 34);
    const auto *ksi_38 = buffer.data(ksi + 38);
    const auto *ksi_49 = buffer.data(ksi + 49);
    const auto *ksi_55 = buffer.data(ksi + 55);
    const auto *ksi_56 = buffer.data(ksi + 56);
    const auto *ksi_58 = buffer.data(ksi + 58);
    const auto *ksi_61 = buffer.data(ksi + 61);
    const auto *ksi_64 = buffer.data(ksi + 64);
    const auto *ksi_65 = buffer.data(ksi + 65);
    const auto *ksi_68 = buffer.data(ksi + 68);
    const auto *ksi_69 = buffer.data(ksi + 69);
    const auto *ksi_70 = buffer.data(ksi + 70);
    const auto *ksi_79 = buffer.data(ksi + 79);
    const auto *ksi_80 = buffer.data(ksi + 80);
    const auto *ksi_81 = buffer.data(ksi + 81);
    const auto *ksi_82 = buffer.data(ksi + 82);
    const auto *ksi_83 = buffer.data(ksi + 83);
    const auto *ksi_84 = buffer.data(ksi + 84);
    const auto *ksi_89 = buffer.data(ksi + 89);
    const auto *ksi_93 = buffer.data(ksi + 93);
    const auto *ksi_98 = buffer.data(ksi + 98);
    const auto *ksi_105 = buffer.data(ksi + 105);
    const auto *ksi_107 = buffer.data(ksi + 107);
    const auto *ksi_108 = buffer.data(ksi + 108);
    const auto *ksi_109 = buffer.data(ksi + 109);
    const auto *ksi_110 = buffer.data(ksi + 110);
    const auto *ksi_111 = buffer.data(ksi + 111);
    const auto *ksi_133 = buffer.data(ksi + 133);
    const auto *ksi_134 = buffer.data(ksi + 134);
    const auto *ksi_135 = buffer.data(ksi + 135);
    const auto *ksi_136 = buffer.data(ksi + 136);
    const auto *ksi_137 = buffer.data(ksi + 137);
    const auto *ksi_138 = buffer.data(ksi + 138);
    const auto *ksi_139 = buffer.data(ksi + 139);
    const auto *ksi_140 = buffer.data(ksi + 140);
    const auto *ksi_145 = buffer.data(ksi + 145);
    const auto *ksi_149 = buffer.data(ksi + 149);
    const auto *ksi_154 = buffer.data(ksi + 154);
    const auto *ksi_160 = buffer.data(ksi + 160);
    const auto *ksi_161 = buffer.data(ksi + 161);
    const auto *ksi_162 = buffer.data(ksi + 162);
    const auto *ksi_163 = buffer.data(ksi + 163);
    const auto *ksi_164 = buffer.data(ksi + 164);
    const auto *ksi_165 = buffer.data(ksi + 165);
    const auto *ksi_167 = buffer.data(ksi + 167);
    const auto *ksi_168 = buffer.data(ksi + 168);
    const auto *ksi_171 = buffer.data(ksi + 171);
    const auto *ksi_174 = buffer.data(ksi + 174);
    const auto *ksi_178 = buffer.data(ksi + 178);
    const auto *ksi_183 = buffer.data(ksi + 183);
    const auto *ksi_189 = buffer.data(ksi + 189);
    const auto *ksi_191 = buffer.data(ksi + 191);
    const auto *ksi_192 = buffer.data(ksi + 192);
    const auto *ksi_193 = buffer.data(ksi + 193);
    const auto *ksi_194 = buffer.data(ksi + 194);
    const auto *ksi_195 = buffer.data(ksi + 195);

    const auto *ksk1_39 = buffer.data(ksk1 + 39);
    const auto *ksk1_42 = buffer.data(ksk1 + 42);
    const auto *ksk1_46 = buffer.data(ksk1 + 46);
    const auto *ksk1_51 = buffer.data(ksk1 + 51);
    const auto *ksk1_64 = buffer.data(ksk1 + 64);
    const auto *ksk1_72 = buffer.data(ksk1 + 72);
    const auto *ksk1_77 = buffer.data(ksk1 + 77);
    const auto *ksk1_81 = buffer.data(ksk1 + 81);
    const auto *ksk1_84 = buffer.data(ksk1 + 84);
    const auto *ksk1_86 = buffer.data(ksk1 + 86);
    const auto *ksk1_89 = buffer.data(ksk1 + 89);
    const auto *ksk1_90 = buffer.data(ksk1 + 90);
    const auto *ksk1_92 = buffer.data(ksk1 + 92);
    const auto *ksk1_107 = buffer.data(ksk1 + 107);

    const auto *lsh0_78 = buffer.data(lsh0 + 78);
    const auto *lsh0_79 = buffer.data(lsh0 + 79);
    const auto *lsh0_80 = buffer.data(lsh0 + 80);
    const auto *lsh0_81 = buffer.data(lsh0 + 81);
    const auto *lsh0_83 = buffer.data(lsh0 + 83);
    const auto *lsh0_101 = buffer.data(lsh0 + 101);
    const auto *lsh0_102 = buffer.data(lsh0 + 102);
    const auto *lsh0_103 = buffer.data(lsh0 + 103);
    const auto *lsh0_104 = buffer.data(lsh0 + 104);
    const auto *lsh0_105 = buffer.data(lsh0 + 105);
    const auto *lsh0_106 = buffer.data(lsh0 + 106);
    const auto *lsh0_107 = buffer.data(lsh0 + 107);
    const auto *lsh0_108 = buffer.data(lsh0 + 108);
    const auto *lsh0_109 = buffer.data(lsh0 + 109);
    const auto *lsh0_110 = buffer.data(lsh0 + 110);
    const auto *lsh0_111 = buffer.data(lsh0 + 111);
    const auto *lsh0_112 = buffer.data(lsh0 + 112);
    const auto *lsh0_113 = buffer.data(lsh0 + 113);
    const auto *lsh0_114 = buffer.data(lsh0 + 114);
    const auto *lsh0_119 = buffer.data(lsh0 + 119);
    const auto *lsh0_120 = buffer.data(lsh0 + 120);
    const auto *lsh0_121 = buffer.data(lsh0 + 121);
    const auto *lsh0_122 = buffer.data(lsh0 + 122);
    const auto *lsh0_123 = buffer.data(lsh0 + 123);
    const auto *lsh0_124 = buffer.data(lsh0 + 124);
    const auto *lsh0_125 = buffer.data(lsh0 + 125);
    const auto *lsh0_126 = buffer.data(lsh0 + 126);
    const auto *lsh0_128 = buffer.data(lsh0 + 128);
    const auto *lsh0_129 = buffer.data(lsh0 + 129);
    const auto *lsh0_131 = buffer.data(lsh0 + 131);
    const auto *lsh0_132 = buffer.data(lsh0 + 132);
    const auto *lsh0_133 = buffer.data(lsh0 + 133);
    const auto *lsh0_135 = buffer.data(lsh0 + 135);
    const auto *lsh0_136 = buffer.data(lsh0 + 136);
    const auto *lsh0_141 = buffer.data(lsh0 + 141);
    const auto *lsh0_142 = buffer.data(lsh0 + 142);
    const auto *lsh0_143 = buffer.data(lsh0 + 143);
    const auto *lsh0_144 = buffer.data(lsh0 + 144);

    const auto *lsh1_78 = buffer.data(lsh1 + 78);
    const auto *lsh1_79 = buffer.data(lsh1 + 79);
    const auto *lsh1_80 = buffer.data(lsh1 + 80);
    const auto *lsh1_81 = buffer.data(lsh1 + 81);
    const auto *lsh1_83 = buffer.data(lsh1 + 83);
    const auto *lsh1_101 = buffer.data(lsh1 + 101);
    const auto *lsh1_102 = buffer.data(lsh1 + 102);
    const auto *lsh1_103 = buffer.data(lsh1 + 103);
    const auto *lsh1_104 = buffer.data(lsh1 + 104);
    const auto *lsh1_105 = buffer.data(lsh1 + 105);
    const auto *lsh1_106 = buffer.data(lsh1 + 106);
    const auto *lsh1_107 = buffer.data(lsh1 + 107);
    const auto *lsh1_108 = buffer.data(lsh1 + 108);
    const auto *lsh1_109 = buffer.data(lsh1 + 109);
    const auto *lsh1_110 = buffer.data(lsh1 + 110);
    const auto *lsh1_111 = buffer.data(lsh1 + 111);
    const auto *lsh1_112 = buffer.data(lsh1 + 112);
    const auto *lsh1_113 = buffer.data(lsh1 + 113);
    const auto *lsh1_114 = buffer.data(lsh1 + 114);
    const auto *lsh1_119 = buffer.data(lsh1 + 119);
    const auto *lsh1_120 = buffer.data(lsh1 + 120);
    const auto *lsh1_121 = buffer.data(lsh1 + 121);
    const auto *lsh1_122 = buffer.data(lsh1 + 122);
    const auto *lsh1_123 = buffer.data(lsh1 + 123);
    const auto *lsh1_124 = buffer.data(lsh1 + 124);
    const auto *lsh1_125 = buffer.data(lsh1 + 125);
    const auto *lsh1_126 = buffer.data(lsh1 + 126);
    const auto *lsh1_128 = buffer.data(lsh1 + 128);
    const auto *lsh1_129 = buffer.data(lsh1 + 129);
    const auto *lsh1_131 = buffer.data(lsh1 + 131);
    const auto *lsh1_132 = buffer.data(lsh1 + 132);
    const auto *lsh1_133 = buffer.data(lsh1 + 133);
    const auto *lsh1_135 = buffer.data(lsh1 + 135);
    const auto *lsh1_136 = buffer.data(lsh1 + 136);
    const auto *lsh1_141 = buffer.data(lsh1 + 141);
    const auto *lsh1_142 = buffer.data(lsh1 + 142);
    const auto *lsh1_143 = buffer.data(lsh1 + 143);
    const auto *lsh1_144 = buffer.data(lsh1 + 144);

    const auto *lsi_99 = buffer.data(lsi + 99);
    const auto *lsi_105 = buffer.data(lsi + 105);
    const auto *lsi_106 = buffer.data(lsi + 106);
    const auto *lsi_107 = buffer.data(lsi + 107);
    const auto *lsi_108 = buffer.data(lsi + 108);
    const auto *lsi_109 = buffer.data(lsi + 109);
    const auto *lsi_110 = buffer.data(lsi + 110);
    const auto *lsi_111 = buffer.data(lsi + 111);
    const auto *lsi_112 = buffer.data(lsi + 112);
    const auto *lsi_114 = buffer.data(lsi + 114);
    const auto *lsi_115 = buffer.data(lsi + 115);
    const auto *lsi_117 = buffer.data(lsi + 117);
    const auto *lsi_118 = buffer.data(lsi + 118);
    const auto *lsi_121 = buffer.data(lsi + 121);
    const auto *lsi_122 = buffer.data(lsi + 122);
    const auto *lsi_126 = buffer.data(lsi + 126);
    const auto *lsi_133 = buffer.data(lsi + 133);
    const auto *lsi_134 = buffer.data(lsi + 134);
    const auto *lsi_135 = buffer.data(lsi + 135);
    const auto *lsi_136 = buffer.data(lsi + 136);
    const auto *lsi_137 = buffer.data(lsi + 137);
    const auto *lsi_138 = buffer.data(lsi + 138);
    const auto *lsi_139 = buffer.data(lsi + 139);
    const auto *lsi_140 = buffer.data(lsi + 140);
    const auto *lsi_141 = buffer.data(lsi + 141);
    const auto *lsi_142 = buffer.data(lsi + 142);
    const auto *lsi_143 = buffer.data(lsi + 143);
    const auto *lsi_144 = buffer.data(lsi + 144);
    const auto *lsi_145 = buffer.data(lsi + 145);
    const auto *lsi_146 = buffer.data(lsi + 146);
    const auto *lsi_147 = buffer.data(lsi + 147);
    const auto *lsi_148 = buffer.data(lsi + 148);
    const auto *lsi_149 = buffer.data(lsi + 149);
    const auto *lsi_150 = buffer.data(lsi + 150);
    const auto *lsi_151 = buffer.data(lsi + 151);
    const auto *lsi_152 = buffer.data(lsi + 152);
    const auto *lsi_153 = buffer.data(lsi + 153);
    const auto *lsi_154 = buffer.data(lsi + 154);
    const auto *lsi_160 = buffer.data(lsi + 160);
    const auto *lsi_161 = buffer.data(lsi + 161);
    const auto *lsi_162 = buffer.data(lsi + 162);
    const auto *lsi_163 = buffer.data(lsi + 163);
    const auto *lsi_164 = buffer.data(lsi + 164);
    const auto *lsi_165 = buffer.data(lsi + 165);
    const auto *lsi_166 = buffer.data(lsi + 166);
    const auto *lsi_167 = buffer.data(lsi + 167);
    const auto *lsi_168 = buffer.data(lsi + 168);
    const auto *lsi_169 = buffer.data(lsi + 169);
    const auto *lsi_170 = buffer.data(lsi + 170);
    const auto *lsi_171 = buffer.data(lsi + 171);
    const auto *lsi_173 = buffer.data(lsi + 173);
    const auto *lsi_174 = buffer.data(lsi + 174);
    const auto *lsi_175 = buffer.data(lsi + 175);
    const auto *lsi_177 = buffer.data(lsi + 177);
    const auto *lsi_178 = buffer.data(lsi + 178);
    const auto *lsi_179 = buffer.data(lsi + 179);
    const auto *lsi_180 = buffer.data(lsi + 180);
    const auto *lsi_182 = buffer.data(lsi + 182);
    const auto *lsi_183 = buffer.data(lsi + 183);
    const auto *lsi_189 = buffer.data(lsi + 189);
    const auto *lsi_190 = buffer.data(lsi + 190);
    const auto *lsi_191 = buffer.data(lsi + 191);
    const auto *lsi_192 = buffer.data(lsi + 192);
    const auto *lsi_193 = buffer.data(lsi + 193);
    const auto *lsi_194 = buffer.data(lsi + 194);
    const auto *lsi_195 = buffer.data(lsi + 195);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pc_x, pc_z, ksi_105, ksi_107, \
                         ksi_108, ksi_109, lsi_99, lsi_105, lsi_107, lsi_108, \
                         lsi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_21 * ksi_105[k]
                   + f_3 * pc_x[k] * lsi_105[k];

        t_130[k] = f_3 * pc_z[k] * lsi_99[k];

        t_131[k] = f_21 * ksi_107[k]
                   + f_3 * pc_x[k] * lsi_107[k];

        t_132[k] = f_21 * ksi_108[k]
                   + f_3 * pc_x[k] * lsi_108[k];

        t_133[k] = f_21 * ksi_109[k]
                   + f_3 * pc_x[k] * lsi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, ksi_49, ksi_110, \
                         ksi_111, lsh0_78, lsh1_78, lsi_105, lsi_110, \
                         lsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_21 * ksi_110[k]
                   + f_3 * pc_x[k] * lsi_110[k];

        t_135[k] = f_21 * ksi_111[k]
                   + f_3 * pc_x[k] * lsi_111[k];

        t_136[k] = f_14 * ksi_49[k]
                   + f_1 * lsh0_78[k]
                   - f_2 * lsh1_78[k]
                   + f_3 * pc_y[k] * lsi_105[k];

        t_137[k] = f_3 * pc_z[k] * lsi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_z, lsh0_78, lsh0_79, lsh0_80, lsh1_78, \
                         lsh1_79, lsh1_80, lsi_106, lsi_107, lsi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_4 * lsh0_78[k]
                   - f_5 * lsh1_78[k]
                   + f_3 * pc_z[k] * lsi_106[k];

        t_139[k] = f_6 * lsh0_79[k]
                   - f_7 * lsh1_79[k]
                   + f_3 * pc_z[k] * lsi_107[k];

        t_140[k] = f_8 * lsh0_80[k]
                   - f_9 * lsh1_80[k]
                   + f_3 * pc_z[k] * lsi_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pc_y, pc_z, ksk0_72, ksi_55, \
                         ksk1_72, lsh0_81, lsh0_83, lsh1_81, lsh1_83, lsi_109, \
                         lsi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_10 * lsh0_81[k]
                   - f_11 * lsh1_81[k]
                   + f_3 * pc_z[k] * lsi_109[k];

        t_142[k] = f_14 * ksi_55[k]
                   + f_3 * pc_y[k] * lsi_111[k];

        t_143[k] = f_1 * lsh0_83[k]
                   - f_2 * lsh1_83[k]
                   + f_3 * pc_z[k] * lsi_111[k];

        t_144[k] = pa_y[k] * ksk0_72[k]
                   - f_12 * pc_y[k] * ksk1_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_z, pc_y, pc_z, ksk0_39, ksi_28, \
                         ksi_56, ksi_58, ksk1_39, lsi_112, lsi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * ksi_56[k]
                   + f_3 * pc_y[k] * lsi_112[k];

        t_146[k] = f_13 * ksi_28[k]
                   + f_3 * pc_z[k] * lsi_112[k];

        t_147[k] = pa_z[k] * ksk0_39[k]
                   - f_12 * pc_z[k] * ksk1_39[k];

        t_148[k] = f_13 * ksi_58[k]
                   + f_3 * pc_y[k] * lsi_114[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pa_z, pc_y, pc_z, ksk0_42, ksk0_77, \
                         ksi_31, ksi_61, ksk1_42, ksk1_77, lsi_115, \
                         lsi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * ksk0_77[k]
                   - f_12 * pc_y[k] * ksk1_77[k];

        t_150[k] = pa_z[k] * ksk0_42[k]
                   - f_12 * pc_z[k] * ksk1_42[k];

        t_151[k] = f_13 * ksi_31[k]
                   + f_3 * pc_z[k] * lsi_115[k];

        t_152[k] = f_13 * ksi_61[k]
                   + f_3 * pc_y[k] * lsi_117[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_y, pa_z, pc_y, pc_z, ksk0_46, ksk0_81, \
                         ksi_34, ksk1_46, ksk1_81, lsi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_y[k] * ksk0_81[k]
                   - f_12 * pc_y[k] * ksk1_81[k];

        t_154[k] = pa_z[k] * ksk0_46[k]
                   - f_12 * pc_z[k] * ksk1_46[k];

        t_155[k] = f_13 * ksi_34[k]
                   + f_3 * pc_z[k] * lsi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_y, pc_y, ksk0_84, ksk0_86, ksi_64, ksi_65, \
                         ksk1_84, ksk1_86, lsi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_y[k] * ksk0_84[k]
                   + f_14 * ksi_64[k]
                   - f_12 * pc_y[k] * ksk1_84[k];

        t_157[k] = f_13 * ksi_65[k]
                   + f_3 * pc_y[k] * lsi_121[k];

        t_158[k] = pa_y[k] * ksk0_86[k]
                   - f_12 * pc_y[k] * ksk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pa_z, pc_y, pc_z, ksk0_51, ksk0_89, \
                         ksi_38, ksi_68, ksk1_51, ksk1_89, lsi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * ksk0_51[k]
                   - f_12 * pc_z[k] * ksk1_51[k];

        t_160[k] = f_13 * ksi_38[k]
                   + f_3 * pc_z[k] * lsi_122[k];

        t_161[k] = pa_y[k] * ksk0_89[k]
                   + f_15 * ksi_68[k]
                   - f_12 * pc_y[k] * ksk1_89[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_y, pc_x, pc_y, ksk0_90, ksk0_92, \
                         ksi_69, ksi_70, ksi_133, ksk1_90, ksk1_92, lsi_126, \
                         lsi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_y[k] * ksk0_90[k]
                   + f_14 * ksi_69[k]
                   - f_12 * pc_y[k] * ksk1_90[k];

        t_163[k] = f_13 * ksi_70[k]
                   + f_3 * pc_y[k] * lsi_126[k];

        t_164[k] = pa_y[k] * ksk0_92[k]
                   - f_12 * pc_y[k] * ksk1_92[k];

        t_165[k] = f_21 * ksi_133[k]
                   + f_3 * pc_x[k] * lsi_133[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, ksi_134, ksi_135, ksi_136, \
                         ksi_137, ksi_138, lsi_134, lsi_135, lsi_136, lsi_137, \
                         lsi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_21 * ksi_134[k]
                   + f_3 * pc_x[k] * lsi_134[k];

        t_167[k] = f_21 * ksi_135[k]
                   + f_3 * pc_x[k] * lsi_135[k];

        t_168[k] = f_21 * ksi_136[k]
                   + f_3 * pc_x[k] * lsi_136[k];

        t_169[k] = f_21 * ksi_137[k]
                   + f_3 * pc_x[k] * lsi_137[k];

        t_170[k] = f_21 * ksi_138[k]
                   + f_3 * pc_x[k] * lsi_138[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_z, pc_x, pc_z, ksk0_64, ksi_49, ksi_139, \
                         ksk1_64, lsi_133, lsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_21 * ksi_139[k]
                   + f_3 * pc_x[k] * lsi_139[k];

        t_172[k] = pa_z[k] * ksk0_64[k]
                   - f_12 * pc_z[k] * ksk1_64[k];

        t_173[k] = f_13 * ksi_49[k]
                   + f_3 * pc_z[k] * lsi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, ksi_79, ksi_80, ksi_81, lsh0_101, \
                         lsh0_102, lsh0_103, lsh1_101, lsh1_102, lsh1_103, lsi_135, lsi_136, \
                         lsi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * ksi_79[k]
                   + f_10 * lsh0_101[k]
                   - f_11 * lsh1_101[k]
                   + f_3 * pc_y[k] * lsi_135[k];

        t_175[k] = f_13 * ksi_80[k]
                   + f_8 * lsh0_102[k]
                   - f_9 * lsh1_102[k]
                   + f_3 * pc_y[k] * lsi_136[k];

        t_176[k] = f_13 * ksi_81[k]
                   + f_6 * lsh0_103[k]
                   - f_7 * lsh1_103[k]
                   + f_3 * pc_y[k] * lsi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, ksk0_107, ksi_82, ksi_83, ksk1_107, \
                         lsh0_104, lsh1_104, lsi_138, lsi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * ksi_82[k]
                   + f_4 * lsh0_104[k]
                   - f_5 * lsh1_104[k]
                   + f_3 * pc_y[k] * lsi_138[k];

        t_178[k] = f_13 * ksi_83[k]
                   + f_3 * pc_y[k] * lsi_139[k];

        t_179[k] = pa_y[k] * ksk0_107[k]
                   - f_12 * pc_y[k] * ksk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, ksi_56, ksi_140, \
                         lsh0_105, lsh1_105, lsi_140, lsi_141, \
                         lsi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_21 * ksi_140[k]
                   + f_1 * lsh0_105[k]
                   - f_2 * lsh1_105[k]
                   + f_3 * pc_x[k] * lsi_140[k];

        t_181[k] = f_3 * pc_y[k] * lsi_140[k];

        t_182[k] = f_14 * ksi_56[k]
                   + f_3 * pc_z[k] * lsi_140[k];

        t_183[k] = f_4 * lsh0_105[k]
                   - f_5 * lsh1_105[k]
                   + f_3 * pc_y[k] * lsi_141[k];

        t_184[k] = f_3 * pc_y[k] * lsi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, ksi_145, lsh0_106, lsh0_107, \
                         lsh0_110, lsh1_106, lsh1_107, lsh1_110, lsi_143, lsi_144, \
                         lsi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_21 * ksi_145[k]
                   + f_10 * lsh0_110[k]
                   - f_11 * lsh1_110[k]
                   + f_3 * pc_x[k] * lsi_145[k];

        t_186[k] = f_6 * lsh0_106[k]
                   - f_7 * lsh1_106[k]
                   + f_3 * pc_y[k] * lsi_143[k];

        t_187[k] = f_4 * lsh0_107[k]
                   - f_5 * lsh1_107[k]
                   + f_3 * pc_y[k] * lsi_144[k];

        t_188[k] = f_3 * pc_y[k] * lsi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, ksi_149, lsh0_108, lsh0_109, \
                         lsh0_114, lsh1_108, lsh1_109, lsh1_114, lsi_146, lsi_147, \
                         lsi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_21 * ksi_149[k]
                   + f_8 * lsh0_114[k]
                   - f_9 * lsh1_114[k]
                   + f_3 * pc_x[k] * lsi_149[k];

        t_190[k] = f_8 * lsh0_108[k]
                   - f_9 * lsh1_108[k]
                   + f_3 * pc_y[k] * lsi_146[k];

        t_191[k] = f_6 * lsh0_109[k]
                   - f_7 * lsh1_109[k]
                   + f_3 * pc_y[k] * lsi_147[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pc_x, pc_y, ksi_154, lsh0_110, lsh0_119, \
                         lsh1_110, lsh1_119, lsi_148, lsi_149, \
                         lsi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_4 * lsh0_110[k]
                   - f_5 * lsh1_110[k]
                   + f_3 * pc_y[k] * lsi_148[k];

        t_193[k] = f_3 * pc_y[k] * lsi_149[k];

        t_194[k] = f_21 * ksi_154[k]
                   + f_6 * lsh0_119[k]
                   - f_7 * lsh1_119[k]
                   + f_3 * pc_x[k] * lsi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_y, lsh0_111, lsh0_112, lsh0_113, lsh1_111, \
                         lsh1_112, lsh1_113, lsi_150, lsi_151, \
                         lsi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_10 * lsh0_111[k]
                   - f_11 * lsh1_111[k]
                   + f_3 * pc_y[k] * lsi_150[k];

        t_196[k] = f_8 * lsh0_112[k]
                   - f_9 * lsh1_112[k]
                   + f_3 * pc_y[k] * lsi_151[k];

        t_197[k] = f_6 * lsh0_113[k]
                   - f_7 * lsh1_113[k]
                   + f_3 * pc_y[k] * lsi_152[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, ksi_160, ksi_161, lsh0_114, \
                         lsh0_125, lsh1_114, lsh1_125, lsi_153, lsi_154, lsi_160, \
                         lsi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * lsh0_114[k]
                   - f_5 * lsh1_114[k]
                   + f_3 * pc_y[k] * lsi_153[k];

        t_199[k] = f_3 * pc_y[k] * lsi_154[k];

        t_200[k] = f_21 * ksi_160[k]
                   + f_4 * lsh0_125[k]
                   - f_5 * lsh1_125[k]
                   + f_3 * pc_x[k] * lsi_160[k];

        t_201[k] = f_21 * ksi_161[k]
                   + f_3 * pc_x[k] * lsi_161[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pc_x, pc_y, ksi_162, ksi_163, \
                         ksi_164, ksi_165, lsi_160, lsi_162, lsi_163, lsi_164, \
                         lsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_21 * ksi_162[k]
                   + f_3 * pc_x[k] * lsi_162[k];

        t_203[k] = f_21 * ksi_163[k]
                   + f_3 * pc_x[k] * lsi_163[k];

        t_204[k] = f_21 * ksi_164[k]
                   + f_3 * pc_x[k] * lsi_164[k];

        t_205[k] = f_21 * ksi_165[k]
                   + f_3 * pc_x[k] * lsi_165[k];

        t_206[k] = f_3 * pc_y[k] * lsi_160[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pc_x, pc_y, ksi_167, lsh0_120, lsh0_121, \
                         lsh1_120, lsh1_121, lsi_161, lsi_162, \
                         lsi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_21 * ksi_167[k]
                   + f_3 * pc_x[k] * lsi_167[k];

        t_208[k] = f_1 * lsh0_120[k]
                   - f_2 * lsh1_120[k]
                   + f_3 * pc_y[k] * lsi_161[k];

        t_209[k] = f_19 * lsh0_121[k]
                   - f_20 * lsh1_121[k]
                   + f_3 * pc_y[k] * lsi_162[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pc_y, lsh0_122, lsh0_123, lsh0_124, lsh1_122, \
                         lsh1_123, lsh1_124, lsi_163, lsi_164, \
                         lsi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * lsh0_122[k]
                   - f_11 * lsh1_122[k]
                   + f_3 * pc_y[k] * lsi_163[k];

        t_211[k] = f_8 * lsh0_123[k]
                   - f_9 * lsh1_123[k]
                   + f_3 * pc_y[k] * lsi_164[k];

        t_212[k] = f_6 * lsh0_124[k]
                   - f_7 * lsh1_124[k]
                   + f_3 * pc_y[k] * lsi_165[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, pc_z, ksi_83, ksi_168, \
                         lsh0_125, lsh0_126, lsh1_125, lsh1_126, lsi_166, lsi_167, \
                         lsi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * lsh0_125[k]
                   - f_5 * lsh1_125[k]
                   + f_3 * pc_y[k] * lsi_166[k];

        t_214[k] = f_3 * pc_y[k] * lsi_167[k];

        t_215[k] = f_14 * ksi_83[k]
                   + f_1 * lsh0_125[k]
                   - f_2 * lsh1_125[k]
                   + f_3 * pc_z[k] * lsi_167[k];

        t_216[k] = f_17 * ksi_168[k]
                   + f_1 * lsh0_126[k]
                   - f_2 * lsh1_126[k]
                   + f_3 * pc_x[k] * lsi_168[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, ksi_84, ksi_171, \
                         lsh0_129, lsh1_129, lsi_168, lsi_169, \
                         lsi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * ksi_84[k]
                   + f_3 * pc_y[k] * lsi_168[k];

        t_218[k] = f_3 * pc_z[k] * lsi_168[k];

        t_219[k] = f_17 * ksi_171[k]
                   + f_10 * lsh0_129[k]
                   - f_11 * lsh1_129[k]
                   + f_3 * pc_x[k] * lsi_171[k];

        t_220[k] = f_3 * pc_z[k] * lsi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pc_x, pc_z, ksi_174, lsh0_126, lsh0_132, \
                         lsh1_126, lsh1_132, lsi_170, lsi_171, \
                         lsi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_4 * lsh0_126[k]
                   - f_5 * lsh1_126[k]
                   + f_3 * pc_z[k] * lsi_170[k];

        t_222[k] = f_17 * ksi_174[k]
                   + f_8 * lsh0_132[k]
                   - f_9 * lsh1_132[k]
                   + f_3 * pc_x[k] * lsi_174[k];

        t_223[k] = f_3 * pc_z[k] * lsi_171[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pc_x, pc_y, pc_z, ksi_89, ksi_178, \
                         lsh0_128, lsh0_136, lsh1_128, lsh1_136, lsi_173, lsi_174, \
                         lsi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_15 * ksi_89[k]
                   + f_3 * pc_y[k] * lsi_173[k];

        t_225[k] = f_6 * lsh0_128[k]
                   - f_7 * lsh1_128[k]
                   + f_3 * pc_z[k] * lsi_173[k];

        t_226[k] = f_17 * ksi_178[k]
                   + f_6 * lsh0_136[k]
                   - f_7 * lsh1_136[k]
                   + f_3 * pc_x[k] * lsi_178[k];

        t_227[k] = f_3 * pc_z[k] * lsi_174[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_y, pc_z, ksi_93, lsh0_129, lsh0_131, \
                         lsh1_129, lsh1_131, lsi_175, lsi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * lsh0_129[k]
                   - f_5 * lsh1_129[k]
                   + f_3 * pc_z[k] * lsi_175[k];

        t_229[k] = f_15 * ksi_93[k]
                   + f_3 * pc_y[k] * lsi_177[k];

        t_230[k] = f_8 * lsh0_131[k]
                   - f_9 * lsh1_131[k]
                   + f_3 * pc_z[k] * lsi_177[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pc_x, pc_z, ksi_183, lsh0_132, lsh0_141, \
                         lsh1_132, lsh1_141, lsi_178, lsi_179, \
                         lsi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_17 * ksi_183[k]
                   + f_4 * lsh0_141[k]
                   - f_5 * lsh1_141[k]
                   + f_3 * pc_x[k] * lsi_183[k];

        t_232[k] = f_3 * pc_z[k] * lsi_178[k];

        t_233[k] = f_4 * lsh0_132[k]
                   - f_5 * lsh1_132[k]
                   + f_3 * pc_z[k] * lsi_179[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, ksi_98, ksi_189, \
                         lsh0_133, lsh0_135, lsh1_133, lsh1_135, lsi_180, lsi_182, \
                         lsi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_6 * lsh0_133[k]
                   - f_7 * lsh1_133[k]
                   + f_3 * pc_z[k] * lsi_180[k];

        t_235[k] = f_15 * ksi_98[k]
                   + f_3 * pc_y[k] * lsi_182[k];

        t_236[k] = f_10 * lsh0_135[k]
                   - f_11 * lsh1_135[k]
                   + f_3 * pc_z[k] * lsi_182[k];

        t_237[k] = f_17 * ksi_189[k]
                   + f_3 * pc_x[k] * lsi_189[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pc_x, pc_z, ksi_191, ksi_192, \
                         ksi_193, ksi_194, lsi_183, lsi_191, lsi_192, lsi_193, \
                         lsi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_3 * pc_z[k] * lsi_183[k];

        t_239[k] = f_17 * ksi_191[k]
                   + f_3 * pc_x[k] * lsi_191[k];

        t_240[k] = f_17 * ksi_192[k]
                   + f_3 * pc_x[k] * lsi_192[k];

        t_241[k] = f_17 * ksi_193[k]
                   + f_3 * pc_x[k] * lsi_193[k];

        t_242[k] = f_17 * ksi_194[k]
                   + f_3 * pc_x[k] * lsi_194[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pc_x, pc_y, pc_z, ksi_105, ksi_195, \
                         lsh0_141, lsh1_141, lsi_189, lsi_190, \
                         lsi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_17 * ksi_195[k]
                   + f_3 * pc_x[k] * lsi_195[k];

        t_244[k] = f_15 * ksi_105[k]
                   + f_1 * lsh0_141[k]
                   - f_2 * lsh1_141[k]
                   + f_3 * pc_y[k] * lsi_189[k];

        t_245[k] = f_3 * pc_z[k] * lsi_189[k];

        t_246[k] = f_4 * lsh0_141[k]
                   - f_5 * lsh1_141[k]
                   + f_3 * pc_z[k] * lsi_190[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_z, lsh0_142, lsh0_143, lsh0_144, lsh1_142, \
                         lsh1_143, lsh1_144, lsi_191, lsi_192, \
                         lsi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * lsh0_142[k]
                   - f_7 * lsh1_142[k]
                   + f_3 * pc_z[k] * lsi_191[k];

        t_248[k] = f_8 * lsh0_143[k]
                   - f_9 * lsh1_143[k]
                   + f_3 * pc_z[k] * lsi_192[k];

        t_249[k] = f_10 * lsh0_144[k]
                   - f_11 * lsh1_144[k]
                   + f_3 * pc_z[k] * lsi_193[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsh0, const size_t lsh1,
                                                          const size_t lsi, const size_t ncols,
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
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_108 = buffer.data(ksk0 + 108);
    const auto *ksk0_111 = buffer.data(ksk0 + 111);
    const auto *ksk0_114 = buffer.data(ksk0 + 114);
    const auto *ksk0_118 = buffer.data(ksk0 + 118);
    const auto *ksk0_120 = buffer.data(ksk0 + 120);
    const auto *ksk0_123 = buffer.data(ksk0 + 123);
    const auto *ksk0_125 = buffer.data(ksk0 + 125);
    const auto *ksk0_126 = buffer.data(ksk0 + 126);
    const auto *ksk0_136 = buffer.data(ksk0 + 136);
    const auto *ksk0_180 = buffer.data(ksk0 + 180);
    const auto *ksk0_183 = buffer.data(ksk0 + 183);
    const auto *ksk0_185 = buffer.data(ksk0 + 185);
    const auto *ksk0_186 = buffer.data(ksk0 + 186);
    const auto *ksk0_189 = buffer.data(ksk0 + 189);
    const auto *ksk0_190 = buffer.data(ksk0 + 190);
    const auto *ksk0_192 = buffer.data(ksk0 + 192);
    const auto *ksk0_194 = buffer.data(ksk0 + 194);
    const auto *ksk0_195 = buffer.data(ksk0 + 195);
    const auto *ksk0_197 = buffer.data(ksk0 + 197);
    const auto *ksk0_198 = buffer.data(ksk0 + 198);
    const auto *ksk0_200 = buffer.data(ksk0 + 200);
    const auto *ksk0_215 = buffer.data(ksk0 + 215);

    const auto *ksi_84 = buffer.data(ksi + 84);
    const auto *ksi_87 = buffer.data(ksi + 87);
    const auto *ksi_90 = buffer.data(ksi + 90);
    const auto *ksi_91 = buffer.data(ksi + 91);
    const auto *ksi_94 = buffer.data(ksi + 94);
    const auto *ksi_95 = buffer.data(ksi + 95);
    const auto *ksi_96 = buffer.data(ksi + 96);
    const auto *ksi_105 = buffer.data(ksi + 105);
    const auto *ksi_111 = buffer.data(ksi + 111);
    const auto *ksi_112 = buffer.data(ksi + 112);
    const auto *ksi_114 = buffer.data(ksi + 114);
    const auto *ksi_115 = buffer.data(ksi + 115);
    const auto *ksi_117 = buffer.data(ksi + 117);
    const auto *ksi_118 = buffer.data(ksi + 118);
    const auto *ksi_121 = buffer.data(ksi + 121);
    const auto *ksi_122 = buffer.data(ksi + 122);
    const auto *ksi_126 = buffer.data(ksi + 126);
    const auto *ksi_133 = buffer.data(ksi + 133);
    const auto *ksi_135 = buffer.data(ksi + 135);
    const auto *ksi_136 = buffer.data(ksi + 136);
    const auto *ksi_137 = buffer.data(ksi + 137);
    const auto *ksi_138 = buffer.data(ksi + 138);
    const auto *ksi_139 = buffer.data(ksi + 139);
    const auto *ksi_140 = buffer.data(ksi + 140);
    const auto *ksi_141 = buffer.data(ksi + 141);
    const auto *ksi_142 = buffer.data(ksi + 142);
    const auto *ksi_143 = buffer.data(ksi + 143);
    const auto *ksi_145 = buffer.data(ksi + 145);
    const auto *ksi_146 = buffer.data(ksi + 146);
    const auto *ksi_148 = buffer.data(ksi + 148);
    const auto *ksi_149 = buffer.data(ksi + 149);
    const auto *ksi_150 = buffer.data(ksi + 150);
    const auto *ksi_152 = buffer.data(ksi + 152);
    const auto *ksi_153 = buffer.data(ksi + 153);
    const auto *ksi_154 = buffer.data(ksi + 154);
    const auto *ksi_161 = buffer.data(ksi + 161);
    const auto *ksi_163 = buffer.data(ksi + 163);
    const auto *ksi_164 = buffer.data(ksi + 164);
    const auto *ksi_165 = buffer.data(ksi + 165);
    const auto *ksi_166 = buffer.data(ksi + 166);
    const auto *ksi_167 = buffer.data(ksi + 167);
    const auto *ksi_168 = buffer.data(ksi + 168);
    const auto *ksi_201 = buffer.data(ksi + 201);
    const auto *ksi_205 = buffer.data(ksi + 205);
    const auto *ksi_210 = buffer.data(ksi + 210);
    const auto *ksi_216 = buffer.data(ksi + 216);
    const auto *ksi_217 = buffer.data(ksi + 217);
    const auto *ksi_218 = buffer.data(ksi + 218);
    const auto *ksi_219 = buffer.data(ksi + 219);
    const auto *ksi_220 = buffer.data(ksi + 220);
    const auto *ksi_221 = buffer.data(ksi + 221);
    const auto *ksi_222 = buffer.data(ksi + 222);
    const auto *ksi_223 = buffer.data(ksi + 223);
    const auto *ksi_245 = buffer.data(ksi + 245);
    const auto *ksi_246 = buffer.data(ksi + 246);
    const auto *ksi_247 = buffer.data(ksi + 247);
    const auto *ksi_248 = buffer.data(ksi + 248);
    const auto *ksi_249 = buffer.data(ksi + 249);
    const auto *ksi_250 = buffer.data(ksi + 250);
    const auto *ksi_251 = buffer.data(ksi + 251);
    const auto *ksi_252 = buffer.data(ksi + 252);
    const auto *ksi_257 = buffer.data(ksi + 257);
    const auto *ksi_261 = buffer.data(ksi + 261);
    const auto *ksi_266 = buffer.data(ksi + 266);
    const auto *ksi_272 = buffer.data(ksi + 272);
    const auto *ksi_273 = buffer.data(ksi + 273);
    const auto *ksi_274 = buffer.data(ksi + 274);
    const auto *ksi_275 = buffer.data(ksi + 275);
    const auto *ksi_276 = buffer.data(ksi + 276);
    const auto *ksi_277 = buffer.data(ksi + 277);
    const auto *ksi_279 = buffer.data(ksi + 279);
    const auto *ksi_280 = buffer.data(ksi + 280);
    const auto *ksi_283 = buffer.data(ksi + 283);
    const auto *ksi_286 = buffer.data(ksi + 286);

    const auto *ksk1_108 = buffer.data(ksk1 + 108);
    const auto *ksk1_111 = buffer.data(ksk1 + 111);
    const auto *ksk1_114 = buffer.data(ksk1 + 114);
    const auto *ksk1_118 = buffer.data(ksk1 + 118);
    const auto *ksk1_120 = buffer.data(ksk1 + 120);
    const auto *ksk1_123 = buffer.data(ksk1 + 123);
    const auto *ksk1_125 = buffer.data(ksk1 + 125);
    const auto *ksk1_126 = buffer.data(ksk1 + 126);
    const auto *ksk1_136 = buffer.data(ksk1 + 136);
    const auto *ksk1_180 = buffer.data(ksk1 + 180);
    const auto *ksk1_183 = buffer.data(ksk1 + 183);
    const auto *ksk1_185 = buffer.data(ksk1 + 185);
    const auto *ksk1_186 = buffer.data(ksk1 + 186);
    const auto *ksk1_189 = buffer.data(ksk1 + 189);
    const auto *ksk1_190 = buffer.data(ksk1 + 190);
    const auto *ksk1_192 = buffer.data(ksk1 + 192);
    const auto *ksk1_194 = buffer.data(ksk1 + 194);
    const auto *ksk1_195 = buffer.data(ksk1 + 195);
    const auto *ksk1_197 = buffer.data(ksk1 + 197);
    const auto *ksk1_198 = buffer.data(ksk1 + 198);
    const auto *ksk1_200 = buffer.data(ksk1 + 200);
    const auto *ksk1_215 = buffer.data(ksk1 + 215);

    const auto *lsh0_146 = buffer.data(lsh0 + 146);
    const auto *lsh0_152 = buffer.data(lsh0 + 152);
    const auto *lsh0_156 = buffer.data(lsh0 + 156);
    const auto *lsh0_161 = buffer.data(lsh0 + 161);
    const auto *lsh0_164 = buffer.data(lsh0 + 164);
    const auto *lsh0_165 = buffer.data(lsh0 + 165);
    const auto *lsh0_166 = buffer.data(lsh0 + 166);
    const auto *lsh0_167 = buffer.data(lsh0 + 167);
    const auto *lsh0_183 = buffer.data(lsh0 + 183);
    const auto *lsh0_185 = buffer.data(lsh0 + 185);
    const auto *lsh0_186 = buffer.data(lsh0 + 186);
    const auto *lsh0_187 = buffer.data(lsh0 + 187);
    const auto *lsh0_188 = buffer.data(lsh0 + 188);
    const auto *lsh0_189 = buffer.data(lsh0 + 189);
    const auto *lsh0_190 = buffer.data(lsh0 + 190);
    const auto *lsh0_191 = buffer.data(lsh0 + 191);
    const auto *lsh0_192 = buffer.data(lsh0 + 192);
    const auto *lsh0_193 = buffer.data(lsh0 + 193);
    const auto *lsh0_194 = buffer.data(lsh0 + 194);
    const auto *lsh0_195 = buffer.data(lsh0 + 195);
    const auto *lsh0_196 = buffer.data(lsh0 + 196);
    const auto *lsh0_197 = buffer.data(lsh0 + 197);
    const auto *lsh0_198 = buffer.data(lsh0 + 198);
    const auto *lsh0_203 = buffer.data(lsh0 + 203);
    const auto *lsh0_204 = buffer.data(lsh0 + 204);
    const auto *lsh0_205 = buffer.data(lsh0 + 205);
    const auto *lsh0_206 = buffer.data(lsh0 + 206);
    const auto *lsh0_207 = buffer.data(lsh0 + 207);
    const auto *lsh0_208 = buffer.data(lsh0 + 208);
    const auto *lsh0_209 = buffer.data(lsh0 + 209);
    const auto *lsh0_210 = buffer.data(lsh0 + 210);
    const auto *lsh0_213 = buffer.data(lsh0 + 213);
    const auto *lsh0_216 = buffer.data(lsh0 + 216);

    const auto *lsh1_146 = buffer.data(lsh1 + 146);
    const auto *lsh1_152 = buffer.data(lsh1 + 152);
    const auto *lsh1_156 = buffer.data(lsh1 + 156);
    const auto *lsh1_161 = buffer.data(lsh1 + 161);
    const auto *lsh1_164 = buffer.data(lsh1 + 164);
    const auto *lsh1_165 = buffer.data(lsh1 + 165);
    const auto *lsh1_166 = buffer.data(lsh1 + 166);
    const auto *lsh1_167 = buffer.data(lsh1 + 167);
    const auto *lsh1_183 = buffer.data(lsh1 + 183);
    const auto *lsh1_185 = buffer.data(lsh1 + 185);
    const auto *lsh1_186 = buffer.data(lsh1 + 186);
    const auto *lsh1_187 = buffer.data(lsh1 + 187);
    const auto *lsh1_188 = buffer.data(lsh1 + 188);
    const auto *lsh1_189 = buffer.data(lsh1 + 189);
    const auto *lsh1_190 = buffer.data(lsh1 + 190);
    const auto *lsh1_191 = buffer.data(lsh1 + 191);
    const auto *lsh1_192 = buffer.data(lsh1 + 192);
    const auto *lsh1_193 = buffer.data(lsh1 + 193);
    const auto *lsh1_194 = buffer.data(lsh1 + 194);
    const auto *lsh1_195 = buffer.data(lsh1 + 195);
    const auto *lsh1_196 = buffer.data(lsh1 + 196);
    const auto *lsh1_197 = buffer.data(lsh1 + 197);
    const auto *lsh1_198 = buffer.data(lsh1 + 198);
    const auto *lsh1_203 = buffer.data(lsh1 + 203);
    const auto *lsh1_204 = buffer.data(lsh1 + 204);
    const auto *lsh1_205 = buffer.data(lsh1 + 205);
    const auto *lsh1_206 = buffer.data(lsh1 + 206);
    const auto *lsh1_207 = buffer.data(lsh1 + 207);
    const auto *lsh1_208 = buffer.data(lsh1 + 208);
    const auto *lsh1_209 = buffer.data(lsh1 + 209);
    const auto *lsh1_210 = buffer.data(lsh1 + 210);
    const auto *lsh1_213 = buffer.data(lsh1 + 213);
    const auto *lsh1_216 = buffer.data(lsh1 + 216);

    const auto *lsi_195 = buffer.data(lsi + 195);
    const auto *lsi_196 = buffer.data(lsi + 196);
    const auto *lsi_198 = buffer.data(lsi + 198);
    const auto *lsi_199 = buffer.data(lsi + 199);
    const auto *lsi_201 = buffer.data(lsi + 201);
    const auto *lsi_202 = buffer.data(lsi + 202);
    const auto *lsi_205 = buffer.data(lsi + 205);
    const auto *lsi_206 = buffer.data(lsi + 206);
    const auto *lsi_210 = buffer.data(lsi + 210);
    const auto *lsi_216 = buffer.data(lsi + 216);
    const auto *lsi_217 = buffer.data(lsi + 217);
    const auto *lsi_218 = buffer.data(lsi + 218);
    const auto *lsi_219 = buffer.data(lsi + 219);
    const auto *lsi_220 = buffer.data(lsi + 220);
    const auto *lsi_221 = buffer.data(lsi + 221);
    const auto *lsi_222 = buffer.data(lsi + 222);
    const auto *lsi_223 = buffer.data(lsi + 223);
    const auto *lsi_224 = buffer.data(lsi + 224);
    const auto *lsi_226 = buffer.data(lsi + 226);
    const auto *lsi_227 = buffer.data(lsi + 227);
    const auto *lsi_229 = buffer.data(lsi + 229);
    const auto *lsi_230 = buffer.data(lsi + 230);
    const auto *lsi_233 = buffer.data(lsi + 233);
    const auto *lsi_234 = buffer.data(lsi + 234);
    const auto *lsi_238 = buffer.data(lsi + 238);
    const auto *lsi_245 = buffer.data(lsi + 245);
    const auto *lsi_246 = buffer.data(lsi + 246);
    const auto *lsi_247 = buffer.data(lsi + 247);
    const auto *lsi_248 = buffer.data(lsi + 248);
    const auto *lsi_249 = buffer.data(lsi + 249);
    const auto *lsi_250 = buffer.data(lsi + 250);
    const auto *lsi_251 = buffer.data(lsi + 251);
    const auto *lsi_252 = buffer.data(lsi + 252);
    const auto *lsi_253 = buffer.data(lsi + 253);
    const auto *lsi_254 = buffer.data(lsi + 254);
    const auto *lsi_255 = buffer.data(lsi + 255);
    const auto *lsi_256 = buffer.data(lsi + 256);
    const auto *lsi_257 = buffer.data(lsi + 257);
    const auto *lsi_258 = buffer.data(lsi + 258);
    const auto *lsi_259 = buffer.data(lsi + 259);
    const auto *lsi_260 = buffer.data(lsi + 260);
    const auto *lsi_261 = buffer.data(lsi + 261);
    const auto *lsi_262 = buffer.data(lsi + 262);
    const auto *lsi_263 = buffer.data(lsi + 263);
    const auto *lsi_264 = buffer.data(lsi + 264);
    const auto *lsi_265 = buffer.data(lsi + 265);
    const auto *lsi_266 = buffer.data(lsi + 266);
    const auto *lsi_272 = buffer.data(lsi + 272);
    const auto *lsi_273 = buffer.data(lsi + 273);
    const auto *lsi_274 = buffer.data(lsi + 274);
    const auto *lsi_275 = buffer.data(lsi + 275);
    const auto *lsi_276 = buffer.data(lsi + 276);
    const auto *lsi_277 = buffer.data(lsi + 277);
    const auto *lsi_278 = buffer.data(lsi + 278);
    const auto *lsi_279 = buffer.data(lsi + 279);
    const auto *lsi_280 = buffer.data(lsi + 280);
    const auto *lsi_281 = buffer.data(lsi + 281);
    const auto *lsi_282 = buffer.data(lsi + 282);
    const auto *lsi_283 = buffer.data(lsi + 283);
    const auto *lsi_286 = buffer.data(lsi + 286);

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_z, pc_y, pc_z, ksk0_108, ksi_111, \
                         ksi_112, ksk1_108, lsh0_146, lsh1_146, lsi_195, \
                         lsi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * ksi_111[k]
                   + f_3 * pc_y[k] * lsi_195[k];

        t_251[k] = f_1 * lsh0_146[k]
                   - f_2 * lsh1_146[k]
                   + f_3 * pc_z[k] * lsi_195[k];

        t_252[k] = pa_z[k] * ksk0_108[k]
                   - f_12 * pc_z[k] * ksk1_108[k];

        t_253[k] = f_14 * ksi_112[k]
                   + f_3 * pc_y[k] * lsi_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pa_z, pc_y, pc_z, ksk0_111, ksi_84, ksi_114, \
                         ksk1_111, lsi_196, lsi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * ksi_84[k]
                   + f_3 * pc_z[k] * lsi_196[k];

        t_255[k] = pa_z[k] * ksk0_111[k]
                   - f_12 * pc_z[k] * ksk1_111[k];

        t_256[k] = f_14 * ksi_114[k]
                   + f_3 * pc_y[k] * lsi_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_z, pc_x, pc_z, ksk0_114, ksi_87, ksi_201, \
                         ksk1_114, lsh0_152, lsh1_152, lsi_199, \
                         lsi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_17 * ksi_201[k]
                   + f_10 * lsh0_152[k]
                   - f_11 * lsh1_152[k]
                   + f_3 * pc_x[k] * lsi_201[k];

        t_258[k] = pa_z[k] * ksk0_114[k]
                   - f_12 * pc_z[k] * ksk1_114[k];

        t_259[k] = f_13 * ksi_87[k]
                   + f_3 * pc_z[k] * lsi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_z, pc_x, pc_y, pc_z, ksk0_118, ksi_117, \
                         ksi_205, ksk1_118, lsh0_156, lsh1_156, lsi_201, \
                         lsi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * ksi_117[k]
                   + f_3 * pc_y[k] * lsi_201[k];

        t_261[k] = f_17 * ksi_205[k]
                   + f_8 * lsh0_156[k]
                   - f_9 * lsh1_156[k]
                   + f_3 * pc_x[k] * lsi_205[k];

        t_262[k] = pa_z[k] * ksk0_118[k]
                   - f_12 * pc_z[k] * ksk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pa_z, pc_y, pc_z, ksk0_120, ksi_90, ksi_91, \
                         ksi_121, ksk1_120, lsi_202, lsi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * ksi_90[k]
                   + f_3 * pc_z[k] * lsi_202[k];

        t_264[k] = pa_z[k] * ksk0_120[k]
                   + f_14 * ksi_91[k]
                   - f_12 * pc_z[k] * ksk1_120[k];

        t_265[k] = f_14 * ksi_121[k]
                   + f_3 * pc_y[k] * lsi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_z, pc_x, pc_z, ksk0_123, ksi_94, ksi_210, \
                         ksk1_123, lsh0_161, lsh1_161, lsi_206, \
                         lsi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_17 * ksi_210[k]
                   + f_6 * lsh0_161[k]
                   - f_7 * lsh1_161[k]
                   + f_3 * pc_x[k] * lsi_210[k];

        t_267[k] = pa_z[k] * ksk0_123[k]
                   - f_12 * pc_z[k] * ksk1_123[k];

        t_268[k] = f_13 * ksi_94[k]
                   + f_3 * pc_z[k] * lsi_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_z, pc_y, pc_z, ksk0_125, ksk0_126, ksi_95, \
                         ksi_96, ksi_126, ksk1_125, ksk1_126, lsi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_z[k] * ksk0_125[k]
                   + f_14 * ksi_95[k]
                   - f_12 * pc_z[k] * ksk1_125[k];

        t_270[k] = pa_z[k] * ksk0_126[k]
                   + f_15 * ksi_96[k]
                   - f_12 * pc_z[k] * ksk1_126[k];

        t_271[k] = f_14 * ksi_126[k]
                   + f_3 * pc_y[k] * lsi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, ksi_216, ksi_217, ksi_218, ksi_219, \
                         lsh0_167, lsh1_167, lsi_216, lsi_217, lsi_218, \
                         lsi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_17 * ksi_216[k]
                   + f_4 * lsh0_167[k]
                   - f_5 * lsh1_167[k]
                   + f_3 * pc_x[k] * lsi_216[k];

        t_273[k] = f_17 * ksi_217[k]
                   + f_3 * pc_x[k] * lsi_217[k];

        t_274[k] = f_17 * ksi_218[k]
                   + f_3 * pc_x[k] * lsi_218[k];

        t_275[k] = f_17 * ksi_219[k]
                   + f_3 * pc_x[k] * lsi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, ksi_220, ksi_221, ksi_222, ksi_223, \
                         lsi_220, lsi_221, lsi_222, lsi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_17 * ksi_220[k]
                   + f_3 * pc_x[k] * lsi_220[k];

        t_277[k] = f_17 * ksi_221[k]
                   + f_3 * pc_x[k] * lsi_221[k];

        t_278[k] = f_17 * ksi_222[k]
                   + f_3 * pc_x[k] * lsi_222[k];

        t_279[k] = f_17 * ksi_223[k]
                   + f_3 * pc_x[k] * lsi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_z, pc_y, pc_z, ksk0_136, ksi_105, ksi_135, \
                         ksk1_136, lsh0_164, lsh1_164, lsi_217, \
                         lsi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * ksk0_136[k]
                   - f_12 * pc_z[k] * ksk1_136[k];

        t_281[k] = f_13 * ksi_105[k]
                   + f_3 * pc_z[k] * lsi_217[k];

        t_282[k] = f_14 * ksi_135[k]
                   + f_10 * lsh0_164[k]
                   - f_11 * lsh1_164[k]
                   + f_3 * pc_y[k] * lsi_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, ksi_136, ksi_137, ksi_138, lsh0_165, \
                         lsh0_166, lsh0_167, lsh1_165, lsh1_166, lsh1_167, lsi_220, lsi_221, \
                         lsi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * ksi_136[k]
                   + f_8 * lsh0_165[k]
                   - f_9 * lsh1_165[k]
                   + f_3 * pc_y[k] * lsi_220[k];

        t_284[k] = f_14 * ksi_137[k]
                   + f_6 * lsh0_166[k]
                   - f_7 * lsh1_166[k]
                   + f_3 * pc_y[k] * lsi_221[k];

        t_285[k] = f_14 * ksi_138[k]
                   + f_4 * lsh0_167[k]
                   - f_5 * lsh1_167[k]
                   + f_3 * pc_y[k] * lsi_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pc_y, pc_z, ksk0_180, ksi_111, \
                         ksi_139, ksi_140, ksk1_180, lsh0_167, lsh1_167, lsi_223, \
                         lsi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * ksi_139[k]
                   + f_3 * pc_y[k] * lsi_223[k];

        t_287[k] = f_13 * ksi_111[k]
                   + f_1 * lsh0_167[k]
                   - f_2 * lsh1_167[k]
                   + f_3 * pc_z[k] * lsi_223[k];

        t_288[k] = pa_y[k] * ksk0_180[k]
                   - f_12 * pc_y[k] * ksk1_180[k];

        t_289[k] = f_13 * ksi_140[k]
                   + f_3 * pc_y[k] * lsi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pc_y, pc_z, ksk0_183, ksk0_185, \
                         ksi_112, ksi_141, ksi_142, ksk1_183, ksk1_185, lsi_224, \
                         lsi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * ksi_112[k]
                   + f_3 * pc_z[k] * lsi_224[k];

        t_291[k] = pa_y[k] * ksk0_183[k]
                   + f_14 * ksi_141[k]
                   - f_12 * pc_y[k] * ksk1_183[k];

        t_292[k] = f_13 * ksi_142[k]
                   + f_3 * pc_y[k] * lsi_226[k];

        t_293[k] = pa_y[k] * ksk0_185[k]
                   - f_12 * pc_y[k] * ksk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pc_y, pc_z, ksk0_186, ksk0_189, \
                         ksi_115, ksi_143, ksi_145, ksk1_186, ksk1_189, lsi_227, \
                         lsi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pa_y[k] * ksk0_186[k]
                   + f_15 * ksi_143[k]
                   - f_12 * pc_y[k] * ksk1_186[k];

        t_295[k] = f_14 * ksi_115[k]
                   + f_3 * pc_z[k] * lsi_227[k];

        t_296[k] = f_13 * ksi_145[k]
                   + f_3 * pc_y[k] * lsi_229[k];

        t_297[k] = pa_y[k] * ksk0_189[k]
                   - f_12 * pc_y[k] * ksk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pa_y, pc_y, pc_z, ksk0_190, ksk0_192, ksi_118, \
                         ksi_146, ksi_148, ksk1_190, ksk1_192, \
                         lsi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pa_y[k] * ksk0_190[k]
                   + f_16 * ksi_146[k]
                   - f_12 * pc_y[k] * ksk1_190[k];

        t_299[k] = f_14 * ksi_118[k]
                   + f_3 * pc_z[k] * lsi_230[k];

        t_300[k] = pa_y[k] * ksk0_192[k]
                   + f_14 * ksi_148[k]
                   - f_12 * pc_y[k] * ksk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_y, pc_y, pc_z, ksk0_194, ksk0_195, \
                         ksi_122, ksi_149, ksi_150, ksk1_194, ksk1_195, lsi_233, \
                         lsi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * ksi_149[k]
                   + f_3 * pc_y[k] * lsi_233[k];

        t_302[k] = pa_y[k] * ksk0_194[k]
                   - f_12 * pc_y[k] * ksk1_194[k];

        t_303[k] = pa_y[k] * ksk0_195[k]
                   + f_17 * ksi_150[k]
                   - f_12 * pc_y[k] * ksk1_195[k];

        t_304[k] = f_14 * ksi_122[k]
                   + f_3 * pc_z[k] * lsi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_y, pc_y, ksk0_197, ksk0_198, ksk0_200, \
                         ksi_152, ksi_153, ksi_154, ksk1_197, ksk1_198, ksk1_200, \
                         lsi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_y[k] * ksk0_197[k]
                   + f_15 * ksi_152[k]
                   - f_12 * pc_y[k] * ksk1_197[k];

        t_306[k] = pa_y[k] * ksk0_198[k]
                   + f_14 * ksi_153[k]
                   - f_12 * pc_y[k] * ksk1_198[k];

        t_307[k] = f_13 * ksi_154[k]
                   + f_3 * pc_y[k] * lsi_238[k];

        t_308[k] = pa_y[k] * ksk0_200[k]
                   - f_12 * pc_y[k] * ksk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, ksi_245, ksi_246, ksi_247, \
                         ksi_248, ksi_249, lsi_245, lsi_246, lsi_247, lsi_248, \
                         lsi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_17 * ksi_245[k]
                   + f_3 * pc_x[k] * lsi_245[k];

        t_310[k] = f_17 * ksi_246[k]
                   + f_3 * pc_x[k] * lsi_246[k];

        t_311[k] = f_17 * ksi_247[k]
                   + f_3 * pc_x[k] * lsi_247[k];

        t_312[k] = f_17 * ksi_248[k]
                   + f_3 * pc_x[k] * lsi_248[k];

        t_313[k] = f_17 * ksi_249[k]
                   + f_3 * pc_x[k] * lsi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, ksi_133, ksi_161, \
                         ksi_250, ksi_251, lsh0_183, lsh1_183, lsi_245, lsi_250, \
                         lsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_17 * ksi_250[k]
                   + f_3 * pc_x[k] * lsi_250[k];

        t_315[k] = f_17 * ksi_251[k]
                   + f_3 * pc_x[k] * lsi_251[k];

        t_316[k] = f_13 * ksi_161[k]
                   + f_1 * lsh0_183[k]
                   - f_2 * lsh1_183[k]
                   + f_3 * pc_y[k] * lsi_245[k];

        t_317[k] = f_14 * ksi_133[k]
                   + f_3 * pc_z[k] * lsi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, ksi_163, ksi_164, ksi_165, lsh0_185, \
                         lsh0_186, lsh0_187, lsh1_185, lsh1_186, lsh1_187, lsi_247, lsi_248, \
                         lsi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * ksi_163[k]
                   + f_10 * lsh0_185[k]
                   - f_11 * lsh1_185[k]
                   + f_3 * pc_y[k] * lsi_247[k];

        t_319[k] = f_13 * ksi_164[k]
                   + f_8 * lsh0_186[k]
                   - f_9 * lsh1_186[k]
                   + f_3 * pc_y[k] * lsi_248[k];

        t_320[k] = f_13 * ksi_165[k]
                   + f_6 * lsh0_187[k]
                   - f_7 * lsh1_187[k]
                   + f_3 * pc_y[k] * lsi_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_y, pc_y, ksk0_215, ksi_166, ksi_167, \
                         ksk1_215, lsh0_188, lsh1_188, lsi_250, \
                         lsi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * ksi_166[k]
                   + f_4 * lsh0_188[k]
                   - f_5 * lsh1_188[k]
                   + f_3 * pc_y[k] * lsi_250[k];

        t_322[k] = f_13 * ksi_167[k]
                   + f_3 * pc_y[k] * lsi_251[k];

        t_323[k] = pa_y[k] * ksk0_215[k]
                   - f_12 * pc_y[k] * ksk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, ksi_140, \
                         ksi_252, lsh0_189, lsh1_189, lsi_252, lsi_253, \
                         lsi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_17 * ksi_252[k]
                   + f_1 * lsh0_189[k]
                   - f_2 * lsh1_189[k]
                   + f_3 * pc_x[k] * lsi_252[k];

        t_325[k] = f_3 * pc_y[k] * lsi_252[k];

        t_326[k] = f_15 * ksi_140[k]
                   + f_3 * pc_z[k] * lsi_252[k];

        t_327[k] = f_4 * lsh0_189[k]
                   - f_5 * lsh1_189[k]
                   + f_3 * pc_y[k] * lsi_253[k];

        t_328[k] = f_3 * pc_y[k] * lsi_254[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pc_x, pc_y, ksi_257, lsh0_190, lsh0_191, \
                         lsh0_194, lsh1_190, lsh1_191, lsh1_194, lsi_255, lsi_256, \
                         lsi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_17 * ksi_257[k]
                   + f_10 * lsh0_194[k]
                   - f_11 * lsh1_194[k]
                   + f_3 * pc_x[k] * lsi_257[k];

        t_330[k] = f_6 * lsh0_190[k]
                   - f_7 * lsh1_190[k]
                   + f_3 * pc_y[k] * lsi_255[k];

        t_331[k] = f_4 * lsh0_191[k]
                   - f_5 * lsh1_191[k]
                   + f_3 * pc_y[k] * lsi_256[k];

        t_332[k] = f_3 * pc_y[k] * lsi_257[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_x, pc_y, ksi_261, lsh0_192, lsh0_193, \
                         lsh0_198, lsh1_192, lsh1_193, lsh1_198, lsi_258, lsi_259, \
                         lsi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_17 * ksi_261[k]
                   + f_8 * lsh0_198[k]
                   - f_9 * lsh1_198[k]
                   + f_3 * pc_x[k] * lsi_261[k];

        t_334[k] = f_8 * lsh0_192[k]
                   - f_9 * lsh1_192[k]
                   + f_3 * pc_y[k] * lsi_258[k];

        t_335[k] = f_6 * lsh0_193[k]
                   - f_7 * lsh1_193[k]
                   + f_3 * pc_y[k] * lsi_259[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, ksi_266, lsh0_194, lsh0_203, \
                         lsh1_194, lsh1_203, lsi_260, lsi_261, \
                         lsi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_4 * lsh0_194[k]
                   - f_5 * lsh1_194[k]
                   + f_3 * pc_y[k] * lsi_260[k];

        t_337[k] = f_3 * pc_y[k] * lsi_261[k];

        t_338[k] = f_17 * ksi_266[k]
                   + f_6 * lsh0_203[k]
                   - f_7 * lsh1_203[k]
                   + f_3 * pc_x[k] * lsi_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_y, lsh0_195, lsh0_196, lsh0_197, lsh1_195, \
                         lsh1_196, lsh1_197, lsi_262, lsi_263, \
                         lsi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_10 * lsh0_195[k]
                   - f_11 * lsh1_195[k]
                   + f_3 * pc_y[k] * lsi_262[k];

        t_340[k] = f_8 * lsh0_196[k]
                   - f_9 * lsh1_196[k]
                   + f_3 * pc_y[k] * lsi_263[k];

        t_341[k] = f_6 * lsh0_197[k]
                   - f_7 * lsh1_197[k]
                   + f_3 * pc_y[k] * lsi_264[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_y, ksi_272, ksi_273, lsh0_198, \
                         lsh0_209, lsh1_198, lsh1_209, lsi_265, lsi_266, lsi_272, \
                         lsi_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_4 * lsh0_198[k]
                   - f_5 * lsh1_198[k]
                   + f_3 * pc_y[k] * lsi_265[k];

        t_343[k] = f_3 * pc_y[k] * lsi_266[k];

        t_344[k] = f_17 * ksi_272[k]
                   + f_4 * lsh0_209[k]
                   - f_5 * lsh1_209[k]
                   + f_3 * pc_x[k] * lsi_272[k];

        t_345[k] = f_17 * ksi_273[k]
                   + f_3 * pc_x[k] * lsi_273[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, pc_x, pc_y, ksi_274, ksi_275, \
                         ksi_276, ksi_277, lsi_272, lsi_274, lsi_275, lsi_276, \
                         lsi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_17 * ksi_274[k]
                   + f_3 * pc_x[k] * lsi_274[k];

        t_347[k] = f_17 * ksi_275[k]
                   + f_3 * pc_x[k] * lsi_275[k];

        t_348[k] = f_17 * ksi_276[k]
                   + f_3 * pc_x[k] * lsi_276[k];

        t_349[k] = f_17 * ksi_277[k]
                   + f_3 * pc_x[k] * lsi_277[k];

        t_350[k] = f_3 * pc_y[k] * lsi_272[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pc_x, pc_y, ksi_279, lsh0_204, lsh0_205, \
                         lsh1_204, lsh1_205, lsi_273, lsi_274, \
                         lsi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_17 * ksi_279[k]
                   + f_3 * pc_x[k] * lsi_279[k];

        t_352[k] = f_1 * lsh0_204[k]
                   - f_2 * lsh1_204[k]
                   + f_3 * pc_y[k] * lsi_273[k];

        t_353[k] = f_19 * lsh0_205[k]
                   - f_20 * lsh1_205[k]
                   + f_3 * pc_y[k] * lsi_274[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_y, lsh0_206, lsh0_207, lsh0_208, lsh1_206, \
                         lsh1_207, lsh1_208, lsi_275, lsi_276, \
                         lsi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_10 * lsh0_206[k]
                   - f_11 * lsh1_206[k]
                   + f_3 * pc_y[k] * lsi_275[k];

        t_355[k] = f_8 * lsh0_207[k]
                   - f_9 * lsh1_207[k]
                   + f_3 * pc_y[k] * lsi_276[k];

        t_356[k] = f_6 * lsh0_208[k]
                   - f_7 * lsh1_208[k]
                   + f_3 * pc_y[k] * lsi_277[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pc_x, pc_y, pc_z, ksi_167, ksi_280, \
                         lsh0_209, lsh0_210, lsh1_209, lsh1_210, lsi_278, lsi_279, \
                         lsi_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_4 * lsh0_209[k]
                   - f_5 * lsh1_209[k]
                   + f_3 * pc_y[k] * lsi_278[k];

        t_358[k] = f_3 * pc_y[k] * lsi_279[k];

        t_359[k] = f_15 * ksi_167[k]
                   + f_1 * lsh0_209[k]
                   - f_2 * lsh1_209[k]
                   + f_3 * pc_z[k] * lsi_279[k];

        t_360[k] = f_16 * ksi_280[k]
                   + f_1 * lsh0_210[k]
                   - f_2 * lsh1_210[k]
                   + f_3 * pc_x[k] * lsi_280[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, ksi_168, ksi_283, \
                         lsh0_213, lsh1_213, lsi_280, lsi_281, \
                         lsi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_16 * ksi_168[k]
                   + f_3 * pc_y[k] * lsi_280[k];

        t_362[k] = f_3 * pc_z[k] * lsi_280[k];

        t_363[k] = f_16 * ksi_283[k]
                   + f_10 * lsh0_213[k]
                   - f_11 * lsh1_213[k]
                   + f_3 * pc_x[k] * lsi_283[k];

        t_364[k] = f_3 * pc_z[k] * lsi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, pc_x, pc_z, ksi_286, lsh0_210, lsh0_216, \
                         lsh1_210, lsh1_216, lsi_282, lsi_283, \
                         lsi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_4 * lsh0_210[k]
                   - f_5 * lsh1_210[k]
                   + f_3 * pc_z[k] * lsi_282[k];

        t_366[k] = f_16 * ksi_286[k]
                   + f_8 * lsh0_216[k]
                   - f_9 * lsh1_216[k]
                   + f_3 * pc_x[k] * lsi_286[k];

        t_367[k] = f_3 * pc_z[k] * lsi_283[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsh0, const size_t lsh1,
                                                          const size_t lsi, const size_t ncols,
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
    const auto f_16 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_216 = buffer.data(ksk0 + 216);
    const auto *ksk0_219 = buffer.data(ksk0 + 219);
    const auto *ksk0_222 = buffer.data(ksk0 + 222);
    const auto *ksk0_226 = buffer.data(ksk0 + 226);
    const auto *ksk0_228 = buffer.data(ksk0 + 228);
    const auto *ksk0_231 = buffer.data(ksk0 + 231);
    const auto *ksk0_233 = buffer.data(ksk0 + 233);
    const auto *ksk0_234 = buffer.data(ksk0 + 234);
    const auto *ksk0_244 = buffer.data(ksk0 + 244);
    const auto *ksk0_324 = buffer.data(ksk0 + 324);
    const auto *ksk0_327 = buffer.data(ksk0 + 327);
    const auto *ksk0_329 = buffer.data(ksk0 + 329);
    const auto *ksk0_330 = buffer.data(ksk0 + 330);
    const auto *ksk0_333 = buffer.data(ksk0 + 333);
    const auto *ksk0_334 = buffer.data(ksk0 + 334);
    const auto *ksk0_336 = buffer.data(ksk0 + 336);

    const auto *ksi_168 = buffer.data(ksi + 168);
    const auto *ksi_171 = buffer.data(ksi + 171);
    const auto *ksi_173 = buffer.data(ksi + 173);
    const auto *ksi_174 = buffer.data(ksi + 174);
    const auto *ksi_175 = buffer.data(ksi + 175);
    const auto *ksi_177 = buffer.data(ksi + 177);
    const auto *ksi_178 = buffer.data(ksi + 178);
    const auto *ksi_179 = buffer.data(ksi + 179);
    const auto *ksi_180 = buffer.data(ksi + 180);
    const auto *ksi_182 = buffer.data(ksi + 182);
    const auto *ksi_189 = buffer.data(ksi + 189);
    const auto *ksi_195 = buffer.data(ksi + 195);
    const auto *ksi_196 = buffer.data(ksi + 196);
    const auto *ksi_198 = buffer.data(ksi + 198);
    const auto *ksi_199 = buffer.data(ksi + 199);
    const auto *ksi_201 = buffer.data(ksi + 201);
    const auto *ksi_202 = buffer.data(ksi + 202);
    const auto *ksi_205 = buffer.data(ksi + 205);
    const auto *ksi_206 = buffer.data(ksi + 206);
    const auto *ksi_210 = buffer.data(ksi + 210);
    const auto *ksi_217 = buffer.data(ksi + 217);
    const auto *ksi_219 = buffer.data(ksi + 219);
    const auto *ksi_220 = buffer.data(ksi + 220);
    const auto *ksi_221 = buffer.data(ksi + 221);
    const auto *ksi_222 = buffer.data(ksi + 222);
    const auto *ksi_223 = buffer.data(ksi + 223);
    const auto *ksi_224 = buffer.data(ksi + 224);
    const auto *ksi_226 = buffer.data(ksi + 226);
    const auto *ksi_227 = buffer.data(ksi + 227);
    const auto *ksi_229 = buffer.data(ksi + 229);
    const auto *ksi_230 = buffer.data(ksi + 230);
    const auto *ksi_233 = buffer.data(ksi + 233);
    const auto *ksi_238 = buffer.data(ksi + 238);
    const auto *ksi_245 = buffer.data(ksi + 245);
    const auto *ksi_247 = buffer.data(ksi + 247);
    const auto *ksi_248 = buffer.data(ksi + 248);
    const auto *ksi_249 = buffer.data(ksi + 249);
    const auto *ksi_250 = buffer.data(ksi + 250);
    const auto *ksi_251 = buffer.data(ksi + 251);
    const auto *ksi_252 = buffer.data(ksi + 252);
    const auto *ksi_253 = buffer.data(ksi + 253);
    const auto *ksi_254 = buffer.data(ksi + 254);
    const auto *ksi_255 = buffer.data(ksi + 255);
    const auto *ksi_257 = buffer.data(ksi + 257);
    const auto *ksi_258 = buffer.data(ksi + 258);
    const auto *ksi_260 = buffer.data(ksi + 260);
    const auto *ksi_290 = buffer.data(ksi + 290);
    const auto *ksi_295 = buffer.data(ksi + 295);
    const auto *ksi_301 = buffer.data(ksi + 301);
    const auto *ksi_303 = buffer.data(ksi + 303);
    const auto *ksi_304 = buffer.data(ksi + 304);
    const auto *ksi_305 = buffer.data(ksi + 305);
    const auto *ksi_306 = buffer.data(ksi + 306);
    const auto *ksi_307 = buffer.data(ksi + 307);
    const auto *ksi_313 = buffer.data(ksi + 313);
    const auto *ksi_317 = buffer.data(ksi + 317);
    const auto *ksi_322 = buffer.data(ksi + 322);
    const auto *ksi_328 = buffer.data(ksi + 328);
    const auto *ksi_329 = buffer.data(ksi + 329);
    const auto *ksi_330 = buffer.data(ksi + 330);
    const auto *ksi_331 = buffer.data(ksi + 331);
    const auto *ksi_332 = buffer.data(ksi + 332);
    const auto *ksi_333 = buffer.data(ksi + 333);
    const auto *ksi_334 = buffer.data(ksi + 334);
    const auto *ksi_335 = buffer.data(ksi + 335);
    const auto *ksi_336 = buffer.data(ksi + 336);
    const auto *ksi_339 = buffer.data(ksi + 339);
    const auto *ksi_341 = buffer.data(ksi + 341);
    const auto *ksi_342 = buffer.data(ksi + 342);
    const auto *ksi_345 = buffer.data(ksi + 345);
    const auto *ksi_346 = buffer.data(ksi + 346);
    const auto *ksi_348 = buffer.data(ksi + 348);
    const auto *ksi_350 = buffer.data(ksi + 350);
    const auto *ksi_351 = buffer.data(ksi + 351);
    const auto *ksi_353 = buffer.data(ksi + 353);
    const auto *ksi_354 = buffer.data(ksi + 354);
    const auto *ksi_356 = buffer.data(ksi + 356);
    const auto *ksi_357 = buffer.data(ksi + 357);
    const auto *ksi_358 = buffer.data(ksi + 358);
    const auto *ksi_359 = buffer.data(ksi + 359);
    const auto *ksi_360 = buffer.data(ksi + 360);
    const auto *ksi_361 = buffer.data(ksi + 361);
    const auto *ksi_362 = buffer.data(ksi + 362);
    const auto *ksi_363 = buffer.data(ksi + 363);

    const auto *ksk1_216 = buffer.data(ksk1 + 216);
    const auto *ksk1_219 = buffer.data(ksk1 + 219);
    const auto *ksk1_222 = buffer.data(ksk1 + 222);
    const auto *ksk1_226 = buffer.data(ksk1 + 226);
    const auto *ksk1_228 = buffer.data(ksk1 + 228);
    const auto *ksk1_231 = buffer.data(ksk1 + 231);
    const auto *ksk1_233 = buffer.data(ksk1 + 233);
    const auto *ksk1_234 = buffer.data(ksk1 + 234);
    const auto *ksk1_244 = buffer.data(ksk1 + 244);
    const auto *ksk1_324 = buffer.data(ksk1 + 324);
    const auto *ksk1_327 = buffer.data(ksk1 + 327);
    const auto *ksk1_329 = buffer.data(ksk1 + 329);
    const auto *ksk1_330 = buffer.data(ksk1 + 330);
    const auto *ksk1_333 = buffer.data(ksk1 + 333);
    const auto *ksk1_334 = buffer.data(ksk1 + 334);
    const auto *ksk1_336 = buffer.data(ksk1 + 336);

    const auto *lsh0_212 = buffer.data(lsh0 + 212);
    const auto *lsh0_213 = buffer.data(lsh0 + 213);
    const auto *lsh0_215 = buffer.data(lsh0 + 215);
    const auto *lsh0_216 = buffer.data(lsh0 + 216);
    const auto *lsh0_217 = buffer.data(lsh0 + 217);
    const auto *lsh0_219 = buffer.data(lsh0 + 219);
    const auto *lsh0_220 = buffer.data(lsh0 + 220);
    const auto *lsh0_225 = buffer.data(lsh0 + 225);
    const auto *lsh0_226 = buffer.data(lsh0 + 226);
    const auto *lsh0_227 = buffer.data(lsh0 + 227);
    const auto *lsh0_228 = buffer.data(lsh0 + 228);
    const auto *lsh0_230 = buffer.data(lsh0 + 230);
    const auto *lsh0_236 = buffer.data(lsh0 + 236);
    const auto *lsh0_240 = buffer.data(lsh0 + 240);
    const auto *lsh0_245 = buffer.data(lsh0 + 245);
    const auto *lsh0_248 = buffer.data(lsh0 + 248);
    const auto *lsh0_249 = buffer.data(lsh0 + 249);
    const auto *lsh0_250 = buffer.data(lsh0 + 250);
    const auto *lsh0_251 = buffer.data(lsh0 + 251);
    const auto *lsh0_252 = buffer.data(lsh0 + 252);
    const auto *lsh0_255 = buffer.data(lsh0 + 255);
    const auto *lsh0_257 = buffer.data(lsh0 + 257);
    const auto *lsh0_258 = buffer.data(lsh0 + 258);
    const auto *lsh0_261 = buffer.data(lsh0 + 261);
    const auto *lsh0_262 = buffer.data(lsh0 + 262);
    const auto *lsh0_264 = buffer.data(lsh0 + 264);
    const auto *lsh0_266 = buffer.data(lsh0 + 266);
    const auto *lsh0_267 = buffer.data(lsh0 + 267);
    const auto *lsh0_269 = buffer.data(lsh0 + 269);
    const auto *lsh0_270 = buffer.data(lsh0 + 270);
    const auto *lsh0_271 = buffer.data(lsh0 + 271);
    const auto *lsh0_272 = buffer.data(lsh0 + 272);

    const auto *lsh1_212 = buffer.data(lsh1 + 212);
    const auto *lsh1_213 = buffer.data(lsh1 + 213);
    const auto *lsh1_215 = buffer.data(lsh1 + 215);
    const auto *lsh1_216 = buffer.data(lsh1 + 216);
    const auto *lsh1_217 = buffer.data(lsh1 + 217);
    const auto *lsh1_219 = buffer.data(lsh1 + 219);
    const auto *lsh1_220 = buffer.data(lsh1 + 220);
    const auto *lsh1_225 = buffer.data(lsh1 + 225);
    const auto *lsh1_226 = buffer.data(lsh1 + 226);
    const auto *lsh1_227 = buffer.data(lsh1 + 227);
    const auto *lsh1_228 = buffer.data(lsh1 + 228);
    const auto *lsh1_230 = buffer.data(lsh1 + 230);
    const auto *lsh1_236 = buffer.data(lsh1 + 236);
    const auto *lsh1_240 = buffer.data(lsh1 + 240);
    const auto *lsh1_245 = buffer.data(lsh1 + 245);
    const auto *lsh1_248 = buffer.data(lsh1 + 248);
    const auto *lsh1_249 = buffer.data(lsh1 + 249);
    const auto *lsh1_250 = buffer.data(lsh1 + 250);
    const auto *lsh1_251 = buffer.data(lsh1 + 251);
    const auto *lsh1_252 = buffer.data(lsh1 + 252);
    const auto *lsh1_255 = buffer.data(lsh1 + 255);
    const auto *lsh1_257 = buffer.data(lsh1 + 257);
    const auto *lsh1_258 = buffer.data(lsh1 + 258);
    const auto *lsh1_261 = buffer.data(lsh1 + 261);
    const auto *lsh1_262 = buffer.data(lsh1 + 262);
    const auto *lsh1_264 = buffer.data(lsh1 + 264);
    const auto *lsh1_266 = buffer.data(lsh1 + 266);
    const auto *lsh1_267 = buffer.data(lsh1 + 267);
    const auto *lsh1_269 = buffer.data(lsh1 + 269);
    const auto *lsh1_270 = buffer.data(lsh1 + 270);
    const auto *lsh1_271 = buffer.data(lsh1 + 271);
    const auto *lsh1_272 = buffer.data(lsh1 + 272);

    const auto *lsi_285 = buffer.data(lsi + 285);
    const auto *lsi_286 = buffer.data(lsi + 286);
    const auto *lsi_287 = buffer.data(lsi + 287);
    const auto *lsi_289 = buffer.data(lsi + 289);
    const auto *lsi_290 = buffer.data(lsi + 290);
    const auto *lsi_291 = buffer.data(lsi + 291);
    const auto *lsi_292 = buffer.data(lsi + 292);
    const auto *lsi_294 = buffer.data(lsi + 294);
    const auto *lsi_295 = buffer.data(lsi + 295);
    const auto *lsi_301 = buffer.data(lsi + 301);
    const auto *lsi_302 = buffer.data(lsi + 302);
    const auto *lsi_303 = buffer.data(lsi + 303);
    const auto *lsi_304 = buffer.data(lsi + 304);
    const auto *lsi_305 = buffer.data(lsi + 305);
    const auto *lsi_306 = buffer.data(lsi + 306);
    const auto *lsi_307 = buffer.data(lsi + 307);
    const auto *lsi_308 = buffer.data(lsi + 308);
    const auto *lsi_310 = buffer.data(lsi + 310);
    const auto *lsi_311 = buffer.data(lsi + 311);
    const auto *lsi_313 = buffer.data(lsi + 313);
    const auto *lsi_314 = buffer.data(lsi + 314);
    const auto *lsi_317 = buffer.data(lsi + 317);
    const auto *lsi_318 = buffer.data(lsi + 318);
    const auto *lsi_322 = buffer.data(lsi + 322);
    const auto *lsi_328 = buffer.data(lsi + 328);
    const auto *lsi_329 = buffer.data(lsi + 329);
    const auto *lsi_330 = buffer.data(lsi + 330);
    const auto *lsi_331 = buffer.data(lsi + 331);
    const auto *lsi_332 = buffer.data(lsi + 332);
    const auto *lsi_333 = buffer.data(lsi + 333);
    const auto *lsi_334 = buffer.data(lsi + 334);
    const auto *lsi_335 = buffer.data(lsi + 335);
    const auto *lsi_336 = buffer.data(lsi + 336);
    const auto *lsi_338 = buffer.data(lsi + 338);
    const auto *lsi_339 = buffer.data(lsi + 339);
    const auto *lsi_341 = buffer.data(lsi + 341);
    const auto *lsi_342 = buffer.data(lsi + 342);
    const auto *lsi_345 = buffer.data(lsi + 345);
    const auto *lsi_346 = buffer.data(lsi + 346);
    const auto *lsi_348 = buffer.data(lsi + 348);
    const auto *lsi_350 = buffer.data(lsi + 350);
    const auto *lsi_351 = buffer.data(lsi + 351);
    const auto *lsi_353 = buffer.data(lsi + 353);
    const auto *lsi_354 = buffer.data(lsi + 354);
    const auto *lsi_356 = buffer.data(lsi + 356);
    const auto *lsi_357 = buffer.data(lsi + 357);
    const auto *lsi_358 = buffer.data(lsi + 358);
    const auto *lsi_359 = buffer.data(lsi + 359);
    const auto *lsi_360 = buffer.data(lsi + 360);
    const auto *lsi_361 = buffer.data(lsi + 361);
    const auto *lsi_362 = buffer.data(lsi + 362);
    const auto *lsi_363 = buffer.data(lsi + 363);
    const auto *lsi_364 = buffer.data(lsi + 364);
    const auto *lsi_366 = buffer.data(lsi + 366);
    const auto *lsi_367 = buffer.data(lsi + 367);
    const auto *lsi_369 = buffer.data(lsi + 369);
    const auto *lsi_370 = buffer.data(lsi + 370);

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pc_x, pc_y, pc_z, ksi_173, ksi_290, \
                         lsh0_212, lsh0_220, lsh1_212, lsh1_220, lsi_285, lsi_286, \
                         lsi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_16 * ksi_173[k]
                   + f_3 * pc_y[k] * lsi_285[k];

        t_369[k] = f_6 * lsh0_212[k]
                   - f_7 * lsh1_212[k]
                   + f_3 * pc_z[k] * lsi_285[k];

        t_370[k] = f_16 * ksi_290[k]
                   + f_6 * lsh0_220[k]
                   - f_7 * lsh1_220[k]
                   + f_3 * pc_x[k] * lsi_290[k];

        t_371[k] = f_3 * pc_z[k] * lsi_286[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, ksi_177, lsh0_213, lsh0_215, \
                         lsh1_213, lsh1_215, lsi_287, lsi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_4 * lsh0_213[k]
                   - f_5 * lsh1_213[k]
                   + f_3 * pc_z[k] * lsi_287[k];

        t_373[k] = f_16 * ksi_177[k]
                   + f_3 * pc_y[k] * lsi_289[k];

        t_374[k] = f_8 * lsh0_215[k]
                   - f_9 * lsh1_215[k]
                   + f_3 * pc_z[k] * lsi_289[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_z, ksi_295, lsh0_216, lsh0_225, \
                         lsh1_216, lsh1_225, lsi_290, lsi_291, \
                         lsi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_16 * ksi_295[k]
                   + f_4 * lsh0_225[k]
                   - f_5 * lsh1_225[k]
                   + f_3 * pc_x[k] * lsi_295[k];

        t_376[k] = f_3 * pc_z[k] * lsi_290[k];

        t_377[k] = f_4 * lsh0_216[k]
                   - f_5 * lsh1_216[k]
                   + f_3 * pc_z[k] * lsi_291[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, pc_y, pc_z, ksi_182, ksi_301, \
                         lsh0_217, lsh0_219, lsh1_217, lsh1_219, lsi_292, lsi_294, \
                         lsi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_6 * lsh0_217[k]
                   - f_7 * lsh1_217[k]
                   + f_3 * pc_z[k] * lsi_292[k];

        t_379[k] = f_16 * ksi_182[k]
                   + f_3 * pc_y[k] * lsi_294[k];

        t_380[k] = f_10 * lsh0_219[k]
                   - f_11 * lsh1_219[k]
                   + f_3 * pc_z[k] * lsi_294[k];

        t_381[k] = f_16 * ksi_301[k]
                   + f_3 * pc_x[k] * lsi_301[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_z, ksi_303, ksi_304, \
                         ksi_305, ksi_306, lsi_295, lsi_303, lsi_304, lsi_305, \
                         lsi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_3 * pc_z[k] * lsi_295[k];

        t_383[k] = f_16 * ksi_303[k]
                   + f_3 * pc_x[k] * lsi_303[k];

        t_384[k] = f_16 * ksi_304[k]
                   + f_3 * pc_x[k] * lsi_304[k];

        t_385[k] = f_16 * ksi_305[k]
                   + f_3 * pc_x[k] * lsi_305[k];

        t_386[k] = f_16 * ksi_306[k]
                   + f_3 * pc_x[k] * lsi_306[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pc_x, pc_y, pc_z, ksi_189, ksi_307, \
                         lsh0_225, lsh1_225, lsi_301, lsi_302, \
                         lsi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_16 * ksi_307[k]
                   + f_3 * pc_x[k] * lsi_307[k];

        t_388[k] = f_16 * ksi_189[k]
                   + f_1 * lsh0_225[k]
                   - f_2 * lsh1_225[k]
                   + f_3 * pc_y[k] * lsi_301[k];

        t_389[k] = f_3 * pc_z[k] * lsi_301[k];

        t_390[k] = f_4 * lsh0_225[k]
                   - f_5 * lsh1_225[k]
                   + f_3 * pc_z[k] * lsi_302[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_z, lsh0_226, lsh0_227, lsh0_228, lsh1_226, \
                         lsh1_227, lsh1_228, lsi_303, lsi_304, \
                         lsi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_6 * lsh0_226[k]
                   - f_7 * lsh1_226[k]
                   + f_3 * pc_z[k] * lsi_303[k];

        t_392[k] = f_8 * lsh0_227[k]
                   - f_9 * lsh1_227[k]
                   + f_3 * pc_z[k] * lsi_304[k];

        t_393[k] = f_10 * lsh0_228[k]
                   - f_11 * lsh1_228[k]
                   + f_3 * pc_z[k] * lsi_305[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pa_z, pc_y, pc_z, ksk0_216, ksi_195, \
                         ksi_196, ksk1_216, lsh0_230, lsh1_230, lsi_307, \
                         lsi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * ksi_195[k]
                   + f_3 * pc_y[k] * lsi_307[k];

        t_395[k] = f_1 * lsh0_230[k]
                   - f_2 * lsh1_230[k]
                   + f_3 * pc_z[k] * lsi_307[k];

        t_396[k] = pa_z[k] * ksk0_216[k]
                   - f_12 * pc_z[k] * ksk1_216[k];

        t_397[k] = f_15 * ksi_196[k]
                   + f_3 * pc_y[k] * lsi_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_z, pc_y, pc_z, ksk0_219, ksi_168, ksi_198, \
                         ksk1_219, lsi_308, lsi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * ksi_168[k]
                   + f_3 * pc_z[k] * lsi_308[k];

        t_399[k] = pa_z[k] * ksk0_219[k]
                   - f_12 * pc_z[k] * ksk1_219[k];

        t_400[k] = f_15 * ksi_198[k]
                   + f_3 * pc_y[k] * lsi_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_z, pc_x, pc_z, ksk0_222, ksi_171, ksi_313, \
                         ksk1_222, lsh0_236, lsh1_236, lsi_311, \
                         lsi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_16 * ksi_313[k]
                   + f_10 * lsh0_236[k]
                   - f_11 * lsh1_236[k]
                   + f_3 * pc_x[k] * lsi_313[k];

        t_402[k] = pa_z[k] * ksk0_222[k]
                   - f_12 * pc_z[k] * ksk1_222[k];

        t_403[k] = f_13 * ksi_171[k]
                   + f_3 * pc_z[k] * lsi_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_z, pc_x, pc_y, pc_z, ksk0_226, ksi_201, \
                         ksi_317, ksk1_226, lsh0_240, lsh1_240, lsi_313, \
                         lsi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * ksi_201[k]
                   + f_3 * pc_y[k] * lsi_313[k];

        t_405[k] = f_16 * ksi_317[k]
                   + f_8 * lsh0_240[k]
                   - f_9 * lsh1_240[k]
                   + f_3 * pc_x[k] * lsi_317[k];

        t_406[k] = pa_z[k] * ksk0_226[k]
                   - f_12 * pc_z[k] * ksk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pa_z, pc_y, pc_z, ksk0_228, ksi_174, ksi_175, \
                         ksi_205, ksk1_228, lsi_314, lsi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * ksi_174[k]
                   + f_3 * pc_z[k] * lsi_314[k];

        t_408[k] = pa_z[k] * ksk0_228[k]
                   + f_14 * ksi_175[k]
                   - f_12 * pc_z[k] * ksk1_228[k];

        t_409[k] = f_15 * ksi_205[k]
                   + f_3 * pc_y[k] * lsi_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pa_z, pc_x, pc_z, ksk0_231, ksi_178, ksi_322, \
                         ksk1_231, lsh0_245, lsh1_245, lsi_318, \
                         lsi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_16 * ksi_322[k]
                   + f_6 * lsh0_245[k]
                   - f_7 * lsh1_245[k]
                   + f_3 * pc_x[k] * lsi_322[k];

        t_411[k] = pa_z[k] * ksk0_231[k]
                   - f_12 * pc_z[k] * ksk1_231[k];

        t_412[k] = f_13 * ksi_178[k]
                   + f_3 * pc_z[k] * lsi_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pa_z, pc_y, pc_z, ksk0_233, ksk0_234, ksi_179, \
                         ksi_180, ksi_210, ksk1_233, ksk1_234, \
                         lsi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pa_z[k] * ksk0_233[k]
                   + f_14 * ksi_179[k]
                   - f_12 * pc_z[k] * ksk1_233[k];

        t_414[k] = pa_z[k] * ksk0_234[k]
                   + f_15 * ksi_180[k]
                   - f_12 * pc_z[k] * ksk1_234[k];

        t_415[k] = f_15 * ksi_210[k]
                   + f_3 * pc_y[k] * lsi_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, ksi_328, ksi_329, ksi_330, ksi_331, \
                         lsh0_251, lsh1_251, lsi_328, lsi_329, lsi_330, \
                         lsi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_16 * ksi_328[k]
                   + f_4 * lsh0_251[k]
                   - f_5 * lsh1_251[k]
                   + f_3 * pc_x[k] * lsi_328[k];

        t_417[k] = f_16 * ksi_329[k]
                   + f_3 * pc_x[k] * lsi_329[k];

        t_418[k] = f_16 * ksi_330[k]
                   + f_3 * pc_x[k] * lsi_330[k];

        t_419[k] = f_16 * ksi_331[k]
                   + f_3 * pc_x[k] * lsi_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, ksi_332, ksi_333, ksi_334, ksi_335, \
                         lsi_332, lsi_333, lsi_334, lsi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_16 * ksi_332[k]
                   + f_3 * pc_x[k] * lsi_332[k];

        t_421[k] = f_16 * ksi_333[k]
                   + f_3 * pc_x[k] * lsi_333[k];

        t_422[k] = f_16 * ksi_334[k]
                   + f_3 * pc_x[k] * lsi_334[k];

        t_423[k] = f_16 * ksi_335[k]
                   + f_3 * pc_x[k] * lsi_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pa_z, pc_y, pc_z, ksk0_244, ksi_189, ksi_219, \
                         ksk1_244, lsh0_248, lsh1_248, lsi_329, \
                         lsi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pa_z[k] * ksk0_244[k]
                   - f_12 * pc_z[k] * ksk1_244[k];

        t_425[k] = f_13 * ksi_189[k]
                   + f_3 * pc_z[k] * lsi_329[k];

        t_426[k] = f_15 * ksi_219[k]
                   + f_10 * lsh0_248[k]
                   - f_11 * lsh1_248[k]
                   + f_3 * pc_y[k] * lsi_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, ksi_220, ksi_221, ksi_222, lsh0_249, \
                         lsh0_250, lsh0_251, lsh1_249, lsh1_250, lsh1_251, lsi_332, lsi_333, \
                         lsi_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * ksi_220[k]
                   + f_8 * lsh0_249[k]
                   - f_9 * lsh1_249[k]
                   + f_3 * pc_y[k] * lsi_332[k];

        t_428[k] = f_15 * ksi_221[k]
                   + f_6 * lsh0_250[k]
                   - f_7 * lsh1_250[k]
                   + f_3 * pc_y[k] * lsi_333[k];

        t_429[k] = f_15 * ksi_222[k]
                   + f_4 * lsh0_251[k]
                   - f_5 * lsh1_251[k]
                   + f_3 * pc_y[k] * lsi_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, ksi_195, ksi_223, ksi_336, \
                         lsh0_251, lsh0_252, lsh1_251, lsh1_252, lsi_335, \
                         lsi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * ksi_223[k]
                   + f_3 * pc_y[k] * lsi_335[k];

        t_431[k] = f_13 * ksi_195[k]
                   + f_1 * lsh0_251[k]
                   - f_2 * lsh1_251[k]
                   + f_3 * pc_z[k] * lsi_335[k];

        t_432[k] = f_16 * ksi_336[k]
                   + f_1 * lsh0_252[k]
                   - f_2 * lsh1_252[k]
                   + f_3 * pc_x[k] * lsi_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, ksi_196, ksi_224, \
                         ksi_226, ksi_339, lsh0_255, lsh1_255, lsi_336, lsi_338, \
                         lsi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * ksi_224[k]
                   + f_3 * pc_y[k] * lsi_336[k];

        t_434[k] = f_14 * ksi_196[k]
                   + f_3 * pc_z[k] * lsi_336[k];

        t_435[k] = f_16 * ksi_339[k]
                   + f_10 * lsh0_255[k]
                   - f_11 * lsh1_255[k]
                   + f_3 * pc_x[k] * lsi_339[k];

        t_436[k] = f_14 * ksi_226[k]
                   + f_3 * pc_y[k] * lsi_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, ksi_199, ksi_341, ksi_342, lsh0_257, \
                         lsh0_258, lsh1_257, lsh1_258, lsi_339, lsi_341, \
                         lsi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_16 * ksi_341[k]
                   + f_10 * lsh0_257[k]
                   - f_11 * lsh1_257[k]
                   + f_3 * pc_x[k] * lsi_341[k];

        t_438[k] = f_16 * ksi_342[k]
                   + f_8 * lsh0_258[k]
                   - f_9 * lsh1_258[k]
                   + f_3 * pc_x[k] * lsi_342[k];

        t_439[k] = f_14 * ksi_199[k]
                   + f_3 * pc_z[k] * lsi_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, ksi_229, ksi_345, ksi_346, lsh0_261, \
                         lsh0_262, lsh1_261, lsh1_262, lsi_341, lsi_345, \
                         lsi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * ksi_229[k]
                   + f_3 * pc_y[k] * lsi_341[k];

        t_441[k] = f_16 * ksi_345[k]
                   + f_8 * lsh0_261[k]
                   - f_9 * lsh1_261[k]
                   + f_3 * pc_x[k] * lsi_345[k];

        t_442[k] = f_16 * ksi_346[k]
                   + f_6 * lsh0_262[k]
                   - f_7 * lsh1_262[k]
                   + f_3 * pc_x[k] * lsi_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, ksi_202, ksi_233, ksi_348, \
                         lsh0_264, lsh1_264, lsi_342, lsi_345, \
                         lsi_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * ksi_202[k]
                   + f_3 * pc_z[k] * lsi_342[k];

        t_444[k] = f_16 * ksi_348[k]
                   + f_6 * lsh0_264[k]
                   - f_7 * lsh1_264[k]
                   + f_3 * pc_x[k] * lsi_348[k];

        t_445[k] = f_14 * ksi_233[k]
                   + f_3 * pc_y[k] * lsi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, ksi_206, ksi_350, ksi_351, lsh0_266, \
                         lsh0_267, lsh1_266, lsh1_267, lsi_346, lsi_350, \
                         lsi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_16 * ksi_350[k]
                   + f_6 * lsh0_266[k]
                   - f_7 * lsh1_266[k]
                   + f_3 * pc_x[k] * lsi_350[k];

        t_447[k] = f_16 * ksi_351[k]
                   + f_4 * lsh0_267[k]
                   - f_5 * lsh1_267[k]
                   + f_3 * pc_x[k] * lsi_351[k];

        t_448[k] = f_14 * ksi_206[k]
                   + f_3 * pc_z[k] * lsi_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, ksi_238, ksi_353, ksi_354, lsh0_269, \
                         lsh0_270, lsh1_269, lsh1_270, lsi_350, lsi_353, \
                         lsi_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_16 * ksi_353[k]
                   + f_4 * lsh0_269[k]
                   - f_5 * lsh1_269[k]
                   + f_3 * pc_x[k] * lsi_353[k];

        t_450[k] = f_16 * ksi_354[k]
                   + f_4 * lsh0_270[k]
                   - f_5 * lsh1_270[k]
                   + f_3 * pc_x[k] * lsi_354[k];

        t_451[k] = f_14 * ksi_238[k]
                   + f_3 * pc_y[k] * lsi_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, ksi_356, ksi_357, ksi_358, ksi_359, \
                         lsh0_272, lsh1_272, lsi_356, lsi_357, lsi_358, \
                         lsi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_16 * ksi_356[k]
                   + f_4 * lsh0_272[k]
                   - f_5 * lsh1_272[k]
                   + f_3 * pc_x[k] * lsi_356[k];

        t_453[k] = f_16 * ksi_357[k]
                   + f_3 * pc_x[k] * lsi_357[k];

        t_454[k] = f_16 * ksi_358[k]
                   + f_3 * pc_x[k] * lsi_358[k];

        t_455[k] = f_16 * ksi_359[k]
                   + f_3 * pc_x[k] * lsi_359[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, ksi_360, ksi_361, ksi_362, ksi_363, \
                         lsi_360, lsi_361, lsi_362, lsi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_16 * ksi_360[k]
                   + f_3 * pc_x[k] * lsi_360[k];

        t_457[k] = f_16 * ksi_361[k]
                   + f_3 * pc_x[k] * lsi_361[k];

        t_458[k] = f_16 * ksi_362[k]
                   + f_3 * pc_x[k] * lsi_362[k];

        t_459[k] = f_16 * ksi_363[k]
                   + f_3 * pc_x[k] * lsi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, ksi_217, ksi_245, ksi_247, lsh0_267, \
                         lsh0_269, lsh1_267, lsh1_269, lsi_357, \
                         lsi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * ksi_245[k]
                   + f_1 * lsh0_267[k]
                   - f_2 * lsh1_267[k]
                   + f_3 * pc_y[k] * lsi_357[k];

        t_461[k] = f_14 * ksi_217[k]
                   + f_3 * pc_z[k] * lsi_357[k];

        t_462[k] = f_14 * ksi_247[k]
                   + f_10 * lsh0_269[k]
                   - f_11 * lsh1_269[k]
                   + f_3 * pc_y[k] * lsi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, ksi_248, ksi_249, ksi_250, lsh0_270, \
                         lsh0_271, lsh0_272, lsh1_270, lsh1_271, lsh1_272, lsi_360, lsi_361, \
                         lsi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * ksi_248[k]
                   + f_8 * lsh0_270[k]
                   - f_9 * lsh1_270[k]
                   + f_3 * pc_y[k] * lsi_360[k];

        t_464[k] = f_14 * ksi_249[k]
                   + f_6 * lsh0_271[k]
                   - f_7 * lsh1_271[k]
                   + f_3 * pc_y[k] * lsi_361[k];

        t_465[k] = f_14 * ksi_250[k]
                   + f_4 * lsh0_272[k]
                   - f_5 * lsh1_272[k]
                   + f_3 * pc_y[k] * lsi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_y, pc_y, pc_z, ksk0_324, ksi_223, \
                         ksi_251, ksi_252, ksk1_324, lsh0_272, lsh1_272, lsi_363, \
                         lsi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * ksi_251[k]
                   + f_3 * pc_y[k] * lsi_363[k];

        t_467[k] = f_14 * ksi_223[k]
                   + f_1 * lsh0_272[k]
                   - f_2 * lsh1_272[k]
                   + f_3 * pc_z[k] * lsi_363[k];

        t_468[k] = pa_y[k] * ksk0_324[k]
                   - f_12 * pc_y[k] * ksk1_324[k];

        t_469[k] = f_13 * ksi_252[k]
                   + f_3 * pc_y[k] * lsi_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_y, pc_y, pc_z, ksk0_327, ksk0_329, \
                         ksi_224, ksi_253, ksi_254, ksk1_327, ksk1_329, lsi_364, \
                         lsi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * ksi_224[k]
                   + f_3 * pc_z[k] * lsi_364[k];

        t_471[k] = pa_y[k] * ksk0_327[k]
                   + f_14 * ksi_253[k]
                   - f_12 * pc_y[k] * ksk1_327[k];

        t_472[k] = f_13 * ksi_254[k]
                   + f_3 * pc_y[k] * lsi_366[k];

        t_473[k] = pa_y[k] * ksk0_329[k]
                   - f_12 * pc_y[k] * ksk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_y, pc_y, pc_z, ksk0_330, ksk0_333, \
                         ksi_227, ksi_255, ksi_257, ksk1_330, ksk1_333, lsi_367, \
                         lsi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pa_y[k] * ksk0_330[k]
                   + f_15 * ksi_255[k]
                   - f_12 * pc_y[k] * ksk1_330[k];

        t_475[k] = f_15 * ksi_227[k]
                   + f_3 * pc_z[k] * lsi_367[k];

        t_476[k] = f_13 * ksi_257[k]
                   + f_3 * pc_y[k] * lsi_369[k];

        t_477[k] = pa_y[k] * ksk0_333[k]
                   - f_12 * pc_y[k] * ksk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pa_y, pc_y, pc_z, ksk0_334, ksk0_336, ksi_230, \
                         ksi_258, ksi_260, ksk1_334, ksk1_336, \
                         lsi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pa_y[k] * ksk0_334[k]
                   + f_16 * ksi_258[k]
                   - f_12 * pc_y[k] * ksk1_334[k];

        t_479[k] = f_15 * ksi_230[k]
                   + f_3 * pc_z[k] * lsi_370[k];

        t_480[k] = pa_y[k] * ksk0_336[k]
                   + f_14 * ksi_260[k]
                   - f_12 * pc_y[k] * ksk1_336[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsh0, const size_t lsh1,
                                                          const size_t lsi, const size_t ncols,
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
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);

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
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_338 = buffer.data(ksk0 + 338);
    const auto *ksk0_339 = buffer.data(ksk0 + 339);
    const auto *ksk0_341 = buffer.data(ksk0 + 341);
    const auto *ksk0_342 = buffer.data(ksk0 + 342);
    const auto *ksk0_344 = buffer.data(ksk0 + 344);
    const auto *ksk0_359 = buffer.data(ksk0 + 359);
    const auto *ksk0_360 = buffer.data(ksk0 + 360);
    const auto *ksk0_363 = buffer.data(ksk0 + 363);
    const auto *ksk0_366 = buffer.data(ksk0 + 366);
    const auto *ksk0_370 = buffer.data(ksk0 + 370);
    const auto *ksk0_372 = buffer.data(ksk0 + 372);
    const auto *ksk0_375 = buffer.data(ksk0 + 375);
    const auto *ksk0_377 = buffer.data(ksk0 + 377);
    const auto *ksk0_378 = buffer.data(ksk0 + 378);

    const auto *ksi_234 = buffer.data(ksi + 234);
    const auto *ksi_245 = buffer.data(ksi + 245);
    const auto *ksi_252 = buffer.data(ksi + 252);
    const auto *ksi_261 = buffer.data(ksi + 261);
    const auto *ksi_262 = buffer.data(ksi + 262);
    const auto *ksi_264 = buffer.data(ksi + 264);
    const auto *ksi_265 = buffer.data(ksi + 265);
    const auto *ksi_266 = buffer.data(ksi + 266);
    const auto *ksi_273 = buffer.data(ksi + 273);
    const auto *ksi_275 = buffer.data(ksi + 275);
    const auto *ksi_276 = buffer.data(ksi + 276);
    const auto *ksi_277 = buffer.data(ksi + 277);
    const auto *ksi_278 = buffer.data(ksi + 278);
    const auto *ksi_279 = buffer.data(ksi + 279);
    const auto *ksi_280 = buffer.data(ksi + 280);
    const auto *ksi_283 = buffer.data(ksi + 283);
    const auto *ksi_285 = buffer.data(ksi + 285);
    const auto *ksi_286 = buffer.data(ksi + 286);
    const auto *ksi_287 = buffer.data(ksi + 287);
    const auto *ksi_289 = buffer.data(ksi + 289);
    const auto *ksi_290 = buffer.data(ksi + 290);
    const auto *ksi_291 = buffer.data(ksi + 291);
    const auto *ksi_292 = buffer.data(ksi + 292);
    const auto *ksi_294 = buffer.data(ksi + 294);
    const auto *ksi_301 = buffer.data(ksi + 301);
    const auto *ksi_307 = buffer.data(ksi + 307);
    const auto *ksi_308 = buffer.data(ksi + 308);
    const auto *ksi_310 = buffer.data(ksi + 310);
    const auto *ksi_313 = buffer.data(ksi + 313);
    const auto *ksi_317 = buffer.data(ksi + 317);
    const auto *ksi_322 = buffer.data(ksi + 322);
    const auto *ksi_385 = buffer.data(ksi + 385);
    const auto *ksi_386 = buffer.data(ksi + 386);
    const auto *ksi_387 = buffer.data(ksi + 387);
    const auto *ksi_388 = buffer.data(ksi + 388);
    const auto *ksi_389 = buffer.data(ksi + 389);
    const auto *ksi_390 = buffer.data(ksi + 390);
    const auto *ksi_391 = buffer.data(ksi + 391);
    const auto *ksi_392 = buffer.data(ksi + 392);
    const auto *ksi_397 = buffer.data(ksi + 397);
    const auto *ksi_401 = buffer.data(ksi + 401);
    const auto *ksi_406 = buffer.data(ksi + 406);
    const auto *ksi_412 = buffer.data(ksi + 412);
    const auto *ksi_413 = buffer.data(ksi + 413);
    const auto *ksi_414 = buffer.data(ksi + 414);
    const auto *ksi_415 = buffer.data(ksi + 415);
    const auto *ksi_416 = buffer.data(ksi + 416);
    const auto *ksi_417 = buffer.data(ksi + 417);
    const auto *ksi_419 = buffer.data(ksi + 419);
    const auto *ksi_420 = buffer.data(ksi + 420);
    const auto *ksi_423 = buffer.data(ksi + 423);
    const auto *ksi_426 = buffer.data(ksi + 426);
    const auto *ksi_430 = buffer.data(ksi + 430);
    const auto *ksi_435 = buffer.data(ksi + 435);
    const auto *ksi_441 = buffer.data(ksi + 441);
    const auto *ksi_443 = buffer.data(ksi + 443);
    const auto *ksi_444 = buffer.data(ksi + 444);
    const auto *ksi_445 = buffer.data(ksi + 445);
    const auto *ksi_446 = buffer.data(ksi + 446);
    const auto *ksi_447 = buffer.data(ksi + 447);
    const auto *ksi_453 = buffer.data(ksi + 453);
    const auto *ksi_457 = buffer.data(ksi + 457);
    const auto *ksi_462 = buffer.data(ksi + 462);
    const auto *ksi_468 = buffer.data(ksi + 468);
    const auto *ksi_469 = buffer.data(ksi + 469);
    const auto *ksi_470 = buffer.data(ksi + 470);
    const auto *ksi_471 = buffer.data(ksi + 471);

    const auto *ksk1_338 = buffer.data(ksk1 + 338);
    const auto *ksk1_339 = buffer.data(ksk1 + 339);
    const auto *ksk1_341 = buffer.data(ksk1 + 341);
    const auto *ksk1_342 = buffer.data(ksk1 + 342);
    const auto *ksk1_344 = buffer.data(ksk1 + 344);
    const auto *ksk1_359 = buffer.data(ksk1 + 359);
    const auto *ksk1_360 = buffer.data(ksk1 + 360);
    const auto *ksk1_363 = buffer.data(ksk1 + 363);
    const auto *ksk1_366 = buffer.data(ksk1 + 366);
    const auto *ksk1_370 = buffer.data(ksk1 + 370);
    const auto *ksk1_372 = buffer.data(ksk1 + 372);
    const auto *ksk1_375 = buffer.data(ksk1 + 375);
    const auto *ksk1_377 = buffer.data(ksk1 + 377);
    const auto *ksk1_378 = buffer.data(ksk1 + 378);

    const auto *lsh0_288 = buffer.data(lsh0 + 288);
    const auto *lsh0_290 = buffer.data(lsh0 + 290);
    const auto *lsh0_291 = buffer.data(lsh0 + 291);
    const auto *lsh0_292 = buffer.data(lsh0 + 292);
    const auto *lsh0_293 = buffer.data(lsh0 + 293);
    const auto *lsh0_294 = buffer.data(lsh0 + 294);
    const auto *lsh0_295 = buffer.data(lsh0 + 295);
    const auto *lsh0_296 = buffer.data(lsh0 + 296);
    const auto *lsh0_297 = buffer.data(lsh0 + 297);
    const auto *lsh0_298 = buffer.data(lsh0 + 298);
    const auto *lsh0_299 = buffer.data(lsh0 + 299);
    const auto *lsh0_300 = buffer.data(lsh0 + 300);
    const auto *lsh0_301 = buffer.data(lsh0 + 301);
    const auto *lsh0_302 = buffer.data(lsh0 + 302);
    const auto *lsh0_303 = buffer.data(lsh0 + 303);
    const auto *lsh0_308 = buffer.data(lsh0 + 308);
    const auto *lsh0_309 = buffer.data(lsh0 + 309);
    const auto *lsh0_310 = buffer.data(lsh0 + 310);
    const auto *lsh0_311 = buffer.data(lsh0 + 311);
    const auto *lsh0_312 = buffer.data(lsh0 + 312);
    const auto *lsh0_313 = buffer.data(lsh0 + 313);
    const auto *lsh0_314 = buffer.data(lsh0 + 314);
    const auto *lsh0_315 = buffer.data(lsh0 + 315);
    const auto *lsh0_317 = buffer.data(lsh0 + 317);
    const auto *lsh0_318 = buffer.data(lsh0 + 318);
    const auto *lsh0_320 = buffer.data(lsh0 + 320);
    const auto *lsh0_321 = buffer.data(lsh0 + 321);
    const auto *lsh0_322 = buffer.data(lsh0 + 322);
    const auto *lsh0_324 = buffer.data(lsh0 + 324);
    const auto *lsh0_325 = buffer.data(lsh0 + 325);
    const auto *lsh0_330 = buffer.data(lsh0 + 330);
    const auto *lsh0_331 = buffer.data(lsh0 + 331);
    const auto *lsh0_332 = buffer.data(lsh0 + 332);
    const auto *lsh0_333 = buffer.data(lsh0 + 333);
    const auto *lsh0_335 = buffer.data(lsh0 + 335);
    const auto *lsh0_341 = buffer.data(lsh0 + 341);
    const auto *lsh0_345 = buffer.data(lsh0 + 345);
    const auto *lsh0_350 = buffer.data(lsh0 + 350);
    const auto *lsh0_356 = buffer.data(lsh0 + 356);

    const auto *lsh1_288 = buffer.data(lsh1 + 288);
    const auto *lsh1_290 = buffer.data(lsh1 + 290);
    const auto *lsh1_291 = buffer.data(lsh1 + 291);
    const auto *lsh1_292 = buffer.data(lsh1 + 292);
    const auto *lsh1_293 = buffer.data(lsh1 + 293);
    const auto *lsh1_294 = buffer.data(lsh1 + 294);
    const auto *lsh1_295 = buffer.data(lsh1 + 295);
    const auto *lsh1_296 = buffer.data(lsh1 + 296);
    const auto *lsh1_297 = buffer.data(lsh1 + 297);
    const auto *lsh1_298 = buffer.data(lsh1 + 298);
    const auto *lsh1_299 = buffer.data(lsh1 + 299);
    const auto *lsh1_300 = buffer.data(lsh1 + 300);
    const auto *lsh1_301 = buffer.data(lsh1 + 301);
    const auto *lsh1_302 = buffer.data(lsh1 + 302);
    const auto *lsh1_303 = buffer.data(lsh1 + 303);
    const auto *lsh1_308 = buffer.data(lsh1 + 308);
    const auto *lsh1_309 = buffer.data(lsh1 + 309);
    const auto *lsh1_310 = buffer.data(lsh1 + 310);
    const auto *lsh1_311 = buffer.data(lsh1 + 311);
    const auto *lsh1_312 = buffer.data(lsh1 + 312);
    const auto *lsh1_313 = buffer.data(lsh1 + 313);
    const auto *lsh1_314 = buffer.data(lsh1 + 314);
    const auto *lsh1_315 = buffer.data(lsh1 + 315);
    const auto *lsh1_317 = buffer.data(lsh1 + 317);
    const auto *lsh1_318 = buffer.data(lsh1 + 318);
    const auto *lsh1_320 = buffer.data(lsh1 + 320);
    const auto *lsh1_321 = buffer.data(lsh1 + 321);
    const auto *lsh1_322 = buffer.data(lsh1 + 322);
    const auto *lsh1_324 = buffer.data(lsh1 + 324);
    const auto *lsh1_325 = buffer.data(lsh1 + 325);
    const auto *lsh1_330 = buffer.data(lsh1 + 330);
    const auto *lsh1_331 = buffer.data(lsh1 + 331);
    const auto *lsh1_332 = buffer.data(lsh1 + 332);
    const auto *lsh1_333 = buffer.data(lsh1 + 333);
    const auto *lsh1_335 = buffer.data(lsh1 + 335);
    const auto *lsh1_341 = buffer.data(lsh1 + 341);
    const auto *lsh1_345 = buffer.data(lsh1 + 345);
    const auto *lsh1_350 = buffer.data(lsh1 + 350);
    const auto *lsh1_356 = buffer.data(lsh1 + 356);

    const auto *lsi_373 = buffer.data(lsi + 373);
    const auto *lsi_374 = buffer.data(lsi + 374);
    const auto *lsi_378 = buffer.data(lsi + 378);
    const auto *lsi_385 = buffer.data(lsi + 385);
    const auto *lsi_386 = buffer.data(lsi + 386);
    const auto *lsi_387 = buffer.data(lsi + 387);
    const auto *lsi_388 = buffer.data(lsi + 388);
    const auto *lsi_389 = buffer.data(lsi + 389);
    const auto *lsi_390 = buffer.data(lsi + 390);
    const auto *lsi_391 = buffer.data(lsi + 391);
    const auto *lsi_392 = buffer.data(lsi + 392);
    const auto *lsi_393 = buffer.data(lsi + 393);
    const auto *lsi_394 = buffer.data(lsi + 394);
    const auto *lsi_395 = buffer.data(lsi + 395);
    const auto *lsi_396 = buffer.data(lsi + 396);
    const auto *lsi_397 = buffer.data(lsi + 397);
    const auto *lsi_398 = buffer.data(lsi + 398);
    const auto *lsi_399 = buffer.data(lsi + 399);
    const auto *lsi_400 = buffer.data(lsi + 400);
    const auto *lsi_401 = buffer.data(lsi + 401);
    const auto *lsi_402 = buffer.data(lsi + 402);
    const auto *lsi_403 = buffer.data(lsi + 403);
    const auto *lsi_404 = buffer.data(lsi + 404);
    const auto *lsi_405 = buffer.data(lsi + 405);
    const auto *lsi_406 = buffer.data(lsi + 406);
    const auto *lsi_412 = buffer.data(lsi + 412);
    const auto *lsi_413 = buffer.data(lsi + 413);
    const auto *lsi_414 = buffer.data(lsi + 414);
    const auto *lsi_415 = buffer.data(lsi + 415);
    const auto *lsi_416 = buffer.data(lsi + 416);
    const auto *lsi_417 = buffer.data(lsi + 417);
    const auto *lsi_418 = buffer.data(lsi + 418);
    const auto *lsi_419 = buffer.data(lsi + 419);
    const auto *lsi_420 = buffer.data(lsi + 420);
    const auto *lsi_421 = buffer.data(lsi + 421);
    const auto *lsi_422 = buffer.data(lsi + 422);
    const auto *lsi_423 = buffer.data(lsi + 423);
    const auto *lsi_425 = buffer.data(lsi + 425);
    const auto *lsi_426 = buffer.data(lsi + 426);
    const auto *lsi_427 = buffer.data(lsi + 427);
    const auto *lsi_429 = buffer.data(lsi + 429);
    const auto *lsi_430 = buffer.data(lsi + 430);
    const auto *lsi_431 = buffer.data(lsi + 431);
    const auto *lsi_432 = buffer.data(lsi + 432);
    const auto *lsi_434 = buffer.data(lsi + 434);
    const auto *lsi_435 = buffer.data(lsi + 435);
    const auto *lsi_441 = buffer.data(lsi + 441);
    const auto *lsi_442 = buffer.data(lsi + 442);
    const auto *lsi_443 = buffer.data(lsi + 443);
    const auto *lsi_444 = buffer.data(lsi + 444);
    const auto *lsi_445 = buffer.data(lsi + 445);
    const auto *lsi_446 = buffer.data(lsi + 446);
    const auto *lsi_447 = buffer.data(lsi + 447);
    const auto *lsi_448 = buffer.data(lsi + 448);
    const auto *lsi_450 = buffer.data(lsi + 450);
    const auto *lsi_451 = buffer.data(lsi + 451);
    const auto *lsi_453 = buffer.data(lsi + 453);
    const auto *lsi_454 = buffer.data(lsi + 454);
    const auto *lsi_457 = buffer.data(lsi + 457);
    const auto *lsi_458 = buffer.data(lsi + 458);
    const auto *lsi_462 = buffer.data(lsi + 462);
    const auto *lsi_468 = buffer.data(lsi + 468);
    const auto *lsi_469 = buffer.data(lsi + 469);
    const auto *lsi_470 = buffer.data(lsi + 470);
    const auto *lsi_471 = buffer.data(lsi + 471);

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pa_y, pc_y, pc_z, ksk0_338, ksk0_339, \
                         ksi_234, ksi_261, ksi_262, ksk1_338, ksk1_339, lsi_373, \
                         lsi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * ksi_261[k]
                   + f_3 * pc_y[k] * lsi_373[k];

        t_482[k] = pa_y[k] * ksk0_338[k]
                   - f_12 * pc_y[k] * ksk1_338[k];

        t_483[k] = pa_y[k] * ksk0_339[k]
                   + f_17 * ksi_262[k]
                   - f_12 * pc_y[k] * ksk1_339[k];

        t_484[k] = f_15 * ksi_234[k]
                   + f_3 * pc_z[k] * lsi_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pa_y, pc_y, ksk0_341, ksk0_342, ksk0_344, \
                         ksi_264, ksi_265, ksi_266, ksk1_341, ksk1_342, ksk1_344, \
                         lsi_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pa_y[k] * ksk0_341[k]
                   + f_15 * ksi_264[k]
                   - f_12 * pc_y[k] * ksk1_341[k];

        t_486[k] = pa_y[k] * ksk0_342[k]
                   + f_14 * ksi_265[k]
                   - f_12 * pc_y[k] * ksk1_342[k];

        t_487[k] = f_13 * ksi_266[k]
                   + f_3 * pc_y[k] * lsi_378[k];

        t_488[k] = pa_y[k] * ksk0_344[k]
                   - f_12 * pc_y[k] * ksk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, ksi_385, ksi_386, ksi_387, \
                         ksi_388, ksi_389, lsi_385, lsi_386, lsi_387, lsi_388, \
                         lsi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_16 * ksi_385[k]
                   + f_3 * pc_x[k] * lsi_385[k];

        t_490[k] = f_16 * ksi_386[k]
                   + f_3 * pc_x[k] * lsi_386[k];

        t_491[k] = f_16 * ksi_387[k]
                   + f_3 * pc_x[k] * lsi_387[k];

        t_492[k] = f_16 * ksi_388[k]
                   + f_3 * pc_x[k] * lsi_388[k];

        t_493[k] = f_16 * ksi_389[k]
                   + f_3 * pc_x[k] * lsi_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, ksi_245, ksi_273, \
                         ksi_390, ksi_391, lsh0_288, lsh1_288, lsi_385, lsi_390, \
                         lsi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_16 * ksi_390[k]
                   + f_3 * pc_x[k] * lsi_390[k];

        t_495[k] = f_16 * ksi_391[k]
                   + f_3 * pc_x[k] * lsi_391[k];

        t_496[k] = f_13 * ksi_273[k]
                   + f_1 * lsh0_288[k]
                   - f_2 * lsh1_288[k]
                   + f_3 * pc_y[k] * lsi_385[k];

        t_497[k] = f_15 * ksi_245[k]
                   + f_3 * pc_z[k] * lsi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, ksi_275, ksi_276, ksi_277, lsh0_290, \
                         lsh0_291, lsh0_292, lsh1_290, lsh1_291, lsh1_292, lsi_387, lsi_388, \
                         lsi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * ksi_275[k]
                   + f_10 * lsh0_290[k]
                   - f_11 * lsh1_290[k]
                   + f_3 * pc_y[k] * lsi_387[k];

        t_499[k] = f_13 * ksi_276[k]
                   + f_8 * lsh0_291[k]
                   - f_9 * lsh1_291[k]
                   + f_3 * pc_y[k] * lsi_388[k];

        t_500[k] = f_13 * ksi_277[k]
                   + f_6 * lsh0_292[k]
                   - f_7 * lsh1_292[k]
                   + f_3 * pc_y[k] * lsi_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pa_y, pc_y, ksk0_359, ksi_278, ksi_279, \
                         ksk1_359, lsh0_293, lsh1_293, lsi_390, \
                         lsi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * ksi_278[k]
                   + f_4 * lsh0_293[k]
                   - f_5 * lsh1_293[k]
                   + f_3 * pc_y[k] * lsi_390[k];

        t_502[k] = f_13 * ksi_279[k]
                   + f_3 * pc_y[k] * lsi_391[k];

        t_503[k] = pa_y[k] * ksk0_359[k]
                   - f_12 * pc_y[k] * ksk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, ksi_252, \
                         ksi_392, lsh0_294, lsh1_294, lsi_392, lsi_393, \
                         lsi_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_16 * ksi_392[k]
                   + f_1 * lsh0_294[k]
                   - f_2 * lsh1_294[k]
                   + f_3 * pc_x[k] * lsi_392[k];

        t_505[k] = f_3 * pc_y[k] * lsi_392[k];

        t_506[k] = f_16 * ksi_252[k]
                   + f_3 * pc_z[k] * lsi_392[k];

        t_507[k] = f_4 * lsh0_294[k]
                   - f_5 * lsh1_294[k]
                   + f_3 * pc_y[k] * lsi_393[k];

        t_508[k] = f_3 * pc_y[k] * lsi_394[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, pc_x, pc_y, ksi_397, lsh0_295, lsh0_296, \
                         lsh0_299, lsh1_295, lsh1_296, lsh1_299, lsi_395, lsi_396, \
                         lsi_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_16 * ksi_397[k]
                   + f_10 * lsh0_299[k]
                   - f_11 * lsh1_299[k]
                   + f_3 * pc_x[k] * lsi_397[k];

        t_510[k] = f_6 * lsh0_295[k]
                   - f_7 * lsh1_295[k]
                   + f_3 * pc_y[k] * lsi_395[k];

        t_511[k] = f_4 * lsh0_296[k]
                   - f_5 * lsh1_296[k]
                   + f_3 * pc_y[k] * lsi_396[k];

        t_512[k] = f_3 * pc_y[k] * lsi_397[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_y, ksi_401, lsh0_297, lsh0_298, \
                         lsh0_303, lsh1_297, lsh1_298, lsh1_303, lsi_398, lsi_399, \
                         lsi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_16 * ksi_401[k]
                   + f_8 * lsh0_303[k]
                   - f_9 * lsh1_303[k]
                   + f_3 * pc_x[k] * lsi_401[k];

        t_514[k] = f_8 * lsh0_297[k]
                   - f_9 * lsh1_297[k]
                   + f_3 * pc_y[k] * lsi_398[k];

        t_515[k] = f_6 * lsh0_298[k]
                   - f_7 * lsh1_298[k]
                   + f_3 * pc_y[k] * lsi_399[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pc_x, pc_y, ksi_406, lsh0_299, lsh0_308, \
                         lsh1_299, lsh1_308, lsi_400, lsi_401, \
                         lsi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_4 * lsh0_299[k]
                   - f_5 * lsh1_299[k]
                   + f_3 * pc_y[k] * lsi_400[k];

        t_517[k] = f_3 * pc_y[k] * lsi_401[k];

        t_518[k] = f_16 * ksi_406[k]
                   + f_6 * lsh0_308[k]
                   - f_7 * lsh1_308[k]
                   + f_3 * pc_x[k] * lsi_406[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_y, lsh0_300, lsh0_301, lsh0_302, lsh1_300, \
                         lsh1_301, lsh1_302, lsi_402, lsi_403, \
                         lsi_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_10 * lsh0_300[k]
                   - f_11 * lsh1_300[k]
                   + f_3 * pc_y[k] * lsi_402[k];

        t_520[k] = f_8 * lsh0_301[k]
                   - f_9 * lsh1_301[k]
                   + f_3 * pc_y[k] * lsi_403[k];

        t_521[k] = f_6 * lsh0_302[k]
                   - f_7 * lsh1_302[k]
                   + f_3 * pc_y[k] * lsi_404[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, ksi_412, ksi_413, lsh0_303, \
                         lsh0_314, lsh1_303, lsh1_314, lsi_405, lsi_406, lsi_412, \
                         lsi_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_4 * lsh0_303[k]
                   - f_5 * lsh1_303[k]
                   + f_3 * pc_y[k] * lsi_405[k];

        t_523[k] = f_3 * pc_y[k] * lsi_406[k];

        t_524[k] = f_16 * ksi_412[k]
                   + f_4 * lsh0_314[k]
                   - f_5 * lsh1_314[k]
                   + f_3 * pc_x[k] * lsi_412[k];

        t_525[k] = f_16 * ksi_413[k]
                   + f_3 * pc_x[k] * lsi_413[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, pc_x, pc_y, ksi_414, ksi_415, \
                         ksi_416, ksi_417, lsi_412, lsi_414, lsi_415, lsi_416, \
                         lsi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_16 * ksi_414[k]
                   + f_3 * pc_x[k] * lsi_414[k];

        t_527[k] = f_16 * ksi_415[k]
                   + f_3 * pc_x[k] * lsi_415[k];

        t_528[k] = f_16 * ksi_416[k]
                   + f_3 * pc_x[k] * lsi_416[k];

        t_529[k] = f_16 * ksi_417[k]
                   + f_3 * pc_x[k] * lsi_417[k];

        t_530[k] = f_3 * pc_y[k] * lsi_412[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pc_x, pc_y, ksi_419, lsh0_309, lsh0_310, \
                         lsh1_309, lsh1_310, lsi_413, lsi_414, \
                         lsi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_16 * ksi_419[k]
                   + f_3 * pc_x[k] * lsi_419[k];

        t_532[k] = f_1 * lsh0_309[k]
                   - f_2 * lsh1_309[k]
                   + f_3 * pc_y[k] * lsi_413[k];

        t_533[k] = f_19 * lsh0_310[k]
                   - f_20 * lsh1_310[k]
                   + f_3 * pc_y[k] * lsi_414[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_y, lsh0_311, lsh0_312, lsh0_313, lsh1_311, \
                         lsh1_312, lsh1_313, lsi_415, lsi_416, \
                         lsi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_10 * lsh0_311[k]
                   - f_11 * lsh1_311[k]
                   + f_3 * pc_y[k] * lsi_415[k];

        t_535[k] = f_8 * lsh0_312[k]
                   - f_9 * lsh1_312[k]
                   + f_3 * pc_y[k] * lsi_416[k];

        t_536[k] = f_6 * lsh0_313[k]
                   - f_7 * lsh1_313[k]
                   + f_3 * pc_y[k] * lsi_417[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, ksi_279, ksi_420, \
                         lsh0_314, lsh0_315, lsh1_314, lsh1_315, lsi_418, lsi_419, \
                         lsi_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * lsh0_314[k]
                   - f_5 * lsh1_314[k]
                   + f_3 * pc_y[k] * lsi_418[k];

        t_538[k] = f_3 * pc_y[k] * lsi_419[k];

        t_539[k] = f_16 * ksi_279[k]
                   + f_1 * lsh0_314[k]
                   - f_2 * lsh1_314[k]
                   + f_3 * pc_z[k] * lsi_419[k];

        t_540[k] = f_15 * ksi_420[k]
                   + f_1 * lsh0_315[k]
                   - f_2 * lsh1_315[k]
                   + f_3 * pc_x[k] * lsi_420[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, ksi_280, ksi_423, \
                         lsh0_318, lsh1_318, lsi_420, lsi_421, \
                         lsi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_17 * ksi_280[k]
                   + f_3 * pc_y[k] * lsi_420[k];

        t_542[k] = f_3 * pc_z[k] * lsi_420[k];

        t_543[k] = f_15 * ksi_423[k]
                   + f_10 * lsh0_318[k]
                   - f_11 * lsh1_318[k]
                   + f_3 * pc_x[k] * lsi_423[k];

        t_544[k] = f_3 * pc_z[k] * lsi_421[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, ksi_426, lsh0_315, lsh0_321, \
                         lsh1_315, lsh1_321, lsi_422, lsi_423, \
                         lsi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * lsh0_315[k]
                   - f_5 * lsh1_315[k]
                   + f_3 * pc_z[k] * lsi_422[k];

        t_546[k] = f_15 * ksi_426[k]
                   + f_8 * lsh0_321[k]
                   - f_9 * lsh1_321[k]
                   + f_3 * pc_x[k] * lsi_426[k];

        t_547[k] = f_3 * pc_z[k] * lsi_423[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pc_x, pc_y, pc_z, ksi_285, ksi_430, \
                         lsh0_317, lsh0_325, lsh1_317, lsh1_325, lsi_425, lsi_426, \
                         lsi_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_17 * ksi_285[k]
                   + f_3 * pc_y[k] * lsi_425[k];

        t_549[k] = f_6 * lsh0_317[k]
                   - f_7 * lsh1_317[k]
                   + f_3 * pc_z[k] * lsi_425[k];

        t_550[k] = f_15 * ksi_430[k]
                   + f_6 * lsh0_325[k]
                   - f_7 * lsh1_325[k]
                   + f_3 * pc_x[k] * lsi_430[k];

        t_551[k] = f_3 * pc_z[k] * lsi_426[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pc_y, pc_z, ksi_289, lsh0_318, lsh0_320, \
                         lsh1_318, lsh1_320, lsi_427, lsi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_4 * lsh0_318[k]
                   - f_5 * lsh1_318[k]
                   + f_3 * pc_z[k] * lsi_427[k];

        t_553[k] = f_17 * ksi_289[k]
                   + f_3 * pc_y[k] * lsi_429[k];

        t_554[k] = f_8 * lsh0_320[k]
                   - f_9 * lsh1_320[k]
                   + f_3 * pc_z[k] * lsi_429[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pc_x, pc_z, ksi_435, lsh0_321, lsh0_330, \
                         lsh1_321, lsh1_330, lsi_430, lsi_431, \
                         lsi_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_15 * ksi_435[k]
                   + f_4 * lsh0_330[k]
                   - f_5 * lsh1_330[k]
                   + f_3 * pc_x[k] * lsi_435[k];

        t_556[k] = f_3 * pc_z[k] * lsi_430[k];

        t_557[k] = f_4 * lsh0_321[k]
                   - f_5 * lsh1_321[k]
                   + f_3 * pc_z[k] * lsi_431[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pc_x, pc_y, pc_z, ksi_294, ksi_441, \
                         lsh0_322, lsh0_324, lsh1_322, lsh1_324, lsi_432, lsi_434, \
                         lsi_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_6 * lsh0_322[k]
                   - f_7 * lsh1_322[k]
                   + f_3 * pc_z[k] * lsi_432[k];

        t_559[k] = f_17 * ksi_294[k]
                   + f_3 * pc_y[k] * lsi_434[k];

        t_560[k] = f_10 * lsh0_324[k]
                   - f_11 * lsh1_324[k]
                   + f_3 * pc_z[k] * lsi_434[k];

        t_561[k] = f_15 * ksi_441[k]
                   + f_3 * pc_x[k] * lsi_441[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, pc_x, pc_z, ksi_443, ksi_444, \
                         ksi_445, ksi_446, lsi_435, lsi_443, lsi_444, lsi_445, \
                         lsi_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_3 * pc_z[k] * lsi_435[k];

        t_563[k] = f_15 * ksi_443[k]
                   + f_3 * pc_x[k] * lsi_443[k];

        t_564[k] = f_15 * ksi_444[k]
                   + f_3 * pc_x[k] * lsi_444[k];

        t_565[k] = f_15 * ksi_445[k]
                   + f_3 * pc_x[k] * lsi_445[k];

        t_566[k] = f_15 * ksi_446[k]
                   + f_3 * pc_x[k] * lsi_446[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, ksi_301, ksi_447, \
                         lsh0_330, lsh1_330, lsi_441, lsi_442, \
                         lsi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_15 * ksi_447[k]
                   + f_3 * pc_x[k] * lsi_447[k];

        t_568[k] = f_17 * ksi_301[k]
                   + f_1 * lsh0_330[k]
                   - f_2 * lsh1_330[k]
                   + f_3 * pc_y[k] * lsi_441[k];

        t_569[k] = f_3 * pc_z[k] * lsi_441[k];

        t_570[k] = f_4 * lsh0_330[k]
                   - f_5 * lsh1_330[k]
                   + f_3 * pc_z[k] * lsi_442[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_z, lsh0_331, lsh0_332, lsh0_333, lsh1_331, \
                         lsh1_332, lsh1_333, lsi_443, lsi_444, \
                         lsi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_6 * lsh0_331[k]
                   - f_7 * lsh1_331[k]
                   + f_3 * pc_z[k] * lsi_443[k];

        t_572[k] = f_8 * lsh0_332[k]
                   - f_9 * lsh1_332[k]
                   + f_3 * pc_z[k] * lsi_444[k];

        t_573[k] = f_10 * lsh0_333[k]
                   - f_11 * lsh1_333[k]
                   + f_3 * pc_z[k] * lsi_445[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_z, pc_y, pc_z, ksk0_360, ksi_307, \
                         ksi_308, ksk1_360, lsh0_335, lsh1_335, lsi_447, \
                         lsi_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * ksi_307[k]
                   + f_3 * pc_y[k] * lsi_447[k];

        t_575[k] = f_1 * lsh0_335[k]
                   - f_2 * lsh1_335[k]
                   + f_3 * pc_z[k] * lsi_447[k];

        t_576[k] = pa_z[k] * ksk0_360[k]
                   - f_12 * pc_z[k] * ksk1_360[k];

        t_577[k] = f_16 * ksi_308[k]
                   + f_3 * pc_y[k] * lsi_448[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pa_z, pc_y, pc_z, ksk0_363, ksi_280, ksi_310, \
                         ksk1_363, lsi_448, lsi_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * ksi_280[k]
                   + f_3 * pc_z[k] * lsi_448[k];

        t_579[k] = pa_z[k] * ksk0_363[k]
                   - f_12 * pc_z[k] * ksk1_363[k];

        t_580[k] = f_16 * ksi_310[k]
                   + f_3 * pc_y[k] * lsi_450[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pa_z, pc_x, pc_z, ksk0_366, ksi_283, ksi_453, \
                         ksk1_366, lsh0_341, lsh1_341, lsi_451, \
                         lsi_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_15 * ksi_453[k]
                   + f_10 * lsh0_341[k]
                   - f_11 * lsh1_341[k]
                   + f_3 * pc_x[k] * lsi_453[k];

        t_582[k] = pa_z[k] * ksk0_366[k]
                   - f_12 * pc_z[k] * ksk1_366[k];

        t_583[k] = f_13 * ksi_283[k]
                   + f_3 * pc_z[k] * lsi_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pa_z, pc_x, pc_y, pc_z, ksk0_370, ksi_313, \
                         ksi_457, ksk1_370, lsh0_345, lsh1_345, lsi_453, \
                         lsi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * ksi_313[k]
                   + f_3 * pc_y[k] * lsi_453[k];

        t_585[k] = f_15 * ksi_457[k]
                   + f_8 * lsh0_345[k]
                   - f_9 * lsh1_345[k]
                   + f_3 * pc_x[k] * lsi_457[k];

        t_586[k] = pa_z[k] * ksk0_370[k]
                   - f_12 * pc_z[k] * ksk1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pa_z, pc_y, pc_z, ksk0_372, ksi_286, ksi_287, \
                         ksi_317, ksk1_372, lsi_454, lsi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_13 * ksi_286[k]
                   + f_3 * pc_z[k] * lsi_454[k];

        t_588[k] = pa_z[k] * ksk0_372[k]
                   + f_14 * ksi_287[k]
                   - f_12 * pc_z[k] * ksk1_372[k];

        t_589[k] = f_16 * ksi_317[k]
                   + f_3 * pc_y[k] * lsi_457[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pa_z, pc_x, pc_z, ksk0_375, ksi_290, ksi_462, \
                         ksk1_375, lsh0_350, lsh1_350, lsi_458, \
                         lsi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_15 * ksi_462[k]
                   + f_6 * lsh0_350[k]
                   - f_7 * lsh1_350[k]
                   + f_3 * pc_x[k] * lsi_462[k];

        t_591[k] = pa_z[k] * ksk0_375[k]
                   - f_12 * pc_z[k] * ksk1_375[k];

        t_592[k] = f_13 * ksi_290[k]
                   + f_3 * pc_z[k] * lsi_458[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pa_z, pc_y, pc_z, ksk0_377, ksk0_378, ksi_291, \
                         ksi_292, ksi_322, ksk1_377, ksk1_378, \
                         lsi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pa_z[k] * ksk0_377[k]
                   + f_14 * ksi_291[k]
                   - f_12 * pc_z[k] * ksk1_377[k];

        t_594[k] = pa_z[k] * ksk0_378[k]
                   + f_15 * ksi_292[k]
                   - f_12 * pc_z[k] * ksk1_378[k];

        t_595[k] = f_16 * ksi_322[k]
                   + f_3 * pc_y[k] * lsi_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, ksi_468, ksi_469, ksi_470, ksi_471, \
                         lsh0_356, lsh1_356, lsi_468, lsi_469, lsi_470, \
                         lsi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_15 * ksi_468[k]
                   + f_4 * lsh0_356[k]
                   - f_5 * lsh1_356[k]
                   + f_3 * pc_x[k] * lsi_468[k];

        t_597[k] = f_15 * ksi_469[k]
                   + f_3 * pc_x[k] * lsi_469[k];

        t_598[k] = f_15 * ksi_470[k]
                   + f_3 * pc_x[k] * lsi_470[k];

        t_599[k] = f_15 * ksi_471[k]
                   + f_3 * pc_x[k] * lsi_471[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsh0, const size_t lsh1,
                                                          const size_t lsi, const size_t ncols,
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
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_388 = buffer.data(ksk0 + 388);
    const auto *ksk0_504 = buffer.data(ksk0 + 504);
    const auto *ksk0_507 = buffer.data(ksk0 + 507);
    const auto *ksk0_509 = buffer.data(ksk0 + 509);
    const auto *ksk0_510 = buffer.data(ksk0 + 510);
    const auto *ksk0_513 = buffer.data(ksk0 + 513);
    const auto *ksk0_514 = buffer.data(ksk0 + 514);
    const auto *ksk0_516 = buffer.data(ksk0 + 516);
    const auto *ksk0_518 = buffer.data(ksk0 + 518);
    const auto *ksk0_519 = buffer.data(ksk0 + 519);
    const auto *ksk0_521 = buffer.data(ksk0 + 521);
    const auto *ksk0_522 = buffer.data(ksk0 + 522);
    const auto *ksk0_524 = buffer.data(ksk0 + 524);

    const auto *ksi_301 = buffer.data(ksi + 301);
    const auto *ksi_307 = buffer.data(ksi + 307);
    const auto *ksi_308 = buffer.data(ksi + 308);
    const auto *ksi_311 = buffer.data(ksi + 311);
    const auto *ksi_314 = buffer.data(ksi + 314);
    const auto *ksi_318 = buffer.data(ksi + 318);
    const auto *ksi_329 = buffer.data(ksi + 329);
    const auto *ksi_331 = buffer.data(ksi + 331);
    const auto *ksi_332 = buffer.data(ksi + 332);
    const auto *ksi_333 = buffer.data(ksi + 333);
    const auto *ksi_334 = buffer.data(ksi + 334);
    const auto *ksi_335 = buffer.data(ksi + 335);
    const auto *ksi_336 = buffer.data(ksi + 336);
    const auto *ksi_338 = buffer.data(ksi + 338);
    const auto *ksi_339 = buffer.data(ksi + 339);
    const auto *ksi_341 = buffer.data(ksi + 341);
    const auto *ksi_342 = buffer.data(ksi + 342);
    const auto *ksi_345 = buffer.data(ksi + 345);
    const auto *ksi_346 = buffer.data(ksi + 346);
    const auto *ksi_350 = buffer.data(ksi + 350);
    const auto *ksi_357 = buffer.data(ksi + 357);
    const auto *ksi_359 = buffer.data(ksi + 359);
    const auto *ksi_360 = buffer.data(ksi + 360);
    const auto *ksi_361 = buffer.data(ksi + 361);
    const auto *ksi_362 = buffer.data(ksi + 362);
    const auto *ksi_363 = buffer.data(ksi + 363);
    const auto *ksi_364 = buffer.data(ksi + 364);
    const auto *ksi_366 = buffer.data(ksi + 366);
    const auto *ksi_367 = buffer.data(ksi + 367);
    const auto *ksi_369 = buffer.data(ksi + 369);
    const auto *ksi_370 = buffer.data(ksi + 370);
    const auto *ksi_373 = buffer.data(ksi + 373);
    const auto *ksi_374 = buffer.data(ksi + 374);
    const auto *ksi_378 = buffer.data(ksi + 378);
    const auto *ksi_385 = buffer.data(ksi + 385);
    const auto *ksi_387 = buffer.data(ksi + 387);
    const auto *ksi_388 = buffer.data(ksi + 388);
    const auto *ksi_389 = buffer.data(ksi + 389);
    const auto *ksi_390 = buffer.data(ksi + 390);
    const auto *ksi_391 = buffer.data(ksi + 391);
    const auto *ksi_392 = buffer.data(ksi + 392);
    const auto *ksi_393 = buffer.data(ksi + 393);
    const auto *ksi_394 = buffer.data(ksi + 394);
    const auto *ksi_395 = buffer.data(ksi + 395);
    const auto *ksi_397 = buffer.data(ksi + 397);
    const auto *ksi_398 = buffer.data(ksi + 398);
    const auto *ksi_400 = buffer.data(ksi + 400);
    const auto *ksi_401 = buffer.data(ksi + 401);
    const auto *ksi_402 = buffer.data(ksi + 402);
    const auto *ksi_404 = buffer.data(ksi + 404);
    const auto *ksi_405 = buffer.data(ksi + 405);
    const auto *ksi_406 = buffer.data(ksi + 406);
    const auto *ksi_472 = buffer.data(ksi + 472);
    const auto *ksi_473 = buffer.data(ksi + 473);
    const auto *ksi_474 = buffer.data(ksi + 474);
    const auto *ksi_475 = buffer.data(ksi + 475);
    const auto *ksi_476 = buffer.data(ksi + 476);
    const auto *ksi_479 = buffer.data(ksi + 479);
    const auto *ksi_481 = buffer.data(ksi + 481);
    const auto *ksi_482 = buffer.data(ksi + 482);
    const auto *ksi_485 = buffer.data(ksi + 485);
    const auto *ksi_486 = buffer.data(ksi + 486);
    const auto *ksi_488 = buffer.data(ksi + 488);
    const auto *ksi_490 = buffer.data(ksi + 490);
    const auto *ksi_491 = buffer.data(ksi + 491);
    const auto *ksi_493 = buffer.data(ksi + 493);
    const auto *ksi_494 = buffer.data(ksi + 494);
    const auto *ksi_496 = buffer.data(ksi + 496);
    const auto *ksi_497 = buffer.data(ksi + 497);
    const auto *ksi_498 = buffer.data(ksi + 498);
    const auto *ksi_499 = buffer.data(ksi + 499);
    const auto *ksi_500 = buffer.data(ksi + 500);
    const auto *ksi_501 = buffer.data(ksi + 501);
    const auto *ksi_502 = buffer.data(ksi + 502);
    const auto *ksi_503 = buffer.data(ksi + 503);
    const auto *ksi_504 = buffer.data(ksi + 504);
    const auto *ksi_507 = buffer.data(ksi + 507);
    const auto *ksi_509 = buffer.data(ksi + 509);
    const auto *ksi_510 = buffer.data(ksi + 510);
    const auto *ksi_513 = buffer.data(ksi + 513);
    const auto *ksi_514 = buffer.data(ksi + 514);
    const auto *ksi_516 = buffer.data(ksi + 516);
    const auto *ksi_518 = buffer.data(ksi + 518);
    const auto *ksi_519 = buffer.data(ksi + 519);
    const auto *ksi_521 = buffer.data(ksi + 521);
    const auto *ksi_522 = buffer.data(ksi + 522);
    const auto *ksi_524 = buffer.data(ksi + 524);
    const auto *ksi_525 = buffer.data(ksi + 525);
    const auto *ksi_526 = buffer.data(ksi + 526);
    const auto *ksi_527 = buffer.data(ksi + 527);
    const auto *ksi_528 = buffer.data(ksi + 528);
    const auto *ksi_529 = buffer.data(ksi + 529);
    const auto *ksi_530 = buffer.data(ksi + 530);
    const auto *ksi_531 = buffer.data(ksi + 531);
    const auto *ksi_553 = buffer.data(ksi + 553);
    const auto *ksi_554 = buffer.data(ksi + 554);
    const auto *ksi_555 = buffer.data(ksi + 555);
    const auto *ksi_556 = buffer.data(ksi + 556);
    const auto *ksi_557 = buffer.data(ksi + 557);

    const auto *ksk1_388 = buffer.data(ksk1 + 388);
    const auto *ksk1_504 = buffer.data(ksk1 + 504);
    const auto *ksk1_507 = buffer.data(ksk1 + 507);
    const auto *ksk1_509 = buffer.data(ksk1 + 509);
    const auto *ksk1_510 = buffer.data(ksk1 + 510);
    const auto *ksk1_513 = buffer.data(ksk1 + 513);
    const auto *ksk1_514 = buffer.data(ksk1 + 514);
    const auto *ksk1_516 = buffer.data(ksk1 + 516);
    const auto *ksk1_518 = buffer.data(ksk1 + 518);
    const auto *ksk1_519 = buffer.data(ksk1 + 519);
    const auto *ksk1_521 = buffer.data(ksk1 + 521);
    const auto *ksk1_522 = buffer.data(ksk1 + 522);
    const auto *ksk1_524 = buffer.data(ksk1 + 524);

    const auto *lsh0_353 = buffer.data(lsh0 + 353);
    const auto *lsh0_354 = buffer.data(lsh0 + 354);
    const auto *lsh0_355 = buffer.data(lsh0 + 355);
    const auto *lsh0_356 = buffer.data(lsh0 + 356);
    const auto *lsh0_357 = buffer.data(lsh0 + 357);
    const auto *lsh0_360 = buffer.data(lsh0 + 360);
    const auto *lsh0_362 = buffer.data(lsh0 + 362);
    const auto *lsh0_363 = buffer.data(lsh0 + 363);
    const auto *lsh0_366 = buffer.data(lsh0 + 366);
    const auto *lsh0_367 = buffer.data(lsh0 + 367);
    const auto *lsh0_369 = buffer.data(lsh0 + 369);
    const auto *lsh0_371 = buffer.data(lsh0 + 371);
    const auto *lsh0_372 = buffer.data(lsh0 + 372);
    const auto *lsh0_374 = buffer.data(lsh0 + 374);
    const auto *lsh0_375 = buffer.data(lsh0 + 375);
    const auto *lsh0_376 = buffer.data(lsh0 + 376);
    const auto *lsh0_377 = buffer.data(lsh0 + 377);
    const auto *lsh0_378 = buffer.data(lsh0 + 378);
    const auto *lsh0_381 = buffer.data(lsh0 + 381);
    const auto *lsh0_383 = buffer.data(lsh0 + 383);
    const auto *lsh0_384 = buffer.data(lsh0 + 384);
    const auto *lsh0_387 = buffer.data(lsh0 + 387);
    const auto *lsh0_388 = buffer.data(lsh0 + 388);
    const auto *lsh0_390 = buffer.data(lsh0 + 390);
    const auto *lsh0_392 = buffer.data(lsh0 + 392);
    const auto *lsh0_393 = buffer.data(lsh0 + 393);
    const auto *lsh0_395 = buffer.data(lsh0 + 395);
    const auto *lsh0_396 = buffer.data(lsh0 + 396);
    const auto *lsh0_397 = buffer.data(lsh0 + 397);
    const auto *lsh0_398 = buffer.data(lsh0 + 398);

    const auto *lsh1_353 = buffer.data(lsh1 + 353);
    const auto *lsh1_354 = buffer.data(lsh1 + 354);
    const auto *lsh1_355 = buffer.data(lsh1 + 355);
    const auto *lsh1_356 = buffer.data(lsh1 + 356);
    const auto *lsh1_357 = buffer.data(lsh1 + 357);
    const auto *lsh1_360 = buffer.data(lsh1 + 360);
    const auto *lsh1_362 = buffer.data(lsh1 + 362);
    const auto *lsh1_363 = buffer.data(lsh1 + 363);
    const auto *lsh1_366 = buffer.data(lsh1 + 366);
    const auto *lsh1_367 = buffer.data(lsh1 + 367);
    const auto *lsh1_369 = buffer.data(lsh1 + 369);
    const auto *lsh1_371 = buffer.data(lsh1 + 371);
    const auto *lsh1_372 = buffer.data(lsh1 + 372);
    const auto *lsh1_374 = buffer.data(lsh1 + 374);
    const auto *lsh1_375 = buffer.data(lsh1 + 375);
    const auto *lsh1_376 = buffer.data(lsh1 + 376);
    const auto *lsh1_377 = buffer.data(lsh1 + 377);
    const auto *lsh1_378 = buffer.data(lsh1 + 378);
    const auto *lsh1_381 = buffer.data(lsh1 + 381);
    const auto *lsh1_383 = buffer.data(lsh1 + 383);
    const auto *lsh1_384 = buffer.data(lsh1 + 384);
    const auto *lsh1_387 = buffer.data(lsh1 + 387);
    const auto *lsh1_388 = buffer.data(lsh1 + 388);
    const auto *lsh1_390 = buffer.data(lsh1 + 390);
    const auto *lsh1_392 = buffer.data(lsh1 + 392);
    const auto *lsh1_393 = buffer.data(lsh1 + 393);
    const auto *lsh1_395 = buffer.data(lsh1 + 395);
    const auto *lsh1_396 = buffer.data(lsh1 + 396);
    const auto *lsh1_397 = buffer.data(lsh1 + 397);
    const auto *lsh1_398 = buffer.data(lsh1 + 398);

    const auto *lsi_469 = buffer.data(lsi + 469);
    const auto *lsi_471 = buffer.data(lsi + 471);
    const auto *lsi_472 = buffer.data(lsi + 472);
    const auto *lsi_473 = buffer.data(lsi + 473);
    const auto *lsi_474 = buffer.data(lsi + 474);
    const auto *lsi_475 = buffer.data(lsi + 475);
    const auto *lsi_476 = buffer.data(lsi + 476);
    const auto *lsi_478 = buffer.data(lsi + 478);
    const auto *lsi_479 = buffer.data(lsi + 479);
    const auto *lsi_481 = buffer.data(lsi + 481);
    const auto *lsi_482 = buffer.data(lsi + 482);
    const auto *lsi_485 = buffer.data(lsi + 485);
    const auto *lsi_486 = buffer.data(lsi + 486);
    const auto *lsi_488 = buffer.data(lsi + 488);
    const auto *lsi_490 = buffer.data(lsi + 490);
    const auto *lsi_491 = buffer.data(lsi + 491);
    const auto *lsi_493 = buffer.data(lsi + 493);
    const auto *lsi_494 = buffer.data(lsi + 494);
    const auto *lsi_496 = buffer.data(lsi + 496);
    const auto *lsi_497 = buffer.data(lsi + 497);
    const auto *lsi_498 = buffer.data(lsi + 498);
    const auto *lsi_499 = buffer.data(lsi + 499);
    const auto *lsi_500 = buffer.data(lsi + 500);
    const auto *lsi_501 = buffer.data(lsi + 501);
    const auto *lsi_502 = buffer.data(lsi + 502);
    const auto *lsi_503 = buffer.data(lsi + 503);
    const auto *lsi_504 = buffer.data(lsi + 504);
    const auto *lsi_506 = buffer.data(lsi + 506);
    const auto *lsi_507 = buffer.data(lsi + 507);
    const auto *lsi_509 = buffer.data(lsi + 509);
    const auto *lsi_510 = buffer.data(lsi + 510);
    const auto *lsi_513 = buffer.data(lsi + 513);
    const auto *lsi_514 = buffer.data(lsi + 514);
    const auto *lsi_516 = buffer.data(lsi + 516);
    const auto *lsi_518 = buffer.data(lsi + 518);
    const auto *lsi_519 = buffer.data(lsi + 519);
    const auto *lsi_521 = buffer.data(lsi + 521);
    const auto *lsi_522 = buffer.data(lsi + 522);
    const auto *lsi_524 = buffer.data(lsi + 524);
    const auto *lsi_525 = buffer.data(lsi + 525);
    const auto *lsi_526 = buffer.data(lsi + 526);
    const auto *lsi_527 = buffer.data(lsi + 527);
    const auto *lsi_528 = buffer.data(lsi + 528);
    const auto *lsi_529 = buffer.data(lsi + 529);
    const auto *lsi_530 = buffer.data(lsi + 530);
    const auto *lsi_531 = buffer.data(lsi + 531);
    const auto *lsi_532 = buffer.data(lsi + 532);
    const auto *lsi_534 = buffer.data(lsi + 534);
    const auto *lsi_535 = buffer.data(lsi + 535);
    const auto *lsi_537 = buffer.data(lsi + 537);
    const auto *lsi_538 = buffer.data(lsi + 538);
    const auto *lsi_541 = buffer.data(lsi + 541);
    const auto *lsi_542 = buffer.data(lsi + 542);
    const auto *lsi_546 = buffer.data(lsi + 546);
    const auto *lsi_553 = buffer.data(lsi + 553);
    const auto *lsi_554 = buffer.data(lsi + 554);
    const auto *lsi_555 = buffer.data(lsi + 555);
    const auto *lsi_556 = buffer.data(lsi + 556);
    const auto *lsi_557 = buffer.data(lsi + 557);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, ksi_472, ksi_473, ksi_474, ksi_475, \
                         lsi_472, lsi_473, lsi_474, lsi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_15 * ksi_472[k]
                   + f_3 * pc_x[k] * lsi_472[k];

        t_601[k] = f_15 * ksi_473[k]
                   + f_3 * pc_x[k] * lsi_473[k];

        t_602[k] = f_15 * ksi_474[k]
                   + f_3 * pc_x[k] * lsi_474[k];

        t_603[k] = f_15 * ksi_475[k]
                   + f_3 * pc_x[k] * lsi_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pa_z, pc_y, pc_z, ksk0_388, ksi_301, ksi_331, \
                         ksk1_388, lsh0_353, lsh1_353, lsi_469, \
                         lsi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pa_z[k] * ksk0_388[k]
                   - f_12 * pc_z[k] * ksk1_388[k];

        t_605[k] = f_13 * ksi_301[k]
                   + f_3 * pc_z[k] * lsi_469[k];

        t_606[k] = f_16 * ksi_331[k]
                   + f_10 * lsh0_353[k]
                   - f_11 * lsh1_353[k]
                   + f_3 * pc_y[k] * lsi_471[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, ksi_332, ksi_333, ksi_334, lsh0_354, \
                         lsh0_355, lsh0_356, lsh1_354, lsh1_355, lsh1_356, lsi_472, lsi_473, \
                         lsi_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_16 * ksi_332[k]
                   + f_8 * lsh0_354[k]
                   - f_9 * lsh1_354[k]
                   + f_3 * pc_y[k] * lsi_472[k];

        t_608[k] = f_16 * ksi_333[k]
                   + f_6 * lsh0_355[k]
                   - f_7 * lsh1_355[k]
                   + f_3 * pc_y[k] * lsi_473[k];

        t_609[k] = f_16 * ksi_334[k]
                   + f_4 * lsh0_356[k]
                   - f_5 * lsh1_356[k]
                   + f_3 * pc_y[k] * lsi_474[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_x, pc_y, pc_z, ksi_307, ksi_335, ksi_476, \
                         lsh0_356, lsh0_357, lsh1_356, lsh1_357, lsi_475, \
                         lsi_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * ksi_335[k]
                   + f_3 * pc_y[k] * lsi_475[k];

        t_611[k] = f_13 * ksi_307[k]
                   + f_1 * lsh0_356[k]
                   - f_2 * lsh1_356[k]
                   + f_3 * pc_z[k] * lsi_475[k];

        t_612[k] = f_15 * ksi_476[k]
                   + f_1 * lsh0_357[k]
                   - f_2 * lsh1_357[k]
                   + f_3 * pc_x[k] * lsi_476[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pc_x, pc_y, pc_z, ksi_308, ksi_336, \
                         ksi_338, ksi_479, lsh0_360, lsh1_360, lsi_476, lsi_478, \
                         lsi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_15 * ksi_336[k]
                   + f_3 * pc_y[k] * lsi_476[k];

        t_614[k] = f_14 * ksi_308[k]
                   + f_3 * pc_z[k] * lsi_476[k];

        t_615[k] = f_15 * ksi_479[k]
                   + f_10 * lsh0_360[k]
                   - f_11 * lsh1_360[k]
                   + f_3 * pc_x[k] * lsi_479[k];

        t_616[k] = f_15 * ksi_338[k]
                   + f_3 * pc_y[k] * lsi_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pc_x, pc_z, ksi_311, ksi_481, ksi_482, lsh0_362, \
                         lsh0_363, lsh1_362, lsh1_363, lsi_479, lsi_481, \
                         lsi_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_15 * ksi_481[k]
                   + f_10 * lsh0_362[k]
                   - f_11 * lsh1_362[k]
                   + f_3 * pc_x[k] * lsi_481[k];

        t_618[k] = f_15 * ksi_482[k]
                   + f_8 * lsh0_363[k]
                   - f_9 * lsh1_363[k]
                   + f_3 * pc_x[k] * lsi_482[k];

        t_619[k] = f_14 * ksi_311[k]
                   + f_3 * pc_z[k] * lsi_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pc_x, pc_y, ksi_341, ksi_485, ksi_486, lsh0_366, \
                         lsh0_367, lsh1_366, lsh1_367, lsi_481, lsi_485, \
                         lsi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * ksi_341[k]
                   + f_3 * pc_y[k] * lsi_481[k];

        t_621[k] = f_15 * ksi_485[k]
                   + f_8 * lsh0_366[k]
                   - f_9 * lsh1_366[k]
                   + f_3 * pc_x[k] * lsi_485[k];

        t_622[k] = f_15 * ksi_486[k]
                   + f_6 * lsh0_367[k]
                   - f_7 * lsh1_367[k]
                   + f_3 * pc_x[k] * lsi_486[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, ksi_314, ksi_345, ksi_488, \
                         lsh0_369, lsh1_369, lsi_482, lsi_485, \
                         lsi_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * ksi_314[k]
                   + f_3 * pc_z[k] * lsi_482[k];

        t_624[k] = f_15 * ksi_488[k]
                   + f_6 * lsh0_369[k]
                   - f_7 * lsh1_369[k]
                   + f_3 * pc_x[k] * lsi_488[k];

        t_625[k] = f_15 * ksi_345[k]
                   + f_3 * pc_y[k] * lsi_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_x, pc_z, ksi_318, ksi_490, ksi_491, lsh0_371, \
                         lsh0_372, lsh1_371, lsh1_372, lsi_486, lsi_490, \
                         lsi_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * ksi_490[k]
                   + f_6 * lsh0_371[k]
                   - f_7 * lsh1_371[k]
                   + f_3 * pc_x[k] * lsi_490[k];

        t_627[k] = f_15 * ksi_491[k]
                   + f_4 * lsh0_372[k]
                   - f_5 * lsh1_372[k]
                   + f_3 * pc_x[k] * lsi_491[k];

        t_628[k] = f_14 * ksi_318[k]
                   + f_3 * pc_z[k] * lsi_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, ksi_350, ksi_493, ksi_494, lsh0_374, \
                         lsh0_375, lsh1_374, lsh1_375, lsi_490, lsi_493, \
                         lsi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_15 * ksi_493[k]
                   + f_4 * lsh0_374[k]
                   - f_5 * lsh1_374[k]
                   + f_3 * pc_x[k] * lsi_493[k];

        t_630[k] = f_15 * ksi_494[k]
                   + f_4 * lsh0_375[k]
                   - f_5 * lsh1_375[k]
                   + f_3 * pc_x[k] * lsi_494[k];

        t_631[k] = f_15 * ksi_350[k]
                   + f_3 * pc_y[k] * lsi_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, ksi_496, ksi_497, ksi_498, ksi_499, \
                         lsh0_377, lsh1_377, lsi_496, lsi_497, lsi_498, \
                         lsi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_15 * ksi_496[k]
                   + f_4 * lsh0_377[k]
                   - f_5 * lsh1_377[k]
                   + f_3 * pc_x[k] * lsi_496[k];

        t_633[k] = f_15 * ksi_497[k]
                   + f_3 * pc_x[k] * lsi_497[k];

        t_634[k] = f_15 * ksi_498[k]
                   + f_3 * pc_x[k] * lsi_498[k];

        t_635[k] = f_15 * ksi_499[k]
                   + f_3 * pc_x[k] * lsi_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, ksi_500, ksi_501, ksi_502, ksi_503, \
                         lsi_500, lsi_501, lsi_502, lsi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_15 * ksi_500[k]
                   + f_3 * pc_x[k] * lsi_500[k];

        t_637[k] = f_15 * ksi_501[k]
                   + f_3 * pc_x[k] * lsi_501[k];

        t_638[k] = f_15 * ksi_502[k]
                   + f_3 * pc_x[k] * lsi_502[k];

        t_639[k] = f_15 * ksi_503[k]
                   + f_3 * pc_x[k] * lsi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, ksi_329, ksi_357, ksi_359, lsh0_372, \
                         lsh0_374, lsh1_372, lsh1_374, lsi_497, \
                         lsi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * ksi_357[k]
                   + f_1 * lsh0_372[k]
                   - f_2 * lsh1_372[k]
                   + f_3 * pc_y[k] * lsi_497[k];

        t_641[k] = f_14 * ksi_329[k]
                   + f_3 * pc_z[k] * lsi_497[k];

        t_642[k] = f_15 * ksi_359[k]
                   + f_10 * lsh0_374[k]
                   - f_11 * lsh1_374[k]
                   + f_3 * pc_y[k] * lsi_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, ksi_360, ksi_361, ksi_362, lsh0_375, \
                         lsh0_376, lsh0_377, lsh1_375, lsh1_376, lsh1_377, lsi_500, lsi_501, \
                         lsi_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * ksi_360[k]
                   + f_8 * lsh0_375[k]
                   - f_9 * lsh1_375[k]
                   + f_3 * pc_y[k] * lsi_500[k];

        t_644[k] = f_15 * ksi_361[k]
                   + f_6 * lsh0_376[k]
                   - f_7 * lsh1_376[k]
                   + f_3 * pc_y[k] * lsi_501[k];

        t_645[k] = f_15 * ksi_362[k]
                   + f_4 * lsh0_377[k]
                   - f_5 * lsh1_377[k]
                   + f_3 * pc_y[k] * lsi_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, ksi_335, ksi_363, ksi_504, \
                         lsh0_377, lsh0_378, lsh1_377, lsh1_378, lsi_503, \
                         lsi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * ksi_363[k]
                   + f_3 * pc_y[k] * lsi_503[k];

        t_647[k] = f_14 * ksi_335[k]
                   + f_1 * lsh0_377[k]
                   - f_2 * lsh1_377[k]
                   + f_3 * pc_z[k] * lsi_503[k];

        t_648[k] = f_15 * ksi_504[k]
                   + f_1 * lsh0_378[k]
                   - f_2 * lsh1_378[k]
                   + f_3 * pc_x[k] * lsi_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, ksi_336, ksi_364, \
                         ksi_366, ksi_507, lsh0_381, lsh1_381, lsi_504, lsi_506, \
                         lsi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_14 * ksi_364[k]
                   + f_3 * pc_y[k] * lsi_504[k];

        t_650[k] = f_15 * ksi_336[k]
                   + f_3 * pc_z[k] * lsi_504[k];

        t_651[k] = f_15 * ksi_507[k]
                   + f_10 * lsh0_381[k]
                   - f_11 * lsh1_381[k]
                   + f_3 * pc_x[k] * lsi_507[k];

        t_652[k] = f_14 * ksi_366[k]
                   + f_3 * pc_y[k] * lsi_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, ksi_339, ksi_509, ksi_510, lsh0_383, \
                         lsh0_384, lsh1_383, lsh1_384, lsi_507, lsi_509, \
                         lsi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_15 * ksi_509[k]
                   + f_10 * lsh0_383[k]
                   - f_11 * lsh1_383[k]
                   + f_3 * pc_x[k] * lsi_509[k];

        t_654[k] = f_15 * ksi_510[k]
                   + f_8 * lsh0_384[k]
                   - f_9 * lsh1_384[k]
                   + f_3 * pc_x[k] * lsi_510[k];

        t_655[k] = f_15 * ksi_339[k]
                   + f_3 * pc_z[k] * lsi_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, ksi_369, ksi_513, ksi_514, lsh0_387, \
                         lsh0_388, lsh1_387, lsh1_388, lsi_509, lsi_513, \
                         lsi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * ksi_369[k]
                   + f_3 * pc_y[k] * lsi_509[k];

        t_657[k] = f_15 * ksi_513[k]
                   + f_8 * lsh0_387[k]
                   - f_9 * lsh1_387[k]
                   + f_3 * pc_x[k] * lsi_513[k];

        t_658[k] = f_15 * ksi_514[k]
                   + f_6 * lsh0_388[k]
                   - f_7 * lsh1_388[k]
                   + f_3 * pc_x[k] * lsi_514[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pc_x, pc_y, pc_z, ksi_342, ksi_373, ksi_516, \
                         lsh0_390, lsh1_390, lsi_510, lsi_513, \
                         lsi_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * ksi_342[k]
                   + f_3 * pc_z[k] * lsi_510[k];

        t_660[k] = f_15 * ksi_516[k]
                   + f_6 * lsh0_390[k]
                   - f_7 * lsh1_390[k]
                   + f_3 * pc_x[k] * lsi_516[k];

        t_661[k] = f_14 * ksi_373[k]
                   + f_3 * pc_y[k] * lsi_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_z, ksi_346, ksi_518, ksi_519, lsh0_392, \
                         lsh0_393, lsh1_392, lsh1_393, lsi_514, lsi_518, \
                         lsi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_15 * ksi_518[k]
                   + f_6 * lsh0_392[k]
                   - f_7 * lsh1_392[k]
                   + f_3 * pc_x[k] * lsi_518[k];

        t_663[k] = f_15 * ksi_519[k]
                   + f_4 * lsh0_393[k]
                   - f_5 * lsh1_393[k]
                   + f_3 * pc_x[k] * lsi_519[k];

        t_664[k] = f_15 * ksi_346[k]
                   + f_3 * pc_z[k] * lsi_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_x, pc_y, ksi_378, ksi_521, ksi_522, lsh0_395, \
                         lsh0_396, lsh1_395, lsh1_396, lsi_518, lsi_521, \
                         lsi_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_15 * ksi_521[k]
                   + f_4 * lsh0_395[k]
                   - f_5 * lsh1_395[k]
                   + f_3 * pc_x[k] * lsi_521[k];

        t_666[k] = f_15 * ksi_522[k]
                   + f_4 * lsh0_396[k]
                   - f_5 * lsh1_396[k]
                   + f_3 * pc_x[k] * lsi_522[k];

        t_667[k] = f_14 * ksi_378[k]
                   + f_3 * pc_y[k] * lsi_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pc_x, ksi_524, ksi_525, ksi_526, ksi_527, \
                         lsh0_398, lsh1_398, lsi_524, lsi_525, lsi_526, \
                         lsi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_15 * ksi_524[k]
                   + f_4 * lsh0_398[k]
                   - f_5 * lsh1_398[k]
                   + f_3 * pc_x[k] * lsi_524[k];

        t_669[k] = f_15 * ksi_525[k]
                   + f_3 * pc_x[k] * lsi_525[k];

        t_670[k] = f_15 * ksi_526[k]
                   + f_3 * pc_x[k] * lsi_526[k];

        t_671[k] = f_15 * ksi_527[k]
                   + f_3 * pc_x[k] * lsi_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, ksi_528, ksi_529, ksi_530, ksi_531, \
                         lsi_528, lsi_529, lsi_530, lsi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_15 * ksi_528[k]
                   + f_3 * pc_x[k] * lsi_528[k];

        t_673[k] = f_15 * ksi_529[k]
                   + f_3 * pc_x[k] * lsi_529[k];

        t_674[k] = f_15 * ksi_530[k]
                   + f_3 * pc_x[k] * lsi_530[k];

        t_675[k] = f_15 * ksi_531[k]
                   + f_3 * pc_x[k] * lsi_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, ksi_357, ksi_385, ksi_387, lsh0_393, \
                         lsh0_395, lsh1_393, lsh1_395, lsi_525, \
                         lsi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * ksi_385[k]
                   + f_1 * lsh0_393[k]
                   - f_2 * lsh1_393[k]
                   + f_3 * pc_y[k] * lsi_525[k];

        t_677[k] = f_15 * ksi_357[k]
                   + f_3 * pc_z[k] * lsi_525[k];

        t_678[k] = f_14 * ksi_387[k]
                   + f_10 * lsh0_395[k]
                   - f_11 * lsh1_395[k]
                   + f_3 * pc_y[k] * lsi_527[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, ksi_388, ksi_389, ksi_390, lsh0_396, \
                         lsh0_397, lsh0_398, lsh1_396, lsh1_397, lsh1_398, lsi_528, lsi_529, \
                         lsi_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * ksi_388[k]
                   + f_8 * lsh0_396[k]
                   - f_9 * lsh1_396[k]
                   + f_3 * pc_y[k] * lsi_528[k];

        t_680[k] = f_14 * ksi_389[k]
                   + f_6 * lsh0_397[k]
                   - f_7 * lsh1_397[k]
                   + f_3 * pc_y[k] * lsi_529[k];

        t_681[k] = f_14 * ksi_390[k]
                   + f_4 * lsh0_398[k]
                   - f_5 * lsh1_398[k]
                   + f_3 * pc_y[k] * lsi_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pa_y, pc_y, pc_z, ksk0_504, ksi_363, \
                         ksi_391, ksi_392, ksk1_504, lsh0_398, lsh1_398, lsi_531, \
                         lsi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * ksi_391[k]
                   + f_3 * pc_y[k] * lsi_531[k];

        t_683[k] = f_15 * ksi_363[k]
                   + f_1 * lsh0_398[k]
                   - f_2 * lsh1_398[k]
                   + f_3 * pc_z[k] * lsi_531[k];

        t_684[k] = pa_y[k] * ksk0_504[k]
                   - f_12 * pc_y[k] * ksk1_504[k];

        t_685[k] = f_13 * ksi_392[k]
                   + f_3 * pc_y[k] * lsi_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pa_y, pc_y, pc_z, ksk0_507, ksk0_509, \
                         ksi_364, ksi_393, ksi_394, ksk1_507, ksk1_509, lsi_532, \
                         lsi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * ksi_364[k]
                   + f_3 * pc_z[k] * lsi_532[k];

        t_687[k] = pa_y[k] * ksk0_507[k]
                   + f_14 * ksi_393[k]
                   - f_12 * pc_y[k] * ksk1_507[k];

        t_688[k] = f_13 * ksi_394[k]
                   + f_3 * pc_y[k] * lsi_534[k];

        t_689[k] = pa_y[k] * ksk0_509[k]
                   - f_12 * pc_y[k] * ksk1_509[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pa_y, pc_y, pc_z, ksk0_510, ksk0_513, \
                         ksi_367, ksi_395, ksi_397, ksk1_510, ksk1_513, lsi_535, \
                         lsi_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pa_y[k] * ksk0_510[k]
                   + f_15 * ksi_395[k]
                   - f_12 * pc_y[k] * ksk1_510[k];

        t_691[k] = f_16 * ksi_367[k]
                   + f_3 * pc_z[k] * lsi_535[k];

        t_692[k] = f_13 * ksi_397[k]
                   + f_3 * pc_y[k] * lsi_537[k];

        t_693[k] = pa_y[k] * ksk0_513[k]
                   - f_12 * pc_y[k] * ksk1_513[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pa_y, pc_y, pc_z, ksk0_514, ksk0_516, ksi_370, \
                         ksi_398, ksi_400, ksk1_514, ksk1_516, \
                         lsi_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pa_y[k] * ksk0_514[k]
                   + f_16 * ksi_398[k]
                   - f_12 * pc_y[k] * ksk1_514[k];

        t_695[k] = f_16 * ksi_370[k]
                   + f_3 * pc_z[k] * lsi_538[k];

        t_696[k] = pa_y[k] * ksk0_516[k]
                   + f_14 * ksi_400[k]
                   - f_12 * pc_y[k] * ksk1_516[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, pa_y, pc_y, pc_z, ksk0_518, ksk0_519, \
                         ksi_374, ksi_401, ksi_402, ksk1_518, ksk1_519, lsi_541, \
                         lsi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * ksi_401[k]
                   + f_3 * pc_y[k] * lsi_541[k];

        t_698[k] = pa_y[k] * ksk0_518[k]
                   - f_12 * pc_y[k] * ksk1_518[k];

        t_699[k] = pa_y[k] * ksk0_519[k]
                   + f_17 * ksi_402[k]
                   - f_12 * pc_y[k] * ksk1_519[k];

        t_700[k] = f_16 * ksi_374[k]
                   + f_3 * pc_z[k] * lsi_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pa_y, pc_y, ksk0_521, ksk0_522, ksk0_524, \
                         ksi_404, ksi_405, ksi_406, ksk1_521, ksk1_522, ksk1_524, \
                         lsi_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pa_y[k] * ksk0_521[k]
                   + f_15 * ksi_404[k]
                   - f_12 * pc_y[k] * ksk1_521[k];

        t_702[k] = pa_y[k] * ksk0_522[k]
                   + f_14 * ksi_405[k]
                   - f_12 * pc_y[k] * ksk1_522[k];

        t_703[k] = f_13 * ksi_406[k]
                   + f_3 * pc_y[k] * lsi_546[k];

        t_704[k] = pa_y[k] * ksk0_524[k]
                   - f_12 * pc_y[k] * ksk1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, ksi_553, ksi_554, ksi_555, \
                         ksi_556, ksi_557, lsi_553, lsi_554, lsi_555, lsi_556, \
                         lsi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_15 * ksi_553[k]
                   + f_3 * pc_x[k] * lsi_553[k];

        t_706[k] = f_15 * ksi_554[k]
                   + f_3 * pc_x[k] * lsi_554[k];

        t_707[k] = f_15 * ksi_555[k]
                   + f_3 * pc_x[k] * lsi_555[k];

        t_708[k] = f_15 * ksi_556[k]
                   + f_3 * pc_x[k] * lsi_556[k];

        t_709[k] = f_15 * ksi_557[k]
                   + f_3 * pc_x[k] * lsi_557[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsh0, const size_t lsh1,
                                                          const size_t lsi, const size_t ncols,
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
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.0 / q;

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
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_539 = buffer.data(ksk0 + 539);
    const auto *ksk0_540 = buffer.data(ksk0 + 540);
    const auto *ksk0_543 = buffer.data(ksk0 + 543);
    const auto *ksk0_546 = buffer.data(ksk0 + 546);
    const auto *ksk0_550 = buffer.data(ksk0 + 550);
    const auto *ksk0_552 = buffer.data(ksk0 + 552);
    const auto *ksk0_555 = buffer.data(ksk0 + 555);
    const auto *ksk0_557 = buffer.data(ksk0 + 557);
    const auto *ksk0_558 = buffer.data(ksk0 + 558);
    const auto *ksk0_568 = buffer.data(ksk0 + 568);

    const auto *ksi_385 = buffer.data(ksi + 385);
    const auto *ksi_392 = buffer.data(ksi + 392);
    const auto *ksi_413 = buffer.data(ksi + 413);
    const auto *ksi_415 = buffer.data(ksi + 415);
    const auto *ksi_416 = buffer.data(ksi + 416);
    const auto *ksi_417 = buffer.data(ksi + 417);
    const auto *ksi_418 = buffer.data(ksi + 418);
    const auto *ksi_419 = buffer.data(ksi + 419);
    const auto *ksi_420 = buffer.data(ksi + 420);
    const auto *ksi_423 = buffer.data(ksi + 423);
    const auto *ksi_425 = buffer.data(ksi + 425);
    const auto *ksi_426 = buffer.data(ksi + 426);
    const auto *ksi_427 = buffer.data(ksi + 427);
    const auto *ksi_429 = buffer.data(ksi + 429);
    const auto *ksi_430 = buffer.data(ksi + 430);
    const auto *ksi_431 = buffer.data(ksi + 431);
    const auto *ksi_432 = buffer.data(ksi + 432);
    const auto *ksi_434 = buffer.data(ksi + 434);
    const auto *ksi_441 = buffer.data(ksi + 441);
    const auto *ksi_447 = buffer.data(ksi + 447);
    const auto *ksi_448 = buffer.data(ksi + 448);
    const auto *ksi_450 = buffer.data(ksi + 450);
    const auto *ksi_453 = buffer.data(ksi + 453);
    const auto *ksi_457 = buffer.data(ksi + 457);
    const auto *ksi_462 = buffer.data(ksi + 462);
    const auto *ksi_471 = buffer.data(ksi + 471);
    const auto *ksi_472 = buffer.data(ksi + 472);
    const auto *ksi_473 = buffer.data(ksi + 473);
    const auto *ksi_474 = buffer.data(ksi + 474);
    const auto *ksi_558 = buffer.data(ksi + 558);
    const auto *ksi_559 = buffer.data(ksi + 559);
    const auto *ksi_560 = buffer.data(ksi + 560);
    const auto *ksi_565 = buffer.data(ksi + 565);
    const auto *ksi_569 = buffer.data(ksi + 569);
    const auto *ksi_574 = buffer.data(ksi + 574);
    const auto *ksi_580 = buffer.data(ksi + 580);
    const auto *ksi_581 = buffer.data(ksi + 581);
    const auto *ksi_582 = buffer.data(ksi + 582);
    const auto *ksi_583 = buffer.data(ksi + 583);
    const auto *ksi_584 = buffer.data(ksi + 584);
    const auto *ksi_585 = buffer.data(ksi + 585);
    const auto *ksi_587 = buffer.data(ksi + 587);
    const auto *ksi_588 = buffer.data(ksi + 588);
    const auto *ksi_591 = buffer.data(ksi + 591);
    const auto *ksi_594 = buffer.data(ksi + 594);
    const auto *ksi_598 = buffer.data(ksi + 598);
    const auto *ksi_603 = buffer.data(ksi + 603);
    const auto *ksi_609 = buffer.data(ksi + 609);
    const auto *ksi_611 = buffer.data(ksi + 611);
    const auto *ksi_612 = buffer.data(ksi + 612);
    const auto *ksi_613 = buffer.data(ksi + 613);
    const auto *ksi_614 = buffer.data(ksi + 614);
    const auto *ksi_615 = buffer.data(ksi + 615);
    const auto *ksi_621 = buffer.data(ksi + 621);
    const auto *ksi_625 = buffer.data(ksi + 625);
    const auto *ksi_630 = buffer.data(ksi + 630);
    const auto *ksi_636 = buffer.data(ksi + 636);
    const auto *ksi_637 = buffer.data(ksi + 637);
    const auto *ksi_638 = buffer.data(ksi + 638);
    const auto *ksi_639 = buffer.data(ksi + 639);
    const auto *ksi_640 = buffer.data(ksi + 640);
    const auto *ksi_641 = buffer.data(ksi + 641);
    const auto *ksi_642 = buffer.data(ksi + 642);
    const auto *ksi_643 = buffer.data(ksi + 643);

    const auto *ksk1_539 = buffer.data(ksk1 + 539);
    const auto *ksk1_540 = buffer.data(ksk1 + 540);
    const auto *ksk1_543 = buffer.data(ksk1 + 543);
    const auto *ksk1_546 = buffer.data(ksk1 + 546);
    const auto *ksk1_550 = buffer.data(ksk1 + 550);
    const auto *ksk1_552 = buffer.data(ksk1 + 552);
    const auto *ksk1_555 = buffer.data(ksk1 + 555);
    const auto *ksk1_557 = buffer.data(ksk1 + 557);
    const auto *ksk1_558 = buffer.data(ksk1 + 558);
    const auto *ksk1_568 = buffer.data(ksk1 + 568);

    const auto *lsh0_414 = buffer.data(lsh0 + 414);
    const auto *lsh0_416 = buffer.data(lsh0 + 416);
    const auto *lsh0_417 = buffer.data(lsh0 + 417);
    const auto *lsh0_418 = buffer.data(lsh0 + 418);
    const auto *lsh0_419 = buffer.data(lsh0 + 419);
    const auto *lsh0_420 = buffer.data(lsh0 + 420);
    const auto *lsh0_421 = buffer.data(lsh0 + 421);
    const auto *lsh0_422 = buffer.data(lsh0 + 422);
    const auto *lsh0_423 = buffer.data(lsh0 + 423);
    const auto *lsh0_424 = buffer.data(lsh0 + 424);
    const auto *lsh0_425 = buffer.data(lsh0 + 425);
    const auto *lsh0_426 = buffer.data(lsh0 + 426);
    const auto *lsh0_427 = buffer.data(lsh0 + 427);
    const auto *lsh0_428 = buffer.data(lsh0 + 428);
    const auto *lsh0_429 = buffer.data(lsh0 + 429);
    const auto *lsh0_434 = buffer.data(lsh0 + 434);
    const auto *lsh0_435 = buffer.data(lsh0 + 435);
    const auto *lsh0_436 = buffer.data(lsh0 + 436);
    const auto *lsh0_437 = buffer.data(lsh0 + 437);
    const auto *lsh0_438 = buffer.data(lsh0 + 438);
    const auto *lsh0_439 = buffer.data(lsh0 + 439);
    const auto *lsh0_440 = buffer.data(lsh0 + 440);
    const auto *lsh0_441 = buffer.data(lsh0 + 441);
    const auto *lsh0_443 = buffer.data(lsh0 + 443);
    const auto *lsh0_444 = buffer.data(lsh0 + 444);
    const auto *lsh0_446 = buffer.data(lsh0 + 446);
    const auto *lsh0_447 = buffer.data(lsh0 + 447);
    const auto *lsh0_448 = buffer.data(lsh0 + 448);
    const auto *lsh0_450 = buffer.data(lsh0 + 450);
    const auto *lsh0_451 = buffer.data(lsh0 + 451);
    const auto *lsh0_456 = buffer.data(lsh0 + 456);
    const auto *lsh0_457 = buffer.data(lsh0 + 457);
    const auto *lsh0_458 = buffer.data(lsh0 + 458);
    const auto *lsh0_459 = buffer.data(lsh0 + 459);
    const auto *lsh0_461 = buffer.data(lsh0 + 461);
    const auto *lsh0_467 = buffer.data(lsh0 + 467);
    const auto *lsh0_471 = buffer.data(lsh0 + 471);
    const auto *lsh0_476 = buffer.data(lsh0 + 476);
    const auto *lsh0_479 = buffer.data(lsh0 + 479);
    const auto *lsh0_480 = buffer.data(lsh0 + 480);
    const auto *lsh0_481 = buffer.data(lsh0 + 481);
    const auto *lsh0_482 = buffer.data(lsh0 + 482);

    const auto *lsh1_414 = buffer.data(lsh1 + 414);
    const auto *lsh1_416 = buffer.data(lsh1 + 416);
    const auto *lsh1_417 = buffer.data(lsh1 + 417);
    const auto *lsh1_418 = buffer.data(lsh1 + 418);
    const auto *lsh1_419 = buffer.data(lsh1 + 419);
    const auto *lsh1_420 = buffer.data(lsh1 + 420);
    const auto *lsh1_421 = buffer.data(lsh1 + 421);
    const auto *lsh1_422 = buffer.data(lsh1 + 422);
    const auto *lsh1_423 = buffer.data(lsh1 + 423);
    const auto *lsh1_424 = buffer.data(lsh1 + 424);
    const auto *lsh1_425 = buffer.data(lsh1 + 425);
    const auto *lsh1_426 = buffer.data(lsh1 + 426);
    const auto *lsh1_427 = buffer.data(lsh1 + 427);
    const auto *lsh1_428 = buffer.data(lsh1 + 428);
    const auto *lsh1_429 = buffer.data(lsh1 + 429);
    const auto *lsh1_434 = buffer.data(lsh1 + 434);
    const auto *lsh1_435 = buffer.data(lsh1 + 435);
    const auto *lsh1_436 = buffer.data(lsh1 + 436);
    const auto *lsh1_437 = buffer.data(lsh1 + 437);
    const auto *lsh1_438 = buffer.data(lsh1 + 438);
    const auto *lsh1_439 = buffer.data(lsh1 + 439);
    const auto *lsh1_440 = buffer.data(lsh1 + 440);
    const auto *lsh1_441 = buffer.data(lsh1 + 441);
    const auto *lsh1_443 = buffer.data(lsh1 + 443);
    const auto *lsh1_444 = buffer.data(lsh1 + 444);
    const auto *lsh1_446 = buffer.data(lsh1 + 446);
    const auto *lsh1_447 = buffer.data(lsh1 + 447);
    const auto *lsh1_448 = buffer.data(lsh1 + 448);
    const auto *lsh1_450 = buffer.data(lsh1 + 450);
    const auto *lsh1_451 = buffer.data(lsh1 + 451);
    const auto *lsh1_456 = buffer.data(lsh1 + 456);
    const auto *lsh1_457 = buffer.data(lsh1 + 457);
    const auto *lsh1_458 = buffer.data(lsh1 + 458);
    const auto *lsh1_459 = buffer.data(lsh1 + 459);
    const auto *lsh1_461 = buffer.data(lsh1 + 461);
    const auto *lsh1_467 = buffer.data(lsh1 + 467);
    const auto *lsh1_471 = buffer.data(lsh1 + 471);
    const auto *lsh1_476 = buffer.data(lsh1 + 476);
    const auto *lsh1_479 = buffer.data(lsh1 + 479);
    const auto *lsh1_480 = buffer.data(lsh1 + 480);
    const auto *lsh1_481 = buffer.data(lsh1 + 481);
    const auto *lsh1_482 = buffer.data(lsh1 + 482);

    const auto *lsi_553 = buffer.data(lsi + 553);
    const auto *lsi_555 = buffer.data(lsi + 555);
    const auto *lsi_556 = buffer.data(lsi + 556);
    const auto *lsi_557 = buffer.data(lsi + 557);
    const auto *lsi_558 = buffer.data(lsi + 558);
    const auto *lsi_559 = buffer.data(lsi + 559);
    const auto *lsi_560 = buffer.data(lsi + 560);
    const auto *lsi_561 = buffer.data(lsi + 561);
    const auto *lsi_562 = buffer.data(lsi + 562);
    const auto *lsi_563 = buffer.data(lsi + 563);
    const auto *lsi_564 = buffer.data(lsi + 564);
    const auto *lsi_565 = buffer.data(lsi + 565);
    const auto *lsi_566 = buffer.data(lsi + 566);
    const auto *lsi_567 = buffer.data(lsi + 567);
    const auto *lsi_568 = buffer.data(lsi + 568);
    const auto *lsi_569 = buffer.data(lsi + 569);
    const auto *lsi_570 = buffer.data(lsi + 570);
    const auto *lsi_571 = buffer.data(lsi + 571);
    const auto *lsi_572 = buffer.data(lsi + 572);
    const auto *lsi_573 = buffer.data(lsi + 573);
    const auto *lsi_574 = buffer.data(lsi + 574);
    const auto *lsi_580 = buffer.data(lsi + 580);
    const auto *lsi_581 = buffer.data(lsi + 581);
    const auto *lsi_582 = buffer.data(lsi + 582);
    const auto *lsi_583 = buffer.data(lsi + 583);
    const auto *lsi_584 = buffer.data(lsi + 584);
    const auto *lsi_585 = buffer.data(lsi + 585);
    const auto *lsi_586 = buffer.data(lsi + 586);
    const auto *lsi_587 = buffer.data(lsi + 587);
    const auto *lsi_588 = buffer.data(lsi + 588);
    const auto *lsi_589 = buffer.data(lsi + 589);
    const auto *lsi_590 = buffer.data(lsi + 590);
    const auto *lsi_591 = buffer.data(lsi + 591);
    const auto *lsi_593 = buffer.data(lsi + 593);
    const auto *lsi_594 = buffer.data(lsi + 594);
    const auto *lsi_595 = buffer.data(lsi + 595);
    const auto *lsi_597 = buffer.data(lsi + 597);
    const auto *lsi_598 = buffer.data(lsi + 598);
    const auto *lsi_599 = buffer.data(lsi + 599);
    const auto *lsi_600 = buffer.data(lsi + 600);
    const auto *lsi_602 = buffer.data(lsi + 602);
    const auto *lsi_603 = buffer.data(lsi + 603);
    const auto *lsi_609 = buffer.data(lsi + 609);
    const auto *lsi_610 = buffer.data(lsi + 610);
    const auto *lsi_611 = buffer.data(lsi + 611);
    const auto *lsi_612 = buffer.data(lsi + 612);
    const auto *lsi_613 = buffer.data(lsi + 613);
    const auto *lsi_614 = buffer.data(lsi + 614);
    const auto *lsi_615 = buffer.data(lsi + 615);
    const auto *lsi_616 = buffer.data(lsi + 616);
    const auto *lsi_618 = buffer.data(lsi + 618);
    const auto *lsi_619 = buffer.data(lsi + 619);
    const auto *lsi_621 = buffer.data(lsi + 621);
    const auto *lsi_622 = buffer.data(lsi + 622);
    const auto *lsi_625 = buffer.data(lsi + 625);
    const auto *lsi_626 = buffer.data(lsi + 626);
    const auto *lsi_630 = buffer.data(lsi + 630);
    const auto *lsi_636 = buffer.data(lsi + 636);
    const auto *lsi_637 = buffer.data(lsi + 637);
    const auto *lsi_638 = buffer.data(lsi + 638);
    const auto *lsi_639 = buffer.data(lsi + 639);
    const auto *lsi_640 = buffer.data(lsi + 640);
    const auto *lsi_641 = buffer.data(lsi + 641);
    const auto *lsi_642 = buffer.data(lsi + 642);
    const auto *lsi_643 = buffer.data(lsi + 643);

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, pc_y, pc_z, ksi_385, ksi_413, \
                         ksi_558, ksi_559, lsh0_414, lsh1_414, lsi_553, lsi_558, \
                         lsi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_15 * ksi_558[k]
                   + f_3 * pc_x[k] * lsi_558[k];

        t_711[k] = f_15 * ksi_559[k]
                   + f_3 * pc_x[k] * lsi_559[k];

        t_712[k] = f_13 * ksi_413[k]
                   + f_1 * lsh0_414[k]
                   - f_2 * lsh1_414[k]
                   + f_3 * pc_y[k] * lsi_553[k];

        t_713[k] = f_16 * ksi_385[k]
                   + f_3 * pc_z[k] * lsi_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, ksi_415, ksi_416, ksi_417, lsh0_416, \
                         lsh0_417, lsh0_418, lsh1_416, lsh1_417, lsh1_418, lsi_555, lsi_556, \
                         lsi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_13 * ksi_415[k]
                   + f_10 * lsh0_416[k]
                   - f_11 * lsh1_416[k]
                   + f_3 * pc_y[k] * lsi_555[k];

        t_715[k] = f_13 * ksi_416[k]
                   + f_8 * lsh0_417[k]
                   - f_9 * lsh1_417[k]
                   + f_3 * pc_y[k] * lsi_556[k];

        t_716[k] = f_13 * ksi_417[k]
                   + f_6 * lsh0_418[k]
                   - f_7 * lsh1_418[k]
                   + f_3 * pc_y[k] * lsi_557[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pa_y, pc_y, ksk0_539, ksi_418, ksi_419, \
                         ksk1_539, lsh0_419, lsh1_419, lsi_558, \
                         lsi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_13 * ksi_418[k]
                   + f_4 * lsh0_419[k]
                   - f_5 * lsh1_419[k]
                   + f_3 * pc_y[k] * lsi_558[k];

        t_718[k] = f_13 * ksi_419[k]
                   + f_3 * pc_y[k] * lsi_559[k];

        t_719[k] = pa_y[k] * ksk0_539[k]
                   - f_12 * pc_y[k] * ksk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, pc_x, pc_y, pc_z, ksi_392, \
                         ksi_560, lsh0_420, lsh1_420, lsi_560, lsi_561, \
                         lsi_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_15 * ksi_560[k]
                   + f_1 * lsh0_420[k]
                   - f_2 * lsh1_420[k]
                   + f_3 * pc_x[k] * lsi_560[k];

        t_721[k] = f_3 * pc_y[k] * lsi_560[k];

        t_722[k] = f_17 * ksi_392[k]
                   + f_3 * pc_z[k] * lsi_560[k];

        t_723[k] = f_4 * lsh0_420[k]
                   - f_5 * lsh1_420[k]
                   + f_3 * pc_y[k] * lsi_561[k];

        t_724[k] = f_3 * pc_y[k] * lsi_562[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pc_x, pc_y, ksi_565, lsh0_421, lsh0_422, \
                         lsh0_425, lsh1_421, lsh1_422, lsh1_425, lsi_563, lsi_564, \
                         lsi_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_15 * ksi_565[k]
                   + f_10 * lsh0_425[k]
                   - f_11 * lsh1_425[k]
                   + f_3 * pc_x[k] * lsi_565[k];

        t_726[k] = f_6 * lsh0_421[k]
                   - f_7 * lsh1_421[k]
                   + f_3 * pc_y[k] * lsi_563[k];

        t_727[k] = f_4 * lsh0_422[k]
                   - f_5 * lsh1_422[k]
                   + f_3 * pc_y[k] * lsi_564[k];

        t_728[k] = f_3 * pc_y[k] * lsi_565[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_y, ksi_569, lsh0_423, lsh0_424, \
                         lsh0_429, lsh1_423, lsh1_424, lsh1_429, lsi_566, lsi_567, \
                         lsi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_15 * ksi_569[k]
                   + f_8 * lsh0_429[k]
                   - f_9 * lsh1_429[k]
                   + f_3 * pc_x[k] * lsi_569[k];

        t_730[k] = f_8 * lsh0_423[k]
                   - f_9 * lsh1_423[k]
                   + f_3 * pc_y[k] * lsi_566[k];

        t_731[k] = f_6 * lsh0_424[k]
                   - f_7 * lsh1_424[k]
                   + f_3 * pc_y[k] * lsi_567[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pc_x, pc_y, ksi_574, lsh0_425, lsh0_434, \
                         lsh1_425, lsh1_434, lsi_568, lsi_569, \
                         lsi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_4 * lsh0_425[k]
                   - f_5 * lsh1_425[k]
                   + f_3 * pc_y[k] * lsi_568[k];

        t_733[k] = f_3 * pc_y[k] * lsi_569[k];

        t_734[k] = f_15 * ksi_574[k]
                   + f_6 * lsh0_434[k]
                   - f_7 * lsh1_434[k]
                   + f_3 * pc_x[k] * lsi_574[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pc_y, lsh0_426, lsh0_427, lsh0_428, lsh1_426, \
                         lsh1_427, lsh1_428, lsi_570, lsi_571, \
                         lsi_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_10 * lsh0_426[k]
                   - f_11 * lsh1_426[k]
                   + f_3 * pc_y[k] * lsi_570[k];

        t_736[k] = f_8 * lsh0_427[k]
                   - f_9 * lsh1_427[k]
                   + f_3 * pc_y[k] * lsi_571[k];

        t_737[k] = f_6 * lsh0_428[k]
                   - f_7 * lsh1_428[k]
                   + f_3 * pc_y[k] * lsi_572[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, pc_x, pc_y, ksi_580, ksi_581, lsh0_429, \
                         lsh0_440, lsh1_429, lsh1_440, lsi_573, lsi_574, lsi_580, \
                         lsi_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_4 * lsh0_429[k]
                   - f_5 * lsh1_429[k]
                   + f_3 * pc_y[k] * lsi_573[k];

        t_739[k] = f_3 * pc_y[k] * lsi_574[k];

        t_740[k] = f_15 * ksi_580[k]
                   + f_4 * lsh0_440[k]
                   - f_5 * lsh1_440[k]
                   + f_3 * pc_x[k] * lsi_580[k];

        t_741[k] = f_15 * ksi_581[k]
                   + f_3 * pc_x[k] * lsi_581[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, pc_x, pc_y, ksi_582, ksi_583, \
                         ksi_584, ksi_585, lsi_580, lsi_582, lsi_583, lsi_584, \
                         lsi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_15 * ksi_582[k]
                   + f_3 * pc_x[k] * lsi_582[k];

        t_743[k] = f_15 * ksi_583[k]
                   + f_3 * pc_x[k] * lsi_583[k];

        t_744[k] = f_15 * ksi_584[k]
                   + f_3 * pc_x[k] * lsi_584[k];

        t_745[k] = f_15 * ksi_585[k]
                   + f_3 * pc_x[k] * lsi_585[k];

        t_746[k] = f_3 * pc_y[k] * lsi_580[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pc_x, pc_y, ksi_587, lsh0_435, lsh0_436, \
                         lsh1_435, lsh1_436, lsi_581, lsi_582, \
                         lsi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_15 * ksi_587[k]
                   + f_3 * pc_x[k] * lsi_587[k];

        t_748[k] = f_1 * lsh0_435[k]
                   - f_2 * lsh1_435[k]
                   + f_3 * pc_y[k] * lsi_581[k];

        t_749[k] = f_19 * lsh0_436[k]
                   - f_20 * lsh1_436[k]
                   + f_3 * pc_y[k] * lsi_582[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_y, lsh0_437, lsh0_438, lsh0_439, lsh1_437, \
                         lsh1_438, lsh1_439, lsi_583, lsi_584, \
                         lsi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_10 * lsh0_437[k]
                   - f_11 * lsh1_437[k]
                   + f_3 * pc_y[k] * lsi_583[k];

        t_751[k] = f_8 * lsh0_438[k]
                   - f_9 * lsh1_438[k]
                   + f_3 * pc_y[k] * lsi_584[k];

        t_752[k] = f_6 * lsh0_439[k]
                   - f_7 * lsh1_439[k]
                   + f_3 * pc_y[k] * lsi_585[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, pc_x, pc_y, pc_z, ksi_419, ksi_588, \
                         lsh0_440, lsh0_441, lsh1_440, lsh1_441, lsi_586, lsi_587, \
                         lsi_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_4 * lsh0_440[k]
                   - f_5 * lsh1_440[k]
                   + f_3 * pc_y[k] * lsi_586[k];

        t_754[k] = f_3 * pc_y[k] * lsi_587[k];

        t_755[k] = f_17 * ksi_419[k]
                   + f_1 * lsh0_440[k]
                   - f_2 * lsh1_440[k]
                   + f_3 * pc_z[k] * lsi_587[k];

        t_756[k] = f_14 * ksi_588[k]
                   + f_1 * lsh0_441[k]
                   - f_2 * lsh1_441[k]
                   + f_3 * pc_x[k] * lsi_588[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pc_x, pc_y, pc_z, ksi_420, ksi_591, \
                         lsh0_444, lsh1_444, lsi_588, lsi_589, \
                         lsi_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_21 * ksi_420[k]
                   + f_3 * pc_y[k] * lsi_588[k];

        t_758[k] = f_3 * pc_z[k] * lsi_588[k];

        t_759[k] = f_14 * ksi_591[k]
                   + f_10 * lsh0_444[k]
                   - f_11 * lsh1_444[k]
                   + f_3 * pc_x[k] * lsi_591[k];

        t_760[k] = f_3 * pc_z[k] * lsi_589[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_x, pc_z, ksi_594, lsh0_441, lsh0_447, \
                         lsh1_441, lsh1_447, lsi_590, lsi_591, \
                         lsi_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_4 * lsh0_441[k]
                   - f_5 * lsh1_441[k]
                   + f_3 * pc_z[k] * lsi_590[k];

        t_762[k] = f_14 * ksi_594[k]
                   + f_8 * lsh0_447[k]
                   - f_9 * lsh1_447[k]
                   + f_3 * pc_x[k] * lsi_594[k];

        t_763[k] = f_3 * pc_z[k] * lsi_591[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, ksi_425, ksi_598, \
                         lsh0_443, lsh0_451, lsh1_443, lsh1_451, lsi_593, lsi_594, \
                         lsi_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_21 * ksi_425[k]
                   + f_3 * pc_y[k] * lsi_593[k];

        t_765[k] = f_6 * lsh0_443[k]
                   - f_7 * lsh1_443[k]
                   + f_3 * pc_z[k] * lsi_593[k];

        t_766[k] = f_14 * ksi_598[k]
                   + f_6 * lsh0_451[k]
                   - f_7 * lsh1_451[k]
                   + f_3 * pc_x[k] * lsi_598[k];

        t_767[k] = f_3 * pc_z[k] * lsi_594[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_y, pc_z, ksi_429, lsh0_444, lsh0_446, \
                         lsh1_444, lsh1_446, lsi_595, lsi_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * lsh0_444[k]
                   - f_5 * lsh1_444[k]
                   + f_3 * pc_z[k] * lsi_595[k];

        t_769[k] = f_21 * ksi_429[k]
                   + f_3 * pc_y[k] * lsi_597[k];

        t_770[k] = f_8 * lsh0_446[k]
                   - f_9 * lsh1_446[k]
                   + f_3 * pc_z[k] * lsi_597[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, pc_x, pc_z, ksi_603, lsh0_447, lsh0_456, \
                         lsh1_447, lsh1_456, lsi_598, lsi_599, \
                         lsi_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_14 * ksi_603[k]
                   + f_4 * lsh0_456[k]
                   - f_5 * lsh1_456[k]
                   + f_3 * pc_x[k] * lsi_603[k];

        t_772[k] = f_3 * pc_z[k] * lsi_598[k];

        t_773[k] = f_4 * lsh0_447[k]
                   - f_5 * lsh1_447[k]
                   + f_3 * pc_z[k] * lsi_599[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pc_x, pc_y, pc_z, ksi_434, ksi_609, \
                         lsh0_448, lsh0_450, lsh1_448, lsh1_450, lsi_600, lsi_602, \
                         lsi_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_6 * lsh0_448[k]
                   - f_7 * lsh1_448[k]
                   + f_3 * pc_z[k] * lsi_600[k];

        t_775[k] = f_21 * ksi_434[k]
                   + f_3 * pc_y[k] * lsi_602[k];

        t_776[k] = f_10 * lsh0_450[k]
                   - f_11 * lsh1_450[k]
                   + f_3 * pc_z[k] * lsi_602[k];

        t_777[k] = f_14 * ksi_609[k]
                   + f_3 * pc_x[k] * lsi_609[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, pc_x, pc_z, ksi_611, ksi_612, \
                         ksi_613, ksi_614, lsi_603, lsi_611, lsi_612, lsi_613, \
                         lsi_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_3 * pc_z[k] * lsi_603[k];

        t_779[k] = f_14 * ksi_611[k]
                   + f_3 * pc_x[k] * lsi_611[k];

        t_780[k] = f_14 * ksi_612[k]
                   + f_3 * pc_x[k] * lsi_612[k];

        t_781[k] = f_14 * ksi_613[k]
                   + f_3 * pc_x[k] * lsi_613[k];

        t_782[k] = f_14 * ksi_614[k]
                   + f_3 * pc_x[k] * lsi_614[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, t_786, pc_x, pc_y, pc_z, ksi_441, ksi_615, \
                         lsh0_456, lsh1_456, lsi_609, lsi_610, \
                         lsi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_14 * ksi_615[k]
                   + f_3 * pc_x[k] * lsi_615[k];

        t_784[k] = f_21 * ksi_441[k]
                   + f_1 * lsh0_456[k]
                   - f_2 * lsh1_456[k]
                   + f_3 * pc_y[k] * lsi_609[k];

        t_785[k] = f_3 * pc_z[k] * lsi_609[k];

        t_786[k] = f_4 * lsh0_456[k]
                   - f_5 * lsh1_456[k]
                   + f_3 * pc_z[k] * lsi_610[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_z, lsh0_457, lsh0_458, lsh0_459, lsh1_457, \
                         lsh1_458, lsh1_459, lsi_611, lsi_612, \
                         lsi_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_6 * lsh0_457[k]
                   - f_7 * lsh1_457[k]
                   + f_3 * pc_z[k] * lsi_611[k];

        t_788[k] = f_8 * lsh0_458[k]
                   - f_9 * lsh1_458[k]
                   + f_3 * pc_z[k] * lsi_612[k];

        t_789[k] = f_10 * lsh0_459[k]
                   - f_11 * lsh1_459[k]
                   + f_3 * pc_z[k] * lsi_613[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_z, pc_y, pc_z, ksk0_540, ksi_447, \
                         ksi_448, ksk1_540, lsh0_461, lsh1_461, lsi_615, \
                         lsi_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_21 * ksi_447[k]
                   + f_3 * pc_y[k] * lsi_615[k];

        t_791[k] = f_1 * lsh0_461[k]
                   - f_2 * lsh1_461[k]
                   + f_3 * pc_z[k] * lsi_615[k];

        t_792[k] = pa_z[k] * ksk0_540[k]
                   - f_12 * pc_z[k] * ksk1_540[k];

        t_793[k] = f_17 * ksi_448[k]
                   + f_3 * pc_y[k] * lsi_616[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pa_z, pc_y, pc_z, ksk0_543, ksi_420, ksi_450, \
                         ksk1_543, lsi_616, lsi_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_13 * ksi_420[k]
                   + f_3 * pc_z[k] * lsi_616[k];

        t_795[k] = pa_z[k] * ksk0_543[k]
                   - f_12 * pc_z[k] * ksk1_543[k];

        t_796[k] = f_17 * ksi_450[k]
                   + f_3 * pc_y[k] * lsi_618[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pa_z, pc_x, pc_z, ksk0_546, ksi_423, ksi_621, \
                         ksk1_546, lsh0_467, lsh1_467, lsi_619, \
                         lsi_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_14 * ksi_621[k]
                   + f_10 * lsh0_467[k]
                   - f_11 * lsh1_467[k]
                   + f_3 * pc_x[k] * lsi_621[k];

        t_798[k] = pa_z[k] * ksk0_546[k]
                   - f_12 * pc_z[k] * ksk1_546[k];

        t_799[k] = f_13 * ksi_423[k]
                   + f_3 * pc_z[k] * lsi_619[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pa_z, pc_x, pc_y, pc_z, ksk0_550, ksi_453, \
                         ksi_625, ksk1_550, lsh0_471, lsh1_471, lsi_621, \
                         lsi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_17 * ksi_453[k]
                   + f_3 * pc_y[k] * lsi_621[k];

        t_801[k] = f_14 * ksi_625[k]
                   + f_8 * lsh0_471[k]
                   - f_9 * lsh1_471[k]
                   + f_3 * pc_x[k] * lsi_625[k];

        t_802[k] = pa_z[k] * ksk0_550[k]
                   - f_12 * pc_z[k] * ksk1_550[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pa_z, pc_y, pc_z, ksk0_552, ksi_426, ksi_427, \
                         ksi_457, ksk1_552, lsi_622, lsi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_13 * ksi_426[k]
                   + f_3 * pc_z[k] * lsi_622[k];

        t_804[k] = pa_z[k] * ksk0_552[k]
                   + f_14 * ksi_427[k]
                   - f_12 * pc_z[k] * ksk1_552[k];

        t_805[k] = f_17 * ksi_457[k]
                   + f_3 * pc_y[k] * lsi_625[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pa_z, pc_x, pc_z, ksk0_555, ksi_430, ksi_630, \
                         ksk1_555, lsh0_476, lsh1_476, lsi_626, \
                         lsi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_14 * ksi_630[k]
                   + f_6 * lsh0_476[k]
                   - f_7 * lsh1_476[k]
                   + f_3 * pc_x[k] * lsi_630[k];

        t_807[k] = pa_z[k] * ksk0_555[k]
                   - f_12 * pc_z[k] * ksk1_555[k];

        t_808[k] = f_13 * ksi_430[k]
                   + f_3 * pc_z[k] * lsi_626[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pa_z, pc_y, pc_z, ksk0_557, ksk0_558, ksi_431, \
                         ksi_432, ksi_462, ksk1_557, ksk1_558, \
                         lsi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pa_z[k] * ksk0_557[k]
                   + f_14 * ksi_431[k]
                   - f_12 * pc_z[k] * ksk1_557[k];

        t_810[k] = pa_z[k] * ksk0_558[k]
                   + f_15 * ksi_432[k]
                   - f_12 * pc_z[k] * ksk1_558[k];

        t_811[k] = f_17 * ksi_462[k]
                   + f_3 * pc_y[k] * lsi_630[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pc_x, ksi_636, ksi_637, ksi_638, ksi_639, \
                         lsh0_482, lsh1_482, lsi_636, lsi_637, lsi_638, \
                         lsi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_14 * ksi_636[k]
                   + f_4 * lsh0_482[k]
                   - f_5 * lsh1_482[k]
                   + f_3 * pc_x[k] * lsi_636[k];

        t_813[k] = f_14 * ksi_637[k]
                   + f_3 * pc_x[k] * lsi_637[k];

        t_814[k] = f_14 * ksi_638[k]
                   + f_3 * pc_x[k] * lsi_638[k];

        t_815[k] = f_14 * ksi_639[k]
                   + f_3 * pc_x[k] * lsi_639[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, pc_x, ksi_640, ksi_641, ksi_642, ksi_643, \
                         lsi_640, lsi_641, lsi_642, lsi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_14 * ksi_640[k]
                   + f_3 * pc_x[k] * lsi_640[k];

        t_817[k] = f_14 * ksi_641[k]
                   + f_3 * pc_x[k] * lsi_641[k];

        t_818[k] = f_14 * ksi_642[k]
                   + f_3 * pc_x[k] * lsi_642[k];

        t_819[k] = f_14 * ksi_643[k]
                   + f_3 * pc_x[k] * lsi_643[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, pa_z, pc_y, pc_z, ksk0_568, ksi_441, ksi_471, \
                         ksk1_568, lsh0_479, lsh1_479, lsi_637, \
                         lsi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = pa_z[k] * ksk0_568[k]
                   - f_12 * pc_z[k] * ksk1_568[k];

        t_821[k] = f_13 * ksi_441[k]
                   + f_3 * pc_z[k] * lsi_637[k];

        t_822[k] = f_17 * ksi_471[k]
                   + f_10 * lsh0_479[k]
                   - f_11 * lsh1_479[k]
                   + f_3 * pc_y[k] * lsi_639[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pc_y, ksi_472, ksi_473, ksi_474, lsh0_480, \
                         lsh0_481, lsh0_482, lsh1_480, lsh1_481, lsh1_482, lsi_640, lsi_641, \
                         lsi_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_17 * ksi_472[k]
                   + f_8 * lsh0_480[k]
                   - f_9 * lsh1_480[k]
                   + f_3 * pc_y[k] * lsi_640[k];

        t_824[k] = f_17 * ksi_473[k]
                   + f_6 * lsh0_481[k]
                   - f_7 * lsh1_481[k]
                   + f_3 * pc_y[k] * lsi_641[k];

        t_825[k] = f_17 * ksi_474[k]
                   + f_4 * lsh0_482[k]
                   - f_5 * lsh1_482[k]
                   + f_3 * pc_y[k] * lsi_642[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t ksi, const size_t lsh0,
                                                          const size_t lsh1, const size_t lsi,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
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
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi_447 = buffer.data(ksi + 447);
    const auto *ksi_448 = buffer.data(ksi + 448);
    const auto *ksi_451 = buffer.data(ksi + 451);
    const auto *ksi_454 = buffer.data(ksi + 454);
    const auto *ksi_458 = buffer.data(ksi + 458);
    const auto *ksi_469 = buffer.data(ksi + 469);
    const auto *ksi_475 = buffer.data(ksi + 475);
    const auto *ksi_476 = buffer.data(ksi + 476);
    const auto *ksi_478 = buffer.data(ksi + 478);
    const auto *ksi_479 = buffer.data(ksi + 479);
    const auto *ksi_481 = buffer.data(ksi + 481);
    const auto *ksi_482 = buffer.data(ksi + 482);
    const auto *ksi_485 = buffer.data(ksi + 485);
    const auto *ksi_486 = buffer.data(ksi + 486);
    const auto *ksi_490 = buffer.data(ksi + 490);
    const auto *ksi_497 = buffer.data(ksi + 497);
    const auto *ksi_499 = buffer.data(ksi + 499);
    const auto *ksi_500 = buffer.data(ksi + 500);
    const auto *ksi_501 = buffer.data(ksi + 501);
    const auto *ksi_502 = buffer.data(ksi + 502);
    const auto *ksi_503 = buffer.data(ksi + 503);
    const auto *ksi_504 = buffer.data(ksi + 504);
    const auto *ksi_506 = buffer.data(ksi + 506);
    const auto *ksi_507 = buffer.data(ksi + 507);
    const auto *ksi_509 = buffer.data(ksi + 509);
    const auto *ksi_510 = buffer.data(ksi + 510);
    const auto *ksi_513 = buffer.data(ksi + 513);
    const auto *ksi_514 = buffer.data(ksi + 514);
    const auto *ksi_518 = buffer.data(ksi + 518);
    const auto *ksi_525 = buffer.data(ksi + 525);
    const auto *ksi_527 = buffer.data(ksi + 527);
    const auto *ksi_528 = buffer.data(ksi + 528);
    const auto *ksi_529 = buffer.data(ksi + 529);
    const auto *ksi_530 = buffer.data(ksi + 530);
    const auto *ksi_531 = buffer.data(ksi + 531);
    const auto *ksi_532 = buffer.data(ksi + 532);
    const auto *ksi_534 = buffer.data(ksi + 534);
    const auto *ksi_537 = buffer.data(ksi + 537);
    const auto *ksi_541 = buffer.data(ksi + 541);
    const auto *ksi_546 = buffer.data(ksi + 546);
    const auto *ksi_553 = buffer.data(ksi + 553);
    const auto *ksi_555 = buffer.data(ksi + 555);
    const auto *ksi_644 = buffer.data(ksi + 644);
    const auto *ksi_647 = buffer.data(ksi + 647);
    const auto *ksi_649 = buffer.data(ksi + 649);
    const auto *ksi_650 = buffer.data(ksi + 650);
    const auto *ksi_653 = buffer.data(ksi + 653);
    const auto *ksi_654 = buffer.data(ksi + 654);
    const auto *ksi_656 = buffer.data(ksi + 656);
    const auto *ksi_658 = buffer.data(ksi + 658);
    const auto *ksi_659 = buffer.data(ksi + 659);
    const auto *ksi_661 = buffer.data(ksi + 661);
    const auto *ksi_662 = buffer.data(ksi + 662);
    const auto *ksi_664 = buffer.data(ksi + 664);
    const auto *ksi_665 = buffer.data(ksi + 665);
    const auto *ksi_666 = buffer.data(ksi + 666);
    const auto *ksi_667 = buffer.data(ksi + 667);
    const auto *ksi_668 = buffer.data(ksi + 668);
    const auto *ksi_669 = buffer.data(ksi + 669);
    const auto *ksi_670 = buffer.data(ksi + 670);
    const auto *ksi_671 = buffer.data(ksi + 671);
    const auto *ksi_672 = buffer.data(ksi + 672);
    const auto *ksi_675 = buffer.data(ksi + 675);
    const auto *ksi_677 = buffer.data(ksi + 677);
    const auto *ksi_678 = buffer.data(ksi + 678);
    const auto *ksi_681 = buffer.data(ksi + 681);
    const auto *ksi_682 = buffer.data(ksi + 682);
    const auto *ksi_684 = buffer.data(ksi + 684);
    const auto *ksi_686 = buffer.data(ksi + 686);
    const auto *ksi_687 = buffer.data(ksi + 687);
    const auto *ksi_689 = buffer.data(ksi + 689);
    const auto *ksi_690 = buffer.data(ksi + 690);
    const auto *ksi_692 = buffer.data(ksi + 692);
    const auto *ksi_693 = buffer.data(ksi + 693);
    const auto *ksi_694 = buffer.data(ksi + 694);
    const auto *ksi_695 = buffer.data(ksi + 695);
    const auto *ksi_696 = buffer.data(ksi + 696);
    const auto *ksi_697 = buffer.data(ksi + 697);
    const auto *ksi_698 = buffer.data(ksi + 698);
    const auto *ksi_699 = buffer.data(ksi + 699);
    const auto *ksi_700 = buffer.data(ksi + 700);
    const auto *ksi_703 = buffer.data(ksi + 703);
    const auto *ksi_705 = buffer.data(ksi + 705);
    const auto *ksi_706 = buffer.data(ksi + 706);
    const auto *ksi_709 = buffer.data(ksi + 709);
    const auto *ksi_710 = buffer.data(ksi + 710);
    const auto *ksi_712 = buffer.data(ksi + 712);
    const auto *ksi_714 = buffer.data(ksi + 714);
    const auto *ksi_715 = buffer.data(ksi + 715);
    const auto *ksi_717 = buffer.data(ksi + 717);
    const auto *ksi_718 = buffer.data(ksi + 718);
    const auto *ksi_720 = buffer.data(ksi + 720);
    const auto *ksi_721 = buffer.data(ksi + 721);
    const auto *ksi_722 = buffer.data(ksi + 722);
    const auto *ksi_723 = buffer.data(ksi + 723);
    const auto *ksi_724 = buffer.data(ksi + 724);
    const auto *ksi_725 = buffer.data(ksi + 725);
    const auto *ksi_726 = buffer.data(ksi + 726);
    const auto *ksi_727 = buffer.data(ksi + 727);

    const auto *lsh0_482 = buffer.data(lsh0 + 482);
    const auto *lsh0_483 = buffer.data(lsh0 + 483);
    const auto *lsh0_486 = buffer.data(lsh0 + 486);
    const auto *lsh0_488 = buffer.data(lsh0 + 488);
    const auto *lsh0_489 = buffer.data(lsh0 + 489);
    const auto *lsh0_492 = buffer.data(lsh0 + 492);
    const auto *lsh0_493 = buffer.data(lsh0 + 493);
    const auto *lsh0_495 = buffer.data(lsh0 + 495);
    const auto *lsh0_497 = buffer.data(lsh0 + 497);
    const auto *lsh0_498 = buffer.data(lsh0 + 498);
    const auto *lsh0_500 = buffer.data(lsh0 + 500);
    const auto *lsh0_501 = buffer.data(lsh0 + 501);
    const auto *lsh0_502 = buffer.data(lsh0 + 502);
    const auto *lsh0_503 = buffer.data(lsh0 + 503);
    const auto *lsh0_504 = buffer.data(lsh0 + 504);
    const auto *lsh0_507 = buffer.data(lsh0 + 507);
    const auto *lsh0_509 = buffer.data(lsh0 + 509);
    const auto *lsh0_510 = buffer.data(lsh0 + 510);
    const auto *lsh0_513 = buffer.data(lsh0 + 513);
    const auto *lsh0_514 = buffer.data(lsh0 + 514);
    const auto *lsh0_516 = buffer.data(lsh0 + 516);
    const auto *lsh0_518 = buffer.data(lsh0 + 518);
    const auto *lsh0_519 = buffer.data(lsh0 + 519);
    const auto *lsh0_521 = buffer.data(lsh0 + 521);
    const auto *lsh0_522 = buffer.data(lsh0 + 522);
    const auto *lsh0_523 = buffer.data(lsh0 + 523);
    const auto *lsh0_524 = buffer.data(lsh0 + 524);
    const auto *lsh0_525 = buffer.data(lsh0 + 525);
    const auto *lsh0_528 = buffer.data(lsh0 + 528);
    const auto *lsh0_530 = buffer.data(lsh0 + 530);
    const auto *lsh0_531 = buffer.data(lsh0 + 531);
    const auto *lsh0_534 = buffer.data(lsh0 + 534);
    const auto *lsh0_535 = buffer.data(lsh0 + 535);
    const auto *lsh0_537 = buffer.data(lsh0 + 537);
    const auto *lsh0_539 = buffer.data(lsh0 + 539);
    const auto *lsh0_540 = buffer.data(lsh0 + 540);
    const auto *lsh0_542 = buffer.data(lsh0 + 542);
    const auto *lsh0_543 = buffer.data(lsh0 + 543);
    const auto *lsh0_545 = buffer.data(lsh0 + 545);

    const auto *lsh1_482 = buffer.data(lsh1 + 482);
    const auto *lsh1_483 = buffer.data(lsh1 + 483);
    const auto *lsh1_486 = buffer.data(lsh1 + 486);
    const auto *lsh1_488 = buffer.data(lsh1 + 488);
    const auto *lsh1_489 = buffer.data(lsh1 + 489);
    const auto *lsh1_492 = buffer.data(lsh1 + 492);
    const auto *lsh1_493 = buffer.data(lsh1 + 493);
    const auto *lsh1_495 = buffer.data(lsh1 + 495);
    const auto *lsh1_497 = buffer.data(lsh1 + 497);
    const auto *lsh1_498 = buffer.data(lsh1 + 498);
    const auto *lsh1_500 = buffer.data(lsh1 + 500);
    const auto *lsh1_501 = buffer.data(lsh1 + 501);
    const auto *lsh1_502 = buffer.data(lsh1 + 502);
    const auto *lsh1_503 = buffer.data(lsh1 + 503);
    const auto *lsh1_504 = buffer.data(lsh1 + 504);
    const auto *lsh1_507 = buffer.data(lsh1 + 507);
    const auto *lsh1_509 = buffer.data(lsh1 + 509);
    const auto *lsh1_510 = buffer.data(lsh1 + 510);
    const auto *lsh1_513 = buffer.data(lsh1 + 513);
    const auto *lsh1_514 = buffer.data(lsh1 + 514);
    const auto *lsh1_516 = buffer.data(lsh1 + 516);
    const auto *lsh1_518 = buffer.data(lsh1 + 518);
    const auto *lsh1_519 = buffer.data(lsh1 + 519);
    const auto *lsh1_521 = buffer.data(lsh1 + 521);
    const auto *lsh1_522 = buffer.data(lsh1 + 522);
    const auto *lsh1_523 = buffer.data(lsh1 + 523);
    const auto *lsh1_524 = buffer.data(lsh1 + 524);
    const auto *lsh1_525 = buffer.data(lsh1 + 525);
    const auto *lsh1_528 = buffer.data(lsh1 + 528);
    const auto *lsh1_530 = buffer.data(lsh1 + 530);
    const auto *lsh1_531 = buffer.data(lsh1 + 531);
    const auto *lsh1_534 = buffer.data(lsh1 + 534);
    const auto *lsh1_535 = buffer.data(lsh1 + 535);
    const auto *lsh1_537 = buffer.data(lsh1 + 537);
    const auto *lsh1_539 = buffer.data(lsh1 + 539);
    const auto *lsh1_540 = buffer.data(lsh1 + 540);
    const auto *lsh1_542 = buffer.data(lsh1 + 542);
    const auto *lsh1_543 = buffer.data(lsh1 + 543);
    const auto *lsh1_545 = buffer.data(lsh1 + 545);

    const auto *lsi_643 = buffer.data(lsi + 643);
    const auto *lsi_644 = buffer.data(lsi + 644);
    const auto *lsi_646 = buffer.data(lsi + 646);
    const auto *lsi_647 = buffer.data(lsi + 647);
    const auto *lsi_649 = buffer.data(lsi + 649);
    const auto *lsi_650 = buffer.data(lsi + 650);
    const auto *lsi_653 = buffer.data(lsi + 653);
    const auto *lsi_654 = buffer.data(lsi + 654);
    const auto *lsi_656 = buffer.data(lsi + 656);
    const auto *lsi_658 = buffer.data(lsi + 658);
    const auto *lsi_659 = buffer.data(lsi + 659);
    const auto *lsi_661 = buffer.data(lsi + 661);
    const auto *lsi_662 = buffer.data(lsi + 662);
    const auto *lsi_664 = buffer.data(lsi + 664);
    const auto *lsi_665 = buffer.data(lsi + 665);
    const auto *lsi_666 = buffer.data(lsi + 666);
    const auto *lsi_667 = buffer.data(lsi + 667);
    const auto *lsi_668 = buffer.data(lsi + 668);
    const auto *lsi_669 = buffer.data(lsi + 669);
    const auto *lsi_670 = buffer.data(lsi + 670);
    const auto *lsi_671 = buffer.data(lsi + 671);
    const auto *lsi_672 = buffer.data(lsi + 672);
    const auto *lsi_674 = buffer.data(lsi + 674);
    const auto *lsi_675 = buffer.data(lsi + 675);
    const auto *lsi_677 = buffer.data(lsi + 677);
    const auto *lsi_678 = buffer.data(lsi + 678);
    const auto *lsi_681 = buffer.data(lsi + 681);
    const auto *lsi_682 = buffer.data(lsi + 682);
    const auto *lsi_684 = buffer.data(lsi + 684);
    const auto *lsi_686 = buffer.data(lsi + 686);
    const auto *lsi_687 = buffer.data(lsi + 687);
    const auto *lsi_689 = buffer.data(lsi + 689);
    const auto *lsi_690 = buffer.data(lsi + 690);
    const auto *lsi_692 = buffer.data(lsi + 692);
    const auto *lsi_693 = buffer.data(lsi + 693);
    const auto *lsi_694 = buffer.data(lsi + 694);
    const auto *lsi_695 = buffer.data(lsi + 695);
    const auto *lsi_696 = buffer.data(lsi + 696);
    const auto *lsi_697 = buffer.data(lsi + 697);
    const auto *lsi_698 = buffer.data(lsi + 698);
    const auto *lsi_699 = buffer.data(lsi + 699);
    const auto *lsi_700 = buffer.data(lsi + 700);
    const auto *lsi_702 = buffer.data(lsi + 702);
    const auto *lsi_703 = buffer.data(lsi + 703);
    const auto *lsi_705 = buffer.data(lsi + 705);
    const auto *lsi_706 = buffer.data(lsi + 706);
    const auto *lsi_709 = buffer.data(lsi + 709);
    const auto *lsi_710 = buffer.data(lsi + 710);
    const auto *lsi_712 = buffer.data(lsi + 712);
    const auto *lsi_714 = buffer.data(lsi + 714);
    const auto *lsi_715 = buffer.data(lsi + 715);
    const auto *lsi_717 = buffer.data(lsi + 717);
    const auto *lsi_718 = buffer.data(lsi + 718);
    const auto *lsi_720 = buffer.data(lsi + 720);
    const auto *lsi_721 = buffer.data(lsi + 721);
    const auto *lsi_722 = buffer.data(lsi + 722);
    const auto *lsi_723 = buffer.data(lsi + 723);
    const auto *lsi_724 = buffer.data(lsi + 724);
    const auto *lsi_725 = buffer.data(lsi + 725);
    const auto *lsi_726 = buffer.data(lsi + 726);
    const auto *lsi_727 = buffer.data(lsi + 727);

#pragma omp simd aligned(t_826, t_827, t_828, pc_x, pc_y, pc_z, ksi_447, ksi_475, ksi_644, \
                         lsh0_482, lsh0_483, lsh1_482, lsh1_483, lsi_643, \
                         lsi_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_17 * ksi_475[k]
                   + f_3 * pc_y[k] * lsi_643[k];

        t_827[k] = f_13 * ksi_447[k]
                   + f_1 * lsh0_482[k]
                   - f_2 * lsh1_482[k]
                   + f_3 * pc_z[k] * lsi_643[k];

        t_828[k] = f_14 * ksi_644[k]
                   + f_1 * lsh0_483[k]
                   - f_2 * lsh1_483[k]
                   + f_3 * pc_x[k] * lsi_644[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pc_x, pc_y, pc_z, ksi_448, ksi_476, \
                         ksi_478, ksi_647, lsh0_486, lsh1_486, lsi_644, lsi_646, \
                         lsi_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_16 * ksi_476[k]
                   + f_3 * pc_y[k] * lsi_644[k];

        t_830[k] = f_14 * ksi_448[k]
                   + f_3 * pc_z[k] * lsi_644[k];

        t_831[k] = f_14 * ksi_647[k]
                   + f_10 * lsh0_486[k]
                   - f_11 * lsh1_486[k]
                   + f_3 * pc_x[k] * lsi_647[k];

        t_832[k] = f_16 * ksi_478[k]
                   + f_3 * pc_y[k] * lsi_646[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, pc_z, ksi_451, ksi_649, ksi_650, lsh0_488, \
                         lsh0_489, lsh1_488, lsh1_489, lsi_647, lsi_649, \
                         lsi_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_14 * ksi_649[k]
                   + f_10 * lsh0_488[k]
                   - f_11 * lsh1_488[k]
                   + f_3 * pc_x[k] * lsi_649[k];

        t_834[k] = f_14 * ksi_650[k]
                   + f_8 * lsh0_489[k]
                   - f_9 * lsh1_489[k]
                   + f_3 * pc_x[k] * lsi_650[k];

        t_835[k] = f_14 * ksi_451[k]
                   + f_3 * pc_z[k] * lsi_647[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pc_x, pc_y, ksi_481, ksi_653, ksi_654, lsh0_492, \
                         lsh0_493, lsh1_492, lsh1_493, lsi_649, lsi_653, \
                         lsi_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * ksi_481[k]
                   + f_3 * pc_y[k] * lsi_649[k];

        t_837[k] = f_14 * ksi_653[k]
                   + f_8 * lsh0_492[k]
                   - f_9 * lsh1_492[k]
                   + f_3 * pc_x[k] * lsi_653[k];

        t_838[k] = f_14 * ksi_654[k]
                   + f_6 * lsh0_493[k]
                   - f_7 * lsh1_493[k]
                   + f_3 * pc_x[k] * lsi_654[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pc_x, pc_y, pc_z, ksi_454, ksi_485, ksi_656, \
                         lsh0_495, lsh1_495, lsi_650, lsi_653, \
                         lsi_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_14 * ksi_454[k]
                   + f_3 * pc_z[k] * lsi_650[k];

        t_840[k] = f_14 * ksi_656[k]
                   + f_6 * lsh0_495[k]
                   - f_7 * lsh1_495[k]
                   + f_3 * pc_x[k] * lsi_656[k];

        t_841[k] = f_16 * ksi_485[k]
                   + f_3 * pc_y[k] * lsi_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, pc_z, ksi_458, ksi_658, ksi_659, lsh0_497, \
                         lsh0_498, lsh1_497, lsh1_498, lsi_654, lsi_658, \
                         lsi_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_14 * ksi_658[k]
                   + f_6 * lsh0_497[k]
                   - f_7 * lsh1_497[k]
                   + f_3 * pc_x[k] * lsi_658[k];

        t_843[k] = f_14 * ksi_659[k]
                   + f_4 * lsh0_498[k]
                   - f_5 * lsh1_498[k]
                   + f_3 * pc_x[k] * lsi_659[k];

        t_844[k] = f_14 * ksi_458[k]
                   + f_3 * pc_z[k] * lsi_654[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, ksi_490, ksi_661, ksi_662, lsh0_500, \
                         lsh0_501, lsh1_500, lsh1_501, lsi_658, lsi_661, \
                         lsi_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_14 * ksi_661[k]
                   + f_4 * lsh0_500[k]
                   - f_5 * lsh1_500[k]
                   + f_3 * pc_x[k] * lsi_661[k];

        t_846[k] = f_14 * ksi_662[k]
                   + f_4 * lsh0_501[k]
                   - f_5 * lsh1_501[k]
                   + f_3 * pc_x[k] * lsi_662[k];

        t_847[k] = f_16 * ksi_490[k]
                   + f_3 * pc_y[k] * lsi_658[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, ksi_664, ksi_665, ksi_666, ksi_667, \
                         lsh0_503, lsh1_503, lsi_664, lsi_665, lsi_666, \
                         lsi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_14 * ksi_664[k]
                   + f_4 * lsh0_503[k]
                   - f_5 * lsh1_503[k]
                   + f_3 * pc_x[k] * lsi_664[k];

        t_849[k] = f_14 * ksi_665[k]
                   + f_3 * pc_x[k] * lsi_665[k];

        t_850[k] = f_14 * ksi_666[k]
                   + f_3 * pc_x[k] * lsi_666[k];

        t_851[k] = f_14 * ksi_667[k]
                   + f_3 * pc_x[k] * lsi_667[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, ksi_668, ksi_669, ksi_670, ksi_671, \
                         lsi_668, lsi_669, lsi_670, lsi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_14 * ksi_668[k]
                   + f_3 * pc_x[k] * lsi_668[k];

        t_853[k] = f_14 * ksi_669[k]
                   + f_3 * pc_x[k] * lsi_669[k];

        t_854[k] = f_14 * ksi_670[k]
                   + f_3 * pc_x[k] * lsi_670[k];

        t_855[k] = f_14 * ksi_671[k]
                   + f_3 * pc_x[k] * lsi_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, ksi_469, ksi_497, ksi_499, lsh0_498, \
                         lsh0_500, lsh1_498, lsh1_500, lsi_665, \
                         lsi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * ksi_497[k]
                   + f_1 * lsh0_498[k]
                   - f_2 * lsh1_498[k]
                   + f_3 * pc_y[k] * lsi_665[k];

        t_857[k] = f_14 * ksi_469[k]
                   + f_3 * pc_z[k] * lsi_665[k];

        t_858[k] = f_16 * ksi_499[k]
                   + f_10 * lsh0_500[k]
                   - f_11 * lsh1_500[k]
                   + f_3 * pc_y[k] * lsi_667[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_y, ksi_500, ksi_501, ksi_502, lsh0_501, \
                         lsh0_502, lsh0_503, lsh1_501, lsh1_502, lsh1_503, lsi_668, lsi_669, \
                         lsi_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * ksi_500[k]
                   + f_8 * lsh0_501[k]
                   - f_9 * lsh1_501[k]
                   + f_3 * pc_y[k] * lsi_668[k];

        t_860[k] = f_16 * ksi_501[k]
                   + f_6 * lsh0_502[k]
                   - f_7 * lsh1_502[k]
                   + f_3 * pc_y[k] * lsi_669[k];

        t_861[k] = f_16 * ksi_502[k]
                   + f_4 * lsh0_503[k]
                   - f_5 * lsh1_503[k]
                   + f_3 * pc_y[k] * lsi_670[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pc_x, pc_y, pc_z, ksi_475, ksi_503, ksi_672, \
                         lsh0_503, lsh0_504, lsh1_503, lsh1_504, lsi_671, \
                         lsi_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * ksi_503[k]
                   + f_3 * pc_y[k] * lsi_671[k];

        t_863[k] = f_14 * ksi_475[k]
                   + f_1 * lsh0_503[k]
                   - f_2 * lsh1_503[k]
                   + f_3 * pc_z[k] * lsi_671[k];

        t_864[k] = f_14 * ksi_672[k]
                   + f_1 * lsh0_504[k]
                   - f_2 * lsh1_504[k]
                   + f_3 * pc_x[k] * lsi_672[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pc_x, pc_y, pc_z, ksi_476, ksi_504, \
                         ksi_506, ksi_675, lsh0_507, lsh1_507, lsi_672, lsi_674, \
                         lsi_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_15 * ksi_504[k]
                   + f_3 * pc_y[k] * lsi_672[k];

        t_866[k] = f_15 * ksi_476[k]
                   + f_3 * pc_z[k] * lsi_672[k];

        t_867[k] = f_14 * ksi_675[k]
                   + f_10 * lsh0_507[k]
                   - f_11 * lsh1_507[k]
                   + f_3 * pc_x[k] * lsi_675[k];

        t_868[k] = f_15 * ksi_506[k]
                   + f_3 * pc_y[k] * lsi_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pc_x, pc_z, ksi_479, ksi_677, ksi_678, lsh0_509, \
                         lsh0_510, lsh1_509, lsh1_510, lsi_675, lsi_677, \
                         lsi_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * ksi_677[k]
                   + f_10 * lsh0_509[k]
                   - f_11 * lsh1_509[k]
                   + f_3 * pc_x[k] * lsi_677[k];

        t_870[k] = f_14 * ksi_678[k]
                   + f_8 * lsh0_510[k]
                   - f_9 * lsh1_510[k]
                   + f_3 * pc_x[k] * lsi_678[k];

        t_871[k] = f_15 * ksi_479[k]
                   + f_3 * pc_z[k] * lsi_675[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pc_x, pc_y, ksi_509, ksi_681, ksi_682, lsh0_513, \
                         lsh0_514, lsh1_513, lsh1_514, lsi_677, lsi_681, \
                         lsi_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_15 * ksi_509[k]
                   + f_3 * pc_y[k] * lsi_677[k];

        t_873[k] = f_14 * ksi_681[k]
                   + f_8 * lsh0_513[k]
                   - f_9 * lsh1_513[k]
                   + f_3 * pc_x[k] * lsi_681[k];

        t_874[k] = f_14 * ksi_682[k]
                   + f_6 * lsh0_514[k]
                   - f_7 * lsh1_514[k]
                   + f_3 * pc_x[k] * lsi_682[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pc_x, pc_y, pc_z, ksi_482, ksi_513, ksi_684, \
                         lsh0_516, lsh1_516, lsi_678, lsi_681, \
                         lsi_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_15 * ksi_482[k]
                   + f_3 * pc_z[k] * lsi_678[k];

        t_876[k] = f_14 * ksi_684[k]
                   + f_6 * lsh0_516[k]
                   - f_7 * lsh1_516[k]
                   + f_3 * pc_x[k] * lsi_684[k];

        t_877[k] = f_15 * ksi_513[k]
                   + f_3 * pc_y[k] * lsi_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_x, pc_z, ksi_486, ksi_686, ksi_687, lsh0_518, \
                         lsh0_519, lsh1_518, lsh1_519, lsi_682, lsi_686, \
                         lsi_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_14 * ksi_686[k]
                   + f_6 * lsh0_518[k]
                   - f_7 * lsh1_518[k]
                   + f_3 * pc_x[k] * lsi_686[k];

        t_879[k] = f_14 * ksi_687[k]
                   + f_4 * lsh0_519[k]
                   - f_5 * lsh1_519[k]
                   + f_3 * pc_x[k] * lsi_687[k];

        t_880[k] = f_15 * ksi_486[k]
                   + f_3 * pc_z[k] * lsi_682[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pc_x, pc_y, ksi_518, ksi_689, ksi_690, lsh0_521, \
                         lsh0_522, lsh1_521, lsh1_522, lsi_686, lsi_689, \
                         lsi_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_14 * ksi_689[k]
                   + f_4 * lsh0_521[k]
                   - f_5 * lsh1_521[k]
                   + f_3 * pc_x[k] * lsi_689[k];

        t_882[k] = f_14 * ksi_690[k]
                   + f_4 * lsh0_522[k]
                   - f_5 * lsh1_522[k]
                   + f_3 * pc_x[k] * lsi_690[k];

        t_883[k] = f_15 * ksi_518[k]
                   + f_3 * pc_y[k] * lsi_686[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, pc_x, ksi_692, ksi_693, ksi_694, ksi_695, \
                         lsh0_524, lsh1_524, lsi_692, lsi_693, lsi_694, \
                         lsi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_14 * ksi_692[k]
                   + f_4 * lsh0_524[k]
                   - f_5 * lsh1_524[k]
                   + f_3 * pc_x[k] * lsi_692[k];

        t_885[k] = f_14 * ksi_693[k]
                   + f_3 * pc_x[k] * lsi_693[k];

        t_886[k] = f_14 * ksi_694[k]
                   + f_3 * pc_x[k] * lsi_694[k];

        t_887[k] = f_14 * ksi_695[k]
                   + f_3 * pc_x[k] * lsi_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pc_x, ksi_696, ksi_697, ksi_698, ksi_699, \
                         lsi_696, lsi_697, lsi_698, lsi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_14 * ksi_696[k]
                   + f_3 * pc_x[k] * lsi_696[k];

        t_889[k] = f_14 * ksi_697[k]
                   + f_3 * pc_x[k] * lsi_697[k];

        t_890[k] = f_14 * ksi_698[k]
                   + f_3 * pc_x[k] * lsi_698[k];

        t_891[k] = f_14 * ksi_699[k]
                   + f_3 * pc_x[k] * lsi_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pc_y, pc_z, ksi_497, ksi_525, ksi_527, lsh0_519, \
                         lsh0_521, lsh1_519, lsh1_521, lsi_693, \
                         lsi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_15 * ksi_525[k]
                   + f_1 * lsh0_519[k]
                   - f_2 * lsh1_519[k]
                   + f_3 * pc_y[k] * lsi_693[k];

        t_893[k] = f_15 * ksi_497[k]
                   + f_3 * pc_z[k] * lsi_693[k];

        t_894[k] = f_15 * ksi_527[k]
                   + f_10 * lsh0_521[k]
                   - f_11 * lsh1_521[k]
                   + f_3 * pc_y[k] * lsi_695[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pc_y, ksi_528, ksi_529, ksi_530, lsh0_522, \
                         lsh0_523, lsh0_524, lsh1_522, lsh1_523, lsh1_524, lsi_696, lsi_697, \
                         lsi_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_15 * ksi_528[k]
                   + f_8 * lsh0_522[k]
                   - f_9 * lsh1_522[k]
                   + f_3 * pc_y[k] * lsi_696[k];

        t_896[k] = f_15 * ksi_529[k]
                   + f_6 * lsh0_523[k]
                   - f_7 * lsh1_523[k]
                   + f_3 * pc_y[k] * lsi_697[k];

        t_897[k] = f_15 * ksi_530[k]
                   + f_4 * lsh0_524[k]
                   - f_5 * lsh1_524[k]
                   + f_3 * pc_y[k] * lsi_698[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pc_x, pc_y, pc_z, ksi_503, ksi_531, ksi_700, \
                         lsh0_524, lsh0_525, lsh1_524, lsh1_525, lsi_699, \
                         lsi_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * ksi_531[k]
                   + f_3 * pc_y[k] * lsi_699[k];

        t_899[k] = f_15 * ksi_503[k]
                   + f_1 * lsh0_524[k]
                   - f_2 * lsh1_524[k]
                   + f_3 * pc_z[k] * lsi_699[k];

        t_900[k] = f_14 * ksi_700[k]
                   + f_1 * lsh0_525[k]
                   - f_2 * lsh1_525[k]
                   + f_3 * pc_x[k] * lsi_700[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pc_x, pc_y, pc_z, ksi_504, ksi_532, \
                         ksi_534, ksi_703, lsh0_528, lsh1_528, lsi_700, lsi_702, \
                         lsi_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_14 * ksi_532[k]
                   + f_3 * pc_y[k] * lsi_700[k];

        t_902[k] = f_16 * ksi_504[k]
                   + f_3 * pc_z[k] * lsi_700[k];

        t_903[k] = f_14 * ksi_703[k]
                   + f_10 * lsh0_528[k]
                   - f_11 * lsh1_528[k]
                   + f_3 * pc_x[k] * lsi_703[k];

        t_904[k] = f_14 * ksi_534[k]
                   + f_3 * pc_y[k] * lsi_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, ksi_507, ksi_705, ksi_706, lsh0_530, \
                         lsh0_531, lsh1_530, lsh1_531, lsi_703, lsi_705, \
                         lsi_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_14 * ksi_705[k]
                   + f_10 * lsh0_530[k]
                   - f_11 * lsh1_530[k]
                   + f_3 * pc_x[k] * lsi_705[k];

        t_906[k] = f_14 * ksi_706[k]
                   + f_8 * lsh0_531[k]
                   - f_9 * lsh1_531[k]
                   + f_3 * pc_x[k] * lsi_706[k];

        t_907[k] = f_16 * ksi_507[k]
                   + f_3 * pc_z[k] * lsi_703[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, ksi_537, ksi_709, ksi_710, lsh0_534, \
                         lsh0_535, lsh1_534, lsh1_535, lsi_705, lsi_709, \
                         lsi_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * ksi_537[k]
                   + f_3 * pc_y[k] * lsi_705[k];

        t_909[k] = f_14 * ksi_709[k]
                   + f_8 * lsh0_534[k]
                   - f_9 * lsh1_534[k]
                   + f_3 * pc_x[k] * lsi_709[k];

        t_910[k] = f_14 * ksi_710[k]
                   + f_6 * lsh0_535[k]
                   - f_7 * lsh1_535[k]
                   + f_3 * pc_x[k] * lsi_710[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pc_x, pc_y, pc_z, ksi_510, ksi_541, ksi_712, \
                         lsh0_537, lsh1_537, lsi_706, lsi_709, \
                         lsi_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_16 * ksi_510[k]
                   + f_3 * pc_z[k] * lsi_706[k];

        t_912[k] = f_14 * ksi_712[k]
                   + f_6 * lsh0_537[k]
                   - f_7 * lsh1_537[k]
                   + f_3 * pc_x[k] * lsi_712[k];

        t_913[k] = f_14 * ksi_541[k]
                   + f_3 * pc_y[k] * lsi_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pc_x, pc_z, ksi_514, ksi_714, ksi_715, lsh0_539, \
                         lsh0_540, lsh1_539, lsh1_540, lsi_710, lsi_714, \
                         lsi_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_14 * ksi_714[k]
                   + f_6 * lsh0_539[k]
                   - f_7 * lsh1_539[k]
                   + f_3 * pc_x[k] * lsi_714[k];

        t_915[k] = f_14 * ksi_715[k]
                   + f_4 * lsh0_540[k]
                   - f_5 * lsh1_540[k]
                   + f_3 * pc_x[k] * lsi_715[k];

        t_916[k] = f_16 * ksi_514[k]
                   + f_3 * pc_z[k] * lsi_710[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pc_x, pc_y, ksi_546, ksi_717, ksi_718, lsh0_542, \
                         lsh0_543, lsh1_542, lsh1_543, lsi_714, lsi_717, \
                         lsi_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_14 * ksi_717[k]
                   + f_4 * lsh0_542[k]
                   - f_5 * lsh1_542[k]
                   + f_3 * pc_x[k] * lsi_717[k];

        t_918[k] = f_14 * ksi_718[k]
                   + f_4 * lsh0_543[k]
                   - f_5 * lsh1_543[k]
                   + f_3 * pc_x[k] * lsi_718[k];

        t_919[k] = f_14 * ksi_546[k]
                   + f_3 * pc_y[k] * lsi_714[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pc_x, ksi_720, ksi_721, ksi_722, ksi_723, \
                         lsh0_545, lsh1_545, lsi_720, lsi_721, lsi_722, \
                         lsi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_14 * ksi_720[k]
                   + f_4 * lsh0_545[k]
                   - f_5 * lsh1_545[k]
                   + f_3 * pc_x[k] * lsi_720[k];

        t_921[k] = f_14 * ksi_721[k]
                   + f_3 * pc_x[k] * lsi_721[k];

        t_922[k] = f_14 * ksi_722[k]
                   + f_3 * pc_x[k] * lsi_722[k];

        t_923[k] = f_14 * ksi_723[k]
                   + f_3 * pc_x[k] * lsi_723[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, ksi_724, ksi_725, ksi_726, ksi_727, \
                         lsi_724, lsi_725, lsi_726, lsi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_14 * ksi_724[k]
                   + f_3 * pc_x[k] * lsi_724[k];

        t_925[k] = f_14 * ksi_725[k]
                   + f_3 * pc_x[k] * lsi_725[k];

        t_926[k] = f_14 * ksi_726[k]
                   + f_3 * pc_x[k] * lsi_726[k];

        t_927[k] = f_14 * ksi_727[k]
                   + f_3 * pc_x[k] * lsi_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_y, pc_z, ksi_525, ksi_553, ksi_555, lsh0_540, \
                         lsh0_542, lsh1_540, lsh1_542, lsi_721, \
                         lsi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_14 * ksi_553[k]
                   + f_1 * lsh0_540[k]
                   - f_2 * lsh1_540[k]
                   + f_3 * pc_y[k] * lsi_721[k];

        t_929[k] = f_16 * ksi_525[k]
                   + f_3 * pc_z[k] * lsi_721[k];

        t_930[k] = f_14 * ksi_555[k]
                   + f_10 * lsh0_542[k]
                   - f_11 * lsh1_542[k]
                   + f_3 * pc_y[k] * lsi_723[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsh0, const size_t lsh1,
                                                          const size_t lsi, const size_t ncols,
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
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.0 / q;

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
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_720 = buffer.data(ksk0 + 720);
    const auto *ksk0_723 = buffer.data(ksk0 + 723);
    const auto *ksk0_725 = buffer.data(ksk0 + 725);
    const auto *ksk0_726 = buffer.data(ksk0 + 726);
    const auto *ksk0_729 = buffer.data(ksk0 + 729);
    const auto *ksk0_730 = buffer.data(ksk0 + 730);
    const auto *ksk0_732 = buffer.data(ksk0 + 732);
    const auto *ksk0_734 = buffer.data(ksk0 + 734);
    const auto *ksk0_735 = buffer.data(ksk0 + 735);
    const auto *ksk0_737 = buffer.data(ksk0 + 737);
    const auto *ksk0_738 = buffer.data(ksk0 + 738);
    const auto *ksk0_740 = buffer.data(ksk0 + 740);
    const auto *ksk0_755 = buffer.data(ksk0 + 755);
    const auto *ksk0_756 = buffer.data(ksk0 + 756);
    const auto *ksk0_759 = buffer.data(ksk0 + 759);
    const auto *ksk0_1008 = buffer.data(ksk0 + 1008);
    const auto *ksk0_1011 = buffer.data(ksk0 + 1011);
    const auto *ksk0_1014 = buffer.data(ksk0 + 1014);
    const auto *ksk0_1018 = buffer.data(ksk0 + 1018);
    const auto *ksk0_1023 = buffer.data(ksk0 + 1023);
    const auto *ksk0_1036 = buffer.data(ksk0 + 1036);
    const auto *ksk0_1038 = buffer.data(ksk0 + 1038);
    const auto *ksk0_1039 = buffer.data(ksk0 + 1039);
    const auto *ksk0_1040 = buffer.data(ksk0 + 1040);
    const auto *ksk0_1041 = buffer.data(ksk0 + 1041);
    const auto *ksk0_1043 = buffer.data(ksk0 + 1043);
    const auto *ksk0_1049 = buffer.data(ksk0 + 1049);

    const auto *ksi_531 = buffer.data(ksi + 531);
    const auto *ksi_532 = buffer.data(ksi + 532);
    const auto *ksi_535 = buffer.data(ksi + 535);
    const auto *ksi_538 = buffer.data(ksi + 538);
    const auto *ksi_542 = buffer.data(ksi + 542);
    const auto *ksi_553 = buffer.data(ksi + 553);
    const auto *ksi_556 = buffer.data(ksi + 556);
    const auto *ksi_557 = buffer.data(ksi + 557);
    const auto *ksi_558 = buffer.data(ksi + 558);
    const auto *ksi_559 = buffer.data(ksi + 559);
    const auto *ksi_560 = buffer.data(ksi + 560);
    const auto *ksi_561 = buffer.data(ksi + 561);
    const auto *ksi_562 = buffer.data(ksi + 562);
    const auto *ksi_563 = buffer.data(ksi + 563);
    const auto *ksi_565 = buffer.data(ksi + 565);
    const auto *ksi_566 = buffer.data(ksi + 566);
    const auto *ksi_568 = buffer.data(ksi + 568);
    const auto *ksi_569 = buffer.data(ksi + 569);
    const auto *ksi_570 = buffer.data(ksi + 570);
    const auto *ksi_572 = buffer.data(ksi + 572);
    const auto *ksi_573 = buffer.data(ksi + 573);
    const auto *ksi_574 = buffer.data(ksi + 574);
    const auto *ksi_581 = buffer.data(ksi + 581);
    const auto *ksi_583 = buffer.data(ksi + 583);
    const auto *ksi_584 = buffer.data(ksi + 584);
    const auto *ksi_585 = buffer.data(ksi + 585);
    const auto *ksi_586 = buffer.data(ksi + 586);
    const auto *ksi_587 = buffer.data(ksi + 587);
    const auto *ksi_588 = buffer.data(ksi + 588);
    const auto *ksi_593 = buffer.data(ksi + 593);
    const auto *ksi_597 = buffer.data(ksi + 597);
    const auto *ksi_602 = buffer.data(ksi + 602);
    const auto *ksi_615 = buffer.data(ksi + 615);
    const auto *ksi_616 = buffer.data(ksi + 616);
    const auto *ksi_618 = buffer.data(ksi + 618);
    const auto *ksi_749 = buffer.data(ksi + 749);
    const auto *ksi_750 = buffer.data(ksi + 750);
    const auto *ksi_751 = buffer.data(ksi + 751);
    const auto *ksi_752 = buffer.data(ksi + 752);
    const auto *ksi_753 = buffer.data(ksi + 753);
    const auto *ksi_754 = buffer.data(ksi + 754);
    const auto *ksi_755 = buffer.data(ksi + 755);
    const auto *ksi_756 = buffer.data(ksi + 756);
    const auto *ksi_761 = buffer.data(ksi + 761);
    const auto *ksi_765 = buffer.data(ksi + 765);
    const auto *ksi_770 = buffer.data(ksi + 770);
    const auto *ksi_776 = buffer.data(ksi + 776);
    const auto *ksi_777 = buffer.data(ksi + 777);
    const auto *ksi_778 = buffer.data(ksi + 778);
    const auto *ksi_779 = buffer.data(ksi + 779);
    const auto *ksi_780 = buffer.data(ksi + 780);
    const auto *ksi_781 = buffer.data(ksi + 781);
    const auto *ksi_783 = buffer.data(ksi + 783);
    const auto *ksi_784 = buffer.data(ksi + 784);
    const auto *ksi_787 = buffer.data(ksi + 787);
    const auto *ksi_790 = buffer.data(ksi + 790);
    const auto *ksi_794 = buffer.data(ksi + 794);
    const auto *ksi_799 = buffer.data(ksi + 799);
    const auto *ksi_805 = buffer.data(ksi + 805);
    const auto *ksi_807 = buffer.data(ksi + 807);
    const auto *ksi_808 = buffer.data(ksi + 808);
    const auto *ksi_809 = buffer.data(ksi + 809);
    const auto *ksi_810 = buffer.data(ksi + 810);
    const auto *ksi_811 = buffer.data(ksi + 811);
    const auto *ksi_817 = buffer.data(ksi + 817);

    const auto *ksk1_720 = buffer.data(ksk1 + 720);
    const auto *ksk1_723 = buffer.data(ksk1 + 723);
    const auto *ksk1_725 = buffer.data(ksk1 + 725);
    const auto *ksk1_726 = buffer.data(ksk1 + 726);
    const auto *ksk1_729 = buffer.data(ksk1 + 729);
    const auto *ksk1_730 = buffer.data(ksk1 + 730);
    const auto *ksk1_732 = buffer.data(ksk1 + 732);
    const auto *ksk1_734 = buffer.data(ksk1 + 734);
    const auto *ksk1_735 = buffer.data(ksk1 + 735);
    const auto *ksk1_737 = buffer.data(ksk1 + 737);
    const auto *ksk1_738 = buffer.data(ksk1 + 738);
    const auto *ksk1_740 = buffer.data(ksk1 + 740);
    const auto *ksk1_755 = buffer.data(ksk1 + 755);
    const auto *ksk1_756 = buffer.data(ksk1 + 756);
    const auto *ksk1_759 = buffer.data(ksk1 + 759);
    const auto *ksk1_1008 = buffer.data(ksk1 + 1008);
    const auto *ksk1_1011 = buffer.data(ksk1 + 1011);
    const auto *ksk1_1014 = buffer.data(ksk1 + 1014);
    const auto *ksk1_1018 = buffer.data(ksk1 + 1018);
    const auto *ksk1_1023 = buffer.data(ksk1 + 1023);
    const auto *ksk1_1036 = buffer.data(ksk1 + 1036);
    const auto *ksk1_1038 = buffer.data(ksk1 + 1038);
    const auto *ksk1_1039 = buffer.data(ksk1 + 1039);
    const auto *ksk1_1040 = buffer.data(ksk1 + 1040);
    const auto *ksk1_1041 = buffer.data(ksk1 + 1041);
    const auto *ksk1_1043 = buffer.data(ksk1 + 1043);
    const auto *ksk1_1049 = buffer.data(ksk1 + 1049);

    const auto *lsh0_543 = buffer.data(lsh0 + 543);
    const auto *lsh0_544 = buffer.data(lsh0 + 544);
    const auto *lsh0_545 = buffer.data(lsh0 + 545);
    const auto *lsh0_561 = buffer.data(lsh0 + 561);
    const auto *lsh0_563 = buffer.data(lsh0 + 563);
    const auto *lsh0_564 = buffer.data(lsh0 + 564);
    const auto *lsh0_565 = buffer.data(lsh0 + 565);
    const auto *lsh0_566 = buffer.data(lsh0 + 566);
    const auto *lsh0_567 = buffer.data(lsh0 + 567);
    const auto *lsh0_568 = buffer.data(lsh0 + 568);
    const auto *lsh0_569 = buffer.data(lsh0 + 569);
    const auto *lsh0_570 = buffer.data(lsh0 + 570);
    const auto *lsh0_571 = buffer.data(lsh0 + 571);
    const auto *lsh0_572 = buffer.data(lsh0 + 572);
    const auto *lsh0_573 = buffer.data(lsh0 + 573);
    const auto *lsh0_574 = buffer.data(lsh0 + 574);
    const auto *lsh0_575 = buffer.data(lsh0 + 575);
    const auto *lsh0_576 = buffer.data(lsh0 + 576);
    const auto *lsh0_581 = buffer.data(lsh0 + 581);
    const auto *lsh0_582 = buffer.data(lsh0 + 582);
    const auto *lsh0_583 = buffer.data(lsh0 + 583);
    const auto *lsh0_584 = buffer.data(lsh0 + 584);
    const auto *lsh0_585 = buffer.data(lsh0 + 585);
    const auto *lsh0_586 = buffer.data(lsh0 + 586);
    const auto *lsh0_587 = buffer.data(lsh0 + 587);
    const auto *lsh0_588 = buffer.data(lsh0 + 588);
    const auto *lsh0_590 = buffer.data(lsh0 + 590);
    const auto *lsh0_591 = buffer.data(lsh0 + 591);
    const auto *lsh0_593 = buffer.data(lsh0 + 593);
    const auto *lsh0_594 = buffer.data(lsh0 + 594);
    const auto *lsh0_595 = buffer.data(lsh0 + 595);
    const auto *lsh0_597 = buffer.data(lsh0 + 597);

    const auto *lsh1_543 = buffer.data(lsh1 + 543);
    const auto *lsh1_544 = buffer.data(lsh1 + 544);
    const auto *lsh1_545 = buffer.data(lsh1 + 545);
    const auto *lsh1_561 = buffer.data(lsh1 + 561);
    const auto *lsh1_563 = buffer.data(lsh1 + 563);
    const auto *lsh1_564 = buffer.data(lsh1 + 564);
    const auto *lsh1_565 = buffer.data(lsh1 + 565);
    const auto *lsh1_566 = buffer.data(lsh1 + 566);
    const auto *lsh1_567 = buffer.data(lsh1 + 567);
    const auto *lsh1_568 = buffer.data(lsh1 + 568);
    const auto *lsh1_569 = buffer.data(lsh1 + 569);
    const auto *lsh1_570 = buffer.data(lsh1 + 570);
    const auto *lsh1_571 = buffer.data(lsh1 + 571);
    const auto *lsh1_572 = buffer.data(lsh1 + 572);
    const auto *lsh1_573 = buffer.data(lsh1 + 573);
    const auto *lsh1_574 = buffer.data(lsh1 + 574);
    const auto *lsh1_575 = buffer.data(lsh1 + 575);
    const auto *lsh1_576 = buffer.data(lsh1 + 576);
    const auto *lsh1_581 = buffer.data(lsh1 + 581);
    const auto *lsh1_582 = buffer.data(lsh1 + 582);
    const auto *lsh1_583 = buffer.data(lsh1 + 583);
    const auto *lsh1_584 = buffer.data(lsh1 + 584);
    const auto *lsh1_585 = buffer.data(lsh1 + 585);
    const auto *lsh1_586 = buffer.data(lsh1 + 586);
    const auto *lsh1_587 = buffer.data(lsh1 + 587);
    const auto *lsh1_588 = buffer.data(lsh1 + 588);
    const auto *lsh1_590 = buffer.data(lsh1 + 590);
    const auto *lsh1_591 = buffer.data(lsh1 + 591);
    const auto *lsh1_593 = buffer.data(lsh1 + 593);
    const auto *lsh1_594 = buffer.data(lsh1 + 594);
    const auto *lsh1_595 = buffer.data(lsh1 + 595);
    const auto *lsh1_597 = buffer.data(lsh1 + 597);

    const auto *lsi_724 = buffer.data(lsi + 724);
    const auto *lsi_725 = buffer.data(lsi + 725);
    const auto *lsi_726 = buffer.data(lsi + 726);
    const auto *lsi_727 = buffer.data(lsi + 727);
    const auto *lsi_728 = buffer.data(lsi + 728);
    const auto *lsi_730 = buffer.data(lsi + 730);
    const auto *lsi_731 = buffer.data(lsi + 731);
    const auto *lsi_733 = buffer.data(lsi + 733);
    const auto *lsi_734 = buffer.data(lsi + 734);
    const auto *lsi_737 = buffer.data(lsi + 737);
    const auto *lsi_738 = buffer.data(lsi + 738);
    const auto *lsi_742 = buffer.data(lsi + 742);
    const auto *lsi_749 = buffer.data(lsi + 749);
    const auto *lsi_750 = buffer.data(lsi + 750);
    const auto *lsi_751 = buffer.data(lsi + 751);
    const auto *lsi_752 = buffer.data(lsi + 752);
    const auto *lsi_753 = buffer.data(lsi + 753);
    const auto *lsi_754 = buffer.data(lsi + 754);
    const auto *lsi_755 = buffer.data(lsi + 755);
    const auto *lsi_756 = buffer.data(lsi + 756);
    const auto *lsi_757 = buffer.data(lsi + 757);
    const auto *lsi_758 = buffer.data(lsi + 758);
    const auto *lsi_759 = buffer.data(lsi + 759);
    const auto *lsi_760 = buffer.data(lsi + 760);
    const auto *lsi_761 = buffer.data(lsi + 761);
    const auto *lsi_762 = buffer.data(lsi + 762);
    const auto *lsi_763 = buffer.data(lsi + 763);
    const auto *lsi_764 = buffer.data(lsi + 764);
    const auto *lsi_765 = buffer.data(lsi + 765);
    const auto *lsi_766 = buffer.data(lsi + 766);
    const auto *lsi_767 = buffer.data(lsi + 767);
    const auto *lsi_768 = buffer.data(lsi + 768);
    const auto *lsi_769 = buffer.data(lsi + 769);
    const auto *lsi_770 = buffer.data(lsi + 770);
    const auto *lsi_776 = buffer.data(lsi + 776);
    const auto *lsi_777 = buffer.data(lsi + 777);
    const auto *lsi_778 = buffer.data(lsi + 778);
    const auto *lsi_779 = buffer.data(lsi + 779);
    const auto *lsi_780 = buffer.data(lsi + 780);
    const auto *lsi_781 = buffer.data(lsi + 781);
    const auto *lsi_782 = buffer.data(lsi + 782);
    const auto *lsi_783 = buffer.data(lsi + 783);
    const auto *lsi_784 = buffer.data(lsi + 784);
    const auto *lsi_785 = buffer.data(lsi + 785);
    const auto *lsi_786 = buffer.data(lsi + 786);
    const auto *lsi_787 = buffer.data(lsi + 787);
    const auto *lsi_789 = buffer.data(lsi + 789);
    const auto *lsi_790 = buffer.data(lsi + 790);
    const auto *lsi_791 = buffer.data(lsi + 791);
    const auto *lsi_793 = buffer.data(lsi + 793);
    const auto *lsi_794 = buffer.data(lsi + 794);
    const auto *lsi_795 = buffer.data(lsi + 795);
    const auto *lsi_796 = buffer.data(lsi + 796);
    const auto *lsi_798 = buffer.data(lsi + 798);
    const auto *lsi_799 = buffer.data(lsi + 799);
    const auto *lsi_805 = buffer.data(lsi + 805);
    const auto *lsi_807 = buffer.data(lsi + 807);
    const auto *lsi_808 = buffer.data(lsi + 808);
    const auto *lsi_809 = buffer.data(lsi + 809);
    const auto *lsi_810 = buffer.data(lsi + 810);
    const auto *lsi_811 = buffer.data(lsi + 811);
    const auto *lsi_812 = buffer.data(lsi + 812);
    const auto *lsi_814 = buffer.data(lsi + 814);

#pragma omp simd aligned(t_931, t_932, t_933, pc_y, ksi_556, ksi_557, ksi_558, lsh0_543, \
                         lsh0_544, lsh0_545, lsh1_543, lsh1_544, lsh1_545, lsi_724, lsi_725, \
                         lsi_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_14 * ksi_556[k]
                   + f_8 * lsh0_543[k]
                   - f_9 * lsh1_543[k]
                   + f_3 * pc_y[k] * lsi_724[k];

        t_932[k] = f_14 * ksi_557[k]
                   + f_6 * lsh0_544[k]
                   - f_7 * lsh1_544[k]
                   + f_3 * pc_y[k] * lsi_725[k];

        t_933[k] = f_14 * ksi_558[k]
                   + f_4 * lsh0_545[k]
                   - f_5 * lsh1_545[k]
                   + f_3 * pc_y[k] * lsi_726[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pa_y, pc_y, pc_z, ksk0_720, ksi_531, \
                         ksi_559, ksi_560, ksk1_720, lsh0_545, lsh1_545, lsi_727, \
                         lsi_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_14 * ksi_559[k]
                   + f_3 * pc_y[k] * lsi_727[k];

        t_935[k] = f_16 * ksi_531[k]
                   + f_1 * lsh0_545[k]
                   - f_2 * lsh1_545[k]
                   + f_3 * pc_z[k] * lsi_727[k];

        t_936[k] = pa_y[k] * ksk0_720[k]
                   - f_12 * pc_y[k] * ksk1_720[k];

        t_937[k] = f_13 * ksi_560[k]
                   + f_3 * pc_y[k] * lsi_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pa_y, pc_y, pc_z, ksk0_723, ksk0_725, \
                         ksi_532, ksi_561, ksi_562, ksk1_723, ksk1_725, lsi_728, \
                         lsi_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_17 * ksi_532[k]
                   + f_3 * pc_z[k] * lsi_728[k];

        t_939[k] = pa_y[k] * ksk0_723[k]
                   + f_14 * ksi_561[k]
                   - f_12 * pc_y[k] * ksk1_723[k];

        t_940[k] = f_13 * ksi_562[k]
                   + f_3 * pc_y[k] * lsi_730[k];

        t_941[k] = pa_y[k] * ksk0_725[k]
                   - f_12 * pc_y[k] * ksk1_725[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pa_y, pc_y, pc_z, ksk0_726, ksk0_729, \
                         ksi_535, ksi_563, ksi_565, ksk1_726, ksk1_729, lsi_731, \
                         lsi_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = pa_y[k] * ksk0_726[k]
                   + f_15 * ksi_563[k]
                   - f_12 * pc_y[k] * ksk1_726[k];

        t_943[k] = f_17 * ksi_535[k]
                   + f_3 * pc_z[k] * lsi_731[k];

        t_944[k] = f_13 * ksi_565[k]
                   + f_3 * pc_y[k] * lsi_733[k];

        t_945[k] = pa_y[k] * ksk0_729[k]
                   - f_12 * pc_y[k] * ksk1_729[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pa_y, pc_y, pc_z, ksk0_730, ksk0_732, ksi_538, \
                         ksi_566, ksi_568, ksk1_730, ksk1_732, \
                         lsi_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = pa_y[k] * ksk0_730[k]
                   + f_16 * ksi_566[k]
                   - f_12 * pc_y[k] * ksk1_730[k];

        t_947[k] = f_17 * ksi_538[k]
                   + f_3 * pc_z[k] * lsi_734[k];

        t_948[k] = pa_y[k] * ksk0_732[k]
                   + f_14 * ksi_568[k]
                   - f_12 * pc_y[k] * ksk1_732[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_y, pc_y, pc_z, ksk0_734, ksk0_735, \
                         ksi_542, ksi_569, ksi_570, ksk1_734, ksk1_735, lsi_737, \
                         lsi_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_13 * ksi_569[k]
                   + f_3 * pc_y[k] * lsi_737[k];

        t_950[k] = pa_y[k] * ksk0_734[k]
                   - f_12 * pc_y[k] * ksk1_734[k];

        t_951[k] = pa_y[k] * ksk0_735[k]
                   + f_17 * ksi_570[k]
                   - f_12 * pc_y[k] * ksk1_735[k];

        t_952[k] = f_17 * ksi_542[k]
                   + f_3 * pc_z[k] * lsi_738[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pa_y, pc_y, ksk0_737, ksk0_738, ksk0_740, \
                         ksi_572, ksi_573, ksi_574, ksk1_737, ksk1_738, ksk1_740, \
                         lsi_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = pa_y[k] * ksk0_737[k]
                   + f_15 * ksi_572[k]
                   - f_12 * pc_y[k] * ksk1_737[k];

        t_954[k] = pa_y[k] * ksk0_738[k]
                   + f_14 * ksi_573[k]
                   - f_12 * pc_y[k] * ksk1_738[k];

        t_955[k] = f_13 * ksi_574[k]
                   + f_3 * pc_y[k] * lsi_742[k];

        t_956[k] = pa_y[k] * ksk0_740[k]
                   - f_12 * pc_y[k] * ksk1_740[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, pc_x, ksi_749, ksi_750, ksi_751, \
                         ksi_752, ksi_753, lsi_749, lsi_750, lsi_751, lsi_752, \
                         lsi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_14 * ksi_749[k]
                   + f_3 * pc_x[k] * lsi_749[k];

        t_958[k] = f_14 * ksi_750[k]
                   + f_3 * pc_x[k] * lsi_750[k];

        t_959[k] = f_14 * ksi_751[k]
                   + f_3 * pc_x[k] * lsi_751[k];

        t_960[k] = f_14 * ksi_752[k]
                   + f_3 * pc_x[k] * lsi_752[k];

        t_961[k] = f_14 * ksi_753[k]
                   + f_3 * pc_x[k] * lsi_753[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pc_x, pc_y, pc_z, ksi_553, ksi_581, \
                         ksi_754, ksi_755, lsh0_561, lsh1_561, lsi_749, lsi_754, \
                         lsi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_14 * ksi_754[k]
                   + f_3 * pc_x[k] * lsi_754[k];

        t_963[k] = f_14 * ksi_755[k]
                   + f_3 * pc_x[k] * lsi_755[k];

        t_964[k] = f_13 * ksi_581[k]
                   + f_1 * lsh0_561[k]
                   - f_2 * lsh1_561[k]
                   + f_3 * pc_y[k] * lsi_749[k];

        t_965[k] = f_17 * ksi_553[k]
                   + f_3 * pc_z[k] * lsi_749[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_y, ksi_583, ksi_584, ksi_585, lsh0_563, \
                         lsh0_564, lsh0_565, lsh1_563, lsh1_564, lsh1_565, lsi_751, lsi_752, \
                         lsi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_13 * ksi_583[k]
                   + f_10 * lsh0_563[k]
                   - f_11 * lsh1_563[k]
                   + f_3 * pc_y[k] * lsi_751[k];

        t_967[k] = f_13 * ksi_584[k]
                   + f_8 * lsh0_564[k]
                   - f_9 * lsh1_564[k]
                   + f_3 * pc_y[k] * lsi_752[k];

        t_968[k] = f_13 * ksi_585[k]
                   + f_6 * lsh0_565[k]
                   - f_7 * lsh1_565[k]
                   + f_3 * pc_y[k] * lsi_753[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pa_y, pc_y, ksk0_755, ksi_586, ksi_587, \
                         ksk1_755, lsh0_566, lsh1_566, lsi_754, \
                         lsi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_13 * ksi_586[k]
                   + f_4 * lsh0_566[k]
                   - f_5 * lsh1_566[k]
                   + f_3 * pc_y[k] * lsi_754[k];

        t_970[k] = f_13 * ksi_587[k]
                   + f_3 * pc_y[k] * lsi_755[k];

        t_971[k] = pa_y[k] * ksk0_755[k]
                   - f_12 * pc_y[k] * ksk1_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, t_976, pc_x, pc_y, pc_z, ksi_560, \
                         ksi_756, lsh0_567, lsh1_567, lsi_756, lsi_757, \
                         lsi_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_14 * ksi_756[k]
                   + f_1 * lsh0_567[k]
                   - f_2 * lsh1_567[k]
                   + f_3 * pc_x[k] * lsi_756[k];

        t_973[k] = f_3 * pc_y[k] * lsi_756[k];

        t_974[k] = f_21 * ksi_560[k]
                   + f_3 * pc_z[k] * lsi_756[k];

        t_975[k] = f_4 * lsh0_567[k]
                   - f_5 * lsh1_567[k]
                   + f_3 * pc_y[k] * lsi_757[k];

        t_976[k] = f_3 * pc_y[k] * lsi_758[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, pc_x, pc_y, ksi_761, lsh0_568, lsh0_569, \
                         lsh0_572, lsh1_568, lsh1_569, lsh1_572, lsi_759, lsi_760, \
                         lsi_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_14 * ksi_761[k]
                   + f_10 * lsh0_572[k]
                   - f_11 * lsh1_572[k]
                   + f_3 * pc_x[k] * lsi_761[k];

        t_978[k] = f_6 * lsh0_568[k]
                   - f_7 * lsh1_568[k]
                   + f_3 * pc_y[k] * lsi_759[k];

        t_979[k] = f_4 * lsh0_569[k]
                   - f_5 * lsh1_569[k]
                   + f_3 * pc_y[k] * lsi_760[k];

        t_980[k] = f_3 * pc_y[k] * lsi_761[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_x, pc_y, ksi_765, lsh0_570, lsh0_571, \
                         lsh0_576, lsh1_570, lsh1_571, lsh1_576, lsi_762, lsi_763, \
                         lsi_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_14 * ksi_765[k]
                   + f_8 * lsh0_576[k]
                   - f_9 * lsh1_576[k]
                   + f_3 * pc_x[k] * lsi_765[k];

        t_982[k] = f_8 * lsh0_570[k]
                   - f_9 * lsh1_570[k]
                   + f_3 * pc_y[k] * lsi_762[k];

        t_983[k] = f_6 * lsh0_571[k]
                   - f_7 * lsh1_571[k]
                   + f_3 * pc_y[k] * lsi_763[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pc_x, pc_y, ksi_770, lsh0_572, lsh0_581, \
                         lsh1_572, lsh1_581, lsi_764, lsi_765, \
                         lsi_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_4 * lsh0_572[k]
                   - f_5 * lsh1_572[k]
                   + f_3 * pc_y[k] * lsi_764[k];

        t_985[k] = f_3 * pc_y[k] * lsi_765[k];

        t_986[k] = f_14 * ksi_770[k]
                   + f_6 * lsh0_581[k]
                   - f_7 * lsh1_581[k]
                   + f_3 * pc_x[k] * lsi_770[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, pc_y, lsh0_573, lsh0_574, lsh0_575, lsh1_573, \
                         lsh1_574, lsh1_575, lsi_766, lsi_767, \
                         lsi_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_10 * lsh0_573[k]
                   - f_11 * lsh1_573[k]
                   + f_3 * pc_y[k] * lsi_766[k];

        t_988[k] = f_8 * lsh0_574[k]
                   - f_9 * lsh1_574[k]
                   + f_3 * pc_y[k] * lsi_767[k];

        t_989[k] = f_6 * lsh0_575[k]
                   - f_7 * lsh1_575[k]
                   + f_3 * pc_y[k] * lsi_768[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, pc_x, pc_y, ksi_776, ksi_777, lsh0_576, \
                         lsh0_587, lsh1_576, lsh1_587, lsi_769, lsi_770, lsi_776, \
                         lsi_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = f_4 * lsh0_576[k]
                   - f_5 * lsh1_576[k]
                   + f_3 * pc_y[k] * lsi_769[k];

        t_991[k] = f_3 * pc_y[k] * lsi_770[k];

        t_992[k] = f_14 * ksi_776[k]
                   + f_4 * lsh0_587[k]
                   - f_5 * lsh1_587[k]
                   + f_3 * pc_x[k] * lsi_776[k];

        t_993[k] = f_14 * ksi_777[k]
                   + f_3 * pc_x[k] * lsi_777[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, t_997, t_998, pc_x, pc_y, ksi_778, ksi_779, \
                         ksi_780, ksi_781, lsi_776, lsi_778, lsi_779, lsi_780, \
                         lsi_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_14 * ksi_778[k]
                   + f_3 * pc_x[k] * lsi_778[k];

        t_995[k] = f_14 * ksi_779[k]
                   + f_3 * pc_x[k] * lsi_779[k];

        t_996[k] = f_14 * ksi_780[k]
                   + f_3 * pc_x[k] * lsi_780[k];

        t_997[k] = f_14 * ksi_781[k]
                   + f_3 * pc_x[k] * lsi_781[k];

        t_998[k] = f_3 * pc_y[k] * lsi_776[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, pc_x, pc_y, ksi_783, lsh0_582, lsh0_583, \
                         lsh1_582, lsh1_583, lsi_777, lsi_778, \
                         lsi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_14 * ksi_783[k]
                   + f_3 * pc_x[k] * lsi_783[k];

        t_1000[k] = f_1 * lsh0_582[k]
                    - f_2 * lsh1_582[k]
                    + f_3 * pc_y[k] * lsi_777[k];

        t_1001[k] = f_19 * lsh0_583[k]
                    - f_20 * lsh1_583[k]
                    + f_3 * pc_y[k] * lsi_778[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, lsh0_584, lsh0_585, lsh0_586, lsh1_584, \
                         lsh1_585, lsh1_586, lsi_779, lsi_780, \
                         lsi_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_10 * lsh0_584[k]
                    - f_11 * lsh1_584[k]
                    + f_3 * pc_y[k] * lsi_779[k];

        t_1003[k] = f_8 * lsh0_585[k]
                    - f_9 * lsh1_585[k]
                    + f_3 * pc_y[k] * lsi_780[k];

        t_1004[k] = f_6 * lsh0_586[k]
                    - f_7 * lsh1_586[k]
                    + f_3 * pc_y[k] * lsi_781[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, t_1008, pa_x, pc_x, pc_y, pc_z, ksk0_1008, \
                         ksi_587, ksi_784, ksk1_1008, lsh0_587, lsh1_587, lsi_782, \
                         lsi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_4 * lsh0_587[k]
                    - f_5 * lsh1_587[k]
                    + f_3 * pc_y[k] * lsi_782[k];

        t_1006[k] = f_3 * pc_y[k] * lsi_783[k];

        t_1007[k] = f_21 * ksi_587[k]
                    + f_1 * lsh0_587[k]
                    - f_2 * lsh1_587[k]
                    + f_3 * pc_z[k] * lsi_783[k];

        t_1008[k] = pa_x[k] * ksk0_1008[k]
                    + f_18 * ksi_784[k]
                    - f_12 * pc_x[k] * ksk1_1008[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, pa_x, pc_x, pc_y, pc_z, ksk0_1011, \
                         ksi_588, ksi_787, ksk1_1011, lsi_784, \
                         lsi_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_18 * ksi_588[k]
                    + f_3 * pc_y[k] * lsi_784[k];

        t_1010[k] = f_3 * pc_z[k] * lsi_784[k];

        t_1011[k] = pa_x[k] * ksk0_1011[k]
                    + f_17 * ksi_787[k]
                    - f_12 * pc_x[k] * ksk1_1011[k];

        t_1012[k] = f_3 * pc_z[k] * lsi_785[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pa_x, pc_x, pc_z, ksk0_1014, ksi_790, \
                         ksk1_1014, lsh0_588, lsh1_588, lsi_786, \
                         lsi_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_4 * lsh0_588[k]
                    - f_5 * lsh1_588[k]
                    + f_3 * pc_z[k] * lsi_786[k];

        t_1014[k] = pa_x[k] * ksk0_1014[k]
                    + f_16 * ksi_790[k]
                    - f_12 * pc_x[k] * ksk1_1014[k];

        t_1015[k] = f_3 * pc_z[k] * lsi_787[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pa_x, pc_x, pc_y, pc_z, ksk0_1018, \
                         ksi_593, ksi_794, ksk1_1018, lsh0_590, lsh1_590, lsi_789, \
                         lsi_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_18 * ksi_593[k]
                    + f_3 * pc_y[k] * lsi_789[k];

        t_1017[k] = f_6 * lsh0_590[k]
                    - f_7 * lsh1_590[k]
                    + f_3 * pc_z[k] * lsi_789[k];

        t_1018[k] = pa_x[k] * ksk0_1018[k]
                    + f_15 * ksi_794[k]
                    - f_12 * pc_x[k] * ksk1_1018[k];

        t_1019[k] = f_3 * pc_z[k] * lsi_790[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, pc_y, pc_z, ksi_597, lsh0_591, lsh0_593, \
                         lsh1_591, lsh1_593, lsi_791, lsi_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_4 * lsh0_591[k]
                    - f_5 * lsh1_591[k]
                    + f_3 * pc_z[k] * lsi_791[k];

        t_1021[k] = f_18 * ksi_597[k]
                    + f_3 * pc_y[k] * lsi_793[k];

        t_1022[k] = f_8 * lsh0_593[k]
                    - f_9 * lsh1_593[k]
                    + f_3 * pc_z[k] * lsi_793[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, pa_x, pc_x, pc_z, ksk0_1023, ksi_799, \
                         ksk1_1023, lsh0_594, lsh1_594, lsi_794, \
                         lsi_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = pa_x[k] * ksk0_1023[k]
                    + f_14 * ksi_799[k]
                    - f_12 * pc_x[k] * ksk1_1023[k];

        t_1024[k] = f_3 * pc_z[k] * lsi_794[k];

        t_1025[k] = f_4 * lsh0_594[k]
                    - f_5 * lsh1_594[k]
                    + f_3 * pc_z[k] * lsi_795[k];
    }

#pragma omp simd aligned(t_1026, t_1027, t_1028, t_1029, pc_x, pc_y, pc_z, ksi_602, ksi_805, \
                         lsh0_595, lsh0_597, lsh1_595, lsh1_597, lsi_796, lsi_798, \
                         lsi_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1026[k] = f_6 * lsh0_595[k]
                    - f_7 * lsh1_595[k]
                    + f_3 * pc_z[k] * lsi_796[k];

        t_1027[k] = f_18 * ksi_602[k]
                    + f_3 * pc_y[k] * lsi_798[k];

        t_1028[k] = f_10 * lsh0_597[k]
                    - f_11 * lsh1_597[k]
                    + f_3 * pc_z[k] * lsi_798[k];

        t_1029[k] = f_13 * ksi_805[k]
                    + f_3 * pc_x[k] * lsi_805[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, t_1034, pc_x, pc_z, ksi_807, ksi_808, \
                         ksi_809, ksi_810, lsi_799, lsi_807, lsi_808, lsi_809, \
                         lsi_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = f_3 * pc_z[k] * lsi_799[k];

        t_1031[k] = f_13 * ksi_807[k]
                    + f_3 * pc_x[k] * lsi_807[k];

        t_1032[k] = f_13 * ksi_808[k]
                    + f_3 * pc_x[k] * lsi_808[k];

        t_1033[k] = f_13 * ksi_809[k]
                    + f_3 * pc_x[k] * lsi_809[k];

        t_1034[k] = f_13 * ksi_810[k]
                    + f_3 * pc_x[k] * lsi_810[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, pa_x, pc_x, pc_z, ksk0_1036, \
                         ksk0_1038, ksi_811, ksk1_1036, ksk1_1038, lsi_805, \
                         lsi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = f_13 * ksi_811[k]
                    + f_3 * pc_x[k] * lsi_811[k];

        t_1036[k] = pa_x[k] * ksk0_1036[k]
                    - f_12 * pc_x[k] * ksk1_1036[k];

        t_1037[k] = f_3 * pc_z[k] * lsi_805[k];

        t_1038[k] = pa_x[k] * ksk0_1038[k]
                    - f_12 * pc_x[k] * ksk1_1038[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, t_1042, pa_x, pc_x, pc_y, ksk0_1039, \
                         ksk0_1040, ksk0_1041, ksi_615, ksk1_1039, ksk1_1040, ksk1_1041, \
                         lsi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = pa_x[k] * ksk0_1039[k]
                    - f_12 * pc_x[k] * ksk1_1039[k];

        t_1040[k] = pa_x[k] * ksk0_1040[k]
                    - f_12 * pc_x[k] * ksk1_1040[k];

        t_1041[k] = pa_x[k] * ksk0_1041[k]
                    - f_12 * pc_x[k] * ksk1_1041[k];

        t_1042[k] = f_18 * ksi_615[k]
                    + f_3 * pc_y[k] * lsi_811[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, t_1046, pa_x, pa_z, pc_x, pc_y, pc_z, \
                         ksk0_756, ksk0_1043, ksi_588, ksi_616, ksk1_756, ksk1_1043, \
                         lsi_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = pa_x[k] * ksk0_1043[k]
                    - f_12 * pc_x[k] * ksk1_1043[k];

        t_1044[k] = pa_z[k] * ksk0_756[k]
                    - f_12 * pc_z[k] * ksk1_756[k];

        t_1045[k] = f_21 * ksi_616[k]
                    + f_3 * pc_y[k] * lsi_812[k];

        t_1046[k] = f_13 * ksi_588[k]
                    + f_3 * pc_z[k] * lsi_812[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pa_x, pa_z, pc_x, pc_y, pc_z, ksk0_759, \
                         ksk0_1049, ksi_618, ksi_817, ksk1_759, ksk1_1049, \
                         lsi_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = pa_z[k] * ksk0_759[k]
                    - f_12 * pc_z[k] * ksk1_759[k];

        t_1048[k] = f_21 * ksi_618[k]
                    + f_3 * pc_y[k] * lsi_814[k];

        t_1049[k] = pa_x[k] * ksk0_1049[k]
                    + f_17 * ksi_817[k]
                    - f_12 * pc_x[k] * ksk1_1049[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksk0,
                                                          const size_t ksi, const size_t ksk1,
                                                          const size_t lsi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_762 = buffer.data(ksk0 + 762);
    const auto *ksk0_766 = buffer.data(ksk0 + 766);
    const auto *ksk0_771 = buffer.data(ksk0 + 771);
    const auto *ksk0_1053 = buffer.data(ksk0 + 1053);
    const auto *ksk0_1056 = buffer.data(ksk0 + 1056);
    const auto *ksk0_1058 = buffer.data(ksk0 + 1058);
    const auto *ksk0_1061 = buffer.data(ksk0 + 1061);
    const auto *ksk0_1062 = buffer.data(ksk0 + 1062);
    const auto *ksk0_1064 = buffer.data(ksk0 + 1064);
    const auto *ksk0_1072 = buffer.data(ksk0 + 1072);
    const auto *ksk0_1074 = buffer.data(ksk0 + 1074);
    const auto *ksk0_1075 = buffer.data(ksk0 + 1075);
    const auto *ksk0_1076 = buffer.data(ksk0 + 1076);
    const auto *ksk0_1077 = buffer.data(ksk0 + 1077);
    const auto *ksk0_1079 = buffer.data(ksk0 + 1079);
    const auto *ksk0_1080 = buffer.data(ksk0 + 1080);
    const auto *ksk0_1083 = buffer.data(ksk0 + 1083);
    const auto *ksk0_1085 = buffer.data(ksk0 + 1085);
    const auto *ksk0_1086 = buffer.data(ksk0 + 1086);
    const auto *ksk0_1089 = buffer.data(ksk0 + 1089);
    const auto *ksk0_1090 = buffer.data(ksk0 + 1090);
    const auto *ksk0_1092 = buffer.data(ksk0 + 1092);
    const auto *ksk0_1094 = buffer.data(ksk0 + 1094);
    const auto *ksk0_1095 = buffer.data(ksk0 + 1095);
    const auto *ksk0_1097 = buffer.data(ksk0 + 1097);
    const auto *ksk0_1098 = buffer.data(ksk0 + 1098);
    const auto *ksk0_1100 = buffer.data(ksk0 + 1100);
    const auto *ksk0_1108 = buffer.data(ksk0 + 1108);
    const auto *ksk0_1110 = buffer.data(ksk0 + 1110);
    const auto *ksk0_1111 = buffer.data(ksk0 + 1111);
    const auto *ksk0_1112 = buffer.data(ksk0 + 1112);
    const auto *ksk0_1113 = buffer.data(ksk0 + 1113);
    const auto *ksk0_1115 = buffer.data(ksk0 + 1115);
    const auto *ksk0_1116 = buffer.data(ksk0 + 1116);
    const auto *ksk0_1119 = buffer.data(ksk0 + 1119);
    const auto *ksk0_1121 = buffer.data(ksk0 + 1121);
    const auto *ksk0_1122 = buffer.data(ksk0 + 1122);
    const auto *ksk0_1125 = buffer.data(ksk0 + 1125);
    const auto *ksk0_1126 = buffer.data(ksk0 + 1126);
    const auto *ksk0_1128 = buffer.data(ksk0 + 1128);
    const auto *ksk0_1130 = buffer.data(ksk0 + 1130);
    const auto *ksk0_1131 = buffer.data(ksk0 + 1131);
    const auto *ksk0_1133 = buffer.data(ksk0 + 1133);
    const auto *ksk0_1134 = buffer.data(ksk0 + 1134);
    const auto *ksk0_1136 = buffer.data(ksk0 + 1136);
    const auto *ksk0_1144 = buffer.data(ksk0 + 1144);
    const auto *ksk0_1146 = buffer.data(ksk0 + 1146);
    const auto *ksk0_1147 = buffer.data(ksk0 + 1147);
    const auto *ksk0_1148 = buffer.data(ksk0 + 1148);
    const auto *ksk0_1149 = buffer.data(ksk0 + 1149);
    const auto *ksk0_1151 = buffer.data(ksk0 + 1151);
    const auto *ksk0_1152 = buffer.data(ksk0 + 1152);
    const auto *ksk0_1155 = buffer.data(ksk0 + 1155);
    const auto *ksk0_1157 = buffer.data(ksk0 + 1157);
    const auto *ksk0_1158 = buffer.data(ksk0 + 1158);
    const auto *ksk0_1161 = buffer.data(ksk0 + 1161);
    const auto *ksk0_1162 = buffer.data(ksk0 + 1162);
    const auto *ksk0_1164 = buffer.data(ksk0 + 1164);
    const auto *ksk0_1166 = buffer.data(ksk0 + 1166);
    const auto *ksk0_1167 = buffer.data(ksk0 + 1167);

    const auto *ksi_591 = buffer.data(ksi + 591);
    const auto *ksi_594 = buffer.data(ksi + 594);
    const auto *ksi_598 = buffer.data(ksi + 598);
    const auto *ksi_609 = buffer.data(ksi + 609);
    const auto *ksi_616 = buffer.data(ksi + 616);
    const auto *ksi_619 = buffer.data(ksi + 619);
    const auto *ksi_621 = buffer.data(ksi + 621);
    const auto *ksi_622 = buffer.data(ksi + 622);
    const auto *ksi_625 = buffer.data(ksi + 625);
    const auto *ksi_626 = buffer.data(ksi + 626);
    const auto *ksi_630 = buffer.data(ksi + 630);
    const auto *ksi_637 = buffer.data(ksi + 637);
    const auto *ksi_643 = buffer.data(ksi + 643);
    const auto *ksi_644 = buffer.data(ksi + 644);
    const auto *ksi_646 = buffer.data(ksi + 646);
    const auto *ksi_647 = buffer.data(ksi + 647);
    const auto *ksi_649 = buffer.data(ksi + 649);
    const auto *ksi_650 = buffer.data(ksi + 650);
    const auto *ksi_653 = buffer.data(ksi + 653);
    const auto *ksi_654 = buffer.data(ksi + 654);
    const auto *ksi_658 = buffer.data(ksi + 658);
    const auto *ksi_665 = buffer.data(ksi + 665);
    const auto *ksi_671 = buffer.data(ksi + 671);
    const auto *ksi_672 = buffer.data(ksi + 672);
    const auto *ksi_674 = buffer.data(ksi + 674);
    const auto *ksi_675 = buffer.data(ksi + 675);
    const auto *ksi_677 = buffer.data(ksi + 677);
    const auto *ksi_678 = buffer.data(ksi + 678);
    const auto *ksi_681 = buffer.data(ksi + 681);
    const auto *ksi_682 = buffer.data(ksi + 682);
    const auto *ksi_686 = buffer.data(ksi + 686);
    const auto *ksi_699 = buffer.data(ksi + 699);
    const auto *ksi_700 = buffer.data(ksi + 700);
    const auto *ksi_702 = buffer.data(ksi + 702);
    const auto *ksi_705 = buffer.data(ksi + 705);
    const auto *ksi_709 = buffer.data(ksi + 709);
    const auto *ksi_821 = buffer.data(ksi + 821);
    const auto *ksi_824 = buffer.data(ksi + 824);
    const auto *ksi_826 = buffer.data(ksi + 826);
    const auto *ksi_829 = buffer.data(ksi + 829);
    const auto *ksi_830 = buffer.data(ksi + 830);
    const auto *ksi_832 = buffer.data(ksi + 832);
    const auto *ksi_833 = buffer.data(ksi + 833);
    const auto *ksi_834 = buffer.data(ksi + 834);
    const auto *ksi_835 = buffer.data(ksi + 835);
    const auto *ksi_836 = buffer.data(ksi + 836);
    const auto *ksi_837 = buffer.data(ksi + 837);
    const auto *ksi_838 = buffer.data(ksi + 838);
    const auto *ksi_839 = buffer.data(ksi + 839);
    const auto *ksi_840 = buffer.data(ksi + 840);
    const auto *ksi_843 = buffer.data(ksi + 843);
    const auto *ksi_845 = buffer.data(ksi + 845);
    const auto *ksi_846 = buffer.data(ksi + 846);
    const auto *ksi_849 = buffer.data(ksi + 849);
    const auto *ksi_850 = buffer.data(ksi + 850);
    const auto *ksi_852 = buffer.data(ksi + 852);
    const auto *ksi_854 = buffer.data(ksi + 854);
    const auto *ksi_855 = buffer.data(ksi + 855);
    const auto *ksi_857 = buffer.data(ksi + 857);
    const auto *ksi_858 = buffer.data(ksi + 858);
    const auto *ksi_860 = buffer.data(ksi + 860);
    const auto *ksi_861 = buffer.data(ksi + 861);
    const auto *ksi_862 = buffer.data(ksi + 862);
    const auto *ksi_863 = buffer.data(ksi + 863);
    const auto *ksi_864 = buffer.data(ksi + 864);
    const auto *ksi_865 = buffer.data(ksi + 865);
    const auto *ksi_866 = buffer.data(ksi + 866);
    const auto *ksi_867 = buffer.data(ksi + 867);
    const auto *ksi_868 = buffer.data(ksi + 868);
    const auto *ksi_871 = buffer.data(ksi + 871);
    const auto *ksi_873 = buffer.data(ksi + 873);
    const auto *ksi_874 = buffer.data(ksi + 874);
    const auto *ksi_877 = buffer.data(ksi + 877);
    const auto *ksi_878 = buffer.data(ksi + 878);
    const auto *ksi_880 = buffer.data(ksi + 880);
    const auto *ksi_882 = buffer.data(ksi + 882);
    const auto *ksi_883 = buffer.data(ksi + 883);
    const auto *ksi_885 = buffer.data(ksi + 885);
    const auto *ksi_886 = buffer.data(ksi + 886);
    const auto *ksi_888 = buffer.data(ksi + 888);
    const auto *ksi_889 = buffer.data(ksi + 889);
    const auto *ksi_890 = buffer.data(ksi + 890);
    const auto *ksi_891 = buffer.data(ksi + 891);
    const auto *ksi_892 = buffer.data(ksi + 892);
    const auto *ksi_893 = buffer.data(ksi + 893);
    const auto *ksi_894 = buffer.data(ksi + 894);
    const auto *ksi_895 = buffer.data(ksi + 895);
    const auto *ksi_896 = buffer.data(ksi + 896);
    const auto *ksi_899 = buffer.data(ksi + 899);
    const auto *ksi_901 = buffer.data(ksi + 901);
    const auto *ksi_902 = buffer.data(ksi + 902);
    const auto *ksi_905 = buffer.data(ksi + 905);
    const auto *ksi_906 = buffer.data(ksi + 906);
    const auto *ksi_908 = buffer.data(ksi + 908);
    const auto *ksi_910 = buffer.data(ksi + 910);
    const auto *ksi_911 = buffer.data(ksi + 911);

    const auto *ksk1_762 = buffer.data(ksk1 + 762);
    const auto *ksk1_766 = buffer.data(ksk1 + 766);
    const auto *ksk1_771 = buffer.data(ksk1 + 771);
    const auto *ksk1_1053 = buffer.data(ksk1 + 1053);
    const auto *ksk1_1056 = buffer.data(ksk1 + 1056);
    const auto *ksk1_1058 = buffer.data(ksk1 + 1058);
    const auto *ksk1_1061 = buffer.data(ksk1 + 1061);
    const auto *ksk1_1062 = buffer.data(ksk1 + 1062);
    const auto *ksk1_1064 = buffer.data(ksk1 + 1064);
    const auto *ksk1_1072 = buffer.data(ksk1 + 1072);
    const auto *ksk1_1074 = buffer.data(ksk1 + 1074);
    const auto *ksk1_1075 = buffer.data(ksk1 + 1075);
    const auto *ksk1_1076 = buffer.data(ksk1 + 1076);
    const auto *ksk1_1077 = buffer.data(ksk1 + 1077);
    const auto *ksk1_1079 = buffer.data(ksk1 + 1079);
    const auto *ksk1_1080 = buffer.data(ksk1 + 1080);
    const auto *ksk1_1083 = buffer.data(ksk1 + 1083);
    const auto *ksk1_1085 = buffer.data(ksk1 + 1085);
    const auto *ksk1_1086 = buffer.data(ksk1 + 1086);
    const auto *ksk1_1089 = buffer.data(ksk1 + 1089);
    const auto *ksk1_1090 = buffer.data(ksk1 + 1090);
    const auto *ksk1_1092 = buffer.data(ksk1 + 1092);
    const auto *ksk1_1094 = buffer.data(ksk1 + 1094);
    const auto *ksk1_1095 = buffer.data(ksk1 + 1095);
    const auto *ksk1_1097 = buffer.data(ksk1 + 1097);
    const auto *ksk1_1098 = buffer.data(ksk1 + 1098);
    const auto *ksk1_1100 = buffer.data(ksk1 + 1100);
    const auto *ksk1_1108 = buffer.data(ksk1 + 1108);
    const auto *ksk1_1110 = buffer.data(ksk1 + 1110);
    const auto *ksk1_1111 = buffer.data(ksk1 + 1111);
    const auto *ksk1_1112 = buffer.data(ksk1 + 1112);
    const auto *ksk1_1113 = buffer.data(ksk1 + 1113);
    const auto *ksk1_1115 = buffer.data(ksk1 + 1115);
    const auto *ksk1_1116 = buffer.data(ksk1 + 1116);
    const auto *ksk1_1119 = buffer.data(ksk1 + 1119);
    const auto *ksk1_1121 = buffer.data(ksk1 + 1121);
    const auto *ksk1_1122 = buffer.data(ksk1 + 1122);
    const auto *ksk1_1125 = buffer.data(ksk1 + 1125);
    const auto *ksk1_1126 = buffer.data(ksk1 + 1126);
    const auto *ksk1_1128 = buffer.data(ksk1 + 1128);
    const auto *ksk1_1130 = buffer.data(ksk1 + 1130);
    const auto *ksk1_1131 = buffer.data(ksk1 + 1131);
    const auto *ksk1_1133 = buffer.data(ksk1 + 1133);
    const auto *ksk1_1134 = buffer.data(ksk1 + 1134);
    const auto *ksk1_1136 = buffer.data(ksk1 + 1136);
    const auto *ksk1_1144 = buffer.data(ksk1 + 1144);
    const auto *ksk1_1146 = buffer.data(ksk1 + 1146);
    const auto *ksk1_1147 = buffer.data(ksk1 + 1147);
    const auto *ksk1_1148 = buffer.data(ksk1 + 1148);
    const auto *ksk1_1149 = buffer.data(ksk1 + 1149);
    const auto *ksk1_1151 = buffer.data(ksk1 + 1151);
    const auto *ksk1_1152 = buffer.data(ksk1 + 1152);
    const auto *ksk1_1155 = buffer.data(ksk1 + 1155);
    const auto *ksk1_1157 = buffer.data(ksk1 + 1157);
    const auto *ksk1_1158 = buffer.data(ksk1 + 1158);
    const auto *ksk1_1161 = buffer.data(ksk1 + 1161);
    const auto *ksk1_1162 = buffer.data(ksk1 + 1162);
    const auto *ksk1_1164 = buffer.data(ksk1 + 1164);
    const auto *ksk1_1166 = buffer.data(ksk1 + 1166);
    const auto *ksk1_1167 = buffer.data(ksk1 + 1167);

    const auto *lsi_815 = buffer.data(lsi + 815);
    const auto *lsi_817 = buffer.data(lsi + 817);
    const auto *lsi_818 = buffer.data(lsi + 818);
    const auto *lsi_821 = buffer.data(lsi + 821);
    const auto *lsi_822 = buffer.data(lsi + 822);
    const auto *lsi_826 = buffer.data(lsi + 826);
    const auto *lsi_833 = buffer.data(lsi + 833);
    const auto *lsi_834 = buffer.data(lsi + 834);
    const auto *lsi_835 = buffer.data(lsi + 835);
    const auto *lsi_836 = buffer.data(lsi + 836);
    const auto *lsi_837 = buffer.data(lsi + 837);
    const auto *lsi_838 = buffer.data(lsi + 838);
    const auto *lsi_839 = buffer.data(lsi + 839);
    const auto *lsi_840 = buffer.data(lsi + 840);
    const auto *lsi_842 = buffer.data(lsi + 842);
    const auto *lsi_843 = buffer.data(lsi + 843);
    const auto *lsi_845 = buffer.data(lsi + 845);
    const auto *lsi_846 = buffer.data(lsi + 846);
    const auto *lsi_849 = buffer.data(lsi + 849);
    const auto *lsi_850 = buffer.data(lsi + 850);
    const auto *lsi_854 = buffer.data(lsi + 854);
    const auto *lsi_861 = buffer.data(lsi + 861);
    const auto *lsi_862 = buffer.data(lsi + 862);
    const auto *lsi_863 = buffer.data(lsi + 863);
    const auto *lsi_864 = buffer.data(lsi + 864);
    const auto *lsi_865 = buffer.data(lsi + 865);
    const auto *lsi_866 = buffer.data(lsi + 866);
    const auto *lsi_867 = buffer.data(lsi + 867);
    const auto *lsi_868 = buffer.data(lsi + 868);
    const auto *lsi_870 = buffer.data(lsi + 870);
    const auto *lsi_871 = buffer.data(lsi + 871);
    const auto *lsi_873 = buffer.data(lsi + 873);
    const auto *lsi_874 = buffer.data(lsi + 874);
    const auto *lsi_877 = buffer.data(lsi + 877);
    const auto *lsi_878 = buffer.data(lsi + 878);
    const auto *lsi_882 = buffer.data(lsi + 882);
    const auto *lsi_889 = buffer.data(lsi + 889);
    const auto *lsi_890 = buffer.data(lsi + 890);
    const auto *lsi_891 = buffer.data(lsi + 891);
    const auto *lsi_892 = buffer.data(lsi + 892);
    const auto *lsi_893 = buffer.data(lsi + 893);
    const auto *lsi_894 = buffer.data(lsi + 894);
    const auto *lsi_895 = buffer.data(lsi + 895);
    const auto *lsi_896 = buffer.data(lsi + 896);
    const auto *lsi_898 = buffer.data(lsi + 898);
    const auto *lsi_899 = buffer.data(lsi + 899);
    const auto *lsi_901 = buffer.data(lsi + 901);
    const auto *lsi_902 = buffer.data(lsi + 902);
    const auto *lsi_905 = buffer.data(lsi + 905);
    const auto *lsi_906 = buffer.data(lsi + 906);

#pragma omp simd aligned(t_1050, t_1051, t_1052, pa_z, pc_y, pc_z, ksk0_762, ksi_591, ksi_621, \
                         ksk1_762, lsi_815, lsi_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = pa_z[k] * ksk0_762[k]
                    - f_12 * pc_z[k] * ksk1_762[k];

        t_1051[k] = f_13 * ksi_591[k]
                    + f_3 * pc_z[k] * lsi_815[k];

        t_1052[k] = f_21 * ksi_621[k]
                    + f_3 * pc_y[k] * lsi_817[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, pa_x, pa_z, pc_x, pc_z, ksk0_766, ksk0_1053, \
                         ksi_594, ksi_821, ksk1_766, ksk1_1053, \
                         lsi_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = pa_x[k] * ksk0_1053[k]
                    + f_16 * ksi_821[k]
                    - f_12 * pc_x[k] * ksk1_1053[k];

        t_1054[k] = pa_z[k] * ksk0_766[k]
                    - f_12 * pc_z[k] * ksk1_766[k];

        t_1055[k] = f_13 * ksi_594[k]
                    + f_3 * pc_z[k] * lsi_818[k];
    }

#pragma omp simd aligned(t_1056, t_1057, t_1058, pa_x, pc_x, pc_y, ksk0_1056, ksk0_1058, \
                         ksi_625, ksi_824, ksi_826, ksk1_1056, ksk1_1058, \
                         lsi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1056[k] = pa_x[k] * ksk0_1056[k]
                    + f_15 * ksi_824[k]
                    - f_12 * pc_x[k] * ksk1_1056[k];

        t_1057[k] = f_21 * ksi_625[k]
                    + f_3 * pc_y[k] * lsi_821[k];

        t_1058[k] = pa_x[k] * ksk0_1058[k]
                    + f_15 * ksi_826[k]
                    - f_12 * pc_x[k] * ksk1_1058[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, pa_x, pa_z, pc_x, pc_z, ksk0_771, ksk0_1061, \
                         ksi_598, ksi_829, ksk1_771, ksk1_1061, \
                         lsi_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = pa_z[k] * ksk0_771[k]
                    - f_12 * pc_z[k] * ksk1_771[k];

        t_1060[k] = f_13 * ksi_598[k]
                    + f_3 * pc_z[k] * lsi_822[k];

        t_1061[k] = pa_x[k] * ksk0_1061[k]
                    + f_14 * ksi_829[k]
                    - f_12 * pc_x[k] * ksk1_1061[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, pa_x, pc_x, pc_y, ksk0_1062, ksk0_1064, \
                         ksi_630, ksi_830, ksi_832, ksk1_1062, ksk1_1064, \
                         lsi_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = pa_x[k] * ksk0_1062[k]
                    + f_14 * ksi_830[k]
                    - f_12 * pc_x[k] * ksk1_1062[k];

        t_1063[k] = f_21 * ksi_630[k]
                    + f_3 * pc_y[k] * lsi_826[k];

        t_1064[k] = pa_x[k] * ksk0_1064[k]
                    + f_14 * ksi_832[k]
                    - f_12 * pc_x[k] * ksk1_1064[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, pc_x, ksi_833, ksi_834, \
                         ksi_835, ksi_836, ksi_837, lsi_833, lsi_834, lsi_835, lsi_836, \
                         lsi_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_13 * ksi_833[k]
                    + f_3 * pc_x[k] * lsi_833[k];

        t_1066[k] = f_13 * ksi_834[k]
                    + f_3 * pc_x[k] * lsi_834[k];

        t_1067[k] = f_13 * ksi_835[k]
                    + f_3 * pc_x[k] * lsi_835[k];

        t_1068[k] = f_13 * ksi_836[k]
                    + f_3 * pc_x[k] * lsi_836[k];

        t_1069[k] = f_13 * ksi_837[k]
                    + f_3 * pc_x[k] * lsi_837[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, t_1073, pa_x, pc_x, pc_z, ksk0_1072, ksi_609, \
                         ksi_838, ksi_839, ksk1_1072, lsi_833, lsi_838, \
                         lsi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_13 * ksi_838[k]
                    + f_3 * pc_x[k] * lsi_838[k];

        t_1071[k] = f_13 * ksi_839[k]
                    + f_3 * pc_x[k] * lsi_839[k];

        t_1072[k] = pa_x[k] * ksk0_1072[k]
                    - f_12 * pc_x[k] * ksk1_1072[k];

        t_1073[k] = f_13 * ksi_609[k]
                    + f_3 * pc_z[k] * lsi_833[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, t_1077, pa_x, pc_x, ksk0_1074, ksk0_1075, \
                         ksk0_1076, ksk0_1077, ksk1_1074, ksk1_1075, ksk1_1076, \
                         ksk1_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = pa_x[k] * ksk0_1074[k]
                    - f_12 * pc_x[k] * ksk1_1074[k];

        t_1075[k] = pa_x[k] * ksk0_1075[k]
                    - f_12 * pc_x[k] * ksk1_1075[k];

        t_1076[k] = pa_x[k] * ksk0_1076[k]
                    - f_12 * pc_x[k] * ksk1_1076[k];

        t_1077[k] = pa_x[k] * ksk0_1077[k]
                    - f_12 * pc_x[k] * ksk1_1077[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, t_1081, pa_x, pc_x, pc_y, ksk0_1079, \
                         ksk0_1080, ksi_643, ksi_644, ksi_840, ksk1_1079, ksk1_1080, lsi_839, \
                         lsi_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = f_21 * ksi_643[k]
                    + f_3 * pc_y[k] * lsi_839[k];

        t_1079[k] = pa_x[k] * ksk0_1079[k]
                    - f_12 * pc_x[k] * ksk1_1079[k];

        t_1080[k] = pa_x[k] * ksk0_1080[k]
                    + f_18 * ksi_840[k]
                    - f_12 * pc_x[k] * ksk1_1080[k];

        t_1081[k] = f_17 * ksi_644[k]
                    + f_3 * pc_y[k] * lsi_840[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, pa_x, pc_x, pc_y, pc_z, ksk0_1083, ksi_616, \
                         ksi_646, ksi_843, ksk1_1083, lsi_840, \
                         lsi_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_14 * ksi_616[k]
                    + f_3 * pc_z[k] * lsi_840[k];

        t_1083[k] = pa_x[k] * ksk0_1083[k]
                    + f_17 * ksi_843[k]
                    - f_12 * pc_x[k] * ksk1_1083[k];

        t_1084[k] = f_17 * ksi_646[k]
                    + f_3 * pc_y[k] * lsi_842[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pa_x, pc_x, pc_z, ksk0_1085, ksk0_1086, \
                         ksi_619, ksi_845, ksi_846, ksk1_1085, ksk1_1086, \
                         lsi_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = pa_x[k] * ksk0_1085[k]
                    + f_17 * ksi_845[k]
                    - f_12 * pc_x[k] * ksk1_1085[k];

        t_1086[k] = pa_x[k] * ksk0_1086[k]
                    + f_16 * ksi_846[k]
                    - f_12 * pc_x[k] * ksk1_1086[k];

        t_1087[k] = f_14 * ksi_619[k]
                    + f_3 * pc_z[k] * lsi_843[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pa_x, pc_x, pc_y, ksk0_1089, ksk0_1090, \
                         ksi_649, ksi_849, ksi_850, ksk1_1089, ksk1_1090, \
                         lsi_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * ksi_649[k]
                    + f_3 * pc_y[k] * lsi_845[k];

        t_1089[k] = pa_x[k] * ksk0_1089[k]
                    + f_16 * ksi_849[k]
                    - f_12 * pc_x[k] * ksk1_1089[k];

        t_1090[k] = pa_x[k] * ksk0_1090[k]
                    + f_15 * ksi_850[k]
                    - f_12 * pc_x[k] * ksk1_1090[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pa_x, pc_x, pc_y, pc_z, ksk0_1092, ksi_622, \
                         ksi_653, ksi_852, ksk1_1092, lsi_846, \
                         lsi_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_14 * ksi_622[k]
                    + f_3 * pc_z[k] * lsi_846[k];

        t_1092[k] = pa_x[k] * ksk0_1092[k]
                    + f_15 * ksi_852[k]
                    - f_12 * pc_x[k] * ksk1_1092[k];

        t_1093[k] = f_17 * ksi_653[k]
                    + f_3 * pc_y[k] * lsi_849[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pa_x, pc_x, pc_z, ksk0_1094, ksk0_1095, \
                         ksi_626, ksi_854, ksi_855, ksk1_1094, ksk1_1095, \
                         lsi_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = pa_x[k] * ksk0_1094[k]
                    + f_15 * ksi_854[k]
                    - f_12 * pc_x[k] * ksk1_1094[k];

        t_1095[k] = pa_x[k] * ksk0_1095[k]
                    + f_14 * ksi_855[k]
                    - f_12 * pc_x[k] * ksk1_1095[k];

        t_1096[k] = f_14 * ksi_626[k]
                    + f_3 * pc_z[k] * lsi_850[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pa_x, pc_x, pc_y, ksk0_1097, ksk0_1098, \
                         ksi_658, ksi_857, ksi_858, ksk1_1097, ksk1_1098, \
                         lsi_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = pa_x[k] * ksk0_1097[k]
                    + f_14 * ksi_857[k]
                    - f_12 * pc_x[k] * ksk1_1097[k];

        t_1098[k] = pa_x[k] * ksk0_1098[k]
                    + f_14 * ksi_858[k]
                    - f_12 * pc_x[k] * ksk1_1098[k];

        t_1099[k] = f_17 * ksi_658[k]
                    + f_3 * pc_y[k] * lsi_854[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pa_x, pc_x, ksk0_1100, ksi_860, \
                         ksi_861, ksi_862, ksi_863, ksk1_1100, lsi_861, lsi_862, \
                         lsi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = pa_x[k] * ksk0_1100[k]
                    + f_14 * ksi_860[k]
                    - f_12 * pc_x[k] * ksk1_1100[k];

        t_1101[k] = f_13 * ksi_861[k]
                    + f_3 * pc_x[k] * lsi_861[k];

        t_1102[k] = f_13 * ksi_862[k]
                    + f_3 * pc_x[k] * lsi_862[k];

        t_1103[k] = f_13 * ksi_863[k]
                    + f_3 * pc_x[k] * lsi_863[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, ksi_864, ksi_865, ksi_866, \
                         ksi_867, lsi_864, lsi_865, lsi_866, lsi_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_13 * ksi_864[k]
                    + f_3 * pc_x[k] * lsi_864[k];

        t_1105[k] = f_13 * ksi_865[k]
                    + f_3 * pc_x[k] * lsi_865[k];

        t_1106[k] = f_13 * ksi_866[k]
                    + f_3 * pc_x[k] * lsi_866[k];

        t_1107[k] = f_13 * ksi_867[k]
                    + f_3 * pc_x[k] * lsi_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, t_1111, pa_x, pc_x, pc_z, ksk0_1108, \
                         ksk0_1110, ksk0_1111, ksi_637, ksk1_1108, ksk1_1110, ksk1_1111, \
                         lsi_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = pa_x[k] * ksk0_1108[k]
                    - f_12 * pc_x[k] * ksk1_1108[k];

        t_1109[k] = f_14 * ksi_637[k]
                    + f_3 * pc_z[k] * lsi_861[k];

        t_1110[k] = pa_x[k] * ksk0_1110[k]
                    - f_12 * pc_x[k] * ksk1_1110[k];

        t_1111[k] = pa_x[k] * ksk0_1111[k]
                    - f_12 * pc_x[k] * ksk1_1111[k];
    }

#pragma omp simd aligned(t_1112, t_1113, t_1114, t_1115, pa_x, pc_x, pc_y, ksk0_1112, \
                         ksk0_1113, ksk0_1115, ksi_671, ksk1_1112, ksk1_1113, ksk1_1115, \
                         lsi_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1112[k] = pa_x[k] * ksk0_1112[k]
                    - f_12 * pc_x[k] * ksk1_1112[k];

        t_1113[k] = pa_x[k] * ksk0_1113[k]
                    - f_12 * pc_x[k] * ksk1_1113[k];

        t_1114[k] = f_17 * ksi_671[k]
                    + f_3 * pc_y[k] * lsi_867[k];

        t_1115[k] = pa_x[k] * ksk0_1115[k]
                    - f_12 * pc_x[k] * ksk1_1115[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, pa_x, pc_x, pc_y, pc_z, ksk0_1116, ksi_644, \
                         ksi_672, ksi_868, ksk1_1116, lsi_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = pa_x[k] * ksk0_1116[k]
                    + f_18 * ksi_868[k]
                    - f_12 * pc_x[k] * ksk1_1116[k];

        t_1117[k] = f_16 * ksi_672[k]
                    + f_3 * pc_y[k] * lsi_868[k];

        t_1118[k] = f_15 * ksi_644[k]
                    + f_3 * pc_z[k] * lsi_868[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, pa_x, pc_x, pc_y, ksk0_1119, ksk0_1121, \
                         ksi_674, ksi_871, ksi_873, ksk1_1119, ksk1_1121, \
                         lsi_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = pa_x[k] * ksk0_1119[k]
                    + f_17 * ksi_871[k]
                    - f_12 * pc_x[k] * ksk1_1119[k];

        t_1120[k] = f_16 * ksi_674[k]
                    + f_3 * pc_y[k] * lsi_870[k];

        t_1121[k] = pa_x[k] * ksk0_1121[k]
                    + f_17 * ksi_873[k]
                    - f_12 * pc_x[k] * ksk1_1121[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, pa_x, pc_x, pc_y, pc_z, ksk0_1122, ksi_647, \
                         ksi_677, ksi_874, ksk1_1122, lsi_871, \
                         lsi_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = pa_x[k] * ksk0_1122[k]
                    + f_16 * ksi_874[k]
                    - f_12 * pc_x[k] * ksk1_1122[k];

        t_1123[k] = f_15 * ksi_647[k]
                    + f_3 * pc_z[k] * lsi_871[k];

        t_1124[k] = f_16 * ksi_677[k]
                    + f_3 * pc_y[k] * lsi_873[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, pa_x, pc_x, pc_z, ksk0_1125, ksk0_1126, \
                         ksi_650, ksi_877, ksi_878, ksk1_1125, ksk1_1126, \
                         lsi_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = pa_x[k] * ksk0_1125[k]
                    + f_16 * ksi_877[k]
                    - f_12 * pc_x[k] * ksk1_1125[k];

        t_1126[k] = pa_x[k] * ksk0_1126[k]
                    + f_15 * ksi_878[k]
                    - f_12 * pc_x[k] * ksk1_1126[k];

        t_1127[k] = f_15 * ksi_650[k]
                    + f_3 * pc_z[k] * lsi_874[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pa_x, pc_x, pc_y, ksk0_1128, ksk0_1130, \
                         ksi_681, ksi_880, ksi_882, ksk1_1128, ksk1_1130, \
                         lsi_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = pa_x[k] * ksk0_1128[k]
                    + f_15 * ksi_880[k]
                    - f_12 * pc_x[k] * ksk1_1128[k];

        t_1129[k] = f_16 * ksi_681[k]
                    + f_3 * pc_y[k] * lsi_877[k];

        t_1130[k] = pa_x[k] * ksk0_1130[k]
                    + f_15 * ksi_882[k]
                    - f_12 * pc_x[k] * ksk1_1130[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pa_x, pc_x, pc_z, ksk0_1131, ksk0_1133, \
                         ksi_654, ksi_883, ksi_885, ksk1_1131, ksk1_1133, \
                         lsi_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = pa_x[k] * ksk0_1131[k]
                    + f_14 * ksi_883[k]
                    - f_12 * pc_x[k] * ksk1_1131[k];

        t_1132[k] = f_15 * ksi_654[k]
                    + f_3 * pc_z[k] * lsi_878[k];

        t_1133[k] = pa_x[k] * ksk0_1133[k]
                    + f_14 * ksi_885[k]
                    - f_12 * pc_x[k] * ksk1_1133[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, pa_x, pc_x, pc_y, ksk0_1134, ksk0_1136, \
                         ksi_686, ksi_886, ksi_888, ksk1_1134, ksk1_1136, \
                         lsi_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = pa_x[k] * ksk0_1134[k]
                    + f_14 * ksi_886[k]
                    - f_12 * pc_x[k] * ksk1_1134[k];

        t_1135[k] = f_16 * ksi_686[k]
                    + f_3 * pc_y[k] * lsi_882[k];

        t_1136[k] = pa_x[k] * ksk0_1136[k]
                    + f_14 * ksi_888[k]
                    - f_12 * pc_x[k] * ksk1_1136[k];
    }

#pragma omp simd aligned(t_1137, t_1138, t_1139, t_1140, t_1141, pc_x, ksi_889, ksi_890, \
                         ksi_891, ksi_892, ksi_893, lsi_889, lsi_890, lsi_891, lsi_892, \
                         lsi_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = f_13 * ksi_889[k]
                    + f_3 * pc_x[k] * lsi_889[k];

        t_1138[k] = f_13 * ksi_890[k]
                    + f_3 * pc_x[k] * lsi_890[k];

        t_1139[k] = f_13 * ksi_891[k]
                    + f_3 * pc_x[k] * lsi_891[k];

        t_1140[k] = f_13 * ksi_892[k]
                    + f_3 * pc_x[k] * lsi_892[k];

        t_1141[k] = f_13 * ksi_893[k]
                    + f_3 * pc_x[k] * lsi_893[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, t_1145, pa_x, pc_x, pc_z, ksk0_1144, ksi_665, \
                         ksi_894, ksi_895, ksk1_1144, lsi_889, lsi_894, \
                         lsi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_13 * ksi_894[k]
                    + f_3 * pc_x[k] * lsi_894[k];

        t_1143[k] = f_13 * ksi_895[k]
                    + f_3 * pc_x[k] * lsi_895[k];

        t_1144[k] = pa_x[k] * ksk0_1144[k]
                    - f_12 * pc_x[k] * ksk1_1144[k];

        t_1145[k] = f_15 * ksi_665[k]
                    + f_3 * pc_z[k] * lsi_889[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, t_1149, pa_x, pc_x, ksk0_1146, ksk0_1147, \
                         ksk0_1148, ksk0_1149, ksk1_1146, ksk1_1147, ksk1_1148, \
                         ksk1_1149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = pa_x[k] * ksk0_1146[k]
                    - f_12 * pc_x[k] * ksk1_1146[k];

        t_1147[k] = pa_x[k] * ksk0_1147[k]
                    - f_12 * pc_x[k] * ksk1_1147[k];

        t_1148[k] = pa_x[k] * ksk0_1148[k]
                    - f_12 * pc_x[k] * ksk1_1148[k];

        t_1149[k] = pa_x[k] * ksk0_1149[k]
                    - f_12 * pc_x[k] * ksk1_1149[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, t_1153, pa_x, pc_x, pc_y, ksk0_1151, \
                         ksk0_1152, ksi_699, ksi_700, ksi_896, ksk1_1151, ksk1_1152, lsi_895, \
                         lsi_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_16 * ksi_699[k]
                    + f_3 * pc_y[k] * lsi_895[k];

        t_1151[k] = pa_x[k] * ksk0_1151[k]
                    - f_12 * pc_x[k] * ksk1_1151[k];

        t_1152[k] = pa_x[k] * ksk0_1152[k]
                    + f_18 * ksi_896[k]
                    - f_12 * pc_x[k] * ksk1_1152[k];

        t_1153[k] = f_15 * ksi_700[k]
                    + f_3 * pc_y[k] * lsi_896[k];
    }

#pragma omp simd aligned(t_1154, t_1155, t_1156, pa_x, pc_x, pc_y, pc_z, ksk0_1155, ksi_672, \
                         ksi_702, ksi_899, ksk1_1155, lsi_896, \
                         lsi_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1154[k] = f_16 * ksi_672[k]
                    + f_3 * pc_z[k] * lsi_896[k];

        t_1155[k] = pa_x[k] * ksk0_1155[k]
                    + f_17 * ksi_899[k]
                    - f_12 * pc_x[k] * ksk1_1155[k];

        t_1156[k] = f_15 * ksi_702[k]
                    + f_3 * pc_y[k] * lsi_898[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pa_x, pc_x, pc_z, ksk0_1157, ksk0_1158, \
                         ksi_675, ksi_901, ksi_902, ksk1_1157, ksk1_1158, \
                         lsi_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = pa_x[k] * ksk0_1157[k]
                    + f_17 * ksi_901[k]
                    - f_12 * pc_x[k] * ksk1_1157[k];

        t_1158[k] = pa_x[k] * ksk0_1158[k]
                    + f_16 * ksi_902[k]
                    - f_12 * pc_x[k] * ksk1_1158[k];

        t_1159[k] = f_16 * ksi_675[k]
                    + f_3 * pc_z[k] * lsi_899[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pa_x, pc_x, pc_y, ksk0_1161, ksk0_1162, \
                         ksi_705, ksi_905, ksi_906, ksk1_1161, ksk1_1162, \
                         lsi_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_15 * ksi_705[k]
                    + f_3 * pc_y[k] * lsi_901[k];

        t_1161[k] = pa_x[k] * ksk0_1161[k]
                    + f_16 * ksi_905[k]
                    - f_12 * pc_x[k] * ksk1_1161[k];

        t_1162[k] = pa_x[k] * ksk0_1162[k]
                    + f_15 * ksi_906[k]
                    - f_12 * pc_x[k] * ksk1_1162[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pa_x, pc_x, pc_y, pc_z, ksk0_1164, ksi_678, \
                         ksi_709, ksi_908, ksk1_1164, lsi_902, \
                         lsi_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * ksi_678[k]
                    + f_3 * pc_z[k] * lsi_902[k];

        t_1164[k] = pa_x[k] * ksk0_1164[k]
                    + f_15 * ksi_908[k]
                    - f_12 * pc_x[k] * ksk1_1164[k];

        t_1165[k] = f_15 * ksi_709[k]
                    + f_3 * pc_y[k] * lsi_905[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pa_x, pc_x, pc_z, ksk0_1166, ksk0_1167, \
                         ksi_682, ksi_910, ksi_911, ksk1_1166, ksk1_1167, \
                         lsi_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = pa_x[k] * ksk0_1166[k]
                    + f_15 * ksi_910[k]
                    - f_12 * pc_x[k] * ksk1_1166[k];

        t_1167[k] = pa_x[k] * ksk0_1167[k]
                    + f_14 * ksi_911[k]
                    - f_12 * pc_x[k] * ksk1_1167[k];

        t_1168[k] = f_16 * ksi_682[k]
                    + f_3 * pc_z[k] * lsi_906[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksk0,
                                                           const size_t ksi, const size_t ksk1,
                                                           const size_t lsh0, const size_t lsh1,
                                                           const size_t lsi, const size_t ncols,
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
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *ksk0_972 = buffer.data(ksk0 + 972);
    const auto *ksk0_977 = buffer.data(ksk0 + 977);
    const auto *ksk0_981 = buffer.data(ksk0 + 981);
    const auto *ksk0_986 = buffer.data(ksk0 + 986);
    const auto *ksk0_992 = buffer.data(ksk0 + 992);
    const auto *ksk0_1169 = buffer.data(ksk0 + 1169);
    const auto *ksk0_1170 = buffer.data(ksk0 + 1170);
    const auto *ksk0_1172 = buffer.data(ksk0 + 1172);
    const auto *ksk0_1180 = buffer.data(ksk0 + 1180);
    const auto *ksk0_1182 = buffer.data(ksk0 + 1182);
    const auto *ksk0_1183 = buffer.data(ksk0 + 1183);
    const auto *ksk0_1184 = buffer.data(ksk0 + 1184);
    const auto *ksk0_1185 = buffer.data(ksk0 + 1185);
    const auto *ksk0_1187 = buffer.data(ksk0 + 1187);
    const auto *ksk0_1188 = buffer.data(ksk0 + 1188);
    const auto *ksk0_1191 = buffer.data(ksk0 + 1191);
    const auto *ksk0_1193 = buffer.data(ksk0 + 1193);
    const auto *ksk0_1194 = buffer.data(ksk0 + 1194);
    const auto *ksk0_1197 = buffer.data(ksk0 + 1197);
    const auto *ksk0_1198 = buffer.data(ksk0 + 1198);
    const auto *ksk0_1200 = buffer.data(ksk0 + 1200);
    const auto *ksk0_1202 = buffer.data(ksk0 + 1202);
    const auto *ksk0_1203 = buffer.data(ksk0 + 1203);
    const auto *ksk0_1205 = buffer.data(ksk0 + 1205);
    const auto *ksk0_1206 = buffer.data(ksk0 + 1206);
    const auto *ksk0_1208 = buffer.data(ksk0 + 1208);
    const auto *ksk0_1216 = buffer.data(ksk0 + 1216);
    const auto *ksk0_1218 = buffer.data(ksk0 + 1218);
    const auto *ksk0_1219 = buffer.data(ksk0 + 1219);
    const auto *ksk0_1220 = buffer.data(ksk0 + 1220);
    const auto *ksk0_1221 = buffer.data(ksk0 + 1221);
    const auto *ksk0_1223 = buffer.data(ksk0 + 1223);
    const auto *ksk0_1227 = buffer.data(ksk0 + 1227);
    const auto *ksk0_1230 = buffer.data(ksk0 + 1230);
    const auto *ksk0_1234 = buffer.data(ksk0 + 1234);
    const auto *ksk0_1236 = buffer.data(ksk0 + 1236);
    const auto *ksk0_1239 = buffer.data(ksk0 + 1239);
    const auto *ksk0_1241 = buffer.data(ksk0 + 1241);
    const auto *ksk0_1242 = buffer.data(ksk0 + 1242);
    const auto *ksk0_1252 = buffer.data(ksk0 + 1252);
    const auto *ksk0_1254 = buffer.data(ksk0 + 1254);
    const auto *ksk0_1255 = buffer.data(ksk0 + 1255);
    const auto *ksk0_1256 = buffer.data(ksk0 + 1256);
    const auto *ksk0_1257 = buffer.data(ksk0 + 1257);
    const auto *ksk0_1259 = buffer.data(ksk0 + 1259);
    const auto *ksk0_1260 = buffer.data(ksk0 + 1260);
    const auto *ksk0_1265 = buffer.data(ksk0 + 1265);
    const auto *ksk0_1269 = buffer.data(ksk0 + 1269);
    const auto *ksk0_1274 = buffer.data(ksk0 + 1274);
    const auto *ksk0_1280 = buffer.data(ksk0 + 1280);

    const auto *ksi_693 = buffer.data(ksi + 693);
    const auto *ksi_700 = buffer.data(ksi + 700);
    const auto *ksi_703 = buffer.data(ksi + 703);
    const auto *ksi_706 = buffer.data(ksi + 706);
    const auto *ksi_710 = buffer.data(ksi + 710);
    const auto *ksi_714 = buffer.data(ksi + 714);
    const auto *ksi_721 = buffer.data(ksi + 721);
    const auto *ksi_727 = buffer.data(ksi + 727);
    const auto *ksi_728 = buffer.data(ksi + 728);
    const auto *ksi_730 = buffer.data(ksi + 730);
    const auto *ksi_731 = buffer.data(ksi + 731);
    const auto *ksi_733 = buffer.data(ksi + 733);
    const auto *ksi_734 = buffer.data(ksi + 734);
    const auto *ksi_737 = buffer.data(ksi + 737);
    const auto *ksi_738 = buffer.data(ksi + 738);
    const auto *ksi_742 = buffer.data(ksi + 742);
    const auto *ksi_749 = buffer.data(ksi + 749);
    const auto *ksi_755 = buffer.data(ksi + 755);
    const auto *ksi_756 = buffer.data(ksi + 756);
    const auto *ksi_758 = buffer.data(ksi + 758);
    const auto *ksi_761 = buffer.data(ksi + 761);
    const auto *ksi_765 = buffer.data(ksi + 765);
    const auto *ksi_770 = buffer.data(ksi + 770);
    const auto *ksi_783 = buffer.data(ksi + 783);
    const auto *ksi_913 = buffer.data(ksi + 913);
    const auto *ksi_914 = buffer.data(ksi + 914);
    const auto *ksi_916 = buffer.data(ksi + 916);
    const auto *ksi_917 = buffer.data(ksi + 917);
    const auto *ksi_918 = buffer.data(ksi + 918);
    const auto *ksi_919 = buffer.data(ksi + 919);
    const auto *ksi_920 = buffer.data(ksi + 920);
    const auto *ksi_921 = buffer.data(ksi + 921);
    const auto *ksi_922 = buffer.data(ksi + 922);
    const auto *ksi_923 = buffer.data(ksi + 923);
    const auto *ksi_924 = buffer.data(ksi + 924);
    const auto *ksi_927 = buffer.data(ksi + 927);
    const auto *ksi_929 = buffer.data(ksi + 929);
    const auto *ksi_930 = buffer.data(ksi + 930);
    const auto *ksi_933 = buffer.data(ksi + 933);
    const auto *ksi_934 = buffer.data(ksi + 934);
    const auto *ksi_936 = buffer.data(ksi + 936);
    const auto *ksi_938 = buffer.data(ksi + 938);
    const auto *ksi_939 = buffer.data(ksi + 939);
    const auto *ksi_941 = buffer.data(ksi + 941);
    const auto *ksi_942 = buffer.data(ksi + 942);
    const auto *ksi_944 = buffer.data(ksi + 944);
    const auto *ksi_945 = buffer.data(ksi + 945);
    const auto *ksi_946 = buffer.data(ksi + 946);
    const auto *ksi_947 = buffer.data(ksi + 947);
    const auto *ksi_948 = buffer.data(ksi + 948);
    const auto *ksi_949 = buffer.data(ksi + 949);
    const auto *ksi_950 = buffer.data(ksi + 950);
    const auto *ksi_951 = buffer.data(ksi + 951);
    const auto *ksi_955 = buffer.data(ksi + 955);
    const auto *ksi_958 = buffer.data(ksi + 958);
    const auto *ksi_962 = buffer.data(ksi + 962);
    const auto *ksi_964 = buffer.data(ksi + 964);
    const auto *ksi_967 = buffer.data(ksi + 967);
    const auto *ksi_969 = buffer.data(ksi + 969);
    const auto *ksi_970 = buffer.data(ksi + 970);
    const auto *ksi_973 = buffer.data(ksi + 973);
    const auto *ksi_974 = buffer.data(ksi + 974);
    const auto *ksi_975 = buffer.data(ksi + 975);
    const auto *ksi_976 = buffer.data(ksi + 976);
    const auto *ksi_977 = buffer.data(ksi + 977);
    const auto *ksi_978 = buffer.data(ksi + 978);
    const auto *ksi_979 = buffer.data(ksi + 979);
    const auto *ksi_980 = buffer.data(ksi + 980);
    const auto *ksi_985 = buffer.data(ksi + 985);
    const auto *ksi_989 = buffer.data(ksi + 989);
    const auto *ksi_994 = buffer.data(ksi + 994);
    const auto *ksi_1000 = buffer.data(ksi + 1000);
    const auto *ksi_1001 = buffer.data(ksi + 1001);
    const auto *ksi_1002 = buffer.data(ksi + 1002);
    const auto *ksi_1003 = buffer.data(ksi + 1003);
    const auto *ksi_1004 = buffer.data(ksi + 1004);
    const auto *ksi_1005 = buffer.data(ksi + 1005);
    const auto *ksi_1007 = buffer.data(ksi + 1007);

    const auto *ksk1_972 = buffer.data(ksk1 + 972);
    const auto *ksk1_977 = buffer.data(ksk1 + 977);
    const auto *ksk1_981 = buffer.data(ksk1 + 981);
    const auto *ksk1_986 = buffer.data(ksk1 + 986);
    const auto *ksk1_992 = buffer.data(ksk1 + 992);
    const auto *ksk1_1169 = buffer.data(ksk1 + 1169);
    const auto *ksk1_1170 = buffer.data(ksk1 + 1170);
    const auto *ksk1_1172 = buffer.data(ksk1 + 1172);
    const auto *ksk1_1180 = buffer.data(ksk1 + 1180);
    const auto *ksk1_1182 = buffer.data(ksk1 + 1182);
    const auto *ksk1_1183 = buffer.data(ksk1 + 1183);
    const auto *ksk1_1184 = buffer.data(ksk1 + 1184);
    const auto *ksk1_1185 = buffer.data(ksk1 + 1185);
    const auto *ksk1_1187 = buffer.data(ksk1 + 1187);
    const auto *ksk1_1188 = buffer.data(ksk1 + 1188);
    const auto *ksk1_1191 = buffer.data(ksk1 + 1191);
    const auto *ksk1_1193 = buffer.data(ksk1 + 1193);
    const auto *ksk1_1194 = buffer.data(ksk1 + 1194);
    const auto *ksk1_1197 = buffer.data(ksk1 + 1197);
    const auto *ksk1_1198 = buffer.data(ksk1 + 1198);
    const auto *ksk1_1200 = buffer.data(ksk1 + 1200);
    const auto *ksk1_1202 = buffer.data(ksk1 + 1202);
    const auto *ksk1_1203 = buffer.data(ksk1 + 1203);
    const auto *ksk1_1205 = buffer.data(ksk1 + 1205);
    const auto *ksk1_1206 = buffer.data(ksk1 + 1206);
    const auto *ksk1_1208 = buffer.data(ksk1 + 1208);
    const auto *ksk1_1216 = buffer.data(ksk1 + 1216);
    const auto *ksk1_1218 = buffer.data(ksk1 + 1218);
    const auto *ksk1_1219 = buffer.data(ksk1 + 1219);
    const auto *ksk1_1220 = buffer.data(ksk1 + 1220);
    const auto *ksk1_1221 = buffer.data(ksk1 + 1221);
    const auto *ksk1_1223 = buffer.data(ksk1 + 1223);
    const auto *ksk1_1227 = buffer.data(ksk1 + 1227);
    const auto *ksk1_1230 = buffer.data(ksk1 + 1230);
    const auto *ksk1_1234 = buffer.data(ksk1 + 1234);
    const auto *ksk1_1236 = buffer.data(ksk1 + 1236);
    const auto *ksk1_1239 = buffer.data(ksk1 + 1239);
    const auto *ksk1_1241 = buffer.data(ksk1 + 1241);
    const auto *ksk1_1242 = buffer.data(ksk1 + 1242);
    const auto *ksk1_1252 = buffer.data(ksk1 + 1252);
    const auto *ksk1_1254 = buffer.data(ksk1 + 1254);
    const auto *ksk1_1255 = buffer.data(ksk1 + 1255);
    const auto *ksk1_1256 = buffer.data(ksk1 + 1256);
    const auto *ksk1_1257 = buffer.data(ksk1 + 1257);
    const auto *ksk1_1259 = buffer.data(ksk1 + 1259);
    const auto *ksk1_1260 = buffer.data(ksk1 + 1260);
    const auto *ksk1_1265 = buffer.data(ksk1 + 1265);
    const auto *ksk1_1269 = buffer.data(ksk1 + 1269);
    const auto *ksk1_1274 = buffer.data(ksk1 + 1274);
    const auto *ksk1_1280 = buffer.data(ksk1 + 1280);

    const auto *lsh0_735 = buffer.data(lsh0 + 735);
    const auto *lsh0_736 = buffer.data(lsh0 + 736);
    const auto *lsh0_737 = buffer.data(lsh0 + 737);
    const auto *lsh0_738 = buffer.data(lsh0 + 738);
    const auto *lsh0_739 = buffer.data(lsh0 + 739);
    const auto *lsh0_740 = buffer.data(lsh0 + 740);
    const auto *lsh0_741 = buffer.data(lsh0 + 741);
    const auto *lsh0_742 = buffer.data(lsh0 + 742);
    const auto *lsh0_743 = buffer.data(lsh0 + 743);
    const auto *lsh0_744 = buffer.data(lsh0 + 744);

    const auto *lsh1_735 = buffer.data(lsh1 + 735);
    const auto *lsh1_736 = buffer.data(lsh1 + 736);
    const auto *lsh1_737 = buffer.data(lsh1 + 737);
    const auto *lsh1_738 = buffer.data(lsh1 + 738);
    const auto *lsh1_739 = buffer.data(lsh1 + 739);
    const auto *lsh1_740 = buffer.data(lsh1 + 740);
    const auto *lsh1_741 = buffer.data(lsh1 + 741);
    const auto *lsh1_742 = buffer.data(lsh1 + 742);
    const auto *lsh1_743 = buffer.data(lsh1 + 743);
    const auto *lsh1_744 = buffer.data(lsh1 + 744);

    const auto *lsi_910 = buffer.data(lsi + 910);
    const auto *lsi_917 = buffer.data(lsi + 917);
    const auto *lsi_918 = buffer.data(lsi + 918);
    const auto *lsi_919 = buffer.data(lsi + 919);
    const auto *lsi_920 = buffer.data(lsi + 920);
    const auto *lsi_921 = buffer.data(lsi + 921);
    const auto *lsi_922 = buffer.data(lsi + 922);
    const auto *lsi_923 = buffer.data(lsi + 923);
    const auto *lsi_924 = buffer.data(lsi + 924);
    const auto *lsi_926 = buffer.data(lsi + 926);
    const auto *lsi_927 = buffer.data(lsi + 927);
    const auto *lsi_929 = buffer.data(lsi + 929);
    const auto *lsi_930 = buffer.data(lsi + 930);
    const auto *lsi_933 = buffer.data(lsi + 933);
    const auto *lsi_934 = buffer.data(lsi + 934);
    const auto *lsi_938 = buffer.data(lsi + 938);
    const auto *lsi_945 = buffer.data(lsi + 945);
    const auto *lsi_946 = buffer.data(lsi + 946);
    const auto *lsi_947 = buffer.data(lsi + 947);
    const auto *lsi_948 = buffer.data(lsi + 948);
    const auto *lsi_949 = buffer.data(lsi + 949);
    const auto *lsi_950 = buffer.data(lsi + 950);
    const auto *lsi_951 = buffer.data(lsi + 951);
    const auto *lsi_952 = buffer.data(lsi + 952);
    const auto *lsi_954 = buffer.data(lsi + 954);
    const auto *lsi_955 = buffer.data(lsi + 955);
    const auto *lsi_957 = buffer.data(lsi + 957);
    const auto *lsi_958 = buffer.data(lsi + 958);
    const auto *lsi_961 = buffer.data(lsi + 961);
    const auto *lsi_962 = buffer.data(lsi + 962);
    const auto *lsi_966 = buffer.data(lsi + 966);
    const auto *lsi_973 = buffer.data(lsi + 973);
    const auto *lsi_974 = buffer.data(lsi + 974);
    const auto *lsi_975 = buffer.data(lsi + 975);
    const auto *lsi_976 = buffer.data(lsi + 976);
    const auto *lsi_977 = buffer.data(lsi + 977);
    const auto *lsi_978 = buffer.data(lsi + 978);
    const auto *lsi_979 = buffer.data(lsi + 979);
    const auto *lsi_980 = buffer.data(lsi + 980);
    const auto *lsi_981 = buffer.data(lsi + 981);
    const auto *lsi_982 = buffer.data(lsi + 982);
    const auto *lsi_983 = buffer.data(lsi + 983);
    const auto *lsi_984 = buffer.data(lsi + 984);
    const auto *lsi_985 = buffer.data(lsi + 985);
    const auto *lsi_986 = buffer.data(lsi + 986);
    const auto *lsi_987 = buffer.data(lsi + 987);
    const auto *lsi_988 = buffer.data(lsi + 988);
    const auto *lsi_989 = buffer.data(lsi + 989);
    const auto *lsi_990 = buffer.data(lsi + 990);
    const auto *lsi_991 = buffer.data(lsi + 991);
    const auto *lsi_992 = buffer.data(lsi + 992);
    const auto *lsi_993 = buffer.data(lsi + 993);
    const auto *lsi_994 = buffer.data(lsi + 994);
    const auto *lsi_1000 = buffer.data(lsi + 1000);
    const auto *lsi_1001 = buffer.data(lsi + 1001);
    const auto *lsi_1002 = buffer.data(lsi + 1002);
    const auto *lsi_1003 = buffer.data(lsi + 1003);
    const auto *lsi_1004 = buffer.data(lsi + 1004);
    const auto *lsi_1005 = buffer.data(lsi + 1005);
    const auto *lsi_1007 = buffer.data(lsi + 1007);

#pragma omp simd aligned(t_1169, t_1170, t_1171, pa_x, pc_x, pc_y, ksk0_1169, ksk0_1170, \
                         ksi_714, ksi_913, ksi_914, ksk1_1169, ksk1_1170, \
                         lsi_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = pa_x[k] * ksk0_1169[k]
                    + f_14 * ksi_913[k]
                    - f_12 * pc_x[k] * ksk1_1169[k];

        t_1170[k] = pa_x[k] * ksk0_1170[k]
                    + f_14 * ksi_914[k]
                    - f_12 * pc_x[k] * ksk1_1170[k];

        t_1171[k] = f_15 * ksi_714[k]
                    + f_3 * pc_y[k] * lsi_910[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, pa_x, pc_x, ksk0_1172, ksi_916, \
                         ksi_917, ksi_918, ksi_919, ksk1_1172, lsi_917, lsi_918, \
                         lsi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = pa_x[k] * ksk0_1172[k]
                    + f_14 * ksi_916[k]
                    - f_12 * pc_x[k] * ksk1_1172[k];

        t_1173[k] = f_13 * ksi_917[k]
                    + f_3 * pc_x[k] * lsi_917[k];

        t_1174[k] = f_13 * ksi_918[k]
                    + f_3 * pc_x[k] * lsi_918[k];

        t_1175[k] = f_13 * ksi_919[k]
                    + f_3 * pc_x[k] * lsi_919[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, t_1179, pc_x, ksi_920, ksi_921, ksi_922, \
                         ksi_923, lsi_920, lsi_921, lsi_922, lsi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_13 * ksi_920[k]
                    + f_3 * pc_x[k] * lsi_920[k];

        t_1177[k] = f_13 * ksi_921[k]
                    + f_3 * pc_x[k] * lsi_921[k];

        t_1178[k] = f_13 * ksi_922[k]
                    + f_3 * pc_x[k] * lsi_922[k];

        t_1179[k] = f_13 * ksi_923[k]
                    + f_3 * pc_x[k] * lsi_923[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, t_1183, pa_x, pc_x, pc_z, ksk0_1180, \
                         ksk0_1182, ksk0_1183, ksi_693, ksk1_1180, ksk1_1182, ksk1_1183, \
                         lsi_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = pa_x[k] * ksk0_1180[k]
                    - f_12 * pc_x[k] * ksk1_1180[k];

        t_1181[k] = f_16 * ksi_693[k]
                    + f_3 * pc_z[k] * lsi_917[k];

        t_1182[k] = pa_x[k] * ksk0_1182[k]
                    - f_12 * pc_x[k] * ksk1_1182[k];

        t_1183[k] = pa_x[k] * ksk0_1183[k]
                    - f_12 * pc_x[k] * ksk1_1183[k];
    }

#pragma omp simd aligned(t_1184, t_1185, t_1186, t_1187, pa_x, pc_x, pc_y, ksk0_1184, \
                         ksk0_1185, ksk0_1187, ksi_727, ksk1_1184, ksk1_1185, ksk1_1187, \
                         lsi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1184[k] = pa_x[k] * ksk0_1184[k]
                    - f_12 * pc_x[k] * ksk1_1184[k];

        t_1185[k] = pa_x[k] * ksk0_1185[k]
                    - f_12 * pc_x[k] * ksk1_1185[k];

        t_1186[k] = f_15 * ksi_727[k]
                    + f_3 * pc_y[k] * lsi_923[k];

        t_1187[k] = pa_x[k] * ksk0_1187[k]
                    - f_12 * pc_x[k] * ksk1_1187[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, pa_x, pc_x, pc_y, pc_z, ksk0_1188, ksi_700, \
                         ksi_728, ksi_924, ksk1_1188, lsi_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = pa_x[k] * ksk0_1188[k]
                    + f_18 * ksi_924[k]
                    - f_12 * pc_x[k] * ksk1_1188[k];

        t_1189[k] = f_14 * ksi_728[k]
                    + f_3 * pc_y[k] * lsi_924[k];

        t_1190[k] = f_17 * ksi_700[k]
                    + f_3 * pc_z[k] * lsi_924[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, pa_x, pc_x, pc_y, ksk0_1191, ksk0_1193, \
                         ksi_730, ksi_927, ksi_929, ksk1_1191, ksk1_1193, \
                         lsi_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = pa_x[k] * ksk0_1191[k]
                    + f_17 * ksi_927[k]
                    - f_12 * pc_x[k] * ksk1_1191[k];

        t_1192[k] = f_14 * ksi_730[k]
                    + f_3 * pc_y[k] * lsi_926[k];

        t_1193[k] = pa_x[k] * ksk0_1193[k]
                    + f_17 * ksi_929[k]
                    - f_12 * pc_x[k] * ksk1_1193[k];
    }

#pragma omp simd aligned(t_1194, t_1195, t_1196, pa_x, pc_x, pc_y, pc_z, ksk0_1194, ksi_703, \
                         ksi_733, ksi_930, ksk1_1194, lsi_927, \
                         lsi_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1194[k] = pa_x[k] * ksk0_1194[k]
                    + f_16 * ksi_930[k]
                    - f_12 * pc_x[k] * ksk1_1194[k];

        t_1195[k] = f_17 * ksi_703[k]
                    + f_3 * pc_z[k] * lsi_927[k];

        t_1196[k] = f_14 * ksi_733[k]
                    + f_3 * pc_y[k] * lsi_929[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, pa_x, pc_x, pc_z, ksk0_1197, ksk0_1198, \
                         ksi_706, ksi_933, ksi_934, ksk1_1197, ksk1_1198, \
                         lsi_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = pa_x[k] * ksk0_1197[k]
                    + f_16 * ksi_933[k]
                    - f_12 * pc_x[k] * ksk1_1197[k];

        t_1198[k] = pa_x[k] * ksk0_1198[k]
                    + f_15 * ksi_934[k]
                    - f_12 * pc_x[k] * ksk1_1198[k];

        t_1199[k] = f_17 * ksi_706[k]
                    + f_3 * pc_z[k] * lsi_930[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, pa_x, pc_x, pc_y, ksk0_1200, ksk0_1202, \
                         ksi_737, ksi_936, ksi_938, ksk1_1200, ksk1_1202, \
                         lsi_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = pa_x[k] * ksk0_1200[k]
                    + f_15 * ksi_936[k]
                    - f_12 * pc_x[k] * ksk1_1200[k];

        t_1201[k] = f_14 * ksi_737[k]
                    + f_3 * pc_y[k] * lsi_933[k];

        t_1202[k] = pa_x[k] * ksk0_1202[k]
                    + f_15 * ksi_938[k]
                    - f_12 * pc_x[k] * ksk1_1202[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, pa_x, pc_x, pc_z, ksk0_1203, ksk0_1205, \
                         ksi_710, ksi_939, ksi_941, ksk1_1203, ksk1_1205, \
                         lsi_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = pa_x[k] * ksk0_1203[k]
                    + f_14 * ksi_939[k]
                    - f_12 * pc_x[k] * ksk1_1203[k];

        t_1204[k] = f_17 * ksi_710[k]
                    + f_3 * pc_z[k] * lsi_934[k];

        t_1205[k] = pa_x[k] * ksk0_1205[k]
                    + f_14 * ksi_941[k]
                    - f_12 * pc_x[k] * ksk1_1205[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, pa_x, pc_x, pc_y, ksk0_1206, ksk0_1208, \
                         ksi_742, ksi_942, ksi_944, ksk1_1206, ksk1_1208, \
                         lsi_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = pa_x[k] * ksk0_1206[k]
                    + f_14 * ksi_942[k]
                    - f_12 * pc_x[k] * ksk1_1206[k];

        t_1207[k] = f_14 * ksi_742[k]
                    + f_3 * pc_y[k] * lsi_938[k];

        t_1208[k] = pa_x[k] * ksk0_1208[k]
                    + f_14 * ksi_944[k]
                    - f_12 * pc_x[k] * ksk1_1208[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, t_1212, t_1213, pc_x, ksi_945, ksi_946, \
                         ksi_947, ksi_948, ksi_949, lsi_945, lsi_946, lsi_947, lsi_948, \
                         lsi_949 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = f_13 * ksi_945[k]
                    + f_3 * pc_x[k] * lsi_945[k];

        t_1210[k] = f_13 * ksi_946[k]
                    + f_3 * pc_x[k] * lsi_946[k];

        t_1211[k] = f_13 * ksi_947[k]
                    + f_3 * pc_x[k] * lsi_947[k];

        t_1212[k] = f_13 * ksi_948[k]
                    + f_3 * pc_x[k] * lsi_948[k];

        t_1213[k] = f_13 * ksi_949[k]
                    + f_3 * pc_x[k] * lsi_949[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pa_x, pc_x, pc_z, ksk0_1216, ksi_721, \
                         ksi_950, ksi_951, ksk1_1216, lsi_945, lsi_950, \
                         lsi_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = f_13 * ksi_950[k]
                    + f_3 * pc_x[k] * lsi_950[k];

        t_1215[k] = f_13 * ksi_951[k]
                    + f_3 * pc_x[k] * lsi_951[k];

        t_1216[k] = pa_x[k] * ksk0_1216[k]
                    - f_12 * pc_x[k] * ksk1_1216[k];

        t_1217[k] = f_17 * ksi_721[k]
                    + f_3 * pc_z[k] * lsi_945[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, t_1221, pa_x, pc_x, ksk0_1218, ksk0_1219, \
                         ksk0_1220, ksk0_1221, ksk1_1218, ksk1_1219, ksk1_1220, \
                         ksk1_1221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = pa_x[k] * ksk0_1218[k]
                    - f_12 * pc_x[k] * ksk1_1218[k];

        t_1219[k] = pa_x[k] * ksk0_1219[k]
                    - f_12 * pc_x[k] * ksk1_1219[k];

        t_1220[k] = pa_x[k] * ksk0_1220[k]
                    - f_12 * pc_x[k] * ksk1_1220[k];

        t_1221[k] = pa_x[k] * ksk0_1221[k]
                    - f_12 * pc_x[k] * ksk1_1221[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pa_x, pa_y, pc_x, pc_y, ksk0_972, \
                         ksk0_1223, ksi_755, ksi_756, ksk1_972, ksk1_1223, lsi_951, \
                         lsi_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_14 * ksi_755[k]
                    + f_3 * pc_y[k] * lsi_951[k];

        t_1223[k] = pa_x[k] * ksk0_1223[k]
                    - f_12 * pc_x[k] * ksk1_1223[k];

        t_1224[k] = pa_y[k] * ksk0_972[k]
                    - f_12 * pc_y[k] * ksk1_972[k];

        t_1225[k] = f_13 * ksi_756[k]
                    + f_3 * pc_y[k] * lsi_952[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, pa_x, pc_x, pc_y, pc_z, ksk0_1227, ksi_728, \
                         ksi_758, ksi_955, ksk1_1227, lsi_952, \
                         lsi_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_21 * ksi_728[k]
                    + f_3 * pc_z[k] * lsi_952[k];

        t_1227[k] = pa_x[k] * ksk0_1227[k]
                    + f_17 * ksi_955[k]
                    - f_12 * pc_x[k] * ksk1_1227[k];

        t_1228[k] = f_13 * ksi_758[k]
                    + f_3 * pc_y[k] * lsi_954[k];
    }

#pragma omp simd aligned(t_1229, t_1230, t_1231, pa_x, pa_y, pc_x, pc_y, pc_z, ksk0_977, \
                         ksk0_1230, ksi_731, ksi_958, ksk1_977, ksk1_1230, \
                         lsi_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1229[k] = pa_y[k] * ksk0_977[k]
                    - f_12 * pc_y[k] * ksk1_977[k];

        t_1230[k] = pa_x[k] * ksk0_1230[k]
                    + f_16 * ksi_958[k]
                    - f_12 * pc_x[k] * ksk1_1230[k];

        t_1231[k] = f_21 * ksi_731[k]
                    + f_3 * pc_z[k] * lsi_955[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, pa_x, pa_y, pc_x, pc_y, ksk0_981, ksk0_1234, \
                         ksi_761, ksi_962, ksk1_981, ksk1_1234, \
                         lsi_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = f_13 * ksi_761[k]
                    + f_3 * pc_y[k] * lsi_957[k];

        t_1233[k] = pa_y[k] * ksk0_981[k]
                    - f_12 * pc_y[k] * ksk1_981[k];

        t_1234[k] = pa_x[k] * ksk0_1234[k]
                    + f_15 * ksi_962[k]
                    - f_12 * pc_x[k] * ksk1_1234[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pa_x, pc_x, pc_y, pc_z, ksk0_1236, ksi_734, \
                         ksi_765, ksi_964, ksk1_1236, lsi_958, \
                         lsi_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_21 * ksi_734[k]
                    + f_3 * pc_z[k] * lsi_958[k];

        t_1236[k] = pa_x[k] * ksk0_1236[k]
                    + f_15 * ksi_964[k]
                    - f_12 * pc_x[k] * ksk1_1236[k];

        t_1237[k] = f_13 * ksi_765[k]
                    + f_3 * pc_y[k] * lsi_961[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, pa_x, pa_y, pc_x, pc_y, pc_z, ksk0_986, \
                         ksk0_1239, ksi_738, ksi_967, ksk1_986, ksk1_1239, \
                         lsi_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = pa_y[k] * ksk0_986[k]
                    - f_12 * pc_y[k] * ksk1_986[k];

        t_1239[k] = pa_x[k] * ksk0_1239[k]
                    + f_14 * ksi_967[k]
                    - f_12 * pc_x[k] * ksk1_1239[k];

        t_1240[k] = f_21 * ksi_738[k]
                    + f_3 * pc_z[k] * lsi_962[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pa_x, pc_x, pc_y, ksk0_1241, ksk0_1242, \
                         ksi_770, ksi_969, ksi_970, ksk1_1241, ksk1_1242, \
                         lsi_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = pa_x[k] * ksk0_1241[k]
                    + f_14 * ksi_969[k]
                    - f_12 * pc_x[k] * ksk1_1241[k];

        t_1242[k] = pa_x[k] * ksk0_1242[k]
                    + f_14 * ksi_970[k]
                    - f_12 * pc_x[k] * ksk1_1242[k];

        t_1243[k] = f_13 * ksi_770[k]
                    + f_3 * pc_y[k] * lsi_966[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, pa_y, pc_x, pc_y, ksk0_992, ksi_973, \
                         ksi_974, ksi_975, ksk1_992, lsi_973, lsi_974, \
                         lsi_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = pa_y[k] * ksk0_992[k]
                    - f_12 * pc_y[k] * ksk1_992[k];

        t_1245[k] = f_13 * ksi_973[k]
                    + f_3 * pc_x[k] * lsi_973[k];

        t_1246[k] = f_13 * ksi_974[k]
                    + f_3 * pc_x[k] * lsi_974[k];

        t_1247[k] = f_13 * ksi_975[k]
                    + f_3 * pc_x[k] * lsi_975[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, pc_x, ksi_976, ksi_977, ksi_978, \
                         ksi_979, lsi_976, lsi_977, lsi_978, lsi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = f_13 * ksi_976[k]
                    + f_3 * pc_x[k] * lsi_976[k];

        t_1249[k] = f_13 * ksi_977[k]
                    + f_3 * pc_x[k] * lsi_977[k];

        t_1250[k] = f_13 * ksi_978[k]
                    + f_3 * pc_x[k] * lsi_978[k];

        t_1251[k] = f_13 * ksi_979[k]
                    + f_3 * pc_x[k] * lsi_979[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, pa_x, pc_x, pc_z, ksk0_1252, \
                         ksk0_1254, ksk0_1255, ksi_749, ksk1_1252, ksk1_1254, ksk1_1255, \
                         lsi_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = pa_x[k] * ksk0_1252[k]
                    - f_12 * pc_x[k] * ksk1_1252[k];

        t_1253[k] = f_21 * ksi_749[k]
                    + f_3 * pc_z[k] * lsi_973[k];

        t_1254[k] = pa_x[k] * ksk0_1254[k]
                    - f_12 * pc_x[k] * ksk1_1254[k];

        t_1255[k] = pa_x[k] * ksk0_1255[k]
                    - f_12 * pc_x[k] * ksk1_1255[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pa_x, pc_x, pc_y, ksk0_1256, \
                         ksk0_1257, ksk0_1259, ksi_783, ksk1_1256, ksk1_1257, ksk1_1259, \
                         lsi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = pa_x[k] * ksk0_1256[k]
                    - f_12 * pc_x[k] * ksk1_1256[k];

        t_1257[k] = pa_x[k] * ksk0_1257[k]
                    - f_12 * pc_x[k] * ksk1_1257[k];

        t_1258[k] = f_13 * ksi_783[k]
                    + f_3 * pc_y[k] * lsi_979[k];

        t_1259[k] = pa_x[k] * ksk0_1259[k]
                    - f_12 * pc_x[k] * ksk1_1259[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pa_x, pc_x, pc_y, pc_z, ksk0_1260, \
                         ksi_756, ksi_980, ksk1_1260, lsh0_735, lsh1_735, lsi_980, \
                         lsi_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = pa_x[k] * ksk0_1260[k]
                    + f_18 * ksi_980[k]
                    - f_12 * pc_x[k] * ksk1_1260[k];

        t_1261[k] = f_3 * pc_y[k] * lsi_980[k];

        t_1262[k] = f_18 * ksi_756[k]
                    + f_3 * pc_z[k] * lsi_980[k];

        t_1263[k] = f_4 * lsh0_735[k]
                    - f_5 * lsh1_735[k]
                    + f_3 * pc_y[k] * lsi_981[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, pa_x, pc_x, pc_y, ksk0_1265, ksi_985, \
                         ksk1_1265, lsh0_736, lsh1_736, lsi_982, \
                         lsi_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_3 * pc_y[k] * lsi_982[k];

        t_1265[k] = pa_x[k] * ksk0_1265[k]
                    + f_17 * ksi_985[k]
                    - f_12 * pc_x[k] * ksk1_1265[k];

        t_1266[k] = f_6 * lsh0_736[k]
                    - f_7 * lsh1_736[k]
                    + f_3 * pc_y[k] * lsi_983[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, pa_x, pc_x, pc_y, ksk0_1269, ksi_989, \
                         ksk1_1269, lsh0_737, lsh1_737, lsi_984, \
                         lsi_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_4 * lsh0_737[k]
                    - f_5 * lsh1_737[k]
                    + f_3 * pc_y[k] * lsi_984[k];

        t_1268[k] = f_3 * pc_y[k] * lsi_985[k];

        t_1269[k] = pa_x[k] * ksk0_1269[k]
                    + f_16 * ksi_989[k]
                    - f_12 * pc_x[k] * ksk1_1269[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, t_1273, pc_y, lsh0_738, lsh0_739, lsh0_740, \
                         lsh1_738, lsh1_739, lsh1_740, lsi_986, lsi_987, lsi_988, \
                         lsi_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = f_8 * lsh0_738[k]
                    - f_9 * lsh1_738[k]
                    + f_3 * pc_y[k] * lsi_986[k];

        t_1271[k] = f_6 * lsh0_739[k]
                    - f_7 * lsh1_739[k]
                    + f_3 * pc_y[k] * lsi_987[k];

        t_1272[k] = f_4 * lsh0_740[k]
                    - f_5 * lsh1_740[k]
                    + f_3 * pc_y[k] * lsi_988[k];

        t_1273[k] = f_3 * pc_y[k] * lsi_989[k];
    }

#pragma omp simd aligned(t_1274, t_1275, t_1276, pa_x, pc_x, pc_y, ksk0_1274, ksi_994, \
                         ksk1_1274, lsh0_741, lsh0_742, lsh1_741, lsh1_742, lsi_990, \
                         lsi_991 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1274[k] = pa_x[k] * ksk0_1274[k]
                    + f_15 * ksi_994[k]
                    - f_12 * pc_x[k] * ksk1_1274[k];

        t_1275[k] = f_10 * lsh0_741[k]
                    - f_11 * lsh1_741[k]
                    + f_3 * pc_y[k] * lsi_990[k];

        t_1276[k] = f_8 * lsh0_742[k]
                    - f_9 * lsh1_742[k]
                    + f_3 * pc_y[k] * lsi_991[k];
    }

#pragma omp simd aligned(t_1277, t_1278, t_1279, pc_y, lsh0_743, lsh0_744, lsh1_743, lsh1_744, \
                         lsi_992, lsi_993, lsi_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = f_6 * lsh0_743[k]
                    - f_7 * lsh1_743[k]
                    + f_3 * pc_y[k] * lsi_992[k];

        t_1278[k] = f_4 * lsh0_744[k]
                    - f_5 * lsh1_744[k]
                    + f_3 * pc_y[k] * lsi_993[k];

        t_1279[k] = f_3 * pc_y[k] * lsi_994[k];
    }

#pragma omp simd aligned(t_1280, t_1281, t_1282, t_1283, pa_x, pc_x, ksk0_1280, ksi_1000, \
                         ksi_1001, ksi_1002, ksi_1003, ksk1_1280, lsi_1001, lsi_1002, \
                         lsi_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1280[k] = pa_x[k] * ksk0_1280[k]
                    + f_14 * ksi_1000[k]
                    - f_12 * pc_x[k] * ksk1_1280[k];

        t_1281[k] = f_13 * ksi_1001[k]
                    + f_3 * pc_x[k] * lsi_1001[k];

        t_1282[k] = f_13 * ksi_1002[k]
                    + f_3 * pc_x[k] * lsi_1002[k];

        t_1283[k] = f_13 * ksi_1003[k]
                    + f_3 * pc_x[k] * lsi_1003[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, t_1287, pc_x, pc_y, ksi_1004, ksi_1005, \
                         ksi_1007, lsi_1000, lsi_1004, lsi_1005, \
                         lsi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_13 * ksi_1004[k]
                    + f_3 * pc_x[k] * lsi_1004[k];

        t_1285[k] = f_13 * ksi_1005[k]
                    + f_3 * pc_x[k] * lsi_1005[k];

        t_1286[k] = f_3 * pc_y[k] * lsi_1000[k];

        t_1287[k] = f_13 * ksi_1007[k]
                    + f_3 * pc_x[k] * lsi_1007[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksk0,
                                                           const size_t ksi, const size_t ksk1,
                                                           const size_t lsh0, const size_t lsh1,
                                                           const size_t lsi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.0 / q;

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
    auto *t_1404 = buffer.data(target + 1404);
    auto *t_1405 = buffer.data(target + 1405);
    auto *t_1406 = buffer.data(target + 1406);
    auto *t_1407 = buffer.data(target + 1407);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_1008 = buffer.data(ksk0 + 1008);
    const auto *ksk0_1009 = buffer.data(ksk0 + 1009);
    const auto *ksk0_1011 = buffer.data(ksk0 + 1011);
    const auto *ksk0_1014 = buffer.data(ksk0 + 1014);
    const auto *ksk0_1018 = buffer.data(ksk0 + 1018);
    const auto *ksk0_1023 = buffer.data(ksk0 + 1023);
    const auto *ksk0_1036 = buffer.data(ksk0 + 1036);
    const auto *ksk0_1038 = buffer.data(ksk0 + 1038);
    const auto *ksk0_1039 = buffer.data(ksk0 + 1039);
    const auto *ksk0_1040 = buffer.data(ksk0 + 1040);
    const auto *ksk0_1041 = buffer.data(ksk0 + 1041);
    const auto *ksk0_1288 = buffer.data(ksk0 + 1288);
    const auto *ksk0_1289 = buffer.data(ksk0 + 1289);
    const auto *ksk0_1290 = buffer.data(ksk0 + 1290);
    const auto *ksk0_1291 = buffer.data(ksk0 + 1291);
    const auto *ksk0_1292 = buffer.data(ksk0 + 1292);
    const auto *ksk0_1293 = buffer.data(ksk0 + 1293);
    const auto *ksk0_1295 = buffer.data(ksk0 + 1295);

    const auto *ksi_805 = buffer.data(ksi + 805);
    const auto *ksi_806 = buffer.data(ksi + 806);
    const auto *ksi_807 = buffer.data(ksi + 807);
    const auto *ksi_808 = buffer.data(ksi + 808);
    const auto *ksi_809 = buffer.data(ksi + 809);
    const auto *ksi_811 = buffer.data(ksi + 811);
    const auto *ksi_833 = buffer.data(ksi + 833);
    const auto *ksi_839 = buffer.data(ksi + 839);
    const auto *ksi_861 = buffer.data(ksi + 861);
    const auto *ksi_863 = buffer.data(ksi + 863);
    const auto *ksi_864 = buffer.data(ksi + 864);
    const auto *ksi_865 = buffer.data(ksi + 865);
    const auto *ksi_866 = buffer.data(ksi + 866);
    const auto *ksi_867 = buffer.data(ksi + 867);

    const auto *ksk1_1008 = buffer.data(ksk1 + 1008);
    const auto *ksk1_1009 = buffer.data(ksk1 + 1009);
    const auto *ksk1_1011 = buffer.data(ksk1 + 1011);
    const auto *ksk1_1014 = buffer.data(ksk1 + 1014);
    const auto *ksk1_1018 = buffer.data(ksk1 + 1018);
    const auto *ksk1_1023 = buffer.data(ksk1 + 1023);
    const auto *ksk1_1036 = buffer.data(ksk1 + 1036);
    const auto *ksk1_1038 = buffer.data(ksk1 + 1038);
    const auto *ksk1_1039 = buffer.data(ksk1 + 1039);
    const auto *ksk1_1040 = buffer.data(ksk1 + 1040);
    const auto *ksk1_1041 = buffer.data(ksk1 + 1041);
    const auto *ksk1_1288 = buffer.data(ksk1 + 1288);
    const auto *ksk1_1289 = buffer.data(ksk1 + 1289);
    const auto *ksk1_1290 = buffer.data(ksk1 + 1290);
    const auto *ksk1_1291 = buffer.data(ksk1 + 1291);
    const auto *ksk1_1292 = buffer.data(ksk1 + 1292);
    const auto *ksk1_1293 = buffer.data(ksk1 + 1293);
    const auto *ksk1_1295 = buffer.data(ksk1 + 1295);

    const auto *lsh0_756 = buffer.data(lsh0 + 756);
    const auto *lsh0_757 = buffer.data(lsh0 + 757);
    const auto *lsh0_759 = buffer.data(lsh0 + 759);
    const auto *lsh0_761 = buffer.data(lsh0 + 761);
    const auto *lsh0_762 = buffer.data(lsh0 + 762);
    const auto *lsh0_764 = buffer.data(lsh0 + 764);
    const auto *lsh0_765 = buffer.data(lsh0 + 765);
    const auto *lsh0_766 = buffer.data(lsh0 + 766);
    const auto *lsh0_768 = buffer.data(lsh0 + 768);
    const auto *lsh0_769 = buffer.data(lsh0 + 769);
    const auto *lsh0_770 = buffer.data(lsh0 + 770);
    const auto *lsh0_771 = buffer.data(lsh0 + 771);
    const auto *lsh0_772 = buffer.data(lsh0 + 772);
    const auto *lsh0_773 = buffer.data(lsh0 + 773);
    const auto *lsh0_774 = buffer.data(lsh0 + 774);
    const auto *lsh0_775 = buffer.data(lsh0 + 775);
    const auto *lsh0_776 = buffer.data(lsh0 + 776);
    const auto *lsh0_779 = buffer.data(lsh0 + 779);
    const auto *lsh0_781 = buffer.data(lsh0 + 781);
    const auto *lsh0_782 = buffer.data(lsh0 + 782);
    const auto *lsh0_784 = buffer.data(lsh0 + 784);
    const auto *lsh0_785 = buffer.data(lsh0 + 785);
    const auto *lsh0_786 = buffer.data(lsh0 + 786);
    const auto *lsh0_788 = buffer.data(lsh0 + 788);
    const auto *lsh0_789 = buffer.data(lsh0 + 789);
    const auto *lsh0_790 = buffer.data(lsh0 + 790);
    const auto *lsh0_791 = buffer.data(lsh0 + 791);
    const auto *lsh0_793 = buffer.data(lsh0 + 793);
    const auto *lsh0_794 = buffer.data(lsh0 + 794);
    const auto *lsh0_795 = buffer.data(lsh0 + 795);
    const auto *lsh0_796 = buffer.data(lsh0 + 796);
    const auto *lsh0_797 = buffer.data(lsh0 + 797);
    const auto *lsh0_798 = buffer.data(lsh0 + 798);
    const auto *lsh0_799 = buffer.data(lsh0 + 799);
    const auto *lsh0_800 = buffer.data(lsh0 + 800);
    const auto *lsh0_801 = buffer.data(lsh0 + 801);
    const auto *lsh0_802 = buffer.data(lsh0 + 802);
    const auto *lsh0_803 = buffer.data(lsh0 + 803);
    const auto *lsh0_804 = buffer.data(lsh0 + 804);
    const auto *lsh0_805 = buffer.data(lsh0 + 805);
    const auto *lsh0_806 = buffer.data(lsh0 + 806);
    const auto *lsh0_807 = buffer.data(lsh0 + 807);
    const auto *lsh0_808 = buffer.data(lsh0 + 808);
    const auto *lsh0_809 = buffer.data(lsh0 + 809);
    const auto *lsh0_810 = buffer.data(lsh0 + 810);
    const auto *lsh0_811 = buffer.data(lsh0 + 811);
    const auto *lsh0_812 = buffer.data(lsh0 + 812);
    const auto *lsh0_813 = buffer.data(lsh0 + 813);
    const auto *lsh0_814 = buffer.data(lsh0 + 814);
    const auto *lsh0_815 = buffer.data(lsh0 + 815);
    const auto *lsh0_816 = buffer.data(lsh0 + 816);
    const auto *lsh0_817 = buffer.data(lsh0 + 817);
    const auto *lsh0_818 = buffer.data(lsh0 + 818);
    const auto *lsh0_819 = buffer.data(lsh0 + 819);
    const auto *lsh0_820 = buffer.data(lsh0 + 820);
    const auto *lsh0_821 = buffer.data(lsh0 + 821);
    const auto *lsh0_822 = buffer.data(lsh0 + 822);

    const auto *lsh1_756 = buffer.data(lsh1 + 756);
    const auto *lsh1_757 = buffer.data(lsh1 + 757);
    const auto *lsh1_759 = buffer.data(lsh1 + 759);
    const auto *lsh1_761 = buffer.data(lsh1 + 761);
    const auto *lsh1_762 = buffer.data(lsh1 + 762);
    const auto *lsh1_764 = buffer.data(lsh1 + 764);
    const auto *lsh1_765 = buffer.data(lsh1 + 765);
    const auto *lsh1_766 = buffer.data(lsh1 + 766);
    const auto *lsh1_768 = buffer.data(lsh1 + 768);
    const auto *lsh1_769 = buffer.data(lsh1 + 769);
    const auto *lsh1_770 = buffer.data(lsh1 + 770);
    const auto *lsh1_771 = buffer.data(lsh1 + 771);
    const auto *lsh1_772 = buffer.data(lsh1 + 772);
    const auto *lsh1_773 = buffer.data(lsh1 + 773);
    const auto *lsh1_774 = buffer.data(lsh1 + 774);
    const auto *lsh1_775 = buffer.data(lsh1 + 775);
    const auto *lsh1_776 = buffer.data(lsh1 + 776);
    const auto *lsh1_779 = buffer.data(lsh1 + 779);
    const auto *lsh1_781 = buffer.data(lsh1 + 781);
    const auto *lsh1_782 = buffer.data(lsh1 + 782);
    const auto *lsh1_784 = buffer.data(lsh1 + 784);
    const auto *lsh1_785 = buffer.data(lsh1 + 785);
    const auto *lsh1_786 = buffer.data(lsh1 + 786);
    const auto *lsh1_788 = buffer.data(lsh1 + 788);
    const auto *lsh1_789 = buffer.data(lsh1 + 789);
    const auto *lsh1_790 = buffer.data(lsh1 + 790);
    const auto *lsh1_791 = buffer.data(lsh1 + 791);
    const auto *lsh1_793 = buffer.data(lsh1 + 793);
    const auto *lsh1_794 = buffer.data(lsh1 + 794);
    const auto *lsh1_795 = buffer.data(lsh1 + 795);
    const auto *lsh1_796 = buffer.data(lsh1 + 796);
    const auto *lsh1_797 = buffer.data(lsh1 + 797);
    const auto *lsh1_798 = buffer.data(lsh1 + 798);
    const auto *lsh1_799 = buffer.data(lsh1 + 799);
    const auto *lsh1_800 = buffer.data(lsh1 + 800);
    const auto *lsh1_801 = buffer.data(lsh1 + 801);
    const auto *lsh1_802 = buffer.data(lsh1 + 802);
    const auto *lsh1_803 = buffer.data(lsh1 + 803);
    const auto *lsh1_804 = buffer.data(lsh1 + 804);
    const auto *lsh1_805 = buffer.data(lsh1 + 805);
    const auto *lsh1_806 = buffer.data(lsh1 + 806);
    const auto *lsh1_807 = buffer.data(lsh1 + 807);
    const auto *lsh1_808 = buffer.data(lsh1 + 808);
    const auto *lsh1_809 = buffer.data(lsh1 + 809);
    const auto *lsh1_810 = buffer.data(lsh1 + 810);
    const auto *lsh1_811 = buffer.data(lsh1 + 811);
    const auto *lsh1_812 = buffer.data(lsh1 + 812);
    const auto *lsh1_813 = buffer.data(lsh1 + 813);
    const auto *lsh1_814 = buffer.data(lsh1 + 814);
    const auto *lsh1_815 = buffer.data(lsh1 + 815);
    const auto *lsh1_816 = buffer.data(lsh1 + 816);
    const auto *lsh1_817 = buffer.data(lsh1 + 817);
    const auto *lsh1_818 = buffer.data(lsh1 + 818);
    const auto *lsh1_819 = buffer.data(lsh1 + 819);
    const auto *lsh1_820 = buffer.data(lsh1 + 820);
    const auto *lsh1_821 = buffer.data(lsh1 + 821);
    const auto *lsh1_822 = buffer.data(lsh1 + 822);

    const auto *lsi_1007 = buffer.data(lsi + 1007);
    const auto *lsi_1008 = buffer.data(lsi + 1008);
    const auto *lsi_1009 = buffer.data(lsi + 1009);
    const auto *lsi_1011 = buffer.data(lsi + 1011);
    const auto *lsi_1013 = buffer.data(lsi + 1013);
    const auto *lsi_1014 = buffer.data(lsi + 1014);
    const auto *lsi_1016 = buffer.data(lsi + 1016);
    const auto *lsi_1017 = buffer.data(lsi + 1017);
    const auto *lsi_1018 = buffer.data(lsi + 1018);
    const auto *lsi_1020 = buffer.data(lsi + 1020);
    const auto *lsi_1021 = buffer.data(lsi + 1021);
    const auto *lsi_1022 = buffer.data(lsi + 1022);
    const auto *lsi_1023 = buffer.data(lsi + 1023);
    const auto *lsi_1025 = buffer.data(lsi + 1025);
    const auto *lsi_1026 = buffer.data(lsi + 1026);
    const auto *lsi_1027 = buffer.data(lsi + 1027);
    const auto *lsi_1028 = buffer.data(lsi + 1028);
    const auto *lsi_1029 = buffer.data(lsi + 1029);
    const auto *lsi_1030 = buffer.data(lsi + 1030);
    const auto *lsi_1031 = buffer.data(lsi + 1031);
    const auto *lsi_1032 = buffer.data(lsi + 1032);
    const auto *lsi_1033 = buffer.data(lsi + 1033);
    const auto *lsi_1034 = buffer.data(lsi + 1034);
    const auto *lsi_1035 = buffer.data(lsi + 1035);
    const auto *lsi_1038 = buffer.data(lsi + 1038);
    const auto *lsi_1040 = buffer.data(lsi + 1040);
    const auto *lsi_1041 = buffer.data(lsi + 1041);
    const auto *lsi_1043 = buffer.data(lsi + 1043);
    const auto *lsi_1044 = buffer.data(lsi + 1044);
    const auto *lsi_1045 = buffer.data(lsi + 1045);
    const auto *lsi_1047 = buffer.data(lsi + 1047);
    const auto *lsi_1048 = buffer.data(lsi + 1048);
    const auto *lsi_1049 = buffer.data(lsi + 1049);
    const auto *lsi_1050 = buffer.data(lsi + 1050);
    const auto *lsi_1052 = buffer.data(lsi + 1052);
    const auto *lsi_1053 = buffer.data(lsi + 1053);
    const auto *lsi_1054 = buffer.data(lsi + 1054);
    const auto *lsi_1055 = buffer.data(lsi + 1055);
    const auto *lsi_1056 = buffer.data(lsi + 1056);
    const auto *lsi_1057 = buffer.data(lsi + 1057);
    const auto *lsi_1058 = buffer.data(lsi + 1058);
    const auto *lsi_1059 = buffer.data(lsi + 1059);
    const auto *lsi_1060 = buffer.data(lsi + 1060);
    const auto *lsi_1061 = buffer.data(lsi + 1061);
    const auto *lsi_1062 = buffer.data(lsi + 1062);
    const auto *lsi_1063 = buffer.data(lsi + 1063);
    const auto *lsi_1064 = buffer.data(lsi + 1064);
    const auto *lsi_1065 = buffer.data(lsi + 1065);
    const auto *lsi_1066 = buffer.data(lsi + 1066);
    const auto *lsi_1067 = buffer.data(lsi + 1067);
    const auto *lsi_1068 = buffer.data(lsi + 1068);
    const auto *lsi_1069 = buffer.data(lsi + 1069);
    const auto *lsi_1070 = buffer.data(lsi + 1070);
    const auto *lsi_1071 = buffer.data(lsi + 1071);
    const auto *lsi_1072 = buffer.data(lsi + 1072);
    const auto *lsi_1073 = buffer.data(lsi + 1073);
    const auto *lsi_1074 = buffer.data(lsi + 1074);
    const auto *lsi_1075 = buffer.data(lsi + 1075);
    const auto *lsi_1076 = buffer.data(lsi + 1076);
    const auto *lsi_1077 = buffer.data(lsi + 1077);
    const auto *lsi_1078 = buffer.data(lsi + 1078);
    const auto *lsi_1079 = buffer.data(lsi + 1079);
    const auto *lsi_1080 = buffer.data(lsi + 1080);
    const auto *lsi_1081 = buffer.data(lsi + 1081);
    const auto *lsi_1082 = buffer.data(lsi + 1082);
    const auto *lsi_1083 = buffer.data(lsi + 1083);
    const auto *lsi_1084 = buffer.data(lsi + 1084);
    const auto *lsi_1085 = buffer.data(lsi + 1085);
    const auto *lsi_1086 = buffer.data(lsi + 1086);
    const auto *lsi_1087 = buffer.data(lsi + 1087);
    const auto *lsi_1088 = buffer.data(lsi + 1088);
    const auto *lsi_1089 = buffer.data(lsi + 1089);
    const auto *lsi_1090 = buffer.data(lsi + 1090);
    const auto *lsi_1091 = buffer.data(lsi + 1091);
    const auto *lsi_1092 = buffer.data(lsi + 1092);
    const auto *lsi_1093 = buffer.data(lsi + 1093);
    const auto *lsi_1094 = buffer.data(lsi + 1094);
    const auto *lsi_1095 = buffer.data(lsi + 1095);

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pa_x, pc_x, ksk0_1288, ksk0_1289, \
                         ksk0_1290, ksk0_1291, ksk1_1288, ksk1_1289, ksk1_1290, \
                         ksk1_1291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = pa_x[k] * ksk0_1288[k]
                    - f_12 * pc_x[k] * ksk1_1288[k];

        t_1289[k] = pa_x[k] * ksk0_1289[k]
                    - f_12 * pc_x[k] * ksk1_1289[k];

        t_1290[k] = pa_x[k] * ksk0_1290[k]
                    - f_12 * pc_x[k] * ksk1_1290[k];

        t_1291[k] = pa_x[k] * ksk0_1291[k]
                    - f_12 * pc_x[k] * ksk1_1291[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, t_1295, pa_x, pc_x, pc_y, ksk0_1292, \
                         ksk0_1293, ksk0_1295, ksk1_1292, ksk1_1293, ksk1_1295, \
                         lsi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = pa_x[k] * ksk0_1292[k]
                    - f_12 * pc_x[k] * ksk1_1292[k];

        t_1293[k] = pa_x[k] * ksk0_1293[k]
                    - f_12 * pc_x[k] * ksk1_1293[k];

        t_1294[k] = f_3 * pc_y[k] * lsi_1007[k];

        t_1295[k] = pa_x[k] * ksk0_1295[k]
                    - f_12 * pc_x[k] * ksk1_1295[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, t_1299, t_1300, pc_x, pc_z, lsh0_756, \
                         lsh0_757, lsh0_759, lsh1_756, lsh1_757, lsh1_759, lsi_1008, lsi_1009, \
                         lsi_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_1 * lsh0_756[k]
                    - f_2 * lsh1_756[k]
                    + f_3 * pc_x[k] * lsi_1008[k];

        t_1297[k] = f_19 * lsh0_757[k]
                    - f_20 * lsh1_757[k]
                    + f_3 * pc_x[k] * lsi_1009[k];

        t_1298[k] = f_3 * pc_z[k] * lsi_1008[k];

        t_1299[k] = f_10 * lsh0_759[k]
                    - f_11 * lsh1_759[k]
                    + f_3 * pc_x[k] * lsi_1011[k];

        t_1300[k] = f_3 * pc_z[k] * lsi_1009[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, pc_x, pc_z, lsh0_761, lsh0_762, \
                         lsh0_764, lsh1_761, lsh1_762, lsh1_764, lsi_1011, lsi_1013, lsi_1014, \
                         lsi_1016 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_10 * lsh0_761[k]
                    - f_11 * lsh1_761[k]
                    + f_3 * pc_x[k] * lsi_1013[k];

        t_1302[k] = f_8 * lsh0_762[k]
                    - f_9 * lsh1_762[k]
                    + f_3 * pc_x[k] * lsi_1014[k];

        t_1303[k] = f_3 * pc_z[k] * lsi_1011[k];

        t_1304[k] = f_8 * lsh0_764[k]
                    - f_9 * lsh1_764[k]
                    + f_3 * pc_x[k] * lsi_1016[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, t_1308, pc_x, pc_z, lsh0_765, lsh0_766, \
                         lsh0_768, lsh1_765, lsh1_766, lsh1_768, lsi_1014, lsi_1017, lsi_1018, \
                         lsi_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = f_8 * lsh0_765[k]
                    - f_9 * lsh1_765[k]
                    + f_3 * pc_x[k] * lsi_1017[k];

        t_1306[k] = f_6 * lsh0_766[k]
                    - f_7 * lsh1_766[k]
                    + f_3 * pc_x[k] * lsi_1018[k];

        t_1307[k] = f_3 * pc_z[k] * lsi_1014[k];

        t_1308[k] = f_6 * lsh0_768[k]
                    - f_7 * lsh1_768[k]
                    + f_3 * pc_x[k] * lsi_1020[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, t_1312, pc_x, pc_z, lsh0_769, lsh0_770, \
                         lsh0_771, lsh1_769, lsh1_770, lsh1_771, lsi_1018, lsi_1021, lsi_1022, \
                         lsi_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = f_6 * lsh0_769[k]
                    - f_7 * lsh1_769[k]
                    + f_3 * pc_x[k] * lsi_1021[k];

        t_1310[k] = f_6 * lsh0_770[k]
                    - f_7 * lsh1_770[k]
                    + f_3 * pc_x[k] * lsi_1022[k];

        t_1311[k] = f_4 * lsh0_771[k]
                    - f_5 * lsh1_771[k]
                    + f_3 * pc_x[k] * lsi_1023[k];

        t_1312[k] = f_3 * pc_z[k] * lsi_1018[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, pc_x, lsh0_773, lsh0_774, lsh0_775, lsh1_773, \
                         lsh1_774, lsh1_775, lsi_1025, lsi_1026, \
                         lsi_1027 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_4 * lsh0_773[k]
                    - f_5 * lsh1_773[k]
                    + f_3 * pc_x[k] * lsi_1025[k];

        t_1314[k] = f_4 * lsh0_774[k]
                    - f_5 * lsh1_774[k]
                    + f_3 * pc_x[k] * lsi_1026[k];

        t_1315[k] = f_4 * lsh0_775[k]
                    - f_5 * lsh1_775[k]
                    + f_3 * pc_x[k] * lsi_1027[k];
    }

#pragma omp simd aligned(t_1316, t_1317, t_1318, t_1319, t_1320, t_1321, pc_x, lsh0_776, \
                         lsh1_776, lsi_1028, lsi_1029, lsi_1030, lsi_1031, lsi_1032, \
                         lsi_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1316[k] = f_4 * lsh0_776[k]
                    - f_5 * lsh1_776[k]
                    + f_3 * pc_x[k] * lsi_1028[k];

        t_1317[k] = f_3 * pc_x[k] * lsi_1029[k];

        t_1318[k] = f_3 * pc_x[k] * lsi_1030[k];

        t_1319[k] = f_3 * pc_x[k] * lsi_1031[k];

        t_1320[k] = f_3 * pc_x[k] * lsi_1032[k];

        t_1321[k] = f_3 * pc_x[k] * lsi_1033[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, t_1325, t_1326, pc_x, pc_y, pc_z, ksi_805, \
                         lsh0_771, lsh1_771, lsi_1029, lsi_1030, lsi_1034, \
                         lsi_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = f_3 * pc_x[k] * lsi_1034[k];

        t_1323[k] = f_3 * pc_x[k] * lsi_1035[k];

        t_1324[k] = f_0 * ksi_805[k]
                    + f_1 * lsh0_771[k]
                    - f_2 * lsh1_771[k]
                    + f_3 * pc_y[k] * lsi_1029[k];

        t_1325[k] = f_3 * pc_z[k] * lsi_1029[k];

        t_1326[k] = f_4 * lsh0_771[k]
                    - f_5 * lsh1_771[k]
                    + f_3 * pc_z[k] * lsi_1030[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, pc_z, lsh0_772, lsh0_773, lsh0_774, lsh1_772, \
                         lsh1_773, lsh1_774, lsi_1031, lsi_1032, \
                         lsi_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = f_6 * lsh0_772[k]
                    - f_7 * lsh1_772[k]
                    + f_3 * pc_z[k] * lsi_1031[k];

        t_1328[k] = f_8 * lsh0_773[k]
                    - f_9 * lsh1_773[k]
                    + f_3 * pc_z[k] * lsi_1032[k];

        t_1329[k] = f_10 * lsh0_774[k]
                    - f_11 * lsh1_774[k]
                    + f_3 * pc_z[k] * lsi_1033[k];
    }

#pragma omp simd aligned(t_1330, t_1331, t_1332, t_1333, pa_z, pc_y, pc_z, ksk0_1008, \
                         ksk0_1009, ksi_811, ksk1_1008, ksk1_1009, lsh0_776, lsh1_776, \
                         lsi_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1330[k] = f_0 * ksi_811[k]
                    + f_3 * pc_y[k] * lsi_1035[k];

        t_1331[k] = f_1 * lsh0_776[k]
                    - f_2 * lsh1_776[k]
                    + f_3 * pc_z[k] * lsi_1035[k];

        t_1332[k] = pa_z[k] * ksk0_1008[k]
                    - f_12 * pc_z[k] * ksk1_1008[k];

        t_1333[k] = pa_z[k] * ksk0_1009[k]
                    - f_12 * pc_z[k] * ksk1_1009[k];
    }

#pragma omp simd aligned(t_1334, t_1335, t_1336, pa_z, pc_x, pc_z, ksk0_1011, ksk1_1011, \
                         lsh0_779, lsh0_781, lsh1_779, lsh1_781, lsi_1038, \
                         lsi_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1334[k] = f_19 * lsh0_779[k]
                    - f_20 * lsh1_779[k]
                    + f_3 * pc_x[k] * lsi_1038[k];

        t_1335[k] = pa_z[k] * ksk0_1011[k]
                    - f_12 * pc_z[k] * ksk1_1011[k];

        t_1336[k] = f_10 * lsh0_781[k]
                    - f_11 * lsh1_781[k]
                    + f_3 * pc_x[k] * lsi_1040[k];
    }

#pragma omp simd aligned(t_1337, t_1338, t_1339, pa_z, pc_x, pc_z, ksk0_1014, ksk1_1014, \
                         lsh0_782, lsh0_784, lsh1_782, lsh1_784, lsi_1041, \
                         lsi_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1337[k] = f_10 * lsh0_782[k]
                    - f_11 * lsh1_782[k]
                    + f_3 * pc_x[k] * lsi_1041[k];

        t_1338[k] = pa_z[k] * ksk0_1014[k]
                    - f_12 * pc_z[k] * ksk1_1014[k];

        t_1339[k] = f_8 * lsh0_784[k]
                    - f_9 * lsh1_784[k]
                    + f_3 * pc_x[k] * lsi_1043[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, pa_z, pc_x, pc_z, ksk0_1018, ksk1_1018, \
                         lsh0_785, lsh0_786, lsh1_785, lsh1_786, lsi_1044, \
                         lsi_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_8 * lsh0_785[k]
                    - f_9 * lsh1_785[k]
                    + f_3 * pc_x[k] * lsi_1044[k];

        t_1341[k] = f_8 * lsh0_786[k]
                    - f_9 * lsh1_786[k]
                    + f_3 * pc_x[k] * lsi_1045[k];

        t_1342[k] = pa_z[k] * ksk0_1018[k]
                    - f_12 * pc_z[k] * ksk1_1018[k];
    }

#pragma omp simd aligned(t_1343, t_1344, t_1345, pc_x, lsh0_788, lsh0_789, lsh0_790, lsh1_788, \
                         lsh1_789, lsh1_790, lsi_1047, lsi_1048, \
                         lsi_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1343[k] = f_6 * lsh0_788[k]
                    - f_7 * lsh1_788[k]
                    + f_3 * pc_x[k] * lsi_1047[k];

        t_1344[k] = f_6 * lsh0_789[k]
                    - f_7 * lsh1_789[k]
                    + f_3 * pc_x[k] * lsi_1048[k];

        t_1345[k] = f_6 * lsh0_790[k]
                    - f_7 * lsh1_790[k]
                    + f_3 * pc_x[k] * lsi_1049[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pa_z, pc_x, pc_z, ksk0_1023, ksk1_1023, \
                         lsh0_791, lsh0_793, lsh1_791, lsh1_793, lsi_1050, \
                         lsi_1052 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_6 * lsh0_791[k]
                    - f_7 * lsh1_791[k]
                    + f_3 * pc_x[k] * lsi_1050[k];

        t_1347[k] = pa_z[k] * ksk0_1023[k]
                    - f_12 * pc_z[k] * ksk1_1023[k];

        t_1348[k] = f_4 * lsh0_793[k]
                    - f_5 * lsh1_793[k]
                    + f_3 * pc_x[k] * lsi_1052[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pc_x, lsh0_794, lsh0_795, lsh0_796, lsh1_794, \
                         lsh1_795, lsh1_796, lsi_1053, lsi_1054, \
                         lsi_1055 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = f_4 * lsh0_794[k]
                    - f_5 * lsh1_794[k]
                    + f_3 * pc_x[k] * lsi_1053[k];

        t_1350[k] = f_4 * lsh0_795[k]
                    - f_5 * lsh1_795[k]
                    + f_3 * pc_x[k] * lsi_1054[k];

        t_1351[k] = f_4 * lsh0_796[k]
                    - f_5 * lsh1_796[k]
                    + f_3 * pc_x[k] * lsi_1055[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, t_1355, t_1356, t_1357, pc_x, lsh0_797, \
                         lsh1_797, lsi_1056, lsi_1057, lsi_1058, lsi_1059, lsi_1060, \
                         lsi_1061 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_4 * lsh0_797[k]
                    - f_5 * lsh1_797[k]
                    + f_3 * pc_x[k] * lsi_1056[k];

        t_1353[k] = f_3 * pc_x[k] * lsi_1057[k];

        t_1354[k] = f_3 * pc_x[k] * lsi_1058[k];

        t_1355[k] = f_3 * pc_x[k] * lsi_1059[k];

        t_1356[k] = f_3 * pc_x[k] * lsi_1060[k];

        t_1357[k] = f_3 * pc_x[k] * lsi_1061[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, t_1361, pa_z, pc_x, pc_z, ksk0_1036, ksi_805, \
                         ksk1_1036, lsi_1057, lsi_1062, lsi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_3 * pc_x[k] * lsi_1062[k];

        t_1359[k] = f_3 * pc_x[k] * lsi_1063[k];

        t_1360[k] = pa_z[k] * ksk0_1036[k]
                    - f_12 * pc_z[k] * ksk1_1036[k];

        t_1361[k] = f_13 * ksi_805[k]
                    + f_3 * pc_z[k] * lsi_1057[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, pa_z, pc_z, ksk0_1038, ksk0_1039, ksk0_1040, \
                         ksi_806, ksi_807, ksi_808, ksk1_1038, ksk1_1039, \
                         ksk1_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = pa_z[k] * ksk0_1038[k]
                    + f_14 * ksi_806[k]
                    - f_12 * pc_z[k] * ksk1_1038[k];

        t_1363[k] = pa_z[k] * ksk0_1039[k]
                    + f_15 * ksi_807[k]
                    - f_12 * pc_z[k] * ksk1_1039[k];

        t_1364[k] = pa_z[k] * ksk0_1040[k]
                    + f_16 * ksi_808[k]
                    - f_12 * pc_z[k] * ksk1_1040[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, pa_z, pc_y, pc_z, ksk0_1041, ksi_809, \
                         ksi_811, ksi_839, ksk1_1041, lsh0_797, lsh1_797, \
                         lsi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = pa_z[k] * ksk0_1041[k]
                    + f_17 * ksi_809[k]
                    - f_12 * pc_z[k] * ksk1_1041[k];

        t_1366[k] = f_18 * ksi_839[k]
                    + f_3 * pc_y[k] * lsi_1063[k];

        t_1367[k] = f_13 * ksi_811[k]
                    + f_1 * lsh0_797[k]
                    - f_2 * lsh1_797[k]
                    + f_3 * pc_z[k] * lsi_1063[k];
    }

#pragma omp simd aligned(t_1368, t_1369, t_1370, pc_x, lsh0_798, lsh0_799, lsh0_800, lsh1_798, \
                         lsh1_799, lsh1_800, lsi_1064, lsi_1065, \
                         lsi_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1368[k] = f_1 * lsh0_798[k]
                    - f_2 * lsh1_798[k]
                    + f_3 * pc_x[k] * lsi_1064[k];

        t_1369[k] = f_19 * lsh0_799[k]
                    - f_20 * lsh1_799[k]
                    + f_3 * pc_x[k] * lsi_1065[k];

        t_1370[k] = f_19 * lsh0_800[k]
                    - f_20 * lsh1_800[k]
                    + f_3 * pc_x[k] * lsi_1066[k];
    }

#pragma omp simd aligned(t_1371, t_1372, t_1373, pc_x, lsh0_801, lsh0_802, lsh0_803, lsh1_801, \
                         lsh1_802, lsh1_803, lsi_1067, lsi_1068, \
                         lsi_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1371[k] = f_10 * lsh0_801[k]
                    - f_11 * lsh1_801[k]
                    + f_3 * pc_x[k] * lsi_1067[k];

        t_1372[k] = f_10 * lsh0_802[k]
                    - f_11 * lsh1_802[k]
                    + f_3 * pc_x[k] * lsi_1068[k];

        t_1373[k] = f_10 * lsh0_803[k]
                    - f_11 * lsh1_803[k]
                    + f_3 * pc_x[k] * lsi_1069[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, pc_x, lsh0_804, lsh0_805, lsh0_806, lsh1_804, \
                         lsh1_805, lsh1_806, lsi_1070, lsi_1071, \
                         lsi_1072 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_8 * lsh0_804[k]
                    - f_9 * lsh1_804[k]
                    + f_3 * pc_x[k] * lsi_1070[k];

        t_1375[k] = f_8 * lsh0_805[k]
                    - f_9 * lsh1_805[k]
                    + f_3 * pc_x[k] * lsi_1071[k];

        t_1376[k] = f_8 * lsh0_806[k]
                    - f_9 * lsh1_806[k]
                    + f_3 * pc_x[k] * lsi_1072[k];
    }

#pragma omp simd aligned(t_1377, t_1378, t_1379, pc_x, lsh0_807, lsh0_808, lsh0_809, lsh1_807, \
                         lsh1_808, lsh1_809, lsi_1073, lsi_1074, \
                         lsi_1075 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1377[k] = f_8 * lsh0_807[k]
                    - f_9 * lsh1_807[k]
                    + f_3 * pc_x[k] * lsi_1073[k];

        t_1378[k] = f_6 * lsh0_808[k]
                    - f_7 * lsh1_808[k]
                    + f_3 * pc_x[k] * lsi_1074[k];

        t_1379[k] = f_6 * lsh0_809[k]
                    - f_7 * lsh1_809[k]
                    + f_3 * pc_x[k] * lsi_1075[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, pc_x, lsh0_810, lsh0_811, lsh0_812, lsh1_810, \
                         lsh1_811, lsh1_812, lsi_1076, lsi_1077, \
                         lsi_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_6 * lsh0_810[k]
                    - f_7 * lsh1_810[k]
                    + f_3 * pc_x[k] * lsi_1076[k];

        t_1381[k] = f_6 * lsh0_811[k]
                    - f_7 * lsh1_811[k]
                    + f_3 * pc_x[k] * lsi_1077[k];

        t_1382[k] = f_6 * lsh0_812[k]
                    - f_7 * lsh1_812[k]
                    + f_3 * pc_x[k] * lsi_1078[k];
    }

#pragma omp simd aligned(t_1383, t_1384, t_1385, pc_x, lsh0_813, lsh0_814, lsh0_815, lsh1_813, \
                         lsh1_814, lsh1_815, lsi_1079, lsi_1080, \
                         lsi_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1383[k] = f_4 * lsh0_813[k]
                    - f_5 * lsh1_813[k]
                    + f_3 * pc_x[k] * lsi_1079[k];

        t_1384[k] = f_4 * lsh0_814[k]
                    - f_5 * lsh1_814[k]
                    + f_3 * pc_x[k] * lsi_1080[k];

        t_1385[k] = f_4 * lsh0_815[k]
                    - f_5 * lsh1_815[k]
                    + f_3 * pc_x[k] * lsi_1081[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, t_1389, pc_x, lsh0_816, lsh0_817, lsh0_818, \
                         lsh1_816, lsh1_817, lsh1_818, lsi_1082, lsi_1083, lsi_1084, \
                         lsi_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = f_4 * lsh0_816[k]
                    - f_5 * lsh1_816[k]
                    + f_3 * pc_x[k] * lsi_1082[k];

        t_1387[k] = f_4 * lsh0_817[k]
                    - f_5 * lsh1_817[k]
                    + f_3 * pc_x[k] * lsi_1083[k];

        t_1388[k] = f_4 * lsh0_818[k]
                    - f_5 * lsh1_818[k]
                    + f_3 * pc_x[k] * lsi_1084[k];

        t_1389[k] = f_3 * pc_x[k] * lsi_1085[k];
    }

#pragma omp simd aligned(t_1390, t_1391, t_1392, t_1393, t_1394, t_1395, pc_x, lsi_1086, \
                         lsi_1087, lsi_1088, lsi_1089, lsi_1090, \
                         lsi_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1390[k] = f_3 * pc_x[k] * lsi_1086[k];

        t_1391[k] = f_3 * pc_x[k] * lsi_1087[k];

        t_1392[k] = f_3 * pc_x[k] * lsi_1088[k];

        t_1393[k] = f_3 * pc_x[k] * lsi_1089[k];

        t_1394[k] = f_3 * pc_x[k] * lsi_1090[k];

        t_1395[k] = f_3 * pc_x[k] * lsi_1091[k];
    }

#pragma omp simd aligned(t_1396, t_1397, t_1398, pc_y, pc_z, ksi_833, ksi_861, ksi_863, \
                         lsh0_813, lsh0_815, lsh1_813, lsh1_815, lsi_1085, \
                         lsi_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1396[k] = f_21 * ksi_861[k]
                    + f_1 * lsh0_813[k]
                    - f_2 * lsh1_813[k]
                    + f_3 * pc_y[k] * lsi_1085[k];

        t_1397[k] = f_14 * ksi_833[k]
                    + f_3 * pc_z[k] * lsi_1085[k];

        t_1398[k] = f_21 * ksi_863[k]
                    + f_10 * lsh0_815[k]
                    - f_11 * lsh1_815[k]
                    + f_3 * pc_y[k] * lsi_1087[k];
    }

#pragma omp simd aligned(t_1399, t_1400, t_1401, pc_y, ksi_864, ksi_865, ksi_866, lsh0_816, \
                         lsh0_817, lsh0_818, lsh1_816, lsh1_817, lsh1_818, lsi_1088, lsi_1089, \
                         lsi_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1399[k] = f_21 * ksi_864[k]
                    + f_8 * lsh0_816[k]
                    - f_9 * lsh1_816[k]
                    + f_3 * pc_y[k] * lsi_1088[k];

        t_1400[k] = f_21 * ksi_865[k]
                    + f_6 * lsh0_817[k]
                    - f_7 * lsh1_817[k]
                    + f_3 * pc_y[k] * lsi_1089[k];

        t_1401[k] = f_21 * ksi_866[k]
                    + f_4 * lsh0_818[k]
                    - f_5 * lsh1_818[k]
                    + f_3 * pc_y[k] * lsi_1090[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, pc_x, pc_y, pc_z, ksi_839, ksi_867, lsh0_818, \
                         lsh0_819, lsh1_818, lsh1_819, lsi_1091, \
                         lsi_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_21 * ksi_867[k]
                    + f_3 * pc_y[k] * lsi_1091[k];

        t_1403[k] = f_14 * ksi_839[k]
                    + f_1 * lsh0_818[k]
                    - f_2 * lsh1_818[k]
                    + f_3 * pc_z[k] * lsi_1091[k];

        t_1404[k] = f_1 * lsh0_819[k]
                    - f_2 * lsh1_819[k]
                    + f_3 * pc_x[k] * lsi_1092[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, pc_x, lsh0_820, lsh0_821, lsh0_822, lsh1_820, \
                         lsh1_821, lsh1_822, lsi_1093, lsi_1094, \
                         lsi_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = f_19 * lsh0_820[k]
                    - f_20 * lsh1_820[k]
                    + f_3 * pc_x[k] * lsi_1093[k];

        t_1406[k] = f_19 * lsh0_821[k]
                    - f_20 * lsh1_821[k]
                    + f_3 * pc_x[k] * lsi_1094[k];

        t_1407[k] = f_10 * lsh0_822[k]
                    - f_11 * lsh1_822[k]
                    + f_3 * pc_x[k] * lsi_1095[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t ksi, const size_t lsh0,
                                                           const size_t lsh1, const size_t lsi,
                                                           const size_t ncols,
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
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);

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
    auto *t_1521 = buffer.data(target + 1521);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi_861 = buffer.data(ksi + 861);
    const auto *ksi_867 = buffer.data(ksi + 867);
    const auto *ksi_889 = buffer.data(ksi + 889);
    const auto *ksi_891 = buffer.data(ksi + 891);
    const auto *ksi_892 = buffer.data(ksi + 892);
    const auto *ksi_893 = buffer.data(ksi + 893);
    const auto *ksi_894 = buffer.data(ksi + 894);
    const auto *ksi_895 = buffer.data(ksi + 895);
    const auto *ksi_917 = buffer.data(ksi + 917);
    const auto *ksi_919 = buffer.data(ksi + 919);
    const auto *ksi_920 = buffer.data(ksi + 920);
    const auto *ksi_921 = buffer.data(ksi + 921);
    const auto *ksi_922 = buffer.data(ksi + 922);
    const auto *ksi_923 = buffer.data(ksi + 923);
    const auto *ksi_945 = buffer.data(ksi + 945);
    const auto *ksi_947 = buffer.data(ksi + 947);
    const auto *ksi_948 = buffer.data(ksi + 948);
    const auto *ksi_949 = buffer.data(ksi + 949);
    const auto *ksi_950 = buffer.data(ksi + 950);
    const auto *ksi_951 = buffer.data(ksi + 951);

    const auto *lsh0_823 = buffer.data(lsh0 + 823);
    const auto *lsh0_824 = buffer.data(lsh0 + 824);
    const auto *lsh0_825 = buffer.data(lsh0 + 825);
    const auto *lsh0_826 = buffer.data(lsh0 + 826);
    const auto *lsh0_827 = buffer.data(lsh0 + 827);
    const auto *lsh0_828 = buffer.data(lsh0 + 828);
    const auto *lsh0_829 = buffer.data(lsh0 + 829);
    const auto *lsh0_830 = buffer.data(lsh0 + 830);
    const auto *lsh0_831 = buffer.data(lsh0 + 831);
    const auto *lsh0_832 = buffer.data(lsh0 + 832);
    const auto *lsh0_833 = buffer.data(lsh0 + 833);
    const auto *lsh0_834 = buffer.data(lsh0 + 834);
    const auto *lsh0_835 = buffer.data(lsh0 + 835);
    const auto *lsh0_836 = buffer.data(lsh0 + 836);
    const auto *lsh0_837 = buffer.data(lsh0 + 837);
    const auto *lsh0_838 = buffer.data(lsh0 + 838);
    const auto *lsh0_839 = buffer.data(lsh0 + 839);
    const auto *lsh0_840 = buffer.data(lsh0 + 840);
    const auto *lsh0_841 = buffer.data(lsh0 + 841);
    const auto *lsh0_842 = buffer.data(lsh0 + 842);
    const auto *lsh0_843 = buffer.data(lsh0 + 843);
    const auto *lsh0_844 = buffer.data(lsh0 + 844);
    const auto *lsh0_845 = buffer.data(lsh0 + 845);
    const auto *lsh0_846 = buffer.data(lsh0 + 846);
    const auto *lsh0_847 = buffer.data(lsh0 + 847);
    const auto *lsh0_848 = buffer.data(lsh0 + 848);
    const auto *lsh0_849 = buffer.data(lsh0 + 849);
    const auto *lsh0_850 = buffer.data(lsh0 + 850);
    const auto *lsh0_851 = buffer.data(lsh0 + 851);
    const auto *lsh0_852 = buffer.data(lsh0 + 852);
    const auto *lsh0_853 = buffer.data(lsh0 + 853);
    const auto *lsh0_854 = buffer.data(lsh0 + 854);
    const auto *lsh0_855 = buffer.data(lsh0 + 855);
    const auto *lsh0_856 = buffer.data(lsh0 + 856);
    const auto *lsh0_857 = buffer.data(lsh0 + 857);
    const auto *lsh0_858 = buffer.data(lsh0 + 858);
    const auto *lsh0_859 = buffer.data(lsh0 + 859);
    const auto *lsh0_860 = buffer.data(lsh0 + 860);
    const auto *lsh0_861 = buffer.data(lsh0 + 861);
    const auto *lsh0_862 = buffer.data(lsh0 + 862);
    const auto *lsh0_863 = buffer.data(lsh0 + 863);
    const auto *lsh0_864 = buffer.data(lsh0 + 864);
    const auto *lsh0_865 = buffer.data(lsh0 + 865);
    const auto *lsh0_866 = buffer.data(lsh0 + 866);
    const auto *lsh0_867 = buffer.data(lsh0 + 867);
    const auto *lsh0_868 = buffer.data(lsh0 + 868);
    const auto *lsh0_869 = buffer.data(lsh0 + 869);
    const auto *lsh0_870 = buffer.data(lsh0 + 870);
    const auto *lsh0_871 = buffer.data(lsh0 + 871);
    const auto *lsh0_872 = buffer.data(lsh0 + 872);
    const auto *lsh0_873 = buffer.data(lsh0 + 873);
    const auto *lsh0_874 = buffer.data(lsh0 + 874);
    const auto *lsh0_875 = buffer.data(lsh0 + 875);
    const auto *lsh0_876 = buffer.data(lsh0 + 876);
    const auto *lsh0_877 = buffer.data(lsh0 + 877);
    const auto *lsh0_878 = buffer.data(lsh0 + 878);
    const auto *lsh0_879 = buffer.data(lsh0 + 879);
    const auto *lsh0_880 = buffer.data(lsh0 + 880);
    const auto *lsh0_881 = buffer.data(lsh0 + 881);
    const auto *lsh0_882 = buffer.data(lsh0 + 882);
    const auto *lsh0_883 = buffer.data(lsh0 + 883);
    const auto *lsh0_884 = buffer.data(lsh0 + 884);
    const auto *lsh0_885 = buffer.data(lsh0 + 885);
    const auto *lsh0_886 = buffer.data(lsh0 + 886);
    const auto *lsh0_887 = buffer.data(lsh0 + 887);
    const auto *lsh0_888 = buffer.data(lsh0 + 888);
    const auto *lsh0_889 = buffer.data(lsh0 + 889);
    const auto *lsh0_890 = buffer.data(lsh0 + 890);
    const auto *lsh0_891 = buffer.data(lsh0 + 891);

    const auto *lsh1_823 = buffer.data(lsh1 + 823);
    const auto *lsh1_824 = buffer.data(lsh1 + 824);
    const auto *lsh1_825 = buffer.data(lsh1 + 825);
    const auto *lsh1_826 = buffer.data(lsh1 + 826);
    const auto *lsh1_827 = buffer.data(lsh1 + 827);
    const auto *lsh1_828 = buffer.data(lsh1 + 828);
    const auto *lsh1_829 = buffer.data(lsh1 + 829);
    const auto *lsh1_830 = buffer.data(lsh1 + 830);
    const auto *lsh1_831 = buffer.data(lsh1 + 831);
    const auto *lsh1_832 = buffer.data(lsh1 + 832);
    const auto *lsh1_833 = buffer.data(lsh1 + 833);
    const auto *lsh1_834 = buffer.data(lsh1 + 834);
    const auto *lsh1_835 = buffer.data(lsh1 + 835);
    const auto *lsh1_836 = buffer.data(lsh1 + 836);
    const auto *lsh1_837 = buffer.data(lsh1 + 837);
    const auto *lsh1_838 = buffer.data(lsh1 + 838);
    const auto *lsh1_839 = buffer.data(lsh1 + 839);
    const auto *lsh1_840 = buffer.data(lsh1 + 840);
    const auto *lsh1_841 = buffer.data(lsh1 + 841);
    const auto *lsh1_842 = buffer.data(lsh1 + 842);
    const auto *lsh1_843 = buffer.data(lsh1 + 843);
    const auto *lsh1_844 = buffer.data(lsh1 + 844);
    const auto *lsh1_845 = buffer.data(lsh1 + 845);
    const auto *lsh1_846 = buffer.data(lsh1 + 846);
    const auto *lsh1_847 = buffer.data(lsh1 + 847);
    const auto *lsh1_848 = buffer.data(lsh1 + 848);
    const auto *lsh1_849 = buffer.data(lsh1 + 849);
    const auto *lsh1_850 = buffer.data(lsh1 + 850);
    const auto *lsh1_851 = buffer.data(lsh1 + 851);
    const auto *lsh1_852 = buffer.data(lsh1 + 852);
    const auto *lsh1_853 = buffer.data(lsh1 + 853);
    const auto *lsh1_854 = buffer.data(lsh1 + 854);
    const auto *lsh1_855 = buffer.data(lsh1 + 855);
    const auto *lsh1_856 = buffer.data(lsh1 + 856);
    const auto *lsh1_857 = buffer.data(lsh1 + 857);
    const auto *lsh1_858 = buffer.data(lsh1 + 858);
    const auto *lsh1_859 = buffer.data(lsh1 + 859);
    const auto *lsh1_860 = buffer.data(lsh1 + 860);
    const auto *lsh1_861 = buffer.data(lsh1 + 861);
    const auto *lsh1_862 = buffer.data(lsh1 + 862);
    const auto *lsh1_863 = buffer.data(lsh1 + 863);
    const auto *lsh1_864 = buffer.data(lsh1 + 864);
    const auto *lsh1_865 = buffer.data(lsh1 + 865);
    const auto *lsh1_866 = buffer.data(lsh1 + 866);
    const auto *lsh1_867 = buffer.data(lsh1 + 867);
    const auto *lsh1_868 = buffer.data(lsh1 + 868);
    const auto *lsh1_869 = buffer.data(lsh1 + 869);
    const auto *lsh1_870 = buffer.data(lsh1 + 870);
    const auto *lsh1_871 = buffer.data(lsh1 + 871);
    const auto *lsh1_872 = buffer.data(lsh1 + 872);
    const auto *lsh1_873 = buffer.data(lsh1 + 873);
    const auto *lsh1_874 = buffer.data(lsh1 + 874);
    const auto *lsh1_875 = buffer.data(lsh1 + 875);
    const auto *lsh1_876 = buffer.data(lsh1 + 876);
    const auto *lsh1_877 = buffer.data(lsh1 + 877);
    const auto *lsh1_878 = buffer.data(lsh1 + 878);
    const auto *lsh1_879 = buffer.data(lsh1 + 879);
    const auto *lsh1_880 = buffer.data(lsh1 + 880);
    const auto *lsh1_881 = buffer.data(lsh1 + 881);
    const auto *lsh1_882 = buffer.data(lsh1 + 882);
    const auto *lsh1_883 = buffer.data(lsh1 + 883);
    const auto *lsh1_884 = buffer.data(lsh1 + 884);
    const auto *lsh1_885 = buffer.data(lsh1 + 885);
    const auto *lsh1_886 = buffer.data(lsh1 + 886);
    const auto *lsh1_887 = buffer.data(lsh1 + 887);
    const auto *lsh1_888 = buffer.data(lsh1 + 888);
    const auto *lsh1_889 = buffer.data(lsh1 + 889);
    const auto *lsh1_890 = buffer.data(lsh1 + 890);
    const auto *lsh1_891 = buffer.data(lsh1 + 891);

    const auto *lsi_1096 = buffer.data(lsi + 1096);
    const auto *lsi_1097 = buffer.data(lsi + 1097);
    const auto *lsi_1098 = buffer.data(lsi + 1098);
    const auto *lsi_1099 = buffer.data(lsi + 1099);
    const auto *lsi_1100 = buffer.data(lsi + 1100);
    const auto *lsi_1101 = buffer.data(lsi + 1101);
    const auto *lsi_1102 = buffer.data(lsi + 1102);
    const auto *lsi_1103 = buffer.data(lsi + 1103);
    const auto *lsi_1104 = buffer.data(lsi + 1104);
    const auto *lsi_1105 = buffer.data(lsi + 1105);
    const auto *lsi_1106 = buffer.data(lsi + 1106);
    const auto *lsi_1107 = buffer.data(lsi + 1107);
    const auto *lsi_1108 = buffer.data(lsi + 1108);
    const auto *lsi_1109 = buffer.data(lsi + 1109);
    const auto *lsi_1110 = buffer.data(lsi + 1110);
    const auto *lsi_1111 = buffer.data(lsi + 1111);
    const auto *lsi_1112 = buffer.data(lsi + 1112);
    const auto *lsi_1113 = buffer.data(lsi + 1113);
    const auto *lsi_1114 = buffer.data(lsi + 1114);
    const auto *lsi_1115 = buffer.data(lsi + 1115);
    const auto *lsi_1116 = buffer.data(lsi + 1116);
    const auto *lsi_1117 = buffer.data(lsi + 1117);
    const auto *lsi_1118 = buffer.data(lsi + 1118);
    const auto *lsi_1119 = buffer.data(lsi + 1119);
    const auto *lsi_1120 = buffer.data(lsi + 1120);
    const auto *lsi_1121 = buffer.data(lsi + 1121);
    const auto *lsi_1122 = buffer.data(lsi + 1122);
    const auto *lsi_1123 = buffer.data(lsi + 1123);
    const auto *lsi_1124 = buffer.data(lsi + 1124);
    const auto *lsi_1125 = buffer.data(lsi + 1125);
    const auto *lsi_1126 = buffer.data(lsi + 1126);
    const auto *lsi_1127 = buffer.data(lsi + 1127);
    const auto *lsi_1128 = buffer.data(lsi + 1128);
    const auto *lsi_1129 = buffer.data(lsi + 1129);
    const auto *lsi_1130 = buffer.data(lsi + 1130);
    const auto *lsi_1131 = buffer.data(lsi + 1131);
    const auto *lsi_1132 = buffer.data(lsi + 1132);
    const auto *lsi_1133 = buffer.data(lsi + 1133);
    const auto *lsi_1134 = buffer.data(lsi + 1134);
    const auto *lsi_1135 = buffer.data(lsi + 1135);
    const auto *lsi_1136 = buffer.data(lsi + 1136);
    const auto *lsi_1137 = buffer.data(lsi + 1137);
    const auto *lsi_1138 = buffer.data(lsi + 1138);
    const auto *lsi_1139 = buffer.data(lsi + 1139);
    const auto *lsi_1140 = buffer.data(lsi + 1140);
    const auto *lsi_1141 = buffer.data(lsi + 1141);
    const auto *lsi_1142 = buffer.data(lsi + 1142);
    const auto *lsi_1143 = buffer.data(lsi + 1143);
    const auto *lsi_1144 = buffer.data(lsi + 1144);
    const auto *lsi_1145 = buffer.data(lsi + 1145);
    const auto *lsi_1146 = buffer.data(lsi + 1146);
    const auto *lsi_1147 = buffer.data(lsi + 1147);
    const auto *lsi_1148 = buffer.data(lsi + 1148);
    const auto *lsi_1149 = buffer.data(lsi + 1149);
    const auto *lsi_1150 = buffer.data(lsi + 1150);
    const auto *lsi_1151 = buffer.data(lsi + 1151);
    const auto *lsi_1152 = buffer.data(lsi + 1152);
    const auto *lsi_1153 = buffer.data(lsi + 1153);
    const auto *lsi_1154 = buffer.data(lsi + 1154);
    const auto *lsi_1155 = buffer.data(lsi + 1155);
    const auto *lsi_1156 = buffer.data(lsi + 1156);
    const auto *lsi_1157 = buffer.data(lsi + 1157);
    const auto *lsi_1158 = buffer.data(lsi + 1158);
    const auto *lsi_1159 = buffer.data(lsi + 1159);
    const auto *lsi_1160 = buffer.data(lsi + 1160);
    const auto *lsi_1161 = buffer.data(lsi + 1161);
    const auto *lsi_1162 = buffer.data(lsi + 1162);
    const auto *lsi_1163 = buffer.data(lsi + 1163);
    const auto *lsi_1164 = buffer.data(lsi + 1164);
    const auto *lsi_1165 = buffer.data(lsi + 1165);
    const auto *lsi_1166 = buffer.data(lsi + 1166);
    const auto *lsi_1167 = buffer.data(lsi + 1167);
    const auto *lsi_1168 = buffer.data(lsi + 1168);
    const auto *lsi_1169 = buffer.data(lsi + 1169);
    const auto *lsi_1170 = buffer.data(lsi + 1170);
    const auto *lsi_1171 = buffer.data(lsi + 1171);
    const auto *lsi_1172 = buffer.data(lsi + 1172);
    const auto *lsi_1173 = buffer.data(lsi + 1173);
    const auto *lsi_1174 = buffer.data(lsi + 1174);
    const auto *lsi_1175 = buffer.data(lsi + 1175);
    const auto *lsi_1176 = buffer.data(lsi + 1176);
    const auto *lsi_1177 = buffer.data(lsi + 1177);
    const auto *lsi_1178 = buffer.data(lsi + 1178);
    const auto *lsi_1179 = buffer.data(lsi + 1179);
    const auto *lsi_1180 = buffer.data(lsi + 1180);
    const auto *lsi_1181 = buffer.data(lsi + 1181);
    const auto *lsi_1182 = buffer.data(lsi + 1182);
    const auto *lsi_1183 = buffer.data(lsi + 1183);
    const auto *lsi_1184 = buffer.data(lsi + 1184);
    const auto *lsi_1185 = buffer.data(lsi + 1185);

#pragma omp simd aligned(t_1408, t_1409, t_1410, pc_x, lsh0_823, lsh0_824, lsh0_825, lsh1_823, \
                         lsh1_824, lsh1_825, lsi_1096, lsi_1097, \
                         lsi_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1408[k] = f_10 * lsh0_823[k]
                    - f_11 * lsh1_823[k]
                    + f_3 * pc_x[k] * lsi_1096[k];

        t_1409[k] = f_10 * lsh0_824[k]
                    - f_11 * lsh1_824[k]
                    + f_3 * pc_x[k] * lsi_1097[k];

        t_1410[k] = f_8 * lsh0_825[k]
                    - f_9 * lsh1_825[k]
                    + f_3 * pc_x[k] * lsi_1098[k];
    }

#pragma omp simd aligned(t_1411, t_1412, t_1413, pc_x, lsh0_826, lsh0_827, lsh0_828, lsh1_826, \
                         lsh1_827, lsh1_828, lsi_1099, lsi_1100, \
                         lsi_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1411[k] = f_8 * lsh0_826[k]
                    - f_9 * lsh1_826[k]
                    + f_3 * pc_x[k] * lsi_1099[k];

        t_1412[k] = f_8 * lsh0_827[k]
                    - f_9 * lsh1_827[k]
                    + f_3 * pc_x[k] * lsi_1100[k];

        t_1413[k] = f_8 * lsh0_828[k]
                    - f_9 * lsh1_828[k]
                    + f_3 * pc_x[k] * lsi_1101[k];
    }

#pragma omp simd aligned(t_1414, t_1415, t_1416, pc_x, lsh0_829, lsh0_830, lsh0_831, lsh1_829, \
                         lsh1_830, lsh1_831, lsi_1102, lsi_1103, \
                         lsi_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1414[k] = f_6 * lsh0_829[k]
                    - f_7 * lsh1_829[k]
                    + f_3 * pc_x[k] * lsi_1102[k];

        t_1415[k] = f_6 * lsh0_830[k]
                    - f_7 * lsh1_830[k]
                    + f_3 * pc_x[k] * lsi_1103[k];

        t_1416[k] = f_6 * lsh0_831[k]
                    - f_7 * lsh1_831[k]
                    + f_3 * pc_x[k] * lsi_1104[k];
    }

#pragma omp simd aligned(t_1417, t_1418, t_1419, pc_x, lsh0_832, lsh0_833, lsh0_834, lsh1_832, \
                         lsh1_833, lsh1_834, lsi_1105, lsi_1106, \
                         lsi_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1417[k] = f_6 * lsh0_832[k]
                    - f_7 * lsh1_832[k]
                    + f_3 * pc_x[k] * lsi_1105[k];

        t_1418[k] = f_6 * lsh0_833[k]
                    - f_7 * lsh1_833[k]
                    + f_3 * pc_x[k] * lsi_1106[k];

        t_1419[k] = f_4 * lsh0_834[k]
                    - f_5 * lsh1_834[k]
                    + f_3 * pc_x[k] * lsi_1107[k];
    }

#pragma omp simd aligned(t_1420, t_1421, t_1422, pc_x, lsh0_835, lsh0_836, lsh0_837, lsh1_835, \
                         lsh1_836, lsh1_837, lsi_1108, lsi_1109, \
                         lsi_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1420[k] = f_4 * lsh0_835[k]
                    - f_5 * lsh1_835[k]
                    + f_3 * pc_x[k] * lsi_1108[k];

        t_1421[k] = f_4 * lsh0_836[k]
                    - f_5 * lsh1_836[k]
                    + f_3 * pc_x[k] * lsi_1109[k];

        t_1422[k] = f_4 * lsh0_837[k]
                    - f_5 * lsh1_837[k]
                    + f_3 * pc_x[k] * lsi_1110[k];
    }

#pragma omp simd aligned(t_1423, t_1424, t_1425, t_1426, t_1427, pc_x, lsh0_838, lsh0_839, \
                         lsh1_838, lsh1_839, lsi_1111, lsi_1112, lsi_1113, lsi_1114, \
                         lsi_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1423[k] = f_4 * lsh0_838[k]
                    - f_5 * lsh1_838[k]
                    + f_3 * pc_x[k] * lsi_1111[k];

        t_1424[k] = f_4 * lsh0_839[k]
                    - f_5 * lsh1_839[k]
                    + f_3 * pc_x[k] * lsi_1112[k];

        t_1425[k] = f_3 * pc_x[k] * lsi_1113[k];

        t_1426[k] = f_3 * pc_x[k] * lsi_1114[k];

        t_1427[k] = f_3 * pc_x[k] * lsi_1115[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, t_1431, t_1432, pc_x, pc_y, ksi_889, \
                         lsh0_834, lsh1_834, lsi_1113, lsi_1116, lsi_1117, lsi_1118, \
                         lsi_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_3 * pc_x[k] * lsi_1116[k];

        t_1429[k] = f_3 * pc_x[k] * lsi_1117[k];

        t_1430[k] = f_3 * pc_x[k] * lsi_1118[k];

        t_1431[k] = f_3 * pc_x[k] * lsi_1119[k];

        t_1432[k] = f_17 * ksi_889[k]
                    + f_1 * lsh0_834[k]
                    - f_2 * lsh1_834[k]
                    + f_3 * pc_y[k] * lsi_1113[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pc_y, pc_z, ksi_861, ksi_891, ksi_892, \
                         lsh0_836, lsh0_837, lsh1_836, lsh1_837, lsi_1113, lsi_1115, \
                         lsi_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_15 * ksi_861[k]
                    + f_3 * pc_z[k] * lsi_1113[k];

        t_1434[k] = f_17 * ksi_891[k]
                    + f_10 * lsh0_836[k]
                    - f_11 * lsh1_836[k]
                    + f_3 * pc_y[k] * lsi_1115[k];

        t_1435[k] = f_17 * ksi_892[k]
                    + f_8 * lsh0_837[k]
                    - f_9 * lsh1_837[k]
                    + f_3 * pc_y[k] * lsi_1116[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, pc_y, ksi_893, ksi_894, ksi_895, lsh0_838, \
                         lsh0_839, lsh1_838, lsh1_839, lsi_1117, lsi_1118, \
                         lsi_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_17 * ksi_893[k]
                    + f_6 * lsh0_838[k]
                    - f_7 * lsh1_838[k]
                    + f_3 * pc_y[k] * lsi_1117[k];

        t_1437[k] = f_17 * ksi_894[k]
                    + f_4 * lsh0_839[k]
                    - f_5 * lsh1_839[k]
                    + f_3 * pc_y[k] * lsi_1118[k];

        t_1438[k] = f_17 * ksi_895[k]
                    + f_3 * pc_y[k] * lsi_1119[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, pc_x, pc_z, ksi_867, lsh0_839, lsh0_840, \
                         lsh0_841, lsh1_839, lsh1_840, lsh1_841, lsi_1119, lsi_1120, \
                         lsi_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = f_15 * ksi_867[k]
                    + f_1 * lsh0_839[k]
                    - f_2 * lsh1_839[k]
                    + f_3 * pc_z[k] * lsi_1119[k];

        t_1440[k] = f_1 * lsh0_840[k]
                    - f_2 * lsh1_840[k]
                    + f_3 * pc_x[k] * lsi_1120[k];

        t_1441[k] = f_19 * lsh0_841[k]
                    - f_20 * lsh1_841[k]
                    + f_3 * pc_x[k] * lsi_1121[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, pc_x, lsh0_842, lsh0_843, lsh0_844, lsh1_842, \
                         lsh1_843, lsh1_844, lsi_1122, lsi_1123, \
                         lsi_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_19 * lsh0_842[k]
                    - f_20 * lsh1_842[k]
                    + f_3 * pc_x[k] * lsi_1122[k];

        t_1443[k] = f_10 * lsh0_843[k]
                    - f_11 * lsh1_843[k]
                    + f_3 * pc_x[k] * lsi_1123[k];

        t_1444[k] = f_10 * lsh0_844[k]
                    - f_11 * lsh1_844[k]
                    + f_3 * pc_x[k] * lsi_1124[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pc_x, lsh0_845, lsh0_846, lsh0_847, lsh1_845, \
                         lsh1_846, lsh1_847, lsi_1125, lsi_1126, \
                         lsi_1127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_10 * lsh0_845[k]
                    - f_11 * lsh1_845[k]
                    + f_3 * pc_x[k] * lsi_1125[k];

        t_1446[k] = f_8 * lsh0_846[k]
                    - f_9 * lsh1_846[k]
                    + f_3 * pc_x[k] * lsi_1126[k];

        t_1447[k] = f_8 * lsh0_847[k]
                    - f_9 * lsh1_847[k]
                    + f_3 * pc_x[k] * lsi_1127[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, lsh0_848, lsh0_849, lsh0_850, lsh1_848, \
                         lsh1_849, lsh1_850, lsi_1128, lsi_1129, \
                         lsi_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_8 * lsh0_848[k]
                    - f_9 * lsh1_848[k]
                    + f_3 * pc_x[k] * lsi_1128[k];

        t_1449[k] = f_8 * lsh0_849[k]
                    - f_9 * lsh1_849[k]
                    + f_3 * pc_x[k] * lsi_1129[k];

        t_1450[k] = f_6 * lsh0_850[k]
                    - f_7 * lsh1_850[k]
                    + f_3 * pc_x[k] * lsi_1130[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_x, lsh0_851, lsh0_852, lsh0_853, lsh1_851, \
                         lsh1_852, lsh1_853, lsi_1131, lsi_1132, \
                         lsi_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_6 * lsh0_851[k]
                    - f_7 * lsh1_851[k]
                    + f_3 * pc_x[k] * lsi_1131[k];

        t_1452[k] = f_6 * lsh0_852[k]
                    - f_7 * lsh1_852[k]
                    + f_3 * pc_x[k] * lsi_1132[k];

        t_1453[k] = f_6 * lsh0_853[k]
                    - f_7 * lsh1_853[k]
                    + f_3 * pc_x[k] * lsi_1133[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, lsh0_854, lsh0_855, lsh0_856, lsh1_854, \
                         lsh1_855, lsh1_856, lsi_1134, lsi_1135, \
                         lsi_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_6 * lsh0_854[k]
                    - f_7 * lsh1_854[k]
                    + f_3 * pc_x[k] * lsi_1134[k];

        t_1455[k] = f_4 * lsh0_855[k]
                    - f_5 * lsh1_855[k]
                    + f_3 * pc_x[k] * lsi_1135[k];

        t_1456[k] = f_4 * lsh0_856[k]
                    - f_5 * lsh1_856[k]
                    + f_3 * pc_x[k] * lsi_1136[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pc_x, lsh0_857, lsh0_858, lsh0_859, lsh1_857, \
                         lsh1_858, lsh1_859, lsi_1137, lsi_1138, \
                         lsi_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_4 * lsh0_857[k]
                    - f_5 * lsh1_857[k]
                    + f_3 * pc_x[k] * lsi_1137[k];

        t_1458[k] = f_4 * lsh0_858[k]
                    - f_5 * lsh1_858[k]
                    + f_3 * pc_x[k] * lsi_1138[k];

        t_1459[k] = f_4 * lsh0_859[k]
                    - f_5 * lsh1_859[k]
                    + f_3 * pc_x[k] * lsi_1139[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, t_1463, t_1464, t_1465, pc_x, lsh0_860, \
                         lsh1_860, lsi_1140, lsi_1141, lsi_1142, lsi_1143, lsi_1144, \
                         lsi_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_4 * lsh0_860[k]
                    - f_5 * lsh1_860[k]
                    + f_3 * pc_x[k] * lsi_1140[k];

        t_1461[k] = f_3 * pc_x[k] * lsi_1141[k];

        t_1462[k] = f_3 * pc_x[k] * lsi_1142[k];

        t_1463[k] = f_3 * pc_x[k] * lsi_1143[k];

        t_1464[k] = f_3 * pc_x[k] * lsi_1144[k];

        t_1465[k] = f_3 * pc_x[k] * lsi_1145[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, t_1469, pc_x, pc_y, pc_z, ksi_889, ksi_917, \
                         lsh0_855, lsh1_855, lsi_1141, lsi_1146, \
                         lsi_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_3 * pc_x[k] * lsi_1146[k];

        t_1467[k] = f_3 * pc_x[k] * lsi_1147[k];

        t_1468[k] = f_16 * ksi_917[k]
                    + f_1 * lsh0_855[k]
                    - f_2 * lsh1_855[k]
                    + f_3 * pc_y[k] * lsi_1141[k];

        t_1469[k] = f_16 * ksi_889[k]
                    + f_3 * pc_z[k] * lsi_1141[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, pc_y, ksi_919, ksi_920, ksi_921, lsh0_857, \
                         lsh0_858, lsh0_859, lsh1_857, lsh1_858, lsh1_859, lsi_1143, lsi_1144, \
                         lsi_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = f_16 * ksi_919[k]
                    + f_10 * lsh0_857[k]
                    - f_11 * lsh1_857[k]
                    + f_3 * pc_y[k] * lsi_1143[k];

        t_1471[k] = f_16 * ksi_920[k]
                    + f_8 * lsh0_858[k]
                    - f_9 * lsh1_858[k]
                    + f_3 * pc_y[k] * lsi_1144[k];

        t_1472[k] = f_16 * ksi_921[k]
                    + f_6 * lsh0_859[k]
                    - f_7 * lsh1_859[k]
                    + f_3 * pc_y[k] * lsi_1145[k];
    }

#pragma omp simd aligned(t_1473, t_1474, t_1475, pc_y, pc_z, ksi_895, ksi_922, ksi_923, \
                         lsh0_860, lsh1_860, lsi_1146, lsi_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1473[k] = f_16 * ksi_922[k]
                    + f_4 * lsh0_860[k]
                    - f_5 * lsh1_860[k]
                    + f_3 * pc_y[k] * lsi_1146[k];

        t_1474[k] = f_16 * ksi_923[k]
                    + f_3 * pc_y[k] * lsi_1147[k];

        t_1475[k] = f_16 * ksi_895[k]
                    + f_1 * lsh0_860[k]
                    - f_2 * lsh1_860[k]
                    + f_3 * pc_z[k] * lsi_1147[k];
    }

#pragma omp simd aligned(t_1476, t_1477, t_1478, pc_x, lsh0_861, lsh0_862, lsh0_863, lsh1_861, \
                         lsh1_862, lsh1_863, lsi_1148, lsi_1149, \
                         lsi_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1476[k] = f_1 * lsh0_861[k]
                    - f_2 * lsh1_861[k]
                    + f_3 * pc_x[k] * lsi_1148[k];

        t_1477[k] = f_19 * lsh0_862[k]
                    - f_20 * lsh1_862[k]
                    + f_3 * pc_x[k] * lsi_1149[k];

        t_1478[k] = f_19 * lsh0_863[k]
                    - f_20 * lsh1_863[k]
                    + f_3 * pc_x[k] * lsi_1150[k];
    }

#pragma omp simd aligned(t_1479, t_1480, t_1481, pc_x, lsh0_864, lsh0_865, lsh0_866, lsh1_864, \
                         lsh1_865, lsh1_866, lsi_1151, lsi_1152, \
                         lsi_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1479[k] = f_10 * lsh0_864[k]
                    - f_11 * lsh1_864[k]
                    + f_3 * pc_x[k] * lsi_1151[k];

        t_1480[k] = f_10 * lsh0_865[k]
                    - f_11 * lsh1_865[k]
                    + f_3 * pc_x[k] * lsi_1152[k];

        t_1481[k] = f_10 * lsh0_866[k]
                    - f_11 * lsh1_866[k]
                    + f_3 * pc_x[k] * lsi_1153[k];
    }

#pragma omp simd aligned(t_1482, t_1483, t_1484, pc_x, lsh0_867, lsh0_868, lsh0_869, lsh1_867, \
                         lsh1_868, lsh1_869, lsi_1154, lsi_1155, \
                         lsi_1156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1482[k] = f_8 * lsh0_867[k]
                    - f_9 * lsh1_867[k]
                    + f_3 * pc_x[k] * lsi_1154[k];

        t_1483[k] = f_8 * lsh0_868[k]
                    - f_9 * lsh1_868[k]
                    + f_3 * pc_x[k] * lsi_1155[k];

        t_1484[k] = f_8 * lsh0_869[k]
                    - f_9 * lsh1_869[k]
                    + f_3 * pc_x[k] * lsi_1156[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, pc_x, lsh0_870, lsh0_871, lsh0_872, lsh1_870, \
                         lsh1_871, lsh1_872, lsi_1157, lsi_1158, \
                         lsi_1159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = f_8 * lsh0_870[k]
                    - f_9 * lsh1_870[k]
                    + f_3 * pc_x[k] * lsi_1157[k];

        t_1486[k] = f_6 * lsh0_871[k]
                    - f_7 * lsh1_871[k]
                    + f_3 * pc_x[k] * lsi_1158[k];

        t_1487[k] = f_6 * lsh0_872[k]
                    - f_7 * lsh1_872[k]
                    + f_3 * pc_x[k] * lsi_1159[k];
    }

#pragma omp simd aligned(t_1488, t_1489, t_1490, pc_x, lsh0_873, lsh0_874, lsh0_875, lsh1_873, \
                         lsh1_874, lsh1_875, lsi_1160, lsi_1161, \
                         lsi_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1488[k] = f_6 * lsh0_873[k]
                    - f_7 * lsh1_873[k]
                    + f_3 * pc_x[k] * lsi_1160[k];

        t_1489[k] = f_6 * lsh0_874[k]
                    - f_7 * lsh1_874[k]
                    + f_3 * pc_x[k] * lsi_1161[k];

        t_1490[k] = f_6 * lsh0_875[k]
                    - f_7 * lsh1_875[k]
                    + f_3 * pc_x[k] * lsi_1162[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, pc_x, lsh0_876, lsh0_877, lsh0_878, lsh1_876, \
                         lsh1_877, lsh1_878, lsi_1163, lsi_1164, \
                         lsi_1165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = f_4 * lsh0_876[k]
                    - f_5 * lsh1_876[k]
                    + f_3 * pc_x[k] * lsi_1163[k];

        t_1492[k] = f_4 * lsh0_877[k]
                    - f_5 * lsh1_877[k]
                    + f_3 * pc_x[k] * lsi_1164[k];

        t_1493[k] = f_4 * lsh0_878[k]
                    - f_5 * lsh1_878[k]
                    + f_3 * pc_x[k] * lsi_1165[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, t_1497, pc_x, lsh0_879, lsh0_880, lsh0_881, \
                         lsh1_879, lsh1_880, lsh1_881, lsi_1166, lsi_1167, lsi_1168, \
                         lsi_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = f_4 * lsh0_879[k]
                    - f_5 * lsh1_879[k]
                    + f_3 * pc_x[k] * lsi_1166[k];

        t_1495[k] = f_4 * lsh0_880[k]
                    - f_5 * lsh1_880[k]
                    + f_3 * pc_x[k] * lsi_1167[k];

        t_1496[k] = f_4 * lsh0_881[k]
                    - f_5 * lsh1_881[k]
                    + f_3 * pc_x[k] * lsi_1168[k];

        t_1497[k] = f_3 * pc_x[k] * lsi_1169[k];
    }

#pragma omp simd aligned(t_1498, t_1499, t_1500, t_1501, t_1502, t_1503, pc_x, lsi_1170, \
                         lsi_1171, lsi_1172, lsi_1173, lsi_1174, \
                         lsi_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1498[k] = f_3 * pc_x[k] * lsi_1170[k];

        t_1499[k] = f_3 * pc_x[k] * lsi_1171[k];

        t_1500[k] = f_3 * pc_x[k] * lsi_1172[k];

        t_1501[k] = f_3 * pc_x[k] * lsi_1173[k];

        t_1502[k] = f_3 * pc_x[k] * lsi_1174[k];

        t_1503[k] = f_3 * pc_x[k] * lsi_1175[k];
    }

#pragma omp simd aligned(t_1504, t_1505, t_1506, pc_y, pc_z, ksi_917, ksi_945, ksi_947, \
                         lsh0_876, lsh0_878, lsh1_876, lsh1_878, lsi_1169, \
                         lsi_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1504[k] = f_15 * ksi_945[k]
                    + f_1 * lsh0_876[k]
                    - f_2 * lsh1_876[k]
                    + f_3 * pc_y[k] * lsi_1169[k];

        t_1505[k] = f_17 * ksi_917[k]
                    + f_3 * pc_z[k] * lsi_1169[k];

        t_1506[k] = f_15 * ksi_947[k]
                    + f_10 * lsh0_878[k]
                    - f_11 * lsh1_878[k]
                    + f_3 * pc_y[k] * lsi_1171[k];
    }

#pragma omp simd aligned(t_1507, t_1508, t_1509, pc_y, ksi_948, ksi_949, ksi_950, lsh0_879, \
                         lsh0_880, lsh0_881, lsh1_879, lsh1_880, lsh1_881, lsi_1172, lsi_1173, \
                         lsi_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1507[k] = f_15 * ksi_948[k]
                    + f_8 * lsh0_879[k]
                    - f_9 * lsh1_879[k]
                    + f_3 * pc_y[k] * lsi_1172[k];

        t_1508[k] = f_15 * ksi_949[k]
                    + f_6 * lsh0_880[k]
                    - f_7 * lsh1_880[k]
                    + f_3 * pc_y[k] * lsi_1173[k];

        t_1509[k] = f_15 * ksi_950[k]
                    + f_4 * lsh0_881[k]
                    - f_5 * lsh1_881[k]
                    + f_3 * pc_y[k] * lsi_1174[k];
    }

#pragma omp simd aligned(t_1510, t_1511, t_1512, pc_x, pc_y, pc_z, ksi_923, ksi_951, lsh0_881, \
                         lsh0_882, lsh1_881, lsh1_882, lsi_1175, \
                         lsi_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1510[k] = f_15 * ksi_951[k]
                    + f_3 * pc_y[k] * lsi_1175[k];

        t_1511[k] = f_17 * ksi_923[k]
                    + f_1 * lsh0_881[k]
                    - f_2 * lsh1_881[k]
                    + f_3 * pc_z[k] * lsi_1175[k];

        t_1512[k] = f_1 * lsh0_882[k]
                    - f_2 * lsh1_882[k]
                    + f_3 * pc_x[k] * lsi_1176[k];
    }

#pragma omp simd aligned(t_1513, t_1514, t_1515, pc_x, lsh0_883, lsh0_884, lsh0_885, lsh1_883, \
                         lsh1_884, lsh1_885, lsi_1177, lsi_1178, \
                         lsi_1179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1513[k] = f_19 * lsh0_883[k]
                    - f_20 * lsh1_883[k]
                    + f_3 * pc_x[k] * lsi_1177[k];

        t_1514[k] = f_19 * lsh0_884[k]
                    - f_20 * lsh1_884[k]
                    + f_3 * pc_x[k] * lsi_1178[k];

        t_1515[k] = f_10 * lsh0_885[k]
                    - f_11 * lsh1_885[k]
                    + f_3 * pc_x[k] * lsi_1179[k];
    }

#pragma omp simd aligned(t_1516, t_1517, t_1518, pc_x, lsh0_886, lsh0_887, lsh0_888, lsh1_886, \
                         lsh1_887, lsh1_888, lsi_1180, lsi_1181, \
                         lsi_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1516[k] = f_10 * lsh0_886[k]
                    - f_11 * lsh1_886[k]
                    + f_3 * pc_x[k] * lsi_1180[k];

        t_1517[k] = f_10 * lsh0_887[k]
                    - f_11 * lsh1_887[k]
                    + f_3 * pc_x[k] * lsi_1181[k];

        t_1518[k] = f_8 * lsh0_888[k]
                    - f_9 * lsh1_888[k]
                    + f_3 * pc_x[k] * lsi_1182[k];
    }

#pragma omp simd aligned(t_1519, t_1520, t_1521, pc_x, lsh0_889, lsh0_890, lsh0_891, lsh1_889, \
                         lsh1_890, lsh1_891, lsi_1183, lsi_1184, \
                         lsi_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1519[k] = f_8 * lsh0_889[k]
                    - f_9 * lsh1_889[k]
                    + f_3 * pc_x[k] * lsi_1183[k];

        t_1520[k] = f_8 * lsh0_890[k]
                    - f_9 * lsh1_890[k]
                    + f_3 * pc_x[k] * lsi_1184[k];

        t_1521[k] = f_8 * lsh0_891[k]
                    - f_9 * lsh1_891[k]
                    + f_3 * pc_x[k] * lsi_1185[k];
    }
}

static auto
compute_prim_lsk_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksk0,
                                                           const size_t ksi, const size_t ksk1,
                                                           const size_t lsh0, const size_t lsh1,
                                                           const size_t lsi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_17 = 2.5 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 2.5 / gamma;
    const auto f_20 = 2.5 * p / (gamma * q);
    const auto f_21 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk0_1260 = buffer.data(ksk0 + 1260);
    const auto *ksk0_1262 = buffer.data(ksk0 + 1262);
    const auto *ksk0_1265 = buffer.data(ksk0 + 1265);
    const auto *ksk0_1269 = buffer.data(ksk0 + 1269);
    const auto *ksk0_1274 = buffer.data(ksk0 + 1274);
    const auto *ksk0_1280 = buffer.data(ksk0 + 1280);
    const auto *ksk0_1288 = buffer.data(ksk0 + 1288);
    const auto *ksk0_1290 = buffer.data(ksk0 + 1290);
    const auto *ksk0_1291 = buffer.data(ksk0 + 1291);
    const auto *ksk0_1292 = buffer.data(ksk0 + 1292);
    const auto *ksk0_1293 = buffer.data(ksk0 + 1293);
    const auto *ksk0_1295 = buffer.data(ksk0 + 1295);

    const auto *ksi_945 = buffer.data(ksi + 945);
    const auto *ksi_951 = buffer.data(ksi + 951);
    const auto *ksi_973 = buffer.data(ksi + 973);
    const auto *ksi_975 = buffer.data(ksi + 975);
    const auto *ksi_976 = buffer.data(ksi + 976);
    const auto *ksi_977 = buffer.data(ksi + 977);
    const auto *ksi_978 = buffer.data(ksi + 978);
    const auto *ksi_979 = buffer.data(ksi + 979);
    const auto *ksi_1001 = buffer.data(ksi + 1001);
    const auto *ksi_1003 = buffer.data(ksi + 1003);
    const auto *ksi_1004 = buffer.data(ksi + 1004);
    const auto *ksi_1005 = buffer.data(ksi + 1005);
    const auto *ksi_1006 = buffer.data(ksi + 1006);
    const auto *ksi_1007 = buffer.data(ksi + 1007);

    const auto *ksk1_1260 = buffer.data(ksk1 + 1260);
    const auto *ksk1_1262 = buffer.data(ksk1 + 1262);
    const auto *ksk1_1265 = buffer.data(ksk1 + 1265);
    const auto *ksk1_1269 = buffer.data(ksk1 + 1269);
    const auto *ksk1_1274 = buffer.data(ksk1 + 1274);
    const auto *ksk1_1280 = buffer.data(ksk1 + 1280);
    const auto *ksk1_1288 = buffer.data(ksk1 + 1288);
    const auto *ksk1_1290 = buffer.data(ksk1 + 1290);
    const auto *ksk1_1291 = buffer.data(ksk1 + 1291);
    const auto *ksk1_1292 = buffer.data(ksk1 + 1292);
    const auto *ksk1_1293 = buffer.data(ksk1 + 1293);
    const auto *ksk1_1295 = buffer.data(ksk1 + 1295);

    const auto *lsh0_892 = buffer.data(lsh0 + 892);
    const auto *lsh0_893 = buffer.data(lsh0 + 893);
    const auto *lsh0_894 = buffer.data(lsh0 + 894);
    const auto *lsh0_895 = buffer.data(lsh0 + 895);
    const auto *lsh0_896 = buffer.data(lsh0 + 896);
    const auto *lsh0_897 = buffer.data(lsh0 + 897);
    const auto *lsh0_898 = buffer.data(lsh0 + 898);
    const auto *lsh0_899 = buffer.data(lsh0 + 899);
    const auto *lsh0_900 = buffer.data(lsh0 + 900);
    const auto *lsh0_901 = buffer.data(lsh0 + 901);
    const auto *lsh0_902 = buffer.data(lsh0 + 902);
    const auto *lsh0_904 = buffer.data(lsh0 + 904);
    const auto *lsh0_906 = buffer.data(lsh0 + 906);
    const auto *lsh0_907 = buffer.data(lsh0 + 907);
    const auto *lsh0_909 = buffer.data(lsh0 + 909);
    const auto *lsh0_910 = buffer.data(lsh0 + 910);
    const auto *lsh0_911 = buffer.data(lsh0 + 911);
    const auto *lsh0_913 = buffer.data(lsh0 + 913);
    const auto *lsh0_914 = buffer.data(lsh0 + 914);
    const auto *lsh0_915 = buffer.data(lsh0 + 915);
    const auto *lsh0_916 = buffer.data(lsh0 + 916);
    const auto *lsh0_918 = buffer.data(lsh0 + 918);
    const auto *lsh0_919 = buffer.data(lsh0 + 919);
    const auto *lsh0_920 = buffer.data(lsh0 + 920);
    const auto *lsh0_921 = buffer.data(lsh0 + 921);
    const auto *lsh0_922 = buffer.data(lsh0 + 922);
    const auto *lsh0_924 = buffer.data(lsh0 + 924);
    const auto *lsh0_926 = buffer.data(lsh0 + 926);
    const auto *lsh0_927 = buffer.data(lsh0 + 927);
    const auto *lsh0_929 = buffer.data(lsh0 + 929);
    const auto *lsh0_930 = buffer.data(lsh0 + 930);
    const auto *lsh0_931 = buffer.data(lsh0 + 931);
    const auto *lsh0_933 = buffer.data(lsh0 + 933);
    const auto *lsh0_934 = buffer.data(lsh0 + 934);
    const auto *lsh0_935 = buffer.data(lsh0 + 935);
    const auto *lsh0_936 = buffer.data(lsh0 + 936);
    const auto *lsh0_938 = buffer.data(lsh0 + 938);
    const auto *lsh0_939 = buffer.data(lsh0 + 939);
    const auto *lsh0_940 = buffer.data(lsh0 + 940);
    const auto *lsh0_941 = buffer.data(lsh0 + 941);
    const auto *lsh0_942 = buffer.data(lsh0 + 942);
    const auto *lsh0_943 = buffer.data(lsh0 + 943);
    const auto *lsh0_944 = buffer.data(lsh0 + 944);

    const auto *lsh1_892 = buffer.data(lsh1 + 892);
    const auto *lsh1_893 = buffer.data(lsh1 + 893);
    const auto *lsh1_894 = buffer.data(lsh1 + 894);
    const auto *lsh1_895 = buffer.data(lsh1 + 895);
    const auto *lsh1_896 = buffer.data(lsh1 + 896);
    const auto *lsh1_897 = buffer.data(lsh1 + 897);
    const auto *lsh1_898 = buffer.data(lsh1 + 898);
    const auto *lsh1_899 = buffer.data(lsh1 + 899);
    const auto *lsh1_900 = buffer.data(lsh1 + 900);
    const auto *lsh1_901 = buffer.data(lsh1 + 901);
    const auto *lsh1_902 = buffer.data(lsh1 + 902);
    const auto *lsh1_904 = buffer.data(lsh1 + 904);
    const auto *lsh1_906 = buffer.data(lsh1 + 906);
    const auto *lsh1_907 = buffer.data(lsh1 + 907);
    const auto *lsh1_909 = buffer.data(lsh1 + 909);
    const auto *lsh1_910 = buffer.data(lsh1 + 910);
    const auto *lsh1_911 = buffer.data(lsh1 + 911);
    const auto *lsh1_913 = buffer.data(lsh1 + 913);
    const auto *lsh1_914 = buffer.data(lsh1 + 914);
    const auto *lsh1_915 = buffer.data(lsh1 + 915);
    const auto *lsh1_916 = buffer.data(lsh1 + 916);
    const auto *lsh1_918 = buffer.data(lsh1 + 918);
    const auto *lsh1_919 = buffer.data(lsh1 + 919);
    const auto *lsh1_920 = buffer.data(lsh1 + 920);
    const auto *lsh1_921 = buffer.data(lsh1 + 921);
    const auto *lsh1_922 = buffer.data(lsh1 + 922);
    const auto *lsh1_924 = buffer.data(lsh1 + 924);
    const auto *lsh1_926 = buffer.data(lsh1 + 926);
    const auto *lsh1_927 = buffer.data(lsh1 + 927);
    const auto *lsh1_929 = buffer.data(lsh1 + 929);
    const auto *lsh1_930 = buffer.data(lsh1 + 930);
    const auto *lsh1_931 = buffer.data(lsh1 + 931);
    const auto *lsh1_933 = buffer.data(lsh1 + 933);
    const auto *lsh1_934 = buffer.data(lsh1 + 934);
    const auto *lsh1_935 = buffer.data(lsh1 + 935);
    const auto *lsh1_936 = buffer.data(lsh1 + 936);
    const auto *lsh1_938 = buffer.data(lsh1 + 938);
    const auto *lsh1_939 = buffer.data(lsh1 + 939);
    const auto *lsh1_940 = buffer.data(lsh1 + 940);
    const auto *lsh1_941 = buffer.data(lsh1 + 941);
    const auto *lsh1_942 = buffer.data(lsh1 + 942);
    const auto *lsh1_943 = buffer.data(lsh1 + 943);
    const auto *lsh1_944 = buffer.data(lsh1 + 944);

    const auto *lsi_1186 = buffer.data(lsi + 1186);
    const auto *lsi_1187 = buffer.data(lsi + 1187);
    const auto *lsi_1188 = buffer.data(lsi + 1188);
    const auto *lsi_1189 = buffer.data(lsi + 1189);
    const auto *lsi_1190 = buffer.data(lsi + 1190);
    const auto *lsi_1191 = buffer.data(lsi + 1191);
    const auto *lsi_1192 = buffer.data(lsi + 1192);
    const auto *lsi_1193 = buffer.data(lsi + 1193);
    const auto *lsi_1194 = buffer.data(lsi + 1194);
    const auto *lsi_1195 = buffer.data(lsi + 1195);
    const auto *lsi_1196 = buffer.data(lsi + 1196);
    const auto *lsi_1197 = buffer.data(lsi + 1197);
    const auto *lsi_1198 = buffer.data(lsi + 1198);
    const auto *lsi_1199 = buffer.data(lsi + 1199);
    const auto *lsi_1200 = buffer.data(lsi + 1200);
    const auto *lsi_1201 = buffer.data(lsi + 1201);
    const auto *lsi_1202 = buffer.data(lsi + 1202);
    const auto *lsi_1203 = buffer.data(lsi + 1203);
    const auto *lsi_1205 = buffer.data(lsi + 1205);
    const auto *lsi_1207 = buffer.data(lsi + 1207);
    const auto *lsi_1208 = buffer.data(lsi + 1208);
    const auto *lsi_1210 = buffer.data(lsi + 1210);
    const auto *lsi_1211 = buffer.data(lsi + 1211);
    const auto *lsi_1212 = buffer.data(lsi + 1212);
    const auto *lsi_1214 = buffer.data(lsi + 1214);
    const auto *lsi_1215 = buffer.data(lsi + 1215);
    const auto *lsi_1216 = buffer.data(lsi + 1216);
    const auto *lsi_1217 = buffer.data(lsi + 1217);
    const auto *lsi_1219 = buffer.data(lsi + 1219);
    const auto *lsi_1220 = buffer.data(lsi + 1220);
    const auto *lsi_1221 = buffer.data(lsi + 1221);
    const auto *lsi_1222 = buffer.data(lsi + 1222);
    const auto *lsi_1223 = buffer.data(lsi + 1223);
    const auto *lsi_1225 = buffer.data(lsi + 1225);
    const auto *lsi_1226 = buffer.data(lsi + 1226);
    const auto *lsi_1227 = buffer.data(lsi + 1227);
    const auto *lsi_1228 = buffer.data(lsi + 1228);
    const auto *lsi_1229 = buffer.data(lsi + 1229);
    const auto *lsi_1230 = buffer.data(lsi + 1230);
    const auto *lsi_1231 = buffer.data(lsi + 1231);
    const auto *lsi_1232 = buffer.data(lsi + 1232);
    const auto *lsi_1234 = buffer.data(lsi + 1234);
    const auto *lsi_1235 = buffer.data(lsi + 1235);
    const auto *lsi_1237 = buffer.data(lsi + 1237);
    const auto *lsi_1238 = buffer.data(lsi + 1238);
    const auto *lsi_1239 = buffer.data(lsi + 1239);
    const auto *lsi_1241 = buffer.data(lsi + 1241);
    const auto *lsi_1242 = buffer.data(lsi + 1242);
    const auto *lsi_1243 = buffer.data(lsi + 1243);
    const auto *lsi_1244 = buffer.data(lsi + 1244);
    const auto *lsi_1246 = buffer.data(lsi + 1246);
    const auto *lsi_1247 = buffer.data(lsi + 1247);
    const auto *lsi_1248 = buffer.data(lsi + 1248);
    const auto *lsi_1249 = buffer.data(lsi + 1249);
    const auto *lsi_1250 = buffer.data(lsi + 1250);
    const auto *lsi_1252 = buffer.data(lsi + 1252);
    const auto *lsi_1253 = buffer.data(lsi + 1253);
    const auto *lsi_1254 = buffer.data(lsi + 1254);
    const auto *lsi_1255 = buffer.data(lsi + 1255);
    const auto *lsi_1256 = buffer.data(lsi + 1256);
    const auto *lsi_1257 = buffer.data(lsi + 1257);
    const auto *lsi_1258 = buffer.data(lsi + 1258);
    const auto *lsi_1259 = buffer.data(lsi + 1259);

#pragma omp simd aligned(t_1522, t_1523, t_1524, pc_x, lsh0_892, lsh0_893, lsh0_894, lsh1_892, \
                         lsh1_893, lsh1_894, lsi_1186, lsi_1187, \
                         lsi_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1522[k] = f_6 * lsh0_892[k]
                    - f_7 * lsh1_892[k]
                    + f_3 * pc_x[k] * lsi_1186[k];

        t_1523[k] = f_6 * lsh0_893[k]
                    - f_7 * lsh1_893[k]
                    + f_3 * pc_x[k] * lsi_1187[k];

        t_1524[k] = f_6 * lsh0_894[k]
                    - f_7 * lsh1_894[k]
                    + f_3 * pc_x[k] * lsi_1188[k];
    }

#pragma omp simd aligned(t_1525, t_1526, t_1527, pc_x, lsh0_895, lsh0_896, lsh0_897, lsh1_895, \
                         lsh1_896, lsh1_897, lsi_1189, lsi_1190, \
                         lsi_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1525[k] = f_6 * lsh0_895[k]
                    - f_7 * lsh1_895[k]
                    + f_3 * pc_x[k] * lsi_1189[k];

        t_1526[k] = f_6 * lsh0_896[k]
                    - f_7 * lsh1_896[k]
                    + f_3 * pc_x[k] * lsi_1190[k];

        t_1527[k] = f_4 * lsh0_897[k]
                    - f_5 * lsh1_897[k]
                    + f_3 * pc_x[k] * lsi_1191[k];
    }

#pragma omp simd aligned(t_1528, t_1529, t_1530, pc_x, lsh0_898, lsh0_899, lsh0_900, lsh1_898, \
                         lsh1_899, lsh1_900, lsi_1192, lsi_1193, \
                         lsi_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1528[k] = f_4 * lsh0_898[k]
                    - f_5 * lsh1_898[k]
                    + f_3 * pc_x[k] * lsi_1192[k];

        t_1529[k] = f_4 * lsh0_899[k]
                    - f_5 * lsh1_899[k]
                    + f_3 * pc_x[k] * lsi_1193[k];

        t_1530[k] = f_4 * lsh0_900[k]
                    - f_5 * lsh1_900[k]
                    + f_3 * pc_x[k] * lsi_1194[k];
    }

#pragma omp simd aligned(t_1531, t_1532, t_1533, t_1534, t_1535, pc_x, lsh0_901, lsh0_902, \
                         lsh1_901, lsh1_902, lsi_1195, lsi_1196, lsi_1197, lsi_1198, \
                         lsi_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1531[k] = f_4 * lsh0_901[k]
                    - f_5 * lsh1_901[k]
                    + f_3 * pc_x[k] * lsi_1195[k];

        t_1532[k] = f_4 * lsh0_902[k]
                    - f_5 * lsh1_902[k]
                    + f_3 * pc_x[k] * lsi_1196[k];

        t_1533[k] = f_3 * pc_x[k] * lsi_1197[k];

        t_1534[k] = f_3 * pc_x[k] * lsi_1198[k];

        t_1535[k] = f_3 * pc_x[k] * lsi_1199[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, t_1539, t_1540, pc_x, pc_y, ksi_973, \
                         lsh0_897, lsh1_897, lsi_1197, lsi_1200, lsi_1201, lsi_1202, \
                         lsi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_3 * pc_x[k] * lsi_1200[k];

        t_1537[k] = f_3 * pc_x[k] * lsi_1201[k];

        t_1538[k] = f_3 * pc_x[k] * lsi_1202[k];

        t_1539[k] = f_3 * pc_x[k] * lsi_1203[k];

        t_1540[k] = f_14 * ksi_973[k]
                    + f_1 * lsh0_897[k]
                    - f_2 * lsh1_897[k]
                    + f_3 * pc_y[k] * lsi_1197[k];
    }

#pragma omp simd aligned(t_1541, t_1542, t_1543, pc_y, pc_z, ksi_945, ksi_975, ksi_976, \
                         lsh0_899, lsh0_900, lsh1_899, lsh1_900, lsi_1197, lsi_1199, \
                         lsi_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1541[k] = f_21 * ksi_945[k]
                    + f_3 * pc_z[k] * lsi_1197[k];

        t_1542[k] = f_14 * ksi_975[k]
                    + f_10 * lsh0_899[k]
                    - f_11 * lsh1_899[k]
                    + f_3 * pc_y[k] * lsi_1199[k];

        t_1543[k] = f_14 * ksi_976[k]
                    + f_8 * lsh0_900[k]
                    - f_9 * lsh1_900[k]
                    + f_3 * pc_y[k] * lsi_1200[k];
    }

#pragma omp simd aligned(t_1544, t_1545, t_1546, pc_y, ksi_977, ksi_978, ksi_979, lsh0_901, \
                         lsh0_902, lsh1_901, lsh1_902, lsi_1201, lsi_1202, \
                         lsi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1544[k] = f_14 * ksi_977[k]
                    + f_6 * lsh0_901[k]
                    - f_7 * lsh1_901[k]
                    + f_3 * pc_y[k] * lsi_1201[k];

        t_1545[k] = f_14 * ksi_978[k]
                    + f_4 * lsh0_902[k]
                    - f_5 * lsh1_902[k]
                    + f_3 * pc_y[k] * lsi_1202[k];

        t_1546[k] = f_14 * ksi_979[k]
                    + f_3 * pc_y[k] * lsi_1203[k];
    }

#pragma omp simd aligned(t_1547, t_1548, t_1549, pa_y, pc_x, pc_y, pc_z, ksk0_1260, ksi_951, \
                         ksk1_1260, lsh0_902, lsh0_904, lsh1_902, lsh1_904, lsi_1203, \
                         lsi_1205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1547[k] = f_21 * ksi_951[k]
                    + f_1 * lsh0_902[k]
                    - f_2 * lsh1_902[k]
                    + f_3 * pc_z[k] * lsi_1203[k];

        t_1548[k] = pa_y[k] * ksk0_1260[k]
                    - f_12 * pc_y[k] * ksk1_1260[k];

        t_1549[k] = f_19 * lsh0_904[k]
                    - f_20 * lsh1_904[k]
                    + f_3 * pc_x[k] * lsi_1205[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, pa_y, pc_x, pc_y, ksk0_1262, ksk1_1262, \
                         lsh0_906, lsh0_907, lsh1_906, lsh1_907, lsi_1207, \
                         lsi_1208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = pa_y[k] * ksk0_1262[k]
                    - f_12 * pc_y[k] * ksk1_1262[k];

        t_1551[k] = f_10 * lsh0_906[k]
                    - f_11 * lsh1_906[k]
                    + f_3 * pc_x[k] * lsi_1207[k];

        t_1552[k] = f_10 * lsh0_907[k]
                    - f_11 * lsh1_907[k]
                    + f_3 * pc_x[k] * lsi_1208[k];
    }

#pragma omp simd aligned(t_1553, t_1554, t_1555, pa_y, pc_x, pc_y, ksk0_1265, ksk1_1265, \
                         lsh0_909, lsh0_910, lsh1_909, lsh1_910, lsi_1210, \
                         lsi_1211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1553[k] = pa_y[k] * ksk0_1265[k]
                    - f_12 * pc_y[k] * ksk1_1265[k];

        t_1554[k] = f_8 * lsh0_909[k]
                    - f_9 * lsh1_909[k]
                    + f_3 * pc_x[k] * lsi_1210[k];

        t_1555[k] = f_8 * lsh0_910[k]
                    - f_9 * lsh1_910[k]
                    + f_3 * pc_x[k] * lsi_1211[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, pa_y, pc_x, pc_y, ksk0_1269, ksk1_1269, \
                         lsh0_911, lsh0_913, lsh1_911, lsh1_913, lsi_1212, \
                         lsi_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_8 * lsh0_911[k]
                    - f_9 * lsh1_911[k]
                    + f_3 * pc_x[k] * lsi_1212[k];

        t_1557[k] = pa_y[k] * ksk0_1269[k]
                    - f_12 * pc_y[k] * ksk1_1269[k];

        t_1558[k] = f_6 * lsh0_913[k]
                    - f_7 * lsh1_913[k]
                    + f_3 * pc_x[k] * lsi_1214[k];
    }

#pragma omp simd aligned(t_1559, t_1560, t_1561, pc_x, lsh0_914, lsh0_915, lsh0_916, lsh1_914, \
                         lsh1_915, lsh1_916, lsi_1215, lsi_1216, \
                         lsi_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1559[k] = f_6 * lsh0_914[k]
                    - f_7 * lsh1_914[k]
                    + f_3 * pc_x[k] * lsi_1215[k];

        t_1560[k] = f_6 * lsh0_915[k]
                    - f_7 * lsh1_915[k]
                    + f_3 * pc_x[k] * lsi_1216[k];

        t_1561[k] = f_6 * lsh0_916[k]
                    - f_7 * lsh1_916[k]
                    + f_3 * pc_x[k] * lsi_1217[k];
    }

#pragma omp simd aligned(t_1562, t_1563, t_1564, pa_y, pc_x, pc_y, ksk0_1274, ksk1_1274, \
                         lsh0_918, lsh0_919, lsh1_918, lsh1_919, lsi_1219, \
                         lsi_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1562[k] = pa_y[k] * ksk0_1274[k]
                    - f_12 * pc_y[k] * ksk1_1274[k];

        t_1563[k] = f_4 * lsh0_918[k]
                    - f_5 * lsh1_918[k]
                    + f_3 * pc_x[k] * lsi_1219[k];

        t_1564[k] = f_4 * lsh0_919[k]
                    - f_5 * lsh1_919[k]
                    + f_3 * pc_x[k] * lsi_1220[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, pc_x, lsh0_920, lsh0_921, lsh0_922, lsh1_920, \
                         lsh1_921, lsh1_922, lsi_1221, lsi_1222, \
                         lsi_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = f_4 * lsh0_920[k]
                    - f_5 * lsh1_920[k]
                    + f_3 * pc_x[k] * lsi_1221[k];

        t_1566[k] = f_4 * lsh0_921[k]
                    - f_5 * lsh1_921[k]
                    + f_3 * pc_x[k] * lsi_1222[k];

        t_1567[k] = f_4 * lsh0_922[k]
                    - f_5 * lsh1_922[k]
                    + f_3 * pc_x[k] * lsi_1223[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, t_1571, t_1572, t_1573, pa_y, pc_x, pc_y, \
                         ksk0_1280, ksk1_1280, lsi_1225, lsi_1226, lsi_1227, lsi_1228, \
                         lsi_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = pa_y[k] * ksk0_1280[k]
                    - f_12 * pc_y[k] * ksk1_1280[k];

        t_1569[k] = f_3 * pc_x[k] * lsi_1225[k];

        t_1570[k] = f_3 * pc_x[k] * lsi_1226[k];

        t_1571[k] = f_3 * pc_x[k] * lsi_1227[k];

        t_1572[k] = f_3 * pc_x[k] * lsi_1228[k];

        t_1573[k] = f_3 * pc_x[k] * lsi_1229[k];
    }

#pragma omp simd aligned(t_1574, t_1575, t_1576, t_1577, pa_y, pc_x, pc_y, pc_z, ksk0_1288, \
                         ksi_973, ksi_1001, ksk1_1288, lsi_1225, lsi_1230, \
                         lsi_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = f_3 * pc_x[k] * lsi_1230[k];

        t_1575[k] = f_3 * pc_x[k] * lsi_1231[k];

        t_1576[k] = pa_y[k] * ksk0_1288[k]
                    + f_18 * ksi_1001[k]
                    - f_12 * pc_y[k] * ksk1_1288[k];

        t_1577[k] = f_18 * ksi_973[k]
                    + f_3 * pc_z[k] * lsi_1225[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pa_y, pc_y, ksk0_1290, ksk0_1291, ksk0_1292, \
                         ksi_1003, ksi_1004, ksi_1005, ksk1_1290, ksk1_1291, \
                         ksk1_1292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = pa_y[k] * ksk0_1290[k]
                    + f_17 * ksi_1003[k]
                    - f_12 * pc_y[k] * ksk1_1290[k];

        t_1579[k] = pa_y[k] * ksk0_1291[k]
                    + f_16 * ksi_1004[k]
                    - f_12 * pc_y[k] * ksk1_1291[k];

        t_1580[k] = pa_y[k] * ksk0_1292[k]
                    + f_15 * ksi_1005[k]
                    - f_12 * pc_y[k] * ksk1_1292[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, pa_y, pc_y, ksk0_1293, ksk0_1295, ksi_1006, \
                         ksi_1007, ksk1_1293, ksk1_1295, lsi_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = pa_y[k] * ksk0_1293[k]
                    + f_14 * ksi_1006[k]
                    - f_12 * pc_y[k] * ksk1_1293[k];

        t_1582[k] = f_13 * ksi_1007[k]
                    + f_3 * pc_y[k] * lsi_1231[k];

        t_1583[k] = pa_y[k] * ksk0_1295[k]
                    - f_12 * pc_y[k] * ksk1_1295[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, t_1587, t_1588, pc_x, pc_y, lsh0_924, \
                         lsh0_926, lsh0_927, lsh1_924, lsh1_926, lsh1_927, lsi_1232, lsi_1234, \
                         lsi_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = f_1 * lsh0_924[k]
                    - f_2 * lsh1_924[k]
                    + f_3 * pc_x[k] * lsi_1232[k];

        t_1585[k] = f_3 * pc_y[k] * lsi_1232[k];

        t_1586[k] = f_19 * lsh0_926[k]
                    - f_20 * lsh1_926[k]
                    + f_3 * pc_x[k] * lsi_1234[k];

        t_1587[k] = f_10 * lsh0_927[k]
                    - f_11 * lsh1_927[k]
                    + f_3 * pc_x[k] * lsi_1235[k];

        t_1588[k] = f_3 * pc_y[k] * lsi_1234[k];
    }

#pragma omp simd aligned(t_1589, t_1590, t_1591, t_1592, pc_x, pc_y, lsh0_929, lsh0_930, \
                         lsh0_931, lsh1_929, lsh1_930, lsh1_931, lsi_1237, lsi_1238, \
                         lsi_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1589[k] = f_10 * lsh0_929[k]
                    - f_11 * lsh1_929[k]
                    + f_3 * pc_x[k] * lsi_1237[k];

        t_1590[k] = f_8 * lsh0_930[k]
                    - f_9 * lsh1_930[k]
                    + f_3 * pc_x[k] * lsi_1238[k];

        t_1591[k] = f_8 * lsh0_931[k]
                    - f_9 * lsh1_931[k]
                    + f_3 * pc_x[k] * lsi_1239[k];

        t_1592[k] = f_3 * pc_y[k] * lsi_1237[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, pc_x, lsh0_933, lsh0_934, lsh0_935, lsh1_933, \
                         lsh1_934, lsh1_935, lsi_1241, lsi_1242, \
                         lsi_1243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = f_8 * lsh0_933[k]
                    - f_9 * lsh1_933[k]
                    + f_3 * pc_x[k] * lsi_1241[k];

        t_1594[k] = f_6 * lsh0_934[k]
                    - f_7 * lsh1_934[k]
                    + f_3 * pc_x[k] * lsi_1242[k];

        t_1595[k] = f_6 * lsh0_935[k]
                    - f_7 * lsh1_935[k]
                    + f_3 * pc_x[k] * lsi_1243[k];
    }

#pragma omp simd aligned(t_1596, t_1597, t_1598, t_1599, pc_x, pc_y, lsh0_936, lsh0_938, \
                         lsh0_939, lsh1_936, lsh1_938, lsh1_939, lsi_1241, lsi_1244, lsi_1246, \
                         lsi_1247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = f_6 * lsh0_936[k]
                    - f_7 * lsh1_936[k]
                    + f_3 * pc_x[k] * lsi_1244[k];

        t_1597[k] = f_3 * pc_y[k] * lsi_1241[k];

        t_1598[k] = f_6 * lsh0_938[k]
                    - f_7 * lsh1_938[k]
                    + f_3 * pc_x[k] * lsi_1246[k];

        t_1599[k] = f_4 * lsh0_939[k]
                    - f_5 * lsh1_939[k]
                    + f_3 * pc_x[k] * lsi_1247[k];
    }

#pragma omp simd aligned(t_1600, t_1601, t_1602, t_1603, pc_x, pc_y, lsh0_940, lsh0_941, \
                         lsh0_942, lsh1_940, lsh1_941, lsh1_942, lsi_1246, lsi_1248, lsi_1249, \
                         lsi_1250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1600[k] = f_4 * lsh0_940[k]
                    - f_5 * lsh1_940[k]
                    + f_3 * pc_x[k] * lsi_1248[k];

        t_1601[k] = f_4 * lsh0_941[k]
                    - f_5 * lsh1_941[k]
                    + f_3 * pc_x[k] * lsi_1249[k];

        t_1602[k] = f_4 * lsh0_942[k]
                    - f_5 * lsh1_942[k]
                    + f_3 * pc_x[k] * lsi_1250[k];

        t_1603[k] = f_3 * pc_y[k] * lsi_1246[k];
    }

#pragma omp simd aligned(t_1604, t_1605, t_1606, t_1607, t_1608, t_1609, pc_x, lsh0_944, \
                         lsh1_944, lsi_1252, lsi_1253, lsi_1254, lsi_1255, lsi_1256, \
                         lsi_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1604[k] = f_4 * lsh0_944[k]
                    - f_5 * lsh1_944[k]
                    + f_3 * pc_x[k] * lsi_1252[k];

        t_1605[k] = f_3 * pc_x[k] * lsi_1253[k];

        t_1606[k] = f_3 * pc_x[k] * lsi_1254[k];

        t_1607[k] = f_3 * pc_x[k] * lsi_1255[k];

        t_1608[k] = f_3 * pc_x[k] * lsi_1256[k];

        t_1609[k] = f_3 * pc_x[k] * lsi_1257[k];
    }

#pragma omp simd aligned(t_1610, t_1611, t_1612, t_1613, pc_x, pc_y, lsh0_939, lsh0_940, \
                         lsh1_939, lsh1_940, lsi_1253, lsi_1254, lsi_1258, \
                         lsi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1610[k] = f_3 * pc_x[k] * lsi_1258[k];

        t_1611[k] = f_3 * pc_x[k] * lsi_1259[k];

        t_1612[k] = f_1 * lsh0_939[k]
                    - f_2 * lsh1_939[k]
                    + f_3 * pc_y[k] * lsi_1253[k];

        t_1613[k] = f_19 * lsh0_940[k]
                    - f_20 * lsh1_940[k]
                    + f_3 * pc_y[k] * lsi_1254[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, pc_y, lsh0_941, lsh0_942, lsh0_943, lsh1_941, \
                         lsh1_942, lsh1_943, lsi_1255, lsi_1256, \
                         lsi_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = f_10 * lsh0_941[k]
                    - f_11 * lsh1_941[k]
                    + f_3 * pc_y[k] * lsi_1255[k];

        t_1615[k] = f_8 * lsh0_942[k]
                    - f_9 * lsh1_942[k]
                    + f_3 * pc_y[k] * lsi_1256[k];

        t_1616[k] = f_6 * lsh0_943[k]
                    - f_7 * lsh1_943[k]
                    + f_3 * pc_y[k] * lsi_1257[k];
    }

#pragma omp simd aligned(t_1617, t_1618, t_1619, pc_y, pc_z, ksi_1007, lsh0_944, lsh1_944, \
                         lsi_1258, lsi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = f_4 * lsh0_944[k]
                    - f_5 * lsh1_944[k]
                    + f_3 * pc_y[k] * lsi_1258[k];

        t_1618[k] = f_3 * pc_y[k] * lsi_1259[k];

        t_1619[k] = f_0 * ksi_1007[k]
                    + f_1 * lsh0_944[k]
                    - f_2 * lsh1_944[k]
                    + f_3 * pc_z[k] * lsi_1259[k];
    }
}

auto
compute_prim_lsk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ksk0, const size_t ksi,
                                                   const size_t ksk1, const size_t lsh0,
                                                   const size_t lsh1, const size_t lsi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_lsk_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsh0, lsh1, lsi, ncols,
                                                              gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsh0, lsh1, lsi, ncols,
                                                              gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsh0, lsh1, lsi, ncols,
                                                              gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsh0, lsh1, lsi, ncols,
                                                              gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsh0, lsh1, lsi, ncols,
                                                              gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsh0, lsh1, lsi, ncols,
                                                              gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsh0, lsh1, lsi, ncols,
                                                              gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece7(buffer, target, pc, ksi, lsh0,
                                                              lsh1, lsi, ncols, gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsh0, lsh1, lsi, ncols,
                                                              gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, ksk0, ksi,
                                                              ksk1, lsi, ncols, gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, ksk0,
                                                               ksi, ksk1, lsh0, lsh1, lsi,
                                                               ncols, gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, ksk0,
                                                               ksi, ksk1, lsh0, lsh1, lsi,
                                                               ncols, gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece12(buffer, target, pc, ksi, lsh0,
                                                               lsh1, lsi, ncols, gamma, p, q);

    compute_prim_lsk_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, ksk0,
                                                               ksi, ksk1, lsh0, lsh1, lsi,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
