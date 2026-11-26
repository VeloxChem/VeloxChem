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


#include "SimdThreeCenterElectronRepulsionVrrRecLSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *ksl0_0 = buffer.data(ksl0 + 0);
    const auto *ksl0_3 = buffer.data(ksl0 + 3);
    const auto *ksl0_5 = buffer.data(ksl0 + 5);
    const auto *ksl0_6 = buffer.data(ksl0 + 6);
    const auto *ksl0_9 = buffer.data(ksl0 + 9);
    const auto *ksl0_10 = buffer.data(ksl0 + 10);
    const auto *ksl0_14 = buffer.data(ksl0 + 14);
    const auto *ksl0_15 = buffer.data(ksl0 + 15);
    const auto *ksl0_20 = buffer.data(ksl0 + 20);
    const auto *ksl0_21 = buffer.data(ksl0 + 21);
    const auto *ksl0_27 = buffer.data(ksl0 + 27);
    const auto *ksl0_36 = buffer.data(ksl0 + 36);
    const auto *ksl0_44 = buffer.data(ksl0 + 44);

    const auto *ksk_0 = buffer.data(ksk + 0);
    const auto *ksk_1 = buffer.data(ksk + 1);
    const auto *ksk_2 = buffer.data(ksk + 2);
    const auto *ksk_3 = buffer.data(ksk + 3);
    const auto *ksk_5 = buffer.data(ksk + 5);
    const auto *ksk_6 = buffer.data(ksk + 6);
    const auto *ksk_9 = buffer.data(ksk + 9);
    const auto *ksk_10 = buffer.data(ksk + 10);
    const auto *ksk_14 = buffer.data(ksk + 14);
    const auto *ksk_15 = buffer.data(ksk + 15);
    const auto *ksk_20 = buffer.data(ksk + 20);
    const auto *ksk_28 = buffer.data(ksk + 28);
    const auto *ksk_30 = buffer.data(ksk + 30);
    const auto *ksk_31 = buffer.data(ksk + 31);
    const auto *ksk_32 = buffer.data(ksk + 32);
    const auto *ksk_33 = buffer.data(ksk + 33);
    const auto *ksk_35 = buffer.data(ksk + 35);
    const auto *ksk_64 = buffer.data(ksk + 64);
    const auto *ksk_66 = buffer.data(ksk + 66);
    const auto *ksk_67 = buffer.data(ksk + 67);
    const auto *ksk_68 = buffer.data(ksk + 68);
    const auto *ksk_69 = buffer.data(ksk + 69);
    const auto *ksk_70 = buffer.data(ksk + 70);
    const auto *ksk_71 = buffer.data(ksk + 71);
    const auto *ksk_100 = buffer.data(ksk + 100);
    const auto *ksk_101 = buffer.data(ksk + 101);
    const auto *ksk_102 = buffer.data(ksk + 102);
    const auto *ksk_103 = buffer.data(ksk + 103);
    const auto *ksk_104 = buffer.data(ksk + 104);
    const auto *ksk_105 = buffer.data(ksk + 105);
    const auto *ksk_107 = buffer.data(ksk + 107);

    const auto *ksl1_0 = buffer.data(ksl1 + 0);
    const auto *ksl1_3 = buffer.data(ksl1 + 3);
    const auto *ksl1_5 = buffer.data(ksl1 + 5);
    const auto *ksl1_6 = buffer.data(ksl1 + 6);
    const auto *ksl1_9 = buffer.data(ksl1 + 9);
    const auto *ksl1_10 = buffer.data(ksl1 + 10);
    const auto *ksl1_14 = buffer.data(ksl1 + 14);
    const auto *ksl1_15 = buffer.data(ksl1 + 15);
    const auto *ksl1_20 = buffer.data(ksl1 + 20);
    const auto *ksl1_21 = buffer.data(ksl1 + 21);
    const auto *ksl1_27 = buffer.data(ksl1 + 27);
    const auto *ksl1_36 = buffer.data(ksl1 + 36);
    const auto *ksl1_44 = buffer.data(ksl1 + 44);

    const auto *lsi0_0 = buffer.data(lsi0 + 0);
    const auto *lsi0_1 = buffer.data(lsi0 + 1);
    const auto *lsi0_2 = buffer.data(lsi0 + 2);
    const auto *lsi0_3 = buffer.data(lsi0 + 3);
    const auto *lsi0_5 = buffer.data(lsi0 + 5);
    const auto *lsi0_6 = buffer.data(lsi0 + 6);
    const auto *lsi0_8 = buffer.data(lsi0 + 8);
    const auto *lsi0_9 = buffer.data(lsi0 + 9);
    const auto *lsi0_10 = buffer.data(lsi0 + 10);
    const auto *lsi0_12 = buffer.data(lsi0 + 12);
    const auto *lsi0_13 = buffer.data(lsi0 + 13);
    const auto *lsi0_14 = buffer.data(lsi0 + 14);
    const auto *lsi0_21 = buffer.data(lsi0 + 21);
    const auto *lsi0_23 = buffer.data(lsi0 + 23);
    const auto *lsi0_24 = buffer.data(lsi0 + 24);
    const auto *lsi0_25 = buffer.data(lsi0 + 25);
    const auto *lsi0_26 = buffer.data(lsi0 + 26);
    const auto *lsi0_27 = buffer.data(lsi0 + 27);
    const auto *lsi0_31 = buffer.data(lsi0 + 31);
    const auto *lsi0_34 = buffer.data(lsi0 + 34);
    const auto *lsi0_35 = buffer.data(lsi0 + 35);
    const auto *lsi0_38 = buffer.data(lsi0 + 38);
    const auto *lsi0_39 = buffer.data(lsi0 + 39);
    const auto *lsi0_40 = buffer.data(lsi0 + 40);
    const auto *lsi0_49 = buffer.data(lsi0 + 49);
    const auto *lsi0_50 = buffer.data(lsi0 + 50);
    const auto *lsi0_51 = buffer.data(lsi0 + 51);
    const auto *lsi0_52 = buffer.data(lsi0 + 52);
    const auto *lsi0_53 = buffer.data(lsi0 + 53);
    const auto *lsi0_58 = buffer.data(lsi0 + 58);
    const auto *lsi0_60 = buffer.data(lsi0 + 60);
    const auto *lsi0_61 = buffer.data(lsi0 + 61);
    const auto *lsi0_63 = buffer.data(lsi0 + 63);
    const auto *lsi0_64 = buffer.data(lsi0 + 64);
    const auto *lsi0_65 = buffer.data(lsi0 + 65);
    const auto *lsi0_67 = buffer.data(lsi0 + 67);
    const auto *lsi0_68 = buffer.data(lsi0 + 68);
    const auto *lsi0_69 = buffer.data(lsi0 + 69);
    const auto *lsi0_70 = buffer.data(lsi0 + 70);
    const auto *lsi0_78 = buffer.data(lsi0 + 78);
    const auto *lsi0_79 = buffer.data(lsi0 + 79);

    const auto *lsi1_0 = buffer.data(lsi1 + 0);
    const auto *lsi1_1 = buffer.data(lsi1 + 1);
    const auto *lsi1_2 = buffer.data(lsi1 + 2);
    const auto *lsi1_3 = buffer.data(lsi1 + 3);
    const auto *lsi1_5 = buffer.data(lsi1 + 5);
    const auto *lsi1_6 = buffer.data(lsi1 + 6);
    const auto *lsi1_8 = buffer.data(lsi1 + 8);
    const auto *lsi1_9 = buffer.data(lsi1 + 9);
    const auto *lsi1_10 = buffer.data(lsi1 + 10);
    const auto *lsi1_12 = buffer.data(lsi1 + 12);
    const auto *lsi1_13 = buffer.data(lsi1 + 13);
    const auto *lsi1_14 = buffer.data(lsi1 + 14);
    const auto *lsi1_21 = buffer.data(lsi1 + 21);
    const auto *lsi1_23 = buffer.data(lsi1 + 23);
    const auto *lsi1_24 = buffer.data(lsi1 + 24);
    const auto *lsi1_25 = buffer.data(lsi1 + 25);
    const auto *lsi1_26 = buffer.data(lsi1 + 26);
    const auto *lsi1_27 = buffer.data(lsi1 + 27);
    const auto *lsi1_31 = buffer.data(lsi1 + 31);
    const auto *lsi1_34 = buffer.data(lsi1 + 34);
    const auto *lsi1_35 = buffer.data(lsi1 + 35);
    const auto *lsi1_38 = buffer.data(lsi1 + 38);
    const auto *lsi1_39 = buffer.data(lsi1 + 39);
    const auto *lsi1_40 = buffer.data(lsi1 + 40);
    const auto *lsi1_49 = buffer.data(lsi1 + 49);
    const auto *lsi1_50 = buffer.data(lsi1 + 50);
    const auto *lsi1_51 = buffer.data(lsi1 + 51);
    const auto *lsi1_52 = buffer.data(lsi1 + 52);
    const auto *lsi1_53 = buffer.data(lsi1 + 53);
    const auto *lsi1_58 = buffer.data(lsi1 + 58);
    const auto *lsi1_60 = buffer.data(lsi1 + 60);
    const auto *lsi1_61 = buffer.data(lsi1 + 61);
    const auto *lsi1_63 = buffer.data(lsi1 + 63);
    const auto *lsi1_64 = buffer.data(lsi1 + 64);
    const auto *lsi1_65 = buffer.data(lsi1 + 65);
    const auto *lsi1_67 = buffer.data(lsi1 + 67);
    const auto *lsi1_68 = buffer.data(lsi1 + 68);
    const auto *lsi1_69 = buffer.data(lsi1 + 69);
    const auto *lsi1_70 = buffer.data(lsi1 + 70);
    const auto *lsi1_78 = buffer.data(lsi1 + 78);
    const auto *lsi1_79 = buffer.data(lsi1 + 79);

    const auto *lsk_0 = buffer.data(lsk + 0);
    const auto *lsk_1 = buffer.data(lsk + 1);
    const auto *lsk_2 = buffer.data(lsk + 2);
    const auto *lsk_3 = buffer.data(lsk + 3);
    const auto *lsk_5 = buffer.data(lsk + 5);
    const auto *lsk_6 = buffer.data(lsk + 6);
    const auto *lsk_8 = buffer.data(lsk + 8);
    const auto *lsk_9 = buffer.data(lsk + 9);
    const auto *lsk_10 = buffer.data(lsk + 10);
    const auto *lsk_12 = buffer.data(lsk + 12);
    const auto *lsk_13 = buffer.data(lsk + 13);
    const auto *lsk_14 = buffer.data(lsk + 14);
    const auto *lsk_15 = buffer.data(lsk + 15);
    const auto *lsk_17 = buffer.data(lsk + 17);
    const auto *lsk_18 = buffer.data(lsk + 18);
    const auto *lsk_19 = buffer.data(lsk + 19);
    const auto *lsk_20 = buffer.data(lsk + 20);
    const auto *lsk_21 = buffer.data(lsk + 21);
    const auto *lsk_27 = buffer.data(lsk + 27);
    const auto *lsk_28 = buffer.data(lsk + 28);
    const auto *lsk_30 = buffer.data(lsk + 30);
    const auto *lsk_31 = buffer.data(lsk + 31);
    const auto *lsk_32 = buffer.data(lsk + 32);
    const auto *lsk_33 = buffer.data(lsk + 33);
    const auto *lsk_34 = buffer.data(lsk + 34);
    const auto *lsk_35 = buffer.data(lsk + 35);
    const auto *lsk_36 = buffer.data(lsk + 36);
    const auto *lsk_37 = buffer.data(lsk + 37);
    const auto *lsk_39 = buffer.data(lsk + 39);
    const auto *lsk_41 = buffer.data(lsk + 41);
    const auto *lsk_42 = buffer.data(lsk + 42);
    const auto *lsk_43 = buffer.data(lsk + 43);
    const auto *lsk_45 = buffer.data(lsk + 45);
    const auto *lsk_46 = buffer.data(lsk + 46);
    const auto *lsk_47 = buffer.data(lsk + 47);
    const auto *lsk_48 = buffer.data(lsk + 48);
    const auto *lsk_50 = buffer.data(lsk + 50);
    const auto *lsk_51 = buffer.data(lsk + 51);
    const auto *lsk_52 = buffer.data(lsk + 52);
    const auto *lsk_53 = buffer.data(lsk + 53);
    const auto *lsk_54 = buffer.data(lsk + 54);
    const auto *lsk_56 = buffer.data(lsk + 56);
    const auto *lsk_57 = buffer.data(lsk + 57);
    const auto *lsk_64 = buffer.data(lsk + 64);
    const auto *lsk_65 = buffer.data(lsk + 65);
    const auto *lsk_66 = buffer.data(lsk + 66);
    const auto *lsk_67 = buffer.data(lsk + 67);
    const auto *lsk_68 = buffer.data(lsk + 68);
    const auto *lsk_69 = buffer.data(lsk + 69);
    const auto *lsk_70 = buffer.data(lsk + 70);
    const auto *lsk_71 = buffer.data(lsk + 71);
    const auto *lsk_72 = buffer.data(lsk + 72);
    const auto *lsk_74 = buffer.data(lsk + 74);
    const auto *lsk_76 = buffer.data(lsk + 76);
    const auto *lsk_77 = buffer.data(lsk + 77);
    const auto *lsk_79 = buffer.data(lsk + 79);
    const auto *lsk_80 = buffer.data(lsk + 80);
    const auto *lsk_81 = buffer.data(lsk + 81);
    const auto *lsk_83 = buffer.data(lsk + 83);
    const auto *lsk_84 = buffer.data(lsk + 84);
    const auto *lsk_85 = buffer.data(lsk + 85);
    const auto *lsk_86 = buffer.data(lsk + 86);
    const auto *lsk_88 = buffer.data(lsk + 88);
    const auto *lsk_89 = buffer.data(lsk + 89);
    const auto *lsk_90 = buffer.data(lsk + 90);
    const auto *lsk_91 = buffer.data(lsk + 91);
    const auto *lsk_92 = buffer.data(lsk + 92);
    const auto *lsk_99 = buffer.data(lsk + 99);
    const auto *lsk_100 = buffer.data(lsk + 100);
    const auto *lsk_101 = buffer.data(lsk + 101);
    const auto *lsk_102 = buffer.data(lsk + 102);
    const auto *lsk_103 = buffer.data(lsk + 103);
    const auto *lsk_104 = buffer.data(lsk + 104);
    const auto *lsk_105 = buffer.data(lsk + 105);
    const auto *lsk_107 = buffer.data(lsk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ksk_0, lsi0_0, \
                         lsi1_0, lsk_0, lsk_1, lsk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ksk_0[k]
                 + f_1 * lsi0_0[k]
                 - f_2 * lsi1_0[k]
                 + f_3 * pc_x[k] * lsk_0[k];

        t_1[k] = f_3 * pc_y[k] * lsk_0[k];

        t_2[k] = f_3 * pc_z[k] * lsk_0[k];

        t_3[k] = f_4 * lsi0_0[k]
                 - f_5 * lsi1_0[k]
                 + f_3 * pc_y[k] * lsk_1[k];

        t_4[k] = f_3 * pc_y[k] * lsk_2[k];

        t_5[k] = f_4 * lsi0_0[k]
                 - f_5 * lsi1_0[k]
                 + f_3 * pc_z[k] * lsk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, lsi0_1, lsi0_2, lsi0_3, lsi1_1, \
                         lsi1_2, lsi1_3, lsk_3, lsk_5, lsk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * lsi0_1[k]
                 - f_7 * lsi1_1[k]
                 + f_3 * pc_y[k] * lsk_3[k];

        t_7[k] = f_3 * pc_z[k] * lsk_3[k];

        t_8[k] = f_3 * pc_y[k] * lsk_5[k];

        t_9[k] = f_6 * lsi0_2[k]
                 - f_7 * lsi1_2[k]
                 + f_3 * pc_z[k] * lsk_5[k];

        t_10[k] = f_8 * lsi0_3[k]
                  - f_9 * lsi1_3[k]
                  + f_3 * pc_y[k] * lsk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, lsi0_5, lsi0_6, \
                         lsi1_5, lsi1_6, lsk_6, lsk_8, lsk_9, lsk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * lsk_6[k];

        t_12[k] = f_4 * lsi0_5[k]
                  - f_5 * lsi1_5[k]
                  + f_3 * pc_y[k] * lsk_8[k];

        t_13[k] = f_3 * pc_y[k] * lsk_9[k];

        t_14[k] = f_8 * lsi0_5[k]
                  - f_9 * lsi1_5[k]
                  + f_3 * pc_z[k] * lsk_9[k];

        t_15[k] = f_10 * lsi0_6[k]
                  - f_11 * lsi1_6[k]
                  + f_3 * pc_y[k] * lsk_10[k];

        t_16[k] = f_3 * pc_z[k] * lsk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, lsi0_8, lsi0_9, lsi1_8, lsi1_9, \
                         lsk_12, lsk_13, lsk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * lsi0_8[k]
                  - f_7 * lsi1_8[k]
                  + f_3 * pc_y[k] * lsk_12[k];

        t_18[k] = f_4 * lsi0_9[k]
                  - f_5 * lsi1_9[k]
                  + f_3 * pc_y[k] * lsk_13[k];

        t_19[k] = f_3 * pc_y[k] * lsk_14[k];

        t_20[k] = f_10 * lsi0_9[k]
                  - f_11 * lsi1_9[k]
                  + f_3 * pc_z[k] * lsk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, lsi0_10, lsi0_12, lsi0_13, \
                         lsi1_10, lsi1_12, lsi1_13, lsk_15, lsk_17, \
                         lsk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_12 * lsi0_10[k]
                  - f_13 * lsi1_10[k]
                  + f_3 * pc_y[k] * lsk_15[k];

        t_22[k] = f_3 * pc_z[k] * lsk_15[k];

        t_23[k] = f_8 * lsi0_12[k]
                  - f_9 * lsi1_12[k]
                  + f_3 * pc_y[k] * lsk_17[k];

        t_24[k] = f_6 * lsi0_13[k]
                  - f_7 * lsi1_13[k]
                  + f_3 * pc_y[k] * lsk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, ksk_28, lsi0_14, \
                         lsi1_14, lsk_19, lsk_20, lsk_21, lsk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * lsi0_14[k]
                  - f_5 * lsi1_14[k]
                  + f_3 * pc_y[k] * lsk_19[k];

        t_26[k] = f_3 * pc_y[k] * lsk_20[k];

        t_27[k] = f_12 * lsi0_14[k]
                  - f_13 * lsi1_14[k]
                  + f_3 * pc_z[k] * lsk_20[k];

        t_28[k] = f_0 * ksk_28[k]
                  + f_3 * pc_x[k] * lsk_28[k];

        t_29[k] = f_3 * pc_z[k] * lsk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, ksk_30, ksk_31, ksk_32, \
                         ksk_33, lsk_27, lsk_30, lsk_31, lsk_32, \
                         lsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * ksk_30[k]
                  + f_3 * pc_x[k] * lsk_30[k];

        t_31[k] = f_0 * ksk_31[k]
                  + f_3 * pc_x[k] * lsk_31[k];

        t_32[k] = f_0 * ksk_32[k]
                  + f_3 * pc_x[k] * lsk_32[k];

        t_33[k] = f_0 * ksk_33[k]
                  + f_3 * pc_x[k] * lsk_33[k];

        t_34[k] = f_3 * pc_y[k] * lsk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, ksk_35, lsi0_21, lsi0_23, \
                         lsi1_21, lsi1_23, lsk_28, lsk_30, lsk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * ksk_35[k]
                  + f_3 * pc_x[k] * lsk_35[k];

        t_36[k] = f_1 * lsi0_21[k]
                  - f_2 * lsi1_21[k]
                  + f_3 * pc_y[k] * lsk_28[k];

        t_37[k] = f_3 * pc_z[k] * lsk_28[k];

        t_38[k] = f_12 * lsi0_23[k]
                  - f_13 * lsi1_23[k]
                  + f_3 * pc_y[k] * lsk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pc_y, lsi0_24, lsi0_25, lsi0_26, lsi1_24, lsi1_25, \
                         lsi1_26, lsk_31, lsk_32, lsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * lsi0_24[k]
                  - f_11 * lsi1_24[k]
                  + f_3 * pc_y[k] * lsk_31[k];

        t_40[k] = f_8 * lsi0_25[k]
                  - f_9 * lsi1_25[k]
                  + f_3 * pc_y[k] * lsk_32[k];

        t_41[k] = f_6 * lsi0_26[k]
                  - f_7 * lsi1_26[k]
                  + f_3 * pc_y[k] * lsk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_y, pc_z, ksl0_0, ksk_0, \
                         ksl1_0, lsi0_27, lsi1_27, lsk_34, lsk_35, \
                         lsk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * lsi0_27[k]
                  - f_5 * lsi1_27[k]
                  + f_3 * pc_y[k] * lsk_34[k];

        t_43[k] = f_3 * pc_y[k] * lsk_35[k];

        t_44[k] = f_1 * lsi0_27[k]
                  - f_2 * lsi1_27[k]
                  + f_3 * pc_z[k] * lsk_35[k];

        t_45[k] = pa_y[k] * ksl0_0[k]
                  - f_14 * pc_y[k] * ksl1_0[k];

        t_46[k] = f_15 * ksk_0[k]
                  + f_3 * pc_y[k] * lsk_36[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_y, pc_z, ksl0_3, ksl0_5, ksk_1, \
                         ksl1_3, ksl1_5, lsk_36, lsk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * lsk_36[k];

        t_48[k] = pa_y[k] * ksl0_3[k]
                  + f_16 * ksk_1[k]
                  - f_14 * pc_y[k] * ksl1_3[k];

        t_49[k] = f_3 * pc_z[k] * lsk_37[k];

        t_50[k] = pa_y[k] * ksl0_5[k]
                  - f_14 * pc_y[k] * ksl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pc_y, pc_z, ksl0_6, ksl0_9, ksk_3, \
                         ksk_5, ksl1_6, ksl1_9, lsk_39, lsk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * ksl0_6[k]
                  + f_17 * ksk_3[k]
                  - f_14 * pc_y[k] * ksl1_6[k];

        t_52[k] = f_3 * pc_z[k] * lsk_39[k];

        t_53[k] = f_15 * ksk_5[k]
                  + f_3 * pc_y[k] * lsk_41[k];

        t_54[k] = pa_y[k] * ksl0_9[k]
                  - f_14 * pc_y[k] * ksl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pc_y, pc_z, ksl0_10, ksk_6, ksk_9, \
                         ksl1_10, lsi0_31, lsi1_31, lsk_42, lsk_43, \
                         lsk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * ksl0_10[k]
                  + f_18 * ksk_6[k]
                  - f_14 * pc_y[k] * ksl1_10[k];

        t_56[k] = f_3 * pc_z[k] * lsk_42[k];

        t_57[k] = f_4 * lsi0_31[k]
                  - f_5 * lsi1_31[k]
                  + f_3 * pc_z[k] * lsk_43[k];

        t_58[k] = f_15 * ksk_9[k]
                  + f_3 * pc_y[k] * lsk_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, ksl0_14, ksl0_15, ksk_10, \
                         ksl1_14, ksl1_15, lsi0_34, lsi1_34, lsk_46, \
                         lsk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * ksl0_14[k]
                  - f_14 * pc_y[k] * ksl1_14[k];

        t_60[k] = pa_y[k] * ksl0_15[k]
                  + f_19 * ksk_10[k]
                  - f_14 * pc_y[k] * ksl1_15[k];

        t_61[k] = f_3 * pc_z[k] * lsk_46[k];

        t_62[k] = f_4 * lsi0_34[k]
                  - f_5 * lsi1_34[k]
                  + f_3 * pc_z[k] * lsk_47[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pc_y, pc_z, ksl0_20, ksk_14, ksl1_20, \
                         lsi0_35, lsi1_35, lsk_48, lsk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_6 * lsi0_35[k]
                  - f_7 * lsi1_35[k]
                  + f_3 * pc_z[k] * lsk_48[k];

        t_64[k] = f_15 * ksk_14[k]
                  + f_3 * pc_y[k] * lsk_50[k];

        t_65[k] = pa_y[k] * ksl0_20[k]
                  - f_14 * pc_y[k] * ksl1_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pc_y, pc_z, ksl0_21, ksk_15, ksl1_21, \
                         lsi0_38, lsi1_38, lsk_51, lsk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * ksl0_21[k]
                  + f_20 * ksk_15[k]
                  - f_14 * pc_y[k] * ksl1_21[k];

        t_67[k] = f_3 * pc_z[k] * lsk_51[k];

        t_68[k] = f_4 * lsi0_38[k]
                  - f_5 * lsi1_38[k]
                  + f_3 * pc_z[k] * lsk_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_y, pc_z, ksk_20, lsi0_39, lsi0_40, lsi1_39, \
                         lsi1_40, lsk_53, lsk_54, lsk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * lsi0_39[k]
                  - f_7 * lsi1_39[k]
                  + f_3 * pc_z[k] * lsk_53[k];

        t_70[k] = f_8 * lsi0_40[k]
                  - f_9 * lsi1_40[k]
                  + f_3 * pc_z[k] * lsk_54[k];

        t_71[k] = f_15 * ksk_20[k]
                  + f_3 * pc_y[k] * lsk_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, pc_z, ksl0_27, ksk_64, \
                         ksk_66, ksl1_27, lsk_57, lsk_64, lsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * ksl0_27[k]
                  - f_14 * pc_y[k] * ksl1_27[k];

        t_73[k] = f_21 * ksk_64[k]
                  + f_3 * pc_x[k] * lsk_64[k];

        t_74[k] = f_3 * pc_z[k] * lsk_57[k];

        t_75[k] = f_21 * ksk_66[k]
                  + f_3 * pc_x[k] * lsk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, ksk_67, ksk_68, ksk_69, ksk_70, \
                         ksk_71, lsk_67, lsk_68, lsk_69, lsk_70, \
                         lsk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_21 * ksk_67[k]
                  + f_3 * pc_x[k] * lsk_67[k];

        t_77[k] = f_21 * ksk_68[k]
                  + f_3 * pc_x[k] * lsk_68[k];

        t_78[k] = f_21 * ksk_69[k]
                  + f_3 * pc_x[k] * lsk_69[k];

        t_79[k] = f_21 * ksk_70[k]
                  + f_3 * pc_x[k] * lsk_70[k];

        t_80[k] = f_21 * ksk_71[k]
                  + f_3 * pc_x[k] * lsk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_y, pc_z, ksk_28, lsi0_49, lsi0_50, \
                         lsi1_49, lsi1_50, lsk_64, lsk_65, lsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_15 * ksk_28[k]
                  + f_1 * lsi0_49[k]
                  - f_2 * lsi1_49[k]
                  + f_3 * pc_y[k] * lsk_64[k];

        t_82[k] = f_3 * pc_z[k] * lsk_64[k];

        t_83[k] = f_4 * lsi0_49[k]
                  - f_5 * lsi1_49[k]
                  + f_3 * pc_z[k] * lsk_65[k];

        t_84[k] = f_6 * lsi0_50[k]
                  - f_7 * lsi1_50[k]
                  + f_3 * pc_z[k] * lsk_66[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_z, lsi0_51, lsi0_52, lsi0_53, lsi1_51, lsi1_52, \
                         lsi1_53, lsk_67, lsk_68, lsk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_8 * lsi0_51[k]
                  - f_9 * lsi1_51[k]
                  + f_3 * pc_z[k] * lsk_67[k];

        t_86[k] = f_10 * lsi0_52[k]
                  - f_11 * lsi1_52[k]
                  + f_3 * pc_z[k] * lsk_68[k];

        t_87[k] = f_12 * lsi0_53[k]
                  - f_13 * lsi1_53[k]
                  + f_3 * pc_z[k] * lsk_69[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pc_y, pc_z, ksl0_0, ksl0_44, \
                         ksk_35, ksl1_0, ksl1_44, lsk_71, lsk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * ksk_35[k]
                  + f_3 * pc_y[k] * lsk_71[k];

        t_89[k] = pa_y[k] * ksl0_44[k]
                  - f_14 * pc_y[k] * ksl1_44[k];

        t_90[k] = pa_z[k] * ksl0_0[k]
                  - f_14 * pc_z[k] * ksl1_0[k];

        t_91[k] = f_3 * pc_y[k] * lsk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_y, pc_z, ksl0_3, ksl0_5, ksk_0, \
                         ksk_2, ksl1_3, ksl1_5, lsk_72, lsk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * ksk_0[k]
                  + f_3 * pc_z[k] * lsk_72[k];

        t_93[k] = pa_z[k] * ksl0_3[k]
                  - f_14 * pc_z[k] * ksl1_3[k];

        t_94[k] = f_3 * pc_y[k] * lsk_74[k];

        t_95[k] = pa_z[k] * ksl0_5[k]
                  + f_16 * ksk_2[k]
                  - f_14 * pc_z[k] * ksl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pc_y, pc_z, ksl0_6, ksl0_9, ksk_5, \
                         ksl1_6, ksl1_9, lsi0_58, lsi1_58, lsk_76, \
                         lsk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * ksl0_6[k]
                  - f_14 * pc_z[k] * ksl1_6[k];

        t_97[k] = f_4 * lsi0_58[k]
                  - f_5 * lsi1_58[k]
                  + f_3 * pc_y[k] * lsk_76[k];

        t_98[k] = f_3 * pc_y[k] * lsk_77[k];

        t_99[k] = pa_z[k] * ksl0_9[k]
                  + f_17 * ksk_5[k]
                  - f_14 * pc_z[k] * ksl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pc_y, pc_z, ksl0_10, ksl1_10, \
                         lsi0_60, lsi0_61, lsi1_60, lsi1_61, lsk_79, lsk_80, \
                         lsk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * ksl0_10[k]
                   - f_14 * pc_z[k] * ksl1_10[k];

        t_101[k] = f_6 * lsi0_60[k]
                   - f_7 * lsi1_60[k]
                   + f_3 * pc_y[k] * lsk_79[k];

        t_102[k] = f_4 * lsi0_61[k]
                   - f_5 * lsi1_61[k]
                   + f_3 * pc_y[k] * lsk_80[k];

        t_103[k] = f_3 * pc_y[k] * lsk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pc_y, pc_z, ksl0_14, ksl0_15, ksk_9, \
                         ksl1_14, ksl1_15, lsi0_63, lsi1_63, lsk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_z[k] * ksl0_14[k]
                   + f_18 * ksk_9[k]
                   - f_14 * pc_z[k] * ksl1_14[k];

        t_105[k] = pa_z[k] * ksl0_15[k]
                   - f_14 * pc_z[k] * ksl1_15[k];

        t_106[k] = f_8 * lsi0_63[k]
                   - f_9 * lsi1_63[k]
                   + f_3 * pc_y[k] * lsk_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, lsi0_64, lsi0_65, lsi1_64, lsi1_65, \
                         lsk_84, lsk_85, lsk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * lsi0_64[k]
                   - f_7 * lsi1_64[k]
                   + f_3 * pc_y[k] * lsk_84[k];

        t_108[k] = f_4 * lsi0_65[k]
                   - f_5 * lsi1_65[k]
                   + f_3 * pc_y[k] * lsk_85[k];

        t_109[k] = f_3 * pc_y[k] * lsk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pc_y, pc_z, ksl0_20, ksl0_21, ksk_14, \
                         ksl1_20, ksl1_21, lsi0_67, lsi1_67, lsk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * ksl0_20[k]
                   + f_19 * ksk_14[k]
                   - f_14 * pc_z[k] * ksl1_20[k];

        t_111[k] = pa_z[k] * ksl0_21[k]
                   - f_14 * pc_z[k] * ksl1_21[k];

        t_112[k] = f_10 * lsi0_67[k]
                   - f_11 * lsi1_67[k]
                   + f_3 * pc_y[k] * lsk_88[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_y, lsi0_68, lsi0_69, lsi0_70, lsi1_68, \
                         lsi1_69, lsi1_70, lsk_89, lsk_90, lsk_91, \
                         lsk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * lsi0_68[k]
                   - f_9 * lsi1_68[k]
                   + f_3 * pc_y[k] * lsk_89[k];

        t_114[k] = f_6 * lsi0_69[k]
                   - f_7 * lsi1_69[k]
                   + f_3 * pc_y[k] * lsk_90[k];

        t_115[k] = f_4 * lsi0_70[k]
                   - f_5 * lsi1_70[k]
                   + f_3 * pc_y[k] * lsk_91[k];

        t_116[k] = f_3 * pc_y[k] * lsk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pc_x, pc_z, ksl0_27, ksk_20, \
                         ksk_100, ksk_101, ksk_102, ksl1_27, lsk_100, lsk_101, \
                         lsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_z[k] * ksl0_27[k]
                   + f_20 * ksk_20[k]
                   - f_14 * pc_z[k] * ksl1_27[k];

        t_118[k] = f_21 * ksk_100[k]
                   + f_3 * pc_x[k] * lsk_100[k];

        t_119[k] = f_21 * ksk_101[k]
                   + f_3 * pc_x[k] * lsk_101[k];

        t_120[k] = f_21 * ksk_102[k]
                   + f_3 * pc_x[k] * lsk_102[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, ksk_103, ksk_104, \
                         ksk_105, ksk_107, lsk_99, lsk_103, lsk_104, lsk_105, \
                         lsk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_21 * ksk_103[k]
                   + f_3 * pc_x[k] * lsk_103[k];

        t_122[k] = f_21 * ksk_104[k]
                   + f_3 * pc_x[k] * lsk_104[k];

        t_123[k] = f_21 * ksk_105[k]
                   + f_3 * pc_x[k] * lsk_105[k];

        t_124[k] = f_3 * pc_y[k] * lsk_99[k];

        t_125[k] = f_21 * ksk_107[k]
                   + f_3 * pc_x[k] * lsk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_z, pc_y, pc_z, ksl0_36, ksl1_36, lsi0_78, \
                         lsi0_79, lsi1_78, lsi1_79, lsk_101, lsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * ksl0_36[k]
                   - f_14 * pc_z[k] * ksl1_36[k];

        t_127[k] = f_22 * lsi0_78[k]
                   - f_23 * lsi1_78[k]
                   + f_3 * pc_y[k] * lsk_101[k];

        t_128[k] = f_12 * lsi0_79[k]
                   - f_13 * lsi1_79[k]
                   + f_3 * pc_y[k] * lsk_102[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_48 = buffer.data(ksl0 + 48);
    const auto *ksl0_51 = buffer.data(ksl0 + 51);
    const auto *ksl0_55 = buffer.data(ksl0 + 55);
    const auto *ksl0_60 = buffer.data(ksl0 + 60);
    const auto *ksl0_66 = buffer.data(ksl0 + 66);
    const auto *ksl0_81 = buffer.data(ksl0 + 81);
    const auto *ksl0_90 = buffer.data(ksl0 + 90);
    const auto *ksl0_95 = buffer.data(ksl0 + 95);
    const auto *ksl0_99 = buffer.data(ksl0 + 99);
    const auto *ksl0_102 = buffer.data(ksl0 + 102);
    const auto *ksl0_104 = buffer.data(ksl0 + 104);
    const auto *ksl0_107 = buffer.data(ksl0 + 107);
    const auto *ksl0_108 = buffer.data(ksl0 + 108);
    const auto *ksl0_110 = buffer.data(ksl0 + 110);
    const auto *ksl0_113 = buffer.data(ksl0 + 113);
    const auto *ksl0_114 = buffer.data(ksl0 + 114);
    const auto *ksl0_115 = buffer.data(ksl0 + 115);
    const auto *ksl0_117 = buffer.data(ksl0 + 117);
    const auto *ksl0_134 = buffer.data(ksl0 + 134);

    const auto *ksk_35 = buffer.data(ksk + 35);
    const auto *ksk_36 = buffer.data(ksk + 36);
    const auto *ksk_39 = buffer.data(ksk + 39);
    const auto *ksk_41 = buffer.data(ksk + 41);
    const auto *ksk_42 = buffer.data(ksk + 42);
    const auto *ksk_45 = buffer.data(ksk + 45);
    const auto *ksk_46 = buffer.data(ksk + 46);
    const auto *ksk_50 = buffer.data(ksk + 50);
    const auto *ksk_51 = buffer.data(ksk + 51);
    const auto *ksk_56 = buffer.data(ksk + 56);
    const auto *ksk_64 = buffer.data(ksk + 64);
    const auto *ksk_71 = buffer.data(ksk + 71);
    const auto *ksk_72 = buffer.data(ksk + 72);
    const auto *ksk_74 = buffer.data(ksk + 74);
    const auto *ksk_77 = buffer.data(ksk + 77);
    const auto *ksk_80 = buffer.data(ksk + 80);
    const auto *ksk_81 = buffer.data(ksk + 81);
    const auto *ksk_84 = buffer.data(ksk + 84);
    const auto *ksk_85 = buffer.data(ksk + 85);
    const auto *ksk_86 = buffer.data(ksk + 86);
    const auto *ksk_89 = buffer.data(ksk + 89);
    const auto *ksk_90 = buffer.data(ksk + 90);
    const auto *ksk_91 = buffer.data(ksk + 91);
    const auto *ksk_92 = buffer.data(ksk + 92);
    const auto *ksk_102 = buffer.data(ksk + 102);
    const auto *ksk_103 = buffer.data(ksk + 103);
    const auto *ksk_104 = buffer.data(ksk + 104);
    const auto *ksk_105 = buffer.data(ksk + 105);
    const auto *ksk_106 = buffer.data(ksk + 106);
    const auto *ksk_107 = buffer.data(ksk + 107);
    const auto *ksk_108 = buffer.data(ksk + 108);
    const auto *ksk_111 = buffer.data(ksk + 111);
    const auto *ksk_114 = buffer.data(ksk + 114);
    const auto *ksk_118 = buffer.data(ksk + 118);
    const auto *ksk_123 = buffer.data(ksk + 123);
    const auto *ksk_129 = buffer.data(ksk + 129);
    const auto *ksk_136 = buffer.data(ksk + 136);
    const auto *ksk_138 = buffer.data(ksk + 138);
    const auto *ksk_139 = buffer.data(ksk + 139);
    const auto *ksk_140 = buffer.data(ksk + 140);
    const auto *ksk_141 = buffer.data(ksk + 141);
    const auto *ksk_142 = buffer.data(ksk + 142);
    const auto *ksk_143 = buffer.data(ksk + 143);
    const auto *ksk_172 = buffer.data(ksk + 172);
    const auto *ksk_173 = buffer.data(ksk + 173);
    const auto *ksk_174 = buffer.data(ksk + 174);
    const auto *ksk_175 = buffer.data(ksk + 175);
    const auto *ksk_176 = buffer.data(ksk + 176);
    const auto *ksk_177 = buffer.data(ksk + 177);
    const auto *ksk_178 = buffer.data(ksk + 178);
    const auto *ksk_179 = buffer.data(ksk + 179);
    const auto *ksk_180 = buffer.data(ksk + 180);
    const auto *ksk_185 = buffer.data(ksk + 185);
    const auto *ksk_189 = buffer.data(ksk + 189);
    const auto *ksk_194 = buffer.data(ksk + 194);
    const auto *ksk_200 = buffer.data(ksk + 200);

    const auto *ksl1_48 = buffer.data(ksl1 + 48);
    const auto *ksl1_51 = buffer.data(ksl1 + 51);
    const auto *ksl1_55 = buffer.data(ksl1 + 55);
    const auto *ksl1_60 = buffer.data(ksl1 + 60);
    const auto *ksl1_66 = buffer.data(ksl1 + 66);
    const auto *ksl1_81 = buffer.data(ksl1 + 81);
    const auto *ksl1_90 = buffer.data(ksl1 + 90);
    const auto *ksl1_95 = buffer.data(ksl1 + 95);
    const auto *ksl1_99 = buffer.data(ksl1 + 99);
    const auto *ksl1_102 = buffer.data(ksl1 + 102);
    const auto *ksl1_104 = buffer.data(ksl1 + 104);
    const auto *ksl1_107 = buffer.data(ksl1 + 107);
    const auto *ksl1_108 = buffer.data(ksl1 + 108);
    const auto *ksl1_110 = buffer.data(ksl1 + 110);
    const auto *ksl1_113 = buffer.data(ksl1 + 113);
    const auto *ksl1_114 = buffer.data(ksl1 + 114);
    const auto *ksl1_115 = buffer.data(ksl1 + 115);
    const auto *ksl1_117 = buffer.data(ksl1 + 117);
    const auto *ksl1_134 = buffer.data(ksl1 + 134);

    const auto *lsi0_80 = buffer.data(lsi0 + 80);
    const auto *lsi0_81 = buffer.data(lsi0 + 81);
    const auto *lsi0_82 = buffer.data(lsi0 + 82);
    const auto *lsi0_83 = buffer.data(lsi0 + 83);
    const auto *lsi0_84 = buffer.data(lsi0 + 84);
    const auto *lsi0_86 = buffer.data(lsi0 + 86);
    const auto *lsi0_87 = buffer.data(lsi0 + 87);
    const auto *lsi0_89 = buffer.data(lsi0 + 89);
    const auto *lsi0_90 = buffer.data(lsi0 + 90);
    const auto *lsi0_91 = buffer.data(lsi0 + 91);
    const auto *lsi0_93 = buffer.data(lsi0 + 93);
    const auto *lsi0_94 = buffer.data(lsi0 + 94);
    const auto *lsi0_95 = buffer.data(lsi0 + 95);
    const auto *lsi0_96 = buffer.data(lsi0 + 96);
    const auto *lsi0_98 = buffer.data(lsi0 + 98);
    const auto *lsi0_99 = buffer.data(lsi0 + 99);
    const auto *lsi0_105 = buffer.data(lsi0 + 105);
    const auto *lsi0_106 = buffer.data(lsi0 + 106);
    const auto *lsi0_107 = buffer.data(lsi0 + 107);
    const auto *lsi0_108 = buffer.data(lsi0 + 108);
    const auto *lsi0_109 = buffer.data(lsi0 + 109);
    const auto *lsi0_111 = buffer.data(lsi0 + 111);
    const auto *lsi0_135 = buffer.data(lsi0 + 135);
    const auto *lsi0_136 = buffer.data(lsi0 + 136);
    const auto *lsi0_137 = buffer.data(lsi0 + 137);
    const auto *lsi0_138 = buffer.data(lsi0 + 138);
    const auto *lsi0_139 = buffer.data(lsi0 + 139);
    const auto *lsi0_140 = buffer.data(lsi0 + 140);
    const auto *lsi0_141 = buffer.data(lsi0 + 141);
    const auto *lsi0_142 = buffer.data(lsi0 + 142);
    const auto *lsi0_143 = buffer.data(lsi0 + 143);
    const auto *lsi0_144 = buffer.data(lsi0 + 144);
    const auto *lsi0_145 = buffer.data(lsi0 + 145);
    const auto *lsi0_146 = buffer.data(lsi0 + 146);
    const auto *lsi0_147 = buffer.data(lsi0 + 147);
    const auto *lsi0_148 = buffer.data(lsi0 + 148);
    const auto *lsi0_149 = buffer.data(lsi0 + 149);
    const auto *lsi0_154 = buffer.data(lsi0 + 154);
    const auto *lsi0_160 = buffer.data(lsi0 + 160);

    const auto *lsi1_80 = buffer.data(lsi1 + 80);
    const auto *lsi1_81 = buffer.data(lsi1 + 81);
    const auto *lsi1_82 = buffer.data(lsi1 + 82);
    const auto *lsi1_83 = buffer.data(lsi1 + 83);
    const auto *lsi1_84 = buffer.data(lsi1 + 84);
    const auto *lsi1_86 = buffer.data(lsi1 + 86);
    const auto *lsi1_87 = buffer.data(lsi1 + 87);
    const auto *lsi1_89 = buffer.data(lsi1 + 89);
    const auto *lsi1_90 = buffer.data(lsi1 + 90);
    const auto *lsi1_91 = buffer.data(lsi1 + 91);
    const auto *lsi1_93 = buffer.data(lsi1 + 93);
    const auto *lsi1_94 = buffer.data(lsi1 + 94);
    const auto *lsi1_95 = buffer.data(lsi1 + 95);
    const auto *lsi1_96 = buffer.data(lsi1 + 96);
    const auto *lsi1_98 = buffer.data(lsi1 + 98);
    const auto *lsi1_99 = buffer.data(lsi1 + 99);
    const auto *lsi1_105 = buffer.data(lsi1 + 105);
    const auto *lsi1_106 = buffer.data(lsi1 + 106);
    const auto *lsi1_107 = buffer.data(lsi1 + 107);
    const auto *lsi1_108 = buffer.data(lsi1 + 108);
    const auto *lsi1_109 = buffer.data(lsi1 + 109);
    const auto *lsi1_111 = buffer.data(lsi1 + 111);
    const auto *lsi1_135 = buffer.data(lsi1 + 135);
    const auto *lsi1_136 = buffer.data(lsi1 + 136);
    const auto *lsi1_137 = buffer.data(lsi1 + 137);
    const auto *lsi1_138 = buffer.data(lsi1 + 138);
    const auto *lsi1_139 = buffer.data(lsi1 + 139);
    const auto *lsi1_140 = buffer.data(lsi1 + 140);
    const auto *lsi1_141 = buffer.data(lsi1 + 141);
    const auto *lsi1_142 = buffer.data(lsi1 + 142);
    const auto *lsi1_143 = buffer.data(lsi1 + 143);
    const auto *lsi1_144 = buffer.data(lsi1 + 144);
    const auto *lsi1_145 = buffer.data(lsi1 + 145);
    const auto *lsi1_146 = buffer.data(lsi1 + 146);
    const auto *lsi1_147 = buffer.data(lsi1 + 147);
    const auto *lsi1_148 = buffer.data(lsi1 + 148);
    const auto *lsi1_149 = buffer.data(lsi1 + 149);
    const auto *lsi1_154 = buffer.data(lsi1 + 154);
    const auto *lsi1_160 = buffer.data(lsi1 + 160);

    const auto *lsk_103 = buffer.data(lsk + 103);
    const auto *lsk_104 = buffer.data(lsk + 104);
    const auto *lsk_105 = buffer.data(lsk + 105);
    const auto *lsk_106 = buffer.data(lsk + 106);
    const auto *lsk_107 = buffer.data(lsk + 107);
    const auto *lsk_108 = buffer.data(lsk + 108);
    const auto *lsk_109 = buffer.data(lsk + 109);
    const auto *lsk_110 = buffer.data(lsk + 110);
    const auto *lsk_111 = buffer.data(lsk + 111);
    const auto *lsk_113 = buffer.data(lsk + 113);
    const auto *lsk_114 = buffer.data(lsk + 114);
    const auto *lsk_115 = buffer.data(lsk + 115);
    const auto *lsk_117 = buffer.data(lsk + 117);
    const auto *lsk_118 = buffer.data(lsk + 118);
    const auto *lsk_119 = buffer.data(lsk + 119);
    const auto *lsk_120 = buffer.data(lsk + 120);
    const auto *lsk_122 = buffer.data(lsk + 122);
    const auto *lsk_123 = buffer.data(lsk + 123);
    const auto *lsk_124 = buffer.data(lsk + 124);
    const auto *lsk_125 = buffer.data(lsk + 125);
    const auto *lsk_126 = buffer.data(lsk + 126);
    const auto *lsk_128 = buffer.data(lsk + 128);
    const auto *lsk_129 = buffer.data(lsk + 129);
    const auto *lsk_136 = buffer.data(lsk + 136);
    const auto *lsk_137 = buffer.data(lsk + 137);
    const auto *lsk_138 = buffer.data(lsk + 138);
    const auto *lsk_139 = buffer.data(lsk + 139);
    const auto *lsk_140 = buffer.data(lsk + 140);
    const auto *lsk_141 = buffer.data(lsk + 141);
    const auto *lsk_142 = buffer.data(lsk + 142);
    const auto *lsk_143 = buffer.data(lsk + 143);
    const auto *lsk_144 = buffer.data(lsk + 144);
    const auto *lsk_146 = buffer.data(lsk + 146);
    const auto *lsk_147 = buffer.data(lsk + 147);
    const auto *lsk_149 = buffer.data(lsk + 149);
    const auto *lsk_150 = buffer.data(lsk + 150);
    const auto *lsk_153 = buffer.data(lsk + 153);
    const auto *lsk_154 = buffer.data(lsk + 154);
    const auto *lsk_158 = buffer.data(lsk + 158);
    const auto *lsk_159 = buffer.data(lsk + 159);
    const auto *lsk_164 = buffer.data(lsk + 164);
    const auto *lsk_172 = buffer.data(lsk + 172);
    const auto *lsk_173 = buffer.data(lsk + 173);
    const auto *lsk_174 = buffer.data(lsk + 174);
    const auto *lsk_175 = buffer.data(lsk + 175);
    const auto *lsk_176 = buffer.data(lsk + 176);
    const auto *lsk_177 = buffer.data(lsk + 177);
    const auto *lsk_178 = buffer.data(lsk + 178);
    const auto *lsk_179 = buffer.data(lsk + 179);
    const auto *lsk_180 = buffer.data(lsk + 180);
    const auto *lsk_181 = buffer.data(lsk + 181);
    const auto *lsk_182 = buffer.data(lsk + 182);
    const auto *lsk_183 = buffer.data(lsk + 183);
    const auto *lsk_184 = buffer.data(lsk + 184);
    const auto *lsk_185 = buffer.data(lsk + 185);
    const auto *lsk_186 = buffer.data(lsk + 186);
    const auto *lsk_187 = buffer.data(lsk + 187);
    const auto *lsk_188 = buffer.data(lsk + 188);
    const auto *lsk_189 = buffer.data(lsk + 189);
    const auto *lsk_190 = buffer.data(lsk + 190);
    const auto *lsk_191 = buffer.data(lsk + 191);
    const auto *lsk_192 = buffer.data(lsk + 192);
    const auto *lsk_193 = buffer.data(lsk + 193);
    const auto *lsk_194 = buffer.data(lsk + 194);
    const auto *lsk_200 = buffer.data(lsk + 200);

#pragma omp simd aligned(t_129, t_130, t_131, pc_y, lsi0_80, lsi0_81, lsi0_82, lsi1_80, \
                         lsi1_81, lsi1_82, lsk_103, lsk_104, lsk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * lsi0_80[k]
                   - f_11 * lsi1_80[k]
                   + f_3 * pc_y[k] * lsk_103[k];

        t_130[k] = f_8 * lsi0_81[k]
                   - f_9 * lsi1_81[k]
                   + f_3 * pc_y[k] * lsk_104[k];

        t_131[k] = f_6 * lsi0_82[k]
                   - f_7 * lsi1_82[k]
                   + f_3 * pc_y[k] * lsk_105[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, pc_y, pc_z, ksk_35, ksk_108, \
                         lsi0_83, lsi0_84, lsi1_83, lsi1_84, lsk_106, lsk_107, \
                         lsk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * lsi0_83[k]
                   - f_5 * lsi1_83[k]
                   + f_3 * pc_y[k] * lsk_106[k];

        t_133[k] = f_3 * pc_y[k] * lsk_107[k];

        t_134[k] = f_15 * ksk_35[k]
                   + f_1 * lsi0_83[k]
                   - f_2 * lsi1_83[k]
                   + f_3 * pc_z[k] * lsk_107[k];

        t_135[k] = f_20 * ksk_108[k]
                   + f_1 * lsi0_84[k]
                   - f_2 * lsi1_84[k]
                   + f_3 * pc_x[k] * lsk_108[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, ksk_36, ksk_111, \
                         lsi0_87, lsi1_87, lsk_108, lsk_109, lsk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_16 * ksk_36[k]
                   + f_3 * pc_y[k] * lsk_108[k];

        t_137[k] = f_3 * pc_z[k] * lsk_108[k];

        t_138[k] = f_20 * ksk_111[k]
                   + f_12 * lsi0_87[k]
                   - f_13 * lsi1_87[k]
                   + f_3 * pc_x[k] * lsk_111[k];

        t_139[k] = f_3 * pc_z[k] * lsk_109[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pc_x, pc_z, ksk_114, lsi0_84, lsi0_90, lsi1_84, \
                         lsi1_90, lsk_110, lsk_111, lsk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * lsi0_84[k]
                   - f_5 * lsi1_84[k]
                   + f_3 * pc_z[k] * lsk_110[k];

        t_141[k] = f_20 * ksk_114[k]
                   + f_10 * lsi0_90[k]
                   - f_11 * lsi1_90[k]
                   + f_3 * pc_x[k] * lsk_114[k];

        t_142[k] = f_3 * pc_z[k] * lsk_111[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, ksk_41, ksk_118, \
                         lsi0_86, lsi0_94, lsi1_86, lsi1_94, lsk_113, lsk_114, \
                         lsk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * ksk_41[k]
                   + f_3 * pc_y[k] * lsk_113[k];

        t_144[k] = f_6 * lsi0_86[k]
                   - f_7 * lsi1_86[k]
                   + f_3 * pc_z[k] * lsk_113[k];

        t_145[k] = f_20 * ksk_118[k]
                   + f_8 * lsi0_94[k]
                   - f_9 * lsi1_94[k]
                   + f_3 * pc_x[k] * lsk_118[k];

        t_146[k] = f_3 * pc_z[k] * lsk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_y, pc_z, ksk_45, lsi0_87, lsi0_89, lsi1_87, \
                         lsi1_89, lsk_115, lsk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * lsi0_87[k]
                   - f_5 * lsi1_87[k]
                   + f_3 * pc_z[k] * lsk_115[k];

        t_148[k] = f_16 * ksk_45[k]
                   + f_3 * pc_y[k] * lsk_117[k];

        t_149[k] = f_8 * lsi0_89[k]
                   - f_9 * lsi1_89[k]
                   + f_3 * pc_z[k] * lsk_117[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, ksk_123, lsi0_90, lsi0_99, lsi1_90, \
                         lsi1_99, lsk_118, lsk_119, lsk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_20 * ksk_123[k]
                   + f_6 * lsi0_99[k]
                   - f_7 * lsi1_99[k]
                   + f_3 * pc_x[k] * lsk_123[k];

        t_151[k] = f_3 * pc_z[k] * lsk_118[k];

        t_152[k] = f_4 * lsi0_90[k]
                   - f_5 * lsi1_90[k]
                   + f_3 * pc_z[k] * lsk_119[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_y, pc_z, ksk_50, lsi0_91, lsi0_93, lsi1_91, \
                         lsi1_93, lsk_120, lsk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * lsi0_91[k]
                   - f_7 * lsi1_91[k]
                   + f_3 * pc_z[k] * lsk_120[k];

        t_154[k] = f_16 * ksk_50[k]
                   + f_3 * pc_y[k] * lsk_122[k];

        t_155[k] = f_10 * lsi0_93[k]
                   - f_11 * lsi1_93[k]
                   + f_3 * pc_z[k] * lsk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, ksk_129, lsi0_94, lsi0_105, lsi1_94, \
                         lsi1_105, lsk_123, lsk_124, lsk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_20 * ksk_129[k]
                   + f_4 * lsi0_105[k]
                   - f_5 * lsi1_105[k]
                   + f_3 * pc_x[k] * lsk_129[k];

        t_157[k] = f_3 * pc_z[k] * lsk_123[k];

        t_158[k] = f_4 * lsi0_94[k]
                   - f_5 * lsi1_94[k]
                   + f_3 * pc_z[k] * lsk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pc_y, pc_z, ksk_56, lsi0_95, lsi0_96, \
                         lsi0_98, lsi1_95, lsi1_96, lsi1_98, lsk_125, lsk_126, \
                         lsk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_6 * lsi0_95[k]
                   - f_7 * lsi1_95[k]
                   + f_3 * pc_z[k] * lsk_125[k];

        t_160[k] = f_8 * lsi0_96[k]
                   - f_9 * lsi1_96[k]
                   + f_3 * pc_z[k] * lsk_126[k];

        t_161[k] = f_16 * ksk_56[k]
                   + f_3 * pc_y[k] * lsk_128[k];

        t_162[k] = f_12 * lsi0_98[k]
                   - f_13 * lsi1_98[k]
                   + f_3 * pc_z[k] * lsk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pc_x, pc_z, ksk_136, ksk_138, \
                         ksk_139, ksk_140, lsk_129, lsk_136, lsk_138, lsk_139, \
                         lsk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_20 * ksk_136[k]
                   + f_3 * pc_x[k] * lsk_136[k];

        t_164[k] = f_3 * pc_z[k] * lsk_129[k];

        t_165[k] = f_20 * ksk_138[k]
                   + f_3 * pc_x[k] * lsk_138[k];

        t_166[k] = f_20 * ksk_139[k]
                   + f_3 * pc_x[k] * lsk_139[k];

        t_167[k] = f_20 * ksk_140[k]
                   + f_3 * pc_x[k] * lsk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, ksk_64, ksk_141, ksk_142, \
                         ksk_143, lsi0_105, lsi1_105, lsk_136, lsk_141, lsk_142, \
                         lsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_20 * ksk_141[k]
                   + f_3 * pc_x[k] * lsk_141[k];

        t_169[k] = f_20 * ksk_142[k]
                   + f_3 * pc_x[k] * lsk_142[k];

        t_170[k] = f_20 * ksk_143[k]
                   + f_3 * pc_x[k] * lsk_143[k];

        t_171[k] = f_16 * ksk_64[k]
                   + f_1 * lsi0_105[k]
                   - f_2 * lsi1_105[k]
                   + f_3 * pc_y[k] * lsk_136[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_z, lsi0_105, lsi0_106, lsi0_107, \
                         lsi1_105, lsi1_106, lsi1_107, lsk_136, lsk_137, lsk_138, \
                         lsk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * lsk_136[k];

        t_173[k] = f_4 * lsi0_105[k]
                   - f_5 * lsi1_105[k]
                   + f_3 * pc_z[k] * lsk_137[k];

        t_174[k] = f_6 * lsi0_106[k]
                   - f_7 * lsi1_106[k]
                   + f_3 * pc_z[k] * lsk_138[k];

        t_175[k] = f_8 * lsi0_107[k]
                   - f_9 * lsi1_107[k]
                   + f_3 * pc_z[k] * lsk_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_y, pc_z, ksk_71, lsi0_108, lsi0_109, \
                         lsi0_111, lsi1_108, lsi1_109, lsi1_111, lsk_140, lsk_141, \
                         lsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * lsi0_108[k]
                   - f_11 * lsi1_108[k]
                   + f_3 * pc_z[k] * lsk_140[k];

        t_177[k] = f_12 * lsi0_109[k]
                   - f_13 * lsi1_109[k]
                   + f_3 * pc_z[k] * lsk_141[k];

        t_178[k] = f_16 * ksk_71[k]
                   + f_3 * pc_y[k] * lsk_143[k];

        t_179[k] = f_1 * lsi0_111[k]
                   - f_2 * lsi1_111[k]
                   + f_3 * pc_z[k] * lsk_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, ksl0_48, ksl0_90, \
                         ksk_36, ksk_72, ksl1_48, ksl1_90, lsk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * ksl0_90[k]
                   - f_14 * pc_y[k] * ksl1_90[k];

        t_181[k] = f_15 * ksk_72[k]
                   + f_3 * pc_y[k] * lsk_144[k];

        t_182[k] = f_15 * ksk_36[k]
                   + f_3 * pc_z[k] * lsk_144[k];

        t_183[k] = pa_z[k] * ksl0_48[k]
                   - f_14 * pc_z[k] * ksl1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pa_z, pc_y, pc_z, ksl0_51, ksl0_95, \
                         ksk_39, ksk_74, ksl1_51, ksl1_95, lsk_146, \
                         lsk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_15 * ksk_74[k]
                   + f_3 * pc_y[k] * lsk_146[k];

        t_185[k] = pa_y[k] * ksl0_95[k]
                   - f_14 * pc_y[k] * ksl1_95[k];

        t_186[k] = pa_z[k] * ksl0_51[k]
                   - f_14 * pc_z[k] * ksl1_51[k];

        t_187[k] = f_15 * ksk_39[k]
                   + f_3 * pc_z[k] * lsk_147[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pc_y, pc_z, ksl0_55, ksl0_99, \
                         ksk_42, ksk_77, ksl1_55, ksl1_99, lsk_149, \
                         lsk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_15 * ksk_77[k]
                   + f_3 * pc_y[k] * lsk_149[k];

        t_189[k] = pa_y[k] * ksl0_99[k]
                   - f_14 * pc_y[k] * ksl1_99[k];

        t_190[k] = pa_z[k] * ksl0_55[k]
                   - f_14 * pc_z[k] * ksl1_55[k];

        t_191[k] = f_15 * ksk_42[k]
                   + f_3 * pc_z[k] * lsk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pc_y, ksl0_102, ksl0_104, ksk_80, ksk_81, \
                         ksl1_102, ksl1_104, lsk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_y[k] * ksl0_102[k]
                   + f_16 * ksk_80[k]
                   - f_14 * pc_y[k] * ksl1_102[k];

        t_193[k] = f_15 * ksk_81[k]
                   + f_3 * pc_y[k] * lsk_153[k];

        t_194[k] = pa_y[k] * ksl0_104[k]
                   - f_14 * pc_y[k] * ksl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_y, pa_z, pc_y, pc_z, ksl0_60, ksl0_107, \
                         ksk_46, ksk_84, ksl1_60, ksl1_107, lsk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * ksl0_60[k]
                   - f_14 * pc_z[k] * ksl1_60[k];

        t_196[k] = f_15 * ksk_46[k]
                   + f_3 * pc_z[k] * lsk_154[k];

        t_197[k] = pa_y[k] * ksl0_107[k]
                   + f_17 * ksk_84[k]
                   - f_14 * pc_y[k] * ksl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_y, pc_y, ksl0_108, ksl0_110, ksk_85, ksk_86, \
                         ksl1_108, ksl1_110, lsk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_y[k] * ksl0_108[k]
                   + f_16 * ksk_85[k]
                   - f_14 * pc_y[k] * ksl1_108[k];

        t_199[k] = f_15 * ksk_86[k]
                   + f_3 * pc_y[k] * lsk_158[k];

        t_200[k] = pa_y[k] * ksl0_110[k]
                   - f_14 * pc_y[k] * ksl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_y, pa_z, pc_y, pc_z, ksl0_66, ksl0_113, \
                         ksk_51, ksk_89, ksl1_66, ksl1_113, lsk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_z[k] * ksl0_66[k]
                   - f_14 * pc_z[k] * ksl1_66[k];

        t_202[k] = f_15 * ksk_51[k]
                   + f_3 * pc_z[k] * lsk_159[k];

        t_203[k] = pa_y[k] * ksl0_113[k]
                   + f_18 * ksk_89[k]
                   - f_14 * pc_y[k] * ksl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_y, ksl0_114, ksl0_115, ksl0_117, \
                         ksk_90, ksk_91, ksk_92, ksl1_114, ksl1_115, ksl1_117, \
                         lsk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * ksl0_114[k]
                   + f_17 * ksk_90[k]
                   - f_14 * pc_y[k] * ksl1_114[k];

        t_205[k] = pa_y[k] * ksl0_115[k]
                   + f_16 * ksk_91[k]
                   - f_14 * pc_y[k] * ksl1_115[k];

        t_206[k] = f_15 * ksk_92[k]
                   + f_3 * pc_y[k] * lsk_164[k];

        t_207[k] = pa_y[k] * ksl0_117[k]
                   - f_14 * pc_y[k] * ksl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, ksk_172, ksk_173, ksk_174, \
                         ksk_175, ksk_176, lsk_172, lsk_173, lsk_174, lsk_175, \
                         lsk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_20 * ksk_172[k]
                   + f_3 * pc_x[k] * lsk_172[k];

        t_209[k] = f_20 * ksk_173[k]
                   + f_3 * pc_x[k] * lsk_173[k];

        t_210[k] = f_20 * ksk_174[k]
                   + f_3 * pc_x[k] * lsk_174[k];

        t_211[k] = f_20 * ksk_175[k]
                   + f_3 * pc_x[k] * lsk_175[k];

        t_212[k] = f_20 * ksk_176[k]
                   + f_3 * pc_x[k] * lsk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_z, pc_x, pc_z, ksl0_81, ksk_177, \
                         ksk_178, ksk_179, ksl1_81, lsk_177, lsk_178, \
                         lsk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_20 * ksk_177[k]
                   + f_3 * pc_x[k] * lsk_177[k];

        t_214[k] = f_20 * ksk_178[k]
                   + f_3 * pc_x[k] * lsk_178[k];

        t_215[k] = f_20 * ksk_179[k]
                   + f_3 * pc_x[k] * lsk_179[k];

        t_216[k] = pa_z[k] * ksl0_81[k]
                   - f_14 * pc_z[k] * ksl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, ksk_64, ksk_102, ksk_103, lsi0_135, \
                         lsi0_136, lsi1_135, lsi1_136, lsk_172, lsk_174, \
                         lsk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * ksk_64[k]
                   + f_3 * pc_z[k] * lsk_172[k];

        t_218[k] = f_15 * ksk_102[k]
                   + f_12 * lsi0_135[k]
                   - f_13 * lsi1_135[k]
                   + f_3 * pc_y[k] * lsk_174[k];

        t_219[k] = f_15 * ksk_103[k]
                   + f_10 * lsi0_136[k]
                   - f_11 * lsi1_136[k]
                   + f_3 * pc_y[k] * lsk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, ksk_104, ksk_105, ksk_106, lsi0_137, \
                         lsi0_138, lsi0_139, lsi1_137, lsi1_138, lsi1_139, lsk_176, lsk_177, \
                         lsk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * ksk_104[k]
                   + f_8 * lsi0_137[k]
                   - f_9 * lsi1_137[k]
                   + f_3 * pc_y[k] * lsk_176[k];

        t_221[k] = f_15 * ksk_105[k]
                   + f_6 * lsi0_138[k]
                   - f_7 * lsi1_138[k]
                   + f_3 * pc_y[k] * lsk_177[k];

        t_222[k] = f_15 * ksk_106[k]
                   + f_4 * lsi0_139[k]
                   - f_5 * lsi1_139[k]
                   + f_3 * pc_y[k] * lsk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_y, pc_x, pc_y, ksl0_134, ksk_107, \
                         ksk_180, ksl1_134, lsi0_140, lsi1_140, lsk_179, \
                         lsk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * ksk_107[k]
                   + f_3 * pc_y[k] * lsk_179[k];

        t_224[k] = pa_y[k] * ksl0_134[k]
                   - f_14 * pc_y[k] * ksl1_134[k];

        t_225[k] = f_20 * ksk_180[k]
                   + f_1 * lsi0_140[k]
                   - f_2 * lsi1_140[k]
                   + f_3 * pc_x[k] * lsk_180[k];

        t_226[k] = f_3 * pc_y[k] * lsk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_y, pc_z, ksk_72, lsi0_140, lsi1_140, lsk_180, \
                         lsk_181, lsk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * ksk_72[k]
                   + f_3 * pc_z[k] * lsk_180[k];

        t_228[k] = f_4 * lsi0_140[k]
                   - f_5 * lsi1_140[k]
                   + f_3 * pc_y[k] * lsk_181[k];

        t_229[k] = f_3 * pc_y[k] * lsk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, ksk_185, lsi0_141, lsi0_142, \
                         lsi0_145, lsi1_141, lsi1_142, lsi1_145, lsk_183, lsk_184, \
                         lsk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_20 * ksk_185[k]
                   + f_12 * lsi0_145[k]
                   - f_13 * lsi1_145[k]
                   + f_3 * pc_x[k] * lsk_185[k];

        t_231[k] = f_6 * lsi0_141[k]
                   - f_7 * lsi1_141[k]
                   + f_3 * pc_y[k] * lsk_183[k];

        t_232[k] = f_4 * lsi0_142[k]
                   - f_5 * lsi1_142[k]
                   + f_3 * pc_y[k] * lsk_184[k];

        t_233[k] = f_3 * pc_y[k] * lsk_185[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, ksk_189, lsi0_143, lsi0_144, \
                         lsi0_149, lsi1_143, lsi1_144, lsi1_149, lsk_186, lsk_187, \
                         lsk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_20 * ksk_189[k]
                   + f_10 * lsi0_149[k]
                   - f_11 * lsi1_149[k]
                   + f_3 * pc_x[k] * lsk_189[k];

        t_235[k] = f_8 * lsi0_143[k]
                   - f_9 * lsi1_143[k]
                   + f_3 * pc_y[k] * lsk_186[k];

        t_236[k] = f_6 * lsi0_144[k]
                   - f_7 * lsi1_144[k]
                   + f_3 * pc_y[k] * lsk_187[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_x, pc_y, ksk_194, lsi0_145, lsi0_154, \
                         lsi1_145, lsi1_154, lsk_188, lsk_189, \
                         lsk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * lsi0_145[k]
                   - f_5 * lsi1_145[k]
                   + f_3 * pc_y[k] * lsk_188[k];

        t_238[k] = f_3 * pc_y[k] * lsk_189[k];

        t_239[k] = f_20 * ksk_194[k]
                   + f_8 * lsi0_154[k]
                   - f_9 * lsi1_154[k]
                   + f_3 * pc_x[k] * lsk_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pc_y, lsi0_146, lsi0_147, lsi0_148, lsi1_146, \
                         lsi1_147, lsi1_148, lsk_190, lsk_191, \
                         lsk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * lsi0_146[k]
                   - f_11 * lsi1_146[k]
                   + f_3 * pc_y[k] * lsk_190[k];

        t_241[k] = f_8 * lsi0_147[k]
                   - f_9 * lsi1_147[k]
                   + f_3 * pc_y[k] * lsk_191[k];

        t_242[k] = f_6 * lsi0_148[k]
                   - f_7 * lsi1_148[k]
                   + f_3 * pc_y[k] * lsk_192[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pc_x, pc_y, ksk_200, lsi0_149, lsi0_160, \
                         lsi1_149, lsi1_160, lsk_193, lsk_194, \
                         lsk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_4 * lsi0_149[k]
                   - f_5 * lsi1_149[k]
                   + f_3 * pc_y[k] * lsk_193[k];

        t_244[k] = f_3 * pc_y[k] * lsk_194[k];

        t_245[k] = f_20 * ksk_200[k]
                   + f_6 * lsi0_160[k]
                   - f_7 * lsi1_160[k]
                   + f_3 * pc_x[k] * lsk_200[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_135 = buffer.data(ksl0 + 135);
    const auto *ksl0_138 = buffer.data(ksl0 + 138);
    const auto *ksl0_141 = buffer.data(ksl0 + 141);
    const auto *ksl0_145 = buffer.data(ksl0 + 145);
    const auto *ksl0_147 = buffer.data(ksl0 + 147);
    const auto *ksl0_150 = buffer.data(ksl0 + 150);
    const auto *ksl0_152 = buffer.data(ksl0 + 152);
    const auto *ksl0_153 = buffer.data(ksl0 + 153);
    const auto *ksl0_156 = buffer.data(ksl0 + 156);
    const auto *ksl0_158 = buffer.data(ksl0 + 158);
    const auto *ksl0_159 = buffer.data(ksl0 + 159);
    const auto *ksl0_160 = buffer.data(ksl0 + 160);
    const auto *ksl0_171 = buffer.data(ksl0 + 171);
    const auto *ksl0_225 = buffer.data(ksl0 + 225);

    const auto *ksk_107 = buffer.data(ksk + 107);
    const auto *ksk_108 = buffer.data(ksk + 108);
    const auto *ksk_111 = buffer.data(ksk + 111);
    const auto *ksk_113 = buffer.data(ksk + 113);
    const auto *ksk_114 = buffer.data(ksk + 114);
    const auto *ksk_115 = buffer.data(ksk + 115);
    const auto *ksk_117 = buffer.data(ksk + 117);
    const auto *ksk_118 = buffer.data(ksk + 118);
    const auto *ksk_119 = buffer.data(ksk + 119);
    const auto *ksk_120 = buffer.data(ksk + 120);
    const auto *ksk_122 = buffer.data(ksk + 122);
    const auto *ksk_123 = buffer.data(ksk + 123);
    const auto *ksk_124 = buffer.data(ksk + 124);
    const auto *ksk_125 = buffer.data(ksk + 125);
    const auto *ksk_126 = buffer.data(ksk + 126);
    const auto *ksk_128 = buffer.data(ksk + 128);
    const auto *ksk_136 = buffer.data(ksk + 136);
    const auto *ksk_143 = buffer.data(ksk + 143);
    const auto *ksk_144 = buffer.data(ksk + 144);
    const auto *ksk_146 = buffer.data(ksk + 146);
    const auto *ksk_149 = buffer.data(ksk + 149);
    const auto *ksk_153 = buffer.data(ksk + 153);
    const auto *ksk_158 = buffer.data(ksk + 158);
    const auto *ksk_164 = buffer.data(ksk + 164);
    const auto *ksk_174 = buffer.data(ksk + 174);
    const auto *ksk_175 = buffer.data(ksk + 175);
    const auto *ksk_176 = buffer.data(ksk + 176);
    const auto *ksk_177 = buffer.data(ksk + 177);
    const auto *ksk_178 = buffer.data(ksk + 178);
    const auto *ksk_179 = buffer.data(ksk + 179);
    const auto *ksk_180 = buffer.data(ksk + 180);
    const auto *ksk_207 = buffer.data(ksk + 207);
    const auto *ksk_208 = buffer.data(ksk + 208);
    const auto *ksk_209 = buffer.data(ksk + 209);
    const auto *ksk_210 = buffer.data(ksk + 210);
    const auto *ksk_211 = buffer.data(ksk + 211);
    const auto *ksk_212 = buffer.data(ksk + 212);
    const auto *ksk_213 = buffer.data(ksk + 213);
    const auto *ksk_215 = buffer.data(ksk + 215);
    const auto *ksk_216 = buffer.data(ksk + 216);
    const auto *ksk_219 = buffer.data(ksk + 219);
    const auto *ksk_222 = buffer.data(ksk + 222);
    const auto *ksk_226 = buffer.data(ksk + 226);
    const auto *ksk_231 = buffer.data(ksk + 231);
    const auto *ksk_237 = buffer.data(ksk + 237);
    const auto *ksk_244 = buffer.data(ksk + 244);
    const auto *ksk_246 = buffer.data(ksk + 246);
    const auto *ksk_247 = buffer.data(ksk + 247);
    const auto *ksk_248 = buffer.data(ksk + 248);
    const auto *ksk_249 = buffer.data(ksk + 249);
    const auto *ksk_250 = buffer.data(ksk + 250);
    const auto *ksk_251 = buffer.data(ksk + 251);
    const auto *ksk_257 = buffer.data(ksk + 257);
    const auto *ksk_261 = buffer.data(ksk + 261);
    const auto *ksk_266 = buffer.data(ksk + 266);
    const auto *ksk_272 = buffer.data(ksk + 272);
    const auto *ksk_279 = buffer.data(ksk + 279);
    const auto *ksk_280 = buffer.data(ksk + 280);
    const auto *ksk_281 = buffer.data(ksk + 281);
    const auto *ksk_282 = buffer.data(ksk + 282);
    const auto *ksk_283 = buffer.data(ksk + 283);
    const auto *ksk_284 = buffer.data(ksk + 284);
    const auto *ksk_285 = buffer.data(ksk + 285);
    const auto *ksk_286 = buffer.data(ksk + 286);
    const auto *ksk_287 = buffer.data(ksk + 287);

    const auto *ksl1_135 = buffer.data(ksl1 + 135);
    const auto *ksl1_138 = buffer.data(ksl1 + 138);
    const auto *ksl1_141 = buffer.data(ksl1 + 141);
    const auto *ksl1_145 = buffer.data(ksl1 + 145);
    const auto *ksl1_147 = buffer.data(ksl1 + 147);
    const auto *ksl1_150 = buffer.data(ksl1 + 150);
    const auto *ksl1_152 = buffer.data(ksl1 + 152);
    const auto *ksl1_153 = buffer.data(ksl1 + 153);
    const auto *ksl1_156 = buffer.data(ksl1 + 156);
    const auto *ksl1_158 = buffer.data(ksl1 + 158);
    const auto *ksl1_159 = buffer.data(ksl1 + 159);
    const auto *ksl1_160 = buffer.data(ksl1 + 160);
    const auto *ksl1_171 = buffer.data(ksl1 + 171);
    const auto *ksl1_225 = buffer.data(ksl1 + 225);

    const auto *lsi0_150 = buffer.data(lsi0 + 150);
    const auto *lsi0_151 = buffer.data(lsi0 + 151);
    const auto *lsi0_152 = buffer.data(lsi0 + 152);
    const auto *lsi0_153 = buffer.data(lsi0 + 153);
    const auto *lsi0_154 = buffer.data(lsi0 + 154);
    const auto *lsi0_161 = buffer.data(lsi0 + 161);
    const auto *lsi0_162 = buffer.data(lsi0 + 162);
    const auto *lsi0_163 = buffer.data(lsi0 + 163);
    const auto *lsi0_164 = buffer.data(lsi0 + 164);
    const auto *lsi0_165 = buffer.data(lsi0 + 165);
    const auto *lsi0_166 = buffer.data(lsi0 + 166);
    const auto *lsi0_167 = buffer.data(lsi0 + 167);
    const auto *lsi0_168 = buffer.data(lsi0 + 168);
    const auto *lsi0_170 = buffer.data(lsi0 + 170);
    const auto *lsi0_171 = buffer.data(lsi0 + 171);
    const auto *lsi0_173 = buffer.data(lsi0 + 173);
    const auto *lsi0_174 = buffer.data(lsi0 + 174);
    const auto *lsi0_175 = buffer.data(lsi0 + 175);
    const auto *lsi0_177 = buffer.data(lsi0 + 177);
    const auto *lsi0_178 = buffer.data(lsi0 + 178);
    const auto *lsi0_179 = buffer.data(lsi0 + 179);
    const auto *lsi0_180 = buffer.data(lsi0 + 180);
    const auto *lsi0_182 = buffer.data(lsi0 + 182);
    const auto *lsi0_183 = buffer.data(lsi0 + 183);
    const auto *lsi0_189 = buffer.data(lsi0 + 189);
    const auto *lsi0_190 = buffer.data(lsi0 + 190);
    const auto *lsi0_191 = buffer.data(lsi0 + 191);
    const auto *lsi0_192 = buffer.data(lsi0 + 192);
    const auto *lsi0_193 = buffer.data(lsi0 + 193);
    const auto *lsi0_195 = buffer.data(lsi0 + 195);
    const auto *lsi0_201 = buffer.data(lsi0 + 201);
    const auto *lsi0_205 = buffer.data(lsi0 + 205);
    const auto *lsi0_210 = buffer.data(lsi0 + 210);
    const auto *lsi0_216 = buffer.data(lsi0 + 216);
    const auto *lsi0_219 = buffer.data(lsi0 + 219);
    const auto *lsi0_220 = buffer.data(lsi0 + 220);
    const auto *lsi0_221 = buffer.data(lsi0 + 221);
    const auto *lsi0_222 = buffer.data(lsi0 + 222);
    const auto *lsi0_223 = buffer.data(lsi0 + 223);

    const auto *lsi1_150 = buffer.data(lsi1 + 150);
    const auto *lsi1_151 = buffer.data(lsi1 + 151);
    const auto *lsi1_152 = buffer.data(lsi1 + 152);
    const auto *lsi1_153 = buffer.data(lsi1 + 153);
    const auto *lsi1_154 = buffer.data(lsi1 + 154);
    const auto *lsi1_161 = buffer.data(lsi1 + 161);
    const auto *lsi1_162 = buffer.data(lsi1 + 162);
    const auto *lsi1_163 = buffer.data(lsi1 + 163);
    const auto *lsi1_164 = buffer.data(lsi1 + 164);
    const auto *lsi1_165 = buffer.data(lsi1 + 165);
    const auto *lsi1_166 = buffer.data(lsi1 + 166);
    const auto *lsi1_167 = buffer.data(lsi1 + 167);
    const auto *lsi1_168 = buffer.data(lsi1 + 168);
    const auto *lsi1_170 = buffer.data(lsi1 + 170);
    const auto *lsi1_171 = buffer.data(lsi1 + 171);
    const auto *lsi1_173 = buffer.data(lsi1 + 173);
    const auto *lsi1_174 = buffer.data(lsi1 + 174);
    const auto *lsi1_175 = buffer.data(lsi1 + 175);
    const auto *lsi1_177 = buffer.data(lsi1 + 177);
    const auto *lsi1_178 = buffer.data(lsi1 + 178);
    const auto *lsi1_179 = buffer.data(lsi1 + 179);
    const auto *lsi1_180 = buffer.data(lsi1 + 180);
    const auto *lsi1_182 = buffer.data(lsi1 + 182);
    const auto *lsi1_183 = buffer.data(lsi1 + 183);
    const auto *lsi1_189 = buffer.data(lsi1 + 189);
    const auto *lsi1_190 = buffer.data(lsi1 + 190);
    const auto *lsi1_191 = buffer.data(lsi1 + 191);
    const auto *lsi1_192 = buffer.data(lsi1 + 192);
    const auto *lsi1_193 = buffer.data(lsi1 + 193);
    const auto *lsi1_195 = buffer.data(lsi1 + 195);
    const auto *lsi1_201 = buffer.data(lsi1 + 201);
    const auto *lsi1_205 = buffer.data(lsi1 + 205);
    const auto *lsi1_210 = buffer.data(lsi1 + 210);
    const auto *lsi1_216 = buffer.data(lsi1 + 216);
    const auto *lsi1_219 = buffer.data(lsi1 + 219);
    const auto *lsi1_220 = buffer.data(lsi1 + 220);
    const auto *lsi1_221 = buffer.data(lsi1 + 221);
    const auto *lsi1_222 = buffer.data(lsi1 + 222);
    const auto *lsi1_223 = buffer.data(lsi1 + 223);

    const auto *lsk_195 = buffer.data(lsk + 195);
    const auto *lsk_196 = buffer.data(lsk + 196);
    const auto *lsk_197 = buffer.data(lsk + 197);
    const auto *lsk_198 = buffer.data(lsk + 198);
    const auto *lsk_199 = buffer.data(lsk + 199);
    const auto *lsk_200 = buffer.data(lsk + 200);
    const auto *lsk_207 = buffer.data(lsk + 207);
    const auto *lsk_208 = buffer.data(lsk + 208);
    const auto *lsk_209 = buffer.data(lsk + 209);
    const auto *lsk_210 = buffer.data(lsk + 210);
    const auto *lsk_211 = buffer.data(lsk + 211);
    const auto *lsk_212 = buffer.data(lsk + 212);
    const auto *lsk_213 = buffer.data(lsk + 213);
    const auto *lsk_214 = buffer.data(lsk + 214);
    const auto *lsk_215 = buffer.data(lsk + 215);
    const auto *lsk_216 = buffer.data(lsk + 216);
    const auto *lsk_217 = buffer.data(lsk + 217);
    const auto *lsk_218 = buffer.data(lsk + 218);
    const auto *lsk_219 = buffer.data(lsk + 219);
    const auto *lsk_221 = buffer.data(lsk + 221);
    const auto *lsk_222 = buffer.data(lsk + 222);
    const auto *lsk_223 = buffer.data(lsk + 223);
    const auto *lsk_225 = buffer.data(lsk + 225);
    const auto *lsk_226 = buffer.data(lsk + 226);
    const auto *lsk_227 = buffer.data(lsk + 227);
    const auto *lsk_228 = buffer.data(lsk + 228);
    const auto *lsk_230 = buffer.data(lsk + 230);
    const auto *lsk_231 = buffer.data(lsk + 231);
    const auto *lsk_232 = buffer.data(lsk + 232);
    const auto *lsk_233 = buffer.data(lsk + 233);
    const auto *lsk_234 = buffer.data(lsk + 234);
    const auto *lsk_236 = buffer.data(lsk + 236);
    const auto *lsk_237 = buffer.data(lsk + 237);
    const auto *lsk_244 = buffer.data(lsk + 244);
    const auto *lsk_245 = buffer.data(lsk + 245);
    const auto *lsk_246 = buffer.data(lsk + 246);
    const auto *lsk_247 = buffer.data(lsk + 247);
    const auto *lsk_248 = buffer.data(lsk + 248);
    const auto *lsk_249 = buffer.data(lsk + 249);
    const auto *lsk_250 = buffer.data(lsk + 250);
    const auto *lsk_251 = buffer.data(lsk + 251);
    const auto *lsk_252 = buffer.data(lsk + 252);
    const auto *lsk_254 = buffer.data(lsk + 254);
    const auto *lsk_255 = buffer.data(lsk + 255);
    const auto *lsk_257 = buffer.data(lsk + 257);
    const auto *lsk_258 = buffer.data(lsk + 258);
    const auto *lsk_261 = buffer.data(lsk + 261);
    const auto *lsk_262 = buffer.data(lsk + 262);
    const auto *lsk_266 = buffer.data(lsk + 266);
    const auto *lsk_267 = buffer.data(lsk + 267);
    const auto *lsk_272 = buffer.data(lsk + 272);
    const auto *lsk_279 = buffer.data(lsk + 279);
    const auto *lsk_280 = buffer.data(lsk + 280);
    const auto *lsk_281 = buffer.data(lsk + 281);
    const auto *lsk_282 = buffer.data(lsk + 282);
    const auto *lsk_283 = buffer.data(lsk + 283);
    const auto *lsk_284 = buffer.data(lsk + 284);
    const auto *lsk_285 = buffer.data(lsk + 285);
    const auto *lsk_286 = buffer.data(lsk + 286);
    const auto *lsk_287 = buffer.data(lsk + 287);
    const auto *lsk_288 = buffer.data(lsk + 288);

#pragma omp simd aligned(t_246, t_247, t_248, pc_y, lsi0_150, lsi0_151, lsi0_152, lsi1_150, \
                         lsi1_151, lsi1_152, lsk_195, lsk_196, \
                         lsk_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * lsi0_150[k]
                   - f_13 * lsi1_150[k]
                   + f_3 * pc_y[k] * lsk_195[k];

        t_247[k] = f_10 * lsi0_151[k]
                   - f_11 * lsi1_151[k]
                   + f_3 * pc_y[k] * lsk_196[k];

        t_248[k] = f_8 * lsi0_152[k]
                   - f_9 * lsi1_152[k]
                   + f_3 * pc_y[k] * lsk_197[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, lsi0_153, lsi0_154, lsi1_153, lsi1_154, \
                         lsk_198, lsk_199, lsk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * lsi0_153[k]
                   - f_7 * lsi1_153[k]
                   + f_3 * pc_y[k] * lsk_198[k];

        t_250[k] = f_4 * lsi0_154[k]
                   - f_5 * lsi1_154[k]
                   + f_3 * pc_y[k] * lsk_199[k];

        t_251[k] = f_3 * pc_y[k] * lsk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, ksk_207, ksk_208, ksk_209, ksk_210, \
                         lsi0_167, lsi1_167, lsk_207, lsk_208, lsk_209, \
                         lsk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_20 * ksk_207[k]
                   + f_4 * lsi0_167[k]
                   - f_5 * lsi1_167[k]
                   + f_3 * pc_x[k] * lsk_207[k];

        t_253[k] = f_20 * ksk_208[k]
                   + f_3 * pc_x[k] * lsk_208[k];

        t_254[k] = f_20 * ksk_209[k]
                   + f_3 * pc_x[k] * lsk_209[k];

        t_255[k] = f_20 * ksk_210[k]
                   + f_3 * pc_x[k] * lsk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pc_x, pc_y, ksk_211, ksk_212, \
                         ksk_213, ksk_215, lsk_207, lsk_211, lsk_212, lsk_213, \
                         lsk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_20 * ksk_211[k]
                   + f_3 * pc_x[k] * lsk_211[k];

        t_257[k] = f_20 * ksk_212[k]
                   + f_3 * pc_x[k] * lsk_212[k];

        t_258[k] = f_20 * ksk_213[k]
                   + f_3 * pc_x[k] * lsk_213[k];

        t_259[k] = f_3 * pc_y[k] * lsk_207[k];

        t_260[k] = f_20 * ksk_215[k]
                   + f_3 * pc_x[k] * lsk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_y, lsi0_161, lsi0_162, lsi0_163, lsi1_161, \
                         lsi1_162, lsi1_163, lsk_208, lsk_209, \
                         lsk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * lsi0_161[k]
                   - f_2 * lsi1_161[k]
                   + f_3 * pc_y[k] * lsk_208[k];

        t_262[k] = f_22 * lsi0_162[k]
                   - f_23 * lsi1_162[k]
                   + f_3 * pc_y[k] * lsk_209[k];

        t_263[k] = f_12 * lsi0_163[k]
                   - f_13 * lsi1_163[k]
                   + f_3 * pc_y[k] * lsk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, lsi0_164, lsi0_165, lsi0_166, lsi1_164, \
                         lsi1_165, lsi1_166, lsk_211, lsk_212, \
                         lsk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * lsi0_164[k]
                   - f_11 * lsi1_164[k]
                   + f_3 * pc_y[k] * lsk_211[k];

        t_265[k] = f_8 * lsi0_165[k]
                   - f_9 * lsi1_165[k]
                   + f_3 * pc_y[k] * lsk_212[k];

        t_266[k] = f_6 * lsi0_166[k]
                   - f_7 * lsi1_166[k]
                   + f_3 * pc_y[k] * lsk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, ksk_107, ksk_216, \
                         lsi0_167, lsi0_168, lsi1_167, lsi1_168, lsk_214, lsk_215, \
                         lsk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_4 * lsi0_167[k]
                   - f_5 * lsi1_167[k]
                   + f_3 * pc_y[k] * lsk_214[k];

        t_268[k] = f_3 * pc_y[k] * lsk_215[k];

        t_269[k] = f_16 * ksk_107[k]
                   + f_1 * lsi0_167[k]
                   - f_2 * lsi1_167[k]
                   + f_3 * pc_z[k] * lsk_215[k];

        t_270[k] = f_19 * ksk_216[k]
                   + f_1 * lsi0_168[k]
                   - f_2 * lsi1_168[k]
                   + f_3 * pc_x[k] * lsk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, ksk_108, ksk_219, \
                         lsi0_171, lsi1_171, lsk_216, lsk_217, \
                         lsk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * ksk_108[k]
                   + f_3 * pc_y[k] * lsk_216[k];

        t_272[k] = f_3 * pc_z[k] * lsk_216[k];

        t_273[k] = f_19 * ksk_219[k]
                   + f_12 * lsi0_171[k]
                   - f_13 * lsi1_171[k]
                   + f_3 * pc_x[k] * lsk_219[k];

        t_274[k] = f_3 * pc_z[k] * lsk_217[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, ksk_222, lsi0_168, lsi0_174, \
                         lsi1_168, lsi1_174, lsk_218, lsk_219, \
                         lsk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_4 * lsi0_168[k]
                   - f_5 * lsi1_168[k]
                   + f_3 * pc_z[k] * lsk_218[k];

        t_276[k] = f_19 * ksk_222[k]
                   + f_10 * lsi0_174[k]
                   - f_11 * lsi1_174[k]
                   + f_3 * pc_x[k] * lsk_222[k];

        t_277[k] = f_3 * pc_z[k] * lsk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pc_x, pc_y, pc_z, ksk_113, ksk_226, \
                         lsi0_170, lsi0_178, lsi1_170, lsi1_178, lsk_221, lsk_222, \
                         lsk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * ksk_113[k]
                   + f_3 * pc_y[k] * lsk_221[k];

        t_279[k] = f_6 * lsi0_170[k]
                   - f_7 * lsi1_170[k]
                   + f_3 * pc_z[k] * lsk_221[k];

        t_280[k] = f_19 * ksk_226[k]
                   + f_8 * lsi0_178[k]
                   - f_9 * lsi1_178[k]
                   + f_3 * pc_x[k] * lsk_226[k];

        t_281[k] = f_3 * pc_z[k] * lsk_222[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pc_y, pc_z, ksk_117, lsi0_171, lsi0_173, \
                         lsi1_171, lsi1_173, lsk_223, lsk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_4 * lsi0_171[k]
                   - f_5 * lsi1_171[k]
                   + f_3 * pc_z[k] * lsk_223[k];

        t_283[k] = f_17 * ksk_117[k]
                   + f_3 * pc_y[k] * lsk_225[k];

        t_284[k] = f_8 * lsi0_173[k]
                   - f_9 * lsi1_173[k]
                   + f_3 * pc_z[k] * lsk_225[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pc_x, pc_z, ksk_231, lsi0_174, lsi0_183, \
                         lsi1_174, lsi1_183, lsk_226, lsk_227, \
                         lsk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_19 * ksk_231[k]
                   + f_6 * lsi0_183[k]
                   - f_7 * lsi1_183[k]
                   + f_3 * pc_x[k] * lsk_231[k];

        t_286[k] = f_3 * pc_z[k] * lsk_226[k];

        t_287[k] = f_4 * lsi0_174[k]
                   - f_5 * lsi1_174[k]
                   + f_3 * pc_z[k] * lsk_227[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, ksk_122, lsi0_175, lsi0_177, \
                         lsi1_175, lsi1_177, lsk_228, lsk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_6 * lsi0_175[k]
                   - f_7 * lsi1_175[k]
                   + f_3 * pc_z[k] * lsk_228[k];

        t_289[k] = f_17 * ksk_122[k]
                   + f_3 * pc_y[k] * lsk_230[k];

        t_290[k] = f_10 * lsi0_177[k]
                   - f_11 * lsi1_177[k]
                   + f_3 * pc_z[k] * lsk_230[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pc_x, pc_z, ksk_237, lsi0_178, lsi0_189, \
                         lsi1_178, lsi1_189, lsk_231, lsk_232, \
                         lsk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_19 * ksk_237[k]
                   + f_4 * lsi0_189[k]
                   - f_5 * lsi1_189[k]
                   + f_3 * pc_x[k] * lsk_237[k];

        t_292[k] = f_3 * pc_z[k] * lsk_231[k];

        t_293[k] = f_4 * lsi0_178[k]
                   - f_5 * lsi1_178[k]
                   + f_3 * pc_z[k] * lsk_232[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_y, pc_z, ksk_128, lsi0_179, lsi0_180, \
                         lsi0_182, lsi1_179, lsi1_180, lsi1_182, lsk_233, lsk_234, \
                         lsk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * lsi0_179[k]
                   - f_7 * lsi1_179[k]
                   + f_3 * pc_z[k] * lsk_233[k];

        t_295[k] = f_8 * lsi0_180[k]
                   - f_9 * lsi1_180[k]
                   + f_3 * pc_z[k] * lsk_234[k];

        t_296[k] = f_17 * ksk_128[k]
                   + f_3 * pc_y[k] * lsk_236[k];

        t_297[k] = f_12 * lsi0_182[k]
                   - f_13 * lsi1_182[k]
                   + f_3 * pc_z[k] * lsk_236[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pc_x, pc_z, ksk_244, ksk_246, \
                         ksk_247, ksk_248, lsk_237, lsk_244, lsk_246, lsk_247, \
                         lsk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_19 * ksk_244[k]
                   + f_3 * pc_x[k] * lsk_244[k];

        t_299[k] = f_3 * pc_z[k] * lsk_237[k];

        t_300[k] = f_19 * ksk_246[k]
                   + f_3 * pc_x[k] * lsk_246[k];

        t_301[k] = f_19 * ksk_247[k]
                   + f_3 * pc_x[k] * lsk_247[k];

        t_302[k] = f_19 * ksk_248[k]
                   + f_3 * pc_x[k] * lsk_248[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, pc_y, ksk_136, ksk_249, ksk_250, \
                         ksk_251, lsi0_189, lsi1_189, lsk_244, lsk_249, lsk_250, \
                         lsk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_19 * ksk_249[k]
                   + f_3 * pc_x[k] * lsk_249[k];

        t_304[k] = f_19 * ksk_250[k]
                   + f_3 * pc_x[k] * lsk_250[k];

        t_305[k] = f_19 * ksk_251[k]
                   + f_3 * pc_x[k] * lsk_251[k];

        t_306[k] = f_17 * ksk_136[k]
                   + f_1 * lsi0_189[k]
                   - f_2 * lsi1_189[k]
                   + f_3 * pc_y[k] * lsk_244[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_z, lsi0_189, lsi0_190, lsi0_191, \
                         lsi1_189, lsi1_190, lsi1_191, lsk_244, lsk_245, lsk_246, \
                         lsk_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_z[k] * lsk_244[k];

        t_308[k] = f_4 * lsi0_189[k]
                   - f_5 * lsi1_189[k]
                   + f_3 * pc_z[k] * lsk_245[k];

        t_309[k] = f_6 * lsi0_190[k]
                   - f_7 * lsi1_190[k]
                   + f_3 * pc_z[k] * lsk_246[k];

        t_310[k] = f_8 * lsi0_191[k]
                   - f_9 * lsi1_191[k]
                   + f_3 * pc_z[k] * lsk_247[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, ksk_143, lsi0_192, lsi0_193, \
                         lsi0_195, lsi1_192, lsi1_193, lsi1_195, lsk_248, lsk_249, \
                         lsk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_10 * lsi0_192[k]
                   - f_11 * lsi1_192[k]
                   + f_3 * pc_z[k] * lsk_248[k];

        t_312[k] = f_12 * lsi0_193[k]
                   - f_13 * lsi1_193[k]
                   + f_3 * pc_z[k] * lsk_249[k];

        t_313[k] = f_17 * ksk_143[k]
                   + f_3 * pc_y[k] * lsk_251[k];

        t_314[k] = f_1 * lsi0_195[k]
                   - f_2 * lsi1_195[k]
                   + f_3 * pc_z[k] * lsk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_z, pc_y, pc_z, ksl0_135, ksl0_138, \
                         ksk_108, ksk_144, ksl1_135, ksl1_138, \
                         lsk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * ksl0_135[k]
                   - f_14 * pc_z[k] * ksl1_135[k];

        t_316[k] = f_16 * ksk_144[k]
                   + f_3 * pc_y[k] * lsk_252[k];

        t_317[k] = f_15 * ksk_108[k]
                   + f_3 * pc_z[k] * lsk_252[k];

        t_318[k] = pa_z[k] * ksl0_138[k]
                   - f_14 * pc_z[k] * ksl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_x, pc_y, pc_z, ksl0_141, ksk_146, \
                         ksk_257, ksl1_141, lsi0_201, lsi1_201, lsk_254, \
                         lsk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * ksk_146[k]
                   + f_3 * pc_y[k] * lsk_254[k];

        t_320[k] = f_19 * ksk_257[k]
                   + f_12 * lsi0_201[k]
                   - f_13 * lsi1_201[k]
                   + f_3 * pc_x[k] * lsk_257[k];

        t_321[k] = pa_z[k] * ksl0_141[k]
                   - f_14 * pc_z[k] * ksl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, ksk_111, ksk_149, ksk_261, \
                         lsi0_205, lsi1_205, lsk_255, lsk_257, \
                         lsk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * ksk_111[k]
                   + f_3 * pc_z[k] * lsk_255[k];

        t_323[k] = f_16 * ksk_149[k]
                   + f_3 * pc_y[k] * lsk_257[k];

        t_324[k] = f_19 * ksk_261[k]
                   + f_10 * lsi0_205[k]
                   - f_11 * lsi1_205[k]
                   + f_3 * pc_x[k] * lsk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pa_z, pc_y, pc_z, ksl0_145, ksl0_147, \
                         ksk_114, ksk_115, ksk_153, ksl1_145, ksl1_147, lsk_258, \
                         lsk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_z[k] * ksl0_145[k]
                   - f_14 * pc_z[k] * ksl1_145[k];

        t_326[k] = f_15 * ksk_114[k]
                   + f_3 * pc_z[k] * lsk_258[k];

        t_327[k] = pa_z[k] * ksl0_147[k]
                   + f_16 * ksk_115[k]
                   - f_14 * pc_z[k] * ksl1_147[k];

        t_328[k] = f_16 * ksk_153[k]
                   + f_3 * pc_y[k] * lsk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_z, pc_x, pc_z, ksl0_150, ksk_118, ksk_266, \
                         ksl1_150, lsi0_210, lsi1_210, lsk_262, \
                         lsk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_19 * ksk_266[k]
                   + f_8 * lsi0_210[k]
                   - f_9 * lsi1_210[k]
                   + f_3 * pc_x[k] * lsk_266[k];

        t_330[k] = pa_z[k] * ksl0_150[k]
                   - f_14 * pc_z[k] * ksl1_150[k];

        t_331[k] = f_15 * ksk_118[k]
                   + f_3 * pc_z[k] * lsk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pa_z, pc_y, pc_z, ksl0_152, ksl0_153, ksk_119, \
                         ksk_120, ksk_158, ksl1_152, ksl1_153, \
                         lsk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pa_z[k] * ksl0_152[k]
                   + f_16 * ksk_119[k]
                   - f_14 * pc_z[k] * ksl1_152[k];

        t_333[k] = pa_z[k] * ksl0_153[k]
                   + f_17 * ksk_120[k]
                   - f_14 * pc_z[k] * ksl1_153[k];

        t_334[k] = f_16 * ksk_158[k]
                   + f_3 * pc_y[k] * lsk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pa_z, pc_x, pc_z, ksl0_156, ksk_123, ksk_272, \
                         ksl1_156, lsi0_216, lsi1_216, lsk_267, \
                         lsk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_19 * ksk_272[k]
                   + f_6 * lsi0_216[k]
                   - f_7 * lsi1_216[k]
                   + f_3 * pc_x[k] * lsk_272[k];

        t_336[k] = pa_z[k] * ksl0_156[k]
                   - f_14 * pc_z[k] * ksl1_156[k];

        t_337[k] = f_15 * ksk_123[k]
                   + f_3 * pc_z[k] * lsk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pa_z, pc_z, ksl0_158, ksl0_159, ksl0_160, \
                         ksk_124, ksk_125, ksk_126, ksl1_158, ksl1_159, \
                         ksl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pa_z[k] * ksl0_158[k]
                   + f_16 * ksk_124[k]
                   - f_14 * pc_z[k] * ksl1_158[k];

        t_339[k] = pa_z[k] * ksl0_159[k]
                   + f_17 * ksk_125[k]
                   - f_14 * pc_z[k] * ksl1_159[k];

        t_340[k] = pa_z[k] * ksl0_160[k]
                   + f_18 * ksk_126[k]
                   - f_14 * pc_z[k] * ksl1_160[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, ksk_164, ksk_279, ksk_280, \
                         ksk_281, lsi0_223, lsi1_223, lsk_272, lsk_279, lsk_280, \
                         lsk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * ksk_164[k]
                   + f_3 * pc_y[k] * lsk_272[k];

        t_342[k] = f_19 * ksk_279[k]
                   + f_4 * lsi0_223[k]
                   - f_5 * lsi1_223[k]
                   + f_3 * pc_x[k] * lsk_279[k];

        t_343[k] = f_19 * ksk_280[k]
                   + f_3 * pc_x[k] * lsk_280[k];

        t_344[k] = f_19 * ksk_281[k]
                   + f_3 * pc_x[k] * lsk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, ksk_282, ksk_283, ksk_284, \
                         ksk_285, ksk_286, lsk_282, lsk_283, lsk_284, lsk_285, \
                         lsk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_19 * ksk_282[k]
                   + f_3 * pc_x[k] * lsk_282[k];

        t_346[k] = f_19 * ksk_283[k]
                   + f_3 * pc_x[k] * lsk_283[k];

        t_347[k] = f_19 * ksk_284[k]
                   + f_3 * pc_x[k] * lsk_284[k];

        t_348[k] = f_19 * ksk_285[k]
                   + f_3 * pc_x[k] * lsk_285[k];

        t_349[k] = f_19 * ksk_286[k]
                   + f_3 * pc_x[k] * lsk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pa_z, pc_x, pc_z, ksl0_171, ksk_136, ksk_287, \
                         ksl1_171, lsk_280, lsk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_19 * ksk_287[k]
                   + f_3 * pc_x[k] * lsk_287[k];

        t_351[k] = pa_z[k] * ksl0_171[k]
                   - f_14 * pc_z[k] * ksl1_171[k];

        t_352[k] = f_15 * ksk_136[k]
                   + f_3 * pc_z[k] * lsk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, ksk_174, ksk_175, ksk_176, lsi0_219, \
                         lsi0_220, lsi0_221, lsi1_219, lsi1_220, lsi1_221, lsk_282, lsk_283, \
                         lsk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * ksk_174[k]
                   + f_12 * lsi0_219[k]
                   - f_13 * lsi1_219[k]
                   + f_3 * pc_y[k] * lsk_282[k];

        t_354[k] = f_16 * ksk_175[k]
                   + f_10 * lsi0_220[k]
                   - f_11 * lsi1_220[k]
                   + f_3 * pc_y[k] * lsk_283[k];

        t_355[k] = f_16 * ksk_176[k]
                   + f_8 * lsi0_221[k]
                   - f_9 * lsi1_221[k]
                   + f_3 * pc_y[k] * lsk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, ksk_177, ksk_178, ksk_179, lsi0_222, \
                         lsi0_223, lsi1_222, lsi1_223, lsk_285, lsk_286, \
                         lsk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * ksk_177[k]
                   + f_6 * lsi0_222[k]
                   - f_7 * lsi1_222[k]
                   + f_3 * pc_y[k] * lsk_285[k];

        t_357[k] = f_16 * ksk_178[k]
                   + f_4 * lsi0_223[k]
                   - f_5 * lsi1_223[k]
                   + f_3 * pc_y[k] * lsk_286[k];

        t_358[k] = f_16 * ksk_179[k]
                   + f_3 * pc_y[k] * lsk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_y, pc_y, pc_z, ksl0_225, ksk_143, \
                         ksk_144, ksk_180, ksl1_225, lsi0_223, lsi1_223, lsk_287, \
                         lsk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * ksk_143[k]
                   + f_1 * lsi0_223[k]
                   - f_2 * lsi1_223[k]
                   + f_3 * pc_z[k] * lsk_287[k];

        t_360[k] = pa_y[k] * ksl0_225[k]
                   - f_14 * pc_y[k] * ksl1_225[k];

        t_361[k] = f_15 * ksk_180[k]
                   + f_3 * pc_y[k] * lsk_288[k];

        t_362[k] = f_16 * ksk_144[k]
                   + f_3 * pc_z[k] * lsk_288[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_228 = buffer.data(ksl0 + 228);
    const auto *ksl0_230 = buffer.data(ksl0 + 230);
    const auto *ksl0_231 = buffer.data(ksl0 + 231);
    const auto *ksl0_234 = buffer.data(ksl0 + 234);
    const auto *ksl0_235 = buffer.data(ksl0 + 235);
    const auto *ksl0_237 = buffer.data(ksl0 + 237);
    const auto *ksl0_239 = buffer.data(ksl0 + 239);
    const auto *ksl0_240 = buffer.data(ksl0 + 240);
    const auto *ksl0_242 = buffer.data(ksl0 + 242);
    const auto *ksl0_243 = buffer.data(ksl0 + 243);
    const auto *ksl0_245 = buffer.data(ksl0 + 245);
    const auto *ksl0_246 = buffer.data(ksl0 + 246);
    const auto *ksl0_248 = buffer.data(ksl0 + 248);
    const auto *ksl0_249 = buffer.data(ksl0 + 249);
    const auto *ksl0_250 = buffer.data(ksl0 + 250);
    const auto *ksl0_252 = buffer.data(ksl0 + 252);
    const auto *ksl0_269 = buffer.data(ksl0 + 269);

    const auto *ksk_147 = buffer.data(ksk + 147);
    const auto *ksk_150 = buffer.data(ksk + 150);
    const auto *ksk_154 = buffer.data(ksk + 154);
    const auto *ksk_159 = buffer.data(ksk + 159);
    const auto *ksk_172 = buffer.data(ksk + 172);
    const auto *ksk_180 = buffer.data(ksk + 180);
    const auto *ksk_181 = buffer.data(ksk + 181);
    const auto *ksk_182 = buffer.data(ksk + 182);
    const auto *ksk_183 = buffer.data(ksk + 183);
    const auto *ksk_185 = buffer.data(ksk + 185);
    const auto *ksk_186 = buffer.data(ksk + 186);
    const auto *ksk_188 = buffer.data(ksk + 188);
    const auto *ksk_189 = buffer.data(ksk + 189);
    const auto *ksk_190 = buffer.data(ksk + 190);
    const auto *ksk_192 = buffer.data(ksk + 192);
    const auto *ksk_193 = buffer.data(ksk + 193);
    const auto *ksk_194 = buffer.data(ksk + 194);
    const auto *ksk_195 = buffer.data(ksk + 195);
    const auto *ksk_197 = buffer.data(ksk + 197);
    const auto *ksk_198 = buffer.data(ksk + 198);
    const auto *ksk_199 = buffer.data(ksk + 199);
    const auto *ksk_200 = buffer.data(ksk + 200);
    const auto *ksk_208 = buffer.data(ksk + 208);
    const auto *ksk_210 = buffer.data(ksk + 210);
    const auto *ksk_211 = buffer.data(ksk + 211);
    const auto *ksk_212 = buffer.data(ksk + 212);
    const auto *ksk_213 = buffer.data(ksk + 213);
    const auto *ksk_214 = buffer.data(ksk + 214);
    const auto *ksk_215 = buffer.data(ksk + 215);
    const auto *ksk_216 = buffer.data(ksk + 216);
    const auto *ksk_221 = buffer.data(ksk + 221);
    const auto *ksk_225 = buffer.data(ksk + 225);
    const auto *ksk_230 = buffer.data(ksk + 230);
    const auto *ksk_236 = buffer.data(ksk + 236);
    const auto *ksk_316 = buffer.data(ksk + 316);
    const auto *ksk_317 = buffer.data(ksk + 317);
    const auto *ksk_318 = buffer.data(ksk + 318);
    const auto *ksk_319 = buffer.data(ksk + 319);
    const auto *ksk_320 = buffer.data(ksk + 320);
    const auto *ksk_321 = buffer.data(ksk + 321);
    const auto *ksk_322 = buffer.data(ksk + 322);
    const auto *ksk_323 = buffer.data(ksk + 323);
    const auto *ksk_324 = buffer.data(ksk + 324);
    const auto *ksk_329 = buffer.data(ksk + 329);
    const auto *ksk_333 = buffer.data(ksk + 333);
    const auto *ksk_338 = buffer.data(ksk + 338);
    const auto *ksk_344 = buffer.data(ksk + 344);
    const auto *ksk_351 = buffer.data(ksk + 351);
    const auto *ksk_352 = buffer.data(ksk + 352);
    const auto *ksk_353 = buffer.data(ksk + 353);
    const auto *ksk_354 = buffer.data(ksk + 354);
    const auto *ksk_355 = buffer.data(ksk + 355);
    const auto *ksk_356 = buffer.data(ksk + 356);
    const auto *ksk_357 = buffer.data(ksk + 357);
    const auto *ksk_359 = buffer.data(ksk + 359);
    const auto *ksk_360 = buffer.data(ksk + 360);
    const auto *ksk_363 = buffer.data(ksk + 363);
    const auto *ksk_366 = buffer.data(ksk + 366);
    const auto *ksk_370 = buffer.data(ksk + 370);
    const auto *ksk_375 = buffer.data(ksk + 375);
    const auto *ksk_381 = buffer.data(ksk + 381);

    const auto *ksl1_228 = buffer.data(ksl1 + 228);
    const auto *ksl1_230 = buffer.data(ksl1 + 230);
    const auto *ksl1_231 = buffer.data(ksl1 + 231);
    const auto *ksl1_234 = buffer.data(ksl1 + 234);
    const auto *ksl1_235 = buffer.data(ksl1 + 235);
    const auto *ksl1_237 = buffer.data(ksl1 + 237);
    const auto *ksl1_239 = buffer.data(ksl1 + 239);
    const auto *ksl1_240 = buffer.data(ksl1 + 240);
    const auto *ksl1_242 = buffer.data(ksl1 + 242);
    const auto *ksl1_243 = buffer.data(ksl1 + 243);
    const auto *ksl1_245 = buffer.data(ksl1 + 245);
    const auto *ksl1_246 = buffer.data(ksl1 + 246);
    const auto *ksl1_248 = buffer.data(ksl1 + 248);
    const auto *ksl1_249 = buffer.data(ksl1 + 249);
    const auto *ksl1_250 = buffer.data(ksl1 + 250);
    const auto *ksl1_252 = buffer.data(ksl1 + 252);
    const auto *ksl1_269 = buffer.data(ksl1 + 269);

    const auto *lsi0_245 = buffer.data(lsi0 + 245);
    const auto *lsi0_247 = buffer.data(lsi0 + 247);
    const auto *lsi0_248 = buffer.data(lsi0 + 248);
    const auto *lsi0_249 = buffer.data(lsi0 + 249);
    const auto *lsi0_250 = buffer.data(lsi0 + 250);
    const auto *lsi0_251 = buffer.data(lsi0 + 251);
    const auto *lsi0_252 = buffer.data(lsi0 + 252);
    const auto *lsi0_253 = buffer.data(lsi0 + 253);
    const auto *lsi0_254 = buffer.data(lsi0 + 254);
    const auto *lsi0_255 = buffer.data(lsi0 + 255);
    const auto *lsi0_256 = buffer.data(lsi0 + 256);
    const auto *lsi0_257 = buffer.data(lsi0 + 257);
    const auto *lsi0_258 = buffer.data(lsi0 + 258);
    const auto *lsi0_259 = buffer.data(lsi0 + 259);
    const auto *lsi0_260 = buffer.data(lsi0 + 260);
    const auto *lsi0_261 = buffer.data(lsi0 + 261);
    const auto *lsi0_262 = buffer.data(lsi0 + 262);
    const auto *lsi0_263 = buffer.data(lsi0 + 263);
    const auto *lsi0_264 = buffer.data(lsi0 + 264);
    const auto *lsi0_265 = buffer.data(lsi0 + 265);
    const auto *lsi0_266 = buffer.data(lsi0 + 266);
    const auto *lsi0_272 = buffer.data(lsi0 + 272);
    const auto *lsi0_273 = buffer.data(lsi0 + 273);
    const auto *lsi0_274 = buffer.data(lsi0 + 274);
    const auto *lsi0_275 = buffer.data(lsi0 + 275);
    const auto *lsi0_276 = buffer.data(lsi0 + 276);
    const auto *lsi0_277 = buffer.data(lsi0 + 277);
    const auto *lsi0_278 = buffer.data(lsi0 + 278);
    const auto *lsi0_279 = buffer.data(lsi0 + 279);
    const auto *lsi0_280 = buffer.data(lsi0 + 280);
    const auto *lsi0_282 = buffer.data(lsi0 + 282);
    const auto *lsi0_283 = buffer.data(lsi0 + 283);
    const auto *lsi0_285 = buffer.data(lsi0 + 285);
    const auto *lsi0_286 = buffer.data(lsi0 + 286);
    const auto *lsi0_287 = buffer.data(lsi0 + 287);
    const auto *lsi0_289 = buffer.data(lsi0 + 289);
    const auto *lsi0_290 = buffer.data(lsi0 + 290);
    const auto *lsi0_291 = buffer.data(lsi0 + 291);
    const auto *lsi0_292 = buffer.data(lsi0 + 292);
    const auto *lsi0_294 = buffer.data(lsi0 + 294);
    const auto *lsi0_295 = buffer.data(lsi0 + 295);
    const auto *lsi0_301 = buffer.data(lsi0 + 301);

    const auto *lsi1_245 = buffer.data(lsi1 + 245);
    const auto *lsi1_247 = buffer.data(lsi1 + 247);
    const auto *lsi1_248 = buffer.data(lsi1 + 248);
    const auto *lsi1_249 = buffer.data(lsi1 + 249);
    const auto *lsi1_250 = buffer.data(lsi1 + 250);
    const auto *lsi1_251 = buffer.data(lsi1 + 251);
    const auto *lsi1_252 = buffer.data(lsi1 + 252);
    const auto *lsi1_253 = buffer.data(lsi1 + 253);
    const auto *lsi1_254 = buffer.data(lsi1 + 254);
    const auto *lsi1_255 = buffer.data(lsi1 + 255);
    const auto *lsi1_256 = buffer.data(lsi1 + 256);
    const auto *lsi1_257 = buffer.data(lsi1 + 257);
    const auto *lsi1_258 = buffer.data(lsi1 + 258);
    const auto *lsi1_259 = buffer.data(lsi1 + 259);
    const auto *lsi1_260 = buffer.data(lsi1 + 260);
    const auto *lsi1_261 = buffer.data(lsi1 + 261);
    const auto *lsi1_262 = buffer.data(lsi1 + 262);
    const auto *lsi1_263 = buffer.data(lsi1 + 263);
    const auto *lsi1_264 = buffer.data(lsi1 + 264);
    const auto *lsi1_265 = buffer.data(lsi1 + 265);
    const auto *lsi1_266 = buffer.data(lsi1 + 266);
    const auto *lsi1_272 = buffer.data(lsi1 + 272);
    const auto *lsi1_273 = buffer.data(lsi1 + 273);
    const auto *lsi1_274 = buffer.data(lsi1 + 274);
    const auto *lsi1_275 = buffer.data(lsi1 + 275);
    const auto *lsi1_276 = buffer.data(lsi1 + 276);
    const auto *lsi1_277 = buffer.data(lsi1 + 277);
    const auto *lsi1_278 = buffer.data(lsi1 + 278);
    const auto *lsi1_279 = buffer.data(lsi1 + 279);
    const auto *lsi1_280 = buffer.data(lsi1 + 280);
    const auto *lsi1_282 = buffer.data(lsi1 + 282);
    const auto *lsi1_283 = buffer.data(lsi1 + 283);
    const auto *lsi1_285 = buffer.data(lsi1 + 285);
    const auto *lsi1_286 = buffer.data(lsi1 + 286);
    const auto *lsi1_287 = buffer.data(lsi1 + 287);
    const auto *lsi1_289 = buffer.data(lsi1 + 289);
    const auto *lsi1_290 = buffer.data(lsi1 + 290);
    const auto *lsi1_291 = buffer.data(lsi1 + 291);
    const auto *lsi1_292 = buffer.data(lsi1 + 292);
    const auto *lsi1_294 = buffer.data(lsi1 + 294);
    const auto *lsi1_295 = buffer.data(lsi1 + 295);
    const auto *lsi1_301 = buffer.data(lsi1 + 301);

    const auto *lsk_290 = buffer.data(lsk + 290);
    const auto *lsk_291 = buffer.data(lsk + 291);
    const auto *lsk_293 = buffer.data(lsk + 293);
    const auto *lsk_294 = buffer.data(lsk + 294);
    const auto *lsk_297 = buffer.data(lsk + 297);
    const auto *lsk_298 = buffer.data(lsk + 298);
    const auto *lsk_302 = buffer.data(lsk + 302);
    const auto *lsk_303 = buffer.data(lsk + 303);
    const auto *lsk_308 = buffer.data(lsk + 308);
    const auto *lsk_316 = buffer.data(lsk + 316);
    const auto *lsk_317 = buffer.data(lsk + 317);
    const auto *lsk_318 = buffer.data(lsk + 318);
    const auto *lsk_319 = buffer.data(lsk + 319);
    const auto *lsk_320 = buffer.data(lsk + 320);
    const auto *lsk_321 = buffer.data(lsk + 321);
    const auto *lsk_322 = buffer.data(lsk + 322);
    const auto *lsk_323 = buffer.data(lsk + 323);
    const auto *lsk_324 = buffer.data(lsk + 324);
    const auto *lsk_325 = buffer.data(lsk + 325);
    const auto *lsk_326 = buffer.data(lsk + 326);
    const auto *lsk_327 = buffer.data(lsk + 327);
    const auto *lsk_328 = buffer.data(lsk + 328);
    const auto *lsk_329 = buffer.data(lsk + 329);
    const auto *lsk_330 = buffer.data(lsk + 330);
    const auto *lsk_331 = buffer.data(lsk + 331);
    const auto *lsk_332 = buffer.data(lsk + 332);
    const auto *lsk_333 = buffer.data(lsk + 333);
    const auto *lsk_334 = buffer.data(lsk + 334);
    const auto *lsk_335 = buffer.data(lsk + 335);
    const auto *lsk_336 = buffer.data(lsk + 336);
    const auto *lsk_337 = buffer.data(lsk + 337);
    const auto *lsk_338 = buffer.data(lsk + 338);
    const auto *lsk_339 = buffer.data(lsk + 339);
    const auto *lsk_340 = buffer.data(lsk + 340);
    const auto *lsk_341 = buffer.data(lsk + 341);
    const auto *lsk_342 = buffer.data(lsk + 342);
    const auto *lsk_343 = buffer.data(lsk + 343);
    const auto *lsk_344 = buffer.data(lsk + 344);
    const auto *lsk_351 = buffer.data(lsk + 351);
    const auto *lsk_352 = buffer.data(lsk + 352);
    const auto *lsk_353 = buffer.data(lsk + 353);
    const auto *lsk_354 = buffer.data(lsk + 354);
    const auto *lsk_355 = buffer.data(lsk + 355);
    const auto *lsk_356 = buffer.data(lsk + 356);
    const auto *lsk_357 = buffer.data(lsk + 357);
    const auto *lsk_358 = buffer.data(lsk + 358);
    const auto *lsk_359 = buffer.data(lsk + 359);
    const auto *lsk_360 = buffer.data(lsk + 360);
    const auto *lsk_361 = buffer.data(lsk + 361);
    const auto *lsk_362 = buffer.data(lsk + 362);
    const auto *lsk_363 = buffer.data(lsk + 363);
    const auto *lsk_365 = buffer.data(lsk + 365);
    const auto *lsk_366 = buffer.data(lsk + 366);
    const auto *lsk_367 = buffer.data(lsk + 367);
    const auto *lsk_369 = buffer.data(lsk + 369);
    const auto *lsk_370 = buffer.data(lsk + 370);
    const auto *lsk_371 = buffer.data(lsk + 371);
    const auto *lsk_372 = buffer.data(lsk + 372);
    const auto *lsk_374 = buffer.data(lsk + 374);
    const auto *lsk_375 = buffer.data(lsk + 375);
    const auto *lsk_376 = buffer.data(lsk + 376);
    const auto *lsk_377 = buffer.data(lsk + 377);
    const auto *lsk_378 = buffer.data(lsk + 378);
    const auto *lsk_380 = buffer.data(lsk + 380);
    const auto *lsk_381 = buffer.data(lsk + 381);

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_y, pc_y, ksl0_228, ksl0_230, ksl0_231, \
                         ksk_181, ksk_182, ksk_183, ksl1_228, ksl1_230, ksl1_231, \
                         lsk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_y[k] * ksl0_228[k]
                   + f_16 * ksk_181[k]
                   - f_14 * pc_y[k] * ksl1_228[k];

        t_364[k] = f_15 * ksk_182[k]
                   + f_3 * pc_y[k] * lsk_290[k];

        t_365[k] = pa_y[k] * ksl0_230[k]
                   - f_14 * pc_y[k] * ksl1_230[k];

        t_366[k] = pa_y[k] * ksl0_231[k]
                   + f_17 * ksk_183[k]
                   - f_14 * pc_y[k] * ksl1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_y, pc_y, pc_z, ksl0_234, ksl0_235, \
                         ksk_147, ksk_185, ksk_186, ksl1_234, ksl1_235, lsk_291, \
                         lsk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * ksk_147[k]
                   + f_3 * pc_z[k] * lsk_291[k];

        t_368[k] = f_15 * ksk_185[k]
                   + f_3 * pc_y[k] * lsk_293[k];

        t_369[k] = pa_y[k] * ksl0_234[k]
                   - f_14 * pc_y[k] * ksl1_234[k];

        t_370[k] = pa_y[k] * ksl0_235[k]
                   + f_18 * ksk_186[k]
                   - f_14 * pc_y[k] * ksl1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pc_y, pc_z, ksl0_237, ksl0_239, \
                         ksk_150, ksk_188, ksk_189, ksl1_237, ksl1_239, lsk_294, \
                         lsk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * ksk_150[k]
                   + f_3 * pc_z[k] * lsk_294[k];

        t_372[k] = pa_y[k] * ksl0_237[k]
                   + f_16 * ksk_188[k]
                   - f_14 * pc_y[k] * ksl1_237[k];

        t_373[k] = f_15 * ksk_189[k]
                   + f_3 * pc_y[k] * lsk_297[k];

        t_374[k] = pa_y[k] * ksl0_239[k]
                   - f_14 * pc_y[k] * ksl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_y, pc_y, pc_z, ksl0_240, ksl0_242, ksk_154, \
                         ksk_190, ksk_192, ksl1_240, ksl1_242, \
                         lsk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pa_y[k] * ksl0_240[k]
                   + f_19 * ksk_190[k]
                   - f_14 * pc_y[k] * ksl1_240[k];

        t_376[k] = f_16 * ksk_154[k]
                   + f_3 * pc_z[k] * lsk_298[k];

        t_377[k] = pa_y[k] * ksl0_242[k]
                   + f_17 * ksk_192[k]
                   - f_14 * pc_y[k] * ksl1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_y, pc_y, ksl0_243, ksl0_245, ksl0_246, \
                         ksk_193, ksk_194, ksk_195, ksl1_243, ksl1_245, ksl1_246, \
                         lsk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * ksl0_243[k]
                   + f_16 * ksk_193[k]
                   - f_14 * pc_y[k] * ksl1_243[k];

        t_379[k] = f_15 * ksk_194[k]
                   + f_3 * pc_y[k] * lsk_302[k];

        t_380[k] = pa_y[k] * ksl0_245[k]
                   - f_14 * pc_y[k] * ksl1_245[k];

        t_381[k] = pa_y[k] * ksl0_246[k]
                   + f_20 * ksk_195[k]
                   - f_14 * pc_y[k] * ksl1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pa_y, pc_y, pc_z, ksl0_248, ksl0_249, ksk_159, \
                         ksk_197, ksk_198, ksl1_248, ksl1_249, \
                         lsk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * ksk_159[k]
                   + f_3 * pc_z[k] * lsk_303[k];

        t_383[k] = pa_y[k] * ksl0_248[k]
                   + f_18 * ksk_197[k]
                   - f_14 * pc_y[k] * ksl1_248[k];

        t_384[k] = pa_y[k] * ksl0_249[k]
                   + f_17 * ksk_198[k]
                   - f_14 * pc_y[k] * ksl1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_y, pc_x, pc_y, ksl0_250, ksl0_252, \
                         ksk_199, ksk_200, ksk_316, ksl1_250, ksl1_252, lsk_308, \
                         lsk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pa_y[k] * ksl0_250[k]
                   + f_16 * ksk_199[k]
                   - f_14 * pc_y[k] * ksl1_250[k];

        t_386[k] = f_15 * ksk_200[k]
                   + f_3 * pc_y[k] * lsk_308[k];

        t_387[k] = pa_y[k] * ksl0_252[k]
                   - f_14 * pc_y[k] * ksl1_252[k];

        t_388[k] = f_19 * ksk_316[k]
                   + f_3 * pc_x[k] * lsk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, ksk_317, ksk_318, ksk_319, \
                         ksk_320, ksk_321, lsk_317, lsk_318, lsk_319, lsk_320, \
                         lsk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_19 * ksk_317[k]
                   + f_3 * pc_x[k] * lsk_317[k];

        t_390[k] = f_19 * ksk_318[k]
                   + f_3 * pc_x[k] * lsk_318[k];

        t_391[k] = f_19 * ksk_319[k]
                   + f_3 * pc_x[k] * lsk_319[k];

        t_392[k] = f_19 * ksk_320[k]
                   + f_3 * pc_x[k] * lsk_320[k];

        t_393[k] = f_19 * ksk_321[k]
                   + f_3 * pc_x[k] * lsk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, ksk_172, ksk_208, \
                         ksk_322, ksk_323, lsi0_245, lsi1_245, lsk_316, lsk_322, \
                         lsk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_19 * ksk_322[k]
                   + f_3 * pc_x[k] * lsk_322[k];

        t_395[k] = f_19 * ksk_323[k]
                   + f_3 * pc_x[k] * lsk_323[k];

        t_396[k] = f_15 * ksk_208[k]
                   + f_1 * lsi0_245[k]
                   - f_2 * lsi1_245[k]
                   + f_3 * pc_y[k] * lsk_316[k];

        t_397[k] = f_16 * ksk_172[k]
                   + f_3 * pc_z[k] * lsk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, ksk_210, ksk_211, ksk_212, lsi0_247, \
                         lsi0_248, lsi0_249, lsi1_247, lsi1_248, lsi1_249, lsk_318, lsk_319, \
                         lsk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * ksk_210[k]
                   + f_12 * lsi0_247[k]
                   - f_13 * lsi1_247[k]
                   + f_3 * pc_y[k] * lsk_318[k];

        t_399[k] = f_15 * ksk_211[k]
                   + f_10 * lsi0_248[k]
                   - f_11 * lsi1_248[k]
                   + f_3 * pc_y[k] * lsk_319[k];

        t_400[k] = f_15 * ksk_212[k]
                   + f_8 * lsi0_249[k]
                   - f_9 * lsi1_249[k]
                   + f_3 * pc_y[k] * lsk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, ksk_213, ksk_214, ksk_215, lsi0_250, \
                         lsi0_251, lsi1_250, lsi1_251, lsk_321, lsk_322, \
                         lsk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * ksk_213[k]
                   + f_6 * lsi0_250[k]
                   - f_7 * lsi1_250[k]
                   + f_3 * pc_y[k] * lsk_321[k];

        t_402[k] = f_15 * ksk_214[k]
                   + f_4 * lsi0_251[k]
                   - f_5 * lsi1_251[k]
                   + f_3 * pc_y[k] * lsk_322[k];

        t_403[k] = f_15 * ksk_215[k]
                   + f_3 * pc_y[k] * lsk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_y, pc_x, pc_y, pc_z, ksl0_269, \
                         ksk_180, ksk_324, ksl1_269, lsi0_252, lsi1_252, \
                         lsk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * ksl0_269[k]
                   - f_14 * pc_y[k] * ksl1_269[k];

        t_405[k] = f_19 * ksk_324[k]
                   + f_1 * lsi0_252[k]
                   - f_2 * lsi1_252[k]
                   + f_3 * pc_x[k] * lsk_324[k];

        t_406[k] = f_3 * pc_y[k] * lsk_324[k];

        t_407[k] = f_17 * ksk_180[k]
                   + f_3 * pc_z[k] * lsk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, ksk_329, lsi0_252, lsi0_257, \
                         lsi1_252, lsi1_257, lsk_325, lsk_326, \
                         lsk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_4 * lsi0_252[k]
                   - f_5 * lsi1_252[k]
                   + f_3 * pc_y[k] * lsk_325[k];

        t_409[k] = f_3 * pc_y[k] * lsk_326[k];

        t_410[k] = f_19 * ksk_329[k]
                   + f_12 * lsi0_257[k]
                   - f_13 * lsi1_257[k]
                   + f_3 * pc_x[k] * lsk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_y, lsi0_253, lsi0_254, lsi1_253, lsi1_254, \
                         lsk_327, lsk_328, lsk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_6 * lsi0_253[k]
                   - f_7 * lsi1_253[k]
                   + f_3 * pc_y[k] * lsk_327[k];

        t_412[k] = f_4 * lsi0_254[k]
                   - f_5 * lsi1_254[k]
                   + f_3 * pc_y[k] * lsk_328[k];

        t_413[k] = f_3 * pc_y[k] * lsk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, ksk_333, lsi0_255, lsi0_256, \
                         lsi0_261, lsi1_255, lsi1_256, lsi1_261, lsk_330, lsk_331, \
                         lsk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_19 * ksk_333[k]
                   + f_10 * lsi0_261[k]
                   - f_11 * lsi1_261[k]
                   + f_3 * pc_x[k] * lsk_333[k];

        t_415[k] = f_8 * lsi0_255[k]
                   - f_9 * lsi1_255[k]
                   + f_3 * pc_y[k] * lsk_330[k];

        t_416[k] = f_6 * lsi0_256[k]
                   - f_7 * lsi1_256[k]
                   + f_3 * pc_y[k] * lsk_331[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, ksk_338, lsi0_257, lsi0_266, \
                         lsi1_257, lsi1_266, lsk_332, lsk_333, \
                         lsk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * lsi0_257[k]
                   - f_5 * lsi1_257[k]
                   + f_3 * pc_y[k] * lsk_332[k];

        t_418[k] = f_3 * pc_y[k] * lsk_333[k];

        t_419[k] = f_19 * ksk_338[k]
                   + f_8 * lsi0_266[k]
                   - f_9 * lsi1_266[k]
                   + f_3 * pc_x[k] * lsk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_y, lsi0_258, lsi0_259, lsi0_260, lsi1_258, \
                         lsi1_259, lsi1_260, lsk_334, lsk_335, \
                         lsk_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_10 * lsi0_258[k]
                   - f_11 * lsi1_258[k]
                   + f_3 * pc_y[k] * lsk_334[k];

        t_421[k] = f_8 * lsi0_259[k]
                   - f_9 * lsi1_259[k]
                   + f_3 * pc_y[k] * lsk_335[k];

        t_422[k] = f_6 * lsi0_260[k]
                   - f_7 * lsi1_260[k]
                   + f_3 * pc_y[k] * lsk_336[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, ksk_344, lsi0_261, lsi0_272, \
                         lsi1_261, lsi1_272, lsk_337, lsk_338, \
                         lsk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * lsi0_261[k]
                   - f_5 * lsi1_261[k]
                   + f_3 * pc_y[k] * lsk_337[k];

        t_424[k] = f_3 * pc_y[k] * lsk_338[k];

        t_425[k] = f_19 * ksk_344[k]
                   + f_6 * lsi0_272[k]
                   - f_7 * lsi1_272[k]
                   + f_3 * pc_x[k] * lsk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_y, lsi0_262, lsi0_263, lsi0_264, lsi1_262, \
                         lsi1_263, lsi1_264, lsk_339, lsk_340, \
                         lsk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_12 * lsi0_262[k]
                   - f_13 * lsi1_262[k]
                   + f_3 * pc_y[k] * lsk_339[k];

        t_427[k] = f_10 * lsi0_263[k]
                   - f_11 * lsi1_263[k]
                   + f_3 * pc_y[k] * lsk_340[k];

        t_428[k] = f_8 * lsi0_264[k]
                   - f_9 * lsi1_264[k]
                   + f_3 * pc_y[k] * lsk_341[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_y, lsi0_265, lsi0_266, lsi1_265, lsi1_266, \
                         lsk_342, lsk_343, lsk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_6 * lsi0_265[k]
                   - f_7 * lsi1_265[k]
                   + f_3 * pc_y[k] * lsk_342[k];

        t_430[k] = f_4 * lsi0_266[k]
                   - f_5 * lsi1_266[k]
                   + f_3 * pc_y[k] * lsk_343[k];

        t_431[k] = f_3 * pc_y[k] * lsk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, ksk_351, ksk_352, ksk_353, ksk_354, \
                         lsi0_279, lsi1_279, lsk_351, lsk_352, lsk_353, \
                         lsk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_19 * ksk_351[k]
                   + f_4 * lsi0_279[k]
                   - f_5 * lsi1_279[k]
                   + f_3 * pc_x[k] * lsk_351[k];

        t_433[k] = f_19 * ksk_352[k]
                   + f_3 * pc_x[k] * lsk_352[k];

        t_434[k] = f_19 * ksk_353[k]
                   + f_3 * pc_x[k] * lsk_353[k];

        t_435[k] = f_19 * ksk_354[k]
                   + f_3 * pc_x[k] * lsk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_y, ksk_355, ksk_356, \
                         ksk_357, ksk_359, lsk_351, lsk_355, lsk_356, lsk_357, \
                         lsk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_19 * ksk_355[k]
                   + f_3 * pc_x[k] * lsk_355[k];

        t_437[k] = f_19 * ksk_356[k]
                   + f_3 * pc_x[k] * lsk_356[k];

        t_438[k] = f_19 * ksk_357[k]
                   + f_3 * pc_x[k] * lsk_357[k];

        t_439[k] = f_3 * pc_y[k] * lsk_351[k];

        t_440[k] = f_19 * ksk_359[k]
                   + f_3 * pc_x[k] * lsk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, lsi0_273, lsi0_274, lsi0_275, lsi1_273, \
                         lsi1_274, lsi1_275, lsk_352, lsk_353, \
                         lsk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * lsi0_273[k]
                   - f_2 * lsi1_273[k]
                   + f_3 * pc_y[k] * lsk_352[k];

        t_442[k] = f_22 * lsi0_274[k]
                   - f_23 * lsi1_274[k]
                   + f_3 * pc_y[k] * lsk_353[k];

        t_443[k] = f_12 * lsi0_275[k]
                   - f_13 * lsi1_275[k]
                   + f_3 * pc_y[k] * lsk_354[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_y, lsi0_276, lsi0_277, lsi0_278, lsi1_276, \
                         lsi1_277, lsi1_278, lsk_355, lsk_356, \
                         lsk_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_10 * lsi0_276[k]
                   - f_11 * lsi1_276[k]
                   + f_3 * pc_y[k] * lsk_355[k];

        t_445[k] = f_8 * lsi0_277[k]
                   - f_9 * lsi1_277[k]
                   + f_3 * pc_y[k] * lsk_356[k];

        t_446[k] = f_6 * lsi0_278[k]
                   - f_7 * lsi1_278[k]
                   + f_3 * pc_y[k] * lsk_357[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pc_x, pc_y, pc_z, ksk_215, ksk_360, \
                         lsi0_279, lsi0_280, lsi1_279, lsi1_280, lsk_358, lsk_359, \
                         lsk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_4 * lsi0_279[k]
                   - f_5 * lsi1_279[k]
                   + f_3 * pc_y[k] * lsk_358[k];

        t_448[k] = f_3 * pc_y[k] * lsk_359[k];

        t_449[k] = f_17 * ksk_215[k]
                   + f_1 * lsi0_279[k]
                   - f_2 * lsi1_279[k]
                   + f_3 * pc_z[k] * lsk_359[k];

        t_450[k] = f_18 * ksk_360[k]
                   + f_1 * lsi0_280[k]
                   - f_2 * lsi1_280[k]
                   + f_3 * pc_x[k] * lsk_360[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pc_x, pc_y, pc_z, ksk_216, ksk_363, \
                         lsi0_283, lsi1_283, lsk_360, lsk_361, \
                         lsk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_18 * ksk_216[k]
                   + f_3 * pc_y[k] * lsk_360[k];

        t_452[k] = f_3 * pc_z[k] * lsk_360[k];

        t_453[k] = f_18 * ksk_363[k]
                   + f_12 * lsi0_283[k]
                   - f_13 * lsi1_283[k]
                   + f_3 * pc_x[k] * lsk_363[k];

        t_454[k] = f_3 * pc_z[k] * lsk_361[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_z, ksk_366, lsi0_280, lsi0_286, \
                         lsi1_280, lsi1_286, lsk_362, lsk_363, \
                         lsk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_4 * lsi0_280[k]
                   - f_5 * lsi1_280[k]
                   + f_3 * pc_z[k] * lsk_362[k];

        t_456[k] = f_18 * ksk_366[k]
                   + f_10 * lsi0_286[k]
                   - f_11 * lsi1_286[k]
                   + f_3 * pc_x[k] * lsk_366[k];

        t_457[k] = f_3 * pc_z[k] * lsk_363[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pc_x, pc_y, pc_z, ksk_221, ksk_370, \
                         lsi0_282, lsi0_290, lsi1_282, lsi1_290, lsk_365, lsk_366, \
                         lsk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_18 * ksk_221[k]
                   + f_3 * pc_y[k] * lsk_365[k];

        t_459[k] = f_6 * lsi0_282[k]
                   - f_7 * lsi1_282[k]
                   + f_3 * pc_z[k] * lsk_365[k];

        t_460[k] = f_18 * ksk_370[k]
                   + f_8 * lsi0_290[k]
                   - f_9 * lsi1_290[k]
                   + f_3 * pc_x[k] * lsk_370[k];

        t_461[k] = f_3 * pc_z[k] * lsk_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, ksk_225, lsi0_283, lsi0_285, \
                         lsi1_283, lsi1_285, lsk_367, lsk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_4 * lsi0_283[k]
                   - f_5 * lsi1_283[k]
                   + f_3 * pc_z[k] * lsk_367[k];

        t_463[k] = f_18 * ksk_225[k]
                   + f_3 * pc_y[k] * lsk_369[k];

        t_464[k] = f_8 * lsi0_285[k]
                   - f_9 * lsi1_285[k]
                   + f_3 * pc_z[k] * lsk_369[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_z, ksk_375, lsi0_286, lsi0_295, \
                         lsi1_286, lsi1_295, lsk_370, lsk_371, \
                         lsk_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_18 * ksk_375[k]
                   + f_6 * lsi0_295[k]
                   - f_7 * lsi1_295[k]
                   + f_3 * pc_x[k] * lsk_375[k];

        t_466[k] = f_3 * pc_z[k] * lsk_370[k];

        t_467[k] = f_4 * lsi0_286[k]
                   - f_5 * lsi1_286[k]
                   + f_3 * pc_z[k] * lsk_371[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_y, pc_z, ksk_230, lsi0_287, lsi0_289, \
                         lsi1_287, lsi1_289, lsk_372, lsk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_6 * lsi0_287[k]
                   - f_7 * lsi1_287[k]
                   + f_3 * pc_z[k] * lsk_372[k];

        t_469[k] = f_18 * ksk_230[k]
                   + f_3 * pc_y[k] * lsk_374[k];

        t_470[k] = f_10 * lsi0_289[k]
                   - f_11 * lsi1_289[k]
                   + f_3 * pc_z[k] * lsk_374[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_x, pc_z, ksk_381, lsi0_290, lsi0_301, \
                         lsi1_290, lsi1_301, lsk_375, lsk_376, \
                         lsk_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_18 * ksk_381[k]
                   + f_4 * lsi0_301[k]
                   - f_5 * lsi1_301[k]
                   + f_3 * pc_x[k] * lsk_381[k];

        t_472[k] = f_3 * pc_z[k] * lsk_375[k];

        t_473[k] = f_4 * lsi0_290[k]
                   - f_5 * lsi1_290[k]
                   + f_3 * pc_z[k] * lsk_376[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_y, pc_z, ksk_236, lsi0_291, lsi0_292, \
                         lsi0_294, lsi1_291, lsi1_292, lsi1_294, lsk_377, lsk_378, \
                         lsk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_6 * lsi0_291[k]
                   - f_7 * lsi1_291[k]
                   + f_3 * pc_z[k] * lsk_377[k];

        t_475[k] = f_8 * lsi0_292[k]
                   - f_9 * lsi1_292[k]
                   + f_3 * pc_z[k] * lsk_378[k];

        t_476[k] = f_18 * ksk_236[k]
                   + f_3 * pc_y[k] * lsk_380[k];

        t_477[k] = f_12 * lsi0_294[k]
                   - f_13 * lsi1_294[k]
                   + f_3 * pc_z[k] * lsk_380[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_270 = buffer.data(ksl0 + 270);
    const auto *ksl0_273 = buffer.data(ksl0 + 273);
    const auto *ksl0_276 = buffer.data(ksl0 + 276);
    const auto *ksl0_280 = buffer.data(ksl0 + 280);
    const auto *ksl0_282 = buffer.data(ksl0 + 282);
    const auto *ksl0_285 = buffer.data(ksl0 + 285);
    const auto *ksl0_287 = buffer.data(ksl0 + 287);
    const auto *ksl0_288 = buffer.data(ksl0 + 288);
    const auto *ksl0_291 = buffer.data(ksl0 + 291);
    const auto *ksl0_293 = buffer.data(ksl0 + 293);
    const auto *ksl0_294 = buffer.data(ksl0 + 294);
    const auto *ksl0_295 = buffer.data(ksl0 + 295);
    const auto *ksl0_306 = buffer.data(ksl0 + 306);
    const auto *ksl0_405 = buffer.data(ksl0 + 405);

    const auto *ksk_216 = buffer.data(ksk + 216);
    const auto *ksk_219 = buffer.data(ksk + 219);
    const auto *ksk_222 = buffer.data(ksk + 222);
    const auto *ksk_223 = buffer.data(ksk + 223);
    const auto *ksk_226 = buffer.data(ksk + 226);
    const auto *ksk_227 = buffer.data(ksk + 227);
    const auto *ksk_228 = buffer.data(ksk + 228);
    const auto *ksk_231 = buffer.data(ksk + 231);
    const auto *ksk_232 = buffer.data(ksk + 232);
    const auto *ksk_233 = buffer.data(ksk + 233);
    const auto *ksk_234 = buffer.data(ksk + 234);
    const auto *ksk_244 = buffer.data(ksk + 244);
    const auto *ksk_251 = buffer.data(ksk + 251);
    const auto *ksk_252 = buffer.data(ksk + 252);
    const auto *ksk_254 = buffer.data(ksk + 254);
    const auto *ksk_255 = buffer.data(ksk + 255);
    const auto *ksk_257 = buffer.data(ksk + 257);
    const auto *ksk_258 = buffer.data(ksk + 258);
    const auto *ksk_261 = buffer.data(ksk + 261);
    const auto *ksk_262 = buffer.data(ksk + 262);
    const auto *ksk_266 = buffer.data(ksk + 266);
    const auto *ksk_267 = buffer.data(ksk + 267);
    const auto *ksk_272 = buffer.data(ksk + 272);
    const auto *ksk_280 = buffer.data(ksk + 280);
    const auto *ksk_282 = buffer.data(ksk + 282);
    const auto *ksk_283 = buffer.data(ksk + 283);
    const auto *ksk_284 = buffer.data(ksk + 284);
    const auto *ksk_285 = buffer.data(ksk + 285);
    const auto *ksk_286 = buffer.data(ksk + 286);
    const auto *ksk_287 = buffer.data(ksk + 287);
    const auto *ksk_288 = buffer.data(ksk + 288);
    const auto *ksk_290 = buffer.data(ksk + 290);
    const auto *ksk_293 = buffer.data(ksk + 293);
    const auto *ksk_297 = buffer.data(ksk + 297);
    const auto *ksk_302 = buffer.data(ksk + 302);
    const auto *ksk_308 = buffer.data(ksk + 308);
    const auto *ksk_316 = buffer.data(ksk + 316);
    const auto *ksk_318 = buffer.data(ksk + 318);
    const auto *ksk_319 = buffer.data(ksk + 319);
    const auto *ksk_320 = buffer.data(ksk + 320);
    const auto *ksk_321 = buffer.data(ksk + 321);
    const auto *ksk_322 = buffer.data(ksk + 322);
    const auto *ksk_323 = buffer.data(ksk + 323);
    const auto *ksk_324 = buffer.data(ksk + 324);
    const auto *ksk_388 = buffer.data(ksk + 388);
    const auto *ksk_390 = buffer.data(ksk + 390);
    const auto *ksk_391 = buffer.data(ksk + 391);
    const auto *ksk_392 = buffer.data(ksk + 392);
    const auto *ksk_393 = buffer.data(ksk + 393);
    const auto *ksk_394 = buffer.data(ksk + 394);
    const auto *ksk_395 = buffer.data(ksk + 395);
    const auto *ksk_401 = buffer.data(ksk + 401);
    const auto *ksk_405 = buffer.data(ksk + 405);
    const auto *ksk_410 = buffer.data(ksk + 410);
    const auto *ksk_416 = buffer.data(ksk + 416);
    const auto *ksk_423 = buffer.data(ksk + 423);
    const auto *ksk_424 = buffer.data(ksk + 424);
    const auto *ksk_425 = buffer.data(ksk + 425);
    const auto *ksk_426 = buffer.data(ksk + 426);
    const auto *ksk_427 = buffer.data(ksk + 427);
    const auto *ksk_428 = buffer.data(ksk + 428);
    const auto *ksk_429 = buffer.data(ksk + 429);
    const auto *ksk_430 = buffer.data(ksk + 430);
    const auto *ksk_431 = buffer.data(ksk + 431);
    const auto *ksk_432 = buffer.data(ksk + 432);
    const auto *ksk_435 = buffer.data(ksk + 435);
    const auto *ksk_437 = buffer.data(ksk + 437);
    const auto *ksk_438 = buffer.data(ksk + 438);
    const auto *ksk_441 = buffer.data(ksk + 441);
    const auto *ksk_442 = buffer.data(ksk + 442);
    const auto *ksk_444 = buffer.data(ksk + 444);
    const auto *ksk_446 = buffer.data(ksk + 446);
    const auto *ksk_447 = buffer.data(ksk + 447);
    const auto *ksk_449 = buffer.data(ksk + 449);
    const auto *ksk_450 = buffer.data(ksk + 450);
    const auto *ksk_452 = buffer.data(ksk + 452);
    const auto *ksk_453 = buffer.data(ksk + 453);
    const auto *ksk_455 = buffer.data(ksk + 455);
    const auto *ksk_456 = buffer.data(ksk + 456);
    const auto *ksk_457 = buffer.data(ksk + 457);
    const auto *ksk_459 = buffer.data(ksk + 459);
    const auto *ksk_460 = buffer.data(ksk + 460);
    const auto *ksk_461 = buffer.data(ksk + 461);
    const auto *ksk_462 = buffer.data(ksk + 462);
    const auto *ksk_463 = buffer.data(ksk + 463);
    const auto *ksk_464 = buffer.data(ksk + 464);
    const auto *ksk_465 = buffer.data(ksk + 465);
    const auto *ksk_466 = buffer.data(ksk + 466);
    const auto *ksk_467 = buffer.data(ksk + 467);

    const auto *ksl1_270 = buffer.data(ksl1 + 270);
    const auto *ksl1_273 = buffer.data(ksl1 + 273);
    const auto *ksl1_276 = buffer.data(ksl1 + 276);
    const auto *ksl1_280 = buffer.data(ksl1 + 280);
    const auto *ksl1_282 = buffer.data(ksl1 + 282);
    const auto *ksl1_285 = buffer.data(ksl1 + 285);
    const auto *ksl1_287 = buffer.data(ksl1 + 287);
    const auto *ksl1_288 = buffer.data(ksl1 + 288);
    const auto *ksl1_291 = buffer.data(ksl1 + 291);
    const auto *ksl1_293 = buffer.data(ksl1 + 293);
    const auto *ksl1_294 = buffer.data(ksl1 + 294);
    const auto *ksl1_295 = buffer.data(ksl1 + 295);
    const auto *ksl1_306 = buffer.data(ksl1 + 306);
    const auto *ksl1_405 = buffer.data(ksl1 + 405);

    const auto *lsi0_301 = buffer.data(lsi0 + 301);
    const auto *lsi0_302 = buffer.data(lsi0 + 302);
    const auto *lsi0_303 = buffer.data(lsi0 + 303);
    const auto *lsi0_304 = buffer.data(lsi0 + 304);
    const auto *lsi0_305 = buffer.data(lsi0 + 305);
    const auto *lsi0_307 = buffer.data(lsi0 + 307);
    const auto *lsi0_313 = buffer.data(lsi0 + 313);
    const auto *lsi0_317 = buffer.data(lsi0 + 317);
    const auto *lsi0_322 = buffer.data(lsi0 + 322);
    const auto *lsi0_328 = buffer.data(lsi0 + 328);
    const auto *lsi0_331 = buffer.data(lsi0 + 331);
    const auto *lsi0_332 = buffer.data(lsi0 + 332);
    const auto *lsi0_333 = buffer.data(lsi0 + 333);
    const auto *lsi0_334 = buffer.data(lsi0 + 334);
    const auto *lsi0_335 = buffer.data(lsi0 + 335);
    const auto *lsi0_336 = buffer.data(lsi0 + 336);
    const auto *lsi0_339 = buffer.data(lsi0 + 339);
    const auto *lsi0_341 = buffer.data(lsi0 + 341);
    const auto *lsi0_342 = buffer.data(lsi0 + 342);
    const auto *lsi0_345 = buffer.data(lsi0 + 345);
    const auto *lsi0_346 = buffer.data(lsi0 + 346);
    const auto *lsi0_348 = buffer.data(lsi0 + 348);
    const auto *lsi0_350 = buffer.data(lsi0 + 350);
    const auto *lsi0_351 = buffer.data(lsi0 + 351);
    const auto *lsi0_353 = buffer.data(lsi0 + 353);
    const auto *lsi0_354 = buffer.data(lsi0 + 354);
    const auto *lsi0_356 = buffer.data(lsi0 + 356);
    const auto *lsi0_357 = buffer.data(lsi0 + 357);
    const auto *lsi0_359 = buffer.data(lsi0 + 359);
    const auto *lsi0_360 = buffer.data(lsi0 + 360);
    const auto *lsi0_361 = buffer.data(lsi0 + 361);
    const auto *lsi0_362 = buffer.data(lsi0 + 362);
    const auto *lsi0_363 = buffer.data(lsi0 + 363);

    const auto *lsi1_301 = buffer.data(lsi1 + 301);
    const auto *lsi1_302 = buffer.data(lsi1 + 302);
    const auto *lsi1_303 = buffer.data(lsi1 + 303);
    const auto *lsi1_304 = buffer.data(lsi1 + 304);
    const auto *lsi1_305 = buffer.data(lsi1 + 305);
    const auto *lsi1_307 = buffer.data(lsi1 + 307);
    const auto *lsi1_313 = buffer.data(lsi1 + 313);
    const auto *lsi1_317 = buffer.data(lsi1 + 317);
    const auto *lsi1_322 = buffer.data(lsi1 + 322);
    const auto *lsi1_328 = buffer.data(lsi1 + 328);
    const auto *lsi1_331 = buffer.data(lsi1 + 331);
    const auto *lsi1_332 = buffer.data(lsi1 + 332);
    const auto *lsi1_333 = buffer.data(lsi1 + 333);
    const auto *lsi1_334 = buffer.data(lsi1 + 334);
    const auto *lsi1_335 = buffer.data(lsi1 + 335);
    const auto *lsi1_336 = buffer.data(lsi1 + 336);
    const auto *lsi1_339 = buffer.data(lsi1 + 339);
    const auto *lsi1_341 = buffer.data(lsi1 + 341);
    const auto *lsi1_342 = buffer.data(lsi1 + 342);
    const auto *lsi1_345 = buffer.data(lsi1 + 345);
    const auto *lsi1_346 = buffer.data(lsi1 + 346);
    const auto *lsi1_348 = buffer.data(lsi1 + 348);
    const auto *lsi1_350 = buffer.data(lsi1 + 350);
    const auto *lsi1_351 = buffer.data(lsi1 + 351);
    const auto *lsi1_353 = buffer.data(lsi1 + 353);
    const auto *lsi1_354 = buffer.data(lsi1 + 354);
    const auto *lsi1_356 = buffer.data(lsi1 + 356);
    const auto *lsi1_357 = buffer.data(lsi1 + 357);
    const auto *lsi1_359 = buffer.data(lsi1 + 359);
    const auto *lsi1_360 = buffer.data(lsi1 + 360);
    const auto *lsi1_361 = buffer.data(lsi1 + 361);
    const auto *lsi1_362 = buffer.data(lsi1 + 362);
    const auto *lsi1_363 = buffer.data(lsi1 + 363);

    const auto *lsk_381 = buffer.data(lsk + 381);
    const auto *lsk_388 = buffer.data(lsk + 388);
    const auto *lsk_389 = buffer.data(lsk + 389);
    const auto *lsk_390 = buffer.data(lsk + 390);
    const auto *lsk_391 = buffer.data(lsk + 391);
    const auto *lsk_392 = buffer.data(lsk + 392);
    const auto *lsk_393 = buffer.data(lsk + 393);
    const auto *lsk_394 = buffer.data(lsk + 394);
    const auto *lsk_395 = buffer.data(lsk + 395);
    const auto *lsk_396 = buffer.data(lsk + 396);
    const auto *lsk_398 = buffer.data(lsk + 398);
    const auto *lsk_399 = buffer.data(lsk + 399);
    const auto *lsk_401 = buffer.data(lsk + 401);
    const auto *lsk_402 = buffer.data(lsk + 402);
    const auto *lsk_405 = buffer.data(lsk + 405);
    const auto *lsk_406 = buffer.data(lsk + 406);
    const auto *lsk_410 = buffer.data(lsk + 410);
    const auto *lsk_411 = buffer.data(lsk + 411);
    const auto *lsk_416 = buffer.data(lsk + 416);
    const auto *lsk_423 = buffer.data(lsk + 423);
    const auto *lsk_424 = buffer.data(lsk + 424);
    const auto *lsk_425 = buffer.data(lsk + 425);
    const auto *lsk_426 = buffer.data(lsk + 426);
    const auto *lsk_427 = buffer.data(lsk + 427);
    const auto *lsk_428 = buffer.data(lsk + 428);
    const auto *lsk_429 = buffer.data(lsk + 429);
    const auto *lsk_430 = buffer.data(lsk + 430);
    const auto *lsk_431 = buffer.data(lsk + 431);
    const auto *lsk_432 = buffer.data(lsk + 432);
    const auto *lsk_434 = buffer.data(lsk + 434);
    const auto *lsk_435 = buffer.data(lsk + 435);
    const auto *lsk_437 = buffer.data(lsk + 437);
    const auto *lsk_438 = buffer.data(lsk + 438);
    const auto *lsk_441 = buffer.data(lsk + 441);
    const auto *lsk_442 = buffer.data(lsk + 442);
    const auto *lsk_444 = buffer.data(lsk + 444);
    const auto *lsk_446 = buffer.data(lsk + 446);
    const auto *lsk_447 = buffer.data(lsk + 447);
    const auto *lsk_449 = buffer.data(lsk + 449);
    const auto *lsk_450 = buffer.data(lsk + 450);
    const auto *lsk_452 = buffer.data(lsk + 452);
    const auto *lsk_453 = buffer.data(lsk + 453);
    const auto *lsk_455 = buffer.data(lsk + 455);
    const auto *lsk_456 = buffer.data(lsk + 456);
    const auto *lsk_457 = buffer.data(lsk + 457);
    const auto *lsk_459 = buffer.data(lsk + 459);
    const auto *lsk_460 = buffer.data(lsk + 460);
    const auto *lsk_461 = buffer.data(lsk + 461);
    const auto *lsk_462 = buffer.data(lsk + 462);
    const auto *lsk_463 = buffer.data(lsk + 463);
    const auto *lsk_464 = buffer.data(lsk + 464);
    const auto *lsk_465 = buffer.data(lsk + 465);
    const auto *lsk_466 = buffer.data(lsk + 466);
    const auto *lsk_467 = buffer.data(lsk + 467);
    const auto *lsk_468 = buffer.data(lsk + 468);

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, pc_x, pc_z, ksk_388, ksk_390, \
                         ksk_391, ksk_392, lsk_381, lsk_388, lsk_390, lsk_391, \
                         lsk_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_18 * ksk_388[k]
                   + f_3 * pc_x[k] * lsk_388[k];

        t_479[k] = f_3 * pc_z[k] * lsk_381[k];

        t_480[k] = f_18 * ksk_390[k]
                   + f_3 * pc_x[k] * lsk_390[k];

        t_481[k] = f_18 * ksk_391[k]
                   + f_3 * pc_x[k] * lsk_391[k];

        t_482[k] = f_18 * ksk_392[k]
                   + f_3 * pc_x[k] * lsk_392[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pc_x, pc_y, ksk_244, ksk_393, ksk_394, \
                         ksk_395, lsi0_301, lsi1_301, lsk_388, lsk_393, lsk_394, \
                         lsk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_18 * ksk_393[k]
                   + f_3 * pc_x[k] * lsk_393[k];

        t_484[k] = f_18 * ksk_394[k]
                   + f_3 * pc_x[k] * lsk_394[k];

        t_485[k] = f_18 * ksk_395[k]
                   + f_3 * pc_x[k] * lsk_395[k];

        t_486[k] = f_18 * ksk_244[k]
                   + f_1 * lsi0_301[k]
                   - f_2 * lsi1_301[k]
                   + f_3 * pc_y[k] * lsk_388[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pc_z, lsi0_301, lsi0_302, lsi0_303, \
                         lsi1_301, lsi1_302, lsi1_303, lsk_388, lsk_389, lsk_390, \
                         lsk_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_3 * pc_z[k] * lsk_388[k];

        t_488[k] = f_4 * lsi0_301[k]
                   - f_5 * lsi1_301[k]
                   + f_3 * pc_z[k] * lsk_389[k];

        t_489[k] = f_6 * lsi0_302[k]
                   - f_7 * lsi1_302[k]
                   + f_3 * pc_z[k] * lsk_390[k];

        t_490[k] = f_8 * lsi0_303[k]
                   - f_9 * lsi1_303[k]
                   + f_3 * pc_z[k] * lsk_391[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_y, pc_z, ksk_251, lsi0_304, lsi0_305, \
                         lsi0_307, lsi1_304, lsi1_305, lsi1_307, lsk_392, lsk_393, \
                         lsk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_10 * lsi0_304[k]
                   - f_11 * lsi1_304[k]
                   + f_3 * pc_z[k] * lsk_392[k];

        t_492[k] = f_12 * lsi0_305[k]
                   - f_13 * lsi1_305[k]
                   + f_3 * pc_z[k] * lsk_393[k];

        t_493[k] = f_18 * ksk_251[k]
                   + f_3 * pc_y[k] * lsk_395[k];

        t_494[k] = f_1 * lsi0_307[k]
                   - f_2 * lsi1_307[k]
                   + f_3 * pc_z[k] * lsk_395[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pa_z, pc_y, pc_z, ksl0_270, ksl0_273, \
                         ksk_216, ksk_252, ksl1_270, ksl1_273, \
                         lsk_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * ksl0_270[k]
                   - f_14 * pc_z[k] * ksl1_270[k];

        t_496[k] = f_17 * ksk_252[k]
                   + f_3 * pc_y[k] * lsk_396[k];

        t_497[k] = f_15 * ksk_216[k]
                   + f_3 * pc_z[k] * lsk_396[k];

        t_498[k] = pa_z[k] * ksl0_273[k]
                   - f_14 * pc_z[k] * ksl1_273[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pa_z, pc_x, pc_y, pc_z, ksl0_276, ksk_254, \
                         ksk_401, ksl1_276, lsi0_313, lsi1_313, lsk_398, \
                         lsk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_17 * ksk_254[k]
                   + f_3 * pc_y[k] * lsk_398[k];

        t_500[k] = f_18 * ksk_401[k]
                   + f_12 * lsi0_313[k]
                   - f_13 * lsi1_313[k]
                   + f_3 * pc_x[k] * lsk_401[k];

        t_501[k] = pa_z[k] * ksl0_276[k]
                   - f_14 * pc_z[k] * ksl1_276[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, ksk_219, ksk_257, ksk_405, \
                         lsi0_317, lsi1_317, lsk_399, lsk_401, \
                         lsk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_15 * ksk_219[k]
                   + f_3 * pc_z[k] * lsk_399[k];

        t_503[k] = f_17 * ksk_257[k]
                   + f_3 * pc_y[k] * lsk_401[k];

        t_504[k] = f_18 * ksk_405[k]
                   + f_10 * lsi0_317[k]
                   - f_11 * lsi1_317[k]
                   + f_3 * pc_x[k] * lsk_405[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pa_z, pc_y, pc_z, ksl0_280, ksl0_282, \
                         ksk_222, ksk_223, ksk_261, ksl1_280, ksl1_282, lsk_402, \
                         lsk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pa_z[k] * ksl0_280[k]
                   - f_14 * pc_z[k] * ksl1_280[k];

        t_506[k] = f_15 * ksk_222[k]
                   + f_3 * pc_z[k] * lsk_402[k];

        t_507[k] = pa_z[k] * ksl0_282[k]
                   + f_16 * ksk_223[k]
                   - f_14 * pc_z[k] * ksl1_282[k];

        t_508[k] = f_17 * ksk_261[k]
                   + f_3 * pc_y[k] * lsk_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_z, pc_x, pc_z, ksl0_285, ksk_226, ksk_410, \
                         ksl1_285, lsi0_322, lsi1_322, lsk_406, \
                         lsk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_18 * ksk_410[k]
                   + f_8 * lsi0_322[k]
                   - f_9 * lsi1_322[k]
                   + f_3 * pc_x[k] * lsk_410[k];

        t_510[k] = pa_z[k] * ksl0_285[k]
                   - f_14 * pc_z[k] * ksl1_285[k];

        t_511[k] = f_15 * ksk_226[k]
                   + f_3 * pc_z[k] * lsk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pa_z, pc_y, pc_z, ksl0_287, ksl0_288, ksk_227, \
                         ksk_228, ksk_266, ksl1_287, ksl1_288, \
                         lsk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pa_z[k] * ksl0_287[k]
                   + f_16 * ksk_227[k]
                   - f_14 * pc_z[k] * ksl1_287[k];

        t_513[k] = pa_z[k] * ksl0_288[k]
                   + f_17 * ksk_228[k]
                   - f_14 * pc_z[k] * ksl1_288[k];

        t_514[k] = f_17 * ksk_266[k]
                   + f_3 * pc_y[k] * lsk_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pa_z, pc_x, pc_z, ksl0_291, ksk_231, ksk_416, \
                         ksl1_291, lsi0_328, lsi1_328, lsk_411, \
                         lsk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_18 * ksk_416[k]
                   + f_6 * lsi0_328[k]
                   - f_7 * lsi1_328[k]
                   + f_3 * pc_x[k] * lsk_416[k];

        t_516[k] = pa_z[k] * ksl0_291[k]
                   - f_14 * pc_z[k] * ksl1_291[k];

        t_517[k] = f_15 * ksk_231[k]
                   + f_3 * pc_z[k] * lsk_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pa_z, pc_z, ksl0_293, ksl0_294, ksl0_295, \
                         ksk_232, ksk_233, ksk_234, ksl1_293, ksl1_294, \
                         ksl1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pa_z[k] * ksl0_293[k]
                   + f_16 * ksk_232[k]
                   - f_14 * pc_z[k] * ksl1_293[k];

        t_519[k] = pa_z[k] * ksl0_294[k]
                   + f_17 * ksk_233[k]
                   - f_14 * pc_z[k] * ksl1_294[k];

        t_520[k] = pa_z[k] * ksl0_295[k]
                   + f_18 * ksk_234[k]
                   - f_14 * pc_z[k] * ksl1_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, ksk_272, ksk_423, ksk_424, \
                         ksk_425, lsi0_335, lsi1_335, lsk_416, lsk_423, lsk_424, \
                         lsk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * ksk_272[k]
                   + f_3 * pc_y[k] * lsk_416[k];

        t_522[k] = f_18 * ksk_423[k]
                   + f_4 * lsi0_335[k]
                   - f_5 * lsi1_335[k]
                   + f_3 * pc_x[k] * lsk_423[k];

        t_523[k] = f_18 * ksk_424[k]
                   + f_3 * pc_x[k] * lsk_424[k];

        t_524[k] = f_18 * ksk_425[k]
                   + f_3 * pc_x[k] * lsk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, ksk_426, ksk_427, ksk_428, \
                         ksk_429, ksk_430, lsk_426, lsk_427, lsk_428, lsk_429, \
                         lsk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_18 * ksk_426[k]
                   + f_3 * pc_x[k] * lsk_426[k];

        t_526[k] = f_18 * ksk_427[k]
                   + f_3 * pc_x[k] * lsk_427[k];

        t_527[k] = f_18 * ksk_428[k]
                   + f_3 * pc_x[k] * lsk_428[k];

        t_528[k] = f_18 * ksk_429[k]
                   + f_3 * pc_x[k] * lsk_429[k];

        t_529[k] = f_18 * ksk_430[k]
                   + f_3 * pc_x[k] * lsk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pa_z, pc_x, pc_z, ksl0_306, ksk_244, ksk_431, \
                         ksl1_306, lsk_424, lsk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_18 * ksk_431[k]
                   + f_3 * pc_x[k] * lsk_431[k];

        t_531[k] = pa_z[k] * ksl0_306[k]
                   - f_14 * pc_z[k] * ksl1_306[k];

        t_532[k] = f_15 * ksk_244[k]
                   + f_3 * pc_z[k] * lsk_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pc_y, ksk_282, ksk_283, ksk_284, lsi0_331, \
                         lsi0_332, lsi0_333, lsi1_331, lsi1_332, lsi1_333, lsk_426, lsk_427, \
                         lsk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * ksk_282[k]
                   + f_12 * lsi0_331[k]
                   - f_13 * lsi1_331[k]
                   + f_3 * pc_y[k] * lsk_426[k];

        t_534[k] = f_17 * ksk_283[k]
                   + f_10 * lsi0_332[k]
                   - f_11 * lsi1_332[k]
                   + f_3 * pc_y[k] * lsk_427[k];

        t_535[k] = f_17 * ksk_284[k]
                   + f_8 * lsi0_333[k]
                   - f_9 * lsi1_333[k]
                   + f_3 * pc_y[k] * lsk_428[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pc_y, ksk_285, ksk_286, ksk_287, lsi0_334, \
                         lsi0_335, lsi1_334, lsi1_335, lsk_429, lsk_430, \
                         lsk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_17 * ksk_285[k]
                   + f_6 * lsi0_334[k]
                   - f_7 * lsi1_334[k]
                   + f_3 * pc_y[k] * lsk_429[k];

        t_537[k] = f_17 * ksk_286[k]
                   + f_4 * lsi0_335[k]
                   - f_5 * lsi1_335[k]
                   + f_3 * pc_y[k] * lsk_430[k];

        t_538[k] = f_17 * ksk_287[k]
                   + f_3 * pc_y[k] * lsk_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pc_x, pc_y, pc_z, ksk_251, ksk_288, ksk_432, \
                         lsi0_335, lsi0_336, lsi1_335, lsi1_336, lsk_431, \
                         lsk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * ksk_251[k]
                   + f_1 * lsi0_335[k]
                   - f_2 * lsi1_335[k]
                   + f_3 * pc_z[k] * lsk_431[k];

        t_540[k] = f_18 * ksk_432[k]
                   + f_1 * lsi0_336[k]
                   - f_2 * lsi1_336[k]
                   + f_3 * pc_x[k] * lsk_432[k];

        t_541[k] = f_16 * ksk_288[k]
                   + f_3 * pc_y[k] * lsk_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, pc_y, pc_z, ksk_252, ksk_290, ksk_435, \
                         lsi0_339, lsi1_339, lsk_432, lsk_434, \
                         lsk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * ksk_252[k]
                   + f_3 * pc_z[k] * lsk_432[k];

        t_543[k] = f_18 * ksk_435[k]
                   + f_12 * lsi0_339[k]
                   - f_13 * lsi1_339[k]
                   + f_3 * pc_x[k] * lsk_435[k];

        t_544[k] = f_16 * ksk_290[k]
                   + f_3 * pc_y[k] * lsk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, ksk_255, ksk_437, ksk_438, lsi0_341, \
                         lsi0_342, lsi1_341, lsi1_342, lsk_435, lsk_437, \
                         lsk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_18 * ksk_437[k]
                   + f_12 * lsi0_341[k]
                   - f_13 * lsi1_341[k]
                   + f_3 * pc_x[k] * lsk_437[k];

        t_546[k] = f_18 * ksk_438[k]
                   + f_10 * lsi0_342[k]
                   - f_11 * lsi1_342[k]
                   + f_3 * pc_x[k] * lsk_438[k];

        t_547[k] = f_16 * ksk_255[k]
                   + f_3 * pc_z[k] * lsk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pc_x, pc_y, ksk_293, ksk_441, ksk_442, lsi0_345, \
                         lsi0_346, lsi1_345, lsi1_346, lsk_437, lsk_441, \
                         lsk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * ksk_293[k]
                   + f_3 * pc_y[k] * lsk_437[k];

        t_549[k] = f_18 * ksk_441[k]
                   + f_10 * lsi0_345[k]
                   - f_11 * lsi1_345[k]
                   + f_3 * pc_x[k] * lsk_441[k];

        t_550[k] = f_18 * ksk_442[k]
                   + f_8 * lsi0_346[k]
                   - f_9 * lsi1_346[k]
                   + f_3 * pc_x[k] * lsk_442[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pc_x, pc_y, pc_z, ksk_258, ksk_297, ksk_444, \
                         lsi0_348, lsi1_348, lsk_438, lsk_441, \
                         lsk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * ksk_258[k]
                   + f_3 * pc_z[k] * lsk_438[k];

        t_552[k] = f_18 * ksk_444[k]
                   + f_8 * lsi0_348[k]
                   - f_9 * lsi1_348[k]
                   + f_3 * pc_x[k] * lsk_444[k];

        t_553[k] = f_16 * ksk_297[k]
                   + f_3 * pc_y[k] * lsk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_x, pc_z, ksk_262, ksk_446, ksk_447, lsi0_350, \
                         lsi0_351, lsi1_350, lsi1_351, lsk_442, lsk_446, \
                         lsk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_18 * ksk_446[k]
                   + f_8 * lsi0_350[k]
                   - f_9 * lsi1_350[k]
                   + f_3 * pc_x[k] * lsk_446[k];

        t_555[k] = f_18 * ksk_447[k]
                   + f_6 * lsi0_351[k]
                   - f_7 * lsi1_351[k]
                   + f_3 * pc_x[k] * lsk_447[k];

        t_556[k] = f_16 * ksk_262[k]
                   + f_3 * pc_z[k] * lsk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, pc_y, ksk_302, ksk_449, ksk_450, lsi0_353, \
                         lsi0_354, lsi1_353, lsi1_354, lsk_446, lsk_449, \
                         lsk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_18 * ksk_449[k]
                   + f_6 * lsi0_353[k]
                   - f_7 * lsi1_353[k]
                   + f_3 * pc_x[k] * lsk_449[k];

        t_558[k] = f_18 * ksk_450[k]
                   + f_6 * lsi0_354[k]
                   - f_7 * lsi1_354[k]
                   + f_3 * pc_x[k] * lsk_450[k];

        t_559[k] = f_16 * ksk_302[k]
                   + f_3 * pc_y[k] * lsk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pc_x, pc_z, ksk_267, ksk_452, ksk_453, lsi0_356, \
                         lsi0_357, lsi1_356, lsi1_357, lsk_447, lsk_452, \
                         lsk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_18 * ksk_452[k]
                   + f_6 * lsi0_356[k]
                   - f_7 * lsi1_356[k]
                   + f_3 * pc_x[k] * lsk_452[k];

        t_561[k] = f_18 * ksk_453[k]
                   + f_4 * lsi0_357[k]
                   - f_5 * lsi1_357[k]
                   + f_3 * pc_x[k] * lsk_453[k];

        t_562[k] = f_16 * ksk_267[k]
                   + f_3 * pc_z[k] * lsk_447[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, ksk_455, ksk_456, ksk_457, lsi0_359, \
                         lsi0_360, lsi0_361, lsi1_359, lsi1_360, lsi1_361, lsk_455, lsk_456, \
                         lsk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_18 * ksk_455[k]
                   + f_4 * lsi0_359[k]
                   - f_5 * lsi1_359[k]
                   + f_3 * pc_x[k] * lsk_455[k];

        t_564[k] = f_18 * ksk_456[k]
                   + f_4 * lsi0_360[k]
                   - f_5 * lsi1_360[k]
                   + f_3 * pc_x[k] * lsk_456[k];

        t_565[k] = f_18 * ksk_457[k]
                   + f_4 * lsi0_361[k]
                   - f_5 * lsi1_361[k]
                   + f_3 * pc_x[k] * lsk_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pc_x, pc_y, ksk_308, ksk_459, ksk_460, \
                         ksk_461, lsi0_363, lsi1_363, lsk_452, lsk_459, lsk_460, \
                         lsk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * ksk_308[k]
                   + f_3 * pc_y[k] * lsk_452[k];

        t_567[k] = f_18 * ksk_459[k]
                   + f_4 * lsi0_363[k]
                   - f_5 * lsi1_363[k]
                   + f_3 * pc_x[k] * lsk_459[k];

        t_568[k] = f_18 * ksk_460[k]
                   + f_3 * pc_x[k] * lsk_460[k];

        t_569[k] = f_18 * ksk_461[k]
                   + f_3 * pc_x[k] * lsk_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, ksk_462, ksk_463, ksk_464, \
                         ksk_465, ksk_466, lsk_462, lsk_463, lsk_464, lsk_465, \
                         lsk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_18 * ksk_462[k]
                   + f_3 * pc_x[k] * lsk_462[k];

        t_571[k] = f_18 * ksk_463[k]
                   + f_3 * pc_x[k] * lsk_463[k];

        t_572[k] = f_18 * ksk_464[k]
                   + f_3 * pc_x[k] * lsk_464[k];

        t_573[k] = f_18 * ksk_465[k]
                   + f_3 * pc_x[k] * lsk_465[k];

        t_574[k] = f_18 * ksk_466[k]
                   + f_3 * pc_x[k] * lsk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, pc_z, ksk_280, ksk_316, ksk_467, \
                         lsi0_357, lsi1_357, lsk_460, lsk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_18 * ksk_467[k]
                   + f_3 * pc_x[k] * lsk_467[k];

        t_576[k] = f_16 * ksk_316[k]
                   + f_1 * lsi0_357[k]
                   - f_2 * lsi1_357[k]
                   + f_3 * pc_y[k] * lsk_460[k];

        t_577[k] = f_16 * ksk_280[k]
                   + f_3 * pc_z[k] * lsk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, ksk_318, ksk_319, ksk_320, lsi0_359, \
                         lsi0_360, lsi0_361, lsi1_359, lsi1_360, lsi1_361, lsk_462, lsk_463, \
                         lsk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * ksk_318[k]
                   + f_12 * lsi0_359[k]
                   - f_13 * lsi1_359[k]
                   + f_3 * pc_y[k] * lsk_462[k];

        t_579[k] = f_16 * ksk_319[k]
                   + f_10 * lsi0_360[k]
                   - f_11 * lsi1_360[k]
                   + f_3 * pc_y[k] * lsk_463[k];

        t_580[k] = f_16 * ksk_320[k]
                   + f_8 * lsi0_361[k]
                   - f_9 * lsi1_361[k]
                   + f_3 * pc_y[k] * lsk_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, ksk_321, ksk_322, ksk_323, lsi0_362, \
                         lsi0_363, lsi1_362, lsi1_363, lsk_465, lsk_466, \
                         lsk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * ksk_321[k]
                   + f_6 * lsi0_362[k]
                   - f_7 * lsi1_362[k]
                   + f_3 * pc_y[k] * lsk_465[k];

        t_582[k] = f_16 * ksk_322[k]
                   + f_4 * lsi0_363[k]
                   - f_5 * lsi1_363[k]
                   + f_3 * pc_y[k] * lsk_466[k];

        t_583[k] = f_16 * ksk_323[k]
                   + f_3 * pc_y[k] * lsk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_y, pc_y, pc_z, ksl0_405, ksk_287, \
                         ksk_288, ksk_324, ksl1_405, lsi0_363, lsi1_363, lsk_467, \
                         lsk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * ksk_287[k]
                   + f_1 * lsi0_363[k]
                   - f_2 * lsi1_363[k]
                   + f_3 * pc_z[k] * lsk_467[k];

        t_585[k] = pa_y[k] * ksl0_405[k]
                   - f_14 * pc_y[k] * ksl1_405[k];

        t_586[k] = f_15 * ksk_324[k]
                   + f_3 * pc_y[k] * lsk_468[k];

        t_587[k] = f_17 * ksk_288[k]
                   + f_3 * pc_z[k] * lsk_468[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_408 = buffer.data(ksl0 + 408);
    const auto *ksl0_410 = buffer.data(ksl0 + 410);
    const auto *ksl0_411 = buffer.data(ksl0 + 411);
    const auto *ksl0_414 = buffer.data(ksl0 + 414);
    const auto *ksl0_415 = buffer.data(ksl0 + 415);
    const auto *ksl0_417 = buffer.data(ksl0 + 417);
    const auto *ksl0_419 = buffer.data(ksl0 + 419);
    const auto *ksl0_420 = buffer.data(ksl0 + 420);
    const auto *ksl0_422 = buffer.data(ksl0 + 422);
    const auto *ksl0_423 = buffer.data(ksl0 + 423);
    const auto *ksl0_425 = buffer.data(ksl0 + 425);
    const auto *ksl0_426 = buffer.data(ksl0 + 426);
    const auto *ksl0_428 = buffer.data(ksl0 + 428);
    const auto *ksl0_429 = buffer.data(ksl0 + 429);
    const auto *ksl0_430 = buffer.data(ksl0 + 430);
    const auto *ksl0_432 = buffer.data(ksl0 + 432);
    const auto *ksl0_449 = buffer.data(ksl0 + 449);

    const auto *ksk_291 = buffer.data(ksk + 291);
    const auto *ksk_294 = buffer.data(ksk + 294);
    const auto *ksk_298 = buffer.data(ksk + 298);
    const auto *ksk_303 = buffer.data(ksk + 303);
    const auto *ksk_316 = buffer.data(ksk + 316);
    const auto *ksk_324 = buffer.data(ksk + 324);
    const auto *ksk_325 = buffer.data(ksk + 325);
    const auto *ksk_326 = buffer.data(ksk + 326);
    const auto *ksk_327 = buffer.data(ksk + 327);
    const auto *ksk_329 = buffer.data(ksk + 329);
    const auto *ksk_330 = buffer.data(ksk + 330);
    const auto *ksk_332 = buffer.data(ksk + 332);
    const auto *ksk_333 = buffer.data(ksk + 333);
    const auto *ksk_334 = buffer.data(ksk + 334);
    const auto *ksk_336 = buffer.data(ksk + 336);
    const auto *ksk_337 = buffer.data(ksk + 337);
    const auto *ksk_338 = buffer.data(ksk + 338);
    const auto *ksk_339 = buffer.data(ksk + 339);
    const auto *ksk_341 = buffer.data(ksk + 341);
    const auto *ksk_342 = buffer.data(ksk + 342);
    const auto *ksk_343 = buffer.data(ksk + 343);
    const auto *ksk_344 = buffer.data(ksk + 344);
    const auto *ksk_352 = buffer.data(ksk + 352);
    const auto *ksk_354 = buffer.data(ksk + 354);
    const auto *ksk_355 = buffer.data(ksk + 355);
    const auto *ksk_356 = buffer.data(ksk + 356);
    const auto *ksk_357 = buffer.data(ksk + 357);
    const auto *ksk_358 = buffer.data(ksk + 358);
    const auto *ksk_359 = buffer.data(ksk + 359);
    const auto *ksk_360 = buffer.data(ksk + 360);
    const auto *ksk_365 = buffer.data(ksk + 365);
    const auto *ksk_369 = buffer.data(ksk + 369);
    const auto *ksk_374 = buffer.data(ksk + 374);
    const auto *ksk_380 = buffer.data(ksk + 380);
    const auto *ksk_496 = buffer.data(ksk + 496);
    const auto *ksk_497 = buffer.data(ksk + 497);
    const auto *ksk_498 = buffer.data(ksk + 498);
    const auto *ksk_499 = buffer.data(ksk + 499);
    const auto *ksk_500 = buffer.data(ksk + 500);
    const auto *ksk_501 = buffer.data(ksk + 501);
    const auto *ksk_502 = buffer.data(ksk + 502);
    const auto *ksk_503 = buffer.data(ksk + 503);
    const auto *ksk_504 = buffer.data(ksk + 504);
    const auto *ksk_509 = buffer.data(ksk + 509);
    const auto *ksk_513 = buffer.data(ksk + 513);
    const auto *ksk_518 = buffer.data(ksk + 518);
    const auto *ksk_524 = buffer.data(ksk + 524);
    const auto *ksk_531 = buffer.data(ksk + 531);
    const auto *ksk_532 = buffer.data(ksk + 532);
    const auto *ksk_533 = buffer.data(ksk + 533);
    const auto *ksk_534 = buffer.data(ksk + 534);
    const auto *ksk_535 = buffer.data(ksk + 535);
    const auto *ksk_536 = buffer.data(ksk + 536);
    const auto *ksk_537 = buffer.data(ksk + 537);
    const auto *ksk_539 = buffer.data(ksk + 539);
    const auto *ksk_540 = buffer.data(ksk + 540);
    const auto *ksk_543 = buffer.data(ksk + 543);
    const auto *ksk_546 = buffer.data(ksk + 546);
    const auto *ksk_550 = buffer.data(ksk + 550);
    const auto *ksk_555 = buffer.data(ksk + 555);
    const auto *ksk_561 = buffer.data(ksk + 561);

    const auto *ksl1_408 = buffer.data(ksl1 + 408);
    const auto *ksl1_410 = buffer.data(ksl1 + 410);
    const auto *ksl1_411 = buffer.data(ksl1 + 411);
    const auto *ksl1_414 = buffer.data(ksl1 + 414);
    const auto *ksl1_415 = buffer.data(ksl1 + 415);
    const auto *ksl1_417 = buffer.data(ksl1 + 417);
    const auto *ksl1_419 = buffer.data(ksl1 + 419);
    const auto *ksl1_420 = buffer.data(ksl1 + 420);
    const auto *ksl1_422 = buffer.data(ksl1 + 422);
    const auto *ksl1_423 = buffer.data(ksl1 + 423);
    const auto *ksl1_425 = buffer.data(ksl1 + 425);
    const auto *ksl1_426 = buffer.data(ksl1 + 426);
    const auto *ksl1_428 = buffer.data(ksl1 + 428);
    const auto *ksl1_429 = buffer.data(ksl1 + 429);
    const auto *ksl1_430 = buffer.data(ksl1 + 430);
    const auto *ksl1_432 = buffer.data(ksl1 + 432);
    const auto *ksl1_449 = buffer.data(ksl1 + 449);

    const auto *lsi0_385 = buffer.data(lsi0 + 385);
    const auto *lsi0_387 = buffer.data(lsi0 + 387);
    const auto *lsi0_388 = buffer.data(lsi0 + 388);
    const auto *lsi0_389 = buffer.data(lsi0 + 389);
    const auto *lsi0_390 = buffer.data(lsi0 + 390);
    const auto *lsi0_391 = buffer.data(lsi0 + 391);
    const auto *lsi0_392 = buffer.data(lsi0 + 392);
    const auto *lsi0_393 = buffer.data(lsi0 + 393);
    const auto *lsi0_394 = buffer.data(lsi0 + 394);
    const auto *lsi0_395 = buffer.data(lsi0 + 395);
    const auto *lsi0_396 = buffer.data(lsi0 + 396);
    const auto *lsi0_397 = buffer.data(lsi0 + 397);
    const auto *lsi0_398 = buffer.data(lsi0 + 398);
    const auto *lsi0_399 = buffer.data(lsi0 + 399);
    const auto *lsi0_400 = buffer.data(lsi0 + 400);
    const auto *lsi0_401 = buffer.data(lsi0 + 401);
    const auto *lsi0_402 = buffer.data(lsi0 + 402);
    const auto *lsi0_403 = buffer.data(lsi0 + 403);
    const auto *lsi0_404 = buffer.data(lsi0 + 404);
    const auto *lsi0_405 = buffer.data(lsi0 + 405);
    const auto *lsi0_406 = buffer.data(lsi0 + 406);
    const auto *lsi0_412 = buffer.data(lsi0 + 412);
    const auto *lsi0_413 = buffer.data(lsi0 + 413);
    const auto *lsi0_414 = buffer.data(lsi0 + 414);
    const auto *lsi0_415 = buffer.data(lsi0 + 415);
    const auto *lsi0_416 = buffer.data(lsi0 + 416);
    const auto *lsi0_417 = buffer.data(lsi0 + 417);
    const auto *lsi0_418 = buffer.data(lsi0 + 418);
    const auto *lsi0_419 = buffer.data(lsi0 + 419);
    const auto *lsi0_420 = buffer.data(lsi0 + 420);
    const auto *lsi0_422 = buffer.data(lsi0 + 422);
    const auto *lsi0_423 = buffer.data(lsi0 + 423);
    const auto *lsi0_425 = buffer.data(lsi0 + 425);
    const auto *lsi0_426 = buffer.data(lsi0 + 426);
    const auto *lsi0_427 = buffer.data(lsi0 + 427);
    const auto *lsi0_429 = buffer.data(lsi0 + 429);
    const auto *lsi0_430 = buffer.data(lsi0 + 430);
    const auto *lsi0_431 = buffer.data(lsi0 + 431);
    const auto *lsi0_432 = buffer.data(lsi0 + 432);
    const auto *lsi0_434 = buffer.data(lsi0 + 434);
    const auto *lsi0_435 = buffer.data(lsi0 + 435);
    const auto *lsi0_441 = buffer.data(lsi0 + 441);

    const auto *lsi1_385 = buffer.data(lsi1 + 385);
    const auto *lsi1_387 = buffer.data(lsi1 + 387);
    const auto *lsi1_388 = buffer.data(lsi1 + 388);
    const auto *lsi1_389 = buffer.data(lsi1 + 389);
    const auto *lsi1_390 = buffer.data(lsi1 + 390);
    const auto *lsi1_391 = buffer.data(lsi1 + 391);
    const auto *lsi1_392 = buffer.data(lsi1 + 392);
    const auto *lsi1_393 = buffer.data(lsi1 + 393);
    const auto *lsi1_394 = buffer.data(lsi1 + 394);
    const auto *lsi1_395 = buffer.data(lsi1 + 395);
    const auto *lsi1_396 = buffer.data(lsi1 + 396);
    const auto *lsi1_397 = buffer.data(lsi1 + 397);
    const auto *lsi1_398 = buffer.data(lsi1 + 398);
    const auto *lsi1_399 = buffer.data(lsi1 + 399);
    const auto *lsi1_400 = buffer.data(lsi1 + 400);
    const auto *lsi1_401 = buffer.data(lsi1 + 401);
    const auto *lsi1_402 = buffer.data(lsi1 + 402);
    const auto *lsi1_403 = buffer.data(lsi1 + 403);
    const auto *lsi1_404 = buffer.data(lsi1 + 404);
    const auto *lsi1_405 = buffer.data(lsi1 + 405);
    const auto *lsi1_406 = buffer.data(lsi1 + 406);
    const auto *lsi1_412 = buffer.data(lsi1 + 412);
    const auto *lsi1_413 = buffer.data(lsi1 + 413);
    const auto *lsi1_414 = buffer.data(lsi1 + 414);
    const auto *lsi1_415 = buffer.data(lsi1 + 415);
    const auto *lsi1_416 = buffer.data(lsi1 + 416);
    const auto *lsi1_417 = buffer.data(lsi1 + 417);
    const auto *lsi1_418 = buffer.data(lsi1 + 418);
    const auto *lsi1_419 = buffer.data(lsi1 + 419);
    const auto *lsi1_420 = buffer.data(lsi1 + 420);
    const auto *lsi1_422 = buffer.data(lsi1 + 422);
    const auto *lsi1_423 = buffer.data(lsi1 + 423);
    const auto *lsi1_425 = buffer.data(lsi1 + 425);
    const auto *lsi1_426 = buffer.data(lsi1 + 426);
    const auto *lsi1_427 = buffer.data(lsi1 + 427);
    const auto *lsi1_429 = buffer.data(lsi1 + 429);
    const auto *lsi1_430 = buffer.data(lsi1 + 430);
    const auto *lsi1_431 = buffer.data(lsi1 + 431);
    const auto *lsi1_432 = buffer.data(lsi1 + 432);
    const auto *lsi1_434 = buffer.data(lsi1 + 434);
    const auto *lsi1_435 = buffer.data(lsi1 + 435);
    const auto *lsi1_441 = buffer.data(lsi1 + 441);

    const auto *lsk_470 = buffer.data(lsk + 470);
    const auto *lsk_471 = buffer.data(lsk + 471);
    const auto *lsk_473 = buffer.data(lsk + 473);
    const auto *lsk_474 = buffer.data(lsk + 474);
    const auto *lsk_477 = buffer.data(lsk + 477);
    const auto *lsk_478 = buffer.data(lsk + 478);
    const auto *lsk_482 = buffer.data(lsk + 482);
    const auto *lsk_483 = buffer.data(lsk + 483);
    const auto *lsk_488 = buffer.data(lsk + 488);
    const auto *lsk_496 = buffer.data(lsk + 496);
    const auto *lsk_497 = buffer.data(lsk + 497);
    const auto *lsk_498 = buffer.data(lsk + 498);
    const auto *lsk_499 = buffer.data(lsk + 499);
    const auto *lsk_500 = buffer.data(lsk + 500);
    const auto *lsk_501 = buffer.data(lsk + 501);
    const auto *lsk_502 = buffer.data(lsk + 502);
    const auto *lsk_503 = buffer.data(lsk + 503);
    const auto *lsk_504 = buffer.data(lsk + 504);
    const auto *lsk_505 = buffer.data(lsk + 505);
    const auto *lsk_506 = buffer.data(lsk + 506);
    const auto *lsk_507 = buffer.data(lsk + 507);
    const auto *lsk_508 = buffer.data(lsk + 508);
    const auto *lsk_509 = buffer.data(lsk + 509);
    const auto *lsk_510 = buffer.data(lsk + 510);
    const auto *lsk_511 = buffer.data(lsk + 511);
    const auto *lsk_512 = buffer.data(lsk + 512);
    const auto *lsk_513 = buffer.data(lsk + 513);
    const auto *lsk_514 = buffer.data(lsk + 514);
    const auto *lsk_515 = buffer.data(lsk + 515);
    const auto *lsk_516 = buffer.data(lsk + 516);
    const auto *lsk_517 = buffer.data(lsk + 517);
    const auto *lsk_518 = buffer.data(lsk + 518);
    const auto *lsk_519 = buffer.data(lsk + 519);
    const auto *lsk_520 = buffer.data(lsk + 520);
    const auto *lsk_521 = buffer.data(lsk + 521);
    const auto *lsk_522 = buffer.data(lsk + 522);
    const auto *lsk_523 = buffer.data(lsk + 523);
    const auto *lsk_524 = buffer.data(lsk + 524);
    const auto *lsk_531 = buffer.data(lsk + 531);
    const auto *lsk_532 = buffer.data(lsk + 532);
    const auto *lsk_533 = buffer.data(lsk + 533);
    const auto *lsk_534 = buffer.data(lsk + 534);
    const auto *lsk_535 = buffer.data(lsk + 535);
    const auto *lsk_536 = buffer.data(lsk + 536);
    const auto *lsk_537 = buffer.data(lsk + 537);
    const auto *lsk_538 = buffer.data(lsk + 538);
    const auto *lsk_539 = buffer.data(lsk + 539);
    const auto *lsk_540 = buffer.data(lsk + 540);
    const auto *lsk_541 = buffer.data(lsk + 541);
    const auto *lsk_542 = buffer.data(lsk + 542);
    const auto *lsk_543 = buffer.data(lsk + 543);
    const auto *lsk_545 = buffer.data(lsk + 545);
    const auto *lsk_546 = buffer.data(lsk + 546);
    const auto *lsk_547 = buffer.data(lsk + 547);
    const auto *lsk_549 = buffer.data(lsk + 549);
    const auto *lsk_550 = buffer.data(lsk + 550);
    const auto *lsk_551 = buffer.data(lsk + 551);
    const auto *lsk_552 = buffer.data(lsk + 552);
    const auto *lsk_554 = buffer.data(lsk + 554);
    const auto *lsk_555 = buffer.data(lsk + 555);
    const auto *lsk_556 = buffer.data(lsk + 556);
    const auto *lsk_557 = buffer.data(lsk + 557);
    const auto *lsk_558 = buffer.data(lsk + 558);
    const auto *lsk_560 = buffer.data(lsk + 560);
    const auto *lsk_561 = buffer.data(lsk + 561);

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_y, pc_y, ksl0_408, ksl0_410, ksl0_411, \
                         ksk_325, ksk_326, ksk_327, ksl1_408, ksl1_410, ksl1_411, \
                         lsk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_y[k] * ksl0_408[k]
                   + f_16 * ksk_325[k]
                   - f_14 * pc_y[k] * ksl1_408[k];

        t_589[k] = f_15 * ksk_326[k]
                   + f_3 * pc_y[k] * lsk_470[k];

        t_590[k] = pa_y[k] * ksl0_410[k]
                   - f_14 * pc_y[k] * ksl1_410[k];

        t_591[k] = pa_y[k] * ksl0_411[k]
                   + f_17 * ksk_327[k]
                   - f_14 * pc_y[k] * ksl1_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_y, pc_y, pc_z, ksl0_414, ksl0_415, \
                         ksk_291, ksk_329, ksk_330, ksl1_414, ksl1_415, lsk_471, \
                         lsk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * ksk_291[k]
                   + f_3 * pc_z[k] * lsk_471[k];

        t_593[k] = f_15 * ksk_329[k]
                   + f_3 * pc_y[k] * lsk_473[k];

        t_594[k] = pa_y[k] * ksl0_414[k]
                   - f_14 * pc_y[k] * ksl1_414[k];

        t_595[k] = pa_y[k] * ksl0_415[k]
                   + f_18 * ksk_330[k]
                   - f_14 * pc_y[k] * ksl1_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pa_y, pc_y, pc_z, ksl0_417, ksl0_419, \
                         ksk_294, ksk_332, ksk_333, ksl1_417, ksl1_419, lsk_474, \
                         lsk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * ksk_294[k]
                   + f_3 * pc_z[k] * lsk_474[k];

        t_597[k] = pa_y[k] * ksl0_417[k]
                   + f_16 * ksk_332[k]
                   - f_14 * pc_y[k] * ksl1_417[k];

        t_598[k] = f_15 * ksk_333[k]
                   + f_3 * pc_y[k] * lsk_477[k];

        t_599[k] = pa_y[k] * ksl0_419[k]
                   - f_14 * pc_y[k] * ksl1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pa_y, pc_y, pc_z, ksl0_420, ksl0_422, ksk_298, \
                         ksk_334, ksk_336, ksl1_420, ksl1_422, \
                         lsk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pa_y[k] * ksl0_420[k]
                   + f_19 * ksk_334[k]
                   - f_14 * pc_y[k] * ksl1_420[k];

        t_601[k] = f_17 * ksk_298[k]
                   + f_3 * pc_z[k] * lsk_478[k];

        t_602[k] = pa_y[k] * ksl0_422[k]
                   + f_17 * ksk_336[k]
                   - f_14 * pc_y[k] * ksl1_422[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pa_y, pc_y, ksl0_423, ksl0_425, ksl0_426, \
                         ksk_337, ksk_338, ksk_339, ksl1_423, ksl1_425, ksl1_426, \
                         lsk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pa_y[k] * ksl0_423[k]
                   + f_16 * ksk_337[k]
                   - f_14 * pc_y[k] * ksl1_423[k];

        t_604[k] = f_15 * ksk_338[k]
                   + f_3 * pc_y[k] * lsk_482[k];

        t_605[k] = pa_y[k] * ksl0_425[k]
                   - f_14 * pc_y[k] * ksl1_425[k];

        t_606[k] = pa_y[k] * ksl0_426[k]
                   + f_20 * ksk_339[k]
                   - f_14 * pc_y[k] * ksl1_426[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pa_y, pc_y, pc_z, ksl0_428, ksl0_429, ksk_303, \
                         ksk_341, ksk_342, ksl1_428, ksl1_429, \
                         lsk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_17 * ksk_303[k]
                   + f_3 * pc_z[k] * lsk_483[k];

        t_608[k] = pa_y[k] * ksl0_428[k]
                   + f_18 * ksk_341[k]
                   - f_14 * pc_y[k] * ksl1_428[k];

        t_609[k] = pa_y[k] * ksl0_429[k]
                   + f_17 * ksk_342[k]
                   - f_14 * pc_y[k] * ksl1_429[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pa_y, pc_x, pc_y, ksl0_430, ksl0_432, \
                         ksk_343, ksk_344, ksk_496, ksl1_430, ksl1_432, lsk_488, \
                         lsk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pa_y[k] * ksl0_430[k]
                   + f_16 * ksk_343[k]
                   - f_14 * pc_y[k] * ksl1_430[k];

        t_611[k] = f_15 * ksk_344[k]
                   + f_3 * pc_y[k] * lsk_488[k];

        t_612[k] = pa_y[k] * ksl0_432[k]
                   - f_14 * pc_y[k] * ksl1_432[k];

        t_613[k] = f_18 * ksk_496[k]
                   + f_3 * pc_x[k] * lsk_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, ksk_497, ksk_498, ksk_499, \
                         ksk_500, ksk_501, lsk_497, lsk_498, lsk_499, lsk_500, \
                         lsk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_18 * ksk_497[k]
                   + f_3 * pc_x[k] * lsk_497[k];

        t_615[k] = f_18 * ksk_498[k]
                   + f_3 * pc_x[k] * lsk_498[k];

        t_616[k] = f_18 * ksk_499[k]
                   + f_3 * pc_x[k] * lsk_499[k];

        t_617[k] = f_18 * ksk_500[k]
                   + f_3 * pc_x[k] * lsk_500[k];

        t_618[k] = f_18 * ksk_501[k]
                   + f_3 * pc_x[k] * lsk_501[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, ksk_316, ksk_352, \
                         ksk_502, ksk_503, lsi0_385, lsi1_385, lsk_496, lsk_502, \
                         lsk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_18 * ksk_502[k]
                   + f_3 * pc_x[k] * lsk_502[k];

        t_620[k] = f_18 * ksk_503[k]
                   + f_3 * pc_x[k] * lsk_503[k];

        t_621[k] = f_15 * ksk_352[k]
                   + f_1 * lsi0_385[k]
                   - f_2 * lsi1_385[k]
                   + f_3 * pc_y[k] * lsk_496[k];

        t_622[k] = f_17 * ksk_316[k]
                   + f_3 * pc_z[k] * lsk_496[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_y, ksk_354, ksk_355, ksk_356, lsi0_387, \
                         lsi0_388, lsi0_389, lsi1_387, lsi1_388, lsi1_389, lsk_498, lsk_499, \
                         lsk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_15 * ksk_354[k]
                   + f_12 * lsi0_387[k]
                   - f_13 * lsi1_387[k]
                   + f_3 * pc_y[k] * lsk_498[k];

        t_624[k] = f_15 * ksk_355[k]
                   + f_10 * lsi0_388[k]
                   - f_11 * lsi1_388[k]
                   + f_3 * pc_y[k] * lsk_499[k];

        t_625[k] = f_15 * ksk_356[k]
                   + f_8 * lsi0_389[k]
                   - f_9 * lsi1_389[k]
                   + f_3 * pc_y[k] * lsk_500[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, ksk_357, ksk_358, ksk_359, lsi0_390, \
                         lsi0_391, lsi1_390, lsi1_391, lsk_501, lsk_502, \
                         lsk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * ksk_357[k]
                   + f_6 * lsi0_390[k]
                   - f_7 * lsi1_390[k]
                   + f_3 * pc_y[k] * lsk_501[k];

        t_627[k] = f_15 * ksk_358[k]
                   + f_4 * lsi0_391[k]
                   - f_5 * lsi1_391[k]
                   + f_3 * pc_y[k] * lsk_502[k];

        t_628[k] = f_15 * ksk_359[k]
                   + f_3 * pc_y[k] * lsk_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pa_y, pc_x, pc_y, pc_z, ksl0_449, \
                         ksk_324, ksk_504, ksl1_449, lsi0_392, lsi1_392, \
                         lsk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pa_y[k] * ksl0_449[k]
                   - f_14 * pc_y[k] * ksl1_449[k];

        t_630[k] = f_18 * ksk_504[k]
                   + f_1 * lsi0_392[k]
                   - f_2 * lsi1_392[k]
                   + f_3 * pc_x[k] * lsk_504[k];

        t_631[k] = f_3 * pc_y[k] * lsk_504[k];

        t_632[k] = f_18 * ksk_324[k]
                   + f_3 * pc_z[k] * lsk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, ksk_509, lsi0_392, lsi0_397, \
                         lsi1_392, lsi1_397, lsk_505, lsk_506, \
                         lsk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_4 * lsi0_392[k]
                   - f_5 * lsi1_392[k]
                   + f_3 * pc_y[k] * lsk_505[k];

        t_634[k] = f_3 * pc_y[k] * lsk_506[k];

        t_635[k] = f_18 * ksk_509[k]
                   + f_12 * lsi0_397[k]
                   - f_13 * lsi1_397[k]
                   + f_3 * pc_x[k] * lsk_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_y, lsi0_393, lsi0_394, lsi1_393, lsi1_394, \
                         lsk_507, lsk_508, lsk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_6 * lsi0_393[k]
                   - f_7 * lsi1_393[k]
                   + f_3 * pc_y[k] * lsk_507[k];

        t_637[k] = f_4 * lsi0_394[k]
                   - f_5 * lsi1_394[k]
                   + f_3 * pc_y[k] * lsk_508[k];

        t_638[k] = f_3 * pc_y[k] * lsk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_y, ksk_513, lsi0_395, lsi0_396, \
                         lsi0_401, lsi1_395, lsi1_396, lsi1_401, lsk_510, lsk_511, \
                         lsk_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_18 * ksk_513[k]
                   + f_10 * lsi0_401[k]
                   - f_11 * lsi1_401[k]
                   + f_3 * pc_x[k] * lsk_513[k];

        t_640[k] = f_8 * lsi0_395[k]
                   - f_9 * lsi1_395[k]
                   + f_3 * pc_y[k] * lsk_510[k];

        t_641[k] = f_6 * lsi0_396[k]
                   - f_7 * lsi1_396[k]
                   + f_3 * pc_y[k] * lsk_511[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, ksk_518, lsi0_397, lsi0_406, \
                         lsi1_397, lsi1_406, lsk_512, lsk_513, \
                         lsk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_4 * lsi0_397[k]
                   - f_5 * lsi1_397[k]
                   + f_3 * pc_y[k] * lsk_512[k];

        t_643[k] = f_3 * pc_y[k] * lsk_513[k];

        t_644[k] = f_18 * ksk_518[k]
                   + f_8 * lsi0_406[k]
                   - f_9 * lsi1_406[k]
                   + f_3 * pc_x[k] * lsk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_y, lsi0_398, lsi0_399, lsi0_400, lsi1_398, \
                         lsi1_399, lsi1_400, lsk_514, lsk_515, \
                         lsk_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_10 * lsi0_398[k]
                   - f_11 * lsi1_398[k]
                   + f_3 * pc_y[k] * lsk_514[k];

        t_646[k] = f_8 * lsi0_399[k]
                   - f_9 * lsi1_399[k]
                   + f_3 * pc_y[k] * lsk_515[k];

        t_647[k] = f_6 * lsi0_400[k]
                   - f_7 * lsi1_400[k]
                   + f_3 * pc_y[k] * lsk_516[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, pc_y, ksk_524, lsi0_401, lsi0_412, \
                         lsi1_401, lsi1_412, lsk_517, lsk_518, \
                         lsk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_4 * lsi0_401[k]
                   - f_5 * lsi1_401[k]
                   + f_3 * pc_y[k] * lsk_517[k];

        t_649[k] = f_3 * pc_y[k] * lsk_518[k];

        t_650[k] = f_18 * ksk_524[k]
                   + f_6 * lsi0_412[k]
                   - f_7 * lsi1_412[k]
                   + f_3 * pc_x[k] * lsk_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_y, lsi0_402, lsi0_403, lsi0_404, lsi1_402, \
                         lsi1_403, lsi1_404, lsk_519, lsk_520, \
                         lsk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_12 * lsi0_402[k]
                   - f_13 * lsi1_402[k]
                   + f_3 * pc_y[k] * lsk_519[k];

        t_652[k] = f_10 * lsi0_403[k]
                   - f_11 * lsi1_403[k]
                   + f_3 * pc_y[k] * lsk_520[k];

        t_653[k] = f_8 * lsi0_404[k]
                   - f_9 * lsi1_404[k]
                   + f_3 * pc_y[k] * lsk_521[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_y, lsi0_405, lsi0_406, lsi1_405, lsi1_406, \
                         lsk_522, lsk_523, lsk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_6 * lsi0_405[k]
                   - f_7 * lsi1_405[k]
                   + f_3 * pc_y[k] * lsk_522[k];

        t_655[k] = f_4 * lsi0_406[k]
                   - f_5 * lsi1_406[k]
                   + f_3 * pc_y[k] * lsk_523[k];

        t_656[k] = f_3 * pc_y[k] * lsk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pc_x, ksk_531, ksk_532, ksk_533, ksk_534, \
                         lsi0_419, lsi1_419, lsk_531, lsk_532, lsk_533, \
                         lsk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_18 * ksk_531[k]
                   + f_4 * lsi0_419[k]
                   - f_5 * lsi1_419[k]
                   + f_3 * pc_x[k] * lsk_531[k];

        t_658[k] = f_18 * ksk_532[k]
                   + f_3 * pc_x[k] * lsk_532[k];

        t_659[k] = f_18 * ksk_533[k]
                   + f_3 * pc_x[k] * lsk_533[k];

        t_660[k] = f_18 * ksk_534[k]
                   + f_3 * pc_x[k] * lsk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, pc_y, ksk_535, ksk_536, \
                         ksk_537, ksk_539, lsk_531, lsk_535, lsk_536, lsk_537, \
                         lsk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_18 * ksk_535[k]
                   + f_3 * pc_x[k] * lsk_535[k];

        t_662[k] = f_18 * ksk_536[k]
                   + f_3 * pc_x[k] * lsk_536[k];

        t_663[k] = f_18 * ksk_537[k]
                   + f_3 * pc_x[k] * lsk_537[k];

        t_664[k] = f_3 * pc_y[k] * lsk_531[k];

        t_665[k] = f_18 * ksk_539[k]
                   + f_3 * pc_x[k] * lsk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, pc_y, lsi0_413, lsi0_414, lsi0_415, lsi1_413, \
                         lsi1_414, lsi1_415, lsk_532, lsk_533, \
                         lsk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * lsi0_413[k]
                   - f_2 * lsi1_413[k]
                   + f_3 * pc_y[k] * lsk_532[k];

        t_667[k] = f_22 * lsi0_414[k]
                   - f_23 * lsi1_414[k]
                   + f_3 * pc_y[k] * lsk_533[k];

        t_668[k] = f_12 * lsi0_415[k]
                   - f_13 * lsi1_415[k]
                   + f_3 * pc_y[k] * lsk_534[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pc_y, lsi0_416, lsi0_417, lsi0_418, lsi1_416, \
                         lsi1_417, lsi1_418, lsk_535, lsk_536, \
                         lsk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_10 * lsi0_416[k]
                   - f_11 * lsi1_416[k]
                   + f_3 * pc_y[k] * lsk_535[k];

        t_670[k] = f_8 * lsi0_417[k]
                   - f_9 * lsi1_417[k]
                   + f_3 * pc_y[k] * lsk_536[k];

        t_671[k] = f_6 * lsi0_418[k]
                   - f_7 * lsi1_418[k]
                   + f_3 * pc_y[k] * lsk_537[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, pc_y, pc_z, ksk_359, ksk_540, \
                         lsi0_419, lsi0_420, lsi1_419, lsi1_420, lsk_538, lsk_539, \
                         lsk_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_4 * lsi0_419[k]
                   - f_5 * lsi1_419[k]
                   + f_3 * pc_y[k] * lsk_538[k];

        t_673[k] = f_3 * pc_y[k] * lsk_539[k];

        t_674[k] = f_18 * ksk_359[k]
                   + f_1 * lsi0_419[k]
                   - f_2 * lsi1_419[k]
                   + f_3 * pc_z[k] * lsk_539[k];

        t_675[k] = f_17 * ksk_540[k]
                   + f_1 * lsi0_420[k]
                   - f_2 * lsi1_420[k]
                   + f_3 * pc_x[k] * lsk_540[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pc_x, pc_y, pc_z, ksk_360, ksk_543, \
                         lsi0_423, lsi1_423, lsk_540, lsk_541, \
                         lsk_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_19 * ksk_360[k]
                   + f_3 * pc_y[k] * lsk_540[k];

        t_677[k] = f_3 * pc_z[k] * lsk_540[k];

        t_678[k] = f_17 * ksk_543[k]
                   + f_12 * lsi0_423[k]
                   - f_13 * lsi1_423[k]
                   + f_3 * pc_x[k] * lsk_543[k];

        t_679[k] = f_3 * pc_z[k] * lsk_541[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_z, ksk_546, lsi0_420, lsi0_426, \
                         lsi1_420, lsi1_426, lsk_542, lsk_543, \
                         lsk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * lsi0_420[k]
                   - f_5 * lsi1_420[k]
                   + f_3 * pc_z[k] * lsk_542[k];

        t_681[k] = f_17 * ksk_546[k]
                   + f_10 * lsi0_426[k]
                   - f_11 * lsi1_426[k]
                   + f_3 * pc_x[k] * lsk_546[k];

        t_682[k] = f_3 * pc_z[k] * lsk_543[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, pc_x, pc_y, pc_z, ksk_365, ksk_550, \
                         lsi0_422, lsi0_430, lsi1_422, lsi1_430, lsk_545, lsk_546, \
                         lsk_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_19 * ksk_365[k]
                   + f_3 * pc_y[k] * lsk_545[k];

        t_684[k] = f_6 * lsi0_422[k]
                   - f_7 * lsi1_422[k]
                   + f_3 * pc_z[k] * lsk_545[k];

        t_685[k] = f_17 * ksk_550[k]
                   + f_8 * lsi0_430[k]
                   - f_9 * lsi1_430[k]
                   + f_3 * pc_x[k] * lsk_550[k];

        t_686[k] = f_3 * pc_z[k] * lsk_546[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pc_y, pc_z, ksk_369, lsi0_423, lsi0_425, \
                         lsi1_423, lsi1_425, lsk_547, lsk_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_4 * lsi0_423[k]
                   - f_5 * lsi1_423[k]
                   + f_3 * pc_z[k] * lsk_547[k];

        t_688[k] = f_19 * ksk_369[k]
                   + f_3 * pc_y[k] * lsk_549[k];

        t_689[k] = f_8 * lsi0_425[k]
                   - f_9 * lsi1_425[k]
                   + f_3 * pc_z[k] * lsk_549[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, pc_x, pc_z, ksk_555, lsi0_426, lsi0_435, \
                         lsi1_426, lsi1_435, lsk_550, lsk_551, \
                         lsk_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_17 * ksk_555[k]
                   + f_6 * lsi0_435[k]
                   - f_7 * lsi1_435[k]
                   + f_3 * pc_x[k] * lsk_555[k];

        t_691[k] = f_3 * pc_z[k] * lsk_550[k];

        t_692[k] = f_4 * lsi0_426[k]
                   - f_5 * lsi1_426[k]
                   + f_3 * pc_z[k] * lsk_551[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pc_y, pc_z, ksk_374, lsi0_427, lsi0_429, \
                         lsi1_427, lsi1_429, lsk_552, lsk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_6 * lsi0_427[k]
                   - f_7 * lsi1_427[k]
                   + f_3 * pc_z[k] * lsk_552[k];

        t_694[k] = f_19 * ksk_374[k]
                   + f_3 * pc_y[k] * lsk_554[k];

        t_695[k] = f_10 * lsi0_429[k]
                   - f_11 * lsi1_429[k]
                   + f_3 * pc_z[k] * lsk_554[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pc_x, pc_z, ksk_561, lsi0_430, lsi0_441, \
                         lsi1_430, lsi1_441, lsk_555, lsk_556, \
                         lsk_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_17 * ksk_561[k]
                   + f_4 * lsi0_441[k]
                   - f_5 * lsi1_441[k]
                   + f_3 * pc_x[k] * lsk_561[k];

        t_697[k] = f_3 * pc_z[k] * lsk_555[k];

        t_698[k] = f_4 * lsi0_430[k]
                   - f_5 * lsi1_430[k]
                   + f_3 * pc_z[k] * lsk_556[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pc_y, pc_z, ksk_380, lsi0_431, lsi0_432, \
                         lsi0_434, lsi1_431, lsi1_432, lsi1_434, lsk_557, lsk_558, \
                         lsk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_6 * lsi0_431[k]
                   - f_7 * lsi1_431[k]
                   + f_3 * pc_z[k] * lsk_557[k];

        t_700[k] = f_8 * lsi0_432[k]
                   - f_9 * lsi1_432[k]
                   + f_3 * pc_z[k] * lsk_558[k];

        t_701[k] = f_19 * ksk_380[k]
                   + f_3 * pc_y[k] * lsk_560[k];

        t_702[k] = f_12 * lsi0_434[k]
                   - f_13 * lsi1_434[k]
                   + f_3 * pc_z[k] * lsk_560[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_450 = buffer.data(ksl0 + 450);
    const auto *ksl0_453 = buffer.data(ksl0 + 453);
    const auto *ksl0_456 = buffer.data(ksl0 + 456);
    const auto *ksl0_460 = buffer.data(ksl0 + 460);
    const auto *ksl0_462 = buffer.data(ksl0 + 462);
    const auto *ksl0_465 = buffer.data(ksl0 + 465);
    const auto *ksl0_467 = buffer.data(ksl0 + 467);
    const auto *ksl0_468 = buffer.data(ksl0 + 468);
    const auto *ksl0_471 = buffer.data(ksl0 + 471);
    const auto *ksl0_473 = buffer.data(ksl0 + 473);
    const auto *ksl0_474 = buffer.data(ksl0 + 474);
    const auto *ksl0_475 = buffer.data(ksl0 + 475);
    const auto *ksl0_486 = buffer.data(ksl0 + 486);

    const auto *ksk_360 = buffer.data(ksk + 360);
    const auto *ksk_363 = buffer.data(ksk + 363);
    const auto *ksk_366 = buffer.data(ksk + 366);
    const auto *ksk_367 = buffer.data(ksk + 367);
    const auto *ksk_370 = buffer.data(ksk + 370);
    const auto *ksk_371 = buffer.data(ksk + 371);
    const auto *ksk_372 = buffer.data(ksk + 372);
    const auto *ksk_375 = buffer.data(ksk + 375);
    const auto *ksk_376 = buffer.data(ksk + 376);
    const auto *ksk_377 = buffer.data(ksk + 377);
    const auto *ksk_378 = buffer.data(ksk + 378);
    const auto *ksk_388 = buffer.data(ksk + 388);
    const auto *ksk_395 = buffer.data(ksk + 395);
    const auto *ksk_396 = buffer.data(ksk + 396);
    const auto *ksk_398 = buffer.data(ksk + 398);
    const auto *ksk_399 = buffer.data(ksk + 399);
    const auto *ksk_401 = buffer.data(ksk + 401);
    const auto *ksk_402 = buffer.data(ksk + 402);
    const auto *ksk_405 = buffer.data(ksk + 405);
    const auto *ksk_406 = buffer.data(ksk + 406);
    const auto *ksk_410 = buffer.data(ksk + 410);
    const auto *ksk_411 = buffer.data(ksk + 411);
    const auto *ksk_416 = buffer.data(ksk + 416);
    const auto *ksk_424 = buffer.data(ksk + 424);
    const auto *ksk_426 = buffer.data(ksk + 426);
    const auto *ksk_427 = buffer.data(ksk + 427);
    const auto *ksk_428 = buffer.data(ksk + 428);
    const auto *ksk_429 = buffer.data(ksk + 429);
    const auto *ksk_430 = buffer.data(ksk + 430);
    const auto *ksk_431 = buffer.data(ksk + 431);
    const auto *ksk_432 = buffer.data(ksk + 432);
    const auto *ksk_434 = buffer.data(ksk + 434);
    const auto *ksk_437 = buffer.data(ksk + 437);
    const auto *ksk_441 = buffer.data(ksk + 441);
    const auto *ksk_446 = buffer.data(ksk + 446);
    const auto *ksk_452 = buffer.data(ksk + 452);
    const auto *ksk_460 = buffer.data(ksk + 460);
    const auto *ksk_462 = buffer.data(ksk + 462);
    const auto *ksk_463 = buffer.data(ksk + 463);
    const auto *ksk_464 = buffer.data(ksk + 464);
    const auto *ksk_465 = buffer.data(ksk + 465);
    const auto *ksk_466 = buffer.data(ksk + 466);
    const auto *ksk_467 = buffer.data(ksk + 467);
    const auto *ksk_468 = buffer.data(ksk + 468);
    const auto *ksk_568 = buffer.data(ksk + 568);
    const auto *ksk_570 = buffer.data(ksk + 570);
    const auto *ksk_571 = buffer.data(ksk + 571);
    const auto *ksk_572 = buffer.data(ksk + 572);
    const auto *ksk_573 = buffer.data(ksk + 573);
    const auto *ksk_574 = buffer.data(ksk + 574);
    const auto *ksk_575 = buffer.data(ksk + 575);
    const auto *ksk_581 = buffer.data(ksk + 581);
    const auto *ksk_585 = buffer.data(ksk + 585);
    const auto *ksk_590 = buffer.data(ksk + 590);
    const auto *ksk_596 = buffer.data(ksk + 596);
    const auto *ksk_603 = buffer.data(ksk + 603);
    const auto *ksk_604 = buffer.data(ksk + 604);
    const auto *ksk_605 = buffer.data(ksk + 605);
    const auto *ksk_606 = buffer.data(ksk + 606);
    const auto *ksk_607 = buffer.data(ksk + 607);
    const auto *ksk_608 = buffer.data(ksk + 608);
    const auto *ksk_609 = buffer.data(ksk + 609);
    const auto *ksk_610 = buffer.data(ksk + 610);
    const auto *ksk_611 = buffer.data(ksk + 611);
    const auto *ksk_612 = buffer.data(ksk + 612);
    const auto *ksk_615 = buffer.data(ksk + 615);
    const auto *ksk_617 = buffer.data(ksk + 617);
    const auto *ksk_618 = buffer.data(ksk + 618);
    const auto *ksk_621 = buffer.data(ksk + 621);
    const auto *ksk_622 = buffer.data(ksk + 622);
    const auto *ksk_624 = buffer.data(ksk + 624);
    const auto *ksk_626 = buffer.data(ksk + 626);
    const auto *ksk_627 = buffer.data(ksk + 627);
    const auto *ksk_629 = buffer.data(ksk + 629);
    const auto *ksk_630 = buffer.data(ksk + 630);
    const auto *ksk_632 = buffer.data(ksk + 632);
    const auto *ksk_633 = buffer.data(ksk + 633);
    const auto *ksk_635 = buffer.data(ksk + 635);
    const auto *ksk_636 = buffer.data(ksk + 636);
    const auto *ksk_637 = buffer.data(ksk + 637);
    const auto *ksk_639 = buffer.data(ksk + 639);
    const auto *ksk_640 = buffer.data(ksk + 640);
    const auto *ksk_641 = buffer.data(ksk + 641);
    const auto *ksk_642 = buffer.data(ksk + 642);
    const auto *ksk_643 = buffer.data(ksk + 643);
    const auto *ksk_644 = buffer.data(ksk + 644);
    const auto *ksk_645 = buffer.data(ksk + 645);
    const auto *ksk_646 = buffer.data(ksk + 646);
    const auto *ksk_647 = buffer.data(ksk + 647);
    const auto *ksk_648 = buffer.data(ksk + 648);

    const auto *ksl1_450 = buffer.data(ksl1 + 450);
    const auto *ksl1_453 = buffer.data(ksl1 + 453);
    const auto *ksl1_456 = buffer.data(ksl1 + 456);
    const auto *ksl1_460 = buffer.data(ksl1 + 460);
    const auto *ksl1_462 = buffer.data(ksl1 + 462);
    const auto *ksl1_465 = buffer.data(ksl1 + 465);
    const auto *ksl1_467 = buffer.data(ksl1 + 467);
    const auto *ksl1_468 = buffer.data(ksl1 + 468);
    const auto *ksl1_471 = buffer.data(ksl1 + 471);
    const auto *ksl1_473 = buffer.data(ksl1 + 473);
    const auto *ksl1_474 = buffer.data(ksl1 + 474);
    const auto *ksl1_475 = buffer.data(ksl1 + 475);
    const auto *ksl1_486 = buffer.data(ksl1 + 486);

    const auto *lsi0_441 = buffer.data(lsi0 + 441);
    const auto *lsi0_442 = buffer.data(lsi0 + 442);
    const auto *lsi0_443 = buffer.data(lsi0 + 443);
    const auto *lsi0_444 = buffer.data(lsi0 + 444);
    const auto *lsi0_445 = buffer.data(lsi0 + 445);
    const auto *lsi0_447 = buffer.data(lsi0 + 447);
    const auto *lsi0_453 = buffer.data(lsi0 + 453);
    const auto *lsi0_457 = buffer.data(lsi0 + 457);
    const auto *lsi0_462 = buffer.data(lsi0 + 462);
    const auto *lsi0_468 = buffer.data(lsi0 + 468);
    const auto *lsi0_471 = buffer.data(lsi0 + 471);
    const auto *lsi0_472 = buffer.data(lsi0 + 472);
    const auto *lsi0_473 = buffer.data(lsi0 + 473);
    const auto *lsi0_474 = buffer.data(lsi0 + 474);
    const auto *lsi0_475 = buffer.data(lsi0 + 475);
    const auto *lsi0_476 = buffer.data(lsi0 + 476);
    const auto *lsi0_479 = buffer.data(lsi0 + 479);
    const auto *lsi0_481 = buffer.data(lsi0 + 481);
    const auto *lsi0_482 = buffer.data(lsi0 + 482);
    const auto *lsi0_485 = buffer.data(lsi0 + 485);
    const auto *lsi0_486 = buffer.data(lsi0 + 486);
    const auto *lsi0_488 = buffer.data(lsi0 + 488);
    const auto *lsi0_490 = buffer.data(lsi0 + 490);
    const auto *lsi0_491 = buffer.data(lsi0 + 491);
    const auto *lsi0_493 = buffer.data(lsi0 + 493);
    const auto *lsi0_494 = buffer.data(lsi0 + 494);
    const auto *lsi0_496 = buffer.data(lsi0 + 496);
    const auto *lsi0_497 = buffer.data(lsi0 + 497);
    const auto *lsi0_499 = buffer.data(lsi0 + 499);
    const auto *lsi0_500 = buffer.data(lsi0 + 500);
    const auto *lsi0_501 = buffer.data(lsi0 + 501);
    const auto *lsi0_502 = buffer.data(lsi0 + 502);
    const auto *lsi0_503 = buffer.data(lsi0 + 503);
    const auto *lsi0_504 = buffer.data(lsi0 + 504);

    const auto *lsi1_441 = buffer.data(lsi1 + 441);
    const auto *lsi1_442 = buffer.data(lsi1 + 442);
    const auto *lsi1_443 = buffer.data(lsi1 + 443);
    const auto *lsi1_444 = buffer.data(lsi1 + 444);
    const auto *lsi1_445 = buffer.data(lsi1 + 445);
    const auto *lsi1_447 = buffer.data(lsi1 + 447);
    const auto *lsi1_453 = buffer.data(lsi1 + 453);
    const auto *lsi1_457 = buffer.data(lsi1 + 457);
    const auto *lsi1_462 = buffer.data(lsi1 + 462);
    const auto *lsi1_468 = buffer.data(lsi1 + 468);
    const auto *lsi1_471 = buffer.data(lsi1 + 471);
    const auto *lsi1_472 = buffer.data(lsi1 + 472);
    const auto *lsi1_473 = buffer.data(lsi1 + 473);
    const auto *lsi1_474 = buffer.data(lsi1 + 474);
    const auto *lsi1_475 = buffer.data(lsi1 + 475);
    const auto *lsi1_476 = buffer.data(lsi1 + 476);
    const auto *lsi1_479 = buffer.data(lsi1 + 479);
    const auto *lsi1_481 = buffer.data(lsi1 + 481);
    const auto *lsi1_482 = buffer.data(lsi1 + 482);
    const auto *lsi1_485 = buffer.data(lsi1 + 485);
    const auto *lsi1_486 = buffer.data(lsi1 + 486);
    const auto *lsi1_488 = buffer.data(lsi1 + 488);
    const auto *lsi1_490 = buffer.data(lsi1 + 490);
    const auto *lsi1_491 = buffer.data(lsi1 + 491);
    const auto *lsi1_493 = buffer.data(lsi1 + 493);
    const auto *lsi1_494 = buffer.data(lsi1 + 494);
    const auto *lsi1_496 = buffer.data(lsi1 + 496);
    const auto *lsi1_497 = buffer.data(lsi1 + 497);
    const auto *lsi1_499 = buffer.data(lsi1 + 499);
    const auto *lsi1_500 = buffer.data(lsi1 + 500);
    const auto *lsi1_501 = buffer.data(lsi1 + 501);
    const auto *lsi1_502 = buffer.data(lsi1 + 502);
    const auto *lsi1_503 = buffer.data(lsi1 + 503);
    const auto *lsi1_504 = buffer.data(lsi1 + 504);

    const auto *lsk_561 = buffer.data(lsk + 561);
    const auto *lsk_568 = buffer.data(lsk + 568);
    const auto *lsk_569 = buffer.data(lsk + 569);
    const auto *lsk_570 = buffer.data(lsk + 570);
    const auto *lsk_571 = buffer.data(lsk + 571);
    const auto *lsk_572 = buffer.data(lsk + 572);
    const auto *lsk_573 = buffer.data(lsk + 573);
    const auto *lsk_574 = buffer.data(lsk + 574);
    const auto *lsk_575 = buffer.data(lsk + 575);
    const auto *lsk_576 = buffer.data(lsk + 576);
    const auto *lsk_578 = buffer.data(lsk + 578);
    const auto *lsk_579 = buffer.data(lsk + 579);
    const auto *lsk_581 = buffer.data(lsk + 581);
    const auto *lsk_582 = buffer.data(lsk + 582);
    const auto *lsk_585 = buffer.data(lsk + 585);
    const auto *lsk_586 = buffer.data(lsk + 586);
    const auto *lsk_590 = buffer.data(lsk + 590);
    const auto *lsk_591 = buffer.data(lsk + 591);
    const auto *lsk_596 = buffer.data(lsk + 596);
    const auto *lsk_603 = buffer.data(lsk + 603);
    const auto *lsk_604 = buffer.data(lsk + 604);
    const auto *lsk_605 = buffer.data(lsk + 605);
    const auto *lsk_606 = buffer.data(lsk + 606);
    const auto *lsk_607 = buffer.data(lsk + 607);
    const auto *lsk_608 = buffer.data(lsk + 608);
    const auto *lsk_609 = buffer.data(lsk + 609);
    const auto *lsk_610 = buffer.data(lsk + 610);
    const auto *lsk_611 = buffer.data(lsk + 611);
    const auto *lsk_612 = buffer.data(lsk + 612);
    const auto *lsk_614 = buffer.data(lsk + 614);
    const auto *lsk_615 = buffer.data(lsk + 615);
    const auto *lsk_617 = buffer.data(lsk + 617);
    const auto *lsk_618 = buffer.data(lsk + 618);
    const auto *lsk_621 = buffer.data(lsk + 621);
    const auto *lsk_622 = buffer.data(lsk + 622);
    const auto *lsk_624 = buffer.data(lsk + 624);
    const auto *lsk_626 = buffer.data(lsk + 626);
    const auto *lsk_627 = buffer.data(lsk + 627);
    const auto *lsk_629 = buffer.data(lsk + 629);
    const auto *lsk_630 = buffer.data(lsk + 630);
    const auto *lsk_632 = buffer.data(lsk + 632);
    const auto *lsk_633 = buffer.data(lsk + 633);
    const auto *lsk_635 = buffer.data(lsk + 635);
    const auto *lsk_636 = buffer.data(lsk + 636);
    const auto *lsk_637 = buffer.data(lsk + 637);
    const auto *lsk_639 = buffer.data(lsk + 639);
    const auto *lsk_640 = buffer.data(lsk + 640);
    const auto *lsk_641 = buffer.data(lsk + 641);
    const auto *lsk_642 = buffer.data(lsk + 642);
    const auto *lsk_643 = buffer.data(lsk + 643);
    const auto *lsk_644 = buffer.data(lsk + 644);
    const auto *lsk_645 = buffer.data(lsk + 645);
    const auto *lsk_646 = buffer.data(lsk + 646);
    const auto *lsk_647 = buffer.data(lsk + 647);
    const auto *lsk_648 = buffer.data(lsk + 648);

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, pc_x, pc_z, ksk_568, ksk_570, \
                         ksk_571, ksk_572, lsk_561, lsk_568, lsk_570, lsk_571, \
                         lsk_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_17 * ksk_568[k]
                   + f_3 * pc_x[k] * lsk_568[k];

        t_704[k] = f_3 * pc_z[k] * lsk_561[k];

        t_705[k] = f_17 * ksk_570[k]
                   + f_3 * pc_x[k] * lsk_570[k];

        t_706[k] = f_17 * ksk_571[k]
                   + f_3 * pc_x[k] * lsk_571[k];

        t_707[k] = f_17 * ksk_572[k]
                   + f_3 * pc_x[k] * lsk_572[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, pc_x, pc_y, ksk_388, ksk_573, ksk_574, \
                         ksk_575, lsi0_441, lsi1_441, lsk_568, lsk_573, lsk_574, \
                         lsk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_17 * ksk_573[k]
                   + f_3 * pc_x[k] * lsk_573[k];

        t_709[k] = f_17 * ksk_574[k]
                   + f_3 * pc_x[k] * lsk_574[k];

        t_710[k] = f_17 * ksk_575[k]
                   + f_3 * pc_x[k] * lsk_575[k];

        t_711[k] = f_19 * ksk_388[k]
                   + f_1 * lsi0_441[k]
                   - f_2 * lsi1_441[k]
                   + f_3 * pc_y[k] * lsk_568[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pc_z, lsi0_441, lsi0_442, lsi0_443, \
                         lsi1_441, lsi1_442, lsi1_443, lsk_568, lsk_569, lsk_570, \
                         lsk_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_3 * pc_z[k] * lsk_568[k];

        t_713[k] = f_4 * lsi0_441[k]
                   - f_5 * lsi1_441[k]
                   + f_3 * pc_z[k] * lsk_569[k];

        t_714[k] = f_6 * lsi0_442[k]
                   - f_7 * lsi1_442[k]
                   + f_3 * pc_z[k] * lsk_570[k];

        t_715[k] = f_8 * lsi0_443[k]
                   - f_9 * lsi1_443[k]
                   + f_3 * pc_z[k] * lsk_571[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pc_y, pc_z, ksk_395, lsi0_444, lsi0_445, \
                         lsi0_447, lsi1_444, lsi1_445, lsi1_447, lsk_572, lsk_573, \
                         lsk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_10 * lsi0_444[k]
                   - f_11 * lsi1_444[k]
                   + f_3 * pc_z[k] * lsk_572[k];

        t_717[k] = f_12 * lsi0_445[k]
                   - f_13 * lsi1_445[k]
                   + f_3 * pc_z[k] * lsk_573[k];

        t_718[k] = f_19 * ksk_395[k]
                   + f_3 * pc_y[k] * lsk_575[k];

        t_719[k] = f_1 * lsi0_447[k]
                   - f_2 * lsi1_447[k]
                   + f_3 * pc_z[k] * lsk_575[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_z, pc_y, pc_z, ksl0_450, ksl0_453, \
                         ksk_360, ksk_396, ksl1_450, ksl1_453, \
                         lsk_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_z[k] * ksl0_450[k]
                   - f_14 * pc_z[k] * ksl1_450[k];

        t_721[k] = f_18 * ksk_396[k]
                   + f_3 * pc_y[k] * lsk_576[k];

        t_722[k] = f_15 * ksk_360[k]
                   + f_3 * pc_z[k] * lsk_576[k];

        t_723[k] = pa_z[k] * ksl0_453[k]
                   - f_14 * pc_z[k] * ksl1_453[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pa_z, pc_x, pc_y, pc_z, ksl0_456, ksk_398, \
                         ksk_581, ksl1_456, lsi0_453, lsi1_453, lsk_578, \
                         lsk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_18 * ksk_398[k]
                   + f_3 * pc_y[k] * lsk_578[k];

        t_725[k] = f_17 * ksk_581[k]
                   + f_12 * lsi0_453[k]
                   - f_13 * lsi1_453[k]
                   + f_3 * pc_x[k] * lsk_581[k];

        t_726[k] = pa_z[k] * ksl0_456[k]
                   - f_14 * pc_z[k] * ksl1_456[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_x, pc_y, pc_z, ksk_363, ksk_401, ksk_585, \
                         lsi0_457, lsi1_457, lsk_579, lsk_581, \
                         lsk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_15 * ksk_363[k]
                   + f_3 * pc_z[k] * lsk_579[k];

        t_728[k] = f_18 * ksk_401[k]
                   + f_3 * pc_y[k] * lsk_581[k];

        t_729[k] = f_17 * ksk_585[k]
                   + f_10 * lsi0_457[k]
                   - f_11 * lsi1_457[k]
                   + f_3 * pc_x[k] * lsk_585[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_z, pc_y, pc_z, ksl0_460, ksl0_462, \
                         ksk_366, ksk_367, ksk_405, ksl1_460, ksl1_462, lsk_582, \
                         lsk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = pa_z[k] * ksl0_460[k]
                   - f_14 * pc_z[k] * ksl1_460[k];

        t_731[k] = f_15 * ksk_366[k]
                   + f_3 * pc_z[k] * lsk_582[k];

        t_732[k] = pa_z[k] * ksl0_462[k]
                   + f_16 * ksk_367[k]
                   - f_14 * pc_z[k] * ksl1_462[k];

        t_733[k] = f_18 * ksk_405[k]
                   + f_3 * pc_y[k] * lsk_585[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pa_z, pc_x, pc_z, ksl0_465, ksk_370, ksk_590, \
                         ksl1_465, lsi0_462, lsi1_462, lsk_586, \
                         lsk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_17 * ksk_590[k]
                   + f_8 * lsi0_462[k]
                   - f_9 * lsi1_462[k]
                   + f_3 * pc_x[k] * lsk_590[k];

        t_735[k] = pa_z[k] * ksl0_465[k]
                   - f_14 * pc_z[k] * ksl1_465[k];

        t_736[k] = f_15 * ksk_370[k]
                   + f_3 * pc_z[k] * lsk_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pa_z, pc_y, pc_z, ksl0_467, ksl0_468, ksk_371, \
                         ksk_372, ksk_410, ksl1_467, ksl1_468, \
                         lsk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pa_z[k] * ksl0_467[k]
                   + f_16 * ksk_371[k]
                   - f_14 * pc_z[k] * ksl1_467[k];

        t_738[k] = pa_z[k] * ksl0_468[k]
                   + f_17 * ksk_372[k]
                   - f_14 * pc_z[k] * ksl1_468[k];

        t_739[k] = f_18 * ksk_410[k]
                   + f_3 * pc_y[k] * lsk_590[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pa_z, pc_x, pc_z, ksl0_471, ksk_375, ksk_596, \
                         ksl1_471, lsi0_468, lsi1_468, lsk_591, \
                         lsk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_17 * ksk_596[k]
                   + f_6 * lsi0_468[k]
                   - f_7 * lsi1_468[k]
                   + f_3 * pc_x[k] * lsk_596[k];

        t_741[k] = pa_z[k] * ksl0_471[k]
                   - f_14 * pc_z[k] * ksl1_471[k];

        t_742[k] = f_15 * ksk_375[k]
                   + f_3 * pc_z[k] * lsk_591[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, pa_z, pc_z, ksl0_473, ksl0_474, ksl0_475, \
                         ksk_376, ksk_377, ksk_378, ksl1_473, ksl1_474, \
                         ksl1_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = pa_z[k] * ksl0_473[k]
                   + f_16 * ksk_376[k]
                   - f_14 * pc_z[k] * ksl1_473[k];

        t_744[k] = pa_z[k] * ksl0_474[k]
                   + f_17 * ksk_377[k]
                   - f_14 * pc_z[k] * ksl1_474[k];

        t_745[k] = pa_z[k] * ksl0_475[k]
                   + f_18 * ksk_378[k]
                   - f_14 * pc_z[k] * ksl1_475[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, ksk_416, ksk_603, ksk_604, \
                         ksk_605, lsi0_475, lsi1_475, lsk_596, lsk_603, lsk_604, \
                         lsk_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_18 * ksk_416[k]
                   + f_3 * pc_y[k] * lsk_596[k];

        t_747[k] = f_17 * ksk_603[k]
                   + f_4 * lsi0_475[k]
                   - f_5 * lsi1_475[k]
                   + f_3 * pc_x[k] * lsk_603[k];

        t_748[k] = f_17 * ksk_604[k]
                   + f_3 * pc_x[k] * lsk_604[k];

        t_749[k] = f_17 * ksk_605[k]
                   + f_3 * pc_x[k] * lsk_605[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, pc_x, ksk_606, ksk_607, ksk_608, \
                         ksk_609, ksk_610, lsk_606, lsk_607, lsk_608, lsk_609, \
                         lsk_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_17 * ksk_606[k]
                   + f_3 * pc_x[k] * lsk_606[k];

        t_751[k] = f_17 * ksk_607[k]
                   + f_3 * pc_x[k] * lsk_607[k];

        t_752[k] = f_17 * ksk_608[k]
                   + f_3 * pc_x[k] * lsk_608[k];

        t_753[k] = f_17 * ksk_609[k]
                   + f_3 * pc_x[k] * lsk_609[k];

        t_754[k] = f_17 * ksk_610[k]
                   + f_3 * pc_x[k] * lsk_610[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, pa_z, pc_x, pc_z, ksl0_486, ksk_388, ksk_611, \
                         ksl1_486, lsk_604, lsk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_17 * ksk_611[k]
                   + f_3 * pc_x[k] * lsk_611[k];

        t_756[k] = pa_z[k] * ksl0_486[k]
                   - f_14 * pc_z[k] * ksl1_486[k];

        t_757[k] = f_15 * ksk_388[k]
                   + f_3 * pc_z[k] * lsk_604[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_y, ksk_426, ksk_427, ksk_428, lsi0_471, \
                         lsi0_472, lsi0_473, lsi1_471, lsi1_472, lsi1_473, lsk_606, lsk_607, \
                         lsk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_18 * ksk_426[k]
                   + f_12 * lsi0_471[k]
                   - f_13 * lsi1_471[k]
                   + f_3 * pc_y[k] * lsk_606[k];

        t_759[k] = f_18 * ksk_427[k]
                   + f_10 * lsi0_472[k]
                   - f_11 * lsi1_472[k]
                   + f_3 * pc_y[k] * lsk_607[k];

        t_760[k] = f_18 * ksk_428[k]
                   + f_8 * lsi0_473[k]
                   - f_9 * lsi1_473[k]
                   + f_3 * pc_y[k] * lsk_608[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_y, ksk_429, ksk_430, ksk_431, lsi0_474, \
                         lsi0_475, lsi1_474, lsi1_475, lsk_609, lsk_610, \
                         lsk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_18 * ksk_429[k]
                   + f_6 * lsi0_474[k]
                   - f_7 * lsi1_474[k]
                   + f_3 * pc_y[k] * lsk_609[k];

        t_762[k] = f_18 * ksk_430[k]
                   + f_4 * lsi0_475[k]
                   - f_5 * lsi1_475[k]
                   + f_3 * pc_y[k] * lsk_610[k];

        t_763[k] = f_18 * ksk_431[k]
                   + f_3 * pc_y[k] * lsk_611[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pc_x, pc_y, pc_z, ksk_395, ksk_432, ksk_612, \
                         lsi0_475, lsi0_476, lsi1_475, lsi1_476, lsk_611, \
                         lsk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_15 * ksk_395[k]
                   + f_1 * lsi0_475[k]
                   - f_2 * lsi1_475[k]
                   + f_3 * pc_z[k] * lsk_611[k];

        t_765[k] = f_17 * ksk_612[k]
                   + f_1 * lsi0_476[k]
                   - f_2 * lsi1_476[k]
                   + f_3 * pc_x[k] * lsk_612[k];

        t_766[k] = f_17 * ksk_432[k]
                   + f_3 * pc_y[k] * lsk_612[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_x, pc_y, pc_z, ksk_396, ksk_434, ksk_615, \
                         lsi0_479, lsi1_479, lsk_612, lsk_614, \
                         lsk_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_16 * ksk_396[k]
                   + f_3 * pc_z[k] * lsk_612[k];

        t_768[k] = f_17 * ksk_615[k]
                   + f_12 * lsi0_479[k]
                   - f_13 * lsi1_479[k]
                   + f_3 * pc_x[k] * lsk_615[k];

        t_769[k] = f_17 * ksk_434[k]
                   + f_3 * pc_y[k] * lsk_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pc_x, pc_z, ksk_399, ksk_617, ksk_618, lsi0_481, \
                         lsi0_482, lsi1_481, lsi1_482, lsk_615, lsk_617, \
                         lsk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_17 * ksk_617[k]
                   + f_12 * lsi0_481[k]
                   - f_13 * lsi1_481[k]
                   + f_3 * pc_x[k] * lsk_617[k];

        t_771[k] = f_17 * ksk_618[k]
                   + f_10 * lsi0_482[k]
                   - f_11 * lsi1_482[k]
                   + f_3 * pc_x[k] * lsk_618[k];

        t_772[k] = f_16 * ksk_399[k]
                   + f_3 * pc_z[k] * lsk_615[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pc_x, pc_y, ksk_437, ksk_621, ksk_622, lsi0_485, \
                         lsi0_486, lsi1_485, lsi1_486, lsk_617, lsk_621, \
                         lsk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_17 * ksk_437[k]
                   + f_3 * pc_y[k] * lsk_617[k];

        t_774[k] = f_17 * ksk_621[k]
                   + f_10 * lsi0_485[k]
                   - f_11 * lsi1_485[k]
                   + f_3 * pc_x[k] * lsk_621[k];

        t_775[k] = f_17 * ksk_622[k]
                   + f_8 * lsi0_486[k]
                   - f_9 * lsi1_486[k]
                   + f_3 * pc_x[k] * lsk_622[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pc_x, pc_y, pc_z, ksk_402, ksk_441, ksk_624, \
                         lsi0_488, lsi1_488, lsk_618, lsk_621, \
                         lsk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_16 * ksk_402[k]
                   + f_3 * pc_z[k] * lsk_618[k];

        t_777[k] = f_17 * ksk_624[k]
                   + f_8 * lsi0_488[k]
                   - f_9 * lsi1_488[k]
                   + f_3 * pc_x[k] * lsk_624[k];

        t_778[k] = f_17 * ksk_441[k]
                   + f_3 * pc_y[k] * lsk_621[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pc_x, pc_z, ksk_406, ksk_626, ksk_627, lsi0_490, \
                         lsi0_491, lsi1_490, lsi1_491, lsk_622, lsk_626, \
                         lsk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_17 * ksk_626[k]
                   + f_8 * lsi0_490[k]
                   - f_9 * lsi1_490[k]
                   + f_3 * pc_x[k] * lsk_626[k];

        t_780[k] = f_17 * ksk_627[k]
                   + f_6 * lsi0_491[k]
                   - f_7 * lsi1_491[k]
                   + f_3 * pc_x[k] * lsk_627[k];

        t_781[k] = f_16 * ksk_406[k]
                   + f_3 * pc_z[k] * lsk_622[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pc_x, pc_y, ksk_446, ksk_629, ksk_630, lsi0_493, \
                         lsi0_494, lsi1_493, lsi1_494, lsk_626, lsk_629, \
                         lsk_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_17 * ksk_629[k]
                   + f_6 * lsi0_493[k]
                   - f_7 * lsi1_493[k]
                   + f_3 * pc_x[k] * lsk_629[k];

        t_783[k] = f_17 * ksk_630[k]
                   + f_6 * lsi0_494[k]
                   - f_7 * lsi1_494[k]
                   + f_3 * pc_x[k] * lsk_630[k];

        t_784[k] = f_17 * ksk_446[k]
                   + f_3 * pc_y[k] * lsk_626[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pc_x, pc_z, ksk_411, ksk_632, ksk_633, lsi0_496, \
                         lsi0_497, lsi1_496, lsi1_497, lsk_627, lsk_632, \
                         lsk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_17 * ksk_632[k]
                   + f_6 * lsi0_496[k]
                   - f_7 * lsi1_496[k]
                   + f_3 * pc_x[k] * lsk_632[k];

        t_786[k] = f_17 * ksk_633[k]
                   + f_4 * lsi0_497[k]
                   - f_5 * lsi1_497[k]
                   + f_3 * pc_x[k] * lsk_633[k];

        t_787[k] = f_16 * ksk_411[k]
                   + f_3 * pc_z[k] * lsk_627[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, ksk_635, ksk_636, ksk_637, lsi0_499, \
                         lsi0_500, lsi0_501, lsi1_499, lsi1_500, lsi1_501, lsk_635, lsk_636, \
                         lsk_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_17 * ksk_635[k]
                   + f_4 * lsi0_499[k]
                   - f_5 * lsi1_499[k]
                   + f_3 * pc_x[k] * lsk_635[k];

        t_789[k] = f_17 * ksk_636[k]
                   + f_4 * lsi0_500[k]
                   - f_5 * lsi1_500[k]
                   + f_3 * pc_x[k] * lsk_636[k];

        t_790[k] = f_17 * ksk_637[k]
                   + f_4 * lsi0_501[k]
                   - f_5 * lsi1_501[k]
                   + f_3 * pc_x[k] * lsk_637[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pc_x, pc_y, ksk_452, ksk_639, ksk_640, \
                         ksk_641, lsi0_503, lsi1_503, lsk_632, lsk_639, lsk_640, \
                         lsk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_17 * ksk_452[k]
                   + f_3 * pc_y[k] * lsk_632[k];

        t_792[k] = f_17 * ksk_639[k]
                   + f_4 * lsi0_503[k]
                   - f_5 * lsi1_503[k]
                   + f_3 * pc_x[k] * lsk_639[k];

        t_793[k] = f_17 * ksk_640[k]
                   + f_3 * pc_x[k] * lsk_640[k];

        t_794[k] = f_17 * ksk_641[k]
                   + f_3 * pc_x[k] * lsk_641[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, pc_x, ksk_642, ksk_643, ksk_644, \
                         ksk_645, ksk_646, lsk_642, lsk_643, lsk_644, lsk_645, \
                         lsk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_17 * ksk_642[k]
                   + f_3 * pc_x[k] * lsk_642[k];

        t_796[k] = f_17 * ksk_643[k]
                   + f_3 * pc_x[k] * lsk_643[k];

        t_797[k] = f_17 * ksk_644[k]
                   + f_3 * pc_x[k] * lsk_644[k];

        t_798[k] = f_17 * ksk_645[k]
                   + f_3 * pc_x[k] * lsk_645[k];

        t_799[k] = f_17 * ksk_646[k]
                   + f_3 * pc_x[k] * lsk_646[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, ksk_424, ksk_460, ksk_647, \
                         lsi0_497, lsi1_497, lsk_640, lsk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_17 * ksk_647[k]
                   + f_3 * pc_x[k] * lsk_647[k];

        t_801[k] = f_17 * ksk_460[k]
                   + f_1 * lsi0_497[k]
                   - f_2 * lsi1_497[k]
                   + f_3 * pc_y[k] * lsk_640[k];

        t_802[k] = f_16 * ksk_424[k]
                   + f_3 * pc_z[k] * lsk_640[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_y, ksk_462, ksk_463, ksk_464, lsi0_499, \
                         lsi0_500, lsi0_501, lsi1_499, lsi1_500, lsi1_501, lsk_642, lsk_643, \
                         lsk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_17 * ksk_462[k]
                   + f_12 * lsi0_499[k]
                   - f_13 * lsi1_499[k]
                   + f_3 * pc_y[k] * lsk_642[k];

        t_804[k] = f_17 * ksk_463[k]
                   + f_10 * lsi0_500[k]
                   - f_11 * lsi1_500[k]
                   + f_3 * pc_y[k] * lsk_643[k];

        t_805[k] = f_17 * ksk_464[k]
                   + f_8 * lsi0_501[k]
                   - f_9 * lsi1_501[k]
                   + f_3 * pc_y[k] * lsk_644[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pc_y, ksk_465, ksk_466, ksk_467, lsi0_502, \
                         lsi0_503, lsi1_502, lsi1_503, lsk_645, lsk_646, \
                         lsk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * ksk_465[k]
                   + f_6 * lsi0_502[k]
                   - f_7 * lsi1_502[k]
                   + f_3 * pc_y[k] * lsk_645[k];

        t_807[k] = f_17 * ksk_466[k]
                   + f_4 * lsi0_503[k]
                   - f_5 * lsi1_503[k]
                   + f_3 * pc_y[k] * lsk_646[k];

        t_808[k] = f_17 * ksk_467[k]
                   + f_3 * pc_y[k] * lsk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pc_x, pc_y, pc_z, ksk_431, ksk_468, ksk_648, \
                         lsi0_503, lsi0_504, lsi1_503, lsi1_504, lsk_647, \
                         lsk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_16 * ksk_431[k]
                   + f_1 * lsi0_503[k]
                   - f_2 * lsi1_503[k]
                   + f_3 * pc_z[k] * lsk_647[k];

        t_810[k] = f_17 * ksk_648[k]
                   + f_1 * lsi0_504[k]
                   - f_2 * lsi1_504[k]
                   + f_3 * pc_x[k] * lsk_648[k];

        t_811[k] = f_16 * ksk_468[k]
                   + f_3 * pc_y[k] * lsk_648[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_630 = buffer.data(ksl0 + 630);
    const auto *ksl0_633 = buffer.data(ksl0 + 633);
    const auto *ksl0_635 = buffer.data(ksl0 + 635);
    const auto *ksl0_636 = buffer.data(ksl0 + 636);
    const auto *ksl0_639 = buffer.data(ksl0 + 639);
    const auto *ksl0_640 = buffer.data(ksl0 + 640);
    const auto *ksl0_642 = buffer.data(ksl0 + 642);
    const auto *ksl0_644 = buffer.data(ksl0 + 644);
    const auto *ksl0_645 = buffer.data(ksl0 + 645);
    const auto *ksl0_647 = buffer.data(ksl0 + 647);
    const auto *ksl0_648 = buffer.data(ksl0 + 648);
    const auto *ksl0_650 = buffer.data(ksl0 + 650);
    const auto *ksl0_651 = buffer.data(ksl0 + 651);
    const auto *ksl0_653 = buffer.data(ksl0 + 653);
    const auto *ksl0_654 = buffer.data(ksl0 + 654);
    const auto *ksl0_655 = buffer.data(ksl0 + 655);
    const auto *ksl0_657 = buffer.data(ksl0 + 657);
    const auto *ksl0_674 = buffer.data(ksl0 + 674);

    const auto *ksk_432 = buffer.data(ksk + 432);
    const auto *ksk_435 = buffer.data(ksk + 435);
    const auto *ksk_438 = buffer.data(ksk + 438);
    const auto *ksk_442 = buffer.data(ksk + 442);
    const auto *ksk_447 = buffer.data(ksk + 447);
    const auto *ksk_460 = buffer.data(ksk + 460);
    const auto *ksk_467 = buffer.data(ksk + 467);
    const auto *ksk_468 = buffer.data(ksk + 468);
    const auto *ksk_470 = buffer.data(ksk + 470);
    const auto *ksk_471 = buffer.data(ksk + 471);
    const auto *ksk_473 = buffer.data(ksk + 473);
    const auto *ksk_474 = buffer.data(ksk + 474);
    const auto *ksk_477 = buffer.data(ksk + 477);
    const auto *ksk_478 = buffer.data(ksk + 478);
    const auto *ksk_482 = buffer.data(ksk + 482);
    const auto *ksk_483 = buffer.data(ksk + 483);
    const auto *ksk_488 = buffer.data(ksk + 488);
    const auto *ksk_496 = buffer.data(ksk + 496);
    const auto *ksk_498 = buffer.data(ksk + 498);
    const auto *ksk_499 = buffer.data(ksk + 499);
    const auto *ksk_500 = buffer.data(ksk + 500);
    const auto *ksk_501 = buffer.data(ksk + 501);
    const auto *ksk_502 = buffer.data(ksk + 502);
    const auto *ksk_503 = buffer.data(ksk + 503);
    const auto *ksk_504 = buffer.data(ksk + 504);
    const auto *ksk_505 = buffer.data(ksk + 505);
    const auto *ksk_506 = buffer.data(ksk + 506);
    const auto *ksk_507 = buffer.data(ksk + 507);
    const auto *ksk_509 = buffer.data(ksk + 509);
    const auto *ksk_510 = buffer.data(ksk + 510);
    const auto *ksk_512 = buffer.data(ksk + 512);
    const auto *ksk_513 = buffer.data(ksk + 513);
    const auto *ksk_514 = buffer.data(ksk + 514);
    const auto *ksk_516 = buffer.data(ksk + 516);
    const auto *ksk_517 = buffer.data(ksk + 517);
    const auto *ksk_518 = buffer.data(ksk + 518);
    const auto *ksk_519 = buffer.data(ksk + 519);
    const auto *ksk_521 = buffer.data(ksk + 521);
    const auto *ksk_522 = buffer.data(ksk + 522);
    const auto *ksk_523 = buffer.data(ksk + 523);
    const auto *ksk_524 = buffer.data(ksk + 524);
    const auto *ksk_532 = buffer.data(ksk + 532);
    const auto *ksk_534 = buffer.data(ksk + 534);
    const auto *ksk_535 = buffer.data(ksk + 535);
    const auto *ksk_536 = buffer.data(ksk + 536);
    const auto *ksk_537 = buffer.data(ksk + 537);
    const auto *ksk_538 = buffer.data(ksk + 538);
    const auto *ksk_539 = buffer.data(ksk + 539);
    const auto *ksk_651 = buffer.data(ksk + 651);
    const auto *ksk_653 = buffer.data(ksk + 653);
    const auto *ksk_654 = buffer.data(ksk + 654);
    const auto *ksk_657 = buffer.data(ksk + 657);
    const auto *ksk_658 = buffer.data(ksk + 658);
    const auto *ksk_660 = buffer.data(ksk + 660);
    const auto *ksk_662 = buffer.data(ksk + 662);
    const auto *ksk_663 = buffer.data(ksk + 663);
    const auto *ksk_665 = buffer.data(ksk + 665);
    const auto *ksk_666 = buffer.data(ksk + 666);
    const auto *ksk_668 = buffer.data(ksk + 668);
    const auto *ksk_669 = buffer.data(ksk + 669);
    const auto *ksk_671 = buffer.data(ksk + 671);
    const auto *ksk_672 = buffer.data(ksk + 672);
    const auto *ksk_673 = buffer.data(ksk + 673);
    const auto *ksk_675 = buffer.data(ksk + 675);
    const auto *ksk_676 = buffer.data(ksk + 676);
    const auto *ksk_677 = buffer.data(ksk + 677);
    const auto *ksk_678 = buffer.data(ksk + 678);
    const auto *ksk_679 = buffer.data(ksk + 679);
    const auto *ksk_680 = buffer.data(ksk + 680);
    const auto *ksk_681 = buffer.data(ksk + 681);
    const auto *ksk_682 = buffer.data(ksk + 682);
    const auto *ksk_683 = buffer.data(ksk + 683);
    const auto *ksk_712 = buffer.data(ksk + 712);
    const auto *ksk_713 = buffer.data(ksk + 713);
    const auto *ksk_714 = buffer.data(ksk + 714);
    const auto *ksk_715 = buffer.data(ksk + 715);
    const auto *ksk_716 = buffer.data(ksk + 716);
    const auto *ksk_717 = buffer.data(ksk + 717);
    const auto *ksk_718 = buffer.data(ksk + 718);
    const auto *ksk_719 = buffer.data(ksk + 719);
    const auto *ksk_720 = buffer.data(ksk + 720);
    const auto *ksk_725 = buffer.data(ksk + 725);
    const auto *ksk_729 = buffer.data(ksk + 729);
    const auto *ksk_734 = buffer.data(ksk + 734);
    const auto *ksk_740 = buffer.data(ksk + 740);

    const auto *ksl1_630 = buffer.data(ksl1 + 630);
    const auto *ksl1_633 = buffer.data(ksl1 + 633);
    const auto *ksl1_635 = buffer.data(ksl1 + 635);
    const auto *ksl1_636 = buffer.data(ksl1 + 636);
    const auto *ksl1_639 = buffer.data(ksl1 + 639);
    const auto *ksl1_640 = buffer.data(ksl1 + 640);
    const auto *ksl1_642 = buffer.data(ksl1 + 642);
    const auto *ksl1_644 = buffer.data(ksl1 + 644);
    const auto *ksl1_645 = buffer.data(ksl1 + 645);
    const auto *ksl1_647 = buffer.data(ksl1 + 647);
    const auto *ksl1_648 = buffer.data(ksl1 + 648);
    const auto *ksl1_650 = buffer.data(ksl1 + 650);
    const auto *ksl1_651 = buffer.data(ksl1 + 651);
    const auto *ksl1_653 = buffer.data(ksl1 + 653);
    const auto *ksl1_654 = buffer.data(ksl1 + 654);
    const auto *ksl1_655 = buffer.data(ksl1 + 655);
    const auto *ksl1_657 = buffer.data(ksl1 + 657);
    const auto *ksl1_674 = buffer.data(ksl1 + 674);

    const auto *lsi0_507 = buffer.data(lsi0 + 507);
    const auto *lsi0_509 = buffer.data(lsi0 + 509);
    const auto *lsi0_510 = buffer.data(lsi0 + 510);
    const auto *lsi0_513 = buffer.data(lsi0 + 513);
    const auto *lsi0_514 = buffer.data(lsi0 + 514);
    const auto *lsi0_516 = buffer.data(lsi0 + 516);
    const auto *lsi0_518 = buffer.data(lsi0 + 518);
    const auto *lsi0_519 = buffer.data(lsi0 + 519);
    const auto *lsi0_521 = buffer.data(lsi0 + 521);
    const auto *lsi0_522 = buffer.data(lsi0 + 522);
    const auto *lsi0_524 = buffer.data(lsi0 + 524);
    const auto *lsi0_525 = buffer.data(lsi0 + 525);
    const auto *lsi0_527 = buffer.data(lsi0 + 527);
    const auto *lsi0_528 = buffer.data(lsi0 + 528);
    const auto *lsi0_529 = buffer.data(lsi0 + 529);
    const auto *lsi0_530 = buffer.data(lsi0 + 530);
    const auto *lsi0_531 = buffer.data(lsi0 + 531);
    const auto *lsi0_553 = buffer.data(lsi0 + 553);
    const auto *lsi0_555 = buffer.data(lsi0 + 555);
    const auto *lsi0_556 = buffer.data(lsi0 + 556);
    const auto *lsi0_557 = buffer.data(lsi0 + 557);
    const auto *lsi0_558 = buffer.data(lsi0 + 558);
    const auto *lsi0_559 = buffer.data(lsi0 + 559);
    const auto *lsi0_560 = buffer.data(lsi0 + 560);
    const auto *lsi0_561 = buffer.data(lsi0 + 561);
    const auto *lsi0_562 = buffer.data(lsi0 + 562);
    const auto *lsi0_563 = buffer.data(lsi0 + 563);
    const auto *lsi0_564 = buffer.data(lsi0 + 564);
    const auto *lsi0_565 = buffer.data(lsi0 + 565);
    const auto *lsi0_566 = buffer.data(lsi0 + 566);
    const auto *lsi0_567 = buffer.data(lsi0 + 567);
    const auto *lsi0_568 = buffer.data(lsi0 + 568);
    const auto *lsi0_569 = buffer.data(lsi0 + 569);
    const auto *lsi0_574 = buffer.data(lsi0 + 574);
    const auto *lsi0_580 = buffer.data(lsi0 + 580);

    const auto *lsi1_507 = buffer.data(lsi1 + 507);
    const auto *lsi1_509 = buffer.data(lsi1 + 509);
    const auto *lsi1_510 = buffer.data(lsi1 + 510);
    const auto *lsi1_513 = buffer.data(lsi1 + 513);
    const auto *lsi1_514 = buffer.data(lsi1 + 514);
    const auto *lsi1_516 = buffer.data(lsi1 + 516);
    const auto *lsi1_518 = buffer.data(lsi1 + 518);
    const auto *lsi1_519 = buffer.data(lsi1 + 519);
    const auto *lsi1_521 = buffer.data(lsi1 + 521);
    const auto *lsi1_522 = buffer.data(lsi1 + 522);
    const auto *lsi1_524 = buffer.data(lsi1 + 524);
    const auto *lsi1_525 = buffer.data(lsi1 + 525);
    const auto *lsi1_527 = buffer.data(lsi1 + 527);
    const auto *lsi1_528 = buffer.data(lsi1 + 528);
    const auto *lsi1_529 = buffer.data(lsi1 + 529);
    const auto *lsi1_530 = buffer.data(lsi1 + 530);
    const auto *lsi1_531 = buffer.data(lsi1 + 531);
    const auto *lsi1_553 = buffer.data(lsi1 + 553);
    const auto *lsi1_555 = buffer.data(lsi1 + 555);
    const auto *lsi1_556 = buffer.data(lsi1 + 556);
    const auto *lsi1_557 = buffer.data(lsi1 + 557);
    const auto *lsi1_558 = buffer.data(lsi1 + 558);
    const auto *lsi1_559 = buffer.data(lsi1 + 559);
    const auto *lsi1_560 = buffer.data(lsi1 + 560);
    const auto *lsi1_561 = buffer.data(lsi1 + 561);
    const auto *lsi1_562 = buffer.data(lsi1 + 562);
    const auto *lsi1_563 = buffer.data(lsi1 + 563);
    const auto *lsi1_564 = buffer.data(lsi1 + 564);
    const auto *lsi1_565 = buffer.data(lsi1 + 565);
    const auto *lsi1_566 = buffer.data(lsi1 + 566);
    const auto *lsi1_567 = buffer.data(lsi1 + 567);
    const auto *lsi1_568 = buffer.data(lsi1 + 568);
    const auto *lsi1_569 = buffer.data(lsi1 + 569);
    const auto *lsi1_574 = buffer.data(lsi1 + 574);
    const auto *lsi1_580 = buffer.data(lsi1 + 580);

    const auto *lsk_648 = buffer.data(lsk + 648);
    const auto *lsk_650 = buffer.data(lsk + 650);
    const auto *lsk_651 = buffer.data(lsk + 651);
    const auto *lsk_653 = buffer.data(lsk + 653);
    const auto *lsk_654 = buffer.data(lsk + 654);
    const auto *lsk_657 = buffer.data(lsk + 657);
    const auto *lsk_658 = buffer.data(lsk + 658);
    const auto *lsk_660 = buffer.data(lsk + 660);
    const auto *lsk_662 = buffer.data(lsk + 662);
    const auto *lsk_663 = buffer.data(lsk + 663);
    const auto *lsk_665 = buffer.data(lsk + 665);
    const auto *lsk_666 = buffer.data(lsk + 666);
    const auto *lsk_668 = buffer.data(lsk + 668);
    const auto *lsk_669 = buffer.data(lsk + 669);
    const auto *lsk_671 = buffer.data(lsk + 671);
    const auto *lsk_672 = buffer.data(lsk + 672);
    const auto *lsk_673 = buffer.data(lsk + 673);
    const auto *lsk_675 = buffer.data(lsk + 675);
    const auto *lsk_676 = buffer.data(lsk + 676);
    const auto *lsk_677 = buffer.data(lsk + 677);
    const auto *lsk_678 = buffer.data(lsk + 678);
    const auto *lsk_679 = buffer.data(lsk + 679);
    const auto *lsk_680 = buffer.data(lsk + 680);
    const auto *lsk_681 = buffer.data(lsk + 681);
    const auto *lsk_682 = buffer.data(lsk + 682);
    const auto *lsk_683 = buffer.data(lsk + 683);
    const auto *lsk_684 = buffer.data(lsk + 684);
    const auto *lsk_686 = buffer.data(lsk + 686);
    const auto *lsk_687 = buffer.data(lsk + 687);
    const auto *lsk_689 = buffer.data(lsk + 689);
    const auto *lsk_690 = buffer.data(lsk + 690);
    const auto *lsk_693 = buffer.data(lsk + 693);
    const auto *lsk_694 = buffer.data(lsk + 694);
    const auto *lsk_698 = buffer.data(lsk + 698);
    const auto *lsk_699 = buffer.data(lsk + 699);
    const auto *lsk_704 = buffer.data(lsk + 704);
    const auto *lsk_712 = buffer.data(lsk + 712);
    const auto *lsk_713 = buffer.data(lsk + 713);
    const auto *lsk_714 = buffer.data(lsk + 714);
    const auto *lsk_715 = buffer.data(lsk + 715);
    const auto *lsk_716 = buffer.data(lsk + 716);
    const auto *lsk_717 = buffer.data(lsk + 717);
    const auto *lsk_718 = buffer.data(lsk + 718);
    const auto *lsk_719 = buffer.data(lsk + 719);
    const auto *lsk_720 = buffer.data(lsk + 720);
    const auto *lsk_721 = buffer.data(lsk + 721);
    const auto *lsk_722 = buffer.data(lsk + 722);
    const auto *lsk_723 = buffer.data(lsk + 723);
    const auto *lsk_724 = buffer.data(lsk + 724);
    const auto *lsk_725 = buffer.data(lsk + 725);
    const auto *lsk_726 = buffer.data(lsk + 726);
    const auto *lsk_727 = buffer.data(lsk + 727);
    const auto *lsk_728 = buffer.data(lsk + 728);
    const auto *lsk_729 = buffer.data(lsk + 729);
    const auto *lsk_730 = buffer.data(lsk + 730);
    const auto *lsk_731 = buffer.data(lsk + 731);
    const auto *lsk_732 = buffer.data(lsk + 732);
    const auto *lsk_733 = buffer.data(lsk + 733);
    const auto *lsk_734 = buffer.data(lsk + 734);
    const auto *lsk_740 = buffer.data(lsk + 740);

#pragma omp simd aligned(t_812, t_813, t_814, pc_x, pc_y, pc_z, ksk_432, ksk_470, ksk_651, \
                         lsi0_507, lsi1_507, lsk_648, lsk_650, \
                         lsk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_17 * ksk_432[k]
                   + f_3 * pc_z[k] * lsk_648[k];

        t_813[k] = f_17 * ksk_651[k]
                   + f_12 * lsi0_507[k]
                   - f_13 * lsi1_507[k]
                   + f_3 * pc_x[k] * lsk_651[k];

        t_814[k] = f_16 * ksk_470[k]
                   + f_3 * pc_y[k] * lsk_650[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pc_x, pc_z, ksk_435, ksk_653, ksk_654, lsi0_509, \
                         lsi0_510, lsi1_509, lsi1_510, lsk_651, lsk_653, \
                         lsk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_17 * ksk_653[k]
                   + f_12 * lsi0_509[k]
                   - f_13 * lsi1_509[k]
                   + f_3 * pc_x[k] * lsk_653[k];

        t_816[k] = f_17 * ksk_654[k]
                   + f_10 * lsi0_510[k]
                   - f_11 * lsi1_510[k]
                   + f_3 * pc_x[k] * lsk_654[k];

        t_817[k] = f_17 * ksk_435[k]
                   + f_3 * pc_z[k] * lsk_651[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pc_x, pc_y, ksk_473, ksk_657, ksk_658, lsi0_513, \
                         lsi0_514, lsi1_513, lsi1_514, lsk_653, lsk_657, \
                         lsk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_16 * ksk_473[k]
                   + f_3 * pc_y[k] * lsk_653[k];

        t_819[k] = f_17 * ksk_657[k]
                   + f_10 * lsi0_513[k]
                   - f_11 * lsi1_513[k]
                   + f_3 * pc_x[k] * lsk_657[k];

        t_820[k] = f_17 * ksk_658[k]
                   + f_8 * lsi0_514[k]
                   - f_9 * lsi1_514[k]
                   + f_3 * pc_x[k] * lsk_658[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pc_x, pc_y, pc_z, ksk_438, ksk_477, ksk_660, \
                         lsi0_516, lsi1_516, lsk_654, lsk_657, \
                         lsk_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_17 * ksk_438[k]
                   + f_3 * pc_z[k] * lsk_654[k];

        t_822[k] = f_17 * ksk_660[k]
                   + f_8 * lsi0_516[k]
                   - f_9 * lsi1_516[k]
                   + f_3 * pc_x[k] * lsk_660[k];

        t_823[k] = f_16 * ksk_477[k]
                   + f_3 * pc_y[k] * lsk_657[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, ksk_442, ksk_662, ksk_663, lsi0_518, \
                         lsi0_519, lsi1_518, lsi1_519, lsk_658, lsk_662, \
                         lsk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_17 * ksk_662[k]
                   + f_8 * lsi0_518[k]
                   - f_9 * lsi1_518[k]
                   + f_3 * pc_x[k] * lsk_662[k];

        t_825[k] = f_17 * ksk_663[k]
                   + f_6 * lsi0_519[k]
                   - f_7 * lsi1_519[k]
                   + f_3 * pc_x[k] * lsk_663[k];

        t_826[k] = f_17 * ksk_442[k]
                   + f_3 * pc_z[k] * lsk_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, pc_x, pc_y, ksk_482, ksk_665, ksk_666, lsi0_521, \
                         lsi0_522, lsi1_521, lsi1_522, lsk_662, lsk_665, \
                         lsk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_17 * ksk_665[k]
                   + f_6 * lsi0_521[k]
                   - f_7 * lsi1_521[k]
                   + f_3 * pc_x[k] * lsk_665[k];

        t_828[k] = f_17 * ksk_666[k]
                   + f_6 * lsi0_522[k]
                   - f_7 * lsi1_522[k]
                   + f_3 * pc_x[k] * lsk_666[k];

        t_829[k] = f_16 * ksk_482[k]
                   + f_3 * pc_y[k] * lsk_662[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, pc_x, pc_z, ksk_447, ksk_668, ksk_669, lsi0_524, \
                         lsi0_525, lsi1_524, lsi1_525, lsk_663, lsk_668, \
                         lsk_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_17 * ksk_668[k]
                   + f_6 * lsi0_524[k]
                   - f_7 * lsi1_524[k]
                   + f_3 * pc_x[k] * lsk_668[k];

        t_831[k] = f_17 * ksk_669[k]
                   + f_4 * lsi0_525[k]
                   - f_5 * lsi1_525[k]
                   + f_3 * pc_x[k] * lsk_669[k];

        t_832[k] = f_17 * ksk_447[k]
                   + f_3 * pc_z[k] * lsk_663[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, ksk_671, ksk_672, ksk_673, lsi0_527, \
                         lsi0_528, lsi0_529, lsi1_527, lsi1_528, lsi1_529, lsk_671, lsk_672, \
                         lsk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_17 * ksk_671[k]
                   + f_4 * lsi0_527[k]
                   - f_5 * lsi1_527[k]
                   + f_3 * pc_x[k] * lsk_671[k];

        t_834[k] = f_17 * ksk_672[k]
                   + f_4 * lsi0_528[k]
                   - f_5 * lsi1_528[k]
                   + f_3 * pc_x[k] * lsk_672[k];

        t_835[k] = f_17 * ksk_673[k]
                   + f_4 * lsi0_529[k]
                   - f_5 * lsi1_529[k]
                   + f_3 * pc_x[k] * lsk_673[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pc_x, pc_y, ksk_488, ksk_675, ksk_676, \
                         ksk_677, lsi0_531, lsi1_531, lsk_668, lsk_675, lsk_676, \
                         lsk_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * ksk_488[k]
                   + f_3 * pc_y[k] * lsk_668[k];

        t_837[k] = f_17 * ksk_675[k]
                   + f_4 * lsi0_531[k]
                   - f_5 * lsi1_531[k]
                   + f_3 * pc_x[k] * lsk_675[k];

        t_838[k] = f_17 * ksk_676[k]
                   + f_3 * pc_x[k] * lsk_676[k];

        t_839[k] = f_17 * ksk_677[k]
                   + f_3 * pc_x[k] * lsk_677[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, pc_x, ksk_678, ksk_679, ksk_680, \
                         ksk_681, ksk_682, lsk_678, lsk_679, lsk_680, lsk_681, \
                         lsk_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_17 * ksk_678[k]
                   + f_3 * pc_x[k] * lsk_678[k];

        t_841[k] = f_17 * ksk_679[k]
                   + f_3 * pc_x[k] * lsk_679[k];

        t_842[k] = f_17 * ksk_680[k]
                   + f_3 * pc_x[k] * lsk_680[k];

        t_843[k] = f_17 * ksk_681[k]
                   + f_3 * pc_x[k] * lsk_681[k];

        t_844[k] = f_17 * ksk_682[k]
                   + f_3 * pc_x[k] * lsk_682[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, pc_z, ksk_460, ksk_496, ksk_683, \
                         lsi0_525, lsi1_525, lsk_676, lsk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_17 * ksk_683[k]
                   + f_3 * pc_x[k] * lsk_683[k];

        t_846[k] = f_16 * ksk_496[k]
                   + f_1 * lsi0_525[k]
                   - f_2 * lsi1_525[k]
                   + f_3 * pc_y[k] * lsk_676[k];

        t_847[k] = f_17 * ksk_460[k]
                   + f_3 * pc_z[k] * lsk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pc_y, ksk_498, ksk_499, ksk_500, lsi0_527, \
                         lsi0_528, lsi0_529, lsi1_527, lsi1_528, lsi1_529, lsk_678, lsk_679, \
                         lsk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * ksk_498[k]
                   + f_12 * lsi0_527[k]
                   - f_13 * lsi1_527[k]
                   + f_3 * pc_y[k] * lsk_678[k];

        t_849[k] = f_16 * ksk_499[k]
                   + f_10 * lsi0_528[k]
                   - f_11 * lsi1_528[k]
                   + f_3 * pc_y[k] * lsk_679[k];

        t_850[k] = f_16 * ksk_500[k]
                   + f_8 * lsi0_529[k]
                   - f_9 * lsi1_529[k]
                   + f_3 * pc_y[k] * lsk_680[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pc_y, ksk_501, ksk_502, ksk_503, lsi0_530, \
                         lsi0_531, lsi1_530, lsi1_531, lsk_681, lsk_682, \
                         lsk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_16 * ksk_501[k]
                   + f_6 * lsi0_530[k]
                   - f_7 * lsi1_530[k]
                   + f_3 * pc_y[k] * lsk_681[k];

        t_852[k] = f_16 * ksk_502[k]
                   + f_4 * lsi0_531[k]
                   - f_5 * lsi1_531[k]
                   + f_3 * pc_y[k] * lsk_682[k];

        t_853[k] = f_16 * ksk_503[k]
                   + f_3 * pc_y[k] * lsk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pa_y, pc_y, pc_z, ksl0_630, ksk_467, \
                         ksk_468, ksk_504, ksl1_630, lsi0_531, lsi1_531, lsk_683, \
                         lsk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_17 * ksk_467[k]
                   + f_1 * lsi0_531[k]
                   - f_2 * lsi1_531[k]
                   + f_3 * pc_z[k] * lsk_683[k];

        t_855[k] = pa_y[k] * ksl0_630[k]
                   - f_14 * pc_y[k] * ksl1_630[k];

        t_856[k] = f_15 * ksk_504[k]
                   + f_3 * pc_y[k] * lsk_684[k];

        t_857[k] = f_18 * ksk_468[k]
                   + f_3 * pc_z[k] * lsk_684[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pa_y, pc_y, ksl0_633, ksl0_635, ksl0_636, \
                         ksk_505, ksk_506, ksk_507, ksl1_633, ksl1_635, ksl1_636, \
                         lsk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pa_y[k] * ksl0_633[k]
                   + f_16 * ksk_505[k]
                   - f_14 * pc_y[k] * ksl1_633[k];

        t_859[k] = f_15 * ksk_506[k]
                   + f_3 * pc_y[k] * lsk_686[k];

        t_860[k] = pa_y[k] * ksl0_635[k]
                   - f_14 * pc_y[k] * ksl1_635[k];

        t_861[k] = pa_y[k] * ksl0_636[k]
                   + f_17 * ksk_507[k]
                   - f_14 * pc_y[k] * ksl1_636[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pa_y, pc_y, pc_z, ksl0_639, ksl0_640, \
                         ksk_471, ksk_509, ksk_510, ksl1_639, ksl1_640, lsk_687, \
                         lsk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_18 * ksk_471[k]
                   + f_3 * pc_z[k] * lsk_687[k];

        t_863[k] = f_15 * ksk_509[k]
                   + f_3 * pc_y[k] * lsk_689[k];

        t_864[k] = pa_y[k] * ksl0_639[k]
                   - f_14 * pc_y[k] * ksl1_639[k];

        t_865[k] = pa_y[k] * ksl0_640[k]
                   + f_18 * ksk_510[k]
                   - f_14 * pc_y[k] * ksl1_640[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pa_y, pc_y, pc_z, ksl0_642, ksl0_644, \
                         ksk_474, ksk_512, ksk_513, ksl1_642, ksl1_644, lsk_690, \
                         lsk_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * ksk_474[k]
                   + f_3 * pc_z[k] * lsk_690[k];

        t_867[k] = pa_y[k] * ksl0_642[k]
                   + f_16 * ksk_512[k]
                   - f_14 * pc_y[k] * ksl1_642[k];

        t_868[k] = f_15 * ksk_513[k]
                   + f_3 * pc_y[k] * lsk_693[k];

        t_869[k] = pa_y[k] * ksl0_644[k]
                   - f_14 * pc_y[k] * ksl1_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pa_y, pc_y, pc_z, ksl0_645, ksl0_647, ksk_478, \
                         ksk_514, ksk_516, ksl1_645, ksl1_647, \
                         lsk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pa_y[k] * ksl0_645[k]
                   + f_19 * ksk_514[k]
                   - f_14 * pc_y[k] * ksl1_645[k];

        t_871[k] = f_18 * ksk_478[k]
                   + f_3 * pc_z[k] * lsk_694[k];

        t_872[k] = pa_y[k] * ksl0_647[k]
                   + f_17 * ksk_516[k]
                   - f_14 * pc_y[k] * ksl1_647[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pa_y, pc_y, ksl0_648, ksl0_650, ksl0_651, \
                         ksk_517, ksk_518, ksk_519, ksl1_648, ksl1_650, ksl1_651, \
                         lsk_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pa_y[k] * ksl0_648[k]
                   + f_16 * ksk_517[k]
                   - f_14 * pc_y[k] * ksl1_648[k];

        t_874[k] = f_15 * ksk_518[k]
                   + f_3 * pc_y[k] * lsk_698[k];

        t_875[k] = pa_y[k] * ksl0_650[k]
                   - f_14 * pc_y[k] * ksl1_650[k];

        t_876[k] = pa_y[k] * ksl0_651[k]
                   + f_20 * ksk_519[k]
                   - f_14 * pc_y[k] * ksl1_651[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pa_y, pc_y, pc_z, ksl0_653, ksl0_654, ksk_483, \
                         ksk_521, ksk_522, ksl1_653, ksl1_654, \
                         lsk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_18 * ksk_483[k]
                   + f_3 * pc_z[k] * lsk_699[k];

        t_878[k] = pa_y[k] * ksl0_653[k]
                   + f_18 * ksk_521[k]
                   - f_14 * pc_y[k] * ksl1_653[k];

        t_879[k] = pa_y[k] * ksl0_654[k]
                   + f_17 * ksk_522[k]
                   - f_14 * pc_y[k] * ksl1_654[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pa_y, pc_x, pc_y, ksl0_655, ksl0_657, \
                         ksk_523, ksk_524, ksk_712, ksl1_655, ksl1_657, lsk_704, \
                         lsk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = pa_y[k] * ksl0_655[k]
                   + f_16 * ksk_523[k]
                   - f_14 * pc_y[k] * ksl1_655[k];

        t_881[k] = f_15 * ksk_524[k]
                   + f_3 * pc_y[k] * lsk_704[k];

        t_882[k] = pa_y[k] * ksl0_657[k]
                   - f_14 * pc_y[k] * ksl1_657[k];

        t_883[k] = f_17 * ksk_712[k]
                   + f_3 * pc_x[k] * lsk_712[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, pc_x, ksk_713, ksk_714, ksk_715, \
                         ksk_716, ksk_717, lsk_713, lsk_714, lsk_715, lsk_716, \
                         lsk_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_17 * ksk_713[k]
                   + f_3 * pc_x[k] * lsk_713[k];

        t_885[k] = f_17 * ksk_714[k]
                   + f_3 * pc_x[k] * lsk_714[k];

        t_886[k] = f_17 * ksk_715[k]
                   + f_3 * pc_x[k] * lsk_715[k];

        t_887[k] = f_17 * ksk_716[k]
                   + f_3 * pc_x[k] * lsk_716[k];

        t_888[k] = f_17 * ksk_717[k]
                   + f_3 * pc_x[k] * lsk_717[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, t_892, pc_x, pc_y, pc_z, ksk_496, ksk_532, \
                         ksk_718, ksk_719, lsi0_553, lsi1_553, lsk_712, lsk_718, \
                         lsk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_17 * ksk_718[k]
                   + f_3 * pc_x[k] * lsk_718[k];

        t_890[k] = f_17 * ksk_719[k]
                   + f_3 * pc_x[k] * lsk_719[k];

        t_891[k] = f_15 * ksk_532[k]
                   + f_1 * lsi0_553[k]
                   - f_2 * lsi1_553[k]
                   + f_3 * pc_y[k] * lsk_712[k];

        t_892[k] = f_18 * ksk_496[k]
                   + f_3 * pc_z[k] * lsk_712[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, ksk_534, ksk_535, ksk_536, lsi0_555, \
                         lsi0_556, lsi0_557, lsi1_555, lsi1_556, lsi1_557, lsk_714, lsk_715, \
                         lsk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_15 * ksk_534[k]
                   + f_12 * lsi0_555[k]
                   - f_13 * lsi1_555[k]
                   + f_3 * pc_y[k] * lsk_714[k];

        t_894[k] = f_15 * ksk_535[k]
                   + f_10 * lsi0_556[k]
                   - f_11 * lsi1_556[k]
                   + f_3 * pc_y[k] * lsk_715[k];

        t_895[k] = f_15 * ksk_536[k]
                   + f_8 * lsi0_557[k]
                   - f_9 * lsi1_557[k]
                   + f_3 * pc_y[k] * lsk_716[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_y, ksk_537, ksk_538, ksk_539, lsi0_558, \
                         lsi0_559, lsi1_558, lsi1_559, lsk_717, lsk_718, \
                         lsk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_15 * ksk_537[k]
                   + f_6 * lsi0_558[k]
                   - f_7 * lsi1_558[k]
                   + f_3 * pc_y[k] * lsk_717[k];

        t_897[k] = f_15 * ksk_538[k]
                   + f_4 * lsi0_559[k]
                   - f_5 * lsi1_559[k]
                   + f_3 * pc_y[k] * lsk_718[k];

        t_898[k] = f_15 * ksk_539[k]
                   + f_3 * pc_y[k] * lsk_719[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pa_y, pc_x, pc_y, pc_z, ksl0_674, \
                         ksk_504, ksk_720, ksl1_674, lsi0_560, lsi1_560, \
                         lsk_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pa_y[k] * ksl0_674[k]
                   - f_14 * pc_y[k] * ksl1_674[k];

        t_900[k] = f_17 * ksk_720[k]
                   + f_1 * lsi0_560[k]
                   - f_2 * lsi1_560[k]
                   + f_3 * pc_x[k] * lsk_720[k];

        t_901[k] = f_3 * pc_y[k] * lsk_720[k];

        t_902[k] = f_19 * ksk_504[k]
                   + f_3 * pc_z[k] * lsk_720[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pc_x, pc_y, ksk_725, lsi0_560, lsi0_565, \
                         lsi1_560, lsi1_565, lsk_721, lsk_722, \
                         lsk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_4 * lsi0_560[k]
                   - f_5 * lsi1_560[k]
                   + f_3 * pc_y[k] * lsk_721[k];

        t_904[k] = f_3 * pc_y[k] * lsk_722[k];

        t_905[k] = f_17 * ksk_725[k]
                   + f_12 * lsi0_565[k]
                   - f_13 * lsi1_565[k]
                   + f_3 * pc_x[k] * lsk_725[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pc_y, lsi0_561, lsi0_562, lsi1_561, lsi1_562, \
                         lsk_723, lsk_724, lsk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_6 * lsi0_561[k]
                   - f_7 * lsi1_561[k]
                   + f_3 * pc_y[k] * lsk_723[k];

        t_907[k] = f_4 * lsi0_562[k]
                   - f_5 * lsi1_562[k]
                   + f_3 * pc_y[k] * lsk_724[k];

        t_908[k] = f_3 * pc_y[k] * lsk_725[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, pc_y, ksk_729, lsi0_563, lsi0_564, \
                         lsi0_569, lsi1_563, lsi1_564, lsi1_569, lsk_726, lsk_727, \
                         lsk_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_17 * ksk_729[k]
                   + f_10 * lsi0_569[k]
                   - f_11 * lsi1_569[k]
                   + f_3 * pc_x[k] * lsk_729[k];

        t_910[k] = f_8 * lsi0_563[k]
                   - f_9 * lsi1_563[k]
                   + f_3 * pc_y[k] * lsk_726[k];

        t_911[k] = f_6 * lsi0_564[k]
                   - f_7 * lsi1_564[k]
                   + f_3 * pc_y[k] * lsk_727[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_x, pc_y, ksk_734, lsi0_565, lsi0_574, \
                         lsi1_565, lsi1_574, lsk_728, lsk_729, \
                         lsk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_4 * lsi0_565[k]
                   - f_5 * lsi1_565[k]
                   + f_3 * pc_y[k] * lsk_728[k];

        t_913[k] = f_3 * pc_y[k] * lsk_729[k];

        t_914[k] = f_17 * ksk_734[k]
                   + f_8 * lsi0_574[k]
                   - f_9 * lsi1_574[k]
                   + f_3 * pc_x[k] * lsk_734[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_y, lsi0_566, lsi0_567, lsi0_568, lsi1_566, \
                         lsi1_567, lsi1_568, lsk_730, lsk_731, \
                         lsk_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_10 * lsi0_566[k]
                   - f_11 * lsi1_566[k]
                   + f_3 * pc_y[k] * lsk_730[k];

        t_916[k] = f_8 * lsi0_567[k]
                   - f_9 * lsi1_567[k]
                   + f_3 * pc_y[k] * lsk_731[k];

        t_917[k] = f_6 * lsi0_568[k]
                   - f_7 * lsi1_568[k]
                   + f_3 * pc_y[k] * lsk_732[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_x, pc_y, ksk_740, lsi0_569, lsi0_580, \
                         lsi1_569, lsi1_580, lsk_733, lsk_734, \
                         lsk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_4 * lsi0_569[k]
                   - f_5 * lsi1_569[k]
                   + f_3 * pc_y[k] * lsk_733[k];

        t_919[k] = f_3 * pc_y[k] * lsk_734[k];

        t_920[k] = f_17 * ksk_740[k]
                   + f_6 * lsi0_580[k]
                   - f_7 * lsi1_580[k]
                   + f_3 * pc_x[k] * lsk_740[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksl0,
                                                          const size_t ksk, const size_t ksl1,
                                                          const size_t lsi0, const size_t lsi1,
                                                          const size_t lsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_675 = buffer.data(ksl0 + 675);
    const auto *ksl0_678 = buffer.data(ksl0 + 678);
    const auto *ksl0_681 = buffer.data(ksl0 + 681);
    const auto *ksl0_685 = buffer.data(ksl0 + 685);
    const auto *ksl0_687 = buffer.data(ksl0 + 687);
    const auto *ksl0_690 = buffer.data(ksl0 + 690);
    const auto *ksl0_692 = buffer.data(ksl0 + 692);
    const auto *ksl0_693 = buffer.data(ksl0 + 693);
    const auto *ksl0_696 = buffer.data(ksl0 + 696);
    const auto *ksl0_698 = buffer.data(ksl0 + 698);
    const auto *ksl0_699 = buffer.data(ksl0 + 699);
    const auto *ksl0_700 = buffer.data(ksl0 + 700);
    const auto *ksl0_711 = buffer.data(ksl0 + 711);

    const auto *ksk_539 = buffer.data(ksk + 539);
    const auto *ksk_540 = buffer.data(ksk + 540);
    const auto *ksk_543 = buffer.data(ksk + 543);
    const auto *ksk_545 = buffer.data(ksk + 545);
    const auto *ksk_546 = buffer.data(ksk + 546);
    const auto *ksk_547 = buffer.data(ksk + 547);
    const auto *ksk_549 = buffer.data(ksk + 549);
    const auto *ksk_550 = buffer.data(ksk + 550);
    const auto *ksk_551 = buffer.data(ksk + 551);
    const auto *ksk_552 = buffer.data(ksk + 552);
    const auto *ksk_554 = buffer.data(ksk + 554);
    const auto *ksk_555 = buffer.data(ksk + 555);
    const auto *ksk_556 = buffer.data(ksk + 556);
    const auto *ksk_557 = buffer.data(ksk + 557);
    const auto *ksk_558 = buffer.data(ksk + 558);
    const auto *ksk_560 = buffer.data(ksk + 560);
    const auto *ksk_568 = buffer.data(ksk + 568);
    const auto *ksk_575 = buffer.data(ksk + 575);
    const auto *ksk_576 = buffer.data(ksk + 576);
    const auto *ksk_578 = buffer.data(ksk + 578);
    const auto *ksk_581 = buffer.data(ksk + 581);
    const auto *ksk_585 = buffer.data(ksk + 585);
    const auto *ksk_590 = buffer.data(ksk + 590);
    const auto *ksk_596 = buffer.data(ksk + 596);
    const auto *ksk_606 = buffer.data(ksk + 606);
    const auto *ksk_607 = buffer.data(ksk + 607);
    const auto *ksk_608 = buffer.data(ksk + 608);
    const auto *ksk_609 = buffer.data(ksk + 609);
    const auto *ksk_610 = buffer.data(ksk + 610);
    const auto *ksk_611 = buffer.data(ksk + 611);
    const auto *ksk_612 = buffer.data(ksk + 612);
    const auto *ksk_747 = buffer.data(ksk + 747);
    const auto *ksk_748 = buffer.data(ksk + 748);
    const auto *ksk_749 = buffer.data(ksk + 749);
    const auto *ksk_750 = buffer.data(ksk + 750);
    const auto *ksk_751 = buffer.data(ksk + 751);
    const auto *ksk_752 = buffer.data(ksk + 752);
    const auto *ksk_753 = buffer.data(ksk + 753);
    const auto *ksk_755 = buffer.data(ksk + 755);
    const auto *ksk_756 = buffer.data(ksk + 756);
    const auto *ksk_759 = buffer.data(ksk + 759);
    const auto *ksk_762 = buffer.data(ksk + 762);
    const auto *ksk_766 = buffer.data(ksk + 766);
    const auto *ksk_771 = buffer.data(ksk + 771);
    const auto *ksk_777 = buffer.data(ksk + 777);
    const auto *ksk_784 = buffer.data(ksk + 784);
    const auto *ksk_786 = buffer.data(ksk + 786);
    const auto *ksk_787 = buffer.data(ksk + 787);
    const auto *ksk_788 = buffer.data(ksk + 788);
    const auto *ksk_789 = buffer.data(ksk + 789);
    const auto *ksk_790 = buffer.data(ksk + 790);
    const auto *ksk_791 = buffer.data(ksk + 791);
    const auto *ksk_797 = buffer.data(ksk + 797);
    const auto *ksk_801 = buffer.data(ksk + 801);
    const auto *ksk_806 = buffer.data(ksk + 806);
    const auto *ksk_812 = buffer.data(ksk + 812);
    const auto *ksk_819 = buffer.data(ksk + 819);
    const auto *ksk_820 = buffer.data(ksk + 820);
    const auto *ksk_821 = buffer.data(ksk + 821);
    const auto *ksk_822 = buffer.data(ksk + 822);
    const auto *ksk_823 = buffer.data(ksk + 823);
    const auto *ksk_824 = buffer.data(ksk + 824);
    const auto *ksk_825 = buffer.data(ksk + 825);
    const auto *ksk_826 = buffer.data(ksk + 826);
    const auto *ksk_827 = buffer.data(ksk + 827);
    const auto *ksk_828 = buffer.data(ksk + 828);

    const auto *ksl1_675 = buffer.data(ksl1 + 675);
    const auto *ksl1_678 = buffer.data(ksl1 + 678);
    const auto *ksl1_681 = buffer.data(ksl1 + 681);
    const auto *ksl1_685 = buffer.data(ksl1 + 685);
    const auto *ksl1_687 = buffer.data(ksl1 + 687);
    const auto *ksl1_690 = buffer.data(ksl1 + 690);
    const auto *ksl1_692 = buffer.data(ksl1 + 692);
    const auto *ksl1_693 = buffer.data(ksl1 + 693);
    const auto *ksl1_696 = buffer.data(ksl1 + 696);
    const auto *ksl1_698 = buffer.data(ksl1 + 698);
    const auto *ksl1_699 = buffer.data(ksl1 + 699);
    const auto *ksl1_700 = buffer.data(ksl1 + 700);
    const auto *ksl1_711 = buffer.data(ksl1 + 711);

    const auto *lsi0_570 = buffer.data(lsi0 + 570);
    const auto *lsi0_571 = buffer.data(lsi0 + 571);
    const auto *lsi0_572 = buffer.data(lsi0 + 572);
    const auto *lsi0_573 = buffer.data(lsi0 + 573);
    const auto *lsi0_574 = buffer.data(lsi0 + 574);
    const auto *lsi0_581 = buffer.data(lsi0 + 581);
    const auto *lsi0_582 = buffer.data(lsi0 + 582);
    const auto *lsi0_583 = buffer.data(lsi0 + 583);
    const auto *lsi0_584 = buffer.data(lsi0 + 584);
    const auto *lsi0_585 = buffer.data(lsi0 + 585);
    const auto *lsi0_586 = buffer.data(lsi0 + 586);
    const auto *lsi0_587 = buffer.data(lsi0 + 587);
    const auto *lsi0_588 = buffer.data(lsi0 + 588);
    const auto *lsi0_590 = buffer.data(lsi0 + 590);
    const auto *lsi0_591 = buffer.data(lsi0 + 591);
    const auto *lsi0_593 = buffer.data(lsi0 + 593);
    const auto *lsi0_594 = buffer.data(lsi0 + 594);
    const auto *lsi0_595 = buffer.data(lsi0 + 595);
    const auto *lsi0_597 = buffer.data(lsi0 + 597);
    const auto *lsi0_598 = buffer.data(lsi0 + 598);
    const auto *lsi0_599 = buffer.data(lsi0 + 599);
    const auto *lsi0_600 = buffer.data(lsi0 + 600);
    const auto *lsi0_602 = buffer.data(lsi0 + 602);
    const auto *lsi0_603 = buffer.data(lsi0 + 603);
    const auto *lsi0_609 = buffer.data(lsi0 + 609);
    const auto *lsi0_610 = buffer.data(lsi0 + 610);
    const auto *lsi0_611 = buffer.data(lsi0 + 611);
    const auto *lsi0_612 = buffer.data(lsi0 + 612);
    const auto *lsi0_613 = buffer.data(lsi0 + 613);
    const auto *lsi0_615 = buffer.data(lsi0 + 615);
    const auto *lsi0_621 = buffer.data(lsi0 + 621);
    const auto *lsi0_625 = buffer.data(lsi0 + 625);
    const auto *lsi0_630 = buffer.data(lsi0 + 630);
    const auto *lsi0_636 = buffer.data(lsi0 + 636);
    const auto *lsi0_639 = buffer.data(lsi0 + 639);
    const auto *lsi0_640 = buffer.data(lsi0 + 640);
    const auto *lsi0_641 = buffer.data(lsi0 + 641);
    const auto *lsi0_642 = buffer.data(lsi0 + 642);
    const auto *lsi0_643 = buffer.data(lsi0 + 643);
    const auto *lsi0_644 = buffer.data(lsi0 + 644);

    const auto *lsi1_570 = buffer.data(lsi1 + 570);
    const auto *lsi1_571 = buffer.data(lsi1 + 571);
    const auto *lsi1_572 = buffer.data(lsi1 + 572);
    const auto *lsi1_573 = buffer.data(lsi1 + 573);
    const auto *lsi1_574 = buffer.data(lsi1 + 574);
    const auto *lsi1_581 = buffer.data(lsi1 + 581);
    const auto *lsi1_582 = buffer.data(lsi1 + 582);
    const auto *lsi1_583 = buffer.data(lsi1 + 583);
    const auto *lsi1_584 = buffer.data(lsi1 + 584);
    const auto *lsi1_585 = buffer.data(lsi1 + 585);
    const auto *lsi1_586 = buffer.data(lsi1 + 586);
    const auto *lsi1_587 = buffer.data(lsi1 + 587);
    const auto *lsi1_588 = buffer.data(lsi1 + 588);
    const auto *lsi1_590 = buffer.data(lsi1 + 590);
    const auto *lsi1_591 = buffer.data(lsi1 + 591);
    const auto *lsi1_593 = buffer.data(lsi1 + 593);
    const auto *lsi1_594 = buffer.data(lsi1 + 594);
    const auto *lsi1_595 = buffer.data(lsi1 + 595);
    const auto *lsi1_597 = buffer.data(lsi1 + 597);
    const auto *lsi1_598 = buffer.data(lsi1 + 598);
    const auto *lsi1_599 = buffer.data(lsi1 + 599);
    const auto *lsi1_600 = buffer.data(lsi1 + 600);
    const auto *lsi1_602 = buffer.data(lsi1 + 602);
    const auto *lsi1_603 = buffer.data(lsi1 + 603);
    const auto *lsi1_609 = buffer.data(lsi1 + 609);
    const auto *lsi1_610 = buffer.data(lsi1 + 610);
    const auto *lsi1_611 = buffer.data(lsi1 + 611);
    const auto *lsi1_612 = buffer.data(lsi1 + 612);
    const auto *lsi1_613 = buffer.data(lsi1 + 613);
    const auto *lsi1_615 = buffer.data(lsi1 + 615);
    const auto *lsi1_621 = buffer.data(lsi1 + 621);
    const auto *lsi1_625 = buffer.data(lsi1 + 625);
    const auto *lsi1_630 = buffer.data(lsi1 + 630);
    const auto *lsi1_636 = buffer.data(lsi1 + 636);
    const auto *lsi1_639 = buffer.data(lsi1 + 639);
    const auto *lsi1_640 = buffer.data(lsi1 + 640);
    const auto *lsi1_641 = buffer.data(lsi1 + 641);
    const auto *lsi1_642 = buffer.data(lsi1 + 642);
    const auto *lsi1_643 = buffer.data(lsi1 + 643);
    const auto *lsi1_644 = buffer.data(lsi1 + 644);

    const auto *lsk_735 = buffer.data(lsk + 735);
    const auto *lsk_736 = buffer.data(lsk + 736);
    const auto *lsk_737 = buffer.data(lsk + 737);
    const auto *lsk_738 = buffer.data(lsk + 738);
    const auto *lsk_739 = buffer.data(lsk + 739);
    const auto *lsk_740 = buffer.data(lsk + 740);
    const auto *lsk_747 = buffer.data(lsk + 747);
    const auto *lsk_748 = buffer.data(lsk + 748);
    const auto *lsk_749 = buffer.data(lsk + 749);
    const auto *lsk_750 = buffer.data(lsk + 750);
    const auto *lsk_751 = buffer.data(lsk + 751);
    const auto *lsk_752 = buffer.data(lsk + 752);
    const auto *lsk_753 = buffer.data(lsk + 753);
    const auto *lsk_754 = buffer.data(lsk + 754);
    const auto *lsk_755 = buffer.data(lsk + 755);
    const auto *lsk_756 = buffer.data(lsk + 756);
    const auto *lsk_757 = buffer.data(lsk + 757);
    const auto *lsk_758 = buffer.data(lsk + 758);
    const auto *lsk_759 = buffer.data(lsk + 759);
    const auto *lsk_761 = buffer.data(lsk + 761);
    const auto *lsk_762 = buffer.data(lsk + 762);
    const auto *lsk_763 = buffer.data(lsk + 763);
    const auto *lsk_765 = buffer.data(lsk + 765);
    const auto *lsk_766 = buffer.data(lsk + 766);
    const auto *lsk_767 = buffer.data(lsk + 767);
    const auto *lsk_768 = buffer.data(lsk + 768);
    const auto *lsk_770 = buffer.data(lsk + 770);
    const auto *lsk_771 = buffer.data(lsk + 771);
    const auto *lsk_772 = buffer.data(lsk + 772);
    const auto *lsk_773 = buffer.data(lsk + 773);
    const auto *lsk_774 = buffer.data(lsk + 774);
    const auto *lsk_776 = buffer.data(lsk + 776);
    const auto *lsk_777 = buffer.data(lsk + 777);
    const auto *lsk_784 = buffer.data(lsk + 784);
    const auto *lsk_785 = buffer.data(lsk + 785);
    const auto *lsk_786 = buffer.data(lsk + 786);
    const auto *lsk_787 = buffer.data(lsk + 787);
    const auto *lsk_788 = buffer.data(lsk + 788);
    const auto *lsk_789 = buffer.data(lsk + 789);
    const auto *lsk_790 = buffer.data(lsk + 790);
    const auto *lsk_791 = buffer.data(lsk + 791);
    const auto *lsk_792 = buffer.data(lsk + 792);
    const auto *lsk_794 = buffer.data(lsk + 794);
    const auto *lsk_795 = buffer.data(lsk + 795);
    const auto *lsk_797 = buffer.data(lsk + 797);
    const auto *lsk_798 = buffer.data(lsk + 798);
    const auto *lsk_801 = buffer.data(lsk + 801);
    const auto *lsk_802 = buffer.data(lsk + 802);
    const auto *lsk_806 = buffer.data(lsk + 806);
    const auto *lsk_807 = buffer.data(lsk + 807);
    const auto *lsk_812 = buffer.data(lsk + 812);
    const auto *lsk_819 = buffer.data(lsk + 819);
    const auto *lsk_820 = buffer.data(lsk + 820);
    const auto *lsk_821 = buffer.data(lsk + 821);
    const auto *lsk_822 = buffer.data(lsk + 822);
    const auto *lsk_823 = buffer.data(lsk + 823);
    const auto *lsk_824 = buffer.data(lsk + 824);
    const auto *lsk_825 = buffer.data(lsk + 825);
    const auto *lsk_826 = buffer.data(lsk + 826);
    const auto *lsk_827 = buffer.data(lsk + 827);
    const auto *lsk_828 = buffer.data(lsk + 828);

#pragma omp simd aligned(t_921, t_922, t_923, pc_y, lsi0_570, lsi0_571, lsi0_572, lsi1_570, \
                         lsi1_571, lsi1_572, lsk_735, lsk_736, \
                         lsk_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_12 * lsi0_570[k]
                   - f_13 * lsi1_570[k]
                   + f_3 * pc_y[k] * lsk_735[k];

        t_922[k] = f_10 * lsi0_571[k]
                   - f_11 * lsi1_571[k]
                   + f_3 * pc_y[k] * lsk_736[k];

        t_923[k] = f_8 * lsi0_572[k]
                   - f_9 * lsi1_572[k]
                   + f_3 * pc_y[k] * lsk_737[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_y, lsi0_573, lsi0_574, lsi1_573, lsi1_574, \
                         lsk_738, lsk_739, lsk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_6 * lsi0_573[k]
                   - f_7 * lsi1_573[k]
                   + f_3 * pc_y[k] * lsk_738[k];

        t_925[k] = f_4 * lsi0_574[k]
                   - f_5 * lsi1_574[k]
                   + f_3 * pc_y[k] * lsk_739[k];

        t_926[k] = f_3 * pc_y[k] * lsk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pc_x, ksk_747, ksk_748, ksk_749, ksk_750, \
                         lsi0_587, lsi1_587, lsk_747, lsk_748, lsk_749, \
                         lsk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_17 * ksk_747[k]
                   + f_4 * lsi0_587[k]
                   - f_5 * lsi1_587[k]
                   + f_3 * pc_x[k] * lsk_747[k];

        t_928[k] = f_17 * ksk_748[k]
                   + f_3 * pc_x[k] * lsk_748[k];

        t_929[k] = f_17 * ksk_749[k]
                   + f_3 * pc_x[k] * lsk_749[k];

        t_930[k] = f_17 * ksk_750[k]
                   + f_3 * pc_x[k] * lsk_750[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pc_x, pc_y, ksk_751, ksk_752, \
                         ksk_753, ksk_755, lsk_747, lsk_751, lsk_752, lsk_753, \
                         lsk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_17 * ksk_751[k]
                   + f_3 * pc_x[k] * lsk_751[k];

        t_932[k] = f_17 * ksk_752[k]
                   + f_3 * pc_x[k] * lsk_752[k];

        t_933[k] = f_17 * ksk_753[k]
                   + f_3 * pc_x[k] * lsk_753[k];

        t_934[k] = f_3 * pc_y[k] * lsk_747[k];

        t_935[k] = f_17 * ksk_755[k]
                   + f_3 * pc_x[k] * lsk_755[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, pc_y, lsi0_581, lsi0_582, lsi0_583, lsi1_581, \
                         lsi1_582, lsi1_583, lsk_748, lsk_749, \
                         lsk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_1 * lsi0_581[k]
                   - f_2 * lsi1_581[k]
                   + f_3 * pc_y[k] * lsk_748[k];

        t_937[k] = f_22 * lsi0_582[k]
                   - f_23 * lsi1_582[k]
                   + f_3 * pc_y[k] * lsk_749[k];

        t_938[k] = f_12 * lsi0_583[k]
                   - f_13 * lsi1_583[k]
                   + f_3 * pc_y[k] * lsk_750[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, pc_y, lsi0_584, lsi0_585, lsi0_586, lsi1_584, \
                         lsi1_585, lsi1_586, lsk_751, lsk_752, \
                         lsk_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_10 * lsi0_584[k]
                   - f_11 * lsi1_584[k]
                   + f_3 * pc_y[k] * lsk_751[k];

        t_940[k] = f_8 * lsi0_585[k]
                   - f_9 * lsi1_585[k]
                   + f_3 * pc_y[k] * lsk_752[k];

        t_941[k] = f_6 * lsi0_586[k]
                   - f_7 * lsi1_586[k]
                   + f_3 * pc_y[k] * lsk_753[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, pc_z, ksk_539, ksk_756, \
                         lsi0_587, lsi0_588, lsi1_587, lsi1_588, lsk_754, lsk_755, \
                         lsk_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_4 * lsi0_587[k]
                   - f_5 * lsi1_587[k]
                   + f_3 * pc_y[k] * lsk_754[k];

        t_943[k] = f_3 * pc_y[k] * lsk_755[k];

        t_944[k] = f_19 * ksk_539[k]
                   + f_1 * lsi0_587[k]
                   - f_2 * lsi1_587[k]
                   + f_3 * pc_z[k] * lsk_755[k];

        t_945[k] = f_16 * ksk_756[k]
                   + f_1 * lsi0_588[k]
                   - f_2 * lsi1_588[k]
                   + f_3 * pc_x[k] * lsk_756[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, t_949, pc_x, pc_y, pc_z, ksk_540, ksk_759, \
                         lsi0_591, lsi1_591, lsk_756, lsk_757, \
                         lsk_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_20 * ksk_540[k]
                   + f_3 * pc_y[k] * lsk_756[k];

        t_947[k] = f_3 * pc_z[k] * lsk_756[k];

        t_948[k] = f_16 * ksk_759[k]
                   + f_12 * lsi0_591[k]
                   - f_13 * lsi1_591[k]
                   + f_3 * pc_x[k] * lsk_759[k];

        t_949[k] = f_3 * pc_z[k] * lsk_757[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pc_x, pc_z, ksk_762, lsi0_588, lsi0_594, \
                         lsi1_588, lsi1_594, lsk_758, lsk_759, \
                         lsk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_4 * lsi0_588[k]
                   - f_5 * lsi1_588[k]
                   + f_3 * pc_z[k] * lsk_758[k];

        t_951[k] = f_16 * ksk_762[k]
                   + f_10 * lsi0_594[k]
                   - f_11 * lsi1_594[k]
                   + f_3 * pc_x[k] * lsk_762[k];

        t_952[k] = f_3 * pc_z[k] * lsk_759[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pc_x, pc_y, pc_z, ksk_545, ksk_766, \
                         lsi0_590, lsi0_598, lsi1_590, lsi1_598, lsk_761, lsk_762, \
                         lsk_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_20 * ksk_545[k]
                   + f_3 * pc_y[k] * lsk_761[k];

        t_954[k] = f_6 * lsi0_590[k]
                   - f_7 * lsi1_590[k]
                   + f_3 * pc_z[k] * lsk_761[k];

        t_955[k] = f_16 * ksk_766[k]
                   + f_8 * lsi0_598[k]
                   - f_9 * lsi1_598[k]
                   + f_3 * pc_x[k] * lsk_766[k];

        t_956[k] = f_3 * pc_z[k] * lsk_762[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, pc_y, pc_z, ksk_549, lsi0_591, lsi0_593, \
                         lsi1_591, lsi1_593, lsk_763, lsk_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_4 * lsi0_591[k]
                   - f_5 * lsi1_591[k]
                   + f_3 * pc_z[k] * lsk_763[k];

        t_958[k] = f_20 * ksk_549[k]
                   + f_3 * pc_y[k] * lsk_765[k];

        t_959[k] = f_8 * lsi0_593[k]
                   - f_9 * lsi1_593[k]
                   + f_3 * pc_z[k] * lsk_765[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, pc_x, pc_z, ksk_771, lsi0_594, lsi0_603, \
                         lsi1_594, lsi1_603, lsk_766, lsk_767, \
                         lsk_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_16 * ksk_771[k]
                   + f_6 * lsi0_603[k]
                   - f_7 * lsi1_603[k]
                   + f_3 * pc_x[k] * lsk_771[k];

        t_961[k] = f_3 * pc_z[k] * lsk_766[k];

        t_962[k] = f_4 * lsi0_594[k]
                   - f_5 * lsi1_594[k]
                   + f_3 * pc_z[k] * lsk_767[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, pc_y, pc_z, ksk_554, lsi0_595, lsi0_597, \
                         lsi1_595, lsi1_597, lsk_768, lsk_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = f_6 * lsi0_595[k]
                   - f_7 * lsi1_595[k]
                   + f_3 * pc_z[k] * lsk_768[k];

        t_964[k] = f_20 * ksk_554[k]
                   + f_3 * pc_y[k] * lsk_770[k];

        t_965[k] = f_10 * lsi0_597[k]
                   - f_11 * lsi1_597[k]
                   + f_3 * pc_z[k] * lsk_770[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_x, pc_z, ksk_777, lsi0_598, lsi0_609, \
                         lsi1_598, lsi1_609, lsk_771, lsk_772, \
                         lsk_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_16 * ksk_777[k]
                   + f_4 * lsi0_609[k]
                   - f_5 * lsi1_609[k]
                   + f_3 * pc_x[k] * lsk_777[k];

        t_967[k] = f_3 * pc_z[k] * lsk_771[k];

        t_968[k] = f_4 * lsi0_598[k]
                   - f_5 * lsi1_598[k]
                   + f_3 * pc_z[k] * lsk_772[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, pc_y, pc_z, ksk_560, lsi0_599, lsi0_600, \
                         lsi0_602, lsi1_599, lsi1_600, lsi1_602, lsk_773, lsk_774, \
                         lsk_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_6 * lsi0_599[k]
                   - f_7 * lsi1_599[k]
                   + f_3 * pc_z[k] * lsk_773[k];

        t_970[k] = f_8 * lsi0_600[k]
                   - f_9 * lsi1_600[k]
                   + f_3 * pc_z[k] * lsk_774[k];

        t_971[k] = f_20 * ksk_560[k]
                   + f_3 * pc_y[k] * lsk_776[k];

        t_972[k] = f_12 * lsi0_602[k]
                   - f_13 * lsi1_602[k]
                   + f_3 * pc_z[k] * lsk_776[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, t_977, pc_x, pc_z, ksk_784, ksk_786, \
                         ksk_787, ksk_788, lsk_777, lsk_784, lsk_786, lsk_787, \
                         lsk_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_16 * ksk_784[k]
                   + f_3 * pc_x[k] * lsk_784[k];

        t_974[k] = f_3 * pc_z[k] * lsk_777[k];

        t_975[k] = f_16 * ksk_786[k]
                   + f_3 * pc_x[k] * lsk_786[k];

        t_976[k] = f_16 * ksk_787[k]
                   + f_3 * pc_x[k] * lsk_787[k];

        t_977[k] = f_16 * ksk_788[k]
                   + f_3 * pc_x[k] * lsk_788[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, pc_x, pc_y, ksk_568, ksk_789, ksk_790, \
                         ksk_791, lsi0_609, lsi1_609, lsk_784, lsk_789, lsk_790, \
                         lsk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_16 * ksk_789[k]
                   + f_3 * pc_x[k] * lsk_789[k];

        t_979[k] = f_16 * ksk_790[k]
                   + f_3 * pc_x[k] * lsk_790[k];

        t_980[k] = f_16 * ksk_791[k]
                   + f_3 * pc_x[k] * lsk_791[k];

        t_981[k] = f_20 * ksk_568[k]
                   + f_1 * lsi0_609[k]
                   - f_2 * lsi1_609[k]
                   + f_3 * pc_y[k] * lsk_784[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, t_985, pc_z, lsi0_609, lsi0_610, lsi0_611, \
                         lsi1_609, lsi1_610, lsi1_611, lsk_784, lsk_785, lsk_786, \
                         lsk_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_3 * pc_z[k] * lsk_784[k];

        t_983[k] = f_4 * lsi0_609[k]
                   - f_5 * lsi1_609[k]
                   + f_3 * pc_z[k] * lsk_785[k];

        t_984[k] = f_6 * lsi0_610[k]
                   - f_7 * lsi1_610[k]
                   + f_3 * pc_z[k] * lsk_786[k];

        t_985[k] = f_8 * lsi0_611[k]
                   - f_9 * lsi1_611[k]
                   + f_3 * pc_z[k] * lsk_787[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, t_989, pc_y, pc_z, ksk_575, lsi0_612, lsi0_613, \
                         lsi0_615, lsi1_612, lsi1_613, lsi1_615, lsk_788, lsk_789, \
                         lsk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_10 * lsi0_612[k]
                   - f_11 * lsi1_612[k]
                   + f_3 * pc_z[k] * lsk_788[k];

        t_987[k] = f_12 * lsi0_613[k]
                   - f_13 * lsi1_613[k]
                   + f_3 * pc_z[k] * lsk_789[k];

        t_988[k] = f_20 * ksk_575[k]
                   + f_3 * pc_y[k] * lsk_791[k];

        t_989[k] = f_1 * lsi0_615[k]
                   - f_2 * lsi1_615[k]
                   + f_3 * pc_z[k] * lsk_791[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, pa_z, pc_y, pc_z, ksl0_675, ksl0_678, \
                         ksk_540, ksk_576, ksl1_675, ksl1_678, \
                         lsk_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = pa_z[k] * ksl0_675[k]
                   - f_14 * pc_z[k] * ksl1_675[k];

        t_991[k] = f_19 * ksk_576[k]
                   + f_3 * pc_y[k] * lsk_792[k];

        t_992[k] = f_15 * ksk_540[k]
                   + f_3 * pc_z[k] * lsk_792[k];

        t_993[k] = pa_z[k] * ksl0_678[k]
                   - f_14 * pc_z[k] * ksl1_678[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, pa_z, pc_x, pc_y, pc_z, ksl0_681, ksk_578, \
                         ksk_797, ksl1_681, lsi0_621, lsi1_621, lsk_794, \
                         lsk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_19 * ksk_578[k]
                   + f_3 * pc_y[k] * lsk_794[k];

        t_995[k] = f_16 * ksk_797[k]
                   + f_12 * lsi0_621[k]
                   - f_13 * lsi1_621[k]
                   + f_3 * pc_x[k] * lsk_797[k];

        t_996[k] = pa_z[k] * ksl0_681[k]
                   - f_14 * pc_z[k] * ksl1_681[k];
    }

#pragma omp simd aligned(t_997, t_998, t_999, pc_x, pc_y, pc_z, ksk_543, ksk_581, ksk_801, \
                         lsi0_625, lsi1_625, lsk_795, lsk_797, \
                         lsk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = f_15 * ksk_543[k]
                   + f_3 * pc_z[k] * lsk_795[k];

        t_998[k] = f_19 * ksk_581[k]
                   + f_3 * pc_y[k] * lsk_797[k];

        t_999[k] = f_16 * ksk_801[k]
                   + f_10 * lsi0_625[k]
                   - f_11 * lsi1_625[k]
                   + f_3 * pc_x[k] * lsk_801[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pa_z, pc_y, pc_z, ksl0_685, ksl0_687, \
                         ksk_546, ksk_547, ksk_585, ksl1_685, ksl1_687, lsk_798, \
                         lsk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = pa_z[k] * ksl0_685[k]
                    - f_14 * pc_z[k] * ksl1_685[k];

        t_1001[k] = f_15 * ksk_546[k]
                    + f_3 * pc_z[k] * lsk_798[k];

        t_1002[k] = pa_z[k] * ksl0_687[k]
                    + f_16 * ksk_547[k]
                    - f_14 * pc_z[k] * ksl1_687[k];

        t_1003[k] = f_19 * ksk_585[k]
                    + f_3 * pc_y[k] * lsk_801[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, pa_z, pc_x, pc_z, ksl0_690, ksk_550, ksk_806, \
                         ksl1_690, lsi0_630, lsi1_630, lsk_802, \
                         lsk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_16 * ksk_806[k]
                    + f_8 * lsi0_630[k]
                    - f_9 * lsi1_630[k]
                    + f_3 * pc_x[k] * lsk_806[k];

        t_1005[k] = pa_z[k] * ksl0_690[k]
                    - f_14 * pc_z[k] * ksl1_690[k];

        t_1006[k] = f_15 * ksk_550[k]
                    + f_3 * pc_z[k] * lsk_802[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, pa_z, pc_y, pc_z, ksl0_692, ksl0_693, \
                         ksk_551, ksk_552, ksk_590, ksl1_692, ksl1_693, \
                         lsk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = pa_z[k] * ksl0_692[k]
                    + f_16 * ksk_551[k]
                    - f_14 * pc_z[k] * ksl1_692[k];

        t_1008[k] = pa_z[k] * ksl0_693[k]
                    + f_17 * ksk_552[k]
                    - f_14 * pc_z[k] * ksl1_693[k];

        t_1009[k] = f_19 * ksk_590[k]
                    + f_3 * pc_y[k] * lsk_806[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, pa_z, pc_x, pc_z, ksl0_696, ksk_555, ksk_812, \
                         ksl1_696, lsi0_636, lsi1_636, lsk_807, \
                         lsk_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_16 * ksk_812[k]
                    + f_6 * lsi0_636[k]
                    - f_7 * lsi1_636[k]
                    + f_3 * pc_x[k] * lsk_812[k];

        t_1011[k] = pa_z[k] * ksl0_696[k]
                    - f_14 * pc_z[k] * ksl1_696[k];

        t_1012[k] = f_15 * ksk_555[k]
                    + f_3 * pc_z[k] * lsk_807[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pa_z, pc_z, ksl0_698, ksl0_699, ksl0_700, \
                         ksk_556, ksk_557, ksk_558, ksl1_698, ksl1_699, \
                         ksl1_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = pa_z[k] * ksl0_698[k]
                    + f_16 * ksk_556[k]
                    - f_14 * pc_z[k] * ksl1_698[k];

        t_1014[k] = pa_z[k] * ksl0_699[k]
                    + f_17 * ksk_557[k]
                    - f_14 * pc_z[k] * ksl1_699[k];

        t_1015[k] = pa_z[k] * ksl0_700[k]
                    + f_18 * ksk_558[k]
                    - f_14 * pc_z[k] * ksl1_700[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, ksk_596, ksk_819, \
                         ksk_820, ksk_821, lsi0_643, lsi1_643, lsk_812, lsk_819, lsk_820, \
                         lsk_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_19 * ksk_596[k]
                    + f_3 * pc_y[k] * lsk_812[k];

        t_1017[k] = f_16 * ksk_819[k]
                    + f_4 * lsi0_643[k]
                    - f_5 * lsi1_643[k]
                    + f_3 * pc_x[k] * lsk_819[k];

        t_1018[k] = f_16 * ksk_820[k]
                    + f_3 * pc_x[k] * lsk_820[k];

        t_1019[k] = f_16 * ksk_821[k]
                    + f_3 * pc_x[k] * lsk_821[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, t_1024, pc_x, ksk_822, ksk_823, \
                         ksk_824, ksk_825, ksk_826, lsk_822, lsk_823, lsk_824, lsk_825, \
                         lsk_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_16 * ksk_822[k]
                    + f_3 * pc_x[k] * lsk_822[k];

        t_1021[k] = f_16 * ksk_823[k]
                    + f_3 * pc_x[k] * lsk_823[k];

        t_1022[k] = f_16 * ksk_824[k]
                    + f_3 * pc_x[k] * lsk_824[k];

        t_1023[k] = f_16 * ksk_825[k]
                    + f_3 * pc_x[k] * lsk_825[k];

        t_1024[k] = f_16 * ksk_826[k]
                    + f_3 * pc_x[k] * lsk_826[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, pa_z, pc_x, pc_z, ksl0_711, ksk_568, ksk_827, \
                         ksl1_711, lsk_820, lsk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_16 * ksk_827[k]
                    + f_3 * pc_x[k] * lsk_827[k];

        t_1026[k] = pa_z[k] * ksl0_711[k]
                    - f_14 * pc_z[k] * ksl1_711[k];

        t_1027[k] = f_15 * ksk_568[k]
                    + f_3 * pc_z[k] * lsk_820[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, pc_y, ksk_606, ksk_607, ksk_608, lsi0_639, \
                         lsi0_640, lsi0_641, lsi1_639, lsi1_640, lsi1_641, lsk_822, lsk_823, \
                         lsk_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_19 * ksk_606[k]
                    + f_12 * lsi0_639[k]
                    - f_13 * lsi1_639[k]
                    + f_3 * pc_y[k] * lsk_822[k];

        t_1029[k] = f_19 * ksk_607[k]
                    + f_10 * lsi0_640[k]
                    - f_11 * lsi1_640[k]
                    + f_3 * pc_y[k] * lsk_823[k];

        t_1030[k] = f_19 * ksk_608[k]
                    + f_8 * lsi0_641[k]
                    - f_9 * lsi1_641[k]
                    + f_3 * pc_y[k] * lsk_824[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, pc_y, ksk_609, ksk_610, ksk_611, lsi0_642, \
                         lsi0_643, lsi1_642, lsi1_643, lsk_825, lsk_826, \
                         lsk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_19 * ksk_609[k]
                    + f_6 * lsi0_642[k]
                    - f_7 * lsi1_642[k]
                    + f_3 * pc_y[k] * lsk_825[k];

        t_1032[k] = f_19 * ksk_610[k]
                    + f_4 * lsi0_643[k]
                    - f_5 * lsi1_643[k]
                    + f_3 * pc_y[k] * lsk_826[k];

        t_1033[k] = f_19 * ksk_611[k]
                    + f_3 * pc_y[k] * lsk_827[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, pc_x, pc_y, pc_z, ksk_575, ksk_612, ksk_828, \
                         lsi0_643, lsi0_644, lsi1_643, lsi1_644, lsk_827, \
                         lsk_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_15 * ksk_575[k]
                    + f_1 * lsi0_643[k]
                    - f_2 * lsi1_643[k]
                    + f_3 * pc_z[k] * lsk_827[k];

        t_1035[k] = f_16 * ksk_828[k]
                    + f_1 * lsi0_644[k]
                    - f_2 * lsi1_644[k]
                    + f_3 * pc_x[k] * lsk_828[k];

        t_1036[k] = f_18 * ksk_612[k]
                    + f_3 * pc_y[k] * lsk_828[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t ksk, const size_t lsi0,
                                                          const size_t lsi1, const size_t lsk,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk_576 = buffer.data(ksk + 576);
    const auto *ksk_579 = buffer.data(ksk + 579);
    const auto *ksk_582 = buffer.data(ksk + 582);
    const auto *ksk_586 = buffer.data(ksk + 586);
    const auto *ksk_591 = buffer.data(ksk + 591);
    const auto *ksk_604 = buffer.data(ksk + 604);
    const auto *ksk_611 = buffer.data(ksk + 611);
    const auto *ksk_612 = buffer.data(ksk + 612);
    const auto *ksk_614 = buffer.data(ksk + 614);
    const auto *ksk_615 = buffer.data(ksk + 615);
    const auto *ksk_617 = buffer.data(ksk + 617);
    const auto *ksk_618 = buffer.data(ksk + 618);
    const auto *ksk_621 = buffer.data(ksk + 621);
    const auto *ksk_622 = buffer.data(ksk + 622);
    const auto *ksk_626 = buffer.data(ksk + 626);
    const auto *ksk_627 = buffer.data(ksk + 627);
    const auto *ksk_632 = buffer.data(ksk + 632);
    const auto *ksk_640 = buffer.data(ksk + 640);
    const auto *ksk_642 = buffer.data(ksk + 642);
    const auto *ksk_643 = buffer.data(ksk + 643);
    const auto *ksk_644 = buffer.data(ksk + 644);
    const auto *ksk_645 = buffer.data(ksk + 645);
    const auto *ksk_646 = buffer.data(ksk + 646);
    const auto *ksk_647 = buffer.data(ksk + 647);
    const auto *ksk_648 = buffer.data(ksk + 648);
    const auto *ksk_650 = buffer.data(ksk + 650);
    const auto *ksk_651 = buffer.data(ksk + 651);
    const auto *ksk_653 = buffer.data(ksk + 653);
    const auto *ksk_654 = buffer.data(ksk + 654);
    const auto *ksk_657 = buffer.data(ksk + 657);
    const auto *ksk_662 = buffer.data(ksk + 662);
    const auto *ksk_668 = buffer.data(ksk + 668);
    const auto *ksk_676 = buffer.data(ksk + 676);
    const auto *ksk_678 = buffer.data(ksk + 678);
    const auto *ksk_679 = buffer.data(ksk + 679);
    const auto *ksk_680 = buffer.data(ksk + 680);
    const auto *ksk_681 = buffer.data(ksk + 681);
    const auto *ksk_682 = buffer.data(ksk + 682);
    const auto *ksk_683 = buffer.data(ksk + 683);
    const auto *ksk_684 = buffer.data(ksk + 684);
    const auto *ksk_686 = buffer.data(ksk + 686);
    const auto *ksk_689 = buffer.data(ksk + 689);
    const auto *ksk_693 = buffer.data(ksk + 693);
    const auto *ksk_831 = buffer.data(ksk + 831);
    const auto *ksk_833 = buffer.data(ksk + 833);
    const auto *ksk_834 = buffer.data(ksk + 834);
    const auto *ksk_837 = buffer.data(ksk + 837);
    const auto *ksk_838 = buffer.data(ksk + 838);
    const auto *ksk_840 = buffer.data(ksk + 840);
    const auto *ksk_842 = buffer.data(ksk + 842);
    const auto *ksk_843 = buffer.data(ksk + 843);
    const auto *ksk_845 = buffer.data(ksk + 845);
    const auto *ksk_846 = buffer.data(ksk + 846);
    const auto *ksk_848 = buffer.data(ksk + 848);
    const auto *ksk_849 = buffer.data(ksk + 849);
    const auto *ksk_851 = buffer.data(ksk + 851);
    const auto *ksk_852 = buffer.data(ksk + 852);
    const auto *ksk_853 = buffer.data(ksk + 853);
    const auto *ksk_855 = buffer.data(ksk + 855);
    const auto *ksk_856 = buffer.data(ksk + 856);
    const auto *ksk_857 = buffer.data(ksk + 857);
    const auto *ksk_858 = buffer.data(ksk + 858);
    const auto *ksk_859 = buffer.data(ksk + 859);
    const auto *ksk_860 = buffer.data(ksk + 860);
    const auto *ksk_861 = buffer.data(ksk + 861);
    const auto *ksk_862 = buffer.data(ksk + 862);
    const auto *ksk_863 = buffer.data(ksk + 863);
    const auto *ksk_864 = buffer.data(ksk + 864);
    const auto *ksk_867 = buffer.data(ksk + 867);
    const auto *ksk_869 = buffer.data(ksk + 869);
    const auto *ksk_870 = buffer.data(ksk + 870);
    const auto *ksk_873 = buffer.data(ksk + 873);
    const auto *ksk_874 = buffer.data(ksk + 874);
    const auto *ksk_876 = buffer.data(ksk + 876);
    const auto *ksk_878 = buffer.data(ksk + 878);
    const auto *ksk_879 = buffer.data(ksk + 879);
    const auto *ksk_881 = buffer.data(ksk + 881);
    const auto *ksk_882 = buffer.data(ksk + 882);
    const auto *ksk_884 = buffer.data(ksk + 884);
    const auto *ksk_885 = buffer.data(ksk + 885);
    const auto *ksk_887 = buffer.data(ksk + 887);
    const auto *ksk_888 = buffer.data(ksk + 888);
    const auto *ksk_889 = buffer.data(ksk + 889);
    const auto *ksk_891 = buffer.data(ksk + 891);
    const auto *ksk_892 = buffer.data(ksk + 892);
    const auto *ksk_893 = buffer.data(ksk + 893);
    const auto *ksk_894 = buffer.data(ksk + 894);
    const auto *ksk_895 = buffer.data(ksk + 895);
    const auto *ksk_896 = buffer.data(ksk + 896);
    const auto *ksk_897 = buffer.data(ksk + 897);
    const auto *ksk_898 = buffer.data(ksk + 898);
    const auto *ksk_899 = buffer.data(ksk + 899);
    const auto *ksk_900 = buffer.data(ksk + 900);
    const auto *ksk_903 = buffer.data(ksk + 903);
    const auto *ksk_905 = buffer.data(ksk + 905);
    const auto *ksk_906 = buffer.data(ksk + 906);
    const auto *ksk_909 = buffer.data(ksk + 909);
    const auto *ksk_910 = buffer.data(ksk + 910);
    const auto *ksk_912 = buffer.data(ksk + 912);

    const auto *lsi0_647 = buffer.data(lsi0 + 647);
    const auto *lsi0_649 = buffer.data(lsi0 + 649);
    const auto *lsi0_650 = buffer.data(lsi0 + 650);
    const auto *lsi0_653 = buffer.data(lsi0 + 653);
    const auto *lsi0_654 = buffer.data(lsi0 + 654);
    const auto *lsi0_656 = buffer.data(lsi0 + 656);
    const auto *lsi0_658 = buffer.data(lsi0 + 658);
    const auto *lsi0_659 = buffer.data(lsi0 + 659);
    const auto *lsi0_661 = buffer.data(lsi0 + 661);
    const auto *lsi0_662 = buffer.data(lsi0 + 662);
    const auto *lsi0_664 = buffer.data(lsi0 + 664);
    const auto *lsi0_665 = buffer.data(lsi0 + 665);
    const auto *lsi0_667 = buffer.data(lsi0 + 667);
    const auto *lsi0_668 = buffer.data(lsi0 + 668);
    const auto *lsi0_669 = buffer.data(lsi0 + 669);
    const auto *lsi0_670 = buffer.data(lsi0 + 670);
    const auto *lsi0_671 = buffer.data(lsi0 + 671);
    const auto *lsi0_672 = buffer.data(lsi0 + 672);
    const auto *lsi0_675 = buffer.data(lsi0 + 675);
    const auto *lsi0_677 = buffer.data(lsi0 + 677);
    const auto *lsi0_678 = buffer.data(lsi0 + 678);
    const auto *lsi0_681 = buffer.data(lsi0 + 681);
    const auto *lsi0_682 = buffer.data(lsi0 + 682);
    const auto *lsi0_684 = buffer.data(lsi0 + 684);
    const auto *lsi0_686 = buffer.data(lsi0 + 686);
    const auto *lsi0_687 = buffer.data(lsi0 + 687);
    const auto *lsi0_689 = buffer.data(lsi0 + 689);
    const auto *lsi0_690 = buffer.data(lsi0 + 690);
    const auto *lsi0_692 = buffer.data(lsi0 + 692);
    const auto *lsi0_693 = buffer.data(lsi0 + 693);
    const auto *lsi0_695 = buffer.data(lsi0 + 695);
    const auto *lsi0_696 = buffer.data(lsi0 + 696);
    const auto *lsi0_697 = buffer.data(lsi0 + 697);
    const auto *lsi0_698 = buffer.data(lsi0 + 698);
    const auto *lsi0_699 = buffer.data(lsi0 + 699);
    const auto *lsi0_700 = buffer.data(lsi0 + 700);
    const auto *lsi0_703 = buffer.data(lsi0 + 703);
    const auto *lsi0_705 = buffer.data(lsi0 + 705);
    const auto *lsi0_706 = buffer.data(lsi0 + 706);
    const auto *lsi0_709 = buffer.data(lsi0 + 709);
    const auto *lsi0_710 = buffer.data(lsi0 + 710);
    const auto *lsi0_712 = buffer.data(lsi0 + 712);

    const auto *lsi1_647 = buffer.data(lsi1 + 647);
    const auto *lsi1_649 = buffer.data(lsi1 + 649);
    const auto *lsi1_650 = buffer.data(lsi1 + 650);
    const auto *lsi1_653 = buffer.data(lsi1 + 653);
    const auto *lsi1_654 = buffer.data(lsi1 + 654);
    const auto *lsi1_656 = buffer.data(lsi1 + 656);
    const auto *lsi1_658 = buffer.data(lsi1 + 658);
    const auto *lsi1_659 = buffer.data(lsi1 + 659);
    const auto *lsi1_661 = buffer.data(lsi1 + 661);
    const auto *lsi1_662 = buffer.data(lsi1 + 662);
    const auto *lsi1_664 = buffer.data(lsi1 + 664);
    const auto *lsi1_665 = buffer.data(lsi1 + 665);
    const auto *lsi1_667 = buffer.data(lsi1 + 667);
    const auto *lsi1_668 = buffer.data(lsi1 + 668);
    const auto *lsi1_669 = buffer.data(lsi1 + 669);
    const auto *lsi1_670 = buffer.data(lsi1 + 670);
    const auto *lsi1_671 = buffer.data(lsi1 + 671);
    const auto *lsi1_672 = buffer.data(lsi1 + 672);
    const auto *lsi1_675 = buffer.data(lsi1 + 675);
    const auto *lsi1_677 = buffer.data(lsi1 + 677);
    const auto *lsi1_678 = buffer.data(lsi1 + 678);
    const auto *lsi1_681 = buffer.data(lsi1 + 681);
    const auto *lsi1_682 = buffer.data(lsi1 + 682);
    const auto *lsi1_684 = buffer.data(lsi1 + 684);
    const auto *lsi1_686 = buffer.data(lsi1 + 686);
    const auto *lsi1_687 = buffer.data(lsi1 + 687);
    const auto *lsi1_689 = buffer.data(lsi1 + 689);
    const auto *lsi1_690 = buffer.data(lsi1 + 690);
    const auto *lsi1_692 = buffer.data(lsi1 + 692);
    const auto *lsi1_693 = buffer.data(lsi1 + 693);
    const auto *lsi1_695 = buffer.data(lsi1 + 695);
    const auto *lsi1_696 = buffer.data(lsi1 + 696);
    const auto *lsi1_697 = buffer.data(lsi1 + 697);
    const auto *lsi1_698 = buffer.data(lsi1 + 698);
    const auto *lsi1_699 = buffer.data(lsi1 + 699);
    const auto *lsi1_700 = buffer.data(lsi1 + 700);
    const auto *lsi1_703 = buffer.data(lsi1 + 703);
    const auto *lsi1_705 = buffer.data(lsi1 + 705);
    const auto *lsi1_706 = buffer.data(lsi1 + 706);
    const auto *lsi1_709 = buffer.data(lsi1 + 709);
    const auto *lsi1_710 = buffer.data(lsi1 + 710);
    const auto *lsi1_712 = buffer.data(lsi1 + 712);

    const auto *lsk_828 = buffer.data(lsk + 828);
    const auto *lsk_830 = buffer.data(lsk + 830);
    const auto *lsk_831 = buffer.data(lsk + 831);
    const auto *lsk_833 = buffer.data(lsk + 833);
    const auto *lsk_834 = buffer.data(lsk + 834);
    const auto *lsk_837 = buffer.data(lsk + 837);
    const auto *lsk_838 = buffer.data(lsk + 838);
    const auto *lsk_840 = buffer.data(lsk + 840);
    const auto *lsk_842 = buffer.data(lsk + 842);
    const auto *lsk_843 = buffer.data(lsk + 843);
    const auto *lsk_845 = buffer.data(lsk + 845);
    const auto *lsk_846 = buffer.data(lsk + 846);
    const auto *lsk_848 = buffer.data(lsk + 848);
    const auto *lsk_849 = buffer.data(lsk + 849);
    const auto *lsk_851 = buffer.data(lsk + 851);
    const auto *lsk_852 = buffer.data(lsk + 852);
    const auto *lsk_853 = buffer.data(lsk + 853);
    const auto *lsk_855 = buffer.data(lsk + 855);
    const auto *lsk_856 = buffer.data(lsk + 856);
    const auto *lsk_857 = buffer.data(lsk + 857);
    const auto *lsk_858 = buffer.data(lsk + 858);
    const auto *lsk_859 = buffer.data(lsk + 859);
    const auto *lsk_860 = buffer.data(lsk + 860);
    const auto *lsk_861 = buffer.data(lsk + 861);
    const auto *lsk_862 = buffer.data(lsk + 862);
    const auto *lsk_863 = buffer.data(lsk + 863);
    const auto *lsk_864 = buffer.data(lsk + 864);
    const auto *lsk_866 = buffer.data(lsk + 866);
    const auto *lsk_867 = buffer.data(lsk + 867);
    const auto *lsk_869 = buffer.data(lsk + 869);
    const auto *lsk_870 = buffer.data(lsk + 870);
    const auto *lsk_873 = buffer.data(lsk + 873);
    const auto *lsk_874 = buffer.data(lsk + 874);
    const auto *lsk_876 = buffer.data(lsk + 876);
    const auto *lsk_878 = buffer.data(lsk + 878);
    const auto *lsk_879 = buffer.data(lsk + 879);
    const auto *lsk_881 = buffer.data(lsk + 881);
    const auto *lsk_882 = buffer.data(lsk + 882);
    const auto *lsk_884 = buffer.data(lsk + 884);
    const auto *lsk_885 = buffer.data(lsk + 885);
    const auto *lsk_887 = buffer.data(lsk + 887);
    const auto *lsk_888 = buffer.data(lsk + 888);
    const auto *lsk_889 = buffer.data(lsk + 889);
    const auto *lsk_891 = buffer.data(lsk + 891);
    const auto *lsk_892 = buffer.data(lsk + 892);
    const auto *lsk_893 = buffer.data(lsk + 893);
    const auto *lsk_894 = buffer.data(lsk + 894);
    const auto *lsk_895 = buffer.data(lsk + 895);
    const auto *lsk_896 = buffer.data(lsk + 896);
    const auto *lsk_897 = buffer.data(lsk + 897);
    const auto *lsk_898 = buffer.data(lsk + 898);
    const auto *lsk_899 = buffer.data(lsk + 899);
    const auto *lsk_900 = buffer.data(lsk + 900);
    const auto *lsk_902 = buffer.data(lsk + 902);
    const auto *lsk_903 = buffer.data(lsk + 903);
    const auto *lsk_905 = buffer.data(lsk + 905);
    const auto *lsk_906 = buffer.data(lsk + 906);
    const auto *lsk_909 = buffer.data(lsk + 909);
    const auto *lsk_910 = buffer.data(lsk + 910);
    const auto *lsk_912 = buffer.data(lsk + 912);

#pragma omp simd aligned(t_1037, t_1038, t_1039, pc_x, pc_y, pc_z, ksk_576, ksk_614, ksk_831, \
                         lsi0_647, lsi1_647, lsk_828, lsk_830, \
                         lsk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_16 * ksk_576[k]
                    + f_3 * pc_z[k] * lsk_828[k];

        t_1038[k] = f_16 * ksk_831[k]
                    + f_12 * lsi0_647[k]
                    - f_13 * lsi1_647[k]
                    + f_3 * pc_x[k] * lsk_831[k];

        t_1039[k] = f_18 * ksk_614[k]
                    + f_3 * pc_y[k] * lsk_830[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pc_x, pc_z, ksk_579, ksk_833, ksk_834, \
                         lsi0_649, lsi0_650, lsi1_649, lsi1_650, lsk_831, lsk_833, \
                         lsk_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_16 * ksk_833[k]
                    + f_12 * lsi0_649[k]
                    - f_13 * lsi1_649[k]
                    + f_3 * pc_x[k] * lsk_833[k];

        t_1041[k] = f_16 * ksk_834[k]
                    + f_10 * lsi0_650[k]
                    - f_11 * lsi1_650[k]
                    + f_3 * pc_x[k] * lsk_834[k];

        t_1042[k] = f_16 * ksk_579[k]
                    + f_3 * pc_z[k] * lsk_831[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pc_x, pc_y, ksk_617, ksk_837, ksk_838, \
                         lsi0_653, lsi0_654, lsi1_653, lsi1_654, lsk_833, lsk_837, \
                         lsk_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_18 * ksk_617[k]
                    + f_3 * pc_y[k] * lsk_833[k];

        t_1044[k] = f_16 * ksk_837[k]
                    + f_10 * lsi0_653[k]
                    - f_11 * lsi1_653[k]
                    + f_3 * pc_x[k] * lsk_837[k];

        t_1045[k] = f_16 * ksk_838[k]
                    + f_8 * lsi0_654[k]
                    - f_9 * lsi1_654[k]
                    + f_3 * pc_x[k] * lsk_838[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pc_x, pc_y, pc_z, ksk_582, ksk_621, ksk_840, \
                         lsi0_656, lsi1_656, lsk_834, lsk_837, \
                         lsk_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_16 * ksk_582[k]
                    + f_3 * pc_z[k] * lsk_834[k];

        t_1047[k] = f_16 * ksk_840[k]
                    + f_8 * lsi0_656[k]
                    - f_9 * lsi1_656[k]
                    + f_3 * pc_x[k] * lsk_840[k];

        t_1048[k] = f_18 * ksk_621[k]
                    + f_3 * pc_y[k] * lsk_837[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pc_x, pc_z, ksk_586, ksk_842, ksk_843, \
                         lsi0_658, lsi0_659, lsi1_658, lsi1_659, lsk_838, lsk_842, \
                         lsk_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_16 * ksk_842[k]
                    + f_8 * lsi0_658[k]
                    - f_9 * lsi1_658[k]
                    + f_3 * pc_x[k] * lsk_842[k];

        t_1050[k] = f_16 * ksk_843[k]
                    + f_6 * lsi0_659[k]
                    - f_7 * lsi1_659[k]
                    + f_3 * pc_x[k] * lsk_843[k];

        t_1051[k] = f_16 * ksk_586[k]
                    + f_3 * pc_z[k] * lsk_838[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pc_x, pc_y, ksk_626, ksk_845, ksk_846, \
                         lsi0_661, lsi0_662, lsi1_661, lsi1_662, lsk_842, lsk_845, \
                         lsk_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_16 * ksk_845[k]
                    + f_6 * lsi0_661[k]
                    - f_7 * lsi1_661[k]
                    + f_3 * pc_x[k] * lsk_845[k];

        t_1053[k] = f_16 * ksk_846[k]
                    + f_6 * lsi0_662[k]
                    - f_7 * lsi1_662[k]
                    + f_3 * pc_x[k] * lsk_846[k];

        t_1054[k] = f_18 * ksk_626[k]
                    + f_3 * pc_y[k] * lsk_842[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, pc_z, ksk_591, ksk_848, ksk_849, \
                         lsi0_664, lsi0_665, lsi1_664, lsi1_665, lsk_843, lsk_848, \
                         lsk_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_16 * ksk_848[k]
                    + f_6 * lsi0_664[k]
                    - f_7 * lsi1_664[k]
                    + f_3 * pc_x[k] * lsk_848[k];

        t_1056[k] = f_16 * ksk_849[k]
                    + f_4 * lsi0_665[k]
                    - f_5 * lsi1_665[k]
                    + f_3 * pc_x[k] * lsk_849[k];

        t_1057[k] = f_16 * ksk_591[k]
                    + f_3 * pc_z[k] * lsk_843[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_x, ksk_851, ksk_852, ksk_853, lsi0_667, \
                         lsi0_668, lsi0_669, lsi1_667, lsi1_668, lsi1_669, lsk_851, lsk_852, \
                         lsk_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_16 * ksk_851[k]
                    + f_4 * lsi0_667[k]
                    - f_5 * lsi1_667[k]
                    + f_3 * pc_x[k] * lsk_851[k];

        t_1059[k] = f_16 * ksk_852[k]
                    + f_4 * lsi0_668[k]
                    - f_5 * lsi1_668[k]
                    + f_3 * pc_x[k] * lsk_852[k];

        t_1060[k] = f_16 * ksk_853[k]
                    + f_4 * lsi0_669[k]
                    - f_5 * lsi1_669[k]
                    + f_3 * pc_x[k] * lsk_853[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, pc_x, pc_y, ksk_632, ksk_855, \
                         ksk_856, ksk_857, lsi0_671, lsi1_671, lsk_848, lsk_855, lsk_856, \
                         lsk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_18 * ksk_632[k]
                    + f_3 * pc_y[k] * lsk_848[k];

        t_1062[k] = f_16 * ksk_855[k]
                    + f_4 * lsi0_671[k]
                    - f_5 * lsi1_671[k]
                    + f_3 * pc_x[k] * lsk_855[k];

        t_1063[k] = f_16 * ksk_856[k]
                    + f_3 * pc_x[k] * lsk_856[k];

        t_1064[k] = f_16 * ksk_857[k]
                    + f_3 * pc_x[k] * lsk_857[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, pc_x, ksk_858, ksk_859, \
                         ksk_860, ksk_861, ksk_862, lsk_858, lsk_859, lsk_860, lsk_861, \
                         lsk_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_16 * ksk_858[k]
                    + f_3 * pc_x[k] * lsk_858[k];

        t_1066[k] = f_16 * ksk_859[k]
                    + f_3 * pc_x[k] * lsk_859[k];

        t_1067[k] = f_16 * ksk_860[k]
                    + f_3 * pc_x[k] * lsk_860[k];

        t_1068[k] = f_16 * ksk_861[k]
                    + f_3 * pc_x[k] * lsk_861[k];

        t_1069[k] = f_16 * ksk_862[k]
                    + f_3 * pc_x[k] * lsk_862[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, ksk_604, ksk_640, ksk_863, \
                         lsi0_665, lsi1_665, lsk_856, lsk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_16 * ksk_863[k]
                    + f_3 * pc_x[k] * lsk_863[k];

        t_1071[k] = f_18 * ksk_640[k]
                    + f_1 * lsi0_665[k]
                    - f_2 * lsi1_665[k]
                    + f_3 * pc_y[k] * lsk_856[k];

        t_1072[k] = f_16 * ksk_604[k]
                    + f_3 * pc_z[k] * lsk_856[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_y, ksk_642, ksk_643, ksk_644, lsi0_667, \
                         lsi0_668, lsi0_669, lsi1_667, lsi1_668, lsi1_669, lsk_858, lsk_859, \
                         lsk_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_18 * ksk_642[k]
                    + f_12 * lsi0_667[k]
                    - f_13 * lsi1_667[k]
                    + f_3 * pc_y[k] * lsk_858[k];

        t_1074[k] = f_18 * ksk_643[k]
                    + f_10 * lsi0_668[k]
                    - f_11 * lsi1_668[k]
                    + f_3 * pc_y[k] * lsk_859[k];

        t_1075[k] = f_18 * ksk_644[k]
                    + f_8 * lsi0_669[k]
                    - f_9 * lsi1_669[k]
                    + f_3 * pc_y[k] * lsk_860[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_y, ksk_645, ksk_646, ksk_647, lsi0_670, \
                         lsi0_671, lsi1_670, lsi1_671, lsk_861, lsk_862, \
                         lsk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_18 * ksk_645[k]
                    + f_6 * lsi0_670[k]
                    - f_7 * lsi1_670[k]
                    + f_3 * pc_y[k] * lsk_861[k];

        t_1077[k] = f_18 * ksk_646[k]
                    + f_4 * lsi0_671[k]
                    - f_5 * lsi1_671[k]
                    + f_3 * pc_y[k] * lsk_862[k];

        t_1078[k] = f_18 * ksk_647[k]
                    + f_3 * pc_y[k] * lsk_863[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, pc_x, pc_y, pc_z, ksk_611, ksk_648, ksk_864, \
                         lsi0_671, lsi0_672, lsi1_671, lsi1_672, lsk_863, \
                         lsk_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_16 * ksk_611[k]
                    + f_1 * lsi0_671[k]
                    - f_2 * lsi1_671[k]
                    + f_3 * pc_z[k] * lsk_863[k];

        t_1080[k] = f_16 * ksk_864[k]
                    + f_1 * lsi0_672[k]
                    - f_2 * lsi1_672[k]
                    + f_3 * pc_x[k] * lsk_864[k];

        t_1081[k] = f_17 * ksk_648[k]
                    + f_3 * pc_y[k] * lsk_864[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, pc_x, pc_y, pc_z, ksk_612, ksk_650, ksk_867, \
                         lsi0_675, lsi1_675, lsk_864, lsk_866, \
                         lsk_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_17 * ksk_612[k]
                    + f_3 * pc_z[k] * lsk_864[k];

        t_1083[k] = f_16 * ksk_867[k]
                    + f_12 * lsi0_675[k]
                    - f_13 * lsi1_675[k]
                    + f_3 * pc_x[k] * lsk_867[k];

        t_1084[k] = f_17 * ksk_650[k]
                    + f_3 * pc_y[k] * lsk_866[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_x, pc_z, ksk_615, ksk_869, ksk_870, \
                         lsi0_677, lsi0_678, lsi1_677, lsi1_678, lsk_867, lsk_869, \
                         lsk_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_16 * ksk_869[k]
                    + f_12 * lsi0_677[k]
                    - f_13 * lsi1_677[k]
                    + f_3 * pc_x[k] * lsk_869[k];

        t_1086[k] = f_16 * ksk_870[k]
                    + f_10 * lsi0_678[k]
                    - f_11 * lsi1_678[k]
                    + f_3 * pc_x[k] * lsk_870[k];

        t_1087[k] = f_17 * ksk_615[k]
                    + f_3 * pc_z[k] * lsk_867[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_x, pc_y, ksk_653, ksk_873, ksk_874, \
                         lsi0_681, lsi0_682, lsi1_681, lsi1_682, lsk_869, lsk_873, \
                         lsk_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * ksk_653[k]
                    + f_3 * pc_y[k] * lsk_869[k];

        t_1089[k] = f_16 * ksk_873[k]
                    + f_10 * lsi0_681[k]
                    - f_11 * lsi1_681[k]
                    + f_3 * pc_x[k] * lsk_873[k];

        t_1090[k] = f_16 * ksk_874[k]
                    + f_8 * lsi0_682[k]
                    - f_9 * lsi1_682[k]
                    + f_3 * pc_x[k] * lsk_874[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_y, pc_z, ksk_618, ksk_657, ksk_876, \
                         lsi0_684, lsi1_684, lsk_870, lsk_873, \
                         lsk_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_17 * ksk_618[k]
                    + f_3 * pc_z[k] * lsk_870[k];

        t_1092[k] = f_16 * ksk_876[k]
                    + f_8 * lsi0_684[k]
                    - f_9 * lsi1_684[k]
                    + f_3 * pc_x[k] * lsk_876[k];

        t_1093[k] = f_17 * ksk_657[k]
                    + f_3 * pc_y[k] * lsk_873[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, pc_z, ksk_622, ksk_878, ksk_879, \
                         lsi0_686, lsi0_687, lsi1_686, lsi1_687, lsk_874, lsk_878, \
                         lsk_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_16 * ksk_878[k]
                    + f_8 * lsi0_686[k]
                    - f_9 * lsi1_686[k]
                    + f_3 * pc_x[k] * lsk_878[k];

        t_1095[k] = f_16 * ksk_879[k]
                    + f_6 * lsi0_687[k]
                    - f_7 * lsi1_687[k]
                    + f_3 * pc_x[k] * lsk_879[k];

        t_1096[k] = f_17 * ksk_622[k]
                    + f_3 * pc_z[k] * lsk_874[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_y, ksk_662, ksk_881, ksk_882, \
                         lsi0_689, lsi0_690, lsi1_689, lsi1_690, lsk_878, lsk_881, \
                         lsk_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_16 * ksk_881[k]
                    + f_6 * lsi0_689[k]
                    - f_7 * lsi1_689[k]
                    + f_3 * pc_x[k] * lsk_881[k];

        t_1098[k] = f_16 * ksk_882[k]
                    + f_6 * lsi0_690[k]
                    - f_7 * lsi1_690[k]
                    + f_3 * pc_x[k] * lsk_882[k];

        t_1099[k] = f_17 * ksk_662[k]
                    + f_3 * pc_y[k] * lsk_878[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_z, ksk_627, ksk_884, ksk_885, \
                         lsi0_692, lsi0_693, lsi1_692, lsi1_693, lsk_879, lsk_884, \
                         lsk_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_16 * ksk_884[k]
                    + f_6 * lsi0_692[k]
                    - f_7 * lsi1_692[k]
                    + f_3 * pc_x[k] * lsk_884[k];

        t_1101[k] = f_16 * ksk_885[k]
                    + f_4 * lsi0_693[k]
                    - f_5 * lsi1_693[k]
                    + f_3 * pc_x[k] * lsk_885[k];

        t_1102[k] = f_17 * ksk_627[k]
                    + f_3 * pc_z[k] * lsk_879[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, ksk_887, ksk_888, ksk_889, lsi0_695, \
                         lsi0_696, lsi0_697, lsi1_695, lsi1_696, lsi1_697, lsk_887, lsk_888, \
                         lsk_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_16 * ksk_887[k]
                    + f_4 * lsi0_695[k]
                    - f_5 * lsi1_695[k]
                    + f_3 * pc_x[k] * lsk_887[k];

        t_1104[k] = f_16 * ksk_888[k]
                    + f_4 * lsi0_696[k]
                    - f_5 * lsi1_696[k]
                    + f_3 * pc_x[k] * lsk_888[k];

        t_1105[k] = f_16 * ksk_889[k]
                    + f_4 * lsi0_697[k]
                    - f_5 * lsi1_697[k]
                    + f_3 * pc_x[k] * lsk_889[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, pc_y, ksk_668, ksk_891, \
                         ksk_892, ksk_893, lsi0_699, lsi1_699, lsk_884, lsk_891, lsk_892, \
                         lsk_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_17 * ksk_668[k]
                    + f_3 * pc_y[k] * lsk_884[k];

        t_1107[k] = f_16 * ksk_891[k]
                    + f_4 * lsi0_699[k]
                    - f_5 * lsi1_699[k]
                    + f_3 * pc_x[k] * lsk_891[k];

        t_1108[k] = f_16 * ksk_892[k]
                    + f_3 * pc_x[k] * lsk_892[k];

        t_1109[k] = f_16 * ksk_893[k]
                    + f_3 * pc_x[k] * lsk_893[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, pc_x, ksk_894, ksk_895, \
                         ksk_896, ksk_897, ksk_898, lsk_894, lsk_895, lsk_896, lsk_897, \
                         lsk_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_16 * ksk_894[k]
                    + f_3 * pc_x[k] * lsk_894[k];

        t_1111[k] = f_16 * ksk_895[k]
                    + f_3 * pc_x[k] * lsk_895[k];

        t_1112[k] = f_16 * ksk_896[k]
                    + f_3 * pc_x[k] * lsk_896[k];

        t_1113[k] = f_16 * ksk_897[k]
                    + f_3 * pc_x[k] * lsk_897[k];

        t_1114[k] = f_16 * ksk_898[k]
                    + f_3 * pc_x[k] * lsk_898[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, pc_x, pc_y, pc_z, ksk_640, ksk_676, ksk_899, \
                         lsi0_693, lsi1_693, lsk_892, lsk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_16 * ksk_899[k]
                    + f_3 * pc_x[k] * lsk_899[k];

        t_1116[k] = f_17 * ksk_676[k]
                    + f_1 * lsi0_693[k]
                    - f_2 * lsi1_693[k]
                    + f_3 * pc_y[k] * lsk_892[k];

        t_1117[k] = f_17 * ksk_640[k]
                    + f_3 * pc_z[k] * lsk_892[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pc_y, ksk_678, ksk_679, ksk_680, lsi0_695, \
                         lsi0_696, lsi0_697, lsi1_695, lsi1_696, lsi1_697, lsk_894, lsk_895, \
                         lsk_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_17 * ksk_678[k]
                    + f_12 * lsi0_695[k]
                    - f_13 * lsi1_695[k]
                    + f_3 * pc_y[k] * lsk_894[k];

        t_1119[k] = f_17 * ksk_679[k]
                    + f_10 * lsi0_696[k]
                    - f_11 * lsi1_696[k]
                    + f_3 * pc_y[k] * lsk_895[k];

        t_1120[k] = f_17 * ksk_680[k]
                    + f_8 * lsi0_697[k]
                    - f_9 * lsi1_697[k]
                    + f_3 * pc_y[k] * lsk_896[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_y, ksk_681, ksk_682, ksk_683, lsi0_698, \
                         lsi0_699, lsi1_698, lsi1_699, lsk_897, lsk_898, \
                         lsk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_17 * ksk_681[k]
                    + f_6 * lsi0_698[k]
                    - f_7 * lsi1_698[k]
                    + f_3 * pc_y[k] * lsk_897[k];

        t_1122[k] = f_17 * ksk_682[k]
                    + f_4 * lsi0_699[k]
                    - f_5 * lsi1_699[k]
                    + f_3 * pc_y[k] * lsk_898[k];

        t_1123[k] = f_17 * ksk_683[k]
                    + f_3 * pc_y[k] * lsk_899[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_y, pc_z, ksk_647, ksk_684, ksk_900, \
                         lsi0_699, lsi0_700, lsi1_699, lsi1_700, lsk_899, \
                         lsk_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_17 * ksk_647[k]
                    + f_1 * lsi0_699[k]
                    - f_2 * lsi1_699[k]
                    + f_3 * pc_z[k] * lsk_899[k];

        t_1125[k] = f_16 * ksk_900[k]
                    + f_1 * lsi0_700[k]
                    - f_2 * lsi1_700[k]
                    + f_3 * pc_x[k] * lsk_900[k];

        t_1126[k] = f_16 * ksk_684[k]
                    + f_3 * pc_y[k] * lsk_900[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, pc_y, pc_z, ksk_648, ksk_686, ksk_903, \
                         lsi0_703, lsi1_703, lsk_900, lsk_902, \
                         lsk_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_18 * ksk_648[k]
                    + f_3 * pc_z[k] * lsk_900[k];

        t_1128[k] = f_16 * ksk_903[k]
                    + f_12 * lsi0_703[k]
                    - f_13 * lsi1_703[k]
                    + f_3 * pc_x[k] * lsk_903[k];

        t_1129[k] = f_16 * ksk_686[k]
                    + f_3 * pc_y[k] * lsk_902[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, pc_z, ksk_651, ksk_905, ksk_906, \
                         lsi0_705, lsi0_706, lsi1_705, lsi1_706, lsk_903, lsk_905, \
                         lsk_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_16 * ksk_905[k]
                    + f_12 * lsi0_705[k]
                    - f_13 * lsi1_705[k]
                    + f_3 * pc_x[k] * lsk_905[k];

        t_1131[k] = f_16 * ksk_906[k]
                    + f_10 * lsi0_706[k]
                    - f_11 * lsi1_706[k]
                    + f_3 * pc_x[k] * lsk_906[k];

        t_1132[k] = f_18 * ksk_651[k]
                    + f_3 * pc_z[k] * lsk_903[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, pc_y, ksk_689, ksk_909, ksk_910, \
                         lsi0_709, lsi0_710, lsi1_709, lsi1_710, lsk_905, lsk_909, \
                         lsk_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_16 * ksk_689[k]
                    + f_3 * pc_y[k] * lsk_905[k];

        t_1134[k] = f_16 * ksk_909[k]
                    + f_10 * lsi0_709[k]
                    - f_11 * lsi1_709[k]
                    + f_3 * pc_x[k] * lsk_909[k];

        t_1135[k] = f_16 * ksk_910[k]
                    + f_8 * lsi0_710[k]
                    - f_9 * lsi1_710[k]
                    + f_3 * pc_x[k] * lsk_910[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, pc_x, pc_y, pc_z, ksk_654, ksk_693, ksk_912, \
                         lsi0_712, lsi1_712, lsk_906, lsk_909, \
                         lsk_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_18 * ksk_654[k]
                    + f_3 * pc_z[k] * lsk_906[k];

        t_1137[k] = f_16 * ksk_912[k]
                    + f_8 * lsi0_712[k]
                    - f_9 * lsi1_712[k]
                    + f_3 * pc_x[k] * lsk_912[k];

        t_1138[k] = f_16 * ksk_693[k]
                    + f_3 * pc_y[k] * lsk_909[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksl0,
                                                           const size_t ksk, const size_t ksl1,
                                                           const size_t lsi0, const size_t lsi1,
                                                           const size_t lsk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_900 = buffer.data(ksl0 + 900);
    const auto *ksl0_903 = buffer.data(ksl0 + 903);
    const auto *ksl0_905 = buffer.data(ksl0 + 905);
    const auto *ksl0_906 = buffer.data(ksl0 + 906);
    const auto *ksl0_909 = buffer.data(ksl0 + 909);
    const auto *ksl0_910 = buffer.data(ksl0 + 910);
    const auto *ksl0_912 = buffer.data(ksl0 + 912);
    const auto *ksl0_914 = buffer.data(ksl0 + 914);
    const auto *ksl0_915 = buffer.data(ksl0 + 915);
    const auto *ksl0_917 = buffer.data(ksl0 + 917);
    const auto *ksl0_918 = buffer.data(ksl0 + 918);
    const auto *ksl0_920 = buffer.data(ksl0 + 920);
    const auto *ksl0_921 = buffer.data(ksl0 + 921);
    const auto *ksl0_923 = buffer.data(ksl0 + 923);
    const auto *ksl0_924 = buffer.data(ksl0 + 924);
    const auto *ksl0_925 = buffer.data(ksl0 + 925);
    const auto *ksl0_927 = buffer.data(ksl0 + 927);
    const auto *ksl0_944 = buffer.data(ksl0 + 944);

    const auto *ksk_658 = buffer.data(ksk + 658);
    const auto *ksk_663 = buffer.data(ksk + 663);
    const auto *ksk_676 = buffer.data(ksk + 676);
    const auto *ksk_683 = buffer.data(ksk + 683);
    const auto *ksk_684 = buffer.data(ksk + 684);
    const auto *ksk_687 = buffer.data(ksk + 687);
    const auto *ksk_690 = buffer.data(ksk + 690);
    const auto *ksk_694 = buffer.data(ksk + 694);
    const auto *ksk_698 = buffer.data(ksk + 698);
    const auto *ksk_699 = buffer.data(ksk + 699);
    const auto *ksk_704 = buffer.data(ksk + 704);
    const auto *ksk_712 = buffer.data(ksk + 712);
    const auto *ksk_714 = buffer.data(ksk + 714);
    const auto *ksk_715 = buffer.data(ksk + 715);
    const auto *ksk_716 = buffer.data(ksk + 716);
    const auto *ksk_717 = buffer.data(ksk + 717);
    const auto *ksk_718 = buffer.data(ksk + 718);
    const auto *ksk_719 = buffer.data(ksk + 719);
    const auto *ksk_720 = buffer.data(ksk + 720);
    const auto *ksk_721 = buffer.data(ksk + 721);
    const auto *ksk_722 = buffer.data(ksk + 722);
    const auto *ksk_723 = buffer.data(ksk + 723);
    const auto *ksk_725 = buffer.data(ksk + 725);
    const auto *ksk_726 = buffer.data(ksk + 726);
    const auto *ksk_728 = buffer.data(ksk + 728);
    const auto *ksk_729 = buffer.data(ksk + 729);
    const auto *ksk_730 = buffer.data(ksk + 730);
    const auto *ksk_732 = buffer.data(ksk + 732);
    const auto *ksk_733 = buffer.data(ksk + 733);
    const auto *ksk_734 = buffer.data(ksk + 734);
    const auto *ksk_735 = buffer.data(ksk + 735);
    const auto *ksk_737 = buffer.data(ksk + 737);
    const auto *ksk_738 = buffer.data(ksk + 738);
    const auto *ksk_739 = buffer.data(ksk + 739);
    const auto *ksk_740 = buffer.data(ksk + 740);
    const auto *ksk_748 = buffer.data(ksk + 748);
    const auto *ksk_750 = buffer.data(ksk + 750);
    const auto *ksk_751 = buffer.data(ksk + 751);
    const auto *ksk_752 = buffer.data(ksk + 752);
    const auto *ksk_753 = buffer.data(ksk + 753);
    const auto *ksk_754 = buffer.data(ksk + 754);
    const auto *ksk_755 = buffer.data(ksk + 755);
    const auto *ksk_914 = buffer.data(ksk + 914);
    const auto *ksk_915 = buffer.data(ksk + 915);
    const auto *ksk_917 = buffer.data(ksk + 917);
    const auto *ksk_918 = buffer.data(ksk + 918);
    const auto *ksk_920 = buffer.data(ksk + 920);
    const auto *ksk_921 = buffer.data(ksk + 921);
    const auto *ksk_923 = buffer.data(ksk + 923);
    const auto *ksk_924 = buffer.data(ksk + 924);
    const auto *ksk_925 = buffer.data(ksk + 925);
    const auto *ksk_927 = buffer.data(ksk + 927);
    const auto *ksk_928 = buffer.data(ksk + 928);
    const auto *ksk_929 = buffer.data(ksk + 929);
    const auto *ksk_930 = buffer.data(ksk + 930);
    const auto *ksk_931 = buffer.data(ksk + 931);
    const auto *ksk_932 = buffer.data(ksk + 932);
    const auto *ksk_933 = buffer.data(ksk + 933);
    const auto *ksk_934 = buffer.data(ksk + 934);
    const auto *ksk_935 = buffer.data(ksk + 935);
    const auto *ksk_964 = buffer.data(ksk + 964);
    const auto *ksk_965 = buffer.data(ksk + 965);
    const auto *ksk_966 = buffer.data(ksk + 966);
    const auto *ksk_967 = buffer.data(ksk + 967);
    const auto *ksk_968 = buffer.data(ksk + 968);
    const auto *ksk_969 = buffer.data(ksk + 969);
    const auto *ksk_970 = buffer.data(ksk + 970);
    const auto *ksk_971 = buffer.data(ksk + 971);
    const auto *ksk_972 = buffer.data(ksk + 972);
    const auto *ksk_977 = buffer.data(ksk + 977);
    const auto *ksk_981 = buffer.data(ksk + 981);
    const auto *ksk_986 = buffer.data(ksk + 986);
    const auto *ksk_992 = buffer.data(ksk + 992);
    const auto *ksk_999 = buffer.data(ksk + 999);
    const auto *ksk_1000 = buffer.data(ksk + 1000);
    const auto *ksk_1001 = buffer.data(ksk + 1001);
    const auto *ksk_1002 = buffer.data(ksk + 1002);
    const auto *ksk_1003 = buffer.data(ksk + 1003);
    const auto *ksk_1004 = buffer.data(ksk + 1004);
    const auto *ksk_1005 = buffer.data(ksk + 1005);
    const auto *ksk_1007 = buffer.data(ksk + 1007);

    const auto *ksl1_900 = buffer.data(ksl1 + 900);
    const auto *ksl1_903 = buffer.data(ksl1 + 903);
    const auto *ksl1_905 = buffer.data(ksl1 + 905);
    const auto *ksl1_906 = buffer.data(ksl1 + 906);
    const auto *ksl1_909 = buffer.data(ksl1 + 909);
    const auto *ksl1_910 = buffer.data(ksl1 + 910);
    const auto *ksl1_912 = buffer.data(ksl1 + 912);
    const auto *ksl1_914 = buffer.data(ksl1 + 914);
    const auto *ksl1_915 = buffer.data(ksl1 + 915);
    const auto *ksl1_917 = buffer.data(ksl1 + 917);
    const auto *ksl1_918 = buffer.data(ksl1 + 918);
    const auto *ksl1_920 = buffer.data(ksl1 + 920);
    const auto *ksl1_921 = buffer.data(ksl1 + 921);
    const auto *ksl1_923 = buffer.data(ksl1 + 923);
    const auto *ksl1_924 = buffer.data(ksl1 + 924);
    const auto *ksl1_925 = buffer.data(ksl1 + 925);
    const auto *ksl1_927 = buffer.data(ksl1 + 927);
    const auto *ksl1_944 = buffer.data(ksl1 + 944);

    const auto *lsi0_714 = buffer.data(lsi0 + 714);
    const auto *lsi0_715 = buffer.data(lsi0 + 715);
    const auto *lsi0_717 = buffer.data(lsi0 + 717);
    const auto *lsi0_718 = buffer.data(lsi0 + 718);
    const auto *lsi0_720 = buffer.data(lsi0 + 720);
    const auto *lsi0_721 = buffer.data(lsi0 + 721);
    const auto *lsi0_723 = buffer.data(lsi0 + 723);
    const auto *lsi0_724 = buffer.data(lsi0 + 724);
    const auto *lsi0_725 = buffer.data(lsi0 + 725);
    const auto *lsi0_726 = buffer.data(lsi0 + 726);
    const auto *lsi0_727 = buffer.data(lsi0 + 727);
    const auto *lsi0_749 = buffer.data(lsi0 + 749);
    const auto *lsi0_751 = buffer.data(lsi0 + 751);
    const auto *lsi0_752 = buffer.data(lsi0 + 752);
    const auto *lsi0_753 = buffer.data(lsi0 + 753);
    const auto *lsi0_754 = buffer.data(lsi0 + 754);
    const auto *lsi0_755 = buffer.data(lsi0 + 755);
    const auto *lsi0_756 = buffer.data(lsi0 + 756);
    const auto *lsi0_757 = buffer.data(lsi0 + 757);
    const auto *lsi0_758 = buffer.data(lsi0 + 758);
    const auto *lsi0_759 = buffer.data(lsi0 + 759);
    const auto *lsi0_760 = buffer.data(lsi0 + 760);
    const auto *lsi0_761 = buffer.data(lsi0 + 761);
    const auto *lsi0_762 = buffer.data(lsi0 + 762);
    const auto *lsi0_763 = buffer.data(lsi0 + 763);
    const auto *lsi0_764 = buffer.data(lsi0 + 764);
    const auto *lsi0_765 = buffer.data(lsi0 + 765);
    const auto *lsi0_766 = buffer.data(lsi0 + 766);
    const auto *lsi0_767 = buffer.data(lsi0 + 767);
    const auto *lsi0_768 = buffer.data(lsi0 + 768);
    const auto *lsi0_769 = buffer.data(lsi0 + 769);
    const auto *lsi0_770 = buffer.data(lsi0 + 770);
    const auto *lsi0_776 = buffer.data(lsi0 + 776);
    const auto *lsi0_783 = buffer.data(lsi0 + 783);

    const auto *lsi1_714 = buffer.data(lsi1 + 714);
    const auto *lsi1_715 = buffer.data(lsi1 + 715);
    const auto *lsi1_717 = buffer.data(lsi1 + 717);
    const auto *lsi1_718 = buffer.data(lsi1 + 718);
    const auto *lsi1_720 = buffer.data(lsi1 + 720);
    const auto *lsi1_721 = buffer.data(lsi1 + 721);
    const auto *lsi1_723 = buffer.data(lsi1 + 723);
    const auto *lsi1_724 = buffer.data(lsi1 + 724);
    const auto *lsi1_725 = buffer.data(lsi1 + 725);
    const auto *lsi1_726 = buffer.data(lsi1 + 726);
    const auto *lsi1_727 = buffer.data(lsi1 + 727);
    const auto *lsi1_749 = buffer.data(lsi1 + 749);
    const auto *lsi1_751 = buffer.data(lsi1 + 751);
    const auto *lsi1_752 = buffer.data(lsi1 + 752);
    const auto *lsi1_753 = buffer.data(lsi1 + 753);
    const auto *lsi1_754 = buffer.data(lsi1 + 754);
    const auto *lsi1_755 = buffer.data(lsi1 + 755);
    const auto *lsi1_756 = buffer.data(lsi1 + 756);
    const auto *lsi1_757 = buffer.data(lsi1 + 757);
    const auto *lsi1_758 = buffer.data(lsi1 + 758);
    const auto *lsi1_759 = buffer.data(lsi1 + 759);
    const auto *lsi1_760 = buffer.data(lsi1 + 760);
    const auto *lsi1_761 = buffer.data(lsi1 + 761);
    const auto *lsi1_762 = buffer.data(lsi1 + 762);
    const auto *lsi1_763 = buffer.data(lsi1 + 763);
    const auto *lsi1_764 = buffer.data(lsi1 + 764);
    const auto *lsi1_765 = buffer.data(lsi1 + 765);
    const auto *lsi1_766 = buffer.data(lsi1 + 766);
    const auto *lsi1_767 = buffer.data(lsi1 + 767);
    const auto *lsi1_768 = buffer.data(lsi1 + 768);
    const auto *lsi1_769 = buffer.data(lsi1 + 769);
    const auto *lsi1_770 = buffer.data(lsi1 + 770);
    const auto *lsi1_776 = buffer.data(lsi1 + 776);
    const auto *lsi1_783 = buffer.data(lsi1 + 783);

    const auto *lsk_910 = buffer.data(lsk + 910);
    const auto *lsk_914 = buffer.data(lsk + 914);
    const auto *lsk_915 = buffer.data(lsk + 915);
    const auto *lsk_917 = buffer.data(lsk + 917);
    const auto *lsk_918 = buffer.data(lsk + 918);
    const auto *lsk_920 = buffer.data(lsk + 920);
    const auto *lsk_921 = buffer.data(lsk + 921);
    const auto *lsk_923 = buffer.data(lsk + 923);
    const auto *lsk_924 = buffer.data(lsk + 924);
    const auto *lsk_925 = buffer.data(lsk + 925);
    const auto *lsk_927 = buffer.data(lsk + 927);
    const auto *lsk_928 = buffer.data(lsk + 928);
    const auto *lsk_929 = buffer.data(lsk + 929);
    const auto *lsk_930 = buffer.data(lsk + 930);
    const auto *lsk_931 = buffer.data(lsk + 931);
    const auto *lsk_932 = buffer.data(lsk + 932);
    const auto *lsk_933 = buffer.data(lsk + 933);
    const auto *lsk_934 = buffer.data(lsk + 934);
    const auto *lsk_935 = buffer.data(lsk + 935);
    const auto *lsk_936 = buffer.data(lsk + 936);
    const auto *lsk_938 = buffer.data(lsk + 938);
    const auto *lsk_939 = buffer.data(lsk + 939);
    const auto *lsk_941 = buffer.data(lsk + 941);
    const auto *lsk_942 = buffer.data(lsk + 942);
    const auto *lsk_945 = buffer.data(lsk + 945);
    const auto *lsk_946 = buffer.data(lsk + 946);
    const auto *lsk_950 = buffer.data(lsk + 950);
    const auto *lsk_951 = buffer.data(lsk + 951);
    const auto *lsk_956 = buffer.data(lsk + 956);
    const auto *lsk_964 = buffer.data(lsk + 964);
    const auto *lsk_965 = buffer.data(lsk + 965);
    const auto *lsk_966 = buffer.data(lsk + 966);
    const auto *lsk_967 = buffer.data(lsk + 967);
    const auto *lsk_968 = buffer.data(lsk + 968);
    const auto *lsk_969 = buffer.data(lsk + 969);
    const auto *lsk_970 = buffer.data(lsk + 970);
    const auto *lsk_971 = buffer.data(lsk + 971);
    const auto *lsk_972 = buffer.data(lsk + 972);
    const auto *lsk_973 = buffer.data(lsk + 973);
    const auto *lsk_974 = buffer.data(lsk + 974);
    const auto *lsk_975 = buffer.data(lsk + 975);
    const auto *lsk_976 = buffer.data(lsk + 976);
    const auto *lsk_977 = buffer.data(lsk + 977);
    const auto *lsk_978 = buffer.data(lsk + 978);
    const auto *lsk_979 = buffer.data(lsk + 979);
    const auto *lsk_980 = buffer.data(lsk + 980);
    const auto *lsk_981 = buffer.data(lsk + 981);
    const auto *lsk_982 = buffer.data(lsk + 982);
    const auto *lsk_983 = buffer.data(lsk + 983);
    const auto *lsk_984 = buffer.data(lsk + 984);
    const auto *lsk_985 = buffer.data(lsk + 985);
    const auto *lsk_986 = buffer.data(lsk + 986);
    const auto *lsk_987 = buffer.data(lsk + 987);
    const auto *lsk_988 = buffer.data(lsk + 988);
    const auto *lsk_989 = buffer.data(lsk + 989);
    const auto *lsk_990 = buffer.data(lsk + 990);
    const auto *lsk_991 = buffer.data(lsk + 991);
    const auto *lsk_992 = buffer.data(lsk + 992);
    const auto *lsk_999 = buffer.data(lsk + 999);
    const auto *lsk_1000 = buffer.data(lsk + 1000);
    const auto *lsk_1001 = buffer.data(lsk + 1001);
    const auto *lsk_1002 = buffer.data(lsk + 1002);
    const auto *lsk_1003 = buffer.data(lsk + 1003);
    const auto *lsk_1004 = buffer.data(lsk + 1004);
    const auto *lsk_1005 = buffer.data(lsk + 1005);
    const auto *lsk_1007 = buffer.data(lsk + 1007);

#pragma omp simd aligned(t_1139, t_1140, t_1141, pc_x, pc_z, ksk_658, ksk_914, ksk_915, \
                         lsi0_714, lsi0_715, lsi1_714, lsi1_715, lsk_910, lsk_914, \
                         lsk_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_16 * ksk_914[k]
                    + f_8 * lsi0_714[k]
                    - f_9 * lsi1_714[k]
                    + f_3 * pc_x[k] * lsk_914[k];

        t_1140[k] = f_16 * ksk_915[k]
                    + f_6 * lsi0_715[k]
                    - f_7 * lsi1_715[k]
                    + f_3 * pc_x[k] * lsk_915[k];

        t_1141[k] = f_18 * ksk_658[k]
                    + f_3 * pc_z[k] * lsk_910[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, pc_x, pc_y, ksk_698, ksk_917, ksk_918, \
                         lsi0_717, lsi0_718, lsi1_717, lsi1_718, lsk_914, lsk_917, \
                         lsk_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_16 * ksk_917[k]
                    + f_6 * lsi0_717[k]
                    - f_7 * lsi1_717[k]
                    + f_3 * pc_x[k] * lsk_917[k];

        t_1143[k] = f_16 * ksk_918[k]
                    + f_6 * lsi0_718[k]
                    - f_7 * lsi1_718[k]
                    + f_3 * pc_x[k] * lsk_918[k];

        t_1144[k] = f_16 * ksk_698[k]
                    + f_3 * pc_y[k] * lsk_914[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, pc_x, pc_z, ksk_663, ksk_920, ksk_921, \
                         lsi0_720, lsi0_721, lsi1_720, lsi1_721, lsk_915, lsk_920, \
                         lsk_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_16 * ksk_920[k]
                    + f_6 * lsi0_720[k]
                    - f_7 * lsi1_720[k]
                    + f_3 * pc_x[k] * lsk_920[k];

        t_1146[k] = f_16 * ksk_921[k]
                    + f_4 * lsi0_721[k]
                    - f_5 * lsi1_721[k]
                    + f_3 * pc_x[k] * lsk_921[k];

        t_1147[k] = f_18 * ksk_663[k]
                    + f_3 * pc_z[k] * lsk_915[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, pc_x, ksk_923, ksk_924, ksk_925, lsi0_723, \
                         lsi0_724, lsi0_725, lsi1_723, lsi1_724, lsi1_725, lsk_923, lsk_924, \
                         lsk_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = f_16 * ksk_923[k]
                    + f_4 * lsi0_723[k]
                    - f_5 * lsi1_723[k]
                    + f_3 * pc_x[k] * lsk_923[k];

        t_1149[k] = f_16 * ksk_924[k]
                    + f_4 * lsi0_724[k]
                    - f_5 * lsi1_724[k]
                    + f_3 * pc_x[k] * lsk_924[k];

        t_1150[k] = f_16 * ksk_925[k]
                    + f_4 * lsi0_725[k]
                    - f_5 * lsi1_725[k]
                    + f_3 * pc_x[k] * lsk_925[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, pc_x, pc_y, ksk_704, ksk_927, \
                         ksk_928, ksk_929, lsi0_727, lsi1_727, lsk_920, lsk_927, lsk_928, \
                         lsk_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_16 * ksk_704[k]
                    + f_3 * pc_y[k] * lsk_920[k];

        t_1152[k] = f_16 * ksk_927[k]
                    + f_4 * lsi0_727[k]
                    - f_5 * lsi1_727[k]
                    + f_3 * pc_x[k] * lsk_927[k];

        t_1153[k] = f_16 * ksk_928[k]
                    + f_3 * pc_x[k] * lsk_928[k];

        t_1154[k] = f_16 * ksk_929[k]
                    + f_3 * pc_x[k] * lsk_929[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, pc_x, ksk_930, ksk_931, \
                         ksk_932, ksk_933, ksk_934, lsk_930, lsk_931, lsk_932, lsk_933, \
                         lsk_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_16 * ksk_930[k]
                    + f_3 * pc_x[k] * lsk_930[k];

        t_1156[k] = f_16 * ksk_931[k]
                    + f_3 * pc_x[k] * lsk_931[k];

        t_1157[k] = f_16 * ksk_932[k]
                    + f_3 * pc_x[k] * lsk_932[k];

        t_1158[k] = f_16 * ksk_933[k]
                    + f_3 * pc_x[k] * lsk_933[k];

        t_1159[k] = f_16 * ksk_934[k]
                    + f_3 * pc_x[k] * lsk_934[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, pc_y, pc_z, ksk_676, ksk_712, ksk_935, \
                         lsi0_721, lsi1_721, lsk_928, lsk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_16 * ksk_935[k]
                    + f_3 * pc_x[k] * lsk_935[k];

        t_1161[k] = f_16 * ksk_712[k]
                    + f_1 * lsi0_721[k]
                    - f_2 * lsi1_721[k]
                    + f_3 * pc_y[k] * lsk_928[k];

        t_1162[k] = f_18 * ksk_676[k]
                    + f_3 * pc_z[k] * lsk_928[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_y, ksk_714, ksk_715, ksk_716, lsi0_723, \
                         lsi0_724, lsi0_725, lsi1_723, lsi1_724, lsi1_725, lsk_930, lsk_931, \
                         lsk_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * ksk_714[k]
                    + f_12 * lsi0_723[k]
                    - f_13 * lsi1_723[k]
                    + f_3 * pc_y[k] * lsk_930[k];

        t_1164[k] = f_16 * ksk_715[k]
                    + f_10 * lsi0_724[k]
                    - f_11 * lsi1_724[k]
                    + f_3 * pc_y[k] * lsk_931[k];

        t_1165[k] = f_16 * ksk_716[k]
                    + f_8 * lsi0_725[k]
                    - f_9 * lsi1_725[k]
                    + f_3 * pc_y[k] * lsk_932[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_y, ksk_717, ksk_718, ksk_719, lsi0_726, \
                         lsi0_727, lsi1_726, lsi1_727, lsk_933, lsk_934, \
                         lsk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_16 * ksk_717[k]
                    + f_6 * lsi0_726[k]
                    - f_7 * lsi1_726[k]
                    + f_3 * pc_y[k] * lsk_933[k];

        t_1167[k] = f_16 * ksk_718[k]
                    + f_4 * lsi0_727[k]
                    - f_5 * lsi1_727[k]
                    + f_3 * pc_y[k] * lsk_934[k];

        t_1168[k] = f_16 * ksk_719[k]
                    + f_3 * pc_y[k] * lsk_935[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, pa_y, pc_y, pc_z, ksl0_900, ksk_683, \
                         ksk_684, ksk_720, ksl1_900, lsi0_727, lsi1_727, lsk_935, \
                         lsk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_18 * ksk_683[k]
                    + f_1 * lsi0_727[k]
                    - f_2 * lsi1_727[k]
                    + f_3 * pc_z[k] * lsk_935[k];

        t_1170[k] = pa_y[k] * ksl0_900[k]
                    - f_14 * pc_y[k] * ksl1_900[k];

        t_1171[k] = f_15 * ksk_720[k]
                    + f_3 * pc_y[k] * lsk_936[k];

        t_1172[k] = f_19 * ksk_684[k]
                    + f_3 * pc_z[k] * lsk_936[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, t_1176, pa_y, pc_y, ksl0_903, ksl0_905, \
                         ksl0_906, ksk_721, ksk_722, ksk_723, ksl1_903, ksl1_905, ksl1_906, \
                         lsk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = pa_y[k] * ksl0_903[k]
                    + f_16 * ksk_721[k]
                    - f_14 * pc_y[k] * ksl1_903[k];

        t_1174[k] = f_15 * ksk_722[k]
                    + f_3 * pc_y[k] * lsk_938[k];

        t_1175[k] = pa_y[k] * ksl0_905[k]
                    - f_14 * pc_y[k] * ksl1_905[k];

        t_1176[k] = pa_y[k] * ksl0_906[k]
                    + f_17 * ksk_723[k]
                    - f_14 * pc_y[k] * ksl1_906[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, pa_y, pc_y, pc_z, ksl0_909, ksl0_910, \
                         ksk_687, ksk_725, ksk_726, ksl1_909, ksl1_910, lsk_939, \
                         lsk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_19 * ksk_687[k]
                    + f_3 * pc_z[k] * lsk_939[k];

        t_1178[k] = f_15 * ksk_725[k]
                    + f_3 * pc_y[k] * lsk_941[k];

        t_1179[k] = pa_y[k] * ksl0_909[k]
                    - f_14 * pc_y[k] * ksl1_909[k];

        t_1180[k] = pa_y[k] * ksl0_910[k]
                    + f_18 * ksk_726[k]
                    - f_14 * pc_y[k] * ksl1_910[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, t_1184, pa_y, pc_y, pc_z, ksl0_912, ksl0_914, \
                         ksk_690, ksk_728, ksk_729, ksl1_912, ksl1_914, lsk_942, \
                         lsk_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_19 * ksk_690[k]
                    + f_3 * pc_z[k] * lsk_942[k];

        t_1182[k] = pa_y[k] * ksl0_912[k]
                    + f_16 * ksk_728[k]
                    - f_14 * pc_y[k] * ksl1_912[k];

        t_1183[k] = f_15 * ksk_729[k]
                    + f_3 * pc_y[k] * lsk_945[k];

        t_1184[k] = pa_y[k] * ksl0_914[k]
                    - f_14 * pc_y[k] * ksl1_914[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pa_y, pc_y, pc_z, ksl0_915, ksl0_917, \
                         ksk_694, ksk_730, ksk_732, ksl1_915, ksl1_917, \
                         lsk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = pa_y[k] * ksl0_915[k]
                    + f_19 * ksk_730[k]
                    - f_14 * pc_y[k] * ksl1_915[k];

        t_1186[k] = f_19 * ksk_694[k]
                    + f_3 * pc_z[k] * lsk_946[k];

        t_1187[k] = pa_y[k] * ksl0_917[k]
                    + f_17 * ksk_732[k]
                    - f_14 * pc_y[k] * ksl1_917[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, t_1191, pa_y, pc_y, ksl0_918, ksl0_920, \
                         ksl0_921, ksk_733, ksk_734, ksk_735, ksl1_918, ksl1_920, ksl1_921, \
                         lsk_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = pa_y[k] * ksl0_918[k]
                    + f_16 * ksk_733[k]
                    - f_14 * pc_y[k] * ksl1_918[k];

        t_1189[k] = f_15 * ksk_734[k]
                    + f_3 * pc_y[k] * lsk_950[k];

        t_1190[k] = pa_y[k] * ksl0_920[k]
                    - f_14 * pc_y[k] * ksl1_920[k];

        t_1191[k] = pa_y[k] * ksl0_921[k]
                    + f_20 * ksk_735[k]
                    - f_14 * pc_y[k] * ksl1_921[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, pa_y, pc_y, pc_z, ksl0_923, ksl0_924, \
                         ksk_699, ksk_737, ksk_738, ksl1_923, ksl1_924, \
                         lsk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = f_19 * ksk_699[k]
                    + f_3 * pc_z[k] * lsk_951[k];

        t_1193[k] = pa_y[k] * ksl0_923[k]
                    + f_18 * ksk_737[k]
                    - f_14 * pc_y[k] * ksl1_923[k];

        t_1194[k] = pa_y[k] * ksl0_924[k]
                    + f_17 * ksk_738[k]
                    - f_14 * pc_y[k] * ksl1_924[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, pa_y, pc_x, pc_y, ksl0_925, ksl0_927, \
                         ksk_739, ksk_740, ksk_964, ksl1_925, ksl1_927, lsk_956, \
                         lsk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = pa_y[k] * ksl0_925[k]
                    + f_16 * ksk_739[k]
                    - f_14 * pc_y[k] * ksl1_925[k];

        t_1196[k] = f_15 * ksk_740[k]
                    + f_3 * pc_y[k] * lsk_956[k];

        t_1197[k] = pa_y[k] * ksl0_927[k]
                    - f_14 * pc_y[k] * ksl1_927[k];

        t_1198[k] = f_16 * ksk_964[k]
                    + f_3 * pc_x[k] * lsk_964[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, t_1202, t_1203, pc_x, ksk_965, ksk_966, \
                         ksk_967, ksk_968, ksk_969, lsk_965, lsk_966, lsk_967, lsk_968, \
                         lsk_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_16 * ksk_965[k]
                    + f_3 * pc_x[k] * lsk_965[k];

        t_1200[k] = f_16 * ksk_966[k]
                    + f_3 * pc_x[k] * lsk_966[k];

        t_1201[k] = f_16 * ksk_967[k]
                    + f_3 * pc_x[k] * lsk_967[k];

        t_1202[k] = f_16 * ksk_968[k]
                    + f_3 * pc_x[k] * lsk_968[k];

        t_1203[k] = f_16 * ksk_969[k]
                    + f_3 * pc_x[k] * lsk_969[k];
    }

#pragma omp simd aligned(t_1204, t_1205, t_1206, t_1207, pc_x, pc_y, pc_z, ksk_712, ksk_748, \
                         ksk_970, ksk_971, lsi0_749, lsi1_749, lsk_964, lsk_970, \
                         lsk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1204[k] = f_16 * ksk_970[k]
                    + f_3 * pc_x[k] * lsk_970[k];

        t_1205[k] = f_16 * ksk_971[k]
                    + f_3 * pc_x[k] * lsk_971[k];

        t_1206[k] = f_15 * ksk_748[k]
                    + f_1 * lsi0_749[k]
                    - f_2 * lsi1_749[k]
                    + f_3 * pc_y[k] * lsk_964[k];

        t_1207[k] = f_19 * ksk_712[k]
                    + f_3 * pc_z[k] * lsk_964[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, pc_y, ksk_750, ksk_751, ksk_752, lsi0_751, \
                         lsi0_752, lsi0_753, lsi1_751, lsi1_752, lsi1_753, lsk_966, lsk_967, \
                         lsk_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_15 * ksk_750[k]
                    + f_12 * lsi0_751[k]
                    - f_13 * lsi1_751[k]
                    + f_3 * pc_y[k] * lsk_966[k];

        t_1209[k] = f_15 * ksk_751[k]
                    + f_10 * lsi0_752[k]
                    - f_11 * lsi1_752[k]
                    + f_3 * pc_y[k] * lsk_967[k];

        t_1210[k] = f_15 * ksk_752[k]
                    + f_8 * lsi0_753[k]
                    - f_9 * lsi1_753[k]
                    + f_3 * pc_y[k] * lsk_968[k];
    }

#pragma omp simd aligned(t_1211, t_1212, t_1213, pc_y, ksk_753, ksk_754, ksk_755, lsi0_754, \
                         lsi0_755, lsi1_754, lsi1_755, lsk_969, lsk_970, \
                         lsk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = f_15 * ksk_753[k]
                    + f_6 * lsi0_754[k]
                    - f_7 * lsi1_754[k]
                    + f_3 * pc_y[k] * lsk_969[k];

        t_1212[k] = f_15 * ksk_754[k]
                    + f_4 * lsi0_755[k]
                    - f_5 * lsi1_755[k]
                    + f_3 * pc_y[k] * lsk_970[k];

        t_1213[k] = f_15 * ksk_755[k]
                    + f_3 * pc_y[k] * lsk_971[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pa_y, pc_x, pc_y, pc_z, ksl0_944, \
                         ksk_720, ksk_972, ksl1_944, lsi0_756, lsi1_756, \
                         lsk_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pa_y[k] * ksl0_944[k]
                    - f_14 * pc_y[k] * ksl1_944[k];

        t_1215[k] = f_16 * ksk_972[k]
                    + f_1 * lsi0_756[k]
                    - f_2 * lsi1_756[k]
                    + f_3 * pc_x[k] * lsk_972[k];

        t_1216[k] = f_3 * pc_y[k] * lsk_972[k];

        t_1217[k] = f_20 * ksk_720[k]
                    + f_3 * pc_z[k] * lsk_972[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, pc_x, pc_y, ksk_977, lsi0_756, lsi0_761, \
                         lsi1_756, lsi1_761, lsk_973, lsk_974, \
                         lsk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = f_4 * lsi0_756[k]
                    - f_5 * lsi1_756[k]
                    + f_3 * pc_y[k] * lsk_973[k];

        t_1219[k] = f_3 * pc_y[k] * lsk_974[k];

        t_1220[k] = f_16 * ksk_977[k]
                    + f_12 * lsi0_761[k]
                    - f_13 * lsi1_761[k]
                    + f_3 * pc_x[k] * lsk_977[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, pc_y, lsi0_757, lsi0_758, lsi1_757, lsi1_758, \
                         lsk_975, lsk_976, lsk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_6 * lsi0_757[k]
                    - f_7 * lsi1_757[k]
                    + f_3 * pc_y[k] * lsk_975[k];

        t_1222[k] = f_4 * lsi0_758[k]
                    - f_5 * lsi1_758[k]
                    + f_3 * pc_y[k] * lsk_976[k];

        t_1223[k] = f_3 * pc_y[k] * lsk_977[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pc_x, pc_y, ksk_981, lsi0_759, lsi0_760, \
                         lsi0_765, lsi1_759, lsi1_760, lsi1_765, lsk_978, lsk_979, \
                         lsk_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = f_16 * ksk_981[k]
                    + f_10 * lsi0_765[k]
                    - f_11 * lsi1_765[k]
                    + f_3 * pc_x[k] * lsk_981[k];

        t_1225[k] = f_8 * lsi0_759[k]
                    - f_9 * lsi1_759[k]
                    + f_3 * pc_y[k] * lsk_978[k];

        t_1226[k] = f_6 * lsi0_760[k]
                    - f_7 * lsi1_760[k]
                    + f_3 * pc_y[k] * lsk_979[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pc_x, pc_y, ksk_986, lsi0_761, lsi0_770, \
                         lsi1_761, lsi1_770, lsk_980, lsk_981, \
                         lsk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_4 * lsi0_761[k]
                    - f_5 * lsi1_761[k]
                    + f_3 * pc_y[k] * lsk_980[k];

        t_1228[k] = f_3 * pc_y[k] * lsk_981[k];

        t_1229[k] = f_16 * ksk_986[k]
                    + f_8 * lsi0_770[k]
                    - f_9 * lsi1_770[k]
                    + f_3 * pc_x[k] * lsk_986[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, pc_y, lsi0_762, lsi0_763, lsi0_764, lsi1_762, \
                         lsi1_763, lsi1_764, lsk_982, lsk_983, \
                         lsk_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_10 * lsi0_762[k]
                    - f_11 * lsi1_762[k]
                    + f_3 * pc_y[k] * lsk_982[k];

        t_1231[k] = f_8 * lsi0_763[k]
                    - f_9 * lsi1_763[k]
                    + f_3 * pc_y[k] * lsk_983[k];

        t_1232[k] = f_6 * lsi0_764[k]
                    - f_7 * lsi1_764[k]
                    + f_3 * pc_y[k] * lsk_984[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, pc_x, pc_y, ksk_992, lsi0_765, lsi0_776, \
                         lsi1_765, lsi1_776, lsk_985, lsk_986, \
                         lsk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = f_4 * lsi0_765[k]
                    - f_5 * lsi1_765[k]
                    + f_3 * pc_y[k] * lsk_985[k];

        t_1234[k] = f_3 * pc_y[k] * lsk_986[k];

        t_1235[k] = f_16 * ksk_992[k]
                    + f_6 * lsi0_776[k]
                    - f_7 * lsi1_776[k]
                    + f_3 * pc_x[k] * lsk_992[k];
    }

#pragma omp simd aligned(t_1236, t_1237, t_1238, pc_y, lsi0_766, lsi0_767, lsi0_768, lsi1_766, \
                         lsi1_767, lsi1_768, lsk_987, lsk_988, \
                         lsk_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1236[k] = f_12 * lsi0_766[k]
                    - f_13 * lsi1_766[k]
                    + f_3 * pc_y[k] * lsk_987[k];

        t_1237[k] = f_10 * lsi0_767[k]
                    - f_11 * lsi1_767[k]
                    + f_3 * pc_y[k] * lsk_988[k];

        t_1238[k] = f_8 * lsi0_768[k]
                    - f_9 * lsi1_768[k]
                    + f_3 * pc_y[k] * lsk_989[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, pc_y, lsi0_769, lsi0_770, lsi1_769, lsi1_770, \
                         lsk_990, lsk_991, lsk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_6 * lsi0_769[k]
                    - f_7 * lsi1_769[k]
                    + f_3 * pc_y[k] * lsk_990[k];

        t_1240[k] = f_4 * lsi0_770[k]
                    - f_5 * lsi1_770[k]
                    + f_3 * pc_y[k] * lsk_991[k];

        t_1241[k] = f_3 * pc_y[k] * lsk_992[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, pc_x, ksk_999, ksk_1000, ksk_1001, \
                         ksk_1002, lsi0_783, lsi1_783, lsk_999, lsk_1000, lsk_1001, \
                         lsk_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = f_16 * ksk_999[k]
                    + f_4 * lsi0_783[k]
                    - f_5 * lsi1_783[k]
                    + f_3 * pc_x[k] * lsk_999[k];

        t_1243[k] = f_16 * ksk_1000[k]
                    + f_3 * pc_x[k] * lsk_1000[k];

        t_1244[k] = f_16 * ksk_1001[k]
                    + f_3 * pc_x[k] * lsk_1001[k];

        t_1245[k] = f_16 * ksk_1002[k]
                    + f_3 * pc_x[k] * lsk_1002[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, t_1250, pc_x, pc_y, ksk_1003, \
                         ksk_1004, ksk_1005, ksk_1007, lsk_999, lsk_1003, lsk_1004, lsk_1005, \
                         lsk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = f_16 * ksk_1003[k]
                    + f_3 * pc_x[k] * lsk_1003[k];

        t_1247[k] = f_16 * ksk_1004[k]
                    + f_3 * pc_x[k] * lsk_1004[k];

        t_1248[k] = f_16 * ksk_1005[k]
                    + f_3 * pc_x[k] * lsk_1005[k];

        t_1249[k] = f_3 * pc_y[k] * lsk_999[k];

        t_1250[k] = f_16 * ksk_1007[k]
                    + f_3 * pc_x[k] * lsk_1007[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksl0,
                                                           const size_t ksk, const size_t ksl1,
                                                           const size_t lsi0, const size_t lsi1,
                                                           const size_t lsk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_945 = buffer.data(ksl0 + 945);
    const auto *ksl0_948 = buffer.data(ksl0 + 948);
    const auto *ksl0_951 = buffer.data(ksl0 + 951);
    const auto *ksl0_955 = buffer.data(ksl0 + 955);
    const auto *ksl0_960 = buffer.data(ksl0 + 960);
    const auto *ksl0_966 = buffer.data(ksl0 + 966);
    const auto *ksl0_1260 = buffer.data(ksl0 + 1260);
    const auto *ksl0_1263 = buffer.data(ksl0 + 1263);
    const auto *ksl0_1266 = buffer.data(ksl0 + 1266);
    const auto *ksl0_1270 = buffer.data(ksl0 + 1270);
    const auto *ksl0_1275 = buffer.data(ksl0 + 1275);
    const auto *ksl0_1281 = buffer.data(ksl0 + 1281);
    const auto *ksl0_1296 = buffer.data(ksl0 + 1296);
    const auto *ksl0_1298 = buffer.data(ksl0 + 1298);
    const auto *ksl0_1299 = buffer.data(ksl0 + 1299);
    const auto *ksl0_1300 = buffer.data(ksl0 + 1300);
    const auto *ksl0_1301 = buffer.data(ksl0 + 1301);
    const auto *ksl0_1302 = buffer.data(ksl0 + 1302);
    const auto *ksl0_1304 = buffer.data(ksl0 + 1304);
    const auto *ksl0_1310 = buffer.data(ksl0 + 1310);
    const auto *ksl0_1314 = buffer.data(ksl0 + 1314);
    const auto *ksl0_1317 = buffer.data(ksl0 + 1317);
    const auto *ksl0_1319 = buffer.data(ksl0 + 1319);
    const auto *ksl0_1322 = buffer.data(ksl0 + 1322);
    const auto *ksl0_1323 = buffer.data(ksl0 + 1323);
    const auto *ksl0_1325 = buffer.data(ksl0 + 1325);
    const auto *ksl0_1328 = buffer.data(ksl0 + 1328);
    const auto *ksl0_1329 = buffer.data(ksl0 + 1329);
    const auto *ksl0_1330 = buffer.data(ksl0 + 1330);
    const auto *ksl0_1332 = buffer.data(ksl0 + 1332);
    const auto *ksl0_1341 = buffer.data(ksl0 + 1341);
    const auto *ksl0_1343 = buffer.data(ksl0 + 1343);
    const auto *ksl0_1344 = buffer.data(ksl0 + 1344);
    const auto *ksl0_1345 = buffer.data(ksl0 + 1345);
    const auto *ksl0_1346 = buffer.data(ksl0 + 1346);
    const auto *ksl0_1347 = buffer.data(ksl0 + 1347);
    const auto *ksl0_1349 = buffer.data(ksl0 + 1349);
    const auto *ksl0_1350 = buffer.data(ksl0 + 1350);
    const auto *ksl0_1353 = buffer.data(ksl0 + 1353);
    const auto *ksl0_1355 = buffer.data(ksl0 + 1355);
    const auto *ksl0_1356 = buffer.data(ksl0 + 1356);
    const auto *ksl0_1359 = buffer.data(ksl0 + 1359);
    const auto *ksl0_1360 = buffer.data(ksl0 + 1360);
    const auto *ksl0_1362 = buffer.data(ksl0 + 1362);
    const auto *ksl0_1364 = buffer.data(ksl0 + 1364);
    const auto *ksl0_1365 = buffer.data(ksl0 + 1365);
    const auto *ksl0_1367 = buffer.data(ksl0 + 1367);

    const auto *ksk_755 = buffer.data(ksk + 755);
    const auto *ksk_756 = buffer.data(ksk + 756);
    const auto *ksk_759 = buffer.data(ksk + 759);
    const auto *ksk_761 = buffer.data(ksk + 761);
    const auto *ksk_762 = buffer.data(ksk + 762);
    const auto *ksk_765 = buffer.data(ksk + 765);
    const auto *ksk_766 = buffer.data(ksk + 766);
    const auto *ksk_770 = buffer.data(ksk + 770);
    const auto *ksk_771 = buffer.data(ksk + 771);
    const auto *ksk_776 = buffer.data(ksk + 776);
    const auto *ksk_784 = buffer.data(ksk + 784);
    const auto *ksk_791 = buffer.data(ksk + 791);
    const auto *ksk_792 = buffer.data(ksk + 792);
    const auto *ksk_794 = buffer.data(ksk + 794);
    const auto *ksk_795 = buffer.data(ksk + 795);
    const auto *ksk_797 = buffer.data(ksk + 797);
    const auto *ksk_798 = buffer.data(ksk + 798);
    const auto *ksk_801 = buffer.data(ksk + 801);
    const auto *ksk_802 = buffer.data(ksk + 802);
    const auto *ksk_806 = buffer.data(ksk + 806);
    const auto *ksk_812 = buffer.data(ksk + 812);
    const auto *ksk_827 = buffer.data(ksk + 827);
    const auto *ksk_828 = buffer.data(ksk + 828);
    const auto *ksk_830 = buffer.data(ksk + 830);
    const auto *ksk_833 = buffer.data(ksk + 833);
    const auto *ksk_837 = buffer.data(ksk + 837);
    const auto *ksk_1008 = buffer.data(ksk + 1008);
    const auto *ksk_1011 = buffer.data(ksk + 1011);
    const auto *ksk_1014 = buffer.data(ksk + 1014);
    const auto *ksk_1018 = buffer.data(ksk + 1018);
    const auto *ksk_1023 = buffer.data(ksk + 1023);
    const auto *ksk_1029 = buffer.data(ksk + 1029);
    const auto *ksk_1036 = buffer.data(ksk + 1036);
    const auto *ksk_1038 = buffer.data(ksk + 1038);
    const auto *ksk_1039 = buffer.data(ksk + 1039);
    const auto *ksk_1040 = buffer.data(ksk + 1040);
    const auto *ksk_1041 = buffer.data(ksk + 1041);
    const auto *ksk_1042 = buffer.data(ksk + 1042);
    const auto *ksk_1043 = buffer.data(ksk + 1043);
    const auto *ksk_1049 = buffer.data(ksk + 1049);
    const auto *ksk_1053 = buffer.data(ksk + 1053);
    const auto *ksk_1056 = buffer.data(ksk + 1056);
    const auto *ksk_1058 = buffer.data(ksk + 1058);
    const auto *ksk_1061 = buffer.data(ksk + 1061);
    const auto *ksk_1062 = buffer.data(ksk + 1062);
    const auto *ksk_1064 = buffer.data(ksk + 1064);
    const auto *ksk_1067 = buffer.data(ksk + 1067);
    const auto *ksk_1068 = buffer.data(ksk + 1068);
    const auto *ksk_1069 = buffer.data(ksk + 1069);
    const auto *ksk_1071 = buffer.data(ksk + 1071);
    const auto *ksk_1072 = buffer.data(ksk + 1072);
    const auto *ksk_1073 = buffer.data(ksk + 1073);
    const auto *ksk_1074 = buffer.data(ksk + 1074);
    const auto *ksk_1075 = buffer.data(ksk + 1075);
    const auto *ksk_1076 = buffer.data(ksk + 1076);
    const auto *ksk_1077 = buffer.data(ksk + 1077);
    const auto *ksk_1078 = buffer.data(ksk + 1078);
    const auto *ksk_1079 = buffer.data(ksk + 1079);
    const auto *ksk_1080 = buffer.data(ksk + 1080);
    const auto *ksk_1083 = buffer.data(ksk + 1083);
    const auto *ksk_1085 = buffer.data(ksk + 1085);
    const auto *ksk_1086 = buffer.data(ksk + 1086);
    const auto *ksk_1089 = buffer.data(ksk + 1089);
    const auto *ksk_1090 = buffer.data(ksk + 1090);
    const auto *ksk_1092 = buffer.data(ksk + 1092);
    const auto *ksk_1094 = buffer.data(ksk + 1094);
    const auto *ksk_1095 = buffer.data(ksk + 1095);
    const auto *ksk_1097 = buffer.data(ksk + 1097);

    const auto *ksl1_945 = buffer.data(ksl1 + 945);
    const auto *ksl1_948 = buffer.data(ksl1 + 948);
    const auto *ksl1_951 = buffer.data(ksl1 + 951);
    const auto *ksl1_955 = buffer.data(ksl1 + 955);
    const auto *ksl1_960 = buffer.data(ksl1 + 960);
    const auto *ksl1_966 = buffer.data(ksl1 + 966);
    const auto *ksl1_1260 = buffer.data(ksl1 + 1260);
    const auto *ksl1_1263 = buffer.data(ksl1 + 1263);
    const auto *ksl1_1266 = buffer.data(ksl1 + 1266);
    const auto *ksl1_1270 = buffer.data(ksl1 + 1270);
    const auto *ksl1_1275 = buffer.data(ksl1 + 1275);
    const auto *ksl1_1281 = buffer.data(ksl1 + 1281);
    const auto *ksl1_1296 = buffer.data(ksl1 + 1296);
    const auto *ksl1_1298 = buffer.data(ksl1 + 1298);
    const auto *ksl1_1299 = buffer.data(ksl1 + 1299);
    const auto *ksl1_1300 = buffer.data(ksl1 + 1300);
    const auto *ksl1_1301 = buffer.data(ksl1 + 1301);
    const auto *ksl1_1302 = buffer.data(ksl1 + 1302);
    const auto *ksl1_1304 = buffer.data(ksl1 + 1304);
    const auto *ksl1_1310 = buffer.data(ksl1 + 1310);
    const auto *ksl1_1314 = buffer.data(ksl1 + 1314);
    const auto *ksl1_1317 = buffer.data(ksl1 + 1317);
    const auto *ksl1_1319 = buffer.data(ksl1 + 1319);
    const auto *ksl1_1322 = buffer.data(ksl1 + 1322);
    const auto *ksl1_1323 = buffer.data(ksl1 + 1323);
    const auto *ksl1_1325 = buffer.data(ksl1 + 1325);
    const auto *ksl1_1328 = buffer.data(ksl1 + 1328);
    const auto *ksl1_1329 = buffer.data(ksl1 + 1329);
    const auto *ksl1_1330 = buffer.data(ksl1 + 1330);
    const auto *ksl1_1332 = buffer.data(ksl1 + 1332);
    const auto *ksl1_1341 = buffer.data(ksl1 + 1341);
    const auto *ksl1_1343 = buffer.data(ksl1 + 1343);
    const auto *ksl1_1344 = buffer.data(ksl1 + 1344);
    const auto *ksl1_1345 = buffer.data(ksl1 + 1345);
    const auto *ksl1_1346 = buffer.data(ksl1 + 1346);
    const auto *ksl1_1347 = buffer.data(ksl1 + 1347);
    const auto *ksl1_1349 = buffer.data(ksl1 + 1349);
    const auto *ksl1_1350 = buffer.data(ksl1 + 1350);
    const auto *ksl1_1353 = buffer.data(ksl1 + 1353);
    const auto *ksl1_1355 = buffer.data(ksl1 + 1355);
    const auto *ksl1_1356 = buffer.data(ksl1 + 1356);
    const auto *ksl1_1359 = buffer.data(ksl1 + 1359);
    const auto *ksl1_1360 = buffer.data(ksl1 + 1360);
    const auto *ksl1_1362 = buffer.data(ksl1 + 1362);
    const auto *ksl1_1364 = buffer.data(ksl1 + 1364);
    const auto *ksl1_1365 = buffer.data(ksl1 + 1365);
    const auto *ksl1_1367 = buffer.data(ksl1 + 1367);

    const auto *lsi0_777 = buffer.data(lsi0 + 777);
    const auto *lsi0_778 = buffer.data(lsi0 + 778);
    const auto *lsi0_779 = buffer.data(lsi0 + 779);
    const auto *lsi0_780 = buffer.data(lsi0 + 780);
    const auto *lsi0_781 = buffer.data(lsi0 + 781);
    const auto *lsi0_782 = buffer.data(lsi0 + 782);
    const auto *lsi0_783 = buffer.data(lsi0 + 783);
    const auto *lsi0_784 = buffer.data(lsi0 + 784);
    const auto *lsi0_786 = buffer.data(lsi0 + 786);
    const auto *lsi0_787 = buffer.data(lsi0 + 787);
    const auto *lsi0_789 = buffer.data(lsi0 + 789);
    const auto *lsi0_790 = buffer.data(lsi0 + 790);
    const auto *lsi0_791 = buffer.data(lsi0 + 791);
    const auto *lsi0_793 = buffer.data(lsi0 + 793);
    const auto *lsi0_794 = buffer.data(lsi0 + 794);
    const auto *lsi0_795 = buffer.data(lsi0 + 795);
    const auto *lsi0_796 = buffer.data(lsi0 + 796);
    const auto *lsi0_798 = buffer.data(lsi0 + 798);

    const auto *lsi1_777 = buffer.data(lsi1 + 777);
    const auto *lsi1_778 = buffer.data(lsi1 + 778);
    const auto *lsi1_779 = buffer.data(lsi1 + 779);
    const auto *lsi1_780 = buffer.data(lsi1 + 780);
    const auto *lsi1_781 = buffer.data(lsi1 + 781);
    const auto *lsi1_782 = buffer.data(lsi1 + 782);
    const auto *lsi1_783 = buffer.data(lsi1 + 783);
    const auto *lsi1_784 = buffer.data(lsi1 + 784);
    const auto *lsi1_786 = buffer.data(lsi1 + 786);
    const auto *lsi1_787 = buffer.data(lsi1 + 787);
    const auto *lsi1_789 = buffer.data(lsi1 + 789);
    const auto *lsi1_790 = buffer.data(lsi1 + 790);
    const auto *lsi1_791 = buffer.data(lsi1 + 791);
    const auto *lsi1_793 = buffer.data(lsi1 + 793);
    const auto *lsi1_794 = buffer.data(lsi1 + 794);
    const auto *lsi1_795 = buffer.data(lsi1 + 795);
    const auto *lsi1_796 = buffer.data(lsi1 + 796);
    const auto *lsi1_798 = buffer.data(lsi1 + 798);

    const auto *lsk_1000 = buffer.data(lsk + 1000);
    const auto *lsk_1001 = buffer.data(lsk + 1001);
    const auto *lsk_1002 = buffer.data(lsk + 1002);
    const auto *lsk_1003 = buffer.data(lsk + 1003);
    const auto *lsk_1004 = buffer.data(lsk + 1004);
    const auto *lsk_1005 = buffer.data(lsk + 1005);
    const auto *lsk_1006 = buffer.data(lsk + 1006);
    const auto *lsk_1007 = buffer.data(lsk + 1007);
    const auto *lsk_1008 = buffer.data(lsk + 1008);
    const auto *lsk_1009 = buffer.data(lsk + 1009);
    const auto *lsk_1010 = buffer.data(lsk + 1010);
    const auto *lsk_1011 = buffer.data(lsk + 1011);
    const auto *lsk_1013 = buffer.data(lsk + 1013);
    const auto *lsk_1014 = buffer.data(lsk + 1014);
    const auto *lsk_1015 = buffer.data(lsk + 1015);
    const auto *lsk_1017 = buffer.data(lsk + 1017);
    const auto *lsk_1018 = buffer.data(lsk + 1018);
    const auto *lsk_1019 = buffer.data(lsk + 1019);
    const auto *lsk_1020 = buffer.data(lsk + 1020);
    const auto *lsk_1022 = buffer.data(lsk + 1022);
    const auto *lsk_1023 = buffer.data(lsk + 1023);
    const auto *lsk_1024 = buffer.data(lsk + 1024);
    const auto *lsk_1025 = buffer.data(lsk + 1025);
    const auto *lsk_1026 = buffer.data(lsk + 1026);
    const auto *lsk_1028 = buffer.data(lsk + 1028);
    const auto *lsk_1029 = buffer.data(lsk + 1029);
    const auto *lsk_1036 = buffer.data(lsk + 1036);
    const auto *lsk_1038 = buffer.data(lsk + 1038);
    const auto *lsk_1039 = buffer.data(lsk + 1039);
    const auto *lsk_1040 = buffer.data(lsk + 1040);
    const auto *lsk_1041 = buffer.data(lsk + 1041);
    const auto *lsk_1042 = buffer.data(lsk + 1042);
    const auto *lsk_1043 = buffer.data(lsk + 1043);
    const auto *lsk_1044 = buffer.data(lsk + 1044);
    const auto *lsk_1046 = buffer.data(lsk + 1046);
    const auto *lsk_1047 = buffer.data(lsk + 1047);
    const auto *lsk_1049 = buffer.data(lsk + 1049);
    const auto *lsk_1050 = buffer.data(lsk + 1050);
    const auto *lsk_1053 = buffer.data(lsk + 1053);
    const auto *lsk_1054 = buffer.data(lsk + 1054);
    const auto *lsk_1058 = buffer.data(lsk + 1058);
    const auto *lsk_1059 = buffer.data(lsk + 1059);
    const auto *lsk_1064 = buffer.data(lsk + 1064);
    const auto *lsk_1072 = buffer.data(lsk + 1072);
    const auto *lsk_1073 = buffer.data(lsk + 1073);
    const auto *lsk_1074 = buffer.data(lsk + 1074);
    const auto *lsk_1075 = buffer.data(lsk + 1075);
    const auto *lsk_1076 = buffer.data(lsk + 1076);
    const auto *lsk_1077 = buffer.data(lsk + 1077);
    const auto *lsk_1078 = buffer.data(lsk + 1078);
    const auto *lsk_1079 = buffer.data(lsk + 1079);
    const auto *lsk_1080 = buffer.data(lsk + 1080);
    const auto *lsk_1082 = buffer.data(lsk + 1082);
    const auto *lsk_1083 = buffer.data(lsk + 1083);
    const auto *lsk_1085 = buffer.data(lsk + 1085);
    const auto *lsk_1086 = buffer.data(lsk + 1086);
    const auto *lsk_1089 = buffer.data(lsk + 1089);
    const auto *lsk_1090 = buffer.data(lsk + 1090);

#pragma omp simd aligned(t_1251, t_1252, t_1253, pc_y, lsi0_777, lsi0_778, lsi0_779, lsi1_777, \
                         lsi1_778, lsi1_779, lsk_1000, lsk_1001, \
                         lsk_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_1 * lsi0_777[k]
                    - f_2 * lsi1_777[k]
                    + f_3 * pc_y[k] * lsk_1000[k];

        t_1252[k] = f_22 * lsi0_778[k]
                    - f_23 * lsi1_778[k]
                    + f_3 * pc_y[k] * lsk_1001[k];

        t_1253[k] = f_12 * lsi0_779[k]
                    - f_13 * lsi1_779[k]
                    + f_3 * pc_y[k] * lsk_1002[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, lsi0_780, lsi0_781, lsi0_782, lsi1_780, \
                         lsi1_781, lsi1_782, lsk_1003, lsk_1004, \
                         lsk_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_10 * lsi0_780[k]
                    - f_11 * lsi1_780[k]
                    + f_3 * pc_y[k] * lsk_1003[k];

        t_1255[k] = f_8 * lsi0_781[k]
                    - f_9 * lsi1_781[k]
                    + f_3 * pc_y[k] * lsk_1004[k];

        t_1256[k] = f_6 * lsi0_782[k]
                    - f_7 * lsi1_782[k]
                    + f_3 * pc_y[k] * lsk_1005[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, t_1260, pa_x, pc_x, pc_y, pc_z, ksl0_1260, \
                         ksk_755, ksk_1008, ksl1_1260, lsi0_783, lsi1_783, lsk_1006, \
                         lsk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_4 * lsi0_783[k]
                    - f_5 * lsi1_783[k]
                    + f_3 * pc_y[k] * lsk_1006[k];

        t_1258[k] = f_3 * pc_y[k] * lsk_1007[k];

        t_1259[k] = f_20 * ksk_755[k]
                    + f_1 * lsi0_783[k]
                    - f_2 * lsi1_783[k]
                    + f_3 * pc_z[k] * lsk_1007[k];

        t_1260[k] = pa_x[k] * ksl0_1260[k]
                    + f_0 * ksk_1008[k]
                    - f_14 * pc_x[k] * ksl1_1260[k];
    }

#pragma omp simd aligned(t_1261, t_1262, t_1263, t_1264, pa_x, pc_x, pc_y, pc_z, ksl0_1263, \
                         ksk_756, ksk_1011, ksl1_1263, lsk_1008, \
                         lsk_1009 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1261[k] = f_21 * ksk_756[k]
                    + f_3 * pc_y[k] * lsk_1008[k];

        t_1262[k] = f_3 * pc_z[k] * lsk_1008[k];

        t_1263[k] = pa_x[k] * ksl0_1263[k]
                    + f_20 * ksk_1011[k]
                    - f_14 * pc_x[k] * ksl1_1263[k];

        t_1264[k] = f_3 * pc_z[k] * lsk_1009[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, pa_x, pc_x, pc_z, ksl0_1266, ksk_1014, \
                         ksl1_1266, lsi0_784, lsi1_784, lsk_1010, \
                         lsk_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_4 * lsi0_784[k]
                    - f_5 * lsi1_784[k]
                    + f_3 * pc_z[k] * lsk_1010[k];

        t_1266[k] = pa_x[k] * ksl0_1266[k]
                    + f_19 * ksk_1014[k]
                    - f_14 * pc_x[k] * ksl1_1266[k];

        t_1267[k] = f_3 * pc_z[k] * lsk_1011[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, t_1271, pa_x, pc_x, pc_y, pc_z, ksl0_1270, \
                         ksk_761, ksk_1018, ksl1_1270, lsi0_786, lsi1_786, lsk_1013, \
                         lsk_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_21 * ksk_761[k]
                    + f_3 * pc_y[k] * lsk_1013[k];

        t_1269[k] = f_6 * lsi0_786[k]
                    - f_7 * lsi1_786[k]
                    + f_3 * pc_z[k] * lsk_1013[k];

        t_1270[k] = pa_x[k] * ksl0_1270[k]
                    + f_18 * ksk_1018[k]
                    - f_14 * pc_x[k] * ksl1_1270[k];

        t_1271[k] = f_3 * pc_z[k] * lsk_1014[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, pc_y, pc_z, ksk_765, lsi0_787, lsi0_789, \
                         lsi1_787, lsi1_789, lsk_1015, lsk_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_4 * lsi0_787[k]
                    - f_5 * lsi1_787[k]
                    + f_3 * pc_z[k] * lsk_1015[k];

        t_1273[k] = f_21 * ksk_765[k]
                    + f_3 * pc_y[k] * lsk_1017[k];

        t_1274[k] = f_8 * lsi0_789[k]
                    - f_9 * lsi1_789[k]
                    + f_3 * pc_z[k] * lsk_1017[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, pa_x, pc_x, pc_z, ksl0_1275, ksk_1023, \
                         ksl1_1275, lsi0_790, lsi1_790, lsk_1018, \
                         lsk_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = pa_x[k] * ksl0_1275[k]
                    + f_17 * ksk_1023[k]
                    - f_14 * pc_x[k] * ksl1_1275[k];

        t_1276[k] = f_3 * pc_z[k] * lsk_1018[k];

        t_1277[k] = f_4 * lsi0_790[k]
                    - f_5 * lsi1_790[k]
                    + f_3 * pc_z[k] * lsk_1019[k];
    }

#pragma omp simd aligned(t_1278, t_1279, t_1280, pc_y, pc_z, ksk_770, lsi0_791, lsi0_793, \
                         lsi1_791, lsi1_793, lsk_1020, lsk_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1278[k] = f_6 * lsi0_791[k]
                    - f_7 * lsi1_791[k]
                    + f_3 * pc_z[k] * lsk_1020[k];

        t_1279[k] = f_21 * ksk_770[k]
                    + f_3 * pc_y[k] * lsk_1022[k];

        t_1280[k] = f_10 * lsi0_793[k]
                    - f_11 * lsi1_793[k]
                    + f_3 * pc_z[k] * lsk_1022[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, pa_x, pc_x, pc_z, ksl0_1281, ksk_1029, \
                         ksl1_1281, lsi0_794, lsi1_794, lsk_1023, \
                         lsk_1024 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = pa_x[k] * ksl0_1281[k]
                    + f_16 * ksk_1029[k]
                    - f_14 * pc_x[k] * ksl1_1281[k];

        t_1282[k] = f_3 * pc_z[k] * lsk_1023[k];

        t_1283[k] = f_4 * lsi0_794[k]
                    - f_5 * lsi1_794[k]
                    + f_3 * pc_z[k] * lsk_1024[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, t_1287, pc_y, pc_z, ksk_776, lsi0_795, \
                         lsi0_796, lsi0_798, lsi1_795, lsi1_796, lsi1_798, lsk_1025, lsk_1026, \
                         lsk_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_6 * lsi0_795[k]
                    - f_7 * lsi1_795[k]
                    + f_3 * pc_z[k] * lsk_1025[k];

        t_1285[k] = f_8 * lsi0_796[k]
                    - f_9 * lsi1_796[k]
                    + f_3 * pc_z[k] * lsk_1026[k];

        t_1286[k] = f_21 * ksk_776[k]
                    + f_3 * pc_y[k] * lsk_1028[k];

        t_1287[k] = f_12 * lsi0_798[k]
                    - f_13 * lsi1_798[k]
                    + f_3 * pc_z[k] * lsk_1028[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, t_1292, pc_x, pc_z, ksk_1036, \
                         ksk_1038, ksk_1039, ksk_1040, lsk_1029, lsk_1036, lsk_1038, lsk_1039, \
                         lsk_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = f_15 * ksk_1036[k]
                    + f_3 * pc_x[k] * lsk_1036[k];

        t_1289[k] = f_3 * pc_z[k] * lsk_1029[k];

        t_1290[k] = f_15 * ksk_1038[k]
                    + f_3 * pc_x[k] * lsk_1038[k];

        t_1291[k] = f_15 * ksk_1039[k]
                    + f_3 * pc_x[k] * lsk_1039[k];

        t_1292[k] = f_15 * ksk_1040[k]
                    + f_3 * pc_x[k] * lsk_1040[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, t_1296, pa_x, pc_x, ksl0_1296, ksk_1041, \
                         ksk_1042, ksk_1043, ksl1_1296, lsk_1041, lsk_1042, \
                         lsk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_15 * ksk_1041[k]
                    + f_3 * pc_x[k] * lsk_1041[k];

        t_1294[k] = f_15 * ksk_1042[k]
                    + f_3 * pc_x[k] * lsk_1042[k];

        t_1295[k] = f_15 * ksk_1043[k]
                    + f_3 * pc_x[k] * lsk_1043[k];

        t_1296[k] = pa_x[k] * ksl0_1296[k]
                    - f_14 * pc_x[k] * ksl1_1296[k];
    }

#pragma omp simd aligned(t_1297, t_1298, t_1299, t_1300, pa_x, pc_x, pc_z, ksl0_1298, \
                         ksl0_1299, ksl0_1300, ksl1_1298, ksl1_1299, ksl1_1300, \
                         lsk_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1297[k] = f_3 * pc_z[k] * lsk_1036[k];

        t_1298[k] = pa_x[k] * ksl0_1298[k]
                    - f_14 * pc_x[k] * ksl1_1298[k];

        t_1299[k] = pa_x[k] * ksl0_1299[k]
                    - f_14 * pc_x[k] * ksl1_1299[k];

        t_1300[k] = pa_x[k] * ksl0_1300[k]
                    - f_14 * pc_x[k] * ksl1_1300[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, pa_x, pc_x, pc_y, ksl0_1301, \
                         ksl0_1302, ksl0_1304, ksk_791, ksl1_1301, ksl1_1302, ksl1_1304, \
                         lsk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = pa_x[k] * ksl0_1301[k]
                    - f_14 * pc_x[k] * ksl1_1301[k];

        t_1302[k] = pa_x[k] * ksl0_1302[k]
                    - f_14 * pc_x[k] * ksl1_1302[k];

        t_1303[k] = f_21 * ksk_791[k]
                    + f_3 * pc_y[k] * lsk_1043[k];

        t_1304[k] = pa_x[k] * ksl0_1304[k]
                    - f_14 * pc_x[k] * ksl1_1304[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, t_1308, pa_z, pc_y, pc_z, ksl0_945, ksl0_948, \
                         ksk_756, ksk_792, ksl1_945, ksl1_948, \
                         lsk_1044 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = pa_z[k] * ksl0_945[k]
                    - f_14 * pc_z[k] * ksl1_945[k];

        t_1306[k] = f_20 * ksk_792[k]
                    + f_3 * pc_y[k] * lsk_1044[k];

        t_1307[k] = f_15 * ksk_756[k]
                    + f_3 * pc_z[k] * lsk_1044[k];

        t_1308[k] = pa_z[k] * ksl0_948[k]
                    - f_14 * pc_z[k] * ksl1_948[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, pa_x, pa_z, pc_x, pc_y, pc_z, ksl0_951, \
                         ksl0_1310, ksk_794, ksk_1049, ksl1_951, ksl1_1310, \
                         lsk_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = f_20 * ksk_794[k]
                    + f_3 * pc_y[k] * lsk_1046[k];

        t_1310[k] = pa_x[k] * ksl0_1310[k]
                    + f_20 * ksk_1049[k]
                    - f_14 * pc_x[k] * ksl1_1310[k];

        t_1311[k] = pa_z[k] * ksl0_951[k]
                    - f_14 * pc_z[k] * ksl1_951[k];
    }

#pragma omp simd aligned(t_1312, t_1313, t_1314, pa_x, pc_x, pc_y, pc_z, ksl0_1314, ksk_759, \
                         ksk_797, ksk_1053, ksl1_1314, lsk_1047, \
                         lsk_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1312[k] = f_15 * ksk_759[k]
                    + f_3 * pc_z[k] * lsk_1047[k];

        t_1313[k] = f_20 * ksk_797[k]
                    + f_3 * pc_y[k] * lsk_1049[k];

        t_1314[k] = pa_x[k] * ksl0_1314[k]
                    + f_19 * ksk_1053[k]
                    - f_14 * pc_x[k] * ksl1_1314[k];
    }

#pragma omp simd aligned(t_1315, t_1316, t_1317, pa_x, pa_z, pc_x, pc_z, ksl0_955, ksl0_1317, \
                         ksk_762, ksk_1056, ksl1_955, ksl1_1317, \
                         lsk_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1315[k] = pa_z[k] * ksl0_955[k]
                    - f_14 * pc_z[k] * ksl1_955[k];

        t_1316[k] = f_15 * ksk_762[k]
                    + f_3 * pc_z[k] * lsk_1050[k];

        t_1317[k] = pa_x[k] * ksl0_1317[k]
                    + f_18 * ksk_1056[k]
                    - f_14 * pc_x[k] * ksl1_1317[k];
    }

#pragma omp simd aligned(t_1318, t_1319, t_1320, pa_x, pa_z, pc_x, pc_y, pc_z, ksl0_960, \
                         ksl0_1319, ksk_801, ksk_1058, ksl1_960, ksl1_1319, \
                         lsk_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1318[k] = f_20 * ksk_801[k]
                    + f_3 * pc_y[k] * lsk_1053[k];

        t_1319[k] = pa_x[k] * ksl0_1319[k]
                    + f_18 * ksk_1058[k]
                    - f_14 * pc_x[k] * ksl1_1319[k];

        t_1320[k] = pa_z[k] * ksl0_960[k]
                    - f_14 * pc_z[k] * ksl1_960[k];
    }

#pragma omp simd aligned(t_1321, t_1322, t_1323, pa_x, pc_x, pc_z, ksl0_1322, ksl0_1323, \
                         ksk_766, ksk_1061, ksk_1062, ksl1_1322, ksl1_1323, \
                         lsk_1054 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1321[k] = f_15 * ksk_766[k]
                    + f_3 * pc_z[k] * lsk_1054[k];

        t_1322[k] = pa_x[k] * ksl0_1322[k]
                    + f_17 * ksk_1061[k]
                    - f_14 * pc_x[k] * ksl1_1322[k];

        t_1323[k] = pa_x[k] * ksl0_1323[k]
                    + f_17 * ksk_1062[k]
                    - f_14 * pc_x[k] * ksl1_1323[k];
    }

#pragma omp simd aligned(t_1324, t_1325, t_1326, pa_x, pa_z, pc_x, pc_y, pc_z, ksl0_966, \
                         ksl0_1325, ksk_806, ksk_1064, ksl1_966, ksl1_1325, \
                         lsk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1324[k] = f_20 * ksk_806[k]
                    + f_3 * pc_y[k] * lsk_1058[k];

        t_1325[k] = pa_x[k] * ksl0_1325[k]
                    + f_17 * ksk_1064[k]
                    - f_14 * pc_x[k] * ksl1_1325[k];

        t_1326[k] = pa_z[k] * ksl0_966[k]
                    - f_14 * pc_z[k] * ksl1_966[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, pa_x, pc_x, pc_z, ksl0_1328, ksl0_1329, \
                         ksk_771, ksk_1067, ksk_1068, ksl1_1328, ksl1_1329, \
                         lsk_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = f_15 * ksk_771[k]
                    + f_3 * pc_z[k] * lsk_1059[k];

        t_1328[k] = pa_x[k] * ksl0_1328[k]
                    + f_16 * ksk_1067[k]
                    - f_14 * pc_x[k] * ksl1_1328[k];

        t_1329[k] = pa_x[k] * ksl0_1329[k]
                    + f_16 * ksk_1068[k]
                    - f_14 * pc_x[k] * ksl1_1329[k];
    }

#pragma omp simd aligned(t_1330, t_1331, t_1332, pa_x, pc_x, pc_y, ksl0_1330, ksl0_1332, \
                         ksk_812, ksk_1069, ksk_1071, ksl1_1330, ksl1_1332, \
                         lsk_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1330[k] = pa_x[k] * ksl0_1330[k]
                    + f_16 * ksk_1069[k]
                    - f_14 * pc_x[k] * ksl1_1330[k];

        t_1331[k] = f_20 * ksk_812[k]
                    + f_3 * pc_y[k] * lsk_1064[k];

        t_1332[k] = pa_x[k] * ksl0_1332[k]
                    + f_16 * ksk_1071[k]
                    - f_14 * pc_x[k] * ksl1_1332[k];
    }

#pragma omp simd aligned(t_1333, t_1334, t_1335, t_1336, t_1337, pc_x, ksk_1072, ksk_1073, \
                         ksk_1074, ksk_1075, ksk_1076, lsk_1072, lsk_1073, lsk_1074, lsk_1075, \
                         lsk_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1333[k] = f_15 * ksk_1072[k]
                    + f_3 * pc_x[k] * lsk_1072[k];

        t_1334[k] = f_15 * ksk_1073[k]
                    + f_3 * pc_x[k] * lsk_1073[k];

        t_1335[k] = f_15 * ksk_1074[k]
                    + f_3 * pc_x[k] * lsk_1074[k];

        t_1336[k] = f_15 * ksk_1075[k]
                    + f_3 * pc_x[k] * lsk_1075[k];

        t_1337[k] = f_15 * ksk_1076[k]
                    + f_3 * pc_x[k] * lsk_1076[k];
    }

#pragma omp simd aligned(t_1338, t_1339, t_1340, t_1341, pa_x, pc_x, ksl0_1341, ksk_1077, \
                         ksk_1078, ksk_1079, ksl1_1341, lsk_1077, lsk_1078, \
                         lsk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1338[k] = f_15 * ksk_1077[k]
                    + f_3 * pc_x[k] * lsk_1077[k];

        t_1339[k] = f_15 * ksk_1078[k]
                    + f_3 * pc_x[k] * lsk_1078[k];

        t_1340[k] = f_15 * ksk_1079[k]
                    + f_3 * pc_x[k] * lsk_1079[k];

        t_1341[k] = pa_x[k] * ksl0_1341[k]
                    - f_14 * pc_x[k] * ksl1_1341[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, t_1345, pa_x, pc_x, pc_z, ksl0_1343, \
                         ksl0_1344, ksl0_1345, ksk_784, ksl1_1343, ksl1_1344, ksl1_1345, \
                         lsk_1072 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_15 * ksk_784[k]
                    + f_3 * pc_z[k] * lsk_1072[k];

        t_1343[k] = pa_x[k] * ksl0_1343[k]
                    - f_14 * pc_x[k] * ksl1_1343[k];

        t_1344[k] = pa_x[k] * ksl0_1344[k]
                    - f_14 * pc_x[k] * ksl1_1344[k];

        t_1345[k] = pa_x[k] * ksl0_1345[k]
                    - f_14 * pc_x[k] * ksl1_1345[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, t_1349, pa_x, pc_x, pc_y, ksl0_1346, \
                         ksl0_1347, ksl0_1349, ksk_827, ksl1_1346, ksl1_1347, ksl1_1349, \
                         lsk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = pa_x[k] * ksl0_1346[k]
                    - f_14 * pc_x[k] * ksl1_1346[k];

        t_1347[k] = pa_x[k] * ksl0_1347[k]
                    - f_14 * pc_x[k] * ksl1_1347[k];

        t_1348[k] = f_20 * ksk_827[k]
                    + f_3 * pc_y[k] * lsk_1079[k];

        t_1349[k] = pa_x[k] * ksl0_1349[k]
                    - f_14 * pc_x[k] * ksl1_1349[k];
    }

#pragma omp simd aligned(t_1350, t_1351, t_1352, pa_x, pc_x, pc_y, pc_z, ksl0_1350, ksk_792, \
                         ksk_828, ksk_1080, ksl1_1350, lsk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1350[k] = pa_x[k] * ksl0_1350[k]
                    + f_0 * ksk_1080[k]
                    - f_14 * pc_x[k] * ksl1_1350[k];

        t_1351[k] = f_19 * ksk_828[k]
                    + f_3 * pc_y[k] * lsk_1080[k];

        t_1352[k] = f_16 * ksk_792[k]
                    + f_3 * pc_z[k] * lsk_1080[k];
    }

#pragma omp simd aligned(t_1353, t_1354, t_1355, pa_x, pc_x, pc_y, ksl0_1353, ksl0_1355, \
                         ksk_830, ksk_1083, ksk_1085, ksl1_1353, ksl1_1355, \
                         lsk_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1353[k] = pa_x[k] * ksl0_1353[k]
                    + f_20 * ksk_1083[k]
                    - f_14 * pc_x[k] * ksl1_1353[k];

        t_1354[k] = f_19 * ksk_830[k]
                    + f_3 * pc_y[k] * lsk_1082[k];

        t_1355[k] = pa_x[k] * ksl0_1355[k]
                    + f_20 * ksk_1085[k]
                    - f_14 * pc_x[k] * ksl1_1355[k];
    }

#pragma omp simd aligned(t_1356, t_1357, t_1358, pa_x, pc_x, pc_y, pc_z, ksl0_1356, ksk_795, \
                         ksk_833, ksk_1086, ksl1_1356, lsk_1083, \
                         lsk_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1356[k] = pa_x[k] * ksl0_1356[k]
                    + f_19 * ksk_1086[k]
                    - f_14 * pc_x[k] * ksl1_1356[k];

        t_1357[k] = f_16 * ksk_795[k]
                    + f_3 * pc_z[k] * lsk_1083[k];

        t_1358[k] = f_19 * ksk_833[k]
                    + f_3 * pc_y[k] * lsk_1085[k];
    }

#pragma omp simd aligned(t_1359, t_1360, t_1361, pa_x, pc_x, pc_z, ksl0_1359, ksl0_1360, \
                         ksk_798, ksk_1089, ksk_1090, ksl1_1359, ksl1_1360, \
                         lsk_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1359[k] = pa_x[k] * ksl0_1359[k]
                    + f_19 * ksk_1089[k]
                    - f_14 * pc_x[k] * ksl1_1359[k];

        t_1360[k] = pa_x[k] * ksl0_1360[k]
                    + f_18 * ksk_1090[k]
                    - f_14 * pc_x[k] * ksl1_1360[k];

        t_1361[k] = f_16 * ksk_798[k]
                    + f_3 * pc_z[k] * lsk_1086[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, pa_x, pc_x, pc_y, ksl0_1362, ksl0_1364, \
                         ksk_837, ksk_1092, ksk_1094, ksl1_1362, ksl1_1364, \
                         lsk_1089 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = pa_x[k] * ksl0_1362[k]
                    + f_18 * ksk_1092[k]
                    - f_14 * pc_x[k] * ksl1_1362[k];

        t_1363[k] = f_19 * ksk_837[k]
                    + f_3 * pc_y[k] * lsk_1089[k];

        t_1364[k] = pa_x[k] * ksl0_1364[k]
                    + f_18 * ksk_1094[k]
                    - f_14 * pc_x[k] * ksl1_1364[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, pa_x, pc_x, pc_z, ksl0_1365, ksl0_1367, \
                         ksk_802, ksk_1095, ksk_1097, ksl1_1365, ksl1_1367, \
                         lsk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = pa_x[k] * ksl0_1365[k]
                    + f_17 * ksk_1095[k]
                    - f_14 * pc_x[k] * ksl1_1365[k];

        t_1366[k] = f_16 * ksk_802[k]
                    + f_3 * pc_z[k] * lsk_1090[k];

        t_1367[k] = pa_x[k] * ksl0_1367[k]
                    + f_17 * ksk_1097[k]
                    - f_14 * pc_x[k] * ksl1_1367[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksl0,
                                                           const size_t ksk, const size_t ksl1,
                                                           const size_t lsk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_1368 = buffer.data(ksl0 + 1368);
    const auto *ksl0_1370 = buffer.data(ksl0 + 1370);
    const auto *ksl0_1371 = buffer.data(ksl0 + 1371);
    const auto *ksl0_1373 = buffer.data(ksl0 + 1373);
    const auto *ksl0_1374 = buffer.data(ksl0 + 1374);
    const auto *ksl0_1375 = buffer.data(ksl0 + 1375);
    const auto *ksl0_1377 = buffer.data(ksl0 + 1377);
    const auto *ksl0_1386 = buffer.data(ksl0 + 1386);
    const auto *ksl0_1388 = buffer.data(ksl0 + 1388);
    const auto *ksl0_1389 = buffer.data(ksl0 + 1389);
    const auto *ksl0_1390 = buffer.data(ksl0 + 1390);
    const auto *ksl0_1391 = buffer.data(ksl0 + 1391);
    const auto *ksl0_1392 = buffer.data(ksl0 + 1392);
    const auto *ksl0_1394 = buffer.data(ksl0 + 1394);
    const auto *ksl0_1395 = buffer.data(ksl0 + 1395);
    const auto *ksl0_1398 = buffer.data(ksl0 + 1398);
    const auto *ksl0_1400 = buffer.data(ksl0 + 1400);
    const auto *ksl0_1401 = buffer.data(ksl0 + 1401);
    const auto *ksl0_1404 = buffer.data(ksl0 + 1404);
    const auto *ksl0_1405 = buffer.data(ksl0 + 1405);
    const auto *ksl0_1407 = buffer.data(ksl0 + 1407);
    const auto *ksl0_1409 = buffer.data(ksl0 + 1409);
    const auto *ksl0_1410 = buffer.data(ksl0 + 1410);
    const auto *ksl0_1412 = buffer.data(ksl0 + 1412);
    const auto *ksl0_1413 = buffer.data(ksl0 + 1413);
    const auto *ksl0_1415 = buffer.data(ksl0 + 1415);
    const auto *ksl0_1416 = buffer.data(ksl0 + 1416);
    const auto *ksl0_1418 = buffer.data(ksl0 + 1418);
    const auto *ksl0_1419 = buffer.data(ksl0 + 1419);
    const auto *ksl0_1420 = buffer.data(ksl0 + 1420);
    const auto *ksl0_1422 = buffer.data(ksl0 + 1422);
    const auto *ksl0_1431 = buffer.data(ksl0 + 1431);
    const auto *ksl0_1433 = buffer.data(ksl0 + 1433);
    const auto *ksl0_1434 = buffer.data(ksl0 + 1434);
    const auto *ksl0_1435 = buffer.data(ksl0 + 1435);
    const auto *ksl0_1436 = buffer.data(ksl0 + 1436);
    const auto *ksl0_1437 = buffer.data(ksl0 + 1437);
    const auto *ksl0_1439 = buffer.data(ksl0 + 1439);
    const auto *ksl0_1440 = buffer.data(ksl0 + 1440);
    const auto *ksl0_1443 = buffer.data(ksl0 + 1443);
    const auto *ksl0_1445 = buffer.data(ksl0 + 1445);
    const auto *ksl0_1446 = buffer.data(ksl0 + 1446);
    const auto *ksl0_1449 = buffer.data(ksl0 + 1449);
    const auto *ksl0_1450 = buffer.data(ksl0 + 1450);
    const auto *ksl0_1452 = buffer.data(ksl0 + 1452);
    const auto *ksl0_1454 = buffer.data(ksl0 + 1454);
    const auto *ksl0_1455 = buffer.data(ksl0 + 1455);
    const auto *ksl0_1457 = buffer.data(ksl0 + 1457);
    const auto *ksl0_1458 = buffer.data(ksl0 + 1458);
    const auto *ksl0_1460 = buffer.data(ksl0 + 1460);
    const auto *ksl0_1461 = buffer.data(ksl0 + 1461);
    const auto *ksl0_1463 = buffer.data(ksl0 + 1463);
    const auto *ksl0_1464 = buffer.data(ksl0 + 1464);
    const auto *ksl0_1465 = buffer.data(ksl0 + 1465);
    const auto *ksl0_1467 = buffer.data(ksl0 + 1467);
    const auto *ksl0_1476 = buffer.data(ksl0 + 1476);
    const auto *ksl0_1478 = buffer.data(ksl0 + 1478);
    const auto *ksl0_1479 = buffer.data(ksl0 + 1479);
    const auto *ksl0_1480 = buffer.data(ksl0 + 1480);
    const auto *ksl0_1481 = buffer.data(ksl0 + 1481);
    const auto *ksl0_1482 = buffer.data(ksl0 + 1482);

    const auto *ksk_807 = buffer.data(ksk + 807);
    const auto *ksk_820 = buffer.data(ksk + 820);
    const auto *ksk_828 = buffer.data(ksk + 828);
    const auto *ksk_831 = buffer.data(ksk + 831);
    const auto *ksk_834 = buffer.data(ksk + 834);
    const auto *ksk_838 = buffer.data(ksk + 838);
    const auto *ksk_842 = buffer.data(ksk + 842);
    const auto *ksk_843 = buffer.data(ksk + 843);
    const auto *ksk_848 = buffer.data(ksk + 848);
    const auto *ksk_856 = buffer.data(ksk + 856);
    const auto *ksk_863 = buffer.data(ksk + 863);
    const auto *ksk_864 = buffer.data(ksk + 864);
    const auto *ksk_866 = buffer.data(ksk + 866);
    const auto *ksk_867 = buffer.data(ksk + 867);
    const auto *ksk_869 = buffer.data(ksk + 869);
    const auto *ksk_870 = buffer.data(ksk + 870);
    const auto *ksk_873 = buffer.data(ksk + 873);
    const auto *ksk_874 = buffer.data(ksk + 874);
    const auto *ksk_878 = buffer.data(ksk + 878);
    const auto *ksk_879 = buffer.data(ksk + 879);
    const auto *ksk_884 = buffer.data(ksk + 884);
    const auto *ksk_892 = buffer.data(ksk + 892);
    const auto *ksk_899 = buffer.data(ksk + 899);
    const auto *ksk_900 = buffer.data(ksk + 900);
    const auto *ksk_902 = buffer.data(ksk + 902);
    const auto *ksk_905 = buffer.data(ksk + 905);
    const auto *ksk_909 = buffer.data(ksk + 909);
    const auto *ksk_914 = buffer.data(ksk + 914);
    const auto *ksk_920 = buffer.data(ksk + 920);
    const auto *ksk_935 = buffer.data(ksk + 935);
    const auto *ksk_1098 = buffer.data(ksk + 1098);
    const auto *ksk_1100 = buffer.data(ksk + 1100);
    const auto *ksk_1101 = buffer.data(ksk + 1101);
    const auto *ksk_1103 = buffer.data(ksk + 1103);
    const auto *ksk_1104 = buffer.data(ksk + 1104);
    const auto *ksk_1105 = buffer.data(ksk + 1105);
    const auto *ksk_1107 = buffer.data(ksk + 1107);
    const auto *ksk_1108 = buffer.data(ksk + 1108);
    const auto *ksk_1109 = buffer.data(ksk + 1109);
    const auto *ksk_1110 = buffer.data(ksk + 1110);
    const auto *ksk_1111 = buffer.data(ksk + 1111);
    const auto *ksk_1112 = buffer.data(ksk + 1112);
    const auto *ksk_1113 = buffer.data(ksk + 1113);
    const auto *ksk_1114 = buffer.data(ksk + 1114);
    const auto *ksk_1115 = buffer.data(ksk + 1115);
    const auto *ksk_1116 = buffer.data(ksk + 1116);
    const auto *ksk_1119 = buffer.data(ksk + 1119);
    const auto *ksk_1121 = buffer.data(ksk + 1121);
    const auto *ksk_1122 = buffer.data(ksk + 1122);
    const auto *ksk_1125 = buffer.data(ksk + 1125);
    const auto *ksk_1126 = buffer.data(ksk + 1126);
    const auto *ksk_1128 = buffer.data(ksk + 1128);
    const auto *ksk_1130 = buffer.data(ksk + 1130);
    const auto *ksk_1131 = buffer.data(ksk + 1131);
    const auto *ksk_1133 = buffer.data(ksk + 1133);
    const auto *ksk_1134 = buffer.data(ksk + 1134);
    const auto *ksk_1136 = buffer.data(ksk + 1136);
    const auto *ksk_1137 = buffer.data(ksk + 1137);
    const auto *ksk_1139 = buffer.data(ksk + 1139);
    const auto *ksk_1140 = buffer.data(ksk + 1140);
    const auto *ksk_1141 = buffer.data(ksk + 1141);
    const auto *ksk_1143 = buffer.data(ksk + 1143);
    const auto *ksk_1144 = buffer.data(ksk + 1144);
    const auto *ksk_1145 = buffer.data(ksk + 1145);
    const auto *ksk_1146 = buffer.data(ksk + 1146);
    const auto *ksk_1147 = buffer.data(ksk + 1147);
    const auto *ksk_1148 = buffer.data(ksk + 1148);
    const auto *ksk_1149 = buffer.data(ksk + 1149);
    const auto *ksk_1150 = buffer.data(ksk + 1150);
    const auto *ksk_1151 = buffer.data(ksk + 1151);
    const auto *ksk_1152 = buffer.data(ksk + 1152);
    const auto *ksk_1155 = buffer.data(ksk + 1155);
    const auto *ksk_1157 = buffer.data(ksk + 1157);
    const auto *ksk_1158 = buffer.data(ksk + 1158);
    const auto *ksk_1161 = buffer.data(ksk + 1161);
    const auto *ksk_1162 = buffer.data(ksk + 1162);
    const auto *ksk_1164 = buffer.data(ksk + 1164);
    const auto *ksk_1166 = buffer.data(ksk + 1166);
    const auto *ksk_1167 = buffer.data(ksk + 1167);
    const auto *ksk_1169 = buffer.data(ksk + 1169);
    const auto *ksk_1170 = buffer.data(ksk + 1170);
    const auto *ksk_1172 = buffer.data(ksk + 1172);
    const auto *ksk_1173 = buffer.data(ksk + 1173);
    const auto *ksk_1175 = buffer.data(ksk + 1175);
    const auto *ksk_1176 = buffer.data(ksk + 1176);
    const auto *ksk_1177 = buffer.data(ksk + 1177);
    const auto *ksk_1179 = buffer.data(ksk + 1179);
    const auto *ksk_1180 = buffer.data(ksk + 1180);
    const auto *ksk_1181 = buffer.data(ksk + 1181);
    const auto *ksk_1182 = buffer.data(ksk + 1182);
    const auto *ksk_1183 = buffer.data(ksk + 1183);
    const auto *ksk_1184 = buffer.data(ksk + 1184);
    const auto *ksk_1185 = buffer.data(ksk + 1185);
    const auto *ksk_1186 = buffer.data(ksk + 1186);
    const auto *ksk_1187 = buffer.data(ksk + 1187);

    const auto *ksl1_1368 = buffer.data(ksl1 + 1368);
    const auto *ksl1_1370 = buffer.data(ksl1 + 1370);
    const auto *ksl1_1371 = buffer.data(ksl1 + 1371);
    const auto *ksl1_1373 = buffer.data(ksl1 + 1373);
    const auto *ksl1_1374 = buffer.data(ksl1 + 1374);
    const auto *ksl1_1375 = buffer.data(ksl1 + 1375);
    const auto *ksl1_1377 = buffer.data(ksl1 + 1377);
    const auto *ksl1_1386 = buffer.data(ksl1 + 1386);
    const auto *ksl1_1388 = buffer.data(ksl1 + 1388);
    const auto *ksl1_1389 = buffer.data(ksl1 + 1389);
    const auto *ksl1_1390 = buffer.data(ksl1 + 1390);
    const auto *ksl1_1391 = buffer.data(ksl1 + 1391);
    const auto *ksl1_1392 = buffer.data(ksl1 + 1392);
    const auto *ksl1_1394 = buffer.data(ksl1 + 1394);
    const auto *ksl1_1395 = buffer.data(ksl1 + 1395);
    const auto *ksl1_1398 = buffer.data(ksl1 + 1398);
    const auto *ksl1_1400 = buffer.data(ksl1 + 1400);
    const auto *ksl1_1401 = buffer.data(ksl1 + 1401);
    const auto *ksl1_1404 = buffer.data(ksl1 + 1404);
    const auto *ksl1_1405 = buffer.data(ksl1 + 1405);
    const auto *ksl1_1407 = buffer.data(ksl1 + 1407);
    const auto *ksl1_1409 = buffer.data(ksl1 + 1409);
    const auto *ksl1_1410 = buffer.data(ksl1 + 1410);
    const auto *ksl1_1412 = buffer.data(ksl1 + 1412);
    const auto *ksl1_1413 = buffer.data(ksl1 + 1413);
    const auto *ksl1_1415 = buffer.data(ksl1 + 1415);
    const auto *ksl1_1416 = buffer.data(ksl1 + 1416);
    const auto *ksl1_1418 = buffer.data(ksl1 + 1418);
    const auto *ksl1_1419 = buffer.data(ksl1 + 1419);
    const auto *ksl1_1420 = buffer.data(ksl1 + 1420);
    const auto *ksl1_1422 = buffer.data(ksl1 + 1422);
    const auto *ksl1_1431 = buffer.data(ksl1 + 1431);
    const auto *ksl1_1433 = buffer.data(ksl1 + 1433);
    const auto *ksl1_1434 = buffer.data(ksl1 + 1434);
    const auto *ksl1_1435 = buffer.data(ksl1 + 1435);
    const auto *ksl1_1436 = buffer.data(ksl1 + 1436);
    const auto *ksl1_1437 = buffer.data(ksl1 + 1437);
    const auto *ksl1_1439 = buffer.data(ksl1 + 1439);
    const auto *ksl1_1440 = buffer.data(ksl1 + 1440);
    const auto *ksl1_1443 = buffer.data(ksl1 + 1443);
    const auto *ksl1_1445 = buffer.data(ksl1 + 1445);
    const auto *ksl1_1446 = buffer.data(ksl1 + 1446);
    const auto *ksl1_1449 = buffer.data(ksl1 + 1449);
    const auto *ksl1_1450 = buffer.data(ksl1 + 1450);
    const auto *ksl1_1452 = buffer.data(ksl1 + 1452);
    const auto *ksl1_1454 = buffer.data(ksl1 + 1454);
    const auto *ksl1_1455 = buffer.data(ksl1 + 1455);
    const auto *ksl1_1457 = buffer.data(ksl1 + 1457);
    const auto *ksl1_1458 = buffer.data(ksl1 + 1458);
    const auto *ksl1_1460 = buffer.data(ksl1 + 1460);
    const auto *ksl1_1461 = buffer.data(ksl1 + 1461);
    const auto *ksl1_1463 = buffer.data(ksl1 + 1463);
    const auto *ksl1_1464 = buffer.data(ksl1 + 1464);
    const auto *ksl1_1465 = buffer.data(ksl1 + 1465);
    const auto *ksl1_1467 = buffer.data(ksl1 + 1467);
    const auto *ksl1_1476 = buffer.data(ksl1 + 1476);
    const auto *ksl1_1478 = buffer.data(ksl1 + 1478);
    const auto *ksl1_1479 = buffer.data(ksl1 + 1479);
    const auto *ksl1_1480 = buffer.data(ksl1 + 1480);
    const auto *ksl1_1481 = buffer.data(ksl1 + 1481);
    const auto *ksl1_1482 = buffer.data(ksl1 + 1482);

    const auto *lsk_1094 = buffer.data(lsk + 1094);
    const auto *lsk_1095 = buffer.data(lsk + 1095);
    const auto *lsk_1100 = buffer.data(lsk + 1100);
    const auto *lsk_1108 = buffer.data(lsk + 1108);
    const auto *lsk_1109 = buffer.data(lsk + 1109);
    const auto *lsk_1110 = buffer.data(lsk + 1110);
    const auto *lsk_1111 = buffer.data(lsk + 1111);
    const auto *lsk_1112 = buffer.data(lsk + 1112);
    const auto *lsk_1113 = buffer.data(lsk + 1113);
    const auto *lsk_1114 = buffer.data(lsk + 1114);
    const auto *lsk_1115 = buffer.data(lsk + 1115);
    const auto *lsk_1116 = buffer.data(lsk + 1116);
    const auto *lsk_1118 = buffer.data(lsk + 1118);
    const auto *lsk_1119 = buffer.data(lsk + 1119);
    const auto *lsk_1121 = buffer.data(lsk + 1121);
    const auto *lsk_1122 = buffer.data(lsk + 1122);
    const auto *lsk_1125 = buffer.data(lsk + 1125);
    const auto *lsk_1126 = buffer.data(lsk + 1126);
    const auto *lsk_1130 = buffer.data(lsk + 1130);
    const auto *lsk_1131 = buffer.data(lsk + 1131);
    const auto *lsk_1136 = buffer.data(lsk + 1136);
    const auto *lsk_1144 = buffer.data(lsk + 1144);
    const auto *lsk_1145 = buffer.data(lsk + 1145);
    const auto *lsk_1146 = buffer.data(lsk + 1146);
    const auto *lsk_1147 = buffer.data(lsk + 1147);
    const auto *lsk_1148 = buffer.data(lsk + 1148);
    const auto *lsk_1149 = buffer.data(lsk + 1149);
    const auto *lsk_1150 = buffer.data(lsk + 1150);
    const auto *lsk_1151 = buffer.data(lsk + 1151);
    const auto *lsk_1152 = buffer.data(lsk + 1152);
    const auto *lsk_1154 = buffer.data(lsk + 1154);
    const auto *lsk_1155 = buffer.data(lsk + 1155);
    const auto *lsk_1157 = buffer.data(lsk + 1157);
    const auto *lsk_1158 = buffer.data(lsk + 1158);
    const auto *lsk_1161 = buffer.data(lsk + 1161);
    const auto *lsk_1162 = buffer.data(lsk + 1162);
    const auto *lsk_1166 = buffer.data(lsk + 1166);
    const auto *lsk_1167 = buffer.data(lsk + 1167);
    const auto *lsk_1172 = buffer.data(lsk + 1172);
    const auto *lsk_1180 = buffer.data(lsk + 1180);
    const auto *lsk_1181 = buffer.data(lsk + 1181);
    const auto *lsk_1182 = buffer.data(lsk + 1182);
    const auto *lsk_1183 = buffer.data(lsk + 1183);
    const auto *lsk_1184 = buffer.data(lsk + 1184);
    const auto *lsk_1185 = buffer.data(lsk + 1185);
    const auto *lsk_1186 = buffer.data(lsk + 1186);
    const auto *lsk_1187 = buffer.data(lsk + 1187);

#pragma omp simd aligned(t_1368, t_1369, t_1370, pa_x, pc_x, pc_y, ksl0_1368, ksl0_1370, \
                         ksk_842, ksk_1098, ksk_1100, ksl1_1368, ksl1_1370, \
                         lsk_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1368[k] = pa_x[k] * ksl0_1368[k]
                    + f_17 * ksk_1098[k]
                    - f_14 * pc_x[k] * ksl1_1368[k];

        t_1369[k] = f_19 * ksk_842[k]
                    + f_3 * pc_y[k] * lsk_1094[k];

        t_1370[k] = pa_x[k] * ksl0_1370[k]
                    + f_17 * ksk_1100[k]
                    - f_14 * pc_x[k] * ksl1_1370[k];
    }

#pragma omp simd aligned(t_1371, t_1372, t_1373, pa_x, pc_x, pc_z, ksl0_1371, ksl0_1373, \
                         ksk_807, ksk_1101, ksk_1103, ksl1_1371, ksl1_1373, \
                         lsk_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1371[k] = pa_x[k] * ksl0_1371[k]
                    + f_16 * ksk_1101[k]
                    - f_14 * pc_x[k] * ksl1_1371[k];

        t_1372[k] = f_16 * ksk_807[k]
                    + f_3 * pc_z[k] * lsk_1095[k];

        t_1373[k] = pa_x[k] * ksl0_1373[k]
                    + f_16 * ksk_1103[k]
                    - f_14 * pc_x[k] * ksl1_1373[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, pa_x, pc_x, pc_y, ksl0_1374, ksl0_1375, \
                         ksk_848, ksk_1104, ksk_1105, ksl1_1374, ksl1_1375, \
                         lsk_1100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = pa_x[k] * ksl0_1374[k]
                    + f_16 * ksk_1104[k]
                    - f_14 * pc_x[k] * ksl1_1374[k];

        t_1375[k] = pa_x[k] * ksl0_1375[k]
                    + f_16 * ksk_1105[k]
                    - f_14 * pc_x[k] * ksl1_1375[k];

        t_1376[k] = f_19 * ksk_848[k]
                    + f_3 * pc_y[k] * lsk_1100[k];
    }

#pragma omp simd aligned(t_1377, t_1378, t_1379, t_1380, pa_x, pc_x, ksl0_1377, ksk_1107, \
                         ksk_1108, ksk_1109, ksk_1110, ksl1_1377, lsk_1108, lsk_1109, \
                         lsk_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1377[k] = pa_x[k] * ksl0_1377[k]
                    + f_16 * ksk_1107[k]
                    - f_14 * pc_x[k] * ksl1_1377[k];

        t_1378[k] = f_15 * ksk_1108[k]
                    + f_3 * pc_x[k] * lsk_1108[k];

        t_1379[k] = f_15 * ksk_1109[k]
                    + f_3 * pc_x[k] * lsk_1109[k];

        t_1380[k] = f_15 * ksk_1110[k]
                    + f_3 * pc_x[k] * lsk_1110[k];
    }

#pragma omp simd aligned(t_1381, t_1382, t_1383, t_1384, t_1385, pc_x, ksk_1111, ksk_1112, \
                         ksk_1113, ksk_1114, ksk_1115, lsk_1111, lsk_1112, lsk_1113, lsk_1114, \
                         lsk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1381[k] = f_15 * ksk_1111[k]
                    + f_3 * pc_x[k] * lsk_1111[k];

        t_1382[k] = f_15 * ksk_1112[k]
                    + f_3 * pc_x[k] * lsk_1112[k];

        t_1383[k] = f_15 * ksk_1113[k]
                    + f_3 * pc_x[k] * lsk_1113[k];

        t_1384[k] = f_15 * ksk_1114[k]
                    + f_3 * pc_x[k] * lsk_1114[k];

        t_1385[k] = f_15 * ksk_1115[k]
                    + f_3 * pc_x[k] * lsk_1115[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, t_1389, pa_x, pc_x, pc_z, ksl0_1386, \
                         ksl0_1388, ksl0_1389, ksk_820, ksl1_1386, ksl1_1388, ksl1_1389, \
                         lsk_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = pa_x[k] * ksl0_1386[k]
                    - f_14 * pc_x[k] * ksl1_1386[k];

        t_1387[k] = f_16 * ksk_820[k]
                    + f_3 * pc_z[k] * lsk_1108[k];

        t_1388[k] = pa_x[k] * ksl0_1388[k]
                    - f_14 * pc_x[k] * ksl1_1388[k];

        t_1389[k] = pa_x[k] * ksl0_1389[k]
                    - f_14 * pc_x[k] * ksl1_1389[k];
    }

#pragma omp simd aligned(t_1390, t_1391, t_1392, t_1393, pa_x, pc_x, pc_y, ksl0_1390, \
                         ksl0_1391, ksl0_1392, ksk_863, ksl1_1390, ksl1_1391, ksl1_1392, \
                         lsk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1390[k] = pa_x[k] * ksl0_1390[k]
                    - f_14 * pc_x[k] * ksl1_1390[k];

        t_1391[k] = pa_x[k] * ksl0_1391[k]
                    - f_14 * pc_x[k] * ksl1_1391[k];

        t_1392[k] = pa_x[k] * ksl0_1392[k]
                    - f_14 * pc_x[k] * ksl1_1392[k];

        t_1393[k] = f_19 * ksk_863[k]
                    + f_3 * pc_y[k] * lsk_1115[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, t_1397, pa_x, pc_x, pc_y, pc_z, ksl0_1394, \
                         ksl0_1395, ksk_828, ksk_864, ksk_1116, ksl1_1394, ksl1_1395, \
                         lsk_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = pa_x[k] * ksl0_1394[k]
                    - f_14 * pc_x[k] * ksl1_1394[k];

        t_1395[k] = pa_x[k] * ksl0_1395[k]
                    + f_0 * ksk_1116[k]
                    - f_14 * pc_x[k] * ksl1_1395[k];

        t_1396[k] = f_18 * ksk_864[k]
                    + f_3 * pc_y[k] * lsk_1116[k];

        t_1397[k] = f_17 * ksk_828[k]
                    + f_3 * pc_z[k] * lsk_1116[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, pa_x, pc_x, pc_y, ksl0_1398, ksl0_1400, \
                         ksk_866, ksk_1119, ksk_1121, ksl1_1398, ksl1_1400, \
                         lsk_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = pa_x[k] * ksl0_1398[k]
                    + f_20 * ksk_1119[k]
                    - f_14 * pc_x[k] * ksl1_1398[k];

        t_1399[k] = f_18 * ksk_866[k]
                    + f_3 * pc_y[k] * lsk_1118[k];

        t_1400[k] = pa_x[k] * ksl0_1400[k]
                    + f_20 * ksk_1121[k]
                    - f_14 * pc_x[k] * ksl1_1400[k];
    }

#pragma omp simd aligned(t_1401, t_1402, t_1403, pa_x, pc_x, pc_y, pc_z, ksl0_1401, ksk_831, \
                         ksk_869, ksk_1122, ksl1_1401, lsk_1119, \
                         lsk_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1401[k] = pa_x[k] * ksl0_1401[k]
                    + f_19 * ksk_1122[k]
                    - f_14 * pc_x[k] * ksl1_1401[k];

        t_1402[k] = f_17 * ksk_831[k]
                    + f_3 * pc_z[k] * lsk_1119[k];

        t_1403[k] = f_18 * ksk_869[k]
                    + f_3 * pc_y[k] * lsk_1121[k];
    }

#pragma omp simd aligned(t_1404, t_1405, t_1406, pa_x, pc_x, pc_z, ksl0_1404, ksl0_1405, \
                         ksk_834, ksk_1125, ksk_1126, ksl1_1404, ksl1_1405, \
                         lsk_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1404[k] = pa_x[k] * ksl0_1404[k]
                    + f_19 * ksk_1125[k]
                    - f_14 * pc_x[k] * ksl1_1404[k];

        t_1405[k] = pa_x[k] * ksl0_1405[k]
                    + f_18 * ksk_1126[k]
                    - f_14 * pc_x[k] * ksl1_1405[k];

        t_1406[k] = f_17 * ksk_834[k]
                    + f_3 * pc_z[k] * lsk_1122[k];
    }

#pragma omp simd aligned(t_1407, t_1408, t_1409, pa_x, pc_x, pc_y, ksl0_1407, ksl0_1409, \
                         ksk_873, ksk_1128, ksk_1130, ksl1_1407, ksl1_1409, \
                         lsk_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1407[k] = pa_x[k] * ksl0_1407[k]
                    + f_18 * ksk_1128[k]
                    - f_14 * pc_x[k] * ksl1_1407[k];

        t_1408[k] = f_18 * ksk_873[k]
                    + f_3 * pc_y[k] * lsk_1125[k];

        t_1409[k] = pa_x[k] * ksl0_1409[k]
                    + f_18 * ksk_1130[k]
                    - f_14 * pc_x[k] * ksl1_1409[k];
    }

#pragma omp simd aligned(t_1410, t_1411, t_1412, pa_x, pc_x, pc_z, ksl0_1410, ksl0_1412, \
                         ksk_838, ksk_1131, ksk_1133, ksl1_1410, ksl1_1412, \
                         lsk_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1410[k] = pa_x[k] * ksl0_1410[k]
                    + f_17 * ksk_1131[k]
                    - f_14 * pc_x[k] * ksl1_1410[k];

        t_1411[k] = f_17 * ksk_838[k]
                    + f_3 * pc_z[k] * lsk_1126[k];

        t_1412[k] = pa_x[k] * ksl0_1412[k]
                    + f_17 * ksk_1133[k]
                    - f_14 * pc_x[k] * ksl1_1412[k];
    }

#pragma omp simd aligned(t_1413, t_1414, t_1415, pa_x, pc_x, pc_y, ksl0_1413, ksl0_1415, \
                         ksk_878, ksk_1134, ksk_1136, ksl1_1413, ksl1_1415, \
                         lsk_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1413[k] = pa_x[k] * ksl0_1413[k]
                    + f_17 * ksk_1134[k]
                    - f_14 * pc_x[k] * ksl1_1413[k];

        t_1414[k] = f_18 * ksk_878[k]
                    + f_3 * pc_y[k] * lsk_1130[k];

        t_1415[k] = pa_x[k] * ksl0_1415[k]
                    + f_17 * ksk_1136[k]
                    - f_14 * pc_x[k] * ksl1_1415[k];
    }

#pragma omp simd aligned(t_1416, t_1417, t_1418, pa_x, pc_x, pc_z, ksl0_1416, ksl0_1418, \
                         ksk_843, ksk_1137, ksk_1139, ksl1_1416, ksl1_1418, \
                         lsk_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1416[k] = pa_x[k] * ksl0_1416[k]
                    + f_16 * ksk_1137[k]
                    - f_14 * pc_x[k] * ksl1_1416[k];

        t_1417[k] = f_17 * ksk_843[k]
                    + f_3 * pc_z[k] * lsk_1131[k];

        t_1418[k] = pa_x[k] * ksl0_1418[k]
                    + f_16 * ksk_1139[k]
                    - f_14 * pc_x[k] * ksl1_1418[k];
    }

#pragma omp simd aligned(t_1419, t_1420, t_1421, pa_x, pc_x, pc_y, ksl0_1419, ksl0_1420, \
                         ksk_884, ksk_1140, ksk_1141, ksl1_1419, ksl1_1420, \
                         lsk_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1419[k] = pa_x[k] * ksl0_1419[k]
                    + f_16 * ksk_1140[k]
                    - f_14 * pc_x[k] * ksl1_1419[k];

        t_1420[k] = pa_x[k] * ksl0_1420[k]
                    + f_16 * ksk_1141[k]
                    - f_14 * pc_x[k] * ksl1_1420[k];

        t_1421[k] = f_18 * ksk_884[k]
                    + f_3 * pc_y[k] * lsk_1136[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, t_1425, pa_x, pc_x, ksl0_1422, ksk_1143, \
                         ksk_1144, ksk_1145, ksk_1146, ksl1_1422, lsk_1144, lsk_1145, \
                         lsk_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = pa_x[k] * ksl0_1422[k]
                    + f_16 * ksk_1143[k]
                    - f_14 * pc_x[k] * ksl1_1422[k];

        t_1423[k] = f_15 * ksk_1144[k]
                    + f_3 * pc_x[k] * lsk_1144[k];

        t_1424[k] = f_15 * ksk_1145[k]
                    + f_3 * pc_x[k] * lsk_1145[k];

        t_1425[k] = f_15 * ksk_1146[k]
                    + f_3 * pc_x[k] * lsk_1146[k];
    }

#pragma omp simd aligned(t_1426, t_1427, t_1428, t_1429, t_1430, pc_x, ksk_1147, ksk_1148, \
                         ksk_1149, ksk_1150, ksk_1151, lsk_1147, lsk_1148, lsk_1149, lsk_1150, \
                         lsk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1426[k] = f_15 * ksk_1147[k]
                    + f_3 * pc_x[k] * lsk_1147[k];

        t_1427[k] = f_15 * ksk_1148[k]
                    + f_3 * pc_x[k] * lsk_1148[k];

        t_1428[k] = f_15 * ksk_1149[k]
                    + f_3 * pc_x[k] * lsk_1149[k];

        t_1429[k] = f_15 * ksk_1150[k]
                    + f_3 * pc_x[k] * lsk_1150[k];

        t_1430[k] = f_15 * ksk_1151[k]
                    + f_3 * pc_x[k] * lsk_1151[k];
    }

#pragma omp simd aligned(t_1431, t_1432, t_1433, t_1434, pa_x, pc_x, pc_z, ksl0_1431, \
                         ksl0_1433, ksl0_1434, ksk_856, ksl1_1431, ksl1_1433, ksl1_1434, \
                         lsk_1144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1431[k] = pa_x[k] * ksl0_1431[k]
                    - f_14 * pc_x[k] * ksl1_1431[k];

        t_1432[k] = f_17 * ksk_856[k]
                    + f_3 * pc_z[k] * lsk_1144[k];

        t_1433[k] = pa_x[k] * ksl0_1433[k]
                    - f_14 * pc_x[k] * ksl1_1433[k];

        t_1434[k] = pa_x[k] * ksl0_1434[k]
                    - f_14 * pc_x[k] * ksl1_1434[k];
    }

#pragma omp simd aligned(t_1435, t_1436, t_1437, t_1438, pa_x, pc_x, pc_y, ksl0_1435, \
                         ksl0_1436, ksl0_1437, ksk_899, ksl1_1435, ksl1_1436, ksl1_1437, \
                         lsk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1435[k] = pa_x[k] * ksl0_1435[k]
                    - f_14 * pc_x[k] * ksl1_1435[k];

        t_1436[k] = pa_x[k] * ksl0_1436[k]
                    - f_14 * pc_x[k] * ksl1_1436[k];

        t_1437[k] = pa_x[k] * ksl0_1437[k]
                    - f_14 * pc_x[k] * ksl1_1437[k];

        t_1438[k] = f_18 * ksk_899[k]
                    + f_3 * pc_y[k] * lsk_1151[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, t_1442, pa_x, pc_x, pc_y, pc_z, ksl0_1439, \
                         ksl0_1440, ksk_864, ksk_900, ksk_1152, ksl1_1439, ksl1_1440, \
                         lsk_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = pa_x[k] * ksl0_1439[k]
                    - f_14 * pc_x[k] * ksl1_1439[k];

        t_1440[k] = pa_x[k] * ksl0_1440[k]
                    + f_0 * ksk_1152[k]
                    - f_14 * pc_x[k] * ksl1_1440[k];

        t_1441[k] = f_17 * ksk_900[k]
                    + f_3 * pc_y[k] * lsk_1152[k];

        t_1442[k] = f_18 * ksk_864[k]
                    + f_3 * pc_z[k] * lsk_1152[k];
    }

#pragma omp simd aligned(t_1443, t_1444, t_1445, pa_x, pc_x, pc_y, ksl0_1443, ksl0_1445, \
                         ksk_902, ksk_1155, ksk_1157, ksl1_1443, ksl1_1445, \
                         lsk_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1443[k] = pa_x[k] * ksl0_1443[k]
                    + f_20 * ksk_1155[k]
                    - f_14 * pc_x[k] * ksl1_1443[k];

        t_1444[k] = f_17 * ksk_902[k]
                    + f_3 * pc_y[k] * lsk_1154[k];

        t_1445[k] = pa_x[k] * ksl0_1445[k]
                    + f_20 * ksk_1157[k]
                    - f_14 * pc_x[k] * ksl1_1445[k];
    }

#pragma omp simd aligned(t_1446, t_1447, t_1448, pa_x, pc_x, pc_y, pc_z, ksl0_1446, ksk_867, \
                         ksk_905, ksk_1158, ksl1_1446, lsk_1155, \
                         lsk_1157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1446[k] = pa_x[k] * ksl0_1446[k]
                    + f_19 * ksk_1158[k]
                    - f_14 * pc_x[k] * ksl1_1446[k];

        t_1447[k] = f_18 * ksk_867[k]
                    + f_3 * pc_z[k] * lsk_1155[k];

        t_1448[k] = f_17 * ksk_905[k]
                    + f_3 * pc_y[k] * lsk_1157[k];
    }

#pragma omp simd aligned(t_1449, t_1450, t_1451, pa_x, pc_x, pc_z, ksl0_1449, ksl0_1450, \
                         ksk_870, ksk_1161, ksk_1162, ksl1_1449, ksl1_1450, \
                         lsk_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1449[k] = pa_x[k] * ksl0_1449[k]
                    + f_19 * ksk_1161[k]
                    - f_14 * pc_x[k] * ksl1_1449[k];

        t_1450[k] = pa_x[k] * ksl0_1450[k]
                    + f_18 * ksk_1162[k]
                    - f_14 * pc_x[k] * ksl1_1450[k];

        t_1451[k] = f_18 * ksk_870[k]
                    + f_3 * pc_z[k] * lsk_1158[k];
    }

#pragma omp simd aligned(t_1452, t_1453, t_1454, pa_x, pc_x, pc_y, ksl0_1452, ksl0_1454, \
                         ksk_909, ksk_1164, ksk_1166, ksl1_1452, ksl1_1454, \
                         lsk_1161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1452[k] = pa_x[k] * ksl0_1452[k]
                    + f_18 * ksk_1164[k]
                    - f_14 * pc_x[k] * ksl1_1452[k];

        t_1453[k] = f_17 * ksk_909[k]
                    + f_3 * pc_y[k] * lsk_1161[k];

        t_1454[k] = pa_x[k] * ksl0_1454[k]
                    + f_18 * ksk_1166[k]
                    - f_14 * pc_x[k] * ksl1_1454[k];
    }

#pragma omp simd aligned(t_1455, t_1456, t_1457, pa_x, pc_x, pc_z, ksl0_1455, ksl0_1457, \
                         ksk_874, ksk_1167, ksk_1169, ksl1_1455, ksl1_1457, \
                         lsk_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1455[k] = pa_x[k] * ksl0_1455[k]
                    + f_17 * ksk_1167[k]
                    - f_14 * pc_x[k] * ksl1_1455[k];

        t_1456[k] = f_18 * ksk_874[k]
                    + f_3 * pc_z[k] * lsk_1162[k];

        t_1457[k] = pa_x[k] * ksl0_1457[k]
                    + f_17 * ksk_1169[k]
                    - f_14 * pc_x[k] * ksl1_1457[k];
    }

#pragma omp simd aligned(t_1458, t_1459, t_1460, pa_x, pc_x, pc_y, ksl0_1458, ksl0_1460, \
                         ksk_914, ksk_1170, ksk_1172, ksl1_1458, ksl1_1460, \
                         lsk_1166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1458[k] = pa_x[k] * ksl0_1458[k]
                    + f_17 * ksk_1170[k]
                    - f_14 * pc_x[k] * ksl1_1458[k];

        t_1459[k] = f_17 * ksk_914[k]
                    + f_3 * pc_y[k] * lsk_1166[k];

        t_1460[k] = pa_x[k] * ksl0_1460[k]
                    + f_17 * ksk_1172[k]
                    - f_14 * pc_x[k] * ksl1_1460[k];
    }

#pragma omp simd aligned(t_1461, t_1462, t_1463, pa_x, pc_x, pc_z, ksl0_1461, ksl0_1463, \
                         ksk_879, ksk_1173, ksk_1175, ksl1_1461, ksl1_1463, \
                         lsk_1167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1461[k] = pa_x[k] * ksl0_1461[k]
                    + f_16 * ksk_1173[k]
                    - f_14 * pc_x[k] * ksl1_1461[k];

        t_1462[k] = f_18 * ksk_879[k]
                    + f_3 * pc_z[k] * lsk_1167[k];

        t_1463[k] = pa_x[k] * ksl0_1463[k]
                    + f_16 * ksk_1175[k]
                    - f_14 * pc_x[k] * ksl1_1463[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, pa_x, pc_x, pc_y, ksl0_1464, ksl0_1465, \
                         ksk_920, ksk_1176, ksk_1177, ksl1_1464, ksl1_1465, \
                         lsk_1172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = pa_x[k] * ksl0_1464[k]
                    + f_16 * ksk_1176[k]
                    - f_14 * pc_x[k] * ksl1_1464[k];

        t_1465[k] = pa_x[k] * ksl0_1465[k]
                    + f_16 * ksk_1177[k]
                    - f_14 * pc_x[k] * ksl1_1465[k];

        t_1466[k] = f_17 * ksk_920[k]
                    + f_3 * pc_y[k] * lsk_1172[k];
    }

#pragma omp simd aligned(t_1467, t_1468, t_1469, t_1470, pa_x, pc_x, ksl0_1467, ksk_1179, \
                         ksk_1180, ksk_1181, ksk_1182, ksl1_1467, lsk_1180, lsk_1181, \
                         lsk_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1467[k] = pa_x[k] * ksl0_1467[k]
                    + f_16 * ksk_1179[k]
                    - f_14 * pc_x[k] * ksl1_1467[k];

        t_1468[k] = f_15 * ksk_1180[k]
                    + f_3 * pc_x[k] * lsk_1180[k];

        t_1469[k] = f_15 * ksk_1181[k]
                    + f_3 * pc_x[k] * lsk_1181[k];

        t_1470[k] = f_15 * ksk_1182[k]
                    + f_3 * pc_x[k] * lsk_1182[k];
    }

#pragma omp simd aligned(t_1471, t_1472, t_1473, t_1474, t_1475, pc_x, ksk_1183, ksk_1184, \
                         ksk_1185, ksk_1186, ksk_1187, lsk_1183, lsk_1184, lsk_1185, lsk_1186, \
                         lsk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1471[k] = f_15 * ksk_1183[k]
                    + f_3 * pc_x[k] * lsk_1183[k];

        t_1472[k] = f_15 * ksk_1184[k]
                    + f_3 * pc_x[k] * lsk_1184[k];

        t_1473[k] = f_15 * ksk_1185[k]
                    + f_3 * pc_x[k] * lsk_1185[k];

        t_1474[k] = f_15 * ksk_1186[k]
                    + f_3 * pc_x[k] * lsk_1186[k];

        t_1475[k] = f_15 * ksk_1187[k]
                    + f_3 * pc_x[k] * lsk_1187[k];
    }

#pragma omp simd aligned(t_1476, t_1477, t_1478, t_1479, pa_x, pc_x, pc_z, ksl0_1476, \
                         ksl0_1478, ksl0_1479, ksk_892, ksl1_1476, ksl1_1478, ksl1_1479, \
                         lsk_1180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1476[k] = pa_x[k] * ksl0_1476[k]
                    - f_14 * pc_x[k] * ksl1_1476[k];

        t_1477[k] = f_18 * ksk_892[k]
                    + f_3 * pc_z[k] * lsk_1180[k];

        t_1478[k] = pa_x[k] * ksl0_1478[k]
                    - f_14 * pc_x[k] * ksl1_1478[k];

        t_1479[k] = pa_x[k] * ksl0_1479[k]
                    - f_14 * pc_x[k] * ksl1_1479[k];
    }

#pragma omp simd aligned(t_1480, t_1481, t_1482, t_1483, pa_x, pc_x, pc_y, ksl0_1480, \
                         ksl0_1481, ksl0_1482, ksk_935, ksl1_1480, ksl1_1481, ksl1_1482, \
                         lsk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1480[k] = pa_x[k] * ksl0_1480[k]
                    - f_14 * pc_x[k] * ksl1_1480[k];

        t_1481[k] = pa_x[k] * ksl0_1481[k]
                    - f_14 * pc_x[k] * ksl1_1481[k];

        t_1482[k] = pa_x[k] * ksl0_1482[k]
                    - f_14 * pc_x[k] * ksl1_1482[k];

        t_1483[k] = f_17 * ksk_935[k]
                    + f_3 * pc_y[k] * lsk_1187[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksl0,
                                                           const size_t ksk, const size_t ksl1,
                                                           const size_t lsi0, const size_t lsi1,
                                                           const size_t lsk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 3.5 / q;

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
    auto *t_1596 = buffer.data(target + 1596);
    auto *t_1597 = buffer.data(target + 1597);
    auto *t_1598 = buffer.data(target + 1598);
    auto *t_1599 = buffer.data(target + 1599);
    auto *t_1600 = buffer.data(target + 1600);
    auto *t_1601 = buffer.data(target + 1601);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_1215 = buffer.data(ksl0 + 1215);
    const auto *ksl0_1220 = buffer.data(ksl0 + 1220);
    const auto *ksl0_1224 = buffer.data(ksl0 + 1224);
    const auto *ksl0_1229 = buffer.data(ksl0 + 1229);
    const auto *ksl0_1235 = buffer.data(ksl0 + 1235);
    const auto *ksl0_1242 = buffer.data(ksl0 + 1242);
    const auto *ksl0_1484 = buffer.data(ksl0 + 1484);
    const auto *ksl0_1485 = buffer.data(ksl0 + 1485);
    const auto *ksl0_1488 = buffer.data(ksl0 + 1488);
    const auto *ksl0_1490 = buffer.data(ksl0 + 1490);
    const auto *ksl0_1491 = buffer.data(ksl0 + 1491);
    const auto *ksl0_1494 = buffer.data(ksl0 + 1494);
    const auto *ksl0_1495 = buffer.data(ksl0 + 1495);
    const auto *ksl0_1497 = buffer.data(ksl0 + 1497);
    const auto *ksl0_1499 = buffer.data(ksl0 + 1499);
    const auto *ksl0_1500 = buffer.data(ksl0 + 1500);
    const auto *ksl0_1502 = buffer.data(ksl0 + 1502);
    const auto *ksl0_1503 = buffer.data(ksl0 + 1503);
    const auto *ksl0_1505 = buffer.data(ksl0 + 1505);
    const auto *ksl0_1506 = buffer.data(ksl0 + 1506);
    const auto *ksl0_1508 = buffer.data(ksl0 + 1508);
    const auto *ksl0_1509 = buffer.data(ksl0 + 1509);
    const auto *ksl0_1510 = buffer.data(ksl0 + 1510);
    const auto *ksl0_1512 = buffer.data(ksl0 + 1512);
    const auto *ksl0_1521 = buffer.data(ksl0 + 1521);
    const auto *ksl0_1523 = buffer.data(ksl0 + 1523);
    const auto *ksl0_1524 = buffer.data(ksl0 + 1524);
    const auto *ksl0_1525 = buffer.data(ksl0 + 1525);
    const auto *ksl0_1526 = buffer.data(ksl0 + 1526);
    const auto *ksl0_1527 = buffer.data(ksl0 + 1527);
    const auto *ksl0_1529 = buffer.data(ksl0 + 1529);
    const auto *ksl0_1533 = buffer.data(ksl0 + 1533);
    const auto *ksl0_1536 = buffer.data(ksl0 + 1536);
    const auto *ksl0_1540 = buffer.data(ksl0 + 1540);
    const auto *ksl0_1542 = buffer.data(ksl0 + 1542);
    const auto *ksl0_1545 = buffer.data(ksl0 + 1545);
    const auto *ksl0_1547 = buffer.data(ksl0 + 1547);
    const auto *ksl0_1548 = buffer.data(ksl0 + 1548);
    const auto *ksl0_1551 = buffer.data(ksl0 + 1551);
    const auto *ksl0_1553 = buffer.data(ksl0 + 1553);
    const auto *ksl0_1554 = buffer.data(ksl0 + 1554);
    const auto *ksl0_1555 = buffer.data(ksl0 + 1555);
    const auto *ksl0_1566 = buffer.data(ksl0 + 1566);
    const auto *ksl0_1568 = buffer.data(ksl0 + 1568);
    const auto *ksl0_1569 = buffer.data(ksl0 + 1569);
    const auto *ksl0_1570 = buffer.data(ksl0 + 1570);
    const auto *ksl0_1571 = buffer.data(ksl0 + 1571);
    const auto *ksl0_1572 = buffer.data(ksl0 + 1572);
    const auto *ksl0_1574 = buffer.data(ksl0 + 1574);
    const auto *ksl0_1575 = buffer.data(ksl0 + 1575);
    const auto *ksl0_1580 = buffer.data(ksl0 + 1580);
    const auto *ksl0_1584 = buffer.data(ksl0 + 1584);
    const auto *ksl0_1589 = buffer.data(ksl0 + 1589);
    const auto *ksl0_1595 = buffer.data(ksl0 + 1595);

    const auto *ksk_900 = buffer.data(ksk + 900);
    const auto *ksk_903 = buffer.data(ksk + 903);
    const auto *ksk_906 = buffer.data(ksk + 906);
    const auto *ksk_910 = buffer.data(ksk + 910);
    const auto *ksk_915 = buffer.data(ksk + 915);
    const auto *ksk_928 = buffer.data(ksk + 928);
    const auto *ksk_936 = buffer.data(ksk + 936);
    const auto *ksk_938 = buffer.data(ksk + 938);
    const auto *ksk_939 = buffer.data(ksk + 939);
    const auto *ksk_941 = buffer.data(ksk + 941);
    const auto *ksk_942 = buffer.data(ksk + 942);
    const auto *ksk_945 = buffer.data(ksk + 945);
    const auto *ksk_946 = buffer.data(ksk + 946);
    const auto *ksk_950 = buffer.data(ksk + 950);
    const auto *ksk_951 = buffer.data(ksk + 951);
    const auto *ksk_956 = buffer.data(ksk + 956);
    const auto *ksk_964 = buffer.data(ksk + 964);
    const auto *ksk_971 = buffer.data(ksk + 971);
    const auto *ksk_972 = buffer.data(ksk + 972);
    const auto *ksk_974 = buffer.data(ksk + 974);
    const auto *ksk_977 = buffer.data(ksk + 977);
    const auto *ksk_981 = buffer.data(ksk + 981);
    const auto *ksk_986 = buffer.data(ksk + 986);
    const auto *ksk_992 = buffer.data(ksk + 992);
    const auto *ksk_1007 = buffer.data(ksk + 1007);
    const auto *ksk_1188 = buffer.data(ksk + 1188);
    const auto *ksk_1191 = buffer.data(ksk + 1191);
    const auto *ksk_1193 = buffer.data(ksk + 1193);
    const auto *ksk_1194 = buffer.data(ksk + 1194);
    const auto *ksk_1197 = buffer.data(ksk + 1197);
    const auto *ksk_1198 = buffer.data(ksk + 1198);
    const auto *ksk_1200 = buffer.data(ksk + 1200);
    const auto *ksk_1202 = buffer.data(ksk + 1202);
    const auto *ksk_1203 = buffer.data(ksk + 1203);
    const auto *ksk_1205 = buffer.data(ksk + 1205);
    const auto *ksk_1206 = buffer.data(ksk + 1206);
    const auto *ksk_1208 = buffer.data(ksk + 1208);
    const auto *ksk_1209 = buffer.data(ksk + 1209);
    const auto *ksk_1211 = buffer.data(ksk + 1211);
    const auto *ksk_1212 = buffer.data(ksk + 1212);
    const auto *ksk_1213 = buffer.data(ksk + 1213);
    const auto *ksk_1215 = buffer.data(ksk + 1215);
    const auto *ksk_1216 = buffer.data(ksk + 1216);
    const auto *ksk_1217 = buffer.data(ksk + 1217);
    const auto *ksk_1218 = buffer.data(ksk + 1218);
    const auto *ksk_1219 = buffer.data(ksk + 1219);
    const auto *ksk_1220 = buffer.data(ksk + 1220);
    const auto *ksk_1221 = buffer.data(ksk + 1221);
    const auto *ksk_1222 = buffer.data(ksk + 1222);
    const auto *ksk_1223 = buffer.data(ksk + 1223);
    const auto *ksk_1227 = buffer.data(ksk + 1227);
    const auto *ksk_1230 = buffer.data(ksk + 1230);
    const auto *ksk_1234 = buffer.data(ksk + 1234);
    const auto *ksk_1236 = buffer.data(ksk + 1236);
    const auto *ksk_1239 = buffer.data(ksk + 1239);
    const auto *ksk_1241 = buffer.data(ksk + 1241);
    const auto *ksk_1242 = buffer.data(ksk + 1242);
    const auto *ksk_1245 = buffer.data(ksk + 1245);
    const auto *ksk_1247 = buffer.data(ksk + 1247);
    const auto *ksk_1248 = buffer.data(ksk + 1248);
    const auto *ksk_1249 = buffer.data(ksk + 1249);
    const auto *ksk_1252 = buffer.data(ksk + 1252);
    const auto *ksk_1253 = buffer.data(ksk + 1253);
    const auto *ksk_1254 = buffer.data(ksk + 1254);
    const auto *ksk_1255 = buffer.data(ksk + 1255);
    const auto *ksk_1256 = buffer.data(ksk + 1256);
    const auto *ksk_1257 = buffer.data(ksk + 1257);
    const auto *ksk_1258 = buffer.data(ksk + 1258);
    const auto *ksk_1259 = buffer.data(ksk + 1259);
    const auto *ksk_1260 = buffer.data(ksk + 1260);
    const auto *ksk_1265 = buffer.data(ksk + 1265);
    const auto *ksk_1269 = buffer.data(ksk + 1269);
    const auto *ksk_1274 = buffer.data(ksk + 1274);
    const auto *ksk_1280 = buffer.data(ksk + 1280);

    const auto *ksl1_1215 = buffer.data(ksl1 + 1215);
    const auto *ksl1_1220 = buffer.data(ksl1 + 1220);
    const auto *ksl1_1224 = buffer.data(ksl1 + 1224);
    const auto *ksl1_1229 = buffer.data(ksl1 + 1229);
    const auto *ksl1_1235 = buffer.data(ksl1 + 1235);
    const auto *ksl1_1242 = buffer.data(ksl1 + 1242);
    const auto *ksl1_1484 = buffer.data(ksl1 + 1484);
    const auto *ksl1_1485 = buffer.data(ksl1 + 1485);
    const auto *ksl1_1488 = buffer.data(ksl1 + 1488);
    const auto *ksl1_1490 = buffer.data(ksl1 + 1490);
    const auto *ksl1_1491 = buffer.data(ksl1 + 1491);
    const auto *ksl1_1494 = buffer.data(ksl1 + 1494);
    const auto *ksl1_1495 = buffer.data(ksl1 + 1495);
    const auto *ksl1_1497 = buffer.data(ksl1 + 1497);
    const auto *ksl1_1499 = buffer.data(ksl1 + 1499);
    const auto *ksl1_1500 = buffer.data(ksl1 + 1500);
    const auto *ksl1_1502 = buffer.data(ksl1 + 1502);
    const auto *ksl1_1503 = buffer.data(ksl1 + 1503);
    const auto *ksl1_1505 = buffer.data(ksl1 + 1505);
    const auto *ksl1_1506 = buffer.data(ksl1 + 1506);
    const auto *ksl1_1508 = buffer.data(ksl1 + 1508);
    const auto *ksl1_1509 = buffer.data(ksl1 + 1509);
    const auto *ksl1_1510 = buffer.data(ksl1 + 1510);
    const auto *ksl1_1512 = buffer.data(ksl1 + 1512);
    const auto *ksl1_1521 = buffer.data(ksl1 + 1521);
    const auto *ksl1_1523 = buffer.data(ksl1 + 1523);
    const auto *ksl1_1524 = buffer.data(ksl1 + 1524);
    const auto *ksl1_1525 = buffer.data(ksl1 + 1525);
    const auto *ksl1_1526 = buffer.data(ksl1 + 1526);
    const auto *ksl1_1527 = buffer.data(ksl1 + 1527);
    const auto *ksl1_1529 = buffer.data(ksl1 + 1529);
    const auto *ksl1_1533 = buffer.data(ksl1 + 1533);
    const auto *ksl1_1536 = buffer.data(ksl1 + 1536);
    const auto *ksl1_1540 = buffer.data(ksl1 + 1540);
    const auto *ksl1_1542 = buffer.data(ksl1 + 1542);
    const auto *ksl1_1545 = buffer.data(ksl1 + 1545);
    const auto *ksl1_1547 = buffer.data(ksl1 + 1547);
    const auto *ksl1_1548 = buffer.data(ksl1 + 1548);
    const auto *ksl1_1551 = buffer.data(ksl1 + 1551);
    const auto *ksl1_1553 = buffer.data(ksl1 + 1553);
    const auto *ksl1_1554 = buffer.data(ksl1 + 1554);
    const auto *ksl1_1555 = buffer.data(ksl1 + 1555);
    const auto *ksl1_1566 = buffer.data(ksl1 + 1566);
    const auto *ksl1_1568 = buffer.data(ksl1 + 1568);
    const auto *ksl1_1569 = buffer.data(ksl1 + 1569);
    const auto *ksl1_1570 = buffer.data(ksl1 + 1570);
    const auto *ksl1_1571 = buffer.data(ksl1 + 1571);
    const auto *ksl1_1572 = buffer.data(ksl1 + 1572);
    const auto *ksl1_1574 = buffer.data(ksl1 + 1574);
    const auto *ksl1_1575 = buffer.data(ksl1 + 1575);
    const auto *ksl1_1580 = buffer.data(ksl1 + 1580);
    const auto *ksl1_1584 = buffer.data(ksl1 + 1584);
    const auto *ksl1_1589 = buffer.data(ksl1 + 1589);
    const auto *ksl1_1595 = buffer.data(ksl1 + 1595);

    const auto *lsi0_980 = buffer.data(lsi0 + 980);
    const auto *lsi0_981 = buffer.data(lsi0 + 981);
    const auto *lsi0_982 = buffer.data(lsi0 + 982);
    const auto *lsi0_983 = buffer.data(lsi0 + 983);
    const auto *lsi0_984 = buffer.data(lsi0 + 984);
    const auto *lsi0_985 = buffer.data(lsi0 + 985);
    const auto *lsi0_986 = buffer.data(lsi0 + 986);
    const auto *lsi0_987 = buffer.data(lsi0 + 987);
    const auto *lsi0_988 = buffer.data(lsi0 + 988);
    const auto *lsi0_989 = buffer.data(lsi0 + 989);
    const auto *lsi0_990 = buffer.data(lsi0 + 990);
    const auto *lsi0_991 = buffer.data(lsi0 + 991);
    const auto *lsi0_992 = buffer.data(lsi0 + 992);
    const auto *lsi0_993 = buffer.data(lsi0 + 993);
    const auto *lsi0_994 = buffer.data(lsi0 + 994);

    const auto *lsi1_980 = buffer.data(lsi1 + 980);
    const auto *lsi1_981 = buffer.data(lsi1 + 981);
    const auto *lsi1_982 = buffer.data(lsi1 + 982);
    const auto *lsi1_983 = buffer.data(lsi1 + 983);
    const auto *lsi1_984 = buffer.data(lsi1 + 984);
    const auto *lsi1_985 = buffer.data(lsi1 + 985);
    const auto *lsi1_986 = buffer.data(lsi1 + 986);
    const auto *lsi1_987 = buffer.data(lsi1 + 987);
    const auto *lsi1_988 = buffer.data(lsi1 + 988);
    const auto *lsi1_989 = buffer.data(lsi1 + 989);
    const auto *lsi1_990 = buffer.data(lsi1 + 990);
    const auto *lsi1_991 = buffer.data(lsi1 + 991);
    const auto *lsi1_992 = buffer.data(lsi1 + 992);
    const auto *lsi1_993 = buffer.data(lsi1 + 993);
    const auto *lsi1_994 = buffer.data(lsi1 + 994);

    const auto *lsk_1188 = buffer.data(lsk + 1188);
    const auto *lsk_1190 = buffer.data(lsk + 1190);
    const auto *lsk_1191 = buffer.data(lsk + 1191);
    const auto *lsk_1193 = buffer.data(lsk + 1193);
    const auto *lsk_1194 = buffer.data(lsk + 1194);
    const auto *lsk_1197 = buffer.data(lsk + 1197);
    const auto *lsk_1198 = buffer.data(lsk + 1198);
    const auto *lsk_1202 = buffer.data(lsk + 1202);
    const auto *lsk_1203 = buffer.data(lsk + 1203);
    const auto *lsk_1208 = buffer.data(lsk + 1208);
    const auto *lsk_1216 = buffer.data(lsk + 1216);
    const auto *lsk_1217 = buffer.data(lsk + 1217);
    const auto *lsk_1218 = buffer.data(lsk + 1218);
    const auto *lsk_1219 = buffer.data(lsk + 1219);
    const auto *lsk_1220 = buffer.data(lsk + 1220);
    const auto *lsk_1221 = buffer.data(lsk + 1221);
    const auto *lsk_1222 = buffer.data(lsk + 1222);
    const auto *lsk_1223 = buffer.data(lsk + 1223);
    const auto *lsk_1224 = buffer.data(lsk + 1224);
    const auto *lsk_1226 = buffer.data(lsk + 1226);
    const auto *lsk_1227 = buffer.data(lsk + 1227);
    const auto *lsk_1229 = buffer.data(lsk + 1229);
    const auto *lsk_1230 = buffer.data(lsk + 1230);
    const auto *lsk_1233 = buffer.data(lsk + 1233);
    const auto *lsk_1234 = buffer.data(lsk + 1234);
    const auto *lsk_1238 = buffer.data(lsk + 1238);
    const auto *lsk_1239 = buffer.data(lsk + 1239);
    const auto *lsk_1244 = buffer.data(lsk + 1244);
    const auto *lsk_1252 = buffer.data(lsk + 1252);
    const auto *lsk_1253 = buffer.data(lsk + 1253);
    const auto *lsk_1254 = buffer.data(lsk + 1254);
    const auto *lsk_1255 = buffer.data(lsk + 1255);
    const auto *lsk_1256 = buffer.data(lsk + 1256);
    const auto *lsk_1257 = buffer.data(lsk + 1257);
    const auto *lsk_1258 = buffer.data(lsk + 1258);
    const auto *lsk_1259 = buffer.data(lsk + 1259);
    const auto *lsk_1260 = buffer.data(lsk + 1260);
    const auto *lsk_1261 = buffer.data(lsk + 1261);
    const auto *lsk_1262 = buffer.data(lsk + 1262);
    const auto *lsk_1263 = buffer.data(lsk + 1263);
    const auto *lsk_1264 = buffer.data(lsk + 1264);
    const auto *lsk_1265 = buffer.data(lsk + 1265);
    const auto *lsk_1266 = buffer.data(lsk + 1266);
    const auto *lsk_1267 = buffer.data(lsk + 1267);
    const auto *lsk_1268 = buffer.data(lsk + 1268);
    const auto *lsk_1269 = buffer.data(lsk + 1269);
    const auto *lsk_1270 = buffer.data(lsk + 1270);
    const auto *lsk_1271 = buffer.data(lsk + 1271);
    const auto *lsk_1272 = buffer.data(lsk + 1272);
    const auto *lsk_1273 = buffer.data(lsk + 1273);
    const auto *lsk_1274 = buffer.data(lsk + 1274);
    const auto *lsk_1275 = buffer.data(lsk + 1275);
    const auto *lsk_1276 = buffer.data(lsk + 1276);
    const auto *lsk_1277 = buffer.data(lsk + 1277);
    const auto *lsk_1278 = buffer.data(lsk + 1278);
    const auto *lsk_1279 = buffer.data(lsk + 1279);
    const auto *lsk_1280 = buffer.data(lsk + 1280);

#pragma omp simd aligned(t_1484, t_1485, t_1486, t_1487, pa_x, pc_x, pc_y, pc_z, ksl0_1484, \
                         ksl0_1485, ksk_900, ksk_936, ksk_1188, ksl1_1484, ksl1_1485, \
                         lsk_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1484[k] = pa_x[k] * ksl0_1484[k]
                    - f_14 * pc_x[k] * ksl1_1484[k];

        t_1485[k] = pa_x[k] * ksl0_1485[k]
                    + f_0 * ksk_1188[k]
                    - f_14 * pc_x[k] * ksl1_1485[k];

        t_1486[k] = f_16 * ksk_936[k]
                    + f_3 * pc_y[k] * lsk_1188[k];

        t_1487[k] = f_19 * ksk_900[k]
                    + f_3 * pc_z[k] * lsk_1188[k];
    }

#pragma omp simd aligned(t_1488, t_1489, t_1490, pa_x, pc_x, pc_y, ksl0_1488, ksl0_1490, \
                         ksk_938, ksk_1191, ksk_1193, ksl1_1488, ksl1_1490, \
                         lsk_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1488[k] = pa_x[k] * ksl0_1488[k]
                    + f_20 * ksk_1191[k]
                    - f_14 * pc_x[k] * ksl1_1488[k];

        t_1489[k] = f_16 * ksk_938[k]
                    + f_3 * pc_y[k] * lsk_1190[k];

        t_1490[k] = pa_x[k] * ksl0_1490[k]
                    + f_20 * ksk_1193[k]
                    - f_14 * pc_x[k] * ksl1_1490[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, pa_x, pc_x, pc_y, pc_z, ksl0_1491, ksk_903, \
                         ksk_941, ksk_1194, ksl1_1491, lsk_1191, \
                         lsk_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = pa_x[k] * ksl0_1491[k]
                    + f_19 * ksk_1194[k]
                    - f_14 * pc_x[k] * ksl1_1491[k];

        t_1492[k] = f_19 * ksk_903[k]
                    + f_3 * pc_z[k] * lsk_1191[k];

        t_1493[k] = f_16 * ksk_941[k]
                    + f_3 * pc_y[k] * lsk_1193[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, pa_x, pc_x, pc_z, ksl0_1494, ksl0_1495, \
                         ksk_906, ksk_1197, ksk_1198, ksl1_1494, ksl1_1495, \
                         lsk_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = pa_x[k] * ksl0_1494[k]
                    + f_19 * ksk_1197[k]
                    - f_14 * pc_x[k] * ksl1_1494[k];

        t_1495[k] = pa_x[k] * ksl0_1495[k]
                    + f_18 * ksk_1198[k]
                    - f_14 * pc_x[k] * ksl1_1495[k];

        t_1496[k] = f_19 * ksk_906[k]
                    + f_3 * pc_z[k] * lsk_1194[k];
    }

#pragma omp simd aligned(t_1497, t_1498, t_1499, pa_x, pc_x, pc_y, ksl0_1497, ksl0_1499, \
                         ksk_945, ksk_1200, ksk_1202, ksl1_1497, ksl1_1499, \
                         lsk_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1497[k] = pa_x[k] * ksl0_1497[k]
                    + f_18 * ksk_1200[k]
                    - f_14 * pc_x[k] * ksl1_1497[k];

        t_1498[k] = f_16 * ksk_945[k]
                    + f_3 * pc_y[k] * lsk_1197[k];

        t_1499[k] = pa_x[k] * ksl0_1499[k]
                    + f_18 * ksk_1202[k]
                    - f_14 * pc_x[k] * ksl1_1499[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, pa_x, pc_x, pc_z, ksl0_1500, ksl0_1502, \
                         ksk_910, ksk_1203, ksk_1205, ksl1_1500, ksl1_1502, \
                         lsk_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = pa_x[k] * ksl0_1500[k]
                    + f_17 * ksk_1203[k]
                    - f_14 * pc_x[k] * ksl1_1500[k];

        t_1501[k] = f_19 * ksk_910[k]
                    + f_3 * pc_z[k] * lsk_1198[k];

        t_1502[k] = pa_x[k] * ksl0_1502[k]
                    + f_17 * ksk_1205[k]
                    - f_14 * pc_x[k] * ksl1_1502[k];
    }

#pragma omp simd aligned(t_1503, t_1504, t_1505, pa_x, pc_x, pc_y, ksl0_1503, ksl0_1505, \
                         ksk_950, ksk_1206, ksk_1208, ksl1_1503, ksl1_1505, \
                         lsk_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1503[k] = pa_x[k] * ksl0_1503[k]
                    + f_17 * ksk_1206[k]
                    - f_14 * pc_x[k] * ksl1_1503[k];

        t_1504[k] = f_16 * ksk_950[k]
                    + f_3 * pc_y[k] * lsk_1202[k];

        t_1505[k] = pa_x[k] * ksl0_1505[k]
                    + f_17 * ksk_1208[k]
                    - f_14 * pc_x[k] * ksl1_1505[k];
    }

#pragma omp simd aligned(t_1506, t_1507, t_1508, pa_x, pc_x, pc_z, ksl0_1506, ksl0_1508, \
                         ksk_915, ksk_1209, ksk_1211, ksl1_1506, ksl1_1508, \
                         lsk_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1506[k] = pa_x[k] * ksl0_1506[k]
                    + f_16 * ksk_1209[k]
                    - f_14 * pc_x[k] * ksl1_1506[k];

        t_1507[k] = f_19 * ksk_915[k]
                    + f_3 * pc_z[k] * lsk_1203[k];

        t_1508[k] = pa_x[k] * ksl0_1508[k]
                    + f_16 * ksk_1211[k]
                    - f_14 * pc_x[k] * ksl1_1508[k];
    }

#pragma omp simd aligned(t_1509, t_1510, t_1511, pa_x, pc_x, pc_y, ksl0_1509, ksl0_1510, \
                         ksk_956, ksk_1212, ksk_1213, ksl1_1509, ksl1_1510, \
                         lsk_1208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1509[k] = pa_x[k] * ksl0_1509[k]
                    + f_16 * ksk_1212[k]
                    - f_14 * pc_x[k] * ksl1_1509[k];

        t_1510[k] = pa_x[k] * ksl0_1510[k]
                    + f_16 * ksk_1213[k]
                    - f_14 * pc_x[k] * ksl1_1510[k];

        t_1511[k] = f_16 * ksk_956[k]
                    + f_3 * pc_y[k] * lsk_1208[k];
    }

#pragma omp simd aligned(t_1512, t_1513, t_1514, t_1515, pa_x, pc_x, ksl0_1512, ksk_1215, \
                         ksk_1216, ksk_1217, ksk_1218, ksl1_1512, lsk_1216, lsk_1217, \
                         lsk_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1512[k] = pa_x[k] * ksl0_1512[k]
                    + f_16 * ksk_1215[k]
                    - f_14 * pc_x[k] * ksl1_1512[k];

        t_1513[k] = f_15 * ksk_1216[k]
                    + f_3 * pc_x[k] * lsk_1216[k];

        t_1514[k] = f_15 * ksk_1217[k]
                    + f_3 * pc_x[k] * lsk_1217[k];

        t_1515[k] = f_15 * ksk_1218[k]
                    + f_3 * pc_x[k] * lsk_1218[k];
    }

#pragma omp simd aligned(t_1516, t_1517, t_1518, t_1519, t_1520, pc_x, ksk_1219, ksk_1220, \
                         ksk_1221, ksk_1222, ksk_1223, lsk_1219, lsk_1220, lsk_1221, lsk_1222, \
                         lsk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1516[k] = f_15 * ksk_1219[k]
                    + f_3 * pc_x[k] * lsk_1219[k];

        t_1517[k] = f_15 * ksk_1220[k]
                    + f_3 * pc_x[k] * lsk_1220[k];

        t_1518[k] = f_15 * ksk_1221[k]
                    + f_3 * pc_x[k] * lsk_1221[k];

        t_1519[k] = f_15 * ksk_1222[k]
                    + f_3 * pc_x[k] * lsk_1222[k];

        t_1520[k] = f_15 * ksk_1223[k]
                    + f_3 * pc_x[k] * lsk_1223[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, t_1524, pa_x, pc_x, pc_z, ksl0_1521, \
                         ksl0_1523, ksl0_1524, ksk_928, ksl1_1521, ksl1_1523, ksl1_1524, \
                         lsk_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = pa_x[k] * ksl0_1521[k]
                    - f_14 * pc_x[k] * ksl1_1521[k];

        t_1522[k] = f_19 * ksk_928[k]
                    + f_3 * pc_z[k] * lsk_1216[k];

        t_1523[k] = pa_x[k] * ksl0_1523[k]
                    - f_14 * pc_x[k] * ksl1_1523[k];

        t_1524[k] = pa_x[k] * ksl0_1524[k]
                    - f_14 * pc_x[k] * ksl1_1524[k];
    }

#pragma omp simd aligned(t_1525, t_1526, t_1527, t_1528, pa_x, pc_x, pc_y, ksl0_1525, \
                         ksl0_1526, ksl0_1527, ksk_971, ksl1_1525, ksl1_1526, ksl1_1527, \
                         lsk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1525[k] = pa_x[k] * ksl0_1525[k]
                    - f_14 * pc_x[k] * ksl1_1525[k];

        t_1526[k] = pa_x[k] * ksl0_1526[k]
                    - f_14 * pc_x[k] * ksl1_1526[k];

        t_1527[k] = pa_x[k] * ksl0_1527[k]
                    - f_14 * pc_x[k] * ksl1_1527[k];

        t_1528[k] = f_16 * ksk_971[k]
                    + f_3 * pc_y[k] * lsk_1223[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, t_1532, pa_x, pa_y, pc_x, pc_y, pc_z, \
                         ksl0_1215, ksl0_1529, ksk_936, ksk_972, ksl1_1215, ksl1_1529, \
                         lsk_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = pa_x[k] * ksl0_1529[k]
                    - f_14 * pc_x[k] * ksl1_1529[k];

        t_1530[k] = pa_y[k] * ksl0_1215[k]
                    - f_14 * pc_y[k] * ksl1_1215[k];

        t_1531[k] = f_15 * ksk_972[k]
                    + f_3 * pc_y[k] * lsk_1224[k];

        t_1532[k] = f_20 * ksk_936[k]
                    + f_3 * pc_z[k] * lsk_1224[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, pa_x, pa_y, pc_x, pc_y, ksl0_1220, ksl0_1533, \
                         ksk_974, ksk_1227, ksl1_1220, ksl1_1533, \
                         lsk_1226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = pa_x[k] * ksl0_1533[k]
                    + f_20 * ksk_1227[k]
                    - f_14 * pc_x[k] * ksl1_1533[k];

        t_1534[k] = f_15 * ksk_974[k]
                    + f_3 * pc_y[k] * lsk_1226[k];

        t_1535[k] = pa_y[k] * ksl0_1220[k]
                    - f_14 * pc_y[k] * ksl1_1220[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, pa_x, pc_x, pc_y, pc_z, ksl0_1536, ksk_939, \
                         ksk_977, ksk_1230, ksl1_1536, lsk_1227, \
                         lsk_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = pa_x[k] * ksl0_1536[k]
                    + f_19 * ksk_1230[k]
                    - f_14 * pc_x[k] * ksl1_1536[k];

        t_1537[k] = f_20 * ksk_939[k]
                    + f_3 * pc_z[k] * lsk_1227[k];

        t_1538[k] = f_15 * ksk_977[k]
                    + f_3 * pc_y[k] * lsk_1229[k];
    }

#pragma omp simd aligned(t_1539, t_1540, t_1541, pa_x, pa_y, pc_x, pc_y, pc_z, ksl0_1224, \
                         ksl0_1540, ksk_942, ksk_1234, ksl1_1224, ksl1_1540, \
                         lsk_1230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1539[k] = pa_y[k] * ksl0_1224[k]
                    - f_14 * pc_y[k] * ksl1_1224[k];

        t_1540[k] = pa_x[k] * ksl0_1540[k]
                    + f_18 * ksk_1234[k]
                    - f_14 * pc_x[k] * ksl1_1540[k];

        t_1541[k] = f_20 * ksk_942[k]
                    + f_3 * pc_z[k] * lsk_1230[k];
    }

#pragma omp simd aligned(t_1542, t_1543, t_1544, pa_x, pa_y, pc_x, pc_y, ksl0_1229, ksl0_1542, \
                         ksk_981, ksk_1236, ksl1_1229, ksl1_1542, \
                         lsk_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1542[k] = pa_x[k] * ksl0_1542[k]
                    + f_18 * ksk_1236[k]
                    - f_14 * pc_x[k] * ksl1_1542[k];

        t_1543[k] = f_15 * ksk_981[k]
                    + f_3 * pc_y[k] * lsk_1233[k];

        t_1544[k] = pa_y[k] * ksl0_1229[k]
                    - f_14 * pc_y[k] * ksl1_1229[k];
    }

#pragma omp simd aligned(t_1545, t_1546, t_1547, pa_x, pc_x, pc_z, ksl0_1545, ksl0_1547, \
                         ksk_946, ksk_1239, ksk_1241, ksl1_1545, ksl1_1547, \
                         lsk_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1545[k] = pa_x[k] * ksl0_1545[k]
                    + f_17 * ksk_1239[k]
                    - f_14 * pc_x[k] * ksl1_1545[k];

        t_1546[k] = f_20 * ksk_946[k]
                    + f_3 * pc_z[k] * lsk_1234[k];

        t_1547[k] = pa_x[k] * ksl0_1547[k]
                    + f_17 * ksk_1241[k]
                    - f_14 * pc_x[k] * ksl1_1547[k];
    }

#pragma omp simd aligned(t_1548, t_1549, t_1550, pa_x, pa_y, pc_x, pc_y, ksl0_1235, ksl0_1548, \
                         ksk_986, ksk_1242, ksl1_1235, ksl1_1548, \
                         lsk_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1548[k] = pa_x[k] * ksl0_1548[k]
                    + f_17 * ksk_1242[k]
                    - f_14 * pc_x[k] * ksl1_1548[k];

        t_1549[k] = f_15 * ksk_986[k]
                    + f_3 * pc_y[k] * lsk_1238[k];

        t_1550[k] = pa_y[k] * ksl0_1235[k]
                    - f_14 * pc_y[k] * ksl1_1235[k];
    }

#pragma omp simd aligned(t_1551, t_1552, t_1553, pa_x, pc_x, pc_z, ksl0_1551, ksl0_1553, \
                         ksk_951, ksk_1245, ksk_1247, ksl1_1551, ksl1_1553, \
                         lsk_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1551[k] = pa_x[k] * ksl0_1551[k]
                    + f_16 * ksk_1245[k]
                    - f_14 * pc_x[k] * ksl1_1551[k];

        t_1552[k] = f_20 * ksk_951[k]
                    + f_3 * pc_z[k] * lsk_1239[k];

        t_1553[k] = pa_x[k] * ksl0_1553[k]
                    + f_16 * ksk_1247[k]
                    - f_14 * pc_x[k] * ksl1_1553[k];
    }

#pragma omp simd aligned(t_1554, t_1555, t_1556, pa_x, pc_x, pc_y, ksl0_1554, ksl0_1555, \
                         ksk_992, ksk_1248, ksk_1249, ksl1_1554, ksl1_1555, \
                         lsk_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1554[k] = pa_x[k] * ksl0_1554[k]
                    + f_16 * ksk_1248[k]
                    - f_14 * pc_x[k] * ksl1_1554[k];

        t_1555[k] = pa_x[k] * ksl0_1555[k]
                    + f_16 * ksk_1249[k]
                    - f_14 * pc_x[k] * ksl1_1555[k];

        t_1556[k] = f_15 * ksk_992[k]
                    + f_3 * pc_y[k] * lsk_1244[k];
    }

#pragma omp simd aligned(t_1557, t_1558, t_1559, t_1560, pa_y, pc_x, pc_y, ksl0_1242, \
                         ksk_1252, ksk_1253, ksk_1254, ksl1_1242, lsk_1252, lsk_1253, \
                         lsk_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1557[k] = pa_y[k] * ksl0_1242[k]
                    - f_14 * pc_y[k] * ksl1_1242[k];

        t_1558[k] = f_15 * ksk_1252[k]
                    + f_3 * pc_x[k] * lsk_1252[k];

        t_1559[k] = f_15 * ksk_1253[k]
                    + f_3 * pc_x[k] * lsk_1253[k];

        t_1560[k] = f_15 * ksk_1254[k]
                    + f_3 * pc_x[k] * lsk_1254[k];
    }

#pragma omp simd aligned(t_1561, t_1562, t_1563, t_1564, t_1565, pc_x, ksk_1255, ksk_1256, \
                         ksk_1257, ksk_1258, ksk_1259, lsk_1255, lsk_1256, lsk_1257, lsk_1258, \
                         lsk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1561[k] = f_15 * ksk_1255[k]
                    + f_3 * pc_x[k] * lsk_1255[k];

        t_1562[k] = f_15 * ksk_1256[k]
                    + f_3 * pc_x[k] * lsk_1256[k];

        t_1563[k] = f_15 * ksk_1257[k]
                    + f_3 * pc_x[k] * lsk_1257[k];

        t_1564[k] = f_15 * ksk_1258[k]
                    + f_3 * pc_x[k] * lsk_1258[k];

        t_1565[k] = f_15 * ksk_1259[k]
                    + f_3 * pc_x[k] * lsk_1259[k];
    }

#pragma omp simd aligned(t_1566, t_1567, t_1568, t_1569, pa_x, pc_x, pc_z, ksl0_1566, \
                         ksl0_1568, ksl0_1569, ksk_964, ksl1_1566, ksl1_1568, ksl1_1569, \
                         lsk_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1566[k] = pa_x[k] * ksl0_1566[k]
                    - f_14 * pc_x[k] * ksl1_1566[k];

        t_1567[k] = f_20 * ksk_964[k]
                    + f_3 * pc_z[k] * lsk_1252[k];

        t_1568[k] = pa_x[k] * ksl0_1568[k]
                    - f_14 * pc_x[k] * ksl1_1568[k];

        t_1569[k] = pa_x[k] * ksl0_1569[k]
                    - f_14 * pc_x[k] * ksl1_1569[k];
    }

#pragma omp simd aligned(t_1570, t_1571, t_1572, t_1573, pa_x, pc_x, pc_y, ksl0_1570, \
                         ksl0_1571, ksl0_1572, ksk_1007, ksl1_1570, ksl1_1571, ksl1_1572, \
                         lsk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1570[k] = pa_x[k] * ksl0_1570[k]
                    - f_14 * pc_x[k] * ksl1_1570[k];

        t_1571[k] = pa_x[k] * ksl0_1571[k]
                    - f_14 * pc_x[k] * ksl1_1571[k];

        t_1572[k] = pa_x[k] * ksl0_1572[k]
                    - f_14 * pc_x[k] * ksl1_1572[k];

        t_1573[k] = f_15 * ksk_1007[k]
                    + f_3 * pc_y[k] * lsk_1259[k];
    }

#pragma omp simd aligned(t_1574, t_1575, t_1576, t_1577, pa_x, pc_x, pc_y, pc_z, ksl0_1574, \
                         ksl0_1575, ksk_972, ksk_1260, ksl1_1574, ksl1_1575, \
                         lsk_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = pa_x[k] * ksl0_1574[k]
                    - f_14 * pc_x[k] * ksl1_1574[k];

        t_1575[k] = pa_x[k] * ksl0_1575[k]
                    + f_0 * ksk_1260[k]
                    - f_14 * pc_x[k] * ksl1_1575[k];

        t_1576[k] = f_3 * pc_y[k] * lsk_1260[k];

        t_1577[k] = f_21 * ksk_972[k]
                    + f_3 * pc_z[k] * lsk_1260[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pa_x, pc_x, pc_y, ksl0_1580, ksk_1265, \
                         ksl1_1580, lsi0_980, lsi1_980, lsk_1261, \
                         lsk_1262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = f_4 * lsi0_980[k]
                    - f_5 * lsi1_980[k]
                    + f_3 * pc_y[k] * lsk_1261[k];

        t_1579[k] = f_3 * pc_y[k] * lsk_1262[k];

        t_1580[k] = pa_x[k] * ksl0_1580[k]
                    + f_20 * ksk_1265[k]
                    - f_14 * pc_x[k] * ksl1_1580[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, pc_y, lsi0_981, lsi0_982, lsi1_981, lsi1_982, \
                         lsk_1263, lsk_1264, lsk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_6 * lsi0_981[k]
                    - f_7 * lsi1_981[k]
                    + f_3 * pc_y[k] * lsk_1263[k];

        t_1582[k] = f_4 * lsi0_982[k]
                    - f_5 * lsi1_982[k]
                    + f_3 * pc_y[k] * lsk_1264[k];

        t_1583[k] = f_3 * pc_y[k] * lsk_1265[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, pa_x, pc_x, pc_y, ksl0_1584, ksk_1269, \
                         ksl1_1584, lsi0_983, lsi0_984, lsi1_983, lsi1_984, lsk_1266, \
                         lsk_1267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = pa_x[k] * ksl0_1584[k]
                    + f_19 * ksk_1269[k]
                    - f_14 * pc_x[k] * ksl1_1584[k];

        t_1585[k] = f_8 * lsi0_983[k]
                    - f_9 * lsi1_983[k]
                    + f_3 * pc_y[k] * lsk_1266[k];

        t_1586[k] = f_6 * lsi0_984[k]
                    - f_7 * lsi1_984[k]
                    + f_3 * pc_y[k] * lsk_1267[k];
    }

#pragma omp simd aligned(t_1587, t_1588, t_1589, pa_x, pc_x, pc_y, ksl0_1589, ksk_1274, \
                         ksl1_1589, lsi0_985, lsi1_985, lsk_1268, \
                         lsk_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1587[k] = f_4 * lsi0_985[k]
                    - f_5 * lsi1_985[k]
                    + f_3 * pc_y[k] * lsk_1268[k];

        t_1588[k] = f_3 * pc_y[k] * lsk_1269[k];

        t_1589[k] = pa_x[k] * ksl0_1589[k]
                    + f_18 * ksk_1274[k]
                    - f_14 * pc_x[k] * ksl1_1589[k];
    }

#pragma omp simd aligned(t_1590, t_1591, t_1592, pc_y, lsi0_986, lsi0_987, lsi0_988, lsi1_986, \
                         lsi1_987, lsi1_988, lsk_1270, lsk_1271, \
                         lsk_1272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1590[k] = f_10 * lsi0_986[k]
                    - f_11 * lsi1_986[k]
                    + f_3 * pc_y[k] * lsk_1270[k];

        t_1591[k] = f_8 * lsi0_987[k]
                    - f_9 * lsi1_987[k]
                    + f_3 * pc_y[k] * lsk_1271[k];

        t_1592[k] = f_6 * lsi0_988[k]
                    - f_7 * lsi1_988[k]
                    + f_3 * pc_y[k] * lsk_1272[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, pa_x, pc_x, pc_y, ksl0_1595, ksk_1280, \
                         ksl1_1595, lsi0_989, lsi1_989, lsk_1273, \
                         lsk_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = f_4 * lsi0_989[k]
                    - f_5 * lsi1_989[k]
                    + f_3 * pc_y[k] * lsk_1273[k];

        t_1594[k] = f_3 * pc_y[k] * lsk_1274[k];

        t_1595[k] = pa_x[k] * ksl0_1595[k]
                    + f_17 * ksk_1280[k]
                    - f_14 * pc_x[k] * ksl1_1595[k];
    }

#pragma omp simd aligned(t_1596, t_1597, t_1598, pc_y, lsi0_990, lsi0_991, lsi0_992, lsi1_990, \
                         lsi1_991, lsi1_992, lsk_1275, lsk_1276, \
                         lsk_1277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = f_12 * lsi0_990[k]
                    - f_13 * lsi1_990[k]
                    + f_3 * pc_y[k] * lsk_1275[k];

        t_1597[k] = f_10 * lsi0_991[k]
                    - f_11 * lsi1_991[k]
                    + f_3 * pc_y[k] * lsk_1276[k];

        t_1598[k] = f_8 * lsi0_992[k]
                    - f_9 * lsi1_992[k]
                    + f_3 * pc_y[k] * lsk_1277[k];
    }

#pragma omp simd aligned(t_1599, t_1600, t_1601, pc_y, lsi0_993, lsi0_994, lsi1_993, lsi1_994, \
                         lsk_1278, lsk_1279, lsk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1599[k] = f_6 * lsi0_993[k]
                    - f_7 * lsi1_993[k]
                    + f_3 * pc_y[k] * lsk_1278[k];

        t_1600[k] = f_4 * lsi0_994[k]
                    - f_5 * lsi1_994[k]
                    + f_3 * pc_y[k] * lsk_1279[k];

        t_1601[k] = f_3 * pc_y[k] * lsk_1280[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksl0,
                                                           const size_t ksk, const size_t ksl1,
                                                           const size_t lsi0, const size_t lsi1,
                                                           const size_t lsk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_1260 = buffer.data(ksl0 + 1260);
    const auto *ksl0_1261 = buffer.data(ksl0 + 1261);
    const auto *ksl0_1263 = buffer.data(ksl0 + 1263);
    const auto *ksl0_1266 = buffer.data(ksl0 + 1266);
    const auto *ksl0_1270 = buffer.data(ksl0 + 1270);
    const auto *ksl0_1275 = buffer.data(ksl0 + 1275);
    const auto *ksl0_1281 = buffer.data(ksl0 + 1281);
    const auto *ksl0_1296 = buffer.data(ksl0 + 1296);
    const auto *ksl0_1298 = buffer.data(ksl0 + 1298);
    const auto *ksl0_1299 = buffer.data(ksl0 + 1299);
    const auto *ksl0_1300 = buffer.data(ksl0 + 1300);
    const auto *ksl0_1301 = buffer.data(ksl0 + 1301);
    const auto *ksl0_1302 = buffer.data(ksl0 + 1302);
    const auto *ksl0_1602 = buffer.data(ksl0 + 1602);
    const auto *ksl0_1611 = buffer.data(ksl0 + 1611);
    const auto *ksl0_1612 = buffer.data(ksl0 + 1612);
    const auto *ksl0_1613 = buffer.data(ksl0 + 1613);
    const auto *ksl0_1614 = buffer.data(ksl0 + 1614);
    const auto *ksl0_1615 = buffer.data(ksl0 + 1615);
    const auto *ksl0_1616 = buffer.data(ksl0 + 1616);
    const auto *ksl0_1617 = buffer.data(ksl0 + 1617);
    const auto *ksl0_1619 = buffer.data(ksl0 + 1619);

    const auto *ksk_1036 = buffer.data(ksk + 1036);
    const auto *ksk_1037 = buffer.data(ksk + 1037);
    const auto *ksk_1038 = buffer.data(ksk + 1038);
    const auto *ksk_1039 = buffer.data(ksk + 1039);
    const auto *ksk_1040 = buffer.data(ksk + 1040);
    const auto *ksk_1041 = buffer.data(ksk + 1041);
    const auto *ksk_1043 = buffer.data(ksk + 1043);
    const auto *ksk_1079 = buffer.data(ksk + 1079);
    const auto *ksk_1287 = buffer.data(ksk + 1287);
    const auto *ksk_1288 = buffer.data(ksk + 1288);
    const auto *ksk_1289 = buffer.data(ksk + 1289);
    const auto *ksk_1290 = buffer.data(ksk + 1290);
    const auto *ksk_1291 = buffer.data(ksk + 1291);
    const auto *ksk_1292 = buffer.data(ksk + 1292);
    const auto *ksk_1293 = buffer.data(ksk + 1293);
    const auto *ksk_1295 = buffer.data(ksk + 1295);

    const auto *ksl1_1260 = buffer.data(ksl1 + 1260);
    const auto *ksl1_1261 = buffer.data(ksl1 + 1261);
    const auto *ksl1_1263 = buffer.data(ksl1 + 1263);
    const auto *ksl1_1266 = buffer.data(ksl1 + 1266);
    const auto *ksl1_1270 = buffer.data(ksl1 + 1270);
    const auto *ksl1_1275 = buffer.data(ksl1 + 1275);
    const auto *ksl1_1281 = buffer.data(ksl1 + 1281);
    const auto *ksl1_1296 = buffer.data(ksl1 + 1296);
    const auto *ksl1_1298 = buffer.data(ksl1 + 1298);
    const auto *ksl1_1299 = buffer.data(ksl1 + 1299);
    const auto *ksl1_1300 = buffer.data(ksl1 + 1300);
    const auto *ksl1_1301 = buffer.data(ksl1 + 1301);
    const auto *ksl1_1302 = buffer.data(ksl1 + 1302);
    const auto *ksl1_1602 = buffer.data(ksl1 + 1602);
    const auto *ksl1_1611 = buffer.data(ksl1 + 1611);
    const auto *ksl1_1612 = buffer.data(ksl1 + 1612);
    const auto *ksl1_1613 = buffer.data(ksl1 + 1613);
    const auto *ksl1_1614 = buffer.data(ksl1 + 1614);
    const auto *ksl1_1615 = buffer.data(ksl1 + 1615);
    const auto *ksl1_1616 = buffer.data(ksl1 + 1616);
    const auto *ksl1_1617 = buffer.data(ksl1 + 1617);
    const auto *ksl1_1619 = buffer.data(ksl1 + 1619);

    const auto *lsi0_1008 = buffer.data(lsi0 + 1008);
    const auto *lsi0_1009 = buffer.data(lsi0 + 1009);
    const auto *lsi0_1011 = buffer.data(lsi0 + 1011);
    const auto *lsi0_1013 = buffer.data(lsi0 + 1013);
    const auto *lsi0_1014 = buffer.data(lsi0 + 1014);
    const auto *lsi0_1016 = buffer.data(lsi0 + 1016);
    const auto *lsi0_1017 = buffer.data(lsi0 + 1017);
    const auto *lsi0_1018 = buffer.data(lsi0 + 1018);
    const auto *lsi0_1020 = buffer.data(lsi0 + 1020);
    const auto *lsi0_1021 = buffer.data(lsi0 + 1021);
    const auto *lsi0_1022 = buffer.data(lsi0 + 1022);
    const auto *lsi0_1023 = buffer.data(lsi0 + 1023);
    const auto *lsi0_1025 = buffer.data(lsi0 + 1025);
    const auto *lsi0_1026 = buffer.data(lsi0 + 1026);
    const auto *lsi0_1027 = buffer.data(lsi0 + 1027);
    const auto *lsi0_1028 = buffer.data(lsi0 + 1028);
    const auto *lsi0_1029 = buffer.data(lsi0 + 1029);
    const auto *lsi0_1030 = buffer.data(lsi0 + 1030);
    const auto *lsi0_1031 = buffer.data(lsi0 + 1031);
    const auto *lsi0_1032 = buffer.data(lsi0 + 1032);
    const auto *lsi0_1033 = buffer.data(lsi0 + 1033);
    const auto *lsi0_1034 = buffer.data(lsi0 + 1034);
    const auto *lsi0_1035 = buffer.data(lsi0 + 1035);
    const auto *lsi0_1038 = buffer.data(lsi0 + 1038);
    const auto *lsi0_1040 = buffer.data(lsi0 + 1040);
    const auto *lsi0_1041 = buffer.data(lsi0 + 1041);
    const auto *lsi0_1043 = buffer.data(lsi0 + 1043);
    const auto *lsi0_1044 = buffer.data(lsi0 + 1044);
    const auto *lsi0_1045 = buffer.data(lsi0 + 1045);
    const auto *lsi0_1047 = buffer.data(lsi0 + 1047);
    const auto *lsi0_1048 = buffer.data(lsi0 + 1048);
    const auto *lsi0_1049 = buffer.data(lsi0 + 1049);
    const auto *lsi0_1050 = buffer.data(lsi0 + 1050);
    const auto *lsi0_1052 = buffer.data(lsi0 + 1052);
    const auto *lsi0_1053 = buffer.data(lsi0 + 1053);
    const auto *lsi0_1054 = buffer.data(lsi0 + 1054);
    const auto *lsi0_1055 = buffer.data(lsi0 + 1055);
    const auto *lsi0_1056 = buffer.data(lsi0 + 1056);
    const auto *lsi0_1058 = buffer.data(lsi0 + 1058);
    const auto *lsi0_1059 = buffer.data(lsi0 + 1059);
    const auto *lsi0_1060 = buffer.data(lsi0 + 1060);
    const auto *lsi0_1061 = buffer.data(lsi0 + 1061);
    const auto *lsi0_1062 = buffer.data(lsi0 + 1062);
    const auto *lsi0_1063 = buffer.data(lsi0 + 1063);
    const auto *lsi0_1064 = buffer.data(lsi0 + 1064);
    const auto *lsi0_1065 = buffer.data(lsi0 + 1065);
    const auto *lsi0_1066 = buffer.data(lsi0 + 1066);
    const auto *lsi0_1067 = buffer.data(lsi0 + 1067);
    const auto *lsi0_1068 = buffer.data(lsi0 + 1068);
    const auto *lsi0_1069 = buffer.data(lsi0 + 1069);
    const auto *lsi0_1070 = buffer.data(lsi0 + 1070);
    const auto *lsi0_1071 = buffer.data(lsi0 + 1071);
    const auto *lsi0_1072 = buffer.data(lsi0 + 1072);
    const auto *lsi0_1073 = buffer.data(lsi0 + 1073);
    const auto *lsi0_1074 = buffer.data(lsi0 + 1074);
    const auto *lsi0_1075 = buffer.data(lsi0 + 1075);
    const auto *lsi0_1076 = buffer.data(lsi0 + 1076);

    const auto *lsi1_1008 = buffer.data(lsi1 + 1008);
    const auto *lsi1_1009 = buffer.data(lsi1 + 1009);
    const auto *lsi1_1011 = buffer.data(lsi1 + 1011);
    const auto *lsi1_1013 = buffer.data(lsi1 + 1013);
    const auto *lsi1_1014 = buffer.data(lsi1 + 1014);
    const auto *lsi1_1016 = buffer.data(lsi1 + 1016);
    const auto *lsi1_1017 = buffer.data(lsi1 + 1017);
    const auto *lsi1_1018 = buffer.data(lsi1 + 1018);
    const auto *lsi1_1020 = buffer.data(lsi1 + 1020);
    const auto *lsi1_1021 = buffer.data(lsi1 + 1021);
    const auto *lsi1_1022 = buffer.data(lsi1 + 1022);
    const auto *lsi1_1023 = buffer.data(lsi1 + 1023);
    const auto *lsi1_1025 = buffer.data(lsi1 + 1025);
    const auto *lsi1_1026 = buffer.data(lsi1 + 1026);
    const auto *lsi1_1027 = buffer.data(lsi1 + 1027);
    const auto *lsi1_1028 = buffer.data(lsi1 + 1028);
    const auto *lsi1_1029 = buffer.data(lsi1 + 1029);
    const auto *lsi1_1030 = buffer.data(lsi1 + 1030);
    const auto *lsi1_1031 = buffer.data(lsi1 + 1031);
    const auto *lsi1_1032 = buffer.data(lsi1 + 1032);
    const auto *lsi1_1033 = buffer.data(lsi1 + 1033);
    const auto *lsi1_1034 = buffer.data(lsi1 + 1034);
    const auto *lsi1_1035 = buffer.data(lsi1 + 1035);
    const auto *lsi1_1038 = buffer.data(lsi1 + 1038);
    const auto *lsi1_1040 = buffer.data(lsi1 + 1040);
    const auto *lsi1_1041 = buffer.data(lsi1 + 1041);
    const auto *lsi1_1043 = buffer.data(lsi1 + 1043);
    const auto *lsi1_1044 = buffer.data(lsi1 + 1044);
    const auto *lsi1_1045 = buffer.data(lsi1 + 1045);
    const auto *lsi1_1047 = buffer.data(lsi1 + 1047);
    const auto *lsi1_1048 = buffer.data(lsi1 + 1048);
    const auto *lsi1_1049 = buffer.data(lsi1 + 1049);
    const auto *lsi1_1050 = buffer.data(lsi1 + 1050);
    const auto *lsi1_1052 = buffer.data(lsi1 + 1052);
    const auto *lsi1_1053 = buffer.data(lsi1 + 1053);
    const auto *lsi1_1054 = buffer.data(lsi1 + 1054);
    const auto *lsi1_1055 = buffer.data(lsi1 + 1055);
    const auto *lsi1_1056 = buffer.data(lsi1 + 1056);
    const auto *lsi1_1058 = buffer.data(lsi1 + 1058);
    const auto *lsi1_1059 = buffer.data(lsi1 + 1059);
    const auto *lsi1_1060 = buffer.data(lsi1 + 1060);
    const auto *lsi1_1061 = buffer.data(lsi1 + 1061);
    const auto *lsi1_1062 = buffer.data(lsi1 + 1062);
    const auto *lsi1_1063 = buffer.data(lsi1 + 1063);
    const auto *lsi1_1064 = buffer.data(lsi1 + 1064);
    const auto *lsi1_1065 = buffer.data(lsi1 + 1065);
    const auto *lsi1_1066 = buffer.data(lsi1 + 1066);
    const auto *lsi1_1067 = buffer.data(lsi1 + 1067);
    const auto *lsi1_1068 = buffer.data(lsi1 + 1068);
    const auto *lsi1_1069 = buffer.data(lsi1 + 1069);
    const auto *lsi1_1070 = buffer.data(lsi1 + 1070);
    const auto *lsi1_1071 = buffer.data(lsi1 + 1071);
    const auto *lsi1_1072 = buffer.data(lsi1 + 1072);
    const auto *lsi1_1073 = buffer.data(lsi1 + 1073);
    const auto *lsi1_1074 = buffer.data(lsi1 + 1074);
    const auto *lsi1_1075 = buffer.data(lsi1 + 1075);
    const auto *lsi1_1076 = buffer.data(lsi1 + 1076);

    const auto *lsk_1287 = buffer.data(lsk + 1287);
    const auto *lsk_1288 = buffer.data(lsk + 1288);
    const auto *lsk_1289 = buffer.data(lsk + 1289);
    const auto *lsk_1290 = buffer.data(lsk + 1290);
    const auto *lsk_1291 = buffer.data(lsk + 1291);
    const auto *lsk_1292 = buffer.data(lsk + 1292);
    const auto *lsk_1293 = buffer.data(lsk + 1293);
    const auto *lsk_1295 = buffer.data(lsk + 1295);
    const auto *lsk_1296 = buffer.data(lsk + 1296);
    const auto *lsk_1297 = buffer.data(lsk + 1297);
    const auto *lsk_1299 = buffer.data(lsk + 1299);
    const auto *lsk_1301 = buffer.data(lsk + 1301);
    const auto *lsk_1302 = buffer.data(lsk + 1302);
    const auto *lsk_1304 = buffer.data(lsk + 1304);
    const auto *lsk_1305 = buffer.data(lsk + 1305);
    const auto *lsk_1306 = buffer.data(lsk + 1306);
    const auto *lsk_1308 = buffer.data(lsk + 1308);
    const auto *lsk_1309 = buffer.data(lsk + 1309);
    const auto *lsk_1310 = buffer.data(lsk + 1310);
    const auto *lsk_1311 = buffer.data(lsk + 1311);
    const auto *lsk_1313 = buffer.data(lsk + 1313);
    const auto *lsk_1314 = buffer.data(lsk + 1314);
    const auto *lsk_1315 = buffer.data(lsk + 1315);
    const auto *lsk_1316 = buffer.data(lsk + 1316);
    const auto *lsk_1317 = buffer.data(lsk + 1317);
    const auto *lsk_1319 = buffer.data(lsk + 1319);
    const auto *lsk_1320 = buffer.data(lsk + 1320);
    const auto *lsk_1321 = buffer.data(lsk + 1321);
    const auto *lsk_1322 = buffer.data(lsk + 1322);
    const auto *lsk_1323 = buffer.data(lsk + 1323);
    const auto *lsk_1324 = buffer.data(lsk + 1324);
    const auto *lsk_1325 = buffer.data(lsk + 1325);
    const auto *lsk_1326 = buffer.data(lsk + 1326);
    const auto *lsk_1327 = buffer.data(lsk + 1327);
    const auto *lsk_1328 = buffer.data(lsk + 1328);
    const auto *lsk_1329 = buffer.data(lsk + 1329);
    const auto *lsk_1330 = buffer.data(lsk + 1330);
    const auto *lsk_1331 = buffer.data(lsk + 1331);
    const auto *lsk_1334 = buffer.data(lsk + 1334);
    const auto *lsk_1336 = buffer.data(lsk + 1336);
    const auto *lsk_1337 = buffer.data(lsk + 1337);
    const auto *lsk_1339 = buffer.data(lsk + 1339);
    const auto *lsk_1340 = buffer.data(lsk + 1340);
    const auto *lsk_1341 = buffer.data(lsk + 1341);
    const auto *lsk_1343 = buffer.data(lsk + 1343);
    const auto *lsk_1344 = buffer.data(lsk + 1344);
    const auto *lsk_1345 = buffer.data(lsk + 1345);
    const auto *lsk_1346 = buffer.data(lsk + 1346);
    const auto *lsk_1348 = buffer.data(lsk + 1348);
    const auto *lsk_1349 = buffer.data(lsk + 1349);
    const auto *lsk_1350 = buffer.data(lsk + 1350);
    const auto *lsk_1351 = buffer.data(lsk + 1351);
    const auto *lsk_1352 = buffer.data(lsk + 1352);
    const auto *lsk_1354 = buffer.data(lsk + 1354);
    const auto *lsk_1355 = buffer.data(lsk + 1355);
    const auto *lsk_1356 = buffer.data(lsk + 1356);
    const auto *lsk_1357 = buffer.data(lsk + 1357);
    const auto *lsk_1358 = buffer.data(lsk + 1358);
    const auto *lsk_1359 = buffer.data(lsk + 1359);
    const auto *lsk_1360 = buffer.data(lsk + 1360);
    const auto *lsk_1361 = buffer.data(lsk + 1361);
    const auto *lsk_1362 = buffer.data(lsk + 1362);
    const auto *lsk_1363 = buffer.data(lsk + 1363);
    const auto *lsk_1364 = buffer.data(lsk + 1364);
    const auto *lsk_1365 = buffer.data(lsk + 1365);
    const auto *lsk_1366 = buffer.data(lsk + 1366);
    const auto *lsk_1367 = buffer.data(lsk + 1367);
    const auto *lsk_1368 = buffer.data(lsk + 1368);
    const auto *lsk_1369 = buffer.data(lsk + 1369);
    const auto *lsk_1370 = buffer.data(lsk + 1370);
    const auto *lsk_1371 = buffer.data(lsk + 1371);
    const auto *lsk_1372 = buffer.data(lsk + 1372);
    const auto *lsk_1373 = buffer.data(lsk + 1373);
    const auto *lsk_1374 = buffer.data(lsk + 1374);
    const auto *lsk_1375 = buffer.data(lsk + 1375);
    const auto *lsk_1376 = buffer.data(lsk + 1376);
    const auto *lsk_1377 = buffer.data(lsk + 1377);
    const auto *lsk_1378 = buffer.data(lsk + 1378);
    const auto *lsk_1379 = buffer.data(lsk + 1379);
    const auto *lsk_1380 = buffer.data(lsk + 1380);

#pragma omp simd aligned(t_1602, t_1603, t_1604, t_1605, pa_x, pc_x, ksl0_1602, ksk_1287, \
                         ksk_1288, ksk_1289, ksk_1290, ksl1_1602, lsk_1288, lsk_1289, \
                         lsk_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1602[k] = pa_x[k] * ksl0_1602[k]
                    + f_16 * ksk_1287[k]
                    - f_14 * pc_x[k] * ksl1_1602[k];

        t_1603[k] = f_15 * ksk_1288[k]
                    + f_3 * pc_x[k] * lsk_1288[k];

        t_1604[k] = f_15 * ksk_1289[k]
                    + f_3 * pc_x[k] * lsk_1289[k];

        t_1605[k] = f_15 * ksk_1290[k]
                    + f_3 * pc_x[k] * lsk_1290[k];
    }

#pragma omp simd aligned(t_1606, t_1607, t_1608, t_1609, t_1610, pc_x, pc_y, ksk_1291, \
                         ksk_1292, ksk_1293, ksk_1295, lsk_1287, lsk_1291, lsk_1292, lsk_1293, \
                         lsk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1606[k] = f_15 * ksk_1291[k]
                    + f_3 * pc_x[k] * lsk_1291[k];

        t_1607[k] = f_15 * ksk_1292[k]
                    + f_3 * pc_x[k] * lsk_1292[k];

        t_1608[k] = f_15 * ksk_1293[k]
                    + f_3 * pc_x[k] * lsk_1293[k];

        t_1609[k] = f_3 * pc_y[k] * lsk_1287[k];

        t_1610[k] = f_15 * ksk_1295[k]
                    + f_3 * pc_x[k] * lsk_1295[k];
    }

#pragma omp simd aligned(t_1611, t_1612, t_1613, t_1614, pa_x, pc_x, ksl0_1611, ksl0_1612, \
                         ksl0_1613, ksl0_1614, ksl1_1611, ksl1_1612, ksl1_1613, \
                         ksl1_1614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1611[k] = pa_x[k] * ksl0_1611[k]
                    - f_14 * pc_x[k] * ksl1_1611[k];

        t_1612[k] = pa_x[k] * ksl0_1612[k]
                    - f_14 * pc_x[k] * ksl1_1612[k];

        t_1613[k] = pa_x[k] * ksl0_1613[k]
                    - f_14 * pc_x[k] * ksl1_1613[k];

        t_1614[k] = pa_x[k] * ksl0_1614[k]
                    - f_14 * pc_x[k] * ksl1_1614[k];
    }

#pragma omp simd aligned(t_1615, t_1616, t_1617, t_1618, pa_x, pc_x, pc_y, ksl0_1615, \
                         ksl0_1616, ksl0_1617, ksl1_1615, ksl1_1616, ksl1_1617, \
                         lsk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1615[k] = pa_x[k] * ksl0_1615[k]
                    - f_14 * pc_x[k] * ksl1_1615[k];

        t_1616[k] = pa_x[k] * ksl0_1616[k]
                    - f_14 * pc_x[k] * ksl1_1616[k];

        t_1617[k] = pa_x[k] * ksl0_1617[k]
                    - f_14 * pc_x[k] * ksl1_1617[k];

        t_1618[k] = f_3 * pc_y[k] * lsk_1295[k];
    }

#pragma omp simd aligned(t_1619, t_1620, t_1621, t_1622, pa_x, pc_x, pc_z, ksl0_1619, \
                         ksl1_1619, lsi0_1008, lsi0_1009, lsi1_1008, lsi1_1009, lsk_1296, \
                         lsk_1297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1619[k] = pa_x[k] * ksl0_1619[k]
                    - f_14 * pc_x[k] * ksl1_1619[k];

        t_1620[k] = f_1 * lsi0_1008[k]
                    - f_2 * lsi1_1008[k]
                    + f_3 * pc_x[k] * lsk_1296[k];

        t_1621[k] = f_22 * lsi0_1009[k]
                    - f_23 * lsi1_1009[k]
                    + f_3 * pc_x[k] * lsk_1297[k];

        t_1622[k] = f_3 * pc_z[k] * lsk_1296[k];
    }

#pragma omp simd aligned(t_1623, t_1624, t_1625, t_1626, pc_x, pc_z, lsi0_1011, lsi0_1013, \
                         lsi0_1014, lsi1_1011, lsi1_1013, lsi1_1014, lsk_1297, lsk_1299, \
                         lsk_1301, lsk_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1623[k] = f_12 * lsi0_1011[k]
                    - f_13 * lsi1_1011[k]
                    + f_3 * pc_x[k] * lsk_1299[k];

        t_1624[k] = f_3 * pc_z[k] * lsk_1297[k];

        t_1625[k] = f_12 * lsi0_1013[k]
                    - f_13 * lsi1_1013[k]
                    + f_3 * pc_x[k] * lsk_1301[k];

        t_1626[k] = f_10 * lsi0_1014[k]
                    - f_11 * lsi1_1014[k]
                    + f_3 * pc_x[k] * lsk_1302[k];
    }

#pragma omp simd aligned(t_1627, t_1628, t_1629, t_1630, pc_x, pc_z, lsi0_1016, lsi0_1017, \
                         lsi0_1018, lsi1_1016, lsi1_1017, lsi1_1018, lsk_1299, lsk_1304, \
                         lsk_1305, lsk_1306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1627[k] = f_3 * pc_z[k] * lsk_1299[k];

        t_1628[k] = f_10 * lsi0_1016[k]
                    - f_11 * lsi1_1016[k]
                    + f_3 * pc_x[k] * lsk_1304[k];

        t_1629[k] = f_10 * lsi0_1017[k]
                    - f_11 * lsi1_1017[k]
                    + f_3 * pc_x[k] * lsk_1305[k];

        t_1630[k] = f_8 * lsi0_1018[k]
                    - f_9 * lsi1_1018[k]
                    + f_3 * pc_x[k] * lsk_1306[k];
    }

#pragma omp simd aligned(t_1631, t_1632, t_1633, t_1634, pc_x, pc_z, lsi0_1020, lsi0_1021, \
                         lsi0_1022, lsi1_1020, lsi1_1021, lsi1_1022, lsk_1302, lsk_1308, \
                         lsk_1309, lsk_1310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1631[k] = f_3 * pc_z[k] * lsk_1302[k];

        t_1632[k] = f_8 * lsi0_1020[k]
                    - f_9 * lsi1_1020[k]
                    + f_3 * pc_x[k] * lsk_1308[k];

        t_1633[k] = f_8 * lsi0_1021[k]
                    - f_9 * lsi1_1021[k]
                    + f_3 * pc_x[k] * lsk_1309[k];

        t_1634[k] = f_8 * lsi0_1022[k]
                    - f_9 * lsi1_1022[k]
                    + f_3 * pc_x[k] * lsk_1310[k];
    }

#pragma omp simd aligned(t_1635, t_1636, t_1637, t_1638, pc_x, pc_z, lsi0_1023, lsi0_1025, \
                         lsi0_1026, lsi1_1023, lsi1_1025, lsi1_1026, lsk_1306, lsk_1311, \
                         lsk_1313, lsk_1314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1635[k] = f_6 * lsi0_1023[k]
                    - f_7 * lsi1_1023[k]
                    + f_3 * pc_x[k] * lsk_1311[k];

        t_1636[k] = f_3 * pc_z[k] * lsk_1306[k];

        t_1637[k] = f_6 * lsi0_1025[k]
                    - f_7 * lsi1_1025[k]
                    + f_3 * pc_x[k] * lsk_1313[k];

        t_1638[k] = f_6 * lsi0_1026[k]
                    - f_7 * lsi1_1026[k]
                    + f_3 * pc_x[k] * lsk_1314[k];
    }

#pragma omp simd aligned(t_1639, t_1640, t_1641, t_1642, pc_x, pc_z, lsi0_1027, lsi0_1028, \
                         lsi0_1029, lsi1_1027, lsi1_1028, lsi1_1029, lsk_1311, lsk_1315, \
                         lsk_1316, lsk_1317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1639[k] = f_6 * lsi0_1027[k]
                    - f_7 * lsi1_1027[k]
                    + f_3 * pc_x[k] * lsk_1315[k];

        t_1640[k] = f_6 * lsi0_1028[k]
                    - f_7 * lsi1_1028[k]
                    + f_3 * pc_x[k] * lsk_1316[k];

        t_1641[k] = f_4 * lsi0_1029[k]
                    - f_5 * lsi1_1029[k]
                    + f_3 * pc_x[k] * lsk_1317[k];

        t_1642[k] = f_3 * pc_z[k] * lsk_1311[k];
    }

#pragma omp simd aligned(t_1643, t_1644, t_1645, pc_x, lsi0_1031, lsi0_1032, lsi0_1033, \
                         lsi1_1031, lsi1_1032, lsi1_1033, lsk_1319, lsk_1320, \
                         lsk_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1643[k] = f_4 * lsi0_1031[k]
                    - f_5 * lsi1_1031[k]
                    + f_3 * pc_x[k] * lsk_1319[k];

        t_1644[k] = f_4 * lsi0_1032[k]
                    - f_5 * lsi1_1032[k]
                    + f_3 * pc_x[k] * lsk_1320[k];

        t_1645[k] = f_4 * lsi0_1033[k]
                    - f_5 * lsi1_1033[k]
                    + f_3 * pc_x[k] * lsk_1321[k];
    }

#pragma omp simd aligned(t_1646, t_1647, t_1648, t_1649, t_1650, pc_x, lsi0_1034, lsi0_1035, \
                         lsi1_1034, lsi1_1035, lsk_1322, lsk_1323, lsk_1324, lsk_1325, \
                         lsk_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1646[k] = f_4 * lsi0_1034[k]
                    - f_5 * lsi1_1034[k]
                    + f_3 * pc_x[k] * lsk_1322[k];

        t_1647[k] = f_4 * lsi0_1035[k]
                    - f_5 * lsi1_1035[k]
                    + f_3 * pc_x[k] * lsk_1323[k];

        t_1648[k] = f_3 * pc_x[k] * lsk_1324[k];

        t_1649[k] = f_3 * pc_x[k] * lsk_1325[k];

        t_1650[k] = f_3 * pc_x[k] * lsk_1326[k];
    }

#pragma omp simd aligned(t_1651, t_1652, t_1653, t_1654, t_1655, pc_x, lsk_1327, lsk_1328, \
                         lsk_1329, lsk_1330, lsk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1651[k] = f_3 * pc_x[k] * lsk_1327[k];

        t_1652[k] = f_3 * pc_x[k] * lsk_1328[k];

        t_1653[k] = f_3 * pc_x[k] * lsk_1329[k];

        t_1654[k] = f_3 * pc_x[k] * lsk_1330[k];

        t_1655[k] = f_3 * pc_x[k] * lsk_1331[k];
    }

#pragma omp simd aligned(t_1656, t_1657, t_1658, t_1659, pc_y, pc_z, ksk_1036, lsi0_1029, \
                         lsi0_1030, lsi1_1029, lsi1_1030, lsk_1324, lsk_1325, \
                         lsk_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1656[k] = f_0 * ksk_1036[k]
                    + f_1 * lsi0_1029[k]
                    - f_2 * lsi1_1029[k]
                    + f_3 * pc_y[k] * lsk_1324[k];

        t_1657[k] = f_3 * pc_z[k] * lsk_1324[k];

        t_1658[k] = f_4 * lsi0_1029[k]
                    - f_5 * lsi1_1029[k]
                    + f_3 * pc_z[k] * lsk_1325[k];

        t_1659[k] = f_6 * lsi0_1030[k]
                    - f_7 * lsi1_1030[k]
                    + f_3 * pc_z[k] * lsk_1326[k];
    }

#pragma omp simd aligned(t_1660, t_1661, t_1662, pc_z, lsi0_1031, lsi0_1032, lsi0_1033, \
                         lsi1_1031, lsi1_1032, lsi1_1033, lsk_1327, lsk_1328, \
                         lsk_1329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1660[k] = f_8 * lsi0_1031[k]
                    - f_9 * lsi1_1031[k]
                    + f_3 * pc_z[k] * lsk_1327[k];

        t_1661[k] = f_10 * lsi0_1032[k]
                    - f_11 * lsi1_1032[k]
                    + f_3 * pc_z[k] * lsk_1328[k];

        t_1662[k] = f_12 * lsi0_1033[k]
                    - f_13 * lsi1_1033[k]
                    + f_3 * pc_z[k] * lsk_1329[k];
    }

#pragma omp simd aligned(t_1663, t_1664, t_1665, t_1666, pa_z, pc_y, pc_z, ksl0_1260, \
                         ksl0_1261, ksk_1043, ksl1_1260, ksl1_1261, lsi0_1035, lsi1_1035, \
                         lsk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1663[k] = f_0 * ksk_1043[k]
                    + f_3 * pc_y[k] * lsk_1331[k];

        t_1664[k] = f_1 * lsi0_1035[k]
                    - f_2 * lsi1_1035[k]
                    + f_3 * pc_z[k] * lsk_1331[k];

        t_1665[k] = pa_z[k] * ksl0_1260[k]
                    - f_14 * pc_z[k] * ksl1_1260[k];

        t_1666[k] = pa_z[k] * ksl0_1261[k]
                    - f_14 * pc_z[k] * ksl1_1261[k];
    }

#pragma omp simd aligned(t_1667, t_1668, t_1669, pa_z, pc_x, pc_z, ksl0_1263, ksl1_1263, \
                         lsi0_1038, lsi0_1040, lsi1_1038, lsi1_1040, lsk_1334, \
                         lsk_1336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1667[k] = f_22 * lsi0_1038[k]
                    - f_23 * lsi1_1038[k]
                    + f_3 * pc_x[k] * lsk_1334[k];

        t_1668[k] = pa_z[k] * ksl0_1263[k]
                    - f_14 * pc_z[k] * ksl1_1263[k];

        t_1669[k] = f_12 * lsi0_1040[k]
                    - f_13 * lsi1_1040[k]
                    + f_3 * pc_x[k] * lsk_1336[k];
    }

#pragma omp simd aligned(t_1670, t_1671, t_1672, pa_z, pc_x, pc_z, ksl0_1266, ksl1_1266, \
                         lsi0_1041, lsi0_1043, lsi1_1041, lsi1_1043, lsk_1337, \
                         lsk_1339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1670[k] = f_12 * lsi0_1041[k]
                    - f_13 * lsi1_1041[k]
                    + f_3 * pc_x[k] * lsk_1337[k];

        t_1671[k] = pa_z[k] * ksl0_1266[k]
                    - f_14 * pc_z[k] * ksl1_1266[k];

        t_1672[k] = f_10 * lsi0_1043[k]
                    - f_11 * lsi1_1043[k]
                    + f_3 * pc_x[k] * lsk_1339[k];
    }

#pragma omp simd aligned(t_1673, t_1674, t_1675, pa_z, pc_x, pc_z, ksl0_1270, ksl1_1270, \
                         lsi0_1044, lsi0_1045, lsi1_1044, lsi1_1045, lsk_1340, \
                         lsk_1341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1673[k] = f_10 * lsi0_1044[k]
                    - f_11 * lsi1_1044[k]
                    + f_3 * pc_x[k] * lsk_1340[k];

        t_1674[k] = f_10 * lsi0_1045[k]
                    - f_11 * lsi1_1045[k]
                    + f_3 * pc_x[k] * lsk_1341[k];

        t_1675[k] = pa_z[k] * ksl0_1270[k]
                    - f_14 * pc_z[k] * ksl1_1270[k];
    }

#pragma omp simd aligned(t_1676, t_1677, t_1678, pc_x, lsi0_1047, lsi0_1048, lsi0_1049, \
                         lsi1_1047, lsi1_1048, lsi1_1049, lsk_1343, lsk_1344, \
                         lsk_1345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1676[k] = f_8 * lsi0_1047[k]
                    - f_9 * lsi1_1047[k]
                    + f_3 * pc_x[k] * lsk_1343[k];

        t_1677[k] = f_8 * lsi0_1048[k]
                    - f_9 * lsi1_1048[k]
                    + f_3 * pc_x[k] * lsk_1344[k];

        t_1678[k] = f_8 * lsi0_1049[k]
                    - f_9 * lsi1_1049[k]
                    + f_3 * pc_x[k] * lsk_1345[k];
    }

#pragma omp simd aligned(t_1679, t_1680, t_1681, pa_z, pc_x, pc_z, ksl0_1275, ksl1_1275, \
                         lsi0_1050, lsi0_1052, lsi1_1050, lsi1_1052, lsk_1346, \
                         lsk_1348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1679[k] = f_8 * lsi0_1050[k]
                    - f_9 * lsi1_1050[k]
                    + f_3 * pc_x[k] * lsk_1346[k];

        t_1680[k] = pa_z[k] * ksl0_1275[k]
                    - f_14 * pc_z[k] * ksl1_1275[k];

        t_1681[k] = f_6 * lsi0_1052[k]
                    - f_7 * lsi1_1052[k]
                    + f_3 * pc_x[k] * lsk_1348[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, pc_x, lsi0_1053, lsi0_1054, lsi0_1055, \
                         lsi1_1053, lsi1_1054, lsi1_1055, lsk_1349, lsk_1350, \
                         lsk_1351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = f_6 * lsi0_1053[k]
                    - f_7 * lsi1_1053[k]
                    + f_3 * pc_x[k] * lsk_1349[k];

        t_1683[k] = f_6 * lsi0_1054[k]
                    - f_7 * lsi1_1054[k]
                    + f_3 * pc_x[k] * lsk_1350[k];

        t_1684[k] = f_6 * lsi0_1055[k]
                    - f_7 * lsi1_1055[k]
                    + f_3 * pc_x[k] * lsk_1351[k];
    }

#pragma omp simd aligned(t_1685, t_1686, t_1687, pa_z, pc_x, pc_z, ksl0_1281, ksl1_1281, \
                         lsi0_1056, lsi0_1058, lsi1_1056, lsi1_1058, lsk_1352, \
                         lsk_1354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1685[k] = f_6 * lsi0_1056[k]
                    - f_7 * lsi1_1056[k]
                    + f_3 * pc_x[k] * lsk_1352[k];

        t_1686[k] = pa_z[k] * ksl0_1281[k]
                    - f_14 * pc_z[k] * ksl1_1281[k];

        t_1687[k] = f_4 * lsi0_1058[k]
                    - f_5 * lsi1_1058[k]
                    + f_3 * pc_x[k] * lsk_1354[k];
    }

#pragma omp simd aligned(t_1688, t_1689, t_1690, pc_x, lsi0_1059, lsi0_1060, lsi0_1061, \
                         lsi1_1059, lsi1_1060, lsi1_1061, lsk_1355, lsk_1356, \
                         lsk_1357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1688[k] = f_4 * lsi0_1059[k]
                    - f_5 * lsi1_1059[k]
                    + f_3 * pc_x[k] * lsk_1355[k];

        t_1689[k] = f_4 * lsi0_1060[k]
                    - f_5 * lsi1_1060[k]
                    + f_3 * pc_x[k] * lsk_1356[k];

        t_1690[k] = f_4 * lsi0_1061[k]
                    - f_5 * lsi1_1061[k]
                    + f_3 * pc_x[k] * lsk_1357[k];
    }

#pragma omp simd aligned(t_1691, t_1692, t_1693, t_1694, t_1695, pc_x, lsi0_1062, lsi0_1063, \
                         lsi1_1062, lsi1_1063, lsk_1358, lsk_1359, lsk_1360, lsk_1361, \
                         lsk_1362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1691[k] = f_4 * lsi0_1062[k]
                    - f_5 * lsi1_1062[k]
                    + f_3 * pc_x[k] * lsk_1358[k];

        t_1692[k] = f_4 * lsi0_1063[k]
                    - f_5 * lsi1_1063[k]
                    + f_3 * pc_x[k] * lsk_1359[k];

        t_1693[k] = f_3 * pc_x[k] * lsk_1360[k];

        t_1694[k] = f_3 * pc_x[k] * lsk_1361[k];

        t_1695[k] = f_3 * pc_x[k] * lsk_1362[k];
    }

#pragma omp simd aligned(t_1696, t_1697, t_1698, t_1699, t_1700, t_1701, pa_z, pc_x, pc_z, \
                         ksl0_1296, ksl1_1296, lsk_1363, lsk_1364, lsk_1365, lsk_1366, \
                         lsk_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1696[k] = f_3 * pc_x[k] * lsk_1363[k];

        t_1697[k] = f_3 * pc_x[k] * lsk_1364[k];

        t_1698[k] = f_3 * pc_x[k] * lsk_1365[k];

        t_1699[k] = f_3 * pc_x[k] * lsk_1366[k];

        t_1700[k] = f_3 * pc_x[k] * lsk_1367[k];

        t_1701[k] = pa_z[k] * ksl0_1296[k]
                    - f_14 * pc_z[k] * ksl1_1296[k];
    }

#pragma omp simd aligned(t_1702, t_1703, t_1704, pa_z, pc_z, ksl0_1298, ksl0_1299, ksk_1036, \
                         ksk_1037, ksk_1038, ksl1_1298, ksl1_1299, \
                         lsk_1360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1702[k] = f_15 * ksk_1036[k]
                    + f_3 * pc_z[k] * lsk_1360[k];

        t_1703[k] = pa_z[k] * ksl0_1298[k]
                    + f_16 * ksk_1037[k]
                    - f_14 * pc_z[k] * ksl1_1298[k];

        t_1704[k] = pa_z[k] * ksl0_1299[k]
                    + f_17 * ksk_1038[k]
                    - f_14 * pc_z[k] * ksl1_1299[k];
    }

#pragma omp simd aligned(t_1705, t_1706, t_1707, pa_z, pc_z, ksl0_1300, ksl0_1301, ksl0_1302, \
                         ksk_1039, ksk_1040, ksk_1041, ksl1_1300, ksl1_1301, \
                         ksl1_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1705[k] = pa_z[k] * ksl0_1300[k]
                    + f_18 * ksk_1039[k]
                    - f_14 * pc_z[k] * ksl1_1300[k];

        t_1706[k] = pa_z[k] * ksl0_1301[k]
                    + f_19 * ksk_1040[k]
                    - f_14 * pc_z[k] * ksl1_1301[k];

        t_1707[k] = pa_z[k] * ksl0_1302[k]
                    + f_20 * ksk_1041[k]
                    - f_14 * pc_z[k] * ksl1_1302[k];
    }

#pragma omp simd aligned(t_1708, t_1709, t_1710, pc_x, pc_y, pc_z, ksk_1043, ksk_1079, \
                         lsi0_1063, lsi0_1064, lsi1_1063, lsi1_1064, lsk_1367, \
                         lsk_1368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1708[k] = f_21 * ksk_1079[k]
                    + f_3 * pc_y[k] * lsk_1367[k];

        t_1709[k] = f_15 * ksk_1043[k]
                    + f_1 * lsi0_1063[k]
                    - f_2 * lsi1_1063[k]
                    + f_3 * pc_z[k] * lsk_1367[k];

        t_1710[k] = f_1 * lsi0_1064[k]
                    - f_2 * lsi1_1064[k]
                    + f_3 * pc_x[k] * lsk_1368[k];
    }

#pragma omp simd aligned(t_1711, t_1712, t_1713, pc_x, lsi0_1065, lsi0_1066, lsi0_1067, \
                         lsi1_1065, lsi1_1066, lsi1_1067, lsk_1369, lsk_1370, \
                         lsk_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1711[k] = f_22 * lsi0_1065[k]
                    - f_23 * lsi1_1065[k]
                    + f_3 * pc_x[k] * lsk_1369[k];

        t_1712[k] = f_22 * lsi0_1066[k]
                    - f_23 * lsi1_1066[k]
                    + f_3 * pc_x[k] * lsk_1370[k];

        t_1713[k] = f_12 * lsi0_1067[k]
                    - f_13 * lsi1_1067[k]
                    + f_3 * pc_x[k] * lsk_1371[k];
    }

#pragma omp simd aligned(t_1714, t_1715, t_1716, pc_x, lsi0_1068, lsi0_1069, lsi0_1070, \
                         lsi1_1068, lsi1_1069, lsi1_1070, lsk_1372, lsk_1373, \
                         lsk_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1714[k] = f_12 * lsi0_1068[k]
                    - f_13 * lsi1_1068[k]
                    + f_3 * pc_x[k] * lsk_1372[k];

        t_1715[k] = f_12 * lsi0_1069[k]
                    - f_13 * lsi1_1069[k]
                    + f_3 * pc_x[k] * lsk_1373[k];

        t_1716[k] = f_10 * lsi0_1070[k]
                    - f_11 * lsi1_1070[k]
                    + f_3 * pc_x[k] * lsk_1374[k];
    }

#pragma omp simd aligned(t_1717, t_1718, t_1719, pc_x, lsi0_1071, lsi0_1072, lsi0_1073, \
                         lsi1_1071, lsi1_1072, lsi1_1073, lsk_1375, lsk_1376, \
                         lsk_1377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1717[k] = f_10 * lsi0_1071[k]
                    - f_11 * lsi1_1071[k]
                    + f_3 * pc_x[k] * lsk_1375[k];

        t_1718[k] = f_10 * lsi0_1072[k]
                    - f_11 * lsi1_1072[k]
                    + f_3 * pc_x[k] * lsk_1376[k];

        t_1719[k] = f_10 * lsi0_1073[k]
                    - f_11 * lsi1_1073[k]
                    + f_3 * pc_x[k] * lsk_1377[k];
    }

#pragma omp simd aligned(t_1720, t_1721, t_1722, pc_x, lsi0_1074, lsi0_1075, lsi0_1076, \
                         lsi1_1074, lsi1_1075, lsi1_1076, lsk_1378, lsk_1379, \
                         lsk_1380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1720[k] = f_8 * lsi0_1074[k]
                    - f_9 * lsi1_1074[k]
                    + f_3 * pc_x[k] * lsk_1378[k];

        t_1721[k] = f_8 * lsi0_1075[k]
                    - f_9 * lsi1_1075[k]
                    + f_3 * pc_x[k] * lsk_1379[k];

        t_1722[k] = f_8 * lsi0_1076[k]
                    - f_9 * lsi1_1076[k]
                    + f_3 * pc_x[k] * lsk_1380[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t ksk, const size_t lsi0,
                                                           const size_t lsi1, const size_t lsk,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksk_1072 = buffer.data(ksk + 1072);
    const auto *ksk_1079 = buffer.data(ksk + 1079);
    const auto *ksk_1108 = buffer.data(ksk + 1108);
    const auto *ksk_1110 = buffer.data(ksk + 1110);
    const auto *ksk_1111 = buffer.data(ksk + 1111);
    const auto *ksk_1112 = buffer.data(ksk + 1112);
    const auto *ksk_1113 = buffer.data(ksk + 1113);
    const auto *ksk_1114 = buffer.data(ksk + 1114);
    const auto *ksk_1115 = buffer.data(ksk + 1115);
    const auto *ksk_1144 = buffer.data(ksk + 1144);
    const auto *ksk_1146 = buffer.data(ksk + 1146);
    const auto *ksk_1147 = buffer.data(ksk + 1147);
    const auto *ksk_1148 = buffer.data(ksk + 1148);
    const auto *ksk_1149 = buffer.data(ksk + 1149);
    const auto *ksk_1150 = buffer.data(ksk + 1150);
    const auto *ksk_1151 = buffer.data(ksk + 1151);
    const auto *ksk_1180 = buffer.data(ksk + 1180);
    const auto *ksk_1182 = buffer.data(ksk + 1182);

    const auto *lsi0_1077 = buffer.data(lsi0 + 1077);
    const auto *lsi0_1078 = buffer.data(lsi0 + 1078);
    const auto *lsi0_1079 = buffer.data(lsi0 + 1079);
    const auto *lsi0_1080 = buffer.data(lsi0 + 1080);
    const auto *lsi0_1081 = buffer.data(lsi0 + 1081);
    const auto *lsi0_1082 = buffer.data(lsi0 + 1082);
    const auto *lsi0_1083 = buffer.data(lsi0 + 1083);
    const auto *lsi0_1084 = buffer.data(lsi0 + 1084);
    const auto *lsi0_1085 = buffer.data(lsi0 + 1085);
    const auto *lsi0_1086 = buffer.data(lsi0 + 1086);
    const auto *lsi0_1087 = buffer.data(lsi0 + 1087);
    const auto *lsi0_1088 = buffer.data(lsi0 + 1088);
    const auto *lsi0_1089 = buffer.data(lsi0 + 1089);
    const auto *lsi0_1090 = buffer.data(lsi0 + 1090);
    const auto *lsi0_1091 = buffer.data(lsi0 + 1091);
    const auto *lsi0_1092 = buffer.data(lsi0 + 1092);
    const auto *lsi0_1093 = buffer.data(lsi0 + 1093);
    const auto *lsi0_1094 = buffer.data(lsi0 + 1094);
    const auto *lsi0_1095 = buffer.data(lsi0 + 1095);
    const auto *lsi0_1096 = buffer.data(lsi0 + 1096);
    const auto *lsi0_1097 = buffer.data(lsi0 + 1097);
    const auto *lsi0_1098 = buffer.data(lsi0 + 1098);
    const auto *lsi0_1099 = buffer.data(lsi0 + 1099);
    const auto *lsi0_1100 = buffer.data(lsi0 + 1100);
    const auto *lsi0_1101 = buffer.data(lsi0 + 1101);
    const auto *lsi0_1102 = buffer.data(lsi0 + 1102);
    const auto *lsi0_1103 = buffer.data(lsi0 + 1103);
    const auto *lsi0_1104 = buffer.data(lsi0 + 1104);
    const auto *lsi0_1105 = buffer.data(lsi0 + 1105);
    const auto *lsi0_1106 = buffer.data(lsi0 + 1106);
    const auto *lsi0_1107 = buffer.data(lsi0 + 1107);
    const auto *lsi0_1108 = buffer.data(lsi0 + 1108);
    const auto *lsi0_1109 = buffer.data(lsi0 + 1109);
    const auto *lsi0_1110 = buffer.data(lsi0 + 1110);
    const auto *lsi0_1111 = buffer.data(lsi0 + 1111);
    const auto *lsi0_1112 = buffer.data(lsi0 + 1112);
    const auto *lsi0_1113 = buffer.data(lsi0 + 1113);
    const auto *lsi0_1114 = buffer.data(lsi0 + 1114);
    const auto *lsi0_1115 = buffer.data(lsi0 + 1115);
    const auto *lsi0_1116 = buffer.data(lsi0 + 1116);
    const auto *lsi0_1117 = buffer.data(lsi0 + 1117);
    const auto *lsi0_1118 = buffer.data(lsi0 + 1118);
    const auto *lsi0_1119 = buffer.data(lsi0 + 1119);
    const auto *lsi0_1120 = buffer.data(lsi0 + 1120);
    const auto *lsi0_1121 = buffer.data(lsi0 + 1121);
    const auto *lsi0_1122 = buffer.data(lsi0 + 1122);
    const auto *lsi0_1123 = buffer.data(lsi0 + 1123);
    const auto *lsi0_1124 = buffer.data(lsi0 + 1124);
    const auto *lsi0_1125 = buffer.data(lsi0 + 1125);
    const auto *lsi0_1126 = buffer.data(lsi0 + 1126);
    const auto *lsi0_1127 = buffer.data(lsi0 + 1127);
    const auto *lsi0_1128 = buffer.data(lsi0 + 1128);
    const auto *lsi0_1129 = buffer.data(lsi0 + 1129);
    const auto *lsi0_1130 = buffer.data(lsi0 + 1130);
    const auto *lsi0_1131 = buffer.data(lsi0 + 1131);
    const auto *lsi0_1132 = buffer.data(lsi0 + 1132);
    const auto *lsi0_1133 = buffer.data(lsi0 + 1133);
    const auto *lsi0_1134 = buffer.data(lsi0 + 1134);
    const auto *lsi0_1135 = buffer.data(lsi0 + 1135);
    const auto *lsi0_1136 = buffer.data(lsi0 + 1136);
    const auto *lsi0_1137 = buffer.data(lsi0 + 1137);
    const auto *lsi0_1138 = buffer.data(lsi0 + 1138);
    const auto *lsi0_1139 = buffer.data(lsi0 + 1139);
    const auto *lsi0_1140 = buffer.data(lsi0 + 1140);
    const auto *lsi0_1141 = buffer.data(lsi0 + 1141);
    const auto *lsi0_1142 = buffer.data(lsi0 + 1142);
    const auto *lsi0_1143 = buffer.data(lsi0 + 1143);
    const auto *lsi0_1144 = buffer.data(lsi0 + 1144);
    const auto *lsi0_1145 = buffer.data(lsi0 + 1145);
    const auto *lsi0_1146 = buffer.data(lsi0 + 1146);
    const auto *lsi0_1147 = buffer.data(lsi0 + 1147);

    const auto *lsi1_1077 = buffer.data(lsi1 + 1077);
    const auto *lsi1_1078 = buffer.data(lsi1 + 1078);
    const auto *lsi1_1079 = buffer.data(lsi1 + 1079);
    const auto *lsi1_1080 = buffer.data(lsi1 + 1080);
    const auto *lsi1_1081 = buffer.data(lsi1 + 1081);
    const auto *lsi1_1082 = buffer.data(lsi1 + 1082);
    const auto *lsi1_1083 = buffer.data(lsi1 + 1083);
    const auto *lsi1_1084 = buffer.data(lsi1 + 1084);
    const auto *lsi1_1085 = buffer.data(lsi1 + 1085);
    const auto *lsi1_1086 = buffer.data(lsi1 + 1086);
    const auto *lsi1_1087 = buffer.data(lsi1 + 1087);
    const auto *lsi1_1088 = buffer.data(lsi1 + 1088);
    const auto *lsi1_1089 = buffer.data(lsi1 + 1089);
    const auto *lsi1_1090 = buffer.data(lsi1 + 1090);
    const auto *lsi1_1091 = buffer.data(lsi1 + 1091);
    const auto *lsi1_1092 = buffer.data(lsi1 + 1092);
    const auto *lsi1_1093 = buffer.data(lsi1 + 1093);
    const auto *lsi1_1094 = buffer.data(lsi1 + 1094);
    const auto *lsi1_1095 = buffer.data(lsi1 + 1095);
    const auto *lsi1_1096 = buffer.data(lsi1 + 1096);
    const auto *lsi1_1097 = buffer.data(lsi1 + 1097);
    const auto *lsi1_1098 = buffer.data(lsi1 + 1098);
    const auto *lsi1_1099 = buffer.data(lsi1 + 1099);
    const auto *lsi1_1100 = buffer.data(lsi1 + 1100);
    const auto *lsi1_1101 = buffer.data(lsi1 + 1101);
    const auto *lsi1_1102 = buffer.data(lsi1 + 1102);
    const auto *lsi1_1103 = buffer.data(lsi1 + 1103);
    const auto *lsi1_1104 = buffer.data(lsi1 + 1104);
    const auto *lsi1_1105 = buffer.data(lsi1 + 1105);
    const auto *lsi1_1106 = buffer.data(lsi1 + 1106);
    const auto *lsi1_1107 = buffer.data(lsi1 + 1107);
    const auto *lsi1_1108 = buffer.data(lsi1 + 1108);
    const auto *lsi1_1109 = buffer.data(lsi1 + 1109);
    const auto *lsi1_1110 = buffer.data(lsi1 + 1110);
    const auto *lsi1_1111 = buffer.data(lsi1 + 1111);
    const auto *lsi1_1112 = buffer.data(lsi1 + 1112);
    const auto *lsi1_1113 = buffer.data(lsi1 + 1113);
    const auto *lsi1_1114 = buffer.data(lsi1 + 1114);
    const auto *lsi1_1115 = buffer.data(lsi1 + 1115);
    const auto *lsi1_1116 = buffer.data(lsi1 + 1116);
    const auto *lsi1_1117 = buffer.data(lsi1 + 1117);
    const auto *lsi1_1118 = buffer.data(lsi1 + 1118);
    const auto *lsi1_1119 = buffer.data(lsi1 + 1119);
    const auto *lsi1_1120 = buffer.data(lsi1 + 1120);
    const auto *lsi1_1121 = buffer.data(lsi1 + 1121);
    const auto *lsi1_1122 = buffer.data(lsi1 + 1122);
    const auto *lsi1_1123 = buffer.data(lsi1 + 1123);
    const auto *lsi1_1124 = buffer.data(lsi1 + 1124);
    const auto *lsi1_1125 = buffer.data(lsi1 + 1125);
    const auto *lsi1_1126 = buffer.data(lsi1 + 1126);
    const auto *lsi1_1127 = buffer.data(lsi1 + 1127);
    const auto *lsi1_1128 = buffer.data(lsi1 + 1128);
    const auto *lsi1_1129 = buffer.data(lsi1 + 1129);
    const auto *lsi1_1130 = buffer.data(lsi1 + 1130);
    const auto *lsi1_1131 = buffer.data(lsi1 + 1131);
    const auto *lsi1_1132 = buffer.data(lsi1 + 1132);
    const auto *lsi1_1133 = buffer.data(lsi1 + 1133);
    const auto *lsi1_1134 = buffer.data(lsi1 + 1134);
    const auto *lsi1_1135 = buffer.data(lsi1 + 1135);
    const auto *lsi1_1136 = buffer.data(lsi1 + 1136);
    const auto *lsi1_1137 = buffer.data(lsi1 + 1137);
    const auto *lsi1_1138 = buffer.data(lsi1 + 1138);
    const auto *lsi1_1139 = buffer.data(lsi1 + 1139);
    const auto *lsi1_1140 = buffer.data(lsi1 + 1140);
    const auto *lsi1_1141 = buffer.data(lsi1 + 1141);
    const auto *lsi1_1142 = buffer.data(lsi1 + 1142);
    const auto *lsi1_1143 = buffer.data(lsi1 + 1143);
    const auto *lsi1_1144 = buffer.data(lsi1 + 1144);
    const auto *lsi1_1145 = buffer.data(lsi1 + 1145);
    const auto *lsi1_1146 = buffer.data(lsi1 + 1146);
    const auto *lsi1_1147 = buffer.data(lsi1 + 1147);

    const auto *lsk_1381 = buffer.data(lsk + 1381);
    const auto *lsk_1382 = buffer.data(lsk + 1382);
    const auto *lsk_1383 = buffer.data(lsk + 1383);
    const auto *lsk_1384 = buffer.data(lsk + 1384);
    const auto *lsk_1385 = buffer.data(lsk + 1385);
    const auto *lsk_1386 = buffer.data(lsk + 1386);
    const auto *lsk_1387 = buffer.data(lsk + 1387);
    const auto *lsk_1388 = buffer.data(lsk + 1388);
    const auto *lsk_1389 = buffer.data(lsk + 1389);
    const auto *lsk_1390 = buffer.data(lsk + 1390);
    const auto *lsk_1391 = buffer.data(lsk + 1391);
    const auto *lsk_1392 = buffer.data(lsk + 1392);
    const auto *lsk_1393 = buffer.data(lsk + 1393);
    const auto *lsk_1394 = buffer.data(lsk + 1394);
    const auto *lsk_1395 = buffer.data(lsk + 1395);
    const auto *lsk_1396 = buffer.data(lsk + 1396);
    const auto *lsk_1397 = buffer.data(lsk + 1397);
    const auto *lsk_1398 = buffer.data(lsk + 1398);
    const auto *lsk_1399 = buffer.data(lsk + 1399);
    const auto *lsk_1400 = buffer.data(lsk + 1400);
    const auto *lsk_1401 = buffer.data(lsk + 1401);
    const auto *lsk_1402 = buffer.data(lsk + 1402);
    const auto *lsk_1403 = buffer.data(lsk + 1403);
    const auto *lsk_1404 = buffer.data(lsk + 1404);
    const auto *lsk_1405 = buffer.data(lsk + 1405);
    const auto *lsk_1406 = buffer.data(lsk + 1406);
    const auto *lsk_1407 = buffer.data(lsk + 1407);
    const auto *lsk_1408 = buffer.data(lsk + 1408);
    const auto *lsk_1409 = buffer.data(lsk + 1409);
    const auto *lsk_1410 = buffer.data(lsk + 1410);
    const auto *lsk_1411 = buffer.data(lsk + 1411);
    const auto *lsk_1412 = buffer.data(lsk + 1412);
    const auto *lsk_1413 = buffer.data(lsk + 1413);
    const auto *lsk_1414 = buffer.data(lsk + 1414);
    const auto *lsk_1415 = buffer.data(lsk + 1415);
    const auto *lsk_1416 = buffer.data(lsk + 1416);
    const auto *lsk_1417 = buffer.data(lsk + 1417);
    const auto *lsk_1418 = buffer.data(lsk + 1418);
    const auto *lsk_1419 = buffer.data(lsk + 1419);
    const auto *lsk_1420 = buffer.data(lsk + 1420);
    const auto *lsk_1421 = buffer.data(lsk + 1421);
    const auto *lsk_1422 = buffer.data(lsk + 1422);
    const auto *lsk_1423 = buffer.data(lsk + 1423);
    const auto *lsk_1424 = buffer.data(lsk + 1424);
    const auto *lsk_1425 = buffer.data(lsk + 1425);
    const auto *lsk_1426 = buffer.data(lsk + 1426);
    const auto *lsk_1427 = buffer.data(lsk + 1427);
    const auto *lsk_1428 = buffer.data(lsk + 1428);
    const auto *lsk_1429 = buffer.data(lsk + 1429);
    const auto *lsk_1430 = buffer.data(lsk + 1430);
    const auto *lsk_1431 = buffer.data(lsk + 1431);
    const auto *lsk_1432 = buffer.data(lsk + 1432);
    const auto *lsk_1433 = buffer.data(lsk + 1433);
    const auto *lsk_1434 = buffer.data(lsk + 1434);
    const auto *lsk_1435 = buffer.data(lsk + 1435);
    const auto *lsk_1436 = buffer.data(lsk + 1436);
    const auto *lsk_1437 = buffer.data(lsk + 1437);
    const auto *lsk_1438 = buffer.data(lsk + 1438);
    const auto *lsk_1439 = buffer.data(lsk + 1439);
    const auto *lsk_1440 = buffer.data(lsk + 1440);
    const auto *lsk_1441 = buffer.data(lsk + 1441);
    const auto *lsk_1442 = buffer.data(lsk + 1442);
    const auto *lsk_1443 = buffer.data(lsk + 1443);
    const auto *lsk_1444 = buffer.data(lsk + 1444);
    const auto *lsk_1445 = buffer.data(lsk + 1445);
    const auto *lsk_1446 = buffer.data(lsk + 1446);
    const auto *lsk_1447 = buffer.data(lsk + 1447);
    const auto *lsk_1448 = buffer.data(lsk + 1448);
    const auto *lsk_1449 = buffer.data(lsk + 1449);
    const auto *lsk_1450 = buffer.data(lsk + 1450);
    const auto *lsk_1451 = buffer.data(lsk + 1451);
    const auto *lsk_1452 = buffer.data(lsk + 1452);
    const auto *lsk_1453 = buffer.data(lsk + 1453);
    const auto *lsk_1454 = buffer.data(lsk + 1454);
    const auto *lsk_1455 = buffer.data(lsk + 1455);
    const auto *lsk_1456 = buffer.data(lsk + 1456);
    const auto *lsk_1457 = buffer.data(lsk + 1457);
    const auto *lsk_1458 = buffer.data(lsk + 1458);
    const auto *lsk_1459 = buffer.data(lsk + 1459);
    const auto *lsk_1460 = buffer.data(lsk + 1460);
    const auto *lsk_1461 = buffer.data(lsk + 1461);
    const auto *lsk_1462 = buffer.data(lsk + 1462);
    const auto *lsk_1463 = buffer.data(lsk + 1463);
    const auto *lsk_1464 = buffer.data(lsk + 1464);
    const auto *lsk_1465 = buffer.data(lsk + 1465);
    const auto *lsk_1466 = buffer.data(lsk + 1466);
    const auto *lsk_1467 = buffer.data(lsk + 1467);
    const auto *lsk_1468 = buffer.data(lsk + 1468);
    const auto *lsk_1469 = buffer.data(lsk + 1469);
    const auto *lsk_1470 = buffer.data(lsk + 1470);
    const auto *lsk_1471 = buffer.data(lsk + 1471);
    const auto *lsk_1472 = buffer.data(lsk + 1472);
    const auto *lsk_1473 = buffer.data(lsk + 1473);
    const auto *lsk_1474 = buffer.data(lsk + 1474);
    const auto *lsk_1475 = buffer.data(lsk + 1475);

#pragma omp simd aligned(t_1723, t_1724, t_1725, pc_x, lsi0_1077, lsi0_1078, lsi0_1079, \
                         lsi1_1077, lsi1_1078, lsi1_1079, lsk_1381, lsk_1382, \
                         lsk_1383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1723[k] = f_8 * lsi0_1077[k]
                    - f_9 * lsi1_1077[k]
                    + f_3 * pc_x[k] * lsk_1381[k];

        t_1724[k] = f_8 * lsi0_1078[k]
                    - f_9 * lsi1_1078[k]
                    + f_3 * pc_x[k] * lsk_1382[k];

        t_1725[k] = f_6 * lsi0_1079[k]
                    - f_7 * lsi1_1079[k]
                    + f_3 * pc_x[k] * lsk_1383[k];
    }

#pragma omp simd aligned(t_1726, t_1727, t_1728, pc_x, lsi0_1080, lsi0_1081, lsi0_1082, \
                         lsi1_1080, lsi1_1081, lsi1_1082, lsk_1384, lsk_1385, \
                         lsk_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1726[k] = f_6 * lsi0_1080[k]
                    - f_7 * lsi1_1080[k]
                    + f_3 * pc_x[k] * lsk_1384[k];

        t_1727[k] = f_6 * lsi0_1081[k]
                    - f_7 * lsi1_1081[k]
                    + f_3 * pc_x[k] * lsk_1385[k];

        t_1728[k] = f_6 * lsi0_1082[k]
                    - f_7 * lsi1_1082[k]
                    + f_3 * pc_x[k] * lsk_1386[k];
    }

#pragma omp simd aligned(t_1729, t_1730, t_1731, pc_x, lsi0_1083, lsi0_1084, lsi0_1085, \
                         lsi1_1083, lsi1_1084, lsi1_1085, lsk_1387, lsk_1388, \
                         lsk_1389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1729[k] = f_6 * lsi0_1083[k]
                    - f_7 * lsi1_1083[k]
                    + f_3 * pc_x[k] * lsk_1387[k];

        t_1730[k] = f_6 * lsi0_1084[k]
                    - f_7 * lsi1_1084[k]
                    + f_3 * pc_x[k] * lsk_1388[k];

        t_1731[k] = f_4 * lsi0_1085[k]
                    - f_5 * lsi1_1085[k]
                    + f_3 * pc_x[k] * lsk_1389[k];
    }

#pragma omp simd aligned(t_1732, t_1733, t_1734, pc_x, lsi0_1086, lsi0_1087, lsi0_1088, \
                         lsi1_1086, lsi1_1087, lsi1_1088, lsk_1390, lsk_1391, \
                         lsk_1392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1732[k] = f_4 * lsi0_1086[k]
                    - f_5 * lsi1_1086[k]
                    + f_3 * pc_x[k] * lsk_1390[k];

        t_1733[k] = f_4 * lsi0_1087[k]
                    - f_5 * lsi1_1087[k]
                    + f_3 * pc_x[k] * lsk_1391[k];

        t_1734[k] = f_4 * lsi0_1088[k]
                    - f_5 * lsi1_1088[k]
                    + f_3 * pc_x[k] * lsk_1392[k];
    }

#pragma omp simd aligned(t_1735, t_1736, t_1737, t_1738, pc_x, lsi0_1089, lsi0_1090, \
                         lsi0_1091, lsi1_1089, lsi1_1090, lsi1_1091, lsk_1393, lsk_1394, \
                         lsk_1395, lsk_1396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1735[k] = f_4 * lsi0_1089[k]
                    - f_5 * lsi1_1089[k]
                    + f_3 * pc_x[k] * lsk_1393[k];

        t_1736[k] = f_4 * lsi0_1090[k]
                    - f_5 * lsi1_1090[k]
                    + f_3 * pc_x[k] * lsk_1394[k];

        t_1737[k] = f_4 * lsi0_1091[k]
                    - f_5 * lsi1_1091[k]
                    + f_3 * pc_x[k] * lsk_1395[k];

        t_1738[k] = f_3 * pc_x[k] * lsk_1396[k];
    }

#pragma omp simd aligned(t_1739, t_1740, t_1741, t_1742, t_1743, t_1744, t_1745, pc_x, \
                         lsk_1397, lsk_1398, lsk_1399, lsk_1400, lsk_1401, lsk_1402, \
                         lsk_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1739[k] = f_3 * pc_x[k] * lsk_1397[k];

        t_1740[k] = f_3 * pc_x[k] * lsk_1398[k];

        t_1741[k] = f_3 * pc_x[k] * lsk_1399[k];

        t_1742[k] = f_3 * pc_x[k] * lsk_1400[k];

        t_1743[k] = f_3 * pc_x[k] * lsk_1401[k];

        t_1744[k] = f_3 * pc_x[k] * lsk_1402[k];

        t_1745[k] = f_3 * pc_x[k] * lsk_1403[k];
    }

#pragma omp simd aligned(t_1746, t_1747, t_1748, pc_y, pc_z, ksk_1072, ksk_1108, ksk_1110, \
                         lsi0_1085, lsi0_1087, lsi1_1085, lsi1_1087, lsk_1396, \
                         lsk_1398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1746[k] = f_20 * ksk_1108[k]
                    + f_1 * lsi0_1085[k]
                    - f_2 * lsi1_1085[k]
                    + f_3 * pc_y[k] * lsk_1396[k];

        t_1747[k] = f_16 * ksk_1072[k]
                    + f_3 * pc_z[k] * lsk_1396[k];

        t_1748[k] = f_20 * ksk_1110[k]
                    + f_12 * lsi0_1087[k]
                    - f_13 * lsi1_1087[k]
                    + f_3 * pc_y[k] * lsk_1398[k];
    }

#pragma omp simd aligned(t_1749, t_1750, t_1751, pc_y, ksk_1111, ksk_1112, ksk_1113, \
                         lsi0_1088, lsi0_1089, lsi0_1090, lsi1_1088, lsi1_1089, lsi1_1090, \
                         lsk_1399, lsk_1400, lsk_1401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1749[k] = f_20 * ksk_1111[k]
                    + f_10 * lsi0_1088[k]
                    - f_11 * lsi1_1088[k]
                    + f_3 * pc_y[k] * lsk_1399[k];

        t_1750[k] = f_20 * ksk_1112[k]
                    + f_8 * lsi0_1089[k]
                    - f_9 * lsi1_1089[k]
                    + f_3 * pc_y[k] * lsk_1400[k];

        t_1751[k] = f_20 * ksk_1113[k]
                    + f_6 * lsi0_1090[k]
                    - f_7 * lsi1_1090[k]
                    + f_3 * pc_y[k] * lsk_1401[k];
    }

#pragma omp simd aligned(t_1752, t_1753, t_1754, pc_y, pc_z, ksk_1079, ksk_1114, ksk_1115, \
                         lsi0_1091, lsi1_1091, lsk_1402, lsk_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1752[k] = f_20 * ksk_1114[k]
                    + f_4 * lsi0_1091[k]
                    - f_5 * lsi1_1091[k]
                    + f_3 * pc_y[k] * lsk_1402[k];

        t_1753[k] = f_20 * ksk_1115[k]
                    + f_3 * pc_y[k] * lsk_1403[k];

        t_1754[k] = f_16 * ksk_1079[k]
                    + f_1 * lsi0_1091[k]
                    - f_2 * lsi1_1091[k]
                    + f_3 * pc_z[k] * lsk_1403[k];
    }

#pragma omp simd aligned(t_1755, t_1756, t_1757, pc_x, lsi0_1092, lsi0_1093, lsi0_1094, \
                         lsi1_1092, lsi1_1093, lsi1_1094, lsk_1404, lsk_1405, \
                         lsk_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1755[k] = f_1 * lsi0_1092[k]
                    - f_2 * lsi1_1092[k]
                    + f_3 * pc_x[k] * lsk_1404[k];

        t_1756[k] = f_22 * lsi0_1093[k]
                    - f_23 * lsi1_1093[k]
                    + f_3 * pc_x[k] * lsk_1405[k];

        t_1757[k] = f_22 * lsi0_1094[k]
                    - f_23 * lsi1_1094[k]
                    + f_3 * pc_x[k] * lsk_1406[k];
    }

#pragma omp simd aligned(t_1758, t_1759, t_1760, pc_x, lsi0_1095, lsi0_1096, lsi0_1097, \
                         lsi1_1095, lsi1_1096, lsi1_1097, lsk_1407, lsk_1408, \
                         lsk_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1758[k] = f_12 * lsi0_1095[k]
                    - f_13 * lsi1_1095[k]
                    + f_3 * pc_x[k] * lsk_1407[k];

        t_1759[k] = f_12 * lsi0_1096[k]
                    - f_13 * lsi1_1096[k]
                    + f_3 * pc_x[k] * lsk_1408[k];

        t_1760[k] = f_12 * lsi0_1097[k]
                    - f_13 * lsi1_1097[k]
                    + f_3 * pc_x[k] * lsk_1409[k];
    }

#pragma omp simd aligned(t_1761, t_1762, t_1763, pc_x, lsi0_1098, lsi0_1099, lsi0_1100, \
                         lsi1_1098, lsi1_1099, lsi1_1100, lsk_1410, lsk_1411, \
                         lsk_1412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1761[k] = f_10 * lsi0_1098[k]
                    - f_11 * lsi1_1098[k]
                    + f_3 * pc_x[k] * lsk_1410[k];

        t_1762[k] = f_10 * lsi0_1099[k]
                    - f_11 * lsi1_1099[k]
                    + f_3 * pc_x[k] * lsk_1411[k];

        t_1763[k] = f_10 * lsi0_1100[k]
                    - f_11 * lsi1_1100[k]
                    + f_3 * pc_x[k] * lsk_1412[k];
    }

#pragma omp simd aligned(t_1764, t_1765, t_1766, pc_x, lsi0_1101, lsi0_1102, lsi0_1103, \
                         lsi1_1101, lsi1_1102, lsi1_1103, lsk_1413, lsk_1414, \
                         lsk_1415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1764[k] = f_10 * lsi0_1101[k]
                    - f_11 * lsi1_1101[k]
                    + f_3 * pc_x[k] * lsk_1413[k];

        t_1765[k] = f_8 * lsi0_1102[k]
                    - f_9 * lsi1_1102[k]
                    + f_3 * pc_x[k] * lsk_1414[k];

        t_1766[k] = f_8 * lsi0_1103[k]
                    - f_9 * lsi1_1103[k]
                    + f_3 * pc_x[k] * lsk_1415[k];
    }

#pragma omp simd aligned(t_1767, t_1768, t_1769, pc_x, lsi0_1104, lsi0_1105, lsi0_1106, \
                         lsi1_1104, lsi1_1105, lsi1_1106, lsk_1416, lsk_1417, \
                         lsk_1418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1767[k] = f_8 * lsi0_1104[k]
                    - f_9 * lsi1_1104[k]
                    + f_3 * pc_x[k] * lsk_1416[k];

        t_1768[k] = f_8 * lsi0_1105[k]
                    - f_9 * lsi1_1105[k]
                    + f_3 * pc_x[k] * lsk_1417[k];

        t_1769[k] = f_8 * lsi0_1106[k]
                    - f_9 * lsi1_1106[k]
                    + f_3 * pc_x[k] * lsk_1418[k];
    }

#pragma omp simd aligned(t_1770, t_1771, t_1772, pc_x, lsi0_1107, lsi0_1108, lsi0_1109, \
                         lsi1_1107, lsi1_1108, lsi1_1109, lsk_1419, lsk_1420, \
                         lsk_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1770[k] = f_6 * lsi0_1107[k]
                    - f_7 * lsi1_1107[k]
                    + f_3 * pc_x[k] * lsk_1419[k];

        t_1771[k] = f_6 * lsi0_1108[k]
                    - f_7 * lsi1_1108[k]
                    + f_3 * pc_x[k] * lsk_1420[k];

        t_1772[k] = f_6 * lsi0_1109[k]
                    - f_7 * lsi1_1109[k]
                    + f_3 * pc_x[k] * lsk_1421[k];
    }

#pragma omp simd aligned(t_1773, t_1774, t_1775, pc_x, lsi0_1110, lsi0_1111, lsi0_1112, \
                         lsi1_1110, lsi1_1111, lsi1_1112, lsk_1422, lsk_1423, \
                         lsk_1424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1773[k] = f_6 * lsi0_1110[k]
                    - f_7 * lsi1_1110[k]
                    + f_3 * pc_x[k] * lsk_1422[k];

        t_1774[k] = f_6 * lsi0_1111[k]
                    - f_7 * lsi1_1111[k]
                    + f_3 * pc_x[k] * lsk_1423[k];

        t_1775[k] = f_6 * lsi0_1112[k]
                    - f_7 * lsi1_1112[k]
                    + f_3 * pc_x[k] * lsk_1424[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, pc_x, lsi0_1113, lsi0_1114, lsi0_1115, \
                         lsi1_1113, lsi1_1114, lsi1_1115, lsk_1425, lsk_1426, \
                         lsk_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = f_4 * lsi0_1113[k]
                    - f_5 * lsi1_1113[k]
                    + f_3 * pc_x[k] * lsk_1425[k];

        t_1777[k] = f_4 * lsi0_1114[k]
                    - f_5 * lsi1_1114[k]
                    + f_3 * pc_x[k] * lsk_1426[k];

        t_1778[k] = f_4 * lsi0_1115[k]
                    - f_5 * lsi1_1115[k]
                    + f_3 * pc_x[k] * lsk_1427[k];
    }

#pragma omp simd aligned(t_1779, t_1780, t_1781, pc_x, lsi0_1116, lsi0_1117, lsi0_1118, \
                         lsi1_1116, lsi1_1117, lsi1_1118, lsk_1428, lsk_1429, \
                         lsk_1430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1779[k] = f_4 * lsi0_1116[k]
                    - f_5 * lsi1_1116[k]
                    + f_3 * pc_x[k] * lsk_1428[k];

        t_1780[k] = f_4 * lsi0_1117[k]
                    - f_5 * lsi1_1117[k]
                    + f_3 * pc_x[k] * lsk_1429[k];

        t_1781[k] = f_4 * lsi0_1118[k]
                    - f_5 * lsi1_1118[k]
                    + f_3 * pc_x[k] * lsk_1430[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, t_1785, t_1786, t_1787, pc_x, lsi0_1119, \
                         lsi1_1119, lsk_1431, lsk_1432, lsk_1433, lsk_1434, lsk_1435, \
                         lsk_1436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = f_4 * lsi0_1119[k]
                    - f_5 * lsi1_1119[k]
                    + f_3 * pc_x[k] * lsk_1431[k];

        t_1783[k] = f_3 * pc_x[k] * lsk_1432[k];

        t_1784[k] = f_3 * pc_x[k] * lsk_1433[k];

        t_1785[k] = f_3 * pc_x[k] * lsk_1434[k];

        t_1786[k] = f_3 * pc_x[k] * lsk_1435[k];

        t_1787[k] = f_3 * pc_x[k] * lsk_1436[k];
    }

#pragma omp simd aligned(t_1788, t_1789, t_1790, t_1791, t_1792, pc_x, pc_y, pc_z, ksk_1108, \
                         ksk_1144, lsi0_1113, lsi1_1113, lsk_1432, lsk_1437, lsk_1438, \
                         lsk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1788[k] = f_3 * pc_x[k] * lsk_1437[k];

        t_1789[k] = f_3 * pc_x[k] * lsk_1438[k];

        t_1790[k] = f_3 * pc_x[k] * lsk_1439[k];

        t_1791[k] = f_19 * ksk_1144[k]
                    + f_1 * lsi0_1113[k]
                    - f_2 * lsi1_1113[k]
                    + f_3 * pc_y[k] * lsk_1432[k];

        t_1792[k] = f_17 * ksk_1108[k]
                    + f_3 * pc_z[k] * lsk_1432[k];
    }

#pragma omp simd aligned(t_1793, t_1794, t_1795, pc_y, ksk_1146, ksk_1147, ksk_1148, \
                         lsi0_1115, lsi0_1116, lsi0_1117, lsi1_1115, lsi1_1116, lsi1_1117, \
                         lsk_1434, lsk_1435, lsk_1436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1793[k] = f_19 * ksk_1146[k]
                    + f_12 * lsi0_1115[k]
                    - f_13 * lsi1_1115[k]
                    + f_3 * pc_y[k] * lsk_1434[k];

        t_1794[k] = f_19 * ksk_1147[k]
                    + f_10 * lsi0_1116[k]
                    - f_11 * lsi1_1116[k]
                    + f_3 * pc_y[k] * lsk_1435[k];

        t_1795[k] = f_19 * ksk_1148[k]
                    + f_8 * lsi0_1117[k]
                    - f_9 * lsi1_1117[k]
                    + f_3 * pc_y[k] * lsk_1436[k];
    }

#pragma omp simd aligned(t_1796, t_1797, t_1798, pc_y, ksk_1149, ksk_1150, ksk_1151, \
                         lsi0_1118, lsi0_1119, lsi1_1118, lsi1_1119, lsk_1437, lsk_1438, \
                         lsk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1796[k] = f_19 * ksk_1149[k]
                    + f_6 * lsi0_1118[k]
                    - f_7 * lsi1_1118[k]
                    + f_3 * pc_y[k] * lsk_1437[k];

        t_1797[k] = f_19 * ksk_1150[k]
                    + f_4 * lsi0_1119[k]
                    - f_5 * lsi1_1119[k]
                    + f_3 * pc_y[k] * lsk_1438[k];

        t_1798[k] = f_19 * ksk_1151[k]
                    + f_3 * pc_y[k] * lsk_1439[k];
    }

#pragma omp simd aligned(t_1799, t_1800, t_1801, pc_x, pc_z, ksk_1115, lsi0_1119, lsi0_1120, \
                         lsi0_1121, lsi1_1119, lsi1_1120, lsi1_1121, lsk_1439, lsk_1440, \
                         lsk_1441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1799[k] = f_17 * ksk_1115[k]
                    + f_1 * lsi0_1119[k]
                    - f_2 * lsi1_1119[k]
                    + f_3 * pc_z[k] * lsk_1439[k];

        t_1800[k] = f_1 * lsi0_1120[k]
                    - f_2 * lsi1_1120[k]
                    + f_3 * pc_x[k] * lsk_1440[k];

        t_1801[k] = f_22 * lsi0_1121[k]
                    - f_23 * lsi1_1121[k]
                    + f_3 * pc_x[k] * lsk_1441[k];
    }

#pragma omp simd aligned(t_1802, t_1803, t_1804, pc_x, lsi0_1122, lsi0_1123, lsi0_1124, \
                         lsi1_1122, lsi1_1123, lsi1_1124, lsk_1442, lsk_1443, \
                         lsk_1444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1802[k] = f_22 * lsi0_1122[k]
                    - f_23 * lsi1_1122[k]
                    + f_3 * pc_x[k] * lsk_1442[k];

        t_1803[k] = f_12 * lsi0_1123[k]
                    - f_13 * lsi1_1123[k]
                    + f_3 * pc_x[k] * lsk_1443[k];

        t_1804[k] = f_12 * lsi0_1124[k]
                    - f_13 * lsi1_1124[k]
                    + f_3 * pc_x[k] * lsk_1444[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, pc_x, lsi0_1125, lsi0_1126, lsi0_1127, \
                         lsi1_1125, lsi1_1126, lsi1_1127, lsk_1445, lsk_1446, \
                         lsk_1447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = f_12 * lsi0_1125[k]
                    - f_13 * lsi1_1125[k]
                    + f_3 * pc_x[k] * lsk_1445[k];

        t_1806[k] = f_10 * lsi0_1126[k]
                    - f_11 * lsi1_1126[k]
                    + f_3 * pc_x[k] * lsk_1446[k];

        t_1807[k] = f_10 * lsi0_1127[k]
                    - f_11 * lsi1_1127[k]
                    + f_3 * pc_x[k] * lsk_1447[k];
    }

#pragma omp simd aligned(t_1808, t_1809, t_1810, pc_x, lsi0_1128, lsi0_1129, lsi0_1130, \
                         lsi1_1128, lsi1_1129, lsi1_1130, lsk_1448, lsk_1449, \
                         lsk_1450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1808[k] = f_10 * lsi0_1128[k]
                    - f_11 * lsi1_1128[k]
                    + f_3 * pc_x[k] * lsk_1448[k];

        t_1809[k] = f_10 * lsi0_1129[k]
                    - f_11 * lsi1_1129[k]
                    + f_3 * pc_x[k] * lsk_1449[k];

        t_1810[k] = f_8 * lsi0_1130[k]
                    - f_9 * lsi1_1130[k]
                    + f_3 * pc_x[k] * lsk_1450[k];
    }

#pragma omp simd aligned(t_1811, t_1812, t_1813, pc_x, lsi0_1131, lsi0_1132, lsi0_1133, \
                         lsi1_1131, lsi1_1132, lsi1_1133, lsk_1451, lsk_1452, \
                         lsk_1453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1811[k] = f_8 * lsi0_1131[k]
                    - f_9 * lsi1_1131[k]
                    + f_3 * pc_x[k] * lsk_1451[k];

        t_1812[k] = f_8 * lsi0_1132[k]
                    - f_9 * lsi1_1132[k]
                    + f_3 * pc_x[k] * lsk_1452[k];

        t_1813[k] = f_8 * lsi0_1133[k]
                    - f_9 * lsi1_1133[k]
                    + f_3 * pc_x[k] * lsk_1453[k];
    }

#pragma omp simd aligned(t_1814, t_1815, t_1816, pc_x, lsi0_1134, lsi0_1135, lsi0_1136, \
                         lsi1_1134, lsi1_1135, lsi1_1136, lsk_1454, lsk_1455, \
                         lsk_1456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1814[k] = f_8 * lsi0_1134[k]
                    - f_9 * lsi1_1134[k]
                    + f_3 * pc_x[k] * lsk_1454[k];

        t_1815[k] = f_6 * lsi0_1135[k]
                    - f_7 * lsi1_1135[k]
                    + f_3 * pc_x[k] * lsk_1455[k];

        t_1816[k] = f_6 * lsi0_1136[k]
                    - f_7 * lsi1_1136[k]
                    + f_3 * pc_x[k] * lsk_1456[k];
    }

#pragma omp simd aligned(t_1817, t_1818, t_1819, pc_x, lsi0_1137, lsi0_1138, lsi0_1139, \
                         lsi1_1137, lsi1_1138, lsi1_1139, lsk_1457, lsk_1458, \
                         lsk_1459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1817[k] = f_6 * lsi0_1137[k]
                    - f_7 * lsi1_1137[k]
                    + f_3 * pc_x[k] * lsk_1457[k];

        t_1818[k] = f_6 * lsi0_1138[k]
                    - f_7 * lsi1_1138[k]
                    + f_3 * pc_x[k] * lsk_1458[k];

        t_1819[k] = f_6 * lsi0_1139[k]
                    - f_7 * lsi1_1139[k]
                    + f_3 * pc_x[k] * lsk_1459[k];
    }

#pragma omp simd aligned(t_1820, t_1821, t_1822, pc_x, lsi0_1140, lsi0_1141, lsi0_1142, \
                         lsi1_1140, lsi1_1141, lsi1_1142, lsk_1460, lsk_1461, \
                         lsk_1462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1820[k] = f_6 * lsi0_1140[k]
                    - f_7 * lsi1_1140[k]
                    + f_3 * pc_x[k] * lsk_1460[k];

        t_1821[k] = f_4 * lsi0_1141[k]
                    - f_5 * lsi1_1141[k]
                    + f_3 * pc_x[k] * lsk_1461[k];

        t_1822[k] = f_4 * lsi0_1142[k]
                    - f_5 * lsi1_1142[k]
                    + f_3 * pc_x[k] * lsk_1462[k];
    }

#pragma omp simd aligned(t_1823, t_1824, t_1825, pc_x, lsi0_1143, lsi0_1144, lsi0_1145, \
                         lsi1_1143, lsi1_1144, lsi1_1145, lsk_1463, lsk_1464, \
                         lsk_1465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1823[k] = f_4 * lsi0_1143[k]
                    - f_5 * lsi1_1143[k]
                    + f_3 * pc_x[k] * lsk_1463[k];

        t_1824[k] = f_4 * lsi0_1144[k]
                    - f_5 * lsi1_1144[k]
                    + f_3 * pc_x[k] * lsk_1464[k];

        t_1825[k] = f_4 * lsi0_1145[k]
                    - f_5 * lsi1_1145[k]
                    + f_3 * pc_x[k] * lsk_1465[k];
    }

#pragma omp simd aligned(t_1826, t_1827, t_1828, t_1829, t_1830, pc_x, lsi0_1146, lsi0_1147, \
                         lsi1_1146, lsi1_1147, lsk_1466, lsk_1467, lsk_1468, lsk_1469, \
                         lsk_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1826[k] = f_4 * lsi0_1146[k]
                    - f_5 * lsi1_1146[k]
                    + f_3 * pc_x[k] * lsk_1466[k];

        t_1827[k] = f_4 * lsi0_1147[k]
                    - f_5 * lsi1_1147[k]
                    + f_3 * pc_x[k] * lsk_1467[k];

        t_1828[k] = f_3 * pc_x[k] * lsk_1468[k];

        t_1829[k] = f_3 * pc_x[k] * lsk_1469[k];

        t_1830[k] = f_3 * pc_x[k] * lsk_1470[k];
    }

#pragma omp simd aligned(t_1831, t_1832, t_1833, t_1834, t_1835, pc_x, lsk_1471, lsk_1472, \
                         lsk_1473, lsk_1474, lsk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1831[k] = f_3 * pc_x[k] * lsk_1471[k];

        t_1832[k] = f_3 * pc_x[k] * lsk_1472[k];

        t_1833[k] = f_3 * pc_x[k] * lsk_1473[k];

        t_1834[k] = f_3 * pc_x[k] * lsk_1474[k];

        t_1835[k] = f_3 * pc_x[k] * lsk_1475[k];
    }

#pragma omp simd aligned(t_1836, t_1837, t_1838, pc_y, pc_z, ksk_1144, ksk_1180, ksk_1182, \
                         lsi0_1141, lsi0_1143, lsi1_1141, lsi1_1143, lsk_1468, \
                         lsk_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1836[k] = f_18 * ksk_1180[k]
                    + f_1 * lsi0_1141[k]
                    - f_2 * lsi1_1141[k]
                    + f_3 * pc_y[k] * lsk_1468[k];

        t_1837[k] = f_18 * ksk_1144[k]
                    + f_3 * pc_z[k] * lsk_1468[k];

        t_1838[k] = f_18 * ksk_1182[k]
                    + f_12 * lsi0_1143[k]
                    - f_13 * lsi1_1143[k]
                    + f_3 * pc_y[k] * lsk_1470[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksl0,
                                                           const size_t ksk, const size_t ksl1,
                                                           const size_t lsi0, const size_t lsi1,
                                                           const size_t lsk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_1575 = buffer.data(ksl0 + 1575);
    const auto *ksl0_1577 = buffer.data(ksl0 + 1577);
    const auto *ksl0_1580 = buffer.data(ksl0 + 1580);
    const auto *ksl0_1584 = buffer.data(ksl0 + 1584);
    const auto *ksl0_1589 = buffer.data(ksl0 + 1589);

    const auto *ksk_1151 = buffer.data(ksk + 1151);
    const auto *ksk_1180 = buffer.data(ksk + 1180);
    const auto *ksk_1183 = buffer.data(ksk + 1183);
    const auto *ksk_1184 = buffer.data(ksk + 1184);
    const auto *ksk_1185 = buffer.data(ksk + 1185);
    const auto *ksk_1186 = buffer.data(ksk + 1186);
    const auto *ksk_1187 = buffer.data(ksk + 1187);
    const auto *ksk_1216 = buffer.data(ksk + 1216);
    const auto *ksk_1218 = buffer.data(ksk + 1218);
    const auto *ksk_1219 = buffer.data(ksk + 1219);
    const auto *ksk_1220 = buffer.data(ksk + 1220);
    const auto *ksk_1221 = buffer.data(ksk + 1221);
    const auto *ksk_1222 = buffer.data(ksk + 1222);
    const auto *ksk_1223 = buffer.data(ksk + 1223);
    const auto *ksk_1252 = buffer.data(ksk + 1252);
    const auto *ksk_1254 = buffer.data(ksk + 1254);
    const auto *ksk_1255 = buffer.data(ksk + 1255);
    const auto *ksk_1256 = buffer.data(ksk + 1256);
    const auto *ksk_1257 = buffer.data(ksk + 1257);
    const auto *ksk_1258 = buffer.data(ksk + 1258);
    const auto *ksk_1259 = buffer.data(ksk + 1259);

    const auto *ksl1_1575 = buffer.data(ksl1 + 1575);
    const auto *ksl1_1577 = buffer.data(ksl1 + 1577);
    const auto *ksl1_1580 = buffer.data(ksl1 + 1580);
    const auto *ksl1_1584 = buffer.data(ksl1 + 1584);
    const auto *ksl1_1589 = buffer.data(ksl1 + 1589);

    const auto *lsi0_1144 = buffer.data(lsi0 + 1144);
    const auto *lsi0_1145 = buffer.data(lsi0 + 1145);
    const auto *lsi0_1146 = buffer.data(lsi0 + 1146);
    const auto *lsi0_1147 = buffer.data(lsi0 + 1147);
    const auto *lsi0_1148 = buffer.data(lsi0 + 1148);
    const auto *lsi0_1149 = buffer.data(lsi0 + 1149);
    const auto *lsi0_1150 = buffer.data(lsi0 + 1150);
    const auto *lsi0_1151 = buffer.data(lsi0 + 1151);
    const auto *lsi0_1152 = buffer.data(lsi0 + 1152);
    const auto *lsi0_1153 = buffer.data(lsi0 + 1153);
    const auto *lsi0_1154 = buffer.data(lsi0 + 1154);
    const auto *lsi0_1155 = buffer.data(lsi0 + 1155);
    const auto *lsi0_1156 = buffer.data(lsi0 + 1156);
    const auto *lsi0_1157 = buffer.data(lsi0 + 1157);
    const auto *lsi0_1158 = buffer.data(lsi0 + 1158);
    const auto *lsi0_1159 = buffer.data(lsi0 + 1159);
    const auto *lsi0_1160 = buffer.data(lsi0 + 1160);
    const auto *lsi0_1161 = buffer.data(lsi0 + 1161);
    const auto *lsi0_1162 = buffer.data(lsi0 + 1162);
    const auto *lsi0_1163 = buffer.data(lsi0 + 1163);
    const auto *lsi0_1164 = buffer.data(lsi0 + 1164);
    const auto *lsi0_1165 = buffer.data(lsi0 + 1165);
    const auto *lsi0_1166 = buffer.data(lsi0 + 1166);
    const auto *lsi0_1167 = buffer.data(lsi0 + 1167);
    const auto *lsi0_1168 = buffer.data(lsi0 + 1168);
    const auto *lsi0_1169 = buffer.data(lsi0 + 1169);
    const auto *lsi0_1170 = buffer.data(lsi0 + 1170);
    const auto *lsi0_1171 = buffer.data(lsi0 + 1171);
    const auto *lsi0_1172 = buffer.data(lsi0 + 1172);
    const auto *lsi0_1173 = buffer.data(lsi0 + 1173);
    const auto *lsi0_1174 = buffer.data(lsi0 + 1174);
    const auto *lsi0_1175 = buffer.data(lsi0 + 1175);
    const auto *lsi0_1176 = buffer.data(lsi0 + 1176);
    const auto *lsi0_1177 = buffer.data(lsi0 + 1177);
    const auto *lsi0_1178 = buffer.data(lsi0 + 1178);
    const auto *lsi0_1179 = buffer.data(lsi0 + 1179);
    const auto *lsi0_1180 = buffer.data(lsi0 + 1180);
    const auto *lsi0_1181 = buffer.data(lsi0 + 1181);
    const auto *lsi0_1182 = buffer.data(lsi0 + 1182);
    const auto *lsi0_1183 = buffer.data(lsi0 + 1183);
    const auto *lsi0_1184 = buffer.data(lsi0 + 1184);
    const auto *lsi0_1185 = buffer.data(lsi0 + 1185);
    const auto *lsi0_1186 = buffer.data(lsi0 + 1186);
    const auto *lsi0_1187 = buffer.data(lsi0 + 1187);
    const auto *lsi0_1188 = buffer.data(lsi0 + 1188);
    const auto *lsi0_1189 = buffer.data(lsi0 + 1189);
    const auto *lsi0_1190 = buffer.data(lsi0 + 1190);
    const auto *lsi0_1191 = buffer.data(lsi0 + 1191);
    const auto *lsi0_1192 = buffer.data(lsi0 + 1192);
    const auto *lsi0_1193 = buffer.data(lsi0 + 1193);
    const auto *lsi0_1194 = buffer.data(lsi0 + 1194);
    const auto *lsi0_1195 = buffer.data(lsi0 + 1195);
    const auto *lsi0_1196 = buffer.data(lsi0 + 1196);
    const auto *lsi0_1197 = buffer.data(lsi0 + 1197);
    const auto *lsi0_1198 = buffer.data(lsi0 + 1198);
    const auto *lsi0_1199 = buffer.data(lsi0 + 1199);
    const auto *lsi0_1200 = buffer.data(lsi0 + 1200);
    const auto *lsi0_1201 = buffer.data(lsi0 + 1201);
    const auto *lsi0_1202 = buffer.data(lsi0 + 1202);
    const auto *lsi0_1203 = buffer.data(lsi0 + 1203);
    const auto *lsi0_1205 = buffer.data(lsi0 + 1205);
    const auto *lsi0_1207 = buffer.data(lsi0 + 1207);
    const auto *lsi0_1208 = buffer.data(lsi0 + 1208);
    const auto *lsi0_1210 = buffer.data(lsi0 + 1210);
    const auto *lsi0_1211 = buffer.data(lsi0 + 1211);
    const auto *lsi0_1212 = buffer.data(lsi0 + 1212);
    const auto *lsi0_1214 = buffer.data(lsi0 + 1214);
    const auto *lsi0_1215 = buffer.data(lsi0 + 1215);
    const auto *lsi0_1216 = buffer.data(lsi0 + 1216);
    const auto *lsi0_1217 = buffer.data(lsi0 + 1217);
    const auto *lsi0_1219 = buffer.data(lsi0 + 1219);

    const auto *lsi1_1144 = buffer.data(lsi1 + 1144);
    const auto *lsi1_1145 = buffer.data(lsi1 + 1145);
    const auto *lsi1_1146 = buffer.data(lsi1 + 1146);
    const auto *lsi1_1147 = buffer.data(lsi1 + 1147);
    const auto *lsi1_1148 = buffer.data(lsi1 + 1148);
    const auto *lsi1_1149 = buffer.data(lsi1 + 1149);
    const auto *lsi1_1150 = buffer.data(lsi1 + 1150);
    const auto *lsi1_1151 = buffer.data(lsi1 + 1151);
    const auto *lsi1_1152 = buffer.data(lsi1 + 1152);
    const auto *lsi1_1153 = buffer.data(lsi1 + 1153);
    const auto *lsi1_1154 = buffer.data(lsi1 + 1154);
    const auto *lsi1_1155 = buffer.data(lsi1 + 1155);
    const auto *lsi1_1156 = buffer.data(lsi1 + 1156);
    const auto *lsi1_1157 = buffer.data(lsi1 + 1157);
    const auto *lsi1_1158 = buffer.data(lsi1 + 1158);
    const auto *lsi1_1159 = buffer.data(lsi1 + 1159);
    const auto *lsi1_1160 = buffer.data(lsi1 + 1160);
    const auto *lsi1_1161 = buffer.data(lsi1 + 1161);
    const auto *lsi1_1162 = buffer.data(lsi1 + 1162);
    const auto *lsi1_1163 = buffer.data(lsi1 + 1163);
    const auto *lsi1_1164 = buffer.data(lsi1 + 1164);
    const auto *lsi1_1165 = buffer.data(lsi1 + 1165);
    const auto *lsi1_1166 = buffer.data(lsi1 + 1166);
    const auto *lsi1_1167 = buffer.data(lsi1 + 1167);
    const auto *lsi1_1168 = buffer.data(lsi1 + 1168);
    const auto *lsi1_1169 = buffer.data(lsi1 + 1169);
    const auto *lsi1_1170 = buffer.data(lsi1 + 1170);
    const auto *lsi1_1171 = buffer.data(lsi1 + 1171);
    const auto *lsi1_1172 = buffer.data(lsi1 + 1172);
    const auto *lsi1_1173 = buffer.data(lsi1 + 1173);
    const auto *lsi1_1174 = buffer.data(lsi1 + 1174);
    const auto *lsi1_1175 = buffer.data(lsi1 + 1175);
    const auto *lsi1_1176 = buffer.data(lsi1 + 1176);
    const auto *lsi1_1177 = buffer.data(lsi1 + 1177);
    const auto *lsi1_1178 = buffer.data(lsi1 + 1178);
    const auto *lsi1_1179 = buffer.data(lsi1 + 1179);
    const auto *lsi1_1180 = buffer.data(lsi1 + 1180);
    const auto *lsi1_1181 = buffer.data(lsi1 + 1181);
    const auto *lsi1_1182 = buffer.data(lsi1 + 1182);
    const auto *lsi1_1183 = buffer.data(lsi1 + 1183);
    const auto *lsi1_1184 = buffer.data(lsi1 + 1184);
    const auto *lsi1_1185 = buffer.data(lsi1 + 1185);
    const auto *lsi1_1186 = buffer.data(lsi1 + 1186);
    const auto *lsi1_1187 = buffer.data(lsi1 + 1187);
    const auto *lsi1_1188 = buffer.data(lsi1 + 1188);
    const auto *lsi1_1189 = buffer.data(lsi1 + 1189);
    const auto *lsi1_1190 = buffer.data(lsi1 + 1190);
    const auto *lsi1_1191 = buffer.data(lsi1 + 1191);
    const auto *lsi1_1192 = buffer.data(lsi1 + 1192);
    const auto *lsi1_1193 = buffer.data(lsi1 + 1193);
    const auto *lsi1_1194 = buffer.data(lsi1 + 1194);
    const auto *lsi1_1195 = buffer.data(lsi1 + 1195);
    const auto *lsi1_1196 = buffer.data(lsi1 + 1196);
    const auto *lsi1_1197 = buffer.data(lsi1 + 1197);
    const auto *lsi1_1198 = buffer.data(lsi1 + 1198);
    const auto *lsi1_1199 = buffer.data(lsi1 + 1199);
    const auto *lsi1_1200 = buffer.data(lsi1 + 1200);
    const auto *lsi1_1201 = buffer.data(lsi1 + 1201);
    const auto *lsi1_1202 = buffer.data(lsi1 + 1202);
    const auto *lsi1_1203 = buffer.data(lsi1 + 1203);
    const auto *lsi1_1205 = buffer.data(lsi1 + 1205);
    const auto *lsi1_1207 = buffer.data(lsi1 + 1207);
    const auto *lsi1_1208 = buffer.data(lsi1 + 1208);
    const auto *lsi1_1210 = buffer.data(lsi1 + 1210);
    const auto *lsi1_1211 = buffer.data(lsi1 + 1211);
    const auto *lsi1_1212 = buffer.data(lsi1 + 1212);
    const auto *lsi1_1214 = buffer.data(lsi1 + 1214);
    const auto *lsi1_1215 = buffer.data(lsi1 + 1215);
    const auto *lsi1_1216 = buffer.data(lsi1 + 1216);
    const auto *lsi1_1217 = buffer.data(lsi1 + 1217);
    const auto *lsi1_1219 = buffer.data(lsi1 + 1219);

    const auto *lsk_1471 = buffer.data(lsk + 1471);
    const auto *lsk_1472 = buffer.data(lsk + 1472);
    const auto *lsk_1473 = buffer.data(lsk + 1473);
    const auto *lsk_1474 = buffer.data(lsk + 1474);
    const auto *lsk_1475 = buffer.data(lsk + 1475);
    const auto *lsk_1476 = buffer.data(lsk + 1476);
    const auto *lsk_1477 = buffer.data(lsk + 1477);
    const auto *lsk_1478 = buffer.data(lsk + 1478);
    const auto *lsk_1479 = buffer.data(lsk + 1479);
    const auto *lsk_1480 = buffer.data(lsk + 1480);
    const auto *lsk_1481 = buffer.data(lsk + 1481);
    const auto *lsk_1482 = buffer.data(lsk + 1482);
    const auto *lsk_1483 = buffer.data(lsk + 1483);
    const auto *lsk_1484 = buffer.data(lsk + 1484);
    const auto *lsk_1485 = buffer.data(lsk + 1485);
    const auto *lsk_1486 = buffer.data(lsk + 1486);
    const auto *lsk_1487 = buffer.data(lsk + 1487);
    const auto *lsk_1488 = buffer.data(lsk + 1488);
    const auto *lsk_1489 = buffer.data(lsk + 1489);
    const auto *lsk_1490 = buffer.data(lsk + 1490);
    const auto *lsk_1491 = buffer.data(lsk + 1491);
    const auto *lsk_1492 = buffer.data(lsk + 1492);
    const auto *lsk_1493 = buffer.data(lsk + 1493);
    const auto *lsk_1494 = buffer.data(lsk + 1494);
    const auto *lsk_1495 = buffer.data(lsk + 1495);
    const auto *lsk_1496 = buffer.data(lsk + 1496);
    const auto *lsk_1497 = buffer.data(lsk + 1497);
    const auto *lsk_1498 = buffer.data(lsk + 1498);
    const auto *lsk_1499 = buffer.data(lsk + 1499);
    const auto *lsk_1500 = buffer.data(lsk + 1500);
    const auto *lsk_1501 = buffer.data(lsk + 1501);
    const auto *lsk_1502 = buffer.data(lsk + 1502);
    const auto *lsk_1503 = buffer.data(lsk + 1503);
    const auto *lsk_1504 = buffer.data(lsk + 1504);
    const auto *lsk_1505 = buffer.data(lsk + 1505);
    const auto *lsk_1506 = buffer.data(lsk + 1506);
    const auto *lsk_1507 = buffer.data(lsk + 1507);
    const auto *lsk_1508 = buffer.data(lsk + 1508);
    const auto *lsk_1509 = buffer.data(lsk + 1509);
    const auto *lsk_1510 = buffer.data(lsk + 1510);
    const auto *lsk_1511 = buffer.data(lsk + 1511);
    const auto *lsk_1512 = buffer.data(lsk + 1512);
    const auto *lsk_1513 = buffer.data(lsk + 1513);
    const auto *lsk_1514 = buffer.data(lsk + 1514);
    const auto *lsk_1515 = buffer.data(lsk + 1515);
    const auto *lsk_1516 = buffer.data(lsk + 1516);
    const auto *lsk_1517 = buffer.data(lsk + 1517);
    const auto *lsk_1518 = buffer.data(lsk + 1518);
    const auto *lsk_1519 = buffer.data(lsk + 1519);
    const auto *lsk_1520 = buffer.data(lsk + 1520);
    const auto *lsk_1521 = buffer.data(lsk + 1521);
    const auto *lsk_1522 = buffer.data(lsk + 1522);
    const auto *lsk_1523 = buffer.data(lsk + 1523);
    const auto *lsk_1524 = buffer.data(lsk + 1524);
    const auto *lsk_1525 = buffer.data(lsk + 1525);
    const auto *lsk_1526 = buffer.data(lsk + 1526);
    const auto *lsk_1527 = buffer.data(lsk + 1527);
    const auto *lsk_1528 = buffer.data(lsk + 1528);
    const auto *lsk_1529 = buffer.data(lsk + 1529);
    const auto *lsk_1530 = buffer.data(lsk + 1530);
    const auto *lsk_1531 = buffer.data(lsk + 1531);
    const auto *lsk_1532 = buffer.data(lsk + 1532);
    const auto *lsk_1533 = buffer.data(lsk + 1533);
    const auto *lsk_1534 = buffer.data(lsk + 1534);
    const auto *lsk_1535 = buffer.data(lsk + 1535);
    const auto *lsk_1536 = buffer.data(lsk + 1536);
    const auto *lsk_1537 = buffer.data(lsk + 1537);
    const auto *lsk_1538 = buffer.data(lsk + 1538);
    const auto *lsk_1539 = buffer.data(lsk + 1539);
    const auto *lsk_1540 = buffer.data(lsk + 1540);
    const auto *lsk_1541 = buffer.data(lsk + 1541);
    const auto *lsk_1542 = buffer.data(lsk + 1542);
    const auto *lsk_1543 = buffer.data(lsk + 1543);
    const auto *lsk_1544 = buffer.data(lsk + 1544);
    const auto *lsk_1545 = buffer.data(lsk + 1545);
    const auto *lsk_1546 = buffer.data(lsk + 1546);
    const auto *lsk_1547 = buffer.data(lsk + 1547);
    const auto *lsk_1549 = buffer.data(lsk + 1549);
    const auto *lsk_1551 = buffer.data(lsk + 1551);
    const auto *lsk_1552 = buffer.data(lsk + 1552);
    const auto *lsk_1554 = buffer.data(lsk + 1554);
    const auto *lsk_1555 = buffer.data(lsk + 1555);
    const auto *lsk_1556 = buffer.data(lsk + 1556);
    const auto *lsk_1558 = buffer.data(lsk + 1558);
    const auto *lsk_1559 = buffer.data(lsk + 1559);
    const auto *lsk_1560 = buffer.data(lsk + 1560);
    const auto *lsk_1561 = buffer.data(lsk + 1561);
    const auto *lsk_1563 = buffer.data(lsk + 1563);

#pragma omp simd aligned(t_1839, t_1840, t_1841, pc_y, ksk_1183, ksk_1184, ksk_1185, \
                         lsi0_1144, lsi0_1145, lsi0_1146, lsi1_1144, lsi1_1145, lsi1_1146, \
                         lsk_1471, lsk_1472, lsk_1473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1839[k] = f_18 * ksk_1183[k]
                    + f_10 * lsi0_1144[k]
                    - f_11 * lsi1_1144[k]
                    + f_3 * pc_y[k] * lsk_1471[k];

        t_1840[k] = f_18 * ksk_1184[k]
                    + f_8 * lsi0_1145[k]
                    - f_9 * lsi1_1145[k]
                    + f_3 * pc_y[k] * lsk_1472[k];

        t_1841[k] = f_18 * ksk_1185[k]
                    + f_6 * lsi0_1146[k]
                    - f_7 * lsi1_1146[k]
                    + f_3 * pc_y[k] * lsk_1473[k];
    }

#pragma omp simd aligned(t_1842, t_1843, t_1844, pc_y, pc_z, ksk_1151, ksk_1186, ksk_1187, \
                         lsi0_1147, lsi1_1147, lsk_1474, lsk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1842[k] = f_18 * ksk_1186[k]
                    + f_4 * lsi0_1147[k]
                    - f_5 * lsi1_1147[k]
                    + f_3 * pc_y[k] * lsk_1474[k];

        t_1843[k] = f_18 * ksk_1187[k]
                    + f_3 * pc_y[k] * lsk_1475[k];

        t_1844[k] = f_18 * ksk_1151[k]
                    + f_1 * lsi0_1147[k]
                    - f_2 * lsi1_1147[k]
                    + f_3 * pc_z[k] * lsk_1475[k];
    }

#pragma omp simd aligned(t_1845, t_1846, t_1847, pc_x, lsi0_1148, lsi0_1149, lsi0_1150, \
                         lsi1_1148, lsi1_1149, lsi1_1150, lsk_1476, lsk_1477, \
                         lsk_1478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1845[k] = f_1 * lsi0_1148[k]
                    - f_2 * lsi1_1148[k]
                    + f_3 * pc_x[k] * lsk_1476[k];

        t_1846[k] = f_22 * lsi0_1149[k]
                    - f_23 * lsi1_1149[k]
                    + f_3 * pc_x[k] * lsk_1477[k];

        t_1847[k] = f_22 * lsi0_1150[k]
                    - f_23 * lsi1_1150[k]
                    + f_3 * pc_x[k] * lsk_1478[k];
    }

#pragma omp simd aligned(t_1848, t_1849, t_1850, pc_x, lsi0_1151, lsi0_1152, lsi0_1153, \
                         lsi1_1151, lsi1_1152, lsi1_1153, lsk_1479, lsk_1480, \
                         lsk_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1848[k] = f_12 * lsi0_1151[k]
                    - f_13 * lsi1_1151[k]
                    + f_3 * pc_x[k] * lsk_1479[k];

        t_1849[k] = f_12 * lsi0_1152[k]
                    - f_13 * lsi1_1152[k]
                    + f_3 * pc_x[k] * lsk_1480[k];

        t_1850[k] = f_12 * lsi0_1153[k]
                    - f_13 * lsi1_1153[k]
                    + f_3 * pc_x[k] * lsk_1481[k];
    }

#pragma omp simd aligned(t_1851, t_1852, t_1853, pc_x, lsi0_1154, lsi0_1155, lsi0_1156, \
                         lsi1_1154, lsi1_1155, lsi1_1156, lsk_1482, lsk_1483, \
                         lsk_1484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1851[k] = f_10 * lsi0_1154[k]
                    - f_11 * lsi1_1154[k]
                    + f_3 * pc_x[k] * lsk_1482[k];

        t_1852[k] = f_10 * lsi0_1155[k]
                    - f_11 * lsi1_1155[k]
                    + f_3 * pc_x[k] * lsk_1483[k];

        t_1853[k] = f_10 * lsi0_1156[k]
                    - f_11 * lsi1_1156[k]
                    + f_3 * pc_x[k] * lsk_1484[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, pc_x, lsi0_1157, lsi0_1158, lsi0_1159, \
                         lsi1_1157, lsi1_1158, lsi1_1159, lsk_1485, lsk_1486, \
                         lsk_1487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = f_10 * lsi0_1157[k]
                    - f_11 * lsi1_1157[k]
                    + f_3 * pc_x[k] * lsk_1485[k];

        t_1855[k] = f_8 * lsi0_1158[k]
                    - f_9 * lsi1_1158[k]
                    + f_3 * pc_x[k] * lsk_1486[k];

        t_1856[k] = f_8 * lsi0_1159[k]
                    - f_9 * lsi1_1159[k]
                    + f_3 * pc_x[k] * lsk_1487[k];
    }

#pragma omp simd aligned(t_1857, t_1858, t_1859, pc_x, lsi0_1160, lsi0_1161, lsi0_1162, \
                         lsi1_1160, lsi1_1161, lsi1_1162, lsk_1488, lsk_1489, \
                         lsk_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1857[k] = f_8 * lsi0_1160[k]
                    - f_9 * lsi1_1160[k]
                    + f_3 * pc_x[k] * lsk_1488[k];

        t_1858[k] = f_8 * lsi0_1161[k]
                    - f_9 * lsi1_1161[k]
                    + f_3 * pc_x[k] * lsk_1489[k];

        t_1859[k] = f_8 * lsi0_1162[k]
                    - f_9 * lsi1_1162[k]
                    + f_3 * pc_x[k] * lsk_1490[k];
    }

#pragma omp simd aligned(t_1860, t_1861, t_1862, pc_x, lsi0_1163, lsi0_1164, lsi0_1165, \
                         lsi1_1163, lsi1_1164, lsi1_1165, lsk_1491, lsk_1492, \
                         lsk_1493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1860[k] = f_6 * lsi0_1163[k]
                    - f_7 * lsi1_1163[k]
                    + f_3 * pc_x[k] * lsk_1491[k];

        t_1861[k] = f_6 * lsi0_1164[k]
                    - f_7 * lsi1_1164[k]
                    + f_3 * pc_x[k] * lsk_1492[k];

        t_1862[k] = f_6 * lsi0_1165[k]
                    - f_7 * lsi1_1165[k]
                    + f_3 * pc_x[k] * lsk_1493[k];
    }

#pragma omp simd aligned(t_1863, t_1864, t_1865, pc_x, lsi0_1166, lsi0_1167, lsi0_1168, \
                         lsi1_1166, lsi1_1167, lsi1_1168, lsk_1494, lsk_1495, \
                         lsk_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1863[k] = f_6 * lsi0_1166[k]
                    - f_7 * lsi1_1166[k]
                    + f_3 * pc_x[k] * lsk_1494[k];

        t_1864[k] = f_6 * lsi0_1167[k]
                    - f_7 * lsi1_1167[k]
                    + f_3 * pc_x[k] * lsk_1495[k];

        t_1865[k] = f_6 * lsi0_1168[k]
                    - f_7 * lsi1_1168[k]
                    + f_3 * pc_x[k] * lsk_1496[k];
    }

#pragma omp simd aligned(t_1866, t_1867, t_1868, pc_x, lsi0_1169, lsi0_1170, lsi0_1171, \
                         lsi1_1169, lsi1_1170, lsi1_1171, lsk_1497, lsk_1498, \
                         lsk_1499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1866[k] = f_4 * lsi0_1169[k]
                    - f_5 * lsi1_1169[k]
                    + f_3 * pc_x[k] * lsk_1497[k];

        t_1867[k] = f_4 * lsi0_1170[k]
                    - f_5 * lsi1_1170[k]
                    + f_3 * pc_x[k] * lsk_1498[k];

        t_1868[k] = f_4 * lsi0_1171[k]
                    - f_5 * lsi1_1171[k]
                    + f_3 * pc_x[k] * lsk_1499[k];
    }

#pragma omp simd aligned(t_1869, t_1870, t_1871, pc_x, lsi0_1172, lsi0_1173, lsi0_1174, \
                         lsi1_1172, lsi1_1173, lsi1_1174, lsk_1500, lsk_1501, \
                         lsk_1502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1869[k] = f_4 * lsi0_1172[k]
                    - f_5 * lsi1_1172[k]
                    + f_3 * pc_x[k] * lsk_1500[k];

        t_1870[k] = f_4 * lsi0_1173[k]
                    - f_5 * lsi1_1173[k]
                    + f_3 * pc_x[k] * lsk_1501[k];

        t_1871[k] = f_4 * lsi0_1174[k]
                    - f_5 * lsi1_1174[k]
                    + f_3 * pc_x[k] * lsk_1502[k];
    }

#pragma omp simd aligned(t_1872, t_1873, t_1874, t_1875, t_1876, t_1877, pc_x, lsi0_1175, \
                         lsi1_1175, lsk_1503, lsk_1504, lsk_1505, lsk_1506, lsk_1507, \
                         lsk_1508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1872[k] = f_4 * lsi0_1175[k]
                    - f_5 * lsi1_1175[k]
                    + f_3 * pc_x[k] * lsk_1503[k];

        t_1873[k] = f_3 * pc_x[k] * lsk_1504[k];

        t_1874[k] = f_3 * pc_x[k] * lsk_1505[k];

        t_1875[k] = f_3 * pc_x[k] * lsk_1506[k];

        t_1876[k] = f_3 * pc_x[k] * lsk_1507[k];

        t_1877[k] = f_3 * pc_x[k] * lsk_1508[k];
    }

#pragma omp simd aligned(t_1878, t_1879, t_1880, t_1881, t_1882, pc_x, pc_y, pc_z, ksk_1180, \
                         ksk_1216, lsi0_1169, lsi1_1169, lsk_1504, lsk_1509, lsk_1510, \
                         lsk_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1878[k] = f_3 * pc_x[k] * lsk_1509[k];

        t_1879[k] = f_3 * pc_x[k] * lsk_1510[k];

        t_1880[k] = f_3 * pc_x[k] * lsk_1511[k];

        t_1881[k] = f_17 * ksk_1216[k]
                    + f_1 * lsi0_1169[k]
                    - f_2 * lsi1_1169[k]
                    + f_3 * pc_y[k] * lsk_1504[k];

        t_1882[k] = f_19 * ksk_1180[k]
                    + f_3 * pc_z[k] * lsk_1504[k];
    }

#pragma omp simd aligned(t_1883, t_1884, t_1885, pc_y, ksk_1218, ksk_1219, ksk_1220, \
                         lsi0_1171, lsi0_1172, lsi0_1173, lsi1_1171, lsi1_1172, lsi1_1173, \
                         lsk_1506, lsk_1507, lsk_1508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1883[k] = f_17 * ksk_1218[k]
                    + f_12 * lsi0_1171[k]
                    - f_13 * lsi1_1171[k]
                    + f_3 * pc_y[k] * lsk_1506[k];

        t_1884[k] = f_17 * ksk_1219[k]
                    + f_10 * lsi0_1172[k]
                    - f_11 * lsi1_1172[k]
                    + f_3 * pc_y[k] * lsk_1507[k];

        t_1885[k] = f_17 * ksk_1220[k]
                    + f_8 * lsi0_1173[k]
                    - f_9 * lsi1_1173[k]
                    + f_3 * pc_y[k] * lsk_1508[k];
    }

#pragma omp simd aligned(t_1886, t_1887, t_1888, pc_y, ksk_1221, ksk_1222, ksk_1223, \
                         lsi0_1174, lsi0_1175, lsi1_1174, lsi1_1175, lsk_1509, lsk_1510, \
                         lsk_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1886[k] = f_17 * ksk_1221[k]
                    + f_6 * lsi0_1174[k]
                    - f_7 * lsi1_1174[k]
                    + f_3 * pc_y[k] * lsk_1509[k];

        t_1887[k] = f_17 * ksk_1222[k]
                    + f_4 * lsi0_1175[k]
                    - f_5 * lsi1_1175[k]
                    + f_3 * pc_y[k] * lsk_1510[k];

        t_1888[k] = f_17 * ksk_1223[k]
                    + f_3 * pc_y[k] * lsk_1511[k];
    }

#pragma omp simd aligned(t_1889, t_1890, t_1891, pc_x, pc_z, ksk_1187, lsi0_1175, lsi0_1176, \
                         lsi0_1177, lsi1_1175, lsi1_1176, lsi1_1177, lsk_1511, lsk_1512, \
                         lsk_1513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1889[k] = f_19 * ksk_1187[k]
                    + f_1 * lsi0_1175[k]
                    - f_2 * lsi1_1175[k]
                    + f_3 * pc_z[k] * lsk_1511[k];

        t_1890[k] = f_1 * lsi0_1176[k]
                    - f_2 * lsi1_1176[k]
                    + f_3 * pc_x[k] * lsk_1512[k];

        t_1891[k] = f_22 * lsi0_1177[k]
                    - f_23 * lsi1_1177[k]
                    + f_3 * pc_x[k] * lsk_1513[k];
    }

#pragma omp simd aligned(t_1892, t_1893, t_1894, pc_x, lsi0_1178, lsi0_1179, lsi0_1180, \
                         lsi1_1178, lsi1_1179, lsi1_1180, lsk_1514, lsk_1515, \
                         lsk_1516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1892[k] = f_22 * lsi0_1178[k]
                    - f_23 * lsi1_1178[k]
                    + f_3 * pc_x[k] * lsk_1514[k];

        t_1893[k] = f_12 * lsi0_1179[k]
                    - f_13 * lsi1_1179[k]
                    + f_3 * pc_x[k] * lsk_1515[k];

        t_1894[k] = f_12 * lsi0_1180[k]
                    - f_13 * lsi1_1180[k]
                    + f_3 * pc_x[k] * lsk_1516[k];
    }

#pragma omp simd aligned(t_1895, t_1896, t_1897, pc_x, lsi0_1181, lsi0_1182, lsi0_1183, \
                         lsi1_1181, lsi1_1182, lsi1_1183, lsk_1517, lsk_1518, \
                         lsk_1519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1895[k] = f_12 * lsi0_1181[k]
                    - f_13 * lsi1_1181[k]
                    + f_3 * pc_x[k] * lsk_1517[k];

        t_1896[k] = f_10 * lsi0_1182[k]
                    - f_11 * lsi1_1182[k]
                    + f_3 * pc_x[k] * lsk_1518[k];

        t_1897[k] = f_10 * lsi0_1183[k]
                    - f_11 * lsi1_1183[k]
                    + f_3 * pc_x[k] * lsk_1519[k];
    }

#pragma omp simd aligned(t_1898, t_1899, t_1900, pc_x, lsi0_1184, lsi0_1185, lsi0_1186, \
                         lsi1_1184, lsi1_1185, lsi1_1186, lsk_1520, lsk_1521, \
                         lsk_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1898[k] = f_10 * lsi0_1184[k]
                    - f_11 * lsi1_1184[k]
                    + f_3 * pc_x[k] * lsk_1520[k];

        t_1899[k] = f_10 * lsi0_1185[k]
                    - f_11 * lsi1_1185[k]
                    + f_3 * pc_x[k] * lsk_1521[k];

        t_1900[k] = f_8 * lsi0_1186[k]
                    - f_9 * lsi1_1186[k]
                    + f_3 * pc_x[k] * lsk_1522[k];
    }

#pragma omp simd aligned(t_1901, t_1902, t_1903, pc_x, lsi0_1187, lsi0_1188, lsi0_1189, \
                         lsi1_1187, lsi1_1188, lsi1_1189, lsk_1523, lsk_1524, \
                         lsk_1525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1901[k] = f_8 * lsi0_1187[k]
                    - f_9 * lsi1_1187[k]
                    + f_3 * pc_x[k] * lsk_1523[k];

        t_1902[k] = f_8 * lsi0_1188[k]
                    - f_9 * lsi1_1188[k]
                    + f_3 * pc_x[k] * lsk_1524[k];

        t_1903[k] = f_8 * lsi0_1189[k]
                    - f_9 * lsi1_1189[k]
                    + f_3 * pc_x[k] * lsk_1525[k];
    }

#pragma omp simd aligned(t_1904, t_1905, t_1906, pc_x, lsi0_1190, lsi0_1191, lsi0_1192, \
                         lsi1_1190, lsi1_1191, lsi1_1192, lsk_1526, lsk_1527, \
                         lsk_1528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1904[k] = f_8 * lsi0_1190[k]
                    - f_9 * lsi1_1190[k]
                    + f_3 * pc_x[k] * lsk_1526[k];

        t_1905[k] = f_6 * lsi0_1191[k]
                    - f_7 * lsi1_1191[k]
                    + f_3 * pc_x[k] * lsk_1527[k];

        t_1906[k] = f_6 * lsi0_1192[k]
                    - f_7 * lsi1_1192[k]
                    + f_3 * pc_x[k] * lsk_1528[k];
    }

#pragma omp simd aligned(t_1907, t_1908, t_1909, pc_x, lsi0_1193, lsi0_1194, lsi0_1195, \
                         lsi1_1193, lsi1_1194, lsi1_1195, lsk_1529, lsk_1530, \
                         lsk_1531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1907[k] = f_6 * lsi0_1193[k]
                    - f_7 * lsi1_1193[k]
                    + f_3 * pc_x[k] * lsk_1529[k];

        t_1908[k] = f_6 * lsi0_1194[k]
                    - f_7 * lsi1_1194[k]
                    + f_3 * pc_x[k] * lsk_1530[k];

        t_1909[k] = f_6 * lsi0_1195[k]
                    - f_7 * lsi1_1195[k]
                    + f_3 * pc_x[k] * lsk_1531[k];
    }

#pragma omp simd aligned(t_1910, t_1911, t_1912, pc_x, lsi0_1196, lsi0_1197, lsi0_1198, \
                         lsi1_1196, lsi1_1197, lsi1_1198, lsk_1532, lsk_1533, \
                         lsk_1534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1910[k] = f_6 * lsi0_1196[k]
                    - f_7 * lsi1_1196[k]
                    + f_3 * pc_x[k] * lsk_1532[k];

        t_1911[k] = f_4 * lsi0_1197[k]
                    - f_5 * lsi1_1197[k]
                    + f_3 * pc_x[k] * lsk_1533[k];

        t_1912[k] = f_4 * lsi0_1198[k]
                    - f_5 * lsi1_1198[k]
                    + f_3 * pc_x[k] * lsk_1534[k];
    }

#pragma omp simd aligned(t_1913, t_1914, t_1915, pc_x, lsi0_1199, lsi0_1200, lsi0_1201, \
                         lsi1_1199, lsi1_1200, lsi1_1201, lsk_1535, lsk_1536, \
                         lsk_1537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1913[k] = f_4 * lsi0_1199[k]
                    - f_5 * lsi1_1199[k]
                    + f_3 * pc_x[k] * lsk_1535[k];

        t_1914[k] = f_4 * lsi0_1200[k]
                    - f_5 * lsi1_1200[k]
                    + f_3 * pc_x[k] * lsk_1536[k];

        t_1915[k] = f_4 * lsi0_1201[k]
                    - f_5 * lsi1_1201[k]
                    + f_3 * pc_x[k] * lsk_1537[k];
    }

#pragma omp simd aligned(t_1916, t_1917, t_1918, t_1919, t_1920, pc_x, lsi0_1202, lsi0_1203, \
                         lsi1_1202, lsi1_1203, lsk_1538, lsk_1539, lsk_1540, lsk_1541, \
                         lsk_1542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1916[k] = f_4 * lsi0_1202[k]
                    - f_5 * lsi1_1202[k]
                    + f_3 * pc_x[k] * lsk_1538[k];

        t_1917[k] = f_4 * lsi0_1203[k]
                    - f_5 * lsi1_1203[k]
                    + f_3 * pc_x[k] * lsk_1539[k];

        t_1918[k] = f_3 * pc_x[k] * lsk_1540[k];

        t_1919[k] = f_3 * pc_x[k] * lsk_1541[k];

        t_1920[k] = f_3 * pc_x[k] * lsk_1542[k];
    }

#pragma omp simd aligned(t_1921, t_1922, t_1923, t_1924, t_1925, pc_x, lsk_1543, lsk_1544, \
                         lsk_1545, lsk_1546, lsk_1547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1921[k] = f_3 * pc_x[k] * lsk_1543[k];

        t_1922[k] = f_3 * pc_x[k] * lsk_1544[k];

        t_1923[k] = f_3 * pc_x[k] * lsk_1545[k];

        t_1924[k] = f_3 * pc_x[k] * lsk_1546[k];

        t_1925[k] = f_3 * pc_x[k] * lsk_1547[k];
    }

#pragma omp simd aligned(t_1926, t_1927, t_1928, pc_y, pc_z, ksk_1216, ksk_1252, ksk_1254, \
                         lsi0_1197, lsi0_1199, lsi1_1197, lsi1_1199, lsk_1540, \
                         lsk_1542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1926[k] = f_16 * ksk_1252[k]
                    + f_1 * lsi0_1197[k]
                    - f_2 * lsi1_1197[k]
                    + f_3 * pc_y[k] * lsk_1540[k];

        t_1927[k] = f_20 * ksk_1216[k]
                    + f_3 * pc_z[k] * lsk_1540[k];

        t_1928[k] = f_16 * ksk_1254[k]
                    + f_12 * lsi0_1199[k]
                    - f_13 * lsi1_1199[k]
                    + f_3 * pc_y[k] * lsk_1542[k];
    }

#pragma omp simd aligned(t_1929, t_1930, t_1931, pc_y, ksk_1255, ksk_1256, ksk_1257, \
                         lsi0_1200, lsi0_1201, lsi0_1202, lsi1_1200, lsi1_1201, lsi1_1202, \
                         lsk_1543, lsk_1544, lsk_1545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1929[k] = f_16 * ksk_1255[k]
                    + f_10 * lsi0_1200[k]
                    - f_11 * lsi1_1200[k]
                    + f_3 * pc_y[k] * lsk_1543[k];

        t_1930[k] = f_16 * ksk_1256[k]
                    + f_8 * lsi0_1201[k]
                    - f_9 * lsi1_1201[k]
                    + f_3 * pc_y[k] * lsk_1544[k];

        t_1931[k] = f_16 * ksk_1257[k]
                    + f_6 * lsi0_1202[k]
                    - f_7 * lsi1_1202[k]
                    + f_3 * pc_y[k] * lsk_1545[k];
    }

#pragma omp simd aligned(t_1932, t_1933, t_1934, t_1935, pa_y, pc_y, pc_z, ksl0_1575, \
                         ksk_1223, ksk_1258, ksk_1259, ksl1_1575, lsi0_1203, lsi1_1203, \
                         lsk_1546, lsk_1547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1932[k] = f_16 * ksk_1258[k]
                    + f_4 * lsi0_1203[k]
                    - f_5 * lsi1_1203[k]
                    + f_3 * pc_y[k] * lsk_1546[k];

        t_1933[k] = f_16 * ksk_1259[k]
                    + f_3 * pc_y[k] * lsk_1547[k];

        t_1934[k] = f_20 * ksk_1223[k]
                    + f_1 * lsi0_1203[k]
                    - f_2 * lsi1_1203[k]
                    + f_3 * pc_z[k] * lsk_1547[k];

        t_1935[k] = pa_y[k] * ksl0_1575[k]
                    - f_14 * pc_y[k] * ksl1_1575[k];
    }

#pragma omp simd aligned(t_1936, t_1937, t_1938, pa_y, pc_x, pc_y, ksl0_1577, ksl1_1577, \
                         lsi0_1205, lsi0_1207, lsi1_1205, lsi1_1207, lsk_1549, \
                         lsk_1551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1936[k] = f_22 * lsi0_1205[k]
                    - f_23 * lsi1_1205[k]
                    + f_3 * pc_x[k] * lsk_1549[k];

        t_1937[k] = pa_y[k] * ksl0_1577[k]
                    - f_14 * pc_y[k] * ksl1_1577[k];

        t_1938[k] = f_12 * lsi0_1207[k]
                    - f_13 * lsi1_1207[k]
                    + f_3 * pc_x[k] * lsk_1551[k];
    }

#pragma omp simd aligned(t_1939, t_1940, t_1941, pa_y, pc_x, pc_y, ksl0_1580, ksl1_1580, \
                         lsi0_1208, lsi0_1210, lsi1_1208, lsi1_1210, lsk_1552, \
                         lsk_1554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1939[k] = f_12 * lsi0_1208[k]
                    - f_13 * lsi1_1208[k]
                    + f_3 * pc_x[k] * lsk_1552[k];

        t_1940[k] = pa_y[k] * ksl0_1580[k]
                    - f_14 * pc_y[k] * ksl1_1580[k];

        t_1941[k] = f_10 * lsi0_1210[k]
                    - f_11 * lsi1_1210[k]
                    + f_3 * pc_x[k] * lsk_1554[k];
    }

#pragma omp simd aligned(t_1942, t_1943, t_1944, pa_y, pc_x, pc_y, ksl0_1584, ksl1_1584, \
                         lsi0_1211, lsi0_1212, lsi1_1211, lsi1_1212, lsk_1555, \
                         lsk_1556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1942[k] = f_10 * lsi0_1211[k]
                    - f_11 * lsi1_1211[k]
                    + f_3 * pc_x[k] * lsk_1555[k];

        t_1943[k] = f_10 * lsi0_1212[k]
                    - f_11 * lsi1_1212[k]
                    + f_3 * pc_x[k] * lsk_1556[k];

        t_1944[k] = pa_y[k] * ksl0_1584[k]
                    - f_14 * pc_y[k] * ksl1_1584[k];
    }

#pragma omp simd aligned(t_1945, t_1946, t_1947, pc_x, lsi0_1214, lsi0_1215, lsi0_1216, \
                         lsi1_1214, lsi1_1215, lsi1_1216, lsk_1558, lsk_1559, \
                         lsk_1560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1945[k] = f_8 * lsi0_1214[k]
                    - f_9 * lsi1_1214[k]
                    + f_3 * pc_x[k] * lsk_1558[k];

        t_1946[k] = f_8 * lsi0_1215[k]
                    - f_9 * lsi1_1215[k]
                    + f_3 * pc_x[k] * lsk_1559[k];

        t_1947[k] = f_8 * lsi0_1216[k]
                    - f_9 * lsi1_1216[k]
                    + f_3 * pc_x[k] * lsk_1560[k];
    }

#pragma omp simd aligned(t_1948, t_1949, t_1950, pa_y, pc_x, pc_y, ksl0_1589, ksl1_1589, \
                         lsi0_1217, lsi0_1219, lsi1_1217, lsi1_1219, lsk_1561, \
                         lsk_1563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1948[k] = f_8 * lsi0_1217[k]
                    - f_9 * lsi1_1217[k]
                    + f_3 * pc_x[k] * lsk_1561[k];

        t_1949[k] = pa_y[k] * ksl0_1589[k]
                    - f_14 * pc_y[k] * ksl1_1589[k];

        t_1950[k] = f_6 * lsi0_1219[k]
                    - f_7 * lsi1_1219[k]
                    + f_3 * pc_x[k] * lsk_1563[k];
    }
}

static auto
compute_prim_lsl_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksl0,
                                                           const size_t ksk, const size_t ksl1,
                                                           const size_t lsi0, const size_t lsi1,
                                                           const size_t lsk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 3.5 / q;
    const auto f_22 = 3.0 / gamma;
    const auto f_23 = 3.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksl0_1595 = buffer.data(ksl0 + 1595);
    const auto *ksl0_1602 = buffer.data(ksl0 + 1602);
    const auto *ksl0_1611 = buffer.data(ksl0 + 1611);
    const auto *ksl0_1613 = buffer.data(ksl0 + 1613);
    const auto *ksl0_1614 = buffer.data(ksl0 + 1614);
    const auto *ksl0_1615 = buffer.data(ksl0 + 1615);
    const auto *ksl0_1616 = buffer.data(ksl0 + 1616);
    const auto *ksl0_1617 = buffer.data(ksl0 + 1617);
    const auto *ksl0_1619 = buffer.data(ksl0 + 1619);

    const auto *ksk_1252 = buffer.data(ksk + 1252);
    const auto *ksk_1288 = buffer.data(ksk + 1288);
    const auto *ksk_1290 = buffer.data(ksk + 1290);
    const auto *ksk_1291 = buffer.data(ksk + 1291);
    const auto *ksk_1292 = buffer.data(ksk + 1292);
    const auto *ksk_1293 = buffer.data(ksk + 1293);
    const auto *ksk_1294 = buffer.data(ksk + 1294);
    const auto *ksk_1295 = buffer.data(ksk + 1295);

    const auto *ksl1_1595 = buffer.data(ksl1 + 1595);
    const auto *ksl1_1602 = buffer.data(ksl1 + 1602);
    const auto *ksl1_1611 = buffer.data(ksl1 + 1611);
    const auto *ksl1_1613 = buffer.data(ksl1 + 1613);
    const auto *ksl1_1614 = buffer.data(ksl1 + 1614);
    const auto *ksl1_1615 = buffer.data(ksl1 + 1615);
    const auto *ksl1_1616 = buffer.data(ksl1 + 1616);
    const auto *ksl1_1617 = buffer.data(ksl1 + 1617);
    const auto *ksl1_1619 = buffer.data(ksl1 + 1619);

    const auto *lsi0_1220 = buffer.data(lsi0 + 1220);
    const auto *lsi0_1221 = buffer.data(lsi0 + 1221);
    const auto *lsi0_1222 = buffer.data(lsi0 + 1222);
    const auto *lsi0_1223 = buffer.data(lsi0 + 1223);
    const auto *lsi0_1225 = buffer.data(lsi0 + 1225);
    const auto *lsi0_1226 = buffer.data(lsi0 + 1226);
    const auto *lsi0_1227 = buffer.data(lsi0 + 1227);
    const auto *lsi0_1228 = buffer.data(lsi0 + 1228);
    const auto *lsi0_1229 = buffer.data(lsi0 + 1229);
    const auto *lsi0_1230 = buffer.data(lsi0 + 1230);
    const auto *lsi0_1232 = buffer.data(lsi0 + 1232);
    const auto *lsi0_1234 = buffer.data(lsi0 + 1234);
    const auto *lsi0_1235 = buffer.data(lsi0 + 1235);
    const auto *lsi0_1237 = buffer.data(lsi0 + 1237);
    const auto *lsi0_1238 = buffer.data(lsi0 + 1238);
    const auto *lsi0_1239 = buffer.data(lsi0 + 1239);
    const auto *lsi0_1241 = buffer.data(lsi0 + 1241);
    const auto *lsi0_1242 = buffer.data(lsi0 + 1242);
    const auto *lsi0_1243 = buffer.data(lsi0 + 1243);
    const auto *lsi0_1244 = buffer.data(lsi0 + 1244);
    const auto *lsi0_1246 = buffer.data(lsi0 + 1246);
    const auto *lsi0_1247 = buffer.data(lsi0 + 1247);
    const auto *lsi0_1248 = buffer.data(lsi0 + 1248);
    const auto *lsi0_1249 = buffer.data(lsi0 + 1249);
    const auto *lsi0_1250 = buffer.data(lsi0 + 1250);
    const auto *lsi0_1252 = buffer.data(lsi0 + 1252);
    const auto *lsi0_1253 = buffer.data(lsi0 + 1253);
    const auto *lsi0_1254 = buffer.data(lsi0 + 1254);
    const auto *lsi0_1255 = buffer.data(lsi0 + 1255);
    const auto *lsi0_1256 = buffer.data(lsi0 + 1256);
    const auto *lsi0_1257 = buffer.data(lsi0 + 1257);
    const auto *lsi0_1258 = buffer.data(lsi0 + 1258);
    const auto *lsi0_1259 = buffer.data(lsi0 + 1259);

    const auto *lsi1_1220 = buffer.data(lsi1 + 1220);
    const auto *lsi1_1221 = buffer.data(lsi1 + 1221);
    const auto *lsi1_1222 = buffer.data(lsi1 + 1222);
    const auto *lsi1_1223 = buffer.data(lsi1 + 1223);
    const auto *lsi1_1225 = buffer.data(lsi1 + 1225);
    const auto *lsi1_1226 = buffer.data(lsi1 + 1226);
    const auto *lsi1_1227 = buffer.data(lsi1 + 1227);
    const auto *lsi1_1228 = buffer.data(lsi1 + 1228);
    const auto *lsi1_1229 = buffer.data(lsi1 + 1229);
    const auto *lsi1_1230 = buffer.data(lsi1 + 1230);
    const auto *lsi1_1232 = buffer.data(lsi1 + 1232);
    const auto *lsi1_1234 = buffer.data(lsi1 + 1234);
    const auto *lsi1_1235 = buffer.data(lsi1 + 1235);
    const auto *lsi1_1237 = buffer.data(lsi1 + 1237);
    const auto *lsi1_1238 = buffer.data(lsi1 + 1238);
    const auto *lsi1_1239 = buffer.data(lsi1 + 1239);
    const auto *lsi1_1241 = buffer.data(lsi1 + 1241);
    const auto *lsi1_1242 = buffer.data(lsi1 + 1242);
    const auto *lsi1_1243 = buffer.data(lsi1 + 1243);
    const auto *lsi1_1244 = buffer.data(lsi1 + 1244);
    const auto *lsi1_1246 = buffer.data(lsi1 + 1246);
    const auto *lsi1_1247 = buffer.data(lsi1 + 1247);
    const auto *lsi1_1248 = buffer.data(lsi1 + 1248);
    const auto *lsi1_1249 = buffer.data(lsi1 + 1249);
    const auto *lsi1_1250 = buffer.data(lsi1 + 1250);
    const auto *lsi1_1252 = buffer.data(lsi1 + 1252);
    const auto *lsi1_1253 = buffer.data(lsi1 + 1253);
    const auto *lsi1_1254 = buffer.data(lsi1 + 1254);
    const auto *lsi1_1255 = buffer.data(lsi1 + 1255);
    const auto *lsi1_1256 = buffer.data(lsi1 + 1256);
    const auto *lsi1_1257 = buffer.data(lsi1 + 1257);
    const auto *lsi1_1258 = buffer.data(lsi1 + 1258);
    const auto *lsi1_1259 = buffer.data(lsi1 + 1259);

    const auto *lsk_1564 = buffer.data(lsk + 1564);
    const auto *lsk_1565 = buffer.data(lsk + 1565);
    const auto *lsk_1566 = buffer.data(lsk + 1566);
    const auto *lsk_1567 = buffer.data(lsk + 1567);
    const auto *lsk_1569 = buffer.data(lsk + 1569);
    const auto *lsk_1570 = buffer.data(lsk + 1570);
    const auto *lsk_1571 = buffer.data(lsk + 1571);
    const auto *lsk_1572 = buffer.data(lsk + 1572);
    const auto *lsk_1573 = buffer.data(lsk + 1573);
    const auto *lsk_1574 = buffer.data(lsk + 1574);
    const auto *lsk_1576 = buffer.data(lsk + 1576);
    const auto *lsk_1577 = buffer.data(lsk + 1577);
    const auto *lsk_1578 = buffer.data(lsk + 1578);
    const auto *lsk_1579 = buffer.data(lsk + 1579);
    const auto *lsk_1580 = buffer.data(lsk + 1580);
    const auto *lsk_1581 = buffer.data(lsk + 1581);
    const auto *lsk_1582 = buffer.data(lsk + 1582);
    const auto *lsk_1583 = buffer.data(lsk + 1583);
    const auto *lsk_1584 = buffer.data(lsk + 1584);
    const auto *lsk_1586 = buffer.data(lsk + 1586);
    const auto *lsk_1587 = buffer.data(lsk + 1587);
    const auto *lsk_1589 = buffer.data(lsk + 1589);
    const auto *lsk_1590 = buffer.data(lsk + 1590);
    const auto *lsk_1591 = buffer.data(lsk + 1591);
    const auto *lsk_1593 = buffer.data(lsk + 1593);
    const auto *lsk_1594 = buffer.data(lsk + 1594);
    const auto *lsk_1595 = buffer.data(lsk + 1595);
    const auto *lsk_1596 = buffer.data(lsk + 1596);
    const auto *lsk_1598 = buffer.data(lsk + 1598);
    const auto *lsk_1599 = buffer.data(lsk + 1599);
    const auto *lsk_1600 = buffer.data(lsk + 1600);
    const auto *lsk_1601 = buffer.data(lsk + 1601);
    const auto *lsk_1602 = buffer.data(lsk + 1602);
    const auto *lsk_1604 = buffer.data(lsk + 1604);
    const auto *lsk_1605 = buffer.data(lsk + 1605);
    const auto *lsk_1606 = buffer.data(lsk + 1606);
    const auto *lsk_1607 = buffer.data(lsk + 1607);
    const auto *lsk_1608 = buffer.data(lsk + 1608);
    const auto *lsk_1609 = buffer.data(lsk + 1609);
    const auto *lsk_1611 = buffer.data(lsk + 1611);
    const auto *lsk_1612 = buffer.data(lsk + 1612);
    const auto *lsk_1613 = buffer.data(lsk + 1613);
    const auto *lsk_1614 = buffer.data(lsk + 1614);
    const auto *lsk_1615 = buffer.data(lsk + 1615);
    const auto *lsk_1616 = buffer.data(lsk + 1616);
    const auto *lsk_1617 = buffer.data(lsk + 1617);
    const auto *lsk_1618 = buffer.data(lsk + 1618);
    const auto *lsk_1619 = buffer.data(lsk + 1619);

#pragma omp simd aligned(t_1951, t_1952, t_1953, pc_x, lsi0_1220, lsi0_1221, lsi0_1222, \
                         lsi1_1220, lsi1_1221, lsi1_1222, lsk_1564, lsk_1565, \
                         lsk_1566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1951[k] = f_6 * lsi0_1220[k]
                    - f_7 * lsi1_1220[k]
                    + f_3 * pc_x[k] * lsk_1564[k];

        t_1952[k] = f_6 * lsi0_1221[k]
                    - f_7 * lsi1_1221[k]
                    + f_3 * pc_x[k] * lsk_1565[k];

        t_1953[k] = f_6 * lsi0_1222[k]
                    - f_7 * lsi1_1222[k]
                    + f_3 * pc_x[k] * lsk_1566[k];
    }

#pragma omp simd aligned(t_1954, t_1955, t_1956, pa_y, pc_x, pc_y, ksl0_1595, ksl1_1595, \
                         lsi0_1223, lsi0_1225, lsi1_1223, lsi1_1225, lsk_1567, \
                         lsk_1569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1954[k] = f_6 * lsi0_1223[k]
                    - f_7 * lsi1_1223[k]
                    + f_3 * pc_x[k] * lsk_1567[k];

        t_1955[k] = pa_y[k] * ksl0_1595[k]
                    - f_14 * pc_y[k] * ksl1_1595[k];

        t_1956[k] = f_4 * lsi0_1225[k]
                    - f_5 * lsi1_1225[k]
                    + f_3 * pc_x[k] * lsk_1569[k];
    }

#pragma omp simd aligned(t_1957, t_1958, t_1959, pc_x, lsi0_1226, lsi0_1227, lsi0_1228, \
                         lsi1_1226, lsi1_1227, lsi1_1228, lsk_1570, lsk_1571, \
                         lsk_1572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1957[k] = f_4 * lsi0_1226[k]
                    - f_5 * lsi1_1226[k]
                    + f_3 * pc_x[k] * lsk_1570[k];

        t_1958[k] = f_4 * lsi0_1227[k]
                    - f_5 * lsi1_1227[k]
                    + f_3 * pc_x[k] * lsk_1571[k];

        t_1959[k] = f_4 * lsi0_1228[k]
                    - f_5 * lsi1_1228[k]
                    + f_3 * pc_x[k] * lsk_1572[k];
    }

#pragma omp simd aligned(t_1960, t_1961, t_1962, t_1963, pa_y, pc_x, pc_y, ksl0_1602, \
                         ksl1_1602, lsi0_1229, lsi0_1230, lsi1_1229, lsi1_1230, lsk_1573, \
                         lsk_1574, lsk_1576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1960[k] = f_4 * lsi0_1229[k]
                    - f_5 * lsi1_1229[k]
                    + f_3 * pc_x[k] * lsk_1573[k];

        t_1961[k] = f_4 * lsi0_1230[k]
                    - f_5 * lsi1_1230[k]
                    + f_3 * pc_x[k] * lsk_1574[k];

        t_1962[k] = pa_y[k] * ksl0_1602[k]
                    - f_14 * pc_y[k] * ksl1_1602[k];

        t_1963[k] = f_3 * pc_x[k] * lsk_1576[k];
    }

#pragma omp simd aligned(t_1964, t_1965, t_1966, t_1967, t_1968, t_1969, t_1970, pc_x, \
                         lsk_1577, lsk_1578, lsk_1579, lsk_1580, lsk_1581, lsk_1582, \
                         lsk_1583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1964[k] = f_3 * pc_x[k] * lsk_1577[k];

        t_1965[k] = f_3 * pc_x[k] * lsk_1578[k];

        t_1966[k] = f_3 * pc_x[k] * lsk_1579[k];

        t_1967[k] = f_3 * pc_x[k] * lsk_1580[k];

        t_1968[k] = f_3 * pc_x[k] * lsk_1581[k];

        t_1969[k] = f_3 * pc_x[k] * lsk_1582[k];

        t_1970[k] = f_3 * pc_x[k] * lsk_1583[k];
    }

#pragma omp simd aligned(t_1971, t_1972, t_1973, pa_y, pc_y, pc_z, ksl0_1611, ksl0_1613, \
                         ksk_1252, ksk_1288, ksk_1290, ksl1_1611, ksl1_1613, \
                         lsk_1576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1971[k] = pa_y[k] * ksl0_1611[k]
                    + f_0 * ksk_1288[k]
                    - f_14 * pc_y[k] * ksl1_1611[k];

        t_1972[k] = f_21 * ksk_1252[k]
                    + f_3 * pc_z[k] * lsk_1576[k];

        t_1973[k] = pa_y[k] * ksl0_1613[k]
                    + f_20 * ksk_1290[k]
                    - f_14 * pc_y[k] * ksl1_1613[k];
    }

#pragma omp simd aligned(t_1974, t_1975, t_1976, pa_y, pc_y, ksl0_1614, ksl0_1615, ksl0_1616, \
                         ksk_1291, ksk_1292, ksk_1293, ksl1_1614, ksl1_1615, \
                         ksl1_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1974[k] = pa_y[k] * ksl0_1614[k]
                    + f_19 * ksk_1291[k]
                    - f_14 * pc_y[k] * ksl1_1614[k];

        t_1975[k] = pa_y[k] * ksl0_1615[k]
                    + f_18 * ksk_1292[k]
                    - f_14 * pc_y[k] * ksl1_1615[k];

        t_1976[k] = pa_y[k] * ksl0_1616[k]
                    + f_17 * ksk_1293[k]
                    - f_14 * pc_y[k] * ksl1_1616[k];
    }

#pragma omp simd aligned(t_1977, t_1978, t_1979, pa_y, pc_y, ksl0_1617, ksl0_1619, ksk_1294, \
                         ksk_1295, ksl1_1617, ksl1_1619, lsk_1583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1977[k] = pa_y[k] * ksl0_1617[k]
                    + f_16 * ksk_1294[k]
                    - f_14 * pc_y[k] * ksl1_1617[k];

        t_1978[k] = f_15 * ksk_1295[k]
                    + f_3 * pc_y[k] * lsk_1583[k];

        t_1979[k] = pa_y[k] * ksl0_1619[k]
                    - f_14 * pc_y[k] * ksl1_1619[k];
    }

#pragma omp simd aligned(t_1980, t_1981, t_1982, t_1983, t_1984, pc_x, pc_y, lsi0_1232, \
                         lsi0_1234, lsi0_1235, lsi1_1232, lsi1_1234, lsi1_1235, lsk_1584, \
                         lsk_1586, lsk_1587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1980[k] = f_1 * lsi0_1232[k]
                    - f_2 * lsi1_1232[k]
                    + f_3 * pc_x[k] * lsk_1584[k];

        t_1981[k] = f_3 * pc_y[k] * lsk_1584[k];

        t_1982[k] = f_22 * lsi0_1234[k]
                    - f_23 * lsi1_1234[k]
                    + f_3 * pc_x[k] * lsk_1586[k];

        t_1983[k] = f_12 * lsi0_1235[k]
                    - f_13 * lsi1_1235[k]
                    + f_3 * pc_x[k] * lsk_1587[k];

        t_1984[k] = f_3 * pc_y[k] * lsk_1586[k];
    }

#pragma omp simd aligned(t_1985, t_1986, t_1987, t_1988, pc_x, pc_y, lsi0_1237, lsi0_1238, \
                         lsi0_1239, lsi1_1237, lsi1_1238, lsi1_1239, lsk_1589, lsk_1590, \
                         lsk_1591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1985[k] = f_12 * lsi0_1237[k]
                    - f_13 * lsi1_1237[k]
                    + f_3 * pc_x[k] * lsk_1589[k];

        t_1986[k] = f_10 * lsi0_1238[k]
                    - f_11 * lsi1_1238[k]
                    + f_3 * pc_x[k] * lsk_1590[k];

        t_1987[k] = f_10 * lsi0_1239[k]
                    - f_11 * lsi1_1239[k]
                    + f_3 * pc_x[k] * lsk_1591[k];

        t_1988[k] = f_3 * pc_y[k] * lsk_1589[k];
    }

#pragma omp simd aligned(t_1989, t_1990, t_1991, pc_x, lsi0_1241, lsi0_1242, lsi0_1243, \
                         lsi1_1241, lsi1_1242, lsi1_1243, lsk_1593, lsk_1594, \
                         lsk_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1989[k] = f_10 * lsi0_1241[k]
                    - f_11 * lsi1_1241[k]
                    + f_3 * pc_x[k] * lsk_1593[k];

        t_1990[k] = f_8 * lsi0_1242[k]
                    - f_9 * lsi1_1242[k]
                    + f_3 * pc_x[k] * lsk_1594[k];

        t_1991[k] = f_8 * lsi0_1243[k]
                    - f_9 * lsi1_1243[k]
                    + f_3 * pc_x[k] * lsk_1595[k];
    }

#pragma omp simd aligned(t_1992, t_1993, t_1994, t_1995, pc_x, pc_y, lsi0_1244, lsi0_1246, \
                         lsi0_1247, lsi1_1244, lsi1_1246, lsi1_1247, lsk_1593, lsk_1596, \
                         lsk_1598, lsk_1599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1992[k] = f_8 * lsi0_1244[k]
                    - f_9 * lsi1_1244[k]
                    + f_3 * pc_x[k] * lsk_1596[k];

        t_1993[k] = f_3 * pc_y[k] * lsk_1593[k];

        t_1994[k] = f_8 * lsi0_1246[k]
                    - f_9 * lsi1_1246[k]
                    + f_3 * pc_x[k] * lsk_1598[k];

        t_1995[k] = f_6 * lsi0_1247[k]
                    - f_7 * lsi1_1247[k]
                    + f_3 * pc_x[k] * lsk_1599[k];
    }

#pragma omp simd aligned(t_1996, t_1997, t_1998, t_1999, pc_x, pc_y, lsi0_1248, lsi0_1249, \
                         lsi0_1250, lsi1_1248, lsi1_1249, lsi1_1250, lsk_1598, lsk_1600, \
                         lsk_1601, lsk_1602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1996[k] = f_6 * lsi0_1248[k]
                    - f_7 * lsi1_1248[k]
                    + f_3 * pc_x[k] * lsk_1600[k];

        t_1997[k] = f_6 * lsi0_1249[k]
                    - f_7 * lsi1_1249[k]
                    + f_3 * pc_x[k] * lsk_1601[k];

        t_1998[k] = f_6 * lsi0_1250[k]
                    - f_7 * lsi1_1250[k]
                    + f_3 * pc_x[k] * lsk_1602[k];

        t_1999[k] = f_3 * pc_y[k] * lsk_1598[k];
    }

#pragma omp simd aligned(t_2000, t_2001, t_2002, pc_x, lsi0_1252, lsi0_1253, lsi0_1254, \
                         lsi1_1252, lsi1_1253, lsi1_1254, lsk_1604, lsk_1605, \
                         lsk_1606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2000[k] = f_6 * lsi0_1252[k]
                    - f_7 * lsi1_1252[k]
                    + f_3 * pc_x[k] * lsk_1604[k];

        t_2001[k] = f_4 * lsi0_1253[k]
                    - f_5 * lsi1_1253[k]
                    + f_3 * pc_x[k] * lsk_1605[k];

        t_2002[k] = f_4 * lsi0_1254[k]
                    - f_5 * lsi1_1254[k]
                    + f_3 * pc_x[k] * lsk_1606[k];
    }

#pragma omp simd aligned(t_2003, t_2004, t_2005, t_2006, pc_x, pc_y, lsi0_1255, lsi0_1256, \
                         lsi0_1257, lsi1_1255, lsi1_1256, lsi1_1257, lsk_1604, lsk_1607, \
                         lsk_1608, lsk_1609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2003[k] = f_4 * lsi0_1255[k]
                    - f_5 * lsi1_1255[k]
                    + f_3 * pc_x[k] * lsk_1607[k];

        t_2004[k] = f_4 * lsi0_1256[k]
                    - f_5 * lsi1_1256[k]
                    + f_3 * pc_x[k] * lsk_1608[k];

        t_2005[k] = f_4 * lsi0_1257[k]
                    - f_5 * lsi1_1257[k]
                    + f_3 * pc_x[k] * lsk_1609[k];

        t_2006[k] = f_3 * pc_y[k] * lsk_1604[k];
    }

#pragma omp simd aligned(t_2007, t_2008, t_2009, t_2010, t_2011, t_2012, pc_x, lsi0_1259, \
                         lsi1_1259, lsk_1611, lsk_1612, lsk_1613, lsk_1614, lsk_1615, \
                         lsk_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2007[k] = f_4 * lsi0_1259[k]
                    - f_5 * lsi1_1259[k]
                    + f_3 * pc_x[k] * lsk_1611[k];

        t_2008[k] = f_3 * pc_x[k] * lsk_1612[k];

        t_2009[k] = f_3 * pc_x[k] * lsk_1613[k];

        t_2010[k] = f_3 * pc_x[k] * lsk_1614[k];

        t_2011[k] = f_3 * pc_x[k] * lsk_1615[k];

        t_2012[k] = f_3 * pc_x[k] * lsk_1616[k];
    }

#pragma omp simd aligned(t_2013, t_2014, t_2015, t_2016, t_2017, pc_x, pc_y, lsi0_1253, \
                         lsi0_1254, lsi1_1253, lsi1_1254, lsk_1612, lsk_1613, lsk_1617, \
                         lsk_1618, lsk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2013[k] = f_3 * pc_x[k] * lsk_1617[k];

        t_2014[k] = f_3 * pc_x[k] * lsk_1618[k];

        t_2015[k] = f_3 * pc_x[k] * lsk_1619[k];

        t_2016[k] = f_1 * lsi0_1253[k]
                    - f_2 * lsi1_1253[k]
                    + f_3 * pc_y[k] * lsk_1612[k];

        t_2017[k] = f_22 * lsi0_1254[k]
                    - f_23 * lsi1_1254[k]
                    + f_3 * pc_y[k] * lsk_1613[k];
    }

#pragma omp simd aligned(t_2018, t_2019, t_2020, pc_y, lsi0_1255, lsi0_1256, lsi0_1257, \
                         lsi1_1255, lsi1_1256, lsi1_1257, lsk_1614, lsk_1615, \
                         lsk_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2018[k] = f_12 * lsi0_1255[k]
                    - f_13 * lsi1_1255[k]
                    + f_3 * pc_y[k] * lsk_1614[k];

        t_2019[k] = f_10 * lsi0_1256[k]
                    - f_11 * lsi1_1256[k]
                    + f_3 * pc_y[k] * lsk_1615[k];

        t_2020[k] = f_8 * lsi0_1257[k]
                    - f_9 * lsi1_1257[k]
                    + f_3 * pc_y[k] * lsk_1616[k];
    }

#pragma omp simd aligned(t_2021, t_2022, t_2023, t_2024, pc_y, pc_z, ksk_1295, lsi0_1258, \
                         lsi0_1259, lsi1_1258, lsi1_1259, lsk_1617, lsk_1618, \
                         lsk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2021[k] = f_6 * lsi0_1258[k]
                    - f_7 * lsi1_1258[k]
                    + f_3 * pc_y[k] * lsk_1617[k];

        t_2022[k] = f_4 * lsi0_1259[k]
                    - f_5 * lsi1_1259[k]
                    + f_3 * pc_y[k] * lsk_1618[k];

        t_2023[k] = f_3 * pc_y[k] * lsk_1619[k];

        t_2024[k] = f_0 * ksk_1295[k]
                    + f_1 * lsi0_1259[k]
                    - f_2 * lsi1_1259[k]
                    + f_3 * pc_z[k] * lsk_1619[k];
    }
}

auto
compute_prim_lsl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ksl0, const size_t ksk,
                                                   const size_t ksl1, const size_t lsi0,
                                                   const size_t lsi1, const size_t lsk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_lsl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, ksl0, ksk,
                                                              ksl1, lsi0, lsi1, lsk, ncols,
                                                              gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece9(buffer, target, pc, ksk, lsi0,
                                                              lsi1, lsk, ncols, gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, ksl0,
                                                               ksk, ksl1, lsi0, lsi1, lsk,
                                                               ncols, gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, ksl0,
                                                               ksk, ksl1, lsi0, lsi1, lsk,
                                                               ncols, gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece12(buffer, target, pa, pc, ksl0,
                                                               ksk, ksl1, lsk, ncols, gamma, p,
                                                               q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, ksl0,
                                                               ksk, ksl1, lsi0, lsi1, lsk,
                                                               ncols, gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece14(buffer, target, pa, pc, ksl0,
                                                               ksk, ksl1, lsi0, lsi1, lsk,
                                                               ncols, gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece15(buffer, target, pc, ksk, lsi0,
                                                               lsi1, lsk, ncols, gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece16(buffer, target, pa, pc, ksl0,
                                                               ksk, ksl1, lsi0, lsi1, lsk,
                                                               ncols, gamma, p, q);

    compute_prim_lsl_three_center_electron_repulsion_0_piece17(buffer, target, pa, pc, ksl0,
                                                               ksk, ksl1, lsi0, lsi1, lsk,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
