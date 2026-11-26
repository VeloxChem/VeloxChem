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


#include "SimdThreeCenterElectronRepulsionVrrRecHSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_hsl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsl0,
                                                          const size_t gsk, const size_t gsl1,
                                                          const size_t hsi0, const size_t hsi1,
                                                          const size_t hsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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

    const auto *gsl0_0 = buffer.data(gsl0 + 0);
    const auto *gsl0_3 = buffer.data(gsl0 + 3);
    const auto *gsl0_5 = buffer.data(gsl0 + 5);
    const auto *gsl0_6 = buffer.data(gsl0 + 6);
    const auto *gsl0_9 = buffer.data(gsl0 + 9);
    const auto *gsl0_10 = buffer.data(gsl0 + 10);
    const auto *gsl0_14 = buffer.data(gsl0 + 14);
    const auto *gsl0_15 = buffer.data(gsl0 + 15);
    const auto *gsl0_20 = buffer.data(gsl0 + 20);
    const auto *gsl0_21 = buffer.data(gsl0 + 21);
    const auto *gsl0_27 = buffer.data(gsl0 + 27);
    const auto *gsl0_36 = buffer.data(gsl0 + 36);
    const auto *gsl0_44 = buffer.data(gsl0 + 44);

    const auto *gsk_0 = buffer.data(gsk + 0);
    const auto *gsk_1 = buffer.data(gsk + 1);
    const auto *gsk_2 = buffer.data(gsk + 2);
    const auto *gsk_3 = buffer.data(gsk + 3);
    const auto *gsk_5 = buffer.data(gsk + 5);
    const auto *gsk_6 = buffer.data(gsk + 6);
    const auto *gsk_9 = buffer.data(gsk + 9);
    const auto *gsk_10 = buffer.data(gsk + 10);
    const auto *gsk_14 = buffer.data(gsk + 14);
    const auto *gsk_15 = buffer.data(gsk + 15);
    const auto *gsk_20 = buffer.data(gsk + 20);
    const auto *gsk_28 = buffer.data(gsk + 28);
    const auto *gsk_30 = buffer.data(gsk + 30);
    const auto *gsk_31 = buffer.data(gsk + 31);
    const auto *gsk_32 = buffer.data(gsk + 32);
    const auto *gsk_33 = buffer.data(gsk + 33);
    const auto *gsk_35 = buffer.data(gsk + 35);
    const auto *gsk_64 = buffer.data(gsk + 64);
    const auto *gsk_66 = buffer.data(gsk + 66);
    const auto *gsk_67 = buffer.data(gsk + 67);
    const auto *gsk_68 = buffer.data(gsk + 68);
    const auto *gsk_69 = buffer.data(gsk + 69);
    const auto *gsk_70 = buffer.data(gsk + 70);
    const auto *gsk_71 = buffer.data(gsk + 71);
    const auto *gsk_100 = buffer.data(gsk + 100);
    const auto *gsk_101 = buffer.data(gsk + 101);
    const auto *gsk_102 = buffer.data(gsk + 102);
    const auto *gsk_103 = buffer.data(gsk + 103);
    const auto *gsk_104 = buffer.data(gsk + 104);
    const auto *gsk_105 = buffer.data(gsk + 105);
    const auto *gsk_107 = buffer.data(gsk + 107);

    const auto *gsl1_0 = buffer.data(gsl1 + 0);
    const auto *gsl1_3 = buffer.data(gsl1 + 3);
    const auto *gsl1_5 = buffer.data(gsl1 + 5);
    const auto *gsl1_6 = buffer.data(gsl1 + 6);
    const auto *gsl1_9 = buffer.data(gsl1 + 9);
    const auto *gsl1_10 = buffer.data(gsl1 + 10);
    const auto *gsl1_14 = buffer.data(gsl1 + 14);
    const auto *gsl1_15 = buffer.data(gsl1 + 15);
    const auto *gsl1_20 = buffer.data(gsl1 + 20);
    const auto *gsl1_21 = buffer.data(gsl1 + 21);
    const auto *gsl1_27 = buffer.data(gsl1 + 27);
    const auto *gsl1_36 = buffer.data(gsl1 + 36);
    const auto *gsl1_44 = buffer.data(gsl1 + 44);

    const auto *hsi0_0 = buffer.data(hsi0 + 0);
    const auto *hsi0_1 = buffer.data(hsi0 + 1);
    const auto *hsi0_2 = buffer.data(hsi0 + 2);
    const auto *hsi0_3 = buffer.data(hsi0 + 3);
    const auto *hsi0_5 = buffer.data(hsi0 + 5);
    const auto *hsi0_6 = buffer.data(hsi0 + 6);
    const auto *hsi0_8 = buffer.data(hsi0 + 8);
    const auto *hsi0_9 = buffer.data(hsi0 + 9);
    const auto *hsi0_10 = buffer.data(hsi0 + 10);
    const auto *hsi0_12 = buffer.data(hsi0 + 12);
    const auto *hsi0_13 = buffer.data(hsi0 + 13);
    const auto *hsi0_14 = buffer.data(hsi0 + 14);
    const auto *hsi0_21 = buffer.data(hsi0 + 21);
    const auto *hsi0_23 = buffer.data(hsi0 + 23);
    const auto *hsi0_24 = buffer.data(hsi0 + 24);
    const auto *hsi0_25 = buffer.data(hsi0 + 25);
    const auto *hsi0_26 = buffer.data(hsi0 + 26);
    const auto *hsi0_27 = buffer.data(hsi0 + 27);
    const auto *hsi0_31 = buffer.data(hsi0 + 31);
    const auto *hsi0_34 = buffer.data(hsi0 + 34);
    const auto *hsi0_35 = buffer.data(hsi0 + 35);
    const auto *hsi0_38 = buffer.data(hsi0 + 38);
    const auto *hsi0_39 = buffer.data(hsi0 + 39);
    const auto *hsi0_40 = buffer.data(hsi0 + 40);
    const auto *hsi0_49 = buffer.data(hsi0 + 49);
    const auto *hsi0_50 = buffer.data(hsi0 + 50);
    const auto *hsi0_51 = buffer.data(hsi0 + 51);
    const auto *hsi0_52 = buffer.data(hsi0 + 52);
    const auto *hsi0_53 = buffer.data(hsi0 + 53);
    const auto *hsi0_58 = buffer.data(hsi0 + 58);
    const auto *hsi0_60 = buffer.data(hsi0 + 60);
    const auto *hsi0_61 = buffer.data(hsi0 + 61);
    const auto *hsi0_63 = buffer.data(hsi0 + 63);
    const auto *hsi0_64 = buffer.data(hsi0 + 64);
    const auto *hsi0_65 = buffer.data(hsi0 + 65);
    const auto *hsi0_67 = buffer.data(hsi0 + 67);
    const auto *hsi0_68 = buffer.data(hsi0 + 68);
    const auto *hsi0_69 = buffer.data(hsi0 + 69);
    const auto *hsi0_70 = buffer.data(hsi0 + 70);
    const auto *hsi0_78 = buffer.data(hsi0 + 78);
    const auto *hsi0_79 = buffer.data(hsi0 + 79);

    const auto *hsi1_0 = buffer.data(hsi1 + 0);
    const auto *hsi1_1 = buffer.data(hsi1 + 1);
    const auto *hsi1_2 = buffer.data(hsi1 + 2);
    const auto *hsi1_3 = buffer.data(hsi1 + 3);
    const auto *hsi1_5 = buffer.data(hsi1 + 5);
    const auto *hsi1_6 = buffer.data(hsi1 + 6);
    const auto *hsi1_8 = buffer.data(hsi1 + 8);
    const auto *hsi1_9 = buffer.data(hsi1 + 9);
    const auto *hsi1_10 = buffer.data(hsi1 + 10);
    const auto *hsi1_12 = buffer.data(hsi1 + 12);
    const auto *hsi1_13 = buffer.data(hsi1 + 13);
    const auto *hsi1_14 = buffer.data(hsi1 + 14);
    const auto *hsi1_21 = buffer.data(hsi1 + 21);
    const auto *hsi1_23 = buffer.data(hsi1 + 23);
    const auto *hsi1_24 = buffer.data(hsi1 + 24);
    const auto *hsi1_25 = buffer.data(hsi1 + 25);
    const auto *hsi1_26 = buffer.data(hsi1 + 26);
    const auto *hsi1_27 = buffer.data(hsi1 + 27);
    const auto *hsi1_31 = buffer.data(hsi1 + 31);
    const auto *hsi1_34 = buffer.data(hsi1 + 34);
    const auto *hsi1_35 = buffer.data(hsi1 + 35);
    const auto *hsi1_38 = buffer.data(hsi1 + 38);
    const auto *hsi1_39 = buffer.data(hsi1 + 39);
    const auto *hsi1_40 = buffer.data(hsi1 + 40);
    const auto *hsi1_49 = buffer.data(hsi1 + 49);
    const auto *hsi1_50 = buffer.data(hsi1 + 50);
    const auto *hsi1_51 = buffer.data(hsi1 + 51);
    const auto *hsi1_52 = buffer.data(hsi1 + 52);
    const auto *hsi1_53 = buffer.data(hsi1 + 53);
    const auto *hsi1_58 = buffer.data(hsi1 + 58);
    const auto *hsi1_60 = buffer.data(hsi1 + 60);
    const auto *hsi1_61 = buffer.data(hsi1 + 61);
    const auto *hsi1_63 = buffer.data(hsi1 + 63);
    const auto *hsi1_64 = buffer.data(hsi1 + 64);
    const auto *hsi1_65 = buffer.data(hsi1 + 65);
    const auto *hsi1_67 = buffer.data(hsi1 + 67);
    const auto *hsi1_68 = buffer.data(hsi1 + 68);
    const auto *hsi1_69 = buffer.data(hsi1 + 69);
    const auto *hsi1_70 = buffer.data(hsi1 + 70);
    const auto *hsi1_78 = buffer.data(hsi1 + 78);
    const auto *hsi1_79 = buffer.data(hsi1 + 79);

    const auto *hsk_0 = buffer.data(hsk + 0);
    const auto *hsk_1 = buffer.data(hsk + 1);
    const auto *hsk_2 = buffer.data(hsk + 2);
    const auto *hsk_3 = buffer.data(hsk + 3);
    const auto *hsk_5 = buffer.data(hsk + 5);
    const auto *hsk_6 = buffer.data(hsk + 6);
    const auto *hsk_8 = buffer.data(hsk + 8);
    const auto *hsk_9 = buffer.data(hsk + 9);
    const auto *hsk_10 = buffer.data(hsk + 10);
    const auto *hsk_12 = buffer.data(hsk + 12);
    const auto *hsk_13 = buffer.data(hsk + 13);
    const auto *hsk_14 = buffer.data(hsk + 14);
    const auto *hsk_15 = buffer.data(hsk + 15);
    const auto *hsk_17 = buffer.data(hsk + 17);
    const auto *hsk_18 = buffer.data(hsk + 18);
    const auto *hsk_19 = buffer.data(hsk + 19);
    const auto *hsk_20 = buffer.data(hsk + 20);
    const auto *hsk_21 = buffer.data(hsk + 21);
    const auto *hsk_27 = buffer.data(hsk + 27);
    const auto *hsk_28 = buffer.data(hsk + 28);
    const auto *hsk_30 = buffer.data(hsk + 30);
    const auto *hsk_31 = buffer.data(hsk + 31);
    const auto *hsk_32 = buffer.data(hsk + 32);
    const auto *hsk_33 = buffer.data(hsk + 33);
    const auto *hsk_34 = buffer.data(hsk + 34);
    const auto *hsk_35 = buffer.data(hsk + 35);
    const auto *hsk_36 = buffer.data(hsk + 36);
    const auto *hsk_37 = buffer.data(hsk + 37);
    const auto *hsk_39 = buffer.data(hsk + 39);
    const auto *hsk_41 = buffer.data(hsk + 41);
    const auto *hsk_42 = buffer.data(hsk + 42);
    const auto *hsk_43 = buffer.data(hsk + 43);
    const auto *hsk_45 = buffer.data(hsk + 45);
    const auto *hsk_46 = buffer.data(hsk + 46);
    const auto *hsk_47 = buffer.data(hsk + 47);
    const auto *hsk_48 = buffer.data(hsk + 48);
    const auto *hsk_50 = buffer.data(hsk + 50);
    const auto *hsk_51 = buffer.data(hsk + 51);
    const auto *hsk_52 = buffer.data(hsk + 52);
    const auto *hsk_53 = buffer.data(hsk + 53);
    const auto *hsk_54 = buffer.data(hsk + 54);
    const auto *hsk_56 = buffer.data(hsk + 56);
    const auto *hsk_57 = buffer.data(hsk + 57);
    const auto *hsk_64 = buffer.data(hsk + 64);
    const auto *hsk_65 = buffer.data(hsk + 65);
    const auto *hsk_66 = buffer.data(hsk + 66);
    const auto *hsk_67 = buffer.data(hsk + 67);
    const auto *hsk_68 = buffer.data(hsk + 68);
    const auto *hsk_69 = buffer.data(hsk + 69);
    const auto *hsk_70 = buffer.data(hsk + 70);
    const auto *hsk_71 = buffer.data(hsk + 71);
    const auto *hsk_72 = buffer.data(hsk + 72);
    const auto *hsk_74 = buffer.data(hsk + 74);
    const auto *hsk_76 = buffer.data(hsk + 76);
    const auto *hsk_77 = buffer.data(hsk + 77);
    const auto *hsk_79 = buffer.data(hsk + 79);
    const auto *hsk_80 = buffer.data(hsk + 80);
    const auto *hsk_81 = buffer.data(hsk + 81);
    const auto *hsk_83 = buffer.data(hsk + 83);
    const auto *hsk_84 = buffer.data(hsk + 84);
    const auto *hsk_85 = buffer.data(hsk + 85);
    const auto *hsk_86 = buffer.data(hsk + 86);
    const auto *hsk_88 = buffer.data(hsk + 88);
    const auto *hsk_89 = buffer.data(hsk + 89);
    const auto *hsk_90 = buffer.data(hsk + 90);
    const auto *hsk_91 = buffer.data(hsk + 91);
    const auto *hsk_92 = buffer.data(hsk + 92);
    const auto *hsk_99 = buffer.data(hsk + 99);
    const auto *hsk_100 = buffer.data(hsk + 100);
    const auto *hsk_101 = buffer.data(hsk + 101);
    const auto *hsk_102 = buffer.data(hsk + 102);
    const auto *hsk_103 = buffer.data(hsk + 103);
    const auto *hsk_104 = buffer.data(hsk + 104);
    const auto *hsk_105 = buffer.data(hsk + 105);
    const auto *hsk_107 = buffer.data(hsk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, gsk_0, hsi0_0, \
                         hsi1_0, hsk_0, hsk_1, hsk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gsk_0[k]
                 + f_1 * hsi0_0[k]
                 - f_2 * hsi1_0[k]
                 + f_3 * pc_x[k] * hsk_0[k];

        t_1[k] = f_3 * pc_y[k] * hsk_0[k];

        t_2[k] = f_3 * pc_z[k] * hsk_0[k];

        t_3[k] = f_4 * hsi0_0[k]
                 - f_5 * hsi1_0[k]
                 + f_3 * pc_y[k] * hsk_1[k];

        t_4[k] = f_3 * pc_y[k] * hsk_2[k];

        t_5[k] = f_4 * hsi0_0[k]
                 - f_5 * hsi1_0[k]
                 + f_3 * pc_z[k] * hsk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, hsi0_1, hsi0_2, hsi0_3, hsi1_1, \
                         hsi1_2, hsi1_3, hsk_3, hsk_5, hsk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hsi0_1[k]
                 - f_7 * hsi1_1[k]
                 + f_3 * pc_y[k] * hsk_3[k];

        t_7[k] = f_3 * pc_z[k] * hsk_3[k];

        t_8[k] = f_3 * pc_y[k] * hsk_5[k];

        t_9[k] = f_6 * hsi0_2[k]
                 - f_7 * hsi1_2[k]
                 + f_3 * pc_z[k] * hsk_5[k];

        t_10[k] = f_8 * hsi0_3[k]
                  - f_9 * hsi1_3[k]
                  + f_3 * pc_y[k] * hsk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, hsi0_5, hsi0_6, \
                         hsi1_5, hsi1_6, hsk_6, hsk_8, hsk_9, hsk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * hsk_6[k];

        t_12[k] = f_4 * hsi0_5[k]
                  - f_5 * hsi1_5[k]
                  + f_3 * pc_y[k] * hsk_8[k];

        t_13[k] = f_3 * pc_y[k] * hsk_9[k];

        t_14[k] = f_8 * hsi0_5[k]
                  - f_9 * hsi1_5[k]
                  + f_3 * pc_z[k] * hsk_9[k];

        t_15[k] = f_10 * hsi0_6[k]
                  - f_11 * hsi1_6[k]
                  + f_3 * pc_y[k] * hsk_10[k];

        t_16[k] = f_3 * pc_z[k] * hsk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, hsi0_8, hsi0_9, hsi1_8, hsi1_9, \
                         hsk_12, hsk_13, hsk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * hsi0_8[k]
                  - f_7 * hsi1_8[k]
                  + f_3 * pc_y[k] * hsk_12[k];

        t_18[k] = f_4 * hsi0_9[k]
                  - f_5 * hsi1_9[k]
                  + f_3 * pc_y[k] * hsk_13[k];

        t_19[k] = f_3 * pc_y[k] * hsk_14[k];

        t_20[k] = f_10 * hsi0_9[k]
                  - f_11 * hsi1_9[k]
                  + f_3 * pc_z[k] * hsk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, hsi0_10, hsi0_12, hsi0_13, \
                         hsi1_10, hsi1_12, hsi1_13, hsk_15, hsk_17, \
                         hsk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_12 * hsi0_10[k]
                  - f_13 * hsi1_10[k]
                  + f_3 * pc_y[k] * hsk_15[k];

        t_22[k] = f_3 * pc_z[k] * hsk_15[k];

        t_23[k] = f_8 * hsi0_12[k]
                  - f_9 * hsi1_12[k]
                  + f_3 * pc_y[k] * hsk_17[k];

        t_24[k] = f_6 * hsi0_13[k]
                  - f_7 * hsi1_13[k]
                  + f_3 * pc_y[k] * hsk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, gsk_28, hsi0_14, \
                         hsi1_14, hsk_19, hsk_20, hsk_21, hsk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * hsi0_14[k]
                  - f_5 * hsi1_14[k]
                  + f_3 * pc_y[k] * hsk_19[k];

        t_26[k] = f_3 * pc_y[k] * hsk_20[k];

        t_27[k] = f_12 * hsi0_14[k]
                  - f_13 * hsi1_14[k]
                  + f_3 * pc_z[k] * hsk_20[k];

        t_28[k] = f_0 * gsk_28[k]
                  + f_3 * pc_x[k] * hsk_28[k];

        t_29[k] = f_3 * pc_z[k] * hsk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, gsk_30, gsk_31, gsk_32, \
                         gsk_33, hsk_27, hsk_30, hsk_31, hsk_32, \
                         hsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * gsk_30[k]
                  + f_3 * pc_x[k] * hsk_30[k];

        t_31[k] = f_0 * gsk_31[k]
                  + f_3 * pc_x[k] * hsk_31[k];

        t_32[k] = f_0 * gsk_32[k]
                  + f_3 * pc_x[k] * hsk_32[k];

        t_33[k] = f_0 * gsk_33[k]
                  + f_3 * pc_x[k] * hsk_33[k];

        t_34[k] = f_3 * pc_y[k] * hsk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, gsk_35, hsi0_21, hsi0_23, \
                         hsi1_21, hsi1_23, hsk_28, hsk_30, hsk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * gsk_35[k]
                  + f_3 * pc_x[k] * hsk_35[k];

        t_36[k] = f_1 * hsi0_21[k]
                  - f_2 * hsi1_21[k]
                  + f_3 * pc_y[k] * hsk_28[k];

        t_37[k] = f_3 * pc_z[k] * hsk_28[k];

        t_38[k] = f_12 * hsi0_23[k]
                  - f_13 * hsi1_23[k]
                  + f_3 * pc_y[k] * hsk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pc_y, hsi0_24, hsi0_25, hsi0_26, hsi1_24, hsi1_25, \
                         hsi1_26, hsk_31, hsk_32, hsk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * hsi0_24[k]
                  - f_11 * hsi1_24[k]
                  + f_3 * pc_y[k] * hsk_31[k];

        t_40[k] = f_8 * hsi0_25[k]
                  - f_9 * hsi1_25[k]
                  + f_3 * pc_y[k] * hsk_32[k];

        t_41[k] = f_6 * hsi0_26[k]
                  - f_7 * hsi1_26[k]
                  + f_3 * pc_y[k] * hsk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_y, pc_z, gsl0_0, gsk_0, \
                         gsl1_0, hsi0_27, hsi1_27, hsk_34, hsk_35, \
                         hsk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * hsi0_27[k]
                  - f_5 * hsi1_27[k]
                  + f_3 * pc_y[k] * hsk_34[k];

        t_43[k] = f_3 * pc_y[k] * hsk_35[k];

        t_44[k] = f_1 * hsi0_27[k]
                  - f_2 * hsi1_27[k]
                  + f_3 * pc_z[k] * hsk_35[k];

        t_45[k] = pa_y[k] * gsl0_0[k]
                  - f_14 * pc_y[k] * gsl1_0[k];

        t_46[k] = f_15 * gsk_0[k]
                  + f_3 * pc_y[k] * hsk_36[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_y, pc_z, gsl0_3, gsl0_5, gsk_1, \
                         gsl1_3, gsl1_5, hsk_36, hsk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * hsk_36[k];

        t_48[k] = pa_y[k] * gsl0_3[k]
                  + f_16 * gsk_1[k]
                  - f_14 * pc_y[k] * gsl1_3[k];

        t_49[k] = f_3 * pc_z[k] * hsk_37[k];

        t_50[k] = pa_y[k] * gsl0_5[k]
                  - f_14 * pc_y[k] * gsl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pc_y, pc_z, gsl0_6, gsl0_9, gsk_3, \
                         gsk_5, gsl1_6, gsl1_9, hsk_39, hsk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * gsl0_6[k]
                  + f_17 * gsk_3[k]
                  - f_14 * pc_y[k] * gsl1_6[k];

        t_52[k] = f_3 * pc_z[k] * hsk_39[k];

        t_53[k] = f_15 * gsk_5[k]
                  + f_3 * pc_y[k] * hsk_41[k];

        t_54[k] = pa_y[k] * gsl0_9[k]
                  - f_14 * pc_y[k] * gsl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pc_y, pc_z, gsl0_10, gsk_6, gsk_9, \
                         gsl1_10, hsi0_31, hsi1_31, hsk_42, hsk_43, \
                         hsk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * gsl0_10[k]
                  + f_18 * gsk_6[k]
                  - f_14 * pc_y[k] * gsl1_10[k];

        t_56[k] = f_3 * pc_z[k] * hsk_42[k];

        t_57[k] = f_4 * hsi0_31[k]
                  - f_5 * hsi1_31[k]
                  + f_3 * pc_z[k] * hsk_43[k];

        t_58[k] = f_15 * gsk_9[k]
                  + f_3 * pc_y[k] * hsk_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, gsl0_14, gsl0_15, gsk_10, \
                         gsl1_14, gsl1_15, hsi0_34, hsi1_34, hsk_46, \
                         hsk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * gsl0_14[k]
                  - f_14 * pc_y[k] * gsl1_14[k];

        t_60[k] = pa_y[k] * gsl0_15[k]
                  + f_0 * gsk_10[k]
                  - f_14 * pc_y[k] * gsl1_15[k];

        t_61[k] = f_3 * pc_z[k] * hsk_46[k];

        t_62[k] = f_4 * hsi0_34[k]
                  - f_5 * hsi1_34[k]
                  + f_3 * pc_z[k] * hsk_47[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pc_y, pc_z, gsl0_20, gsk_14, gsl1_20, \
                         hsi0_35, hsi1_35, hsk_48, hsk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_6 * hsi0_35[k]
                  - f_7 * hsi1_35[k]
                  + f_3 * pc_z[k] * hsk_48[k];

        t_64[k] = f_15 * gsk_14[k]
                  + f_3 * pc_y[k] * hsk_50[k];

        t_65[k] = pa_y[k] * gsl0_20[k]
                  - f_14 * pc_y[k] * gsl1_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pc_y, pc_z, gsl0_21, gsk_15, gsl1_21, \
                         hsi0_38, hsi1_38, hsk_51, hsk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * gsl0_21[k]
                  + f_19 * gsk_15[k]
                  - f_14 * pc_y[k] * gsl1_21[k];

        t_67[k] = f_3 * pc_z[k] * hsk_51[k];

        t_68[k] = f_4 * hsi0_38[k]
                  - f_5 * hsi1_38[k]
                  + f_3 * pc_z[k] * hsk_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_y, pc_z, gsk_20, hsi0_39, hsi0_40, hsi1_39, \
                         hsi1_40, hsk_53, hsk_54, hsk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * hsi0_39[k]
                  - f_7 * hsi1_39[k]
                  + f_3 * pc_z[k] * hsk_53[k];

        t_70[k] = f_8 * hsi0_40[k]
                  - f_9 * hsi1_40[k]
                  + f_3 * pc_z[k] * hsk_54[k];

        t_71[k] = f_15 * gsk_20[k]
                  + f_3 * pc_y[k] * hsk_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, pc_z, gsl0_27, gsk_64, \
                         gsk_66, gsl1_27, hsk_57, hsk_64, hsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * gsl0_27[k]
                  - f_14 * pc_y[k] * gsl1_27[k];

        t_73[k] = f_18 * gsk_64[k]
                  + f_3 * pc_x[k] * hsk_64[k];

        t_74[k] = f_3 * pc_z[k] * hsk_57[k];

        t_75[k] = f_18 * gsk_66[k]
                  + f_3 * pc_x[k] * hsk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, gsk_67, gsk_68, gsk_69, gsk_70, \
                         gsk_71, hsk_67, hsk_68, hsk_69, hsk_70, \
                         hsk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_18 * gsk_67[k]
                  + f_3 * pc_x[k] * hsk_67[k];

        t_77[k] = f_18 * gsk_68[k]
                  + f_3 * pc_x[k] * hsk_68[k];

        t_78[k] = f_18 * gsk_69[k]
                  + f_3 * pc_x[k] * hsk_69[k];

        t_79[k] = f_18 * gsk_70[k]
                  + f_3 * pc_x[k] * hsk_70[k];

        t_80[k] = f_18 * gsk_71[k]
                  + f_3 * pc_x[k] * hsk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_y, pc_z, gsk_28, hsi0_49, hsi0_50, \
                         hsi1_49, hsi1_50, hsk_64, hsk_65, hsk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_15 * gsk_28[k]
                  + f_1 * hsi0_49[k]
                  - f_2 * hsi1_49[k]
                  + f_3 * pc_y[k] * hsk_64[k];

        t_82[k] = f_3 * pc_z[k] * hsk_64[k];

        t_83[k] = f_4 * hsi0_49[k]
                  - f_5 * hsi1_49[k]
                  + f_3 * pc_z[k] * hsk_65[k];

        t_84[k] = f_6 * hsi0_50[k]
                  - f_7 * hsi1_50[k]
                  + f_3 * pc_z[k] * hsk_66[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_z, hsi0_51, hsi0_52, hsi0_53, hsi1_51, hsi1_52, \
                         hsi1_53, hsk_67, hsk_68, hsk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_8 * hsi0_51[k]
                  - f_9 * hsi1_51[k]
                  + f_3 * pc_z[k] * hsk_67[k];

        t_86[k] = f_10 * hsi0_52[k]
                  - f_11 * hsi1_52[k]
                  + f_3 * pc_z[k] * hsk_68[k];

        t_87[k] = f_12 * hsi0_53[k]
                  - f_13 * hsi1_53[k]
                  + f_3 * pc_z[k] * hsk_69[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pc_y, pc_z, gsl0_0, gsl0_44, \
                         gsk_35, gsl1_0, gsl1_44, hsk_71, hsk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * gsk_35[k]
                  + f_3 * pc_y[k] * hsk_71[k];

        t_89[k] = pa_y[k] * gsl0_44[k]
                  - f_14 * pc_y[k] * gsl1_44[k];

        t_90[k] = pa_z[k] * gsl0_0[k]
                  - f_14 * pc_z[k] * gsl1_0[k];

        t_91[k] = f_3 * pc_y[k] * hsk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_y, pc_z, gsl0_3, gsl0_5, gsk_0, \
                         gsk_2, gsl1_3, gsl1_5, hsk_72, hsk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * gsk_0[k]
                  + f_3 * pc_z[k] * hsk_72[k];

        t_93[k] = pa_z[k] * gsl0_3[k]
                  - f_14 * pc_z[k] * gsl1_3[k];

        t_94[k] = f_3 * pc_y[k] * hsk_74[k];

        t_95[k] = pa_z[k] * gsl0_5[k]
                  + f_16 * gsk_2[k]
                  - f_14 * pc_z[k] * gsl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pc_y, pc_z, gsl0_6, gsl0_9, gsk_5, \
                         gsl1_6, gsl1_9, hsi0_58, hsi1_58, hsk_76, \
                         hsk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * gsl0_6[k]
                  - f_14 * pc_z[k] * gsl1_6[k];

        t_97[k] = f_4 * hsi0_58[k]
                  - f_5 * hsi1_58[k]
                  + f_3 * pc_y[k] * hsk_76[k];

        t_98[k] = f_3 * pc_y[k] * hsk_77[k];

        t_99[k] = pa_z[k] * gsl0_9[k]
                  + f_17 * gsk_5[k]
                  - f_14 * pc_z[k] * gsl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pc_y, pc_z, gsl0_10, gsl1_10, \
                         hsi0_60, hsi0_61, hsi1_60, hsi1_61, hsk_79, hsk_80, \
                         hsk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * gsl0_10[k]
                   - f_14 * pc_z[k] * gsl1_10[k];

        t_101[k] = f_6 * hsi0_60[k]
                   - f_7 * hsi1_60[k]
                   + f_3 * pc_y[k] * hsk_79[k];

        t_102[k] = f_4 * hsi0_61[k]
                   - f_5 * hsi1_61[k]
                   + f_3 * pc_y[k] * hsk_80[k];

        t_103[k] = f_3 * pc_y[k] * hsk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pc_y, pc_z, gsl0_14, gsl0_15, gsk_9, \
                         gsl1_14, gsl1_15, hsi0_63, hsi1_63, hsk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_z[k] * gsl0_14[k]
                   + f_18 * gsk_9[k]
                   - f_14 * pc_z[k] * gsl1_14[k];

        t_105[k] = pa_z[k] * gsl0_15[k]
                   - f_14 * pc_z[k] * gsl1_15[k];

        t_106[k] = f_8 * hsi0_63[k]
                   - f_9 * hsi1_63[k]
                   + f_3 * pc_y[k] * hsk_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, hsi0_64, hsi0_65, hsi1_64, hsi1_65, \
                         hsk_84, hsk_85, hsk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * hsi0_64[k]
                   - f_7 * hsi1_64[k]
                   + f_3 * pc_y[k] * hsk_84[k];

        t_108[k] = f_4 * hsi0_65[k]
                   - f_5 * hsi1_65[k]
                   + f_3 * pc_y[k] * hsk_85[k];

        t_109[k] = f_3 * pc_y[k] * hsk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pc_y, pc_z, gsl0_20, gsl0_21, gsk_14, \
                         gsl1_20, gsl1_21, hsi0_67, hsi1_67, hsk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * gsl0_20[k]
                   + f_0 * gsk_14[k]
                   - f_14 * pc_z[k] * gsl1_20[k];

        t_111[k] = pa_z[k] * gsl0_21[k]
                   - f_14 * pc_z[k] * gsl1_21[k];

        t_112[k] = f_10 * hsi0_67[k]
                   - f_11 * hsi1_67[k]
                   + f_3 * pc_y[k] * hsk_88[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_y, hsi0_68, hsi0_69, hsi0_70, hsi1_68, \
                         hsi1_69, hsi1_70, hsk_89, hsk_90, hsk_91, \
                         hsk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * hsi0_68[k]
                   - f_9 * hsi1_68[k]
                   + f_3 * pc_y[k] * hsk_89[k];

        t_114[k] = f_6 * hsi0_69[k]
                   - f_7 * hsi1_69[k]
                   + f_3 * pc_y[k] * hsk_90[k];

        t_115[k] = f_4 * hsi0_70[k]
                   - f_5 * hsi1_70[k]
                   + f_3 * pc_y[k] * hsk_91[k];

        t_116[k] = f_3 * pc_y[k] * hsk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pc_x, pc_z, gsl0_27, gsk_20, \
                         gsk_100, gsk_101, gsk_102, gsl1_27, hsk_100, hsk_101, \
                         hsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_z[k] * gsl0_27[k]
                   + f_19 * gsk_20[k]
                   - f_14 * pc_z[k] * gsl1_27[k];

        t_118[k] = f_18 * gsk_100[k]
                   + f_3 * pc_x[k] * hsk_100[k];

        t_119[k] = f_18 * gsk_101[k]
                   + f_3 * pc_x[k] * hsk_101[k];

        t_120[k] = f_18 * gsk_102[k]
                   + f_3 * pc_x[k] * hsk_102[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, gsk_103, gsk_104, \
                         gsk_105, gsk_107, hsk_99, hsk_103, hsk_104, hsk_105, \
                         hsk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_18 * gsk_103[k]
                   + f_3 * pc_x[k] * hsk_103[k];

        t_122[k] = f_18 * gsk_104[k]
                   + f_3 * pc_x[k] * hsk_104[k];

        t_123[k] = f_18 * gsk_105[k]
                   + f_3 * pc_x[k] * hsk_105[k];

        t_124[k] = f_3 * pc_y[k] * hsk_99[k];

        t_125[k] = f_18 * gsk_107[k]
                   + f_3 * pc_x[k] * hsk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_z, pc_y, pc_z, gsl0_36, gsl1_36, hsi0_78, \
                         hsi0_79, hsi1_78, hsi1_79, hsk_101, hsk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * gsl0_36[k]
                   - f_14 * pc_z[k] * gsl1_36[k];

        t_127[k] = f_20 * hsi0_78[k]
                   - f_21 * hsi1_78[k]
                   + f_3 * pc_y[k] * hsk_101[k];

        t_128[k] = f_12 * hsi0_79[k]
                   - f_13 * hsi1_79[k]
                   + f_3 * pc_y[k] * hsk_102[k];
    }
}

static auto
compute_prim_hsl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsl0,
                                                          const size_t gsk, const size_t gsl1,
                                                          const size_t hsi0, const size_t hsi1,
                                                          const size_t hsk, const size_t ncols,
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

    const auto *gsl0_48 = buffer.data(gsl0 + 48);
    const auto *gsl0_51 = buffer.data(gsl0 + 51);
    const auto *gsl0_55 = buffer.data(gsl0 + 55);
    const auto *gsl0_60 = buffer.data(gsl0 + 60);
    const auto *gsl0_66 = buffer.data(gsl0 + 66);
    const auto *gsl0_81 = buffer.data(gsl0 + 81);
    const auto *gsl0_90 = buffer.data(gsl0 + 90);
    const auto *gsl0_95 = buffer.data(gsl0 + 95);
    const auto *gsl0_99 = buffer.data(gsl0 + 99);
    const auto *gsl0_102 = buffer.data(gsl0 + 102);
    const auto *gsl0_104 = buffer.data(gsl0 + 104);
    const auto *gsl0_107 = buffer.data(gsl0 + 107);
    const auto *gsl0_108 = buffer.data(gsl0 + 108);
    const auto *gsl0_110 = buffer.data(gsl0 + 110);
    const auto *gsl0_113 = buffer.data(gsl0 + 113);
    const auto *gsl0_114 = buffer.data(gsl0 + 114);
    const auto *gsl0_115 = buffer.data(gsl0 + 115);
    const auto *gsl0_117 = buffer.data(gsl0 + 117);
    const auto *gsl0_134 = buffer.data(gsl0 + 134);

    const auto *gsk_35 = buffer.data(gsk + 35);
    const auto *gsk_36 = buffer.data(gsk + 36);
    const auto *gsk_39 = buffer.data(gsk + 39);
    const auto *gsk_41 = buffer.data(gsk + 41);
    const auto *gsk_42 = buffer.data(gsk + 42);
    const auto *gsk_45 = buffer.data(gsk + 45);
    const auto *gsk_46 = buffer.data(gsk + 46);
    const auto *gsk_50 = buffer.data(gsk + 50);
    const auto *gsk_51 = buffer.data(gsk + 51);
    const auto *gsk_56 = buffer.data(gsk + 56);
    const auto *gsk_64 = buffer.data(gsk + 64);
    const auto *gsk_71 = buffer.data(gsk + 71);
    const auto *gsk_72 = buffer.data(gsk + 72);
    const auto *gsk_74 = buffer.data(gsk + 74);
    const auto *gsk_77 = buffer.data(gsk + 77);
    const auto *gsk_80 = buffer.data(gsk + 80);
    const auto *gsk_81 = buffer.data(gsk + 81);
    const auto *gsk_84 = buffer.data(gsk + 84);
    const auto *gsk_85 = buffer.data(gsk + 85);
    const auto *gsk_86 = buffer.data(gsk + 86);
    const auto *gsk_89 = buffer.data(gsk + 89);
    const auto *gsk_90 = buffer.data(gsk + 90);
    const auto *gsk_91 = buffer.data(gsk + 91);
    const auto *gsk_92 = buffer.data(gsk + 92);
    const auto *gsk_102 = buffer.data(gsk + 102);
    const auto *gsk_103 = buffer.data(gsk + 103);
    const auto *gsk_104 = buffer.data(gsk + 104);
    const auto *gsk_105 = buffer.data(gsk + 105);
    const auto *gsk_106 = buffer.data(gsk + 106);
    const auto *gsk_107 = buffer.data(gsk + 107);
    const auto *gsk_108 = buffer.data(gsk + 108);
    const auto *gsk_111 = buffer.data(gsk + 111);
    const auto *gsk_114 = buffer.data(gsk + 114);
    const auto *gsk_118 = buffer.data(gsk + 118);
    const auto *gsk_123 = buffer.data(gsk + 123);
    const auto *gsk_129 = buffer.data(gsk + 129);
    const auto *gsk_136 = buffer.data(gsk + 136);
    const auto *gsk_138 = buffer.data(gsk + 138);
    const auto *gsk_139 = buffer.data(gsk + 139);
    const auto *gsk_140 = buffer.data(gsk + 140);
    const auto *gsk_141 = buffer.data(gsk + 141);
    const auto *gsk_142 = buffer.data(gsk + 142);
    const auto *gsk_143 = buffer.data(gsk + 143);
    const auto *gsk_172 = buffer.data(gsk + 172);
    const auto *gsk_173 = buffer.data(gsk + 173);
    const auto *gsk_174 = buffer.data(gsk + 174);
    const auto *gsk_175 = buffer.data(gsk + 175);
    const auto *gsk_176 = buffer.data(gsk + 176);
    const auto *gsk_177 = buffer.data(gsk + 177);
    const auto *gsk_178 = buffer.data(gsk + 178);
    const auto *gsk_179 = buffer.data(gsk + 179);
    const auto *gsk_180 = buffer.data(gsk + 180);
    const auto *gsk_185 = buffer.data(gsk + 185);
    const auto *gsk_189 = buffer.data(gsk + 189);
    const auto *gsk_194 = buffer.data(gsk + 194);
    const auto *gsk_200 = buffer.data(gsk + 200);

    const auto *gsl1_48 = buffer.data(gsl1 + 48);
    const auto *gsl1_51 = buffer.data(gsl1 + 51);
    const auto *gsl1_55 = buffer.data(gsl1 + 55);
    const auto *gsl1_60 = buffer.data(gsl1 + 60);
    const auto *gsl1_66 = buffer.data(gsl1 + 66);
    const auto *gsl1_81 = buffer.data(gsl1 + 81);
    const auto *gsl1_90 = buffer.data(gsl1 + 90);
    const auto *gsl1_95 = buffer.data(gsl1 + 95);
    const auto *gsl1_99 = buffer.data(gsl1 + 99);
    const auto *gsl1_102 = buffer.data(gsl1 + 102);
    const auto *gsl1_104 = buffer.data(gsl1 + 104);
    const auto *gsl1_107 = buffer.data(gsl1 + 107);
    const auto *gsl1_108 = buffer.data(gsl1 + 108);
    const auto *gsl1_110 = buffer.data(gsl1 + 110);
    const auto *gsl1_113 = buffer.data(gsl1 + 113);
    const auto *gsl1_114 = buffer.data(gsl1 + 114);
    const auto *gsl1_115 = buffer.data(gsl1 + 115);
    const auto *gsl1_117 = buffer.data(gsl1 + 117);
    const auto *gsl1_134 = buffer.data(gsl1 + 134);

    const auto *hsi0_80 = buffer.data(hsi0 + 80);
    const auto *hsi0_81 = buffer.data(hsi0 + 81);
    const auto *hsi0_82 = buffer.data(hsi0 + 82);
    const auto *hsi0_83 = buffer.data(hsi0 + 83);
    const auto *hsi0_84 = buffer.data(hsi0 + 84);
    const auto *hsi0_86 = buffer.data(hsi0 + 86);
    const auto *hsi0_87 = buffer.data(hsi0 + 87);
    const auto *hsi0_89 = buffer.data(hsi0 + 89);
    const auto *hsi0_90 = buffer.data(hsi0 + 90);
    const auto *hsi0_91 = buffer.data(hsi0 + 91);
    const auto *hsi0_93 = buffer.data(hsi0 + 93);
    const auto *hsi0_94 = buffer.data(hsi0 + 94);
    const auto *hsi0_95 = buffer.data(hsi0 + 95);
    const auto *hsi0_96 = buffer.data(hsi0 + 96);
    const auto *hsi0_98 = buffer.data(hsi0 + 98);
    const auto *hsi0_99 = buffer.data(hsi0 + 99);
    const auto *hsi0_105 = buffer.data(hsi0 + 105);
    const auto *hsi0_106 = buffer.data(hsi0 + 106);
    const auto *hsi0_107 = buffer.data(hsi0 + 107);
    const auto *hsi0_108 = buffer.data(hsi0 + 108);
    const auto *hsi0_109 = buffer.data(hsi0 + 109);
    const auto *hsi0_111 = buffer.data(hsi0 + 111);
    const auto *hsi0_135 = buffer.data(hsi0 + 135);
    const auto *hsi0_136 = buffer.data(hsi0 + 136);
    const auto *hsi0_137 = buffer.data(hsi0 + 137);
    const auto *hsi0_138 = buffer.data(hsi0 + 138);
    const auto *hsi0_139 = buffer.data(hsi0 + 139);
    const auto *hsi0_140 = buffer.data(hsi0 + 140);
    const auto *hsi0_141 = buffer.data(hsi0 + 141);
    const auto *hsi0_142 = buffer.data(hsi0 + 142);
    const auto *hsi0_143 = buffer.data(hsi0 + 143);
    const auto *hsi0_144 = buffer.data(hsi0 + 144);
    const auto *hsi0_145 = buffer.data(hsi0 + 145);
    const auto *hsi0_146 = buffer.data(hsi0 + 146);
    const auto *hsi0_147 = buffer.data(hsi0 + 147);
    const auto *hsi0_148 = buffer.data(hsi0 + 148);
    const auto *hsi0_149 = buffer.data(hsi0 + 149);
    const auto *hsi0_154 = buffer.data(hsi0 + 154);
    const auto *hsi0_160 = buffer.data(hsi0 + 160);

    const auto *hsi1_80 = buffer.data(hsi1 + 80);
    const auto *hsi1_81 = buffer.data(hsi1 + 81);
    const auto *hsi1_82 = buffer.data(hsi1 + 82);
    const auto *hsi1_83 = buffer.data(hsi1 + 83);
    const auto *hsi1_84 = buffer.data(hsi1 + 84);
    const auto *hsi1_86 = buffer.data(hsi1 + 86);
    const auto *hsi1_87 = buffer.data(hsi1 + 87);
    const auto *hsi1_89 = buffer.data(hsi1 + 89);
    const auto *hsi1_90 = buffer.data(hsi1 + 90);
    const auto *hsi1_91 = buffer.data(hsi1 + 91);
    const auto *hsi1_93 = buffer.data(hsi1 + 93);
    const auto *hsi1_94 = buffer.data(hsi1 + 94);
    const auto *hsi1_95 = buffer.data(hsi1 + 95);
    const auto *hsi1_96 = buffer.data(hsi1 + 96);
    const auto *hsi1_98 = buffer.data(hsi1 + 98);
    const auto *hsi1_99 = buffer.data(hsi1 + 99);
    const auto *hsi1_105 = buffer.data(hsi1 + 105);
    const auto *hsi1_106 = buffer.data(hsi1 + 106);
    const auto *hsi1_107 = buffer.data(hsi1 + 107);
    const auto *hsi1_108 = buffer.data(hsi1 + 108);
    const auto *hsi1_109 = buffer.data(hsi1 + 109);
    const auto *hsi1_111 = buffer.data(hsi1 + 111);
    const auto *hsi1_135 = buffer.data(hsi1 + 135);
    const auto *hsi1_136 = buffer.data(hsi1 + 136);
    const auto *hsi1_137 = buffer.data(hsi1 + 137);
    const auto *hsi1_138 = buffer.data(hsi1 + 138);
    const auto *hsi1_139 = buffer.data(hsi1 + 139);
    const auto *hsi1_140 = buffer.data(hsi1 + 140);
    const auto *hsi1_141 = buffer.data(hsi1 + 141);
    const auto *hsi1_142 = buffer.data(hsi1 + 142);
    const auto *hsi1_143 = buffer.data(hsi1 + 143);
    const auto *hsi1_144 = buffer.data(hsi1 + 144);
    const auto *hsi1_145 = buffer.data(hsi1 + 145);
    const auto *hsi1_146 = buffer.data(hsi1 + 146);
    const auto *hsi1_147 = buffer.data(hsi1 + 147);
    const auto *hsi1_148 = buffer.data(hsi1 + 148);
    const auto *hsi1_149 = buffer.data(hsi1 + 149);
    const auto *hsi1_154 = buffer.data(hsi1 + 154);
    const auto *hsi1_160 = buffer.data(hsi1 + 160);

    const auto *hsk_103 = buffer.data(hsk + 103);
    const auto *hsk_104 = buffer.data(hsk + 104);
    const auto *hsk_105 = buffer.data(hsk + 105);
    const auto *hsk_106 = buffer.data(hsk + 106);
    const auto *hsk_107 = buffer.data(hsk + 107);
    const auto *hsk_108 = buffer.data(hsk + 108);
    const auto *hsk_109 = buffer.data(hsk + 109);
    const auto *hsk_110 = buffer.data(hsk + 110);
    const auto *hsk_111 = buffer.data(hsk + 111);
    const auto *hsk_113 = buffer.data(hsk + 113);
    const auto *hsk_114 = buffer.data(hsk + 114);
    const auto *hsk_115 = buffer.data(hsk + 115);
    const auto *hsk_117 = buffer.data(hsk + 117);
    const auto *hsk_118 = buffer.data(hsk + 118);
    const auto *hsk_119 = buffer.data(hsk + 119);
    const auto *hsk_120 = buffer.data(hsk + 120);
    const auto *hsk_122 = buffer.data(hsk + 122);
    const auto *hsk_123 = buffer.data(hsk + 123);
    const auto *hsk_124 = buffer.data(hsk + 124);
    const auto *hsk_125 = buffer.data(hsk + 125);
    const auto *hsk_126 = buffer.data(hsk + 126);
    const auto *hsk_128 = buffer.data(hsk + 128);
    const auto *hsk_129 = buffer.data(hsk + 129);
    const auto *hsk_136 = buffer.data(hsk + 136);
    const auto *hsk_137 = buffer.data(hsk + 137);
    const auto *hsk_138 = buffer.data(hsk + 138);
    const auto *hsk_139 = buffer.data(hsk + 139);
    const auto *hsk_140 = buffer.data(hsk + 140);
    const auto *hsk_141 = buffer.data(hsk + 141);
    const auto *hsk_142 = buffer.data(hsk + 142);
    const auto *hsk_143 = buffer.data(hsk + 143);
    const auto *hsk_144 = buffer.data(hsk + 144);
    const auto *hsk_146 = buffer.data(hsk + 146);
    const auto *hsk_147 = buffer.data(hsk + 147);
    const auto *hsk_149 = buffer.data(hsk + 149);
    const auto *hsk_150 = buffer.data(hsk + 150);
    const auto *hsk_153 = buffer.data(hsk + 153);
    const auto *hsk_154 = buffer.data(hsk + 154);
    const auto *hsk_158 = buffer.data(hsk + 158);
    const auto *hsk_159 = buffer.data(hsk + 159);
    const auto *hsk_164 = buffer.data(hsk + 164);
    const auto *hsk_172 = buffer.data(hsk + 172);
    const auto *hsk_173 = buffer.data(hsk + 173);
    const auto *hsk_174 = buffer.data(hsk + 174);
    const auto *hsk_175 = buffer.data(hsk + 175);
    const auto *hsk_176 = buffer.data(hsk + 176);
    const auto *hsk_177 = buffer.data(hsk + 177);
    const auto *hsk_178 = buffer.data(hsk + 178);
    const auto *hsk_179 = buffer.data(hsk + 179);
    const auto *hsk_180 = buffer.data(hsk + 180);
    const auto *hsk_181 = buffer.data(hsk + 181);
    const auto *hsk_182 = buffer.data(hsk + 182);
    const auto *hsk_183 = buffer.data(hsk + 183);
    const auto *hsk_184 = buffer.data(hsk + 184);
    const auto *hsk_185 = buffer.data(hsk + 185);
    const auto *hsk_186 = buffer.data(hsk + 186);
    const auto *hsk_187 = buffer.data(hsk + 187);
    const auto *hsk_188 = buffer.data(hsk + 188);
    const auto *hsk_189 = buffer.data(hsk + 189);
    const auto *hsk_190 = buffer.data(hsk + 190);
    const auto *hsk_191 = buffer.data(hsk + 191);
    const auto *hsk_192 = buffer.data(hsk + 192);
    const auto *hsk_193 = buffer.data(hsk + 193);
    const auto *hsk_194 = buffer.data(hsk + 194);
    const auto *hsk_200 = buffer.data(hsk + 200);

#pragma omp simd aligned(t_129, t_130, t_131, pc_y, hsi0_80, hsi0_81, hsi0_82, hsi1_80, \
                         hsi1_81, hsi1_82, hsk_103, hsk_104, hsk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * hsi0_80[k]
                   - f_11 * hsi1_80[k]
                   + f_3 * pc_y[k] * hsk_103[k];

        t_130[k] = f_8 * hsi0_81[k]
                   - f_9 * hsi1_81[k]
                   + f_3 * pc_y[k] * hsk_104[k];

        t_131[k] = f_6 * hsi0_82[k]
                   - f_7 * hsi1_82[k]
                   + f_3 * pc_y[k] * hsk_105[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, pc_y, pc_z, gsk_35, gsk_108, \
                         hsi0_83, hsi0_84, hsi1_83, hsi1_84, hsk_106, hsk_107, \
                         hsk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * hsi0_83[k]
                   - f_5 * hsi1_83[k]
                   + f_3 * pc_y[k] * hsk_106[k];

        t_133[k] = f_3 * pc_y[k] * hsk_107[k];

        t_134[k] = f_15 * gsk_35[k]
                   + f_1 * hsi0_83[k]
                   - f_2 * hsi1_83[k]
                   + f_3 * pc_z[k] * hsk_107[k];

        t_135[k] = f_17 * gsk_108[k]
                   + f_1 * hsi0_84[k]
                   - f_2 * hsi1_84[k]
                   + f_3 * pc_x[k] * hsk_108[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, gsk_36, gsk_111, \
                         hsi0_87, hsi1_87, hsk_108, hsk_109, hsk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_16 * gsk_36[k]
                   + f_3 * pc_y[k] * hsk_108[k];

        t_137[k] = f_3 * pc_z[k] * hsk_108[k];

        t_138[k] = f_17 * gsk_111[k]
                   + f_12 * hsi0_87[k]
                   - f_13 * hsi1_87[k]
                   + f_3 * pc_x[k] * hsk_111[k];

        t_139[k] = f_3 * pc_z[k] * hsk_109[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pc_x, pc_z, gsk_114, hsi0_84, hsi0_90, hsi1_84, \
                         hsi1_90, hsk_110, hsk_111, hsk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * hsi0_84[k]
                   - f_5 * hsi1_84[k]
                   + f_3 * pc_z[k] * hsk_110[k];

        t_141[k] = f_17 * gsk_114[k]
                   + f_10 * hsi0_90[k]
                   - f_11 * hsi1_90[k]
                   + f_3 * pc_x[k] * hsk_114[k];

        t_142[k] = f_3 * pc_z[k] * hsk_111[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, gsk_41, gsk_118, \
                         hsi0_86, hsi0_94, hsi1_86, hsi1_94, hsk_113, hsk_114, \
                         hsk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * gsk_41[k]
                   + f_3 * pc_y[k] * hsk_113[k];

        t_144[k] = f_6 * hsi0_86[k]
                   - f_7 * hsi1_86[k]
                   + f_3 * pc_z[k] * hsk_113[k];

        t_145[k] = f_17 * gsk_118[k]
                   + f_8 * hsi0_94[k]
                   - f_9 * hsi1_94[k]
                   + f_3 * pc_x[k] * hsk_118[k];

        t_146[k] = f_3 * pc_z[k] * hsk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_y, pc_z, gsk_45, hsi0_87, hsi0_89, hsi1_87, \
                         hsi1_89, hsk_115, hsk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * hsi0_87[k]
                   - f_5 * hsi1_87[k]
                   + f_3 * pc_z[k] * hsk_115[k];

        t_148[k] = f_16 * gsk_45[k]
                   + f_3 * pc_y[k] * hsk_117[k];

        t_149[k] = f_8 * hsi0_89[k]
                   - f_9 * hsi1_89[k]
                   + f_3 * pc_z[k] * hsk_117[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, gsk_123, hsi0_90, hsi0_99, hsi1_90, \
                         hsi1_99, hsk_118, hsk_119, hsk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_17 * gsk_123[k]
                   + f_6 * hsi0_99[k]
                   - f_7 * hsi1_99[k]
                   + f_3 * pc_x[k] * hsk_123[k];

        t_151[k] = f_3 * pc_z[k] * hsk_118[k];

        t_152[k] = f_4 * hsi0_90[k]
                   - f_5 * hsi1_90[k]
                   + f_3 * pc_z[k] * hsk_119[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_y, pc_z, gsk_50, hsi0_91, hsi0_93, hsi1_91, \
                         hsi1_93, hsk_120, hsk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * hsi0_91[k]
                   - f_7 * hsi1_91[k]
                   + f_3 * pc_z[k] * hsk_120[k];

        t_154[k] = f_16 * gsk_50[k]
                   + f_3 * pc_y[k] * hsk_122[k];

        t_155[k] = f_10 * hsi0_93[k]
                   - f_11 * hsi1_93[k]
                   + f_3 * pc_z[k] * hsk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, gsk_129, hsi0_94, hsi0_105, hsi1_94, \
                         hsi1_105, hsk_123, hsk_124, hsk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_17 * gsk_129[k]
                   + f_4 * hsi0_105[k]
                   - f_5 * hsi1_105[k]
                   + f_3 * pc_x[k] * hsk_129[k];

        t_157[k] = f_3 * pc_z[k] * hsk_123[k];

        t_158[k] = f_4 * hsi0_94[k]
                   - f_5 * hsi1_94[k]
                   + f_3 * pc_z[k] * hsk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pc_y, pc_z, gsk_56, hsi0_95, hsi0_96, \
                         hsi0_98, hsi1_95, hsi1_96, hsi1_98, hsk_125, hsk_126, \
                         hsk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_6 * hsi0_95[k]
                   - f_7 * hsi1_95[k]
                   + f_3 * pc_z[k] * hsk_125[k];

        t_160[k] = f_8 * hsi0_96[k]
                   - f_9 * hsi1_96[k]
                   + f_3 * pc_z[k] * hsk_126[k];

        t_161[k] = f_16 * gsk_56[k]
                   + f_3 * pc_y[k] * hsk_128[k];

        t_162[k] = f_12 * hsi0_98[k]
                   - f_13 * hsi1_98[k]
                   + f_3 * pc_z[k] * hsk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pc_x, pc_z, gsk_136, gsk_138, \
                         gsk_139, gsk_140, hsk_129, hsk_136, hsk_138, hsk_139, \
                         hsk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_17 * gsk_136[k]
                   + f_3 * pc_x[k] * hsk_136[k];

        t_164[k] = f_3 * pc_z[k] * hsk_129[k];

        t_165[k] = f_17 * gsk_138[k]
                   + f_3 * pc_x[k] * hsk_138[k];

        t_166[k] = f_17 * gsk_139[k]
                   + f_3 * pc_x[k] * hsk_139[k];

        t_167[k] = f_17 * gsk_140[k]
                   + f_3 * pc_x[k] * hsk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, gsk_64, gsk_141, gsk_142, \
                         gsk_143, hsi0_105, hsi1_105, hsk_136, hsk_141, hsk_142, \
                         hsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_17 * gsk_141[k]
                   + f_3 * pc_x[k] * hsk_141[k];

        t_169[k] = f_17 * gsk_142[k]
                   + f_3 * pc_x[k] * hsk_142[k];

        t_170[k] = f_17 * gsk_143[k]
                   + f_3 * pc_x[k] * hsk_143[k];

        t_171[k] = f_16 * gsk_64[k]
                   + f_1 * hsi0_105[k]
                   - f_2 * hsi1_105[k]
                   + f_3 * pc_y[k] * hsk_136[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_z, hsi0_105, hsi0_106, hsi0_107, \
                         hsi1_105, hsi1_106, hsi1_107, hsk_136, hsk_137, hsk_138, \
                         hsk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * hsk_136[k];

        t_173[k] = f_4 * hsi0_105[k]
                   - f_5 * hsi1_105[k]
                   + f_3 * pc_z[k] * hsk_137[k];

        t_174[k] = f_6 * hsi0_106[k]
                   - f_7 * hsi1_106[k]
                   + f_3 * pc_z[k] * hsk_138[k];

        t_175[k] = f_8 * hsi0_107[k]
                   - f_9 * hsi1_107[k]
                   + f_3 * pc_z[k] * hsk_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_y, pc_z, gsk_71, hsi0_108, hsi0_109, \
                         hsi0_111, hsi1_108, hsi1_109, hsi1_111, hsk_140, hsk_141, \
                         hsk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * hsi0_108[k]
                   - f_11 * hsi1_108[k]
                   + f_3 * pc_z[k] * hsk_140[k];

        t_177[k] = f_12 * hsi0_109[k]
                   - f_13 * hsi1_109[k]
                   + f_3 * pc_z[k] * hsk_141[k];

        t_178[k] = f_16 * gsk_71[k]
                   + f_3 * pc_y[k] * hsk_143[k];

        t_179[k] = f_1 * hsi0_111[k]
                   - f_2 * hsi1_111[k]
                   + f_3 * pc_z[k] * hsk_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, gsl0_48, gsl0_90, \
                         gsk_36, gsk_72, gsl1_48, gsl1_90, hsk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * gsl0_90[k]
                   - f_14 * pc_y[k] * gsl1_90[k];

        t_181[k] = f_15 * gsk_72[k]
                   + f_3 * pc_y[k] * hsk_144[k];

        t_182[k] = f_15 * gsk_36[k]
                   + f_3 * pc_z[k] * hsk_144[k];

        t_183[k] = pa_z[k] * gsl0_48[k]
                   - f_14 * pc_z[k] * gsl1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pa_z, pc_y, pc_z, gsl0_51, gsl0_95, \
                         gsk_39, gsk_74, gsl1_51, gsl1_95, hsk_146, \
                         hsk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_15 * gsk_74[k]
                   + f_3 * pc_y[k] * hsk_146[k];

        t_185[k] = pa_y[k] * gsl0_95[k]
                   - f_14 * pc_y[k] * gsl1_95[k];

        t_186[k] = pa_z[k] * gsl0_51[k]
                   - f_14 * pc_z[k] * gsl1_51[k];

        t_187[k] = f_15 * gsk_39[k]
                   + f_3 * pc_z[k] * hsk_147[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pc_y, pc_z, gsl0_55, gsl0_99, \
                         gsk_42, gsk_77, gsl1_55, gsl1_99, hsk_149, \
                         hsk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_15 * gsk_77[k]
                   + f_3 * pc_y[k] * hsk_149[k];

        t_189[k] = pa_y[k] * gsl0_99[k]
                   - f_14 * pc_y[k] * gsl1_99[k];

        t_190[k] = pa_z[k] * gsl0_55[k]
                   - f_14 * pc_z[k] * gsl1_55[k];

        t_191[k] = f_15 * gsk_42[k]
                   + f_3 * pc_z[k] * hsk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pc_y, gsl0_102, gsl0_104, gsk_80, gsk_81, \
                         gsl1_102, gsl1_104, hsk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_y[k] * gsl0_102[k]
                   + f_16 * gsk_80[k]
                   - f_14 * pc_y[k] * gsl1_102[k];

        t_193[k] = f_15 * gsk_81[k]
                   + f_3 * pc_y[k] * hsk_153[k];

        t_194[k] = pa_y[k] * gsl0_104[k]
                   - f_14 * pc_y[k] * gsl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_y, pa_z, pc_y, pc_z, gsl0_60, gsl0_107, \
                         gsk_46, gsk_84, gsl1_60, gsl1_107, hsk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * gsl0_60[k]
                   - f_14 * pc_z[k] * gsl1_60[k];

        t_196[k] = f_15 * gsk_46[k]
                   + f_3 * pc_z[k] * hsk_154[k];

        t_197[k] = pa_y[k] * gsl0_107[k]
                   + f_17 * gsk_84[k]
                   - f_14 * pc_y[k] * gsl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_y, pc_y, gsl0_108, gsl0_110, gsk_85, gsk_86, \
                         gsl1_108, gsl1_110, hsk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_y[k] * gsl0_108[k]
                   + f_16 * gsk_85[k]
                   - f_14 * pc_y[k] * gsl1_108[k];

        t_199[k] = f_15 * gsk_86[k]
                   + f_3 * pc_y[k] * hsk_158[k];

        t_200[k] = pa_y[k] * gsl0_110[k]
                   - f_14 * pc_y[k] * gsl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_y, pa_z, pc_y, pc_z, gsl0_66, gsl0_113, \
                         gsk_51, gsk_89, gsl1_66, gsl1_113, hsk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_z[k] * gsl0_66[k]
                   - f_14 * pc_z[k] * gsl1_66[k];

        t_202[k] = f_15 * gsk_51[k]
                   + f_3 * pc_z[k] * hsk_159[k];

        t_203[k] = pa_y[k] * gsl0_113[k]
                   + f_18 * gsk_89[k]
                   - f_14 * pc_y[k] * gsl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_y, gsl0_114, gsl0_115, gsl0_117, \
                         gsk_90, gsk_91, gsk_92, gsl1_114, gsl1_115, gsl1_117, \
                         hsk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * gsl0_114[k]
                   + f_17 * gsk_90[k]
                   - f_14 * pc_y[k] * gsl1_114[k];

        t_205[k] = pa_y[k] * gsl0_115[k]
                   + f_16 * gsk_91[k]
                   - f_14 * pc_y[k] * gsl1_115[k];

        t_206[k] = f_15 * gsk_92[k]
                   + f_3 * pc_y[k] * hsk_164[k];

        t_207[k] = pa_y[k] * gsl0_117[k]
                   - f_14 * pc_y[k] * gsl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, gsk_172, gsk_173, gsk_174, \
                         gsk_175, gsk_176, hsk_172, hsk_173, hsk_174, hsk_175, \
                         hsk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_17 * gsk_172[k]
                   + f_3 * pc_x[k] * hsk_172[k];

        t_209[k] = f_17 * gsk_173[k]
                   + f_3 * pc_x[k] * hsk_173[k];

        t_210[k] = f_17 * gsk_174[k]
                   + f_3 * pc_x[k] * hsk_174[k];

        t_211[k] = f_17 * gsk_175[k]
                   + f_3 * pc_x[k] * hsk_175[k];

        t_212[k] = f_17 * gsk_176[k]
                   + f_3 * pc_x[k] * hsk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_z, pc_x, pc_z, gsl0_81, gsk_177, \
                         gsk_178, gsk_179, gsl1_81, hsk_177, hsk_178, \
                         hsk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_17 * gsk_177[k]
                   + f_3 * pc_x[k] * hsk_177[k];

        t_214[k] = f_17 * gsk_178[k]
                   + f_3 * pc_x[k] * hsk_178[k];

        t_215[k] = f_17 * gsk_179[k]
                   + f_3 * pc_x[k] * hsk_179[k];

        t_216[k] = pa_z[k] * gsl0_81[k]
                   - f_14 * pc_z[k] * gsl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, gsk_64, gsk_102, gsk_103, hsi0_135, \
                         hsi0_136, hsi1_135, hsi1_136, hsk_172, hsk_174, \
                         hsk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * gsk_64[k]
                   + f_3 * pc_z[k] * hsk_172[k];

        t_218[k] = f_15 * gsk_102[k]
                   + f_12 * hsi0_135[k]
                   - f_13 * hsi1_135[k]
                   + f_3 * pc_y[k] * hsk_174[k];

        t_219[k] = f_15 * gsk_103[k]
                   + f_10 * hsi0_136[k]
                   - f_11 * hsi1_136[k]
                   + f_3 * pc_y[k] * hsk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, gsk_104, gsk_105, gsk_106, hsi0_137, \
                         hsi0_138, hsi0_139, hsi1_137, hsi1_138, hsi1_139, hsk_176, hsk_177, \
                         hsk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * gsk_104[k]
                   + f_8 * hsi0_137[k]
                   - f_9 * hsi1_137[k]
                   + f_3 * pc_y[k] * hsk_176[k];

        t_221[k] = f_15 * gsk_105[k]
                   + f_6 * hsi0_138[k]
                   - f_7 * hsi1_138[k]
                   + f_3 * pc_y[k] * hsk_177[k];

        t_222[k] = f_15 * gsk_106[k]
                   + f_4 * hsi0_139[k]
                   - f_5 * hsi1_139[k]
                   + f_3 * pc_y[k] * hsk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_y, pc_x, pc_y, gsl0_134, gsk_107, \
                         gsk_180, gsl1_134, hsi0_140, hsi1_140, hsk_179, \
                         hsk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * gsk_107[k]
                   + f_3 * pc_y[k] * hsk_179[k];

        t_224[k] = pa_y[k] * gsl0_134[k]
                   - f_14 * pc_y[k] * gsl1_134[k];

        t_225[k] = f_17 * gsk_180[k]
                   + f_1 * hsi0_140[k]
                   - f_2 * hsi1_140[k]
                   + f_3 * pc_x[k] * hsk_180[k];

        t_226[k] = f_3 * pc_y[k] * hsk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_y, pc_z, gsk_72, hsi0_140, hsi1_140, hsk_180, \
                         hsk_181, hsk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * gsk_72[k]
                   + f_3 * pc_z[k] * hsk_180[k];

        t_228[k] = f_4 * hsi0_140[k]
                   - f_5 * hsi1_140[k]
                   + f_3 * pc_y[k] * hsk_181[k];

        t_229[k] = f_3 * pc_y[k] * hsk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, gsk_185, hsi0_141, hsi0_142, \
                         hsi0_145, hsi1_141, hsi1_142, hsi1_145, hsk_183, hsk_184, \
                         hsk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_17 * gsk_185[k]
                   + f_12 * hsi0_145[k]
                   - f_13 * hsi1_145[k]
                   + f_3 * pc_x[k] * hsk_185[k];

        t_231[k] = f_6 * hsi0_141[k]
                   - f_7 * hsi1_141[k]
                   + f_3 * pc_y[k] * hsk_183[k];

        t_232[k] = f_4 * hsi0_142[k]
                   - f_5 * hsi1_142[k]
                   + f_3 * pc_y[k] * hsk_184[k];

        t_233[k] = f_3 * pc_y[k] * hsk_185[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, gsk_189, hsi0_143, hsi0_144, \
                         hsi0_149, hsi1_143, hsi1_144, hsi1_149, hsk_186, hsk_187, \
                         hsk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_17 * gsk_189[k]
                   + f_10 * hsi0_149[k]
                   - f_11 * hsi1_149[k]
                   + f_3 * pc_x[k] * hsk_189[k];

        t_235[k] = f_8 * hsi0_143[k]
                   - f_9 * hsi1_143[k]
                   + f_3 * pc_y[k] * hsk_186[k];

        t_236[k] = f_6 * hsi0_144[k]
                   - f_7 * hsi1_144[k]
                   + f_3 * pc_y[k] * hsk_187[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_x, pc_y, gsk_194, hsi0_145, hsi0_154, \
                         hsi1_145, hsi1_154, hsk_188, hsk_189, \
                         hsk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * hsi0_145[k]
                   - f_5 * hsi1_145[k]
                   + f_3 * pc_y[k] * hsk_188[k];

        t_238[k] = f_3 * pc_y[k] * hsk_189[k];

        t_239[k] = f_17 * gsk_194[k]
                   + f_8 * hsi0_154[k]
                   - f_9 * hsi1_154[k]
                   + f_3 * pc_x[k] * hsk_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pc_y, hsi0_146, hsi0_147, hsi0_148, hsi1_146, \
                         hsi1_147, hsi1_148, hsk_190, hsk_191, \
                         hsk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * hsi0_146[k]
                   - f_11 * hsi1_146[k]
                   + f_3 * pc_y[k] * hsk_190[k];

        t_241[k] = f_8 * hsi0_147[k]
                   - f_9 * hsi1_147[k]
                   + f_3 * pc_y[k] * hsk_191[k];

        t_242[k] = f_6 * hsi0_148[k]
                   - f_7 * hsi1_148[k]
                   + f_3 * pc_y[k] * hsk_192[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pc_x, pc_y, gsk_200, hsi0_149, hsi0_160, \
                         hsi1_149, hsi1_160, hsk_193, hsk_194, \
                         hsk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_4 * hsi0_149[k]
                   - f_5 * hsi1_149[k]
                   + f_3 * pc_y[k] * hsk_193[k];

        t_244[k] = f_3 * pc_y[k] * hsk_194[k];

        t_245[k] = f_17 * gsk_200[k]
                   + f_6 * hsi0_160[k]
                   - f_7 * hsi1_160[k]
                   + f_3 * pc_x[k] * hsk_200[k];
    }
}

static auto
compute_prim_hsl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsl0,
                                                          const size_t gsk, const size_t gsl1,
                                                          const size_t hsi0, const size_t hsi1,
                                                          const size_t hsk, const size_t ncols,
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
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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

    const auto *gsl0_135 = buffer.data(gsl0 + 135);
    const auto *gsl0_138 = buffer.data(gsl0 + 138);
    const auto *gsl0_141 = buffer.data(gsl0 + 141);
    const auto *gsl0_145 = buffer.data(gsl0 + 145);
    const auto *gsl0_147 = buffer.data(gsl0 + 147);
    const auto *gsl0_150 = buffer.data(gsl0 + 150);
    const auto *gsl0_152 = buffer.data(gsl0 + 152);
    const auto *gsl0_153 = buffer.data(gsl0 + 153);
    const auto *gsl0_156 = buffer.data(gsl0 + 156);
    const auto *gsl0_158 = buffer.data(gsl0 + 158);
    const auto *gsl0_159 = buffer.data(gsl0 + 159);
    const auto *gsl0_160 = buffer.data(gsl0 + 160);
    const auto *gsl0_171 = buffer.data(gsl0 + 171);
    const auto *gsl0_225 = buffer.data(gsl0 + 225);

    const auto *gsk_107 = buffer.data(gsk + 107);
    const auto *gsk_108 = buffer.data(gsk + 108);
    const auto *gsk_111 = buffer.data(gsk + 111);
    const auto *gsk_113 = buffer.data(gsk + 113);
    const auto *gsk_114 = buffer.data(gsk + 114);
    const auto *gsk_115 = buffer.data(gsk + 115);
    const auto *gsk_117 = buffer.data(gsk + 117);
    const auto *gsk_118 = buffer.data(gsk + 118);
    const auto *gsk_119 = buffer.data(gsk + 119);
    const auto *gsk_120 = buffer.data(gsk + 120);
    const auto *gsk_122 = buffer.data(gsk + 122);
    const auto *gsk_123 = buffer.data(gsk + 123);
    const auto *gsk_124 = buffer.data(gsk + 124);
    const auto *gsk_125 = buffer.data(gsk + 125);
    const auto *gsk_126 = buffer.data(gsk + 126);
    const auto *gsk_128 = buffer.data(gsk + 128);
    const auto *gsk_136 = buffer.data(gsk + 136);
    const auto *gsk_143 = buffer.data(gsk + 143);
    const auto *gsk_144 = buffer.data(gsk + 144);
    const auto *gsk_146 = buffer.data(gsk + 146);
    const auto *gsk_149 = buffer.data(gsk + 149);
    const auto *gsk_153 = buffer.data(gsk + 153);
    const auto *gsk_158 = buffer.data(gsk + 158);
    const auto *gsk_164 = buffer.data(gsk + 164);
    const auto *gsk_174 = buffer.data(gsk + 174);
    const auto *gsk_175 = buffer.data(gsk + 175);
    const auto *gsk_176 = buffer.data(gsk + 176);
    const auto *gsk_177 = buffer.data(gsk + 177);
    const auto *gsk_178 = buffer.data(gsk + 178);
    const auto *gsk_179 = buffer.data(gsk + 179);
    const auto *gsk_180 = buffer.data(gsk + 180);
    const auto *gsk_207 = buffer.data(gsk + 207);
    const auto *gsk_208 = buffer.data(gsk + 208);
    const auto *gsk_209 = buffer.data(gsk + 209);
    const auto *gsk_210 = buffer.data(gsk + 210);
    const auto *gsk_211 = buffer.data(gsk + 211);
    const auto *gsk_212 = buffer.data(gsk + 212);
    const auto *gsk_213 = buffer.data(gsk + 213);
    const auto *gsk_215 = buffer.data(gsk + 215);
    const auto *gsk_216 = buffer.data(gsk + 216);
    const auto *gsk_219 = buffer.data(gsk + 219);
    const auto *gsk_222 = buffer.data(gsk + 222);
    const auto *gsk_226 = buffer.data(gsk + 226);
    const auto *gsk_231 = buffer.data(gsk + 231);
    const auto *gsk_237 = buffer.data(gsk + 237);
    const auto *gsk_244 = buffer.data(gsk + 244);
    const auto *gsk_246 = buffer.data(gsk + 246);
    const auto *gsk_247 = buffer.data(gsk + 247);
    const auto *gsk_248 = buffer.data(gsk + 248);
    const auto *gsk_249 = buffer.data(gsk + 249);
    const auto *gsk_250 = buffer.data(gsk + 250);
    const auto *gsk_251 = buffer.data(gsk + 251);
    const auto *gsk_257 = buffer.data(gsk + 257);
    const auto *gsk_261 = buffer.data(gsk + 261);
    const auto *gsk_266 = buffer.data(gsk + 266);
    const auto *gsk_272 = buffer.data(gsk + 272);
    const auto *gsk_279 = buffer.data(gsk + 279);
    const auto *gsk_280 = buffer.data(gsk + 280);
    const auto *gsk_281 = buffer.data(gsk + 281);
    const auto *gsk_282 = buffer.data(gsk + 282);
    const auto *gsk_283 = buffer.data(gsk + 283);
    const auto *gsk_284 = buffer.data(gsk + 284);
    const auto *gsk_285 = buffer.data(gsk + 285);
    const auto *gsk_286 = buffer.data(gsk + 286);
    const auto *gsk_287 = buffer.data(gsk + 287);

    const auto *gsl1_135 = buffer.data(gsl1 + 135);
    const auto *gsl1_138 = buffer.data(gsl1 + 138);
    const auto *gsl1_141 = buffer.data(gsl1 + 141);
    const auto *gsl1_145 = buffer.data(gsl1 + 145);
    const auto *gsl1_147 = buffer.data(gsl1 + 147);
    const auto *gsl1_150 = buffer.data(gsl1 + 150);
    const auto *gsl1_152 = buffer.data(gsl1 + 152);
    const auto *gsl1_153 = buffer.data(gsl1 + 153);
    const auto *gsl1_156 = buffer.data(gsl1 + 156);
    const auto *gsl1_158 = buffer.data(gsl1 + 158);
    const auto *gsl1_159 = buffer.data(gsl1 + 159);
    const auto *gsl1_160 = buffer.data(gsl1 + 160);
    const auto *gsl1_171 = buffer.data(gsl1 + 171);
    const auto *gsl1_225 = buffer.data(gsl1 + 225);

    const auto *hsi0_150 = buffer.data(hsi0 + 150);
    const auto *hsi0_151 = buffer.data(hsi0 + 151);
    const auto *hsi0_152 = buffer.data(hsi0 + 152);
    const auto *hsi0_153 = buffer.data(hsi0 + 153);
    const auto *hsi0_154 = buffer.data(hsi0 + 154);
    const auto *hsi0_161 = buffer.data(hsi0 + 161);
    const auto *hsi0_162 = buffer.data(hsi0 + 162);
    const auto *hsi0_163 = buffer.data(hsi0 + 163);
    const auto *hsi0_164 = buffer.data(hsi0 + 164);
    const auto *hsi0_165 = buffer.data(hsi0 + 165);
    const auto *hsi0_166 = buffer.data(hsi0 + 166);
    const auto *hsi0_167 = buffer.data(hsi0 + 167);
    const auto *hsi0_168 = buffer.data(hsi0 + 168);
    const auto *hsi0_170 = buffer.data(hsi0 + 170);
    const auto *hsi0_171 = buffer.data(hsi0 + 171);
    const auto *hsi0_173 = buffer.data(hsi0 + 173);
    const auto *hsi0_174 = buffer.data(hsi0 + 174);
    const auto *hsi0_175 = buffer.data(hsi0 + 175);
    const auto *hsi0_177 = buffer.data(hsi0 + 177);
    const auto *hsi0_178 = buffer.data(hsi0 + 178);
    const auto *hsi0_179 = buffer.data(hsi0 + 179);
    const auto *hsi0_180 = buffer.data(hsi0 + 180);
    const auto *hsi0_182 = buffer.data(hsi0 + 182);
    const auto *hsi0_183 = buffer.data(hsi0 + 183);
    const auto *hsi0_189 = buffer.data(hsi0 + 189);
    const auto *hsi0_190 = buffer.data(hsi0 + 190);
    const auto *hsi0_191 = buffer.data(hsi0 + 191);
    const auto *hsi0_192 = buffer.data(hsi0 + 192);
    const auto *hsi0_193 = buffer.data(hsi0 + 193);
    const auto *hsi0_195 = buffer.data(hsi0 + 195);
    const auto *hsi0_201 = buffer.data(hsi0 + 201);
    const auto *hsi0_205 = buffer.data(hsi0 + 205);
    const auto *hsi0_210 = buffer.data(hsi0 + 210);
    const auto *hsi0_216 = buffer.data(hsi0 + 216);
    const auto *hsi0_219 = buffer.data(hsi0 + 219);
    const auto *hsi0_220 = buffer.data(hsi0 + 220);
    const auto *hsi0_221 = buffer.data(hsi0 + 221);
    const auto *hsi0_222 = buffer.data(hsi0 + 222);
    const auto *hsi0_223 = buffer.data(hsi0 + 223);

    const auto *hsi1_150 = buffer.data(hsi1 + 150);
    const auto *hsi1_151 = buffer.data(hsi1 + 151);
    const auto *hsi1_152 = buffer.data(hsi1 + 152);
    const auto *hsi1_153 = buffer.data(hsi1 + 153);
    const auto *hsi1_154 = buffer.data(hsi1 + 154);
    const auto *hsi1_161 = buffer.data(hsi1 + 161);
    const auto *hsi1_162 = buffer.data(hsi1 + 162);
    const auto *hsi1_163 = buffer.data(hsi1 + 163);
    const auto *hsi1_164 = buffer.data(hsi1 + 164);
    const auto *hsi1_165 = buffer.data(hsi1 + 165);
    const auto *hsi1_166 = buffer.data(hsi1 + 166);
    const auto *hsi1_167 = buffer.data(hsi1 + 167);
    const auto *hsi1_168 = buffer.data(hsi1 + 168);
    const auto *hsi1_170 = buffer.data(hsi1 + 170);
    const auto *hsi1_171 = buffer.data(hsi1 + 171);
    const auto *hsi1_173 = buffer.data(hsi1 + 173);
    const auto *hsi1_174 = buffer.data(hsi1 + 174);
    const auto *hsi1_175 = buffer.data(hsi1 + 175);
    const auto *hsi1_177 = buffer.data(hsi1 + 177);
    const auto *hsi1_178 = buffer.data(hsi1 + 178);
    const auto *hsi1_179 = buffer.data(hsi1 + 179);
    const auto *hsi1_180 = buffer.data(hsi1 + 180);
    const auto *hsi1_182 = buffer.data(hsi1 + 182);
    const auto *hsi1_183 = buffer.data(hsi1 + 183);
    const auto *hsi1_189 = buffer.data(hsi1 + 189);
    const auto *hsi1_190 = buffer.data(hsi1 + 190);
    const auto *hsi1_191 = buffer.data(hsi1 + 191);
    const auto *hsi1_192 = buffer.data(hsi1 + 192);
    const auto *hsi1_193 = buffer.data(hsi1 + 193);
    const auto *hsi1_195 = buffer.data(hsi1 + 195);
    const auto *hsi1_201 = buffer.data(hsi1 + 201);
    const auto *hsi1_205 = buffer.data(hsi1 + 205);
    const auto *hsi1_210 = buffer.data(hsi1 + 210);
    const auto *hsi1_216 = buffer.data(hsi1 + 216);
    const auto *hsi1_219 = buffer.data(hsi1 + 219);
    const auto *hsi1_220 = buffer.data(hsi1 + 220);
    const auto *hsi1_221 = buffer.data(hsi1 + 221);
    const auto *hsi1_222 = buffer.data(hsi1 + 222);
    const auto *hsi1_223 = buffer.data(hsi1 + 223);

    const auto *hsk_195 = buffer.data(hsk + 195);
    const auto *hsk_196 = buffer.data(hsk + 196);
    const auto *hsk_197 = buffer.data(hsk + 197);
    const auto *hsk_198 = buffer.data(hsk + 198);
    const auto *hsk_199 = buffer.data(hsk + 199);
    const auto *hsk_200 = buffer.data(hsk + 200);
    const auto *hsk_207 = buffer.data(hsk + 207);
    const auto *hsk_208 = buffer.data(hsk + 208);
    const auto *hsk_209 = buffer.data(hsk + 209);
    const auto *hsk_210 = buffer.data(hsk + 210);
    const auto *hsk_211 = buffer.data(hsk + 211);
    const auto *hsk_212 = buffer.data(hsk + 212);
    const auto *hsk_213 = buffer.data(hsk + 213);
    const auto *hsk_214 = buffer.data(hsk + 214);
    const auto *hsk_215 = buffer.data(hsk + 215);
    const auto *hsk_216 = buffer.data(hsk + 216);
    const auto *hsk_217 = buffer.data(hsk + 217);
    const auto *hsk_218 = buffer.data(hsk + 218);
    const auto *hsk_219 = buffer.data(hsk + 219);
    const auto *hsk_221 = buffer.data(hsk + 221);
    const auto *hsk_222 = buffer.data(hsk + 222);
    const auto *hsk_223 = buffer.data(hsk + 223);
    const auto *hsk_225 = buffer.data(hsk + 225);
    const auto *hsk_226 = buffer.data(hsk + 226);
    const auto *hsk_227 = buffer.data(hsk + 227);
    const auto *hsk_228 = buffer.data(hsk + 228);
    const auto *hsk_230 = buffer.data(hsk + 230);
    const auto *hsk_231 = buffer.data(hsk + 231);
    const auto *hsk_232 = buffer.data(hsk + 232);
    const auto *hsk_233 = buffer.data(hsk + 233);
    const auto *hsk_234 = buffer.data(hsk + 234);
    const auto *hsk_236 = buffer.data(hsk + 236);
    const auto *hsk_237 = buffer.data(hsk + 237);
    const auto *hsk_244 = buffer.data(hsk + 244);
    const auto *hsk_245 = buffer.data(hsk + 245);
    const auto *hsk_246 = buffer.data(hsk + 246);
    const auto *hsk_247 = buffer.data(hsk + 247);
    const auto *hsk_248 = buffer.data(hsk + 248);
    const auto *hsk_249 = buffer.data(hsk + 249);
    const auto *hsk_250 = buffer.data(hsk + 250);
    const auto *hsk_251 = buffer.data(hsk + 251);
    const auto *hsk_252 = buffer.data(hsk + 252);
    const auto *hsk_254 = buffer.data(hsk + 254);
    const auto *hsk_255 = buffer.data(hsk + 255);
    const auto *hsk_257 = buffer.data(hsk + 257);
    const auto *hsk_258 = buffer.data(hsk + 258);
    const auto *hsk_261 = buffer.data(hsk + 261);
    const auto *hsk_262 = buffer.data(hsk + 262);
    const auto *hsk_266 = buffer.data(hsk + 266);
    const auto *hsk_267 = buffer.data(hsk + 267);
    const auto *hsk_272 = buffer.data(hsk + 272);
    const auto *hsk_279 = buffer.data(hsk + 279);
    const auto *hsk_280 = buffer.data(hsk + 280);
    const auto *hsk_281 = buffer.data(hsk + 281);
    const auto *hsk_282 = buffer.data(hsk + 282);
    const auto *hsk_283 = buffer.data(hsk + 283);
    const auto *hsk_284 = buffer.data(hsk + 284);
    const auto *hsk_285 = buffer.data(hsk + 285);
    const auto *hsk_286 = buffer.data(hsk + 286);
    const auto *hsk_287 = buffer.data(hsk + 287);
    const auto *hsk_288 = buffer.data(hsk + 288);

#pragma omp simd aligned(t_246, t_247, t_248, pc_y, hsi0_150, hsi0_151, hsi0_152, hsi1_150, \
                         hsi1_151, hsi1_152, hsk_195, hsk_196, \
                         hsk_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * hsi0_150[k]
                   - f_13 * hsi1_150[k]
                   + f_3 * pc_y[k] * hsk_195[k];

        t_247[k] = f_10 * hsi0_151[k]
                   - f_11 * hsi1_151[k]
                   + f_3 * pc_y[k] * hsk_196[k];

        t_248[k] = f_8 * hsi0_152[k]
                   - f_9 * hsi1_152[k]
                   + f_3 * pc_y[k] * hsk_197[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, hsi0_153, hsi0_154, hsi1_153, hsi1_154, \
                         hsk_198, hsk_199, hsk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * hsi0_153[k]
                   - f_7 * hsi1_153[k]
                   + f_3 * pc_y[k] * hsk_198[k];

        t_250[k] = f_4 * hsi0_154[k]
                   - f_5 * hsi1_154[k]
                   + f_3 * pc_y[k] * hsk_199[k];

        t_251[k] = f_3 * pc_y[k] * hsk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, gsk_207, gsk_208, gsk_209, gsk_210, \
                         hsi0_167, hsi1_167, hsk_207, hsk_208, hsk_209, \
                         hsk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_17 * gsk_207[k]
                   + f_4 * hsi0_167[k]
                   - f_5 * hsi1_167[k]
                   + f_3 * pc_x[k] * hsk_207[k];

        t_253[k] = f_17 * gsk_208[k]
                   + f_3 * pc_x[k] * hsk_208[k];

        t_254[k] = f_17 * gsk_209[k]
                   + f_3 * pc_x[k] * hsk_209[k];

        t_255[k] = f_17 * gsk_210[k]
                   + f_3 * pc_x[k] * hsk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pc_x, pc_y, gsk_211, gsk_212, \
                         gsk_213, gsk_215, hsk_207, hsk_211, hsk_212, hsk_213, \
                         hsk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_17 * gsk_211[k]
                   + f_3 * pc_x[k] * hsk_211[k];

        t_257[k] = f_17 * gsk_212[k]
                   + f_3 * pc_x[k] * hsk_212[k];

        t_258[k] = f_17 * gsk_213[k]
                   + f_3 * pc_x[k] * hsk_213[k];

        t_259[k] = f_3 * pc_y[k] * hsk_207[k];

        t_260[k] = f_17 * gsk_215[k]
                   + f_3 * pc_x[k] * hsk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_y, hsi0_161, hsi0_162, hsi0_163, hsi1_161, \
                         hsi1_162, hsi1_163, hsk_208, hsk_209, \
                         hsk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * hsi0_161[k]
                   - f_2 * hsi1_161[k]
                   + f_3 * pc_y[k] * hsk_208[k];

        t_262[k] = f_20 * hsi0_162[k]
                   - f_21 * hsi1_162[k]
                   + f_3 * pc_y[k] * hsk_209[k];

        t_263[k] = f_12 * hsi0_163[k]
                   - f_13 * hsi1_163[k]
                   + f_3 * pc_y[k] * hsk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, hsi0_164, hsi0_165, hsi0_166, hsi1_164, \
                         hsi1_165, hsi1_166, hsk_211, hsk_212, \
                         hsk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * hsi0_164[k]
                   - f_11 * hsi1_164[k]
                   + f_3 * pc_y[k] * hsk_211[k];

        t_265[k] = f_8 * hsi0_165[k]
                   - f_9 * hsi1_165[k]
                   + f_3 * pc_y[k] * hsk_212[k];

        t_266[k] = f_6 * hsi0_166[k]
                   - f_7 * hsi1_166[k]
                   + f_3 * pc_y[k] * hsk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, gsk_107, gsk_216, \
                         hsi0_167, hsi0_168, hsi1_167, hsi1_168, hsk_214, hsk_215, \
                         hsk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_4 * hsi0_167[k]
                   - f_5 * hsi1_167[k]
                   + f_3 * pc_y[k] * hsk_214[k];

        t_268[k] = f_3 * pc_y[k] * hsk_215[k];

        t_269[k] = f_16 * gsk_107[k]
                   + f_1 * hsi0_167[k]
                   - f_2 * hsi1_167[k]
                   + f_3 * pc_z[k] * hsk_215[k];

        t_270[k] = f_16 * gsk_216[k]
                   + f_1 * hsi0_168[k]
                   - f_2 * hsi1_168[k]
                   + f_3 * pc_x[k] * hsk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, gsk_108, gsk_219, \
                         hsi0_171, hsi1_171, hsk_216, hsk_217, \
                         hsk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * gsk_108[k]
                   + f_3 * pc_y[k] * hsk_216[k];

        t_272[k] = f_3 * pc_z[k] * hsk_216[k];

        t_273[k] = f_16 * gsk_219[k]
                   + f_12 * hsi0_171[k]
                   - f_13 * hsi1_171[k]
                   + f_3 * pc_x[k] * hsk_219[k];

        t_274[k] = f_3 * pc_z[k] * hsk_217[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, gsk_222, hsi0_168, hsi0_174, \
                         hsi1_168, hsi1_174, hsk_218, hsk_219, \
                         hsk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_4 * hsi0_168[k]
                   - f_5 * hsi1_168[k]
                   + f_3 * pc_z[k] * hsk_218[k];

        t_276[k] = f_16 * gsk_222[k]
                   + f_10 * hsi0_174[k]
                   - f_11 * hsi1_174[k]
                   + f_3 * pc_x[k] * hsk_222[k];

        t_277[k] = f_3 * pc_z[k] * hsk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pc_x, pc_y, pc_z, gsk_113, gsk_226, \
                         hsi0_170, hsi0_178, hsi1_170, hsi1_178, hsk_221, hsk_222, \
                         hsk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * gsk_113[k]
                   + f_3 * pc_y[k] * hsk_221[k];

        t_279[k] = f_6 * hsi0_170[k]
                   - f_7 * hsi1_170[k]
                   + f_3 * pc_z[k] * hsk_221[k];

        t_280[k] = f_16 * gsk_226[k]
                   + f_8 * hsi0_178[k]
                   - f_9 * hsi1_178[k]
                   + f_3 * pc_x[k] * hsk_226[k];

        t_281[k] = f_3 * pc_z[k] * hsk_222[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pc_y, pc_z, gsk_117, hsi0_171, hsi0_173, \
                         hsi1_171, hsi1_173, hsk_223, hsk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_4 * hsi0_171[k]
                   - f_5 * hsi1_171[k]
                   + f_3 * pc_z[k] * hsk_223[k];

        t_283[k] = f_17 * gsk_117[k]
                   + f_3 * pc_y[k] * hsk_225[k];

        t_284[k] = f_8 * hsi0_173[k]
                   - f_9 * hsi1_173[k]
                   + f_3 * pc_z[k] * hsk_225[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pc_x, pc_z, gsk_231, hsi0_174, hsi0_183, \
                         hsi1_174, hsi1_183, hsk_226, hsk_227, \
                         hsk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_16 * gsk_231[k]
                   + f_6 * hsi0_183[k]
                   - f_7 * hsi1_183[k]
                   + f_3 * pc_x[k] * hsk_231[k];

        t_286[k] = f_3 * pc_z[k] * hsk_226[k];

        t_287[k] = f_4 * hsi0_174[k]
                   - f_5 * hsi1_174[k]
                   + f_3 * pc_z[k] * hsk_227[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, gsk_122, hsi0_175, hsi0_177, \
                         hsi1_175, hsi1_177, hsk_228, hsk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_6 * hsi0_175[k]
                   - f_7 * hsi1_175[k]
                   + f_3 * pc_z[k] * hsk_228[k];

        t_289[k] = f_17 * gsk_122[k]
                   + f_3 * pc_y[k] * hsk_230[k];

        t_290[k] = f_10 * hsi0_177[k]
                   - f_11 * hsi1_177[k]
                   + f_3 * pc_z[k] * hsk_230[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pc_x, pc_z, gsk_237, hsi0_178, hsi0_189, \
                         hsi1_178, hsi1_189, hsk_231, hsk_232, \
                         hsk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_16 * gsk_237[k]
                   + f_4 * hsi0_189[k]
                   - f_5 * hsi1_189[k]
                   + f_3 * pc_x[k] * hsk_237[k];

        t_292[k] = f_3 * pc_z[k] * hsk_231[k];

        t_293[k] = f_4 * hsi0_178[k]
                   - f_5 * hsi1_178[k]
                   + f_3 * pc_z[k] * hsk_232[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_y, pc_z, gsk_128, hsi0_179, hsi0_180, \
                         hsi0_182, hsi1_179, hsi1_180, hsi1_182, hsk_233, hsk_234, \
                         hsk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * hsi0_179[k]
                   - f_7 * hsi1_179[k]
                   + f_3 * pc_z[k] * hsk_233[k];

        t_295[k] = f_8 * hsi0_180[k]
                   - f_9 * hsi1_180[k]
                   + f_3 * pc_z[k] * hsk_234[k];

        t_296[k] = f_17 * gsk_128[k]
                   + f_3 * pc_y[k] * hsk_236[k];

        t_297[k] = f_12 * hsi0_182[k]
                   - f_13 * hsi1_182[k]
                   + f_3 * pc_z[k] * hsk_236[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pc_x, pc_z, gsk_244, gsk_246, \
                         gsk_247, gsk_248, hsk_237, hsk_244, hsk_246, hsk_247, \
                         hsk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_16 * gsk_244[k]
                   + f_3 * pc_x[k] * hsk_244[k];

        t_299[k] = f_3 * pc_z[k] * hsk_237[k];

        t_300[k] = f_16 * gsk_246[k]
                   + f_3 * pc_x[k] * hsk_246[k];

        t_301[k] = f_16 * gsk_247[k]
                   + f_3 * pc_x[k] * hsk_247[k];

        t_302[k] = f_16 * gsk_248[k]
                   + f_3 * pc_x[k] * hsk_248[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, pc_y, gsk_136, gsk_249, gsk_250, \
                         gsk_251, hsi0_189, hsi1_189, hsk_244, hsk_249, hsk_250, \
                         hsk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_16 * gsk_249[k]
                   + f_3 * pc_x[k] * hsk_249[k];

        t_304[k] = f_16 * gsk_250[k]
                   + f_3 * pc_x[k] * hsk_250[k];

        t_305[k] = f_16 * gsk_251[k]
                   + f_3 * pc_x[k] * hsk_251[k];

        t_306[k] = f_17 * gsk_136[k]
                   + f_1 * hsi0_189[k]
                   - f_2 * hsi1_189[k]
                   + f_3 * pc_y[k] * hsk_244[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_z, hsi0_189, hsi0_190, hsi0_191, \
                         hsi1_189, hsi1_190, hsi1_191, hsk_244, hsk_245, hsk_246, \
                         hsk_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_z[k] * hsk_244[k];

        t_308[k] = f_4 * hsi0_189[k]
                   - f_5 * hsi1_189[k]
                   + f_3 * pc_z[k] * hsk_245[k];

        t_309[k] = f_6 * hsi0_190[k]
                   - f_7 * hsi1_190[k]
                   + f_3 * pc_z[k] * hsk_246[k];

        t_310[k] = f_8 * hsi0_191[k]
                   - f_9 * hsi1_191[k]
                   + f_3 * pc_z[k] * hsk_247[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, gsk_143, hsi0_192, hsi0_193, \
                         hsi0_195, hsi1_192, hsi1_193, hsi1_195, hsk_248, hsk_249, \
                         hsk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_10 * hsi0_192[k]
                   - f_11 * hsi1_192[k]
                   + f_3 * pc_z[k] * hsk_248[k];

        t_312[k] = f_12 * hsi0_193[k]
                   - f_13 * hsi1_193[k]
                   + f_3 * pc_z[k] * hsk_249[k];

        t_313[k] = f_17 * gsk_143[k]
                   + f_3 * pc_y[k] * hsk_251[k];

        t_314[k] = f_1 * hsi0_195[k]
                   - f_2 * hsi1_195[k]
                   + f_3 * pc_z[k] * hsk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_z, pc_y, pc_z, gsl0_135, gsl0_138, \
                         gsk_108, gsk_144, gsl1_135, gsl1_138, \
                         hsk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * gsl0_135[k]
                   - f_14 * pc_z[k] * gsl1_135[k];

        t_316[k] = f_16 * gsk_144[k]
                   + f_3 * pc_y[k] * hsk_252[k];

        t_317[k] = f_15 * gsk_108[k]
                   + f_3 * pc_z[k] * hsk_252[k];

        t_318[k] = pa_z[k] * gsl0_138[k]
                   - f_14 * pc_z[k] * gsl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_x, pc_y, pc_z, gsl0_141, gsk_146, \
                         gsk_257, gsl1_141, hsi0_201, hsi1_201, hsk_254, \
                         hsk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * gsk_146[k]
                   + f_3 * pc_y[k] * hsk_254[k];

        t_320[k] = f_16 * gsk_257[k]
                   + f_12 * hsi0_201[k]
                   - f_13 * hsi1_201[k]
                   + f_3 * pc_x[k] * hsk_257[k];

        t_321[k] = pa_z[k] * gsl0_141[k]
                   - f_14 * pc_z[k] * gsl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, gsk_111, gsk_149, gsk_261, \
                         hsi0_205, hsi1_205, hsk_255, hsk_257, \
                         hsk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * gsk_111[k]
                   + f_3 * pc_z[k] * hsk_255[k];

        t_323[k] = f_16 * gsk_149[k]
                   + f_3 * pc_y[k] * hsk_257[k];

        t_324[k] = f_16 * gsk_261[k]
                   + f_10 * hsi0_205[k]
                   - f_11 * hsi1_205[k]
                   + f_3 * pc_x[k] * hsk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pa_z, pc_y, pc_z, gsl0_145, gsl0_147, \
                         gsk_114, gsk_115, gsk_153, gsl1_145, gsl1_147, hsk_258, \
                         hsk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_z[k] * gsl0_145[k]
                   - f_14 * pc_z[k] * gsl1_145[k];

        t_326[k] = f_15 * gsk_114[k]
                   + f_3 * pc_z[k] * hsk_258[k];

        t_327[k] = pa_z[k] * gsl0_147[k]
                   + f_16 * gsk_115[k]
                   - f_14 * pc_z[k] * gsl1_147[k];

        t_328[k] = f_16 * gsk_153[k]
                   + f_3 * pc_y[k] * hsk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_z, pc_x, pc_z, gsl0_150, gsk_118, gsk_266, \
                         gsl1_150, hsi0_210, hsi1_210, hsk_262, \
                         hsk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_16 * gsk_266[k]
                   + f_8 * hsi0_210[k]
                   - f_9 * hsi1_210[k]
                   + f_3 * pc_x[k] * hsk_266[k];

        t_330[k] = pa_z[k] * gsl0_150[k]
                   - f_14 * pc_z[k] * gsl1_150[k];

        t_331[k] = f_15 * gsk_118[k]
                   + f_3 * pc_z[k] * hsk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pa_z, pc_y, pc_z, gsl0_152, gsl0_153, gsk_119, \
                         gsk_120, gsk_158, gsl1_152, gsl1_153, \
                         hsk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pa_z[k] * gsl0_152[k]
                   + f_16 * gsk_119[k]
                   - f_14 * pc_z[k] * gsl1_152[k];

        t_333[k] = pa_z[k] * gsl0_153[k]
                   + f_17 * gsk_120[k]
                   - f_14 * pc_z[k] * gsl1_153[k];

        t_334[k] = f_16 * gsk_158[k]
                   + f_3 * pc_y[k] * hsk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pa_z, pc_x, pc_z, gsl0_156, gsk_123, gsk_272, \
                         gsl1_156, hsi0_216, hsi1_216, hsk_267, \
                         hsk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_16 * gsk_272[k]
                   + f_6 * hsi0_216[k]
                   - f_7 * hsi1_216[k]
                   + f_3 * pc_x[k] * hsk_272[k];

        t_336[k] = pa_z[k] * gsl0_156[k]
                   - f_14 * pc_z[k] * gsl1_156[k];

        t_337[k] = f_15 * gsk_123[k]
                   + f_3 * pc_z[k] * hsk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pa_z, pc_z, gsl0_158, gsl0_159, gsl0_160, \
                         gsk_124, gsk_125, gsk_126, gsl1_158, gsl1_159, \
                         gsl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pa_z[k] * gsl0_158[k]
                   + f_16 * gsk_124[k]
                   - f_14 * pc_z[k] * gsl1_158[k];

        t_339[k] = pa_z[k] * gsl0_159[k]
                   + f_17 * gsk_125[k]
                   - f_14 * pc_z[k] * gsl1_159[k];

        t_340[k] = pa_z[k] * gsl0_160[k]
                   + f_18 * gsk_126[k]
                   - f_14 * pc_z[k] * gsl1_160[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, gsk_164, gsk_279, gsk_280, \
                         gsk_281, hsi0_223, hsi1_223, hsk_272, hsk_279, hsk_280, \
                         hsk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * gsk_164[k]
                   + f_3 * pc_y[k] * hsk_272[k];

        t_342[k] = f_16 * gsk_279[k]
                   + f_4 * hsi0_223[k]
                   - f_5 * hsi1_223[k]
                   + f_3 * pc_x[k] * hsk_279[k];

        t_343[k] = f_16 * gsk_280[k]
                   + f_3 * pc_x[k] * hsk_280[k];

        t_344[k] = f_16 * gsk_281[k]
                   + f_3 * pc_x[k] * hsk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, gsk_282, gsk_283, gsk_284, \
                         gsk_285, gsk_286, hsk_282, hsk_283, hsk_284, hsk_285, \
                         hsk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_16 * gsk_282[k]
                   + f_3 * pc_x[k] * hsk_282[k];

        t_346[k] = f_16 * gsk_283[k]
                   + f_3 * pc_x[k] * hsk_283[k];

        t_347[k] = f_16 * gsk_284[k]
                   + f_3 * pc_x[k] * hsk_284[k];

        t_348[k] = f_16 * gsk_285[k]
                   + f_3 * pc_x[k] * hsk_285[k];

        t_349[k] = f_16 * gsk_286[k]
                   + f_3 * pc_x[k] * hsk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pa_z, pc_x, pc_z, gsl0_171, gsk_136, gsk_287, \
                         gsl1_171, hsk_280, hsk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_16 * gsk_287[k]
                   + f_3 * pc_x[k] * hsk_287[k];

        t_351[k] = pa_z[k] * gsl0_171[k]
                   - f_14 * pc_z[k] * gsl1_171[k];

        t_352[k] = f_15 * gsk_136[k]
                   + f_3 * pc_z[k] * hsk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, gsk_174, gsk_175, gsk_176, hsi0_219, \
                         hsi0_220, hsi0_221, hsi1_219, hsi1_220, hsi1_221, hsk_282, hsk_283, \
                         hsk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * gsk_174[k]
                   + f_12 * hsi0_219[k]
                   - f_13 * hsi1_219[k]
                   + f_3 * pc_y[k] * hsk_282[k];

        t_354[k] = f_16 * gsk_175[k]
                   + f_10 * hsi0_220[k]
                   - f_11 * hsi1_220[k]
                   + f_3 * pc_y[k] * hsk_283[k];

        t_355[k] = f_16 * gsk_176[k]
                   + f_8 * hsi0_221[k]
                   - f_9 * hsi1_221[k]
                   + f_3 * pc_y[k] * hsk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, gsk_177, gsk_178, gsk_179, hsi0_222, \
                         hsi0_223, hsi1_222, hsi1_223, hsk_285, hsk_286, \
                         hsk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * gsk_177[k]
                   + f_6 * hsi0_222[k]
                   - f_7 * hsi1_222[k]
                   + f_3 * pc_y[k] * hsk_285[k];

        t_357[k] = f_16 * gsk_178[k]
                   + f_4 * hsi0_223[k]
                   - f_5 * hsi1_223[k]
                   + f_3 * pc_y[k] * hsk_286[k];

        t_358[k] = f_16 * gsk_179[k]
                   + f_3 * pc_y[k] * hsk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_y, pc_y, pc_z, gsl0_225, gsk_143, \
                         gsk_144, gsk_180, gsl1_225, hsi0_223, hsi1_223, hsk_287, \
                         hsk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * gsk_143[k]
                   + f_1 * hsi0_223[k]
                   - f_2 * hsi1_223[k]
                   + f_3 * pc_z[k] * hsk_287[k];

        t_360[k] = pa_y[k] * gsl0_225[k]
                   - f_14 * pc_y[k] * gsl1_225[k];

        t_361[k] = f_15 * gsk_180[k]
                   + f_3 * pc_y[k] * hsk_288[k];

        t_362[k] = f_16 * gsk_144[k]
                   + f_3 * pc_z[k] * hsk_288[k];
    }
}

static auto
compute_prim_hsl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsl0,
                                                          const size_t gsk, const size_t gsl1,
                                                          const size_t hsi0, const size_t hsi1,
                                                          const size_t hsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsl0_228 = buffer.data(gsl0 + 228);
    const auto *gsl0_230 = buffer.data(gsl0 + 230);
    const auto *gsl0_231 = buffer.data(gsl0 + 231);
    const auto *gsl0_234 = buffer.data(gsl0 + 234);
    const auto *gsl0_235 = buffer.data(gsl0 + 235);
    const auto *gsl0_237 = buffer.data(gsl0 + 237);
    const auto *gsl0_239 = buffer.data(gsl0 + 239);
    const auto *gsl0_240 = buffer.data(gsl0 + 240);
    const auto *gsl0_242 = buffer.data(gsl0 + 242);
    const auto *gsl0_243 = buffer.data(gsl0 + 243);
    const auto *gsl0_245 = buffer.data(gsl0 + 245);
    const auto *gsl0_246 = buffer.data(gsl0 + 246);
    const auto *gsl0_248 = buffer.data(gsl0 + 248);
    const auto *gsl0_249 = buffer.data(gsl0 + 249);
    const auto *gsl0_250 = buffer.data(gsl0 + 250);
    const auto *gsl0_252 = buffer.data(gsl0 + 252);
    const auto *gsl0_269 = buffer.data(gsl0 + 269);
    const auto *gsl0_450 = buffer.data(gsl0 + 450);
    const auto *gsl0_453 = buffer.data(gsl0 + 453);
    const auto *gsl0_456 = buffer.data(gsl0 + 456);
    const auto *gsl0_460 = buffer.data(gsl0 + 460);
    const auto *gsl0_465 = buffer.data(gsl0 + 465);
    const auto *gsl0_471 = buffer.data(gsl0 + 471);

    const auto *gsk_147 = buffer.data(gsk + 147);
    const auto *gsk_150 = buffer.data(gsk + 150);
    const auto *gsk_154 = buffer.data(gsk + 154);
    const auto *gsk_159 = buffer.data(gsk + 159);
    const auto *gsk_172 = buffer.data(gsk + 172);
    const auto *gsk_180 = buffer.data(gsk + 180);
    const auto *gsk_181 = buffer.data(gsk + 181);
    const auto *gsk_182 = buffer.data(gsk + 182);
    const auto *gsk_183 = buffer.data(gsk + 183);
    const auto *gsk_185 = buffer.data(gsk + 185);
    const auto *gsk_186 = buffer.data(gsk + 186);
    const auto *gsk_188 = buffer.data(gsk + 188);
    const auto *gsk_189 = buffer.data(gsk + 189);
    const auto *gsk_190 = buffer.data(gsk + 190);
    const auto *gsk_192 = buffer.data(gsk + 192);
    const auto *gsk_193 = buffer.data(gsk + 193);
    const auto *gsk_194 = buffer.data(gsk + 194);
    const auto *gsk_195 = buffer.data(gsk + 195);
    const auto *gsk_197 = buffer.data(gsk + 197);
    const auto *gsk_198 = buffer.data(gsk + 198);
    const auto *gsk_199 = buffer.data(gsk + 199);
    const auto *gsk_200 = buffer.data(gsk + 200);
    const auto *gsk_208 = buffer.data(gsk + 208);
    const auto *gsk_210 = buffer.data(gsk + 210);
    const auto *gsk_211 = buffer.data(gsk + 211);
    const auto *gsk_212 = buffer.data(gsk + 212);
    const auto *gsk_213 = buffer.data(gsk + 213);
    const auto *gsk_214 = buffer.data(gsk + 214);
    const auto *gsk_215 = buffer.data(gsk + 215);
    const auto *gsk_216 = buffer.data(gsk + 216);
    const auto *gsk_221 = buffer.data(gsk + 221);
    const auto *gsk_225 = buffer.data(gsk + 225);
    const auto *gsk_230 = buffer.data(gsk + 230);
    const auto *gsk_236 = buffer.data(gsk + 236);
    const auto *gsk_316 = buffer.data(gsk + 316);
    const auto *gsk_317 = buffer.data(gsk + 317);
    const auto *gsk_318 = buffer.data(gsk + 318);
    const auto *gsk_319 = buffer.data(gsk + 319);
    const auto *gsk_320 = buffer.data(gsk + 320);
    const auto *gsk_321 = buffer.data(gsk + 321);
    const auto *gsk_322 = buffer.data(gsk + 322);
    const auto *gsk_323 = buffer.data(gsk + 323);
    const auto *gsk_324 = buffer.data(gsk + 324);
    const auto *gsk_329 = buffer.data(gsk + 329);
    const auto *gsk_333 = buffer.data(gsk + 333);
    const auto *gsk_338 = buffer.data(gsk + 338);
    const auto *gsk_344 = buffer.data(gsk + 344);
    const auto *gsk_351 = buffer.data(gsk + 351);
    const auto *gsk_352 = buffer.data(gsk + 352);
    const auto *gsk_353 = buffer.data(gsk + 353);
    const auto *gsk_354 = buffer.data(gsk + 354);
    const auto *gsk_355 = buffer.data(gsk + 355);
    const auto *gsk_356 = buffer.data(gsk + 356);
    const auto *gsk_357 = buffer.data(gsk + 357);
    const auto *gsk_359 = buffer.data(gsk + 359);
    const auto *gsk_360 = buffer.data(gsk + 360);
    const auto *gsk_363 = buffer.data(gsk + 363);
    const auto *gsk_366 = buffer.data(gsk + 366);
    const auto *gsk_370 = buffer.data(gsk + 370);
    const auto *gsk_375 = buffer.data(gsk + 375);
    const auto *gsk_381 = buffer.data(gsk + 381);

    const auto *gsl1_228 = buffer.data(gsl1 + 228);
    const auto *gsl1_230 = buffer.data(gsl1 + 230);
    const auto *gsl1_231 = buffer.data(gsl1 + 231);
    const auto *gsl1_234 = buffer.data(gsl1 + 234);
    const auto *gsl1_235 = buffer.data(gsl1 + 235);
    const auto *gsl1_237 = buffer.data(gsl1 + 237);
    const auto *gsl1_239 = buffer.data(gsl1 + 239);
    const auto *gsl1_240 = buffer.data(gsl1 + 240);
    const auto *gsl1_242 = buffer.data(gsl1 + 242);
    const auto *gsl1_243 = buffer.data(gsl1 + 243);
    const auto *gsl1_245 = buffer.data(gsl1 + 245);
    const auto *gsl1_246 = buffer.data(gsl1 + 246);
    const auto *gsl1_248 = buffer.data(gsl1 + 248);
    const auto *gsl1_249 = buffer.data(gsl1 + 249);
    const auto *gsl1_250 = buffer.data(gsl1 + 250);
    const auto *gsl1_252 = buffer.data(gsl1 + 252);
    const auto *gsl1_269 = buffer.data(gsl1 + 269);
    const auto *gsl1_450 = buffer.data(gsl1 + 450);
    const auto *gsl1_453 = buffer.data(gsl1 + 453);
    const auto *gsl1_456 = buffer.data(gsl1 + 456);
    const auto *gsl1_460 = buffer.data(gsl1 + 460);
    const auto *gsl1_465 = buffer.data(gsl1 + 465);
    const auto *gsl1_471 = buffer.data(gsl1 + 471);

    const auto *hsi0_245 = buffer.data(hsi0 + 245);
    const auto *hsi0_247 = buffer.data(hsi0 + 247);
    const auto *hsi0_248 = buffer.data(hsi0 + 248);
    const auto *hsi0_249 = buffer.data(hsi0 + 249);
    const auto *hsi0_250 = buffer.data(hsi0 + 250);
    const auto *hsi0_251 = buffer.data(hsi0 + 251);
    const auto *hsi0_252 = buffer.data(hsi0 + 252);
    const auto *hsi0_253 = buffer.data(hsi0 + 253);
    const auto *hsi0_254 = buffer.data(hsi0 + 254);
    const auto *hsi0_255 = buffer.data(hsi0 + 255);
    const auto *hsi0_256 = buffer.data(hsi0 + 256);
    const auto *hsi0_257 = buffer.data(hsi0 + 257);
    const auto *hsi0_258 = buffer.data(hsi0 + 258);
    const auto *hsi0_259 = buffer.data(hsi0 + 259);
    const auto *hsi0_260 = buffer.data(hsi0 + 260);
    const auto *hsi0_261 = buffer.data(hsi0 + 261);
    const auto *hsi0_262 = buffer.data(hsi0 + 262);
    const auto *hsi0_263 = buffer.data(hsi0 + 263);
    const auto *hsi0_264 = buffer.data(hsi0 + 264);
    const auto *hsi0_265 = buffer.data(hsi0 + 265);
    const auto *hsi0_266 = buffer.data(hsi0 + 266);
    const auto *hsi0_272 = buffer.data(hsi0 + 272);
    const auto *hsi0_273 = buffer.data(hsi0 + 273);
    const auto *hsi0_274 = buffer.data(hsi0 + 274);
    const auto *hsi0_275 = buffer.data(hsi0 + 275);
    const auto *hsi0_276 = buffer.data(hsi0 + 276);
    const auto *hsi0_277 = buffer.data(hsi0 + 277);
    const auto *hsi0_278 = buffer.data(hsi0 + 278);
    const auto *hsi0_279 = buffer.data(hsi0 + 279);
    const auto *hsi0_280 = buffer.data(hsi0 + 280);
    const auto *hsi0_282 = buffer.data(hsi0 + 282);
    const auto *hsi0_283 = buffer.data(hsi0 + 283);
    const auto *hsi0_285 = buffer.data(hsi0 + 285);
    const auto *hsi0_286 = buffer.data(hsi0 + 286);
    const auto *hsi0_287 = buffer.data(hsi0 + 287);
    const auto *hsi0_289 = buffer.data(hsi0 + 289);
    const auto *hsi0_290 = buffer.data(hsi0 + 290);
    const auto *hsi0_291 = buffer.data(hsi0 + 291);
    const auto *hsi0_292 = buffer.data(hsi0 + 292);
    const auto *hsi0_294 = buffer.data(hsi0 + 294);

    const auto *hsi1_245 = buffer.data(hsi1 + 245);
    const auto *hsi1_247 = buffer.data(hsi1 + 247);
    const auto *hsi1_248 = buffer.data(hsi1 + 248);
    const auto *hsi1_249 = buffer.data(hsi1 + 249);
    const auto *hsi1_250 = buffer.data(hsi1 + 250);
    const auto *hsi1_251 = buffer.data(hsi1 + 251);
    const auto *hsi1_252 = buffer.data(hsi1 + 252);
    const auto *hsi1_253 = buffer.data(hsi1 + 253);
    const auto *hsi1_254 = buffer.data(hsi1 + 254);
    const auto *hsi1_255 = buffer.data(hsi1 + 255);
    const auto *hsi1_256 = buffer.data(hsi1 + 256);
    const auto *hsi1_257 = buffer.data(hsi1 + 257);
    const auto *hsi1_258 = buffer.data(hsi1 + 258);
    const auto *hsi1_259 = buffer.data(hsi1 + 259);
    const auto *hsi1_260 = buffer.data(hsi1 + 260);
    const auto *hsi1_261 = buffer.data(hsi1 + 261);
    const auto *hsi1_262 = buffer.data(hsi1 + 262);
    const auto *hsi1_263 = buffer.data(hsi1 + 263);
    const auto *hsi1_264 = buffer.data(hsi1 + 264);
    const auto *hsi1_265 = buffer.data(hsi1 + 265);
    const auto *hsi1_266 = buffer.data(hsi1 + 266);
    const auto *hsi1_272 = buffer.data(hsi1 + 272);
    const auto *hsi1_273 = buffer.data(hsi1 + 273);
    const auto *hsi1_274 = buffer.data(hsi1 + 274);
    const auto *hsi1_275 = buffer.data(hsi1 + 275);
    const auto *hsi1_276 = buffer.data(hsi1 + 276);
    const auto *hsi1_277 = buffer.data(hsi1 + 277);
    const auto *hsi1_278 = buffer.data(hsi1 + 278);
    const auto *hsi1_279 = buffer.data(hsi1 + 279);
    const auto *hsi1_280 = buffer.data(hsi1 + 280);
    const auto *hsi1_282 = buffer.data(hsi1 + 282);
    const auto *hsi1_283 = buffer.data(hsi1 + 283);
    const auto *hsi1_285 = buffer.data(hsi1 + 285);
    const auto *hsi1_286 = buffer.data(hsi1 + 286);
    const auto *hsi1_287 = buffer.data(hsi1 + 287);
    const auto *hsi1_289 = buffer.data(hsi1 + 289);
    const auto *hsi1_290 = buffer.data(hsi1 + 290);
    const auto *hsi1_291 = buffer.data(hsi1 + 291);
    const auto *hsi1_292 = buffer.data(hsi1 + 292);
    const auto *hsi1_294 = buffer.data(hsi1 + 294);

    const auto *hsk_290 = buffer.data(hsk + 290);
    const auto *hsk_291 = buffer.data(hsk + 291);
    const auto *hsk_293 = buffer.data(hsk + 293);
    const auto *hsk_294 = buffer.data(hsk + 294);
    const auto *hsk_297 = buffer.data(hsk + 297);
    const auto *hsk_298 = buffer.data(hsk + 298);
    const auto *hsk_302 = buffer.data(hsk + 302);
    const auto *hsk_303 = buffer.data(hsk + 303);
    const auto *hsk_308 = buffer.data(hsk + 308);
    const auto *hsk_316 = buffer.data(hsk + 316);
    const auto *hsk_317 = buffer.data(hsk + 317);
    const auto *hsk_318 = buffer.data(hsk + 318);
    const auto *hsk_319 = buffer.data(hsk + 319);
    const auto *hsk_320 = buffer.data(hsk + 320);
    const auto *hsk_321 = buffer.data(hsk + 321);
    const auto *hsk_322 = buffer.data(hsk + 322);
    const auto *hsk_323 = buffer.data(hsk + 323);
    const auto *hsk_324 = buffer.data(hsk + 324);
    const auto *hsk_325 = buffer.data(hsk + 325);
    const auto *hsk_326 = buffer.data(hsk + 326);
    const auto *hsk_327 = buffer.data(hsk + 327);
    const auto *hsk_328 = buffer.data(hsk + 328);
    const auto *hsk_329 = buffer.data(hsk + 329);
    const auto *hsk_330 = buffer.data(hsk + 330);
    const auto *hsk_331 = buffer.data(hsk + 331);
    const auto *hsk_332 = buffer.data(hsk + 332);
    const auto *hsk_333 = buffer.data(hsk + 333);
    const auto *hsk_334 = buffer.data(hsk + 334);
    const auto *hsk_335 = buffer.data(hsk + 335);
    const auto *hsk_336 = buffer.data(hsk + 336);
    const auto *hsk_337 = buffer.data(hsk + 337);
    const auto *hsk_338 = buffer.data(hsk + 338);
    const auto *hsk_339 = buffer.data(hsk + 339);
    const auto *hsk_340 = buffer.data(hsk + 340);
    const auto *hsk_341 = buffer.data(hsk + 341);
    const auto *hsk_342 = buffer.data(hsk + 342);
    const auto *hsk_343 = buffer.data(hsk + 343);
    const auto *hsk_344 = buffer.data(hsk + 344);
    const auto *hsk_351 = buffer.data(hsk + 351);
    const auto *hsk_352 = buffer.data(hsk + 352);
    const auto *hsk_353 = buffer.data(hsk + 353);
    const auto *hsk_354 = buffer.data(hsk + 354);
    const auto *hsk_355 = buffer.data(hsk + 355);
    const auto *hsk_356 = buffer.data(hsk + 356);
    const auto *hsk_357 = buffer.data(hsk + 357);
    const auto *hsk_358 = buffer.data(hsk + 358);
    const auto *hsk_359 = buffer.data(hsk + 359);
    const auto *hsk_360 = buffer.data(hsk + 360);
    const auto *hsk_361 = buffer.data(hsk + 361);
    const auto *hsk_362 = buffer.data(hsk + 362);
    const auto *hsk_363 = buffer.data(hsk + 363);
    const auto *hsk_365 = buffer.data(hsk + 365);
    const auto *hsk_366 = buffer.data(hsk + 366);
    const auto *hsk_367 = buffer.data(hsk + 367);
    const auto *hsk_369 = buffer.data(hsk + 369);
    const auto *hsk_370 = buffer.data(hsk + 370);
    const auto *hsk_371 = buffer.data(hsk + 371);
    const auto *hsk_372 = buffer.data(hsk + 372);
    const auto *hsk_374 = buffer.data(hsk + 374);
    const auto *hsk_375 = buffer.data(hsk + 375);
    const auto *hsk_376 = buffer.data(hsk + 376);
    const auto *hsk_377 = buffer.data(hsk + 377);
    const auto *hsk_378 = buffer.data(hsk + 378);
    const auto *hsk_380 = buffer.data(hsk + 380);

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_y, pc_y, gsl0_228, gsl0_230, gsl0_231, \
                         gsk_181, gsk_182, gsk_183, gsl1_228, gsl1_230, gsl1_231, \
                         hsk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_y[k] * gsl0_228[k]
                   + f_16 * gsk_181[k]
                   - f_14 * pc_y[k] * gsl1_228[k];

        t_364[k] = f_15 * gsk_182[k]
                   + f_3 * pc_y[k] * hsk_290[k];

        t_365[k] = pa_y[k] * gsl0_230[k]
                   - f_14 * pc_y[k] * gsl1_230[k];

        t_366[k] = pa_y[k] * gsl0_231[k]
                   + f_17 * gsk_183[k]
                   - f_14 * pc_y[k] * gsl1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_y, pc_y, pc_z, gsl0_234, gsl0_235, \
                         gsk_147, gsk_185, gsk_186, gsl1_234, gsl1_235, hsk_291, \
                         hsk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * gsk_147[k]
                   + f_3 * pc_z[k] * hsk_291[k];

        t_368[k] = f_15 * gsk_185[k]
                   + f_3 * pc_y[k] * hsk_293[k];

        t_369[k] = pa_y[k] * gsl0_234[k]
                   - f_14 * pc_y[k] * gsl1_234[k];

        t_370[k] = pa_y[k] * gsl0_235[k]
                   + f_18 * gsk_186[k]
                   - f_14 * pc_y[k] * gsl1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pc_y, pc_z, gsl0_237, gsl0_239, \
                         gsk_150, gsk_188, gsk_189, gsl1_237, gsl1_239, hsk_294, \
                         hsk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * gsk_150[k]
                   + f_3 * pc_z[k] * hsk_294[k];

        t_372[k] = pa_y[k] * gsl0_237[k]
                   + f_16 * gsk_188[k]
                   - f_14 * pc_y[k] * gsl1_237[k];

        t_373[k] = f_15 * gsk_189[k]
                   + f_3 * pc_y[k] * hsk_297[k];

        t_374[k] = pa_y[k] * gsl0_239[k]
                   - f_14 * pc_y[k] * gsl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_y, pc_y, pc_z, gsl0_240, gsl0_242, gsk_154, \
                         gsk_190, gsk_192, gsl1_240, gsl1_242, \
                         hsk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pa_y[k] * gsl0_240[k]
                   + f_0 * gsk_190[k]
                   - f_14 * pc_y[k] * gsl1_240[k];

        t_376[k] = f_16 * gsk_154[k]
                   + f_3 * pc_z[k] * hsk_298[k];

        t_377[k] = pa_y[k] * gsl0_242[k]
                   + f_17 * gsk_192[k]
                   - f_14 * pc_y[k] * gsl1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_y, pc_y, gsl0_243, gsl0_245, gsl0_246, \
                         gsk_193, gsk_194, gsk_195, gsl1_243, gsl1_245, gsl1_246, \
                         hsk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * gsl0_243[k]
                   + f_16 * gsk_193[k]
                   - f_14 * pc_y[k] * gsl1_243[k];

        t_379[k] = f_15 * gsk_194[k]
                   + f_3 * pc_y[k] * hsk_302[k];

        t_380[k] = pa_y[k] * gsl0_245[k]
                   - f_14 * pc_y[k] * gsl1_245[k];

        t_381[k] = pa_y[k] * gsl0_246[k]
                   + f_19 * gsk_195[k]
                   - f_14 * pc_y[k] * gsl1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pa_y, pc_y, pc_z, gsl0_248, gsl0_249, gsk_159, \
                         gsk_197, gsk_198, gsl1_248, gsl1_249, \
                         hsk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * gsk_159[k]
                   + f_3 * pc_z[k] * hsk_303[k];

        t_383[k] = pa_y[k] * gsl0_248[k]
                   + f_18 * gsk_197[k]
                   - f_14 * pc_y[k] * gsl1_248[k];

        t_384[k] = pa_y[k] * gsl0_249[k]
                   + f_17 * gsk_198[k]
                   - f_14 * pc_y[k] * gsl1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_y, pc_x, pc_y, gsl0_250, gsl0_252, \
                         gsk_199, gsk_200, gsk_316, gsl1_250, gsl1_252, hsk_308, \
                         hsk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pa_y[k] * gsl0_250[k]
                   + f_16 * gsk_199[k]
                   - f_14 * pc_y[k] * gsl1_250[k];

        t_386[k] = f_15 * gsk_200[k]
                   + f_3 * pc_y[k] * hsk_308[k];

        t_387[k] = pa_y[k] * gsl0_252[k]
                   - f_14 * pc_y[k] * gsl1_252[k];

        t_388[k] = f_16 * gsk_316[k]
                   + f_3 * pc_x[k] * hsk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, gsk_317, gsk_318, gsk_319, \
                         gsk_320, gsk_321, hsk_317, hsk_318, hsk_319, hsk_320, \
                         hsk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_16 * gsk_317[k]
                   + f_3 * pc_x[k] * hsk_317[k];

        t_390[k] = f_16 * gsk_318[k]
                   + f_3 * pc_x[k] * hsk_318[k];

        t_391[k] = f_16 * gsk_319[k]
                   + f_3 * pc_x[k] * hsk_319[k];

        t_392[k] = f_16 * gsk_320[k]
                   + f_3 * pc_x[k] * hsk_320[k];

        t_393[k] = f_16 * gsk_321[k]
                   + f_3 * pc_x[k] * hsk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, gsk_172, gsk_208, \
                         gsk_322, gsk_323, hsi0_245, hsi1_245, hsk_316, hsk_322, \
                         hsk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * gsk_322[k]
                   + f_3 * pc_x[k] * hsk_322[k];

        t_395[k] = f_16 * gsk_323[k]
                   + f_3 * pc_x[k] * hsk_323[k];

        t_396[k] = f_15 * gsk_208[k]
                   + f_1 * hsi0_245[k]
                   - f_2 * hsi1_245[k]
                   + f_3 * pc_y[k] * hsk_316[k];

        t_397[k] = f_16 * gsk_172[k]
                   + f_3 * pc_z[k] * hsk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, gsk_210, gsk_211, gsk_212, hsi0_247, \
                         hsi0_248, hsi0_249, hsi1_247, hsi1_248, hsi1_249, hsk_318, hsk_319, \
                         hsk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * gsk_210[k]
                   + f_12 * hsi0_247[k]
                   - f_13 * hsi1_247[k]
                   + f_3 * pc_y[k] * hsk_318[k];

        t_399[k] = f_15 * gsk_211[k]
                   + f_10 * hsi0_248[k]
                   - f_11 * hsi1_248[k]
                   + f_3 * pc_y[k] * hsk_319[k];

        t_400[k] = f_15 * gsk_212[k]
                   + f_8 * hsi0_249[k]
                   - f_9 * hsi1_249[k]
                   + f_3 * pc_y[k] * hsk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, gsk_213, gsk_214, gsk_215, hsi0_250, \
                         hsi0_251, hsi1_250, hsi1_251, hsk_321, hsk_322, \
                         hsk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * gsk_213[k]
                   + f_6 * hsi0_250[k]
                   - f_7 * hsi1_250[k]
                   + f_3 * pc_y[k] * hsk_321[k];

        t_402[k] = f_15 * gsk_214[k]
                   + f_4 * hsi0_251[k]
                   - f_5 * hsi1_251[k]
                   + f_3 * pc_y[k] * hsk_322[k];

        t_403[k] = f_15 * gsk_215[k]
                   + f_3 * pc_y[k] * hsk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_y, pc_x, pc_y, pc_z, gsl0_269, \
                         gsk_180, gsk_324, gsl1_269, hsi0_252, hsi1_252, \
                         hsk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * gsl0_269[k]
                   - f_14 * pc_y[k] * gsl1_269[k];

        t_405[k] = f_16 * gsk_324[k]
                   + f_1 * hsi0_252[k]
                   - f_2 * hsi1_252[k]
                   + f_3 * pc_x[k] * hsk_324[k];

        t_406[k] = f_3 * pc_y[k] * hsk_324[k];

        t_407[k] = f_17 * gsk_180[k]
                   + f_3 * pc_z[k] * hsk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, gsk_329, hsi0_252, hsi0_257, \
                         hsi1_252, hsi1_257, hsk_325, hsk_326, \
                         hsk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_4 * hsi0_252[k]
                   - f_5 * hsi1_252[k]
                   + f_3 * pc_y[k] * hsk_325[k];

        t_409[k] = f_3 * pc_y[k] * hsk_326[k];

        t_410[k] = f_16 * gsk_329[k]
                   + f_12 * hsi0_257[k]
                   - f_13 * hsi1_257[k]
                   + f_3 * pc_x[k] * hsk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_y, hsi0_253, hsi0_254, hsi1_253, hsi1_254, \
                         hsk_327, hsk_328, hsk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_6 * hsi0_253[k]
                   - f_7 * hsi1_253[k]
                   + f_3 * pc_y[k] * hsk_327[k];

        t_412[k] = f_4 * hsi0_254[k]
                   - f_5 * hsi1_254[k]
                   + f_3 * pc_y[k] * hsk_328[k];

        t_413[k] = f_3 * pc_y[k] * hsk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, gsk_333, hsi0_255, hsi0_256, \
                         hsi0_261, hsi1_255, hsi1_256, hsi1_261, hsk_330, hsk_331, \
                         hsk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_16 * gsk_333[k]
                   + f_10 * hsi0_261[k]
                   - f_11 * hsi1_261[k]
                   + f_3 * pc_x[k] * hsk_333[k];

        t_415[k] = f_8 * hsi0_255[k]
                   - f_9 * hsi1_255[k]
                   + f_3 * pc_y[k] * hsk_330[k];

        t_416[k] = f_6 * hsi0_256[k]
                   - f_7 * hsi1_256[k]
                   + f_3 * pc_y[k] * hsk_331[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, gsk_338, hsi0_257, hsi0_266, \
                         hsi1_257, hsi1_266, hsk_332, hsk_333, \
                         hsk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * hsi0_257[k]
                   - f_5 * hsi1_257[k]
                   + f_3 * pc_y[k] * hsk_332[k];

        t_418[k] = f_3 * pc_y[k] * hsk_333[k];

        t_419[k] = f_16 * gsk_338[k]
                   + f_8 * hsi0_266[k]
                   - f_9 * hsi1_266[k]
                   + f_3 * pc_x[k] * hsk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_y, hsi0_258, hsi0_259, hsi0_260, hsi1_258, \
                         hsi1_259, hsi1_260, hsk_334, hsk_335, \
                         hsk_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_10 * hsi0_258[k]
                   - f_11 * hsi1_258[k]
                   + f_3 * pc_y[k] * hsk_334[k];

        t_421[k] = f_8 * hsi0_259[k]
                   - f_9 * hsi1_259[k]
                   + f_3 * pc_y[k] * hsk_335[k];

        t_422[k] = f_6 * hsi0_260[k]
                   - f_7 * hsi1_260[k]
                   + f_3 * pc_y[k] * hsk_336[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, gsk_344, hsi0_261, hsi0_272, \
                         hsi1_261, hsi1_272, hsk_337, hsk_338, \
                         hsk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * hsi0_261[k]
                   - f_5 * hsi1_261[k]
                   + f_3 * pc_y[k] * hsk_337[k];

        t_424[k] = f_3 * pc_y[k] * hsk_338[k];

        t_425[k] = f_16 * gsk_344[k]
                   + f_6 * hsi0_272[k]
                   - f_7 * hsi1_272[k]
                   + f_3 * pc_x[k] * hsk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_y, hsi0_262, hsi0_263, hsi0_264, hsi1_262, \
                         hsi1_263, hsi1_264, hsk_339, hsk_340, \
                         hsk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_12 * hsi0_262[k]
                   - f_13 * hsi1_262[k]
                   + f_3 * pc_y[k] * hsk_339[k];

        t_427[k] = f_10 * hsi0_263[k]
                   - f_11 * hsi1_263[k]
                   + f_3 * pc_y[k] * hsk_340[k];

        t_428[k] = f_8 * hsi0_264[k]
                   - f_9 * hsi1_264[k]
                   + f_3 * pc_y[k] * hsk_341[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_y, hsi0_265, hsi0_266, hsi1_265, hsi1_266, \
                         hsk_342, hsk_343, hsk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_6 * hsi0_265[k]
                   - f_7 * hsi1_265[k]
                   + f_3 * pc_y[k] * hsk_342[k];

        t_430[k] = f_4 * hsi0_266[k]
                   - f_5 * hsi1_266[k]
                   + f_3 * pc_y[k] * hsk_343[k];

        t_431[k] = f_3 * pc_y[k] * hsk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, gsk_351, gsk_352, gsk_353, gsk_354, \
                         hsi0_279, hsi1_279, hsk_351, hsk_352, hsk_353, \
                         hsk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_16 * gsk_351[k]
                   + f_4 * hsi0_279[k]
                   - f_5 * hsi1_279[k]
                   + f_3 * pc_x[k] * hsk_351[k];

        t_433[k] = f_16 * gsk_352[k]
                   + f_3 * pc_x[k] * hsk_352[k];

        t_434[k] = f_16 * gsk_353[k]
                   + f_3 * pc_x[k] * hsk_353[k];

        t_435[k] = f_16 * gsk_354[k]
                   + f_3 * pc_x[k] * hsk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_y, gsk_355, gsk_356, \
                         gsk_357, gsk_359, hsk_351, hsk_355, hsk_356, hsk_357, \
                         hsk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_16 * gsk_355[k]
                   + f_3 * pc_x[k] * hsk_355[k];

        t_437[k] = f_16 * gsk_356[k]
                   + f_3 * pc_x[k] * hsk_356[k];

        t_438[k] = f_16 * gsk_357[k]
                   + f_3 * pc_x[k] * hsk_357[k];

        t_439[k] = f_3 * pc_y[k] * hsk_351[k];

        t_440[k] = f_16 * gsk_359[k]
                   + f_3 * pc_x[k] * hsk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, hsi0_273, hsi0_274, hsi0_275, hsi1_273, \
                         hsi1_274, hsi1_275, hsk_352, hsk_353, \
                         hsk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * hsi0_273[k]
                   - f_2 * hsi1_273[k]
                   + f_3 * pc_y[k] * hsk_352[k];

        t_442[k] = f_20 * hsi0_274[k]
                   - f_21 * hsi1_274[k]
                   + f_3 * pc_y[k] * hsk_353[k];

        t_443[k] = f_12 * hsi0_275[k]
                   - f_13 * hsi1_275[k]
                   + f_3 * pc_y[k] * hsk_354[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_y, hsi0_276, hsi0_277, hsi0_278, hsi1_276, \
                         hsi1_277, hsi1_278, hsk_355, hsk_356, \
                         hsk_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_10 * hsi0_276[k]
                   - f_11 * hsi1_276[k]
                   + f_3 * pc_y[k] * hsk_355[k];

        t_445[k] = f_8 * hsi0_277[k]
                   - f_9 * hsi1_277[k]
                   + f_3 * pc_y[k] * hsk_356[k];

        t_446[k] = f_6 * hsi0_278[k]
                   - f_7 * hsi1_278[k]
                   + f_3 * pc_y[k] * hsk_357[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pa_x, pc_x, pc_y, pc_z, gsl0_450, \
                         gsk_215, gsk_360, gsl1_450, hsi0_279, hsi1_279, hsk_358, \
                         hsk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_4 * hsi0_279[k]
                   - f_5 * hsi1_279[k]
                   + f_3 * pc_y[k] * hsk_358[k];

        t_448[k] = f_3 * pc_y[k] * hsk_359[k];

        t_449[k] = f_17 * gsk_215[k]
                   + f_1 * hsi0_279[k]
                   - f_2 * hsi1_279[k]
                   + f_3 * pc_z[k] * hsk_359[k];

        t_450[k] = pa_x[k] * gsl0_450[k]
                   + f_22 * gsk_360[k]
                   - f_14 * pc_x[k] * gsl1_450[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pa_x, pc_x, pc_y, pc_z, gsl0_453, \
                         gsk_216, gsk_363, gsl1_453, hsk_360, hsk_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_18 * gsk_216[k]
                   + f_3 * pc_y[k] * hsk_360[k];

        t_452[k] = f_3 * pc_z[k] * hsk_360[k];

        t_453[k] = pa_x[k] * gsl0_453[k]
                   + f_19 * gsk_363[k]
                   - f_14 * pc_x[k] * gsl1_453[k];

        t_454[k] = f_3 * pc_z[k] * hsk_361[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pa_x, pc_x, pc_z, gsl0_456, gsk_366, gsl1_456, \
                         hsi0_280, hsi1_280, hsk_362, hsk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_4 * hsi0_280[k]
                   - f_5 * hsi1_280[k]
                   + f_3 * pc_z[k] * hsk_362[k];

        t_456[k] = pa_x[k] * gsl0_456[k]
                   + f_0 * gsk_366[k]
                   - f_14 * pc_x[k] * gsl1_456[k];

        t_457[k] = f_3 * pc_z[k] * hsk_363[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pa_x, pc_x, pc_y, pc_z, gsl0_460, \
                         gsk_221, gsk_370, gsl1_460, hsi0_282, hsi1_282, hsk_365, \
                         hsk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_18 * gsk_221[k]
                   + f_3 * pc_y[k] * hsk_365[k];

        t_459[k] = f_6 * hsi0_282[k]
                   - f_7 * hsi1_282[k]
                   + f_3 * pc_z[k] * hsk_365[k];

        t_460[k] = pa_x[k] * gsl0_460[k]
                   + f_18 * gsk_370[k]
                   - f_14 * pc_x[k] * gsl1_460[k];

        t_461[k] = f_3 * pc_z[k] * hsk_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, gsk_225, hsi0_283, hsi0_285, \
                         hsi1_283, hsi1_285, hsk_367, hsk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_4 * hsi0_283[k]
                   - f_5 * hsi1_283[k]
                   + f_3 * pc_z[k] * hsk_367[k];

        t_463[k] = f_18 * gsk_225[k]
                   + f_3 * pc_y[k] * hsk_369[k];

        t_464[k] = f_8 * hsi0_285[k]
                   - f_9 * hsi1_285[k]
                   + f_3 * pc_z[k] * hsk_369[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pa_x, pc_x, pc_z, gsl0_465, gsk_375, gsl1_465, \
                         hsi0_286, hsi1_286, hsk_370, hsk_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_x[k] * gsl0_465[k]
                   + f_17 * gsk_375[k]
                   - f_14 * pc_x[k] * gsl1_465[k];

        t_466[k] = f_3 * pc_z[k] * hsk_370[k];

        t_467[k] = f_4 * hsi0_286[k]
                   - f_5 * hsi1_286[k]
                   + f_3 * pc_z[k] * hsk_371[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_y, pc_z, gsk_230, hsi0_287, hsi0_289, \
                         hsi1_287, hsi1_289, hsk_372, hsk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_6 * hsi0_287[k]
                   - f_7 * hsi1_287[k]
                   + f_3 * pc_z[k] * hsk_372[k];

        t_469[k] = f_18 * gsk_230[k]
                   + f_3 * pc_y[k] * hsk_374[k];

        t_470[k] = f_10 * hsi0_289[k]
                   - f_11 * hsi1_289[k]
                   + f_3 * pc_z[k] * hsk_374[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pa_x, pc_x, pc_z, gsl0_471, gsk_381, gsl1_471, \
                         hsi0_290, hsi1_290, hsk_375, hsk_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pa_x[k] * gsl0_471[k]
                   + f_16 * gsk_381[k]
                   - f_14 * pc_x[k] * gsl1_471[k];

        t_472[k] = f_3 * pc_z[k] * hsk_375[k];

        t_473[k] = f_4 * hsi0_290[k]
                   - f_5 * hsi1_290[k]
                   + f_3 * pc_z[k] * hsk_376[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_y, pc_z, gsk_236, hsi0_291, hsi0_292, \
                         hsi0_294, hsi1_291, hsi1_292, hsi1_294, hsk_377, hsk_378, \
                         hsk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_6 * hsi0_291[k]
                   - f_7 * hsi1_291[k]
                   + f_3 * pc_z[k] * hsk_377[k];

        t_475[k] = f_8 * hsi0_292[k]
                   - f_9 * hsi1_292[k]
                   + f_3 * pc_z[k] * hsk_378[k];

        t_476[k] = f_18 * gsk_236[k]
                   + f_3 * pc_y[k] * hsk_380[k];

        t_477[k] = f_12 * hsi0_294[k]
                   - f_13 * hsi1_294[k]
                   + f_3 * pc_z[k] * hsk_380[k];
    }
}

static auto
compute_prim_hsl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsl0,
                                                          const size_t gsk, const size_t gsl1,
                                                          const size_t hsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_22 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsl0_270 = buffer.data(gsl0 + 270);
    const auto *gsl0_273 = buffer.data(gsl0 + 273);
    const auto *gsl0_276 = buffer.data(gsl0 + 276);
    const auto *gsl0_280 = buffer.data(gsl0 + 280);
    const auto *gsl0_285 = buffer.data(gsl0 + 285);
    const auto *gsl0_291 = buffer.data(gsl0 + 291);
    const auto *gsl0_405 = buffer.data(gsl0 + 405);
    const auto *gsl0_410 = buffer.data(gsl0 + 410);
    const auto *gsl0_414 = buffer.data(gsl0 + 414);
    const auto *gsl0_486 = buffer.data(gsl0 + 486);
    const auto *gsl0_488 = buffer.data(gsl0 + 488);
    const auto *gsl0_489 = buffer.data(gsl0 + 489);
    const auto *gsl0_490 = buffer.data(gsl0 + 490);
    const auto *gsl0_491 = buffer.data(gsl0 + 491);
    const auto *gsl0_492 = buffer.data(gsl0 + 492);
    const auto *gsl0_494 = buffer.data(gsl0 + 494);
    const auto *gsl0_500 = buffer.data(gsl0 + 500);
    const auto *gsl0_504 = buffer.data(gsl0 + 504);
    const auto *gsl0_507 = buffer.data(gsl0 + 507);
    const auto *gsl0_509 = buffer.data(gsl0 + 509);
    const auto *gsl0_512 = buffer.data(gsl0 + 512);
    const auto *gsl0_513 = buffer.data(gsl0 + 513);
    const auto *gsl0_515 = buffer.data(gsl0 + 515);
    const auto *gsl0_518 = buffer.data(gsl0 + 518);
    const auto *gsl0_519 = buffer.data(gsl0 + 519);
    const auto *gsl0_520 = buffer.data(gsl0 + 520);
    const auto *gsl0_522 = buffer.data(gsl0 + 522);
    const auto *gsl0_531 = buffer.data(gsl0 + 531);
    const auto *gsl0_533 = buffer.data(gsl0 + 533);
    const auto *gsl0_534 = buffer.data(gsl0 + 534);
    const auto *gsl0_535 = buffer.data(gsl0 + 535);
    const auto *gsl0_536 = buffer.data(gsl0 + 536);
    const auto *gsl0_537 = buffer.data(gsl0 + 537);
    const auto *gsl0_539 = buffer.data(gsl0 + 539);
    const auto *gsl0_540 = buffer.data(gsl0 + 540);
    const auto *gsl0_543 = buffer.data(gsl0 + 543);
    const auto *gsl0_545 = buffer.data(gsl0 + 545);
    const auto *gsl0_546 = buffer.data(gsl0 + 546);
    const auto *gsl0_549 = buffer.data(gsl0 + 549);
    const auto *gsl0_550 = buffer.data(gsl0 + 550);
    const auto *gsl0_552 = buffer.data(gsl0 + 552);
    const auto *gsl0_554 = buffer.data(gsl0 + 554);
    const auto *gsl0_555 = buffer.data(gsl0 + 555);
    const auto *gsl0_557 = buffer.data(gsl0 + 557);
    const auto *gsl0_558 = buffer.data(gsl0 + 558);
    const auto *gsl0_560 = buffer.data(gsl0 + 560);
    const auto *gsl0_561 = buffer.data(gsl0 + 561);
    const auto *gsl0_563 = buffer.data(gsl0 + 563);
    const auto *gsl0_564 = buffer.data(gsl0 + 564);
    const auto *gsl0_565 = buffer.data(gsl0 + 565);
    const auto *gsl0_567 = buffer.data(gsl0 + 567);
    const auto *gsl0_576 = buffer.data(gsl0 + 576);
    const auto *gsl0_578 = buffer.data(gsl0 + 578);
    const auto *gsl0_579 = buffer.data(gsl0 + 579);
    const auto *gsl0_580 = buffer.data(gsl0 + 580);
    const auto *gsl0_581 = buffer.data(gsl0 + 581);
    const auto *gsl0_582 = buffer.data(gsl0 + 582);
    const auto *gsl0_584 = buffer.data(gsl0 + 584);
    const auto *gsl0_588 = buffer.data(gsl0 + 588);
    const auto *gsl0_591 = buffer.data(gsl0 + 591);
    const auto *gsl0_595 = buffer.data(gsl0 + 595);

    const auto *gsk_216 = buffer.data(gsk + 216);
    const auto *gsk_219 = buffer.data(gsk + 219);
    const auto *gsk_222 = buffer.data(gsk + 222);
    const auto *gsk_226 = buffer.data(gsk + 226);
    const auto *gsk_231 = buffer.data(gsk + 231);
    const auto *gsk_244 = buffer.data(gsk + 244);
    const auto *gsk_251 = buffer.data(gsk + 251);
    const auto *gsk_252 = buffer.data(gsk + 252);
    const auto *gsk_254 = buffer.data(gsk + 254);
    const auto *gsk_255 = buffer.data(gsk + 255);
    const auto *gsk_257 = buffer.data(gsk + 257);
    const auto *gsk_258 = buffer.data(gsk + 258);
    const auto *gsk_261 = buffer.data(gsk + 261);
    const auto *gsk_262 = buffer.data(gsk + 262);
    const auto *gsk_266 = buffer.data(gsk + 266);
    const auto *gsk_267 = buffer.data(gsk + 267);
    const auto *gsk_272 = buffer.data(gsk + 272);
    const auto *gsk_280 = buffer.data(gsk + 280);
    const auto *gsk_287 = buffer.data(gsk + 287);
    const auto *gsk_288 = buffer.data(gsk + 288);
    const auto *gsk_290 = buffer.data(gsk + 290);
    const auto *gsk_291 = buffer.data(gsk + 291);
    const auto *gsk_293 = buffer.data(gsk + 293);
    const auto *gsk_294 = buffer.data(gsk + 294);
    const auto *gsk_297 = buffer.data(gsk + 297);
    const auto *gsk_302 = buffer.data(gsk + 302);
    const auto *gsk_308 = buffer.data(gsk + 308);
    const auto *gsk_323 = buffer.data(gsk + 323);
    const auto *gsk_324 = buffer.data(gsk + 324);
    const auto *gsk_326 = buffer.data(gsk + 326);
    const auto *gsk_329 = buffer.data(gsk + 329);
    const auto *gsk_388 = buffer.data(gsk + 388);
    const auto *gsk_390 = buffer.data(gsk + 390);
    const auto *gsk_391 = buffer.data(gsk + 391);
    const auto *gsk_392 = buffer.data(gsk + 392);
    const auto *gsk_393 = buffer.data(gsk + 393);
    const auto *gsk_394 = buffer.data(gsk + 394);
    const auto *gsk_395 = buffer.data(gsk + 395);
    const auto *gsk_401 = buffer.data(gsk + 401);
    const auto *gsk_405 = buffer.data(gsk + 405);
    const auto *gsk_408 = buffer.data(gsk + 408);
    const auto *gsk_410 = buffer.data(gsk + 410);
    const auto *gsk_413 = buffer.data(gsk + 413);
    const auto *gsk_414 = buffer.data(gsk + 414);
    const auto *gsk_416 = buffer.data(gsk + 416);
    const auto *gsk_419 = buffer.data(gsk + 419);
    const auto *gsk_420 = buffer.data(gsk + 420);
    const auto *gsk_421 = buffer.data(gsk + 421);
    const auto *gsk_423 = buffer.data(gsk + 423);
    const auto *gsk_424 = buffer.data(gsk + 424);
    const auto *gsk_425 = buffer.data(gsk + 425);
    const auto *gsk_426 = buffer.data(gsk + 426);
    const auto *gsk_427 = buffer.data(gsk + 427);
    const auto *gsk_428 = buffer.data(gsk + 428);
    const auto *gsk_429 = buffer.data(gsk + 429);
    const auto *gsk_430 = buffer.data(gsk + 430);
    const auto *gsk_431 = buffer.data(gsk + 431);
    const auto *gsk_432 = buffer.data(gsk + 432);
    const auto *gsk_435 = buffer.data(gsk + 435);
    const auto *gsk_437 = buffer.data(gsk + 437);
    const auto *gsk_438 = buffer.data(gsk + 438);
    const auto *gsk_441 = buffer.data(gsk + 441);
    const auto *gsk_442 = buffer.data(gsk + 442);
    const auto *gsk_444 = buffer.data(gsk + 444);
    const auto *gsk_446 = buffer.data(gsk + 446);
    const auto *gsk_447 = buffer.data(gsk + 447);
    const auto *gsk_449 = buffer.data(gsk + 449);
    const auto *gsk_450 = buffer.data(gsk + 450);
    const auto *gsk_452 = buffer.data(gsk + 452);
    const auto *gsk_453 = buffer.data(gsk + 453);
    const auto *gsk_455 = buffer.data(gsk + 455);
    const auto *gsk_456 = buffer.data(gsk + 456);
    const auto *gsk_457 = buffer.data(gsk + 457);
    const auto *gsk_459 = buffer.data(gsk + 459);
    const auto *gsk_460 = buffer.data(gsk + 460);
    const auto *gsk_461 = buffer.data(gsk + 461);
    const auto *gsk_462 = buffer.data(gsk + 462);
    const auto *gsk_463 = buffer.data(gsk + 463);
    const auto *gsk_464 = buffer.data(gsk + 464);
    const auto *gsk_465 = buffer.data(gsk + 465);
    const auto *gsk_466 = buffer.data(gsk + 466);
    const auto *gsk_467 = buffer.data(gsk + 467);
    const auto *gsk_471 = buffer.data(gsk + 471);
    const auto *gsk_474 = buffer.data(gsk + 474);
    const auto *gsk_478 = buffer.data(gsk + 478);

    const auto *gsl1_270 = buffer.data(gsl1 + 270);
    const auto *gsl1_273 = buffer.data(gsl1 + 273);
    const auto *gsl1_276 = buffer.data(gsl1 + 276);
    const auto *gsl1_280 = buffer.data(gsl1 + 280);
    const auto *gsl1_285 = buffer.data(gsl1 + 285);
    const auto *gsl1_291 = buffer.data(gsl1 + 291);
    const auto *gsl1_405 = buffer.data(gsl1 + 405);
    const auto *gsl1_410 = buffer.data(gsl1 + 410);
    const auto *gsl1_414 = buffer.data(gsl1 + 414);
    const auto *gsl1_486 = buffer.data(gsl1 + 486);
    const auto *gsl1_488 = buffer.data(gsl1 + 488);
    const auto *gsl1_489 = buffer.data(gsl1 + 489);
    const auto *gsl1_490 = buffer.data(gsl1 + 490);
    const auto *gsl1_491 = buffer.data(gsl1 + 491);
    const auto *gsl1_492 = buffer.data(gsl1 + 492);
    const auto *gsl1_494 = buffer.data(gsl1 + 494);
    const auto *gsl1_500 = buffer.data(gsl1 + 500);
    const auto *gsl1_504 = buffer.data(gsl1 + 504);
    const auto *gsl1_507 = buffer.data(gsl1 + 507);
    const auto *gsl1_509 = buffer.data(gsl1 + 509);
    const auto *gsl1_512 = buffer.data(gsl1 + 512);
    const auto *gsl1_513 = buffer.data(gsl1 + 513);
    const auto *gsl1_515 = buffer.data(gsl1 + 515);
    const auto *gsl1_518 = buffer.data(gsl1 + 518);
    const auto *gsl1_519 = buffer.data(gsl1 + 519);
    const auto *gsl1_520 = buffer.data(gsl1 + 520);
    const auto *gsl1_522 = buffer.data(gsl1 + 522);
    const auto *gsl1_531 = buffer.data(gsl1 + 531);
    const auto *gsl1_533 = buffer.data(gsl1 + 533);
    const auto *gsl1_534 = buffer.data(gsl1 + 534);
    const auto *gsl1_535 = buffer.data(gsl1 + 535);
    const auto *gsl1_536 = buffer.data(gsl1 + 536);
    const auto *gsl1_537 = buffer.data(gsl1 + 537);
    const auto *gsl1_539 = buffer.data(gsl1 + 539);
    const auto *gsl1_540 = buffer.data(gsl1 + 540);
    const auto *gsl1_543 = buffer.data(gsl1 + 543);
    const auto *gsl1_545 = buffer.data(gsl1 + 545);
    const auto *gsl1_546 = buffer.data(gsl1 + 546);
    const auto *gsl1_549 = buffer.data(gsl1 + 549);
    const auto *gsl1_550 = buffer.data(gsl1 + 550);
    const auto *gsl1_552 = buffer.data(gsl1 + 552);
    const auto *gsl1_554 = buffer.data(gsl1 + 554);
    const auto *gsl1_555 = buffer.data(gsl1 + 555);
    const auto *gsl1_557 = buffer.data(gsl1 + 557);
    const auto *gsl1_558 = buffer.data(gsl1 + 558);
    const auto *gsl1_560 = buffer.data(gsl1 + 560);
    const auto *gsl1_561 = buffer.data(gsl1 + 561);
    const auto *gsl1_563 = buffer.data(gsl1 + 563);
    const auto *gsl1_564 = buffer.data(gsl1 + 564);
    const auto *gsl1_565 = buffer.data(gsl1 + 565);
    const auto *gsl1_567 = buffer.data(gsl1 + 567);
    const auto *gsl1_576 = buffer.data(gsl1 + 576);
    const auto *gsl1_578 = buffer.data(gsl1 + 578);
    const auto *gsl1_579 = buffer.data(gsl1 + 579);
    const auto *gsl1_580 = buffer.data(gsl1 + 580);
    const auto *gsl1_581 = buffer.data(gsl1 + 581);
    const auto *gsl1_582 = buffer.data(gsl1 + 582);
    const auto *gsl1_584 = buffer.data(gsl1 + 584);
    const auto *gsl1_588 = buffer.data(gsl1 + 588);
    const auto *gsl1_591 = buffer.data(gsl1 + 591);
    const auto *gsl1_595 = buffer.data(gsl1 + 595);

    const auto *hsk_381 = buffer.data(hsk + 381);
    const auto *hsk_388 = buffer.data(hsk + 388);
    const auto *hsk_390 = buffer.data(hsk + 390);
    const auto *hsk_391 = buffer.data(hsk + 391);
    const auto *hsk_392 = buffer.data(hsk + 392);
    const auto *hsk_393 = buffer.data(hsk + 393);
    const auto *hsk_394 = buffer.data(hsk + 394);
    const auto *hsk_395 = buffer.data(hsk + 395);
    const auto *hsk_396 = buffer.data(hsk + 396);
    const auto *hsk_398 = buffer.data(hsk + 398);
    const auto *hsk_399 = buffer.data(hsk + 399);
    const auto *hsk_401 = buffer.data(hsk + 401);
    const auto *hsk_402 = buffer.data(hsk + 402);
    const auto *hsk_405 = buffer.data(hsk + 405);
    const auto *hsk_406 = buffer.data(hsk + 406);
    const auto *hsk_410 = buffer.data(hsk + 410);
    const auto *hsk_411 = buffer.data(hsk + 411);
    const auto *hsk_416 = buffer.data(hsk + 416);
    const auto *hsk_424 = buffer.data(hsk + 424);
    const auto *hsk_425 = buffer.data(hsk + 425);
    const auto *hsk_426 = buffer.data(hsk + 426);
    const auto *hsk_427 = buffer.data(hsk + 427);
    const auto *hsk_428 = buffer.data(hsk + 428);
    const auto *hsk_429 = buffer.data(hsk + 429);
    const auto *hsk_430 = buffer.data(hsk + 430);
    const auto *hsk_431 = buffer.data(hsk + 431);
    const auto *hsk_432 = buffer.data(hsk + 432);
    const auto *hsk_434 = buffer.data(hsk + 434);
    const auto *hsk_435 = buffer.data(hsk + 435);
    const auto *hsk_437 = buffer.data(hsk + 437);
    const auto *hsk_438 = buffer.data(hsk + 438);
    const auto *hsk_441 = buffer.data(hsk + 441);
    const auto *hsk_442 = buffer.data(hsk + 442);
    const auto *hsk_446 = buffer.data(hsk + 446);
    const auto *hsk_447 = buffer.data(hsk + 447);
    const auto *hsk_452 = buffer.data(hsk + 452);
    const auto *hsk_460 = buffer.data(hsk + 460);
    const auto *hsk_461 = buffer.data(hsk + 461);
    const auto *hsk_462 = buffer.data(hsk + 462);
    const auto *hsk_463 = buffer.data(hsk + 463);
    const auto *hsk_464 = buffer.data(hsk + 464);
    const auto *hsk_465 = buffer.data(hsk + 465);
    const auto *hsk_466 = buffer.data(hsk + 466);
    const auto *hsk_467 = buffer.data(hsk + 467);
    const auto *hsk_468 = buffer.data(hsk + 468);
    const auto *hsk_470 = buffer.data(hsk + 470);
    const auto *hsk_471 = buffer.data(hsk + 471);
    const auto *hsk_473 = buffer.data(hsk + 473);
    const auto *hsk_474 = buffer.data(hsk + 474);

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, pc_x, pc_z, gsk_388, gsk_390, \
                         gsk_391, gsk_392, hsk_381, hsk_388, hsk_390, hsk_391, \
                         hsk_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_15 * gsk_388[k]
                   + f_3 * pc_x[k] * hsk_388[k];

        t_479[k] = f_3 * pc_z[k] * hsk_381[k];

        t_480[k] = f_15 * gsk_390[k]
                   + f_3 * pc_x[k] * hsk_390[k];

        t_481[k] = f_15 * gsk_391[k]
                   + f_3 * pc_x[k] * hsk_391[k];

        t_482[k] = f_15 * gsk_392[k]
                   + f_3 * pc_x[k] * hsk_392[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pa_x, pc_x, gsl0_486, gsk_393, gsk_394, \
                         gsk_395, gsl1_486, hsk_393, hsk_394, hsk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_15 * gsk_393[k]
                   + f_3 * pc_x[k] * hsk_393[k];

        t_484[k] = f_15 * gsk_394[k]
                   + f_3 * pc_x[k] * hsk_394[k];

        t_485[k] = f_15 * gsk_395[k]
                   + f_3 * pc_x[k] * hsk_395[k];

        t_486[k] = pa_x[k] * gsl0_486[k]
                   - f_14 * pc_x[k] * gsl1_486[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pa_x, pc_x, pc_z, gsl0_488, gsl0_489, \
                         gsl0_490, gsl1_488, gsl1_489, gsl1_490, \
                         hsk_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_3 * pc_z[k] * hsk_388[k];

        t_488[k] = pa_x[k] * gsl0_488[k]
                   - f_14 * pc_x[k] * gsl1_488[k];

        t_489[k] = pa_x[k] * gsl0_489[k]
                   - f_14 * pc_x[k] * gsl1_489[k];

        t_490[k] = pa_x[k] * gsl0_490[k]
                   - f_14 * pc_x[k] * gsl1_490[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pa_x, pc_x, pc_y, gsl0_491, gsl0_492, \
                         gsl0_494, gsk_251, gsl1_491, gsl1_492, gsl1_494, \
                         hsk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = pa_x[k] * gsl0_491[k]
                   - f_14 * pc_x[k] * gsl1_491[k];

        t_492[k] = pa_x[k] * gsl0_492[k]
                   - f_14 * pc_x[k] * gsl1_492[k];

        t_493[k] = f_18 * gsk_251[k]
                   + f_3 * pc_y[k] * hsk_395[k];

        t_494[k] = pa_x[k] * gsl0_494[k]
                   - f_14 * pc_x[k] * gsl1_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pa_z, pc_y, pc_z, gsl0_270, gsl0_273, \
                         gsk_216, gsk_252, gsl1_270, gsl1_273, \
                         hsk_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * gsl0_270[k]
                   - f_14 * pc_z[k] * gsl1_270[k];

        t_496[k] = f_17 * gsk_252[k]
                   + f_3 * pc_y[k] * hsk_396[k];

        t_497[k] = f_15 * gsk_216[k]
                   + f_3 * pc_z[k] * hsk_396[k];

        t_498[k] = pa_z[k] * gsl0_273[k]
                   - f_14 * pc_z[k] * gsl1_273[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pa_x, pa_z, pc_x, pc_y, pc_z, gsl0_276, \
                         gsl0_500, gsk_254, gsk_401, gsl1_276, gsl1_500, \
                         hsk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_17 * gsk_254[k]
                   + f_3 * pc_y[k] * hsk_398[k];

        t_500[k] = pa_x[k] * gsl0_500[k]
                   + f_19 * gsk_401[k]
                   - f_14 * pc_x[k] * gsl1_500[k];

        t_501[k] = pa_z[k] * gsl0_276[k]
                   - f_14 * pc_z[k] * gsl1_276[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pa_x, pc_x, pc_y, pc_z, gsl0_504, gsk_219, \
                         gsk_257, gsk_405, gsl1_504, hsk_399, hsk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_15 * gsk_219[k]
                   + f_3 * pc_z[k] * hsk_399[k];

        t_503[k] = f_17 * gsk_257[k]
                   + f_3 * pc_y[k] * hsk_401[k];

        t_504[k] = pa_x[k] * gsl0_504[k]
                   + f_0 * gsk_405[k]
                   - f_14 * pc_x[k] * gsl1_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pa_x, pa_z, pc_x, pc_z, gsl0_280, gsl0_507, \
                         gsk_222, gsk_408, gsl1_280, gsl1_507, \
                         hsk_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pa_z[k] * gsl0_280[k]
                   - f_14 * pc_z[k] * gsl1_280[k];

        t_506[k] = f_15 * gsk_222[k]
                   + f_3 * pc_z[k] * hsk_402[k];

        t_507[k] = pa_x[k] * gsl0_507[k]
                   + f_18 * gsk_408[k]
                   - f_14 * pc_x[k] * gsl1_507[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pa_x, pa_z, pc_x, pc_y, pc_z, gsl0_285, \
                         gsl0_509, gsk_261, gsk_410, gsl1_285, gsl1_509, \
                         hsk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_17 * gsk_261[k]
                   + f_3 * pc_y[k] * hsk_405[k];

        t_509[k] = pa_x[k] * gsl0_509[k]
                   + f_18 * gsk_410[k]
                   - f_14 * pc_x[k] * gsl1_509[k];

        t_510[k] = pa_z[k] * gsl0_285[k]
                   - f_14 * pc_z[k] * gsl1_285[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pa_x, pc_x, pc_z, gsl0_512, gsl0_513, gsk_226, \
                         gsk_413, gsk_414, gsl1_512, gsl1_513, \
                         hsk_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_15 * gsk_226[k]
                   + f_3 * pc_z[k] * hsk_406[k];

        t_512[k] = pa_x[k] * gsl0_512[k]
                   + f_17 * gsk_413[k]
                   - f_14 * pc_x[k] * gsl1_512[k];

        t_513[k] = pa_x[k] * gsl0_513[k]
                   + f_17 * gsk_414[k]
                   - f_14 * pc_x[k] * gsl1_513[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pa_x, pa_z, pc_x, pc_y, pc_z, gsl0_291, \
                         gsl0_515, gsk_266, gsk_416, gsl1_291, gsl1_515, \
                         hsk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_17 * gsk_266[k]
                   + f_3 * pc_y[k] * hsk_410[k];

        t_515[k] = pa_x[k] * gsl0_515[k]
                   + f_17 * gsk_416[k]
                   - f_14 * pc_x[k] * gsl1_515[k];

        t_516[k] = pa_z[k] * gsl0_291[k]
                   - f_14 * pc_z[k] * gsl1_291[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pa_x, pc_x, pc_z, gsl0_518, gsl0_519, gsk_231, \
                         gsk_419, gsk_420, gsl1_518, gsl1_519, \
                         hsk_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_15 * gsk_231[k]
                   + f_3 * pc_z[k] * hsk_411[k];

        t_518[k] = pa_x[k] * gsl0_518[k]
                   + f_16 * gsk_419[k]
                   - f_14 * pc_x[k] * gsl1_518[k];

        t_519[k] = pa_x[k] * gsl0_519[k]
                   + f_16 * gsk_420[k]
                   - f_14 * pc_x[k] * gsl1_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pa_x, pc_x, pc_y, gsl0_520, gsl0_522, gsk_272, \
                         gsk_421, gsk_423, gsl1_520, gsl1_522, \
                         hsk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = pa_x[k] * gsl0_520[k]
                   + f_16 * gsk_421[k]
                   - f_14 * pc_x[k] * gsl1_520[k];

        t_521[k] = f_17 * gsk_272[k]
                   + f_3 * pc_y[k] * hsk_416[k];

        t_522[k] = pa_x[k] * gsl0_522[k]
                   + f_16 * gsk_423[k]
                   - f_14 * pc_x[k] * gsl1_522[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, t_527, pc_x, gsk_424, gsk_425, gsk_426, \
                         gsk_427, gsk_428, hsk_424, hsk_425, hsk_426, hsk_427, \
                         hsk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_15 * gsk_424[k]
                   + f_3 * pc_x[k] * hsk_424[k];

        t_524[k] = f_15 * gsk_425[k]
                   + f_3 * pc_x[k] * hsk_425[k];

        t_525[k] = f_15 * gsk_426[k]
                   + f_3 * pc_x[k] * hsk_426[k];

        t_526[k] = f_15 * gsk_427[k]
                   + f_3 * pc_x[k] * hsk_427[k];

        t_527[k] = f_15 * gsk_428[k]
                   + f_3 * pc_x[k] * hsk_428[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, pa_x, pc_x, gsl0_531, gsk_429, gsk_430, \
                         gsk_431, gsl1_531, hsk_429, hsk_430, hsk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_15 * gsk_429[k]
                   + f_3 * pc_x[k] * hsk_429[k];

        t_529[k] = f_15 * gsk_430[k]
                   + f_3 * pc_x[k] * hsk_430[k];

        t_530[k] = f_15 * gsk_431[k]
                   + f_3 * pc_x[k] * hsk_431[k];

        t_531[k] = pa_x[k] * gsl0_531[k]
                   - f_14 * pc_x[k] * gsl1_531[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pa_x, pc_x, pc_z, gsl0_533, gsl0_534, \
                         gsl0_535, gsk_244, gsl1_533, gsl1_534, gsl1_535, \
                         hsk_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_15 * gsk_244[k]
                   + f_3 * pc_z[k] * hsk_424[k];

        t_533[k] = pa_x[k] * gsl0_533[k]
                   - f_14 * pc_x[k] * gsl1_533[k];

        t_534[k] = pa_x[k] * gsl0_534[k]
                   - f_14 * pc_x[k] * gsl1_534[k];

        t_535[k] = pa_x[k] * gsl0_535[k]
                   - f_14 * pc_x[k] * gsl1_535[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pa_x, pc_x, pc_y, gsl0_536, gsl0_537, \
                         gsl0_539, gsk_287, gsl1_536, gsl1_537, gsl1_539, \
                         hsk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = pa_x[k] * gsl0_536[k]
                   - f_14 * pc_x[k] * gsl1_536[k];

        t_537[k] = pa_x[k] * gsl0_537[k]
                   - f_14 * pc_x[k] * gsl1_537[k];

        t_538[k] = f_17 * gsk_287[k]
                   + f_3 * pc_y[k] * hsk_431[k];

        t_539[k] = pa_x[k] * gsl0_539[k]
                   - f_14 * pc_x[k] * gsl1_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pa_x, pc_x, pc_y, pc_z, gsl0_540, gsk_252, \
                         gsk_288, gsk_432, gsl1_540, hsk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = pa_x[k] * gsl0_540[k]
                   + f_22 * gsk_432[k]
                   - f_14 * pc_x[k] * gsl1_540[k];

        t_541[k] = f_16 * gsk_288[k]
                   + f_3 * pc_y[k] * hsk_432[k];

        t_542[k] = f_16 * gsk_252[k]
                   + f_3 * pc_z[k] * hsk_432[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pa_x, pc_x, pc_y, gsl0_543, gsl0_545, gsk_290, \
                         gsk_435, gsk_437, gsl1_543, gsl1_545, \
                         hsk_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = pa_x[k] * gsl0_543[k]
                   + f_19 * gsk_435[k]
                   - f_14 * pc_x[k] * gsl1_543[k];

        t_544[k] = f_16 * gsk_290[k]
                   + f_3 * pc_y[k] * hsk_434[k];

        t_545[k] = pa_x[k] * gsl0_545[k]
                   + f_19 * gsk_437[k]
                   - f_14 * pc_x[k] * gsl1_545[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, pa_x, pc_x, pc_y, pc_z, gsl0_546, gsk_255, \
                         gsk_293, gsk_438, gsl1_546, hsk_435, hsk_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_x[k] * gsl0_546[k]
                   + f_0 * gsk_438[k]
                   - f_14 * pc_x[k] * gsl1_546[k];

        t_547[k] = f_16 * gsk_255[k]
                   + f_3 * pc_z[k] * hsk_435[k];

        t_548[k] = f_16 * gsk_293[k]
                   + f_3 * pc_y[k] * hsk_437[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pa_x, pc_x, pc_z, gsl0_549, gsl0_550, gsk_258, \
                         gsk_441, gsk_442, gsl1_549, gsl1_550, \
                         hsk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = pa_x[k] * gsl0_549[k]
                   + f_0 * gsk_441[k]
                   - f_14 * pc_x[k] * gsl1_549[k];

        t_550[k] = pa_x[k] * gsl0_550[k]
                   + f_18 * gsk_442[k]
                   - f_14 * pc_x[k] * gsl1_550[k];

        t_551[k] = f_16 * gsk_258[k]
                   + f_3 * pc_z[k] * hsk_438[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pa_x, pc_x, pc_y, gsl0_552, gsl0_554, gsk_297, \
                         gsk_444, gsk_446, gsl1_552, gsl1_554, \
                         hsk_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_x[k] * gsl0_552[k]
                   + f_18 * gsk_444[k]
                   - f_14 * pc_x[k] * gsl1_552[k];

        t_553[k] = f_16 * gsk_297[k]
                   + f_3 * pc_y[k] * hsk_441[k];

        t_554[k] = pa_x[k] * gsl0_554[k]
                   + f_18 * gsk_446[k]
                   - f_14 * pc_x[k] * gsl1_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pa_x, pc_x, pc_z, gsl0_555, gsl0_557, gsk_262, \
                         gsk_447, gsk_449, gsl1_555, gsl1_557, \
                         hsk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pa_x[k] * gsl0_555[k]
                   + f_17 * gsk_447[k]
                   - f_14 * pc_x[k] * gsl1_555[k];

        t_556[k] = f_16 * gsk_262[k]
                   + f_3 * pc_z[k] * hsk_442[k];

        t_557[k] = pa_x[k] * gsl0_557[k]
                   + f_17 * gsk_449[k]
                   - f_14 * pc_x[k] * gsl1_557[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, pa_x, pc_x, pc_y, gsl0_558, gsl0_560, gsk_302, \
                         gsk_450, gsk_452, gsl1_558, gsl1_560, \
                         hsk_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = pa_x[k] * gsl0_558[k]
                   + f_17 * gsk_450[k]
                   - f_14 * pc_x[k] * gsl1_558[k];

        t_559[k] = f_16 * gsk_302[k]
                   + f_3 * pc_y[k] * hsk_446[k];

        t_560[k] = pa_x[k] * gsl0_560[k]
                   + f_17 * gsk_452[k]
                   - f_14 * pc_x[k] * gsl1_560[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pa_x, pc_x, pc_z, gsl0_561, gsl0_563, gsk_267, \
                         gsk_453, gsk_455, gsl1_561, gsl1_563, \
                         hsk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = pa_x[k] * gsl0_561[k]
                   + f_16 * gsk_453[k]
                   - f_14 * pc_x[k] * gsl1_561[k];

        t_562[k] = f_16 * gsk_267[k]
                   + f_3 * pc_z[k] * hsk_447[k];

        t_563[k] = pa_x[k] * gsl0_563[k]
                   + f_16 * gsk_455[k]
                   - f_14 * pc_x[k] * gsl1_563[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_x, pc_x, pc_y, gsl0_564, gsl0_565, gsk_308, \
                         gsk_456, gsk_457, gsl1_564, gsl1_565, \
                         hsk_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = pa_x[k] * gsl0_564[k]
                   + f_16 * gsk_456[k]
                   - f_14 * pc_x[k] * gsl1_564[k];

        t_565[k] = pa_x[k] * gsl0_565[k]
                   + f_16 * gsk_457[k]
                   - f_14 * pc_x[k] * gsl1_565[k];

        t_566[k] = f_16 * gsk_308[k]
                   + f_3 * pc_y[k] * hsk_452[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pa_x, pc_x, gsl0_567, gsk_459, gsk_460, \
                         gsk_461, gsk_462, gsl1_567, hsk_460, hsk_461, \
                         hsk_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = pa_x[k] * gsl0_567[k]
                   + f_16 * gsk_459[k]
                   - f_14 * pc_x[k] * gsl1_567[k];

        t_568[k] = f_15 * gsk_460[k]
                   + f_3 * pc_x[k] * hsk_460[k];

        t_569[k] = f_15 * gsk_461[k]
                   + f_3 * pc_x[k] * hsk_461[k];

        t_570[k] = f_15 * gsk_462[k]
                   + f_3 * pc_x[k] * hsk_462[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, t_575, pc_x, gsk_463, gsk_464, gsk_465, \
                         gsk_466, gsk_467, hsk_463, hsk_464, hsk_465, hsk_466, \
                         hsk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_15 * gsk_463[k]
                   + f_3 * pc_x[k] * hsk_463[k];

        t_572[k] = f_15 * gsk_464[k]
                   + f_3 * pc_x[k] * hsk_464[k];

        t_573[k] = f_15 * gsk_465[k]
                   + f_3 * pc_x[k] * hsk_465[k];

        t_574[k] = f_15 * gsk_466[k]
                   + f_3 * pc_x[k] * hsk_466[k];

        t_575[k] = f_15 * gsk_467[k]
                   + f_3 * pc_x[k] * hsk_467[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pa_x, pc_x, pc_z, gsl0_576, gsl0_578, \
                         gsl0_579, gsk_280, gsl1_576, gsl1_578, gsl1_579, \
                         hsk_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = pa_x[k] * gsl0_576[k]
                   - f_14 * pc_x[k] * gsl1_576[k];

        t_577[k] = f_16 * gsk_280[k]
                   + f_3 * pc_z[k] * hsk_460[k];

        t_578[k] = pa_x[k] * gsl0_578[k]
                   - f_14 * pc_x[k] * gsl1_578[k];

        t_579[k] = pa_x[k] * gsl0_579[k]
                   - f_14 * pc_x[k] * gsl1_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pa_x, pc_x, pc_y, gsl0_580, gsl0_581, \
                         gsl0_582, gsk_323, gsl1_580, gsl1_581, gsl1_582, \
                         hsk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = pa_x[k] * gsl0_580[k]
                   - f_14 * pc_x[k] * gsl1_580[k];

        t_581[k] = pa_x[k] * gsl0_581[k]
                   - f_14 * pc_x[k] * gsl1_581[k];

        t_582[k] = pa_x[k] * gsl0_582[k]
                   - f_14 * pc_x[k] * gsl1_582[k];

        t_583[k] = f_16 * gsk_323[k]
                   + f_3 * pc_y[k] * hsk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_x, pa_y, pc_x, pc_y, pc_z, gsl0_405, \
                         gsl0_584, gsk_288, gsk_324, gsl1_405, gsl1_584, \
                         hsk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = pa_x[k] * gsl0_584[k]
                   - f_14 * pc_x[k] * gsl1_584[k];

        t_585[k] = pa_y[k] * gsl0_405[k]
                   - f_14 * pc_y[k] * gsl1_405[k];

        t_586[k] = f_15 * gsk_324[k]
                   + f_3 * pc_y[k] * hsk_468[k];

        t_587[k] = f_17 * gsk_288[k]
                   + f_3 * pc_z[k] * hsk_468[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, pa_x, pa_y, pc_x, pc_y, gsl0_410, gsl0_588, \
                         gsk_326, gsk_471, gsl1_410, gsl1_588, \
                         hsk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_x[k] * gsl0_588[k]
                   + f_19 * gsk_471[k]
                   - f_14 * pc_x[k] * gsl1_588[k];

        t_589[k] = f_15 * gsk_326[k]
                   + f_3 * pc_y[k] * hsk_470[k];

        t_590[k] = pa_y[k] * gsl0_410[k]
                   - f_14 * pc_y[k] * gsl1_410[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, pa_x, pc_x, pc_y, pc_z, gsl0_591, gsk_291, \
                         gsk_329, gsk_474, gsl1_591, hsk_471, hsk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = pa_x[k] * gsl0_591[k]
                   + f_0 * gsk_474[k]
                   - f_14 * pc_x[k] * gsl1_591[k];

        t_592[k] = f_17 * gsk_291[k]
                   + f_3 * pc_z[k] * hsk_471[k];

        t_593[k] = f_15 * gsk_329[k]
                   + f_3 * pc_y[k] * hsk_473[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pa_x, pa_y, pc_x, pc_y, pc_z, gsl0_414, \
                         gsl0_595, gsk_294, gsk_478, gsl1_414, gsl1_595, \
                         hsk_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = pa_y[k] * gsl0_414[k]
                   - f_14 * pc_y[k] * gsl1_414[k];

        t_595[k] = pa_x[k] * gsl0_595[k]
                   + f_18 * gsk_478[k]
                   - f_14 * pc_x[k] * gsl1_595[k];

        t_596[k] = f_17 * gsk_294[k]
                   + f_3 * pc_z[k] * hsk_474[k];
    }
}

static auto
compute_prim_hsl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsl0,
                                                          const size_t gsk, const size_t gsl1,
                                                          const size_t hsi0, const size_t hsi1,
                                                          const size_t hsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsl0_419 = buffer.data(gsl0 + 419);
    const auto *gsl0_425 = buffer.data(gsl0 + 425);
    const auto *gsl0_432 = buffer.data(gsl0 + 432);
    const auto *gsl0_450 = buffer.data(gsl0 + 450);
    const auto *gsl0_451 = buffer.data(gsl0 + 451);
    const auto *gsl0_597 = buffer.data(gsl0 + 597);
    const auto *gsl0_600 = buffer.data(gsl0 + 600);
    const auto *gsl0_602 = buffer.data(gsl0 + 602);
    const auto *gsl0_603 = buffer.data(gsl0 + 603);
    const auto *gsl0_606 = buffer.data(gsl0 + 606);
    const auto *gsl0_608 = buffer.data(gsl0 + 608);
    const auto *gsl0_609 = buffer.data(gsl0 + 609);
    const auto *gsl0_610 = buffer.data(gsl0 + 610);
    const auto *gsl0_621 = buffer.data(gsl0 + 621);
    const auto *gsl0_623 = buffer.data(gsl0 + 623);
    const auto *gsl0_624 = buffer.data(gsl0 + 624);
    const auto *gsl0_625 = buffer.data(gsl0 + 625);
    const auto *gsl0_626 = buffer.data(gsl0 + 626);
    const auto *gsl0_627 = buffer.data(gsl0 + 627);
    const auto *gsl0_629 = buffer.data(gsl0 + 629);
    const auto *gsl0_630 = buffer.data(gsl0 + 630);
    const auto *gsl0_635 = buffer.data(gsl0 + 635);
    const auto *gsl0_639 = buffer.data(gsl0 + 639);
    const auto *gsl0_644 = buffer.data(gsl0 + 644);
    const auto *gsl0_650 = buffer.data(gsl0 + 650);
    const auto *gsl0_657 = buffer.data(gsl0 + 657);
    const auto *gsl0_666 = buffer.data(gsl0 + 666);
    const auto *gsl0_667 = buffer.data(gsl0 + 667);
    const auto *gsl0_668 = buffer.data(gsl0 + 668);
    const auto *gsl0_669 = buffer.data(gsl0 + 669);
    const auto *gsl0_670 = buffer.data(gsl0 + 670);
    const auto *gsl0_671 = buffer.data(gsl0 + 671);
    const auto *gsl0_672 = buffer.data(gsl0 + 672);
    const auto *gsl0_674 = buffer.data(gsl0 + 674);

    const auto *gsk_298 = buffer.data(gsk + 298);
    const auto *gsk_303 = buffer.data(gsk + 303);
    const auto *gsk_316 = buffer.data(gsk + 316);
    const auto *gsk_324 = buffer.data(gsk + 324);
    const auto *gsk_333 = buffer.data(gsk + 333);
    const auto *gsk_338 = buffer.data(gsk + 338);
    const auto *gsk_344 = buffer.data(gsk + 344);
    const auto *gsk_359 = buffer.data(gsk + 359);
    const auto *gsk_388 = buffer.data(gsk + 388);
    const auto *gsk_395 = buffer.data(gsk + 395);
    const auto *gsk_480 = buffer.data(gsk + 480);
    const auto *gsk_483 = buffer.data(gsk + 483);
    const auto *gsk_485 = buffer.data(gsk + 485);
    const auto *gsk_486 = buffer.data(gsk + 486);
    const auto *gsk_489 = buffer.data(gsk + 489);
    const auto *gsk_491 = buffer.data(gsk + 491);
    const auto *gsk_492 = buffer.data(gsk + 492);
    const auto *gsk_493 = buffer.data(gsk + 493);
    const auto *gsk_496 = buffer.data(gsk + 496);
    const auto *gsk_497 = buffer.data(gsk + 497);
    const auto *gsk_498 = buffer.data(gsk + 498);
    const auto *gsk_499 = buffer.data(gsk + 499);
    const auto *gsk_500 = buffer.data(gsk + 500);
    const auto *gsk_501 = buffer.data(gsk + 501);
    const auto *gsk_502 = buffer.data(gsk + 502);
    const auto *gsk_503 = buffer.data(gsk + 503);
    const auto *gsk_504 = buffer.data(gsk + 504);
    const auto *gsk_509 = buffer.data(gsk + 509);
    const auto *gsk_513 = buffer.data(gsk + 513);
    const auto *gsk_518 = buffer.data(gsk + 518);
    const auto *gsk_524 = buffer.data(gsk + 524);
    const auto *gsk_531 = buffer.data(gsk + 531);
    const auto *gsk_532 = buffer.data(gsk + 532);
    const auto *gsk_533 = buffer.data(gsk + 533);
    const auto *gsk_534 = buffer.data(gsk + 534);
    const auto *gsk_535 = buffer.data(gsk + 535);
    const auto *gsk_536 = buffer.data(gsk + 536);
    const auto *gsk_537 = buffer.data(gsk + 537);
    const auto *gsk_539 = buffer.data(gsk + 539);

    const auto *gsl1_419 = buffer.data(gsl1 + 419);
    const auto *gsl1_425 = buffer.data(gsl1 + 425);
    const auto *gsl1_432 = buffer.data(gsl1 + 432);
    const auto *gsl1_450 = buffer.data(gsl1 + 450);
    const auto *gsl1_451 = buffer.data(gsl1 + 451);
    const auto *gsl1_597 = buffer.data(gsl1 + 597);
    const auto *gsl1_600 = buffer.data(gsl1 + 600);
    const auto *gsl1_602 = buffer.data(gsl1 + 602);
    const auto *gsl1_603 = buffer.data(gsl1 + 603);
    const auto *gsl1_606 = buffer.data(gsl1 + 606);
    const auto *gsl1_608 = buffer.data(gsl1 + 608);
    const auto *gsl1_609 = buffer.data(gsl1 + 609);
    const auto *gsl1_610 = buffer.data(gsl1 + 610);
    const auto *gsl1_621 = buffer.data(gsl1 + 621);
    const auto *gsl1_623 = buffer.data(gsl1 + 623);
    const auto *gsl1_624 = buffer.data(gsl1 + 624);
    const auto *gsl1_625 = buffer.data(gsl1 + 625);
    const auto *gsl1_626 = buffer.data(gsl1 + 626);
    const auto *gsl1_627 = buffer.data(gsl1 + 627);
    const auto *gsl1_629 = buffer.data(gsl1 + 629);
    const auto *gsl1_630 = buffer.data(gsl1 + 630);
    const auto *gsl1_635 = buffer.data(gsl1 + 635);
    const auto *gsl1_639 = buffer.data(gsl1 + 639);
    const auto *gsl1_644 = buffer.data(gsl1 + 644);
    const auto *gsl1_650 = buffer.data(gsl1 + 650);
    const auto *gsl1_657 = buffer.data(gsl1 + 657);
    const auto *gsl1_666 = buffer.data(gsl1 + 666);
    const auto *gsl1_667 = buffer.data(gsl1 + 667);
    const auto *gsl1_668 = buffer.data(gsl1 + 668);
    const auto *gsl1_669 = buffer.data(gsl1 + 669);
    const auto *gsl1_670 = buffer.data(gsl1 + 670);
    const auto *gsl1_671 = buffer.data(gsl1 + 671);
    const auto *gsl1_672 = buffer.data(gsl1 + 672);
    const auto *gsl1_674 = buffer.data(gsl1 + 674);

    const auto *hsi0_392 = buffer.data(hsi0 + 392);
    const auto *hsi0_393 = buffer.data(hsi0 + 393);
    const auto *hsi0_394 = buffer.data(hsi0 + 394);
    const auto *hsi0_395 = buffer.data(hsi0 + 395);
    const auto *hsi0_396 = buffer.data(hsi0 + 396);
    const auto *hsi0_397 = buffer.data(hsi0 + 397);
    const auto *hsi0_398 = buffer.data(hsi0 + 398);
    const auto *hsi0_399 = buffer.data(hsi0 + 399);
    const auto *hsi0_400 = buffer.data(hsi0 + 400);
    const auto *hsi0_401 = buffer.data(hsi0 + 401);
    const auto *hsi0_402 = buffer.data(hsi0 + 402);
    const auto *hsi0_403 = buffer.data(hsi0 + 403);
    const auto *hsi0_404 = buffer.data(hsi0 + 404);
    const auto *hsi0_405 = buffer.data(hsi0 + 405);
    const auto *hsi0_406 = buffer.data(hsi0 + 406);
    const auto *hsi0_420 = buffer.data(hsi0 + 420);
    const auto *hsi0_421 = buffer.data(hsi0 + 421);
    const auto *hsi0_423 = buffer.data(hsi0 + 423);
    const auto *hsi0_425 = buffer.data(hsi0 + 425);
    const auto *hsi0_426 = buffer.data(hsi0 + 426);
    const auto *hsi0_428 = buffer.data(hsi0 + 428);
    const auto *hsi0_429 = buffer.data(hsi0 + 429);
    const auto *hsi0_430 = buffer.data(hsi0 + 430);
    const auto *hsi0_432 = buffer.data(hsi0 + 432);
    const auto *hsi0_433 = buffer.data(hsi0 + 433);
    const auto *hsi0_434 = buffer.data(hsi0 + 434);
    const auto *hsi0_435 = buffer.data(hsi0 + 435);
    const auto *hsi0_437 = buffer.data(hsi0 + 437);
    const auto *hsi0_438 = buffer.data(hsi0 + 438);
    const auto *hsi0_439 = buffer.data(hsi0 + 439);
    const auto *hsi0_440 = buffer.data(hsi0 + 440);
    const auto *hsi0_441 = buffer.data(hsi0 + 441);
    const auto *hsi0_442 = buffer.data(hsi0 + 442);
    const auto *hsi0_443 = buffer.data(hsi0 + 443);
    const auto *hsi0_444 = buffer.data(hsi0 + 444);
    const auto *hsi0_445 = buffer.data(hsi0 + 445);
    const auto *hsi0_446 = buffer.data(hsi0 + 446);
    const auto *hsi0_447 = buffer.data(hsi0 + 447);

    const auto *hsi1_392 = buffer.data(hsi1 + 392);
    const auto *hsi1_393 = buffer.data(hsi1 + 393);
    const auto *hsi1_394 = buffer.data(hsi1 + 394);
    const auto *hsi1_395 = buffer.data(hsi1 + 395);
    const auto *hsi1_396 = buffer.data(hsi1 + 396);
    const auto *hsi1_397 = buffer.data(hsi1 + 397);
    const auto *hsi1_398 = buffer.data(hsi1 + 398);
    const auto *hsi1_399 = buffer.data(hsi1 + 399);
    const auto *hsi1_400 = buffer.data(hsi1 + 400);
    const auto *hsi1_401 = buffer.data(hsi1 + 401);
    const auto *hsi1_402 = buffer.data(hsi1 + 402);
    const auto *hsi1_403 = buffer.data(hsi1 + 403);
    const auto *hsi1_404 = buffer.data(hsi1 + 404);
    const auto *hsi1_405 = buffer.data(hsi1 + 405);
    const auto *hsi1_406 = buffer.data(hsi1 + 406);
    const auto *hsi1_420 = buffer.data(hsi1 + 420);
    const auto *hsi1_421 = buffer.data(hsi1 + 421);
    const auto *hsi1_423 = buffer.data(hsi1 + 423);
    const auto *hsi1_425 = buffer.data(hsi1 + 425);
    const auto *hsi1_426 = buffer.data(hsi1 + 426);
    const auto *hsi1_428 = buffer.data(hsi1 + 428);
    const auto *hsi1_429 = buffer.data(hsi1 + 429);
    const auto *hsi1_430 = buffer.data(hsi1 + 430);
    const auto *hsi1_432 = buffer.data(hsi1 + 432);
    const auto *hsi1_433 = buffer.data(hsi1 + 433);
    const auto *hsi1_434 = buffer.data(hsi1 + 434);
    const auto *hsi1_435 = buffer.data(hsi1 + 435);
    const auto *hsi1_437 = buffer.data(hsi1 + 437);
    const auto *hsi1_438 = buffer.data(hsi1 + 438);
    const auto *hsi1_439 = buffer.data(hsi1 + 439);
    const auto *hsi1_440 = buffer.data(hsi1 + 440);
    const auto *hsi1_441 = buffer.data(hsi1 + 441);
    const auto *hsi1_442 = buffer.data(hsi1 + 442);
    const auto *hsi1_443 = buffer.data(hsi1 + 443);
    const auto *hsi1_444 = buffer.data(hsi1 + 444);
    const auto *hsi1_445 = buffer.data(hsi1 + 445);
    const auto *hsi1_446 = buffer.data(hsi1 + 446);
    const auto *hsi1_447 = buffer.data(hsi1 + 447);

    const auto *hsk_477 = buffer.data(hsk + 477);
    const auto *hsk_478 = buffer.data(hsk + 478);
    const auto *hsk_482 = buffer.data(hsk + 482);
    const auto *hsk_483 = buffer.data(hsk + 483);
    const auto *hsk_488 = buffer.data(hsk + 488);
    const auto *hsk_496 = buffer.data(hsk + 496);
    const auto *hsk_497 = buffer.data(hsk + 497);
    const auto *hsk_498 = buffer.data(hsk + 498);
    const auto *hsk_499 = buffer.data(hsk + 499);
    const auto *hsk_500 = buffer.data(hsk + 500);
    const auto *hsk_501 = buffer.data(hsk + 501);
    const auto *hsk_502 = buffer.data(hsk + 502);
    const auto *hsk_503 = buffer.data(hsk + 503);
    const auto *hsk_504 = buffer.data(hsk + 504);
    const auto *hsk_505 = buffer.data(hsk + 505);
    const auto *hsk_506 = buffer.data(hsk + 506);
    const auto *hsk_507 = buffer.data(hsk + 507);
    const auto *hsk_508 = buffer.data(hsk + 508);
    const auto *hsk_509 = buffer.data(hsk + 509);
    const auto *hsk_510 = buffer.data(hsk + 510);
    const auto *hsk_511 = buffer.data(hsk + 511);
    const auto *hsk_512 = buffer.data(hsk + 512);
    const auto *hsk_513 = buffer.data(hsk + 513);
    const auto *hsk_514 = buffer.data(hsk + 514);
    const auto *hsk_515 = buffer.data(hsk + 515);
    const auto *hsk_516 = buffer.data(hsk + 516);
    const auto *hsk_517 = buffer.data(hsk + 517);
    const auto *hsk_518 = buffer.data(hsk + 518);
    const auto *hsk_519 = buffer.data(hsk + 519);
    const auto *hsk_520 = buffer.data(hsk + 520);
    const auto *hsk_521 = buffer.data(hsk + 521);
    const auto *hsk_522 = buffer.data(hsk + 522);
    const auto *hsk_523 = buffer.data(hsk + 523);
    const auto *hsk_524 = buffer.data(hsk + 524);
    const auto *hsk_531 = buffer.data(hsk + 531);
    const auto *hsk_532 = buffer.data(hsk + 532);
    const auto *hsk_533 = buffer.data(hsk + 533);
    const auto *hsk_534 = buffer.data(hsk + 534);
    const auto *hsk_535 = buffer.data(hsk + 535);
    const auto *hsk_536 = buffer.data(hsk + 536);
    const auto *hsk_537 = buffer.data(hsk + 537);
    const auto *hsk_539 = buffer.data(hsk + 539);
    const auto *hsk_540 = buffer.data(hsk + 540);
    const auto *hsk_541 = buffer.data(hsk + 541);
    const auto *hsk_543 = buffer.data(hsk + 543);
    const auto *hsk_545 = buffer.data(hsk + 545);
    const auto *hsk_546 = buffer.data(hsk + 546);
    const auto *hsk_548 = buffer.data(hsk + 548);
    const auto *hsk_549 = buffer.data(hsk + 549);
    const auto *hsk_550 = buffer.data(hsk + 550);
    const auto *hsk_552 = buffer.data(hsk + 552);
    const auto *hsk_553 = buffer.data(hsk + 553);
    const auto *hsk_554 = buffer.data(hsk + 554);
    const auto *hsk_555 = buffer.data(hsk + 555);
    const auto *hsk_557 = buffer.data(hsk + 557);
    const auto *hsk_558 = buffer.data(hsk + 558);
    const auto *hsk_559 = buffer.data(hsk + 559);
    const auto *hsk_560 = buffer.data(hsk + 560);
    const auto *hsk_561 = buffer.data(hsk + 561);
    const auto *hsk_563 = buffer.data(hsk + 563);
    const auto *hsk_564 = buffer.data(hsk + 564);
    const auto *hsk_565 = buffer.data(hsk + 565);
    const auto *hsk_566 = buffer.data(hsk + 566);
    const auto *hsk_567 = buffer.data(hsk + 567);
    const auto *hsk_568 = buffer.data(hsk + 568);
    const auto *hsk_569 = buffer.data(hsk + 569);
    const auto *hsk_570 = buffer.data(hsk + 570);
    const auto *hsk_571 = buffer.data(hsk + 571);
    const auto *hsk_572 = buffer.data(hsk + 572);
    const auto *hsk_573 = buffer.data(hsk + 573);
    const auto *hsk_574 = buffer.data(hsk + 574);
    const auto *hsk_575 = buffer.data(hsk + 575);

#pragma omp simd aligned(t_597, t_598, t_599, pa_x, pa_y, pc_x, pc_y, gsl0_419, gsl0_597, \
                         gsk_333, gsk_480, gsl1_419, gsl1_597, \
                         hsk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = pa_x[k] * gsl0_597[k]
                   + f_18 * gsk_480[k]
                   - f_14 * pc_x[k] * gsl1_597[k];

        t_598[k] = f_15 * gsk_333[k]
                   + f_3 * pc_y[k] * hsk_477[k];

        t_599[k] = pa_y[k] * gsl0_419[k]
                   - f_14 * pc_y[k] * gsl1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pa_x, pc_x, pc_z, gsl0_600, gsl0_602, gsk_298, \
                         gsk_483, gsk_485, gsl1_600, gsl1_602, \
                         hsk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pa_x[k] * gsl0_600[k]
                   + f_17 * gsk_483[k]
                   - f_14 * pc_x[k] * gsl1_600[k];

        t_601[k] = f_17 * gsk_298[k]
                   + f_3 * pc_z[k] * hsk_478[k];

        t_602[k] = pa_x[k] * gsl0_602[k]
                   + f_17 * gsk_485[k]
                   - f_14 * pc_x[k] * gsl1_602[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, pa_x, pa_y, pc_x, pc_y, gsl0_425, gsl0_603, \
                         gsk_338, gsk_486, gsl1_425, gsl1_603, \
                         hsk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pa_x[k] * gsl0_603[k]
                   + f_17 * gsk_486[k]
                   - f_14 * pc_x[k] * gsl1_603[k];

        t_604[k] = f_15 * gsk_338[k]
                   + f_3 * pc_y[k] * hsk_482[k];

        t_605[k] = pa_y[k] * gsl0_425[k]
                   - f_14 * pc_y[k] * gsl1_425[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, pa_x, pc_x, pc_z, gsl0_606, gsl0_608, gsk_303, \
                         gsk_489, gsk_491, gsl1_606, gsl1_608, \
                         hsk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = pa_x[k] * gsl0_606[k]
                   + f_16 * gsk_489[k]
                   - f_14 * pc_x[k] * gsl1_606[k];

        t_607[k] = f_17 * gsk_303[k]
                   + f_3 * pc_z[k] * hsk_483[k];

        t_608[k] = pa_x[k] * gsl0_608[k]
                   + f_16 * gsk_491[k]
                   - f_14 * pc_x[k] * gsl1_608[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pa_x, pc_x, pc_y, gsl0_609, gsl0_610, gsk_344, \
                         gsk_492, gsk_493, gsl1_609, gsl1_610, \
                         hsk_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = pa_x[k] * gsl0_609[k]
                   + f_16 * gsk_492[k]
                   - f_14 * pc_x[k] * gsl1_609[k];

        t_610[k] = pa_x[k] * gsl0_610[k]
                   + f_16 * gsk_493[k]
                   - f_14 * pc_x[k] * gsl1_610[k];

        t_611[k] = f_15 * gsk_344[k]
                   + f_3 * pc_y[k] * hsk_488[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pa_y, pc_x, pc_y, gsl0_432, gsk_496, \
                         gsk_497, gsk_498, gsl1_432, hsk_496, hsk_497, \
                         hsk_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pa_y[k] * gsl0_432[k]
                   - f_14 * pc_y[k] * gsl1_432[k];

        t_613[k] = f_15 * gsk_496[k]
                   + f_3 * pc_x[k] * hsk_496[k];

        t_614[k] = f_15 * gsk_497[k]
                   + f_3 * pc_x[k] * hsk_497[k];

        t_615[k] = f_15 * gsk_498[k]
                   + f_3 * pc_x[k] * hsk_498[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, t_620, pc_x, gsk_499, gsk_500, gsk_501, \
                         gsk_502, gsk_503, hsk_499, hsk_500, hsk_501, hsk_502, \
                         hsk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = f_15 * gsk_499[k]
                   + f_3 * pc_x[k] * hsk_499[k];

        t_617[k] = f_15 * gsk_500[k]
                   + f_3 * pc_x[k] * hsk_500[k];

        t_618[k] = f_15 * gsk_501[k]
                   + f_3 * pc_x[k] * hsk_501[k];

        t_619[k] = f_15 * gsk_502[k]
                   + f_3 * pc_x[k] * hsk_502[k];

        t_620[k] = f_15 * gsk_503[k]
                   + f_3 * pc_x[k] * hsk_503[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, pa_x, pc_x, pc_z, gsl0_621, gsl0_623, \
                         gsl0_624, gsk_316, gsl1_621, gsl1_623, gsl1_624, \
                         hsk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = pa_x[k] * gsl0_621[k]
                   - f_14 * pc_x[k] * gsl1_621[k];

        t_622[k] = f_17 * gsk_316[k]
                   + f_3 * pc_z[k] * hsk_496[k];

        t_623[k] = pa_x[k] * gsl0_623[k]
                   - f_14 * pc_x[k] * gsl1_623[k];

        t_624[k] = pa_x[k] * gsl0_624[k]
                   - f_14 * pc_x[k] * gsl1_624[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, pa_x, pc_x, pc_y, gsl0_625, gsl0_626, \
                         gsl0_627, gsk_359, gsl1_625, gsl1_626, gsl1_627, \
                         hsk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = pa_x[k] * gsl0_625[k]
                   - f_14 * pc_x[k] * gsl1_625[k];

        t_626[k] = pa_x[k] * gsl0_626[k]
                   - f_14 * pc_x[k] * gsl1_626[k];

        t_627[k] = pa_x[k] * gsl0_627[k]
                   - f_14 * pc_x[k] * gsl1_627[k];

        t_628[k] = f_15 * gsk_359[k]
                   + f_3 * pc_y[k] * hsk_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pa_x, pc_x, pc_y, pc_z, gsl0_629, \
                         gsl0_630, gsk_324, gsk_504, gsl1_629, gsl1_630, \
                         hsk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pa_x[k] * gsl0_629[k]
                   - f_14 * pc_x[k] * gsl1_629[k];

        t_630[k] = pa_x[k] * gsl0_630[k]
                   + f_22 * gsk_504[k]
                   - f_14 * pc_x[k] * gsl1_630[k];

        t_631[k] = f_3 * pc_y[k] * hsk_504[k];

        t_632[k] = f_18 * gsk_324[k]
                   + f_3 * pc_z[k] * hsk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pa_x, pc_x, pc_y, gsl0_635, gsk_509, gsl1_635, \
                         hsi0_392, hsi1_392, hsk_505, hsk_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_4 * hsi0_392[k]
                   - f_5 * hsi1_392[k]
                   + f_3 * pc_y[k] * hsk_505[k];

        t_634[k] = f_3 * pc_y[k] * hsk_506[k];

        t_635[k] = pa_x[k] * gsl0_635[k]
                   + f_19 * gsk_509[k]
                   - f_14 * pc_x[k] * gsl1_635[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_y, hsi0_393, hsi0_394, hsi1_393, hsi1_394, \
                         hsk_507, hsk_508, hsk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_6 * hsi0_393[k]
                   - f_7 * hsi1_393[k]
                   + f_3 * pc_y[k] * hsk_507[k];

        t_637[k] = f_4 * hsi0_394[k]
                   - f_5 * hsi1_394[k]
                   + f_3 * pc_y[k] * hsk_508[k];

        t_638[k] = f_3 * pc_y[k] * hsk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pa_x, pc_x, pc_y, gsl0_639, gsk_513, gsl1_639, \
                         hsi0_395, hsi0_396, hsi1_395, hsi1_396, hsk_510, \
                         hsk_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pa_x[k] * gsl0_639[k]
                   + f_0 * gsk_513[k]
                   - f_14 * pc_x[k] * gsl1_639[k];

        t_640[k] = f_8 * hsi0_395[k]
                   - f_9 * hsi1_395[k]
                   + f_3 * pc_y[k] * hsk_510[k];

        t_641[k] = f_6 * hsi0_396[k]
                   - f_7 * hsi1_396[k]
                   + f_3 * pc_y[k] * hsk_511[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pa_x, pc_x, pc_y, gsl0_644, gsk_518, gsl1_644, \
                         hsi0_397, hsi1_397, hsk_512, hsk_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_4 * hsi0_397[k]
                   - f_5 * hsi1_397[k]
                   + f_3 * pc_y[k] * hsk_512[k];

        t_643[k] = f_3 * pc_y[k] * hsk_513[k];

        t_644[k] = pa_x[k] * gsl0_644[k]
                   + f_18 * gsk_518[k]
                   - f_14 * pc_x[k] * gsl1_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_y, hsi0_398, hsi0_399, hsi0_400, hsi1_398, \
                         hsi1_399, hsi1_400, hsk_514, hsk_515, \
                         hsk_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_10 * hsi0_398[k]
                   - f_11 * hsi1_398[k]
                   + f_3 * pc_y[k] * hsk_514[k];

        t_646[k] = f_8 * hsi0_399[k]
                   - f_9 * hsi1_399[k]
                   + f_3 * pc_y[k] * hsk_515[k];

        t_647[k] = f_6 * hsi0_400[k]
                   - f_7 * hsi1_400[k]
                   + f_3 * pc_y[k] * hsk_516[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pa_x, pc_x, pc_y, gsl0_650, gsk_524, gsl1_650, \
                         hsi0_401, hsi1_401, hsk_517, hsk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_4 * hsi0_401[k]
                   - f_5 * hsi1_401[k]
                   + f_3 * pc_y[k] * hsk_517[k];

        t_649[k] = f_3 * pc_y[k] * hsk_518[k];

        t_650[k] = pa_x[k] * gsl0_650[k]
                   + f_17 * gsk_524[k]
                   - f_14 * pc_x[k] * gsl1_650[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_y, hsi0_402, hsi0_403, hsi0_404, hsi1_402, \
                         hsi1_403, hsi1_404, hsk_519, hsk_520, \
                         hsk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_12 * hsi0_402[k]
                   - f_13 * hsi1_402[k]
                   + f_3 * pc_y[k] * hsk_519[k];

        t_652[k] = f_10 * hsi0_403[k]
                   - f_11 * hsi1_403[k]
                   + f_3 * pc_y[k] * hsk_520[k];

        t_653[k] = f_8 * hsi0_404[k]
                   - f_9 * hsi1_404[k]
                   + f_3 * pc_y[k] * hsk_521[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_y, hsi0_405, hsi0_406, hsi1_405, hsi1_406, \
                         hsk_522, hsk_523, hsk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_6 * hsi0_405[k]
                   - f_7 * hsi1_405[k]
                   + f_3 * pc_y[k] * hsk_522[k];

        t_655[k] = f_4 * hsi0_406[k]
                   - f_5 * hsi1_406[k]
                   + f_3 * pc_y[k] * hsk_523[k];

        t_656[k] = f_3 * pc_y[k] * hsk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pa_x, pc_x, gsl0_657, gsk_531, gsk_532, \
                         gsk_533, gsk_534, gsl1_657, hsk_532, hsk_533, \
                         hsk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = pa_x[k] * gsl0_657[k]
                   + f_16 * gsk_531[k]
                   - f_14 * pc_x[k] * gsl1_657[k];

        t_658[k] = f_15 * gsk_532[k]
                   + f_3 * pc_x[k] * hsk_532[k];

        t_659[k] = f_15 * gsk_533[k]
                   + f_3 * pc_x[k] * hsk_533[k];

        t_660[k] = f_15 * gsk_534[k]
                   + f_3 * pc_x[k] * hsk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, pc_y, gsk_535, gsk_536, \
                         gsk_537, gsk_539, hsk_531, hsk_535, hsk_536, hsk_537, \
                         hsk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_15 * gsk_535[k]
                   + f_3 * pc_x[k] * hsk_535[k];

        t_662[k] = f_15 * gsk_536[k]
                   + f_3 * pc_x[k] * hsk_536[k];

        t_663[k] = f_15 * gsk_537[k]
                   + f_3 * pc_x[k] * hsk_537[k];

        t_664[k] = f_3 * pc_y[k] * hsk_531[k];

        t_665[k] = f_15 * gsk_539[k]
                   + f_3 * pc_x[k] * hsk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pa_x, pc_x, gsl0_666, gsl0_667, gsl0_668, \
                         gsl0_669, gsl1_666, gsl1_667, gsl1_668, \
                         gsl1_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = pa_x[k] * gsl0_666[k]
                   - f_14 * pc_x[k] * gsl1_666[k];

        t_667[k] = pa_x[k] * gsl0_667[k]
                   - f_14 * pc_x[k] * gsl1_667[k];

        t_668[k] = pa_x[k] * gsl0_668[k]
                   - f_14 * pc_x[k] * gsl1_668[k];

        t_669[k] = pa_x[k] * gsl0_669[k]
                   - f_14 * pc_x[k] * gsl1_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pa_x, pc_x, pc_y, gsl0_670, gsl0_671, \
                         gsl0_672, gsl1_670, gsl1_671, gsl1_672, \
                         hsk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = pa_x[k] * gsl0_670[k]
                   - f_14 * pc_x[k] * gsl1_670[k];

        t_671[k] = pa_x[k] * gsl0_671[k]
                   - f_14 * pc_x[k] * gsl1_671[k];

        t_672[k] = pa_x[k] * gsl0_672[k]
                   - f_14 * pc_x[k] * gsl1_672[k];

        t_673[k] = f_3 * pc_y[k] * hsk_539[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pa_x, pc_x, pc_z, gsl0_674, gsl1_674, \
                         hsi0_420, hsi0_421, hsi1_420, hsi1_421, hsk_540, \
                         hsk_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = pa_x[k] * gsl0_674[k]
                   - f_14 * pc_x[k] * gsl1_674[k];

        t_675[k] = f_1 * hsi0_420[k]
                   - f_2 * hsi1_420[k]
                   + f_3 * pc_x[k] * hsk_540[k];

        t_676[k] = f_20 * hsi0_421[k]
                   - f_21 * hsi1_421[k]
                   + f_3 * pc_x[k] * hsk_541[k];

        t_677[k] = f_3 * pc_z[k] * hsk_540[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, t_681, pc_x, pc_z, hsi0_423, hsi0_425, hsi0_426, \
                         hsi1_423, hsi1_425, hsi1_426, hsk_541, hsk_543, hsk_545, \
                         hsk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_12 * hsi0_423[k]
                   - f_13 * hsi1_423[k]
                   + f_3 * pc_x[k] * hsk_543[k];

        t_679[k] = f_3 * pc_z[k] * hsk_541[k];

        t_680[k] = f_12 * hsi0_425[k]
                   - f_13 * hsi1_425[k]
                   + f_3 * pc_x[k] * hsk_545[k];

        t_681[k] = f_10 * hsi0_426[k]
                   - f_11 * hsi1_426[k]
                   + f_3 * pc_x[k] * hsk_546[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pc_x, pc_z, hsi0_428, hsi0_429, hsi0_430, \
                         hsi1_428, hsi1_429, hsi1_430, hsk_543, hsk_548, hsk_549, \
                         hsk_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_3 * pc_z[k] * hsk_543[k];

        t_683[k] = f_10 * hsi0_428[k]
                   - f_11 * hsi1_428[k]
                   + f_3 * pc_x[k] * hsk_548[k];

        t_684[k] = f_10 * hsi0_429[k]
                   - f_11 * hsi1_429[k]
                   + f_3 * pc_x[k] * hsk_549[k];

        t_685[k] = f_8 * hsi0_430[k]
                   - f_9 * hsi1_430[k]
                   + f_3 * pc_x[k] * hsk_550[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, pc_z, hsi0_432, hsi0_433, hsi0_434, \
                         hsi1_432, hsi1_433, hsi1_434, hsk_546, hsk_552, hsk_553, \
                         hsk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_3 * pc_z[k] * hsk_546[k];

        t_687[k] = f_8 * hsi0_432[k]
                   - f_9 * hsi1_432[k]
                   + f_3 * pc_x[k] * hsk_552[k];

        t_688[k] = f_8 * hsi0_433[k]
                   - f_9 * hsi1_433[k]
                   + f_3 * pc_x[k] * hsk_553[k];

        t_689[k] = f_8 * hsi0_434[k]
                   - f_9 * hsi1_434[k]
                   + f_3 * pc_x[k] * hsk_554[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_z, hsi0_435, hsi0_437, hsi0_438, \
                         hsi1_435, hsi1_437, hsi1_438, hsk_550, hsk_555, hsk_557, \
                         hsk_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_6 * hsi0_435[k]
                   - f_7 * hsi1_435[k]
                   + f_3 * pc_x[k] * hsk_555[k];

        t_691[k] = f_3 * pc_z[k] * hsk_550[k];

        t_692[k] = f_6 * hsi0_437[k]
                   - f_7 * hsi1_437[k]
                   + f_3 * pc_x[k] * hsk_557[k];

        t_693[k] = f_6 * hsi0_438[k]
                   - f_7 * hsi1_438[k]
                   + f_3 * pc_x[k] * hsk_558[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pc_x, pc_z, hsi0_439, hsi0_440, hsi0_441, \
                         hsi1_439, hsi1_440, hsi1_441, hsk_555, hsk_559, hsk_560, \
                         hsk_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_6 * hsi0_439[k]
                   - f_7 * hsi1_439[k]
                   + f_3 * pc_x[k] * hsk_559[k];

        t_695[k] = f_6 * hsi0_440[k]
                   - f_7 * hsi1_440[k]
                   + f_3 * pc_x[k] * hsk_560[k];

        t_696[k] = f_4 * hsi0_441[k]
                   - f_5 * hsi1_441[k]
                   + f_3 * pc_x[k] * hsk_561[k];

        t_697[k] = f_3 * pc_z[k] * hsk_555[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, hsi0_443, hsi0_444, hsi0_445, hsi1_443, \
                         hsi1_444, hsi1_445, hsk_563, hsk_564, \
                         hsk_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_4 * hsi0_443[k]
                   - f_5 * hsi1_443[k]
                   + f_3 * pc_x[k] * hsk_563[k];

        t_699[k] = f_4 * hsi0_444[k]
                   - f_5 * hsi1_444[k]
                   + f_3 * pc_x[k] * hsk_564[k];

        t_700[k] = f_4 * hsi0_445[k]
                   - f_5 * hsi1_445[k]
                   + f_3 * pc_x[k] * hsk_565[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, t_705, pc_x, hsi0_446, hsi0_447, \
                         hsi1_446, hsi1_447, hsk_566, hsk_567, hsk_568, hsk_569, \
                         hsk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_4 * hsi0_446[k]
                   - f_5 * hsi1_446[k]
                   + f_3 * pc_x[k] * hsk_566[k];

        t_702[k] = f_4 * hsi0_447[k]
                   - f_5 * hsi1_447[k]
                   + f_3 * pc_x[k] * hsk_567[k];

        t_703[k] = f_3 * pc_x[k] * hsk_568[k];

        t_704[k] = f_3 * pc_x[k] * hsk_569[k];

        t_705[k] = f_3 * pc_x[k] * hsk_570[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, t_710, pc_x, hsk_571, hsk_572, hsk_573, \
                         hsk_574, hsk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_3 * pc_x[k] * hsk_571[k];

        t_707[k] = f_3 * pc_x[k] * hsk_572[k];

        t_708[k] = f_3 * pc_x[k] * hsk_573[k];

        t_709[k] = f_3 * pc_x[k] * hsk_574[k];

        t_710[k] = f_3 * pc_x[k] * hsk_575[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, t_714, pc_y, pc_z, gsk_388, hsi0_441, hsi0_442, \
                         hsi1_441, hsi1_442, hsk_568, hsk_569, \
                         hsk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_0 * gsk_388[k]
                   + f_1 * hsi0_441[k]
                   - f_2 * hsi1_441[k]
                   + f_3 * pc_y[k] * hsk_568[k];

        t_712[k] = f_3 * pc_z[k] * hsk_568[k];

        t_713[k] = f_4 * hsi0_441[k]
                   - f_5 * hsi1_441[k]
                   + f_3 * pc_z[k] * hsk_569[k];

        t_714[k] = f_6 * hsi0_442[k]
                   - f_7 * hsi1_442[k]
                   + f_3 * pc_z[k] * hsk_570[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, pc_z, hsi0_443, hsi0_444, hsi0_445, hsi1_443, \
                         hsi1_444, hsi1_445, hsk_571, hsk_572, \
                         hsk_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_8 * hsi0_443[k]
                   - f_9 * hsi1_443[k]
                   + f_3 * pc_z[k] * hsk_571[k];

        t_716[k] = f_10 * hsi0_444[k]
                   - f_11 * hsi1_444[k]
                   + f_3 * pc_z[k] * hsk_572[k];

        t_717[k] = f_12 * hsi0_445[k]
                   - f_13 * hsi1_445[k]
                   + f_3 * pc_z[k] * hsk_573[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pa_z, pc_y, pc_z, gsl0_450, gsl0_451, \
                         gsk_395, gsl1_450, gsl1_451, hsi0_447, hsi1_447, \
                         hsk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_0 * gsk_395[k]
                   + f_3 * pc_y[k] * hsk_575[k];

        t_719[k] = f_1 * hsi0_447[k]
                   - f_2 * hsi1_447[k]
                   + f_3 * pc_z[k] * hsk_575[k];

        t_720[k] = pa_z[k] * gsl0_450[k]
                   - f_14 * pc_z[k] * gsl1_450[k];

        t_721[k] = pa_z[k] * gsl0_451[k]
                   - f_14 * pc_z[k] * gsl1_451[k];
    }
}

static auto
compute_prim_hsl_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsl0,
                                                          const size_t gsk, const size_t gsl1,
                                                          const size_t hsi0, const size_t hsi1,
                                                          const size_t hsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsl0_453 = buffer.data(gsl0 + 453);
    const auto *gsl0_456 = buffer.data(gsl0 + 456);
    const auto *gsl0_460 = buffer.data(gsl0 + 460);
    const auto *gsl0_465 = buffer.data(gsl0 + 465);
    const auto *gsl0_471 = buffer.data(gsl0 + 471);
    const auto *gsl0_486 = buffer.data(gsl0 + 486);
    const auto *gsl0_488 = buffer.data(gsl0 + 488);
    const auto *gsl0_489 = buffer.data(gsl0 + 489);
    const auto *gsl0_490 = buffer.data(gsl0 + 490);
    const auto *gsl0_491 = buffer.data(gsl0 + 491);
    const auto *gsl0_492 = buffer.data(gsl0 + 492);

    const auto *gsk_388 = buffer.data(gsk + 388);
    const auto *gsk_389 = buffer.data(gsk + 389);
    const auto *gsk_390 = buffer.data(gsk + 390);
    const auto *gsk_391 = buffer.data(gsk + 391);
    const auto *gsk_392 = buffer.data(gsk + 392);
    const auto *gsk_393 = buffer.data(gsk + 393);
    const auto *gsk_395 = buffer.data(gsk + 395);
    const auto *gsk_424 = buffer.data(gsk + 424);
    const auto *gsk_431 = buffer.data(gsk + 431);
    const auto *gsk_460 = buffer.data(gsk + 460);
    const auto *gsk_462 = buffer.data(gsk + 462);
    const auto *gsk_463 = buffer.data(gsk + 463);
    const auto *gsk_464 = buffer.data(gsk + 464);
    const auto *gsk_465 = buffer.data(gsk + 465);
    const auto *gsk_466 = buffer.data(gsk + 466);
    const auto *gsk_467 = buffer.data(gsk + 467);

    const auto *gsl1_453 = buffer.data(gsl1 + 453);
    const auto *gsl1_456 = buffer.data(gsl1 + 456);
    const auto *gsl1_460 = buffer.data(gsl1 + 460);
    const auto *gsl1_465 = buffer.data(gsl1 + 465);
    const auto *gsl1_471 = buffer.data(gsl1 + 471);
    const auto *gsl1_486 = buffer.data(gsl1 + 486);
    const auto *gsl1_488 = buffer.data(gsl1 + 488);
    const auto *gsl1_489 = buffer.data(gsl1 + 489);
    const auto *gsl1_490 = buffer.data(gsl1 + 490);
    const auto *gsl1_491 = buffer.data(gsl1 + 491);
    const auto *gsl1_492 = buffer.data(gsl1 + 492);

    const auto *hsi0_450 = buffer.data(hsi0 + 450);
    const auto *hsi0_452 = buffer.data(hsi0 + 452);
    const auto *hsi0_453 = buffer.data(hsi0 + 453);
    const auto *hsi0_455 = buffer.data(hsi0 + 455);
    const auto *hsi0_456 = buffer.data(hsi0 + 456);
    const auto *hsi0_457 = buffer.data(hsi0 + 457);
    const auto *hsi0_459 = buffer.data(hsi0 + 459);
    const auto *hsi0_460 = buffer.data(hsi0 + 460);
    const auto *hsi0_461 = buffer.data(hsi0 + 461);
    const auto *hsi0_462 = buffer.data(hsi0 + 462);
    const auto *hsi0_464 = buffer.data(hsi0 + 464);
    const auto *hsi0_465 = buffer.data(hsi0 + 465);
    const auto *hsi0_466 = buffer.data(hsi0 + 466);
    const auto *hsi0_467 = buffer.data(hsi0 + 467);
    const auto *hsi0_468 = buffer.data(hsi0 + 468);
    const auto *hsi0_470 = buffer.data(hsi0 + 470);
    const auto *hsi0_471 = buffer.data(hsi0 + 471);
    const auto *hsi0_472 = buffer.data(hsi0 + 472);
    const auto *hsi0_473 = buffer.data(hsi0 + 473);
    const auto *hsi0_474 = buffer.data(hsi0 + 474);
    const auto *hsi0_475 = buffer.data(hsi0 + 475);
    const auto *hsi0_476 = buffer.data(hsi0 + 476);
    const auto *hsi0_477 = buffer.data(hsi0 + 477);
    const auto *hsi0_478 = buffer.data(hsi0 + 478);
    const auto *hsi0_479 = buffer.data(hsi0 + 479);
    const auto *hsi0_480 = buffer.data(hsi0 + 480);
    const auto *hsi0_481 = buffer.data(hsi0 + 481);
    const auto *hsi0_482 = buffer.data(hsi0 + 482);
    const auto *hsi0_483 = buffer.data(hsi0 + 483);
    const auto *hsi0_484 = buffer.data(hsi0 + 484);
    const auto *hsi0_485 = buffer.data(hsi0 + 485);
    const auto *hsi0_486 = buffer.data(hsi0 + 486);
    const auto *hsi0_487 = buffer.data(hsi0 + 487);
    const auto *hsi0_488 = buffer.data(hsi0 + 488);
    const auto *hsi0_489 = buffer.data(hsi0 + 489);
    const auto *hsi0_490 = buffer.data(hsi0 + 490);
    const auto *hsi0_491 = buffer.data(hsi0 + 491);
    const auto *hsi0_492 = buffer.data(hsi0 + 492);
    const auto *hsi0_493 = buffer.data(hsi0 + 493);
    const auto *hsi0_494 = buffer.data(hsi0 + 494);
    const auto *hsi0_495 = buffer.data(hsi0 + 495);
    const auto *hsi0_496 = buffer.data(hsi0 + 496);
    const auto *hsi0_497 = buffer.data(hsi0 + 497);
    const auto *hsi0_498 = buffer.data(hsi0 + 498);
    const auto *hsi0_499 = buffer.data(hsi0 + 499);
    const auto *hsi0_500 = buffer.data(hsi0 + 500);
    const auto *hsi0_501 = buffer.data(hsi0 + 501);
    const auto *hsi0_502 = buffer.data(hsi0 + 502);
    const auto *hsi0_503 = buffer.data(hsi0 + 503);
    const auto *hsi0_504 = buffer.data(hsi0 + 504);
    const auto *hsi0_505 = buffer.data(hsi0 + 505);
    const auto *hsi0_506 = buffer.data(hsi0 + 506);
    const auto *hsi0_507 = buffer.data(hsi0 + 507);
    const auto *hsi0_508 = buffer.data(hsi0 + 508);
    const auto *hsi0_509 = buffer.data(hsi0 + 509);
    const auto *hsi0_510 = buffer.data(hsi0 + 510);
    const auto *hsi0_511 = buffer.data(hsi0 + 511);
    const auto *hsi0_512 = buffer.data(hsi0 + 512);
    const auto *hsi0_513 = buffer.data(hsi0 + 513);
    const auto *hsi0_514 = buffer.data(hsi0 + 514);
    const auto *hsi0_515 = buffer.data(hsi0 + 515);
    const auto *hsi0_516 = buffer.data(hsi0 + 516);
    const auto *hsi0_517 = buffer.data(hsi0 + 517);
    const auto *hsi0_518 = buffer.data(hsi0 + 518);
    const auto *hsi0_519 = buffer.data(hsi0 + 519);
    const auto *hsi0_520 = buffer.data(hsi0 + 520);
    const auto *hsi0_521 = buffer.data(hsi0 + 521);
    const auto *hsi0_522 = buffer.data(hsi0 + 522);
    const auto *hsi0_523 = buffer.data(hsi0 + 523);
    const auto *hsi0_524 = buffer.data(hsi0 + 524);
    const auto *hsi0_525 = buffer.data(hsi0 + 525);
    const auto *hsi0_526 = buffer.data(hsi0 + 526);
    const auto *hsi0_527 = buffer.data(hsi0 + 527);

    const auto *hsi1_450 = buffer.data(hsi1 + 450);
    const auto *hsi1_452 = buffer.data(hsi1 + 452);
    const auto *hsi1_453 = buffer.data(hsi1 + 453);
    const auto *hsi1_455 = buffer.data(hsi1 + 455);
    const auto *hsi1_456 = buffer.data(hsi1 + 456);
    const auto *hsi1_457 = buffer.data(hsi1 + 457);
    const auto *hsi1_459 = buffer.data(hsi1 + 459);
    const auto *hsi1_460 = buffer.data(hsi1 + 460);
    const auto *hsi1_461 = buffer.data(hsi1 + 461);
    const auto *hsi1_462 = buffer.data(hsi1 + 462);
    const auto *hsi1_464 = buffer.data(hsi1 + 464);
    const auto *hsi1_465 = buffer.data(hsi1 + 465);
    const auto *hsi1_466 = buffer.data(hsi1 + 466);
    const auto *hsi1_467 = buffer.data(hsi1 + 467);
    const auto *hsi1_468 = buffer.data(hsi1 + 468);
    const auto *hsi1_470 = buffer.data(hsi1 + 470);
    const auto *hsi1_471 = buffer.data(hsi1 + 471);
    const auto *hsi1_472 = buffer.data(hsi1 + 472);
    const auto *hsi1_473 = buffer.data(hsi1 + 473);
    const auto *hsi1_474 = buffer.data(hsi1 + 474);
    const auto *hsi1_475 = buffer.data(hsi1 + 475);
    const auto *hsi1_476 = buffer.data(hsi1 + 476);
    const auto *hsi1_477 = buffer.data(hsi1 + 477);
    const auto *hsi1_478 = buffer.data(hsi1 + 478);
    const auto *hsi1_479 = buffer.data(hsi1 + 479);
    const auto *hsi1_480 = buffer.data(hsi1 + 480);
    const auto *hsi1_481 = buffer.data(hsi1 + 481);
    const auto *hsi1_482 = buffer.data(hsi1 + 482);
    const auto *hsi1_483 = buffer.data(hsi1 + 483);
    const auto *hsi1_484 = buffer.data(hsi1 + 484);
    const auto *hsi1_485 = buffer.data(hsi1 + 485);
    const auto *hsi1_486 = buffer.data(hsi1 + 486);
    const auto *hsi1_487 = buffer.data(hsi1 + 487);
    const auto *hsi1_488 = buffer.data(hsi1 + 488);
    const auto *hsi1_489 = buffer.data(hsi1 + 489);
    const auto *hsi1_490 = buffer.data(hsi1 + 490);
    const auto *hsi1_491 = buffer.data(hsi1 + 491);
    const auto *hsi1_492 = buffer.data(hsi1 + 492);
    const auto *hsi1_493 = buffer.data(hsi1 + 493);
    const auto *hsi1_494 = buffer.data(hsi1 + 494);
    const auto *hsi1_495 = buffer.data(hsi1 + 495);
    const auto *hsi1_496 = buffer.data(hsi1 + 496);
    const auto *hsi1_497 = buffer.data(hsi1 + 497);
    const auto *hsi1_498 = buffer.data(hsi1 + 498);
    const auto *hsi1_499 = buffer.data(hsi1 + 499);
    const auto *hsi1_500 = buffer.data(hsi1 + 500);
    const auto *hsi1_501 = buffer.data(hsi1 + 501);
    const auto *hsi1_502 = buffer.data(hsi1 + 502);
    const auto *hsi1_503 = buffer.data(hsi1 + 503);
    const auto *hsi1_504 = buffer.data(hsi1 + 504);
    const auto *hsi1_505 = buffer.data(hsi1 + 505);
    const auto *hsi1_506 = buffer.data(hsi1 + 506);
    const auto *hsi1_507 = buffer.data(hsi1 + 507);
    const auto *hsi1_508 = buffer.data(hsi1 + 508);
    const auto *hsi1_509 = buffer.data(hsi1 + 509);
    const auto *hsi1_510 = buffer.data(hsi1 + 510);
    const auto *hsi1_511 = buffer.data(hsi1 + 511);
    const auto *hsi1_512 = buffer.data(hsi1 + 512);
    const auto *hsi1_513 = buffer.data(hsi1 + 513);
    const auto *hsi1_514 = buffer.data(hsi1 + 514);
    const auto *hsi1_515 = buffer.data(hsi1 + 515);
    const auto *hsi1_516 = buffer.data(hsi1 + 516);
    const auto *hsi1_517 = buffer.data(hsi1 + 517);
    const auto *hsi1_518 = buffer.data(hsi1 + 518);
    const auto *hsi1_519 = buffer.data(hsi1 + 519);
    const auto *hsi1_520 = buffer.data(hsi1 + 520);
    const auto *hsi1_521 = buffer.data(hsi1 + 521);
    const auto *hsi1_522 = buffer.data(hsi1 + 522);
    const auto *hsi1_523 = buffer.data(hsi1 + 523);
    const auto *hsi1_524 = buffer.data(hsi1 + 524);
    const auto *hsi1_525 = buffer.data(hsi1 + 525);
    const auto *hsi1_526 = buffer.data(hsi1 + 526);
    const auto *hsi1_527 = buffer.data(hsi1 + 527);

    const auto *hsk_578 = buffer.data(hsk + 578);
    const auto *hsk_580 = buffer.data(hsk + 580);
    const auto *hsk_581 = buffer.data(hsk + 581);
    const auto *hsk_583 = buffer.data(hsk + 583);
    const auto *hsk_584 = buffer.data(hsk + 584);
    const auto *hsk_585 = buffer.data(hsk + 585);
    const auto *hsk_587 = buffer.data(hsk + 587);
    const auto *hsk_588 = buffer.data(hsk + 588);
    const auto *hsk_589 = buffer.data(hsk + 589);
    const auto *hsk_590 = buffer.data(hsk + 590);
    const auto *hsk_592 = buffer.data(hsk + 592);
    const auto *hsk_593 = buffer.data(hsk + 593);
    const auto *hsk_594 = buffer.data(hsk + 594);
    const auto *hsk_595 = buffer.data(hsk + 595);
    const auto *hsk_596 = buffer.data(hsk + 596);
    const auto *hsk_598 = buffer.data(hsk + 598);
    const auto *hsk_599 = buffer.data(hsk + 599);
    const auto *hsk_600 = buffer.data(hsk + 600);
    const auto *hsk_601 = buffer.data(hsk + 601);
    const auto *hsk_602 = buffer.data(hsk + 602);
    const auto *hsk_603 = buffer.data(hsk + 603);
    const auto *hsk_604 = buffer.data(hsk + 604);
    const auto *hsk_605 = buffer.data(hsk + 605);
    const auto *hsk_606 = buffer.data(hsk + 606);
    const auto *hsk_607 = buffer.data(hsk + 607);
    const auto *hsk_608 = buffer.data(hsk + 608);
    const auto *hsk_609 = buffer.data(hsk + 609);
    const auto *hsk_610 = buffer.data(hsk + 610);
    const auto *hsk_611 = buffer.data(hsk + 611);
    const auto *hsk_612 = buffer.data(hsk + 612);
    const auto *hsk_613 = buffer.data(hsk + 613);
    const auto *hsk_614 = buffer.data(hsk + 614);
    const auto *hsk_615 = buffer.data(hsk + 615);
    const auto *hsk_616 = buffer.data(hsk + 616);
    const auto *hsk_617 = buffer.data(hsk + 617);
    const auto *hsk_618 = buffer.data(hsk + 618);
    const auto *hsk_619 = buffer.data(hsk + 619);
    const auto *hsk_620 = buffer.data(hsk + 620);
    const auto *hsk_621 = buffer.data(hsk + 621);
    const auto *hsk_622 = buffer.data(hsk + 622);
    const auto *hsk_623 = buffer.data(hsk + 623);
    const auto *hsk_624 = buffer.data(hsk + 624);
    const auto *hsk_625 = buffer.data(hsk + 625);
    const auto *hsk_626 = buffer.data(hsk + 626);
    const auto *hsk_627 = buffer.data(hsk + 627);
    const auto *hsk_628 = buffer.data(hsk + 628);
    const auto *hsk_629 = buffer.data(hsk + 629);
    const auto *hsk_630 = buffer.data(hsk + 630);
    const auto *hsk_631 = buffer.data(hsk + 631);
    const auto *hsk_632 = buffer.data(hsk + 632);
    const auto *hsk_633 = buffer.data(hsk + 633);
    const auto *hsk_634 = buffer.data(hsk + 634);
    const auto *hsk_635 = buffer.data(hsk + 635);
    const auto *hsk_636 = buffer.data(hsk + 636);
    const auto *hsk_637 = buffer.data(hsk + 637);
    const auto *hsk_638 = buffer.data(hsk + 638);
    const auto *hsk_639 = buffer.data(hsk + 639);
    const auto *hsk_640 = buffer.data(hsk + 640);
    const auto *hsk_641 = buffer.data(hsk + 641);
    const auto *hsk_642 = buffer.data(hsk + 642);
    const auto *hsk_643 = buffer.data(hsk + 643);
    const auto *hsk_644 = buffer.data(hsk + 644);
    const auto *hsk_645 = buffer.data(hsk + 645);
    const auto *hsk_646 = buffer.data(hsk + 646);
    const auto *hsk_647 = buffer.data(hsk + 647);
    const auto *hsk_648 = buffer.data(hsk + 648);
    const auto *hsk_649 = buffer.data(hsk + 649);
    const auto *hsk_650 = buffer.data(hsk + 650);
    const auto *hsk_651 = buffer.data(hsk + 651);
    const auto *hsk_652 = buffer.data(hsk + 652);
    const auto *hsk_653 = buffer.data(hsk + 653);
    const auto *hsk_654 = buffer.data(hsk + 654);
    const auto *hsk_655 = buffer.data(hsk + 655);
    const auto *hsk_656 = buffer.data(hsk + 656);
    const auto *hsk_657 = buffer.data(hsk + 657);
    const auto *hsk_658 = buffer.data(hsk + 658);
    const auto *hsk_659 = buffer.data(hsk + 659);
    const auto *hsk_660 = buffer.data(hsk + 660);
    const auto *hsk_661 = buffer.data(hsk + 661);
    const auto *hsk_662 = buffer.data(hsk + 662);
    const auto *hsk_663 = buffer.data(hsk + 663);
    const auto *hsk_664 = buffer.data(hsk + 664);
    const auto *hsk_665 = buffer.data(hsk + 665);
    const auto *hsk_666 = buffer.data(hsk + 666);
    const auto *hsk_667 = buffer.data(hsk + 667);
    const auto *hsk_668 = buffer.data(hsk + 668);
    const auto *hsk_669 = buffer.data(hsk + 669);
    const auto *hsk_670 = buffer.data(hsk + 670);
    const auto *hsk_671 = buffer.data(hsk + 671);

#pragma omp simd aligned(t_722, t_723, t_724, pa_z, pc_x, pc_z, gsl0_453, gsl1_453, hsi0_450, \
                         hsi0_452, hsi1_450, hsi1_452, hsk_578, \
                         hsk_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_20 * hsi0_450[k]
                   - f_21 * hsi1_450[k]
                   + f_3 * pc_x[k] * hsk_578[k];

        t_723[k] = pa_z[k] * gsl0_453[k]
                   - f_14 * pc_z[k] * gsl1_453[k];

        t_724[k] = f_12 * hsi0_452[k]
                   - f_13 * hsi1_452[k]
                   + f_3 * pc_x[k] * hsk_580[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, pa_z, pc_x, pc_z, gsl0_456, gsl1_456, hsi0_453, \
                         hsi0_455, hsi1_453, hsi1_455, hsk_581, \
                         hsk_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_12 * hsi0_453[k]
                   - f_13 * hsi1_453[k]
                   + f_3 * pc_x[k] * hsk_581[k];

        t_726[k] = pa_z[k] * gsl0_456[k]
                   - f_14 * pc_z[k] * gsl1_456[k];

        t_727[k] = f_10 * hsi0_455[k]
                   - f_11 * hsi1_455[k]
                   + f_3 * pc_x[k] * hsk_583[k];
    }

#pragma omp simd aligned(t_728, t_729, t_730, pa_z, pc_x, pc_z, gsl0_460, gsl1_460, hsi0_456, \
                         hsi0_457, hsi1_456, hsi1_457, hsk_584, \
                         hsk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_728[k] = f_10 * hsi0_456[k]
                   - f_11 * hsi1_456[k]
                   + f_3 * pc_x[k] * hsk_584[k];

        t_729[k] = f_10 * hsi0_457[k]
                   - f_11 * hsi1_457[k]
                   + f_3 * pc_x[k] * hsk_585[k];

        t_730[k] = pa_z[k] * gsl0_460[k]
                   - f_14 * pc_z[k] * gsl1_460[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, pc_x, hsi0_459, hsi0_460, hsi0_461, hsi1_459, \
                         hsi1_460, hsi1_461, hsk_587, hsk_588, \
                         hsk_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_8 * hsi0_459[k]
                   - f_9 * hsi1_459[k]
                   + f_3 * pc_x[k] * hsk_587[k];

        t_732[k] = f_8 * hsi0_460[k]
                   - f_9 * hsi1_460[k]
                   + f_3 * pc_x[k] * hsk_588[k];

        t_733[k] = f_8 * hsi0_461[k]
                   - f_9 * hsi1_461[k]
                   + f_3 * pc_x[k] * hsk_589[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pa_z, pc_x, pc_z, gsl0_465, gsl1_465, hsi0_462, \
                         hsi0_464, hsi1_462, hsi1_464, hsk_590, \
                         hsk_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_8 * hsi0_462[k]
                   - f_9 * hsi1_462[k]
                   + f_3 * pc_x[k] * hsk_590[k];

        t_735[k] = pa_z[k] * gsl0_465[k]
                   - f_14 * pc_z[k] * gsl1_465[k];

        t_736[k] = f_6 * hsi0_464[k]
                   - f_7 * hsi1_464[k]
                   + f_3 * pc_x[k] * hsk_592[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pc_x, hsi0_465, hsi0_466, hsi0_467, hsi1_465, \
                         hsi1_466, hsi1_467, hsk_593, hsk_594, \
                         hsk_595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_6 * hsi0_465[k]
                   - f_7 * hsi1_465[k]
                   + f_3 * pc_x[k] * hsk_593[k];

        t_738[k] = f_6 * hsi0_466[k]
                   - f_7 * hsi1_466[k]
                   + f_3 * pc_x[k] * hsk_594[k];

        t_739[k] = f_6 * hsi0_467[k]
                   - f_7 * hsi1_467[k]
                   + f_3 * pc_x[k] * hsk_595[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pa_z, pc_x, pc_z, gsl0_471, gsl1_471, hsi0_468, \
                         hsi0_470, hsi1_468, hsi1_470, hsk_596, \
                         hsk_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_6 * hsi0_468[k]
                   - f_7 * hsi1_468[k]
                   + f_3 * pc_x[k] * hsk_596[k];

        t_741[k] = pa_z[k] * gsl0_471[k]
                   - f_14 * pc_z[k] * gsl1_471[k];

        t_742[k] = f_4 * hsi0_470[k]
                   - f_5 * hsi1_470[k]
                   + f_3 * pc_x[k] * hsk_598[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, pc_x, hsi0_471, hsi0_472, hsi0_473, hsi1_471, \
                         hsi1_472, hsi1_473, hsk_599, hsk_600, \
                         hsk_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_4 * hsi0_471[k]
                   - f_5 * hsi1_471[k]
                   + f_3 * pc_x[k] * hsk_599[k];

        t_744[k] = f_4 * hsi0_472[k]
                   - f_5 * hsi1_472[k]
                   + f_3 * pc_x[k] * hsk_600[k];

        t_745[k] = f_4 * hsi0_473[k]
                   - f_5 * hsi1_473[k]
                   + f_3 * pc_x[k] * hsk_601[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, t_750, pc_x, hsi0_474, hsi0_475, \
                         hsi1_474, hsi1_475, hsk_602, hsk_603, hsk_604, hsk_605, \
                         hsk_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_4 * hsi0_474[k]
                   - f_5 * hsi1_474[k]
                   + f_3 * pc_x[k] * hsk_602[k];

        t_747[k] = f_4 * hsi0_475[k]
                   - f_5 * hsi1_475[k]
                   + f_3 * pc_x[k] * hsk_603[k];

        t_748[k] = f_3 * pc_x[k] * hsk_604[k];

        t_749[k] = f_3 * pc_x[k] * hsk_605[k];

        t_750[k] = f_3 * pc_x[k] * hsk_606[k];
    }

#pragma omp simd aligned(t_751, t_752, t_753, t_754, t_755, t_756, pa_z, pc_x, pc_z, gsl0_486, \
                         gsl1_486, hsk_607, hsk_608, hsk_609, hsk_610, \
                         hsk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_751[k] = f_3 * pc_x[k] * hsk_607[k];

        t_752[k] = f_3 * pc_x[k] * hsk_608[k];

        t_753[k] = f_3 * pc_x[k] * hsk_609[k];

        t_754[k] = f_3 * pc_x[k] * hsk_610[k];

        t_755[k] = f_3 * pc_x[k] * hsk_611[k];

        t_756[k] = pa_z[k] * gsl0_486[k]
                   - f_14 * pc_z[k] * gsl1_486[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, pa_z, pc_z, gsl0_488, gsl0_489, gsk_388, \
                         gsk_389, gsk_390, gsl1_488, gsl1_489, \
                         hsk_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_15 * gsk_388[k]
                   + f_3 * pc_z[k] * hsk_604[k];

        t_758[k] = pa_z[k] * gsl0_488[k]
                   + f_16 * gsk_389[k]
                   - f_14 * pc_z[k] * gsl1_488[k];

        t_759[k] = pa_z[k] * gsl0_489[k]
                   + f_17 * gsk_390[k]
                   - f_14 * pc_z[k] * gsl1_489[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, pa_z, pc_z, gsl0_490, gsl0_491, gsl0_492, \
                         gsk_391, gsk_392, gsk_393, gsl1_490, gsl1_491, \
                         gsl1_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = pa_z[k] * gsl0_490[k]
                   + f_18 * gsk_391[k]
                   - f_14 * pc_z[k] * gsl1_490[k];

        t_761[k] = pa_z[k] * gsl0_491[k]
                   + f_0 * gsk_392[k]
                   - f_14 * pc_z[k] * gsl1_491[k];

        t_762[k] = pa_z[k] * gsl0_492[k]
                   + f_19 * gsk_393[k]
                   - f_14 * pc_z[k] * gsl1_492[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, pc_x, pc_y, pc_z, gsk_395, gsk_431, hsi0_475, \
                         hsi0_476, hsi1_475, hsi1_476, hsk_611, \
                         hsk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_18 * gsk_431[k]
                   + f_3 * pc_y[k] * hsk_611[k];

        t_764[k] = f_15 * gsk_395[k]
                   + f_1 * hsi0_475[k]
                   - f_2 * hsi1_475[k]
                   + f_3 * pc_z[k] * hsk_611[k];

        t_765[k] = f_1 * hsi0_476[k]
                   - f_2 * hsi1_476[k]
                   + f_3 * pc_x[k] * hsk_612[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, pc_x, hsi0_477, hsi0_478, hsi0_479, hsi1_477, \
                         hsi1_478, hsi1_479, hsk_613, hsk_614, \
                         hsk_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_20 * hsi0_477[k]
                   - f_21 * hsi1_477[k]
                   + f_3 * pc_x[k] * hsk_613[k];

        t_767[k] = f_20 * hsi0_478[k]
                   - f_21 * hsi1_478[k]
                   + f_3 * pc_x[k] * hsk_614[k];

        t_768[k] = f_12 * hsi0_479[k]
                   - f_13 * hsi1_479[k]
                   + f_3 * pc_x[k] * hsk_615[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, pc_x, hsi0_480, hsi0_481, hsi0_482, hsi1_480, \
                         hsi1_481, hsi1_482, hsk_616, hsk_617, \
                         hsk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_12 * hsi0_480[k]
                   - f_13 * hsi1_480[k]
                   + f_3 * pc_x[k] * hsk_616[k];

        t_770[k] = f_12 * hsi0_481[k]
                   - f_13 * hsi1_481[k]
                   + f_3 * pc_x[k] * hsk_617[k];

        t_771[k] = f_10 * hsi0_482[k]
                   - f_11 * hsi1_482[k]
                   + f_3 * pc_x[k] * hsk_618[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pc_x, hsi0_483, hsi0_484, hsi0_485, hsi1_483, \
                         hsi1_484, hsi1_485, hsk_619, hsk_620, \
                         hsk_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_10 * hsi0_483[k]
                   - f_11 * hsi1_483[k]
                   + f_3 * pc_x[k] * hsk_619[k];

        t_773[k] = f_10 * hsi0_484[k]
                   - f_11 * hsi1_484[k]
                   + f_3 * pc_x[k] * hsk_620[k];

        t_774[k] = f_10 * hsi0_485[k]
                   - f_11 * hsi1_485[k]
                   + f_3 * pc_x[k] * hsk_621[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, pc_x, hsi0_486, hsi0_487, hsi0_488, hsi1_486, \
                         hsi1_487, hsi1_488, hsk_622, hsk_623, \
                         hsk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_8 * hsi0_486[k]
                   - f_9 * hsi1_486[k]
                   + f_3 * pc_x[k] * hsk_622[k];

        t_776[k] = f_8 * hsi0_487[k]
                   - f_9 * hsi1_487[k]
                   + f_3 * pc_x[k] * hsk_623[k];

        t_777[k] = f_8 * hsi0_488[k]
                   - f_9 * hsi1_488[k]
                   + f_3 * pc_x[k] * hsk_624[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, pc_x, hsi0_489, hsi0_490, hsi0_491, hsi1_489, \
                         hsi1_490, hsi1_491, hsk_625, hsk_626, \
                         hsk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_8 * hsi0_489[k]
                   - f_9 * hsi1_489[k]
                   + f_3 * pc_x[k] * hsk_625[k];

        t_779[k] = f_8 * hsi0_490[k]
                   - f_9 * hsi1_490[k]
                   + f_3 * pc_x[k] * hsk_626[k];

        t_780[k] = f_6 * hsi0_491[k]
                   - f_7 * hsi1_491[k]
                   + f_3 * pc_x[k] * hsk_627[k];
    }

#pragma omp simd aligned(t_781, t_782, t_783, pc_x, hsi0_492, hsi0_493, hsi0_494, hsi1_492, \
                         hsi1_493, hsi1_494, hsk_628, hsk_629, \
                         hsk_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_781[k] = f_6 * hsi0_492[k]
                   - f_7 * hsi1_492[k]
                   + f_3 * pc_x[k] * hsk_628[k];

        t_782[k] = f_6 * hsi0_493[k]
                   - f_7 * hsi1_493[k]
                   + f_3 * pc_x[k] * hsk_629[k];

        t_783[k] = f_6 * hsi0_494[k]
                   - f_7 * hsi1_494[k]
                   + f_3 * pc_x[k] * hsk_630[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, pc_x, hsi0_495, hsi0_496, hsi0_497, hsi1_495, \
                         hsi1_496, hsi1_497, hsk_631, hsk_632, \
                         hsk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_6 * hsi0_495[k]
                   - f_7 * hsi1_495[k]
                   + f_3 * pc_x[k] * hsk_631[k];

        t_785[k] = f_6 * hsi0_496[k]
                   - f_7 * hsi1_496[k]
                   + f_3 * pc_x[k] * hsk_632[k];

        t_786[k] = f_4 * hsi0_497[k]
                   - f_5 * hsi1_497[k]
                   + f_3 * pc_x[k] * hsk_633[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_x, hsi0_498, hsi0_499, hsi0_500, hsi1_498, \
                         hsi1_499, hsi1_500, hsk_634, hsk_635, \
                         hsk_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_4 * hsi0_498[k]
                   - f_5 * hsi1_498[k]
                   + f_3 * pc_x[k] * hsk_634[k];

        t_788[k] = f_4 * hsi0_499[k]
                   - f_5 * hsi1_499[k]
                   + f_3 * pc_x[k] * hsk_635[k];

        t_789[k] = f_4 * hsi0_500[k]
                   - f_5 * hsi1_500[k]
                   + f_3 * pc_x[k] * hsk_636[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pc_x, hsi0_501, hsi0_502, hsi0_503, \
                         hsi1_501, hsi1_502, hsi1_503, hsk_637, hsk_638, hsk_639, \
                         hsk_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_4 * hsi0_501[k]
                   - f_5 * hsi1_501[k]
                   + f_3 * pc_x[k] * hsk_637[k];

        t_791[k] = f_4 * hsi0_502[k]
                   - f_5 * hsi1_502[k]
                   + f_3 * pc_x[k] * hsk_638[k];

        t_792[k] = f_4 * hsi0_503[k]
                   - f_5 * hsi1_503[k]
                   + f_3 * pc_x[k] * hsk_639[k];

        t_793[k] = f_3 * pc_x[k] * hsk_640[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, t_799, t_800, pc_x, hsk_641, \
                         hsk_642, hsk_643, hsk_644, hsk_645, hsk_646, \
                         hsk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_3 * pc_x[k] * hsk_641[k];

        t_795[k] = f_3 * pc_x[k] * hsk_642[k];

        t_796[k] = f_3 * pc_x[k] * hsk_643[k];

        t_797[k] = f_3 * pc_x[k] * hsk_644[k];

        t_798[k] = f_3 * pc_x[k] * hsk_645[k];

        t_799[k] = f_3 * pc_x[k] * hsk_646[k];

        t_800[k] = f_3 * pc_x[k] * hsk_647[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, pc_y, pc_z, gsk_424, gsk_460, gsk_462, hsi0_497, \
                         hsi0_499, hsi1_497, hsi1_499, hsk_640, \
                         hsk_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_17 * gsk_460[k]
                   + f_1 * hsi0_497[k]
                   - f_2 * hsi1_497[k]
                   + f_3 * pc_y[k] * hsk_640[k];

        t_802[k] = f_16 * gsk_424[k]
                   + f_3 * pc_z[k] * hsk_640[k];

        t_803[k] = f_17 * gsk_462[k]
                   + f_12 * hsi0_499[k]
                   - f_13 * hsi1_499[k]
                   + f_3 * pc_y[k] * hsk_642[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pc_y, gsk_463, gsk_464, gsk_465, hsi0_500, \
                         hsi0_501, hsi0_502, hsi1_500, hsi1_501, hsi1_502, hsk_643, hsk_644, \
                         hsk_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_17 * gsk_463[k]
                   + f_10 * hsi0_500[k]
                   - f_11 * hsi1_500[k]
                   + f_3 * pc_y[k] * hsk_643[k];

        t_805[k] = f_17 * gsk_464[k]
                   + f_8 * hsi0_501[k]
                   - f_9 * hsi1_501[k]
                   + f_3 * pc_y[k] * hsk_644[k];

        t_806[k] = f_17 * gsk_465[k]
                   + f_6 * hsi0_502[k]
                   - f_7 * hsi1_502[k]
                   + f_3 * pc_y[k] * hsk_645[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, pc_y, pc_z, gsk_431, gsk_466, gsk_467, hsi0_503, \
                         hsi1_503, hsk_646, hsk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = f_17 * gsk_466[k]
                   + f_4 * hsi0_503[k]
                   - f_5 * hsi1_503[k]
                   + f_3 * pc_y[k] * hsk_646[k];

        t_808[k] = f_17 * gsk_467[k]
                   + f_3 * pc_y[k] * hsk_647[k];

        t_809[k] = f_16 * gsk_431[k]
                   + f_1 * hsi0_503[k]
                   - f_2 * hsi1_503[k]
                   + f_3 * pc_z[k] * hsk_647[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, pc_x, hsi0_504, hsi0_505, hsi0_506, hsi1_504, \
                         hsi1_505, hsi1_506, hsk_648, hsk_649, \
                         hsk_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_1 * hsi0_504[k]
                   - f_2 * hsi1_504[k]
                   + f_3 * pc_x[k] * hsk_648[k];

        t_811[k] = f_20 * hsi0_505[k]
                   - f_21 * hsi1_505[k]
                   + f_3 * pc_x[k] * hsk_649[k];

        t_812[k] = f_20 * hsi0_506[k]
                   - f_21 * hsi1_506[k]
                   + f_3 * pc_x[k] * hsk_650[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, pc_x, hsi0_507, hsi0_508, hsi0_509, hsi1_507, \
                         hsi1_508, hsi1_509, hsk_651, hsk_652, \
                         hsk_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_12 * hsi0_507[k]
                   - f_13 * hsi1_507[k]
                   + f_3 * pc_x[k] * hsk_651[k];

        t_814[k] = f_12 * hsi0_508[k]
                   - f_13 * hsi1_508[k]
                   + f_3 * pc_x[k] * hsk_652[k];

        t_815[k] = f_12 * hsi0_509[k]
                   - f_13 * hsi1_509[k]
                   + f_3 * pc_x[k] * hsk_653[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pc_x, hsi0_510, hsi0_511, hsi0_512, hsi1_510, \
                         hsi1_511, hsi1_512, hsk_654, hsk_655, \
                         hsk_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_10 * hsi0_510[k]
                   - f_11 * hsi1_510[k]
                   + f_3 * pc_x[k] * hsk_654[k];

        t_817[k] = f_10 * hsi0_511[k]
                   - f_11 * hsi1_511[k]
                   + f_3 * pc_x[k] * hsk_655[k];

        t_818[k] = f_10 * hsi0_512[k]
                   - f_11 * hsi1_512[k]
                   + f_3 * pc_x[k] * hsk_656[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, hsi0_513, hsi0_514, hsi0_515, hsi1_513, \
                         hsi1_514, hsi1_515, hsk_657, hsk_658, \
                         hsk_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_10 * hsi0_513[k]
                   - f_11 * hsi1_513[k]
                   + f_3 * pc_x[k] * hsk_657[k];

        t_820[k] = f_8 * hsi0_514[k]
                   - f_9 * hsi1_514[k]
                   + f_3 * pc_x[k] * hsk_658[k];

        t_821[k] = f_8 * hsi0_515[k]
                   - f_9 * hsi1_515[k]
                   + f_3 * pc_x[k] * hsk_659[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pc_x, hsi0_516, hsi0_517, hsi0_518, hsi1_516, \
                         hsi1_517, hsi1_518, hsk_660, hsk_661, \
                         hsk_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_8 * hsi0_516[k]
                   - f_9 * hsi1_516[k]
                   + f_3 * pc_x[k] * hsk_660[k];

        t_823[k] = f_8 * hsi0_517[k]
                   - f_9 * hsi1_517[k]
                   + f_3 * pc_x[k] * hsk_661[k];

        t_824[k] = f_8 * hsi0_518[k]
                   - f_9 * hsi1_518[k]
                   + f_3 * pc_x[k] * hsk_662[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, pc_x, hsi0_519, hsi0_520, hsi0_521, hsi1_519, \
                         hsi1_520, hsi1_521, hsk_663, hsk_664, \
                         hsk_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = f_6 * hsi0_519[k]
                   - f_7 * hsi1_519[k]
                   + f_3 * pc_x[k] * hsk_663[k];

        t_826[k] = f_6 * hsi0_520[k]
                   - f_7 * hsi1_520[k]
                   + f_3 * pc_x[k] * hsk_664[k];

        t_827[k] = f_6 * hsi0_521[k]
                   - f_7 * hsi1_521[k]
                   + f_3 * pc_x[k] * hsk_665[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, pc_x, hsi0_522, hsi0_523, hsi0_524, hsi1_522, \
                         hsi1_523, hsi1_524, hsk_666, hsk_667, \
                         hsk_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = f_6 * hsi0_522[k]
                   - f_7 * hsi1_522[k]
                   + f_3 * pc_x[k] * hsk_666[k];

        t_829[k] = f_6 * hsi0_523[k]
                   - f_7 * hsi1_523[k]
                   + f_3 * pc_x[k] * hsk_667[k];

        t_830[k] = f_6 * hsi0_524[k]
                   - f_7 * hsi1_524[k]
                   + f_3 * pc_x[k] * hsk_668[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, pc_x, hsi0_525, hsi0_526, hsi0_527, hsi1_525, \
                         hsi1_526, hsi1_527, hsk_669, hsk_670, \
                         hsk_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_4 * hsi0_525[k]
                   - f_5 * hsi1_525[k]
                   + f_3 * pc_x[k] * hsk_669[k];

        t_832[k] = f_4 * hsi0_526[k]
                   - f_5 * hsi1_526[k]
                   + f_3 * pc_x[k] * hsk_670[k];

        t_833[k] = f_4 * hsi0_527[k]
                   - f_5 * hsi1_527[k]
                   + f_3 * pc_x[k] * hsk_671[k];
    }
}

static auto
compute_prim_hsl_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsl0,
                                                          const size_t gsk, const size_t gsl1,
                                                          const size_t hsi0, const size_t hsi1,
                                                          const size_t hsk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsl0_630 = buffer.data(gsl0 + 630);
    const auto *gsl0_632 = buffer.data(gsl0 + 632);
    const auto *gsl0_635 = buffer.data(gsl0 + 635);
    const auto *gsl0_639 = buffer.data(gsl0 + 639);
    const auto *gsl0_644 = buffer.data(gsl0 + 644);
    const auto *gsl0_650 = buffer.data(gsl0 + 650);
    const auto *gsl0_657 = buffer.data(gsl0 + 657);
    const auto *gsl0_666 = buffer.data(gsl0 + 666);
    const auto *gsl0_668 = buffer.data(gsl0 + 668);
    const auto *gsl0_669 = buffer.data(gsl0 + 669);
    const auto *gsl0_670 = buffer.data(gsl0 + 670);
    const auto *gsl0_671 = buffer.data(gsl0 + 671);
    const auto *gsl0_672 = buffer.data(gsl0 + 672);
    const auto *gsl0_674 = buffer.data(gsl0 + 674);

    const auto *gsk_460 = buffer.data(gsk + 460);
    const auto *gsk_467 = buffer.data(gsk + 467);
    const auto *gsk_496 = buffer.data(gsk + 496);
    const auto *gsk_498 = buffer.data(gsk + 498);
    const auto *gsk_499 = buffer.data(gsk + 499);
    const auto *gsk_500 = buffer.data(gsk + 500);
    const auto *gsk_501 = buffer.data(gsk + 501);
    const auto *gsk_502 = buffer.data(gsk + 502);
    const auto *gsk_503 = buffer.data(gsk + 503);
    const auto *gsk_532 = buffer.data(gsk + 532);
    const auto *gsk_534 = buffer.data(gsk + 534);
    const auto *gsk_535 = buffer.data(gsk + 535);
    const auto *gsk_536 = buffer.data(gsk + 536);
    const auto *gsk_537 = buffer.data(gsk + 537);
    const auto *gsk_538 = buffer.data(gsk + 538);
    const auto *gsk_539 = buffer.data(gsk + 539);

    const auto *gsl1_630 = buffer.data(gsl1 + 630);
    const auto *gsl1_632 = buffer.data(gsl1 + 632);
    const auto *gsl1_635 = buffer.data(gsl1 + 635);
    const auto *gsl1_639 = buffer.data(gsl1 + 639);
    const auto *gsl1_644 = buffer.data(gsl1 + 644);
    const auto *gsl1_650 = buffer.data(gsl1 + 650);
    const auto *gsl1_657 = buffer.data(gsl1 + 657);
    const auto *gsl1_666 = buffer.data(gsl1 + 666);
    const auto *gsl1_668 = buffer.data(gsl1 + 668);
    const auto *gsl1_669 = buffer.data(gsl1 + 669);
    const auto *gsl1_670 = buffer.data(gsl1 + 670);
    const auto *gsl1_671 = buffer.data(gsl1 + 671);
    const auto *gsl1_672 = buffer.data(gsl1 + 672);
    const auto *gsl1_674 = buffer.data(gsl1 + 674);

    const auto *hsi0_525 = buffer.data(hsi0 + 525);
    const auto *hsi0_527 = buffer.data(hsi0 + 527);
    const auto *hsi0_528 = buffer.data(hsi0 + 528);
    const auto *hsi0_529 = buffer.data(hsi0 + 529);
    const auto *hsi0_530 = buffer.data(hsi0 + 530);
    const auto *hsi0_531 = buffer.data(hsi0 + 531);
    const auto *hsi0_533 = buffer.data(hsi0 + 533);
    const auto *hsi0_535 = buffer.data(hsi0 + 535);
    const auto *hsi0_536 = buffer.data(hsi0 + 536);
    const auto *hsi0_538 = buffer.data(hsi0 + 538);
    const auto *hsi0_539 = buffer.data(hsi0 + 539);
    const auto *hsi0_540 = buffer.data(hsi0 + 540);
    const auto *hsi0_542 = buffer.data(hsi0 + 542);
    const auto *hsi0_543 = buffer.data(hsi0 + 543);
    const auto *hsi0_544 = buffer.data(hsi0 + 544);
    const auto *hsi0_545 = buffer.data(hsi0 + 545);
    const auto *hsi0_547 = buffer.data(hsi0 + 547);
    const auto *hsi0_548 = buffer.data(hsi0 + 548);
    const auto *hsi0_549 = buffer.data(hsi0 + 549);
    const auto *hsi0_550 = buffer.data(hsi0 + 550);
    const auto *hsi0_551 = buffer.data(hsi0 + 551);
    const auto *hsi0_553 = buffer.data(hsi0 + 553);
    const auto *hsi0_554 = buffer.data(hsi0 + 554);
    const auto *hsi0_555 = buffer.data(hsi0 + 555);
    const auto *hsi0_556 = buffer.data(hsi0 + 556);
    const auto *hsi0_557 = buffer.data(hsi0 + 557);
    const auto *hsi0_558 = buffer.data(hsi0 + 558);
    const auto *hsi0_560 = buffer.data(hsi0 + 560);
    const auto *hsi0_562 = buffer.data(hsi0 + 562);
    const auto *hsi0_563 = buffer.data(hsi0 + 563);
    const auto *hsi0_565 = buffer.data(hsi0 + 565);
    const auto *hsi0_566 = buffer.data(hsi0 + 566);
    const auto *hsi0_567 = buffer.data(hsi0 + 567);
    const auto *hsi0_569 = buffer.data(hsi0 + 569);
    const auto *hsi0_570 = buffer.data(hsi0 + 570);
    const auto *hsi0_571 = buffer.data(hsi0 + 571);
    const auto *hsi0_572 = buffer.data(hsi0 + 572);
    const auto *hsi0_574 = buffer.data(hsi0 + 574);
    const auto *hsi0_575 = buffer.data(hsi0 + 575);
    const auto *hsi0_576 = buffer.data(hsi0 + 576);
    const auto *hsi0_577 = buffer.data(hsi0 + 577);
    const auto *hsi0_578 = buffer.data(hsi0 + 578);
    const auto *hsi0_580 = buffer.data(hsi0 + 580);
    const auto *hsi0_581 = buffer.data(hsi0 + 581);
    const auto *hsi0_582 = buffer.data(hsi0 + 582);
    const auto *hsi0_583 = buffer.data(hsi0 + 583);
    const auto *hsi0_584 = buffer.data(hsi0 + 584);
    const auto *hsi0_585 = buffer.data(hsi0 + 585);
    const auto *hsi0_586 = buffer.data(hsi0 + 586);
    const auto *hsi0_587 = buffer.data(hsi0 + 587);

    const auto *hsi1_525 = buffer.data(hsi1 + 525);
    const auto *hsi1_527 = buffer.data(hsi1 + 527);
    const auto *hsi1_528 = buffer.data(hsi1 + 528);
    const auto *hsi1_529 = buffer.data(hsi1 + 529);
    const auto *hsi1_530 = buffer.data(hsi1 + 530);
    const auto *hsi1_531 = buffer.data(hsi1 + 531);
    const auto *hsi1_533 = buffer.data(hsi1 + 533);
    const auto *hsi1_535 = buffer.data(hsi1 + 535);
    const auto *hsi1_536 = buffer.data(hsi1 + 536);
    const auto *hsi1_538 = buffer.data(hsi1 + 538);
    const auto *hsi1_539 = buffer.data(hsi1 + 539);
    const auto *hsi1_540 = buffer.data(hsi1 + 540);
    const auto *hsi1_542 = buffer.data(hsi1 + 542);
    const auto *hsi1_543 = buffer.data(hsi1 + 543);
    const auto *hsi1_544 = buffer.data(hsi1 + 544);
    const auto *hsi1_545 = buffer.data(hsi1 + 545);
    const auto *hsi1_547 = buffer.data(hsi1 + 547);
    const auto *hsi1_548 = buffer.data(hsi1 + 548);
    const auto *hsi1_549 = buffer.data(hsi1 + 549);
    const auto *hsi1_550 = buffer.data(hsi1 + 550);
    const auto *hsi1_551 = buffer.data(hsi1 + 551);
    const auto *hsi1_553 = buffer.data(hsi1 + 553);
    const auto *hsi1_554 = buffer.data(hsi1 + 554);
    const auto *hsi1_555 = buffer.data(hsi1 + 555);
    const auto *hsi1_556 = buffer.data(hsi1 + 556);
    const auto *hsi1_557 = buffer.data(hsi1 + 557);
    const auto *hsi1_558 = buffer.data(hsi1 + 558);
    const auto *hsi1_560 = buffer.data(hsi1 + 560);
    const auto *hsi1_562 = buffer.data(hsi1 + 562);
    const auto *hsi1_563 = buffer.data(hsi1 + 563);
    const auto *hsi1_565 = buffer.data(hsi1 + 565);
    const auto *hsi1_566 = buffer.data(hsi1 + 566);
    const auto *hsi1_567 = buffer.data(hsi1 + 567);
    const auto *hsi1_569 = buffer.data(hsi1 + 569);
    const auto *hsi1_570 = buffer.data(hsi1 + 570);
    const auto *hsi1_571 = buffer.data(hsi1 + 571);
    const auto *hsi1_572 = buffer.data(hsi1 + 572);
    const auto *hsi1_574 = buffer.data(hsi1 + 574);
    const auto *hsi1_575 = buffer.data(hsi1 + 575);
    const auto *hsi1_576 = buffer.data(hsi1 + 576);
    const auto *hsi1_577 = buffer.data(hsi1 + 577);
    const auto *hsi1_578 = buffer.data(hsi1 + 578);
    const auto *hsi1_580 = buffer.data(hsi1 + 580);
    const auto *hsi1_581 = buffer.data(hsi1 + 581);
    const auto *hsi1_582 = buffer.data(hsi1 + 582);
    const auto *hsi1_583 = buffer.data(hsi1 + 583);
    const auto *hsi1_584 = buffer.data(hsi1 + 584);
    const auto *hsi1_585 = buffer.data(hsi1 + 585);
    const auto *hsi1_586 = buffer.data(hsi1 + 586);
    const auto *hsi1_587 = buffer.data(hsi1 + 587);

    const auto *hsk_672 = buffer.data(hsk + 672);
    const auto *hsk_673 = buffer.data(hsk + 673);
    const auto *hsk_674 = buffer.data(hsk + 674);
    const auto *hsk_675 = buffer.data(hsk + 675);
    const auto *hsk_676 = buffer.data(hsk + 676);
    const auto *hsk_677 = buffer.data(hsk + 677);
    const auto *hsk_678 = buffer.data(hsk + 678);
    const auto *hsk_679 = buffer.data(hsk + 679);
    const auto *hsk_680 = buffer.data(hsk + 680);
    const auto *hsk_681 = buffer.data(hsk + 681);
    const auto *hsk_682 = buffer.data(hsk + 682);
    const auto *hsk_683 = buffer.data(hsk + 683);
    const auto *hsk_685 = buffer.data(hsk + 685);
    const auto *hsk_687 = buffer.data(hsk + 687);
    const auto *hsk_688 = buffer.data(hsk + 688);
    const auto *hsk_690 = buffer.data(hsk + 690);
    const auto *hsk_691 = buffer.data(hsk + 691);
    const auto *hsk_692 = buffer.data(hsk + 692);
    const auto *hsk_694 = buffer.data(hsk + 694);
    const auto *hsk_695 = buffer.data(hsk + 695);
    const auto *hsk_696 = buffer.data(hsk + 696);
    const auto *hsk_697 = buffer.data(hsk + 697);
    const auto *hsk_699 = buffer.data(hsk + 699);
    const auto *hsk_700 = buffer.data(hsk + 700);
    const auto *hsk_701 = buffer.data(hsk + 701);
    const auto *hsk_702 = buffer.data(hsk + 702);
    const auto *hsk_703 = buffer.data(hsk + 703);
    const auto *hsk_705 = buffer.data(hsk + 705);
    const auto *hsk_706 = buffer.data(hsk + 706);
    const auto *hsk_707 = buffer.data(hsk + 707);
    const auto *hsk_708 = buffer.data(hsk + 708);
    const auto *hsk_709 = buffer.data(hsk + 709);
    const auto *hsk_710 = buffer.data(hsk + 710);
    const auto *hsk_712 = buffer.data(hsk + 712);
    const auto *hsk_713 = buffer.data(hsk + 713);
    const auto *hsk_714 = buffer.data(hsk + 714);
    const auto *hsk_715 = buffer.data(hsk + 715);
    const auto *hsk_716 = buffer.data(hsk + 716);
    const auto *hsk_717 = buffer.data(hsk + 717);
    const auto *hsk_718 = buffer.data(hsk + 718);
    const auto *hsk_719 = buffer.data(hsk + 719);
    const auto *hsk_720 = buffer.data(hsk + 720);
    const auto *hsk_722 = buffer.data(hsk + 722);
    const auto *hsk_723 = buffer.data(hsk + 723);
    const auto *hsk_725 = buffer.data(hsk + 725);
    const auto *hsk_726 = buffer.data(hsk + 726);
    const auto *hsk_727 = buffer.data(hsk + 727);
    const auto *hsk_729 = buffer.data(hsk + 729);
    const auto *hsk_730 = buffer.data(hsk + 730);
    const auto *hsk_731 = buffer.data(hsk + 731);
    const auto *hsk_732 = buffer.data(hsk + 732);
    const auto *hsk_734 = buffer.data(hsk + 734);
    const auto *hsk_735 = buffer.data(hsk + 735);
    const auto *hsk_736 = buffer.data(hsk + 736);
    const auto *hsk_737 = buffer.data(hsk + 737);
    const auto *hsk_738 = buffer.data(hsk + 738);
    const auto *hsk_740 = buffer.data(hsk + 740);
    const auto *hsk_741 = buffer.data(hsk + 741);
    const auto *hsk_742 = buffer.data(hsk + 742);
    const auto *hsk_743 = buffer.data(hsk + 743);
    const auto *hsk_744 = buffer.data(hsk + 744);
    const auto *hsk_745 = buffer.data(hsk + 745);
    const auto *hsk_747 = buffer.data(hsk + 747);
    const auto *hsk_748 = buffer.data(hsk + 748);
    const auto *hsk_749 = buffer.data(hsk + 749);
    const auto *hsk_750 = buffer.data(hsk + 750);
    const auto *hsk_751 = buffer.data(hsk + 751);
    const auto *hsk_752 = buffer.data(hsk + 752);
    const auto *hsk_753 = buffer.data(hsk + 753);
    const auto *hsk_754 = buffer.data(hsk + 754);
    const auto *hsk_755 = buffer.data(hsk + 755);

#pragma omp simd aligned(t_834, t_835, t_836, pc_x, hsi0_528, hsi0_529, hsi0_530, hsi1_528, \
                         hsi1_529, hsi1_530, hsk_672, hsk_673, \
                         hsk_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_4 * hsi0_528[k]
                   - f_5 * hsi1_528[k]
                   + f_3 * pc_x[k] * hsk_672[k];

        t_835[k] = f_4 * hsi0_529[k]
                   - f_5 * hsi1_529[k]
                   + f_3 * pc_x[k] * hsk_673[k];

        t_836[k] = f_4 * hsi0_530[k]
                   - f_5 * hsi1_530[k]
                   + f_3 * pc_x[k] * hsk_674[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, t_841, t_842, pc_x, hsi0_531, hsi1_531, \
                         hsk_675, hsk_676, hsk_677, hsk_678, hsk_679, \
                         hsk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_4 * hsi0_531[k]
                   - f_5 * hsi1_531[k]
                   + f_3 * pc_x[k] * hsk_675[k];

        t_838[k] = f_3 * pc_x[k] * hsk_676[k];

        t_839[k] = f_3 * pc_x[k] * hsk_677[k];

        t_840[k] = f_3 * pc_x[k] * hsk_678[k];

        t_841[k] = f_3 * pc_x[k] * hsk_679[k];

        t_842[k] = f_3 * pc_x[k] * hsk_680[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, t_846, t_847, pc_x, pc_y, pc_z, gsk_460, \
                         gsk_496, hsi0_525, hsi1_525, hsk_676, hsk_681, hsk_682, \
                         hsk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_3 * pc_x[k] * hsk_681[k];

        t_844[k] = f_3 * pc_x[k] * hsk_682[k];

        t_845[k] = f_3 * pc_x[k] * hsk_683[k];

        t_846[k] = f_16 * gsk_496[k]
                   + f_1 * hsi0_525[k]
                   - f_2 * hsi1_525[k]
                   + f_3 * pc_y[k] * hsk_676[k];

        t_847[k] = f_17 * gsk_460[k]
                   + f_3 * pc_z[k] * hsk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pc_y, gsk_498, gsk_499, gsk_500, hsi0_527, \
                         hsi0_528, hsi0_529, hsi1_527, hsi1_528, hsi1_529, hsk_678, hsk_679, \
                         hsk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * gsk_498[k]
                   + f_12 * hsi0_527[k]
                   - f_13 * hsi1_527[k]
                   + f_3 * pc_y[k] * hsk_678[k];

        t_849[k] = f_16 * gsk_499[k]
                   + f_10 * hsi0_528[k]
                   - f_11 * hsi1_528[k]
                   + f_3 * pc_y[k] * hsk_679[k];

        t_850[k] = f_16 * gsk_500[k]
                   + f_8 * hsi0_529[k]
                   - f_9 * hsi1_529[k]
                   + f_3 * pc_y[k] * hsk_680[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pc_y, gsk_501, gsk_502, gsk_503, hsi0_530, \
                         hsi0_531, hsi1_530, hsi1_531, hsk_681, hsk_682, \
                         hsk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_16 * gsk_501[k]
                   + f_6 * hsi0_530[k]
                   - f_7 * hsi1_530[k]
                   + f_3 * pc_y[k] * hsk_681[k];

        t_852[k] = f_16 * gsk_502[k]
                   + f_4 * hsi0_531[k]
                   - f_5 * hsi1_531[k]
                   + f_3 * pc_y[k] * hsk_682[k];

        t_853[k] = f_16 * gsk_503[k]
                   + f_3 * pc_y[k] * hsk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, pa_y, pc_x, pc_y, pc_z, gsl0_630, gsk_467, \
                         gsl1_630, hsi0_531, hsi0_533, hsi1_531, hsi1_533, hsk_683, \
                         hsk_685 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_17 * gsk_467[k]
                   + f_1 * hsi0_531[k]
                   - f_2 * hsi1_531[k]
                   + f_3 * pc_z[k] * hsk_683[k];

        t_855[k] = pa_y[k] * gsl0_630[k]
                   - f_14 * pc_y[k] * gsl1_630[k];

        t_856[k] = f_20 * hsi0_533[k]
                   - f_21 * hsi1_533[k]
                   + f_3 * pc_x[k] * hsk_685[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, pa_y, pc_x, pc_y, gsl0_632, gsl1_632, hsi0_535, \
                         hsi0_536, hsi1_535, hsi1_536, hsk_687, \
                         hsk_688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = pa_y[k] * gsl0_632[k]
                   - f_14 * pc_y[k] * gsl1_632[k];

        t_858[k] = f_12 * hsi0_535[k]
                   - f_13 * hsi1_535[k]
                   + f_3 * pc_x[k] * hsk_687[k];

        t_859[k] = f_12 * hsi0_536[k]
                   - f_13 * hsi1_536[k]
                   + f_3 * pc_x[k] * hsk_688[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pa_y, pc_x, pc_y, gsl0_635, gsl1_635, hsi0_538, \
                         hsi0_539, hsi1_538, hsi1_539, hsk_690, \
                         hsk_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = pa_y[k] * gsl0_635[k]
                   - f_14 * pc_y[k] * gsl1_635[k];

        t_861[k] = f_10 * hsi0_538[k]
                   - f_11 * hsi1_538[k]
                   + f_3 * pc_x[k] * hsk_690[k];

        t_862[k] = f_10 * hsi0_539[k]
                   - f_11 * hsi1_539[k]
                   + f_3 * pc_x[k] * hsk_691[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pa_y, pc_x, pc_y, gsl0_639, gsl1_639, hsi0_540, \
                         hsi0_542, hsi1_540, hsi1_542, hsk_692, \
                         hsk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_10 * hsi0_540[k]
                   - f_11 * hsi1_540[k]
                   + f_3 * pc_x[k] * hsk_692[k];

        t_864[k] = pa_y[k] * gsl0_639[k]
                   - f_14 * pc_y[k] * gsl1_639[k];

        t_865[k] = f_8 * hsi0_542[k]
                   - f_9 * hsi1_542[k]
                   + f_3 * pc_x[k] * hsk_694[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, hsi0_543, hsi0_544, hsi0_545, hsi1_543, \
                         hsi1_544, hsi1_545, hsk_695, hsk_696, \
                         hsk_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_8 * hsi0_543[k]
                   - f_9 * hsi1_543[k]
                   + f_3 * pc_x[k] * hsk_695[k];

        t_867[k] = f_8 * hsi0_544[k]
                   - f_9 * hsi1_544[k]
                   + f_3 * pc_x[k] * hsk_696[k];

        t_868[k] = f_8 * hsi0_545[k]
                   - f_9 * hsi1_545[k]
                   + f_3 * pc_x[k] * hsk_697[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pa_y, pc_x, pc_y, gsl0_644, gsl1_644, hsi0_547, \
                         hsi0_548, hsi1_547, hsi1_548, hsk_699, \
                         hsk_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = pa_y[k] * gsl0_644[k]
                   - f_14 * pc_y[k] * gsl1_644[k];

        t_870[k] = f_6 * hsi0_547[k]
                   - f_7 * hsi1_547[k]
                   + f_3 * pc_x[k] * hsk_699[k];

        t_871[k] = f_6 * hsi0_548[k]
                   - f_7 * hsi1_548[k]
                   + f_3 * pc_x[k] * hsk_700[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pc_x, hsi0_549, hsi0_550, hsi0_551, hsi1_549, \
                         hsi1_550, hsi1_551, hsk_701, hsk_702, \
                         hsk_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_6 * hsi0_549[k]
                   - f_7 * hsi1_549[k]
                   + f_3 * pc_x[k] * hsk_701[k];

        t_873[k] = f_6 * hsi0_550[k]
                   - f_7 * hsi1_550[k]
                   + f_3 * pc_x[k] * hsk_702[k];

        t_874[k] = f_6 * hsi0_551[k]
                   - f_7 * hsi1_551[k]
                   + f_3 * pc_x[k] * hsk_703[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pa_y, pc_x, pc_y, gsl0_650, gsl1_650, hsi0_553, \
                         hsi0_554, hsi1_553, hsi1_554, hsk_705, \
                         hsk_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = pa_y[k] * gsl0_650[k]
                   - f_14 * pc_y[k] * gsl1_650[k];

        t_876[k] = f_4 * hsi0_553[k]
                   - f_5 * hsi1_553[k]
                   + f_3 * pc_x[k] * hsk_705[k];

        t_877[k] = f_4 * hsi0_554[k]
                   - f_5 * hsi1_554[k]
                   + f_3 * pc_x[k] * hsk_706[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_x, hsi0_555, hsi0_556, hsi0_557, hsi1_555, \
                         hsi1_556, hsi1_557, hsk_707, hsk_708, \
                         hsk_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_4 * hsi0_555[k]
                   - f_5 * hsi1_555[k]
                   + f_3 * pc_x[k] * hsk_707[k];

        t_879[k] = f_4 * hsi0_556[k]
                   - f_5 * hsi1_556[k]
                   + f_3 * pc_x[k] * hsk_708[k];

        t_880[k] = f_4 * hsi0_557[k]
                   - f_5 * hsi1_557[k]
                   + f_3 * pc_x[k] * hsk_709[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, t_884, t_885, pa_y, pc_x, pc_y, gsl0_657, \
                         gsl1_657, hsi0_558, hsi1_558, hsk_710, hsk_712, hsk_713, \
                         hsk_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_4 * hsi0_558[k]
                   - f_5 * hsi1_558[k]
                   + f_3 * pc_x[k] * hsk_710[k];

        t_882[k] = pa_y[k] * gsl0_657[k]
                   - f_14 * pc_y[k] * gsl1_657[k];

        t_883[k] = f_3 * pc_x[k] * hsk_712[k];

        t_884[k] = f_3 * pc_x[k] * hsk_713[k];

        t_885[k] = f_3 * pc_x[k] * hsk_714[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, t_890, pc_x, hsk_715, hsk_716, hsk_717, \
                         hsk_718, hsk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_3 * pc_x[k] * hsk_715[k];

        t_887[k] = f_3 * pc_x[k] * hsk_716[k];

        t_888[k] = f_3 * pc_x[k] * hsk_717[k];

        t_889[k] = f_3 * pc_x[k] * hsk_718[k];

        t_890[k] = f_3 * pc_x[k] * hsk_719[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, pa_y, pc_y, pc_z, gsl0_666, gsl0_668, gsk_496, \
                         gsk_532, gsk_534, gsl1_666, gsl1_668, \
                         hsk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = pa_y[k] * gsl0_666[k]
                   + f_22 * gsk_532[k]
                   - f_14 * pc_y[k] * gsl1_666[k];

        t_892[k] = f_18 * gsk_496[k]
                   + f_3 * pc_z[k] * hsk_712[k];

        t_893[k] = pa_y[k] * gsl0_668[k]
                   + f_19 * gsk_534[k]
                   - f_14 * pc_y[k] * gsl1_668[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, pa_y, pc_y, gsl0_669, gsl0_670, gsl0_671, \
                         gsk_535, gsk_536, gsk_537, gsl1_669, gsl1_670, \
                         gsl1_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = pa_y[k] * gsl0_669[k]
                   + f_0 * gsk_535[k]
                   - f_14 * pc_y[k] * gsl1_669[k];

        t_895[k] = pa_y[k] * gsl0_670[k]
                   + f_18 * gsk_536[k]
                   - f_14 * pc_y[k] * gsl1_670[k];

        t_896[k] = pa_y[k] * gsl0_671[k]
                   + f_17 * gsk_537[k]
                   - f_14 * pc_y[k] * gsl1_671[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pa_y, pc_y, gsl0_672, gsl0_674, gsk_538, \
                         gsk_539, gsl1_672, gsl1_674, hsk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = pa_y[k] * gsl0_672[k]
                   + f_16 * gsk_538[k]
                   - f_14 * pc_y[k] * gsl1_672[k];

        t_898[k] = f_15 * gsk_539[k]
                   + f_3 * pc_y[k] * hsk_719[k];

        t_899[k] = pa_y[k] * gsl0_674[k]
                   - f_14 * pc_y[k] * gsl1_674[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, pc_x, pc_y, hsi0_560, hsi0_562, \
                         hsi0_563, hsi1_560, hsi1_562, hsi1_563, hsk_720, hsk_722, \
                         hsk_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_1 * hsi0_560[k]
                   - f_2 * hsi1_560[k]
                   + f_3 * pc_x[k] * hsk_720[k];

        t_901[k] = f_3 * pc_y[k] * hsk_720[k];

        t_902[k] = f_20 * hsi0_562[k]
                   - f_21 * hsi1_562[k]
                   + f_3 * pc_x[k] * hsk_722[k];

        t_903[k] = f_12 * hsi0_563[k]
                   - f_13 * hsi1_563[k]
                   + f_3 * pc_x[k] * hsk_723[k];

        t_904[k] = f_3 * pc_y[k] * hsk_722[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pc_x, pc_y, hsi0_565, hsi0_566, hsi0_567, \
                         hsi1_565, hsi1_566, hsi1_567, hsk_725, hsk_726, \
                         hsk_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_12 * hsi0_565[k]
                   - f_13 * hsi1_565[k]
                   + f_3 * pc_x[k] * hsk_725[k];

        t_906[k] = f_10 * hsi0_566[k]
                   - f_11 * hsi1_566[k]
                   + f_3 * pc_x[k] * hsk_726[k];

        t_907[k] = f_10 * hsi0_567[k]
                   - f_11 * hsi1_567[k]
                   + f_3 * pc_x[k] * hsk_727[k];

        t_908[k] = f_3 * pc_y[k] * hsk_725[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, hsi0_569, hsi0_570, hsi0_571, hsi1_569, \
                         hsi1_570, hsi1_571, hsk_729, hsk_730, \
                         hsk_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_10 * hsi0_569[k]
                   - f_11 * hsi1_569[k]
                   + f_3 * pc_x[k] * hsk_729[k];

        t_910[k] = f_8 * hsi0_570[k]
                   - f_9 * hsi1_570[k]
                   + f_3 * pc_x[k] * hsk_730[k];

        t_911[k] = f_8 * hsi0_571[k]
                   - f_9 * hsi1_571[k]
                   + f_3 * pc_x[k] * hsk_731[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, t_915, pc_x, pc_y, hsi0_572, hsi0_574, hsi0_575, \
                         hsi1_572, hsi1_574, hsi1_575, hsk_729, hsk_732, hsk_734, \
                         hsk_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_8 * hsi0_572[k]
                   - f_9 * hsi1_572[k]
                   + f_3 * pc_x[k] * hsk_732[k];

        t_913[k] = f_3 * pc_y[k] * hsk_729[k];

        t_914[k] = f_8 * hsi0_574[k]
                   - f_9 * hsi1_574[k]
                   + f_3 * pc_x[k] * hsk_734[k];

        t_915[k] = f_6 * hsi0_575[k]
                   - f_7 * hsi1_575[k]
                   + f_3 * pc_x[k] * hsk_735[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, t_919, pc_x, pc_y, hsi0_576, hsi0_577, hsi0_578, \
                         hsi1_576, hsi1_577, hsi1_578, hsk_734, hsk_736, hsk_737, \
                         hsk_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_6 * hsi0_576[k]
                   - f_7 * hsi1_576[k]
                   + f_3 * pc_x[k] * hsk_736[k];

        t_917[k] = f_6 * hsi0_577[k]
                   - f_7 * hsi1_577[k]
                   + f_3 * pc_x[k] * hsk_737[k];

        t_918[k] = f_6 * hsi0_578[k]
                   - f_7 * hsi1_578[k]
                   + f_3 * pc_x[k] * hsk_738[k];

        t_919[k] = f_3 * pc_y[k] * hsk_734[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, pc_x, hsi0_580, hsi0_581, hsi0_582, hsi1_580, \
                         hsi1_581, hsi1_582, hsk_740, hsk_741, \
                         hsk_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_6 * hsi0_580[k]
                   - f_7 * hsi1_580[k]
                   + f_3 * pc_x[k] * hsk_740[k];

        t_921[k] = f_4 * hsi0_581[k]
                   - f_5 * hsi1_581[k]
                   + f_3 * pc_x[k] * hsk_741[k];

        t_922[k] = f_4 * hsi0_582[k]
                   - f_5 * hsi1_582[k]
                   + f_3 * pc_x[k] * hsk_742[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, t_926, pc_x, pc_y, hsi0_583, hsi0_584, hsi0_585, \
                         hsi1_583, hsi1_584, hsi1_585, hsk_740, hsk_743, hsk_744, \
                         hsk_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = f_4 * hsi0_583[k]
                   - f_5 * hsi1_583[k]
                   + f_3 * pc_x[k] * hsk_743[k];

        t_924[k] = f_4 * hsi0_584[k]
                   - f_5 * hsi1_584[k]
                   + f_3 * pc_x[k] * hsk_744[k];

        t_925[k] = f_4 * hsi0_585[k]
                   - f_5 * hsi1_585[k]
                   + f_3 * pc_x[k] * hsk_745[k];

        t_926[k] = f_3 * pc_y[k] * hsk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, t_931, t_932, pc_x, hsi0_587, hsi1_587, \
                         hsk_747, hsk_748, hsk_749, hsk_750, hsk_751, \
                         hsk_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_4 * hsi0_587[k]
                   - f_5 * hsi1_587[k]
                   + f_3 * pc_x[k] * hsk_747[k];

        t_928[k] = f_3 * pc_x[k] * hsk_748[k];

        t_929[k] = f_3 * pc_x[k] * hsk_749[k];

        t_930[k] = f_3 * pc_x[k] * hsk_750[k];

        t_931[k] = f_3 * pc_x[k] * hsk_751[k];

        t_932[k] = f_3 * pc_x[k] * hsk_752[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, t_937, pc_x, pc_y, hsi0_581, hsi0_582, \
                         hsi1_581, hsi1_582, hsk_748, hsk_749, hsk_753, hsk_754, \
                         hsk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_3 * pc_x[k] * hsk_753[k];

        t_934[k] = f_3 * pc_x[k] * hsk_754[k];

        t_935[k] = f_3 * pc_x[k] * hsk_755[k];

        t_936[k] = f_1 * hsi0_581[k]
                   - f_2 * hsi1_581[k]
                   + f_3 * pc_y[k] * hsk_748[k];

        t_937[k] = f_20 * hsi0_582[k]
                   - f_21 * hsi1_582[k]
                   + f_3 * pc_y[k] * hsk_749[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, pc_y, hsi0_583, hsi0_584, hsi0_585, hsi1_583, \
                         hsi1_584, hsi1_585, hsk_750, hsk_751, \
                         hsk_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_12 * hsi0_583[k]
                   - f_13 * hsi1_583[k]
                   + f_3 * pc_y[k] * hsk_750[k];

        t_939[k] = f_10 * hsi0_584[k]
                   - f_11 * hsi1_584[k]
                   + f_3 * pc_y[k] * hsk_751[k];

        t_940[k] = f_8 * hsi0_585[k]
                   - f_9 * hsi1_585[k]
                   + f_3 * pc_y[k] * hsk_752[k];
    }

#pragma omp simd aligned(t_941, t_942, t_943, t_944, pc_y, pc_z, gsk_539, hsi0_586, hsi0_587, \
                         hsi1_586, hsi1_587, hsk_753, hsk_754, \
                         hsk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_941[k] = f_6 * hsi0_586[k]
                   - f_7 * hsi1_586[k]
                   + f_3 * pc_y[k] * hsk_753[k];

        t_942[k] = f_4 * hsi0_587[k]
                   - f_5 * hsi1_587[k]
                   + f_3 * pc_y[k] * hsk_754[k];

        t_943[k] = f_3 * pc_y[k] * hsk_755[k];

        t_944[k] = f_0 * gsk_539[k]
                   + f_1 * hsi0_587[k]
                   - f_2 * hsi1_587[k]
                   + f_3 * pc_z[k] * hsk_755[k];
    }
}

auto
compute_prim_hsl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t gsl0, const size_t gsk,
                                                   const size_t gsl1, const size_t hsi0,
                                                   const size_t hsi1, const size_t hsk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_hsl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, gsl0, gsk,
                                                              gsl1, hsi0, hsi1, hsk, ncols,
                                                              gamma, p, q);

    compute_prim_hsl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, gsl0, gsk,
                                                              gsl1, hsi0, hsi1, hsk, ncols,
                                                              gamma, p, q);

    compute_prim_hsl_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, gsl0, gsk,
                                                              gsl1, hsi0, hsi1, hsk, ncols,
                                                              gamma, p, q);

    compute_prim_hsl_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, gsl0, gsk,
                                                              gsl1, hsi0, hsi1, hsk, ncols,
                                                              gamma, p, q);

    compute_prim_hsl_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, gsl0, gsk,
                                                              gsl1, hsk, ncols, gamma, p, q);

    compute_prim_hsl_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, gsl0, gsk,
                                                              gsl1, hsi0, hsi1, hsk, ncols,
                                                              gamma, p, q);

    compute_prim_hsl_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, gsl0, gsk,
                                                              gsl1, hsi0, hsi1, hsk, ncols,
                                                              gamma, p, q);

    compute_prim_hsl_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, gsl0, gsk,
                                                              gsl1, hsi0, hsi1, hsk, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
