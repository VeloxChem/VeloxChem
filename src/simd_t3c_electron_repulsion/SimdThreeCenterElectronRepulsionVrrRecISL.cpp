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


#include "SimdThreeCenterElectronRepulsionVrrRecISL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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

    const auto *hsl0_0 = buffer.data(hsl0 + 0);
    const auto *hsl0_3 = buffer.data(hsl0 + 3);
    const auto *hsl0_5 = buffer.data(hsl0 + 5);
    const auto *hsl0_6 = buffer.data(hsl0 + 6);
    const auto *hsl0_9 = buffer.data(hsl0 + 9);
    const auto *hsl0_10 = buffer.data(hsl0 + 10);
    const auto *hsl0_14 = buffer.data(hsl0 + 14);
    const auto *hsl0_15 = buffer.data(hsl0 + 15);
    const auto *hsl0_20 = buffer.data(hsl0 + 20);
    const auto *hsl0_21 = buffer.data(hsl0 + 21);
    const auto *hsl0_27 = buffer.data(hsl0 + 27);
    const auto *hsl0_36 = buffer.data(hsl0 + 36);
    const auto *hsl0_44 = buffer.data(hsl0 + 44);

    const auto *hsk_0 = buffer.data(hsk + 0);
    const auto *hsk_1 = buffer.data(hsk + 1);
    const auto *hsk_2 = buffer.data(hsk + 2);
    const auto *hsk_3 = buffer.data(hsk + 3);
    const auto *hsk_5 = buffer.data(hsk + 5);
    const auto *hsk_6 = buffer.data(hsk + 6);
    const auto *hsk_9 = buffer.data(hsk + 9);
    const auto *hsk_10 = buffer.data(hsk + 10);
    const auto *hsk_14 = buffer.data(hsk + 14);
    const auto *hsk_15 = buffer.data(hsk + 15);
    const auto *hsk_20 = buffer.data(hsk + 20);
    const auto *hsk_28 = buffer.data(hsk + 28);
    const auto *hsk_30 = buffer.data(hsk + 30);
    const auto *hsk_31 = buffer.data(hsk + 31);
    const auto *hsk_32 = buffer.data(hsk + 32);
    const auto *hsk_33 = buffer.data(hsk + 33);
    const auto *hsk_35 = buffer.data(hsk + 35);
    const auto *hsk_64 = buffer.data(hsk + 64);
    const auto *hsk_66 = buffer.data(hsk + 66);
    const auto *hsk_67 = buffer.data(hsk + 67);
    const auto *hsk_68 = buffer.data(hsk + 68);
    const auto *hsk_69 = buffer.data(hsk + 69);
    const auto *hsk_70 = buffer.data(hsk + 70);
    const auto *hsk_71 = buffer.data(hsk + 71);
    const auto *hsk_100 = buffer.data(hsk + 100);
    const auto *hsk_101 = buffer.data(hsk + 101);
    const auto *hsk_102 = buffer.data(hsk + 102);
    const auto *hsk_103 = buffer.data(hsk + 103);
    const auto *hsk_104 = buffer.data(hsk + 104);
    const auto *hsk_105 = buffer.data(hsk + 105);
    const auto *hsk_107 = buffer.data(hsk + 107);

    const auto *hsl1_0 = buffer.data(hsl1 + 0);
    const auto *hsl1_3 = buffer.data(hsl1 + 3);
    const auto *hsl1_5 = buffer.data(hsl1 + 5);
    const auto *hsl1_6 = buffer.data(hsl1 + 6);
    const auto *hsl1_9 = buffer.data(hsl1 + 9);
    const auto *hsl1_10 = buffer.data(hsl1 + 10);
    const auto *hsl1_14 = buffer.data(hsl1 + 14);
    const auto *hsl1_15 = buffer.data(hsl1 + 15);
    const auto *hsl1_20 = buffer.data(hsl1 + 20);
    const auto *hsl1_21 = buffer.data(hsl1 + 21);
    const auto *hsl1_27 = buffer.data(hsl1 + 27);
    const auto *hsl1_36 = buffer.data(hsl1 + 36);
    const auto *hsl1_44 = buffer.data(hsl1 + 44);

    const auto *isi0_0 = buffer.data(isi0 + 0);
    const auto *isi0_1 = buffer.data(isi0 + 1);
    const auto *isi0_2 = buffer.data(isi0 + 2);
    const auto *isi0_3 = buffer.data(isi0 + 3);
    const auto *isi0_5 = buffer.data(isi0 + 5);
    const auto *isi0_6 = buffer.data(isi0 + 6);
    const auto *isi0_8 = buffer.data(isi0 + 8);
    const auto *isi0_9 = buffer.data(isi0 + 9);
    const auto *isi0_10 = buffer.data(isi0 + 10);
    const auto *isi0_12 = buffer.data(isi0 + 12);
    const auto *isi0_13 = buffer.data(isi0 + 13);
    const auto *isi0_14 = buffer.data(isi0 + 14);
    const auto *isi0_21 = buffer.data(isi0 + 21);
    const auto *isi0_23 = buffer.data(isi0 + 23);
    const auto *isi0_24 = buffer.data(isi0 + 24);
    const auto *isi0_25 = buffer.data(isi0 + 25);
    const auto *isi0_26 = buffer.data(isi0 + 26);
    const auto *isi0_27 = buffer.data(isi0 + 27);
    const auto *isi0_31 = buffer.data(isi0 + 31);
    const auto *isi0_34 = buffer.data(isi0 + 34);
    const auto *isi0_35 = buffer.data(isi0 + 35);
    const auto *isi0_38 = buffer.data(isi0 + 38);
    const auto *isi0_39 = buffer.data(isi0 + 39);
    const auto *isi0_40 = buffer.data(isi0 + 40);
    const auto *isi0_49 = buffer.data(isi0 + 49);
    const auto *isi0_50 = buffer.data(isi0 + 50);
    const auto *isi0_51 = buffer.data(isi0 + 51);
    const auto *isi0_52 = buffer.data(isi0 + 52);
    const auto *isi0_53 = buffer.data(isi0 + 53);
    const auto *isi0_58 = buffer.data(isi0 + 58);
    const auto *isi0_60 = buffer.data(isi0 + 60);
    const auto *isi0_61 = buffer.data(isi0 + 61);
    const auto *isi0_63 = buffer.data(isi0 + 63);
    const auto *isi0_64 = buffer.data(isi0 + 64);
    const auto *isi0_65 = buffer.data(isi0 + 65);
    const auto *isi0_67 = buffer.data(isi0 + 67);
    const auto *isi0_68 = buffer.data(isi0 + 68);
    const auto *isi0_69 = buffer.data(isi0 + 69);
    const auto *isi0_70 = buffer.data(isi0 + 70);
    const auto *isi0_78 = buffer.data(isi0 + 78);
    const auto *isi0_79 = buffer.data(isi0 + 79);

    const auto *isi1_0 = buffer.data(isi1 + 0);
    const auto *isi1_1 = buffer.data(isi1 + 1);
    const auto *isi1_2 = buffer.data(isi1 + 2);
    const auto *isi1_3 = buffer.data(isi1 + 3);
    const auto *isi1_5 = buffer.data(isi1 + 5);
    const auto *isi1_6 = buffer.data(isi1 + 6);
    const auto *isi1_8 = buffer.data(isi1 + 8);
    const auto *isi1_9 = buffer.data(isi1 + 9);
    const auto *isi1_10 = buffer.data(isi1 + 10);
    const auto *isi1_12 = buffer.data(isi1 + 12);
    const auto *isi1_13 = buffer.data(isi1 + 13);
    const auto *isi1_14 = buffer.data(isi1 + 14);
    const auto *isi1_21 = buffer.data(isi1 + 21);
    const auto *isi1_23 = buffer.data(isi1 + 23);
    const auto *isi1_24 = buffer.data(isi1 + 24);
    const auto *isi1_25 = buffer.data(isi1 + 25);
    const auto *isi1_26 = buffer.data(isi1 + 26);
    const auto *isi1_27 = buffer.data(isi1 + 27);
    const auto *isi1_31 = buffer.data(isi1 + 31);
    const auto *isi1_34 = buffer.data(isi1 + 34);
    const auto *isi1_35 = buffer.data(isi1 + 35);
    const auto *isi1_38 = buffer.data(isi1 + 38);
    const auto *isi1_39 = buffer.data(isi1 + 39);
    const auto *isi1_40 = buffer.data(isi1 + 40);
    const auto *isi1_49 = buffer.data(isi1 + 49);
    const auto *isi1_50 = buffer.data(isi1 + 50);
    const auto *isi1_51 = buffer.data(isi1 + 51);
    const auto *isi1_52 = buffer.data(isi1 + 52);
    const auto *isi1_53 = buffer.data(isi1 + 53);
    const auto *isi1_58 = buffer.data(isi1 + 58);
    const auto *isi1_60 = buffer.data(isi1 + 60);
    const auto *isi1_61 = buffer.data(isi1 + 61);
    const auto *isi1_63 = buffer.data(isi1 + 63);
    const auto *isi1_64 = buffer.data(isi1 + 64);
    const auto *isi1_65 = buffer.data(isi1 + 65);
    const auto *isi1_67 = buffer.data(isi1 + 67);
    const auto *isi1_68 = buffer.data(isi1 + 68);
    const auto *isi1_69 = buffer.data(isi1 + 69);
    const auto *isi1_70 = buffer.data(isi1 + 70);
    const auto *isi1_78 = buffer.data(isi1 + 78);
    const auto *isi1_79 = buffer.data(isi1 + 79);

    const auto *isk_0 = buffer.data(isk + 0);
    const auto *isk_1 = buffer.data(isk + 1);
    const auto *isk_2 = buffer.data(isk + 2);
    const auto *isk_3 = buffer.data(isk + 3);
    const auto *isk_5 = buffer.data(isk + 5);
    const auto *isk_6 = buffer.data(isk + 6);
    const auto *isk_8 = buffer.data(isk + 8);
    const auto *isk_9 = buffer.data(isk + 9);
    const auto *isk_10 = buffer.data(isk + 10);
    const auto *isk_12 = buffer.data(isk + 12);
    const auto *isk_13 = buffer.data(isk + 13);
    const auto *isk_14 = buffer.data(isk + 14);
    const auto *isk_15 = buffer.data(isk + 15);
    const auto *isk_17 = buffer.data(isk + 17);
    const auto *isk_18 = buffer.data(isk + 18);
    const auto *isk_19 = buffer.data(isk + 19);
    const auto *isk_20 = buffer.data(isk + 20);
    const auto *isk_21 = buffer.data(isk + 21);
    const auto *isk_27 = buffer.data(isk + 27);
    const auto *isk_28 = buffer.data(isk + 28);
    const auto *isk_30 = buffer.data(isk + 30);
    const auto *isk_31 = buffer.data(isk + 31);
    const auto *isk_32 = buffer.data(isk + 32);
    const auto *isk_33 = buffer.data(isk + 33);
    const auto *isk_34 = buffer.data(isk + 34);
    const auto *isk_35 = buffer.data(isk + 35);
    const auto *isk_36 = buffer.data(isk + 36);
    const auto *isk_37 = buffer.data(isk + 37);
    const auto *isk_39 = buffer.data(isk + 39);
    const auto *isk_41 = buffer.data(isk + 41);
    const auto *isk_42 = buffer.data(isk + 42);
    const auto *isk_43 = buffer.data(isk + 43);
    const auto *isk_45 = buffer.data(isk + 45);
    const auto *isk_46 = buffer.data(isk + 46);
    const auto *isk_47 = buffer.data(isk + 47);
    const auto *isk_48 = buffer.data(isk + 48);
    const auto *isk_50 = buffer.data(isk + 50);
    const auto *isk_51 = buffer.data(isk + 51);
    const auto *isk_52 = buffer.data(isk + 52);
    const auto *isk_53 = buffer.data(isk + 53);
    const auto *isk_54 = buffer.data(isk + 54);
    const auto *isk_56 = buffer.data(isk + 56);
    const auto *isk_57 = buffer.data(isk + 57);
    const auto *isk_64 = buffer.data(isk + 64);
    const auto *isk_65 = buffer.data(isk + 65);
    const auto *isk_66 = buffer.data(isk + 66);
    const auto *isk_67 = buffer.data(isk + 67);
    const auto *isk_68 = buffer.data(isk + 68);
    const auto *isk_69 = buffer.data(isk + 69);
    const auto *isk_70 = buffer.data(isk + 70);
    const auto *isk_71 = buffer.data(isk + 71);
    const auto *isk_72 = buffer.data(isk + 72);
    const auto *isk_74 = buffer.data(isk + 74);
    const auto *isk_76 = buffer.data(isk + 76);
    const auto *isk_77 = buffer.data(isk + 77);
    const auto *isk_79 = buffer.data(isk + 79);
    const auto *isk_80 = buffer.data(isk + 80);
    const auto *isk_81 = buffer.data(isk + 81);
    const auto *isk_83 = buffer.data(isk + 83);
    const auto *isk_84 = buffer.data(isk + 84);
    const auto *isk_85 = buffer.data(isk + 85);
    const auto *isk_86 = buffer.data(isk + 86);
    const auto *isk_88 = buffer.data(isk + 88);
    const auto *isk_89 = buffer.data(isk + 89);
    const auto *isk_90 = buffer.data(isk + 90);
    const auto *isk_91 = buffer.data(isk + 91);
    const auto *isk_92 = buffer.data(isk + 92);
    const auto *isk_99 = buffer.data(isk + 99);
    const auto *isk_100 = buffer.data(isk + 100);
    const auto *isk_101 = buffer.data(isk + 101);
    const auto *isk_102 = buffer.data(isk + 102);
    const auto *isk_103 = buffer.data(isk + 103);
    const auto *isk_104 = buffer.data(isk + 104);
    const auto *isk_105 = buffer.data(isk + 105);
    const auto *isk_107 = buffer.data(isk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, hsk_0, isi0_0, \
                         isi1_0, isk_0, isk_1, isk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hsk_0[k]
                 + f_1 * isi0_0[k]
                 - f_2 * isi1_0[k]
                 + f_3 * pc_x[k] * isk_0[k];

        t_1[k] = f_3 * pc_y[k] * isk_0[k];

        t_2[k] = f_3 * pc_z[k] * isk_0[k];

        t_3[k] = f_4 * isi0_0[k]
                 - f_5 * isi1_0[k]
                 + f_3 * pc_y[k] * isk_1[k];

        t_4[k] = f_3 * pc_y[k] * isk_2[k];

        t_5[k] = f_4 * isi0_0[k]
                 - f_5 * isi1_0[k]
                 + f_3 * pc_z[k] * isk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, isi0_1, isi0_2, isi0_3, isi1_1, \
                         isi1_2, isi1_3, isk_3, isk_5, isk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * isi0_1[k]
                 - f_7 * isi1_1[k]
                 + f_3 * pc_y[k] * isk_3[k];

        t_7[k] = f_3 * pc_z[k] * isk_3[k];

        t_8[k] = f_3 * pc_y[k] * isk_5[k];

        t_9[k] = f_6 * isi0_2[k]
                 - f_7 * isi1_2[k]
                 + f_3 * pc_z[k] * isk_5[k];

        t_10[k] = f_8 * isi0_3[k]
                  - f_9 * isi1_3[k]
                  + f_3 * pc_y[k] * isk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, isi0_5, isi0_6, \
                         isi1_5, isi1_6, isk_6, isk_8, isk_9, isk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * isk_6[k];

        t_12[k] = f_4 * isi0_5[k]
                  - f_5 * isi1_5[k]
                  + f_3 * pc_y[k] * isk_8[k];

        t_13[k] = f_3 * pc_y[k] * isk_9[k];

        t_14[k] = f_8 * isi0_5[k]
                  - f_9 * isi1_5[k]
                  + f_3 * pc_z[k] * isk_9[k];

        t_15[k] = f_10 * isi0_6[k]
                  - f_11 * isi1_6[k]
                  + f_3 * pc_y[k] * isk_10[k];

        t_16[k] = f_3 * pc_z[k] * isk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, isi0_8, isi0_9, isi1_8, isi1_9, \
                         isk_12, isk_13, isk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * isi0_8[k]
                  - f_7 * isi1_8[k]
                  + f_3 * pc_y[k] * isk_12[k];

        t_18[k] = f_4 * isi0_9[k]
                  - f_5 * isi1_9[k]
                  + f_3 * pc_y[k] * isk_13[k];

        t_19[k] = f_3 * pc_y[k] * isk_14[k];

        t_20[k] = f_10 * isi0_9[k]
                  - f_11 * isi1_9[k]
                  + f_3 * pc_z[k] * isk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, isi0_10, isi0_12, isi0_13, \
                         isi1_10, isi1_12, isi1_13, isk_15, isk_17, \
                         isk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_12 * isi0_10[k]
                  - f_13 * isi1_10[k]
                  + f_3 * pc_y[k] * isk_15[k];

        t_22[k] = f_3 * pc_z[k] * isk_15[k];

        t_23[k] = f_8 * isi0_12[k]
                  - f_9 * isi1_12[k]
                  + f_3 * pc_y[k] * isk_17[k];

        t_24[k] = f_6 * isi0_13[k]
                  - f_7 * isi1_13[k]
                  + f_3 * pc_y[k] * isk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, hsk_28, isi0_14, \
                         isi1_14, isk_19, isk_20, isk_21, isk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * isi0_14[k]
                  - f_5 * isi1_14[k]
                  + f_3 * pc_y[k] * isk_19[k];

        t_26[k] = f_3 * pc_y[k] * isk_20[k];

        t_27[k] = f_12 * isi0_14[k]
                  - f_13 * isi1_14[k]
                  + f_3 * pc_z[k] * isk_20[k];

        t_28[k] = f_0 * hsk_28[k]
                  + f_3 * pc_x[k] * isk_28[k];

        t_29[k] = f_3 * pc_z[k] * isk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, hsk_30, hsk_31, hsk_32, \
                         hsk_33, isk_27, isk_30, isk_31, isk_32, \
                         isk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * hsk_30[k]
                  + f_3 * pc_x[k] * isk_30[k];

        t_31[k] = f_0 * hsk_31[k]
                  + f_3 * pc_x[k] * isk_31[k];

        t_32[k] = f_0 * hsk_32[k]
                  + f_3 * pc_x[k] * isk_32[k];

        t_33[k] = f_0 * hsk_33[k]
                  + f_3 * pc_x[k] * isk_33[k];

        t_34[k] = f_3 * pc_y[k] * isk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, hsk_35, isi0_21, isi0_23, \
                         isi1_21, isi1_23, isk_28, isk_30, isk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * hsk_35[k]
                  + f_3 * pc_x[k] * isk_35[k];

        t_36[k] = f_1 * isi0_21[k]
                  - f_2 * isi1_21[k]
                  + f_3 * pc_y[k] * isk_28[k];

        t_37[k] = f_3 * pc_z[k] * isk_28[k];

        t_38[k] = f_12 * isi0_23[k]
                  - f_13 * isi1_23[k]
                  + f_3 * pc_y[k] * isk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pc_y, isi0_24, isi0_25, isi0_26, isi1_24, isi1_25, \
                         isi1_26, isk_31, isk_32, isk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * isi0_24[k]
                  - f_11 * isi1_24[k]
                  + f_3 * pc_y[k] * isk_31[k];

        t_40[k] = f_8 * isi0_25[k]
                  - f_9 * isi1_25[k]
                  + f_3 * pc_y[k] * isk_32[k];

        t_41[k] = f_6 * isi0_26[k]
                  - f_7 * isi1_26[k]
                  + f_3 * pc_y[k] * isk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_y, pc_z, hsl0_0, hsk_0, \
                         hsl1_0, isi0_27, isi1_27, isk_34, isk_35, \
                         isk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * isi0_27[k]
                  - f_5 * isi1_27[k]
                  + f_3 * pc_y[k] * isk_34[k];

        t_43[k] = f_3 * pc_y[k] * isk_35[k];

        t_44[k] = f_1 * isi0_27[k]
                  - f_2 * isi1_27[k]
                  + f_3 * pc_z[k] * isk_35[k];

        t_45[k] = pa_y[k] * hsl0_0[k]
                  - f_14 * pc_y[k] * hsl1_0[k];

        t_46[k] = f_15 * hsk_0[k]
                  + f_3 * pc_y[k] * isk_36[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_y, pc_z, hsl0_3, hsl0_5, hsk_1, \
                         hsl1_3, hsl1_5, isk_36, isk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * isk_36[k];

        t_48[k] = pa_y[k] * hsl0_3[k]
                  + f_16 * hsk_1[k]
                  - f_14 * pc_y[k] * hsl1_3[k];

        t_49[k] = f_3 * pc_z[k] * isk_37[k];

        t_50[k] = pa_y[k] * hsl0_5[k]
                  - f_14 * pc_y[k] * hsl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pc_y, pc_z, hsl0_6, hsl0_9, hsk_3, \
                         hsk_5, hsl1_6, hsl1_9, isk_39, isk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * hsl0_6[k]
                  + f_17 * hsk_3[k]
                  - f_14 * pc_y[k] * hsl1_6[k];

        t_52[k] = f_3 * pc_z[k] * isk_39[k];

        t_53[k] = f_15 * hsk_5[k]
                  + f_3 * pc_y[k] * isk_41[k];

        t_54[k] = pa_y[k] * hsl0_9[k]
                  - f_14 * pc_y[k] * hsl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pc_y, pc_z, hsl0_10, hsk_6, hsk_9, \
                         hsl1_10, isi0_31, isi1_31, isk_42, isk_43, \
                         isk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * hsl0_10[k]
                  + f_18 * hsk_6[k]
                  - f_14 * pc_y[k] * hsl1_10[k];

        t_56[k] = f_3 * pc_z[k] * isk_42[k];

        t_57[k] = f_4 * isi0_31[k]
                  - f_5 * isi1_31[k]
                  + f_3 * pc_z[k] * isk_43[k];

        t_58[k] = f_15 * hsk_9[k]
                  + f_3 * pc_y[k] * isk_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, hsl0_14, hsl0_15, hsk_10, \
                         hsl1_14, hsl1_15, isi0_34, isi1_34, isk_46, \
                         isk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * hsl0_14[k]
                  - f_14 * pc_y[k] * hsl1_14[k];

        t_60[k] = pa_y[k] * hsl0_15[k]
                  + f_19 * hsk_10[k]
                  - f_14 * pc_y[k] * hsl1_15[k];

        t_61[k] = f_3 * pc_z[k] * isk_46[k];

        t_62[k] = f_4 * isi0_34[k]
                  - f_5 * isi1_34[k]
                  + f_3 * pc_z[k] * isk_47[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pc_y, pc_z, hsl0_20, hsk_14, hsl1_20, \
                         isi0_35, isi1_35, isk_48, isk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_6 * isi0_35[k]
                  - f_7 * isi1_35[k]
                  + f_3 * pc_z[k] * isk_48[k];

        t_64[k] = f_15 * hsk_14[k]
                  + f_3 * pc_y[k] * isk_50[k];

        t_65[k] = pa_y[k] * hsl0_20[k]
                  - f_14 * pc_y[k] * hsl1_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pc_y, pc_z, hsl0_21, hsk_15, hsl1_21, \
                         isi0_38, isi1_38, isk_51, isk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * hsl0_21[k]
                  + f_0 * hsk_15[k]
                  - f_14 * pc_y[k] * hsl1_21[k];

        t_67[k] = f_3 * pc_z[k] * isk_51[k];

        t_68[k] = f_4 * isi0_38[k]
                  - f_5 * isi1_38[k]
                  + f_3 * pc_z[k] * isk_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_y, pc_z, hsk_20, isi0_39, isi0_40, isi1_39, \
                         isi1_40, isk_53, isk_54, isk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * isi0_39[k]
                  - f_7 * isi1_39[k]
                  + f_3 * pc_z[k] * isk_53[k];

        t_70[k] = f_8 * isi0_40[k]
                  - f_9 * isi1_40[k]
                  + f_3 * pc_z[k] * isk_54[k];

        t_71[k] = f_15 * hsk_20[k]
                  + f_3 * pc_y[k] * isk_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, pc_z, hsl0_27, hsk_64, \
                         hsk_66, hsl1_27, isk_57, isk_64, isk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * hsl0_27[k]
                  - f_14 * pc_y[k] * hsl1_27[k];

        t_73[k] = f_19 * hsk_64[k]
                  + f_3 * pc_x[k] * isk_64[k];

        t_74[k] = f_3 * pc_z[k] * isk_57[k];

        t_75[k] = f_19 * hsk_66[k]
                  + f_3 * pc_x[k] * isk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, hsk_67, hsk_68, hsk_69, hsk_70, \
                         hsk_71, isk_67, isk_68, isk_69, isk_70, \
                         isk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_19 * hsk_67[k]
                  + f_3 * pc_x[k] * isk_67[k];

        t_77[k] = f_19 * hsk_68[k]
                  + f_3 * pc_x[k] * isk_68[k];

        t_78[k] = f_19 * hsk_69[k]
                  + f_3 * pc_x[k] * isk_69[k];

        t_79[k] = f_19 * hsk_70[k]
                  + f_3 * pc_x[k] * isk_70[k];

        t_80[k] = f_19 * hsk_71[k]
                  + f_3 * pc_x[k] * isk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_y, pc_z, hsk_28, isi0_49, isi0_50, \
                         isi1_49, isi1_50, isk_64, isk_65, isk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_15 * hsk_28[k]
                  + f_1 * isi0_49[k]
                  - f_2 * isi1_49[k]
                  + f_3 * pc_y[k] * isk_64[k];

        t_82[k] = f_3 * pc_z[k] * isk_64[k];

        t_83[k] = f_4 * isi0_49[k]
                  - f_5 * isi1_49[k]
                  + f_3 * pc_z[k] * isk_65[k];

        t_84[k] = f_6 * isi0_50[k]
                  - f_7 * isi1_50[k]
                  + f_3 * pc_z[k] * isk_66[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_z, isi0_51, isi0_52, isi0_53, isi1_51, isi1_52, \
                         isi1_53, isk_67, isk_68, isk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_8 * isi0_51[k]
                  - f_9 * isi1_51[k]
                  + f_3 * pc_z[k] * isk_67[k];

        t_86[k] = f_10 * isi0_52[k]
                  - f_11 * isi1_52[k]
                  + f_3 * pc_z[k] * isk_68[k];

        t_87[k] = f_12 * isi0_53[k]
                  - f_13 * isi1_53[k]
                  + f_3 * pc_z[k] * isk_69[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pc_y, pc_z, hsl0_0, hsl0_44, \
                         hsk_35, hsl1_0, hsl1_44, isk_71, isk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * hsk_35[k]
                  + f_3 * pc_y[k] * isk_71[k];

        t_89[k] = pa_y[k] * hsl0_44[k]
                  - f_14 * pc_y[k] * hsl1_44[k];

        t_90[k] = pa_z[k] * hsl0_0[k]
                  - f_14 * pc_z[k] * hsl1_0[k];

        t_91[k] = f_3 * pc_y[k] * isk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_y, pc_z, hsl0_3, hsl0_5, hsk_0, \
                         hsk_2, hsl1_3, hsl1_5, isk_72, isk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * hsk_0[k]
                  + f_3 * pc_z[k] * isk_72[k];

        t_93[k] = pa_z[k] * hsl0_3[k]
                  - f_14 * pc_z[k] * hsl1_3[k];

        t_94[k] = f_3 * pc_y[k] * isk_74[k];

        t_95[k] = pa_z[k] * hsl0_5[k]
                  + f_16 * hsk_2[k]
                  - f_14 * pc_z[k] * hsl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pc_y, pc_z, hsl0_6, hsl0_9, hsk_5, \
                         hsl1_6, hsl1_9, isi0_58, isi1_58, isk_76, \
                         isk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * hsl0_6[k]
                  - f_14 * pc_z[k] * hsl1_6[k];

        t_97[k] = f_4 * isi0_58[k]
                  - f_5 * isi1_58[k]
                  + f_3 * pc_y[k] * isk_76[k];

        t_98[k] = f_3 * pc_y[k] * isk_77[k];

        t_99[k] = pa_z[k] * hsl0_9[k]
                  + f_17 * hsk_5[k]
                  - f_14 * pc_z[k] * hsl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pc_y, pc_z, hsl0_10, hsl1_10, \
                         isi0_60, isi0_61, isi1_60, isi1_61, isk_79, isk_80, \
                         isk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * hsl0_10[k]
                   - f_14 * pc_z[k] * hsl1_10[k];

        t_101[k] = f_6 * isi0_60[k]
                   - f_7 * isi1_60[k]
                   + f_3 * pc_y[k] * isk_79[k];

        t_102[k] = f_4 * isi0_61[k]
                   - f_5 * isi1_61[k]
                   + f_3 * pc_y[k] * isk_80[k];

        t_103[k] = f_3 * pc_y[k] * isk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pc_y, pc_z, hsl0_14, hsl0_15, hsk_9, \
                         hsl1_14, hsl1_15, isi0_63, isi1_63, isk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_z[k] * hsl0_14[k]
                   + f_18 * hsk_9[k]
                   - f_14 * pc_z[k] * hsl1_14[k];

        t_105[k] = pa_z[k] * hsl0_15[k]
                   - f_14 * pc_z[k] * hsl1_15[k];

        t_106[k] = f_8 * isi0_63[k]
                   - f_9 * isi1_63[k]
                   + f_3 * pc_y[k] * isk_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, isi0_64, isi0_65, isi1_64, isi1_65, \
                         isk_84, isk_85, isk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * isi0_64[k]
                   - f_7 * isi1_64[k]
                   + f_3 * pc_y[k] * isk_84[k];

        t_108[k] = f_4 * isi0_65[k]
                   - f_5 * isi1_65[k]
                   + f_3 * pc_y[k] * isk_85[k];

        t_109[k] = f_3 * pc_y[k] * isk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pc_y, pc_z, hsl0_20, hsl0_21, hsk_14, \
                         hsl1_20, hsl1_21, isi0_67, isi1_67, isk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * hsl0_20[k]
                   + f_19 * hsk_14[k]
                   - f_14 * pc_z[k] * hsl1_20[k];

        t_111[k] = pa_z[k] * hsl0_21[k]
                   - f_14 * pc_z[k] * hsl1_21[k];

        t_112[k] = f_10 * isi0_67[k]
                   - f_11 * isi1_67[k]
                   + f_3 * pc_y[k] * isk_88[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_y, isi0_68, isi0_69, isi0_70, isi1_68, \
                         isi1_69, isi1_70, isk_89, isk_90, isk_91, \
                         isk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * isi0_68[k]
                   - f_9 * isi1_68[k]
                   + f_3 * pc_y[k] * isk_89[k];

        t_114[k] = f_6 * isi0_69[k]
                   - f_7 * isi1_69[k]
                   + f_3 * pc_y[k] * isk_90[k];

        t_115[k] = f_4 * isi0_70[k]
                   - f_5 * isi1_70[k]
                   + f_3 * pc_y[k] * isk_91[k];

        t_116[k] = f_3 * pc_y[k] * isk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pc_x, pc_z, hsl0_27, hsk_20, \
                         hsk_100, hsk_101, hsk_102, hsl1_27, isk_100, isk_101, \
                         isk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_z[k] * hsl0_27[k]
                   + f_0 * hsk_20[k]
                   - f_14 * pc_z[k] * hsl1_27[k];

        t_118[k] = f_19 * hsk_100[k]
                   + f_3 * pc_x[k] * isk_100[k];

        t_119[k] = f_19 * hsk_101[k]
                   + f_3 * pc_x[k] * isk_101[k];

        t_120[k] = f_19 * hsk_102[k]
                   + f_3 * pc_x[k] * isk_102[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, hsk_103, hsk_104, \
                         hsk_105, hsk_107, isk_99, isk_103, isk_104, isk_105, \
                         isk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_19 * hsk_103[k]
                   + f_3 * pc_x[k] * isk_103[k];

        t_122[k] = f_19 * hsk_104[k]
                   + f_3 * pc_x[k] * isk_104[k];

        t_123[k] = f_19 * hsk_105[k]
                   + f_3 * pc_x[k] * isk_105[k];

        t_124[k] = f_3 * pc_y[k] * isk_99[k];

        t_125[k] = f_19 * hsk_107[k]
                   + f_3 * pc_x[k] * isk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_z, pc_y, pc_z, hsl0_36, hsl1_36, isi0_78, \
                         isi0_79, isi1_78, isi1_79, isk_101, isk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * hsl0_36[k]
                   - f_14 * pc_z[k] * hsl1_36[k];

        t_127[k] = f_20 * isi0_78[k]
                   - f_21 * isi1_78[k]
                   + f_3 * pc_y[k] * isk_101[k];

        t_128[k] = f_12 * isi0_79[k]
                   - f_13 * isi1_79[k]
                   + f_3 * pc_y[k] * isk_102[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
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

    const auto *hsl0_48 = buffer.data(hsl0 + 48);
    const auto *hsl0_51 = buffer.data(hsl0 + 51);
    const auto *hsl0_55 = buffer.data(hsl0 + 55);
    const auto *hsl0_60 = buffer.data(hsl0 + 60);
    const auto *hsl0_66 = buffer.data(hsl0 + 66);
    const auto *hsl0_81 = buffer.data(hsl0 + 81);
    const auto *hsl0_90 = buffer.data(hsl0 + 90);
    const auto *hsl0_95 = buffer.data(hsl0 + 95);
    const auto *hsl0_99 = buffer.data(hsl0 + 99);
    const auto *hsl0_102 = buffer.data(hsl0 + 102);
    const auto *hsl0_104 = buffer.data(hsl0 + 104);
    const auto *hsl0_107 = buffer.data(hsl0 + 107);
    const auto *hsl0_108 = buffer.data(hsl0 + 108);
    const auto *hsl0_110 = buffer.data(hsl0 + 110);
    const auto *hsl0_113 = buffer.data(hsl0 + 113);
    const auto *hsl0_114 = buffer.data(hsl0 + 114);
    const auto *hsl0_115 = buffer.data(hsl0 + 115);
    const auto *hsl0_117 = buffer.data(hsl0 + 117);
    const auto *hsl0_134 = buffer.data(hsl0 + 134);

    const auto *hsk_35 = buffer.data(hsk + 35);
    const auto *hsk_36 = buffer.data(hsk + 36);
    const auto *hsk_39 = buffer.data(hsk + 39);
    const auto *hsk_41 = buffer.data(hsk + 41);
    const auto *hsk_42 = buffer.data(hsk + 42);
    const auto *hsk_45 = buffer.data(hsk + 45);
    const auto *hsk_46 = buffer.data(hsk + 46);
    const auto *hsk_50 = buffer.data(hsk + 50);
    const auto *hsk_51 = buffer.data(hsk + 51);
    const auto *hsk_56 = buffer.data(hsk + 56);
    const auto *hsk_64 = buffer.data(hsk + 64);
    const auto *hsk_71 = buffer.data(hsk + 71);
    const auto *hsk_72 = buffer.data(hsk + 72);
    const auto *hsk_74 = buffer.data(hsk + 74);
    const auto *hsk_77 = buffer.data(hsk + 77);
    const auto *hsk_80 = buffer.data(hsk + 80);
    const auto *hsk_81 = buffer.data(hsk + 81);
    const auto *hsk_84 = buffer.data(hsk + 84);
    const auto *hsk_85 = buffer.data(hsk + 85);
    const auto *hsk_86 = buffer.data(hsk + 86);
    const auto *hsk_89 = buffer.data(hsk + 89);
    const auto *hsk_90 = buffer.data(hsk + 90);
    const auto *hsk_91 = buffer.data(hsk + 91);
    const auto *hsk_92 = buffer.data(hsk + 92);
    const auto *hsk_102 = buffer.data(hsk + 102);
    const auto *hsk_103 = buffer.data(hsk + 103);
    const auto *hsk_104 = buffer.data(hsk + 104);
    const auto *hsk_105 = buffer.data(hsk + 105);
    const auto *hsk_106 = buffer.data(hsk + 106);
    const auto *hsk_107 = buffer.data(hsk + 107);
    const auto *hsk_108 = buffer.data(hsk + 108);
    const auto *hsk_111 = buffer.data(hsk + 111);
    const auto *hsk_114 = buffer.data(hsk + 114);
    const auto *hsk_118 = buffer.data(hsk + 118);
    const auto *hsk_123 = buffer.data(hsk + 123);
    const auto *hsk_129 = buffer.data(hsk + 129);
    const auto *hsk_136 = buffer.data(hsk + 136);
    const auto *hsk_138 = buffer.data(hsk + 138);
    const auto *hsk_139 = buffer.data(hsk + 139);
    const auto *hsk_140 = buffer.data(hsk + 140);
    const auto *hsk_141 = buffer.data(hsk + 141);
    const auto *hsk_142 = buffer.data(hsk + 142);
    const auto *hsk_143 = buffer.data(hsk + 143);
    const auto *hsk_172 = buffer.data(hsk + 172);
    const auto *hsk_173 = buffer.data(hsk + 173);
    const auto *hsk_174 = buffer.data(hsk + 174);
    const auto *hsk_175 = buffer.data(hsk + 175);
    const auto *hsk_176 = buffer.data(hsk + 176);
    const auto *hsk_177 = buffer.data(hsk + 177);
    const auto *hsk_178 = buffer.data(hsk + 178);
    const auto *hsk_179 = buffer.data(hsk + 179);
    const auto *hsk_180 = buffer.data(hsk + 180);
    const auto *hsk_185 = buffer.data(hsk + 185);
    const auto *hsk_189 = buffer.data(hsk + 189);
    const auto *hsk_194 = buffer.data(hsk + 194);
    const auto *hsk_200 = buffer.data(hsk + 200);

    const auto *hsl1_48 = buffer.data(hsl1 + 48);
    const auto *hsl1_51 = buffer.data(hsl1 + 51);
    const auto *hsl1_55 = buffer.data(hsl1 + 55);
    const auto *hsl1_60 = buffer.data(hsl1 + 60);
    const auto *hsl1_66 = buffer.data(hsl1 + 66);
    const auto *hsl1_81 = buffer.data(hsl1 + 81);
    const auto *hsl1_90 = buffer.data(hsl1 + 90);
    const auto *hsl1_95 = buffer.data(hsl1 + 95);
    const auto *hsl1_99 = buffer.data(hsl1 + 99);
    const auto *hsl1_102 = buffer.data(hsl1 + 102);
    const auto *hsl1_104 = buffer.data(hsl1 + 104);
    const auto *hsl1_107 = buffer.data(hsl1 + 107);
    const auto *hsl1_108 = buffer.data(hsl1 + 108);
    const auto *hsl1_110 = buffer.data(hsl1 + 110);
    const auto *hsl1_113 = buffer.data(hsl1 + 113);
    const auto *hsl1_114 = buffer.data(hsl1 + 114);
    const auto *hsl1_115 = buffer.data(hsl1 + 115);
    const auto *hsl1_117 = buffer.data(hsl1 + 117);
    const auto *hsl1_134 = buffer.data(hsl1 + 134);

    const auto *isi0_80 = buffer.data(isi0 + 80);
    const auto *isi0_81 = buffer.data(isi0 + 81);
    const auto *isi0_82 = buffer.data(isi0 + 82);
    const auto *isi0_83 = buffer.data(isi0 + 83);
    const auto *isi0_84 = buffer.data(isi0 + 84);
    const auto *isi0_86 = buffer.data(isi0 + 86);
    const auto *isi0_87 = buffer.data(isi0 + 87);
    const auto *isi0_89 = buffer.data(isi0 + 89);
    const auto *isi0_90 = buffer.data(isi0 + 90);
    const auto *isi0_91 = buffer.data(isi0 + 91);
    const auto *isi0_93 = buffer.data(isi0 + 93);
    const auto *isi0_94 = buffer.data(isi0 + 94);
    const auto *isi0_95 = buffer.data(isi0 + 95);
    const auto *isi0_96 = buffer.data(isi0 + 96);
    const auto *isi0_98 = buffer.data(isi0 + 98);
    const auto *isi0_99 = buffer.data(isi0 + 99);
    const auto *isi0_105 = buffer.data(isi0 + 105);
    const auto *isi0_106 = buffer.data(isi0 + 106);
    const auto *isi0_107 = buffer.data(isi0 + 107);
    const auto *isi0_108 = buffer.data(isi0 + 108);
    const auto *isi0_109 = buffer.data(isi0 + 109);
    const auto *isi0_111 = buffer.data(isi0 + 111);
    const auto *isi0_135 = buffer.data(isi0 + 135);
    const auto *isi0_136 = buffer.data(isi0 + 136);
    const auto *isi0_137 = buffer.data(isi0 + 137);
    const auto *isi0_138 = buffer.data(isi0 + 138);
    const auto *isi0_139 = buffer.data(isi0 + 139);
    const auto *isi0_140 = buffer.data(isi0 + 140);
    const auto *isi0_141 = buffer.data(isi0 + 141);
    const auto *isi0_142 = buffer.data(isi0 + 142);
    const auto *isi0_143 = buffer.data(isi0 + 143);
    const auto *isi0_144 = buffer.data(isi0 + 144);
    const auto *isi0_145 = buffer.data(isi0 + 145);
    const auto *isi0_146 = buffer.data(isi0 + 146);
    const auto *isi0_147 = buffer.data(isi0 + 147);
    const auto *isi0_148 = buffer.data(isi0 + 148);
    const auto *isi0_149 = buffer.data(isi0 + 149);
    const auto *isi0_154 = buffer.data(isi0 + 154);
    const auto *isi0_160 = buffer.data(isi0 + 160);

    const auto *isi1_80 = buffer.data(isi1 + 80);
    const auto *isi1_81 = buffer.data(isi1 + 81);
    const auto *isi1_82 = buffer.data(isi1 + 82);
    const auto *isi1_83 = buffer.data(isi1 + 83);
    const auto *isi1_84 = buffer.data(isi1 + 84);
    const auto *isi1_86 = buffer.data(isi1 + 86);
    const auto *isi1_87 = buffer.data(isi1 + 87);
    const auto *isi1_89 = buffer.data(isi1 + 89);
    const auto *isi1_90 = buffer.data(isi1 + 90);
    const auto *isi1_91 = buffer.data(isi1 + 91);
    const auto *isi1_93 = buffer.data(isi1 + 93);
    const auto *isi1_94 = buffer.data(isi1 + 94);
    const auto *isi1_95 = buffer.data(isi1 + 95);
    const auto *isi1_96 = buffer.data(isi1 + 96);
    const auto *isi1_98 = buffer.data(isi1 + 98);
    const auto *isi1_99 = buffer.data(isi1 + 99);
    const auto *isi1_105 = buffer.data(isi1 + 105);
    const auto *isi1_106 = buffer.data(isi1 + 106);
    const auto *isi1_107 = buffer.data(isi1 + 107);
    const auto *isi1_108 = buffer.data(isi1 + 108);
    const auto *isi1_109 = buffer.data(isi1 + 109);
    const auto *isi1_111 = buffer.data(isi1 + 111);
    const auto *isi1_135 = buffer.data(isi1 + 135);
    const auto *isi1_136 = buffer.data(isi1 + 136);
    const auto *isi1_137 = buffer.data(isi1 + 137);
    const auto *isi1_138 = buffer.data(isi1 + 138);
    const auto *isi1_139 = buffer.data(isi1 + 139);
    const auto *isi1_140 = buffer.data(isi1 + 140);
    const auto *isi1_141 = buffer.data(isi1 + 141);
    const auto *isi1_142 = buffer.data(isi1 + 142);
    const auto *isi1_143 = buffer.data(isi1 + 143);
    const auto *isi1_144 = buffer.data(isi1 + 144);
    const auto *isi1_145 = buffer.data(isi1 + 145);
    const auto *isi1_146 = buffer.data(isi1 + 146);
    const auto *isi1_147 = buffer.data(isi1 + 147);
    const auto *isi1_148 = buffer.data(isi1 + 148);
    const auto *isi1_149 = buffer.data(isi1 + 149);
    const auto *isi1_154 = buffer.data(isi1 + 154);
    const auto *isi1_160 = buffer.data(isi1 + 160);

    const auto *isk_103 = buffer.data(isk + 103);
    const auto *isk_104 = buffer.data(isk + 104);
    const auto *isk_105 = buffer.data(isk + 105);
    const auto *isk_106 = buffer.data(isk + 106);
    const auto *isk_107 = buffer.data(isk + 107);
    const auto *isk_108 = buffer.data(isk + 108);
    const auto *isk_109 = buffer.data(isk + 109);
    const auto *isk_110 = buffer.data(isk + 110);
    const auto *isk_111 = buffer.data(isk + 111);
    const auto *isk_113 = buffer.data(isk + 113);
    const auto *isk_114 = buffer.data(isk + 114);
    const auto *isk_115 = buffer.data(isk + 115);
    const auto *isk_117 = buffer.data(isk + 117);
    const auto *isk_118 = buffer.data(isk + 118);
    const auto *isk_119 = buffer.data(isk + 119);
    const auto *isk_120 = buffer.data(isk + 120);
    const auto *isk_122 = buffer.data(isk + 122);
    const auto *isk_123 = buffer.data(isk + 123);
    const auto *isk_124 = buffer.data(isk + 124);
    const auto *isk_125 = buffer.data(isk + 125);
    const auto *isk_126 = buffer.data(isk + 126);
    const auto *isk_128 = buffer.data(isk + 128);
    const auto *isk_129 = buffer.data(isk + 129);
    const auto *isk_136 = buffer.data(isk + 136);
    const auto *isk_137 = buffer.data(isk + 137);
    const auto *isk_138 = buffer.data(isk + 138);
    const auto *isk_139 = buffer.data(isk + 139);
    const auto *isk_140 = buffer.data(isk + 140);
    const auto *isk_141 = buffer.data(isk + 141);
    const auto *isk_142 = buffer.data(isk + 142);
    const auto *isk_143 = buffer.data(isk + 143);
    const auto *isk_144 = buffer.data(isk + 144);
    const auto *isk_146 = buffer.data(isk + 146);
    const auto *isk_147 = buffer.data(isk + 147);
    const auto *isk_149 = buffer.data(isk + 149);
    const auto *isk_150 = buffer.data(isk + 150);
    const auto *isk_153 = buffer.data(isk + 153);
    const auto *isk_154 = buffer.data(isk + 154);
    const auto *isk_158 = buffer.data(isk + 158);
    const auto *isk_159 = buffer.data(isk + 159);
    const auto *isk_164 = buffer.data(isk + 164);
    const auto *isk_172 = buffer.data(isk + 172);
    const auto *isk_173 = buffer.data(isk + 173);
    const auto *isk_174 = buffer.data(isk + 174);
    const auto *isk_175 = buffer.data(isk + 175);
    const auto *isk_176 = buffer.data(isk + 176);
    const auto *isk_177 = buffer.data(isk + 177);
    const auto *isk_178 = buffer.data(isk + 178);
    const auto *isk_179 = buffer.data(isk + 179);
    const auto *isk_180 = buffer.data(isk + 180);
    const auto *isk_181 = buffer.data(isk + 181);
    const auto *isk_182 = buffer.data(isk + 182);
    const auto *isk_183 = buffer.data(isk + 183);
    const auto *isk_184 = buffer.data(isk + 184);
    const auto *isk_185 = buffer.data(isk + 185);
    const auto *isk_186 = buffer.data(isk + 186);
    const auto *isk_187 = buffer.data(isk + 187);
    const auto *isk_188 = buffer.data(isk + 188);
    const auto *isk_189 = buffer.data(isk + 189);
    const auto *isk_190 = buffer.data(isk + 190);
    const auto *isk_191 = buffer.data(isk + 191);
    const auto *isk_192 = buffer.data(isk + 192);
    const auto *isk_193 = buffer.data(isk + 193);
    const auto *isk_194 = buffer.data(isk + 194);
    const auto *isk_200 = buffer.data(isk + 200);

#pragma omp simd aligned(t_129, t_130, t_131, pc_y, isi0_80, isi0_81, isi0_82, isi1_80, \
                         isi1_81, isi1_82, isk_103, isk_104, isk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * isi0_80[k]
                   - f_11 * isi1_80[k]
                   + f_3 * pc_y[k] * isk_103[k];

        t_130[k] = f_8 * isi0_81[k]
                   - f_9 * isi1_81[k]
                   + f_3 * pc_y[k] * isk_104[k];

        t_131[k] = f_6 * isi0_82[k]
                   - f_7 * isi1_82[k]
                   + f_3 * pc_y[k] * isk_105[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, pc_y, pc_z, hsk_35, hsk_108, \
                         isi0_83, isi0_84, isi1_83, isi1_84, isk_106, isk_107, \
                         isk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * isi0_83[k]
                   - f_5 * isi1_83[k]
                   + f_3 * pc_y[k] * isk_106[k];

        t_133[k] = f_3 * pc_y[k] * isk_107[k];

        t_134[k] = f_15 * hsk_35[k]
                   + f_1 * isi0_83[k]
                   - f_2 * isi1_83[k]
                   + f_3 * pc_z[k] * isk_107[k];

        t_135[k] = f_18 * hsk_108[k]
                   + f_1 * isi0_84[k]
                   - f_2 * isi1_84[k]
                   + f_3 * pc_x[k] * isk_108[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, hsk_36, hsk_111, \
                         isi0_87, isi1_87, isk_108, isk_109, isk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_16 * hsk_36[k]
                   + f_3 * pc_y[k] * isk_108[k];

        t_137[k] = f_3 * pc_z[k] * isk_108[k];

        t_138[k] = f_18 * hsk_111[k]
                   + f_12 * isi0_87[k]
                   - f_13 * isi1_87[k]
                   + f_3 * pc_x[k] * isk_111[k];

        t_139[k] = f_3 * pc_z[k] * isk_109[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pc_x, pc_z, hsk_114, isi0_84, isi0_90, isi1_84, \
                         isi1_90, isk_110, isk_111, isk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * isi0_84[k]
                   - f_5 * isi1_84[k]
                   + f_3 * pc_z[k] * isk_110[k];

        t_141[k] = f_18 * hsk_114[k]
                   + f_10 * isi0_90[k]
                   - f_11 * isi1_90[k]
                   + f_3 * pc_x[k] * isk_114[k];

        t_142[k] = f_3 * pc_z[k] * isk_111[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, hsk_41, hsk_118, \
                         isi0_86, isi0_94, isi1_86, isi1_94, isk_113, isk_114, \
                         isk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * hsk_41[k]
                   + f_3 * pc_y[k] * isk_113[k];

        t_144[k] = f_6 * isi0_86[k]
                   - f_7 * isi1_86[k]
                   + f_3 * pc_z[k] * isk_113[k];

        t_145[k] = f_18 * hsk_118[k]
                   + f_8 * isi0_94[k]
                   - f_9 * isi1_94[k]
                   + f_3 * pc_x[k] * isk_118[k];

        t_146[k] = f_3 * pc_z[k] * isk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_y, pc_z, hsk_45, isi0_87, isi0_89, isi1_87, \
                         isi1_89, isk_115, isk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * isi0_87[k]
                   - f_5 * isi1_87[k]
                   + f_3 * pc_z[k] * isk_115[k];

        t_148[k] = f_16 * hsk_45[k]
                   + f_3 * pc_y[k] * isk_117[k];

        t_149[k] = f_8 * isi0_89[k]
                   - f_9 * isi1_89[k]
                   + f_3 * pc_z[k] * isk_117[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, hsk_123, isi0_90, isi0_99, isi1_90, \
                         isi1_99, isk_118, isk_119, isk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_18 * hsk_123[k]
                   + f_6 * isi0_99[k]
                   - f_7 * isi1_99[k]
                   + f_3 * pc_x[k] * isk_123[k];

        t_151[k] = f_3 * pc_z[k] * isk_118[k];

        t_152[k] = f_4 * isi0_90[k]
                   - f_5 * isi1_90[k]
                   + f_3 * pc_z[k] * isk_119[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_y, pc_z, hsk_50, isi0_91, isi0_93, isi1_91, \
                         isi1_93, isk_120, isk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * isi0_91[k]
                   - f_7 * isi1_91[k]
                   + f_3 * pc_z[k] * isk_120[k];

        t_154[k] = f_16 * hsk_50[k]
                   + f_3 * pc_y[k] * isk_122[k];

        t_155[k] = f_10 * isi0_93[k]
                   - f_11 * isi1_93[k]
                   + f_3 * pc_z[k] * isk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, hsk_129, isi0_94, isi0_105, isi1_94, \
                         isi1_105, isk_123, isk_124, isk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_18 * hsk_129[k]
                   + f_4 * isi0_105[k]
                   - f_5 * isi1_105[k]
                   + f_3 * pc_x[k] * isk_129[k];

        t_157[k] = f_3 * pc_z[k] * isk_123[k];

        t_158[k] = f_4 * isi0_94[k]
                   - f_5 * isi1_94[k]
                   + f_3 * pc_z[k] * isk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pc_y, pc_z, hsk_56, isi0_95, isi0_96, \
                         isi0_98, isi1_95, isi1_96, isi1_98, isk_125, isk_126, \
                         isk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_6 * isi0_95[k]
                   - f_7 * isi1_95[k]
                   + f_3 * pc_z[k] * isk_125[k];

        t_160[k] = f_8 * isi0_96[k]
                   - f_9 * isi1_96[k]
                   + f_3 * pc_z[k] * isk_126[k];

        t_161[k] = f_16 * hsk_56[k]
                   + f_3 * pc_y[k] * isk_128[k];

        t_162[k] = f_12 * isi0_98[k]
                   - f_13 * isi1_98[k]
                   + f_3 * pc_z[k] * isk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pc_x, pc_z, hsk_136, hsk_138, \
                         hsk_139, hsk_140, isk_129, isk_136, isk_138, isk_139, \
                         isk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_18 * hsk_136[k]
                   + f_3 * pc_x[k] * isk_136[k];

        t_164[k] = f_3 * pc_z[k] * isk_129[k];

        t_165[k] = f_18 * hsk_138[k]
                   + f_3 * pc_x[k] * isk_138[k];

        t_166[k] = f_18 * hsk_139[k]
                   + f_3 * pc_x[k] * isk_139[k];

        t_167[k] = f_18 * hsk_140[k]
                   + f_3 * pc_x[k] * isk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, hsk_64, hsk_141, hsk_142, \
                         hsk_143, isi0_105, isi1_105, isk_136, isk_141, isk_142, \
                         isk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_18 * hsk_141[k]
                   + f_3 * pc_x[k] * isk_141[k];

        t_169[k] = f_18 * hsk_142[k]
                   + f_3 * pc_x[k] * isk_142[k];

        t_170[k] = f_18 * hsk_143[k]
                   + f_3 * pc_x[k] * isk_143[k];

        t_171[k] = f_16 * hsk_64[k]
                   + f_1 * isi0_105[k]
                   - f_2 * isi1_105[k]
                   + f_3 * pc_y[k] * isk_136[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_z, isi0_105, isi0_106, isi0_107, \
                         isi1_105, isi1_106, isi1_107, isk_136, isk_137, isk_138, \
                         isk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * isk_136[k];

        t_173[k] = f_4 * isi0_105[k]
                   - f_5 * isi1_105[k]
                   + f_3 * pc_z[k] * isk_137[k];

        t_174[k] = f_6 * isi0_106[k]
                   - f_7 * isi1_106[k]
                   + f_3 * pc_z[k] * isk_138[k];

        t_175[k] = f_8 * isi0_107[k]
                   - f_9 * isi1_107[k]
                   + f_3 * pc_z[k] * isk_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_y, pc_z, hsk_71, isi0_108, isi0_109, \
                         isi0_111, isi1_108, isi1_109, isi1_111, isk_140, isk_141, \
                         isk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * isi0_108[k]
                   - f_11 * isi1_108[k]
                   + f_3 * pc_z[k] * isk_140[k];

        t_177[k] = f_12 * isi0_109[k]
                   - f_13 * isi1_109[k]
                   + f_3 * pc_z[k] * isk_141[k];

        t_178[k] = f_16 * hsk_71[k]
                   + f_3 * pc_y[k] * isk_143[k];

        t_179[k] = f_1 * isi0_111[k]
                   - f_2 * isi1_111[k]
                   + f_3 * pc_z[k] * isk_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, hsl0_48, hsl0_90, \
                         hsk_36, hsk_72, hsl1_48, hsl1_90, isk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * hsl0_90[k]
                   - f_14 * pc_y[k] * hsl1_90[k];

        t_181[k] = f_15 * hsk_72[k]
                   + f_3 * pc_y[k] * isk_144[k];

        t_182[k] = f_15 * hsk_36[k]
                   + f_3 * pc_z[k] * isk_144[k];

        t_183[k] = pa_z[k] * hsl0_48[k]
                   - f_14 * pc_z[k] * hsl1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pa_z, pc_y, pc_z, hsl0_51, hsl0_95, \
                         hsk_39, hsk_74, hsl1_51, hsl1_95, isk_146, \
                         isk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_15 * hsk_74[k]
                   + f_3 * pc_y[k] * isk_146[k];

        t_185[k] = pa_y[k] * hsl0_95[k]
                   - f_14 * pc_y[k] * hsl1_95[k];

        t_186[k] = pa_z[k] * hsl0_51[k]
                   - f_14 * pc_z[k] * hsl1_51[k];

        t_187[k] = f_15 * hsk_39[k]
                   + f_3 * pc_z[k] * isk_147[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pc_y, pc_z, hsl0_55, hsl0_99, \
                         hsk_42, hsk_77, hsl1_55, hsl1_99, isk_149, \
                         isk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_15 * hsk_77[k]
                   + f_3 * pc_y[k] * isk_149[k];

        t_189[k] = pa_y[k] * hsl0_99[k]
                   - f_14 * pc_y[k] * hsl1_99[k];

        t_190[k] = pa_z[k] * hsl0_55[k]
                   - f_14 * pc_z[k] * hsl1_55[k];

        t_191[k] = f_15 * hsk_42[k]
                   + f_3 * pc_z[k] * isk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pc_y, hsl0_102, hsl0_104, hsk_80, hsk_81, \
                         hsl1_102, hsl1_104, isk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_y[k] * hsl0_102[k]
                   + f_16 * hsk_80[k]
                   - f_14 * pc_y[k] * hsl1_102[k];

        t_193[k] = f_15 * hsk_81[k]
                   + f_3 * pc_y[k] * isk_153[k];

        t_194[k] = pa_y[k] * hsl0_104[k]
                   - f_14 * pc_y[k] * hsl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_y, pa_z, pc_y, pc_z, hsl0_60, hsl0_107, \
                         hsk_46, hsk_84, hsl1_60, hsl1_107, isk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * hsl0_60[k]
                   - f_14 * pc_z[k] * hsl1_60[k];

        t_196[k] = f_15 * hsk_46[k]
                   + f_3 * pc_z[k] * isk_154[k];

        t_197[k] = pa_y[k] * hsl0_107[k]
                   + f_17 * hsk_84[k]
                   - f_14 * pc_y[k] * hsl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_y, pc_y, hsl0_108, hsl0_110, hsk_85, hsk_86, \
                         hsl1_108, hsl1_110, isk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_y[k] * hsl0_108[k]
                   + f_16 * hsk_85[k]
                   - f_14 * pc_y[k] * hsl1_108[k];

        t_199[k] = f_15 * hsk_86[k]
                   + f_3 * pc_y[k] * isk_158[k];

        t_200[k] = pa_y[k] * hsl0_110[k]
                   - f_14 * pc_y[k] * hsl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_y, pa_z, pc_y, pc_z, hsl0_66, hsl0_113, \
                         hsk_51, hsk_89, hsl1_66, hsl1_113, isk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_z[k] * hsl0_66[k]
                   - f_14 * pc_z[k] * hsl1_66[k];

        t_202[k] = f_15 * hsk_51[k]
                   + f_3 * pc_z[k] * isk_159[k];

        t_203[k] = pa_y[k] * hsl0_113[k]
                   + f_18 * hsk_89[k]
                   - f_14 * pc_y[k] * hsl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_y, hsl0_114, hsl0_115, hsl0_117, \
                         hsk_90, hsk_91, hsk_92, hsl1_114, hsl1_115, hsl1_117, \
                         isk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * hsl0_114[k]
                   + f_17 * hsk_90[k]
                   - f_14 * pc_y[k] * hsl1_114[k];

        t_205[k] = pa_y[k] * hsl0_115[k]
                   + f_16 * hsk_91[k]
                   - f_14 * pc_y[k] * hsl1_115[k];

        t_206[k] = f_15 * hsk_92[k]
                   + f_3 * pc_y[k] * isk_164[k];

        t_207[k] = pa_y[k] * hsl0_117[k]
                   - f_14 * pc_y[k] * hsl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, hsk_172, hsk_173, hsk_174, \
                         hsk_175, hsk_176, isk_172, isk_173, isk_174, isk_175, \
                         isk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_18 * hsk_172[k]
                   + f_3 * pc_x[k] * isk_172[k];

        t_209[k] = f_18 * hsk_173[k]
                   + f_3 * pc_x[k] * isk_173[k];

        t_210[k] = f_18 * hsk_174[k]
                   + f_3 * pc_x[k] * isk_174[k];

        t_211[k] = f_18 * hsk_175[k]
                   + f_3 * pc_x[k] * isk_175[k];

        t_212[k] = f_18 * hsk_176[k]
                   + f_3 * pc_x[k] * isk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_z, pc_x, pc_z, hsl0_81, hsk_177, \
                         hsk_178, hsk_179, hsl1_81, isk_177, isk_178, \
                         isk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_18 * hsk_177[k]
                   + f_3 * pc_x[k] * isk_177[k];

        t_214[k] = f_18 * hsk_178[k]
                   + f_3 * pc_x[k] * isk_178[k];

        t_215[k] = f_18 * hsk_179[k]
                   + f_3 * pc_x[k] * isk_179[k];

        t_216[k] = pa_z[k] * hsl0_81[k]
                   - f_14 * pc_z[k] * hsl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, hsk_64, hsk_102, hsk_103, isi0_135, \
                         isi0_136, isi1_135, isi1_136, isk_172, isk_174, \
                         isk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * hsk_64[k]
                   + f_3 * pc_z[k] * isk_172[k];

        t_218[k] = f_15 * hsk_102[k]
                   + f_12 * isi0_135[k]
                   - f_13 * isi1_135[k]
                   + f_3 * pc_y[k] * isk_174[k];

        t_219[k] = f_15 * hsk_103[k]
                   + f_10 * isi0_136[k]
                   - f_11 * isi1_136[k]
                   + f_3 * pc_y[k] * isk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, hsk_104, hsk_105, hsk_106, isi0_137, \
                         isi0_138, isi0_139, isi1_137, isi1_138, isi1_139, isk_176, isk_177, \
                         isk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * hsk_104[k]
                   + f_8 * isi0_137[k]
                   - f_9 * isi1_137[k]
                   + f_3 * pc_y[k] * isk_176[k];

        t_221[k] = f_15 * hsk_105[k]
                   + f_6 * isi0_138[k]
                   - f_7 * isi1_138[k]
                   + f_3 * pc_y[k] * isk_177[k];

        t_222[k] = f_15 * hsk_106[k]
                   + f_4 * isi0_139[k]
                   - f_5 * isi1_139[k]
                   + f_3 * pc_y[k] * isk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_y, pc_x, pc_y, hsl0_134, hsk_107, \
                         hsk_180, hsl1_134, isi0_140, isi1_140, isk_179, \
                         isk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * hsk_107[k]
                   + f_3 * pc_y[k] * isk_179[k];

        t_224[k] = pa_y[k] * hsl0_134[k]
                   - f_14 * pc_y[k] * hsl1_134[k];

        t_225[k] = f_18 * hsk_180[k]
                   + f_1 * isi0_140[k]
                   - f_2 * isi1_140[k]
                   + f_3 * pc_x[k] * isk_180[k];

        t_226[k] = f_3 * pc_y[k] * isk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_y, pc_z, hsk_72, isi0_140, isi1_140, isk_180, \
                         isk_181, isk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * hsk_72[k]
                   + f_3 * pc_z[k] * isk_180[k];

        t_228[k] = f_4 * isi0_140[k]
                   - f_5 * isi1_140[k]
                   + f_3 * pc_y[k] * isk_181[k];

        t_229[k] = f_3 * pc_y[k] * isk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, hsk_185, isi0_141, isi0_142, \
                         isi0_145, isi1_141, isi1_142, isi1_145, isk_183, isk_184, \
                         isk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_18 * hsk_185[k]
                   + f_12 * isi0_145[k]
                   - f_13 * isi1_145[k]
                   + f_3 * pc_x[k] * isk_185[k];

        t_231[k] = f_6 * isi0_141[k]
                   - f_7 * isi1_141[k]
                   + f_3 * pc_y[k] * isk_183[k];

        t_232[k] = f_4 * isi0_142[k]
                   - f_5 * isi1_142[k]
                   + f_3 * pc_y[k] * isk_184[k];

        t_233[k] = f_3 * pc_y[k] * isk_185[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, hsk_189, isi0_143, isi0_144, \
                         isi0_149, isi1_143, isi1_144, isi1_149, isk_186, isk_187, \
                         isk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_18 * hsk_189[k]
                   + f_10 * isi0_149[k]
                   - f_11 * isi1_149[k]
                   + f_3 * pc_x[k] * isk_189[k];

        t_235[k] = f_8 * isi0_143[k]
                   - f_9 * isi1_143[k]
                   + f_3 * pc_y[k] * isk_186[k];

        t_236[k] = f_6 * isi0_144[k]
                   - f_7 * isi1_144[k]
                   + f_3 * pc_y[k] * isk_187[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_x, pc_y, hsk_194, isi0_145, isi0_154, \
                         isi1_145, isi1_154, isk_188, isk_189, \
                         isk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * isi0_145[k]
                   - f_5 * isi1_145[k]
                   + f_3 * pc_y[k] * isk_188[k];

        t_238[k] = f_3 * pc_y[k] * isk_189[k];

        t_239[k] = f_18 * hsk_194[k]
                   + f_8 * isi0_154[k]
                   - f_9 * isi1_154[k]
                   + f_3 * pc_x[k] * isk_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pc_y, isi0_146, isi0_147, isi0_148, isi1_146, \
                         isi1_147, isi1_148, isk_190, isk_191, \
                         isk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * isi0_146[k]
                   - f_11 * isi1_146[k]
                   + f_3 * pc_y[k] * isk_190[k];

        t_241[k] = f_8 * isi0_147[k]
                   - f_9 * isi1_147[k]
                   + f_3 * pc_y[k] * isk_191[k];

        t_242[k] = f_6 * isi0_148[k]
                   - f_7 * isi1_148[k]
                   + f_3 * pc_y[k] * isk_192[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pc_x, pc_y, hsk_200, isi0_149, isi0_160, \
                         isi1_149, isi1_160, isk_193, isk_194, \
                         isk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_4 * isi0_149[k]
                   - f_5 * isi1_149[k]
                   + f_3 * pc_y[k] * isk_193[k];

        t_244[k] = f_3 * pc_y[k] * isk_194[k];

        t_245[k] = f_18 * hsk_200[k]
                   + f_6 * isi0_160[k]
                   - f_7 * isi1_160[k]
                   + f_3 * pc_x[k] * isk_200[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
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

    const auto *hsl0_135 = buffer.data(hsl0 + 135);
    const auto *hsl0_138 = buffer.data(hsl0 + 138);
    const auto *hsl0_141 = buffer.data(hsl0 + 141);
    const auto *hsl0_145 = buffer.data(hsl0 + 145);
    const auto *hsl0_147 = buffer.data(hsl0 + 147);
    const auto *hsl0_150 = buffer.data(hsl0 + 150);
    const auto *hsl0_152 = buffer.data(hsl0 + 152);
    const auto *hsl0_153 = buffer.data(hsl0 + 153);
    const auto *hsl0_156 = buffer.data(hsl0 + 156);
    const auto *hsl0_158 = buffer.data(hsl0 + 158);
    const auto *hsl0_159 = buffer.data(hsl0 + 159);
    const auto *hsl0_160 = buffer.data(hsl0 + 160);
    const auto *hsl0_171 = buffer.data(hsl0 + 171);
    const auto *hsl0_225 = buffer.data(hsl0 + 225);

    const auto *hsk_107 = buffer.data(hsk + 107);
    const auto *hsk_108 = buffer.data(hsk + 108);
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
    const auto *hsk_136 = buffer.data(hsk + 136);
    const auto *hsk_143 = buffer.data(hsk + 143);
    const auto *hsk_144 = buffer.data(hsk + 144);
    const auto *hsk_146 = buffer.data(hsk + 146);
    const auto *hsk_149 = buffer.data(hsk + 149);
    const auto *hsk_153 = buffer.data(hsk + 153);
    const auto *hsk_158 = buffer.data(hsk + 158);
    const auto *hsk_164 = buffer.data(hsk + 164);
    const auto *hsk_174 = buffer.data(hsk + 174);
    const auto *hsk_175 = buffer.data(hsk + 175);
    const auto *hsk_176 = buffer.data(hsk + 176);
    const auto *hsk_177 = buffer.data(hsk + 177);
    const auto *hsk_178 = buffer.data(hsk + 178);
    const auto *hsk_179 = buffer.data(hsk + 179);
    const auto *hsk_180 = buffer.data(hsk + 180);
    const auto *hsk_207 = buffer.data(hsk + 207);
    const auto *hsk_208 = buffer.data(hsk + 208);
    const auto *hsk_209 = buffer.data(hsk + 209);
    const auto *hsk_210 = buffer.data(hsk + 210);
    const auto *hsk_211 = buffer.data(hsk + 211);
    const auto *hsk_212 = buffer.data(hsk + 212);
    const auto *hsk_213 = buffer.data(hsk + 213);
    const auto *hsk_215 = buffer.data(hsk + 215);
    const auto *hsk_216 = buffer.data(hsk + 216);
    const auto *hsk_219 = buffer.data(hsk + 219);
    const auto *hsk_222 = buffer.data(hsk + 222);
    const auto *hsk_226 = buffer.data(hsk + 226);
    const auto *hsk_231 = buffer.data(hsk + 231);
    const auto *hsk_237 = buffer.data(hsk + 237);
    const auto *hsk_244 = buffer.data(hsk + 244);
    const auto *hsk_246 = buffer.data(hsk + 246);
    const auto *hsk_247 = buffer.data(hsk + 247);
    const auto *hsk_248 = buffer.data(hsk + 248);
    const auto *hsk_249 = buffer.data(hsk + 249);
    const auto *hsk_250 = buffer.data(hsk + 250);
    const auto *hsk_251 = buffer.data(hsk + 251);
    const auto *hsk_257 = buffer.data(hsk + 257);
    const auto *hsk_261 = buffer.data(hsk + 261);
    const auto *hsk_266 = buffer.data(hsk + 266);
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

    const auto *hsl1_135 = buffer.data(hsl1 + 135);
    const auto *hsl1_138 = buffer.data(hsl1 + 138);
    const auto *hsl1_141 = buffer.data(hsl1 + 141);
    const auto *hsl1_145 = buffer.data(hsl1 + 145);
    const auto *hsl1_147 = buffer.data(hsl1 + 147);
    const auto *hsl1_150 = buffer.data(hsl1 + 150);
    const auto *hsl1_152 = buffer.data(hsl1 + 152);
    const auto *hsl1_153 = buffer.data(hsl1 + 153);
    const auto *hsl1_156 = buffer.data(hsl1 + 156);
    const auto *hsl1_158 = buffer.data(hsl1 + 158);
    const auto *hsl1_159 = buffer.data(hsl1 + 159);
    const auto *hsl1_160 = buffer.data(hsl1 + 160);
    const auto *hsl1_171 = buffer.data(hsl1 + 171);
    const auto *hsl1_225 = buffer.data(hsl1 + 225);

    const auto *isi0_150 = buffer.data(isi0 + 150);
    const auto *isi0_151 = buffer.data(isi0 + 151);
    const auto *isi0_152 = buffer.data(isi0 + 152);
    const auto *isi0_153 = buffer.data(isi0 + 153);
    const auto *isi0_154 = buffer.data(isi0 + 154);
    const auto *isi0_161 = buffer.data(isi0 + 161);
    const auto *isi0_162 = buffer.data(isi0 + 162);
    const auto *isi0_163 = buffer.data(isi0 + 163);
    const auto *isi0_164 = buffer.data(isi0 + 164);
    const auto *isi0_165 = buffer.data(isi0 + 165);
    const auto *isi0_166 = buffer.data(isi0 + 166);
    const auto *isi0_167 = buffer.data(isi0 + 167);
    const auto *isi0_168 = buffer.data(isi0 + 168);
    const auto *isi0_170 = buffer.data(isi0 + 170);
    const auto *isi0_171 = buffer.data(isi0 + 171);
    const auto *isi0_173 = buffer.data(isi0 + 173);
    const auto *isi0_174 = buffer.data(isi0 + 174);
    const auto *isi0_175 = buffer.data(isi0 + 175);
    const auto *isi0_177 = buffer.data(isi0 + 177);
    const auto *isi0_178 = buffer.data(isi0 + 178);
    const auto *isi0_179 = buffer.data(isi0 + 179);
    const auto *isi0_180 = buffer.data(isi0 + 180);
    const auto *isi0_182 = buffer.data(isi0 + 182);
    const auto *isi0_183 = buffer.data(isi0 + 183);
    const auto *isi0_189 = buffer.data(isi0 + 189);
    const auto *isi0_190 = buffer.data(isi0 + 190);
    const auto *isi0_191 = buffer.data(isi0 + 191);
    const auto *isi0_192 = buffer.data(isi0 + 192);
    const auto *isi0_193 = buffer.data(isi0 + 193);
    const auto *isi0_195 = buffer.data(isi0 + 195);
    const auto *isi0_201 = buffer.data(isi0 + 201);
    const auto *isi0_205 = buffer.data(isi0 + 205);
    const auto *isi0_210 = buffer.data(isi0 + 210);
    const auto *isi0_216 = buffer.data(isi0 + 216);
    const auto *isi0_219 = buffer.data(isi0 + 219);
    const auto *isi0_220 = buffer.data(isi0 + 220);
    const auto *isi0_221 = buffer.data(isi0 + 221);
    const auto *isi0_222 = buffer.data(isi0 + 222);
    const auto *isi0_223 = buffer.data(isi0 + 223);

    const auto *isi1_150 = buffer.data(isi1 + 150);
    const auto *isi1_151 = buffer.data(isi1 + 151);
    const auto *isi1_152 = buffer.data(isi1 + 152);
    const auto *isi1_153 = buffer.data(isi1 + 153);
    const auto *isi1_154 = buffer.data(isi1 + 154);
    const auto *isi1_161 = buffer.data(isi1 + 161);
    const auto *isi1_162 = buffer.data(isi1 + 162);
    const auto *isi1_163 = buffer.data(isi1 + 163);
    const auto *isi1_164 = buffer.data(isi1 + 164);
    const auto *isi1_165 = buffer.data(isi1 + 165);
    const auto *isi1_166 = buffer.data(isi1 + 166);
    const auto *isi1_167 = buffer.data(isi1 + 167);
    const auto *isi1_168 = buffer.data(isi1 + 168);
    const auto *isi1_170 = buffer.data(isi1 + 170);
    const auto *isi1_171 = buffer.data(isi1 + 171);
    const auto *isi1_173 = buffer.data(isi1 + 173);
    const auto *isi1_174 = buffer.data(isi1 + 174);
    const auto *isi1_175 = buffer.data(isi1 + 175);
    const auto *isi1_177 = buffer.data(isi1 + 177);
    const auto *isi1_178 = buffer.data(isi1 + 178);
    const auto *isi1_179 = buffer.data(isi1 + 179);
    const auto *isi1_180 = buffer.data(isi1 + 180);
    const auto *isi1_182 = buffer.data(isi1 + 182);
    const auto *isi1_183 = buffer.data(isi1 + 183);
    const auto *isi1_189 = buffer.data(isi1 + 189);
    const auto *isi1_190 = buffer.data(isi1 + 190);
    const auto *isi1_191 = buffer.data(isi1 + 191);
    const auto *isi1_192 = buffer.data(isi1 + 192);
    const auto *isi1_193 = buffer.data(isi1 + 193);
    const auto *isi1_195 = buffer.data(isi1 + 195);
    const auto *isi1_201 = buffer.data(isi1 + 201);
    const auto *isi1_205 = buffer.data(isi1 + 205);
    const auto *isi1_210 = buffer.data(isi1 + 210);
    const auto *isi1_216 = buffer.data(isi1 + 216);
    const auto *isi1_219 = buffer.data(isi1 + 219);
    const auto *isi1_220 = buffer.data(isi1 + 220);
    const auto *isi1_221 = buffer.data(isi1 + 221);
    const auto *isi1_222 = buffer.data(isi1 + 222);
    const auto *isi1_223 = buffer.data(isi1 + 223);

    const auto *isk_195 = buffer.data(isk + 195);
    const auto *isk_196 = buffer.data(isk + 196);
    const auto *isk_197 = buffer.data(isk + 197);
    const auto *isk_198 = buffer.data(isk + 198);
    const auto *isk_199 = buffer.data(isk + 199);
    const auto *isk_200 = buffer.data(isk + 200);
    const auto *isk_207 = buffer.data(isk + 207);
    const auto *isk_208 = buffer.data(isk + 208);
    const auto *isk_209 = buffer.data(isk + 209);
    const auto *isk_210 = buffer.data(isk + 210);
    const auto *isk_211 = buffer.data(isk + 211);
    const auto *isk_212 = buffer.data(isk + 212);
    const auto *isk_213 = buffer.data(isk + 213);
    const auto *isk_214 = buffer.data(isk + 214);
    const auto *isk_215 = buffer.data(isk + 215);
    const auto *isk_216 = buffer.data(isk + 216);
    const auto *isk_217 = buffer.data(isk + 217);
    const auto *isk_218 = buffer.data(isk + 218);
    const auto *isk_219 = buffer.data(isk + 219);
    const auto *isk_221 = buffer.data(isk + 221);
    const auto *isk_222 = buffer.data(isk + 222);
    const auto *isk_223 = buffer.data(isk + 223);
    const auto *isk_225 = buffer.data(isk + 225);
    const auto *isk_226 = buffer.data(isk + 226);
    const auto *isk_227 = buffer.data(isk + 227);
    const auto *isk_228 = buffer.data(isk + 228);
    const auto *isk_230 = buffer.data(isk + 230);
    const auto *isk_231 = buffer.data(isk + 231);
    const auto *isk_232 = buffer.data(isk + 232);
    const auto *isk_233 = buffer.data(isk + 233);
    const auto *isk_234 = buffer.data(isk + 234);
    const auto *isk_236 = buffer.data(isk + 236);
    const auto *isk_237 = buffer.data(isk + 237);
    const auto *isk_244 = buffer.data(isk + 244);
    const auto *isk_245 = buffer.data(isk + 245);
    const auto *isk_246 = buffer.data(isk + 246);
    const auto *isk_247 = buffer.data(isk + 247);
    const auto *isk_248 = buffer.data(isk + 248);
    const auto *isk_249 = buffer.data(isk + 249);
    const auto *isk_250 = buffer.data(isk + 250);
    const auto *isk_251 = buffer.data(isk + 251);
    const auto *isk_252 = buffer.data(isk + 252);
    const auto *isk_254 = buffer.data(isk + 254);
    const auto *isk_255 = buffer.data(isk + 255);
    const auto *isk_257 = buffer.data(isk + 257);
    const auto *isk_258 = buffer.data(isk + 258);
    const auto *isk_261 = buffer.data(isk + 261);
    const auto *isk_262 = buffer.data(isk + 262);
    const auto *isk_266 = buffer.data(isk + 266);
    const auto *isk_267 = buffer.data(isk + 267);
    const auto *isk_272 = buffer.data(isk + 272);
    const auto *isk_279 = buffer.data(isk + 279);
    const auto *isk_280 = buffer.data(isk + 280);
    const auto *isk_281 = buffer.data(isk + 281);
    const auto *isk_282 = buffer.data(isk + 282);
    const auto *isk_283 = buffer.data(isk + 283);
    const auto *isk_284 = buffer.data(isk + 284);
    const auto *isk_285 = buffer.data(isk + 285);
    const auto *isk_286 = buffer.data(isk + 286);
    const auto *isk_287 = buffer.data(isk + 287);
    const auto *isk_288 = buffer.data(isk + 288);

#pragma omp simd aligned(t_246, t_247, t_248, pc_y, isi0_150, isi0_151, isi0_152, isi1_150, \
                         isi1_151, isi1_152, isk_195, isk_196, \
                         isk_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * isi0_150[k]
                   - f_13 * isi1_150[k]
                   + f_3 * pc_y[k] * isk_195[k];

        t_247[k] = f_10 * isi0_151[k]
                   - f_11 * isi1_151[k]
                   + f_3 * pc_y[k] * isk_196[k];

        t_248[k] = f_8 * isi0_152[k]
                   - f_9 * isi1_152[k]
                   + f_3 * pc_y[k] * isk_197[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, isi0_153, isi0_154, isi1_153, isi1_154, \
                         isk_198, isk_199, isk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * isi0_153[k]
                   - f_7 * isi1_153[k]
                   + f_3 * pc_y[k] * isk_198[k];

        t_250[k] = f_4 * isi0_154[k]
                   - f_5 * isi1_154[k]
                   + f_3 * pc_y[k] * isk_199[k];

        t_251[k] = f_3 * pc_y[k] * isk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, hsk_207, hsk_208, hsk_209, hsk_210, \
                         isi0_167, isi1_167, isk_207, isk_208, isk_209, \
                         isk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_18 * hsk_207[k]
                   + f_4 * isi0_167[k]
                   - f_5 * isi1_167[k]
                   + f_3 * pc_x[k] * isk_207[k];

        t_253[k] = f_18 * hsk_208[k]
                   + f_3 * pc_x[k] * isk_208[k];

        t_254[k] = f_18 * hsk_209[k]
                   + f_3 * pc_x[k] * isk_209[k];

        t_255[k] = f_18 * hsk_210[k]
                   + f_3 * pc_x[k] * isk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pc_x, pc_y, hsk_211, hsk_212, \
                         hsk_213, hsk_215, isk_207, isk_211, isk_212, isk_213, \
                         isk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_18 * hsk_211[k]
                   + f_3 * pc_x[k] * isk_211[k];

        t_257[k] = f_18 * hsk_212[k]
                   + f_3 * pc_x[k] * isk_212[k];

        t_258[k] = f_18 * hsk_213[k]
                   + f_3 * pc_x[k] * isk_213[k];

        t_259[k] = f_3 * pc_y[k] * isk_207[k];

        t_260[k] = f_18 * hsk_215[k]
                   + f_3 * pc_x[k] * isk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_y, isi0_161, isi0_162, isi0_163, isi1_161, \
                         isi1_162, isi1_163, isk_208, isk_209, \
                         isk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * isi0_161[k]
                   - f_2 * isi1_161[k]
                   + f_3 * pc_y[k] * isk_208[k];

        t_262[k] = f_20 * isi0_162[k]
                   - f_21 * isi1_162[k]
                   + f_3 * pc_y[k] * isk_209[k];

        t_263[k] = f_12 * isi0_163[k]
                   - f_13 * isi1_163[k]
                   + f_3 * pc_y[k] * isk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, isi0_164, isi0_165, isi0_166, isi1_164, \
                         isi1_165, isi1_166, isk_211, isk_212, \
                         isk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * isi0_164[k]
                   - f_11 * isi1_164[k]
                   + f_3 * pc_y[k] * isk_211[k];

        t_265[k] = f_8 * isi0_165[k]
                   - f_9 * isi1_165[k]
                   + f_3 * pc_y[k] * isk_212[k];

        t_266[k] = f_6 * isi0_166[k]
                   - f_7 * isi1_166[k]
                   + f_3 * pc_y[k] * isk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, hsk_107, hsk_216, \
                         isi0_167, isi0_168, isi1_167, isi1_168, isk_214, isk_215, \
                         isk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_4 * isi0_167[k]
                   - f_5 * isi1_167[k]
                   + f_3 * pc_y[k] * isk_214[k];

        t_268[k] = f_3 * pc_y[k] * isk_215[k];

        t_269[k] = f_16 * hsk_107[k]
                   + f_1 * isi0_167[k]
                   - f_2 * isi1_167[k]
                   + f_3 * pc_z[k] * isk_215[k];

        t_270[k] = f_17 * hsk_216[k]
                   + f_1 * isi0_168[k]
                   - f_2 * isi1_168[k]
                   + f_3 * pc_x[k] * isk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, hsk_108, hsk_219, \
                         isi0_171, isi1_171, isk_216, isk_217, \
                         isk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * hsk_108[k]
                   + f_3 * pc_y[k] * isk_216[k];

        t_272[k] = f_3 * pc_z[k] * isk_216[k];

        t_273[k] = f_17 * hsk_219[k]
                   + f_12 * isi0_171[k]
                   - f_13 * isi1_171[k]
                   + f_3 * pc_x[k] * isk_219[k];

        t_274[k] = f_3 * pc_z[k] * isk_217[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, hsk_222, isi0_168, isi0_174, \
                         isi1_168, isi1_174, isk_218, isk_219, \
                         isk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_4 * isi0_168[k]
                   - f_5 * isi1_168[k]
                   + f_3 * pc_z[k] * isk_218[k];

        t_276[k] = f_17 * hsk_222[k]
                   + f_10 * isi0_174[k]
                   - f_11 * isi1_174[k]
                   + f_3 * pc_x[k] * isk_222[k];

        t_277[k] = f_3 * pc_z[k] * isk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pc_x, pc_y, pc_z, hsk_113, hsk_226, \
                         isi0_170, isi0_178, isi1_170, isi1_178, isk_221, isk_222, \
                         isk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * hsk_113[k]
                   + f_3 * pc_y[k] * isk_221[k];

        t_279[k] = f_6 * isi0_170[k]
                   - f_7 * isi1_170[k]
                   + f_3 * pc_z[k] * isk_221[k];

        t_280[k] = f_17 * hsk_226[k]
                   + f_8 * isi0_178[k]
                   - f_9 * isi1_178[k]
                   + f_3 * pc_x[k] * isk_226[k];

        t_281[k] = f_3 * pc_z[k] * isk_222[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pc_y, pc_z, hsk_117, isi0_171, isi0_173, \
                         isi1_171, isi1_173, isk_223, isk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_4 * isi0_171[k]
                   - f_5 * isi1_171[k]
                   + f_3 * pc_z[k] * isk_223[k];

        t_283[k] = f_17 * hsk_117[k]
                   + f_3 * pc_y[k] * isk_225[k];

        t_284[k] = f_8 * isi0_173[k]
                   - f_9 * isi1_173[k]
                   + f_3 * pc_z[k] * isk_225[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pc_x, pc_z, hsk_231, isi0_174, isi0_183, \
                         isi1_174, isi1_183, isk_226, isk_227, \
                         isk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_17 * hsk_231[k]
                   + f_6 * isi0_183[k]
                   - f_7 * isi1_183[k]
                   + f_3 * pc_x[k] * isk_231[k];

        t_286[k] = f_3 * pc_z[k] * isk_226[k];

        t_287[k] = f_4 * isi0_174[k]
                   - f_5 * isi1_174[k]
                   + f_3 * pc_z[k] * isk_227[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, hsk_122, isi0_175, isi0_177, \
                         isi1_175, isi1_177, isk_228, isk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_6 * isi0_175[k]
                   - f_7 * isi1_175[k]
                   + f_3 * pc_z[k] * isk_228[k];

        t_289[k] = f_17 * hsk_122[k]
                   + f_3 * pc_y[k] * isk_230[k];

        t_290[k] = f_10 * isi0_177[k]
                   - f_11 * isi1_177[k]
                   + f_3 * pc_z[k] * isk_230[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pc_x, pc_z, hsk_237, isi0_178, isi0_189, \
                         isi1_178, isi1_189, isk_231, isk_232, \
                         isk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_17 * hsk_237[k]
                   + f_4 * isi0_189[k]
                   - f_5 * isi1_189[k]
                   + f_3 * pc_x[k] * isk_237[k];

        t_292[k] = f_3 * pc_z[k] * isk_231[k];

        t_293[k] = f_4 * isi0_178[k]
                   - f_5 * isi1_178[k]
                   + f_3 * pc_z[k] * isk_232[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_y, pc_z, hsk_128, isi0_179, isi0_180, \
                         isi0_182, isi1_179, isi1_180, isi1_182, isk_233, isk_234, \
                         isk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * isi0_179[k]
                   - f_7 * isi1_179[k]
                   + f_3 * pc_z[k] * isk_233[k];

        t_295[k] = f_8 * isi0_180[k]
                   - f_9 * isi1_180[k]
                   + f_3 * pc_z[k] * isk_234[k];

        t_296[k] = f_17 * hsk_128[k]
                   + f_3 * pc_y[k] * isk_236[k];

        t_297[k] = f_12 * isi0_182[k]
                   - f_13 * isi1_182[k]
                   + f_3 * pc_z[k] * isk_236[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pc_x, pc_z, hsk_244, hsk_246, \
                         hsk_247, hsk_248, isk_237, isk_244, isk_246, isk_247, \
                         isk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_17 * hsk_244[k]
                   + f_3 * pc_x[k] * isk_244[k];

        t_299[k] = f_3 * pc_z[k] * isk_237[k];

        t_300[k] = f_17 * hsk_246[k]
                   + f_3 * pc_x[k] * isk_246[k];

        t_301[k] = f_17 * hsk_247[k]
                   + f_3 * pc_x[k] * isk_247[k];

        t_302[k] = f_17 * hsk_248[k]
                   + f_3 * pc_x[k] * isk_248[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, pc_y, hsk_136, hsk_249, hsk_250, \
                         hsk_251, isi0_189, isi1_189, isk_244, isk_249, isk_250, \
                         isk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_17 * hsk_249[k]
                   + f_3 * pc_x[k] * isk_249[k];

        t_304[k] = f_17 * hsk_250[k]
                   + f_3 * pc_x[k] * isk_250[k];

        t_305[k] = f_17 * hsk_251[k]
                   + f_3 * pc_x[k] * isk_251[k];

        t_306[k] = f_17 * hsk_136[k]
                   + f_1 * isi0_189[k]
                   - f_2 * isi1_189[k]
                   + f_3 * pc_y[k] * isk_244[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_z, isi0_189, isi0_190, isi0_191, \
                         isi1_189, isi1_190, isi1_191, isk_244, isk_245, isk_246, \
                         isk_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_z[k] * isk_244[k];

        t_308[k] = f_4 * isi0_189[k]
                   - f_5 * isi1_189[k]
                   + f_3 * pc_z[k] * isk_245[k];

        t_309[k] = f_6 * isi0_190[k]
                   - f_7 * isi1_190[k]
                   + f_3 * pc_z[k] * isk_246[k];

        t_310[k] = f_8 * isi0_191[k]
                   - f_9 * isi1_191[k]
                   + f_3 * pc_z[k] * isk_247[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, hsk_143, isi0_192, isi0_193, \
                         isi0_195, isi1_192, isi1_193, isi1_195, isk_248, isk_249, \
                         isk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_10 * isi0_192[k]
                   - f_11 * isi1_192[k]
                   + f_3 * pc_z[k] * isk_248[k];

        t_312[k] = f_12 * isi0_193[k]
                   - f_13 * isi1_193[k]
                   + f_3 * pc_z[k] * isk_249[k];

        t_313[k] = f_17 * hsk_143[k]
                   + f_3 * pc_y[k] * isk_251[k];

        t_314[k] = f_1 * isi0_195[k]
                   - f_2 * isi1_195[k]
                   + f_3 * pc_z[k] * isk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_z, pc_y, pc_z, hsl0_135, hsl0_138, \
                         hsk_108, hsk_144, hsl1_135, hsl1_138, \
                         isk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * hsl0_135[k]
                   - f_14 * pc_z[k] * hsl1_135[k];

        t_316[k] = f_16 * hsk_144[k]
                   + f_3 * pc_y[k] * isk_252[k];

        t_317[k] = f_15 * hsk_108[k]
                   + f_3 * pc_z[k] * isk_252[k];

        t_318[k] = pa_z[k] * hsl0_138[k]
                   - f_14 * pc_z[k] * hsl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_x, pc_y, pc_z, hsl0_141, hsk_146, \
                         hsk_257, hsl1_141, isi0_201, isi1_201, isk_254, \
                         isk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * hsk_146[k]
                   + f_3 * pc_y[k] * isk_254[k];

        t_320[k] = f_17 * hsk_257[k]
                   + f_12 * isi0_201[k]
                   - f_13 * isi1_201[k]
                   + f_3 * pc_x[k] * isk_257[k];

        t_321[k] = pa_z[k] * hsl0_141[k]
                   - f_14 * pc_z[k] * hsl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, hsk_111, hsk_149, hsk_261, \
                         isi0_205, isi1_205, isk_255, isk_257, \
                         isk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * hsk_111[k]
                   + f_3 * pc_z[k] * isk_255[k];

        t_323[k] = f_16 * hsk_149[k]
                   + f_3 * pc_y[k] * isk_257[k];

        t_324[k] = f_17 * hsk_261[k]
                   + f_10 * isi0_205[k]
                   - f_11 * isi1_205[k]
                   + f_3 * pc_x[k] * isk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pa_z, pc_y, pc_z, hsl0_145, hsl0_147, \
                         hsk_114, hsk_115, hsk_153, hsl1_145, hsl1_147, isk_258, \
                         isk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_z[k] * hsl0_145[k]
                   - f_14 * pc_z[k] * hsl1_145[k];

        t_326[k] = f_15 * hsk_114[k]
                   + f_3 * pc_z[k] * isk_258[k];

        t_327[k] = pa_z[k] * hsl0_147[k]
                   + f_16 * hsk_115[k]
                   - f_14 * pc_z[k] * hsl1_147[k];

        t_328[k] = f_16 * hsk_153[k]
                   + f_3 * pc_y[k] * isk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_z, pc_x, pc_z, hsl0_150, hsk_118, hsk_266, \
                         hsl1_150, isi0_210, isi1_210, isk_262, \
                         isk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_17 * hsk_266[k]
                   + f_8 * isi0_210[k]
                   - f_9 * isi1_210[k]
                   + f_3 * pc_x[k] * isk_266[k];

        t_330[k] = pa_z[k] * hsl0_150[k]
                   - f_14 * pc_z[k] * hsl1_150[k];

        t_331[k] = f_15 * hsk_118[k]
                   + f_3 * pc_z[k] * isk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pa_z, pc_y, pc_z, hsl0_152, hsl0_153, hsk_119, \
                         hsk_120, hsk_158, hsl1_152, hsl1_153, \
                         isk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pa_z[k] * hsl0_152[k]
                   + f_16 * hsk_119[k]
                   - f_14 * pc_z[k] * hsl1_152[k];

        t_333[k] = pa_z[k] * hsl0_153[k]
                   + f_17 * hsk_120[k]
                   - f_14 * pc_z[k] * hsl1_153[k];

        t_334[k] = f_16 * hsk_158[k]
                   + f_3 * pc_y[k] * isk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pa_z, pc_x, pc_z, hsl0_156, hsk_123, hsk_272, \
                         hsl1_156, isi0_216, isi1_216, isk_267, \
                         isk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_17 * hsk_272[k]
                   + f_6 * isi0_216[k]
                   - f_7 * isi1_216[k]
                   + f_3 * pc_x[k] * isk_272[k];

        t_336[k] = pa_z[k] * hsl0_156[k]
                   - f_14 * pc_z[k] * hsl1_156[k];

        t_337[k] = f_15 * hsk_123[k]
                   + f_3 * pc_z[k] * isk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pa_z, pc_z, hsl0_158, hsl0_159, hsl0_160, \
                         hsk_124, hsk_125, hsk_126, hsl1_158, hsl1_159, \
                         hsl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pa_z[k] * hsl0_158[k]
                   + f_16 * hsk_124[k]
                   - f_14 * pc_z[k] * hsl1_158[k];

        t_339[k] = pa_z[k] * hsl0_159[k]
                   + f_17 * hsk_125[k]
                   - f_14 * pc_z[k] * hsl1_159[k];

        t_340[k] = pa_z[k] * hsl0_160[k]
                   + f_18 * hsk_126[k]
                   - f_14 * pc_z[k] * hsl1_160[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, hsk_164, hsk_279, hsk_280, \
                         hsk_281, isi0_223, isi1_223, isk_272, isk_279, isk_280, \
                         isk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * hsk_164[k]
                   + f_3 * pc_y[k] * isk_272[k];

        t_342[k] = f_17 * hsk_279[k]
                   + f_4 * isi0_223[k]
                   - f_5 * isi1_223[k]
                   + f_3 * pc_x[k] * isk_279[k];

        t_343[k] = f_17 * hsk_280[k]
                   + f_3 * pc_x[k] * isk_280[k];

        t_344[k] = f_17 * hsk_281[k]
                   + f_3 * pc_x[k] * isk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, hsk_282, hsk_283, hsk_284, \
                         hsk_285, hsk_286, isk_282, isk_283, isk_284, isk_285, \
                         isk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_17 * hsk_282[k]
                   + f_3 * pc_x[k] * isk_282[k];

        t_346[k] = f_17 * hsk_283[k]
                   + f_3 * pc_x[k] * isk_283[k];

        t_347[k] = f_17 * hsk_284[k]
                   + f_3 * pc_x[k] * isk_284[k];

        t_348[k] = f_17 * hsk_285[k]
                   + f_3 * pc_x[k] * isk_285[k];

        t_349[k] = f_17 * hsk_286[k]
                   + f_3 * pc_x[k] * isk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pa_z, pc_x, pc_z, hsl0_171, hsk_136, hsk_287, \
                         hsl1_171, isk_280, isk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_17 * hsk_287[k]
                   + f_3 * pc_x[k] * isk_287[k];

        t_351[k] = pa_z[k] * hsl0_171[k]
                   - f_14 * pc_z[k] * hsl1_171[k];

        t_352[k] = f_15 * hsk_136[k]
                   + f_3 * pc_z[k] * isk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, hsk_174, hsk_175, hsk_176, isi0_219, \
                         isi0_220, isi0_221, isi1_219, isi1_220, isi1_221, isk_282, isk_283, \
                         isk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * hsk_174[k]
                   + f_12 * isi0_219[k]
                   - f_13 * isi1_219[k]
                   + f_3 * pc_y[k] * isk_282[k];

        t_354[k] = f_16 * hsk_175[k]
                   + f_10 * isi0_220[k]
                   - f_11 * isi1_220[k]
                   + f_3 * pc_y[k] * isk_283[k];

        t_355[k] = f_16 * hsk_176[k]
                   + f_8 * isi0_221[k]
                   - f_9 * isi1_221[k]
                   + f_3 * pc_y[k] * isk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, hsk_177, hsk_178, hsk_179, isi0_222, \
                         isi0_223, isi1_222, isi1_223, isk_285, isk_286, \
                         isk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * hsk_177[k]
                   + f_6 * isi0_222[k]
                   - f_7 * isi1_222[k]
                   + f_3 * pc_y[k] * isk_285[k];

        t_357[k] = f_16 * hsk_178[k]
                   + f_4 * isi0_223[k]
                   - f_5 * isi1_223[k]
                   + f_3 * pc_y[k] * isk_286[k];

        t_358[k] = f_16 * hsk_179[k]
                   + f_3 * pc_y[k] * isk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_y, pc_y, pc_z, hsl0_225, hsk_143, \
                         hsk_144, hsk_180, hsl1_225, isi0_223, isi1_223, isk_287, \
                         isk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * hsk_143[k]
                   + f_1 * isi0_223[k]
                   - f_2 * isi1_223[k]
                   + f_3 * pc_z[k] * isk_287[k];

        t_360[k] = pa_y[k] * hsl0_225[k]
                   - f_14 * pc_y[k] * hsl1_225[k];

        t_361[k] = f_15 * hsk_180[k]
                   + f_3 * pc_y[k] * isk_288[k];

        t_362[k] = f_16 * hsk_144[k]
                   + f_3 * pc_z[k] * isk_288[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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

    const auto *hsl0_228 = buffer.data(hsl0 + 228);
    const auto *hsl0_230 = buffer.data(hsl0 + 230);
    const auto *hsl0_231 = buffer.data(hsl0 + 231);
    const auto *hsl0_234 = buffer.data(hsl0 + 234);
    const auto *hsl0_235 = buffer.data(hsl0 + 235);
    const auto *hsl0_237 = buffer.data(hsl0 + 237);
    const auto *hsl0_239 = buffer.data(hsl0 + 239);
    const auto *hsl0_240 = buffer.data(hsl0 + 240);
    const auto *hsl0_242 = buffer.data(hsl0 + 242);
    const auto *hsl0_243 = buffer.data(hsl0 + 243);
    const auto *hsl0_245 = buffer.data(hsl0 + 245);
    const auto *hsl0_246 = buffer.data(hsl0 + 246);
    const auto *hsl0_248 = buffer.data(hsl0 + 248);
    const auto *hsl0_249 = buffer.data(hsl0 + 249);
    const auto *hsl0_250 = buffer.data(hsl0 + 250);
    const auto *hsl0_252 = buffer.data(hsl0 + 252);
    const auto *hsl0_269 = buffer.data(hsl0 + 269);

    const auto *hsk_147 = buffer.data(hsk + 147);
    const auto *hsk_150 = buffer.data(hsk + 150);
    const auto *hsk_154 = buffer.data(hsk + 154);
    const auto *hsk_159 = buffer.data(hsk + 159);
    const auto *hsk_172 = buffer.data(hsk + 172);
    const auto *hsk_180 = buffer.data(hsk + 180);
    const auto *hsk_181 = buffer.data(hsk + 181);
    const auto *hsk_182 = buffer.data(hsk + 182);
    const auto *hsk_183 = buffer.data(hsk + 183);
    const auto *hsk_185 = buffer.data(hsk + 185);
    const auto *hsk_186 = buffer.data(hsk + 186);
    const auto *hsk_188 = buffer.data(hsk + 188);
    const auto *hsk_189 = buffer.data(hsk + 189);
    const auto *hsk_190 = buffer.data(hsk + 190);
    const auto *hsk_192 = buffer.data(hsk + 192);
    const auto *hsk_193 = buffer.data(hsk + 193);
    const auto *hsk_194 = buffer.data(hsk + 194);
    const auto *hsk_195 = buffer.data(hsk + 195);
    const auto *hsk_197 = buffer.data(hsk + 197);
    const auto *hsk_198 = buffer.data(hsk + 198);
    const auto *hsk_199 = buffer.data(hsk + 199);
    const auto *hsk_200 = buffer.data(hsk + 200);
    const auto *hsk_208 = buffer.data(hsk + 208);
    const auto *hsk_210 = buffer.data(hsk + 210);
    const auto *hsk_211 = buffer.data(hsk + 211);
    const auto *hsk_212 = buffer.data(hsk + 212);
    const auto *hsk_213 = buffer.data(hsk + 213);
    const auto *hsk_214 = buffer.data(hsk + 214);
    const auto *hsk_215 = buffer.data(hsk + 215);
    const auto *hsk_216 = buffer.data(hsk + 216);
    const auto *hsk_221 = buffer.data(hsk + 221);
    const auto *hsk_225 = buffer.data(hsk + 225);
    const auto *hsk_230 = buffer.data(hsk + 230);
    const auto *hsk_236 = buffer.data(hsk + 236);
    const auto *hsk_316 = buffer.data(hsk + 316);
    const auto *hsk_317 = buffer.data(hsk + 317);
    const auto *hsk_318 = buffer.data(hsk + 318);
    const auto *hsk_319 = buffer.data(hsk + 319);
    const auto *hsk_320 = buffer.data(hsk + 320);
    const auto *hsk_321 = buffer.data(hsk + 321);
    const auto *hsk_322 = buffer.data(hsk + 322);
    const auto *hsk_323 = buffer.data(hsk + 323);
    const auto *hsk_324 = buffer.data(hsk + 324);
    const auto *hsk_329 = buffer.data(hsk + 329);
    const auto *hsk_333 = buffer.data(hsk + 333);
    const auto *hsk_338 = buffer.data(hsk + 338);
    const auto *hsk_344 = buffer.data(hsk + 344);
    const auto *hsk_351 = buffer.data(hsk + 351);
    const auto *hsk_352 = buffer.data(hsk + 352);
    const auto *hsk_353 = buffer.data(hsk + 353);
    const auto *hsk_354 = buffer.data(hsk + 354);
    const auto *hsk_355 = buffer.data(hsk + 355);
    const auto *hsk_356 = buffer.data(hsk + 356);
    const auto *hsk_357 = buffer.data(hsk + 357);
    const auto *hsk_359 = buffer.data(hsk + 359);
    const auto *hsk_360 = buffer.data(hsk + 360);
    const auto *hsk_363 = buffer.data(hsk + 363);
    const auto *hsk_366 = buffer.data(hsk + 366);
    const auto *hsk_370 = buffer.data(hsk + 370);
    const auto *hsk_375 = buffer.data(hsk + 375);
    const auto *hsk_381 = buffer.data(hsk + 381);

    const auto *hsl1_228 = buffer.data(hsl1 + 228);
    const auto *hsl1_230 = buffer.data(hsl1 + 230);
    const auto *hsl1_231 = buffer.data(hsl1 + 231);
    const auto *hsl1_234 = buffer.data(hsl1 + 234);
    const auto *hsl1_235 = buffer.data(hsl1 + 235);
    const auto *hsl1_237 = buffer.data(hsl1 + 237);
    const auto *hsl1_239 = buffer.data(hsl1 + 239);
    const auto *hsl1_240 = buffer.data(hsl1 + 240);
    const auto *hsl1_242 = buffer.data(hsl1 + 242);
    const auto *hsl1_243 = buffer.data(hsl1 + 243);
    const auto *hsl1_245 = buffer.data(hsl1 + 245);
    const auto *hsl1_246 = buffer.data(hsl1 + 246);
    const auto *hsl1_248 = buffer.data(hsl1 + 248);
    const auto *hsl1_249 = buffer.data(hsl1 + 249);
    const auto *hsl1_250 = buffer.data(hsl1 + 250);
    const auto *hsl1_252 = buffer.data(hsl1 + 252);
    const auto *hsl1_269 = buffer.data(hsl1 + 269);

    const auto *isi0_245 = buffer.data(isi0 + 245);
    const auto *isi0_247 = buffer.data(isi0 + 247);
    const auto *isi0_248 = buffer.data(isi0 + 248);
    const auto *isi0_249 = buffer.data(isi0 + 249);
    const auto *isi0_250 = buffer.data(isi0 + 250);
    const auto *isi0_251 = buffer.data(isi0 + 251);
    const auto *isi0_252 = buffer.data(isi0 + 252);
    const auto *isi0_253 = buffer.data(isi0 + 253);
    const auto *isi0_254 = buffer.data(isi0 + 254);
    const auto *isi0_255 = buffer.data(isi0 + 255);
    const auto *isi0_256 = buffer.data(isi0 + 256);
    const auto *isi0_257 = buffer.data(isi0 + 257);
    const auto *isi0_258 = buffer.data(isi0 + 258);
    const auto *isi0_259 = buffer.data(isi0 + 259);
    const auto *isi0_260 = buffer.data(isi0 + 260);
    const auto *isi0_261 = buffer.data(isi0 + 261);
    const auto *isi0_262 = buffer.data(isi0 + 262);
    const auto *isi0_263 = buffer.data(isi0 + 263);
    const auto *isi0_264 = buffer.data(isi0 + 264);
    const auto *isi0_265 = buffer.data(isi0 + 265);
    const auto *isi0_266 = buffer.data(isi0 + 266);
    const auto *isi0_272 = buffer.data(isi0 + 272);
    const auto *isi0_273 = buffer.data(isi0 + 273);
    const auto *isi0_274 = buffer.data(isi0 + 274);
    const auto *isi0_275 = buffer.data(isi0 + 275);
    const auto *isi0_276 = buffer.data(isi0 + 276);
    const auto *isi0_277 = buffer.data(isi0 + 277);
    const auto *isi0_278 = buffer.data(isi0 + 278);
    const auto *isi0_279 = buffer.data(isi0 + 279);
    const auto *isi0_280 = buffer.data(isi0 + 280);
    const auto *isi0_282 = buffer.data(isi0 + 282);
    const auto *isi0_283 = buffer.data(isi0 + 283);
    const auto *isi0_285 = buffer.data(isi0 + 285);
    const auto *isi0_286 = buffer.data(isi0 + 286);
    const auto *isi0_287 = buffer.data(isi0 + 287);
    const auto *isi0_289 = buffer.data(isi0 + 289);
    const auto *isi0_290 = buffer.data(isi0 + 290);
    const auto *isi0_291 = buffer.data(isi0 + 291);
    const auto *isi0_292 = buffer.data(isi0 + 292);
    const auto *isi0_294 = buffer.data(isi0 + 294);
    const auto *isi0_295 = buffer.data(isi0 + 295);
    const auto *isi0_301 = buffer.data(isi0 + 301);

    const auto *isi1_245 = buffer.data(isi1 + 245);
    const auto *isi1_247 = buffer.data(isi1 + 247);
    const auto *isi1_248 = buffer.data(isi1 + 248);
    const auto *isi1_249 = buffer.data(isi1 + 249);
    const auto *isi1_250 = buffer.data(isi1 + 250);
    const auto *isi1_251 = buffer.data(isi1 + 251);
    const auto *isi1_252 = buffer.data(isi1 + 252);
    const auto *isi1_253 = buffer.data(isi1 + 253);
    const auto *isi1_254 = buffer.data(isi1 + 254);
    const auto *isi1_255 = buffer.data(isi1 + 255);
    const auto *isi1_256 = buffer.data(isi1 + 256);
    const auto *isi1_257 = buffer.data(isi1 + 257);
    const auto *isi1_258 = buffer.data(isi1 + 258);
    const auto *isi1_259 = buffer.data(isi1 + 259);
    const auto *isi1_260 = buffer.data(isi1 + 260);
    const auto *isi1_261 = buffer.data(isi1 + 261);
    const auto *isi1_262 = buffer.data(isi1 + 262);
    const auto *isi1_263 = buffer.data(isi1 + 263);
    const auto *isi1_264 = buffer.data(isi1 + 264);
    const auto *isi1_265 = buffer.data(isi1 + 265);
    const auto *isi1_266 = buffer.data(isi1 + 266);
    const auto *isi1_272 = buffer.data(isi1 + 272);
    const auto *isi1_273 = buffer.data(isi1 + 273);
    const auto *isi1_274 = buffer.data(isi1 + 274);
    const auto *isi1_275 = buffer.data(isi1 + 275);
    const auto *isi1_276 = buffer.data(isi1 + 276);
    const auto *isi1_277 = buffer.data(isi1 + 277);
    const auto *isi1_278 = buffer.data(isi1 + 278);
    const auto *isi1_279 = buffer.data(isi1 + 279);
    const auto *isi1_280 = buffer.data(isi1 + 280);
    const auto *isi1_282 = buffer.data(isi1 + 282);
    const auto *isi1_283 = buffer.data(isi1 + 283);
    const auto *isi1_285 = buffer.data(isi1 + 285);
    const auto *isi1_286 = buffer.data(isi1 + 286);
    const auto *isi1_287 = buffer.data(isi1 + 287);
    const auto *isi1_289 = buffer.data(isi1 + 289);
    const auto *isi1_290 = buffer.data(isi1 + 290);
    const auto *isi1_291 = buffer.data(isi1 + 291);
    const auto *isi1_292 = buffer.data(isi1 + 292);
    const auto *isi1_294 = buffer.data(isi1 + 294);
    const auto *isi1_295 = buffer.data(isi1 + 295);
    const auto *isi1_301 = buffer.data(isi1 + 301);

    const auto *isk_290 = buffer.data(isk + 290);
    const auto *isk_291 = buffer.data(isk + 291);
    const auto *isk_293 = buffer.data(isk + 293);
    const auto *isk_294 = buffer.data(isk + 294);
    const auto *isk_297 = buffer.data(isk + 297);
    const auto *isk_298 = buffer.data(isk + 298);
    const auto *isk_302 = buffer.data(isk + 302);
    const auto *isk_303 = buffer.data(isk + 303);
    const auto *isk_308 = buffer.data(isk + 308);
    const auto *isk_316 = buffer.data(isk + 316);
    const auto *isk_317 = buffer.data(isk + 317);
    const auto *isk_318 = buffer.data(isk + 318);
    const auto *isk_319 = buffer.data(isk + 319);
    const auto *isk_320 = buffer.data(isk + 320);
    const auto *isk_321 = buffer.data(isk + 321);
    const auto *isk_322 = buffer.data(isk + 322);
    const auto *isk_323 = buffer.data(isk + 323);
    const auto *isk_324 = buffer.data(isk + 324);
    const auto *isk_325 = buffer.data(isk + 325);
    const auto *isk_326 = buffer.data(isk + 326);
    const auto *isk_327 = buffer.data(isk + 327);
    const auto *isk_328 = buffer.data(isk + 328);
    const auto *isk_329 = buffer.data(isk + 329);
    const auto *isk_330 = buffer.data(isk + 330);
    const auto *isk_331 = buffer.data(isk + 331);
    const auto *isk_332 = buffer.data(isk + 332);
    const auto *isk_333 = buffer.data(isk + 333);
    const auto *isk_334 = buffer.data(isk + 334);
    const auto *isk_335 = buffer.data(isk + 335);
    const auto *isk_336 = buffer.data(isk + 336);
    const auto *isk_337 = buffer.data(isk + 337);
    const auto *isk_338 = buffer.data(isk + 338);
    const auto *isk_339 = buffer.data(isk + 339);
    const auto *isk_340 = buffer.data(isk + 340);
    const auto *isk_341 = buffer.data(isk + 341);
    const auto *isk_342 = buffer.data(isk + 342);
    const auto *isk_343 = buffer.data(isk + 343);
    const auto *isk_344 = buffer.data(isk + 344);
    const auto *isk_351 = buffer.data(isk + 351);
    const auto *isk_352 = buffer.data(isk + 352);
    const auto *isk_353 = buffer.data(isk + 353);
    const auto *isk_354 = buffer.data(isk + 354);
    const auto *isk_355 = buffer.data(isk + 355);
    const auto *isk_356 = buffer.data(isk + 356);
    const auto *isk_357 = buffer.data(isk + 357);
    const auto *isk_358 = buffer.data(isk + 358);
    const auto *isk_359 = buffer.data(isk + 359);
    const auto *isk_360 = buffer.data(isk + 360);
    const auto *isk_361 = buffer.data(isk + 361);
    const auto *isk_362 = buffer.data(isk + 362);
    const auto *isk_363 = buffer.data(isk + 363);
    const auto *isk_365 = buffer.data(isk + 365);
    const auto *isk_366 = buffer.data(isk + 366);
    const auto *isk_367 = buffer.data(isk + 367);
    const auto *isk_369 = buffer.data(isk + 369);
    const auto *isk_370 = buffer.data(isk + 370);
    const auto *isk_371 = buffer.data(isk + 371);
    const auto *isk_372 = buffer.data(isk + 372);
    const auto *isk_374 = buffer.data(isk + 374);
    const auto *isk_375 = buffer.data(isk + 375);
    const auto *isk_376 = buffer.data(isk + 376);
    const auto *isk_377 = buffer.data(isk + 377);
    const auto *isk_378 = buffer.data(isk + 378);
    const auto *isk_380 = buffer.data(isk + 380);
    const auto *isk_381 = buffer.data(isk + 381);

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_y, pc_y, hsl0_228, hsl0_230, hsl0_231, \
                         hsk_181, hsk_182, hsk_183, hsl1_228, hsl1_230, hsl1_231, \
                         isk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_y[k] * hsl0_228[k]
                   + f_16 * hsk_181[k]
                   - f_14 * pc_y[k] * hsl1_228[k];

        t_364[k] = f_15 * hsk_182[k]
                   + f_3 * pc_y[k] * isk_290[k];

        t_365[k] = pa_y[k] * hsl0_230[k]
                   - f_14 * pc_y[k] * hsl1_230[k];

        t_366[k] = pa_y[k] * hsl0_231[k]
                   + f_17 * hsk_183[k]
                   - f_14 * pc_y[k] * hsl1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_y, pc_y, pc_z, hsl0_234, hsl0_235, \
                         hsk_147, hsk_185, hsk_186, hsl1_234, hsl1_235, isk_291, \
                         isk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * hsk_147[k]
                   + f_3 * pc_z[k] * isk_291[k];

        t_368[k] = f_15 * hsk_185[k]
                   + f_3 * pc_y[k] * isk_293[k];

        t_369[k] = pa_y[k] * hsl0_234[k]
                   - f_14 * pc_y[k] * hsl1_234[k];

        t_370[k] = pa_y[k] * hsl0_235[k]
                   + f_18 * hsk_186[k]
                   - f_14 * pc_y[k] * hsl1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pc_y, pc_z, hsl0_237, hsl0_239, \
                         hsk_150, hsk_188, hsk_189, hsl1_237, hsl1_239, isk_294, \
                         isk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * hsk_150[k]
                   + f_3 * pc_z[k] * isk_294[k];

        t_372[k] = pa_y[k] * hsl0_237[k]
                   + f_16 * hsk_188[k]
                   - f_14 * pc_y[k] * hsl1_237[k];

        t_373[k] = f_15 * hsk_189[k]
                   + f_3 * pc_y[k] * isk_297[k];

        t_374[k] = pa_y[k] * hsl0_239[k]
                   - f_14 * pc_y[k] * hsl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_y, pc_y, pc_z, hsl0_240, hsl0_242, hsk_154, \
                         hsk_190, hsk_192, hsl1_240, hsl1_242, \
                         isk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pa_y[k] * hsl0_240[k]
                   + f_19 * hsk_190[k]
                   - f_14 * pc_y[k] * hsl1_240[k];

        t_376[k] = f_16 * hsk_154[k]
                   + f_3 * pc_z[k] * isk_298[k];

        t_377[k] = pa_y[k] * hsl0_242[k]
                   + f_17 * hsk_192[k]
                   - f_14 * pc_y[k] * hsl1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_y, pc_y, hsl0_243, hsl0_245, hsl0_246, \
                         hsk_193, hsk_194, hsk_195, hsl1_243, hsl1_245, hsl1_246, \
                         isk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * hsl0_243[k]
                   + f_16 * hsk_193[k]
                   - f_14 * pc_y[k] * hsl1_243[k];

        t_379[k] = f_15 * hsk_194[k]
                   + f_3 * pc_y[k] * isk_302[k];

        t_380[k] = pa_y[k] * hsl0_245[k]
                   - f_14 * pc_y[k] * hsl1_245[k];

        t_381[k] = pa_y[k] * hsl0_246[k]
                   + f_0 * hsk_195[k]
                   - f_14 * pc_y[k] * hsl1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pa_y, pc_y, pc_z, hsl0_248, hsl0_249, hsk_159, \
                         hsk_197, hsk_198, hsl1_248, hsl1_249, \
                         isk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * hsk_159[k]
                   + f_3 * pc_z[k] * isk_303[k];

        t_383[k] = pa_y[k] * hsl0_248[k]
                   + f_18 * hsk_197[k]
                   - f_14 * pc_y[k] * hsl1_248[k];

        t_384[k] = pa_y[k] * hsl0_249[k]
                   + f_17 * hsk_198[k]
                   - f_14 * pc_y[k] * hsl1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_y, pc_x, pc_y, hsl0_250, hsl0_252, \
                         hsk_199, hsk_200, hsk_316, hsl1_250, hsl1_252, isk_308, \
                         isk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pa_y[k] * hsl0_250[k]
                   + f_16 * hsk_199[k]
                   - f_14 * pc_y[k] * hsl1_250[k];

        t_386[k] = f_15 * hsk_200[k]
                   + f_3 * pc_y[k] * isk_308[k];

        t_387[k] = pa_y[k] * hsl0_252[k]
                   - f_14 * pc_y[k] * hsl1_252[k];

        t_388[k] = f_17 * hsk_316[k]
                   + f_3 * pc_x[k] * isk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, hsk_317, hsk_318, hsk_319, \
                         hsk_320, hsk_321, isk_317, isk_318, isk_319, isk_320, \
                         isk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_17 * hsk_317[k]
                   + f_3 * pc_x[k] * isk_317[k];

        t_390[k] = f_17 * hsk_318[k]
                   + f_3 * pc_x[k] * isk_318[k];

        t_391[k] = f_17 * hsk_319[k]
                   + f_3 * pc_x[k] * isk_319[k];

        t_392[k] = f_17 * hsk_320[k]
                   + f_3 * pc_x[k] * isk_320[k];

        t_393[k] = f_17 * hsk_321[k]
                   + f_3 * pc_x[k] * isk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, hsk_172, hsk_208, \
                         hsk_322, hsk_323, isi0_245, isi1_245, isk_316, isk_322, \
                         isk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_17 * hsk_322[k]
                   + f_3 * pc_x[k] * isk_322[k];

        t_395[k] = f_17 * hsk_323[k]
                   + f_3 * pc_x[k] * isk_323[k];

        t_396[k] = f_15 * hsk_208[k]
                   + f_1 * isi0_245[k]
                   - f_2 * isi1_245[k]
                   + f_3 * pc_y[k] * isk_316[k];

        t_397[k] = f_16 * hsk_172[k]
                   + f_3 * pc_z[k] * isk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, hsk_210, hsk_211, hsk_212, isi0_247, \
                         isi0_248, isi0_249, isi1_247, isi1_248, isi1_249, isk_318, isk_319, \
                         isk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * hsk_210[k]
                   + f_12 * isi0_247[k]
                   - f_13 * isi1_247[k]
                   + f_3 * pc_y[k] * isk_318[k];

        t_399[k] = f_15 * hsk_211[k]
                   + f_10 * isi0_248[k]
                   - f_11 * isi1_248[k]
                   + f_3 * pc_y[k] * isk_319[k];

        t_400[k] = f_15 * hsk_212[k]
                   + f_8 * isi0_249[k]
                   - f_9 * isi1_249[k]
                   + f_3 * pc_y[k] * isk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, hsk_213, hsk_214, hsk_215, isi0_250, \
                         isi0_251, isi1_250, isi1_251, isk_321, isk_322, \
                         isk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * hsk_213[k]
                   + f_6 * isi0_250[k]
                   - f_7 * isi1_250[k]
                   + f_3 * pc_y[k] * isk_321[k];

        t_402[k] = f_15 * hsk_214[k]
                   + f_4 * isi0_251[k]
                   - f_5 * isi1_251[k]
                   + f_3 * pc_y[k] * isk_322[k];

        t_403[k] = f_15 * hsk_215[k]
                   + f_3 * pc_y[k] * isk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_y, pc_x, pc_y, pc_z, hsl0_269, \
                         hsk_180, hsk_324, hsl1_269, isi0_252, isi1_252, \
                         isk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * hsl0_269[k]
                   - f_14 * pc_y[k] * hsl1_269[k];

        t_405[k] = f_17 * hsk_324[k]
                   + f_1 * isi0_252[k]
                   - f_2 * isi1_252[k]
                   + f_3 * pc_x[k] * isk_324[k];

        t_406[k] = f_3 * pc_y[k] * isk_324[k];

        t_407[k] = f_17 * hsk_180[k]
                   + f_3 * pc_z[k] * isk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, hsk_329, isi0_252, isi0_257, \
                         isi1_252, isi1_257, isk_325, isk_326, \
                         isk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_4 * isi0_252[k]
                   - f_5 * isi1_252[k]
                   + f_3 * pc_y[k] * isk_325[k];

        t_409[k] = f_3 * pc_y[k] * isk_326[k];

        t_410[k] = f_17 * hsk_329[k]
                   + f_12 * isi0_257[k]
                   - f_13 * isi1_257[k]
                   + f_3 * pc_x[k] * isk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_y, isi0_253, isi0_254, isi1_253, isi1_254, \
                         isk_327, isk_328, isk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_6 * isi0_253[k]
                   - f_7 * isi1_253[k]
                   + f_3 * pc_y[k] * isk_327[k];

        t_412[k] = f_4 * isi0_254[k]
                   - f_5 * isi1_254[k]
                   + f_3 * pc_y[k] * isk_328[k];

        t_413[k] = f_3 * pc_y[k] * isk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, hsk_333, isi0_255, isi0_256, \
                         isi0_261, isi1_255, isi1_256, isi1_261, isk_330, isk_331, \
                         isk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_17 * hsk_333[k]
                   + f_10 * isi0_261[k]
                   - f_11 * isi1_261[k]
                   + f_3 * pc_x[k] * isk_333[k];

        t_415[k] = f_8 * isi0_255[k]
                   - f_9 * isi1_255[k]
                   + f_3 * pc_y[k] * isk_330[k];

        t_416[k] = f_6 * isi0_256[k]
                   - f_7 * isi1_256[k]
                   + f_3 * pc_y[k] * isk_331[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, hsk_338, isi0_257, isi0_266, \
                         isi1_257, isi1_266, isk_332, isk_333, \
                         isk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * isi0_257[k]
                   - f_5 * isi1_257[k]
                   + f_3 * pc_y[k] * isk_332[k];

        t_418[k] = f_3 * pc_y[k] * isk_333[k];

        t_419[k] = f_17 * hsk_338[k]
                   + f_8 * isi0_266[k]
                   - f_9 * isi1_266[k]
                   + f_3 * pc_x[k] * isk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_y, isi0_258, isi0_259, isi0_260, isi1_258, \
                         isi1_259, isi1_260, isk_334, isk_335, \
                         isk_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_10 * isi0_258[k]
                   - f_11 * isi1_258[k]
                   + f_3 * pc_y[k] * isk_334[k];

        t_421[k] = f_8 * isi0_259[k]
                   - f_9 * isi1_259[k]
                   + f_3 * pc_y[k] * isk_335[k];

        t_422[k] = f_6 * isi0_260[k]
                   - f_7 * isi1_260[k]
                   + f_3 * pc_y[k] * isk_336[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, hsk_344, isi0_261, isi0_272, \
                         isi1_261, isi1_272, isk_337, isk_338, \
                         isk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * isi0_261[k]
                   - f_5 * isi1_261[k]
                   + f_3 * pc_y[k] * isk_337[k];

        t_424[k] = f_3 * pc_y[k] * isk_338[k];

        t_425[k] = f_17 * hsk_344[k]
                   + f_6 * isi0_272[k]
                   - f_7 * isi1_272[k]
                   + f_3 * pc_x[k] * isk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_y, isi0_262, isi0_263, isi0_264, isi1_262, \
                         isi1_263, isi1_264, isk_339, isk_340, \
                         isk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_12 * isi0_262[k]
                   - f_13 * isi1_262[k]
                   + f_3 * pc_y[k] * isk_339[k];

        t_427[k] = f_10 * isi0_263[k]
                   - f_11 * isi1_263[k]
                   + f_3 * pc_y[k] * isk_340[k];

        t_428[k] = f_8 * isi0_264[k]
                   - f_9 * isi1_264[k]
                   + f_3 * pc_y[k] * isk_341[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_y, isi0_265, isi0_266, isi1_265, isi1_266, \
                         isk_342, isk_343, isk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_6 * isi0_265[k]
                   - f_7 * isi1_265[k]
                   + f_3 * pc_y[k] * isk_342[k];

        t_430[k] = f_4 * isi0_266[k]
                   - f_5 * isi1_266[k]
                   + f_3 * pc_y[k] * isk_343[k];

        t_431[k] = f_3 * pc_y[k] * isk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, hsk_351, hsk_352, hsk_353, hsk_354, \
                         isi0_279, isi1_279, isk_351, isk_352, isk_353, \
                         isk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_17 * hsk_351[k]
                   + f_4 * isi0_279[k]
                   - f_5 * isi1_279[k]
                   + f_3 * pc_x[k] * isk_351[k];

        t_433[k] = f_17 * hsk_352[k]
                   + f_3 * pc_x[k] * isk_352[k];

        t_434[k] = f_17 * hsk_353[k]
                   + f_3 * pc_x[k] * isk_353[k];

        t_435[k] = f_17 * hsk_354[k]
                   + f_3 * pc_x[k] * isk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_y, hsk_355, hsk_356, \
                         hsk_357, hsk_359, isk_351, isk_355, isk_356, isk_357, \
                         isk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_17 * hsk_355[k]
                   + f_3 * pc_x[k] * isk_355[k];

        t_437[k] = f_17 * hsk_356[k]
                   + f_3 * pc_x[k] * isk_356[k];

        t_438[k] = f_17 * hsk_357[k]
                   + f_3 * pc_x[k] * isk_357[k];

        t_439[k] = f_3 * pc_y[k] * isk_351[k];

        t_440[k] = f_17 * hsk_359[k]
                   + f_3 * pc_x[k] * isk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, isi0_273, isi0_274, isi0_275, isi1_273, \
                         isi1_274, isi1_275, isk_352, isk_353, \
                         isk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * isi0_273[k]
                   - f_2 * isi1_273[k]
                   + f_3 * pc_y[k] * isk_352[k];

        t_442[k] = f_20 * isi0_274[k]
                   - f_21 * isi1_274[k]
                   + f_3 * pc_y[k] * isk_353[k];

        t_443[k] = f_12 * isi0_275[k]
                   - f_13 * isi1_275[k]
                   + f_3 * pc_y[k] * isk_354[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_y, isi0_276, isi0_277, isi0_278, isi1_276, \
                         isi1_277, isi1_278, isk_355, isk_356, \
                         isk_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_10 * isi0_276[k]
                   - f_11 * isi1_276[k]
                   + f_3 * pc_y[k] * isk_355[k];

        t_445[k] = f_8 * isi0_277[k]
                   - f_9 * isi1_277[k]
                   + f_3 * pc_y[k] * isk_356[k];

        t_446[k] = f_6 * isi0_278[k]
                   - f_7 * isi1_278[k]
                   + f_3 * pc_y[k] * isk_357[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pc_x, pc_y, pc_z, hsk_215, hsk_360, \
                         isi0_279, isi0_280, isi1_279, isi1_280, isk_358, isk_359, \
                         isk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_4 * isi0_279[k]
                   - f_5 * isi1_279[k]
                   + f_3 * pc_y[k] * isk_358[k];

        t_448[k] = f_3 * pc_y[k] * isk_359[k];

        t_449[k] = f_17 * hsk_215[k]
                   + f_1 * isi0_279[k]
                   - f_2 * isi1_279[k]
                   + f_3 * pc_z[k] * isk_359[k];

        t_450[k] = f_16 * hsk_360[k]
                   + f_1 * isi0_280[k]
                   - f_2 * isi1_280[k]
                   + f_3 * pc_x[k] * isk_360[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pc_x, pc_y, pc_z, hsk_216, hsk_363, \
                         isi0_283, isi1_283, isk_360, isk_361, \
                         isk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_18 * hsk_216[k]
                   + f_3 * pc_y[k] * isk_360[k];

        t_452[k] = f_3 * pc_z[k] * isk_360[k];

        t_453[k] = f_16 * hsk_363[k]
                   + f_12 * isi0_283[k]
                   - f_13 * isi1_283[k]
                   + f_3 * pc_x[k] * isk_363[k];

        t_454[k] = f_3 * pc_z[k] * isk_361[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_z, hsk_366, isi0_280, isi0_286, \
                         isi1_280, isi1_286, isk_362, isk_363, \
                         isk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_4 * isi0_280[k]
                   - f_5 * isi1_280[k]
                   + f_3 * pc_z[k] * isk_362[k];

        t_456[k] = f_16 * hsk_366[k]
                   + f_10 * isi0_286[k]
                   - f_11 * isi1_286[k]
                   + f_3 * pc_x[k] * isk_366[k];

        t_457[k] = f_3 * pc_z[k] * isk_363[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pc_x, pc_y, pc_z, hsk_221, hsk_370, \
                         isi0_282, isi0_290, isi1_282, isi1_290, isk_365, isk_366, \
                         isk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_18 * hsk_221[k]
                   + f_3 * pc_y[k] * isk_365[k];

        t_459[k] = f_6 * isi0_282[k]
                   - f_7 * isi1_282[k]
                   + f_3 * pc_z[k] * isk_365[k];

        t_460[k] = f_16 * hsk_370[k]
                   + f_8 * isi0_290[k]
                   - f_9 * isi1_290[k]
                   + f_3 * pc_x[k] * isk_370[k];

        t_461[k] = f_3 * pc_z[k] * isk_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, hsk_225, isi0_283, isi0_285, \
                         isi1_283, isi1_285, isk_367, isk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_4 * isi0_283[k]
                   - f_5 * isi1_283[k]
                   + f_3 * pc_z[k] * isk_367[k];

        t_463[k] = f_18 * hsk_225[k]
                   + f_3 * pc_y[k] * isk_369[k];

        t_464[k] = f_8 * isi0_285[k]
                   - f_9 * isi1_285[k]
                   + f_3 * pc_z[k] * isk_369[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_z, hsk_375, isi0_286, isi0_295, \
                         isi1_286, isi1_295, isk_370, isk_371, \
                         isk_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_16 * hsk_375[k]
                   + f_6 * isi0_295[k]
                   - f_7 * isi1_295[k]
                   + f_3 * pc_x[k] * isk_375[k];

        t_466[k] = f_3 * pc_z[k] * isk_370[k];

        t_467[k] = f_4 * isi0_286[k]
                   - f_5 * isi1_286[k]
                   + f_3 * pc_z[k] * isk_371[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_y, pc_z, hsk_230, isi0_287, isi0_289, \
                         isi1_287, isi1_289, isk_372, isk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_6 * isi0_287[k]
                   - f_7 * isi1_287[k]
                   + f_3 * pc_z[k] * isk_372[k];

        t_469[k] = f_18 * hsk_230[k]
                   + f_3 * pc_y[k] * isk_374[k];

        t_470[k] = f_10 * isi0_289[k]
                   - f_11 * isi1_289[k]
                   + f_3 * pc_z[k] * isk_374[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_x, pc_z, hsk_381, isi0_290, isi0_301, \
                         isi1_290, isi1_301, isk_375, isk_376, \
                         isk_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_16 * hsk_381[k]
                   + f_4 * isi0_301[k]
                   - f_5 * isi1_301[k]
                   + f_3 * pc_x[k] * isk_381[k];

        t_472[k] = f_3 * pc_z[k] * isk_375[k];

        t_473[k] = f_4 * isi0_290[k]
                   - f_5 * isi1_290[k]
                   + f_3 * pc_z[k] * isk_376[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_y, pc_z, hsk_236, isi0_291, isi0_292, \
                         isi0_294, isi1_291, isi1_292, isi1_294, isk_377, isk_378, \
                         isk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_6 * isi0_291[k]
                   - f_7 * isi1_291[k]
                   + f_3 * pc_z[k] * isk_377[k];

        t_475[k] = f_8 * isi0_292[k]
                   - f_9 * isi1_292[k]
                   + f_3 * pc_z[k] * isk_378[k];

        t_476[k] = f_18 * hsk_236[k]
                   + f_3 * pc_y[k] * isk_380[k];

        t_477[k] = f_12 * isi0_294[k]
                   - f_13 * isi1_294[k]
                   + f_3 * pc_z[k] * isk_380[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
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

    const auto *hsl0_270 = buffer.data(hsl0 + 270);
    const auto *hsl0_273 = buffer.data(hsl0 + 273);
    const auto *hsl0_276 = buffer.data(hsl0 + 276);
    const auto *hsl0_280 = buffer.data(hsl0 + 280);
    const auto *hsl0_282 = buffer.data(hsl0 + 282);
    const auto *hsl0_285 = buffer.data(hsl0 + 285);
    const auto *hsl0_287 = buffer.data(hsl0 + 287);
    const auto *hsl0_288 = buffer.data(hsl0 + 288);
    const auto *hsl0_291 = buffer.data(hsl0 + 291);
    const auto *hsl0_293 = buffer.data(hsl0 + 293);
    const auto *hsl0_294 = buffer.data(hsl0 + 294);
    const auto *hsl0_295 = buffer.data(hsl0 + 295);
    const auto *hsl0_306 = buffer.data(hsl0 + 306);
    const auto *hsl0_405 = buffer.data(hsl0 + 405);

    const auto *hsk_216 = buffer.data(hsk + 216);
    const auto *hsk_219 = buffer.data(hsk + 219);
    const auto *hsk_222 = buffer.data(hsk + 222);
    const auto *hsk_223 = buffer.data(hsk + 223);
    const auto *hsk_226 = buffer.data(hsk + 226);
    const auto *hsk_227 = buffer.data(hsk + 227);
    const auto *hsk_228 = buffer.data(hsk + 228);
    const auto *hsk_231 = buffer.data(hsk + 231);
    const auto *hsk_232 = buffer.data(hsk + 232);
    const auto *hsk_233 = buffer.data(hsk + 233);
    const auto *hsk_234 = buffer.data(hsk + 234);
    const auto *hsk_244 = buffer.data(hsk + 244);
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
    const auto *hsk_280 = buffer.data(hsk + 280);
    const auto *hsk_282 = buffer.data(hsk + 282);
    const auto *hsk_283 = buffer.data(hsk + 283);
    const auto *hsk_284 = buffer.data(hsk + 284);
    const auto *hsk_285 = buffer.data(hsk + 285);
    const auto *hsk_286 = buffer.data(hsk + 286);
    const auto *hsk_287 = buffer.data(hsk + 287);
    const auto *hsk_288 = buffer.data(hsk + 288);
    const auto *hsk_290 = buffer.data(hsk + 290);
    const auto *hsk_293 = buffer.data(hsk + 293);
    const auto *hsk_297 = buffer.data(hsk + 297);
    const auto *hsk_302 = buffer.data(hsk + 302);
    const auto *hsk_308 = buffer.data(hsk + 308);
    const auto *hsk_316 = buffer.data(hsk + 316);
    const auto *hsk_318 = buffer.data(hsk + 318);
    const auto *hsk_319 = buffer.data(hsk + 319);
    const auto *hsk_320 = buffer.data(hsk + 320);
    const auto *hsk_321 = buffer.data(hsk + 321);
    const auto *hsk_322 = buffer.data(hsk + 322);
    const auto *hsk_323 = buffer.data(hsk + 323);
    const auto *hsk_324 = buffer.data(hsk + 324);
    const auto *hsk_388 = buffer.data(hsk + 388);
    const auto *hsk_390 = buffer.data(hsk + 390);
    const auto *hsk_391 = buffer.data(hsk + 391);
    const auto *hsk_392 = buffer.data(hsk + 392);
    const auto *hsk_393 = buffer.data(hsk + 393);
    const auto *hsk_394 = buffer.data(hsk + 394);
    const auto *hsk_395 = buffer.data(hsk + 395);
    const auto *hsk_401 = buffer.data(hsk + 401);
    const auto *hsk_405 = buffer.data(hsk + 405);
    const auto *hsk_410 = buffer.data(hsk + 410);
    const auto *hsk_416 = buffer.data(hsk + 416);
    const auto *hsk_423 = buffer.data(hsk + 423);
    const auto *hsk_424 = buffer.data(hsk + 424);
    const auto *hsk_425 = buffer.data(hsk + 425);
    const auto *hsk_426 = buffer.data(hsk + 426);
    const auto *hsk_427 = buffer.data(hsk + 427);
    const auto *hsk_428 = buffer.data(hsk + 428);
    const auto *hsk_429 = buffer.data(hsk + 429);
    const auto *hsk_430 = buffer.data(hsk + 430);
    const auto *hsk_431 = buffer.data(hsk + 431);
    const auto *hsk_432 = buffer.data(hsk + 432);
    const auto *hsk_435 = buffer.data(hsk + 435);
    const auto *hsk_437 = buffer.data(hsk + 437);
    const auto *hsk_438 = buffer.data(hsk + 438);
    const auto *hsk_441 = buffer.data(hsk + 441);
    const auto *hsk_442 = buffer.data(hsk + 442);
    const auto *hsk_444 = buffer.data(hsk + 444);
    const auto *hsk_446 = buffer.data(hsk + 446);
    const auto *hsk_447 = buffer.data(hsk + 447);
    const auto *hsk_449 = buffer.data(hsk + 449);
    const auto *hsk_450 = buffer.data(hsk + 450);
    const auto *hsk_452 = buffer.data(hsk + 452);
    const auto *hsk_453 = buffer.data(hsk + 453);
    const auto *hsk_455 = buffer.data(hsk + 455);
    const auto *hsk_456 = buffer.data(hsk + 456);
    const auto *hsk_457 = buffer.data(hsk + 457);
    const auto *hsk_459 = buffer.data(hsk + 459);
    const auto *hsk_460 = buffer.data(hsk + 460);
    const auto *hsk_461 = buffer.data(hsk + 461);
    const auto *hsk_462 = buffer.data(hsk + 462);
    const auto *hsk_463 = buffer.data(hsk + 463);
    const auto *hsk_464 = buffer.data(hsk + 464);
    const auto *hsk_465 = buffer.data(hsk + 465);
    const auto *hsk_466 = buffer.data(hsk + 466);
    const auto *hsk_467 = buffer.data(hsk + 467);

    const auto *hsl1_270 = buffer.data(hsl1 + 270);
    const auto *hsl1_273 = buffer.data(hsl1 + 273);
    const auto *hsl1_276 = buffer.data(hsl1 + 276);
    const auto *hsl1_280 = buffer.data(hsl1 + 280);
    const auto *hsl1_282 = buffer.data(hsl1 + 282);
    const auto *hsl1_285 = buffer.data(hsl1 + 285);
    const auto *hsl1_287 = buffer.data(hsl1 + 287);
    const auto *hsl1_288 = buffer.data(hsl1 + 288);
    const auto *hsl1_291 = buffer.data(hsl1 + 291);
    const auto *hsl1_293 = buffer.data(hsl1 + 293);
    const auto *hsl1_294 = buffer.data(hsl1 + 294);
    const auto *hsl1_295 = buffer.data(hsl1 + 295);
    const auto *hsl1_306 = buffer.data(hsl1 + 306);
    const auto *hsl1_405 = buffer.data(hsl1 + 405);

    const auto *isi0_301 = buffer.data(isi0 + 301);
    const auto *isi0_302 = buffer.data(isi0 + 302);
    const auto *isi0_303 = buffer.data(isi0 + 303);
    const auto *isi0_304 = buffer.data(isi0 + 304);
    const auto *isi0_305 = buffer.data(isi0 + 305);
    const auto *isi0_307 = buffer.data(isi0 + 307);
    const auto *isi0_313 = buffer.data(isi0 + 313);
    const auto *isi0_317 = buffer.data(isi0 + 317);
    const auto *isi0_322 = buffer.data(isi0 + 322);
    const auto *isi0_328 = buffer.data(isi0 + 328);
    const auto *isi0_331 = buffer.data(isi0 + 331);
    const auto *isi0_332 = buffer.data(isi0 + 332);
    const auto *isi0_333 = buffer.data(isi0 + 333);
    const auto *isi0_334 = buffer.data(isi0 + 334);
    const auto *isi0_335 = buffer.data(isi0 + 335);
    const auto *isi0_336 = buffer.data(isi0 + 336);
    const auto *isi0_339 = buffer.data(isi0 + 339);
    const auto *isi0_341 = buffer.data(isi0 + 341);
    const auto *isi0_342 = buffer.data(isi0 + 342);
    const auto *isi0_345 = buffer.data(isi0 + 345);
    const auto *isi0_346 = buffer.data(isi0 + 346);
    const auto *isi0_348 = buffer.data(isi0 + 348);
    const auto *isi0_350 = buffer.data(isi0 + 350);
    const auto *isi0_351 = buffer.data(isi0 + 351);
    const auto *isi0_353 = buffer.data(isi0 + 353);
    const auto *isi0_354 = buffer.data(isi0 + 354);
    const auto *isi0_356 = buffer.data(isi0 + 356);
    const auto *isi0_357 = buffer.data(isi0 + 357);
    const auto *isi0_359 = buffer.data(isi0 + 359);
    const auto *isi0_360 = buffer.data(isi0 + 360);
    const auto *isi0_361 = buffer.data(isi0 + 361);
    const auto *isi0_362 = buffer.data(isi0 + 362);
    const auto *isi0_363 = buffer.data(isi0 + 363);

    const auto *isi1_301 = buffer.data(isi1 + 301);
    const auto *isi1_302 = buffer.data(isi1 + 302);
    const auto *isi1_303 = buffer.data(isi1 + 303);
    const auto *isi1_304 = buffer.data(isi1 + 304);
    const auto *isi1_305 = buffer.data(isi1 + 305);
    const auto *isi1_307 = buffer.data(isi1 + 307);
    const auto *isi1_313 = buffer.data(isi1 + 313);
    const auto *isi1_317 = buffer.data(isi1 + 317);
    const auto *isi1_322 = buffer.data(isi1 + 322);
    const auto *isi1_328 = buffer.data(isi1 + 328);
    const auto *isi1_331 = buffer.data(isi1 + 331);
    const auto *isi1_332 = buffer.data(isi1 + 332);
    const auto *isi1_333 = buffer.data(isi1 + 333);
    const auto *isi1_334 = buffer.data(isi1 + 334);
    const auto *isi1_335 = buffer.data(isi1 + 335);
    const auto *isi1_336 = buffer.data(isi1 + 336);
    const auto *isi1_339 = buffer.data(isi1 + 339);
    const auto *isi1_341 = buffer.data(isi1 + 341);
    const auto *isi1_342 = buffer.data(isi1 + 342);
    const auto *isi1_345 = buffer.data(isi1 + 345);
    const auto *isi1_346 = buffer.data(isi1 + 346);
    const auto *isi1_348 = buffer.data(isi1 + 348);
    const auto *isi1_350 = buffer.data(isi1 + 350);
    const auto *isi1_351 = buffer.data(isi1 + 351);
    const auto *isi1_353 = buffer.data(isi1 + 353);
    const auto *isi1_354 = buffer.data(isi1 + 354);
    const auto *isi1_356 = buffer.data(isi1 + 356);
    const auto *isi1_357 = buffer.data(isi1 + 357);
    const auto *isi1_359 = buffer.data(isi1 + 359);
    const auto *isi1_360 = buffer.data(isi1 + 360);
    const auto *isi1_361 = buffer.data(isi1 + 361);
    const auto *isi1_362 = buffer.data(isi1 + 362);
    const auto *isi1_363 = buffer.data(isi1 + 363);

    const auto *isk_381 = buffer.data(isk + 381);
    const auto *isk_388 = buffer.data(isk + 388);
    const auto *isk_389 = buffer.data(isk + 389);
    const auto *isk_390 = buffer.data(isk + 390);
    const auto *isk_391 = buffer.data(isk + 391);
    const auto *isk_392 = buffer.data(isk + 392);
    const auto *isk_393 = buffer.data(isk + 393);
    const auto *isk_394 = buffer.data(isk + 394);
    const auto *isk_395 = buffer.data(isk + 395);
    const auto *isk_396 = buffer.data(isk + 396);
    const auto *isk_398 = buffer.data(isk + 398);
    const auto *isk_399 = buffer.data(isk + 399);
    const auto *isk_401 = buffer.data(isk + 401);
    const auto *isk_402 = buffer.data(isk + 402);
    const auto *isk_405 = buffer.data(isk + 405);
    const auto *isk_406 = buffer.data(isk + 406);
    const auto *isk_410 = buffer.data(isk + 410);
    const auto *isk_411 = buffer.data(isk + 411);
    const auto *isk_416 = buffer.data(isk + 416);
    const auto *isk_423 = buffer.data(isk + 423);
    const auto *isk_424 = buffer.data(isk + 424);
    const auto *isk_425 = buffer.data(isk + 425);
    const auto *isk_426 = buffer.data(isk + 426);
    const auto *isk_427 = buffer.data(isk + 427);
    const auto *isk_428 = buffer.data(isk + 428);
    const auto *isk_429 = buffer.data(isk + 429);
    const auto *isk_430 = buffer.data(isk + 430);
    const auto *isk_431 = buffer.data(isk + 431);
    const auto *isk_432 = buffer.data(isk + 432);
    const auto *isk_434 = buffer.data(isk + 434);
    const auto *isk_435 = buffer.data(isk + 435);
    const auto *isk_437 = buffer.data(isk + 437);
    const auto *isk_438 = buffer.data(isk + 438);
    const auto *isk_441 = buffer.data(isk + 441);
    const auto *isk_442 = buffer.data(isk + 442);
    const auto *isk_444 = buffer.data(isk + 444);
    const auto *isk_446 = buffer.data(isk + 446);
    const auto *isk_447 = buffer.data(isk + 447);
    const auto *isk_449 = buffer.data(isk + 449);
    const auto *isk_450 = buffer.data(isk + 450);
    const auto *isk_452 = buffer.data(isk + 452);
    const auto *isk_453 = buffer.data(isk + 453);
    const auto *isk_455 = buffer.data(isk + 455);
    const auto *isk_456 = buffer.data(isk + 456);
    const auto *isk_457 = buffer.data(isk + 457);
    const auto *isk_459 = buffer.data(isk + 459);
    const auto *isk_460 = buffer.data(isk + 460);
    const auto *isk_461 = buffer.data(isk + 461);
    const auto *isk_462 = buffer.data(isk + 462);
    const auto *isk_463 = buffer.data(isk + 463);
    const auto *isk_464 = buffer.data(isk + 464);
    const auto *isk_465 = buffer.data(isk + 465);
    const auto *isk_466 = buffer.data(isk + 466);
    const auto *isk_467 = buffer.data(isk + 467);
    const auto *isk_468 = buffer.data(isk + 468);

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, pc_x, pc_z, hsk_388, hsk_390, \
                         hsk_391, hsk_392, isk_381, isk_388, isk_390, isk_391, \
                         isk_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_16 * hsk_388[k]
                   + f_3 * pc_x[k] * isk_388[k];

        t_479[k] = f_3 * pc_z[k] * isk_381[k];

        t_480[k] = f_16 * hsk_390[k]
                   + f_3 * pc_x[k] * isk_390[k];

        t_481[k] = f_16 * hsk_391[k]
                   + f_3 * pc_x[k] * isk_391[k];

        t_482[k] = f_16 * hsk_392[k]
                   + f_3 * pc_x[k] * isk_392[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pc_x, pc_y, hsk_244, hsk_393, hsk_394, \
                         hsk_395, isi0_301, isi1_301, isk_388, isk_393, isk_394, \
                         isk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_16 * hsk_393[k]
                   + f_3 * pc_x[k] * isk_393[k];

        t_484[k] = f_16 * hsk_394[k]
                   + f_3 * pc_x[k] * isk_394[k];

        t_485[k] = f_16 * hsk_395[k]
                   + f_3 * pc_x[k] * isk_395[k];

        t_486[k] = f_18 * hsk_244[k]
                   + f_1 * isi0_301[k]
                   - f_2 * isi1_301[k]
                   + f_3 * pc_y[k] * isk_388[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pc_z, isi0_301, isi0_302, isi0_303, \
                         isi1_301, isi1_302, isi1_303, isk_388, isk_389, isk_390, \
                         isk_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_3 * pc_z[k] * isk_388[k];

        t_488[k] = f_4 * isi0_301[k]
                   - f_5 * isi1_301[k]
                   + f_3 * pc_z[k] * isk_389[k];

        t_489[k] = f_6 * isi0_302[k]
                   - f_7 * isi1_302[k]
                   + f_3 * pc_z[k] * isk_390[k];

        t_490[k] = f_8 * isi0_303[k]
                   - f_9 * isi1_303[k]
                   + f_3 * pc_z[k] * isk_391[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_y, pc_z, hsk_251, isi0_304, isi0_305, \
                         isi0_307, isi1_304, isi1_305, isi1_307, isk_392, isk_393, \
                         isk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_10 * isi0_304[k]
                   - f_11 * isi1_304[k]
                   + f_3 * pc_z[k] * isk_392[k];

        t_492[k] = f_12 * isi0_305[k]
                   - f_13 * isi1_305[k]
                   + f_3 * pc_z[k] * isk_393[k];

        t_493[k] = f_18 * hsk_251[k]
                   + f_3 * pc_y[k] * isk_395[k];

        t_494[k] = f_1 * isi0_307[k]
                   - f_2 * isi1_307[k]
                   + f_3 * pc_z[k] * isk_395[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pa_z, pc_y, pc_z, hsl0_270, hsl0_273, \
                         hsk_216, hsk_252, hsl1_270, hsl1_273, \
                         isk_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * hsl0_270[k]
                   - f_14 * pc_z[k] * hsl1_270[k];

        t_496[k] = f_17 * hsk_252[k]
                   + f_3 * pc_y[k] * isk_396[k];

        t_497[k] = f_15 * hsk_216[k]
                   + f_3 * pc_z[k] * isk_396[k];

        t_498[k] = pa_z[k] * hsl0_273[k]
                   - f_14 * pc_z[k] * hsl1_273[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pa_z, pc_x, pc_y, pc_z, hsl0_276, hsk_254, \
                         hsk_401, hsl1_276, isi0_313, isi1_313, isk_398, \
                         isk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_17 * hsk_254[k]
                   + f_3 * pc_y[k] * isk_398[k];

        t_500[k] = f_16 * hsk_401[k]
                   + f_12 * isi0_313[k]
                   - f_13 * isi1_313[k]
                   + f_3 * pc_x[k] * isk_401[k];

        t_501[k] = pa_z[k] * hsl0_276[k]
                   - f_14 * pc_z[k] * hsl1_276[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, hsk_219, hsk_257, hsk_405, \
                         isi0_317, isi1_317, isk_399, isk_401, \
                         isk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_15 * hsk_219[k]
                   + f_3 * pc_z[k] * isk_399[k];

        t_503[k] = f_17 * hsk_257[k]
                   + f_3 * pc_y[k] * isk_401[k];

        t_504[k] = f_16 * hsk_405[k]
                   + f_10 * isi0_317[k]
                   - f_11 * isi1_317[k]
                   + f_3 * pc_x[k] * isk_405[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pa_z, pc_y, pc_z, hsl0_280, hsl0_282, \
                         hsk_222, hsk_223, hsk_261, hsl1_280, hsl1_282, isk_402, \
                         isk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pa_z[k] * hsl0_280[k]
                   - f_14 * pc_z[k] * hsl1_280[k];

        t_506[k] = f_15 * hsk_222[k]
                   + f_3 * pc_z[k] * isk_402[k];

        t_507[k] = pa_z[k] * hsl0_282[k]
                   + f_16 * hsk_223[k]
                   - f_14 * pc_z[k] * hsl1_282[k];

        t_508[k] = f_17 * hsk_261[k]
                   + f_3 * pc_y[k] * isk_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_z, pc_x, pc_z, hsl0_285, hsk_226, hsk_410, \
                         hsl1_285, isi0_322, isi1_322, isk_406, \
                         isk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_16 * hsk_410[k]
                   + f_8 * isi0_322[k]
                   - f_9 * isi1_322[k]
                   + f_3 * pc_x[k] * isk_410[k];

        t_510[k] = pa_z[k] * hsl0_285[k]
                   - f_14 * pc_z[k] * hsl1_285[k];

        t_511[k] = f_15 * hsk_226[k]
                   + f_3 * pc_z[k] * isk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pa_z, pc_y, pc_z, hsl0_287, hsl0_288, hsk_227, \
                         hsk_228, hsk_266, hsl1_287, hsl1_288, \
                         isk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pa_z[k] * hsl0_287[k]
                   + f_16 * hsk_227[k]
                   - f_14 * pc_z[k] * hsl1_287[k];

        t_513[k] = pa_z[k] * hsl0_288[k]
                   + f_17 * hsk_228[k]
                   - f_14 * pc_z[k] * hsl1_288[k];

        t_514[k] = f_17 * hsk_266[k]
                   + f_3 * pc_y[k] * isk_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pa_z, pc_x, pc_z, hsl0_291, hsk_231, hsk_416, \
                         hsl1_291, isi0_328, isi1_328, isk_411, \
                         isk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_16 * hsk_416[k]
                   + f_6 * isi0_328[k]
                   - f_7 * isi1_328[k]
                   + f_3 * pc_x[k] * isk_416[k];

        t_516[k] = pa_z[k] * hsl0_291[k]
                   - f_14 * pc_z[k] * hsl1_291[k];

        t_517[k] = f_15 * hsk_231[k]
                   + f_3 * pc_z[k] * isk_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pa_z, pc_z, hsl0_293, hsl0_294, hsl0_295, \
                         hsk_232, hsk_233, hsk_234, hsl1_293, hsl1_294, \
                         hsl1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pa_z[k] * hsl0_293[k]
                   + f_16 * hsk_232[k]
                   - f_14 * pc_z[k] * hsl1_293[k];

        t_519[k] = pa_z[k] * hsl0_294[k]
                   + f_17 * hsk_233[k]
                   - f_14 * pc_z[k] * hsl1_294[k];

        t_520[k] = pa_z[k] * hsl0_295[k]
                   + f_18 * hsk_234[k]
                   - f_14 * pc_z[k] * hsl1_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, hsk_272, hsk_423, hsk_424, \
                         hsk_425, isi0_335, isi1_335, isk_416, isk_423, isk_424, \
                         isk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * hsk_272[k]
                   + f_3 * pc_y[k] * isk_416[k];

        t_522[k] = f_16 * hsk_423[k]
                   + f_4 * isi0_335[k]
                   - f_5 * isi1_335[k]
                   + f_3 * pc_x[k] * isk_423[k];

        t_523[k] = f_16 * hsk_424[k]
                   + f_3 * pc_x[k] * isk_424[k];

        t_524[k] = f_16 * hsk_425[k]
                   + f_3 * pc_x[k] * isk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, hsk_426, hsk_427, hsk_428, \
                         hsk_429, hsk_430, isk_426, isk_427, isk_428, isk_429, \
                         isk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_16 * hsk_426[k]
                   + f_3 * pc_x[k] * isk_426[k];

        t_526[k] = f_16 * hsk_427[k]
                   + f_3 * pc_x[k] * isk_427[k];

        t_527[k] = f_16 * hsk_428[k]
                   + f_3 * pc_x[k] * isk_428[k];

        t_528[k] = f_16 * hsk_429[k]
                   + f_3 * pc_x[k] * isk_429[k];

        t_529[k] = f_16 * hsk_430[k]
                   + f_3 * pc_x[k] * isk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pa_z, pc_x, pc_z, hsl0_306, hsk_244, hsk_431, \
                         hsl1_306, isk_424, isk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_16 * hsk_431[k]
                   + f_3 * pc_x[k] * isk_431[k];

        t_531[k] = pa_z[k] * hsl0_306[k]
                   - f_14 * pc_z[k] * hsl1_306[k];

        t_532[k] = f_15 * hsk_244[k]
                   + f_3 * pc_z[k] * isk_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pc_y, hsk_282, hsk_283, hsk_284, isi0_331, \
                         isi0_332, isi0_333, isi1_331, isi1_332, isi1_333, isk_426, isk_427, \
                         isk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * hsk_282[k]
                   + f_12 * isi0_331[k]
                   - f_13 * isi1_331[k]
                   + f_3 * pc_y[k] * isk_426[k];

        t_534[k] = f_17 * hsk_283[k]
                   + f_10 * isi0_332[k]
                   - f_11 * isi1_332[k]
                   + f_3 * pc_y[k] * isk_427[k];

        t_535[k] = f_17 * hsk_284[k]
                   + f_8 * isi0_333[k]
                   - f_9 * isi1_333[k]
                   + f_3 * pc_y[k] * isk_428[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pc_y, hsk_285, hsk_286, hsk_287, isi0_334, \
                         isi0_335, isi1_334, isi1_335, isk_429, isk_430, \
                         isk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_17 * hsk_285[k]
                   + f_6 * isi0_334[k]
                   - f_7 * isi1_334[k]
                   + f_3 * pc_y[k] * isk_429[k];

        t_537[k] = f_17 * hsk_286[k]
                   + f_4 * isi0_335[k]
                   - f_5 * isi1_335[k]
                   + f_3 * pc_y[k] * isk_430[k];

        t_538[k] = f_17 * hsk_287[k]
                   + f_3 * pc_y[k] * isk_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pc_x, pc_y, pc_z, hsk_251, hsk_288, hsk_432, \
                         isi0_335, isi0_336, isi1_335, isi1_336, isk_431, \
                         isk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * hsk_251[k]
                   + f_1 * isi0_335[k]
                   - f_2 * isi1_335[k]
                   + f_3 * pc_z[k] * isk_431[k];

        t_540[k] = f_16 * hsk_432[k]
                   + f_1 * isi0_336[k]
                   - f_2 * isi1_336[k]
                   + f_3 * pc_x[k] * isk_432[k];

        t_541[k] = f_16 * hsk_288[k]
                   + f_3 * pc_y[k] * isk_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, pc_y, pc_z, hsk_252, hsk_290, hsk_435, \
                         isi0_339, isi1_339, isk_432, isk_434, \
                         isk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * hsk_252[k]
                   + f_3 * pc_z[k] * isk_432[k];

        t_543[k] = f_16 * hsk_435[k]
                   + f_12 * isi0_339[k]
                   - f_13 * isi1_339[k]
                   + f_3 * pc_x[k] * isk_435[k];

        t_544[k] = f_16 * hsk_290[k]
                   + f_3 * pc_y[k] * isk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, hsk_255, hsk_437, hsk_438, isi0_341, \
                         isi0_342, isi1_341, isi1_342, isk_435, isk_437, \
                         isk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_16 * hsk_437[k]
                   + f_12 * isi0_341[k]
                   - f_13 * isi1_341[k]
                   + f_3 * pc_x[k] * isk_437[k];

        t_546[k] = f_16 * hsk_438[k]
                   + f_10 * isi0_342[k]
                   - f_11 * isi1_342[k]
                   + f_3 * pc_x[k] * isk_438[k];

        t_547[k] = f_16 * hsk_255[k]
                   + f_3 * pc_z[k] * isk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pc_x, pc_y, hsk_293, hsk_441, hsk_442, isi0_345, \
                         isi0_346, isi1_345, isi1_346, isk_437, isk_441, \
                         isk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * hsk_293[k]
                   + f_3 * pc_y[k] * isk_437[k];

        t_549[k] = f_16 * hsk_441[k]
                   + f_10 * isi0_345[k]
                   - f_11 * isi1_345[k]
                   + f_3 * pc_x[k] * isk_441[k];

        t_550[k] = f_16 * hsk_442[k]
                   + f_8 * isi0_346[k]
                   - f_9 * isi1_346[k]
                   + f_3 * pc_x[k] * isk_442[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pc_x, pc_y, pc_z, hsk_258, hsk_297, hsk_444, \
                         isi0_348, isi1_348, isk_438, isk_441, \
                         isk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * hsk_258[k]
                   + f_3 * pc_z[k] * isk_438[k];

        t_552[k] = f_16 * hsk_444[k]
                   + f_8 * isi0_348[k]
                   - f_9 * isi1_348[k]
                   + f_3 * pc_x[k] * isk_444[k];

        t_553[k] = f_16 * hsk_297[k]
                   + f_3 * pc_y[k] * isk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_x, pc_z, hsk_262, hsk_446, hsk_447, isi0_350, \
                         isi0_351, isi1_350, isi1_351, isk_442, isk_446, \
                         isk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_16 * hsk_446[k]
                   + f_8 * isi0_350[k]
                   - f_9 * isi1_350[k]
                   + f_3 * pc_x[k] * isk_446[k];

        t_555[k] = f_16 * hsk_447[k]
                   + f_6 * isi0_351[k]
                   - f_7 * isi1_351[k]
                   + f_3 * pc_x[k] * isk_447[k];

        t_556[k] = f_16 * hsk_262[k]
                   + f_3 * pc_z[k] * isk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, pc_y, hsk_302, hsk_449, hsk_450, isi0_353, \
                         isi0_354, isi1_353, isi1_354, isk_446, isk_449, \
                         isk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_16 * hsk_449[k]
                   + f_6 * isi0_353[k]
                   - f_7 * isi1_353[k]
                   + f_3 * pc_x[k] * isk_449[k];

        t_558[k] = f_16 * hsk_450[k]
                   + f_6 * isi0_354[k]
                   - f_7 * isi1_354[k]
                   + f_3 * pc_x[k] * isk_450[k];

        t_559[k] = f_16 * hsk_302[k]
                   + f_3 * pc_y[k] * isk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pc_x, pc_z, hsk_267, hsk_452, hsk_453, isi0_356, \
                         isi0_357, isi1_356, isi1_357, isk_447, isk_452, \
                         isk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_16 * hsk_452[k]
                   + f_6 * isi0_356[k]
                   - f_7 * isi1_356[k]
                   + f_3 * pc_x[k] * isk_452[k];

        t_561[k] = f_16 * hsk_453[k]
                   + f_4 * isi0_357[k]
                   - f_5 * isi1_357[k]
                   + f_3 * pc_x[k] * isk_453[k];

        t_562[k] = f_16 * hsk_267[k]
                   + f_3 * pc_z[k] * isk_447[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, hsk_455, hsk_456, hsk_457, isi0_359, \
                         isi0_360, isi0_361, isi1_359, isi1_360, isi1_361, isk_455, isk_456, \
                         isk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_16 * hsk_455[k]
                   + f_4 * isi0_359[k]
                   - f_5 * isi1_359[k]
                   + f_3 * pc_x[k] * isk_455[k];

        t_564[k] = f_16 * hsk_456[k]
                   + f_4 * isi0_360[k]
                   - f_5 * isi1_360[k]
                   + f_3 * pc_x[k] * isk_456[k];

        t_565[k] = f_16 * hsk_457[k]
                   + f_4 * isi0_361[k]
                   - f_5 * isi1_361[k]
                   + f_3 * pc_x[k] * isk_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pc_x, pc_y, hsk_308, hsk_459, hsk_460, \
                         hsk_461, isi0_363, isi1_363, isk_452, isk_459, isk_460, \
                         isk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * hsk_308[k]
                   + f_3 * pc_y[k] * isk_452[k];

        t_567[k] = f_16 * hsk_459[k]
                   + f_4 * isi0_363[k]
                   - f_5 * isi1_363[k]
                   + f_3 * pc_x[k] * isk_459[k];

        t_568[k] = f_16 * hsk_460[k]
                   + f_3 * pc_x[k] * isk_460[k];

        t_569[k] = f_16 * hsk_461[k]
                   + f_3 * pc_x[k] * isk_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, hsk_462, hsk_463, hsk_464, \
                         hsk_465, hsk_466, isk_462, isk_463, isk_464, isk_465, \
                         isk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_16 * hsk_462[k]
                   + f_3 * pc_x[k] * isk_462[k];

        t_571[k] = f_16 * hsk_463[k]
                   + f_3 * pc_x[k] * isk_463[k];

        t_572[k] = f_16 * hsk_464[k]
                   + f_3 * pc_x[k] * isk_464[k];

        t_573[k] = f_16 * hsk_465[k]
                   + f_3 * pc_x[k] * isk_465[k];

        t_574[k] = f_16 * hsk_466[k]
                   + f_3 * pc_x[k] * isk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, pc_z, hsk_280, hsk_316, hsk_467, \
                         isi0_357, isi1_357, isk_460, isk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_16 * hsk_467[k]
                   + f_3 * pc_x[k] * isk_467[k];

        t_576[k] = f_16 * hsk_316[k]
                   + f_1 * isi0_357[k]
                   - f_2 * isi1_357[k]
                   + f_3 * pc_y[k] * isk_460[k];

        t_577[k] = f_16 * hsk_280[k]
                   + f_3 * pc_z[k] * isk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, hsk_318, hsk_319, hsk_320, isi0_359, \
                         isi0_360, isi0_361, isi1_359, isi1_360, isi1_361, isk_462, isk_463, \
                         isk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * hsk_318[k]
                   + f_12 * isi0_359[k]
                   - f_13 * isi1_359[k]
                   + f_3 * pc_y[k] * isk_462[k];

        t_579[k] = f_16 * hsk_319[k]
                   + f_10 * isi0_360[k]
                   - f_11 * isi1_360[k]
                   + f_3 * pc_y[k] * isk_463[k];

        t_580[k] = f_16 * hsk_320[k]
                   + f_8 * isi0_361[k]
                   - f_9 * isi1_361[k]
                   + f_3 * pc_y[k] * isk_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, hsk_321, hsk_322, hsk_323, isi0_362, \
                         isi0_363, isi1_362, isi1_363, isk_465, isk_466, \
                         isk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * hsk_321[k]
                   + f_6 * isi0_362[k]
                   - f_7 * isi1_362[k]
                   + f_3 * pc_y[k] * isk_465[k];

        t_582[k] = f_16 * hsk_322[k]
                   + f_4 * isi0_363[k]
                   - f_5 * isi1_363[k]
                   + f_3 * pc_y[k] * isk_466[k];

        t_583[k] = f_16 * hsk_323[k]
                   + f_3 * pc_y[k] * isk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_y, pc_y, pc_z, hsl0_405, hsk_287, \
                         hsk_288, hsk_324, hsl1_405, isi0_363, isi1_363, isk_467, \
                         isk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * hsk_287[k]
                   + f_1 * isi0_363[k]
                   - f_2 * isi1_363[k]
                   + f_3 * pc_z[k] * isk_467[k];

        t_585[k] = pa_y[k] * hsl0_405[k]
                   - f_14 * pc_y[k] * hsl1_405[k];

        t_586[k] = f_15 * hsk_324[k]
                   + f_3 * pc_y[k] * isk_468[k];

        t_587[k] = f_17 * hsk_288[k]
                   + f_3 * pc_z[k] * isk_468[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsl0_408 = buffer.data(hsl0 + 408);
    const auto *hsl0_410 = buffer.data(hsl0 + 410);
    const auto *hsl0_411 = buffer.data(hsl0 + 411);
    const auto *hsl0_414 = buffer.data(hsl0 + 414);
    const auto *hsl0_415 = buffer.data(hsl0 + 415);
    const auto *hsl0_417 = buffer.data(hsl0 + 417);
    const auto *hsl0_419 = buffer.data(hsl0 + 419);
    const auto *hsl0_420 = buffer.data(hsl0 + 420);
    const auto *hsl0_422 = buffer.data(hsl0 + 422);
    const auto *hsl0_423 = buffer.data(hsl0 + 423);
    const auto *hsl0_425 = buffer.data(hsl0 + 425);
    const auto *hsl0_426 = buffer.data(hsl0 + 426);
    const auto *hsl0_428 = buffer.data(hsl0 + 428);
    const auto *hsl0_429 = buffer.data(hsl0 + 429);
    const auto *hsl0_430 = buffer.data(hsl0 + 430);
    const auto *hsl0_432 = buffer.data(hsl0 + 432);
    const auto *hsl0_449 = buffer.data(hsl0 + 449);
    const auto *hsl0_675 = buffer.data(hsl0 + 675);
    const auto *hsl0_678 = buffer.data(hsl0 + 678);
    const auto *hsl0_681 = buffer.data(hsl0 + 681);
    const auto *hsl0_685 = buffer.data(hsl0 + 685);
    const auto *hsl0_690 = buffer.data(hsl0 + 690);
    const auto *hsl0_696 = buffer.data(hsl0 + 696);

    const auto *hsk_291 = buffer.data(hsk + 291);
    const auto *hsk_294 = buffer.data(hsk + 294);
    const auto *hsk_298 = buffer.data(hsk + 298);
    const auto *hsk_303 = buffer.data(hsk + 303);
    const auto *hsk_316 = buffer.data(hsk + 316);
    const auto *hsk_324 = buffer.data(hsk + 324);
    const auto *hsk_325 = buffer.data(hsk + 325);
    const auto *hsk_326 = buffer.data(hsk + 326);
    const auto *hsk_327 = buffer.data(hsk + 327);
    const auto *hsk_329 = buffer.data(hsk + 329);
    const auto *hsk_330 = buffer.data(hsk + 330);
    const auto *hsk_332 = buffer.data(hsk + 332);
    const auto *hsk_333 = buffer.data(hsk + 333);
    const auto *hsk_334 = buffer.data(hsk + 334);
    const auto *hsk_336 = buffer.data(hsk + 336);
    const auto *hsk_337 = buffer.data(hsk + 337);
    const auto *hsk_338 = buffer.data(hsk + 338);
    const auto *hsk_339 = buffer.data(hsk + 339);
    const auto *hsk_341 = buffer.data(hsk + 341);
    const auto *hsk_342 = buffer.data(hsk + 342);
    const auto *hsk_343 = buffer.data(hsk + 343);
    const auto *hsk_344 = buffer.data(hsk + 344);
    const auto *hsk_352 = buffer.data(hsk + 352);
    const auto *hsk_354 = buffer.data(hsk + 354);
    const auto *hsk_355 = buffer.data(hsk + 355);
    const auto *hsk_356 = buffer.data(hsk + 356);
    const auto *hsk_357 = buffer.data(hsk + 357);
    const auto *hsk_358 = buffer.data(hsk + 358);
    const auto *hsk_359 = buffer.data(hsk + 359);
    const auto *hsk_360 = buffer.data(hsk + 360);
    const auto *hsk_365 = buffer.data(hsk + 365);
    const auto *hsk_369 = buffer.data(hsk + 369);
    const auto *hsk_374 = buffer.data(hsk + 374);
    const auto *hsk_380 = buffer.data(hsk + 380);
    const auto *hsk_496 = buffer.data(hsk + 496);
    const auto *hsk_497 = buffer.data(hsk + 497);
    const auto *hsk_498 = buffer.data(hsk + 498);
    const auto *hsk_499 = buffer.data(hsk + 499);
    const auto *hsk_500 = buffer.data(hsk + 500);
    const auto *hsk_501 = buffer.data(hsk + 501);
    const auto *hsk_502 = buffer.data(hsk + 502);
    const auto *hsk_503 = buffer.data(hsk + 503);
    const auto *hsk_504 = buffer.data(hsk + 504);
    const auto *hsk_509 = buffer.data(hsk + 509);
    const auto *hsk_513 = buffer.data(hsk + 513);
    const auto *hsk_518 = buffer.data(hsk + 518);
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
    const auto *hsk_543 = buffer.data(hsk + 543);
    const auto *hsk_546 = buffer.data(hsk + 546);
    const auto *hsk_550 = buffer.data(hsk + 550);
    const auto *hsk_555 = buffer.data(hsk + 555);
    const auto *hsk_561 = buffer.data(hsk + 561);

    const auto *hsl1_408 = buffer.data(hsl1 + 408);
    const auto *hsl1_410 = buffer.data(hsl1 + 410);
    const auto *hsl1_411 = buffer.data(hsl1 + 411);
    const auto *hsl1_414 = buffer.data(hsl1 + 414);
    const auto *hsl1_415 = buffer.data(hsl1 + 415);
    const auto *hsl1_417 = buffer.data(hsl1 + 417);
    const auto *hsl1_419 = buffer.data(hsl1 + 419);
    const auto *hsl1_420 = buffer.data(hsl1 + 420);
    const auto *hsl1_422 = buffer.data(hsl1 + 422);
    const auto *hsl1_423 = buffer.data(hsl1 + 423);
    const auto *hsl1_425 = buffer.data(hsl1 + 425);
    const auto *hsl1_426 = buffer.data(hsl1 + 426);
    const auto *hsl1_428 = buffer.data(hsl1 + 428);
    const auto *hsl1_429 = buffer.data(hsl1 + 429);
    const auto *hsl1_430 = buffer.data(hsl1 + 430);
    const auto *hsl1_432 = buffer.data(hsl1 + 432);
    const auto *hsl1_449 = buffer.data(hsl1 + 449);
    const auto *hsl1_675 = buffer.data(hsl1 + 675);
    const auto *hsl1_678 = buffer.data(hsl1 + 678);
    const auto *hsl1_681 = buffer.data(hsl1 + 681);
    const auto *hsl1_685 = buffer.data(hsl1 + 685);
    const auto *hsl1_690 = buffer.data(hsl1 + 690);
    const auto *hsl1_696 = buffer.data(hsl1 + 696);

    const auto *isi0_385 = buffer.data(isi0 + 385);
    const auto *isi0_387 = buffer.data(isi0 + 387);
    const auto *isi0_388 = buffer.data(isi0 + 388);
    const auto *isi0_389 = buffer.data(isi0 + 389);
    const auto *isi0_390 = buffer.data(isi0 + 390);
    const auto *isi0_391 = buffer.data(isi0 + 391);
    const auto *isi0_392 = buffer.data(isi0 + 392);
    const auto *isi0_393 = buffer.data(isi0 + 393);
    const auto *isi0_394 = buffer.data(isi0 + 394);
    const auto *isi0_395 = buffer.data(isi0 + 395);
    const auto *isi0_396 = buffer.data(isi0 + 396);
    const auto *isi0_397 = buffer.data(isi0 + 397);
    const auto *isi0_398 = buffer.data(isi0 + 398);
    const auto *isi0_399 = buffer.data(isi0 + 399);
    const auto *isi0_400 = buffer.data(isi0 + 400);
    const auto *isi0_401 = buffer.data(isi0 + 401);
    const auto *isi0_402 = buffer.data(isi0 + 402);
    const auto *isi0_403 = buffer.data(isi0 + 403);
    const auto *isi0_404 = buffer.data(isi0 + 404);
    const auto *isi0_405 = buffer.data(isi0 + 405);
    const auto *isi0_406 = buffer.data(isi0 + 406);
    const auto *isi0_412 = buffer.data(isi0 + 412);
    const auto *isi0_413 = buffer.data(isi0 + 413);
    const auto *isi0_414 = buffer.data(isi0 + 414);
    const auto *isi0_415 = buffer.data(isi0 + 415);
    const auto *isi0_416 = buffer.data(isi0 + 416);
    const auto *isi0_417 = buffer.data(isi0 + 417);
    const auto *isi0_418 = buffer.data(isi0 + 418);
    const auto *isi0_419 = buffer.data(isi0 + 419);
    const auto *isi0_420 = buffer.data(isi0 + 420);
    const auto *isi0_422 = buffer.data(isi0 + 422);
    const auto *isi0_423 = buffer.data(isi0 + 423);
    const auto *isi0_425 = buffer.data(isi0 + 425);
    const auto *isi0_426 = buffer.data(isi0 + 426);
    const auto *isi0_427 = buffer.data(isi0 + 427);
    const auto *isi0_429 = buffer.data(isi0 + 429);
    const auto *isi0_430 = buffer.data(isi0 + 430);
    const auto *isi0_431 = buffer.data(isi0 + 431);
    const auto *isi0_432 = buffer.data(isi0 + 432);
    const auto *isi0_434 = buffer.data(isi0 + 434);

    const auto *isi1_385 = buffer.data(isi1 + 385);
    const auto *isi1_387 = buffer.data(isi1 + 387);
    const auto *isi1_388 = buffer.data(isi1 + 388);
    const auto *isi1_389 = buffer.data(isi1 + 389);
    const auto *isi1_390 = buffer.data(isi1 + 390);
    const auto *isi1_391 = buffer.data(isi1 + 391);
    const auto *isi1_392 = buffer.data(isi1 + 392);
    const auto *isi1_393 = buffer.data(isi1 + 393);
    const auto *isi1_394 = buffer.data(isi1 + 394);
    const auto *isi1_395 = buffer.data(isi1 + 395);
    const auto *isi1_396 = buffer.data(isi1 + 396);
    const auto *isi1_397 = buffer.data(isi1 + 397);
    const auto *isi1_398 = buffer.data(isi1 + 398);
    const auto *isi1_399 = buffer.data(isi1 + 399);
    const auto *isi1_400 = buffer.data(isi1 + 400);
    const auto *isi1_401 = buffer.data(isi1 + 401);
    const auto *isi1_402 = buffer.data(isi1 + 402);
    const auto *isi1_403 = buffer.data(isi1 + 403);
    const auto *isi1_404 = buffer.data(isi1 + 404);
    const auto *isi1_405 = buffer.data(isi1 + 405);
    const auto *isi1_406 = buffer.data(isi1 + 406);
    const auto *isi1_412 = buffer.data(isi1 + 412);
    const auto *isi1_413 = buffer.data(isi1 + 413);
    const auto *isi1_414 = buffer.data(isi1 + 414);
    const auto *isi1_415 = buffer.data(isi1 + 415);
    const auto *isi1_416 = buffer.data(isi1 + 416);
    const auto *isi1_417 = buffer.data(isi1 + 417);
    const auto *isi1_418 = buffer.data(isi1 + 418);
    const auto *isi1_419 = buffer.data(isi1 + 419);
    const auto *isi1_420 = buffer.data(isi1 + 420);
    const auto *isi1_422 = buffer.data(isi1 + 422);
    const auto *isi1_423 = buffer.data(isi1 + 423);
    const auto *isi1_425 = buffer.data(isi1 + 425);
    const auto *isi1_426 = buffer.data(isi1 + 426);
    const auto *isi1_427 = buffer.data(isi1 + 427);
    const auto *isi1_429 = buffer.data(isi1 + 429);
    const auto *isi1_430 = buffer.data(isi1 + 430);
    const auto *isi1_431 = buffer.data(isi1 + 431);
    const auto *isi1_432 = buffer.data(isi1 + 432);
    const auto *isi1_434 = buffer.data(isi1 + 434);

    const auto *isk_470 = buffer.data(isk + 470);
    const auto *isk_471 = buffer.data(isk + 471);
    const auto *isk_473 = buffer.data(isk + 473);
    const auto *isk_474 = buffer.data(isk + 474);
    const auto *isk_477 = buffer.data(isk + 477);
    const auto *isk_478 = buffer.data(isk + 478);
    const auto *isk_482 = buffer.data(isk + 482);
    const auto *isk_483 = buffer.data(isk + 483);
    const auto *isk_488 = buffer.data(isk + 488);
    const auto *isk_496 = buffer.data(isk + 496);
    const auto *isk_497 = buffer.data(isk + 497);
    const auto *isk_498 = buffer.data(isk + 498);
    const auto *isk_499 = buffer.data(isk + 499);
    const auto *isk_500 = buffer.data(isk + 500);
    const auto *isk_501 = buffer.data(isk + 501);
    const auto *isk_502 = buffer.data(isk + 502);
    const auto *isk_503 = buffer.data(isk + 503);
    const auto *isk_504 = buffer.data(isk + 504);
    const auto *isk_505 = buffer.data(isk + 505);
    const auto *isk_506 = buffer.data(isk + 506);
    const auto *isk_507 = buffer.data(isk + 507);
    const auto *isk_508 = buffer.data(isk + 508);
    const auto *isk_509 = buffer.data(isk + 509);
    const auto *isk_510 = buffer.data(isk + 510);
    const auto *isk_511 = buffer.data(isk + 511);
    const auto *isk_512 = buffer.data(isk + 512);
    const auto *isk_513 = buffer.data(isk + 513);
    const auto *isk_514 = buffer.data(isk + 514);
    const auto *isk_515 = buffer.data(isk + 515);
    const auto *isk_516 = buffer.data(isk + 516);
    const auto *isk_517 = buffer.data(isk + 517);
    const auto *isk_518 = buffer.data(isk + 518);
    const auto *isk_519 = buffer.data(isk + 519);
    const auto *isk_520 = buffer.data(isk + 520);
    const auto *isk_521 = buffer.data(isk + 521);
    const auto *isk_522 = buffer.data(isk + 522);
    const auto *isk_523 = buffer.data(isk + 523);
    const auto *isk_524 = buffer.data(isk + 524);
    const auto *isk_531 = buffer.data(isk + 531);
    const auto *isk_532 = buffer.data(isk + 532);
    const auto *isk_533 = buffer.data(isk + 533);
    const auto *isk_534 = buffer.data(isk + 534);
    const auto *isk_535 = buffer.data(isk + 535);
    const auto *isk_536 = buffer.data(isk + 536);
    const auto *isk_537 = buffer.data(isk + 537);
    const auto *isk_538 = buffer.data(isk + 538);
    const auto *isk_539 = buffer.data(isk + 539);
    const auto *isk_540 = buffer.data(isk + 540);
    const auto *isk_541 = buffer.data(isk + 541);
    const auto *isk_542 = buffer.data(isk + 542);
    const auto *isk_543 = buffer.data(isk + 543);
    const auto *isk_545 = buffer.data(isk + 545);
    const auto *isk_546 = buffer.data(isk + 546);
    const auto *isk_547 = buffer.data(isk + 547);
    const auto *isk_549 = buffer.data(isk + 549);
    const auto *isk_550 = buffer.data(isk + 550);
    const auto *isk_551 = buffer.data(isk + 551);
    const auto *isk_552 = buffer.data(isk + 552);
    const auto *isk_554 = buffer.data(isk + 554);
    const auto *isk_555 = buffer.data(isk + 555);
    const auto *isk_556 = buffer.data(isk + 556);
    const auto *isk_557 = buffer.data(isk + 557);
    const auto *isk_558 = buffer.data(isk + 558);
    const auto *isk_560 = buffer.data(isk + 560);

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_y, pc_y, hsl0_408, hsl0_410, hsl0_411, \
                         hsk_325, hsk_326, hsk_327, hsl1_408, hsl1_410, hsl1_411, \
                         isk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_y[k] * hsl0_408[k]
                   + f_16 * hsk_325[k]
                   - f_14 * pc_y[k] * hsl1_408[k];

        t_589[k] = f_15 * hsk_326[k]
                   + f_3 * pc_y[k] * isk_470[k];

        t_590[k] = pa_y[k] * hsl0_410[k]
                   - f_14 * pc_y[k] * hsl1_410[k];

        t_591[k] = pa_y[k] * hsl0_411[k]
                   + f_17 * hsk_327[k]
                   - f_14 * pc_y[k] * hsl1_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_y, pc_y, pc_z, hsl0_414, hsl0_415, \
                         hsk_291, hsk_329, hsk_330, hsl1_414, hsl1_415, isk_471, \
                         isk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * hsk_291[k]
                   + f_3 * pc_z[k] * isk_471[k];

        t_593[k] = f_15 * hsk_329[k]
                   + f_3 * pc_y[k] * isk_473[k];

        t_594[k] = pa_y[k] * hsl0_414[k]
                   - f_14 * pc_y[k] * hsl1_414[k];

        t_595[k] = pa_y[k] * hsl0_415[k]
                   + f_18 * hsk_330[k]
                   - f_14 * pc_y[k] * hsl1_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pa_y, pc_y, pc_z, hsl0_417, hsl0_419, \
                         hsk_294, hsk_332, hsk_333, hsl1_417, hsl1_419, isk_474, \
                         isk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * hsk_294[k]
                   + f_3 * pc_z[k] * isk_474[k];

        t_597[k] = pa_y[k] * hsl0_417[k]
                   + f_16 * hsk_332[k]
                   - f_14 * pc_y[k] * hsl1_417[k];

        t_598[k] = f_15 * hsk_333[k]
                   + f_3 * pc_y[k] * isk_477[k];

        t_599[k] = pa_y[k] * hsl0_419[k]
                   - f_14 * pc_y[k] * hsl1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pa_y, pc_y, pc_z, hsl0_420, hsl0_422, hsk_298, \
                         hsk_334, hsk_336, hsl1_420, hsl1_422, \
                         isk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pa_y[k] * hsl0_420[k]
                   + f_19 * hsk_334[k]
                   - f_14 * pc_y[k] * hsl1_420[k];

        t_601[k] = f_17 * hsk_298[k]
                   + f_3 * pc_z[k] * isk_478[k];

        t_602[k] = pa_y[k] * hsl0_422[k]
                   + f_17 * hsk_336[k]
                   - f_14 * pc_y[k] * hsl1_422[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pa_y, pc_y, hsl0_423, hsl0_425, hsl0_426, \
                         hsk_337, hsk_338, hsk_339, hsl1_423, hsl1_425, hsl1_426, \
                         isk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pa_y[k] * hsl0_423[k]
                   + f_16 * hsk_337[k]
                   - f_14 * pc_y[k] * hsl1_423[k];

        t_604[k] = f_15 * hsk_338[k]
                   + f_3 * pc_y[k] * isk_482[k];

        t_605[k] = pa_y[k] * hsl0_425[k]
                   - f_14 * pc_y[k] * hsl1_425[k];

        t_606[k] = pa_y[k] * hsl0_426[k]
                   + f_0 * hsk_339[k]
                   - f_14 * pc_y[k] * hsl1_426[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pa_y, pc_y, pc_z, hsl0_428, hsl0_429, hsk_303, \
                         hsk_341, hsk_342, hsl1_428, hsl1_429, \
                         isk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_17 * hsk_303[k]
                   + f_3 * pc_z[k] * isk_483[k];

        t_608[k] = pa_y[k] * hsl0_428[k]
                   + f_18 * hsk_341[k]
                   - f_14 * pc_y[k] * hsl1_428[k];

        t_609[k] = pa_y[k] * hsl0_429[k]
                   + f_17 * hsk_342[k]
                   - f_14 * pc_y[k] * hsl1_429[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pa_y, pc_x, pc_y, hsl0_430, hsl0_432, \
                         hsk_343, hsk_344, hsk_496, hsl1_430, hsl1_432, isk_488, \
                         isk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pa_y[k] * hsl0_430[k]
                   + f_16 * hsk_343[k]
                   - f_14 * pc_y[k] * hsl1_430[k];

        t_611[k] = f_15 * hsk_344[k]
                   + f_3 * pc_y[k] * isk_488[k];

        t_612[k] = pa_y[k] * hsl0_432[k]
                   - f_14 * pc_y[k] * hsl1_432[k];

        t_613[k] = f_16 * hsk_496[k]
                   + f_3 * pc_x[k] * isk_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, hsk_497, hsk_498, hsk_499, \
                         hsk_500, hsk_501, isk_497, isk_498, isk_499, isk_500, \
                         isk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_16 * hsk_497[k]
                   + f_3 * pc_x[k] * isk_497[k];

        t_615[k] = f_16 * hsk_498[k]
                   + f_3 * pc_x[k] * isk_498[k];

        t_616[k] = f_16 * hsk_499[k]
                   + f_3 * pc_x[k] * isk_499[k];

        t_617[k] = f_16 * hsk_500[k]
                   + f_3 * pc_x[k] * isk_500[k];

        t_618[k] = f_16 * hsk_501[k]
                   + f_3 * pc_x[k] * isk_501[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, hsk_316, hsk_352, \
                         hsk_502, hsk_503, isi0_385, isi1_385, isk_496, isk_502, \
                         isk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_16 * hsk_502[k]
                   + f_3 * pc_x[k] * isk_502[k];

        t_620[k] = f_16 * hsk_503[k]
                   + f_3 * pc_x[k] * isk_503[k];

        t_621[k] = f_15 * hsk_352[k]
                   + f_1 * isi0_385[k]
                   - f_2 * isi1_385[k]
                   + f_3 * pc_y[k] * isk_496[k];

        t_622[k] = f_17 * hsk_316[k]
                   + f_3 * pc_z[k] * isk_496[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_y, hsk_354, hsk_355, hsk_356, isi0_387, \
                         isi0_388, isi0_389, isi1_387, isi1_388, isi1_389, isk_498, isk_499, \
                         isk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_15 * hsk_354[k]
                   + f_12 * isi0_387[k]
                   - f_13 * isi1_387[k]
                   + f_3 * pc_y[k] * isk_498[k];

        t_624[k] = f_15 * hsk_355[k]
                   + f_10 * isi0_388[k]
                   - f_11 * isi1_388[k]
                   + f_3 * pc_y[k] * isk_499[k];

        t_625[k] = f_15 * hsk_356[k]
                   + f_8 * isi0_389[k]
                   - f_9 * isi1_389[k]
                   + f_3 * pc_y[k] * isk_500[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, hsk_357, hsk_358, hsk_359, isi0_390, \
                         isi0_391, isi1_390, isi1_391, isk_501, isk_502, \
                         isk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * hsk_357[k]
                   + f_6 * isi0_390[k]
                   - f_7 * isi1_390[k]
                   + f_3 * pc_y[k] * isk_501[k];

        t_627[k] = f_15 * hsk_358[k]
                   + f_4 * isi0_391[k]
                   - f_5 * isi1_391[k]
                   + f_3 * pc_y[k] * isk_502[k];

        t_628[k] = f_15 * hsk_359[k]
                   + f_3 * pc_y[k] * isk_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pa_y, pc_x, pc_y, pc_z, hsl0_449, \
                         hsk_324, hsk_504, hsl1_449, isi0_392, isi1_392, \
                         isk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pa_y[k] * hsl0_449[k]
                   - f_14 * pc_y[k] * hsl1_449[k];

        t_630[k] = f_16 * hsk_504[k]
                   + f_1 * isi0_392[k]
                   - f_2 * isi1_392[k]
                   + f_3 * pc_x[k] * isk_504[k];

        t_631[k] = f_3 * pc_y[k] * isk_504[k];

        t_632[k] = f_18 * hsk_324[k]
                   + f_3 * pc_z[k] * isk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, hsk_509, isi0_392, isi0_397, \
                         isi1_392, isi1_397, isk_505, isk_506, \
                         isk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_4 * isi0_392[k]
                   - f_5 * isi1_392[k]
                   + f_3 * pc_y[k] * isk_505[k];

        t_634[k] = f_3 * pc_y[k] * isk_506[k];

        t_635[k] = f_16 * hsk_509[k]
                   + f_12 * isi0_397[k]
                   - f_13 * isi1_397[k]
                   + f_3 * pc_x[k] * isk_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_y, isi0_393, isi0_394, isi1_393, isi1_394, \
                         isk_507, isk_508, isk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_6 * isi0_393[k]
                   - f_7 * isi1_393[k]
                   + f_3 * pc_y[k] * isk_507[k];

        t_637[k] = f_4 * isi0_394[k]
                   - f_5 * isi1_394[k]
                   + f_3 * pc_y[k] * isk_508[k];

        t_638[k] = f_3 * pc_y[k] * isk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_y, hsk_513, isi0_395, isi0_396, \
                         isi0_401, isi1_395, isi1_396, isi1_401, isk_510, isk_511, \
                         isk_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_16 * hsk_513[k]
                   + f_10 * isi0_401[k]
                   - f_11 * isi1_401[k]
                   + f_3 * pc_x[k] * isk_513[k];

        t_640[k] = f_8 * isi0_395[k]
                   - f_9 * isi1_395[k]
                   + f_3 * pc_y[k] * isk_510[k];

        t_641[k] = f_6 * isi0_396[k]
                   - f_7 * isi1_396[k]
                   + f_3 * pc_y[k] * isk_511[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, hsk_518, isi0_397, isi0_406, \
                         isi1_397, isi1_406, isk_512, isk_513, \
                         isk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_4 * isi0_397[k]
                   - f_5 * isi1_397[k]
                   + f_3 * pc_y[k] * isk_512[k];

        t_643[k] = f_3 * pc_y[k] * isk_513[k];

        t_644[k] = f_16 * hsk_518[k]
                   + f_8 * isi0_406[k]
                   - f_9 * isi1_406[k]
                   + f_3 * pc_x[k] * isk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_y, isi0_398, isi0_399, isi0_400, isi1_398, \
                         isi1_399, isi1_400, isk_514, isk_515, \
                         isk_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_10 * isi0_398[k]
                   - f_11 * isi1_398[k]
                   + f_3 * pc_y[k] * isk_514[k];

        t_646[k] = f_8 * isi0_399[k]
                   - f_9 * isi1_399[k]
                   + f_3 * pc_y[k] * isk_515[k];

        t_647[k] = f_6 * isi0_400[k]
                   - f_7 * isi1_400[k]
                   + f_3 * pc_y[k] * isk_516[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, pc_y, hsk_524, isi0_401, isi0_412, \
                         isi1_401, isi1_412, isk_517, isk_518, \
                         isk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_4 * isi0_401[k]
                   - f_5 * isi1_401[k]
                   + f_3 * pc_y[k] * isk_517[k];

        t_649[k] = f_3 * pc_y[k] * isk_518[k];

        t_650[k] = f_16 * hsk_524[k]
                   + f_6 * isi0_412[k]
                   - f_7 * isi1_412[k]
                   + f_3 * pc_x[k] * isk_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_y, isi0_402, isi0_403, isi0_404, isi1_402, \
                         isi1_403, isi1_404, isk_519, isk_520, \
                         isk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_12 * isi0_402[k]
                   - f_13 * isi1_402[k]
                   + f_3 * pc_y[k] * isk_519[k];

        t_652[k] = f_10 * isi0_403[k]
                   - f_11 * isi1_403[k]
                   + f_3 * pc_y[k] * isk_520[k];

        t_653[k] = f_8 * isi0_404[k]
                   - f_9 * isi1_404[k]
                   + f_3 * pc_y[k] * isk_521[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_y, isi0_405, isi0_406, isi1_405, isi1_406, \
                         isk_522, isk_523, isk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_6 * isi0_405[k]
                   - f_7 * isi1_405[k]
                   + f_3 * pc_y[k] * isk_522[k];

        t_655[k] = f_4 * isi0_406[k]
                   - f_5 * isi1_406[k]
                   + f_3 * pc_y[k] * isk_523[k];

        t_656[k] = f_3 * pc_y[k] * isk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pc_x, hsk_531, hsk_532, hsk_533, hsk_534, \
                         isi0_419, isi1_419, isk_531, isk_532, isk_533, \
                         isk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_16 * hsk_531[k]
                   + f_4 * isi0_419[k]
                   - f_5 * isi1_419[k]
                   + f_3 * pc_x[k] * isk_531[k];

        t_658[k] = f_16 * hsk_532[k]
                   + f_3 * pc_x[k] * isk_532[k];

        t_659[k] = f_16 * hsk_533[k]
                   + f_3 * pc_x[k] * isk_533[k];

        t_660[k] = f_16 * hsk_534[k]
                   + f_3 * pc_x[k] * isk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, pc_y, hsk_535, hsk_536, \
                         hsk_537, hsk_539, isk_531, isk_535, isk_536, isk_537, \
                         isk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_16 * hsk_535[k]
                   + f_3 * pc_x[k] * isk_535[k];

        t_662[k] = f_16 * hsk_536[k]
                   + f_3 * pc_x[k] * isk_536[k];

        t_663[k] = f_16 * hsk_537[k]
                   + f_3 * pc_x[k] * isk_537[k];

        t_664[k] = f_3 * pc_y[k] * isk_531[k];

        t_665[k] = f_16 * hsk_539[k]
                   + f_3 * pc_x[k] * isk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, pc_y, isi0_413, isi0_414, isi0_415, isi1_413, \
                         isi1_414, isi1_415, isk_532, isk_533, \
                         isk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * isi0_413[k]
                   - f_2 * isi1_413[k]
                   + f_3 * pc_y[k] * isk_532[k];

        t_667[k] = f_20 * isi0_414[k]
                   - f_21 * isi1_414[k]
                   + f_3 * pc_y[k] * isk_533[k];

        t_668[k] = f_12 * isi0_415[k]
                   - f_13 * isi1_415[k]
                   + f_3 * pc_y[k] * isk_534[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pc_y, isi0_416, isi0_417, isi0_418, isi1_416, \
                         isi1_417, isi1_418, isk_535, isk_536, \
                         isk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_10 * isi0_416[k]
                   - f_11 * isi1_416[k]
                   + f_3 * pc_y[k] * isk_535[k];

        t_670[k] = f_8 * isi0_417[k]
                   - f_9 * isi1_417[k]
                   + f_3 * pc_y[k] * isk_536[k];

        t_671[k] = f_6 * isi0_418[k]
                   - f_7 * isi1_418[k]
                   + f_3 * pc_y[k] * isk_537[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pa_x, pc_x, pc_y, pc_z, hsl0_675, \
                         hsk_359, hsk_540, hsl1_675, isi0_419, isi1_419, isk_538, \
                         isk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_4 * isi0_419[k]
                   - f_5 * isi1_419[k]
                   + f_3 * pc_y[k] * isk_538[k];

        t_673[k] = f_3 * pc_y[k] * isk_539[k];

        t_674[k] = f_18 * hsk_359[k]
                   + f_1 * isi0_419[k]
                   - f_2 * isi1_419[k]
                   + f_3 * pc_z[k] * isk_539[k];

        t_675[k] = pa_x[k] * hsl0_675[k]
                   + f_22 * hsk_540[k]
                   - f_14 * pc_x[k] * hsl1_675[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pa_x, pc_x, pc_y, pc_z, hsl0_678, \
                         hsk_360, hsk_543, hsl1_678, isk_540, isk_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_19 * hsk_360[k]
                   + f_3 * pc_y[k] * isk_540[k];

        t_677[k] = f_3 * pc_z[k] * isk_540[k];

        t_678[k] = pa_x[k] * hsl0_678[k]
                   + f_0 * hsk_543[k]
                   - f_14 * pc_x[k] * hsl1_678[k];

        t_679[k] = f_3 * pc_z[k] * isk_541[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pa_x, pc_x, pc_z, hsl0_681, hsk_546, hsl1_681, \
                         isi0_420, isi1_420, isk_542, isk_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * isi0_420[k]
                   - f_5 * isi1_420[k]
                   + f_3 * pc_z[k] * isk_542[k];

        t_681[k] = pa_x[k] * hsl0_681[k]
                   + f_19 * hsk_546[k]
                   - f_14 * pc_x[k] * hsl1_681[k];

        t_682[k] = f_3 * pc_z[k] * isk_543[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, pa_x, pc_x, pc_y, pc_z, hsl0_685, \
                         hsk_365, hsk_550, hsl1_685, isi0_422, isi1_422, isk_545, \
                         isk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_19 * hsk_365[k]
                   + f_3 * pc_y[k] * isk_545[k];

        t_684[k] = f_6 * isi0_422[k]
                   - f_7 * isi1_422[k]
                   + f_3 * pc_z[k] * isk_545[k];

        t_685[k] = pa_x[k] * hsl0_685[k]
                   + f_18 * hsk_550[k]
                   - f_14 * pc_x[k] * hsl1_685[k];

        t_686[k] = f_3 * pc_z[k] * isk_546[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pc_y, pc_z, hsk_369, isi0_423, isi0_425, \
                         isi1_423, isi1_425, isk_547, isk_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_4 * isi0_423[k]
                   - f_5 * isi1_423[k]
                   + f_3 * pc_z[k] * isk_547[k];

        t_688[k] = f_19 * hsk_369[k]
                   + f_3 * pc_y[k] * isk_549[k];

        t_689[k] = f_8 * isi0_425[k]
                   - f_9 * isi1_425[k]
                   + f_3 * pc_z[k] * isk_549[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, pa_x, pc_x, pc_z, hsl0_690, hsk_555, hsl1_690, \
                         isi0_426, isi1_426, isk_550, isk_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pa_x[k] * hsl0_690[k]
                   + f_17 * hsk_555[k]
                   - f_14 * pc_x[k] * hsl1_690[k];

        t_691[k] = f_3 * pc_z[k] * isk_550[k];

        t_692[k] = f_4 * isi0_426[k]
                   - f_5 * isi1_426[k]
                   + f_3 * pc_z[k] * isk_551[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pc_y, pc_z, hsk_374, isi0_427, isi0_429, \
                         isi1_427, isi1_429, isk_552, isk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_6 * isi0_427[k]
                   - f_7 * isi1_427[k]
                   + f_3 * pc_z[k] * isk_552[k];

        t_694[k] = f_19 * hsk_374[k]
                   + f_3 * pc_y[k] * isk_554[k];

        t_695[k] = f_10 * isi0_429[k]
                   - f_11 * isi1_429[k]
                   + f_3 * pc_z[k] * isk_554[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pa_x, pc_x, pc_z, hsl0_696, hsk_561, hsl1_696, \
                         isi0_430, isi1_430, isk_555, isk_556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = pa_x[k] * hsl0_696[k]
                   + f_16 * hsk_561[k]
                   - f_14 * pc_x[k] * hsl1_696[k];

        t_697[k] = f_3 * pc_z[k] * isk_555[k];

        t_698[k] = f_4 * isi0_430[k]
                   - f_5 * isi1_430[k]
                   + f_3 * pc_z[k] * isk_556[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pc_y, pc_z, hsk_380, isi0_431, isi0_432, \
                         isi0_434, isi1_431, isi1_432, isi1_434, isk_557, isk_558, \
                         isk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_6 * isi0_431[k]
                   - f_7 * isi1_431[k]
                   + f_3 * pc_z[k] * isk_557[k];

        t_700[k] = f_8 * isi0_432[k]
                   - f_9 * isi1_432[k]
                   + f_3 * pc_z[k] * isk_558[k];

        t_701[k] = f_19 * hsk_380[k]
                   + f_3 * pc_y[k] * isk_560[k];

        t_702[k] = f_12 * isi0_434[k]
                   - f_13 * isi1_434[k]
                   + f_3 * pc_z[k] * isk_560[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_22 = 4.0 / q;

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
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsl0_450 = buffer.data(hsl0 + 450);
    const auto *hsl0_453 = buffer.data(hsl0 + 453);
    const auto *hsl0_456 = buffer.data(hsl0 + 456);
    const auto *hsl0_460 = buffer.data(hsl0 + 460);
    const auto *hsl0_465 = buffer.data(hsl0 + 465);
    const auto *hsl0_471 = buffer.data(hsl0 + 471);
    const auto *hsl0_711 = buffer.data(hsl0 + 711);
    const auto *hsl0_713 = buffer.data(hsl0 + 713);
    const auto *hsl0_714 = buffer.data(hsl0 + 714);
    const auto *hsl0_715 = buffer.data(hsl0 + 715);
    const auto *hsl0_716 = buffer.data(hsl0 + 716);
    const auto *hsl0_717 = buffer.data(hsl0 + 717);
    const auto *hsl0_719 = buffer.data(hsl0 + 719);
    const auto *hsl0_725 = buffer.data(hsl0 + 725);
    const auto *hsl0_729 = buffer.data(hsl0 + 729);
    const auto *hsl0_732 = buffer.data(hsl0 + 732);
    const auto *hsl0_734 = buffer.data(hsl0 + 734);
    const auto *hsl0_737 = buffer.data(hsl0 + 737);
    const auto *hsl0_738 = buffer.data(hsl0 + 738);
    const auto *hsl0_740 = buffer.data(hsl0 + 740);
    const auto *hsl0_743 = buffer.data(hsl0 + 743);
    const auto *hsl0_744 = buffer.data(hsl0 + 744);
    const auto *hsl0_745 = buffer.data(hsl0 + 745);
    const auto *hsl0_747 = buffer.data(hsl0 + 747);
    const auto *hsl0_756 = buffer.data(hsl0 + 756);
    const auto *hsl0_758 = buffer.data(hsl0 + 758);
    const auto *hsl0_759 = buffer.data(hsl0 + 759);
    const auto *hsl0_760 = buffer.data(hsl0 + 760);
    const auto *hsl0_761 = buffer.data(hsl0 + 761);
    const auto *hsl0_762 = buffer.data(hsl0 + 762);
    const auto *hsl0_764 = buffer.data(hsl0 + 764);
    const auto *hsl0_765 = buffer.data(hsl0 + 765);
    const auto *hsl0_768 = buffer.data(hsl0 + 768);
    const auto *hsl0_770 = buffer.data(hsl0 + 770);
    const auto *hsl0_771 = buffer.data(hsl0 + 771);
    const auto *hsl0_774 = buffer.data(hsl0 + 774);
    const auto *hsl0_775 = buffer.data(hsl0 + 775);
    const auto *hsl0_777 = buffer.data(hsl0 + 777);
    const auto *hsl0_779 = buffer.data(hsl0 + 779);
    const auto *hsl0_780 = buffer.data(hsl0 + 780);
    const auto *hsl0_782 = buffer.data(hsl0 + 782);
    const auto *hsl0_783 = buffer.data(hsl0 + 783);
    const auto *hsl0_785 = buffer.data(hsl0 + 785);
    const auto *hsl0_786 = buffer.data(hsl0 + 786);
    const auto *hsl0_788 = buffer.data(hsl0 + 788);
    const auto *hsl0_789 = buffer.data(hsl0 + 789);
    const auto *hsl0_790 = buffer.data(hsl0 + 790);
    const auto *hsl0_792 = buffer.data(hsl0 + 792);
    const auto *hsl0_801 = buffer.data(hsl0 + 801);
    const auto *hsl0_803 = buffer.data(hsl0 + 803);
    const auto *hsl0_804 = buffer.data(hsl0 + 804);
    const auto *hsl0_805 = buffer.data(hsl0 + 805);
    const auto *hsl0_806 = buffer.data(hsl0 + 806);
    const auto *hsl0_807 = buffer.data(hsl0 + 807);
    const auto *hsl0_809 = buffer.data(hsl0 + 809);
    const auto *hsl0_810 = buffer.data(hsl0 + 810);
    const auto *hsl0_813 = buffer.data(hsl0 + 813);
    const auto *hsl0_815 = buffer.data(hsl0 + 815);
    const auto *hsl0_816 = buffer.data(hsl0 + 816);
    const auto *hsl0_819 = buffer.data(hsl0 + 819);
    const auto *hsl0_820 = buffer.data(hsl0 + 820);

    const auto *hsk_360 = buffer.data(hsk + 360);
    const auto *hsk_363 = buffer.data(hsk + 363);
    const auto *hsk_366 = buffer.data(hsk + 366);
    const auto *hsk_370 = buffer.data(hsk + 370);
    const auto *hsk_375 = buffer.data(hsk + 375);
    const auto *hsk_388 = buffer.data(hsk + 388);
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
    const auto *hsk_431 = buffer.data(hsk + 431);
    const auto *hsk_432 = buffer.data(hsk + 432);
    const auto *hsk_434 = buffer.data(hsk + 434);
    const auto *hsk_435 = buffer.data(hsk + 435);
    const auto *hsk_437 = buffer.data(hsk + 437);
    const auto *hsk_438 = buffer.data(hsk + 438);
    const auto *hsk_441 = buffer.data(hsk + 441);
    const auto *hsk_446 = buffer.data(hsk + 446);
    const auto *hsk_452 = buffer.data(hsk + 452);
    const auto *hsk_467 = buffer.data(hsk + 467);
    const auto *hsk_468 = buffer.data(hsk + 468);
    const auto *hsk_470 = buffer.data(hsk + 470);
    const auto *hsk_473 = buffer.data(hsk + 473);
    const auto *hsk_568 = buffer.data(hsk + 568);
    const auto *hsk_570 = buffer.data(hsk + 570);
    const auto *hsk_571 = buffer.data(hsk + 571);
    const auto *hsk_572 = buffer.data(hsk + 572);
    const auto *hsk_573 = buffer.data(hsk + 573);
    const auto *hsk_574 = buffer.data(hsk + 574);
    const auto *hsk_575 = buffer.data(hsk + 575);
    const auto *hsk_581 = buffer.data(hsk + 581);
    const auto *hsk_585 = buffer.data(hsk + 585);
    const auto *hsk_588 = buffer.data(hsk + 588);
    const auto *hsk_590 = buffer.data(hsk + 590);
    const auto *hsk_593 = buffer.data(hsk + 593);
    const auto *hsk_594 = buffer.data(hsk + 594);
    const auto *hsk_596 = buffer.data(hsk + 596);
    const auto *hsk_599 = buffer.data(hsk + 599);
    const auto *hsk_600 = buffer.data(hsk + 600);
    const auto *hsk_601 = buffer.data(hsk + 601);
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
    const auto *hsk_615 = buffer.data(hsk + 615);
    const auto *hsk_617 = buffer.data(hsk + 617);
    const auto *hsk_618 = buffer.data(hsk + 618);
    const auto *hsk_621 = buffer.data(hsk + 621);
    const auto *hsk_622 = buffer.data(hsk + 622);
    const auto *hsk_624 = buffer.data(hsk + 624);
    const auto *hsk_626 = buffer.data(hsk + 626);
    const auto *hsk_627 = buffer.data(hsk + 627);
    const auto *hsk_629 = buffer.data(hsk + 629);
    const auto *hsk_630 = buffer.data(hsk + 630);
    const auto *hsk_632 = buffer.data(hsk + 632);
    const auto *hsk_633 = buffer.data(hsk + 633);
    const auto *hsk_635 = buffer.data(hsk + 635);
    const auto *hsk_636 = buffer.data(hsk + 636);
    const auto *hsk_637 = buffer.data(hsk + 637);
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
    const auto *hsk_651 = buffer.data(hsk + 651);
    const auto *hsk_653 = buffer.data(hsk + 653);
    const auto *hsk_654 = buffer.data(hsk + 654);
    const auto *hsk_657 = buffer.data(hsk + 657);
    const auto *hsk_658 = buffer.data(hsk + 658);

    const auto *hsl1_450 = buffer.data(hsl1 + 450);
    const auto *hsl1_453 = buffer.data(hsl1 + 453);
    const auto *hsl1_456 = buffer.data(hsl1 + 456);
    const auto *hsl1_460 = buffer.data(hsl1 + 460);
    const auto *hsl1_465 = buffer.data(hsl1 + 465);
    const auto *hsl1_471 = buffer.data(hsl1 + 471);
    const auto *hsl1_711 = buffer.data(hsl1 + 711);
    const auto *hsl1_713 = buffer.data(hsl1 + 713);
    const auto *hsl1_714 = buffer.data(hsl1 + 714);
    const auto *hsl1_715 = buffer.data(hsl1 + 715);
    const auto *hsl1_716 = buffer.data(hsl1 + 716);
    const auto *hsl1_717 = buffer.data(hsl1 + 717);
    const auto *hsl1_719 = buffer.data(hsl1 + 719);
    const auto *hsl1_725 = buffer.data(hsl1 + 725);
    const auto *hsl1_729 = buffer.data(hsl1 + 729);
    const auto *hsl1_732 = buffer.data(hsl1 + 732);
    const auto *hsl1_734 = buffer.data(hsl1 + 734);
    const auto *hsl1_737 = buffer.data(hsl1 + 737);
    const auto *hsl1_738 = buffer.data(hsl1 + 738);
    const auto *hsl1_740 = buffer.data(hsl1 + 740);
    const auto *hsl1_743 = buffer.data(hsl1 + 743);
    const auto *hsl1_744 = buffer.data(hsl1 + 744);
    const auto *hsl1_745 = buffer.data(hsl1 + 745);
    const auto *hsl1_747 = buffer.data(hsl1 + 747);
    const auto *hsl1_756 = buffer.data(hsl1 + 756);
    const auto *hsl1_758 = buffer.data(hsl1 + 758);
    const auto *hsl1_759 = buffer.data(hsl1 + 759);
    const auto *hsl1_760 = buffer.data(hsl1 + 760);
    const auto *hsl1_761 = buffer.data(hsl1 + 761);
    const auto *hsl1_762 = buffer.data(hsl1 + 762);
    const auto *hsl1_764 = buffer.data(hsl1 + 764);
    const auto *hsl1_765 = buffer.data(hsl1 + 765);
    const auto *hsl1_768 = buffer.data(hsl1 + 768);
    const auto *hsl1_770 = buffer.data(hsl1 + 770);
    const auto *hsl1_771 = buffer.data(hsl1 + 771);
    const auto *hsl1_774 = buffer.data(hsl1 + 774);
    const auto *hsl1_775 = buffer.data(hsl1 + 775);
    const auto *hsl1_777 = buffer.data(hsl1 + 777);
    const auto *hsl1_779 = buffer.data(hsl1 + 779);
    const auto *hsl1_780 = buffer.data(hsl1 + 780);
    const auto *hsl1_782 = buffer.data(hsl1 + 782);
    const auto *hsl1_783 = buffer.data(hsl1 + 783);
    const auto *hsl1_785 = buffer.data(hsl1 + 785);
    const auto *hsl1_786 = buffer.data(hsl1 + 786);
    const auto *hsl1_788 = buffer.data(hsl1 + 788);
    const auto *hsl1_789 = buffer.data(hsl1 + 789);
    const auto *hsl1_790 = buffer.data(hsl1 + 790);
    const auto *hsl1_792 = buffer.data(hsl1 + 792);
    const auto *hsl1_801 = buffer.data(hsl1 + 801);
    const auto *hsl1_803 = buffer.data(hsl1 + 803);
    const auto *hsl1_804 = buffer.data(hsl1 + 804);
    const auto *hsl1_805 = buffer.data(hsl1 + 805);
    const auto *hsl1_806 = buffer.data(hsl1 + 806);
    const auto *hsl1_807 = buffer.data(hsl1 + 807);
    const auto *hsl1_809 = buffer.data(hsl1 + 809);
    const auto *hsl1_810 = buffer.data(hsl1 + 810);
    const auto *hsl1_813 = buffer.data(hsl1 + 813);
    const auto *hsl1_815 = buffer.data(hsl1 + 815);
    const auto *hsl1_816 = buffer.data(hsl1 + 816);
    const auto *hsl1_819 = buffer.data(hsl1 + 819);
    const auto *hsl1_820 = buffer.data(hsl1 + 820);

    const auto *isk_561 = buffer.data(isk + 561);
    const auto *isk_568 = buffer.data(isk + 568);
    const auto *isk_570 = buffer.data(isk + 570);
    const auto *isk_571 = buffer.data(isk + 571);
    const auto *isk_572 = buffer.data(isk + 572);
    const auto *isk_573 = buffer.data(isk + 573);
    const auto *isk_574 = buffer.data(isk + 574);
    const auto *isk_575 = buffer.data(isk + 575);
    const auto *isk_576 = buffer.data(isk + 576);
    const auto *isk_578 = buffer.data(isk + 578);
    const auto *isk_579 = buffer.data(isk + 579);
    const auto *isk_581 = buffer.data(isk + 581);
    const auto *isk_582 = buffer.data(isk + 582);
    const auto *isk_585 = buffer.data(isk + 585);
    const auto *isk_586 = buffer.data(isk + 586);
    const auto *isk_590 = buffer.data(isk + 590);
    const auto *isk_591 = buffer.data(isk + 591);
    const auto *isk_596 = buffer.data(isk + 596);
    const auto *isk_604 = buffer.data(isk + 604);
    const auto *isk_605 = buffer.data(isk + 605);
    const auto *isk_606 = buffer.data(isk + 606);
    const auto *isk_607 = buffer.data(isk + 607);
    const auto *isk_608 = buffer.data(isk + 608);
    const auto *isk_609 = buffer.data(isk + 609);
    const auto *isk_610 = buffer.data(isk + 610);
    const auto *isk_611 = buffer.data(isk + 611);
    const auto *isk_612 = buffer.data(isk + 612);
    const auto *isk_614 = buffer.data(isk + 614);
    const auto *isk_615 = buffer.data(isk + 615);
    const auto *isk_617 = buffer.data(isk + 617);
    const auto *isk_618 = buffer.data(isk + 618);
    const auto *isk_621 = buffer.data(isk + 621);
    const auto *isk_622 = buffer.data(isk + 622);
    const auto *isk_626 = buffer.data(isk + 626);
    const auto *isk_627 = buffer.data(isk + 627);
    const auto *isk_632 = buffer.data(isk + 632);
    const auto *isk_640 = buffer.data(isk + 640);
    const auto *isk_641 = buffer.data(isk + 641);
    const auto *isk_642 = buffer.data(isk + 642);
    const auto *isk_643 = buffer.data(isk + 643);
    const auto *isk_644 = buffer.data(isk + 644);
    const auto *isk_645 = buffer.data(isk + 645);
    const auto *isk_646 = buffer.data(isk + 646);
    const auto *isk_647 = buffer.data(isk + 647);
    const auto *isk_648 = buffer.data(isk + 648);
    const auto *isk_650 = buffer.data(isk + 650);
    const auto *isk_651 = buffer.data(isk + 651);
    const auto *isk_653 = buffer.data(isk + 653);
    const auto *isk_654 = buffer.data(isk + 654);

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, pc_x, pc_z, hsk_568, hsk_570, \
                         hsk_571, hsk_572, isk_561, isk_568, isk_570, isk_571, \
                         isk_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_15 * hsk_568[k]
                   + f_3 * pc_x[k] * isk_568[k];

        t_704[k] = f_3 * pc_z[k] * isk_561[k];

        t_705[k] = f_15 * hsk_570[k]
                   + f_3 * pc_x[k] * isk_570[k];

        t_706[k] = f_15 * hsk_571[k]
                   + f_3 * pc_x[k] * isk_571[k];

        t_707[k] = f_15 * hsk_572[k]
                   + f_3 * pc_x[k] * isk_572[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, pa_x, pc_x, hsl0_711, hsk_573, hsk_574, \
                         hsk_575, hsl1_711, isk_573, isk_574, isk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_15 * hsk_573[k]
                   + f_3 * pc_x[k] * isk_573[k];

        t_709[k] = f_15 * hsk_574[k]
                   + f_3 * pc_x[k] * isk_574[k];

        t_710[k] = f_15 * hsk_575[k]
                   + f_3 * pc_x[k] * isk_575[k];

        t_711[k] = pa_x[k] * hsl0_711[k]
                   - f_14 * pc_x[k] * hsl1_711[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_x, pc_x, pc_z, hsl0_713, hsl0_714, \
                         hsl0_715, hsl1_713, hsl1_714, hsl1_715, \
                         isk_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_3 * pc_z[k] * isk_568[k];

        t_713[k] = pa_x[k] * hsl0_713[k]
                   - f_14 * pc_x[k] * hsl1_713[k];

        t_714[k] = pa_x[k] * hsl0_714[k]
                   - f_14 * pc_x[k] * hsl1_714[k];

        t_715[k] = pa_x[k] * hsl0_715[k]
                   - f_14 * pc_x[k] * hsl1_715[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pa_x, pc_x, pc_y, hsl0_716, hsl0_717, \
                         hsl0_719, hsk_395, hsl1_716, hsl1_717, hsl1_719, \
                         isk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = pa_x[k] * hsl0_716[k]
                   - f_14 * pc_x[k] * hsl1_716[k];

        t_717[k] = pa_x[k] * hsl0_717[k]
                   - f_14 * pc_x[k] * hsl1_717[k];

        t_718[k] = f_19 * hsk_395[k]
                   + f_3 * pc_y[k] * isk_575[k];

        t_719[k] = pa_x[k] * hsl0_719[k]
                   - f_14 * pc_x[k] * hsl1_719[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_z, pc_y, pc_z, hsl0_450, hsl0_453, \
                         hsk_360, hsk_396, hsl1_450, hsl1_453, \
                         isk_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_z[k] * hsl0_450[k]
                   - f_14 * pc_z[k] * hsl1_450[k];

        t_721[k] = f_18 * hsk_396[k]
                   + f_3 * pc_y[k] * isk_576[k];

        t_722[k] = f_15 * hsk_360[k]
                   + f_3 * pc_z[k] * isk_576[k];

        t_723[k] = pa_z[k] * hsl0_453[k]
                   - f_14 * pc_z[k] * hsl1_453[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pa_x, pa_z, pc_x, pc_y, pc_z, hsl0_456, \
                         hsl0_725, hsk_398, hsk_581, hsl1_456, hsl1_725, \
                         isk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_18 * hsk_398[k]
                   + f_3 * pc_y[k] * isk_578[k];

        t_725[k] = pa_x[k] * hsl0_725[k]
                   + f_0 * hsk_581[k]
                   - f_14 * pc_x[k] * hsl1_725[k];

        t_726[k] = pa_z[k] * hsl0_456[k]
                   - f_14 * pc_z[k] * hsl1_456[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pa_x, pc_x, pc_y, pc_z, hsl0_729, hsk_363, \
                         hsk_401, hsk_585, hsl1_729, isk_579, isk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_15 * hsk_363[k]
                   + f_3 * pc_z[k] * isk_579[k];

        t_728[k] = f_18 * hsk_401[k]
                   + f_3 * pc_y[k] * isk_581[k];

        t_729[k] = pa_x[k] * hsl0_729[k]
                   + f_19 * hsk_585[k]
                   - f_14 * pc_x[k] * hsl1_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pa_x, pa_z, pc_x, pc_z, hsl0_460, hsl0_732, \
                         hsk_366, hsk_588, hsl1_460, hsl1_732, \
                         isk_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = pa_z[k] * hsl0_460[k]
                   - f_14 * pc_z[k] * hsl1_460[k];

        t_731[k] = f_15 * hsk_366[k]
                   + f_3 * pc_z[k] * isk_582[k];

        t_732[k] = pa_x[k] * hsl0_732[k]
                   + f_18 * hsk_588[k]
                   - f_14 * pc_x[k] * hsl1_732[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pa_x, pa_z, pc_x, pc_y, pc_z, hsl0_465, \
                         hsl0_734, hsk_405, hsk_590, hsl1_465, hsl1_734, \
                         isk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_18 * hsk_405[k]
                   + f_3 * pc_y[k] * isk_585[k];

        t_734[k] = pa_x[k] * hsl0_734[k]
                   + f_18 * hsk_590[k]
                   - f_14 * pc_x[k] * hsl1_734[k];

        t_735[k] = pa_z[k] * hsl0_465[k]
                   - f_14 * pc_z[k] * hsl1_465[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pa_x, pc_x, pc_z, hsl0_737, hsl0_738, hsk_370, \
                         hsk_593, hsk_594, hsl1_737, hsl1_738, \
                         isk_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_15 * hsk_370[k]
                   + f_3 * pc_z[k] * isk_586[k];

        t_737[k] = pa_x[k] * hsl0_737[k]
                   + f_17 * hsk_593[k]
                   - f_14 * pc_x[k] * hsl1_737[k];

        t_738[k] = pa_x[k] * hsl0_738[k]
                   + f_17 * hsk_594[k]
                   - f_14 * pc_x[k] * hsl1_738[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pa_x, pa_z, pc_x, pc_y, pc_z, hsl0_471, \
                         hsl0_740, hsk_410, hsk_596, hsl1_471, hsl1_740, \
                         isk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_18 * hsk_410[k]
                   + f_3 * pc_y[k] * isk_590[k];

        t_740[k] = pa_x[k] * hsl0_740[k]
                   + f_17 * hsk_596[k]
                   - f_14 * pc_x[k] * hsl1_740[k];

        t_741[k] = pa_z[k] * hsl0_471[k]
                   - f_14 * pc_z[k] * hsl1_471[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, pa_x, pc_x, pc_z, hsl0_743, hsl0_744, hsk_375, \
                         hsk_599, hsk_600, hsl1_743, hsl1_744, \
                         isk_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_15 * hsk_375[k]
                   + f_3 * pc_z[k] * isk_591[k];

        t_743[k] = pa_x[k] * hsl0_743[k]
                   + f_16 * hsk_599[k]
                   - f_14 * pc_x[k] * hsl1_743[k];

        t_744[k] = pa_x[k] * hsl0_744[k]
                   + f_16 * hsk_600[k]
                   - f_14 * pc_x[k] * hsl1_744[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, pa_x, pc_x, pc_y, hsl0_745, hsl0_747, hsk_416, \
                         hsk_601, hsk_603, hsl1_745, hsl1_747, \
                         isk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = pa_x[k] * hsl0_745[k]
                   + f_16 * hsk_601[k]
                   - f_14 * pc_x[k] * hsl1_745[k];

        t_746[k] = f_18 * hsk_416[k]
                   + f_3 * pc_y[k] * isk_596[k];

        t_747[k] = pa_x[k] * hsl0_747[k]
                   + f_16 * hsk_603[k]
                   - f_14 * pc_x[k] * hsl1_747[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, t_752, pc_x, hsk_604, hsk_605, hsk_606, \
                         hsk_607, hsk_608, isk_604, isk_605, isk_606, isk_607, \
                         isk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_15 * hsk_604[k]
                   + f_3 * pc_x[k] * isk_604[k];

        t_749[k] = f_15 * hsk_605[k]
                   + f_3 * pc_x[k] * isk_605[k];

        t_750[k] = f_15 * hsk_606[k]
                   + f_3 * pc_x[k] * isk_606[k];

        t_751[k] = f_15 * hsk_607[k]
                   + f_3 * pc_x[k] * isk_607[k];

        t_752[k] = f_15 * hsk_608[k]
                   + f_3 * pc_x[k] * isk_608[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, pa_x, pc_x, hsl0_756, hsk_609, hsk_610, \
                         hsk_611, hsl1_756, isk_609, isk_610, isk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_15 * hsk_609[k]
                   + f_3 * pc_x[k] * isk_609[k];

        t_754[k] = f_15 * hsk_610[k]
                   + f_3 * pc_x[k] * isk_610[k];

        t_755[k] = f_15 * hsk_611[k]
                   + f_3 * pc_x[k] * isk_611[k];

        t_756[k] = pa_x[k] * hsl0_756[k]
                   - f_14 * pc_x[k] * hsl1_756[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pa_x, pc_x, pc_z, hsl0_758, hsl0_759, \
                         hsl0_760, hsk_388, hsl1_758, hsl1_759, hsl1_760, \
                         isk_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_15 * hsk_388[k]
                   + f_3 * pc_z[k] * isk_604[k];

        t_758[k] = pa_x[k] * hsl0_758[k]
                   - f_14 * pc_x[k] * hsl1_758[k];

        t_759[k] = pa_x[k] * hsl0_759[k]
                   - f_14 * pc_x[k] * hsl1_759[k];

        t_760[k] = pa_x[k] * hsl0_760[k]
                   - f_14 * pc_x[k] * hsl1_760[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, pa_x, pc_x, pc_y, hsl0_761, hsl0_762, \
                         hsl0_764, hsk_431, hsl1_761, hsl1_762, hsl1_764, \
                         isk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = pa_x[k] * hsl0_761[k]
                   - f_14 * pc_x[k] * hsl1_761[k];

        t_762[k] = pa_x[k] * hsl0_762[k]
                   - f_14 * pc_x[k] * hsl1_762[k];

        t_763[k] = f_18 * hsk_431[k]
                   + f_3 * pc_y[k] * isk_611[k];

        t_764[k] = pa_x[k] * hsl0_764[k]
                   - f_14 * pc_x[k] * hsl1_764[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pa_x, pc_x, pc_y, pc_z, hsl0_765, hsk_396, \
                         hsk_432, hsk_612, hsl1_765, isk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = pa_x[k] * hsl0_765[k]
                   + f_22 * hsk_612[k]
                   - f_14 * pc_x[k] * hsl1_765[k];

        t_766[k] = f_17 * hsk_432[k]
                   + f_3 * pc_y[k] * isk_612[k];

        t_767[k] = f_16 * hsk_396[k]
                   + f_3 * pc_z[k] * isk_612[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pa_x, pc_x, pc_y, hsl0_768, hsl0_770, hsk_434, \
                         hsk_615, hsk_617, hsl1_768, hsl1_770, \
                         isk_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = pa_x[k] * hsl0_768[k]
                   + f_0 * hsk_615[k]
                   - f_14 * pc_x[k] * hsl1_768[k];

        t_769[k] = f_17 * hsk_434[k]
                   + f_3 * pc_y[k] * isk_614[k];

        t_770[k] = pa_x[k] * hsl0_770[k]
                   + f_0 * hsk_617[k]
                   - f_14 * pc_x[k] * hsl1_770[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, pa_x, pc_x, pc_y, pc_z, hsl0_771, hsk_399, \
                         hsk_437, hsk_618, hsl1_771, isk_615, isk_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = pa_x[k] * hsl0_771[k]
                   + f_19 * hsk_618[k]
                   - f_14 * pc_x[k] * hsl1_771[k];

        t_772[k] = f_16 * hsk_399[k]
                   + f_3 * pc_z[k] * isk_615[k];

        t_773[k] = f_17 * hsk_437[k]
                   + f_3 * pc_y[k] * isk_617[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, pa_x, pc_x, pc_z, hsl0_774, hsl0_775, hsk_402, \
                         hsk_621, hsk_622, hsl1_774, hsl1_775, \
                         isk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = pa_x[k] * hsl0_774[k]
                   + f_19 * hsk_621[k]
                   - f_14 * pc_x[k] * hsl1_774[k];

        t_775[k] = pa_x[k] * hsl0_775[k]
                   + f_18 * hsk_622[k]
                   - f_14 * pc_x[k] * hsl1_775[k];

        t_776[k] = f_16 * hsk_402[k]
                   + f_3 * pc_z[k] * isk_618[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pa_x, pc_x, pc_y, hsl0_777, hsl0_779, hsk_441, \
                         hsk_624, hsk_626, hsl1_777, hsl1_779, \
                         isk_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = pa_x[k] * hsl0_777[k]
                   + f_18 * hsk_624[k]
                   - f_14 * pc_x[k] * hsl1_777[k];

        t_778[k] = f_17 * hsk_441[k]
                   + f_3 * pc_y[k] * isk_621[k];

        t_779[k] = pa_x[k] * hsl0_779[k]
                   + f_18 * hsk_626[k]
                   - f_14 * pc_x[k] * hsl1_779[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pa_x, pc_x, pc_z, hsl0_780, hsl0_782, hsk_406, \
                         hsk_627, hsk_629, hsl1_780, hsl1_782, \
                         isk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = pa_x[k] * hsl0_780[k]
                   + f_17 * hsk_627[k]
                   - f_14 * pc_x[k] * hsl1_780[k];

        t_781[k] = f_16 * hsk_406[k]
                   + f_3 * pc_z[k] * isk_622[k];

        t_782[k] = pa_x[k] * hsl0_782[k]
                   + f_17 * hsk_629[k]
                   - f_14 * pc_x[k] * hsl1_782[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pa_x, pc_x, pc_y, hsl0_783, hsl0_785, hsk_446, \
                         hsk_630, hsk_632, hsl1_783, hsl1_785, \
                         isk_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = pa_x[k] * hsl0_783[k]
                   + f_17 * hsk_630[k]
                   - f_14 * pc_x[k] * hsl1_783[k];

        t_784[k] = f_17 * hsk_446[k]
                   + f_3 * pc_y[k] * isk_626[k];

        t_785[k] = pa_x[k] * hsl0_785[k]
                   + f_17 * hsk_632[k]
                   - f_14 * pc_x[k] * hsl1_785[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, pa_x, pc_x, pc_z, hsl0_786, hsl0_788, hsk_411, \
                         hsk_633, hsk_635, hsl1_786, hsl1_788, \
                         isk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = pa_x[k] * hsl0_786[k]
                   + f_16 * hsk_633[k]
                   - f_14 * pc_x[k] * hsl1_786[k];

        t_787[k] = f_16 * hsk_411[k]
                   + f_3 * pc_z[k] * isk_627[k];

        t_788[k] = pa_x[k] * hsl0_788[k]
                   + f_16 * hsk_635[k]
                   - f_14 * pc_x[k] * hsl1_788[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, pa_x, pc_x, pc_y, hsl0_789, hsl0_790, hsk_452, \
                         hsk_636, hsk_637, hsl1_789, hsl1_790, \
                         isk_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = pa_x[k] * hsl0_789[k]
                   + f_16 * hsk_636[k]
                   - f_14 * pc_x[k] * hsl1_789[k];

        t_790[k] = pa_x[k] * hsl0_790[k]
                   + f_16 * hsk_637[k]
                   - f_14 * pc_x[k] * hsl1_790[k];

        t_791[k] = f_17 * hsk_452[k]
                   + f_3 * pc_y[k] * isk_632[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pa_x, pc_x, hsl0_792, hsk_639, hsk_640, \
                         hsk_641, hsk_642, hsl1_792, isk_640, isk_641, \
                         isk_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = pa_x[k] * hsl0_792[k]
                   + f_16 * hsk_639[k]
                   - f_14 * pc_x[k] * hsl1_792[k];

        t_793[k] = f_15 * hsk_640[k]
                   + f_3 * pc_x[k] * isk_640[k];

        t_794[k] = f_15 * hsk_641[k]
                   + f_3 * pc_x[k] * isk_641[k];

        t_795[k] = f_15 * hsk_642[k]
                   + f_3 * pc_x[k] * isk_642[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, t_800, pc_x, hsk_643, hsk_644, hsk_645, \
                         hsk_646, hsk_647, isk_643, isk_644, isk_645, isk_646, \
                         isk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_15 * hsk_643[k]
                   + f_3 * pc_x[k] * isk_643[k];

        t_797[k] = f_15 * hsk_644[k]
                   + f_3 * pc_x[k] * isk_644[k];

        t_798[k] = f_15 * hsk_645[k]
                   + f_3 * pc_x[k] * isk_645[k];

        t_799[k] = f_15 * hsk_646[k]
                   + f_3 * pc_x[k] * isk_646[k];

        t_800[k] = f_15 * hsk_647[k]
                   + f_3 * pc_x[k] * isk_647[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, t_804, pa_x, pc_x, pc_z, hsl0_801, hsl0_803, \
                         hsl0_804, hsk_424, hsl1_801, hsl1_803, hsl1_804, \
                         isk_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = pa_x[k] * hsl0_801[k]
                   - f_14 * pc_x[k] * hsl1_801[k];

        t_802[k] = f_16 * hsk_424[k]
                   + f_3 * pc_z[k] * isk_640[k];

        t_803[k] = pa_x[k] * hsl0_803[k]
                   - f_14 * pc_x[k] * hsl1_803[k];

        t_804[k] = pa_x[k] * hsl0_804[k]
                   - f_14 * pc_x[k] * hsl1_804[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, pa_x, pc_x, pc_y, hsl0_805, hsl0_806, \
                         hsl0_807, hsk_467, hsl1_805, hsl1_806, hsl1_807, \
                         isk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = pa_x[k] * hsl0_805[k]
                   - f_14 * pc_x[k] * hsl1_805[k];

        t_806[k] = pa_x[k] * hsl0_806[k]
                   - f_14 * pc_x[k] * hsl1_806[k];

        t_807[k] = pa_x[k] * hsl0_807[k]
                   - f_14 * pc_x[k] * hsl1_807[k];

        t_808[k] = f_17 * hsk_467[k]
                   + f_3 * pc_y[k] * isk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, pa_x, pc_x, pc_y, pc_z, hsl0_809, \
                         hsl0_810, hsk_432, hsk_468, hsk_648, hsl1_809, hsl1_810, \
                         isk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pa_x[k] * hsl0_809[k]
                   - f_14 * pc_x[k] * hsl1_809[k];

        t_810[k] = pa_x[k] * hsl0_810[k]
                   + f_22 * hsk_648[k]
                   - f_14 * pc_x[k] * hsl1_810[k];

        t_811[k] = f_16 * hsk_468[k]
                   + f_3 * pc_y[k] * isk_648[k];

        t_812[k] = f_17 * hsk_432[k]
                   + f_3 * pc_z[k] * isk_648[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, pa_x, pc_x, pc_y, hsl0_813, hsl0_815, hsk_470, \
                         hsk_651, hsk_653, hsl1_813, hsl1_815, \
                         isk_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = pa_x[k] * hsl0_813[k]
                   + f_0 * hsk_651[k]
                   - f_14 * pc_x[k] * hsl1_813[k];

        t_814[k] = f_16 * hsk_470[k]
                   + f_3 * pc_y[k] * isk_650[k];

        t_815[k] = pa_x[k] * hsl0_815[k]
                   + f_0 * hsk_653[k]
                   - f_14 * pc_x[k] * hsl1_815[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pa_x, pc_x, pc_y, pc_z, hsl0_816, hsk_435, \
                         hsk_473, hsk_654, hsl1_816, isk_651, isk_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = pa_x[k] * hsl0_816[k]
                   + f_19 * hsk_654[k]
                   - f_14 * pc_x[k] * hsl1_816[k];

        t_817[k] = f_17 * hsk_435[k]
                   + f_3 * pc_z[k] * isk_651[k];

        t_818[k] = f_16 * hsk_473[k]
                   + f_3 * pc_y[k] * isk_653[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pa_x, pc_x, pc_z, hsl0_819, hsl0_820, hsk_438, \
                         hsk_657, hsk_658, hsl1_819, hsl1_820, \
                         isk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = pa_x[k] * hsl0_819[k]
                   + f_19 * hsk_657[k]
                   - f_14 * pc_x[k] * hsl1_819[k];

        t_820[k] = pa_x[k] * hsl0_820[k]
                   + f_18 * hsk_658[k]
                   - f_14 * pc_x[k] * hsl1_820[k];

        t_821[k] = f_17 * hsk_438[k]
                   + f_3 * pc_z[k] * isk_654[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_22 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsl0_630 = buffer.data(hsl0 + 630);
    const auto *hsl0_635 = buffer.data(hsl0 + 635);
    const auto *hsl0_639 = buffer.data(hsl0 + 639);
    const auto *hsl0_644 = buffer.data(hsl0 + 644);
    const auto *hsl0_650 = buffer.data(hsl0 + 650);
    const auto *hsl0_657 = buffer.data(hsl0 + 657);
    const auto *hsl0_822 = buffer.data(hsl0 + 822);
    const auto *hsl0_824 = buffer.data(hsl0 + 824);
    const auto *hsl0_825 = buffer.data(hsl0 + 825);
    const auto *hsl0_827 = buffer.data(hsl0 + 827);
    const auto *hsl0_828 = buffer.data(hsl0 + 828);
    const auto *hsl0_830 = buffer.data(hsl0 + 830);
    const auto *hsl0_831 = buffer.data(hsl0 + 831);
    const auto *hsl0_833 = buffer.data(hsl0 + 833);
    const auto *hsl0_834 = buffer.data(hsl0 + 834);
    const auto *hsl0_835 = buffer.data(hsl0 + 835);
    const auto *hsl0_837 = buffer.data(hsl0 + 837);
    const auto *hsl0_846 = buffer.data(hsl0 + 846);
    const auto *hsl0_848 = buffer.data(hsl0 + 848);
    const auto *hsl0_849 = buffer.data(hsl0 + 849);
    const auto *hsl0_850 = buffer.data(hsl0 + 850);
    const auto *hsl0_851 = buffer.data(hsl0 + 851);
    const auto *hsl0_852 = buffer.data(hsl0 + 852);
    const auto *hsl0_854 = buffer.data(hsl0 + 854);
    const auto *hsl0_858 = buffer.data(hsl0 + 858);
    const auto *hsl0_861 = buffer.data(hsl0 + 861);
    const auto *hsl0_865 = buffer.data(hsl0 + 865);
    const auto *hsl0_867 = buffer.data(hsl0 + 867);
    const auto *hsl0_870 = buffer.data(hsl0 + 870);
    const auto *hsl0_872 = buffer.data(hsl0 + 872);
    const auto *hsl0_873 = buffer.data(hsl0 + 873);
    const auto *hsl0_876 = buffer.data(hsl0 + 876);
    const auto *hsl0_878 = buffer.data(hsl0 + 878);
    const auto *hsl0_879 = buffer.data(hsl0 + 879);
    const auto *hsl0_880 = buffer.data(hsl0 + 880);
    const auto *hsl0_891 = buffer.data(hsl0 + 891);
    const auto *hsl0_893 = buffer.data(hsl0 + 893);
    const auto *hsl0_894 = buffer.data(hsl0 + 894);
    const auto *hsl0_895 = buffer.data(hsl0 + 895);
    const auto *hsl0_896 = buffer.data(hsl0 + 896);
    const auto *hsl0_897 = buffer.data(hsl0 + 897);
    const auto *hsl0_899 = buffer.data(hsl0 + 899);
    const auto *hsl0_900 = buffer.data(hsl0 + 900);
    const auto *hsl0_905 = buffer.data(hsl0 + 905);
    const auto *hsl0_909 = buffer.data(hsl0 + 909);
    const auto *hsl0_914 = buffer.data(hsl0 + 914);
    const auto *hsl0_920 = buffer.data(hsl0 + 920);
    const auto *hsl0_927 = buffer.data(hsl0 + 927);
    const auto *hsl0_936 = buffer.data(hsl0 + 936);
    const auto *hsl0_937 = buffer.data(hsl0 + 937);
    const auto *hsl0_938 = buffer.data(hsl0 + 938);
    const auto *hsl0_939 = buffer.data(hsl0 + 939);

    const auto *hsk_442 = buffer.data(hsk + 442);
    const auto *hsk_447 = buffer.data(hsk + 447);
    const auto *hsk_460 = buffer.data(hsk + 460);
    const auto *hsk_468 = buffer.data(hsk + 468);
    const auto *hsk_471 = buffer.data(hsk + 471);
    const auto *hsk_474 = buffer.data(hsk + 474);
    const auto *hsk_477 = buffer.data(hsk + 477);
    const auto *hsk_478 = buffer.data(hsk + 478);
    const auto *hsk_482 = buffer.data(hsk + 482);
    const auto *hsk_483 = buffer.data(hsk + 483);
    const auto *hsk_488 = buffer.data(hsk + 488);
    const auto *hsk_496 = buffer.data(hsk + 496);
    const auto *hsk_503 = buffer.data(hsk + 503);
    const auto *hsk_504 = buffer.data(hsk + 504);
    const auto *hsk_506 = buffer.data(hsk + 506);
    const auto *hsk_509 = buffer.data(hsk + 509);
    const auto *hsk_513 = buffer.data(hsk + 513);
    const auto *hsk_518 = buffer.data(hsk + 518);
    const auto *hsk_524 = buffer.data(hsk + 524);
    const auto *hsk_539 = buffer.data(hsk + 539);
    const auto *hsk_660 = buffer.data(hsk + 660);
    const auto *hsk_662 = buffer.data(hsk + 662);
    const auto *hsk_663 = buffer.data(hsk + 663);
    const auto *hsk_665 = buffer.data(hsk + 665);
    const auto *hsk_666 = buffer.data(hsk + 666);
    const auto *hsk_668 = buffer.data(hsk + 668);
    const auto *hsk_669 = buffer.data(hsk + 669);
    const auto *hsk_671 = buffer.data(hsk + 671);
    const auto *hsk_672 = buffer.data(hsk + 672);
    const auto *hsk_673 = buffer.data(hsk + 673);
    const auto *hsk_675 = buffer.data(hsk + 675);
    const auto *hsk_676 = buffer.data(hsk + 676);
    const auto *hsk_677 = buffer.data(hsk + 677);
    const auto *hsk_678 = buffer.data(hsk + 678);
    const auto *hsk_679 = buffer.data(hsk + 679);
    const auto *hsk_680 = buffer.data(hsk + 680);
    const auto *hsk_681 = buffer.data(hsk + 681);
    const auto *hsk_682 = buffer.data(hsk + 682);
    const auto *hsk_683 = buffer.data(hsk + 683);
    const auto *hsk_687 = buffer.data(hsk + 687);
    const auto *hsk_690 = buffer.data(hsk + 690);
    const auto *hsk_694 = buffer.data(hsk + 694);
    const auto *hsk_696 = buffer.data(hsk + 696);
    const auto *hsk_699 = buffer.data(hsk + 699);
    const auto *hsk_701 = buffer.data(hsk + 701);
    const auto *hsk_702 = buffer.data(hsk + 702);
    const auto *hsk_705 = buffer.data(hsk + 705);
    const auto *hsk_707 = buffer.data(hsk + 707);
    const auto *hsk_708 = buffer.data(hsk + 708);
    const auto *hsk_709 = buffer.data(hsk + 709);
    const auto *hsk_712 = buffer.data(hsk + 712);
    const auto *hsk_713 = buffer.data(hsk + 713);
    const auto *hsk_714 = buffer.data(hsk + 714);
    const auto *hsk_715 = buffer.data(hsk + 715);
    const auto *hsk_716 = buffer.data(hsk + 716);
    const auto *hsk_717 = buffer.data(hsk + 717);
    const auto *hsk_718 = buffer.data(hsk + 718);
    const auto *hsk_719 = buffer.data(hsk + 719);
    const auto *hsk_720 = buffer.data(hsk + 720);
    const auto *hsk_725 = buffer.data(hsk + 725);
    const auto *hsk_729 = buffer.data(hsk + 729);
    const auto *hsk_734 = buffer.data(hsk + 734);
    const auto *hsk_740 = buffer.data(hsk + 740);
    const auto *hsk_747 = buffer.data(hsk + 747);
    const auto *hsk_748 = buffer.data(hsk + 748);
    const auto *hsk_749 = buffer.data(hsk + 749);
    const auto *hsk_750 = buffer.data(hsk + 750);
    const auto *hsk_751 = buffer.data(hsk + 751);
    const auto *hsk_752 = buffer.data(hsk + 752);
    const auto *hsk_753 = buffer.data(hsk + 753);
    const auto *hsk_755 = buffer.data(hsk + 755);

    const auto *hsl1_630 = buffer.data(hsl1 + 630);
    const auto *hsl1_635 = buffer.data(hsl1 + 635);
    const auto *hsl1_639 = buffer.data(hsl1 + 639);
    const auto *hsl1_644 = buffer.data(hsl1 + 644);
    const auto *hsl1_650 = buffer.data(hsl1 + 650);
    const auto *hsl1_657 = buffer.data(hsl1 + 657);
    const auto *hsl1_822 = buffer.data(hsl1 + 822);
    const auto *hsl1_824 = buffer.data(hsl1 + 824);
    const auto *hsl1_825 = buffer.data(hsl1 + 825);
    const auto *hsl1_827 = buffer.data(hsl1 + 827);
    const auto *hsl1_828 = buffer.data(hsl1 + 828);
    const auto *hsl1_830 = buffer.data(hsl1 + 830);
    const auto *hsl1_831 = buffer.data(hsl1 + 831);
    const auto *hsl1_833 = buffer.data(hsl1 + 833);
    const auto *hsl1_834 = buffer.data(hsl1 + 834);
    const auto *hsl1_835 = buffer.data(hsl1 + 835);
    const auto *hsl1_837 = buffer.data(hsl1 + 837);
    const auto *hsl1_846 = buffer.data(hsl1 + 846);
    const auto *hsl1_848 = buffer.data(hsl1 + 848);
    const auto *hsl1_849 = buffer.data(hsl1 + 849);
    const auto *hsl1_850 = buffer.data(hsl1 + 850);
    const auto *hsl1_851 = buffer.data(hsl1 + 851);
    const auto *hsl1_852 = buffer.data(hsl1 + 852);
    const auto *hsl1_854 = buffer.data(hsl1 + 854);
    const auto *hsl1_858 = buffer.data(hsl1 + 858);
    const auto *hsl1_861 = buffer.data(hsl1 + 861);
    const auto *hsl1_865 = buffer.data(hsl1 + 865);
    const auto *hsl1_867 = buffer.data(hsl1 + 867);
    const auto *hsl1_870 = buffer.data(hsl1 + 870);
    const auto *hsl1_872 = buffer.data(hsl1 + 872);
    const auto *hsl1_873 = buffer.data(hsl1 + 873);
    const auto *hsl1_876 = buffer.data(hsl1 + 876);
    const auto *hsl1_878 = buffer.data(hsl1 + 878);
    const auto *hsl1_879 = buffer.data(hsl1 + 879);
    const auto *hsl1_880 = buffer.data(hsl1 + 880);
    const auto *hsl1_891 = buffer.data(hsl1 + 891);
    const auto *hsl1_893 = buffer.data(hsl1 + 893);
    const auto *hsl1_894 = buffer.data(hsl1 + 894);
    const auto *hsl1_895 = buffer.data(hsl1 + 895);
    const auto *hsl1_896 = buffer.data(hsl1 + 896);
    const auto *hsl1_897 = buffer.data(hsl1 + 897);
    const auto *hsl1_899 = buffer.data(hsl1 + 899);
    const auto *hsl1_900 = buffer.data(hsl1 + 900);
    const auto *hsl1_905 = buffer.data(hsl1 + 905);
    const auto *hsl1_909 = buffer.data(hsl1 + 909);
    const auto *hsl1_914 = buffer.data(hsl1 + 914);
    const auto *hsl1_920 = buffer.data(hsl1 + 920);
    const auto *hsl1_927 = buffer.data(hsl1 + 927);
    const auto *hsl1_936 = buffer.data(hsl1 + 936);
    const auto *hsl1_937 = buffer.data(hsl1 + 937);
    const auto *hsl1_938 = buffer.data(hsl1 + 938);
    const auto *hsl1_939 = buffer.data(hsl1 + 939);

    const auto *isi0_560 = buffer.data(isi0 + 560);
    const auto *isi0_561 = buffer.data(isi0 + 561);
    const auto *isi0_562 = buffer.data(isi0 + 562);
    const auto *isi0_563 = buffer.data(isi0 + 563);
    const auto *isi0_564 = buffer.data(isi0 + 564);
    const auto *isi0_565 = buffer.data(isi0 + 565);
    const auto *isi0_566 = buffer.data(isi0 + 566);
    const auto *isi0_567 = buffer.data(isi0 + 567);
    const auto *isi0_568 = buffer.data(isi0 + 568);
    const auto *isi0_569 = buffer.data(isi0 + 569);
    const auto *isi0_570 = buffer.data(isi0 + 570);
    const auto *isi0_571 = buffer.data(isi0 + 571);
    const auto *isi0_572 = buffer.data(isi0 + 572);
    const auto *isi0_573 = buffer.data(isi0 + 573);
    const auto *isi0_574 = buffer.data(isi0 + 574);

    const auto *isi1_560 = buffer.data(isi1 + 560);
    const auto *isi1_561 = buffer.data(isi1 + 561);
    const auto *isi1_562 = buffer.data(isi1 + 562);
    const auto *isi1_563 = buffer.data(isi1 + 563);
    const auto *isi1_564 = buffer.data(isi1 + 564);
    const auto *isi1_565 = buffer.data(isi1 + 565);
    const auto *isi1_566 = buffer.data(isi1 + 566);
    const auto *isi1_567 = buffer.data(isi1 + 567);
    const auto *isi1_568 = buffer.data(isi1 + 568);
    const auto *isi1_569 = buffer.data(isi1 + 569);
    const auto *isi1_570 = buffer.data(isi1 + 570);
    const auto *isi1_571 = buffer.data(isi1 + 571);
    const auto *isi1_572 = buffer.data(isi1 + 572);
    const auto *isi1_573 = buffer.data(isi1 + 573);
    const auto *isi1_574 = buffer.data(isi1 + 574);

    const auto *isk_657 = buffer.data(isk + 657);
    const auto *isk_658 = buffer.data(isk + 658);
    const auto *isk_662 = buffer.data(isk + 662);
    const auto *isk_663 = buffer.data(isk + 663);
    const auto *isk_668 = buffer.data(isk + 668);
    const auto *isk_676 = buffer.data(isk + 676);
    const auto *isk_677 = buffer.data(isk + 677);
    const auto *isk_678 = buffer.data(isk + 678);
    const auto *isk_679 = buffer.data(isk + 679);
    const auto *isk_680 = buffer.data(isk + 680);
    const auto *isk_681 = buffer.data(isk + 681);
    const auto *isk_682 = buffer.data(isk + 682);
    const auto *isk_683 = buffer.data(isk + 683);
    const auto *isk_684 = buffer.data(isk + 684);
    const auto *isk_686 = buffer.data(isk + 686);
    const auto *isk_687 = buffer.data(isk + 687);
    const auto *isk_689 = buffer.data(isk + 689);
    const auto *isk_690 = buffer.data(isk + 690);
    const auto *isk_693 = buffer.data(isk + 693);
    const auto *isk_694 = buffer.data(isk + 694);
    const auto *isk_698 = buffer.data(isk + 698);
    const auto *isk_699 = buffer.data(isk + 699);
    const auto *isk_704 = buffer.data(isk + 704);
    const auto *isk_712 = buffer.data(isk + 712);
    const auto *isk_713 = buffer.data(isk + 713);
    const auto *isk_714 = buffer.data(isk + 714);
    const auto *isk_715 = buffer.data(isk + 715);
    const auto *isk_716 = buffer.data(isk + 716);
    const auto *isk_717 = buffer.data(isk + 717);
    const auto *isk_718 = buffer.data(isk + 718);
    const auto *isk_719 = buffer.data(isk + 719);
    const auto *isk_720 = buffer.data(isk + 720);
    const auto *isk_721 = buffer.data(isk + 721);
    const auto *isk_722 = buffer.data(isk + 722);
    const auto *isk_723 = buffer.data(isk + 723);
    const auto *isk_724 = buffer.data(isk + 724);
    const auto *isk_725 = buffer.data(isk + 725);
    const auto *isk_726 = buffer.data(isk + 726);
    const auto *isk_727 = buffer.data(isk + 727);
    const auto *isk_728 = buffer.data(isk + 728);
    const auto *isk_729 = buffer.data(isk + 729);
    const auto *isk_730 = buffer.data(isk + 730);
    const auto *isk_731 = buffer.data(isk + 731);
    const auto *isk_732 = buffer.data(isk + 732);
    const auto *isk_733 = buffer.data(isk + 733);
    const auto *isk_734 = buffer.data(isk + 734);
    const auto *isk_735 = buffer.data(isk + 735);
    const auto *isk_736 = buffer.data(isk + 736);
    const auto *isk_737 = buffer.data(isk + 737);
    const auto *isk_738 = buffer.data(isk + 738);
    const auto *isk_739 = buffer.data(isk + 739);
    const auto *isk_740 = buffer.data(isk + 740);
    const auto *isk_747 = buffer.data(isk + 747);
    const auto *isk_748 = buffer.data(isk + 748);
    const auto *isk_749 = buffer.data(isk + 749);
    const auto *isk_750 = buffer.data(isk + 750);
    const auto *isk_751 = buffer.data(isk + 751);
    const auto *isk_752 = buffer.data(isk + 752);
    const auto *isk_753 = buffer.data(isk + 753);
    const auto *isk_755 = buffer.data(isk + 755);

#pragma omp simd aligned(t_822, t_823, t_824, pa_x, pc_x, pc_y, hsl0_822, hsl0_824, hsk_477, \
                         hsk_660, hsk_662, hsl1_822, hsl1_824, \
                         isk_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pa_x[k] * hsl0_822[k]
                   + f_18 * hsk_660[k]
                   - f_14 * pc_x[k] * hsl1_822[k];

        t_823[k] = f_16 * hsk_477[k]
                   + f_3 * pc_y[k] * isk_657[k];

        t_824[k] = pa_x[k] * hsl0_824[k]
                   + f_18 * hsk_662[k]
                   - f_14 * pc_x[k] * hsl1_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, pa_x, pc_x, pc_z, hsl0_825, hsl0_827, hsk_442, \
                         hsk_663, hsk_665, hsl1_825, hsl1_827, \
                         isk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = pa_x[k] * hsl0_825[k]
                   + f_17 * hsk_663[k]
                   - f_14 * pc_x[k] * hsl1_825[k];

        t_826[k] = f_17 * hsk_442[k]
                   + f_3 * pc_z[k] * isk_658[k];

        t_827[k] = pa_x[k] * hsl0_827[k]
                   + f_17 * hsk_665[k]
                   - f_14 * pc_x[k] * hsl1_827[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, pa_x, pc_x, pc_y, hsl0_828, hsl0_830, hsk_482, \
                         hsk_666, hsk_668, hsl1_828, hsl1_830, \
                         isk_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = pa_x[k] * hsl0_828[k]
                   + f_17 * hsk_666[k]
                   - f_14 * pc_x[k] * hsl1_828[k];

        t_829[k] = f_16 * hsk_482[k]
                   + f_3 * pc_y[k] * isk_662[k];

        t_830[k] = pa_x[k] * hsl0_830[k]
                   + f_17 * hsk_668[k]
                   - f_14 * pc_x[k] * hsl1_830[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, pa_x, pc_x, pc_z, hsl0_831, hsl0_833, hsk_447, \
                         hsk_669, hsk_671, hsl1_831, hsl1_833, \
                         isk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = pa_x[k] * hsl0_831[k]
                   + f_16 * hsk_669[k]
                   - f_14 * pc_x[k] * hsl1_831[k];

        t_832[k] = f_17 * hsk_447[k]
                   + f_3 * pc_z[k] * isk_663[k];

        t_833[k] = pa_x[k] * hsl0_833[k]
                   + f_16 * hsk_671[k]
                   - f_14 * pc_x[k] * hsl1_833[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pa_x, pc_x, pc_y, hsl0_834, hsl0_835, hsk_488, \
                         hsk_672, hsk_673, hsl1_834, hsl1_835, \
                         isk_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = pa_x[k] * hsl0_834[k]
                   + f_16 * hsk_672[k]
                   - f_14 * pc_x[k] * hsl1_834[k];

        t_835[k] = pa_x[k] * hsl0_835[k]
                   + f_16 * hsk_673[k]
                   - f_14 * pc_x[k] * hsl1_835[k];

        t_836[k] = f_16 * hsk_488[k]
                   + f_3 * pc_y[k] * isk_668[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, pa_x, pc_x, hsl0_837, hsk_675, hsk_676, \
                         hsk_677, hsk_678, hsl1_837, isk_676, isk_677, \
                         isk_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = pa_x[k] * hsl0_837[k]
                   + f_16 * hsk_675[k]
                   - f_14 * pc_x[k] * hsl1_837[k];

        t_838[k] = f_15 * hsk_676[k]
                   + f_3 * pc_x[k] * isk_676[k];

        t_839[k] = f_15 * hsk_677[k]
                   + f_3 * pc_x[k] * isk_677[k];

        t_840[k] = f_15 * hsk_678[k]
                   + f_3 * pc_x[k] * isk_678[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, t_845, pc_x, hsk_679, hsk_680, hsk_681, \
                         hsk_682, hsk_683, isk_679, isk_680, isk_681, isk_682, \
                         isk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_15 * hsk_679[k]
                   + f_3 * pc_x[k] * isk_679[k];

        t_842[k] = f_15 * hsk_680[k]
                   + f_3 * pc_x[k] * isk_680[k];

        t_843[k] = f_15 * hsk_681[k]
                   + f_3 * pc_x[k] * isk_681[k];

        t_844[k] = f_15 * hsk_682[k]
                   + f_3 * pc_x[k] * isk_682[k];

        t_845[k] = f_15 * hsk_683[k]
                   + f_3 * pc_x[k] * isk_683[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, t_849, pa_x, pc_x, pc_z, hsl0_846, hsl0_848, \
                         hsl0_849, hsk_460, hsl1_846, hsl1_848, hsl1_849, \
                         isk_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = pa_x[k] * hsl0_846[k]
                   - f_14 * pc_x[k] * hsl1_846[k];

        t_847[k] = f_17 * hsk_460[k]
                   + f_3 * pc_z[k] * isk_676[k];

        t_848[k] = pa_x[k] * hsl0_848[k]
                   - f_14 * pc_x[k] * hsl1_848[k];

        t_849[k] = pa_x[k] * hsl0_849[k]
                   - f_14 * pc_x[k] * hsl1_849[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, pa_x, pc_x, pc_y, hsl0_850, hsl0_851, \
                         hsl0_852, hsk_503, hsl1_850, hsl1_851, hsl1_852, \
                         isk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = pa_x[k] * hsl0_850[k]
                   - f_14 * pc_x[k] * hsl1_850[k];

        t_851[k] = pa_x[k] * hsl0_851[k]
                   - f_14 * pc_x[k] * hsl1_851[k];

        t_852[k] = pa_x[k] * hsl0_852[k]
                   - f_14 * pc_x[k] * hsl1_852[k];

        t_853[k] = f_16 * hsk_503[k]
                   + f_3 * pc_y[k] * isk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pa_x, pa_y, pc_x, pc_y, pc_z, hsl0_630, \
                         hsl0_854, hsk_468, hsk_504, hsl1_630, hsl1_854, \
                         isk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = pa_x[k] * hsl0_854[k]
                   - f_14 * pc_x[k] * hsl1_854[k];

        t_855[k] = pa_y[k] * hsl0_630[k]
                   - f_14 * pc_y[k] * hsl1_630[k];

        t_856[k] = f_15 * hsk_504[k]
                   + f_3 * pc_y[k] * isk_684[k];

        t_857[k] = f_18 * hsk_468[k]
                   + f_3 * pc_z[k] * isk_684[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, pa_x, pa_y, pc_x, pc_y, hsl0_635, hsl0_858, \
                         hsk_506, hsk_687, hsl1_635, hsl1_858, \
                         isk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pa_x[k] * hsl0_858[k]
                   + f_0 * hsk_687[k]
                   - f_14 * pc_x[k] * hsl1_858[k];

        t_859[k] = f_15 * hsk_506[k]
                   + f_3 * pc_y[k] * isk_686[k];

        t_860[k] = pa_y[k] * hsl0_635[k]
                   - f_14 * pc_y[k] * hsl1_635[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, pa_x, pc_x, pc_y, pc_z, hsl0_861, hsk_471, \
                         hsk_509, hsk_690, hsl1_861, isk_687, isk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = pa_x[k] * hsl0_861[k]
                   + f_19 * hsk_690[k]
                   - f_14 * pc_x[k] * hsl1_861[k];

        t_862[k] = f_18 * hsk_471[k]
                   + f_3 * pc_z[k] * isk_687[k];

        t_863[k] = f_15 * hsk_509[k]
                   + f_3 * pc_y[k] * isk_689[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, pa_x, pa_y, pc_x, pc_y, pc_z, hsl0_639, \
                         hsl0_865, hsk_474, hsk_694, hsl1_639, hsl1_865, \
                         isk_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = pa_y[k] * hsl0_639[k]
                   - f_14 * pc_y[k] * hsl1_639[k];

        t_865[k] = pa_x[k] * hsl0_865[k]
                   + f_18 * hsk_694[k]
                   - f_14 * pc_x[k] * hsl1_865[k];

        t_866[k] = f_18 * hsk_474[k]
                   + f_3 * pc_z[k] * isk_690[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, pa_x, pa_y, pc_x, pc_y, hsl0_644, hsl0_867, \
                         hsk_513, hsk_696, hsl1_644, hsl1_867, \
                         isk_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = pa_x[k] * hsl0_867[k]
                   + f_18 * hsk_696[k]
                   - f_14 * pc_x[k] * hsl1_867[k];

        t_868[k] = f_15 * hsk_513[k]
                   + f_3 * pc_y[k] * isk_693[k];

        t_869[k] = pa_y[k] * hsl0_644[k]
                   - f_14 * pc_y[k] * hsl1_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pa_x, pc_x, pc_z, hsl0_870, hsl0_872, hsk_478, \
                         hsk_699, hsk_701, hsl1_870, hsl1_872, \
                         isk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pa_x[k] * hsl0_870[k]
                   + f_17 * hsk_699[k]
                   - f_14 * pc_x[k] * hsl1_870[k];

        t_871[k] = f_18 * hsk_478[k]
                   + f_3 * pc_z[k] * isk_694[k];

        t_872[k] = pa_x[k] * hsl0_872[k]
                   + f_17 * hsk_701[k]
                   - f_14 * pc_x[k] * hsl1_872[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pa_x, pa_y, pc_x, pc_y, hsl0_650, hsl0_873, \
                         hsk_518, hsk_702, hsl1_650, hsl1_873, \
                         isk_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pa_x[k] * hsl0_873[k]
                   + f_17 * hsk_702[k]
                   - f_14 * pc_x[k] * hsl1_873[k];

        t_874[k] = f_15 * hsk_518[k]
                   + f_3 * pc_y[k] * isk_698[k];

        t_875[k] = pa_y[k] * hsl0_650[k]
                   - f_14 * pc_y[k] * hsl1_650[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pa_x, pc_x, pc_z, hsl0_876, hsl0_878, hsk_483, \
                         hsk_705, hsk_707, hsl1_876, hsl1_878, \
                         isk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = pa_x[k] * hsl0_876[k]
                   + f_16 * hsk_705[k]
                   - f_14 * pc_x[k] * hsl1_876[k];

        t_877[k] = f_18 * hsk_483[k]
                   + f_3 * pc_z[k] * isk_699[k];

        t_878[k] = pa_x[k] * hsl0_878[k]
                   + f_16 * hsk_707[k]
                   - f_14 * pc_x[k] * hsl1_878[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pa_x, pc_x, pc_y, hsl0_879, hsl0_880, hsk_524, \
                         hsk_708, hsk_709, hsl1_879, hsl1_880, \
                         isk_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = pa_x[k] * hsl0_879[k]
                   + f_16 * hsk_708[k]
                   - f_14 * pc_x[k] * hsl1_879[k];

        t_880[k] = pa_x[k] * hsl0_880[k]
                   + f_16 * hsk_709[k]
                   - f_14 * pc_x[k] * hsl1_880[k];

        t_881[k] = f_15 * hsk_524[k]
                   + f_3 * pc_y[k] * isk_704[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pa_y, pc_x, pc_y, hsl0_657, hsk_712, \
                         hsk_713, hsk_714, hsl1_657, isk_712, isk_713, \
                         isk_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pa_y[k] * hsl0_657[k]
                   - f_14 * pc_y[k] * hsl1_657[k];

        t_883[k] = f_15 * hsk_712[k]
                   + f_3 * pc_x[k] * isk_712[k];

        t_884[k] = f_15 * hsk_713[k]
                   + f_3 * pc_x[k] * isk_713[k];

        t_885[k] = f_15 * hsk_714[k]
                   + f_3 * pc_x[k] * isk_714[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, t_890, pc_x, hsk_715, hsk_716, hsk_717, \
                         hsk_718, hsk_719, isk_715, isk_716, isk_717, isk_718, \
                         isk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_15 * hsk_715[k]
                   + f_3 * pc_x[k] * isk_715[k];

        t_887[k] = f_15 * hsk_716[k]
                   + f_3 * pc_x[k] * isk_716[k];

        t_888[k] = f_15 * hsk_717[k]
                   + f_3 * pc_x[k] * isk_717[k];

        t_889[k] = f_15 * hsk_718[k]
                   + f_3 * pc_x[k] * isk_718[k];

        t_890[k] = f_15 * hsk_719[k]
                   + f_3 * pc_x[k] * isk_719[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, pa_x, pc_x, pc_z, hsl0_891, hsl0_893, \
                         hsl0_894, hsk_496, hsl1_891, hsl1_893, hsl1_894, \
                         isk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = pa_x[k] * hsl0_891[k]
                   - f_14 * pc_x[k] * hsl1_891[k];

        t_892[k] = f_18 * hsk_496[k]
                   + f_3 * pc_z[k] * isk_712[k];

        t_893[k] = pa_x[k] * hsl0_893[k]
                   - f_14 * pc_x[k] * hsl1_893[k];

        t_894[k] = pa_x[k] * hsl0_894[k]
                   - f_14 * pc_x[k] * hsl1_894[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, pa_x, pc_x, pc_y, hsl0_895, hsl0_896, \
                         hsl0_897, hsk_539, hsl1_895, hsl1_896, hsl1_897, \
                         isk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = pa_x[k] * hsl0_895[k]
                   - f_14 * pc_x[k] * hsl1_895[k];

        t_896[k] = pa_x[k] * hsl0_896[k]
                   - f_14 * pc_x[k] * hsl1_896[k];

        t_897[k] = pa_x[k] * hsl0_897[k]
                   - f_14 * pc_x[k] * hsl1_897[k];

        t_898[k] = f_15 * hsk_539[k]
                   + f_3 * pc_y[k] * isk_719[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pa_x, pc_x, pc_y, pc_z, hsl0_899, \
                         hsl0_900, hsk_504, hsk_720, hsl1_899, hsl1_900, \
                         isk_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pa_x[k] * hsl0_899[k]
                   - f_14 * pc_x[k] * hsl1_899[k];

        t_900[k] = pa_x[k] * hsl0_900[k]
                   + f_22 * hsk_720[k]
                   - f_14 * pc_x[k] * hsl1_900[k];

        t_901[k] = f_3 * pc_y[k] * isk_720[k];

        t_902[k] = f_19 * hsk_504[k]
                   + f_3 * pc_z[k] * isk_720[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pa_x, pc_x, pc_y, hsl0_905, hsk_725, hsl1_905, \
                         isi0_560, isi1_560, isk_721, isk_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_4 * isi0_560[k]
                   - f_5 * isi1_560[k]
                   + f_3 * pc_y[k] * isk_721[k];

        t_904[k] = f_3 * pc_y[k] * isk_722[k];

        t_905[k] = pa_x[k] * hsl0_905[k]
                   + f_0 * hsk_725[k]
                   - f_14 * pc_x[k] * hsl1_905[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pc_y, isi0_561, isi0_562, isi1_561, isi1_562, \
                         isk_723, isk_724, isk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_6 * isi0_561[k]
                   - f_7 * isi1_561[k]
                   + f_3 * pc_y[k] * isk_723[k];

        t_907[k] = f_4 * isi0_562[k]
                   - f_5 * isi1_562[k]
                   + f_3 * pc_y[k] * isk_724[k];

        t_908[k] = f_3 * pc_y[k] * isk_725[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pa_x, pc_x, pc_y, hsl0_909, hsk_729, hsl1_909, \
                         isi0_563, isi0_564, isi1_563, isi1_564, isk_726, \
                         isk_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pa_x[k] * hsl0_909[k]
                   + f_19 * hsk_729[k]
                   - f_14 * pc_x[k] * hsl1_909[k];

        t_910[k] = f_8 * isi0_563[k]
                   - f_9 * isi1_563[k]
                   + f_3 * pc_y[k] * isk_726[k];

        t_911[k] = f_6 * isi0_564[k]
                   - f_7 * isi1_564[k]
                   + f_3 * pc_y[k] * isk_727[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pa_x, pc_x, pc_y, hsl0_914, hsk_734, hsl1_914, \
                         isi0_565, isi1_565, isk_728, isk_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_4 * isi0_565[k]
                   - f_5 * isi1_565[k]
                   + f_3 * pc_y[k] * isk_728[k];

        t_913[k] = f_3 * pc_y[k] * isk_729[k];

        t_914[k] = pa_x[k] * hsl0_914[k]
                   + f_18 * hsk_734[k]
                   - f_14 * pc_x[k] * hsl1_914[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_y, isi0_566, isi0_567, isi0_568, isi1_566, \
                         isi1_567, isi1_568, isk_730, isk_731, \
                         isk_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_10 * isi0_566[k]
                   - f_11 * isi1_566[k]
                   + f_3 * pc_y[k] * isk_730[k];

        t_916[k] = f_8 * isi0_567[k]
                   - f_9 * isi1_567[k]
                   + f_3 * pc_y[k] * isk_731[k];

        t_917[k] = f_6 * isi0_568[k]
                   - f_7 * isi1_568[k]
                   + f_3 * pc_y[k] * isk_732[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pa_x, pc_x, pc_y, hsl0_920, hsk_740, hsl1_920, \
                         isi0_569, isi1_569, isk_733, isk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_4 * isi0_569[k]
                   - f_5 * isi1_569[k]
                   + f_3 * pc_y[k] * isk_733[k];

        t_919[k] = f_3 * pc_y[k] * isk_734[k];

        t_920[k] = pa_x[k] * hsl0_920[k]
                   + f_17 * hsk_740[k]
                   - f_14 * pc_x[k] * hsl1_920[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pc_y, isi0_570, isi0_571, isi0_572, isi1_570, \
                         isi1_571, isi1_572, isk_735, isk_736, \
                         isk_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_12 * isi0_570[k]
                   - f_13 * isi1_570[k]
                   + f_3 * pc_y[k] * isk_735[k];

        t_922[k] = f_10 * isi0_571[k]
                   - f_11 * isi1_571[k]
                   + f_3 * pc_y[k] * isk_736[k];

        t_923[k] = f_8 * isi0_572[k]
                   - f_9 * isi1_572[k]
                   + f_3 * pc_y[k] * isk_737[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_y, isi0_573, isi0_574, isi1_573, isi1_574, \
                         isk_738, isk_739, isk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_6 * isi0_573[k]
                   - f_7 * isi1_573[k]
                   + f_3 * pc_y[k] * isk_738[k];

        t_925[k] = f_4 * isi0_574[k]
                   - f_5 * isi1_574[k]
                   + f_3 * pc_y[k] * isk_739[k];

        t_926[k] = f_3 * pc_y[k] * isk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pa_x, pc_x, hsl0_927, hsk_747, hsk_748, \
                         hsk_749, hsk_750, hsl1_927, isk_748, isk_749, \
                         isk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = pa_x[k] * hsl0_927[k]
                   + f_16 * hsk_747[k]
                   - f_14 * pc_x[k] * hsl1_927[k];

        t_928[k] = f_15 * hsk_748[k]
                   + f_3 * pc_x[k] * isk_748[k];

        t_929[k] = f_15 * hsk_749[k]
                   + f_3 * pc_x[k] * isk_749[k];

        t_930[k] = f_15 * hsk_750[k]
                   + f_3 * pc_x[k] * isk_750[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pc_x, pc_y, hsk_751, hsk_752, \
                         hsk_753, hsk_755, isk_747, isk_751, isk_752, isk_753, \
                         isk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_15 * hsk_751[k]
                   + f_3 * pc_x[k] * isk_751[k];

        t_932[k] = f_15 * hsk_752[k]
                   + f_3 * pc_x[k] * isk_752[k];

        t_933[k] = f_15 * hsk_753[k]
                   + f_3 * pc_x[k] * isk_753[k];

        t_934[k] = f_3 * pc_y[k] * isk_747[k];

        t_935[k] = f_15 * hsk_755[k]
                   + f_3 * pc_x[k] * isk_755[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, pa_x, pc_x, hsl0_936, hsl0_937, hsl0_938, \
                         hsl0_939, hsl1_936, hsl1_937, hsl1_938, \
                         hsl1_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = pa_x[k] * hsl0_936[k]
                   - f_14 * pc_x[k] * hsl1_936[k];

        t_937[k] = pa_x[k] * hsl0_937[k]
                   - f_14 * pc_x[k] * hsl1_937[k];

        t_938[k] = pa_x[k] * hsl0_938[k]
                   - f_14 * pc_x[k] * hsl1_938[k];

        t_939[k] = pa_x[k] * hsl0_939[k]
                   - f_14 * pc_x[k] * hsl1_939[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsl0_675 = buffer.data(hsl0 + 675);
    const auto *hsl0_676 = buffer.data(hsl0 + 676);
    const auto *hsl0_678 = buffer.data(hsl0 + 678);
    const auto *hsl0_681 = buffer.data(hsl0 + 681);
    const auto *hsl0_685 = buffer.data(hsl0 + 685);
    const auto *hsl0_690 = buffer.data(hsl0 + 690);
    const auto *hsl0_696 = buffer.data(hsl0 + 696);
    const auto *hsl0_711 = buffer.data(hsl0 + 711);
    const auto *hsl0_713 = buffer.data(hsl0 + 713);
    const auto *hsl0_714 = buffer.data(hsl0 + 714);
    const auto *hsl0_715 = buffer.data(hsl0 + 715);
    const auto *hsl0_716 = buffer.data(hsl0 + 716);
    const auto *hsl0_717 = buffer.data(hsl0 + 717);
    const auto *hsl0_940 = buffer.data(hsl0 + 940);
    const auto *hsl0_941 = buffer.data(hsl0 + 941);
    const auto *hsl0_942 = buffer.data(hsl0 + 942);
    const auto *hsl0_944 = buffer.data(hsl0 + 944);

    const auto *hsk_568 = buffer.data(hsk + 568);
    const auto *hsk_569 = buffer.data(hsk + 569);
    const auto *hsk_570 = buffer.data(hsk + 570);
    const auto *hsk_571 = buffer.data(hsk + 571);
    const auto *hsk_572 = buffer.data(hsk + 572);
    const auto *hsk_573 = buffer.data(hsk + 573);
    const auto *hsk_575 = buffer.data(hsk + 575);
    const auto *hsk_611 = buffer.data(hsk + 611);

    const auto *hsl1_675 = buffer.data(hsl1 + 675);
    const auto *hsl1_676 = buffer.data(hsl1 + 676);
    const auto *hsl1_678 = buffer.data(hsl1 + 678);
    const auto *hsl1_681 = buffer.data(hsl1 + 681);
    const auto *hsl1_685 = buffer.data(hsl1 + 685);
    const auto *hsl1_690 = buffer.data(hsl1 + 690);
    const auto *hsl1_696 = buffer.data(hsl1 + 696);
    const auto *hsl1_711 = buffer.data(hsl1 + 711);
    const auto *hsl1_713 = buffer.data(hsl1 + 713);
    const auto *hsl1_714 = buffer.data(hsl1 + 714);
    const auto *hsl1_715 = buffer.data(hsl1 + 715);
    const auto *hsl1_716 = buffer.data(hsl1 + 716);
    const auto *hsl1_717 = buffer.data(hsl1 + 717);
    const auto *hsl1_940 = buffer.data(hsl1 + 940);
    const auto *hsl1_941 = buffer.data(hsl1 + 941);
    const auto *hsl1_942 = buffer.data(hsl1 + 942);
    const auto *hsl1_944 = buffer.data(hsl1 + 944);

    const auto *isi0_588 = buffer.data(isi0 + 588);
    const auto *isi0_589 = buffer.data(isi0 + 589);
    const auto *isi0_591 = buffer.data(isi0 + 591);
    const auto *isi0_593 = buffer.data(isi0 + 593);
    const auto *isi0_594 = buffer.data(isi0 + 594);
    const auto *isi0_596 = buffer.data(isi0 + 596);
    const auto *isi0_597 = buffer.data(isi0 + 597);
    const auto *isi0_598 = buffer.data(isi0 + 598);
    const auto *isi0_600 = buffer.data(isi0 + 600);
    const auto *isi0_601 = buffer.data(isi0 + 601);
    const auto *isi0_602 = buffer.data(isi0 + 602);
    const auto *isi0_603 = buffer.data(isi0 + 603);
    const auto *isi0_605 = buffer.data(isi0 + 605);
    const auto *isi0_606 = buffer.data(isi0 + 606);
    const auto *isi0_607 = buffer.data(isi0 + 607);
    const auto *isi0_608 = buffer.data(isi0 + 608);
    const auto *isi0_609 = buffer.data(isi0 + 609);
    const auto *isi0_610 = buffer.data(isi0 + 610);
    const auto *isi0_611 = buffer.data(isi0 + 611);
    const auto *isi0_612 = buffer.data(isi0 + 612);
    const auto *isi0_613 = buffer.data(isi0 + 613);
    const auto *isi0_614 = buffer.data(isi0 + 614);
    const auto *isi0_615 = buffer.data(isi0 + 615);
    const auto *isi0_618 = buffer.data(isi0 + 618);
    const auto *isi0_620 = buffer.data(isi0 + 620);
    const auto *isi0_621 = buffer.data(isi0 + 621);
    const auto *isi0_623 = buffer.data(isi0 + 623);
    const auto *isi0_624 = buffer.data(isi0 + 624);
    const auto *isi0_625 = buffer.data(isi0 + 625);
    const auto *isi0_627 = buffer.data(isi0 + 627);
    const auto *isi0_628 = buffer.data(isi0 + 628);
    const auto *isi0_629 = buffer.data(isi0 + 629);
    const auto *isi0_630 = buffer.data(isi0 + 630);
    const auto *isi0_632 = buffer.data(isi0 + 632);
    const auto *isi0_633 = buffer.data(isi0 + 633);
    const auto *isi0_634 = buffer.data(isi0 + 634);
    const auto *isi0_635 = buffer.data(isi0 + 635);
    const auto *isi0_636 = buffer.data(isi0 + 636);
    const auto *isi0_638 = buffer.data(isi0 + 638);
    const auto *isi0_639 = buffer.data(isi0 + 639);
    const auto *isi0_640 = buffer.data(isi0 + 640);
    const auto *isi0_641 = buffer.data(isi0 + 641);
    const auto *isi0_642 = buffer.data(isi0 + 642);
    const auto *isi0_643 = buffer.data(isi0 + 643);
    const auto *isi0_644 = buffer.data(isi0 + 644);
    const auto *isi0_645 = buffer.data(isi0 + 645);
    const auto *isi0_646 = buffer.data(isi0 + 646);
    const auto *isi0_647 = buffer.data(isi0 + 647);
    const auto *isi0_648 = buffer.data(isi0 + 648);
    const auto *isi0_649 = buffer.data(isi0 + 649);
    const auto *isi0_650 = buffer.data(isi0 + 650);
    const auto *isi0_651 = buffer.data(isi0 + 651);
    const auto *isi0_652 = buffer.data(isi0 + 652);
    const auto *isi0_653 = buffer.data(isi0 + 653);
    const auto *isi0_654 = buffer.data(isi0 + 654);
    const auto *isi0_655 = buffer.data(isi0 + 655);
    const auto *isi0_656 = buffer.data(isi0 + 656);
    const auto *isi0_657 = buffer.data(isi0 + 657);
    const auto *isi0_658 = buffer.data(isi0 + 658);
    const auto *isi0_659 = buffer.data(isi0 + 659);
    const auto *isi0_660 = buffer.data(isi0 + 660);
    const auto *isi0_661 = buffer.data(isi0 + 661);
    const auto *isi0_662 = buffer.data(isi0 + 662);
    const auto *isi0_663 = buffer.data(isi0 + 663);
    const auto *isi0_664 = buffer.data(isi0 + 664);
    const auto *isi0_665 = buffer.data(isi0 + 665);
    const auto *isi0_666 = buffer.data(isi0 + 666);
    const auto *isi0_667 = buffer.data(isi0 + 667);
    const auto *isi0_668 = buffer.data(isi0 + 668);

    const auto *isi1_588 = buffer.data(isi1 + 588);
    const auto *isi1_589 = buffer.data(isi1 + 589);
    const auto *isi1_591 = buffer.data(isi1 + 591);
    const auto *isi1_593 = buffer.data(isi1 + 593);
    const auto *isi1_594 = buffer.data(isi1 + 594);
    const auto *isi1_596 = buffer.data(isi1 + 596);
    const auto *isi1_597 = buffer.data(isi1 + 597);
    const auto *isi1_598 = buffer.data(isi1 + 598);
    const auto *isi1_600 = buffer.data(isi1 + 600);
    const auto *isi1_601 = buffer.data(isi1 + 601);
    const auto *isi1_602 = buffer.data(isi1 + 602);
    const auto *isi1_603 = buffer.data(isi1 + 603);
    const auto *isi1_605 = buffer.data(isi1 + 605);
    const auto *isi1_606 = buffer.data(isi1 + 606);
    const auto *isi1_607 = buffer.data(isi1 + 607);
    const auto *isi1_608 = buffer.data(isi1 + 608);
    const auto *isi1_609 = buffer.data(isi1 + 609);
    const auto *isi1_610 = buffer.data(isi1 + 610);
    const auto *isi1_611 = buffer.data(isi1 + 611);
    const auto *isi1_612 = buffer.data(isi1 + 612);
    const auto *isi1_613 = buffer.data(isi1 + 613);
    const auto *isi1_614 = buffer.data(isi1 + 614);
    const auto *isi1_615 = buffer.data(isi1 + 615);
    const auto *isi1_618 = buffer.data(isi1 + 618);
    const auto *isi1_620 = buffer.data(isi1 + 620);
    const auto *isi1_621 = buffer.data(isi1 + 621);
    const auto *isi1_623 = buffer.data(isi1 + 623);
    const auto *isi1_624 = buffer.data(isi1 + 624);
    const auto *isi1_625 = buffer.data(isi1 + 625);
    const auto *isi1_627 = buffer.data(isi1 + 627);
    const auto *isi1_628 = buffer.data(isi1 + 628);
    const auto *isi1_629 = buffer.data(isi1 + 629);
    const auto *isi1_630 = buffer.data(isi1 + 630);
    const auto *isi1_632 = buffer.data(isi1 + 632);
    const auto *isi1_633 = buffer.data(isi1 + 633);
    const auto *isi1_634 = buffer.data(isi1 + 634);
    const auto *isi1_635 = buffer.data(isi1 + 635);
    const auto *isi1_636 = buffer.data(isi1 + 636);
    const auto *isi1_638 = buffer.data(isi1 + 638);
    const auto *isi1_639 = buffer.data(isi1 + 639);
    const auto *isi1_640 = buffer.data(isi1 + 640);
    const auto *isi1_641 = buffer.data(isi1 + 641);
    const auto *isi1_642 = buffer.data(isi1 + 642);
    const auto *isi1_643 = buffer.data(isi1 + 643);
    const auto *isi1_644 = buffer.data(isi1 + 644);
    const auto *isi1_645 = buffer.data(isi1 + 645);
    const auto *isi1_646 = buffer.data(isi1 + 646);
    const auto *isi1_647 = buffer.data(isi1 + 647);
    const auto *isi1_648 = buffer.data(isi1 + 648);
    const auto *isi1_649 = buffer.data(isi1 + 649);
    const auto *isi1_650 = buffer.data(isi1 + 650);
    const auto *isi1_651 = buffer.data(isi1 + 651);
    const auto *isi1_652 = buffer.data(isi1 + 652);
    const auto *isi1_653 = buffer.data(isi1 + 653);
    const auto *isi1_654 = buffer.data(isi1 + 654);
    const auto *isi1_655 = buffer.data(isi1 + 655);
    const auto *isi1_656 = buffer.data(isi1 + 656);
    const auto *isi1_657 = buffer.data(isi1 + 657);
    const auto *isi1_658 = buffer.data(isi1 + 658);
    const auto *isi1_659 = buffer.data(isi1 + 659);
    const auto *isi1_660 = buffer.data(isi1 + 660);
    const auto *isi1_661 = buffer.data(isi1 + 661);
    const auto *isi1_662 = buffer.data(isi1 + 662);
    const auto *isi1_663 = buffer.data(isi1 + 663);
    const auto *isi1_664 = buffer.data(isi1 + 664);
    const auto *isi1_665 = buffer.data(isi1 + 665);
    const auto *isi1_666 = buffer.data(isi1 + 666);
    const auto *isi1_667 = buffer.data(isi1 + 667);
    const auto *isi1_668 = buffer.data(isi1 + 668);

    const auto *isk_755 = buffer.data(isk + 755);
    const auto *isk_756 = buffer.data(isk + 756);
    const auto *isk_757 = buffer.data(isk + 757);
    const auto *isk_759 = buffer.data(isk + 759);
    const auto *isk_761 = buffer.data(isk + 761);
    const auto *isk_762 = buffer.data(isk + 762);
    const auto *isk_764 = buffer.data(isk + 764);
    const auto *isk_765 = buffer.data(isk + 765);
    const auto *isk_766 = buffer.data(isk + 766);
    const auto *isk_768 = buffer.data(isk + 768);
    const auto *isk_769 = buffer.data(isk + 769);
    const auto *isk_770 = buffer.data(isk + 770);
    const auto *isk_771 = buffer.data(isk + 771);
    const auto *isk_773 = buffer.data(isk + 773);
    const auto *isk_774 = buffer.data(isk + 774);
    const auto *isk_775 = buffer.data(isk + 775);
    const auto *isk_776 = buffer.data(isk + 776);
    const auto *isk_777 = buffer.data(isk + 777);
    const auto *isk_779 = buffer.data(isk + 779);
    const auto *isk_780 = buffer.data(isk + 780);
    const auto *isk_781 = buffer.data(isk + 781);
    const auto *isk_782 = buffer.data(isk + 782);
    const auto *isk_783 = buffer.data(isk + 783);
    const auto *isk_784 = buffer.data(isk + 784);
    const auto *isk_785 = buffer.data(isk + 785);
    const auto *isk_786 = buffer.data(isk + 786);
    const auto *isk_787 = buffer.data(isk + 787);
    const auto *isk_788 = buffer.data(isk + 788);
    const auto *isk_789 = buffer.data(isk + 789);
    const auto *isk_790 = buffer.data(isk + 790);
    const auto *isk_791 = buffer.data(isk + 791);
    const auto *isk_794 = buffer.data(isk + 794);
    const auto *isk_796 = buffer.data(isk + 796);
    const auto *isk_797 = buffer.data(isk + 797);
    const auto *isk_799 = buffer.data(isk + 799);
    const auto *isk_800 = buffer.data(isk + 800);
    const auto *isk_801 = buffer.data(isk + 801);
    const auto *isk_803 = buffer.data(isk + 803);
    const auto *isk_804 = buffer.data(isk + 804);
    const auto *isk_805 = buffer.data(isk + 805);
    const auto *isk_806 = buffer.data(isk + 806);
    const auto *isk_808 = buffer.data(isk + 808);
    const auto *isk_809 = buffer.data(isk + 809);
    const auto *isk_810 = buffer.data(isk + 810);
    const auto *isk_811 = buffer.data(isk + 811);
    const auto *isk_812 = buffer.data(isk + 812);
    const auto *isk_814 = buffer.data(isk + 814);
    const auto *isk_815 = buffer.data(isk + 815);
    const auto *isk_816 = buffer.data(isk + 816);
    const auto *isk_817 = buffer.data(isk + 817);
    const auto *isk_818 = buffer.data(isk + 818);
    const auto *isk_819 = buffer.data(isk + 819);
    const auto *isk_820 = buffer.data(isk + 820);
    const auto *isk_821 = buffer.data(isk + 821);
    const auto *isk_822 = buffer.data(isk + 822);
    const auto *isk_823 = buffer.data(isk + 823);
    const auto *isk_824 = buffer.data(isk + 824);
    const auto *isk_825 = buffer.data(isk + 825);
    const auto *isk_826 = buffer.data(isk + 826);
    const auto *isk_827 = buffer.data(isk + 827);
    const auto *isk_828 = buffer.data(isk + 828);
    const auto *isk_829 = buffer.data(isk + 829);
    const auto *isk_830 = buffer.data(isk + 830);
    const auto *isk_831 = buffer.data(isk + 831);
    const auto *isk_832 = buffer.data(isk + 832);
    const auto *isk_833 = buffer.data(isk + 833);
    const auto *isk_834 = buffer.data(isk + 834);
    const auto *isk_835 = buffer.data(isk + 835);
    const auto *isk_836 = buffer.data(isk + 836);
    const auto *isk_837 = buffer.data(isk + 837);
    const auto *isk_838 = buffer.data(isk + 838);
    const auto *isk_839 = buffer.data(isk + 839);
    const auto *isk_840 = buffer.data(isk + 840);
    const auto *isk_841 = buffer.data(isk + 841);
    const auto *isk_842 = buffer.data(isk + 842);
    const auto *isk_843 = buffer.data(isk + 843);
    const auto *isk_844 = buffer.data(isk + 844);
    const auto *isk_845 = buffer.data(isk + 845);
    const auto *isk_846 = buffer.data(isk + 846);
    const auto *isk_847 = buffer.data(isk + 847);
    const auto *isk_848 = buffer.data(isk + 848);
    const auto *isk_849 = buffer.data(isk + 849);
    const auto *isk_850 = buffer.data(isk + 850);
    const auto *isk_851 = buffer.data(isk + 851);
    const auto *isk_852 = buffer.data(isk + 852);

#pragma omp simd aligned(t_940, t_941, t_942, t_943, pa_x, pc_x, pc_y, hsl0_940, hsl0_941, \
                         hsl0_942, hsl1_940, hsl1_941, hsl1_942, \
                         isk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = pa_x[k] * hsl0_940[k]
                   - f_14 * pc_x[k] * hsl1_940[k];

        t_941[k] = pa_x[k] * hsl0_941[k]
                   - f_14 * pc_x[k] * hsl1_941[k];

        t_942[k] = pa_x[k] * hsl0_942[k]
                   - f_14 * pc_x[k] * hsl1_942[k];

        t_943[k] = f_3 * pc_y[k] * isk_755[k];
    }

#pragma omp simd aligned(t_944, t_945, t_946, t_947, pa_x, pc_x, pc_z, hsl0_944, hsl1_944, \
                         isi0_588, isi0_589, isi1_588, isi1_589, isk_756, \
                         isk_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = pa_x[k] * hsl0_944[k]
                   - f_14 * pc_x[k] * hsl1_944[k];

        t_945[k] = f_1 * isi0_588[k]
                   - f_2 * isi1_588[k]
                   + f_3 * pc_x[k] * isk_756[k];

        t_946[k] = f_20 * isi0_589[k]
                   - f_21 * isi1_589[k]
                   + f_3 * pc_x[k] * isk_757[k];

        t_947[k] = f_3 * pc_z[k] * isk_756[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, t_951, pc_x, pc_z, isi0_591, isi0_593, isi0_594, \
                         isi1_591, isi1_593, isi1_594, isk_757, isk_759, isk_761, \
                         isk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = f_12 * isi0_591[k]
                   - f_13 * isi1_591[k]
                   + f_3 * pc_x[k] * isk_759[k];

        t_949[k] = f_3 * pc_z[k] * isk_757[k];

        t_950[k] = f_12 * isi0_593[k]
                   - f_13 * isi1_593[k]
                   + f_3 * pc_x[k] * isk_761[k];

        t_951[k] = f_10 * isi0_594[k]
                   - f_11 * isi1_594[k]
                   + f_3 * pc_x[k] * isk_762[k];
    }

#pragma omp simd aligned(t_952, t_953, t_954, t_955, pc_x, pc_z, isi0_596, isi0_597, isi0_598, \
                         isi1_596, isi1_597, isi1_598, isk_759, isk_764, isk_765, \
                         isk_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_952[k] = f_3 * pc_z[k] * isk_759[k];

        t_953[k] = f_10 * isi0_596[k]
                   - f_11 * isi1_596[k]
                   + f_3 * pc_x[k] * isk_764[k];

        t_954[k] = f_10 * isi0_597[k]
                   - f_11 * isi1_597[k]
                   + f_3 * pc_x[k] * isk_765[k];

        t_955[k] = f_8 * isi0_598[k]
                   - f_9 * isi1_598[k]
                   + f_3 * pc_x[k] * isk_766[k];
    }

#pragma omp simd aligned(t_956, t_957, t_958, t_959, pc_x, pc_z, isi0_600, isi0_601, isi0_602, \
                         isi1_600, isi1_601, isi1_602, isk_762, isk_768, isk_769, \
                         isk_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_956[k] = f_3 * pc_z[k] * isk_762[k];

        t_957[k] = f_8 * isi0_600[k]
                   - f_9 * isi1_600[k]
                   + f_3 * pc_x[k] * isk_768[k];

        t_958[k] = f_8 * isi0_601[k]
                   - f_9 * isi1_601[k]
                   + f_3 * pc_x[k] * isk_769[k];

        t_959[k] = f_8 * isi0_602[k]
                   - f_9 * isi1_602[k]
                   + f_3 * pc_x[k] * isk_770[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, pc_x, pc_z, isi0_603, isi0_605, isi0_606, \
                         isi1_603, isi1_605, isi1_606, isk_766, isk_771, isk_773, \
                         isk_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_6 * isi0_603[k]
                   - f_7 * isi1_603[k]
                   + f_3 * pc_x[k] * isk_771[k];

        t_961[k] = f_3 * pc_z[k] * isk_766[k];

        t_962[k] = f_6 * isi0_605[k]
                   - f_7 * isi1_605[k]
                   + f_3 * pc_x[k] * isk_773[k];

        t_963[k] = f_6 * isi0_606[k]
                   - f_7 * isi1_606[k]
                   + f_3 * pc_x[k] * isk_774[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pc_x, pc_z, isi0_607, isi0_608, isi0_609, \
                         isi1_607, isi1_608, isi1_609, isk_771, isk_775, isk_776, \
                         isk_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = f_6 * isi0_607[k]
                   - f_7 * isi1_607[k]
                   + f_3 * pc_x[k] * isk_775[k];

        t_965[k] = f_6 * isi0_608[k]
                   - f_7 * isi1_608[k]
                   + f_3 * pc_x[k] * isk_776[k];

        t_966[k] = f_4 * isi0_609[k]
                   - f_5 * isi1_609[k]
                   + f_3 * pc_x[k] * isk_777[k];

        t_967[k] = f_3 * pc_z[k] * isk_771[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, pc_x, isi0_611, isi0_612, isi0_613, isi1_611, \
                         isi1_612, isi1_613, isk_779, isk_780, \
                         isk_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_4 * isi0_611[k]
                   - f_5 * isi1_611[k]
                   + f_3 * pc_x[k] * isk_779[k];

        t_969[k] = f_4 * isi0_612[k]
                   - f_5 * isi1_612[k]
                   + f_3 * pc_x[k] * isk_780[k];

        t_970[k] = f_4 * isi0_613[k]
                   - f_5 * isi1_613[k]
                   + f_3 * pc_x[k] * isk_781[k];
    }

#pragma omp simd aligned(t_971, t_972, t_973, t_974, t_975, pc_x, isi0_614, isi0_615, \
                         isi1_614, isi1_615, isk_782, isk_783, isk_784, isk_785, \
                         isk_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_971[k] = f_4 * isi0_614[k]
                   - f_5 * isi1_614[k]
                   + f_3 * pc_x[k] * isk_782[k];

        t_972[k] = f_4 * isi0_615[k]
                   - f_5 * isi1_615[k]
                   + f_3 * pc_x[k] * isk_783[k];

        t_973[k] = f_3 * pc_x[k] * isk_784[k];

        t_974[k] = f_3 * pc_x[k] * isk_785[k];

        t_975[k] = f_3 * pc_x[k] * isk_786[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, t_979, t_980, pc_x, isk_787, isk_788, isk_789, \
                         isk_790, isk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_3 * pc_x[k] * isk_787[k];

        t_977[k] = f_3 * pc_x[k] * isk_788[k];

        t_978[k] = f_3 * pc_x[k] * isk_789[k];

        t_979[k] = f_3 * pc_x[k] * isk_790[k];

        t_980[k] = f_3 * pc_x[k] * isk_791[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, t_984, pc_y, pc_z, hsk_568, isi0_609, isi0_610, \
                         isi1_609, isi1_610, isk_784, isk_785, \
                         isk_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_0 * hsk_568[k]
                   + f_1 * isi0_609[k]
                   - f_2 * isi1_609[k]
                   + f_3 * pc_y[k] * isk_784[k];

        t_982[k] = f_3 * pc_z[k] * isk_784[k];

        t_983[k] = f_4 * isi0_609[k]
                   - f_5 * isi1_609[k]
                   + f_3 * pc_z[k] * isk_785[k];

        t_984[k] = f_6 * isi0_610[k]
                   - f_7 * isi1_610[k]
                   + f_3 * pc_z[k] * isk_786[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pc_z, isi0_611, isi0_612, isi0_613, isi1_611, \
                         isi1_612, isi1_613, isk_787, isk_788, \
                         isk_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_8 * isi0_611[k]
                   - f_9 * isi1_611[k]
                   + f_3 * pc_z[k] * isk_787[k];

        t_986[k] = f_10 * isi0_612[k]
                   - f_11 * isi1_612[k]
                   + f_3 * pc_z[k] * isk_788[k];

        t_987[k] = f_12 * isi0_613[k]
                   - f_13 * isi1_613[k]
                   + f_3 * pc_z[k] * isk_789[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, pa_z, pc_y, pc_z, hsl0_675, hsl0_676, \
                         hsk_575, hsl1_675, hsl1_676, isi0_615, isi1_615, \
                         isk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_0 * hsk_575[k]
                   + f_3 * pc_y[k] * isk_791[k];

        t_989[k] = f_1 * isi0_615[k]
                   - f_2 * isi1_615[k]
                   + f_3 * pc_z[k] * isk_791[k];

        t_990[k] = pa_z[k] * hsl0_675[k]
                   - f_14 * pc_z[k] * hsl1_675[k];

        t_991[k] = pa_z[k] * hsl0_676[k]
                   - f_14 * pc_z[k] * hsl1_676[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, pa_z, pc_x, pc_z, hsl0_678, hsl1_678, isi0_618, \
                         isi0_620, isi1_618, isi1_620, isk_794, \
                         isk_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_20 * isi0_618[k]
                   - f_21 * isi1_618[k]
                   + f_3 * pc_x[k] * isk_794[k];

        t_993[k] = pa_z[k] * hsl0_678[k]
                   - f_14 * pc_z[k] * hsl1_678[k];

        t_994[k] = f_12 * isi0_620[k]
                   - f_13 * isi1_620[k]
                   + f_3 * pc_x[k] * isk_796[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, pa_z, pc_x, pc_z, hsl0_681, hsl1_681, isi0_621, \
                         isi0_623, isi1_621, isi1_623, isk_797, \
                         isk_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_12 * isi0_621[k]
                   - f_13 * isi1_621[k]
                   + f_3 * pc_x[k] * isk_797[k];

        t_996[k] = pa_z[k] * hsl0_681[k]
                   - f_14 * pc_z[k] * hsl1_681[k];

        t_997[k] = f_10 * isi0_623[k]
                   - f_11 * isi1_623[k]
                   + f_3 * pc_x[k] * isk_799[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, pa_z, pc_x, pc_z, hsl0_685, hsl1_685, isi0_624, \
                         isi0_625, isi1_624, isi1_625, isk_800, \
                         isk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_10 * isi0_624[k]
                   - f_11 * isi1_624[k]
                   + f_3 * pc_x[k] * isk_800[k];

        t_999[k] = f_10 * isi0_625[k]
                   - f_11 * isi1_625[k]
                   + f_3 * pc_x[k] * isk_801[k];

        t_1000[k] = pa_z[k] * hsl0_685[k]
                    - f_14 * pc_z[k] * hsl1_685[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, pc_x, isi0_627, isi0_628, isi0_629, isi1_627, \
                         isi1_628, isi1_629, isk_803, isk_804, \
                         isk_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_8 * isi0_627[k]
                    - f_9 * isi1_627[k]
                    + f_3 * pc_x[k] * isk_803[k];

        t_1002[k] = f_8 * isi0_628[k]
                    - f_9 * isi1_628[k]
                    + f_3 * pc_x[k] * isk_804[k];

        t_1003[k] = f_8 * isi0_629[k]
                    - f_9 * isi1_629[k]
                    + f_3 * pc_x[k] * isk_805[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, pa_z, pc_x, pc_z, hsl0_690, hsl1_690, \
                         isi0_630, isi0_632, isi1_630, isi1_632, isk_806, \
                         isk_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_8 * isi0_630[k]
                    - f_9 * isi1_630[k]
                    + f_3 * pc_x[k] * isk_806[k];

        t_1005[k] = pa_z[k] * hsl0_690[k]
                    - f_14 * pc_z[k] * hsl1_690[k];

        t_1006[k] = f_6 * isi0_632[k]
                    - f_7 * isi1_632[k]
                    + f_3 * pc_x[k] * isk_808[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, pc_x, isi0_633, isi0_634, isi0_635, isi1_633, \
                         isi1_634, isi1_635, isk_809, isk_810, \
                         isk_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = f_6 * isi0_633[k]
                    - f_7 * isi1_633[k]
                    + f_3 * pc_x[k] * isk_809[k];

        t_1008[k] = f_6 * isi0_634[k]
                    - f_7 * isi1_634[k]
                    + f_3 * pc_x[k] * isk_810[k];

        t_1009[k] = f_6 * isi0_635[k]
                    - f_7 * isi1_635[k]
                    + f_3 * pc_x[k] * isk_811[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, pa_z, pc_x, pc_z, hsl0_696, hsl1_696, \
                         isi0_636, isi0_638, isi1_636, isi1_638, isk_812, \
                         isk_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_6 * isi0_636[k]
                    - f_7 * isi1_636[k]
                    + f_3 * pc_x[k] * isk_812[k];

        t_1011[k] = pa_z[k] * hsl0_696[k]
                    - f_14 * pc_z[k] * hsl1_696[k];

        t_1012[k] = f_4 * isi0_638[k]
                    - f_5 * isi1_638[k]
                    + f_3 * pc_x[k] * isk_814[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pc_x, isi0_639, isi0_640, isi0_641, isi1_639, \
                         isi1_640, isi1_641, isk_815, isk_816, \
                         isk_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_4 * isi0_639[k]
                    - f_5 * isi1_639[k]
                    + f_3 * pc_x[k] * isk_815[k];

        t_1014[k] = f_4 * isi0_640[k]
                    - f_5 * isi1_640[k]
                    + f_3 * pc_x[k] * isk_816[k];

        t_1015[k] = f_4 * isi0_641[k]
                    - f_5 * isi1_641[k]
                    + f_3 * pc_x[k] * isk_817[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, t_1020, pc_x, isi0_642, isi0_643, \
                         isi1_642, isi1_643, isk_818, isk_819, isk_820, isk_821, \
                         isk_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_4 * isi0_642[k]
                    - f_5 * isi1_642[k]
                    + f_3 * pc_x[k] * isk_818[k];

        t_1017[k] = f_4 * isi0_643[k]
                    - f_5 * isi1_643[k]
                    + f_3 * pc_x[k] * isk_819[k];

        t_1018[k] = f_3 * pc_x[k] * isk_820[k];

        t_1019[k] = f_3 * pc_x[k] * isk_821[k];

        t_1020[k] = f_3 * pc_x[k] * isk_822[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, t_1025, t_1026, pa_z, pc_x, pc_z, \
                         hsl0_711, hsl1_711, isk_823, isk_824, isk_825, isk_826, \
                         isk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_3 * pc_x[k] * isk_823[k];

        t_1022[k] = f_3 * pc_x[k] * isk_824[k];

        t_1023[k] = f_3 * pc_x[k] * isk_825[k];

        t_1024[k] = f_3 * pc_x[k] * isk_826[k];

        t_1025[k] = f_3 * pc_x[k] * isk_827[k];

        t_1026[k] = pa_z[k] * hsl0_711[k]
                    - f_14 * pc_z[k] * hsl1_711[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, pa_z, pc_z, hsl0_713, hsl0_714, hsk_568, \
                         hsk_569, hsk_570, hsl1_713, hsl1_714, \
                         isk_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_15 * hsk_568[k]
                    + f_3 * pc_z[k] * isk_820[k];

        t_1028[k] = pa_z[k] * hsl0_713[k]
                    + f_16 * hsk_569[k]
                    - f_14 * pc_z[k] * hsl1_713[k];

        t_1029[k] = pa_z[k] * hsl0_714[k]
                    + f_17 * hsk_570[k]
                    - f_14 * pc_z[k] * hsl1_714[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, pa_z, pc_z, hsl0_715, hsl0_716, hsl0_717, \
                         hsk_571, hsk_572, hsk_573, hsl1_715, hsl1_716, \
                         hsl1_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = pa_z[k] * hsl0_715[k]
                    + f_18 * hsk_571[k]
                    - f_14 * pc_z[k] * hsl1_715[k];

        t_1031[k] = pa_z[k] * hsl0_716[k]
                    + f_19 * hsk_572[k]
                    - f_14 * pc_z[k] * hsl1_716[k];

        t_1032[k] = pa_z[k] * hsl0_717[k]
                    + f_0 * hsk_573[k]
                    - f_14 * pc_z[k] * hsl1_717[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, pc_x, pc_y, pc_z, hsk_575, hsk_611, isi0_643, \
                         isi0_644, isi1_643, isi1_644, isk_827, \
                         isk_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_19 * hsk_611[k]
                    + f_3 * pc_y[k] * isk_827[k];

        t_1034[k] = f_15 * hsk_575[k]
                    + f_1 * isi0_643[k]
                    - f_2 * isi1_643[k]
                    + f_3 * pc_z[k] * isk_827[k];

        t_1035[k] = f_1 * isi0_644[k]
                    - f_2 * isi1_644[k]
                    + f_3 * pc_x[k] * isk_828[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, pc_x, isi0_645, isi0_646, isi0_647, isi1_645, \
                         isi1_646, isi1_647, isk_829, isk_830, \
                         isk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = f_20 * isi0_645[k]
                    - f_21 * isi1_645[k]
                    + f_3 * pc_x[k] * isk_829[k];

        t_1037[k] = f_20 * isi0_646[k]
                    - f_21 * isi1_646[k]
                    + f_3 * pc_x[k] * isk_830[k];

        t_1038[k] = f_12 * isi0_647[k]
                    - f_13 * isi1_647[k]
                    + f_3 * pc_x[k] * isk_831[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pc_x, isi0_648, isi0_649, isi0_650, isi1_648, \
                         isi1_649, isi1_650, isk_832, isk_833, \
                         isk_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_12 * isi0_648[k]
                    - f_13 * isi1_648[k]
                    + f_3 * pc_x[k] * isk_832[k];

        t_1040[k] = f_12 * isi0_649[k]
                    - f_13 * isi1_649[k]
                    + f_3 * pc_x[k] * isk_833[k];

        t_1041[k] = f_10 * isi0_650[k]
                    - f_11 * isi1_650[k]
                    + f_3 * pc_x[k] * isk_834[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, pc_x, isi0_651, isi0_652, isi0_653, isi1_651, \
                         isi1_652, isi1_653, isk_835, isk_836, \
                         isk_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_10 * isi0_651[k]
                    - f_11 * isi1_651[k]
                    + f_3 * pc_x[k] * isk_835[k];

        t_1043[k] = f_10 * isi0_652[k]
                    - f_11 * isi1_652[k]
                    + f_3 * pc_x[k] * isk_836[k];

        t_1044[k] = f_10 * isi0_653[k]
                    - f_11 * isi1_653[k]
                    + f_3 * pc_x[k] * isk_837[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, pc_x, isi0_654, isi0_655, isi0_656, isi1_654, \
                         isi1_655, isi1_656, isk_838, isk_839, \
                         isk_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = f_8 * isi0_654[k]
                    - f_9 * isi1_654[k]
                    + f_3 * pc_x[k] * isk_838[k];

        t_1046[k] = f_8 * isi0_655[k]
                    - f_9 * isi1_655[k]
                    + f_3 * pc_x[k] * isk_839[k];

        t_1047[k] = f_8 * isi0_656[k]
                    - f_9 * isi1_656[k]
                    + f_3 * pc_x[k] * isk_840[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, pc_x, isi0_657, isi0_658, isi0_659, isi1_657, \
                         isi1_658, isi1_659, isk_841, isk_842, \
                         isk_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_8 * isi0_657[k]
                    - f_9 * isi1_657[k]
                    + f_3 * pc_x[k] * isk_841[k];

        t_1049[k] = f_8 * isi0_658[k]
                    - f_9 * isi1_658[k]
                    + f_3 * pc_x[k] * isk_842[k];

        t_1050[k] = f_6 * isi0_659[k]
                    - f_7 * isi1_659[k]
                    + f_3 * pc_x[k] * isk_843[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, pc_x, isi0_660, isi0_661, isi0_662, isi1_660, \
                         isi1_661, isi1_662, isk_844, isk_845, \
                         isk_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_6 * isi0_660[k]
                    - f_7 * isi1_660[k]
                    + f_3 * pc_x[k] * isk_844[k];

        t_1052[k] = f_6 * isi0_661[k]
                    - f_7 * isi1_661[k]
                    + f_3 * pc_x[k] * isk_845[k];

        t_1053[k] = f_6 * isi0_662[k]
                    - f_7 * isi1_662[k]
                    + f_3 * pc_x[k] * isk_846[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, pc_x, isi0_663, isi0_664, isi0_665, isi1_663, \
                         isi1_664, isi1_665, isk_847, isk_848, \
                         isk_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_6 * isi0_663[k]
                    - f_7 * isi1_663[k]
                    + f_3 * pc_x[k] * isk_847[k];

        t_1055[k] = f_6 * isi0_664[k]
                    - f_7 * isi1_664[k]
                    + f_3 * pc_x[k] * isk_848[k];

        t_1056[k] = f_4 * isi0_665[k]
                    - f_5 * isi1_665[k]
                    + f_3 * pc_x[k] * isk_849[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, pc_x, isi0_666, isi0_667, isi0_668, isi1_666, \
                         isi1_667, isi1_668, isk_850, isk_851, \
                         isk_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_4 * isi0_666[k]
                    - f_5 * isi1_666[k]
                    + f_3 * pc_x[k] * isk_850[k];

        t_1058[k] = f_4 * isi0_667[k]
                    - f_5 * isi1_667[k]
                    + f_3 * pc_x[k] * isk_851[k];

        t_1059[k] = f_4 * isi0_668[k]
                    - f_5 * isi1_668[k]
                    + f_3 * pc_x[k] * isk_852[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsl0,
                                                          const size_t hsk, const size_t hsl1,
                                                          const size_t isi0, const size_t isi1,
                                                          const size_t isk, const size_t ncols,
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
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);

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
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsl0_900 = buffer.data(hsl0 + 900);
    const auto *hsl0_902 = buffer.data(hsl0 + 902);

    const auto *hsk_604 = buffer.data(hsk + 604);
    const auto *hsk_611 = buffer.data(hsk + 611);
    const auto *hsk_640 = buffer.data(hsk + 640);
    const auto *hsk_642 = buffer.data(hsk + 642);
    const auto *hsk_643 = buffer.data(hsk + 643);
    const auto *hsk_644 = buffer.data(hsk + 644);
    const auto *hsk_645 = buffer.data(hsk + 645);
    const auto *hsk_646 = buffer.data(hsk + 646);
    const auto *hsk_647 = buffer.data(hsk + 647);
    const auto *hsk_676 = buffer.data(hsk + 676);
    const auto *hsk_678 = buffer.data(hsk + 678);
    const auto *hsk_679 = buffer.data(hsk + 679);
    const auto *hsk_680 = buffer.data(hsk + 680);
    const auto *hsk_681 = buffer.data(hsk + 681);
    const auto *hsk_682 = buffer.data(hsk + 682);
    const auto *hsk_683 = buffer.data(hsk + 683);
    const auto *hsk_712 = buffer.data(hsk + 712);
    const auto *hsk_714 = buffer.data(hsk + 714);
    const auto *hsk_715 = buffer.data(hsk + 715);
    const auto *hsk_716 = buffer.data(hsk + 716);
    const auto *hsk_717 = buffer.data(hsk + 717);
    const auto *hsk_718 = buffer.data(hsk + 718);
    const auto *hsk_719 = buffer.data(hsk + 719);

    const auto *hsl1_900 = buffer.data(hsl1 + 900);
    const auto *hsl1_902 = buffer.data(hsl1 + 902);

    const auto *isi0_665 = buffer.data(isi0 + 665);
    const auto *isi0_667 = buffer.data(isi0 + 667);
    const auto *isi0_668 = buffer.data(isi0 + 668);
    const auto *isi0_669 = buffer.data(isi0 + 669);
    const auto *isi0_670 = buffer.data(isi0 + 670);
    const auto *isi0_671 = buffer.data(isi0 + 671);
    const auto *isi0_672 = buffer.data(isi0 + 672);
    const auto *isi0_673 = buffer.data(isi0 + 673);
    const auto *isi0_674 = buffer.data(isi0 + 674);
    const auto *isi0_675 = buffer.data(isi0 + 675);
    const auto *isi0_676 = buffer.data(isi0 + 676);
    const auto *isi0_677 = buffer.data(isi0 + 677);
    const auto *isi0_678 = buffer.data(isi0 + 678);
    const auto *isi0_679 = buffer.data(isi0 + 679);
    const auto *isi0_680 = buffer.data(isi0 + 680);
    const auto *isi0_681 = buffer.data(isi0 + 681);
    const auto *isi0_682 = buffer.data(isi0 + 682);
    const auto *isi0_683 = buffer.data(isi0 + 683);
    const auto *isi0_684 = buffer.data(isi0 + 684);
    const auto *isi0_685 = buffer.data(isi0 + 685);
    const auto *isi0_686 = buffer.data(isi0 + 686);
    const auto *isi0_687 = buffer.data(isi0 + 687);
    const auto *isi0_688 = buffer.data(isi0 + 688);
    const auto *isi0_689 = buffer.data(isi0 + 689);
    const auto *isi0_690 = buffer.data(isi0 + 690);
    const auto *isi0_691 = buffer.data(isi0 + 691);
    const auto *isi0_692 = buffer.data(isi0 + 692);
    const auto *isi0_693 = buffer.data(isi0 + 693);
    const auto *isi0_694 = buffer.data(isi0 + 694);
    const auto *isi0_695 = buffer.data(isi0 + 695);
    const auto *isi0_696 = buffer.data(isi0 + 696);
    const auto *isi0_697 = buffer.data(isi0 + 697);
    const auto *isi0_698 = buffer.data(isi0 + 698);
    const auto *isi0_699 = buffer.data(isi0 + 699);
    const auto *isi0_700 = buffer.data(isi0 + 700);
    const auto *isi0_701 = buffer.data(isi0 + 701);
    const auto *isi0_702 = buffer.data(isi0 + 702);
    const auto *isi0_703 = buffer.data(isi0 + 703);
    const auto *isi0_704 = buffer.data(isi0 + 704);
    const auto *isi0_705 = buffer.data(isi0 + 705);
    const auto *isi0_706 = buffer.data(isi0 + 706);
    const auto *isi0_707 = buffer.data(isi0 + 707);
    const auto *isi0_708 = buffer.data(isi0 + 708);
    const auto *isi0_709 = buffer.data(isi0 + 709);
    const auto *isi0_710 = buffer.data(isi0 + 710);
    const auto *isi0_711 = buffer.data(isi0 + 711);
    const auto *isi0_712 = buffer.data(isi0 + 712);
    const auto *isi0_713 = buffer.data(isi0 + 713);
    const auto *isi0_714 = buffer.data(isi0 + 714);
    const auto *isi0_715 = buffer.data(isi0 + 715);
    const auto *isi0_716 = buffer.data(isi0 + 716);
    const auto *isi0_717 = buffer.data(isi0 + 717);
    const auto *isi0_718 = buffer.data(isi0 + 718);
    const auto *isi0_719 = buffer.data(isi0 + 719);
    const auto *isi0_720 = buffer.data(isi0 + 720);
    const auto *isi0_721 = buffer.data(isi0 + 721);
    const auto *isi0_722 = buffer.data(isi0 + 722);
    const auto *isi0_723 = buffer.data(isi0 + 723);
    const auto *isi0_724 = buffer.data(isi0 + 724);
    const auto *isi0_725 = buffer.data(isi0 + 725);
    const auto *isi0_726 = buffer.data(isi0 + 726);
    const auto *isi0_727 = buffer.data(isi0 + 727);
    const auto *isi0_729 = buffer.data(isi0 + 729);
    const auto *isi0_731 = buffer.data(isi0 + 731);

    const auto *isi1_665 = buffer.data(isi1 + 665);
    const auto *isi1_667 = buffer.data(isi1 + 667);
    const auto *isi1_668 = buffer.data(isi1 + 668);
    const auto *isi1_669 = buffer.data(isi1 + 669);
    const auto *isi1_670 = buffer.data(isi1 + 670);
    const auto *isi1_671 = buffer.data(isi1 + 671);
    const auto *isi1_672 = buffer.data(isi1 + 672);
    const auto *isi1_673 = buffer.data(isi1 + 673);
    const auto *isi1_674 = buffer.data(isi1 + 674);
    const auto *isi1_675 = buffer.data(isi1 + 675);
    const auto *isi1_676 = buffer.data(isi1 + 676);
    const auto *isi1_677 = buffer.data(isi1 + 677);
    const auto *isi1_678 = buffer.data(isi1 + 678);
    const auto *isi1_679 = buffer.data(isi1 + 679);
    const auto *isi1_680 = buffer.data(isi1 + 680);
    const auto *isi1_681 = buffer.data(isi1 + 681);
    const auto *isi1_682 = buffer.data(isi1 + 682);
    const auto *isi1_683 = buffer.data(isi1 + 683);
    const auto *isi1_684 = buffer.data(isi1 + 684);
    const auto *isi1_685 = buffer.data(isi1 + 685);
    const auto *isi1_686 = buffer.data(isi1 + 686);
    const auto *isi1_687 = buffer.data(isi1 + 687);
    const auto *isi1_688 = buffer.data(isi1 + 688);
    const auto *isi1_689 = buffer.data(isi1 + 689);
    const auto *isi1_690 = buffer.data(isi1 + 690);
    const auto *isi1_691 = buffer.data(isi1 + 691);
    const auto *isi1_692 = buffer.data(isi1 + 692);
    const auto *isi1_693 = buffer.data(isi1 + 693);
    const auto *isi1_694 = buffer.data(isi1 + 694);
    const auto *isi1_695 = buffer.data(isi1 + 695);
    const auto *isi1_696 = buffer.data(isi1 + 696);
    const auto *isi1_697 = buffer.data(isi1 + 697);
    const auto *isi1_698 = buffer.data(isi1 + 698);
    const auto *isi1_699 = buffer.data(isi1 + 699);
    const auto *isi1_700 = buffer.data(isi1 + 700);
    const auto *isi1_701 = buffer.data(isi1 + 701);
    const auto *isi1_702 = buffer.data(isi1 + 702);
    const auto *isi1_703 = buffer.data(isi1 + 703);
    const auto *isi1_704 = buffer.data(isi1 + 704);
    const auto *isi1_705 = buffer.data(isi1 + 705);
    const auto *isi1_706 = buffer.data(isi1 + 706);
    const auto *isi1_707 = buffer.data(isi1 + 707);
    const auto *isi1_708 = buffer.data(isi1 + 708);
    const auto *isi1_709 = buffer.data(isi1 + 709);
    const auto *isi1_710 = buffer.data(isi1 + 710);
    const auto *isi1_711 = buffer.data(isi1 + 711);
    const auto *isi1_712 = buffer.data(isi1 + 712);
    const auto *isi1_713 = buffer.data(isi1 + 713);
    const auto *isi1_714 = buffer.data(isi1 + 714);
    const auto *isi1_715 = buffer.data(isi1 + 715);
    const auto *isi1_716 = buffer.data(isi1 + 716);
    const auto *isi1_717 = buffer.data(isi1 + 717);
    const auto *isi1_718 = buffer.data(isi1 + 718);
    const auto *isi1_719 = buffer.data(isi1 + 719);
    const auto *isi1_720 = buffer.data(isi1 + 720);
    const auto *isi1_721 = buffer.data(isi1 + 721);
    const auto *isi1_722 = buffer.data(isi1 + 722);
    const auto *isi1_723 = buffer.data(isi1 + 723);
    const auto *isi1_724 = buffer.data(isi1 + 724);
    const auto *isi1_725 = buffer.data(isi1 + 725);
    const auto *isi1_726 = buffer.data(isi1 + 726);
    const auto *isi1_727 = buffer.data(isi1 + 727);
    const auto *isi1_729 = buffer.data(isi1 + 729);
    const auto *isi1_731 = buffer.data(isi1 + 731);

    const auto *isk_853 = buffer.data(isk + 853);
    const auto *isk_854 = buffer.data(isk + 854);
    const auto *isk_855 = buffer.data(isk + 855);
    const auto *isk_856 = buffer.data(isk + 856);
    const auto *isk_857 = buffer.data(isk + 857);
    const auto *isk_858 = buffer.data(isk + 858);
    const auto *isk_859 = buffer.data(isk + 859);
    const auto *isk_860 = buffer.data(isk + 860);
    const auto *isk_861 = buffer.data(isk + 861);
    const auto *isk_862 = buffer.data(isk + 862);
    const auto *isk_863 = buffer.data(isk + 863);
    const auto *isk_864 = buffer.data(isk + 864);
    const auto *isk_865 = buffer.data(isk + 865);
    const auto *isk_866 = buffer.data(isk + 866);
    const auto *isk_867 = buffer.data(isk + 867);
    const auto *isk_868 = buffer.data(isk + 868);
    const auto *isk_869 = buffer.data(isk + 869);
    const auto *isk_870 = buffer.data(isk + 870);
    const auto *isk_871 = buffer.data(isk + 871);
    const auto *isk_872 = buffer.data(isk + 872);
    const auto *isk_873 = buffer.data(isk + 873);
    const auto *isk_874 = buffer.data(isk + 874);
    const auto *isk_875 = buffer.data(isk + 875);
    const auto *isk_876 = buffer.data(isk + 876);
    const auto *isk_877 = buffer.data(isk + 877);
    const auto *isk_878 = buffer.data(isk + 878);
    const auto *isk_879 = buffer.data(isk + 879);
    const auto *isk_880 = buffer.data(isk + 880);
    const auto *isk_881 = buffer.data(isk + 881);
    const auto *isk_882 = buffer.data(isk + 882);
    const auto *isk_883 = buffer.data(isk + 883);
    const auto *isk_884 = buffer.data(isk + 884);
    const auto *isk_885 = buffer.data(isk + 885);
    const auto *isk_886 = buffer.data(isk + 886);
    const auto *isk_887 = buffer.data(isk + 887);
    const auto *isk_888 = buffer.data(isk + 888);
    const auto *isk_889 = buffer.data(isk + 889);
    const auto *isk_890 = buffer.data(isk + 890);
    const auto *isk_891 = buffer.data(isk + 891);
    const auto *isk_892 = buffer.data(isk + 892);
    const auto *isk_893 = buffer.data(isk + 893);
    const auto *isk_894 = buffer.data(isk + 894);
    const auto *isk_895 = buffer.data(isk + 895);
    const auto *isk_896 = buffer.data(isk + 896);
    const auto *isk_897 = buffer.data(isk + 897);
    const auto *isk_898 = buffer.data(isk + 898);
    const auto *isk_899 = buffer.data(isk + 899);
    const auto *isk_900 = buffer.data(isk + 900);
    const auto *isk_901 = buffer.data(isk + 901);
    const auto *isk_902 = buffer.data(isk + 902);
    const auto *isk_903 = buffer.data(isk + 903);
    const auto *isk_904 = buffer.data(isk + 904);
    const auto *isk_905 = buffer.data(isk + 905);
    const auto *isk_906 = buffer.data(isk + 906);
    const auto *isk_907 = buffer.data(isk + 907);
    const auto *isk_908 = buffer.data(isk + 908);
    const auto *isk_909 = buffer.data(isk + 909);
    const auto *isk_910 = buffer.data(isk + 910);
    const auto *isk_911 = buffer.data(isk + 911);
    const auto *isk_912 = buffer.data(isk + 912);
    const auto *isk_913 = buffer.data(isk + 913);
    const auto *isk_914 = buffer.data(isk + 914);
    const auto *isk_915 = buffer.data(isk + 915);
    const auto *isk_916 = buffer.data(isk + 916);
    const auto *isk_917 = buffer.data(isk + 917);
    const auto *isk_918 = buffer.data(isk + 918);
    const auto *isk_919 = buffer.data(isk + 919);
    const auto *isk_920 = buffer.data(isk + 920);
    const auto *isk_921 = buffer.data(isk + 921);
    const auto *isk_922 = buffer.data(isk + 922);
    const auto *isk_923 = buffer.data(isk + 923);
    const auto *isk_924 = buffer.data(isk + 924);
    const auto *isk_925 = buffer.data(isk + 925);
    const auto *isk_926 = buffer.data(isk + 926);
    const auto *isk_927 = buffer.data(isk + 927);
    const auto *isk_928 = buffer.data(isk + 928);
    const auto *isk_929 = buffer.data(isk + 929);
    const auto *isk_930 = buffer.data(isk + 930);
    const auto *isk_931 = buffer.data(isk + 931);
    const auto *isk_932 = buffer.data(isk + 932);
    const auto *isk_933 = buffer.data(isk + 933);
    const auto *isk_934 = buffer.data(isk + 934);
    const auto *isk_935 = buffer.data(isk + 935);
    const auto *isk_937 = buffer.data(isk + 937);
    const auto *isk_939 = buffer.data(isk + 939);

#pragma omp simd aligned(t_1060, t_1061, t_1062, t_1063, pc_x, isi0_669, isi0_670, isi0_671, \
                         isi1_669, isi1_670, isi1_671, isk_853, isk_854, isk_855, \
                         isk_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1060[k] = f_4 * isi0_669[k]
                    - f_5 * isi1_669[k]
                    + f_3 * pc_x[k] * isk_853[k];

        t_1061[k] = f_4 * isi0_670[k]
                    - f_5 * isi1_670[k]
                    + f_3 * pc_x[k] * isk_854[k];

        t_1062[k] = f_4 * isi0_671[k]
                    - f_5 * isi1_671[k]
                    + f_3 * pc_x[k] * isk_855[k];

        t_1063[k] = f_3 * pc_x[k] * isk_856[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, t_1067, t_1068, t_1069, t_1070, pc_x, \
                         isk_857, isk_858, isk_859, isk_860, isk_861, isk_862, \
                         isk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_3 * pc_x[k] * isk_857[k];

        t_1065[k] = f_3 * pc_x[k] * isk_858[k];

        t_1066[k] = f_3 * pc_x[k] * isk_859[k];

        t_1067[k] = f_3 * pc_x[k] * isk_860[k];

        t_1068[k] = f_3 * pc_x[k] * isk_861[k];

        t_1069[k] = f_3 * pc_x[k] * isk_862[k];

        t_1070[k] = f_3 * pc_x[k] * isk_863[k];
    }

#pragma omp simd aligned(t_1071, t_1072, t_1073, pc_y, pc_z, hsk_604, hsk_640, hsk_642, \
                         isi0_665, isi0_667, isi1_665, isi1_667, isk_856, \
                         isk_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1071[k] = f_18 * hsk_640[k]
                    + f_1 * isi0_665[k]
                    - f_2 * isi1_665[k]
                    + f_3 * pc_y[k] * isk_856[k];

        t_1072[k] = f_16 * hsk_604[k]
                    + f_3 * pc_z[k] * isk_856[k];

        t_1073[k] = f_18 * hsk_642[k]
                    + f_12 * isi0_667[k]
                    - f_13 * isi1_667[k]
                    + f_3 * pc_y[k] * isk_858[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, pc_y, hsk_643, hsk_644, hsk_645, isi0_668, \
                         isi0_669, isi0_670, isi1_668, isi1_669, isi1_670, isk_859, isk_860, \
                         isk_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = f_18 * hsk_643[k]
                    + f_10 * isi0_668[k]
                    - f_11 * isi1_668[k]
                    + f_3 * pc_y[k] * isk_859[k];

        t_1075[k] = f_18 * hsk_644[k]
                    + f_8 * isi0_669[k]
                    - f_9 * isi1_669[k]
                    + f_3 * pc_y[k] * isk_860[k];

        t_1076[k] = f_18 * hsk_645[k]
                    + f_6 * isi0_670[k]
                    - f_7 * isi1_670[k]
                    + f_3 * pc_y[k] * isk_861[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, pc_y, pc_z, hsk_611, hsk_646, hsk_647, \
                         isi0_671, isi1_671, isk_862, isk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = f_18 * hsk_646[k]
                    + f_4 * isi0_671[k]
                    - f_5 * isi1_671[k]
                    + f_3 * pc_y[k] * isk_862[k];

        t_1078[k] = f_18 * hsk_647[k]
                    + f_3 * pc_y[k] * isk_863[k];

        t_1079[k] = f_16 * hsk_611[k]
                    + f_1 * isi0_671[k]
                    - f_2 * isi1_671[k]
                    + f_3 * pc_z[k] * isk_863[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, pc_x, isi0_672, isi0_673, isi0_674, isi1_672, \
                         isi1_673, isi1_674, isk_864, isk_865, \
                         isk_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_1 * isi0_672[k]
                    - f_2 * isi1_672[k]
                    + f_3 * pc_x[k] * isk_864[k];

        t_1081[k] = f_20 * isi0_673[k]
                    - f_21 * isi1_673[k]
                    + f_3 * pc_x[k] * isk_865[k];

        t_1082[k] = f_20 * isi0_674[k]
                    - f_21 * isi1_674[k]
                    + f_3 * pc_x[k] * isk_866[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, pc_x, isi0_675, isi0_676, isi0_677, isi1_675, \
                         isi1_676, isi1_677, isk_867, isk_868, \
                         isk_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_12 * isi0_675[k]
                    - f_13 * isi1_675[k]
                    + f_3 * pc_x[k] * isk_867[k];

        t_1084[k] = f_12 * isi0_676[k]
                    - f_13 * isi1_676[k]
                    + f_3 * pc_x[k] * isk_868[k];

        t_1085[k] = f_12 * isi0_677[k]
                    - f_13 * isi1_677[k]
                    + f_3 * pc_x[k] * isk_869[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pc_x, isi0_678, isi0_679, isi0_680, isi1_678, \
                         isi1_679, isi1_680, isk_870, isk_871, \
                         isk_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_10 * isi0_678[k]
                    - f_11 * isi1_678[k]
                    + f_3 * pc_x[k] * isk_870[k];

        t_1087[k] = f_10 * isi0_679[k]
                    - f_11 * isi1_679[k]
                    + f_3 * pc_x[k] * isk_871[k];

        t_1088[k] = f_10 * isi0_680[k]
                    - f_11 * isi1_680[k]
                    + f_3 * pc_x[k] * isk_872[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pc_x, isi0_681, isi0_682, isi0_683, isi1_681, \
                         isi1_682, isi1_683, isk_873, isk_874, \
                         isk_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_10 * isi0_681[k]
                    - f_11 * isi1_681[k]
                    + f_3 * pc_x[k] * isk_873[k];

        t_1090[k] = f_8 * isi0_682[k]
                    - f_9 * isi1_682[k]
                    + f_3 * pc_x[k] * isk_874[k];

        t_1091[k] = f_8 * isi0_683[k]
                    - f_9 * isi1_683[k]
                    + f_3 * pc_x[k] * isk_875[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pc_x, isi0_684, isi0_685, isi0_686, isi1_684, \
                         isi1_685, isi1_686, isk_876, isk_877, \
                         isk_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_8 * isi0_684[k]
                    - f_9 * isi1_684[k]
                    + f_3 * pc_x[k] * isk_876[k];

        t_1093[k] = f_8 * isi0_685[k]
                    - f_9 * isi1_685[k]
                    + f_3 * pc_x[k] * isk_877[k];

        t_1094[k] = f_8 * isi0_686[k]
                    - f_9 * isi1_686[k]
                    + f_3 * pc_x[k] * isk_878[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pc_x, isi0_687, isi0_688, isi0_689, isi1_687, \
                         isi1_688, isi1_689, isk_879, isk_880, \
                         isk_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = f_6 * isi0_687[k]
                    - f_7 * isi1_687[k]
                    + f_3 * pc_x[k] * isk_879[k];

        t_1096[k] = f_6 * isi0_688[k]
                    - f_7 * isi1_688[k]
                    + f_3 * pc_x[k] * isk_880[k];

        t_1097[k] = f_6 * isi0_689[k]
                    - f_7 * isi1_689[k]
                    + f_3 * pc_x[k] * isk_881[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pc_x, isi0_690, isi0_691, isi0_692, isi1_690, \
                         isi1_691, isi1_692, isk_882, isk_883, \
                         isk_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = f_6 * isi0_690[k]
                    - f_7 * isi1_690[k]
                    + f_3 * pc_x[k] * isk_882[k];

        t_1099[k] = f_6 * isi0_691[k]
                    - f_7 * isi1_691[k]
                    + f_3 * pc_x[k] * isk_883[k];

        t_1100[k] = f_6 * isi0_692[k]
                    - f_7 * isi1_692[k]
                    + f_3 * pc_x[k] * isk_884[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, pc_x, isi0_693, isi0_694, isi0_695, isi1_693, \
                         isi1_694, isi1_695, isk_885, isk_886, \
                         isk_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = f_4 * isi0_693[k]
                    - f_5 * isi1_693[k]
                    + f_3 * pc_x[k] * isk_885[k];

        t_1102[k] = f_4 * isi0_694[k]
                    - f_5 * isi1_694[k]
                    + f_3 * pc_x[k] * isk_886[k];

        t_1103[k] = f_4 * isi0_695[k]
                    - f_5 * isi1_695[k]
                    + f_3 * pc_x[k] * isk_887[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, pc_x, isi0_696, isi0_697, isi0_698, isi1_696, \
                         isi1_697, isi1_698, isk_888, isk_889, \
                         isk_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_4 * isi0_696[k]
                    - f_5 * isi1_696[k]
                    + f_3 * pc_x[k] * isk_888[k];

        t_1105[k] = f_4 * isi0_697[k]
                    - f_5 * isi1_697[k]
                    + f_3 * pc_x[k] * isk_889[k];

        t_1106[k] = f_4 * isi0_698[k]
                    - f_5 * isi1_698[k]
                    + f_3 * pc_x[k] * isk_890[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, t_1110, t_1111, t_1112, pc_x, isi0_699, \
                         isi1_699, isk_891, isk_892, isk_893, isk_894, isk_895, \
                         isk_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = f_4 * isi0_699[k]
                    - f_5 * isi1_699[k]
                    + f_3 * pc_x[k] * isk_891[k];

        t_1108[k] = f_3 * pc_x[k] * isk_892[k];

        t_1109[k] = f_3 * pc_x[k] * isk_893[k];

        t_1110[k] = f_3 * pc_x[k] * isk_894[k];

        t_1111[k] = f_3 * pc_x[k] * isk_895[k];

        t_1112[k] = f_3 * pc_x[k] * isk_896[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, t_1116, t_1117, pc_x, pc_y, pc_z, hsk_640, \
                         hsk_676, isi0_693, isi1_693, isk_892, isk_897, isk_898, \
                         isk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = f_3 * pc_x[k] * isk_897[k];

        t_1114[k] = f_3 * pc_x[k] * isk_898[k];

        t_1115[k] = f_3 * pc_x[k] * isk_899[k];

        t_1116[k] = f_17 * hsk_676[k]
                    + f_1 * isi0_693[k]
                    - f_2 * isi1_693[k]
                    + f_3 * pc_y[k] * isk_892[k];

        t_1117[k] = f_17 * hsk_640[k]
                    + f_3 * pc_z[k] * isk_892[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pc_y, hsk_678, hsk_679, hsk_680, isi0_695, \
                         isi0_696, isi0_697, isi1_695, isi1_696, isi1_697, isk_894, isk_895, \
                         isk_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_17 * hsk_678[k]
                    + f_12 * isi0_695[k]
                    - f_13 * isi1_695[k]
                    + f_3 * pc_y[k] * isk_894[k];

        t_1119[k] = f_17 * hsk_679[k]
                    + f_10 * isi0_696[k]
                    - f_11 * isi1_696[k]
                    + f_3 * pc_y[k] * isk_895[k];

        t_1120[k] = f_17 * hsk_680[k]
                    + f_8 * isi0_697[k]
                    - f_9 * isi1_697[k]
                    + f_3 * pc_y[k] * isk_896[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_y, hsk_681, hsk_682, hsk_683, isi0_698, \
                         isi0_699, isi1_698, isi1_699, isk_897, isk_898, \
                         isk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_17 * hsk_681[k]
                    + f_6 * isi0_698[k]
                    - f_7 * isi1_698[k]
                    + f_3 * pc_y[k] * isk_897[k];

        t_1122[k] = f_17 * hsk_682[k]
                    + f_4 * isi0_699[k]
                    - f_5 * isi1_699[k]
                    + f_3 * pc_y[k] * isk_898[k];

        t_1123[k] = f_17 * hsk_683[k]
                    + f_3 * pc_y[k] * isk_899[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_z, hsk_647, isi0_699, isi0_700, \
                         isi0_701, isi1_699, isi1_700, isi1_701, isk_899, isk_900, \
                         isk_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_17 * hsk_647[k]
                    + f_1 * isi0_699[k]
                    - f_2 * isi1_699[k]
                    + f_3 * pc_z[k] * isk_899[k];

        t_1125[k] = f_1 * isi0_700[k]
                    - f_2 * isi1_700[k]
                    + f_3 * pc_x[k] * isk_900[k];

        t_1126[k] = f_20 * isi0_701[k]
                    - f_21 * isi1_701[k]
                    + f_3 * pc_x[k] * isk_901[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, isi0_702, isi0_703, isi0_704, isi1_702, \
                         isi1_703, isi1_704, isk_902, isk_903, \
                         isk_904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_20 * isi0_702[k]
                    - f_21 * isi1_702[k]
                    + f_3 * pc_x[k] * isk_902[k];

        t_1128[k] = f_12 * isi0_703[k]
                    - f_13 * isi1_703[k]
                    + f_3 * pc_x[k] * isk_903[k];

        t_1129[k] = f_12 * isi0_704[k]
                    - f_13 * isi1_704[k]
                    + f_3 * pc_x[k] * isk_904[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, isi0_705, isi0_706, isi0_707, isi1_705, \
                         isi1_706, isi1_707, isk_905, isk_906, \
                         isk_907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_12 * isi0_705[k]
                    - f_13 * isi1_705[k]
                    + f_3 * pc_x[k] * isk_905[k];

        t_1131[k] = f_10 * isi0_706[k]
                    - f_11 * isi1_706[k]
                    + f_3 * pc_x[k] * isk_906[k];

        t_1132[k] = f_10 * isi0_707[k]
                    - f_11 * isi1_707[k]
                    + f_3 * pc_x[k] * isk_907[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, isi0_708, isi0_709, isi0_710, isi1_708, \
                         isi1_709, isi1_710, isk_908, isk_909, \
                         isk_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_10 * isi0_708[k]
                    - f_11 * isi1_708[k]
                    + f_3 * pc_x[k] * isk_908[k];

        t_1134[k] = f_10 * isi0_709[k]
                    - f_11 * isi1_709[k]
                    + f_3 * pc_x[k] * isk_909[k];

        t_1135[k] = f_8 * isi0_710[k]
                    - f_9 * isi1_710[k]
                    + f_3 * pc_x[k] * isk_910[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, pc_x, isi0_711, isi0_712, isi0_713, isi1_711, \
                         isi1_712, isi1_713, isk_911, isk_912, \
                         isk_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_8 * isi0_711[k]
                    - f_9 * isi1_711[k]
                    + f_3 * pc_x[k] * isk_911[k];

        t_1137[k] = f_8 * isi0_712[k]
                    - f_9 * isi1_712[k]
                    + f_3 * pc_x[k] * isk_912[k];

        t_1138[k] = f_8 * isi0_713[k]
                    - f_9 * isi1_713[k]
                    + f_3 * pc_x[k] * isk_913[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, pc_x, isi0_714, isi0_715, isi0_716, isi1_714, \
                         isi1_715, isi1_716, isk_914, isk_915, \
                         isk_916 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_8 * isi0_714[k]
                    - f_9 * isi1_714[k]
                    + f_3 * pc_x[k] * isk_914[k];

        t_1140[k] = f_6 * isi0_715[k]
                    - f_7 * isi1_715[k]
                    + f_3 * pc_x[k] * isk_915[k];

        t_1141[k] = f_6 * isi0_716[k]
                    - f_7 * isi1_716[k]
                    + f_3 * pc_x[k] * isk_916[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, pc_x, isi0_717, isi0_718, isi0_719, isi1_717, \
                         isi1_718, isi1_719, isk_917, isk_918, \
                         isk_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_6 * isi0_717[k]
                    - f_7 * isi1_717[k]
                    + f_3 * pc_x[k] * isk_917[k];

        t_1143[k] = f_6 * isi0_718[k]
                    - f_7 * isi1_718[k]
                    + f_3 * pc_x[k] * isk_918[k];

        t_1144[k] = f_6 * isi0_719[k]
                    - f_7 * isi1_719[k]
                    + f_3 * pc_x[k] * isk_919[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, pc_x, isi0_720, isi0_721, isi0_722, isi1_720, \
                         isi1_721, isi1_722, isk_920, isk_921, \
                         isk_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_6 * isi0_720[k]
                    - f_7 * isi1_720[k]
                    + f_3 * pc_x[k] * isk_920[k];

        t_1146[k] = f_4 * isi0_721[k]
                    - f_5 * isi1_721[k]
                    + f_3 * pc_x[k] * isk_921[k];

        t_1147[k] = f_4 * isi0_722[k]
                    - f_5 * isi1_722[k]
                    + f_3 * pc_x[k] * isk_922[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, pc_x, isi0_723, isi0_724, isi0_725, isi1_723, \
                         isi1_724, isi1_725, isk_923, isk_924, \
                         isk_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = f_4 * isi0_723[k]
                    - f_5 * isi1_723[k]
                    + f_3 * pc_x[k] * isk_923[k];

        t_1149[k] = f_4 * isi0_724[k]
                    - f_5 * isi1_724[k]
                    + f_3 * pc_x[k] * isk_924[k];

        t_1150[k] = f_4 * isi0_725[k]
                    - f_5 * isi1_725[k]
                    + f_3 * pc_x[k] * isk_925[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, t_1155, pc_x, isi0_726, isi0_727, \
                         isi1_726, isi1_727, isk_926, isk_927, isk_928, isk_929, \
                         isk_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_4 * isi0_726[k]
                    - f_5 * isi1_726[k]
                    + f_3 * pc_x[k] * isk_926[k];

        t_1152[k] = f_4 * isi0_727[k]
                    - f_5 * isi1_727[k]
                    + f_3 * pc_x[k] * isk_927[k];

        t_1153[k] = f_3 * pc_x[k] * isk_928[k];

        t_1154[k] = f_3 * pc_x[k] * isk_929[k];

        t_1155[k] = f_3 * pc_x[k] * isk_930[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, t_1159, t_1160, pc_x, isk_931, isk_932, \
                         isk_933, isk_934, isk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_3 * pc_x[k] * isk_931[k];

        t_1157[k] = f_3 * pc_x[k] * isk_932[k];

        t_1158[k] = f_3 * pc_x[k] * isk_933[k];

        t_1159[k] = f_3 * pc_x[k] * isk_934[k];

        t_1160[k] = f_3 * pc_x[k] * isk_935[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, pc_y, pc_z, hsk_676, hsk_712, hsk_714, \
                         isi0_721, isi0_723, isi1_721, isi1_723, isk_928, \
                         isk_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_16 * hsk_712[k]
                    + f_1 * isi0_721[k]
                    - f_2 * isi1_721[k]
                    + f_3 * pc_y[k] * isk_928[k];

        t_1162[k] = f_18 * hsk_676[k]
                    + f_3 * pc_z[k] * isk_928[k];

        t_1163[k] = f_16 * hsk_714[k]
                    + f_12 * isi0_723[k]
                    - f_13 * isi1_723[k]
                    + f_3 * pc_y[k] * isk_930[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, pc_y, hsk_715, hsk_716, hsk_717, isi0_724, \
                         isi0_725, isi0_726, isi1_724, isi1_725, isi1_726, isk_931, isk_932, \
                         isk_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = f_16 * hsk_715[k]
                    + f_10 * isi0_724[k]
                    - f_11 * isi1_724[k]
                    + f_3 * pc_y[k] * isk_931[k];

        t_1165[k] = f_16 * hsk_716[k]
                    + f_8 * isi0_725[k]
                    - f_9 * isi1_725[k]
                    + f_3 * pc_y[k] * isk_932[k];

        t_1166[k] = f_16 * hsk_717[k]
                    + f_6 * isi0_726[k]
                    - f_7 * isi1_726[k]
                    + f_3 * pc_y[k] * isk_933[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, t_1170, pa_y, pc_y, pc_z, hsl0_900, hsk_683, \
                         hsk_718, hsk_719, hsl1_900, isi0_727, isi1_727, isk_934, \
                         isk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = f_16 * hsk_718[k]
                    + f_4 * isi0_727[k]
                    - f_5 * isi1_727[k]
                    + f_3 * pc_y[k] * isk_934[k];

        t_1168[k] = f_16 * hsk_719[k]
                    + f_3 * pc_y[k] * isk_935[k];

        t_1169[k] = f_18 * hsk_683[k]
                    + f_1 * isi0_727[k]
                    - f_2 * isi1_727[k]
                    + f_3 * pc_z[k] * isk_935[k];

        t_1170[k] = pa_y[k] * hsl0_900[k]
                    - f_14 * pc_y[k] * hsl1_900[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, pa_y, pc_x, pc_y, hsl0_902, hsl1_902, \
                         isi0_729, isi0_731, isi1_729, isi1_731, isk_937, \
                         isk_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_20 * isi0_729[k]
                    - f_21 * isi1_729[k]
                    + f_3 * pc_x[k] * isk_937[k];

        t_1172[k] = pa_y[k] * hsl0_902[k]
                    - f_14 * pc_y[k] * hsl1_902[k];

        t_1173[k] = f_12 * isi0_731[k]
                    - f_13 * isi1_731[k]
                    + f_3 * pc_x[k] * isk_939[k];
    }
}

static auto
compute_prim_isl_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t hsl0,
                                                           const size_t hsk, const size_t hsl1,
                                                           const size_t isi0, const size_t isi1,
                                                           const size_t isk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_20 = 3.0 / gamma;
    const auto f_21 = 3.0 * p / (gamma * q);
    const auto f_22 = 4.0 / q;

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

    const auto *hsl0_905 = buffer.data(hsl0 + 905);
    const auto *hsl0_909 = buffer.data(hsl0 + 909);
    const auto *hsl0_914 = buffer.data(hsl0 + 914);
    const auto *hsl0_920 = buffer.data(hsl0 + 920);
    const auto *hsl0_927 = buffer.data(hsl0 + 927);
    const auto *hsl0_936 = buffer.data(hsl0 + 936);
    const auto *hsl0_938 = buffer.data(hsl0 + 938);
    const auto *hsl0_939 = buffer.data(hsl0 + 939);
    const auto *hsl0_940 = buffer.data(hsl0 + 940);
    const auto *hsl0_941 = buffer.data(hsl0 + 941);
    const auto *hsl0_942 = buffer.data(hsl0 + 942);
    const auto *hsl0_944 = buffer.data(hsl0 + 944);

    const auto *hsk_712 = buffer.data(hsk + 712);
    const auto *hsk_748 = buffer.data(hsk + 748);
    const auto *hsk_750 = buffer.data(hsk + 750);
    const auto *hsk_751 = buffer.data(hsk + 751);
    const auto *hsk_752 = buffer.data(hsk + 752);
    const auto *hsk_753 = buffer.data(hsk + 753);
    const auto *hsk_754 = buffer.data(hsk + 754);
    const auto *hsk_755 = buffer.data(hsk + 755);

    const auto *hsl1_905 = buffer.data(hsl1 + 905);
    const auto *hsl1_909 = buffer.data(hsl1 + 909);
    const auto *hsl1_914 = buffer.data(hsl1 + 914);
    const auto *hsl1_920 = buffer.data(hsl1 + 920);
    const auto *hsl1_927 = buffer.data(hsl1 + 927);
    const auto *hsl1_936 = buffer.data(hsl1 + 936);
    const auto *hsl1_938 = buffer.data(hsl1 + 938);
    const auto *hsl1_939 = buffer.data(hsl1 + 939);
    const auto *hsl1_940 = buffer.data(hsl1 + 940);
    const auto *hsl1_941 = buffer.data(hsl1 + 941);
    const auto *hsl1_942 = buffer.data(hsl1 + 942);
    const auto *hsl1_944 = buffer.data(hsl1 + 944);

    const auto *isi0_732 = buffer.data(isi0 + 732);
    const auto *isi0_734 = buffer.data(isi0 + 734);
    const auto *isi0_735 = buffer.data(isi0 + 735);
    const auto *isi0_736 = buffer.data(isi0 + 736);
    const auto *isi0_738 = buffer.data(isi0 + 738);
    const auto *isi0_739 = buffer.data(isi0 + 739);
    const auto *isi0_740 = buffer.data(isi0 + 740);
    const auto *isi0_741 = buffer.data(isi0 + 741);
    const auto *isi0_743 = buffer.data(isi0 + 743);
    const auto *isi0_744 = buffer.data(isi0 + 744);
    const auto *isi0_745 = buffer.data(isi0 + 745);
    const auto *isi0_746 = buffer.data(isi0 + 746);
    const auto *isi0_747 = buffer.data(isi0 + 747);
    const auto *isi0_749 = buffer.data(isi0 + 749);
    const auto *isi0_750 = buffer.data(isi0 + 750);
    const auto *isi0_751 = buffer.data(isi0 + 751);
    const auto *isi0_752 = buffer.data(isi0 + 752);
    const auto *isi0_753 = buffer.data(isi0 + 753);
    const auto *isi0_754 = buffer.data(isi0 + 754);
    const auto *isi0_756 = buffer.data(isi0 + 756);
    const auto *isi0_758 = buffer.data(isi0 + 758);
    const auto *isi0_759 = buffer.data(isi0 + 759);
    const auto *isi0_761 = buffer.data(isi0 + 761);
    const auto *isi0_762 = buffer.data(isi0 + 762);
    const auto *isi0_763 = buffer.data(isi0 + 763);
    const auto *isi0_765 = buffer.data(isi0 + 765);
    const auto *isi0_766 = buffer.data(isi0 + 766);
    const auto *isi0_767 = buffer.data(isi0 + 767);
    const auto *isi0_768 = buffer.data(isi0 + 768);
    const auto *isi0_770 = buffer.data(isi0 + 770);
    const auto *isi0_771 = buffer.data(isi0 + 771);
    const auto *isi0_772 = buffer.data(isi0 + 772);
    const auto *isi0_773 = buffer.data(isi0 + 773);
    const auto *isi0_774 = buffer.data(isi0 + 774);
    const auto *isi0_776 = buffer.data(isi0 + 776);
    const auto *isi0_777 = buffer.data(isi0 + 777);
    const auto *isi0_778 = buffer.data(isi0 + 778);
    const auto *isi0_779 = buffer.data(isi0 + 779);
    const auto *isi0_780 = buffer.data(isi0 + 780);
    const auto *isi0_781 = buffer.data(isi0 + 781);
    const auto *isi0_782 = buffer.data(isi0 + 782);
    const auto *isi0_783 = buffer.data(isi0 + 783);

    const auto *isi1_732 = buffer.data(isi1 + 732);
    const auto *isi1_734 = buffer.data(isi1 + 734);
    const auto *isi1_735 = buffer.data(isi1 + 735);
    const auto *isi1_736 = buffer.data(isi1 + 736);
    const auto *isi1_738 = buffer.data(isi1 + 738);
    const auto *isi1_739 = buffer.data(isi1 + 739);
    const auto *isi1_740 = buffer.data(isi1 + 740);
    const auto *isi1_741 = buffer.data(isi1 + 741);
    const auto *isi1_743 = buffer.data(isi1 + 743);
    const auto *isi1_744 = buffer.data(isi1 + 744);
    const auto *isi1_745 = buffer.data(isi1 + 745);
    const auto *isi1_746 = buffer.data(isi1 + 746);
    const auto *isi1_747 = buffer.data(isi1 + 747);
    const auto *isi1_749 = buffer.data(isi1 + 749);
    const auto *isi1_750 = buffer.data(isi1 + 750);
    const auto *isi1_751 = buffer.data(isi1 + 751);
    const auto *isi1_752 = buffer.data(isi1 + 752);
    const auto *isi1_753 = buffer.data(isi1 + 753);
    const auto *isi1_754 = buffer.data(isi1 + 754);
    const auto *isi1_756 = buffer.data(isi1 + 756);
    const auto *isi1_758 = buffer.data(isi1 + 758);
    const auto *isi1_759 = buffer.data(isi1 + 759);
    const auto *isi1_761 = buffer.data(isi1 + 761);
    const auto *isi1_762 = buffer.data(isi1 + 762);
    const auto *isi1_763 = buffer.data(isi1 + 763);
    const auto *isi1_765 = buffer.data(isi1 + 765);
    const auto *isi1_766 = buffer.data(isi1 + 766);
    const auto *isi1_767 = buffer.data(isi1 + 767);
    const auto *isi1_768 = buffer.data(isi1 + 768);
    const auto *isi1_770 = buffer.data(isi1 + 770);
    const auto *isi1_771 = buffer.data(isi1 + 771);
    const auto *isi1_772 = buffer.data(isi1 + 772);
    const auto *isi1_773 = buffer.data(isi1 + 773);
    const auto *isi1_774 = buffer.data(isi1 + 774);
    const auto *isi1_776 = buffer.data(isi1 + 776);
    const auto *isi1_777 = buffer.data(isi1 + 777);
    const auto *isi1_778 = buffer.data(isi1 + 778);
    const auto *isi1_779 = buffer.data(isi1 + 779);
    const auto *isi1_780 = buffer.data(isi1 + 780);
    const auto *isi1_781 = buffer.data(isi1 + 781);
    const auto *isi1_782 = buffer.data(isi1 + 782);
    const auto *isi1_783 = buffer.data(isi1 + 783);

    const auto *isk_940 = buffer.data(isk + 940);
    const auto *isk_942 = buffer.data(isk + 942);
    const auto *isk_943 = buffer.data(isk + 943);
    const auto *isk_944 = buffer.data(isk + 944);
    const auto *isk_946 = buffer.data(isk + 946);
    const auto *isk_947 = buffer.data(isk + 947);
    const auto *isk_948 = buffer.data(isk + 948);
    const auto *isk_949 = buffer.data(isk + 949);
    const auto *isk_951 = buffer.data(isk + 951);
    const auto *isk_952 = buffer.data(isk + 952);
    const auto *isk_953 = buffer.data(isk + 953);
    const auto *isk_954 = buffer.data(isk + 954);
    const auto *isk_955 = buffer.data(isk + 955);
    const auto *isk_957 = buffer.data(isk + 957);
    const auto *isk_958 = buffer.data(isk + 958);
    const auto *isk_959 = buffer.data(isk + 959);
    const auto *isk_960 = buffer.data(isk + 960);
    const auto *isk_961 = buffer.data(isk + 961);
    const auto *isk_962 = buffer.data(isk + 962);
    const auto *isk_964 = buffer.data(isk + 964);
    const auto *isk_965 = buffer.data(isk + 965);
    const auto *isk_966 = buffer.data(isk + 966);
    const auto *isk_967 = buffer.data(isk + 967);
    const auto *isk_968 = buffer.data(isk + 968);
    const auto *isk_969 = buffer.data(isk + 969);
    const auto *isk_970 = buffer.data(isk + 970);
    const auto *isk_971 = buffer.data(isk + 971);
    const auto *isk_972 = buffer.data(isk + 972);
    const auto *isk_974 = buffer.data(isk + 974);
    const auto *isk_975 = buffer.data(isk + 975);
    const auto *isk_977 = buffer.data(isk + 977);
    const auto *isk_978 = buffer.data(isk + 978);
    const auto *isk_979 = buffer.data(isk + 979);
    const auto *isk_981 = buffer.data(isk + 981);
    const auto *isk_982 = buffer.data(isk + 982);
    const auto *isk_983 = buffer.data(isk + 983);
    const auto *isk_984 = buffer.data(isk + 984);
    const auto *isk_986 = buffer.data(isk + 986);
    const auto *isk_987 = buffer.data(isk + 987);
    const auto *isk_988 = buffer.data(isk + 988);
    const auto *isk_989 = buffer.data(isk + 989);
    const auto *isk_990 = buffer.data(isk + 990);
    const auto *isk_992 = buffer.data(isk + 992);
    const auto *isk_993 = buffer.data(isk + 993);
    const auto *isk_994 = buffer.data(isk + 994);
    const auto *isk_995 = buffer.data(isk + 995);
    const auto *isk_996 = buffer.data(isk + 996);
    const auto *isk_997 = buffer.data(isk + 997);
    const auto *isk_999 = buffer.data(isk + 999);
    const auto *isk_1000 = buffer.data(isk + 1000);
    const auto *isk_1001 = buffer.data(isk + 1001);
    const auto *isk_1002 = buffer.data(isk + 1002);
    const auto *isk_1003 = buffer.data(isk + 1003);
    const auto *isk_1004 = buffer.data(isk + 1004);
    const auto *isk_1005 = buffer.data(isk + 1005);
    const auto *isk_1006 = buffer.data(isk + 1006);
    const auto *isk_1007 = buffer.data(isk + 1007);

#pragma omp simd aligned(t_1174, t_1175, t_1176, pa_y, pc_x, pc_y, hsl0_905, hsl1_905, \
                         isi0_732, isi0_734, isi1_732, isi1_734, isk_940, \
                         isk_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1174[k] = f_12 * isi0_732[k]
                    - f_13 * isi1_732[k]
                    + f_3 * pc_x[k] * isk_940[k];

        t_1175[k] = pa_y[k] * hsl0_905[k]
                    - f_14 * pc_y[k] * hsl1_905[k];

        t_1176[k] = f_10 * isi0_734[k]
                    - f_11 * isi1_734[k]
                    + f_3 * pc_x[k] * isk_942[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, pa_y, pc_x, pc_y, hsl0_909, hsl1_909, \
                         isi0_735, isi0_736, isi1_735, isi1_736, isk_943, \
                         isk_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_10 * isi0_735[k]
                    - f_11 * isi1_735[k]
                    + f_3 * pc_x[k] * isk_943[k];

        t_1178[k] = f_10 * isi0_736[k]
                    - f_11 * isi1_736[k]
                    + f_3 * pc_x[k] * isk_944[k];

        t_1179[k] = pa_y[k] * hsl0_909[k]
                    - f_14 * pc_y[k] * hsl1_909[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pc_x, isi0_738, isi0_739, isi0_740, isi1_738, \
                         isi1_739, isi1_740, isk_946, isk_947, \
                         isk_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_8 * isi0_738[k]
                    - f_9 * isi1_738[k]
                    + f_3 * pc_x[k] * isk_946[k];

        t_1181[k] = f_8 * isi0_739[k]
                    - f_9 * isi1_739[k]
                    + f_3 * pc_x[k] * isk_947[k];

        t_1182[k] = f_8 * isi0_740[k]
                    - f_9 * isi1_740[k]
                    + f_3 * pc_x[k] * isk_948[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, pa_y, pc_x, pc_y, hsl0_914, hsl1_914, \
                         isi0_741, isi0_743, isi1_741, isi1_743, isk_949, \
                         isk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_8 * isi0_741[k]
                    - f_9 * isi1_741[k]
                    + f_3 * pc_x[k] * isk_949[k];

        t_1184[k] = pa_y[k] * hsl0_914[k]
                    - f_14 * pc_y[k] * hsl1_914[k];

        t_1185[k] = f_6 * isi0_743[k]
                    - f_7 * isi1_743[k]
                    + f_3 * pc_x[k] * isk_951[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, pc_x, isi0_744, isi0_745, isi0_746, isi1_744, \
                         isi1_745, isi1_746, isk_952, isk_953, \
                         isk_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_6 * isi0_744[k]
                    - f_7 * isi1_744[k]
                    + f_3 * pc_x[k] * isk_952[k];

        t_1187[k] = f_6 * isi0_745[k]
                    - f_7 * isi1_745[k]
                    + f_3 * pc_x[k] * isk_953[k];

        t_1188[k] = f_6 * isi0_746[k]
                    - f_7 * isi1_746[k]
                    + f_3 * pc_x[k] * isk_954[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, pa_y, pc_x, pc_y, hsl0_920, hsl1_920, \
                         isi0_747, isi0_749, isi1_747, isi1_749, isk_955, \
                         isk_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_6 * isi0_747[k]
                    - f_7 * isi1_747[k]
                    + f_3 * pc_x[k] * isk_955[k];

        t_1190[k] = pa_y[k] * hsl0_920[k]
                    - f_14 * pc_y[k] * hsl1_920[k];

        t_1191[k] = f_4 * isi0_749[k]
                    - f_5 * isi1_749[k]
                    + f_3 * pc_x[k] * isk_957[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, pc_x, isi0_750, isi0_751, isi0_752, isi1_750, \
                         isi1_751, isi1_752, isk_958, isk_959, \
                         isk_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = f_4 * isi0_750[k]
                    - f_5 * isi1_750[k]
                    + f_3 * pc_x[k] * isk_958[k];

        t_1193[k] = f_4 * isi0_751[k]
                    - f_5 * isi1_751[k]
                    + f_3 * pc_x[k] * isk_959[k];

        t_1194[k] = f_4 * isi0_752[k]
                    - f_5 * isi1_752[k]
                    + f_3 * pc_x[k] * isk_960[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, pa_y, pc_x, pc_y, hsl0_927, hsl1_927, \
                         isi0_753, isi0_754, isi1_753, isi1_754, isk_961, isk_962, \
                         isk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = f_4 * isi0_753[k]
                    - f_5 * isi1_753[k]
                    + f_3 * pc_x[k] * isk_961[k];

        t_1196[k] = f_4 * isi0_754[k]
                    - f_5 * isi1_754[k]
                    + f_3 * pc_x[k] * isk_962[k];

        t_1197[k] = pa_y[k] * hsl0_927[k]
                    - f_14 * pc_y[k] * hsl1_927[k];

        t_1198[k] = f_3 * pc_x[k] * isk_964[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, t_1202, t_1203, t_1204, t_1205, pc_x, \
                         isk_965, isk_966, isk_967, isk_968, isk_969, isk_970, \
                         isk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_3 * pc_x[k] * isk_965[k];

        t_1200[k] = f_3 * pc_x[k] * isk_966[k];

        t_1201[k] = f_3 * pc_x[k] * isk_967[k];

        t_1202[k] = f_3 * pc_x[k] * isk_968[k];

        t_1203[k] = f_3 * pc_x[k] * isk_969[k];

        t_1204[k] = f_3 * pc_x[k] * isk_970[k];

        t_1205[k] = f_3 * pc_x[k] * isk_971[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, pa_y, pc_y, pc_z, hsl0_936, hsl0_938, \
                         hsk_712, hsk_748, hsk_750, hsl1_936, hsl1_938, \
                         isk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = pa_y[k] * hsl0_936[k]
                    + f_22 * hsk_748[k]
                    - f_14 * pc_y[k] * hsl1_936[k];

        t_1207[k] = f_19 * hsk_712[k]
                    + f_3 * pc_z[k] * isk_964[k];

        t_1208[k] = pa_y[k] * hsl0_938[k]
                    + f_0 * hsk_750[k]
                    - f_14 * pc_y[k] * hsl1_938[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, pa_y, pc_y, hsl0_939, hsl0_940, hsl0_941, \
                         hsk_751, hsk_752, hsk_753, hsl1_939, hsl1_940, \
                         hsl1_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = pa_y[k] * hsl0_939[k]
                    + f_19 * hsk_751[k]
                    - f_14 * pc_y[k] * hsl1_939[k];

        t_1210[k] = pa_y[k] * hsl0_940[k]
                    + f_18 * hsk_752[k]
                    - f_14 * pc_y[k] * hsl1_940[k];

        t_1211[k] = pa_y[k] * hsl0_941[k]
                    + f_17 * hsk_753[k]
                    - f_14 * pc_y[k] * hsl1_941[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, pa_y, pc_y, hsl0_942, hsl0_944, hsk_754, \
                         hsk_755, hsl1_942, hsl1_944, isk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = pa_y[k] * hsl0_942[k]
                    + f_16 * hsk_754[k]
                    - f_14 * pc_y[k] * hsl1_942[k];

        t_1213[k] = f_15 * hsk_755[k]
                    + f_3 * pc_y[k] * isk_971[k];

        t_1214[k] = pa_y[k] * hsl0_944[k]
                    - f_14 * pc_y[k] * hsl1_944[k];
    }

#pragma omp simd aligned(t_1215, t_1216, t_1217, t_1218, t_1219, pc_x, pc_y, isi0_756, \
                         isi0_758, isi0_759, isi1_756, isi1_758, isi1_759, isk_972, isk_974, \
                         isk_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1215[k] = f_1 * isi0_756[k]
                    - f_2 * isi1_756[k]
                    + f_3 * pc_x[k] * isk_972[k];

        t_1216[k] = f_3 * pc_y[k] * isk_972[k];

        t_1217[k] = f_20 * isi0_758[k]
                    - f_21 * isi1_758[k]
                    + f_3 * pc_x[k] * isk_974[k];

        t_1218[k] = f_12 * isi0_759[k]
                    - f_13 * isi1_759[k]
                    + f_3 * pc_x[k] * isk_975[k];

        t_1219[k] = f_3 * pc_y[k] * isk_974[k];
    }

#pragma omp simd aligned(t_1220, t_1221, t_1222, t_1223, pc_x, pc_y, isi0_761, isi0_762, \
                         isi0_763, isi1_761, isi1_762, isi1_763, isk_977, isk_978, \
                         isk_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1220[k] = f_12 * isi0_761[k]
                    - f_13 * isi1_761[k]
                    + f_3 * pc_x[k] * isk_977[k];

        t_1221[k] = f_10 * isi0_762[k]
                    - f_11 * isi1_762[k]
                    + f_3 * pc_x[k] * isk_978[k];

        t_1222[k] = f_10 * isi0_763[k]
                    - f_11 * isi1_763[k]
                    + f_3 * pc_x[k] * isk_979[k];

        t_1223[k] = f_3 * pc_y[k] * isk_977[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pc_x, isi0_765, isi0_766, isi0_767, isi1_765, \
                         isi1_766, isi1_767, isk_981, isk_982, \
                         isk_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = f_10 * isi0_765[k]
                    - f_11 * isi1_765[k]
                    + f_3 * pc_x[k] * isk_981[k];

        t_1225[k] = f_8 * isi0_766[k]
                    - f_9 * isi1_766[k]
                    + f_3 * pc_x[k] * isk_982[k];

        t_1226[k] = f_8 * isi0_767[k]
                    - f_9 * isi1_767[k]
                    + f_3 * pc_x[k] * isk_983[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, t_1230, pc_x, pc_y, isi0_768, isi0_770, \
                         isi0_771, isi1_768, isi1_770, isi1_771, isk_981, isk_984, isk_986, \
                         isk_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_8 * isi0_768[k]
                    - f_9 * isi1_768[k]
                    + f_3 * pc_x[k] * isk_984[k];

        t_1228[k] = f_3 * pc_y[k] * isk_981[k];

        t_1229[k] = f_8 * isi0_770[k]
                    - f_9 * isi1_770[k]
                    + f_3 * pc_x[k] * isk_986[k];

        t_1230[k] = f_6 * isi0_771[k]
                    - f_7 * isi1_771[k]
                    + f_3 * pc_x[k] * isk_987[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pc_x, pc_y, isi0_772, isi0_773, \
                         isi0_774, isi1_772, isi1_773, isi1_774, isk_986, isk_988, isk_989, \
                         isk_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = f_6 * isi0_772[k]
                    - f_7 * isi1_772[k]
                    + f_3 * pc_x[k] * isk_988[k];

        t_1232[k] = f_6 * isi0_773[k]
                    - f_7 * isi1_773[k]
                    + f_3 * pc_x[k] * isk_989[k];

        t_1233[k] = f_6 * isi0_774[k]
                    - f_7 * isi1_774[k]
                    + f_3 * pc_x[k] * isk_990[k];

        t_1234[k] = f_3 * pc_y[k] * isk_986[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pc_x, isi0_776, isi0_777, isi0_778, isi1_776, \
                         isi1_777, isi1_778, isk_992, isk_993, \
                         isk_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_6 * isi0_776[k]
                    - f_7 * isi1_776[k]
                    + f_3 * pc_x[k] * isk_992[k];

        t_1236[k] = f_4 * isi0_777[k]
                    - f_5 * isi1_777[k]
                    + f_3 * pc_x[k] * isk_993[k];

        t_1237[k] = f_4 * isi0_778[k]
                    - f_5 * isi1_778[k]
                    + f_3 * pc_x[k] * isk_994[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, t_1241, pc_x, pc_y, isi0_779, isi0_780, \
                         isi0_781, isi1_779, isi1_780, isi1_781, isk_992, isk_995, isk_996, \
                         isk_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_4 * isi0_779[k]
                    - f_5 * isi1_779[k]
                    + f_3 * pc_x[k] * isk_995[k];

        t_1239[k] = f_4 * isi0_780[k]
                    - f_5 * isi1_780[k]
                    + f_3 * pc_x[k] * isk_996[k];

        t_1240[k] = f_4 * isi0_781[k]
                    - f_5 * isi1_781[k]
                    + f_3 * pc_x[k] * isk_997[k];

        t_1241[k] = f_3 * pc_y[k] * isk_992[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, t_1246, t_1247, pc_x, isi0_783, \
                         isi1_783, isk_999, isk_1000, isk_1001, isk_1002, isk_1003, \
                         isk_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = f_4 * isi0_783[k]
                    - f_5 * isi1_783[k]
                    + f_3 * pc_x[k] * isk_999[k];

        t_1243[k] = f_3 * pc_x[k] * isk_1000[k];

        t_1244[k] = f_3 * pc_x[k] * isk_1001[k];

        t_1245[k] = f_3 * pc_x[k] * isk_1002[k];

        t_1246[k] = f_3 * pc_x[k] * isk_1003[k];

        t_1247[k] = f_3 * pc_x[k] * isk_1004[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, t_1252, pc_x, pc_y, isi0_777, \
                         isi0_778, isi1_777, isi1_778, isk_1000, isk_1001, isk_1005, isk_1006, \
                         isk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = f_3 * pc_x[k] * isk_1005[k];

        t_1249[k] = f_3 * pc_x[k] * isk_1006[k];

        t_1250[k] = f_3 * pc_x[k] * isk_1007[k];

        t_1251[k] = f_1 * isi0_777[k]
                    - f_2 * isi1_777[k]
                    + f_3 * pc_y[k] * isk_1000[k];

        t_1252[k] = f_20 * isi0_778[k]
                    - f_21 * isi1_778[k]
                    + f_3 * pc_y[k] * isk_1001[k];
    }

#pragma omp simd aligned(t_1253, t_1254, t_1255, pc_y, isi0_779, isi0_780, isi0_781, isi1_779, \
                         isi1_780, isi1_781, isk_1002, isk_1003, \
                         isk_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1253[k] = f_12 * isi0_779[k]
                    - f_13 * isi1_779[k]
                    + f_3 * pc_y[k] * isk_1002[k];

        t_1254[k] = f_10 * isi0_780[k]
                    - f_11 * isi1_780[k]
                    + f_3 * pc_y[k] * isk_1003[k];

        t_1255[k] = f_8 * isi0_781[k]
                    - f_9 * isi1_781[k]
                    + f_3 * pc_y[k] * isk_1004[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pc_y, pc_z, hsk_755, isi0_782, \
                         isi0_783, isi1_782, isi1_783, isk_1005, isk_1006, \
                         isk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_6 * isi0_782[k]
                    - f_7 * isi1_782[k]
                    + f_3 * pc_y[k] * isk_1005[k];

        t_1257[k] = f_4 * isi0_783[k]
                    - f_5 * isi1_783[k]
                    + f_3 * pc_y[k] * isk_1006[k];

        t_1258[k] = f_3 * pc_y[k] * isk_1007[k];

        t_1259[k] = f_0 * hsk_755[k]
                    + f_1 * isi0_783[k]
                    - f_2 * isi1_783[k]
                    + f_3 * pc_z[k] * isk_1007[k];
    }
}

auto
compute_prim_isl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t hsl0, const size_t hsk,
                                                   const size_t hsl1, const size_t isi0,
                                                   const size_t isi1, const size_t isk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_isl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isk, ncols, gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, hsl0, hsk,
                                                              hsl1, isi0, isi1, isk, ncols,
                                                              gamma, p, q);

    compute_prim_isl_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, hsl0,
                                                               hsk, hsl1, isi0, isi1, isk,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
