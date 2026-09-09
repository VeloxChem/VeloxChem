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


#include "SimdThreeCenterElectronRepulsionVrrRecKSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);

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

    const auto *isl0_0 = buffer.data(isl0 + 0);
    const auto *isl0_3 = buffer.data(isl0 + 3);
    const auto *isl0_5 = buffer.data(isl0 + 5);
    const auto *isl0_6 = buffer.data(isl0 + 6);
    const auto *isl0_9 = buffer.data(isl0 + 9);
    const auto *isl0_10 = buffer.data(isl0 + 10);
    const auto *isl0_14 = buffer.data(isl0 + 14);
    const auto *isl0_15 = buffer.data(isl0 + 15);
    const auto *isl0_20 = buffer.data(isl0 + 20);
    const auto *isl0_21 = buffer.data(isl0 + 21);
    const auto *isl0_27 = buffer.data(isl0 + 27);
    const auto *isl0_36 = buffer.data(isl0 + 36);
    const auto *isl0_44 = buffer.data(isl0 + 44);

    const auto *isk_0 = buffer.data(isk + 0);
    const auto *isk_1 = buffer.data(isk + 1);
    const auto *isk_2 = buffer.data(isk + 2);
    const auto *isk_3 = buffer.data(isk + 3);
    const auto *isk_5 = buffer.data(isk + 5);
    const auto *isk_6 = buffer.data(isk + 6);
    const auto *isk_9 = buffer.data(isk + 9);
    const auto *isk_10 = buffer.data(isk + 10);
    const auto *isk_14 = buffer.data(isk + 14);
    const auto *isk_15 = buffer.data(isk + 15);
    const auto *isk_20 = buffer.data(isk + 20);
    const auto *isk_28 = buffer.data(isk + 28);
    const auto *isk_30 = buffer.data(isk + 30);
    const auto *isk_31 = buffer.data(isk + 31);
    const auto *isk_32 = buffer.data(isk + 32);
    const auto *isk_33 = buffer.data(isk + 33);
    const auto *isk_35 = buffer.data(isk + 35);
    const auto *isk_64 = buffer.data(isk + 64);
    const auto *isk_66 = buffer.data(isk + 66);
    const auto *isk_67 = buffer.data(isk + 67);
    const auto *isk_68 = buffer.data(isk + 68);
    const auto *isk_69 = buffer.data(isk + 69);
    const auto *isk_70 = buffer.data(isk + 70);
    const auto *isk_71 = buffer.data(isk + 71);
    const auto *isk_100 = buffer.data(isk + 100);
    const auto *isk_101 = buffer.data(isk + 101);
    const auto *isk_102 = buffer.data(isk + 102);
    const auto *isk_103 = buffer.data(isk + 103);
    const auto *isk_104 = buffer.data(isk + 104);
    const auto *isk_105 = buffer.data(isk + 105);
    const auto *isk_107 = buffer.data(isk + 107);

    const auto *isl1_0 = buffer.data(isl1 + 0);
    const auto *isl1_3 = buffer.data(isl1 + 3);
    const auto *isl1_5 = buffer.data(isl1 + 5);
    const auto *isl1_6 = buffer.data(isl1 + 6);
    const auto *isl1_9 = buffer.data(isl1 + 9);
    const auto *isl1_10 = buffer.data(isl1 + 10);
    const auto *isl1_14 = buffer.data(isl1 + 14);
    const auto *isl1_15 = buffer.data(isl1 + 15);
    const auto *isl1_20 = buffer.data(isl1 + 20);
    const auto *isl1_21 = buffer.data(isl1 + 21);
    const auto *isl1_27 = buffer.data(isl1 + 27);
    const auto *isl1_36 = buffer.data(isl1 + 36);
    const auto *isl1_44 = buffer.data(isl1 + 44);

    const auto *ksi0_0 = buffer.data(ksi0 + 0);
    const auto *ksi0_1 = buffer.data(ksi0 + 1);
    const auto *ksi0_2 = buffer.data(ksi0 + 2);
    const auto *ksi0_3 = buffer.data(ksi0 + 3);
    const auto *ksi0_5 = buffer.data(ksi0 + 5);
    const auto *ksi0_6 = buffer.data(ksi0 + 6);
    const auto *ksi0_8 = buffer.data(ksi0 + 8);
    const auto *ksi0_9 = buffer.data(ksi0 + 9);
    const auto *ksi0_10 = buffer.data(ksi0 + 10);
    const auto *ksi0_12 = buffer.data(ksi0 + 12);
    const auto *ksi0_13 = buffer.data(ksi0 + 13);
    const auto *ksi0_14 = buffer.data(ksi0 + 14);
    const auto *ksi0_21 = buffer.data(ksi0 + 21);
    const auto *ksi0_23 = buffer.data(ksi0 + 23);
    const auto *ksi0_24 = buffer.data(ksi0 + 24);
    const auto *ksi0_25 = buffer.data(ksi0 + 25);
    const auto *ksi0_26 = buffer.data(ksi0 + 26);
    const auto *ksi0_27 = buffer.data(ksi0 + 27);
    const auto *ksi0_31 = buffer.data(ksi0 + 31);
    const auto *ksi0_34 = buffer.data(ksi0 + 34);
    const auto *ksi0_35 = buffer.data(ksi0 + 35);
    const auto *ksi0_38 = buffer.data(ksi0 + 38);
    const auto *ksi0_39 = buffer.data(ksi0 + 39);
    const auto *ksi0_40 = buffer.data(ksi0 + 40);
    const auto *ksi0_49 = buffer.data(ksi0 + 49);
    const auto *ksi0_50 = buffer.data(ksi0 + 50);
    const auto *ksi0_51 = buffer.data(ksi0 + 51);
    const auto *ksi0_52 = buffer.data(ksi0 + 52);
    const auto *ksi0_53 = buffer.data(ksi0 + 53);
    const auto *ksi0_58 = buffer.data(ksi0 + 58);
    const auto *ksi0_60 = buffer.data(ksi0 + 60);
    const auto *ksi0_61 = buffer.data(ksi0 + 61);
    const auto *ksi0_63 = buffer.data(ksi0 + 63);
    const auto *ksi0_64 = buffer.data(ksi0 + 64);
    const auto *ksi0_65 = buffer.data(ksi0 + 65);
    const auto *ksi0_67 = buffer.data(ksi0 + 67);
    const auto *ksi0_68 = buffer.data(ksi0 + 68);
    const auto *ksi0_69 = buffer.data(ksi0 + 69);
    const auto *ksi0_70 = buffer.data(ksi0 + 70);
    const auto *ksi0_78 = buffer.data(ksi0 + 78);
    const auto *ksi0_79 = buffer.data(ksi0 + 79);

    const auto *ksi1_0 = buffer.data(ksi1 + 0);
    const auto *ksi1_1 = buffer.data(ksi1 + 1);
    const auto *ksi1_2 = buffer.data(ksi1 + 2);
    const auto *ksi1_3 = buffer.data(ksi1 + 3);
    const auto *ksi1_5 = buffer.data(ksi1 + 5);
    const auto *ksi1_6 = buffer.data(ksi1 + 6);
    const auto *ksi1_8 = buffer.data(ksi1 + 8);
    const auto *ksi1_9 = buffer.data(ksi1 + 9);
    const auto *ksi1_10 = buffer.data(ksi1 + 10);
    const auto *ksi1_12 = buffer.data(ksi1 + 12);
    const auto *ksi1_13 = buffer.data(ksi1 + 13);
    const auto *ksi1_14 = buffer.data(ksi1 + 14);
    const auto *ksi1_21 = buffer.data(ksi1 + 21);
    const auto *ksi1_23 = buffer.data(ksi1 + 23);
    const auto *ksi1_24 = buffer.data(ksi1 + 24);
    const auto *ksi1_25 = buffer.data(ksi1 + 25);
    const auto *ksi1_26 = buffer.data(ksi1 + 26);
    const auto *ksi1_27 = buffer.data(ksi1 + 27);
    const auto *ksi1_31 = buffer.data(ksi1 + 31);
    const auto *ksi1_34 = buffer.data(ksi1 + 34);
    const auto *ksi1_35 = buffer.data(ksi1 + 35);
    const auto *ksi1_38 = buffer.data(ksi1 + 38);
    const auto *ksi1_39 = buffer.data(ksi1 + 39);
    const auto *ksi1_40 = buffer.data(ksi1 + 40);
    const auto *ksi1_49 = buffer.data(ksi1 + 49);
    const auto *ksi1_50 = buffer.data(ksi1 + 50);
    const auto *ksi1_51 = buffer.data(ksi1 + 51);
    const auto *ksi1_52 = buffer.data(ksi1 + 52);
    const auto *ksi1_53 = buffer.data(ksi1 + 53);
    const auto *ksi1_58 = buffer.data(ksi1 + 58);
    const auto *ksi1_60 = buffer.data(ksi1 + 60);
    const auto *ksi1_61 = buffer.data(ksi1 + 61);
    const auto *ksi1_63 = buffer.data(ksi1 + 63);
    const auto *ksi1_64 = buffer.data(ksi1 + 64);
    const auto *ksi1_65 = buffer.data(ksi1 + 65);
    const auto *ksi1_67 = buffer.data(ksi1 + 67);
    const auto *ksi1_68 = buffer.data(ksi1 + 68);
    const auto *ksi1_69 = buffer.data(ksi1 + 69);
    const auto *ksi1_70 = buffer.data(ksi1 + 70);
    const auto *ksi1_78 = buffer.data(ksi1 + 78);
    const auto *ksi1_79 = buffer.data(ksi1 + 79);

    const auto *ksk_0 = buffer.data(ksk + 0);
    const auto *ksk_1 = buffer.data(ksk + 1);
    const auto *ksk_2 = buffer.data(ksk + 2);
    const auto *ksk_3 = buffer.data(ksk + 3);
    const auto *ksk_5 = buffer.data(ksk + 5);
    const auto *ksk_6 = buffer.data(ksk + 6);
    const auto *ksk_8 = buffer.data(ksk + 8);
    const auto *ksk_9 = buffer.data(ksk + 9);
    const auto *ksk_10 = buffer.data(ksk + 10);
    const auto *ksk_12 = buffer.data(ksk + 12);
    const auto *ksk_13 = buffer.data(ksk + 13);
    const auto *ksk_14 = buffer.data(ksk + 14);
    const auto *ksk_15 = buffer.data(ksk + 15);
    const auto *ksk_17 = buffer.data(ksk + 17);
    const auto *ksk_18 = buffer.data(ksk + 18);
    const auto *ksk_19 = buffer.data(ksk + 19);
    const auto *ksk_20 = buffer.data(ksk + 20);
    const auto *ksk_21 = buffer.data(ksk + 21);
    const auto *ksk_27 = buffer.data(ksk + 27);
    const auto *ksk_28 = buffer.data(ksk + 28);
    const auto *ksk_30 = buffer.data(ksk + 30);
    const auto *ksk_31 = buffer.data(ksk + 31);
    const auto *ksk_32 = buffer.data(ksk + 32);
    const auto *ksk_33 = buffer.data(ksk + 33);
    const auto *ksk_34 = buffer.data(ksk + 34);
    const auto *ksk_35 = buffer.data(ksk + 35);
    const auto *ksk_36 = buffer.data(ksk + 36);
    const auto *ksk_37 = buffer.data(ksk + 37);
    const auto *ksk_39 = buffer.data(ksk + 39);
    const auto *ksk_41 = buffer.data(ksk + 41);
    const auto *ksk_42 = buffer.data(ksk + 42);
    const auto *ksk_43 = buffer.data(ksk + 43);
    const auto *ksk_45 = buffer.data(ksk + 45);
    const auto *ksk_46 = buffer.data(ksk + 46);
    const auto *ksk_47 = buffer.data(ksk + 47);
    const auto *ksk_48 = buffer.data(ksk + 48);
    const auto *ksk_50 = buffer.data(ksk + 50);
    const auto *ksk_51 = buffer.data(ksk + 51);
    const auto *ksk_52 = buffer.data(ksk + 52);
    const auto *ksk_53 = buffer.data(ksk + 53);
    const auto *ksk_54 = buffer.data(ksk + 54);
    const auto *ksk_56 = buffer.data(ksk + 56);
    const auto *ksk_57 = buffer.data(ksk + 57);
    const auto *ksk_64 = buffer.data(ksk + 64);
    const auto *ksk_65 = buffer.data(ksk + 65);
    const auto *ksk_66 = buffer.data(ksk + 66);
    const auto *ksk_67 = buffer.data(ksk + 67);
    const auto *ksk_68 = buffer.data(ksk + 68);
    const auto *ksk_69 = buffer.data(ksk + 69);
    const auto *ksk_70 = buffer.data(ksk + 70);
    const auto *ksk_71 = buffer.data(ksk + 71);
    const auto *ksk_72 = buffer.data(ksk + 72);
    const auto *ksk_74 = buffer.data(ksk + 74);
    const auto *ksk_76 = buffer.data(ksk + 76);
    const auto *ksk_77 = buffer.data(ksk + 77);
    const auto *ksk_79 = buffer.data(ksk + 79);
    const auto *ksk_80 = buffer.data(ksk + 80);
    const auto *ksk_81 = buffer.data(ksk + 81);
    const auto *ksk_83 = buffer.data(ksk + 83);
    const auto *ksk_84 = buffer.data(ksk + 84);
    const auto *ksk_85 = buffer.data(ksk + 85);
    const auto *ksk_86 = buffer.data(ksk + 86);
    const auto *ksk_88 = buffer.data(ksk + 88);
    const auto *ksk_89 = buffer.data(ksk + 89);
    const auto *ksk_90 = buffer.data(ksk + 90);
    const auto *ksk_91 = buffer.data(ksk + 91);
    const auto *ksk_92 = buffer.data(ksk + 92);
    const auto *ksk_99 = buffer.data(ksk + 99);
    const auto *ksk_100 = buffer.data(ksk + 100);
    const auto *ksk_101 = buffer.data(ksk + 101);
    const auto *ksk_102 = buffer.data(ksk + 102);
    const auto *ksk_103 = buffer.data(ksk + 103);
    const auto *ksk_104 = buffer.data(ksk + 104);
    const auto *ksk_105 = buffer.data(ksk + 105);
    const auto *ksk_107 = buffer.data(ksk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, isk_0, ksi0_0, \
                         ksi1_0, ksk_0, ksk_1, ksk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * isk_0[k]
                 + f_1 * ksi0_0[k]
                 - f_2 * ksi1_0[k]
                 + f_3 * pc_x[k] * ksk_0[k];

        t_1[k] = f_3 * pc_y[k] * ksk_0[k];

        t_2[k] = f_3 * pc_z[k] * ksk_0[k];

        t_3[k] = f_4 * ksi0_0[k]
                 - f_5 * ksi1_0[k]
                 + f_3 * pc_y[k] * ksk_1[k];

        t_4[k] = f_3 * pc_y[k] * ksk_2[k];

        t_5[k] = f_4 * ksi0_0[k]
                 - f_5 * ksi1_0[k]
                 + f_3 * pc_z[k] * ksk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, ksi0_1, ksi0_2, ksi0_3, ksi1_1, \
                         ksi1_2, ksi1_3, ksk_3, ksk_5, ksk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * ksi0_1[k]
                 - f_7 * ksi1_1[k]
                 + f_3 * pc_y[k] * ksk_3[k];

        t_7[k] = f_3 * pc_z[k] * ksk_3[k];

        t_8[k] = f_3 * pc_y[k] * ksk_5[k];

        t_9[k] = f_6 * ksi0_2[k]
                 - f_7 * ksi1_2[k]
                 + f_3 * pc_z[k] * ksk_5[k];

        t_10[k] = f_8 * ksi0_3[k]
                  - f_9 * ksi1_3[k]
                  + f_3 * pc_y[k] * ksk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pc_y, pc_z, ksi0_5, ksi0_6, \
                         ksi1_5, ksi1_6, ksk_6, ksk_8, ksk_9, ksk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * ksk_6[k];

        t_12[k] = f_4 * ksi0_5[k]
                  - f_5 * ksi1_5[k]
                  + f_3 * pc_y[k] * ksk_8[k];

        t_13[k] = f_3 * pc_y[k] * ksk_9[k];

        t_14[k] = f_8 * ksi0_5[k]
                  - f_9 * ksi1_5[k]
                  + f_3 * pc_z[k] * ksk_9[k];

        t_15[k] = f_10 * ksi0_6[k]
                  - f_11 * ksi1_6[k]
                  + f_3 * pc_y[k] * ksk_10[k];

        t_16[k] = f_3 * pc_z[k] * ksk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_y, pc_z, ksi0_8, ksi0_9, ksi1_8, ksi1_9, \
                         ksk_12, ksk_13, ksk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * ksi0_8[k]
                  - f_7 * ksi1_8[k]
                  + f_3 * pc_y[k] * ksk_12[k];

        t_18[k] = f_4 * ksi0_9[k]
                  - f_5 * ksi1_9[k]
                  + f_3 * pc_y[k] * ksk_13[k];

        t_19[k] = f_3 * pc_y[k] * ksk_14[k];

        t_20[k] = f_10 * ksi0_9[k]
                  - f_11 * ksi1_9[k]
                  + f_3 * pc_z[k] * ksk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, ksi0_10, ksi0_12, ksi0_13, \
                         ksi1_10, ksi1_12, ksi1_13, ksk_15, ksk_17, \
                         ksk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_12 * ksi0_10[k]
                  - f_13 * ksi1_10[k]
                  + f_3 * pc_y[k] * ksk_15[k];

        t_22[k] = f_3 * pc_z[k] * ksk_15[k];

        t_23[k] = f_8 * ksi0_12[k]
                  - f_9 * ksi1_12[k]
                  + f_3 * pc_y[k] * ksk_17[k];

        t_24[k] = f_6 * ksi0_13[k]
                  - f_7 * ksi1_13[k]
                  + f_3 * pc_y[k] * ksk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, pc_z, isk_28, ksi0_14, \
                         ksi1_14, ksk_19, ksk_20, ksk_21, ksk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * ksi0_14[k]
                  - f_5 * ksi1_14[k]
                  + f_3 * pc_y[k] * ksk_19[k];

        t_26[k] = f_3 * pc_y[k] * ksk_20[k];

        t_27[k] = f_12 * ksi0_14[k]
                  - f_13 * ksi1_14[k]
                  + f_3 * pc_z[k] * ksk_20[k];

        t_28[k] = f_0 * isk_28[k]
                  + f_3 * pc_x[k] * ksk_28[k];

        t_29[k] = f_3 * pc_z[k] * ksk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, isk_30, isk_31, isk_32, \
                         isk_33, ksk_27, ksk_30, ksk_31, ksk_32, \
                         ksk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * isk_30[k]
                  + f_3 * pc_x[k] * ksk_30[k];

        t_31[k] = f_0 * isk_31[k]
                  + f_3 * pc_x[k] * ksk_31[k];

        t_32[k] = f_0 * isk_32[k]
                  + f_3 * pc_x[k] * ksk_32[k];

        t_33[k] = f_0 * isk_33[k]
                  + f_3 * pc_x[k] * ksk_33[k];

        t_34[k] = f_3 * pc_y[k] * ksk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, isk_35, ksi0_21, ksi0_23, \
                         ksi1_21, ksi1_23, ksk_28, ksk_30, ksk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * isk_35[k]
                  + f_3 * pc_x[k] * ksk_35[k];

        t_36[k] = f_1 * ksi0_21[k]
                  - f_2 * ksi1_21[k]
                  + f_3 * pc_y[k] * ksk_28[k];

        t_37[k] = f_3 * pc_z[k] * ksk_28[k];

        t_38[k] = f_12 * ksi0_23[k]
                  - f_13 * ksi1_23[k]
                  + f_3 * pc_y[k] * ksk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pc_y, ksi0_24, ksi0_25, ksi0_26, ksi1_24, ksi1_25, \
                         ksi1_26, ksk_31, ksk_32, ksk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * ksi0_24[k]
                  - f_11 * ksi1_24[k]
                  + f_3 * pc_y[k] * ksk_31[k];

        t_40[k] = f_8 * ksi0_25[k]
                  - f_9 * ksi1_25[k]
                  + f_3 * pc_y[k] * ksk_32[k];

        t_41[k] = f_6 * ksi0_26[k]
                  - f_7 * ksi1_26[k]
                  + f_3 * pc_y[k] * ksk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pc_y, pc_z, isl0_0, isk_0, \
                         isl1_0, ksi0_27, ksi1_27, ksk_34, ksk_35, \
                         ksk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * ksi0_27[k]
                  - f_5 * ksi1_27[k]
                  + f_3 * pc_y[k] * ksk_34[k];

        t_43[k] = f_3 * pc_y[k] * ksk_35[k];

        t_44[k] = f_1 * ksi0_27[k]
                  - f_2 * ksi1_27[k]
                  + f_3 * pc_z[k] * ksk_35[k];

        t_45[k] = pa_y[k] * isl0_0[k]
                  - f_14 * pc_y[k] * isl1_0[k];

        t_46[k] = f_15 * isk_0[k]
                  + f_3 * pc_y[k] * ksk_36[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_y, pc_z, isl0_3, isl0_5, isk_1, \
                         isl1_3, isl1_5, ksk_36, ksk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * ksk_36[k];

        t_48[k] = pa_y[k] * isl0_3[k]
                  + f_16 * isk_1[k]
                  - f_14 * pc_y[k] * isl1_3[k];

        t_49[k] = f_3 * pc_z[k] * ksk_37[k];

        t_50[k] = pa_y[k] * isl0_5[k]
                  - f_14 * pc_y[k] * isl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pc_y, pc_z, isl0_6, isl0_9, isk_3, \
                         isk_5, isl1_6, isl1_9, ksk_39, ksk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * isl0_6[k]
                  + f_17 * isk_3[k]
                  - f_14 * pc_y[k] * isl1_6[k];

        t_52[k] = f_3 * pc_z[k] * ksk_39[k];

        t_53[k] = f_15 * isk_5[k]
                  + f_3 * pc_y[k] * ksk_41[k];

        t_54[k] = pa_y[k] * isl0_9[k]
                  - f_14 * pc_y[k] * isl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pc_y, pc_z, isl0_10, isk_6, isk_9, \
                         isl1_10, ksi0_31, ksi1_31, ksk_42, ksk_43, \
                         ksk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * isl0_10[k]
                  + f_18 * isk_6[k]
                  - f_14 * pc_y[k] * isl1_10[k];

        t_56[k] = f_3 * pc_z[k] * ksk_42[k];

        t_57[k] = f_4 * ksi0_31[k]
                  - f_5 * ksi1_31[k]
                  + f_3 * pc_z[k] * ksk_43[k];

        t_58[k] = f_15 * isk_9[k]
                  + f_3 * pc_y[k] * ksk_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, isl0_14, isl0_15, isk_10, \
                         isl1_14, isl1_15, ksi0_34, ksi1_34, ksk_46, \
                         ksk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * isl0_14[k]
                  - f_14 * pc_y[k] * isl1_14[k];

        t_60[k] = pa_y[k] * isl0_15[k]
                  + f_19 * isk_10[k]
                  - f_14 * pc_y[k] * isl1_15[k];

        t_61[k] = f_3 * pc_z[k] * ksk_46[k];

        t_62[k] = f_4 * ksi0_34[k]
                  - f_5 * ksi1_34[k]
                  + f_3 * pc_z[k] * ksk_47[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pc_y, pc_z, isl0_20, isk_14, isl1_20, \
                         ksi0_35, ksi1_35, ksk_48, ksk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_6 * ksi0_35[k]
                  - f_7 * ksi1_35[k]
                  + f_3 * pc_z[k] * ksk_48[k];

        t_64[k] = f_15 * isk_14[k]
                  + f_3 * pc_y[k] * ksk_50[k];

        t_65[k] = pa_y[k] * isl0_20[k]
                  - f_14 * pc_y[k] * isl1_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pc_y, pc_z, isl0_21, isk_15, isl1_21, \
                         ksi0_38, ksi1_38, ksk_51, ksk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * isl0_21[k]
                  + f_20 * isk_15[k]
                  - f_14 * pc_y[k] * isl1_21[k];

        t_67[k] = f_3 * pc_z[k] * ksk_51[k];

        t_68[k] = f_4 * ksi0_38[k]
                  - f_5 * ksi1_38[k]
                  + f_3 * pc_z[k] * ksk_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_y, pc_z, isk_20, ksi0_39, ksi0_40, ksi1_39, \
                         ksi1_40, ksk_53, ksk_54, ksk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * ksi0_39[k]
                  - f_7 * ksi1_39[k]
                  + f_3 * pc_z[k] * ksk_53[k];

        t_70[k] = f_8 * ksi0_40[k]
                  - f_9 * ksi1_40[k]
                  + f_3 * pc_z[k] * ksk_54[k];

        t_71[k] = f_15 * isk_20[k]
                  + f_3 * pc_y[k] * ksk_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, pc_z, isl0_27, isk_64, \
                         isk_66, isl1_27, ksk_57, ksk_64, ksk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * isl0_27[k]
                  - f_14 * pc_y[k] * isl1_27[k];

        t_73[k] = f_20 * isk_64[k]
                  + f_3 * pc_x[k] * ksk_64[k];

        t_74[k] = f_3 * pc_z[k] * ksk_57[k];

        t_75[k] = f_20 * isk_66[k]
                  + f_3 * pc_x[k] * ksk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, isk_67, isk_68, isk_69, isk_70, \
                         isk_71, ksk_67, ksk_68, ksk_69, ksk_70, \
                         ksk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_20 * isk_67[k]
                  + f_3 * pc_x[k] * ksk_67[k];

        t_77[k] = f_20 * isk_68[k]
                  + f_3 * pc_x[k] * ksk_68[k];

        t_78[k] = f_20 * isk_69[k]
                  + f_3 * pc_x[k] * ksk_69[k];

        t_79[k] = f_20 * isk_70[k]
                  + f_3 * pc_x[k] * ksk_70[k];

        t_80[k] = f_20 * isk_71[k]
                  + f_3 * pc_x[k] * ksk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_y, pc_z, isk_28, ksi0_49, ksi0_50, \
                         ksi1_49, ksi1_50, ksk_64, ksk_65, ksk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_15 * isk_28[k]
                  + f_1 * ksi0_49[k]
                  - f_2 * ksi1_49[k]
                  + f_3 * pc_y[k] * ksk_64[k];

        t_82[k] = f_3 * pc_z[k] * ksk_64[k];

        t_83[k] = f_4 * ksi0_49[k]
                  - f_5 * ksi1_49[k]
                  + f_3 * pc_z[k] * ksk_65[k];

        t_84[k] = f_6 * ksi0_50[k]
                  - f_7 * ksi1_50[k]
                  + f_3 * pc_z[k] * ksk_66[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_z, ksi0_51, ksi0_52, ksi0_53, ksi1_51, ksi1_52, \
                         ksi1_53, ksk_67, ksk_68, ksk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_8 * ksi0_51[k]
                  - f_9 * ksi1_51[k]
                  + f_3 * pc_z[k] * ksk_67[k];

        t_86[k] = f_10 * ksi0_52[k]
                  - f_11 * ksi1_52[k]
                  + f_3 * pc_z[k] * ksk_68[k];

        t_87[k] = f_12 * ksi0_53[k]
                  - f_13 * ksi1_53[k]
                  + f_3 * pc_z[k] * ksk_69[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pc_y, pc_z, isl0_0, isl0_44, \
                         isk_35, isl1_0, isl1_44, ksk_71, ksk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * isk_35[k]
                  + f_3 * pc_y[k] * ksk_71[k];

        t_89[k] = pa_y[k] * isl0_44[k]
                  - f_14 * pc_y[k] * isl1_44[k];

        t_90[k] = pa_z[k] * isl0_0[k]
                  - f_14 * pc_z[k] * isl1_0[k];

        t_91[k] = f_3 * pc_y[k] * ksk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_y, pc_z, isl0_3, isl0_5, isk_0, \
                         isk_2, isl1_3, isl1_5, ksk_72, ksk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * isk_0[k]
                  + f_3 * pc_z[k] * ksk_72[k];

        t_93[k] = pa_z[k] * isl0_3[k]
                  - f_14 * pc_z[k] * isl1_3[k];

        t_94[k] = f_3 * pc_y[k] * ksk_74[k];

        t_95[k] = pa_z[k] * isl0_5[k]
                  + f_16 * isk_2[k]
                  - f_14 * pc_z[k] * isl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pc_y, pc_z, isl0_6, isl0_9, isk_5, \
                         isl1_6, isl1_9, ksi0_58, ksi1_58, ksk_76, \
                         ksk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * isl0_6[k]
                  - f_14 * pc_z[k] * isl1_6[k];

        t_97[k] = f_4 * ksi0_58[k]
                  - f_5 * ksi1_58[k]
                  + f_3 * pc_y[k] * ksk_76[k];

        t_98[k] = f_3 * pc_y[k] * ksk_77[k];

        t_99[k] = pa_z[k] * isl0_9[k]
                  + f_17 * isk_5[k]
                  - f_14 * pc_z[k] * isl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pc_y, pc_z, isl0_10, isl1_10, \
                         ksi0_60, ksi0_61, ksi1_60, ksi1_61, ksk_79, ksk_80, \
                         ksk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * isl0_10[k]
                   - f_14 * pc_z[k] * isl1_10[k];

        t_101[k] = f_6 * ksi0_60[k]
                   - f_7 * ksi1_60[k]
                   + f_3 * pc_y[k] * ksk_79[k];

        t_102[k] = f_4 * ksi0_61[k]
                   - f_5 * ksi1_61[k]
                   + f_3 * pc_y[k] * ksk_80[k];

        t_103[k] = f_3 * pc_y[k] * ksk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pc_y, pc_z, isl0_14, isl0_15, isk_9, \
                         isl1_14, isl1_15, ksi0_63, ksi1_63, ksk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_z[k] * isl0_14[k]
                   + f_18 * isk_9[k]
                   - f_14 * pc_z[k] * isl1_14[k];

        t_105[k] = pa_z[k] * isl0_15[k]
                   - f_14 * pc_z[k] * isl1_15[k];

        t_106[k] = f_8 * ksi0_63[k]
                   - f_9 * ksi1_63[k]
                   + f_3 * pc_y[k] * ksk_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, ksi0_64, ksi0_65, ksi1_64, ksi1_65, \
                         ksk_84, ksk_85, ksk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * ksi0_64[k]
                   - f_7 * ksi1_64[k]
                   + f_3 * pc_y[k] * ksk_84[k];

        t_108[k] = f_4 * ksi0_65[k]
                   - f_5 * ksi1_65[k]
                   + f_3 * pc_y[k] * ksk_85[k];

        t_109[k] = f_3 * pc_y[k] * ksk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pc_y, pc_z, isl0_20, isl0_21, isk_14, \
                         isl1_20, isl1_21, ksi0_67, ksi1_67, ksk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * isl0_20[k]
                   + f_19 * isk_14[k]
                   - f_14 * pc_z[k] * isl1_20[k];

        t_111[k] = pa_z[k] * isl0_21[k]
                   - f_14 * pc_z[k] * isl1_21[k];

        t_112[k] = f_10 * ksi0_67[k]
                   - f_11 * ksi1_67[k]
                   + f_3 * pc_y[k] * ksk_88[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_y, ksi0_68, ksi0_69, ksi0_70, ksi1_68, \
                         ksi1_69, ksi1_70, ksk_89, ksk_90, ksk_91, \
                         ksk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * ksi0_68[k]
                   - f_9 * ksi1_68[k]
                   + f_3 * pc_y[k] * ksk_89[k];

        t_114[k] = f_6 * ksi0_69[k]
                   - f_7 * ksi1_69[k]
                   + f_3 * pc_y[k] * ksk_90[k];

        t_115[k] = f_4 * ksi0_70[k]
                   - f_5 * ksi1_70[k]
                   + f_3 * pc_y[k] * ksk_91[k];

        t_116[k] = f_3 * pc_y[k] * ksk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pc_x, pc_z, isl0_27, isk_20, \
                         isk_100, isk_101, isk_102, isl1_27, ksk_100, ksk_101, \
                         ksk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_z[k] * isl0_27[k]
                   + f_20 * isk_20[k]
                   - f_14 * pc_z[k] * isl1_27[k];

        t_118[k] = f_20 * isk_100[k]
                   + f_3 * pc_x[k] * ksk_100[k];

        t_119[k] = f_20 * isk_101[k]
                   + f_3 * pc_x[k] * ksk_101[k];

        t_120[k] = f_20 * isk_102[k]
                   + f_3 * pc_x[k] * ksk_102[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, isk_103, isk_104, \
                         isk_105, isk_107, ksk_99, ksk_103, ksk_104, ksk_105, \
                         ksk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_20 * isk_103[k]
                   + f_3 * pc_x[k] * ksk_103[k];

        t_122[k] = f_20 * isk_104[k]
                   + f_3 * pc_x[k] * ksk_104[k];

        t_123[k] = f_20 * isk_105[k]
                   + f_3 * pc_x[k] * ksk_105[k];

        t_124[k] = f_3 * pc_y[k] * ksk_99[k];

        t_125[k] = f_20 * isk_107[k]
                   + f_3 * pc_x[k] * ksk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_z, pc_y, pc_z, isl0_36, isl1_36, ksi0_78, \
                         ksi0_79, ksi1_78, ksi1_79, ksk_101, ksk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * isl0_36[k]
                   - f_14 * pc_z[k] * isl1_36[k];

        t_127[k] = f_21 * ksi0_78[k]
                   - f_22 * ksi1_78[k]
                   + f_3 * pc_y[k] * ksk_101[k];

        t_128[k] = f_12 * ksi0_79[k]
                   - f_13 * ksi1_79[k]
                   + f_3 * pc_y[k] * ksk_102[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
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

    const auto *isl0_48 = buffer.data(isl0 + 48);
    const auto *isl0_51 = buffer.data(isl0 + 51);
    const auto *isl0_55 = buffer.data(isl0 + 55);
    const auto *isl0_60 = buffer.data(isl0 + 60);
    const auto *isl0_66 = buffer.data(isl0 + 66);
    const auto *isl0_81 = buffer.data(isl0 + 81);
    const auto *isl0_90 = buffer.data(isl0 + 90);
    const auto *isl0_95 = buffer.data(isl0 + 95);
    const auto *isl0_99 = buffer.data(isl0 + 99);
    const auto *isl0_102 = buffer.data(isl0 + 102);
    const auto *isl0_104 = buffer.data(isl0 + 104);
    const auto *isl0_107 = buffer.data(isl0 + 107);
    const auto *isl0_108 = buffer.data(isl0 + 108);
    const auto *isl0_110 = buffer.data(isl0 + 110);
    const auto *isl0_113 = buffer.data(isl0 + 113);
    const auto *isl0_114 = buffer.data(isl0 + 114);
    const auto *isl0_115 = buffer.data(isl0 + 115);
    const auto *isl0_117 = buffer.data(isl0 + 117);
    const auto *isl0_134 = buffer.data(isl0 + 134);

    const auto *isk_35 = buffer.data(isk + 35);
    const auto *isk_36 = buffer.data(isk + 36);
    const auto *isk_39 = buffer.data(isk + 39);
    const auto *isk_41 = buffer.data(isk + 41);
    const auto *isk_42 = buffer.data(isk + 42);
    const auto *isk_45 = buffer.data(isk + 45);
    const auto *isk_46 = buffer.data(isk + 46);
    const auto *isk_50 = buffer.data(isk + 50);
    const auto *isk_51 = buffer.data(isk + 51);
    const auto *isk_56 = buffer.data(isk + 56);
    const auto *isk_64 = buffer.data(isk + 64);
    const auto *isk_71 = buffer.data(isk + 71);
    const auto *isk_72 = buffer.data(isk + 72);
    const auto *isk_74 = buffer.data(isk + 74);
    const auto *isk_77 = buffer.data(isk + 77);
    const auto *isk_80 = buffer.data(isk + 80);
    const auto *isk_81 = buffer.data(isk + 81);
    const auto *isk_84 = buffer.data(isk + 84);
    const auto *isk_85 = buffer.data(isk + 85);
    const auto *isk_86 = buffer.data(isk + 86);
    const auto *isk_89 = buffer.data(isk + 89);
    const auto *isk_90 = buffer.data(isk + 90);
    const auto *isk_91 = buffer.data(isk + 91);
    const auto *isk_92 = buffer.data(isk + 92);
    const auto *isk_102 = buffer.data(isk + 102);
    const auto *isk_103 = buffer.data(isk + 103);
    const auto *isk_104 = buffer.data(isk + 104);
    const auto *isk_105 = buffer.data(isk + 105);
    const auto *isk_106 = buffer.data(isk + 106);
    const auto *isk_107 = buffer.data(isk + 107);
    const auto *isk_108 = buffer.data(isk + 108);
    const auto *isk_111 = buffer.data(isk + 111);
    const auto *isk_114 = buffer.data(isk + 114);
    const auto *isk_118 = buffer.data(isk + 118);
    const auto *isk_123 = buffer.data(isk + 123);
    const auto *isk_129 = buffer.data(isk + 129);
    const auto *isk_136 = buffer.data(isk + 136);
    const auto *isk_138 = buffer.data(isk + 138);
    const auto *isk_139 = buffer.data(isk + 139);
    const auto *isk_140 = buffer.data(isk + 140);
    const auto *isk_141 = buffer.data(isk + 141);
    const auto *isk_142 = buffer.data(isk + 142);
    const auto *isk_143 = buffer.data(isk + 143);
    const auto *isk_172 = buffer.data(isk + 172);
    const auto *isk_173 = buffer.data(isk + 173);
    const auto *isk_174 = buffer.data(isk + 174);
    const auto *isk_175 = buffer.data(isk + 175);
    const auto *isk_176 = buffer.data(isk + 176);
    const auto *isk_177 = buffer.data(isk + 177);
    const auto *isk_178 = buffer.data(isk + 178);
    const auto *isk_179 = buffer.data(isk + 179);
    const auto *isk_180 = buffer.data(isk + 180);
    const auto *isk_185 = buffer.data(isk + 185);
    const auto *isk_189 = buffer.data(isk + 189);
    const auto *isk_194 = buffer.data(isk + 194);
    const auto *isk_200 = buffer.data(isk + 200);

    const auto *isl1_48 = buffer.data(isl1 + 48);
    const auto *isl1_51 = buffer.data(isl1 + 51);
    const auto *isl1_55 = buffer.data(isl1 + 55);
    const auto *isl1_60 = buffer.data(isl1 + 60);
    const auto *isl1_66 = buffer.data(isl1 + 66);
    const auto *isl1_81 = buffer.data(isl1 + 81);
    const auto *isl1_90 = buffer.data(isl1 + 90);
    const auto *isl1_95 = buffer.data(isl1 + 95);
    const auto *isl1_99 = buffer.data(isl1 + 99);
    const auto *isl1_102 = buffer.data(isl1 + 102);
    const auto *isl1_104 = buffer.data(isl1 + 104);
    const auto *isl1_107 = buffer.data(isl1 + 107);
    const auto *isl1_108 = buffer.data(isl1 + 108);
    const auto *isl1_110 = buffer.data(isl1 + 110);
    const auto *isl1_113 = buffer.data(isl1 + 113);
    const auto *isl1_114 = buffer.data(isl1 + 114);
    const auto *isl1_115 = buffer.data(isl1 + 115);
    const auto *isl1_117 = buffer.data(isl1 + 117);
    const auto *isl1_134 = buffer.data(isl1 + 134);

    const auto *ksi0_80 = buffer.data(ksi0 + 80);
    const auto *ksi0_81 = buffer.data(ksi0 + 81);
    const auto *ksi0_82 = buffer.data(ksi0 + 82);
    const auto *ksi0_83 = buffer.data(ksi0 + 83);
    const auto *ksi0_84 = buffer.data(ksi0 + 84);
    const auto *ksi0_86 = buffer.data(ksi0 + 86);
    const auto *ksi0_87 = buffer.data(ksi0 + 87);
    const auto *ksi0_89 = buffer.data(ksi0 + 89);
    const auto *ksi0_90 = buffer.data(ksi0 + 90);
    const auto *ksi0_91 = buffer.data(ksi0 + 91);
    const auto *ksi0_93 = buffer.data(ksi0 + 93);
    const auto *ksi0_94 = buffer.data(ksi0 + 94);
    const auto *ksi0_95 = buffer.data(ksi0 + 95);
    const auto *ksi0_96 = buffer.data(ksi0 + 96);
    const auto *ksi0_98 = buffer.data(ksi0 + 98);
    const auto *ksi0_99 = buffer.data(ksi0 + 99);
    const auto *ksi0_105 = buffer.data(ksi0 + 105);
    const auto *ksi0_106 = buffer.data(ksi0 + 106);
    const auto *ksi0_107 = buffer.data(ksi0 + 107);
    const auto *ksi0_108 = buffer.data(ksi0 + 108);
    const auto *ksi0_109 = buffer.data(ksi0 + 109);
    const auto *ksi0_111 = buffer.data(ksi0 + 111);
    const auto *ksi0_135 = buffer.data(ksi0 + 135);
    const auto *ksi0_136 = buffer.data(ksi0 + 136);
    const auto *ksi0_137 = buffer.data(ksi0 + 137);
    const auto *ksi0_138 = buffer.data(ksi0 + 138);
    const auto *ksi0_139 = buffer.data(ksi0 + 139);
    const auto *ksi0_140 = buffer.data(ksi0 + 140);
    const auto *ksi0_141 = buffer.data(ksi0 + 141);
    const auto *ksi0_142 = buffer.data(ksi0 + 142);
    const auto *ksi0_143 = buffer.data(ksi0 + 143);
    const auto *ksi0_144 = buffer.data(ksi0 + 144);
    const auto *ksi0_145 = buffer.data(ksi0 + 145);
    const auto *ksi0_146 = buffer.data(ksi0 + 146);
    const auto *ksi0_147 = buffer.data(ksi0 + 147);
    const auto *ksi0_148 = buffer.data(ksi0 + 148);
    const auto *ksi0_149 = buffer.data(ksi0 + 149);
    const auto *ksi0_154 = buffer.data(ksi0 + 154);
    const auto *ksi0_160 = buffer.data(ksi0 + 160);

    const auto *ksi1_80 = buffer.data(ksi1 + 80);
    const auto *ksi1_81 = buffer.data(ksi1 + 81);
    const auto *ksi1_82 = buffer.data(ksi1 + 82);
    const auto *ksi1_83 = buffer.data(ksi1 + 83);
    const auto *ksi1_84 = buffer.data(ksi1 + 84);
    const auto *ksi1_86 = buffer.data(ksi1 + 86);
    const auto *ksi1_87 = buffer.data(ksi1 + 87);
    const auto *ksi1_89 = buffer.data(ksi1 + 89);
    const auto *ksi1_90 = buffer.data(ksi1 + 90);
    const auto *ksi1_91 = buffer.data(ksi1 + 91);
    const auto *ksi1_93 = buffer.data(ksi1 + 93);
    const auto *ksi1_94 = buffer.data(ksi1 + 94);
    const auto *ksi1_95 = buffer.data(ksi1 + 95);
    const auto *ksi1_96 = buffer.data(ksi1 + 96);
    const auto *ksi1_98 = buffer.data(ksi1 + 98);
    const auto *ksi1_99 = buffer.data(ksi1 + 99);
    const auto *ksi1_105 = buffer.data(ksi1 + 105);
    const auto *ksi1_106 = buffer.data(ksi1 + 106);
    const auto *ksi1_107 = buffer.data(ksi1 + 107);
    const auto *ksi1_108 = buffer.data(ksi1 + 108);
    const auto *ksi1_109 = buffer.data(ksi1 + 109);
    const auto *ksi1_111 = buffer.data(ksi1 + 111);
    const auto *ksi1_135 = buffer.data(ksi1 + 135);
    const auto *ksi1_136 = buffer.data(ksi1 + 136);
    const auto *ksi1_137 = buffer.data(ksi1 + 137);
    const auto *ksi1_138 = buffer.data(ksi1 + 138);
    const auto *ksi1_139 = buffer.data(ksi1 + 139);
    const auto *ksi1_140 = buffer.data(ksi1 + 140);
    const auto *ksi1_141 = buffer.data(ksi1 + 141);
    const auto *ksi1_142 = buffer.data(ksi1 + 142);
    const auto *ksi1_143 = buffer.data(ksi1 + 143);
    const auto *ksi1_144 = buffer.data(ksi1 + 144);
    const auto *ksi1_145 = buffer.data(ksi1 + 145);
    const auto *ksi1_146 = buffer.data(ksi1 + 146);
    const auto *ksi1_147 = buffer.data(ksi1 + 147);
    const auto *ksi1_148 = buffer.data(ksi1 + 148);
    const auto *ksi1_149 = buffer.data(ksi1 + 149);
    const auto *ksi1_154 = buffer.data(ksi1 + 154);
    const auto *ksi1_160 = buffer.data(ksi1 + 160);

    const auto *ksk_103 = buffer.data(ksk + 103);
    const auto *ksk_104 = buffer.data(ksk + 104);
    const auto *ksk_105 = buffer.data(ksk + 105);
    const auto *ksk_106 = buffer.data(ksk + 106);
    const auto *ksk_107 = buffer.data(ksk + 107);
    const auto *ksk_108 = buffer.data(ksk + 108);
    const auto *ksk_109 = buffer.data(ksk + 109);
    const auto *ksk_110 = buffer.data(ksk + 110);
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
    const auto *ksk_129 = buffer.data(ksk + 129);
    const auto *ksk_136 = buffer.data(ksk + 136);
    const auto *ksk_137 = buffer.data(ksk + 137);
    const auto *ksk_138 = buffer.data(ksk + 138);
    const auto *ksk_139 = buffer.data(ksk + 139);
    const auto *ksk_140 = buffer.data(ksk + 140);
    const auto *ksk_141 = buffer.data(ksk + 141);
    const auto *ksk_142 = buffer.data(ksk + 142);
    const auto *ksk_143 = buffer.data(ksk + 143);
    const auto *ksk_144 = buffer.data(ksk + 144);
    const auto *ksk_146 = buffer.data(ksk + 146);
    const auto *ksk_147 = buffer.data(ksk + 147);
    const auto *ksk_149 = buffer.data(ksk + 149);
    const auto *ksk_150 = buffer.data(ksk + 150);
    const auto *ksk_153 = buffer.data(ksk + 153);
    const auto *ksk_154 = buffer.data(ksk + 154);
    const auto *ksk_158 = buffer.data(ksk + 158);
    const auto *ksk_159 = buffer.data(ksk + 159);
    const auto *ksk_164 = buffer.data(ksk + 164);
    const auto *ksk_172 = buffer.data(ksk + 172);
    const auto *ksk_173 = buffer.data(ksk + 173);
    const auto *ksk_174 = buffer.data(ksk + 174);
    const auto *ksk_175 = buffer.data(ksk + 175);
    const auto *ksk_176 = buffer.data(ksk + 176);
    const auto *ksk_177 = buffer.data(ksk + 177);
    const auto *ksk_178 = buffer.data(ksk + 178);
    const auto *ksk_179 = buffer.data(ksk + 179);
    const auto *ksk_180 = buffer.data(ksk + 180);
    const auto *ksk_181 = buffer.data(ksk + 181);
    const auto *ksk_182 = buffer.data(ksk + 182);
    const auto *ksk_183 = buffer.data(ksk + 183);
    const auto *ksk_184 = buffer.data(ksk + 184);
    const auto *ksk_185 = buffer.data(ksk + 185);
    const auto *ksk_186 = buffer.data(ksk + 186);
    const auto *ksk_187 = buffer.data(ksk + 187);
    const auto *ksk_188 = buffer.data(ksk + 188);
    const auto *ksk_189 = buffer.data(ksk + 189);
    const auto *ksk_190 = buffer.data(ksk + 190);
    const auto *ksk_191 = buffer.data(ksk + 191);
    const auto *ksk_192 = buffer.data(ksk + 192);
    const auto *ksk_193 = buffer.data(ksk + 193);
    const auto *ksk_194 = buffer.data(ksk + 194);
    const auto *ksk_200 = buffer.data(ksk + 200);

#pragma omp simd aligned(t_129, t_130, t_131, pc_y, ksi0_80, ksi0_81, ksi0_82, ksi1_80, \
                         ksi1_81, ksi1_82, ksk_103, ksk_104, ksk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * ksi0_80[k]
                   - f_11 * ksi1_80[k]
                   + f_3 * pc_y[k] * ksk_103[k];

        t_130[k] = f_8 * ksi0_81[k]
                   - f_9 * ksi1_81[k]
                   + f_3 * pc_y[k] * ksk_104[k];

        t_131[k] = f_6 * ksi0_82[k]
                   - f_7 * ksi1_82[k]
                   + f_3 * pc_y[k] * ksk_105[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, pc_y, pc_z, isk_35, isk_108, \
                         ksi0_83, ksi0_84, ksi1_83, ksi1_84, ksk_106, ksk_107, \
                         ksk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * ksi0_83[k]
                   - f_5 * ksi1_83[k]
                   + f_3 * pc_y[k] * ksk_106[k];

        t_133[k] = f_3 * pc_y[k] * ksk_107[k];

        t_134[k] = f_15 * isk_35[k]
                   + f_1 * ksi0_83[k]
                   - f_2 * ksi1_83[k]
                   + f_3 * pc_z[k] * ksk_107[k];

        t_135[k] = f_19 * isk_108[k]
                   + f_1 * ksi0_84[k]
                   - f_2 * ksi1_84[k]
                   + f_3 * pc_x[k] * ksk_108[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, isk_36, isk_111, \
                         ksi0_87, ksi1_87, ksk_108, ksk_109, ksk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_16 * isk_36[k]
                   + f_3 * pc_y[k] * ksk_108[k];

        t_137[k] = f_3 * pc_z[k] * ksk_108[k];

        t_138[k] = f_19 * isk_111[k]
                   + f_12 * ksi0_87[k]
                   - f_13 * ksi1_87[k]
                   + f_3 * pc_x[k] * ksk_111[k];

        t_139[k] = f_3 * pc_z[k] * ksk_109[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pc_x, pc_z, isk_114, ksi0_84, ksi0_90, ksi1_84, \
                         ksi1_90, ksk_110, ksk_111, ksk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * ksi0_84[k]
                   - f_5 * ksi1_84[k]
                   + f_3 * pc_z[k] * ksk_110[k];

        t_141[k] = f_19 * isk_114[k]
                   + f_10 * ksi0_90[k]
                   - f_11 * ksi1_90[k]
                   + f_3 * pc_x[k] * ksk_114[k];

        t_142[k] = f_3 * pc_z[k] * ksk_111[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, isk_41, isk_118, \
                         ksi0_86, ksi0_94, ksi1_86, ksi1_94, ksk_113, ksk_114, \
                         ksk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * isk_41[k]
                   + f_3 * pc_y[k] * ksk_113[k];

        t_144[k] = f_6 * ksi0_86[k]
                   - f_7 * ksi1_86[k]
                   + f_3 * pc_z[k] * ksk_113[k];

        t_145[k] = f_19 * isk_118[k]
                   + f_8 * ksi0_94[k]
                   - f_9 * ksi1_94[k]
                   + f_3 * pc_x[k] * ksk_118[k];

        t_146[k] = f_3 * pc_z[k] * ksk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_y, pc_z, isk_45, ksi0_87, ksi0_89, ksi1_87, \
                         ksi1_89, ksk_115, ksk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * ksi0_87[k]
                   - f_5 * ksi1_87[k]
                   + f_3 * pc_z[k] * ksk_115[k];

        t_148[k] = f_16 * isk_45[k]
                   + f_3 * pc_y[k] * ksk_117[k];

        t_149[k] = f_8 * ksi0_89[k]
                   - f_9 * ksi1_89[k]
                   + f_3 * pc_z[k] * ksk_117[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, isk_123, ksi0_90, ksi0_99, ksi1_90, \
                         ksi1_99, ksk_118, ksk_119, ksk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_19 * isk_123[k]
                   + f_6 * ksi0_99[k]
                   - f_7 * ksi1_99[k]
                   + f_3 * pc_x[k] * ksk_123[k];

        t_151[k] = f_3 * pc_z[k] * ksk_118[k];

        t_152[k] = f_4 * ksi0_90[k]
                   - f_5 * ksi1_90[k]
                   + f_3 * pc_z[k] * ksk_119[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_y, pc_z, isk_50, ksi0_91, ksi0_93, ksi1_91, \
                         ksi1_93, ksk_120, ksk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * ksi0_91[k]
                   - f_7 * ksi1_91[k]
                   + f_3 * pc_z[k] * ksk_120[k];

        t_154[k] = f_16 * isk_50[k]
                   + f_3 * pc_y[k] * ksk_122[k];

        t_155[k] = f_10 * ksi0_93[k]
                   - f_11 * ksi1_93[k]
                   + f_3 * pc_z[k] * ksk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, isk_129, ksi0_94, ksi0_105, ksi1_94, \
                         ksi1_105, ksk_123, ksk_124, ksk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_19 * isk_129[k]
                   + f_4 * ksi0_105[k]
                   - f_5 * ksi1_105[k]
                   + f_3 * pc_x[k] * ksk_129[k];

        t_157[k] = f_3 * pc_z[k] * ksk_123[k];

        t_158[k] = f_4 * ksi0_94[k]
                   - f_5 * ksi1_94[k]
                   + f_3 * pc_z[k] * ksk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pc_y, pc_z, isk_56, ksi0_95, ksi0_96, \
                         ksi0_98, ksi1_95, ksi1_96, ksi1_98, ksk_125, ksk_126, \
                         ksk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_6 * ksi0_95[k]
                   - f_7 * ksi1_95[k]
                   + f_3 * pc_z[k] * ksk_125[k];

        t_160[k] = f_8 * ksi0_96[k]
                   - f_9 * ksi1_96[k]
                   + f_3 * pc_z[k] * ksk_126[k];

        t_161[k] = f_16 * isk_56[k]
                   + f_3 * pc_y[k] * ksk_128[k];

        t_162[k] = f_12 * ksi0_98[k]
                   - f_13 * ksi1_98[k]
                   + f_3 * pc_z[k] * ksk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pc_x, pc_z, isk_136, isk_138, \
                         isk_139, isk_140, ksk_129, ksk_136, ksk_138, ksk_139, \
                         ksk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_19 * isk_136[k]
                   + f_3 * pc_x[k] * ksk_136[k];

        t_164[k] = f_3 * pc_z[k] * ksk_129[k];

        t_165[k] = f_19 * isk_138[k]
                   + f_3 * pc_x[k] * ksk_138[k];

        t_166[k] = f_19 * isk_139[k]
                   + f_3 * pc_x[k] * ksk_139[k];

        t_167[k] = f_19 * isk_140[k]
                   + f_3 * pc_x[k] * ksk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, isk_64, isk_141, isk_142, \
                         isk_143, ksi0_105, ksi1_105, ksk_136, ksk_141, ksk_142, \
                         ksk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * isk_141[k]
                   + f_3 * pc_x[k] * ksk_141[k];

        t_169[k] = f_19 * isk_142[k]
                   + f_3 * pc_x[k] * ksk_142[k];

        t_170[k] = f_19 * isk_143[k]
                   + f_3 * pc_x[k] * ksk_143[k];

        t_171[k] = f_16 * isk_64[k]
                   + f_1 * ksi0_105[k]
                   - f_2 * ksi1_105[k]
                   + f_3 * pc_y[k] * ksk_136[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_z, ksi0_105, ksi0_106, ksi0_107, \
                         ksi1_105, ksi1_106, ksi1_107, ksk_136, ksk_137, ksk_138, \
                         ksk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * ksk_136[k];

        t_173[k] = f_4 * ksi0_105[k]
                   - f_5 * ksi1_105[k]
                   + f_3 * pc_z[k] * ksk_137[k];

        t_174[k] = f_6 * ksi0_106[k]
                   - f_7 * ksi1_106[k]
                   + f_3 * pc_z[k] * ksk_138[k];

        t_175[k] = f_8 * ksi0_107[k]
                   - f_9 * ksi1_107[k]
                   + f_3 * pc_z[k] * ksk_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_y, pc_z, isk_71, ksi0_108, ksi0_109, \
                         ksi0_111, ksi1_108, ksi1_109, ksi1_111, ksk_140, ksk_141, \
                         ksk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * ksi0_108[k]
                   - f_11 * ksi1_108[k]
                   + f_3 * pc_z[k] * ksk_140[k];

        t_177[k] = f_12 * ksi0_109[k]
                   - f_13 * ksi1_109[k]
                   + f_3 * pc_z[k] * ksk_141[k];

        t_178[k] = f_16 * isk_71[k]
                   + f_3 * pc_y[k] * ksk_143[k];

        t_179[k] = f_1 * ksi0_111[k]
                   - f_2 * ksi1_111[k]
                   + f_3 * pc_z[k] * ksk_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, isl0_48, isl0_90, \
                         isk_36, isk_72, isl1_48, isl1_90, ksk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * isl0_90[k]
                   - f_14 * pc_y[k] * isl1_90[k];

        t_181[k] = f_15 * isk_72[k]
                   + f_3 * pc_y[k] * ksk_144[k];

        t_182[k] = f_15 * isk_36[k]
                   + f_3 * pc_z[k] * ksk_144[k];

        t_183[k] = pa_z[k] * isl0_48[k]
                   - f_14 * pc_z[k] * isl1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pa_z, pc_y, pc_z, isl0_51, isl0_95, \
                         isk_39, isk_74, isl1_51, isl1_95, ksk_146, \
                         ksk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_15 * isk_74[k]
                   + f_3 * pc_y[k] * ksk_146[k];

        t_185[k] = pa_y[k] * isl0_95[k]
                   - f_14 * pc_y[k] * isl1_95[k];

        t_186[k] = pa_z[k] * isl0_51[k]
                   - f_14 * pc_z[k] * isl1_51[k];

        t_187[k] = f_15 * isk_39[k]
                   + f_3 * pc_z[k] * ksk_147[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pc_y, pc_z, isl0_55, isl0_99, \
                         isk_42, isk_77, isl1_55, isl1_99, ksk_149, \
                         ksk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_15 * isk_77[k]
                   + f_3 * pc_y[k] * ksk_149[k];

        t_189[k] = pa_y[k] * isl0_99[k]
                   - f_14 * pc_y[k] * isl1_99[k];

        t_190[k] = pa_z[k] * isl0_55[k]
                   - f_14 * pc_z[k] * isl1_55[k];

        t_191[k] = f_15 * isk_42[k]
                   + f_3 * pc_z[k] * ksk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pc_y, isl0_102, isl0_104, isk_80, isk_81, \
                         isl1_102, isl1_104, ksk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_y[k] * isl0_102[k]
                   + f_16 * isk_80[k]
                   - f_14 * pc_y[k] * isl1_102[k];

        t_193[k] = f_15 * isk_81[k]
                   + f_3 * pc_y[k] * ksk_153[k];

        t_194[k] = pa_y[k] * isl0_104[k]
                   - f_14 * pc_y[k] * isl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_y, pa_z, pc_y, pc_z, isl0_60, isl0_107, \
                         isk_46, isk_84, isl1_60, isl1_107, ksk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * isl0_60[k]
                   - f_14 * pc_z[k] * isl1_60[k];

        t_196[k] = f_15 * isk_46[k]
                   + f_3 * pc_z[k] * ksk_154[k];

        t_197[k] = pa_y[k] * isl0_107[k]
                   + f_17 * isk_84[k]
                   - f_14 * pc_y[k] * isl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_y, pc_y, isl0_108, isl0_110, isk_85, isk_86, \
                         isl1_108, isl1_110, ksk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_y[k] * isl0_108[k]
                   + f_16 * isk_85[k]
                   - f_14 * pc_y[k] * isl1_108[k];

        t_199[k] = f_15 * isk_86[k]
                   + f_3 * pc_y[k] * ksk_158[k];

        t_200[k] = pa_y[k] * isl0_110[k]
                   - f_14 * pc_y[k] * isl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_y, pa_z, pc_y, pc_z, isl0_66, isl0_113, \
                         isk_51, isk_89, isl1_66, isl1_113, ksk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_z[k] * isl0_66[k]
                   - f_14 * pc_z[k] * isl1_66[k];

        t_202[k] = f_15 * isk_51[k]
                   + f_3 * pc_z[k] * ksk_159[k];

        t_203[k] = pa_y[k] * isl0_113[k]
                   + f_18 * isk_89[k]
                   - f_14 * pc_y[k] * isl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_y, isl0_114, isl0_115, isl0_117, \
                         isk_90, isk_91, isk_92, isl1_114, isl1_115, isl1_117, \
                         ksk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * isl0_114[k]
                   + f_17 * isk_90[k]
                   - f_14 * pc_y[k] * isl1_114[k];

        t_205[k] = pa_y[k] * isl0_115[k]
                   + f_16 * isk_91[k]
                   - f_14 * pc_y[k] * isl1_115[k];

        t_206[k] = f_15 * isk_92[k]
                   + f_3 * pc_y[k] * ksk_164[k];

        t_207[k] = pa_y[k] * isl0_117[k]
                   - f_14 * pc_y[k] * isl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, isk_172, isk_173, isk_174, \
                         isk_175, isk_176, ksk_172, ksk_173, ksk_174, ksk_175, \
                         ksk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_19 * isk_172[k]
                   + f_3 * pc_x[k] * ksk_172[k];

        t_209[k] = f_19 * isk_173[k]
                   + f_3 * pc_x[k] * ksk_173[k];

        t_210[k] = f_19 * isk_174[k]
                   + f_3 * pc_x[k] * ksk_174[k];

        t_211[k] = f_19 * isk_175[k]
                   + f_3 * pc_x[k] * ksk_175[k];

        t_212[k] = f_19 * isk_176[k]
                   + f_3 * pc_x[k] * ksk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_z, pc_x, pc_z, isl0_81, isk_177, \
                         isk_178, isk_179, isl1_81, ksk_177, ksk_178, \
                         ksk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_19 * isk_177[k]
                   + f_3 * pc_x[k] * ksk_177[k];

        t_214[k] = f_19 * isk_178[k]
                   + f_3 * pc_x[k] * ksk_178[k];

        t_215[k] = f_19 * isk_179[k]
                   + f_3 * pc_x[k] * ksk_179[k];

        t_216[k] = pa_z[k] * isl0_81[k]
                   - f_14 * pc_z[k] * isl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, isk_64, isk_102, isk_103, ksi0_135, \
                         ksi0_136, ksi1_135, ksi1_136, ksk_172, ksk_174, \
                         ksk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * isk_64[k]
                   + f_3 * pc_z[k] * ksk_172[k];

        t_218[k] = f_15 * isk_102[k]
                   + f_12 * ksi0_135[k]
                   - f_13 * ksi1_135[k]
                   + f_3 * pc_y[k] * ksk_174[k];

        t_219[k] = f_15 * isk_103[k]
                   + f_10 * ksi0_136[k]
                   - f_11 * ksi1_136[k]
                   + f_3 * pc_y[k] * ksk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, isk_104, isk_105, isk_106, ksi0_137, \
                         ksi0_138, ksi0_139, ksi1_137, ksi1_138, ksi1_139, ksk_176, ksk_177, \
                         ksk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * isk_104[k]
                   + f_8 * ksi0_137[k]
                   - f_9 * ksi1_137[k]
                   + f_3 * pc_y[k] * ksk_176[k];

        t_221[k] = f_15 * isk_105[k]
                   + f_6 * ksi0_138[k]
                   - f_7 * ksi1_138[k]
                   + f_3 * pc_y[k] * ksk_177[k];

        t_222[k] = f_15 * isk_106[k]
                   + f_4 * ksi0_139[k]
                   - f_5 * ksi1_139[k]
                   + f_3 * pc_y[k] * ksk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_y, pc_x, pc_y, isl0_134, isk_107, \
                         isk_180, isl1_134, ksi0_140, ksi1_140, ksk_179, \
                         ksk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * isk_107[k]
                   + f_3 * pc_y[k] * ksk_179[k];

        t_224[k] = pa_y[k] * isl0_134[k]
                   - f_14 * pc_y[k] * isl1_134[k];

        t_225[k] = f_19 * isk_180[k]
                   + f_1 * ksi0_140[k]
                   - f_2 * ksi1_140[k]
                   + f_3 * pc_x[k] * ksk_180[k];

        t_226[k] = f_3 * pc_y[k] * ksk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_y, pc_z, isk_72, ksi0_140, ksi1_140, ksk_180, \
                         ksk_181, ksk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * isk_72[k]
                   + f_3 * pc_z[k] * ksk_180[k];

        t_228[k] = f_4 * ksi0_140[k]
                   - f_5 * ksi1_140[k]
                   + f_3 * pc_y[k] * ksk_181[k];

        t_229[k] = f_3 * pc_y[k] * ksk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, isk_185, ksi0_141, ksi0_142, \
                         ksi0_145, ksi1_141, ksi1_142, ksi1_145, ksk_183, ksk_184, \
                         ksk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_19 * isk_185[k]
                   + f_12 * ksi0_145[k]
                   - f_13 * ksi1_145[k]
                   + f_3 * pc_x[k] * ksk_185[k];

        t_231[k] = f_6 * ksi0_141[k]
                   - f_7 * ksi1_141[k]
                   + f_3 * pc_y[k] * ksk_183[k];

        t_232[k] = f_4 * ksi0_142[k]
                   - f_5 * ksi1_142[k]
                   + f_3 * pc_y[k] * ksk_184[k];

        t_233[k] = f_3 * pc_y[k] * ksk_185[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, isk_189, ksi0_143, ksi0_144, \
                         ksi0_149, ksi1_143, ksi1_144, ksi1_149, ksk_186, ksk_187, \
                         ksk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_19 * isk_189[k]
                   + f_10 * ksi0_149[k]
                   - f_11 * ksi1_149[k]
                   + f_3 * pc_x[k] * ksk_189[k];

        t_235[k] = f_8 * ksi0_143[k]
                   - f_9 * ksi1_143[k]
                   + f_3 * pc_y[k] * ksk_186[k];

        t_236[k] = f_6 * ksi0_144[k]
                   - f_7 * ksi1_144[k]
                   + f_3 * pc_y[k] * ksk_187[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_x, pc_y, isk_194, ksi0_145, ksi0_154, \
                         ksi1_145, ksi1_154, ksk_188, ksk_189, \
                         ksk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * ksi0_145[k]
                   - f_5 * ksi1_145[k]
                   + f_3 * pc_y[k] * ksk_188[k];

        t_238[k] = f_3 * pc_y[k] * ksk_189[k];

        t_239[k] = f_19 * isk_194[k]
                   + f_8 * ksi0_154[k]
                   - f_9 * ksi1_154[k]
                   + f_3 * pc_x[k] * ksk_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pc_y, ksi0_146, ksi0_147, ksi0_148, ksi1_146, \
                         ksi1_147, ksi1_148, ksk_190, ksk_191, \
                         ksk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * ksi0_146[k]
                   - f_11 * ksi1_146[k]
                   + f_3 * pc_y[k] * ksk_190[k];

        t_241[k] = f_8 * ksi0_147[k]
                   - f_9 * ksi1_147[k]
                   + f_3 * pc_y[k] * ksk_191[k];

        t_242[k] = f_6 * ksi0_148[k]
                   - f_7 * ksi1_148[k]
                   + f_3 * pc_y[k] * ksk_192[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pc_x, pc_y, isk_200, ksi0_149, ksi0_160, \
                         ksi1_149, ksi1_160, ksk_193, ksk_194, \
                         ksk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_4 * ksi0_149[k]
                   - f_5 * ksi1_149[k]
                   + f_3 * pc_y[k] * ksk_193[k];

        t_244[k] = f_3 * pc_y[k] * ksk_194[k];

        t_245[k] = f_19 * isk_200[k]
                   + f_6 * ksi0_160[k]
                   - f_7 * ksi1_160[k]
                   + f_3 * pc_x[k] * ksk_200[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
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
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);

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

    const auto *isl0_135 = buffer.data(isl0 + 135);
    const auto *isl0_138 = buffer.data(isl0 + 138);
    const auto *isl0_141 = buffer.data(isl0 + 141);
    const auto *isl0_145 = buffer.data(isl0 + 145);
    const auto *isl0_147 = buffer.data(isl0 + 147);
    const auto *isl0_150 = buffer.data(isl0 + 150);
    const auto *isl0_152 = buffer.data(isl0 + 152);
    const auto *isl0_153 = buffer.data(isl0 + 153);
    const auto *isl0_156 = buffer.data(isl0 + 156);
    const auto *isl0_158 = buffer.data(isl0 + 158);
    const auto *isl0_159 = buffer.data(isl0 + 159);
    const auto *isl0_160 = buffer.data(isl0 + 160);
    const auto *isl0_171 = buffer.data(isl0 + 171);
    const auto *isl0_225 = buffer.data(isl0 + 225);

    const auto *isk_107 = buffer.data(isk + 107);
    const auto *isk_108 = buffer.data(isk + 108);
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
    const auto *isk_136 = buffer.data(isk + 136);
    const auto *isk_143 = buffer.data(isk + 143);
    const auto *isk_144 = buffer.data(isk + 144);
    const auto *isk_146 = buffer.data(isk + 146);
    const auto *isk_149 = buffer.data(isk + 149);
    const auto *isk_153 = buffer.data(isk + 153);
    const auto *isk_158 = buffer.data(isk + 158);
    const auto *isk_164 = buffer.data(isk + 164);
    const auto *isk_174 = buffer.data(isk + 174);
    const auto *isk_175 = buffer.data(isk + 175);
    const auto *isk_176 = buffer.data(isk + 176);
    const auto *isk_177 = buffer.data(isk + 177);
    const auto *isk_178 = buffer.data(isk + 178);
    const auto *isk_179 = buffer.data(isk + 179);
    const auto *isk_180 = buffer.data(isk + 180);
    const auto *isk_207 = buffer.data(isk + 207);
    const auto *isk_208 = buffer.data(isk + 208);
    const auto *isk_209 = buffer.data(isk + 209);
    const auto *isk_210 = buffer.data(isk + 210);
    const auto *isk_211 = buffer.data(isk + 211);
    const auto *isk_212 = buffer.data(isk + 212);
    const auto *isk_213 = buffer.data(isk + 213);
    const auto *isk_215 = buffer.data(isk + 215);
    const auto *isk_216 = buffer.data(isk + 216);
    const auto *isk_219 = buffer.data(isk + 219);
    const auto *isk_222 = buffer.data(isk + 222);
    const auto *isk_226 = buffer.data(isk + 226);
    const auto *isk_231 = buffer.data(isk + 231);
    const auto *isk_237 = buffer.data(isk + 237);
    const auto *isk_244 = buffer.data(isk + 244);
    const auto *isk_246 = buffer.data(isk + 246);
    const auto *isk_247 = buffer.data(isk + 247);
    const auto *isk_248 = buffer.data(isk + 248);
    const auto *isk_249 = buffer.data(isk + 249);
    const auto *isk_250 = buffer.data(isk + 250);
    const auto *isk_251 = buffer.data(isk + 251);
    const auto *isk_257 = buffer.data(isk + 257);
    const auto *isk_261 = buffer.data(isk + 261);
    const auto *isk_266 = buffer.data(isk + 266);
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

    const auto *isl1_135 = buffer.data(isl1 + 135);
    const auto *isl1_138 = buffer.data(isl1 + 138);
    const auto *isl1_141 = buffer.data(isl1 + 141);
    const auto *isl1_145 = buffer.data(isl1 + 145);
    const auto *isl1_147 = buffer.data(isl1 + 147);
    const auto *isl1_150 = buffer.data(isl1 + 150);
    const auto *isl1_152 = buffer.data(isl1 + 152);
    const auto *isl1_153 = buffer.data(isl1 + 153);
    const auto *isl1_156 = buffer.data(isl1 + 156);
    const auto *isl1_158 = buffer.data(isl1 + 158);
    const auto *isl1_159 = buffer.data(isl1 + 159);
    const auto *isl1_160 = buffer.data(isl1 + 160);
    const auto *isl1_171 = buffer.data(isl1 + 171);
    const auto *isl1_225 = buffer.data(isl1 + 225);

    const auto *ksi0_150 = buffer.data(ksi0 + 150);
    const auto *ksi0_151 = buffer.data(ksi0 + 151);
    const auto *ksi0_152 = buffer.data(ksi0 + 152);
    const auto *ksi0_153 = buffer.data(ksi0 + 153);
    const auto *ksi0_154 = buffer.data(ksi0 + 154);
    const auto *ksi0_161 = buffer.data(ksi0 + 161);
    const auto *ksi0_162 = buffer.data(ksi0 + 162);
    const auto *ksi0_163 = buffer.data(ksi0 + 163);
    const auto *ksi0_164 = buffer.data(ksi0 + 164);
    const auto *ksi0_165 = buffer.data(ksi0 + 165);
    const auto *ksi0_166 = buffer.data(ksi0 + 166);
    const auto *ksi0_167 = buffer.data(ksi0 + 167);
    const auto *ksi0_168 = buffer.data(ksi0 + 168);
    const auto *ksi0_170 = buffer.data(ksi0 + 170);
    const auto *ksi0_171 = buffer.data(ksi0 + 171);
    const auto *ksi0_173 = buffer.data(ksi0 + 173);
    const auto *ksi0_174 = buffer.data(ksi0 + 174);
    const auto *ksi0_175 = buffer.data(ksi0 + 175);
    const auto *ksi0_177 = buffer.data(ksi0 + 177);
    const auto *ksi0_178 = buffer.data(ksi0 + 178);
    const auto *ksi0_179 = buffer.data(ksi0 + 179);
    const auto *ksi0_180 = buffer.data(ksi0 + 180);
    const auto *ksi0_182 = buffer.data(ksi0 + 182);
    const auto *ksi0_183 = buffer.data(ksi0 + 183);
    const auto *ksi0_189 = buffer.data(ksi0 + 189);
    const auto *ksi0_190 = buffer.data(ksi0 + 190);
    const auto *ksi0_191 = buffer.data(ksi0 + 191);
    const auto *ksi0_192 = buffer.data(ksi0 + 192);
    const auto *ksi0_193 = buffer.data(ksi0 + 193);
    const auto *ksi0_195 = buffer.data(ksi0 + 195);
    const auto *ksi0_201 = buffer.data(ksi0 + 201);
    const auto *ksi0_205 = buffer.data(ksi0 + 205);
    const auto *ksi0_210 = buffer.data(ksi0 + 210);
    const auto *ksi0_216 = buffer.data(ksi0 + 216);
    const auto *ksi0_219 = buffer.data(ksi0 + 219);
    const auto *ksi0_220 = buffer.data(ksi0 + 220);
    const auto *ksi0_221 = buffer.data(ksi0 + 221);
    const auto *ksi0_222 = buffer.data(ksi0 + 222);
    const auto *ksi0_223 = buffer.data(ksi0 + 223);

    const auto *ksi1_150 = buffer.data(ksi1 + 150);
    const auto *ksi1_151 = buffer.data(ksi1 + 151);
    const auto *ksi1_152 = buffer.data(ksi1 + 152);
    const auto *ksi1_153 = buffer.data(ksi1 + 153);
    const auto *ksi1_154 = buffer.data(ksi1 + 154);
    const auto *ksi1_161 = buffer.data(ksi1 + 161);
    const auto *ksi1_162 = buffer.data(ksi1 + 162);
    const auto *ksi1_163 = buffer.data(ksi1 + 163);
    const auto *ksi1_164 = buffer.data(ksi1 + 164);
    const auto *ksi1_165 = buffer.data(ksi1 + 165);
    const auto *ksi1_166 = buffer.data(ksi1 + 166);
    const auto *ksi1_167 = buffer.data(ksi1 + 167);
    const auto *ksi1_168 = buffer.data(ksi1 + 168);
    const auto *ksi1_170 = buffer.data(ksi1 + 170);
    const auto *ksi1_171 = buffer.data(ksi1 + 171);
    const auto *ksi1_173 = buffer.data(ksi1 + 173);
    const auto *ksi1_174 = buffer.data(ksi1 + 174);
    const auto *ksi1_175 = buffer.data(ksi1 + 175);
    const auto *ksi1_177 = buffer.data(ksi1 + 177);
    const auto *ksi1_178 = buffer.data(ksi1 + 178);
    const auto *ksi1_179 = buffer.data(ksi1 + 179);
    const auto *ksi1_180 = buffer.data(ksi1 + 180);
    const auto *ksi1_182 = buffer.data(ksi1 + 182);
    const auto *ksi1_183 = buffer.data(ksi1 + 183);
    const auto *ksi1_189 = buffer.data(ksi1 + 189);
    const auto *ksi1_190 = buffer.data(ksi1 + 190);
    const auto *ksi1_191 = buffer.data(ksi1 + 191);
    const auto *ksi1_192 = buffer.data(ksi1 + 192);
    const auto *ksi1_193 = buffer.data(ksi1 + 193);
    const auto *ksi1_195 = buffer.data(ksi1 + 195);
    const auto *ksi1_201 = buffer.data(ksi1 + 201);
    const auto *ksi1_205 = buffer.data(ksi1 + 205);
    const auto *ksi1_210 = buffer.data(ksi1 + 210);
    const auto *ksi1_216 = buffer.data(ksi1 + 216);
    const auto *ksi1_219 = buffer.data(ksi1 + 219);
    const auto *ksi1_220 = buffer.data(ksi1 + 220);
    const auto *ksi1_221 = buffer.data(ksi1 + 221);
    const auto *ksi1_222 = buffer.data(ksi1 + 222);
    const auto *ksi1_223 = buffer.data(ksi1 + 223);

    const auto *ksk_195 = buffer.data(ksk + 195);
    const auto *ksk_196 = buffer.data(ksk + 196);
    const auto *ksk_197 = buffer.data(ksk + 197);
    const auto *ksk_198 = buffer.data(ksk + 198);
    const auto *ksk_199 = buffer.data(ksk + 199);
    const auto *ksk_200 = buffer.data(ksk + 200);
    const auto *ksk_207 = buffer.data(ksk + 207);
    const auto *ksk_208 = buffer.data(ksk + 208);
    const auto *ksk_209 = buffer.data(ksk + 209);
    const auto *ksk_210 = buffer.data(ksk + 210);
    const auto *ksk_211 = buffer.data(ksk + 211);
    const auto *ksk_212 = buffer.data(ksk + 212);
    const auto *ksk_213 = buffer.data(ksk + 213);
    const auto *ksk_214 = buffer.data(ksk + 214);
    const auto *ksk_215 = buffer.data(ksk + 215);
    const auto *ksk_216 = buffer.data(ksk + 216);
    const auto *ksk_217 = buffer.data(ksk + 217);
    const auto *ksk_218 = buffer.data(ksk + 218);
    const auto *ksk_219 = buffer.data(ksk + 219);
    const auto *ksk_221 = buffer.data(ksk + 221);
    const auto *ksk_222 = buffer.data(ksk + 222);
    const auto *ksk_223 = buffer.data(ksk + 223);
    const auto *ksk_225 = buffer.data(ksk + 225);
    const auto *ksk_226 = buffer.data(ksk + 226);
    const auto *ksk_227 = buffer.data(ksk + 227);
    const auto *ksk_228 = buffer.data(ksk + 228);
    const auto *ksk_230 = buffer.data(ksk + 230);
    const auto *ksk_231 = buffer.data(ksk + 231);
    const auto *ksk_232 = buffer.data(ksk + 232);
    const auto *ksk_233 = buffer.data(ksk + 233);
    const auto *ksk_234 = buffer.data(ksk + 234);
    const auto *ksk_236 = buffer.data(ksk + 236);
    const auto *ksk_237 = buffer.data(ksk + 237);
    const auto *ksk_244 = buffer.data(ksk + 244);
    const auto *ksk_245 = buffer.data(ksk + 245);
    const auto *ksk_246 = buffer.data(ksk + 246);
    const auto *ksk_247 = buffer.data(ksk + 247);
    const auto *ksk_248 = buffer.data(ksk + 248);
    const auto *ksk_249 = buffer.data(ksk + 249);
    const auto *ksk_250 = buffer.data(ksk + 250);
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
    const auto *ksk_279 = buffer.data(ksk + 279);
    const auto *ksk_280 = buffer.data(ksk + 280);
    const auto *ksk_281 = buffer.data(ksk + 281);
    const auto *ksk_282 = buffer.data(ksk + 282);
    const auto *ksk_283 = buffer.data(ksk + 283);
    const auto *ksk_284 = buffer.data(ksk + 284);
    const auto *ksk_285 = buffer.data(ksk + 285);
    const auto *ksk_286 = buffer.data(ksk + 286);
    const auto *ksk_287 = buffer.data(ksk + 287);
    const auto *ksk_288 = buffer.data(ksk + 288);

#pragma omp simd aligned(t_246, t_247, t_248, pc_y, ksi0_150, ksi0_151, ksi0_152, ksi1_150, \
                         ksi1_151, ksi1_152, ksk_195, ksk_196, \
                         ksk_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * ksi0_150[k]
                   - f_13 * ksi1_150[k]
                   + f_3 * pc_y[k] * ksk_195[k];

        t_247[k] = f_10 * ksi0_151[k]
                   - f_11 * ksi1_151[k]
                   + f_3 * pc_y[k] * ksk_196[k];

        t_248[k] = f_8 * ksi0_152[k]
                   - f_9 * ksi1_152[k]
                   + f_3 * pc_y[k] * ksk_197[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, ksi0_153, ksi0_154, ksi1_153, ksi1_154, \
                         ksk_198, ksk_199, ksk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * ksi0_153[k]
                   - f_7 * ksi1_153[k]
                   + f_3 * pc_y[k] * ksk_198[k];

        t_250[k] = f_4 * ksi0_154[k]
                   - f_5 * ksi1_154[k]
                   + f_3 * pc_y[k] * ksk_199[k];

        t_251[k] = f_3 * pc_y[k] * ksk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, isk_207, isk_208, isk_209, isk_210, \
                         ksi0_167, ksi1_167, ksk_207, ksk_208, ksk_209, \
                         ksk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_19 * isk_207[k]
                   + f_4 * ksi0_167[k]
                   - f_5 * ksi1_167[k]
                   + f_3 * pc_x[k] * ksk_207[k];

        t_253[k] = f_19 * isk_208[k]
                   + f_3 * pc_x[k] * ksk_208[k];

        t_254[k] = f_19 * isk_209[k]
                   + f_3 * pc_x[k] * ksk_209[k];

        t_255[k] = f_19 * isk_210[k]
                   + f_3 * pc_x[k] * ksk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pc_x, pc_y, isk_211, isk_212, \
                         isk_213, isk_215, ksk_207, ksk_211, ksk_212, ksk_213, \
                         ksk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_19 * isk_211[k]
                   + f_3 * pc_x[k] * ksk_211[k];

        t_257[k] = f_19 * isk_212[k]
                   + f_3 * pc_x[k] * ksk_212[k];

        t_258[k] = f_19 * isk_213[k]
                   + f_3 * pc_x[k] * ksk_213[k];

        t_259[k] = f_3 * pc_y[k] * ksk_207[k];

        t_260[k] = f_19 * isk_215[k]
                   + f_3 * pc_x[k] * ksk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_y, ksi0_161, ksi0_162, ksi0_163, ksi1_161, \
                         ksi1_162, ksi1_163, ksk_208, ksk_209, \
                         ksk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * ksi0_161[k]
                   - f_2 * ksi1_161[k]
                   + f_3 * pc_y[k] * ksk_208[k];

        t_262[k] = f_21 * ksi0_162[k]
                   - f_22 * ksi1_162[k]
                   + f_3 * pc_y[k] * ksk_209[k];

        t_263[k] = f_12 * ksi0_163[k]
                   - f_13 * ksi1_163[k]
                   + f_3 * pc_y[k] * ksk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, ksi0_164, ksi0_165, ksi0_166, ksi1_164, \
                         ksi1_165, ksi1_166, ksk_211, ksk_212, \
                         ksk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * ksi0_164[k]
                   - f_11 * ksi1_164[k]
                   + f_3 * pc_y[k] * ksk_211[k];

        t_265[k] = f_8 * ksi0_165[k]
                   - f_9 * ksi1_165[k]
                   + f_3 * pc_y[k] * ksk_212[k];

        t_266[k] = f_6 * ksi0_166[k]
                   - f_7 * ksi1_166[k]
                   + f_3 * pc_y[k] * ksk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, isk_107, isk_216, \
                         ksi0_167, ksi0_168, ksi1_167, ksi1_168, ksk_214, ksk_215, \
                         ksk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_4 * ksi0_167[k]
                   - f_5 * ksi1_167[k]
                   + f_3 * pc_y[k] * ksk_214[k];

        t_268[k] = f_3 * pc_y[k] * ksk_215[k];

        t_269[k] = f_16 * isk_107[k]
                   + f_1 * ksi0_167[k]
                   - f_2 * ksi1_167[k]
                   + f_3 * pc_z[k] * ksk_215[k];

        t_270[k] = f_18 * isk_216[k]
                   + f_1 * ksi0_168[k]
                   - f_2 * ksi1_168[k]
                   + f_3 * pc_x[k] * ksk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, isk_108, isk_219, \
                         ksi0_171, ksi1_171, ksk_216, ksk_217, \
                         ksk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * isk_108[k]
                   + f_3 * pc_y[k] * ksk_216[k];

        t_272[k] = f_3 * pc_z[k] * ksk_216[k];

        t_273[k] = f_18 * isk_219[k]
                   + f_12 * ksi0_171[k]
                   - f_13 * ksi1_171[k]
                   + f_3 * pc_x[k] * ksk_219[k];

        t_274[k] = f_3 * pc_z[k] * ksk_217[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, isk_222, ksi0_168, ksi0_174, \
                         ksi1_168, ksi1_174, ksk_218, ksk_219, \
                         ksk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_4 * ksi0_168[k]
                   - f_5 * ksi1_168[k]
                   + f_3 * pc_z[k] * ksk_218[k];

        t_276[k] = f_18 * isk_222[k]
                   + f_10 * ksi0_174[k]
                   - f_11 * ksi1_174[k]
                   + f_3 * pc_x[k] * ksk_222[k];

        t_277[k] = f_3 * pc_z[k] * ksk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pc_x, pc_y, pc_z, isk_113, isk_226, \
                         ksi0_170, ksi0_178, ksi1_170, ksi1_178, ksk_221, ksk_222, \
                         ksk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * isk_113[k]
                   + f_3 * pc_y[k] * ksk_221[k];

        t_279[k] = f_6 * ksi0_170[k]
                   - f_7 * ksi1_170[k]
                   + f_3 * pc_z[k] * ksk_221[k];

        t_280[k] = f_18 * isk_226[k]
                   + f_8 * ksi0_178[k]
                   - f_9 * ksi1_178[k]
                   + f_3 * pc_x[k] * ksk_226[k];

        t_281[k] = f_3 * pc_z[k] * ksk_222[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pc_y, pc_z, isk_117, ksi0_171, ksi0_173, \
                         ksi1_171, ksi1_173, ksk_223, ksk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_4 * ksi0_171[k]
                   - f_5 * ksi1_171[k]
                   + f_3 * pc_z[k] * ksk_223[k];

        t_283[k] = f_17 * isk_117[k]
                   + f_3 * pc_y[k] * ksk_225[k];

        t_284[k] = f_8 * ksi0_173[k]
                   - f_9 * ksi1_173[k]
                   + f_3 * pc_z[k] * ksk_225[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pc_x, pc_z, isk_231, ksi0_174, ksi0_183, \
                         ksi1_174, ksi1_183, ksk_226, ksk_227, \
                         ksk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_18 * isk_231[k]
                   + f_6 * ksi0_183[k]
                   - f_7 * ksi1_183[k]
                   + f_3 * pc_x[k] * ksk_231[k];

        t_286[k] = f_3 * pc_z[k] * ksk_226[k];

        t_287[k] = f_4 * ksi0_174[k]
                   - f_5 * ksi1_174[k]
                   + f_3 * pc_z[k] * ksk_227[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, isk_122, ksi0_175, ksi0_177, \
                         ksi1_175, ksi1_177, ksk_228, ksk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_6 * ksi0_175[k]
                   - f_7 * ksi1_175[k]
                   + f_3 * pc_z[k] * ksk_228[k];

        t_289[k] = f_17 * isk_122[k]
                   + f_3 * pc_y[k] * ksk_230[k];

        t_290[k] = f_10 * ksi0_177[k]
                   - f_11 * ksi1_177[k]
                   + f_3 * pc_z[k] * ksk_230[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pc_x, pc_z, isk_237, ksi0_178, ksi0_189, \
                         ksi1_178, ksi1_189, ksk_231, ksk_232, \
                         ksk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_18 * isk_237[k]
                   + f_4 * ksi0_189[k]
                   - f_5 * ksi1_189[k]
                   + f_3 * pc_x[k] * ksk_237[k];

        t_292[k] = f_3 * pc_z[k] * ksk_231[k];

        t_293[k] = f_4 * ksi0_178[k]
                   - f_5 * ksi1_178[k]
                   + f_3 * pc_z[k] * ksk_232[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_y, pc_z, isk_128, ksi0_179, ksi0_180, \
                         ksi0_182, ksi1_179, ksi1_180, ksi1_182, ksk_233, ksk_234, \
                         ksk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * ksi0_179[k]
                   - f_7 * ksi1_179[k]
                   + f_3 * pc_z[k] * ksk_233[k];

        t_295[k] = f_8 * ksi0_180[k]
                   - f_9 * ksi1_180[k]
                   + f_3 * pc_z[k] * ksk_234[k];

        t_296[k] = f_17 * isk_128[k]
                   + f_3 * pc_y[k] * ksk_236[k];

        t_297[k] = f_12 * ksi0_182[k]
                   - f_13 * ksi1_182[k]
                   + f_3 * pc_z[k] * ksk_236[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pc_x, pc_z, isk_244, isk_246, \
                         isk_247, isk_248, ksk_237, ksk_244, ksk_246, ksk_247, \
                         ksk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_18 * isk_244[k]
                   + f_3 * pc_x[k] * ksk_244[k];

        t_299[k] = f_3 * pc_z[k] * ksk_237[k];

        t_300[k] = f_18 * isk_246[k]
                   + f_3 * pc_x[k] * ksk_246[k];

        t_301[k] = f_18 * isk_247[k]
                   + f_3 * pc_x[k] * ksk_247[k];

        t_302[k] = f_18 * isk_248[k]
                   + f_3 * pc_x[k] * ksk_248[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, pc_y, isk_136, isk_249, isk_250, \
                         isk_251, ksi0_189, ksi1_189, ksk_244, ksk_249, ksk_250, \
                         ksk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_18 * isk_249[k]
                   + f_3 * pc_x[k] * ksk_249[k];

        t_304[k] = f_18 * isk_250[k]
                   + f_3 * pc_x[k] * ksk_250[k];

        t_305[k] = f_18 * isk_251[k]
                   + f_3 * pc_x[k] * ksk_251[k];

        t_306[k] = f_17 * isk_136[k]
                   + f_1 * ksi0_189[k]
                   - f_2 * ksi1_189[k]
                   + f_3 * pc_y[k] * ksk_244[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_z, ksi0_189, ksi0_190, ksi0_191, \
                         ksi1_189, ksi1_190, ksi1_191, ksk_244, ksk_245, ksk_246, \
                         ksk_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_z[k] * ksk_244[k];

        t_308[k] = f_4 * ksi0_189[k]
                   - f_5 * ksi1_189[k]
                   + f_3 * pc_z[k] * ksk_245[k];

        t_309[k] = f_6 * ksi0_190[k]
                   - f_7 * ksi1_190[k]
                   + f_3 * pc_z[k] * ksk_246[k];

        t_310[k] = f_8 * ksi0_191[k]
                   - f_9 * ksi1_191[k]
                   + f_3 * pc_z[k] * ksk_247[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, isk_143, ksi0_192, ksi0_193, \
                         ksi0_195, ksi1_192, ksi1_193, ksi1_195, ksk_248, ksk_249, \
                         ksk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_10 * ksi0_192[k]
                   - f_11 * ksi1_192[k]
                   + f_3 * pc_z[k] * ksk_248[k];

        t_312[k] = f_12 * ksi0_193[k]
                   - f_13 * ksi1_193[k]
                   + f_3 * pc_z[k] * ksk_249[k];

        t_313[k] = f_17 * isk_143[k]
                   + f_3 * pc_y[k] * ksk_251[k];

        t_314[k] = f_1 * ksi0_195[k]
                   - f_2 * ksi1_195[k]
                   + f_3 * pc_z[k] * ksk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_z, pc_y, pc_z, isl0_135, isl0_138, \
                         isk_108, isk_144, isl1_135, isl1_138, \
                         ksk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * isl0_135[k]
                   - f_14 * pc_z[k] * isl1_135[k];

        t_316[k] = f_16 * isk_144[k]
                   + f_3 * pc_y[k] * ksk_252[k];

        t_317[k] = f_15 * isk_108[k]
                   + f_3 * pc_z[k] * ksk_252[k];

        t_318[k] = pa_z[k] * isl0_138[k]
                   - f_14 * pc_z[k] * isl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_x, pc_y, pc_z, isl0_141, isk_146, \
                         isk_257, isl1_141, ksi0_201, ksi1_201, ksk_254, \
                         ksk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * isk_146[k]
                   + f_3 * pc_y[k] * ksk_254[k];

        t_320[k] = f_18 * isk_257[k]
                   + f_12 * ksi0_201[k]
                   - f_13 * ksi1_201[k]
                   + f_3 * pc_x[k] * ksk_257[k];

        t_321[k] = pa_z[k] * isl0_141[k]
                   - f_14 * pc_z[k] * isl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, isk_111, isk_149, isk_261, \
                         ksi0_205, ksi1_205, ksk_255, ksk_257, \
                         ksk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * isk_111[k]
                   + f_3 * pc_z[k] * ksk_255[k];

        t_323[k] = f_16 * isk_149[k]
                   + f_3 * pc_y[k] * ksk_257[k];

        t_324[k] = f_18 * isk_261[k]
                   + f_10 * ksi0_205[k]
                   - f_11 * ksi1_205[k]
                   + f_3 * pc_x[k] * ksk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pa_z, pc_y, pc_z, isl0_145, isl0_147, \
                         isk_114, isk_115, isk_153, isl1_145, isl1_147, ksk_258, \
                         ksk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_z[k] * isl0_145[k]
                   - f_14 * pc_z[k] * isl1_145[k];

        t_326[k] = f_15 * isk_114[k]
                   + f_3 * pc_z[k] * ksk_258[k];

        t_327[k] = pa_z[k] * isl0_147[k]
                   + f_16 * isk_115[k]
                   - f_14 * pc_z[k] * isl1_147[k];

        t_328[k] = f_16 * isk_153[k]
                   + f_3 * pc_y[k] * ksk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_z, pc_x, pc_z, isl0_150, isk_118, isk_266, \
                         isl1_150, ksi0_210, ksi1_210, ksk_262, \
                         ksk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_18 * isk_266[k]
                   + f_8 * ksi0_210[k]
                   - f_9 * ksi1_210[k]
                   + f_3 * pc_x[k] * ksk_266[k];

        t_330[k] = pa_z[k] * isl0_150[k]
                   - f_14 * pc_z[k] * isl1_150[k];

        t_331[k] = f_15 * isk_118[k]
                   + f_3 * pc_z[k] * ksk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pa_z, pc_y, pc_z, isl0_152, isl0_153, isk_119, \
                         isk_120, isk_158, isl1_152, isl1_153, \
                         ksk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pa_z[k] * isl0_152[k]
                   + f_16 * isk_119[k]
                   - f_14 * pc_z[k] * isl1_152[k];

        t_333[k] = pa_z[k] * isl0_153[k]
                   + f_17 * isk_120[k]
                   - f_14 * pc_z[k] * isl1_153[k];

        t_334[k] = f_16 * isk_158[k]
                   + f_3 * pc_y[k] * ksk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pa_z, pc_x, pc_z, isl0_156, isk_123, isk_272, \
                         isl1_156, ksi0_216, ksi1_216, ksk_267, \
                         ksk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_18 * isk_272[k]
                   + f_6 * ksi0_216[k]
                   - f_7 * ksi1_216[k]
                   + f_3 * pc_x[k] * ksk_272[k];

        t_336[k] = pa_z[k] * isl0_156[k]
                   - f_14 * pc_z[k] * isl1_156[k];

        t_337[k] = f_15 * isk_123[k]
                   + f_3 * pc_z[k] * ksk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pa_z, pc_z, isl0_158, isl0_159, isl0_160, \
                         isk_124, isk_125, isk_126, isl1_158, isl1_159, \
                         isl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pa_z[k] * isl0_158[k]
                   + f_16 * isk_124[k]
                   - f_14 * pc_z[k] * isl1_158[k];

        t_339[k] = pa_z[k] * isl0_159[k]
                   + f_17 * isk_125[k]
                   - f_14 * pc_z[k] * isl1_159[k];

        t_340[k] = pa_z[k] * isl0_160[k]
                   + f_18 * isk_126[k]
                   - f_14 * pc_z[k] * isl1_160[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, isk_164, isk_279, isk_280, \
                         isk_281, ksi0_223, ksi1_223, ksk_272, ksk_279, ksk_280, \
                         ksk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * isk_164[k]
                   + f_3 * pc_y[k] * ksk_272[k];

        t_342[k] = f_18 * isk_279[k]
                   + f_4 * ksi0_223[k]
                   - f_5 * ksi1_223[k]
                   + f_3 * pc_x[k] * ksk_279[k];

        t_343[k] = f_18 * isk_280[k]
                   + f_3 * pc_x[k] * ksk_280[k];

        t_344[k] = f_18 * isk_281[k]
                   + f_3 * pc_x[k] * ksk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, isk_282, isk_283, isk_284, \
                         isk_285, isk_286, ksk_282, ksk_283, ksk_284, ksk_285, \
                         ksk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_18 * isk_282[k]
                   + f_3 * pc_x[k] * ksk_282[k];

        t_346[k] = f_18 * isk_283[k]
                   + f_3 * pc_x[k] * ksk_283[k];

        t_347[k] = f_18 * isk_284[k]
                   + f_3 * pc_x[k] * ksk_284[k];

        t_348[k] = f_18 * isk_285[k]
                   + f_3 * pc_x[k] * ksk_285[k];

        t_349[k] = f_18 * isk_286[k]
                   + f_3 * pc_x[k] * ksk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pa_z, pc_x, pc_z, isl0_171, isk_136, isk_287, \
                         isl1_171, ksk_280, ksk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_18 * isk_287[k]
                   + f_3 * pc_x[k] * ksk_287[k];

        t_351[k] = pa_z[k] * isl0_171[k]
                   - f_14 * pc_z[k] * isl1_171[k];

        t_352[k] = f_15 * isk_136[k]
                   + f_3 * pc_z[k] * ksk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, isk_174, isk_175, isk_176, ksi0_219, \
                         ksi0_220, ksi0_221, ksi1_219, ksi1_220, ksi1_221, ksk_282, ksk_283, \
                         ksk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * isk_174[k]
                   + f_12 * ksi0_219[k]
                   - f_13 * ksi1_219[k]
                   + f_3 * pc_y[k] * ksk_282[k];

        t_354[k] = f_16 * isk_175[k]
                   + f_10 * ksi0_220[k]
                   - f_11 * ksi1_220[k]
                   + f_3 * pc_y[k] * ksk_283[k];

        t_355[k] = f_16 * isk_176[k]
                   + f_8 * ksi0_221[k]
                   - f_9 * ksi1_221[k]
                   + f_3 * pc_y[k] * ksk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, isk_177, isk_178, isk_179, ksi0_222, \
                         ksi0_223, ksi1_222, ksi1_223, ksk_285, ksk_286, \
                         ksk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * isk_177[k]
                   + f_6 * ksi0_222[k]
                   - f_7 * ksi1_222[k]
                   + f_3 * pc_y[k] * ksk_285[k];

        t_357[k] = f_16 * isk_178[k]
                   + f_4 * ksi0_223[k]
                   - f_5 * ksi1_223[k]
                   + f_3 * pc_y[k] * ksk_286[k];

        t_358[k] = f_16 * isk_179[k]
                   + f_3 * pc_y[k] * ksk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_y, pc_y, pc_z, isl0_225, isk_143, \
                         isk_144, isk_180, isl1_225, ksi0_223, ksi1_223, ksk_287, \
                         ksk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * isk_143[k]
                   + f_1 * ksi0_223[k]
                   - f_2 * ksi1_223[k]
                   + f_3 * pc_z[k] * ksk_287[k];

        t_360[k] = pa_y[k] * isl0_225[k]
                   - f_14 * pc_y[k] * isl1_225[k];

        t_361[k] = f_15 * isk_180[k]
                   + f_3 * pc_y[k] * ksk_288[k];

        t_362[k] = f_16 * isk_144[k]
                   + f_3 * pc_z[k] * ksk_288[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
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
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);

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

    const auto *isl0_228 = buffer.data(isl0 + 228);
    const auto *isl0_230 = buffer.data(isl0 + 230);
    const auto *isl0_231 = buffer.data(isl0 + 231);
    const auto *isl0_234 = buffer.data(isl0 + 234);
    const auto *isl0_235 = buffer.data(isl0 + 235);
    const auto *isl0_237 = buffer.data(isl0 + 237);
    const auto *isl0_239 = buffer.data(isl0 + 239);
    const auto *isl0_240 = buffer.data(isl0 + 240);
    const auto *isl0_242 = buffer.data(isl0 + 242);
    const auto *isl0_243 = buffer.data(isl0 + 243);
    const auto *isl0_245 = buffer.data(isl0 + 245);
    const auto *isl0_246 = buffer.data(isl0 + 246);
    const auto *isl0_248 = buffer.data(isl0 + 248);
    const auto *isl0_249 = buffer.data(isl0 + 249);
    const auto *isl0_250 = buffer.data(isl0 + 250);
    const auto *isl0_252 = buffer.data(isl0 + 252);
    const auto *isl0_269 = buffer.data(isl0 + 269);

    const auto *isk_147 = buffer.data(isk + 147);
    const auto *isk_150 = buffer.data(isk + 150);
    const auto *isk_154 = buffer.data(isk + 154);
    const auto *isk_159 = buffer.data(isk + 159);
    const auto *isk_172 = buffer.data(isk + 172);
    const auto *isk_180 = buffer.data(isk + 180);
    const auto *isk_181 = buffer.data(isk + 181);
    const auto *isk_182 = buffer.data(isk + 182);
    const auto *isk_183 = buffer.data(isk + 183);
    const auto *isk_185 = buffer.data(isk + 185);
    const auto *isk_186 = buffer.data(isk + 186);
    const auto *isk_188 = buffer.data(isk + 188);
    const auto *isk_189 = buffer.data(isk + 189);
    const auto *isk_190 = buffer.data(isk + 190);
    const auto *isk_192 = buffer.data(isk + 192);
    const auto *isk_193 = buffer.data(isk + 193);
    const auto *isk_194 = buffer.data(isk + 194);
    const auto *isk_195 = buffer.data(isk + 195);
    const auto *isk_197 = buffer.data(isk + 197);
    const auto *isk_198 = buffer.data(isk + 198);
    const auto *isk_199 = buffer.data(isk + 199);
    const auto *isk_200 = buffer.data(isk + 200);
    const auto *isk_208 = buffer.data(isk + 208);
    const auto *isk_210 = buffer.data(isk + 210);
    const auto *isk_211 = buffer.data(isk + 211);
    const auto *isk_212 = buffer.data(isk + 212);
    const auto *isk_213 = buffer.data(isk + 213);
    const auto *isk_214 = buffer.data(isk + 214);
    const auto *isk_215 = buffer.data(isk + 215);
    const auto *isk_216 = buffer.data(isk + 216);
    const auto *isk_221 = buffer.data(isk + 221);
    const auto *isk_225 = buffer.data(isk + 225);
    const auto *isk_230 = buffer.data(isk + 230);
    const auto *isk_236 = buffer.data(isk + 236);
    const auto *isk_316 = buffer.data(isk + 316);
    const auto *isk_317 = buffer.data(isk + 317);
    const auto *isk_318 = buffer.data(isk + 318);
    const auto *isk_319 = buffer.data(isk + 319);
    const auto *isk_320 = buffer.data(isk + 320);
    const auto *isk_321 = buffer.data(isk + 321);
    const auto *isk_322 = buffer.data(isk + 322);
    const auto *isk_323 = buffer.data(isk + 323);
    const auto *isk_324 = buffer.data(isk + 324);
    const auto *isk_329 = buffer.data(isk + 329);
    const auto *isk_333 = buffer.data(isk + 333);
    const auto *isk_338 = buffer.data(isk + 338);
    const auto *isk_344 = buffer.data(isk + 344);
    const auto *isk_351 = buffer.data(isk + 351);
    const auto *isk_352 = buffer.data(isk + 352);
    const auto *isk_353 = buffer.data(isk + 353);
    const auto *isk_354 = buffer.data(isk + 354);
    const auto *isk_355 = buffer.data(isk + 355);
    const auto *isk_356 = buffer.data(isk + 356);
    const auto *isk_357 = buffer.data(isk + 357);
    const auto *isk_359 = buffer.data(isk + 359);
    const auto *isk_360 = buffer.data(isk + 360);
    const auto *isk_363 = buffer.data(isk + 363);
    const auto *isk_366 = buffer.data(isk + 366);
    const auto *isk_370 = buffer.data(isk + 370);
    const auto *isk_375 = buffer.data(isk + 375);
    const auto *isk_381 = buffer.data(isk + 381);

    const auto *isl1_228 = buffer.data(isl1 + 228);
    const auto *isl1_230 = buffer.data(isl1 + 230);
    const auto *isl1_231 = buffer.data(isl1 + 231);
    const auto *isl1_234 = buffer.data(isl1 + 234);
    const auto *isl1_235 = buffer.data(isl1 + 235);
    const auto *isl1_237 = buffer.data(isl1 + 237);
    const auto *isl1_239 = buffer.data(isl1 + 239);
    const auto *isl1_240 = buffer.data(isl1 + 240);
    const auto *isl1_242 = buffer.data(isl1 + 242);
    const auto *isl1_243 = buffer.data(isl1 + 243);
    const auto *isl1_245 = buffer.data(isl1 + 245);
    const auto *isl1_246 = buffer.data(isl1 + 246);
    const auto *isl1_248 = buffer.data(isl1 + 248);
    const auto *isl1_249 = buffer.data(isl1 + 249);
    const auto *isl1_250 = buffer.data(isl1 + 250);
    const auto *isl1_252 = buffer.data(isl1 + 252);
    const auto *isl1_269 = buffer.data(isl1 + 269);

    const auto *ksi0_245 = buffer.data(ksi0 + 245);
    const auto *ksi0_247 = buffer.data(ksi0 + 247);
    const auto *ksi0_248 = buffer.data(ksi0 + 248);
    const auto *ksi0_249 = buffer.data(ksi0 + 249);
    const auto *ksi0_250 = buffer.data(ksi0 + 250);
    const auto *ksi0_251 = buffer.data(ksi0 + 251);
    const auto *ksi0_252 = buffer.data(ksi0 + 252);
    const auto *ksi0_253 = buffer.data(ksi0 + 253);
    const auto *ksi0_254 = buffer.data(ksi0 + 254);
    const auto *ksi0_255 = buffer.data(ksi0 + 255);
    const auto *ksi0_256 = buffer.data(ksi0 + 256);
    const auto *ksi0_257 = buffer.data(ksi0 + 257);
    const auto *ksi0_258 = buffer.data(ksi0 + 258);
    const auto *ksi0_259 = buffer.data(ksi0 + 259);
    const auto *ksi0_260 = buffer.data(ksi0 + 260);
    const auto *ksi0_261 = buffer.data(ksi0 + 261);
    const auto *ksi0_262 = buffer.data(ksi0 + 262);
    const auto *ksi0_263 = buffer.data(ksi0 + 263);
    const auto *ksi0_264 = buffer.data(ksi0 + 264);
    const auto *ksi0_265 = buffer.data(ksi0 + 265);
    const auto *ksi0_266 = buffer.data(ksi0 + 266);
    const auto *ksi0_272 = buffer.data(ksi0 + 272);
    const auto *ksi0_273 = buffer.data(ksi0 + 273);
    const auto *ksi0_274 = buffer.data(ksi0 + 274);
    const auto *ksi0_275 = buffer.data(ksi0 + 275);
    const auto *ksi0_276 = buffer.data(ksi0 + 276);
    const auto *ksi0_277 = buffer.data(ksi0 + 277);
    const auto *ksi0_278 = buffer.data(ksi0 + 278);
    const auto *ksi0_279 = buffer.data(ksi0 + 279);
    const auto *ksi0_280 = buffer.data(ksi0 + 280);
    const auto *ksi0_282 = buffer.data(ksi0 + 282);
    const auto *ksi0_283 = buffer.data(ksi0 + 283);
    const auto *ksi0_285 = buffer.data(ksi0 + 285);
    const auto *ksi0_286 = buffer.data(ksi0 + 286);
    const auto *ksi0_287 = buffer.data(ksi0 + 287);
    const auto *ksi0_289 = buffer.data(ksi0 + 289);
    const auto *ksi0_290 = buffer.data(ksi0 + 290);
    const auto *ksi0_291 = buffer.data(ksi0 + 291);
    const auto *ksi0_292 = buffer.data(ksi0 + 292);
    const auto *ksi0_294 = buffer.data(ksi0 + 294);
    const auto *ksi0_295 = buffer.data(ksi0 + 295);
    const auto *ksi0_301 = buffer.data(ksi0 + 301);

    const auto *ksi1_245 = buffer.data(ksi1 + 245);
    const auto *ksi1_247 = buffer.data(ksi1 + 247);
    const auto *ksi1_248 = buffer.data(ksi1 + 248);
    const auto *ksi1_249 = buffer.data(ksi1 + 249);
    const auto *ksi1_250 = buffer.data(ksi1 + 250);
    const auto *ksi1_251 = buffer.data(ksi1 + 251);
    const auto *ksi1_252 = buffer.data(ksi1 + 252);
    const auto *ksi1_253 = buffer.data(ksi1 + 253);
    const auto *ksi1_254 = buffer.data(ksi1 + 254);
    const auto *ksi1_255 = buffer.data(ksi1 + 255);
    const auto *ksi1_256 = buffer.data(ksi1 + 256);
    const auto *ksi1_257 = buffer.data(ksi1 + 257);
    const auto *ksi1_258 = buffer.data(ksi1 + 258);
    const auto *ksi1_259 = buffer.data(ksi1 + 259);
    const auto *ksi1_260 = buffer.data(ksi1 + 260);
    const auto *ksi1_261 = buffer.data(ksi1 + 261);
    const auto *ksi1_262 = buffer.data(ksi1 + 262);
    const auto *ksi1_263 = buffer.data(ksi1 + 263);
    const auto *ksi1_264 = buffer.data(ksi1 + 264);
    const auto *ksi1_265 = buffer.data(ksi1 + 265);
    const auto *ksi1_266 = buffer.data(ksi1 + 266);
    const auto *ksi1_272 = buffer.data(ksi1 + 272);
    const auto *ksi1_273 = buffer.data(ksi1 + 273);
    const auto *ksi1_274 = buffer.data(ksi1 + 274);
    const auto *ksi1_275 = buffer.data(ksi1 + 275);
    const auto *ksi1_276 = buffer.data(ksi1 + 276);
    const auto *ksi1_277 = buffer.data(ksi1 + 277);
    const auto *ksi1_278 = buffer.data(ksi1 + 278);
    const auto *ksi1_279 = buffer.data(ksi1 + 279);
    const auto *ksi1_280 = buffer.data(ksi1 + 280);
    const auto *ksi1_282 = buffer.data(ksi1 + 282);
    const auto *ksi1_283 = buffer.data(ksi1 + 283);
    const auto *ksi1_285 = buffer.data(ksi1 + 285);
    const auto *ksi1_286 = buffer.data(ksi1 + 286);
    const auto *ksi1_287 = buffer.data(ksi1 + 287);
    const auto *ksi1_289 = buffer.data(ksi1 + 289);
    const auto *ksi1_290 = buffer.data(ksi1 + 290);
    const auto *ksi1_291 = buffer.data(ksi1 + 291);
    const auto *ksi1_292 = buffer.data(ksi1 + 292);
    const auto *ksi1_294 = buffer.data(ksi1 + 294);
    const auto *ksi1_295 = buffer.data(ksi1 + 295);
    const auto *ksi1_301 = buffer.data(ksi1 + 301);

    const auto *ksk_290 = buffer.data(ksk + 290);
    const auto *ksk_291 = buffer.data(ksk + 291);
    const auto *ksk_293 = buffer.data(ksk + 293);
    const auto *ksk_294 = buffer.data(ksk + 294);
    const auto *ksk_297 = buffer.data(ksk + 297);
    const auto *ksk_298 = buffer.data(ksk + 298);
    const auto *ksk_302 = buffer.data(ksk + 302);
    const auto *ksk_303 = buffer.data(ksk + 303);
    const auto *ksk_308 = buffer.data(ksk + 308);
    const auto *ksk_316 = buffer.data(ksk + 316);
    const auto *ksk_317 = buffer.data(ksk + 317);
    const auto *ksk_318 = buffer.data(ksk + 318);
    const auto *ksk_319 = buffer.data(ksk + 319);
    const auto *ksk_320 = buffer.data(ksk + 320);
    const auto *ksk_321 = buffer.data(ksk + 321);
    const auto *ksk_322 = buffer.data(ksk + 322);
    const auto *ksk_323 = buffer.data(ksk + 323);
    const auto *ksk_324 = buffer.data(ksk + 324);
    const auto *ksk_325 = buffer.data(ksk + 325);
    const auto *ksk_326 = buffer.data(ksk + 326);
    const auto *ksk_327 = buffer.data(ksk + 327);
    const auto *ksk_328 = buffer.data(ksk + 328);
    const auto *ksk_329 = buffer.data(ksk + 329);
    const auto *ksk_330 = buffer.data(ksk + 330);
    const auto *ksk_331 = buffer.data(ksk + 331);
    const auto *ksk_332 = buffer.data(ksk + 332);
    const auto *ksk_333 = buffer.data(ksk + 333);
    const auto *ksk_334 = buffer.data(ksk + 334);
    const auto *ksk_335 = buffer.data(ksk + 335);
    const auto *ksk_336 = buffer.data(ksk + 336);
    const auto *ksk_337 = buffer.data(ksk + 337);
    const auto *ksk_338 = buffer.data(ksk + 338);
    const auto *ksk_339 = buffer.data(ksk + 339);
    const auto *ksk_340 = buffer.data(ksk + 340);
    const auto *ksk_341 = buffer.data(ksk + 341);
    const auto *ksk_342 = buffer.data(ksk + 342);
    const auto *ksk_343 = buffer.data(ksk + 343);
    const auto *ksk_344 = buffer.data(ksk + 344);
    const auto *ksk_351 = buffer.data(ksk + 351);
    const auto *ksk_352 = buffer.data(ksk + 352);
    const auto *ksk_353 = buffer.data(ksk + 353);
    const auto *ksk_354 = buffer.data(ksk + 354);
    const auto *ksk_355 = buffer.data(ksk + 355);
    const auto *ksk_356 = buffer.data(ksk + 356);
    const auto *ksk_357 = buffer.data(ksk + 357);
    const auto *ksk_358 = buffer.data(ksk + 358);
    const auto *ksk_359 = buffer.data(ksk + 359);
    const auto *ksk_360 = buffer.data(ksk + 360);
    const auto *ksk_361 = buffer.data(ksk + 361);
    const auto *ksk_362 = buffer.data(ksk + 362);
    const auto *ksk_363 = buffer.data(ksk + 363);
    const auto *ksk_365 = buffer.data(ksk + 365);
    const auto *ksk_366 = buffer.data(ksk + 366);
    const auto *ksk_367 = buffer.data(ksk + 367);
    const auto *ksk_369 = buffer.data(ksk + 369);
    const auto *ksk_370 = buffer.data(ksk + 370);
    const auto *ksk_371 = buffer.data(ksk + 371);
    const auto *ksk_372 = buffer.data(ksk + 372);
    const auto *ksk_374 = buffer.data(ksk + 374);
    const auto *ksk_375 = buffer.data(ksk + 375);
    const auto *ksk_376 = buffer.data(ksk + 376);
    const auto *ksk_377 = buffer.data(ksk + 377);
    const auto *ksk_378 = buffer.data(ksk + 378);
    const auto *ksk_380 = buffer.data(ksk + 380);
    const auto *ksk_381 = buffer.data(ksk + 381);

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_y, pc_y, isl0_228, isl0_230, isl0_231, \
                         isk_181, isk_182, isk_183, isl1_228, isl1_230, isl1_231, \
                         ksk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_y[k] * isl0_228[k]
                   + f_16 * isk_181[k]
                   - f_14 * pc_y[k] * isl1_228[k];

        t_364[k] = f_15 * isk_182[k]
                   + f_3 * pc_y[k] * ksk_290[k];

        t_365[k] = pa_y[k] * isl0_230[k]
                   - f_14 * pc_y[k] * isl1_230[k];

        t_366[k] = pa_y[k] * isl0_231[k]
                   + f_17 * isk_183[k]
                   - f_14 * pc_y[k] * isl1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_y, pc_y, pc_z, isl0_234, isl0_235, \
                         isk_147, isk_185, isk_186, isl1_234, isl1_235, ksk_291, \
                         ksk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * isk_147[k]
                   + f_3 * pc_z[k] * ksk_291[k];

        t_368[k] = f_15 * isk_185[k]
                   + f_3 * pc_y[k] * ksk_293[k];

        t_369[k] = pa_y[k] * isl0_234[k]
                   - f_14 * pc_y[k] * isl1_234[k];

        t_370[k] = pa_y[k] * isl0_235[k]
                   + f_18 * isk_186[k]
                   - f_14 * pc_y[k] * isl1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pc_y, pc_z, isl0_237, isl0_239, \
                         isk_150, isk_188, isk_189, isl1_237, isl1_239, ksk_294, \
                         ksk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * isk_150[k]
                   + f_3 * pc_z[k] * ksk_294[k];

        t_372[k] = pa_y[k] * isl0_237[k]
                   + f_16 * isk_188[k]
                   - f_14 * pc_y[k] * isl1_237[k];

        t_373[k] = f_15 * isk_189[k]
                   + f_3 * pc_y[k] * ksk_297[k];

        t_374[k] = pa_y[k] * isl0_239[k]
                   - f_14 * pc_y[k] * isl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_y, pc_y, pc_z, isl0_240, isl0_242, isk_154, \
                         isk_190, isk_192, isl1_240, isl1_242, \
                         ksk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pa_y[k] * isl0_240[k]
                   + f_19 * isk_190[k]
                   - f_14 * pc_y[k] * isl1_240[k];

        t_376[k] = f_16 * isk_154[k]
                   + f_3 * pc_z[k] * ksk_298[k];

        t_377[k] = pa_y[k] * isl0_242[k]
                   + f_17 * isk_192[k]
                   - f_14 * pc_y[k] * isl1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_y, pc_y, isl0_243, isl0_245, isl0_246, \
                         isk_193, isk_194, isk_195, isl1_243, isl1_245, isl1_246, \
                         ksk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * isl0_243[k]
                   + f_16 * isk_193[k]
                   - f_14 * pc_y[k] * isl1_243[k];

        t_379[k] = f_15 * isk_194[k]
                   + f_3 * pc_y[k] * ksk_302[k];

        t_380[k] = pa_y[k] * isl0_245[k]
                   - f_14 * pc_y[k] * isl1_245[k];

        t_381[k] = pa_y[k] * isl0_246[k]
                   + f_20 * isk_195[k]
                   - f_14 * pc_y[k] * isl1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pa_y, pc_y, pc_z, isl0_248, isl0_249, isk_159, \
                         isk_197, isk_198, isl1_248, isl1_249, \
                         ksk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * isk_159[k]
                   + f_3 * pc_z[k] * ksk_303[k];

        t_383[k] = pa_y[k] * isl0_248[k]
                   + f_18 * isk_197[k]
                   - f_14 * pc_y[k] * isl1_248[k];

        t_384[k] = pa_y[k] * isl0_249[k]
                   + f_17 * isk_198[k]
                   - f_14 * pc_y[k] * isl1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_y, pc_x, pc_y, isl0_250, isl0_252, \
                         isk_199, isk_200, isk_316, isl1_250, isl1_252, ksk_308, \
                         ksk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pa_y[k] * isl0_250[k]
                   + f_16 * isk_199[k]
                   - f_14 * pc_y[k] * isl1_250[k];

        t_386[k] = f_15 * isk_200[k]
                   + f_3 * pc_y[k] * ksk_308[k];

        t_387[k] = pa_y[k] * isl0_252[k]
                   - f_14 * pc_y[k] * isl1_252[k];

        t_388[k] = f_18 * isk_316[k]
                   + f_3 * pc_x[k] * ksk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, isk_317, isk_318, isk_319, \
                         isk_320, isk_321, ksk_317, ksk_318, ksk_319, ksk_320, \
                         ksk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_18 * isk_317[k]
                   + f_3 * pc_x[k] * ksk_317[k];

        t_390[k] = f_18 * isk_318[k]
                   + f_3 * pc_x[k] * ksk_318[k];

        t_391[k] = f_18 * isk_319[k]
                   + f_3 * pc_x[k] * ksk_319[k];

        t_392[k] = f_18 * isk_320[k]
                   + f_3 * pc_x[k] * ksk_320[k];

        t_393[k] = f_18 * isk_321[k]
                   + f_3 * pc_x[k] * ksk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, isk_172, isk_208, \
                         isk_322, isk_323, ksi0_245, ksi1_245, ksk_316, ksk_322, \
                         ksk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_18 * isk_322[k]
                   + f_3 * pc_x[k] * ksk_322[k];

        t_395[k] = f_18 * isk_323[k]
                   + f_3 * pc_x[k] * ksk_323[k];

        t_396[k] = f_15 * isk_208[k]
                   + f_1 * ksi0_245[k]
                   - f_2 * ksi1_245[k]
                   + f_3 * pc_y[k] * ksk_316[k];

        t_397[k] = f_16 * isk_172[k]
                   + f_3 * pc_z[k] * ksk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, isk_210, isk_211, isk_212, ksi0_247, \
                         ksi0_248, ksi0_249, ksi1_247, ksi1_248, ksi1_249, ksk_318, ksk_319, \
                         ksk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * isk_210[k]
                   + f_12 * ksi0_247[k]
                   - f_13 * ksi1_247[k]
                   + f_3 * pc_y[k] * ksk_318[k];

        t_399[k] = f_15 * isk_211[k]
                   + f_10 * ksi0_248[k]
                   - f_11 * ksi1_248[k]
                   + f_3 * pc_y[k] * ksk_319[k];

        t_400[k] = f_15 * isk_212[k]
                   + f_8 * ksi0_249[k]
                   - f_9 * ksi1_249[k]
                   + f_3 * pc_y[k] * ksk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, isk_213, isk_214, isk_215, ksi0_250, \
                         ksi0_251, ksi1_250, ksi1_251, ksk_321, ksk_322, \
                         ksk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * isk_213[k]
                   + f_6 * ksi0_250[k]
                   - f_7 * ksi1_250[k]
                   + f_3 * pc_y[k] * ksk_321[k];

        t_402[k] = f_15 * isk_214[k]
                   + f_4 * ksi0_251[k]
                   - f_5 * ksi1_251[k]
                   + f_3 * pc_y[k] * ksk_322[k];

        t_403[k] = f_15 * isk_215[k]
                   + f_3 * pc_y[k] * ksk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_y, pc_x, pc_y, pc_z, isl0_269, \
                         isk_180, isk_324, isl1_269, ksi0_252, ksi1_252, \
                         ksk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * isl0_269[k]
                   - f_14 * pc_y[k] * isl1_269[k];

        t_405[k] = f_18 * isk_324[k]
                   + f_1 * ksi0_252[k]
                   - f_2 * ksi1_252[k]
                   + f_3 * pc_x[k] * ksk_324[k];

        t_406[k] = f_3 * pc_y[k] * ksk_324[k];

        t_407[k] = f_17 * isk_180[k]
                   + f_3 * pc_z[k] * ksk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, isk_329, ksi0_252, ksi0_257, \
                         ksi1_252, ksi1_257, ksk_325, ksk_326, \
                         ksk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_4 * ksi0_252[k]
                   - f_5 * ksi1_252[k]
                   + f_3 * pc_y[k] * ksk_325[k];

        t_409[k] = f_3 * pc_y[k] * ksk_326[k];

        t_410[k] = f_18 * isk_329[k]
                   + f_12 * ksi0_257[k]
                   - f_13 * ksi1_257[k]
                   + f_3 * pc_x[k] * ksk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_y, ksi0_253, ksi0_254, ksi1_253, ksi1_254, \
                         ksk_327, ksk_328, ksk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_6 * ksi0_253[k]
                   - f_7 * ksi1_253[k]
                   + f_3 * pc_y[k] * ksk_327[k];

        t_412[k] = f_4 * ksi0_254[k]
                   - f_5 * ksi1_254[k]
                   + f_3 * pc_y[k] * ksk_328[k];

        t_413[k] = f_3 * pc_y[k] * ksk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, isk_333, ksi0_255, ksi0_256, \
                         ksi0_261, ksi1_255, ksi1_256, ksi1_261, ksk_330, ksk_331, \
                         ksk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_18 * isk_333[k]
                   + f_10 * ksi0_261[k]
                   - f_11 * ksi1_261[k]
                   + f_3 * pc_x[k] * ksk_333[k];

        t_415[k] = f_8 * ksi0_255[k]
                   - f_9 * ksi1_255[k]
                   + f_3 * pc_y[k] * ksk_330[k];

        t_416[k] = f_6 * ksi0_256[k]
                   - f_7 * ksi1_256[k]
                   + f_3 * pc_y[k] * ksk_331[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, isk_338, ksi0_257, ksi0_266, \
                         ksi1_257, ksi1_266, ksk_332, ksk_333, \
                         ksk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * ksi0_257[k]
                   - f_5 * ksi1_257[k]
                   + f_3 * pc_y[k] * ksk_332[k];

        t_418[k] = f_3 * pc_y[k] * ksk_333[k];

        t_419[k] = f_18 * isk_338[k]
                   + f_8 * ksi0_266[k]
                   - f_9 * ksi1_266[k]
                   + f_3 * pc_x[k] * ksk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_y, ksi0_258, ksi0_259, ksi0_260, ksi1_258, \
                         ksi1_259, ksi1_260, ksk_334, ksk_335, \
                         ksk_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_10 * ksi0_258[k]
                   - f_11 * ksi1_258[k]
                   + f_3 * pc_y[k] * ksk_334[k];

        t_421[k] = f_8 * ksi0_259[k]
                   - f_9 * ksi1_259[k]
                   + f_3 * pc_y[k] * ksk_335[k];

        t_422[k] = f_6 * ksi0_260[k]
                   - f_7 * ksi1_260[k]
                   + f_3 * pc_y[k] * ksk_336[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, isk_344, ksi0_261, ksi0_272, \
                         ksi1_261, ksi1_272, ksk_337, ksk_338, \
                         ksk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * ksi0_261[k]
                   - f_5 * ksi1_261[k]
                   + f_3 * pc_y[k] * ksk_337[k];

        t_424[k] = f_3 * pc_y[k] * ksk_338[k];

        t_425[k] = f_18 * isk_344[k]
                   + f_6 * ksi0_272[k]
                   - f_7 * ksi1_272[k]
                   + f_3 * pc_x[k] * ksk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_y, ksi0_262, ksi0_263, ksi0_264, ksi1_262, \
                         ksi1_263, ksi1_264, ksk_339, ksk_340, \
                         ksk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_12 * ksi0_262[k]
                   - f_13 * ksi1_262[k]
                   + f_3 * pc_y[k] * ksk_339[k];

        t_427[k] = f_10 * ksi0_263[k]
                   - f_11 * ksi1_263[k]
                   + f_3 * pc_y[k] * ksk_340[k];

        t_428[k] = f_8 * ksi0_264[k]
                   - f_9 * ksi1_264[k]
                   + f_3 * pc_y[k] * ksk_341[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_y, ksi0_265, ksi0_266, ksi1_265, ksi1_266, \
                         ksk_342, ksk_343, ksk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_6 * ksi0_265[k]
                   - f_7 * ksi1_265[k]
                   + f_3 * pc_y[k] * ksk_342[k];

        t_430[k] = f_4 * ksi0_266[k]
                   - f_5 * ksi1_266[k]
                   + f_3 * pc_y[k] * ksk_343[k];

        t_431[k] = f_3 * pc_y[k] * ksk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, isk_351, isk_352, isk_353, isk_354, \
                         ksi0_279, ksi1_279, ksk_351, ksk_352, ksk_353, \
                         ksk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_18 * isk_351[k]
                   + f_4 * ksi0_279[k]
                   - f_5 * ksi1_279[k]
                   + f_3 * pc_x[k] * ksk_351[k];

        t_433[k] = f_18 * isk_352[k]
                   + f_3 * pc_x[k] * ksk_352[k];

        t_434[k] = f_18 * isk_353[k]
                   + f_3 * pc_x[k] * ksk_353[k];

        t_435[k] = f_18 * isk_354[k]
                   + f_3 * pc_x[k] * ksk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_y, isk_355, isk_356, \
                         isk_357, isk_359, ksk_351, ksk_355, ksk_356, ksk_357, \
                         ksk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_18 * isk_355[k]
                   + f_3 * pc_x[k] * ksk_355[k];

        t_437[k] = f_18 * isk_356[k]
                   + f_3 * pc_x[k] * ksk_356[k];

        t_438[k] = f_18 * isk_357[k]
                   + f_3 * pc_x[k] * ksk_357[k];

        t_439[k] = f_3 * pc_y[k] * ksk_351[k];

        t_440[k] = f_18 * isk_359[k]
                   + f_3 * pc_x[k] * ksk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, ksi0_273, ksi0_274, ksi0_275, ksi1_273, \
                         ksi1_274, ksi1_275, ksk_352, ksk_353, \
                         ksk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * ksi0_273[k]
                   - f_2 * ksi1_273[k]
                   + f_3 * pc_y[k] * ksk_352[k];

        t_442[k] = f_21 * ksi0_274[k]
                   - f_22 * ksi1_274[k]
                   + f_3 * pc_y[k] * ksk_353[k];

        t_443[k] = f_12 * ksi0_275[k]
                   - f_13 * ksi1_275[k]
                   + f_3 * pc_y[k] * ksk_354[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_y, ksi0_276, ksi0_277, ksi0_278, ksi1_276, \
                         ksi1_277, ksi1_278, ksk_355, ksk_356, \
                         ksk_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_10 * ksi0_276[k]
                   - f_11 * ksi1_276[k]
                   + f_3 * pc_y[k] * ksk_355[k];

        t_445[k] = f_8 * ksi0_277[k]
                   - f_9 * ksi1_277[k]
                   + f_3 * pc_y[k] * ksk_356[k];

        t_446[k] = f_6 * ksi0_278[k]
                   - f_7 * ksi1_278[k]
                   + f_3 * pc_y[k] * ksk_357[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pc_x, pc_y, pc_z, isk_215, isk_360, \
                         ksi0_279, ksi0_280, ksi1_279, ksi1_280, ksk_358, ksk_359, \
                         ksk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_4 * ksi0_279[k]
                   - f_5 * ksi1_279[k]
                   + f_3 * pc_y[k] * ksk_358[k];

        t_448[k] = f_3 * pc_y[k] * ksk_359[k];

        t_449[k] = f_17 * isk_215[k]
                   + f_1 * ksi0_279[k]
                   - f_2 * ksi1_279[k]
                   + f_3 * pc_z[k] * ksk_359[k];

        t_450[k] = f_17 * isk_360[k]
                   + f_1 * ksi0_280[k]
                   - f_2 * ksi1_280[k]
                   + f_3 * pc_x[k] * ksk_360[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pc_x, pc_y, pc_z, isk_216, isk_363, \
                         ksi0_283, ksi1_283, ksk_360, ksk_361, \
                         ksk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_18 * isk_216[k]
                   + f_3 * pc_y[k] * ksk_360[k];

        t_452[k] = f_3 * pc_z[k] * ksk_360[k];

        t_453[k] = f_17 * isk_363[k]
                   + f_12 * ksi0_283[k]
                   - f_13 * ksi1_283[k]
                   + f_3 * pc_x[k] * ksk_363[k];

        t_454[k] = f_3 * pc_z[k] * ksk_361[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_z, isk_366, ksi0_280, ksi0_286, \
                         ksi1_280, ksi1_286, ksk_362, ksk_363, \
                         ksk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_4 * ksi0_280[k]
                   - f_5 * ksi1_280[k]
                   + f_3 * pc_z[k] * ksk_362[k];

        t_456[k] = f_17 * isk_366[k]
                   + f_10 * ksi0_286[k]
                   - f_11 * ksi1_286[k]
                   + f_3 * pc_x[k] * ksk_366[k];

        t_457[k] = f_3 * pc_z[k] * ksk_363[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pc_x, pc_y, pc_z, isk_221, isk_370, \
                         ksi0_282, ksi0_290, ksi1_282, ksi1_290, ksk_365, ksk_366, \
                         ksk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_18 * isk_221[k]
                   + f_3 * pc_y[k] * ksk_365[k];

        t_459[k] = f_6 * ksi0_282[k]
                   - f_7 * ksi1_282[k]
                   + f_3 * pc_z[k] * ksk_365[k];

        t_460[k] = f_17 * isk_370[k]
                   + f_8 * ksi0_290[k]
                   - f_9 * ksi1_290[k]
                   + f_3 * pc_x[k] * ksk_370[k];

        t_461[k] = f_3 * pc_z[k] * ksk_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, isk_225, ksi0_283, ksi0_285, \
                         ksi1_283, ksi1_285, ksk_367, ksk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_4 * ksi0_283[k]
                   - f_5 * ksi1_283[k]
                   + f_3 * pc_z[k] * ksk_367[k];

        t_463[k] = f_18 * isk_225[k]
                   + f_3 * pc_y[k] * ksk_369[k];

        t_464[k] = f_8 * ksi0_285[k]
                   - f_9 * ksi1_285[k]
                   + f_3 * pc_z[k] * ksk_369[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_z, isk_375, ksi0_286, ksi0_295, \
                         ksi1_286, ksi1_295, ksk_370, ksk_371, \
                         ksk_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_17 * isk_375[k]
                   + f_6 * ksi0_295[k]
                   - f_7 * ksi1_295[k]
                   + f_3 * pc_x[k] * ksk_375[k];

        t_466[k] = f_3 * pc_z[k] * ksk_370[k];

        t_467[k] = f_4 * ksi0_286[k]
                   - f_5 * ksi1_286[k]
                   + f_3 * pc_z[k] * ksk_371[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_y, pc_z, isk_230, ksi0_287, ksi0_289, \
                         ksi1_287, ksi1_289, ksk_372, ksk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_6 * ksi0_287[k]
                   - f_7 * ksi1_287[k]
                   + f_3 * pc_z[k] * ksk_372[k];

        t_469[k] = f_18 * isk_230[k]
                   + f_3 * pc_y[k] * ksk_374[k];

        t_470[k] = f_10 * ksi0_289[k]
                   - f_11 * ksi1_289[k]
                   + f_3 * pc_z[k] * ksk_374[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_x, pc_z, isk_381, ksi0_290, ksi0_301, \
                         ksi1_290, ksi1_301, ksk_375, ksk_376, \
                         ksk_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_17 * isk_381[k]
                   + f_4 * ksi0_301[k]
                   - f_5 * ksi1_301[k]
                   + f_3 * pc_x[k] * ksk_381[k];

        t_472[k] = f_3 * pc_z[k] * ksk_375[k];

        t_473[k] = f_4 * ksi0_290[k]
                   - f_5 * ksi1_290[k]
                   + f_3 * pc_z[k] * ksk_376[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_y, pc_z, isk_236, ksi0_291, ksi0_292, \
                         ksi0_294, ksi1_291, ksi1_292, ksi1_294, ksk_377, ksk_378, \
                         ksk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_6 * ksi0_291[k]
                   - f_7 * ksi1_291[k]
                   + f_3 * pc_z[k] * ksk_377[k];

        t_475[k] = f_8 * ksi0_292[k]
                   - f_9 * ksi1_292[k]
                   + f_3 * pc_z[k] * ksk_378[k];

        t_476[k] = f_18 * isk_236[k]
                   + f_3 * pc_y[k] * ksk_380[k];

        t_477[k] = f_12 * ksi0_294[k]
                   - f_13 * ksi1_294[k]
                   + f_3 * pc_z[k] * ksk_380[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
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

    const auto *isl0_270 = buffer.data(isl0 + 270);
    const auto *isl0_273 = buffer.data(isl0 + 273);
    const auto *isl0_276 = buffer.data(isl0 + 276);
    const auto *isl0_280 = buffer.data(isl0 + 280);
    const auto *isl0_282 = buffer.data(isl0 + 282);
    const auto *isl0_285 = buffer.data(isl0 + 285);
    const auto *isl0_287 = buffer.data(isl0 + 287);
    const auto *isl0_288 = buffer.data(isl0 + 288);
    const auto *isl0_291 = buffer.data(isl0 + 291);
    const auto *isl0_293 = buffer.data(isl0 + 293);
    const auto *isl0_294 = buffer.data(isl0 + 294);
    const auto *isl0_295 = buffer.data(isl0 + 295);
    const auto *isl0_306 = buffer.data(isl0 + 306);
    const auto *isl0_405 = buffer.data(isl0 + 405);

    const auto *isk_216 = buffer.data(isk + 216);
    const auto *isk_219 = buffer.data(isk + 219);
    const auto *isk_222 = buffer.data(isk + 222);
    const auto *isk_223 = buffer.data(isk + 223);
    const auto *isk_226 = buffer.data(isk + 226);
    const auto *isk_227 = buffer.data(isk + 227);
    const auto *isk_228 = buffer.data(isk + 228);
    const auto *isk_231 = buffer.data(isk + 231);
    const auto *isk_232 = buffer.data(isk + 232);
    const auto *isk_233 = buffer.data(isk + 233);
    const auto *isk_234 = buffer.data(isk + 234);
    const auto *isk_244 = buffer.data(isk + 244);
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
    const auto *isk_280 = buffer.data(isk + 280);
    const auto *isk_282 = buffer.data(isk + 282);
    const auto *isk_283 = buffer.data(isk + 283);
    const auto *isk_284 = buffer.data(isk + 284);
    const auto *isk_285 = buffer.data(isk + 285);
    const auto *isk_286 = buffer.data(isk + 286);
    const auto *isk_287 = buffer.data(isk + 287);
    const auto *isk_288 = buffer.data(isk + 288);
    const auto *isk_290 = buffer.data(isk + 290);
    const auto *isk_293 = buffer.data(isk + 293);
    const auto *isk_297 = buffer.data(isk + 297);
    const auto *isk_302 = buffer.data(isk + 302);
    const auto *isk_308 = buffer.data(isk + 308);
    const auto *isk_316 = buffer.data(isk + 316);
    const auto *isk_318 = buffer.data(isk + 318);
    const auto *isk_319 = buffer.data(isk + 319);
    const auto *isk_320 = buffer.data(isk + 320);
    const auto *isk_321 = buffer.data(isk + 321);
    const auto *isk_322 = buffer.data(isk + 322);
    const auto *isk_323 = buffer.data(isk + 323);
    const auto *isk_324 = buffer.data(isk + 324);
    const auto *isk_388 = buffer.data(isk + 388);
    const auto *isk_390 = buffer.data(isk + 390);
    const auto *isk_391 = buffer.data(isk + 391);
    const auto *isk_392 = buffer.data(isk + 392);
    const auto *isk_393 = buffer.data(isk + 393);
    const auto *isk_394 = buffer.data(isk + 394);
    const auto *isk_395 = buffer.data(isk + 395);
    const auto *isk_401 = buffer.data(isk + 401);
    const auto *isk_405 = buffer.data(isk + 405);
    const auto *isk_410 = buffer.data(isk + 410);
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

    const auto *isl1_270 = buffer.data(isl1 + 270);
    const auto *isl1_273 = buffer.data(isl1 + 273);
    const auto *isl1_276 = buffer.data(isl1 + 276);
    const auto *isl1_280 = buffer.data(isl1 + 280);
    const auto *isl1_282 = buffer.data(isl1 + 282);
    const auto *isl1_285 = buffer.data(isl1 + 285);
    const auto *isl1_287 = buffer.data(isl1 + 287);
    const auto *isl1_288 = buffer.data(isl1 + 288);
    const auto *isl1_291 = buffer.data(isl1 + 291);
    const auto *isl1_293 = buffer.data(isl1 + 293);
    const auto *isl1_294 = buffer.data(isl1 + 294);
    const auto *isl1_295 = buffer.data(isl1 + 295);
    const auto *isl1_306 = buffer.data(isl1 + 306);
    const auto *isl1_405 = buffer.data(isl1 + 405);

    const auto *ksi0_301 = buffer.data(ksi0 + 301);
    const auto *ksi0_302 = buffer.data(ksi0 + 302);
    const auto *ksi0_303 = buffer.data(ksi0 + 303);
    const auto *ksi0_304 = buffer.data(ksi0 + 304);
    const auto *ksi0_305 = buffer.data(ksi0 + 305);
    const auto *ksi0_307 = buffer.data(ksi0 + 307);
    const auto *ksi0_313 = buffer.data(ksi0 + 313);
    const auto *ksi0_317 = buffer.data(ksi0 + 317);
    const auto *ksi0_322 = buffer.data(ksi0 + 322);
    const auto *ksi0_328 = buffer.data(ksi0 + 328);
    const auto *ksi0_331 = buffer.data(ksi0 + 331);
    const auto *ksi0_332 = buffer.data(ksi0 + 332);
    const auto *ksi0_333 = buffer.data(ksi0 + 333);
    const auto *ksi0_334 = buffer.data(ksi0 + 334);
    const auto *ksi0_335 = buffer.data(ksi0 + 335);
    const auto *ksi0_336 = buffer.data(ksi0 + 336);
    const auto *ksi0_339 = buffer.data(ksi0 + 339);
    const auto *ksi0_341 = buffer.data(ksi0 + 341);
    const auto *ksi0_342 = buffer.data(ksi0 + 342);
    const auto *ksi0_345 = buffer.data(ksi0 + 345);
    const auto *ksi0_346 = buffer.data(ksi0 + 346);
    const auto *ksi0_348 = buffer.data(ksi0 + 348);
    const auto *ksi0_350 = buffer.data(ksi0 + 350);
    const auto *ksi0_351 = buffer.data(ksi0 + 351);
    const auto *ksi0_353 = buffer.data(ksi0 + 353);
    const auto *ksi0_354 = buffer.data(ksi0 + 354);
    const auto *ksi0_356 = buffer.data(ksi0 + 356);
    const auto *ksi0_357 = buffer.data(ksi0 + 357);
    const auto *ksi0_359 = buffer.data(ksi0 + 359);
    const auto *ksi0_360 = buffer.data(ksi0 + 360);
    const auto *ksi0_361 = buffer.data(ksi0 + 361);
    const auto *ksi0_362 = buffer.data(ksi0 + 362);
    const auto *ksi0_363 = buffer.data(ksi0 + 363);

    const auto *ksi1_301 = buffer.data(ksi1 + 301);
    const auto *ksi1_302 = buffer.data(ksi1 + 302);
    const auto *ksi1_303 = buffer.data(ksi1 + 303);
    const auto *ksi1_304 = buffer.data(ksi1 + 304);
    const auto *ksi1_305 = buffer.data(ksi1 + 305);
    const auto *ksi1_307 = buffer.data(ksi1 + 307);
    const auto *ksi1_313 = buffer.data(ksi1 + 313);
    const auto *ksi1_317 = buffer.data(ksi1 + 317);
    const auto *ksi1_322 = buffer.data(ksi1 + 322);
    const auto *ksi1_328 = buffer.data(ksi1 + 328);
    const auto *ksi1_331 = buffer.data(ksi1 + 331);
    const auto *ksi1_332 = buffer.data(ksi1 + 332);
    const auto *ksi1_333 = buffer.data(ksi1 + 333);
    const auto *ksi1_334 = buffer.data(ksi1 + 334);
    const auto *ksi1_335 = buffer.data(ksi1 + 335);
    const auto *ksi1_336 = buffer.data(ksi1 + 336);
    const auto *ksi1_339 = buffer.data(ksi1 + 339);
    const auto *ksi1_341 = buffer.data(ksi1 + 341);
    const auto *ksi1_342 = buffer.data(ksi1 + 342);
    const auto *ksi1_345 = buffer.data(ksi1 + 345);
    const auto *ksi1_346 = buffer.data(ksi1 + 346);
    const auto *ksi1_348 = buffer.data(ksi1 + 348);
    const auto *ksi1_350 = buffer.data(ksi1 + 350);
    const auto *ksi1_351 = buffer.data(ksi1 + 351);
    const auto *ksi1_353 = buffer.data(ksi1 + 353);
    const auto *ksi1_354 = buffer.data(ksi1 + 354);
    const auto *ksi1_356 = buffer.data(ksi1 + 356);
    const auto *ksi1_357 = buffer.data(ksi1 + 357);
    const auto *ksi1_359 = buffer.data(ksi1 + 359);
    const auto *ksi1_360 = buffer.data(ksi1 + 360);
    const auto *ksi1_361 = buffer.data(ksi1 + 361);
    const auto *ksi1_362 = buffer.data(ksi1 + 362);
    const auto *ksi1_363 = buffer.data(ksi1 + 363);

    const auto *ksk_381 = buffer.data(ksk + 381);
    const auto *ksk_388 = buffer.data(ksk + 388);
    const auto *ksk_389 = buffer.data(ksk + 389);
    const auto *ksk_390 = buffer.data(ksk + 390);
    const auto *ksk_391 = buffer.data(ksk + 391);
    const auto *ksk_392 = buffer.data(ksk + 392);
    const auto *ksk_393 = buffer.data(ksk + 393);
    const auto *ksk_394 = buffer.data(ksk + 394);
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
    const auto *ksk_434 = buffer.data(ksk + 434);
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
    const auto *ksk_468 = buffer.data(ksk + 468);

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, pc_x, pc_z, isk_388, isk_390, \
                         isk_391, isk_392, ksk_381, ksk_388, ksk_390, ksk_391, \
                         ksk_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_17 * isk_388[k]
                   + f_3 * pc_x[k] * ksk_388[k];

        t_479[k] = f_3 * pc_z[k] * ksk_381[k];

        t_480[k] = f_17 * isk_390[k]
                   + f_3 * pc_x[k] * ksk_390[k];

        t_481[k] = f_17 * isk_391[k]
                   + f_3 * pc_x[k] * ksk_391[k];

        t_482[k] = f_17 * isk_392[k]
                   + f_3 * pc_x[k] * ksk_392[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pc_x, pc_y, isk_244, isk_393, isk_394, \
                         isk_395, ksi0_301, ksi1_301, ksk_388, ksk_393, ksk_394, \
                         ksk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_17 * isk_393[k]
                   + f_3 * pc_x[k] * ksk_393[k];

        t_484[k] = f_17 * isk_394[k]
                   + f_3 * pc_x[k] * ksk_394[k];

        t_485[k] = f_17 * isk_395[k]
                   + f_3 * pc_x[k] * ksk_395[k];

        t_486[k] = f_18 * isk_244[k]
                   + f_1 * ksi0_301[k]
                   - f_2 * ksi1_301[k]
                   + f_3 * pc_y[k] * ksk_388[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pc_z, ksi0_301, ksi0_302, ksi0_303, \
                         ksi1_301, ksi1_302, ksi1_303, ksk_388, ksk_389, ksk_390, \
                         ksk_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_3 * pc_z[k] * ksk_388[k];

        t_488[k] = f_4 * ksi0_301[k]
                   - f_5 * ksi1_301[k]
                   + f_3 * pc_z[k] * ksk_389[k];

        t_489[k] = f_6 * ksi0_302[k]
                   - f_7 * ksi1_302[k]
                   + f_3 * pc_z[k] * ksk_390[k];

        t_490[k] = f_8 * ksi0_303[k]
                   - f_9 * ksi1_303[k]
                   + f_3 * pc_z[k] * ksk_391[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_y, pc_z, isk_251, ksi0_304, ksi0_305, \
                         ksi0_307, ksi1_304, ksi1_305, ksi1_307, ksk_392, ksk_393, \
                         ksk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_10 * ksi0_304[k]
                   - f_11 * ksi1_304[k]
                   + f_3 * pc_z[k] * ksk_392[k];

        t_492[k] = f_12 * ksi0_305[k]
                   - f_13 * ksi1_305[k]
                   + f_3 * pc_z[k] * ksk_393[k];

        t_493[k] = f_18 * isk_251[k]
                   + f_3 * pc_y[k] * ksk_395[k];

        t_494[k] = f_1 * ksi0_307[k]
                   - f_2 * ksi1_307[k]
                   + f_3 * pc_z[k] * ksk_395[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pa_z, pc_y, pc_z, isl0_270, isl0_273, \
                         isk_216, isk_252, isl1_270, isl1_273, \
                         ksk_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * isl0_270[k]
                   - f_14 * pc_z[k] * isl1_270[k];

        t_496[k] = f_17 * isk_252[k]
                   + f_3 * pc_y[k] * ksk_396[k];

        t_497[k] = f_15 * isk_216[k]
                   + f_3 * pc_z[k] * ksk_396[k];

        t_498[k] = pa_z[k] * isl0_273[k]
                   - f_14 * pc_z[k] * isl1_273[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pa_z, pc_x, pc_y, pc_z, isl0_276, isk_254, \
                         isk_401, isl1_276, ksi0_313, ksi1_313, ksk_398, \
                         ksk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_17 * isk_254[k]
                   + f_3 * pc_y[k] * ksk_398[k];

        t_500[k] = f_17 * isk_401[k]
                   + f_12 * ksi0_313[k]
                   - f_13 * ksi1_313[k]
                   + f_3 * pc_x[k] * ksk_401[k];

        t_501[k] = pa_z[k] * isl0_276[k]
                   - f_14 * pc_z[k] * isl1_276[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, isk_219, isk_257, isk_405, \
                         ksi0_317, ksi1_317, ksk_399, ksk_401, \
                         ksk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_15 * isk_219[k]
                   + f_3 * pc_z[k] * ksk_399[k];

        t_503[k] = f_17 * isk_257[k]
                   + f_3 * pc_y[k] * ksk_401[k];

        t_504[k] = f_17 * isk_405[k]
                   + f_10 * ksi0_317[k]
                   - f_11 * ksi1_317[k]
                   + f_3 * pc_x[k] * ksk_405[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pa_z, pc_y, pc_z, isl0_280, isl0_282, \
                         isk_222, isk_223, isk_261, isl1_280, isl1_282, ksk_402, \
                         ksk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pa_z[k] * isl0_280[k]
                   - f_14 * pc_z[k] * isl1_280[k];

        t_506[k] = f_15 * isk_222[k]
                   + f_3 * pc_z[k] * ksk_402[k];

        t_507[k] = pa_z[k] * isl0_282[k]
                   + f_16 * isk_223[k]
                   - f_14 * pc_z[k] * isl1_282[k];

        t_508[k] = f_17 * isk_261[k]
                   + f_3 * pc_y[k] * ksk_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_z, pc_x, pc_z, isl0_285, isk_226, isk_410, \
                         isl1_285, ksi0_322, ksi1_322, ksk_406, \
                         ksk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_17 * isk_410[k]
                   + f_8 * ksi0_322[k]
                   - f_9 * ksi1_322[k]
                   + f_3 * pc_x[k] * ksk_410[k];

        t_510[k] = pa_z[k] * isl0_285[k]
                   - f_14 * pc_z[k] * isl1_285[k];

        t_511[k] = f_15 * isk_226[k]
                   + f_3 * pc_z[k] * ksk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pa_z, pc_y, pc_z, isl0_287, isl0_288, isk_227, \
                         isk_228, isk_266, isl1_287, isl1_288, \
                         ksk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pa_z[k] * isl0_287[k]
                   + f_16 * isk_227[k]
                   - f_14 * pc_z[k] * isl1_287[k];

        t_513[k] = pa_z[k] * isl0_288[k]
                   + f_17 * isk_228[k]
                   - f_14 * pc_z[k] * isl1_288[k];

        t_514[k] = f_17 * isk_266[k]
                   + f_3 * pc_y[k] * ksk_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pa_z, pc_x, pc_z, isl0_291, isk_231, isk_416, \
                         isl1_291, ksi0_328, ksi1_328, ksk_411, \
                         ksk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_17 * isk_416[k]
                   + f_6 * ksi0_328[k]
                   - f_7 * ksi1_328[k]
                   + f_3 * pc_x[k] * ksk_416[k];

        t_516[k] = pa_z[k] * isl0_291[k]
                   - f_14 * pc_z[k] * isl1_291[k];

        t_517[k] = f_15 * isk_231[k]
                   + f_3 * pc_z[k] * ksk_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pa_z, pc_z, isl0_293, isl0_294, isl0_295, \
                         isk_232, isk_233, isk_234, isl1_293, isl1_294, \
                         isl1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pa_z[k] * isl0_293[k]
                   + f_16 * isk_232[k]
                   - f_14 * pc_z[k] * isl1_293[k];

        t_519[k] = pa_z[k] * isl0_294[k]
                   + f_17 * isk_233[k]
                   - f_14 * pc_z[k] * isl1_294[k];

        t_520[k] = pa_z[k] * isl0_295[k]
                   + f_18 * isk_234[k]
                   - f_14 * pc_z[k] * isl1_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, isk_272, isk_423, isk_424, \
                         isk_425, ksi0_335, ksi1_335, ksk_416, ksk_423, ksk_424, \
                         ksk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * isk_272[k]
                   + f_3 * pc_y[k] * ksk_416[k];

        t_522[k] = f_17 * isk_423[k]
                   + f_4 * ksi0_335[k]
                   - f_5 * ksi1_335[k]
                   + f_3 * pc_x[k] * ksk_423[k];

        t_523[k] = f_17 * isk_424[k]
                   + f_3 * pc_x[k] * ksk_424[k];

        t_524[k] = f_17 * isk_425[k]
                   + f_3 * pc_x[k] * ksk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, isk_426, isk_427, isk_428, \
                         isk_429, isk_430, ksk_426, ksk_427, ksk_428, ksk_429, \
                         ksk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_17 * isk_426[k]
                   + f_3 * pc_x[k] * ksk_426[k];

        t_526[k] = f_17 * isk_427[k]
                   + f_3 * pc_x[k] * ksk_427[k];

        t_527[k] = f_17 * isk_428[k]
                   + f_3 * pc_x[k] * ksk_428[k];

        t_528[k] = f_17 * isk_429[k]
                   + f_3 * pc_x[k] * ksk_429[k];

        t_529[k] = f_17 * isk_430[k]
                   + f_3 * pc_x[k] * ksk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pa_z, pc_x, pc_z, isl0_306, isk_244, isk_431, \
                         isl1_306, ksk_424, ksk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_17 * isk_431[k]
                   + f_3 * pc_x[k] * ksk_431[k];

        t_531[k] = pa_z[k] * isl0_306[k]
                   - f_14 * pc_z[k] * isl1_306[k];

        t_532[k] = f_15 * isk_244[k]
                   + f_3 * pc_z[k] * ksk_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pc_y, isk_282, isk_283, isk_284, ksi0_331, \
                         ksi0_332, ksi0_333, ksi1_331, ksi1_332, ksi1_333, ksk_426, ksk_427, \
                         ksk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * isk_282[k]
                   + f_12 * ksi0_331[k]
                   - f_13 * ksi1_331[k]
                   + f_3 * pc_y[k] * ksk_426[k];

        t_534[k] = f_17 * isk_283[k]
                   + f_10 * ksi0_332[k]
                   - f_11 * ksi1_332[k]
                   + f_3 * pc_y[k] * ksk_427[k];

        t_535[k] = f_17 * isk_284[k]
                   + f_8 * ksi0_333[k]
                   - f_9 * ksi1_333[k]
                   + f_3 * pc_y[k] * ksk_428[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pc_y, isk_285, isk_286, isk_287, ksi0_334, \
                         ksi0_335, ksi1_334, ksi1_335, ksk_429, ksk_430, \
                         ksk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_17 * isk_285[k]
                   + f_6 * ksi0_334[k]
                   - f_7 * ksi1_334[k]
                   + f_3 * pc_y[k] * ksk_429[k];

        t_537[k] = f_17 * isk_286[k]
                   + f_4 * ksi0_335[k]
                   - f_5 * ksi1_335[k]
                   + f_3 * pc_y[k] * ksk_430[k];

        t_538[k] = f_17 * isk_287[k]
                   + f_3 * pc_y[k] * ksk_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pc_x, pc_y, pc_z, isk_251, isk_288, isk_432, \
                         ksi0_335, ksi0_336, ksi1_335, ksi1_336, ksk_431, \
                         ksk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * isk_251[k]
                   + f_1 * ksi0_335[k]
                   - f_2 * ksi1_335[k]
                   + f_3 * pc_z[k] * ksk_431[k];

        t_540[k] = f_17 * isk_432[k]
                   + f_1 * ksi0_336[k]
                   - f_2 * ksi1_336[k]
                   + f_3 * pc_x[k] * ksk_432[k];

        t_541[k] = f_16 * isk_288[k]
                   + f_3 * pc_y[k] * ksk_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, pc_y, pc_z, isk_252, isk_290, isk_435, \
                         ksi0_339, ksi1_339, ksk_432, ksk_434, \
                         ksk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * isk_252[k]
                   + f_3 * pc_z[k] * ksk_432[k];

        t_543[k] = f_17 * isk_435[k]
                   + f_12 * ksi0_339[k]
                   - f_13 * ksi1_339[k]
                   + f_3 * pc_x[k] * ksk_435[k];

        t_544[k] = f_16 * isk_290[k]
                   + f_3 * pc_y[k] * ksk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, isk_255, isk_437, isk_438, ksi0_341, \
                         ksi0_342, ksi1_341, ksi1_342, ksk_435, ksk_437, \
                         ksk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_17 * isk_437[k]
                   + f_12 * ksi0_341[k]
                   - f_13 * ksi1_341[k]
                   + f_3 * pc_x[k] * ksk_437[k];

        t_546[k] = f_17 * isk_438[k]
                   + f_10 * ksi0_342[k]
                   - f_11 * ksi1_342[k]
                   + f_3 * pc_x[k] * ksk_438[k];

        t_547[k] = f_16 * isk_255[k]
                   + f_3 * pc_z[k] * ksk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pc_x, pc_y, isk_293, isk_441, isk_442, ksi0_345, \
                         ksi0_346, ksi1_345, ksi1_346, ksk_437, ksk_441, \
                         ksk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * isk_293[k]
                   + f_3 * pc_y[k] * ksk_437[k];

        t_549[k] = f_17 * isk_441[k]
                   + f_10 * ksi0_345[k]
                   - f_11 * ksi1_345[k]
                   + f_3 * pc_x[k] * ksk_441[k];

        t_550[k] = f_17 * isk_442[k]
                   + f_8 * ksi0_346[k]
                   - f_9 * ksi1_346[k]
                   + f_3 * pc_x[k] * ksk_442[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pc_x, pc_y, pc_z, isk_258, isk_297, isk_444, \
                         ksi0_348, ksi1_348, ksk_438, ksk_441, \
                         ksk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * isk_258[k]
                   + f_3 * pc_z[k] * ksk_438[k];

        t_552[k] = f_17 * isk_444[k]
                   + f_8 * ksi0_348[k]
                   - f_9 * ksi1_348[k]
                   + f_3 * pc_x[k] * ksk_444[k];

        t_553[k] = f_16 * isk_297[k]
                   + f_3 * pc_y[k] * ksk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_x, pc_z, isk_262, isk_446, isk_447, ksi0_350, \
                         ksi0_351, ksi1_350, ksi1_351, ksk_442, ksk_446, \
                         ksk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_17 * isk_446[k]
                   + f_8 * ksi0_350[k]
                   - f_9 * ksi1_350[k]
                   + f_3 * pc_x[k] * ksk_446[k];

        t_555[k] = f_17 * isk_447[k]
                   + f_6 * ksi0_351[k]
                   - f_7 * ksi1_351[k]
                   + f_3 * pc_x[k] * ksk_447[k];

        t_556[k] = f_16 * isk_262[k]
                   + f_3 * pc_z[k] * ksk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, pc_y, isk_302, isk_449, isk_450, ksi0_353, \
                         ksi0_354, ksi1_353, ksi1_354, ksk_446, ksk_449, \
                         ksk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_17 * isk_449[k]
                   + f_6 * ksi0_353[k]
                   - f_7 * ksi1_353[k]
                   + f_3 * pc_x[k] * ksk_449[k];

        t_558[k] = f_17 * isk_450[k]
                   + f_6 * ksi0_354[k]
                   - f_7 * ksi1_354[k]
                   + f_3 * pc_x[k] * ksk_450[k];

        t_559[k] = f_16 * isk_302[k]
                   + f_3 * pc_y[k] * ksk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pc_x, pc_z, isk_267, isk_452, isk_453, ksi0_356, \
                         ksi0_357, ksi1_356, ksi1_357, ksk_447, ksk_452, \
                         ksk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_17 * isk_452[k]
                   + f_6 * ksi0_356[k]
                   - f_7 * ksi1_356[k]
                   + f_3 * pc_x[k] * ksk_452[k];

        t_561[k] = f_17 * isk_453[k]
                   + f_4 * ksi0_357[k]
                   - f_5 * ksi1_357[k]
                   + f_3 * pc_x[k] * ksk_453[k];

        t_562[k] = f_16 * isk_267[k]
                   + f_3 * pc_z[k] * ksk_447[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, isk_455, isk_456, isk_457, ksi0_359, \
                         ksi0_360, ksi0_361, ksi1_359, ksi1_360, ksi1_361, ksk_455, ksk_456, \
                         ksk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_17 * isk_455[k]
                   + f_4 * ksi0_359[k]
                   - f_5 * ksi1_359[k]
                   + f_3 * pc_x[k] * ksk_455[k];

        t_564[k] = f_17 * isk_456[k]
                   + f_4 * ksi0_360[k]
                   - f_5 * ksi1_360[k]
                   + f_3 * pc_x[k] * ksk_456[k];

        t_565[k] = f_17 * isk_457[k]
                   + f_4 * ksi0_361[k]
                   - f_5 * ksi1_361[k]
                   + f_3 * pc_x[k] * ksk_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pc_x, pc_y, isk_308, isk_459, isk_460, \
                         isk_461, ksi0_363, ksi1_363, ksk_452, ksk_459, ksk_460, \
                         ksk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * isk_308[k]
                   + f_3 * pc_y[k] * ksk_452[k];

        t_567[k] = f_17 * isk_459[k]
                   + f_4 * ksi0_363[k]
                   - f_5 * ksi1_363[k]
                   + f_3 * pc_x[k] * ksk_459[k];

        t_568[k] = f_17 * isk_460[k]
                   + f_3 * pc_x[k] * ksk_460[k];

        t_569[k] = f_17 * isk_461[k]
                   + f_3 * pc_x[k] * ksk_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, isk_462, isk_463, isk_464, \
                         isk_465, isk_466, ksk_462, ksk_463, ksk_464, ksk_465, \
                         ksk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_17 * isk_462[k]
                   + f_3 * pc_x[k] * ksk_462[k];

        t_571[k] = f_17 * isk_463[k]
                   + f_3 * pc_x[k] * ksk_463[k];

        t_572[k] = f_17 * isk_464[k]
                   + f_3 * pc_x[k] * ksk_464[k];

        t_573[k] = f_17 * isk_465[k]
                   + f_3 * pc_x[k] * ksk_465[k];

        t_574[k] = f_17 * isk_466[k]
                   + f_3 * pc_x[k] * ksk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, pc_z, isk_280, isk_316, isk_467, \
                         ksi0_357, ksi1_357, ksk_460, ksk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_17 * isk_467[k]
                   + f_3 * pc_x[k] * ksk_467[k];

        t_576[k] = f_16 * isk_316[k]
                   + f_1 * ksi0_357[k]
                   - f_2 * ksi1_357[k]
                   + f_3 * pc_y[k] * ksk_460[k];

        t_577[k] = f_16 * isk_280[k]
                   + f_3 * pc_z[k] * ksk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, isk_318, isk_319, isk_320, ksi0_359, \
                         ksi0_360, ksi0_361, ksi1_359, ksi1_360, ksi1_361, ksk_462, ksk_463, \
                         ksk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * isk_318[k]
                   + f_12 * ksi0_359[k]
                   - f_13 * ksi1_359[k]
                   + f_3 * pc_y[k] * ksk_462[k];

        t_579[k] = f_16 * isk_319[k]
                   + f_10 * ksi0_360[k]
                   - f_11 * ksi1_360[k]
                   + f_3 * pc_y[k] * ksk_463[k];

        t_580[k] = f_16 * isk_320[k]
                   + f_8 * ksi0_361[k]
                   - f_9 * ksi1_361[k]
                   + f_3 * pc_y[k] * ksk_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, isk_321, isk_322, isk_323, ksi0_362, \
                         ksi0_363, ksi1_362, ksi1_363, ksk_465, ksk_466, \
                         ksk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * isk_321[k]
                   + f_6 * ksi0_362[k]
                   - f_7 * ksi1_362[k]
                   + f_3 * pc_y[k] * ksk_465[k];

        t_582[k] = f_16 * isk_322[k]
                   + f_4 * ksi0_363[k]
                   - f_5 * ksi1_363[k]
                   + f_3 * pc_y[k] * ksk_466[k];

        t_583[k] = f_16 * isk_323[k]
                   + f_3 * pc_y[k] * ksk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_y, pc_y, pc_z, isl0_405, isk_287, \
                         isk_288, isk_324, isl1_405, ksi0_363, ksi1_363, ksk_467, \
                         ksk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * isk_287[k]
                   + f_1 * ksi0_363[k]
                   - f_2 * ksi1_363[k]
                   + f_3 * pc_z[k] * ksk_467[k];

        t_585[k] = pa_y[k] * isl0_405[k]
                   - f_14 * pc_y[k] * isl1_405[k];

        t_586[k] = f_15 * isk_324[k]
                   + f_3 * pc_y[k] * ksk_468[k];

        t_587[k] = f_17 * isk_288[k]
                   + f_3 * pc_z[k] * ksk_468[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
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
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);

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

    const auto *isl0_408 = buffer.data(isl0 + 408);
    const auto *isl0_410 = buffer.data(isl0 + 410);
    const auto *isl0_411 = buffer.data(isl0 + 411);
    const auto *isl0_414 = buffer.data(isl0 + 414);
    const auto *isl0_415 = buffer.data(isl0 + 415);
    const auto *isl0_417 = buffer.data(isl0 + 417);
    const auto *isl0_419 = buffer.data(isl0 + 419);
    const auto *isl0_420 = buffer.data(isl0 + 420);
    const auto *isl0_422 = buffer.data(isl0 + 422);
    const auto *isl0_423 = buffer.data(isl0 + 423);
    const auto *isl0_425 = buffer.data(isl0 + 425);
    const auto *isl0_426 = buffer.data(isl0 + 426);
    const auto *isl0_428 = buffer.data(isl0 + 428);
    const auto *isl0_429 = buffer.data(isl0 + 429);
    const auto *isl0_430 = buffer.data(isl0 + 430);
    const auto *isl0_432 = buffer.data(isl0 + 432);
    const auto *isl0_449 = buffer.data(isl0 + 449);

    const auto *isk_291 = buffer.data(isk + 291);
    const auto *isk_294 = buffer.data(isk + 294);
    const auto *isk_298 = buffer.data(isk + 298);
    const auto *isk_303 = buffer.data(isk + 303);
    const auto *isk_316 = buffer.data(isk + 316);
    const auto *isk_324 = buffer.data(isk + 324);
    const auto *isk_325 = buffer.data(isk + 325);
    const auto *isk_326 = buffer.data(isk + 326);
    const auto *isk_327 = buffer.data(isk + 327);
    const auto *isk_329 = buffer.data(isk + 329);
    const auto *isk_330 = buffer.data(isk + 330);
    const auto *isk_332 = buffer.data(isk + 332);
    const auto *isk_333 = buffer.data(isk + 333);
    const auto *isk_334 = buffer.data(isk + 334);
    const auto *isk_336 = buffer.data(isk + 336);
    const auto *isk_337 = buffer.data(isk + 337);
    const auto *isk_338 = buffer.data(isk + 338);
    const auto *isk_339 = buffer.data(isk + 339);
    const auto *isk_341 = buffer.data(isk + 341);
    const auto *isk_342 = buffer.data(isk + 342);
    const auto *isk_343 = buffer.data(isk + 343);
    const auto *isk_344 = buffer.data(isk + 344);
    const auto *isk_352 = buffer.data(isk + 352);
    const auto *isk_354 = buffer.data(isk + 354);
    const auto *isk_355 = buffer.data(isk + 355);
    const auto *isk_356 = buffer.data(isk + 356);
    const auto *isk_357 = buffer.data(isk + 357);
    const auto *isk_358 = buffer.data(isk + 358);
    const auto *isk_359 = buffer.data(isk + 359);
    const auto *isk_360 = buffer.data(isk + 360);
    const auto *isk_365 = buffer.data(isk + 365);
    const auto *isk_369 = buffer.data(isk + 369);
    const auto *isk_374 = buffer.data(isk + 374);
    const auto *isk_380 = buffer.data(isk + 380);
    const auto *isk_496 = buffer.data(isk + 496);
    const auto *isk_497 = buffer.data(isk + 497);
    const auto *isk_498 = buffer.data(isk + 498);
    const auto *isk_499 = buffer.data(isk + 499);
    const auto *isk_500 = buffer.data(isk + 500);
    const auto *isk_501 = buffer.data(isk + 501);
    const auto *isk_502 = buffer.data(isk + 502);
    const auto *isk_503 = buffer.data(isk + 503);
    const auto *isk_504 = buffer.data(isk + 504);
    const auto *isk_509 = buffer.data(isk + 509);
    const auto *isk_513 = buffer.data(isk + 513);
    const auto *isk_518 = buffer.data(isk + 518);
    const auto *isk_524 = buffer.data(isk + 524);
    const auto *isk_531 = buffer.data(isk + 531);
    const auto *isk_532 = buffer.data(isk + 532);
    const auto *isk_533 = buffer.data(isk + 533);
    const auto *isk_534 = buffer.data(isk + 534);
    const auto *isk_535 = buffer.data(isk + 535);
    const auto *isk_536 = buffer.data(isk + 536);
    const auto *isk_537 = buffer.data(isk + 537);
    const auto *isk_539 = buffer.data(isk + 539);
    const auto *isk_540 = buffer.data(isk + 540);
    const auto *isk_543 = buffer.data(isk + 543);
    const auto *isk_546 = buffer.data(isk + 546);
    const auto *isk_550 = buffer.data(isk + 550);
    const auto *isk_555 = buffer.data(isk + 555);
    const auto *isk_561 = buffer.data(isk + 561);

    const auto *isl1_408 = buffer.data(isl1 + 408);
    const auto *isl1_410 = buffer.data(isl1 + 410);
    const auto *isl1_411 = buffer.data(isl1 + 411);
    const auto *isl1_414 = buffer.data(isl1 + 414);
    const auto *isl1_415 = buffer.data(isl1 + 415);
    const auto *isl1_417 = buffer.data(isl1 + 417);
    const auto *isl1_419 = buffer.data(isl1 + 419);
    const auto *isl1_420 = buffer.data(isl1 + 420);
    const auto *isl1_422 = buffer.data(isl1 + 422);
    const auto *isl1_423 = buffer.data(isl1 + 423);
    const auto *isl1_425 = buffer.data(isl1 + 425);
    const auto *isl1_426 = buffer.data(isl1 + 426);
    const auto *isl1_428 = buffer.data(isl1 + 428);
    const auto *isl1_429 = buffer.data(isl1 + 429);
    const auto *isl1_430 = buffer.data(isl1 + 430);
    const auto *isl1_432 = buffer.data(isl1 + 432);
    const auto *isl1_449 = buffer.data(isl1 + 449);

    const auto *ksi0_385 = buffer.data(ksi0 + 385);
    const auto *ksi0_387 = buffer.data(ksi0 + 387);
    const auto *ksi0_388 = buffer.data(ksi0 + 388);
    const auto *ksi0_389 = buffer.data(ksi0 + 389);
    const auto *ksi0_390 = buffer.data(ksi0 + 390);
    const auto *ksi0_391 = buffer.data(ksi0 + 391);
    const auto *ksi0_392 = buffer.data(ksi0 + 392);
    const auto *ksi0_393 = buffer.data(ksi0 + 393);
    const auto *ksi0_394 = buffer.data(ksi0 + 394);
    const auto *ksi0_395 = buffer.data(ksi0 + 395);
    const auto *ksi0_396 = buffer.data(ksi0 + 396);
    const auto *ksi0_397 = buffer.data(ksi0 + 397);
    const auto *ksi0_398 = buffer.data(ksi0 + 398);
    const auto *ksi0_399 = buffer.data(ksi0 + 399);
    const auto *ksi0_400 = buffer.data(ksi0 + 400);
    const auto *ksi0_401 = buffer.data(ksi0 + 401);
    const auto *ksi0_402 = buffer.data(ksi0 + 402);
    const auto *ksi0_403 = buffer.data(ksi0 + 403);
    const auto *ksi0_404 = buffer.data(ksi0 + 404);
    const auto *ksi0_405 = buffer.data(ksi0 + 405);
    const auto *ksi0_406 = buffer.data(ksi0 + 406);
    const auto *ksi0_412 = buffer.data(ksi0 + 412);
    const auto *ksi0_413 = buffer.data(ksi0 + 413);
    const auto *ksi0_414 = buffer.data(ksi0 + 414);
    const auto *ksi0_415 = buffer.data(ksi0 + 415);
    const auto *ksi0_416 = buffer.data(ksi0 + 416);
    const auto *ksi0_417 = buffer.data(ksi0 + 417);
    const auto *ksi0_418 = buffer.data(ksi0 + 418);
    const auto *ksi0_419 = buffer.data(ksi0 + 419);
    const auto *ksi0_420 = buffer.data(ksi0 + 420);
    const auto *ksi0_422 = buffer.data(ksi0 + 422);
    const auto *ksi0_423 = buffer.data(ksi0 + 423);
    const auto *ksi0_425 = buffer.data(ksi0 + 425);
    const auto *ksi0_426 = buffer.data(ksi0 + 426);
    const auto *ksi0_427 = buffer.data(ksi0 + 427);
    const auto *ksi0_429 = buffer.data(ksi0 + 429);
    const auto *ksi0_430 = buffer.data(ksi0 + 430);
    const auto *ksi0_431 = buffer.data(ksi0 + 431);
    const auto *ksi0_432 = buffer.data(ksi0 + 432);
    const auto *ksi0_434 = buffer.data(ksi0 + 434);
    const auto *ksi0_435 = buffer.data(ksi0 + 435);
    const auto *ksi0_441 = buffer.data(ksi0 + 441);

    const auto *ksi1_385 = buffer.data(ksi1 + 385);
    const auto *ksi1_387 = buffer.data(ksi1 + 387);
    const auto *ksi1_388 = buffer.data(ksi1 + 388);
    const auto *ksi1_389 = buffer.data(ksi1 + 389);
    const auto *ksi1_390 = buffer.data(ksi1 + 390);
    const auto *ksi1_391 = buffer.data(ksi1 + 391);
    const auto *ksi1_392 = buffer.data(ksi1 + 392);
    const auto *ksi1_393 = buffer.data(ksi1 + 393);
    const auto *ksi1_394 = buffer.data(ksi1 + 394);
    const auto *ksi1_395 = buffer.data(ksi1 + 395);
    const auto *ksi1_396 = buffer.data(ksi1 + 396);
    const auto *ksi1_397 = buffer.data(ksi1 + 397);
    const auto *ksi1_398 = buffer.data(ksi1 + 398);
    const auto *ksi1_399 = buffer.data(ksi1 + 399);
    const auto *ksi1_400 = buffer.data(ksi1 + 400);
    const auto *ksi1_401 = buffer.data(ksi1 + 401);
    const auto *ksi1_402 = buffer.data(ksi1 + 402);
    const auto *ksi1_403 = buffer.data(ksi1 + 403);
    const auto *ksi1_404 = buffer.data(ksi1 + 404);
    const auto *ksi1_405 = buffer.data(ksi1 + 405);
    const auto *ksi1_406 = buffer.data(ksi1 + 406);
    const auto *ksi1_412 = buffer.data(ksi1 + 412);
    const auto *ksi1_413 = buffer.data(ksi1 + 413);
    const auto *ksi1_414 = buffer.data(ksi1 + 414);
    const auto *ksi1_415 = buffer.data(ksi1 + 415);
    const auto *ksi1_416 = buffer.data(ksi1 + 416);
    const auto *ksi1_417 = buffer.data(ksi1 + 417);
    const auto *ksi1_418 = buffer.data(ksi1 + 418);
    const auto *ksi1_419 = buffer.data(ksi1 + 419);
    const auto *ksi1_420 = buffer.data(ksi1 + 420);
    const auto *ksi1_422 = buffer.data(ksi1 + 422);
    const auto *ksi1_423 = buffer.data(ksi1 + 423);
    const auto *ksi1_425 = buffer.data(ksi1 + 425);
    const auto *ksi1_426 = buffer.data(ksi1 + 426);
    const auto *ksi1_427 = buffer.data(ksi1 + 427);
    const auto *ksi1_429 = buffer.data(ksi1 + 429);
    const auto *ksi1_430 = buffer.data(ksi1 + 430);
    const auto *ksi1_431 = buffer.data(ksi1 + 431);
    const auto *ksi1_432 = buffer.data(ksi1 + 432);
    const auto *ksi1_434 = buffer.data(ksi1 + 434);
    const auto *ksi1_435 = buffer.data(ksi1 + 435);
    const auto *ksi1_441 = buffer.data(ksi1 + 441);

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
    const auto *ksk_497 = buffer.data(ksk + 497);
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
    const auto *ksk_508 = buffer.data(ksk + 508);
    const auto *ksk_509 = buffer.data(ksk + 509);
    const auto *ksk_510 = buffer.data(ksk + 510);
    const auto *ksk_511 = buffer.data(ksk + 511);
    const auto *ksk_512 = buffer.data(ksk + 512);
    const auto *ksk_513 = buffer.data(ksk + 513);
    const auto *ksk_514 = buffer.data(ksk + 514);
    const auto *ksk_515 = buffer.data(ksk + 515);
    const auto *ksk_516 = buffer.data(ksk + 516);
    const auto *ksk_517 = buffer.data(ksk + 517);
    const auto *ksk_518 = buffer.data(ksk + 518);
    const auto *ksk_519 = buffer.data(ksk + 519);
    const auto *ksk_520 = buffer.data(ksk + 520);
    const auto *ksk_521 = buffer.data(ksk + 521);
    const auto *ksk_522 = buffer.data(ksk + 522);
    const auto *ksk_523 = buffer.data(ksk + 523);
    const auto *ksk_524 = buffer.data(ksk + 524);
    const auto *ksk_531 = buffer.data(ksk + 531);
    const auto *ksk_532 = buffer.data(ksk + 532);
    const auto *ksk_533 = buffer.data(ksk + 533);
    const auto *ksk_534 = buffer.data(ksk + 534);
    const auto *ksk_535 = buffer.data(ksk + 535);
    const auto *ksk_536 = buffer.data(ksk + 536);
    const auto *ksk_537 = buffer.data(ksk + 537);
    const auto *ksk_538 = buffer.data(ksk + 538);
    const auto *ksk_539 = buffer.data(ksk + 539);
    const auto *ksk_540 = buffer.data(ksk + 540);
    const auto *ksk_541 = buffer.data(ksk + 541);
    const auto *ksk_542 = buffer.data(ksk + 542);
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
    const auto *ksk_561 = buffer.data(ksk + 561);

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_y, pc_y, isl0_408, isl0_410, isl0_411, \
                         isk_325, isk_326, isk_327, isl1_408, isl1_410, isl1_411, \
                         ksk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_y[k] * isl0_408[k]
                   + f_16 * isk_325[k]
                   - f_14 * pc_y[k] * isl1_408[k];

        t_589[k] = f_15 * isk_326[k]
                   + f_3 * pc_y[k] * ksk_470[k];

        t_590[k] = pa_y[k] * isl0_410[k]
                   - f_14 * pc_y[k] * isl1_410[k];

        t_591[k] = pa_y[k] * isl0_411[k]
                   + f_17 * isk_327[k]
                   - f_14 * pc_y[k] * isl1_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_y, pc_y, pc_z, isl0_414, isl0_415, \
                         isk_291, isk_329, isk_330, isl1_414, isl1_415, ksk_471, \
                         ksk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * isk_291[k]
                   + f_3 * pc_z[k] * ksk_471[k];

        t_593[k] = f_15 * isk_329[k]
                   + f_3 * pc_y[k] * ksk_473[k];

        t_594[k] = pa_y[k] * isl0_414[k]
                   - f_14 * pc_y[k] * isl1_414[k];

        t_595[k] = pa_y[k] * isl0_415[k]
                   + f_18 * isk_330[k]
                   - f_14 * pc_y[k] * isl1_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pa_y, pc_y, pc_z, isl0_417, isl0_419, \
                         isk_294, isk_332, isk_333, isl1_417, isl1_419, ksk_474, \
                         ksk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * isk_294[k]
                   + f_3 * pc_z[k] * ksk_474[k];

        t_597[k] = pa_y[k] * isl0_417[k]
                   + f_16 * isk_332[k]
                   - f_14 * pc_y[k] * isl1_417[k];

        t_598[k] = f_15 * isk_333[k]
                   + f_3 * pc_y[k] * ksk_477[k];

        t_599[k] = pa_y[k] * isl0_419[k]
                   - f_14 * pc_y[k] * isl1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pa_y, pc_y, pc_z, isl0_420, isl0_422, isk_298, \
                         isk_334, isk_336, isl1_420, isl1_422, \
                         ksk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pa_y[k] * isl0_420[k]
                   + f_19 * isk_334[k]
                   - f_14 * pc_y[k] * isl1_420[k];

        t_601[k] = f_17 * isk_298[k]
                   + f_3 * pc_z[k] * ksk_478[k];

        t_602[k] = pa_y[k] * isl0_422[k]
                   + f_17 * isk_336[k]
                   - f_14 * pc_y[k] * isl1_422[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pa_y, pc_y, isl0_423, isl0_425, isl0_426, \
                         isk_337, isk_338, isk_339, isl1_423, isl1_425, isl1_426, \
                         ksk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pa_y[k] * isl0_423[k]
                   + f_16 * isk_337[k]
                   - f_14 * pc_y[k] * isl1_423[k];

        t_604[k] = f_15 * isk_338[k]
                   + f_3 * pc_y[k] * ksk_482[k];

        t_605[k] = pa_y[k] * isl0_425[k]
                   - f_14 * pc_y[k] * isl1_425[k];

        t_606[k] = pa_y[k] * isl0_426[k]
                   + f_20 * isk_339[k]
                   - f_14 * pc_y[k] * isl1_426[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pa_y, pc_y, pc_z, isl0_428, isl0_429, isk_303, \
                         isk_341, isk_342, isl1_428, isl1_429, \
                         ksk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_17 * isk_303[k]
                   + f_3 * pc_z[k] * ksk_483[k];

        t_608[k] = pa_y[k] * isl0_428[k]
                   + f_18 * isk_341[k]
                   - f_14 * pc_y[k] * isl1_428[k];

        t_609[k] = pa_y[k] * isl0_429[k]
                   + f_17 * isk_342[k]
                   - f_14 * pc_y[k] * isl1_429[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pa_y, pc_x, pc_y, isl0_430, isl0_432, \
                         isk_343, isk_344, isk_496, isl1_430, isl1_432, ksk_488, \
                         ksk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pa_y[k] * isl0_430[k]
                   + f_16 * isk_343[k]
                   - f_14 * pc_y[k] * isl1_430[k];

        t_611[k] = f_15 * isk_344[k]
                   + f_3 * pc_y[k] * ksk_488[k];

        t_612[k] = pa_y[k] * isl0_432[k]
                   - f_14 * pc_y[k] * isl1_432[k];

        t_613[k] = f_17 * isk_496[k]
                   + f_3 * pc_x[k] * ksk_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, isk_497, isk_498, isk_499, \
                         isk_500, isk_501, ksk_497, ksk_498, ksk_499, ksk_500, \
                         ksk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_17 * isk_497[k]
                   + f_3 * pc_x[k] * ksk_497[k];

        t_615[k] = f_17 * isk_498[k]
                   + f_3 * pc_x[k] * ksk_498[k];

        t_616[k] = f_17 * isk_499[k]
                   + f_3 * pc_x[k] * ksk_499[k];

        t_617[k] = f_17 * isk_500[k]
                   + f_3 * pc_x[k] * ksk_500[k];

        t_618[k] = f_17 * isk_501[k]
                   + f_3 * pc_x[k] * ksk_501[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, isk_316, isk_352, \
                         isk_502, isk_503, ksi0_385, ksi1_385, ksk_496, ksk_502, \
                         ksk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_17 * isk_502[k]
                   + f_3 * pc_x[k] * ksk_502[k];

        t_620[k] = f_17 * isk_503[k]
                   + f_3 * pc_x[k] * ksk_503[k];

        t_621[k] = f_15 * isk_352[k]
                   + f_1 * ksi0_385[k]
                   - f_2 * ksi1_385[k]
                   + f_3 * pc_y[k] * ksk_496[k];

        t_622[k] = f_17 * isk_316[k]
                   + f_3 * pc_z[k] * ksk_496[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_y, isk_354, isk_355, isk_356, ksi0_387, \
                         ksi0_388, ksi0_389, ksi1_387, ksi1_388, ksi1_389, ksk_498, ksk_499, \
                         ksk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_15 * isk_354[k]
                   + f_12 * ksi0_387[k]
                   - f_13 * ksi1_387[k]
                   + f_3 * pc_y[k] * ksk_498[k];

        t_624[k] = f_15 * isk_355[k]
                   + f_10 * ksi0_388[k]
                   - f_11 * ksi1_388[k]
                   + f_3 * pc_y[k] * ksk_499[k];

        t_625[k] = f_15 * isk_356[k]
                   + f_8 * ksi0_389[k]
                   - f_9 * ksi1_389[k]
                   + f_3 * pc_y[k] * ksk_500[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, isk_357, isk_358, isk_359, ksi0_390, \
                         ksi0_391, ksi1_390, ksi1_391, ksk_501, ksk_502, \
                         ksk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * isk_357[k]
                   + f_6 * ksi0_390[k]
                   - f_7 * ksi1_390[k]
                   + f_3 * pc_y[k] * ksk_501[k];

        t_627[k] = f_15 * isk_358[k]
                   + f_4 * ksi0_391[k]
                   - f_5 * ksi1_391[k]
                   + f_3 * pc_y[k] * ksk_502[k];

        t_628[k] = f_15 * isk_359[k]
                   + f_3 * pc_y[k] * ksk_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pa_y, pc_x, pc_y, pc_z, isl0_449, \
                         isk_324, isk_504, isl1_449, ksi0_392, ksi1_392, \
                         ksk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pa_y[k] * isl0_449[k]
                   - f_14 * pc_y[k] * isl1_449[k];

        t_630[k] = f_17 * isk_504[k]
                   + f_1 * ksi0_392[k]
                   - f_2 * ksi1_392[k]
                   + f_3 * pc_x[k] * ksk_504[k];

        t_631[k] = f_3 * pc_y[k] * ksk_504[k];

        t_632[k] = f_18 * isk_324[k]
                   + f_3 * pc_z[k] * ksk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, isk_509, ksi0_392, ksi0_397, \
                         ksi1_392, ksi1_397, ksk_505, ksk_506, \
                         ksk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_4 * ksi0_392[k]
                   - f_5 * ksi1_392[k]
                   + f_3 * pc_y[k] * ksk_505[k];

        t_634[k] = f_3 * pc_y[k] * ksk_506[k];

        t_635[k] = f_17 * isk_509[k]
                   + f_12 * ksi0_397[k]
                   - f_13 * ksi1_397[k]
                   + f_3 * pc_x[k] * ksk_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_y, ksi0_393, ksi0_394, ksi1_393, ksi1_394, \
                         ksk_507, ksk_508, ksk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_6 * ksi0_393[k]
                   - f_7 * ksi1_393[k]
                   + f_3 * pc_y[k] * ksk_507[k];

        t_637[k] = f_4 * ksi0_394[k]
                   - f_5 * ksi1_394[k]
                   + f_3 * pc_y[k] * ksk_508[k];

        t_638[k] = f_3 * pc_y[k] * ksk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_y, isk_513, ksi0_395, ksi0_396, \
                         ksi0_401, ksi1_395, ksi1_396, ksi1_401, ksk_510, ksk_511, \
                         ksk_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_17 * isk_513[k]
                   + f_10 * ksi0_401[k]
                   - f_11 * ksi1_401[k]
                   + f_3 * pc_x[k] * ksk_513[k];

        t_640[k] = f_8 * ksi0_395[k]
                   - f_9 * ksi1_395[k]
                   + f_3 * pc_y[k] * ksk_510[k];

        t_641[k] = f_6 * ksi0_396[k]
                   - f_7 * ksi1_396[k]
                   + f_3 * pc_y[k] * ksk_511[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, isk_518, ksi0_397, ksi0_406, \
                         ksi1_397, ksi1_406, ksk_512, ksk_513, \
                         ksk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_4 * ksi0_397[k]
                   - f_5 * ksi1_397[k]
                   + f_3 * pc_y[k] * ksk_512[k];

        t_643[k] = f_3 * pc_y[k] * ksk_513[k];

        t_644[k] = f_17 * isk_518[k]
                   + f_8 * ksi0_406[k]
                   - f_9 * ksi1_406[k]
                   + f_3 * pc_x[k] * ksk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_y, ksi0_398, ksi0_399, ksi0_400, ksi1_398, \
                         ksi1_399, ksi1_400, ksk_514, ksk_515, \
                         ksk_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_10 * ksi0_398[k]
                   - f_11 * ksi1_398[k]
                   + f_3 * pc_y[k] * ksk_514[k];

        t_646[k] = f_8 * ksi0_399[k]
                   - f_9 * ksi1_399[k]
                   + f_3 * pc_y[k] * ksk_515[k];

        t_647[k] = f_6 * ksi0_400[k]
                   - f_7 * ksi1_400[k]
                   + f_3 * pc_y[k] * ksk_516[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, pc_y, isk_524, ksi0_401, ksi0_412, \
                         ksi1_401, ksi1_412, ksk_517, ksk_518, \
                         ksk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_4 * ksi0_401[k]
                   - f_5 * ksi1_401[k]
                   + f_3 * pc_y[k] * ksk_517[k];

        t_649[k] = f_3 * pc_y[k] * ksk_518[k];

        t_650[k] = f_17 * isk_524[k]
                   + f_6 * ksi0_412[k]
                   - f_7 * ksi1_412[k]
                   + f_3 * pc_x[k] * ksk_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_y, ksi0_402, ksi0_403, ksi0_404, ksi1_402, \
                         ksi1_403, ksi1_404, ksk_519, ksk_520, \
                         ksk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_12 * ksi0_402[k]
                   - f_13 * ksi1_402[k]
                   + f_3 * pc_y[k] * ksk_519[k];

        t_652[k] = f_10 * ksi0_403[k]
                   - f_11 * ksi1_403[k]
                   + f_3 * pc_y[k] * ksk_520[k];

        t_653[k] = f_8 * ksi0_404[k]
                   - f_9 * ksi1_404[k]
                   + f_3 * pc_y[k] * ksk_521[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_y, ksi0_405, ksi0_406, ksi1_405, ksi1_406, \
                         ksk_522, ksk_523, ksk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_6 * ksi0_405[k]
                   - f_7 * ksi1_405[k]
                   + f_3 * pc_y[k] * ksk_522[k];

        t_655[k] = f_4 * ksi0_406[k]
                   - f_5 * ksi1_406[k]
                   + f_3 * pc_y[k] * ksk_523[k];

        t_656[k] = f_3 * pc_y[k] * ksk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pc_x, isk_531, isk_532, isk_533, isk_534, \
                         ksi0_419, ksi1_419, ksk_531, ksk_532, ksk_533, \
                         ksk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_17 * isk_531[k]
                   + f_4 * ksi0_419[k]
                   - f_5 * ksi1_419[k]
                   + f_3 * pc_x[k] * ksk_531[k];

        t_658[k] = f_17 * isk_532[k]
                   + f_3 * pc_x[k] * ksk_532[k];

        t_659[k] = f_17 * isk_533[k]
                   + f_3 * pc_x[k] * ksk_533[k];

        t_660[k] = f_17 * isk_534[k]
                   + f_3 * pc_x[k] * ksk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, pc_y, isk_535, isk_536, \
                         isk_537, isk_539, ksk_531, ksk_535, ksk_536, ksk_537, \
                         ksk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_17 * isk_535[k]
                   + f_3 * pc_x[k] * ksk_535[k];

        t_662[k] = f_17 * isk_536[k]
                   + f_3 * pc_x[k] * ksk_536[k];

        t_663[k] = f_17 * isk_537[k]
                   + f_3 * pc_x[k] * ksk_537[k];

        t_664[k] = f_3 * pc_y[k] * ksk_531[k];

        t_665[k] = f_17 * isk_539[k]
                   + f_3 * pc_x[k] * ksk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, pc_y, ksi0_413, ksi0_414, ksi0_415, ksi1_413, \
                         ksi1_414, ksi1_415, ksk_532, ksk_533, \
                         ksk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * ksi0_413[k]
                   - f_2 * ksi1_413[k]
                   + f_3 * pc_y[k] * ksk_532[k];

        t_667[k] = f_21 * ksi0_414[k]
                   - f_22 * ksi1_414[k]
                   + f_3 * pc_y[k] * ksk_533[k];

        t_668[k] = f_12 * ksi0_415[k]
                   - f_13 * ksi1_415[k]
                   + f_3 * pc_y[k] * ksk_534[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pc_y, ksi0_416, ksi0_417, ksi0_418, ksi1_416, \
                         ksi1_417, ksi1_418, ksk_535, ksk_536, \
                         ksk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_10 * ksi0_416[k]
                   - f_11 * ksi1_416[k]
                   + f_3 * pc_y[k] * ksk_535[k];

        t_670[k] = f_8 * ksi0_417[k]
                   - f_9 * ksi1_417[k]
                   + f_3 * pc_y[k] * ksk_536[k];

        t_671[k] = f_6 * ksi0_418[k]
                   - f_7 * ksi1_418[k]
                   + f_3 * pc_y[k] * ksk_537[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, pc_y, pc_z, isk_359, isk_540, \
                         ksi0_419, ksi0_420, ksi1_419, ksi1_420, ksk_538, ksk_539, \
                         ksk_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_4 * ksi0_419[k]
                   - f_5 * ksi1_419[k]
                   + f_3 * pc_y[k] * ksk_538[k];

        t_673[k] = f_3 * pc_y[k] * ksk_539[k];

        t_674[k] = f_18 * isk_359[k]
                   + f_1 * ksi0_419[k]
                   - f_2 * ksi1_419[k]
                   + f_3 * pc_z[k] * ksk_539[k];

        t_675[k] = f_16 * isk_540[k]
                   + f_1 * ksi0_420[k]
                   - f_2 * ksi1_420[k]
                   + f_3 * pc_x[k] * ksk_540[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pc_x, pc_y, pc_z, isk_360, isk_543, \
                         ksi0_423, ksi1_423, ksk_540, ksk_541, \
                         ksk_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_19 * isk_360[k]
                   + f_3 * pc_y[k] * ksk_540[k];

        t_677[k] = f_3 * pc_z[k] * ksk_540[k];

        t_678[k] = f_16 * isk_543[k]
                   + f_12 * ksi0_423[k]
                   - f_13 * ksi1_423[k]
                   + f_3 * pc_x[k] * ksk_543[k];

        t_679[k] = f_3 * pc_z[k] * ksk_541[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_z, isk_546, ksi0_420, ksi0_426, \
                         ksi1_420, ksi1_426, ksk_542, ksk_543, \
                         ksk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * ksi0_420[k]
                   - f_5 * ksi1_420[k]
                   + f_3 * pc_z[k] * ksk_542[k];

        t_681[k] = f_16 * isk_546[k]
                   + f_10 * ksi0_426[k]
                   - f_11 * ksi1_426[k]
                   + f_3 * pc_x[k] * ksk_546[k];

        t_682[k] = f_3 * pc_z[k] * ksk_543[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, pc_x, pc_y, pc_z, isk_365, isk_550, \
                         ksi0_422, ksi0_430, ksi1_422, ksi1_430, ksk_545, ksk_546, \
                         ksk_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_19 * isk_365[k]
                   + f_3 * pc_y[k] * ksk_545[k];

        t_684[k] = f_6 * ksi0_422[k]
                   - f_7 * ksi1_422[k]
                   + f_3 * pc_z[k] * ksk_545[k];

        t_685[k] = f_16 * isk_550[k]
                   + f_8 * ksi0_430[k]
                   - f_9 * ksi1_430[k]
                   + f_3 * pc_x[k] * ksk_550[k];

        t_686[k] = f_3 * pc_z[k] * ksk_546[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pc_y, pc_z, isk_369, ksi0_423, ksi0_425, \
                         ksi1_423, ksi1_425, ksk_547, ksk_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_4 * ksi0_423[k]
                   - f_5 * ksi1_423[k]
                   + f_3 * pc_z[k] * ksk_547[k];

        t_688[k] = f_19 * isk_369[k]
                   + f_3 * pc_y[k] * ksk_549[k];

        t_689[k] = f_8 * ksi0_425[k]
                   - f_9 * ksi1_425[k]
                   + f_3 * pc_z[k] * ksk_549[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, pc_x, pc_z, isk_555, ksi0_426, ksi0_435, \
                         ksi1_426, ksi1_435, ksk_550, ksk_551, \
                         ksk_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_16 * isk_555[k]
                   + f_6 * ksi0_435[k]
                   - f_7 * ksi1_435[k]
                   + f_3 * pc_x[k] * ksk_555[k];

        t_691[k] = f_3 * pc_z[k] * ksk_550[k];

        t_692[k] = f_4 * ksi0_426[k]
                   - f_5 * ksi1_426[k]
                   + f_3 * pc_z[k] * ksk_551[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pc_y, pc_z, isk_374, ksi0_427, ksi0_429, \
                         ksi1_427, ksi1_429, ksk_552, ksk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_6 * ksi0_427[k]
                   - f_7 * ksi1_427[k]
                   + f_3 * pc_z[k] * ksk_552[k];

        t_694[k] = f_19 * isk_374[k]
                   + f_3 * pc_y[k] * ksk_554[k];

        t_695[k] = f_10 * ksi0_429[k]
                   - f_11 * ksi1_429[k]
                   + f_3 * pc_z[k] * ksk_554[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pc_x, pc_z, isk_561, ksi0_430, ksi0_441, \
                         ksi1_430, ksi1_441, ksk_555, ksk_556, \
                         ksk_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_16 * isk_561[k]
                   + f_4 * ksi0_441[k]
                   - f_5 * ksi1_441[k]
                   + f_3 * pc_x[k] * ksk_561[k];

        t_697[k] = f_3 * pc_z[k] * ksk_555[k];

        t_698[k] = f_4 * ksi0_430[k]
                   - f_5 * ksi1_430[k]
                   + f_3 * pc_z[k] * ksk_556[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pc_y, pc_z, isk_380, ksi0_431, ksi0_432, \
                         ksi0_434, ksi1_431, ksi1_432, ksi1_434, ksk_557, ksk_558, \
                         ksk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_6 * ksi0_431[k]
                   - f_7 * ksi1_431[k]
                   + f_3 * pc_z[k] * ksk_557[k];

        t_700[k] = f_8 * ksi0_432[k]
                   - f_9 * ksi1_432[k]
                   + f_3 * pc_z[k] * ksk_558[k];

        t_701[k] = f_19 * isk_380[k]
                   + f_3 * pc_y[k] * ksk_560[k];

        t_702[k] = f_12 * ksi0_434[k]
                   - f_13 * ksi1_434[k]
                   + f_3 * pc_z[k] * ksk_560[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
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

    const auto *isl0_450 = buffer.data(isl0 + 450);
    const auto *isl0_453 = buffer.data(isl0 + 453);
    const auto *isl0_456 = buffer.data(isl0 + 456);
    const auto *isl0_460 = buffer.data(isl0 + 460);
    const auto *isl0_462 = buffer.data(isl0 + 462);
    const auto *isl0_465 = buffer.data(isl0 + 465);
    const auto *isl0_467 = buffer.data(isl0 + 467);
    const auto *isl0_468 = buffer.data(isl0 + 468);
    const auto *isl0_471 = buffer.data(isl0 + 471);
    const auto *isl0_473 = buffer.data(isl0 + 473);
    const auto *isl0_474 = buffer.data(isl0 + 474);
    const auto *isl0_475 = buffer.data(isl0 + 475);
    const auto *isl0_486 = buffer.data(isl0 + 486);

    const auto *isk_360 = buffer.data(isk + 360);
    const auto *isk_363 = buffer.data(isk + 363);
    const auto *isk_366 = buffer.data(isk + 366);
    const auto *isk_367 = buffer.data(isk + 367);
    const auto *isk_370 = buffer.data(isk + 370);
    const auto *isk_371 = buffer.data(isk + 371);
    const auto *isk_372 = buffer.data(isk + 372);
    const auto *isk_375 = buffer.data(isk + 375);
    const auto *isk_376 = buffer.data(isk + 376);
    const auto *isk_377 = buffer.data(isk + 377);
    const auto *isk_378 = buffer.data(isk + 378);
    const auto *isk_388 = buffer.data(isk + 388);
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
    const auto *isk_424 = buffer.data(isk + 424);
    const auto *isk_426 = buffer.data(isk + 426);
    const auto *isk_427 = buffer.data(isk + 427);
    const auto *isk_428 = buffer.data(isk + 428);
    const auto *isk_429 = buffer.data(isk + 429);
    const auto *isk_430 = buffer.data(isk + 430);
    const auto *isk_431 = buffer.data(isk + 431);
    const auto *isk_432 = buffer.data(isk + 432);
    const auto *isk_434 = buffer.data(isk + 434);
    const auto *isk_437 = buffer.data(isk + 437);
    const auto *isk_441 = buffer.data(isk + 441);
    const auto *isk_446 = buffer.data(isk + 446);
    const auto *isk_452 = buffer.data(isk + 452);
    const auto *isk_460 = buffer.data(isk + 460);
    const auto *isk_462 = buffer.data(isk + 462);
    const auto *isk_463 = buffer.data(isk + 463);
    const auto *isk_464 = buffer.data(isk + 464);
    const auto *isk_465 = buffer.data(isk + 465);
    const auto *isk_466 = buffer.data(isk + 466);
    const auto *isk_467 = buffer.data(isk + 467);
    const auto *isk_468 = buffer.data(isk + 468);
    const auto *isk_568 = buffer.data(isk + 568);
    const auto *isk_570 = buffer.data(isk + 570);
    const auto *isk_571 = buffer.data(isk + 571);
    const auto *isk_572 = buffer.data(isk + 572);
    const auto *isk_573 = buffer.data(isk + 573);
    const auto *isk_574 = buffer.data(isk + 574);
    const auto *isk_575 = buffer.data(isk + 575);
    const auto *isk_581 = buffer.data(isk + 581);
    const auto *isk_585 = buffer.data(isk + 585);
    const auto *isk_590 = buffer.data(isk + 590);
    const auto *isk_596 = buffer.data(isk + 596);
    const auto *isk_603 = buffer.data(isk + 603);
    const auto *isk_604 = buffer.data(isk + 604);
    const auto *isk_605 = buffer.data(isk + 605);
    const auto *isk_606 = buffer.data(isk + 606);
    const auto *isk_607 = buffer.data(isk + 607);
    const auto *isk_608 = buffer.data(isk + 608);
    const auto *isk_609 = buffer.data(isk + 609);
    const auto *isk_610 = buffer.data(isk + 610);
    const auto *isk_611 = buffer.data(isk + 611);
    const auto *isk_612 = buffer.data(isk + 612);
    const auto *isk_615 = buffer.data(isk + 615);
    const auto *isk_617 = buffer.data(isk + 617);
    const auto *isk_618 = buffer.data(isk + 618);
    const auto *isk_621 = buffer.data(isk + 621);
    const auto *isk_622 = buffer.data(isk + 622);
    const auto *isk_624 = buffer.data(isk + 624);
    const auto *isk_626 = buffer.data(isk + 626);
    const auto *isk_627 = buffer.data(isk + 627);
    const auto *isk_629 = buffer.data(isk + 629);
    const auto *isk_630 = buffer.data(isk + 630);
    const auto *isk_632 = buffer.data(isk + 632);
    const auto *isk_633 = buffer.data(isk + 633);
    const auto *isk_635 = buffer.data(isk + 635);
    const auto *isk_636 = buffer.data(isk + 636);
    const auto *isk_637 = buffer.data(isk + 637);
    const auto *isk_639 = buffer.data(isk + 639);
    const auto *isk_640 = buffer.data(isk + 640);
    const auto *isk_641 = buffer.data(isk + 641);
    const auto *isk_642 = buffer.data(isk + 642);
    const auto *isk_643 = buffer.data(isk + 643);
    const auto *isk_644 = buffer.data(isk + 644);
    const auto *isk_645 = buffer.data(isk + 645);
    const auto *isk_646 = buffer.data(isk + 646);
    const auto *isk_647 = buffer.data(isk + 647);
    const auto *isk_648 = buffer.data(isk + 648);

    const auto *isl1_450 = buffer.data(isl1 + 450);
    const auto *isl1_453 = buffer.data(isl1 + 453);
    const auto *isl1_456 = buffer.data(isl1 + 456);
    const auto *isl1_460 = buffer.data(isl1 + 460);
    const auto *isl1_462 = buffer.data(isl1 + 462);
    const auto *isl1_465 = buffer.data(isl1 + 465);
    const auto *isl1_467 = buffer.data(isl1 + 467);
    const auto *isl1_468 = buffer.data(isl1 + 468);
    const auto *isl1_471 = buffer.data(isl1 + 471);
    const auto *isl1_473 = buffer.data(isl1 + 473);
    const auto *isl1_474 = buffer.data(isl1 + 474);
    const auto *isl1_475 = buffer.data(isl1 + 475);
    const auto *isl1_486 = buffer.data(isl1 + 486);

    const auto *ksi0_441 = buffer.data(ksi0 + 441);
    const auto *ksi0_442 = buffer.data(ksi0 + 442);
    const auto *ksi0_443 = buffer.data(ksi0 + 443);
    const auto *ksi0_444 = buffer.data(ksi0 + 444);
    const auto *ksi0_445 = buffer.data(ksi0 + 445);
    const auto *ksi0_447 = buffer.data(ksi0 + 447);
    const auto *ksi0_453 = buffer.data(ksi0 + 453);
    const auto *ksi0_457 = buffer.data(ksi0 + 457);
    const auto *ksi0_462 = buffer.data(ksi0 + 462);
    const auto *ksi0_468 = buffer.data(ksi0 + 468);
    const auto *ksi0_471 = buffer.data(ksi0 + 471);
    const auto *ksi0_472 = buffer.data(ksi0 + 472);
    const auto *ksi0_473 = buffer.data(ksi0 + 473);
    const auto *ksi0_474 = buffer.data(ksi0 + 474);
    const auto *ksi0_475 = buffer.data(ksi0 + 475);
    const auto *ksi0_476 = buffer.data(ksi0 + 476);
    const auto *ksi0_479 = buffer.data(ksi0 + 479);
    const auto *ksi0_481 = buffer.data(ksi0 + 481);
    const auto *ksi0_482 = buffer.data(ksi0 + 482);
    const auto *ksi0_485 = buffer.data(ksi0 + 485);
    const auto *ksi0_486 = buffer.data(ksi0 + 486);
    const auto *ksi0_488 = buffer.data(ksi0 + 488);
    const auto *ksi0_490 = buffer.data(ksi0 + 490);
    const auto *ksi0_491 = buffer.data(ksi0 + 491);
    const auto *ksi0_493 = buffer.data(ksi0 + 493);
    const auto *ksi0_494 = buffer.data(ksi0 + 494);
    const auto *ksi0_496 = buffer.data(ksi0 + 496);
    const auto *ksi0_497 = buffer.data(ksi0 + 497);
    const auto *ksi0_499 = buffer.data(ksi0 + 499);
    const auto *ksi0_500 = buffer.data(ksi0 + 500);
    const auto *ksi0_501 = buffer.data(ksi0 + 501);
    const auto *ksi0_502 = buffer.data(ksi0 + 502);
    const auto *ksi0_503 = buffer.data(ksi0 + 503);
    const auto *ksi0_504 = buffer.data(ksi0 + 504);

    const auto *ksi1_441 = buffer.data(ksi1 + 441);
    const auto *ksi1_442 = buffer.data(ksi1 + 442);
    const auto *ksi1_443 = buffer.data(ksi1 + 443);
    const auto *ksi1_444 = buffer.data(ksi1 + 444);
    const auto *ksi1_445 = buffer.data(ksi1 + 445);
    const auto *ksi1_447 = buffer.data(ksi1 + 447);
    const auto *ksi1_453 = buffer.data(ksi1 + 453);
    const auto *ksi1_457 = buffer.data(ksi1 + 457);
    const auto *ksi1_462 = buffer.data(ksi1 + 462);
    const auto *ksi1_468 = buffer.data(ksi1 + 468);
    const auto *ksi1_471 = buffer.data(ksi1 + 471);
    const auto *ksi1_472 = buffer.data(ksi1 + 472);
    const auto *ksi1_473 = buffer.data(ksi1 + 473);
    const auto *ksi1_474 = buffer.data(ksi1 + 474);
    const auto *ksi1_475 = buffer.data(ksi1 + 475);
    const auto *ksi1_476 = buffer.data(ksi1 + 476);
    const auto *ksi1_479 = buffer.data(ksi1 + 479);
    const auto *ksi1_481 = buffer.data(ksi1 + 481);
    const auto *ksi1_482 = buffer.data(ksi1 + 482);
    const auto *ksi1_485 = buffer.data(ksi1 + 485);
    const auto *ksi1_486 = buffer.data(ksi1 + 486);
    const auto *ksi1_488 = buffer.data(ksi1 + 488);
    const auto *ksi1_490 = buffer.data(ksi1 + 490);
    const auto *ksi1_491 = buffer.data(ksi1 + 491);
    const auto *ksi1_493 = buffer.data(ksi1 + 493);
    const auto *ksi1_494 = buffer.data(ksi1 + 494);
    const auto *ksi1_496 = buffer.data(ksi1 + 496);
    const auto *ksi1_497 = buffer.data(ksi1 + 497);
    const auto *ksi1_499 = buffer.data(ksi1 + 499);
    const auto *ksi1_500 = buffer.data(ksi1 + 500);
    const auto *ksi1_501 = buffer.data(ksi1 + 501);
    const auto *ksi1_502 = buffer.data(ksi1 + 502);
    const auto *ksi1_503 = buffer.data(ksi1 + 503);
    const auto *ksi1_504 = buffer.data(ksi1 + 504);

    const auto *ksk_561 = buffer.data(ksk + 561);
    const auto *ksk_568 = buffer.data(ksk + 568);
    const auto *ksk_569 = buffer.data(ksk + 569);
    const auto *ksk_570 = buffer.data(ksk + 570);
    const auto *ksk_571 = buffer.data(ksk + 571);
    const auto *ksk_572 = buffer.data(ksk + 572);
    const auto *ksk_573 = buffer.data(ksk + 573);
    const auto *ksk_574 = buffer.data(ksk + 574);
    const auto *ksk_575 = buffer.data(ksk + 575);
    const auto *ksk_576 = buffer.data(ksk + 576);
    const auto *ksk_578 = buffer.data(ksk + 578);
    const auto *ksk_579 = buffer.data(ksk + 579);
    const auto *ksk_581 = buffer.data(ksk + 581);
    const auto *ksk_582 = buffer.data(ksk + 582);
    const auto *ksk_585 = buffer.data(ksk + 585);
    const auto *ksk_586 = buffer.data(ksk + 586);
    const auto *ksk_590 = buffer.data(ksk + 590);
    const auto *ksk_591 = buffer.data(ksk + 591);
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
    const auto *ksk_614 = buffer.data(ksk + 614);
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

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, pc_x, pc_z, isk_568, isk_570, \
                         isk_571, isk_572, ksk_561, ksk_568, ksk_570, ksk_571, \
                         ksk_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_16 * isk_568[k]
                   + f_3 * pc_x[k] * ksk_568[k];

        t_704[k] = f_3 * pc_z[k] * ksk_561[k];

        t_705[k] = f_16 * isk_570[k]
                   + f_3 * pc_x[k] * ksk_570[k];

        t_706[k] = f_16 * isk_571[k]
                   + f_3 * pc_x[k] * ksk_571[k];

        t_707[k] = f_16 * isk_572[k]
                   + f_3 * pc_x[k] * ksk_572[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, pc_x, pc_y, isk_388, isk_573, isk_574, \
                         isk_575, ksi0_441, ksi1_441, ksk_568, ksk_573, ksk_574, \
                         ksk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_16 * isk_573[k]
                   + f_3 * pc_x[k] * ksk_573[k];

        t_709[k] = f_16 * isk_574[k]
                   + f_3 * pc_x[k] * ksk_574[k];

        t_710[k] = f_16 * isk_575[k]
                   + f_3 * pc_x[k] * ksk_575[k];

        t_711[k] = f_19 * isk_388[k]
                   + f_1 * ksi0_441[k]
                   - f_2 * ksi1_441[k]
                   + f_3 * pc_y[k] * ksk_568[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pc_z, ksi0_441, ksi0_442, ksi0_443, \
                         ksi1_441, ksi1_442, ksi1_443, ksk_568, ksk_569, ksk_570, \
                         ksk_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_3 * pc_z[k] * ksk_568[k];

        t_713[k] = f_4 * ksi0_441[k]
                   - f_5 * ksi1_441[k]
                   + f_3 * pc_z[k] * ksk_569[k];

        t_714[k] = f_6 * ksi0_442[k]
                   - f_7 * ksi1_442[k]
                   + f_3 * pc_z[k] * ksk_570[k];

        t_715[k] = f_8 * ksi0_443[k]
                   - f_9 * ksi1_443[k]
                   + f_3 * pc_z[k] * ksk_571[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pc_y, pc_z, isk_395, ksi0_444, ksi0_445, \
                         ksi0_447, ksi1_444, ksi1_445, ksi1_447, ksk_572, ksk_573, \
                         ksk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_10 * ksi0_444[k]
                   - f_11 * ksi1_444[k]
                   + f_3 * pc_z[k] * ksk_572[k];

        t_717[k] = f_12 * ksi0_445[k]
                   - f_13 * ksi1_445[k]
                   + f_3 * pc_z[k] * ksk_573[k];

        t_718[k] = f_19 * isk_395[k]
                   + f_3 * pc_y[k] * ksk_575[k];

        t_719[k] = f_1 * ksi0_447[k]
                   - f_2 * ksi1_447[k]
                   + f_3 * pc_z[k] * ksk_575[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_z, pc_y, pc_z, isl0_450, isl0_453, \
                         isk_360, isk_396, isl1_450, isl1_453, \
                         ksk_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_z[k] * isl0_450[k]
                   - f_14 * pc_z[k] * isl1_450[k];

        t_721[k] = f_18 * isk_396[k]
                   + f_3 * pc_y[k] * ksk_576[k];

        t_722[k] = f_15 * isk_360[k]
                   + f_3 * pc_z[k] * ksk_576[k];

        t_723[k] = pa_z[k] * isl0_453[k]
                   - f_14 * pc_z[k] * isl1_453[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pa_z, pc_x, pc_y, pc_z, isl0_456, isk_398, \
                         isk_581, isl1_456, ksi0_453, ksi1_453, ksk_578, \
                         ksk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_18 * isk_398[k]
                   + f_3 * pc_y[k] * ksk_578[k];

        t_725[k] = f_16 * isk_581[k]
                   + f_12 * ksi0_453[k]
                   - f_13 * ksi1_453[k]
                   + f_3 * pc_x[k] * ksk_581[k];

        t_726[k] = pa_z[k] * isl0_456[k]
                   - f_14 * pc_z[k] * isl1_456[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_x, pc_y, pc_z, isk_363, isk_401, isk_585, \
                         ksi0_457, ksi1_457, ksk_579, ksk_581, \
                         ksk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_15 * isk_363[k]
                   + f_3 * pc_z[k] * ksk_579[k];

        t_728[k] = f_18 * isk_401[k]
                   + f_3 * pc_y[k] * ksk_581[k];

        t_729[k] = f_16 * isk_585[k]
                   + f_10 * ksi0_457[k]
                   - f_11 * ksi1_457[k]
                   + f_3 * pc_x[k] * ksk_585[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_z, pc_y, pc_z, isl0_460, isl0_462, \
                         isk_366, isk_367, isk_405, isl1_460, isl1_462, ksk_582, \
                         ksk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = pa_z[k] * isl0_460[k]
                   - f_14 * pc_z[k] * isl1_460[k];

        t_731[k] = f_15 * isk_366[k]
                   + f_3 * pc_z[k] * ksk_582[k];

        t_732[k] = pa_z[k] * isl0_462[k]
                   + f_16 * isk_367[k]
                   - f_14 * pc_z[k] * isl1_462[k];

        t_733[k] = f_18 * isk_405[k]
                   + f_3 * pc_y[k] * ksk_585[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pa_z, pc_x, pc_z, isl0_465, isk_370, isk_590, \
                         isl1_465, ksi0_462, ksi1_462, ksk_586, \
                         ksk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_16 * isk_590[k]
                   + f_8 * ksi0_462[k]
                   - f_9 * ksi1_462[k]
                   + f_3 * pc_x[k] * ksk_590[k];

        t_735[k] = pa_z[k] * isl0_465[k]
                   - f_14 * pc_z[k] * isl1_465[k];

        t_736[k] = f_15 * isk_370[k]
                   + f_3 * pc_z[k] * ksk_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pa_z, pc_y, pc_z, isl0_467, isl0_468, isk_371, \
                         isk_372, isk_410, isl1_467, isl1_468, \
                         ksk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pa_z[k] * isl0_467[k]
                   + f_16 * isk_371[k]
                   - f_14 * pc_z[k] * isl1_467[k];

        t_738[k] = pa_z[k] * isl0_468[k]
                   + f_17 * isk_372[k]
                   - f_14 * pc_z[k] * isl1_468[k];

        t_739[k] = f_18 * isk_410[k]
                   + f_3 * pc_y[k] * ksk_590[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pa_z, pc_x, pc_z, isl0_471, isk_375, isk_596, \
                         isl1_471, ksi0_468, ksi1_468, ksk_591, \
                         ksk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_16 * isk_596[k]
                   + f_6 * ksi0_468[k]
                   - f_7 * ksi1_468[k]
                   + f_3 * pc_x[k] * ksk_596[k];

        t_741[k] = pa_z[k] * isl0_471[k]
                   - f_14 * pc_z[k] * isl1_471[k];

        t_742[k] = f_15 * isk_375[k]
                   + f_3 * pc_z[k] * ksk_591[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, pa_z, pc_z, isl0_473, isl0_474, isl0_475, \
                         isk_376, isk_377, isk_378, isl1_473, isl1_474, \
                         isl1_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = pa_z[k] * isl0_473[k]
                   + f_16 * isk_376[k]
                   - f_14 * pc_z[k] * isl1_473[k];

        t_744[k] = pa_z[k] * isl0_474[k]
                   + f_17 * isk_377[k]
                   - f_14 * pc_z[k] * isl1_474[k];

        t_745[k] = pa_z[k] * isl0_475[k]
                   + f_18 * isk_378[k]
                   - f_14 * pc_z[k] * isl1_475[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, isk_416, isk_603, isk_604, \
                         isk_605, ksi0_475, ksi1_475, ksk_596, ksk_603, ksk_604, \
                         ksk_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_18 * isk_416[k]
                   + f_3 * pc_y[k] * ksk_596[k];

        t_747[k] = f_16 * isk_603[k]
                   + f_4 * ksi0_475[k]
                   - f_5 * ksi1_475[k]
                   + f_3 * pc_x[k] * ksk_603[k];

        t_748[k] = f_16 * isk_604[k]
                   + f_3 * pc_x[k] * ksk_604[k];

        t_749[k] = f_16 * isk_605[k]
                   + f_3 * pc_x[k] * ksk_605[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, pc_x, isk_606, isk_607, isk_608, \
                         isk_609, isk_610, ksk_606, ksk_607, ksk_608, ksk_609, \
                         ksk_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_16 * isk_606[k]
                   + f_3 * pc_x[k] * ksk_606[k];

        t_751[k] = f_16 * isk_607[k]
                   + f_3 * pc_x[k] * ksk_607[k];

        t_752[k] = f_16 * isk_608[k]
                   + f_3 * pc_x[k] * ksk_608[k];

        t_753[k] = f_16 * isk_609[k]
                   + f_3 * pc_x[k] * ksk_609[k];

        t_754[k] = f_16 * isk_610[k]
                   + f_3 * pc_x[k] * ksk_610[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, pa_z, pc_x, pc_z, isl0_486, isk_388, isk_611, \
                         isl1_486, ksk_604, ksk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_16 * isk_611[k]
                   + f_3 * pc_x[k] * ksk_611[k];

        t_756[k] = pa_z[k] * isl0_486[k]
                   - f_14 * pc_z[k] * isl1_486[k];

        t_757[k] = f_15 * isk_388[k]
                   + f_3 * pc_z[k] * ksk_604[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_y, isk_426, isk_427, isk_428, ksi0_471, \
                         ksi0_472, ksi0_473, ksi1_471, ksi1_472, ksi1_473, ksk_606, ksk_607, \
                         ksk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_18 * isk_426[k]
                   + f_12 * ksi0_471[k]
                   - f_13 * ksi1_471[k]
                   + f_3 * pc_y[k] * ksk_606[k];

        t_759[k] = f_18 * isk_427[k]
                   + f_10 * ksi0_472[k]
                   - f_11 * ksi1_472[k]
                   + f_3 * pc_y[k] * ksk_607[k];

        t_760[k] = f_18 * isk_428[k]
                   + f_8 * ksi0_473[k]
                   - f_9 * ksi1_473[k]
                   + f_3 * pc_y[k] * ksk_608[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_y, isk_429, isk_430, isk_431, ksi0_474, \
                         ksi0_475, ksi1_474, ksi1_475, ksk_609, ksk_610, \
                         ksk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_18 * isk_429[k]
                   + f_6 * ksi0_474[k]
                   - f_7 * ksi1_474[k]
                   + f_3 * pc_y[k] * ksk_609[k];

        t_762[k] = f_18 * isk_430[k]
                   + f_4 * ksi0_475[k]
                   - f_5 * ksi1_475[k]
                   + f_3 * pc_y[k] * ksk_610[k];

        t_763[k] = f_18 * isk_431[k]
                   + f_3 * pc_y[k] * ksk_611[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pc_x, pc_y, pc_z, isk_395, isk_432, isk_612, \
                         ksi0_475, ksi0_476, ksi1_475, ksi1_476, ksk_611, \
                         ksk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_15 * isk_395[k]
                   + f_1 * ksi0_475[k]
                   - f_2 * ksi1_475[k]
                   + f_3 * pc_z[k] * ksk_611[k];

        t_765[k] = f_16 * isk_612[k]
                   + f_1 * ksi0_476[k]
                   - f_2 * ksi1_476[k]
                   + f_3 * pc_x[k] * ksk_612[k];

        t_766[k] = f_17 * isk_432[k]
                   + f_3 * pc_y[k] * ksk_612[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_x, pc_y, pc_z, isk_396, isk_434, isk_615, \
                         ksi0_479, ksi1_479, ksk_612, ksk_614, \
                         ksk_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_16 * isk_396[k]
                   + f_3 * pc_z[k] * ksk_612[k];

        t_768[k] = f_16 * isk_615[k]
                   + f_12 * ksi0_479[k]
                   - f_13 * ksi1_479[k]
                   + f_3 * pc_x[k] * ksk_615[k];

        t_769[k] = f_17 * isk_434[k]
                   + f_3 * pc_y[k] * ksk_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pc_x, pc_z, isk_399, isk_617, isk_618, ksi0_481, \
                         ksi0_482, ksi1_481, ksi1_482, ksk_615, ksk_617, \
                         ksk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_16 * isk_617[k]
                   + f_12 * ksi0_481[k]
                   - f_13 * ksi1_481[k]
                   + f_3 * pc_x[k] * ksk_617[k];

        t_771[k] = f_16 * isk_618[k]
                   + f_10 * ksi0_482[k]
                   - f_11 * ksi1_482[k]
                   + f_3 * pc_x[k] * ksk_618[k];

        t_772[k] = f_16 * isk_399[k]
                   + f_3 * pc_z[k] * ksk_615[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pc_x, pc_y, isk_437, isk_621, isk_622, ksi0_485, \
                         ksi0_486, ksi1_485, ksi1_486, ksk_617, ksk_621, \
                         ksk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_17 * isk_437[k]
                   + f_3 * pc_y[k] * ksk_617[k];

        t_774[k] = f_16 * isk_621[k]
                   + f_10 * ksi0_485[k]
                   - f_11 * ksi1_485[k]
                   + f_3 * pc_x[k] * ksk_621[k];

        t_775[k] = f_16 * isk_622[k]
                   + f_8 * ksi0_486[k]
                   - f_9 * ksi1_486[k]
                   + f_3 * pc_x[k] * ksk_622[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pc_x, pc_y, pc_z, isk_402, isk_441, isk_624, \
                         ksi0_488, ksi1_488, ksk_618, ksk_621, \
                         ksk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_16 * isk_402[k]
                   + f_3 * pc_z[k] * ksk_618[k];

        t_777[k] = f_16 * isk_624[k]
                   + f_8 * ksi0_488[k]
                   - f_9 * ksi1_488[k]
                   + f_3 * pc_x[k] * ksk_624[k];

        t_778[k] = f_17 * isk_441[k]
                   + f_3 * pc_y[k] * ksk_621[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pc_x, pc_z, isk_406, isk_626, isk_627, ksi0_490, \
                         ksi0_491, ksi1_490, ksi1_491, ksk_622, ksk_626, \
                         ksk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_16 * isk_626[k]
                   + f_8 * ksi0_490[k]
                   - f_9 * ksi1_490[k]
                   + f_3 * pc_x[k] * ksk_626[k];

        t_780[k] = f_16 * isk_627[k]
                   + f_6 * ksi0_491[k]
                   - f_7 * ksi1_491[k]
                   + f_3 * pc_x[k] * ksk_627[k];

        t_781[k] = f_16 * isk_406[k]
                   + f_3 * pc_z[k] * ksk_622[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pc_x, pc_y, isk_446, isk_629, isk_630, ksi0_493, \
                         ksi0_494, ksi1_493, ksi1_494, ksk_626, ksk_629, \
                         ksk_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_16 * isk_629[k]
                   + f_6 * ksi0_493[k]
                   - f_7 * ksi1_493[k]
                   + f_3 * pc_x[k] * ksk_629[k];

        t_783[k] = f_16 * isk_630[k]
                   + f_6 * ksi0_494[k]
                   - f_7 * ksi1_494[k]
                   + f_3 * pc_x[k] * ksk_630[k];

        t_784[k] = f_17 * isk_446[k]
                   + f_3 * pc_y[k] * ksk_626[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pc_x, pc_z, isk_411, isk_632, isk_633, ksi0_496, \
                         ksi0_497, ksi1_496, ksi1_497, ksk_627, ksk_632, \
                         ksk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_16 * isk_632[k]
                   + f_6 * ksi0_496[k]
                   - f_7 * ksi1_496[k]
                   + f_3 * pc_x[k] * ksk_632[k];

        t_786[k] = f_16 * isk_633[k]
                   + f_4 * ksi0_497[k]
                   - f_5 * ksi1_497[k]
                   + f_3 * pc_x[k] * ksk_633[k];

        t_787[k] = f_16 * isk_411[k]
                   + f_3 * pc_z[k] * ksk_627[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, isk_635, isk_636, isk_637, ksi0_499, \
                         ksi0_500, ksi0_501, ksi1_499, ksi1_500, ksi1_501, ksk_635, ksk_636, \
                         ksk_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_16 * isk_635[k]
                   + f_4 * ksi0_499[k]
                   - f_5 * ksi1_499[k]
                   + f_3 * pc_x[k] * ksk_635[k];

        t_789[k] = f_16 * isk_636[k]
                   + f_4 * ksi0_500[k]
                   - f_5 * ksi1_500[k]
                   + f_3 * pc_x[k] * ksk_636[k];

        t_790[k] = f_16 * isk_637[k]
                   + f_4 * ksi0_501[k]
                   - f_5 * ksi1_501[k]
                   + f_3 * pc_x[k] * ksk_637[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pc_x, pc_y, isk_452, isk_639, isk_640, \
                         isk_641, ksi0_503, ksi1_503, ksk_632, ksk_639, ksk_640, \
                         ksk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_17 * isk_452[k]
                   + f_3 * pc_y[k] * ksk_632[k];

        t_792[k] = f_16 * isk_639[k]
                   + f_4 * ksi0_503[k]
                   - f_5 * ksi1_503[k]
                   + f_3 * pc_x[k] * ksk_639[k];

        t_793[k] = f_16 * isk_640[k]
                   + f_3 * pc_x[k] * ksk_640[k];

        t_794[k] = f_16 * isk_641[k]
                   + f_3 * pc_x[k] * ksk_641[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, pc_x, isk_642, isk_643, isk_644, \
                         isk_645, isk_646, ksk_642, ksk_643, ksk_644, ksk_645, \
                         ksk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_16 * isk_642[k]
                   + f_3 * pc_x[k] * ksk_642[k];

        t_796[k] = f_16 * isk_643[k]
                   + f_3 * pc_x[k] * ksk_643[k];

        t_797[k] = f_16 * isk_644[k]
                   + f_3 * pc_x[k] * ksk_644[k];

        t_798[k] = f_16 * isk_645[k]
                   + f_3 * pc_x[k] * ksk_645[k];

        t_799[k] = f_16 * isk_646[k]
                   + f_3 * pc_x[k] * ksk_646[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, isk_424, isk_460, isk_647, \
                         ksi0_497, ksi1_497, ksk_640, ksk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_16 * isk_647[k]
                   + f_3 * pc_x[k] * ksk_647[k];

        t_801[k] = f_17 * isk_460[k]
                   + f_1 * ksi0_497[k]
                   - f_2 * ksi1_497[k]
                   + f_3 * pc_y[k] * ksk_640[k];

        t_802[k] = f_16 * isk_424[k]
                   + f_3 * pc_z[k] * ksk_640[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_y, isk_462, isk_463, isk_464, ksi0_499, \
                         ksi0_500, ksi0_501, ksi1_499, ksi1_500, ksi1_501, ksk_642, ksk_643, \
                         ksk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_17 * isk_462[k]
                   + f_12 * ksi0_499[k]
                   - f_13 * ksi1_499[k]
                   + f_3 * pc_y[k] * ksk_642[k];

        t_804[k] = f_17 * isk_463[k]
                   + f_10 * ksi0_500[k]
                   - f_11 * ksi1_500[k]
                   + f_3 * pc_y[k] * ksk_643[k];

        t_805[k] = f_17 * isk_464[k]
                   + f_8 * ksi0_501[k]
                   - f_9 * ksi1_501[k]
                   + f_3 * pc_y[k] * ksk_644[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pc_y, isk_465, isk_466, isk_467, ksi0_502, \
                         ksi0_503, ksi1_502, ksi1_503, ksk_645, ksk_646, \
                         ksk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * isk_465[k]
                   + f_6 * ksi0_502[k]
                   - f_7 * ksi1_502[k]
                   + f_3 * pc_y[k] * ksk_645[k];

        t_807[k] = f_17 * isk_466[k]
                   + f_4 * ksi0_503[k]
                   - f_5 * ksi1_503[k]
                   + f_3 * pc_y[k] * ksk_646[k];

        t_808[k] = f_17 * isk_467[k]
                   + f_3 * pc_y[k] * ksk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pc_x, pc_y, pc_z, isk_431, isk_468, isk_648, \
                         ksi0_503, ksi0_504, ksi1_503, ksi1_504, ksk_647, \
                         ksk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_16 * isk_431[k]
                   + f_1 * ksi0_503[k]
                   - f_2 * ksi1_503[k]
                   + f_3 * pc_z[k] * ksk_647[k];

        t_810[k] = f_16 * isk_648[k]
                   + f_1 * ksi0_504[k]
                   - f_2 * ksi1_504[k]
                   + f_3 * pc_x[k] * ksk_648[k];

        t_811[k] = f_16 * isk_468[k]
                   + f_3 * pc_y[k] * ksk_648[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
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

    const auto *isl0_630 = buffer.data(isl0 + 630);
    const auto *isl0_633 = buffer.data(isl0 + 633);
    const auto *isl0_635 = buffer.data(isl0 + 635);
    const auto *isl0_636 = buffer.data(isl0 + 636);
    const auto *isl0_639 = buffer.data(isl0 + 639);
    const auto *isl0_640 = buffer.data(isl0 + 640);
    const auto *isl0_642 = buffer.data(isl0 + 642);
    const auto *isl0_644 = buffer.data(isl0 + 644);
    const auto *isl0_645 = buffer.data(isl0 + 645);
    const auto *isl0_647 = buffer.data(isl0 + 647);
    const auto *isl0_648 = buffer.data(isl0 + 648);
    const auto *isl0_650 = buffer.data(isl0 + 650);
    const auto *isl0_651 = buffer.data(isl0 + 651);
    const auto *isl0_653 = buffer.data(isl0 + 653);
    const auto *isl0_654 = buffer.data(isl0 + 654);
    const auto *isl0_655 = buffer.data(isl0 + 655);
    const auto *isl0_657 = buffer.data(isl0 + 657);
    const auto *isl0_674 = buffer.data(isl0 + 674);

    const auto *isk_432 = buffer.data(isk + 432);
    const auto *isk_435 = buffer.data(isk + 435);
    const auto *isk_438 = buffer.data(isk + 438);
    const auto *isk_442 = buffer.data(isk + 442);
    const auto *isk_447 = buffer.data(isk + 447);
    const auto *isk_460 = buffer.data(isk + 460);
    const auto *isk_467 = buffer.data(isk + 467);
    const auto *isk_468 = buffer.data(isk + 468);
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
    const auto *isk_509 = buffer.data(isk + 509);
    const auto *isk_510 = buffer.data(isk + 510);
    const auto *isk_512 = buffer.data(isk + 512);
    const auto *isk_513 = buffer.data(isk + 513);
    const auto *isk_514 = buffer.data(isk + 514);
    const auto *isk_516 = buffer.data(isk + 516);
    const auto *isk_517 = buffer.data(isk + 517);
    const auto *isk_518 = buffer.data(isk + 518);
    const auto *isk_519 = buffer.data(isk + 519);
    const auto *isk_521 = buffer.data(isk + 521);
    const auto *isk_522 = buffer.data(isk + 522);
    const auto *isk_523 = buffer.data(isk + 523);
    const auto *isk_524 = buffer.data(isk + 524);
    const auto *isk_532 = buffer.data(isk + 532);
    const auto *isk_534 = buffer.data(isk + 534);
    const auto *isk_535 = buffer.data(isk + 535);
    const auto *isk_536 = buffer.data(isk + 536);
    const auto *isk_537 = buffer.data(isk + 537);
    const auto *isk_538 = buffer.data(isk + 538);
    const auto *isk_539 = buffer.data(isk + 539);
    const auto *isk_651 = buffer.data(isk + 651);
    const auto *isk_653 = buffer.data(isk + 653);
    const auto *isk_654 = buffer.data(isk + 654);
    const auto *isk_657 = buffer.data(isk + 657);
    const auto *isk_658 = buffer.data(isk + 658);
    const auto *isk_660 = buffer.data(isk + 660);
    const auto *isk_662 = buffer.data(isk + 662);
    const auto *isk_663 = buffer.data(isk + 663);
    const auto *isk_665 = buffer.data(isk + 665);
    const auto *isk_666 = buffer.data(isk + 666);
    const auto *isk_668 = buffer.data(isk + 668);
    const auto *isk_669 = buffer.data(isk + 669);
    const auto *isk_671 = buffer.data(isk + 671);
    const auto *isk_672 = buffer.data(isk + 672);
    const auto *isk_673 = buffer.data(isk + 673);
    const auto *isk_675 = buffer.data(isk + 675);
    const auto *isk_676 = buffer.data(isk + 676);
    const auto *isk_677 = buffer.data(isk + 677);
    const auto *isk_678 = buffer.data(isk + 678);
    const auto *isk_679 = buffer.data(isk + 679);
    const auto *isk_680 = buffer.data(isk + 680);
    const auto *isk_681 = buffer.data(isk + 681);
    const auto *isk_682 = buffer.data(isk + 682);
    const auto *isk_683 = buffer.data(isk + 683);
    const auto *isk_712 = buffer.data(isk + 712);
    const auto *isk_713 = buffer.data(isk + 713);
    const auto *isk_714 = buffer.data(isk + 714);
    const auto *isk_715 = buffer.data(isk + 715);
    const auto *isk_716 = buffer.data(isk + 716);
    const auto *isk_717 = buffer.data(isk + 717);
    const auto *isk_718 = buffer.data(isk + 718);
    const auto *isk_719 = buffer.data(isk + 719);
    const auto *isk_720 = buffer.data(isk + 720);
    const auto *isk_725 = buffer.data(isk + 725);
    const auto *isk_729 = buffer.data(isk + 729);
    const auto *isk_734 = buffer.data(isk + 734);
    const auto *isk_740 = buffer.data(isk + 740);

    const auto *isl1_630 = buffer.data(isl1 + 630);
    const auto *isl1_633 = buffer.data(isl1 + 633);
    const auto *isl1_635 = buffer.data(isl1 + 635);
    const auto *isl1_636 = buffer.data(isl1 + 636);
    const auto *isl1_639 = buffer.data(isl1 + 639);
    const auto *isl1_640 = buffer.data(isl1 + 640);
    const auto *isl1_642 = buffer.data(isl1 + 642);
    const auto *isl1_644 = buffer.data(isl1 + 644);
    const auto *isl1_645 = buffer.data(isl1 + 645);
    const auto *isl1_647 = buffer.data(isl1 + 647);
    const auto *isl1_648 = buffer.data(isl1 + 648);
    const auto *isl1_650 = buffer.data(isl1 + 650);
    const auto *isl1_651 = buffer.data(isl1 + 651);
    const auto *isl1_653 = buffer.data(isl1 + 653);
    const auto *isl1_654 = buffer.data(isl1 + 654);
    const auto *isl1_655 = buffer.data(isl1 + 655);
    const auto *isl1_657 = buffer.data(isl1 + 657);
    const auto *isl1_674 = buffer.data(isl1 + 674);

    const auto *ksi0_507 = buffer.data(ksi0 + 507);
    const auto *ksi0_509 = buffer.data(ksi0 + 509);
    const auto *ksi0_510 = buffer.data(ksi0 + 510);
    const auto *ksi0_513 = buffer.data(ksi0 + 513);
    const auto *ksi0_514 = buffer.data(ksi0 + 514);
    const auto *ksi0_516 = buffer.data(ksi0 + 516);
    const auto *ksi0_518 = buffer.data(ksi0 + 518);
    const auto *ksi0_519 = buffer.data(ksi0 + 519);
    const auto *ksi0_521 = buffer.data(ksi0 + 521);
    const auto *ksi0_522 = buffer.data(ksi0 + 522);
    const auto *ksi0_524 = buffer.data(ksi0 + 524);
    const auto *ksi0_525 = buffer.data(ksi0 + 525);
    const auto *ksi0_527 = buffer.data(ksi0 + 527);
    const auto *ksi0_528 = buffer.data(ksi0 + 528);
    const auto *ksi0_529 = buffer.data(ksi0 + 529);
    const auto *ksi0_530 = buffer.data(ksi0 + 530);
    const auto *ksi0_531 = buffer.data(ksi0 + 531);
    const auto *ksi0_553 = buffer.data(ksi0 + 553);
    const auto *ksi0_555 = buffer.data(ksi0 + 555);
    const auto *ksi0_556 = buffer.data(ksi0 + 556);
    const auto *ksi0_557 = buffer.data(ksi0 + 557);
    const auto *ksi0_558 = buffer.data(ksi0 + 558);
    const auto *ksi0_559 = buffer.data(ksi0 + 559);
    const auto *ksi0_560 = buffer.data(ksi0 + 560);
    const auto *ksi0_561 = buffer.data(ksi0 + 561);
    const auto *ksi0_562 = buffer.data(ksi0 + 562);
    const auto *ksi0_563 = buffer.data(ksi0 + 563);
    const auto *ksi0_564 = buffer.data(ksi0 + 564);
    const auto *ksi0_565 = buffer.data(ksi0 + 565);
    const auto *ksi0_566 = buffer.data(ksi0 + 566);
    const auto *ksi0_567 = buffer.data(ksi0 + 567);
    const auto *ksi0_568 = buffer.data(ksi0 + 568);
    const auto *ksi0_569 = buffer.data(ksi0 + 569);
    const auto *ksi0_574 = buffer.data(ksi0 + 574);
    const auto *ksi0_580 = buffer.data(ksi0 + 580);

    const auto *ksi1_507 = buffer.data(ksi1 + 507);
    const auto *ksi1_509 = buffer.data(ksi1 + 509);
    const auto *ksi1_510 = buffer.data(ksi1 + 510);
    const auto *ksi1_513 = buffer.data(ksi1 + 513);
    const auto *ksi1_514 = buffer.data(ksi1 + 514);
    const auto *ksi1_516 = buffer.data(ksi1 + 516);
    const auto *ksi1_518 = buffer.data(ksi1 + 518);
    const auto *ksi1_519 = buffer.data(ksi1 + 519);
    const auto *ksi1_521 = buffer.data(ksi1 + 521);
    const auto *ksi1_522 = buffer.data(ksi1 + 522);
    const auto *ksi1_524 = buffer.data(ksi1 + 524);
    const auto *ksi1_525 = buffer.data(ksi1 + 525);
    const auto *ksi1_527 = buffer.data(ksi1 + 527);
    const auto *ksi1_528 = buffer.data(ksi1 + 528);
    const auto *ksi1_529 = buffer.data(ksi1 + 529);
    const auto *ksi1_530 = buffer.data(ksi1 + 530);
    const auto *ksi1_531 = buffer.data(ksi1 + 531);
    const auto *ksi1_553 = buffer.data(ksi1 + 553);
    const auto *ksi1_555 = buffer.data(ksi1 + 555);
    const auto *ksi1_556 = buffer.data(ksi1 + 556);
    const auto *ksi1_557 = buffer.data(ksi1 + 557);
    const auto *ksi1_558 = buffer.data(ksi1 + 558);
    const auto *ksi1_559 = buffer.data(ksi1 + 559);
    const auto *ksi1_560 = buffer.data(ksi1 + 560);
    const auto *ksi1_561 = buffer.data(ksi1 + 561);
    const auto *ksi1_562 = buffer.data(ksi1 + 562);
    const auto *ksi1_563 = buffer.data(ksi1 + 563);
    const auto *ksi1_564 = buffer.data(ksi1 + 564);
    const auto *ksi1_565 = buffer.data(ksi1 + 565);
    const auto *ksi1_566 = buffer.data(ksi1 + 566);
    const auto *ksi1_567 = buffer.data(ksi1 + 567);
    const auto *ksi1_568 = buffer.data(ksi1 + 568);
    const auto *ksi1_569 = buffer.data(ksi1 + 569);
    const auto *ksi1_574 = buffer.data(ksi1 + 574);
    const auto *ksi1_580 = buffer.data(ksi1 + 580);

    const auto *ksk_648 = buffer.data(ksk + 648);
    const auto *ksk_650 = buffer.data(ksk + 650);
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
    const auto *ksk_684 = buffer.data(ksk + 684);
    const auto *ksk_686 = buffer.data(ksk + 686);
    const auto *ksk_687 = buffer.data(ksk + 687);
    const auto *ksk_689 = buffer.data(ksk + 689);
    const auto *ksk_690 = buffer.data(ksk + 690);
    const auto *ksk_693 = buffer.data(ksk + 693);
    const auto *ksk_694 = buffer.data(ksk + 694);
    const auto *ksk_698 = buffer.data(ksk + 698);
    const auto *ksk_699 = buffer.data(ksk + 699);
    const auto *ksk_704 = buffer.data(ksk + 704);
    const auto *ksk_712 = buffer.data(ksk + 712);
    const auto *ksk_713 = buffer.data(ksk + 713);
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
    const auto *ksk_724 = buffer.data(ksk + 724);
    const auto *ksk_725 = buffer.data(ksk + 725);
    const auto *ksk_726 = buffer.data(ksk + 726);
    const auto *ksk_727 = buffer.data(ksk + 727);
    const auto *ksk_728 = buffer.data(ksk + 728);
    const auto *ksk_729 = buffer.data(ksk + 729);
    const auto *ksk_730 = buffer.data(ksk + 730);
    const auto *ksk_731 = buffer.data(ksk + 731);
    const auto *ksk_732 = buffer.data(ksk + 732);
    const auto *ksk_733 = buffer.data(ksk + 733);
    const auto *ksk_734 = buffer.data(ksk + 734);
    const auto *ksk_740 = buffer.data(ksk + 740);

#pragma omp simd aligned(t_812, t_813, t_814, pc_x, pc_y, pc_z, isk_432, isk_470, isk_651, \
                         ksi0_507, ksi1_507, ksk_648, ksk_650, \
                         ksk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_17 * isk_432[k]
                   + f_3 * pc_z[k] * ksk_648[k];

        t_813[k] = f_16 * isk_651[k]
                   + f_12 * ksi0_507[k]
                   - f_13 * ksi1_507[k]
                   + f_3 * pc_x[k] * ksk_651[k];

        t_814[k] = f_16 * isk_470[k]
                   + f_3 * pc_y[k] * ksk_650[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pc_x, pc_z, isk_435, isk_653, isk_654, ksi0_509, \
                         ksi0_510, ksi1_509, ksi1_510, ksk_651, ksk_653, \
                         ksk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_16 * isk_653[k]
                   + f_12 * ksi0_509[k]
                   - f_13 * ksi1_509[k]
                   + f_3 * pc_x[k] * ksk_653[k];

        t_816[k] = f_16 * isk_654[k]
                   + f_10 * ksi0_510[k]
                   - f_11 * ksi1_510[k]
                   + f_3 * pc_x[k] * ksk_654[k];

        t_817[k] = f_17 * isk_435[k]
                   + f_3 * pc_z[k] * ksk_651[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pc_x, pc_y, isk_473, isk_657, isk_658, ksi0_513, \
                         ksi0_514, ksi1_513, ksi1_514, ksk_653, ksk_657, \
                         ksk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_16 * isk_473[k]
                   + f_3 * pc_y[k] * ksk_653[k];

        t_819[k] = f_16 * isk_657[k]
                   + f_10 * ksi0_513[k]
                   - f_11 * ksi1_513[k]
                   + f_3 * pc_x[k] * ksk_657[k];

        t_820[k] = f_16 * isk_658[k]
                   + f_8 * ksi0_514[k]
                   - f_9 * ksi1_514[k]
                   + f_3 * pc_x[k] * ksk_658[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pc_x, pc_y, pc_z, isk_438, isk_477, isk_660, \
                         ksi0_516, ksi1_516, ksk_654, ksk_657, \
                         ksk_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_17 * isk_438[k]
                   + f_3 * pc_z[k] * ksk_654[k];

        t_822[k] = f_16 * isk_660[k]
                   + f_8 * ksi0_516[k]
                   - f_9 * ksi1_516[k]
                   + f_3 * pc_x[k] * ksk_660[k];

        t_823[k] = f_16 * isk_477[k]
                   + f_3 * pc_y[k] * ksk_657[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, isk_442, isk_662, isk_663, ksi0_518, \
                         ksi0_519, ksi1_518, ksi1_519, ksk_658, ksk_662, \
                         ksk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_16 * isk_662[k]
                   + f_8 * ksi0_518[k]
                   - f_9 * ksi1_518[k]
                   + f_3 * pc_x[k] * ksk_662[k];

        t_825[k] = f_16 * isk_663[k]
                   + f_6 * ksi0_519[k]
                   - f_7 * ksi1_519[k]
                   + f_3 * pc_x[k] * ksk_663[k];

        t_826[k] = f_17 * isk_442[k]
                   + f_3 * pc_z[k] * ksk_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, pc_x, pc_y, isk_482, isk_665, isk_666, ksi0_521, \
                         ksi0_522, ksi1_521, ksi1_522, ksk_662, ksk_665, \
                         ksk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_16 * isk_665[k]
                   + f_6 * ksi0_521[k]
                   - f_7 * ksi1_521[k]
                   + f_3 * pc_x[k] * ksk_665[k];

        t_828[k] = f_16 * isk_666[k]
                   + f_6 * ksi0_522[k]
                   - f_7 * ksi1_522[k]
                   + f_3 * pc_x[k] * ksk_666[k];

        t_829[k] = f_16 * isk_482[k]
                   + f_3 * pc_y[k] * ksk_662[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, pc_x, pc_z, isk_447, isk_668, isk_669, ksi0_524, \
                         ksi0_525, ksi1_524, ksi1_525, ksk_663, ksk_668, \
                         ksk_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_16 * isk_668[k]
                   + f_6 * ksi0_524[k]
                   - f_7 * ksi1_524[k]
                   + f_3 * pc_x[k] * ksk_668[k];

        t_831[k] = f_16 * isk_669[k]
                   + f_4 * ksi0_525[k]
                   - f_5 * ksi1_525[k]
                   + f_3 * pc_x[k] * ksk_669[k];

        t_832[k] = f_17 * isk_447[k]
                   + f_3 * pc_z[k] * ksk_663[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, isk_671, isk_672, isk_673, ksi0_527, \
                         ksi0_528, ksi0_529, ksi1_527, ksi1_528, ksi1_529, ksk_671, ksk_672, \
                         ksk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_16 * isk_671[k]
                   + f_4 * ksi0_527[k]
                   - f_5 * ksi1_527[k]
                   + f_3 * pc_x[k] * ksk_671[k];

        t_834[k] = f_16 * isk_672[k]
                   + f_4 * ksi0_528[k]
                   - f_5 * ksi1_528[k]
                   + f_3 * pc_x[k] * ksk_672[k];

        t_835[k] = f_16 * isk_673[k]
                   + f_4 * ksi0_529[k]
                   - f_5 * ksi1_529[k]
                   + f_3 * pc_x[k] * ksk_673[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pc_x, pc_y, isk_488, isk_675, isk_676, \
                         isk_677, ksi0_531, ksi1_531, ksk_668, ksk_675, ksk_676, \
                         ksk_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * isk_488[k]
                   + f_3 * pc_y[k] * ksk_668[k];

        t_837[k] = f_16 * isk_675[k]
                   + f_4 * ksi0_531[k]
                   - f_5 * ksi1_531[k]
                   + f_3 * pc_x[k] * ksk_675[k];

        t_838[k] = f_16 * isk_676[k]
                   + f_3 * pc_x[k] * ksk_676[k];

        t_839[k] = f_16 * isk_677[k]
                   + f_3 * pc_x[k] * ksk_677[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, pc_x, isk_678, isk_679, isk_680, \
                         isk_681, isk_682, ksk_678, ksk_679, ksk_680, ksk_681, \
                         ksk_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_16 * isk_678[k]
                   + f_3 * pc_x[k] * ksk_678[k];

        t_841[k] = f_16 * isk_679[k]
                   + f_3 * pc_x[k] * ksk_679[k];

        t_842[k] = f_16 * isk_680[k]
                   + f_3 * pc_x[k] * ksk_680[k];

        t_843[k] = f_16 * isk_681[k]
                   + f_3 * pc_x[k] * ksk_681[k];

        t_844[k] = f_16 * isk_682[k]
                   + f_3 * pc_x[k] * ksk_682[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, pc_z, isk_460, isk_496, isk_683, \
                         ksi0_525, ksi1_525, ksk_676, ksk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_16 * isk_683[k]
                   + f_3 * pc_x[k] * ksk_683[k];

        t_846[k] = f_16 * isk_496[k]
                   + f_1 * ksi0_525[k]
                   - f_2 * ksi1_525[k]
                   + f_3 * pc_y[k] * ksk_676[k];

        t_847[k] = f_17 * isk_460[k]
                   + f_3 * pc_z[k] * ksk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pc_y, isk_498, isk_499, isk_500, ksi0_527, \
                         ksi0_528, ksi0_529, ksi1_527, ksi1_528, ksi1_529, ksk_678, ksk_679, \
                         ksk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * isk_498[k]
                   + f_12 * ksi0_527[k]
                   - f_13 * ksi1_527[k]
                   + f_3 * pc_y[k] * ksk_678[k];

        t_849[k] = f_16 * isk_499[k]
                   + f_10 * ksi0_528[k]
                   - f_11 * ksi1_528[k]
                   + f_3 * pc_y[k] * ksk_679[k];

        t_850[k] = f_16 * isk_500[k]
                   + f_8 * ksi0_529[k]
                   - f_9 * ksi1_529[k]
                   + f_3 * pc_y[k] * ksk_680[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pc_y, isk_501, isk_502, isk_503, ksi0_530, \
                         ksi0_531, ksi1_530, ksi1_531, ksk_681, ksk_682, \
                         ksk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_16 * isk_501[k]
                   + f_6 * ksi0_530[k]
                   - f_7 * ksi1_530[k]
                   + f_3 * pc_y[k] * ksk_681[k];

        t_852[k] = f_16 * isk_502[k]
                   + f_4 * ksi0_531[k]
                   - f_5 * ksi1_531[k]
                   + f_3 * pc_y[k] * ksk_682[k];

        t_853[k] = f_16 * isk_503[k]
                   + f_3 * pc_y[k] * ksk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pa_y, pc_y, pc_z, isl0_630, isk_467, \
                         isk_468, isk_504, isl1_630, ksi0_531, ksi1_531, ksk_683, \
                         ksk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_17 * isk_467[k]
                   + f_1 * ksi0_531[k]
                   - f_2 * ksi1_531[k]
                   + f_3 * pc_z[k] * ksk_683[k];

        t_855[k] = pa_y[k] * isl0_630[k]
                   - f_14 * pc_y[k] * isl1_630[k];

        t_856[k] = f_15 * isk_504[k]
                   + f_3 * pc_y[k] * ksk_684[k];

        t_857[k] = f_18 * isk_468[k]
                   + f_3 * pc_z[k] * ksk_684[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pa_y, pc_y, isl0_633, isl0_635, isl0_636, \
                         isk_505, isk_506, isk_507, isl1_633, isl1_635, isl1_636, \
                         ksk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pa_y[k] * isl0_633[k]
                   + f_16 * isk_505[k]
                   - f_14 * pc_y[k] * isl1_633[k];

        t_859[k] = f_15 * isk_506[k]
                   + f_3 * pc_y[k] * ksk_686[k];

        t_860[k] = pa_y[k] * isl0_635[k]
                   - f_14 * pc_y[k] * isl1_635[k];

        t_861[k] = pa_y[k] * isl0_636[k]
                   + f_17 * isk_507[k]
                   - f_14 * pc_y[k] * isl1_636[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pa_y, pc_y, pc_z, isl0_639, isl0_640, \
                         isk_471, isk_509, isk_510, isl1_639, isl1_640, ksk_687, \
                         ksk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_18 * isk_471[k]
                   + f_3 * pc_z[k] * ksk_687[k];

        t_863[k] = f_15 * isk_509[k]
                   + f_3 * pc_y[k] * ksk_689[k];

        t_864[k] = pa_y[k] * isl0_639[k]
                   - f_14 * pc_y[k] * isl1_639[k];

        t_865[k] = pa_y[k] * isl0_640[k]
                   + f_18 * isk_510[k]
                   - f_14 * pc_y[k] * isl1_640[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pa_y, pc_y, pc_z, isl0_642, isl0_644, \
                         isk_474, isk_512, isk_513, isl1_642, isl1_644, ksk_690, \
                         ksk_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * isk_474[k]
                   + f_3 * pc_z[k] * ksk_690[k];

        t_867[k] = pa_y[k] * isl0_642[k]
                   + f_16 * isk_512[k]
                   - f_14 * pc_y[k] * isl1_642[k];

        t_868[k] = f_15 * isk_513[k]
                   + f_3 * pc_y[k] * ksk_693[k];

        t_869[k] = pa_y[k] * isl0_644[k]
                   - f_14 * pc_y[k] * isl1_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pa_y, pc_y, pc_z, isl0_645, isl0_647, isk_478, \
                         isk_514, isk_516, isl1_645, isl1_647, \
                         ksk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pa_y[k] * isl0_645[k]
                   + f_19 * isk_514[k]
                   - f_14 * pc_y[k] * isl1_645[k];

        t_871[k] = f_18 * isk_478[k]
                   + f_3 * pc_z[k] * ksk_694[k];

        t_872[k] = pa_y[k] * isl0_647[k]
                   + f_17 * isk_516[k]
                   - f_14 * pc_y[k] * isl1_647[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pa_y, pc_y, isl0_648, isl0_650, isl0_651, \
                         isk_517, isk_518, isk_519, isl1_648, isl1_650, isl1_651, \
                         ksk_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pa_y[k] * isl0_648[k]
                   + f_16 * isk_517[k]
                   - f_14 * pc_y[k] * isl1_648[k];

        t_874[k] = f_15 * isk_518[k]
                   + f_3 * pc_y[k] * ksk_698[k];

        t_875[k] = pa_y[k] * isl0_650[k]
                   - f_14 * pc_y[k] * isl1_650[k];

        t_876[k] = pa_y[k] * isl0_651[k]
                   + f_20 * isk_519[k]
                   - f_14 * pc_y[k] * isl1_651[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pa_y, pc_y, pc_z, isl0_653, isl0_654, isk_483, \
                         isk_521, isk_522, isl1_653, isl1_654, \
                         ksk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_18 * isk_483[k]
                   + f_3 * pc_z[k] * ksk_699[k];

        t_878[k] = pa_y[k] * isl0_653[k]
                   + f_18 * isk_521[k]
                   - f_14 * pc_y[k] * isl1_653[k];

        t_879[k] = pa_y[k] * isl0_654[k]
                   + f_17 * isk_522[k]
                   - f_14 * pc_y[k] * isl1_654[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pa_y, pc_x, pc_y, isl0_655, isl0_657, \
                         isk_523, isk_524, isk_712, isl1_655, isl1_657, ksk_704, \
                         ksk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = pa_y[k] * isl0_655[k]
                   + f_16 * isk_523[k]
                   - f_14 * pc_y[k] * isl1_655[k];

        t_881[k] = f_15 * isk_524[k]
                   + f_3 * pc_y[k] * ksk_704[k];

        t_882[k] = pa_y[k] * isl0_657[k]
                   - f_14 * pc_y[k] * isl1_657[k];

        t_883[k] = f_16 * isk_712[k]
                   + f_3 * pc_x[k] * ksk_712[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, pc_x, isk_713, isk_714, isk_715, \
                         isk_716, isk_717, ksk_713, ksk_714, ksk_715, ksk_716, \
                         ksk_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_16 * isk_713[k]
                   + f_3 * pc_x[k] * ksk_713[k];

        t_885[k] = f_16 * isk_714[k]
                   + f_3 * pc_x[k] * ksk_714[k];

        t_886[k] = f_16 * isk_715[k]
                   + f_3 * pc_x[k] * ksk_715[k];

        t_887[k] = f_16 * isk_716[k]
                   + f_3 * pc_x[k] * ksk_716[k];

        t_888[k] = f_16 * isk_717[k]
                   + f_3 * pc_x[k] * ksk_717[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, t_892, pc_x, pc_y, pc_z, isk_496, isk_532, \
                         isk_718, isk_719, ksi0_553, ksi1_553, ksk_712, ksk_718, \
                         ksk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_16 * isk_718[k]
                   + f_3 * pc_x[k] * ksk_718[k];

        t_890[k] = f_16 * isk_719[k]
                   + f_3 * pc_x[k] * ksk_719[k];

        t_891[k] = f_15 * isk_532[k]
                   + f_1 * ksi0_553[k]
                   - f_2 * ksi1_553[k]
                   + f_3 * pc_y[k] * ksk_712[k];

        t_892[k] = f_18 * isk_496[k]
                   + f_3 * pc_z[k] * ksk_712[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, isk_534, isk_535, isk_536, ksi0_555, \
                         ksi0_556, ksi0_557, ksi1_555, ksi1_556, ksi1_557, ksk_714, ksk_715, \
                         ksk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_15 * isk_534[k]
                   + f_12 * ksi0_555[k]
                   - f_13 * ksi1_555[k]
                   + f_3 * pc_y[k] * ksk_714[k];

        t_894[k] = f_15 * isk_535[k]
                   + f_10 * ksi0_556[k]
                   - f_11 * ksi1_556[k]
                   + f_3 * pc_y[k] * ksk_715[k];

        t_895[k] = f_15 * isk_536[k]
                   + f_8 * ksi0_557[k]
                   - f_9 * ksi1_557[k]
                   + f_3 * pc_y[k] * ksk_716[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_y, isk_537, isk_538, isk_539, ksi0_558, \
                         ksi0_559, ksi1_558, ksi1_559, ksk_717, ksk_718, \
                         ksk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_15 * isk_537[k]
                   + f_6 * ksi0_558[k]
                   - f_7 * ksi1_558[k]
                   + f_3 * pc_y[k] * ksk_717[k];

        t_897[k] = f_15 * isk_538[k]
                   + f_4 * ksi0_559[k]
                   - f_5 * ksi1_559[k]
                   + f_3 * pc_y[k] * ksk_718[k];

        t_898[k] = f_15 * isk_539[k]
                   + f_3 * pc_y[k] * ksk_719[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pa_y, pc_x, pc_y, pc_z, isl0_674, \
                         isk_504, isk_720, isl1_674, ksi0_560, ksi1_560, \
                         ksk_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pa_y[k] * isl0_674[k]
                   - f_14 * pc_y[k] * isl1_674[k];

        t_900[k] = f_16 * isk_720[k]
                   + f_1 * ksi0_560[k]
                   - f_2 * ksi1_560[k]
                   + f_3 * pc_x[k] * ksk_720[k];

        t_901[k] = f_3 * pc_y[k] * ksk_720[k];

        t_902[k] = f_19 * isk_504[k]
                   + f_3 * pc_z[k] * ksk_720[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pc_x, pc_y, isk_725, ksi0_560, ksi0_565, \
                         ksi1_560, ksi1_565, ksk_721, ksk_722, \
                         ksk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_4 * ksi0_560[k]
                   - f_5 * ksi1_560[k]
                   + f_3 * pc_y[k] * ksk_721[k];

        t_904[k] = f_3 * pc_y[k] * ksk_722[k];

        t_905[k] = f_16 * isk_725[k]
                   + f_12 * ksi0_565[k]
                   - f_13 * ksi1_565[k]
                   + f_3 * pc_x[k] * ksk_725[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pc_y, ksi0_561, ksi0_562, ksi1_561, ksi1_562, \
                         ksk_723, ksk_724, ksk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_6 * ksi0_561[k]
                   - f_7 * ksi1_561[k]
                   + f_3 * pc_y[k] * ksk_723[k];

        t_907[k] = f_4 * ksi0_562[k]
                   - f_5 * ksi1_562[k]
                   + f_3 * pc_y[k] * ksk_724[k];

        t_908[k] = f_3 * pc_y[k] * ksk_725[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, pc_y, isk_729, ksi0_563, ksi0_564, \
                         ksi0_569, ksi1_563, ksi1_564, ksi1_569, ksk_726, ksk_727, \
                         ksk_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_16 * isk_729[k]
                   + f_10 * ksi0_569[k]
                   - f_11 * ksi1_569[k]
                   + f_3 * pc_x[k] * ksk_729[k];

        t_910[k] = f_8 * ksi0_563[k]
                   - f_9 * ksi1_563[k]
                   + f_3 * pc_y[k] * ksk_726[k];

        t_911[k] = f_6 * ksi0_564[k]
                   - f_7 * ksi1_564[k]
                   + f_3 * pc_y[k] * ksk_727[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_x, pc_y, isk_734, ksi0_565, ksi0_574, \
                         ksi1_565, ksi1_574, ksk_728, ksk_729, \
                         ksk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_4 * ksi0_565[k]
                   - f_5 * ksi1_565[k]
                   + f_3 * pc_y[k] * ksk_728[k];

        t_913[k] = f_3 * pc_y[k] * ksk_729[k];

        t_914[k] = f_16 * isk_734[k]
                   + f_8 * ksi0_574[k]
                   - f_9 * ksi1_574[k]
                   + f_3 * pc_x[k] * ksk_734[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_y, ksi0_566, ksi0_567, ksi0_568, ksi1_566, \
                         ksi1_567, ksi1_568, ksk_730, ksk_731, \
                         ksk_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_10 * ksi0_566[k]
                   - f_11 * ksi1_566[k]
                   + f_3 * pc_y[k] * ksk_730[k];

        t_916[k] = f_8 * ksi0_567[k]
                   - f_9 * ksi1_567[k]
                   + f_3 * pc_y[k] * ksk_731[k];

        t_917[k] = f_6 * ksi0_568[k]
                   - f_7 * ksi1_568[k]
                   + f_3 * pc_y[k] * ksk_732[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_x, pc_y, isk_740, ksi0_569, ksi0_580, \
                         ksi1_569, ksi1_580, ksk_733, ksk_734, \
                         ksk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_4 * ksi0_569[k]
                   - f_5 * ksi1_569[k]
                   + f_3 * pc_y[k] * ksk_733[k];

        t_919[k] = f_3 * pc_y[k] * ksk_734[k];

        t_920[k] = f_16 * isk_740[k]
                   + f_6 * ksi0_580[k]
                   - f_7 * ksi1_580[k]
                   + f_3 * pc_x[k] * ksk_740[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksi0, const size_t ksi1,
                                                          const size_t ksk, const size_t ncols,
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
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);
    const auto f_23 = 4.0 / q;

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
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isl0_675 = buffer.data(isl0 + 675);
    const auto *isl0_678 = buffer.data(isl0 + 678);
    const auto *isl0_681 = buffer.data(isl0 + 681);
    const auto *isl0_685 = buffer.data(isl0 + 685);
    const auto *isl0_690 = buffer.data(isl0 + 690);
    const auto *isl0_696 = buffer.data(isl0 + 696);
    const auto *isl0_945 = buffer.data(isl0 + 945);
    const auto *isl0_948 = buffer.data(isl0 + 948);
    const auto *isl0_951 = buffer.data(isl0 + 951);
    const auto *isl0_955 = buffer.data(isl0 + 955);
    const auto *isl0_960 = buffer.data(isl0 + 960);
    const auto *isl0_966 = buffer.data(isl0 + 966);
    const auto *isl0_981 = buffer.data(isl0 + 981);
    const auto *isl0_983 = buffer.data(isl0 + 983);
    const auto *isl0_984 = buffer.data(isl0 + 984);
    const auto *isl0_985 = buffer.data(isl0 + 985);
    const auto *isl0_986 = buffer.data(isl0 + 986);
    const auto *isl0_987 = buffer.data(isl0 + 987);
    const auto *isl0_989 = buffer.data(isl0 + 989);
    const auto *isl0_995 = buffer.data(isl0 + 995);
    const auto *isl0_999 = buffer.data(isl0 + 999);
    const auto *isl0_1002 = buffer.data(isl0 + 1002);
    const auto *isl0_1004 = buffer.data(isl0 + 1004);
    const auto *isl0_1007 = buffer.data(isl0 + 1007);
    const auto *isl0_1008 = buffer.data(isl0 + 1008);
    const auto *isl0_1010 = buffer.data(isl0 + 1010);
    const auto *isl0_1013 = buffer.data(isl0 + 1013);
    const auto *isl0_1014 = buffer.data(isl0 + 1014);
    const auto *isl0_1015 = buffer.data(isl0 + 1015);
    const auto *isl0_1017 = buffer.data(isl0 + 1017);
    const auto *isl0_1026 = buffer.data(isl0 + 1026);
    const auto *isl0_1028 = buffer.data(isl0 + 1028);
    const auto *isl0_1029 = buffer.data(isl0 + 1029);
    const auto *isl0_1030 = buffer.data(isl0 + 1030);
    const auto *isl0_1031 = buffer.data(isl0 + 1031);
    const auto *isl0_1032 = buffer.data(isl0 + 1032);
    const auto *isl0_1034 = buffer.data(isl0 + 1034);
    const auto *isl0_1035 = buffer.data(isl0 + 1035);
    const auto *isl0_1038 = buffer.data(isl0 + 1038);
    const auto *isl0_1040 = buffer.data(isl0 + 1040);

    const auto *isk_539 = buffer.data(isk + 539);
    const auto *isk_540 = buffer.data(isk + 540);
    const auto *isk_543 = buffer.data(isk + 543);
    const auto *isk_545 = buffer.data(isk + 545);
    const auto *isk_546 = buffer.data(isk + 546);
    const auto *isk_549 = buffer.data(isk + 549);
    const auto *isk_550 = buffer.data(isk + 550);
    const auto *isk_554 = buffer.data(isk + 554);
    const auto *isk_555 = buffer.data(isk + 555);
    const auto *isk_560 = buffer.data(isk + 560);
    const auto *isk_568 = buffer.data(isk + 568);
    const auto *isk_575 = buffer.data(isk + 575);
    const auto *isk_576 = buffer.data(isk + 576);
    const auto *isk_578 = buffer.data(isk + 578);
    const auto *isk_581 = buffer.data(isk + 581);
    const auto *isk_585 = buffer.data(isk + 585);
    const auto *isk_590 = buffer.data(isk + 590);
    const auto *isk_596 = buffer.data(isk + 596);
    const auto *isk_611 = buffer.data(isk + 611);
    const auto *isk_612 = buffer.data(isk + 612);
    const auto *isk_614 = buffer.data(isk + 614);
    const auto *isk_747 = buffer.data(isk + 747);
    const auto *isk_748 = buffer.data(isk + 748);
    const auto *isk_749 = buffer.data(isk + 749);
    const auto *isk_750 = buffer.data(isk + 750);
    const auto *isk_751 = buffer.data(isk + 751);
    const auto *isk_752 = buffer.data(isk + 752);
    const auto *isk_753 = buffer.data(isk + 753);
    const auto *isk_755 = buffer.data(isk + 755);
    const auto *isk_756 = buffer.data(isk + 756);
    const auto *isk_759 = buffer.data(isk + 759);
    const auto *isk_762 = buffer.data(isk + 762);
    const auto *isk_766 = buffer.data(isk + 766);
    const auto *isk_771 = buffer.data(isk + 771);
    const auto *isk_777 = buffer.data(isk + 777);
    const auto *isk_784 = buffer.data(isk + 784);
    const auto *isk_786 = buffer.data(isk + 786);
    const auto *isk_787 = buffer.data(isk + 787);
    const auto *isk_788 = buffer.data(isk + 788);
    const auto *isk_789 = buffer.data(isk + 789);
    const auto *isk_790 = buffer.data(isk + 790);
    const auto *isk_791 = buffer.data(isk + 791);
    const auto *isk_797 = buffer.data(isk + 797);
    const auto *isk_801 = buffer.data(isk + 801);
    const auto *isk_804 = buffer.data(isk + 804);
    const auto *isk_806 = buffer.data(isk + 806);
    const auto *isk_809 = buffer.data(isk + 809);
    const auto *isk_810 = buffer.data(isk + 810);
    const auto *isk_812 = buffer.data(isk + 812);
    const auto *isk_815 = buffer.data(isk + 815);
    const auto *isk_816 = buffer.data(isk + 816);
    const auto *isk_817 = buffer.data(isk + 817);
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
    const auto *isk_831 = buffer.data(isk + 831);
    const auto *isk_833 = buffer.data(isk + 833);

    const auto *isl1_675 = buffer.data(isl1 + 675);
    const auto *isl1_678 = buffer.data(isl1 + 678);
    const auto *isl1_681 = buffer.data(isl1 + 681);
    const auto *isl1_685 = buffer.data(isl1 + 685);
    const auto *isl1_690 = buffer.data(isl1 + 690);
    const auto *isl1_696 = buffer.data(isl1 + 696);
    const auto *isl1_945 = buffer.data(isl1 + 945);
    const auto *isl1_948 = buffer.data(isl1 + 948);
    const auto *isl1_951 = buffer.data(isl1 + 951);
    const auto *isl1_955 = buffer.data(isl1 + 955);
    const auto *isl1_960 = buffer.data(isl1 + 960);
    const auto *isl1_966 = buffer.data(isl1 + 966);
    const auto *isl1_981 = buffer.data(isl1 + 981);
    const auto *isl1_983 = buffer.data(isl1 + 983);
    const auto *isl1_984 = buffer.data(isl1 + 984);
    const auto *isl1_985 = buffer.data(isl1 + 985);
    const auto *isl1_986 = buffer.data(isl1 + 986);
    const auto *isl1_987 = buffer.data(isl1 + 987);
    const auto *isl1_989 = buffer.data(isl1 + 989);
    const auto *isl1_995 = buffer.data(isl1 + 995);
    const auto *isl1_999 = buffer.data(isl1 + 999);
    const auto *isl1_1002 = buffer.data(isl1 + 1002);
    const auto *isl1_1004 = buffer.data(isl1 + 1004);
    const auto *isl1_1007 = buffer.data(isl1 + 1007);
    const auto *isl1_1008 = buffer.data(isl1 + 1008);
    const auto *isl1_1010 = buffer.data(isl1 + 1010);
    const auto *isl1_1013 = buffer.data(isl1 + 1013);
    const auto *isl1_1014 = buffer.data(isl1 + 1014);
    const auto *isl1_1015 = buffer.data(isl1 + 1015);
    const auto *isl1_1017 = buffer.data(isl1 + 1017);
    const auto *isl1_1026 = buffer.data(isl1 + 1026);
    const auto *isl1_1028 = buffer.data(isl1 + 1028);
    const auto *isl1_1029 = buffer.data(isl1 + 1029);
    const auto *isl1_1030 = buffer.data(isl1 + 1030);
    const auto *isl1_1031 = buffer.data(isl1 + 1031);
    const auto *isl1_1032 = buffer.data(isl1 + 1032);
    const auto *isl1_1034 = buffer.data(isl1 + 1034);
    const auto *isl1_1035 = buffer.data(isl1 + 1035);
    const auto *isl1_1038 = buffer.data(isl1 + 1038);
    const auto *isl1_1040 = buffer.data(isl1 + 1040);

    const auto *ksi0_570 = buffer.data(ksi0 + 570);
    const auto *ksi0_571 = buffer.data(ksi0 + 571);
    const auto *ksi0_572 = buffer.data(ksi0 + 572);
    const auto *ksi0_573 = buffer.data(ksi0 + 573);
    const auto *ksi0_574 = buffer.data(ksi0 + 574);
    const auto *ksi0_581 = buffer.data(ksi0 + 581);
    const auto *ksi0_582 = buffer.data(ksi0 + 582);
    const auto *ksi0_583 = buffer.data(ksi0 + 583);
    const auto *ksi0_584 = buffer.data(ksi0 + 584);
    const auto *ksi0_585 = buffer.data(ksi0 + 585);
    const auto *ksi0_586 = buffer.data(ksi0 + 586);
    const auto *ksi0_587 = buffer.data(ksi0 + 587);
    const auto *ksi0_588 = buffer.data(ksi0 + 588);
    const auto *ksi0_590 = buffer.data(ksi0 + 590);
    const auto *ksi0_591 = buffer.data(ksi0 + 591);
    const auto *ksi0_593 = buffer.data(ksi0 + 593);
    const auto *ksi0_594 = buffer.data(ksi0 + 594);
    const auto *ksi0_595 = buffer.data(ksi0 + 595);
    const auto *ksi0_597 = buffer.data(ksi0 + 597);
    const auto *ksi0_598 = buffer.data(ksi0 + 598);
    const auto *ksi0_599 = buffer.data(ksi0 + 599);
    const auto *ksi0_600 = buffer.data(ksi0 + 600);
    const auto *ksi0_602 = buffer.data(ksi0 + 602);

    const auto *ksi1_570 = buffer.data(ksi1 + 570);
    const auto *ksi1_571 = buffer.data(ksi1 + 571);
    const auto *ksi1_572 = buffer.data(ksi1 + 572);
    const auto *ksi1_573 = buffer.data(ksi1 + 573);
    const auto *ksi1_574 = buffer.data(ksi1 + 574);
    const auto *ksi1_581 = buffer.data(ksi1 + 581);
    const auto *ksi1_582 = buffer.data(ksi1 + 582);
    const auto *ksi1_583 = buffer.data(ksi1 + 583);
    const auto *ksi1_584 = buffer.data(ksi1 + 584);
    const auto *ksi1_585 = buffer.data(ksi1 + 585);
    const auto *ksi1_586 = buffer.data(ksi1 + 586);
    const auto *ksi1_587 = buffer.data(ksi1 + 587);
    const auto *ksi1_588 = buffer.data(ksi1 + 588);
    const auto *ksi1_590 = buffer.data(ksi1 + 590);
    const auto *ksi1_591 = buffer.data(ksi1 + 591);
    const auto *ksi1_593 = buffer.data(ksi1 + 593);
    const auto *ksi1_594 = buffer.data(ksi1 + 594);
    const auto *ksi1_595 = buffer.data(ksi1 + 595);
    const auto *ksi1_597 = buffer.data(ksi1 + 597);
    const auto *ksi1_598 = buffer.data(ksi1 + 598);
    const auto *ksi1_599 = buffer.data(ksi1 + 599);
    const auto *ksi1_600 = buffer.data(ksi1 + 600);
    const auto *ksi1_602 = buffer.data(ksi1 + 602);

    const auto *ksk_735 = buffer.data(ksk + 735);
    const auto *ksk_736 = buffer.data(ksk + 736);
    const auto *ksk_737 = buffer.data(ksk + 737);
    const auto *ksk_738 = buffer.data(ksk + 738);
    const auto *ksk_739 = buffer.data(ksk + 739);
    const auto *ksk_740 = buffer.data(ksk + 740);
    const auto *ksk_747 = buffer.data(ksk + 747);
    const auto *ksk_748 = buffer.data(ksk + 748);
    const auto *ksk_749 = buffer.data(ksk + 749);
    const auto *ksk_750 = buffer.data(ksk + 750);
    const auto *ksk_751 = buffer.data(ksk + 751);
    const auto *ksk_752 = buffer.data(ksk + 752);
    const auto *ksk_753 = buffer.data(ksk + 753);
    const auto *ksk_754 = buffer.data(ksk + 754);
    const auto *ksk_755 = buffer.data(ksk + 755);
    const auto *ksk_756 = buffer.data(ksk + 756);
    const auto *ksk_757 = buffer.data(ksk + 757);
    const auto *ksk_758 = buffer.data(ksk + 758);
    const auto *ksk_759 = buffer.data(ksk + 759);
    const auto *ksk_761 = buffer.data(ksk + 761);
    const auto *ksk_762 = buffer.data(ksk + 762);
    const auto *ksk_763 = buffer.data(ksk + 763);
    const auto *ksk_765 = buffer.data(ksk + 765);
    const auto *ksk_766 = buffer.data(ksk + 766);
    const auto *ksk_767 = buffer.data(ksk + 767);
    const auto *ksk_768 = buffer.data(ksk + 768);
    const auto *ksk_770 = buffer.data(ksk + 770);
    const auto *ksk_771 = buffer.data(ksk + 771);
    const auto *ksk_772 = buffer.data(ksk + 772);
    const auto *ksk_773 = buffer.data(ksk + 773);
    const auto *ksk_774 = buffer.data(ksk + 774);
    const auto *ksk_776 = buffer.data(ksk + 776);
    const auto *ksk_777 = buffer.data(ksk + 777);
    const auto *ksk_784 = buffer.data(ksk + 784);
    const auto *ksk_786 = buffer.data(ksk + 786);
    const auto *ksk_787 = buffer.data(ksk + 787);
    const auto *ksk_788 = buffer.data(ksk + 788);
    const auto *ksk_789 = buffer.data(ksk + 789);
    const auto *ksk_790 = buffer.data(ksk + 790);
    const auto *ksk_791 = buffer.data(ksk + 791);
    const auto *ksk_792 = buffer.data(ksk + 792);
    const auto *ksk_794 = buffer.data(ksk + 794);
    const auto *ksk_795 = buffer.data(ksk + 795);
    const auto *ksk_797 = buffer.data(ksk + 797);
    const auto *ksk_798 = buffer.data(ksk + 798);
    const auto *ksk_801 = buffer.data(ksk + 801);
    const auto *ksk_802 = buffer.data(ksk + 802);
    const auto *ksk_806 = buffer.data(ksk + 806);
    const auto *ksk_807 = buffer.data(ksk + 807);
    const auto *ksk_812 = buffer.data(ksk + 812);
    const auto *ksk_820 = buffer.data(ksk + 820);
    const auto *ksk_821 = buffer.data(ksk + 821);
    const auto *ksk_822 = buffer.data(ksk + 822);
    const auto *ksk_823 = buffer.data(ksk + 823);
    const auto *ksk_824 = buffer.data(ksk + 824);
    const auto *ksk_825 = buffer.data(ksk + 825);
    const auto *ksk_826 = buffer.data(ksk + 826);
    const auto *ksk_827 = buffer.data(ksk + 827);
    const auto *ksk_828 = buffer.data(ksk + 828);
    const auto *ksk_830 = buffer.data(ksk + 830);

#pragma omp simd aligned(t_921, t_922, t_923, pc_y, ksi0_570, ksi0_571, ksi0_572, ksi1_570, \
                         ksi1_571, ksi1_572, ksk_735, ksk_736, \
                         ksk_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_12 * ksi0_570[k]
                   - f_13 * ksi1_570[k]
                   + f_3 * pc_y[k] * ksk_735[k];

        t_922[k] = f_10 * ksi0_571[k]
                   - f_11 * ksi1_571[k]
                   + f_3 * pc_y[k] * ksk_736[k];

        t_923[k] = f_8 * ksi0_572[k]
                   - f_9 * ksi1_572[k]
                   + f_3 * pc_y[k] * ksk_737[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_y, ksi0_573, ksi0_574, ksi1_573, ksi1_574, \
                         ksk_738, ksk_739, ksk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_6 * ksi0_573[k]
                   - f_7 * ksi1_573[k]
                   + f_3 * pc_y[k] * ksk_738[k];

        t_925[k] = f_4 * ksi0_574[k]
                   - f_5 * ksi1_574[k]
                   + f_3 * pc_y[k] * ksk_739[k];

        t_926[k] = f_3 * pc_y[k] * ksk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pc_x, isk_747, isk_748, isk_749, isk_750, \
                         ksi0_587, ksi1_587, ksk_747, ksk_748, ksk_749, \
                         ksk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_16 * isk_747[k]
                   + f_4 * ksi0_587[k]
                   - f_5 * ksi1_587[k]
                   + f_3 * pc_x[k] * ksk_747[k];

        t_928[k] = f_16 * isk_748[k]
                   + f_3 * pc_x[k] * ksk_748[k];

        t_929[k] = f_16 * isk_749[k]
                   + f_3 * pc_x[k] * ksk_749[k];

        t_930[k] = f_16 * isk_750[k]
                   + f_3 * pc_x[k] * ksk_750[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pc_x, pc_y, isk_751, isk_752, \
                         isk_753, isk_755, ksk_747, ksk_751, ksk_752, ksk_753, \
                         ksk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_16 * isk_751[k]
                   + f_3 * pc_x[k] * ksk_751[k];

        t_932[k] = f_16 * isk_752[k]
                   + f_3 * pc_x[k] * ksk_752[k];

        t_933[k] = f_16 * isk_753[k]
                   + f_3 * pc_x[k] * ksk_753[k];

        t_934[k] = f_3 * pc_y[k] * ksk_747[k];

        t_935[k] = f_16 * isk_755[k]
                   + f_3 * pc_x[k] * ksk_755[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, pc_y, ksi0_581, ksi0_582, ksi0_583, ksi1_581, \
                         ksi1_582, ksi1_583, ksk_748, ksk_749, \
                         ksk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_1 * ksi0_581[k]
                   - f_2 * ksi1_581[k]
                   + f_3 * pc_y[k] * ksk_748[k];

        t_937[k] = f_21 * ksi0_582[k]
                   - f_22 * ksi1_582[k]
                   + f_3 * pc_y[k] * ksk_749[k];

        t_938[k] = f_12 * ksi0_583[k]
                   - f_13 * ksi1_583[k]
                   + f_3 * pc_y[k] * ksk_750[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, pc_y, ksi0_584, ksi0_585, ksi0_586, ksi1_584, \
                         ksi1_585, ksi1_586, ksk_751, ksk_752, \
                         ksk_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_10 * ksi0_584[k]
                   - f_11 * ksi1_584[k]
                   + f_3 * pc_y[k] * ksk_751[k];

        t_940[k] = f_8 * ksi0_585[k]
                   - f_9 * ksi1_585[k]
                   + f_3 * pc_y[k] * ksk_752[k];

        t_941[k] = f_6 * ksi0_586[k]
                   - f_7 * ksi1_586[k]
                   + f_3 * pc_y[k] * ksk_753[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pa_x, pc_x, pc_y, pc_z, isl0_945, \
                         isk_539, isk_756, isl1_945, ksi0_587, ksi1_587, ksk_754, \
                         ksk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_4 * ksi0_587[k]
                   - f_5 * ksi1_587[k]
                   + f_3 * pc_y[k] * ksk_754[k];

        t_943[k] = f_3 * pc_y[k] * ksk_755[k];

        t_944[k] = f_19 * isk_539[k]
                   + f_1 * ksi0_587[k]
                   - f_2 * ksi1_587[k]
                   + f_3 * pc_z[k] * ksk_755[k];

        t_945[k] = pa_x[k] * isl0_945[k]
                   + f_23 * isk_756[k]
                   - f_14 * pc_x[k] * isl1_945[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, t_949, pa_x, pc_x, pc_y, pc_z, isl0_948, \
                         isk_540, isk_759, isl1_948, ksk_756, ksk_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_20 * isk_540[k]
                   + f_3 * pc_y[k] * ksk_756[k];

        t_947[k] = f_3 * pc_z[k] * ksk_756[k];

        t_948[k] = pa_x[k] * isl0_948[k]
                   + f_20 * isk_759[k]
                   - f_14 * pc_x[k] * isl1_948[k];

        t_949[k] = f_3 * pc_z[k] * ksk_757[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pa_x, pc_x, pc_z, isl0_951, isk_762, isl1_951, \
                         ksi0_588, ksi1_588, ksk_758, ksk_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_4 * ksi0_588[k]
                   - f_5 * ksi1_588[k]
                   + f_3 * pc_z[k] * ksk_758[k];

        t_951[k] = pa_x[k] * isl0_951[k]
                   + f_19 * isk_762[k]
                   - f_14 * pc_x[k] * isl1_951[k];

        t_952[k] = f_3 * pc_z[k] * ksk_759[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pa_x, pc_x, pc_y, pc_z, isl0_955, \
                         isk_545, isk_766, isl1_955, ksi0_590, ksi1_590, ksk_761, \
                         ksk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_20 * isk_545[k]
                   + f_3 * pc_y[k] * ksk_761[k];

        t_954[k] = f_6 * ksi0_590[k]
                   - f_7 * ksi1_590[k]
                   + f_3 * pc_z[k] * ksk_761[k];

        t_955[k] = pa_x[k] * isl0_955[k]
                   + f_18 * isk_766[k]
                   - f_14 * pc_x[k] * isl1_955[k];

        t_956[k] = f_3 * pc_z[k] * ksk_762[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, pc_y, pc_z, isk_549, ksi0_591, ksi0_593, \
                         ksi1_591, ksi1_593, ksk_763, ksk_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_4 * ksi0_591[k]
                   - f_5 * ksi1_591[k]
                   + f_3 * pc_z[k] * ksk_763[k];

        t_958[k] = f_20 * isk_549[k]
                   + f_3 * pc_y[k] * ksk_765[k];

        t_959[k] = f_8 * ksi0_593[k]
                   - f_9 * ksi1_593[k]
                   + f_3 * pc_z[k] * ksk_765[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, pa_x, pc_x, pc_z, isl0_960, isk_771, isl1_960, \
                         ksi0_594, ksi1_594, ksk_766, ksk_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = pa_x[k] * isl0_960[k]
                   + f_17 * isk_771[k]
                   - f_14 * pc_x[k] * isl1_960[k];

        t_961[k] = f_3 * pc_z[k] * ksk_766[k];

        t_962[k] = f_4 * ksi0_594[k]
                   - f_5 * ksi1_594[k]
                   + f_3 * pc_z[k] * ksk_767[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, pc_y, pc_z, isk_554, ksi0_595, ksi0_597, \
                         ksi1_595, ksi1_597, ksk_768, ksk_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = f_6 * ksi0_595[k]
                   - f_7 * ksi1_595[k]
                   + f_3 * pc_z[k] * ksk_768[k];

        t_964[k] = f_20 * isk_554[k]
                   + f_3 * pc_y[k] * ksk_770[k];

        t_965[k] = f_10 * ksi0_597[k]
                   - f_11 * ksi1_597[k]
                   + f_3 * pc_z[k] * ksk_770[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pa_x, pc_x, pc_z, isl0_966, isk_777, isl1_966, \
                         ksi0_598, ksi1_598, ksk_771, ksk_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = pa_x[k] * isl0_966[k]
                   + f_16 * isk_777[k]
                   - f_14 * pc_x[k] * isl1_966[k];

        t_967[k] = f_3 * pc_z[k] * ksk_771[k];

        t_968[k] = f_4 * ksi0_598[k]
                   - f_5 * ksi1_598[k]
                   + f_3 * pc_z[k] * ksk_772[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, pc_y, pc_z, isk_560, ksi0_599, ksi0_600, \
                         ksi0_602, ksi1_599, ksi1_600, ksi1_602, ksk_773, ksk_774, \
                         ksk_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_6 * ksi0_599[k]
                   - f_7 * ksi1_599[k]
                   + f_3 * pc_z[k] * ksk_773[k];

        t_970[k] = f_8 * ksi0_600[k]
                   - f_9 * ksi1_600[k]
                   + f_3 * pc_z[k] * ksk_774[k];

        t_971[k] = f_20 * isk_560[k]
                   + f_3 * pc_y[k] * ksk_776[k];

        t_972[k] = f_12 * ksi0_602[k]
                   - f_13 * ksi1_602[k]
                   + f_3 * pc_z[k] * ksk_776[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, t_977, pc_x, pc_z, isk_784, isk_786, \
                         isk_787, isk_788, ksk_777, ksk_784, ksk_786, ksk_787, \
                         ksk_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_15 * isk_784[k]
                   + f_3 * pc_x[k] * ksk_784[k];

        t_974[k] = f_3 * pc_z[k] * ksk_777[k];

        t_975[k] = f_15 * isk_786[k]
                   + f_3 * pc_x[k] * ksk_786[k];

        t_976[k] = f_15 * isk_787[k]
                   + f_3 * pc_x[k] * ksk_787[k];

        t_977[k] = f_15 * isk_788[k]
                   + f_3 * pc_x[k] * ksk_788[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, pa_x, pc_x, isl0_981, isk_789, isk_790, \
                         isk_791, isl1_981, ksk_789, ksk_790, ksk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_15 * isk_789[k]
                   + f_3 * pc_x[k] * ksk_789[k];

        t_979[k] = f_15 * isk_790[k]
                   + f_3 * pc_x[k] * ksk_790[k];

        t_980[k] = f_15 * isk_791[k]
                   + f_3 * pc_x[k] * ksk_791[k];

        t_981[k] = pa_x[k] * isl0_981[k]
                   - f_14 * pc_x[k] * isl1_981[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, t_985, pa_x, pc_x, pc_z, isl0_983, isl0_984, \
                         isl0_985, isl1_983, isl1_984, isl1_985, \
                         ksk_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_3 * pc_z[k] * ksk_784[k];

        t_983[k] = pa_x[k] * isl0_983[k]
                   - f_14 * pc_x[k] * isl1_983[k];

        t_984[k] = pa_x[k] * isl0_984[k]
                   - f_14 * pc_x[k] * isl1_984[k];

        t_985[k] = pa_x[k] * isl0_985[k]
                   - f_14 * pc_x[k] * isl1_985[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, t_989, pa_x, pc_x, pc_y, isl0_986, isl0_987, \
                         isl0_989, isk_575, isl1_986, isl1_987, isl1_989, \
                         ksk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = pa_x[k] * isl0_986[k]
                   - f_14 * pc_x[k] * isl1_986[k];

        t_987[k] = pa_x[k] * isl0_987[k]
                   - f_14 * pc_x[k] * isl1_987[k];

        t_988[k] = f_20 * isk_575[k]
                   + f_3 * pc_y[k] * ksk_791[k];

        t_989[k] = pa_x[k] * isl0_989[k]
                   - f_14 * pc_x[k] * isl1_989[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, pa_z, pc_y, pc_z, isl0_675, isl0_678, \
                         isk_540, isk_576, isl1_675, isl1_678, \
                         ksk_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = pa_z[k] * isl0_675[k]
                   - f_14 * pc_z[k] * isl1_675[k];

        t_991[k] = f_19 * isk_576[k]
                   + f_3 * pc_y[k] * ksk_792[k];

        t_992[k] = f_15 * isk_540[k]
                   + f_3 * pc_z[k] * ksk_792[k];

        t_993[k] = pa_z[k] * isl0_678[k]
                   - f_14 * pc_z[k] * isl1_678[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, pa_x, pa_z, pc_x, pc_y, pc_z, isl0_681, \
                         isl0_995, isk_578, isk_797, isl1_681, isl1_995, \
                         ksk_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_19 * isk_578[k]
                   + f_3 * pc_y[k] * ksk_794[k];

        t_995[k] = pa_x[k] * isl0_995[k]
                   + f_20 * isk_797[k]
                   - f_14 * pc_x[k] * isl1_995[k];

        t_996[k] = pa_z[k] * isl0_681[k]
                   - f_14 * pc_z[k] * isl1_681[k];
    }

#pragma omp simd aligned(t_997, t_998, t_999, pa_x, pc_x, pc_y, pc_z, isl0_999, isk_543, \
                         isk_581, isk_801, isl1_999, ksk_795, ksk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = f_15 * isk_543[k]
                   + f_3 * pc_z[k] * ksk_795[k];

        t_998[k] = f_19 * isk_581[k]
                   + f_3 * pc_y[k] * ksk_797[k];

        t_999[k] = pa_x[k] * isl0_999[k]
                   + f_19 * isk_801[k]
                   - f_14 * pc_x[k] * isl1_999[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, pa_x, pa_z, pc_x, pc_z, isl0_685, isl0_1002, \
                         isk_546, isk_804, isl1_685, isl1_1002, \
                         ksk_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = pa_z[k] * isl0_685[k]
                    - f_14 * pc_z[k] * isl1_685[k];

        t_1001[k] = f_15 * isk_546[k]
                    + f_3 * pc_z[k] * ksk_798[k];

        t_1002[k] = pa_x[k] * isl0_1002[k]
                    + f_18 * isk_804[k]
                    - f_14 * pc_x[k] * isl1_1002[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, pa_x, pa_z, pc_x, pc_y, pc_z, isl0_690, \
                         isl0_1004, isk_585, isk_806, isl1_690, isl1_1004, \
                         ksk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = f_19 * isk_585[k]
                    + f_3 * pc_y[k] * ksk_801[k];

        t_1004[k] = pa_x[k] * isl0_1004[k]
                    + f_18 * isk_806[k]
                    - f_14 * pc_x[k] * isl1_1004[k];

        t_1005[k] = pa_z[k] * isl0_690[k]
                    - f_14 * pc_z[k] * isl1_690[k];
    }

#pragma omp simd aligned(t_1006, t_1007, t_1008, pa_x, pc_x, pc_z, isl0_1007, isl0_1008, \
                         isk_550, isk_809, isk_810, isl1_1007, isl1_1008, \
                         ksk_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1006[k] = f_15 * isk_550[k]
                    + f_3 * pc_z[k] * ksk_802[k];

        t_1007[k] = pa_x[k] * isl0_1007[k]
                    + f_17 * isk_809[k]
                    - f_14 * pc_x[k] * isl1_1007[k];

        t_1008[k] = pa_x[k] * isl0_1008[k]
                    + f_17 * isk_810[k]
                    - f_14 * pc_x[k] * isl1_1008[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, pa_x, pa_z, pc_x, pc_y, pc_z, isl0_696, \
                         isl0_1010, isk_590, isk_812, isl1_696, isl1_1010, \
                         ksk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_19 * isk_590[k]
                    + f_3 * pc_y[k] * ksk_806[k];

        t_1010[k] = pa_x[k] * isl0_1010[k]
                    + f_17 * isk_812[k]
                    - f_14 * pc_x[k] * isl1_1010[k];

        t_1011[k] = pa_z[k] * isl0_696[k]
                    - f_14 * pc_z[k] * isl1_696[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pa_x, pc_x, pc_z, isl0_1013, isl0_1014, \
                         isk_555, isk_815, isk_816, isl1_1013, isl1_1014, \
                         ksk_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_15 * isk_555[k]
                    + f_3 * pc_z[k] * ksk_807[k];

        t_1013[k] = pa_x[k] * isl0_1013[k]
                    + f_16 * isk_815[k]
                    - f_14 * pc_x[k] * isl1_1013[k];

        t_1014[k] = pa_x[k] * isl0_1014[k]
                    + f_16 * isk_816[k]
                    - f_14 * pc_x[k] * isl1_1014[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, pa_x, pc_x, pc_y, isl0_1015, isl0_1017, \
                         isk_596, isk_817, isk_819, isl1_1015, isl1_1017, \
                         ksk_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = pa_x[k] * isl0_1015[k]
                    + f_16 * isk_817[k]
                    - f_14 * pc_x[k] * isl1_1015[k];

        t_1016[k] = f_19 * isk_596[k]
                    + f_3 * pc_y[k] * ksk_812[k];

        t_1017[k] = pa_x[k] * isl0_1017[k]
                    + f_16 * isk_819[k]
                    - f_14 * pc_x[k] * isl1_1017[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, t_1021, t_1022, pc_x, isk_820, isk_821, \
                         isk_822, isk_823, isk_824, ksk_820, ksk_821, ksk_822, ksk_823, \
                         ksk_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_15 * isk_820[k]
                    + f_3 * pc_x[k] * ksk_820[k];

        t_1019[k] = f_15 * isk_821[k]
                    + f_3 * pc_x[k] * ksk_821[k];

        t_1020[k] = f_15 * isk_822[k]
                    + f_3 * pc_x[k] * ksk_822[k];

        t_1021[k] = f_15 * isk_823[k]
                    + f_3 * pc_x[k] * ksk_823[k];

        t_1022[k] = f_15 * isk_824[k]
                    + f_3 * pc_x[k] * ksk_824[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, t_1026, pa_x, pc_x, isl0_1026, isk_825, \
                         isk_826, isk_827, isl1_1026, ksk_825, ksk_826, \
                         ksk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = f_15 * isk_825[k]
                    + f_3 * pc_x[k] * ksk_825[k];

        t_1024[k] = f_15 * isk_826[k]
                    + f_3 * pc_x[k] * ksk_826[k];

        t_1025[k] = f_15 * isk_827[k]
                    + f_3 * pc_x[k] * ksk_827[k];

        t_1026[k] = pa_x[k] * isl0_1026[k]
                    - f_14 * pc_x[k] * isl1_1026[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, pa_x, pc_x, pc_z, isl0_1028, \
                         isl0_1029, isl0_1030, isk_568, isl1_1028, isl1_1029, isl1_1030, \
                         ksk_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_15 * isk_568[k]
                    + f_3 * pc_z[k] * ksk_820[k];

        t_1028[k] = pa_x[k] * isl0_1028[k]
                    - f_14 * pc_x[k] * isl1_1028[k];

        t_1029[k] = pa_x[k] * isl0_1029[k]
                    - f_14 * pc_x[k] * isl1_1029[k];

        t_1030[k] = pa_x[k] * isl0_1030[k]
                    - f_14 * pc_x[k] * isl1_1030[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, t_1034, pa_x, pc_x, pc_y, isl0_1031, \
                         isl0_1032, isl0_1034, isk_611, isl1_1031, isl1_1032, isl1_1034, \
                         ksk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = pa_x[k] * isl0_1031[k]
                    - f_14 * pc_x[k] * isl1_1031[k];

        t_1032[k] = pa_x[k] * isl0_1032[k]
                    - f_14 * pc_x[k] * isl1_1032[k];

        t_1033[k] = f_19 * isk_611[k]
                    + f_3 * pc_y[k] * ksk_827[k];

        t_1034[k] = pa_x[k] * isl0_1034[k]
                    - f_14 * pc_x[k] * isl1_1034[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, pa_x, pc_x, pc_y, pc_z, isl0_1035, isk_576, \
                         isk_612, isk_828, isl1_1035, ksk_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = pa_x[k] * isl0_1035[k]
                    + f_23 * isk_828[k]
                    - f_14 * pc_x[k] * isl1_1035[k];

        t_1036[k] = f_18 * isk_612[k]
                    + f_3 * pc_y[k] * ksk_828[k];

        t_1037[k] = f_16 * isk_576[k]
                    + f_3 * pc_z[k] * ksk_828[k];
    }

#pragma omp simd aligned(t_1038, t_1039, t_1040, pa_x, pc_x, pc_y, isl0_1038, isl0_1040, \
                         isk_614, isk_831, isk_833, isl1_1038, isl1_1040, \
                         ksk_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = pa_x[k] * isl0_1038[k]
                    + f_20 * isk_831[k]
                    - f_14 * pc_x[k] * isl1_1038[k];

        t_1039[k] = f_18 * isk_614[k]
                    + f_3 * pc_y[k] * ksk_830[k];

        t_1040[k] = pa_x[k] * isl0_1040[k]
                    + f_20 * isk_833[k]
                    - f_14 * pc_x[k] * isl1_1040[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isl0,
                                                          const size_t isk, const size_t isl1,
                                                          const size_t ksk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_23 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isl0_1041 = buffer.data(isl0 + 1041);
    const auto *isl0_1044 = buffer.data(isl0 + 1044);
    const auto *isl0_1045 = buffer.data(isl0 + 1045);
    const auto *isl0_1047 = buffer.data(isl0 + 1047);
    const auto *isl0_1049 = buffer.data(isl0 + 1049);
    const auto *isl0_1050 = buffer.data(isl0 + 1050);
    const auto *isl0_1052 = buffer.data(isl0 + 1052);
    const auto *isl0_1053 = buffer.data(isl0 + 1053);
    const auto *isl0_1055 = buffer.data(isl0 + 1055);
    const auto *isl0_1056 = buffer.data(isl0 + 1056);
    const auto *isl0_1058 = buffer.data(isl0 + 1058);
    const auto *isl0_1059 = buffer.data(isl0 + 1059);
    const auto *isl0_1060 = buffer.data(isl0 + 1060);
    const auto *isl0_1062 = buffer.data(isl0 + 1062);
    const auto *isl0_1071 = buffer.data(isl0 + 1071);
    const auto *isl0_1073 = buffer.data(isl0 + 1073);
    const auto *isl0_1074 = buffer.data(isl0 + 1074);
    const auto *isl0_1075 = buffer.data(isl0 + 1075);
    const auto *isl0_1076 = buffer.data(isl0 + 1076);
    const auto *isl0_1077 = buffer.data(isl0 + 1077);
    const auto *isl0_1079 = buffer.data(isl0 + 1079);
    const auto *isl0_1080 = buffer.data(isl0 + 1080);
    const auto *isl0_1083 = buffer.data(isl0 + 1083);
    const auto *isl0_1085 = buffer.data(isl0 + 1085);
    const auto *isl0_1086 = buffer.data(isl0 + 1086);
    const auto *isl0_1089 = buffer.data(isl0 + 1089);
    const auto *isl0_1090 = buffer.data(isl0 + 1090);
    const auto *isl0_1092 = buffer.data(isl0 + 1092);
    const auto *isl0_1094 = buffer.data(isl0 + 1094);
    const auto *isl0_1095 = buffer.data(isl0 + 1095);
    const auto *isl0_1097 = buffer.data(isl0 + 1097);
    const auto *isl0_1098 = buffer.data(isl0 + 1098);
    const auto *isl0_1100 = buffer.data(isl0 + 1100);
    const auto *isl0_1101 = buffer.data(isl0 + 1101);
    const auto *isl0_1103 = buffer.data(isl0 + 1103);
    const auto *isl0_1104 = buffer.data(isl0 + 1104);
    const auto *isl0_1105 = buffer.data(isl0 + 1105);
    const auto *isl0_1107 = buffer.data(isl0 + 1107);
    const auto *isl0_1116 = buffer.data(isl0 + 1116);
    const auto *isl0_1118 = buffer.data(isl0 + 1118);
    const auto *isl0_1119 = buffer.data(isl0 + 1119);
    const auto *isl0_1120 = buffer.data(isl0 + 1120);
    const auto *isl0_1121 = buffer.data(isl0 + 1121);
    const auto *isl0_1122 = buffer.data(isl0 + 1122);
    const auto *isl0_1124 = buffer.data(isl0 + 1124);
    const auto *isl0_1125 = buffer.data(isl0 + 1125);
    const auto *isl0_1128 = buffer.data(isl0 + 1128);
    const auto *isl0_1130 = buffer.data(isl0 + 1130);
    const auto *isl0_1131 = buffer.data(isl0 + 1131);
    const auto *isl0_1134 = buffer.data(isl0 + 1134);
    const auto *isl0_1135 = buffer.data(isl0 + 1135);
    const auto *isl0_1137 = buffer.data(isl0 + 1137);
    const auto *isl0_1139 = buffer.data(isl0 + 1139);
    const auto *isl0_1140 = buffer.data(isl0 + 1140);
    const auto *isl0_1142 = buffer.data(isl0 + 1142);
    const auto *isl0_1143 = buffer.data(isl0 + 1143);
    const auto *isl0_1145 = buffer.data(isl0 + 1145);
    const auto *isl0_1146 = buffer.data(isl0 + 1146);
    const auto *isl0_1148 = buffer.data(isl0 + 1148);
    const auto *isl0_1149 = buffer.data(isl0 + 1149);
    const auto *isl0_1150 = buffer.data(isl0 + 1150);
    const auto *isl0_1152 = buffer.data(isl0 + 1152);

    const auto *isk_579 = buffer.data(isk + 579);
    const auto *isk_582 = buffer.data(isk + 582);
    const auto *isk_586 = buffer.data(isk + 586);
    const auto *isk_591 = buffer.data(isk + 591);
    const auto *isk_604 = buffer.data(isk + 604);
    const auto *isk_612 = buffer.data(isk + 612);
    const auto *isk_615 = buffer.data(isk + 615);
    const auto *isk_617 = buffer.data(isk + 617);
    const auto *isk_618 = buffer.data(isk + 618);
    const auto *isk_621 = buffer.data(isk + 621);
    const auto *isk_622 = buffer.data(isk + 622);
    const auto *isk_626 = buffer.data(isk + 626);
    const auto *isk_627 = buffer.data(isk + 627);
    const auto *isk_632 = buffer.data(isk + 632);
    const auto *isk_640 = buffer.data(isk + 640);
    const auto *isk_647 = buffer.data(isk + 647);
    const auto *isk_648 = buffer.data(isk + 648);
    const auto *isk_650 = buffer.data(isk + 650);
    const auto *isk_651 = buffer.data(isk + 651);
    const auto *isk_653 = buffer.data(isk + 653);
    const auto *isk_654 = buffer.data(isk + 654);
    const auto *isk_657 = buffer.data(isk + 657);
    const auto *isk_658 = buffer.data(isk + 658);
    const auto *isk_662 = buffer.data(isk + 662);
    const auto *isk_663 = buffer.data(isk + 663);
    const auto *isk_668 = buffer.data(isk + 668);
    const auto *isk_683 = buffer.data(isk + 683);
    const auto *isk_684 = buffer.data(isk + 684);
    const auto *isk_686 = buffer.data(isk + 686);
    const auto *isk_689 = buffer.data(isk + 689);
    const auto *isk_693 = buffer.data(isk + 693);
    const auto *isk_698 = buffer.data(isk + 698);
    const auto *isk_704 = buffer.data(isk + 704);
    const auto *isk_834 = buffer.data(isk + 834);
    const auto *isk_837 = buffer.data(isk + 837);
    const auto *isk_838 = buffer.data(isk + 838);
    const auto *isk_840 = buffer.data(isk + 840);
    const auto *isk_842 = buffer.data(isk + 842);
    const auto *isk_843 = buffer.data(isk + 843);
    const auto *isk_845 = buffer.data(isk + 845);
    const auto *isk_846 = buffer.data(isk + 846);
    const auto *isk_848 = buffer.data(isk + 848);
    const auto *isk_849 = buffer.data(isk + 849);
    const auto *isk_851 = buffer.data(isk + 851);
    const auto *isk_852 = buffer.data(isk + 852);
    const auto *isk_853 = buffer.data(isk + 853);
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
    const auto *isk_867 = buffer.data(isk + 867);
    const auto *isk_869 = buffer.data(isk + 869);
    const auto *isk_870 = buffer.data(isk + 870);
    const auto *isk_873 = buffer.data(isk + 873);
    const auto *isk_874 = buffer.data(isk + 874);
    const auto *isk_876 = buffer.data(isk + 876);
    const auto *isk_878 = buffer.data(isk + 878);
    const auto *isk_879 = buffer.data(isk + 879);
    const auto *isk_881 = buffer.data(isk + 881);
    const auto *isk_882 = buffer.data(isk + 882);
    const auto *isk_884 = buffer.data(isk + 884);
    const auto *isk_885 = buffer.data(isk + 885);
    const auto *isk_887 = buffer.data(isk + 887);
    const auto *isk_888 = buffer.data(isk + 888);
    const auto *isk_889 = buffer.data(isk + 889);
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
    const auto *isk_903 = buffer.data(isk + 903);
    const auto *isk_905 = buffer.data(isk + 905);
    const auto *isk_906 = buffer.data(isk + 906);
    const auto *isk_909 = buffer.data(isk + 909);
    const auto *isk_910 = buffer.data(isk + 910);
    const auto *isk_912 = buffer.data(isk + 912);
    const auto *isk_914 = buffer.data(isk + 914);
    const auto *isk_915 = buffer.data(isk + 915);
    const auto *isk_917 = buffer.data(isk + 917);
    const auto *isk_918 = buffer.data(isk + 918);
    const auto *isk_920 = buffer.data(isk + 920);
    const auto *isk_921 = buffer.data(isk + 921);
    const auto *isk_923 = buffer.data(isk + 923);
    const auto *isk_924 = buffer.data(isk + 924);
    const auto *isk_925 = buffer.data(isk + 925);
    const auto *isk_927 = buffer.data(isk + 927);
    const auto *isk_928 = buffer.data(isk + 928);
    const auto *isk_929 = buffer.data(isk + 929);
    const auto *isk_930 = buffer.data(isk + 930);

    const auto *isl1_1041 = buffer.data(isl1 + 1041);
    const auto *isl1_1044 = buffer.data(isl1 + 1044);
    const auto *isl1_1045 = buffer.data(isl1 + 1045);
    const auto *isl1_1047 = buffer.data(isl1 + 1047);
    const auto *isl1_1049 = buffer.data(isl1 + 1049);
    const auto *isl1_1050 = buffer.data(isl1 + 1050);
    const auto *isl1_1052 = buffer.data(isl1 + 1052);
    const auto *isl1_1053 = buffer.data(isl1 + 1053);
    const auto *isl1_1055 = buffer.data(isl1 + 1055);
    const auto *isl1_1056 = buffer.data(isl1 + 1056);
    const auto *isl1_1058 = buffer.data(isl1 + 1058);
    const auto *isl1_1059 = buffer.data(isl1 + 1059);
    const auto *isl1_1060 = buffer.data(isl1 + 1060);
    const auto *isl1_1062 = buffer.data(isl1 + 1062);
    const auto *isl1_1071 = buffer.data(isl1 + 1071);
    const auto *isl1_1073 = buffer.data(isl1 + 1073);
    const auto *isl1_1074 = buffer.data(isl1 + 1074);
    const auto *isl1_1075 = buffer.data(isl1 + 1075);
    const auto *isl1_1076 = buffer.data(isl1 + 1076);
    const auto *isl1_1077 = buffer.data(isl1 + 1077);
    const auto *isl1_1079 = buffer.data(isl1 + 1079);
    const auto *isl1_1080 = buffer.data(isl1 + 1080);
    const auto *isl1_1083 = buffer.data(isl1 + 1083);
    const auto *isl1_1085 = buffer.data(isl1 + 1085);
    const auto *isl1_1086 = buffer.data(isl1 + 1086);
    const auto *isl1_1089 = buffer.data(isl1 + 1089);
    const auto *isl1_1090 = buffer.data(isl1 + 1090);
    const auto *isl1_1092 = buffer.data(isl1 + 1092);
    const auto *isl1_1094 = buffer.data(isl1 + 1094);
    const auto *isl1_1095 = buffer.data(isl1 + 1095);
    const auto *isl1_1097 = buffer.data(isl1 + 1097);
    const auto *isl1_1098 = buffer.data(isl1 + 1098);
    const auto *isl1_1100 = buffer.data(isl1 + 1100);
    const auto *isl1_1101 = buffer.data(isl1 + 1101);
    const auto *isl1_1103 = buffer.data(isl1 + 1103);
    const auto *isl1_1104 = buffer.data(isl1 + 1104);
    const auto *isl1_1105 = buffer.data(isl1 + 1105);
    const auto *isl1_1107 = buffer.data(isl1 + 1107);
    const auto *isl1_1116 = buffer.data(isl1 + 1116);
    const auto *isl1_1118 = buffer.data(isl1 + 1118);
    const auto *isl1_1119 = buffer.data(isl1 + 1119);
    const auto *isl1_1120 = buffer.data(isl1 + 1120);
    const auto *isl1_1121 = buffer.data(isl1 + 1121);
    const auto *isl1_1122 = buffer.data(isl1 + 1122);
    const auto *isl1_1124 = buffer.data(isl1 + 1124);
    const auto *isl1_1125 = buffer.data(isl1 + 1125);
    const auto *isl1_1128 = buffer.data(isl1 + 1128);
    const auto *isl1_1130 = buffer.data(isl1 + 1130);
    const auto *isl1_1131 = buffer.data(isl1 + 1131);
    const auto *isl1_1134 = buffer.data(isl1 + 1134);
    const auto *isl1_1135 = buffer.data(isl1 + 1135);
    const auto *isl1_1137 = buffer.data(isl1 + 1137);
    const auto *isl1_1139 = buffer.data(isl1 + 1139);
    const auto *isl1_1140 = buffer.data(isl1 + 1140);
    const auto *isl1_1142 = buffer.data(isl1 + 1142);
    const auto *isl1_1143 = buffer.data(isl1 + 1143);
    const auto *isl1_1145 = buffer.data(isl1 + 1145);
    const auto *isl1_1146 = buffer.data(isl1 + 1146);
    const auto *isl1_1148 = buffer.data(isl1 + 1148);
    const auto *isl1_1149 = buffer.data(isl1 + 1149);
    const auto *isl1_1150 = buffer.data(isl1 + 1150);
    const auto *isl1_1152 = buffer.data(isl1 + 1152);

    const auto *ksk_831 = buffer.data(ksk + 831);
    const auto *ksk_833 = buffer.data(ksk + 833);
    const auto *ksk_834 = buffer.data(ksk + 834);
    const auto *ksk_837 = buffer.data(ksk + 837);
    const auto *ksk_838 = buffer.data(ksk + 838);
    const auto *ksk_842 = buffer.data(ksk + 842);
    const auto *ksk_843 = buffer.data(ksk + 843);
    const auto *ksk_848 = buffer.data(ksk + 848);
    const auto *ksk_856 = buffer.data(ksk + 856);
    const auto *ksk_857 = buffer.data(ksk + 857);
    const auto *ksk_858 = buffer.data(ksk + 858);
    const auto *ksk_859 = buffer.data(ksk + 859);
    const auto *ksk_860 = buffer.data(ksk + 860);
    const auto *ksk_861 = buffer.data(ksk + 861);
    const auto *ksk_862 = buffer.data(ksk + 862);
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
    const auto *ksk_893 = buffer.data(ksk + 893);
    const auto *ksk_894 = buffer.data(ksk + 894);
    const auto *ksk_895 = buffer.data(ksk + 895);
    const auto *ksk_896 = buffer.data(ksk + 896);
    const auto *ksk_897 = buffer.data(ksk + 897);
    const auto *ksk_898 = buffer.data(ksk + 898);
    const auto *ksk_899 = buffer.data(ksk + 899);
    const auto *ksk_900 = buffer.data(ksk + 900);
    const auto *ksk_902 = buffer.data(ksk + 902);
    const auto *ksk_903 = buffer.data(ksk + 903);
    const auto *ksk_905 = buffer.data(ksk + 905);
    const auto *ksk_906 = buffer.data(ksk + 906);
    const auto *ksk_909 = buffer.data(ksk + 909);
    const auto *ksk_910 = buffer.data(ksk + 910);
    const auto *ksk_914 = buffer.data(ksk + 914);
    const auto *ksk_915 = buffer.data(ksk + 915);
    const auto *ksk_920 = buffer.data(ksk + 920);
    const auto *ksk_928 = buffer.data(ksk + 928);
    const auto *ksk_929 = buffer.data(ksk + 929);
    const auto *ksk_930 = buffer.data(ksk + 930);

#pragma omp simd aligned(t_1041, t_1042, t_1043, pa_x, pc_x, pc_y, pc_z, isl0_1041, isk_579, \
                         isk_617, isk_834, isl1_1041, ksk_831, \
                         ksk_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = pa_x[k] * isl0_1041[k]
                    + f_19 * isk_834[k]
                    - f_14 * pc_x[k] * isl1_1041[k];

        t_1042[k] = f_16 * isk_579[k]
                    + f_3 * pc_z[k] * ksk_831[k];

        t_1043[k] = f_18 * isk_617[k]
                    + f_3 * pc_y[k] * ksk_833[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, pa_x, pc_x, pc_z, isl0_1044, isl0_1045, \
                         isk_582, isk_837, isk_838, isl1_1044, isl1_1045, \
                         ksk_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = pa_x[k] * isl0_1044[k]
                    + f_19 * isk_837[k]
                    - f_14 * pc_x[k] * isl1_1044[k];

        t_1045[k] = pa_x[k] * isl0_1045[k]
                    + f_18 * isk_838[k]
                    - f_14 * pc_x[k] * isl1_1045[k];

        t_1046[k] = f_16 * isk_582[k]
                    + f_3 * pc_z[k] * ksk_834[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pa_x, pc_x, pc_y, isl0_1047, isl0_1049, \
                         isk_621, isk_840, isk_842, isl1_1047, isl1_1049, \
                         ksk_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = pa_x[k] * isl0_1047[k]
                    + f_18 * isk_840[k]
                    - f_14 * pc_x[k] * isl1_1047[k];

        t_1048[k] = f_18 * isk_621[k]
                    + f_3 * pc_y[k] * ksk_837[k];

        t_1049[k] = pa_x[k] * isl0_1049[k]
                    + f_18 * isk_842[k]
                    - f_14 * pc_x[k] * isl1_1049[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, pa_x, pc_x, pc_z, isl0_1050, isl0_1052, \
                         isk_586, isk_843, isk_845, isl1_1050, isl1_1052, \
                         ksk_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = pa_x[k] * isl0_1050[k]
                    + f_17 * isk_843[k]
                    - f_14 * pc_x[k] * isl1_1050[k];

        t_1051[k] = f_16 * isk_586[k]
                    + f_3 * pc_z[k] * ksk_838[k];

        t_1052[k] = pa_x[k] * isl0_1052[k]
                    + f_17 * isk_845[k]
                    - f_14 * pc_x[k] * isl1_1052[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, pa_x, pc_x, pc_y, isl0_1053, isl0_1055, \
                         isk_626, isk_846, isk_848, isl1_1053, isl1_1055, \
                         ksk_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = pa_x[k] * isl0_1053[k]
                    + f_17 * isk_846[k]
                    - f_14 * pc_x[k] * isl1_1053[k];

        t_1054[k] = f_18 * isk_626[k]
                    + f_3 * pc_y[k] * ksk_842[k];

        t_1055[k] = pa_x[k] * isl0_1055[k]
                    + f_17 * isk_848[k]
                    - f_14 * pc_x[k] * isl1_1055[k];
    }

#pragma omp simd aligned(t_1056, t_1057, t_1058, pa_x, pc_x, pc_z, isl0_1056, isl0_1058, \
                         isk_591, isk_849, isk_851, isl1_1056, isl1_1058, \
                         ksk_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1056[k] = pa_x[k] * isl0_1056[k]
                    + f_16 * isk_849[k]
                    - f_14 * pc_x[k] * isl1_1056[k];

        t_1057[k] = f_16 * isk_591[k]
                    + f_3 * pc_z[k] * ksk_843[k];

        t_1058[k] = pa_x[k] * isl0_1058[k]
                    + f_16 * isk_851[k]
                    - f_14 * pc_x[k] * isl1_1058[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, pa_x, pc_x, pc_y, isl0_1059, isl0_1060, \
                         isk_632, isk_852, isk_853, isl1_1059, isl1_1060, \
                         ksk_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = pa_x[k] * isl0_1059[k]
                    + f_16 * isk_852[k]
                    - f_14 * pc_x[k] * isl1_1059[k];

        t_1060[k] = pa_x[k] * isl0_1060[k]
                    + f_16 * isk_853[k]
                    - f_14 * pc_x[k] * isl1_1060[k];

        t_1061[k] = f_18 * isk_632[k]
                    + f_3 * pc_y[k] * ksk_848[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, t_1065, pa_x, pc_x, isl0_1062, isk_855, \
                         isk_856, isk_857, isk_858, isl1_1062, ksk_856, ksk_857, \
                         ksk_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = pa_x[k] * isl0_1062[k]
                    + f_16 * isk_855[k]
                    - f_14 * pc_x[k] * isl1_1062[k];

        t_1063[k] = f_15 * isk_856[k]
                    + f_3 * pc_x[k] * ksk_856[k];

        t_1064[k] = f_15 * isk_857[k]
                    + f_3 * pc_x[k] * ksk_857[k];

        t_1065[k] = f_15 * isk_858[k]
                    + f_3 * pc_x[k] * ksk_858[k];
    }

#pragma omp simd aligned(t_1066, t_1067, t_1068, t_1069, t_1070, pc_x, isk_859, isk_860, \
                         isk_861, isk_862, isk_863, ksk_859, ksk_860, ksk_861, ksk_862, \
                         ksk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1066[k] = f_15 * isk_859[k]
                    + f_3 * pc_x[k] * ksk_859[k];

        t_1067[k] = f_15 * isk_860[k]
                    + f_3 * pc_x[k] * ksk_860[k];

        t_1068[k] = f_15 * isk_861[k]
                    + f_3 * pc_x[k] * ksk_861[k];

        t_1069[k] = f_15 * isk_862[k]
                    + f_3 * pc_x[k] * ksk_862[k];

        t_1070[k] = f_15 * isk_863[k]
                    + f_3 * pc_x[k] * ksk_863[k];
    }

#pragma omp simd aligned(t_1071, t_1072, t_1073, t_1074, pa_x, pc_x, pc_z, isl0_1071, \
                         isl0_1073, isl0_1074, isk_604, isl1_1071, isl1_1073, isl1_1074, \
                         ksk_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1071[k] = pa_x[k] * isl0_1071[k]
                    - f_14 * pc_x[k] * isl1_1071[k];

        t_1072[k] = f_16 * isk_604[k]
                    + f_3 * pc_z[k] * ksk_856[k];

        t_1073[k] = pa_x[k] * isl0_1073[k]
                    - f_14 * pc_x[k] * isl1_1073[k];

        t_1074[k] = pa_x[k] * isl0_1074[k]
                    - f_14 * pc_x[k] * isl1_1074[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, t_1078, pa_x, pc_x, pc_y, isl0_1075, \
                         isl0_1076, isl0_1077, isk_647, isl1_1075, isl1_1076, isl1_1077, \
                         ksk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = pa_x[k] * isl0_1075[k]
                    - f_14 * pc_x[k] * isl1_1075[k];

        t_1076[k] = pa_x[k] * isl0_1076[k]
                    - f_14 * pc_x[k] * isl1_1076[k];

        t_1077[k] = pa_x[k] * isl0_1077[k]
                    - f_14 * pc_x[k] * isl1_1077[k];

        t_1078[k] = f_18 * isk_647[k]
                    + f_3 * pc_y[k] * ksk_863[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, pa_x, pc_x, pc_y, pc_z, isl0_1079, \
                         isl0_1080, isk_612, isk_648, isk_864, isl1_1079, isl1_1080, \
                         ksk_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = pa_x[k] * isl0_1079[k]
                    - f_14 * pc_x[k] * isl1_1079[k];

        t_1080[k] = pa_x[k] * isl0_1080[k]
                    + f_23 * isk_864[k]
                    - f_14 * pc_x[k] * isl1_1080[k];

        t_1081[k] = f_17 * isk_648[k]
                    + f_3 * pc_y[k] * ksk_864[k];

        t_1082[k] = f_17 * isk_612[k]
                    + f_3 * pc_z[k] * ksk_864[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, pa_x, pc_x, pc_y, isl0_1083, isl0_1085, \
                         isk_650, isk_867, isk_869, isl1_1083, isl1_1085, \
                         ksk_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = pa_x[k] * isl0_1083[k]
                    + f_20 * isk_867[k]
                    - f_14 * pc_x[k] * isl1_1083[k];

        t_1084[k] = f_17 * isk_650[k]
                    + f_3 * pc_y[k] * ksk_866[k];

        t_1085[k] = pa_x[k] * isl0_1085[k]
                    + f_20 * isk_869[k]
                    - f_14 * pc_x[k] * isl1_1085[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pa_x, pc_x, pc_y, pc_z, isl0_1086, isk_615, \
                         isk_653, isk_870, isl1_1086, ksk_867, \
                         ksk_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = pa_x[k] * isl0_1086[k]
                    + f_19 * isk_870[k]
                    - f_14 * pc_x[k] * isl1_1086[k];

        t_1087[k] = f_17 * isk_615[k]
                    + f_3 * pc_z[k] * ksk_867[k];

        t_1088[k] = f_17 * isk_653[k]
                    + f_3 * pc_y[k] * ksk_869[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pa_x, pc_x, pc_z, isl0_1089, isl0_1090, \
                         isk_618, isk_873, isk_874, isl1_1089, isl1_1090, \
                         ksk_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = pa_x[k] * isl0_1089[k]
                    + f_19 * isk_873[k]
                    - f_14 * pc_x[k] * isl1_1089[k];

        t_1090[k] = pa_x[k] * isl0_1090[k]
                    + f_18 * isk_874[k]
                    - f_14 * pc_x[k] * isl1_1090[k];

        t_1091[k] = f_17 * isk_618[k]
                    + f_3 * pc_z[k] * ksk_870[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pa_x, pc_x, pc_y, isl0_1092, isl0_1094, \
                         isk_657, isk_876, isk_878, isl1_1092, isl1_1094, \
                         ksk_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = pa_x[k] * isl0_1092[k]
                    + f_18 * isk_876[k]
                    - f_14 * pc_x[k] * isl1_1092[k];

        t_1093[k] = f_17 * isk_657[k]
                    + f_3 * pc_y[k] * ksk_873[k];

        t_1094[k] = pa_x[k] * isl0_1094[k]
                    + f_18 * isk_878[k]
                    - f_14 * pc_x[k] * isl1_1094[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pa_x, pc_x, pc_z, isl0_1095, isl0_1097, \
                         isk_622, isk_879, isk_881, isl1_1095, isl1_1097, \
                         ksk_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = pa_x[k] * isl0_1095[k]
                    + f_17 * isk_879[k]
                    - f_14 * pc_x[k] * isl1_1095[k];

        t_1096[k] = f_17 * isk_622[k]
                    + f_3 * pc_z[k] * ksk_874[k];

        t_1097[k] = pa_x[k] * isl0_1097[k]
                    + f_17 * isk_881[k]
                    - f_14 * pc_x[k] * isl1_1097[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pa_x, pc_x, pc_y, isl0_1098, isl0_1100, \
                         isk_662, isk_882, isk_884, isl1_1098, isl1_1100, \
                         ksk_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = pa_x[k] * isl0_1098[k]
                    + f_17 * isk_882[k]
                    - f_14 * pc_x[k] * isl1_1098[k];

        t_1099[k] = f_17 * isk_662[k]
                    + f_3 * pc_y[k] * ksk_878[k];

        t_1100[k] = pa_x[k] * isl0_1100[k]
                    + f_17 * isk_884[k]
                    - f_14 * pc_x[k] * isl1_1100[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, pa_x, pc_x, pc_z, isl0_1101, isl0_1103, \
                         isk_627, isk_885, isk_887, isl1_1101, isl1_1103, \
                         ksk_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = pa_x[k] * isl0_1101[k]
                    + f_16 * isk_885[k]
                    - f_14 * pc_x[k] * isl1_1101[k];

        t_1102[k] = f_17 * isk_627[k]
                    + f_3 * pc_z[k] * ksk_879[k];

        t_1103[k] = pa_x[k] * isl0_1103[k]
                    + f_16 * isk_887[k]
                    - f_14 * pc_x[k] * isl1_1103[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, pa_x, pc_x, pc_y, isl0_1104, isl0_1105, \
                         isk_668, isk_888, isk_889, isl1_1104, isl1_1105, \
                         ksk_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = pa_x[k] * isl0_1104[k]
                    + f_16 * isk_888[k]
                    - f_14 * pc_x[k] * isl1_1104[k];

        t_1105[k] = pa_x[k] * isl0_1105[k]
                    + f_16 * isk_889[k]
                    - f_14 * pc_x[k] * isl1_1105[k];

        t_1106[k] = f_17 * isk_668[k]
                    + f_3 * pc_y[k] * ksk_884[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, t_1110, pa_x, pc_x, isl0_1107, isk_891, \
                         isk_892, isk_893, isk_894, isl1_1107, ksk_892, ksk_893, \
                         ksk_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = pa_x[k] * isl0_1107[k]
                    + f_16 * isk_891[k]
                    - f_14 * pc_x[k] * isl1_1107[k];

        t_1108[k] = f_15 * isk_892[k]
                    + f_3 * pc_x[k] * ksk_892[k];

        t_1109[k] = f_15 * isk_893[k]
                    + f_3 * pc_x[k] * ksk_893[k];

        t_1110[k] = f_15 * isk_894[k]
                    + f_3 * pc_x[k] * ksk_894[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, t_1114, t_1115, pc_x, isk_895, isk_896, \
                         isk_897, isk_898, isk_899, ksk_895, ksk_896, ksk_897, ksk_898, \
                         ksk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_15 * isk_895[k]
                    + f_3 * pc_x[k] * ksk_895[k];

        t_1112[k] = f_15 * isk_896[k]
                    + f_3 * pc_x[k] * ksk_896[k];

        t_1113[k] = f_15 * isk_897[k]
                    + f_3 * pc_x[k] * ksk_897[k];

        t_1114[k] = f_15 * isk_898[k]
                    + f_3 * pc_x[k] * ksk_898[k];

        t_1115[k] = f_15 * isk_899[k]
                    + f_3 * pc_x[k] * ksk_899[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, t_1119, pa_x, pc_x, pc_z, isl0_1116, \
                         isl0_1118, isl0_1119, isk_640, isl1_1116, isl1_1118, isl1_1119, \
                         ksk_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = pa_x[k] * isl0_1116[k]
                    - f_14 * pc_x[k] * isl1_1116[k];

        t_1117[k] = f_17 * isk_640[k]
                    + f_3 * pc_z[k] * ksk_892[k];

        t_1118[k] = pa_x[k] * isl0_1118[k]
                    - f_14 * pc_x[k] * isl1_1118[k];

        t_1119[k] = pa_x[k] * isl0_1119[k]
                    - f_14 * pc_x[k] * isl1_1119[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, t_1123, pa_x, pc_x, pc_y, isl0_1120, \
                         isl0_1121, isl0_1122, isk_683, isl1_1120, isl1_1121, isl1_1122, \
                         ksk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = pa_x[k] * isl0_1120[k]
                    - f_14 * pc_x[k] * isl1_1120[k];

        t_1121[k] = pa_x[k] * isl0_1121[k]
                    - f_14 * pc_x[k] * isl1_1121[k];

        t_1122[k] = pa_x[k] * isl0_1122[k]
                    - f_14 * pc_x[k] * isl1_1122[k];

        t_1123[k] = f_17 * isk_683[k]
                    + f_3 * pc_y[k] * ksk_899[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, t_1127, pa_x, pc_x, pc_y, pc_z, isl0_1124, \
                         isl0_1125, isk_648, isk_684, isk_900, isl1_1124, isl1_1125, \
                         ksk_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = pa_x[k] * isl0_1124[k]
                    - f_14 * pc_x[k] * isl1_1124[k];

        t_1125[k] = pa_x[k] * isl0_1125[k]
                    + f_23 * isk_900[k]
                    - f_14 * pc_x[k] * isl1_1125[k];

        t_1126[k] = f_16 * isk_684[k]
                    + f_3 * pc_y[k] * ksk_900[k];

        t_1127[k] = f_18 * isk_648[k]
                    + f_3 * pc_z[k] * ksk_900[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pa_x, pc_x, pc_y, isl0_1128, isl0_1130, \
                         isk_686, isk_903, isk_905, isl1_1128, isl1_1130, \
                         ksk_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = pa_x[k] * isl0_1128[k]
                    + f_20 * isk_903[k]
                    - f_14 * pc_x[k] * isl1_1128[k];

        t_1129[k] = f_16 * isk_686[k]
                    + f_3 * pc_y[k] * ksk_902[k];

        t_1130[k] = pa_x[k] * isl0_1130[k]
                    + f_20 * isk_905[k]
                    - f_14 * pc_x[k] * isl1_1130[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pa_x, pc_x, pc_y, pc_z, isl0_1131, isk_651, \
                         isk_689, isk_906, isl1_1131, ksk_903, \
                         ksk_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = pa_x[k] * isl0_1131[k]
                    + f_19 * isk_906[k]
                    - f_14 * pc_x[k] * isl1_1131[k];

        t_1132[k] = f_18 * isk_651[k]
                    + f_3 * pc_z[k] * ksk_903[k];

        t_1133[k] = f_16 * isk_689[k]
                    + f_3 * pc_y[k] * ksk_905[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, pa_x, pc_x, pc_z, isl0_1134, isl0_1135, \
                         isk_654, isk_909, isk_910, isl1_1134, isl1_1135, \
                         ksk_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = pa_x[k] * isl0_1134[k]
                    + f_19 * isk_909[k]
                    - f_14 * pc_x[k] * isl1_1134[k];

        t_1135[k] = pa_x[k] * isl0_1135[k]
                    + f_18 * isk_910[k]
                    - f_14 * pc_x[k] * isl1_1135[k];

        t_1136[k] = f_18 * isk_654[k]
                    + f_3 * pc_z[k] * ksk_906[k];
    }

#pragma omp simd aligned(t_1137, t_1138, t_1139, pa_x, pc_x, pc_y, isl0_1137, isl0_1139, \
                         isk_693, isk_912, isk_914, isl1_1137, isl1_1139, \
                         ksk_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = pa_x[k] * isl0_1137[k]
                    + f_18 * isk_912[k]
                    - f_14 * pc_x[k] * isl1_1137[k];

        t_1138[k] = f_16 * isk_693[k]
                    + f_3 * pc_y[k] * ksk_909[k];

        t_1139[k] = pa_x[k] * isl0_1139[k]
                    + f_18 * isk_914[k]
                    - f_14 * pc_x[k] * isl1_1139[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, pa_x, pc_x, pc_z, isl0_1140, isl0_1142, \
                         isk_658, isk_915, isk_917, isl1_1140, isl1_1142, \
                         ksk_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = pa_x[k] * isl0_1140[k]
                    + f_17 * isk_915[k]
                    - f_14 * pc_x[k] * isl1_1140[k];

        t_1141[k] = f_18 * isk_658[k]
                    + f_3 * pc_z[k] * ksk_910[k];

        t_1142[k] = pa_x[k] * isl0_1142[k]
                    + f_17 * isk_917[k]
                    - f_14 * pc_x[k] * isl1_1142[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, pa_x, pc_x, pc_y, isl0_1143, isl0_1145, \
                         isk_698, isk_918, isk_920, isl1_1143, isl1_1145, \
                         ksk_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = pa_x[k] * isl0_1143[k]
                    + f_17 * isk_918[k]
                    - f_14 * pc_x[k] * isl1_1143[k];

        t_1144[k] = f_16 * isk_698[k]
                    + f_3 * pc_y[k] * ksk_914[k];

        t_1145[k] = pa_x[k] * isl0_1145[k]
                    + f_17 * isk_920[k]
                    - f_14 * pc_x[k] * isl1_1145[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pa_x, pc_x, pc_z, isl0_1146, isl0_1148, \
                         isk_663, isk_921, isk_923, isl1_1146, isl1_1148, \
                         ksk_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = pa_x[k] * isl0_1146[k]
                    + f_16 * isk_921[k]
                    - f_14 * pc_x[k] * isl1_1146[k];

        t_1147[k] = f_18 * isk_663[k]
                    + f_3 * pc_z[k] * ksk_915[k];

        t_1148[k] = pa_x[k] * isl0_1148[k]
                    + f_16 * isk_923[k]
                    - f_14 * pc_x[k] * isl1_1148[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, pa_x, pc_x, pc_y, isl0_1149, isl0_1150, \
                         isk_704, isk_924, isk_925, isl1_1149, isl1_1150, \
                         ksk_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = pa_x[k] * isl0_1149[k]
                    + f_16 * isk_924[k]
                    - f_14 * pc_x[k] * isl1_1149[k];

        t_1150[k] = pa_x[k] * isl0_1150[k]
                    + f_16 * isk_925[k]
                    - f_14 * pc_x[k] * isl1_1150[k];

        t_1151[k] = f_16 * isk_704[k]
                    + f_3 * pc_y[k] * ksk_920[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, t_1155, pa_x, pc_x, isl0_1152, isk_927, \
                         isk_928, isk_929, isk_930, isl1_1152, ksk_928, ksk_929, \
                         ksk_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = pa_x[k] * isl0_1152[k]
                    + f_16 * isk_927[k]
                    - f_14 * pc_x[k] * isl1_1152[k];

        t_1153[k] = f_15 * isk_928[k]
                    + f_3 * pc_x[k] * ksk_928[k];

        t_1154[k] = f_15 * isk_929[k]
                    + f_3 * pc_x[k] * ksk_929[k];

        t_1155[k] = f_15 * isk_930[k]
                    + f_3 * pc_x[k] * ksk_930[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t isl0,
                                                           const size_t isk, const size_t isl1,
                                                           const size_t ksi0, const size_t ksi1,
                                                           const size_t ksk, const size_t ncols,
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
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);
    const auto f_23 = 4.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isl0_900 = buffer.data(isl0 + 900);
    const auto *isl0_905 = buffer.data(isl0 + 905);
    const auto *isl0_909 = buffer.data(isl0 + 909);
    const auto *isl0_914 = buffer.data(isl0 + 914);
    const auto *isl0_920 = buffer.data(isl0 + 920);
    const auto *isl0_927 = buffer.data(isl0 + 927);
    const auto *isl0_1161 = buffer.data(isl0 + 1161);
    const auto *isl0_1163 = buffer.data(isl0 + 1163);
    const auto *isl0_1164 = buffer.data(isl0 + 1164);
    const auto *isl0_1165 = buffer.data(isl0 + 1165);
    const auto *isl0_1166 = buffer.data(isl0 + 1166);
    const auto *isl0_1167 = buffer.data(isl0 + 1167);
    const auto *isl0_1169 = buffer.data(isl0 + 1169);
    const auto *isl0_1173 = buffer.data(isl0 + 1173);
    const auto *isl0_1176 = buffer.data(isl0 + 1176);
    const auto *isl0_1180 = buffer.data(isl0 + 1180);
    const auto *isl0_1182 = buffer.data(isl0 + 1182);
    const auto *isl0_1185 = buffer.data(isl0 + 1185);
    const auto *isl0_1187 = buffer.data(isl0 + 1187);
    const auto *isl0_1188 = buffer.data(isl0 + 1188);
    const auto *isl0_1191 = buffer.data(isl0 + 1191);
    const auto *isl0_1193 = buffer.data(isl0 + 1193);
    const auto *isl0_1194 = buffer.data(isl0 + 1194);
    const auto *isl0_1195 = buffer.data(isl0 + 1195);
    const auto *isl0_1206 = buffer.data(isl0 + 1206);
    const auto *isl0_1208 = buffer.data(isl0 + 1208);
    const auto *isl0_1209 = buffer.data(isl0 + 1209);
    const auto *isl0_1210 = buffer.data(isl0 + 1210);
    const auto *isl0_1211 = buffer.data(isl0 + 1211);
    const auto *isl0_1212 = buffer.data(isl0 + 1212);
    const auto *isl0_1214 = buffer.data(isl0 + 1214);
    const auto *isl0_1215 = buffer.data(isl0 + 1215);
    const auto *isl0_1220 = buffer.data(isl0 + 1220);
    const auto *isl0_1224 = buffer.data(isl0 + 1224);
    const auto *isl0_1229 = buffer.data(isl0 + 1229);
    const auto *isl0_1235 = buffer.data(isl0 + 1235);
    const auto *isl0_1242 = buffer.data(isl0 + 1242);
    const auto *isl0_1251 = buffer.data(isl0 + 1251);
    const auto *isl0_1252 = buffer.data(isl0 + 1252);
    const auto *isl0_1253 = buffer.data(isl0 + 1253);
    const auto *isl0_1254 = buffer.data(isl0 + 1254);
    const auto *isl0_1255 = buffer.data(isl0 + 1255);
    const auto *isl0_1256 = buffer.data(isl0 + 1256);
    const auto *isl0_1257 = buffer.data(isl0 + 1257);
    const auto *isl0_1259 = buffer.data(isl0 + 1259);

    const auto *isk_676 = buffer.data(isk + 676);
    const auto *isk_684 = buffer.data(isk + 684);
    const auto *isk_687 = buffer.data(isk + 687);
    const auto *isk_690 = buffer.data(isk + 690);
    const auto *isk_694 = buffer.data(isk + 694);
    const auto *isk_699 = buffer.data(isk + 699);
    const auto *isk_712 = buffer.data(isk + 712);
    const auto *isk_719 = buffer.data(isk + 719);
    const auto *isk_720 = buffer.data(isk + 720);
    const auto *isk_722 = buffer.data(isk + 722);
    const auto *isk_725 = buffer.data(isk + 725);
    const auto *isk_729 = buffer.data(isk + 729);
    const auto *isk_734 = buffer.data(isk + 734);
    const auto *isk_740 = buffer.data(isk + 740);
    const auto *isk_755 = buffer.data(isk + 755);
    const auto *isk_931 = buffer.data(isk + 931);
    const auto *isk_932 = buffer.data(isk + 932);
    const auto *isk_933 = buffer.data(isk + 933);
    const auto *isk_934 = buffer.data(isk + 934);
    const auto *isk_935 = buffer.data(isk + 935);
    const auto *isk_939 = buffer.data(isk + 939);
    const auto *isk_942 = buffer.data(isk + 942);
    const auto *isk_946 = buffer.data(isk + 946);
    const auto *isk_948 = buffer.data(isk + 948);
    const auto *isk_951 = buffer.data(isk + 951);
    const auto *isk_953 = buffer.data(isk + 953);
    const auto *isk_954 = buffer.data(isk + 954);
    const auto *isk_957 = buffer.data(isk + 957);
    const auto *isk_959 = buffer.data(isk + 959);
    const auto *isk_960 = buffer.data(isk + 960);
    const auto *isk_961 = buffer.data(isk + 961);
    const auto *isk_964 = buffer.data(isk + 964);
    const auto *isk_965 = buffer.data(isk + 965);
    const auto *isk_966 = buffer.data(isk + 966);
    const auto *isk_967 = buffer.data(isk + 967);
    const auto *isk_968 = buffer.data(isk + 968);
    const auto *isk_969 = buffer.data(isk + 969);
    const auto *isk_970 = buffer.data(isk + 970);
    const auto *isk_971 = buffer.data(isk + 971);
    const auto *isk_972 = buffer.data(isk + 972);
    const auto *isk_977 = buffer.data(isk + 977);
    const auto *isk_981 = buffer.data(isk + 981);
    const auto *isk_986 = buffer.data(isk + 986);
    const auto *isk_992 = buffer.data(isk + 992);
    const auto *isk_999 = buffer.data(isk + 999);
    const auto *isk_1000 = buffer.data(isk + 1000);
    const auto *isk_1001 = buffer.data(isk + 1001);
    const auto *isk_1002 = buffer.data(isk + 1002);
    const auto *isk_1003 = buffer.data(isk + 1003);
    const auto *isk_1004 = buffer.data(isk + 1004);
    const auto *isk_1005 = buffer.data(isk + 1005);
    const auto *isk_1007 = buffer.data(isk + 1007);

    const auto *isl1_900 = buffer.data(isl1 + 900);
    const auto *isl1_905 = buffer.data(isl1 + 905);
    const auto *isl1_909 = buffer.data(isl1 + 909);
    const auto *isl1_914 = buffer.data(isl1 + 914);
    const auto *isl1_920 = buffer.data(isl1 + 920);
    const auto *isl1_927 = buffer.data(isl1 + 927);
    const auto *isl1_1161 = buffer.data(isl1 + 1161);
    const auto *isl1_1163 = buffer.data(isl1 + 1163);
    const auto *isl1_1164 = buffer.data(isl1 + 1164);
    const auto *isl1_1165 = buffer.data(isl1 + 1165);
    const auto *isl1_1166 = buffer.data(isl1 + 1166);
    const auto *isl1_1167 = buffer.data(isl1 + 1167);
    const auto *isl1_1169 = buffer.data(isl1 + 1169);
    const auto *isl1_1173 = buffer.data(isl1 + 1173);
    const auto *isl1_1176 = buffer.data(isl1 + 1176);
    const auto *isl1_1180 = buffer.data(isl1 + 1180);
    const auto *isl1_1182 = buffer.data(isl1 + 1182);
    const auto *isl1_1185 = buffer.data(isl1 + 1185);
    const auto *isl1_1187 = buffer.data(isl1 + 1187);
    const auto *isl1_1188 = buffer.data(isl1 + 1188);
    const auto *isl1_1191 = buffer.data(isl1 + 1191);
    const auto *isl1_1193 = buffer.data(isl1 + 1193);
    const auto *isl1_1194 = buffer.data(isl1 + 1194);
    const auto *isl1_1195 = buffer.data(isl1 + 1195);
    const auto *isl1_1206 = buffer.data(isl1 + 1206);
    const auto *isl1_1208 = buffer.data(isl1 + 1208);
    const auto *isl1_1209 = buffer.data(isl1 + 1209);
    const auto *isl1_1210 = buffer.data(isl1 + 1210);
    const auto *isl1_1211 = buffer.data(isl1 + 1211);
    const auto *isl1_1212 = buffer.data(isl1 + 1212);
    const auto *isl1_1214 = buffer.data(isl1 + 1214);
    const auto *isl1_1215 = buffer.data(isl1 + 1215);
    const auto *isl1_1220 = buffer.data(isl1 + 1220);
    const auto *isl1_1224 = buffer.data(isl1 + 1224);
    const auto *isl1_1229 = buffer.data(isl1 + 1229);
    const auto *isl1_1235 = buffer.data(isl1 + 1235);
    const auto *isl1_1242 = buffer.data(isl1 + 1242);
    const auto *isl1_1251 = buffer.data(isl1 + 1251);
    const auto *isl1_1252 = buffer.data(isl1 + 1252);
    const auto *isl1_1253 = buffer.data(isl1 + 1253);
    const auto *isl1_1254 = buffer.data(isl1 + 1254);
    const auto *isl1_1255 = buffer.data(isl1 + 1255);
    const auto *isl1_1256 = buffer.data(isl1 + 1256);
    const auto *isl1_1257 = buffer.data(isl1 + 1257);
    const auto *isl1_1259 = buffer.data(isl1 + 1259);

    const auto *ksi0_756 = buffer.data(ksi0 + 756);
    const auto *ksi0_757 = buffer.data(ksi0 + 757);
    const auto *ksi0_758 = buffer.data(ksi0 + 758);
    const auto *ksi0_759 = buffer.data(ksi0 + 759);
    const auto *ksi0_760 = buffer.data(ksi0 + 760);
    const auto *ksi0_761 = buffer.data(ksi0 + 761);
    const auto *ksi0_762 = buffer.data(ksi0 + 762);
    const auto *ksi0_763 = buffer.data(ksi0 + 763);
    const auto *ksi0_764 = buffer.data(ksi0 + 764);
    const auto *ksi0_765 = buffer.data(ksi0 + 765);
    const auto *ksi0_766 = buffer.data(ksi0 + 766);
    const auto *ksi0_767 = buffer.data(ksi0 + 767);
    const auto *ksi0_768 = buffer.data(ksi0 + 768);
    const auto *ksi0_769 = buffer.data(ksi0 + 769);
    const auto *ksi0_770 = buffer.data(ksi0 + 770);
    const auto *ksi0_784 = buffer.data(ksi0 + 784);
    const auto *ksi0_785 = buffer.data(ksi0 + 785);
    const auto *ksi0_787 = buffer.data(ksi0 + 787);
    const auto *ksi0_789 = buffer.data(ksi0 + 789);
    const auto *ksi0_790 = buffer.data(ksi0 + 790);
    const auto *ksi0_792 = buffer.data(ksi0 + 792);
    const auto *ksi0_793 = buffer.data(ksi0 + 793);
    const auto *ksi0_794 = buffer.data(ksi0 + 794);
    const auto *ksi0_796 = buffer.data(ksi0 + 796);
    const auto *ksi0_797 = buffer.data(ksi0 + 797);
    const auto *ksi0_798 = buffer.data(ksi0 + 798);
    const auto *ksi0_799 = buffer.data(ksi0 + 799);
    const auto *ksi0_801 = buffer.data(ksi0 + 801);
    const auto *ksi0_802 = buffer.data(ksi0 + 802);

    const auto *ksi1_756 = buffer.data(ksi1 + 756);
    const auto *ksi1_757 = buffer.data(ksi1 + 757);
    const auto *ksi1_758 = buffer.data(ksi1 + 758);
    const auto *ksi1_759 = buffer.data(ksi1 + 759);
    const auto *ksi1_760 = buffer.data(ksi1 + 760);
    const auto *ksi1_761 = buffer.data(ksi1 + 761);
    const auto *ksi1_762 = buffer.data(ksi1 + 762);
    const auto *ksi1_763 = buffer.data(ksi1 + 763);
    const auto *ksi1_764 = buffer.data(ksi1 + 764);
    const auto *ksi1_765 = buffer.data(ksi1 + 765);
    const auto *ksi1_766 = buffer.data(ksi1 + 766);
    const auto *ksi1_767 = buffer.data(ksi1 + 767);
    const auto *ksi1_768 = buffer.data(ksi1 + 768);
    const auto *ksi1_769 = buffer.data(ksi1 + 769);
    const auto *ksi1_770 = buffer.data(ksi1 + 770);
    const auto *ksi1_784 = buffer.data(ksi1 + 784);
    const auto *ksi1_785 = buffer.data(ksi1 + 785);
    const auto *ksi1_787 = buffer.data(ksi1 + 787);
    const auto *ksi1_789 = buffer.data(ksi1 + 789);
    const auto *ksi1_790 = buffer.data(ksi1 + 790);
    const auto *ksi1_792 = buffer.data(ksi1 + 792);
    const auto *ksi1_793 = buffer.data(ksi1 + 793);
    const auto *ksi1_794 = buffer.data(ksi1 + 794);
    const auto *ksi1_796 = buffer.data(ksi1 + 796);
    const auto *ksi1_797 = buffer.data(ksi1 + 797);
    const auto *ksi1_798 = buffer.data(ksi1 + 798);
    const auto *ksi1_799 = buffer.data(ksi1 + 799);
    const auto *ksi1_801 = buffer.data(ksi1 + 801);
    const auto *ksi1_802 = buffer.data(ksi1 + 802);

    const auto *ksk_928 = buffer.data(ksk + 928);
    const auto *ksk_931 = buffer.data(ksk + 931);
    const auto *ksk_932 = buffer.data(ksk + 932);
    const auto *ksk_933 = buffer.data(ksk + 933);
    const auto *ksk_934 = buffer.data(ksk + 934);
    const auto *ksk_935 = buffer.data(ksk + 935);
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
    const auto *ksk_965 = buffer.data(ksk + 965);
    const auto *ksk_966 = buffer.data(ksk + 966);
    const auto *ksk_967 = buffer.data(ksk + 967);
    const auto *ksk_968 = buffer.data(ksk + 968);
    const auto *ksk_969 = buffer.data(ksk + 969);
    const auto *ksk_970 = buffer.data(ksk + 970);
    const auto *ksk_971 = buffer.data(ksk + 971);
    const auto *ksk_972 = buffer.data(ksk + 972);
    const auto *ksk_973 = buffer.data(ksk + 973);
    const auto *ksk_974 = buffer.data(ksk + 974);
    const auto *ksk_975 = buffer.data(ksk + 975);
    const auto *ksk_976 = buffer.data(ksk + 976);
    const auto *ksk_977 = buffer.data(ksk + 977);
    const auto *ksk_978 = buffer.data(ksk + 978);
    const auto *ksk_979 = buffer.data(ksk + 979);
    const auto *ksk_980 = buffer.data(ksk + 980);
    const auto *ksk_981 = buffer.data(ksk + 981);
    const auto *ksk_982 = buffer.data(ksk + 982);
    const auto *ksk_983 = buffer.data(ksk + 983);
    const auto *ksk_984 = buffer.data(ksk + 984);
    const auto *ksk_985 = buffer.data(ksk + 985);
    const auto *ksk_986 = buffer.data(ksk + 986);
    const auto *ksk_987 = buffer.data(ksk + 987);
    const auto *ksk_988 = buffer.data(ksk + 988);
    const auto *ksk_989 = buffer.data(ksk + 989);
    const auto *ksk_990 = buffer.data(ksk + 990);
    const auto *ksk_991 = buffer.data(ksk + 991);
    const auto *ksk_992 = buffer.data(ksk + 992);
    const auto *ksk_999 = buffer.data(ksk + 999);
    const auto *ksk_1000 = buffer.data(ksk + 1000);
    const auto *ksk_1001 = buffer.data(ksk + 1001);
    const auto *ksk_1002 = buffer.data(ksk + 1002);
    const auto *ksk_1003 = buffer.data(ksk + 1003);
    const auto *ksk_1004 = buffer.data(ksk + 1004);
    const auto *ksk_1005 = buffer.data(ksk + 1005);
    const auto *ksk_1007 = buffer.data(ksk + 1007);
    const auto *ksk_1008 = buffer.data(ksk + 1008);
    const auto *ksk_1009 = buffer.data(ksk + 1009);
    const auto *ksk_1011 = buffer.data(ksk + 1011);
    const auto *ksk_1013 = buffer.data(ksk + 1013);
    const auto *ksk_1014 = buffer.data(ksk + 1014);
    const auto *ksk_1016 = buffer.data(ksk + 1016);
    const auto *ksk_1017 = buffer.data(ksk + 1017);
    const auto *ksk_1018 = buffer.data(ksk + 1018);
    const auto *ksk_1020 = buffer.data(ksk + 1020);
    const auto *ksk_1021 = buffer.data(ksk + 1021);
    const auto *ksk_1022 = buffer.data(ksk + 1022);
    const auto *ksk_1023 = buffer.data(ksk + 1023);
    const auto *ksk_1025 = buffer.data(ksk + 1025);
    const auto *ksk_1026 = buffer.data(ksk + 1026);

#pragma omp simd aligned(t_1156, t_1157, t_1158, t_1159, t_1160, pc_x, isk_931, isk_932, \
                         isk_933, isk_934, isk_935, ksk_931, ksk_932, ksk_933, ksk_934, \
                         ksk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_15 * isk_931[k]
                    + f_3 * pc_x[k] * ksk_931[k];

        t_1157[k] = f_15 * isk_932[k]
                    + f_3 * pc_x[k] * ksk_932[k];

        t_1158[k] = f_15 * isk_933[k]
                    + f_3 * pc_x[k] * ksk_933[k];

        t_1159[k] = f_15 * isk_934[k]
                    + f_3 * pc_x[k] * ksk_934[k];

        t_1160[k] = f_15 * isk_935[k]
                    + f_3 * pc_x[k] * ksk_935[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, pa_x, pc_x, pc_z, isl0_1161, \
                         isl0_1163, isl0_1164, isk_676, isl1_1161, isl1_1163, isl1_1164, \
                         ksk_928 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = pa_x[k] * isl0_1161[k]
                    - f_14 * pc_x[k] * isl1_1161[k];

        t_1162[k] = f_18 * isk_676[k]
                    + f_3 * pc_z[k] * ksk_928[k];

        t_1163[k] = pa_x[k] * isl0_1163[k]
                    - f_14 * pc_x[k] * isl1_1163[k];

        t_1164[k] = pa_x[k] * isl0_1164[k]
                    - f_14 * pc_x[k] * isl1_1164[k];
    }

#pragma omp simd aligned(t_1165, t_1166, t_1167, t_1168, pa_x, pc_x, pc_y, isl0_1165, \
                         isl0_1166, isl0_1167, isk_719, isl1_1165, isl1_1166, isl1_1167, \
                         ksk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1165[k] = pa_x[k] * isl0_1165[k]
                    - f_14 * pc_x[k] * isl1_1165[k];

        t_1166[k] = pa_x[k] * isl0_1166[k]
                    - f_14 * pc_x[k] * isl1_1166[k];

        t_1167[k] = pa_x[k] * isl0_1167[k]
                    - f_14 * pc_x[k] * isl1_1167[k];

        t_1168[k] = f_16 * isk_719[k]
                    + f_3 * pc_y[k] * ksk_935[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, pa_x, pa_y, pc_x, pc_y, pc_z, \
                         isl0_900, isl0_1169, isk_684, isk_720, isl1_900, isl1_1169, \
                         ksk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = pa_x[k] * isl0_1169[k]
                    - f_14 * pc_x[k] * isl1_1169[k];

        t_1170[k] = pa_y[k] * isl0_900[k]
                    - f_14 * pc_y[k] * isl1_900[k];

        t_1171[k] = f_15 * isk_720[k]
                    + f_3 * pc_y[k] * ksk_936[k];

        t_1172[k] = f_19 * isk_684[k]
                    + f_3 * pc_z[k] * ksk_936[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, pa_x, pa_y, pc_x, pc_y, isl0_905, isl0_1173, \
                         isk_722, isk_939, isl1_905, isl1_1173, \
                         ksk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = pa_x[k] * isl0_1173[k]
                    + f_20 * isk_939[k]
                    - f_14 * pc_x[k] * isl1_1173[k];

        t_1174[k] = f_15 * isk_722[k]
                    + f_3 * pc_y[k] * ksk_938[k];

        t_1175[k] = pa_y[k] * isl0_905[k]
                    - f_14 * pc_y[k] * isl1_905[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, pa_x, pc_x, pc_y, pc_z, isl0_1176, isk_687, \
                         isk_725, isk_942, isl1_1176, ksk_939, \
                         ksk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = pa_x[k] * isl0_1176[k]
                    + f_19 * isk_942[k]
                    - f_14 * pc_x[k] * isl1_1176[k];

        t_1177[k] = f_19 * isk_687[k]
                    + f_3 * pc_z[k] * ksk_939[k];

        t_1178[k] = f_15 * isk_725[k]
                    + f_3 * pc_y[k] * ksk_941[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pa_x, pa_y, pc_x, pc_y, pc_z, isl0_909, \
                         isl0_1180, isk_690, isk_946, isl1_909, isl1_1180, \
                         ksk_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = pa_y[k] * isl0_909[k]
                    - f_14 * pc_y[k] * isl1_909[k];

        t_1180[k] = pa_x[k] * isl0_1180[k]
                    + f_18 * isk_946[k]
                    - f_14 * pc_x[k] * isl1_1180[k];

        t_1181[k] = f_19 * isk_690[k]
                    + f_3 * pc_z[k] * ksk_942[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pa_x, pa_y, pc_x, pc_y, isl0_914, isl0_1182, \
                         isk_729, isk_948, isl1_914, isl1_1182, \
                         ksk_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = pa_x[k] * isl0_1182[k]
                    + f_18 * isk_948[k]
                    - f_14 * pc_x[k] * isl1_1182[k];

        t_1183[k] = f_15 * isk_729[k]
                    + f_3 * pc_y[k] * ksk_945[k];

        t_1184[k] = pa_y[k] * isl0_914[k]
                    - f_14 * pc_y[k] * isl1_914[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pa_x, pc_x, pc_z, isl0_1185, isl0_1187, \
                         isk_694, isk_951, isk_953, isl1_1185, isl1_1187, \
                         ksk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = pa_x[k] * isl0_1185[k]
                    + f_17 * isk_951[k]
                    - f_14 * pc_x[k] * isl1_1185[k];

        t_1186[k] = f_19 * isk_694[k]
                    + f_3 * pc_z[k] * ksk_946[k];

        t_1187[k] = pa_x[k] * isl0_1187[k]
                    + f_17 * isk_953[k]
                    - f_14 * pc_x[k] * isl1_1187[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, pa_x, pa_y, pc_x, pc_y, isl0_920, isl0_1188, \
                         isk_734, isk_954, isl1_920, isl1_1188, \
                         ksk_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = pa_x[k] * isl0_1188[k]
                    + f_17 * isk_954[k]
                    - f_14 * pc_x[k] * isl1_1188[k];

        t_1189[k] = f_15 * isk_734[k]
                    + f_3 * pc_y[k] * ksk_950[k];

        t_1190[k] = pa_y[k] * isl0_920[k]
                    - f_14 * pc_y[k] * isl1_920[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, pa_x, pc_x, pc_z, isl0_1191, isl0_1193, \
                         isk_699, isk_957, isk_959, isl1_1191, isl1_1193, \
                         ksk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = pa_x[k] * isl0_1191[k]
                    + f_16 * isk_957[k]
                    - f_14 * pc_x[k] * isl1_1191[k];

        t_1192[k] = f_19 * isk_699[k]
                    + f_3 * pc_z[k] * ksk_951[k];

        t_1193[k] = pa_x[k] * isl0_1193[k]
                    + f_16 * isk_959[k]
                    - f_14 * pc_x[k] * isl1_1193[k];
    }

#pragma omp simd aligned(t_1194, t_1195, t_1196, pa_x, pc_x, pc_y, isl0_1194, isl0_1195, \
                         isk_740, isk_960, isk_961, isl1_1194, isl1_1195, \
                         ksk_956 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1194[k] = pa_x[k] * isl0_1194[k]
                    + f_16 * isk_960[k]
                    - f_14 * pc_x[k] * isl1_1194[k];

        t_1195[k] = pa_x[k] * isl0_1195[k]
                    + f_16 * isk_961[k]
                    - f_14 * pc_x[k] * isl1_1195[k];

        t_1196[k] = f_15 * isk_740[k]
                    + f_3 * pc_y[k] * ksk_956[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, t_1200, pa_y, pc_x, pc_y, isl0_927, isk_964, \
                         isk_965, isk_966, isl1_927, ksk_964, ksk_965, \
                         ksk_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = pa_y[k] * isl0_927[k]
                    - f_14 * pc_y[k] * isl1_927[k];

        t_1198[k] = f_15 * isk_964[k]
                    + f_3 * pc_x[k] * ksk_964[k];

        t_1199[k] = f_15 * isk_965[k]
                    + f_3 * pc_x[k] * ksk_965[k];

        t_1200[k] = f_15 * isk_966[k]
                    + f_3 * pc_x[k] * ksk_966[k];
    }

#pragma omp simd aligned(t_1201, t_1202, t_1203, t_1204, t_1205, pc_x, isk_967, isk_968, \
                         isk_969, isk_970, isk_971, ksk_967, ksk_968, ksk_969, ksk_970, \
                         ksk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1201[k] = f_15 * isk_967[k]
                    + f_3 * pc_x[k] * ksk_967[k];

        t_1202[k] = f_15 * isk_968[k]
                    + f_3 * pc_x[k] * ksk_968[k];

        t_1203[k] = f_15 * isk_969[k]
                    + f_3 * pc_x[k] * ksk_969[k];

        t_1204[k] = f_15 * isk_970[k]
                    + f_3 * pc_x[k] * ksk_970[k];

        t_1205[k] = f_15 * isk_971[k]
                    + f_3 * pc_x[k] * ksk_971[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pa_x, pc_x, pc_z, isl0_1206, \
                         isl0_1208, isl0_1209, isk_712, isl1_1206, isl1_1208, isl1_1209, \
                         ksk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = pa_x[k] * isl0_1206[k]
                    - f_14 * pc_x[k] * isl1_1206[k];

        t_1207[k] = f_19 * isk_712[k]
                    + f_3 * pc_z[k] * ksk_964[k];

        t_1208[k] = pa_x[k] * isl0_1208[k]
                    - f_14 * pc_x[k] * isl1_1208[k];

        t_1209[k] = pa_x[k] * isl0_1209[k]
                    - f_14 * pc_x[k] * isl1_1209[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pa_x, pc_x, pc_y, isl0_1210, \
                         isl0_1211, isl0_1212, isk_755, isl1_1210, isl1_1211, isl1_1212, \
                         ksk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = pa_x[k] * isl0_1210[k]
                    - f_14 * pc_x[k] * isl1_1210[k];

        t_1211[k] = pa_x[k] * isl0_1211[k]
                    - f_14 * pc_x[k] * isl1_1211[k];

        t_1212[k] = pa_x[k] * isl0_1212[k]
                    - f_14 * pc_x[k] * isl1_1212[k];

        t_1213[k] = f_15 * isk_755[k]
                    + f_3 * pc_y[k] * ksk_971[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pa_x, pc_x, pc_y, pc_z, isl0_1214, \
                         isl0_1215, isk_720, isk_972, isl1_1214, isl1_1215, \
                         ksk_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pa_x[k] * isl0_1214[k]
                    - f_14 * pc_x[k] * isl1_1214[k];

        t_1215[k] = pa_x[k] * isl0_1215[k]
                    + f_23 * isk_972[k]
                    - f_14 * pc_x[k] * isl1_1215[k];

        t_1216[k] = f_3 * pc_y[k] * ksk_972[k];

        t_1217[k] = f_20 * isk_720[k]
                    + f_3 * pc_z[k] * ksk_972[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, pa_x, pc_x, pc_y, isl0_1220, isk_977, \
                         isl1_1220, ksi0_756, ksi1_756, ksk_973, \
                         ksk_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = f_4 * ksi0_756[k]
                    - f_5 * ksi1_756[k]
                    + f_3 * pc_y[k] * ksk_973[k];

        t_1219[k] = f_3 * pc_y[k] * ksk_974[k];

        t_1220[k] = pa_x[k] * isl0_1220[k]
                    + f_20 * isk_977[k]
                    - f_14 * pc_x[k] * isl1_1220[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, pc_y, ksi0_757, ksi0_758, ksi1_757, ksi1_758, \
                         ksk_975, ksk_976, ksk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_6 * ksi0_757[k]
                    - f_7 * ksi1_757[k]
                    + f_3 * pc_y[k] * ksk_975[k];

        t_1222[k] = f_4 * ksi0_758[k]
                    - f_5 * ksi1_758[k]
                    + f_3 * pc_y[k] * ksk_976[k];

        t_1223[k] = f_3 * pc_y[k] * ksk_977[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pa_x, pc_x, pc_y, isl0_1224, isk_981, \
                         isl1_1224, ksi0_759, ksi0_760, ksi1_759, ksi1_760, ksk_978, \
                         ksk_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = pa_x[k] * isl0_1224[k]
                    + f_19 * isk_981[k]
                    - f_14 * pc_x[k] * isl1_1224[k];

        t_1225[k] = f_8 * ksi0_759[k]
                    - f_9 * ksi1_759[k]
                    + f_3 * pc_y[k] * ksk_978[k];

        t_1226[k] = f_6 * ksi0_760[k]
                    - f_7 * ksi1_760[k]
                    + f_3 * pc_y[k] * ksk_979[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pa_x, pc_x, pc_y, isl0_1229, isk_986, \
                         isl1_1229, ksi0_761, ksi1_761, ksk_980, \
                         ksk_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_4 * ksi0_761[k]
                    - f_5 * ksi1_761[k]
                    + f_3 * pc_y[k] * ksk_980[k];

        t_1228[k] = f_3 * pc_y[k] * ksk_981[k];

        t_1229[k] = pa_x[k] * isl0_1229[k]
                    + f_18 * isk_986[k]
                    - f_14 * pc_x[k] * isl1_1229[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, pc_y, ksi0_762, ksi0_763, ksi0_764, ksi1_762, \
                         ksi1_763, ksi1_764, ksk_982, ksk_983, \
                         ksk_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_10 * ksi0_762[k]
                    - f_11 * ksi1_762[k]
                    + f_3 * pc_y[k] * ksk_982[k];

        t_1231[k] = f_8 * ksi0_763[k]
                    - f_9 * ksi1_763[k]
                    + f_3 * pc_y[k] * ksk_983[k];

        t_1232[k] = f_6 * ksi0_764[k]
                    - f_7 * ksi1_764[k]
                    + f_3 * pc_y[k] * ksk_984[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, pa_x, pc_x, pc_y, isl0_1235, isk_992, \
                         isl1_1235, ksi0_765, ksi1_765, ksk_985, \
                         ksk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = f_4 * ksi0_765[k]
                    - f_5 * ksi1_765[k]
                    + f_3 * pc_y[k] * ksk_985[k];

        t_1234[k] = f_3 * pc_y[k] * ksk_986[k];

        t_1235[k] = pa_x[k] * isl0_1235[k]
                    + f_17 * isk_992[k]
                    - f_14 * pc_x[k] * isl1_1235[k];
    }

#pragma omp simd aligned(t_1236, t_1237, t_1238, pc_y, ksi0_766, ksi0_767, ksi0_768, ksi1_766, \
                         ksi1_767, ksi1_768, ksk_987, ksk_988, \
                         ksk_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1236[k] = f_12 * ksi0_766[k]
                    - f_13 * ksi1_766[k]
                    + f_3 * pc_y[k] * ksk_987[k];

        t_1237[k] = f_10 * ksi0_767[k]
                    - f_11 * ksi1_767[k]
                    + f_3 * pc_y[k] * ksk_988[k];

        t_1238[k] = f_8 * ksi0_768[k]
                    - f_9 * ksi1_768[k]
                    + f_3 * pc_y[k] * ksk_989[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, pc_y, ksi0_769, ksi0_770, ksi1_769, ksi1_770, \
                         ksk_990, ksk_991, ksk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_6 * ksi0_769[k]
                    - f_7 * ksi1_769[k]
                    + f_3 * pc_y[k] * ksk_990[k];

        t_1240[k] = f_4 * ksi0_770[k]
                    - f_5 * ksi1_770[k]
                    + f_3 * pc_y[k] * ksk_991[k];

        t_1241[k] = f_3 * pc_y[k] * ksk_992[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, pa_x, pc_x, isl0_1242, isk_999, \
                         isk_1000, isk_1001, isk_1002, isl1_1242, ksk_1000, ksk_1001, \
                         ksk_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = pa_x[k] * isl0_1242[k]
                    + f_16 * isk_999[k]
                    - f_14 * pc_x[k] * isl1_1242[k];

        t_1243[k] = f_15 * isk_1000[k]
                    + f_3 * pc_x[k] * ksk_1000[k];

        t_1244[k] = f_15 * isk_1001[k]
                    + f_3 * pc_x[k] * ksk_1001[k];

        t_1245[k] = f_15 * isk_1002[k]
                    + f_3 * pc_x[k] * ksk_1002[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, t_1250, pc_x, pc_y, isk_1003, \
                         isk_1004, isk_1005, isk_1007, ksk_999, ksk_1003, ksk_1004, ksk_1005, \
                         ksk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = f_15 * isk_1003[k]
                    + f_3 * pc_x[k] * ksk_1003[k];

        t_1247[k] = f_15 * isk_1004[k]
                    + f_3 * pc_x[k] * ksk_1004[k];

        t_1248[k] = f_15 * isk_1005[k]
                    + f_3 * pc_x[k] * ksk_1005[k];

        t_1249[k] = f_3 * pc_y[k] * ksk_999[k];

        t_1250[k] = f_15 * isk_1007[k]
                    + f_3 * pc_x[k] * ksk_1007[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, pa_x, pc_x, isl0_1251, isl0_1252, \
                         isl0_1253, isl0_1254, isl1_1251, isl1_1252, isl1_1253, \
                         isl1_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = pa_x[k] * isl0_1251[k]
                    - f_14 * pc_x[k] * isl1_1251[k];

        t_1252[k] = pa_x[k] * isl0_1252[k]
                    - f_14 * pc_x[k] * isl1_1252[k];

        t_1253[k] = pa_x[k] * isl0_1253[k]
                    - f_14 * pc_x[k] * isl1_1253[k];

        t_1254[k] = pa_x[k] * isl0_1254[k]
                    - f_14 * pc_x[k] * isl1_1254[k];
    }

#pragma omp simd aligned(t_1255, t_1256, t_1257, t_1258, pa_x, pc_x, pc_y, isl0_1255, \
                         isl0_1256, isl0_1257, isl1_1255, isl1_1256, isl1_1257, \
                         ksk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1255[k] = pa_x[k] * isl0_1255[k]
                    - f_14 * pc_x[k] * isl1_1255[k];

        t_1256[k] = pa_x[k] * isl0_1256[k]
                    - f_14 * pc_x[k] * isl1_1256[k];

        t_1257[k] = pa_x[k] * isl0_1257[k]
                    - f_14 * pc_x[k] * isl1_1257[k];

        t_1258[k] = f_3 * pc_y[k] * ksk_1007[k];
    }

#pragma omp simd aligned(t_1259, t_1260, t_1261, t_1262, pa_x, pc_x, pc_z, isl0_1259, \
                         isl1_1259, ksi0_784, ksi0_785, ksi1_784, ksi1_785, ksk_1008, \
                         ksk_1009 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1259[k] = pa_x[k] * isl0_1259[k]
                    - f_14 * pc_x[k] * isl1_1259[k];

        t_1260[k] = f_1 * ksi0_784[k]
                    - f_2 * ksi1_784[k]
                    + f_3 * pc_x[k] * ksk_1008[k];

        t_1261[k] = f_21 * ksi0_785[k]
                    - f_22 * ksi1_785[k]
                    + f_3 * pc_x[k] * ksk_1009[k];

        t_1262[k] = f_3 * pc_z[k] * ksk_1008[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, t_1266, pc_x, pc_z, ksi0_787, ksi0_789, \
                         ksi0_790, ksi1_787, ksi1_789, ksi1_790, ksk_1009, ksk_1011, ksk_1013, \
                         ksk_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = f_12 * ksi0_787[k]
                    - f_13 * ksi1_787[k]
                    + f_3 * pc_x[k] * ksk_1011[k];

        t_1264[k] = f_3 * pc_z[k] * ksk_1009[k];

        t_1265[k] = f_12 * ksi0_789[k]
                    - f_13 * ksi1_789[k]
                    + f_3 * pc_x[k] * ksk_1013[k];

        t_1266[k] = f_10 * ksi0_790[k]
                    - f_11 * ksi1_790[k]
                    + f_3 * pc_x[k] * ksk_1014[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, t_1270, pc_x, pc_z, ksi0_792, ksi0_793, \
                         ksi0_794, ksi1_792, ksi1_793, ksi1_794, ksk_1011, ksk_1016, ksk_1017, \
                         ksk_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_3 * pc_z[k] * ksk_1011[k];

        t_1268[k] = f_10 * ksi0_792[k]
                    - f_11 * ksi1_792[k]
                    + f_3 * pc_x[k] * ksk_1016[k];

        t_1269[k] = f_10 * ksi0_793[k]
                    - f_11 * ksi1_793[k]
                    + f_3 * pc_x[k] * ksk_1017[k];

        t_1270[k] = f_8 * ksi0_794[k]
                    - f_9 * ksi1_794[k]
                    + f_3 * pc_x[k] * ksk_1018[k];
    }

#pragma omp simd aligned(t_1271, t_1272, t_1273, t_1274, pc_x, pc_z, ksi0_796, ksi0_797, \
                         ksi0_798, ksi1_796, ksi1_797, ksi1_798, ksk_1014, ksk_1020, ksk_1021, \
                         ksk_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1271[k] = f_3 * pc_z[k] * ksk_1014[k];

        t_1272[k] = f_8 * ksi0_796[k]
                    - f_9 * ksi1_796[k]
                    + f_3 * pc_x[k] * ksk_1020[k];

        t_1273[k] = f_8 * ksi0_797[k]
                    - f_9 * ksi1_797[k]
                    + f_3 * pc_x[k] * ksk_1021[k];

        t_1274[k] = f_8 * ksi0_798[k]
                    - f_9 * ksi1_798[k]
                    + f_3 * pc_x[k] * ksk_1022[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, t_1278, pc_x, pc_z, ksi0_799, ksi0_801, \
                         ksi0_802, ksi1_799, ksi1_801, ksi1_802, ksk_1018, ksk_1023, ksk_1025, \
                         ksk_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_6 * ksi0_799[k]
                    - f_7 * ksi1_799[k]
                    + f_3 * pc_x[k] * ksk_1023[k];

        t_1276[k] = f_3 * pc_z[k] * ksk_1018[k];

        t_1277[k] = f_6 * ksi0_801[k]
                    - f_7 * ksi1_801[k]
                    + f_3 * pc_x[k] * ksk_1025[k];

        t_1278[k] = f_6 * ksi0_802[k]
                    - f_7 * ksi1_802[k]
                    + f_3 * pc_x[k] * ksk_1026[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t isl0,
                                                           const size_t isk, const size_t isl1,
                                                           const size_t ksi0, const size_t ksi1,
                                                           const size_t ksk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isl0_945 = buffer.data(isl0 + 945);
    const auto *isl0_946 = buffer.data(isl0 + 946);
    const auto *isl0_948 = buffer.data(isl0 + 948);
    const auto *isl0_951 = buffer.data(isl0 + 951);
    const auto *isl0_955 = buffer.data(isl0 + 955);
    const auto *isl0_960 = buffer.data(isl0 + 960);
    const auto *isl0_966 = buffer.data(isl0 + 966);
    const auto *isl0_981 = buffer.data(isl0 + 981);
    const auto *isl0_983 = buffer.data(isl0 + 983);
    const auto *isl0_984 = buffer.data(isl0 + 984);
    const auto *isl0_985 = buffer.data(isl0 + 985);
    const auto *isl0_986 = buffer.data(isl0 + 986);
    const auto *isl0_987 = buffer.data(isl0 + 987);

    const auto *isk_784 = buffer.data(isk + 784);
    const auto *isk_785 = buffer.data(isk + 785);
    const auto *isk_786 = buffer.data(isk + 786);
    const auto *isk_787 = buffer.data(isk + 787);
    const auto *isk_788 = buffer.data(isk + 788);
    const auto *isk_789 = buffer.data(isk + 789);
    const auto *isk_791 = buffer.data(isk + 791);
    const auto *isk_820 = buffer.data(isk + 820);
    const auto *isk_827 = buffer.data(isk + 827);
    const auto *isk_856 = buffer.data(isk + 856);
    const auto *isk_858 = buffer.data(isk + 858);
    const auto *isk_859 = buffer.data(isk + 859);
    const auto *isk_860 = buffer.data(isk + 860);
    const auto *isk_861 = buffer.data(isk + 861);
    const auto *isk_862 = buffer.data(isk + 862);
    const auto *isk_863 = buffer.data(isk + 863);

    const auto *isl1_945 = buffer.data(isl1 + 945);
    const auto *isl1_946 = buffer.data(isl1 + 946);
    const auto *isl1_948 = buffer.data(isl1 + 948);
    const auto *isl1_951 = buffer.data(isl1 + 951);
    const auto *isl1_955 = buffer.data(isl1 + 955);
    const auto *isl1_960 = buffer.data(isl1 + 960);
    const auto *isl1_966 = buffer.data(isl1 + 966);
    const auto *isl1_981 = buffer.data(isl1 + 981);
    const auto *isl1_983 = buffer.data(isl1 + 983);
    const auto *isl1_984 = buffer.data(isl1 + 984);
    const auto *isl1_985 = buffer.data(isl1 + 985);
    const auto *isl1_986 = buffer.data(isl1 + 986);
    const auto *isl1_987 = buffer.data(isl1 + 987);

    const auto *ksi0_803 = buffer.data(ksi0 + 803);
    const auto *ksi0_804 = buffer.data(ksi0 + 804);
    const auto *ksi0_805 = buffer.data(ksi0 + 805);
    const auto *ksi0_806 = buffer.data(ksi0 + 806);
    const auto *ksi0_807 = buffer.data(ksi0 + 807);
    const auto *ksi0_808 = buffer.data(ksi0 + 808);
    const auto *ksi0_809 = buffer.data(ksi0 + 809);
    const auto *ksi0_810 = buffer.data(ksi0 + 810);
    const auto *ksi0_811 = buffer.data(ksi0 + 811);
    const auto *ksi0_814 = buffer.data(ksi0 + 814);
    const auto *ksi0_816 = buffer.data(ksi0 + 816);
    const auto *ksi0_817 = buffer.data(ksi0 + 817);
    const auto *ksi0_819 = buffer.data(ksi0 + 819);
    const auto *ksi0_820 = buffer.data(ksi0 + 820);
    const auto *ksi0_821 = buffer.data(ksi0 + 821);
    const auto *ksi0_823 = buffer.data(ksi0 + 823);
    const auto *ksi0_824 = buffer.data(ksi0 + 824);
    const auto *ksi0_825 = buffer.data(ksi0 + 825);
    const auto *ksi0_826 = buffer.data(ksi0 + 826);
    const auto *ksi0_828 = buffer.data(ksi0 + 828);
    const auto *ksi0_829 = buffer.data(ksi0 + 829);
    const auto *ksi0_830 = buffer.data(ksi0 + 830);
    const auto *ksi0_831 = buffer.data(ksi0 + 831);
    const auto *ksi0_832 = buffer.data(ksi0 + 832);
    const auto *ksi0_834 = buffer.data(ksi0 + 834);
    const auto *ksi0_835 = buffer.data(ksi0 + 835);
    const auto *ksi0_836 = buffer.data(ksi0 + 836);
    const auto *ksi0_837 = buffer.data(ksi0 + 837);
    const auto *ksi0_838 = buffer.data(ksi0 + 838);
    const auto *ksi0_839 = buffer.data(ksi0 + 839);
    const auto *ksi0_840 = buffer.data(ksi0 + 840);
    const auto *ksi0_841 = buffer.data(ksi0 + 841);
    const auto *ksi0_842 = buffer.data(ksi0 + 842);
    const auto *ksi0_843 = buffer.data(ksi0 + 843);
    const auto *ksi0_844 = buffer.data(ksi0 + 844);
    const auto *ksi0_845 = buffer.data(ksi0 + 845);
    const auto *ksi0_846 = buffer.data(ksi0 + 846);
    const auto *ksi0_847 = buffer.data(ksi0 + 847);
    const auto *ksi0_848 = buffer.data(ksi0 + 848);
    const auto *ksi0_849 = buffer.data(ksi0 + 849);
    const auto *ksi0_850 = buffer.data(ksi0 + 850);
    const auto *ksi0_851 = buffer.data(ksi0 + 851);
    const auto *ksi0_852 = buffer.data(ksi0 + 852);
    const auto *ksi0_853 = buffer.data(ksi0 + 853);
    const auto *ksi0_854 = buffer.data(ksi0 + 854);
    const auto *ksi0_855 = buffer.data(ksi0 + 855);
    const auto *ksi0_856 = buffer.data(ksi0 + 856);
    const auto *ksi0_857 = buffer.data(ksi0 + 857);
    const auto *ksi0_858 = buffer.data(ksi0 + 858);
    const auto *ksi0_859 = buffer.data(ksi0 + 859);
    const auto *ksi0_860 = buffer.data(ksi0 + 860);
    const auto *ksi0_861 = buffer.data(ksi0 + 861);
    const auto *ksi0_862 = buffer.data(ksi0 + 862);
    const auto *ksi0_863 = buffer.data(ksi0 + 863);
    const auto *ksi0_864 = buffer.data(ksi0 + 864);
    const auto *ksi0_865 = buffer.data(ksi0 + 865);
    const auto *ksi0_866 = buffer.data(ksi0 + 866);
    const auto *ksi0_867 = buffer.data(ksi0 + 867);
    const auto *ksi0_868 = buffer.data(ksi0 + 868);
    const auto *ksi0_869 = buffer.data(ksi0 + 869);
    const auto *ksi0_870 = buffer.data(ksi0 + 870);

    const auto *ksi1_803 = buffer.data(ksi1 + 803);
    const auto *ksi1_804 = buffer.data(ksi1 + 804);
    const auto *ksi1_805 = buffer.data(ksi1 + 805);
    const auto *ksi1_806 = buffer.data(ksi1 + 806);
    const auto *ksi1_807 = buffer.data(ksi1 + 807);
    const auto *ksi1_808 = buffer.data(ksi1 + 808);
    const auto *ksi1_809 = buffer.data(ksi1 + 809);
    const auto *ksi1_810 = buffer.data(ksi1 + 810);
    const auto *ksi1_811 = buffer.data(ksi1 + 811);
    const auto *ksi1_814 = buffer.data(ksi1 + 814);
    const auto *ksi1_816 = buffer.data(ksi1 + 816);
    const auto *ksi1_817 = buffer.data(ksi1 + 817);
    const auto *ksi1_819 = buffer.data(ksi1 + 819);
    const auto *ksi1_820 = buffer.data(ksi1 + 820);
    const auto *ksi1_821 = buffer.data(ksi1 + 821);
    const auto *ksi1_823 = buffer.data(ksi1 + 823);
    const auto *ksi1_824 = buffer.data(ksi1 + 824);
    const auto *ksi1_825 = buffer.data(ksi1 + 825);
    const auto *ksi1_826 = buffer.data(ksi1 + 826);
    const auto *ksi1_828 = buffer.data(ksi1 + 828);
    const auto *ksi1_829 = buffer.data(ksi1 + 829);
    const auto *ksi1_830 = buffer.data(ksi1 + 830);
    const auto *ksi1_831 = buffer.data(ksi1 + 831);
    const auto *ksi1_832 = buffer.data(ksi1 + 832);
    const auto *ksi1_834 = buffer.data(ksi1 + 834);
    const auto *ksi1_835 = buffer.data(ksi1 + 835);
    const auto *ksi1_836 = buffer.data(ksi1 + 836);
    const auto *ksi1_837 = buffer.data(ksi1 + 837);
    const auto *ksi1_838 = buffer.data(ksi1 + 838);
    const auto *ksi1_839 = buffer.data(ksi1 + 839);
    const auto *ksi1_840 = buffer.data(ksi1 + 840);
    const auto *ksi1_841 = buffer.data(ksi1 + 841);
    const auto *ksi1_842 = buffer.data(ksi1 + 842);
    const auto *ksi1_843 = buffer.data(ksi1 + 843);
    const auto *ksi1_844 = buffer.data(ksi1 + 844);
    const auto *ksi1_845 = buffer.data(ksi1 + 845);
    const auto *ksi1_846 = buffer.data(ksi1 + 846);
    const auto *ksi1_847 = buffer.data(ksi1 + 847);
    const auto *ksi1_848 = buffer.data(ksi1 + 848);
    const auto *ksi1_849 = buffer.data(ksi1 + 849);
    const auto *ksi1_850 = buffer.data(ksi1 + 850);
    const auto *ksi1_851 = buffer.data(ksi1 + 851);
    const auto *ksi1_852 = buffer.data(ksi1 + 852);
    const auto *ksi1_853 = buffer.data(ksi1 + 853);
    const auto *ksi1_854 = buffer.data(ksi1 + 854);
    const auto *ksi1_855 = buffer.data(ksi1 + 855);
    const auto *ksi1_856 = buffer.data(ksi1 + 856);
    const auto *ksi1_857 = buffer.data(ksi1 + 857);
    const auto *ksi1_858 = buffer.data(ksi1 + 858);
    const auto *ksi1_859 = buffer.data(ksi1 + 859);
    const auto *ksi1_860 = buffer.data(ksi1 + 860);
    const auto *ksi1_861 = buffer.data(ksi1 + 861);
    const auto *ksi1_862 = buffer.data(ksi1 + 862);
    const auto *ksi1_863 = buffer.data(ksi1 + 863);
    const auto *ksi1_864 = buffer.data(ksi1 + 864);
    const auto *ksi1_865 = buffer.data(ksi1 + 865);
    const auto *ksi1_866 = buffer.data(ksi1 + 866);
    const auto *ksi1_867 = buffer.data(ksi1 + 867);
    const auto *ksi1_868 = buffer.data(ksi1 + 868);
    const auto *ksi1_869 = buffer.data(ksi1 + 869);
    const auto *ksi1_870 = buffer.data(ksi1 + 870);

    const auto *ksk_1023 = buffer.data(ksk + 1023);
    const auto *ksk_1027 = buffer.data(ksk + 1027);
    const auto *ksk_1028 = buffer.data(ksk + 1028);
    const auto *ksk_1029 = buffer.data(ksk + 1029);
    const auto *ksk_1031 = buffer.data(ksk + 1031);
    const auto *ksk_1032 = buffer.data(ksk + 1032);
    const auto *ksk_1033 = buffer.data(ksk + 1033);
    const auto *ksk_1034 = buffer.data(ksk + 1034);
    const auto *ksk_1035 = buffer.data(ksk + 1035);
    const auto *ksk_1036 = buffer.data(ksk + 1036);
    const auto *ksk_1037 = buffer.data(ksk + 1037);
    const auto *ksk_1038 = buffer.data(ksk + 1038);
    const auto *ksk_1039 = buffer.data(ksk + 1039);
    const auto *ksk_1040 = buffer.data(ksk + 1040);
    const auto *ksk_1041 = buffer.data(ksk + 1041);
    const auto *ksk_1042 = buffer.data(ksk + 1042);
    const auto *ksk_1043 = buffer.data(ksk + 1043);
    const auto *ksk_1046 = buffer.data(ksk + 1046);
    const auto *ksk_1048 = buffer.data(ksk + 1048);
    const auto *ksk_1049 = buffer.data(ksk + 1049);
    const auto *ksk_1051 = buffer.data(ksk + 1051);
    const auto *ksk_1052 = buffer.data(ksk + 1052);
    const auto *ksk_1053 = buffer.data(ksk + 1053);
    const auto *ksk_1055 = buffer.data(ksk + 1055);
    const auto *ksk_1056 = buffer.data(ksk + 1056);
    const auto *ksk_1057 = buffer.data(ksk + 1057);
    const auto *ksk_1058 = buffer.data(ksk + 1058);
    const auto *ksk_1060 = buffer.data(ksk + 1060);
    const auto *ksk_1061 = buffer.data(ksk + 1061);
    const auto *ksk_1062 = buffer.data(ksk + 1062);
    const auto *ksk_1063 = buffer.data(ksk + 1063);
    const auto *ksk_1064 = buffer.data(ksk + 1064);
    const auto *ksk_1066 = buffer.data(ksk + 1066);
    const auto *ksk_1067 = buffer.data(ksk + 1067);
    const auto *ksk_1068 = buffer.data(ksk + 1068);
    const auto *ksk_1069 = buffer.data(ksk + 1069);
    const auto *ksk_1070 = buffer.data(ksk + 1070);
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
    const auto *ksk_1081 = buffer.data(ksk + 1081);
    const auto *ksk_1082 = buffer.data(ksk + 1082);
    const auto *ksk_1083 = buffer.data(ksk + 1083);
    const auto *ksk_1084 = buffer.data(ksk + 1084);
    const auto *ksk_1085 = buffer.data(ksk + 1085);
    const auto *ksk_1086 = buffer.data(ksk + 1086);
    const auto *ksk_1087 = buffer.data(ksk + 1087);
    const auto *ksk_1088 = buffer.data(ksk + 1088);
    const auto *ksk_1089 = buffer.data(ksk + 1089);
    const auto *ksk_1090 = buffer.data(ksk + 1090);
    const auto *ksk_1091 = buffer.data(ksk + 1091);
    const auto *ksk_1092 = buffer.data(ksk + 1092);
    const auto *ksk_1093 = buffer.data(ksk + 1093);
    const auto *ksk_1094 = buffer.data(ksk + 1094);
    const auto *ksk_1095 = buffer.data(ksk + 1095);
    const auto *ksk_1096 = buffer.data(ksk + 1096);
    const auto *ksk_1097 = buffer.data(ksk + 1097);
    const auto *ksk_1098 = buffer.data(ksk + 1098);
    const auto *ksk_1099 = buffer.data(ksk + 1099);
    const auto *ksk_1100 = buffer.data(ksk + 1100);
    const auto *ksk_1101 = buffer.data(ksk + 1101);
    const auto *ksk_1102 = buffer.data(ksk + 1102);
    const auto *ksk_1103 = buffer.data(ksk + 1103);
    const auto *ksk_1104 = buffer.data(ksk + 1104);
    const auto *ksk_1105 = buffer.data(ksk + 1105);
    const auto *ksk_1106 = buffer.data(ksk + 1106);
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
    const auto *ksk_1117 = buffer.data(ksk + 1117);
    const auto *ksk_1118 = buffer.data(ksk + 1118);

#pragma omp simd aligned(t_1279, t_1280, t_1281, t_1282, pc_x, pc_z, ksi0_803, ksi0_804, \
                         ksi0_805, ksi1_803, ksi1_804, ksi1_805, ksk_1023, ksk_1027, ksk_1028, \
                         ksk_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1279[k] = f_6 * ksi0_803[k]
                    - f_7 * ksi1_803[k]
                    + f_3 * pc_x[k] * ksk_1027[k];

        t_1280[k] = f_6 * ksi0_804[k]
                    - f_7 * ksi1_804[k]
                    + f_3 * pc_x[k] * ksk_1028[k];

        t_1281[k] = f_4 * ksi0_805[k]
                    - f_5 * ksi1_805[k]
                    + f_3 * pc_x[k] * ksk_1029[k];

        t_1282[k] = f_3 * pc_z[k] * ksk_1023[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, pc_x, ksi0_807, ksi0_808, ksi0_809, ksi1_807, \
                         ksi1_808, ksi1_809, ksk_1031, ksk_1032, \
                         ksk_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = f_4 * ksi0_807[k]
                    - f_5 * ksi1_807[k]
                    + f_3 * pc_x[k] * ksk_1031[k];

        t_1284[k] = f_4 * ksi0_808[k]
                    - f_5 * ksi1_808[k]
                    + f_3 * pc_x[k] * ksk_1032[k];

        t_1285[k] = f_4 * ksi0_809[k]
                    - f_5 * ksi1_809[k]
                    + f_3 * pc_x[k] * ksk_1033[k];
    }

#pragma omp simd aligned(t_1286, t_1287, t_1288, t_1289, t_1290, pc_x, ksi0_810, ksi0_811, \
                         ksi1_810, ksi1_811, ksk_1034, ksk_1035, ksk_1036, ksk_1037, \
                         ksk_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1286[k] = f_4 * ksi0_810[k]
                    - f_5 * ksi1_810[k]
                    + f_3 * pc_x[k] * ksk_1034[k];

        t_1287[k] = f_4 * ksi0_811[k]
                    - f_5 * ksi1_811[k]
                    + f_3 * pc_x[k] * ksk_1035[k];

        t_1288[k] = f_3 * pc_x[k] * ksk_1036[k];

        t_1289[k] = f_3 * pc_x[k] * ksk_1037[k];

        t_1290[k] = f_3 * pc_x[k] * ksk_1038[k];
    }

#pragma omp simd aligned(t_1291, t_1292, t_1293, t_1294, t_1295, pc_x, ksk_1039, ksk_1040, \
                         ksk_1041, ksk_1042, ksk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1291[k] = f_3 * pc_x[k] * ksk_1039[k];

        t_1292[k] = f_3 * pc_x[k] * ksk_1040[k];

        t_1293[k] = f_3 * pc_x[k] * ksk_1041[k];

        t_1294[k] = f_3 * pc_x[k] * ksk_1042[k];

        t_1295[k] = f_3 * pc_x[k] * ksk_1043[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, t_1299, pc_y, pc_z, isk_784, ksi0_805, \
                         ksi0_806, ksi1_805, ksi1_806, ksk_1036, ksk_1037, \
                         ksk_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_0 * isk_784[k]
                    + f_1 * ksi0_805[k]
                    - f_2 * ksi1_805[k]
                    + f_3 * pc_y[k] * ksk_1036[k];

        t_1297[k] = f_3 * pc_z[k] * ksk_1036[k];

        t_1298[k] = f_4 * ksi0_805[k]
                    - f_5 * ksi1_805[k]
                    + f_3 * pc_z[k] * ksk_1037[k];

        t_1299[k] = f_6 * ksi0_806[k]
                    - f_7 * ksi1_806[k]
                    + f_3 * pc_z[k] * ksk_1038[k];
    }

#pragma omp simd aligned(t_1300, t_1301, t_1302, pc_z, ksi0_807, ksi0_808, ksi0_809, ksi1_807, \
                         ksi1_808, ksi1_809, ksk_1039, ksk_1040, \
                         ksk_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1300[k] = f_8 * ksi0_807[k]
                    - f_9 * ksi1_807[k]
                    + f_3 * pc_z[k] * ksk_1039[k];

        t_1301[k] = f_10 * ksi0_808[k]
                    - f_11 * ksi1_808[k]
                    + f_3 * pc_z[k] * ksk_1040[k];

        t_1302[k] = f_12 * ksi0_809[k]
                    - f_13 * ksi1_809[k]
                    + f_3 * pc_z[k] * ksk_1041[k];
    }

#pragma omp simd aligned(t_1303, t_1304, t_1305, t_1306, pa_z, pc_y, pc_z, isl0_945, isl0_946, \
                         isk_791, isl1_945, isl1_946, ksi0_811, ksi1_811, \
                         ksk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1303[k] = f_0 * isk_791[k]
                    + f_3 * pc_y[k] * ksk_1043[k];

        t_1304[k] = f_1 * ksi0_811[k]
                    - f_2 * ksi1_811[k]
                    + f_3 * pc_z[k] * ksk_1043[k];

        t_1305[k] = pa_z[k] * isl0_945[k]
                    - f_14 * pc_z[k] * isl1_945[k];

        t_1306[k] = pa_z[k] * isl0_946[k]
                    - f_14 * pc_z[k] * isl1_946[k];
    }

#pragma omp simd aligned(t_1307, t_1308, t_1309, pa_z, pc_x, pc_z, isl0_948, isl1_948, \
                         ksi0_814, ksi0_816, ksi1_814, ksi1_816, ksk_1046, \
                         ksk_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1307[k] = f_21 * ksi0_814[k]
                    - f_22 * ksi1_814[k]
                    + f_3 * pc_x[k] * ksk_1046[k];

        t_1308[k] = pa_z[k] * isl0_948[k]
                    - f_14 * pc_z[k] * isl1_948[k];

        t_1309[k] = f_12 * ksi0_816[k]
                    - f_13 * ksi1_816[k]
                    + f_3 * pc_x[k] * ksk_1048[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, pa_z, pc_x, pc_z, isl0_951, isl1_951, \
                         ksi0_817, ksi0_819, ksi1_817, ksi1_819, ksk_1049, \
                         ksk_1051 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_12 * ksi0_817[k]
                    - f_13 * ksi1_817[k]
                    + f_3 * pc_x[k] * ksk_1049[k];

        t_1311[k] = pa_z[k] * isl0_951[k]
                    - f_14 * pc_z[k] * isl1_951[k];

        t_1312[k] = f_10 * ksi0_819[k]
                    - f_11 * ksi1_819[k]
                    + f_3 * pc_x[k] * ksk_1051[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, pa_z, pc_x, pc_z, isl0_955, isl1_955, \
                         ksi0_820, ksi0_821, ksi1_820, ksi1_821, ksk_1052, \
                         ksk_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_10 * ksi0_820[k]
                    - f_11 * ksi1_820[k]
                    + f_3 * pc_x[k] * ksk_1052[k];

        t_1314[k] = f_10 * ksi0_821[k]
                    - f_11 * ksi1_821[k]
                    + f_3 * pc_x[k] * ksk_1053[k];

        t_1315[k] = pa_z[k] * isl0_955[k]
                    - f_14 * pc_z[k] * isl1_955[k];
    }

#pragma omp simd aligned(t_1316, t_1317, t_1318, pc_x, ksi0_823, ksi0_824, ksi0_825, ksi1_823, \
                         ksi1_824, ksi1_825, ksk_1055, ksk_1056, \
                         ksk_1057 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1316[k] = f_8 * ksi0_823[k]
                    - f_9 * ksi1_823[k]
                    + f_3 * pc_x[k] * ksk_1055[k];

        t_1317[k] = f_8 * ksi0_824[k]
                    - f_9 * ksi1_824[k]
                    + f_3 * pc_x[k] * ksk_1056[k];

        t_1318[k] = f_8 * ksi0_825[k]
                    - f_9 * ksi1_825[k]
                    + f_3 * pc_x[k] * ksk_1057[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, pa_z, pc_x, pc_z, isl0_960, isl1_960, \
                         ksi0_826, ksi0_828, ksi1_826, ksi1_828, ksk_1058, \
                         ksk_1060 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_8 * ksi0_826[k]
                    - f_9 * ksi1_826[k]
                    + f_3 * pc_x[k] * ksk_1058[k];

        t_1320[k] = pa_z[k] * isl0_960[k]
                    - f_14 * pc_z[k] * isl1_960[k];

        t_1321[k] = f_6 * ksi0_828[k]
                    - f_7 * ksi1_828[k]
                    + f_3 * pc_x[k] * ksk_1060[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, pc_x, ksi0_829, ksi0_830, ksi0_831, ksi1_829, \
                         ksi1_830, ksi1_831, ksk_1061, ksk_1062, \
                         ksk_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = f_6 * ksi0_829[k]
                    - f_7 * ksi1_829[k]
                    + f_3 * pc_x[k] * ksk_1061[k];

        t_1323[k] = f_6 * ksi0_830[k]
                    - f_7 * ksi1_830[k]
                    + f_3 * pc_x[k] * ksk_1062[k];

        t_1324[k] = f_6 * ksi0_831[k]
                    - f_7 * ksi1_831[k]
                    + f_3 * pc_x[k] * ksk_1063[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pa_z, pc_x, pc_z, isl0_966, isl1_966, \
                         ksi0_832, ksi0_834, ksi1_832, ksi1_834, ksk_1064, \
                         ksk_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_6 * ksi0_832[k]
                    - f_7 * ksi1_832[k]
                    + f_3 * pc_x[k] * ksk_1064[k];

        t_1326[k] = pa_z[k] * isl0_966[k]
                    - f_14 * pc_z[k] * isl1_966[k];

        t_1327[k] = f_4 * ksi0_834[k]
                    - f_5 * ksi1_834[k]
                    + f_3 * pc_x[k] * ksk_1066[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pc_x, ksi0_835, ksi0_836, ksi0_837, ksi1_835, \
                         ksi1_836, ksi1_837, ksk_1067, ksk_1068, \
                         ksk_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_4 * ksi0_835[k]
                    - f_5 * ksi1_835[k]
                    + f_3 * pc_x[k] * ksk_1067[k];

        t_1329[k] = f_4 * ksi0_836[k]
                    - f_5 * ksi1_836[k]
                    + f_3 * pc_x[k] * ksk_1068[k];

        t_1330[k] = f_4 * ksi0_837[k]
                    - f_5 * ksi1_837[k]
                    + f_3 * pc_x[k] * ksk_1069[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, t_1335, pc_x, ksi0_838, ksi0_839, \
                         ksi1_838, ksi1_839, ksk_1070, ksk_1071, ksk_1072, ksk_1073, \
                         ksk_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_4 * ksi0_838[k]
                    - f_5 * ksi1_838[k]
                    + f_3 * pc_x[k] * ksk_1070[k];

        t_1332[k] = f_4 * ksi0_839[k]
                    - f_5 * ksi1_839[k]
                    + f_3 * pc_x[k] * ksk_1071[k];

        t_1333[k] = f_3 * pc_x[k] * ksk_1072[k];

        t_1334[k] = f_3 * pc_x[k] * ksk_1073[k];

        t_1335[k] = f_3 * pc_x[k] * ksk_1074[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, t_1339, t_1340, t_1341, pa_z, pc_x, pc_z, \
                         isl0_981, isl1_981, ksk_1075, ksk_1076, ksk_1077, ksk_1078, \
                         ksk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_3 * pc_x[k] * ksk_1075[k];

        t_1337[k] = f_3 * pc_x[k] * ksk_1076[k];

        t_1338[k] = f_3 * pc_x[k] * ksk_1077[k];

        t_1339[k] = f_3 * pc_x[k] * ksk_1078[k];

        t_1340[k] = f_3 * pc_x[k] * ksk_1079[k];

        t_1341[k] = pa_z[k] * isl0_981[k]
                    - f_14 * pc_z[k] * isl1_981[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, pa_z, pc_z, isl0_983, isl0_984, isk_784, \
                         isk_785, isk_786, isl1_983, isl1_984, \
                         ksk_1072 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_15 * isk_784[k]
                    + f_3 * pc_z[k] * ksk_1072[k];

        t_1343[k] = pa_z[k] * isl0_983[k]
                    + f_16 * isk_785[k]
                    - f_14 * pc_z[k] * isl1_983[k];

        t_1344[k] = pa_z[k] * isl0_984[k]
                    + f_17 * isk_786[k]
                    - f_14 * pc_z[k] * isl1_984[k];
    }

#pragma omp simd aligned(t_1345, t_1346, t_1347, pa_z, pc_z, isl0_985, isl0_986, isl0_987, \
                         isk_787, isk_788, isk_789, isl1_985, isl1_986, \
                         isl1_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1345[k] = pa_z[k] * isl0_985[k]
                    + f_18 * isk_787[k]
                    - f_14 * pc_z[k] * isl1_985[k];

        t_1346[k] = pa_z[k] * isl0_986[k]
                    + f_19 * isk_788[k]
                    - f_14 * pc_z[k] * isl1_986[k];

        t_1347[k] = pa_z[k] * isl0_987[k]
                    + f_20 * isk_789[k]
                    - f_14 * pc_z[k] * isl1_987[k];
    }

#pragma omp simd aligned(t_1348, t_1349, t_1350, pc_x, pc_y, pc_z, isk_791, isk_827, ksi0_839, \
                         ksi0_840, ksi1_839, ksi1_840, ksk_1079, \
                         ksk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1348[k] = f_20 * isk_827[k]
                    + f_3 * pc_y[k] * ksk_1079[k];

        t_1349[k] = f_15 * isk_791[k]
                    + f_1 * ksi0_839[k]
                    - f_2 * ksi1_839[k]
                    + f_3 * pc_z[k] * ksk_1079[k];

        t_1350[k] = f_1 * ksi0_840[k]
                    - f_2 * ksi1_840[k]
                    + f_3 * pc_x[k] * ksk_1080[k];
    }

#pragma omp simd aligned(t_1351, t_1352, t_1353, pc_x, ksi0_841, ksi0_842, ksi0_843, ksi1_841, \
                         ksi1_842, ksi1_843, ksk_1081, ksk_1082, \
                         ksk_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1351[k] = f_21 * ksi0_841[k]
                    - f_22 * ksi1_841[k]
                    + f_3 * pc_x[k] * ksk_1081[k];

        t_1352[k] = f_21 * ksi0_842[k]
                    - f_22 * ksi1_842[k]
                    + f_3 * pc_x[k] * ksk_1082[k];

        t_1353[k] = f_12 * ksi0_843[k]
                    - f_13 * ksi1_843[k]
                    + f_3 * pc_x[k] * ksk_1083[k];
    }

#pragma omp simd aligned(t_1354, t_1355, t_1356, pc_x, ksi0_844, ksi0_845, ksi0_846, ksi1_844, \
                         ksi1_845, ksi1_846, ksk_1084, ksk_1085, \
                         ksk_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1354[k] = f_12 * ksi0_844[k]
                    - f_13 * ksi1_844[k]
                    + f_3 * pc_x[k] * ksk_1084[k];

        t_1355[k] = f_12 * ksi0_845[k]
                    - f_13 * ksi1_845[k]
                    + f_3 * pc_x[k] * ksk_1085[k];

        t_1356[k] = f_10 * ksi0_846[k]
                    - f_11 * ksi1_846[k]
                    + f_3 * pc_x[k] * ksk_1086[k];
    }

#pragma omp simd aligned(t_1357, t_1358, t_1359, pc_x, ksi0_847, ksi0_848, ksi0_849, ksi1_847, \
                         ksi1_848, ksi1_849, ksk_1087, ksk_1088, \
                         ksk_1089 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1357[k] = f_10 * ksi0_847[k]
                    - f_11 * ksi1_847[k]
                    + f_3 * pc_x[k] * ksk_1087[k];

        t_1358[k] = f_10 * ksi0_848[k]
                    - f_11 * ksi1_848[k]
                    + f_3 * pc_x[k] * ksk_1088[k];

        t_1359[k] = f_10 * ksi0_849[k]
                    - f_11 * ksi1_849[k]
                    + f_3 * pc_x[k] * ksk_1089[k];
    }

#pragma omp simd aligned(t_1360, t_1361, t_1362, pc_x, ksi0_850, ksi0_851, ksi0_852, ksi1_850, \
                         ksi1_851, ksi1_852, ksk_1090, ksk_1091, \
                         ksk_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1360[k] = f_8 * ksi0_850[k]
                    - f_9 * ksi1_850[k]
                    + f_3 * pc_x[k] * ksk_1090[k];

        t_1361[k] = f_8 * ksi0_851[k]
                    - f_9 * ksi1_851[k]
                    + f_3 * pc_x[k] * ksk_1091[k];

        t_1362[k] = f_8 * ksi0_852[k]
                    - f_9 * ksi1_852[k]
                    + f_3 * pc_x[k] * ksk_1092[k];
    }

#pragma omp simd aligned(t_1363, t_1364, t_1365, pc_x, ksi0_853, ksi0_854, ksi0_855, ksi1_853, \
                         ksi1_854, ksi1_855, ksk_1093, ksk_1094, \
                         ksk_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1363[k] = f_8 * ksi0_853[k]
                    - f_9 * ksi1_853[k]
                    + f_3 * pc_x[k] * ksk_1093[k];

        t_1364[k] = f_8 * ksi0_854[k]
                    - f_9 * ksi1_854[k]
                    + f_3 * pc_x[k] * ksk_1094[k];

        t_1365[k] = f_6 * ksi0_855[k]
                    - f_7 * ksi1_855[k]
                    + f_3 * pc_x[k] * ksk_1095[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, pc_x, ksi0_856, ksi0_857, ksi0_858, ksi1_856, \
                         ksi1_857, ksi1_858, ksk_1096, ksk_1097, \
                         ksk_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_6 * ksi0_856[k]
                    - f_7 * ksi1_856[k]
                    + f_3 * pc_x[k] * ksk_1096[k];

        t_1367[k] = f_6 * ksi0_857[k]
                    - f_7 * ksi1_857[k]
                    + f_3 * pc_x[k] * ksk_1097[k];

        t_1368[k] = f_6 * ksi0_858[k]
                    - f_7 * ksi1_858[k]
                    + f_3 * pc_x[k] * ksk_1098[k];
    }

#pragma omp simd aligned(t_1369, t_1370, t_1371, pc_x, ksi0_859, ksi0_860, ksi0_861, ksi1_859, \
                         ksi1_860, ksi1_861, ksk_1099, ksk_1100, \
                         ksk_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1369[k] = f_6 * ksi0_859[k]
                    - f_7 * ksi1_859[k]
                    + f_3 * pc_x[k] * ksk_1099[k];

        t_1370[k] = f_6 * ksi0_860[k]
                    - f_7 * ksi1_860[k]
                    + f_3 * pc_x[k] * ksk_1100[k];

        t_1371[k] = f_4 * ksi0_861[k]
                    - f_5 * ksi1_861[k]
                    + f_3 * pc_x[k] * ksk_1101[k];
    }

#pragma omp simd aligned(t_1372, t_1373, t_1374, pc_x, ksi0_862, ksi0_863, ksi0_864, ksi1_862, \
                         ksi1_863, ksi1_864, ksk_1102, ksk_1103, \
                         ksk_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1372[k] = f_4 * ksi0_862[k]
                    - f_5 * ksi1_862[k]
                    + f_3 * pc_x[k] * ksk_1102[k];

        t_1373[k] = f_4 * ksi0_863[k]
                    - f_5 * ksi1_863[k]
                    + f_3 * pc_x[k] * ksk_1103[k];

        t_1374[k] = f_4 * ksi0_864[k]
                    - f_5 * ksi1_864[k]
                    + f_3 * pc_x[k] * ksk_1104[k];
    }

#pragma omp simd aligned(t_1375, t_1376, t_1377, t_1378, pc_x, ksi0_865, ksi0_866, ksi0_867, \
                         ksi1_865, ksi1_866, ksi1_867, ksk_1105, ksk_1106, ksk_1107, \
                         ksk_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1375[k] = f_4 * ksi0_865[k]
                    - f_5 * ksi1_865[k]
                    + f_3 * pc_x[k] * ksk_1105[k];

        t_1376[k] = f_4 * ksi0_866[k]
                    - f_5 * ksi1_866[k]
                    + f_3 * pc_x[k] * ksk_1106[k];

        t_1377[k] = f_4 * ksi0_867[k]
                    - f_5 * ksi1_867[k]
                    + f_3 * pc_x[k] * ksk_1107[k];

        t_1378[k] = f_3 * pc_x[k] * ksk_1108[k];
    }

#pragma omp simd aligned(t_1379, t_1380, t_1381, t_1382, t_1383, t_1384, t_1385, pc_x, \
                         ksk_1109, ksk_1110, ksk_1111, ksk_1112, ksk_1113, ksk_1114, \
                         ksk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1379[k] = f_3 * pc_x[k] * ksk_1109[k];

        t_1380[k] = f_3 * pc_x[k] * ksk_1110[k];

        t_1381[k] = f_3 * pc_x[k] * ksk_1111[k];

        t_1382[k] = f_3 * pc_x[k] * ksk_1112[k];

        t_1383[k] = f_3 * pc_x[k] * ksk_1113[k];

        t_1384[k] = f_3 * pc_x[k] * ksk_1114[k];

        t_1385[k] = f_3 * pc_x[k] * ksk_1115[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, pc_y, pc_z, isk_820, isk_856, isk_858, \
                         ksi0_861, ksi0_863, ksi1_861, ksi1_863, ksk_1108, \
                         ksk_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = f_19 * isk_856[k]
                    + f_1 * ksi0_861[k]
                    - f_2 * ksi1_861[k]
                    + f_3 * pc_y[k] * ksk_1108[k];

        t_1387[k] = f_16 * isk_820[k]
                    + f_3 * pc_z[k] * ksk_1108[k];

        t_1388[k] = f_19 * isk_858[k]
                    + f_12 * ksi0_863[k]
                    - f_13 * ksi1_863[k]
                    + f_3 * pc_y[k] * ksk_1110[k];
    }

#pragma omp simd aligned(t_1389, t_1390, t_1391, pc_y, isk_859, isk_860, isk_861, ksi0_864, \
                         ksi0_865, ksi0_866, ksi1_864, ksi1_865, ksi1_866, ksk_1111, ksk_1112, \
                         ksk_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1389[k] = f_19 * isk_859[k]
                    + f_10 * ksi0_864[k]
                    - f_11 * ksi1_864[k]
                    + f_3 * pc_y[k] * ksk_1111[k];

        t_1390[k] = f_19 * isk_860[k]
                    + f_8 * ksi0_865[k]
                    - f_9 * ksi1_865[k]
                    + f_3 * pc_y[k] * ksk_1112[k];

        t_1391[k] = f_19 * isk_861[k]
                    + f_6 * ksi0_866[k]
                    - f_7 * ksi1_866[k]
                    + f_3 * pc_y[k] * ksk_1113[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, pc_y, pc_z, isk_827, isk_862, isk_863, \
                         ksi0_867, ksi1_867, ksk_1114, ksk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_19 * isk_862[k]
                    + f_4 * ksi0_867[k]
                    - f_5 * ksi1_867[k]
                    + f_3 * pc_y[k] * ksk_1114[k];

        t_1393[k] = f_19 * isk_863[k]
                    + f_3 * pc_y[k] * ksk_1115[k];

        t_1394[k] = f_16 * isk_827[k]
                    + f_1 * ksi0_867[k]
                    - f_2 * ksi1_867[k]
                    + f_3 * pc_z[k] * ksk_1115[k];
    }

#pragma omp simd aligned(t_1395, t_1396, t_1397, pc_x, ksi0_868, ksi0_869, ksi0_870, ksi1_868, \
                         ksi1_869, ksi1_870, ksk_1116, ksk_1117, \
                         ksk_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1395[k] = f_1 * ksi0_868[k]
                    - f_2 * ksi1_868[k]
                    + f_3 * pc_x[k] * ksk_1116[k];

        t_1396[k] = f_21 * ksi0_869[k]
                    - f_22 * ksi1_869[k]
                    + f_3 * pc_x[k] * ksk_1117[k];

        t_1397[k] = f_21 * ksi0_870[k]
                    - f_22 * ksi1_870[k]
                    + f_3 * pc_x[k] * ksk_1118[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t isk, const size_t ksi0,
                                                           const size_t ksi1, const size_t ksk,
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
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isk_856 = buffer.data(isk + 856);
    const auto *isk_863 = buffer.data(isk + 863);
    const auto *isk_892 = buffer.data(isk + 892);
    const auto *isk_894 = buffer.data(isk + 894);
    const auto *isk_895 = buffer.data(isk + 895);
    const auto *isk_896 = buffer.data(isk + 896);
    const auto *isk_897 = buffer.data(isk + 897);
    const auto *isk_898 = buffer.data(isk + 898);
    const auto *isk_899 = buffer.data(isk + 899);
    const auto *isk_928 = buffer.data(isk + 928);
    const auto *isk_930 = buffer.data(isk + 930);
    const auto *isk_931 = buffer.data(isk + 931);
    const auto *isk_932 = buffer.data(isk + 932);
    const auto *isk_933 = buffer.data(isk + 933);
    const auto *isk_934 = buffer.data(isk + 934);
    const auto *isk_935 = buffer.data(isk + 935);

    const auto *ksi0_871 = buffer.data(ksi0 + 871);
    const auto *ksi0_872 = buffer.data(ksi0 + 872);
    const auto *ksi0_873 = buffer.data(ksi0 + 873);
    const auto *ksi0_874 = buffer.data(ksi0 + 874);
    const auto *ksi0_875 = buffer.data(ksi0 + 875);
    const auto *ksi0_876 = buffer.data(ksi0 + 876);
    const auto *ksi0_877 = buffer.data(ksi0 + 877);
    const auto *ksi0_878 = buffer.data(ksi0 + 878);
    const auto *ksi0_879 = buffer.data(ksi0 + 879);
    const auto *ksi0_880 = buffer.data(ksi0 + 880);
    const auto *ksi0_881 = buffer.data(ksi0 + 881);
    const auto *ksi0_882 = buffer.data(ksi0 + 882);
    const auto *ksi0_883 = buffer.data(ksi0 + 883);
    const auto *ksi0_884 = buffer.data(ksi0 + 884);
    const auto *ksi0_885 = buffer.data(ksi0 + 885);
    const auto *ksi0_886 = buffer.data(ksi0 + 886);
    const auto *ksi0_887 = buffer.data(ksi0 + 887);
    const auto *ksi0_888 = buffer.data(ksi0 + 888);
    const auto *ksi0_889 = buffer.data(ksi0 + 889);
    const auto *ksi0_890 = buffer.data(ksi0 + 890);
    const auto *ksi0_891 = buffer.data(ksi0 + 891);
    const auto *ksi0_892 = buffer.data(ksi0 + 892);
    const auto *ksi0_893 = buffer.data(ksi0 + 893);
    const auto *ksi0_894 = buffer.data(ksi0 + 894);
    const auto *ksi0_895 = buffer.data(ksi0 + 895);
    const auto *ksi0_896 = buffer.data(ksi0 + 896);
    const auto *ksi0_897 = buffer.data(ksi0 + 897);
    const auto *ksi0_898 = buffer.data(ksi0 + 898);
    const auto *ksi0_899 = buffer.data(ksi0 + 899);
    const auto *ksi0_900 = buffer.data(ksi0 + 900);
    const auto *ksi0_901 = buffer.data(ksi0 + 901);
    const auto *ksi0_902 = buffer.data(ksi0 + 902);
    const auto *ksi0_903 = buffer.data(ksi0 + 903);
    const auto *ksi0_904 = buffer.data(ksi0 + 904);
    const auto *ksi0_905 = buffer.data(ksi0 + 905);
    const auto *ksi0_906 = buffer.data(ksi0 + 906);
    const auto *ksi0_907 = buffer.data(ksi0 + 907);
    const auto *ksi0_908 = buffer.data(ksi0 + 908);
    const auto *ksi0_909 = buffer.data(ksi0 + 909);
    const auto *ksi0_910 = buffer.data(ksi0 + 910);
    const auto *ksi0_911 = buffer.data(ksi0 + 911);
    const auto *ksi0_912 = buffer.data(ksi0 + 912);
    const auto *ksi0_913 = buffer.data(ksi0 + 913);
    const auto *ksi0_914 = buffer.data(ksi0 + 914);
    const auto *ksi0_915 = buffer.data(ksi0 + 915);
    const auto *ksi0_916 = buffer.data(ksi0 + 916);
    const auto *ksi0_917 = buffer.data(ksi0 + 917);
    const auto *ksi0_918 = buffer.data(ksi0 + 918);
    const auto *ksi0_919 = buffer.data(ksi0 + 919);
    const auto *ksi0_920 = buffer.data(ksi0 + 920);
    const auto *ksi0_921 = buffer.data(ksi0 + 921);
    const auto *ksi0_922 = buffer.data(ksi0 + 922);
    const auto *ksi0_923 = buffer.data(ksi0 + 923);
    const auto *ksi0_924 = buffer.data(ksi0 + 924);
    const auto *ksi0_925 = buffer.data(ksi0 + 925);
    const auto *ksi0_926 = buffer.data(ksi0 + 926);
    const auto *ksi0_927 = buffer.data(ksi0 + 927);
    const auto *ksi0_928 = buffer.data(ksi0 + 928);
    const auto *ksi0_929 = buffer.data(ksi0 + 929);
    const auto *ksi0_930 = buffer.data(ksi0 + 930);
    const auto *ksi0_931 = buffer.data(ksi0 + 931);
    const auto *ksi0_932 = buffer.data(ksi0 + 932);
    const auto *ksi0_933 = buffer.data(ksi0 + 933);
    const auto *ksi0_934 = buffer.data(ksi0 + 934);
    const auto *ksi0_935 = buffer.data(ksi0 + 935);
    const auto *ksi0_936 = buffer.data(ksi0 + 936);
    const auto *ksi0_937 = buffer.data(ksi0 + 937);
    const auto *ksi0_938 = buffer.data(ksi0 + 938);
    const auto *ksi0_939 = buffer.data(ksi0 + 939);
    const auto *ksi0_940 = buffer.data(ksi0 + 940);
    const auto *ksi0_941 = buffer.data(ksi0 + 941);
    const auto *ksi0_942 = buffer.data(ksi0 + 942);
    const auto *ksi0_943 = buffer.data(ksi0 + 943);
    const auto *ksi0_944 = buffer.data(ksi0 + 944);
    const auto *ksi0_945 = buffer.data(ksi0 + 945);
    const auto *ksi0_946 = buffer.data(ksi0 + 946);
    const auto *ksi0_947 = buffer.data(ksi0 + 947);

    const auto *ksi1_871 = buffer.data(ksi1 + 871);
    const auto *ksi1_872 = buffer.data(ksi1 + 872);
    const auto *ksi1_873 = buffer.data(ksi1 + 873);
    const auto *ksi1_874 = buffer.data(ksi1 + 874);
    const auto *ksi1_875 = buffer.data(ksi1 + 875);
    const auto *ksi1_876 = buffer.data(ksi1 + 876);
    const auto *ksi1_877 = buffer.data(ksi1 + 877);
    const auto *ksi1_878 = buffer.data(ksi1 + 878);
    const auto *ksi1_879 = buffer.data(ksi1 + 879);
    const auto *ksi1_880 = buffer.data(ksi1 + 880);
    const auto *ksi1_881 = buffer.data(ksi1 + 881);
    const auto *ksi1_882 = buffer.data(ksi1 + 882);
    const auto *ksi1_883 = buffer.data(ksi1 + 883);
    const auto *ksi1_884 = buffer.data(ksi1 + 884);
    const auto *ksi1_885 = buffer.data(ksi1 + 885);
    const auto *ksi1_886 = buffer.data(ksi1 + 886);
    const auto *ksi1_887 = buffer.data(ksi1 + 887);
    const auto *ksi1_888 = buffer.data(ksi1 + 888);
    const auto *ksi1_889 = buffer.data(ksi1 + 889);
    const auto *ksi1_890 = buffer.data(ksi1 + 890);
    const auto *ksi1_891 = buffer.data(ksi1 + 891);
    const auto *ksi1_892 = buffer.data(ksi1 + 892);
    const auto *ksi1_893 = buffer.data(ksi1 + 893);
    const auto *ksi1_894 = buffer.data(ksi1 + 894);
    const auto *ksi1_895 = buffer.data(ksi1 + 895);
    const auto *ksi1_896 = buffer.data(ksi1 + 896);
    const auto *ksi1_897 = buffer.data(ksi1 + 897);
    const auto *ksi1_898 = buffer.data(ksi1 + 898);
    const auto *ksi1_899 = buffer.data(ksi1 + 899);
    const auto *ksi1_900 = buffer.data(ksi1 + 900);
    const auto *ksi1_901 = buffer.data(ksi1 + 901);
    const auto *ksi1_902 = buffer.data(ksi1 + 902);
    const auto *ksi1_903 = buffer.data(ksi1 + 903);
    const auto *ksi1_904 = buffer.data(ksi1 + 904);
    const auto *ksi1_905 = buffer.data(ksi1 + 905);
    const auto *ksi1_906 = buffer.data(ksi1 + 906);
    const auto *ksi1_907 = buffer.data(ksi1 + 907);
    const auto *ksi1_908 = buffer.data(ksi1 + 908);
    const auto *ksi1_909 = buffer.data(ksi1 + 909);
    const auto *ksi1_910 = buffer.data(ksi1 + 910);
    const auto *ksi1_911 = buffer.data(ksi1 + 911);
    const auto *ksi1_912 = buffer.data(ksi1 + 912);
    const auto *ksi1_913 = buffer.data(ksi1 + 913);
    const auto *ksi1_914 = buffer.data(ksi1 + 914);
    const auto *ksi1_915 = buffer.data(ksi1 + 915);
    const auto *ksi1_916 = buffer.data(ksi1 + 916);
    const auto *ksi1_917 = buffer.data(ksi1 + 917);
    const auto *ksi1_918 = buffer.data(ksi1 + 918);
    const auto *ksi1_919 = buffer.data(ksi1 + 919);
    const auto *ksi1_920 = buffer.data(ksi1 + 920);
    const auto *ksi1_921 = buffer.data(ksi1 + 921);
    const auto *ksi1_922 = buffer.data(ksi1 + 922);
    const auto *ksi1_923 = buffer.data(ksi1 + 923);
    const auto *ksi1_924 = buffer.data(ksi1 + 924);
    const auto *ksi1_925 = buffer.data(ksi1 + 925);
    const auto *ksi1_926 = buffer.data(ksi1 + 926);
    const auto *ksi1_927 = buffer.data(ksi1 + 927);
    const auto *ksi1_928 = buffer.data(ksi1 + 928);
    const auto *ksi1_929 = buffer.data(ksi1 + 929);
    const auto *ksi1_930 = buffer.data(ksi1 + 930);
    const auto *ksi1_931 = buffer.data(ksi1 + 931);
    const auto *ksi1_932 = buffer.data(ksi1 + 932);
    const auto *ksi1_933 = buffer.data(ksi1 + 933);
    const auto *ksi1_934 = buffer.data(ksi1 + 934);
    const auto *ksi1_935 = buffer.data(ksi1 + 935);
    const auto *ksi1_936 = buffer.data(ksi1 + 936);
    const auto *ksi1_937 = buffer.data(ksi1 + 937);
    const auto *ksi1_938 = buffer.data(ksi1 + 938);
    const auto *ksi1_939 = buffer.data(ksi1 + 939);
    const auto *ksi1_940 = buffer.data(ksi1 + 940);
    const auto *ksi1_941 = buffer.data(ksi1 + 941);
    const auto *ksi1_942 = buffer.data(ksi1 + 942);
    const auto *ksi1_943 = buffer.data(ksi1 + 943);
    const auto *ksi1_944 = buffer.data(ksi1 + 944);
    const auto *ksi1_945 = buffer.data(ksi1 + 945);
    const auto *ksi1_946 = buffer.data(ksi1 + 946);
    const auto *ksi1_947 = buffer.data(ksi1 + 947);

    const auto *ksk_1119 = buffer.data(ksk + 1119);
    const auto *ksk_1120 = buffer.data(ksk + 1120);
    const auto *ksk_1121 = buffer.data(ksk + 1121);
    const auto *ksk_1122 = buffer.data(ksk + 1122);
    const auto *ksk_1123 = buffer.data(ksk + 1123);
    const auto *ksk_1124 = buffer.data(ksk + 1124);
    const auto *ksk_1125 = buffer.data(ksk + 1125);
    const auto *ksk_1126 = buffer.data(ksk + 1126);
    const auto *ksk_1127 = buffer.data(ksk + 1127);
    const auto *ksk_1128 = buffer.data(ksk + 1128);
    const auto *ksk_1129 = buffer.data(ksk + 1129);
    const auto *ksk_1130 = buffer.data(ksk + 1130);
    const auto *ksk_1131 = buffer.data(ksk + 1131);
    const auto *ksk_1132 = buffer.data(ksk + 1132);
    const auto *ksk_1133 = buffer.data(ksk + 1133);
    const auto *ksk_1134 = buffer.data(ksk + 1134);
    const auto *ksk_1135 = buffer.data(ksk + 1135);
    const auto *ksk_1136 = buffer.data(ksk + 1136);
    const auto *ksk_1137 = buffer.data(ksk + 1137);
    const auto *ksk_1138 = buffer.data(ksk + 1138);
    const auto *ksk_1139 = buffer.data(ksk + 1139);
    const auto *ksk_1140 = buffer.data(ksk + 1140);
    const auto *ksk_1141 = buffer.data(ksk + 1141);
    const auto *ksk_1142 = buffer.data(ksk + 1142);
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
    const auto *ksk_1153 = buffer.data(ksk + 1153);
    const auto *ksk_1154 = buffer.data(ksk + 1154);
    const auto *ksk_1155 = buffer.data(ksk + 1155);
    const auto *ksk_1156 = buffer.data(ksk + 1156);
    const auto *ksk_1157 = buffer.data(ksk + 1157);
    const auto *ksk_1158 = buffer.data(ksk + 1158);
    const auto *ksk_1159 = buffer.data(ksk + 1159);
    const auto *ksk_1160 = buffer.data(ksk + 1160);
    const auto *ksk_1161 = buffer.data(ksk + 1161);
    const auto *ksk_1162 = buffer.data(ksk + 1162);
    const auto *ksk_1163 = buffer.data(ksk + 1163);
    const auto *ksk_1164 = buffer.data(ksk + 1164);
    const auto *ksk_1165 = buffer.data(ksk + 1165);
    const auto *ksk_1166 = buffer.data(ksk + 1166);
    const auto *ksk_1167 = buffer.data(ksk + 1167);
    const auto *ksk_1168 = buffer.data(ksk + 1168);
    const auto *ksk_1169 = buffer.data(ksk + 1169);
    const auto *ksk_1170 = buffer.data(ksk + 1170);
    const auto *ksk_1171 = buffer.data(ksk + 1171);
    const auto *ksk_1172 = buffer.data(ksk + 1172);
    const auto *ksk_1173 = buffer.data(ksk + 1173);
    const auto *ksk_1174 = buffer.data(ksk + 1174);
    const auto *ksk_1175 = buffer.data(ksk + 1175);
    const auto *ksk_1176 = buffer.data(ksk + 1176);
    const auto *ksk_1177 = buffer.data(ksk + 1177);
    const auto *ksk_1178 = buffer.data(ksk + 1178);
    const auto *ksk_1179 = buffer.data(ksk + 1179);
    const auto *ksk_1180 = buffer.data(ksk + 1180);
    const auto *ksk_1181 = buffer.data(ksk + 1181);
    const auto *ksk_1182 = buffer.data(ksk + 1182);
    const auto *ksk_1183 = buffer.data(ksk + 1183);
    const auto *ksk_1184 = buffer.data(ksk + 1184);
    const auto *ksk_1185 = buffer.data(ksk + 1185);
    const auto *ksk_1186 = buffer.data(ksk + 1186);
    const auto *ksk_1187 = buffer.data(ksk + 1187);
    const auto *ksk_1188 = buffer.data(ksk + 1188);
    const auto *ksk_1189 = buffer.data(ksk + 1189);
    const auto *ksk_1190 = buffer.data(ksk + 1190);
    const auto *ksk_1191 = buffer.data(ksk + 1191);
    const auto *ksk_1192 = buffer.data(ksk + 1192);
    const auto *ksk_1193 = buffer.data(ksk + 1193);
    const auto *ksk_1194 = buffer.data(ksk + 1194);
    const auto *ksk_1195 = buffer.data(ksk + 1195);
    const auto *ksk_1196 = buffer.data(ksk + 1196);
    const auto *ksk_1197 = buffer.data(ksk + 1197);
    const auto *ksk_1198 = buffer.data(ksk + 1198);
    const auto *ksk_1199 = buffer.data(ksk + 1199);
    const auto *ksk_1200 = buffer.data(ksk + 1200);
    const auto *ksk_1201 = buffer.data(ksk + 1201);
    const auto *ksk_1202 = buffer.data(ksk + 1202);
    const auto *ksk_1203 = buffer.data(ksk + 1203);
    const auto *ksk_1204 = buffer.data(ksk + 1204);
    const auto *ksk_1205 = buffer.data(ksk + 1205);
    const auto *ksk_1206 = buffer.data(ksk + 1206);
    const auto *ksk_1207 = buffer.data(ksk + 1207);
    const auto *ksk_1208 = buffer.data(ksk + 1208);
    const auto *ksk_1209 = buffer.data(ksk + 1209);
    const auto *ksk_1210 = buffer.data(ksk + 1210);
    const auto *ksk_1211 = buffer.data(ksk + 1211);

#pragma omp simd aligned(t_1398, t_1399, t_1400, pc_x, ksi0_871, ksi0_872, ksi0_873, ksi1_871, \
                         ksi1_872, ksi1_873, ksk_1119, ksk_1120, \
                         ksk_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_12 * ksi0_871[k]
                    - f_13 * ksi1_871[k]
                    + f_3 * pc_x[k] * ksk_1119[k];

        t_1399[k] = f_12 * ksi0_872[k]
                    - f_13 * ksi1_872[k]
                    + f_3 * pc_x[k] * ksk_1120[k];

        t_1400[k] = f_12 * ksi0_873[k]
                    - f_13 * ksi1_873[k]
                    + f_3 * pc_x[k] * ksk_1121[k];
    }

#pragma omp simd aligned(t_1401, t_1402, t_1403, pc_x, ksi0_874, ksi0_875, ksi0_876, ksi1_874, \
                         ksi1_875, ksi1_876, ksk_1122, ksk_1123, \
                         ksk_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1401[k] = f_10 * ksi0_874[k]
                    - f_11 * ksi1_874[k]
                    + f_3 * pc_x[k] * ksk_1122[k];

        t_1402[k] = f_10 * ksi0_875[k]
                    - f_11 * ksi1_875[k]
                    + f_3 * pc_x[k] * ksk_1123[k];

        t_1403[k] = f_10 * ksi0_876[k]
                    - f_11 * ksi1_876[k]
                    + f_3 * pc_x[k] * ksk_1124[k];
    }

#pragma omp simd aligned(t_1404, t_1405, t_1406, pc_x, ksi0_877, ksi0_878, ksi0_879, ksi1_877, \
                         ksi1_878, ksi1_879, ksk_1125, ksk_1126, \
                         ksk_1127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1404[k] = f_10 * ksi0_877[k]
                    - f_11 * ksi1_877[k]
                    + f_3 * pc_x[k] * ksk_1125[k];

        t_1405[k] = f_8 * ksi0_878[k]
                    - f_9 * ksi1_878[k]
                    + f_3 * pc_x[k] * ksk_1126[k];

        t_1406[k] = f_8 * ksi0_879[k]
                    - f_9 * ksi1_879[k]
                    + f_3 * pc_x[k] * ksk_1127[k];
    }

#pragma omp simd aligned(t_1407, t_1408, t_1409, pc_x, ksi0_880, ksi0_881, ksi0_882, ksi1_880, \
                         ksi1_881, ksi1_882, ksk_1128, ksk_1129, \
                         ksk_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1407[k] = f_8 * ksi0_880[k]
                    - f_9 * ksi1_880[k]
                    + f_3 * pc_x[k] * ksk_1128[k];

        t_1408[k] = f_8 * ksi0_881[k]
                    - f_9 * ksi1_881[k]
                    + f_3 * pc_x[k] * ksk_1129[k];

        t_1409[k] = f_8 * ksi0_882[k]
                    - f_9 * ksi1_882[k]
                    + f_3 * pc_x[k] * ksk_1130[k];
    }

#pragma omp simd aligned(t_1410, t_1411, t_1412, pc_x, ksi0_883, ksi0_884, ksi0_885, ksi1_883, \
                         ksi1_884, ksi1_885, ksk_1131, ksk_1132, \
                         ksk_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1410[k] = f_6 * ksi0_883[k]
                    - f_7 * ksi1_883[k]
                    + f_3 * pc_x[k] * ksk_1131[k];

        t_1411[k] = f_6 * ksi0_884[k]
                    - f_7 * ksi1_884[k]
                    + f_3 * pc_x[k] * ksk_1132[k];

        t_1412[k] = f_6 * ksi0_885[k]
                    - f_7 * ksi1_885[k]
                    + f_3 * pc_x[k] * ksk_1133[k];
    }

#pragma omp simd aligned(t_1413, t_1414, t_1415, pc_x, ksi0_886, ksi0_887, ksi0_888, ksi1_886, \
                         ksi1_887, ksi1_888, ksk_1134, ksk_1135, \
                         ksk_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1413[k] = f_6 * ksi0_886[k]
                    - f_7 * ksi1_886[k]
                    + f_3 * pc_x[k] * ksk_1134[k];

        t_1414[k] = f_6 * ksi0_887[k]
                    - f_7 * ksi1_887[k]
                    + f_3 * pc_x[k] * ksk_1135[k];

        t_1415[k] = f_6 * ksi0_888[k]
                    - f_7 * ksi1_888[k]
                    + f_3 * pc_x[k] * ksk_1136[k];
    }

#pragma omp simd aligned(t_1416, t_1417, t_1418, pc_x, ksi0_889, ksi0_890, ksi0_891, ksi1_889, \
                         ksi1_890, ksi1_891, ksk_1137, ksk_1138, \
                         ksk_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1416[k] = f_4 * ksi0_889[k]
                    - f_5 * ksi1_889[k]
                    + f_3 * pc_x[k] * ksk_1137[k];

        t_1417[k] = f_4 * ksi0_890[k]
                    - f_5 * ksi1_890[k]
                    + f_3 * pc_x[k] * ksk_1138[k];

        t_1418[k] = f_4 * ksi0_891[k]
                    - f_5 * ksi1_891[k]
                    + f_3 * pc_x[k] * ksk_1139[k];
    }

#pragma omp simd aligned(t_1419, t_1420, t_1421, pc_x, ksi0_892, ksi0_893, ksi0_894, ksi1_892, \
                         ksi1_893, ksi1_894, ksk_1140, ksk_1141, \
                         ksk_1142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1419[k] = f_4 * ksi0_892[k]
                    - f_5 * ksi1_892[k]
                    + f_3 * pc_x[k] * ksk_1140[k];

        t_1420[k] = f_4 * ksi0_893[k]
                    - f_5 * ksi1_893[k]
                    + f_3 * pc_x[k] * ksk_1141[k];

        t_1421[k] = f_4 * ksi0_894[k]
                    - f_5 * ksi1_894[k]
                    + f_3 * pc_x[k] * ksk_1142[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, t_1425, t_1426, t_1427, pc_x, ksi0_895, \
                         ksi1_895, ksk_1143, ksk_1144, ksk_1145, ksk_1146, ksk_1147, \
                         ksk_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_4 * ksi0_895[k]
                    - f_5 * ksi1_895[k]
                    + f_3 * pc_x[k] * ksk_1143[k];

        t_1423[k] = f_3 * pc_x[k] * ksk_1144[k];

        t_1424[k] = f_3 * pc_x[k] * ksk_1145[k];

        t_1425[k] = f_3 * pc_x[k] * ksk_1146[k];

        t_1426[k] = f_3 * pc_x[k] * ksk_1147[k];

        t_1427[k] = f_3 * pc_x[k] * ksk_1148[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, t_1431, t_1432, pc_x, pc_y, pc_z, isk_856, \
                         isk_892, ksi0_889, ksi1_889, ksk_1144, ksk_1149, ksk_1150, \
                         ksk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_3 * pc_x[k] * ksk_1149[k];

        t_1429[k] = f_3 * pc_x[k] * ksk_1150[k];

        t_1430[k] = f_3 * pc_x[k] * ksk_1151[k];

        t_1431[k] = f_18 * isk_892[k]
                    + f_1 * ksi0_889[k]
                    - f_2 * ksi1_889[k]
                    + f_3 * pc_y[k] * ksk_1144[k];

        t_1432[k] = f_17 * isk_856[k]
                    + f_3 * pc_z[k] * ksk_1144[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pc_y, isk_894, isk_895, isk_896, ksi0_891, \
                         ksi0_892, ksi0_893, ksi1_891, ksi1_892, ksi1_893, ksk_1146, ksk_1147, \
                         ksk_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_18 * isk_894[k]
                    + f_12 * ksi0_891[k]
                    - f_13 * ksi1_891[k]
                    + f_3 * pc_y[k] * ksk_1146[k];

        t_1434[k] = f_18 * isk_895[k]
                    + f_10 * ksi0_892[k]
                    - f_11 * ksi1_892[k]
                    + f_3 * pc_y[k] * ksk_1147[k];

        t_1435[k] = f_18 * isk_896[k]
                    + f_8 * ksi0_893[k]
                    - f_9 * ksi1_893[k]
                    + f_3 * pc_y[k] * ksk_1148[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, pc_y, isk_897, isk_898, isk_899, ksi0_894, \
                         ksi0_895, ksi1_894, ksi1_895, ksk_1149, ksk_1150, \
                         ksk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_18 * isk_897[k]
                    + f_6 * ksi0_894[k]
                    - f_7 * ksi1_894[k]
                    + f_3 * pc_y[k] * ksk_1149[k];

        t_1437[k] = f_18 * isk_898[k]
                    + f_4 * ksi0_895[k]
                    - f_5 * ksi1_895[k]
                    + f_3 * pc_y[k] * ksk_1150[k];

        t_1438[k] = f_18 * isk_899[k]
                    + f_3 * pc_y[k] * ksk_1151[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, pc_x, pc_z, isk_863, ksi0_895, ksi0_896, \
                         ksi0_897, ksi1_895, ksi1_896, ksi1_897, ksk_1151, ksk_1152, \
                         ksk_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = f_17 * isk_863[k]
                    + f_1 * ksi0_895[k]
                    - f_2 * ksi1_895[k]
                    + f_3 * pc_z[k] * ksk_1151[k];

        t_1440[k] = f_1 * ksi0_896[k]
                    - f_2 * ksi1_896[k]
                    + f_3 * pc_x[k] * ksk_1152[k];

        t_1441[k] = f_21 * ksi0_897[k]
                    - f_22 * ksi1_897[k]
                    + f_3 * pc_x[k] * ksk_1153[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, pc_x, ksi0_898, ksi0_899, ksi0_900, ksi1_898, \
                         ksi1_899, ksi1_900, ksk_1154, ksk_1155, \
                         ksk_1156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_21 * ksi0_898[k]
                    - f_22 * ksi1_898[k]
                    + f_3 * pc_x[k] * ksk_1154[k];

        t_1443[k] = f_12 * ksi0_899[k]
                    - f_13 * ksi1_899[k]
                    + f_3 * pc_x[k] * ksk_1155[k];

        t_1444[k] = f_12 * ksi0_900[k]
                    - f_13 * ksi1_900[k]
                    + f_3 * pc_x[k] * ksk_1156[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pc_x, ksi0_901, ksi0_902, ksi0_903, ksi1_901, \
                         ksi1_902, ksi1_903, ksk_1157, ksk_1158, \
                         ksk_1159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_12 * ksi0_901[k]
                    - f_13 * ksi1_901[k]
                    + f_3 * pc_x[k] * ksk_1157[k];

        t_1446[k] = f_10 * ksi0_902[k]
                    - f_11 * ksi1_902[k]
                    + f_3 * pc_x[k] * ksk_1158[k];

        t_1447[k] = f_10 * ksi0_903[k]
                    - f_11 * ksi1_903[k]
                    + f_3 * pc_x[k] * ksk_1159[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, ksi0_904, ksi0_905, ksi0_906, ksi1_904, \
                         ksi1_905, ksi1_906, ksk_1160, ksk_1161, \
                         ksk_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_10 * ksi0_904[k]
                    - f_11 * ksi1_904[k]
                    + f_3 * pc_x[k] * ksk_1160[k];

        t_1449[k] = f_10 * ksi0_905[k]
                    - f_11 * ksi1_905[k]
                    + f_3 * pc_x[k] * ksk_1161[k];

        t_1450[k] = f_8 * ksi0_906[k]
                    - f_9 * ksi1_906[k]
                    + f_3 * pc_x[k] * ksk_1162[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_x, ksi0_907, ksi0_908, ksi0_909, ksi1_907, \
                         ksi1_908, ksi1_909, ksk_1163, ksk_1164, \
                         ksk_1165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_8 * ksi0_907[k]
                    - f_9 * ksi1_907[k]
                    + f_3 * pc_x[k] * ksk_1163[k];

        t_1452[k] = f_8 * ksi0_908[k]
                    - f_9 * ksi1_908[k]
                    + f_3 * pc_x[k] * ksk_1164[k];

        t_1453[k] = f_8 * ksi0_909[k]
                    - f_9 * ksi1_909[k]
                    + f_3 * pc_x[k] * ksk_1165[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, ksi0_910, ksi0_911, ksi0_912, ksi1_910, \
                         ksi1_911, ksi1_912, ksk_1166, ksk_1167, \
                         ksk_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_8 * ksi0_910[k]
                    - f_9 * ksi1_910[k]
                    + f_3 * pc_x[k] * ksk_1166[k];

        t_1455[k] = f_6 * ksi0_911[k]
                    - f_7 * ksi1_911[k]
                    + f_3 * pc_x[k] * ksk_1167[k];

        t_1456[k] = f_6 * ksi0_912[k]
                    - f_7 * ksi1_912[k]
                    + f_3 * pc_x[k] * ksk_1168[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pc_x, ksi0_913, ksi0_914, ksi0_915, ksi1_913, \
                         ksi1_914, ksi1_915, ksk_1169, ksk_1170, \
                         ksk_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_6 * ksi0_913[k]
                    - f_7 * ksi1_913[k]
                    + f_3 * pc_x[k] * ksk_1169[k];

        t_1458[k] = f_6 * ksi0_914[k]
                    - f_7 * ksi1_914[k]
                    + f_3 * pc_x[k] * ksk_1170[k];

        t_1459[k] = f_6 * ksi0_915[k]
                    - f_7 * ksi1_915[k]
                    + f_3 * pc_x[k] * ksk_1171[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, pc_x, ksi0_916, ksi0_917, ksi0_918, ksi1_916, \
                         ksi1_917, ksi1_918, ksk_1172, ksk_1173, \
                         ksk_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_6 * ksi0_916[k]
                    - f_7 * ksi1_916[k]
                    + f_3 * pc_x[k] * ksk_1172[k];

        t_1461[k] = f_4 * ksi0_917[k]
                    - f_5 * ksi1_917[k]
                    + f_3 * pc_x[k] * ksk_1173[k];

        t_1462[k] = f_4 * ksi0_918[k]
                    - f_5 * ksi1_918[k]
                    + f_3 * pc_x[k] * ksk_1174[k];
    }

#pragma omp simd aligned(t_1463, t_1464, t_1465, pc_x, ksi0_919, ksi0_920, ksi0_921, ksi1_919, \
                         ksi1_920, ksi1_921, ksk_1175, ksk_1176, \
                         ksk_1177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1463[k] = f_4 * ksi0_919[k]
                    - f_5 * ksi1_919[k]
                    + f_3 * pc_x[k] * ksk_1175[k];

        t_1464[k] = f_4 * ksi0_920[k]
                    - f_5 * ksi1_920[k]
                    + f_3 * pc_x[k] * ksk_1176[k];

        t_1465[k] = f_4 * ksi0_921[k]
                    - f_5 * ksi1_921[k]
                    + f_3 * pc_x[k] * ksk_1177[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, t_1469, t_1470, pc_x, ksi0_922, ksi0_923, \
                         ksi1_922, ksi1_923, ksk_1178, ksk_1179, ksk_1180, ksk_1181, \
                         ksk_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_4 * ksi0_922[k]
                    - f_5 * ksi1_922[k]
                    + f_3 * pc_x[k] * ksk_1178[k];

        t_1467[k] = f_4 * ksi0_923[k]
                    - f_5 * ksi1_923[k]
                    + f_3 * pc_x[k] * ksk_1179[k];

        t_1468[k] = f_3 * pc_x[k] * ksk_1180[k];

        t_1469[k] = f_3 * pc_x[k] * ksk_1181[k];

        t_1470[k] = f_3 * pc_x[k] * ksk_1182[k];
    }

#pragma omp simd aligned(t_1471, t_1472, t_1473, t_1474, t_1475, pc_x, ksk_1183, ksk_1184, \
                         ksk_1185, ksk_1186, ksk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1471[k] = f_3 * pc_x[k] * ksk_1183[k];

        t_1472[k] = f_3 * pc_x[k] * ksk_1184[k];

        t_1473[k] = f_3 * pc_x[k] * ksk_1185[k];

        t_1474[k] = f_3 * pc_x[k] * ksk_1186[k];

        t_1475[k] = f_3 * pc_x[k] * ksk_1187[k];
    }

#pragma omp simd aligned(t_1476, t_1477, t_1478, pc_y, pc_z, isk_892, isk_928, isk_930, \
                         ksi0_917, ksi0_919, ksi1_917, ksi1_919, ksk_1180, \
                         ksk_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1476[k] = f_17 * isk_928[k]
                    + f_1 * ksi0_917[k]
                    - f_2 * ksi1_917[k]
                    + f_3 * pc_y[k] * ksk_1180[k];

        t_1477[k] = f_18 * isk_892[k]
                    + f_3 * pc_z[k] * ksk_1180[k];

        t_1478[k] = f_17 * isk_930[k]
                    + f_12 * ksi0_919[k]
                    - f_13 * ksi1_919[k]
                    + f_3 * pc_y[k] * ksk_1182[k];
    }

#pragma omp simd aligned(t_1479, t_1480, t_1481, pc_y, isk_931, isk_932, isk_933, ksi0_920, \
                         ksi0_921, ksi0_922, ksi1_920, ksi1_921, ksi1_922, ksk_1183, ksk_1184, \
                         ksk_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1479[k] = f_17 * isk_931[k]
                    + f_10 * ksi0_920[k]
                    - f_11 * ksi1_920[k]
                    + f_3 * pc_y[k] * ksk_1183[k];

        t_1480[k] = f_17 * isk_932[k]
                    + f_8 * ksi0_921[k]
                    - f_9 * ksi1_921[k]
                    + f_3 * pc_y[k] * ksk_1184[k];

        t_1481[k] = f_17 * isk_933[k]
                    + f_6 * ksi0_922[k]
                    - f_7 * ksi1_922[k]
                    + f_3 * pc_y[k] * ksk_1185[k];
    }

#pragma omp simd aligned(t_1482, t_1483, t_1484, pc_y, pc_z, isk_899, isk_934, isk_935, \
                         ksi0_923, ksi1_923, ksk_1186, ksk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1482[k] = f_17 * isk_934[k]
                    + f_4 * ksi0_923[k]
                    - f_5 * ksi1_923[k]
                    + f_3 * pc_y[k] * ksk_1186[k];

        t_1483[k] = f_17 * isk_935[k]
                    + f_3 * pc_y[k] * ksk_1187[k];

        t_1484[k] = f_18 * isk_899[k]
                    + f_1 * ksi0_923[k]
                    - f_2 * ksi1_923[k]
                    + f_3 * pc_z[k] * ksk_1187[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, pc_x, ksi0_924, ksi0_925, ksi0_926, ksi1_924, \
                         ksi1_925, ksi1_926, ksk_1188, ksk_1189, \
                         ksk_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = f_1 * ksi0_924[k]
                    - f_2 * ksi1_924[k]
                    + f_3 * pc_x[k] * ksk_1188[k];

        t_1486[k] = f_21 * ksi0_925[k]
                    - f_22 * ksi1_925[k]
                    + f_3 * pc_x[k] * ksk_1189[k];

        t_1487[k] = f_21 * ksi0_926[k]
                    - f_22 * ksi1_926[k]
                    + f_3 * pc_x[k] * ksk_1190[k];
    }

#pragma omp simd aligned(t_1488, t_1489, t_1490, pc_x, ksi0_927, ksi0_928, ksi0_929, ksi1_927, \
                         ksi1_928, ksi1_929, ksk_1191, ksk_1192, \
                         ksk_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1488[k] = f_12 * ksi0_927[k]
                    - f_13 * ksi1_927[k]
                    + f_3 * pc_x[k] * ksk_1191[k];

        t_1489[k] = f_12 * ksi0_928[k]
                    - f_13 * ksi1_928[k]
                    + f_3 * pc_x[k] * ksk_1192[k];

        t_1490[k] = f_12 * ksi0_929[k]
                    - f_13 * ksi1_929[k]
                    + f_3 * pc_x[k] * ksk_1193[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, pc_x, ksi0_930, ksi0_931, ksi0_932, ksi1_930, \
                         ksi1_931, ksi1_932, ksk_1194, ksk_1195, \
                         ksk_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = f_10 * ksi0_930[k]
                    - f_11 * ksi1_930[k]
                    + f_3 * pc_x[k] * ksk_1194[k];

        t_1492[k] = f_10 * ksi0_931[k]
                    - f_11 * ksi1_931[k]
                    + f_3 * pc_x[k] * ksk_1195[k];

        t_1493[k] = f_10 * ksi0_932[k]
                    - f_11 * ksi1_932[k]
                    + f_3 * pc_x[k] * ksk_1196[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, pc_x, ksi0_933, ksi0_934, ksi0_935, ksi1_933, \
                         ksi1_934, ksi1_935, ksk_1197, ksk_1198, \
                         ksk_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = f_10 * ksi0_933[k]
                    - f_11 * ksi1_933[k]
                    + f_3 * pc_x[k] * ksk_1197[k];

        t_1495[k] = f_8 * ksi0_934[k]
                    - f_9 * ksi1_934[k]
                    + f_3 * pc_x[k] * ksk_1198[k];

        t_1496[k] = f_8 * ksi0_935[k]
                    - f_9 * ksi1_935[k]
                    + f_3 * pc_x[k] * ksk_1199[k];
    }

#pragma omp simd aligned(t_1497, t_1498, t_1499, pc_x, ksi0_936, ksi0_937, ksi0_938, ksi1_936, \
                         ksi1_937, ksi1_938, ksk_1200, ksk_1201, \
                         ksk_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1497[k] = f_8 * ksi0_936[k]
                    - f_9 * ksi1_936[k]
                    + f_3 * pc_x[k] * ksk_1200[k];

        t_1498[k] = f_8 * ksi0_937[k]
                    - f_9 * ksi1_937[k]
                    + f_3 * pc_x[k] * ksk_1201[k];

        t_1499[k] = f_8 * ksi0_938[k]
                    - f_9 * ksi1_938[k]
                    + f_3 * pc_x[k] * ksk_1202[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, pc_x, ksi0_939, ksi0_940, ksi0_941, ksi1_939, \
                         ksi1_940, ksi1_941, ksk_1203, ksk_1204, \
                         ksk_1205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_6 * ksi0_939[k]
                    - f_7 * ksi1_939[k]
                    + f_3 * pc_x[k] * ksk_1203[k];

        t_1501[k] = f_6 * ksi0_940[k]
                    - f_7 * ksi1_940[k]
                    + f_3 * pc_x[k] * ksk_1204[k];

        t_1502[k] = f_6 * ksi0_941[k]
                    - f_7 * ksi1_941[k]
                    + f_3 * pc_x[k] * ksk_1205[k];
    }

#pragma omp simd aligned(t_1503, t_1504, t_1505, pc_x, ksi0_942, ksi0_943, ksi0_944, ksi1_942, \
                         ksi1_943, ksi1_944, ksk_1206, ksk_1207, \
                         ksk_1208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1503[k] = f_6 * ksi0_942[k]
                    - f_7 * ksi1_942[k]
                    + f_3 * pc_x[k] * ksk_1206[k];

        t_1504[k] = f_6 * ksi0_943[k]
                    - f_7 * ksi1_943[k]
                    + f_3 * pc_x[k] * ksk_1207[k];

        t_1505[k] = f_6 * ksi0_944[k]
                    - f_7 * ksi1_944[k]
                    + f_3 * pc_x[k] * ksk_1208[k];
    }

#pragma omp simd aligned(t_1506, t_1507, t_1508, pc_x, ksi0_945, ksi0_946, ksi0_947, ksi1_945, \
                         ksi1_946, ksi1_947, ksk_1209, ksk_1210, \
                         ksk_1211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1506[k] = f_4 * ksi0_945[k]
                    - f_5 * ksi1_945[k]
                    + f_3 * pc_x[k] * ksk_1209[k];

        t_1507[k] = f_4 * ksi0_946[k]
                    - f_5 * ksi1_946[k]
                    + f_3 * pc_x[k] * ksk_1210[k];

        t_1508[k] = f_4 * ksi0_947[k]
                    - f_5 * ksi1_947[k]
                    + f_3 * pc_x[k] * ksk_1211[k];
    }
}

static auto
compute_prim_ksl_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t isl0,
                                                           const size_t isk, const size_t isl1,
                                                           const size_t ksi0, const size_t ksi1,
                                                           const size_t ksk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_21 = 3.0 / gamma;
    const auto f_22 = 3.0 * p / (gamma * q);
    const auto f_23 = 4.0 / q;

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

    const auto *isl0_1215 = buffer.data(isl0 + 1215);
    const auto *isl0_1217 = buffer.data(isl0 + 1217);
    const auto *isl0_1220 = buffer.data(isl0 + 1220);
    const auto *isl0_1224 = buffer.data(isl0 + 1224);
    const auto *isl0_1229 = buffer.data(isl0 + 1229);
    const auto *isl0_1235 = buffer.data(isl0 + 1235);
    const auto *isl0_1242 = buffer.data(isl0 + 1242);
    const auto *isl0_1251 = buffer.data(isl0 + 1251);
    const auto *isl0_1253 = buffer.data(isl0 + 1253);
    const auto *isl0_1254 = buffer.data(isl0 + 1254);
    const auto *isl0_1255 = buffer.data(isl0 + 1255);
    const auto *isl0_1256 = buffer.data(isl0 + 1256);
    const auto *isl0_1257 = buffer.data(isl0 + 1257);
    const auto *isl0_1259 = buffer.data(isl0 + 1259);

    const auto *isk_928 = buffer.data(isk + 928);
    const auto *isk_935 = buffer.data(isk + 935);
    const auto *isk_964 = buffer.data(isk + 964);
    const auto *isk_966 = buffer.data(isk + 966);
    const auto *isk_967 = buffer.data(isk + 967);
    const auto *isk_968 = buffer.data(isk + 968);
    const auto *isk_969 = buffer.data(isk + 969);
    const auto *isk_970 = buffer.data(isk + 970);
    const auto *isk_971 = buffer.data(isk + 971);
    const auto *isk_1000 = buffer.data(isk + 1000);
    const auto *isk_1002 = buffer.data(isk + 1002);
    const auto *isk_1003 = buffer.data(isk + 1003);
    const auto *isk_1004 = buffer.data(isk + 1004);
    const auto *isk_1005 = buffer.data(isk + 1005);
    const auto *isk_1006 = buffer.data(isk + 1006);
    const auto *isk_1007 = buffer.data(isk + 1007);

    const auto *isl1_1215 = buffer.data(isl1 + 1215);
    const auto *isl1_1217 = buffer.data(isl1 + 1217);
    const auto *isl1_1220 = buffer.data(isl1 + 1220);
    const auto *isl1_1224 = buffer.data(isl1 + 1224);
    const auto *isl1_1229 = buffer.data(isl1 + 1229);
    const auto *isl1_1235 = buffer.data(isl1 + 1235);
    const auto *isl1_1242 = buffer.data(isl1 + 1242);
    const auto *isl1_1251 = buffer.data(isl1 + 1251);
    const auto *isl1_1253 = buffer.data(isl1 + 1253);
    const auto *isl1_1254 = buffer.data(isl1 + 1254);
    const auto *isl1_1255 = buffer.data(isl1 + 1255);
    const auto *isl1_1256 = buffer.data(isl1 + 1256);
    const auto *isl1_1257 = buffer.data(isl1 + 1257);
    const auto *isl1_1259 = buffer.data(isl1 + 1259);

    const auto *ksi0_945 = buffer.data(ksi0 + 945);
    const auto *ksi0_947 = buffer.data(ksi0 + 947);
    const auto *ksi0_948 = buffer.data(ksi0 + 948);
    const auto *ksi0_949 = buffer.data(ksi0 + 949);
    const auto *ksi0_950 = buffer.data(ksi0 + 950);
    const auto *ksi0_951 = buffer.data(ksi0 + 951);
    const auto *ksi0_953 = buffer.data(ksi0 + 953);
    const auto *ksi0_955 = buffer.data(ksi0 + 955);
    const auto *ksi0_956 = buffer.data(ksi0 + 956);
    const auto *ksi0_958 = buffer.data(ksi0 + 958);
    const auto *ksi0_959 = buffer.data(ksi0 + 959);
    const auto *ksi0_960 = buffer.data(ksi0 + 960);
    const auto *ksi0_962 = buffer.data(ksi0 + 962);
    const auto *ksi0_963 = buffer.data(ksi0 + 963);
    const auto *ksi0_964 = buffer.data(ksi0 + 964);
    const auto *ksi0_965 = buffer.data(ksi0 + 965);
    const auto *ksi0_967 = buffer.data(ksi0 + 967);
    const auto *ksi0_968 = buffer.data(ksi0 + 968);
    const auto *ksi0_969 = buffer.data(ksi0 + 969);
    const auto *ksi0_970 = buffer.data(ksi0 + 970);
    const auto *ksi0_971 = buffer.data(ksi0 + 971);
    const auto *ksi0_973 = buffer.data(ksi0 + 973);
    const auto *ksi0_974 = buffer.data(ksi0 + 974);
    const auto *ksi0_975 = buffer.data(ksi0 + 975);
    const auto *ksi0_976 = buffer.data(ksi0 + 976);
    const auto *ksi0_977 = buffer.data(ksi0 + 977);
    const auto *ksi0_978 = buffer.data(ksi0 + 978);
    const auto *ksi0_980 = buffer.data(ksi0 + 980);
    const auto *ksi0_982 = buffer.data(ksi0 + 982);
    const auto *ksi0_983 = buffer.data(ksi0 + 983);
    const auto *ksi0_985 = buffer.data(ksi0 + 985);
    const auto *ksi0_986 = buffer.data(ksi0 + 986);
    const auto *ksi0_987 = buffer.data(ksi0 + 987);
    const auto *ksi0_989 = buffer.data(ksi0 + 989);
    const auto *ksi0_990 = buffer.data(ksi0 + 990);
    const auto *ksi0_991 = buffer.data(ksi0 + 991);
    const auto *ksi0_992 = buffer.data(ksi0 + 992);
    const auto *ksi0_994 = buffer.data(ksi0 + 994);
    const auto *ksi0_995 = buffer.data(ksi0 + 995);
    const auto *ksi0_996 = buffer.data(ksi0 + 996);
    const auto *ksi0_997 = buffer.data(ksi0 + 997);
    const auto *ksi0_998 = buffer.data(ksi0 + 998);
    const auto *ksi0_1000 = buffer.data(ksi0 + 1000);
    const auto *ksi0_1001 = buffer.data(ksi0 + 1001);
    const auto *ksi0_1002 = buffer.data(ksi0 + 1002);
    const auto *ksi0_1003 = buffer.data(ksi0 + 1003);
    const auto *ksi0_1004 = buffer.data(ksi0 + 1004);
    const auto *ksi0_1005 = buffer.data(ksi0 + 1005);
    const auto *ksi0_1006 = buffer.data(ksi0 + 1006);
    const auto *ksi0_1007 = buffer.data(ksi0 + 1007);

    const auto *ksi1_945 = buffer.data(ksi1 + 945);
    const auto *ksi1_947 = buffer.data(ksi1 + 947);
    const auto *ksi1_948 = buffer.data(ksi1 + 948);
    const auto *ksi1_949 = buffer.data(ksi1 + 949);
    const auto *ksi1_950 = buffer.data(ksi1 + 950);
    const auto *ksi1_951 = buffer.data(ksi1 + 951);
    const auto *ksi1_953 = buffer.data(ksi1 + 953);
    const auto *ksi1_955 = buffer.data(ksi1 + 955);
    const auto *ksi1_956 = buffer.data(ksi1 + 956);
    const auto *ksi1_958 = buffer.data(ksi1 + 958);
    const auto *ksi1_959 = buffer.data(ksi1 + 959);
    const auto *ksi1_960 = buffer.data(ksi1 + 960);
    const auto *ksi1_962 = buffer.data(ksi1 + 962);
    const auto *ksi1_963 = buffer.data(ksi1 + 963);
    const auto *ksi1_964 = buffer.data(ksi1 + 964);
    const auto *ksi1_965 = buffer.data(ksi1 + 965);
    const auto *ksi1_967 = buffer.data(ksi1 + 967);
    const auto *ksi1_968 = buffer.data(ksi1 + 968);
    const auto *ksi1_969 = buffer.data(ksi1 + 969);
    const auto *ksi1_970 = buffer.data(ksi1 + 970);
    const auto *ksi1_971 = buffer.data(ksi1 + 971);
    const auto *ksi1_973 = buffer.data(ksi1 + 973);
    const auto *ksi1_974 = buffer.data(ksi1 + 974);
    const auto *ksi1_975 = buffer.data(ksi1 + 975);
    const auto *ksi1_976 = buffer.data(ksi1 + 976);
    const auto *ksi1_977 = buffer.data(ksi1 + 977);
    const auto *ksi1_978 = buffer.data(ksi1 + 978);
    const auto *ksi1_980 = buffer.data(ksi1 + 980);
    const auto *ksi1_982 = buffer.data(ksi1 + 982);
    const auto *ksi1_983 = buffer.data(ksi1 + 983);
    const auto *ksi1_985 = buffer.data(ksi1 + 985);
    const auto *ksi1_986 = buffer.data(ksi1 + 986);
    const auto *ksi1_987 = buffer.data(ksi1 + 987);
    const auto *ksi1_989 = buffer.data(ksi1 + 989);
    const auto *ksi1_990 = buffer.data(ksi1 + 990);
    const auto *ksi1_991 = buffer.data(ksi1 + 991);
    const auto *ksi1_992 = buffer.data(ksi1 + 992);
    const auto *ksi1_994 = buffer.data(ksi1 + 994);
    const auto *ksi1_995 = buffer.data(ksi1 + 995);
    const auto *ksi1_996 = buffer.data(ksi1 + 996);
    const auto *ksi1_997 = buffer.data(ksi1 + 997);
    const auto *ksi1_998 = buffer.data(ksi1 + 998);
    const auto *ksi1_1000 = buffer.data(ksi1 + 1000);
    const auto *ksi1_1001 = buffer.data(ksi1 + 1001);
    const auto *ksi1_1002 = buffer.data(ksi1 + 1002);
    const auto *ksi1_1003 = buffer.data(ksi1 + 1003);
    const auto *ksi1_1004 = buffer.data(ksi1 + 1004);
    const auto *ksi1_1005 = buffer.data(ksi1 + 1005);
    const auto *ksi1_1006 = buffer.data(ksi1 + 1006);
    const auto *ksi1_1007 = buffer.data(ksi1 + 1007);

    const auto *ksk_1212 = buffer.data(ksk + 1212);
    const auto *ksk_1213 = buffer.data(ksk + 1213);
    const auto *ksk_1214 = buffer.data(ksk + 1214);
    const auto *ksk_1215 = buffer.data(ksk + 1215);
    const auto *ksk_1216 = buffer.data(ksk + 1216);
    const auto *ksk_1217 = buffer.data(ksk + 1217);
    const auto *ksk_1218 = buffer.data(ksk + 1218);
    const auto *ksk_1219 = buffer.data(ksk + 1219);
    const auto *ksk_1220 = buffer.data(ksk + 1220);
    const auto *ksk_1221 = buffer.data(ksk + 1221);
    const auto *ksk_1222 = buffer.data(ksk + 1222);
    const auto *ksk_1223 = buffer.data(ksk + 1223);
    const auto *ksk_1225 = buffer.data(ksk + 1225);
    const auto *ksk_1227 = buffer.data(ksk + 1227);
    const auto *ksk_1228 = buffer.data(ksk + 1228);
    const auto *ksk_1230 = buffer.data(ksk + 1230);
    const auto *ksk_1231 = buffer.data(ksk + 1231);
    const auto *ksk_1232 = buffer.data(ksk + 1232);
    const auto *ksk_1234 = buffer.data(ksk + 1234);
    const auto *ksk_1235 = buffer.data(ksk + 1235);
    const auto *ksk_1236 = buffer.data(ksk + 1236);
    const auto *ksk_1237 = buffer.data(ksk + 1237);
    const auto *ksk_1239 = buffer.data(ksk + 1239);
    const auto *ksk_1240 = buffer.data(ksk + 1240);
    const auto *ksk_1241 = buffer.data(ksk + 1241);
    const auto *ksk_1242 = buffer.data(ksk + 1242);
    const auto *ksk_1243 = buffer.data(ksk + 1243);
    const auto *ksk_1245 = buffer.data(ksk + 1245);
    const auto *ksk_1246 = buffer.data(ksk + 1246);
    const auto *ksk_1247 = buffer.data(ksk + 1247);
    const auto *ksk_1248 = buffer.data(ksk + 1248);
    const auto *ksk_1249 = buffer.data(ksk + 1249);
    const auto *ksk_1250 = buffer.data(ksk + 1250);
    const auto *ksk_1252 = buffer.data(ksk + 1252);
    const auto *ksk_1253 = buffer.data(ksk + 1253);
    const auto *ksk_1254 = buffer.data(ksk + 1254);
    const auto *ksk_1255 = buffer.data(ksk + 1255);
    const auto *ksk_1256 = buffer.data(ksk + 1256);
    const auto *ksk_1257 = buffer.data(ksk + 1257);
    const auto *ksk_1258 = buffer.data(ksk + 1258);
    const auto *ksk_1259 = buffer.data(ksk + 1259);
    const auto *ksk_1260 = buffer.data(ksk + 1260);
    const auto *ksk_1262 = buffer.data(ksk + 1262);
    const auto *ksk_1263 = buffer.data(ksk + 1263);
    const auto *ksk_1265 = buffer.data(ksk + 1265);
    const auto *ksk_1266 = buffer.data(ksk + 1266);
    const auto *ksk_1267 = buffer.data(ksk + 1267);
    const auto *ksk_1269 = buffer.data(ksk + 1269);
    const auto *ksk_1270 = buffer.data(ksk + 1270);
    const auto *ksk_1271 = buffer.data(ksk + 1271);
    const auto *ksk_1272 = buffer.data(ksk + 1272);
    const auto *ksk_1274 = buffer.data(ksk + 1274);
    const auto *ksk_1275 = buffer.data(ksk + 1275);
    const auto *ksk_1276 = buffer.data(ksk + 1276);
    const auto *ksk_1277 = buffer.data(ksk + 1277);
    const auto *ksk_1278 = buffer.data(ksk + 1278);
    const auto *ksk_1280 = buffer.data(ksk + 1280);
    const auto *ksk_1281 = buffer.data(ksk + 1281);
    const auto *ksk_1282 = buffer.data(ksk + 1282);
    const auto *ksk_1283 = buffer.data(ksk + 1283);
    const auto *ksk_1284 = buffer.data(ksk + 1284);
    const auto *ksk_1285 = buffer.data(ksk + 1285);
    const auto *ksk_1287 = buffer.data(ksk + 1287);
    const auto *ksk_1288 = buffer.data(ksk + 1288);
    const auto *ksk_1289 = buffer.data(ksk + 1289);
    const auto *ksk_1290 = buffer.data(ksk + 1290);
    const auto *ksk_1291 = buffer.data(ksk + 1291);
    const auto *ksk_1292 = buffer.data(ksk + 1292);
    const auto *ksk_1293 = buffer.data(ksk + 1293);
    const auto *ksk_1294 = buffer.data(ksk + 1294);
    const auto *ksk_1295 = buffer.data(ksk + 1295);

#pragma omp simd aligned(t_1509, t_1510, t_1511, pc_x, ksi0_948, ksi0_949, ksi0_950, ksi1_948, \
                         ksi1_949, ksi1_950, ksk_1212, ksk_1213, \
                         ksk_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1509[k] = f_4 * ksi0_948[k]
                    - f_5 * ksi1_948[k]
                    + f_3 * pc_x[k] * ksk_1212[k];

        t_1510[k] = f_4 * ksi0_949[k]
                    - f_5 * ksi1_949[k]
                    + f_3 * pc_x[k] * ksk_1213[k];

        t_1511[k] = f_4 * ksi0_950[k]
                    - f_5 * ksi1_950[k]
                    + f_3 * pc_x[k] * ksk_1214[k];
    }

#pragma omp simd aligned(t_1512, t_1513, t_1514, t_1515, t_1516, t_1517, pc_x, ksi0_951, \
                         ksi1_951, ksk_1215, ksk_1216, ksk_1217, ksk_1218, ksk_1219, \
                         ksk_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1512[k] = f_4 * ksi0_951[k]
                    - f_5 * ksi1_951[k]
                    + f_3 * pc_x[k] * ksk_1215[k];

        t_1513[k] = f_3 * pc_x[k] * ksk_1216[k];

        t_1514[k] = f_3 * pc_x[k] * ksk_1217[k];

        t_1515[k] = f_3 * pc_x[k] * ksk_1218[k];

        t_1516[k] = f_3 * pc_x[k] * ksk_1219[k];

        t_1517[k] = f_3 * pc_x[k] * ksk_1220[k];
    }

#pragma omp simd aligned(t_1518, t_1519, t_1520, t_1521, t_1522, pc_x, pc_y, pc_z, isk_928, \
                         isk_964, ksi0_945, ksi1_945, ksk_1216, ksk_1221, ksk_1222, \
                         ksk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = f_3 * pc_x[k] * ksk_1221[k];

        t_1519[k] = f_3 * pc_x[k] * ksk_1222[k];

        t_1520[k] = f_3 * pc_x[k] * ksk_1223[k];

        t_1521[k] = f_16 * isk_964[k]
                    + f_1 * ksi0_945[k]
                    - f_2 * ksi1_945[k]
                    + f_3 * pc_y[k] * ksk_1216[k];

        t_1522[k] = f_19 * isk_928[k]
                    + f_3 * pc_z[k] * ksk_1216[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pc_y, isk_966, isk_967, isk_968, ksi0_947, \
                         ksi0_948, ksi0_949, ksi1_947, ksi1_948, ksi1_949, ksk_1218, ksk_1219, \
                         ksk_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_16 * isk_966[k]
                    + f_12 * ksi0_947[k]
                    - f_13 * ksi1_947[k]
                    + f_3 * pc_y[k] * ksk_1218[k];

        t_1524[k] = f_16 * isk_967[k]
                    + f_10 * ksi0_948[k]
                    - f_11 * ksi1_948[k]
                    + f_3 * pc_y[k] * ksk_1219[k];

        t_1525[k] = f_16 * isk_968[k]
                    + f_8 * ksi0_949[k]
                    - f_9 * ksi1_949[k]
                    + f_3 * pc_y[k] * ksk_1220[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, pc_y, isk_969, isk_970, isk_971, ksi0_950, \
                         ksi0_951, ksi1_950, ksi1_951, ksk_1221, ksk_1222, \
                         ksk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_16 * isk_969[k]
                    + f_6 * ksi0_950[k]
                    - f_7 * ksi1_950[k]
                    + f_3 * pc_y[k] * ksk_1221[k];

        t_1527[k] = f_16 * isk_970[k]
                    + f_4 * ksi0_951[k]
                    - f_5 * ksi1_951[k]
                    + f_3 * pc_y[k] * ksk_1222[k];

        t_1528[k] = f_16 * isk_971[k]
                    + f_3 * pc_y[k] * ksk_1223[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, pa_y, pc_x, pc_y, pc_z, isl0_1215, isk_935, \
                         isl1_1215, ksi0_951, ksi0_953, ksi1_951, ksi1_953, ksk_1223, \
                         ksk_1225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_19 * isk_935[k]
                    + f_1 * ksi0_951[k]
                    - f_2 * ksi1_951[k]
                    + f_3 * pc_z[k] * ksk_1223[k];

        t_1530[k] = pa_y[k] * isl0_1215[k]
                    - f_14 * pc_y[k] * isl1_1215[k];

        t_1531[k] = f_21 * ksi0_953[k]
                    - f_22 * ksi1_953[k]
                    + f_3 * pc_x[k] * ksk_1225[k];
    }

#pragma omp simd aligned(t_1532, t_1533, t_1534, pa_y, pc_x, pc_y, isl0_1217, isl1_1217, \
                         ksi0_955, ksi0_956, ksi1_955, ksi1_956, ksk_1227, \
                         ksk_1228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1532[k] = pa_y[k] * isl0_1217[k]
                    - f_14 * pc_y[k] * isl1_1217[k];

        t_1533[k] = f_12 * ksi0_955[k]
                    - f_13 * ksi1_955[k]
                    + f_3 * pc_x[k] * ksk_1227[k];

        t_1534[k] = f_12 * ksi0_956[k]
                    - f_13 * ksi1_956[k]
                    + f_3 * pc_x[k] * ksk_1228[k];
    }

#pragma omp simd aligned(t_1535, t_1536, t_1537, pa_y, pc_x, pc_y, isl0_1220, isl1_1220, \
                         ksi0_958, ksi0_959, ksi1_958, ksi1_959, ksk_1230, \
                         ksk_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1535[k] = pa_y[k] * isl0_1220[k]
                    - f_14 * pc_y[k] * isl1_1220[k];

        t_1536[k] = f_10 * ksi0_958[k]
                    - f_11 * ksi1_958[k]
                    + f_3 * pc_x[k] * ksk_1230[k];

        t_1537[k] = f_10 * ksi0_959[k]
                    - f_11 * ksi1_959[k]
                    + f_3 * pc_x[k] * ksk_1231[k];
    }

#pragma omp simd aligned(t_1538, t_1539, t_1540, pa_y, pc_x, pc_y, isl0_1224, isl1_1224, \
                         ksi0_960, ksi0_962, ksi1_960, ksi1_962, ksk_1232, \
                         ksk_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1538[k] = f_10 * ksi0_960[k]
                    - f_11 * ksi1_960[k]
                    + f_3 * pc_x[k] * ksk_1232[k];

        t_1539[k] = pa_y[k] * isl0_1224[k]
                    - f_14 * pc_y[k] * isl1_1224[k];

        t_1540[k] = f_8 * ksi0_962[k]
                    - f_9 * ksi1_962[k]
                    + f_3 * pc_x[k] * ksk_1234[k];
    }

#pragma omp simd aligned(t_1541, t_1542, t_1543, pc_x, ksi0_963, ksi0_964, ksi0_965, ksi1_963, \
                         ksi1_964, ksi1_965, ksk_1235, ksk_1236, \
                         ksk_1237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1541[k] = f_8 * ksi0_963[k]
                    - f_9 * ksi1_963[k]
                    + f_3 * pc_x[k] * ksk_1235[k];

        t_1542[k] = f_8 * ksi0_964[k]
                    - f_9 * ksi1_964[k]
                    + f_3 * pc_x[k] * ksk_1236[k];

        t_1543[k] = f_8 * ksi0_965[k]
                    - f_9 * ksi1_965[k]
                    + f_3 * pc_x[k] * ksk_1237[k];
    }

#pragma omp simd aligned(t_1544, t_1545, t_1546, pa_y, pc_x, pc_y, isl0_1229, isl1_1229, \
                         ksi0_967, ksi0_968, ksi1_967, ksi1_968, ksk_1239, \
                         ksk_1240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1544[k] = pa_y[k] * isl0_1229[k]
                    - f_14 * pc_y[k] * isl1_1229[k];

        t_1545[k] = f_6 * ksi0_967[k]
                    - f_7 * ksi1_967[k]
                    + f_3 * pc_x[k] * ksk_1239[k];

        t_1546[k] = f_6 * ksi0_968[k]
                    - f_7 * ksi1_968[k]
                    + f_3 * pc_x[k] * ksk_1240[k];
    }

#pragma omp simd aligned(t_1547, t_1548, t_1549, pc_x, ksi0_969, ksi0_970, ksi0_971, ksi1_969, \
                         ksi1_970, ksi1_971, ksk_1241, ksk_1242, \
                         ksk_1243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1547[k] = f_6 * ksi0_969[k]
                    - f_7 * ksi1_969[k]
                    + f_3 * pc_x[k] * ksk_1241[k];

        t_1548[k] = f_6 * ksi0_970[k]
                    - f_7 * ksi1_970[k]
                    + f_3 * pc_x[k] * ksk_1242[k];

        t_1549[k] = f_6 * ksi0_971[k]
                    - f_7 * ksi1_971[k]
                    + f_3 * pc_x[k] * ksk_1243[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, pa_y, pc_x, pc_y, isl0_1235, isl1_1235, \
                         ksi0_973, ksi0_974, ksi1_973, ksi1_974, ksk_1245, \
                         ksk_1246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = pa_y[k] * isl0_1235[k]
                    - f_14 * pc_y[k] * isl1_1235[k];

        t_1551[k] = f_4 * ksi0_973[k]
                    - f_5 * ksi1_973[k]
                    + f_3 * pc_x[k] * ksk_1245[k];

        t_1552[k] = f_4 * ksi0_974[k]
                    - f_5 * ksi1_974[k]
                    + f_3 * pc_x[k] * ksk_1246[k];
    }

#pragma omp simd aligned(t_1553, t_1554, t_1555, pc_x, ksi0_975, ksi0_976, ksi0_977, ksi1_975, \
                         ksi1_976, ksi1_977, ksk_1247, ksk_1248, \
                         ksk_1249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1553[k] = f_4 * ksi0_975[k]
                    - f_5 * ksi1_975[k]
                    + f_3 * pc_x[k] * ksk_1247[k];

        t_1554[k] = f_4 * ksi0_976[k]
                    - f_5 * ksi1_976[k]
                    + f_3 * pc_x[k] * ksk_1248[k];

        t_1555[k] = f_4 * ksi0_977[k]
                    - f_5 * ksi1_977[k]
                    + f_3 * pc_x[k] * ksk_1249[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, t_1559, t_1560, pa_y, pc_x, pc_y, isl0_1242, \
                         isl1_1242, ksi0_978, ksi1_978, ksk_1250, ksk_1252, ksk_1253, \
                         ksk_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_4 * ksi0_978[k]
                    - f_5 * ksi1_978[k]
                    + f_3 * pc_x[k] * ksk_1250[k];

        t_1557[k] = pa_y[k] * isl0_1242[k]
                    - f_14 * pc_y[k] * isl1_1242[k];

        t_1558[k] = f_3 * pc_x[k] * ksk_1252[k];

        t_1559[k] = f_3 * pc_x[k] * ksk_1253[k];

        t_1560[k] = f_3 * pc_x[k] * ksk_1254[k];
    }

#pragma omp simd aligned(t_1561, t_1562, t_1563, t_1564, t_1565, pc_x, ksk_1255, ksk_1256, \
                         ksk_1257, ksk_1258, ksk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1561[k] = f_3 * pc_x[k] * ksk_1255[k];

        t_1562[k] = f_3 * pc_x[k] * ksk_1256[k];

        t_1563[k] = f_3 * pc_x[k] * ksk_1257[k];

        t_1564[k] = f_3 * pc_x[k] * ksk_1258[k];

        t_1565[k] = f_3 * pc_x[k] * ksk_1259[k];
    }

#pragma omp simd aligned(t_1566, t_1567, t_1568, pa_y, pc_y, pc_z, isl0_1251, isl0_1253, \
                         isk_964, isk_1000, isk_1002, isl1_1251, isl1_1253, \
                         ksk_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1566[k] = pa_y[k] * isl0_1251[k]
                    + f_23 * isk_1000[k]
                    - f_14 * pc_y[k] * isl1_1251[k];

        t_1567[k] = f_20 * isk_964[k]
                    + f_3 * pc_z[k] * ksk_1252[k];

        t_1568[k] = pa_y[k] * isl0_1253[k]
                    + f_20 * isk_1002[k]
                    - f_14 * pc_y[k] * isl1_1253[k];
    }

#pragma omp simd aligned(t_1569, t_1570, t_1571, pa_y, pc_y, isl0_1254, isl0_1255, isl0_1256, \
                         isk_1003, isk_1004, isk_1005, isl1_1254, isl1_1255, \
                         isl1_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1569[k] = pa_y[k] * isl0_1254[k]
                    + f_19 * isk_1003[k]
                    - f_14 * pc_y[k] * isl1_1254[k];

        t_1570[k] = pa_y[k] * isl0_1255[k]
                    + f_18 * isk_1004[k]
                    - f_14 * pc_y[k] * isl1_1255[k];

        t_1571[k] = pa_y[k] * isl0_1256[k]
                    + f_17 * isk_1005[k]
                    - f_14 * pc_y[k] * isl1_1256[k];
    }

#pragma omp simd aligned(t_1572, t_1573, t_1574, pa_y, pc_y, isl0_1257, isl0_1259, isk_1006, \
                         isk_1007, isl1_1257, isl1_1259, ksk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1572[k] = pa_y[k] * isl0_1257[k]
                    + f_16 * isk_1006[k]
                    - f_14 * pc_y[k] * isl1_1257[k];

        t_1573[k] = f_15 * isk_1007[k]
                    + f_3 * pc_y[k] * ksk_1259[k];

        t_1574[k] = pa_y[k] * isl0_1259[k]
                    - f_14 * pc_y[k] * isl1_1259[k];
    }

#pragma omp simd aligned(t_1575, t_1576, t_1577, t_1578, t_1579, pc_x, pc_y, ksi0_980, \
                         ksi0_982, ksi0_983, ksi1_980, ksi1_982, ksi1_983, ksk_1260, ksk_1262, \
                         ksk_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1575[k] = f_1 * ksi0_980[k]
                    - f_2 * ksi1_980[k]
                    + f_3 * pc_x[k] * ksk_1260[k];

        t_1576[k] = f_3 * pc_y[k] * ksk_1260[k];

        t_1577[k] = f_21 * ksi0_982[k]
                    - f_22 * ksi1_982[k]
                    + f_3 * pc_x[k] * ksk_1262[k];

        t_1578[k] = f_12 * ksi0_983[k]
                    - f_13 * ksi1_983[k]
                    + f_3 * pc_x[k] * ksk_1263[k];

        t_1579[k] = f_3 * pc_y[k] * ksk_1262[k];
    }

#pragma omp simd aligned(t_1580, t_1581, t_1582, t_1583, pc_x, pc_y, ksi0_985, ksi0_986, \
                         ksi0_987, ksi1_985, ksi1_986, ksi1_987, ksk_1265, ksk_1266, \
                         ksk_1267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1580[k] = f_12 * ksi0_985[k]
                    - f_13 * ksi1_985[k]
                    + f_3 * pc_x[k] * ksk_1265[k];

        t_1581[k] = f_10 * ksi0_986[k]
                    - f_11 * ksi1_986[k]
                    + f_3 * pc_x[k] * ksk_1266[k];

        t_1582[k] = f_10 * ksi0_987[k]
                    - f_11 * ksi1_987[k]
                    + f_3 * pc_x[k] * ksk_1267[k];

        t_1583[k] = f_3 * pc_y[k] * ksk_1265[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, pc_x, ksi0_989, ksi0_990, ksi0_991, ksi1_989, \
                         ksi1_990, ksi1_991, ksk_1269, ksk_1270, \
                         ksk_1271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = f_10 * ksi0_989[k]
                    - f_11 * ksi1_989[k]
                    + f_3 * pc_x[k] * ksk_1269[k];

        t_1585[k] = f_8 * ksi0_990[k]
                    - f_9 * ksi1_990[k]
                    + f_3 * pc_x[k] * ksk_1270[k];

        t_1586[k] = f_8 * ksi0_991[k]
                    - f_9 * ksi1_991[k]
                    + f_3 * pc_x[k] * ksk_1271[k];
    }

#pragma omp simd aligned(t_1587, t_1588, t_1589, t_1590, pc_x, pc_y, ksi0_992, ksi0_994, \
                         ksi0_995, ksi1_992, ksi1_994, ksi1_995, ksk_1269, ksk_1272, ksk_1274, \
                         ksk_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1587[k] = f_8 * ksi0_992[k]
                    - f_9 * ksi1_992[k]
                    + f_3 * pc_x[k] * ksk_1272[k];

        t_1588[k] = f_3 * pc_y[k] * ksk_1269[k];

        t_1589[k] = f_8 * ksi0_994[k]
                    - f_9 * ksi1_994[k]
                    + f_3 * pc_x[k] * ksk_1274[k];

        t_1590[k] = f_6 * ksi0_995[k]
                    - f_7 * ksi1_995[k]
                    + f_3 * pc_x[k] * ksk_1275[k];
    }

#pragma omp simd aligned(t_1591, t_1592, t_1593, t_1594, pc_x, pc_y, ksi0_996, ksi0_997, \
                         ksi0_998, ksi1_996, ksi1_997, ksi1_998, ksk_1274, ksk_1276, ksk_1277, \
                         ksk_1278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1591[k] = f_6 * ksi0_996[k]
                    - f_7 * ksi1_996[k]
                    + f_3 * pc_x[k] * ksk_1276[k];

        t_1592[k] = f_6 * ksi0_997[k]
                    - f_7 * ksi1_997[k]
                    + f_3 * pc_x[k] * ksk_1277[k];

        t_1593[k] = f_6 * ksi0_998[k]
                    - f_7 * ksi1_998[k]
                    + f_3 * pc_x[k] * ksk_1278[k];

        t_1594[k] = f_3 * pc_y[k] * ksk_1274[k];
    }

#pragma omp simd aligned(t_1595, t_1596, t_1597, pc_x, ksi0_1000, ksi0_1001, ksi0_1002, \
                         ksi1_1000, ksi1_1001, ksi1_1002, ksk_1280, ksk_1281, \
                         ksk_1282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1595[k] = f_6 * ksi0_1000[k]
                    - f_7 * ksi1_1000[k]
                    + f_3 * pc_x[k] * ksk_1280[k];

        t_1596[k] = f_4 * ksi0_1001[k]
                    - f_5 * ksi1_1001[k]
                    + f_3 * pc_x[k] * ksk_1281[k];

        t_1597[k] = f_4 * ksi0_1002[k]
                    - f_5 * ksi1_1002[k]
                    + f_3 * pc_x[k] * ksk_1282[k];
    }

#pragma omp simd aligned(t_1598, t_1599, t_1600, t_1601, pc_x, pc_y, ksi0_1003, ksi0_1004, \
                         ksi0_1005, ksi1_1003, ksi1_1004, ksi1_1005, ksk_1280, ksk_1283, \
                         ksk_1284, ksk_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1598[k] = f_4 * ksi0_1003[k]
                    - f_5 * ksi1_1003[k]
                    + f_3 * pc_x[k] * ksk_1283[k];

        t_1599[k] = f_4 * ksi0_1004[k]
                    - f_5 * ksi1_1004[k]
                    + f_3 * pc_x[k] * ksk_1284[k];

        t_1600[k] = f_4 * ksi0_1005[k]
                    - f_5 * ksi1_1005[k]
                    + f_3 * pc_x[k] * ksk_1285[k];

        t_1601[k] = f_3 * pc_y[k] * ksk_1280[k];
    }

#pragma omp simd aligned(t_1602, t_1603, t_1604, t_1605, t_1606, t_1607, pc_x, ksi0_1007, \
                         ksi1_1007, ksk_1287, ksk_1288, ksk_1289, ksk_1290, ksk_1291, \
                         ksk_1292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1602[k] = f_4 * ksi0_1007[k]
                    - f_5 * ksi1_1007[k]
                    + f_3 * pc_x[k] * ksk_1287[k];

        t_1603[k] = f_3 * pc_x[k] * ksk_1288[k];

        t_1604[k] = f_3 * pc_x[k] * ksk_1289[k];

        t_1605[k] = f_3 * pc_x[k] * ksk_1290[k];

        t_1606[k] = f_3 * pc_x[k] * ksk_1291[k];

        t_1607[k] = f_3 * pc_x[k] * ksk_1292[k];
    }

#pragma omp simd aligned(t_1608, t_1609, t_1610, t_1611, t_1612, pc_x, pc_y, ksi0_1001, \
                         ksi0_1002, ksi1_1001, ksi1_1002, ksk_1288, ksk_1289, ksk_1293, \
                         ksk_1294, ksk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1608[k] = f_3 * pc_x[k] * ksk_1293[k];

        t_1609[k] = f_3 * pc_x[k] * ksk_1294[k];

        t_1610[k] = f_3 * pc_x[k] * ksk_1295[k];

        t_1611[k] = f_1 * ksi0_1001[k]
                    - f_2 * ksi1_1001[k]
                    + f_3 * pc_y[k] * ksk_1288[k];

        t_1612[k] = f_21 * ksi0_1002[k]
                    - f_22 * ksi1_1002[k]
                    + f_3 * pc_y[k] * ksk_1289[k];
    }

#pragma omp simd aligned(t_1613, t_1614, t_1615, pc_y, ksi0_1003, ksi0_1004, ksi0_1005, \
                         ksi1_1003, ksi1_1004, ksi1_1005, ksk_1290, ksk_1291, \
                         ksk_1292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1613[k] = f_12 * ksi0_1003[k]
                    - f_13 * ksi1_1003[k]
                    + f_3 * pc_y[k] * ksk_1290[k];

        t_1614[k] = f_10 * ksi0_1004[k]
                    - f_11 * ksi1_1004[k]
                    + f_3 * pc_y[k] * ksk_1291[k];

        t_1615[k] = f_8 * ksi0_1005[k]
                    - f_9 * ksi1_1005[k]
                    + f_3 * pc_y[k] * ksk_1292[k];
    }

#pragma omp simd aligned(t_1616, t_1617, t_1618, t_1619, pc_y, pc_z, isk_1007, ksi0_1006, \
                         ksi0_1007, ksi1_1006, ksi1_1007, ksk_1293, ksk_1294, \
                         ksk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1616[k] = f_6 * ksi0_1006[k]
                    - f_7 * ksi1_1006[k]
                    + f_3 * pc_y[k] * ksk_1293[k];

        t_1617[k] = f_4 * ksi0_1007[k]
                    - f_5 * ksi1_1007[k]
                    + f_3 * pc_y[k] * ksk_1294[k];

        t_1618[k] = f_3 * pc_y[k] * ksk_1295[k];

        t_1619[k] = f_0 * isk_1007[k]
                    + f_1 * ksi0_1007[k]
                    - f_2 * ksi1_1007[k]
                    + f_3 * pc_z[k] * ksk_1295[k];
    }
}

auto
compute_prim_ksl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t isl0, const size_t isk,
                                                   const size_t isl1, const size_t ksi0,
                                                   const size_t ksi1, const size_t ksk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ksl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksi0, ksi1, ksk, ncols,
                                                              gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, isl0, isk,
                                                              isl1, ksk, ncols, gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, isl0,
                                                               isk, isl1, ksi0, ksi1, ksk,
                                                               ncols, gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, isl0,
                                                               isk, isl1, ksi0, ksi1, ksk,
                                                               ncols, gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece12(buffer, target, pc, isk, ksi0,
                                                               ksi1, ksk, ncols, gamma, p, q);

    compute_prim_ksl_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, isl0,
                                                               isk, isl1, ksi0, ksi1, ksk,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
